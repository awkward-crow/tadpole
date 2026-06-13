use std::io::{BufRead, Read};

use crate::{
    distance::{
        CompressedLowerDistanceMatrix, CompressedUpperDistanceMatrix, EuclideanDistanceMatrix,
        SparseDistanceMatrix,
    },
    types::{IndexDiameter, IndexT, ValueT},
};

// ---------------------------------------------------------------------------
// Endian-aware primitive reader (matches C++ read<T>)
// ---------------------------------------------------------------------------

fn read_le<T: Copy + Default>(reader: &mut impl Read) -> Option<T> {
    let mut buf = vec![0u8; std::mem::size_of::<T>()];
    reader.read_exact(&mut buf).ok()?;
    // The C++ code reverses bytes on big-endian. We do the same.
    if cfg!(target_endian = "big") {
        buf.reverse();
    }
    Some(unsafe { std::ptr::read_unaligned(buf.as_ptr() as *const T) })
}

// ---------------------------------------------------------------------------
// Readers
// ---------------------------------------------------------------------------

/// Read comma/whitespace-separated lower-triangular distance matrix.
pub fn read_lower(reader: impl BufRead) -> CompressedLowerDistanceMatrix {
    let mut distances = Vec::new();
    for line in reader.lines() {
        let line = line.unwrap();
        for tok in line.split(|c: char| c == ',' || c.is_whitespace()) {
            let tok = tok.trim();
            if tok.is_empty() {
                continue;
            }
            if let Ok(v) = tok.parse::<ValueT>() {
                distances.push(v);
            }
        }
    }
    CompressedLowerDistanceMatrix::new(distances)
}

pub fn read_upper(reader: impl BufRead) -> CompressedLowerDistanceMatrix {
    let mut distances = Vec::new();
    for line in reader.lines() {
        let line = line.unwrap();
        for tok in line.split(|c: char| c == ',' || c.is_whitespace()) {
            let tok = tok.trim();
            if tok.is_empty() {
                continue;
            }
            if let Ok(v) = tok.parse::<ValueT>() {
                distances.push(v);
            }
        }
    }
    CompressedUpperDistanceMatrix::new(distances).to_lower()
}

/// Full distance matrix — only lower triangular part is read (j < i).
pub fn read_distance_matrix(reader: impl BufRead) -> CompressedLowerDistanceMatrix {
    let mut distances = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line.unwrap();
        let mut j = 0;
        for tok in line.split(|c: char| c == ',' || c.is_whitespace()) {
            let tok = tok.trim();
            if tok.is_empty() {
                continue;
            }
            if j < i {
                if let Ok(v) = tok.parse::<ValueT>() {
                    distances.push(v);
                }
            }
            j += 1;
        }
    }
    CompressedLowerDistanceMatrix::new(distances)
}

pub fn read_point_cloud(reader: impl BufRead) -> EuclideanDistanceMatrix {
    let mut points: Vec<Vec<ValueT>> = Vec::new();
    for line in reader.lines() {
        let line = line.unwrap();
        let point: Vec<ValueT> = line
            .split(|c: char| c == ',' || c.is_whitespace())
            .filter_map(|tok| {
                let tok = tok.trim();
                if tok.is_empty() { None } else { tok.parse().ok() }
            })
            .collect();
        if !point.is_empty() {
            if !points.is_empty() {
                assert_eq!(point.len(), points[0].len(), "point cloud dimension mismatch");
            }
            points.push(point);
        }
    }
    EuclideanDistanceMatrix::new(points)
}

pub fn read_sparse_distance_matrix(reader: impl BufRead) -> SparseDistanceMatrix {
    let mut neighbors: Vec<Vec<IndexDiameter>> = Vec::new();
    let mut num_edges: IndexT = 0;

    for line in reader.lines() {
        let line = line.unwrap();
        let mut parts = line.split(|c: char| c == ',' || c.is_whitespace()).filter(|s| !s.trim().is_empty());
        let i: usize = match parts.next().and_then(|s| s.trim().parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        let j: usize = match parts.next().and_then(|s| s.trim().parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        let value: ValueT = match parts.next().and_then(|s| s.trim().parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        if i != j {
            let needed = i.max(j) + 1;
            if neighbors.len() < needed {
                neighbors.resize(needed, Vec::new());
            }
            neighbors[i].push(IndexDiameter::new(j as IndexT, value));
            neighbors[j].push(IndexDiameter::new(i as IndexT, value));
            num_edges += 1;
        }
    }

    for row in &mut neighbors {
        row.sort();
    }

    SparseDistanceMatrix::new(neighbors, num_edges)
}

pub fn read_dipha(reader: &mut impl Read) -> CompressedLowerDistanceMatrix {
    let magic: i64 = read_le(reader).expect("dipha: short read");
    if magic != 8067171840 {
        eprintln!("input is not a Dipha file (magic number: 8067171840)");
        std::process::exit(-1);
    }
    let file_type: i64 = read_le(reader).expect("dipha: short read");
    if file_type != 7 {
        eprintln!("input is not a Dipha distance matrix (file type: 7)");
        std::process::exit(-1);
    }
    let n: i64 = read_le(reader).expect("dipha: short read");
    let mut distances = Vec::new();
    for i in 0..n {
        for j in 0..n {
            let v: f64 = read_le(reader).expect("dipha: short read");
            if i > j {
                distances.push(v as ValueT);
            }
        }
    }
    CompressedLowerDistanceMatrix::new(distances)
}

pub fn read_binary(reader: &mut impl Read) -> CompressedLowerDistanceMatrix {
    let mut distances = Vec::new();
    loop {
        match read_le::<ValueT>(reader) {
            Some(v) => distances.push(v),
            None => break,
        }
    }
    CompressedLowerDistanceMatrix::new(distances)
}
