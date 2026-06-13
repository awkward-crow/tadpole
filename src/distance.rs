use crate::types::{IndexT, ValueT, IndexDiameter};

// ---------------------------------------------------------------------------
// Trait
// ---------------------------------------------------------------------------

pub trait DistanceMatrix {
    fn get(&self, i: IndexT, j: IndexT) -> ValueT;
    fn size(&self) -> usize;
}

// ---------------------------------------------------------------------------
// CompressedLowerDistanceMatrix
// ---------------------------------------------------------------------------

/// Lower-triangular distance matrix stored in row-major order.
/// distances[i*(i-1)/2 + j] = d(i, j) for j < i.
pub struct CompressedLowerDistanceMatrix {
    pub distances: Vec<ValueT>,
    n: usize,
}

impl CompressedLowerDistanceMatrix {
    pub fn new(distances: Vec<ValueT>) -> Self {
        let n = ((1.0 + (1.0 + 8.0 * distances.len() as f64).sqrt()) / 2.0) as usize;
        assert_eq!(distances.len(), n * (n - 1) / 2);
        CompressedLowerDistanceMatrix { distances, n }
    }

    /// Build from any DistanceMatrix by reading the lower triangle.
    pub fn from_dense<D: DistanceMatrix>(mat: &D) -> Self {
        let n = mat.size();
        let mut distances = Vec::with_capacity(n * (n - 1) / 2);
        for i in 1..n {
            for j in 0..i {
                distances.push(mat.get(i as IndexT, j as IndexT));
            }
        }
        CompressedLowerDistanceMatrix { distances, n }
    }

    /// Row offset for row i (i >= 1).
    #[inline]
    fn row_start(&self, i: usize) -> usize {
        i * (i - 1) / 2
    }
}

impl DistanceMatrix for CompressedLowerDistanceMatrix {
    #[inline]
    fn get(&self, i: IndexT, j: IndexT) -> ValueT {
        let (i, j) = (i as usize, j as usize);
        if i == j {
            0.0
        } else if i < j {
            self.distances[self.row_start(j) + i]
        } else {
            self.distances[self.row_start(i) + j]
        }
    }

    fn size(&self) -> usize {
        self.n
    }
}

// ---------------------------------------------------------------------------
// CompressedUpperDistanceMatrix (converted to lower on construction)
// ---------------------------------------------------------------------------

pub struct CompressedUpperDistanceMatrix {
    distances: Vec<ValueT>,
    n: usize,
}

impl CompressedUpperDistanceMatrix {
    pub fn new(distances: Vec<ValueT>) -> Self {
        let n = ((1.0 + (1.0 + 8.0 * distances.len() as f64).sqrt()) / 2.0) as usize;
        assert_eq!(distances.len(), n * (n - 1) / 2);
        CompressedUpperDistanceMatrix { distances, n }
    }

    /// Convert to lower-triangular by reading the upper triangle and
    /// re-indexing, matching the C++ conversion constructor.
    pub fn to_lower(self) -> CompressedLowerDistanceMatrix {
        let n = self.n;
        let mut lower = vec![0.0_f32; n * (n - 1) / 2];
        // Upper stores d(i,j) for i < j at position: i*(n-1) - i*(i-1)/2 + (j - i - 1)
        // We write d(j,i) = d(i,j) into the lower matrix.
        for i in 0..n {
            for j in (i + 1)..n {
                let upper_idx = i * (2 * n - i - 3) / 2 + j - 1;
                let lower_idx = j * (j - 1) / 2 + i;
                lower[lower_idx] = self.distances[upper_idx];
            }
        }
        CompressedLowerDistanceMatrix { distances: lower, n }
    }
}

// ---------------------------------------------------------------------------
// SparseDistanceMatrix
// ---------------------------------------------------------------------------

pub struct SparseDistanceMatrix {
    pub neighbors: Vec<Vec<IndexDiameter>>,
    pub num_edges: IndexT,
}

impl SparseDistanceMatrix {
    pub fn new(neighbors: Vec<Vec<IndexDiameter>>, num_edges: IndexT) -> Self {
        SparseDistanceMatrix { neighbors, num_edges }
    }

    /// Build from a dense distance matrix, including only edges <= threshold.
    pub fn from_dense<D: DistanceMatrix>(mat: &D, threshold: ValueT) -> Self {
        let n = mat.size();
        let mut neighbors = vec![Vec::new(); n];
        let mut num_edges = 0_i64;
        for i in 0..n {
            for j in 0..n {
                if i != j {
                    let d = mat.get(i as IndexT, j as IndexT);
                    if d <= threshold {
                        num_edges += 1;
                        neighbors[i].push(IndexDiameter::new(j as IndexT, d));
                    }
                }
            }
        }
        SparseDistanceMatrix { neighbors, num_edges }
    }
}

impl DistanceMatrix for SparseDistanceMatrix {
    fn get(&self, i: IndexT, j: IndexT) -> ValueT {
        let row = &self.neighbors[i as usize];
        match row.binary_search_by_key(&j, |nd| nd.index) {
            Ok(pos) => row[pos].diameter,
            Err(_) => f32::INFINITY,
        }
    }

    fn size(&self) -> usize {
        self.neighbors.len()
    }
}

// ---------------------------------------------------------------------------
// EuclideanDistanceMatrix
// ---------------------------------------------------------------------------

pub struct EuclideanDistanceMatrix {
    pub points: Vec<Vec<ValueT>>,
}

impl EuclideanDistanceMatrix {
    pub fn new(points: Vec<Vec<ValueT>>) -> Self {
        if !points.is_empty() {
            let dim = points[0].len();
            for p in &points {
                assert_eq!(p.len(), dim);
            }
        }
        EuclideanDistanceMatrix { points }
    }
}

impl DistanceMatrix for EuclideanDistanceMatrix {
    fn get(&self, i: IndexT, j: IndexT) -> ValueT {
        let (i, j) = (i as usize, j as usize);
        self.points[i]
            .iter()
            .zip(&self.points[j])
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f32>()
            .sqrt()
    }

    fn size(&self) -> usize {
        self.points.len()
    }
}
