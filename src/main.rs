mod binomial;
mod coboundary;
mod distance;
mod io;
mod matrix;
mod ripser;
mod types;
mod union_find;

use types::fmt_float;

use std::{
    fs::File,
    io::{BufReader, Read},
};

use clap::{Parser, ValueEnum};
use distance::{CompressedLowerDistanceMatrix, DistanceMatrix, SparseDistanceMatrix};
use types::ValueT;

#[derive(Debug, Clone, ValueEnum)]
enum FileFormat {
    #[value(name = "lower-distance")]
    LowerDistance,
    #[value(name = "upper-distance")]
    UpperDistance,
    #[value(name = "distance")]
    Distance,
    #[value(name = "point-cloud")]
    PointCloud,
    Dipha,
    Sparse,
    Binary,
}

#[derive(Parser, Debug)]
#[command(about = "tadpole: Vietoris-Rips persistence barcodes (Rust port of ripser)")]
struct Args {
    /// Input file (stdin if omitted)
    filename: Option<String>,

    /// File format
    #[arg(long, default_value = "distance")]
    format: FileFormat,

    /// Maximum homology dimension
    #[arg(long, default_value_t = 1)]
    dim: i64,

    /// Maximum diameter threshold
    #[arg(long)]
    threshold: Option<f32>,

    /// Only show pairs with death/birth ratio > r
    #[arg(long, default_value_t = 1.0)]
    ratio: f32,

    /// Coefficient modulus (must be prime); only available with use_coefficients feature
    #[cfg(feature = "use_coefficients")]
    #[arg(long, default_value_t = 2)]
    modulus: u16,
}

fn main() {
    let args = Args::parse();

    let threshold = args.threshold.unwrap_or(f32::MAX);
    #[cfg(feature = "use_coefficients")]
    let modulus: u16 = args.modulus;
    #[cfg(not(feature = "use_coefficients"))]
    let modulus: u16 = 2;

    // Open input
    let stdin = std::io::stdin();
    let reader: Box<dyn std::io::Read> = match &args.filename {
        Some(path) => Box::new(File::open(path).unwrap_or_else(|_| {
            eprintln!("couldn't open file {}", path);
            std::process::exit(-1);
        })),
        None => Box::new(stdin.lock()),
    };

    match args.format {
        FileFormat::Sparse => {
            let dist = io::read_sparse_distance_matrix(BufReader::new(reader));
            println!(
                "sparse distance matrix with {} points and {}/{} entries",
                dist.size(),
                dist.num_edges,
                (dist.size() * (dist.size() - 1)) / 2
            );
            ripser::Ripser::new(dist, args.dim, threshold, args.ratio, modulus)
                .compute_barcodes();
        }
        FileFormat::PointCloud if threshold < f32::MAX => {
            let eucl = io::read_point_cloud(BufReader::new(reader));
            println!(
                "point cloud with {} points in dimension {}",
                eucl.size(),
                eucl.points.first().map_or(0, |p| p.len())
            );
            let sparse = SparseDistanceMatrix::from_dense(&eucl, threshold);
            ripser::Ripser::new(sparse, args.dim, threshold, args.ratio, modulus)
                .compute_barcodes();
        }
        FileFormat::PointCloud => {
            let eucl = io::read_point_cloud(BufReader::new(reader));
            println!(
                "point cloud with {} points in dimension {}",
                eucl.size(),
                eucl.points.first().map_or(0, |p| p.len())
            );
            let lower = CompressedLowerDistanceMatrix::from_dense(&eucl);
            run_compressed(lower, args.dim, threshold, args.ratio, modulus);
        }
        fmt => {
            let mut buf = Vec::new();
            BufReader::new(reader).read_to_end(&mut buf).unwrap();
            let lower = match fmt {
                FileFormat::LowerDistance => {
                    io::read_lower(BufReader::new(buf.as_slice()))
                }
                FileFormat::UpperDistance => {
                    io::read_upper(BufReader::new(buf.as_slice()))
                }
                FileFormat::Distance => {
                    io::read_distance_matrix(BufReader::new(buf.as_slice()))
                }
                FileFormat::Dipha => io::read_dipha(&mut buf.as_slice()),
                FileFormat::Binary => io::read_binary(&mut buf.as_slice()),
                _ => unreachable!(),
            };
            run_compressed(lower, args.dim, threshold, args.ratio, modulus);
        }
    }
}

fn run_compressed(
    dist: CompressedLowerDistanceMatrix,
    dim_max: i64,
    threshold: ValueT,
    ratio: f32,
    modulus: u16,
) {
    let n = dist.size();
    let mut min = f32::INFINITY;
    let mut max = f32::NEG_INFINITY;
    let mut max_finite = f32::NEG_INFINITY;
    let mut num_edges = 0usize;
    for &d in &dist.distances {
        if d < min { min = d; }
        if d > max { max = d; }
        if d != f32::INFINITY && d > max_finite { max_finite = d; }
        if d <= threshold { num_edges += 1; }
    }
    println!("value range: [{},{}]", fmt_float(min), fmt_float(max_finite));

    if threshold == f32::MAX {
        let mut enclosing_radius = f32::INFINITY;
        for i in 0..n {
            let mut r_i = f32::NEG_INFINITY;
            for j in 0..n {
                let d = dist.get(i as i64, j as i64);
                if d > r_i { r_i = d; }
            }
            if r_i < enclosing_radius { enclosing_radius = r_i; }
        }
        println!(
            "distance matrix with {} points, using threshold at enclosing radius {}",
            n, fmt_float(enclosing_radius)
        );
        ripser::Ripser::new(dist, dim_max, enclosing_radius, ratio, modulus).compute_barcodes();
    } else {
        println!(
            "sparse distance matrix with {} points and {}/{} entries",
            n, num_edges,
            n * (n - 1) / 2
        );
        let sparse = SparseDistanceMatrix::from_dense(&dist, threshold);
        ripser::Ripser::new(sparse, dim_max, threshold, ratio, modulus).compute_barcodes();
    }
}
