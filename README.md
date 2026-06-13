# tadpole

A Rust port of [ripser](https://github.com/Ripser/ripser), which computes
Vietoris-Rips persistence barcodes. 

## Build

```sh
cargo build --release
```

`print_persistence_pairs` is enabled by default (matches ripser's default build).
Optional features:

| Feature | Default | Description |
|---|---|---|
| `print_persistence_pairs` | on | Print barcode intervals to stdout |
| `indicate_progress` | off | Print reduction progress to stderr |
| `use_coefficients` | off | Compute with Z/p coefficients |

## Usage

```sh
tadpole [OPTIONS] [FILENAME]

Options:
  --format <FORMAT>    lower-distance, upper-distance, distance, point-cloud,
                       dipha, binary, sparse (default: distance)
  --dim <DIM>          Maximum homology dimension (default: 1)
  --threshold <T>      Distance threshold (default: enclosing radius)
  --ratio <R>          Skip intervals shorter than ratio * birth (default: 1)
  --modulus <P>        Coefficient modulus, prime (default: 2; use_coefficients feature only)
```

## Examples

[ripser](https://github.com/Ripser/ripser) has a number of example datasets. Copy them to this directory tree,

```sh
git clone https://github.com/Ripser/ripser /tmp/ripser
cp -r /tmp/ripser/examples examples
```

Then, run a small smoke test (after running `cargo build --release`),

```sh
./target/release/tadpole --format lower-distance --dim 2 examples/random16.lower_distance_matrix
```

The following are standard benchmarks,

```sh
./target/release/tadpole --format point-cloud --dim 2 examples/o3_1024.txt
./target/release/tadpole --format lower-distance --dim 2 examples/sphere_3_192.lower_distance_matrix
```

Real-world dataset (cyclooctane conformations)

```sh
./target/release/tadpole --format point-cloud --dim 2 examples/pointsCycloOctane.csv
```

See [benchmark.md](benchmark.md) for timing comparisons against C++ ripser.

## Module layout

```
src/
  main.rs        CLI entry point (clap), dispatch by distance matrix type
  types.rs       Primitive type aliases, Entry/DiameterEntry/DiameterIndex/IndexDiameter
  binomial.rs    BinomialCoeffTable
  distance.rs    DistanceMatrix trait + CompressedLower, CompressedUpper, Sparse, Euclidean
  union_find.rs  UnionFind
  matrix.rs      CompressedSparseMatrix<T>
  coboundary.rs  CoboundaryEnumerator (compressed and sparse variants)
  ripser.rs      Ripser<D> — core algorithm
  io.rs          File readers (lower/upper/distance/point-cloud/dipha/binary/sparse)
```

## Key translation decisions

- C++ templates → Rust generics + traits (`DistanceMatrix`)
- Self-referential coboundary enumerators → pass `&dist` / `&binomial_coeff` into
  `has_next` / `next` explicitly rather than storing back-references
- Pre-allocated enumerator buffers → passed as `&mut` arguments through the call stack
  (reuse allocation without storing them as fields on `Ripser`)
- `entry_t` bitfield → `Entry(i64)` newtype with manual bit ops
- `goto continue_outer` → Rust labeled `continue 'outer`
- `#ifdef` guards → Cargo features (`use_coefficients`, `indicate_progress`,
  `print_persistence_pairs`)
- `mutable` members → `&mut self` methods

## Performance opportunities

Not yet pursued. Listed roughly by expected impact.

### 1. HashMap hasher

`pivot_column_index` uses Rust's default `HashMap`, which uses SipHash — designed
for DoS resistance, not raw speed. C++ ripser uses robin-hood hashing. Swapping in
`rustc-hash` (`FxHashMap`) or `ahash` would likely be the single biggest win, since
the map is hit on every reduction step.

```toml
# Cargo.toml
rustc-hash = "2"
```

```rust
use rustc_hash::FxHashMap;
// replace HashMap<IndexT, usize> with FxHashMap<IndexT, usize>
```

### 2. BinaryHeap allocation per column

`working_reduction_column` and `working_coboundary` are heap-allocated fresh for
every column in `compute_pairs`. They could be fields on `Ripser`, cleared at the
start of each column, saving repeated allocator round-trips on large inputs.

### 3. next_simplices allocation per assemble call

`next_simplices` in `assemble_columns_to_reduce_compressed` is a fresh `Vec`
each call. Same fix as (2): make it a field and clear it.

### 4. Vertex-decode closure duplication

The binary-search vertex-decode closure is copy-pasted into three functions
(`get_zero_pivot_cofacet_compressed`, `assemble_columns_to_reduce_compressed`,
`init_coboundary_and_get_pivot_compressed`). Extracting it to an `#[inline]`
method would reduce code size and make future changes easier, though the compiler
likely already inlines all three copies.
