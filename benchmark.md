# tadpole vs ripser benchmark

Comparing `tadpole` (Rust port) against the reference C++ `ripser` on a range of inputs.
Both produce byte-for-byte identical output. All runs use `--dim 2`.

## Environment

- CPU: 11th Gen Intel Core i7-1165G7 @ 2.80GHz
- OS: Linux 6.19.14-arch1-1 x86_64
- Rust: 1.93.0 (compiled with `cargo build --release`)
- C++ ripser: compiled with GCC 15.2.1 (default ripser Makefile flags)
- Benchmarking tool: hyperfine 1.20.0

## How to reproduce

Build tadpole:

```sh
git clone <repo>
cd tadpole
cargo build --release
```

Build ripser:

```sh
git clone https://github.com/Ripser/ripser
cd ripser
make
```

Run benchmarks:

```sh
# random16 (16 points, lower-distance matrix)
hyperfine --warmup 2 --runs 5 \
  './target/release/tadpole --format lower-distance --dim 2 examples/random16.lower_distance_matrix' \
  'ripser --format lower-distance --dim 2 examples/random16.lower_distance_matrix'

# sphere_3_192 (192 points, lower-distance matrix)
hyperfine --warmup 2 --runs 5 \
  './target/release/tadpole --format lower-distance --dim 2 examples/sphere_3_192.lower_distance_matrix' \
  'ripser --format lower-distance --dim 2 examples/sphere_3_192.lower_distance_matrix'

# o3_1024 (1024 points, point cloud)
hyperfine --warmup 2 --runs 5 \
  './target/release/tadpole --format point-cloud --dim 2 examples/o3_1024.txt' \
  'ripser --format point-cloud --dim 2 examples/o3_1024.txt'
```

## Results

### random16 — 16 points, lower-distance matrix

```
Benchmark 1: tadpole --format lower-distance --dim 2 random16.lower_distance_matrix
  Time (mean ± σ):       1.2 ms ±   0.1 ms    [User: 0.9 ms, System: 0.4 ms]
  Range (min … max):     1.1 ms …   1.4 ms    5 runs

Benchmark 2: ripser --format lower-distance --dim 2 random16.lower_distance_matrix
  Time (mean ± σ):       2.2 ms ±   0.2 ms    [User: 1.4 ms, System: 0.9 ms]
  Range (min … max):     2.0 ms …   2.4 ms    5 runs

Summary: tadpole 1.81x faster
```

> Note: both commands complete in < 5 ms; startup noise dominates. Ratio is not meaningful at this scale.

### sphere_3_192 — 192 points, lower-distance matrix

```
Benchmark 1: tadpole --format lower-distance --dim 2 sphere_3_192.lower_distance_matrix
  Time (mean ± σ):     494.7 ms ±   6.2 ms    [User: 483.2 ms, System: 10.3 ms]
  Range (min … max):   488.1 ms … 502.4 ms    5 runs

Benchmark 2: ripser --format lower-distance --dim 2 sphere_3_192.lower_distance_matrix
  Time (mean ± σ):     550.4 ms ±   8.5 ms    [User: 525.7 ms, System: 23.5 ms]
  Range (min … max):   541.0 ms … 561.9 ms    5 runs

Summary: tadpole 1.11x faster
```

### o3_1024 — 1024 points, point cloud

Single runs (not hyperfine):

```
tadpole: real 10m38.0s  user 9m37.7s  sys 0m59.0s
ripser:  real 12m18.4s  user 10m12.7s sys 2m4.4s
```

Output identical.

## Summary table

| Input        | Points | Format          | tadpole        | ripser         | ratio        |
|--------------|--------|-----------------|----------------|----------------|--------------|
| random16     | 16     | lower-distance  | 1.2 ms         | 2.2 ms         | 1.81x faster (startup-dominated) |
| sphere_3_192 | 192    | lower-distance  | 494.7 ms       | 550.4 ms       | 1.11x faster |
| o3_1024      | 1024   | point-cloud     | 9m37.7s (user) | 10m12.7s (user)| 1.06x faster (single run) |
