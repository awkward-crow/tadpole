use crate::types::{IndexT, check_overflow};

/// Precomputed table of binomial coefficients B[k][n] = C(n, k).
/// Layout mirrors the C++ version: outer index is k, inner index is n.
pub struct BinomialCoeffTable {
    b: Vec<Vec<IndexT>>,
}

impl BinomialCoeffTable {
    pub fn new(n: IndexT, k: IndexT) -> Self {
        let n = n as usize;
        let k = k as usize;
        let mut b = vec![vec![0_i64; n + 1]; k + 1];
        for i in 0..=n {
            b[0][i] = 1;
            for j in 1..i.min(k + 1) {
                b[j][i] = b[j - 1][i - 1] + b[j][i - 1];
            }
            if i <= k {
                b[i][i] = 1;
            }
            check_overflow(b[(i / 2).min(k)][i]);
        }
        BinomialCoeffTable { b }
    }

    /// Returns C(n, k).
    #[inline]
    pub fn get(&self, n: IndexT, k: IndexT) -> IndexT {
        debug_assert!(k >= 0 && (k as usize) < self.b.len());
        debug_assert!(n >= k - 1 && (n as usize) < self.b[k as usize].len());
        self.b[k as usize][n as usize]
    }
}
