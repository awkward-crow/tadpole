/// Coboundary enumerators for compressed-lower and sparse distance matrices.
///
/// In C++, these were inner classes of ripser<DistanceMatrix> holding back-references
/// to the parent.  Here they hold only their own iteration state; the distance matrix
/// and binomial table are passed in as arguments to `set_simplex`, `has_next`, and `next`.
use crate::{
    binomial::BinomialCoeffTable,
    distance::{CompressedLowerDistanceMatrix, DistanceMatrix, SparseDistanceMatrix},
    types::{CoefficientT, DiameterEntry, IndexDiameter, IndexT},
};

// ---------------------------------------------------------------------------
// CompressedCoboundaryEnumerator
// ---------------------------------------------------------------------------

pub struct CompressedCoboundaryEnumerator {
    idx_below: IndexT,
    idx_above: IndexT,
    pub j: IndexT,
    pub k: IndexT,
    vertices: Vec<IndexT>,
    simplex: DiameterEntry,
    modulus: CoefficientT,
}

impl CompressedCoboundaryEnumerator {
    pub fn new(modulus: CoefficientT) -> Self {
        CompressedCoboundaryEnumerator {
            idx_below: 0,
            idx_above: 0,
            j: 0,
            k: 0,
            vertices: Vec::new(),
            simplex: DiameterEntry::invalid(),
            modulus,
        }
    }

    pub fn set_simplex(
        &mut self,
        simplex: DiameterEntry,
        dim: IndexT,
        n: IndexT,
        binomial_coeff: &BinomialCoeffTable,
        get_simplex_vertices: impl Fn(IndexT, IndexT, IndexT, &mut Vec<IndexT>),
    ) {
        self.idx_below = simplex.index();
        self.idx_above = 0;
        self.j = n - 1;
        self.k = dim + 1;
        self.simplex = simplex;
        self.vertices.resize((dim + 1) as usize, 0);
        get_simplex_vertices(simplex.index(), dim, n, &mut self.vertices);
        let _ = binomial_coeff; // stored for has_next / next
    }

    pub fn has_next(&self, all_cofacets: bool, binomial_coeff: &BinomialCoeffTable) -> bool {
        self.j >= self.k
            && (all_cofacets || binomial_coeff.get(self.j, self.k) > self.idx_below)
    }

    #[inline]
    pub fn has_next_compressed(&self, all_cofacets: bool, _dist: &CompressedLowerDistanceMatrix, binomial_coeff: &BinomialCoeffTable) -> bool {
        self.has_next(all_cofacets, binomial_coeff)
    }

    pub fn next(
        &mut self,
        dist: &CompressedLowerDistanceMatrix,
        binomial_coeff: &BinomialCoeffTable,
    ) -> DiameterEntry {
        while binomial_coeff.get(self.j, self.k) <= self.idx_below {
            self.idx_below -= binomial_coeff.get(self.j, self.k);
            self.idx_above += binomial_coeff.get(self.j, self.k + 1);
            self.j -= 1;
            self.k -= 1;
            debug_assert!(self.k != -1);
        }
        let mut cofacet_diameter = self.simplex.diameter;
        for &v in &self.vertices {
            cofacet_diameter = cofacet_diameter.max(dist.get(self.j, v));
        }
        let cofacet_index =
            self.idx_above + binomial_coeff.get(self.j, self.k + 1) + self.idx_below;
        let cofacet_coefficient = if self.k & 1 != 0 {
            (self.modulus - 1) * self.simplex.coefficient() % self.modulus
        } else {
            self.simplex.coefficient() % self.modulus
        };
        self.j -= 1;
        DiameterEntry::new(cofacet_diameter, cofacet_index, cofacet_coefficient)
    }

    #[inline]
    pub fn next_compressed(&mut self, dist: &CompressedLowerDistanceMatrix, binomial_coeff: &BinomialCoeffTable) -> DiameterEntry {
        self.next(dist, binomial_coeff)
    }
}

// ---------------------------------------------------------------------------
// SparseCoboundaryEnumerator
// ---------------------------------------------------------------------------

pub struct SparseCoboundaryEnumerator {
    idx_below: IndexT,
    idx_above: IndexT,
    pub k: IndexT,
    vertices: Vec<IndexT>,
    simplex: DiameterEntry,
    modulus: CoefficientT,
    // Per-vertex iterator positions into dist.neighbors[v]
    neighbor_pos: Vec<usize>,
    neighbor_end: Vec<usize>,
    // Current candidate neighbor
    neighbor: IndexDiameter,
}

impl SparseCoboundaryEnumerator {
    pub fn new(modulus: CoefficientT) -> Self {
        SparseCoboundaryEnumerator {
            idx_below: 0,
            idx_above: 0,
            k: 0,
            vertices: Vec::new(),
            simplex: DiameterEntry::invalid(),
            modulus,
            neighbor_pos: Vec::new(),
            neighbor_end: Vec::new(),
            neighbor: IndexDiameter::new(0, 0.0),
        }
    }

    pub fn set_simplex(
        &mut self,
        simplex: DiameterEntry,
        dim: IndexT,
        dist: &SparseDistanceMatrix,
        binomial_coeff: &BinomialCoeffTable,
        get_simplex_vertices: impl Fn(IndexT, IndexT, IndexT, &mut Vec<IndexT>),
    ) {
        let n = dist.size() as IndexT;
        self.idx_below = simplex.index();
        self.idx_above = 0;
        self.k = dim + 1;
        self.simplex = simplex;
        self.vertices.resize((dim + 1) as usize, 0);
        get_simplex_vertices(simplex.index(), dim, n, &mut self.vertices);

        self.neighbor_pos.resize((dim + 1) as usize, 0);
        self.neighbor_end.resize((dim + 1) as usize, 0);
        for i in 0..=(dim as usize) {
            let v = self.vertices[i] as usize;
            // neighbors are sorted ascending by index; we iterate in reverse
            // (rbegin..rend) to get descending order, matching C++
            self.neighbor_pos[i] = dist.neighbors[v].len(); // points past end; we'll use saturating sub
            self.neighbor_end[i] = 0;
        }
        // Initialise pos to point to last element (rbegin)
        for i in 0..=(dim as usize) {
            let v = self.vertices[i] as usize;
            self.neighbor_pos[i] = dist.neighbors[v].len();
        }
        let _ = binomial_coeff;
    }

    /// Returns true and sets self.neighbor if a valid cofacet candidate exists.
    /// Mirrors C++ `has_next(all_cofacets)`.
    pub fn has_next(
        &mut self,
        all_cofacets: bool,
        dist: &SparseDistanceMatrix,
        binomial_coeff: &BinomialCoeffTable,
    ) -> bool {
        let dim = self.vertices.len();
        // Iterate over neighbors of vertex[0] in descending index order
        'outer: loop {
            if self.neighbor_pos[0] == 0 {
                return false;
            }
            self.neighbor_pos[0] -= 1;
            let v0 = self.vertices[0] as usize;
            self.neighbor = dist.neighbors[v0][self.neighbor_pos[0]];

            // All other vertex neighbor lists must contain the same index
            for idx in 1..dim {
                let v = self.vertices[idx] as usize;
                // Advance this iterator until index <= neighbor.index
                while self.neighbor_pos[idx] > 0 {
                    let candidate = dist.neighbors[v][self.neighbor_pos[idx] - 1];
                    if candidate.index > self.neighbor.index {
                        self.neighbor_pos[idx] -= 1;
                    } else {
                        break;
                    }
                }
                // Check if we've hit the end
                if self.neighbor_pos[idx] == 0 {
                    return false;
                }
                let candidate = dist.neighbors[v][self.neighbor_pos[idx] - 1];
                if candidate.index != self.neighbor.index {
                    continue 'outer;
                }
                // Take max diameter
                if candidate.diameter > self.neighbor.diameter {
                    self.neighbor.diameter = candidate.diameter;
                }
            }

            // Advance k past any vertices larger than neighbor.index
            while self.k > 0
                && self.vertices[(self.k - 1) as usize] > self.neighbor.index
            {
                if !all_cofacets {
                    return false;
                }
                self.idx_below -= binomial_coeff.get(self.vertices[(self.k - 1) as usize], self.k);
                self.idx_above +=
                    binomial_coeff.get(self.vertices[(self.k - 1) as usize], self.k + 1);
                self.k -= 1;
            }
            return true;
        }
    }

    pub fn next(&mut self, binomial_coeff: &BinomialCoeffTable) -> DiameterEntry {
        // Advance iterator[0] (already consumed in has_next)
        let cofacet_diameter = self.simplex.diameter.max(self.neighbor.diameter);
        let cofacet_index =
            self.idx_above + binomial_coeff.get(self.neighbor.index, self.k + 1) + self.idx_below;
        let cofacet_coefficient = if self.k & 1 != 0 {
            (self.modulus - 1) * self.simplex.coefficient() % self.modulus
        } else {
            self.simplex.coefficient() % self.modulus
        };
        DiameterEntry::new(cofacet_diameter, cofacet_index, cofacet_coefficient)
    }

    #[inline]
    pub fn has_next_sparse(&mut self, all_cofacets: bool, dist: &SparseDistanceMatrix, binomial_coeff: &BinomialCoeffTable) -> bool {
        self.has_next(all_cofacets, dist, binomial_coeff)
    }

    #[inline]
    pub fn next_sparse(&mut self, binomial_coeff: &BinomialCoeffTable) -> DiameterEntry {
        self.next(binomial_coeff)
    }
}
