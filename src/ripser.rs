use std::collections::BinaryHeap;

use crate::{
    binomial::BinomialCoeffTable,
    coboundary::{CompressedCoboundaryEnumerator, SparseCoboundaryEnumerator},
    distance::{CompressedLowerDistanceMatrix, DistanceMatrix, SparseDistanceMatrix},
    matrix::CompressedSparseMatrix,
    types::{fmt_float, CoefficientT, DiameterEntry, DiameterIndex, IndexT, ValueT},
    union_find::UnionFind,
};
use std::collections::HashMap;

#[cfg(feature = "use_coefficients")]
fn is_prime(n: CoefficientT) -> bool {
    if n < 2 { return false; }
    if n == 2 { return true; }
    if n & 1 == 0 { return false; }
    let mut p = 3u32;
    while p * p <= n as u32 {
        if n % p as CoefficientT == 0 { return false; }
        p += 2;
    }
    true
}

#[cfg(feature = "use_coefficients")]
fn multiplicative_inverse_vector(m: CoefficientT) -> Vec<CoefficientT> {
    let mut inv = vec![0u16; m as usize];
    if m > 1 {
        inv[1] = 1;
        for a in 2..m as usize {
            inv[a] = (m - (inv[m as usize % a] as u32 * (m as u32 / a as u32) % m as u32) as CoefficientT) % m;
        }
    }
    inv
}

// ---------------------------------------------------------------------------
// Concrete Ripser<D>
// ---------------------------------------------------------------------------

pub struct Ripser<D: DistanceMatrix> {
    pub dist: D,
    pub n: IndexT,
    pub dim_max: IndexT,
    pub threshold: ValueT,
    pub ratio: f32,
    pub modulus: CoefficientT,
    pub binomial_coeff: BinomialCoeffTable,
    #[cfg(feature = "use_coefficients")]
    pub multiplicative_inverse: Vec<CoefficientT>,
    // Scratch buffers (replace C++ mutable members and static locals)
    pub cofacet_entries: Vec<DiameterEntry>,
    pub vertices: Vec<IndexT>,
}

impl<D: DistanceMatrix> Ripser<D> {
    pub fn new(dist: D, dim_max: IndexT, threshold: ValueT, ratio: f32, modulus: CoefficientT) -> Self {
        #[cfg(feature = "use_coefficients")]
        if !is_prime(modulus) {
            eprintln!("modulus {} is not prime", modulus);
            std::process::exit(-1);
        }
        let n = dist.size() as IndexT;
        let dim_max = dim_max.min(n - 2);
        Ripser {
            binomial_coeff: BinomialCoeffTable::new(n, dim_max + 2),
            #[cfg(feature = "use_coefficients")]
            multiplicative_inverse: multiplicative_inverse_vector(modulus),
            dist,
            n,
            dim_max,
            threshold,
            ratio,
            modulus,
            cofacet_entries: Vec::new(),
            vertices: Vec::new(),
        }
    }

    pub fn get_max_vertex(&self, idx: IndexT, k: IndexT, n: IndexT) -> IndexT {
        // Binary search: find largest w >= k-1 such that C(w,k) <= idx
        let mut top = n;
        let bottom = k - 1;
        if self.binomial_coeff.get(top, k) > idx {
            let mut count = top - bottom;
            while count > 0 {
                let step = count >> 1;
                let mid = top - step;
                if self.binomial_coeff.get(mid, k) > idx {
                    top = mid - 1;
                    count -= step + 1;
                } else {
                    count = step;
                }
            }
        }
        top
    }

    pub fn get_edge_index(&self, i: IndexT, j: IndexT) -> IndexT {
        self.binomial_coeff.get(i, 2) + j
    }

    pub fn compute_diameter(&mut self, index: IndexT, dim: IndexT) -> ValueT {
        let mut diam = f32::NEG_INFINITY;
        self.vertices.resize((dim + 1) as usize, 0);
        // get_simplex_vertices fills largest-first
        let n = self.n;
        let bc = &self.binomial_coeff;
        // inline to avoid borrow conflict
        {
            let mut idx = index;
            let mut nn = n - 1;
            for k in (2..=dim + 1).rev() {
                nn = {
                    let mut top = nn;
                    if bc.get(top, k) > idx {
                        let mut count = top - (k - 1);
                        while count > 0 {
                            let step = count >> 1;
                            let mid = top - step;
                            if bc.get(mid, k) > idx {
                                top = mid - 1;
                                count -= step + 1;
                            } else {
                                count = step;
                            }
                        }
                    }
                    top
                };
                self.vertices[(dim + 1 - k) as usize] = nn;
                idx -= bc.get(nn, k);
            }
            self.vertices[dim as usize] = idx;
        }
        for i in 0..=(dim as usize) {
            for j in 0..i {
                diam = diam.max(self.dist.get(self.vertices[i], self.vertices[j]));
            }
        }
        diam
    }
}

// ---------------------------------------------------------------------------
// Ripser<CompressedLowerDistanceMatrix>
// ---------------------------------------------------------------------------

impl Ripser<CompressedLowerDistanceMatrix> {
    fn get_edges_impl(&self) -> Vec<DiameterIndex> {
        let mut edges = Vec::new();
        let total = self.binomial_coeff.get(self.n, 2);
        let mut index = total;
        while index > 0 {
            index -= 1;
            // get vertices of edge `index`
            let v1 = self.get_max_vertex(index, 2, self.n - 1);
            let v0 = index - self.binomial_coeff.get(v1, 2);
            let length = self.dist.get(v0, v1);
            if length <= self.threshold {
                edges.push(DiameterIndex::new(length, index));
            }
        }
        edges
    }

    pub fn compute_barcodes(&mut self) {
        let mut simplices: Vec<DiameterIndex> = Vec::new();
        let mut columns_to_reduce: Vec<DiameterIndex> = Vec::new();

        let mut pivot_cofacet_enum = CompressedCoboundaryEnumerator::new(self.modulus);
        let mut assemble_cofacet_enum = CompressedCoboundaryEnumerator::new(self.modulus);
        let mut init_cofacet_enum = CompressedCoboundaryEnumerator::new(self.modulus);
        let mut add_cofacet_enum = CompressedCoboundaryEnumerator::new(self.modulus);
        let mut apparent_cofacet_enum = CompressedCoboundaryEnumerator::new(self.modulus);

        self.compute_dim_0_pairs_compressed(
            &mut simplices,
            &mut columns_to_reduce,
            &mut pivot_cofacet_enum,
        );

        for dim in 1..=self.dim_max {
            let mut pivot_column_index: HashMap<IndexT, usize> = HashMap::new();
            pivot_column_index.reserve(columns_to_reduce.len());

            self.compute_pairs_compressed(
                &columns_to_reduce,
                &mut pivot_column_index,
                dim,
                &mut init_cofacet_enum,
                &mut add_cofacet_enum,
                &mut apparent_cofacet_enum,
            );

            if dim < self.dim_max {
                self.assemble_columns_to_reduce_compressed(
                    &mut simplices,
                    &mut columns_to_reduce,
                    &pivot_column_index,
                    dim + 1,
                    &mut assemble_cofacet_enum,
                    &mut apparent_cofacet_enum,
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // get_zero_pivot_facet / cofacet / apparent helpers
    // -----------------------------------------------------------------------

    fn get_zero_pivot_facet_inline(&mut self, simplex: DiameterEntry, dim: IndexT) -> DiameterEntry {
        // Iterate over boundary faces using the same logic as simplex_boundary_enumerator
        let idx = simplex.index();
        let diam = simplex.diameter;
        let modulus = self.modulus;
        let coeff = simplex.coefficient();

        let mut idx_below = idx;
        let mut idx_above: IndexT = 0;
        let mut j = self.n - 1;
        let mut k = dim;

        while k >= 0 {
            j = self.get_max_vertex(idx_below, k + 1, j);
            let face_index = idx_above - self.binomial_coeff.get(j, k + 1) + idx_below;
            let face_diam = self.compute_diameter(face_index, dim - 1);
            let face_coeff = if k & 1 != 0 {
                (modulus - 1 + coeff) % modulus
            } else {
                coeff % modulus
            };
            idx_below -= self.binomial_coeff.get(j, k + 1);
            idx_above += self.binomial_coeff.get(j, k);
            k -= 1;

            if face_diam == diam {
                return DiameterEntry::new(face_diam, face_index, face_coeff);
            }
        }
        DiameterEntry::invalid()
    }

    fn get_zero_pivot_cofacet_compressed(
        &mut self,
        simplex: DiameterEntry,
        dim: IndexT,
        cofacets: &mut CompressedCoboundaryEnumerator,
    ) -> DiameterEntry {
        let n = self.n;
        let dist = &self.dist as *const CompressedLowerDistanceMatrix;
        let bc = &self.binomial_coeff as *const BinomialCoeffTable;
        let dist = unsafe { &*dist };
        let bc = unsafe { &*bc };

        let get_verts = |idx: IndexT, d: IndexT, nn: IndexT, out: &mut Vec<IndexT>| {
            // Inline get_simplex_vertices
            out.resize((d + 1) as usize, 0);
            let mut i = idx;
            let mut nn2 = nn - 1;
            for kk in (2..=d + 1).rev() {
                nn2 = {
                    let mut top = nn2;
                    if bc.get(top, kk) > i {
                        let mut count = top - (kk - 1);
                        while count > 0 {
                            let step = count >> 1;
                            let mid = top - step;
                            if bc.get(mid, kk) > i {
                                top = mid - 1;
                                count -= step + 1;
                            } else {
                                count = step;
                            }
                        }
                    }
                    top
                };
                out[(d + 1 - kk) as usize] = nn2;
                i -= bc.get(nn2, kk);
            }
            out[d as usize] = i;
        };

        cofacets.set_simplex(simplex, dim, n, bc, get_verts);
        while cofacets.has_next_compressed(true, dist, bc) {
            let cofacet = cofacets.next_compressed(dist, bc);
            if cofacet.diameter == simplex.diameter {
                return cofacet;
            }
        }
        DiameterEntry::invalid()
    }

    fn get_zero_apparent_facet(
        &mut self,
        simplex: DiameterEntry,
        dim: IndexT,
        cofacet_enum: &mut CompressedCoboundaryEnumerator,
    ) -> DiameterEntry {
        let facet = self.get_zero_pivot_facet_inline(simplex, dim);
        if facet.index() == -1 {
            return DiameterEntry::invalid();
        }
        let cofacet_of_facet = self.get_zero_pivot_cofacet_compressed(facet, dim - 1, cofacet_enum);
        if cofacet_of_facet.index() == simplex.index() {
            facet
        } else {
            DiameterEntry::invalid()
        }
    }

    fn get_zero_apparent_cofacet(
        &mut self,
        simplex: DiameterEntry,
        dim: IndexT,
        cofacet_enum: &mut CompressedCoboundaryEnumerator,
    ) -> DiameterEntry {
        let cofacet = self.get_zero_pivot_cofacet_compressed(simplex, dim, cofacet_enum);
        if cofacet.index() == -1 {
            return DiameterEntry::invalid();
        }
        let facet_of_cofacet = self.get_zero_pivot_facet_inline(cofacet, dim + 1);
        if facet_of_cofacet.index() == simplex.index() {
            cofacet
        } else {
            DiameterEntry::invalid()
        }
    }

    fn is_in_zero_apparent_pair(
        &mut self,
        simplex: DiameterEntry,
        dim: IndexT,
        cofacet_enum: &mut CompressedCoboundaryEnumerator,
    ) -> bool {
        let cofacet = self.get_zero_apparent_cofacet(simplex, dim, cofacet_enum);
        if cofacet.index() != -1 {
            return true;
        }
        let facet = self.get_zero_apparent_facet(simplex, dim, cofacet_enum);
        facet.index() != -1
    }

    // -----------------------------------------------------------------------
    // compute_dim_0_pairs
    // -----------------------------------------------------------------------

    fn compute_dim_0_pairs_compressed(
        &mut self,
        simplices: &mut Vec<DiameterIndex>,
        columns_to_reduce: &mut Vec<DiameterIndex>,
        cofacet_enum: &mut CompressedCoboundaryEnumerator,
    ) {
        #[cfg(feature = "print_persistence_pairs")]
        println!("persistence intervals in dim 0:");

        let mut dset = UnionFind::new(self.n);
        let mut edges = self.get_edges_impl();
        edges.sort_by(|a, b| b.cmp(a)); // descending diameter, ascending index
        // Actually DiameterIndex Ord is already "larger diameter first, smaller index first"
        // For sort in C++ it uses greater_diameter_or_smaller_index as comparator for rbegin/rend sort
        // which means sort ascending and reverse = sort descending by that comparator
        // Let's re-sort properly:
        edges.sort();
        edges.reverse();

        for e in &edges {
            let idx = e.index;
            let v1 = self.get_max_vertex(idx, 2, self.n - 1);
            let v0 = idx - self.binomial_coeff.get(v1, 2);
            let u = dset.find(v0);
            let v = dset.find(v1);
            if u != v {
                #[cfg(feature = "print_persistence_pairs")]
                if e.diameter != 0.0 {
                    println!(" [0,{})", fmt_float(e.diameter));
                }
                dset.link(u, v);
            } else if self.dim_max > 0 {
                // Check for zero-apparent cofacet
                let de = DiameterEntry::from_diameter_index(*e, 1);
                let cofacet = self.get_zero_apparent_cofacet(de, 1, cofacet_enum);
                if cofacet.index() == -1 {
                    columns_to_reduce.push(*e);
                }
            }
        }

        if self.dim_max > 0 {
            columns_to_reduce.reverse();
        }

        // In C++, the "edges" parameter doubles as "simplices" — it's populated
        // with all edges inside compute_dim_0_pairs and passed to
        // assemble_columns_to_reduce as the seed simplices for dim-2+ cofacet enumeration.
        *simplices = edges;

        #[cfg(feature = "print_persistence_pairs")]
        for i in 0..self.n {
            if dset.find(i) == i {
                println!(" [0, )");
            }
        }
    }

    // -----------------------------------------------------------------------
    // assemble_columns_to_reduce
    // -----------------------------------------------------------------------

    fn assemble_columns_to_reduce_compressed(
        &mut self,
        simplices: &mut Vec<DiameterIndex>,
        columns_to_reduce: &mut Vec<DiameterIndex>,
        pivot_column_index: &HashMap<IndexT, usize>,
        dim: IndexT,
        cofacets: &mut CompressedCoboundaryEnumerator,
        apparent_cofacet_enum: &mut CompressedCoboundaryEnumerator,
    ) {
        columns_to_reduce.clear();
        let mut next_simplices: Vec<DiameterIndex> = Vec::new();

        let n = self.n;
        let threshold = self.threshold;
        let dim_max = self.dim_max;

        for &simplex in simplices.iter() {
            let de = DiameterEntry::from_diameter_index(simplex, 1);
            let bc = &self.binomial_coeff as *const BinomialCoeffTable;
            let bc = unsafe { &*bc };
            let dist = &self.dist as *const CompressedLowerDistanceMatrix;
            let dist = unsafe { &*dist };

            cofacets.set_simplex(de, dim - 1, n, bc, |idx, d, nn, out| {
                out.resize((d + 1) as usize, 0);
                let mut i = idx;
                let mut nn2 = nn - 1;
                for kk in (2..=d + 1).rev() {
                    nn2 = {
                        let mut top = nn2;
                        if bc.get(top, kk) > i {
                            let mut count = top - (kk - 1);
                            while count > 0 {
                                let step = count >> 1;
                                let mid = top - step;
                                if bc.get(mid, kk) > i {
                                    top = mid - 1;
                                    count -= step + 1;
                                } else {
                                    count = step;
                                }
                            }
                        }
                        top
                    };
                    out[(d + 1 - kk) as usize] = nn2;
                    i -= bc.get(nn2, kk);
                }
                out[d as usize] = i;
            });

            while cofacets.has_next_compressed(false, dist, bc) {
                let cofacet = cofacets.next_compressed(dist, bc);
                if cofacet.diameter <= threshold {
                    if dim < dim_max {
                        next_simplices.push(DiameterIndex::new(cofacet.diameter, cofacet.index()));
                    }
                    let in_pair = self.is_in_zero_apparent_pair(
                        cofacet, dim, apparent_cofacet_enum
                    );
                    let in_pivot = pivot_column_index.contains_key(&cofacet.index());
                    if !in_pair && !in_pivot {
                        columns_to_reduce.push(DiameterIndex::new(cofacet.diameter, cofacet.index()));
                    }
                }
            }
        }

        if dim < dim_max {
            *simplices = next_simplices;
        }

        // C++ sorts with begin/end + greater_diameter_or_smaller_index → largest first.
        // DiameterIndex::Ord makes larger diameter "less", so sort() alone gives largest first.
        columns_to_reduce.sort();
    }

    // -----------------------------------------------------------------------
    // pop_pivot / get_pivot helpers
    // -----------------------------------------------------------------------

    fn pop_pivot(
        column: &mut BinaryHeap<DiameterEntry>,
        #[cfg(feature = "use_coefficients")] modulus: CoefficientT,
    ) -> DiameterEntry {
        #[cfg(feature = "use_coefficients")]
        {
            let mut pivot = DiameterEntry::invalid();
            while let Some(top) = column.pop() {
                if pivot.coefficient() == 0 {
                    pivot = top;
                } else if pivot.index() != top.index() {
                    column.push(top);
                    return pivot;
                } else {
                    let new_coeff = (pivot.coefficient() + top.coefficient()) % modulus;
                    pivot.set_coefficient(new_coeff);
                }
            }
            if pivot.coefficient() == 0 {
                DiameterEntry::invalid()
            } else {
                pivot
            }
        }
        #[cfg(not(feature = "use_coefficients"))]
        {
            while let Some(pivot) = column.pop() {
                match column.peek() {
                    Some(next) if next.index() == pivot.index() => {
                        column.pop();
                        // cancel: continue loop
                    }
                    _ => return pivot,
                }
            }
            DiameterEntry::invalid()
        }
    }

    fn get_pivot(column: &mut BinaryHeap<DiameterEntry>,
                 #[cfg(feature = "use_coefficients")] modulus: CoefficientT) -> DiameterEntry {
        #[cfg(feature = "use_coefficients")]
        let result = Self::pop_pivot(column, modulus);
        #[cfg(not(feature = "use_coefficients"))]
        let result = Self::pop_pivot(column);
        if result.index() != -1 {
            column.push(result);
        }
        result
    }

    // -----------------------------------------------------------------------
    // compute_pairs
    // -----------------------------------------------------------------------

    fn compute_pairs_compressed(
        &mut self,
        columns_to_reduce: &[DiameterIndex],
        pivot_column_index: &mut HashMap<IndexT, usize>,
        dim: IndexT,
        init_cofacets: &mut CompressedCoboundaryEnumerator,
        add_cofacets: &mut CompressedCoboundaryEnumerator,
        apparent_cofacet_enum: &mut CompressedCoboundaryEnumerator,
    ) {
        #[cfg(feature = "print_persistence_pairs")]
        println!("persistence intervals in dim {}:", dim);

        let mut reduction_matrix: CompressedSparseMatrix<DiameterEntry> =
            CompressedSparseMatrix::new();

        for (index_column_to_reduce, &col_di) in columns_to_reduce.iter().enumerate() {
            let column_to_reduce = DiameterEntry::from_diameter_index(col_di, 1);
            let diameter = col_di.diameter;

            reduction_matrix.append_column();

            let mut working_reduction_column: BinaryHeap<DiameterEntry> = BinaryHeap::new();
            let mut working_coboundary: BinaryHeap<DiameterEntry> = BinaryHeap::new();

            // init_coboundary_and_get_pivot
            let mut pivot = self.init_coboundary_and_get_pivot_compressed(
                column_to_reduce,
                &mut working_coboundary,
                dim,
                pivot_column_index,
                init_cofacets,
                apparent_cofacet_enum,
            );

            loop {
                if pivot.index() != -1 {
                    if let Some(&index_column_to_add) = pivot_column_index.get(&pivot.index()) {
                        // Found existing pivot — add that column
                        let factor = {
                            #[cfg(feature = "use_coefficients")]
                            {
                                let other_coeff = {
                                    // Look up the coefficient of the pivot in the index
                                    // In C++ this is pair->first (the entry_t key which holds coeff)
                                    // We store only the index in our HashMap, not the coefficient.
                                    // For Z/2 this is always 1.
                                    1u16
                                };
                                (self.modulus
                                    - pivot.coefficient()
                                        * self.multiplicative_inverse[other_coeff as usize]
                                        % self.modulus)
                                    % self.modulus
                            }
                            #[cfg(not(feature = "use_coefficients"))]
                            1u16
                        };

                        self.add_coboundary_compressed(
                            &reduction_matrix,
                            columns_to_reduce,
                            index_column_to_add,
                            factor,
                            dim,
                            &mut working_reduction_column,
                            &mut working_coboundary,
                            add_cofacets,
                        );

                        #[cfg(feature = "use_coefficients")]
                        { pivot = Self::get_pivot(&mut working_coboundary, self.modulus); }
                        #[cfg(not(feature = "use_coefficients"))]
                        { pivot = Self::get_pivot(&mut working_coboundary); }
                    } else {
                        // Check for zero apparent facet
                        let e = self.get_zero_apparent_facet(pivot, dim + 1, apparent_cofacet_enum);
                        if e.index() != -1 {
                            let mut e2 = e;
                            e2.set_coefficient((self.modulus - e2.coefficient()) % self.modulus);
                            self.add_simplex_coboundary_compressed(
                                e2,
                                dim,
                                &mut working_reduction_column,
                                &mut working_coboundary,
                                add_cofacets,
                            );
                            #[cfg(feature = "use_coefficients")]
                            { pivot = Self::get_pivot(&mut working_coboundary, self.modulus); }
                            #[cfg(not(feature = "use_coefficients"))]
                            { pivot = Self::get_pivot(&mut working_coboundary); }
                        } else {
                            // New pivot — record persistence pair
                            let death = pivot.diameter;
                            #[cfg(feature = "print_persistence_pairs")]
                            if death > diameter * self.ratio {
                                println!(" [{},{})", fmt_float(diameter), fmt_float(death));
                            }

                            pivot_column_index.insert(pivot.index(), index_column_to_reduce);

                            // Drain working_reduction_column into reduction_matrix
                            loop {
                                #[cfg(feature = "use_coefficients")]
                                let e = Self::pop_pivot(&mut working_reduction_column, self.modulus);
                                #[cfg(not(feature = "use_coefficients"))]
                                let e = Self::pop_pivot(&mut working_reduction_column);
                                if e.index() == -1 { break; }
                                debug_assert!(e.coefficient() > 0);
                                reduction_matrix.push_back(e);
                            }
                            break;
                        }
                    }
                } else {
                    // Pivot is -1 → infinite bar
                    #[cfg(feature = "print_persistence_pairs")]
                    println!(" [{}, )", fmt_float(diameter));
                    break;
                }
            }
        }
    }

    fn init_coboundary_and_get_pivot_compressed(
        &mut self,
        simplex: DiameterEntry,
        working_coboundary: &mut BinaryHeap<DiameterEntry>,
        dim: IndexT,
        pivot_column_index: &HashMap<IndexT, usize>,
        cofacets: &mut CompressedCoboundaryEnumerator,
        apparent_cofacet_enum: &mut CompressedCoboundaryEnumerator,
    ) -> DiameterEntry {
        let mut check_emergent = true;
        self.cofacet_entries.clear();

        let n = self.n;
        let threshold = self.threshold;
        let bc = &self.binomial_coeff as *const BinomialCoeffTable;
        let bc = unsafe { &*bc };
        let dist = &self.dist as *const CompressedLowerDistanceMatrix;
        let dist = unsafe { &*dist };

        cofacets.set_simplex(simplex, dim, n, bc, |idx, d, nn, out| {
            out.resize((d + 1) as usize, 0);
            let mut i = idx;
            let mut nn2 = nn - 1;
            for kk in (2..=d + 1).rev() {
                nn2 = {
                    let mut top = nn2;
                    if bc.get(top, kk) > i {
                        let mut count = top - (kk - 1);
                        while count > 0 {
                            let step = count >> 1;
                            let mid = top - step;
                            if bc.get(mid, kk) > i {
                                top = mid - 1;
                                count -= step + 1;
                            } else {
                                count = step;
                            }
                        }
                    }
                    top
                };
                out[(d + 1 - kk) as usize] = nn2;
                i -= bc.get(nn2, kk);
            }
            out[d as usize] = i;
        });

        while cofacets.has_next_compressed(true, dist, bc) {
            let cofacet = cofacets.next_compressed(dist, bc);
            if cofacet.diameter <= threshold {
                self.cofacet_entries.push(cofacet);
                if check_emergent && cofacet.diameter == simplex.diameter {
                    let not_in_pivot = !pivot_column_index.contains_key(&cofacet.index());
                    let apparent_facet = self.get_zero_apparent_facet(
                        cofacet, dim + 1, apparent_cofacet_enum
                    );
                    if not_in_pivot && apparent_facet.index() == -1 {
                        return cofacet;
                    }
                    check_emergent = false;
                }
            }
        }

        for &cofacet in &self.cofacet_entries {
            working_coboundary.push(cofacet);
        }

        #[cfg(feature = "use_coefficients")]
        return Self::get_pivot(working_coboundary, self.modulus);
        #[cfg(not(feature = "use_coefficients"))]
        return Self::get_pivot(working_coboundary);
    }

    fn add_simplex_coboundary_compressed(
        &self,
        simplex: DiameterEntry,
        dim: IndexT,
        working_reduction_column: &mut BinaryHeap<DiameterEntry>,
        working_coboundary: &mut BinaryHeap<DiameterEntry>,
        cofacets: &mut CompressedCoboundaryEnumerator,
    ) {
        working_reduction_column.push(simplex);
        let n = self.n;
        let threshold = self.threshold;
        let bc = &self.binomial_coeff as *const BinomialCoeffTable;
        let bc_ref = unsafe { &*bc };
        let dist = &self.dist as *const CompressedLowerDistanceMatrix;
        let dist_ref = unsafe { &*dist };

        cofacets.set_simplex(simplex, dim, n, bc_ref, |idx, d, nn, out| {
            out.resize((d + 1) as usize, 0);
            let mut i = idx;
            let mut nn2 = nn - 1;
            for kk in (2..=d + 1).rev() {
                nn2 = {
                    let mut top = nn2;
                    if bc_ref.get(top, kk) > i {
                        let mut count = top - (kk - 1);
                        while count > 0 {
                            let step = count >> 1;
                            let mid = top - step;
                            if bc_ref.get(mid, kk) > i {
                                top = mid - 1;
                                count -= step + 1;
                            } else {
                                count = step;
                            }
                        }
                    }
                    top
                };
                out[(d + 1 - kk) as usize] = nn2;
                i -= bc_ref.get(nn2, kk);
            }
            out[d as usize] = i;
        });

        while cofacets.has_next_compressed(true, dist_ref, bc_ref) {
            let cofacet = cofacets.next_compressed(dist_ref, bc_ref);
            if cofacet.diameter <= threshold {
                working_coboundary.push(cofacet);
            }
        }
    }

    fn add_coboundary_compressed(
        &self,
        reduction_matrix: &CompressedSparseMatrix<DiameterEntry>,
        columns_to_reduce: &[DiameterIndex],
        index_column_to_add: usize,
        factor: CoefficientT,
        dim: IndexT,
        working_reduction_column: &mut BinaryHeap<DiameterEntry>,
        working_coboundary: &mut BinaryHeap<DiameterEntry>,
        cofacets: &mut CompressedCoboundaryEnumerator,
    ) {
        let col = DiameterEntry::from_diameter_index(
            columns_to_reduce[index_column_to_add],
            factor,
        );
        self.add_simplex_coboundary_compressed(
            col, dim, working_reduction_column, working_coboundary, cofacets,
        );
        for &simplex in reduction_matrix.subrange(index_column_to_add as IndexT) {
            let mut s = simplex;
            let new_coeff = s.coefficient() * factor % self.modulus;
            s.set_coefficient(new_coeff);
            self.add_simplex_coboundary_compressed(
                s, dim, working_reduction_column, working_coboundary, cofacets,
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Ripser<SparseDistanceMatrix> — same algorithm, sparse coboundary enumerator
// ---------------------------------------------------------------------------

impl Ripser<SparseDistanceMatrix> {
    fn get_edges_impl(&self) -> Vec<DiameterIndex> {
        let mut edges = Vec::new();
        for i in 0..self.n {
            for nd in &self.dist.neighbors[i as usize] {
                let j = nd.index;
                if i > j {
                    edges.push(DiameterIndex::new(nd.diameter, self.get_edge_index(i, j)));
                }
            }
        }
        edges
    }

    pub fn compute_barcodes(&mut self) {
        let mut simplices: Vec<DiameterIndex> = Vec::new();
        let mut columns_to_reduce: Vec<DiameterIndex> = Vec::new();

        let mut pivot_cofacet_enum = SparseCoboundaryEnumerator::new(self.modulus);
        let mut assemble_cofacet_enum = SparseCoboundaryEnumerator::new(self.modulus);
        let mut init_cofacet_enum = SparseCoboundaryEnumerator::new(self.modulus);
        let mut add_cofacet_enum = SparseCoboundaryEnumerator::new(self.modulus);
        let mut apparent_facet_enum = SparseCoboundaryEnumerator::new(self.modulus);
        let mut apparent_cofacet_enum = SparseCoboundaryEnumerator::new(self.modulus);

        self.compute_dim_0_pairs_sparse(
            &mut simplices,
            &mut columns_to_reduce,
            &mut pivot_cofacet_enum,
        );

        for dim in 1..=self.dim_max {
            let mut pivot_column_index: HashMap<IndexT, usize> = HashMap::new();
            pivot_column_index.reserve(columns_to_reduce.len());

            self.compute_pairs_sparse(
                &columns_to_reduce,
                &mut pivot_column_index,
                dim,
                &mut init_cofacet_enum,
                &mut add_cofacet_enum,
                &mut apparent_facet_enum,
            );

            if dim < self.dim_max {
                self.assemble_columns_to_reduce_sparse(
                    &mut simplices,
                    &mut columns_to_reduce,
                    &pivot_column_index,
                    dim + 1,
                    &mut assemble_cofacet_enum,
                    &mut apparent_facet_enum,
                    &mut apparent_cofacet_enum,
                );
            }
        }
    }

    fn get_zero_pivot_cofacet_sparse(
        &mut self,
        simplex: DiameterEntry,
        dim: IndexT,
        cofacets: &mut SparseCoboundaryEnumerator,
    ) -> DiameterEntry {
        let bc = &self.binomial_coeff as *const BinomialCoeffTable;
        let bc = unsafe { &*bc };
        let dist = &self.dist as *const SparseDistanceMatrix;
        let dist = unsafe { &*dist };
        cofacets.set_simplex(simplex, dim, dist, bc, |idx, d, nn, out| {
            out.resize((d + 1) as usize, 0);
            let mut i = idx;
            let mut nn2 = nn - 1;
            for kk in (2..=d + 1).rev() {
                nn2 = {
                    let mut top = nn2;
                    if bc.get(top, kk) > i {
                        let mut count = top - (kk - 1);
                        while count > 0 {
                            let step = count >> 1;
                            let mid = top - step;
                            if bc.get(mid, kk) > i {
                                top = mid - 1;
                                count -= step + 1;
                            } else {
                                count = step;
                            }
                        }
                    }
                    top
                };
                out[(d + 1 - kk) as usize] = nn2;
                i -= bc.get(nn2, kk);
            }
            out[d as usize] = i;
        });

        while cofacets.has_next_sparse(true, dist, bc) {
            let cofacet = cofacets.next_sparse(bc);
            if cofacet.diameter == simplex.diameter {
                return cofacet;
            }
        }
        DiameterEntry::invalid()
    }

    fn get_zero_apparent_facet_sparse(
        &mut self,
        simplex: DiameterEntry,
        dim: IndexT,
        cofacet_enum: &mut SparseCoboundaryEnumerator,
    ) -> DiameterEntry {
        // Use boundary (facet) enumeration inline
        let facet = self.get_zero_pivot_facet_sparse_inline(simplex, dim);
        if facet.index() == -1 {
            return DiameterEntry::invalid();
        }
        let cofacet = self.get_zero_pivot_cofacet_sparse(facet, dim - 1, cofacet_enum);
        if cofacet.index() == simplex.index() {
            facet
        } else {
            DiameterEntry::invalid()
        }
    }

    fn get_zero_pivot_facet_sparse_inline(&mut self, simplex: DiameterEntry, dim: IndexT) -> DiameterEntry {
        let idx = simplex.index();
        let diam = simplex.diameter;
        let modulus = self.modulus;
        let coeff = simplex.coefficient();

        let mut idx_below = idx;
        let mut idx_above: IndexT = 0;
        let mut j = self.n - 1;
        let mut k = dim;

        while k >= 0 {
            j = self.get_max_vertex(idx_below, k + 1, j);
            let face_index = idx_above - self.binomial_coeff.get(j, k + 1) + idx_below;
            let face_diam = self.compute_diameter(face_index, dim - 1);
            let face_coeff = if k & 1 != 0 {
                (modulus - 1 + coeff) % modulus
            } else {
                coeff % modulus
            };
            idx_below -= self.binomial_coeff.get(j, k + 1);
            idx_above += self.binomial_coeff.get(j, k);
            k -= 1;

            if face_diam == diam {
                return DiameterEntry::new(face_diam, face_index, face_coeff);
            }
        }
        DiameterEntry::invalid()
    }

    fn get_zero_apparent_cofacet_sparse(
        &mut self,
        simplex: DiameterEntry,
        dim: IndexT,
        cofacet_enum: &mut SparseCoboundaryEnumerator,
    ) -> DiameterEntry {
        let cofacet = self.get_zero_pivot_cofacet_sparse(simplex, dim, cofacet_enum);
        if cofacet.index() == -1 {
            return DiameterEntry::invalid();
        }
        let facet = self.get_zero_pivot_facet_sparse_inline(cofacet, dim + 1);
        if facet.index() == simplex.index() {
            cofacet
        } else {
            DiameterEntry::invalid()
        }
    }

    fn is_in_zero_apparent_pair_sparse(
        &mut self,
        simplex: DiameterEntry,
        dim: IndexT,
        cofacet_enum: &mut SparseCoboundaryEnumerator,
        facet_enum: &mut SparseCoboundaryEnumerator,
    ) -> bool {
        let cofacet = self.get_zero_apparent_cofacet_sparse(simplex, dim, cofacet_enum);
        if cofacet.index() != -1 { return true; }
        let facet = self.get_zero_apparent_facet_sparse(simplex, dim, facet_enum);
        facet.index() != -1
    }

    fn compute_dim_0_pairs_sparse(
        &mut self,
        simplices: &mut Vec<DiameterIndex>,
        columns_to_reduce: &mut Vec<DiameterIndex>,
        cofacet_enum: &mut SparseCoboundaryEnumerator,
    ) {
        #[cfg(feature = "print_persistence_pairs")]
        println!("persistence intervals in dim 0:");

        let mut dset = UnionFind::new(self.n);
        let mut edges = self.get_edges_impl();
        edges.sort();
        edges.reverse();

        for e in &edges {
            let idx = e.index;
            let v1 = self.get_max_vertex(idx, 2, self.n - 1);
            let v0 = idx - self.binomial_coeff.get(v1, 2);
            let u = dset.find(v0);
            let v = dset.find(v1);
            if u != v {
                #[cfg(feature = "print_persistence_pairs")]
                if e.diameter != 0.0 {
                    println!(" [0,{})", fmt_float(e.diameter));
                }
                dset.link(u, v);
            } else if self.dim_max > 0 {
                let de = DiameterEntry::from_diameter_index(*e, 1);
                let cofacet = self.get_zero_apparent_cofacet_sparse(de, 1, cofacet_enum);
                if cofacet.index() == -1 {
                    columns_to_reduce.push(*e);
                }
            }
        }

        if self.dim_max > 0 {
            columns_to_reduce.reverse();
        }

        // Populate simplices with all edges for assemble_columns_to_reduce (dim 2+).
        *simplices = edges;

        #[cfg(feature = "print_persistence_pairs")]
        for i in 0..self.n {
            if dset.find(i) == i {
                println!(" [0, )");
            }
        }
    }

    fn assemble_columns_to_reduce_sparse(
        &mut self,
        simplices: &mut Vec<DiameterIndex>,
        columns_to_reduce: &mut Vec<DiameterIndex>,
        pivot_column_index: &HashMap<IndexT, usize>,
        dim: IndexT,
        cofacets: &mut SparseCoboundaryEnumerator,
        apparent_facet_enum: &mut SparseCoboundaryEnumerator,
        apparent_cofacet_enum: &mut SparseCoboundaryEnumerator,
    ) {
        columns_to_reduce.clear();
        let mut next_simplices: Vec<DiameterIndex> = Vec::new();

        let threshold = self.threshold;
        let dim_max = self.dim_max;

        for &simplex in simplices.iter() {
            let de = DiameterEntry::from_diameter_index(simplex, 1);
            let bc = &self.binomial_coeff as *const BinomialCoeffTable;
            let bc = unsafe { &*bc };
            let dist = &self.dist as *const SparseDistanceMatrix;
            let dist = unsafe { &*dist };

            cofacets.set_simplex(de, dim - 1, dist, bc, |idx, d, nn, out| {
                out.resize((d + 1) as usize, 0);
                let mut i = idx;
                let mut nn2 = nn - 1;
                for kk in (2..=d + 1).rev() {
                    nn2 = {
                        let mut top = nn2;
                        if bc.get(top, kk) > i {
                            let mut count = top - (kk - 1);
                            while count > 0 {
                                let step = count >> 1;
                                let mid = top - step;
                                if bc.get(mid, kk) > i {
                                    top = mid - 1;
                                    count -= step + 1;
                                } else {
                                    count = step;
                                }
                            }
                        }
                        top
                    };
                    out[(d + 1 - kk) as usize] = nn2;
                    i -= bc.get(nn2, kk);
                }
                out[d as usize] = i;
            });

            while cofacets.has_next_sparse(false, dist, bc) {
                let cofacet = cofacets.next_sparse(bc);
                if cofacet.diameter <= threshold {
                    if dim < dim_max {
                        next_simplices.push(DiameterIndex::new(cofacet.diameter, cofacet.index()));
                    }
                    let in_pair = self.is_in_zero_apparent_pair_sparse(
                        cofacet, dim, apparent_cofacet_enum, apparent_facet_enum,
                    );
                    let in_pivot = pivot_column_index.contains_key(&cofacet.index());
                    if !in_pair && !in_pivot {
                        columns_to_reduce
                            .push(DiameterIndex::new(cofacet.diameter, cofacet.index()));
                    }
                }
            }
        }

        if dim < dim_max {
            *simplices = next_simplices;
        }

        columns_to_reduce.sort();
        columns_to_reduce.reverse();
    }

    fn pop_pivot(
        column: &mut BinaryHeap<DiameterEntry>,
        #[cfg(feature = "use_coefficients")] modulus: CoefficientT,
    ) -> DiameterEntry {
        #[cfg(not(feature = "use_coefficients"))]
        {
            while let Some(pivot) = column.pop() {
                match column.peek() {
                    Some(next) if next.index() == pivot.index() => { column.pop(); }
                    _ => return pivot,
                }
            }
            DiameterEntry::invalid()
        }
        #[cfg(feature = "use_coefficients")]
        {
            let mut pivot = DiameterEntry::invalid();
            while let Some(top) = column.pop() {
                if pivot.coefficient() == 0 {
                    pivot = top;
                } else if pivot.index() != top.index() {
                    column.push(top);
                    return pivot;
                } else {
                    let c = (pivot.coefficient() + top.coefficient()) % modulus;
                    pivot.set_coefficient(c);
                }
            }
            if pivot.coefficient() == 0 { DiameterEntry::invalid() } else { pivot }
        }
    }

    fn get_pivot(column: &mut BinaryHeap<DiameterEntry>,
                 #[cfg(feature = "use_coefficients")] modulus: CoefficientT) -> DiameterEntry {
        #[cfg(feature = "use_coefficients")]
        let r = Self::pop_pivot(column, modulus);
        #[cfg(not(feature = "use_coefficients"))]
        let r = Self::pop_pivot(column);
        if r.index() != -1 { column.push(r); }
        r
    }

    fn compute_pairs_sparse(
        &mut self,
        columns_to_reduce: &[DiameterIndex],
        pivot_column_index: &mut HashMap<IndexT, usize>,
        dim: IndexT,
        init_cofacets: &mut SparseCoboundaryEnumerator,
        add_cofacets: &mut SparseCoboundaryEnumerator,
        apparent_facet_enum: &mut SparseCoboundaryEnumerator,
    ) {
        #[cfg(feature = "print_persistence_pairs")]
        println!("persistence intervals in dim {}:", dim);

        let mut reduction_matrix: CompressedSparseMatrix<DiameterEntry> =
            CompressedSparseMatrix::new();

        for (index_column_to_reduce, &col_di) in columns_to_reduce.iter().enumerate() {
            let column_to_reduce = DiameterEntry::from_diameter_index(col_di, 1);
            let diameter = col_di.diameter;

            reduction_matrix.append_column();

            let mut working_reduction_column: BinaryHeap<DiameterEntry> = BinaryHeap::new();
            let mut working_coboundary: BinaryHeap<DiameterEntry> = BinaryHeap::new();

            let mut pivot = self.init_coboundary_and_get_pivot_sparse(
                column_to_reduce,
                &mut working_coboundary,
                dim,
                pivot_column_index,
                init_cofacets,
                apparent_facet_enum,
            );

            loop {
                if pivot.index() != -1 {
                    if let Some(&index_column_to_add) = pivot_column_index.get(&pivot.index()) {
                        let factor = 1u16; // Z/2
                        self.add_coboundary_sparse(
                            &reduction_matrix,
                            columns_to_reduce,
                            index_column_to_add,
                            factor,
                            dim,
                            &mut working_reduction_column,
                            &mut working_coboundary,
                            add_cofacets,
                        );
                        #[cfg(feature = "use_coefficients")]
                        { pivot = Self::get_pivot(&mut working_coboundary, self.modulus); }
                        #[cfg(not(feature = "use_coefficients"))]
                        { pivot = Self::get_pivot(&mut working_coboundary); }
                    } else {
                        let e = self.get_zero_apparent_facet_sparse(pivot, dim + 1, apparent_facet_enum);
                        if e.index() != -1 {
                            let mut e2 = e;
                            e2.set_coefficient((self.modulus - e2.coefficient()) % self.modulus);
                            self.add_simplex_coboundary_sparse(
                                e2, dim,
                                &mut working_reduction_column,
                                &mut working_coboundary,
                                add_cofacets,
                            );
                            #[cfg(feature = "use_coefficients")]
                            { pivot = Self::get_pivot(&mut working_coboundary, self.modulus); }
                            #[cfg(not(feature = "use_coefficients"))]
                            { pivot = Self::get_pivot(&mut working_coboundary); }
                        } else {
                            let death = pivot.diameter;
                            #[cfg(feature = "print_persistence_pairs")]
                            if death > diameter * self.ratio {
                                println!(" [{},{})", fmt_float(diameter), fmt_float(death));
                            }
                            pivot_column_index.insert(pivot.index(), index_column_to_reduce);
                            loop {
                                #[cfg(feature = "use_coefficients")]
                                let e = Self::pop_pivot(&mut working_reduction_column, self.modulus);
                                #[cfg(not(feature = "use_coefficients"))]
                                let e = Self::pop_pivot(&mut working_reduction_column);
                                if e.index() == -1 { break; }
                                reduction_matrix.push_back(e);
                            }
                            break;
                        }
                    }
                } else {
                    #[cfg(feature = "print_persistence_pairs")]
                    println!(" [{}, )", fmt_float(diameter));
                    break;
                }
            }
        }
    }

    fn init_coboundary_and_get_pivot_sparse(
        &mut self,
        simplex: DiameterEntry,
        working_coboundary: &mut BinaryHeap<DiameterEntry>,
        dim: IndexT,
        pivot_column_index: &HashMap<IndexT, usize>,
        cofacets: &mut SparseCoboundaryEnumerator,
        apparent_facet_enum: &mut SparseCoboundaryEnumerator,
    ) -> DiameterEntry {
        let mut check_emergent = true;
        self.cofacet_entries.clear();

        let bc = &self.binomial_coeff as *const BinomialCoeffTable;
        let bc = unsafe { &*bc };
        let dist = &self.dist as *const SparseDistanceMatrix;
        let dist = unsafe { &*dist };
        let threshold = self.threshold;

        cofacets.set_simplex(simplex, dim, dist, bc, |idx, d, nn, out| {
            out.resize((d + 1) as usize, 0);
            let mut i = idx;
            let mut nn2 = nn - 1;
            for kk in (2..=d + 1).rev() {
                nn2 = {
                    let mut top = nn2;
                    if bc.get(top, kk) > i {
                        let mut count = top - (kk - 1);
                        while count > 0 {
                            let step = count >> 1;
                            let mid = top - step;
                            if bc.get(mid, kk) > i {
                                top = mid - 1;
                                count -= step + 1;
                            } else {
                                count = step;
                            }
                        }
                    }
                    top
                };
                out[(d + 1 - kk) as usize] = nn2;
                i -= bc.get(nn2, kk);
            }
            out[d as usize] = i;
        });

        while cofacets.has_next_sparse(true, dist, bc) {
            let cofacet = cofacets.next_sparse(bc);
            if cofacet.diameter <= threshold {
                self.cofacet_entries.push(cofacet);
                if check_emergent && cofacet.diameter == simplex.diameter {
                    let not_in_pivot = !pivot_column_index.contains_key(&cofacet.index());
                    let af = self.get_zero_apparent_facet_sparse(cofacet, dim + 1, apparent_facet_enum);
                    if not_in_pivot && af.index() == -1 {
                        return cofacet;
                    }
                    check_emergent = false;
                }
            }
        }

        for &c in &self.cofacet_entries {
            working_coboundary.push(c);
        }
        #[cfg(feature = "use_coefficients")]
        return Self::get_pivot(working_coboundary, self.modulus);
        #[cfg(not(feature = "use_coefficients"))]
        return Self::get_pivot(working_coboundary);
    }

    fn add_simplex_coboundary_sparse(
        &self,
        simplex: DiameterEntry,
        dim: IndexT,
        working_reduction_column: &mut BinaryHeap<DiameterEntry>,
        working_coboundary: &mut BinaryHeap<DiameterEntry>,
        cofacets: &mut SparseCoboundaryEnumerator,
    ) {
        working_reduction_column.push(simplex);
        let bc = &self.binomial_coeff as *const BinomialCoeffTable;
        let bc = unsafe { &*bc };
        let dist = &self.dist as *const SparseDistanceMatrix;
        let dist = unsafe { &*dist };
        let threshold = self.threshold;

        cofacets.set_simplex(simplex, dim, dist, bc, |idx, d, nn, out| {
            out.resize((d + 1) as usize, 0);
            let mut i = idx;
            let mut nn2 = nn - 1;
            for kk in (2..=d + 1).rev() {
                nn2 = {
                    let mut top = nn2;
                    if bc.get(top, kk) > i {
                        let mut count = top - (kk - 1);
                        while count > 0 {
                            let step = count >> 1;
                            let mid = top - step;
                            if bc.get(mid, kk) > i {
                                top = mid - 1;
                                count -= step + 1;
                            } else {
                                count = step;
                            }
                        }
                    }
                    top
                };
                out[(d + 1 - kk) as usize] = nn2;
                i -= bc.get(nn2, kk);
            }
            out[d as usize] = i;
        });

        while cofacets.has_next_sparse(true, dist, bc) {
            let cofacet = cofacets.next_sparse(bc);
            if cofacet.diameter <= threshold {
                working_coboundary.push(cofacet);
            }
        }
    }

    fn add_coboundary_sparse(
        &self,
        reduction_matrix: &CompressedSparseMatrix<DiameterEntry>,
        columns_to_reduce: &[DiameterIndex],
        index_column_to_add: usize,
        factor: CoefficientT,
        dim: IndexT,
        working_reduction_column: &mut BinaryHeap<DiameterEntry>,
        working_coboundary: &mut BinaryHeap<DiameterEntry>,
        cofacets: &mut SparseCoboundaryEnumerator,
    ) {
        let col = DiameterEntry::from_diameter_index(columns_to_reduce[index_column_to_add], factor);
        self.add_simplex_coboundary_sparse(
            col, dim, working_reduction_column, working_coboundary, cofacets,
        );
        for &simplex in reduction_matrix.subrange(index_column_to_add as IndexT) {
            let mut s = simplex;
            let c = s.coefficient() * factor % self.modulus;
            s.set_coefficient(c);
            self.add_simplex_coboundary_sparse(
                s, dim, working_reduction_column, working_coboundary, cofacets,
            );
        }
    }
}
