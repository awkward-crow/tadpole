use crate::types::IndexT;

/// Compressed sparse column matrix — mirrors C++ `compressed_sparse_matrix<T>`.
/// `bounds[i]` is the exclusive end index into `entries` for column i.
pub struct CompressedSparseMatrix<T> {
    bounds: Vec<usize>,
    entries: Vec<T>,
}

impl<T: Clone> CompressedSparseMatrix<T> {
    pub fn new() -> Self {
        CompressedSparseMatrix {
            bounds: Vec::new(),
            entries: Vec::new(),
        }
    }

    pub fn size(&self) -> usize {
        self.bounds.len()
    }

    pub fn append_column(&mut self) {
        self.bounds.push(self.entries.len());
    }

    pub fn push_back(&mut self, e: T) {
        debug_assert!(self.size() > 0);
        self.entries.push(e);
        *self.bounds.last_mut().unwrap() += 1;
    }

    pub fn subrange(&self, index: IndexT) -> &[T] {
        let index = index as usize;
        let start = if index == 0 { 0 } else { self.bounds[index - 1] };
        let end = self.bounds[index];
        &self.entries[start..end]
    }
}
