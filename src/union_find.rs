use crate::types::IndexT;

pub struct UnionFind {
    parent: Vec<IndexT>,
    rank: Vec<u8>,
}

impl UnionFind {
    pub fn new(n: IndexT) -> Self {
        let n = n as usize;
        let parent = (0..n as IndexT).collect();
        UnionFind {
            parent,
            rank: vec![0; n],
        }
    }

    pub fn find(&mut self, mut x: IndexT) -> IndexT {
        // Path halving
        let mut y = x;
        let mut z;
        loop {
            z = self.parent[y as usize];
            if z == y {
                break;
            }
            y = z;
        }
        loop {
            z = self.parent[x as usize];
            if z == y {
                break;
            }
            self.parent[x as usize] = y;
            x = z;
        }
        z
    }

    pub fn link(&mut self, x: IndexT, y: IndexT) {
        let x = self.find(x);
        let y = self.find(y);
        if x == y {
            return;
        }
        if self.rank[x as usize] > self.rank[y as usize] {
            self.parent[y as usize] = x;
        } else {
            self.parent[x as usize] = y;
            if self.rank[x as usize] == self.rank[y as usize] {
                self.rank[y as usize] += 1;
            }
        }
    }
}
