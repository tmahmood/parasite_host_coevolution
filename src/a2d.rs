use std::collections::hash_map::Iter;
use std::collections::HashMap;

pub struct A2D<T: Clone> {
    arr: HashMap<(usize, usize), T>,
    column_len: usize,
}

impl<T: Clone> A2D<T> {
    pub fn new(column_len: usize) -> Self {
        A2D {
            arr: HashMap::new(),
            column_len,
        }
    }

    pub fn insert(&mut self, col: usize, row: usize, item: T) {
        self.arr.insert((row, col), item);
    }

    pub fn get(&self, col: usize, row: usize) -> Option<&T> {
        self.arr.get(&(row, col))
    }

    pub fn remove(&mut self, col: usize, row: usize) -> Option<T> {
        self.arr.remove(&(col, row))
    }

    pub fn iter(&self) -> Iter<'_, (usize, usize), T> {
        self.arr.iter()
    }

    pub fn len(&self) -> usize {
        self.arr.len()
    }
}

impl<T: Clone> From<&Vec<T>> for A2D<T> {
    fn from(vec: &Vec<T>) -> Self {
        let mut a: A2D<T> = A2D::new(vec.len());
        for (i, v) in vec.iter().enumerate() {
            a.insert(i, 0, v.clone());
        }
        a
    }
}
