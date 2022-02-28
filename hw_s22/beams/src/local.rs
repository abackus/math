//! This module deals with procedures that work locally,
//! thus in a single finite element, node, or dof
//! The main purpose is to construct element stiffness matrices and load vectors.

use nalgebra::{Matrix4, Vector4};

/// A degree of freedom is just determined by an index
#[derive(Debug, Copy, Clone, Eq, PartialEq, Hash)]
pub struct Dof {
    index: usize,
}

impl Dof {
    /// Create a dof with some index
    pub fn new(index: usize) -> Dof {
        Dof {
            index,
        }
    }

    /// Return the index of a dof 
    pub fn index(self) -> usize {
        self.index
    }
}

/// We assume that every node has Dofs for the 0th and 1st derivative
// In future homeworks I should replace this struct with 
// pub struct Node<dim, sobolev> {
//    pub coords: [f64; dim],
//    pub dofs: [Dof; sobolev],
// }
// but implementing that will take a lot of work and I'm lazy

#[derive(Debug, Copy, Clone)]
pub struct Node {
    pub coord: f64,
    pub dof0: Dof, // 0th derivative
    pub dof1: Dof, // 1st derivative
}

/// A single finite element.
/// We index the dofs as follows:
/// * 0 - left node, 0th derivative
/// * 1 - left node, 1st derivative
/// * 2 - right node, 0th derivative
/// * 3 - right node, 1st derivative
#[derive(Debug, Copy, Clone)]
pub struct FiniteElement {
    pub left: Node,
    pub right: Node,
}

impl FiniteElement {
    /// Create an element from two nodes 
    pub fn new(left: Node, right: Node) -> FiniteElement {
        FiniteElement {
            left,
            right,
        }
    }

    /// Return the dof of index 0 <= i <= 3
    pub fn eldof(self, i: usize) -> usize {
        match i {
            0 => self.left.dof0.index,
            1 => self.left.dof1.index,
            2 => self.right.dof0.index,
            3 => self.right.dof1.index,
            _ => panic!("Invalid eldof input")
        }
    }

    /// Return the length of an element
    pub fn len(self) -> f64 {
        self.right.coord - self.left.coord
    }

    /// Return the element stiffness matrix for the beam equation
    /// We assume locally constant flexural rigidity `flex`
    pub fn beam_stiffness(self, flex: f64) -> Matrix4<f64> {
        let h = self.len();
        flex * f64::powf(h, -2.0) * Matrix4::new(
            12.0 / h, 6.0, -12.0 / h, 6.0,
            6.0, 4.0 * h, -6.0, 2.0 * h,
            -12.0 / h, -6.0, 12.0 / h, -6.0,
            6.0, 2.0 * h, -6.0, 4.0 * h
        )
    }
    
    /// Return the element load vector for the beam equation
    /// We assume locally constant force `force`
    pub fn beam_load(self, force: f64) -> Vector4<f64> {
        let h = self.len();
        h * force * Vector4::new(
            0.5,
            h / 12.0,
            0.5,
            -h / 12.0,
        )
    }
}