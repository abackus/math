//! This module deals with procedures that work globally,
//! thus in all finite elements at once or over the whole domain
//! The main purpose is to subassemble the stiffness matrix and load vector

use crate::local::{Dof, Node, FiniteElement};
use nalgebra::{DVector, DMatrix};
use std::collections::VecDeque;

/// The domain of the problem consists of all finite elements
#[derive(Debug)]
pub struct Domain {
    pub elts: Vec<FiniteElement>,
    dofcounter: usize,
}

impl Domain {
    /// Create a domain starting at 0 with specified distances between nodes
    pub fn new(spacing: Vec<f64>) -> Domain {
        let mut domain = Domain {
            elts: Vec::with_capacity(spacing.len() + 1),
            dofcounter: 0,
        };
        let mut coord = 0.0;
        let mut curr_nodes = VecDeque::with_capacity(2);
        curr_nodes.push_front(domain.create_node(coord));
        for h in spacing {
            coord += h;
            curr_nodes.push_front(
                domain.create_node(coord)
            );
            domain.create_elt(
                curr_nodes.pop_back().unwrap(),
                *curr_nodes.back().unwrap()
            );
        }
        domain
    }

    /// Create a domain starting at 0 with evenly spaced nodes
    // This isn't the most efficient way to do this (that would be with an Iterator)
    // but I'm lazy
    pub fn evenly_spaced(nodecount: usize, spacing: f64) -> Domain {
        let mut spacing_vector = Vec::with_capacity(nodecount - 1);
        for _ in 1..nodecount {
            spacing_vector.push(spacing);
        }
        Domain::new(spacing_vector)
    }

    /// Create a new node in the domain
    pub fn create_node(&mut self, coord: f64) -> Node {
        let node = Node {
            coord,
            dof0: Dof::new(self.dofcounter),
            dof1: Dof::new(self.dofcounter + 1),
        };
        self.dofcounter += 2;
        node
    }

    /// Create a new element in the domain out of two nodes
    pub fn create_elt(&mut self, left: Node, right: Node) -> FiniteElement {
        let elt = FiniteElement::new(left, right);
        self.elts.push(elt);
        elt
    }

    /// Return the stiffness matrix for the beam equation
    /// We assume constant flexural rigidity `flex`
    pub fn beam_stiffness(&self, flex: f64) -> DMatrix<f64> {
        let mut stiffness = DMatrix::zeros(self.dofcounter, self.dofcounter);
        for elt in &self.elts {
            for i in 0..4 {
                for j in 0..4 {
                    stiffness[(elt.eldof(i), elt.eldof(j))] = elt.beam_stiffness(flex)[(i, j)];
                }
            }
        }
        stiffness
    }

    /// Return the load vector for the beam equation
    /// We assume constant force
    pub fn beam_load(&self, force: f64) -> DVector<f64> {
        let mut load = DVector::zeros(self.dofcounter);
        for elt in &self.elts {
            for i in 0..4 {
                load[elt.eldof(i)] = elt.beam_load(force)[i];
            }
        }
        load
    }
}