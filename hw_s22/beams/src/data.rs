//! This module handles point masses and boundary data
//! It can apply point loads and essential boundary data to the stiffness matrix and load vector

use crate::local::Dof;
use std::collections::{HashMap, HashSet};
use nalgebra::{DVector, DMatrix};

/// A table of point masses
#[derive(Debug)]
pub struct PointMasses {
    dofs: HashMap<Dof, f64>,
}

impl PointMasses {
    /// Create an empty point mass table 
    pub fn new() -> PointMasses {
        PointMasses {
            dofs: HashMap::new()
        }
    }

    /// Insert a new point mass
    pub fn insert(&mut self, dof: Dof, mass: f64) {
        self.dofs.insert(dof, mass);
    }

    /// Mutate a load vector by applying self
    /// This consumes self
    pub fn apply(self, load: &mut DVector<f64>) {
        for (dof, mass) in self.dofs {
            load[dof.index()] += mass;
        }
    }
}

// A set of homogeneous boundary data dofs
#[derive(Debug)]
pub struct EssentialBoundaryData {
    dofs: HashSet<Dof>,
}

impl EssentialBoundaryData {
    /// Create essential boundary data from a list of dofs
    pub fn new(dofs: HashSet<Dof>) -> EssentialBoundaryData {
        EssentialBoundaryData {
            dofs,
        }
    }

    /// Mutate a load vector and stiffness matrix by applying self
    /// This consumes self 
    pub fn apply(self, load: &mut DVector<f64>, stiffness: &mut DMatrix<f64>) {
        for dof in self.dofs {
            let i = dof.index();
            load[i] = 0.0;
            for j in 0..load.len() {
                if i == j {
                    stiffness[(i, i)] = 1.0;
                } else {
                    stiffness[(i, j)] = 0.0;
                    stiffness[(j, i)] = 0.0;
                }
            }
        }
    }
}