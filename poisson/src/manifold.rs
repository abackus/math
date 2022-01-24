//! Discrete manifolds 

use crate::graph::{Graph, Vertex};

/// A manifold is a graph along with an atlas and a dimension
/// An atlas is a vector of charts that cover the manifold
pub struct Manifold {
    graph: Graph,
    atlas: Vec<Chart>,
    dim: usize,
}

/// Each chart is a closure that maps a vector of integral coordinates to a vertex
/// We assume that each vector has length dim
pub type Chart = fn(Vec<i32>) -> Vertex;

/// Two manifolds are equivalent if they have the same dimension and same graph
impl PartialEq for Manifold {
    fn eq(&self, other: &Self) -> bool {
        self.dim == other.dim && self.graph == other.graph
    }
}

impl Manifold {
    /// Returns a discrete flat torus
    /// The returned atlas consists of one chart
    /// The chart maps (x_1, ..., x_dim) to (x_1/scale, ..., x_dim/scale)
    /// # Arguments
    /// scale - each edge has weight 1/scale
    /// lens - the ith side has lens[i] * scale nodes, the dimension is lens.len()
    pub fn flat_torus(scale: usize, lens: Vec<usize>) -> Manifold {
        let dim = lens.len();
        let mut graph = Graph::new();
        let mut atlas = Vec::new();
        // Remove the following
        graph.add_vertex();
        println!("{}", scale);
        // TODO
        Manifold {
            graph,
            atlas,
            dim,
        }
    }
}