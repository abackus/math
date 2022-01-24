//! An undirected weighted graph from which we can create a graph Laplacian

use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};
extern crate nalgebra as na;
use na::{DMatrix, DVector};

/// A single vertex in a graph
#[derive(Copy, Clone, Eq, PartialEq, Hash, PartialOrd, Debug)]
pub struct Vertex {
    index: usize,
}

/// The set of vertices of a graph, which is a lazily evaluted iterator
pub struct VertexSet {
    curr: usize,
    cardinality: usize,
}

/// An unweighted, undirected edge between two vertices
#[derive(Eq, Debug)]
pub struct Edge {
    left: Vertex,
    right: Vertex,
}

/// A graph consists of a set of vertices and a set of edges
/// Every edge must link vertices in the graph
#[derive(PartialEq)]
pub struct Graph {
    edges: HashMap<Edge, f64>,
    vertex_count: usize,
}

/// A function on a graph takes each vertex and returns a float
pub struct Function<'a> {
    pub func: Box<dyn Fn(Vertex) -> f64 + 'a>,
    pub domain: &'a Graph,
}

/// A measure is identified with its Radon-Nikodym derivative, but is assumed *nonnegative*
pub type Measure<'a> = Function<'a>;

impl Edge {
    /// Create an unweighted, undirected edge
    /// Panicks if the edge is from a vertex to itself
    fn new(left: Vertex, right: Vertex) -> Edge {
        if left == right {
            panic!("An edge cannot point from an edge to itself");
        }
        Edge {
            left,
            right,
        }
    }
}

/// The unordered pairs (x, y) and (y, x) are equivalent
impl PartialEq for Edge {
    fn eq(&self, other: &Self) -> bool {
        (self.left == other.left && self.right == other.right) || (self.left == other.right && self.right == other.left)
    }
}

/// The hash of the pairs (x, y) and (y, x) are the same
impl Hash for Edge {
    fn hash<H>(&self, state: &mut H) where H: Hasher {
        if self.left >= self.right {
            self.left.hash(state);
            self.right.hash(state);
        }
        else {
            self.right.hash(state);
            self.left.hash(state);
        }
    }
}

impl Graph {
    /// Creates an empty Graph
    pub fn new() -> Graph {
        Graph {
            edges: HashMap::new(),
            vertex_count: 0,
        }
    }

    /// Create an isolated vertex, add it to the graph, and return it
    pub fn add_vertex(&mut self) -> Vertex {
        let vertex = {
            Vertex {
                index: self.vertex_count,
            }
        };
        self.vertex_count = self.vertex_count + 1;
        vertex
    }

    /// Add m isolated vertices, and return a vector of them
    /// This is faster than calling add_vertex m times
    pub fn add_vertices(&mut self, m: usize) -> Vec<Vertex> {
        let mut v = Vec::with_capacity(m);
        for n in 0..m {
            v.push(
                Vertex {
                    index: self.vertex_count + n,
                }
            );
        }
        self.vertex_count = self.vertex_count + m;
        v
    }

    /// Check if a graph already has a vertex
    pub fn has_vertex(&self, vertex: Vertex) -> bool {
        self.vertex_count > vertex.index
    }

    /// Returns the set of all vertices
    pub fn vertex_set(&self) -> VertexSet {
        VertexSet {
            curr: 0,
            cardinality: self.vertex_count
        }
    }

    /// Return the number of vertices in a graph
    pub fn len(&self) -> usize {
        self.vertex_count
    }

    /// Adds an edge between two vertices already in the graph
    /// Panicks if we refer to a vertex that does not exist
    /// If the edge is already in the graph, it overwrites the current weight
    pub fn add_edge(&mut self, edge: Edge, weight: f64) {
        if !self.has_vertex(edge.left) || !self.has_vertex(edge.right) {
            panic!("We cannot create an edge between vertices that do not exist");
        } else {
            self.edges.insert(edge, weight);
        }
    }

    /// Returns the graph Laplacian matrix
    /// TODO: Make this a sparse matrix?
    pub fn laplacian(&self) -> DMatrix<f64> {
        let d = self.vertex_count;
        let mut m = DMatrix::<f64>::zeros(d, d);
        for (edge, weight) in self.edges.iter() {
            let i = edge.left.index;
            let j = edge.right.index;
            m[(i, j)] = -*weight;
            m[(j, i)] = -*weight;
            m[(i, i)] = m[(i, i)] + *weight;
            m[(j, j)] = m[(j, j)] + *weight;
        }
        m
    }
    
    /// Returns the uniform probability measure
    /// If we add more vertices, this measure becomes incorrectly normalized
    pub fn uniform_measure(&self) -> Measure {
        let mass = 1.0/(self.vertex_count as f64);
        Function {
            domain: self,
            func: Box::new(move |_x| {
                mass
            })
        }
    }
}

impl Iterator for VertexSet {
    type Item = Vertex;

    fn next(&mut self) -> Option<Self::Item> {
        if self.curr < self.cardinality {
            self.curr = self.curr + 1;
            Some(Vertex {
                index: self.curr - 1,
            })
        } else {
            None
        }
    }
}

impl Function<'_> {
    /// Converts a function to a vector indexed by the vertices
    fn to_vector(&self) -> DVector<f64> {
        let mut v = Vec::with_capacity(self.domain.vertex_count);
        for x in self.domain.vertex_set() {
            v.push(self.eval(x));
        }
        return DVector::from_vec(v);
    }

    /// Converts a vector indexed by the vertices to a function
    /// This consumes the input vector
    fn from_vector(g: &Graph, v: DVector<f64>) -> Function {
        Function {
            domain: g,
            func: Box::new(move |x| {
                v[x.index]
            })
        }
    }

    /// Evaluate the function at a vertex
    pub fn eval(&self, x: Vertex) -> f64 {
        if !self.domain.has_vertex(x) {
            panic!("Tried to call a function on a vertex not in its domain");
        }
        (self.func)(x)
    }
    
    /// Returns the support of the function
    pub fn support(&self) -> HashSet<Vertex> {
        let mut s = HashSet::new();
        for x in self.domain.vertex_set() {
            if !(self.eval(x) == 0.0) {
                s.insert(x);
            }
        }
        s
    }

    /// Returns the positive-definite Laplacian of the function
    pub fn laplace(&self) -> Function {
        Function::from_vector(self.domain, self.domain.laplacian() * self.to_vector())
    }
}

impl Measure<'_> {
    /// Returns the Dirac mass at a vertex
    pub fn dirac(domain: &Graph, vertex: Vertex) -> Measure {
        let copy = vertex;
        Measure {
            domain,
            func: Box::new(move |x| {
                if x == copy {
                    1.0
                } else {
                    0.0
                }
            })
        }
    }
    
    /// Compute the integral of a function with respect to self
    pub fn integral(&self, f: &Function) -> f64 {
        if !(self.domain == f.domain) {
            panic!("Tried to integrate function vs measure with different domains")
        }
        let mut s = 0.0;
        for x in self.support() {
            s = s + f.eval(x) * self.eval(x);
        }
        s
    }
    
    /// Returns the norm ||f||_{L^p(mu)}
    /// This method is NOT numerically stable for p large but finite
    pub fn lebesgue_norm(&self, f: &Function, p: f64) -> f64 {
        if p == f64::INFINITY {
            if !(self.domain == f.domain) {
                panic!("Tried to integrate function vs measure with different domains")
            }
            let mut norm = 0.0;
            for x in self.support() {
                let y = f.eval(x).abs();
                if y > norm {
                    norm = y;
                }
            }
            norm
        } else if p > 0.0 {
            self.integral(&Function {
                domain: self.domain,
                func: Box::new(move |x| {
                    f.eval(x).abs().powf(p)
                })
            }).powf(1.0/p)
        } else if p == 0.0 {
            if !(self.domain == f.domain) {
                panic!("Tried to integrate function vs measure with different domains")
            }
            let mut norm = 0.0;
            for x in self.support() {
                if !(f.eval(x) == 0.0) {
                    norm = norm + 1.0
                }
            }
            norm
        } else {
            panic!("p must be nonnegative")
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Try to make a Boolean algebra with two elements.
    #[test]
    fn make_bool_alg2() {
        let mut g = Graph::new();
        let v1 = Vertex { index: 0 };
        let v2 = Vertex { index: 1 };
        assert_eq!(g.add_vertex(), v1);
        assert_eq!(g.add_vertex(), v2);
        g.add_edge(Edge::new(v1, v2), 1.0);
        let mut m = DMatrix::<f64>::zeros(2, 2);
        m[(0, 1)] = -1.0;
        m[(1, 0)] = -1.0;
        m[(1, 1)] = 1.0;
        m[(0, 0)] = 1.0;
        assert_eq!(g.laplacian(), m);
    }

    /// Try to make a Dirac mass on a graph with 100 vertices
    #[test]
    fn make_dirac() {
        let mut g = Graph::new();
        let v = g.add_vertex();
        let w = g.add_vertex();
        assert_ne!(v, w);
        g.add_vertices(98);
        assert_eq!(g.len(), 100);
        let f = Measure::dirac(&g, v);
        assert_eq!(f.eval(v), 1.0);
        assert_eq!(f.eval(w), 0.0);
    } 

    /// The support of a uniform measure should be the whole graph
    #[test]
    fn test_support() {
        let mut g = Graph::new();
        g.add_vertices(100);
        let supp = &g.uniform_measure().support();
        assert_eq!(supp.len(), 100);
    }
}