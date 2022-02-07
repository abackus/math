use std::collections::HashMap;
use nalgebra::DVector;

/// A Point of dimension d is a concrete point (with coordinates) in R^d
#[derive(Copy, Clone)]
pub struct Point {
    pub coords: DVector<f64>,
}

/// A Vertex is an abstract (coordinate-free) point
/// It's better to check if two Vertices are equal than two Points
/// because of floating-point issues
#[derive(Copy, Clone, PartialEq, Eq, Hash)]
pub struct Vertex {
    index: usize,
}

/// A Simplex of dimension d consists of d+1 vertices
/// Simplices are always assumed nonempty!
#[derive(Copy, Clone, PartialEq)]
pub struct Simplex {
    pub vertices: Vec<Vertex>,
}

/// A Complex consists of an ordered set of Simplices, along with an assignment
/// from each Vertex in each Simplex to a Point
/// We always assume that the Points associated to each Simplex are affinely independent
pub struct SimplicialComplex {
    pub simplices: Vec<Simplex>,
    pub vertices: HashMap<Vertex, Point>,
    max_index: usize,
}

impl Point {
    /// Checks if a vector of d + 1 points in R^d is affine-dependent
    pub fn affine_dependent(points: Vec<Point>) -> bool {
        // TODO Implement me
        false
    }

    /// Returns the dimension of the Point
    pub fn dimension(&self) -> usize {
        self.coords.len()
    }

    /// Returns the origin
    pub fn origin(dimension: usize) -> Point {
        Point {
            coords: DVector::<f64>::zeros(dimension),
        }
    }
}

impl Simplex {
    /// Return the simplex with a vertex attached
    pub fn append_vertex(self, vertex: Vertex) -> Simplex { 
        let mut vertices = self.vertices.clone();
        vertices.push(vertex);
        Simplex {
            vertices,
        }
    }

    /// Returns all faces of the simplex
    pub fn boundary(self) -> Vec<Simplex> {
        let d = self.dimension();
        if d <= 0 {
            return vec![];
        }
        let d = d as usize;
        let mut v = Vec::with_capacity(d);
        for i in 0..d {
            // The ith face has the ith vertex deleted
            let mut vertices = Vec::with_capacity(d);
            for j in 0..d + 1 {
                if j != i {
                    vertices.push(self.vertices[j]);
                }
            }
            let face = Simplex {
                vertices,
            };
            v.push(face);
        }
        v
    }

    /// Returns the dimension 
    pub fn dimension(self) -> usize {
        let card = self.vertices.len();
        if card == 0 {
            panic!("Empty simplex");
        }
        card - 1
    }
}

impl SimplicialComplex {
    /// Creates a simplicial complex which consists of a single simplex
    pub fn new(points: Vec<Point>) -> SimplicialComplex {
        let d = points.len() - 1;
        let mut vertices: HashMap<Vertex, Point> = HashMap::new();
        let mut simplex: Vec<Vertex> = Vec::with_capacity(d + 1);
        let mut i = 0;
        if Point::affine_dependent(&points) {
            panic!("Tried to make a degenerate simplex");
        }
        for point in points {
            if point.coords.len() != d {
                panic!("Number of points in simplex doesn't match dimension");
            }
            let vertex = Vertex {
                index: i,
            };
            simplex.push(vertex);
            vertices.insert(vertex, point);
            i += 1;
        }

        SimplicialComplex {
            simplices: vec![Simplex {
                vertices: simplex,
            }],
            vertices,
            max_index: i,
        }
    }

    /// Adds a point as a new vertex 
    pub fn add_point(&mut self, point: Point) -> Vertex {
        // TODO: If this point is already in the Complex, glue it to the old Vertex
        let vertex = Vertex {
            index: self.max_index,
        };
        self.max_index += 1;
        self.vertices.insert(vertex, point);
        vertex
    }

    /// Returns the barycenter of a simplex
    pub fn barycenter(&self, simplex: Simplex) -> Point {
        let d = self.dimension();
        let f = 1.0/((d + 1) as f64);
        let mut barycenter = Point::origin(d);
        for i in 0..d {
            for v in simplex.vertices {
                let p = self.vertex_as_point(v).coords;
                barycenter.coords[i] = f * p[i];
            }
        }
        barycenter
    }

    /// Returns the dimension
    pub fn dimension(&self) -> usize {
        if self.simplices.is_empty() {
            panic!("Empty complex")
        }
        self.simplices[0].dimension()
    }

    // Replaces a Simplex in self with its barycentric subdivision
    pub fn subdivide(&mut self, simplex: Simplex) {
        let capacity = self.simplices.len() + factorial((simplex.dimension() + 1).try_into().unwrap())
        let new_simplices = Vec::with_capacity(capacity);
        for i in 0..self.simplices.len() {
            if self.simplices[i] == simplex {
                for subsimplex in self.subdivide_recursive(simplex) {
                    new_simplices.push(subsimplex);
                }
            } else {
                new_simplices.push(self.simplices[i]);
            }
        }
        self.simplices = new_simplices;
    }

    /// Returns the barycentric subdivision of a Simplex
    /// Helper function for subdivide()
    fn subdivide_recursive(&mut self, simplex: Simplex) -> Vec<Simplex> {
        let d = simplex.dimension();
        if d <= 0 {
            vec![simplex]
        } else {
            let v = Vec::with_capacity(factorial((d + 1).try_into().unwrap()));
            let i = self.simplices.iter().position(|r| r == simplex).unwrap();
            for f in simplex.faces() {
                for g in self.subdivide_recursive(f) {
                    v.push(g.append_vertex(self.barycenter(simplex)));
                }
            }
            v
        }
    }

    /// Return the concrete Point referred to by a Vertex
    pub fn vertex_as_point(&self, vertex: Vertex) -> Point {
        self.vertices.get(vertex)
    }
}

fn factorial(num: u128) -> u128 {
    match num {
        0 => 1,
        1.. => num * factorial(num - 1),
    }
}