use poisson::graph::{Graph, Function, Measure};

fn main() {
    let mut g = Graph::new();
    let v = g.add_vertex();
    g.add_vertices(99);
    let mu = Measure::dirac(&g, v);
    let nu = g.uniform_measure();
    let f = Function {
        domain: &g,
        func: Box::new(|x| {
            15.0 * mu.eval(x)
        }),
    };
    println!("{} {}", mu.lebesgue_norm(&f, 1.0), nu.lebesgue_norm(&f, 1.0));
    println!("{} {}", mu.lebesgue_norm(&f, 0.0), nu.lebesgue_norm(&f, 0.0));
    println!("{} {}", mu.lebesgue_norm(&f, f64::INFINITY), nu.lebesgue_norm(&f, f64::INFINITY));
    println!("{} {}", mu.lebesgue_norm(&f, 10.0), nu.lebesgue_norm(&f, 100.0));
}
