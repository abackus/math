use beams::global::Domain;
use beams::data::{PointMasses, EssentialBoundaryData};

const L: f64 = 5.0;
const A: f64 = 2.0;
const P: f64 = 5.0;
const Q: f64 = 3.0;
const F: f64 = -1.0;
const EI: f64 = 1.0;

fn main() {
    let domain = Domain::new(
         // the ith entry here is the length of the ith interval
         // in the partition
        vec![L/10., L/10., L/10., L/10., L/10.,
            L/10., L/10., L/10., L/10., L/10.,
            A/5., A/5., A/5., A/5., A/5.
        ]
    );
    let mut stiffness = domain.beam_stiffness(EI);
    let mut load = domain.beam_load(F);

    let p_dof = domain.elts[4].right.dof0;
    let q_dof = domain.elts[14].right.dof0;
    
    let mut pointmasses = PointMasses::new();
    pointmasses.insert(
        p_dof,
        P
    );
    pointmasses.insert(
        q_dof,
        Q
    );
    pointmasses.apply(&mut load);

    EssentialBoundaryData::new(
        // each entry in this vec consists of a dof which is set to 0
        vec![
            domain.elts[0].left.dof0,
            domain.elts[9].right.dof0, 
        ].into_iter().collect()
    ).apply(&mut load, &mut stiffness);
    println!("{}", stiffness);
    println!("{}", load);
    println!("{}", stiffness.try_inverse().unwrap() * load);
}
