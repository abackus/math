use beams::global::Domain;
use beams::data::{PointMasses, EssentialBoundaryData};

fn main() {
    let domain = Domain::new(
        vec![0.5, 0.25, 0.25]
    );
    let mut stiffness = domain.beam_stiffness(1.0);
    let mut load = domain.beam_load(1.0);

    let moment_dof = domain.elts[0].right.dof1;
    let pointload_dof = domain.elts[1].right.dof0;
    
    let mut pointmasses = PointMasses::new();
    pointmasses.insert(
        moment_dof,
        0.5
    );
    pointmasses.insert(
        pointload_dof,
        5.0
    );
    pointmasses.apply(&mut load);

    EssentialBoundaryData::new(
        vec![
            domain.elts[0].left.dof0,
            domain.elts[0].left.dof1, 
            domain.elts[0].right.dof0,
            domain.elts[2].right.dof0,
        ].into_iter().collect()
    ).apply(&mut load, &mut stiffness);
    println!("{}", stiffness);
    println!("{}", load);
    println!("{}", stiffness.try_inverse().unwrap() * load);
}
