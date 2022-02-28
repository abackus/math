//! A finite element solver for one-dimensional, linear elliptic PDE
//! We assume that we are working in Sobolev space H_0^1

pub mod local;
pub mod global;
pub mod data;

#[cfg(test)]
mod tests {
    use crate::global::Domain;
    use crate::data::{PointMasses, EssentialBoundaryData};
    use nalgebra::{Vector4, Matrix4}; 

    #[test]
    fn one_element() {
        let domain = Domain::evenly_spaced(2, 1.0);
        let elt = domain.elts[0];
        let good_load = Vector4::new(0.5, 1.0 / 12.0, 0.5, -1.0 / 12.0);
        let good_stiffness = Matrix4::new(
            12.0, 6.0, -12.0, 6.0,
            6.0, 4.0, -6.0, 2.0,
            -12.0, -6.0, 12.0, -6.0,
            6.0, 2.0, -6.0, 4.0
        );
        assert_eq!(elt.beam_load(1.0), good_load);
        assert_eq!(domain.beam_load(1.0), good_load);
        assert_eq!(elt.beam_stiffness(1.0), good_stiffness);
        assert_eq!(domain.beam_stiffness(1.0), good_stiffness);
    }

    #[test]
    fn two_elements() {
        let domain = Domain::evenly_spaced(3, 0.5);
        let stiffness = domain.beam_stiffness(1.0);
        assert_eq!(stiffness[(0, 0)], 96.0);
        assert_eq!(stiffness[(0, 5)], 0.0);
    }

    #[test]
    fn lecture_example() {
        let domain = Domain::new(
            vec![0.5, 0.25, 0.25]
        );
        let mut stiffness = domain.beam_stiffness(1.0);
        let mut load = domain.beam_load(1.0);
        assert_eq!(stiffness[(0, 0)], 96.0);
        assert_eq!(stiffness[(0, 7)], 0.0);

        assert_eq!(domain.elts[0].right.coord, 0.5);
        assert_eq!(domain.elts[1].right.coord, 0.75);
        let moment_dof = domain.elts[0].right.dof1;
        let pointload_dof = domain.elts[1].right.dof0;
        assert_eq!(load[moment_dof.index()], 0.005208333333333333);
        assert_eq!(load[pointload_dof.index()], 0.125);

        let mut pointmasses = PointMasses::new();
        pointmasses.insert(
            moment_dof,
            2.0
        );
        pointmasses.insert(
            pointload_dof,
            5.0
        );
        pointmasses.apply(&mut load);
        assert_eq!(load[moment_dof.index()], 2.005208333333333333);
        assert_eq!(load[pointload_dof.index()], 5.125);

        EssentialBoundaryData::new(
            vec![
                domain.elts[0].left.dof0,
                domain.elts[0].left.dof1, 
                domain.elts[0].right.dof0,
                domain.elts[2].right.dof0,
            ].into_iter().collect()
        ).apply(&mut load, &mut stiffness);
        assert_eq!(load[0], 0.0);
        assert_eq!(stiffness[(0, 0)], 1.0);
        assert_eq!(stiffness[(0, 1)], 0.0);
        assert_eq!(stiffness[(1, 0)], 0.0);
        assert_ne!(stiffness.determinant(), 0.0);
    }
}