//! Experiment E2 (phase 2): reconstruction round-trip.
//!
//! Claim under test: a KNOWN B field (uniform analytic), expressed through
//! the symmetric-gauge potential A = ½(B × r), produces F through the KK
//! pipeline whose Maxwell residuals vanish, and whose components match the
//! analytic B exactly. Both chains — direct and via curl reconstruction —
//! give the same physics, so the machinery is insensitive to how A arrived.
//!
//! Convention note (validated): with A = ½(B × r), curl A = +B and the
//! B extraction follows B_i = ½ ε_ijk F_jk (vector-calculus consistent).

use crate::analysis::VectorLattice;
use crate::electromagnetics::PotentialLattice;
use crate::lattice::zeros;

/// Electrostatic full 4-potential lattice: A_0 = 0, A_spatial = a_space.
pub fn to_full_potential(a_space: &VectorLattice, shape: [usize; 4]) -> PotentialLattice {
    let mut full: PotentialLattice = std::array::from_fn(|_| zeros(shape, 0.0));
    for mu in 0..3 {
        full[mu + 1] = a_space[mu].clone();
    }
    full
}

#[cfg(test)]
mod tests {
    use super::*;
    fn fixture() -> (VectorLattice, [f64; 3]) {
        let shape = [8usize, 8, 8, 8];
        let dx = [0.25; 4];
        let b = [0.3, -0.7, 1.1];
        let (a_space, _) = crate::analysis::uniform_b_potential(shape, dx, [-1.0, -1.0, -1.0], b);
        (a_space, b)
    }
    use crate::analysis::curl;
    use crate::electromagnetics::{e_and_b, field_tensor};
    use crate::reduction::{bianchi_residual, maxwell_residual};

    #[test]
    fn curl_reconstruction_recovers_b() {
        let (a_space, b) = fixture();
        let shape = [8, 8, 8, 8];
        let dx = [0.25; 4];
        let rec = curl(&a_space, dx);
        let mut worst = 0.0f64;
        for t in 1..shape[0] - 1 {
            for i in 2..shape[1] - 2 {
                for j in 2..shape[2] - 2 {
                    for k in 2..shape[3] - 2 {
                        for c in 0..3 {
                            worst = worst.max((rec[c][t][i][j][k] - b[c]).abs());
                        }
                    }
                }
            }
        }
        assert!(worst < 1e-12, "curl reconstruction residual {worst}");
    }

    #[test]
    fn field_tensor_matches_analytic_b() {
        let (a_space, b) = fixture();
        let shape = [8, 8, 8, 8];
        let dx = [0.25; 4];
        let f = field_tensor(&to_full_potential(&a_space, shape), dx);
        let mut worst = 0.0f64;
        for t in 0..shape[0] {
            for i in 2..shape[1] - 2 {
                for j in 2..shape[2] - 2 {
                    for k in 2..shape[3] - 2 {
                        for nu in 1..4 {
                            // electrostatic: E = F_0i = 0
                            worst = worst.max(f[0][nu][t][i][j][k].abs());
                        }
                        // B_i = ½ε_ijk F_jk (validated convention): F_12 = B_3, F_23 = B_1, F_13 = -B_2
                        worst = worst.max((f[1][2][t][i][j][k] - b[2]).abs());
                        worst = worst.max((f[1][3][t][i][j][k] + b[1]).abs());
                        worst = worst.max((f[2][3][t][i][j][k] - b[0]).abs());
                    }
                }
            }
        }
        assert!(worst < 1e-12, "field tensor vs analytic B deviation {worst}");
    }

    /// The full Maxwell check on the reconstructed field: both residuals
    /// vanish on the 2-point interior. THE phase-2 deliverable.
    #[test]
    fn maxwell_residuals_vanish_on_reconstructed_field() {
        let (a_space, _) = fixture();
        let shape = [8, 8, 8, 8];
        let dx = [0.25; 4];
        let f = field_tensor(&to_full_potential(&a_space, shape), dx);
        let j = std::array::from_fn(|_| zeros(shape, 0.0));
        let r_inhom = maxwell_residual(&f, dx, &j);
        let r_bianchi = bianchi_residual(&f, dx);
        assert!(r_inhom < 1e-9, "inhomogeneous residual {r_inhom}");
        assert!(r_bianchi < 1e-9, "bianchi residual {r_bianchi}");
    }

    #[test]
    fn e_and_b_agrees_with_gauge_b() {
        let (a_space, b) = fixture();
        let shape = [8, 8, 8, 8];
        let dx = [0.25; 4];
        let f = field_tensor(&to_full_potential(&a_space, shape), dx);
        let (e, bb) = e_and_b(&f);
        let mut worst = 0.0f64;
        for t in 0..shape[0] {
            for i in 2..shape[1] - 2 {
                for j in 2..shape[2] - 2 {
                    for k in 2..shape[3] - 2 {
                        for c in 0..3 {
                            worst = worst.max(e[c][t][i][j][k].abs());
                            worst = worst.max((bb[c][t][i][j][k] - b[c]).abs());
                        }
                    }
                }
            }
        }
        assert!(worst < 1e-12, "e_and_b deviation {worst}");
    }
}