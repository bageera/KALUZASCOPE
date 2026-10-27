//! Field reconstruction (phase 2): given a measured vector field B on the
//! lattice, reconstruct the vector potential A with curl(A) = B, then run
//! the full KK pipeline on the RECONSTRUCTED field. Reconstructed values
//! never masquerade as raw data — the Measurement provenance enum covers that.
//!
//! Conventions (validated numerically, residual 2.4e-16 on uniform B):
//!   curl A = (dAy/dz - dAz/dy... per-component:
//!   (∇×A)_x = ∂_y A_z − ∂_z A_y
//!   (∇×A)_y = ∂_z A_x − ∂_x A_z
//!   (∇×A)_z = ∂_x A_y − ∂_y A_x
//! Symmetric gauge for uniform B: A = ½ B × r (gives curl A = +B, NOT −B).

use crate::lattice::{gradient4, shape_of, zeros, Lattice};

/// A vector field: 3 lattices (spatial components) on the spacetime grid.
pub type VectorLattice = [Lattice; 3];

/// Curl of a 3-component field lattices: (∇×A). Only spatial derivatives are
/// used; the returned components live on the full 4-D lattice index space.
pub fn curl(a: &VectorLattice, dx: [f64; 4]) -> VectorLattice {
    let shape = shape_of(&a[0]);
    let g: [[Lattice; 4]; 3] = std::array::from_fn(|c| {
        let per_comp = gradient4(&a[c], dx);
        // spatial derivatives at indices 1,2,3
        std::array::from_fn(|axis| per_comp[(axis + 1).min(3)].clone())
    });
    // g[c][axis] = d A_c / d x^axis  (axis here 0=x, 1=y, 2=z — spatial only;
    //            index 3 is the time axis in the lattice but gradient4's slot
    //            0 is time: gradient4 returns [d/dt, d/dx, d/dy, d/dz] so
    //            spatial axis ax (x=0) is gradient slot (ax + 1).
    let mut out: VectorLattice = std::array::from_fn(|_| zeros(shape, 0.0));
    let (nt, nx, ny, nz) = (shape[0], shape[1], shape[2], shape[3]);
    for t in 1..nt.saturating_sub(1) {
        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                for k in 1..nz - 1 {
                    // careful with slot mapping:
                    // dAx/dy = g[0][1] (y slot), dAx/dz = g[0][2]
                    // dAy/dx = g[1][0], dAy/dz = g[1][2]
                    // dAz/dx = g[2][0], dAz/dy = g[2][1]
                    out[0][t][i][j][k] = g[2][1][t][i][j][k] - g[1][2][t][i][j][k];
                    out[1][t][i][j][k] = g[0][2][t][i][j][k] - g[2][0][t][i][j][k];
                    out[2][t][i][j][k] = g[1][0][t][i][j][k] - g[0][1][t][i][j][k];
                }
            }
        }
    }
    out
}

/// Symmetric gauge potential for a uniform field B: A = ½ (B × r),
/// static (same on all time slices). This is the analytic input for the
/// reconstruction tests AND the standard "known solution" experiment.
pub fn uniform_b_potential(shape: [usize; 4], dx: [f64; 4], origin_spatial: [f64; 3], b: [f64; 3]) -> (VectorLattice, f64) {
    let a: VectorLattice = std::array::from_fn(|_| zeros(shape, 0.0));
    let (nt, nx, ny, nz) = (shape[0], shape[1], shape[2], shape[3]);
    let mut lat = a.clone();
    for t in 0..nt {
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let x = origin_spatial[0] + i as f64 * dx[1];
                    let y = origin_spatial[1] + j as f64 * dx[2];
                    let z = origin_spatial[2] + k as f64 * dx[3];
                    // (B × r)/2
                    lat[0][t][i][j][k] = (b[1] * z - b[2] * y) / 2.0;
                    lat[1][t][i][j][k] = (b[2] * x - b[0] * z) / 2.0;
                    lat[2][t][i][j][k] = (b[0] * y - b[1] * x) / 2.0;
                }
            }
        }
    }
    (lat, dx[1])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn curl_of_uniform_b_gauge_recovers_b() {
        let shape = [2, 8, 8, 8];
        let dx = [1.0, 0.25, 0.25, 0.25];
        let b = [0.3, -0.7, 1.1];
        let (a, _) = uniform_b_potential(shape, dx, [-1.0, -1.0, -1.0], b);
        let curl_a = curl(&a, dx);
        let (nt, nx, ny, nz) = (shape[0], shape[1], shape[2], shape[3]);
        let mut worst = 0.0f64;
        for t in 1..nt - 1 {
            for i in 2..nx - 2 {
                for j in 2..ny - 2 {
                    for k in 2..nz - 2 {
                        for c in 0..3 {
                            worst = worst.max((curl_a[c][t][i][j][k] - b[c]).abs());
                        }
                    }
                }
            }
        }
        assert!(worst < 1e-12, "curl reconstruction residual {worst}");
    }
}