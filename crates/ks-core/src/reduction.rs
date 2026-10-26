//! Dimensional reduction: recover 4D Maxwell from the 5D story.
//!
//! The first success criterion (v0.1.0): the KK machinery is consistent with
//! 4D Maxwell on the lattice. Two residual functions do the checking:
//!
//! - `bianchi_residual`: homogeneous Maxwell (d_[a F_bc] = 0). First
//!   derivatives of F — exact on analytically-representable fields, valid on
//!   the full lattice.
//! - `maxwell_residual`: inhomogeneous Maxwell (d_mu F^{mu nu} = j^nu).
//!   SECOND derivatives of A — evaluated only on the 2-point interior since
//!   one-sided boundary stencils have O(dx) error (documented convention).

use crate::lattice::{gradient4, shape_of, Lattice};

const G4: [f64; 4] = [1.0, -1.0, -1.0, -1.0];

/// Inhomogeneous Maxwell residual: max |d_mu F^{mu nu} - j^nu| over the
/// 2-point interior of the lattice and all nu. Central differences only:
/// one-sided boundary stencils carry O(dx) error (documented convention),
/// so the boundary band is excluded. `j` = source lattice per
/// contravariant component (j^0 = rho).
pub fn maxwell_residual(f: &[[Lattice; 4]; 4], dx: [f64; 4], j: &[Lattice; 4]) -> f64 {
    let shape = shape_of(&f[0][1]);
    const MARGIN: usize = 2;
    if shape.iter().any(|&n| n < 2 * MARGIN + 1) {
        panic!("lattice axis smaller than 2*MARGIN+1 = 5: {shape:?}");
    }
    let mut worst = 0.0f64;
    for nu in 0..4 {
        for t in MARGIN..shape[0] - MARGIN { for x in MARGIN..shape[1] - MARGIN {
        for y in MARGIN..shape[2] - MARGIN { for z in MARGIN..shape[3] - MARGIN {
            let mut div = 0.0;
            for mu in 0..4 {
                // d/dx^mu of F^{mu nu}: diagonal metric =>
                // F^{mu nu} = g[mu] g[nu] F_{mu nu}.
                let mut im = [t, x, y, z];
                let mut ip = [t, x, y, z];
                im[mu] -= 1;
                ip[mu] += 1;
                let f_up_ip = G4[mu] * G4[nu] * f[mu][nu][ip[0]][ip[1]][ip[2]][ip[3]];
                let f_up_im = G4[mu] * G4[nu] * f[mu][nu][im[0]][im[1]][im[2]][im[3]];
                div += (f_up_ip - f_up_im) / (2.0 * dx[mu]);
            }
            let residual = (div - j[nu][t][x][y][z]).abs();
            if residual > worst {
                worst = residual;
            }
        }}}}
    }
    worst
}

/// Bianchi identity residual: max |d_a F_bc + d_b F_ca + d_c F_ab| over all
/// a < b < c and lattice points. Uses gradient4 (one-sided at boundary faces
/// is fine — only first derivatives of F are involved).
pub fn bianchi_residual(f: &[[Lattice; 4]; 4], dx: [f64; 4]) -> f64 {
    let shape = shape_of(&f[0][1]);
    let mut worst = 0.0f64;
    for a in 0..4 { for b in 0..4 { for c in 0..4 {
        if a < b && b < c {
            let df_bc = gradient4(&f[b][c], dx);
            let df_ca = gradient4(&f[c][a], dx);
            let df_ab = gradient4(&f[a][b], dx);
            for t in 0..shape[0] { for x in 0..shape[1] { for y in 0..shape[2] { for z in 0..shape[3] {
                let s = df_bc[a][t][x][y][z] + df_ca[b][t][x][y][z] + df_ab[c][t][x][y][z];
                if s.abs() > worst {
                    worst = s.abs();
                }
            }}}}
        }
    }}}
    worst
}

/// Wave operator (c = 1): dd_t^2 - grad^2, central second differences with
/// one-sided at boundaries. Interior accuracy O(dx^2); boundary O(dx).
pub fn wave_operator(field: &Lattice, dx: [f64; 4]) -> Lattice {
    let shape = shape_of(field);
    let mut out = crate::lattice::zeros(shape, 0.0);
    let sgn = [1.0, -1.0, -1.0, -1.0];
    let mut idx = [0usize; 4];
    for i0 in 0..shape[0] { idx[0] = i0;
    for i1 in 0..shape[1] { idx[1] = i1;
    for i2 in 0..shape[2] { idx[2] = i2;
    for i3 in 0..shape[3] { idx[3] = i3;
        let mut acc = 0.0;
        for axis in 0..4 {
            let n = shape[axis];
            let at = |p: [usize; 4]| field[p[0]][p[1]][p[2]][p[3]];
            let mut lo = idx;
            let mut hi = idx;
            if idx[axis] == 0 {
                hi[axis] = 1;
                acc += sgn[axis] * (at(hi) - 2.0 * at(idx) + at(idx)) / dx[axis].powi(2);
            } else if idx[axis] == n - 1 {
                lo[axis] = n - 2;
                acc += sgn[axis] * (at(idx) - 2.0 * at(idx) + at(lo)) / dx[axis].powi(2);
            } else {
                lo[axis] = idx[axis] - 1;
                hi[axis] = idx[axis] + 1;
                acc += sgn[axis] * (at(hi) - 2.0 * at(idx) + at(lo)) / dx[axis].powi(2);
            }
        }
        out[idx[0]][idx[1]][idx[2]][idx[3]] = acc;
    }}}}
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::electromagnetics::field_tensor;
    use crate::lattice::zeros;

    /// Plane wave A_y = sin(w t - k z), uniform across x and y:
    /// B_x = k cos(...), E_y = w cos(...).
    fn plane_wave_a(
        shape: [usize; 4],
        dt: f64,
        dz: f64,
        w: f64,
        k: f64,
    ) -> crate::electromagnetics::PotentialLattice {
        let mut a = std::array::from_fn(|_| zeros(shape, 0.0));
        for t in 0..shape[0] { for x in 0..shape[1] { for y in 0..shape[2] { for z in 0..shape[3] {
            let phase = w * (t as f64) * dt - k * (z as f64) * dz;
            a[2][t][x][y][z] = phase.sin();
        }}}}
        a
    }

    #[test]
    fn bianchi_holds_for_plane_wave_all_points() {
        let shape = [8, 2, 2, 16];
        let (dt, dz, w, k) = (0.1, 0.1, 2.0, 2.0); // w = k, lightlike wave
        let a = plane_wave_a(shape, dt, dz, w, k);
        let f = field_tensor(&a, [dt, 1.0, 1.0, dz]);
        let r = bianchi_residual(&f, [dt, 1.0, 1.0, dz]);
        // F components are exact cosines; d(cos) is exact up to fp noise.
        assert!(r < 1e-9, "bianchi residual {r}");
    }

    #[test]
    fn inhomogeneous_residual_zero_for_sourcefree_lightlike_wave() {
        // w = k: exact solution of sourcefree Maxwell. Interior residual ~ 0.
        let shape = [10, 5, 5, 16];
        let (dt, dz, w, k) = (0.1, 0.1, 2.0, 2.0);
        let a = plane_wave_a(shape, dt, dz, w, k);
        let f = field_tensor(&a, [dt, 1.0, 1.0, dz]);
        let j = std::array::from_fn(|_| zeros(shape, 0.0));
        let r = maxwell_residual(&f, [dt, 1.0, 1.0, dz], &j);
        assert!(r < 1e-9, "sourcefree residual {r}");
    }

    #[test]
    fn inhomogeneous_residual_exposes_wrong_dispersion() {
        // w != k: NOT a sourcefree solution. Residual must be substantial,
        // ~ (w^2 - k^2) |sin| = 5.0 here.
        let shape = [10, 5, 5, 16];
        let (dt, dz, w, k) = (0.1, 0.1, 3.0, 2.0);
        let a = plane_wave_a(shape, dt, dz, w, k);
        let f = field_tensor(&a, [dt, 1.0, 1.0, dz]);
        let j = std::array::from_fn(|_| zeros(shape, 0.0));
        let r = maxwell_residual(&f, [dt, 1.0, 1.0, dz], &j);
        assert!(r > 4.0, "expected O(5) residual for wrong w/k, got {r}");
    }

    #[test]
    fn wave_operator_on_quadratic_in_time() {
        // f = t^2: dd_t^2 f = 2 => box f = 2 (interior).
        let shape = [6, 2, 2, 2];
        let mut f = zeros(shape, 0.0);
        for t in 0..6 { for x in 0..2 { for y in 0..2 { for z in 0..2 {
            f[t][x][y][z] = (t * t) as f64;
        }}}}
        let out = wave_operator(&f, [1.0; 4]);
        assert!((out[3][0][0][0] - 2.0).abs() < 1e-9);
    }
}