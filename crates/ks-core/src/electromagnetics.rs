//! Electromagnetism on the lattice: potential, field tensor, E/B extraction.
//! c = 1, metric signature (+,-,-,-). Covariant potential A_mu, A_0 = phi.
//! F_{mu nu} = d_mu A_nu - d_nu A_mu with F_{0i} = E_i, F_{ij} = -eps_ijk B_k.

use crate::lattice::{gradient4, shape_of, zeros, Lattice};

/// Potential: a[t][x][y][z] = A_mu(mu) at that point.
pub type PotentialLattice = [Lattice; 4];

/// Field tensor: f[mu][nu] is the lattice of F_{mu nu}. Antisymmetric by
/// construction (built as difference of gradients of A).
pub fn field_tensor(a: &PotentialLattice, dx: [f64; 4]) -> [[Lattice; 4]; 4] {
    let d: [[Lattice; 4]; 4] = std::array::from_fn(|mu| gradient4(&a[mu], dx));
    let mut f: [[Lattice; 4]; 4] = Default::default();
    for mu in 0..4 {
        for nu in 0..4 {
            let shape = shape_of(&a[0]);
            let mut out = zeros(shape, 0.0);
            for (i, (t_u, t_v)) in d[nu][mu].iter().zip(&d[mu][nu]).enumerate() {
                for (x, (x_u, x_v)) in t_u.iter().zip(t_v).enumerate() {
                    for (y, (y_u, y_v)) in x_u.iter().zip(x_v).enumerate() {
                        for (z, (z_u, z_v)) in y_u.iter().zip(y_v).enumerate() {
                            // d[nu][mu] = d_mu A_nu ; d[mu][nu] = d_nu A_mu
                            out[i][x][y][z] = z_u - z_v;
                        }
                    }
                }
            }
            f[mu][nu] = out;
        }
    }
    f
}

/// E_i = F_{0i} lattices cloned out, B from F_{12} = -B_3, F_{13} = +B_2,
/// F_{23} = -B_1. Returns ([E; 3], [B; 3]).
pub fn e_and_b(f: &[[Lattice; 4]; 4]) -> ([Lattice; 3], [Lattice; 3]) {
    let e: [Lattice; 3] = std::array::from_fn(|i| f[0][i + 1].clone());
    let shape = shape_of(&f[1][2]);
    let (nt, nx, ny, nz) = (shape[0], shape[1], shape[2], shape[3]);
    let mut b: [Lattice; 3] = std::array::from_fn(|_| zeros(shape, 0.0));
    for t in 0..nt { for x in 0..nx { for y in 0..ny { for z in 0..nz {
        b[0][t][x][y][z] = -f[2][3][t][x][y][z];
        b[1][t][x][y][z] = f[1][3][t][x][y][z];
        b[2][t][x][y][z] = -f[1][2][t][x][y][z];
    }}}}
    (e, b)
}

/// Lorentz invariants per lattice point, computed from E and B directly:
///   S = B^2 - E^2   (half-trace contracted with (+,-,-,-))
///   P = E . B
pub fn invariants(f: &[[Lattice; 4]; 4]) -> (Lattice, Lattice) {
    let shape = shape_of(&f[0][1]);
    let (nt, nx, ny, nz) = (shape[0], shape[1], shape[2], shape[3]);
    let (e, b) = e_and_b(f);
    let mut s = zeros(shape, 0.0);
    let mut p = zeros(shape, 0.0);
    for t in 0..nt { for x in 0..nx { for y in 0..ny { for z in 0..nz {
        let e2: f64 = (0..3).map(|i| e[i][t][x][y][z] * e[i][t][x][y][z]).sum();
        let b2: f64 = (0..3).map(|i| b[i][t][x][y][z] * b[i][t][x][y][z]).sum();
        let edotb: f64 = (0..3).map(|i| e[i][t][x][y][z] * b[i][t][x][y][z]).sum();
        s[t][x][y][z] = b2 - e2;
        p[t][x][y][z] = edotb;
    }}}}
    (s, p)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_potential_has_zero_field() {
        let shape = [2, 3, 3, 3];
        let a: PotentialLattice = std::array::from_fn(|_| zeros(shape, 0.4));
        let f = field_tensor(&a, [1.0; 4]);
        for mu in 0..4 { for nu in 0..4 {
            assert!(f[mu][nu][0][1][1][1].abs() < 1e-12);
        }}
        let (e, b) = e_and_b(&f);
        for i in 0..3 {
            assert!(e[i][0][1][1][1].abs() < 1e-12);
            assert!(b[i][0][1][1][1].abs() < 1e-12);
        }
    }

    #[test]
    fn antisymmetry_and_diagonal_zero() {
        let shape = [2, 4, 2, 4];
        let mut a: PotentialLattice = std::array::from_fn(|mu| zeros(shape, mu as f64 * 0.1));
        for t in 0..2 { for x in 0..4 { for y in 0..2 { for z in 0..4 {
            for mu in 0..4 {
                a[mu][t][x][y][z] =
                    mu as f64 + t as f64 + x as f64 + 2.0 * y as f64 + 3.0 * z as f64;
            }
        }}}}
        let f = field_tensor(&a, [1.0; 4]);
        for mu in 0..4 {
            assert!(f[mu][mu][0][1][1][1].abs() < 1e-12);
            for nu in 0..4 {
                let at = f[mu][nu][0][1][1][1];
                let at2 = f[nu][mu][0][1][1][1];
                assert!((at + at2).abs() < 1e-12, "antisym {mu},{nu}");
            }
        }
    }

    #[test]
    fn invariants_of_pure_e_field() {
        // A_0 = -E_y y with E_y = 2: E = F_{0i}; at interior points F_02 = -d_2 A_0 = +2.
        let shape = [2, 4, 2, 4];
        let mut a: PotentialLattice = std::array::from_fn(|_| zeros(shape, 0.0));
        for t in 0..2 { for x in 0..4 { for y in 0..2 { for z in 0..4 {
            a[0][t][x][y][z] = -2.0 * y as f64;
        }}}}
        let f = field_tensor(&a, [1.0; 4]);
        let (e, b) = e_and_b(&f);
        // interior y = 1 (y=0 is a boundary face, one-sided diff — also fine
        // here since the field is linear, but we test the interior point).
        assert!((e[1][0][1][1][1] - 2.0).abs() < 1e-12, "E_y = {}", e[1][0][1][1][1]);
        assert!(e[0][0][1][1][1].abs() < 1e-12);
        assert!(b.iter().all(|bl| bl[0][1][1][1].abs() < 1e-12));
        let (s, p) = invariants(&f);
        assert!((s[0][1][1][1] + 4.0).abs() < 1e-12, "S = B^2 - E^2 = -4");
        assert!(p[0][1][1][1].abs() < 1e-12);
    }
}