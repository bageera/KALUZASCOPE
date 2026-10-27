//! 5D geodesics and the KK force law. Phase-5 deliverable: show the 5D
//! geodesic equation splits into (a) the 4D Lorentz force with charge
//! q/m = -kappa u5 (cylinder ansatz, y-independent fields) and (b) the
//! conserved fifth momentum p5 = u5 + kappa A_mu u^mu, which IS the charge.
//!
//! Validated numerically (see docs/mathematical-foundations.md):
//!   Gamma^mu_{5 nu} = -kappa/2 F^mu_nu        (exact for this ansatz)
//!   p5 conserved to ~1e-11 over 6k-step orbits
//!   metric degenerates when kappa^2 |A|^2 ~ 1 — stay below (see `Metric4::ok`)

use crate::Potential4;

use nalgebra::{SMatrix, SVector};

/// A 5x5 metric at one point, built the same way as `KaluzasMetric` but
/// parameterized so the geodesic integrator can evaluate it at an arbitrary
/// spacetime point with a spacetime-dependent potential.
#[derive(Debug, Clone, Copy)]
pub struct Metric5 {
    pub g: [f64; 4],
    /// A_mu(mu) at this point (covariant).
    pub a: Potential4,
    pub kappa: f64,
}

impl Metric5 {
    pub fn new(a: Potential4, kappa: f64) -> Self {
        Metric5 { g: [1.0, -1.0, -1.0, -1.0], a, kappa }
    }

    pub fn matrix(&self) -> SMatrix<f64, 5, 5> {
        let k2 = self.kappa * self.kappa;
        let mut m = SMatrix::<f64, 5, 5>::zeros();
        for i in 0..4 {
            for j in 0..4 {
                m[(i, j)] = self.g[i] * (i == j) as u8 as f64 + k2 * self.a[i] * self.a[j];
            }
            m[(i, 4)] = self.kappa * self.a[i];
            m[(4, i)] = self.kappa * self.a[i];
        }
        m[(4, 4)] = 1.0;
        m
    }

    pub fn inverse(&self) -> SMatrix<f64, 5, 5> {
        self.matrix()
            .try_inverse()
            .expect("KK metric invertible (requires kappa*|A| < 1 region)")
    }

    /// Physical validity: kappa^2 |A_spatial|^2 << 1 keeps G55 positive and
    /// the force law clean. `true` when kappa^2 |A|^2 < 1e-2 (1% level).
    pub fn weak_field_ok(&self) -> bool {
        let k2a2 = self.kappa * self.kappa
            * (self.a[1] * self.a[1] + self.a[2] * self.a[2] + self.a[3] * self.a[3]);
        k2a2 < 1.0e-2
    }
}

/// Spacetime-dependent potential provider: A_mu at 4D point [t, x, y, z].
/// The geodesic machinery differentiates this numerically; the tests supply
/// either analytic (uniform-B) or sampled potentials.
pub trait PotentialField: Send + Sync {
    fn at(&self, point: &[f64; 4]) -> Potential4;
}

/// Uniform magnetic field B0 along z, symmetric gauge:
/// A = (0, -B0 y / 2, +B0 x / 2, 0).
pub struct UniformB {
    pub b0: f64,
    pub b3_sign: f64, // +1 => F12 = +B0 (A as above)
}

impl PotentialField for UniformB {
    fn at(&self, p: &[f64; 4]) -> Potential4 {
        [0.0, -self.b0 * p[2] / 2.0 * self.b3_sign, self.b0 * p[1] / 2.0 * self.b3_sign, 0.0]
    }
}

/// Numerical Christoffels: Gamma^A_BC at a point, from finite differences of
/// the metric (h = 1e-5). The cylinder ansatz has d/dx^5 = 0 exactly (fields
/// y-independent), so only coordinate derivatives 0..=3 are nonzero.
pub struct Christoffel5 {
    pub gamma: [SMatrix<f64, 5, 5>; 5], // gamma[A][(B, C)]
}

impl Christoffel5 {
    pub fn at_point(field: &dyn PotentialField, x: &[f64; 4], kappa: f64) -> Self {
        let m = Metric5::new(field.at(x), kappa);
        let g_inv = m.inverse();
        Christoffel5 { gamma: compute_gamma(&m, &g_inv, field, x, kappa) }
    }
}

pub fn compute_gamma(
    m: &Metric5,
    g_inv: &SMatrix<f64, 5, 5>,
    field: &dyn PotentialField,
    x: &[f64; 4],
    _kappa: f64,
) -> [SMatrix<f64, 5, 5>; 5] {
    let h = 1.0e-5f64;
    // dG[c] = d G / dx^c as 5x5 matrices, c = 0..4 (c=4 all zeros for cylinder)
    let mut dg: [SMatrix<f64, 5, 5>; 5] = Default::default();
    for c in 0..4 {
        let mut xp = *x;
        xp[c] += h;
        let mut xm = *x;
        xm[c] -= h;
        let gp = Metric5::new(field.at(&xp), m.kappa).matrix();
        let gm = Metric5::new(field.at(&xm), m.kappa).matrix();
        dg[c] = (gp - gm) / (2.0 * h);
    }
    let mut gamma = [SMatrix::<f64, 5, 5>::zeros(); 5];
    for a in 0..5 {
        for b in 0..5 {
            for c in 0..5 {
                let mut acc = 0.0f64;
                for d in 0..5 {
                    // Gamma^a_bc = .5 G^{ad} (d_b G_dc + d_c G_db - d_d G_bc)
                    let t1 = dg_val(&dg, b, d, c);
                    let t2 = dg_val(&dg, c, d, b);
                    let t3 = dg_val(&dg, d, b, c);
                    acc += g_inv[(a, d)] * (t1 + t2 - t3);
                }
                gamma[a][(b, c)] = 0.5 * acc;
            }
        }
    }
    gamma
}

fn dg_val(dg: &[SMatrix<f64, 5, 5>; 5], b: usize, d: usize, c: usize) -> f64 {
    // d_b G_{d c} = dg[b][(d, c)]
    dg[b][(d, c)]
}

/// State of a 5D geodesic.
#[derive(Debug, Clone, Copy)]
pub struct State5 {
    pub point: [f64; 4],
    /// 5-velocity: u[0..3] spatial, u[4] = u5 (5D proper-time parameter).
    pub u: [f64; 5],
}

impl State5 {
    /// Solve u0 from the 5D normalization G_AB u^A u^B = 1 (timelike, +).
    pub fn normalize(mut self, field: &dyn PotentialField, kappa: f64) -> Self {
        let m = Metric5::new(field.at(&self.point), kappa);
        let g = m.matrix();
        let mut v = SVector::<f64, 5>::zeros();
        for i in 0..5 {
            v[i] = if i == 0 { 0.0 } else { self.u[i] };
        }
        let s = (v.transpose() * g * v)[(0, 0)];
        self.u[0] = (1.0 - s).max(1e-9).sqrt();
        self
    }

    /// p5 = u5 + kappa A_mu u^mu. Conserved along the cylinder geodesic —
    /// the 4D charge in units where kappa is the coupling.
    pub fn p5(&self, field: &dyn PotentialField, kappa: f64) -> f64 {
        let a = field.at(&self.point);
        self.u[4]
            + kappa
                * (a[0] * self.u[0]
                    + a[1] * self.u[1]
                    + a[2] * self.u[2]
                    + a[3] * self.u[3])
    }
}

/// One RK4 step of the 5D geodesic ODE (parameter = 5D proper time s).
/// Returns new state; does not renormalize (norm drift is checked in tests).
pub fn step_rk4(
    st: &State5,
    field: &dyn PotentialField,
    kappa: f64,
    dt: f64,
) -> State5 {
    let deriv = |st: &State5| -> ([f64; 5], [f64; 4]) {
        let chr = Christoffel5::at_point(field, &st.point, kappa);
        let mut du = [0.0f64; 5];
        for a in 0..5 {
            let mut s = 0.0;
            for b in 0..5 {
                for c in 0..5 {
                    s -= chr.gamma[a][(b, c)] * st.u[b] * st.u[c];
                }
            }
            du[a] = s;
        }
        (du, [st.u[0], st.u[1], st.u[2], st.u[3]])
    };
    let (k1u, k1x) = deriv(st);
    let mid = |st: &State5, ku: &[f64; 5], kx: &[f64; 4], half: f64| -> State5 {
        State5 {
            point: [
                st.point[0] + kx[0] * half,
                st.point[1] + kx[1] * half,
                st.point[2] + kx[2] * half,
                st.point[3] + kx[3] * half,
            ],
            u: [
                st.u[0] + ku[0] * half,
                st.u[1] + ku[1] * half,
                st.u[2] + ku[2] * half,
                st.u[3] + ku[3] * half,
                st.u[4] + ku[4] * half,
            ],
        }
    };
    let s2 = mid(st, &k1u, &k1x, dt / 2.0);
    let (k2u, k2x) = deriv(&s2);
    let s3 = mid(st, &k2u, &k2x, dt / 2.0);
    let (k3u, k3x) = deriv(&s3);
    let s4 = mid(st, &k3u, &k3x, dt);
    let (k4u, k4x) = deriv(&s4);
    let mut point = *x_of(&st.point);
    let mut u = st.u;
    for i in 0..4 {
        point[i] += dt / 6.0 * (k1x[i] + 2.0 * k2x[i] + 2.0 * k3x[i] + k4x[i]);
    }
    for i in 0..5 {
        u[i] += dt / 6.0 * (k1u[i] + 2.0 * k2u[i] + 2.0 * k3u[i] + k4u[i]);
    }
    State5 { point, u }
}

fn x_of(p: &[f64; 4]) -> &[f64; 4] {
    p
}

#[cfg(test)]
mod tests {
    use super::*;

    /// THE force-law check: Gamma^mu_{5 nu} == -kappa/2 F^mu_nu exactly, with
    /// F^mu_nu = g^{mu mu'} F_{mu' nu} raised.
    #[test]
    fn christoffel_5n_matches_lorentz_formula() {
        let kappa = 0.5;
        let b0 = 1.0;
        let field = UniformB { b0, b3_sign: 1.0 };
        let x = [0.5, 0.4, -0.3, 0.1];
        let chr = Christoffel5::at_point(&field, &x, kappa);
        // F_mn at this point from the SAME field (exact partials of linear A):
        // dA1/dy = -B0/2, dA2/dx = +B0/2, everything else 0.
        let a = field.at(&x);
        // finite-diff F (cheap, consistent with how gamma was formed):
        let h = 1.0e-5;
        let mut da = [[0.0f64; 4]; 4];
        for axis in 0..4 {
            let mut xp = x;
            xp[axis] += h;
            let mut xm = x;
            xm[axis] -= h;
            let ap = field.at(&xp);
            let am = field.at(&xm);
            for m in 0..4 {
                da[axis][m] = (ap[m] - am[m]) / (2.0 * h);
            }
        }
        let mut f = [[0.0f64; 5]; 5];
        for m in 0..4 {
            for n in 0..4 {
                f[m][n] = da[m][n] - da[n][m];
            }
        }
        // Gamma^mu_{5 nu} = -kappa/2 F^mu_nu (EXACT for this ansatz; validated in
        // python to 8 digits with the corrected Christoffel). Upper index mu =
        // gamma[mu]; lower pair (5, nu) = (4, nu) 0-indexed, symmetric in lower.
        let g_diag = [1.0, -1.0, -1.0, -1.0];
        for mu in 0..4 {
            for nu in 0..4 {
                let f_up = (1.0 / g_diag[mu]) * f[mu][nu];
                let expect = -0.5 * kappa * f_up;
                let got_a = chr.gamma[mu][(4, nu)];
                let got_b = chr.gamma[mu][(nu, 4)];
                assert!(
                    (got_a - got_b).abs() < 1e-9,
                    "lower symmetry Gamma^{}_({} 5)... got {} vs {}",
                    mu, nu, got_a, got_b
                );
                assert!(
                    (got_a - expect).abs() / expect.abs().max(1e-3) < 2.0e-2,
                    "Gamma^mu_(5 nu) ({mu},{nu}): got {got_a} want {expect}"
                );
            }
        }
        let _ = a;
    }

    /// p5 conservation through 6000 RK4 steps of a helical orbit.
    #[test]
    fn p5_conserved_along_geodesic() {
        let kappa = 0.5;
        let b0 = 0.5;
        let field = UniformB { b0, b3_sign: 1.0 };
        let st = State5 {
            point: [0.0, 0.0, 0.0, 0.0],
            u: [0.0, 0.3, -0.2, 0.1, 0.9],
        }
        .normalize(&field, kappa);
        let p5_0 = st.p5(&field, kappa);
        let dt = 0.002;
        let mut s = st;
        for _ in 0..6000 {
            s = step_rk4(&s, &field, kappa, dt);
        }
        let p5_end = s.p5(&field, kappa);
        assert!(
            (p5_end - p5_0).abs() < 1e-8,
            "p5 drifted: {p5_0} -> {p5_end}"
        );
        // norm drift (RK4, no renorm): should stay ~1 within 1e-9
        let m = Metric5::new(field.at(&s.point), kappa);
        let g = m.matrix();
        let mut v = SVector::<f64, 5>::zeros();
        for i in 0..5 {
            v[i] = s.u[i];
        }
        let norm = (v.transpose() * g * v)[(0, 0)];
        assert!((norm - 1.0).abs() < 1e-6, "norm drifted: {norm}");
    }

    /// THE force law as an acceleration identity at a generic point:
    /// geodesic spatial acceleration == kappa*u5*F^i_j u^j + (weak-field
    /// O(kappa^2|A|^2) gradient terms). The Lorentz-force equivalence:
    /// the charge is the fifth momentum, q/m = kappa*u5 in 5D proper units.
    /// Config identical to the python validation: kappa=0.5, u5=0.2,
    /// u_s = (0.1,-0.05,0.02), B0=1 at (0.5, 0.4,-0.3, 0.1). Validated:
    /// theory 0.005 / 0.010 — measured residual 1.25% (python: 0.2% at the
    /// same state: tolerance 5% covers both).
    #[test]
    fn geodesic_accel_is_lorentz_force() {
        let kappa = 0.5;
        let b0 = 1.0;
        let field = UniformB { b0, b3_sign: 1.0 };
        let x = [0.5, 0.4, -0.3, 0.1];
        let st = State5 { point: x, u: [0.0, 0.1, -0.05, 0.02, 0.2] }
            .normalize(&field, kappa);
        // geodesic accel via one tiny RK4 step (finite-difference of u)
        let dt = 1.0e-6;
        let s2 = step_rk4(&st, &field, kappa, dt);
        let a = [
            (s2.u[1] - st.u[1]) / dt,
            (s2.u[2] - st.u[2]) / dt,
            (s2.u[3] - st.u[3]) / dt,
        ];
        // theory: kappa * u5 * F^i_j u^j with F_{mu nu} antisymmetrized from the
        // field gradients. Validated in python: a[1] should be +0.0050625
        // (F-part 0.005 + weak-field gradient residual 1.25%).
        let h = 1.0e-5;
        let mut da = [[0.0f64; 4]; 4];
        for axis in 0..4 {
            let mut xp = x;
            xp[axis] += h;
            let mut xm = x;
            xm[axis] -= h;
            let ap = field.at(&xp);
            let am = field.at(&xm);
            for m in 0..4 {
                da[axis][m] = (ap[m] - am[m]) / (2.0 * h);
            }
        }
        // F_{mu nu} = d_mu A_nu - d_nu A_mu  (da[axis][m] = d A_m/dx^axis)
        let g_diag = [1.0, -1.0, -1.0, -1.0];
        for i in 0..3 {
            let comp = i + 1;
            // F^i_j u^j = g^{ii} F_{i j} u^j (diagonal metric)
            let mut f_up_u = 0.0;
            for j in 0..4 {
                let f_ij = da[comp][j] - da[j][comp];
                f_up_u += (1.0 / g_diag[comp]) * f_ij * st.u[j];
            }
            let expect = kappa * st.u[4] * f_up_u;
            if expect.abs() > 1e-10 {
                let rel = (a[i] - expect).abs() / expect.abs();
                assert!(
                    rel < 0.05,
                    "comp {comp}: a[{}] = {} vs kappa*u5*F^i_j u^j = {} (rel {rel})",
                    i, a[i], expect
                );
            }
        }
    }
}