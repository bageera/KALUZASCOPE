//! The Kaluza–Klein 5D metric built from the 4D metric, the EM potential,
//! and the dilaton. Convention (standard Kaluza 1921 / Jordan-Thiry form,
//! c = 1, G = 0 coupling scaled out of kappa):
//!
//!   G_AB = [[ g_mn + phi^2 kappa^2 A_m A_n,  phi^2 kappa A_m ],
//!           [ phi^2 kappa A_n,               phi^2       ]]
//!
//! A_mu = [phi, A_i] covariant with A_0 = scalar potential.
//! Signature (+,-,-,-,-), fifth coordinate x^5 ~ 2 pi R y, y periodic.

use crate::Potential4;
use nalgebra::SMatrix;

/// Kappa: EM coupling in the metric. Standard result: kappa^2 = 16 pi G / c^4
/// (G = Newton's constant). Keep it an explicit knob, not a hidden global.
pub const KAPPA_DEFAULT: f64 = 1.0;

/// A 4D spatial metric block. For the flat/curved checks we use diagonal
/// metrics; anything more general can go in later.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Metric4 {
    /// diag entries with signature (+1,-1,-1,-1) assumed.
    pub diag: [f64; 4],
}

impl Metric4 {
    pub fn minkowski() -> Self {
        Metric4 { diag: [1.0, -1.0, -1.0, -1.0] }
    }

    pub fn as_matrix(&self) -> SMatrix<f64, 4, 4> {
        let mut m = SMatrix::<f64, 4, 4>::zeros();
        for i in 0..4 { m[(i, i)] = self.diag[i]; }
        m
    }

    pub fn as_flat(&self) -> [f64; 16] {
        let mut m = [0.0; 16];
        for i in 0..4 { m[i * 4 + i] = self.diag[i]; }
        m
    }

    /// g^{mn} (diagonal inverse).
    pub fn inverse(&self) -> [f64; 4] {
        std::array::from_fn(|i| 1.0 / self.diag[i])
    }
}

/// The 5D Kaluza–Klein metric with a constant scalar field phi and constant
/// vector potential (the "zero-mode" truncation). For a full spacetime-
/// dependent field, pass your per-point A_mu / phi and evaluate here.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KaluzasMetric {
    pub base: Metric4,
    pub a: Potential4,
    pub phi: f64,
    pub kappa: f64,
}

impl KaluzasMetric {
    /// Flat spacetime, phi = 1, no potential — pure 5D Minkowski block form.
    pub fn vacuum(kappa: f64) -> Self {
        KaluzasMetric {
            base: Metric4::minkowski(),
            a: [0.0; 4],
            phi: 1.0,
            kappa,
        }
    }

    /// The full 5x5 matrix G_AB in row-major.
    pub fn matrix(&self) -> SMatrix<f64, 5, 5> {
        let g = self.base.as_matrix();
        let k = self.kappa;
        let p2 = self.phi * self.phi;
        let mut m = SMatrix::<f64, 5, 5>::zeros();
        for i in 0..4 {
            for j in 0..4 {
                m[(i, j)] = g[(i, j)] + p2 * k * k * self.a[i] * self.a[j];
            }
        }
        for i in 0..4 {
            m[(i, 4)] = p2 * k * self.a[i];
            m[(4, i)] = p2 * k * self.a[i];
        }
        m[(4, 4)] = p2;
        m
    }

    pub fn as_flat_array(&self) -> [f64; 25] {
        let m = self.matrix();
        std::array::from_fn(|n| m[n])
    }

    /// Determinant via the block formula: det G = det(g + phi^2 k^2 A A^T
    /// - phi^2 k^2 A A^T) ... standard result det G = phi^2 det g. Verify.
    pub fn determinant(&self) -> f64 {
        self.matrix().determinant()
    }

    /// 5D inverse metric G^AB from the matrix.
    pub fn inverse(&self) -> SMatrix<f64, 5, 5> {
        self.matrix().try_inverse().expect("5D KK metric must be invertible for physical parameters")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vacuum_metric_is_5d_minkowski() {
        let m = KaluzasMetric::vacuum(1.0);
        let mat = m.matrix();
        for i in 0..5 {
            assert!((mat[(i, i)] - if i == 4 { 1.0 } else if i == 0 { 1.0 } else { -1.0 }).abs() < 1e-12, "diag {i}");
            for j in 0..5 {
                if i != j {
                    assert!(mat[(i, j)].abs() < 1e-15);
                }
            }
        }
        assert!((m.determinant() - (-1.0)).abs() < 1e-12);
    }

    #[test]
    fn block_determinant_formula() {
        // det G = phi^2 det g for the KK ansatz — verify directly on a
        // generic set of parameters, independent of the shortcut formula.
        let m = KaluzasMetric {
            base: Metric4::minkowski(),
            a: [0.3, -1.2, 0.7, 2.4],
            phi: 1.7,
            kappa: 0.5,
        };
        let g_det = m.base.as_matrix().determinant();
        assert!((m.determinant() - m.phi * m.phi * g_det).abs() < 1e-10);
    }

    #[test]
    fn metric_inverse_roundtrip() {
        let m = KaluzasMetric {
            base: Metric4::minkowski(),
            a: [0.3, -1.2, 0.7, 2.4],
            phi: 1.7,
            kappa: 0.5,
        };
        let g = m.matrix();
        let gi = m.inverse();
        let eye = g * gi;
        for i in 0..5 {
            for j in 0..5 {
                assert!((eye[(i, j)] - if i == j { 1.0 } else { 0.0 }).abs() < 1e-12);
            }
        }
    }
}