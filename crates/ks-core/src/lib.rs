//! Kaluza–Klein 5D metric machinery. See docs/mathematical-foundations.md.

pub mod analysis;
pub mod electromagnetics;
pub mod experiments_e2;
pub mod geodesic;
pub mod lattice;
pub mod metric;
pub mod reduction;

/// Covariant 4-potential A_mu (A_0 = scalar potential, c = 1).
pub type Potential4 = [f64; 4];

/// Generic 4-vector.
pub type Vec4 = [f64; 4];

/// Shared numeric assertion: relative-tolerance closeness with label.
pub fn assert_close(a: f64, b: f64, tol: f64, label: &str) {
    assert!(
        (a - b).abs() <= tol.abs() * (1.0 + a.abs().max(b.abs())),
        "{label}: {a} vs {b} (tol {tol})"
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn assert_close_passes() {
        assert_close(1.0, 1.0 + 1e-12, 1e-9, "same");
        assert_close(100.0, 100.5, 1e-2, "relative");
    }
}