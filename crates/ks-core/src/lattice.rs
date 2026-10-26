//! Lattice utilities: a scalar field on a uniform 4-D grid, indexed
//! f[t][x][y][z]. Central differences; one-sided at the two boundary faces
//! of each axis. Documented caveat: one-sided boundary stencils carry O(dx)
//! error, so any SECOND-derivative residual must be evaluated only on the
//! interior (2-point margin) — see reduction::maxwell_residual.

/// A scalar field on a spacetime lattice. Indexing: f[t][x][y][z].
pub type Lattice = Vec<Vec<Vec<Vec<f64>>>>;

///Allocate a lattice of `shape` filled with `fill`.
pub fn zeros(shape: [usize; 4], fill: f64) -> Lattice {
    vec![vec![vec![vec![fill; shape[3]]; shape[2]]; shape[1]]; shape[0]]
}

pub fn shape_of(field: &Lattice) -> [usize; 4] {
    [
        field.len(),
        field[0].len(),
        field[0][0].len(),
        field[0][0][0].len(),
    ]
}

fn neighbor_range(idx: usize, n: usize) -> (usize, usize, f64) {
    // Returns (lo, hi, divisor): one-sided (stride 1 per dx) at the two
    // boundary faces, central (stride 2 per 2 dx) in the interior.
    if idx == 0 {
        (0, 1, 1.0)
    } else if idx == n - 1 {
        (n - 2, n - 1, 1.0)
    } else {
        (idx - 1, idx + 1, 2.0)
    }
}

/// Central-difference derivative along every axis; one-sided at the two
/// boundary faces of each axis. `dx` = spacings (dt, dx, dy, dz).
/// Returns [df/dt, df/dx, df/dy, df/dz].
pub fn gradient4(field: &Lattice, dx: [f64; 4]) -> [Lattice; 4] {
    let shape = shape_of(field);
    let mut out: [Lattice; 4] = std::array::from_fn(|_| zeros(shape, 0.0));
    let mut idx = [0usize; 4];
    for i0 in 0..shape[0] {
        idx[0] = i0;
        for i1 in 0..shape[1] {
            idx[1] = i1;
            for i2 in 0..shape[2] {
                idx[2] = i2;
                for i3 in 0..shape[3] {
                    idx[3] = i3;
                    for axis in 0..4 {
                        let n = shape[axis];
                        if n < 2 {
                            continue;
                        }
                        let (lo, hi, div) = neighbor_range(idx[axis], n);
                        let mut a = idx;
                        let mut b = idx;
                        a[axis] = lo;
                        b[axis] = hi;
                        let v = (field[b[0]][b[1]][b[2]][b[3]]
                            - field[a[0]][a[1]][a[2]][a[3]])
                            / (div * dx[axis]);
                        out[axis][idx[0]][idx[1]][idx[2]][idx[3]] = v;
                    }
                }
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gradient_of_linear_field_is_constant() {
        let shape = [2, 5, 2, 2];
        let mut f = zeros(shape, 0.0);
        for t in 0..2 { for x in 0..5 { for y in 0..2 { for z in 0..2 {
            f[t][x][y][z] = x as f64;
        }}}}
        let g = gradient4(&f, [1.0; 4]);
        for t in 0..2 { for x in 0..5 { for y in 0..2 { for z in 0..2 {
            assert!((g[0][t][x][y][z]).abs() < 1e-12, "dt");
            assert!((g[1][t][x][y][z] - 1.0).abs() < 1e-12, "dx at {t},{x}");
            assert!((g[2][t][x][y][z]).abs() < 1e-12, "dy");
            assert!((g[3][t][x][y][z]).abs() < 1e-12, "dz");
        }}}}
    }

    #[test]
    fn gradient_of_quadratic_interior_is_exact() {
        let shape = [2, 6, 2, 2];
        let mut f = zeros(shape, 0.0);
        for t in 0..2 { for x in 0..6 { for y in 0..2 { for z in 0..2 {
            f[t][x][y][z] = (x * x) as f64;
        }}}}
        let g = gradient4(&f, [1.0; 4]);
        for x in 1..5 {
            assert!((g[1][0][x][0][0] - 2.0 * x as f64).abs() < 1e-12);
        }
        // boundary faces use one-sided: x=0 => (1-0)/1 = 1; x=5 => (25-16)/1 = 9
        assert!((g[1][0][0][0][0] - 1.0).abs() < 1e-12);
        assert!((g[1][0][5][0][0] - 9.0).abs() < 1e-12);
    }
}