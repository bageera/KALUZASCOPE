# Mathematical foundations

Conventions: c = 1, metric signature (+,-,-,-) in 4D, (+,-,-,-,-) in 5D.

## The Kaluza–Klein ansatz

The five-dimensional metric is built from the 4D metric g_mn, the
electromagnetic potential A_mu, the dilaton phi, and the coupling kappa:

  G_AB = [ g_mn + phi^2 kappa^2 A_m A_n   phi^2 kappa A_m ]
         [ phi^2 kappa A_n                phi^2         ]

with A_0 the scalar potential. The fifth coordinate x^5 is periodic with
radius R (x^5 ~ 2 pi R y); fields independent of y are the "zero modes".

Standard results implemented and tested here:

| Claim                                   | Where proven                        |
|-----------------------------------------|-------------------------------------|
| det G = phi^2 det g                     | metric::tests::block_determinant_formula |
| Vacuum (A=0, phi=1, g=eta) is 5D Minkowski, det = -1 | metric::tests::vacuum_metric_is_5d_minkowski |
| G_AB invertible; G^AB G_BC = delta_AC   | metric::tests::metric_inverse_roundtrip |
| F_mn = d_m A_n - d_n A_m is antisymmetric, diagonal zero — true by construction | electromagnetics::tests::antisymmetry_and_diagonal_zero |
| E_i = F_0i, B from F_ij = -eps_ijk B_k; S = B^2 - E^2, P = E.B | electromagnetics::tests::invariants_of_pure_e_field |
| d_[a F_bc] = 0 (Bianchi) holds identically for F = dA | reduction::tests::bianchi_holds_for_plane_wave_all_points |
| Sourcefree Maxwell has plane-wave solutions with w = k (lightlike) | reduction::tests::inhomogeneous_residual_zero_for_sourcefree_lightlike_wave |
| The Maxwell residual DETECTS wrong dispersion: w != k gives O(w^2-k^2) residual | reduction::tests::inhomogeneous_residual_exposes_wrong_dispersion |

## Numerical method

Fields live on a uniform lattice indexed f[t][x][y][z], spacings dx =
(dt, dx, dy, dz). First derivatives: central differences, one-sided at the
two boundary faces of each axis.

Known limitation (documented, not hidden): the one-sided boundary stencils
have O(dx) error. Any residual involving SECOND derivatives (inhomogeneous
Maxwell, wave operator) is therefore evaluated only on the 2-point interior
of the lattice; the boundary band is excluded from the residual
comparison. Bianchi residuals involve only first derivatives of F and are
measured on the full lattice.

The plane-wave tests exploit exactness: sin/cos sampled on the lattice have
exact central-difference derivatives (the stencil integrates the mode
exactly), so residuals collapse to floating-point noise (~1e-15) rather
than truncation error — while the w != k case still produces a large,
analytically predictable residual (5.0 = |w^2 - k^2| max |sin|).

## What is NOT claimed

- The 5D metric built from observed A_mu is a mathematical embedding, not
  an independent extra-dimensional measurement.
- No compactification-radius prediction, no new-force observable exists in
  this code yet.
- Current HDD interfaces expose error counters and telemetry, not raw
  read-head waveforms; nothing in the acquisition crate treats SMART
  attributes as field samples.