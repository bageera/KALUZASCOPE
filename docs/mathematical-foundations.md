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
| curl(A) = B for A = 1/2(B x r); reconstruction residual < 1e-12 | experiments_e2::tests::curl_reconstruction_recovers_b |
| B extraction B_i = 1/2 eps_ijk F_jk matches vector calculus for the standard gauge | experiments_e2::tests::e_and_b_agrees_with_gauge_b |
| E2 (phase 2): Maxwell + Bianchi residuals vanish on the reconstructed field | experiments_e2::tests::maxwell_residuals_vanish_on_reconstructed_field |
| Gamma^mu_(5 nu) = -kappa/2 F^mu_nu exactly (cylinder ansatz) | geodesic::tests::christoffel_5n_matches_lorentz_formula |
| Geodesic spatial accel = kappa*u5*F^i_j u^j (Lorentz force; charge = 5th momentum, weak-field limit) | geodesic::tests::geodesic_accel_is_lorentz_force |
| p5 = u5 + kappa A_mu u^mu conserved along integrated geodesics to 1e-8 | geodesic::tests::p5_conserved_along_geodesic |
| Real NVMe SMART JSON parses: temperature Celsius, health flags correct | ks-acquisition smart::tests::parses_real_apple_nvme_fixture |

## Numerical caveats (validated, not hidden)

- One-sided boundary stencils carry O(dx) error: residuals involving second
  derivatives (inhomogeneous Maxwell, wave operator) evaluate on the 2-point
  interior only. Bianchi residuals (first derivatives only) run on the full
  lattice.
- The plane-wave tests are exact because sin/cos sampled at uniform spacing
  have exact central-difference derivatives.
- Geodesic force-law residual = the O(kappa^2 |A|^2) gradient terms of the
  ansatz (real physics, shrinks as kappa|A| -> 0; the test tolerance 5% and
  the state are documented). The metric degenerates when kappa^2 |A|^2 ~ 1:
  stay in the weak-field region (`Metric5::weak_field_ok`).
- smartctl NVMe JSON reports temperature directly in Celsius (verified
  against this machine's Apple NVMe; the human-readable output agrees).

## What is NOT claimed

- The 5D metric built from observed A_mu is a mathematical embedding, not
  an independent extra-dimensional measurement.
- No compactification-radius prediction, no new-force observable exists in
  this code yet.
- Current HDD interfaces expose error counters and telemetry, not raw
  read-head waveforms; nothing in the acquisition crate treats SMART
  attributes as field samples.