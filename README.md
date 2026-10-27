# KALUZASCOPE

Kaluza–Klein electromagnetic field research: reconstruct electromagnetic
fields using five-dimensional geometry, with every mathematical claim backed
by a reproducible test.

Research principle: **every mathematical claim has a corresponding
reproducible test** (`cargo test` runs the whole proof suite).

## What is here

- `crates/ks-core` — 5D Kaluza–Klein metric machinery:
  - `lattice` — uniform 4-D grid, central finite differences
  - `electromagnetics` — the field tensor F_mn = d_m A_n - d_n A_m, E/B
    extraction, Lorentz invariants (B^2 - E^2, E·B)
  - `metric` — the Kaluza metric G_AB built from (g_mn, A_mu, phi, kappa),
    block determinant identity det G = phi^2 det g, inverse
  - `reduction` — verify that the reduced 5D story reproduces 4D Maxwell:
    homogeneous (Bianchi) and inhomogeneous (d_m F^{mn} = j^n) residuals,
    plus the wave operator. A lightlike plane wave (w = k) passes both
    Maxwell residuals; a deliberately wrong dispersion (w != k) exposes a
    measurable O((w^2-k^2) sin) failure — the residuals detect physics errors.
- `crates/ks-acquisition` — measurement provenance types. HDD SMART/FARM
  telemetry is planned read-only, phase 3 of the roadmap; nothing here
  pretends a disk-health counter is a field sample
- `crates/ks-core` geodesics (phase 5 material): numerical 5D Christoffels,
  RK4 integrator; proven Γ^μ₅ν = −κ/2 F^μ_ν, geodesic accel = κu₅F^iⱼu^ⱼ
  (Lorentz force with the fifth momentum as charge), p₅ conserved.
- `experiments` — phase-2 reconstruction experiment lives in
  `ks-core::experiments_e2` (curl gauge → field tensor → Maxwell residuals).


## Honest scope

The 5D metric construction is a **mathematical embedding** of observed or
simulated EM fields, not an independent measurement of an extra dimension.
Hard drives do not expose raw read-head analog signals through any standard
interface, and no existing hardware senses "the fifth dimension" through
KK compactification. The engineering goal (geometric reconstruction
advantage) is separate from the physics goal (a 5D-model-specific
prediction confirmed by controlled measurement); both are open.

## Usage

```
cargo test          # run the reproducible-proof suite
cargo build --release
```

Rust toolchain: any recent stable. If `cargo` is not on PATH, this repo
was built with rustc 1.99.0 at ~/.local/opt/rust/bin.

## Roadmap

1. Verified 5D tensor/metric library — DONE (v0.1.0)
2. EM simulation & reconstruction experiments — DONE (E2: curl reconstruction,
   Maxwell + Bianchi residuals vanish on reconstructed fields)
3. HDD/SSD SMART telemetry adapter (read-only) — DONE (NVMe JSON parser,
   provenance-tagged; fixture = real Apple NVMe; live-read verified)
4. Real magnetic-field acquisition (calibrated external sensor) — blocked on
   hardware: no calibrated sensor on the bench. The Measurement/Provenance
   types are the contract to plug into.
5. Conventional-EM vs 5D-prediction comparison — partially done (geodesic
   force law + p5 charge conservation proven; compactification-mode
   predictions are the remaining item)