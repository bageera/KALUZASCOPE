//! ks-acquisition: sensor reads for the experiments. v0.1.0 ships the
//! simulated feed only; the HDD smartctl adapter lands when hardware is in
//! the loop (phase 3 of the roadmap). Measured vs. reconstructed data are
//! NEVER mixed: `Measurement` carries its own provenance.

use ks_core::lattice::Lattice;

/// What a measurement actually is. The distinction matters: an HDD SMART
/// attribute is a device-health counter, not a field sample. The enum makes
/// conflation unrepresentable.
#[derive(Debug, Clone)]
pub enum Provenance {
    /// Known analytic solution, evaluated on the lattice.
    Simulated { model: String },
    /// Read from real hardware (HDD telemetry, sensor probe, ...).
    Measured { device: String },
}

#[derive(Debug, Clone)]
pub struct Measurement {
    pub provenance: Provenance,
    /// Sampled field snapshot [t][x][y][z] (the only channel the pipeline
    /// consumes; the acquisition layer is responsible for its meaning).
    pub field: Lattice,
    pub dx: [f64; 4],
}

impl Measurement {
    /// Simulated measurement from an analytic model.
    pub fn simulated(model: &str, field: Lattice, dx: [f64; 4]) -> Self {
        Measurement {
            provenance: Provenance::Simulated { model: model.to_string() },
            field,
            dx,
        }
    }
}

/// The planned phase-3 adapter. Stub API so the experiment harness can be
/// written against it before hardware exists; the contract: read-only,
/// no firmware writes, provenance always Measured.
pub mod hdd {
    /// One smartctl reading. Only the fields the research can legitimately
    /// use; drives do NOT expose raw head waveforms through any standard
    /// interface (documented in docs/experimental-methodology.md).
    #[derive(Debug, Clone)]
    pub struct SmartSample {
        pub device: String,
        pub temperature_c: f64,
        pub reallocated_sectors: u64,
        pub read_error_rate: f64,
    }

    /// Parse `smartctl -a` output. Not wired to a real device yet: the
    /// parser lands with the first supported drive on the bench.
    pub fn parse_smartctl(_output: &str) -> Option<SmartSample> {
        None // ponytail: stub until phase-3 hardware; parser when real output exists
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ks_core::lattice::zeros;

    #[test]
    fn simulated_measurement_carries_provenance() {
        let m = Measurement::simulated("plane-wave", zeros([2, 2, 2, 2], 0.0), [1.0; 4]);
        match m.provenance {
            Provenance::Simulated { model } => assert_eq!(model, "plane-wave"),
            Provenance::Measured { .. } => panic!("wrong provenance"),
        }
    }
}