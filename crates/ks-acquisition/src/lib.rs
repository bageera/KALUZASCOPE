//! ks-acquisition: sensor reads for the experiments. The HDD/SSD path
//! (smart.rs) parses read-only smartctl output; device-health counters are
//! never treated as field samples — provenance types make conflation
//! unrepresentable.

pub mod smart;

#[cfg(test)]
mod live;

use ks_core::lattice::Lattice;

/// What a measurement actually is.
#[derive(Debug, Clone)]
pub enum Provenance {
    /// Known analytic solution, evaluated on the lattice.
    Simulated { model: String },
    /// Read from real hardware (HDD telemetry, sensor probe, ...).
    Measured { device: String },
    /// Reconstructed from other measured quantities (e.g. potentials from
    /// B-field samples). Carries the source so re-derivation is possible.
    Reconstructed { source: String },
}

#[derive(Debug, Clone)]
pub struct Measurement {
    pub provenance: Provenance,
    /// Sampled field snapshot [t][x][y][z].
    pub field: Lattice,
    pub dx: [f64; 4],
}

impl Measurement {
    pub fn simulated(model: &str, field: Lattice, dx: [f64; 4]) -> Self {
        Measurement {
            provenance: Provenance::Simulated { model: model.to_string() },
            field,
            dx,
        }
    }

    pub fn measured(device: &str, field: Lattice, dx: [f64; 4]) -> Self {
        Measurement {
            provenance: Provenance::Measured { device: device.to_string() },
            field,
            dx,
        }
    }

    pub fn reconstructed(source: &str, field: Lattice, dx: [f64; 4]) -> Self {
        Measurement {
            provenance: Provenance::Reconstructed { source: source.to_string() },
            field,
            dx,
        }
    }

    pub fn is_measured(&self) -> bool {
        matches!(self.provenance, Provenance::Measured { .. })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ks_core::lattice::zeros;

    #[test]
    fn provenance_is_tracked() {
        let m = Measurement::simulated("plane-wave", zeros([2, 2, 2, 2], 0.0), [1.0; 4]);
        assert!(!m.is_measured());
        let m2 = Measurement::measured("/dev/disk0", zeros([2, 2, 2, 2], 0.0), [1.0; 4]);
        assert!(m2.is_measured());
        let m3 = Measurement::reconstructed("curl", zeros([2, 2, 2, 2], 0.0), [1.0; 4]);
        assert!(!m3.is_measured());
    }
}