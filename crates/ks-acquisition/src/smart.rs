//! HDD/SSD SMART telemetry (phase 3): read-only parse of `smartctl -a -j`
//! output into a typed sample. Measured data always carries `Provenance::Measured`.
//! Device-health counters are NOT field samples — see docs/experimental-methodology.md.

//! The parser consumes the JSON emitted by smartctl 7.x. The fixture in
//! `experiments/fixtures/apple_nvme_smart.json` is a real Apple NVMe read.

use serde_json::Value;

/// NVMe SMART/Health log (the fields this project legitimately consumes).
#[derive(Debug, Clone, PartialEq)]
pub struct SmartSample {
    pub device: String,
    pub model: String,
    pub critical_warning_raw: u64,
    pub temperature_c: f64,
    pub available_spare_pct: f64,
    pub available_spare_threshold_pct: f64,
    pub percentage_used: f64,
    pub media_errors: u64,
    pub power_on_hours: f64,
    pub unsafe_shutdowns: u64,
}

impl SmartSample {
    /// True when the drive reports no critical condition and enough spare.
    pub fn healthy(&self) -> bool {
        self.critical_warning_raw == 0
            && self.available_spare_pct >= self.available_spare_threshold_pct
            && self.media_errors == 0
            && self.temperature_c < 70.0
    }
}

/// Parse a smartctl -j document. Returns None if the NVMe health log is
/// missing (e.g. ATA drives have a different schema — not supported yet).
pub fn parse_nvme_smart_json(text: &str, device: &str) -> Option<SmartSample> {
    let doc: Value = parse_json(text)?;
    // smartctl exits nonzero on some read errors but still emits JSON; the
    // health log is what we need.
    let log = doc.get("nvme_smart_health_information_log")?;
    let model = doc
        .pointer("/model_name")
        .and_then(|v| v.as_str())
        .unwrap_or("unknown")
        .to_string();
    Some(SmartSample {
        device: device.to_string(),
        model,
        critical_warning_raw: log.get("critical_warning")?.as_u64()?,
        temperature_c: log.get("temperature")?.as_f64()?,
        available_spare_pct: log.get("available_spare")?.as_f64()?,
        available_spare_threshold_pct: log.get("available_spare_threshold")?.as_f64()?,
        percentage_used: log.get("percentage_used")?.as_f64().unwrap_or(0.0),
        media_errors: log.get("media_errors")?.as_u64()?,
        power_on_hours: log.get("power_on_hours")?.as_f64().unwrap_or(0.0),
        unsafe_shutdowns: log.get("unsafe_shutdowns")?.as_u64().unwrap_or(0),
    })
}

fn parse_json(text: &str) -> Option<Value> {
    serde_json::from_str(text).ok()
}

/// smartctl NVMe JSON reports temperature directly in Celsius (verified
/// against a real Apple NVMe: log says 40, human-readable output says
/// "40 Celsius"). No conversion is applied.
pub fn read_device(device: &str) -> Option<SmartSample> {
    let out = std::process::Command::new("smartctl")
        .arg("-a")
        .arg("-j")
        .arg(device)
        .output()
        .ok()?;
    let text = String::from_utf8_lossy(&out.stdout);
    parse_nvme_smart_json(&text, device)
}

#[cfg(test)]
mod tests {
    use super::*;

    const FIXTURE: &str = include_str!("../../../experiments/fixtures/apple_nvme_smart.json");

    #[test]
    fn parses_real_apple_nvme_fixture() {
        let s = parse_nvme_smart_json(FIXTURE, "/dev/disk0").expect("fixture parse");
        assert_eq!(s.model, "APPLE SSD AP0512Z");
        assert!((s.temperature_c - 40.0).abs() < 0.01, "temp {}", s.temperature_c);
        assert_eq!(s.critical_warning_raw, 0);
        assert_eq!(s.media_errors, 0);
        assert!((s.available_spare_pct - 100.0).abs() < 0.01);
        assert!((s.available_spare_threshold_pct - 99.0).abs() < 0.01);
        assert!(s.healthy());
    }

    #[test]
    fn missing_health_log_is_none() {
        assert!(parse_nvme_smart_json("{\"device\":{}}", "/dev/x").is_none());
        assert!(parse_nvme_smart_json("not json", "/dev/x").is_none());
    }

    #[test]
    fn unhealthy_drive_detected() {
        let doc = r#"{
            "model_name": "TEST DISK",
            "nvme_smart_health_information_log": {
                "critical_warning": 4,
                "temperature": 100,
                "available_spare": 10,
                "available_spare_threshold": 99,
                "percentage_used": 105,
                "media_errors": 17,
                "power_on_hours": 5,
                "unsafe_shutdowns": 1
            }
        }"#;
        let s = parse_nvme_smart_json(doc, "/dev/x").expect("parse");
        assert!(!s.healthy(), "must flag degraded drive");
        assert!((s.temperature_c - 100.0).abs() < 0.2);
    }
}