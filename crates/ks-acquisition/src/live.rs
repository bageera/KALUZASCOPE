// Live phase-3 verification: run the parser against the real NVMe via
// ks-acquisition::smart::read_device. Compile-check + run through the test
// harness (kept as an #[ignore] test so CI never depends on hardware).
//
//   cargo test -p ks-acquisition live_nvme_read -- --ignored --nocapture
#[cfg(test)]
mod live {
    #[test]
    #[ignore]
    fn live_nvme_read() {
        match crate::smart::read_device("/dev/disk0") {
            Some(s) => {
                println!("device:        {}", s.device);
                println!("model:         {}", s.model);
                println!("temperature_c: {}", s.temperature_c);
                println!("spare pct:     {}", s.available_spare_pct);
                println!("pct used:      {}", s.percentage_used);
                println!("media errors:  {}", s.media_errors);
                println!("POH:           {}", s.power_on_hours);
                println!("healthy:       {}", s.healthy());
                assert!(s.temperature_c > 0.0 && s.temperature_c < 120.0);
            }
            None => panic!("smartctl parse returned None on the live drive"),
        }
    }
}