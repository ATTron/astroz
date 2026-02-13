#!/usr/bin/env rust-script
//! ```cargo
//! [dependencies]
//! sgp4 = "2"
//! rayon = "1"
//!
//! [profile.dev]
//! opt-level = 3
//! ```

use std::hint::black_box;
use std::time::Instant;
use rayon::prelude::*;

const ISS_TLE_LINE1: &str = "1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995";
const ISS_TLE_LINE2: &str = "2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123";

const SCENARIOS: &[(&str, usize, f64)] = &[
    ("1 day (minute)", 1440, 1.0),
    ("1 week (minute)", 10080, 1.0),
    ("2 weeks (minute)", 20160, 1.0),
    ("2 weeks (second)", 1209600, 1.0 / 60.0),
    ("1 month (minute)", 43200, 1.0),
];

const ITERATIONS: u32 = 10;

fn main() {
    let elements =
        sgp4::Elements::from_tle(None, ISS_TLE_LINE1.as_bytes(), ISS_TLE_LINE2.as_bytes()).unwrap();
    let constants = sgp4::Constants::from_elements(&elements).unwrap();

    // Warmup
    for i in 0..100 {
        let _ = constants.propagate(sgp4::MinutesSinceEpoch(i as f64));
    }

    println!("\nRust sgp4 Benchmark");
    println!("==================================================");
    println!("Rayon Threads: {}", rayon::current_num_threads());
    println!("\n--- Sequential Propagation ---");

    let mut total_props_per_sec = 0.0;
    for (name, points, step) in SCENARIOS {
        let times: Vec<f64> = (0..*points).map(|i| i as f64 * step).collect();

        let mut total_ns = 0u128;
        for _ in 0..ITERATIONS {
            let start = Instant::now();
            for &t in &times {
                let _ = black_box(constants.propagate(sgp4::MinutesSinceEpoch(t)));
            }
            total_ns += start.elapsed().as_nanos();
        }

        let avg_ms = total_ns as f64 / ITERATIONS as f64 / 1_000_000.0;
        let props_per_sec = *points as f64 / (avg_ms / 1000.0);
        total_props_per_sec += props_per_sec;
        println!(
            "{:<25} {:>10.3} ms  ({:.2} prop/s)",
            name, avg_ms, props_per_sec
        );
    }
    println!("{:<25} {:>17.2} prop/s", "Average", total_props_per_sec / SCENARIOS.len() as f64);

    println!("\n--- Rayon Parallel Propagation ---");

    let mut total_props_per_sec = 0.0;
    for (name, points, step) in SCENARIOS {
        let times: Vec<f64> = (0..*points).map(|i| i as f64 * step).collect();

        let mut total_ns = 0u128;
        for _ in 0..ITERATIONS {
            let start = Instant::now();
            times.par_iter().for_each(|&t| {
                let _ = black_box(constants.propagate(sgp4::MinutesSinceEpoch(t)));
            });
            total_ns += start.elapsed().as_nanos();
        }

        let avg_ms = total_ns as f64 / ITERATIONS as f64 / 1_000_000.0;
        let props_per_sec = *points as f64 / (avg_ms / 1000.0);
        total_props_per_sec += props_per_sec;
        println!(
            "{:<25} {:>10.3} ms  ({:.2} prop/s)",
            name, avg_ms, props_per_sec
        );
    }
    println!("{:<25} {:>17.2} prop/s", "Average", total_props_per_sec / SCENARIOS.len() as f64);
    println!();
}
