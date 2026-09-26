#!/usr/bin/env rust-script
//! SGP4 throughput benchmark for the Rust `sgp4` crate.
//!
//! ```cargo
//! [dependencies]
//! sgp4 = "2"
//! rayon = "1"
//!
//! [profile.dev]
//! opt-level = 3
//! ```

use rayon::prelude::*;
use std::hint::black_box;
use std::time::Instant;

const LINE1: &str = "1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995";
const LINE2: &str = "2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123";

const ITERATIONS: u32 = 10;
const WARMUP: usize = 100;

struct Scenario {
    name: &'static str,
    points: usize,
    step: f64,
}

const STANDARD: &[Scenario] = &[
    Scenario { name: "1 day (minute)", points: 1440, step: 1.0 },
    Scenario { name: "1 week (minute)", points: 10080, step: 1.0 },
    Scenario { name: "2 weeks (minute)", points: 20160, step: 1.0 },
    Scenario { name: "2 weeks (second)", points: 1209600, step: 1.0 / 60.0 },
    Scenario { name: "1 month (minute)", points: 43200, step: 1.0 },
];

const LARGE: &[Scenario] = &[
    Scenario { name: "1 month (second)", points: 2592000, step: 1.0 / 60.0 },
    Scenario { name: "3 months (second)", points: 7776000, step: 1.0 / 60.0 },
    Scenario { name: "1 year (minute)", points: 525600, step: 1.0 },
    Scenario { name: "1 year (second)", points: 31536000, step: 1.0 / 60.0 },
];

fn make_times(points: usize, step: f64) -> Vec<f64> {
    (0..points).map(|i| i as f64 * step).collect()
}

fn propagate_seq(c: &sgp4::Constants, times: &[f64]) {
    for &t in times {
        black_box(c.propagate(sgp4::MinutesSinceEpoch(t)).ok());
    }
}

fn propagate_par(c: &sgp4::Constants, times: &[f64]) {
    times.par_iter().for_each(|&t| {
        black_box(c.propagate(sgp4::MinutesSinceEpoch(t)).ok());
    });
}

fn section(title: &str, scenarios: &[Scenario], c: &sgp4::Constants, run: fn(&sgp4::Constants, &[f64])) {
    println!("\n--- {} ---", title);
    let mut rates = Vec::with_capacity(scenarios.len());
    for s in scenarios {
        let times = make_times(s.points, s.step);
        let start = Instant::now();
        for _ in 0..ITERATIONS {
            run(c, &times);
        }
        let avg = start.elapsed().as_secs_f64() / ITERATIONS as f64;
        let rate = s.points as f64 / avg;
        rates.push(rate);
        println!("{:<25} {:>10.3} ms  ({:.2} prop/s)", s.name, avg * 1e3, rate);
    }
    let mean = rates.iter().sum::<f64>() / rates.len() as f64;
    println!("{:<25} {:>17.2} prop/s", "Average", mean);
}

fn main() {
    let elements = sgp4::Elements::from_tle(None, LINE1.as_bytes(), LINE2.as_bytes())
        .expect("failed to parse TLE");
    let constants = sgp4::Constants::from_elements(&elements).expect("failed to init SGP4");

    let warm = make_times(WARMUP, 1.0);
    propagate_seq(&constants, &warm);
    propagate_par(&constants, &warm);

    println!("\nRust sgp4 Benchmark");
    println!("{}", "=".repeat(50));
    println!("Rayon Threads: {}", rayon::current_num_threads());

    section("Sequential Propagation", STANDARD, &constants, propagate_seq);
    section("Rayon Parallel Propagation", STANDARD, &constants, propagate_par);
    section("Rayon Parallel (Large Workloads)", LARGE, &constants, propagate_par);
}
