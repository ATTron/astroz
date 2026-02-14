#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10,<3.13"
# dependencies = ["astroz>=0.8", "numpy"]
# ///

import time

import numpy as np
from astroz.api import Satrec, WGS72

ISS_TLE_LINE1 = "1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995"
ISS_TLE_LINE2 = "2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123"

SCENARIOS = [
    ("1 day (minute)", 1440, 1.0),
    ("1 week (minute)", 10080, 1.0),
    ("2 weeks (minute)", 20160, 1.0),
    ("2 weeks (second)", 1209600, 1.0 / 60.0),
    ("1 month (minute)", 43200, 1.0),
]

MT_SCENARIOS = [
    ("1 month (second)", 2592000, 1.0 / 60.0),
    ("3 months (second)", 7776000, 1.0 / 60.0),
    ("1 year (minute)", 525600, 1.0),
    ("1 year (second)", 31536000, 1.0 / 60.0),
]

ITERATIONS = 10


def main():
    sat = Satrec.twoline2rv(ISS_TLE_LINE1, ISS_TLE_LINE2, WGS72)

    jd_epoch = sat.jdsatepoch
    fr_epoch = sat.jdsatepochF

    # Warmup
    for i in range(100):
        sat.sgp4(jd_epoch, fr_epoch + i / 1440.0)

    print("\nPython astroz Benchmark")
    print("=" * 50)
    print("\n--- Sequential Propagation ---")

    total_props_per_sec = 0.0
    for name, points, step in SCENARIOS:
        times = [i * step for i in range(points)]

        total_ns = 0
        for _ in range(ITERATIONS):
            start = time.perf_counter_ns()
            for t in times:
                sat.sgp4(jd_epoch, fr_epoch + t / 1440.0)
            total_ns += time.perf_counter_ns() - start

        avg_ms = total_ns / ITERATIONS / 1_000_000
        props_per_sec = points / (avg_ms / 1000)
        total_props_per_sec += props_per_sec
        print(f"{name:<25} {avg_ms:>10.3f} ms  ({props_per_sec:.2f} prop/s)")

    print(f"{'Average':<25} {total_props_per_sec / len(SCENARIOS):>17.2f} prop/s")

    print("\n--- SIMD Batch Propagation ---")

    total_props_per_sec = 0.0
    for name, points, step in SCENARIOS:
        times = np.arange(points, dtype=np.float64) * step
        jd_arr = np.full(points, jd_epoch)
        fr_arr = fr_epoch + times / 1440.0

        total_ns = 0
        for _ in range(ITERATIONS):
            start = time.perf_counter_ns()
            sat.sgp4_array(jd_arr, fr_arr)
            total_ns += time.perf_counter_ns() - start

        avg_ms = total_ns / ITERATIONS / 1_000_000
        props_per_sec = points / (avg_ms / 1000)
        total_props_per_sec += props_per_sec
        print(f"{name:<25} {avg_ms:>10.3f} ms  ({props_per_sec:.2f} prop/s)")

    print(f"{'Average':<25} {total_props_per_sec / len(SCENARIOS):>17.2f} prop/s")

    print("\n--- SIMD Batch Propagation (Large Workloads) ---")

    total_props_per_sec = 0.0
    for name, points, step in MT_SCENARIOS:
        times = np.arange(points, dtype=np.float64) * step
        jd_arr = np.full(points, jd_epoch)
        fr_arr = fr_epoch + times / 1440.0

        total_ns = 0
        for _ in range(ITERATIONS):
            start = time.perf_counter_ns()
            sat.sgp4_array(jd_arr, fr_arr)
            total_ns += time.perf_counter_ns() - start

        avg_ms = total_ns / ITERATIONS / 1_000_000
        props_per_sec = points / (avg_ms / 1000)
        total_props_per_sec += props_per_sec
        print(f"{name:<25} {avg_ms:>10.3f} ms  ({props_per_sec:.2f} prop/s)")

    print(f"{'Average':<25} {total_props_per_sec / len(MT_SCENARIOS):>17.2f} prop/s")
    print()


if __name__ == "__main__":
    main()
