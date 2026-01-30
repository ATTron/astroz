#!/usr/bin/env python3
"""SGP4 benchmark comparing astroz vs python-sgp4."""

import time
import numpy as np

ISS_TLE = """ISS (ZARYA)
1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995
2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123"""

ISS_TLE_LINE1 = "1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995"
ISS_TLE_LINE2 = "2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123"

SCENARIOS = [
    ("1 day (minute)", 1440, 1.0),
    ("1 week (minute)", 10080, 1.0),
    ("2 weeks (minute)", 20160, 1.0),
    ("1 month (minute)", 43200, 1.0),
]

ITERATIONS = 10


def bench_astroz():
    try:
        from astroz import Constellation, propagate
    except ImportError:
        print(
            "astroz not installed... run: zig build python-bindings && pip install -e bindings/python/"
        )
        return None

    # Pre-parse constellation to avoid parsing overhead in benchmark
    c = Constellation(ISS_TLE)

    # Warmup
    for _ in range(10):
        propagate(c, np.arange(100, dtype=np.float64))

    results = {}
    for name, points, step in SCENARIOS:
        times = np.arange(points, dtype=np.float64) * step
        total = 0
        for _ in range(ITERATIONS):
            start = time.perf_counter_ns()
            propagate(c, times)
            total += time.perf_counter_ns() - start
        results[name] = total / ITERATIONS / 1_000_000
    return results


def bench_sgp4():
    from sgp4.api import Satrec, jday

    sat = Satrec.twoline2rv(ISS_TLE_LINE1, ISS_TLE_LINE2)
    jd, fr = jday(2024, 1, 1, 12, 0, 0)
    for i in range(100):  # warmup
        sat.sgp4(jd, fr + i / 1440.0)

    results = {}
    for name, points, step in SCENARIOS:
        jd_arr = np.full(points, jd)
        fr_arr = np.arange(points, dtype=np.float64) * step / 1440.0
        total = 0
        for _ in range(ITERATIONS):
            start = time.perf_counter_ns()
            sat.sgp4_array(jd_arr, fr_arr)
            total += time.perf_counter_ns() - start
        results[name] = total / ITERATIONS / 1_000_000
    return results


def main():
    print("SGP4 Benchmark (ms)\n")
    astroz = bench_astroz()
    sgp4 = bench_sgp4()

    if astroz and sgp4:
        print(f"{'Scenario':<22} {'astroz':>10} {'sgp4':>10} {'speedup':>10}")
        print("-" * 55)
        for name, _, _ in SCENARIOS:
            a, s = astroz[name], sgp4[name]
            print(f"{name:<22} {a:>10.2f} {s:>10.2f} {s / a:>10.2f}x")
    elif sgp4:
        print(f"{'Scenario':<22} {'sgp4':>10}")
        print("-" * 35)
        for name, _, _ in SCENARIOS:
            print(f"{name:<22} {sgp4[name]:>10.2f}")


if __name__ == "__main__":
    main()
