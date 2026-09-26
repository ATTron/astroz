# /// script
# requires-python = ">=3.10"
# dependencies = ["astroz>=0.8", "numpy"]
# ///
"""SGP4 throughput benchmark for astroz through its python-sgp4 compatible API."""

import time

import numpy as np
from astroz.api import WGS72, Satrec

LINE1 = "1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995"
LINE2 = "2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123"

ITERATIONS = 10
WARMUP = 100
MINUTES_PER_DAY = 1440.0

STANDARD = [
    ("1 day (minute)", 1440, 1.0),
    ("1 week (minute)", 10080, 1.0),
    ("2 weeks (minute)", 20160, 1.0),
    ("2 weeks (second)", 1209600, 1.0 / 60.0),
    ("1 month (minute)", 43200, 1.0),
]

LARGE = [
    ("1 month (second)", 2592000, 1.0 / 60.0),
    ("3 months (second)", 7776000, 1.0 / 60.0),
    ("1 year (minute)", 525600, 1.0),
    ("1 year (second)", 31536000, 1.0 / 60.0),
]


def julian_arrays(sat, points, step):
    """Return (jd, fr) arrays for times i*step minutes after the TLE epoch."""
    minutes = np.arange(points, dtype=np.float64) * step
    jd = np.full(points, sat.jdsatepoch, dtype=np.float64)
    fr = sat.jdsatepochF + minutes / MINUTES_PER_DAY
    return jd, fr


def report_row(name, points, elapsed):
    avg = elapsed / ITERATIONS
    rate = points / avg
    print(f"{name:<25} {avg * 1e3:>10.3f} ms  ({rate:.2f} prop/s)", flush=True)
    return rate


def report_average(rates):
    print(f"{'Average':<25} {sum(rates) / len(rates):>17.2f} prop/s", flush=True)


def bench_sequential(sat, scenarios):
    print("\n--- Sequential Propagation ---")
    rates = []
    for name, points, step in scenarios:
        jd, fr = julian_arrays(sat, points, step)
        pairs = list(zip(jd.tolist(), fr.tolist()))
        f = sat.sgp4
        start = time.perf_counter()
        for _ in range(ITERATIONS):
            for j, r in pairs:
                f(j, r)
        rates.append(report_row(name, points, time.perf_counter() - start))
    report_average(rates)


def bench_batch(sat, title, scenarios):
    print(f"\n--- {title} ---")
    rates = []
    for name, points, step in scenarios:
        jd, fr = julian_arrays(sat, points, step)
        start = time.perf_counter()
        for _ in range(ITERATIONS):
            sat.sgp4_array(jd, fr)
        rates.append(report_row(name, points, time.perf_counter() - start))
        del jd, fr
    report_average(rates)


def main():
    sat = Satrec.twoline2rv(LINE1, LINE2, WGS72)

    jd, fr = julian_arrays(sat, WARMUP, 1.0)
    for j, r in zip(jd.tolist(), fr.tolist()):
        sat.sgp4(j, r)
    sat.sgp4_array(jd, fr)

    print("\nPython astroz Benchmark")
    print("=" * 50)

    bench_sequential(sat, STANDARD)
    bench_batch(sat, "SIMD Batch Propagation", STANDARD)
    bench_batch(sat, "SIMD Batch Propagation (Large Workloads)", LARGE)


if __name__ == "__main__":
    main()
