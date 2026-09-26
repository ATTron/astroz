# /// script
# requires-python = ">=3.10"
# dependencies = ["sgp4>=2.22"]
# ///
"""SGP4 throughput benchmark for the reference python-sgp4 package."""

import os
import time
from multiprocessing import Pool

from sgp4.api import WGS72, Satrec

LINE1 = "1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995"
LINE2 = "2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123"

ITERATIONS = 10
WARMUP = 100

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

_worker_sat = None


def make_times(points, step):
    return [i * step for i in range(points)]


def propagate_all(sat, times):
    f = sat.sgp4_tsince
    for t in times:
        f(t)


def init_worker():
    global _worker_sat
    _worker_sat = Satrec.twoline2rv(LINE1, LINE2, WGS72)


def worker_propagate(times):
    propagate_all(_worker_sat, times)
    return len(times)


def split(times, n):
    size = (len(times) + n - 1) // n
    return [times[i : i + size] for i in range(0, len(times), size)]


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
        times = make_times(points, step)
        start = time.perf_counter()
        for _ in range(ITERATIONS):
            propagate_all(sat, times)
        rates.append(report_row(name, points, time.perf_counter() - start))
    report_average(rates)


def bench_pool(pool, workers, title, scenarios):
    print(f"\n--- {title} ---")
    rates = []
    for name, points, step in scenarios:
        chunks = split(make_times(points, step), workers)
        start = time.perf_counter()
        for _ in range(ITERATIONS):
            pool.map(worker_propagate, chunks, chunksize=1)
        rates.append(report_row(name, points, time.perf_counter() - start))
        del chunks
    report_average(rates)


def main():
    workers = os.cpu_count() or 1
    sat = Satrec.twoline2rv(LINE1, LINE2, WGS72)
    propagate_all(sat, make_times(WARMUP, 1.0))

    print("\nPython sgp4 Benchmark")
    print("=" * 50)
    print(f"Workers: {workers}")

    bench_sequential(sat, STANDARD)

    with Pool(workers, initializer=init_worker) as pool:
        pool.map(worker_propagate, split(make_times(WARMUP, 1.0), workers))
        bench_pool(pool, workers, "Multiprocessing Parallel Propagation", STANDARD)
        bench_pool(pool, workers, "Multiprocessing Parallel (Large Workloads)", LARGE)


if __name__ == "__main__":
    main()
