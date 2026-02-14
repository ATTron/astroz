#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = ["sgp4>=2.22"]
# ///

import multiprocessing
import time

from sgp4.api import Satrec, WGS72

ISS_TLE_LINE1 = "1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995"
ISS_TLE_LINE2 = "2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123"

SCENARIOS = [
    ("1 day (minute)", 1440, 1.0),
    ("1 week (minute)", 10080, 1.0),
    ("2 weeks (minute)", 20160, 1.0),
    ("2 weeks (second)", 1209600, 1.0 / 60.0),
    ("1 month (minute)", 43200, 1.0),
]

ITERATIONS = 10

# Worker-local Satrec for multiprocessing
_worker_sat = None


def _init_worker(line1, line2):
    global _worker_sat
    _worker_sat = Satrec.twoline2rv(line1, line2, WGS72)


def _propagate_chunk(times_chunk):
    for t in times_chunk:
        _worker_sat.sgp4_tsince(t)


def main():
    sat = Satrec.twoline2rv(ISS_TLE_LINE1, ISS_TLE_LINE2, WGS72)

    # Warmup
    for i in range(100):
        sat.sgp4_tsince(float(i))

    num_workers = multiprocessing.cpu_count()

    print("\nPython sgp4 Benchmark")
    print("=" * 50)
    print(f"Workers: {num_workers}")
    print("\n--- Sequential Propagation ---")

    total_props_per_sec = 0.0
    for name, points, step in SCENARIOS:
        times = [i * step for i in range(points)]

        total_ns = 0
        for _ in range(ITERATIONS):
            start = time.perf_counter_ns()
            for t in times:
                sat.sgp4_tsince(t)
            total_ns += time.perf_counter_ns() - start

        avg_ms = total_ns / ITERATIONS / 1_000_000
        props_per_sec = points / (avg_ms / 1000)
        total_props_per_sec += props_per_sec
        print(f"{name:<25} {avg_ms:>10.3f} ms  ({props_per_sec:.2f} prop/s)")

    print(f"{'Average':<25} {total_props_per_sec / len(SCENARIOS):>17.2f} prop/s")

    print("\n--- Multiprocessing Parallel Propagation ---")

    total_props_per_sec = 0.0
    with multiprocessing.Pool(
        processes=num_workers,
        initializer=_init_worker,
        initargs=(ISS_TLE_LINE1, ISS_TLE_LINE2),
    ) as pool:
        for name, points, step in SCENARIOS:
            times = [i * step for i in range(points)]

            # Split times into chunks, one per worker
            chunk_size = (len(times) + num_workers - 1) // num_workers
            chunks = [
                times[i : i + chunk_size] for i in range(0, len(times), chunk_size)
            ]

            total_ns = 0
            for _ in range(ITERATIONS):
                start = time.perf_counter_ns()
                pool.map(_propagate_chunk, chunks)
                total_ns += time.perf_counter_ns() - start

            avg_ms = total_ns / ITERATIONS / 1_000_000
            props_per_sec = points / (avg_ms / 1000)
            total_props_per_sec += props_per_sec
            print(f"{name:<25} {avg_ms:>10.3f} ms  ({props_per_sec:.2f} prop/s)")

    print(f"{'Average':<25} {total_props_per_sec / len(SCENARIOS):>17.2f} prop/s")

    MT_SCENARIOS = [
        ("1 month (second)", 2592000, 1.0 / 60.0),
        ("3 months (second)", 7776000, 1.0 / 60.0),
        ("1 year (minute)", 525600, 1.0),
        ("1 year (second)", 31536000, 1.0 / 60.0),
    ]

    print("\n--- Multiprocessing Parallel (Large Workloads) ---")

    total_props_per_sec = 0.0
    with multiprocessing.Pool(
        processes=num_workers,
        initializer=_init_worker,
        initargs=(ISS_TLE_LINE1, ISS_TLE_LINE2),
    ) as pool:
        for name, points, step in MT_SCENARIOS:
            times = [i * step for i in range(points)]

            chunk_size = (len(times) + num_workers - 1) // num_workers
            chunks = [
                times[i : i + chunk_size] for i in range(0, len(times), chunk_size)
            ]

            total_ns = 0
            for _ in range(ITERATIONS):
                start = time.perf_counter_ns()
                pool.map(_propagate_chunk, chunks)
                total_ns += time.perf_counter_ns() - start

            avg_ms = total_ns / ITERATIONS / 1_000_000
            props_per_sec = points / (avg_ms / 1000)
            total_props_per_sec += props_per_sec
            print(f"{name:<25} {avg_ms:>10.3f} ms  ({props_per_sec:.2f} prop/s)")

    print(f"{'Average':<25} {total_props_per_sec / len(MT_SCENARIOS):>17.2f} prop/s")
    print()


if __name__ == "__main__":
    main()
