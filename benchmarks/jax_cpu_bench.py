#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = ["astrojax>=0.3.1", "jax>=0.9.0"]
# ///

import time

import jax
import jax.numpy as jnp

from astrojax import set_dtype
from astrojax.sgp4 import create_sgp4_propagator

# Force CPU
jax.config.update("jax_default_device", jax.devices("cpu")[0])

# Use float64 precision (must be set before any JIT compilation)
set_dtype(jnp.float64)

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
    params, propagate_fn = create_sgp4_propagator(
        ISS_TLE_LINE1, ISS_TLE_LINE2, gravity="wgs72"
    )
    batch_propagate = jax.jit(jax.vmap(propagate_fn))

    # Warmup: 100 scalar propagations
    for i in range(100):
        r, v = propagate_fn(jnp.float64(i))
        r.block_until_ready()

    # Pre-compile all 5 scenario shapes (JAX caches by input shape)
    for _, points, step in SCENARIOS:
        times = jnp.arange(points, dtype=jnp.float64) * step
        r, v = batch_propagate(times)
        r.block_until_ready()

    device = jax.devices("cpu")[0]
    print("\nJAX CPU SGP4 Benchmark")
    print("=" * 50)
    print(f"Device: {device.device_kind}")
    print(f"JAX version: {jax.__version__}")
    print(f"Precision: float64")

    print("\n--- Vectorized Propagation (jit+vmap) ---")

    total_props_per_sec = 0.0
    for name, points, step in SCENARIOS:
        times = jnp.arange(points, dtype=jnp.float64) * step

        total_ns = 0
        for _ in range(ITERATIONS):
            start = time.perf_counter_ns()
            r, v = batch_propagate(times)
            r.block_until_ready()
            total_ns += time.perf_counter_ns() - start

        avg_ms = total_ns / ITERATIONS / 1_000_000
        props_per_sec = points / (avg_ms / 1000)
        total_props_per_sec += props_per_sec
        print(f"{name:<25} {avg_ms:>10.3f} ms  ({props_per_sec:.2f} prop/s)")

    print(f"{'Average':<25} {total_props_per_sec / len(SCENARIOS):>17.2f} prop/s")

    # Pre-compile large workload shapes
    for _, points, step in MT_SCENARIOS:
        times = jnp.arange(points, dtype=jnp.float64) * step
        r, v = batch_propagate(times)
        r.block_until_ready()

    print("\n--- Vectorized Propagation - Large Workloads (jit+vmap) ---")

    total_props_per_sec = 0.0
    for name, points, step in MT_SCENARIOS:
        times = jnp.arange(points, dtype=jnp.float64) * step

        total_ns = 0
        for _ in range(ITERATIONS):
            start = time.perf_counter_ns()
            r, v = batch_propagate(times)
            r.block_until_ready()
            total_ns += time.perf_counter_ns() - start

        avg_ms = total_ns / ITERATIONS / 1_000_000
        props_per_sec = points / (avg_ms / 1000)
        total_props_per_sec += props_per_sec
        print(f"{name:<25} {avg_ms:>10.3f} ms  ({props_per_sec:.2f} prop/s)")

    print(f"{'Average':<25} {total_props_per_sec / len(MT_SCENARIOS):>17.2f} prop/s")
    print()


if __name__ == "__main__":
    main()
