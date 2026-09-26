# /// script
# requires-python = ">=3.10"
# dependencies = ["astrojax>=0.3.1", "jax>=0.9.0"]
# ///
"""SGP4 throughput benchmark for astrojax on the JAX CPU backend."""

import os

os.environ["JAX_PLATFORMS"] = "cpu"

import time

import astrojax
import astrojax.sgp4
import jax
import jax.numpy as jnp

astrojax.set_dtype(jnp.float64)

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


def make_times(points, step):
    return jnp.arange(points, dtype=jnp.float64) * step


def block(result):
    r, v = result
    r.block_until_ready()
    v.block_until_ready()


def report_row(name, points, elapsed):
    avg = elapsed / ITERATIONS
    rate = points / avg
    print(f"{name:<25} {avg * 1e3:>10.3f} ms  ({rate:.2f} prop/s)", flush=True)
    return rate


def report_average(rates):
    print(f"{'Average':<25} {sum(rates) / len(rates):>17.2f} prop/s", flush=True)


def bench(batched, title, scenarios):
    print(f"\n--- {title} ---")
    rates = []
    for name, points, step in scenarios:
        times = make_times(points, step)
        block(batched(times))  # compile for this shape
        start = time.perf_counter()
        for _ in range(ITERATIONS):
            block(batched(times))
        rates.append(report_row(name, points, time.perf_counter() - start))
        del times
    report_average(rates)


def main():
    _, propagate_fn = astrojax.sgp4.create_sgp4_propagator(LINE1, LINE2, gravity="wgs72")
    batched = jax.jit(jax.vmap(propagate_fn))
    block(batched(make_times(WARMUP, 1.0)))

    print("\nJAX CPU SGP4 Benchmark")
    print("=" * 50)
    print(f"Device: {jax.devices('cpu')[0].device_kind}")
    print(f"JAX Version: {jax.__version__}")
    print("Precision: float64")

    bench(batched, "Vectorized Propagation (jit+vmap)", STANDARD)
    bench(batched, "Vectorized Propagation - Large Workloads (jit+vmap)", LARGE)


if __name__ == "__main__":
    main()
