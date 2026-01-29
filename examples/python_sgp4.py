#!/usr/bin/env python3
"""SGP4 orbit propagation with astroz Python bindings."""

import time
import numpy as np
from astroz import Tle, propagate, Constellation

ISS_TLE = """ISS (ZARYA)
1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995
2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123"""


def main():
    # Parse TLE to inspect properties
    tle = Tle(ISS_TLE)
    print(f"Satellite: {tle.satellite_number}")
    print(f"Inclination: {tle.inclination:.2f}°")
    print(f"Eccentricity: {tle.eccentricity:.6f}")

    # Single propagation
    pos = propagate(ISS_TLE, [30.0])
    print(
        f"Position at t=30min: {pos[0, 0, 0]:.1f}, {pos[0, 0, 1]:.1f}, {pos[0, 0, 2]:.1f} km"
    )

    # Batch propagation (1 day at 1-minute intervals)
    n = 1440
    times = np.arange(n, dtype=np.float64)

    # Pre-parse for repeated calls
    c = Constellation(ISS_TLE)

    start = time.perf_counter()
    positions = propagate(c, times)
    elapsed = time.perf_counter() - start

    print(
        f"Batch: {n:,} points in {elapsed * 1000:.1f} ms ({elapsed / n * 1e6:.2f} μs/prop)"
    )

    # Verify altitude (positions shape: n_times, n_sats, 3)
    altitudes = np.linalg.norm(positions[:, 0, :], axis=1) - 6371.0
    print(f"ISS altitude: {altitudes.min():.0f}-{altitudes.max():.0f} km")


if __name__ == "__main__":
    main()
