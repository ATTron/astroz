#!/usr/bin/env python3
"""SGP4 orbit propagation with astroz Python bindings."""

import time
import numpy as np
from astroz import Tle, Sgp4

ISS_TLE = """1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995
2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123"""


def main():
    tle = Tle(ISS_TLE)
    sgp4 = Sgp4(tle)

    # Single propagation
    pos, vel = sgp4.propagate(30.0)
    print(f"Single: pos={pos[0]:.1f}, {pos[1]:.1f}, {pos[2]:.1f} km")

    # Batch propagation (1 day at 1-second intervals)
    n = 86400
    times = np.arange(n, dtype=np.float64) / 60.0
    positions = np.empty((n, 3), dtype=np.float64)
    velocities = np.empty((n, 3), dtype=np.float64)

    start = time.perf_counter()
    sgp4.propagate_into(times, positions, velocities)
    elapsed = time.perf_counter() - start

    print(
        f"Batch:  {n:,} points in {elapsed * 1000:.1f} ms ({elapsed / n * 1e6:.2f} μs/prop)"
    )

    # Verify altitude
    altitudes = np.linalg.norm(positions, axis=1) - 6371.0
    print(f"ISS altitude: {altitudes.min():.0f}-{altitudes.max():.0f} km")


if __name__ == "__main__":
    main()
