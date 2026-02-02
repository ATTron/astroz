#!/usr/bin/env python3
"""SGP4 orbit propagation with astroz Python bindings (python-sgp4 compatible API)."""

import time
import numpy as np
from astroz.api import Satrec, SatrecArray, WGS72

# ISS TLE lines
LINE1 = "1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995"
LINE2 = "2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123"


def main():
    # Create satellite from TLE (python-sgp4 compatible syntax)
    sat = Satrec.twoline2rv(LINE1, LINE2, WGS72)

    print(f"Satellite: {sat.satnum}")
    print(f"Inclination: {np.degrees(sat.inclo):.2f}°")
    print(f"Eccentricity: {sat.ecco:.6f}")

    # Single propagation at epoch + 30 minutes
    jd, fr = sat.jdsatepoch, sat.jdsatepochF + 30.0 / 1440.0
    error, pos, vel = sat.sgp4(jd, fr)
    print(f"Position at t=30min: {pos[0]:.1f}, {pos[1]:.1f}, {pos[2]:.1f} km")

    # Batch propagation (1 day at 1-minute intervals)
    n = 1440
    sat_array = SatrecArray([sat])

    # Create time arrays (JD format)
    base_jd, base_fr = sat.jdsatepoch, sat.jdsatepochF
    jd_arr = np.full(n, base_jd, dtype=np.float64)
    fr_arr = base_fr + np.arange(n, dtype=np.float64) / 1440.0

    start = time.perf_counter()
    errors, positions, _ = sat_array.sgp4(jd_arr, fr_arr, velocities=False)
    elapsed = time.perf_counter() - start

    print(
        f"Batch: {n:,} points in {elapsed * 1000:.1f} ms ({elapsed / n * 1e6:.2f} μs/prop)"
    )

    # Verify altitude (positions shape: n_sats, n_times, 3)
    altitudes = np.linalg.norm(positions[0, :, :], axis=1) - 6371.0
    print(f"ISS altitude: {altitudes.min():.0f}-{altitudes.max():.0f} km")


if __name__ == "__main__":
    main()
