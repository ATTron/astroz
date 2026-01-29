"""Conjunction screening example using astroz.

Demonstrates: load TLEs, screen for close approaches.
"""

import sys
import time
import numpy as np
from astroz import Constellation, screen


def main():
    if len(sys.argv) < 2:
        print("Usage: python conjunction_screening.py <source>")
        print("  source: file path, CelesTrak group (starlink, all), or NORAD ID")
        print("\nExamples:")
        print("  python conjunction_screening.py starlink")
        print("  python conjunction_screening.py satellites.tle")
        return

    source = sys.argv[1]
    threshold_km = 10.0
    times = np.arange(1440, dtype=np.float64)  # 1 day at 1-min resolution

    # Load constellation
    print(f"Loading {source}...")
    t0 = time.perf_counter()
    c = (
        Constellation(norad_id=int(source))
        if source.isdigit()
        else Constellation(source)
    )
    print(f"Loaded {c.num_satellites} satellites in {time.perf_counter() - t0:.3f}s")

    # All-vs-all screening
    print(f"Screening for conjunctions within {threshold_km} km...")
    t0 = time.perf_counter()
    pairs, t_indices = screen(c, times, threshold=threshold_km)
    elapsed = time.perf_counter() - t0
    print(f"Found {len(pairs)} events in {elapsed:.3f}s")

    # Report results
    if len(pairs) == 0:
        print("No conjunctions found.")
        return

    # Group by unique pairs
    from collections import defaultdict

    pair_times = defaultdict(list)
    for (i, j), t in zip(pairs, t_indices):
        pair_times[(min(i, j), max(i, j))].append(t)

    print(f"\n{len(pair_times)} unique satellite pairs with close approaches:")
    for (sat_i, sat_j), times_list in sorted(
        pair_times.items(), key=lambda x: len(x[1]), reverse=True
    )[:20]:
        print(f"  sat {sat_i:>5d} & {sat_j:>5d}: {len(times_list)} events")

    # Single-target screening (faster for tracking one object)
    print("\nSingle-target screening (satellite 0)...")
    t0 = time.perf_counter()
    min_dists, min_t_indices = screen(
        c, np.arange(1440, dtype=np.float64), threshold=50.0, target=0
    )
    elapsed = time.perf_counter() - t0

    close = np.where(np.array(min_dists) < threshold_km)[0]
    print(
        f"Found {len(close)} satellites within {threshold_km} km of satellite 0 ({elapsed:.3f}s)"
    )


if __name__ == "__main__":
    main()
