"""Conjunction screening example using astroz.

Demonstrates: load TLEs, propagate, screen for close approaches.
"""

import sys
import time
from collections import defaultdict
import numpy as np
from astroz import (
    load_constellation,
    propagate_constellation,
    min_distances,
    coarse_screen,
)

COARSE_THRESHOLD_KM = 10.0
FINE_THRESHOLD_KM = 5.0
MAX_COARSE_HITS = 10  # Filter out co-orbital pairs


def main():
    if len(sys.argv) < 2:
        print("Usage: python conjunction_screening.py <source>")
        print("  source: file path, CelesTrak group (starlink, all), or NORAD ID")
        print("\nExamples:")
        print("  python conjunction_screening.py starlink")
        print("  python conjunction_screening.py satellites.tle")
        return

    # Load constellation
    t0 = time.perf_counter()
    source = sys.argv[1]
    constellation = (
        load_constellation(norad_id=int(source))
        if source.isdigit()
        else load_constellation(source)
    )
    num_sats = constellation.num_satellites
    print(f"Loaded {num_sats} satellites in {time.perf_counter() - t0:.3f}s")

    # Phase 1: Coarse propagation (1 day at 1-min resolution)
    t0 = time.perf_counter()
    positions = propagate_constellation(
        constellation,
        np.arange(1440, dtype=np.float64),
        output="teme",
        layout="satellite_major",
    )
    print(
        f"Phase 1: Propagation ({num_sats} x 1440) in {time.perf_counter() - t0:.3f}s"
    )

    # Phase 2: Spatial hash screening
    t0 = time.perf_counter()
    raw_pairs, raw_t_indices = coarse_screen(positions, threshold=COARSE_THRESHOLD_KM)

    # Group by pair and filter out persistent co-orbital pairs
    candidate_pairs = defaultdict(list)
    for (i, j), t_idx in zip(raw_pairs, raw_t_indices):
        candidate_pairs[(min(i, j), max(i, j))].append(t_idx)
    candidate_pairs = {
        k: v for k, v in candidate_pairs.items() if len(v) <= MAX_COARSE_HITS
    }
    print(
        f"Phase 2: {len(candidate_pairs)} candidate pairs in {time.perf_counter() - t0:.3f}s"
    )

    if not candidate_pairs:
        print("No conjunction candidates found.")
        return

    # Phase 3: Get minimum distances for candidate pairs
    t0 = time.perf_counter()
    pair_arr = np.array(list(candidate_pairs.keys()), dtype=np.uint32)
    min_dists, min_times = min_distances(positions, pair_arr)
    print(
        f"Phase 3: Min distances for {len(pair_arr)} pairs in {time.perf_counter() - t0:.3f}s"
    )

    # Report close approaches
    events = [
        (pair_arr[i, 0], pair_arr[i, 1], min_dists[i], min_times[i])
        for i in np.where(min_dists < FINE_THRESHOLD_KM)[0]
    ]
    events.sort(key=lambda e: e[2])

    print(f"\nFound {len(events)} conjunctions (< {FINE_THRESHOLD_KM} km)")
    for i, (sat_i, sat_j, dist, t_min) in enumerate(events[:20]):
        print(
            f"  #{i + 1}: sat {sat_i:>5d} & {sat_j:>5d}  {dist:.3f} km  @ t={t_min} min"
        )


if __name__ == "__main__":
    main()
