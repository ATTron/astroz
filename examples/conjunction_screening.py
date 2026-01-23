"""Conjunction screening example using astroz high-performance APIs.

Demonstrates the full pipeline: TLE parsing, coarse propagation,
spatial hash screening, fine propagation with satellite masks,
and pairwise distance computation.
"""

import sys
import time
from collections import defaultdict
import numpy as np
from astroz import (
    Sgp4Constellation,
    propagate_constellation,
    min_distances,
    coarse_screen,
)

COARSE_THRESHOLD_KM = 10.0
FINE_THRESHOLD_KM = 5.0
FINE_WINDOW_MIN = 5.0
FINE_STEPS = 600
MAX_WINDOW_SPAN = 30.0
MAX_COARSE_HITS = 10


def merge_windows(t_indices):
    """Merge nearby coarse time indices into contiguous windows."""
    sorted_t = sorted(t_indices)
    groups = [[sorted_t[0], sorted_t[0], [sorted_t[0]]]]
    for t in sorted_t[1:]:
        g = groups[-1]
        if t <= g[1] + 2 * FINE_WINDOW_MIN and t - g[0] <= MAX_WINDOW_SPAN:
            g[1] = t
            g[2].append(t)
        else:
            groups.append([t, t, [t]])
    return groups


def main():
    if len(sys.argv) < 2:
        print("Usage: python conjunction_screening.py <tle_file> [reference_jd]")
        return

    with open(sys.argv[1]) as f:
        tle_data = f.read()

    # Phase 0: Parse TLEs
    t0 = time.perf_counter()
    constellation = Sgp4Constellation.from_tle_text(tle_data)
    num_sats = constellation.num_satellites
    print(f"Parsed {num_sats} satellites in {time.perf_counter() - t0:.3f}s")

    epochs = np.array(constellation.epochs, dtype=np.float64)
    ref_jd = float(sys.argv[2]) if len(sys.argv) >= 3 else np.median(epochs)
    epoch_offsets = (ref_jd - epochs) * 1440.0

    # Phase 1: Coarse propagation
    t0 = time.perf_counter()
    base_times = np.arange(1440, dtype=np.float64)
    positions = propagate_constellation(
        constellation,
        base_times,
        epoch_offsets=epoch_offsets,
        output="teme",
        reference_jd=ref_jd,
    )
    print(
        f"Phase 1: Coarse propagation ({num_sats} x 1440) in {time.perf_counter() - t0:.3f}s"
    )

    # Phase 1b: Spatial hash screening
    t0 = time.perf_counter()
    valid_mask = (np.linalg.norm(positions[:, 0, :], axis=1) > 100.0).astype(np.uint8)
    raw_pairs, raw_t_indices = coarse_screen(
        positions,
        threshold=COARSE_THRESHOLD_KM,
        valid_mask=valid_mask,
    )

    candidate_pairs = defaultdict(list)
    for (i, j), t_idx in zip(raw_pairs, raw_t_indices):
        candidate_pairs[(min(i, j), max(i, j))].append(t_idx)

    # Filter persistent co-orbital pairs
    candidate_pairs = {
        k: v for k, v in candidate_pairs.items() if len(v) <= MAX_COARSE_HITS
    }
    print(
        f"Phase 1b: {len(candidate_pairs)} candidate pairs in {time.perf_counter() - t0:.3f}s"
    )

    if not candidate_pairs:
        print("No conjunction candidates found.")
        return

    # Phase 2: Fine propagation per window
    t0 = time.perf_counter()

    pairs_by_t = defaultdict(list)
    for pair, t_indices in candidate_pairs.items():
        for t in t_indices:
            pairs_by_t[t].append(pair)

    merged_groups = merge_windows(pairs_by_t.keys())
    events = []
    mask = np.zeros(num_sats, dtype=np.uint8)

    for g_start, g_end, g_t_indices in merged_groups:
        group_pairs = {p for t in g_t_indices for p in pairs_by_t[t]}
        sat_set = {s for pair in group_pairs for s in pair}

        mask[:] = 0
        for s in sat_set:
            mask[s] = 1

        t_start, t_end = (
            float(g_start) - FINE_WINDOW_MIN,
            float(g_end) + FINE_WINDOW_MIN,
        )
        fine_times = np.linspace(
            t_start,
            t_end,
            max(FINE_STEPS, int((t_end - t_start) * 10)),
            dtype=np.float64,
        )

        fine_pos, fine_vel = propagate_constellation(
            constellation,
            fine_times,
            epoch_offsets=epoch_offsets,
            satellite_mask=mask,
            output="teme",
            reference_jd=ref_jd,
            velocities=True,
        )

        pair_arr = np.array(list(group_pairs), dtype=np.uint32)
        dists, indices, vels = min_distances(fine_pos, pair_arr, fine_vel)

        for k in np.where(dists < FINE_THRESHOLD_KM)[0]:
            events.append(
                {
                    "sat_i": int(pair_arr[k, 0]),
                    "sat_j": int(pair_arr[k, 1]),
                    "dist_km": dists[k],
                    "t_min": fine_times[int(indices[k])],
                    "vrel_km_s": vels[k],
                }
            )

    t_fine = time.perf_counter() - t0
    print(f"Phase 2: Fine screening ({len(merged_groups)} windows) in {t_fine:.3f}s")

    print(f"\nFound {len(events)} conjunction events (< {FINE_THRESHOLD_KM} km)")
    events.sort(key=lambda e: e["dist_km"])
    for i, e in enumerate(events[:20]):
        print(
            f"  #{i + 1}: {e['sat_i']:>5d} & {e['sat_j']:>5d}  "
            f"{e['dist_km']:.3f} km  t={e['t_min']:.2f} min  Vrel={e['vrel_km_s']:.3f} km/s"
        )


if __name__ == "__main__":
    main()
