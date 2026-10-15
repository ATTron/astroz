"""Screen a constellation for close approaches over one day."""

import sys
import time

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
    threshold_km = 5.0

    print(f"Loading {source}...")
    t0 = time.perf_counter()
    c = Constellation(norad_id=int(source)) if source.isdigit() else Constellation(source)
    print(f"Loaded {c.num_satellites} satellites in {time.perf_counter() - t0:.3f}s")

    print(f"Screening one day for approaches within {threshold_km} km...")
    t0 = time.perf_counter()
    events = screen(c, 1440, threshold=threshold_km)
    print(f"Found {len(events.i)} encounters in {time.perf_counter() - t0:.3f}s")

    closest = events.distance.argsort()[:20]
    for k in closest:
        print(
            f"  sat {events.i[k]:>5d} & {events.j[k]:>5d}: {events.distance[k]:6.3f} km "
            f"at t+{events.time[k]:7.2f} min, {events.speed[k]:5.2f} km/s"
        )

    # Screening one satellite against the rest is much cheaper than all pairs
    t0 = time.perf_counter()
    mine = screen(c, 1440, threshold=threshold_km, target=0)
    print(f"\n{len(mine.i)} encounters involving satellite 0 ({time.perf_counter() - t0:.3f}s)")


if __name__ == "__main__":
    main()
