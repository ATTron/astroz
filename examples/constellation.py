"""Load a constellation, propagate it to Earth-fixed and geodetic coordinates, and screen it for close approaches.

Usage: python constellation.py <source>
  source: a CelesTrak group (starlink, stations, active), a NORAD ID, or a TLE / OMM JSON file
"""

import sys
import time

from astroz import Constellation, propagate, screen


def main():
    if len(sys.argv) < 2:
        print(__doc__.split("\n\n")[1])
        return

    source = sys.argv[1]
    c = Constellation(norad_id=int(source)) if source.isdigit() else Constellation(source)
    print(f"{c.num_satellites} satellites")

    # Times are minutes from start_time (default: now); columns follow c.satnums
    ids = c.satnums
    ecef = propagate(c, [0], output="ecef")[0]
    geo = propagate(c, [0], output="geodetic")[0]
    for s in range(min(3, c.num_satellites)):
        (x, y, z), (lat, lon, alt) = ecef[s], geo[s]
        print(f"  {ids[s]:>6}: ({x:9.1f}, {y:9.1f}, {z:9.1f}) km ECEF = {lat:6.2f} deg, {lon:7.2f} deg, {alt:6.1f} km")

    t0 = time.perf_counter()
    events = screen(c, 1440, threshold=5.0)
    print(f"\n{len(events.i)} approaches within 5 km over one day ({time.perf_counter() - t0:.2f}s); closest:")
    for k in events.distance.argsort()[:10]:
        print(f"  {ids[events.i[k]]:>6} & {ids[events.j[k]]:>6}: {events.distance[k]:6.3f} km "
              f"at t+{events.time[k]:7.2f} min, {events.speed[k]:5.2f} km/s")

    # Screening one satellite against the rest is much cheaper than all pairs
    mine = screen(c, 1440, threshold=5.0, target=0)
    print(f"\n{len(mine.i)} of them involve {ids[0]}")


if __name__ == "__main__":
    main()
