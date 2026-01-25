# astroz Python Bindings

High-performance SGP4 satellite propagation, powered by Zig with SIMD acceleration.

**Platforms:** macOS, Linux | **Requires:** Python 3.10+

## Quick Start

```python
from astroz import load_constellation, propagate_constellation
import numpy as np

# Load Starlink satellites and propagate for 1 day
constellation = load_constellation("starlink")
positions = propagate_constellation(constellation, np.arange(1440))
# positions: (1440, num_satellites, 3) in km, ECEF coordinates
```

## Loading TLEs

```python
from astroz import load_constellation

# CelesTrak groups
constellation = load_constellation("starlink")
constellation = load_constellation("iss")
constellation = load_constellation("gps")
constellation = load_constellation("all")  # ~10k active satellites

# By NORAD ID
constellation = load_constellation(norad_id=25544)  # ISS
constellation = load_constellation(norad_id=[25544, 48274])  # Multiple

# Local file or URL
constellation = load_constellation("satellites.tle")
constellation = load_constellation("https://example.com/tles.txt")

# With metadata (name, norad_id, inclination, period, etc.)
constellation, metadata = load_constellation("starlink", with_metadata=True)
for sat in metadata:
    print(f"{sat['name']}: {sat['inclination']:.1f}° inc, {sat['period']:.1f} min period")
```

Groups: `all`, `starlink`, `oneweb`, `planet`, `spire`, `gps`, `glonass`, `galileo`, `beidou`, `stations`/`iss`, `weather`, `geo`

## Propagation

```python
from astroz import load_constellation, propagate_constellation
from datetime import datetime, timezone
import numpy as np

constellation = load_constellation("starlink")

# Simple (defaults: now, ECEF)
positions = propagate_constellation(constellation, np.arange(1440))

# With options
positions = propagate_constellation(
    constellation,
    np.arange(14 * 1440),  # 2 weeks
    start_time=datetime(2024, 6, 1, tzinfo=timezone.utc),
    output="geodetic",  # "ecef" (default), "teme", or "geodetic"
)

# With velocities
positions, velocities = propagate_constellation(
    constellation, np.arange(1440), velocities=True
)
```

## Single Satellite

```python
from astroz import Tle, Sgp4
import numpy as np

tle = Tle("""1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995
2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123""")

sgp4 = Sgp4(tle)
pos, vel = sgp4.propagate(30.0)  # 30 min after epoch
positions, velocities = sgp4.propagate_batch(np.arange(1440))
```

## Collision Screening

```python
from astroz import load_constellation, propagate_constellation, coarse_screen, min_distances
import numpy as np

constellation = load_constellation("starlink")
positions = propagate_constellation(
    constellation, np.arange(1440),
    output="teme", layout="satellite_major"
)

# Find pairs within 10km
pairs, t_indices = coarse_screen(positions, threshold=10.0)

# Get exact minimum distances
pairs_array = np.array(pairs, dtype=np.uint32)
min_dists, min_times = min_distances(positions, pairs_array)
```

## Performance

| Constellation (13,448 sats × 1,440 steps) | Throughput |
|--------------------------------------------|------------|
| 1 thread | 7.7M props/sec |
| 16 threads | 56M props/sec |

Set `ASTROZ_THREADS` to control thread count.

## Building

Requires [Zig](https://ziglang.org/) and Python 3.10+.

```bash
cd bindings/python
pip install -e .
```
