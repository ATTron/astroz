# astroz Python Bindings

High-performance SGP4 satellite propagation, powered by Zig with SIMD acceleration (AVX512/AVX2).

**Platforms:** macOS, Linux | **Requires:** Python 3.10+

## Quick Start

```python
from astroz import Constellation
import numpy as np

# Load and propagate - automatically optimized for 300M+ props/sec
constellation = Constellation("starlink")
positions = constellation.propagate(np.arange(1440))  # 1 day at 1-min intervals
# positions: (1440, num_satellites, 3) in km, ECEF coordinates
```

## Loading Constellations

```python
from astroz import Constellation

# CelesTrak groups
constellation = Constellation("starlink")
constellation = Constellation("iss")
constellation = Constellation("gps")
constellation = Constellation("all")  # ~13k active satellites

# By NORAD ID
constellation = Constellation(norad_id=25544)  # ISS
constellation = Constellation(norad_id=[25544, 48274])  # Multiple

# Local file or URL
constellation = Constellation("satellites.tle")
constellation = Constellation("https://example.com/tles.txt")

# With metadata (name, norad_id, inclination, period, etc.)
constellation = Constellation("starlink", with_metadata=True)
for sat in constellation.metadata:
    print(f"{sat['name']}: {sat['inclination']:.1f}° inc, {sat['period']:.1f} min period")
```

Groups: `all`, `starlink`, `oneweb`, `planet`, `spire`, `gps`, `glonass`, `galileo`, `beidou`, `stations`/`iss`, `weather`, `geo`

## Propagation

```python
from astroz import Constellation
from datetime import datetime, timezone
import numpy as np

constellation = Constellation("starlink")

# Simple (defaults: now UTC, ECEF output)
positions = constellation.propagate(np.arange(1440))

# With options
positions = constellation.propagate(
    np.arange(14 * 1440),  # 2 weeks
    start_time=datetime(2024, 6, 1, tzinfo=timezone.utc),
    output="geodetic",  # "ecef" (default), "teme", or "geodetic"
)

# With velocities
positions, velocities = constellation.propagate(np.arange(1440), velocities=True)
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
from astroz import Constellation, coarse_screen, min_distances
import numpy as np

constellation = Constellation("starlink")
positions = constellation.propagate(np.arange(1440), output="teme")

# Find pairs within 10km
pairs, t_indices = coarse_screen(positions, threshold=10.0)

# Get exact minimum distances
pairs_array = np.array(pairs, dtype=np.uint32)
min_dists, min_times = min_distances(positions, pairs_array)
```

## Performance

| Constellation (13,478 sats × 1,440 steps) | Throughput |
|--------------------------------------------|------------|
| 1 thread | 37.7M props/sec |
| 16 threads | 303M props/sec |

*Benchmarked on AMD Ryzen 7 7840U with AVX512. The `Constellation` class automatically handles buffer reuse and warmup for optimal performance.*

Set `ASTROZ_THREADS` to control thread count (defaults to all cores).

## Advanced: Pre-allocated Buffers

For maximum performance with repeated propagations, pre-allocate output buffers:

```python
from astroz import Constellation
import numpy as np

constellation = Constellation("starlink")
times = np.arange(1440, dtype=np.float64)

# Pre-allocate buffer
out = np.empty((len(times), constellation.num_satellites, 3), dtype=np.float64)
positions = constellation.propagate(times, out=out)
```

Note: The `Constellation` class automatically reuses buffers internally when the output shape matches previous calls.

## Building

Requires [Zig](https://ziglang.org/) and Python 3.10+.

```bash
cd bindings/python
pip install -e .
```
