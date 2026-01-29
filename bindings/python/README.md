# astroz Python Bindings

High-performance SGP4 satellite propagation, powered by Zig with SIMD acceleration (AVX512/AVX2).

**Platforms:** macOS, Linux | **Requires:** Python 3.10+

## Quick Start

```python
from astroz import propagate
import numpy as np

# Load and propagate - automatically optimized for 300M+ props/sec
positions = propagate("starlink", np.arange(1440))  # 1 day at 1-min intervals
# positions: (1440, num_satellites, 3) in km, ECEF coordinates
```

## Loading Sources

```python
from astroz import propagate, Constellation

# CelesTrak groups
positions = propagate("starlink", times)
positions = propagate("iss", times)
positions = propagate("gps", times)
positions = propagate("all", times)  # ~13k active satellites

# By NORAD ID
positions = propagate(None, times, norad_id=25544)  # ISS
positions = propagate(None, times, norad_id=[25544, 48274])  # Multiple

# Local file or URL
positions = propagate("satellites.tle", times)
positions = propagate("https://example.com/tles.txt", times)

# For repeated propagation, pre-parse to avoid overhead
c = Constellation("starlink")
positions1 = propagate(c, times1)
positions2 = propagate(c, times2)
```

Groups: `all`, `starlink`, `oneweb`, `planet`, `spire`, `gps`, `glonass`, `galileo`, `beidou`, `stations`/`iss`, `weather`, `geo`

## Propagation

```python
from astroz import propagate, Constellation
from datetime import datetime, timezone
import numpy as np

times = np.arange(1440)  # 1 day at 1-min intervals

# Simple (defaults: now UTC, ECEF output)
positions = propagate("starlink", times)

# With options
positions = propagate(
    "starlink",
    np.arange(14 * 1440),  # 2 weeks
    start_time=datetime(2024, 6, 1, tzinfo=timezone.utc),
    output="geodetic",  # "ecef" (default), "teme", or "geodetic"
)

# With velocities
positions, velocities = propagate("starlink", times, velocities=True)
```

## Conjunction Screening

```python
from astroz import screen
import numpy as np

times = np.arange(1440)

# Single target (fastest - uses fused propagate+screen, no full position array)
min_dists, min_t_indices = screen("starlink", times, threshold=50.0, target=0)
# Returns per-satellite minimum distance to target and time index

# All-vs-all screening
pairs, t_indices = screen("starlink", times, threshold=10.0)
# Returns all conjunction events within threshold
```

## Performance

| Constellation (13,478 sats x 1,440 steps) | Throughput |
|--------------------------------------------|------------|
| 1 thread | 37.7M props/sec |
| 16 threads | 303M props/sec |

*Benchmarked on AMD Ryzen 7 7840U with AVX512.*

Set `ASTROZ_THREADS` to control thread count (defaults to all cores).

## Building

Requires [Zig](https://ziglang.org/) and Python 3.10+.

```bash
cd bindings/python
pip install -e .
```
