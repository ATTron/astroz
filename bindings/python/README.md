# astroz Python Bindings

High-performance SGP4 satellite propagation, powered by Zig with SIMD acceleration (AVX512/AVX2).

**Platforms:** macOS, Linux | **Requires:** Python 3.10+

## Quick Start

### python-sgp4 Compatible API (Recommended)

Drop-in replacement for [python-sgp4](https://github.com/brandon-rhodes/python-sgp4):

```python
from astroz.api import Satrec, SatrecArray, jday, WGS72
import numpy as np

# Single satellite (same syntax as python-sgp4)
sat = Satrec.twoline2rv(line1, line2, WGS72)
jd, fr = jday(2024, 5, 6, 12, 0, 0.0)
error, position, velocity = sat.sgp4(jd, fr)

# Batch propagation (270-330M props/sec with SIMD)
sat_array = SatrecArray(satrecs)
e, r, v = sat_array.sgp4(jd, fr)  # Scalars or arrays

# Skip velocities for 30% faster propagation
e, r, _ = sat_array.sgp4(jd_array, fr_array, velocities=False)
```

### High-Level API

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

## Migrating from python-sgp4

### Step 1: Change the Import (2x faster)

```python
# Before
from sgp4.api import Satrec, SatrecArray, jday

# After
from astroz.api import Satrec, SatrecArray, jday
```

That's it. Your existing code works unchanged and runs **2x faster**.

### Step 2: Use Batch Methods (5-12x faster)

If you have loops, switch to batch methods:

```python
# Before: loop (1.3M props/sec)
results = []
for jd, fr in zip(jd_list, fr_list):
    e, r, v = sat.sgp4(jd, fr)
    results.append(r)

# After: batch (15M props/sec) - 12x faster
jd_array = np.array(jd_list)
fr_array = np.array(fr_list)
e, r, v = sat.sgp4_array(jd_array, fr_array)
```

### Step 3: Use SatrecArray for Multiple Satellites (100x faster)

```python
# Before: loop over satellites (1.3M props/sec)
for sat in satellites:
    e, r, v = sat.sgp4_array(jd, fr)

# After: batch all satellites (290M props/sec) - 100x faster
sat_array = SatrecArray(satellites)
e, r, v = sat_array.sgp4(jd, fr)
```

### Step 4: Skip Velocities If Not Needed (30% faster)

```python
# If you only need positions
e, r, _ = sat_array.sgp4(jd, fr, velocities=False)
```

### Performance Summary

| Pattern | python-sgp4 | astroz | Speedup |
|---------|-------------|--------|---------|
| `sat.sgp4()` loop | 1.3M/s | 2.5M/s | **2x** |
| `sat.sgp4_array()` | 2.7M/s | 15M/s | **5x** |
| `SatrecArray.sgp4()` | 3M/s | 290M/s | **100x** |

### Compatibility Notes

astroz supports the core python-sgp4 API for satellite propagation. Some rarely-used attributes are not implemented:

- **TLE metadata**: `classification`, `intldesg`, `elnum`, `revnum`, `ephtype`
- **Intermediate elements**: `Om`, `am`, `em`, `im`, `mm`, `nm`, `om` (osculating elements after propagation)
- **Element rates**: `argpdot`, `mdot`, `nodedot`
- **Gravity constants**: `j2`, `j3`, `j4`, `mu`, etc. (available via `astroz.constants`)

For typical satellite tracking and visualization, astroz is a full drop-in replacement.

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
