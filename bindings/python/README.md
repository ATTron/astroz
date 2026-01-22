# astroz Python Bindings

High-performance SGP4 satellite orbit propagation for Python, powered by Zig with SIMD acceleration.

**Supported platforms:** macOS and Linux
**Requires:** Python 3.10+

## Quick Start

```python
from astroz import Tle, Sgp4
import numpy as np

# Parse TLE
tle = Tle("""1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995
2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123""")

# Single propagation
sgp4 = Sgp4(tle)
pos, vel = sgp4.propagate(30.0)  # 30 minutes after TLE epoch

# Batch propagation (convenience method)
times = np.arange(0, 1440, 1.0, dtype=np.float64)  # 1 day, 1-min intervals
positions, velocities = sgp4.propagate_batch(times)

# Or use propagate_into for zero-copy into pre-allocated arrays
positions = np.empty((len(times), 3), dtype=np.float64)
velocities = np.empty((len(times), 3), dtype=np.float64)
sgp4.propagate_into(times, positions, velocities)
```

## Performance

**1.3-2.9x faster** than python-sgp4:

| Scenario | astroz | python-sgp4 | Speedup |
|----------|--------|-------------|---------|
| 2 weeks (second intervals) | 160 ms | 464 ms | **2.9x** |
| 1 month (minute intervals) | 5.9 ms | 16.1 ms | **2.7x** |

## API

### Tle

```python
tle = Tle(tle_string)
tle.satellite_number  # NORAD catalog number
tle.epoch            # Epoch (J2000 seconds)
tle.inclination      # Degrees
tle.eccentricity
tle.mean_motion      # Rev/day
```

### Sgp4

```python
sgp4 = Sgp4(tle, gravity_model=WGS84)  # WGS84 (default) or WGS72

# Single point
pos, vel = sgp4.propagate(tsince)  # tsince in minutes
# Returns ((x,y,z), (vx,vy,vz)) in km and km/s (TEME frame)

# Batch (convenience) - returns allocated arrays
positions, velocities = sgp4.propagate_batch(times)

# Batch (zero-copy) - writes directly to pre-allocated arrays
sgp4.propagate_into(times, positions, velocities)
```

## Building

Requires [Zig](https://ziglang.org/) and Python 3.12+.

```bash
cd bindings/python
pip install -e .
```
