# astroz Python bindings

Python bindings for [astroz](https://github.com/ATTron/astroz), a high-performance astrodynamics library written in Zig.

## Performance

Using the batch numpy API, astroz significantly outperforms other Python SGP4 libraries:

| Scenario | astroz | python-sgp4 | skyfield | vs sgp4 |
|----------|--------|-------------|----------|---------|
| 1 day (minute res) | 0.37 ms | 0.54 ms | 8 ms | **1.5x faster** |
| 1 week (minute res) | 2.6 ms | 3.7 ms | 50 ms | **1.4x faster** |
| 2 weeks (minute res) | 4.8 ms | 7.4 ms | 92 ms | **1.5x faster** |
| 2 weeks (second res) | 280 ms | 460 ms | - | **1.6x faster** |
| 1 month (minute res) | 9.8 ms | 17 ms | 187 ms | **1.7x faster** |

The raw Zig implementation achieves ~5 million propagations/second, outperforming
the Rust [sgp4 crate](https://github.com/neuromorphicsystems/sgp4) by 1.3-2.2x.

## Installation

```bash
# Build the shared library first
cd /path/to/astroz
zig build c-api

# Install Python package
cd python
pip install -e .
```

## Usage

```python
import astroz
from astroz import Tle, Sgp4, orbital, coords

# Parse TLE and propagate with SGP4
tle = Tle("""1 25544U 98067A   24001.50000000  .00016717  00000-0  10270-3 0  9025
2 25544  51.6400 208.9163 0006703  40.6542 114.5435 15.49910020999999""")

sgp4 = Sgp4(tle)
pos, vel = sgp4.propagate(30.0)  # 30 minutes after epoch
print(f"Position: {pos} km")
print(f"Velocity: {vel} km/s")

# Convert to ground coordinates
gmst = coords.julian_to_gmst(2460310.5)
lat, lon, alt = coords.eci_to_geodetic(pos, gmst)
print(f"Lat: {lat:.2f}, Lon: {lon:.2f}, Alt: {alt:.1f} km")

# Orbital mechanics
transfer = orbital.hohmann(398600.5, 6778, 42164)  # LEO to GEO
print(f"Delta-V: {transfer['total_delta_v']:.2f} km/s")
```

## API

### TLE Parsing
- `Tle(tle_string)` - Parse two-line element set
- Properties: `satellite_number`, `epoch`, `inclination`, `eccentricity`, `mean_motion`

### SGP4 Propagation
- `Sgp4(tle, gravity_model=WGS84)` - Initialize propagator
- `propagate(tsince)` - Single propagation (returns tuple)
- `propagate_batch(times)` - Batch propagation (returns list)
- `propagate_batch_np(times)` - Batch with numpy (fastest, returns arrays)

### Coordinates
- `coords.eci_to_ecef(eci, gmst)` - ECI to ECEF
- `coords.ecef_to_geodetic(ecef)` - ECEF to lat/lon/alt
- `coords.eci_to_geodetic(eci, gmst)` - ECI to lat/lon/alt
- `coords.julian_to_gmst(jd)` - Julian date to GMST

### Orbital Mechanics
- `orbital.hohmann(mu, r1, r2)` - Hohmann transfer calculation
- `orbital.velocity(mu, radius, sma=0)` - Orbital velocity
- `orbital.period(mu, sma)` - Orbital period
- `orbital.escape_velocity(mu, radius)` - Escape velocity
