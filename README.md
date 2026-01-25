<h1 align="center">
  <img src="https://repository-images.githubusercontent.com/819657891/7fdb22c8-7171-4b75-9f33-88a62ea67259" width="900" height="600" />
</h1>

[![CI][ci-shd]][ci-url]
[![CD][cd-shd]][cd-url]
[![DC][dc-shd]][dc-url]

## Astronomical and Spacecraft Toolkit Written in Zig

**Featuring the fastest CPU based open source SGP4 propagator**

| Orbital Mechanics | Spacecraft Ops | Astronomy |
|-------------------|----------------|-----------|
| SGP4 propagation | CCSDS packets | FITS parsing |
| TLE parsing | VITA49 packets | WCS coordinates |
| Orbital maneuvers | Attitude determination | Star precession |
| Monte Carlo sims | | Celestial bodies |

### Performance

Sub-meter accuracy validated against reference implementations. Uses SIMD (AVX2/SSE) to process 4 satellites simultaneously, with multithreaded constellation propagation across all available cores.

#### Single Satellite (Python)

| Scenario | astroz | python-sgp4 | Speedup |
|----------|--------|-------------|---------|
| 1 day (minute res) | 0.24 ms | 0.65 ms | **2.7x** |
| 1 week (minute res) | 1.99 ms | 3.41 ms | **1.7x** |
| 2 weeks (minute res) | 3.22 ms | 7.00 ms | **2.2x** |
| 2 weeks (second res) | 144 ms | 438 ms | **3.0x** |
| 1 month (minute res) | 5.29 ms | 14.84 ms | **2.8x** |

#### Multi-Satellite Constellation

| Mode | Throughput |
|------|------------|
| Single-core (SIMD) | **22M props/sec** |
| Multithreaded | **200M+ props/sec** |

Uses SIMD (AVX2/SSE) to process 4 satellites per batch with optional multithreaded time-major iteration. Validated against Vallado AIAA 2006-6753 reference vectors (< 10m position error, < 1µm/s velocity error). Set `ASTROZ_THREADS` environment variable to control thread count (defaults to all available cores).

The [Cesium visualization example](examples/README.md) propagates the entire active satellite catalog (~13,000 satellites) at interactive rates. **[Try the live demo →](https://attron.github.io/astroz-demo/)**

### Python

```bash
pip install astroz
```

```python
from astroz import Tle, Sgp4
import numpy as np

tle = Tle("1 25544U 98067A   24127.82853009 ...\n2 25544  51.6393 ...")
sgp4 = Sgp4(tle)

# Single propagation
pos, vel = sgp4.propagate(30.0)  # 30 min after epoch

# Batch propagation (convenience method)
times = np.arange(1209600, dtype=np.float64) / 60.0  # 2 weeks in minutes
positions, velocities = sgp4.propagate_batch(times)

# Or use propagate_into for zero-copy into pre-allocated arrays
positions = np.empty((len(times), 3), dtype=np.float64)
velocities = np.empty((len(times), 3), dtype=np.float64)
sgp4.propagate_into(times, positions, velocities)

# Multi-satellite constellation (SIMD + multithreaded)
from astroz import load_constellation, propagate_constellation

# Load from CelesTrak group, file, URL, or TLE string
constellation = load_constellation("starlink")  # or load_constellation(norad_id=25544)

# Propagate for 1 day at 1-minute intervals (defaults to current UTC time)
times = np.arange(1440, dtype=np.float64)
positions = propagate_constellation(constellation, times, output="ecef")
# shape: (num_times, num_sats, 3)

# With velocities and custom start time
from datetime import datetime, timezone
positions, velocities = propagate_constellation(
    constellation, times,
    start_time=datetime(2024, 6, 15, tzinfo=timezone.utc),
    output="ecef",
    velocities=True,
)
```

### Usage

- Add `astroz` as a dependency in your `build.zig.zon`.

```sh
zig fetch --save https://github.com/ATTron/astroz/archive/<git_tag_or_commit_hash>.tar.gz
#or
zig fetch --save git+https://github.com/ATTron/astroz/#HEAD
```

- Use `astroz` as a module in your `build.zig`.

```zig
const astroz_dep = b.dependency("astroz", .{
    .target = target,
    .optimize = optimize,
});
const astroz_mod = astroz_dep.module("astroz");
exe.root_module.addImport("astroz", astroz_mod);
```

### Examples

#### Orbital Mechanics

- #### [Planet Transfer & Mission Planning](examples/transfer_propagation.zig)
  Demonstrates interplanetary transfers with mission planning (Hohmann vs Bi-Elliptic comparison) and trajectory propagation.

- #### [Orbit Maneuvers](examples/orbit_maneuvers.zig)
  Comprehensive example showing TLE-based orbit propagation with various maneuver types: impulse, plane change, and phase change.

- #### [Monte Carlo Simulation](examples/simple_monte_carlo.zig)
  Statistical analysis for mission planning with uncertainty.

- #### [SGP4 Propagation](examples/sgp4_propagation.zig)

  Analytical orbit propagation using SGP4 with TLE input. Demonstrates both direct SGP4 usage and the modular propagator interface.

- #### [Cesium Satellite Visualization](examples/README.md) — **[Live Demo](https://attron.github.io/astroz-demo/)**

  Interactive 3D visualization of the entire near-earth satellite catalog (~13,000 satellites) using Cesium. Features multithreaded SGP4 propagation at ~190M props/sec, constellation filtering, search, and satellite tracking.

#### Spacecraft Operations

- #### [Orbit Orientation Determination](examples/simple_spacecraft_orientation.zig)
  Calculate spacecraft attitude and orientation.

#### Telemetry & Data Handling

- #### [Parse Vita49](examples/parse_vita49.zig) / [with Callback](examples/parse_vita49_callback.zig)
  VITA Radio Transport (VRT) packet stream parsing.

- #### [Parse CCSDS](examples/parse_ccsds.zig) / [with File Sync](examples/parse_ccsds_file_sync.zig)
  Parse CCSDS space packet protocol from files.

- #### [Create CCSDS Packet](examples/create_ccsds_packet.zig) / [with Config](examples/create_ccsds_packet_config.zig)
  Generate CCSDS packets for telemetry.

#### Astronomy & Astrometry

- #### [Generate Image from FITS File](examples/parse_fits_file.zig)
  Parse and render FITS astronomical image data.

<img src="test/test.png" width="450" height="400" alt="sample fits image as png"/>

- #### [Precess Star Coordinates](examples/precess_star.zig)
  Calculate stellar precession to a target epoch.

- #### [Calculate WCS from TLE](examples/wcs.zig)
  Compute World Coordinate System values from orbital elements.

<!-- MARKDOWN LINKS -->

[ci-shd]: https://img.shields.io/github/actions/workflow/status/ATTron/astroz/ci.yaml?branch=main&logo=github&label=CI&labelColor=black
[ci-url]: https://github.com/ATTron/astroz/blob/main/.github/workflows/ci.yaml
[cd-shd]: https://img.shields.io/github/actions/workflow/status/ATTron/astroz/cd.yaml?branch=main&logo=github&label=CD&labelColor=black
[cd-url]: https://github.com/ATTron/astroz/blob/main/.github/workflows/cd.yaml
[dc-shd]: https://img.shields.io/badge/click-F6A516?logo=zig&logoColor=F6A516&label=doc&labelColor=black
[dc-url]: https://attron.github.io/astroz
