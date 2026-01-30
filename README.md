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

Sub-meter accuracy validated against reference implementations. Uses SIMD (AVX512/AVX2) to process 8 satellites simultaneously, with multithreaded constellation propagation across all available cores.

#### Single-Threaded (1.2M propagations, single satellite)

| Implementation | Props/sec | Speedup vs python-sgp4 |
|----------------|-----------|------------------------|
| **astroz** | **30.8M** | **11x** |
| Rust sgp4 | 5.1M | 1.8x |
| heyoka | 3.8M | 1.3x |
| satkit | 3.5M | 1.2x |
| python-sgp4 | 2.8M | 1x |

#### Multi-Threaded Constellation (13,478 satellites × 1,440 times)

| Implementation | 1 Thread | 16 Threads |
|----------------|----------|------------|
| **astroz** | 37.7M/s | **303M/s** |
| heyoka | 15.7M/s | 155.6M/s |
| Rust sgp4 (rayon) | 4.4M/s | 47.9M/s |
| satkit | 3.5M/s | 3.5M/s |
| python-sgp4 | 2.7M/s | 2.7M/s |

*Benchmarked on AMD Ryzen 7 7840U (16 threads). All implementations using their optimal configurations (SIMD, pre-allocated outputs, batch mode).*

Uses SIMD (AVX512 for 8-wide, AVX2/SSE for 4-wide) with multithreaded time-major iteration. Validated against Vallado AIAA 2006-6753 reference vectors (< 10m position error, < 1µm/s velocity error). Set `ASTROZ_THREADS` environment variable to control thread count (defaults to all available cores).

The [Cesium visualization example](examples/README.md) propagates the entire active satellite catalog (~13,000 satellites) at interactive rates. **[Try the live demo →](https://attron.github.io/astroz-demo/)**

### Python

```bash
pip install astroz
```

```python
from astroz import propagate, Constellation
import numpy as np

# Load and propagate - automatically optimized for maximum performance
positions = propagate("starlink", np.arange(1440))  # 1 day at 1-min intervals
# shape: (1440, num_satellites, 3) in km, ECEF coordinates

# With options
from datetime import datetime, timezone
positions = propagate(
    "starlink",
    np.arange(1440),
    start_time=datetime(2024, 6, 15, tzinfo=timezone.utc),
    output="geodetic",  # "ecef" (default), "teme", or "geodetic"
)

# With velocities
positions, velocities = propagate("starlink", np.arange(1440), velocities=True)

# For repeated propagation, pre-parse to avoid overhead
c = Constellation("starlink")
positions = propagate(c, np.arange(1440))
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

  Interactive 3D visualization of the entire near-earth satellite catalog (~13,000 satellites) using Cesium. Features multithreaded SGP4 propagation at ~300M props/sec, constellation filtering, search, and satellite tracking.

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
