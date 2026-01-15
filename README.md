<h1 align="center">
  <img src="https://repository-images.githubusercontent.com/819657891/7fdb22c8-7171-4b75-9f33-88a62ea67259" width="900" height="600" />
</h1>

[![CI][ci-shd]][ci-url]
[![CD][cd-shd]][cd-url]
[![DC][dc-shd]][dc-url]


## Astronomical and Spacecraft Toolkit Written in Zig

### Features / Plans

#### Spacecraft

- [x] CCSDS Packets
  - [x] CCSDS Stream Parser
- [x] VITA49 Packets
  - [x] Vita49 Stream Parser
- [x] TLE Support
  - [x] Orbital Propagation
    - [x] RK4 Numerical Integration
    - [x] SGP4 Analytical Propagation
  - [x] Modular Propagator System
    - [x] Composable Force Models (TwoBody, J2, Drag)
    - [x] Swappable Integrators
  - [x] Orbital Maneuvers
    - [x] Impulse Maneuvers
    - [x] Phase Maneuvers
    - [x] Plane Change Maneuvers
- [x] Orientation Determination
- [ ] SPICE Kernel Parsing
    - [x] LSK

#### Astronomical

- [x] Astronomical References
  - [x] J2000 and JD
  - [x] Celestial Bodies
    - [x] Mass
    - [x] Radius
    - [x] Orbital Details
- [x] Astronomical Coordinates
  - [x] Equatorial Coordinate System
  - [x] World Coordinate System
- [x] Astronomical Computation
  - [x] Precession
- [x] Celestial Bodies
- [x] Orbital Mechanics
  - [x] Interplanetary Maneuvers
  - [x] Monte Carlo Simulations
- [x] FITS File Parsing
  - [x] Image Generation
    - [x] MEF Parsing
  - [x] Table Parsing

### SGP4 Performance

Fastest open-source SGP4 implementation available, with sub-meter accuracy

| Implementation | Propagations/sec | vs python-sgp4 |
|----------------|------------------|----------------|
| astroz (Zig) | ~5M/s | 2x faster |
| astroz (Python) | ~4M/s | 1.6x faster |
| Rust sgp4 | ~4M/s | 1.6x faster |
| python-sgp4 | ~2.5M/s | baseline |

### Python Bindings

```bash
pip install astroz
```

```python
from astroz import Tle, Sgp4, orbital, coords

tle = Tle("1 25544U ...")
sgp4 = Sgp4(tle)
pos, vel = sgp4.propagate(30.0)  # 30 min after epoch
```

See [python/README.md](python/README.md) for full API documentation.

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

<img src="assets/earth_mars_transfer.gif" width="300" height="300" alt="Earth-Mars transfer animation"/> <img src="assets/earth_mars_transfer_static.png" width="300" height="300" alt="Earth-Mars transfer trajectory"/> <img src="assets/earth_mars_transfer_analysis.png" width="300" height="300" alt="Earth-Mars transfer analysis"/>

- #### [Orbit Maneuvers](examples/orbit_maneuvers.zig)
  Comprehensive example showing TLE-based orbit propagation with various maneuver types: impulse, plane change, and phase change.

<img src="assets/orbit_prop.gif" width="450" height="400" alt="visualization of orbit prop"/>

- #### [Monte Carlo Simulation](examples/simple_monte_carlo.zig)
  Statistical analysis for mission planning with uncertainty.

<img src="assets/monte_carlo_analysis.png" width="450" height="400" alt="graphs showing monte carlo run analysis"/>

- #### [SGP4 Propagation](examples/sgp4_propagation.zig)

  Analytical orbit propagation using SGP4 with TLE input. Demonstrates both direct SGP4 usage and the modular propagator interface.

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
