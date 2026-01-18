<h1 align="center">
  <img src="https://repository-images.githubusercontent.com/819657891/7fdb22c8-7171-4b75-9f33-88a62ea67259" width="900" height="600" />
</h1>

[![CI][ci-shd]][ci-url]
[![CD][cd-shd]][cd-url]
[![DC][dc-shd]][dc-url]

## Astronomical and Spacecraft Toolkit Written in Zig

**Featuring the fastest open-source SGP4 propagator.**

| Orbital Mechanics | Spacecraft Ops | Astronomy |
|-------------------|----------------|-----------|
| SGP4 propagation | CCSDS packets | FITS parsing |
| TLE parsing | VITA49 packets | WCS coordinates |
| Orbital maneuvers | Attitude determination | Star precession |
| Monte Carlo sims | | Celestial bodies |

### Performance

Sub-meter accuracy validated against reference implementations. Uses SIMD (AVX2/SSE) to process 4 orbits simultaneously.

#### Native (2 weeks @ second resolution)

| Implementation | Props/sec |
|----------------|-----------|
| **astroz (SIMD)** | 10.8M |
| **astroz (scalar)** | 5.3M |
| sgp4 (Rust) | 4.8M |

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
