# ASTROZ  

![Testing](https://github.com/ATTron/astroz/actions/workflows/test_zig.yml/badge.svg) 

<img src="https://repository-images.githubusercontent.com/819657891/291c28ef-4c03-4d0e-bb0c-41d4662867c3" width="100" height="100"/> 

### Astronomical and Spacecraft Toolkit Written in Zig for Zig!  

## Features / Plans

### Spacecraft

- [x] CCSDS Packets
  - [x] CCSDS Stream Parser
- [x] VITA49 Packets
  - [x] Vita49 Stream Parser
- [x] TLE Support
  - [x] Orbital Propagation
    - [x] RK4
  - [ ] Orbital Maneuvers
    - [x] Impulse Maneuvers
    - [ ] Phase Maneuvers
    - [ ] Plane Change Maneuvers

### Astronomical

- [x] Astronomical References
  - [x] J2000 and JD
  - [x] Celestial Bodies
    - [x] Mass
    - [x] Radius
    - [x] Orbital Details
- [x] Astronomical Coordinates
  - [x] Equatorial Coordinate System
- [x] Astronomical Computation
  - [x] Precession
- [x] Celestial Bodies
- [ ] Orbital Mechanics
  - [ ] Interplanetary Maneuvers

### Feature not listed ?

To request a feature, please create an issue for this project and I will try my
best to be responsive.  

## Install

**Please use the master branch of the zig repository as that is what I'm developing against**

The easiest way I've found to get started with dependencies in zig is the following.  

- in your `main.zig` import the dependency.  

```zig
const astroz = @import("astroz");

```

- run `zig fetch --save git+https://github.com/ATTron/astroz/#HEAD`  
- inside `build.zig`  

```zig
const package = b.dependency("astroz", .{
    .target = target,
    .optimize = optimize,
});

const module = package.module("astroz");

exe.root_module.addImport("astroz", module);

b.installArtifact(exe);

```

## Usage

### Examples

#### Parse a TLE

```zig
const std = @import("std");
const astroz = @import("astroz");
const TLE = astroz.tle.TLE;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const test_tle =
        \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
        \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
    ;

    var tle = try TLE.parse(raw_tle, allocator);
    defer tle.deinit();

    tle.output();

}

```

#### Orbit Prop for the next 3 days

```zig
const std = @import("std");
const math = std.math;
const astroz = @import("astroz");
const TLE = astroz.tle.TLE;
const constants = astroz.constants;
const spacecraft = astroz.spacecraft;
const Spacecraft = spacecraft.Spacecraft;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const test_tle =
        \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
        \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
    ;

    var tle = try TLE.parse(test_tle, allocator);
    defer tle.deinit();

    var test_sc = Spacecraft.create("dummy_sc", tle, 300.000, spacecraft.Satellite_Size.Cube, constants.earth, allocator);
    defer test_sc.deinit();

    try test_sc.propagate(
        test_sc.tle.first_line.epoch,
        test_sc.tle.first_line.epoch + 3 * 86400.0, // 3 days worth of orbit predictions
        1,
        null,
    );

    for (test_sc.orbit_predictions.items) |iter| {
        const r = math.sqrt(iter.state[0] * iter.state[0] + iter.state[1] * iter.state[1] + iter.state[2] * iter.state[2]);

        std.debug.print("Next Prediction is: {any}\n", .{r});
    }
}
```

<img src="https://raw.githubusercontent.com/ATTron/astroz/main/assets/orbit_prop.gif" width="450" height="400" alt="visualization of orbit prop"/>

#### Orbit Prop for the next 3 days w/ impulse manuevers

```zig
const std = @import("std");
const math = std.math;
const astroz = @import("astroz");
const TLE = astroz.tle.TLE;
const constants = astroz.constants;
const spacecraft = astroz.spacecraft;
const Spacecraft = spacecraft.Spacecraft;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const test_tle =
        \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
        \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
    ;

    var tle = try TLE.parse(test_tle, allocator);
    defer tle.deinit();

    var test_sc = Spacecraft.create("dummy_sc", tle, 300.000, spacecraft.Satellite_Size.Cube, constants.earth, allocator);
    defer test_sc.deinit();

    const impulses = [_]Impulse{
        .{ .time = 3600.0, .delta_v = .{ 0.05, 0.03, 0.01 } },
        .{ .time = 7200.0, .delta_v = .{ 1.1, -0.05, 0.02 } },
        .{ .time = 10800.0, .delta_v = .{ -0.03, 0.08, -0.01 } },
    };

    try test_sc.propagate(
        test_sc.tle.first_line.epoch,
        test_sc.tle.first_line.epoch + 3 * 86400.0, // 3 days worth of orbit predictions
        1,
        &impulses,
    );

    for (test_sc.orbit_predictions.items) |iter| {
        const r = math.sqrt(iter.state[0] * iter.state[0] + iter.state[1] * iter.state[1] + iter.state[2] * iter.state[2]);

        std.debug.print("Next Prediction is: {any}\n", .{r});
    }
}
```

<img src="https://raw.githubusercontent.com/ATTron/astroz/main/assets/orbit_prop_w_impulse.gif" width="450" height="400" alt="visualization of orbit prop with impulses"/>

#### Setup Vita49 Parser

##### W/ Callback

```zig
const std = @import("std");
const astroz = @import("astroz");
const Vita49 = astroz.vita49.Vita49;
const Parser = astroz.parsers.Parser;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const P = Parser(Vita49);
    const ip = "127.0.0.1".*;
    const port: u16 = 65432;
    var parser = P.new(&ip, port, 1024, allocator);
    defer parser.deinit();
    _ = try parser.start(callback);
}

fn callback(packet: Vita49) void {
    std.debug.print("Packet recieved: {any}", .{packet});
}
```

##### W/O Callback

```zig
const std = @import("std");
const astroz = @import("astroz");
const Vita49 = astroz.vita49.Vita49;
const Parser = astroz.parsers.Parser;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const P = Parser(Vita49);
    const ip = "127.0.0.1".*;
    const port: u16 = 65432;
    var parser = P.new(&ip, port, 1024, allocator);
    defer parser.deinit();
    _ = try parser.start(null);
}

```

#### CCSDS Parser from File w/o file sync

```zig
const std = @import("std");
const astroz = @import("astroz");
const CCSDS = astroz.ccsds.CCSDS;
const Parser = astroz.parsers.Parser;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const file_name = "./test/files/ccsds.bin".*;

    const P = Parser(CCSDS);
    var parser = P.new(null, null, 1024, allocator);
    defer parser.deinit();

    _ = try parser.parse_from_file(&file_name, null, null);

    for (parser.packets.items) |packet| {
        std.log.info("Packets from files: 0x{x}", .{packet.packets});
    }
}
```

#### CCSDS Parser from File w/ file sync

```zig
const std = @import("std");
const astroz = @import("astroz");
const CCSDS = astroz.ccsds.CCSDS;
const Parser = astroz.parsers.Parser;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const file_name = "./test/files/ccsds.bin".*;
    const sync_pattern = .{ 0x78, 0x97, 0xC0, 0x00, 0x00, 0x0A, 0x01, 0x02 };

    const P = Parser(CCSDS);
    var parser = P.new(null, null, 1024, allocator);
    defer parser.deinit();

    _ = try parser.parse_from_file(&file_name, &sync_pattern, null);

    for (parser.packets.items) |packet| {
        std.log.info("Packets from files: 0x{x}", .{packet.packets});
    }
}
```

#### Create a CCSDS Packet

##### W/O Config

```zig
const std = @import("std");
const astroz = @import("astroz");
const CCSDS = astroz.ccsds.CCSDS;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const raw_test_packet: [16]u8 = .{ 0x78, 0x97, 0xC0, 0x00, 0x00, 0x0A, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A };
    const converted_test_packet = CCSDS.new(&raw_test_packet, allocator, null);
    defer converted_test_packet.deinit();

    std.debug.print("CCSDS Packet Created:\n{any}", .{converted_test_packet});
}
```

##### W/ Config

**config.json**

```json
{
  "secondary_header_length": 12
}
```

```zig
const std = @import("std");
const astroz = @import("astroz");
const ccsds = astroz.ccsds;
const CCSDS = ccsds.CCSDS;
const Config = ccsds.Config;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const config_file = try std.fs.cwd().readFileAlloc(allocator, "config.json", 512);
    defer allocator.free(config_file);

    const config = try ccsds.parse_config(config_file, allocator);

    const raw_test_packet: [16]u8 = .{ 0x78, 0x97, 0xC0, 0x00, 0x00, 0x0A, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A };
    const converted_test_packet = CCSDS.new(&raw_test_packet, allocator, config);
    defer converted_test_packet.deinit();

    std.debug.print("\nCCSDS Packet Created:\n{any}", .{converted_test_packet});
}

```

#### precess a star to July 30, 2005

```zig
const std = @import("std");
const astroz = @import("astroz");
const coordinates = astroz.coordinates;
const Datetime = astroz.time.Datetime;

pub fn main() !void {
    const declination = coordinates.Declination.new(40, 10, 10);
    const ra = coordinates.Right_Ascension.new(19, 52, 2);
    const j2000 = coordinates.Equatorial_Coordinate_System.new(declination, ra);

    std.debug.print("Precessed to July 30, 2005:\n{any}", .{j2000.precess(Datetime.new_date(2005, 7, 30))});
}

```


