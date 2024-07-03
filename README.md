# ASTROZ

Astronomical and Spacecraft Toolkit Written in Zig for Zig!  

## Features / Plans

### Spacecraft

- [x] CCSDS Packets
  - [x] CCSDS Stream Parser
- [x] VITA49 Packets
  - [x] Vita49 Stream Parser
- [ ] Orbital Maneuvers
  - [ ] Impulse Maneuvers
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

### Setup Vita49 Parser

#### W/ Callback

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
    _ = try parser.start(callback);
}

fn callback(packet: Vita49) void {
    std.debug.print("Packet recieved: {any}", .{packet});
}
```

#### W/O Callback

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
    _ = try parser.start(null);
}

```

#### Create a CCSDS Packet

##### W/O Config

```zig
const std = @import("std");
const astroz = @import("astroz");
const CCSDS = astroz.ccsds.CCSDS;

pub fn main() !void {
    const raw_test_packet: [16]u8 = .{ 0x78, 0x97, 0xC0, 0x00, 0x00, 0x0A, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A };
    const converted_test_packet = CCSDS.new(&raw_test_packet, null);

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
    const converted_test_packet = CCSDS.new(&raw_test_packet, config);

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
