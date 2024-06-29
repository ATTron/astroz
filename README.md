# ASTROZ

Astronomical and Spacecraft Toolkit Written in Zig for Zig!  

## Features / Plans

### Spacecraft

- [x] CCSDS Packets
- [ ] VITA49 Packets

### Astronomical

- [x] Astronomical Constants
- [x] Astronomical Coordinates
- [x] Astronomical Computation
- [ ] Orbital Mechanics

## Install

**Please use the master branch of the zig repository as that is what I'm developing against**

The easiest way I've found to get started with dependencies in zig is the following.  
- in your `main.zig` import the dependency.
```zig
const astroz = @import("astroz");

...
```
- run `zig fetch --save git+https://github.com/ATTron/astroz/#HEAD`  
- inside `build.zig`  
```zig
...
const package = b.dependency("astroz", .{
    .target = target,
    .optimize = optimize,
});

const module = package.module("astroz");

exe.root_module.addImport("astroz", module);

b.installArtifact(exe);

...
```


## Usage

### Example

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
