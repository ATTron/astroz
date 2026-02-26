//! Unified satellite type that auto-dispatches between SGP4 and SDP4.

const std = @import("std");
const constants = @import("constants.zig");
const Tle = @import("Tle.zig");
const Sgp4 = @import("Sgp4.zig");
const Sdp4 = @import("Sdp4.zig");

const Satellite = @This();

impl: union(enum) {
    sgp4: Sgp4,
    sdp4: Sdp4,
},

pub fn init(tle: Tle, grav: constants.Sgp4GravityModel) Sgp4.Error!Satellite {
    return .{ .impl = .{ .sgp4 = Sgp4.init(tle, grav) catch |err| switch (err) {
        error.DeepSpaceNotSupported => return .{ .impl = .{ .sdp4 = try Sdp4.init(tle, grav) } },
        else => return err,
    } } };
}

pub fn propagate(self: *const Satellite, tsince: f64) Sgp4.Error![2][3]f64 {
    switch (self.impl) {
        inline else => |*s| return s.propagate(tsince),
    }
}

pub fn propagateN(self: *const Satellite, comptime N: usize, times: [N]f64) Sgp4.Error!Sgp4.PosVelArray(N) {
    switch (self.impl) {
        inline else => |*s| return s.propagateN(N, times),
    }
}

pub fn epochJd(self: *const Satellite) f64 {
    return switch (self.impl) {
        .sgp4 => |s| s.elements.epochJd,
        .sdp4 => |s| s.elements.sgp4.epochJd,
    };
}

pub fn isDeepSpace(self: *const Satellite) bool {
    return self.impl == .sdp4;
}

const testing = std.testing;

test "auto-selects sgp4 for near-earth" {
    const tle_str =
        \\1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995
        \\2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123
    ;
    var tle = try Tle.parse(tle_str, testing.allocator);
    defer tle.deinit();

    const sat = try Satellite.init(tle, constants.wgs84);
    try testing.expect(!sat.isDeepSpace());
}

test "auto-selects sdp4 for deep-space" {
    const tle_str =
        \\1 20413U 90005A   24186.00000000  .00000012  00000+0  10000-3 0  9992
        \\2 20413  55.4408  61.4858 0112981 129.5765 231.5553  2.00561730104446
    ;
    var tle = try Tle.parse(tle_str, testing.allocator);
    defer tle.deinit();

    const sat = try Satellite.init(tle, constants.wgs72);
    try testing.expect(sat.isDeepSpace());
}

test "propagate matches direct sgp4" {
    const tle_str =
        \\1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995
        \\2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123
    ;
    var tle = try Tle.parse(tle_str, testing.allocator);
    defer tle.deinit();

    const sat = try Satellite.init(tle, constants.wgs84);
    const sgp4 = try Sgp4.init(tle, constants.wgs84);

    for ([_]f64{ 0, 30, 60, 120 }) |t| {
        const sat_result = try sat.propagate(t);
        const sgp4_result = try sgp4.propagate(t);
        for (0..3) |i| {
            try testing.expectApproxEqAbs(sgp4_result[0][i], sat_result[0][i], 1e-12);
            try testing.expectApproxEqAbs(sgp4_result[1][i], sat_result[1][i], 1e-12);
        }
    }
}

test "propagate matches direct sdp4" {
    const tle_str =
        \\1 20413U 90005A   24186.00000000  .00000012  00000+0  10000-3 0  9992
        \\2 20413  55.4408  61.4858 0112981 129.5765 231.5553  2.00561730104446
    ;
    var tle = try Tle.parse(tle_str, testing.allocator);
    defer tle.deinit();

    const sat = try Satellite.init(tle, constants.wgs72);
    const sdp4 = try Sdp4.init(tle, constants.wgs72);

    for ([_]f64{ 0, 720, 1440 }) |t| {
        const sat_result = try sat.propagate(t);
        const sdp4_result = try sdp4.propagate(t);
        for (0..3) |i| {
            try testing.expectApproxEqAbs(sdp4_result[0][i], sat_result[0][i], 1e-12);
            try testing.expectApproxEqAbs(sdp4_result[1][i], sat_result[1][i], 1e-12);
        }
    }
}

test "propagateN matches scalar" {
    const tle_str =
        \\1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995
        \\2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123
    ;
    var tle = try Tle.parse(tle_str, testing.allocator);
    defer tle.deinit();

    const sat = try Satellite.init(tle, constants.wgs84);
    const times = [4]f64{ 0, 30, 60, 90 };
    const batch = try sat.propagateN(4, times);

    for (times, batch) |t, batch_result| {
        const scalar = try sat.propagate(t);
        for (0..3) |i| {
            try testing.expectApproxEqAbs(scalar[0][i], batch_result[0][i], 1e-4);
            try testing.expectApproxEqAbs(scalar[1][i], batch_result[1][i], 1e-4);
        }
    }
}

test "epochJd" {
    const tle_str =
        \\1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995
        \\2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123
    ;
    var tle = try Tle.parse(tle_str, testing.allocator);
    defer tle.deinit();

    const sat = try Satellite.init(tle, constants.wgs84);
    const sgp4 = try Sgp4.init(tle, constants.wgs84);
    try testing.expectEqual(sgp4.elements.epochJd, sat.epochJd());
}
