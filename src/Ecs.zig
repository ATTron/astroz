//! Equatorial Coordinate System is a commonly used coordinate system in astronomy as it doesnt rely on
//! the position of the viewer. Allowing you to find what you're looking easily without conversion

const std = @import("std");
const Datetime = @import("Datetime.zig");
const constants = @import("constants.zig");
const calculations = @import("calculations.zig");

const Ecs = @This();

declination: Declination,
right_ascension: RightAscension,

pub fn init(declination: Declination, right_ascension: RightAscension) Ecs {
    return .{
        .declination = declination,
        .right_ascension = right_ascension,
    };
}

// TODO: this can 100% be done better
/// Find anything at any date after Jan 1 2000 in the ECS format
pub fn precess(self: Ecs, date: Datetime) Ecs {
    const precess_constants = Precess.init(date);
    const deltas = self.calculateRaDec(precess_constants);

    const total_seconds = self.right_ascension.seconds + deltas.ra;
    const new_seconds = @mod(total_seconds, 60.0);
    const total_minutes = self.right_ascension.minutes + @as(u16, @intFromFloat(total_seconds)) / 60;
    const new_minutes = @mod(total_minutes, 60);
    const new_hours = self.right_ascension.hours + total_minutes / 60;

    const precessed_right_ascension = RightAscension.init(
        new_hours,
        new_minutes,
        new_seconds,
    );

    const total_arcseconds = self.declination.arcseconds + deltas.dec;
    const new_arcseconds = @mod(total_arcseconds, 60.0);
    const total_arcminutes = self.declination.arcminutes + @as(u16, @intFromFloat(total_arcseconds)) / 60;
    const new_arcminutes = @mod(total_arcminutes, 60);
    const new_degrees = self.declination.degrees + total_arcminutes / 60;

    const precessed_declination = Declination.init(
        new_degrees,
        new_arcminutes,
        new_arcseconds,
    );

    return .{
        .declination = precessed_declination,
        .right_ascension = precessed_right_ascension,
    };
}

fn calculateRaDec(self: Ecs, precess_constants: Precess) struct { ra: f64, dec: f64 } {
    const ra_sin = @sin(calculations.degreesToRadians(self.right_ascension.convertToAngular()));
    const dec_tan = @tan(calculations.degreesToRadians(self.declination.convertToAngular()));
    const ra_cos = @cos(calculations.degreesToRadians(self.right_ascension.convertToAngular()));

    const delta_ra = precess_constants.M + (precess_constants.N * ra_sin * dec_tan);
    const delta_dec = precess_constants.N * ra_cos;

    const return_ra = (delta_ra * 3600) / 15;
    const return_dec = delta_dec * 3600;

    return .{ .ra = return_ra, .dec = return_dec };
}

/// Declination portion of ECS.
pub const Declination = struct {
    degrees: u16,
    arcminutes: u16,
    arcseconds: f64,

    pub fn init(degrees: ?u16, arcminutes: ?u16, arcseconds: ?f64) Declination {
        return .{
            .degrees = degrees orelse 0,
            .arcminutes = arcminutes orelse 0,
            .arcseconds = arcseconds orelse 0.0,
        };
    }

    /// Convert to angular degrees, needed to precess
    pub fn convertToAngular(self: Declination) f64 {
        const degrees = @as(f32, @floatFromInt(self.degrees));
        const arcminutes = @as(f32, @floatFromInt(self.arcminutes)) / 60;
        const arcseconds = self.arcseconds / 3600;
        return degrees + arcminutes + arcseconds;
    }
};

/// Right ascension portion of ECS.
pub const RightAscension = struct {
    hours: u16,
    minutes: u16,
    seconds: f64,

    pub fn init(hours: ?u16, minutes: ?u16, seconds: ?f64) RightAscension {
        return .{
            .hours = hours orelse 0,
            .minutes = minutes orelse 0,
            .seconds = seconds orelse 0.0,
        };
    }

    /// Convert to angular degrees, needed to precess
    pub fn convertToAngular(self: RightAscension) f64 {
        const hours = @as(f32, @floatFromInt(self.hours));
        const minutes = @as(f32, @floatFromInt(self.minutes)) / 60.0;
        const seconds = self.seconds / 3600.0;
        return (hours + minutes + seconds) * 15;
    }
};

const Precess = struct {
    datetime: Datetime,
    t: f64,
    T: f64,
    M: f64,
    N: f64,

    pub fn init(datetime: Datetime) Precess {
        const t = @as(f64, @floatFromInt(datetime.year.?)) + calculateDoyPercentage(datetime);
        const T = (t - constants.j2k) / 100.0;
        const M = (1.2812323 * T) + (0.0003879 * std.math.pow(f64, T, 2.0)) + (0.0000101 * std.math.pow(f64, T, 3.0));
        const N = 0.5567530 * T - 0.0001185 * std.math.pow(f64, T, 2.0) - 0.0000116 * std.math.pow(f64, T, 3.0);

        return .{
            .datetime = datetime,
            .t = t,
            .T = T,
            .M = M,
            .N = N,
        };
    }

    fn calculateDoyPercentage(datetime: Datetime) f64 {
        return @as(f64, @floatFromInt(datetime.doy.?)) / @as(f64, @floatFromInt(datetime.days_in_year));
    }
};

test "Equatorial Coordinates" {
    const test_ra = RightAscension.init(19, 50, 47.0);
    const test_dec = Declination.init(8, 52, 6.0);
    const test_coord = Ecs.init(
        test_dec,
        test_ra,
    );
    const angular_ra = test_ra.convertToAngular();
    const angular_dec = test_dec.convertToAngular();
    const precessed_output = test_coord.precess(Datetime.initDate(2005, 7, 30));
    const expected_precessed = Ecs.init(
        Declination.init(8, 52, 57.962516014965246),
        RightAscension.init(19, 51, 3.122949149012854),
    );

    try std.testing.expectEqual(test_coord.right_ascension, test_ra);
    try std.testing.expectEqual(test_coord.declination, test_dec);
    try std.testing.expectEqual(297.6958428700765, angular_ra);
    try std.testing.expectEqual(8.86833346048991, angular_dec);
    try std.testing.expectEqual(expected_precessed, precessed_output);
}

test "Precess" {
    const test_dt = Datetime.initDate(2005, 7, 30);
    const test_pre = Precess.init(test_dt);

    try std.testing.expectEqual(2005.5780821917808, test_pre.t);
    try std.testing.expectEqual(0.055780821917808227, test_pre.T);
    try std.testing.expectEqual(0.07146939946550677, test_pre.M);
    try std.testing.expectEqual(0.03105576921912479, test_pre.N);
}
