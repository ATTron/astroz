//! Equatorial Coordinate System is a commonly used coordinate system in astronomy as it doesnt rely on
//! the position of the viewer. Allowing you to find what you're looking easily without conversion

const std = @import("std");
const Datetime = @import("Datetime.zig");
const constants = @import("constants.zig");
const calculations = @import("calculations.zig");

const EquatorialCoordinateSystem = @This();

declination: Declination,
rightAscension: RightAscension,

pub fn init(declination: Declination, rightAscension: RightAscension) EquatorialCoordinateSystem {
    return .{
        .declination = declination,
        .rightAscension = rightAscension,
    };
}

// TODO: this can 100% be done better
/// Find anything at any date after Jan 1 2000 in the ECS format
pub fn precess(self: EquatorialCoordinateSystem, date: Datetime) EquatorialCoordinateSystem {
    const precessConstants = Precess.init(date);
    const deltas = self.calculateRaDec(precessConstants);

    const totalSeconds = self.rightAscension.seconds + deltas.ra;
    const newSeconds = @mod(totalSeconds, 60.0);
    const totalMinutes = self.rightAscension.minutes + @as(u16, @intFromFloat(totalSeconds)) / 60;
    const newMinutes = @mod(totalMinutes, 60);
    const newHours = self.rightAscension.hours + totalMinutes / 60;

    const precessedRightAscension = RightAscension.init(
        newHours,
        newMinutes,
        newSeconds,
    );

    const totalArcseconds = self.declination.arcseconds + deltas.dec;
    const newArcseconds = @mod(totalArcseconds, 60.0);
    const totalArcminutes = self.declination.arcminutes + @as(u16, @intFromFloat(totalArcseconds)) / 60;
    const newArcminutes = @mod(totalArcminutes, 60);
    const newDegrees = self.declination.degrees + totalArcminutes / 60;

    const precessedDeclination = Declination.init(
        newDegrees,
        newArcminutes,
        newArcseconds,
    );

    return .{
        .declination = precessedDeclination,
        .rightAscension = precessedRightAscension,
    };
}

fn calculateRaDec(self: EquatorialCoordinateSystem, precessConstants: Precess) struct { ra: f64, dec: f64 } {
    const raSin = @sin(calculations.degreesToRadians(self.rightAscension.convertToAngular()));
    const decTan = @tan(calculations.degreesToRadians(self.declination.convertToAngular()));
    const raCos = @cos(calculations.degreesToRadians(self.rightAscension.convertToAngular()));

    const deltaRa = precessConstants.M + (precessConstants.N * raSin * decTan);
    const deltaDec = precessConstants.N * raCos;

    const returnRa = (deltaRa * 3600) / 15;
    const returnDec = deltaDec * 3600;

    return .{ .ra = returnRa, .dec = returnDec };
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

    // TODO: clean up magic numbers
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
        return @as(f64, @floatFromInt(datetime.doy.?)) / @as(f64, @floatFromInt(datetime.daysInYear));
    }
};

test "Equatorial Coordinates" {
    const test_ra = RightAscension.init(19, 50, 47.0);
    const test_dec = Declination.init(8, 52, 6.0);
    const test_coord = EquatorialCoordinateSystem.init(
        test_dec,
        test_ra,
    );
    const angular_ra = test_ra.convertToAngular();
    const angular_dec = test_dec.convertToAngular();
    const precessed_output = test_coord.precess(Datetime.initDate(2005, 7, 30));
    const expected_precessed = EquatorialCoordinateSystem.init(
        Declination.init(8, 52, 57.962516014965246),
        RightAscension.init(19, 51, 3.122949149012854),
    );

    try std.testing.expectEqual(test_coord.rightAscension, test_ra);
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
