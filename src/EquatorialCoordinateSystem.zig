//! Equatorial Coordinate System is a commonly used coordinate system in astronomy as it doesnt rely on
//! the position of the viewer. Allowing you to find what you're looking easily without conversion

const std = @import("std");

const calculations = @import("calculations.zig");
const constants = @import("constants.zig");
const Datetime = @import("Datetime.zig");

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
    const raRad = self.rightAscension.convertToAngular() * constants.deg2rad;
    const decRad = self.declination.convertToAngular() * constants.deg2rad;
    const raSin = @sin(raRad);
    const decTan = @tan(decRad);
    const raCos = @cos(raRad);

    const deltaRa = precessConstants.M + (precessConstants.N * raSin * decTan);
    const deltaDec = precessConstants.N * raCos;

    const returnRa = (deltaRa * constants.arcsecondsPerDegree) / constants.degreesPerHour;
    const returnDec = deltaDec * constants.arcsecondsPerDegree;

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
        const arcminutes = @as(f32, @floatFromInt(self.arcminutes)) / constants.arcminutesPerDegree;
        const arcseconds = self.arcseconds / constants.arcsecondsPerDegree;
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
        const minutes = @as(f32, @floatFromInt(self.minutes)) / constants.minutesPerHour;
        const seconds = self.seconds / constants.secondsPerHour;
        return (hours + minutes + seconds) * constants.degreesPerHour;
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
    const testRa = RightAscension.init(19, 50, 47.0);
    const testDec = Declination.init(8, 52, 6.0);
    const testCoord = EquatorialCoordinateSystem.init(
        testDec,
        testRa,
    );
    const angularRa = testRa.convertToAngular();
    const angularDec = testDec.convertToAngular();
    const precessedOutput = testCoord.precess(Datetime.initDate(2005, 7, 30));
    const expectedPrecessed = EquatorialCoordinateSystem.init(
        Declination.init(8, 52, 57.962516014965246),
        RightAscension.init(19, 51, 3.122949149012854),
    );

    try std.testing.expectEqual(testCoord.rightAscension, testRa);
    try std.testing.expectEqual(testCoord.declination, testDec);
    try std.testing.expectApproxEqAbs(297.6958428700765, angularRa, 1e-4);
    try std.testing.expectApproxEqAbs(8.86833346048991, angularDec, 1e-4);

    // Check precessed coordinates are within tolerance
    try std.testing.expectApproxEqAbs(expectedPrecessed.declination.arcseconds, precessedOutput.declination.arcseconds, 1e-4);
    try std.testing.expectApproxEqAbs(expectedPrecessed.rightAscension.seconds, precessedOutput.rightAscension.seconds, 1e-4);
}

test "Precess" {
    const testDt = Datetime.initDate(2005, 7, 30);
    const testPre = Precess.init(testDt);

    try std.testing.expectEqual(2005.5780821917808, testPre.t);
    try std.testing.expectEqual(0.055780821917808227, testPre.T);
    try std.testing.expectEqual(0.07146939946550677, testPre.M);
    try std.testing.expectEqual(0.03105576921912479, testPre.N);
}
