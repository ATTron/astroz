const std = @import("std");
const time = @import("time.zig");
const constants = @import("constants.zig");
const calculations = @import("calculations.zig");

pub const Equatorial_Coordinate_System = struct {
    declination: Declination,
    right_ascension: Right_Ascension,

    const Self = @This();

    pub fn new(declination: Declination, right_ascension: Right_Ascension) Self {
        return .{ .declination = declination, .right_ascension = right_ascension };
    }

    // this can 100% be done better
    pub fn precess(self: Self, date: time.Datetime) Equatorial_Coordinate_System {
        const precess_constants = Precess.precess(date);
        const deltas = self.calculate_ra_dec(precess_constants);

        std.debug.print("{d}\n{d}", .{ deltas.ra, deltas.dec });

        const total_seconds = self.right_ascension.seconds + deltas.ra;
        const new_seconds = @mod(total_seconds, 60.0);
        const total_minutes = self.right_ascension.minutes + @as(u16, @intFromFloat(total_seconds)) / 60;
        const new_minutes = @mod(total_minutes, 60);
        const new_hours = self.right_ascension.hours + total_minutes / 60;

        const precessed_right_ascension = Right_Ascension.new(new_hours, new_minutes, new_seconds);

        const total_arcseconds = self.declination.arcseconds + deltas.dec;
        const new_arcseconds = @mod(total_arcseconds, 60.0);
        const total_arcminutes = self.declination.arcminutes + @as(u16, @intFromFloat(total_arcseconds)) / 60;
        const new_arcminutes = @mod(total_arcminutes, 60);
        const new_degrees = self.declination.degrees + total_arcminutes / 60;

        const precessed_declination = Declination.new(new_degrees, new_arcminutes, new_arcseconds);

        return .{ .declination = precessed_declination, .right_ascension = precessed_right_ascension };
    }

    fn calculate_ra_dec(self: Self, precess_constants: Precess) struct { ra: f64, dec: f64 } {
        const ra_sin = @sin(calculations.degrees_to_radians(self.right_ascension.convert_to_angular()));
        const dec_tan = @tan(calculations.degrees_to_radians(self.declination.convert_to_angular()));
        const ra_cos = @cos(calculations.degrees_to_radians(self.right_ascension.convert_to_angular()));

        const delta_ra = precess_constants.M + (precess_constants.N * ra_sin * dec_tan);
        const delta_dec = precess_constants.N * ra_cos;

        const return_ra = (delta_ra * 3600) / 15;
        const return_dec = delta_dec * 3600;

        return .{ .ra = return_ra, .dec = return_dec };
    }
};

pub const Declination = struct {
    degrees: u16,
    arcminutes: u16,
    arcseconds: f64,

    const Self = @This();

    pub fn new(degrees: ?u16, arcminutes: ?u16, arcseconds: ?f64) Self {
        return .{ .degrees = degrees orelse 0, .arcminutes = arcminutes orelse 0.0, .arcseconds = arcseconds orelse 0.0 };
    }

    pub fn convert_to_angular(self: Self) f64 {
        const degrees = @as(f32, @floatFromInt(self.degrees));
        const arcminutes = @as(f32, @floatFromInt(self.arcminutes)) / 60;
        const arcseconds = self.arcseconds / 3600;
        return degrees + arcminutes + arcseconds;
    }
};

pub const Right_Ascension = struct {
    hours: u16,
    minutes: u16,
    seconds: f64,

    const Self = @This();

    pub fn new(hours: ?u16, minutes: ?u16, seconds: ?f64) Self {
        return .{ .hours = hours orelse 0, .minutes = minutes orelse 0, .seconds = seconds orelse 0.0 };
    }

    pub fn convert_to_angular(self: Self) f64 {
        const hours = @as(f32, @floatFromInt(self.hours));
        const minutes = @as(f32, @floatFromInt(self.minutes)) / 60.0;
        const seconds = self.seconds / 3600.0;
        return (hours + minutes + seconds) * 15;
    }
};

const Precess = struct {
    datetime: time.Datetime,
    t: f64,
    T: f64,
    M: f64,
    N: f64,

    const Self = @This();

    pub fn precess(date: time.Datetime) Self {
        var pre = Precess{ .datetime = date, .t = 0, .T = 0, .M = 0, .N = 0 };
        pre.t = pre.calculate_t();
        pre.T = pre.calculate_T();
        pre.M = pre.calculate_M();
        pre.N = pre.calculate_N();

        return pre;
    }

    fn calculate_doy_percentage(self: *Self) f64 {
        return @as(f64, @floatFromInt(self.datetime.doy.?)) / @as(f64, @floatFromInt(self.datetime.days_in_year));
    }

    fn calculate_t(self: *Self) f64 {
        return @as(f64, @floatFromInt(self.datetime.year.?)) + self.calculate_doy_percentage();
    }

    fn calculate_T(self: *Self) f64 {
        return (self.calculate_t() - constants.j2k) / 100.0;
    }

    fn calculate_M(self: *Self) f64 {
        return (1.2812323 * self.T) + (0.0003879 * std.math.pow(f64, self.T, 2.0)) + (0.0000101 * std.math.pow(f64, self.T, 3.0));
    }

    fn calculate_N(self: *Self) f64 {
        return 0.5567530 * self.T - 0.0001185 * std.math.pow(f64, self.T, 2.0) - 0.0000116 * std.math.pow(f64, self.T, 3.0);
    }
};

test "Equatorial Coordinates" {
    const test_ra = Right_Ascension.new(19, 50, 47.0);
    const test_dec = Declination.new(8, 52, 6.0);
    const test_coord = Equatorial_Coordinate_System.new(test_dec, test_ra);
    const angular_ra = test_ra.convert_to_angular();
    const angular_dec = test_dec.convert_to_angular();
    const precessed_output = test_coord.precess(time.Datetime.new_date(2005, 7, 30));
    const expected_precessed = Equatorial_Coordinate_System.new(Declination.new(8, 52, 57.962516014965246), Right_Ascension.new(19, 51, 3.122949149012854));

    try std.testing.expectEqual(test_coord.right_ascension, test_ra);
    try std.testing.expectEqual(test_coord.declination, test_dec);
    try std.testing.expectEqual(297.6958428700765, angular_ra);
    try std.testing.expectEqual(8.86833346048991, angular_dec);
    try std.testing.expectEqual(expected_precessed, precessed_output);
}

test "Precess" {
    const test_dt = time.Datetime.new_date(2005, 7, 30);
    const test_pre = Precess.precess(test_dt);

    try std.testing.expectEqual(2005.5780821917808, test_pre.t);
    try std.testing.expectEqual(0.055780821917808227, test_pre.T);
    try std.testing.expectEqual(0.07146939946550677, test_pre.M);
    try std.testing.expectEqual(0.03105576921912479, test_pre.N);
}
