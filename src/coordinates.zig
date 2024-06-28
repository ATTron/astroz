const std = @import("std");
const time = @import("time.zig");
const constants = @import("constants.zig");

pub const Equatorial_Coordinate_System = struct {
    declination: Declination,
    right_ascension: Right_Ascension,

    const Self = @This();

    pub fn new(declination: Declination, right_ascension: Right_Ascension) Self {
        return .{ .declination = declination, .right_ascension = right_ascension };
    }

    pub fn precess(self: Self, date: time.Datetime) void {
        const precess_constants = Precess.precess(date);
        const delta_ra = precess_constants.M + (precess_constants.N * std.math.sin(self.right_ascension.convert_to_angular()) * std.math.tan(self.declination.convert_to_angular()));
        const delta_dec = precess_constants.N * std.math.cos(self.right_ascension.convert_to_angular());

        std.debug.print("\ndelta ra: {d}\ndelta dec: {d}\n", .{ delta_ra, delta_dec });
    }
};

pub const Declination = struct {
    degrees: u8,
    arcminutes: u8,
    arcseconds: f32,

    const Self = @This();

    pub fn new(degrees: ?u8, arcminutes: ?u8, arcseconds: ?f32) Self {
        return .{ .degrees = degrees orelse 0, .arcminutes = arcminutes orelse 0.0, .arcseconds = arcseconds orelse 0.0 };
    }

    pub fn convert_to_angular(self: Self) f32 {
        const degrees = @as(f32, @floatFromInt(self.degrees));
        const arcminutes = @as(f32, @floatFromInt(self.arcminutes)) / 60;
        const arcseconds = self.arcseconds / 3600;
        return degrees + arcminutes + arcseconds;
    }
};

pub const Right_Ascension = struct {
    hours: u8,
    minutes: u8,
    seconds: f32,

    const Self = @This();

    pub fn new(hours: ?u8, minutes: ?u8, seconds: ?f32) Self {
        return .{ .hours = hours orelse 0, .minutes = minutes orelse 0, .seconds = seconds orelse 0.0 };
    }

    pub fn convert_to_angular(self: Self) f32 {
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
        return 1.2812323 * self.T + 0.0003879 * std.math.pow(f64, self.T, 2.0) + 0.0000101 * std.math.pow(f64, self.T, 3.0);
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
    test_coord.precess(time.Datetime.new_date(2005, 6, 30));

    try std.testing.expectEqual(test_coord.right_ascension, test_ra);
    try std.testing.expectEqual(test_coord.declination, test_dec);
    try std.testing.expectEqual(297.69586, angular_ra);
    try std.testing.expectEqual(8.868334, angular_dec);
}

test "Precess" {
    const test_dt = time.Datetime.new_date(2005, 6, 30);
    const test_pre = Precess.precess(test_dt);

    try std.testing.expectEqual(2005.495890410959, test_pre.t);
    try std.testing.expectEqual(0.05495890410958964, test_pre.T);
    try std.testing.expectEqual(0.07041629643906712, test_pre.M);
    try std.testing.expectEqual(0.030598174887084096, test_pre.N);
}
