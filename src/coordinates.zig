const std = @import("std");

pub const Equatorial_Coordinate_System = struct {
    declination: Declination,
    right_ascension: Right_Ascension,

    const Self = @This();

    pub fn new(declination: Declination, right_ascension: Right_Ascension) Self {
        return .{ .declination = declination, .right_ascension = right_ascension };
    }

    pub fn precess(self: Self, date: std.time.Instant) void {
        std.debug.print("HEYO IMPLEMENT ME: {any} : {date}", .{ self, date });
    }
};

pub const Declination = struct {
    degrees: u8,
    arcminutes: u8,
    arcseconds: f32,

    const Self = @This();

    pub fn new(degree: ?u8, arcminutes: ?u8, arcseconds: ?f32) Self {
        return .{ .degree = degree orelse 0, .arcminutes = arcminutes orelse 0.0, .arcseconds = arcseconds orelse 0.0 };
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
};
