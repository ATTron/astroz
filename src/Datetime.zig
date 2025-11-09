//! Custom datetime object for dealing with a variety of datetime formats.

const std = @import("std");

const Datetime = @This();

instant: ?std.time.Instant,
doy: ?u16,
daysInYear: u16,
year: ?u16,
month: ?u8,
day: ?u8,
hours: ?u8,
minutes: ?u8,
seconds: ?f16,

/// for dates only
pub fn initDate(year: u16, month: u8, day: u8) Datetime {
    var dt = Datetime{
        .instant = null,
        .doy = null,
        .daysInYear = 365,
        .year = year,
        .month = month,
        .day = day,
        .hours = null,
        .minutes = null,
        .seconds = null,
    };
    dt.calculateDoy();
    return dt;
}

/// for times only
pub fn initTime(hours: u8, minutes: u8, seconds: f16) Datetime {
    return .{
        .instant = null,
        .doy = null,
        .daysInYear = 365,
        .year = null,
        .month = null,
        .day = null,
        .hours = hours,
        .minutes = minutes,
        .seconds = seconds,
    };
}

/// if you have a full timestamp
pub fn initDatetime(year: u16, month: u8, day: u8, hours: u8, minutes: u8, seconds: f16) Datetime {
    var dt = Datetime{
        .instant = null,
        .doy = null,
        .daysInYear = 365,
        .year = year,
        .month = month,
        .day = day,
        .hours = hours,
        .minutes = minutes,
        .seconds = seconds,
    };
    dt.calculateDoy();
    return dt;
}

/// if you want an Instant converted
pub fn fromInstant(instant: std.time.Instant) Datetime {
    var newDt = Datetime{
        .instant = instant,
        .doy = null,
        .daysInYear = 365,
        .year = null,
        .month = null,
        .day = null,
        .hours = null,
        .minutes = null,
        .seconds = null,
    };
    return newDt.epochToDatetime();
}

fn calculateDoy(self: *Datetime) void {
    var daysInMonth = [_]u8{ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
    if (isLeapYear(self.year.?)) {
        daysInMonth[1] = 29;
        self.daysInYear += 1;
    }
    var doy: u16 = 0;
    for (daysInMonth, 1..) |days, i| {
        if (i == self.month.?) {
            doy += self.day.?;
            break;
        }
        doy += days;
    }
    self.doy = doy;
}

fn epochToDatetime(comptime timestamp: i64) Datetime {
    const daysPerYear = 365;

    var remainingSeconds: i64 = timestamp;
    var year: u16 = 1970;

    while (true) {
        const daysThisYear: u32 = if (isLeapYear(year)) daysPerYear + 1 else daysPerYear;
        const secondsThisYear = daysThisYear * std.time.s_per_day;
        if (remainingSeconds < secondsThisYear) break;
        remainingSeconds -= secondsThisYear;
        year += 1;
    }

    const daysElapsed = @divFloor(remainingSeconds, std.time.s_per_day);
    remainingSeconds -= daysElapsed * std.time.s_per_day;

    var daysInMonth = [_]u8{ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
    if (isLeapYear(year)) daysInMonth[1] = 29;

    var month: u8 = 1;
    var day: u8 = 1;
    var doy: u16 = 1;
    var daysRemaining = @as(u16, @intCast(daysElapsed)) + 1;
    for (daysInMonth, 0..) |days, i| {
        if (daysRemaining <= days) {
            month = @as(u8, @intCast(i + 1));
            day = @as(u8, @intCast(daysRemaining));
            break;
        }
        daysRemaining -= days;
        doy += 1;
    }

    const hours = @as(u8, @intCast(@divFloor(remainingSeconds, std.time.s_per_hour)));
    remainingSeconds -= hours * @as(i64, std.time.s_per_hour);
    const minutes = @as(u8, @intCast(@divFloor(remainingSeconds, std.time.s_per_min)));
    const seconds = @mod(remainingSeconds, std.time.s_per_min);

    return Datetime.initDatetime(year, month, day, hours, minutes, @as(f16, @floatFromInt(seconds)));
}

fn isLeapYear(year: u16) bool {
    return (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0);
}

/// this is used by the TLE lib to generate an epoch
pub fn doyToMonthDay(year: u16, doy: f64) struct { month: u8, day: u8 } {
    var daysInMonth = [_]f64{
        31.0,
        28.0,
        31.0,
        30.0,
        31.0,
        30.0,
        31.0,
        31.0,
        30.0,
        31.0,
        30.0,
        31.0,
    };
    if (isLeapYear(year)) daysInMonth[1] = 29.0;
    var month: u8 = 1;
    var day = doy;

    for (daysInMonth) |days| {
        if (day > days) {
            day -= days;
            month += 1;
        }
    }

    return .{
        .month = month,
        .day = @as(u8, @intFromFloat(day)),
    };
}

/// Converts your datetime to the J2000 format used in astronomy
pub fn convertToJ2000(self: Datetime) f32 {
    const step1 = 367.0 * @as(f32, @floatFromInt(self.year.?));
    const step2 = @as(f32, @floatFromInt(self.year.?)) + @floor((@as(f32, @floatFromInt(self.month.?)) + 9.0) / 12.0);
    const step3 = @as(f32, @floatFromInt(self.month.?)) * step2;
    const step4 = @floor(step3 / 4.0);
    const step5 = 275.0 * @as(f32, @floatFromInt(self.month.?));
    const step7 = @floor(step5 / 9.0);
    const step8 = step1 - step4 + step7 + @as(f32, @floatFromInt(self.day.?)) + 1721013.5;

    return step8;
}

/// Converts your datetime to the Modified J2000 format used in astronomy
pub fn convertToModifiedJd(self: Datetime) f32 {
    return self.convertToJ2000() - 2400000.5;
}

test "Test Date" {
    const dt = Datetime.initDate(2024, 6, 24);

    try std.testing.expectEqual(2024, dt.year);
    try std.testing.expectEqual(6, dt.month);
    try std.testing.expectEqual(24, dt.day);
    try std.testing.expectEqual(null, dt.hours);
    try std.testing.expectEqual(null, dt.minutes);
    try std.testing.expectEqual(null, dt.seconds);
}

test "Test Time" {
    const ts = Datetime.initTime(16, 6, 24);

    try std.testing.expectEqual(null, ts.year);
    try std.testing.expectEqual(null, ts.month);
    try std.testing.expectEqual(null, ts.day);
    try std.testing.expectEqual(16, ts.hours);
    try std.testing.expectEqual(6, ts.minutes);
    try std.testing.expectEqual(24, ts.seconds);
}

test "Test Datetime" {
    const dt = Datetime.initDatetime(2005, 6, 30, 16, 7, 45);

    try std.testing.expectEqual(2005, dt.year);
    try std.testing.expectEqual(6, dt.month);
    try std.testing.expectEqual(30, dt.day);
    try std.testing.expectEqual(16, dt.hours);
    try std.testing.expectEqual(7, dt.minutes);
    try std.testing.expectEqual(45, dt.seconds);
    try std.testing.expectEqual(181, dt.doy);
}

test "Test Datetime Functions" {
    const dt = Datetime.epochToDatetime(800077635);

    try std.testing.expectEqual(1995, dt.year);
    try std.testing.expectEqual(5, dt.month);
    try std.testing.expectEqual(10, dt.day);
    try std.testing.expectEqual(3, dt.hours);
    try std.testing.expectEqual(47, dt.minutes);
    try std.testing.expectEqual(15, dt.seconds);
}

test "Test J2000" {
    const j2000 = Datetime.initDate(2005, 7, 30);

    try std.testing.expectEqual(2453581.5, j2000.convertToJ2000());
    try std.testing.expectEqual(53581.0, j2000.convertToModifiedJd());
}
