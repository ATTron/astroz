//! Custom datetime object for dealing with a variety of datetime formats.

const std = @import("std");
const constants = @import("constants.zig");

const Datetime = @This();

const days_in_month = [_]u8{ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
const days_per_year = 365;

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
        .daysInYear = days_per_year,
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
        .daysInYear = days_per_year,
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
        .daysInYear = days_per_year,
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
        .daysInYear = days_per_year,
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
    var dim = days_in_month;
    if (isLeapYear(self.year.?)) {
        dim[1] = 29;
        self.daysInYear = days_per_year + 1;
    }
    var doy: u16 = 0;
    for (dim, 1..) |days, i| {
        if (i == self.month.?) {
            doy += self.day.?;
            break;
        }
        doy += days;
    }
    self.doy = doy;
}

fn epochToDatetime(timestamp: i64) Datetime {
    var remainingSeconds: i64 = timestamp;
    var year: u16 = 1970;

    while (true) {
        const daysThisYear: u32 = if (isLeapYear(year)) days_per_year + 1 else days_per_year;
        const secondsThisYear = daysThisYear * std.time.s_per_day;
        if (remainingSeconds < secondsThisYear) break;
        remainingSeconds -= secondsThisYear;
        year += 1;
    }

    const daysElapsed = @divFloor(remainingSeconds, std.time.s_per_day);
    remainingSeconds -= daysElapsed * std.time.s_per_day;

    var dim = days_in_month;
    if (isLeapYear(year)) dim[1] = 29;

    var month: u8 = 1;
    var day: u8 = 1;
    var daysRemaining = @as(u16, @intCast(daysElapsed)) + 1;
    for (dim, 0..) |days, i| {
        if (daysRemaining <= days) {
            month = @as(u8, @intCast(i + 1));
            day = @as(u8, @intCast(daysRemaining));
            break;
        }
        daysRemaining -= days;
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
    var month: u8 = 1;
    var day = doy;

    for (days_in_month) |days| {
        const days_f64: f64 = @floatFromInt(if (month == 2 and isLeapYear(year)) days + 1 else days);
        if (day > days_f64) {
            day -= days_f64;
            month += 1;
        } else {
            break;
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

/// Converts to Julian Date (days since Jan 1, 4713 BC)
pub fn toJulianDate(self: Datetime) f64 {
    const y: f64 = @floatFromInt(self.year.?);
    const m: f64 = @floatFromInt(self.month.?);
    const d: f64 = @floatFromInt(self.day.?);

    const a = @floor((14.0 - m) / 12.0);
    const yy = y + 4800.0 - a;
    const mm = m + 12.0 * a - 3.0;

    var jd = d + @floor((153.0 * mm + 2.0) / 5.0) + 365.0 * yy +
        @floor(yy / 4.0) - @floor(yy / 100.0) + @floor(yy / 400.0) - 32045.0;

    // add fractional day from time if present
    if (self.hours != null) {
        const h: f64 = @floatFromInt(self.hours.?);
        const min: f64 = if (self.minutes) |mins| @floatFromInt(mins) else 0.0;
        const sec: f64 = if (self.seconds) |secs| @floatCast(secs) else 0.0;
        jd += (h - 12.0) / constants.hours_per_day + min / constants.minutes_per_day + sec / constants.seconds_per_day;
    }

    return jd;
}

/// initialize from year and fractional day-of-year (TLE epoch format)
pub fn fromYearDoy(year: u16, doy: f64) Datetime {
    const md = doyToMonthDay(year, doy);
    const fractional_day = doy - @floor(doy);
    const total_seconds = fractional_day * constants.seconds_per_day;
    const hours: u8 = @intFromFloat(@floor(total_seconds / constants.seconds_per_hour));
    const minutes: u8 = @intFromFloat(@floor(@mod(total_seconds, constants.seconds_per_hour) / constants.seconds_per_minute));
    const seconds: f16 = @floatCast(@mod(total_seconds, constants.seconds_per_minute));

    return Datetime.initDatetime(year, md.month, md.day, hours, minutes, seconds);
}

/// Convert year and fractional day-of-year directly to Julian Date
pub fn yearDoyToJulianDate(year: u16, doy: f64) f64 {
    const y: f64 = @floatFromInt(year);
    const a = @floor((14.0 - 1.0) / 12.0);
    const yy = y + 4800.0 - a;
    const mm = 1.0 + 12.0 * a - 3.0;
    const jd_jan1 = 1.0 + @floor((153.0 * mm + 2.0) / 5.0) + 365.0 * yy +
        @floor(yy / 4.0) - @floor(yy / 100.0) + @floor(yy / 400.0) - 32045.0;
    return jd_jan1 + doy - 1.0;
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
