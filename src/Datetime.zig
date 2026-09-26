//! Custom datetime object for dealing with a variety of datetime formats.

const std = @import("std");
const constants = @import("constants.zig");

const Datetime = @This();

const daysInMonth = [_]u8{ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
const daysPerYear = 365;

instant: ?std.Io.Timestamp,
doy: ?u16,
daysInYear: u16,
year: ?u16,
month: ?u8,
day: ?u8,
hours: ?u8,
minutes: ?u8,
seconds: ?f64,

/// for dates only
pub fn initDate(year: u16, month: u8, day: u8) Datetime {
    var dt = Datetime{
        .instant = null,
        .doy = null,
        .daysInYear = daysPerYear,
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
pub fn initTime(hours: u8, minutes: u8, seconds: f64) Datetime {
    return .{
        .instant = null,
        .doy = null,
        .daysInYear = daysPerYear,
        .year = null,
        .month = null,
        .day = null,
        .hours = hours,
        .minutes = minutes,
        .seconds = seconds,
    };
}

/// if you have a full timestamp
pub fn initDatetime(year: u16, month: u8, day: u8, hours: u8, minutes: u8, seconds: f64) Datetime {
    var dt = Datetime{
        .instant = null,
        .doy = null,
        .daysInYear = daysPerYear,
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
pub fn fromInstant(instant: std.Io.Timestamp) Datetime {
    var newDt = Datetime{
        .instant = instant,
        .doy = null,
        .daysInYear = daysPerYear,
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
    var dim = daysInMonth;
    if (isLeapYear(self.year.?)) {
        dim[1] = 29;
        self.daysInYear = daysPerYear + 1;
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
        const daysThisYear: u32 = if (isLeapYear(year)) daysPerYear + 1 else daysPerYear;
        const secondsThisYear = daysThisYear * std.time.s_per_day;
        if (remainingSeconds < secondsThisYear) break;
        remainingSeconds -= secondsThisYear;
        year += 1;
    }

    const daysElapsed = @divFloor(remainingSeconds, std.time.s_per_day);
    remainingSeconds -= daysElapsed * std.time.s_per_day;

    var dim = daysInMonth;
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

    return Datetime.initDatetime(year, month, day, hours, minutes, @floatFromInt(seconds));
}

fn isLeapYear(year: u16) bool {
    return (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0);
}

/// Fractional day of year (1.0 = Jan 1 00:00) to calendar month and day
pub fn doyToMonthDay(year: u16, doy: f64) struct { month: u8, day: u8 } {
    var month: u8 = 1;
    var day = doy;

    for (daysInMonth) |days| {
        const daysF64: f64 = @floatFromInt(if (month == 2 and isLeapYear(year)) days + 1 else days);
        // day N runs from N.0 up to N+1, so 31.5 is still Jan 31
        if (day >= daysF64 + 1.0) {
            day -= daysF64;
            month += 1;
        } else {
            break;
        }
    }

    // past Dec 31, keep counting December days (Dec 32, ...), as python-sgp4 does
    if (month == 13) return .{ .month = 12, .day = @as(u8, @trunc(day)) + 31 };
    return .{
        .month = month,
        .day = @as(u8, @trunc(day)),
    };
}

/// Julian date at 00:00 of this date (valid 1901-2099)
pub fn convertToJ2000(self: Datetime) f64 {
    const y: f64 = @floatFromInt(self.year.?);
    const m: f64 = @floatFromInt(self.month.?);
    const d: f64 = @floatFromInt(self.day.?);
    return 367.0 * y - @floor(7.0 * (y + @floor((m + 9.0) / 12.0)) / 4.0) + @floor(275.0 * m / 9.0) + d + 1721013.5;
}

/// Modified Julian date at 00:00 of this date
pub fn convertToModifiedJd(self: Datetime) f64 {
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
        const sec: f64 = self.seconds orelse 0.0;
        jd += (h - 12.0) / constants.hoursPerDay + min / constants.minutesPerDay + sec / constants.secondsPerDay;
    }

    return jd;
}

/// initialize from year and fractional day-of-year (TLE epoch format).
/// Seconds are rounded to the microsecond, the same way python-sgp4's days2mdhms does.
pub fn fromYearDoy(year: u16, doy: f64) Datetime {
    const totalSeconds = roundToMicrosecond(doy * constants.secondsPerDay);
    const totalMinutes = @floor(totalSeconds / constants.secondsPerMinute);
    const seconds = roundToMicrosecond(@mod(totalSeconds, constants.secondsPerMinute));
    const minuteOfDay = @mod(totalMinutes, constants.minutesPerDay);

    const md = doyToMonthDay(year, @floor(totalMinutes / constants.minutesPerDay));
    const hours: u8 = @intFromFloat(@floor(minuteOfDay / 60.0));
    const minutes: u8 = @intFromFloat(@mod(minuteOfDay, 60.0));
    return Datetime.initDatetime(year, md.month, md.day, hours, minutes, seconds);
}

/// Matches Python's round(x, 6): rounds the exact value of `seconds` (not the already
/// rounded product seconds * 1e6) and breaks exact ties to even.
fn roundToMicrosecond(seconds: f64) f64 {
    const scaled = seconds * 1e6;
    const productError = @mulAdd(f64, seconds, 1e6, -scaled);
    var whole = @floor(scaled);
    var frac = (scaled - whole) + productError;
    if (frac < 0) {
        whole -= 1;
        frac += 1;
    } else if (frac >= 1) {
        whole += 1;
        frac -= 1;
    }
    if (frac > 0.5 or (frac == 0.5 and @mod(whole, 2) == 1)) whole += 1;
    return whole / 1e6;
}

/// Convert year and fractional day-of-year directly to Julian Date
/// DOY 1.0 = Jan 1 00:00:00 (midnight), matching TLE epoch convention
pub fn yearDoyToJulianDate(year: u16, doy: f64) f64 {
    const y: f64 = @floatFromInt(year);
    const a = @floor((14.0 - 1.0) / 12.0);
    const yy = y + 4800.0 - a;
    const mm = 1.0 + 12.0 * a - 3.0;
    // jdJan1 is JD at noon on Jan 1; subtract 0.5 to get midnight
    const jdJan1 = 1.0 + @floor((153.0 * mm + 2.0) / 5.0) + 365.0 * yy +
        @floor(yy / 4.0) - @floor(yy / 100.0) + @floor(yy / 400.0) - 32045.0;
    return jdJan1 + doy - 1.5; // -1.5 = -1 for DOY offset, -0.5 for noon→midnight
}

/// python-sgp4 compatible: calendar date and time to a (jd, fr) pair.
/// `jd` is the Julian date at midnight starting the day; `fr` is the fraction of the
/// day, kept separate so it stays precise to well under a microsecond.
pub fn jday(year: u16, month: u8, day: u8, hour: u8, minute: u8, second: f64) struct { jd: f64, fr: f64 } {
    const h: f64 = @floatFromInt(hour);
    const m: f64 = @floatFromInt(minute);
    return .{
        .jd = Datetime.initDate(year, month, day).convertToJ2000(),
        .fr = (second + m * constants.secondsPerMinute + h * constants.secondsPerHour) / constants.secondsPerDay,
    };
}

/// python-sgp4 compatible: fractional doy to (month, day, hour, minute, second)
/// Converts year and fractional day-of-year to calendar components
pub fn days2mdhms(year: u16, days: f64) struct { month: u8, day: u8, hour: u8, minute: u8, second: f64 } {
    const dt = fromYearDoy(year, days);
    return .{
        .month = dt.month.?,
        .day = dt.day.?,
        .hour = dt.hours.?,
        .minute = dt.minutes.?,
        .second = dt.seconds.?,
    };
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
    try std.testing.expectEqual(2461303.5, Datetime.initDate(2026, 9, 20).convertToJ2000());
    try std.testing.expectEqual(2451118.5, Datetime.initDate(1998, 11, 1).convertToJ2000());
    try std.testing.expectEqual(2451544.5, Datetime.initDate(2000, 1, 1).convertToJ2000());
}

test "Test jday" {
    // Test jday: 2019-01-05 04:28:31.5 should give JD ~2458488.686
    const result = Datetime.jday(2019, 1, 5, 4, 28, 31.5);
    // JD at noon should be 2458488.5
    try std.testing.expectEqual(2458488.5, result.jd);
    // Fractional part should be about 0.186... (4:28:31.5 from noon)
    try std.testing.expectApproxEqAbs(0.18647569444444444, result.fr, 1e-9);

    // values from python-sgp4 2.27's sgp4.api.jday
    const cases = [_]struct { [6]f64, f64, f64 }{
        .{ .{ 2026, 9, 1, 0, 3, 59.328 }, 2461284.5, 0.00277 },
        .{ .{ 2020, 2, 11, 13, 57, 0 }, 2458890.5, 0.58125 },
    };
    for (cases) |c| {
        const i = c[0];
        const r = Datetime.jday(@intFromFloat(i[0]), @intFromFloat(i[1]), @intFromFloat(i[2]), @intFromFloat(i[3]), @intFromFloat(i[4]), i[5]);
        try std.testing.expectEqual(c[1], r.jd);
        try std.testing.expectApproxEqAbs(c[2], r.fr, 1e-15);
    }
}

test "Test days2mdhms" {
    // Test days2mdhms: day 5.186475... of 2019 should give Jan 5, 04:28:31.5
    const result = Datetime.days2mdhms(2019, 5.186475694444444);
    try std.testing.expectEqual(1, result.month);
    try std.testing.expectEqual(5, result.day);
    try std.testing.expectEqual(4, result.hour);
    try std.testing.expectEqual(28, result.minute);
    try std.testing.expectApproxEqAbs(31.5, result.second, 0.01);

    // last day of a month must not roll over to day 0 of the next
    const monthEnds = [_]struct { u16, f64, u8, u8 }{
        .{ 2026, 31.5, 1, 31 },   .{ 2026, 32.0, 2, 1 },    .{ 2024, 60.5, 2, 29 },
        .{ 2026, 243.5, 8, 31 },  .{ 2026, 365.5, 12, 31 }, .{ 2026, 366.5, 12, 32 },
        .{ 2024, 367.5, 12, 32 },
    };
    for (monthEnds) |c| {
        const r = Datetime.days2mdhms(c[0], c[1]);
        try std.testing.expectEqual(c[2], r.month);
        try std.testing.expectEqual(c[3], r.day);
    }

    // values from python-sgp4 2.27's sgp4.api.days2mdhms
    const exact = [_]struct { f64, u8, u8, f64 }{
        .{ 244.00277, 0, 3, 59.328 },
        .{ 244.999, 23, 58, 33.6 },
        .{ 244.0 + (3 * 60 + 59.99) / 86400.0, 0, 3, 59.99 },
    };
    for (exact) |c| {
        const r = Datetime.days2mdhms(2026, c[0]);
        try std.testing.expectEqual(c[1], r.hour);
        try std.testing.expectEqual(c[2], r.minute);
        try std.testing.expectApproxEqAbs(c[3], r.second, 1e-9);
    }

    // cases where rounding seconds * 1e6 instead of the exact seconds lands a microsecond off
    const halfMicro = [_]struct { f64, f64 }{
        .{ 352.2400996782234, 44.612198 },
        .{ 302.1565716923206, 27.794216 },
        .{ 106.11012104495948, 34.458284 },
    };
    for (halfMicro) |c| try std.testing.expectEqual(c[1], Datetime.days2mdhms(2026, c[0]).second);

    // every instant of a day comes back within a microsecond, with seconds below 60
    for (0..100_000) |i| {
        const secondOfDay = @as(f64, @floatFromInt(i)) * 0.864;
        const r = Datetime.days2mdhms(2026, 244.0 + secondOfDay / 86400.0);
        const got = @as(f64, @floatFromInt(r.hour)) * 3600 + @as(f64, @floatFromInt(r.minute)) * 60 + r.second;
        try std.testing.expect(r.second < 60);
        try std.testing.expectApproxEqAbs(secondOfDay, got, 1e-6);
    }
}
