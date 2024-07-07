const std = @import("std");

pub const Datetime = struct {
    instant: ?std.time.Instant,
    doy: ?u16,
    days_in_year: u16,
    year: ?u16,
    month: ?u8,
    day: ?u8,
    hours: ?u8,
    minutes: ?u8,
    seconds: ?f16,

    const Self = @This();

    pub fn new_date(year: u16, month: u8, day: u8) Self {
        var dt = Datetime{ .instant = null, .doy = null, .days_in_year = 365, .year = year, .month = month, .day = day, .hours = null, .minutes = null, .seconds = null };
        dt.calculate_doy();
        return dt;
    }

    pub fn new_time(hours: u8, minutes: u8, seconds: f16) Self {
        return .{ .instant = null, .doy = null, .days_in_year = 365, .year = null, .month = null, .day = null, .hours = hours, .minutes = minutes, .seconds = seconds };
    }

    pub fn new_datetime(year: u16, month: u8, day: u8, hours: u8, minutes: u8, seconds: f16) Self {
        var dt = Datetime{ .instant = null, .doy = null, .days_in_year = 365, .year = year, .month = month, .day = day, .hours = hours, .minutes = minutes, .seconds = seconds };
        dt.calculate_doy();
        return dt;
    }

    pub fn from_instant(instant: std.time.Instant) Self {
        var new_dt = Self{ .instant = instant, .doy = null, .days_in_year = 365, .year = null, .month = null, .day = null, .hours = null, .minutes = null, .seconds = null };
        return new_dt.epoch_to_datetime();
    }

    fn calculate_doy(self: *Self) void {
        var days_in_month = [_]u8{ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
        if (isLeapYear(self.year.?)) {
            days_in_month[1] = 29;
            self.days_in_year += 1;
        }
        var doy: u16 = 0;
        for (days_in_month, 1..) |days, i| {
            if (i == self.month.?) {
                doy += self.day.?;
                break;
            }
            doy += days;
        }
        self.doy = doy;
    }

    fn epoch_to_datetime(comptime timestamp: i64) Self {
        const days_per_year = 365;

        var remaining_seconds: i64 = timestamp;
        var year: u16 = 1970;

        while (true) {
            const days_this_year: u32 = if (isLeapYear(year)) days_per_year + 1 else days_per_year;
            const seconds_this_year = days_this_year * std.time.s_per_day;
            if (remaining_seconds < seconds_this_year) break;
            remaining_seconds -= seconds_this_year;
            year += 1;
        }

        const days_elapsed = @divFloor(remaining_seconds, std.time.s_per_day);
        remaining_seconds -= days_elapsed * std.time.s_per_day;

        var days_in_month = [_]u8{ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
        if (isLeapYear(year)) days_in_month[1] = 29;

        var month: u8 = 1;
        var day: u8 = 1;
        var doy: u16 = 1;
        var days_remaining = @as(u16, @intCast(days_elapsed)) + 1;
        for (days_in_month, 0..) |days, i| {
            if (days_remaining <= days) {
                month = @as(u8, @intCast(i + 1));
                day = @as(u8, @intCast(days_remaining));
                break;
            }
            days_remaining -= days;
            doy += 1;
        }

        const hours = @as(u8, @intCast(@divFloor(remaining_seconds, std.time.s_per_hour)));
        remaining_seconds -= hours * @as(i64, std.time.s_per_hour);
        const minutes = @as(u8, @intCast(@divFloor(remaining_seconds, std.time.s_per_min)));
        const seconds = @mod(remaining_seconds, std.time.s_per_min);

        return Datetime.new_datetime(year, month, day, hours, minutes, @as(f16, @floatFromInt(seconds)));
    }

    fn isLeapYear(year: u16) bool {
        return (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0);
    }

    pub fn doy_to_month_day(year: u16, doy: f64) struct { month: u8, day: u8 } {
        var days_in_month = [_]f64{
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
        if (isLeapYear(year)) days_in_month[1] = 29.0;
        var month: u8 = 1;
        var day = doy;

        for (days_in_month) |days| {
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

    pub fn convert_to_j2000(self: Self) f32 {
        const step_1 = 367.0 * @as(f32, @floatFromInt(self.year.?));
        const step_2 = @as(f32, @floatFromInt(self.year.?)) + @floor((@as(f32, @floatFromInt(self.month.?)) + 9.0) / 12.0);
        const step_3 = @as(f32, @floatFromInt(self.month.?)) * step_2;
        const step_4 = @floor(step_3 / 4.0);
        const step_5 = 275.0 * @as(f32, @floatFromInt(self.month.?));
        const step_6 = @floor(step_5 / 9.0);
        const step_7 = step_1 - step_4 + step_6 + @as(f32, @floatFromInt(self.day.?)) + 1721013.5;

        return step_7;
    }

    pub fn convert_to_modified_jd(self: Self) f32 {
        return self.convert_to_j2000() - 2400000.5;
    }
};

test "Test Date" {
    const dt = Datetime.new_date(2024, 6, 24);

    try std.testing.expectEqual(2024, dt.year);
    try std.testing.expectEqual(6, dt.month);
    try std.testing.expectEqual(24, dt.day);
    try std.testing.expectEqual(null, dt.hours);
    try std.testing.expectEqual(null, dt.minutes);
    try std.testing.expectEqual(null, dt.seconds);
}

test "Test Time" {
    const ts = Datetime.new_time(16, 6, 24);

    try std.testing.expectEqual(null, ts.year);
    try std.testing.expectEqual(null, ts.month);
    try std.testing.expectEqual(null, ts.day);
    try std.testing.expectEqual(16, ts.hours);
    try std.testing.expectEqual(6, ts.minutes);
    try std.testing.expectEqual(24, ts.seconds);
}

test "Test Datetime" {
    const dt = Datetime.new_datetime(2005, 6, 30, 16, 7, 45);

    try std.testing.expectEqual(2005, dt.year);
    try std.testing.expectEqual(6, dt.month);
    try std.testing.expectEqual(30, dt.day);
    try std.testing.expectEqual(16, dt.hours);
    try std.testing.expectEqual(7, dt.minutes);
    try std.testing.expectEqual(45, dt.seconds);
    try std.testing.expectEqual(181, dt.doy);
}

test "Test Datetime Functions" {
    const dt = Datetime.epoch_to_datetime(800077635);

    try std.testing.expectEqual(1995, dt.year);
    try std.testing.expectEqual(5, dt.month);
    try std.testing.expectEqual(10, dt.day);
    try std.testing.expectEqual(3, dt.hours);
    try std.testing.expectEqual(47, dt.minutes);
    try std.testing.expectEqual(15, dt.seconds);
}

test "Test J2000" {
    const j2000 = Datetime.new_date(2005, 7, 30);

    try std.testing.expectEqual(2453581.5, j2000.convert_to_j2000());
    try std.testing.expectEqual(53581.0, j2000.convert_to_modified_jd());
}
