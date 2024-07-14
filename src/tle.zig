const std = @import("std");
const DateTime = @import("time.zig").Datetime;

/// TLE Error types
/// Only contains a single error, but may be expanded in the future
pub const TleError = error{BadTleLength};

/// The first line of a TLE which comes from passed in TLE
pub const FirstLine = struct {
    line_number: u8,
    satellite_number: u32,
    classification: u8,
    intl_designator: []u8,
    epoch_year: u16,
    epoch_day: f32,
    epoch: f64,
    first_der_m_motion: f32,
    secnd_der_m_motion: [7]u8,
    bstar_drag: f32,
    ephem_type: u8,
    elem_number: u32,
    checksum: u8,

    pub fn init(line: []const u8, allocator: std.mem.Allocator) !FirstLine {
        if (line.len < 69) {
            return TleError.BadTleLength;
        }

        const intl_year = try std.fmt.parseInt(u16, line[9..11], 10);
        const intl_designator = try std.fmt.parseInt(u32, line[11..14], 10);
        const intl_designator_launch = line[14..17].*;

        const string = try std.fmt.allocPrint(
            allocator,
            "{d}-{d}{s}",
            .{ intl_year, intl_designator, intl_designator_launch },
        );

        const drag = try std.fmt.parseFloat(f32, line[54..59]);
        const leading = try std.fmt.parseInt(i32, line[59..61], 10);
        var bstar_drag: f32 = undefined;
        var power: u32 = undefined;

        if (leading < 0) {
            power = @as(u32, @intCast(-leading));
        } else {
            power = @as(u32, @intCast(leading));
        }
        var result: f32 = 1.0;
        for (0..power) |_| {
            result *= 10.0;
        }
        bstar_drag = @abs(-result / drag);
        const epoch_year = try std.fmt.parseInt(u16, line[18..20], 10);
        const epoch_day = try std.fmt.parseFloat(f32, line[20..32]);
        const epoch = tleEpochToEpoch(epoch_year, epoch_day);

        return .{
            .line_number = line[0],
            .satellite_number = try std.fmt.parseInt(u32, line[2..7], 10),
            .classification = line[7],
            .intl_designator = string,
            .epoch_year = epoch_year,
            .epoch_day = epoch_day,
            .epoch = epoch,
            .first_der_m_motion = try std.fmt.parseFloat(f32, line[34..43]),
            .secnd_der_m_motion = line[45..52].*,
            .bstar_drag = bstar_drag,
            .ephem_type = line[62],
            .elem_number = try std.fmt.parseInt(u32, line[65..68], 10),
            .checksum = line[68],
        };
    }

    fn tleEpochToEpoch(year: u16, doy: f64) f32 {
        const epoch_year = 2000 + year;
        const month_day = DateTime.doyToMonthDay(epoch_year, doy);

        const full_epoch = DateTime.initDate(epoch_year, month_day.month, month_day.day).convertToJ2000();

        return full_epoch;
    }
};

/// Second line that comes from passed in TLE
pub const SecondLine = struct {
    line_number: u8,
    satellite_number: u32,
    inclination: f32,
    right_ascension: f32,
    eccentricity: f32,
    perigee: f32,
    m_anomaly: f32,
    m_motion: f64,
    rev_num: f32,
    checksum: u8,
    clean_line: []const u8,

    pub fn init(line: []const u8, allocator: std.mem.Allocator) !SecondLine {
        if (line.len < 69) {
            return TleError.BadTleLength;
        }

        const clean_line = try std.mem.replaceOwned(u8, allocator, line, " ", "");

        const value = try std.fmt.parseFloat(f32, clean_line[21..28]);
        const divisor = std.math.pow(f32, 10, 7);
        const eccentricity: f32 = @as(f32, value) / divisor;

        return .{
            .line_number = line[0],
            .satellite_number = try std.fmt.parseInt(u32, clean_line[1..6], 10),
            .inclination = try std.fmt.parseFloat(f32, clean_line[6..13]),
            .right_ascension = try std.fmt.parseFloat(f32, clean_line[13..21]),
            .eccentricity = eccentricity,
            .perigee = try std.fmt.parseFloat(f32, clean_line[28..36]),
            .m_anomaly = try std.fmt.parseFloat(f32, clean_line[36..43]),
            .m_motion = try std.fmt.parseFloat(f64, clean_line[43..54]),
            .rev_num = try std.fmt.parseFloat(f32, clean_line[54..58]),
            .checksum = clean_line[58],
            .clean_line = clean_line,
        };
    }
};

/// The proper TLE struct
/// You should only call this when you are parsing a TLE
/// The FirstLine and SecondLine are being built by the parse method
pub const TLE = struct {
    first_line: FirstLine,
    second_line: SecondLine,
    allocator: std.mem.Allocator,

    pub fn parse(tle: []const u8, allocator: std.mem.Allocator) !TLE {
        const delimiter: [1]u8 = .{'\n'};
        var lines = std.mem.splitAny(u8, tle, &delimiter);
        const fl = lines.first();
        const sl = lines.next().?;
        const first_line = try FirstLine.init(fl, allocator);
        const second_line = try SecondLine.init(sl, allocator);
        return .{ .first_line = first_line, .second_line = second_line, .allocator = allocator };
    }

    pub fn deinit(self: *TLE) void {
        self.allocator.free(self.first_line.intl_designator);
        self.allocator.free(self.second_line.clean_line);
    }

    /// Helpful for sanity checking TLE parsing
    pub fn output(self: TLE) void {
        // 1st line
        std.debug.print("line_number: {c}\n", .{self.first_line.line_number});
        std.debug.print("satellite_number: {d}\n", .{self.first_line.satellite_number});
        std.debug.print("classification: {c}\n", .{self.first_line.classification});
        std.debug.print("intl_designator: {s}\n", .{self.first_line.intl_designator});
        std.debug.print("epoch_year: {d}\n", .{self.first_line.epoch_year});
        std.debug.print("epoch_day: {d}\n", .{self.first_line.epoch_day});
        std.debug.print("epoch: {d}\n", .{self.first_line.epoch});
        std.debug.print("first_der_m_motion: {d}\n", .{self.first_line.first_der_m_motion});
        std.debug.print("second_der_m_motion: {s}\n", .{self.first_line.secnd_der_m_motion});
        std.debug.print("bstar_drag: {d}\n", .{self.first_line.bstar_drag});
        std.debug.print("ephem_type: {c}\n", .{self.first_line.ephem_type});
        std.debug.print("elem_number: {d}\n", .{self.first_line.elem_number});
        std.debug.print("checksum: {c}\n", .{self.first_line.checksum});

        // 2nd line
        std.debug.print("\n\n", .{});
        std.debug.print("line_number: {c}\n", .{self.second_line.line_number});
        std.debug.print("satellite_number: {d}\n", .{self.second_line.satellite_number});
        std.debug.print("inclination: {d}\n", .{self.second_line.inclination});
        std.debug.print("right_ascension: {d}\n", .{self.second_line.right_ascension});
        std.debug.print("eccentricity: {d}\n", .{self.second_line.eccentricity});
        std.debug.print("perigee: {d}\n", .{self.second_line.perigee});
        std.debug.print("m_anomaly: {d}\n", .{self.second_line.m_anomaly});
        std.debug.print("m_motion: {d}\n", .{self.second_line.m_motion});
        std.debug.print("rev_num: {d}\n", .{self.second_line.rev_num});
        std.debug.print("checksum: {c}\n", .{self.second_line.checksum});
    }
};

test "test tle parsing" {
    const test_tle =
        \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
        \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
    ;

    var tle = try TLE.parse(test_tle, std.testing.allocator);
    defer tle.deinit();

    try std.testing.expectEqual(55909, tle.first_line.satellite_number);
    try std.testing.expectEqual(0.006211566, tle.first_line.bstar_drag);
    try std.testing.expectEqual(0.00023579, tle.first_line.first_der_m_motion);

    try std.testing.expectEqual(43.9978, tle.second_line.inclination);
    try std.testing.expectEqual(311.8012, tle.second_line.right_ascension);
    try std.testing.expectEqual(0.0011446, tle.second_line.eccentricity);
    try std.testing.expectEqual(15.05761711, tle.second_line.m_motion);
}
