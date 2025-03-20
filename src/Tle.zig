//! The proper TLE struct
//! You should only call this when you are parsing a TLE
//! The FirstLine and SecondLine are being built by the parse method

const std = @import("std");
const DateTime = @import("Datetime.zig");

const Tle = @This();

firstLine: FirstLine,
secondLine: SecondLine,
allocator: std.mem.Allocator,

pub fn parse(tle: []const u8, allocator: std.mem.Allocator) !Tle {
    const delimiter: [1]u8 = .{'\n'};
    var lines = std.mem.splitAny(u8, tle, &delimiter);
    const fl = lines.first();
    const sl = lines.next().?;
    const firstLine = try FirstLine.init(fl, allocator);
    const secondLine = try SecondLine.init(sl, allocator);
    return .{ .firstLine = firstLine, .secondLine = secondLine, .allocator = allocator };
}

pub fn deinit(self: *Tle) void {
    self.allocator.free(self.firstLine.intlDesignator);
    self.allocator.free(self.secondLine.cleanLine);
}

/// Helpful for sanity checking TLE parsing
pub fn output(self: Tle) void {
    // 1st line
    std.debug.print("line_number: {c}\n", .{self.firstLine.lineNumber});
    std.debug.print("satellite_number: {d}\n", .{self.firstLine.satelliteNumber});
    std.debug.print("classification: {c}\n", .{self.firstLine.classification});
    std.debug.print("intl_designator: {s}\n", .{self.firstLine.intlDesignator});
    std.debug.print("epoch_year: {d}\n", .{self.firstLine.epochYear});
    std.debug.print("epoch_day: {d}\n", .{self.firstLine.epochDay});
    std.debug.print("epoch: {d}\n", .{self.firstLine.epoch});
    std.debug.print("first_der_m_motion: {d}\n", .{self.firstLine.firstDerMeanMotion});
    std.debug.print("second_der_m_motion: {s}\n", .{self.firstLine.secondDerMeanMotion});
    std.debug.print("bstar_drag: {d}\n", .{self.firstLine.bstarDrag});
    std.debug.print("ephem_type: {c}\n", .{self.firstLine.ephemType});
    std.debug.print("elem_number: {d}\n", .{self.firstLine.elemNumber});
    std.debug.print("checksum: {c}\n", .{self.firstLine.checksum});

    // 2nd line
    std.debug.print("\n\n", .{});
    std.debug.print("line_number: {c}\n", .{self.secondLine.lineNumber});
    std.debug.print("satellite_number: {d}\n", .{self.secondLine.satelliteNumber});
    std.debug.print("inclination: {d}\n", .{self.secondLine.inclination});
    std.debug.print("right_ascension: {d}\n", .{self.secondLine.rightAscension});
    std.debug.print("eccentricity: {d}\n", .{self.secondLine.eccentricity});
    std.debug.print("perigee: {d}\n", .{self.secondLine.perigee});
    std.debug.print("m_anomaly: {d}\n", .{self.secondLine.mAnomaly});
    std.debug.print("m_motion: {d}\n", .{self.secondLine.mMotion});
    std.debug.print("rev_num: {d}\n", .{self.secondLine.revNum});
    std.debug.print("checksum: {c}\n", .{self.secondLine.checksum});
}

/// TLE Error types
/// Only contains a single error, but may be expanded in the future
pub const Error = error{BadTleLength};

/// The first line of a TLE which comes from passed in TLE
pub const FirstLine = struct {
    lineNumber: u8,
    satelliteNumber: u32,
    classification: u8,
    intlDesignator: []u8,
    epochYear: u16,
    epochDay: f32,
    epoch: f64,
    firstDerMeanMotion: f32,
    secondDerMeanMotion: [7]u8,
    bstarDrag: f32,
    ephemType: u8,
    elemNumber: u32,
    checksum: u8,

    pub fn init(line: []const u8, allocator: std.mem.Allocator) !FirstLine {
        if (line.len < 69) {
            return Error.BadTleLength;
        }

        const intlYear = try std.fmt.parseInt(u16, line[9..11], 10);
        const intlDesignator = try std.fmt.parseInt(u32, line[11..14], 10);
        const intlDesignatorLaunch = line[14..17].*;

        const string = try std.fmt.allocPrint(
            allocator,
            "{d}-{d}{s}",
            .{ intlYear, intlDesignator, intlDesignatorLaunch },
        );

        const drag = try std.fmt.parseFloat(f32, line[54..59]);
        const leading = try std.fmt.parseInt(i32, line[59..61], 10);
        var bstarDrag: f32 = undefined;
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
        bstarDrag = @abs(-result / drag);
        const epochYear = try std.fmt.parseInt(u16, line[18..20], 10);
        const epochDay = try std.fmt.parseFloat(f32, line[20..32]);
        const epoch = tleEpochToEpoch(epochYear, epochDay);

        return .{
            .lineNumber = line[0],
            .satelliteNumber = try std.fmt.parseInt(u32, line[2..7], 10),
            .classification = line[7],
            .intlDesignator = string,
            .epochYear = epochYear,
            .epochDay = epochDay,
            .epoch = epoch,
            .firstDerMeanMotion = try std.fmt.parseFloat(f32, line[34..43]),
            .secondDerMeanMotion = line[45..52].*,
            .bstarDrag = bstarDrag,
            .ephemType = line[62],
            .elemNumber = try std.fmt.parseInt(u32, line[65..68], 10),
            .checksum = line[68],
        };
    }

    fn tleEpochToEpoch(year: u16, doy: f64) f32 {
        const epochYear = 2000 + year;
        const monthDay = DateTime.doyToMonthDay(epochYear, doy);

        const fullEpoch = DateTime.initDate(epochYear, monthDay.month, monthDay.day).convertToJ2000();

        return fullEpoch;
    }
};

/// Second line that comes from passed in TLE
pub const SecondLine = struct {
    lineNumber: u8,
    satelliteNumber: u32,
    inclination: f32,
    rightAscension: f32,
    eccentricity: f32,
    perigee: f32,
    mAnomaly: f32,
    mMotion: f64,
    revNum: f32,
    checksum: u8,
    cleanLine: []const u8,

    pub fn init(line: []const u8, allocator: std.mem.Allocator) !SecondLine {
        if (line.len < 69) {
            return Error.BadTleLength;
        }

        const cleanLine = try std.mem.replaceOwned(u8, allocator, line, " ", "");

        const value = try std.fmt.parseFloat(f32, cleanLine[21..28]);
        const divisor = std.math.pow(f32, 10, 7);
        const eccentricity: f32 = @as(f32, value) / divisor;

        return .{
            .lineNumber = line[0],
            .satelliteNumber = try std.fmt.parseInt(u32, cleanLine[1..6], 10),
            .inclination = try std.fmt.parseFloat(f32, cleanLine[6..13]),
            .rightAscension = try std.fmt.parseFloat(f32, cleanLine[13..21]),
            .eccentricity = eccentricity,
            .perigee = try std.fmt.parseFloat(f32, cleanLine[28..36]),
            .mAnomaly = try std.fmt.parseFloat(f32, cleanLine[36..43]),
            .mMotion = try std.fmt.parseFloat(f64, cleanLine[43..54]),
            .revNum = try std.fmt.parseFloat(f32, cleanLine[54..58]),
            .checksum = cleanLine[58],
            .cleanLine = cleanLine,
        };
    }
};

test "test tle parsing" {
    const test_tle =
        \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
        \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
    ;

    var tle = try Tle.parse(test_tle, std.testing.allocator);
    defer tle.deinit();

    try std.testing.expectEqual(55909, tle.firstLine.satelliteNumber);
    try std.testing.expectEqual(0.006211566, tle.firstLine.bstarDrag);
    try std.testing.expectEqual(0.00023579, tle.firstLine.firstDerMeanMotion);

    try std.testing.expectEqual(43.9978, tle.secondLine.inclination);
    try std.testing.expectEqual(311.8012, tle.secondLine.rightAscension);
    try std.testing.expectEqual(0.0011446, tle.secondLine.eccentricity);
    try std.testing.expectEqual(15.05761711, tle.secondLine.mMotion);
}
