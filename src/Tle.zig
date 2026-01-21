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
    epochDay: f64,
    epoch: f64,
    firstDerMeanMotion: f64,
    secondDerMeanMotion: [7]u8,
    bstarDrag: f64,
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
        errdefer allocator.free(string);

        const mantissa = try std.fmt.parseFloat(f64, line[54..59]);
        const exponent = try std.fmt.parseInt(i32, line[59..61], 10);
        const bstarDrag = (mantissa * 1e-5) * std.math.pow(f64, 10.0, @as(f64, @floatFromInt(exponent)));

        const epochYear = try std.fmt.parseInt(u16, line[18..20], 10);
        const epochDay = try std.fmt.parseFloat(f64, line[20..32]);
        const epoch = tleEpochToEpoch(epochYear, epochDay);

        return .{
            .lineNumber = line[0],
            .satelliteNumber = try std.fmt.parseInt(u32, line[2..7], 10),
            .classification = line[7],
            .intlDesignator = string,
            .epochYear = epochYear,
            .epochDay = epochDay,
            .epoch = epoch,
            .firstDerMeanMotion = try std.fmt.parseFloat(f64, line[34..43]),
            .secondDerMeanMotion = line[45..52].*,
            .bstarDrag = bstarDrag,
            .ephemType = line[62],
            .elemNumber = try std.fmt.parseInt(u32, line[65..68], 10),
            .checksum = line[68],
        };
    }

    fn tleEpochToEpoch(year: u16, doy: f64) f64 {
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
    inclination: f64,
    rightAscension: f64,
    eccentricity: f64,
    perigee: f64,
    mAnomaly: f64,
    mMotion: f64,
    revNum: u32,
    checksum: u8,
    cleanLine: []const u8,

    pub fn init(line: []const u8, allocator: std.mem.Allocator) !SecondLine {
        if (line.len < 69) {
            return Error.BadTleLength;
        }

        // store a copy for potential later use
        const lineCopy = try allocator.dupe(u8, line);

        const eccStr = std.mem.trim(u8, line[26..33], " ");
        const eccValue = try std.fmt.parseFloat(f64, eccStr);
        const eccentricity: f64 = eccValue / 1e7;

        return .{
            .lineNumber = line[0],
            .satelliteNumber = try std.fmt.parseInt(u32, std.mem.trim(u8, line[2..7], " "), 10),
            .inclination = try std.fmt.parseFloat(f64, std.mem.trim(u8, line[8..16], " ")),
            .rightAscension = try std.fmt.parseFloat(f64, std.mem.trim(u8, line[17..25], " ")),
            .eccentricity = eccentricity,
            .perigee = try std.fmt.parseFloat(f64, std.mem.trim(u8, line[34..42], " ")),
            .mAnomaly = try std.fmt.parseFloat(f64, std.mem.trim(u8, line[43..51], " ")),
            .mMotion = try std.fmt.parseFloat(f64, std.mem.trim(u8, line[52..63], " ")),
            .revNum = try std.fmt.parseInt(u32, std.mem.trim(u8, line[63..68], " "), 10),
            .checksum = line[68],
            .cleanLine = lineCopy,
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
    try std.testing.expectApproxEqAbs(@as(f64, 0.0016099), tle.firstLine.bstarDrag, 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 0.00023579), tle.firstLine.firstDerMeanMotion, 1e-10);

    try std.testing.expectApproxEqAbs(@as(f64, 43.9978), tle.secondLine.inclination, 1e-6);
    try std.testing.expectApproxEqAbs(@as(f64, 311.8012), tle.secondLine.rightAscension, 1e-6);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0011446), tle.secondLine.eccentricity, 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 15.05761711), tle.secondLine.mMotion, 1e-10);
}
