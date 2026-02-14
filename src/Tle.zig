//! The proper TLE struct

const std = @import("std");
const DateTime = @import("Datetime.zig");

const Tle = @This();

firstLine: FirstLine,
secondLine: SecondLine,
allocator: std.mem.Allocator,

pub fn parse(tle: []const u8, allocator: std.mem.Allocator) !Tle {
    var lines = std.mem.splitAny(u8, tle, "\n\r");
    var line1: ?[]const u8 = null;

    while (lines.next()) |raw_line| {
        const line = std.mem.trim(u8, raw_line, " \t");
        if (line.len < 69) continue;
        if (line1 == null) {
            line1 = line;
        } else {
            return parseLines(line1.?, line, allocator);
        }
    }
    return Error.BadTleLength;
}

pub fn parseLines(line1: []const u8, line2: []const u8, allocator: std.mem.Allocator) !Tle {
    const firstLine = try FirstLine.init(line1, allocator);
    errdefer allocator.free(firstLine.intlDesignator);
    const secondLine = try SecondLine.init(line2, allocator);
    return .{ .firstLine = firstLine, .secondLine = secondLine, .allocator = allocator };
}

pub const MultiIterator = struct {
    lines: std.mem.SplitIterator(u8, .any),

    pub const LinePair = struct { line1: []const u8, line2: []const u8 };

    pub fn init(text: []const u8) MultiIterator {
        return .{ .lines = std.mem.splitAny(u8, text, "\n\r") };
    }

    pub fn next(self: *MultiIterator) ?LinePair {
        var candidateLine: ?[]const u8 = null;

        while (self.lines.next()) |raw_line| {
            const line = std.mem.trim(u8, raw_line, " \t");
            if (line.len < 69) continue;

            switch (line[0]) {
                '1' => candidateLine = line,
                '2' => {
                    if (candidateLine) |l1| {
                        candidateLine = null;
                        return .{ .line1 = l1, .line2 = line };
                    }
                },
                else => candidateLine = null,
            }
        }
        return null;
    }
};

pub fn deinit(self: *Tle) void {
    self.allocator.free(self.firstLine.intlDesignator);
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

fn trimField(line: []const u8, start: usize, end: usize) []const u8 {
    return std.mem.trim(u8, line[start..end], " ");
}

/// Parse a satellite catalog number, supporting both numeric and Alpha-5 encoding.
/// Alpha-5 uses a leading letter (A=10..Z=35) followed by 4 digits for numbers > 99999.
fn parseSatelliteNumber(field: []const u8) !u32 {
    if (field.len == 0) return error.InvalidCharacter;
    const first = field[0];
    if (first >= 'A' and first <= 'Z') {
        const prefix: u32 = @as(u32, first - 'A') + 10;
        const rest = try std.fmt.parseInt(u32, field[1..], 10);
        return prefix * 10000 + rest;
    }
    return std.fmt.parseInt(u32, field, 10);
}

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
    secondDerMeanMotion: [8]u8,
    bstarDrag: f64,
    ephemType: u8,
    elemNumber: u32,
    checksum: u8,

    pub fn init(line: []const u8, allocator: std.mem.Allocator) !FirstLine {
        if (line.len < 69) {
            return Error.BadTleLength;
        }

        const intlYearField = trimField(line, 9, 11);
        const intlLaunchField = trimField(line, 11, 14);
        const intlPiece = std.mem.trim(u8, line[14..17], " ");

        const string = if (intlYearField.len > 0 and intlLaunchField.len > 0) blk: {
            const intlYear = try std.fmt.parseInt(u16, intlYearField, 10);
            const intlLaunchNum = try std.fmt.parseInt(u32, intlLaunchField, 10);
            break :blk try std.fmt.allocPrint(
                allocator,
                "{d:0>2}-{d:0>3}{s}",
                .{ intlYear, intlLaunchNum, intlPiece },
            );
        } else blk: {
            break :blk try allocator.dupe(u8, "");
        };
        errdefer allocator.free(string);

        const mantissa = try std.fmt.parseFloat(f64, trimField(line, 53, 59));
        const exponent = try std.fmt.parseInt(i32, trimField(line, 59, 61), 10);
        const bstarDrag = (mantissa * 1e-5) * std.math.pow(f64, 10.0, @as(f64, @floatFromInt(exponent)));

        const epochYear = try std.fmt.parseInt(u16, trimField(line, 18, 20), 10);
        const epochDay = try std.fmt.parseFloat(f64, trimField(line, 20, 32));
        const epoch = tleEpochToEpoch(epochYear, epochDay);

        return .{
            .lineNumber = line[0],
            .satelliteNumber = try parseSatelliteNumber(trimField(line, 2, 7)),
            .classification = line[7],
            .intlDesignator = string,
            .epochYear = epochYear,
            .epochDay = epochDay,
            .epoch = epoch,
            .firstDerMeanMotion = try std.fmt.parseFloat(f64, trimField(line, 33, 43)),
            .secondDerMeanMotion = line[44..52].*,
            .bstarDrag = bstarDrag,
            .ephemType = line[62],
            .elemNumber = try std.fmt.parseInt(u32, trimField(line, 64, 68), 10),
            .checksum = line[68],
        };
    }

    fn tleEpochToEpoch(year: u16, doy: f64) f64 {
        const y = 2000 + year;
        const md = DateTime.doyToMonthDay(y, doy);
        return DateTime.initDate(y, md.month, md.day).convertToJ2000();
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

    pub fn init(line: []const u8, _: std.mem.Allocator) !SecondLine {
        if (line.len < 69) {
            return Error.BadTleLength;
        }

        return .{
            .lineNumber = line[0],
            .satelliteNumber = try parseSatelliteNumber(trimField(line, 2, 7)),
            .inclination = try std.fmt.parseFloat(f64, trimField(line, 8, 16)),
            .rightAscension = try std.fmt.parseFloat(f64, trimField(line, 17, 25)),
            .eccentricity = try std.fmt.parseFloat(f64, trimField(line, 26, 33)) / 1e7,
            .perigee = try std.fmt.parseFloat(f64, trimField(line, 34, 42)),
            .mAnomaly = try std.fmt.parseFloat(f64, trimField(line, 43, 51)),
            .mMotion = try std.fmt.parseFloat(f64, trimField(line, 52, 63)),
            .revNum = try std.fmt.parseInt(u32, trimField(line, 63, 68), 10),
            .checksum = line[68],
        };
    }
};

test "parseLines and MultiIterator" {
    const line1 = "1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998";
    const line2 = "2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371";

    // parseLines
    {
        var tle = try Tle.parseLines(line1, line2, std.testing.allocator);
        defer tle.deinit();
        try std.testing.expectEqual(@as(u32, 55909), tle.firstLine.satelliteNumber);
        try std.testing.expectApproxEqAbs(@as(f64, 43.9978), tle.secondLine.inclination, 1e-6);
    }

    // parse handles CRLF, blank lines, leading/trailing whitespace
    {
        const crlf = "  " ++ line1 ++ "  \r\n\r\n  " ++ line2 ++ "  ";
        var tle = try Tle.parse(crlf, std.testing.allocator);
        defer tle.deinit();
        try std.testing.expectEqual(@as(u32, 55909), tle.firstLine.satelliteNumber);
        try std.testing.expectEqual(@as(u8, '8'), tle.firstLine.checksum);
    }

    // parse with missing second line returns error
    try std.testing.expectError(Error.BadTleLength, Tle.parse(line1, std.testing.allocator));

    // MultiIterator: 3-line format with names, orphaned line 1, empty/garbage
    {
        const text =
            "ISS (ZARYA)\n" ++
            "1 25544U 98067A   24187.50000000  .00016717  00000+0  10270-3 0  9993\n" ++
            "2 25544  51.6400 208.5000 0007417  35.0000 325.0000 15.49000000400000\n" ++
            // orphaned line 1 (no matching line 2 before next line 1)
            "1 99999U 00000A   24001.00000000  .00000000  00000+0  00000-0 0  0000\n" ++
            "STARLINK-1234\n" ++
            line1 ++ "\n" ++ line2 ++ "\n";

        var iter = MultiIterator.init(text);

        const pair1 = iter.next().?;
        try std.testing.expect(std.mem.startsWith(u8, pair1.line1[2..7], "25544"));

        const pair2 = iter.next().?;
        try std.testing.expect(std.mem.startsWith(u8, pair2.line1[2..7], "55909"));

        try std.testing.expect(iter.next() == null);
    }

    // Empty/garbage input
    {
        var iter = MultiIterator.init("hello\ngarbage\n");
        try std.testing.expect(iter.next() == null);
    }
}
