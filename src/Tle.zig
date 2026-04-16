//! Orbital element set parsed from TLE or OMM format.

const std = @import("std");
const DateTime = @import("Datetime.zig");

const Tle = @This();

// Identification
satelliteNumber: u32,
classification: u8,
intlDesignator: []u8,

// Epoch
epochYear: u16,
epochDay: f64,
epoch: f64, // J2000 seconds
epochJd: f64, // Julian Date (precomputed)

// Drag / perturbation
firstDerMeanMotion: f64,
bstarDrag: f64,
ephemType: u8,
elemNumber: u32,

// Orbital elements (degrees for angles, rev/day for mean motion)
inclination: f64,
rightAscension: f64,
eccentricity: f64,
perigee: f64,
mAnomaly: f64,
mMotion: f64,
revNum: u32,

allocator: std.mem.Allocator,

// ---------------------------------------------------------------------------
// TLE text parsing
// ---------------------------------------------------------------------------

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
    if (line1.len < 69 or line2.len < 69) return Error.BadTleLength;

    // International designator
    const intlYearField = trimField(line1, 9, 11);
    const intlLaunchField = trimField(line1, 11, 14);
    const intlPiece = std.mem.trim(u8, line1[14..17], " ");

    const intlDesignator = if (intlYearField.len > 0 and intlLaunchField.len > 0) blk: {
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
    errdefer allocator.free(intlDesignator);

    // B* drag
    const mantissa = try std.fmt.parseFloat(f64, trimField(line1, 53, 59));
    const exponent = try std.fmt.parseInt(i32, trimField(line1, 59, 61), 10);
    const bstarDrag = (mantissa * 1e-5) * std.math.pow(f64, 10.0, @as(f64, @floatFromInt(exponent)));

    // Epoch
    const epochYear = try std.fmt.parseInt(u16, trimField(line1, 18, 20), 10);
    const epochDay = try std.fmt.parseFloat(f64, trimField(line1, 20, 32));
    const epoch = tleEpochToJ2000(epochYear, epochDay);
    const epochJd = tleEpochToJd(epochYear, epochDay);

    // Eccentricity: TLE stores as integer, e.g. 0011446 → 0.0011446
    const eccentricity = try std.fmt.parseFloat(f64, trimField(line2, 26, 33)) / 1e7;

    return .{
        .satelliteNumber = try parseSatelliteNumber(trimField(line1, 2, 7)),
        .classification = line1[7],
        .intlDesignator = intlDesignator,
        .epochYear = epochYear,
        .epochDay = epochDay,
        .epoch = epoch,
        .epochJd = epochJd,
        .firstDerMeanMotion = try std.fmt.parseFloat(f64, trimField(line1, 33, 43)),
        .bstarDrag = bstarDrag,
        .ephemType = line1[62],
        .elemNumber = try std.fmt.parseInt(u32, trimField(line1, 64, 68), 10),
        .inclination = try std.fmt.parseFloat(f64, trimField(line2, 8, 16)),
        .rightAscension = try std.fmt.parseFloat(f64, trimField(line2, 17, 25)),
        .eccentricity = eccentricity,
        .perigee = try std.fmt.parseFloat(f64, trimField(line2, 34, 42)),
        .mAnomaly = try std.fmt.parseFloat(f64, trimField(line2, 43, 51)),
        .mMotion = try std.fmt.parseFloat(f64, trimField(line2, 52, 63)),
        .revNum = try std.fmt.parseInt(u32, trimField(line2, 63, 68), 10),
        .allocator = allocator,
    };
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

// ---------------------------------------------------------------------------
// OMM JSON parsing
// ---------------------------------------------------------------------------

pub fn parseOmm(json_text: []const u8, allocator: std.mem.Allocator) !Tle {
    const parsed = try std.json.parseFromSlice(OmmRecord, allocator, json_text, .{
        .ignore_unknown_fields = true,
    });
    defer parsed.deinit();
    return ommRecordToTle(parsed.value, allocator);
}

pub fn parseOmmArray(json_text: []const u8, allocator: std.mem.Allocator) ![]Tle {
    const parsed = try std.json.parseFromSlice([]OmmRecord, allocator, json_text, .{
        .ignore_unknown_fields = true,
    });
    defer parsed.deinit();

    const records = parsed.value;
    const tles = try allocator.alloc(Tle, records.len);
    var count: usize = 0;
    errdefer {
        for (tles[0..count]) |*t| {
            var tle = t.*;
            tle.deinit();
        }
        allocator.free(tles);
    }

    for (records, 0..) |rec, i| {
        tles[i] = try ommRecordToTle(rec, allocator);
        count += 1;
    }
    return tles;
}

const OmmRecord = struct {
    OBJECT_NAME: ?[]const u8 = null,
    OBJECT_ID: ?[]const u8 = null,
    EPOCH: []const u8,
    MEAN_MOTION: f64,
    ECCENTRICITY: f64,
    INCLINATION: f64,
    RA_OF_ASC_NODE: f64,
    ARG_OF_PERICENTER: f64,
    MEAN_ANOMALY: f64,
    EPHEMERIS_TYPE: ?u8 = 0,
    CLASSIFICATION_TYPE: ?[]const u8 = null,
    NORAD_CAT_ID: u32,
    ELEMENT_SET_NO: ?u32 = 0,
    REV_AT_EPOCH: ?u32 = 0,
    BSTAR: f64,
    MEAN_MOTION_DOT: ?f64 = 0,
    MEAN_MOTION_DDOT: ?f64 = 0,
};

fn ommRecordToTle(rec: OmmRecord, allocator: std.mem.Allocator) !Tle {
    // Parse ISO 8601 epoch: "2026-04-15T13:17:52.692576"
    const ep = try parseIso8601Epoch(rec.EPOCH);

    const intlDesignator = if (rec.OBJECT_ID) |id|
        try allocator.dupe(u8, id)
    else
        try allocator.dupe(u8, "");

    return .{
        .satelliteNumber = rec.NORAD_CAT_ID,
        .classification = if (rec.CLASSIFICATION_TYPE) |ct| (if (ct.len > 0) ct[0] else 'U') else 'U',
        .intlDesignator = intlDesignator,
        .epochYear = ep.year,
        .epochDay = ep.doy,
        .epoch = ep.j2000,
        .epochJd = ep.jd,
        .firstDerMeanMotion = rec.MEAN_MOTION_DOT orelse 0,
        .bstarDrag = rec.BSTAR,
        .ephemType = rec.EPHEMERIS_TYPE orelse 0,
        .elemNumber = rec.ELEMENT_SET_NO orelse 0,
        .inclination = rec.INCLINATION,
        .rightAscension = rec.RA_OF_ASC_NODE,
        .eccentricity = rec.ECCENTRICITY,
        .perigee = rec.ARG_OF_PERICENTER,
        .mAnomaly = rec.MEAN_ANOMALY,
        .mMotion = rec.MEAN_MOTION,
        .revNum = rec.REV_AT_EPOCH orelse 0,
        .allocator = allocator,
    };
}

const ParsedEpoch = struct { year: u16, doy: f64, jd: f64, j2000: f64 };

fn parseIso8601Epoch(epoch: []const u8) !ParsedEpoch {
    // "2026-04-15T13:17:52.692576"
    if (epoch.len < 19) return Error.BadTleLength;

    const year = try std.fmt.parseInt(u16, epoch[0..4], 10);
    const month = try std.fmt.parseInt(u8, epoch[5..7], 10);
    const day = try std.fmt.parseInt(u8, epoch[8..10], 10);
    const hour = try std.fmt.parseInt(u8, epoch[11..13], 10);
    const min = try std.fmt.parseInt(u8, epoch[14..16], 10);

    // Seconds may have fractional part; strip trailing 'Z' if present
    const secEnd = if (epoch[epoch.len - 1] == 'Z') epoch.len - 1 else epoch.len;
    const sec = try std.fmt.parseFloat(f64, epoch[17..secEnd]);

    // Day of year (fractional)
    const doy = dayOfYear(year, month, day) +
        (@as(f64, @floatFromInt(hour)) +
            (@as(f64, @floatFromInt(min)) + sec / 60.0) / 60.0) / 24.0;

    // Two-digit year for epochYear field (OMM uses 4-digit, store mod 100)
    const epochYear: u16 = year % 100;
    const jd = DateTime.yearDoyToJulianDate(year, doy);
    const j2000 = tleEpochToJ2000(epochYear, doy);

    return .{ .year = epochYear, .doy = doy, .jd = jd, .j2000 = j2000 };
}

fn dayOfYear(year: u16, month: u8, day: u8) f64 {
    const cumDays = [_]u16{ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 };
    var doy: u16 = cumDays[month - 1] + day;
    if (month > 2 and ((year % 4 == 0 and year % 100 != 0) or (year % 400 == 0))) {
        doy += 1;
    }
    return @floatFromInt(doy);
}

// ---------------------------------------------------------------------------
// Common
// ---------------------------------------------------------------------------

pub fn deinit(self: *Tle) void {
    self.allocator.free(self.intlDesignator);
}

pub fn output(self: Tle) void {
    std.debug.print("satellite_number: {d}\n", .{self.satelliteNumber});
    std.debug.print("classification: {c}\n", .{self.classification});
    std.debug.print("intl_designator: {s}\n", .{self.intlDesignator});
    std.debug.print("epoch_year: {d}\n", .{self.epochYear});
    std.debug.print("epoch_day: {d}\n", .{self.epochDay});
    std.debug.print("epoch_jd: {d}\n", .{self.epochJd});
    std.debug.print("first_der_m_motion: {d}\n", .{self.firstDerMeanMotion});
    std.debug.print("bstar_drag: {d}\n", .{self.bstarDrag});
    std.debug.print("ephem_type: {d}\n", .{self.ephemType});
    std.debug.print("elem_number: {d}\n", .{self.elemNumber});
    std.debug.print("inclination: {d}\n", .{self.inclination});
    std.debug.print("right_ascension: {d}\n", .{self.rightAscension});
    std.debug.print("eccentricity: {d}\n", .{self.eccentricity});
    std.debug.print("perigee: {d}\n", .{self.perigee});
    std.debug.print("m_anomaly: {d}\n", .{self.mAnomaly});
    std.debug.print("m_motion: {d}\n", .{self.mMotion});
    std.debug.print("rev_num: {d}\n", .{self.revNum});
}

pub const Error = error{BadTleLength};

fn trimField(line: []const u8, start: usize, end: usize) []const u8 {
    return std.mem.trim(u8, line[start..end], " ");
}

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

fn tleEpochToJ2000(epochYear: u16, epochDay: f64) f64 {
    const y = 2000 + epochYear;
    const md = DateTime.doyToMonthDay(y, epochDay);
    return DateTime.initDate(y, md.month, md.day).convertToJ2000();
}

fn tleEpochToJd(epochYear: u16, epochDay: f64) f64 {
    const fullYear: u16 = if (epochYear < 57)
        2000 + epochYear
    else
        1900 + epochYear;
    return DateTime.yearDoyToJulianDate(fullYear, epochDay);
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "parseLines and MultiIterator" {
    const line1 = "1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998";
    const line2 = "2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371";

    // parseLines
    {
        var tle = try Tle.parseLines(line1, line2, std.testing.allocator);
        defer tle.deinit();
        try std.testing.expectEqual(@as(u32, 55909), tle.satelliteNumber);
        try std.testing.expectApproxEqAbs(@as(f64, 43.9978), tle.inclination, 1e-6);
    }

    // parse handles CRLF, blank lines, leading/trailing whitespace
    {
        const crlf = "  " ++ line1 ++ "  \r\n\r\n  " ++ line2 ++ "  ";
        var tle = try Tle.parse(crlf, std.testing.allocator);
        defer tle.deinit();
        try std.testing.expectEqual(@as(u32, 55909), tle.satelliteNumber);
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

test "parseOmm" {
    const json =
        \\{"OBJECT_NAME":"ISS (ZARYA)","OBJECT_ID":"1998-067A","EPOCH":"2026-04-15T13:17:52.692576","MEAN_MOTION":15.48924547,"ECCENTRICITY":0.00065833,"INCLINATION":51.6327,"RA_OF_ASC_NODE":250.0746,"ARG_OF_PERICENTER":315.773,"MEAN_ANOMALY":44.2732,"EPHEMERIS_TYPE":0,"CLASSIFICATION_TYPE":"U","NORAD_CAT_ID":25544,"ELEMENT_SET_NO":999,"REV_AT_EPOCH":56204,"BSTAR":0.00010353824,"MEAN_MOTION_DOT":0.00005244,"MEAN_MOTION_DDOT":0}
    ;

    var tle = try Tle.parseOmm(json, std.testing.allocator);
    defer tle.deinit();

    try std.testing.expectEqual(@as(u32, 25544), tle.satelliteNumber);
    try std.testing.expectApproxEqAbs(@as(f64, 51.6327), tle.inclination, 1e-4);
    try std.testing.expectApproxEqAbs(@as(f64, 0.00065833), tle.eccentricity, 1e-8);
    try std.testing.expectApproxEqAbs(@as(f64, 15.48924547), tle.mMotion, 1e-8);
    try std.testing.expectApproxEqAbs(@as(f64, 0.00010353824), tle.bstarDrag, 1e-11);
    try std.testing.expect(std.mem.eql(u8, tle.intlDesignator, "1998-067A"));
}

test "parseOmmArray" {
    const json =
        \\[{"OBJECT_NAME":"ISS","OBJECT_ID":"1998-067A","EPOCH":"2026-04-15T13:17:52.692576","MEAN_MOTION":15.489,"ECCENTRICITY":0.0006,"INCLINATION":51.63,"RA_OF_ASC_NODE":250.07,"ARG_OF_PERICENTER":315.77,"MEAN_ANOMALY":44.27,"NORAD_CAT_ID":25544,"BSTAR":0.0001},{"OBJECT_NAME":"STARLINK","OBJECT_ID":"2023-035B","EPOCH":"2026-04-15T12:00:00.000000","MEAN_MOTION":15.05,"ECCENTRICITY":0.001,"INCLINATION":43.99,"RA_OF_ASC_NODE":311.80,"ARG_OF_PERICENTER":278.62,"MEAN_ANOMALY":81.33,"NORAD_CAT_ID":55909,"BSTAR":0.00016}]
    ;

    const tles = try Tle.parseOmmArray(json, std.testing.allocator);
    defer {
        for (tles) |*t| {
            var tle = t.*;
            tle.deinit();
        }
        std.testing.allocator.free(tles);
    }

    try std.testing.expectEqual(@as(usize, 2), tles.len);
    try std.testing.expectEqual(@as(u32, 25544), tles[0].satelliteNumber);
    try std.testing.expectEqual(@as(u32, 55909), tles[1].satelliteNumber);
}
