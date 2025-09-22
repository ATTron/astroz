//! LSK file structure
const std = @import("std");
const log = std.log;
const ArrayList = std.ArrayList;

pub const Lsk = @This();

const FieldType = enum {
    deltaTA,
    k,
    eb,
    m,
    unknown,
};

pub const LSKError = error{
    InvalidFormat,
    MissingData,
    ParseError,
    OutOfMemory,
};

pub const LeapSecondEntry = struct {
    delta_at: f64,
    jd_utc: f64,
};

delta_t_a: f64,
k: f64,
eb: f64,
m: [2]f64,
leap_seconds: ArrayList(LeapSecondEntry),
allocator: std.mem.Allocator,

pub fn deinit(self: *Lsk) void {
    self.leap_seconds.deinit(self.allocator);
}

pub fn parseLsk(allocator: std.mem.Allocator, content: []const u8) LSKError!Lsk {
    var lsk = Lsk{
        .delta_t_a = 0,
        .k = 0,
        .eb = 0,
        .m = [2]f64{ 0, 0 },
        .leap_seconds = ArrayList(LeapSecondEntry){},
        .allocator = allocator,
    };

    var in_data_section = false;
    var lines = std.mem.splitAny(u8, content, "\n");

    while (lines.next()) |line| {
        const trimmed = std.mem.trim(u8, line, " \t\r\n");

        if (std.mem.eql(u8, trimmed, "\\\\begindata")) {
            in_data_section = true;
            continue;
        }

        if (std.mem.eql(u8, trimmed, "\\\\begintext")) {
            in_data_section = false;
            continue;
        }

        if (!in_data_section or trimmed.len == 0) continue;

        if (std.mem.startsWith(u8, trimmed, "DELTET/DELTA_AT")) {
            try handleDeltaAt(trimmed, &lines, &lsk, allocator, &in_data_section);
        } else {
            try handleSimpleField(trimmed, &lsk);
        }
    }

    return lsk;
}

fn convertSpiceNotation(input: []const u8, buffer: []u8) []const u8 {
    if (std.mem.indexOf(u8, input, "D")) |d_pos| {
        const len = input.len;
        if (len < buffer.len) {
            @memcpy(buffer[0..len], input);
            buffer[d_pos] = 'e';
            return buffer[0..len];
        }
    }
    return input;
}

fn parseFloatFromAssignment(line: []const u8) !f64 {
    if (std.mem.indexOf(u8, line, "=")) |eq_pos| {
        const value_part = std.mem.trim(u8, line[eq_pos + 1 ..], " \t");

        var buffer: [256]u8 = undefined;
        const converted_value = convertSpiceNotation(value_part, &buffer);

        return std.fmt.parseFloat(f64, converted_value);
    }
    return error.ParseError;
}

fn parseMFromAssignment(line: []const u8, m: *[2]f64) !void {
    if (std.mem.indexOf(u8, line, "=")) |eq_pos| {
        const value_part = std.mem.trim(u8, line[eq_pos + 1 ..], " \t()");
        var values = std.mem.splitAny(u8, value_part, " ");

        var count: u32 = 0;
        while (values.next()) |val| {
            const trimmed = std.mem.trim(u8, val, " \t");
            if (trimmed.len == 0) continue;

            var buffer: [256]u8 = undefined;
            const final_val = convertSpiceNotation(trimmed, &buffer);

            if (count < 2) {
                m[count] = try std.fmt.parseFloat(f64, final_val);
                count += 1;
            }
        }
    }
}

fn parseDeltaAtFromAssignment(line: []const u8, leap_seconds: *ArrayList(LeapSecondEntry), allocator: std.mem.Allocator) !void {
    if (std.mem.indexOf(u8, line, "=")) |eq_pos| {
        const value_part = std.mem.trim(u8, line[eq_pos + 1 ..], " \t()");

        var tokens = std.mem.tokenizeAny(u8, value_part, " \t,()");

        while (tokens.next()) |token| {
            if (std.mem.startsWith(u8, token, "@")) continue;

            if (std.fmt.parseFloat(f64, token)) |delta_at| {
                try leap_seconds.append(allocator, LeapSecondEntry{
                    .delta_at = delta_at,
                    .jd_utc = 0.0, // placeholder for now
                });
            } else |_| {
                continue;
            }
        }
    }
}

fn fromPrefix(trimmed_line: []const u8) FieldType {
    if (std.mem.startsWith(u8, trimmed_line, "DELTET/DELTA_T_A")) return .deltaTA;
    if (std.mem.startsWith(u8, trimmed_line, "DELTET/K")) return .k;
    if (std.mem.startsWith(u8, trimmed_line, "DELTET/EB")) return .eb;
    if (std.mem.startsWith(u8, trimmed_line, "DELTET/M")) return .m;
    return .unknown;
}

fn handleSimpleField(trimmed: []const u8, lsk: *Lsk) !void {
    switch (fromPrefix(trimmed)) {
        .deltaTA => lsk.delta_t_a = parseFloatFromAssignment(trimmed) catch return LSKError.ParseError,
        .k => lsk.k = parseFloatFromAssignment(trimmed) catch return LSKError.ParseError,
        .eb => lsk.eb = parseFloatFromAssignment(trimmed) catch return LSKError.ParseError,
        .m => parseMFromAssignment(trimmed, &lsk.m) catch return LSKError.ParseError,
        .unknown => {},
    }
}

fn handleDeltaAt(trimmed: []const u8, lines: *std.mem.SplitIterator(u8, .any), lsk: *Lsk, allocator: std.mem.Allocator, inDataSection: *bool) !void {
    var fullLine: std.ArrayList(u8) = .{};
    defer fullLine.deinit(allocator);

    try fullLine.appendSlice(allocator, trimmed);
    var assignmentProcessed = false;

    while (lines.next()) |continuationLine| {
        const continuationTrimmed = std.mem.trim(u8, continuationLine, " \t\r\n");
        if (continuationTrimmed.len == 0) continue;
        if (continuationLine.len > 0 and (continuationLine[0] == ' ' or continuationLine[0] == '\t')) {
            try fullLine.appendSlice(allocator, " ");
            try fullLine.appendSlice(allocator, continuationTrimmed);
        } else {
            parseDeltaAtFromAssignment(fullLine.items, &lsk.leap_seconds, allocator) catch return LSKError.ParseError;
            assignmentProcessed = true;
            if (std.mem.eql(u8, continuationTrimmed, "\\\\begintext")) {
                inDataSection.* = false;
                break;
            }
            break;
        }
    }
    if (!assignmentProcessed and fullLine.items.len > 0) {
        parseDeltaAtFromAssignment(fullLine.items, &lsk.leap_seconds, allocator) catch return LSKError.ParseError;
    }
}

test "Lsk parser basic functionality" {
    const testing = std.testing;
    const allocator = testing.allocator;

    const sample_lsk =
        \\KPL/Lsk
        \\
        \\LEAPSECONDS KERNEL FILE
        \\
        \\\\begindata
        \\
        \\DELTET/DELTA_T_A       =   32.184
        \\DELTET/K               =    1.657D-3
        \\DELTET/EB              =    1.671D-2
        \\DELTET/M               = (  6.239996D0   1.99096871D-7 )
        \\
        \\DELTET/DELTA_AT        = ( 10,   @1972-JAN-1
        \\                           11,   @1972-JUL-1
        \\                           12,   @1973-JAN-1 )
        \\
        \\\\begintext
    ;

    var lsk = try parseLsk(allocator, sample_lsk);
    defer lsk.deinit();

    try testing.expectEqual(@as(f64, 32.184), lsk.delta_t_a);
    try testing.expectApproxEqRel(@as(f64, 1.657e-3), lsk.k, 1e-10);
    try testing.expectApproxEqRel(@as(f64, 1.671e-2), lsk.eb, 1e-10);
    try testing.expectApproxEqRel(@as(f64, 6.239996), lsk.m[0], 1e-10);
    try testing.expectApproxEqRel(@as(f64, 1.99096871e-7), lsk.m[1], 1e-15);
    try testing.expectEqual(@as(usize, 3), lsk.leap_seconds.items.len);
}

test "Lsk parser error handling" {
    const testing = std.testing;
    const allocator = testing.allocator;

    const invalid_lsk = "invalid content";
    var lsk = try parseLsk(allocator, invalid_lsk);
    defer lsk.deinit();

    try testing.expectEqual(@as(f64, 0), lsk.delta_t_a);
    try testing.expectEqual(@as(usize, 0), lsk.leap_seconds.items.len);
}

test "Lsk parse float from assignment" {
    const testing = std.testing;

    try testing.expectEqual(@as(f64, 32.184), try parseFloatFromAssignment("DELTET/DELTA_T_A = 32.184"));
    try testing.expectApproxEqRel(@as(f64, 1.657e-3), try parseFloatFromAssignment("DELTET/K = 1.657D-3"), 1e-10);
}

test "Lsk parse M values" {
    const testing = std.testing;
    var m: [2]f64 = undefined;

    try parseMFromAssignment("DELTET/M = ( 6.239996D0 1.99096871D-7 )", &m);
    try testing.expectApproxEqRel(@as(f64, 6.239996), m[0], 1e-10);
    try testing.expectApproxEqRel(@as(f64, 1.99096871e-7), m[1], 1e-15);
}
