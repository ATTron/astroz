//! SPICE file structure
const std = @import("std");
const log = std.log;
const ArrayList = std.ArrayList;

pub const Spice = @This();

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

pub const LSK = struct {
    delta_t_a: f64,
    k: f64,
    eb: f64,
    m: [2]f64,
    leap_seconds: ArrayList(LeapSecondEntry),
    allocator: std.mem.Allocator,

    pub fn deinit(self: *LSK) void {
        self.leap_seconds.deinit(self.allocator);
    }
};

lsk: ?LSK,

pub fn new(lsk: LSK) Spice {
    return .{
        .lsk = lsk,
    };
}

pub fn parseLSK(allocator: std.mem.Allocator, content: []const u8) LSKError!LSK {
    var lsk = LSK{
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

        if (std.mem.startsWith(u8, trimmed, "DELTET/DELTA_T_A")) {
            lsk.delta_t_a = parseFloatFromAssignment(trimmed) catch return LSKError.ParseError;
        } else if (std.mem.startsWith(u8, trimmed, "DELTET/K")) {
            lsk.k = parseFloatFromAssignment(trimmed) catch return LSKError.ParseError;
        } else if (std.mem.startsWith(u8, trimmed, "DELTET/EB")) {
            lsk.eb = parseFloatFromAssignment(trimmed) catch return LSKError.ParseError;
        } else if (std.mem.startsWith(u8, trimmed, "DELTET/M")) {
            parseMFromAssignment(trimmed, &lsk.m) catch return LSKError.ParseError;
        } else if (std.mem.startsWith(u8, trimmed, "DELTET/DELTA_AT")) {
            var full_line = std.ArrayList(u8){};
            defer full_line.deinit(allocator);

            try full_line.appendSlice(allocator, trimmed);

            var assignment_processed = false;
            while (lines.next()) |continuation_line| {
                const continuation_trimmed = std.mem.trim(u8, continuation_line, " \t\r\n");
                if (continuation_trimmed.len == 0) continue;

                if (continuation_line.len > 0 and (continuation_line[0] == ' ' or continuation_line[0] == '\t')) {
                    try full_line.appendSlice(allocator, " ");
                    try full_line.appendSlice(allocator, continuation_trimmed);
                } else {
                    parseDeltaAtFromAssignment(full_line.items, &lsk.leap_seconds, allocator) catch return LSKError.ParseError;
                    assignment_processed = true;

                    if (std.mem.eql(u8, continuation_trimmed, "\\\\begintext")) {
                        in_data_section = false;
                        break;
                    }
                    break;
                }
            }

            if (!assignment_processed and full_line.items.len > 0) {
                parseDeltaAtFromAssignment(full_line.items, &lsk.leap_seconds, allocator) catch return LSKError.ParseError;
            }
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

test "LSK parser basic functionality" {
    const testing = std.testing;
    const allocator = testing.allocator;

    const sample_lsk =
        \\KPL/LSK
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

    var lsk = try parseLSK(allocator, sample_lsk);
    defer lsk.deinit();

    try testing.expectEqual(@as(f64, 32.184), lsk.delta_t_a);
    try testing.expectApproxEqRel(@as(f64, 1.657e-3), lsk.k, 1e-10);
    try testing.expectApproxEqRel(@as(f64, 1.671e-2), lsk.eb, 1e-10);
    try testing.expectApproxEqRel(@as(f64, 6.239996), lsk.m[0], 1e-10);
    try testing.expectApproxEqRel(@as(f64, 1.99096871e-7), lsk.m[1], 1e-15);
    try testing.expectEqual(@as(usize, 3), lsk.leap_seconds.items.len);
}

test "LSK parser error handling" {
    const testing = std.testing;
    const allocator = testing.allocator;

    const invalid_lsk = "invalid content";
    var lsk = try parseLSK(allocator, invalid_lsk);
    defer lsk.deinit();

    try testing.expectEqual(@as(f64, 0), lsk.delta_t_a);
    try testing.expectEqual(@as(usize, 0), lsk.leap_seconds.items.len);
}

test "LSK parse float from assignment" {
    const testing = std.testing;

    try testing.expectEqual(@as(f64, 32.184), try parseFloatFromAssignment("DELTET/DELTA_T_A = 32.184"));
    try testing.expectApproxEqRel(@as(f64, 1.657e-3), try parseFloatFromAssignment("DELTET/K = 1.657D-3"), 1e-10);
}

test "LSK parse M values" {
    const testing = std.testing;
    var m: [2]f64 = undefined;

    try parseMFromAssignment("DELTET/M = ( 6.239996D0 1.99096871D-7 )", &m);
    try testing.expectApproxEqRel(@as(f64, 6.239996), m[0], 1e-10);
    try testing.expectApproxEqRel(@as(f64, 1.99096871e-7), m[1], 1e-15);
}
