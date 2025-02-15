//! The Vita49 packet type which will be used by the Parser

const std = @import("std");

const Vita49 = @This();

header: Header,
stream_id: ?u32,
class_id: ?ClassId,
i_timestamp: ?u32,
f_timestamp: ?u64,
payload: []const u8,
trailer: ?Trailer,
allocator: std.mem.Allocator,
end: usize,

// private
raw_data: []const u8,

pub fn init(pl: []const u8, allocator: std.mem.Allocator, config: ?[]const u8) !Vita49 {
    var stream = try allocator.dupe(u8, pl); // dupe the data so we dont ever lose it
    if (config != null) {
        std.log.debug("Config found for Vita49 but this hasn't been implemented yet", .{});
    }
    if (stream.len < 4) {
        return Error.InsufficientData;
    }
    const header = try Header.init(stream[0..4].*);
    var stream_id: ?u32 = undefined;
    var class_id: ?ClassId = undefined;
    var payload: []const u8 = undefined;
    var trailer: ?Trailer = undefined;
    var i_start: usize = 4;
    var f_start: usize = 4;
    var i_timestamp: ?u32 = undefined;
    var f_timestamp: ?u64 = undefined;
    var payload_range = try getPayloadRange(header, true);

    switch (header.packet_type) {
        .signal_w_stream_id,
        .ext_data_w_stream_id,
        .ext_cmd_packet,
        .cmd_packet,
        .ctx_packet,
        .ext_ctx_packet,
        => {
            stream_id = std.mem.readInt(u32, stream[4..8], .little);
            i_start += 4;
            f_start += 4;
        },
        .signal_wo_stream_id,
        .ext_data_wo_stream_id,
        => {
            payload_range = try getPayloadRange(header, false);
        },
    }
    if (header.class_id) {
        class_id = ClassId.init(stream[8..16].*);
        i_start += 8;
        f_start += 8;
    } else {
        class_id = null;
    }
    if (header.trailer) {
        trailer = Trailer.init(stream[payload_range.end..]);
        payload = stream[payload_range.start..payload_range.end];
    } else {
        payload = stream[payload_range.start..payload_range.end];
    }
    if (header.tsi != Tsi.none) {
        var tmp_array: [4]u8 = undefined;
        @memcpy(&tmp_array, stream[i_start .. i_start + 4]);
        i_timestamp = std.mem.readInt(u32, &tmp_array, .little);
        f_start += 4;
    } else {
        i_timestamp = null;
    }
    if (header.tsf != Tsf.none) {
        var tmp_array: [8]u8 = undefined;
        @memcpy(&tmp_array, stream[f_start .. f_start + 8]);
        f_timestamp = std.mem.readInt(u64, &tmp_array, .little);
    } else {
        f_timestamp = null;
    }
    const end_of_data = if (header.trailer) payload_range.end + 4 else payload_range.end;

    _ = allocator.resize(stream, end_of_data);
    return .{
        .header = header,
        .stream_id = stream_id,
        .class_id = class_id,
        .i_timestamp = i_timestamp,
        .f_timestamp = f_timestamp,
        .payload = payload,
        .trailer = trailer,
        .allocator = allocator,
        .end = payload_range.end,
        .raw_data = stream,
    };
}

pub fn deinit(self: *Vita49) void {
    self.allocator.free(self.raw_data);
}

fn getPayloadRange(header: Header, has_stream_id: bool) !struct { start: usize, end: usize } {
    var start: usize = 4;
    var end: usize = @as(usize, (header.packet_size * 4)) - 1;
    if (has_stream_id) {
        start += 4;
    }
    if (header.class_id) {
        start += 8;
    }
    if (header.tsi != Tsi.none) {
        start += 4;
    }
    if (header.tsf != Tsf.none) {
        start += 8;
    }
    if (header.trailer) {
        end -= 4;
    }
    if (start > end) {
        return Error.MalformedPayloadRange;
    }
    return .{ .start = start, .end = end };
}

/// Vita49 Possible error types
pub const Error = error{ MalformedPayloadRange, InsufficientData };

const PacketType = enum(u4) {
    signal_wo_stream_id = 0,
    signal_w_stream_id = 1,
    ext_data_wo_stream_id = 2,
    ext_data_w_stream_id = 3,
    ctx_packet = 4,
    ext_ctx_packet = 5,
    cmd_packet = 6,
    ext_cmd_packet = 7,
};

const Tsf = enum(u2) {
    none = 0,
    sample_count = 1,
    real_time = 2,
    free_running_count = 3,
};

const Tsi = enum(u2) {
    none = 0,
    utc = 1,
    gps = 2,
    other = 3,
};

const Trailer = packed struct {
    enables: u12,
    state: u12,
    e: bool,
    ctx: u7,

    pub fn init(trailer: []const u8) Trailer {
        var tmp_array: [4]u8 = undefined;
        @memcpy(&tmp_array, trailer);
        return @bitCast(tmp_array);
    }
};

pub const Header = packed struct {
    packet_type: PacketType,
    class_id: bool,
    trailer: bool,
    tsi: Tsi,
    tsf: Tsf,
    packet_count: u4,
    packet_size: u16,

    pub fn init(stream: [4]u8) !Header {
        const little_endian_stream = std.mem.readInt(u32, &stream, .little);
        const packet_as_uint = @as(u4, @truncate((little_endian_stream >> 4) & 0xF));
        const tsi = @as(u2, @truncate((little_endian_stream >> 10) & 0x3));
        const tsf = @as(u2, @truncate((little_endian_stream >> 8) & 0x3));

        return .{
            .packet_type = @enumFromInt(packet_as_uint),
            .class_id = ((little_endian_stream >> 5) & 1) == 1,
            .trailer = ((little_endian_stream >> 6) & 1) == 1,
            .tsi = @enumFromInt(tsi),
            .tsf = @enumFromInt(tsf),
            .packet_count = @as(u4, @truncate((little_endian_stream >> 16) & 0xF)),
            .packet_size = @truncate((little_endian_stream >> 16) & 0xFFFF),
        };
    }

    pub fn output(self: Header) void {
        std.log.info("Vita49 Packet Header:\n", .{});
        std.log.info("Packet type: {}\n", .{self.packet_type});
        std.log.info("Class ID: {}\n", .{self.class_id});
        std.log.info("Trailer Present: {}\n", .{self.trailer});
        std.log.info("TSI: {}\n", .{self.tsi});
        std.log.info("TSF: {}\n", .{self.tsf});
        std.log.info("Packet_Count: {}\n", .{self.packet_count});
        std.log.info("Packet_Size: {}\n", .{self.packet_size});
    }
};

pub const ClassId = packed struct {
    reserved: u8,
    oui: u24,
    info_class_code: u16,
    packet_class_code: u16,

    pub fn init(stream: [8]u8) ClassId {
        return .{
            .reserved = stream[0],
            .oui = @as(u24, @bitCast(stream[1..4].*)),
            .info_class_code = @bitCast(stream[4..6].*),
            .packet_class_code = @bitCast(stream[6..8].*),
        };
    }
};

test "Test Vita49 Packet w/o trailer" {
    const vita49_test_packet = [_]u8{
        // Packet header
        0x3A, 0x02, 0x0A, 0x00, // Packet type (0x023A), packet size (10 words), tsi (0), and tsf(2 realtime)
        0x34, 0x12, 0x00, 0x00, // Stream ID (0x1234 in this example)

        // Class ID
        0x00, // Reserved
        0x56, 0x34, 0x12, // OUI (example: 12:34:56)
        0x78, 0x9A, // Information Class Code
        0xBC, 0xDE, // Packet Class Code

        // Timestamp - Fractional-seconds timestamp
        0x80, 0x00,
        0x00, 0x00,
        0x00, 0x00,
        0x00, 0x00,

        // Payload data (example: "Hello, VITA 49!")
        0x48, 0x65,
        0x6C, 0x6C,
        0x6F, 0x2C,
        0x20, 0x56,
        0x49, 0x54,
        0x41, 0x20,
        0x34, 0x39,
        0x21,
    };

    var vita49_packet = try Vita49.init(&vita49_test_packet, std.testing.allocator, null);
    defer vita49_packet.deinit();

    try std.testing.expectEqual(null, vita49_packet.i_timestamp);
    try std.testing.expectEqual(128, vita49_packet.f_timestamp);
    try std.testing.expectEqual(4660, vita49_packet.stream_id.?);
    try std.testing.expectEqual(1193046, vita49_packet.class_id.?.oui);
    try std.testing.expectEqualStrings("Hello, VITA 49!", vita49_packet.payload);
}

test "Test Vita49 Packet w/ trailer" {
    const vita49_test_packet = [_]u8{
        // Packet header
        0x4A, 0x06, 0x0A, 0x00, // Packet type (0x064A) and packet size (9 words)
        0x34, 0x12, 0x00, 0x00, // Stream ID (0x1234 in this example)

        // Timestamp - Integer-seconds timestamp
        0x00, 0x00, 0x00, 0x01,
        // Timestamp - Fractional-seconds timestamp
        0x80, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00,

        // Payload data (example: "Hello, VITA 49!")
        0x48, 0x65, 0x6C, 0x6C,
        0x6F, 0x2C, 0x20, 0x56,
        0x49, 0x54, 0x41, 0x20,
        0x34, 0x39, 0x21,

        // Trailer (4 bytes)
        0xAA,
        0xBB, 0xCC, 0xDD,
    };

    var vita49_packet = try Vita49.init(&vita49_test_packet, std.testing.allocator, null);
    defer vita49_packet.deinit();

    try std.testing.expectEqual(4660, vita49_packet.stream_id);
    try std.testing.expectEqual(null, vita49_packet.class_id);
    try std.testing.expectEqual(16777216, vita49_packet.i_timestamp);
    try std.testing.expectEqual(128, vita49_packet.f_timestamp);
    try std.testing.expectEqualStrings("Hello, VITA 49!", vita49_packet.payload);
}
