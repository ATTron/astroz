const std = @import("std");
const debug = std.debug;

pub const Vita49Error = error{ MalformedPayloadRange, InsufficentData };

pub const Packet_Type = enum(u4) {
    signal_wo_stream_id = 0,
    signal_w_stream_id = 1,
    ext_data_wo_stream_id = 2,
    ext_data_w_stream_id = 3,
    ctx_packet = 4,
    ext_ctx_packet = 5,
    cmd_packet = 6,
    ext_cmd_packet = 7,
};

pub const TSF = enum(u2) {
    none = 0,
    sample_count = 1,
    real_time = 2,
    free_running_count = 3,
};

pub const TSI = enum(u2) {
    none = 0,
    utc = 1,
    gps = 2,
    other = 3,
};

pub const Trailer = packed struct {
    enables: u12,
    state: u12,
    e: bool,
    ctx: u7,

    const Self = @This();

    pub fn new(trailer: []const u8) Self {
        var tmp_array: [4]u8 = undefined;
        @memcpy(&tmp_array, trailer);
        return @bitCast(tmp_array);
    }
};

pub const Header = packed struct {
    packet_type: Packet_Type,
    class_id: bool,
    trailer: bool,
    tsi: TSI,
    tsf: TSF,
    packet_count: u4,
    packet_size: u16,

    const Self = @This();

    pub fn new(stream: [4]u8) !Self {
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

    pub fn output(self: Self) void {
        debug.print("Vita49 Packet Header:\n", .{});
        debug.print("Packet type: {}\n", .{self.packet_type});
        debug.print("Class ID: {}\n", .{self.class_id});
        debug.print("Trailer Present: {}\n", .{self.trailer});
        debug.print("TSI: {}\n", .{self.tsi});
        debug.print("TSF: {}\n", .{self.tsf});
        debug.print("Packet_Count: {}\n", .{self.packet_count});
        debug.print("Packet_Size: {}\n", .{self.packet_size});
    }
};

pub const Class_ID = packed struct {
    reserved: u8,
    oui: u24,
    info_class_code: u16,
    packet_class_code: u16,

    const Self = @This();

    pub fn new(stream: [8]u8) Self {
        return .{
            .reserved = stream[0],
            .oui = @as(u24, @bitCast(stream[1..4].*)),
            .info_class_code = @bitCast(stream[4..6].*),
            .packet_class_code = @bitCast(stream[6..8].*),
        };
    }
};

pub const Vita49 = struct {
    header: Header,
    stream_id: ?u32,
    class_id: ?Class_ID,
    i_timestamp: ?u32,
    f_timestamp: ?u64,
    payload: []const u8,
    trailer: ?Trailer,

    const Self = @This();

    pub fn new(stream: []const u8) !Self {
        if (stream.len < 4) {
            return Vita49Error.InsufficentData;
        }
        const header = try Header.new(stream[0..4].*);
        var stream_id: ?u32 = undefined;
        var class_id: ?Class_ID = undefined;
        var payload: []const u8 = undefined;
        var trailer: ?Trailer = undefined;
        var i_start: usize = 4;
        var f_start: usize = 4;
        var i_timestamp: ?u32 = undefined;
        var f_timestamp: ?u64 = undefined;
        var payload_range = try get_payload_range(header, true);

        switch (header.packet_type) {
            Packet_Type.signal_w_stream_id,
            Packet_Type.ext_data_w_stream_id,
            Packet_Type.ext_cmd_packet,
            Packet_Type.cmd_packet,
            Packet_Type.ctx_packet,
            Packet_Type.ext_ctx_packet,
            => {
                stream_id = std.mem.readInt(u32, stream[4..8], .little);
                i_start += 4;
                f_start += 4;
            },
            Packet_Type.signal_wo_stream_id,
            Packet_Type.ext_data_wo_stream_id,
            => {
                payload_range = try get_payload_range(header, false);
            },
        }
        if (header.class_id) {
            class_id = Class_ID.new(stream[8..16].*);
            i_start += 8;
            f_start += 8;
        } else {
            class_id = null;
        }
        if (header.trailer) {
            trailer = Trailer.new(stream[payload_range.end..]);
            payload = stream[payload_range.start..payload_range.end];
        } else {
            payload = stream[payload_range.start..];
        }
        if (header.tsi != TSI.none) {
            var tmp_array: [4]u8 = undefined;
            @memcpy(&tmp_array, stream[i_start .. i_start + 4]);
            i_timestamp = std.mem.readInt(u32, &tmp_array, .little);
            f_start += 4;
        } else {
            i_timestamp = null;
        }
        if (header.tsf != TSF.none) {
            var tmp_array: [8]u8 = undefined;
            @memcpy(&tmp_array, stream[f_start .. f_start + 8]);
            f_timestamp = std.mem.readInt(u64, &tmp_array, .little);
        } else {
            f_timestamp = null;
        }
        return .{
            .header = header,
            .stream_id = stream_id,
            .class_id = class_id,
            .i_timestamp = i_timestamp,
            .f_timestamp = f_timestamp,
            .payload = payload,
            .trailer = trailer,
        };
    }

    fn get_payload_range(header: Header, has_stream_id: bool) !struct { start: usize, end: usize } {
        var start: usize = 4;
        var end: usize = @as(usize, (header.packet_size * 4)) - 1;
        if (has_stream_id) {
            start += 4;
        }
        if (header.class_id) {
            start += 8;
        }
        if (header.tsi != TSI.none) {
            start += 4;
        }
        if (header.tsf != TSF.none) {
            start += 8;
        }
        if (header.trailer) {
            end -= 4;
        }
        if (start > end) {
            return Vita49Error.MalformedPayloadRange;
        }
        return .{ .start = start, .end = end };
    }
};

pub const Vita49_Parser = struct { stream: []const u8, packets: []const Vita49 };

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

    const vita49_packet = try Vita49.new(&vita49_test_packet);

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

    const vita49_packet = try Vita49.new(&vita49_test_packet);

    try std.testing.expectEqual(4660, vita49_packet.stream_id);
    try std.testing.expectEqual(null, vita49_packet.class_id);
    try std.testing.expectEqual(16777216, vita49_packet.i_timestamp);
    try std.testing.expectEqual(128, vita49_packet.f_timestamp);
    try std.testing.expectEqualStrings("Hello, VITA 49!", vita49_packet.payload);
}