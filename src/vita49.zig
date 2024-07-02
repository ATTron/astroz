const std = @import("std");
const debug = std.debug;

pub const Packet_Type = enum(u4) {
    signal_wo_stream_id = 0,
    signal_w_stream_id = 1,
    ext_data_wo_steam_id = 2,
    ext_data_w_steam_id = 3,
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
    calibrated_time: bool,
    valid_data: bool,
    reference_lock: bool,
    agc_mgc: bool,
    detected_signal: bool,
    spectral_inversion: bool,
    over_range: bool,
    sample_loss: bool,
    user_defined_0: bool,
    user_defined_1: bool,
    user_defined_2: bool,
    user_defined_3: bool,
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

        debug.print("\nshowing entire stream: {x}\n", .{little_endian_stream});
        debug.print("\nshowing packet_type area probably: {x}\n", .{packet_as_uint});
        debug.print("\nshowing class_id area probably: {x}\n", .{(little_endian_stream >> 5) & 1});
        debug.print("\nshowing trailer area probably: {x}\n", .{(little_endian_stream >> 6) & 1});
        debug.print("\nshowing packet count area probably: {x}\n", .{(little_endian_stream >> 16) & 0xF});
        debug.print("\nshowing packet size area probably: {x}\n", .{(little_endian_stream >> 16) & 0xFFFF});
        // debug.print("\nshowing tsi area probably: {x}\n", .{little_endian_stream >> 10});
        // debug.print("\nshowing packet size area probably: {x}\n{d}\n", .{ little_endian_stream >> 16, little_endian_stream >> 16 });

        return .{ .packet_type = @enumFromInt(packet_as_uint), .class_id = ((little_endian_stream >> 5) & 1) == 1, .trailer = ((little_endian_stream >> 6) & 1) == 1, .tsi = @enumFromInt(tsi), .tsf = @enumFromInt(tsf), .packet_count = @as(u4, @truncate((little_endian_stream >> 16) & 0xF)), .packet_size = @truncate((little_endian_stream >> 16) & 0xFFFF) };
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
            // .oui = @as(u24, @truncate(@as(u24, @bitCast(stream[1..4].*)))),
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
            return error.InsufficientData;
        }
        const header = try Header.new(stream[0..4].*);
        debug.print("\n\nshowing header packet size: {d}", .{header.packet_size});
        const packet = stream[4 .. (header.packet_size * 4) - 1]; // get the rest of the packet... need to multiply by 4 to get in bytes
        var stream_id: ?u32 = undefined;
        var class_id: ?Class_ID = undefined;
        switch (header.packet_type) {
            Packet_Type.signal_w_stream_id, Packet_Type.ext_data_w_steam_id => {
                stream_id = std.mem.readInt(u32, packet[4..8], .little);
                class_id = Class_ID.new(packet[8..16].*);
                debug.print("\nshowing class_id: {any}\n", .{class_id});
            },
            Packet_Type.signal_wo_stream_id, Packet_Type.ext_data_wo_steam_id => {
                debug.print("Not currently implemented", .{});
            },
            else => {
                debug.print("Not currently implemented", .{});
            },
        }
        return .{
            .header = header,
            .stream_id = stream_id,
            .class_id = class_id,
            .i_timestamp = null,
            .f_timestamp = null,
            .payload = stream[24..],
            .trailer = null,
        };
    }

    pub fn read_trailer(self: Self) void {
        switch (self.trailer.?) {}
    }
};

pub const Vita49_Parser = struct { stream: []const u8, packets: []const Vita49 };

// test "Test Vita49 Header" {
//     const test_header = [_]u8{ 0x18, 0xC3, 0x00, 0x00 };
//
//     const header = try Header.new(test_header);
//     header.output();
// }

// TODO: Fix this
test "Test Vita49 Packet w/ no trailer" {
    const vita49_test_packet = [_]u8{
        // Packet header
        0x3A, 0x02, 0x0A, 0x00, // Packet type (0x4A02) and packet size (11 words)
        0x00, 0x00, 0x12, 0x34, // Stream ID (0x1234 in this example)

        // Class ID
        0x00, // Reserved
        0x12, 0x34, 0x56, // OUI (example: 12:34:56)
        0x78, 0x9A, // Information Class Code
        0xBC, 0xDE, // Packet Class Code

        // Timestamp - Integer-seconds timestamp
        0x00, 0x00,
        0x00, 0x01,
        // Timestamp - Fractional-seconds timestamp
        0x80, 0x00,
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

    debug.print("\nvita49 packet: {any}\n", .{vita49_packet});
    debug.print("\ndebug output: {s}\n", .{vita49_packet.payload});
    try std.testing.expectEqual(4660, vita49_packet.stream_id.?);
    try std.testing.expectEqual(15715755, vita49_packet.class_id.?.oui);
}
