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

pub const Class_ID = packed struct {
    reserved: u8,
    oui: u24,
    info_class_code: u16,
    packet_class_code: u16,
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
        const packet_as_uint = @as(u4, @truncate((little_endian_stream >> 28) & 0xF));
        const tsi = @as(u2, @truncate((little_endian_stream >> 22) & 0x3));
        const tsf = @as(u2, @truncate((little_endian_stream >> 24) & 0x3));

        debug.print("\nshowing entire stream: {x}\n", .{little_endian_stream});
        debug.print("\nshowing tsi area probably: {x}\n", .{little_endian_stream >> 22});

        return .{ .packet_type = @enumFromInt(packet_as_uint), .class_id = ((little_endian_stream >> 27) & 1) == 1, .trailer = ((little_endian_stream >> 24) & 1) == 1, .tsi = @enumFromInt(tsi), .tsf = @enumFromInt(tsf), .packet_count = @as(u4, @truncate((little_endian_stream >> 14) & 0xF)), .packet_size = @as(u16, @truncate(little_endian_stream & 0xFFFF)) };
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

pub const Vita49 = struct {
    header: Header,
    stream_id: ?u32,
    class_id: ?Class_ID,
    i_timestamp: ?u32,
    f_timestampt: ?u64,
    payload: []const u8,
    trailer: ?Trailer,

    const Self = @This();

    pub fn new(stream: []const u8) !Self {
        const header = try Header.new(stream[0..4]);
        var stream_id: ?32 = undefined;
        switch (header.packet_type) {
            Packet_Type.ctx_packet, Packet_Type.ext_ctx_packet, Packet_Type.signal_w_stream_id, Packet_Type.ext_data_w_steam_id => {
                stream_id = stream[32..63];
            },
            Packet_Type.signal_wo_stream_id, Packet_Type.signal_wo_stream_id, Packet_Type.ext_data_wo_steam_id => {},
            _ => {
                debug.print("Not currently implemented", .{});
            },
        }
    }

    pub fn read_trailer(self: Self) void {
        switch (self.trailer.?) {}
    }
};

pub const Vita49_Parser = struct { stream: []const u8, packets: []const Vita49 };

test "Test Vita49 Header" {
    const test_header = [_]u8{ 0x18, 0xC3, 0x00, 0x00 };

    const header = try Header.new(test_header);
    header.output();
}

// TODO: Fix this
// test "Test Vita49 Packet" {
//     const test_packet = [_]u8{
//         // Packet header
//         0x18, 0xC3, 0x00, 0x00, // Packet type and header size
//         0x00, 0x00, 0x00, 0x01, // Stream ID
//
//         // Class ID
//         0x00, 0x00, 0x00, 0x00, // Reserved
//         0xAB, 0xCD, 0xEF, // OUI (example: AB:CD:EF)
//         0x12, 0x34, // Information Class Code
//         0x56, 0x78, // Packet Class Code
//
//         // Timestamp - Integer-seconds timestamp
//         0x00, 0x00,
//         0x00, 0x00,
//         // Timestamp - Fractional-seconds timestamp
//         0x00, 0x00,
//         0x00, 0x00,
//
//         // Payload data (example: 16 bytes)
//         0x11, 0x22,
//         0x33, 0x44,
//         0x55, 0x66,
//         0x77, 0x88,
//         0x99, 0xAA,
//         0xBB, 0xCC,
//         0xDD, 0xEE,
//         0xFF,
//         0x00,
//
//         // Trailer
//         0x00, 0x00, 0x00, 0x00, // CRC-32C
//     };
//
//     const vita49_packet = Vita49.new(&test_packet);
//
//     debug.print("\nvita49 packet: {any}\n", .{vita49_packet});
// }
