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

        std.debug.print("\nshowing entire stream: {x}\n", .{little_endian_stream});
        std.debug.print("\nshowing tsi area probably: {x}\n", .{little_endian_stream >> 22});

        return .{ .packet_type = @enumFromInt(packet_as_uint), .class_id = ((little_endian_stream >> 27) & 1) == 1, .trailer = ((little_endian_stream >> 24) & 1) == 1, .tsi = @enumFromInt(tsi), .tsf = @enumFromInt(tsf), .packet_count = @as(u4, @truncate((little_endian_stream >> 14) & 0xF)), .packet_size = @as(u16, @truncate(little_endian_stream & 0xFFFF)) };
    }

    pub fn output(self: Self) void {
        std.debug.print("Vita49 Packet Header:\n", .{});
        std.debug.print("Packet type: {}\n", .{self.packet_type});
        std.debug.print("Class ID: {}\n", .{self.class_id});
        std.debug.print("Trailer Present: {}\n", .{self.trailer});
        std.debug.print("TSI: {}\n", .{self.tsi});
        std.debug.print("TSF: {}\n", .{self.tsf});
        std.debug.print("Packet_Count: {}\n", .{self.packet_count});
        std.debug.print("Packet_Size: {}\n", .{self.packet_size});
    }
};

pub const Vita49 = struct {
    header: Header,
    stream_id: ?u32,
    class_id: ?Class_ID,
    i_timestamp: ?u32,
    f_timestampt: ?u64,
    payload: []const u8,
    trailer: ?u32,

    const Self = @This();

    pub fn new(stream: []const u8) Self {
        const header = Header.new(stream[0..31]);
        var stream_id: ?32 = undefined;
        switch (header.packet_type) {
            Packet_Type.ctx_packet, Packet_Type.ext_ctx_packet, Packet_Type.signal_w_stream_id, Packet_Type.ext_data_w_steam_id => {
                stream_id = stream[32..63];
            },
            Packet_Type.signal_wo_stream_id, Packet_Type.signal_wo_stream_id, Packet_Type.ext_data_wo_steam_id => {
                stream_id = null;
            },
            _ => {
                std.debug.print("Not currently implemented", .{});
            },
        }
    }
};

pub const Vita49_Parser = struct { stream: []const u8, packets: []const Vita49 };

test "Test Vita49 Header" {
    const test_packet = [_]u8{ 0x1C, 0x00, 0xC0, 0x18 };

    const header = try Header.new(test_packet);
    header.output();
}
