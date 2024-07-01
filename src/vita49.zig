const std = @import("std");
const debug = std.debug;

pub const Packet_Type = enum(u4) {
    signal_w_stream_id = 0,
    signal_wo_stream_id = 1,
    ext_data_wo_steam_id = 2,
    ext_data_w_steam_id = 3,
    ctx_packet = 4,
    ext_ctx_packet = 5,
    cmd_packet = 6,
    ext_cmd_packet = 7,
};

pub const Header = packed struct {
    packet_type: Packet_Type,
    c: bool,
    indicator: u3,
    tsi: u2,
    tsf: u2,
    packet_count: u4,
    packet_size: u16,

    const Self = @This();

    pub fn new(stream: [4]u8) !Self {
        const little_endian_stream = std.mem.readInt(u32, &stream, .little);

        std.debug.print("\nbig endian = {x}\nlittle endian = {x}\n", .{ stream, little_endian_stream });

        const packet_as_uint = @as(u4, @truncate((little_endian_stream >> 28) & 0xF));
        const c = ((little_endian_stream >> 27) & 1) == 1;
        const indicator = @as(u3, @truncate((little_endian_stream >> 24) & 0x7));
        const tsi = @as(u2, @truncate((little_endian_stream >> 22) & 0x3));
        const tsf = @as(u2, @truncate((little_endian_stream >> 20) & 0x3));
        const packet_count = @as(u4, @truncate((little_endian_stream >> 14) & 0xF));
        const packet_size = @as(u16, @truncate(little_endian_stream & 0xFFFF));

        return .{ .packet_type = @enumFromInt(packet_as_uint), .c = c, .indicator = indicator, .tsi = tsi, .tsf = tsf, .packet_count = packet_count, .packet_size = packet_size };
    }
};

pub const Vita49 = struct {
    header: Header,
    stream_id: ?u32,
    class_id: ?u64,
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

test "Test New Vita49 Header" {
    const test_packet = [_]u8{ 0x1C, 0x00, 0x30, 0x18 };

    const header = try Header.new(test_packet);

    std.debug.print("\nHEADER IS: {any}\n", .{header});
}
