//! The Vita49 packet type which will be used by the Parser

const std = @import("std");

const Vita49 = @This();

header: Header,
streamId: ?u32,
classId: ?ClassId,
iTimestamp: ?u32,
fTimestamp: ?u64,
payload: []const u8,
trailer: ?Trailer,
allocator: std.mem.Allocator,
end: usize,

// private
rawData: []const u8,

pub fn init(pl: []const u8, allocator: std.mem.Allocator, config: ?[]const u8) !Vita49 {
    var stream = try allocator.dupe(u8, pl); // dupe the data so we dont ever lose it
    if (config != null) {
        std.log.debug("Config found for Vita49 but this hasn't been implemented yet", .{});
    }
    if (stream.len < 4) {
        return Error.InsufficientData;
    }
    const header = try Header.init(stream[0..4].*);
    var streamId: ?u32 = undefined;
    var classId: ?ClassId = undefined;
    var payload: []const u8 = undefined;
    var trailer: ?Trailer = undefined;
    var iStart: usize = 4;
    var fStart: usize = 4;
    var iTimestamp: ?u32 = undefined;
    var fTimestamp: ?u64 = undefined;
    var payloadRange = try getPayloadRange(header, true);

    switch (header.packetType) {
        .signal_w_stream_id,
        .ext_data_w_stream_id,
        .ext_cmd_packet,
        .cmd_packet,
        .ctx_packet,
        .ext_ctx_packet,
        => {
            streamId = std.mem.readInt(u32, stream[4..8], .little);
            iStart += 4;
            fStart += 4;
        },
        .signal_wo_stream_id,
        .ext_data_wo_stream_id,
        => {
            payloadRange = try getPayloadRange(header, false);
        },
    }
    if (header.classId) {
        classId = ClassId.init(stream[8..16].*);
        iStart += 8;
        fStart += 8;
    } else {
        classId = null;
    }
    if (header.trailer) {
        trailer = Trailer.init(stream[payloadRange.end..]);
        payload = stream[payloadRange.start..payloadRange.end];
    } else {
        payload = stream[payloadRange.start..payloadRange.end];
    }
    if (header.tsi != Tsi.none) {
        var tmpArray: [4]u8 = undefined;
        @memcpy(&tmpArray, stream[iStart .. iStart + 4]);
        iTimestamp = std.mem.readInt(u32, &tmpArray, .little);
        fStart += 4;
    } else {
        iTimestamp = null;
    }
    if (header.tsf != Tsf.none) {
        var tmpArray: [8]u8 = undefined;
        @memcpy(&tmpArray, stream[fStart .. fStart + 8]);
        fTimestamp = std.mem.readInt(u64, &tmpArray, .little);
    } else {
        fTimestamp = null;
    }
    const endOfData = if (header.trailer) payloadRange.end + 4 else payloadRange.end;

    _ = allocator.resize(stream, endOfData);
    return .{
        .header = header,
        .streamId = streamId,
        .classId = classId,
        .iTimestamp = iTimestamp,
        .fTimestamp = fTimestamp,
        .payload = payload,
        .trailer = trailer,
        .allocator = allocator,
        .end = payloadRange.end,
        .rawData = stream,
    };
}

pub fn deinit(self: *Vita49) void {
    self.allocator.free(self.rawData);
}

fn getPayloadRange(header: Header, hasStreamId: bool) !struct { start: usize, end: usize } {
    var start: usize = 4;
    var end: usize = @as(usize, (header.packetSize * 4)) - 1;
    if (hasStreamId) {
        start += 4;
    }
    if (header.classId) {
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
        var tmpArray: [4]u8 = undefined;
        @memcpy(&tmpArray, trailer);
        return @bitCast(tmpArray);
    }
};

pub const Header = packed struct {
    packetType: PacketType,
    classId: bool,
    trailer: bool,
    tsi: Tsi,
    tsf: Tsf,
    packetCount: u4,
    packetSize: u16,

    pub fn init(stream: [4]u8) !Header {
        const littleEndianStream = std.mem.readInt(u32, &stream, .little);
        const packetAsUint = @as(u4, @truncate((littleEndianStream >> 4) & 0xF));
        const tsi = @as(u2, @truncate((littleEndianStream >> 10) & 0x3));
        const tsf = @as(u2, @truncate((littleEndianStream >> 8) & 0x3));

        return .{
            .packetType = @enumFromInt(packetAsUint),
            .classId = ((littleEndianStream >> 5) & 1) == 1,
            .trailer = ((littleEndianStream >> 6) & 1) == 1,
            .tsi = @enumFromInt(tsi),
            .tsf = @enumFromInt(tsf),
            .packetCount = @as(u4, @truncate((littleEndianStream >> 16) & 0xF)),
            .packetSize = @truncate((littleEndianStream >> 16) & 0xFFFF),
        };
    }

    pub fn output(self: Header) void {
        std.log.info("Vita49 Packet Header:\n", .{});
        std.log.info("Packet type: {}\n", .{self.packetType});
        std.log.info("Class ID: {}\n", .{self.classId});
        std.log.info("Trailer Present: {}\n", .{self.trailer});
        std.log.info("TSI: {}\n", .{self.tsi});
        std.log.info("TSF: {}\n", .{self.tsf});
        std.log.info("Packet_Count: {}\n", .{self.packetCount});
        std.log.info("Packet_Size: {}\n", .{self.packetSize});
    }
};

pub const ClassId = packed struct {
    reserved: u8,
    oui: u24,
    infoClassCode: u16,
    packetClassCode: u16,

    pub fn init(stream: [8]u8) ClassId {
        return .{
            .reserved = stream[0],
            .oui = @as(u24, @bitCast(stream[1..4].*)),
            .infoClassCode = @bitCast(stream[4..6].*),
            .packetClassCode = @bitCast(stream[6..8].*),
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

    try std.testing.expectEqual(null, vita49_packet.iTimestamp);
    try std.testing.expectEqual(128, vita49_packet.fTimestamp);
    try std.testing.expectEqual(4660, vita49_packet.streamId.?);
    try std.testing.expectEqual(1193046, vita49_packet.classId.?.oui);
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

    try std.testing.expectEqual(4660, vita49_packet.streamId);
    try std.testing.expectEqual(null, vita49_packet.classId);
    try std.testing.expectEqual(16777216, vita49_packet.iTimestamp);
    try std.testing.expectEqual(128, vita49_packet.fTimestamp);
    try std.testing.expectEqualStrings("Hello, VITA 49!", vita49_packet.payload);
}
