//! CCSDS Data Structure.

const std = @import("std");

const Ccsds = @This();

header: HeaderMetadata,
primaryHeader: []const u8,
secondaryHeader: ?[]const u8,
packets: []const u8,
rawData: []const u8,
allocator: std.mem.Allocator,

pub fn init(pl: []const u8, allocator: std.mem.Allocator, config: ?Config) !Ccsds {
    var rawPackets = try allocator.dupe(u8, pl);
    const primaryHeader = rawPackets[0..6];
    const version = @as(u3, @truncate((primaryHeader[0] >> 5) & 0x07));
    const packetType = @as(u1, @truncate((primaryHeader[0] >> 4) & 0x01));
    const secondaryHeaderFlag = (primaryHeader[0] >> 3) & 0x01 != 0;
    const apid = @as(u11, @truncate(readBigEndianU16(.{
        primaryHeader[0] & 0x07,
        primaryHeader[1],
    })));
    const sequenceFlag = @as(u2, @truncate((primaryHeader[2] >> 6) & 0x03));
    const packetSequenceCount = @as(u14, @truncate(readBigEndianU16(.{
        primaryHeader[2] & 0x3F,
        primaryHeader[3],
    })));
    const packetSize = readBigEndianU16(.{ primaryHeader[4], primaryHeader[5] });

    var start: u8 = 6;
    const secondaryHeader: ?[]const u8 = if (secondaryHeaderFlag) blk: {
        if (rawPackets.len < 10) {
            std.log.warn("packet length is too short to have a secondary header", .{});
            break :blk null;
        }
        start = if (config != null) config.?.secondaryHeaderLength else 10;
        break :blk rawPackets[6..10];
    } else null;
    const header = HeaderMetadata{
        .version = version,
        .packetType = packetType,
        .secondaryHeaderFlag = secondaryHeaderFlag,
        .apid = apid,
        .sequenceFlag = sequenceFlag,
        .packetSequenceCount = packetSequenceCount,
        .packetSize = packetSize + 1,
    };
    const end = 5 + header.packetSize; // num of header bytes + packet_size
    const packets = rawPackets[start..end];

    _ = allocator.resize(rawPackets, end);

    return .{ .header = header, .primaryHeader = primaryHeader, .secondaryHeader = secondaryHeader, .packets = packets, .rawData = rawPackets, .allocator = allocator };
}

pub fn deinit(self: *Ccsds) void {
    self.allocator.free(self.rawData);
}

/// CCSDS Config Data Structure
/// This functionality is still so-so
/// Probably avoid it for now
pub const Config = struct {
    secondaryHeaderLength: u8,
};

/// The primary headers metadata
pub const HeaderMetadata = packed struct {
    version: u3,
    packetType: u1,
    secondaryHeaderFlag: bool,
    apid: u11,
    sequenceFlag: u2,
    packetSequenceCount: u14,
    packetSize: u16,
};

/// If you choose to use the CCSDS config you need to call this function first to get the Config struct
pub fn parseConfig(configContent: []const u8, allocator: std.mem.Allocator) !Config {
    const configParsed = try std.json.parseFromSlice(Config, allocator, configContent, .{});
    defer configParsed.deinit();

    return .{
        .secondaryHeaderLength = configParsed.value.secondaryHeaderLength,
    };
}

fn readBigEndianU16(bytes: [2]u8) u16 {
    return @as(u16, bytes[0]) << 8 | bytes[1];
}

test "CCSDS Structure Testing w/ config" {
    const test_config =
        \\{"secondaryHeaderLength": 12}
    ;
    const raw_test_packet: [16]u8 = .{ 0x78, 0x97, 0xC0, 0x00, 0x00, 0x0A, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A };
    const test_allocator = std.testing.allocator;

    const config = try parseConfig(test_config, test_allocator);
    var converted_test_packet = try Ccsds.init(&raw_test_packet, std.testing.allocator, config);
    defer converted_test_packet.deinit();

    const packets = .{ 7, 8, 9, 10 };

    try std.testing.expectEqual(3, converted_test_packet.header.version);
    try std.testing.expectEqual(1, converted_test_packet.header.packetType); // telemetry packet type
    try std.testing.expectEqual(true, converted_test_packet.header.secondaryHeaderFlag);
    try std.testing.expectEqual(3, converted_test_packet.header.sequenceFlag);
    try std.testing.expectEqual(0, converted_test_packet.header.packetSequenceCount);
    try std.testing.expectEqual(11, converted_test_packet.header.packetSize);
    try std.testing.expectEqualSlices(u8, &packets, converted_test_packet.packets);
}

test "CCSDS Structure Testing w/o config" {
    const raw_test_packet: [16]u8 = .{ 0x78, 0x97, 0xC0, 0x00, 0x00, 0x0A, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A };
    var converted_test_packet = try Ccsds.init(&raw_test_packet, std.testing.allocator, null);
    defer converted_test_packet.deinit();

    const packets = .{ 5, 6, 7, 8, 9, 10 };

    try std.testing.expectEqual(3, converted_test_packet.header.version);
    try std.testing.expectEqual(1, converted_test_packet.header.packetType); // telemetry packet type
    try std.testing.expectEqual(true, converted_test_packet.header.secondaryHeaderFlag);
    try std.testing.expectEqual(3, converted_test_packet.header.sequenceFlag);
    try std.testing.expectEqual(0, converted_test_packet.header.packetSequenceCount);
    try std.testing.expectEqual(11, converted_test_packet.header.packetSize);
    try std.testing.expectEqualSlices(u8, &packets, converted_test_packet.packets);
}
