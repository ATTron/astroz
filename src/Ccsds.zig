//! CCSDS Data Structure.

const std = @import("std");

const Ccsds = @This();

header: HeaderMetadata,
primary_header: []const u8,
secondary_header: ?[]const u8,
packets: []const u8,
raw_data: []const u8,
allocator: std.mem.Allocator,

pub fn init(pl: []const u8, allocator: std.mem.Allocator, config: ?Config) !Ccsds {
    var raw_packets = try allocator.dupe(u8, pl);
    const primary_header = raw_packets[0..6];
    const version = @as(u3, @truncate((primary_header[0] >> 5) & 0x07));
    const packet_type = @as(u1, @truncate((primary_header[0] >> 4) & 0x01));
    const secondary_header_flag = (primary_header[0] >> 3) & 0x01 != 0;
    const apid = @as(u11, @truncate(readBigEndianU16(.{
        primary_header[0] & 0x07,
        primary_header[1],
    })));
    const sequence_flag = @as(u2, @truncate((primary_header[2] >> 6) & 0x03));
    const packet_sequence_count = @as(u14, @truncate(readBigEndianU16(.{
        primary_header[2] & 0x3F,
        primary_header[3],
    })));
    const packet_size = readBigEndianU16(.{ primary_header[4], primary_header[5] });

    var start: u8 = 6;
    const secondary_header: ?[]const u8 = if (secondary_header_flag) blk: {
        if (raw_packets.len < 10) {
            std.log.warn("packet length is too short to have a secondary header", .{});
            break :blk null;
        }
        start = if (config != null) config.?.secondary_header_length else 10;
        break :blk raw_packets[6..10];
    } else null;
    const header = HeaderMetadata{
        .version = version,
        .packet_type = packet_type,
        .secondary_header_flag = secondary_header_flag,
        .apid = apid,
        .sequence_flag = sequence_flag,
        .packet_sequence_count = packet_sequence_count,
        .packet_size = packet_size + 1,
    };
    const end = 5 + header.packet_size; // num of header bytes + packet_size
    const packets = raw_packets[start..end];

    _ = allocator.resize(raw_packets, end);

    return .{ .header = header, .primary_header = primary_header, .secondary_header = secondary_header, .packets = packets, .raw_data = raw_packets, .allocator = allocator };
}

pub fn deinit(self: *Ccsds) void {
    self.allocator.free(self.raw_data);
}

/// CCSDS Config Data Structure
/// This functionality is still so-so
/// Probably avoid it for now
pub const Config = struct {
    secondary_header_length: u8,
};

/// The primary headers metadata
pub const HeaderMetadata = packed struct {
    version: u3,
    packet_type: u1,
    secondary_header_flag: bool,
    apid: u11,
    sequence_flag: u2,
    packet_sequence_count: u14,
    packet_size: u16,
};

/// If you choose to use the CCSDS config you need to call this function first to get the Config struct
pub fn parseConfig(config_content: []const u8, allocator: std.mem.Allocator) !Config {
    const config_parsed = try std.json.parseFromSlice(Config, allocator, config_content, .{});
    defer config_parsed.deinit();

    return .{
        .secondary_header_length = config_parsed.value.secondary_header_length,
    };
}

fn readBigEndianU16(bytes: [2]u8) u16 {
    return @as(u16, bytes[0]) << 8 | bytes[1];
}

test "CCSDS Structure Testing w/ config" {
    const test_config =
        \\{"secondary_header_length": 12}
    ;
    const raw_test_packet: [16]u8 = .{ 0x78, 0x97, 0xC0, 0x00, 0x00, 0x0A, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A };
    const test_allocator = std.testing.allocator;

    const config = try parseConfig(test_config, test_allocator);
    var converted_test_packet = try Ccsds.init(&raw_test_packet, std.testing.allocator, config);
    defer converted_test_packet.deinit();

    const packets = .{ 7, 8, 9, 10 };

    try std.testing.expectEqual(3, converted_test_packet.header.version);
    try std.testing.expectEqual(1, converted_test_packet.header.packet_type); // telemetry packet type
    try std.testing.expectEqual(true, converted_test_packet.header.secondary_header_flag);
    try std.testing.expectEqual(3, converted_test_packet.header.sequence_flag);
    try std.testing.expectEqual(0, converted_test_packet.header.packet_sequence_count);
    try std.testing.expectEqual(11, converted_test_packet.header.packet_size);
    try std.testing.expectEqualSlices(u8, &packets, converted_test_packet.packets);
}

test "CCSDS Structure Testing w/o config" {
    const raw_test_packet: [16]u8 = .{ 0x78, 0x97, 0xC0, 0x00, 0x00, 0x0A, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A };
    var converted_test_packet = try Ccsds.init(&raw_test_packet, std.testing.allocator, null);
    defer converted_test_packet.deinit();

    const packets = .{ 5, 6, 7, 8, 9, 10 };

    try std.testing.expectEqual(3, converted_test_packet.header.version);
    try std.testing.expectEqual(1, converted_test_packet.header.packet_type); // telemetry packet type
    try std.testing.expectEqual(true, converted_test_packet.header.secondary_header_flag);
    try std.testing.expectEqual(3, converted_test_packet.header.sequence_flag);
    try std.testing.expectEqual(0, converted_test_packet.header.packet_sequence_count);
    try std.testing.expectEqual(11, converted_test_packet.header.packet_size);
    try std.testing.expectEqualSlices(u8, &packets, converted_test_packet.packets);
}
