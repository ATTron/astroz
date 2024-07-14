const std = @import("std");
const astroz = @import("astroz");
const ccsds = astroz.ccsds;
const CCSDS = ccsds.CCSDS;
const Config = ccsds.Config;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const config_file = try std.fs.cwd().readFileAlloc(allocator, "examples/create_ccsds_packet_config.json", 512);
    defer allocator.free(config_file);

    const config = try ccsds.parseConfig(config_file, allocator);

    const raw_test_packet: [16]u8 = .{ 0x78, 0x97, 0xC0, 0x00, 0x00, 0x0A, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A };
    var converted_test_packet = try CCSDS.init(&raw_test_packet, allocator, config);
    defer converted_test_packet.deinit();

    std.debug.print("\nCCSDS Packet Created:\n{any}", .{converted_test_packet});
}
