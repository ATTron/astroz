const std = @import("std");
const astroz = @import("astroz");
const Ccsds = astroz.Ccsds;
const Config = Ccsds.Config;

pub fn main() !void {
    var dbga = std.heap.DebugAllocator(.{}).init;
    defer _ = dbga.deinit();
    const allocator = dbga.allocator();

    const configFile = try std.fs.cwd().readFileAlloc(allocator, "examples/create_ccsds_packet_config.json", 512);
    defer allocator.free(configFile);

    const config = try Ccsds.parseConfig(configFile, allocator);

    const rawTestPacket: [16]u8 = .{ 0x78, 0x97, 0xC0, 0x00, 0x00, 0x0A, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A };
    var convertedTestPacket = try Ccsds.init(&rawTestPacket, allocator, config);
    defer convertedTestPacket.deinit();

    std.debug.print("\nCCSDS Packet Created:\n{any}", .{convertedTestPacket});
}
