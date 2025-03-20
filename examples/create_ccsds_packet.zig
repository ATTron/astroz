const std = @import("std");
const astroz = @import("astroz");
const Ccsds = astroz.Ccsds;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const rawTestPacket: [16]u8 = .{ 0x78, 0x97, 0xC0, 0x00, 0x00, 0x0A, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A };
    var convertedTestPacket = try Ccsds.init(&rawTestPacket, allocator, null);
    defer convertedTestPacket.deinit();

    std.debug.print("CCSDS Packet Created:\n{any}", .{convertedTestPacket});
}
