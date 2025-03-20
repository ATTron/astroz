const std = @import("std");
const astroz = @import("astroz");
const Vita49 = astroz.Vita49;
const Parser = astroz.Parser;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const file_name = "./test/vita49.bin".*;
    const sync_pattern = .{ 0x3A, 0x02, 0x0a, 0x00, 0x34, 0x12, 0x00, 0x00, 0x00, 0x56 };
    const P = Parser(Vita49);
    var parser = try P.init(null, null, 1024, allocator);
    defer parser.deinit();

    _ = try parser.parseFromFile(&file_name, &sync_pattern, callback);
    for (parser.packets.items) |packet| {
        try std.testing.expectEqualStrings("Hello, VITA 49!", packet.payload);
    }
}

fn callback(packet: Vita49) void {
    std.debug.print("Packet received: {any}", .{packet});
}
