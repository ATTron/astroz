const std = @import("std");
const astroz = @import("astroz");
const CCSDS = astroz.ccsds.CCSDS;
const Parser = astroz.parsers.Parser;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const file_name = "./test/ccsds.bin".*;
    const sync_pattern = .{ 0x78, 0x97, 0xC0, 0x00, 0x00, 0x0A, 0x01, 0x02 };

    const P = Parser(CCSDS);
    var parser = try P.init(null, null, 1024, allocator);
    defer parser.deinit();

    _ = try parser.parseFromFile(&file_name, &sync_pattern, null);

    for (parser.packets.items) |packet| {
        std.log.info("Packets from files: 0x{x}", .{packet.packets});
    }
}
