const std = @import("std");
const astroz = @import("astroz");
const CCSDS = astroz.ccsds.CCSDS;
const Parser = astroz.parsers.Parser;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const file_name = "./test/ccsds.bin".*;

    const P = Parser(CCSDS);
    var parser = try P.init(null, null, 1024, allocator);
    defer parser.deinit();

    _ = try parser.parseFromFile(&file_name, null, null);

    for (parser.packets.items) |packet| {
        std.log.info("Packets from files: 0x{x}", .{packet.packets});
    }
}