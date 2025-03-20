const std = @import("std");
const astroz = @import("astroz");
const Ccsds = astroz.Ccsds;
const Parser = astroz.Parser;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const fileName = "./test/ccsds.bin".*;

    const P = Parser(Ccsds);
    var parser = try P.init(null, null, 1024, allocator);
    defer parser.deinit();

    _ = try parser.parseFromFile(&fileName, null, null);

    for (parser.packets.items) |packet| {
        std.log.info("Packets from files: 0x{x}", .{packet.packets});
    }
}
