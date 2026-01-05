const std = @import("std");
const astroz = @import("astroz");
const Ccsds = astroz.Ccsds;
const Parser = astroz.Parser;

pub fn main() !void {
    var dbga = std.heap.DebugAllocator(.{}).init;
    defer _ = dbga.deinit();
    const allocator = dbga.allocator();
    const io = std.Io.Threaded.global_single_threaded.ioBasic();

    const fileName = "./test/ccsds.bin".*;
    const syncPattern = .{ 0x78, 0x97, 0xC0, 0x00, 0x00, 0x0A, 0x01, 0x02 };

    const P = Parser(Ccsds);
    var parser = try P.init(null, null, 1024, io, allocator);
    defer parser.deinit();

    _ = try parser.parseFromFile(&fileName, &syncPattern, null);

    for (parser.packets.items) |packet| {
        std.log.info("Packets from files: 0x{x}", .{packet.packets});
    }
}
