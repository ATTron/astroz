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

    const P = Parser(Ccsds);
    var parser = try P.init(null, null, 1024, io, allocator);
    defer parser.deinit();

    _ = try parser.parseFromFile(&fileName, null, null);

    for (parser.packets.items) |packet| {
        std.log.info("Packets from files: 0x{x}", .{packet.packets});
    }
}
