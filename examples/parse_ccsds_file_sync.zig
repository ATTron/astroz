const std = @import("std");
const astroz = @import("astroz");

const sync = [_]u8{ 0x78, 0x97, 0xC0, 0x00, 0x00, 0x0A, 0x01, 0x02 };

pub fn main(init: std.process.Init) !void {
    var parser = try astroz.Parser(astroz.Ccsds).init(null, null, 1024, init.io, init.gpa);
    defer parser.deinit();

    try parser.parseFromFile("test/ccsds.bin", &sync, null);

    for (parser.packets.items, 0..) |pkt, n| {
        const len = 5 + @as(usize, pkt.header.packetSize);
        std.log.info("packet {d}: {x}", .{ n, pkt.rawData[0..len] });
    }
}
