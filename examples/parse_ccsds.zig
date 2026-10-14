const std = @import("std");
const astroz = @import("astroz");

// First bytes of every packet in the recording
const sync = [_]u8{ 0x78, 0x97, 0xC0, 0x00, 0x00, 0x0A, 0x01, 0x02 };

pub fn main(init: std.process.Init) !void {
    // without a sync pattern, the file is parsed as back-to-back packets
    try parse(init, null);
    // with one, parsing starts at each occurrence of the pattern
    try parse(init, &sync);
}

fn parse(init: std.process.Init, syncPattern: ?[]const u8) !void {
    var parser = try astroz.Parser(astroz.Ccsds).init(null, null, 1024, init.io, init.gpa);
    defer parser.deinit();

    try parser.parseFromFile("test/ccsds.bin", syncPattern, null);

    std.log.info("sync pattern: {s}", .{if (syncPattern == null) "no" else "yes"});
    for (parser.packets.items, 0..) |pkt, n| {
        const len = 5 + @as(usize, pkt.header.packetSize);
        std.log.info("packet {d}: {x}", .{ n, pkt.rawData[0..len] });
    }
}
