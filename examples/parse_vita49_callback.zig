const std = @import("std");
const astroz = @import("astroz");

const sync = [_]u8{ 0x3A, 0x02, 0x0A, 0x00, 0x34, 0x12, 0x00, 0x00, 0x00, 0x56 };
const expected = "Hello, VITA 49!";

fn onPacket(pkt: astroz.Vita49) void {
    std.log.info("callback: stream {?x}, {d} words, payload \"{s}\"", .{ pkt.streamId, pkt.header.packetSize, pkt.payload });
}

pub fn main(init: std.process.Init) !void {
    var parser = try astroz.Parser(astroz.Vita49).init(null, null, 1024, init.io, init.gpa);
    defer parser.deinit();

    try parser.parseFromFile("test/vita49.bin", &sync, onPacket);

    for (parser.packets.items) |pkt| {
        if (!std.mem.eql(u8, pkt.payload, expected)) return error.UnexpectedPayload;
        std.log.info("payload: {s}", .{pkt.payload});
    }
}
