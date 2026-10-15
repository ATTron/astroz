//! Space packet formats: build a CCSDS packet (with and without a config that
//! sets the secondary header), parse a CCSDS recording with and without a sync
//! pattern, and parse a VITA 49 stream with a per-packet callback.

const std = @import("std");
const astroz = @import("astroz");
const Ccsds = astroz.Ccsds;
const Vita49 = astroz.Vita49;

const ccsdsBytes = [_]u8{ 0x78, 0x97, 0xC0, 0x00, 0x00, 0x0A, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A };
const ccsdsSync = ccsdsBytes[0..8];
const vitaSync = [_]u8{ 0x3A, 0x02, 0x0A, 0x00, 0x34, 0x12, 0x00, 0x00, 0x00, 0x56 };

fn onVita49(pkt: Vita49) void {
    std.log.info("  stream {?x}, {d} words, payload \"{s}\"", .{ pkt.streamId, pkt.header.packetSize, pkt.payload });
}

pub fn main(init: std.process.Init) !void {
    const gpa = init.gpa;

    var pkt = try Ccsds.init(&ccsdsBytes, gpa, null);
    defer pkt.deinit();
    std.log.info("CCSDS header: {any}", .{pkt.header});
    std.log.info("  payload: {x}", .{pkt.packets});

    const json = try std.Io.Dir.cwd().readFileAlloc(init.io, "examples/ccsds_config.json", gpa, .limited(4096));
    defer gpa.free(json);
    const config = try Ccsds.parseConfig(json, gpa);
    var configured = try Ccsds.init(&ccsdsBytes, gpa, config);
    defer configured.deinit();
    std.log.info("  with {any}: payload {x}", .{ config, configured.packets });

    // Without a sync pattern the file is read as back-to-back packets; with one,
    // parsing starts at each occurrence of the pattern
    for ([_]?[]const u8{ null, ccsdsSync }) |sync| {
        var parser = try astroz.Parser(Ccsds).init(null, null, 1024, init.io, gpa);
        defer parser.deinit();
        try parser.parseFromFile("test/ccsds.bin", sync, null);
        std.log.info("test/ccsds.bin, {s} sync pattern: {d} packets", .{ if (sync == null) "no" else "a", parser.packets.items.len });
        for (parser.packets.items, 0..) |p, k| {
            std.log.info("  {d}: {x}", .{ k, p.rawData[0 .. 5 + @as(usize, p.header.packetSize)] });
        }
    }

    var vita = try astroz.Parser(Vita49).init(null, null, 1024, init.io, gpa);
    defer vita.deinit();
    std.log.info("test/vita49.bin:", .{});
    try vita.parseFromFile("test/vita49.bin", &vitaSync, onVita49);
}
