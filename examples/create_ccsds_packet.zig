const std = @import("std");
const astroz = @import("astroz");

pub fn main(init: std.process.Init) !void {
    const bytes = [_]u8{ 0x78, 0x97, 0xC0, 0x00, 0x00, 0x0A, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A };

    var pkt = try astroz.Ccsds.init(&bytes, init.gpa, null);
    defer pkt.deinit();

    std.log.info("header: {any}", .{pkt.header});
    std.log.info("payload: {x}", .{pkt.packets});
}
