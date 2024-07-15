const std = @import("std");
const astroz = @import("astroz");
const Vita49 = astroz.Vita49;
const Parser = astroz.Parser;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const P = Parser(Vita49);
    const ip = "127.0.0.1".*;
    const port: u16 = 65432;
    var parser = try P.init(&ip, port, 1024, allocator);
    defer parser.deinit();
    _ = try parser.start(callback);
}

fn callback(packet: Vita49) void {
    std.debug.print("Packet received: {any}", .{packet});
}
