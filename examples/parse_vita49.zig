const std = @import("std");
const astroz = @import("astroz");
const Vita49 = astroz.vita49.Vita49;
const Parser = astroz.parsers.Parser;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const P = Parser(Vita49);
    const ip = "127.0.0.1".*;
    const port: u16 = 65432;
    var parser = try P.init(&ip, port, 1024, allocator);
    defer parser.deinit();
    _ = try parser.start(null);
}
