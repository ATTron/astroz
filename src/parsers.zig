const std = @import("std");
const CCSDS = @import("ccsds.zig").CCSDS;
const Vita49 = @import("vita49.zig").Vita49;

/// Semi-generic function that takes either Vita49 or CCSDS as the Frame type.
pub fn Parser(comptime Frame: type) type {
    return struct {
        ip_address: []const u8,
        port: u16,
        buffer_size: u64 = 1024,
        is_running: bool = false,
        should_stop: bool = false,
        packets: std.ArrayList(Frame),
        allocator: std.mem.Allocator,

        const Self = @This();

        pub fn init(
            ip_address: ?[]const u8,
            port: ?u16,
            buffer_size: u64,
            allocator: std.mem.Allocator,
        ) !Self {
            return .{
                .ip_address = ip_address orelse "127.0.0.1",
                .port = port orelse 65432,
                .buffer_size = buffer_size,
                .packets = std.ArrayList(Frame).init(allocator),
                .allocator = allocator,
            };
        }

        /// If you are using a parser, this will clean up any allocator mess left per packet that is created by the Parser
        pub fn deinit(self: *Self) void {
            for (self.packets.items) |*packet| {
                packet.deinit();
            }
            self.packets.deinit();
        }

        /// Use this if you have a recording you need to parse
        pub fn parseFromFile(
            self: *Self,
            file_name: []const u8,
            sync_pattern: ?[]const u8,
            callback: ?fn (Frame) void,
        ) !void {
            const file = try std.fs.cwd().openFile(file_name, .{});
            defer file.close();

            const stat = try file.stat();
            var file_content = try file.readToEndAlloc(
                self.allocator,
                stat.size,
            );
            defer self.allocator.free(file_content);

            std.log.debug("type of Frame is: {s}", .{@typeName(Frame)});
            if (std.mem.eql(u8, @typeName(Frame), "vita49.Vita49")) {
                if (sync_pattern) |sp| {
                    var i: usize = 0;
                    while (file_content.len > 4) : (i += 1) {
                        if (std.mem.eql(u8, file_content[i .. i + sp.len], sp)) {
                            const new_frame = try Frame.init(file_content[i..], self.allocator, null);
                            try self.packets.append(new_frame);
                            if (callback) |cb| {
                                cb(new_frame);
                            }

                            const skip_lengh = new_frame.header.packet_size * 4;
                            if (skip_lengh > file_content.len - i) {
                                break;
                            }

                            std.mem.copyForwards(u8, file_content[0..], file_content[i..]);

                            const new_alloc_size = file_content.len - i;
                            file_content = self.allocator.realloc(file_content, new_alloc_size) catch |err| {
                                std.log.err("Failed to reallocate memory for parsing: {}", .{err});
                                break;
                            };

                            i = 0;
                        }
                    }
                } else {
                    const new_frame = try Frame.init(file_content, self.allocator, null);
                    try self.packets.append(new_frame);
                    if (callback) |cb| {
                        cb(new_frame);
                    }

                    const skip_lengh = new_frame.header.packet_size * 4;

                    std.mem.copyForwards(u8, file_content[0..], file_content[skip_lengh - 1 ..]);

                    const new_alloc_size = file_content.len - 1;
                    file_content = try self.allocator.realloc(file_content, new_alloc_size);
                }
            } else if (std.mem.eql(u8, @typeName(Frame), "ccsds.CCSDS")) {
                if (sync_pattern) |sp| {
                    var i: usize = 0;
                    while (file_content.len > 4) : (i += 1) {
                        if (std.mem.eql(u8, file_content[i .. i + sp.len], sp)) {
                            const new_frame = try Frame.init(file_content[i..], self.allocator, null);
                            try self.packets.append(new_frame);
                            if (callback) |cb| {
                                cb(new_frame);
                            }

                            const skip_lengh = new_frame.header.packet_size + 6;
                            if (skip_lengh > file_content.len - i) {
                                break;
                            }

                            std.mem.copyForwards(u8, file_content[0..], file_content[i..]);

                            const new_alloc_size = file_content.len - i;
                            file_content = self.allocator.realloc(file_content, new_alloc_size) catch |err| {
                                std.log.err("Failed to reallocate memory for parsing: {}", .{err});
                                break;
                            };

                            i = 0;
                        }
                    }
                } else {
                    const new_frame = try Frame.init(file_content, self.allocator, null);
                    try self.packets.append(new_frame);
                    if (callback) |cb| {
                        cb(new_frame);
                    }

                    const skip_lengh = new_frame.header.packet_size + 6;

                    std.mem.copyForwards(u8, file_content[0..], file_content[skip_lengh - 1 ..]);

                    const new_alloc_size = file_content.len - 1;
                    file_content = try self.allocator.realloc(file_content, new_alloc_size);
                }
            }
        }

        /// This will start the tcp listener and begin parsing as data comes in
        pub fn start(self: *Self, comptime callback: ?fn (Frame) void) !void {
            const addr = try std.net.Address.parseIp4(self.ip_address, self.port);

            const stream = try std.net.tcpConnectToAddress(addr);
            defer stream.close();

            std.log.info("connected to socket successful", .{});

            var incoming_buffer = std.mem.zeroes([1024]u8);
            while (!self.should_stop) {
                _ = try stream.read(&incoming_buffer);
                const new_frame = try Frame.init(&incoming_buffer, self.allocator, null);
                std.log.debug("message recieved: {any}", .{new_frame});
                _ = try self.packets.append(new_frame);
                if (callback != null) {
                    callback.?(new_frame);
                }
            }
        }

        /// Kills the tcp connection
        /// Make sure you clean up
        pub fn stop(self: *Self) void {
            self.should_stop = true;
        }
    };
}

/// this is for running the tests ONLY
fn _run_test_server(parse_type: []const u8) !void {
    const ip_addr = try std.net.Ip4Address.parse("127.0.0.1", 65432);
    const test_host = std.net.Address{ .in = ip_addr };
    var server = try test_host.listen(.{
        .reuse_port = true,
    });
    defer server.deinit();

    const addr = server.listen_address;
    std.log.info("Listening on {}\n", .{addr.getPort()});

    var client = try server.accept();
    defer client.stream.close();

    std.log.info("Connection received! {}\n", .{client.address});

    var _pkt: []const u8 = undefined;
    if (std.mem.eql(u8, parse_type, "vita49")) {
        _pkt = &[_]u8{
            0x3A, 0x02, 0x0A, 0x00, 0x34, 0x12, 0x00, 0x00, 0x00, 0x56, 0x34,
            0x12, 0x78, 0x9A, 0xBC, 0xDE, 0x80, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00, 0x48, 0x65, 0x6C, 0x6C, 0x6F, 0x2C, 0x20, 0x56, 0x49,
            0x54, 0x41, 0x20, 0x34, 0x39, 0x21,
        };
    } else {
        _pkt = &[_]u8{
            0x78, 0x97, 0xC0, 0x00, 0x00, 0x0A, 0x01, 0x02,
            0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A,
        };
    }
    var counter: usize = 0;
    while (counter < 5) {
        _ = try client.stream.writeAll(_pkt);
        std.time.sleep(2 * std.time.ns_per_s);
        counter += 1;
    }
}

/// this is for running tests ONLY
fn _test_callback(packet: Vita49) void {
    std.log.debug("CALLBACK CALLED: {any}", .{packet});
}

test "Vita49 Parse From File w/ sync" {
    const file_name = "./test/vita49.bin".*;
    //3a02 0a00 3412 0000 0056
    const sync_pattern = .{ 0x3A, 0x02, 0x0a, 0x00, 0x34, 0x12, 0x00, 0x00, 0x00, 0x56 };
    const P = Parser(Vita49);
    var parser = try P.init(null, null, 1024, std.testing.allocator);
    defer parser.deinit();

    _ = try parser.parseFromFile(&file_name, &sync_pattern, null);
    for (parser.packets.items) |packet| {
        try std.testing.expectEqualStrings("Hello, VITA 49!", packet.payload);
    }
}

test "Vita49 Parse From File w/o sync" {
    const file_name = "./test/vita49.bin".*;
    //3a02 0a00 3412 0000 0056
    const P = Parser(Vita49);
    var parser = try P.init(null, null, 1024, std.testing.allocator);
    defer parser.deinit();

    _ = try parser.parseFromFile(&file_name, null, null);
    for (parser.packets.items) |packet| {
        try std.testing.expectEqualStrings("Hello, VITA 49!", packet.payload);
    }
}

test "CCSDS Parse From File w/o sync" {
    const file_name = "./test/ccsds.bin".*;
    const P = Parser(CCSDS);
    var parser = try P.init(null, null, 1024, std.testing.allocator);
    defer parser.deinit();

    const packets = .{ 5, 6, 7, 8, 9, 10 };

    _ = try parser.parseFromFile(&file_name, null, null);
    for (parser.packets.items) |packet| {
        try std.testing.expectEqualSlices(u8, &packets, packet.packets);
    }
}

test "CCSDS Parse From File w/ sync" {
    const file_name = "./test/ccsds.bin".*;
    // 7897 c000 000a 0102
    const sync_pattern = .{ 0x78, 0x97, 0xC0, 0x00, 0x00, 0x0A, 0x01, 0x02 };
    const P = Parser(CCSDS);
    var parser = try P.init(null, null, 1024, std.testing.allocator);
    defer parser.deinit();

    const packets = .{ 5, 6, 7, 8, 9, 10 };

    _ = try parser.parseFromFile(&file_name, &sync_pattern, null);
    for (parser.packets.items) |packet| {
        try std.testing.expectEqualSlices(u8, &packets, packet.packets);
    }
}

test "Vita49 Parser Test" {
    const ip = "127.0.0.1".*;
    const port: u16 = 65432;
    const parser = Parser(Vita49);
    var par_test = try parser.init(&ip, port, 1024, std.testing.allocator);
    defer par_test.deinit();

    {
        const t1 = try std.Thread.spawn(.{}, _run_test_server, .{"vita49"});
        defer t1.join();

        std.time.sleep(2 * std.time.ns_per_s);

        const t2 = try std.Thread.spawn(.{}, struct {
            fn run(pt: *parser) !void {
                try pt.start(null);
            }
        }.run, .{&par_test});
        defer t2.join();

        std.time.sleep(10 * std.time.ns_per_s);

        const t3 = try std.Thread.spawn(.{}, struct {
            fn run(pt: *parser) void {
                pt.stop();
            }
        }.run, .{&par_test});
        defer t3.join();
    }
    try std.testing.expectEqual(8, par_test.packets.capacity);
}

test "Vita49 Parser Test w/ Callback" {
    const ip = "127.0.0.1".*;
    const port: u16 = 65432;
    const parser = Parser(Vita49);
    var par_test = try parser.init(&ip, port, 1024, std.testing.allocator);
    defer par_test.deinit();

    {
        const t1 = try std.Thread.spawn(.{}, _run_test_server, .{"vita49"});
        defer t1.join();

        std.time.sleep(2 * std.time.ns_per_s);

        const t2 = try std.Thread.spawn(.{}, struct {
            fn run(pt: *parser) !void {
                try pt.start(_test_callback);
            }
        }.run, .{&par_test});
        defer t2.join();

        std.time.sleep(10 * std.time.ns_per_s);

        const t3 = try std.Thread.spawn(.{}, struct {
            fn run(pt: *parser) void {
                pt.stop();
            }
        }.run, .{&par_test});
        defer t3.join();
    }
    try std.testing.expectEqual(8, par_test.packets.capacity);
}

test "CCSDS Parser Test" {
    const ip = "127.0.0.1".*;
    const port: u16 = 65432;
    const parser = Parser(CCSDS);
    var par_test = try parser.init(&ip, port, 1024, std.testing.allocator);
    defer par_test.deinit();

    {
        const t1 = try std.Thread.spawn(.{}, _run_test_server, .{"ccsds"});
        defer t1.join();

        std.time.sleep(2 * std.time.ns_per_s);

        const t2 = try std.Thread.spawn(.{}, struct {
            fn run(pt: *parser) !void {
                try pt.start(null);
            }
        }.run, .{&par_test});
        defer t2.join();

        std.time.sleep(10 * std.time.ns_per_s);

        const t3 = try std.Thread.spawn(.{}, struct {
            fn run(pt: *parser) void {
                pt.stop();
            }
        }.run, .{&par_test});
        defer t3.join();
    }
    try std.testing.expectEqual(8, par_test.packets.capacity);
}
