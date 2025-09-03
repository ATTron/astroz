const std = @import("std");
const Ccsds = @import("Ccsds.zig");
const Vita49 = @import("Vita49.zig");

/// Semi-generic function that takes either Vita49 or CCSDS as the Frame type.
pub fn Parser(comptime Frame: type) type {
    return struct {
        ipAddress: []const u8,
        port: u16,
        bufferSize: u64 = 1024,
        isRunning: bool = false,
        shouldStop: bool = false,
        packets: std.ArrayList(Frame),
        allocator: std.mem.Allocator,

        const Self = @This();

        pub fn init(
            ipAddress: ?[]const u8,
            port: ?u16,
            bufferSize: u64,
            allocator: std.mem.Allocator,
        ) !Self {
            return .{
                .ipAddress = ipAddress orelse "127.0.0.1",
                .port = port orelse 65432,
                .bufferSize = bufferSize,
                .packets = std.ArrayList(Frame){},
                .allocator = allocator,
            };
        }

        /// Clean up any allocator mess left per packet that is created by the Parser
        pub fn deinit(self: *Self) void {
            for (self.packets.items) |*packet| {
                packet.deinit();
            }
            self.packets.deinit(self.allocator);
        }

        /// Use this if you have a recording you need to parse
        pub fn parseFromFile(
            self: *Self,
            fileName: []const u8,
            syncPattern: ?[]const u8,
            callback: ?fn (Frame) void,
        ) !void {
            var fileContent = try std.fs.cwd().readFileAlloc(fileName, self.allocator, .unlimited);
            defer self.allocator.free(fileContent);

            std.log.debug("type of Frame is: {s}", .{@typeName(Frame)});
            if (std.mem.eql(u8, @typeName(Frame), "Vita49")) {
                if (syncPattern) |sp| {
                    var i: usize = 0;
                    while (fileContent.len > 4) : (i += 1) {
                        if (std.mem.startsWith(u8, fileContent[i..], sp)) {
                            const newFrame = try Frame.init(fileContent[i..], self.allocator, null);
                            try self.packets.append(self.allocator, newFrame);
                            if (callback) |cb| {
                                cb(newFrame);
                            }

                            const skipLength = newFrame.header.packetSize * 4;
                            if (skipLength > fileContent.len - i) {
                                break;
                            }

                            std.mem.copyForwards(u8, fileContent[0..], fileContent[i..]);

                            const newAllocSize = fileContent.len - i;
                            fileContent = self.allocator.realloc(fileContent, newAllocSize) catch |err| {
                                std.log.err("Failed to reallocate memory for parsing: {}", .{err});
                                break;
                            };

                            i = 0;
                        }
                    }
                } else {
                    const newFrame = try Frame.init(fileContent, self.allocator, null);
                    try self.packets.append(self.allocator, newFrame);
                    if (callback) |cb| {
                        cb(newFrame);
                    }

                    const skipLength = newFrame.header.packetSize * 4;

                    std.mem.copyForwards(u8, fileContent[0..], fileContent[skipLength - 1 ..]);

                    const newAllocSize = fileContent.len - 1;
                    fileContent = try self.allocator.realloc(fileContent, newAllocSize);
                }
            } else if (std.mem.eql(u8, @typeName(Frame), "Ccsds")) {
                if (syncPattern) |sp| {
                    var i: usize = 0;
                    while (fileContent.len > 4) : (i += 1) {
                        if (std.mem.startsWith(u8, fileContent[i..], sp)) {
                            const newFrame = try Frame.init(fileContent[i..], self.allocator, null);
                            try self.packets.append(self.allocator, newFrame);
                            if (callback) |cb| {
                                cb(newFrame);
                            }

                            const skipLength = newFrame.header.packetSize + 6;
                            if (skipLength > fileContent.len - i) {
                                break;
                            }

                            std.mem.copyForwards(u8, fileContent[0..], fileContent[i..]);

                            const newAllocSize = fileContent.len - i;
                            fileContent = self.allocator.realloc(fileContent, newAllocSize) catch |err| {
                                std.log.err("Failed to reallocate memory for parsing: {}", .{err});
                                break;
                            };

                            i = 0;
                        }
                    }
                } else {
                    const newFrame = try Frame.init(fileContent, self.allocator, null);
                    try self.packets.append(self.allocator, newFrame);
                    if (callback) |cb| {
                        cb(newFrame);
                    }

                    const skipLength = newFrame.header.packetSize + 6;

                    std.mem.copyForwards(u8, fileContent[0..], fileContent[skipLength - 1 ..]);

                    const newAllocSize = fileContent.len - 1;
                    fileContent = try self.allocator.realloc(fileContent, newAllocSize);
                }
            }
        }

        /// This will start the tcp listener and begin parsing as data comes in
        pub fn start(self: *Self, comptime callback: ?fn (Frame) void) !void {
            const addr = try std.net.Address.parseIp4(self.ipAddress, self.port);

            const stream = try std.net.tcpConnectToAddress(addr);
            defer stream.close();

            std.log.info("connected to socket successful", .{});

            var incomingBuffer = std.mem.zeroes([1024]u8);
            while (!self.shouldStop) {
                _ = try stream.read(&incomingBuffer);
                const newFrame = try Frame.init(&incomingBuffer, self.allocator, null);
                std.log.debug("message received: {any}", .{newFrame});
                _ = try self.packets.append(self.allocator, newFrame);
                if (callback != null) {
                    callback.?(newFrame);
                }
            }
        }

        /// Kills the tcp connection
        /// Make sure you clean up
        pub fn stop(self: *Self) void {
            self.shouldStop = true;
        }
    };
}

fn testRunServer(parse_type: []const u8) !void {
    const ipAddr = try std.net.Ip4Address.parse("127.0.0.1", 65432);
    const testHost = std.net.Address{ .in = ipAddr };
    var server = try testHost.listen(.{
        .reuse_address = true,
    });
    defer server.deinit();

    const addr = server.listen_address;
    std.log.info("Listening on {d}\n", .{addr.getPort()});

    var client = try server.accept();
    defer client.stream.close();

    std.log.info("Connection received! {any}\n", .{client.address});

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
        std.Thread.sleep(2 * std.time.ns_per_s);
        counter += 1;
    }
}

fn testCallback(packet: Vita49) void {
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
    const P = Parser(Ccsds);
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
    const P = Parser(Ccsds);
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
        const t1 = try std.Thread.spawn(.{}, testRunServer, .{"vita49"});
        defer t1.join();

        std.Thread.sleep(2 * std.time.ns_per_s);

        const t2 = try std.Thread.spawn(.{}, struct {
            fn run(pt: *parser) !void {
                try pt.start(null);
            }
        }.run, .{&par_test});
        defer t2.join();

        std.Thread.sleep(10 * std.time.ns_per_s);

        const t3 = try std.Thread.spawn(.{}, struct {
            fn run(pt: *parser) void {
                pt.stop();
            }
        }.run, .{&par_test});
        defer t3.join();
    }
    try std.testing.expectEqual(7, par_test.packets.capacity);
}

test "Vita49 Parser Test w/ Callback" {
    const ip = "127.0.0.1".*;
    const port: u16 = 65432;
    const parser = Parser(Vita49);
    var par_test = try parser.init(&ip, port, 1024, std.testing.allocator);
    defer par_test.deinit();

    {
        const t1 = try std.Thread.spawn(.{}, testRunServer, .{"vita49"});
        defer t1.join();

        std.Thread.sleep(2 * std.time.ns_per_s);

        const t2 = try std.Thread.spawn(.{}, struct {
            fn run(pt: *parser) !void {
                try pt.start(testCallback);
            }
        }.run, .{&par_test});
        defer t2.join();

        std.Thread.sleep(10 * std.time.ns_per_s);

        const t3 = try std.Thread.spawn(.{}, struct {
            fn run(pt: *parser) void {
                pt.stop();
            }
        }.run, .{&par_test});
        defer t3.join();
    }
    try std.testing.expectEqual(7, par_test.packets.capacity);
}

test "CCSDS Parser Test" {
    const ip = "127.0.0.1".*;
    const port: u16 = 65432;
    const parser = Parser(Ccsds);
    var par_test = try parser.init(&ip, port, 1024, std.testing.allocator);
    defer par_test.deinit();

    {
        const t1 = try std.Thread.spawn(.{}, testRunServer, .{"ccsds"});
        defer t1.join();

        std.Thread.sleep(2 * std.time.ns_per_s);

        const t2 = try std.Thread.spawn(.{}, struct {
            fn run(pt: *parser) !void {
                try pt.start(null);
            }
        }.run, .{&par_test});
        defer t2.join();

        std.Thread.sleep(10 * std.time.ns_per_s);

        const t3 = try std.Thread.spawn(.{}, struct {
            fn run(pt: *parser) void {
                pt.stop();
            }
        }.run, .{&par_test});
        defer t3.join();
    }
    try std.testing.expectEqual(7, par_test.packets.capacity);
}
