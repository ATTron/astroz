const std = @import("std");
const net = std.net;
const Vita49 = @import("vita49.zig").Vita49;
const CCSDS = @import("ccsds.zig").CCSDS;

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

        pub fn new(ip_address: []const u8, port: u16, buffer_size: u64, allocator: std.mem.Allocator) Self {
            return .{
                .ip_address = ip_address,
                .port = port,
                .buffer_size = buffer_size,
                .packets = std.ArrayList(Frame).init(allocator),
                .allocator = allocator,
            };
        }

        pub fn deinit(self: *Self) void {
            self.packets.deinit();
        }

        pub fn parse_from_file(self: Self, file_name: []const u8) void {
            const file = try std.fs.cwd().openFile(&file_name, .{});
            defer file.close();

            var buffer = std.io.bufferedReader(file.reader());
            const reader = buffer.reader();

            var line = std.ArrayList(u8).init(self.allocator);
            defer line.deinit();

            const writer = line.writer();
            var line_no: usize = 0;

            while (reader.streamUntilDelimiter(writer, "\n", null)) {
                defer line.clearRetainingCapacity();
                line_no += 1;

                std.log.warn("\n{s}\n", .{line.items});
            } else |err| switch (err) {
                error.EndOfStream => {
                    if (line.items.len > 0) {
                        line_no += 1;
                        std.log.warn("\n{s}\n", .{line.items});
                    }
                },
                else => return err,
            }

            std.log.debug("Total Lines Read: {d}\n", .{line_no});
        }

        pub fn start(self: *Self, comptime callback: ?fn (Frame) void) !void {
            const addr = try net.Address.parseIp4(self.ip_address, self.port);

            const stream = try net.tcpConnectToAddress(addr);
            defer stream.close();

            std.log.info("connected to socket successful", .{});

            var incoming_buffer = std.mem.zeroes([1024]u8);
            while (!self.should_stop) {
                _ = try stream.read(&incoming_buffer);
                const new_frame = try Frame.new(&incoming_buffer, null);
                std.log.debug("message recieved: {any}", .{new_frame});
                _ = try self.packets.append(new_frame);
                if (callback != null) {
                    callback.?(new_frame);
                }
            }
        }

        pub fn stop(self: *Self) void {
            self.should_stop = true;
        }
    };
}

/// this is for running the tests ONLY
fn _run_test_server(parse_type: []const u8) !void {
    const ip_addr = try net.Ip4Address.parse("127.0.0.1", 65432);
    const test_host = net.Address{ .in = ip_addr };
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

test "Vita49 Parser Test" {
    const ip = "127.0.0.1".*;
    const port: u16 = 65432;
    const parser = Parser(Vita49);
    var par_test = parser.new(&ip, port, 1024, std.testing.allocator);
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
    var par_test = parser.new(&ip, port, 1024, std.testing.allocator);
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
    var par_test = parser.new(&ip, port, 1024, std.testing.allocator);
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
