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

        const Self = @This();

        pub fn new(ip_address: []const u8, port: u16, buffer_size: u64) Self {
            return .{
                .ip_address = ip_address,
                .port = port,
                .buffer_size = buffer_size,
            };
        }

        pub fn parse_from_file(self: Self) void {
            self.is_running;
        }

        pub fn start(self: *Self) !void {
            std.log.warn("print something out and keep it please", .{});
            const addr = try net.Address.parseIp4(self.ip_address, self.port);
            // const stream_ip = net.Address{ .in = loopback };

            const stream = try net.tcpConnectToAddress(addr);
            defer stream.close();

            std.log.info("connected to socket successful", .{});

            var incoming_buffer = std.mem.zeroes([1024]u8);
            while (!self.should_stop) {
                _ = try stream.read(&incoming_buffer);
                const new_frame = Frame.new(&incoming_buffer);
                std.log.info("message recieved: {any}", .{new_frame});
            }
        }

        pub fn stop(self: *Self) void {
            self.should_stop = true;
        }
    };
}

/// this is for running the tests ONLY
fn _run_test_server() !void {
    const loopback = try net.Ip4Address.parse("127.0.0.1", 65432);
    const localhost = net.Address{ .in = loopback };
    var server = try localhost.listen(.{
        .reuse_port = true,
    });
    defer server.deinit();

    const addr = server.listen_address;
    std.log.info("Listening on {}\n", .{addr.getPort()});

    var client = try server.accept();
    defer client.stream.close();

    std.log.info("Connection received! {}\n", .{client.address});
    const vita49_pkt = [_]u8{
        0x3A, 0x02, 0x0A, 0x00, 0x34, 0x12, 0x00, 0x00, 0x00, 0x56, 0x34,
        0x12, 0x78, 0x9A, 0xBC, 0xDE, 0x80, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x48, 0x65, 0x6C, 0x6C, 0x6F, 0x2C, 0x20, 0x56, 0x49,
        0x54, 0x41, 0x20, 0x34, 0x39, 0x21,
    };
    var counter: usize = 0;
    while (counter < 5) {
        _ = try client.stream.writeAll(&vita49_pkt);
        std.time.sleep(2 * std.time.ns_per_s);
        counter += 1;
    }
}

test "Vita49 Parser Test" {
    const ip = "127.0.0.1".*;
    const port: u16 = 65432;
    const parser = Parser(Vita49);
    var par_test = parser.new(&ip, port, 1024);

    {
        const t1 = try std.Thread.spawn(.{}, _run_test_server, .{});
        defer t1.join();

        std.time.sleep(2 * std.time.ns_per_s);

        const t2 = try std.Thread.spawn(.{}, struct {
            fn run(pt: *parser) !void {
                try pt.start();
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
}
