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
        io: std.Io,
        allocator: std.mem.Allocator,

        const Self = @This();

        pub fn init(
            ipAddress: ?[]const u8,
            port: ?u16,
            bufferSize: u64,
            io: std.Io,
            allocator: std.mem.Allocator,
        ) !Self {
            return .{
                .ipAddress = ipAddress orelse "127.0.0.1",
                .port = port orelse 65432,
                .bufferSize = bufferSize,
                .packets = std.ArrayList(Frame){},
                .io = io,
                .allocator = allocator,
            };
        }

        /// Clean up any allocator mess left per packet that is created by parser
        pub fn deinit(self: *Self) void {
            for (self.packets.items) |*packet| {
                packet.deinit();
            }
            self.packets.deinit(self.allocator);
        }

        /// Compute total packet byte size from a parsed frame
        fn frameByteSize(frame: Frame) usize {
            if (comptime std.mem.eql(u8, @typeName(Frame), "Vita49"))
                return @as(usize, frame.header.packetSize) * 4
            else
                return @as(usize, frame.header.packetSize) + 6;
        }

        /// Use this if you have a recording you need to parse
        pub fn parseFromFile(
            self: *Self,
            fileName: []const u8,
            syncPattern: ?[]const u8,
            callback: ?fn (Frame) void,
        ) !void {
            var fileContent = try std.Io.Dir.cwd().readFileAlloc(self.io, fileName, self.allocator, .unlimited);
            defer self.allocator.free(fileContent);

            if (syncPattern) |sp| {
                var i: usize = 0;
                while (fileContent.len > 4) : (i += 1) {
                    if (std.mem.startsWith(u8, fileContent[i..], sp)) {
                        const newFrame = try Frame.init(fileContent[i..], self.allocator, null);
                        try self.packets.append(self.allocator, newFrame);
                        if (callback) |cb| cb(newFrame);

                        if (frameByteSize(newFrame) > fileContent.len - i) break;

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
                if (callback) |cb| cb(newFrame);

                const skipLength = frameByteSize(newFrame);
                std.mem.copyForwards(u8, fileContent[0..], fileContent[skipLength - 1 ..]);
                const newAllocSize = fileContent.len - 1;
                fileContent = try self.allocator.realloc(fileContent, newAllocSize);
            }
        }

        /// This will start the tcp listener and begin parsing as data comes in
        pub fn start(self: *Self, io: std.Io, comptime callback: ?fn (Frame) void) !void {
            const addr = try std.Io.net.IpAddress.parseIp4(self.ipAddress, self.port);

            const stream = try std.Io.net.IpAddress.connect(addr, io, .{
                .mode = .raw,
                .protocol = .tcp,
                .timeout = .{ .duration = .{
                    .raw = std.Io.Duration.fromNanoseconds(500 * std.time.ns_per_ms),
                    .clock = .real,
                } },
            });
            defer stream.close(io);

            std.log.info("connected to socket successful", .{});

            var incomingBuffer = std.mem.zeroes([1024]u8);
            var readerBuffer = std.mem.zeroes([1024]u8);
            var reader = stream.reader(io, &readerBuffer);
            while (!self.shouldStop) {
                const bytesRead = try reader.interface.readSliceShort(&incomingBuffer);
                if (bytesRead == 0) continue;

                const newFrame = Frame.init(incomingBuffer[0..bytesRead], self.allocator, null) catch continue;
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

test "Vita49 Parse From File w/ sync" {
    const io = std.Io.Threaded.global_single_threaded.ioBasic();
    const file_name = "./test/vita49.bin".*;
    //3a02 0a00 3412 0000 0056
    const sync_pattern = .{ 0x3A, 0x02, 0x0a, 0x00, 0x34, 0x12, 0x00, 0x00, 0x00, 0x56 };
    const P = Parser(Vita49);
    var parser = try P.init(null, null, 1024, io, std.testing.allocator);
    defer parser.deinit();

    _ = try parser.parseFromFile(&file_name, &sync_pattern, null);
    for (parser.packets.items) |packet| {
        try std.testing.expectEqualStrings("Hello, VITA 49!", packet.payload);
    }
}

test "Vita49 Parse From File w/o sync" {
    const io = std.Io.Threaded.global_single_threaded.ioBasic();
    const file_name = "./test/vita49.bin".*;
    //3a02 0a00 3412 0000 0056
    const P = Parser(Vita49);
    var parser = try P.init(null, null, 1024, io, std.testing.allocator);
    defer parser.deinit();

    _ = try parser.parseFromFile(&file_name, null, null);
    for (parser.packets.items) |packet| {
        try std.testing.expectEqualStrings("Hello, VITA 49!", packet.payload);
    }
}

test "CCSDS Parse From File w/o sync" {
    const io = std.Io.Threaded.global_single_threaded.ioBasic();
    const file_name = "./test/ccsds.bin".*;
    const P = Parser(Ccsds);
    var parser = try P.init(null, null, 1024, io, std.testing.allocator);
    defer parser.deinit();

    const packets = .{ 5, 6, 7, 8, 9, 10 };

    _ = try parser.parseFromFile(&file_name, null, null);
    for (parser.packets.items) |packet| {
        try std.testing.expectEqualSlices(u8, &packets, packet.packets);
    }
}

test "CCSDS Parse From File w/ sync" {
    const io = std.Io.Threaded.global_single_threaded.ioBasic();
    const file_name = "./test/ccsds.bin".*;
    // 7897 c000 000a 0102
    const sync_pattern = .{ 0x78, 0x97, 0xC0, 0x00, 0x00, 0x0A, 0x01, 0x02 };
    const P = Parser(Ccsds);
    var parser = try P.init(null, null, 1024, io, std.testing.allocator);
    defer parser.deinit();

    const packets = .{ 5, 6, 7, 8, 9, 10 };

    _ = try parser.parseFromFile(&file_name, &sync_pattern, null);
    for (parser.packets.items) |packet| {
        try std.testing.expectEqualSlices(u8, &packets, packet.packets);
    }
}
