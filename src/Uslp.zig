const std = @import("std");

pub const Uslp = @This();

primaryHeader: USLPHeader,
insertZone: []const u8,
dataHeader: DataHeader,
TFDZ: []const u8,
OpControl: u32,
FrameError: u16,
rawData: []const u8,
allocator: std.mem.Allocator,

pub fn init(pl: []const u8, allocator: std.mem.Allocator) !Uslp{
    var rawPackets = try allocator.dupe(u8, pl); // dupe so it doesn't go out of scope too early
    const primaryHeader = rawPackets[1..6];
    const headers = fetchHeader(rawPackets);
    const end = 5 + headers.header.packetSize; // num of header bytes + packet_size
    const packets = rawPackets[headers.start..end];

    _ = allocator.resize(rawPackets, end);

    return .{
        .header = headers.header,
        .allocator = allocator,
        .rawData = rawPackets,
    };
}

pub fn deinit(self: *Uslp) !void {
    self.allocator.free(self.rawData);
}

pub fn fetchHeader(rawPackets: u8) struct {header: USLPHeader } {
    const header = rawPackets[0..6];
    const tfv = @as(u4, @truncate((header[0] >> 5) & 0x07));
    std.debug.print("%s", .{tfv});
    // const tfvn = rawPackets[0..3];
    // const scid = rawPackets[3..19];
}

pub const DataHeader = packed struct {};

pub const USLPHeader = packed struct {
    xferFrameVersion: u4,
    scID: u16,
    srcID: u1,
    virtChannel: u6,
    mapID: u4,
    eofFlag: u1,
    frameLength: u16,
    bypassCntrl: u1,
    protocolCmd: u1,
    spare: u2,
    ocf: u1,
    VCFrameCountLength: u3,
    VCFrameCount: u56,
};


test "CCSDS Structure Testing w/o config" {
    const raw_test_packet: [16]u8 = .{ 0x78, 0x97, 0xC0, 0x00, 0x00, 0x0A, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A };
    var converted_test_packet = try Uslp.init(&raw_test_packet, std.testing.allocator, null);
    defer converted_test_packet.deinit();

    try std.testing.expectEqual(3, converted_test_packet.header.version);

}
