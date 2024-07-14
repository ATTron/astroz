const std = @import("std");
const zigimg = @import("zigimg");

const WriteContext = struct {
    file: std.fs.File,
};

pub const FitsError = error{
    InvalidHeader,
    UnsupportedBitpix,
};

pub const FitsHeader = struct {
    bitpix: i32,
    naxis: u8,
    naxes: [2]u64,

    pub fn parse(file: *std.fs.File) !FitsHeader {
        var header: FitsHeader = undefined;
        var buffer: [80]u8 = undefined;

        while (true) {
            _ = try file.read(&buffer);
            const line = buffer[0..80];

            if (std.mem.eql(u8, line[0..8], "END     ")) break;

            if (std.mem.startsWith(u8, line, "BITPIX  =")) {
                header.bitpix = try std.fmt.parseInt(i32, std.mem.trim(u8, line[10..30], " "), 10);
            } else if (std.mem.startsWith(u8, line, "NAXIS   =")) {
                header.naxis = try std.fmt.parseInt(u8, std.mem.trim(u8, line[10..30], " "), 10);
            } else if (std.mem.startsWith(u8, line, "NAXIS1  =")) {
                header.naxes[0] = try std.fmt.parseInt(u64, std.mem.trim(u8, line[10..30], " "), 10);
            } else if (std.mem.startsWith(u8, line, "NAXIS2  =")) {
                header.naxes[1] = try std.fmt.parseInt(u64, std.mem.trim(u8, line[10..30], " "), 10);
            }
        }

        return header;
    }
};

pub const FitsFile = struct {
    header: FitsHeader,
    data: []f32,
    allocator: std.mem.Allocator,

    pub fn init(file_path: []const u8, allocator: std.mem.Allocator) !FitsFile {
        // var file = try std.fs.cwd().openFile(file_path, .{});
        std.log.warn("{s}", .{file_path});
        var file = try std.fs.cwd().openFile("/home/templeton/projects/zig/astroz/test/sample_fits.fits", .{});
        defer file.close();

        const header = try FitsHeader.parse(&file);
        const data = try readImageData(allocator, &file, header);
        return .{
            .header = header,
            .data = data,
            .allocator = allocator,
        };
    }

    fn printRawBytes(data: []const u8, count: usize) void {
        const to_print = @min(data.len, count);
        std.debug.print("First {} bytes: ", .{to_print});
        for (data[0..to_print]) |byte| {
            std.debug.print("{x:0>2} ", .{byte});
        }
        std.debug.print("\n", .{});
    }

    fn readImageData(allocator: std.mem.Allocator, file: *std.fs.File, header: FitsHeader) !struct { raw: []u8, floats: []f32 } {
        const pixel_count = header.naxes[0] * header.naxes[1];
        const bytes_per_pixel = @abs(header.bitpix) / 8;
        const total_bytes = pixel_count * bytes_per_pixel;

        var raw_data = try allocator.alloc(u8, total_bytes);
        errdefer allocator.free(raw_data);

        const bytes_read = try file.readAll(raw_data);
        if (bytes_read != total_bytes) {
            return error.UnexpectedEOF;
        }

        printRawBytes(raw_data, 32); // Print first 32 bytes

        var float_data = try allocator.alloc(f32, pixel_count);
        errdefer allocator.free(float_data);

        // Convert raw bytes to floats
        for (0..pixel_count) |i| {
            const start = i * 4;
            const value = std.mem.readInt(f32, raw_data[start..][0..4], .big);
            float_data[i] = value;
        }

        return .{ .raw = raw_data, .floats = float_data };
    }

    pub fn deinit(self: *FitsFile) void {
        self.allocator.free(self.data);
    }

    pub fn toImage(self: FitsFile, output_path: []const u8) !void {
        const width: u32 = @intCast(self.header.naxes[0]);
        const height: u32 = @intCast(self.header.naxes[1]);

        const stats = analyzeData(self.data);
        std.debug.print("Data stats: min={d}, max={d}, mean={d}\n", .{ stats.min, stats.max, stats.mean });

        var image = try zigimg.Image.create(self.allocator, width, height, .grayscale16);
        defer image.deinit();

        // Use mean and max for scaling
        const low = stats.mean;
        const high = stats.max;
        const range = high - low;

        for (self.data, 0..) |value, i| {
            var normalized = if (range != 0)
                (value - low) / range
            else
                0.5;

            // Clip values
            normalized = std.math.clamp(normalized, 0, 1);

            // Apply curve to enhance contrast
            normalized = std.math.pow(f32, normalized, 0.5);

            const gray_value = @as(u16, @intFromFloat(normalized * 65535.0));
            image.pixels.grayscale16[i] = .{ .value = gray_value };
        }

        try image.writeToFilePath(output_path, .{ .png = .{} });
    }

    fn analyzeData(data: []const f32) struct { min: f32, max: f32, mean: f32 } {
        var min: f32 = std.math.floatMax(f32);
        var max: f32 = std.math.floatMin(f32);
        var sum: f32 = 0;

        for (data) |value| {
            min = @min(min, value);
            max = @max(max, value);
            sum += value;
        }

        return .{
            .min = min,
            .max = max,
            .mean = sum / @as(f32, @floatFromInt(data.len)),
        };
    }
};

test "create sample image from fits" {
    var fits_image = try FitsFile.init("test/sample_fits.fits", std.testing.allocator);
    defer fits_image.deinit();
    std.log.warn("SHOWING FITS FILE INFO: {any}\n", .{fits_image.header});

    try fits_image.toImage("test.png");
}
