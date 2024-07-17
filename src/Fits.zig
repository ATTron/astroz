const std = @import("std");
const zigimg = @import("zigimg");
const cfitsio = @import("cfitsio").c;

const Fits = @This();

allocator: std.mem.Allocator,

pub fn init(allocator: std.mem.Allocator) Fits {
    return .{
        .allocator = allocator,
    };
}

pub fn toImage(self: *Fits, input_path: []const u8, output_path: []const u8, options: StretchOptions) !void {
    var fptr: ?*cfitsio.fitsfile = null;
    var status: c_int = 0;
    _ = cfitsio.fits_open_file(&fptr, input_path.ptr, cfitsio.READONLY, &status);
    defer _ = cfitsio.fits_close_file(fptr, &status);

    if (status != 0) {
        return error.FITSOpenError;
    }

    var naxis: c_int = 0;
    var naxes: [2]c_long = undefined;
    _ = cfitsio.fits_get_img_dim(fptr, &naxis, &status);
    _ = cfitsio.fits_get_img_size(fptr, 2, &naxes, &status);

    const width: usize = @intCast(naxes[0]);
    const height: usize = @intCast(naxes[1]);

    const pixels = try self.allocator.alloc(f32, width * height);
    defer self.allocator.free(pixels);

    const fpixel: c_long = 1;
    const nelem: c_long = @intCast(width * height);
    _ = cfitsio.fits_read_img(fptr, cfitsio.TFLOAT, fpixel, nelem, null, pixels.ptr, null, &status);

    if (status != 0) {
        return error.FITSReadError;
    }

    try self.applyStretch(pixels, width, height, output_path, options);
}

fn applyStretch(self: *Fits, pixels: []f32, width: usize, height: usize, output_path: []const u8, options: StretchOptions) !void {
    const sorted_pixels = try self.allocator.dupe(f32, pixels);
    defer self.allocator.free(sorted_pixels);
    std.sort.heap(f32, sorted_pixels, {}, std.sort.asc(f32));

    const vmin_idx = sorted_pixels.len / 2000;
    const vmax_idx = sorted_pixels.len * 1995 / 2000;
    const vmin = sorted_pixels[vmin_idx];
    const vmax = sorted_pixels[vmax_idx];

    var image = try zigimg.Image.create(self.allocator, width, height, .rgb24);
    defer image.deinit();

    for (pixels, 0..) |pixel, i| {
        const normalized = std.math.clamp((pixel - vmin) / (vmax - vmin), 0, 1);
        const stretched = sineStretch(normalized, options);
        const color = applyColorMap(stretched);

        image.pixels.rgb24[i] = .{
            .r = @intFromFloat(color[0] * 255),
            .g = @intFromFloat(color[1] * 255),
            .b = @intFromFloat(color[2] * 255),
        };
    }

    try image.writeToFilePath(output_path, .{ .png = .{} });
}

fn sineStretch(x: f32, options: StretchOptions) f32 {
    const stretch: f32 = options.stretch;
    const bend: f32 = options.bend;
    return std.math.asinh((x - bend) / stretch) / std.math.asinh((1 - bend) / stretch) * 0.5 + 0.5;
}

fn applyColorMap(value: f32) [3]f32 {
    return .{ value, value, value };
}

pub const StretchOptions = struct {
    stretch: f32 = 0.15,
    bend: f32 = 0.5,
};

pub const FitsError = error{
    InvalidHeader,
    UnsupportedBitpix,
};

test Fits {
    var fits_png = Fits.init(std.testing.allocator);
    try fits_png.toImage("test/sample_fits.fits", "test/test.png", .{});
}
