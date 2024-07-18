const std = @import("std");
const zigimg = @import("zigimg");
const cfitsio = @import("cfitsio").c;

const Fits = @This();

allocator: std.mem.Allocator,
fptr: ?*cfitsio.fitsfile,

pub fn open(file_path: []const u8, allocator: std.mem.Allocator) !Fits {
    std.log.debug("{s}", .{file_path});
    var status: c_int = 0;
    var fptr: ?*cfitsio.fitsfile = null;
    _ = cfitsio.fits_open_file(&fptr, file_path.ptr, cfitsio.READONLY, &status);
    if (status != 0) return FitsError.OpenError;

    return Fits{
        .allocator = allocator,
        .fptr = fptr,
    };
}

pub fn close(self: *Fits) void {
    if (self.fptr) |fptr| {
        var status: c_int = 0;
        _ = cfitsio.fits_close_file(fptr, &status);
        self.fptr = null;
    }
}

pub fn readTable(self: *Fits, hdu: c_int) !void {
    var status: c_int = 0;
    var hdu_type: c_int = undefined;

    _ = cfitsio.fits_movabs_hdu(self.fptr, hdu, &hdu_type, &status);
    if (status != 0) {
        try self.reportError(status);
        return FitsError.ReadError;
    }

    if (hdu_type != cfitsio.BINARY_TBL) {
        std.debug.print("HDU {d} is not a binary table. Type: {d}\n", .{ hdu, hdu_type });
        return FitsError.NotATableError;
    }

    var rows: c_long = undefined;
    var cols: c_int = undefined;
    _ = cfitsio.fits_get_num_rows(self.fptr, &rows, &status);
    _ = cfitsio.fits_get_num_cols(self.fptr, &cols, &status);
    if (status != 0) {
        try self.reportError(status);
        return FitsError.ReadError;
    }

    std.debug.print("Table has {d} rows and {d} columns\n", .{ rows, cols });

    var i: c_int = 1;
    while (i <= cols) : (i += 1) {
        var ttype: [cfitsio.FLEN_VALUE]u8 = undefined;
        var tunit: [cfitsio.FLEN_VALUE]u8 = undefined;
        var typechar: [cfitsio.FLEN_VALUE]u8 = undefined;
        var repeat: c_long = undefined;
        var scale: f64 = undefined;
        var zero: f64 = undefined;
        var nulval: c_long = undefined;
        var tdisp: [cfitsio.FLEN_VALUE]u8 = undefined;

        status = 0;
        _ = cfitsio.fits_get_bcolparms(self.fptr, i, &ttype, &tunit, &typechar, &repeat, &scale, &zero, &nulval, &tdisp, &status);
        if (status != 0) {
            try self.reportError(status);
            return FitsError.ReadError;
        }

        const col_name = std.mem.sliceTo(&ttype, 0);
        const type_str = std.mem.sliceTo(&typechar, 0);
        std.debug.print("Column {d}: name = {s}, type = {s}, repeat = {d}\n", .{ i, col_name, type_str, repeat });
    }
    std.debug.print("Finished processing all columns.\n", .{});
}

fn reportError(self: *Fits, status: c_int) !void {
    _ = self;
    var error_text: [128]u8 = undefined;
    _ = cfitsio.fits_get_errstatus(status, &error_text);
    std.debug.print("CFITSIO Error: {s}\n", .{std.mem.sliceTo(&error_text, 0)});
}

pub fn readImage(self: *Fits, output_path: []const u8, options: StretchOptions) !void {
    const fptr = self.fptr;
    var status: c_int = 0;
    // _ = cfitsio.fits_open_file(&fptr, input_path.ptr, cfitsio.READONLY, &status);
    // defer _ = cfitsio.fits_close_file(fptr, &status);

    // if (status != 0) {
    //     return error.FITSOpenError;
    // }

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

inline fn sineStretch(x: f32, options: StretchOptions) f32 {
    const stretch: f32 = options.stretch;
    const bend: f32 = options.bend;
    return std.math.asinh((x - bend) / stretch) / std.math.asinh((1 - bend) / stretch) * 0.5 + 0.5;
}

inline fn applyColorMap(value: f32) [3]f32 {
    return .{ value, value, value };
}

pub const StretchOptions = struct {
    stretch: f32 = 0.15,
    bend: f32 = 0.5,
};

pub const FitsError = error{
    OpenError,
    ReadError,
    MemoryError,
    UnknownError,
    NotATableError,
};

pub const FitsTableInfo = struct {
    rows: f32,
    columns: f32,
};

test Fits {
    var fits_png = try Fits.open("test/sample_fits.fits", std.testing.allocator);
    defer fits_png.close();
    try fits_png.readImage("test/test.png", .{});
}

test "Read FITS table" {
    std.debug.print("\n--- Starting FITS table test ---\n", .{});

    var fits_file = try Fits.open("./test/table.fits", std.testing.allocator);
    defer fits_file.close();

    std.debug.print("Successfully opened FITS file\n", .{});

    // Try to read the first HDU
    fits_file.readTable(1) catch |err| {
        std.debug.print("Expected error reading HDU 1: {}\n", .{err});
    };

    std.debug.print("Attempting to read HDU 2...\n", .{});
    try fits_file.readTable(2);
    std.debug.print("Successfully read HDU 2\n", .{});
}
