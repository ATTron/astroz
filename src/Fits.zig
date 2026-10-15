const std = @import("std");

const cfitsio = @import("cfitsio").c;
const zignal = @import("zignal");

const Rgba = zignal.Rgba(u8);
const Fits = @This();

allocator: std.mem.Allocator,
io: std.Io,
fptr: ?*cfitsio.fitsfile,
file_path: []const u8,
hdus: []HduInfo,

pub fn open(file_path: []const u8, io: std.Io, allocator: std.mem.Allocator) !Fits {
    var status: c_int = 0;
    var fptr: ?*cfitsio.fitsfile = null;
    _ = cfitsio.fits_open_file(&fptr, file_path.ptr, cfitsio.READONLY, &status);
    if (status != 0) return FitsError.OpenError;

    var num_hdus: c_int = undefined;
    _ = cfitsio.fits_get_num_hdus(fptr, &num_hdus, &status);
    if (status != 0) return FitsError.HeaderReadError;

    const hdus = try allocator.alloc(HduInfo, @intCast(num_hdus));
    errdefer allocator.free(hdus);

    for (hdus, 0..) |*hdu, i| {
        const hdu_num: c_int = @intCast(i + 1);
        _ = cfitsio.fits_movabs_hdu(fptr, hdu_num, null, &status);
        if (status != 0) return FitsError.HeaderReadError;

        var hdu_type: c_int = undefined;
        _ = cfitsio.fits_get_hdu_type(fptr, &hdu_type, &status);
        if (status != 0) return FitsError.HeaderReadError;

        hdu.* = .{
            .number = hdu_num,
            .content_type = switch (hdu_type) {
                cfitsio.IMAGE_HDU => .Image,
                cfitsio.ASCII_TBL => .AsciiTable,
                cfitsio.BINARY_TBL => .BinaryTable,
                else => .Unknown,
            },
        };
    }

    return Fits{
        .allocator = allocator,
        .fptr = fptr,
        .io = io,
        .file_path = file_path,
        .hdus = hdus,
    };
}

pub fn open_and_parse(file_path: []const u8, io: std.Io, allocator: std.mem.Allocator, options: ParseOptions) !Fits {
    var fits_file = try Fits.open(file_path, io, allocator);
    var image_count: usize = 0;
    var ascii_table_count: usize = 0;
    var binary_table_count: usize = 0;

    var output_path_buffer: [std.fs.max_path_bytes]u8 = undefined;
    std.log.debug("Processing HDUs", .{});
    for (fits_file.hdus) |hdu| {
        switch (hdu.content_type) {
            .Image => {
                image_count += 1;
                if (options.createImages) {
                    const filename = try std.fmt.bufPrint(&output_path_buffer, "image_{d}", .{image_count});
                    try fits_file.readImage(hdu.number, filename, options.stretchOptions);
                } else {
                    const filename = try std.fmt.bufPrint(&output_path_buffer, "image_{d}_data", .{image_count});
                    try fits_file.readImageAsTable(hdu.number, filename);
                }
            },
            .AsciiTable => {
                ascii_table_count += 1;
                const filename = try std.fmt.bufPrint(&output_path_buffer, "ascii_table_{d}", .{ascii_table_count});
                std.log.debug("Reading ASCII table", .{});
                try fits_file.readTable(hdu.number, filename);
            },
            .BinaryTable => {
                binary_table_count += 1;
                const filename = try std.fmt.bufPrint(&output_path_buffer, "binary_table_{d}", .{binary_table_count});
                std.log.debug("Reading binary table", .{});
                try fits_file.readTable(hdu.number, filename);
            },
            .Unknown => {
                std.log.warn("Skipping unknown HDU type", .{});
            },
        }
    }
    return fits_file;
}

pub fn close(self: *Fits) void {
    if (self.fptr) |fptr| {
        var status: c_int = 0;
        _ = cfitsio.fits_close_file(fptr, &status);
        self.fptr = null;
    }
    self.allocator.free(self.hdus);
}

pub fn readTable(self: *Fits, hdu_number: c_int, output_path: ?[]const u8) !void {
    if (hdu_number < 1 or hdu_number > self.hdus.len) {
        return FitsError.InvalidHDUError;
    }
    switch (self.hdus[@intCast(hdu_number - 1)].content_type) {
        .AsciiTable, .BinaryTable => {},
        else => return FitsError.NotATableError,
    }
    var status: c_int = 0;
    _ = cfitsio.fits_movabs_hdu(self.fptr, hdu_number, null, &status);
    try check(status);
    var rows: c_long = undefined;
    var cols: c_int = undefined;
    _ = cfitsio.fits_get_num_rows(self.fptr, &rows, &status);
    _ = cfitsio.fits_get_num_cols(self.fptr, &cols, &status);
    try check(status);
    std.log.debug("Table has {d} rows and {d} columns", .{ rows, cols });

    const dir = self.outputDir();
    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const csv_path = try std.fmt.bufPrint(&path_buf, "{s}/{s}.csv", .{ dir, output_path orelse dir });

    const csv_file = try std.Io.Dir.cwd().createFile(self.io, csv_path, .{});
    defer csv_file.close(self.io);
    var out_buf: [4096]u8 = undefined;
    var csv_writer = csv_file.writer(self.io, &out_buf);
    const out = &csv_writer.interface;

    var i: c_int = 1;
    while (i <= cols) : (i += 1) {
        var ttype: [cfitsio.FLEN_VALUE]u8 = undefined;
        var tunit: [cfitsio.FLEN_VALUE]u8 = undefined;
        var typechar: [1]u8 = undefined;
        var repeat: c_long = undefined;
        var scale: f64 = undefined;
        var zero: f64 = undefined;
        var nulval: c_long = undefined;
        var tdisp: [cfitsio.FLEN_VALUE]u8 = undefined;
        status = 0;
        _ = cfitsio.fits_get_bcolparms(self.fptr, i, &ttype, &tunit, &typechar, &repeat, &scale, &zero, &nulval, &tdisp, &status);
        try check(status);
        try out.print("{s}", .{std.mem.sliceTo(&ttype, 0)});
        if (i < cols) {
            try out.writeAll(",");
        }
    }
    try out.writeAll("\n");

    var row_buffer = try self.allocator.alloc(u8, @as(usize, @intCast(cols)) * cfitsio.FLEN_VALUE);
    defer self.allocator.free(row_buffer);

    var row_pointers = try self.allocator.alloc([*c]u8, @as(usize, @intCast(cols)));
    defer self.allocator.free(row_pointers);

    for (0..@as(usize, @intCast(cols))) |j| {
        row_pointers[j] = @ptrCast(&row_buffer[j * cfitsio.FLEN_VALUE]);
    }

    const row_ptr: [*c][*c]u8 = @ptrCast(row_pointers.ptr);

    var row: c_long = 1;
    while (row <= rows) : (row += 1) {
        i = 1;
        while (i <= cols) : (i += 1) {
            _ = cfitsio.fits_read_col_str(self.fptr, i, row, 1, 1, null, row_ptr, null, &status);
            try check(status);
            try out.print("{s}", .{std.mem.sliceTo(row_ptr[0], 0)});
            if (i < cols) {
                try out.writeAll(",");
            }
        }
        try out.writeAll("\n");
    }
    try out.flush();
    std.log.debug("CSV file created: {s}", .{csv_path});
}

/// Prints the CFITSIO message for a non-zero status and returns ReadError.
fn check(status: c_int) FitsError!void {
    if (status == 0) return;
    var error_text: [128]u8 = undefined;
    _ = cfitsio.fits_get_errstatus(status, &error_text);
    std.debug.print("CFITSIO Error: {s}\n", .{std.mem.sliceTo(&error_text, 0)});
    return FitsError.ReadError;
}

const Pixels = struct {
    data: []f32,
    width: usize,
    height: usize,
};

/// Reads a 2D image HDU as f32. Returns null if the image is empty.
fn readPixels(self: *Fits, hdu_number: c_int) !?Pixels {
    if (hdu_number < 1 or hdu_number > self.hdus.len) {
        return FitsError.InvalidHDUError;
    }
    if (self.hdus[@intCast(hdu_number - 1)].content_type != .Image) {
        return FitsError.NotAnImageError;
    }

    var status: c_int = 0;
    _ = cfitsio.fits_movabs_hdu(self.fptr, hdu_number, null, &status);
    try check(status);

    var naxes: [2]c_long = undefined;
    _ = cfitsio.fits_get_img_size(self.fptr, 2, &naxes, &status);
    try check(status);

    if (naxes[0] <= 0 or naxes[1] <= 0) {
        std.log.debug("Image is empty or has invalid dimensions: {d}x{d}", .{ naxes[0], naxes[1] });
        return null;
    }

    const width: usize = @intCast(naxes[0]);
    const height: usize = @intCast(naxes[1]);
    const data = try self.allocator.alloc(f32, width * height);
    errdefer self.allocator.free(data);

    _ = cfitsio.fits_read_img(self.fptr, cfitsio.TFLOAT, 1, @intCast(data.len), null, data.ptr, null, &status);
    try check(status);
    return .{ .data = data, .width = width, .height = height };
}

pub fn readImage(self: *Fits, hdu_number: c_int, output_path: ?[]const u8, options: StretchOptions) !void {
    const pixels = try self.readPixels(hdu_number) orelse return;
    defer self.allocator.free(pixels.data);
    try self.applyStretch(pixels.data, @intCast(pixels.width), @intCast(pixels.height), output_path, options);
}

pub fn readImageAsTable(self: *Fits, hdu_number: c_int, output_path: ?[]const u8) !void {
    const pixels = try self.readPixels(hdu_number) orelse return;
    defer self.allocator.free(pixels.data);

    const dir = self.outputDir();
    var output_path_buffer: [std.fs.max_path_bytes]u8 = undefined;
    const output = if (output_path) |op|
        try std.fmt.bufPrint(&output_path_buffer, "{s}/{s}.csv", .{ dir, op })
    else
        try std.fmt.bufPrint(&output_path_buffer, "{s}/{s}_image_data.csv", .{ dir, dir });

    const csv_file = try std.Io.Dir.cwd().createFile(self.io, output, .{});
    defer csv_file.close(self.io);
    var out_buf: [4096]u8 = undefined;
    var csv_writer = csv_file.writer(self.io, &out_buf);
    const out = &csv_writer.interface;

    try out.writeAll("x,y,value\n");
    for (0..pixels.height) |y| {
        const row = pixels.data[y * pixels.width ..][0..pixels.width];
        for (row, 0..) |value, x| {
            try out.print("{d},{d},{d}\n", .{ x, y, value });
        }
    }
    try out.flush();
    std.log.debug("Image data saved as CSV: {s}", .{output});
}

fn applyStretch(self: *Fits, pixels: []f32, width: u32, height: u32, output_path: ?[]const u8, options: StretchOptions) !void {
    // Clip range: drop the darkest 0.05% and brightest 0.25% of pixels.
    const ordered = try self.allocator.dupe(f32, pixels);
    defer self.allocator.free(ordered);
    std.mem.sort(f32, ordered, {}, std.sort.asc(f32));
    const lo = ordered[ordered.len / 2000];
    const hi = ordered[ordered.len * 1995 / 2000];
    const span = hi - lo;

    var img: zignal.Image(Rgba) = try .init(self.allocator, height, width);
    defer img.deinit(self.allocator);

    for (img.data, pixels) |*dst, raw| {
        const t = std.math.clamp((raw - lo) / span, 0.0, 1.0);
        const gray: u8 = @intFromFloat(sineStretch(t, options) * 255);
        dst.* = .{ .r = gray, .g = gray, .b = gray };
    }

    const dir = self.outputDir();
    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const png_path = try std.fmt.bufPrint(&path_buf, "{s}/{s}.png", .{ dir, output_path orelse dir });
    try img.save(self.io, self.allocator, png_path);
}

inline fn sineStretch(x: f32, options: StretchOptions) f32 {
    const stretch: f32 = options.stretch;
    const bend: f32 = options.bend;
    return std.math.asinh((x - bend) / stretch) / std.math.asinh((1 - bend) / stretch) * 0.5 + 0.5;
}

/// Outputs go in a directory named after the input file's stem, created on demand.
fn outputDir(self: *Fits) []const u8 {
    const stem = std.fs.path.stem(std.fs.path.basename(self.file_path));
    std.Io.Dir.cwd().createDir(self.io, stem, .default_dir) catch |err| switch (err) {
        error.PathAlreadyExists => {},
        else => std.log.warn("Failed to create directory: {s}. Error: {}", .{ stem, err }),
    };
    return stem;
}

pub const StretchOptions = struct {
    stretch: f32 = 0.15,
    bend: f32 = 0.5,
};

pub const ParseOptions = struct {
    createImages: bool = false,
    stretchOptions: StretchOptions = .{},
};

pub const FitsError = error{
    OpenError,
    HeaderReadError,
    ReadError,
    NotATableError,
    NotAnImageError,
    InvalidHDUError,
};

pub const HduInfo = struct {
    number: c_int,
    content_type: FitsContentType,
};

pub const FitsContentType = enum {
    Image,
    AsciiTable,
    BinaryTable,
    Unknown,
};

test Fits {
    const test_file_png = "sample_fits/image_1.png";
    const test_file_table = "sample_fits/ascii_table_1.csv";
    const io = std.testing.io;
    var fits_png = try Fits.open_and_parse("test/sample_fits.fits", io, std.testing.allocator, .{ .createImages = true });
    defer fits_png.close();

    var file = try std.Io.Dir.cwd().openFile(io, test_file_png, .{});

    var stat = try file.stat(io);
    try std.testing.expect(stat.size > 0);

    file = try std.Io.Dir.cwd().openFile(io, test_file_table, .{});
    defer file.close(io);

    stat = try file.stat(io);
    try std.testing.expect(stat.size > 0);
}

test "Read FITS table" {
    const io = std.testing.io;
    const test_file_table = "small/binary_table_1.csv";
    var fits_file = try Fits.open_and_parse("test/small.fits", io, std.testing.allocator, .{});
    defer fits_file.close();

    var file = try std.Io.Dir.cwd().openFile(io, test_file_table, .{});
    defer file.close(io);

    const stat = try file.stat(io);
    try std.testing.expect(stat.size > 0);
}
