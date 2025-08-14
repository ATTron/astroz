const std = @import("std");
const zignal = @import("zignal");
const cfitsio = @import("cfitsio").c;

const Fits = @This();

allocator: std.mem.Allocator,
fptr: ?*cfitsio.fitsfile,
file_path: []const u8,
hdus: []HduInfo,

pub fn open(file_path: []const u8, allocator: std.mem.Allocator) !Fits {
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
        .file_path = file_path,
        .hdus = hdus,
    };
}

pub fn open_and_parse(file_path: []const u8, allocator: std.mem.Allocator, options: ParseOptions) !Fits {
    var fits_file = try Fits.open(file_path, allocator);
    var image_count: usize = 0;
    var ascii_table_count: usize = 0;
    var binary_table_count: usize = 0;

    var output_path_buffer: [std.fs.max_path_bytes]u8 = undefined;
    std.log.debug("Processing HDUs started. . .", .{});
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
                std.log.debug("Starting to read ASCII table...\n", .{});
                try fits_file.readTable(hdu.number, filename);
            },
            .BinaryTable => {
                binary_table_count += 1;
                const filename = try std.fmt.bufPrint(&output_path_buffer, "binary_table_{d}", .{binary_table_count});
                std.log.debug("Starting to read binary table...\n", .{});
                try fits_file.readTable(hdu.number, filename);
            },
            .Unknown => {
                std.log.warn("Unknown Fits HDU found . . . skipping", .{});
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
    const hdu_info = self.hdus[@intCast(hdu_number - 1)];
    switch (hdu_info.content_type) {
        .AsciiTable, .BinaryTable => {},
        else => return FitsError.NotATableError,
    }
    var status: c_int = 0;
    _ = cfitsio.fits_movabs_hdu(self.fptr, hdu_number, null, &status);
    if (status != 0) {
        try self.reportError(status);
        return FitsError.ReadError;
    }
    var rows: c_long = undefined;
    var cols: c_int = undefined;
    _ = cfitsio.fits_get_num_rows(self.fptr, &rows, &status);
    _ = cfitsio.fits_get_num_cols(self.fptr, &cols, &status);
    if (status != 0) {
        try self.reportError(status);
        return FitsError.ReadError;
    }
    std.log.debug("Table has {d} rows and {d} columns\n", .{ rows, cols });
    var file: std.fs.File = undefined;

    const input_path = self.file_path;
    const file_name = std.fs.path.stem(std.fs.path.basename(input_path));

    createDirectoryIfNotExists(file_name) catch |err| {
        std.log.warn("Failed to create directory: {s}. Error: {}", .{ file_name, err });
    };

    var output_path_buffer: [std.fs.max_path_bytes]u8 = undefined;
    if (output_path) |op| {
        const output = try std.fmt.bufPrint(&output_path_buffer, "{s}/{s}.csv", .{ file_name, op });
        file = try std.fs.cwd().createFile(output, .{});
    } else {
        const output = try std.fmt.bufPrint(&output_path_buffer, "{s}/{s}.csv", .{ file_name, file_name });
        file = try std.fs.cwd().createFile(output, .{});
    }

    defer file.close();
    var write_buffer: [4096]u8 = undefined;
    var file_writer = file.writer(&write_buffer);
    var buffered = &file_writer.interface;

    // Write column headers
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
        if (status != 0) {
            try self.reportError(status);
            return FitsError.ReadError;
        }
        const col_name = std.mem.sliceTo(&ttype, 0);
        try buffered.print("{s}", .{col_name});
        if (i < cols) {
            try buffered.writeAll(",");
        }
    }
    try buffered.writeAll("\n");

    var row_buffer = try self.allocator.alloc(u8, @as(usize, @intCast(cols)) * cfitsio.FLEN_VALUE);
    defer self.allocator.free(row_buffer);

    var row_pointers = try self.allocator.alloc([*c]u8, @as(usize, @intCast(cols)));
    defer self.allocator.free(row_pointers);

    for (0..@as(usize, @intCast(cols))) |j| {
        row_pointers[j] = @ptrCast(&row_buffer[j * cfitsio.FLEN_VALUE]);
    }

    const row_ptr: [*c][*c]u8 = @ptrCast(row_pointers.ptr);

    // Write data rows
    var row: c_long = 1;
    while (row <= rows) : (row += 1) {
        i = 1;
        while (i <= cols) : (i += 1) {
            _ = cfitsio.fits_read_col_str(self.fptr, i, row, 1, 1, null, row_ptr, null, &status);
            if (status != 0) {
                try self.reportError(status);
                return FitsError.ReadError;
            }
            const value_str = std.mem.sliceTo(row_ptr[0], 0);
            try buffered.print("{s}", .{value_str});
            if (i < cols) {
                try buffered.writeAll(",");
            }
        }
        try buffered.writeAll("\n");
    }
    try file_writer.interface.flush();
    std.log.debug("Finished processing all columns and rows. CSV file created: {any}\n", .{file});
}

fn reportError(self: *Fits, status: c_int) !void {
    _ = self;
    var error_text: [128]u8 = undefined;
    _ = cfitsio.fits_get_errstatus(status, &error_text);
    std.debug.print("CFITSIO Error: {s}\n", .{std.mem.sliceTo(&error_text, 0)});
}

pub fn readImage(self: *Fits, hdu_number: c_int, output_path: ?[]const u8, options: StretchOptions) !void {
    if (hdu_number < 1 or hdu_number > self.hdus.len) {
        return FitsError.InvalidHDUError;
    }

    const hdu_info = self.hdus[@intCast(hdu_number - 1)];

    if (hdu_info.content_type != .Image) {
        return FitsError.NotAnImageError;
    }

    var status: c_int = 0;

    _ = cfitsio.fits_movabs_hdu(self.fptr, hdu_number, null, &status);
    if (status != 0) {
        try self.reportError(status);
        return FitsError.ReadError;
    }

    var naxis: c_int = 0;
    var naxes: [2]c_long = undefined;
    _ = cfitsio.fits_get_img_dim(self.fptr, &naxis, &status);
    _ = cfitsio.fits_get_img_size(self.fptr, 2, &naxes, &status);
    if (status != 0) {
        try self.reportError(status);
        return FitsError.ReadError;
    }

    // this checks if the image would be empty
    if (naxes[0] > 0 and naxes[1] > 0) {
        const width: usize = @intCast(naxes[0]);
        const height: usize = @intCast(naxes[1]);
        const pixels: []f32 = try self.allocator.alloc(f32, width * height);
        defer self.allocator.free(pixels);

        const fpixel: c_long = 1;
        const nelem: c_long = @intCast(width * height);
        _ = cfitsio.fits_read_img(self.fptr, cfitsio.TFLOAT, fpixel, nelem, null, pixels.ptr, null, &status);
        if (status != 0) {
            try self.reportError(status);
            return FitsError.ReadError;
        }

        try self.applyStretch(pixels, width, height, output_path, options);
    } else {
        std.log.debug("Image is empty or invalid dimensions: {d}x{d}\n", .{ naxes[0], naxes[1] });
    }
}

pub fn readImageAsTable(self: *Fits, hdu_number: c_int, output_path: ?[]const u8) !void {
    if (hdu_number < 1 or hdu_number > self.hdus.len) {
        return FitsError.InvalidHDUError;
    }

    const hdu_info = self.hdus[@intCast(hdu_number - 1)];

    if (hdu_info.content_type != .Image) {
        return FitsError.NotAnImageError;
    }

    var status: c_int = 0;

    _ = cfitsio.fits_movabs_hdu(self.fptr, hdu_number, null, &status);
    if (status != 0) {
        try self.reportError(status);
        return FitsError.ReadError;
    }

    var naxis: c_int = 0;
    var naxes: [2]c_long = undefined;
    _ = cfitsio.fits_get_img_dim(self.fptr, &naxis, &status);
    _ = cfitsio.fits_get_img_size(self.fptr, 2, &naxes, &status);
    if (status != 0) {
        try self.reportError(status);
        return FitsError.ReadError;
    }

    if (naxes[0] > 0 and naxes[1] > 0) {
        const width: usize = @intCast(naxes[0]);
        const height: usize = @intCast(naxes[1]);
        const pixels: []f32 = try self.allocator.alloc(f32, width * height);
        defer self.allocator.free(pixels);

        const fpixel: c_long = 1;
        const nelem: c_long = @intCast(width * height);
        _ = cfitsio.fits_read_img(self.fptr, cfitsio.TFLOAT, fpixel, nelem, null, pixels.ptr, null, &status);
        if (status != 0) {
            try self.reportError(status);
            return FitsError.ReadError;
        }

        const input_path = self.file_path;
        const file_name = std.fs.path.stem(std.fs.path.basename(input_path));

        createDirectoryIfNotExists(file_name) catch |err| {
            std.log.warn("Failed to create directory: {s}. Error: {}", .{ file_name, err });
        };

        var output_path_buffer: [std.fs.max_path_bytes]u8 = undefined;
        const output = if (output_path) |op|
            try std.fmt.bufPrint(&output_path_buffer, "{s}/{s}.csv", .{ file_name, op })
        else
            try std.fmt.bufPrint(&output_path_buffer, "{s}/{s}_image_data.csv", .{ file_name, file_name });

        const file = try std.fs.cwd().createFile(output, .{});
        defer file.close();

        var write_buffer: [4096]u8 = undefined;
        var file_writer = file.writer(&write_buffer);
        var buffered = &file_writer.interface;

        try buffered.writeAll("x,y,value\n");

        for (0..height) |y| {
            for (0..width) |x| {
                const pixel_value = pixels[y * width + x];
                try buffered.print("{d},{d},{d}\n", .{ x, y, pixel_value });
            }
        }

        try file_writer.interface.flush();
        std.log.debug("Image data saved as CSV: {s}\n", .{output});
    } else {
        std.log.debug("Image is empty or invalid dimensions\n", .{});
    }
}

fn applyStretch(self: *Fits, pixels: []f32, width: usize, height: usize, output_path: ?[]const u8, options: StretchOptions) !void {
    const sorted_pixels = try self.allocator.dupe(f32, pixels);
    defer self.allocator.free(sorted_pixels);
    std.sort.heap(f32, sorted_pixels, {}, std.sort.asc(f32));

    const vmin_idx = sorted_pixels.len / 2000;
    const vmax_idx = sorted_pixels.len * 1995 / 2000;
    const vmin = sorted_pixels[vmin_idx];
    const vmax = sorted_pixels[vmax_idx];

    var image: zignal.Image(zignal.Rgb) = try .initAlloc(self.allocator, height, width);
    defer image.deinit(self.allocator);
    std.debug.assert(image.data.len == pixels.len);

    for (pixels, image.data) |pixel, *i| {
        const normalized = std.math.clamp((pixel - vmin) / (vmax - vmin), 0, 1);
        const stretched = sineStretch(normalized, options);
        const color = applyColorMap(stretched);

        i.* = .{
            .r = @intFromFloat(color[0] * 255),
            .g = @intFromFloat(color[1] * 255),
            .b = @intFromFloat(color[2] * 255),
        };
    }

    const input_path = self.file_path;
    const file_name = std.fs.path.stem(std.fs.path.basename(input_path));

    createDirectoryIfNotExists(file_name) catch |err| {
        std.log.warn("Failed to create directory: {s}. Error: {}", .{ file_name, err });
    };

    var output_path_buffer: [std.fs.max_path_bytes]u8 = undefined;
    if (output_path) |op| {
        const output = try std.fmt.bufPrint(&output_path_buffer, "{s}/{s}.png", .{ file_name, op });
        try image.save(self.allocator, output);
    } else {
        const output = try std.fmt.bufPrint(&output_path_buffer, "{s}/{s}.png", .{ file_name, file_name });
        try image.save(self.allocator, output);
    }
}

inline fn sineStretch(x: f32, options: StretchOptions) f32 {
    const stretch: f32 = options.stretch;
    const bend: f32 = options.bend;
    return std.math.asinh((x - bend) / stretch) / std.math.asinh((1 - bend) / stretch) * 0.5 + 0.5;
}

inline fn applyColorMap(value: f32) [3]f32 {
    return .{ value, value, value };
}

fn createDirectoryIfNotExists(dir_name: []const u8) !void {
    std.fs.cwd().makeDir(dir_name) catch |err| {
        switch (err) {
            error.PathAlreadyExists => {
                // Directory already exists, which is fine
                return;
            },
            else => {
                // For any other error, we return it
                return err;
            },
        }
    };
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
    MemoryError,
    UnknownError,
    NotATableError,
    NotAnImageError,
    InvalidHDUError,
    InvalidImageSizeError,
    UnsupportedColumnType,
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

pub const FitsTableInfo = struct {
    rows: f32,
    columns: f32,
};

test Fits {
    const test_file_png = "sample_fits/image_1.png";
    const test_file_table = "sample_fits/ascii_table_1.csv";
    var fits_png = try Fits.open_and_parse("test/sample_fits.fits", std.testing.allocator, .{ .createImages = true });
    defer fits_png.close();

    var file = try std.fs.cwd().openFile(test_file_png, .{});

    var stat = try file.stat();
    try std.testing.expect(stat.size > 0);

    file = try std.fs.cwd().openFile(test_file_table, .{});
    defer file.close();

    stat = try file.stat();
    try std.testing.expect(stat.size > 0);
}

test "Read FITS table" {
    const test_file_table = "small/binary_table_1.csv";
    var fits_file = try Fits.open_and_parse("test/small.fits", std.testing.allocator, .{});
    defer fits_file.close();

    var file = try std.fs.cwd().openFile(test_file_table, .{});
    defer file.close();

    const stat = try file.stat();
    try std.testing.expect(stat.size > 0);
}
