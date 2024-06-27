const std = @import("std");

// Although this function looks imperative, nfbui will be executed by an external
// runner.
pub fn build(b: *std.Build) void {
    // Standard target options allows the person running `zig build` to choose
    // what target to build for. Here we do not override the defaults, which
    // means any target is allowed, and the default is native. Other options
    // for restricting supported target set are available.
    const target = b.standardTargetOptions(.{});

    // Standard optimization options allow the person running `zig build` to select
    // between Debug, ReleaseSafe, ReleaseFast, and ReleaseSmall. Here we do not
    // set a preferred release mode, allowing the user to decide how to optimize.
    const optimize = b.standardOptimizeOption(.{});

    _ = b.addModule("astroz", .{
        .root_source_file = b.path("src/lib.zig"),
        .target = target,
        .optimize = optimize,
    });

    _ = b.addModule("astroz.ccsds", .{
        .root_source_file = b.path("src/ccsds.zig"),
        .target = target,
        .optimize = optimize,
    });

    _ = b.addModule("astroz.constants", .{
        .root_source_file = b.path("src/constants.zig"),
        .target = target,
        .optimize = optimize,
    });

    _ = b.addModule("astroz.coordinates", .{
        .root_source_file = b.path("src/coordinates.zig"),
        .target = target,
        .optimize = optimize,
    });

    const lib_unit_tests = b.addTest(.{
        .root_source_file = b.path("src/ccsds.zig"),
        .target = target,
        .optimize = optimize,
    });

    const run_lib_unit_tests = b.addRunArtifact(lib_unit_tests);

    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_lib_unit_tests.step);
}
