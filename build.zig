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

    _ = b.addModule("astroz.time", .{
        .root_source_file = b.path("src/time.zig"),
        .target = target,
        .optimize = optimize,
    });

    _ = b.addModule("astroz.calculations", .{
        .root_source_file = b.path("src/calculations.zig"),
        .target = target,
        .optimize = optimize,
    });

    _ = b.addModule("astroz.vita49", .{
        .root_source_file = b.path("src/vita49.zig"),
        .target = target,
        .optimize = optimize,
    });

    _ = b.addModule("astroz.parsers", .{
        .root_source_file = b.path("src/parser.zig"),
        .target = target,
        .optimize = optimize,
    });

    _ = b.addModule("astroz.tle", .{
        .root_source_file = b.path("src/tle.zig"),
        .target = target,
        .optimize = optimize,
    });

    _ = b.addModule("astroz.spacecraft", .{
        .root_source_file = b.path("src/spacecraft.zig"),
        .target = target,
        .optimize = optimize,
    });

    const coord_unit_tests = b.addTest(.{
        .root_source_file = b.path("src/coordinates.zig"),
    });

    const run_coord_unit_tests = b.addRunArtifact(coord_unit_tests);

    const ccsds_unit_tests = b.addTest(.{
        .root_source_file = b.path("src/ccsds.zig"),
    });

    const run_ccsds_unit_tests = b.addRunArtifact(ccsds_unit_tests);

    const time_unit_tests = b.addTest(.{
        .root_source_file = b.path("src/time.zig"),
    });

    const run_time_unit_tests = b.addRunArtifact(time_unit_tests);

    const constants_unit_tests = b.addTest(.{
        .root_source_file = b.path("src/constants.zig"),
    });

    const run_constants_unit_tests = b.addRunArtifact(constants_unit_tests);

    const vita49_unit_tests = b.addTest(.{
        .root_source_file = b.path("src/vita49.zig"),
    });

    const run_vita49_unit_tests = b.addRunArtifact(vita49_unit_tests);

    const parsers_unit_tests = b.addTest(.{
        .root_source_file = b.path("src/parsers.zig"),
    });

    const run_parsers_unit_tests = b.addRunArtifact(parsers_unit_tests);

    const tle_unit_tests = b.addTest(.{
        .root_source_file = b.path("src/tle.zig"),
    });

    const run_tle_unit_tests = b.addRunArtifact(tle_unit_tests);

    const spacecraft_unit_tests = b.addTest(.{
        .root_source_file = b.path("src/spacecraft.zig"),
    });

    const run_spacecraft_unit_tests = b.addRunArtifact(spacecraft_unit_tests);

    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_coord_unit_tests.step);
    test_step.dependOn(&run_ccsds_unit_tests.step);
    test_step.dependOn(&run_time_unit_tests.step);
    test_step.dependOn(&run_constants_unit_tests.step);
    test_step.dependOn(&run_vita49_unit_tests.step);
    test_step.dependOn(&run_parsers_unit_tests.step);
    test_step.dependOn(&run_tle_unit_tests.step);
    test_step.dependOn(&run_spacecraft_unit_tests.step);
}
