const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});
    const root_source_file = b.path("src/lib.zig");

    // Module
    const astroz_mod = b.addModule("astroz", .{
        .target = target,
        .optimize = optimize,
        .root_source_file = root_source_file,
    });

    // Library
    const lib_step = b.step("lib", "Install library");

    const lib = b.addStaticLibrary(.{
        .name = "astroz",
        .target = target,
        .optimize = optimize,
        .root_source_file = root_source_file,
    });
    const zigimg_dependency = b.dependency("zigimg", .{
        .target = target,
        .optimize = optimize,
    });

    lib.root_module.addImport("zigimg", zigimg_dependency.module("zigimg"));

    const cfitsio_dep = b.dependency("cfitsio", .{
        .target = target,
        .optimize = optimize,
    });

    lib.root_module.addImport("cfitsio", cfitsio_dep.module("cfitsio"));

    const lib_install = b.addInstallArtifact(lib, .{});
    lib_step.dependOn(&lib_install.step);
    b.default_step.dependOn(lib_step);

    // Documentation
    const doc_step = b.step("doc", "Emit documentation");

    const doc_install = b.addInstallDirectory(.{
        .install_dir = .prefix,
        .install_subdir = "doc",
        .source_dir = lib.getEmittedDocs(),
    });
    doc_step.dependOn(&doc_install.step);
    b.default_step.dependOn(doc_step);

    // Example suite
    const examples_step = b.step("example", "Run example suite");

    inline for (EXAMPLE_NAMES) |EXAMPLE_NAME| {
        const example = b.addExecutable(.{
            .name = EXAMPLE_NAME,
            .target = target,
            .optimize = optimize,
            .root_source_file = b.path("examples/" ++ EXAMPLE_NAME ++ ".zig"),
        });
        example.root_module.addImport("astroz", astroz_mod);

        const example_run = b.addRunArtifact(example);
        examples_step.dependOn(&example_run.step);
    }

    b.default_step.dependOn(examples_step);

    // Test suite
    const tests_step = b.step("test", "Run test suite");

    const tests = b.addTest(.{
        .target = target,
        .root_source_file = root_source_file,
    });

    tests.root_module.addImport("zigimg", zigimg_dependency.module("zigimg"));
    tests.root_module.addImport("cfitsio", cfitsio_dep.module("cfitsio"));

    const tests_run = b.addRunArtifact(tests);
    tests_step.dependOn(&tests_run.step);
    b.default_step.dependOn(tests_step);

    // Formatting checks
    const fmt_step = b.step("fmt", "Run formatting checks");

    const fmt = b.addFmt(.{
        .paths = &.{
            "src/",
            "build.zig",
        },
        .check = true,
    });
    fmt_step.dependOn(&fmt.step);
    b.default_step.dependOn(fmt_step);
}

const EXAMPLE_NAMES = &.{
    "create_ccsds_packet_config",
    "create_ccsds_packet",
    "orbit_phase_change",
    "orbit_plane_change",
    "orbit_prop_impulse",
    "orbit_prop",
    "parse_ccsds_file_sync",
    "parse_ccsds",
    "parse_tle",
    "parse_vita49_callback",
    "parse_vita49",
    "precess_star",
    "simple_spacecraft_orientation",
};
