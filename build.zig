const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});
    const root_source_file = b.path("src/lib.zig");
    const use_llvm = b.option(bool, "use-llvm", "Use Zig's llvm backend");

    // Module
    const astroz_mod = b.addModule("astroz", .{
        .target = target,
        .optimize = optimize,
        .root_source_file = root_source_file,
    });

    // Library
    const lib_step = b.step("lib", "Install library");

    const lib = b.addLibrary(.{
        .name = "astroz",
        .root_module = astroz_mod,
        .use_llvm = use_llvm,
    });

    const zignal_dependency = b.dependency("zignal", .{
        .target = target,
        .optimize = optimize,
    });

    lib.root_module.addImport("zignal", zignal_dependency.module("zignal"));
    astroz_mod.addImport("zignal", zignal_dependency.module("zignal"));

    const cfitsio_dep = b.dependency("cfitsio", .{
        .target = target,
        .optimize = optimize,
    });

    lib.root_module.addImport("cfitsio", cfitsio_dep.module("cfitsio"));
    astroz_mod.addImport("cfitsio", cfitsio_dep.module("cfitsio"));

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
            .root_module = b.createModule(.{
                .target = target,
                .root_source_file = b.path("examples/" ++ EXAMPLE_NAME ++ ".zig"),
                .optimize = optimize,
            }),
        });
        example.root_module.addImport("astroz", astroz_mod);

        const example_run = b.addRunArtifact(example);
        examples_step.dependOn(&example_run.step);
    }

    // Test suite
    const tests_step = b.step("test", "Run test suite");

    const tests = b.addTest(.{
        .root_module = astroz_mod,
    });

    tests.root_module.addImport("zignal", zignal_dependency.module("zignal"));
    tests.root_module.addImport("cfitsio", cfitsio_dep.module("cfitsio"));

    const tests_run = b.addRunArtifact(tests);
    tests_step.dependOn(&tests_run.step);
    // b.default_step.dependOn(tests_step);

    // C API shared library (for Python/FFI bindings)
    const c_api_step = b.step("c-api", "Build C API shared library");

    const c_api_mod = b.createModule(.{
        .target = target,
        .optimize = optimize,
        .root_source_file = b.path("src/c_api/root.zig"),
    });
    c_api_mod.addImport("astroz", astroz_mod);

    const c_api_lib = b.addLibrary(.{
        .linkage = .dynamic,
        .name = "astroz_c",
        .root_module = c_api_mod,
        .use_llvm = use_llvm,
    });

    const c_api_install = b.addInstallArtifact(c_api_lib, .{});
    c_api_step.dependOn(&c_api_install.step);

    // Benchmark
    const bench_step = b.step("bench", "Run SGP4 benchmark");

    const bench = b.addExecutable(.{
        .name = "sgp4_bench",
        .root_module = b.createModule(.{
            .target = target,
            .root_source_file = b.path("benchmarks/zig_sgp4_bench.zig"),
            .optimize = .ReleaseFast,
        }),
    });
    bench.root_module.addImport("astroz", astroz_mod);

    const bench_run = b.addRunArtifact(bench);
    bench_step.dependOn(&bench_run.step);

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
    "orbit_maneuvers",
    "parse_ccsds_file_sync",
    "parse_ccsds",
    "parse_fits_file",
    "parse_vita49_callback",
    "parse_vita49",
    "precess_star",
    "sgp4_propagation",
    "simple_monte_carlo",
    "simple_spacecraft_orientation",
    "transfer_propagation",
    "wcs",
};
