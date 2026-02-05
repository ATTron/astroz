const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});
    const root_source_file = b.path("src/lib.zig");
    const use_llvm = b.option(bool, "use-llvm", "Use Zig's llvm backend");

    // CSPICE configuration
    // Override with: -Dcspice-include=/path/to/cspice/include -Dcspice-lib=/path/to/cspice/lib
    const cspice_include = b.option([]const u8, "cspice-include", "CSPICE include path (containing SpiceUsr.h)");
    const cspice_lib = b.option([]const u8, "cspice-lib", "CSPICE library path (containing libcspice.a)");
    const enable_cspice = b.option(bool, "enable-cspice", "Enable CSPICE support (requires CSPICE to be installed)") orelse false;

    // Build options for conditional compilation
    const build_options = b.addOptions();
    build_options.addOption(bool, "enable_cspice", enable_cspice);

    // Module
    const astroz_mod = b.addModule("astroz", .{
        .target = target,
        .optimize = optimize,
        .root_source_file = root_source_file,
    });

    astroz_mod.addOptions("build_options", build_options);

    // CSPICE setup
    if (enable_cspice) {
        if (cspice_include) |inc| {
            astroz_mod.addIncludePath(.{ .cwd_relative = inc });
        } else {
            astroz_mod.addIncludePath(.{ .cwd_relative = "/usr/include" });
            astroz_mod.addIncludePath(.{ .cwd_relative = "/usr/local/include" });
            astroz_mod.addIncludePath(.{ .cwd_relative = "/usr/local/include/cspice" });
            astroz_mod.addIncludePath(.{ .cwd_relative = "/opt/cspice/include" });
        }

        if (cspice_lib) |lib_path| {
            astroz_mod.addObjectFile(.{ .cwd_relative = lib_path });
        } else {
            // Try AUR location first, then fallback to standard locations
            astroz_mod.addObjectFile(.{ .cwd_relative = "/usr/lib/cspice.a" });
        }
        astroz_mod.link_libc = true;
    }

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

    // Python bindings
    const python_step = b.step("python-bindings", "Build Python native bindings");

    // Create a minimal astroz module for Python bindings (no cfitsio dependency)
    const astroz_python_mod = b.addModule("astroz_python", .{
        .target = target,
        .optimize = optimize,
        .root_source_file = b.path("src/lib.zig"),
    });
    astroz_python_mod.addImport("zignal", zignal_dependency.module("zignal"));

    const python_mod = b.createModule(.{
        .target = target,
        .optimize = optimize,
        .root_source_file = b.path("bindings/python/src/main.zig"),
    });
    python_mod.addImport("astroz", astroz_python_mod);

    // Python configuration for bindings
    // Override with: -Dpython-include=/path -Dpython-lib=python3.X -Dpython-lib-path=/path
    // Example for uv-managed Python 3.12:
    //   zig build python-bindings \
    //     -Dpython-include=$(python3.12 -c "import sysconfig; print(sysconfig.get_path('include'))") \
    //     -Dpython-lib=python3.12 \
    //     -Dpython-lib-path=$(python3.12 -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
    const python_include = b.option([]const u8, "python-include", "Python include path (run: python -c \"import sysconfig; print(sysconfig.get_path('include'))\")");
    const python_lib_name = b.option([]const u8, "python-lib", "Python library name (e.g., python3.12)");
    const python_lib_path = b.option([]const u8, "python-lib-path", "Python library search path");

    if (python_include) |inc| {
        python_mod.addIncludePath(.{ .cwd_relative = inc });
    } else {
        // Try common system locations
        python_mod.addIncludePath(.{ .cwd_relative = "/usr/include/python3.12" });
        python_mod.addIncludePath(.{ .cwd_relative = "/usr/include/python3" });
    }

    const is_macos = target.result.os.tag == .macos;

    // On macOS, don't link against Python library - symbols resolve at load time
    // when the extension is loaded by the Python interpreter. This is the standard
    // approach used by all major Python extension build systems.
    if (!is_macos) {
        if (python_lib_path) |path| {
            python_mod.addLibraryPath(.{ .cwd_relative = path });
        }
        python_mod.linkSystemLibrary(python_lib_name orelse "python3.12", .{});
    }

    python_mod.link_libc = true;

    const python_lib = b.addLibrary(.{
        .linkage = .dynamic,
        .name = "_astroz",
        .root_module = python_mod,
        .use_llvm = use_llvm,
    });

    // On macOS, allow undefined symbols - they resolve when loaded by Python
    if (is_macos) {
        python_lib.linker_allow_shlib_undefined = true;
    }

    const python_install = b.addInstallArtifact(python_lib, .{
        .dest_dir = .{ .override = .{ .custom = "bindings/python/astroz" } },
    });
    python_step.dependOn(&python_install.step);

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
    "maneuver_planning",
    "constellation_phasing",
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
