const std = @import("std");
const oma = @import("oma");

const examples = [_][]const u8{
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
    "propagation",
    "simple_monte_carlo",
    "simple_spacecraft_orientation",
    "spice_propagation",
    "transfer_propagation",
    "wcs",
};

const Kernel = struct { remote: []const u8, local: []const u8 };

const naif_base_url = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/";
const kernel_dir = "data/kernels";
const kernels = [_]Kernel{
    .{ .remote = "lsk/naif0012.tls", .local = "naif0012.tls" },
    .{ .remote = "spk/planets/de440s.bsp", .local = "de440s.bsp" },
    .{ .remote = "pck/pck00011.tpc", .local = "pck00011.tpc" },
    .{ .remote = "pck/gm_de440.tpc", .local = "gm_de440.tpc" },
};

const cspice_include_fallbacks = [_][]const u8{
    "/usr/include",
    "/usr/local/include",
    "/usr/local/include/cspice",
    "/opt/cspice/include",
};

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const use_llvm = b.option(bool, "use-llvm", "Compile with the LLVM backend");
    const enable_cspice = b.option(bool, "enable-cspice", "Link against NAIF CSPICE") orelse false;
    const cspice_include = b.option([]const u8, "cspice-include", "Directory containing the CSPICE headers");
    const cspice_lib = b.option([]const u8, "cspice-lib", "Path to the CSPICE static archive (cspice.a)");
    const python_include = b.option([]const u8, "python-include", "Directory containing Python.h");

    const build_options = b.addOptions();
    build_options.addOption(bool, "enable_cspice", enable_cspice);
    const build_options_mod = build_options.createModule();

    const zignal_mod = b.dependency("zignal", .{ .target = target, .optimize = optimize }).module("zignal");
    const cfitsio_mod = b.dependency("cfitsio", .{ .target = target, .optimize = optimize }).module("cfitsio");
    const oma_dep = b.dependency("oma", .{});
    const oma_mod = oma_dep.module("oma");
    const kernels_src = b.path("src/simdKernels.zig");

    // Core library module
    const astroz = b.addModule("astroz", .{
        .root_source_file = b.path("src/lib.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "zignal", .module = zignal_mod },
            .{ .name = "cfitsio", .module = cfitsio_mod },
            .{ .name = "oma", .module = oma_mod },
            .{ .name = "build_options", .module = build_options_mod },
        },
    });

    if (enable_cspice) {
        if (cspice_include) |dir| {
            astroz.addIncludePath(.{ .cwd_relative = dir });
        } else for (cspice_include_fallbacks) |dir| {
            astroz.addIncludePath(.{ .cwd_relative = dir });
        }
        astroz.addObjectFile(.{ .cwd_relative = cspice_lib orelse "/usr/lib/cspice.a" });
        astroz.link_libc = true;
    }

    // Static library
    const lib = b.addLibrary(.{
        .name = "astroz",
        .root_module = astroz,
        .use_llvm = use_llvm,
    });
    oma.addMultiVersion(oma_dep, lib, .{ .source = kernels_src });
    const lib_step = b.step("lib", "Build the astroz static library");
    lib_step.dependOn(&b.addInstallArtifact(lib, .{}).step);

    // Documentation
    const docs = b.addInstallDirectory(.{
        .source_dir = lib.getEmittedDocs(),
        .install_dir = .prefix,
        .install_subdir = "doc",
    });
    const doc_step = b.step("doc", "Generate API documentation into zig-out/doc");
    doc_step.dependOn(&docs.step);

    // Examples
    const example_step = b.step("example", "Build and run every example");
    const build_examples_step = b.step("build-examples", "Compile every example without running it");
    for (examples) |name| {
        const exe = b.addExecutable(.{
            .name = name,
            .root_module = b.createModule(.{
                .root_source_file = b.path(b.fmt("examples/{s}.zig", .{name})),
                .target = target,
                .optimize = optimize,
                .imports = &.{.{ .name = "astroz", .module = astroz }},
            }),
        });
        build_examples_step.dependOn(&exe.step);
        example_step.dependOn(&b.addRunArtifact(exe).step);
    }

    // Tests
    const tests = b.addTest(.{ .root_module = astroz });
    const test_step = b.step("test", "Run the unit tests");
    test_step.dependOn(&b.addRunArtifact(tests).step);

    // C API
    const c_api = b.addLibrary(.{
        .linkage = .dynamic,
        .name = "astroz_c",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/c_api/root.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "astroz", .module = astroz }},
        }),
        .use_llvm = use_llvm,
    });
    const c_api_step = b.step("c-api", "Build the astroz C shared library");
    c_api_step.dependOn(&b.addInstallArtifact(c_api, .{}).step);

    // Python extension. The module used here leaves out cfitsio.
    //
    // Example with a uv-managed interpreter:
    //   zig build python-bindings -Doptimize=ReleaseFast \
    //     -Dpython-include=$(uv run python -c "import sysconfig; print(sysconfig.get_path('include'))")
    const astroz_python = b.addModule("astroz_python", .{
        .root_source_file = b.path("src/lib.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "zignal", .module = zignal_mod },
            .{ .name = "oma", .module = oma_mod },
            .{ .name = "build_options", .module = build_options_mod },
        },
    });
    const py_mod = b.createModule(.{
        .root_source_file = b.path("bindings/python/src/main.zig"),
        .target = target,
        .optimize = optimize,
        .link_libc = true,
        .imports = &.{.{ .name = "astroz", .module = astroz_python }},
    });
    if (python_include) |dir| {
        py_mod.addIncludePath(.{ .cwd_relative = dir });
    } else {
        py_mod.addIncludePath(.{ .cwd_relative = "/usr/include/python3.12" });
        py_mod.addIncludePath(.{ .cwd_relative = "/usr/include/python3" });
    }
    const py_lib = b.addLibrary(.{
        .linkage = .dynamic,
        .name = "_astroz",
        .root_module = py_mod,
        .use_llvm = use_llvm,
    });
    // Python symbols come from the interpreter at load time; manylinux Pythons
    // have no libpython.so to link against.
    py_lib.linker_allow_shlib_undefined = true;
    oma.addMultiVersion(oma_dep, py_lib, .{ .source = kernels_src, .pic = true });
    const py_step = b.step("python-bindings", "Build the Python extension module");
    py_step.dependOn(&b.addInstallArtifact(py_lib, .{
        .dest_dir = .{ .override = .{ .custom = "bindings/python/astroz" } },
    }).step);

    // WebAssembly module for the JavaScript bindings. Like the Python module it
    // leaves out the C dependencies.
    const wasm_target = b.resolveTargetQuery(.{
        .cpu_arch = .wasm32,
        .os_tag = .freestanding,
        .cpu_features_add = std.Target.wasm.featureSet(&.{ .simd128, .bulk_memory }),
    });
    const astroz_wasm = b.createModule(.{
        .root_source_file = b.path("src/lib.zig"),
        .target = wasm_target,
        .optimize = .ReleaseFast,
        .imports = &.{.{ .name = "build_options", .module = build_options_mod }},
    });
    const wasm = b.addExecutable(.{
        .name = "astroz",
        .root_module = b.createModule(.{
            .root_source_file = b.path("bindings/javascript/zig/root.zig"),
            .target = wasm_target,
            .optimize = .ReleaseFast,
            .strip = true,
            .imports = &.{.{ .name = "astroz", .module = astroz_wasm }},
        }),
    });
    wasm.entry = .disabled;
    wasm.rdynamic = true;
    const wasm_step = b.step("wasm", "Build the WebAssembly module for the JavaScript bindings");
    wasm_step.dependOn(&b.addInstallArtifact(wasm, .{
        .dest_dir = .{ .override = .{ .custom = "bindings/javascript/wasm" } },
    }).step);

    // Benchmark
    const bench = b.addExecutable(.{
        .name = "sgp4_bench",
        .root_module = b.createModule(.{
            .root_source_file = b.path("benchmarks/zig_sgp4_bench.zig"),
            .target = target,
            .optimize = .ReleaseFast,
            .imports = &.{.{ .name = "astroz", .module = astroz }},
        }),
    });
    const bench_step = b.step("bench", "Run the SGP4 benchmark (ReleaseFast)");
    bench_step.dependOn(&b.addRunArtifact(bench).step);

    // Formatting check
    const fmt = b.addFmt(.{ .paths = &.{ "src", "examples", "build.zig" }, .check = true });
    const fmt_step = b.step("fmt", "Check source formatting");
    fmt_step.dependOn(&fmt.step);

    // NAIF kernel download
    const fetch_step = b.step("fetch-kernels", "Download NAIF generic kernels into " ++ kernel_dir);
    const fetch = b.addSystemCommand(&.{ "sh", "-c", fetchKernelsScript() });
    fetch.setCwd(b.path("."));
    fetch_step.dependOn(&fetch.step);

    b.getInstallStep().dependOn(lib_step);
    b.getInstallStep().dependOn(doc_step);
    b.getInstallStep().dependOn(fmt_step);
}

fn fetchKernelsScript() []const u8 {
    comptime var script: []const u8 =
        \\set -e
        \\if ! command -v curl >/dev/null 2>&1; then
        \\  echo "fetch-kernels: 'curl' was not found in PATH." >&2
        \\  echo "Install it first, for example:" >&2
        \\  echo "  Debian/Ubuntu: sudo apt install curl" >&2
        \\  echo "  Fedora:        sudo dnf install curl" >&2
        \\  echo "  Arch:          sudo pacman -S curl" >&2
        \\  echo "  macOS:         brew install curl" >&2
        \\  exit 1
        \\fi
        \\
        \\get() {
        \\  echo "fetching $2"
        \\  if [ -f "$2" ]; then
        \\    curl --fail --location --silent --show-error -z "$2" -o "$2" "$1"
        \\  else
        \\    curl --fail --location --silent --show-error -o "$2" "$1"
        \\  fi
        \\}
        \\
    ++ "mkdir -p " ++ kernel_dir ++ "\n";
    inline for (kernels) |k| {
        script = script ++ "get " ++ naif_base_url ++ k.remote ++ " " ++ kernel_dir ++ "/" ++ k.local ++ "\n";
    }
    return script;
}
