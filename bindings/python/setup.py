"""Setup script for astroz Python bindings using Zig."""

import subprocess
import sys
import sysconfig
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


class ZigBuildExt(build_ext):
    """Custom build_ext that uses Zig to build native Python bindings."""

    def build_extension(self, ext):
        # Get Python configuration
        python_include = sysconfig.get_path("include")
        python_lib_dir = sysconfig.get_config_var("LIBDIR")

        # Get Python library name (e.g., "python3.12")
        py_version = f"{sys.version_info.major}.{sys.version_info.minor}"
        python_lib = f"python{py_version}"

        # Get the output directory
        ext_fullpath = Path(self.get_ext_fullpath(ext.name))
        ext_dir = ext_fullpath.parent
        ext_dir.mkdir(parents=True, exist_ok=True)

        # Project root is 2 levels up from bindings/python/
        project_root = Path(__file__).parent.parent.parent

        # Build with Zig
        cmd = [
            "zig",
            "build",
            "python-bindings",
            f"-Dpython-include={python_include}",
            "-Doptimize=ReleaseFast",
        ]

        # On non-macOS platforms, link against Python library
        # On macOS, Python extensions should NOT link against the Python library -
        # symbols resolve at runtime when loaded by the interpreter
        if sys.platform != "darwin":
            cmd.append(f"-Dpython-lib={python_lib}")
            if python_lib_dir and Path(python_lib_dir).exists():
                cmd.append(f"-Dpython-lib-path={python_lib_dir}")

        print(f"Building with: {' '.join(cmd)}")
        subprocess.check_call(cmd, cwd=project_root)

        # Copy the built library to the expected location
        # Zig outputs to zig-out/bindings/python/astroz/lib_astroz.so
        built_lib = (
            project_root
            / "zig-out"
            / "bindings"
            / "python"
            / "astroz"
            / "lib_astroz.so"
        )

        if not built_lib.exists():
            # Try alternative names
            for name in ["lib_astroz.so", "_astroz.so", "lib_astroz.dylib"]:
                alt = project_root / "zig-out" / "bindings" / "python" / "astroz" / name
                if alt.exists():
                    built_lib = alt
                    break

        if not built_lib.exists():
            # Check in lib directory
            for name in ["lib_astroz.so", "_astroz.so"]:
                alt = project_root / "zig-out" / "lib" / name
                if alt.exists():
                    built_lib = alt
                    break

        if not built_lib.exists():
            raise RuntimeError(
                f"Built library not found. Looked in:\n"
                f"  {project_root / 'zig-out' / 'bindings' / 'python' / 'astroz'}\n"
                f"  {project_root / 'zig-out' / 'lib'}"
            )

        # Copy to target
        import shutil

        target = ext_dir / f"_astroz{sysconfig.get_config_var('EXT_SUFFIX')}"
        print(f"Copying {built_lib} -> {target}")
        shutil.copy2(built_lib, target)


# Dummy extension - actual build is done by ZigBuildExt
ext_modules = [
    Extension(
        "astroz._astroz",
        sources=[],  # No sources - built by Zig
    )
]

setup(
    ext_modules=ext_modules,
    cmdclass={"build_ext": ZigBuildExt},
)
