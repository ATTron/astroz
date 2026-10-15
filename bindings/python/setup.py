"""Builds the astroz extension with `zig build python-bindings`."""

import os
import re
import shutil
import subprocess
import sys
import sysconfig
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


PROJECT_ROOT = Path(__file__).parent.parent.parent


def zon_version():
    """The version from build.zig.zon, the one place it is written."""
    text = (PROJECT_ROOT / "build.zig.zon").read_text()
    return re.search(r'^\s*\.version\s*=\s*"([^"]+)"', text, re.M).group(1)


class ZigBuildExt(build_ext):
    def build_extension(self, ext):
        project_root = PROJECT_ROOT
        cmd = [
            "zig",
            "build",
            "python-bindings",
            f"-Dpython-include={sysconfig.get_path('include')}",
            "-Doptimize=ReleaseFast",
        ]
        # CI sets this (e.g. x86_64-linux-gnu.2.17) so wheels don't pick up the
        # build machine's glibc and CPU
        zig_target = os.environ.get("ASTROZ_ZIG_TARGET")
        if zig_target:
            cmd.append(f"-Dtarget={zig_target}")

        print(f"Building with: {' '.join(cmd)}")
        subprocess.check_call(cmd, cwd=project_root)

        lib_name = "lib_astroz.dylib" if sys.platform == "darwin" else "lib_astroz.so"
        built_lib = project_root / "zig-out" / "bindings" / "python" / "astroz" / lib_name
        if not built_lib.exists():
            raise RuntimeError(f"zig build did not produce {built_lib}")

        target = Path(self.get_ext_fullpath(ext.name))
        target.parent.mkdir(parents=True, exist_ok=True)
        print(f"Copying {built_lib} -> {target}")
        shutil.copy2(built_lib, target)


setup(
    version=zon_version(),
    # No sources: ZigBuildExt does the actual build
    ext_modules=[Extension("astroz._astroz", sources=[])],
    cmdclass={"build_ext": ZigBuildExt},
)
