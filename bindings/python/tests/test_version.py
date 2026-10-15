import re
from importlib.metadata import version
from pathlib import Path

import astroz


def test_version_comes_from_build_zig_zon():
    zon = (Path(__file__).parents[3] / "build.zig.zon").read_text()
    expected = re.search(r'^\s*\.version\s*=\s*"([^"]+)"', zon, re.M).group(1)
    assert astroz.__version__ == expected
    assert version("astroz") == expected
