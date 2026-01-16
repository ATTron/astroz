"""astroz - High-performance astrodynamics library."""

import atexit
from .__version__ import __version__
from ._lib import _lib
from .tle import Tle
from .sgp4 import Sgp4
from . import orbital
from . import coords

# initialize on import
_lib.astroz_init()
atexit.register(_lib.astroz_deinit)


def version() -> str:
    """Return astroz version string."""
    v = _lib.astroz_version()
    return f"{(v >> 16) & 0xFF}.{(v >> 8) & 0xFF}.{v & 0xFF}"


__all__ = ["__version__", "version", "Tle", "Sgp4", "orbital", "coords"]
