"""Type stubs for astroz native bindings."""

from typing import Tuple

import numpy as np
import numpy.typing as npt

# Gravity model constants
WGS84: int
WGS72: int

def version() -> str:
    """Return astroz version string."""

class Tle:
    """Two-Line Element set for satellite orbit data.
    Args:
        tle_string: TLE string (both lines separated by newline)
    """

    def __init__(self, tle_string: str) -> None: ...
    @property
    def satellite_number(self) -> int:
        """NORAD catalog number."""

    @property
    def epoch(self) -> float:
        """Epoch in J2000 seconds."""

    @property
    def inclination(self) -> float:
        """Orbital inclination in degrees."""

    @property
    def eccentricity(self) -> float:
        """Orbital eccentricity."""

    @property
    def mean_motion(self) -> float:
        """Mean motion in revolutions per day."""

class Sgp4:
    """SGP4 orbit propagator with SIMD acceleration.
    Args:
        tle: Tle object
        gravity_model: WGS84 (default) or WGS72
    """

    def __init__(self, tle: Tle, gravity_model: int = WGS84) -> None: ...
    def propagate(
        self, tsince: float
    ) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
        """Propagate to time since TLE epoch.
        Args:
            tsince: Minutes since TLE epoch
        Returns:
            ((x, y, z), (vx, vy, vz)) position [km] and velocity [km/s] in TEME frame
        """

    def propagate_into(
        self,
        times: npt.NDArray[np.float64],
        positions: npt.NDArray[np.float64],
        velocities: npt.NDArray[np.float64],
    ) -> None:
        """Zero-copy batch propagation into pre-allocated arrays.
        Uses SIMD (4-wide) and writes directly to output arrays with no
        Python object creation overhead.
        Args:
            times: float64 array (n,) of minutes since TLE epoch
            positions: float64 array (n, 3) to write positions [km]
            velocities: float64 array (n, 3) to write velocities [km/s]
        Example:
            >>> times = np.linspace(0, 1440, 10000)
            >>> positions = np.empty((len(times), 3), dtype=np.float64)
            >>> velocities = np.empty((len(times), 3), dtype=np.float64)
            >>> sgp4.propagate_into(times, positions, velocities)
        """

    def propagate_batch(
        self,
        times: npt.NDArray[np.float64],
    ) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        """Batch propagation returning position and velocity arrays.
        Convenience wrapper around propagate_into() that allocates and
        returns output arrays.
        Args:
            times: float64 array (n,) of minutes since TLE epoch
        Returns:
            (positions, velocities) tuple of float64 arrays (n, 3)
            positions in km, velocities in km/s, TEME frame
        Example:
            >>> times = np.linspace(0, 1440, 10000)
            >>> positions, velocities = sgp4.propagate_batch(times)
        """
