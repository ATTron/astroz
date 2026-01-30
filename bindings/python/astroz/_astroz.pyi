"""Type stubs for astroz native bindings.

This module contains type stubs for the Zig native extension. These stubs
provide type hints and documentation for IDE support and static analysis.

The main classes are:

- ``Tle``: Parse and access Two-Line Element set data
- ``Sgp4``: Single-satellite SGP4 propagator
- ``Sgp4Constellation``: Multi-satellite constellation propagator with SIMD
"""

from typing import Tuple, Optional, List, Literal

import numpy as np
import numpy.typing as npt

# Gravity model constants
WGS84: int
"""WGS84 gravity model constant (default)."""

WGS72: int
"""WGS72 gravity model constant (legacy)."""

def version() -> str:
    """Return astroz version string.

    Returns
    -------
    str
        Version string in format "major.minor.patch".

    Examples
    --------
    >>> import astroz._astroz as _astroz
    >>> _astroz.version()
    '0.4.4'
    """

def coarse_screen(
    positions: npt.NDArray[np.float64],
    n_satellites: int,
    threshold: float,
    satellite_mask: Optional[npt.NDArray[np.bool_]],
) -> Tuple[npt.NDArray[np.int32], npt.NDArray[np.int32]]:
    """Screen position data for conjunction events.

    Performs all-vs-all distance checking on pre-computed position data
    to identify satellite pairs that come within a threshold distance.

    Parameters
    ----------
    positions : ndarray
        Position array of shape ``(n_times, n_satellites, 3)`` in km.
    n_satellites : int
        Number of satellites in the constellation.
    threshold : float
        Distance threshold in kilometers.
    satellite_mask : ndarray or None
        Boolean mask of shape ``(n_satellites,)`` to exclude satellites,
        or None to include all.

    Returns
    -------
    pairs : ndarray
        Array of shape ``(n_conjunctions, 2)`` with satellite index pairs.
    t_indices : ndarray
        Array of shape ``(n_conjunctions,)`` with time indices.

    See Also
    --------
    Sgp4Constellation.screen_conjunction : Fused propagate+screen for single target.
    """

class Tle:
    """Two-Line Element set for satellite orbit data.

    Parses and provides access to TLE orbital elements.

    Parameters
    ----------
    tle_string : str
        TLE string with both lines separated by newline.

    Attributes
    ----------
    satellite_number : int
        NORAD catalog number (unique satellite identifier).
    epoch : float
        TLE epoch in J2000 seconds.
    inclination : float
        Orbital inclination in degrees.
    eccentricity : float
        Orbital eccentricity (dimensionless, 0-1).
    mean_motion : float
        Mean motion in revolutions per day.

    Examples
    --------
    >>> tle_text = '''ISS (ZARYA)
    ... 1 25544U 98067A   24001.50000000  .00000000  00000-0  00000-0 0    09
    ... 2 25544  51.6400 100.0000 0001000   0.0000   0.0000 15.50000000    05'''
    >>> tle = Tle(tle_text)
    >>> tle.satellite_number
    25544
    >>> tle.inclination
    51.64
    """

    def __init__(self, tle_string: str) -> None: ...
    @property
    def satellite_number(self) -> int:
        """NORAD catalog number (unique satellite identifier)."""

    @property
    def epoch(self) -> float:
        """TLE epoch in J2000 seconds."""

    @property
    def inclination(self) -> float:
        """Orbital inclination in degrees."""

    @property
    def eccentricity(self) -> float:
        """Orbital eccentricity (dimensionless, 0-1)."""

    @property
    def mean_motion(self) -> float:
        """Mean motion in revolutions per day."""

class Sgp4:
    """SGP4 orbit propagator for a single satellite.

    SIMD-accelerated SGP4/SDP4 propagator supporting both near-Earth and
    deep-space satellites.

    Parameters
    ----------
    tle : Tle
        Parsed TLE object containing orbital elements.
    gravity_model : int, default WGS84
        Gravity model to use. Either ``WGS84`` or ``WGS72``.

    Examples
    --------
    >>> tle = Tle(tle_string)
    >>> sgp4 = Sgp4(tle)
    >>> pos, vel = sgp4.propagate(60.0)  # 1 hour after epoch

    See Also
    --------
    Sgp4Constellation : Multi-satellite propagator for constellations.
    """

    def __init__(self, tle: Tle, gravity_model: int = WGS84) -> None: ...
    def propagate(
        self, tsince: float
    ) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
        """Propagate to a single time.

        Parameters
        ----------
        tsince : float
            Minutes since TLE epoch.

        Returns
        -------
        position : tuple of float
            (x, y, z) position in km, TEME frame.
        velocity : tuple of float
            (vx, vy, vz) velocity in km/s, TEME frame.

        Examples
        --------
        >>> pos, vel = sgp4.propagate(60.0)
        >>> x, y, z = pos
        """

    def propagate_into(
        self,
        times: npt.NDArray[np.float64],
        positions: npt.NDArray[np.float64],
        velocities: npt.NDArray[np.float64],
    ) -> None:
        """Batch propagation into pre-allocated arrays.

        Zero-copy propagation using SIMD acceleration (8-wide). Writes
        directly to output arrays with no Python object creation overhead.

        Parameters
        ----------
        times : ndarray
            float64 array of shape ``(n,)`` containing minutes since TLE epoch.
        positions : ndarray
            float64 array of shape ``(n, 3)`` to write positions in km.
        velocities : ndarray
            float64 array of shape ``(n, 3)`` to write velocities in km/s.

        Examples
        --------
        >>> times = np.linspace(0, 1440, 10000)
        >>> positions = np.empty((len(times), 3), dtype=np.float64)
        >>> velocities = np.empty((len(times), 3), dtype=np.float64)
        >>> sgp4.propagate_into(times, positions, velocities)
        """

    def propagate_batch(
        self,
        times: npt.NDArray[np.float64],
    ) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        """Batch propagation returning new arrays.

        Convenience wrapper around ``propagate_into()`` that allocates and
        returns output arrays.

        Parameters
        ----------
        times : ndarray
            float64 array of shape ``(n,)`` containing minutes since TLE epoch.

        Returns
        -------
        positions : ndarray
            float64 array of shape ``(n, 3)`` with positions in km, TEME frame.
        velocities : ndarray
            float64 array of shape ``(n, 3)`` with velocities in km/s, TEME frame.

        Examples
        --------
        >>> times = np.linspace(0, 1440, 10000)
        >>> positions, velocities = sgp4.propagate_batch(times)
        """

class Sgp4Constellation:
    """Multi-satellite constellation propagator with SIMD and threading.

    Efficiently propagates multiple satellites simultaneously using SIMD
    vectorization and multi-threading.

    This class is typically used internally by the high-level ``propagate()``
    and ``screen()`` functions, but can be used directly for advanced use cases.

    Attributes
    ----------
    num_satellites : int
        Number of satellites in the constellation.
    epochs : list of float
        TLE epoch for each satellite as Julian Date.

    Examples
    --------
    >>> constellation = Sgp4Constellation.from_tle_text(tle_text)
    >>> print(f"Loaded {constellation.num_satellites} satellites")

    See Also
    --------
    Sgp4 : Single-satellite propagator.
    """

    @staticmethod
    def from_tle_text(tle_text: str) -> "Sgp4Constellation":
        """Create constellation from TLE text.

        Parameters
        ----------
        tle_text : str
            Multi-line string containing TLE data for multiple satellites.

        Returns
        -------
        Sgp4Constellation
            Constellation object ready for propagation.

        Examples
        --------
        >>> with open("starlink.tle") as f:
        ...     tle_text = f.read()
        >>> constellation = Sgp4Constellation.from_tle_text(tle_text)
        """

    @property
    def num_satellites(self) -> int:
        """Number of satellites in the constellation."""

    @property
    def epochs(self) -> List[float]:
        """TLE epoch for each satellite as Julian Date."""

    def propagate_into(
        self,
        times: npt.NDArray[np.float64],
        positions: npt.NDArray[np.float64],
        velocities: Optional[npt.NDArray[np.float64]],
        *,
        epoch_offsets: npt.NDArray[np.float64],
        satellite_mask: Optional[npt.NDArray[np.bool_]],
        output: Literal["ecef", "teme", "geodetic"],
        reference_jd: float,
        time_major: bool,
    ) -> None:
        """Propagate constellation into pre-allocated arrays.

        Parameters
        ----------
        times : ndarray
            float64 array of time offsets in minutes.
        positions : ndarray
            Output array for positions. Shape depends on ``time_major``.
        velocities : ndarray or None
            Output array for velocities, or None to skip velocity computation.
        epoch_offsets : ndarray
            Per-satellite epoch offsets in minutes.
        satellite_mask : ndarray or None
            Boolean mask to select which satellites to propagate.
        output : {"ecef", "teme", "geodetic"}
            Output coordinate frame.
        reference_jd : float
            Reference Julian Date for time calculations.
        time_major : bool
            If True, output shape is ``(n_times, n_sats, 3)``.
            If False, output shape is ``(n_sats, n_times, 3)``.
        """

    def screen_conjunction(
        self,
        times: npt.NDArray[np.float64],
        target: int,
        threshold: float,
        *,
        epoch_offsets: npt.NDArray[np.float64],
        reference_jd: float,
    ) -> Tuple[List[float], List[int]]:
        """Fused propagate and screen for single-target conjunction detection.

        This is the fastest path for finding close approaches to a specific
        satellite. Combines propagation and distance checking in a single pass.

        Parameters
        ----------
        times : ndarray
            float64 array of time offsets in minutes.
        target : int
            Index of the target satellite to screen against.
        threshold : float
            Distance threshold in kilometers.
        epoch_offsets : ndarray
            Per-satellite epoch offsets in minutes.
        reference_jd : float
            Reference Julian Date for time calculations.

        Returns
        -------
        min_distances : list of float
            Minimum distance to target for each satellite.
        min_t_indices : list of int
            Time index at which minimum distance occurred.

        See Also
        --------
        coarse_screen : All-vs-all screening on pre-computed positions.
        """
