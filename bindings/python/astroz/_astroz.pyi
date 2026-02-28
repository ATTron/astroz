"""Type stubs for astroz native bindings.

This module contains type stubs for the Zig native extension. These stubs
provide type hints and documentation for IDE support and static analysis.

The main classes are:

- ``Tle``: Parse and access Two-Line Element set data
- ``Sgp4``: Single-satellite SGP4 propagator
- ``Sgp4Constellation``: Multi-satellite constellation propagator with SIMD
"""

from typing import Tuple, Optional, List, Literal, Dict

import numpy as np
import numpy.typing as npt

# Gravity model constants
WGS84: int
"""WGS84 gravity model constant (default)."""

WGS72: int
"""WGS72 gravity model constant (legacy)."""

# Physical constants
EARTH_MU: float
"""Earth gravitational parameter (km^3/s^2), WGS84 value 398600.5."""

EARTH_R_EQ: float
"""Earth equatorial radius (km), WGS84 value 6378.137."""

EARTH_J2: float
"""Earth J2 zonal harmonic coefficient, WGS84 value 0.00108262998905."""

SUN_MU: float
"""Sun gravitational parameter (km^3/s^2), 1.32712e11."""

MOON_MU: float
"""Moon gravitational parameter (km^3/s^2), 4902.80."""

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

def jday(
    year: int,
    month: int,
    day: int,
    hour: int,
    minute: int,
    second: float,
) -> Tuple[float, float]:
    """Convert calendar date to Julian date.

    Returns Julian date split into integer (at noon) and fractional parts,
    matching python-sgp4 convention.

    Parameters
    ----------
    year : int
        Year (e.g., 2024).
    month : int
        Month (1-12).
    day : int
        Day of month (1-31).
    hour : int
        Hour (0-23).
    minute : int
        Minute (0-59).
    second : float
        Second (0.0-59.999...).

    Returns
    -------
    jd : float
        Julian date integer part (at noon).
    fr : float
        Fractional part of Julian date.

    Examples
    --------
    >>> jd, fr = jday(2024, 5, 6, 12, 0, 0.0)
    >>> jd + fr  # Full Julian date
    2460437.0
    """

def days2mdhms(
    year: int,
    days: float,
) -> Tuple[int, int, int, int, float]:
    """Convert day of year to month, day, hour, minute, second.

    Parameters
    ----------
    year : int
        Year (used to determine leap year).
    days : float
        Fractional day of year (1.0 = Jan 1 00:00:00).

    Returns
    -------
    month : int
        Month (1-12).
    day : int
        Day of month (1-31).
    hour : int
        Hour (0-23).
    minute : int
        Minute (0-59).
    second : float
        Second (0.0-59.999...).

    Examples
    --------
    >>> days2mdhms(2024, 127.5)  # May 6, noon
    (5, 6, 12, 0, 0.0)
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

def hohmann_transfer(mu: float, r1: float, r2: float) -> Dict[str, float]:
    """Compute Hohmann transfer between two circular orbits.

    Parameters
    ----------
    mu : float
        Gravitational parameter (km^3/s^2).
    r1 : float
        Initial orbit radius (km).
    r2 : float
        Final orbit radius (km).

    Returns
    -------
    dict
        Keys: sma, dv1, dv2, total_dv, transfer_time, transfer_time_days.
    """

def bi_elliptic_transfer(
    mu: float, r1: float, r2: float, r_intermediate: float
) -> Dict[str, float]:
    """Compute bi-elliptic transfer between two circular orbits.

    Parameters
    ----------
    mu : float
        Gravitational parameter (km^3/s^2).
    r1 : float
        Initial orbit radius (km).
    r2 : float
        Final orbit radius (km).
    r_intermediate : float
        Apoapsis radius of the intermediate transfer orbit (km).

    Returns
    -------
    dict
        Keys: sma, dv1, dv2, dv3, total_dv, total_time, total_time_days.
    """

def lambert(
    mu: float,
    r1: Tuple[float, float, float],
    r2: Tuple[float, float, float],
    tof: float,
) -> Dict[str, object]:
    """Solve Lambert's problem for two position vectors and time of flight.

    Parameters
    ----------
    mu : float
        Gravitational parameter (km^3/s^2).
    r1 : tuple of float
        Initial position vector (x, y, z) in km.
    r2 : tuple of float
        Final position vector (x, y, z) in km.
    tof : float
        Time of flight in seconds.

    Returns
    -------
    dict
        Keys: departure_velocity (3-tuple), arrival_velocity (3-tuple),
        transfer_angle, sma, tof.
    """

def orbital_velocity(mu: float, radius: float, sma: Optional[float] = None) -> float:
    """Compute orbital velocity using the vis-viva equation.

    Parameters
    ----------
    mu : float
        Gravitational parameter (km^3/s^2).
    radius : float
        Current orbital radius (km).
    sma : float, optional
        Semi-major axis (km). If None, assumes circular orbit.

    Returns
    -------
    float
        Orbital velocity in km/s.
    """

def orbital_period(mu: float, sma: float) -> float:
    """Compute orbital period from Kepler's third law.

    Parameters
    ----------
    mu : float
        Gravitational parameter (km^3/s^2).
    sma : float
        Semi-major axis (km).

    Returns
    -------
    float
        Orbital period in seconds.
    """

def escape_velocity(mu: float, radius: float) -> float:
    """Compute escape velocity at a given radius.

    Parameters
    ----------
    mu : float
        Gravitational parameter (km^3/s^2).
    radius : float
        Radial distance (km).

    Returns
    -------
    float
        Escape velocity in km/s.
    """

def propagate_numerical(
    state: Tuple[float, float, float, float, float, float],
    t0: float,
    duration: float,
    dt: float,
    mu: float,
    *,
    j2: Optional[float] = None,
    r_eq: Optional[float] = None,
    drag_cd: Optional[float] = None,
    drag_area: Optional[float] = None,
    drag_mass: Optional[float] = None,
    integrator: str = "dp87",
    rtol: float = 1e-9,
    atol: float = 1e-12,
) -> Tuple[List[float], List[Tuple[float, float, float, float, float, float]]]:
    """Numerically propagate an orbit with optional perturbations.

    All positions are in km, velocities in km/s, times in seconds.

    Parameters
    ----------
    state : tuple of float
        Initial state vector [x, y, z, vx, vy, vz] in km and km/s.
    t0 : float
        Start time (s).
    duration : float
        Propagation duration (s).
    dt : float
        Time step (s).
    mu : float
        Central body gravitational parameter (km^3/s^2).
    j2 : float, optional
        J2 zonal harmonic coefficient. Enables J2 perturbation.
        Requires ``r_eq``.
    r_eq : float, optional
        Body equatorial radius (km). Required when ``j2`` or
        ``drag_cd`` is specified.
    drag_cd : float, optional
        Drag coefficient (dimensionless). Enables atmospheric drag
        using an exponential Earth atmosphere model (sea-level density
        1.225 kg/m^3, scale height 7.249 km, cutoff 1500 km).
        Requires ``drag_area``, ``drag_mass``, and ``r_eq``.
    drag_area : float, optional
        Spacecraft cross-sectional area (m^2). Required with ``drag_cd``.
    drag_mass : float, optional
        Spacecraft mass (kg). Required with ``drag_cd``.
    integrator : str, optional
        Integrator: ``"dp87"`` (default, adaptive) or ``"rk4"`` (fixed-step).
    rtol : float, optional
        Relative tolerance for DP87 (default 1e-9).
    atol : float, optional
        Absolute tolerance for DP87 (default 1e-12).

    Returns
    -------
    times : list of float
        Time values (s).
    states : list of tuple
        State vectors (x, y, z, vx, vy, vz) in km and km/s at each time.
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

class Satrec:
    """SGP4 satellite record (python-sgp4 compatible).

    This class provides a drop-in replacement for python-sgp4's Satrec class,
    using the high-performance astroz Zig backend with SIMD acceleration.

    Use ``Satrec.twoline2rv()`` class method to create from TLE lines.

    Attributes
    ----------
    satnum : int
        NORAD catalog number.
    epochyr : int
        Epoch year (2-digit).
    epochdays : float
        Epoch day of year (fractional).
    jdsatepoch : float
        Epoch Julian date (integer part).
    jdsatepochF : float
        Epoch Julian date (fractional part).
    ecco : float
        Eccentricity.
    inclo : float
        Inclination (radians).
    nodeo : float
        Right ascension of ascending node (radians).
    argpo : float
        Argument of perigee (radians).
    mo : float
        Mean anomaly (radians).
    no_kozai : float
        Mean motion (radians/minute, Kozai).
    bstar : float
        B* drag term.
    ndot : float
        First derivative of mean motion (rad/min^2).
    a : float
        Semi-major axis (Earth radii).
    alta : float
        Apoapsis altitude (Earth radii).
    altp : float
        Periapsis altitude (Earth radii).
    error : int
        Last error code (0 = success).
    t : float
        Last propagation time (minutes from epoch).

    Examples
    --------
    >>> from astroz.api import Satrec, jday, WGS72
    >>> line1 = "1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995"
    >>> line2 = "2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123"
    >>> sat = Satrec.twoline2rv(line1, line2, WGS72)
    >>> jd, fr = jday(2024, 5, 6, 19, 52, 0.0)
    >>> error, position, velocity = sat.sgp4(jd, fr)
    >>> print(f"Position: {position}")  # (x, y, z) in km

    See Also
    --------
    jday : Convert calendar date to Julian date.
    days2mdhms : Convert day of year to calendar components.
    """

    @classmethod
    def twoline2rv(
        cls,
        line1: str,
        line2: str,
        whichconst: int = WGS72,
    ) -> "Satrec":
        """Create Satrec from TLE lines.

        Parameters
        ----------
        line1 : str
            First line of TLE (69 characters).
        line2 : str
            Second line of TLE (69 characters).
        whichconst : int, default WGS72
            Gravity model. Either ``WGS72`` or ``WGS84``.

        Returns
        -------
        Satrec
            Initialized satellite record ready for propagation.

        Examples
        --------
        >>> sat = Satrec.twoline2rv(line1, line2, WGS72)
        """

    def sgp4(
        self,
        jd: float,
        fr: float,
    ) -> Tuple[int, Tuple[float, float, float], Tuple[float, float, float]]:
        """Propagate to given Julian date.

        Parameters
        ----------
        jd : float
            Julian date (integer part, from jday()).
        fr : float
            Julian date (fractional part, from jday()).

        Returns
        -------
        error : int
            Error code (0 = success).
        position : tuple of float
            (x, y, z) position in km, TEME frame.
        velocity : tuple of float
            (vx, vy, vz) velocity in km/s, TEME frame.

        Examples
        --------
        >>> jd, fr = jday(2024, 5, 6, 12, 0, 0.0)
        >>> error, position, velocity = sat.sgp4(jd, fr)
        >>> if error == 0:
        ...     print(f"Position: {position}")
        """

    def sgp4_array(
        self,
        jd: npt.NDArray[np.float64],
        fr: npt.NDArray[np.float64],
    ) -> Tuple[
        npt.NDArray[np.uint8],
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
    ]:
        """Propagate to multiple times using SIMD acceleration.

        Parameters
        ----------
        jd : ndarray
            Julian date integer parts, shape (n_times,).
        fr : ndarray
            Julian date fractional parts, shape (n_times,).

        Returns
        -------
        e : ndarray
            Error codes, shape (n_times,). 0 = success.
        r : ndarray
            Positions in km (TEME), shape (n_times, 3).
        v : ndarray
            Velocities in km/s (TEME), shape (n_times, 3).

        Examples
        --------
        >>> jd = np.array([2458826, 2458826, 2458826])
        >>> fr = np.array([0.0001, 0.0002, 0.0003])
        >>> e, r, v = sat.sgp4_array(jd, fr)
        """

    def sgp4_array_into(
        self,
        jd: npt.NDArray[np.float64],
        fr: npt.NDArray[np.float64],
        positions: npt.NDArray[np.float64],
        velocities: npt.NDArray[np.float64],
    ) -> None:
        """SIMD batch propagation writing directly into pre-allocated arrays.

        Uses SIMD batch propagation internally for both SGP4 and SDP4 satellites.

        Parameters
        ----------
        jd : ndarray
            Julian date integer parts, shape (n_times,).
        fr : ndarray
            Julian date fractional parts, shape (n_times,).
        positions : ndarray
            Output array for positions in km (TEME), shape (n_times, 3).
        velocities : ndarray
            Output array for velocities in km/s (TEME), shape (n_times, 3).

        Examples
        --------
        >>> r = np.empty((100, 3), dtype=np.float64)
        >>> v = np.empty((100, 3), dtype=np.float64)
        >>> sat.sgp4_array_into(jd, fr, r, v)
        """

    # Properties
    @property
    def satnum(self) -> int:
        """NORAD catalog number."""

    @property
    def epochyr(self) -> int:
        """Epoch year (2-digit)."""

    @property
    def epochdays(self) -> float:
        """Epoch day of year (fractional)."""

    @property
    def jdsatepoch(self) -> float:
        """Epoch Julian date (integer part)."""

    @property
    def jdsatepochF(self) -> float:
        """Epoch Julian date (fractional part)."""

    @property
    def ecco(self) -> float:
        """Eccentricity."""

    @property
    def inclo(self) -> float:
        """Inclination (radians)."""

    @property
    def nodeo(self) -> float:
        """Right ascension of ascending node (radians)."""

    @property
    def argpo(self) -> float:
        """Argument of perigee (radians)."""

    @property
    def mo(self) -> float:
        """Mean anomaly (radians)."""

    @property
    def no_kozai(self) -> float:
        """Mean motion (radians/minute, Kozai)."""

    @property
    def bstar(self) -> float:
        """B* drag term."""

    @property
    def ndot(self) -> float:
        """First derivative of mean motion (rad/min^2)."""

    @property
    def a(self) -> float:
        """Semi-major axis (Earth radii)."""

    @property
    def alta(self) -> float:
        """Apoapsis altitude (Earth radii)."""

    @property
    def altp(self) -> float:
        """Periapsis altitude (Earth radii)."""

    @property
    def error(self) -> int:
        """Last error code (0 = success)."""

    @property
    def t(self) -> float:
        """Last propagation time (minutes from epoch)."""

    @property
    def is_deep_space(self) -> bool:
        """True if using SDP4 deep-space propagator (period > 225 min)."""

class SatrecArray:
    """Batch SGP4 propagator for multiple satellites (python-sgp4 compatible).

    This class provides SIMD-accelerated batch propagation for multiple satellites,
    achieving 270-330M propagations/second (positions only) or ~200M props/sec
    (with velocities).

    Parameters
    ----------
    satrecs : list of Satrec
        List of Satrec objects to propagate as a batch.

    Attributes
    ----------
    num_satellites : int
        Number of satellites in the array.
    epochs : list of float
        Epoch Julian dates for each satellite.

    Examples
    --------
    >>> from astroz.api import Satrec, SatrecArray, jday, WGS72
    >>> import numpy as np
    >>>
    >>> # Create satellites
    >>> sats = [Satrec.twoline2rv(l1, l2, WGS72) for l1, l2 in tle_pairs]
    >>> sat_array = SatrecArray(sats)
    >>>
    >>> # Propagate to multiple times
    >>> jd = np.full(1440, 2460000.5)
    >>> fr = np.linspace(0, 1, 1440)
    >>> e, r, v = sat_array.sgp4(jd, fr)
    >>>
    >>> # Faster: skip velocity computation
    >>> e, r, _ = sat_array.sgp4(jd, fr, velocities=False)

    See Also
    --------
    Satrec : Single satellite record.
    """

    def __init__(self, satrecs: List[Satrec]) -> None: ...
    @property
    def num_satellites(self) -> int:
        """Number of satellites in the array."""

    @property
    def epochs(self) -> List[float]:
        """Epoch Julian dates for each satellite."""

    def sgp4(
        self,
        jd: float | npt.NDArray[np.float64],
        fr: float | npt.NDArray[np.float64],
        *,
        velocities: bool = True,
    ) -> Tuple[
        npt.NDArray[np.uint8],
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
    ]:
        """Propagate all satellites to the given Julian dates.

        Uses SIMD internally for high throughput.

        Parameters
        ----------
        jd : float or ndarray
            Julian date integer parts. Scalar or array of shape (n_times,).
        fr : float or ndarray
            Julian date fractional parts. Scalar or array of shape (n_times,).
        velocities : bool, default True
            If True, compute velocities. Set False for ~30% faster propagation.

        Returns
        -------
        e : ndarray
            Error codes, shape (n_sats, n_times). 0 = success.
        r : ndarray
            Positions in km (TEME), shape (n_sats, n_times, 3).
        v : ndarray
            Velocities in km/s (TEME), shape (n_sats, n_times, 3).
            All zeros if velocities=False.

        Examples
        --------
        >>> e, r, v = sat_array.sgp4(2458826.5, 0.8625)  # Single time (scalars)
        >>> e, r, v = sat_array.sgp4(jd_arr, fr_arr)  # Multiple times (arrays)
        >>> e, r, _ = sat_array.sgp4(jd, fr, velocities=False)  # Faster
        """
