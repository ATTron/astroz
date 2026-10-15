"""Type stubs for the astroz native extension (``astroz._astroz``).

These describe the native types only. The ``astroz.api.Satrec`` and
``astroz.api.SatrecArray`` wrappers add ``sgp4_array`` and
``SatrecArray.sgp4`` on top of them.
"""

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import numpy.typing as npt

WGS72OLD: int  # python-sgp4 numbering: WGS72OLD 0, WGS72 1, WGS84 2
WGS72: int
WGS84: int

EARTH_MU: float
"""Earth gravitational parameter (km^3/s^2), WGS84."""
EARTH_R_EQ: float
"""Earth equatorial radius (km), WGS84."""
EARTH_J2: float
"""Earth J2 zonal harmonic, WGS84."""
SUN_MU: float
"""Sun gravitational parameter (km^3/s^2)."""
MOON_MU: float
"""Moon gravitational parameter (km^3/s^2)."""

__version__: str

def jday(
    year: int, month: int, day: int, hour: int, minute: int, second: float
) -> Tuple[float, float]:
    """Calendar date to Julian date, split as python-sgp4 does: ``(jd, fr)``."""

def days2mdhms(year: int, days: float) -> Tuple[int, int, int, int, float]:
    """Fractional day of year (1.0 = Jan 1 00:00) to (month, day, hour, minute, second)."""

def screen_candidates(
    positions: npt.NDArray[np.float64],
    velocities: npt.NDArray[np.float64],
    num_sats: int,
    dt: float,
    threshold: float,
    primary: Optional[npt.NDArray[np.bool_]] = None,
) -> bytes:
    """Candidate intervals that may hold a closest approach within ``threshold`` km.

    ``positions`` and ``velocities`` are time-major ``(n_times, num_sats, 3)``
    TEME states on a grid ``dt`` seconds apart. Returns ``(i, j, step)`` triples
    packed as uint32; read them with ``np.frombuffer(..., np.uint32).reshape(-1, 3)``.
    """

def refine_candidates(
    satrecs: Sequence["Satrec"],
    candidates: npt.NDArray[np.uint32],
    jd0: float,
    dt: float,
) -> bytes:
    """Time of closest approach in each ``(i, j, step)`` candidate from ``screen_candidates``.

    ``satrecs[s]`` propagates satellite ``s`` and the grid started at ``jd0``.
    Returns float64 quadruples (time in seconds after ``jd0``, distance km,
    speed km/s, converged 1.0/0.0); read them with
    ``np.frombuffer(...).reshape(-1, 4)``.
    """

def sdp4_batch_propagate_into(
    satrecs: Sequence["Satrec"],
    jd: npt.NDArray[np.float64],
    fr: npt.NDArray[np.float64],
    positions: npt.NDArray[np.float64],
    velocities: npt.NDArray[np.float64],
    output_stride: int = -1,
    sat_offset: int = 0,
) -> None:
    """Threaded SDP4 for deep-space Satrecs, TEME output.

    Writes ``(n_sats, n_times, 3)`` by default. With ``output_stride > 0`` the
    output is time-major ``(n_times, output_stride, 3)`` and satellite ``s``
    goes in column ``sat_offset + s``.
    """

def hohmann_transfer(mu: float, r1: float, r2: float) -> Dict[str, float]:
    """Hohmann transfer between circular orbits of radius r1 and r2 (km).

    Keys: sma, dv1, dv2, total_dv, transfer_time, transfer_time_days.
    """

def bi_elliptic_transfer(
    mu: float, r1: float, r2: float, r_intermediate: float
) -> Dict[str, float]:
    """Bi-elliptic transfer via an intermediate apoapsis radius (km).

    Keys: sma, dv1, dv2, dv3, total_dv, total_time, total_time_days.
    """

def lambert(
    mu: float,
    r1: Tuple[float, float, float],
    r2: Tuple[float, float, float],
    tof: float,
) -> Dict[str, object]:
    """Lambert's problem for position vectors r1, r2 (km) and time of flight (s).

    Keys: departure_velocity, arrival_velocity (3-tuples), transfer_angle, sma, tof.
    """

def orbital_velocity(mu: float, radius: float, sma: Optional[float] = None) -> float:
    """Vis-viva speed (km/s). Circular orbit if ``sma`` is omitted."""

def orbital_period(mu: float, sma: float) -> float:
    """Orbital period (s)."""

def escape_velocity(mu: float, radius: float) -> float:
    """Escape velocity (km/s)."""

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
    """Numerically propagate a state vector (km, km/s; times in s).

    ``j2`` enables J2 and needs ``r_eq``. ``drag_cd`` enables exponential
    atmosphere drag and needs ``drag_area`` (m^2), ``drag_mass`` (kg) and
    ``r_eq``. ``integrator`` is ``"dp87"`` (adaptive, uses rtol/atol) or
    ``"rk4"``. Returns ``(times, states)``.
    """

class Tle:
    """A parsed two-line element set."""

    def __init__(self, tle_string: str) -> None: ...
    @property
    def satellite_number(self) -> int: ...
    @property
    def epoch(self) -> float:
        """Epoch in seconds since J2000."""
    @property
    def inclination(self) -> float:
        """Degrees."""
    @property
    def eccentricity(self) -> float: ...
    @property
    def mean_motion(self) -> float:
        """Revolutions per day."""

class Sgp4Constellation:
    """SIMD SGP4 batches for many near-earth satellites."""

    def __init__(self, tles: Sequence[Tle | Satrec], gravity_model: int = WGS84) -> None: ...
    @classmethod
    def from_tle_text(cls, text: str, gravity_model: int = WGS84) -> "Sgp4Constellation":
        """Parse every TLE in ``text``, skipping deep-space ones."""
    @property
    def num_satellites(self) -> int: ...
    @property
    def num_batches(self) -> int: ...
    @property
    def batch_size(self) -> int: ...
    @property
    def epochs(self) -> List[float]:
        """Epoch of each satellite as a Julian date."""
    def propagate_into(
        self,
        times: npt.NDArray[np.float64],
        positions: npt.NDArray[np.float64],
        velocities: Optional[npt.NDArray[np.float64]] = None,
        epoch_offsets: Optional[npt.NDArray[np.float64]] = None,
        satellite_mask: Optional[npt.NDArray[np.uint8]] = None,
        output: str = "teme",
        reference_jd: float = 0.0,
        time_major: bool = False,
        output_stride: int = -1,
    ) -> None:
        """Propagate every satellite into preallocated arrays.

        ``times`` are minutes; ``epoch_offsets[i]`` is added for satellite i.
        Output is ``(n_sats, n_times, 3)``, or ``(n_times, n_sats, 3)`` with
        ``time_major``. ``output`` is ``"teme"``, ``"ecef"`` or ``"geodetic"``
        (the last two need ``reference_jd``). ``output_stride`` widens the
        satellite axis so other results can share the array.
        """

class Satrec:
    """python-sgp4 compatible satellite record (SGP4 or SDP4)."""

    @classmethod
    def twoline2rv(cls, line1: str, line2: str, whichconst: int = WGS72) -> "Satrec": ...
    @classmethod
    def from_omm(cls, json: str, whichconst: int = WGS72) -> list["Satrec"]:
        """One Satrec per record in CCSDS OMM JSON (one object or an array)."""
    def sgp4init(
        self, whichconst: int, opsmode: str, satnum: int | str, epoch: float, bstar: float, ndot: float,
        nddot: float, ecco: float, argpo: float, inclo: float, mo: float, no_kozai: float, nodeo: float,
    ) -> None:
        """Initialize from mean elements, as python-sgp4: ``epoch`` in days since 1949-12-31 00:00 UT,
        angles in radians, mean motion in rad/min. ``opsmode`` must be ``'i'``."""
    def sgp4(
        self, jd: float, fr: float
    ) -> Tuple[int, Tuple[float, float, float], Tuple[float, float, float]]:
        """``(error, position km, velocity km/s)`` in TEME."""
    def sgp4_array_into(
        self,
        jd: npt.NDArray[np.float64],
        fr: npt.NDArray[np.float64],
        positions: npt.NDArray[np.float64],
        velocities: npt.NDArray[np.float64],
    ) -> None:
        """Fill ``(n, 3)`` arrays; failed times are NaN."""
    @property
    def satnum(self) -> int: ...
    @property
    def satnum_str(self) -> str: ...
    intldesg: str
    classification: str
    elnum: int
    revnum: int
    ephtype: int
    @property
    def epochyr(self) -> int: ...
    @property
    def epochdays(self) -> float: ...
    @property
    def jdsatepoch(self) -> float: ...
    @property
    def jdsatepochF(self) -> float: ...
    @property
    def ecco(self) -> float: ...
    @property
    def inclo(self) -> float:
        """Radians."""
    @property
    def nodeo(self) -> float:
        """Radians."""
    @property
    def argpo(self) -> float:
        """Radians."""
    @property
    def mo(self) -> float:
        """Radians."""
    @property
    def no_kozai(self) -> float:
        """Radians per minute."""
    @property
    def bstar(self) -> float: ...
    @property
    def ndot(self) -> float:
        """Radians per minute^2."""
    @property
    def nddot(self) -> float:
        """Radians per minute^3."""
    @property
    def a(self) -> float:
        """Semi-major axis in Earth radii."""
    @property
    def alta(self) -> float:
        """Apoapsis altitude in Earth radii."""
    @property
    def altp(self) -> float:
        """Periapsis altitude in Earth radii."""
    @property
    def no_unkozai(self) -> float:
        """Brouwer (un-Kozai'd) mean motion, radians per minute."""
    @property
    def method(self) -> str:
        """``'n'`` for SGP4 (near-earth), ``'d'`` for SDP4 (deep space)."""
    @property
    def operationmode(self) -> str: ...
    @property
    def error(self) -> int: ...
    @property
    def t(self) -> float:
        """Minutes since epoch of the last propagation."""
    @property
    def is_deep_space(self) -> bool: ...

class SatrecArray:
    """SIMD SGP4 batches built from near-earth Satrecs."""

    def __init__(self, satrecs: Sequence[Satrec]) -> None: ...
    @property
    def num_satellites(self) -> int: ...
    @property
    def epochs(self) -> List[float]: ...
    def propagate_into(
        self,
        times: npt.NDArray[np.float64],
        positions: npt.NDArray[np.float64],
        velocities: Optional[npt.NDArray[np.float64]] = None,
        epoch_offsets: Optional[npt.NDArray[np.float64]] = None,
    ) -> None:
        """Time-major TEME output, shape ``(n_times, n_sats, 3)``."""
