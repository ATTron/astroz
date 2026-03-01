"""astroz - High-performance astrodynamics library.

A Python library for satellite orbit propagation, conjunction screening,
orbital mechanics, and numerical propagation, powered by a Zig backend
with SIMD acceleration (AVX2/AVX512).

APIs
----
astroz (this module)
    High-level API for quick propagation and screening with automatic
    TLE loading from CelesTrak, files, or URLs.

astroz.api
    python-sgp4 compatible API. Drop-in replacement - just change the import
    for 2-100x speedup. See :mod:`astroz.api` for migration guide.

Quick Start
-----------
**High-level API** (easiest)::

    import astroz
    positions = astroz.propagate("starlink", [0, 60, 120])

**python-sgp4 compatible** (drop-in replacement)::

    from astroz.api import Satrec, jday  # was: from sgp4.api import ...
    sat = Satrec.twoline2rv(line1, line2)
    e, r, v = sat.sgp4(jd, fr)

Main Functions
--------------
propagate
    Propagate satellites to specified times, returning positions (and optionally
    velocities) in ECEF, TEME, or geodetic coordinates.
screen
    Screen a constellation for conjunction events, either all-vs-all or
    single-target mode.

Orbital Mechanics
-----------------
hohmann_transfer
    Compute Hohmann transfer between two circular orbits.
bi_elliptic_transfer
    Compute bi-elliptic transfer between two circular orbits.
lambert
    Solve Lambert's problem for two position vectors and time of flight.
orbital_velocity
    Compute orbital velocity using the vis-viva equation.
orbital_period
    Compute orbital period from Kepler's third law.
escape_velocity
    Compute escape velocity at a given radius.

Numerical Propagation
---------------------
propagate_numerical
    Numerically propagate an orbit with optional perturbations (J2,
    atmospheric drag). Supports RK4 and Dormand-Prince 8(7) integrators.

Constants
---------
EARTH_MU
    Earth gravitational parameter (km^3/s^2), WGS84 value 398600.5.
EARTH_R_EQ
    Earth equatorial radius (km), WGS84 value 6378.137.
EARTH_J2
    Earth J2 zonal harmonic coefficient, WGS84 value 0.00108262998905.
SUN_MU
    Sun gravitational parameter (km^3/s^2), 1.32712e11.
MOON_MU
    Moon gravitational parameter (km^3/s^2), 4902.80.

Classes
-------
Constellation
    Pre-parsed TLE data for efficient repeated propagation/screening.

Data Sources
------------
The library supports multiple TLE data sources:

- **CelesTrak groups**: ``"starlink"``, ``"gps"``, ``"stations"``, etc.
- **NORAD IDs**: ``norad_id=25544`` (ISS) or ``norad_id=[25544, 48274]``
- **URLs**: Any URL returning TLE-format data
- **Local files**: Path to a TLE file
- **Raw TLE text**: Direct TLE string input

Performance
-----------
=================================  ============
Operation                          Throughput
=================================  ============
Multi-satellite batch              290M props/s
Single satellite ``sgp4_array()``  15M props/s
Single ``sgp4()`` call             2.5M props/s
=================================  ============

See Also
--------
- python-sgp4 migration: :mod:`astroz.api`
- Zig documentation: https://attron.github.io/astroz/zig/
- GitHub: https://github.com/ATTron/astroz
"""

import numpy as np
from datetime import datetime, timezone
from pathlib import Path
from ._astroz import (
    Tle,
    Satrec as _Satrec,
    Sgp4Constellation as _Sgp4Constellation,
    coarse_screen as _coarse_screen,
    sdp4_batch_propagate_into as _sdp4_batch_propagate_into,
    hohmann_transfer,
    bi_elliptic_transfer,
    lambert,
    orbital_velocity,
    orbital_period,
    escape_velocity,
    propagate_numerical,
    EARTH_MU,
    EARTH_R_EQ,
    EARTH_J2,
    SUN_MU,
    MOON_MU,
)

_CELESTRAK_ALIASES = {
    "all": "active",
    "iss": "stations",
    "gps": "gps-ops",
    "glonass": "glo-ops",
}


def _fetch_url(url):
    """Fetch URL content."""
    import urllib.request

    req = urllib.request.Request(url, headers={"User-Agent": "astroz"})
    return urllib.request.urlopen(req, timeout=60).read().decode("utf-8")


def _celestrak_url(group=None, norad_id=None):
    """Build CelesTrak URL."""
    if norad_id is not None:
        ids = (
            ",".join(str(i) for i in norad_id)
            if isinstance(norad_id, (list, tuple))
            else str(norad_id)
        )
        return f"https://celestrak.org/NORAD/elements/gp.php?CATNR={ids}&FORMAT=tle"
    group_name = _CELESTRAK_ALIASES.get(group.lower(), group)
    return f"https://celestrak.org/NORAD/elements/gp.php?GROUP={group_name}&FORMAT=tle"


def _load_tle_text(source, norad_id=None):
    """Load TLE text from various sources."""
    if norad_id is not None:
        return _fetch_url(_celestrak_url(norad_id=norad_id))
    if source is None:
        raise ValueError("Must specify 'source' or 'norad_id'")
    if source.startswith(("http://", "https://")):
        return _fetch_url(source)
    if "1 " in source and "2 " in source:
        return source
    if Path(source).exists():
        return Path(source).read_text()
    return _fetch_url(_celestrak_url(group=source))


def _parse_tle_pairs(tle_text):
    """Parse TLE text into list of (line1, line2) pairs."""
    lines = [l.strip() for l in tle_text.strip().splitlines() if l.strip()]
    pairs = []
    i = 0
    while i < len(lines):
        if lines[i].startswith("1 "):
            if i + 1 < len(lines) and lines[i + 1].startswith("2 "):
                pairs.append((lines[i], lines[i + 1]))
                i += 2
                continue
        elif lines[i].startswith("2 "):
            i += 1
            continue
        # Skip name lines
        i += 1
    return pairs


def _start_jd(start_time):
    """Convert start_time to Julian Date, defaulting to current UTC."""
    if start_time is None:
        start_time = datetime.now(timezone.utc)
    return 2440587.5 + (start_time.timestamp() / 86400.0)


def _pad_epochs_for_simd(epochs_arr):
    """Pad epoch array to multiple of 8 for SIMD alignment."""
    n = len(epochs_arr)
    padded = ((n + 7) // 8) * 8
    if padded > n:
        return np.concatenate([epochs_arr, np.full(padded - n, epochs_arr[-1])])
    return epochs_arr


def _compute_epoch_offsets(constellation, start_time):
    """Compute epoch offsets from start_time."""
    start = _start_jd(start_time)
    epochs = np.array(constellation.epochs, dtype=np.float64)
    return (start - epochs) * 1440.0, start


class Constellation:
    """Pre-parsed TLE data for efficient repeated propagation and screening.

    Use this class when you need to propagate or screen the same satellites
    multiple times. It parses TLE data once and caches the result, avoiding
    the overhead of re-parsing on every call.

    Transparently handles mixed SGP4 (near-earth) and SDP4 (deep-space)
    satellites. Near-earth satellites use SIMD batch propagation, deep-space
    satellites use per-satellite SIMD time-major propagation.

    Parameters
    ----------
    source : str, optional
        TLE data source. Can be one of:

        - CelesTrak group name (e.g., ``"starlink"``, ``"gps"``, ``"stations"``)
        - URL to TLE data
        - Local file path
        - Raw TLE text string

    norad_id : int or list of int, optional
        NORAD catalog ID(s) to fetch from CelesTrak. If provided, ``source``
        is ignored.

    Attributes
    ----------
    num_satellites : int
        Number of satellites in the constellation.
    epochs : list of float
        TLE epoch for each satellite as Julian Date.

    Examples
    --------
    Load from CelesTrak group name:

    >>> c = Constellation("starlink")
    >>> c.num_satellites
    6547

    Load specific satellites by NORAD ID:

    >>> iss = Constellation(norad_id=25544)
    >>> iss.num_satellites
    1

    Propagate multiple times efficiently:

    >>> c = Constellation("gps")
    >>> pos1 = propagate(c, [0, 60, 120])
    >>> pos2 = propagate(c, [180, 240, 300])

    See Also
    --------
    propagate : Propagate satellites to specified times.
    screen : Screen for conjunction events.
    """

    def __init__(self, source=None, *, norad_id=None):
        tle_text = _load_tle_text(source, norad_id)
        pairs = _parse_tle_pairs(tle_text)

        # Separate near-earth (SGP4) and deep-space (SDP4) TLEs.
        # Output ordering: SGP4 sats first [0, n_sgp4), SDP4 after [n_sgp4, n_total).
        sgp4_lines = []
        self._sdp4_satrecs = []
        self._total_sats = len(pairs)

        for l1, l2 in pairs:
            sat = _Satrec.twoline2rv(l1, l2)
            if sat.is_deep_space:
                self._sdp4_satrecs.append(sat)
            else:
                sgp4_lines.append(l1)
                sgp4_lines.append(l2)

        self._n_sgp4 = len(sgp4_lines) // 2

        # Build SGP4 constellation from near-earth TLEs only
        self._zig = None
        if sgp4_lines:
            sgp4_text = "\n".join(sgp4_lines)
            self._zig = _Sgp4Constellation.from_tle_text(sgp4_text)

    @property
    def num_satellites(self):
        """Number of satellites in the constellation."""
        return self._total_sats

    @property
    def epochs(self):
        """TLE epoch for each satellite as Julian Date."""
        epochs = []
        if self._zig is not None:
            epochs.extend(self._zig.epochs)
        for sat in self._sdp4_satrecs:
            epochs.append(sat.jdsatepoch + sat.jdsatepochF)
        return epochs


def propagate(
    source, times, *, start_time=None, output="ecef", velocities=False, norad_id=None
):
    """Propagate one or more satellites to specified times.

    Computes satellite positions (and optionally velocities) at the given times
    using SGP4 propagation with SIMD acceleration and multi-threading.

    Parameters
    ----------
    source : Constellation, str, or None
        TLE data source. Can be:

        - ``Constellation`` object (most efficient for repeated calls)
        - CelesTrak group name (e.g., ``"starlink"``, ``"gps"``)
        - URL to TLE data
        - Local file path
        - Raw TLE text string
        - ``None`` if using ``norad_id``

    times : array-like
        Times at which to compute positions, in minutes from ``start_time``.
        Will be converted to a contiguous float64 array.
    start_time : datetime, optional
        Reference time for the ``times`` array. Defaults to current UTC time.
    output : {"ecef", "teme", "geodetic"}, default "ecef"
        Output coordinate frame:

        - ``"ecef"``: Earth-Centered Earth-Fixed (x, y, z in km)
        - ``"teme"``: True Equator Mean Equinox (x, y, z in km)
        - ``"geodetic"``: (latitude in deg, longitude in deg, altitude in km)

    velocities : bool, default False
        If True, also return velocity vectors.
    norad_id : int or list of int, optional
        NORAD catalog ID(s) to fetch from CelesTrak. If provided, ``source``
        is ignored.

    Returns
    -------
    positions : ndarray
        Array of shape ``(n_times, n_satellites, 3)`` containing positions.
        Units depend on ``output`` parameter.
    velocities : ndarray, optional
        Only returned if ``velocities=True``. Array of shape
        ``(n_times, n_satellites, 3)`` containing velocities in km/s.

    Examples
    --------
    Basic propagation with CelesTrak group:

    >>> positions = propagate("starlink", [0, 60, 120])
    >>> positions.shape
    (3, 6547, 3)  # 3 times, 6547 satellites, xyz

    Get positions and velocities:

    >>> pos, vel = propagate("gps", [0, 30, 60], velocities=True)
    >>> vel.shape
    (3, 31, 3)

    Output in geodetic coordinates (lat, lon, alt):

    >>> geo = propagate("stations", [0], output="geodetic")
    >>> lat, lon, alt = geo[0, 0]  # First time, first satellite

    Propagate specific satellite by NORAD ID:

    >>> iss_pos = propagate(None, range(1440), norad_id=25544)

    Efficient repeated propagation with Constellation:

    >>> c = Constellation("starlink")
    >>> pos1 = propagate(c, [0, 60])
    >>> pos2 = propagate(c, [120, 180])

    Notes
    -----
    - Uses SIMD-accelerated SGP4 (8-wide AVX2/AVX-512 when available)
    - Multi-threaded for large constellations
    - For repeated propagation of the same satellites, use a ``Constellation``
      object to avoid re-parsing TLEs

    See Also
    --------
    Constellation : Pre-parsed TLE data for efficient reuse.
    screen : Screen for conjunction events.
    """
    if isinstance(source, Constellation):
        const = source
    else:
        const = Constellation(source, norad_id=norad_id)

    times = np.ascontiguousarray(times, dtype=np.float64)
    n_sats, n_times = const.num_satellites, len(times)
    shape = (n_times, n_sats, 3)
    start = _start_jd(start_time)

    pos = np.empty(shape, dtype=np.float64)
    vel = np.empty(shape, dtype=np.float64) if velocities else None
    n_sgp4 = const._n_sgp4

    # SGP4 batch: positions [0, n_sgp4) — write directly into output
    if const._zig is not None and n_sgp4 > 0:
        sgp4_epochs = np.array(const._zig.epochs, dtype=np.float64)
        sgp4_offsets = (start - _pad_epochs_for_simd(sgp4_epochs)) * 1440.0

        # Write SGP4 results directly into output array.
        # For mixed case, output_stride=n_sats tells Zig to use the full
        # satellite dimension as stride, writing SGP4 at indices [0, n_sgp4).
        const._zig.propagate_into(
            times,
            pos,
            vel,
            epoch_offsets=sgp4_offsets,
            satellite_mask=None,
            output=output,
            reference_jd=start,
            time_major=True,
            output_stride=n_sats,
        )

    # SDP4 batch: positions [n_sgp4, n_total) — threaded in Zig
    n_sdp4 = len(const._sdp4_satrecs)
    if n_sdp4 > 0:
        abs_jd_times = start + times / 1440.0
        common_jd = np.full(
            n_times, np.floor(abs_jd_times[0] - 0.5) + 0.5, dtype=np.float64
        )
        common_fr = np.ascontiguousarray(abs_jd_times - common_jd[0], dtype=np.float64)

        # Write directly into output array using stride — no temp, no copy
        _sdp4_batch_propagate_into(
            const._sdp4_satrecs,
            common_jd,
            common_fr,
            pos,
            vel if velocities else pos,  # vel unused if no velocities
            output_stride=n_sats,
            sat_offset=n_sgp4,
        )

    return (pos, vel) if velocities else pos


def screen(
    source, times, threshold=10.0, *, target=None, start_time=None, norad_id=None
):
    """Screen a constellation for conjunction events.

    Identifies close approaches between satellites within a distance threshold.
    Supports two modes: single-target (fastest) and all-vs-all screening.

    Parameters
    ----------
    source : Constellation, str, or None
        TLE data source. Can be:

        - ``Constellation`` object (most efficient for repeated calls)
        - CelesTrak group name (e.g., ``"starlink"``, ``"gps"``)
        - URL to TLE data
        - Local file path
        - Raw TLE text string
        - ``None`` if using ``norad_id``

    times : array-like
        Times at which to check for conjunctions, in minutes from ``start_time``.
    threshold : float, default 10.0
        Distance threshold in kilometers. Pairs closer than this are reported.
    target : int, optional
        Satellite index for single-target screening. When set, only conjunctions
        with this specific satellite are found. This uses a fused propagate+screen
        algorithm and is significantly faster than all-vs-all mode.
    start_time : datetime, optional
        Reference time for the ``times`` array. Defaults to current UTC time.
    norad_id : int or list of int, optional
        NORAD catalog ID(s) to fetch from CelesTrak. If provided, ``source``
        is ignored.

    Returns
    -------
    For single-target mode (``target`` is set):
        min_distances : ndarray
            Array of shape ``(n_satellites,)`` containing the minimum distance
            to the target satellite over all times.
        min_t_indices : ndarray
            Array of shape ``(n_satellites,)`` containing the time index at
            which the minimum distance occurred.

    For all-vs-all mode (``target`` is None):
        pairs : ndarray
            Array of shape ``(n_conjunctions, 2)`` containing satellite index
            pairs that came within threshold distance.
        t_indices : ndarray
            Array of shape ``(n_conjunctions,)`` containing the time index
            at which each conjunction occurred.

    Examples
    --------
    Single-target screening (fastest):

    >>> min_dist, min_t = screen("starlink", range(1440), threshold=5.0, target=0)
    >>> # Find satellites that came within 5km of satellite 0
    >>> close = np.where(min_dist < 5.0)[0]

    All-vs-all screening:

    >>> pairs, times = screen("starlink", range(1440), threshold=10.0)
    >>> print(f"Found {len(pairs)} conjunction events")

    Screen ISS against all active satellites:

    >>> c = Constellation("active")
    >>> # Find ISS index in constellation, then screen
    >>> min_dist, min_t = screen(c, range(1440), threshold=10.0, target=iss_idx)

    Notes
    -----
    - Single-target mode uses a fused propagate+screen algorithm that is
      significantly faster than all-vs-all mode
    - All-vs-all mode scales as O(n^2) with number of satellites
    - For large constellations, consider using single-target mode with
      specific satellites of interest

    See Also
    --------
    propagate : Propagate satellites to get full position data.
    Constellation : Pre-parsed TLE data for efficient reuse.
    """
    if isinstance(source, Constellation):
        const = source
    else:
        const = Constellation(source, norad_id=norad_id)

    times_arr = np.ascontiguousarray(times, dtype=np.float64)

    # Pure SGP4 constellation — use fast native paths
    if not const._sdp4_satrecs and const._zig is not None:
        constellation = const._zig
        epoch_offsets, start_jd = _compute_epoch_offsets(constellation, start_time)

        if target is not None:
            result = constellation.screen_conjunction(
                times_arr,
                target,
                threshold,
                epoch_offsets=epoch_offsets,
                reference_jd=start_jd,
            )
            return np.array(result[0]), np.array(result[1], dtype=np.uint32)
        else:
            n_sats, n_times = constellation.num_satellites, len(times_arr)
            pos = np.empty((n_times, n_sats, 3), dtype=np.float64)
            constellation.propagate_into(
                times_arr,
                pos,
                None,
                epoch_offsets=epoch_offsets,
                satellite_mask=None,
                output="teme",
                reference_jd=start_jd,
                time_major=True,
            )
            return _coarse_screen(pos, n_sats, threshold, None)

    # Mixed SGP4/SDP4 — propagate all then screen
    pos = propagate(const, times_arr, start_time=start_time, output="teme")
    n_sats = const.num_satellites
    return _coarse_screen(pos, n_sats, threshold, None)


__version__ = "0.7.1"

__all__ = [
    "__version__",
    "Tle",
    "Constellation",
    "propagate",
    "screen",
    "hohmann_transfer",
    "bi_elliptic_transfer",
    "lambert",
    "orbital_velocity",
    "orbital_period",
    "escape_velocity",
    "propagate_numerical",
    "EARTH_MU",
    "EARTH_R_EQ",
    "EARTH_J2",
    "SUN_MU",
    "MOON_MU",
]
