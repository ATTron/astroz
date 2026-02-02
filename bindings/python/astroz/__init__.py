"""astroz - High-performance astrodynamics library.

A Python library for satellite orbit propagation and conjunction screening,
powered by a Zig backend with SIMD acceleration (AVX2/AVX512).

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
    Sgp4Constellation as _Sgp4Constellation,
    coarse_screen as _coarse_screen,
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
    if Path(source).exists():
        return Path(source).read_text()
    if "1 " in source and "2 " in source:
        return source
    return _fetch_url(_celestrak_url(group=source))


def _get_constellation(source, norad_id=None):
    """Get _Sgp4Constellation from source (Constellation object or TLE source)."""
    if isinstance(source, Constellation):
        return source._zig
    tle_text = _load_tle_text(source, norad_id)
    return _Sgp4Constellation.from_tle_text(tle_text)


def _compute_epoch_offsets(constellation, start_time):
    """Compute epoch offsets from start_time."""
    if start_time is None:
        start_time = datetime.now(timezone.utc)
    start_jd = 2440587.5 + (start_time.timestamp() / 86400.0)
    epochs = np.array(constellation.epochs, dtype=np.float64)
    return (start_jd - epochs) * 1440.0, start_jd


class Constellation:
    """Pre-parsed TLE data for efficient repeated propagation and screening.

    Use this class when you need to propagate or screen the same satellites
    multiple times. It parses TLE data once and caches the result, avoiding
    the overhead of re-parsing on every call.

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
        self._zig = _Sgp4Constellation.from_tle_text(tle_text)

    @property
    def num_satellites(self):
        """Number of satellites in the constellation."""
        return self._zig.num_satellites

    @property
    def epochs(self):
        """TLE epoch for each satellite as Julian Date."""
        return self._zig.epochs


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
    constellation = _get_constellation(source, norad_id)

    times = np.ascontiguousarray(times, dtype=np.float64)
    n_sats, n_times = constellation.num_satellites, len(times)
    shape = (n_times, n_sats, 3)

    epoch_offsets, start_jd = _compute_epoch_offsets(constellation, start_time)

    pos = np.empty(shape, dtype=np.float64)
    vel = np.empty(shape, dtype=np.float64) if velocities else None

    constellation.propagate_into(
        times,
        pos,
        vel,
        epoch_offsets=epoch_offsets,
        satellite_mask=None,
        output=output,
        reference_jd=start_jd,
        time_major=True,
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
    constellation = _get_constellation(source, norad_id)

    times = np.ascontiguousarray(times, dtype=np.float64)
    epoch_offsets, start_jd = _compute_epoch_offsets(constellation, start_time)

    if target is not None:
        # Single target - use fused propagate+screen (fastest path)
        result = constellation.screen_conjunction(
            times,
            target,
            threshold,
            epoch_offsets=epoch_offsets,
            reference_jd=start_jd,
        )
        return np.array(result[0]), np.array(result[1], dtype=np.uint32)
    else:
        # All-vs-all - propagate then screen
        n_sats, n_times = constellation.num_satellites, len(times)
        shape = (n_times, n_sats, 3)
        pos = np.empty(shape, dtype=np.float64)

        constellation.propagate_into(
            times,
            pos,
            None,
            epoch_offsets=epoch_offsets,
            satellite_mask=None,
            output="teme",
            reference_jd=start_jd,
            time_major=True,
        )
        return _coarse_screen(pos, n_sats, threshold, None)


__version__ = "0.4.4"

__all__ = [
    "__version__",
    "Tle",
    "Constellation",
    "propagate",
    "screen",
]
