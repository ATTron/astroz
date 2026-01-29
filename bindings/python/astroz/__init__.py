"""astroz - High-performance astrodynamics library."""

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
    """Parsed TLE data for repeated propagation/screening.

    Use this to avoid re-parsing TLEs on every call:
        c = Constellation("starlink")
        pos1 = propagate(c, times1)
        pos2 = propagate(c, times2)
    """

    def __init__(self, source=None, *, norad_id=None):
        tle_text = _load_tle_text(source, norad_id)
        self._zig = _Sgp4Constellation.from_tle_text(tle_text)

    @property
    def num_satellites(self):
        return self._zig.num_satellites

    @property
    def epochs(self):
        return self._zig.epochs


def propagate(
    source, times, *, start_time=None, output="ecef", velocities=False, norad_id=None
):
    """Propagate one or more satellites.

    Args:
        source: Constellation, TLE text string, file path, URL, or CelesTrak group name
        times: array of minutes from start_time
        start_time: datetime (defaults to now UTC)
        output: "ecef" (default), "teme", or "geodetic"
        velocities: if True, also return velocities
        norad_id: NORAD catalog ID or list of IDs (alternative to source)

    Returns:
        positions: array of shape (n_times, n_sats, 3)
        If velocities=True: (positions, velocities) tuple
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
    """Screen constellation for conjunctions.

    Args:
        source: Constellation, TLE text string, file path, URL, or CelesTrak group name
        times: array of minutes from start_time
        threshold: distance threshold in km
        target: satellite index for single-target screening (fastest), or None for all-vs-all
        start_time: datetime (defaults to now UTC)
        norad_id: NORAD catalog ID or list of IDs (alternative to source)

    Returns:
        If target is set (single-target, fused - fastest):
            (min_distances, min_t_indices) - per satellite minimum distance to target
        If target is None (all-vs-all):
            (pairs, t_indices) - all conjunction events
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
