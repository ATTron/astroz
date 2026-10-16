"""Constellation loading: OMM JSON is read directly (not rewritten as TLE
text first), and every record ends up in exactly one place.

The records are CelesTrak data from the gpconf corpus
(https://github.com/hneogy/gp-omm-conformance). SARAMAGO is catalog number
100000, which a TLE cannot hold, and its eccentricity and BSTAR carry more
digits than TLE columns do.
"""

import json
import math
from datetime import datetime, timezone

import numpy as np
import pytest

import astroz
from astroz._astroz import Satrec

SARAMAGO = {
    "OBJECT_NAME": "SARAMAGO", "OBJECT_ID": "2026-067CY", "EPOCH": "2026-07-14T21:45:20.933856",
    "MEAN_MOTION": 15.20467281, "ECCENTRICITY": 0.00055903, "INCLINATION": 97.4593,
    "RA_OF_ASC_NODE": 154.097, "ARG_OF_PERICENTER": 270.5113, "MEAN_ANOMALY": 89.5482,
    "EPHEMERIS_TYPE": 0, "CLASSIFICATION_TYPE": "U", "NORAD_CAT_ID": 100000,
    "ELEMENT_SET_NO": 999, "REV_AT_EPOCH": 1591, "BSTAR": 0.00022159168,
    "MEAN_MOTION_DOT": 4.77e-5, "MEAN_MOTION_DDOT": 0,
}

OBJ_81011 = {
    "OBJECT_NAME": "UNKNOWN", "OBJECT_ID": "", "EPOCH": "2026-04-08T23:58:19.197120",
    "MEAN_MOTION": 14.07928696, "ECCENTRICITY": 0.0134684, "INCLINATION": 99.3661,
    "RA_OF_ASC_NODE": 118.128, "ARG_OF_PERICENTER": 184.7263, "MEAN_ANOMALY": 194.7833,
    "EPHEMERIS_TYPE": 0, "CLASSIFICATION_TYPE": "U", "NORAD_CAT_ID": 81011,
    "ELEMENT_SET_NO": 999, "REV_AT_EPOCH": 1, "BSTAR": 0.0012275,
    "MEAN_MOTION_DOT": 2.218e-5, "MEAN_MOTION_DDOT": 0,
}
TLE_81011 = (
    "1 81011U          26098.99883330  .00002218  00000+0  12275-2 0  9990\n"
    "2 81011  99.3661 118.1280 0134684 184.7263 194.7833 14.07928696    13"
)


def test_six_digit_catalog_number():
    (sat,) = Satrec.from_omm(json.dumps(SARAMAGO))
    assert sat.satnum == 100000


def test_values_keep_full_precision():
    (sat,) = Satrec.from_omm(json.dumps(SARAMAGO))
    assert sat.ecco == 0.00055903
    assert sat.bstar == 0.00022159168


def test_single_object_and_array_agree():
    (a,) = Satrec.from_omm(json.dumps(SARAMAGO))
    (b,) = Satrec.from_omm(json.dumps([SARAMAGO]))
    assert a.sgp4(a.jdsatepoch, a.jdsatepochF) == b.sgp4(b.jdsatepoch, b.jdsatepochF)


def test_constellation_from_mixed_catalog_numbers():
    c = astroz.Constellation(json.dumps([SARAMAGO, OBJ_81011]))
    assert c.num_satellites == 2
    start = datetime(2026, 7, 14, 21, 45, 21, tzinfo=timezone.utc)
    pos = astroz.propagate(c, [0.0, 60.0], start_time=start, output="teme")
    assert pos.shape == (2, 2, 3)
    assert np.isfinite(pos[:, 0]).all()  # SARAMAGO, near its epoch


def test_omm_matches_tle_for_same_object():
    start = datetime(2026, 4, 9, tzinfo=timezone.utc)
    from_omm = astroz.propagate(json.dumps(OBJ_81011), [0.0, 720.0], start_time=start, output="teme")
    from_tle = astroz.propagate(TLE_81011, [0.0, 720.0], start_time=start, output="teme")
    assert np.abs(from_omm - from_tle).max() < 1e-2  # km; the OMM values match the TLE columns


def test_invalid_omm_raises():
    with pytest.raises(ValueError, match="OMM"):
        astroz.Constellation(json.dumps([{"NORAD_CAT_ID": 1}]))


NEAR_EPOCH = datetime(2026, 7, 14, 21, 45, 21, tzinfo=timezone.utc)  # SARAMAGO's epoch


def test_period_near_deep_space_boundary_loads():
    # Mean motion 6.40-6.41 rev/day puts the period right at SGP4's 225-minute
    # deep-space cutoff, where a raw and a corrected period disagree.
    recs = [SARAMAGO] + [dict(SARAMAGO, NORAD_CAT_ID=200000 + i, MEAN_MOTION=6.4 + i * 0.0005, ECCENTRICITY=0.001)
                         for i in range(21)]
    c = astroz.Constellation(json.dumps(recs))
    assert c.num_satellites == len(recs)
    pos = astroz.propagate(c, [0.0], start_time=NEAR_EPOCH, output="teme")
    assert np.isfinite(pos).all()


def test_rejected_record_is_nan_and_does_not_break_the_rest():
    rejected = dict(SARAMAGO, NORAD_CAT_ID=100001, ECCENTRICITY=1.2)  # not an ellipse: SGP4 error 1
    c = astroz.Constellation(json.dumps([SARAMAGO, rejected]))
    assert c.num_satellites == 2
    assert len(c.epochs) == 2
    assert c.satnums == [100000, 100001]  # good one first, rejected one last
    pos, vel = astroz.propagate(c, [0.0, 10.0], start_time=NEAR_EPOCH, output="teme", velocities=True)
    assert np.isfinite(pos[:, 0]).all() and np.isfinite(vel[:, 0]).all()
    assert np.isnan(pos[:, 1]).all() and np.isnan(vel[:, 1]).all()
    events = astroz.screen(c, 10.0, threshold=10.0, start_time=NEAR_EPOCH)
    assert not np.isin(1, events.i) and not np.isin(1, events.j)


def test_propagate_matches_python_sgp4_wgs72():
    sgp4_api = pytest.importorskip("sgp4.api")
    l1, l2 = TLE_81011.splitlines()
    ref = sgp4_api.Satrec.twoline2rv(l1, l2, sgp4_api.WGS72)
    start = datetime(2026, 4, 9, tzinfo=timezone.utc)
    jd, fr = sgp4_api.jday(2026, 4, 9, 0, 0, 0.0)
    for minutes in (0.0, 720.0):
        err, r, _ = ref.sgp4(jd, fr + minutes / 1440.0)
        assert err == 0
        got = astroz.propagate(TLE_81011, [minutes], start_time=start, output="teme")[0, 0]
        assert np.abs(got - np.array(r)).max() < 1e-3  # km


def _shell(n, seed):
    """Random near-circular orbits in a thin shell around 550 km, so plenty of pairs pass close."""
    rng = np.random.default_rng(seed)
    return [dict(SARAMAGO, NORAD_CAT_ID=300000 + k, MEAN_MOTION=15.05 + rng.normal(0, 0.01),
                 ECCENTRICITY=float(rng.uniform(1e-4, 2e-3)), INCLINATION=float(rng.uniform(40, 100)),
                 RA_OF_ASC_NODE=float(rng.uniform(0, 360)), ARG_OF_PERICENTER=float(rng.uniform(0, 360)),
                 MEAN_ANOMALY=float(rng.uniform(0, 360)), BSTAR=1e-5, MEAN_MOTION_DOT=0)
            for k in range(n)]


def _brute_force(c, start, minutes, threshold):
    """Every closest approach within threshold: scan each pair at 1 s, then golden-section
    refine each local minimum with the scalar propagator. Shares no code with screen()."""
    sats = c._satrecs
    jd0 = astroz._start_jd(start)
    jd_day, fr0 = float(math.floor(jd0)), jd0 - math.floor(jd0)
    t = np.arange(0, minutes * 60 + 1, 1.0)
    pos = np.empty((len(sats), len(t), 3))
    for s, sat in enumerate(sats):
        sat.sgp4_array_into(np.full(len(t), jd_day), fr0 + t / 86400, pos[s], np.empty((len(t), 3)))

    def dist(i, j, sec):
        ri, rj, v = np.empty((1, 3)), np.empty((1, 3)), np.empty((1, 3))
        jd, fr = np.array([jd_day]), np.array([fr0 + sec / 86400])
        sats[i].sgp4_array_into(jd, fr, ri, v)
        sats[j].sgp4_array_into(jd, fr, rj, v)
        return float(np.linalg.norm(ri - rj))

    events = {}
    golden = (math.sqrt(5) - 1) / 2
    for i in range(len(sats) - 1):
        d = np.linalg.norm(pos[i][None] - pos[i + 1:], axis=2)
        interior = (d[:, 1:-1] < d[:, :-2]) & (d[:, 1:-1] <= d[:, 2:]) & (d[:, 1:-1] < threshold + 10)
        for jj, k in zip(*np.nonzero(interior)):
            j, lo, hi = i + 1 + jj, t[k], t[k + 2]
            for _ in range(60):
                m1, m2 = hi - golden * (hi - lo), lo + golden * (hi - lo)
                lo, hi = (lo, m2) if dist(i, j, m1) < dist(i, j, m2) else (m1, hi)
            sec = 0.5 * (lo + hi)
            if dist(i, j, sec) <= threshold:
                events[(i, j, round(sec))] = dist(i, j, sec)
    return events


@pytest.mark.parametrize("step", [1.0, 10.0])
def test_screen_finds_every_encounter(step, monkeypatch):
    c = astroz.Constellation(json.dumps(_shell(120, 7)))
    start = datetime(2026, 7, 14, 22, 0, tzinfo=timezone.utc)
    want = _brute_force(c, start, 120, 25.0)
    assert len(want) > 5
    # Tiny chunks, so encounters also land on chunk boundaries
    monkeypatch.setattr(astroz, "_SCREEN_STATE_BUDGET", 120 * 48 * 3)
    ev = astroz.screen(c, 120, threshold=25.0, step=step, start_time=start)
    got = {(int(i), int(j), round(t * 60)): d for i, j, t, d in zip(ev.i, ev.j, ev.time, ev.distance)}
    assert got.keys() == want.keys()
    for k in want:
        # Refinement evaluates several times at once through the SIMD path, which
        # differs from the scalar path by centimetres; astroz allows 1 m for that.
        assert abs(got[k] - want[k]) < 1e-3  # km
    assert np.all(np.diff(ev.time) >= 0) and np.all(ev.speed > 0)


def test_screen_target_limits_pairs():
    c = astroz.Constellation(json.dumps(_shell(120, 7)))
    start = datetime(2026, 7, 14, 22, 0, tzinfo=timezone.utc)
    every = astroz.screen(c, 120, threshold=25.0, start_time=start)
    target = int(every.i[0])
    some = astroz.screen(c, 120, threshold=25.0, target=target, start_time=start)
    involved = (every.i == target) | (every.j == target)
    assert len(some.i) == involved.sum() > 0
    assert np.all((some.i == target) | (some.j == target))


def test_screen_stays_inside_the_window():
    c = astroz.Constellation(json.dumps(_shell(120, 7)))
    start = datetime(2026, 7, 14, 22, 0, tzinfo=timezone.utc)
    every = astroz.screen(c, 120, threshold=25.0, start_time=start)
    cut = float(np.median(every.time))
    # 7-minute steps don't divide the window, so the grid runs past it
    part = astroz.screen(c, cut, threshold=25.0, step=7.0, start_time=start)
    assert len(part.i) == np.sum(every.time <= cut) > 0
    assert np.all(part.time <= cut)
