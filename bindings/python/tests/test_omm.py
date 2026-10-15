"""OMM JSON is read directly, not rewritten as TLE text first.

The records are CelesTrak data from the gpconf corpus
(https://github.com/hneogy/gp-omm-conformance). SARAMAGO is catalog number
100000, which a TLE cannot hold, and its eccentricity and BSTAR carry more
digits than TLE columns do.
"""

import json
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
