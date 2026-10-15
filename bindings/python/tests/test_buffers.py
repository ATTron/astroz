"""Native functions only accept C-contiguous float64 arrays, and check their sizes.

Before, an int64 or float32 array was read as raw float64 bytes and the result
was silently garbage, and short output arrays were written past their end.
"""

import numpy as np
import pytest

from astroz._astroz import Satrec, SatrecArray, Sgp4Constellation, WGS72, sdp4_batch_propagate_into

ISS = (
    "1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995",
    "2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123",
)
GEO = (
    "1 28626U 05008A   06176.46683397 -.00000205  00000-0  10000-3 0  2190",
    "2 28626   0.0019 286.9433 0000335  13.7918  55.6504  1.00270176  4891",
)


def _times(n=4):
    sat = Satrec.twoline2rv(*ISS, WGS72)
    return sat, np.full(n, sat.jdsatepoch), sat.jdsatepochF + np.arange(n) / 1440.0


def test_sgp4_array_into_float64_works():
    sat, jd, fr = _times()
    r, v = np.empty((4, 3)), np.empty((4, 3))
    sat.sgp4_array_into(jd, fr, r, v)
    assert np.isfinite(r).all() and 6000 < np.linalg.norm(r[0]) < 7000


@pytest.mark.parametrize("dtype", [np.int64, np.float32])
def test_sgp4_array_into_rejects_other_dtypes(dtype):
    sat, jd, fr = _times()
    r, v = np.empty((4, 3)), np.empty((4, 3))
    with pytest.raises(TypeError, match="float64"):
        sat.sgp4_array_into(jd.astype(dtype), fr, r, v)
    with pytest.raises(TypeError, match="float64"):
        sat.sgp4_array_into(jd, fr, r.astype(dtype), v)


def test_sgp4_array_into_rejects_non_contiguous_and_short_arrays():
    sat, jd, fr = _times(8)
    r, v = np.empty((8, 3)), np.empty((8, 3))
    with pytest.raises((TypeError, BufferError, ValueError)):
        sat.sgp4_array_into(jd[::2], fr[::2], r, v)
    with pytest.raises(ValueError):
        sat.sgp4_array_into(jd, fr[:4], r, v)
    with pytest.raises(ValueError):
        sat.sgp4_array_into(jd, fr, np.empty((4, 3)), v)


def test_batch_propagate_into_rejects_other_dtypes():
    sat = Satrec.twoline2rv(*ISS, WGS72)
    times = np.arange(5, dtype=np.float64)
    for batch in (SatrecArray([sat]), Sgp4Constellation([sat], gravity_model=WGS72)):
        pos = np.empty((5, 1, 3))
        batch.propagate_into(times, pos, None)
        assert np.isfinite(pos).all()
        with pytest.raises(TypeError, match="float64"):
            batch.propagate_into(np.arange(5), pos, None)
        with pytest.raises(TypeError, match="float64"):
            batch.propagate_into(times, pos.astype(np.float32), None)
        with pytest.raises(TypeError, match="float64"):
            batch.propagate_into(times, pos, None, epoch_offsets=np.zeros(1, dtype=np.int64))


def test_satellite_mask_errors_raise():
    c = Sgp4Constellation([Satrec.twoline2rv(*ISS, WGS72)] * 3, gravity_model=WGS72)
    pos = np.empty((2, 3, 3))
    with pytest.raises(ValueError, match="satellite_mask"):
        c.propagate_into(np.zeros(2), pos, None, satellite_mask=np.ones(1, dtype=bool))
    with pytest.raises(TypeError, match="satellite_mask"):
        c.propagate_into(np.zeros(2), pos, None, satellite_mask=np.ones(3, dtype=np.float64))


def test_sdp4_batch_propagate_into_checks_arrays():
    sat = Satrec.twoline2rv(*GEO, WGS72)
    jd = np.full(3, sat.jdsatepoch)
    fr = sat.jdsatepochF + np.arange(3) / 24.0
    pos, vel = np.empty((1, 3, 3)), np.empty((1, 3, 3))
    sdp4_batch_propagate_into([sat], jd, fr, pos, vel)
    assert np.isfinite(pos).all()
    with pytest.raises(TypeError, match="float64"):
        sdp4_batch_propagate_into([sat], jd.astype(np.int64), fr, pos, vel)
    with pytest.raises(ValueError):
        sdp4_batch_propagate_into([sat], jd, fr, np.empty((1, 2, 3)), vel)
    with pytest.raises(ValueError):
        # time-major output with the column past the stride
        sdp4_batch_propagate_into([sat], jd, fr, np.empty((3, 2, 3)), np.empty((3, 2, 3)), output_stride=2, sat_offset=2)


def test_api_accepts_scalar_jd_and_sliced_arrays():
    sgp4_api = pytest.importorskip("sgp4.api")
    from astroz.api import Satrec as ApiSatrec, SatrecArray as ApiSatrecArray

    sats = [ApiSatrec.twoline2rv(*ISS, WGS72), ApiSatrec.twoline2rv(*GEO, WGS72)]
    ref = sgp4_api.SatrecArray([sgp4_api.Satrec.twoline2rv(*t, sgp4_api.WGS72) for t in (ISS, GEO)])
    jd = 2460437.0
    fr = np.linspace(0.0, 1.0, 10)[::2]  # non-contiguous
    e, r, _ = ApiSatrecArray(sats).sgp4(jd, fr)
    re, rr, _ = ref.sgp4(np.full(len(fr), jd), np.ascontiguousarray(fr))
    ok = (e == 0) & (re == 0)
    assert ok.any()
    assert np.abs(r[ok] - rr[ok]).max() < 1e-3
