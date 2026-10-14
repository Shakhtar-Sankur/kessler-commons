"""Kp lookup and storm tagging."""
import numpy as np

from commons.orbits import epoch_days
from commons.weather import g_scale, load_kp, max_kp


def test_g_scale_matches_noaa():
    assert [g_scale(k) for k in (0, 4.67, 5, 5.67, 6, 7, 8.33, 9)] == [0, 0, 1, 1, 2, 3, 4, 5]


def test_max_kp_over_an_interval():
    recs = [{"time_tag": f"2026-10-10T{h:02d}:00:00", "Kp": k}
            for h, k in [(0, 1), (3, 2), (6, 7), (9, 3), (12, 1)]]
    recs.append({"time_tag": "2026-10-10T15:00:00", "Kp": None})   # dropped
    t, kp = load_kp(recs)
    assert len(t) == 5
    t0 = epoch_days("2026-10-10T05:00:00")
    assert max_kp(t, kp, t0, t0 + 2 / 24, pad_days=0) == 7      # 06:00 inside
    t1 = epoch_days("2026-10-10T13:00:00")
    assert max_kp(t, kp, t1, t1 + 0.5 / 24, pad_days=0) == 1    # nearest sample
    assert max_kp(t, kp, t1 + 5, t1 + 6) is None                  # outside the feed
    assert max_kp(np.array([]), np.array([]), t0, t1) is None
