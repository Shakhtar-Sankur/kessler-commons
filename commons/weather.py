"""Geomagnetic activity, for telling maneuvers apart from storm drag.

A geomagnetic storm heats the upper atmosphere, so low satellites lose altitude faster for a day
or two. That extra drop can look like a small deorbit burn, or it can hide a real one. The planetary
Kp index (NOAA SWPC, every 3 hours, 0-9) says when a storm is on; the NOAA G-scale labels its size
(Kp 5 = G1 ... Kp 9 = G5). This module turns the Kp feed into a lookup over the archive's time axis,
so every detected orbit change can be tagged with the geomagnetic conditions around it.
"""
import json

import numpy as np

from .orbits import epoch_days

KP_URL = "https://services.swpc.noaa.gov/products/noaa-planetary-k-index.json"


def load_kp(records):
    """(epoch days, Kp) from NOAA planetary-k-index records, sorted; non-finite Kp dropped."""
    out = []
    for r in records:
        try:
            out.append((epoch_days(r["time_tag"]), float(r["Kp"])))
        except (KeyError, ValueError, TypeError):
            continue
    out.sort()
    return np.array([t for t, _ in out]), np.array([k for _, k in out])


def load_kp_file(path):
    return load_kp(json.load(open(path)))


def g_scale(kp):
    """NOAA geomagnetic storm scale: 0 below G1, else G1..G5 for Kp 5..9."""
    return 0 if kp < 5 else min(5, int(kp) - 4)


def max_kp(t_kp, kp, t0, t1, pad_days=0.25):
    """Highest Kp seen while a change happened: samples within [t0, t1] widened by pad on each
    side (Kp is 3-hourly, so a short interval may contain no sample). None if the feed does not
    cover the interval at all."""
    if len(t_kp) == 0:
        return None
    lo, hi = min(t0, t1) - pad_days, max(t0, t1) + pad_days
    m = (t_kp >= lo) & (t_kp <= hi)
    if m.any():
        return float(kp[m].max())
    i = int(np.argmin(np.abs(t_kp - 0.5 * (t0 + t1))))
    return float(kp[i]) if abs(t_kp[i] - 0.5 * (t0 + t1)) <= 0.5 else None
