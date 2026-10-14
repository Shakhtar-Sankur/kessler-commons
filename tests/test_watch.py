"""Orbital Watch on a synthetic archive with known maneuvers."""
import csv
import gzip
import math
import os
from datetime import datetime, timedelta

import numpy as np

from commons.orbits import MU, brouwer_sma_km, detect
from commons.watch import element_sets, watch


def mean_motion_for(a_km, ecc=0.0001, incl_deg=53.0):
    """Kozai mean motion (rev/day) whose Brouwer mean semi-major axis is a_km (by bisection)."""
    lo, hi = 10.0, 17.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if brouwer_sma_km(mid / 720.0 * math.pi, ecc, math.radians(incl_deg)) > a_km:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def write_archive(root, sats, start, n_snap=20, every_h=6.0, noise_m=5.0, seed=0):
    rng = np.random.default_rng(seed)
    fields = ["OBJECT_NAME", "NORAD_CAT_ID", "EPOCH", "MEAN_MOTION", "ECCENTRICITY", "INCLINATION"]
    for s in range(n_snap):
        when = start + timedelta(hours=every_h * s)
        d = os.path.join(root, when.strftime("%Y/%m/%d"))
        os.makedirs(d, exist_ok=True)
        with gzip.open(os.path.join(d, f"supgp-test_{when:%H%M}Z.csv.gz"), "wt", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            for norad, (a0, drag_m_per_day, jumps) in sats.items():
                epoch = when - timedelta(minutes=30)
                days = (epoch - start).total_seconds() / 86400.0
                a = a0 * 1e3 - drag_m_per_day * days + sum(dm for td, dm in jumps if days > td)
                a += rng.normal(0, noise_m)
                w.writerow({"OBJECT_NAME": f"SAT-{norad}", "NORAD_CAT_ID": norad,
                            "EPOCH": epoch.strftime("%Y-%m-%dT%H:%M:%S.%f"),
                            "MEAN_MOTION": f"{mean_motion_for(a / 1e3):.10f}",
                            "ECCENTRICITY": "0.0001", "INCLINATION": "53.0"})


def test_finds_the_injected_maneuvers_and_nothing_else(tmp_path):
    start = datetime(2026, 10, 10)
    sats = {
        "1": (6928.0, 20.0, [(2.1, 300.0)]),          # one raise of 300 m on day 2.1
        "2": (6928.0, 20.0, [(1.6, -150.0), (3.6, 200.0)]),
        "3": (6928.0, 20.0, []),                       # drag only
        "4": (6928.0, 60.0, []),                       # strong drag only
    }
    write_archive(str(tmp_path), sats, start)
    sets = element_sets(str(tmp_path))
    assert set(sets["test"]) == set(sats)
    changes, summary = watch(sets)
    found = sorted((c["norad"], round(c["delta_a_m"], -1)) for c in changes)
    assert [n for n, _ in found] == ["1", "2", "2"]
    for c in changes:
        true = {"1": [300.0], "2": [-150.0, 200.0]}[c["norad"]]
        assert min(abs(c["delta_a_m"] - x) for x in true) < 40
    s = summary["test"]
    assert s["satellites"] == 4 and s["changes"] == 3 and s["raises"] == 2 and s["lowerings"] == 1


def test_duplicate_element_sets_across_snapshots_count_once(tmp_path):
    start = datetime(2026, 10, 10)
    write_archive(str(tmp_path), {"9": (6928.0, 20.0, [])}, start, n_snap=8)
    # the same snapshot again under a later file name: identical element sets
    src = next(p for p in tmp_path.rglob("*.csv.gz"))
    dup = src.with_name("supgp-test_2359Z.csv.gz")
    dup.write_bytes(src.read_bytes())
    sets = element_sets(str(tmp_path))
    assert len(sets["test"]["9"]) == 8


def test_detect_needs_enough_points():
    k, s = detect(np.array([0.0, 0.25, 0.5]), np.array([6928.0, 6928.0, 6928.3]))
    assert len(k) == 0
