"""Orbital Watch: every detectable orbit change of every satellite in the archive, by operator.

Reads the archive of operator-derived element sets (CelesTrak SupGP snapshots, one file per
constellation every 6 hours: <root>/YYYY/MM/DD/supgp-<group>_HHMMZ.csv.gz), joins each
satellite's element sets across snapshots into one history, finds the steps in its mean
semi-major axis (commons.orbits.detect), and summarises them per constellation.

    python -m commons.watch --archive ../kessler-archive/data --out out

Outputs (out/, never committed until publication is cleared with the data providers):
  maneuvers.csv   one row per detected orbit change
  summary.json    per constellation: satellites, satellite-days observed, changes, rates, sizes
"""
import argparse
import csv
import gzip
import json
import os
import re
from collections import defaultdict

import numpy as np

from .orbits import detect, history
from .weather import g_scale, load_kp_file, max_kp

STORM_KP = 5.0   # NOAA G1 and above

FILE = re.compile(r"supgp-(?P<group>[a-z0-9-]+)_(?P<hhmm>\d{4})Z\.csv\.gz$")


def snapshots(root):
    """[(group, path)] of every SupGP snapshot under root, oldest first."""
    found = []
    for dirpath, _, files in os.walk(root):
        for f in files:
            m = FILE.match(f)
            if m:
                found.append((os.path.relpath(os.path.join(dirpath, f), root), m.group("group")))
    return [(g, os.path.join(root, p)) for p, g in sorted(found)]


def element_sets(root):
    """{group: {norad: [rows]}} across all snapshots; identical element sets appear once."""
    out = defaultdict(lambda: defaultdict(dict))
    for group, path in snapshots(root):
        with gzip.open(path, "rt") as f:
            for r in csv.DictReader(f):
                out[group][r["NORAD_CAT_ID"]][r["EPOCH"]] = r
    return {g: {n: list(v.values()) for n, v in sats.items()} for g, sats in out.items()}


def watch(sets, min_sets=6, kp=None, **kw):
    """Detected orbit changes and per-group statistics. With kp=(epoch days, Kp), each change is
    tagged with the highest Kp around it, and exposure and change rates are split into storm
    (Kp >= 5) and quiet time."""
    t_kp, v_kp = kp if kp is not None else (np.array([]), np.array([]))
    changes, summary = [], {}
    for group, sats in sorted(sets.items()):
        n_sat, days, sizes = 0, 0.0, []
        exp = {"storm": 0.0, "quiet": 0.0}
        cnt = {"storm": [], "quiet": []}
        for norad, rows in sats.items():
            t, a = history(rows)
            if len(t) < min_sets:
                continue
            n_sat += 1
            days += float(t[-1] - t[0])
            name = rows[0].get("OBJECT_NAME", "")
            for i in range(1, len(t)):          # exposure: time between element sets, by weather
                m = max_kp(t_kp, v_kp, t[i - 1], t[i])
                if m is not None:
                    exp["storm" if m >= STORM_KP else "quiet"] += float(t[i] - t[i - 1])
            for k, s in zip(*detect(t, a, **kw)):
                m = max_kp(t_kp, v_kp, t[k - 1], t[k])
                changes.append({"group": group, "norad": norad, "name": name,
                                "between": [round(float(t[k - 1]), 5), round(float(t[k]), 5)],
                                "delta_a_m": round(float(s), 1), "sma_km": round(float(a[k]), 3),
                                "max_kp": m, "g": None if m is None else g_scale(m)})
                sizes.append(float(s))
                if m is not None:
                    cnt["storm" if m >= STORM_KP else "quiet"].append(float(s))
        sizes = np.array(sizes)
        summary[group] = {
            "satellites": n_sat,
            "satellite_days": round(days, 1),
            "changes": int(len(sizes)),
            "per_satellite_day": round(len(sizes) / days, 4) if days else None,
            "raises": int((sizes > 0).sum()),
            "lowerings": int((sizes < 0).sum()),
            "median_abs_m": round(float(np.median(np.abs(sizes))), 1) if len(sizes) else None,
        }
        if len(t_kp):
            summary[group]["weather"] = {
                w: {"satellite_days": round(exp[w], 1), "changes": len(cnt[w]),
                    "per_satellite_day": round(len(cnt[w]) / exp[w], 4) if exp[w] else None,
                    "lowerings": int(sum(x < 0 for x in cnt[w]))}
                for w in ("quiet", "storm")}
    return changes, summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", required=True)
    ap.add_argument("--out", default="out")
    ap.add_argument("--min-sets", type=int, default=6)
    ap.add_argument("--kp", help="NOAA planetary-k-index.json (tags changes with geomagnetic activity)")
    a = ap.parse_args()
    changes, summary = watch(element_sets(a.archive), a.min_sets,
                             kp=load_kp_file(a.kp) if a.kp else None)
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "maneuvers.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["group", "norad", "name", "after_epoch_days", "before_epoch_days", "delta_a_m", "sma_km",
                    "max_kp", "g_scale"])
        for c in changes:
            w.writerow([c["group"], c["norad"], c["name"], *c["between"], c["delta_a_m"], c["sma_km"],
                        c["max_kp"], c["g"]])
    json.dump(summary, open(os.path.join(a.out, "summary.json"), "w"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
