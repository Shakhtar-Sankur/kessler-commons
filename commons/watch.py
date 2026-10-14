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


def watch(sets, min_sets=6, **kw):
    """Detected orbit changes and per-group statistics."""
    changes, summary = [], {}
    for group, sats in sorted(sets.items()):
        n_sat, days, sizes = 0, 0.0, []
        for norad, rows in sats.items():
            t, a = history(rows)
            if len(t) < min_sets:
                continue
            n_sat += 1
            days += float(t[-1] - t[0])
            name = rows[0].get("OBJECT_NAME", "")
            for k, s in zip(*detect(t, a, **kw)):
                changes.append({"group": group, "norad": norad, "name": name,
                                "between": [round(float(t[k - 1]), 5), round(float(t[k]), 5)],
                                "delta_a_m": round(float(s), 1), "sma_km": round(float(a[k]), 3)})
                sizes.append(float(s))
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
    return changes, summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", required=True)
    ap.add_argument("--out", default="out")
    ap.add_argument("--min-sets", type=int, default=6)
    a = ap.parse_args()
    changes, summary = watch(element_sets(a.archive), a.min_sets)
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "maneuvers.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["group", "norad", "name", "after_epoch_days", "before_epoch_days", "delta_a_m", "sma_km"])
        for c in changes:
            w.writerow([c["group"], c["norad"], c["name"], *c["between"], c["delta_a_m"], c["sma_km"]])
    json.dump(summary, open(os.path.join(a.out, "summary.json"), "w"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
