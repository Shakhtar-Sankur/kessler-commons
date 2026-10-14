# kessler Commons

An open, independent record of what satellites do in orbit, built only on open data.

The core is **Orbital Watch**, which tracks every detectable orbit change of every satellite, for
every major operator, from public element sets. It is meant to become the reference that
operators' own maneuver reports can be checked against. The tools are civilian and for space
safety only. They do no targeting and attribute no intent.

kessler Commons grows out of [kessler](https://github.com/Shakhtar-Sankur/kessler), which provides
the GPU SGP4, the all-vs-all conjunction screening and the single-operator maneuver detector.

## Status

| | What | State |
|---|---|---|
| M1 | Multi-operator maneuver detector over the 6-hourly operator-ephemeris archive | **in progress**: detector ported and tested; runs over the archive privately |
| M1 | Storm-drag separation (NOAA SWPC Kp/F10.7; propulsionless satellites and debris as controls) | next |
| M2 | Validation against operator ephemerides and published maneuver logs | planned |
| M3 | Versioned open catalogue with a DOI | planned, after data rights are cleared |
| M4 | Comparison with operators' published counts; per-storm altitude loss | planned |
| M5 | Public site: operator timelines, storm view | planned |

## Data and publication

Inputs are public: CelesTrak's supplemental GP element sets, which are derived from operators'
own ephemerides, along with CelesTrak's satellite catalogue. TraCSS open data (CC0) will be added
when its public service opens. No Space-Track data is used.

**Results are not published yet.** `out/` is git-ignored. Detections, counts and per-operator
statistics will appear here only after the data providers confirm that publishing derived results
is allowed. Until then, this repository holds the methods and the tests that check them.

## How it works

- **`commons/orbits.py`**
  - Recovers Brouwer's mean semi-major axis from each element set, exactly as SGP4's
    initialisation does.
  - Finds steps in that axis. At every boundary between element sets, it fits a level, a drag
    rate and a step over the neighbouring sets.
  - Counts a step as a maneuver only when two conditions hold:
    - it is above a size floor (50 m);
    - it is more than 6 standard errors of the satellite's own noise.
  - Measures that noise as the median absolute deviation of the satellite's drag-corrected jumps.
- **`commons/watch.py`**
  - Joins each satellite's element sets across all archive snapshots and drops repeats.
  - Runs the detector.
  - Writes `out/maneuvers.csv` (one row per detected change) and `out/summary.json` (per
    constellation: satellites, satellite-days observed, changes, raises and lowerings, median
    size).

```
pip install -r requirements.txt
python -m pytest -q
python -m commons.watch --archive path/to/archive/data --out out
```

The tests build a synthetic archive with known maneuvers, drag and noise, and check two things:
- every injected maneuver is found with its size;
- nothing is found on satellites that only feel drag.

## Wording

Results are reported as "detected from public data" or "could not be verified from public data",
never as claims about an operator's conduct. Every number will come with its detection limits and
the script that produced it.

## License

Apache-2.0.
