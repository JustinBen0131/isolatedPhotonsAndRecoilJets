#!/usr/bin/env python3
"""Calculate photon-trigger sampled luminosity with the same pp/AuAu formula.

For each run: L = MB_live * (photon_scaled / photon_live)
                 * pileup * mb_coverage / sigma_MB.
See PHOTON_LUMINOSITY.md for counter definitions and the PPG12 mapping.
Only Python's standard library is required; no event loop or batch jobs.
"""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def photon_luminosity(mb_live, photon_live, photon_scaled, sigma_barn,
                      *, pileup=1.0, mb_coverage=1.0):
    """Return one run's sampled luminosity in inverse barns.

    live counters include DAQ live-time gating and precede prescaling.
    scaled counters follow prescaling. Their photon ratio is the recorded
    fraction; do not apply another photon prescale correction afterward.

    pileup is an independently supplied MB pileup correction. mb_coverage
    converts the scaler MB count convention to the chosen observed-event
    convention (PPG12: observed MB / scaled MB). It is not a luminosity fit
    or necessarily a fraction <= 1. Defaults of one mean no correction.
    Counters must describe matching acquisition intervals and trigger menus.
    """
    values = (mb_live, photon_live, photon_scaled, sigma_barn, pileup, mb_coverage)
    if not all(math.isfinite(v) for v in values):
        raise ValueError("nonfinite luminosity input")
    if min(mb_live, photon_live, photon_scaled, mb_coverage) < 0 or sigma_barn <= 0 or pileup <= 0:
        raise ValueError("invalid luminosity input")
    if photon_scaled > photon_live:
        raise ValueError("scaled photon count exceeds live count")
    if photon_scaled == 0:
        return 0.0
    if photon_live <= 0:
        raise ValueError("positive scaled count without live count")
    return mb_live * (photon_scaled / photon_live) * pileup * mb_coverage / sigma_barn


def calculate_csv(path, sigma_barn, unit):
    """Read explicitly selected runs; reject missing fields and duplicate runs.

    CSV columns: run,mb_live,photon_live,photon_scaled,pileup,mb_coverage.
    Both factors are mandatory in the CSV so unity assumptions are visible.
    Missing counters must be resolved by the input producer, not silently
    converted to zero. Include known zero-exposure runs with explicit zeros.
    """
    divisor = {"nb^-1": 1e9, "pb^-1": 1e12}[unit]
    if not math.isfinite(sigma_barn) or sigma_barn <= 0:
        raise ValueError("sigma_barn must be finite and positive")
    rows, seen = [], set()
    with Path(path).open(newline="") as stream:
        for row in csv.DictReader(stream):
            run = int(row["run"])
            if run <= 0 or run in seen:
                raise ValueError(f"invalid or duplicate run: {run}")
            seen.add(run)
            inputs = {key: float(row[key]) for key in
                      ("mb_live", "photon_live", "photon_scaled", "pileup", "mb_coverage")}
            value = photon_luminosity(**inputs, sigma_barn=sigma_barn) / divisor
            rows.append({"run": run, **inputs, "luminosity": value})
    if not rows:
        raise ValueError("input contains no runs")
    return {"unit": unit, "sigma_barn": sigma_barn, "runs": len(rows),
            "positive_runs": sum(r["luminosity"] > 0 for r in rows),
            "luminosity": math.fsum(r["luminosity"] for r in rows),
            "input_sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
            "rows": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Per-run counter CSV")
    parser.add_argument("--sigma-barn", type=float, required=True)
    parser.add_argument("--unit", choices=("nb^-1", "pb^-1"), required=True)
    args = parser.parse_args()
    print(json.dumps(calculate_csv(args.input, args.sigma_barn, args.unit), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
