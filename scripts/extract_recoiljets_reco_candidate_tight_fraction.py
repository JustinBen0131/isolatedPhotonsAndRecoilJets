#!/usr/bin/env python3
"""Extract truth-matched reco-candidate tight-ID fraction vs reco photon pT."""

from __future__ import annotations

import argparse
import csv
import ctypes
import math
import re
from pathlib import Path

import ROOT


ROOT.gROOT.SetBatch(True)

COARSE_CENTS = {
    "0_20": ("0_10", "10_20"),
    "20_50": ("20_30", "30_40", "40_50"),
    "50_80": ("50_60", "60_80"),
}


def integral_and_error(hist) -> tuple[float, float]:
    err = ctypes.c_double(0.0)
    val = float(hist.IntegralAndError(1, hist.GetNbinsX(), err))
    return val, float(err.value)


def tight_fraction(tight: float, all_reco: float, e_tight: float, e_all: float) -> tuple[float, float, float]:
    if all_reco <= 0:
        return math.nan, math.nan, 0.0
    non_tight = max(all_reco - tight, 0.0)
    var_tight = e_tight * e_tight
    # The all-reco histogram is tight + non-tight. Infer the disjoint
    # non-tight variance, with a floor for weighted-bin roundoff.
    var_non_tight = max(e_all * e_all - var_tight, 0.0)
    den = tight + non_tight
    value = tight / den
    err = math.sqrt((non_tight * non_tight * var_tight + tight * tight * var_non_tight) / (den * den * den * den))
    return value, err, non_tight


def points(root_path: Path, pt_min: float, pt_max: float, cone: str) -> list[dict]:
    handle = ROOT.TFile.Open(str(root_path), "READ")
    if not handle or handle.IsZombie():
        raise RuntimeError(f"Could not open {root_path}")
    directory = handle.Get("SIM")
    if not directory:
        raise RuntimeError(f"Missing SIM directory in {root_path}")

    all_pattern = re.compile(rf"^h_EisoReco_truthSigMatched_{cone}_pT_([0-9]+)_([0-9]+)_cent_(.+)$")
    keys = directory.GetListOfKeys()
    by_bin: dict[tuple[float, float, str], str] = {}
    for ik in range(keys.GetEntries()):
        name = keys.At(ik).GetName()
        match = all_pattern.match(name)
        if not match:
            continue
        lo = float(match.group(1))
        hi = float(match.group(2))
        fine_cent = match.group(3)
        mid = 0.5 * (lo + hi)
        if mid < pt_min or mid > pt_max:
            continue
        by_bin[(lo, hi, fine_cent)] = name

    rows: list[dict] = []
    for cent, fine_cents in COARSE_CENTS.items():
        available_cents = set(fine_cents)
        if any(fine == cent for _, _, fine in by_bin):
            available_cents.add(cent)
        pt_bins = sorted({(lo, hi) for lo, hi, fine in by_bin if fine in available_cents})
        for lo, hi in pt_bins:
            tight_sum = all_sum = 0.0
            var_tight = var_all = 0.0
            for fine_cent in sorted(available_cents):
                all_name = by_bin.get((lo, hi, fine_cent))
                if not all_name:
                    continue
                tight_name = all_name.replace("h_EisoReco_truthSigMatched_", "h_EisoReco_truthSigMatched_tight_", 1)
                h_all = directory.Get(all_name)
                h_tight = directory.Get(tight_name)
                if not h_all or not h_tight:
                    raise RuntimeError(f"Missing paired histograms for {all_name} / {tight_name}")
                all_val, all_err = integral_and_error(h_all)
                tight_val, tight_err = integral_and_error(h_tight)
                all_sum += all_val
                tight_sum += tight_val
                var_all += all_err * all_err
                var_tight += tight_err * tight_err
            value, err, non_tight = tight_fraction(tight_sum, all_sum, math.sqrt(var_tight), math.sqrt(var_all))
            if not math.isfinite(value):
                continue
            rows.append(
                {
                    "source": "recoiljets_final_root",
                    "centrality": cent,
                    "pt_low": lo,
                    "pt_high": hi,
                    "pt_mid": 0.5 * (lo + hi),
                    "value": value,
                    "error": err,
                    "numerator": tight_sum,
                    "denominator": all_sum,
                    "non_tight": non_tight,
                    "cone": cone,
                    "root_path": str(root_path),
                }
            )
    handle.Close()
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--pt-min", type=float, default=15.0)
    parser.add_argument("--pt-max", type=float, default=35.0)
    parser.add_argument("--cone", default="isoR30")
    args = parser.parse_args()

    rows = points(args.root, args.pt_min, args.pt_max, args.cone)
    if not rows:
        raise RuntimeError("No reco-candidate tight-fraction points were found")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
