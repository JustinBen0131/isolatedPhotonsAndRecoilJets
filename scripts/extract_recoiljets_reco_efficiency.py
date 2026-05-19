#!/usr/bin/env python3
"""Extract truth-photon reco efficiency points from a RecoilJets SIM ROOT file."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import ROOT


ROOT.gROOT.SetBatch(True)

COARSE_CENTS = {
    "0_20": ("0_10", "10_20"),
    "20_50": ("20_30", "30_40", "40_50"),
    "50_80": ("50_60", "60_80"),
}


def get_hist(directory, name: str):
    hist = directory.Get(name)
    if not hist:
        raise RuntimeError(f"Missing histogram {name}")
    return hist


def summed_hist(directory, pattern: str, cent: str):
    direct = directory.Get(pattern.format(cent=cent))
    if direct:
        return direct.Clone(f"{direct.GetName()}_clone_for_{cent}")

    pieces = []
    for fine in COARSE_CENTS[cent]:
        hist = directory.Get(pattern.format(cent=fine))
        if not hist:
            raise RuntimeError(f"Missing histogram {pattern.format(cent=fine)}")
        pieces.append(hist)
    out = pieces[0].Clone(f"{pieces[0].GetName()}_summed_for_{cent}")
    out.SetDirectory(0)
    for hist in pieces[1:]:
        out.Add(hist)
    return out


def points(root_path: Path, pt_min: float, pt_max: float) -> list[dict]:
    handle = ROOT.TFile.Open(str(root_path), "READ")
    if not handle or handle.IsZombie():
        raise RuntimeError(f"Could not open {root_path}")
    directory = handle.Get("SIM")
    if not directory:
        raise RuntimeError(f"Missing SIM directory in {root_path}")

    rows: list[dict] = []
    for cent in COARSE_CENTS:
        h_truth = summed_hist(directory, "h_unfoldTruthPho_pTgamma_cent_{cent}", cent)
        h_miss = summed_hist(directory, "h_unfoldTruthPhoMisses_pTgamma_isoR30_isSliding_cent_{cent}", cent)
        for ib in range(1, h_truth.GetNbinsX() + 1):
            lo = h_truth.GetXaxis().GetBinLowEdge(ib)
            hi = h_truth.GetXaxis().GetBinUpEdge(ib)
            mid = 0.5 * (lo + hi)
            if mid < pt_min or mid > pt_max:
                continue
            miss_bin = h_miss.GetXaxis().FindBin(mid)
            if miss_bin < 1 or miss_bin > h_miss.GetNbinsX():
                continue
            truth = float(h_truth.GetBinContent(ib))
            miss = float(h_miss.GetBinContent(miss_bin))
            if truth <= 0:
                continue
            reco = truth - miss
            eff = reco / truth
            err = 0.0
            if reco >= 0 and miss >= 0:
                e_truth = float(h_truth.GetBinError(ib))
                e_miss = float(h_miss.GetBinError(miss_bin))
                var_miss = e_miss * e_miss
                # The truth histogram is the sum of reconstructed + missed
                # truth photons, so infer the reconstructed variance when the
                # two categories are disjoint. Fall back gracefully for any
                # weighted-bin roundoff where var_truth < var_miss.
                var_reco = max(e_truth * e_truth - var_miss, 0.0)
                denom2 = truth * truth
                err = math.sqrt((miss * miss * var_reco + reco * reco * var_miss) / (denom2 * denom2))
            rows.append(
                {
                    "source": "recoiljets_final_root",
                    "centrality": cent,
                    "pt_low": lo,
                    "pt_high": hi,
                    "pt_mid": mid,
                    "value": eff,
                    "error": err,
                    "numerator": reco,
                    "denominator": truth,
                    "missed": miss,
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
    args = parser.parse_args()

    rows = points(args.root, args.pt_min, args.pt_max)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
