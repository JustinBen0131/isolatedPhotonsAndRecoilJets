#!/usr/bin/env python3
"""Summarize THE-76 local PPG12 A/B/C/D histograms against tree tags."""

from __future__ import annotations

import argparse
from pathlib import Path

import ROOT  # type: ignore


ROOT.gROOT.SetBatch(True)

RECO_BINS = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36]
REGION_HISTS = {
    "A": "SIM/h_tight_iso_cluster_signal_0",
    "B": "SIM/h_tight_noniso_cluster_signal_0",
    "C": "SIM/h_nontight_iso_cluster_signal_0",
    "D": "SIM/h_nontight_noniso_cluster_signal_0",
}


def find_bin(x: float) -> int:
    for i in range(len(RECO_BINS) - 1):
        if RECO_BINS[i] < x < RECO_BINS[i + 1]:
            return i
    return -1


def bin_label(i: int) -> str:
    return f"{RECO_BINS[i]}-{RECO_BINS[i + 1]}"


def safe_div(n: float, d: float) -> float:
    return n / d if d else float("nan")


def hist_region_counts(f: ROOT.TFile) -> dict[str, list[float]]:
    out: dict[str, list[float]] = {}
    for region, name in REGION_HISTS.items():
        h = f.Get(name)
        if not h:
            raise RuntimeError(f"Missing histogram {name}")
        out[region] = [float(h.GetBinContent(i + 1)) for i in range(len(RECO_BINS) - 1)]
    return out


def tree_region_counts(f: ROOT.TFile, weighted: bool) -> dict[str, list[float]]:
    t = f.Get("AuAuPhotonIDTrainingTree")
    if not t:
        raise RuntimeError("Missing AuAuPhotonIDTrainingTree")
    out = {region: [0.0 for _ in range(len(RECO_BINS) - 1)] for region in "ABCD"}
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        if int(getattr(t, "is_signal")) != 1:
            continue
        ib = find_bin(float(getattr(t, "cluster_Et")))
        if ib < 0:
            continue
        tag = int(getattr(t, "ppg12_tight_tag"))
        is_iso = int(getattr(t, "ppg12_is_iso")) == 1
        is_noniso = int(getattr(t, "ppg12_is_noniso")) == 1
        in_window = int(getattr(t, "ppg12_truth_window_pass_r04")) == 1
        weight = float(getattr(t, "event_weight")) if weighted else 1.0
        if tag == 1 and is_iso and in_window:
            out["A"][ib] += weight
        if tag == 1 and is_noniso:
            out["B"][ib] += weight
        if tag == 2 and is_iso:
            out["C"][ib] += weight
        if tag == 2 and is_noniso:
            out["D"][ib] += weight
    return out


def print_summary(label: str, counts: dict[str, list[float]]) -> None:
    print(f"## {label}")
    print("bin,A,B,C,D,B/A,C/A,D/A,CD/A,D_over_CD")
    for i in range(len(RECO_BINS) - 1):
        a = counts["A"][i]
        b = counts["B"][i]
        c = counts["C"][i]
        d = counts["D"][i]
        print(
            f"{bin_label(i)},{a:.9g},{b:.9g},{c:.9g},{d:.9g},"
            f"{safe_div(b, a):.9g},{safe_div(c, a):.9g},"
            f"{safe_div(d, a):.9g},{safe_div(c + d, a):.9g},"
            f"{safe_div(d, c + d):.9g}"
        )
    print()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root_file")
    args = parser.parse_args()
    path = Path(args.root_file)
    f = ROOT.TFile.Open(str(path))
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open {path}")
    print(f"# {path}")
    print_summary("hist", hist_region_counts(f))
    print_summary("tree_unweighted", tree_region_counts(f, weighted=False))
    print_summary("tree_weighted", tree_region_counts(f, weighted=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
