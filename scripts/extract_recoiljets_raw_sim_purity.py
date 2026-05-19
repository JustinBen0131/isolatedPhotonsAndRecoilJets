#!/usr/bin/env python3
"""Extract raw SIM purity points from paired RecoilJets signal/background ROOTs."""

from __future__ import annotations

import argparse
import csv
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


def open_sim(path: Path):
    handle = ROOT.TFile.Open(str(path), "READ")
    if not handle or handle.IsZombie():
        raise RuntimeError(f"Could not open {path}")
    directory = handle.Get("SIM")
    if not directory:
        raise RuntimeError(f"Missing SIM directory in {path}")
    return handle, directory


def hist_integral_and_error(hist) -> tuple[float, float]:
    total = 0.0
    err2 = 0.0
    for ibin in range(0, hist.GetNbinsX() + 2):
        total += float(hist.GetBinContent(ibin))
        err2 += float(hist.GetBinError(ibin)) ** 2
    return total, math.sqrt(err2)


def index_histograms(directory, prefixes: tuple[str, ...], cone: str) -> dict[tuple[str, str, float, float], str]:
    patterns = {
        prefix: re.compile(rf"^{re.escape(prefix)}_{re.escape(cone)}_pT_([0-9]+)_([0-9]+)_cent_(.+)$")
        for prefix in prefixes
    }
    index: dict[tuple[str, str, float, float], str] = {}
    keys = directory.GetListOfKeys()
    for ikey in range(keys.GetEntries()):
        name = keys.At(ikey).GetName()
        for prefix, pattern in patterns.items():
            match = pattern.match(name)
            if not match:
                continue
            lo = float(match.group(1))
            hi = float(match.group(2))
            cent = match.group(3)
            index[(prefix, cent, lo, hi)] = name
    return index


def sum_count_for_bin(directory, index: dict[tuple[str, str, float, float], str], prefix: str, cent: str, lo: float, hi: float) -> tuple[float, float]:
    pieces = (cent,)
    if (prefix, cent, lo, hi) not in index:
        pieces = COARSE_CENTS[cent]

    total = 0.0
    err2 = 0.0
    missing: list[str] = []
    for piece in pieces:
        name = index.get((prefix, piece, lo, hi))
        if not name:
            missing.append(f"{prefix}_pT_{lo:g}_{hi:g}_cent_{piece}")
            continue
        hist = directory.Get(name)
        if not hist:
            missing.append(name)
            continue
        value, error = hist_integral_and_error(hist)
        total += value
        err2 += error * error
    if missing and len(missing) == len(pieces):
        raise RuntimeError(f"Missing all histograms for {prefix}, pT {lo:g}-{hi:g}, cent {cent}: {missing}")
    return total, math.sqrt(err2)


def bins_for_cent(index: dict[tuple[str, str, float, float], str], prefix: str, cent: str, pt_min: float, pt_max: float) -> list[tuple[float, float]]:
    direct = sorted(
        (lo, hi)
        for key_prefix, key_cent, lo, hi in index
        if key_prefix == prefix and key_cent == cent and pt_min <= 0.5 * (lo + hi) <= pt_max
    )
    if direct:
        return direct
    bins: set[tuple[float, float]] = set()
    for fine_cent in COARSE_CENTS[cent]:
        bins.update(
            (lo, hi)
            for key_prefix, key_cent, lo, hi in index
            if key_prefix == prefix and key_cent == fine_cent and pt_min <= 0.5 * (lo + hi) <= pt_max
        )
    return sorted(bins)


def purity_points(signal_root: Path, background_root: Path, pt_min: float, pt_max: float, cone: str) -> list[dict]:
    sig_handle, sig_dir = open_sim(signal_root)
    bkg_handle, bkg_dir = open_sim(background_root)
    signal_prefix = "h_EisoReco_truthSigMatched_tight"
    background_prefix = "h_Eiso_tight"
    sig_index = index_histograms(sig_dir, (signal_prefix,), cone)
    bkg_index = index_histograms(bkg_dir, (background_prefix,), cone)

    rows: list[dict] = []
    try:
        for cent in COARSE_CENTS:
            sig_bins = set(bins_for_cent(sig_index, signal_prefix, cent, pt_min, pt_max))
            bkg_bins = set(bins_for_cent(bkg_index, background_prefix, cent, pt_min, pt_max))
            for lo, hi in sorted(sig_bins & bkg_bins):
                sig, sig_err = sum_count_for_bin(sig_dir, sig_index, signal_prefix, cent, lo, hi)
                bkg, bkg_err = sum_count_for_bin(bkg_dir, bkg_index, background_prefix, cent, lo, hi)
                denom = sig + bkg
                if denom <= 0:
                    continue
                purity = sig / denom
                err = math.sqrt((bkg * bkg * sig_err * sig_err + sig * sig * bkg_err * bkg_err) / (denom**4))
                rows.append(
                    {
                        "source": "recoiljets_final_root",
                        "centrality": cent,
                        "pt_low": lo,
                        "pt_high": hi,
                        "pt_mid": 0.5 * (lo + hi),
                        "value": purity,
                        "error": err,
                        "signal_tight": sig,
                        "signal_tight_error": sig_err,
                        "background_tight": bkg,
                        "background_tight_error": bkg_err,
                        "signal_root": str(signal_root),
                        "background_root": str(background_root),
                    }
                )
    finally:
        sig_handle.Close()
        bkg_handle.Close()
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--signal-root", required=True, type=Path)
    parser.add_argument("--background-root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--pt-min", type=float, default=15.0)
    parser.add_argument("--pt-max", type=float, default=35.0)
    parser.add_argument("--cone", default="isoR30")
    args = parser.parse_args()

    rows = purity_points(args.signal_root, args.background_root, args.pt_min, args.pt_max, args.cone)
    if not rows:
        raise RuntimeError("No purity rows were produced")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
