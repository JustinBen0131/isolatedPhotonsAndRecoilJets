#!/usr/bin/env python3
"""Build split-study compact JSON from validation score_histograms reports."""

from __future__ import annotations
# Keep purpose-folder helpers runnable when invoked directly.
import sys as _codex_sys
from pathlib import Path as _CodexPath
_CODEX_THIS_FILE = _CodexPath(__file__).resolve()
_CODEX_SCRIPTS_DIR = next((p for p in _CODEX_THIS_FILE.parents if p.name == "scripts"), _CODEX_THIS_FILE.parent)
_CODEX_SCRIPTS_DIR_STR = str(_CODEX_SCRIPTS_DIR)
if _CODEX_SCRIPTS_DIR_STR not in _codex_sys.path:
    _codex_sys.path.append(_CODEX_SCRIPTS_DIR_STR)
del _CODEX_THIS_FILE, _CODEX_SCRIPTS_DIR, _CODEX_SCRIPTS_DIR_STR

import argparse
import json
from pathlib import Path

import numpy as np


SPLITS = [
    ("90/10", "train90_test10", 0.90),
    ("50/50", "train50_test50", 0.50),
    ("10/90", "train10_test90", 0.10),
]
CENTS = [("0-20%", "0_20", 0, 20), ("20-50%", "20_50", 20, 50), ("50-80%", "50_80", 50, 80)]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reports", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--report-prefix", default="model_validation_condor_THE8_jet40_centinput_overfitdiag_20260528_")
    ap.add_argument("--report-suffix", default="_quickstat_20260528_1450")
    ap.add_argument("--product", default="centInput_pt1535")
    return ap.parse_args()


def read_summary(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            out[key] = value
    return out


def auc_from_counts(sig_counts: list[int], bkg_counts: list[int]) -> float:
    sig = np.asarray(sig_counts, dtype=float)
    bkg = np.asarray(bkg_counts, dtype=float)
    if sig.sum() <= 0 or bkg.sum() <= 0:
        return float("nan")
    sig = sig / sig.sum()
    bkg = bkg / bkg.sum()
    bkg_less = np.cumsum(bkg) - bkg
    return float(np.sum(sig * (bkg_less + 0.5 * bkg)))


def main() -> None:
    args = parse_args()
    histograms = []
    rows = []
    for split_label, tag, train_frac in SPLITS:
        report = args.reports / f"{args.report_prefix}{tag}{args.report_suffix}"
        score_payload = json.loads((report / "score_histograms.json").read_text())
        summary = read_summary(report / "validation_summary.txt")
        product = score_payload["products"][args.product]
        rows.append(
            {
                "label": split_label,
                "train_fraction": train_frac,
                "test_fraction": 1.0 - train_frac,
                "global_auc": float(summary[f"{args.product}_auc"]),
                "total_entries": int(summary["total_entries"]),
                "signal_entries": int(summary["signal_entries"]),
                "background_entries": int(summary["background_entries"]),
            }
        )
        for cent_label, cent_key, cent_lo, cent_hi in CENTS:
            cent = product["by_centrality"][cent_key]
            sig = cent["signal"]
            bkg = cent["background"]
            histograms.append(
                {
                    "split_label": split_label,
                    "cent_label": cent_label,
                    "cent_lo": cent_lo,
                    "cent_hi": cent_hi,
                    "bin_edges": score_payload["bin_edges"],
                    "signal_density": sig["density"],
                    "background_density": bkg["density"],
                    "signal_entries": int(sig["entries"]),
                    "background_entries": int(bkg["entries"]),
                    "entries": int(sig["entries"]) + int(bkg["entries"]),
                    "auc": auc_from_counts(sig["counts"], bkg["counts"]),
                }
            )
    payload = {
        "schema": "auau_bdt_splitstudy_compact_v1",
        "source": "THE8 Jet40-inclusive quickstat validation reports",
        "score_key": f"score_{args.product}",
        "selection": "15 < cluster_Et < 35 GeV, 0 <= centrality < 80",
        "display_split": "common full-stat score-cache validation rows from Jet12+20+30+40 embedded-inclusive background",
        "random_seed": None,
        "train_fraction_for_display_split": None,
        "val_fraction_for_display_split": None,
        "centrality_bins": [{"label": label, "lo": lo, "hi": hi} for label, _key, lo, hi in CENTS],
        "rows": rows,
        "histograms": histograms,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(args.out)


if __name__ == "__main__":
    main()
