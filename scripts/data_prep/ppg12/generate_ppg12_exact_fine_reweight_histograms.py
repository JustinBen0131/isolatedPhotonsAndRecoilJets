#!/usr/bin/env python3
"""Generate fine ET/eta PPG12-exact reweighting closure histograms from an NPZ cache."""

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
import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--score-cache", required=True, type=Path)
    parser.add_argument("--output", default="-", help="CSV output path, or '-' for stdout")
    parser.add_argument("--et-min", type=float, default=5.0)
    parser.add_argument("--et-max", type=float, default=35.0)
    parser.add_argument("--et-bins", type=int, default=20)
    return parser.parse_args()


def density(values: np.ndarray, bins: np.ndarray, weights: np.ndarray | None = None) -> np.ndarray:
    hist, _ = np.histogram(values, bins=bins, weights=weights)
    total = float(np.sum(hist))
    widths = np.diff(bins)
    if total <= 0.0:
        return np.zeros(len(widths), dtype="float64")
    return hist.astype("float64") / total / widths


def main() -> None:
    args = parse_args()
    repo = Path.cwd()
    sys.path.insert(0, str(repo / "scripts"))
    from train_auau_photon_bdt import compute_ppg12_exact_global_weights

    cache = np.load(args.score_cache, allow_pickle=True)
    frame = pd.DataFrame(
        {
            "is_signal": cache["is_signal"].astype("int32"),
            "cluster_Et": cache["cluster_Et"].astype("float64"),
            "cluster_Eta": cache["cluster_Eta"].astype("float64"),
        }
    )
    weights, _ = compute_ppg12_exact_global_weights(frame, "is_signal")
    labels = frame["is_signal"].to_numpy(dtype="int32")

    specs = [
        ("cluster_Et", np.linspace(args.et_min, args.et_max, args.et_bins + 1)),
        ("cluster_Eta", np.linspace(-0.7, 0.7, 21)),
    ]

    rows: list[dict[str, object]] = []
    for axis, bins in specs:
        values = frame[axis].to_numpy(dtype="float64")
        finite_axis = np.isfinite(values)
        for cls in (0, 1):
            mask = finite_axis & (labels == cls)
            raw_counts, _ = np.histogram(values[mask], bins=bins)
            weighted_sums, _ = np.histogram(values[mask], bins=bins, weights=weights[mask])
            raw_density = density(values[mask], bins)
            weighted_density = density(values[mask], bins, weights=weights[mask])
            for i, (lo, hi) in enumerate(zip(bins[:-1], bins[1:])):
                rows.append(
                    {
                        "axis": axis,
                        "class": cls,
                        "bin_low": float(lo),
                        "bin_high": float(hi),
                        "raw_count": int(raw_counts[i]),
                        "weighted_sum": float(weighted_sums[i]),
                        "raw_density": float(raw_density[i]),
                        "weighted_density": float(weighted_density[i]),
                    }
                )

    fieldnames = [
        "axis",
        "class",
        "bin_low",
        "bin_high",
        "raw_count",
        "weighted_sum",
        "raw_density",
        "weighted_density",
    ]
    if args.output == "-":
        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    else:
        out = Path(args.output)
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        print(out)


if __name__ == "__main__":
    main()
