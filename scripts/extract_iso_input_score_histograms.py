#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


CENT_BINS = [(0.0, 20.0, "0-20%"), (20.0, 50.0, "20-50%"), (50.0, 80.0, "50-80%")]
MODELS = [
    ("Baseline inputs", "baseline_inputs", "score_globalEtCent1535_bdt_noIso_ptCent7"),
    ("Baseline + isolation inputs", "baseline_plus_isolation_inputs", "score_globalEtCent1535_bdt_iso_ptCent7"),
]


def score_cache_paths(path: Path) -> list[Path]:
    if path.is_dir():
        manifest = path / "score_caches.local.list"
        if not manifest.is_file():
            manifest = path / "score_caches.list"
    else:
        manifest = path
    if manifest.suffix == ".npz" and manifest.is_file():
        return [manifest]
    paths = [Path(line.strip()) for line in manifest.read_text().splitlines() if line.strip()]
    if not paths:
        raise SystemExit(f"No score-cache paths in {manifest}")
    return paths


def auc_from_scores(y: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y, dtype=np.int8)
    score = np.asarray(score, dtype=np.float64)
    mask = np.isin(y, [0, 1]) & np.isfinite(score)
    y = y[mask]
    score = score[mask]
    n_sig = int(np.sum(y == 1))
    n_bkg = int(np.sum(y == 0))
    if n_sig == 0 or n_bkg == 0:
        return float("nan")
    order = np.argsort(score, kind="mergesort")
    sorted_score = score[order]
    ranks = np.empty(score.size, dtype=np.float64)
    i = 0
    while i < score.size:
        j = i + 1
        while j < score.size and sorted_score[j] == sorted_score[i]:
            j += 1
        ranks[order[i:j]] = 0.5 * (i + 1 + j)
        i = j
    rank_sum_sig = float(np.sum(ranks[y == 1]))
    return (rank_sum_sig - n_sig * (n_sig + 1) / 2.0) / (n_sig * n_bkg)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    paths = score_cache_paths(args.manifest)
    bins = np.linspace(0.0, 1.0, 51)
    rows: list[dict[str, object]] = []

    for model_label, model_slug, score_key in MODELS:
        for clo, chi, cent_label in CENT_BINS:
            sig_parts: list[np.ndarray] = []
            bkg_parts: list[np.ndarray] = []
            for path in paths:
                with np.load(path, allow_pickle=True) as z:
                    y = np.asarray(z["is_signal"], dtype=np.int8)
                    et = np.asarray(z["cluster_Et"], dtype=np.float32)
                    cent = np.asarray(z["centrality"], dtype=np.float32)
                    score = np.asarray(z[score_key], dtype=np.float32)
                    mask = (
                        np.isin(y, [0, 1])
                        & np.isfinite(et)
                        & np.isfinite(cent)
                        & np.isfinite(score)
                        & (et >= 15.0)
                        & (et < 35.0)
                        & (cent >= clo)
                        & (cent < chi)
                    )
                    if np.any(mask & (y == 1)):
                        sig_parts.append(score[mask & (y == 1)])
                    if np.any(mask & (y == 0)):
                        bkg_parts.append(score[mask & (y == 0)])
            sig = np.concatenate(sig_parts) if sig_parts else np.array([], dtype=np.float32)
            bkg = np.concatenate(bkg_parts) if bkg_parts else np.array([], dtype=np.float32)
            sig_density = np.histogram(sig, bins=bins, density=True)[0] if sig.size else np.zeros(len(bins) - 1)
            bkg_density = np.histogram(bkg, bins=bins, density=True)[0] if bkg.size else np.zeros(len(bins) - 1)
            auc = auc_from_scores(np.r_[np.ones(sig.size, dtype=np.int8), np.zeros(bkg.size, dtype=np.int8)], np.r_[sig, bkg])
            for lo, hi, sd, bd in zip(bins[:-1], bins[1:], sig_density, bkg_density):
                rows.append(
                    {
                        "model_row": model_label,
                        "model_slug": model_slug,
                        "score_key": score_key,
                        "centrality_bin": cent_label,
                        "bin_lo": float(lo),
                        "bin_hi": float(hi),
                        "signal_density": float(sd),
                        "background_density": float(bd),
                        "signal_entries": int(sig.size),
                        "background_entries": int(bkg.size),
                        "auc": float(auc),
                    }
                )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
