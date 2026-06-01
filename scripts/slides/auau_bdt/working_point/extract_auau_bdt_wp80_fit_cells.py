#!/usr/bin/env python3
"""Extract Au+Au BDT WP80 threshold cells and linear ET fits."""

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
import json
import math
from pathlib import Path

import numpy as np


DEFAULT_PT_EDGES = "15,17,19,21,23,25,27,30,35"
DEFAULT_CENT_EDGES = "0,20,50,80"


def parse_edges(spec: str) -> list[float]:
    values = [float(tok) for tok in spec.split(",") if tok.strip()]
    if len(values) < 2 or any(values[i + 1] <= values[i] for i in range(len(values) - 1)):
        raise ValueError(f"Bad edge list: {spec}")
    return values


def score_cache_paths(path: Path) -> list[Path]:
    if path.is_dir():
        manifest = path / "score_caches.list"
    else:
        manifest = path
    if manifest.suffix == ".npz":
        return [manifest]
    paths = [Path(line.strip()) for line in manifest.read_text().splitlines() if line.strip()]
    if not paths:
        raise FileNotFoundError(f"No score caches listed in {manifest}")
    return paths


def threshold_for_target(y: np.ndarray, score: np.ndarray, target: float) -> dict | None:
    valid = np.isfinite(score) & np.isin(y, [0, 1])
    sig = score[valid & (y == 1)]
    bkg = score[valid & (y == 0)]
    if sig.size == 0 or bkg.size == 0:
        return None
    threshold = float(np.quantile(sig, max(0.0, min(1.0, 1.0 - target))))
    return {
        "threshold": threshold,
        "signal_efficiency": float(np.mean(sig > threshold)),
        "background_fake_rate": float(np.mean(bkg > threshold)),
        "signal_entries": int(sig.size),
        "background_entries": int(bkg.size),
    }


def fit_rows(rows: list[dict]) -> dict:
    good = [r for r in rows if r["source"] == "cell" and math.isfinite(float(r["threshold"]))]
    if len(good) < 2:
        return {"slope": math.nan, "intercept": math.nan, "max_abs_residual": math.nan, "n_fit_points": len(good)}
    x = np.array([r["pt_center"] for r in good], dtype=float)
    t = np.array([r["threshold"] for r in good], dtype=float)
    w = np.sqrt(np.array([max(1, r["signal_entries"]) for r in good], dtype=float))
    slope, intercept = np.polyfit(x, t, 1, w=w)
    pred = slope * x + intercept
    return {
        "slope": float(slope),
        "intercept": float(intercept),
        "max_abs_residual": float(np.max(np.abs(t - pred))),
        "n_fit_points": int(len(good)),
    }


def compute(args: argparse.Namespace) -> dict:
    pt_edges = parse_edges(args.pt_edges)
    cent_edges = parse_edges(args.cent_edges)
    score_key = args.score_key
    rows: list[dict] = []
    inclusive_sig: list[np.ndarray] = []
    inclusive_bkg: list[np.ndarray] = []
    by_cell = {
        (clo, chi, plo, phi, cls): []
        for clo, chi in zip(cent_edges[:-1], cent_edges[1:])
        for plo, phi in zip(pt_edges[:-1], pt_edges[1:])
        for cls in (0, 1)
    }
    files_loaded = 0
    rows_loaded = 0

    for cache in score_cache_paths(args.manifest):
        with np.load(cache, allow_pickle=True) as z:
            missing = [key for key in ("is_signal", "cluster_Et", "centrality", score_key) if key not in z.files]
            if missing:
                raise KeyError(f"{cache} missing keys {missing}")
            y = np.asarray(z["is_signal"], dtype=np.int8)
            et = np.asarray(z["cluster_Et"], dtype=np.float32)
            cent = np.asarray(z["centrality"], dtype=np.float32)
            score = np.asarray(z[score_key], dtype=np.float32)
            base = (
                np.isfinite(et)
                & np.isfinite(cent)
                & np.isfinite(score)
                & np.isin(y, [0, 1])
                & (et >= pt_edges[0])
                & (et < pt_edges[-1])
                & (cent >= cent_edges[0])
                & (cent < cent_edges[-1])
            )
            files_loaded += 1
            rows_loaded += int(np.sum(base))
            if not np.any(base):
                continue
            inclusive_sig.append(score[base & (y == 1)])
            inclusive_bkg.append(score[base & (y == 0)])
            for clo, chi in zip(cent_edges[:-1], cent_edges[1:]):
                cmask = base & (cent >= clo) & (cent < chi)
                for plo, phi in zip(pt_edges[:-1], pt_edges[1:]):
                    mask = cmask & (et >= plo) & (et < phi)
                    if not np.any(mask):
                        continue
                    for cls in (0, 1):
                        arr = score[mask & (y == cls)]
                        if arr.size:
                            by_cell[(clo, chi, plo, phi, cls)].append(arr.astype(np.float32, copy=False))

    sig_all = np.concatenate(inclusive_sig) if inclusive_sig else np.array([], dtype=np.float32)
    bkg_all = np.concatenate(inclusive_bkg) if inclusive_bkg else np.array([], dtype=np.float32)
    inclusive = threshold_for_target(
        np.r_[np.ones(sig_all.size, dtype=np.int8), np.zeros(bkg_all.size, dtype=np.int8)],
        np.r_[sig_all, bkg_all],
        args.target,
    )

    for clo, chi in zip(cent_edges[:-1], cent_edges[1:]):
        cent_rows: list[dict] = []
        for plo, phi in zip(pt_edges[:-1], pt_edges[1:]):
            sig = np.concatenate(by_cell[(clo, chi, plo, phi, 1)]) if by_cell[(clo, chi, plo, phi, 1)] else np.array([], dtype=np.float32)
            bkg = np.concatenate(by_cell[(clo, chi, plo, phi, 0)]) if by_cell[(clo, chi, plo, phi, 0)] else np.array([], dtype=np.float32)
            item = threshold_for_target(
                np.r_[np.ones(sig.size, dtype=np.int8), np.zeros(bkg.size, dtype=np.int8)],
                np.r_[sig, bkg],
                args.target,
            )
            if item is None:
                continue
            row = {
                "centrality_min": float(clo),
                "centrality_max": float(chi),
                "centrality_label": f"{int(clo)}-{int(chi)}%",
                "pt_min": float(plo),
                "pt_max": float(phi),
                "pt_center": float(0.5 * (plo + phi)),
                "source": "cell",
                **item,
            }
            rows.append(row)
            cent_rows.append(row)

    fits = {}
    for clo, chi in zip(cent_edges[:-1], cent_edges[1:]):
        key = f"{int(clo)}_{int(chi)}"
        cent_rows = [r for r in rows if r["centrality_min"] == clo and r["centrality_max"] == chi]
        fits[key] = {"centrality_label": f"{int(clo)}-{int(chi)}%", **fit_rows(cent_rows)}

    metadata = {
        "schema": "AUAU_BDT_WP80_ET_FIT_CELLS_V1",
        "manifest": str(args.manifest),
        "score_key": score_key,
        "target_signal_efficiency": args.target,
        "pt_edges": pt_edges,
        "cent_edges": cent_edges,
        "files_loaded": files_loaded,
        "rows_loaded": rows_loaded,
        "source_label": args.source_label,
        "model_label": args.model_label,
        "training_sample": args.training_sample,
        "training_inputs": args.training_inputs,
        "inclusive": inclusive or {},
        "fits": fits,
    }
    return {"metadata": metadata, "cells": rows}


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "centrality_min",
        "centrality_max",
        "centrality_label",
        "pt_min",
        "pt_max",
        "pt_center",
        "threshold",
        "signal_efficiency",
        "background_fake_rate",
        "signal_entries",
        "background_entries",
        "source",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path)
    ap.add_argument("--json-in", type=Path)
    ap.add_argument("--json-out", type=Path)
    ap.add_argument("--csv", type=Path)
    ap.add_argument("--score-key", default="score_globalEtCent1535_bdt_noIso")
    ap.add_argument("--pt-edges", default=DEFAULT_PT_EDGES)
    ap.add_argument("--cent-edges", default=DEFAULT_CENT_EDGES)
    ap.add_argument("--target", type=float, default=0.80)
    ap.add_argument("--source-label", default="")
    ap.add_argument("--model-label", default="global EtCent1535 no-iso BDT")
    ap.add_argument("--training-sample", default="Branch A Jet12+20+30+40")
    ap.add_argument("--training-inputs", default="baseV3E + cent + weta33/wphi33")
    args = ap.parse_args()
    if args.json_in is None and args.manifest is None:
        ap.error("--manifest is required unless --json-in is supplied")
    return args


def main() -> None:
    args = parse_args()
    if args.json_in is not None:
        payload = json.loads(args.json_in.read_text())
    else:
        payload = compute(args)
    if args.json_out is not None:
        text = json.dumps(payload, separators=(",", ":"))
        if str(args.json_out) == "-":
            print(text)
        else:
            args.json_out.parent.mkdir(parents=True, exist_ok=True)
            args.json_out.write_text(text + "\n")
    if args.csv is not None:
        write_csv(args.csv, payload["cells"])
        print(args.csv)


if __name__ == "__main__":
    main()
