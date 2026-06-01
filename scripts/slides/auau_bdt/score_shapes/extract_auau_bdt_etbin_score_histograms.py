#!/usr/bin/env python3
"""Extract ET-binned Au+Au BDT score histograms from validation score caches."""

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
from pathlib import Path

import numpy as np


DEFAULT_ET_EDGES = "15,17,19,21,23,25,27,30,35"


def parse_edges(spec: str) -> list[float]:
    values = [float(tok) for tok in spec.split(",") if tok.strip()]
    if len(values) < 2:
        raise ValueError(f"Need at least two bin edges: {spec}")
    if any(values[i + 1] <= values[i] for i in range(len(values) - 1)):
        raise ValueError(f"Bin edges must be strictly increasing: {spec}")
    return values


def auc_from_scores(y: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y, dtype=np.int8)
    score = np.asarray(score, dtype=np.float64)
    mask = np.isfinite(score) & np.isin(y, [0, 1])
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


def manifest_paths(path: Path) -> list[Path]:
    paths = [Path(line.strip()) for line in path.read_text().splitlines() if line.strip()]
    if not paths:
        raise FileNotFoundError(f"No score caches listed in {path}")
    return paths


def compute_payload(args: argparse.Namespace) -> dict:
    et_edges = parse_edges(args.et_bins)
    score_edges = np.linspace(args.score_min, args.score_max, args.score_bins + 1)
    et_bins = list(zip(et_edges[:-1], et_edges[1:]))
    hist = {(lo, hi, cls): np.zeros(args.score_bins, dtype=np.int64) for lo, hi in et_bins for cls in (0, 1)}
    scores = {(lo, hi, cls): [] for lo, hi in et_bins for cls in (0, 1)}
    rows_loaded = 0
    files_loaded = 0

    for cache in manifest_paths(args.manifest):
        with np.load(cache, allow_pickle=True) as z:
            missing = [key for key in ("is_signal", "centrality", "cluster_Et", args.score_key) if key not in z.files]
            if missing:
                raise KeyError(f"{cache} missing keys {missing}")
            y = np.asarray(z["is_signal"], dtype=np.int8)
            cent = np.asarray(z["centrality"], dtype=np.float32)
            et = np.asarray(z["cluster_Et"], dtype=np.float32)
            score = np.asarray(z[args.score_key], dtype=np.float32)
            base = (
                np.isfinite(cent)
                & np.isfinite(et)
                & np.isfinite(score)
                & np.isin(y, [0, 1])
                & (cent >= args.centrality_min)
                & (cent < args.centrality_max)
                & (et >= et_edges[0])
                & (et < et_edges[-1])
            )
            rows_loaded += int(np.sum(base))
            files_loaded += 1
            if not np.any(base):
                continue
            for lo, hi in et_bins:
                et_mask = base & (et >= lo) & (et < hi)
                if not np.any(et_mask):
                    continue
                for cls in (0, 1):
                    selected = score[et_mask & (y == cls)]
                    if selected.size == 0:
                        continue
                    hist[(lo, hi, cls)] += np.histogram(selected, bins=score_edges)[0]
                    scores[(lo, hi, cls)].append(selected.astype(np.float32, copy=False))

    rows = []
    summary = []
    for lo, hi in et_bins:
        sig = np.concatenate(scores[(lo, hi, 1)]) if scores[(lo, hi, 1)] else np.array([], dtype=np.float32)
        bkg = np.concatenate(scores[(lo, hi, 0)]) if scores[(lo, hi, 0)] else np.array([], dtype=np.float32)
        y_all = np.r_[np.ones(sig.size, dtype=np.int8), np.zeros(bkg.size, dtype=np.int8)]
        score_all = np.r_[sig, bkg]
        auc = auc_from_scores(y_all, score_all)
        sig_entries = int(sig.size)
        bkg_entries = int(bkg.size)
        label = f"{int(lo)}-{int(hi)}"
        summary.append(
            {
                "et_label": label,
                "et_lo": lo,
                "et_hi": hi,
                "entries": sig_entries + bkg_entries,
                "signal_entries": sig_entries,
                "background_entries": bkg_entries,
                "auc": auc,
            }
        )
        for cls, cls_name, entries in ((1, "signal", sig_entries), (0, "background", bkg_entries)):
            counts = hist[(lo, hi, cls)]
            density = counts / max(1, entries) / np.diff(score_edges)
            for bin_low, bin_high, den, cnt in zip(score_edges[:-1], score_edges[1:], density, counts):
                rows.append(
                    {
                        "et_label": label,
                        "et_lo": lo,
                        "et_hi": hi,
                        "class": cls_name,
                        "bin_low": float(bin_low),
                        "bin_high": float(bin_high),
                        "density": float(den),
                        "count": int(cnt),
                        "entries": entries,
                        "auc": auc,
                        "signal_entries": sig_entries,
                        "background_entries": bkg_entries,
                    }
                )

    metadata = {
        "schema": args.schema,
        "manifest": str(args.manifest),
        "score_key": args.score_key,
        "et_bins": [[lo, hi] for lo, hi in et_bins],
        "score_bins": args.score_bins,
        "score_range": [args.score_min, args.score_max],
        "selection": args.selection,
        "source_label": args.source_label,
        "files_loaded": files_loaded,
        "rows_loaded": rows_loaded,
        "summary": summary,
        "slide_title": args.slide_title,
        "slide_subtitle": args.slide_subtitle,
        "internal_subtitle": args.internal_subtitle,
        "model_label": args.model_label,
        "inputs_title": args.inputs_title,
        "inputs_label": args.inputs_label,
        "training_title": args.training_title,
        "training_label": args.training_label,
        "readout_note": args.readout_note.format(rows_loaded=rows_loaded),
    }
    return {"metadata": metadata, "rows": rows}


def write_outputs(payload: dict, csv_path: Path, metadata_path: Path) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    rows = payload["rows"]
    if not rows:
        raise ValueError("No histogram rows to write")
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    metadata_path.write_text(json.dumps(payload["metadata"], indent=2, sort_keys=True) + "\n")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path)
    ap.add_argument("--score-key", default="score_globalEtCent1535_bdt_noIso")
    ap.add_argument("--et-bins", default=DEFAULT_ET_EDGES)
    ap.add_argument("--score-bins", type=int, default=34)
    ap.add_argument("--score-min", type=float, default=0.0)
    ap.add_argument("--score-max", type=float, default=1.0)
    ap.add_argument("--centrality-min", type=float, default=0.0)
    ap.add_argument("--centrality-max", type=float, default=80.0)
    ap.add_argument("--csv", type=Path)
    ap.add_argument("--metadata", type=Path)
    ap.add_argument("--json-out", type=Path)
    ap.add_argument("--json-in", type=Path)
    ap.add_argument("--schema", default="AUAU_GLOBAL_BDT_ETBIN_SCORE_SEPARATION_V2")
    ap.add_argument("--selection", default="15 <= cluster_Et < 35 GeV, 0 <= centrality < 80")
    ap.add_argument("--source-label", default="")
    ap.add_argument(
        "--slide-title",
        default=r"Jet12+20+30+40 BDT keeps strong separation across $15<E_T<35$ GeV",
    )
    ap.add_argument(
        "--slide-subtitle",
        default=r"Branch A $15 < E_T < 35$ GeV Au+Au BDT, evaluated on full score-cache validation; each panel is unit-area signal vs background.",
    )
    ap.add_argument("--internal-subtitle", default="PYTHIA8 Au+Au embedded validation")
    ap.add_argument("--model-label", default=r"global $\bf{EtCent1535}$ no-iso BDT")
    ap.add_argument("--inputs-title", default="Training sample")
    ap.add_argument("--inputs-label", default=r"Branch A $\bf{Jet12+20+30+40}$")
    ap.add_argument("--training-title", default="Training inputs")
    ap.add_argument("--training-label", default="baseV3E + cent + weta33/wphi33")
    ap.add_argument(
        "--readout-note",
        default=(
            r"Full score-cache validation has {rows_loaded:,} scored rows, "
            r"$15 \leq E_T < 35$ GeV and 0-80% centrality; score key is "
            "score_globalEtCent1535_bdt_noIso."
        ),
    )
    args = ap.parse_args()
    if args.json_in is None and args.manifest is None:
        ap.error("--manifest is required unless --json-in is provided")
    if args.json_in is not None and (args.csv is None or args.metadata is None):
        ap.error("--csv and --metadata are required with --json-in")
    if args.json_in is None and args.json_out is None and (args.csv is None or args.metadata is None):
        ap.error("provide --json-out or both --csv and --metadata")
    return args


def main() -> None:
    args = parse_args()
    if args.json_in is not None:
        payload = json.loads(args.json_in.read_text())
    else:
        payload = compute_payload(args)
        if args.json_out is not None:
            text = json.dumps(payload, separators=(",", ":"))
            if str(args.json_out) == "-":
                print(text)
            else:
                args.json_out.parent.mkdir(parents=True, exist_ok=True)
                args.json_out.write_text(text + "\n")
    if args.csv is not None and args.metadata is not None:
        write_outputs(payload, args.csv, args.metadata)
        print(args.csv)
        print(args.metadata)


if __name__ == "__main__":
    main()
