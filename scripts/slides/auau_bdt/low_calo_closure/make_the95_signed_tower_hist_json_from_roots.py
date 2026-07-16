#!/usr/bin/env python3
"""Build THE-95 signed-tower low-calo event histograms from extraction ROOTs.

This helper is meant to run in the sPHENIX/ROOT environment, usually on SDCC,
against the `training_roots.list` produced by the signed-tower diagnostic
extraction. It writes the same compact histogram JSON schema consumed by the
existing Blair kBird and THE58 slide builders, without copying large ROOT files
back to the laptop.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import ROOT


SAMPLES = [
    "run28_embeddedPhoton12",
    "run28_embeddedPhoton20",
    "run28_embeddedJet12",
    "run28_embeddedJet20",
    "run28_embeddedJet30",
    "run28_embeddedJet40",
]


def sample_from_path(path: str) -> str:
    for part in Path(path).parts:
        if part in SAMPLES:
            return part
    for sample in SAMPLES:
        if sample in path:
            return sample
    return "UNKNOWN"


def read_manifest(path: Path, max_files: int | None) -> list[str]:
    files = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if max_files is not None:
        files = files[:max_files]
    if not files:
        raise SystemExit(f"No ROOT files listed in {path}")
    return files


def root_numpy_for_file(path: str, branches: list[str]) -> dict[str, np.ndarray]:
    rdf = ROOT.RDataFrame("AuAuPhotonIDTrainingTree", path)
    arrays = rdf.AsNumpy(branches)
    return {key: np.asarray(value) for key, value in arrays.items()}


def per_file_unique_events(path: str, branches: list[str]) -> tuple[int, np.ndarray, np.ndarray, np.ndarray]:
    arrays = root_numpy_for_file(path, branches)
    n_rows = len(arrays["run"])
    if n_rows == 0:
        return 0, np.array([], dtype="float32"), np.array([], dtype="float32"), np.array([], dtype="int64")

    key = np.empty(n_rows, dtype=[("run", "i4"), ("evt", "i8")])
    key["run"] = arrays["run"].astype("int32", copy=False)
    key["evt"] = arrays["evt"].astype("int64", copy=False)
    unique_idx = np.unique(key, return_index=True)[1]

    cent = arrays["centrality"].astype("float32", copy=False)[unique_idx]
    if "event_calo_log10_total_energy_plus1" in arrays:
        log_calo = arrays["event_calo_log10_total_energy_plus1"].astype("float32", copy=False)[unique_idx]
    else:
        total = arrays["event_calo_total_energy"].astype("float64", copy=False)[unique_idx]
        log_calo = np.log10(np.maximum(total, 0.0) + 1.0).astype("float32")
    evt = arrays["evt"].astype("int64", copy=False)[unique_idx]
    good = np.isfinite(cent) & np.isfinite(log_calo) & (cent >= 0.0) & (cent < 80.0)
    return n_rows, cent[good], log_calo[good], evt[good]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--max-files", type=int, default=None)
    parser.add_argument("--progress-every", type=int, default=100)
    parser.add_argument("--energy-min", type=float, default=2.70)
    parser.add_argument("--energy-max", type=float, default=3.45)
    parser.add_argument("--energy-bin", type=float, default=0.01)
    parser.add_argument(
        "--schema",
        default="THE95_SIGNED_TOWER_LOW_CALO_FULL_COUNT_HISTOGRAMS_V1",
        help="Schema/variant label to record in the compact JSON.",
    )
    parser.add_argument(
        "--cut-variable",
        default="log10(E_CEMC + E_IHCal + E_OHCal + 1), signed finite tower sums",
        help="Human-readable event-calo definition recorded in the JSON.",
    )
    parser.add_argument(
        "--variant-note",
        default="",
        help="Optional exact processing note, e.g. CaloTowerStatus/get_isGood toggles.",
    )
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    files = read_manifest(args.manifest, args.max_files)
    branches = ["run", "evt", "centrality", "event_calo_log10_total_energy_plus1", "event_calo_total_energy"]

    cent_edges = np.arange(0.0, 85.0, 5.0, dtype="float64")
    energy_edges = np.round(np.arange(args.energy_min, args.energy_max + 0.0001, args.energy_bin), 10)
    total_hist = np.zeros((len(cent_edges) - 1, len(energy_edges) - 1), dtype="int64")
    source_total = {sample: np.zeros(len(cent_edges) - 1, dtype="int64") for sample in SAMPLES}
    source_total.setdefault("UNKNOWN", np.zeros(len(cent_edges) - 1, dtype="int64"))

    raw_rows = 0
    unique_events = 0
    low_0_5_counts = {700.0: 0, 725.0: 0, 750.0: 0}
    low_0_5_total = 0

    for idx, path in enumerate(files, 1):
        if idx == 1 or idx % max(args.progress_every, 1) == 0 or idx == len(files):
            print(f"[THE95] reading {idx}/{len(files)} {path}", flush=True)
        sample = sample_from_path(path)
        n_rows, cent, log_calo, _evt = per_file_unique_events(path, branches)
        raw_rows += int(n_rows)
        unique_events += int(len(cent))
        if len(cent) == 0:
            continue
        hist, _, _ = np.histogram2d(cent, log_calo, bins=[cent_edges, energy_edges])
        total_hist += hist.astype("int64")
        cent_bin = np.floor(cent / 5.0).astype("int32")
        valid = (cent_bin >= 0) & (cent_bin < len(cent_edges) - 1)
        source_total[sample] += np.bincount(cent_bin[valid], minlength=len(cent_edges) - 1).astype("int64")

        central = (cent >= 0.0) & (cent < 5.0)
        if np.any(central):
            e_calo = np.power(10.0, log_calo[central].astype("float64")) - 1.0
            low_0_5_total += int(e_calo.size)
            for thr in low_0_5_counts:
                low_0_5_counts[thr] += int(np.sum(e_calo < thr))

    panels = []
    thresholds = []
    for i in range(len(cent_edges) - 1):
        lo = float(cent_edges[i])
        hi = float(cent_edges[i + 1])
        full = total_hist[i].astype(int).tolist()
        total = int(total_hist[i].sum())
        panels.append(
            {
                "cent_lo": lo,
                "cent_hi": hi,
                "threshold": float(energy_edges[0]),
                "event_total": total,
                "event_retained": total,
                "event_removed": 0,
                "event_removed_fraction": 0.0,
                "retained_hist": full,
                "removed_hist": [0 for _ in full],
                "source_event_total": {sample: int(source_total[sample][i]) for sample in sorted(source_total)},
                "source_event_removed": {sample: 0 for sample in sorted(source_total)},
            }
        )
        thresholds.append({"cent_lo": lo, "cent_hi": hi, "threshold": float(energy_edges[0])})

    payload = {
        "schema": args.schema,
        "source_cache_list": str(args.manifest),
        "cache_count": len(files),
        "raw_candidate_rows": raw_rows,
        "event_rows_after_per_cache_dedup": unique_events,
        "centrality_domain": "0 <= centrality < 80",
        "cut_variable": args.cut_variable,
        "variant_note": args.variant_note,
        "energy_edges": [float(x) for x in energy_edges],
        "thresholds": thresholds,
        "panels": panels,
        "low_energy_rate_0_5": {
            str(int(thr)): {
                "count": int(count),
                "total": int(low_0_5_total),
                "fraction": float(count / low_0_5_total) if low_0_5_total else math.nan,
            }
            for thr, count in sorted(low_0_5_counts.items())
        },
    }

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"[THE95] wrote {args.out_json}")
    print(json.dumps(payload["low_energy_rate_0_5"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
