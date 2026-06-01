#!/usr/bin/env python3
"""Summarize local AuAu ML training smoke outputs.

This is a local QA helper for directories pulled with:

  ./scripts/sftp_get_recoiljets_outputs.sh trainingLatest <tight|npb|jetResidual>

It reports whether the expected training trees are present, basic row/class
counts, centrality/radius coverage, and any model summary JSON metrics.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def find_default_base(task: str) -> Path:
    repo = Path(__file__).resolve().parents[1]
    aliases = {
        "tight": "tight",
        "npb": "npb",
        "jetResidual": "jetResidual",
        "jetML": "jetResidual",
    }
    key = aliases.get(task, task)
    return repo / "InputFiles" / "trainingSmoke" / key


def finite_stats(values):
    import numpy as np

    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) == 0:
        return "none"
    return f"n={len(arr)} mean={float(np.mean(arr)):.4g} rms={float(math.sqrt(np.mean(arr * arr))):.4g} min={float(np.min(arr)):.4g} max={float(np.max(arr)):.4g}"


def summarize_root(path: Path):
    try:
        import numpy as np
        import uproot
    except Exception as exc:
        raise SystemExit("This audit helper needs uproot and numpy in the active Python env.") from exc

    print(f"\nROOT: {path}")
    with uproot.open(path) as f:
        keys = set(f.keys())
        if "AuAuPhotonIDTrainingTree;1" in keys or "AuAuPhotonIDTrainingTree" in f:
            t = f["AuAuPhotonIDTrainingTree"]
            print(f"  AuAuPhotonIDTrainingTree entries={t.num_entries}")
            cols = [c for c in ["is_signal", "centrality", "cluster_Et", "event_weight"] if c in t.keys()]
            if cols:
                arr = t.arrays(cols, library="np")
                if "is_signal" in arr:
                    y = arr["is_signal"]
                    print(f"    class balance: signal={int(np.sum(y == 1))} background={int(np.sum(y == 0))}")
                if "centrality" in arr:
                    print(f"    centrality: {finite_stats(arr['centrality'])}")
                if "cluster_Et" in arr:
                    print(f"    cluster_Et: {finite_stats(arr['cluster_Et'])}")
        if "JetResidualMLTrainingTree;1" in keys or "JetResidualMLTrainingTree" in f:
            t = f["JetResidualMLTrainingTree"]
            print(f"  JetResidualMLTrainingTree entries={t.num_entries}")
            cols = [
                c for c in [
                    "is_recoil", "is_matched", "centrality", "R", "rKey",
                    "reco_areaSub_pt", "truth_pt", "delta_pt", "response_ratio",
                ] if c in t.keys()
            ]
            if cols:
                arr = t.arrays(cols, library="np")
                if "is_recoil" in arr:
                    print(f"    recoil rows: {int(np.sum(arr['is_recoil'] == 1))}")
                if "is_matched" in arr:
                    print(f"    matched rows: {int(np.sum(arr['is_matched'] == 1))}")
                if "R" in arr:
                    radii = sorted(set(float(x) for x in arr["R"] if np.isfinite(x)))
                    print(f"    radii: {radii}")
                if "centrality" in arr:
                    print(f"    centrality: {finite_stats(arr['centrality'])}")
                if "delta_pt" in arr:
                    matched = arr.get("is_matched", np.ones(len(arr["delta_pt"]), dtype=bool)) == 1
                    print(f"    delta_pt matched: {finite_stats(arr['delta_pt'][matched])}")
                if "response_ratio" in arr:
                    matched = arr.get("is_matched", np.ones(len(arr["response_ratio"]), dtype=bool)) == 1
                    print(f"    response_ratio matched: {finite_stats(arr['response_ratio'][matched])}")


def summarize_json(path: Path):
    try:
        data = json.loads(path.read_text())
    except Exception:
        return
    print(f"\nJSON: {path}")
    if isinstance(data, list):
        for item in data:
            name = Path(item.get("output_tmva", "model")).name
            auc = item.get("auc", None)
            rows = item.get("n_rows", None)
            sig = item.get("n_signal", None)
            bkg = item.get("n_background", None)
            cent = item.get("cent_range", None)
            print(f"  {name}: cent={cent} rows={rows} sig={sig} bkg={bkg} auc={auc}")
    elif isinstance(data, dict):
        if "metrics" in data:
            metrics = data["metrics"]
            print(f"  model={data.get('tmva_model')}")
            print(f"  rows_loaded={data.get('rows_loaded')} rows_used={data.get('rows_used')}")
            print(
                "  residual RMS: "
                f"area={metrics.get('area_residual_rms')} ml={metrics.get('ml_residual_rms')} "
                f"mean(area)={metrics.get('area_residual_mean')} mean(ml)={metrics.get('ml_residual_mean')}"
            )
        elif "auc" in data:
            print(f"  model={data.get('output_tmva')} rows={data.get('n_rows')} auc={data.get('auc')}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task", nargs="?", default="jetResidual",
                        help="tight, npb, jetResidual, or an explicit directory with --path")
    parser.add_argument("--path", type=Path, default=None)
    args = parser.parse_args()

    base = args.path if args.path else find_default_base(args.task)
    if not base.exists():
        raise SystemExit(f"Path does not exist: {base}")

    print(f"Audit base: {base}")
    roots = sorted(base.rglob("*.root"))
    jsons = sorted(base.rglob("*.json"))
    print(f"ROOT files: {len(roots)}")
    print(f"JSON files: {len(jsons)}")

    for path in roots[:50]:
        summarize_root(path)
    if len(roots) > 50:
        print(f"\n[WARN] Skipped {len(roots) - 50} additional ROOT files")

    for path in jsons:
        summarize_json(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
