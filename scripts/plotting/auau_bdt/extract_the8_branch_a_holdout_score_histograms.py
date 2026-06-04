#!/usr/bin/env python3
"""Extract THE-8 Branch A held-out score histograms from training matrices."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


BRANCHES = [
    ("jet12_20", "Jet12+20"),
    ("jet12_20_30", "Jet12+20+30"),
    ("jet12_20_30_40", "Jet12+20+30+40"),
]
CENTRALITY = [("0_20", "0-20%", 0.0, 20.0), ("20_50", "20-50%", 20.0, 50.0), ("50_80", "50-80%", 50.0, 80.0)]
PRODUCT = "globalEtCent1535_bdt_noIso"
WEIGHT_COL = "__ppg12_exact_training_weight"


def safe_ratio(frame: dict[str, np.ndarray], num: str, den: str) -> np.ndarray:
    n = frame[num].astype("float64", copy=False)
    d = frame[den].astype("float64", copy=False)
    out = np.full(len(n), np.nan, dtype="float64")
    good = np.isfinite(n) & np.isfinite(d) & (np.abs(d) > 1.0e-9)
    out[good] = n[good] / d[good]
    return out


def add_needed_derived(frame: dict[str, np.ndarray]) -> None:
    if "cluster_weta_over_wphi" not in frame:
        frame["cluster_weta_over_wphi"] = safe_ratio(frame, "cluster_weta_cogx", "cluster_wphi_cogx")
    if "cluster_weta33_over_wphi33" not in frame:
        frame["cluster_weta33_over_wphi33"] = safe_ratio(frame, "cluster_weta33_cogx", "cluster_wphi33_cogx")


def load_registry(path: Path) -> dict:
    data = json.loads(path.read_text())
    model = data["models"][0]
    report = model["report"]
    return {
        "features": list(model["features"]),
        "pt_range": model.get("pt_range") or report.get("pt_range"),
        "cent_range": model.get("cent_range") or report.get("cent_range"),
        "split": report["split"],
        "reported_auc": float(report["auc"]),
        "reported_holdout_auc": float(report["overfit_diagnostics"]["holdout_auc"]),
        "reported_holdout_rows": int(report["overfit_diagnostics"]["holdout_rows"]),
        "reported_train_rows": int(report["overfit_diagnostics"]["train_rows"]),
        "model_id": model["model_id"],
        "product": model["product"],
    }


def load_matrix(path: Path, features: list[str]) -> dict[str, np.ndarray]:
    z = np.load(path, allow_pickle=True)
    needed = set(features) | {"is_signal", "centrality", "cluster_Et", "source_sample", WEIGHT_COL}
    base_needed = set()
    for name in needed:
        if name == "cluster_weta_over_wphi":
            base_needed.update(["cluster_weta_cogx", "cluster_wphi_cogx"])
        elif name == "cluster_weta33_over_wphi33":
            base_needed.update(["cluster_weta33_cogx", "cluster_wphi33_cogx"])
        else:
            base_needed.add(name)
    missing = sorted(col for col in base_needed if col not in z.files)
    if missing:
        raise SystemExit(f"{path} is missing columns: {missing}")
    frame = {name: z[name] for name in base_needed}
    add_needed_derived(frame)
    return frame


def filtered_indices(frame: dict[str, np.ndarray], features: list[str], pt_range, cent_range) -> np.ndarray:
    mask = np.ones(len(frame["is_signal"]), dtype=bool)
    if pt_range is not None:
        lo, hi = float(pt_range[0]), float(pt_range[1])
        et = frame["cluster_Et"]
        mask &= np.isfinite(et) & (et >= lo) & (et < hi)
    if cent_range is not None:
        lo, hi = float(cent_range[0]), float(cent_range[1])
        cent = frame["centrality"]
        mask &= np.isfinite(cent) & (cent >= lo) & (cent < hi)
    y = frame["is_signal"].astype("int32", copy=False)
    mask &= np.isin(y, [0, 1])
    for name in features:
        mask &= np.isfinite(frame[name])
    return np.flatnonzero(mask)


def holdout_indices(frame: dict[str, np.ndarray], registry: dict, seed: int) -> np.ndarray:
    from sklearn.model_selection import train_test_split

    idx = filtered_indices(frame, registry["features"], registry["pt_range"], registry["cent_range"])
    y = frame["is_signal"][idx].astype("int32", copy=False)
    _, test_idx = train_test_split(idx, test_size=0.10, random_state=seed, stratify=y)
    if len(test_idx) != registry["reported_holdout_rows"]:
        raise SystemExit(
            "Reconstructed holdout row count mismatch: "
            f"got {len(test_idx)}, expected {registry['reported_holdout_rows']}"
        )
    return np.asarray(test_idx, dtype="int64")


def score_model(model_path: Path, frame: dict[str, np.ndarray], features: list[str], idx: np.ndarray, batch_size: int) -> np.ndarray:
    import xgboost as xgb

    booster = xgb.Booster()
    booster.load_model(str(model_path))
    out = np.empty(len(idx), dtype="float32")
    for start in range(0, len(idx), batch_size):
        stop = min(start + batch_size, len(idx))
        local = idx[start:stop]
        x = np.column_stack([frame[name][local].astype("float32", copy=False) for name in features])
        dmat = xgb.DMatrix(x)
        out[start:stop] = booster.predict(dmat).astype("float32")
    return out


def weighted_auc(y: np.ndarray, score: np.ndarray, weight: np.ndarray) -> float:
    from sklearn.metrics import roc_auc_score

    mask = np.isin(y, [0, 1]) & np.isfinite(score) & np.isfinite(weight) & (weight > 0.0)
    if len(np.unique(y[mask])) != 2:
        return math.nan
    return float(roc_auc_score(y[mask], score[mask], sample_weight=weight[mask]))


def weighted_auc_from_masks(signal_mask: np.ndarray, background_mask: np.ndarray, score: np.ndarray, weight: np.ndarray) -> float:
    y = np.full(len(score), -1, dtype="int32")
    y[background_mask] = 0
    y[signal_mask] = 1
    return weighted_auc(y, score, weight)


def hist_payload(signal_mask: np.ndarray, background_mask: np.ndarray, score: np.ndarray, cent: np.ndarray, bins: np.ndarray) -> dict:
    def one(mask: np.ndarray) -> dict:
        values = score[mask & np.isfinite(score)]
        counts, _ = np.histogram(values, bins=bins)
        width = np.diff(bins)
        density = counts.astype("float64")
        if density.sum() > 0:
            density = density / density.sum() / width
        return {"counts": counts.astype(int).tolist(), "density": density.astype(float).tolist(), "entries": int(len(values))}

    out = {
        "signal": one(signal_mask),
        "background": one(background_mask),
        "by_centrality": {},
    }
    for key, label, lo, hi in CENTRALITY:
        cmask = np.isfinite(cent) & (cent >= lo) & (cent < hi)
        out["by_centrality"][key] = {
            "label": label,
            "signal": one(cmask & signal_mask),
            "background": one(cmask & background_mask),
        }
    return out


PLOT_CLASS_DEFINITION = (
    "Signal MC = source_sample contains embeddedPhoton and is_signal == 1; "
    "Inclusive MC = source_sample contains embeddedJet with no truth-background filter"
)


def source_class_masks(source: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    src = source.astype(str)
    signal_mask = (np.char.find(src, "embeddedPhoton") >= 0) & (y.astype("int32", copy=False) == 1)
    background_mask = np.char.find(src, "embeddedJet") >= 0
    overlap = signal_mask & background_mask
    if np.any(overlap):
        raise SystemExit("Source-defined signal/background masks overlap")
    if not np.any(signal_mask) or not np.any(background_mask):
        raise SystemExit(
            "Source-defined plotting classes require both embeddedPhoton and embeddedJet rows; "
            f"signal={int(signal_mask.sum())} background={int(background_mask.sum())}"
        )
    return signal_mask, background_mask


def branch_record(
    label: str,
    summary: dict,
    source: np.ndarray,
    y: np.ndarray,
    score: np.ndarray,
    cent: np.ndarray,
    weight: np.ndarray,
    bins: np.ndarray,
) -> dict:
    signal_mask, background_mask = source_class_masks(source, y)
    h = hist_payload(signal_mask, background_mask, score, cent, bins)
    auc = weighted_auc_from_masks(signal_mask, background_mask, score, weight)
    total = int(len(score))
    sig = int(np.sum(signal_mask))
    bkg = int(np.sum(background_mask))
    return {
        "label": label,
        "bin_edges": bins.astype(float).tolist(),
        "inclusive": {"signal": h["signal"], "background": h["background"]},
        "by_centrality": h["by_centrality"],
        "auc_rows": [],
        "wp80_rows": [],
        "summary": {
            **summary,
            "status": "READY",
            "total_entries": str(total),
            "scored_entries": str(total),
            "signal_entries": str(sig),
            "background_entries": str(bkg),
            f"{PRODUCT}_auc": f"{auc:.6f}",
            f"{PRODUCT}_finite_score_fraction": f"{float(np.mean(np.isfinite(score))):.6f}",
            "finite_score_fraction": f"{float(np.mean(np.isfinite(score))):.6f}",
            "plot_class_definition": PLOT_CLASS_DEFINITION,
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--remote-source-base", type=Path, default=Path("/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining"))
    ap.add_argument("--remote-model-base", type=Path, default=Path("/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models"))
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--random-seed", type=int, default=13)
    ap.add_argument("--batch-size", type=int, default=200000)
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    bins = np.linspace(0.0, 1.0, 51)
    loaded: dict[str, dict] = {}
    for tag, label in BRANCHES:
        model_dir = args.remote_model_base / f"THE8_branchA_{tag}_global_noiso_bdt_mem24_20260527"
        source_dir = args.remote_source_base / f"THE8_branchA_ladder_{tag}_20260527"
        registry = load_registry(model_dir / "model_registry.json")
        frame = load_matrix(model_dir / "training_matrix.npz", registry["features"])
        holdout = holdout_indices(frame, registry, args.random_seed)
        loaded[tag] = {
            "label": label,
            "model_dir": model_dir,
            "source_dir": source_dir,
            "registry": registry,
            "frame": frame,
            "holdout": holdout,
        }

    own_branches = []
    for tag, label in BRANCHES:
        item = loaded[tag]
        reg = item["registry"]
        idx = item["holdout"]
        score = score_model(
            item["model_dir"] / "auau_tight_bdt_globalEtCent1535_bdt_noIso_tmva.xgb.json",
            item["frame"],
            reg["features"],
            idx,
            args.batch_size,
        )
        y = item["frame"]["is_signal"][idx].astype("int32", copy=False)
        cent = item["frame"]["centrality"][idx].astype("float32", copy=False)
        weight = item["frame"][WEIGHT_COL][idx].astype("float64", copy=False)
        source = item["frame"]["source_sample"][idx].astype(str)
        auc = weighted_auc(y, score, weight)
        if not math.isclose(auc, reg["reported_holdout_auc"], rel_tol=0.0, abs_tol=5.0e-6):
            raise SystemExit(f"{tag} own-holdout AUC mismatch: got {auc}, registry {reg['reported_holdout_auc']}")
        own_branches.append(
            branch_record(
                label,
                {
                "validation_mode": "own_10pct_training_holdout_truth_signal_inclusive_jet",
                    "source": str(item["source_dir"]),
                    "model_dir": str(item["model_dir"]),
                    "model_registry": str(item["model_dir"] / "model_registry.json"),
                    "split_mode": "row",
                    "test_fraction_requested": "0.10",
                    "random_seed": str(args.random_seed),
                    "reported_holdout_auc": f"{reg['reported_holdout_auc']:.9f}",
                },
                source,
                y,
                score,
                cent,
                weight,
                bins,
            )
        )

    common = loaded["jet12_20"]
    common_idx = common["holdout"]
    common_branches = []
    for tag, label in BRANCHES:
        model_item = loaded[tag]
        reg = model_item["registry"]
        score = score_model(
            model_item["model_dir"] / "auau_tight_bdt_globalEtCent1535_bdt_noIso_tmva.xgb.json",
            common["frame"],
            reg["features"],
            common_idx,
            args.batch_size,
        )
        y = common["frame"]["is_signal"][common_idx].astype("int32", copy=False)
        cent = common["frame"]["centrality"][common_idx].astype("float32", copy=False)
        weight = common["frame"][WEIGHT_COL][common_idx].astype("float64", copy=False)
        source = common["frame"]["source_sample"][common_idx].astype(str)
        common_branches.append(
            branch_record(
                label,
                {
                    "validation_mode": "common_jet12_20_10pct_training_holdout_truth_signal_inclusive_jet",
                    "source": str(common["source_dir"]),
                    "model_dir": str(model_item["model_dir"]),
                    "model_registry": str(model_item["model_dir"] / "model_registry.json"),
                    "split_mode": "row",
                    "test_fraction_requested": "0.10",
                    "random_seed": str(args.random_seed),
                    "common_holdout_sample": "Jet12+20",
                },
                source,
                y,
                score,
                cent,
                weight,
                bins,
            )
        )

    outputs = {
        "the8_branch_a_ladder_own_holdout_score_histograms.json": {
            "schema": "THE8_BRANCH_A_LADDER_COMPACT_HOLDOUT_SCORE_HISTOGRAMS_V1",
            "description": "Each BDT scored on its own 10% row holdout from the matching Branch A training matrix.",
            "branches": own_branches,
        },
        "the8_branch_a_ladder_common_jet12_20_holdout_score_histograms.json": {
            "schema": "THE8_BRANCH_A_LADDER_COMPACT_HOLDOUT_SCORE_HISTOGRAMS_V1",
            "description": "All three BDTs scored on the same Jet12+20 10% row holdout from the Branch A training matrix.",
            "branches": common_branches,
        },
    }
    for name, payload in outputs.items():
        path = args.outdir / name
        path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
