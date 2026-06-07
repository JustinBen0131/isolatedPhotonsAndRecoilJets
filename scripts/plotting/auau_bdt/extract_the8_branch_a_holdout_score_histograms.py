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
COMMON_VALIDATION_SAMPLES = [
    ("jet12_20", "Jet12+20"),
    ("jet12_20_30_40", "Jet12+20+30+40"),
]
CENTRALITY = [("0_20", "0-20%", 0.0, 20.0), ("20_50", "20-50%", 20.0, 50.0), ("50_80", "50-80%", 50.0, 80.0)]
PRODUCT = "globalEtCent1535_bdt_noIso"
SCORE_COL = f"score_{PRODUCT}"
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
    mask = np.isin(y, [0, 1]) & np.isfinite(score) & np.isfinite(weight) & (weight > 0.0)
    if len(np.unique(y[mask])) != 2:
        return math.nan
    yy = y[mask].astype("int32", copy=False)
    ss = score[mask].astype("float64", copy=False)
    ww = weight[mask].astype("float64", copy=False)
    try:
        from sklearn.metrics import roc_auc_score

        return float(roc_auc_score(yy, ss, sample_weight=ww))
    except Exception:
        return weighted_auc_rank(yy, ss, ww)


def weighted_auc_rank(y: np.ndarray, score: np.ndarray, weight: np.ndarray) -> float:
    order = np.argsort(score, kind="mergesort")
    y = y[order]
    score = score[order]
    weight = weight[order]
    pos_total = float(np.sum(weight[y == 1]))
    neg_total = float(np.sum(weight[y == 0]))
    if pos_total <= 0.0 or neg_total <= 0.0:
        return math.nan
    wins = 0.0
    neg_below = 0.0
    start = 0
    n = len(score)
    while start < n:
        stop = start + 1
        while stop < n and score[stop] == score[start]:
            stop += 1
        y_group = y[start:stop]
        w_group = weight[start:stop]
        pos_group = float(np.sum(w_group[y_group == 1]))
        neg_group = float(np.sum(w_group[y_group == 0]))
        wins += pos_group * neg_below + 0.5 * pos_group * neg_group
        neg_below += neg_group
        start = stop
    return wins / (pos_total * neg_total)


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
    "Signal MC = source_sample contains embeddedPhoton and is_signal == 1 "
    "(truth-isolated prompt label); "
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
    centrality_summary = {}
    for key, _label, lo, hi in CENTRALITY:
        cmask = np.isfinite(cent) & (cent >= lo) & (cent < hi)
        sig_c = signal_mask & cmask
        bkg_c = background_mask & cmask
        centrality_summary[f"{PRODUCT}_auc_{key}"] = f"{weighted_auc_from_masks(sig_c, bkg_c, score, weight):.6f}"
        centrality_summary[f"signal_entries_{key}"] = str(int(np.sum(sig_c)))
        centrality_summary[f"background_entries_{key}"] = str(int(np.sum(bkg_c)))
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
            **centrality_summary,
            f"{PRODUCT}_finite_score_fraction": f"{float(np.mean(np.isfinite(score))):.6f}",
            "finite_score_fraction": f"{float(np.mean(np.isfinite(score))):.6f}",
            "plot_class_definition": PLOT_CLASS_DEFINITION,
        },
    }


def branch_record_from_score_caches(
    label: str,
    summary: dict,
    model_dir: Path,
    report_dir: Path,
    bins: np.ndarray,
    *,
    matrix_dir: Path | None = None,
    keep_indices: np.ndarray | None = None,
) -> dict:
    validation_matrix_dir = matrix_dir or model_dir
    matrix_path = validation_matrix_dir / "training_matrix.npz"
    cache_list = report_dir / "score_caches.list"
    if cache_list.exists():
        cache_paths = [Path(line.strip()) for line in cache_list.read_text().splitlines() if line.strip()]
    else:
        cache_paths = sorted((report_dir / "score_caches").glob("score_cache_*.npz"))
    if not cache_paths:
        raise SystemExit(f"No score caches found in {report_dir}")

    matrix = np.load(matrix_path, allow_pickle=True)
    missing = [name for name in ["source_sample", "is_signal", "centrality", WEIGHT_COL] if name not in matrix.files]
    if missing:
        raise SystemExit(f"{matrix_path} is missing columns needed for source-defined full-sample plotting: {missing}")
    source_all = matrix["source_sample"].astype(str)
    y_all = matrix["is_signal"].astype("int32", copy=False)
    cent_all = matrix["centrality"].astype("float32", copy=False)
    weight_all = matrix[WEIGHT_COL].astype("float64", copy=False)
    keep_row_mask = None
    if keep_indices is not None:
        keep_row_mask = np.zeros(len(source_all), dtype=bool)
        keep_row_mask[np.asarray(keep_indices, dtype="int64")] = True

    source_parts = []
    y_parts = []
    score_parts = []
    cent_parts = []
    weight_parts = []
    total_entries = 0
    scored_entries = 0
    offset = 0
    for cache_path in cache_paths:
        cache = np.load(cache_path, allow_pickle=True)
        if SCORE_COL not in cache.files:
            raise SystemExit(f"{cache_path} is missing {SCORE_COL}")
        score = cache[SCORE_COL].astype("float32", copy=False)
        n = len(score)
        stop = offset + n
        if stop > len(source_all):
            raise SystemExit(f"Score cache rows exceed training matrix length for {label}: stop={stop} matrix={len(source_all)}")
        y = cache["is_signal"].astype("int32", copy=False)
        if not np.array_equal(y, y_all[offset:stop]):
            raise SystemExit(f"Score cache/training-matrix is_signal order mismatch at {cache_path}")
        source = source_all[offset:stop]
        cent = cent_all[offset:stop]
        weight = weight_all[offset:stop]
        finite = np.isfinite(score)
        row_keep = np.ones(n, dtype=bool) if keep_row_mask is None else keep_row_mask[offset:stop]
        signal = (np.char.find(source, "embeddedPhoton") >= 0) & (y == 1) & finite & row_keep
        background = (np.char.find(source, "embeddedJet") >= 0) & finite & row_keep
        keep = signal | background
        source_parts.append(source[keep])
        y_parts.append(y[keep])
        score_parts.append(score[keep])
        cent_parts.append(cent[keep])
        weight_parts.append(weight[keep])
        total_entries += n
        scored_entries += int(np.sum(finite & row_keep))
        offset = stop
    if offset != len(source_all):
        raise SystemExit(f"Score cache rows do not cover full training matrix for {label}: cache={offset} matrix={len(source_all)}")

    record = branch_record(
        label,
        summary,
        np.concatenate(source_parts),
        np.concatenate(y_parts),
        np.concatenate(score_parts),
        np.concatenate(cent_parts),
        np.concatenate(weight_parts),
        bins,
    )
    selected_entries = int(np.sum(keep_row_mask)) if keep_row_mask is not None else total_entries
    finite_fraction = float(scored_entries / selected_entries) if selected_entries else math.nan
    record["summary"]["total_entries"] = str(selected_entries)
    record["summary"]["scored_entries"] = str(scored_entries)
    record["summary"][f"{PRODUCT}_finite_score_fraction"] = f"{finite_fraction:.6f}"
    record["summary"]["finite_score_fraction"] = f"{finite_fraction:.6f}"
    record["summary"]["score_cache_list"] = str(cache_list)
    record["summary"]["validation_matrix"] = str(matrix_path)
    if keep_row_mask is not None:
        record["summary"]["score_cache_row_filter"] = "source_10pct_holdout_indices"
    return record


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--remote-source-base", type=Path, default=Path("/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining"))
    ap.add_argument("--remote-model-base", type=Path, default=Path("/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models"))
    ap.add_argument("--outdir", type=Path)
    ap.add_argument("--random-seed", type=int, default=13)
    ap.add_argument("--batch-size", type=int, default=200000)
    ap.add_argument(
        "--stdout-full-own",
        action="store_true",
        help="Print only the own full allotted training+holdout sample payload as JSON.",
    )
    ap.add_argument(
        "--stdout-full-common-grid",
        action="store_true",
        help="Print all three BDTs on the two fixed full score-cache validation samples as JSON.",
    )
    ap.add_argument(
        "--stdout-holdout-common-grid-from-scorecache",
        action="store_true",
        help="Print all three BDTs on the two fixed 10% source holdouts using full score-cache reports as score sources.",
    )
    ap.add_argument(
        "--stdout-holdout-three-by-three-direct",
        action="store_true",
        help="Print all three BDTs on all three fixed 10% source holdouts using direct model scoring.",
    )
    args = ap.parse_args()

    if args.outdir is None and not (
        args.stdout_full_own
        or args.stdout_full_common_grid
        or args.stdout_holdout_common_grid_from_scorecache
        or args.stdout_holdout_three_by_three_direct
    ):
        ap.error("--outdir is required unless a stdout extraction mode is set")
    if args.outdir is not None:
        args.outdir.mkdir(parents=True, exist_ok=True)
    bins = np.linspace(0.0, 1.0, 51)
    if args.stdout_full_own:
        full_own_branches = []
        for tag, label in BRANCHES:
            model_dir = args.remote_model_base / f"THE8_branchA_{tag}_global_noiso_bdt_mem24_20260527"
            source_dir = args.remote_source_base / f"THE8_branchA_ladder_{tag}_20260527"
            report_dir = source_dir / "reports" / f"model_validation_condor_THE8_branchA_{tag}_scorecache_fullstat_20260527"
            registry = load_registry(model_dir / "model_registry.json")
            full_own_branches.append(
                branch_record_from_score_caches(
                    label,
                    {
                        "validation_mode": "own_full_allotted_sample_truth_signal_inclusive_jet",
                        "source": str(source_dir),
                        "model_dir": str(model_dir),
                        "model_registry": str(model_dir / "model_registry.json"),
                        "report_dir": str(report_dir),
                        "split_mode": "full_training_matrix",
                        "row_scope": "training_plus_holdout",
                        "test_fraction_requested": "0.10",
                        "random_seed": str(args.random_seed),
                        "reported_train_rows": str(registry["reported_train_rows"]),
                        "reported_holdout_rows": str(registry["reported_holdout_rows"]),
                    },
                    model_dir,
                    report_dir,
                    bins,
                )
            )
        print(
            json.dumps(
                {
                    "schema": "THE8_BRANCH_A_LADDER_COMPACT_FULL_SAMPLE_SCORE_HISTOGRAMS_V1",
                    "description": "Each BDT scored on its own full allotted Branch A training matrix: training rows plus the 10% row holdout.",
                    "plot_class_definition": PLOT_CLASS_DEFINITION,
                    "branches": full_own_branches,
                },
                indent=2,
                sort_keys=True,
            )
        )
        return 0

    if args.stdout_full_common_grid:
        common_samples = []
        model_dirs = {
            tag: args.remote_model_base / f"THE8_branchA_{tag}_global_noiso_bdt_mem24_20260527"
            for tag, _ in BRANCHES
        }
        source_dirs = {
            tag: args.remote_source_base / f"THE8_branchA_ladder_{tag}_20260527"
            for tag, _ in BRANCHES
        }
        common_report_tags = {
            ("jet12_20", "jet12_20"): "THE8_branchA_jet12_20_scorecache_fullstat_20260527",
            ("jet12_20", "jet12_20_30"): "THE8_commonJet12_20_modelJet12_20_30_scorecache_fullstat_20260603",
            ("jet12_20", "jet12_20_30_40"): "THE8_commonJet12_20_modelJet12_20_30_40_scorecache_fullstat_20260603",
            ("jet12_20_30_40", "jet12_20"): "THE8_commonJet12_20_30_40_modelJet12_20_scorecache_fullstat_20260603",
            ("jet12_20_30_40", "jet12_20_30"): "THE8_commonJet12_20_30_40_modelJet12_20_30_scorecache_fullstat_20260603",
            ("jet12_20_30_40", "jet12_20_30_40"): "THE8_branchA_jet12_20_30_40_scorecache_fullstat_20260527",
        }
        for common_tag, common_label in COMMON_VALIDATION_SAMPLES:
            source_dir = source_dirs[common_tag]
            source_matrix_dir = model_dirs[common_tag]
            branches = []
            for model_tag, model_label in BRANCHES:
                model_dir = model_dirs[model_tag]
                report_tag = common_report_tags[(common_tag, model_tag)]
                report_dir = source_dir / "reports" / f"model_validation_condor_{report_tag}"
                branches.append(
                    branch_record_from_score_caches(
                        model_label,
                        {
                            "validation_mode": f"common_{common_tag}_full_scorecache_truth_signal_inclusive_jet",
                            "source": str(source_dir),
                            "model_dir": str(model_dir),
                            "model_registry": str(model_dir / "model_registry.json"),
                            "report_dir": str(report_dir),
                            "split_mode": "full_training_matrix",
                            "row_scope": "training_plus_holdout",
                            "test_fraction_requested": "0.10",
                            "random_seed": str(args.random_seed),
                            "common_validation_sample": common_label,
                            "model_training_sample": model_label,
                        },
                        model_dir,
                        report_dir,
                        bins,
                        matrix_dir=source_matrix_dir,
                    )
                )
            common_samples.append(
                {
                    "validation_sample": common_label,
                    "validation_tag": common_tag,
                    "validation_matrix": str(source_matrix_dir / "training_matrix.npz"),
                    "branches": branches,
                }
            )
        print(
            json.dumps(
                {
                    "schema": "THE8_BRANCH_A_LADDER_COMPACT_FIXED_FULL_SAMPLE_SCORE_HISTOGRAMS_V1",
                    "description": "All three BDTs scored on each fixed full score-cache validation sample: Jet12+20 and Jet12+20+30+40.",
                    "plot_class_definition": PLOT_CLASS_DEFINITION,
                    "common_samples": common_samples,
                },
                indent=2,
                sort_keys=True,
            )
        )
        return 0

    if args.stdout_holdout_three_by_three_direct:
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

        common_samples = []
        for common_tag, common_label in BRANCHES:
            common = loaded[common_tag]
            common_idx = common["holdout"]
            common_branches = []
            for model_tag, model_label in BRANCHES:
                model_item = loaded[model_tag]
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
                        model_label,
                        {
                            "validation_mode": f"common_{common_tag}_10pct_holdout_direct_truth_signal_inclusive_jet",
                            "source": str(common["source_dir"]),
                            "model_dir": str(model_item["model_dir"]),
                            "model_registry": str(model_item["model_dir"] / "model_registry.json"),
                            "split_mode": "row",
                            "row_scope": "source_10pct_holdout_direct_model_scoring",
                            "test_fraction_requested": "0.10",
                            "random_seed": str(args.random_seed),
                            "common_validation_sample": common_label,
                            "model_training_sample": model_label,
                            "reported_holdout_rows": str(common["registry"]["reported_holdout_rows"]),
                        },
                        source,
                        y,
                        score,
                        cent,
                        weight,
                        bins,
                    )
                )
            common_samples.append(
                {
                    "validation_sample": common_label,
                    "validation_tag": common_tag,
                    "validation_matrix": str(common["model_dir"] / "training_matrix.npz"),
                    "row_scope": "source_10pct_holdout_direct_model_scoring",
                    "branches": common_branches,
                }
            )
        print(
            json.dumps(
                {
                    "schema": "THE8_BRANCH_A_LADDER_COMPACT_FIXED_HOLDOUT_3X3_DIRECT_SCORE_HISTOGRAMS_V1",
                    "description": "All three BDTs directly scored on all three fixed 10% source holdouts: Jet12+20, Jet12+20+30, and Jet12+20+30+40.",
                    "plot_class_definition": PLOT_CLASS_DEFINITION,
                    "common_samples": common_samples,
                },
                indent=2,
                sort_keys=True,
            )
        )
        return 0

    if args.stdout_holdout_common_grid_from_scorecache:
        common_samples = []
        model_dirs = {
            tag: args.remote_model_base / f"THE8_branchA_{tag}_global_noiso_bdt_mem24_20260527"
            for tag, _ in BRANCHES
        }
        source_dirs = {
            tag: args.remote_source_base / f"THE8_branchA_ladder_{tag}_20260527"
            for tag, _ in BRANCHES
        }
        common_report_tags = {
            ("jet12_20", "jet12_20"): "THE8_branchA_jet12_20_scorecache_fullstat_20260527",
            ("jet12_20", "jet12_20_30"): "THE8_commonJet12_20_modelJet12_20_30_scorecache_fullstat_20260603",
            ("jet12_20", "jet12_20_30_40"): "THE8_commonJet12_20_modelJet12_20_30_40_scorecache_fullstat_20260603",
            ("jet12_20_30_40", "jet12_20"): "THE8_commonJet12_20_30_40_modelJet12_20_scorecache_fullstat_20260603",
            ("jet12_20_30_40", "jet12_20_30"): "THE8_commonJet12_20_30_40_modelJet12_20_30_scorecache_fullstat_20260603",
            ("jet12_20_30_40", "jet12_20_30_40"): "THE8_branchA_jet12_20_30_40_scorecache_fullstat_20260527",
        }
        holdout_by_source = {}
        registry_by_source = {}
        for common_tag, _ in COMMON_VALIDATION_SAMPLES:
            source_matrix_dir = model_dirs[common_tag]
            registry = load_registry(source_matrix_dir / "model_registry.json")
            frame = load_matrix(source_matrix_dir / "training_matrix.npz", registry["features"])
            holdout_by_source[common_tag] = holdout_indices(frame, registry, args.random_seed)
            registry_by_source[common_tag] = registry
        for common_tag, common_label in COMMON_VALIDATION_SAMPLES:
            source_dir = source_dirs[common_tag]
            source_matrix_dir = model_dirs[common_tag]
            registry = registry_by_source[common_tag]
            branches = []
            for model_tag, model_label in BRANCHES:
                model_dir = model_dirs[model_tag]
                report_tag = common_report_tags[(common_tag, model_tag)]
                report_dir = source_dir / "reports" / f"model_validation_condor_{report_tag}"
                branches.append(
                    branch_record_from_score_caches(
                        model_label,
                        {
                            "validation_mode": f"common_{common_tag}_10pct_holdout_from_full_scorecache_truth_signal_inclusive_jet",
                            "source": str(source_dir),
                            "model_dir": str(model_dir),
                            "model_registry": str(model_dir / "model_registry.json"),
                            "report_dir": str(report_dir),
                            "split_mode": "row",
                            "row_scope": "source_10pct_holdout_from_full_scorecache",
                            "test_fraction_requested": "0.10",
                            "random_seed": str(args.random_seed),
                            "reported_holdout_rows": str(registry["reported_holdout_rows"]),
                            "common_validation_sample": common_label,
                            "model_training_sample": model_label,
                        },
                        model_dir,
                        report_dir,
                        bins,
                        matrix_dir=source_matrix_dir,
                        keep_indices=holdout_by_source[common_tag],
                    )
                )
            common_samples.append(
                {
                    "validation_sample": common_label,
                    "validation_tag": common_tag,
                    "validation_matrix": str(source_matrix_dir / "training_matrix.npz"),
                    "row_scope": "source_10pct_holdout_from_full_scorecache",
                    "branches": branches,
                }
            )
        print(
            json.dumps(
                {
                    "schema": "THE8_BRANCH_A_LADDER_COMPACT_FIXED_HOLDOUT_SCORECACHE_HISTOGRAMS_V1",
                    "description": "All three BDTs scored on each fixed 10% source holdout: Jet12+20 and Jet12+20+30+40. Scores are read from the full score-cache reports.",
                    "plot_class_definition": PLOT_CLASS_DEFINITION,
                    "common_samples": common_samples,
                },
                indent=2,
                sort_keys=True,
            )
        )
        return 0

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

    common_outputs = {}
    for common_tag, common_label in COMMON_VALIDATION_SAMPLES:
        common = loaded[common_tag]
        common_idx = common["holdout"]
        common_branches = []
        mode_tag = common_tag.replace("_", "_")
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
                        "validation_mode": f"common_{mode_tag}_10pct_training_holdout_truth_signal_inclusive_jet",
                        "source": str(common["source_dir"]),
                        "model_dir": str(model_item["model_dir"]),
                        "model_registry": str(model_item["model_dir"] / "model_registry.json"),
                        "split_mode": "row",
                        "test_fraction_requested": "0.10",
                        "random_seed": str(args.random_seed),
                        "common_holdout_sample": common_label,
                    },
                    source,
                    y,
                    score,
                    cent,
                    weight,
                    bins,
                )
            )
        common_outputs[f"the8_branch_a_ladder_common_{common_tag}_holdout_score_histograms.json"] = {
            "schema": "THE8_BRANCH_A_LADDER_COMPACT_HOLDOUT_SCORE_HISTOGRAMS_V1",
            "description": f"All three BDTs scored on the same {common_label} 10% row holdout from the Branch A training matrix.",
            "plot_class_definition": PLOT_CLASS_DEFINITION,
            "branches": common_branches,
        }

    outputs = {
        "the8_branch_a_ladder_own_holdout_score_histograms.json": {
            "schema": "THE8_BRANCH_A_LADDER_COMPACT_HOLDOUT_SCORE_HISTOGRAMS_V1",
            "description": "Each BDT scored on its own 10% row holdout from the matching Branch A training matrix.",
            "plot_class_definition": PLOT_CLASS_DEFINITION,
            "branches": own_branches,
        },
        **common_outputs,
    }
    for name, payload in outputs.items():
        path = args.outdir / name
        path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
