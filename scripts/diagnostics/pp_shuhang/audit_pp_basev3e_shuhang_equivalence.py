#!/usr/bin/env python3
"""Audit pp baseV3E score equivalence against Shuhang/PPG12 TMVA models.

This diagnostic intentionally uses existing extracted ROOT tables and existing
frozen models only. It does not train, submit Condor, or use Shuhang ROOT files
as training inputs.
"""

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
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from train_auau_photon_bdt import add_derived_features, expand_required_columns  # noqa: E402


DEFAULT_TREE = "AuAuPhotonIDTrainingTree"
DEFAULT_PRODUCT = "ppg12_base_v3E_bdt_noIso"
BASEV3E_FEATURES = [
    "cluster_Et",
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "vertexz",
    "cluster_Eta",
    "e11_over_e33",
    "cluster_et1",
    "cluster_et2",
    "cluster_et3",
    "cluster_et4",
    "e32_over_e35",
]
AUDIT_COLUMNS = [
    "cluster_Et",
    "cluster_Eta",
    "npb_score",
    "cluster_prob",
    "e11_over_e33",
    "e32_over_e35",
    "cluster_et1",
    "event_weight",
    "tight_bdt_score",
]
SAMPLES = [
    "run28_photonjet5",
    "run28_photonjet10",
    "run28_photonjet20",
    "run28_jet8",
    "run28_jet12",
    "run28_jet20",
    "run28_jet30",
    "run28_jet40",
    "photonjet5",
    "photonjet10",
    "photonjet20",
    "jet8",
    "jet12",
    "jet20",
    "jet30",
    "jet40",
]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--signal", nargs="+", type=Path, required=True)
    ap.add_argument("--inclusive", nargs="+", type=Path, required=True)
    ap.add_argument("--registry", type=Path, required=True)
    ap.add_argument("--product", default=DEFAULT_PRODUCT)
    ap.add_argument("--tree", default=DEFAULT_TREE)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--pt-range", default="22:28")
    ap.add_argument("--eta-max", type=float, default=0.7)
    ap.add_argument("--npb-min", type=float, default=0.5)
    ap.add_argument("--bins", type=int, default=50)
    ap.add_argument("--step-size", default="100 MB")
    ap.add_argument("--weight-branch", default="event_weight")
    ap.add_argument("--max-files-per-sample", type=int, default=0)
    ap.add_argument("--max-scatter-rows-per-sample", type=int, default=25000)
    ap.add_argument("--random-seed", type=int, default=42)
    ap.add_argument("--shuhang-split-model", type=Path, required=True)
    ap.add_argument("--shuhang-nosplit-model", type=Path, required=True)
    ap.add_argument(
        "--sample-weight-scale",
        action="append",
        default=[],
        metavar="SAMPLE=SCALE",
        help="Optional multiplicative correction applied to event weights for matching path samples.",
    )
    return ap.parse_args()


def parse_range(text: str) -> tuple[float, float]:
    lo_s, hi_s = str(text).split(":", 1)
    lo = float(lo_s)
    hi = float(hi_s)
    if hi <= lo:
        raise SystemExit(f"Bad --pt-range {text!r}")
    return lo, hi


def expand_paths(items: Iterable[Path]) -> list[Path]:
    paths: list[Path] = []
    for item in items:
        text = str(item)
        if text.startswith("@"):
            manifest = Path(text[1:])
            paths.extend(Path(line.strip()) for line in manifest.read_text().splitlines() if line.strip())
        elif item.is_file() and item.suffix in {".list", ".txt"}:
            paths.extend(Path(line.strip()) for line in item.read_text().splitlines() if line.strip())
        elif item.is_dir():
            paths.extend(sorted(item.rglob("*.root")))
        else:
            paths.append(item)
    unique: list[Path] = []
    seen: set[str] = set()
    for path in paths:
        key = str(path)
        if key not in seen:
            unique.append(path)
            seen.add(key)
    return unique


def sample_name(path: Path) -> str:
    text = str(path)
    for name in SAMPLES:
        if name in text:
            return name.replace("run28_", "")
    return "unknown"


def limit_by_sample(paths: list[Path], max_files_per_sample: int) -> list[Path]:
    if max_files_per_sample <= 0:
        return paths
    counts: dict[str, int] = defaultdict(int)
    limited: list[Path] = []
    for path in paths:
        sample = sample_name(path)
        if counts[sample] >= max_files_per_sample:
            continue
        limited.append(path)
        counts[sample] += 1
    return limited


def parse_sample_weight_scales(items: list[str]) -> dict[str, float]:
    out: dict[str, float] = {}
    for item in items:
        if "=" not in item:
            raise SystemExit(f"Bad --sample-weight-scale {item!r}; expected SAMPLE=SCALE")
        sample, value = item.split("=", 1)
        scale = float(value)
        if not math.isfinite(scale) or scale <= 0.0:
            raise SystemExit(f"Bad --sample-weight-scale {item!r}; scale must be positive finite")
        out[sample.strip()] = scale
    return out


def weight_scale_for_path(path: Path, scales: dict[str, float]) -> float:
    text = str(path)
    for sample, scale in scales.items():
        if sample in text or sample.replace("run28_", "") in text:
            return float(scale)
    return 1.0


def load_model(registry_path: Path, product: str):
    import xgboost as xgb

    registry = json.loads(registry_path.read_text())
    chosen = None
    for model in registry.get("models", []):
        names = {str(model.get(key, "")) for key in ("product", "model_id", "name")}
        if product in names or any(product in name for name in names):
            chosen = model
            break
    if chosen is None:
        known = [str(m.get("product") or m.get("model_id")) for m in registry.get("models", [])]
        raise SystemExit(f"Could not find product {product!r}; known={known}")
    features = [str(name) for name in chosen.get("features", [])]
    if not features:
        raise SystemExit(f"Model {product!r} has no feature list")
    xgb_path = chosen.get("output_xgb_json") or chosen.get("report", {}).get("output_xgb_json")
    if not xgb_path:
        raise SystemExit(f"Model {product!r} has no output_xgb_json")
    booster = xgb.Booster()
    booster.load_model(str(xgb_path))
    return booster, features, chosen


def required_columns(our_features: list[str]) -> list[str]:
    cols = set(expand_required_columns(our_features))
    cols.update(BASEV3E_FEATURES)
    cols.update(AUDIT_COLUMNS)
    return sorted(cols)


def column(frame, name: str, default: float = np.nan) -> np.ndarray:
    if name in frame:
        return frame[name].to_numpy(dtype="float64")
    return np.full(len(frame), float(default), dtype="float64")


def event_weights(frame, branch: str, path: Path, scales: dict[str, float]) -> np.ndarray:
    if branch in frame:
        weight = frame[branch].to_numpy(dtype="float64")
        weight = np.where(np.isfinite(weight) & (weight > 0.0), weight, 0.0)
        if float(np.sum(weight)) <= 0.0:
            weight = np.ones(len(frame), dtype="float64")
    else:
        weight = np.ones(len(frame), dtype="float64")
    return weight * weight_scale_for_path(path, scales)


def finite_features(frame, features: list[str]) -> np.ndarray:
    mask = np.ones(len(frame), dtype=bool)
    for feature in features:
        if feature not in frame:
            return np.zeros(len(frame), dtype=bool)
        values = frame[feature].to_numpy(dtype="float64")
        mask &= np.isfinite(values)
    return mask


def fig19_common_mask(frame) -> np.ndarray:
    prob = column(frame, "cluster_prob", 0.5)
    e11e33 = column(frame, "e11_over_e33")
    e32e35 = column(frame, "e32_over_e35")
    et1 = column(frame, "cluster_et1")
    return (
        np.isfinite(prob)
        & (prob > 0.0)
        & (prob < 1.0)
        & np.isfinite(e11e33)
        & (e11e33 > 0.0)
        & (e11e33 < 0.98)
        & np.isfinite(e32e35)
        & (e32e35 > 0.8)
        & (e32e35 < 1.0)
        & np.isfinite(et1)
        & (et1 > 0.6)
        & (et1 < 1.0)
    )


def npb_masks(frame, npb_min: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    npb = column(frame, "npb_score")
    real = np.isfinite(npb) & (npb >= 0.0) & (npb <= 1.0)
    return real, real & (npb > npb_min), real & (npb <= npb_min)


def update_hist(hists: dict[str, np.ndarray], key: str, score: np.ndarray, weight: np.ndarray, edges: np.ndarray):
    if key not in hists:
        hists[key] = np.zeros(len(edges) - 1, dtype="float64")
    if len(score):
        hists[key] += np.histogram(score, bins=edges, weights=weight)[0]


def update_summary(summary: dict, key: tuple[str, str, str], rows: int, weight: float, low: int, high: int):
    item = summary.setdefault(
        key,
        {
            "class": key[0],
            "sample": key[1],
            "selection": key[2],
            "rows": 0,
            "weighted_rows": 0.0,
            "score_lt_0p1": 0,
            "score_gt_0p8": 0,
        },
    )
    item["rows"] += int(rows)
    item["weighted_rows"] += float(weight)
    item["score_lt_0p1"] += int(low)
    item["score_gt_0p8"] += int(high)


def choose_scatter_rows(scoreable: np.ndarray, sample_keep: int, rng: np.random.Generator) -> np.ndarray:
    idx = np.flatnonzero(scoreable)
    if sample_keep <= 0 or len(idx) <= sample_keep:
        return idx
    return rng.choice(idx, size=sample_keep, replace=False)


def tmva_scores(model, matrix: np.ndarray) -> np.ndarray:
    import ROOT

    vec = ROOT.std.vector("float")()
    out = np.empty(matrix.shape[0], dtype="float64")
    for i, row in enumerate(matrix.astype("float32", copy=False)):
        vec.clear()
        for value in row:
            vec.push_back(float(value))
        result = model.Compute(vec)
        out[i] = float(result[0]) if len(result) else np.nan
    return out


def load_tmva(path: Path):
    import ROOT

    if not Path(path).is_file():
        raise SystemExit(f"Missing TMVA model: {path}")
    return ROOT.TMVA.Experimental.RBDT("myBDT", str(path))


def auc_or_nan(y_true: np.ndarray, score: np.ndarray) -> float:
    from sklearn.metrics import roc_auc_score

    mask = np.isfinite(score)
    if len(np.unique(y_true[mask])) < 2:
        return float("nan")
    return float(roc_auc_score(y_true[mask], score[mask]))


def corr_or_nan(a: np.ndarray, b: np.ndarray) -> float:
    mask = np.isfinite(a) & np.isfinite(b)
    if int(np.sum(mask)) < 3:
        return float("nan")
    return float(np.corrcoef(a[mask], b[mask])[0, 1])


def spearman_or_nan(a: np.ndarray, b: np.ndarray) -> float:
    mask = np.isfinite(a) & np.isfinite(b)
    if int(np.sum(mask)) < 3:
        return float("nan")
    aa = a[mask]
    bb = b[mask]
    ar = np.empty_like(aa, dtype="float64")
    br = np.empty_like(bb, dtype="float64")
    ar[np.argsort(aa)] = np.arange(len(aa), dtype="float64")
    br[np.argsort(bb)] = np.arange(len(bb), dtype="float64")
    return float(np.corrcoef(ar, br)[0, 1])


def unit(hist: np.ndarray) -> np.ndarray:
    total = float(np.sum(hist))
    if not math.isfinite(total) or total <= 0.0:
        return hist
    return hist / total


def draw_outputs(outdir: Path, edges: np.ndarray, hists: dict[str, np.ndarray], scatter, metrics: dict):
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve

    centers = 0.5 * (edges[:-1] + edges[1:])

    fig, axes = plt.subplots(1, 2, figsize=(13.333, 7.5), dpi=180)
    ax = axes[0]
    for sel, color in [
        ("inclusive/all_no_pre_no_npb", "#1f5fff"),
        ("inclusive/fig19_common_no_npb", "#3aa76d"),
        ("inclusive/npb_pass", "#7b61ff"),
        ("inclusive/npb_fail", "#e45757"),
    ]:
        if sel in hists:
            ax.step(centers, unit(hists[sel]), where="mid", linewidth=2, label=sel.split("/", 1)[1], color=color)
    ax.set_title("Inclusive score shape by selection")
    ax.set_xlabel("This-analysis BDT score")
    ax.set_ylabel("unit-normalized counts")
    ax.set_xlim(0, 1)
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.2)

    ax = axes[1]
    for key, hist in sorted(hists.items()):
        if key.startswith("inclusive_sample/"):
            ax.step(centers, unit(hist), where="mid", linewidth=1.8, label=key.split("/", 1)[1])
    ax.set_title("Inclusive score shape by jet sample")
    ax.set_xlabel("This-analysis BDT score")
    ax.set_xlim(0, 1)
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.2)
    fig.suptitle("pp baseV3E inclusive selection contract audit, 22<E_T<28 GeV", fontsize=16)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(outdir / "inclusive_selection_breakdown_22_28.png")
    plt.close(fig)

    y = scatter["label"]
    fig, ax = plt.subplots(figsize=(7.4, 6.2), dpi=180)
    for name, color in [("our", "#111111"), ("shuhang_split", "#1f5fff"), ("shuhang_nosplit", "#d62728")]:
        fpr, tpr, _ = roc_curve(y, scatter[name])
        ax.plot(fpr, tpr, linewidth=2.2, label=f"{name} AUC={metrics[name + '_auc']:.4f}", color=color)
    ax.plot([0, 1], [0, 1], linestyle="--", linewidth=1, color="0.5")
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title("Same-row ROC comparison, 22<E_T<28 GeV")
    ax.legend(frameon=False, loc="lower right")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(outdir / "score_equivalence_roc_22_28.png")
    plt.close(fig)

    for name, title in [("shuhang_split", "Shuhang split TMVA"), ("shuhang_nosplit", "Shuhang nosplit TMVA")]:
        fig, ax = plt.subplots(figsize=(7.0, 6.2), dpi=180)
        hb = ax.hexbin(scatter["our"], scatter[name], gridsize=75, bins="log", cmap="viridis", mincnt=1)
        ax.plot([0, 1], [0, 1], color="white", linewidth=1.0, linestyle="--", alpha=0.85)
        ax.set_xlabel("This-analysis frozen XGBoost score")
        ax.set_ylabel(f"{title} score")
        ax.set_title(f"Same-row score correlation: {title}")
        label = (
            f"Pearson={metrics['our_vs_' + name + '_pearson']:.4f}\n"
            f"rank corr={metrics['our_vs_' + name + '_spearman']:.4f}"
        )
        ax.text(0.04, 0.95, label, transform=ax.transAxes, va="top", ha="left", fontsize=11,
                bbox=dict(facecolor="white", alpha=0.85, edgecolor="0.8"))
        fig.colorbar(hb, ax=ax, label="log10(rows)")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.grid(alpha=0.15)
        fig.tight_layout()
        suffix = "nosplit" if name.endswith("nosplit") else "split"
        fig.savefig(outdir / f"our_vs_shuhang_tmva_score_scatter_{suffix}.png")
        plt.close(fig)


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.random_seed)
    signal_paths = limit_by_sample(expand_paths(args.signal), args.max_files_per_sample)
    inclusive_paths = limit_by_sample(expand_paths(args.inclusive), args.max_files_per_sample)
    if not signal_paths or not inclusive_paths:
        raise SystemExit("Missing signal or inclusive inputs")

    import pandas as pd
    import uproot
    import xgboost as xgb

    booster, our_features, model_meta = load_model(args.registry, args.product)
    required = required_columns(our_features)
    scales = parse_sample_weight_scales(args.sample_weight_scale)
    pt_lo, pt_hi = parse_range(args.pt_range)
    edges = np.linspace(0.0, 1.0, args.bins + 1)
    hists: dict[str, np.ndarray] = {}
    summary: dict[tuple[str, str, str], dict] = {}
    file_counts: dict[str, int] = defaultdict(int)
    scatter_parts: dict[str, list[np.ndarray]] = defaultdict(list)
    tmva_split = load_tmva(args.shuhang_split_model)
    tmva_nosplit = load_tmva(args.shuhang_nosplit_model)

    per_sample_scatter_seen: dict[tuple[str, str], int] = defaultdict(int)
    paths_with_class = [(p, "signal") for p in signal_paths] + [(p, "inclusive") for p in inclusive_paths]
    for path, class_name in paths_with_class:
        sample = sample_name(path)
        file_counts[f"{class_name}/{sample}"] += 1
        try:
            with uproot.open(path) as root_file:
                if args.tree not in root_file:
                    continue
                tree = root_file[args.tree]
                keys = set(tree.keys())
                missing = [name for name in expand_required_columns(our_features) if name not in keys]
                missing += [name for name in BASEV3E_FEATURES if name not in keys]
                if missing:
                    raise SystemExit(f"{path} missing required scorer columns: {sorted(set(missing))}")
                have = [name for name in required if name in keys]
                for arrays in tree.iterate(have, library="np", step_size=args.step_size):
                    if not arrays:
                        continue
                    n = len(next(iter(arrays.values())))
                    if n == 0:
                        continue
                    frame = pd.DataFrame({name: arrays[name] for name in have})
                    add_derived_features(frame)
                    et = column(frame, "cluster_Et")
                    eta = column(frame, "cluster_Eta")
                    pt_eta = (
                        np.isfinite(et)
                        & (et >= pt_lo)
                        & (et < pt_hi)
                        & np.isfinite(eta)
                        & (np.abs(eta) < args.eta_max)
                    )
                    scoreable = pt_eta & finite_features(frame, our_features) & finite_features(frame, BASEV3E_FEATURES)
                    if not np.any(scoreable):
                        continue
                    our_x = np.column_stack(
                        [frame[feature].to_numpy(dtype="float32")[scoreable] for feature in our_features]
                    )
                    our_score_all = booster.predict(xgb.DMatrix(our_x)).astype("float64")
                    weight_all = event_weights(frame, args.weight_branch, path, scales)[scoreable]
                    common = fig19_common_mask(frame)[scoreable]
                    npb_real, npb_pass, npb_fail = npb_masks(frame, args.npb_min)
                    npb_real = npb_real[scoreable]
                    npb_pass = npb_pass[scoreable]
                    npb_fail = npb_fail[scoreable]

                    selections = {
                        "all_no_pre_no_npb": np.ones(len(our_score_all), dtype=bool),
                        "fig19_common_no_npb": common,
                        "npb_pass": npb_pass,
                        "npb_fail": npb_fail,
                        "fig19_common_npb_pass": common & npb_pass,
                        "fig19_common_npb_fail": common & npb_fail,
                        "npb_missing_or_sentinel": ~npb_real,
                    }
                    for sel_name, sel_mask in selections.items():
                        if not np.any(sel_mask):
                            continue
                        score = our_score_all[sel_mask]
                        weight = weight_all[sel_mask]
                        key = (class_name, sample, sel_name)
                        update_summary(
                            summary,
                            key,
                            int(np.sum(sel_mask)),
                            float(np.sum(weight)),
                            int(np.sum(score < 0.1)),
                            int(np.sum(score > 0.8)),
                        )
                        update_hist(hists, f"{class_name}/{sel_name}", score, weight, edges)
                        if class_name == "inclusive" and sel_name == "all_no_pre_no_npb":
                            update_hist(hists, f"inclusive_sample/{sample}", score, weight, edges)

                    scatter_budget_key = (class_name, sample)
                    remaining = args.max_scatter_rows_per_sample - per_sample_scatter_seen[scatter_budget_key]
                    if remaining > 0:
                        local_scoreable_idx = choose_scatter_rows(
                            np.ones(len(our_score_all), dtype=bool),
                            min(remaining, len(our_score_all)),
                            rng,
                        )
                        tmva_x = np.column_stack(
                            [
                                frame[feature].to_numpy(dtype="float32")[scoreable][local_scoreable_idx]
                                for feature in BASEV3E_FEATURES
                            ]
                        )
                        split_score = tmva_scores(tmva_split, tmva_x)
                        nosplit_score = tmva_scores(tmva_nosplit, tmva_x)
                        scatter_parts["our"].append(our_score_all[local_scoreable_idx])
                        scatter_parts["shuhang_split"].append(split_score)
                        scatter_parts["shuhang_nosplit"].append(nosplit_score)
                        scatter_parts["label"].append(
                            np.full(len(local_scoreable_idx), 1 if class_name == "signal" else 0, dtype="int8")
                        )
                        scatter_parts["class"].append(
                            np.array([class_name] * len(local_scoreable_idx), dtype=object)
                        )
                        scatter_parts["sample"].append(np.array([sample] * len(local_scoreable_idx), dtype=object))
                        scatter_parts["cluster_Et"].append(column(frame, "cluster_Et")[scoreable][local_scoreable_idx])
                        scatter_parts["npb_score"].append(column(frame, "npb_score")[scoreable][local_scoreable_idx])
                        per_sample_scatter_seen[scatter_budget_key] += len(local_scoreable_idx)
        except OSError as exc:
            print(f"[warn] skipping unreadable file {path}: {exc}", file=sys.stderr)

    if not scatter_parts["our"]:
        raise SystemExit("No same-row scatter rows collected")

    scatter = {key: np.concatenate(parts) for key, parts in scatter_parts.items()}
    metrics = {
        "our_auc": auc_or_nan(scatter["label"], scatter["our"]),
        "shuhang_split_auc": auc_or_nan(scatter["label"], scatter["shuhang_split"]),
        "shuhang_nosplit_auc": auc_or_nan(scatter["label"], scatter["shuhang_nosplit"]),
        "our_vs_shuhang_split_pearson": corr_or_nan(scatter["our"], scatter["shuhang_split"]),
        "our_vs_shuhang_split_spearman": spearman_or_nan(scatter["our"], scatter["shuhang_split"]),
        "our_vs_shuhang_nosplit_pearson": corr_or_nan(scatter["our"], scatter["shuhang_nosplit"]),
        "our_vs_shuhang_nosplit_spearman": spearman_or_nan(scatter["our"], scatter["shuhang_nosplit"]),
        "scatter_rows": int(len(scatter["our"])),
    }
    draw_outputs(args.outdir, edges, hists, scatter, metrics)

    csv_path = args.outdir / "score_contract_audit_22_28.csv"
    with csv_path.open("w", newline="") as handle:
        fieldnames = [
            "class",
            "sample",
            "selection",
            "rows",
            "weighted_rows",
            "score_lt_0p1",
            "score_gt_0p8",
            "frac_score_lt_0p1",
            "frac_score_gt_0p8",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for item in sorted(summary.values(), key=lambda x: (x["class"], x["sample"], x["selection"])):
            rows = max(int(item["rows"]), 1)
            writer.writerow(
                {
                    **item,
                    "frac_score_lt_0p1": item["score_lt_0p1"] / rows,
                    "frac_score_gt_0p8": item["score_gt_0p8"] / rows,
                }
            )

    scatter_csv = args.outdir / "score_equivalence_sample_22_28.csv"
    with scatter_csv.open("w", newline="") as handle:
        fieldnames = [
            "class",
            "sample",
            "cluster_Et",
            "npb_score",
            "our_score",
            "shuhang_split_score",
            "shuhang_nosplit_score",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for i in range(len(scatter["our"])):
            writer.writerow(
                {
                    "class": scatter["class"][i],
                    "sample": scatter["sample"][i],
                    "cluster_Et": float(scatter["cluster_Et"][i]),
                    "npb_score": float(scatter["npb_score"][i]) if np.isfinite(scatter["npb_score"][i]) else "",
                    "our_score": float(scatter["our"][i]),
                    "shuhang_split_score": float(scatter["shuhang_split"][i]),
                    "shuhang_nosplit_score": float(scatter["shuhang_nosplit"][i]),
                }
            )

    summary_json = {
        "schema": "PP_BASEV3E_SHUHANG_MODEL_EQUIVALENCE_AUDIT_V1",
        "pt_range": args.pt_range,
        "eta_max": args.eta_max,
        "product": args.product,
        "registry": str(args.registry),
        "our_model_features": our_features,
        "shuhang_base_v3E_features": BASEV3E_FEATURES,
        "shuhang_split_model": str(args.shuhang_split_model),
        "shuhang_nosplit_model": str(args.shuhang_nosplit_model),
        "file_counts": dict(file_counts),
        "metrics": metrics,
        "outputs": {
            "contract_csv": str(csv_path),
            "scatter_csv": str(scatter_csv),
            "selection_breakdown_png": str(args.outdir / "inclusive_selection_breakdown_22_28.png"),
            "roc_png": str(args.outdir / "score_equivalence_roc_22_28.png"),
            "split_scatter_png": str(args.outdir / "our_vs_shuhang_tmva_score_scatter_split.png"),
            "nosplit_scatter_png": str(args.outdir / "our_vs_shuhang_tmva_score_scatter_nosplit.png"),
        },
        "audit_rows": list(sorted(summary.values(), key=lambda x: (x["class"], x["sample"], x["selection"]))),
        "model_meta": {
            "product": model_meta.get("product"),
            "model_id": model_meta.get("model_id"),
            "output_xgb_json": model_meta.get("output_xgb_json"),
        },
    }
    (args.outdir / "score_contract_audit_22_28.json").write_text(json.dumps(summary_json, indent=2) + "\n")
    print(json.dumps({"ok": True, "outdir": str(args.outdir), "metrics": metrics}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
