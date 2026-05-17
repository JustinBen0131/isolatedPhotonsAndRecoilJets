#!/usr/bin/env python3
"""Plot feature-use diagnostics for exported AuAu BDT+MLP stack artifacts."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_auau_bdt_mlp_stack_roc_overlay as rocutil  # noqa: E402
import train_auau_stacked_bdt_mlp_sweep as stack  # noqa: E402


DEFAULT_STACK_DIR = Path(
    "dataOutput/auauBDTMLPStackFullStat/"
    "stacked_bdt_mlp_targeted_fullstat_ptFine15to35_cent7_full_20260513_123934"
)
DEFAULT_BDT_CACHE = Path(
    "dataOutput/auauMLPDiagnosticPlots/stack_score_separation_fineEt_cent7_20260514/"
    "local_aligned_bdt_score_caches.list"
)
DEFAULT_MLP_CACHE = Path(
    "dataOutput/auauMLPDiagnosticPlots/stack_score_separation_fineEt_cent7_20260514/"
    "local_aligned_mlp_score_caches.list"
)
DEFAULT_OUTDIR = Path("dataOutput/auauMLPDiagnosticPlots/stack_feature_importance_20260515")


DERIVED_FEATURES = {
    "mlp_logit",
    "bdt_is_finite",
    "mlp_is_finite",
    "log_cluster_Et",
    "centrality_scaled",
    "bdt_x_mlp_logit",
    "bdt_x_logEt",
    "mlp_logit_x_logEt",
    "bdt_x_cent",
    "mlp_logit_x_cent",
    "logEt_x_cent",
}


def label_sphenix(fig, y=0.985) -> None:
    fig.text(0.026, y, "sPHENIX", ha="left", va="top", fontsize=15, fontstyle="italic", fontweight="bold")
    fig.text(0.130, y, "Internal", ha="left", va="top", fontsize=15)


def setup_style() -> None:
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "axes.linewidth": 1.1,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def feature_group(name: str) -> str:
    if name == "bdt_score":
        return "BDT score"
    if name in {"mlp_score", "mlp_logit"}:
        return "MLP score"
    if name in {"bdt_is_finite", "mlp_is_finite"}:
        return "finite flags"
    if name in {"cluster_Et", "log_cluster_Et", "centrality", "centrality_scaled"}:
        return "E_T / centrality"
    if name in {"cluster_weta_cogx", "cluster_wphi_cogx", "cluster_weta_over_wphi"}:
        return "full-width shapes"
    if name in {
        "cluster_weta33_cogx",
        "cluster_wphi33_cogx",
        "cluster_weta33_over_wphi33",
        "cluster_weta35_cogx",
        "cluster_wphi53_cogx",
        "cluster_w32",
        "cluster_w52",
        "cluster_w72",
    }:
        return "local/core widths"
    if name.startswith("e") or name.startswith("cluster_et"):
        return "energy sharing"
    if name in {"cluster_Eta", "vertexz"}:
        return "geometry / vertex"
    return "other"


def pretty_feature(name: str) -> str:
    mapping = {
        "bdt_score": "BDT score",
        "mlp_score": "MLP score",
        "mlp_logit": "MLP logit",
        "cluster_Et": r"cluster $E_T$",
        "cluster_weta_cogx": r"$w_\eta$",
        "cluster_wphi_cogx": r"$w_\phi$",
        "cluster_weta33_cogx": r"$w_{\eta,3\times3}$",
        "cluster_wphi33_cogx": r"$w_{\phi,3\times3}$",
        "e11_over_e33": r"$e_{11}/e_{33}$",
        "e32_over_e35": r"$e_{32}/e_{35}$",
    }
    return mapping.get(name, name.replace("cluster_", "").replace("_", " "))


def logistic_importance(model: dict) -> dict[str, float]:
    names = list(model["feature_names"])
    coef = np.abs(np.asarray(model["coef"], dtype="float64"))
    total = float(np.sum(coef))
    if total <= 0.0:
        return {name: 0.0 for name in names}
    return {name: float(val / total) for name, val in zip(names, coef, strict=True)}


def gbm_importance(model: dict) -> dict[str, float]:
    names = list(model["feature_names"])
    counts = np.zeros(len(names), dtype="float64")
    for tree in model.get("trees", []):
        feat = np.asarray(tree.get("feature", []), dtype="int32")
        for idx in feat[feat >= 0]:
            if idx < len(counts):
                counts[idx] += 1.0
    total = float(np.sum(counts))
    if total <= 0.0:
        return {name: 0.0 for name in names}
    return {name: float(val / total) for name, val in zip(names, counts, strict=True)}


def nn_path_importance(model: dict) -> dict[str, float]:
    names = list(model["feature_names"])
    layers = model.get("layers", [])
    if not layers:
        return {name: 0.0 for name in names}
    path = np.abs(np.asarray(layers[0]["weight"], dtype="float64"))
    for layer in layers[1:]:
        path = path @ np.abs(np.asarray(layer["weight"], dtype="float64"))
    vals = np.ravel(path)
    total = float(np.sum(vals))
    if total <= 0.0:
        return {name: 0.0 for name in names}
    return {name: float(val / total) for name, val in zip(names, vals, strict=True)}


def route_importance(artifact: dict) -> dict[str, float]:
    algo = artifact.get("algorithm", "")
    accum: dict[str, list[float]] = defaultdict(list)
    for route in artifact.get("routes", []):
        model = route["model"]
        if algo == "logistic":
            imp = logistic_importance(model)
        elif algo == "gbm":
            imp = gbm_importance(model)
        elif algo == "nn":
            imp = nn_path_importance(model)
        else:
            raise SystemExit(f"Unsupported stack algorithm for importance: {algo}")
        for key, val in imp.items():
            accum[key].append(val)
    return {key: float(np.mean(vals)) for key, vals in accum.items()}


def read_artifacts(stack_dir: Path) -> dict[str, dict]:
    out = {}
    for path in sorted((stack_dir / "artifacts").glob("ptFine15to35_cent7_full_*.json")):
        artifact = json.loads(path.read_text())
        out[artifact["algorithm"]] = artifact
    if not out:
        raise SystemExit(f"No stack artifacts found under {stack_dir / 'artifacts'}")
    return out


def write_rows(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_structural(artifact_imps: dict[str, dict[str, float]], outdir: Path) -> tuple[Path, Path]:
    import matplotlib.pyplot as plt

    alg_labels = {"logistic": "Logistic stack\n|standardized coef|", "gbm": "Tree stack\nsplit-use fraction", "nn": "NN stack\npath-strength"}
    colors = {"BDT score": "#1b4f8a", "MLP score": "#7a3db8", "E_T / centrality": "#3a923a", "full-width shapes": "#b5651d", "local/core widths": "#d95f02", "energy sharing": "#b23a48", "geometry / vertex": "#666666", "finite flags": "#aaaaaa", "other": "#444444"}

    rows = []
    group_rows = []
    for algo, imp in artifact_imps.items():
        for feat, val in sorted(imp.items(), key=lambda item: item[1], reverse=True):
            rows.append({"algorithm": algo, "feature": feat, "group": feature_group(feat), "importance_fraction": val})
        group_sum: dict[str, float] = defaultdict(float)
        for feat, val in imp.items():
            group_sum[feature_group(feat)] += val
        for group, val in sorted(group_sum.items(), key=lambda item: item[1], reverse=True):
            group_rows.append({"algorithm": algo, "group": group, "importance_fraction": val})

    write_rows(outdir / "stack_structural_feature_importance.csv", rows)
    write_rows(outdir / "stack_structural_group_importance.csv", group_rows)

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 6.2))
    for ax, algo in zip(axes, ["logistic", "gbm", "nn"], strict=True):
        top = [r for r in rows if r["algorithm"] == algo][:12]
        top = list(reversed(top))
        vals = [r["importance_fraction"] for r in top]
        labels = [pretty_feature(r["feature"]) for r in top]
        bar_colors = [colors.get(r["group"], "#444444") for r in top]
        ax.barh(np.arange(len(top)), vals, color=bar_colors, edgecolor="white", linewidth=0.8)
        ax.set_yticks(np.arange(len(top)), labels, fontsize=9.5)
        ax.set_title(alg_labels.get(algo, algo), fontsize=13, fontweight="bold")
        ax.set_xlabel("Mean per-route fraction", fontsize=10.5)
        ax.grid(axis="x", color="0.90")
        ax.set_axisbelow(True)
        for i, v in enumerate(vals):
            ax.text(v + max(vals) * 0.015, i, f"{100*v:.1f}%", va="center", fontsize=8.5)
    label_sphenix(fig, y=0.975)
    fig.text(0.54, 0.965, "Feature use inside routed BDT+MLP stackers", ha="center", va="top", fontsize=20, fontweight="bold")
    fig.text(
        0.54,
        0.915,
        "42 local submodels are averaged; NN uses absolute input-to-output path strength, not raw first-layer weights",
        ha="center",
        va="top",
        fontsize=11.5,
        color="0.35",
    )
    fig.tight_layout(rect=[0.02, 0.04, 0.99, 0.88], w_pad=2.1)
    top_path = outdir / "stack_top_feature_importance_by_algorithm.png"
    fig.savefig(top_path, dpi=220)
    plt.close(fig)

    all_groups = ["BDT score", "MLP score", "E_T / centrality", "energy sharing", "local/core widths", "full-width shapes", "geometry / vertex", "finite flags", "other"]
    algos = ["logistic", "gbm", "nn"]
    x = np.arange(len(algos))
    bottom = np.zeros(len(algos), dtype="float64")
    fig, ax = plt.subplots(figsize=(10.8, 6.3))
    group_lookup = {(r["algorithm"], r["group"]): r["importance_fraction"] for r in group_rows}
    for group in all_groups:
        vals = np.array([group_lookup.get((algo, group), 0.0) for algo in algos])
        if np.all(vals == 0.0):
            continue
        ax.bar(x, vals, bottom=bottom, label=group, color=colors.get(group, "#444444"), edgecolor="white", linewidth=0.8)
        bottom += vals
    ax.set_xticks(x, ["Logistic", "Tree / GBM", "Neural net"], fontsize=12)
    ax.set_ylabel("Mean per-route structural fraction", fontsize=12)
    ax.set_ylim(0, 1.0)
    ax.grid(axis="y", color="0.90")
    ax.set_axisbelow(True)
    ax.legend(frameon=False, ncol=3, fontsize=10, loc="upper center", bbox_to_anchor=(0.5, -0.10))
    label_sphenix(fig, y=0.975)
    fig.text(0.54, 0.965, "Grouped information used by each BDT+MLP stacker", ha="center", va="top", fontsize=20, fontweight="bold")
    fig.text(0.54, 0.910, "Structural diagnostic from exported artifacts; route choice already uses fine ET x 7-centrality bins", ha="center", fontsize=11.5, color="0.35")
    fig.tight_layout(rect=[0.04, 0.10, 0.98, 0.86])
    group_path = outdir / "stack_group_importance_by_algorithm.png"
    fig.savefig(group_path, dpi=220)
    plt.close(fig)
    return top_path, group_path


def required_raw_features(artifacts: dict[str, dict]) -> set[str]:
    names = {"is_signal", "cluster_Et", "centrality", "reco_eiso"}
    for artifact in artifacts.values():
        names.update(artifact.get("feature_names", []))
        for route in artifact.get("routes", []):
            names.update(route["model"].get("feature_names", []))
    names.difference_update(DERIVED_FEATURES)
    names.difference_update({"bdt_score", "mlp_score"})
    return names


def load_aligned_frame(args: argparse.Namespace, artifacts: dict[str, dict]) -> dict[str, np.ndarray]:
    mlp_paths = stack.read_manifest(args.mlp_cache)
    bdt_paths = stack.read_manifest(args.bdt_cache)
    if len(mlp_paths) != len(bdt_paths):
        raise SystemExit(f"Cache count mismatch: MLP {len(mlp_paths)} vs BDT {len(bdt_paths)}")
    wanted = required_raw_features(artifacts)
    pieces: dict[str, list[np.ndarray]] = {key: [] for key in sorted(wanted)}
    pieces["mlp_score"] = []
    pieces["bdt_score"] = []
    mlp_key = args.mlp_score
    bdt_key = args.bdt_score
    rows = 0
    for idx, (mlp_path, bdt_path) in enumerate(zip(mlp_paths, bdt_paths, strict=True), 1):
        with np.load(mlp_path, allow_pickle=True) as mlp_cache, np.load(bdt_path, allow_pickle=True) as bdt_cache:
            if idx == 1:
                mlp_key = stack.auto_score_key(mlp_cache, args.mlp_score, "MLP cache")
                bdt_key = stack.auto_score_key(bdt_cache, args.bdt_score, "BDT cache")
            stack.compare_alignment(mlp_cache, bdt_cache, f"shard {idx}")
            n = len(mlp_cache["is_signal"])
            take = n if args.max_rows <= 0 else max(0, min(n, args.max_rows - rows))
            if take <= 0:
                break
            for key in wanted:
                if key in mlp_cache.files:
                    arr = np.asarray(mlp_cache[key])[:take]
                elif key in bdt_cache.files:
                    arr = np.asarray(bdt_cache[key])[:take]
                else:
                    raise SystemExit(f"Missing required feature {key} in aligned caches")
                pieces[key].append(arr)
            pieces["mlp_score"].append(np.asarray(mlp_cache[mlp_key], dtype="float64")[:take])
            pieces["bdt_score"].append(np.asarray(bdt_cache[bdt_key], dtype="float64")[:take])
            rows += take
            if args.max_rows > 0 and rows >= args.max_rows:
                break
    frame = {key: np.concatenate(parts) for key, parts in pieces.items() if parts}
    stack.add_derived_features(frame)
    return frame


def auc_for(frame: dict[str, np.ndarray], score: np.ndarray, pt_range=(15.0, 35.0)) -> float:
    et = np.asarray(frame["cluster_Et"], dtype="float64")
    mask = (et >= pt_range[0]) & (et < pt_range[1]) & np.isfinite(score)
    _, _, auc = rocutil.roc_curve_auc(np.asarray(frame["is_signal"], dtype="int32")[mask], np.asarray(score, dtype="float64")[mask])
    return auc


def copy_frame(frame: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    return {key: np.array(val, copy=True) for key, val in frame.items()}


def permutation_importance(artifacts: dict[str, dict], frame: dict[str, np.ndarray], args: argparse.Namespace, outdir: Path) -> Path:
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(args.random_seed)
    groups = {
        "BDT score": ["bdt_score", "bdt_is_finite"],
        "MLP score": ["mlp_score", "mlp_logit", "mlp_is_finite"],
        "E_T / route": ["cluster_Et", "log_cluster_Et"],
        "centrality route": ["centrality", "centrality_scaled"],
        "full widths": ["cluster_weta_cogx", "cluster_wphi_cogx", "cluster_weta_over_wphi"],
        "local/core widths": [
            "cluster_weta33_cogx",
            "cluster_wphi33_cogx",
            "cluster_weta33_over_wphi33",
            "cluster_weta35_cogx",
            "cluster_wphi53_cogx",
            "cluster_w32",
            "cluster_w52",
            "cluster_w72",
        ],
        "energy sharing": [
            "e11_over_e33",
            "cluster_et1",
            "cluster_et2",
            "cluster_et3",
            "cluster_et4",
            "e32_over_e35",
            "e11_over_e22",
            "e11_over_e13",
            "e11_over_e15",
            "e11_over_e17",
            "e11_over_e31",
            "e11_over_e51",
            "e11_over_e71",
            "e22_over_e33",
            "e22_over_e35",
            "e22_over_e37",
            "e22_over_e53",
        ],
    }
    rows = []
    for algo in ["logistic", "gbm", "nn"]:
        artifact = artifacts.get(algo)
        if artifact is None:
            continue
        base_score = rocutil.route_score(artifact, copy_frame(frame))
        base_auc = auc_for(frame, base_score)
        for group, features in groups.items():
            drops = []
            usable = [feature for feature in features if feature in frame]
            if not usable:
                continue
            for _ in range(args.permutation_repeats):
                trial = copy_frame(frame)
                perm = rng.permutation(len(frame["is_signal"]))
                for feature in usable:
                    trial[feature] = np.asarray(trial[feature])[perm]
                stack.add_derived_features(trial)
                score = rocutil.route_score(artifact, trial)
                drops.append(base_auc - auc_for(trial, score))
            rows.append(
                {
                    "algorithm": algo,
                    "group": group,
                    "base_auc": base_auc,
                    "auc_drop_mean": float(np.mean(drops)),
                    "auc_drop_std": float(np.std(drops)),
                    "features": "+".join(usable),
                }
            )
    write_rows(outdir / "stack_permutation_group_importance.csv", rows)

    algos = ["logistic", "gbm", "nn"]
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.8), sharex=True)
    for ax, algo in zip(axes, algos, strict=True):
        sub = [r for r in rows if r["algorithm"] == algo]
        sub.sort(key=lambda r: r["auc_drop_mean"], reverse=True)
        sub = list(reversed(sub[:8]))
        vals = [r["auc_drop_mean"] for r in sub]
        errs = [r["auc_drop_std"] for r in sub]
        labels = [r["group"] for r in sub]
        ax.barh(np.arange(len(sub)), vals, xerr=errs, color="#1f77b4" if algo == "nn" else "#6a6a6a", alpha=0.88)
        ax.set_yticks(np.arange(len(sub)), labels, fontsize=10)
        base_auc = sub[0]["base_auc"] if sub else math.nan
        title = {"logistic": "Logistic", "gbm": "Tree / GBM", "nn": "Neural net"}[algo]
        ax.set_title(f"{title}\nbase AUC={base_auc:.3f}", fontsize=13, fontweight="bold")
        ax.grid(axis="x", color="0.90")
        ax.set_axisbelow(True)
        for i, val in enumerate(vals):
            ax.text(val + max(max(vals), 0.005) * 0.03, i, f"{val:.3f}", va="center", fontsize=9)
    axes[0].set_xlabel("AUC drop after permutation", fontsize=11)
    axes[1].set_xlabel("AUC drop after permutation", fontsize=11)
    axes[2].set_xlabel("AUC drop after permutation", fontsize=11)
    label_sphenix(fig, y=0.975)
    fig.text(0.54, 0.965, "Data-driven group importance for BDT+MLP stackers", ha="center", va="top", fontsize=20, fontweight="bold")
    fig.text(
        0.54,
        0.910,
        f"Aligned validation cache, first {len(frame['is_signal']):,} rows; larger AUC drop means the stacker relied more on that information",
        ha="center",
        fontsize=11.5,
        color="0.35",
    )
    fig.tight_layout(rect=[0.02, 0.05, 0.99, 0.86], w_pad=2.0)
    path = outdir / "stack_permutation_group_importance.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stack-dir", type=Path, default=DEFAULT_STACK_DIR)
    ap.add_argument("--bdt-cache", type=Path, default=DEFAULT_BDT_CACHE)
    ap.add_argument("--mlp-cache", type=Path, default=DEFAULT_MLP_CACHE)
    ap.add_argument("--mlp-score", default=stack.DEFAULT_MLP_SCORE)
    ap.add_argument("--bdt-score", default=stack.DEFAULT_BDT_SCORE)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    ap.add_argument("--max-rows", type=int, default=300000)
    ap.add_argument("--permutation-repeats", type=int, default=3)
    ap.add_argument("--random-seed", type=int, default=753)
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    setup_style()
    artifacts = read_artifacts(args.stack_dir)
    artifact_imps = {algo: route_importance(artifact) for algo, artifact in artifacts.items()}
    structural_paths = plot_structural(artifact_imps, args.outdir)
    frame = load_aligned_frame(args, artifacts)
    permutation_path = permutation_importance(artifacts, frame, args, args.outdir)
    manifest = {
        "schema": "RJ_AUAU_STACK_FEATURE_IMPORTANCE_V1",
        "stack_dir": str(args.stack_dir),
        "bdt_cache": str(args.bdt_cache),
        "mlp_cache": str(args.mlp_cache),
        "max_rows": args.max_rows,
        "rows_loaded": int(len(frame["is_signal"])),
        "plots": [str(p) for p in (*structural_paths, permutation_path)],
        "caveat": "NN structural path-strength is model-internal. Permutation AUC drop is the data-driven check.",
    }
    (args.outdir / "stack_feature_importance_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    for path in (*structural_paths, permutation_path):
        print(path)


if __name__ == "__main__":
    main()
