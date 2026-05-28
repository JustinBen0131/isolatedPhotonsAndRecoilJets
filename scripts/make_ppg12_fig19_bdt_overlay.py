#!/usr/bin/env python3
"""Make a PPG12 Fig. 19-style pp BDT score overlay.

This diagnostic intentionally differs from the pp table-validation score plots:
PPG12 Fig. 19 compares unit-normalized score distributions for signal MC and
full inclusive MC candidates after photon-ID/NPB-style cuts.  It is not a
closure plot over only signal-labeled rows and background-labeled training rows.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Iterable

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from train_auau_photon_bdt import add_derived_features, expand_required_columns  # noqa: E402


DEFAULT_TREE = "AuAuPhotonIDTrainingTree"
DEFAULT_PRODUCT = "ppg12_base_v3E_bdt_noIso"


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--signal", nargs="+", type=Path, required=True, help="Signal ROOT files or list files.")
    ap.add_argument("--inclusive", nargs="+", type=Path, required=True, help="Inclusive MC ROOT files or list files.")
    ap.add_argument("--registry", type=Path, required=True, help="BDT model_registry.json.")
    ap.add_argument("--product", default=DEFAULT_PRODUCT, help="Product/model id to plot.")
    ap.add_argument("--tree", default=DEFAULT_TREE)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--pt-range", default="18:22")
    ap.add_argument("--eta-max", type=float, default=0.7)
    ap.add_argument("--npb-min", type=float, default=0.5)
    ap.add_argument(
        "--npb-mode",
        choices=["auto", "apply", "skip"],
        default="auto",
        help="auto applies NPB only when real scores are present; -2 sentinel chunks are treated as no-NPB.",
    )
    ap.add_argument("--bins", type=int, default=50)
    ap.add_argument("--step-size", default="100 MB")
    ap.add_argument("--weight-branch", default="event_weight")
    ap.add_argument(
        "--sample-weight-scale",
        action="append",
        default=[],
        metavar="SAMPLE=SCALE",
        help="Optional multiplicative correction applied to the weight branch for paths containing SAMPLE.",
    )
    ap.add_argument("--reference-branch", default="tight_bdt_score")
    ap.add_argument("--shuhang-signal-root", type=Path, default=None)
    ap.add_argument("--shuhang-inclusive-root", type=Path, default=None)
    ap.add_argument("--shuhang-hist", default="h2d_bdt_eta0_pt2_cut1")
    ap.add_argument(
        "--preselection",
        choices=["fig19", "none"],
        default="fig19",
        help=(
            "fig19 applies the historical shower-shape preselection used by this overlay "
            "(e11/e33, e32/e35, et1). none keeps only pt/eta, finite score inputs, "
            "and the requested NPB mode, matching Fig. 13-style no-preselection checks."
        ),
    )
    ap.add_argument(
        "--require-weta-cut",
        action="store_true",
        help="Apply the older reference weta<2 cut. Current May-21 PPG12 IAN preselection leaves this off.",
    )
    ap.add_argument(
        "--require-unfiltered-inclusive",
        action="store_true",
        help=(
            "Fail if the inclusive input appears truth-background filtered. "
            "Use for Shuhang/Fig13-style inclusive-MC overlays, where the blue curve is all inclusive MC."
        ),
    )
    ap.add_argument("--max-files", type=int, default=0, help="Optional smoke-test file cap per sample class.")
    return ap.parse_args()


def parse_sample_weight_scales(items: list[str]) -> dict[str, float]:
    scales: dict[str, float] = {}
    for item in items:
        if "=" not in item:
            raise SystemExit(f"Bad --sample-weight-scale {item!r}; expected SAMPLE=SCALE")
        sample, value = item.split("=", 1)
        sample = sample.strip()
        if not sample:
            raise SystemExit(f"Bad --sample-weight-scale {item!r}; empty sample")
        scale = float(value)
        if not math.isfinite(scale) or scale <= 0.0:
            raise SystemExit(f"Bad --sample-weight-scale {item!r}; scale must be positive finite")
        scales[sample] = scale
    return scales


def weight_scale_for_path(path: Path, scales: dict[str, float]) -> float:
    text = str(path)
    for sample, scale in scales.items():
        if sample in text or sample.replace("run28_", "") in text:
            return float(scale)
    return 1.0


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
            matches = sorted(Path().glob(text)) if any(ch in text for ch in "*?[]") else [item]
            paths.extend(matches)
    unique: list[Path] = []
    seen: set[str] = set()
    for path in paths:
        key = str(path)
        if key not in seen:
            unique.append(path)
            seen.add(key)
    return unique


def parse_range(text: str) -> tuple[float, float]:
    lo_s, hi_s = str(text).split(":", 1)
    lo = float(lo_s)
    hi = float(hi_s)
    if hi <= lo:
        raise SystemExit(f"Bad range: {text}")
    return lo, hi


def load_model(registry_path: Path, product: str):
    try:
        import xgboost as xgb
    except ImportError as exc:
        raise SystemExit("This plot requires xgboost in the active Python environment") from exc

    registry = json.loads(registry_path.read_text())
    models = registry.get("models", [])
    if not isinstance(models, list) or not models:
        raise SystemExit(f"Registry has no models: {registry_path}")
    chosen = None
    for model in models:
        names = {str(model.get(key, "")) for key in ("product", "model_id", "name")}
        if product in names or any(product in name for name in names):
            chosen = model
            break
    if chosen is None:
        known = [str(m.get("product") or m.get("model_id")) for m in models]
        raise SystemExit(f"Could not find product {product!r}. Known products: {known}")
    features = [str(name) for name in chosen.get("features", [])]
    if not features:
        raise SystemExit(f"Model {product} has no feature list in registry")
    xgb_path = chosen.get("output_xgb_json") or chosen.get("report", {}).get("output_xgb_json")
    if not xgb_path:
        raise SystemExit(f"Model {product} has no output_xgb_json")
    xgb_path = Path(xgb_path)
    if not xgb_path.is_file():
        raise SystemExit(f"Missing XGBoost model JSON: {xgb_path}")
    booster = xgb.Booster()
    booster.load_model(str(xgb_path))
    return booster, features, chosen


def required_columns(features: list[str], args: argparse.Namespace) -> list[str]:
    cols = set(expand_required_columns(features))
    cols.update(
        [
            "cluster_Et",
            "cluster_Eta",
            "npb_score",
            "cluster_prob",
            "e11_over_e33",
            "e32_over_e35",
            "cluster_et1",
            "cluster_weta_cogx",
            "is_signal",
            args.weight_branch,
            args.reference_branch,
        ]
    )
    return sorted(cols)


def safe_weight(arrays: dict[str, np.ndarray], branch: str, n: int) -> np.ndarray:
    if branch not in arrays:
        return np.ones(n, dtype="float64")
    weight = np.asarray(arrays[branch], dtype="float64")
    weight = np.where(np.isfinite(weight) & (weight > 0.0), weight, 0.0)
    if np.sum(weight) <= 0.0:
        return np.ones(n, dtype="float64")
    return weight


def column_or_default(frame, name: str, default: float) -> np.ndarray:
    if name in frame:
        return frame[name].to_numpy(dtype="float64")
    return np.full(len(frame), float(default), dtype="float64")


def unit_normalize(hist: np.ndarray) -> np.ndarray:
    total = float(np.sum(hist))
    if total <= 0.0 or not math.isfinite(total):
        return hist
    return hist / total


def score_paths(paths: list[Path], booster, features: list[str], args: argparse.Namespace, label: str):
    import pandas as pd
    import uproot
    import xgboost as xgb

    edges = np.linspace(0.0, 1.0, int(args.bins) + 1)
    hist = np.zeros(len(edges) - 1, dtype="float64")
    ref_hist = np.zeros(len(edges) - 1, dtype="float64")
    counts = {
        "files_seen": 0,
        "files_with_tree": 0,
        "rows_seen": 0,
        "rows_after_cuts": 0,
        "weighted_after_cuts": 0.0,
        "reference_rows_after_cuts": 0,
        "truth_signal_rows_after_cuts": 0,
        "truth_background_rows_after_cuts": 0,
        "npb_chunks_applied": 0,
        "npb_chunks_skipped": 0,
        "cluster_prob_chunks_skipped": 0,
    }
    pt_lo, pt_hi = parse_range(args.pt_range)
    columns = required_columns(features, args)
    sample_scales = parse_sample_weight_scales(args.sample_weight_scale)
    for path in paths[: args.max_files or None]:
        counts["files_seen"] += 1
        try:
            with uproot.open(path) as root_file:
                if args.tree not in root_file:
                    continue
                tree = root_file[args.tree]
                have = [col for col in columns if col in tree.keys()]
                missing_features = [name for name in expand_required_columns(features) if name not in have]
                if missing_features:
                    raise SystemExit(f"{path} is missing model feature columns: {missing_features}")
                counts["files_with_tree"] += 1
                for arrays in tree.iterate(have, library="np", step_size=args.step_size):
                    n = len(next(iter(arrays.values()))) if arrays else 0
                    if n == 0:
                        continue
                    counts["rows_seen"] += int(n)
                    frame = pd.DataFrame({name: arrays[name] for name in have})
                    add_derived_features(frame)
                    et = frame["cluster_Et"].to_numpy(dtype="float64")
                    eta = frame["cluster_Eta"].to_numpy(dtype="float64")
                    npb = frame["npb_score"].to_numpy(dtype="float64")
                    prob = column_or_default(frame, "cluster_prob", 0.5)
                    e11e33 = frame["e11_over_e33"].to_numpy(dtype="float64")
                    e32e35 = frame["e32_over_e35"].to_numpy(dtype="float64")
                    et1 = frame["cluster_et1"].to_numpy(dtype="float64")
                    weta = frame["cluster_weta_cogx"].to_numpy(dtype="float64")
                    if "cluster_prob" not in frame:
                        counts["cluster_prob_chunks_skipped"] += 1
                    mask = (
                        np.isfinite(et)
                        & (et >= pt_lo)
                        & (et < pt_hi)
                        & np.isfinite(eta)
                        & (np.abs(eta) < args.eta_max)
                    )
                    if args.preselection == "fig19":
                        mask &= (
                            np.isfinite(prob)
                            & (prob > 0.0)
                            & (prob < 1.0)
                            & np.isfinite(e11e33)
                            & (e11e33 < 0.98)
                            & (e11e33 > 0.0)
                            & np.isfinite(e32e35)
                            & (e32e35 > 0.8)
                            & (e32e35 < 1.0)
                            & np.isfinite(et1)
                            & (et1 > 0.6)
                            & (et1 < 1.0)
                        )
                    if args.require_weta_cut:
                        mask &= np.isfinite(weta) & (weta < 2.0)
                    npb_real = np.isfinite(npb) & (npb >= 0.0) & (npb <= 1.0)
                    apply_npb = args.npb_mode == "apply" or (args.npb_mode == "auto" and bool(np.any(npb_real)))
                    if apply_npb:
                        mask &= npb_real & (npb > args.npb_min)
                        counts["npb_chunks_applied"] += 1
                    else:
                        counts["npb_chunks_skipped"] += 1
                    for feature in features:
                        values = frame[feature].to_numpy(dtype="float64")
                        mask &= np.isfinite(values)
                    if not np.any(mask):
                        continue
                    x = np.column_stack([frame[feature].to_numpy(dtype="float32")[mask] for feature in features])
                    score = booster.predict(xgb.DMatrix(x)).astype("float64")
                    weight_all = safe_weight(frame, args.weight_branch, n)
                    weight = weight_all[mask] * weight_scale_for_path(path, sample_scales)
                    hist += np.histogram(score, bins=edges, weights=weight)[0]
                    counts["rows_after_cuts"] += int(np.sum(mask))
                    counts["weighted_after_cuts"] += float(np.sum(weight))
                    if "is_signal" in frame:
                        truth_label = frame["is_signal"].to_numpy(dtype="int32")[mask]
                        counts["truth_signal_rows_after_cuts"] += int(np.sum(truth_label == 1))
                        counts["truth_background_rows_after_cuts"] += int(np.sum(truth_label == 0))
                    if args.reference_branch in frame:
                        ref = frame[args.reference_branch].to_numpy(dtype="float64")[mask]
                        ref_mask = np.isfinite(ref) & (ref >= 0.0) & (ref <= 1.0)
                        if np.any(ref_mask):
                            ref_hist += np.histogram(ref[ref_mask], bins=edges, weights=weight[ref_mask])[0]
                            counts["reference_rows_after_cuts"] += int(np.sum(ref_mask))
        except OSError as exc:
            print(f"[warn] skipping unreadable {label} file {path}: {exc}", file=sys.stderr)
    if (
        label == "inclusive"
        and args.require_unfiltered_inclusive
        and counts["rows_after_cuts"] > 0
        and counts["truth_signal_rows_after_cuts"] == 0
    ):
        raise SystemExit(
            "Inclusive input appears truth-background filtered: selected rows have "
            "is_signal==0 only. For Shuhang/Fig13-style inclusive-MC overlays, rerun "
            "extraction with RJ_PP_PHOTONID_PPG12_FILTER=0 and "
            "RJ_PP_PHOTONID_REQUIRE_PRESELECTION=0, or use condorExtractCurrentIANRawOverlay."
        )
    if (
        label == "inclusive"
        and args.preselection == "none"
        and counts["rows_after_cuts"] > 0
        and counts["truth_signal_rows_after_cuts"] == 0
    ):
        print(
            "[fig19][WARN] inclusive rows are all is_signal==0; this is a "
            "truth-background-filtered table, not raw inclusive MC.",
            flush=True,
        )
    return {
        "edges": edges,
        "hist": unit_normalize(hist),
        "raw_hist": hist,
        "reference_hist": unit_normalize(ref_hist),
        "reference_raw_hist": ref_hist,
        "counts": counts,
    }


def load_shuhang_projection(path: Path, hist_name: str, edges: np.ndarray) -> dict:
    import uproot

    if path is None:
        raise SystemExit("Internal error: missing Shuhang ROOT path")
    with uproot.open(path) as root_file:
        if hist_name not in root_file:
            matches = [key.split(";")[0] for key in root_file.keys() if "bdt" in key.lower()]
            raise SystemExit(f"{path} is missing {hist_name}. BDT-like keys include: {matches[:20]}")
        hist = root_file[hist_name]
        values, xedges, yedges = hist.to_numpy(flow=False)
    projection = np.asarray(values, dtype="float64").sum(axis=1)
    centers = 0.5 * (np.asarray(xedges[:-1]) + np.asarray(xedges[1:]))
    selected = np.isfinite(centers) & np.isfinite(projection) & (centers >= edges[0]) & (centers < edges[-1])
    rebinned = np.histogram(centers[selected], bins=edges, weights=projection[selected])[0]
    return {
        "hist": unit_normalize(rebinned),
        "raw_hist": rebinned.tolist(),
        "source_root": str(path),
        "hist_name": hist_name,
        "source_x_edges": np.asarray(xedges).tolist(),
        "source_y_edges": np.asarray(yedges).tolist(),
        "source_integral": float(np.sum(projection)),
        "selected_integral": float(np.sum(rebinned)),
    }


def label_suffix(args: argparse.Namespace) -> str:
    pt_lo, pt_hi = parse_range(args.pt_range)
    if args.npb_mode == "skip":
        npb_text = "w/o NPB cut"
    elif args.npb_mode == "apply":
        npb_text = rf"NPB $>$ {args.npb_min:g}"
    else:
        npb_text = rf"auto NPB $>$ {args.npb_min:g}"
    return rf"${pt_lo:g} < E_T < {pt_hi:g}$ GeV, {npb_text}"


def output_tag(args: argparse.Namespace) -> str:
    pt_lo, pt_hi = parse_range(args.pt_range)
    hist_tag = str(args.shuhang_hist).replace("h2d_bdt_", "").replace("_", "-")
    npb_tag = {"skip": "noNPB", "apply": "withNPB", "auto": "autoNPB"}[args.npb_mode]
    return f"pp_noCent_bdt_{pt_lo:g}_{pt_hi:g}_{npb_tag}_{hist_tag}".replace(".", "p")


def styled_label(ax, args: argparse.Namespace):
    ax.text(0.05, 0.94, "sPHENIX", transform=ax.transAxes, fontsize=13, fontweight="bold", fontstyle="italic")
    ax.text(0.245, 0.94, "Internal", transform=ax.transAxes, fontsize=13)
    ax.text(0.05, 0.875, r"$p$+$p$ $\sqrt{s}=200$ GeV", transform=ax.transAxes, fontsize=12)
    ax.text(0.05, 0.815, r"$|\eta| < 0.7$", transform=ax.transAxes, fontsize=12)
    ax.text(0.05, 0.755, label_suffix(args), transform=ax.transAxes, fontsize=12)


def draw_overlay(signal, inclusive, args, model_meta: dict) -> Path:
    import matplotlib.pyplot as plt

    edges = signal["edges"]
    centers = 0.5 * (edges[:-1] + edges[1:])
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    ax.step(centers, signal["hist"], where="mid", color="#d62728", linewidth=1.8, label="Signal MC")
    ax.step(centers, inclusive["hist"], where="mid", color="#1f5fff", linewidth=1.8, label="Inclusive MC")
    styled_label(ax, args)
    ax.set_xlabel("BDT score", fontsize=14)
    ax.set_ylabel("normalized counts", fontsize=14)
    ax.set_xlim(0.0, 1.0)
    ymax = max(float(np.max(signal["hist"])), float(np.max(inclusive["hist"])), 0.05)
    ax.set_ylim(0.0, ymax * 1.25)
    ax.tick_params(direction="in", top=True, right=True, labelsize=12)
    ax.legend(frameon=False, fontsize=12, loc="upper right")
    ax.grid(False)
    fig.tight_layout()
    tag = output_tag(args)
    path = out / f"ppg12_equiv_{tag}_fullInclusive_unitnorm.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)

    summary = {
        "plot": str(path),
        "normalization": "each curve divided by its own weighted bin-count sum, matching PPG12 scaleToUnit-style normalized counts",
        "product": args.product,
        "model": {
            "model_id": model_meta.get("model_id"),
            "product": model_meta.get("product"),
            "features": model_meta.get("features"),
            "output_xgb_json": model_meta.get("output_xgb_json"),
        },
        "cuts": {
            "pt_range": args.pt_range,
            "abs_eta_lt": args.eta_max,
            "npb_score_gt": args.npb_min,
            "npb_mode": args.npb_mode,
            "npb_note": "auto skips chunks with only sentinel/non-real NPB scores such as -2",
            "preselection": args.preselection,
            "cluster_prob_range": [0.0, 1.0],
            "cluster_weta_cogx_lt": 2.0 if args.require_weta_cut else None,
            "e11_over_e33_lt": 0.98 if args.preselection == "fig19" else None,
            "e32_over_e35_range": [0.8, 1.0] if args.preselection == "fig19" else None,
            "cluster_et1_range": [0.6, 1.0] if args.preselection == "fig19" else None,
        },
        "signal": signal["counts"],
        "inclusive": inclusive["counts"],
        "signal_raw_hist": signal["raw_hist"].tolist(),
        "inclusive_raw_hist": inclusive["raw_hist"].tolist(),
        "bins": edges.tolist(),
    }
    (out / f"ppg12_equiv_{tag}_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return path


def draw_shuhang_overlay(signal, inclusive, shuhang_signal, shuhang_inclusive, args, model_meta: dict) -> Path:
    import matplotlib.pyplot as plt

    edges = signal["edges"]
    centers = 0.5 * (edges[:-1] + edges[1:])
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    ax.step(centers, signal["hist"], where="mid", color="#d62728", linewidth=2.0, label="This analysis signal")
    ax.step(centers, inclusive["hist"], where="mid", color="#1f5fff", linewidth=2.0, label="This analysis inclusive")
    ax.step(centers, shuhang_signal["hist"], where="mid", color="#8b1a1a", linewidth=1.8, linestyle="--", label="PPG12/Shuhang signal")
    ax.step(centers, shuhang_inclusive["hist"], where="mid", color="#173a9a", linewidth=1.8, linestyle="--", label="PPG12/Shuhang inclusive")
    styled_label(ax, args)
    ax.set_xlabel("BDT score", fontsize=14)
    ax.set_ylabel("unit-normalized counts", fontsize=14)
    ax.set_xlim(0.0, 1.0)
    ymax = max(
        float(np.max(signal["hist"])),
        float(np.max(inclusive["hist"])),
        float(np.max(shuhang_signal["hist"])),
        float(np.max(shuhang_inclusive["hist"])),
        0.05,
    )
    ax.set_ylim(0.0, ymax * 1.30)
    ax.tick_params(direction="in", top=True, right=True, labelsize=12)
    ax.legend(frameon=False, fontsize=10.5, loc="upper center", ncol=2)
    ax.grid(False)
    fig.tight_layout()
    tag = output_tag(args)
    path = out / f"pp_currentian_basev3e_bdt_score_overlay_vs_ppg12_{tag}.png"
    fig.savefig(path, dpi=240)
    plt.close(fig)

    summary = {
        "plot": str(path),
        "normalization": "each curve divided by its own bin-count sum after applying its own analysis/reference weights",
        "product": args.product,
        "model": {
            "model_id": model_meta.get("model_id"),
            "product": model_meta.get("product"),
            "features": model_meta.get("features"),
            "output_xgb_json": model_meta.get("output_xgb_json"),
        },
        "cuts": {
            "pt_range": args.pt_range,
            "abs_eta_lt": args.eta_max,
            "npb_score_gt": args.npb_min,
            "npb_mode": args.npb_mode,
            "preselection": args.preselection,
            "e11_over_e33_lt": 0.98 if args.preselection == "fig19" else None,
            "e32_over_e35_range": [0.8, 1.0] if args.preselection == "fig19" else None,
            "cluster_et1_range": [0.6, 1.0] if args.preselection == "fig19" else None,
        },
        "sample_weight_scale": parse_sample_weight_scales(args.sample_weight_scale),
        "signal": signal["counts"],
        "inclusive": inclusive["counts"],
        "shuhang_signal": {key: value for key, value in shuhang_signal.items() if key != "hist"},
        "shuhang_inclusive": {key: value for key, value in shuhang_inclusive.items() if key != "hist"},
        "bins": edges.tolist(),
        "this_analysis_signal_hist": signal["hist"].tolist(),
        "this_analysis_inclusive_hist": inclusive["hist"].tolist(),
        "shuhang_signal_hist": shuhang_signal["hist"].tolist(),
        "shuhang_inclusive_hist": shuhang_inclusive["hist"].tolist(),
    }
    (out / f"pp_currentian_basev3e_bdt_score_overlay_vs_ppg12_{tag}_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    return path


def main() -> int:
    args = parse_args()
    signal_paths = expand_paths(args.signal)
    inclusive_paths = expand_paths(args.inclusive)
    if not signal_paths:
        raise SystemExit("No signal input ROOT files found")
    if not inclusive_paths:
        raise SystemExit("No inclusive input ROOT files found")
    booster, features, meta = load_model(args.registry, args.product)
    print(f"[fig19] product={args.product}")
    print(f"[fig19] features={','.join(features)}")
    print(f"[fig19] signal_files={len(signal_paths)} inclusive_files={len(inclusive_paths)}")
    signal = score_paths(signal_paths, booster, features, args, "signal")
    inclusive = score_paths(inclusive_paths, booster, features, args, "inclusive")
    out = draw_overlay(signal, inclusive, args, meta)
    print(f"[fig19] wrote {out}")
    if args.shuhang_signal_root and args.shuhang_inclusive_root:
        shuhang_signal = load_shuhang_projection(args.shuhang_signal_root, args.shuhang_hist, signal["edges"])
        shuhang_inclusive = load_shuhang_projection(args.shuhang_inclusive_root, args.shuhang_hist, signal["edges"])
        shuhang_out = draw_shuhang_overlay(signal, inclusive, shuhang_signal, shuhang_inclusive, args, meta)
        print(f"[fig19] wrote {shuhang_out}")
    print(f"[fig19] signal selected rows={signal['counts']['rows_after_cuts']} inclusive selected rows={inclusive['counts']['rows_after_cuts']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
