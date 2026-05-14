#!/usr/bin/env python3
"""Promote one AuAu BDT+MLP stacker into a production-ready artifact.

This helper intentionally handles one primary stacker at a time.  It trains
from aligned BDT/MLP score caches, derives a held-out WP80 E_T x centrality
threshold grid, and writes diagnostic plots that must be inspected before the
grid is used in RecoilJets production.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import train_auau_stacked_bdt_mlp_sweep as stack  # noqa: E402
from train_auau_photon_mlp import auc_score, threshold_for_signal_efficiency  # noqa: E402


DEFAULT_VARIANT = "ptFine5to35_cent7_full"
DEFAULT_ALGORITHM = "gbm"
DEFAULT_PT_EDGES = "15,17,19,21,23,25,27,30,35"
DEFAULT_CENT_EDGES = "0,10,20,30,40,50,60,80"


def write_json(path: Path, payload) -> None:
    stack.write_json(path, payload)


def parse_edges(text: str) -> list[float]:
    vals = [float(tok) for tok in text.replace(";", ",").split(",") if tok.strip()]
    if len(vals) < 2:
        raise SystemExit(f"Need at least two bin edges, got: {text}")
    if any(not math.isfinite(x) for x in vals):
        raise SystemExit(f"Non-finite bin edge in: {text}")
    if any(vals[i] >= vals[i + 1] for i in range(len(vals) - 1)):
        raise SystemExit(f"Bin edges must be strictly increasing: {text}")
    return vals


def split_masks(y: np.ndarray, seed: int, train_fraction: float, val_fraction: float) -> dict[str, np.ndarray]:
    masks = stack.stratified_split(y, seed, train_fraction, val_fraction)
    masks["trainval"] = masks["train"] | masks["val"]
    return masks


def find_variant(name: str):
    variants = stack.build_variants(SimpleNamespace(sweep="full", include_controls=True))
    for variant in variants:
        if variant.name == name:
            return variant
    known = ", ".join(v.name for v in variants)
    raise SystemExit(f"Unknown stack variant {name}. Known variants: {known}")


def derive_grid(frame, score, mask, target: float, pt_edges: list[float], cent_edges: list[float]) -> dict:
    y = np.asarray(frame["is_signal"], dtype="int32")
    et = np.asarray(frame["cluster_Et"], dtype="float64")
    cent = np.asarray(frame["centrality"], dtype="float64")
    score = np.asarray(score, dtype="float64")
    eligible = mask & np.isin(y, [0, 1]) & np.isfinite(score) & np.isfinite(et) & np.isfinite(cent)
    eligible &= (et >= pt_edges[0]) & (et < pt_edges[-1]) & (cent >= cent_edges[0]) & (cent < cent_edges[-1])
    inclusive = threshold_for_signal_efficiency(y[eligible], score[eligible], target)
    if inclusive is None:
        raise SystemExit("Cannot derive WP grid: held-out split has insufficient signal/background")

    pt_only = {}
    for ip, (plo, phi) in enumerate(zip(pt_edges[:-1], pt_edges[1:])):
        bmask = eligible & (et >= plo) & (et < phi)
        pt_only[ip] = threshold_for_signal_efficiency(y[bmask], score[bmask], target)
    cent_only = {}
    for ic, (clo, chi) in enumerate(zip(cent_edges[:-1], cent_edges[1:])):
        bmask = eligible & (cent >= clo) & (cent < chi)
        cent_only[ic] = threshold_for_signal_efficiency(y[bmask], score[bmask], target)

    thresholds: list[list[float]] = []
    eff_grid: list[list[float]] = []
    fake_grid: list[list[float]] = []
    cells: list[dict] = []
    for ic, (clo, chi) in enumerate(zip(cent_edges[:-1], cent_edges[1:])):
        thr_row = []
        eff_row = []
        fake_row = []
        for ip, (plo, phi) in enumerate(zip(pt_edges[:-1], pt_edges[1:])):
            bmask = eligible & (cent >= clo) & (cent < chi) & (et >= plo) & (et < phi)
            item = threshold_for_signal_efficiency(y[bmask], score[bmask], target)
            source = "cell"
            if item is None:
                item = cent_only.get(ic)
                source = "centrality_fallback"
            if item is None:
                item = pt_only.get(ip)
                source = "pt_fallback"
            if item is None:
                item = inclusive
                source = "inclusive_fallback"
            threshold = float(item["threshold"])
            sig = score[bmask & (y == 1)]
            bkg = score[bmask & (y == 0)]
            sig_eff = float(np.mean(sig > threshold)) if sig.size else math.nan
            bkg_fake = float(np.mean(bkg > threshold)) if bkg.size else math.nan
            thr_row.append(threshold)
            eff_row.append(sig_eff)
            fake_row.append(bkg_fake)
            cells.append(
                {
                    "centrality_min": float(clo),
                    "centrality_max": float(chi),
                    "pt_min": float(plo),
                    "pt_max": float(phi),
                    "threshold": threshold,
                    "signal_efficiency": sig_eff,
                    "background_fake_rate": bkg_fake,
                    "signal_entries": int(sig.size),
                    "background_entries": int(bkg.size),
                    "source": source,
                }
            )
        thresholds.append(thr_row)
        eff_grid.append(eff_row)
        fake_grid.append(fake_row)

    fit_points = []
    fit_values = []
    fit_weights = []
    for cell in cells:
        if not math.isfinite(cell["threshold"]) or cell["signal_entries"] <= 0:
            continue
        fit_points.append([1.0, 0.5 * (cell["pt_min"] + cell["pt_max"]), 0.5 * (cell["centrality_min"] + cell["centrality_max"])])
        fit_values.append(cell["threshold"])
        fit_weights.append(math.sqrt(max(1, cell["signal_entries"])))
    plane = {"intercept": math.nan, "pt_slope": math.nan, "centrality_slope": math.nan, "max_abs_residual": math.nan}
    residuals = []
    if len(fit_points) >= 3:
        x = np.asarray(fit_points, dtype=float)
        t = np.asarray(fit_values, dtype=float)
        w = np.asarray(fit_weights, dtype=float)
        coeff, *_ = np.linalg.lstsq(x * w[:, None], t * w, rcond=None)
        pred = x @ coeff
        residuals = (t - pred).tolist()
        plane = {
            "intercept": float(coeff[0]),
            "pt_slope": float(coeff[1]),
            "centrality_slope": float(coeff[2]),
            "max_abs_residual": float(np.max(np.abs(t - pred))) if t.size else math.nan,
        }

    finite_eff = np.asarray([c["signal_efficiency"] for c in cells if math.isfinite(c["signal_efficiency"])], dtype=float)
    return {
        "target_signal_efficiency": float(target),
        "mode": "grid2d",
        "pt_min": float(pt_edges[0]),
        "pt_max": float(pt_edges[-1]),
        "max_score": 1.0,
        "pt_edges": [float(x) for x in pt_edges],
        "cent_edges": [float(x) for x in cent_edges],
        "grid_thresholds": thresholds,
        "grid_signal_efficiency": eff_grid,
        "grid_background_fake_rate": fake_grid,
        "cells": cells,
        "inclusive": inclusive,
        "fit_quality": {
            "max_abs_cell_efficiency_error": float(np.max(np.abs(finite_eff - target))) if finite_eff.size else math.nan,
            "mean_abs_cell_efficiency_error": float(np.mean(np.abs(finite_eff - target))) if finite_eff.size else math.nan,
            "min_cell_signal_entries": int(min((c["signal_entries"] for c in cells), default=0)),
            "min_cell_background_entries": int(min((c["background_entries"] for c in cells), default=0)),
            "plane_fit": plane,
            "plane_residuals": residuals,
        },
    }


def runtime_entry(variant: str, wp: dict) -> str:
    pt_edges = ";".join(f"{x:.10g}" for x in wp["pt_edges"])
    cent_edges = ";".join(f"{x:.10g}" for x in wp["cent_edges"])
    thresholds = ";".join(f"{x:.10g}" for row in wp["grid_thresholds"] for x in row)
    return (
        f"{variant}|grid2d|{pt_edges}|{cent_edges}|{thresholds}|"
        f"{wp['pt_min']:.10g}|{wp['pt_max']:.10g}|{wp['max_score']:.10g}"
    )


def write_cell_csv(path: Path, cells: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "centrality_min",
        "centrality_max",
        "pt_min",
        "pt_max",
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
        writer.writerows(cells)


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def make_diagnostics(outdir: Path, wp: dict, target: float, model_name: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outdir.mkdir(parents=True, exist_ok=True)
    thresholds = np.asarray(wp["grid_thresholds"], dtype=float)
    eff = np.asarray(wp["grid_signal_efficiency"], dtype=float)
    fake = np.asarray(wp["grid_background_fake_rate"], dtype=float)
    pt_edges = np.asarray(wp["pt_edges"], dtype=float)
    cent_edges = np.asarray(wp["cent_edges"], dtype=float)
    extent = [pt_edges[0], pt_edges[-1], cent_edges[-1], cent_edges[0]]

    def heat(data, title, cbar, path, cmap="viridis"):
        fig, ax = plt.subplots(figsize=(8.2, 5.2))
        im = ax.imshow(data, aspect="auto", interpolation="nearest", extent=extent, cmap=cmap)
        fig.colorbar(im, ax=ax, label=cbar)
        ax.set_xlabel(r"Photon candidate $E_T$ [GeV]")
        ax.set_ylabel("Centrality [%]")
        ax.set_title(title, pad=24)
        ax.text(0.0, 1.02, r"$\bf{\it{sPHENIX}}$ Internal", transform=ax.transAxes, ha="left", va="bottom", fontsize=12, clip_on=False)
        for ic in range(data.shape[0]):
            for ip in range(data.shape[1]):
                x = 0.5 * (pt_edges[ip] + pt_edges[ip + 1])
                y = 0.5 * (cent_edges[ic] + cent_edges[ic + 1])
                val = data[ic, ip]
                if math.isfinite(float(val)):
                    rgba = im.cmap(im.norm(float(val)))
                    luminance = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
                    text_color = "black" if luminance > 0.58 else "white"
                    ax.text(x, y, f"{val:.3f}", ha="center", va="center", fontsize=7, color=text_color)
        fig.tight_layout()
        fig.savefig(path, dpi=180)
        plt.close(fig)

    heat(thresholds, f"{model_name}: WP{int(round(target*100))} stack threshold grid", "Score threshold", outdir / "stack_wp80_threshold_grid.png")
    heat(eff, f"{model_name}: achieved signal efficiency", "Signal efficiency", outdir / "stack_wp80_signal_efficiency_grid.png", "magma")
    heat(fake, f"{model_name}: background fake rate at WP{int(round(target*100))}", "Background fake rate", outdir / "stack_wp80_fake_rate_grid.png", "cividis")

    plane = wp["fit_quality"].get("plane_fit", {})
    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    centers = 0.5 * (pt_edges[:-1] + pt_edges[1:])
    for ic, (clo, chi) in enumerate(zip(cent_edges[:-1], cent_edges[1:])):
        ax.plot(centers, thresholds[ic], marker="o", lw=1.7, label=f"{clo:g}-{chi:g}%")
        if all(math.isfinite(float(plane.get(k, math.nan))) for k in ("intercept", "pt_slope", "centrality_slope")):
            c = 0.5 * (clo + chi)
            fit = plane["intercept"] + plane["pt_slope"] * centers + plane["centrality_slope"] * c
            ax.plot(centers, fit, ls="--", lw=1.0, alpha=0.65)
    ax.set_xlabel(r"Photon candidate $E_T$ [GeV]")
    ax.set_ylabel(f"Score threshold for {100*target:.0f}% signal efficiency")
    ax.set_title(f"{model_name}: local thresholds with plane-fit overlay", pad=24)
    ax.text(0.0, 1.02, r"$\bf{\it{sPHENIX}}$ Internal", transform=ax.transAxes, ha="left", va="bottom", fontsize=12, clip_on=False)
    ax.legend(ncol=2, fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / "stack_wp80_threshold_vs_et_by_centrality.png", dpi=180)
    plt.close(fig)


def train_and_wp(args: argparse.Namespace) -> None:
    args.outdir.mkdir(parents=True, exist_ok=True)
    sweep_args = SimpleNamespace(
        mlp_cache=args.mlp_cache,
        bdt_cache=args.bdt_cache,
        mlp_score=args.mlp_score,
        bdt_score=args.bdt_score,
        max_rows=args.max_rows,
        max_shards=args.max_shards,
        allow_missing_full_features=args.allow_missing_full_features,
        outdir=args.outdir,
        l2=args.l2,
        linear_backend=args.linear_backend,
        max_linear_steps=args.max_linear_steps,
        gbm_estimators=args.gbm_estimators,
        gbm_learning_rate=args.gbm_learning_rate,
        gbm_max_depth=args.gbm_max_depth,
        gbm_max_leaf_nodes=args.gbm_max_leaf_nodes,
        nn_hidden=args.nn_hidden,
        nn_epochs=args.nn_epochs,
        nn_patience=args.nn_patience,
        nn_batch_size=args.nn_batch_size,
        nn_learning_rate=args.nn_learning_rate,
        nn_l2=args.nn_l2,
    )
    frame, cache_info = stack.load_aligned_caches(sweep_args)
    preflight = stack.preflight(frame, cache_info, sweep_args)
    variant = find_variant(args.variant)
    y = np.asarray(frame["is_signal"], dtype="int32")
    masks = split_masks(y, args.random_seed, args.train_fraction, args.val_fraction)
    feature_names, _missing = stack.feature_names_for_variant(frame, variant, sweep_args)
    x = stack.feature_matrix(frame, feature_names)
    score = np.full(len(y), np.nan, dtype="float64")
    eval_mask = stack.pt_range_mask(frame, variant.pt_lo, variant.pt_hi)
    training_history_rows = []
    if variant.routes:
        artifact_routes = []
        for idx, route in enumerate(variant.routes):
            route_mask = eval_mask & stack.route_mask(frame, route)
            route_train = masks[args.fit_split] & route_mask
            route_val = masks["val"] & route_mask
            route_fit = stack.fit_one(f"{variant.name}_{args.algorithm}_{route.label}", args.algorithm, feature_names, x, y, route_train, sweep_args, args.random_seed + 2000 + idx, val_mask=route_val)
            if route_fit is None:
                continue
            score[route_mask] = route_fit.predict(x[route_mask])
            for row in route_fit.history:
                hist = {"model": f"{variant.name}_{args.algorithm}", "route": route.label, "submodel": f"{variant.name}_{args.algorithm}_{route.label}", "algorithm": args.algorithm}
                hist.update(row)
                training_history_rows.append(hist)
            artifact_routes.append(
                {
                    "label": route.label,
                    "pt_lo": route.pt_lo,
                    "pt_hi": route.pt_hi,
                    "cent_lo": route.cent_lo,
                    "cent_hi": route.cent_hi,
                    "model": route_fit.artifact,
                }
            )
        if not artifact_routes:
            raise SystemExit("No routed stack experts were fitted")
        artifact = {
            "schema": "RJ_AUAU_STACKED_BDT_MLP_ARTIFACT_V1",
            "name": f"{variant.name}_{args.algorithm}",
            "variant": stack.variant_payload(variant),
            "abcd_eligible_isolation_input_free": True,
            "closure_pending": True,
            "uses_runtime_bdt_score": True,
            "uses_runtime_mlp_score": True,
            "algorithm": args.algorithm,
            "feature_names": feature_names,
            "routes": artifact_routes,
        }
    else:
        train_mask = masks[args.fit_split] & eval_mask
        val_mask = masks["val"] & eval_mask
        fitted = stack.fit_one(f"{variant.name}_{args.algorithm}", args.algorithm, feature_names, x, y, train_mask, sweep_args, args.random_seed + 1009, val_mask=val_mask)
        if fitted is None:
            raise SystemExit("Stack model did not fit")
        score[eval_mask] = fitted.predict(x[eval_mask])
        for row in fitted.history:
            hist = {"model": f"{variant.name}_{args.algorithm}", "route": "global", "submodel": f"{variant.name}_{args.algorithm}", "algorithm": args.algorithm}
            hist.update(row)
            training_history_rows.append(hist)
        artifact = {
            "schema": "RJ_AUAU_STACKED_BDT_MLP_ARTIFACT_V1",
            "name": f"{variant.name}_{args.algorithm}",
            "variant": stack.variant_payload(variant),
            "abcd_eligible_isolation_input_free": True,
            "closure_pending": True,
            "uses_runtime_bdt_score": True,
            "uses_runtime_mlp_score": True,
            "algorithm": args.algorithm,
            "feature_names": feature_names,
            "model": fitted.artifact,
        }

    artifact_path = args.outdir / "artifacts" / f"{variant.name}_{args.algorithm}.json"
    write_json(artifact_path, artifact)
    history_path = args.outdir / "stack_training_history.csv"
    write_csv(history_path, training_history_rows)
    np.savez_compressed(
        args.outdir / "stack_scores.npz",
        score=score.astype("float32"),
        is_signal=y.astype("int8"),
        cluster_Et=np.asarray(frame["cluster_Et"], dtype="float32"),
        centrality=np.asarray(frame["centrality"], dtype="float32"),
        reco_eiso=np.asarray(frame.get("reco_eiso", np.full(len(y), np.nan)), dtype="float32"),
        split_train=masks["train"],
        split_val=masks["val"],
        split_test=masks["test"],
    )

    pt_edges = parse_edges(args.wp_pt_edges)
    cent_edges = parse_edges(args.wp_cent_edges)
    wp_mask = masks[args.wp_split]
    wp = derive_grid(frame, score, wp_mask, args.target_signal_efficiency, pt_edges, cent_edges)
    manifest = {
        "schema": "RECOILJETS_AUAU_BDT_MLP_STACK_WORKING_POINTS_V1",
        "target_signal_efficiency": args.target_signal_efficiency,
        "working_point_mode": "centpt",
        "model_name": f"{variant.name}_{args.algorithm}",
        "artifact": str(artifact_path),
        "source_caches": {"mlp_cache": str(args.mlp_cache), "bdt_cache": str(args.bdt_cache)},
        "preflight": preflight,
        "fit_split": args.fit_split,
        "wp_split": args.wp_split,
        "runtime_entries": [runtime_entry("auauBDTMLPStack", wp)],
        "products": {"auauBDTMLPStack": wp},
    }
    wp_path = args.outdir / "stack_working_points_target80.json"
    write_json(wp_path, manifest)
    (args.outdir / "stack_working_points_target80.yaml").write_text(
        "auau_tight_bdt_mlp_stack_working_point_entries: ["
        + ", ".join(json.dumps(x) for x in manifest["runtime_entries"])
        + "]\n"
    )
    write_cell_csv(args.outdir / "stack_working_points_target80_cells.csv", wp["cells"])
    make_diagnostics(args.outdir / "diagnostics", wp, args.target_signal_efficiency, f"{variant.name}_{args.algorithm}")

    eval_finite = eval_mask & np.isfinite(score) & np.isin(y, [0, 1])
    wp_finite = wp_mask & eval_finite
    summary = {
        "status": "READY",
        "model": f"{variant.name}_{args.algorithm}",
        "artifact": str(artifact_path),
        "working_points": str(wp_path),
        "training_history_csv": str(history_path),
        "rows": int(len(y)),
        "eval_rows": int(eval_mask.sum()),
        "finite_eval_fraction": float(np.mean(np.isfinite(score[eval_mask]))) if eval_mask.any() else math.nan,
        "auc_eval_all": auc_score(y[eval_finite], score[eval_finite]),
        "auc_wp_split": auc_score(y[wp_finite], score[wp_finite]),
        "fit_quality": wp["fit_quality"],
        "diagnostics_dir": str(args.outdir / "diagnostics"),
    }
    write_json(args.outdir / "promotion_summary.json", summary)
    print("[stackPromote] READY")
    print(f"[stackPromote] artifact={artifact_path}")
    print(f"[stackPromote] working_points={wp_path}")
    print(f"[stackPromote] diagnostics={args.outdir / 'diagnostics'}")


def replace_top_level_block(lines: list[str], key: str, replacement: list[str]) -> list[str]:
    start = None
    for i, line in enumerate(lines):
        if line.startswith(f"{key}:"):
            start = i
            break
    if start is None:
        return lines + ["\n"] + replacement
    end = start + 1
    while end < len(lines):
        line = lines[end]
        if line and not line.startswith((" ", "\t", "#")) and ":" in line:
            break
        end += 1
    return lines[:start] + replacement + lines[end:]


def make_config(args: argparse.Namespace) -> None:
    wp = json.loads(args.working_points.read_text())
    entry = wp["runtime_entries"][0]
    artifact = args.artifact or wp.get("artifact", "")
    if not artifact:
        raise SystemExit("No stack artifact provided")
    lines = args.template.read_text().splitlines()
    lines = replace_top_level_block(lines, "photon_id_sets", ["photon_id_sets:", "  - [reference, auauBDTMLPStack, auauBDTMLPStackComplement]"])
    replacements = {
        "auau_tight_bdt_mlp_stack_model_file": artifact,
        "auau_tight_bdt_mlp_stack_bdt_mode": args.bdt_mode,
        "auau_tight_bdt_mlp_stack_mlp_mode": args.mlp_mode,
        "auau_tight_bdt_mlp_stack_apply_pt_min": "15.0",
        "auau_tight_bdt_mlp_stack_apply_pt_max": "35.0",
        "auau_tight_bdt_mlp_stack_working_point_entries": "[" + json.dumps(entry) + "]",
    }
    seen = set()
    out_lines = []
    for line in lines:
        key = line.split(":", 1)[0].strip() if ":" in line and not line.startswith((" ", "\t", "#")) else ""
        if key in replacements:
            out_lines.append(f"{key}: {replacements[key]}")
            seen.add(key)
        else:
            out_lines.append(line)
    for key, value in replacements.items():
        if key not in seen:
            out_lines.append(f"{key}: {value}")
    args.out.write_text("\n".join(out_lines) + "\n")
    print(f"[stackPromote] wrote_config={args.out}")


def self_test(args: argparse.Namespace) -> None:
    if args.outdir.exists():
        shutil.rmtree(args.outdir)
    args.outdir.mkdir(parents=True, exist_ok=True)
    mlp_cache, bdt_cache = stack.make_synthetic_cache(args.outdir / "synthetic", 9000, args.random_seed)
    args.mlp_cache = mlp_cache
    args.bdt_cache = bdt_cache
    args.max_rows = 0
    args.max_shards = 0
    train_and_wp(args)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="mode", required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--mlp-cache", type=Path)
    common.add_argument("--bdt-cache", type=Path)
    common.add_argument("--outdir", type=Path, required=True)
    common.add_argument("--mlp-score", default=stack.DEFAULT_MLP_SCORE)
    common.add_argument("--bdt-score", default=stack.DEFAULT_BDT_SCORE)
    common.add_argument("--variant", default=DEFAULT_VARIANT)
    common.add_argument("--algorithm", choices=["gbm", "logistic", "nn"], default=DEFAULT_ALGORITHM)
    common.add_argument("--target-signal-efficiency", type=float, default=0.80)
    common.add_argument("--wp-pt-edges", default=DEFAULT_PT_EDGES)
    common.add_argument("--wp-cent-edges", default=DEFAULT_CENT_EDGES)
    common.add_argument("--fit-split", choices=["train", "trainval"], default="train")
    common.add_argument("--wp-split", choices=["val", "test"], default="val")
    common.add_argument("--train-fraction", type=float, default=0.60)
    common.add_argument("--val-fraction", type=float, default=0.20)
    common.add_argument("--random-seed", type=int, default=24681357)
    common.add_argument("--max-shards", type=int, default=0)
    common.add_argument("--max-rows", type=int, default=0)
    common.add_argument("--allow-missing-full-features", action="store_true")
    common.add_argument("--l2", type=float, default=2.0e-3)
    common.add_argument("--linear-backend", choices=["numpy", "sklearn"], default="numpy")
    common.add_argument("--max-linear-steps", type=int, default=1800)
    common.add_argument("--gbm-estimators", type=int, default=90)
    common.add_argument("--gbm-learning-rate", type=float, default=0.045)
    common.add_argument("--gbm-max-depth", type=int, default=3)
    common.add_argument("--gbm-max-leaf-nodes", type=int, default=8)
    common.add_argument("--nn-hidden", default="64,32")
    common.add_argument("--nn-epochs", type=int, default=180)
    common.add_argument("--nn-patience", type=int, default=24)
    common.add_argument("--nn-batch-size", type=int, default=8192)
    common.add_argument("--nn-learning-rate", type=float, default=1.5e-3)
    common.add_argument("--nn-l2", type=float, default=1.0e-3)
    sub.add_parser("train-wp", parents=[common])
    sub.add_parser("self-test", parents=[common])
    cfg = sub.add_parser("make-config")
    cfg.add_argument("--template", type=Path, required=True)
    cfg.add_argument("--working-points", type=Path, required=True)
    cfg.add_argument("--artifact", default="")
    cfg.add_argument("--out", type=Path, required=True)
    cfg.add_argument("--bdt-mode", default="auauEtFineCent7BDT")
    cfg.add_argument("--mlp-mode", default="auauCentInputBase3x3MLP")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    if args.mode == "self-test":
        self_test(args)
    elif args.mode == "train-wp":
        if not args.mlp_cache or not args.bdt_cache:
            raise SystemExit("train-wp requires --mlp-cache and --bdt-cache")
        train_and_wp(args)
    elif args.mode == "make-config":
        make_config(args)
    else:
        raise SystemExit(f"Unknown mode {args.mode}")


if __name__ == "__main__":
    main()
