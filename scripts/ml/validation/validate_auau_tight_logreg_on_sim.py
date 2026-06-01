#!/usr/bin/env python3
"""Validate AuAu tight-logistic-regression photon-ID artifacts."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from train_auau_photon_logreg import (  # noqa: E402
    PRODUCT_SPECS,
    SCHEMA,
    artifact_score,
    parse_products,
    write_json,
)
from train_auau_photon_mlp import (  # noqa: E402
    add_derived_features,
    auc_score,
    calibration_report,
    expand_input_paths,
    expand_required_columns,
    json_ready,
    limit_paths_by_group,
    load_frame,
    parse_bin_spec,
    threshold_for_signal_efficiency,
)


DEFAULT_PT_BINS = "15,17,19,21,23,25,27,30,35"
DEFAULT_CENT_BINS = "0,10,20,30,40,50,60,80"


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def load_artifact(path: Path) -> dict:
    data = json.loads(path.read_text())
    if data.get("schema") != SCHEMA:
        raise SystemExit(f"Unsupported logreg artifact schema in {path}: {data.get('schema')}")
    return data


def discover_artifacts(model_dir: Path, model_registry: Path | None = None) -> list[tuple[str, Path, dict]]:
    registry = model_registry or (model_dir / "model_registry.json")
    items: list[tuple[str, Path, dict]] = []
    if registry.is_file():
        data = json.loads(registry.read_text())
        for row in data.get("models", []):
            product = row.get("product") or Path(row.get("output_json", "")).stem
            raw = Path(row["output_json"])
            path = raw if raw.is_absolute() else model_dir / raw.name
            items.append((product, path, load_artifact(path)))
    else:
        for path in sorted(model_dir.glob("auau_tight_logreg_*.json")):
            art = load_artifact(path)
            items.append((str(art.get("name") or path.stem), path, art))
    if not items:
        raise SystemExit(f"No logreg artifacts found in {model_dir}")
    return items


def artifact_features(artifact: dict) -> list[str]:
    features: list[str] = []
    for feature in artifact.get("features", []):
        if feature not in features:
            features.append(feature)
    if "model" in artifact:
        for feature in artifact["model"].get("feature_names", []):
            if feature not in features:
                features.append(feature)
    for route in artifact.get("routes", []):
        for feature in route.get("model", {}).get("feature_names", []):
            if feature not in features:
                features.append(feature)
    return features


def validation_bins(text: str | None, default: str) -> list[tuple[float, float]]:
    return parse_bin_spec(text or default)


def score_key(product: str) -> str:
    clean = "".join(ch if ch.isalnum() or ch == "_" else "_" for ch in product)
    return f"score_{clean}"


def load_validation_frame(args, artifacts: list[tuple[str, Path, dict]]):
    required = {"is_signal", "cluster_Et", "centrality", "cluster_Eta", "vertexz", "run", "evt", "reco_eiso"}
    for _, _, artifact in artifacts:
        required.update(expand_required_columns(artifact_features(artifact)))
    paths = expand_input_paths(args.input, args.source, check_exists=False)
    if args.manifest:
        paths = expand_input_paths([Path("@" + str(args.manifest))], None, check_exists=False)
    paths, manifest_report = limit_paths_by_group(paths, args.max_files_per_sample, args.random_seed)
    frame, seen_optional, read_report = load_frame(paths, args.tree, sorted(required), [])
    frame = add_derived_features(frame)
    return frame, {"manifest_report": manifest_report, "seen_optional": seen_optional, "read_report": read_report}


def finite_corr(a, b) -> float:
    a = np.asarray(a, dtype="float64")
    b = np.asarray(b, dtype="float64")
    mask = np.isfinite(a) & np.isfinite(b)
    if int(mask.sum()) < 3:
        return math.nan
    return float(np.corrcoef(a[mask], b[mask])[0, 1])


def bin_metrics(frame, score, bins, axis: str, target: float) -> dict:
    y = frame["is_signal"].to_numpy(dtype="int32")
    values = frame[axis].to_numpy(dtype="float64")
    out = {}
    for lo, hi in bins:
        mask = np.isfinite(values) & (values >= lo) & (values < hi) & np.isfinite(score) & np.isin(y, [0, 1])
        key = f"{lo:g}_{hi:g}".replace(".", "p")
        wp = threshold_for_signal_efficiency(y[mask], score[mask], target)
        out[key] = {
            "lo": float(lo),
            "hi": float(hi),
            "entries": int(mask.sum()),
            "signal_entries": int((mask & (y == 1)).sum()),
            "background_entries": int((mask & (y == 0)).sum()),
            "auc": auc_score(y[mask], score[mask]) if mask.any() else math.nan,
            "wp80_fake": None if wp is None else float(wp["background_fake_rate"]),
            "wp80_threshold": None if wp is None else float(wp["threshold"]),
        }
    return out


def metric_row(product: str, frame, score, pt_bins, cent_bins, target: float) -> dict:
    y = frame["is_signal"].to_numpy(dtype="int32")
    finite = np.isfinite(score) & np.isin(y, [0, 1])
    wp = threshold_for_signal_efficiency(y[finite], score[finite], target)
    pt_report = bin_metrics(frame, score, pt_bins, "cluster_Et", target)
    cent_report = bin_metrics(frame, score, cent_bins, "centrality", target)
    row = {
        "product": product,
        "entries": int(len(frame)),
        "scored_entries": int(finite.sum()),
        "finite_fraction": float(finite.mean()) if len(finite) else math.nan,
        "auc": auc_score(y[finite], score[finite]) if finite.any() else math.nan,
        "wp80_fake": None if wp is None else float(wp["background_fake_rate"]),
        "wp80_threshold": None if wp is None else float(wp["threshold"]),
        "ece": calibration_report(y[finite], score[finite]).get("ece") if finite.any() else math.nan,
        "score_vs_eiso_corr": finite_corr(score, frame["reco_eiso"].to_numpy(dtype="float64")),
        "pt_bins": pt_report,
        "centrality_bins": cent_report,
    }
    for key, item in pt_report.items():
        row[f"pt_{key}_auc"] = item["auc"]
        row[f"pt_{key}_wp80_fake"] = item["wp80_fake"]
    for key, item in cent_report.items():
        row[f"cent_{key}_auc"] = item["auc"]
        row[f"cent_{key}_wp80_fake"] = item["wp80_fake"]
    return row


def derive_grid(frame, score, target: float, pt_bins, cent_bins) -> dict:
    y = frame["is_signal"].to_numpy(dtype="int32")
    et = frame["cluster_Et"].to_numpy(dtype="float64")
    cent = frame["centrality"].to_numpy(dtype="float64")
    finite = np.isfinite(score) & np.isfinite(et) & np.isfinite(cent) & np.isin(y, [0, 1])
    inclusive = threshold_for_signal_efficiency(y[finite], score[finite], target)
    if inclusive is None:
        raise SystemExit("Cannot derive WP grid: no valid inclusive threshold")
    pt_edges = [pt_bins[0][0]] + [hi for _, hi in pt_bins]
    cent_edges = [cent_bins[0][0]] + [hi for _, hi in cent_bins]
    thresholds = []
    eff_grid = []
    fake_grid = []
    cells = []
    for clo, chi in cent_bins:
        thr_row = []
        eff_row = []
        fake_row = []
        for plo, phi in pt_bins:
            mask = finite & (et >= plo) & (et < phi) & (cent >= clo) & (cent < chi)
            wp = threshold_for_signal_efficiency(y[mask], score[mask], target)
            source = "cell"
            if wp is None:
                cmask = finite & (cent >= clo) & (cent < chi)
                wp = threshold_for_signal_efficiency(y[cmask], score[cmask], target)
                source = "centrality_fallback"
            if wp is None:
                pmask = finite & (et >= plo) & (et < phi)
                wp = threshold_for_signal_efficiency(y[pmask], score[pmask], target)
                source = "pt_fallback"
            if wp is None:
                wp = inclusive
                source = "inclusive_fallback"
            threshold = float(wp["threshold"])
            sig = score[mask & (y == 1)]
            bkg = score[mask & (y == 0)]
            sig_eff = float(np.mean(sig > threshold)) if sig.size else math.nan
            fake = float(np.mean(bkg > threshold)) if bkg.size else math.nan
            thr_row.append(threshold)
            eff_row.append(sig_eff)
            fake_row.append(fake)
            cells.append(
                {
                    "centrality_min": float(clo),
                    "centrality_max": float(chi),
                    "pt_min": float(plo),
                    "pt_max": float(phi),
                    "threshold": threshold,
                    "signal_efficiency": sig_eff,
                    "background_fake_rate": fake,
                    "signal_entries": int(sig.size),
                    "background_entries": int(bkg.size),
                    "source": source,
                }
            )
        thresholds.append(thr_row)
        eff_grid.append(eff_row)
        fake_grid.append(fake_row)
    finite_eff = np.asarray([c["signal_efficiency"] for c in cells if math.isfinite(c["signal_efficiency"])])
    return {
        "target_signal_efficiency": float(target),
        "working_point_mode": "centpt",
        "mode": "grid2d",
        "pt_min": float(pt_edges[0]),
        "pt_max": float(pt_edges[-1]),
        "max_score": 1.0,
        "pt_edges": pt_edges,
        "cent_edges": cent_edges,
        "grid_thresholds": thresholds,
        "grid_signal_efficiency": eff_grid,
        "grid_background_fake_rate": fake_grid,
        "cells": cells,
        "fit_quality": {
            "min_cell_signal_entries": int(min((c["signal_entries"] for c in cells), default=0)),
            "min_cell_background_entries": int(min((c["background_entries"] for c in cells), default=0)),
            "max_abs_cell_efficiency_error": float(np.max(np.abs(finite_eff - target))) if finite_eff.size else math.nan,
        },
    }


def runtime_entry(variant: str, wp: dict) -> str:
    pt_edges = ";".join(f"{x:.10g}" for x in wp["pt_edges"])
    cent_edges = ";".join(f"{x:.10g}" for x in wp["cent_edges"])
    thresholds = ";".join(f"{x:.10g}" for row in wp["grid_thresholds"] for x in row)
    return f"{variant}|grid2d|{pt_edges}|{cent_edges}|{thresholds}|{wp['pt_min']:.10g}|{wp['pt_max']:.10g}|{wp['max_score']:.10g}"


def write_working_points(outdir: Path, metrics: dict, frame, scores: dict[str, np.ndarray], target: float, pt_bins, cent_bins):
    products = {}
    for product, score in scores.items():
        grid = derive_grid(frame, score, target, pt_bins, cent_bins)
        entry = runtime_entry("auauTightLogReg", grid)
        products[product] = {"product": product, "variant": "auauTightLogReg", "entry": entry, "working_point": grid}
    payload = {"schema": "RJ_AUAU_TIGHT_LOGREG_WORKING_POINTS_V1", "target_signal_efficiency": target, "products": products}
    label = f"target{int(round(target * 100)):02d}"
    path = outdir / f"logreg_working_points_{label}.json"
    write_json(path, payload)
    quoted = ",".join(json.dumps(v["entry"]) for v in products.values())
    (outdir / f"logreg_working_points_{label}.yaml").write_text(f"auau_tight_logreg_working_point_entries: [{quoted}]\n")
    cell_rows = []
    for product, item in products.items():
        for cell in item["working_point"]["cells"]:
            row = {"product": product}
            row.update(cell)
            cell_rows.append(row)
    write_csv(outdir / f"logreg_working_points_{label}_cells.csv", cell_rows)
    return path


def make_plots(outdir: Path, frame, scores: dict[str, np.ndarray], metrics: dict, wp_path: Path | None) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"[validateAuAuTightLogReg][WARN] matplotlib unavailable; skipping plots: {exc}", file=sys.stderr)
        return
    plot_dir = outdir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    y = frame["is_signal"].to_numpy(dtype="int32")
    for product, score in scores.items():
        finite = np.isfinite(score) & np.isin(y, [0, 1])
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.hist(score[finite & (y == 1)], bins=60, range=(0, 1), histtype="step", density=True, label="Signal")
        ax.hist(score[finite & (y == 0)], bins=60, range=(0, 1), histtype="step", density=True, label="Background")
        ax.set_xlabel("logistic-regression score")
        ax.set_ylabel("Density")
        ax.set_title(f"{product} score separation")
        ax.legend()
        fig.tight_layout()
        fig.savefig(plot_dir / f"{product}_score_separation.png", dpi=180)
        plt.close(fig)

    if wp_path and wp_path.is_file():
        wp = json.loads(wp_path.read_text())
        for product, item in wp.get("products", {}).items():
            grid = item["working_point"]
            for key, title, cmap in [
                ("grid_thresholds", "WP80 threshold", "viridis"),
                ("grid_signal_efficiency", "Achieved signal efficiency", "magma"),
                ("grid_background_fake_rate", "Background fake rate", "cividis"),
            ]:
                data = np.asarray(grid[key], dtype=float)
                fig, ax = plt.subplots(figsize=(8, 5))
                im = ax.imshow(data, aspect="auto", interpolation="nearest", cmap=cmap)
                fig.colorbar(im, ax=ax)
                ax.set_title(f"{product}: {title}")
                ax.set_xlabel("E_T bin index")
                ax.set_ylabel("centrality bin index")
                for ic in range(data.shape[0]):
                    for ip in range(data.shape[1]):
                        ax.text(ip, ic, f"{data[ic, ip]:.3f}", ha="center", va="center", fontsize=7, color="white")
                fig.tight_layout()
                fig.savefig(plot_dir / f"{product}_{key}.png", dpi=180)
                plt.close(fig)


def save_score_cache(path: Path, frame, scores: dict[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    common = {
        "is_signal": frame["is_signal"].to_numpy(dtype="int8"),
        "cluster_Et": frame["cluster_Et"].to_numpy(dtype="float32"),
        "centrality": frame["centrality"].to_numpy(dtype="float32"),
        "cluster_Eta": frame["cluster_Eta"].to_numpy(dtype="float32"),
        "reco_eiso": frame["reco_eiso"].to_numpy(dtype="float32"),
    }
    common.update({score_key(k): v.astype("float32") for k, v in scores.items()})
    np.savez_compressed(path, **common)
    print(f"[validateAuAuTightLogReg] wrote score cache: {path}", flush=True)


def load_score_caches(manifest: Path):
    paths = [Path(x.strip()) for x in manifest.read_text().splitlines() if x.strip()]
    if not paths:
        raise SystemExit(f"No score caches listed in {manifest}")
    pieces: dict[str, list[np.ndarray]] = {}
    for idx, path in enumerate(paths, 1):
        if idx == 1 or idx == len(paths) or idx % 20 == 0:
            print(f"[validateAuAuTightLogReg] merging cache {idx}/{len(paths)}: {path}", flush=True)
        with np.load(path) as cache:
            for key in cache.files:
                pieces.setdefault(key, []).append(np.asarray(cache[key]))
    merged = {key: np.concatenate(vals) for key, vals in pieces.items()}
    import pandas as pd

    frame = pd.DataFrame({k: merged[k] for k in ["is_signal", "cluster_Et", "centrality", "cluster_Eta", "reco_eiso"] if k in merged})
    scores = {key.removeprefix("score_"): val.astype("float64") for key, val in merged.items() if key.startswith("score_")}
    return frame, scores, {"score_caches": len(paths)}


def run_validation(args) -> int:
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    if args.merge_score_caches:
        frame, scores, read_report = load_score_caches(args.merge_score_caches)
        artifacts = []
    else:
        artifacts = discover_artifacts(args.model_dir, args.model_registry)
        frame, read_report = load_validation_frame(args, artifacts)
        if args.score_max_rows and args.score_max_rows > 0 and len(frame) > args.score_max_rows:
            rng = np.random.default_rng(args.random_seed)
            idx = np.sort(rng.choice(np.arange(len(frame)), size=args.score_max_rows, replace=False))
            frame = frame.iloc[idx].copy()
        scores = {product: artifact_score(artifact, frame) for product, _, artifact in artifacts}
        if args.write_score_cache:
            save_score_cache(args.write_score_cache, frame, scores)
            return 0

    pt_bins = validation_bins(args.pt_bins, DEFAULT_PT_BINS)
    cent_bins = validation_bins(args.centrality_bins, DEFAULT_CENT_BINS)
    metrics = {product: metric_row(product, frame, score, pt_bins, cent_bins, args.target_signal_efficiency) for product, score in scores.items()}
    rank_rows = []
    for product, row in metrics.items():
        flat = {k: v for k, v in row.items() if not isinstance(v, dict)}
        rank_rows.append(flat)
    rank_rows.sort(key=lambda r: (r.get("wp80_fake") if r.get("wp80_fake") is not None else math.inf, -(r.get("auc") or -math.inf)))
    write_json(outdir / "validation_metrics.json", {"schema": "RJ_AUAU_TIGHT_LOGREG_VALIDATION_METRICS_V1", "read_report": read_report, "products": metrics})
    write_csv(outdir / "validation_rank_table.csv", rank_rows)
    wp_path = write_working_points(outdir, metrics, frame, scores, args.target_signal_efficiency, pt_bins, cent_bins)
    if not args.no_plots:
        make_plots(outdir, frame, scores, metrics, wp_path)
    cache_note = ""
    if args.merge_score_caches:
        (outdir / "score_caches.list").write_text(args.merge_score_caches.read_text())
        cache_note = f"score_caches={read_report.get('score_caches')}\n"
    summary = outdir / "validation_summary.txt"
    best = rank_rows[0] if rank_rows else {}
    summary.write_text(
        "RECOILJETS_AUAU_TIGHT_LOGREG_SIM_VALIDATION_V1\n"
        "status=READY\n"
        f"outdir={outdir}\n"
        f"{cache_note}"
        f"best_product={best.get('product', '')}\n"
        f"best_auc={best.get('auc', '')}\n"
        f"best_wp80_fake={best.get('wp80_fake', '')}\n"
        f"working_points={wp_path}\n"
    )
    print(f"[validateAuAuTightLogReg] wrote metrics: {outdir / 'validation_metrics.json'}")
    print(f"[validateAuAuTightLogReg] wrote ranking table: {outdir / 'validation_rank_table.csv'}")
    print(f"[validateAuAuTightLogReg] wrote summary: {summary}")
    return 0


def derive_from_report(args) -> int:
    report = args.derive_working_points_from_report
    cache_manifest = report / "score_caches.list"
    if not cache_manifest.is_file():
        raise SystemExit(f"Report has no score_caches.list: {report}")
    args.merge_score_caches = cache_manifest
    args.outdir = report
    return run_validation(args)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", type=Path, default=None)
    ap.add_argument("--input", type=Path, nargs="*", default=[])
    ap.add_argument("--manifest", type=Path, default=None)
    ap.add_argument("--tree", default="AuAuPhotonIDTrainingTree")
    ap.add_argument("--model-dir", type=Path, default=None)
    ap.add_argument("--model-registry", type=Path, default=None)
    ap.add_argument("--outdir", type=Path, default=Path("logreg_validation"))
    ap.add_argument("--write-score-cache", type=Path, default=None)
    ap.add_argument("--merge-score-caches", type=Path, default=None)
    ap.add_argument("--derive-working-points-from-report", type=Path, default=None)
    ap.add_argument("--score-max-rows", type=int, default=300000)
    ap.add_argument("--max-files-per-sample", type=int, default=0)
    ap.add_argument("--pt-bins", default=DEFAULT_PT_BINS)
    ap.add_argument("--centrality-bins", default=DEFAULT_CENT_BINS)
    ap.add_argument("--target-signal-efficiency", type=float, default=0.80)
    ap.add_argument("--random-seed", type=int, default=137)
    ap.add_argument("--no-plots", action="store_true")
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    if args.derive_working_points_from_report:
        return derive_from_report(args)
    if not args.merge_score_caches and not args.model_dir:
        raise SystemExit("--model-dir is required unless --merge-score-caches is used")
    return run_validation(args)


if __name__ == "__main__":
    raise SystemExit(main())
