#!/usr/bin/env python3
"""Validate AuAu tight-MLP JSON artifacts on embedded-sim training trees."""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from train_auau_photon_mlp import (  # noqa: E402
    BASE_AND_3X3_FEATURES,
    BASE_FEATURES,
    MODEL_SPECS,
    WIDTH_RATIO_FEATURES,
    add_derived_features,
    auc_score,
    calibration_report,
    expand_required_columns,
    json_ready,
    limit_paths_by_group,
    load_mlp_artifact,
    path_group_key,
    parse_bin_spec,
    predict_mlp_array,
    threshold_for_signal_efficiency,
    write_json,
)


DIAGNOSTIC_FEATURES = []
for _name in BASE_AND_3X3_FEATURES + WIDTH_RATIO_FEATURES + ["centrality", "reco_eiso"]:
    if _name not in DIAGNOSTIC_FEATURES:
        DIAGNOSTIC_FEATURES.append(_name)

EFFICIENCY_TARGETS = [0.50, 0.70, 0.80, 0.90, 0.95]
DEFAULT_PT_BINS = [(15.0, 20.0), (20.0, 25.0), (25.0, 35.0)]
DEFAULT_CENT_BINS = [(0.0, 20.0), (20.0, 50.0), (50.0, 80.0)]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", type=Path, default=None)
    ap.add_argument("--manifest", type=Path, default=None)
    ap.add_argument("--model-dir", type=Path, default=None)
    ap.add_argument("--model-registry", type=Path, default=None)
    ap.add_argument("--outdir", type=Path, default=None)
    ap.add_argument("--tree", default="AuAuPhotonIDTrainingTree")
    ap.add_argument("--write-score-cache", type=Path, default=None)
    ap.add_argument("--merge-score-caches", type=Path, default=None)
    ap.add_argument("--rescore-score-caches", type=Path, default=None)
    ap.add_argument("--sweep-manifest", type=Path, default=None)
    ap.add_argument("--include-cache-scores", dest="include_cache_scores", action="store_true", default=os.environ.get("RJ_AUAU_TIGHT_MLP_INCLUDE_CACHE_SCORES", "1") != "0")
    ap.add_argument("--drop-cache-scores", dest="include_cache_scores", action="store_false")
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--pt-bins", default=os.environ.get("RJ_AUAU_TIGHT_MLP_VALIDATE_PT_BINS", ""))
    ap.add_argument("--centrality-bins", default=os.environ.get("RJ_AUAU_TIGHT_MLP_VALIDATE_CENTRALITY_BINS", ""))
    ap.add_argument("--score-max-rows", type=int, default=int(os.environ.get("RJ_AUAU_TIGHT_MLP_VALIDATE_SCORE_MAX_ROWS", "300000")))
    ap.add_argument("--max-files-per-sample", type=int, default=int(os.environ.get("RJ_AUAU_TIGHT_MLP_VALIDATE_MAX_FILES_PER_SAMPLE", "0")))
    ap.add_argument("--random-seed", type=int, default=int(os.environ.get("RJ_AUAU_TIGHT_MLP_VALIDATE_RANDOM_SEED", "137")))
    ap.add_argument("--progress-every", type=int, default=int(os.environ.get("RJ_AUAU_TIGHT_MLP_VALIDATE_PROGRESS_EVERY", "500")))
    ap.add_argument("--min-finite-fraction", type=float, default=float(os.environ.get("RJ_AUAU_TIGHT_MLP_VALIDATE_MIN_FINITE_FRACTION", "0.95")))
    ap.add_argument("--min-auc", type=float, default=float(os.environ.get("RJ_AUAU_TIGHT_MLP_VALIDATE_MIN_AUC", "0.55")))
    ap.add_argument("--target-signal-efficiency", type=float, default=float(os.environ.get("RJ_AUAU_MLP_TARGET_SIGNAL_EFF", "0.80")))
    ap.add_argument("--derive-working-points-from-report", type=Path, default=None)
    return ap.parse_args()


def discover_manifest(source: Path) -> list[Path]:
    candidates = [source / "manifests" / "training_roots.list", source / "training_roots.list"]
    for manifest in candidates:
        if manifest.is_file():
            paths = [Path(x.strip()) for x in manifest.read_text().splitlines() if x.strip()]
            if paths:
                return paths
    search_root = source / "extraction" if (source / "extraction").is_dir() else source
    paths = sorted(search_root.rglob("*.root"))
    if not paths:
        raise SystemExit(f"No ROOT files found under {search_root}")
    return paths


def read_manifest(path: Path) -> list[Path]:
    if not path.is_file():
        raise SystemExit(f"Manifest does not exist: {path}")
    paths = [Path(x.strip()) for x in path.read_text().splitlines() if x.strip()]
    if not paths:
        raise SystemExit(f"Manifest is empty: {path}")
    return paths


def validation_pt_bins(args) -> list[tuple[float, float]]:
    return parse_bin_spec(args.pt_bins) or list(DEFAULT_PT_BINS)


def validation_centrality_bins(args) -> list[tuple[float, float]]:
    return parse_bin_spec(args.centrality_bins) or list(DEFAULT_CENT_BINS)


def require_imports(needs_uproot: bool = True):
    try:
        import numpy as np  # noqa: F401
    except Exception as exc:
        raise SystemExit(
            "MLP validation needs numpy. Use "
            "RJ_ML_PYTHON=/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python. "
            f"Import failure: {exc}"
        ) from exc
    if needs_uproot:
        try:
            import uproot  # noqa: F401
        except Exception as exc:
            raise SystemExit(
                "MLP validation needs uproot when reading ROOT inputs. Use "
                "RJ_ML_PYTHON=/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python. "
                f"Import failure: {exc}"
            ) from exc


def product_specs_from_registry(model_dir: Path, registry_path: Path | None = None):
    registry = registry_path if registry_path is not None else model_dir / "model_registry.json"
    if not registry.is_file():
        if registry_path is not None:
            raise SystemExit(f"Requested model registry does not exist: {registry}")
        specs = {}
        for product, spec in MODEL_SPECS.items():
            path = model_dir / spec["filename"]
            if path.is_file():
                specs[product] = {
                    "variant": spec["variant"],
                    "features": spec["features"],
                    "filename": spec["filename"],
                    "model_id": product,
                }
        if not specs:
            raise SystemExit(f"No MLP registry or default model JSON files found in {model_dir}")
        return specs
    data = json.loads(registry.read_text())
    specs = {}
    for spec in data.get("models", []):
        if spec.get("report", {}).get("status") == "skipped":
            continue
        product = spec["product"]
        specs[product] = {
            "variant": spec.get("variant", MODEL_SPECS.get(product, {}).get("variant", product)),
            "features": spec["features"],
            "filename": Path(spec["output_json"]).name,
            "model_id": product,
        }
    if not specs:
        raise SystemExit(f"Registry exists but has no trained MLP reports: {registry}")
    return specs


def load_models(model_dir: Path, product_specs):
    models = {}
    missing = []
    for product, spec in product_specs.items():
        path = model_dir / spec["filename"]
        if not path.is_file() or path.stat().st_size <= 0:
            missing.append(str(path))
            continue
        model = load_mlp_artifact(path)
        artifact_features = list(model.get("features", []))
        if artifact_features != list(spec["features"]):
            raise SystemExit(f"Feature order mismatch for {product}: registry={spec['features']} artifact={artifact_features}")
        models[product] = model
    if missing:
        raise SystemExit("Missing required MLP JSON files:\n  " + "\n  ".join(missing))
    return models


def load_sweep_models(path: Path):
    data = json.loads(path.read_text())
    if data.get("schema") != "RJ_AUAU_TIGHT_MLP_SWEEP_V1":
        raise SystemExit(f"Unsupported sweep manifest schema in {path}: {data.get('schema')}")
    models = {}
    product_specs = {}
    for entry in data.get("variants", []):
        name = entry["name"]
        kind = entry.get("kind", "single")
        product_specs[name] = {
            "variant": entry.get("variant", name),
            "model_id": name,
            "features": [],
            "filename": "",
        }
        if kind == "single":
            artifact = Path(entry["artifact"])
            model = load_mlp_artifact(artifact)
            product_specs[name]["features"] = list(model.get("features", []))
            product_specs[name]["filename"] = artifact.name
            models[name] = model
        elif kind == "routed":
            axis = entry["axis"]
            routes = []
            route_features = []
            for route in entry.get("routes", []):
                artifact = Path(route["artifact"])
                model = load_mlp_artifact(artifact)
                features = list(model.get("features", []))
                route_features.extend(features)
                routes.append(
                    {
                        "lo": float(route["lo"]),
                        "hi": float(route["hi"]),
                        "label": route.get("label", f"{route['lo']}:{route['hi']}"),
                        "artifact": str(artifact),
                        "model": model,
                    }
                )
            if not routes:
                raise SystemExit(f"Sweep variant {name} has no routes")
            product_specs[name]["features"] = sorted(set(route_features))
            product_specs[name]["filename"] = Path(routes[0]["artifact"]).name
            models[name] = {
                "schema": "RJ_AUAU_TIGHT_MLP_ROUTED_V1",
                "axis": axis,
                "routes": routes,
            }
        else:
            raise SystemExit(f"Unsupported sweep variant kind for {name}: {kind}")
    if not models:
        raise SystemExit(f"Sweep manifest has no variants: {path}")
    return models, product_specs


def collect_rows(paths: list[Path], tree_name: str, score_max_rows: int, progress_every: int, feature_names: list[str] | None = None):
    import numpy as np
    import uproot

    base_required = ["is_signal", "cluster_Et", "cluster_Eta", "centrality", "reco_eiso"] + BASE_FEATURES
    requested = list(feature_names or [])
    requested += [name for name in DIAGNOSTIC_FEATURES if name not in requested]
    required = sorted(set(base_required + expand_required_columns(requested)))
    arrays = {name: [] for name in required}
    total_entries = 0
    signal_entries = 0
    background_entries = 0
    missing_tree = []
    missing_branches = {}
    selected_rows = 0
    selected_signal_rows = 0
    selected_background_rows = 0
    selected_groups = {}
    per_file_limit = max(1, score_max_rows // max(1, len(paths))) if score_max_rows > 0 else None
    for idx, path in enumerate(paths, 1):
        if progress_every > 0 and (idx == 1 or idx % progress_every == 0 or idx == len(paths)):
            print(f"[validateAuAuTightMLP] reading {idx}/{len(paths)}: {path}", flush=True)
        try:
            with uproot.open(path) as root_file:
                try:
                    tree = root_file[tree_name]
                except Exception:
                    missing_tree.append(str(path))
                    continue
                keys = set(tree.keys())
                missing = [name for name in required if name not in keys]
                if missing:
                    missing_branches[str(path)] = missing
                    continue
                labels = tree["is_signal"].array(library="np")
                total_entries += int(len(labels))
                signal_entries += int((labels == 1).sum())
                background_entries += int((labels == 0).sum())
                read_n = len(labels) if per_file_limit is None else min(len(labels), per_file_limit)
                if read_n <= 0:
                    continue
                chunk = tree.arrays(required, entry_stop=read_n, library="np")
                for name in required:
                    arrays[name].append(chunk[name])
                selected_rows += read_n
                selected_labels = labels[:read_n]
                selected_signal_rows += int((selected_labels == 1).sum())
                selected_background_rows += int((selected_labels == 0).sum())
                group = path_group_key(path)
                group_report = selected_groups.setdefault(group, {"files": 0, "rows": 0, "signal_rows": 0, "background_rows": 0})
                group_report["files"] += 1
                group_report["rows"] += int(read_n)
                group_report["signal_rows"] += int((selected_labels == 1).sum())
                group_report["background_rows"] += int((selected_labels == 0).sum())
        except Exception as exc:
            missing_tree.append(f"{path}: {exc}")
    frame = {name: np.concatenate(parts) if parts else np.array([], dtype="float32") for name, parts in arrays.items()}
    add_derived_features(frame)
    return frame, {
        "files": len(paths),
        "total_entries": total_entries,
        "signal_entries": signal_entries,
        "background_entries": background_entries,
        "scored_entries": int(selected_rows),
        "scored_signal_entries": int(selected_signal_rows),
        "scored_background_entries": int(selected_background_rows),
        "scored_groups": dict(sorted(selected_groups.items())),
        "missing_tree_files": missing_tree,
        "missing_branches": missing_branches,
    }


def score_rows(frame, models):
    import numpy as np

    n = len(frame["is_signal"])
    scores = {}
    for product, model in models.items():
        out = np.full(n, np.nan, dtype="float32")
        if model.get("schema") == "RJ_AUAU_TIGHT_MLP_ROUTED_V1":
            axis_name = model["axis"]
            if axis_name not in frame:
                raise SystemExit(f"Routed MLP {product} needs missing axis column {axis_name}")
            axis = np.asarray(frame[axis_name], dtype="float64")
            for route in model["routes"]:
                route_model = route["model"]
                features = route_model["features"]
                missing = [feature for feature in features if feature not in frame]
                if missing:
                    raise SystemExit(f"Routed MLP {product}/{route['label']} missing cached feature(s): {missing}")
                vals = np.column_stack([frame[f].astype("float64") for f in features])
                route_mask = np.isfinite(axis) & (axis >= route["lo"]) & (axis < route["hi"])
                finite = route_mask & np.isfinite(vals).all(axis=1)
                if finite.any():
                    out[finite] = predict_mlp_array(route_model, vals[finite]).astype("float32")
        else:
            features = model["features"]
            missing = [feature for feature in features if feature not in frame]
            if missing:
                raise SystemExit(f"MLP {product} missing cached feature(s): {missing}")
            vals = np.column_stack([frame[f].astype("float64") for f in features])
            finite = np.isfinite(vals).all(axis=1)
            if finite.any():
                out[finite] = predict_mlp_array(model, vals[finite]).astype("float32")
        scores[product] = out
    return scores


def write_score_cache(path: Path, frame, scores, counts):
    import numpy as np

    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "is_signal": frame["is_signal"].astype("int32"),
        "counts_json": np.array(json.dumps(json_ready(counts)), dtype=object),
    }
    for name in DIAGNOSTIC_FEATURES:
        if name in frame:
            payload[name] = frame[name].astype("float32")
    for product, score in scores.items():
        payload[f"score_{product}"] = score.astype("float32")
    np.savez_compressed(path, **payload)
    print(f"[validateAuAuTightMLP] wrote score cache: {path}", flush=True)


def load_score_caches(cache_manifest: Path):
    import numpy as np

    cache_paths = read_manifest(cache_manifest)
    frame_parts = {name: [] for name in ["is_signal"] + DIAGNOSTIC_FEATURES}
    score_parts = {}
    counts = {
        "files": 0,
        "total_entries": 0,
        "signal_entries": 0,
        "background_entries": 0,
        "scored_entries": 0,
        "scored_signal_entries": 0,
        "scored_background_entries": 0,
        "missing_tree_files": [],
        "missing_branches": {},
        "score_cache_files": len(cache_paths),
    }
    for idx, path in enumerate(cache_paths, 1):
        print(f"[validateAuAuTightMLP] merging score cache {idx}/{len(cache_paths)}: {path}", flush=True)
        with np.load(path, allow_pickle=True) as data:
            for name in frame_parts:
                if name in data:
                    frame_parts[name].append(data[name])
                else:
                    raise SystemExit(f"Missing diagnostic column {name} in score cache: {path}")
            for key in data.files:
                if key.startswith("score_"):
                    score_parts.setdefault(key[len("score_"):], []).append(data[key])
            cached_counts = json.loads(str(data["counts_json"].item())) if "counts_json" in data else {}
            for key in ("files", "total_entries", "signal_entries", "background_entries", "scored_entries", "scored_signal_entries", "scored_background_entries"):
                counts[key] += int(cached_counts.get(key, 0))
            counts["missing_tree_files"].extend(cached_counts.get("missing_tree_files", []))
            counts["missing_branches"].update(cached_counts.get("missing_branches", {}))
    frame = {name: np.concatenate(parts) if parts else np.array([], dtype="float32") for name, parts in frame_parts.items()}
    scores = {product: np.concatenate(parts) if parts else np.array([], dtype="float32") for product, parts in score_parts.items()}
    counts["scored_entries"] = int(len(frame["is_signal"]))
    return frame, scores, counts


def corrcoef_or_nan(x, y):
    import numpy as np

    xa = np.asarray(x, dtype="float64")
    ya = np.asarray(y, dtype="float64")
    mask = np.isfinite(xa) & np.isfinite(ya)
    if mask.sum() < 3:
        return math.nan
    if np.std(xa[mask]) <= 0.0 or np.std(ya[mask]) <= 0.0:
        return math.nan
    return float(np.corrcoef(xa[mask], ya[mask])[0, 1])


def finite_fraction(score):
    import numpy as np

    return float(np.isfinite(score).mean()) if len(score) else math.nan


def binned_metrics(frame, y, score, bins, axis_name):
    import numpy as np

    axis = np.asarray(frame[axis_name], dtype="float64")
    rows = []
    for lo, hi in bins:
        mask = np.isfinite(axis) & (axis >= lo) & (axis < hi)
        if mask.sum() <= 0:
            continue
        rows.append(
            {
                "lo": float(lo),
                "hi": float(hi),
                "entries": int(mask.sum()),
                "auc": auc_score(y[mask], score[mask]),
                "wp80": threshold_for_signal_efficiency(y[mask], score[mask], 0.80),
            }
        )
    return rows


def build_metrics(frame, scores, product_specs, args):
    import numpy as np

    y = frame["is_signal"].astype("int32")
    metrics = {}
    pt_bins = validation_pt_bins(args)
    cent_bins = validation_centrality_bins(args)
    for product, score in scores.items():
        score = np.asarray(score, dtype="float64")
        metrics[product] = {
            "variant": product_specs.get(product, {}).get("variant", product),
            "finite_fraction": finite_fraction(score),
            "auc": auc_score(y, score),
            "thresholds": {
                f"wp{int(round(target * 100)):03d}": threshold_for_signal_efficiency(y, score, target)
                for target in EFFICIENCY_TARGETS
            },
            "calibration": calibration_report(y, score),
            "correlations": {
                "score_vs_cluster_Et": corrcoef_or_nan(score, frame["cluster_Et"]),
                "score_vs_centrality": corrcoef_or_nan(score, frame["centrality"]),
                "score_vs_reco_eiso": corrcoef_or_nan(score, frame["reco_eiso"]),
            },
            "pt_bins": binned_metrics(frame, y, score, pt_bins, "cluster_Et"),
            "centrality_bins": binned_metrics(frame, y, score, cent_bins, "centrality"),
        }
    return metrics


def derive_working_points(frame, scores, product_specs, target: float, pt_bins: list[tuple[float, float]] | None = None):
    import numpy as np

    y = np.asarray(frame["is_signal"], dtype="int32")
    et = np.asarray(frame["cluster_Et"], dtype="float64")
    pt_bins = pt_bins or list(DEFAULT_PT_BINS)
    manifest = {
        "schema": "RJ_AUAU_TIGHT_MLP_WORKING_POINTS_V1",
        "target_signal_efficiency": float(target),
        "entries": [],
        "runtime_entries": [],
    }
    for product, score in scores.items():
        variant = product_specs.get(product, {}).get("variant", product)
        thresholds = []
        edges = [pt_bins[0][0]]
        bin_reports = []
        for lo, hi in pt_bins:
            mask = np.isfinite(et) & (et >= lo) & (et < hi)
            report = threshold_for_signal_efficiency(y[mask], score[mask], target)
            if report is None:
                thresholds.append(math.nan)
            else:
                thresholds.append(float(report["threshold"]))
            edges.append(hi)
            bin_reports.append({"lo": lo, "hi": hi, "report": report})
        if any(not math.isfinite(t) for t in thresholds):
            global_report = threshold_for_signal_efficiency(y, score, target)
            if global_report is None:
                continue
            runtime = f"{variant}|linear|{global_report['threshold']:.10g}|0|{pt_bins[0][0]:.10g}|{pt_bins[-1][1]:.10g}|1"
        else:
            edge_s = ";".join(f"{x:.10g}" for x in edges)
            thr_s = ";".join(f"{x:.10g}" for x in thresholds)
            runtime = f"{variant}|binned|{edge_s}|{thr_s}|{pt_bins[0][0]:.10g}|{pt_bins[-1][1]:.10g}|1"
        manifest["runtime_entries"].append(runtime)
        manifest["entries"].append(
            {
                "product": product,
                "variant": variant,
                "pt_bins": bin_reports,
                "runtime_entry": runtime,
            }
        )
    return manifest


def write_working_point_outputs(outdir: Path, manifest, target_label: str | None = None):
    label = target_label or f"target{int(round(float(manifest['target_signal_efficiency']) * 100)):02d}"
    json_path = outdir / f"mlp_working_points_{label}.json"
    frag_path = outdir / f"mlp_working_points_{label}.yaml"
    json_path.write_text(json.dumps(json_ready(manifest), indent=2, sort_keys=True) + "\n")
    quoted = ", ".join(json.dumps(x) for x in manifest["runtime_entries"])
    frag_path.write_text(
        "# Runtime MLP working-point entries derived from embedded-sim validation.\n"
        f"auau_tight_mlp_working_point_entries: [{quoted}]\n"
    )
    return json_path, frag_path


def write_curve_exports(outdir: Path, frame, scores, metrics):
    import numpy as np

    curves = outdir / "curves"
    curves.mkdir(parents=True, exist_ok=True)
    y = frame["is_signal"].astype("int32")
    for product, score in scores.items():
        rows = []
        finite = np.isfinite(score) & np.isin(y, [0, 1])
        thresholds = np.unique(np.quantile(score[finite], np.linspace(0, 1, 101))) if finite.any() else []
        for thr in thresholds:
            sig = finite & (y == 1)
            bkg = finite & (y == 0)
            rows.append(
                (
                    float(thr),
                    float(np.mean(score[sig] > thr)) if sig.any() else math.nan,
                    float(np.mean(score[bkg] > thr)) if bkg.any() else math.nan,
                )
            )
        path = curves / f"{product}_roc_scan.csv"
        with path.open("w") as out:
            out.write("threshold,signal_efficiency,background_fake_rate\n")
            for row in rows:
                out.write(",".join(str(x) for x in row) + "\n")


def make_plots(outdir: Path, frame, scores, metrics):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception as exc:
        print(f"[validateAuAuTightMLP][WARN] matplotlib unavailable; skipping plots: {exc}", file=sys.stderr)
        return

    plot_dir = outdir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    y = frame["is_signal"].astype("int32")
    for product, score in scores.items():
        finite = np.isfinite(score) & np.isin(y, [0, 1])
        if finite.sum() == 0:
            continue
        fig, ax = plt.subplots(figsize=(6.5, 5.0))
        ax.hist(score[finite & (y == 1)], bins=60, range=(0, 1), density=True, histtype="step", lw=2, label="signal")
        ax.hist(score[finite & (y == 0)], bins=60, range=(0, 1), density=True, histtype="step", lw=2, label="background")
        ax.set_xlabel("MLP score")
        ax.set_ylabel("Normalized entries")
        ax.set_title(product)
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(plot_dir / f"{product}_score_distribution.png", dpi=160)
        plt.close(fig)

        roc_path = outdir / "curves" / f"{product}_roc_scan.csv"
        if roc_path.is_file():
            data = np.genfromtxt(roc_path, delimiter=",", names=True)
            if data.size:
                fig, ax = plt.subplots(figsize=(5.8, 5.2))
                ax.plot(data["background_fake_rate"], data["signal_efficiency"], lw=2)
                ax.set_xlabel("Background fake rate")
                ax.set_ylabel("Signal efficiency")
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 1)
                ax.grid(alpha=0.25)
                ax.set_title(f"{product} ROC scan")
                fig.tight_layout()
                fig.savefig(plot_dir / f"{product}_roc.png", dpi=160)
                plt.close(fig)


def write_summary(path: Path, status: str, args, counts, metrics, notes: list[str]):
    lines = [
        "RECOILJETS_AUAU_TIGHT_MLP_SIM_VALIDATION_V1",
        f"status={status}",
        f"source={args.source}",
        f"model_dir={args.model_dir}",
        f"report_dir={path.parent}",
        f"scored_entries={counts.get('scored_entries', 0)}",
        f"scored_signal_entries={counts.get('scored_signal_entries', 0)}",
        f"scored_background_entries={counts.get('scored_background_entries', 0)}",
        f"signal_entries={counts.get('signal_entries', 0)}",
        f"background_entries={counts.get('background_entries', 0)}",
    ]
    for product, m in metrics.items():
        wp80 = m.get("thresholds", {}).get("wp080") or {}
        lines.extend(
            [
                f"product.{product}.auc={m.get('auc')}",
                f"product.{product}.finite_fraction={m.get('finite_fraction')}",
                f"product.{product}.wp080_threshold={wp80.get('threshold')}",
                f"product.{product}.wp080_fake_rate={wp80.get('background_fake_rate')}",
                f"product.{product}.score_vs_eiso_corr={m.get('correlations', {}).get('score_vs_reco_eiso')}",
            ]
        )
    for note in notes:
        lines.append(f"note={note}")
    path.write_text("\n".join(lines) + "\n")


def write_ranking_table(outdir: Path, metrics):
    rows = []
    for product, m in metrics.items():
        wp80 = m.get("thresholds", {}).get("wp080") or {}
        row = {
            "product": product,
            "variant": m.get("variant", product),
            "auc": m.get("auc"),
            "wp080_fake_rate": wp80.get("background_fake_rate"),
            "finite_fraction": m.get("finite_fraction"),
            "ece": (m.get("calibration") or {}).get("ece"),
            "score_vs_eiso_corr": (m.get("correlations") or {}).get("score_vs_reco_eiso"),
        }
        for pt in m.get("pt_bins", []):
            label = f"pt_{pt['lo']:g}_{pt['hi']:g}".replace(".", "p")
            row[f"{label}_auc"] = pt.get("auc")
            row[f"{label}_wp080_fake_rate"] = (pt.get("wp80") or {}).get("background_fake_rate")
        rows.append(row)
    keys = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    path = outdir / "validation_rank_table.csv"
    with path.open("w") as out:
        out.write(",".join(keys) + "\n")
        for row in rows:
            out.write(",".join("" if row.get(k) is None else str(row.get(k)) for k in keys) + "\n")
    md = outdir / "validation_rank_table.md"
    metric_keys = [k for k in keys if k not in ("variant",)]
    with md.open("w") as out:
        out.write("| " + " | ".join(metric_keys) + " |\n")
        out.write("| " + " | ".join(["---"] * len(metric_keys)) + " |\n")
        for row in sorted(rows, key=lambda r: (r.get("wp080_fake_rate") is None, r.get("wp080_fake_rate") or math.inf, -(r.get("auc") or -math.inf))):
            out.write("| " + " | ".join("" if row.get(k) is None else str(row.get(k)) for k in metric_keys) + " |\n")
    return path, md


def main() -> int:
    args = parse_args()
    require_imports(needs_uproot=not (args.derive_working_points_from_report or args.merge_score_caches or args.rescore_score_caches))
    if args.derive_working_points_from_report:
        report_dir = args.derive_working_points_from_report
        metrics_path = report_dir / "validation_metrics.json"
        scores_manifest = report_dir / "score_caches.list"
        if not metrics_path.is_file():
            raise SystemExit(f"Missing validation metrics: {metrics_path}")
        if not scores_manifest.is_file():
            raise SystemExit(f"Missing score cache manifest for WP derivation: {scores_manifest}")
        frame, scores, counts = load_score_caches(scores_manifest)
        product_specs = {p: {"variant": (json.loads(metrics_path.read_text()).get("products", {}).get(p, {}).get("variant", p))} for p in scores}
        manifest = derive_working_points(frame, scores, product_specs, args.target_signal_efficiency, validation_pt_bins(args))
        json_path, frag_path = write_working_point_outputs(report_dir, manifest)
        print(f"[validateAuAuTightMLP] wrote {json_path}")
        print(f"[validateAuAuTightMLP] wrote {frag_path}")
        return 0

    if args.outdir is not None:
        outdir = args.outdir
    elif args.source is not None:
        outdir = args.source / "reports" / "mlp_model_validation"
    else:
        outdir = Path("mlp_model_validation")
    outdir.mkdir(parents=True, exist_ok=True)

    if args.rescore_score_caches is not None:
        frame, cached_scores, counts = load_score_caches(args.rescore_score_caches)
        if args.sweep_manifest is not None:
            models, product_specs = load_sweep_models(args.sweep_manifest)
        else:
            if args.model_dir is None:
                raise SystemExit("--model-dir or --sweep-manifest is required with --rescore-score-caches")
            product_specs = product_specs_from_registry(args.model_dir, args.model_registry)
            models = load_models(args.model_dir, product_specs)
        scores = score_rows(frame, models)
        if args.include_cache_scores:
            for product, score in cached_scores.items():
                key = product if product not in scores else f"cache_{product}"
                scores[key] = score
                product_specs[key] = {"variant": MODEL_SPECS.get(product, {}).get("variant", product)}
        (outdir / "score_caches_source.list").write_text(args.rescore_score_caches.read_text())
    elif args.merge_score_caches is not None:
        frame, scores, counts = load_score_caches(args.merge_score_caches)
        product_specs = {product: {"variant": MODEL_SPECS.get(product, {}).get("variant", product)} for product in scores}
    else:
        if args.sweep_manifest is not None:
            models, product_specs = load_sweep_models(args.sweep_manifest)
        else:
            if args.model_dir is None:
                raise SystemExit("--model-dir is required unless --merge-score-caches, --rescore-score-caches, or --sweep-manifest is used")
            product_specs = product_specs_from_registry(args.model_dir, args.model_registry)
            models = load_models(args.model_dir, product_specs)
        paths = read_manifest(args.manifest) if args.manifest else discover_manifest(args.source)
        paths, manifest_report = limit_paths_by_group(paths, args.max_files_per_sample, args.random_seed)
        manifest_report.update(
            {
                "schema": "RJ_AUAU_TIGHT_MLP_VALIDATION_MANIFEST_SUMMARY_V1",
                "source": str(args.source) if args.source else None,
                "tree": args.tree,
                "score_max_rows": int(args.score_max_rows),
                "random_seed": int(args.random_seed),
            }
        )
        write_json(outdir / "validation_manifest_summary.json", manifest_report)
        print(
            "[validateAuAuTightMLP] manifest "
            f"files_before={manifest_report['input_files_before']['total_files']} "
            f"files_after={manifest_report['input_files_after']['total_files']} "
            f"groups={manifest_report['input_files_after']['groups']} "
            f"score_max_rows={args.score_max_rows}",
            flush=True,
        )
        feature_names = sorted(set(name for spec in product_specs.values() for name in spec["features"]))
        frame, counts = collect_rows(paths, args.tree, args.score_max_rows, args.progress_every, feature_names)
        print(
            "[validateAuAuTightMLP] collected "
            f"scored_entries={counts['scored_entries']} "
            f"signal={counts['scored_signal_entries']} "
            f"background={counts['scored_background_entries']}",
            flush=True,
        )
        scores = score_rows(frame, models)

    if args.write_score_cache is not None:
        write_score_cache(args.write_score_cache, frame, scores, counts)
        if args.no_plots:
            return 0
        (outdir / "score_caches.list").write_text(str(args.write_score_cache) + "\n")
    elif args.merge_score_caches is not None:
        (outdir / "score_caches.list").write_text(args.merge_score_caches.read_text())
    else:
        full_cache = outdir / "score_caches" / "score_cache_full.npz"
        write_score_cache(full_cache, frame, scores, counts)
        (outdir / "score_caches.list").write_text(str(full_cache) + "\n")

    metrics = build_metrics(frame, scores, product_specs, args)
    wp_manifest = derive_working_points(frame, scores, product_specs, args.target_signal_efficiency, validation_pt_bins(args))
    wp_json, wp_frag = write_working_point_outputs(outdir, wp_manifest)
    write_curve_exports(outdir, frame, scores, metrics)
    rank_csv, rank_md = write_ranking_table(outdir, metrics)
    if not args.no_plots:
        make_plots(outdir, frame, scores, metrics)

    metrics_path = outdir / "validation_metrics.json"
    metrics_path.write_text(json.dumps(json_ready({"schema": "RJ_AUAU_TIGHT_MLP_VALIDATION_METRICS_V1", "products": metrics}), indent=2, sort_keys=True) + "\n")
    notes = []
    status = "READY"
    for product, m in metrics.items():
        if not math.isfinite(m.get("finite_fraction", math.nan)) or m["finite_fraction"] < args.min_finite_fraction:
            status = "CHECK"
            notes.append(f"{product}: finite fraction below threshold")
        if not math.isfinite(m.get("auc", math.nan)) or m["auc"] < args.min_auc:
            status = "CHECK"
            notes.append(f"{product}: AUC below threshold")
    summary = outdir / "validation_summary.txt"
    write_summary(summary, status, args, counts, metrics, notes)
    print(f"[validateAuAuTightMLP] wrote metrics: {metrics_path}")
    print(f"[validateAuAuTightMLP] wrote ranking table: {rank_csv}")
    print(f"[validateAuAuTightMLP] wrote ranking table markdown: {rank_md}")
    print(f"[validateAuAuTightMLP] wrote working points: {wp_json}")
    print(f"[validateAuAuTightMLP] wrote working-point YAML fragment: {wp_frag}")
    print(f"[validateAuAuTightMLP] wrote summary: {summary}")
    return 0 if status == "READY" else 3


if __name__ == "__main__":
    raise SystemExit(main())
