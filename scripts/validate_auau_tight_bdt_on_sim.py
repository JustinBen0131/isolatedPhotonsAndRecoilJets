#!/usr/bin/env python3
"""Validate trained AuAu tight-BDT TMVA models on embedded-sim training trees.

This is a sidecar QA step: it reads AuAuPhotonIDTrainingTree files produced by
the extract-only workflow, scores rows with the same TMVA/RBDT ROOT files that
production will consume, and writes quick PPG12-style separation diagnostics.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path


BASE_FEATURES = [
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

PRODUCTS = {
    "centINDcontrol": {
        "features": BASE_FEATURES,
        "models": [("all", "auau_tight_bdt_centINDcontrol_allCent_tmva.root")],
    },
    "centAsFeat": {
        "features": BASE_FEATURES + ["centrality"],
        "models": [("all", "auau_tight_bdt_centAsFeat_allCent_tmva.root")],
    },
    "centDepBDTs": {
        "features": BASE_FEATURES,
        "models": [
            ((0.0, 20.0), "auau_tight_bdt_centDepBDTs_cent_000_020_tmva.root"),
            ((20.0, 50.0), "auau_tight_bdt_centDepBDTs_cent_020_050_tmva.root"),
            ((50.0, 80.0), "auau_tight_bdt_centDepBDTs_cent_050_080_tmva.root"),
        ],
    },
}

PLOT_LABELS = {
    "centINDcontrol": "Cluster-shape only",
    "centAsFeat": "Cluster-shape + centrality",
    "centDepBDTs": "Centrality-specific models",
}

PLOT_COLORS = {
    "centINDcontrol": "#1f77b4",
    "centAsFeat": "#ff7f0e",
    "centDepBDTs": "#2ca02c",
}

PT_BINS = [(6, 10), (10, 15), (15, 20), (20, 25), (25, 35)]
CENT_BINS = [(0, 20), (20, 50), (50, 80)]
DIAGNOSTIC_FEATURES = []
for _name in BASE_FEATURES + ["centrality", "reco_eiso"]:
    if _name not in DIAGNOSTIC_FEATURES:
        DIAGNOSTIC_FEATURES.append(_name)
EFFICIENCY_TARGETS = [0.50, 0.70, 0.80, 0.90, 0.95]
FAKE_RATE_TARGETS = [0.01, 0.02, 0.05, 0.10, 0.20]


def range_key(lo, hi) -> str:
    return f"{lo:g}_{hi:g}".replace(".", "p")


def cent_label(lo, hi) -> str:
    return f"{lo:g}-{hi:g}%"


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", type=Path, required=True)
    ap.add_argument("--model-dir", type=Path, required=True)
    ap.add_argument("--manifest", type=Path, default=None, help="Optional ROOT-file manifest for sharded validation.")
    ap.add_argument("--outdir", type=Path, default=None)
    ap.add_argument("--tree", default="AuAuPhotonIDTrainingTree")
    ap.add_argument("--write-score-cache", type=Path, default=None, help="Write scored rows to this .npz cache.")
    ap.add_argument("--merge-score-caches", type=Path, default=None, help="Merge a list of .npz score caches and make final plots.")
    ap.add_argument("--no-plots", action="store_true", help="Skip plot/curve writing; useful for Condor scoring shards.")
    ap.add_argument("--score-max-rows", type=int, default=int(os.environ.get("RJ_AUAU_TIGHT_BDT_VALIDATE_SCORE_MAX_ROWS", "300000")))
    ap.add_argument("--progress-every", type=int, default=int(os.environ.get("RJ_AUAU_TIGHT_BDT_VALIDATE_PROGRESS_EVERY", "500")))
    ap.add_argument("--min-finite-fraction", type=float, default=float(os.environ.get("RJ_AUAU_TIGHT_BDT_VALIDATE_MIN_FINITE_FRACTION", "0.95")))
    ap.add_argument("--min-auc", type=float, default=float(os.environ.get("RJ_AUAU_TIGHT_BDT_VALIDATE_MIN_AUC", "0.55")))
    return ap.parse_args()


def discover_manifest(source: Path) -> list[Path]:
    candidates = [
        source / "manifests" / "training_roots.list",
        source / "training_roots.list",
    ]
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


def require_imports():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # noqa: F401
        import numpy as np  # noqa: F401
        import uproot  # noqa: F401
        from sklearn.metrics import roc_auc_score, roc_curve  # noqa: F401
        import ROOT  # noqa: F401
    except Exception as exc:
        raise SystemExit(
            "Validation needs numpy, uproot, sklearn, matplotlib, and PyROOT. "
            "Use RJ_ML_PYTHON=/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python. "
            f"Import failure: {exc}"
        ) from exc


def load_models(model_dir: Path):
    import ROOT

    loaded = {}
    missing = []
    for product, spec in PRODUCTS.items():
        entries = []
        for selector, filename in spec["models"]:
            path = model_dir / filename
            if not path.is_file() or path.stat().st_size <= 0:
                missing.append(str(path))
                continue
            try:
                entries.append((selector, ROOT.TMVA.Experimental.RBDT("myBDT", str(path))))
            except Exception as exc:
                raise SystemExit(f"Failed to load TMVA model {path}: {exc}") from exc
        loaded[product] = entries
    if missing:
        raise SystemExit("Missing required TMVA model files:\n  " + "\n  ".join(missing))
    return loaded


def cent_selector_mask(cent, selector):
    import numpy as np

    if selector == "all":
        return np.ones(len(cent), dtype=bool)
    lo, hi = selector
    return np.isfinite(cent) & (cent >= lo) & (cent < hi)


def score_rows(frame, models):
    import numpy as np
    import ROOT

    scores = {}
    for product, spec in PRODUCTS.items():
        features = spec["features"]
        out = np.full(len(frame["is_signal"]), np.nan, dtype="float32")
        for selector, model in models[product]:
            mask = cent_selector_mask(frame["centrality"], selector)
            idx = np.flatnonzero(mask)
            if len(idx) == 0:
                continue
            vals = np.column_stack([frame[f][idx].astype("float32") for f in features])
            finite = np.isfinite(vals).all(axis=1)
            idx = idx[finite]
            vals = vals[finite]
            for local_i, row in zip(idx, vals):
                try:
                    xvec = ROOT.std.vector("float")()
                    for x in row:
                        xvec.push_back(float(x))
                    y = model.Compute(xvec)
                    if len(y) and math.isfinite(float(y[0])):
                        out[local_i] = float(y[0])
                except Exception:
                    out[local_i] = np.nan
        scores[product] = out
    return scores


def collect_rows(paths: list[Path], tree_name: str, score_max_rows: int, progress_every: int):
    import numpy as np
    import uproot

    required = sorted(
        set(["is_signal", "cluster_Et", "cluster_Eta", "centrality", "reco_eiso"] + BASE_FEATURES)
    )
    arrays = {name: [] for name in required}
    total_entries = 0
    signal_entries = 0
    background_entries = 0
    missing_tree = []
    missing_branches = {}
    selected_rows = 0
    per_file_score_limit = None
    if score_max_rows > 0:
        per_file_score_limit = max(1, score_max_rows // max(1, len(paths)))

    for idx, path in enumerate(paths, 1):
        if progress_every > 0 and (idx == 1 or idx % progress_every == 0 or idx == len(paths)):
            print(f"[validateAuAuTightBDT] reading {idx}/{len(paths)}: {path}", flush=True)
        try:
            with uproot.open(path) as f:
                try:
                    tree = f[tree_name]
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
                read_n = len(labels)
                if per_file_score_limit is not None:
                    read_n = min(read_n, per_file_score_limit)
                if read_n <= 0:
                    continue
                chunk = tree.arrays(required, entry_stop=read_n, library="np")
                for name in required:
                    arrays[name].append(chunk[name])
                selected_rows += read_n
        except Exception as exc:
            missing_tree.append(f"{path}: {exc}")

    frame = {}
    for name, parts in arrays.items():
        frame[name] = np.concatenate(parts) if parts else np.array([], dtype="float32")
    return frame, {
        "files": len(paths),
        "total_entries": total_entries,
        "signal_entries": signal_entries,
        "background_entries": background_entries,
        "scored_entries": int(len(frame["is_signal"])),
        "missing_tree_files": missing_tree,
        "missing_branches": missing_branches,
    }


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
    print(f"[validateAuAuTightBDT] wrote score cache: {path}", flush=True)


def load_score_caches(cache_manifest: Path):
    import numpy as np

    cache_paths = read_manifest(cache_manifest)
    frame_parts = {name: [] for name in ["is_signal"] + DIAGNOSTIC_FEATURES}
    score_parts = {product: [] for product in PRODUCTS}
    counts = {
        "files": 0,
        "total_entries": 0,
        "signal_entries": 0,
        "background_entries": 0,
        "scored_entries": 0,
        "missing_tree_files": [],
        "missing_branches": {},
        "score_cache_files": len(cache_paths),
    }
    for idx, path in enumerate(cache_paths, 1):
        print(f"[validateAuAuTightBDT] merging score cache {idx}/{len(cache_paths)}: {path}", flush=True)
        with np.load(path, allow_pickle=True) as data:
            for name in frame_parts:
                if name in data:
                    frame_parts[name].append(data[name])
                else:
                    raise SystemExit(
                        f"Missing diagnostic column {name} in score cache: {path}. "
                        "Regenerate validation caches with the updated validator."
                    )
            for product in PRODUCTS:
                key = f"score_{product}"
                if key not in data:
                    raise SystemExit(f"Missing {key} in score cache: {path}")
                score_parts[product].append(data[key])
            cached_counts = json.loads(str(data["counts_json"].item())) if "counts_json" in data else {}
            counts["files"] += int(cached_counts.get("files", 0))
            counts["total_entries"] += int(cached_counts.get("total_entries", 0))
            counts["signal_entries"] += int(cached_counts.get("signal_entries", 0))
            counts["background_entries"] += int(cached_counts.get("background_entries", 0))
            counts["scored_entries"] += int(cached_counts.get("scored_entries", 0))
            counts["missing_tree_files"].extend(cached_counts.get("missing_tree_files", []))
            counts["missing_branches"].update(cached_counts.get("missing_branches", {}))

    frame = {
        name: np.concatenate(parts) if parts else np.array([], dtype="float32")
        for name, parts in frame_parts.items()
    }
    scores = {
        product: np.concatenate(parts) if parts else np.array([], dtype="float32")
        for product, parts in score_parts.items()
    }
    counts["scored_entries"] = int(len(frame["is_signal"]))
    return frame, scores, counts


def auc_for(y, score, weight=None):
    import numpy as np
    from sklearn.metrics import roc_auc_score, roc_curve

    mask = np.isfinite(score) & np.isin(y, [0, 1])
    if mask.sum() < 2 or len(np.unique(y[mask])) < 2:
        return math.nan, None
    auc = float(roc_auc_score(y[mask], score[mask], sample_weight=None if weight is None else weight[mask]))
    fpr, tpr, thresholds = roc_curve(y[mask], score[mask], sample_weight=None if weight is None else weight[mask])
    return auc, (fpr, tpr, thresholds)


def finite_values(values):
    import numpy as np

    arr = np.asarray(values, dtype="float64")
    return arr[np.isfinite(arr)]


def summary_stats(values):
    import numpy as np

    arr = finite_values(values)
    if len(arr) == 0:
        return {
            "entries": 0,
            "mean": math.nan,
            "std": math.nan,
            "min": math.nan,
            "q01": math.nan,
            "q05": math.nan,
            "q10": math.nan,
            "q25": math.nan,
            "median": math.nan,
            "q75": math.nan,
            "q90": math.nan,
            "q95": math.nan,
            "q99": math.nan,
            "max": math.nan,
        }
    qs = np.percentile(arr, [1, 5, 10, 25, 50, 75, 90, 95, 99])
    return {
        "entries": int(len(arr)),
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr)),
        "min": float(np.min(arr)),
        "q01": float(qs[0]),
        "q05": float(qs[1]),
        "q10": float(qs[2]),
        "q25": float(qs[3]),
        "median": float(qs[4]),
        "q75": float(qs[5]),
        "q90": float(qs[6]),
        "q95": float(qs[7]),
        "q99": float(qs[8]),
        "max": float(np.max(arr)),
    }


def corrcoef_or_nan(x, y):
    import numpy as np

    xa = np.asarray(x, dtype="float64")
    ya = np.asarray(y, dtype="float64")
    mask = np.isfinite(xa) & np.isfinite(ya)
    if mask.sum() < 3:
        return math.nan
    if np.std(xa[mask]) <= 0 or np.std(ya[mask]) <= 0:
        return math.nan
    return float(np.corrcoef(xa[mask], ya[mask])[0, 1])


def threshold_report(y, score):
    import numpy as np
    from sklearn.metrics import roc_curve

    mask = np.isfinite(score) & np.isin(y, [0, 1])
    if mask.sum() < 2 or len(np.unique(y[mask])) < 2:
        return {"by_signal_efficiency": [], "by_background_fake_rate": []}
    fpr, tpr, thresholds = roc_curve(y[mask], score[mask])
    thresholds = np.asarray(thresholds, dtype="float64")
    finite_threshold = np.isfinite(thresholds)

    by_eff = []
    for target in EFFICIENCY_TARGETS:
        ok = np.flatnonzero(finite_threshold & (tpr >= target))
        if len(ok) == 0:
            continue
        idx = ok[np.argmin(fpr[ok])]
        by_eff.append(
            {
                "target_signal_efficiency": target,
                "threshold": float(thresholds[idx]),
                "signal_efficiency": float(tpr[idx]),
                "background_fake_rate": float(fpr[idx]),
                "background_rejection": float(1.0 - fpr[idx]),
            }
        )

    by_fake = []
    for target in FAKE_RATE_TARGETS:
        ok = np.flatnonzero(finite_threshold & (fpr <= target))
        if len(ok) == 0:
            continue
        idx = ok[np.argmax(tpr[ok])]
        by_fake.append(
            {
                "target_background_fake_rate": target,
                "threshold": float(thresholds[idx]),
                "signal_efficiency": float(tpr[idx]),
                "background_fake_rate": float(fpr[idx]),
                "background_rejection": float(1.0 - fpr[idx]),
            }
        )
    return {"by_signal_efficiency": by_eff, "by_background_fake_rate": by_fake}


def build_deep_diagnostics(frame, scores):
    import numpy as np

    y = frame["is_signal"].astype("int32")
    cent = frame["centrality"].astype("float64")
    et = frame["cluster_Et"].astype("float64")
    diagnostics = {
        "description": (
            "Offline validation diagnostics for questions about centrality, pT, "
            "threshold choices, feature distributions, and score correlations."
        ),
        "pt_bins": [{"key": range_key(lo, hi), "low": lo, "high": hi} for lo, hi in PT_BINS],
        "centrality_bins": [{"key": range_key(lo, hi), "low": lo, "high": hi} for lo, hi in CENT_BINS],
        "products": {},
        "features": {},
    }

    for product, score in scores.items():
        eligible = np.isin(y, [0, 1]) & np.isfinite(score)
        if product == "centDepBDTs":
            eligible &= np.isfinite(cent) & (cent >= 0.0) & (cent < 80.0)
        product_diag = {
            "thresholds_inclusive": threshold_report(y[eligible], score[eligible]),
            "score_summary": {
                "signal": summary_stats(score[eligible & (y == 1)]),
                "background": summary_stats(score[eligible & (y == 0)]),
            },
            "auc_by_pt": {},
            "auc_by_centrality": {},
            "auc_by_centrality_and_pt": {},
            "score_correlations": {},
        }
        for name in DIAGNOSTIC_FEATURES:
            if name in frame:
                product_diag["score_correlations"][name] = corrcoef_or_nan(score[eligible], frame[name][eligible])

        for lo, hi in PT_BINS:
            key = range_key(lo, hi)
            mask = eligible & np.isfinite(et) & (et >= lo) & (et < hi)
            a, _ = auc_for(y[mask], score[mask])
            product_diag["auc_by_pt"][key] = {
                "auc": a,
                "entries": int(mask.sum()),
                "signal_entries": int((mask & (y == 1)).sum()),
                "background_entries": int((mask & (y == 0)).sum()),
                "thresholds": threshold_report(y[mask], score[mask]),
                "signal_score": summary_stats(score[mask & (y == 1)]),
                "background_score": summary_stats(score[mask & (y == 0)]),
            }

        for clo, chi in CENT_BINS:
            ckey = range_key(clo, chi)
            cmask = eligible & np.isfinite(cent) & (cent >= clo) & (cent < chi)
            a, _ = auc_for(y[cmask], score[cmask])
            product_diag["auc_by_centrality"][ckey] = {
                "auc": a,
                "entries": int(cmask.sum()),
                "signal_entries": int((cmask & (y == 1)).sum()),
                "background_entries": int((cmask & (y == 0)).sum()),
                "thresholds": threshold_report(y[cmask], score[cmask]),
                "signal_score": summary_stats(score[cmask & (y == 1)]),
                "background_score": summary_stats(score[cmask & (y == 0)]),
            }
            product_diag["auc_by_centrality_and_pt"][ckey] = {}
            for plo, phi in PT_BINS:
                pkey = range_key(plo, phi)
                mask = cmask & np.isfinite(et) & (et >= plo) & (et < phi)
                a2, _ = auc_for(y[mask], score[mask])
                product_diag["auc_by_centrality_and_pt"][ckey][pkey] = {
                    "auc": a2,
                    "entries": int(mask.sum()),
                    "signal_entries": int((mask & (y == 1)).sum()),
                    "background_entries": int((mask & (y == 0)).sum()),
                }
        diagnostics["products"][product] = product_diag

    for name in DIAGNOSTIC_FEATURES:
        if name not in frame:
            continue
        values = frame[name]
        feature_diag = {
            "inclusive": {
                "signal": summary_stats(values[y == 1]),
                "background": summary_stats(values[y == 0]),
            },
            "by_centrality": {},
            "by_pt": {},
        }
        for lo, hi in CENT_BINS:
            key = range_key(lo, hi)
            mask = np.isfinite(cent) & (cent >= lo) & (cent < hi)
            feature_diag["by_centrality"][key] = {
                "signal": summary_stats(values[mask & (y == 1)]),
                "background": summary_stats(values[mask & (y == 0)]),
            }
        for lo, hi in PT_BINS:
            key = range_key(lo, hi)
            mask = np.isfinite(et) & (et >= lo) & (et < hi)
            feature_diag["by_pt"][key] = {
                "signal": summary_stats(values[mask & (y == 1)]),
                "background": summary_stats(values[mask & (y == 0)]),
            }
        diagnostics["features"][name] = feature_diag

    return diagnostics


def write_diagnostic_tables(outdir: Path, diagnostics):
    import csv

    auc_path = outdir / "validation_auc_table.csv"
    with auc_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["product", "centrality_bin", "pt_bin", "auc", "entries", "signal_entries", "background_entries"])
        for product, product_diag in diagnostics["products"].items():
            for ckey, cent_diag in product_diag["auc_by_centrality_and_pt"].items():
                for pkey, cell in cent_diag.items():
                    writer.writerow([product, ckey, pkey, cell["auc"], cell["entries"], cell["signal_entries"], cell["background_entries"]])

    threshold_path = outdir / "validation_threshold_table.csv"
    with threshold_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["product", "scope", "target_type", "target", "threshold", "signal_efficiency", "background_fake_rate", "background_rejection"])
        for product, product_diag in diagnostics["products"].items():
            scopes = {"inclusive": product_diag["thresholds_inclusive"]}
            scopes.update({f"cent_{k}": v["thresholds"] for k, v in product_diag["auc_by_centrality"].items()})
            scopes.update({f"pt_{k}": v["thresholds"] for k, v in product_diag["auc_by_pt"].items()})
            for scope, report in scopes.items():
                for row in report["by_signal_efficiency"]:
                    writer.writerow([product, scope, "signal_efficiency", row["target_signal_efficiency"], row["threshold"], row["signal_efficiency"], row["background_fake_rate"], row["background_rejection"]])
                for row in report["by_background_fake_rate"]:
                    writer.writerow([product, scope, "background_fake_rate", row["target_background_fake_rate"], row["threshold"], row["signal_efficiency"], row["background_fake_rate"], row["background_rejection"]])

    feature_path = outdir / "validation_feature_summary.csv"
    with feature_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["feature", "scope", "class", "entries", "mean", "std", "q10", "median", "q90"])
        for feature, feature_diag in diagnostics["features"].items():
            scopes = {"inclusive": feature_diag["inclusive"]}
            scopes.update({f"cent_{k}": v for k, v in feature_diag["by_centrality"].items()})
            scopes.update({f"pt_{k}": v for k, v in feature_diag["by_pt"].items()})
            for scope, classes in scopes.items():
                for class_name, stats in classes.items():
                    writer.writerow([feature, scope, class_name, stats["entries"], stats["mean"], stats["std"], stats["q10"], stats["median"], stats["q90"]])


def write_curve_exports(outdir: Path, frame, scores, metrics):
    """Write the real curve points needed for later ROOT/sPHENIX plotting."""
    import numpy as np
    import ROOT

    y = frame["is_signal"].astype("int32")
    score_bins = np.linspace(0.0, 1.0, 51)
    roc_payload = {
        "description": "ROC points from AuAu tight-BDT validation; no smoothing or refit applied.",
        "x_axis_background_fake_rate": "fraction of background candidates accepted at this score threshold",
        "y_axis_signal_efficiency": "fraction of signal photons accepted at this score threshold",
        "products": {},
    }
    hist_payload = {
        "description": "Binned BDT score shapes from the same validation rows used for ROC/AUC.",
        "bin_edges": score_bins.tolist(),
        "products": {},
    }

    root_path = outdir / "validation_curves.root"
    fout = ROOT.TFile(str(root_path), "RECREATE")
    for product, score in scores.items():
        label = PLOT_LABELS.get(product, product)
        roc = metrics["products"].get(product, {}).get("_roc")
        if roc is not None:
            fpr, tpr, thresholds = roc
            rejection = 1.0 - fpr
            roc_payload["products"][product] = {
                "label": label,
                "color": PLOT_COLORS.get(product, "black"),
                "auc": metrics["products"][product]["auc_inclusive"],
                "background_fake_rate": np.asarray(fpr, dtype=float).tolist(),
                "background_rejection": np.asarray(rejection, dtype=float).tolist(),
                "signal_efficiency": np.asarray(tpr, dtype=float).tolist(),
                "threshold": [
                    None if (not np.isfinite(x)) else float(x)
                    for x in np.asarray(thresholds, dtype=float)
                ],
            }

            g_fake = ROOT.TGraph(len(fpr), np.asarray(fpr, dtype="float64"), np.asarray(tpr, dtype="float64"))
            g_fake.SetName(f"gROC_{product}_signalEfficiency_vs_backgroundFakeRate")
            g_fake.SetTitle(f"{label};Background fake rate;Signal efficiency")
            g_fake.Write()

            g_rej = ROOT.TGraph(len(rejection), np.asarray(rejection, dtype="float64"), np.asarray(tpr, dtype="float64"))
            g_rej.SetName(f"gROC_{product}_signalEfficiency_vs_backgroundRejection")
            g_rej.SetTitle(f"{label};Background rejection;Signal efficiency")
            g_rej.Write()

        finite = np.isfinite(score)
        product_hist = {"label": label, "color": PLOT_COLORS.get(product, "black"), "by_centrality": {}}
        for class_name, class_value in (("background", 0), ("signal", 1)):
            values = score[finite & (y == class_value)]
            counts, _ = np.histogram(values, bins=score_bins)
            width = np.diff(score_bins)
            density = counts.astype(float)
            if density.sum() > 0:
                density = density / density.sum() / width
            product_hist[class_name] = {
                "counts": counts.astype(int).tolist(),
                "density": density.astype(float).tolist(),
                "entries": int(len(values)),
            }

            h = ROOT.TH1D(
                f"hScore_{product}_{class_name}",
                f"{label} {class_name};BDT score;Candidates",
                len(score_bins) - 1,
                score_bins.astype("float64"),
            )
            h.SetDirectory(fout)
            for val in values:
                h.Fill(float(val))
            h.Write()

            h_norm = h.Clone(f"hScore_{product}_{class_name}_unitArea")
            if h_norm.Integral() > 0:
                h_norm.Scale(1.0 / h_norm.Integral(), "width")
            h_norm.SetTitle(f"{label} {class_name};BDT score;Unit-normalized candidates")
            h_norm.Write()

        roc_by_cent = {}
        cent = frame["centrality"].astype("float64")
        for lo, hi in CENT_BINS:
            key = range_key(lo, hi)
            cent_mask = np.isfinite(cent) & (cent >= lo) & (cent < hi)
            auc_cent, roc_cent = auc_for(y[cent_mask], score[cent_mask])
            if roc_cent is not None:
                fpr, tpr, thresholds = roc_cent
                rejection = 1.0 - fpr
                roc_by_cent[key] = {
                    "label": cent_label(lo, hi),
                    "auc": auc_cent,
                    "background_fake_rate": np.asarray(fpr, dtype=float).tolist(),
                    "background_rejection": np.asarray(rejection, dtype=float).tolist(),
                    "signal_efficiency": np.asarray(tpr, dtype=float).tolist(),
                    "threshold": [
                        None if (not np.isfinite(x)) else float(x)
                        for x in np.asarray(thresholds, dtype=float)
                    ],
                }

                g_fake = ROOT.TGraph(len(fpr), np.asarray(fpr, dtype="float64"), np.asarray(tpr, dtype="float64"))
                g_fake.SetName(f"gROC_{product}_cent_{key}_signalEfficiency_vs_backgroundFakeRate")
                g_fake.SetTitle(f"{label} {cent_label(lo, hi)};Background fake rate;Signal efficiency")
                g_fake.Write()

                g_rej = ROOT.TGraph(len(rejection), np.asarray(rejection, dtype="float64"), np.asarray(tpr, dtype="float64"))
                g_rej.SetName(f"gROC_{product}_cent_{key}_signalEfficiency_vs_backgroundRejection")
                g_rej.SetTitle(f"{label} {cent_label(lo, hi)};Background rejection;Signal efficiency")
                g_rej.Write()

            cent_hist = {}
            for class_name, class_value in (("background", 0), ("signal", 1)):
                values = score[finite & cent_mask & (y == class_value)]
                counts, _ = np.histogram(values, bins=score_bins)
                width = np.diff(score_bins)
                density = counts.astype(float)
                if density.sum() > 0:
                    density = density / density.sum() / width
                cent_hist[class_name] = {
                    "counts": counts.astype(int).tolist(),
                    "density": density.astype(float).tolist(),
                    "entries": int(len(values)),
                }

                h = ROOT.TH1D(
                    f"hScore_{product}_cent_{key}_{class_name}",
                    f"{label} {cent_label(lo, hi)} {class_name};BDT score;Candidates",
                    len(score_bins) - 1,
                    score_bins.astype("float64"),
                )
                h.SetDirectory(fout)
                for val in values:
                    h.Fill(float(val))
                h.Write()

                h_norm = h.Clone(f"hScore_{product}_cent_{key}_{class_name}_unitArea")
                if h_norm.Integral() > 0:
                    h_norm.Scale(1.0 / h_norm.Integral(), "width")
                h_norm.SetTitle(f"{label} {cent_label(lo, hi)} {class_name};BDT score;Unit-normalized candidates")
                h_norm.Write()
            product_hist["by_centrality"][key] = cent_hist

        if roc_by_cent:
            roc_payload["products"].setdefault(
                product,
                {"label": label, "color": PLOT_COLORS.get(product, "black")},
            )
            roc_payload["products"][product]["by_centrality"] = roc_by_cent
        roc_by_pt = {}
        et = frame["cluster_Et"].astype("float64")
        for lo, hi in PT_BINS:
            key = range_key(lo, hi)
            pt_mask = np.isfinite(et) & (et >= lo) & (et < hi)
            auc_pt, roc_pt = auc_for(y[pt_mask], score[pt_mask])
            if roc_pt is not None:
                fpr, tpr, thresholds = roc_pt
                rejection = 1.0 - fpr
                roc_by_pt[key] = {
                    "label": f"{lo:g}-{hi:g} GeV",
                    "auc": auc_pt,
                    "background_fake_rate": np.asarray(fpr, dtype=float).tolist(),
                    "background_rejection": np.asarray(rejection, dtype=float).tolist(),
                    "signal_efficiency": np.asarray(tpr, dtype=float).tolist(),
                    "threshold": [
                        None if (not np.isfinite(x)) else float(x)
                        for x in np.asarray(thresholds, dtype=float)
                    ],
                }

                g_fake = ROOT.TGraph(len(fpr), np.asarray(fpr, dtype="float64"), np.asarray(tpr, dtype="float64"))
                g_fake.SetName(f"gROC_{product}_pt_{key}_signalEfficiency_vs_backgroundFakeRate")
                g_fake.SetTitle(f"{label} {lo:g}-{hi:g} GeV;Background fake rate;Signal efficiency")
                g_fake.Write()

                g_rej = ROOT.TGraph(len(rejection), np.asarray(rejection, dtype="float64"), np.asarray(tpr, dtype="float64"))
                g_rej.SetName(f"gROC_{product}_pt_{key}_signalEfficiency_vs_backgroundRejection")
                g_rej.SetTitle(f"{label} {lo:g}-{hi:g} GeV;Background rejection;Signal efficiency")
                g_rej.Write()

            pt_hist = {}
            for class_name, class_value in (("background", 0), ("signal", 1)):
                values = score[finite & pt_mask & (y == class_value)]
                counts, _ = np.histogram(values, bins=score_bins)
                width = np.diff(score_bins)
                density = counts.astype(float)
                if density.sum() > 0:
                    density = density / density.sum() / width
                pt_hist[class_name] = {
                    "counts": counts.astype(int).tolist(),
                    "density": density.astype(float).tolist(),
                    "entries": int(len(values)),
                }
                h = ROOT.TH1D(
                    f"hScore_{product}_pt_{key}_{class_name}",
                    f"{label} {lo:g}-{hi:g} GeV {class_name};BDT score;Candidates",
                    len(score_bins) - 1,
                    score_bins.astype("float64"),
                )
                h.SetDirectory(fout)
                for val in values:
                    h.Fill(float(val))
                h.Write()
                h_norm = h.Clone(f"hScore_{product}_pt_{key}_{class_name}_unitArea")
                if h_norm.Integral() > 0:
                    h_norm.Scale(1.0 / h_norm.Integral(), "width")
                h_norm.SetTitle(f"{label} {lo:g}-{hi:g} GeV {class_name};BDT score;Unit-normalized candidates")
                h_norm.Write()
            product_hist.setdefault("by_pt", {})[key] = pt_hist
        if roc_by_pt:
            roc_payload["products"].setdefault(
                product,
                {"label": label, "color": PLOT_COLORS.get(product, "black")},
            )
            roc_payload["products"][product]["by_pt"] = roc_by_pt
        hist_payload["products"][product] = product_hist

    fout.Close()
    (outdir / "roc_points.json").write_text(json.dumps(json_ready(roc_payload), indent=2, sort_keys=True) + "\n")
    (outdir / "score_histograms.json").write_text(json.dumps(json_ready(hist_payload), indent=2, sort_keys=True) + "\n")


def make_plots(outdir: Path, frame, scores, metrics):
    import matplotlib.pyplot as plt
    import numpy as np

    y = frame["is_signal"].astype("int32")
    et = frame["cluster_Et"].astype("float64")
    cent = frame["centrality"].astype("float64")
    eiso = frame["reco_eiso"].astype("float64")

    colors = {
        "centINDcontrol": "tab:blue",
        "centAsFeat": "tab:orange",
        "centDepBDTs": "tab:green",
    }

    fig, ax = plt.subplots(figsize=(7.0, 6.0))
    for product, score in scores.items():
        auc = metrics["products"][product]["auc_inclusive"]
        roc = metrics["products"][product].get("_roc")
        if roc is None:
            continue
        fpr, tpr, _ = roc
        label = PLOT_LABELS.get(product, product)
        ax.plot(fpr, tpr, lw=2, color=colors[product], label=f"{label} AUC={auc:.3f}")
    ax.plot([0, 1], [0, 1], "k--", lw=1)
    ax.set_xlabel("Background fake rate")
    ax.set_ylabel("Signal efficiency")
    ax.set_title("AuAu tight-BDT ROC on embedded simulation")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(outdir / "roc_inclusive.png", dpi=170)
    plt.close(fig)

    for product, score in scores.items():
        fig, ax = plt.subplots(figsize=(7.0, 6.0))
        any_curve = False
        for lo, hi in CENT_BINS:
            mask = np.isfinite(cent) & (cent >= lo) & (cent < hi)
            auc, roc = auc_for(y[mask], score[mask])
            if roc is None:
                continue
            fpr, tpr, _ = roc
            any_curve = True
            ax.plot(fpr, tpr, lw=2, label=f"{cent_label(lo, hi)} AUC={auc:.3f}")
        if any_curve:
            ax.plot([0, 1], [0, 1], "k--", lw=1)
            ax.set_xlabel("Background fake rate")
            ax.set_ylabel("Signal efficiency")
            ax.set_title(f"{PLOT_LABELS.get(product, product)} ROC by centrality")
            ax.grid(alpha=0.25)
            ax.legend(fontsize=9)
            fig.tight_layout()
            fig.savefig(outdir / f"roc_by_centrality_{product}.png", dpi=170)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(7.0, 6.0))
        any_curve = False
        for lo, hi in PT_BINS:
            mask = np.isfinite(et) & (et >= lo) & (et < hi)
            auc, roc = auc_for(y[mask], score[mask])
            if roc is None:
                continue
            fpr, tpr, _ = roc
            any_curve = True
            ax.plot(fpr, tpr, lw=2, label=f"{lo:g}-{hi:g} GeV AUC={auc:.3f}")
        if any_curve:
            ax.plot([0, 1], [0, 1], "k--", lw=1)
            ax.set_xlabel("Background fake rate")
            ax.set_ylabel("Signal efficiency")
            ax.set_title(f"{PLOT_LABELS.get(product, product)} ROC by photon E$_T$")
            ax.grid(alpha=0.25)
            ax.legend(fontsize=8)
            fig.tight_layout()
            fig.savefig(outdir / f"roc_by_pt_{product}.png", dpi=170)
        plt.close(fig)

        finite = np.isfinite(score)
        fig, ax = plt.subplots(figsize=(7.0, 5.2))
        bins = np.linspace(0, 1, 51)
        ax.hist(score[finite & (y == 0)], bins=bins, histtype="stepfilled", alpha=0.45, density=True, label="background")
        ax.hist(score[finite & (y == 1)], bins=bins, histtype="step", lw=2.2, density=True, label="signal")
        ax.set_yscale("log")
        ax.set_xlabel("BDT score")
        ax.set_ylabel("Normalized candidates")
        ax.set_title(f"{PLOT_LABELS.get(product, product)}: signal/background score separation")
        ax.grid(alpha=0.25)
        ax.legend()
        fig.tight_layout()
        fig.savefig(outdir / f"score_overlay_{product}.png", dpi=170)
        plt.close(fig)

        fig, axes = plt.subplots(1, len(CENT_BINS), figsize=(14.4, 4.4), sharey=True)
        for ax, (lo, hi) in zip(axes, CENT_BINS):
            mask = finite & np.isfinite(cent) & (cent >= lo) & (cent < hi)
            if (mask & (y == 0)).sum() > 0:
                ax.hist(score[mask & (y == 0)], bins=bins, histtype="stepfilled", alpha=0.45, density=True, label="background")
            if (mask & (y == 1)).sum() > 0:
                ax.hist(score[mask & (y == 1)], bins=bins, histtype="step", lw=2.0, density=True, label="signal")
            auc = metrics["products"][product]["auc_by_centrality"].get(range_key(lo, hi), math.nan)
            auc_text = f"AUC={auc:.3f}" if math.isfinite(auc) else "AUC=n/a"
            ax.set_title(f"{cent_label(lo, hi)}  {auc_text}")
            ax.set_xlabel("BDT score")
            ax.set_yscale("log")
            ax.grid(alpha=0.25)
        axes[0].set_ylabel("Normalized candidates")
        axes[-1].legend(fontsize=9)
        fig.suptitle(f"{PLOT_LABELS.get(product, product)} score separation by centrality", y=1.02)
        fig.tight_layout()
        fig.savefig(outdir / f"score_overlay_by_centrality_{product}.png", dpi=170, bbox_inches="tight")
        plt.close(fig)

        for lo, hi in CENT_BINS:
            mask = finite & np.isfinite(cent) & (cent >= lo) & (cent < hi)
            if mask.sum() < 10:
                continue
            fig, ax = plt.subplots(figsize=(7.0, 5.2))
            ax.hist(score[mask & (y == 0)], bins=bins, histtype="stepfilled", alpha=0.45, density=True, label="background")
            ax.hist(score[mask & (y == 1)], bins=bins, histtype="step", lw=2.2, density=True, label="signal")
            ax.set_yscale("log")
            ax.set_xlabel("BDT score")
            ax.set_ylabel("Normalized candidates")
            ax.set_title(f"{PLOT_LABELS.get(product, product)}: {cent_label(lo, hi)} score separation")
            ax.grid(alpha=0.25)
            ax.legend()
            fig.tight_layout()
            fig.savefig(outdir / f"score_overlay_{product}_cent{range_key(lo, hi)}.png", dpi=170)
            plt.close(fig)

        for x, xname, xlabel in ((et, "clusterEt", r"cluster $E_T$ [GeV]"), (eiso, "recoEiso", r"reco $E_T^{iso}$ [GeV]")):
            mask = finite & np.isfinite(x)
            if mask.sum() < 10:
                continue
            fig, ax = plt.subplots(figsize=(7.0, 5.4))
            h = ax.hist2d(x[mask], score[mask], bins=(70, 60), range=((0, 60), (0, 1)), cmap="viridis")
            fig.colorbar(h[3], ax=ax, label="candidates")
            ax.set_xlabel(xlabel)
            ax.set_ylabel("BDT score")
            ax.set_title(f"{PLOT_LABELS.get(product, product)}: BDT score vs {xlabel}")
            fig.tight_layout()
            fig.savefig(outdir / f"score_vs_{xname}_{product}.png", dpi=170)
            plt.close(fig)

        heat = np.full((len(CENT_BINS), len(PT_BINS)), np.nan, dtype="float64")
        for ic, (clo, chi) in enumerate(CENT_BINS):
            for ip, (plo, phi) in enumerate(PT_BINS):
                mask = (
                    np.isfinite(cent)
                    & (cent >= clo)
                    & (cent < chi)
                    & np.isfinite(et)
                    & (et >= plo)
                    & (et < phi)
                )
                heat[ic, ip], _ = auc_for(y[mask], score[mask])
        if np.isfinite(heat).any():
            fig, ax = plt.subplots(figsize=(7.4, 4.8))
            im = ax.imshow(heat, vmin=0.5, vmax=1.0, cmap="viridis", aspect="auto")
            ax.set_xticks(np.arange(len(PT_BINS)), [f"{lo:g}-{hi:g}" for lo, hi in PT_BINS])
            ax.set_yticks(np.arange(len(CENT_BINS)), [cent_label(lo, hi) for lo, hi in CENT_BINS])
            ax.set_xlabel(r"cluster $E_T$ [GeV]")
            ax.set_ylabel("Centrality")
            ax.set_title(f"{PLOT_LABELS.get(product, product)} AUC by centrality and $E_T$")
            for ic in range(len(CENT_BINS)):
                for ip in range(len(PT_BINS)):
                    if np.isfinite(heat[ic, ip]):
                        ax.text(ip, ic, f"{heat[ic, ip]:.3f}", ha="center", va="center", color="white", fontsize=8)
            fig.colorbar(im, ax=ax, label="AUC")
            fig.tight_layout()
            fig.savefig(outdir / f"auc_heatmap_cent_pt_{product}.png", dpi=170)
            plt.close(fig)

    products = list(PRODUCTS)
    x = np.arange(len(products))
    aucs = [metrics["products"][p]["auc_inclusive"] for p in products]
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    ax.bar(x, aucs, color=[colors[p] for p in products])
    ax.set_xticks(x, [PLOT_LABELS.get(p, p) for p in products], rotation=20, ha="right")
    ax.set_ylim(0.45, 1.0)
    ax.set_ylabel("Inclusive AUC")
    ax.set_title("AuAu tight-BDT inclusive AUC")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "auc_summary.png", dpi=170)
    plt.close(fig)


def build_metrics(frame, scores, args):
    import numpy as np

    y = frame["is_signal"].astype("int32")
    et = frame["cluster_Et"].astype("float64")
    cent = frame["centrality"].astype("float64")
    eiso = frame["reco_eiso"].astype("float64")
    metrics = {"products": {}}
    for product, score in scores.items():
        eligible = np.isin(y, [0, 1])
        if product == "centDepBDTs":
            eligible &= np.isfinite(cent) & (cent >= 0.0) & (cent < 80.0)
        finite = eligible & np.isfinite(score)
        auc, roc = auc_for(y[eligible], score[eligible])
        by_pt = {}
        for lo, hi in PT_BINS:
            mask = eligible & (et >= lo) & (et < hi)
            a, _ = auc_for(y[mask], score[mask])
            by_pt[f"{lo:g}_{hi:g}"] = a
        by_cent = {}
        for lo, hi in CENT_BINS:
            mask = eligible & (cent >= lo) & (cent < hi)
            a, _ = auc_for(y[mask], score[mask])
            by_cent[f"{lo:g}_{hi:g}"] = a
        corr_mask = finite & np.isfinite(eiso)
        corr = math.nan
        if corr_mask.sum() >= 3:
            corr = float(np.corrcoef(score[corr_mask], eiso[corr_mask])[0, 1])
        metrics["products"][product] = {
            "auc_inclusive": auc,
            "_roc": roc,
            "auc_by_pt": by_pt,
            "auc_by_centrality": by_cent,
            "eligible_entries": int(eligible.sum()),
            "finite_scores": int(finite.sum()),
            "finite_score_fraction": float(finite.sum() / max(1, eligible.sum())),
            "score_eiso_pearson": corr,
            "signal_score_mean": float(np.nanmean(score[finite & (y == 1)])) if (finite & (y == 1)).any() else math.nan,
            "background_score_mean": float(np.nanmean(score[finite & (y == 0)])) if (finite & (y == 0)).any() else math.nan,
        }
    return metrics


def ensure_product_metrics(metrics):
    defaults = {
        "auc_inclusive": math.nan,
        "_roc": None,
        "auc_by_pt": {f"{lo:g}_{hi:g}": math.nan for lo, hi in PT_BINS},
        "auc_by_centrality": {f"{lo:g}_{hi:g}": math.nan for lo, hi in CENT_BINS},
        "eligible_entries": 0,
        "finite_scores": 0,
        "finite_score_fraction": 0.0,
        "score_eiso_pearson": math.nan,
        "signal_score_mean": math.nan,
        "background_score_mean": math.nan,
    }
    metrics.setdefault("products", {})
    for product in PRODUCTS:
        metrics["products"].setdefault(product, dict(defaults))
    return metrics


def json_ready(obj):
    import numpy as np

    if isinstance(obj, dict):
        return {k: json_ready(v) for k, v in obj.items() if k != "_roc"}
    if isinstance(obj, list):
        return [json_ready(v) for v in obj]
    if isinstance(obj, tuple):
        return [json_ready(v) for v in obj]
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)):
        return None
    return obj


def write_parseable_summary(path: Path, status: str, args, counts, metrics, notes):
    products = metrics["products"]
    finite_fracs = [products[p]["finite_score_fraction"] for p in products]
    with path.open("w") as f:
        f.write("RECOILJETS_AUAU_TIGHT_BDT_SIM_VALIDATION_V1\n")
        f.write(f"status={status}\n")
        f.write(f"source={args.source}\n")
        f.write(f"model_dir={args.model_dir}\n")
        f.write(f"total_entries={counts['total_entries']}\n")
        f.write(f"signal_entries={counts['signal_entries']}\n")
        f.write(f"background_entries={counts['background_entries']}\n")
        f.write(f"scored_entries={counts['scored_entries']}\n")
        for product in PRODUCTS:
            f.write(f"{product}_auc={products[product]['auc_inclusive']:.6g}\n")
            f.write(f"{product}_finite_score_fraction={products[product]['finite_score_fraction']:.6g}\n")
            f.write(f"{product}_score_eiso_pearson={products[product]['score_eiso_pearson']:.6g}\n")
            f.write(f"{product}_signal_score_mean={products[product]['signal_score_mean']:.6g}\n")
            f.write(f"{product}_background_score_mean={products[product]['background_score_mean']:.6g}\n")
        f.write(f"finite_score_fraction={min(finite_fracs) if finite_fracs else 0.0:.6g}\n")
        f.write(f"report_dir={path.parent}\n")
        f.write(f"metrics_json={path.parent / 'validation_metrics.json'}\n")
        f.write(f"deep_diagnostics_json={path.parent / 'validation_deep_diagnostics.json'}\n")
        f.write(f"curves_root={path.parent / 'validation_curves.root'}\n")
        f.write(f"auc_table_csv={path.parent / 'validation_auc_table.csv'}\n")
        f.write(f"threshold_table_csv={path.parent / 'validation_threshold_table.csv'}\n")
        f.write(f"feature_summary_csv={path.parent / 'validation_feature_summary.csv'}\n")
        f.write(f"notes={'; '.join(notes) if notes else 'none'}\n")
        f.write("next_action=If status is READY, inspect the PNG/JSON outputs, then run the constrained data+MC BDT-variant validation.\n")


def main() -> int:
    args = parse_args()
    require_imports()
    import numpy as np

    outdir = args.outdir or (args.source / "reports" / f"model_validation_{os.environ.get('RJ_AUAU_TIGHT_BDT_VALIDATE_STAMP') or __import__('datetime').datetime.now().strftime('%Y%m%d_%H%M%S')}")
    outdir.mkdir(parents=True, exist_ok=True)

    if args.merge_score_caches is not None:
        frame, scores, counts = load_score_caches(args.merge_score_caches)
    else:
        paths = read_manifest(args.manifest) if args.manifest is not None else discover_manifest(args.source)
        models = load_models(args.model_dir)
        frame, counts = collect_rows(paths, args.tree, args.score_max_rows, args.progress_every)
        scores = score_rows(frame, models) if counts["scored_entries"] > 0 else {p: np.array([], dtype="float32") for p in PRODUCTS}
        if args.write_score_cache is not None:
            write_score_cache(args.write_score_cache, frame, scores, counts)

    notes = []
    if counts["missing_tree_files"]:
        notes.append(f"missing_tree_files={len(counts['missing_tree_files'])}")
    if counts["missing_branches"]:
        notes.append(f"missing_branch_files={len(counts['missing_branches'])}")
    if counts["total_entries"] <= 0 or counts["signal_entries"] <= 0 or counts["background_entries"] <= 0:
        notes.append("nonzero signal/background requirement failed")
    if counts["scored_entries"] <= 0:
        notes.append("no rows available for scoring")

    metrics = build_metrics(frame, scores, args) if counts["scored_entries"] > 0 else {"products": {}}
    metrics = ensure_product_metrics(metrics)

    status = "READY"
    if counts["total_entries"] <= 0 or counts["signal_entries"] <= 0 or counts["background_entries"] <= 0:
        status = "CHECK"
    for product, product_metrics in metrics["products"].items():
        auc = product_metrics["auc_inclusive"]
        finite_frac = product_metrics["finite_score_fraction"]
        sig_mean = product_metrics["signal_score_mean"]
        bkg_mean = product_metrics["background_score_mean"]
        if not math.isfinite(auc) or auc < args.min_auc:
            status = "CHECK"
            notes.append(f"{product}_auc_below_threshold")
        if finite_frac < args.min_finite_fraction:
            status = "CHECK"
            notes.append(f"{product}_finite_fraction_below_threshold")
        if math.isfinite(sig_mean) and math.isfinite(bkg_mean) and sig_mean <= bkg_mean:
            status = "CHECK"
            notes.append(f"{product}_signal_mean_not_above_background_mean")
    if set(metrics["products"]) != set(PRODUCTS):
        status = "CHECK"
        notes.append("not_all_products_scored")

    if not args.no_plots:
        diagnostics = build_deep_diagnostics(frame, scores)
        (outdir / "validation_deep_diagnostics.json").write_text(
            json.dumps(json_ready(diagnostics), indent=2, sort_keys=True) + "\n"
        )
        write_diagnostic_tables(outdir, diagnostics)
        make_plots(outdir, frame, scores, metrics)
        write_curve_exports(outdir, frame, scores, metrics)
    metrics_json = json_ready({"counts": counts, **metrics, "status": status, "notes": notes})
    (outdir / "validation_metrics.json").write_text(json.dumps(metrics_json, indent=2, sort_keys=True) + "\n")
    write_parseable_summary(outdir / "validation_summary.txt", status, args, counts, metrics, notes)
    print((outdir / "validation_summary.txt").read_text(), end="")
    return 0 if status == "READY" else 3


if __name__ == "__main__":
    raise SystemExit(main())
