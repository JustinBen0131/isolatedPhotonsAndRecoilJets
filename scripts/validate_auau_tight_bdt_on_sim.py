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
WIDTH_3X3_FEATURES = ["cluster_weta33_cogx", "cluster_wphi33_cogx"]
WIDTH_RATIO_FEATURES = ["cluster_weta_over_wphi", "cluster_weta33_over_wphi33"]
SHAPE_RESIDUAL_FEATURES = [
    "shape_tail_e37_e53_logratio",
    "shape_long_tail_e17_e71_logratio",
    "shape_width_core_full_logratio",
    "shape_core_tail_tension",
    "shape_compactness_gradient",
]
SHAPE_TEMPLATE_FEATURES = [
    "shape_template_diag_chi2",
    "shape_template_max_abs_z",
]
SHAPE_TEMPLATE_BASIS = [
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "cluster_weta33_cogx",
    "cluster_wphi33_cogx",
    "e11_over_e33",
    "e32_over_e35",
    "shape_tail_e37_e53_logratio",
    "shape_width_core_full_logratio",
    "shape_core_tail_tension",
]
SHAPE_TEMPLATE_PT_EDGES = [15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 30.0, 35.0]
SHAPE_TEMPLATE_CENT_EDGES = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0]
ISOLATION_DIAGNOSTIC_FEATURES = [
    "reco_eiso_clip30",
    "reco_eiso_over_cluster_Et",
    "reco_eiso_signed_log1p",
]
ISOLATION_CONE_RAW_FEATURES = [
    "reco_eiso_r30",
    "reco_eiso_r40",
]
EXTENDED_SHOWER_FEATURES = [
    "cluster_weta35_cogx",
    "cluster_wphi53_cogx",
    "cluster_w32",
    "cluster_w52",
    "cluster_w72",
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
]
DERIVED_FEATURE_DEPS = {
    "cluster_weta_over_wphi": ["cluster_weta_cogx", "cluster_wphi_cogx"],
    "cluster_weta33_over_wphi33": ["cluster_weta33_cogx", "cluster_wphi33_cogx"],
    "shape_tail_e37_e53_logratio": ["e22_over_e37", "e22_over_e53"],
    "shape_long_tail_e17_e71_logratio": ["e11_over_e17", "e11_over_e71"],
    "shape_width_core_full_logratio": [
        "cluster_weta_cogx",
        "cluster_wphi_cogx",
        "cluster_weta33_cogx",
        "cluster_wphi33_cogx",
    ],
    "shape_core_tail_tension": ["e11_over_e33", "e22_over_e37", "e22_over_e53"],
    "shape_compactness_gradient": ["e11_over_e33", "e32_over_e35"],
    "shape_template_diag_chi2": [
        "is_signal",
        "cluster_Et",
        "centrality",
        "cluster_weta_cogx",
        "cluster_wphi_cogx",
        "cluster_weta33_cogx",
        "cluster_wphi33_cogx",
        "e11_over_e33",
        "e32_over_e35",
        "e22_over_e37",
        "e22_over_e53",
    ],
    "shape_template_max_abs_z": [
        "is_signal",
        "cluster_Et",
        "centrality",
        "cluster_weta_cogx",
        "cluster_wphi_cogx",
        "cluster_weta33_cogx",
        "cluster_wphi33_cogx",
        "e11_over_e33",
        "e32_over_e35",
        "e22_over_e37",
        "e22_over_e53",
    ],
    "reco_eiso_clip30": ["reco_eiso"],
    "reco_eiso_over_cluster_Et": ["reco_eiso", "cluster_Et"],
    "reco_eiso_signed_log1p": ["reco_eiso"],
}

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
    "centAsFeat3x3_pt5to40": "Centrality input, 3x3 shower widths",
    "centAsFeat_pt15to30": "Base widths, 15-30 GeV",
    "centAsFeat3x3_pt15to30": "3x3 widths, 15-30 GeV",
    "centAsFeatBase3x3_pt15to30": "Base + 3x3 widths, 15-30 GeV",
    "noCent_pt1535": "No centrality input, base + 3x3 widths",
    "centInput_pt1535": "Centrality input, base + 3x3 widths",
    "ptFine_noCent": "Fine E_T bins, no centrality input",
    "ptFine_centInput": "Fine E_T bins, centrality input",
    "ptFine_cent3": "Fine E_T x 3 centrality bins",
    "ptFine_cent7": "Fine E_T x 7 centrality bins",
    "centDepBDTs": "Centrality-specific models",
    "isoBDT_global15to35_EtCent_full": "Isolation-input BDT, global 15-35 GeV",
    "isoBDT_ptFine15to35_cent7_full": "Isolation-input BDT, fine E_T x 7 centrality bins",
    "baseBDT_v3E_withCentrality_w33_E22E37": "Base v3E + centrality + 3x3 widths + E22/E37",
    "baseBDT_v3E_withCentrality_w33_E22E53": "Base v3E + centrality + 3x3 widths + E22/E53",
    "baseBDT_v3E_withCentrality_w33_E22E37_E22E53": "Base v3E + centrality + 3x3 widths + E22/E37 + E22/E53",
    "globalEtCent1535_bdt_eisoR30_ptCent3": "Fine E_T x 3 centrality bins + raw R=0.3 isolation ET",
    "globalEtCent1535_bdt_eisoR30_ptCent7": "Fine E_T x 7 centrality bins + raw R=0.3 isolation ET",
    "globalEtCent1535_bdt_eisoR30R40_ptCent3": "Fine E_T x 3 centrality bins + raw R=0.3 and R=0.4 isolation ET",
    "globalEtCent1535_bdt_eisoR30R40_ptCent7": "Fine E_T x 7 centrality bins + raw R=0.3 and R=0.4 isolation ET",
}

PLOT_COLORS = {
    "centINDcontrol": "#1f77b4",
    "centAsFeat": "#ff7f0e",
    "centDepBDTs": "#2ca02c",
    "noCent_pt1535": "#4D4D4D",
    "centInput_pt1535": "#0072B2",
    "ptFine_noCent": "#999999",
    "ptFine_centInput": "#D55E00",
    "ptFine_cent3": "#009E73",
    "ptFine_cent7": "#CC79A7",
    "isoBDT_global15to35_EtCent_full": "#7E57C2",
    "isoBDT_ptFine15to35_cent7_full": "#C2185B",
    "globalEtCent1535_bdt_eisoR30_ptCent3": "#009E73",
    "globalEtCent1535_bdt_eisoR30_ptCent7": "#0072B2",
    "globalEtCent1535_bdt_eisoR30R40_ptCent3": "#D55E00",
    "globalEtCent1535_bdt_eisoR30R40_ptCent7": "#CC79A7",
}
FALLBACK_COLORS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
]

PT_BINS = [(6, 10), (10, 15), (15, 20), (20, 25), (25, 35)]
CENT_BINS = [(0, 20), (20, 50), (50, 80)]
DIAGNOSTIC_FEATURES = []
for _name in BASE_FEATURES + WIDTH_3X3_FEATURES + EXTENDED_SHOWER_FEATURES + WIDTH_RATIO_FEATURES + SHAPE_RESIDUAL_FEATURES + SHAPE_TEMPLATE_FEATURES + ISOLATION_DIAGNOSTIC_FEATURES + ISOLATION_CONE_RAW_FEATURES + ["centrality", "reco_eiso"]:
    if _name not in DIAGNOSTIC_FEATURES:
        DIAGNOSTIC_FEATURES.append(_name)
MANDATORY_CACHE_COLUMNS = ["is_signal", "cluster_Et", "cluster_Eta", "centrality", "reco_eiso"]
EFFICIENCY_TARGETS = [0.50, 0.70, 0.80, 0.90, 0.95]
FAKE_RATE_TARGETS = [0.01, 0.02, 0.05, 0.10, 0.20]


def expand_required_columns(features: list[str]) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for feature in features:
        for dep in DERIVED_FEATURE_DEPS.get(feature, [feature]):
            if dep not in seen:
                out.append(dep)
                seen.add(dep)
    return out


def add_derived_features(frame):
    import numpy as np

    eps = 1.0e-6

    def safe_ratio(num: str, den: str):
        a = np.asarray(frame[num], dtype="float64")
        b = np.asarray(frame[den], dtype="float64")
        out = np.full(len(a), np.nan, dtype="float32")
        mask = np.isfinite(a) & np.isfinite(b) & (np.abs(b) > 1.0e-12)
        out[mask] = (a[mask] / b[mask]).astype("float32")
        return out

    def safe_log_ratio(num: str, den: str):
        a = np.asarray(frame[num], dtype="float64")
        b = np.asarray(frame[den], dtype="float64")
        out = np.full(len(a), np.nan, dtype="float32")
        mask = np.isfinite(a) & np.isfinite(b) & (a > eps) & (b > eps)
        out[mask] = np.log(a[mask] / b[mask]).astype("float32")
        return out

    def safe_logit(name: str):
        x = np.asarray(frame[name], dtype="float64")
        out = np.full(len(x), np.nan, dtype="float32")
        mask = np.isfinite(x)
        clipped = np.clip(x[mask], eps, 1.0 - eps)
        out[mask] = np.log(clipped / (1.0 - clipped)).astype("float32")
        return out

    if "cluster_weta_over_wphi" not in frame and {"cluster_weta_cogx", "cluster_wphi_cogx"}.issubset(frame):
        frame["cluster_weta_over_wphi"] = safe_ratio("cluster_weta_cogx", "cluster_wphi_cogx")
    if "cluster_weta33_over_wphi33" not in frame and {"cluster_weta33_cogx", "cluster_wphi33_cogx"}.issubset(frame):
        frame["cluster_weta33_over_wphi33"] = safe_ratio("cluster_weta33_cogx", "cluster_wphi33_cogx")
    if "shape_tail_e37_e53_logratio" not in frame and {"e22_over_e37", "e22_over_e53"}.issubset(frame):
        frame["shape_tail_e37_e53_logratio"] = safe_log_ratio("e22_over_e37", "e22_over_e53")
    if "shape_long_tail_e17_e71_logratio" not in frame and {"e11_over_e17", "e11_over_e71"}.issubset(frame):
        frame["shape_long_tail_e17_e71_logratio"] = safe_log_ratio("e11_over_e17", "e11_over_e71")
    if (
        "shape_width_core_full_logratio" not in frame
        and {"cluster_weta_cogx", "cluster_wphi_cogx", "cluster_weta33_cogx", "cluster_wphi33_cogx"}.issubset(frame)
    ):
        core = safe_ratio("cluster_weta33_cogx", "cluster_wphi33_cogx").astype("float64")
        full = safe_ratio("cluster_weta_cogx", "cluster_wphi_cogx").astype("float64")
        out = np.full(len(core), np.nan, dtype="float32")
        mask = np.isfinite(core) & np.isfinite(full) & (core > eps) & (full > eps)
        out[mask] = np.log(core[mask] / full[mask]).astype("float32")
        frame["shape_width_core_full_logratio"] = out
    if "shape_core_tail_tension" not in frame and {"e11_over_e33", "e22_over_e37", "e22_over_e53"}.issubset(frame):
        frame["shape_core_tail_tension"] = (
            safe_logit("e11_over_e33")
            - 0.5 * (safe_logit("e22_over_e37") + safe_logit("e22_over_e53"))
        ).astype("float32")
    if "shape_compactness_gradient" not in frame and {"e11_over_e33", "e32_over_e35"}.issubset(frame):
        frame["shape_compactness_gradient"] = (safe_logit("e11_over_e33") - safe_logit("e32_over_e35")).astype("float32")

    if (
        {"shape_template_diag_chi2", "shape_template_max_abs_z"}.difference(frame)
        and {"is_signal", "cluster_Et", "centrality", *SHAPE_TEMPLATE_BASIS}.issubset(frame)
    ):
        basis = SHAPE_TEMPLATE_BASIS
        values = np.column_stack([np.asarray(frame[name], dtype="float64") for name in basis])
        is_signal = np.asarray(frame["is_signal"], dtype="int32") == 1
        et = np.asarray(frame["cluster_Et"], dtype="float64")
        cent = np.asarray(frame["centrality"], dtype="float64")
        finite = np.isfinite(values).all(axis=1) & np.isfinite(et) & np.isfinite(cent)
        global_mask = is_signal & finite
        if global_mask.sum() >= 50:
            global_mu = np.nanmedian(values[global_mask], axis=0)
            global_mad = np.nanmedian(np.abs(values[global_mask] - global_mu), axis=0)
            global_sigma = np.maximum(1.4826 * global_mad, 1.0e-4)
        else:
            global_mu = np.nanmedian(values[finite], axis=0) if finite.any() else np.zeros(len(basis))
            global_sigma = np.nanstd(values[finite], axis=0) if finite.any() else np.ones(len(basis))
            global_sigma = np.maximum(global_sigma, 1.0e-4)

        chi2 = np.full(len(et), np.nan, dtype="float64")
        max_abs = np.full(len(et), np.nan, dtype="float64")
        for pt_lo, pt_hi in zip(SHAPE_TEMPLATE_PT_EDGES[:-1], SHAPE_TEMPLATE_PT_EDGES[1:]):
            pt_mask = (et >= pt_lo) & (et < pt_hi)
            for cent_lo, cent_hi in zip(SHAPE_TEMPLATE_CENT_EDGES[:-1], SHAPE_TEMPLATE_CENT_EDGES[1:]):
                row_mask = pt_mask & (cent >= cent_lo) & (cent < cent_hi) & finite
                if not row_mask.any():
                    continue
                sig_mask = row_mask & is_signal
                if sig_mask.sum() >= 50:
                    mu = np.nanmedian(values[sig_mask], axis=0)
                    mad = np.nanmedian(np.abs(values[sig_mask] - mu), axis=0)
                    sigma = np.maximum(1.4826 * mad, 1.0e-4)
                else:
                    mu = global_mu
                    sigma = global_sigma
                z = np.clip((values[row_mask] - mu) / sigma, -8.0, 8.0)
                chi2[row_mask] = np.mean(z * z, axis=1)
                max_abs[row_mask] = np.max(np.abs(z), axis=1)
        fallback = finite & ~np.isfinite(chi2)
        if fallback.any():
            z = np.clip((values[fallback] - global_mu) / global_sigma, -8.0, 8.0)
            chi2[fallback] = np.mean(z * z, axis=1)
            max_abs[fallback] = np.max(np.abs(z), axis=1)
        frame["shape_template_diag_chi2"] = chi2.astype("float32")
        frame["shape_template_max_abs_z"] = max_abs.astype("float32")

    cols = set(frame)
    if "reco_eiso_clip30" not in cols and "reco_eiso" in cols:
        reco_eiso = np.asarray(frame["reco_eiso"], dtype="float64")
        frame["reco_eiso_clip30"] = np.where(
            np.isfinite(reco_eiso) & (np.abs(reco_eiso) < 1.0e8),
            np.clip(reco_eiso, -20.0, 30.0),
            np.nan,
        ).astype("float32")
        cols.add("reco_eiso_clip30")
    if "reco_eiso_over_cluster_Et" not in cols and {"reco_eiso", "cluster_Et"}.issubset(cols):
        reco_eiso = np.asarray(frame["reco_eiso"], dtype="float64")
        cluster_et = np.asarray(frame["cluster_Et"], dtype="float64")
        out = np.full(len(reco_eiso), np.nan, dtype="float64")
        mask = (
            np.isfinite(reco_eiso)
            & np.isfinite(cluster_et)
            & (np.abs(reco_eiso) < 1.0e8)
            & (cluster_et > 1.0e-6)
        )
        out[mask] = np.clip(reco_eiso[mask] / cluster_et[mask], -2.0, 3.0)
        frame["reco_eiso_over_cluster_Et"] = out.astype("float32")
        cols.add("reco_eiso_over_cluster_Et")
    if "reco_eiso_signed_log1p" not in cols and "reco_eiso" in cols:
        reco_eiso = np.asarray(frame["reco_eiso"], dtype="float64")
        clipped = np.where(
            np.isfinite(reco_eiso) & (np.abs(reco_eiso) < 1.0e8),
            np.clip(reco_eiso, -20.0, 60.0),
            np.nan,
        )
        frame["reco_eiso_signed_log1p"] = (np.sign(clipped) * np.log1p(np.abs(clipped))).astype("float32")


def range_key(lo, hi) -> str:
    return f"{lo:g}_{hi:g}".replace(".", "p")


def cent_label(lo, hi) -> str:
    return f"{lo:g}-{hi:g}%"


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", type=Path, default=None)
    ap.add_argument("--model-dir", type=Path, default=None)
    ap.add_argument(
        "--model-registry",
        type=Path,
        default=None,
        help="Optional explicit model_registry*.json. If omitted, uses model-dir/model_registry.json when present.",
    )
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
    ap.add_argument("--target-signal-efficiency", type=float, default=float(os.environ.get("RJ_AUAU_BDT_TARGET_SIGNAL_EFF", "0.80")))
    ap.add_argument(
        "--working-point-mode",
        choices=["pt", "centpt"],
        default=os.environ.get("RJ_AUAU_BDT_WP_MODE", "centpt"),
        help="Threshold shape to derive: pt is old E_T-only; centpt is E_T x centrality local target-efficiency bins.",
    )
    ap.add_argument(
        "--working-point-pt-bins",
        default=os.environ.get("RJ_AUAU_BDT_WP_PT_BINS", ""),
        help="Comma-separated E_T edges for target-efficiency working points, e.g. 15,17,19,21,23,25,27,30,35.",
    )
    ap.add_argument(
        "--working-point-cent-bins",
        default=os.environ.get("RJ_AUAU_BDT_WP_CENT_BINS", "0,20,50,80"),
        help="Comma-separated centrality edges for target-efficiency working points, e.g. 0,20,50,80.",
    )
    ap.add_argument(
        "--derive-working-points-from-report",
        type=Path,
        default=None,
        help="Regenerate target-efficiency BDT working-point files from an existing validation report.",
    )
    return ap.parse_args()


def parse_edge_list(text: str, fallback):
    if not text:
        return list(fallback)
    edges = []
    for token in text.replace(":", ",").split(","):
        token = token.strip()
        if not token:
            continue
        edges.append(float(token))
    if len(edges) < 2:
        raise SystemExit(f"Need at least two bin edges, got: {text}")
    for a, b in zip(edges[:-1], edges[1:]):
        if not b > a:
            raise SystemExit(f"Bin edges must be strictly increasing, got: {text}")
    return edges


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


def product_specs_from_registry(model_dir: Path, registry_path: Path | None = None):
    registry = registry_path if registry_path is not None else model_dir / "model_registry.json"
    if not registry.is_file():
        if registry_path is not None:
            raise SystemExit(f"Requested model registry does not exist: {registry}")
        return None
    data = json.loads(registry.read_text())
    grouped = {}
    for spec in data.get("models", []):
        report = spec.get("report")
        if not report:
            continue
        if report.get("status") == "skipped":
            continue
        product = spec["product"]
        grouped.setdefault(product, {"features": spec["features"], "models": []})
        grouped[product]["models"].append(
            {
                "selector": {
                    "cent_range": spec.get("cent_range"),
                    "pt_range": spec.get("pt_range"),
                },
                "filename": Path(spec["output_tmva"]).name,
                "model_id": spec["model_id"],
            }
        )
    if not grouped:
        raise SystemExit(f"Registry exists but has no trained model reports: {registry}")
    return grouped


def active_product_specs(model_dir: Path, registry_path: Path | None = None):
    return product_specs_from_registry(model_dir, registry_path) or {
        product: {
            "features": spec["features"],
            "models": [
                {
                    "selector": {"cent_range": selector if selector != "all" else None, "pt_range": None},
                    "filename": filename,
                    "model_id": product if selector == "all" else f"{product}_{range_key(*selector)}",
                }
                for selector, filename in spec["models"]
            ],
        }
        for product, spec in PRODUCTS.items()
    }


def product_label(product: str) -> str:
    return PLOT_LABELS.get(product, product.replace("_", " "))


def product_color(product: str) -> str:
    if product in PLOT_COLORS:
        return PLOT_COLORS[product]
    return FALLBACK_COLORS[stable_color_index(product)]


def stable_color_index(text: str) -> int:
    return sum(ord(ch) for ch in text) % len(FALLBACK_COLORS)


def load_models(model_dir: Path, product_specs):
    import ROOT

    loaded = {}
    missing = []
    for product, spec in product_specs.items():
        entries = []
        for model_spec in spec["models"]:
            filename = model_spec["filename"]
            path = model_dir / filename
            if not path.is_file() or path.stat().st_size <= 0:
                missing.append(str(path))
                continue
            try:
                entries.append((model_spec["selector"], ROOT.TMVA.Experimental.RBDT("myBDT", str(path))))
            except Exception as exc:
                raise SystemExit(f"Failed to load TMVA model {path}: {exc}") from exc
        loaded[product] = entries
    if missing:
        raise SystemExit("Missing required TMVA model files:\n  " + "\n  ".join(missing))
    return loaded


def selector_mask(frame, selector):
    import numpy as np

    n = len(frame["is_signal"])
    mask = np.ones(n, dtype=bool)
    cent_range = selector.get("cent_range")
    if cent_range is not None:
        lo, hi = cent_range
        cent = frame["centrality"]
        mask &= np.isfinite(cent) & (cent >= lo) & (cent < hi)
    pt_range = selector.get("pt_range")
    if pt_range is not None:
        lo, hi = pt_range
        et = frame["cluster_Et"]
        mask &= np.isfinite(et) & (et >= lo) & (et < hi)
    return mask


def product_eligible_mask(frame, product_spec):
    import numpy as np

    models = product_spec.get("models", []) if product_spec else []
    if not models:
        return np.ones(len(frame["is_signal"]), dtype=bool)
    eligible = np.zeros(len(frame["is_signal"]), dtype=bool)
    for model in models:
        eligible |= selector_mask(frame, model.get("selector", {}))
    return eligible


def score_rows(frame, models, product_specs):
    import numpy as np
    import ROOT

    scores = {}
    for product, spec in product_specs.items():
        features = spec["features"]
        out = np.full(len(frame["is_signal"]), np.nan, dtype="float32")
        for selector, model in models[product]:
            mask = selector_mask(frame, selector)
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


def collect_rows(
    paths: list[Path],
    tree_name: str,
    score_max_rows: int,
    progress_every: int,
    feature_names: list[str] | None = None,
):
    import numpy as np
    import uproot

    # Only scoring inputs and mandatory QA columns may veto a ROOT file. The
    # broader DIAGNOSTIC_FEATURES list includes optional sidecar-only branches
    # such as cone-specific isolation, which are not present in every training
    # extraction and should not force otherwise valid files to be skipped.
    required = sorted(
        set(
            ["is_signal", "cluster_Et", "cluster_Eta", "centrality", "reco_eiso"]
            + BASE_FEATURES
            + expand_required_columns(list(feature_names or []))
        )
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
    add_derived_features(frame)
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
    score_parts = {}
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
            for name in ["is_signal"] + DIAGNOSTIC_FEATURES:
                if name in data:
                    frame_parts[name].append(data[name])
                elif name in MANDATORY_CACHE_COLUMNS:
                    raise SystemExit(
                        f"Missing diagnostic column {name} in score cache: {path}. "
                        "Regenerate validation caches with the updated validator."
                    )
            for key in data.files:
                if not key.startswith("score_"):
                    continue
                product = key[len("score_"):]
                score_parts.setdefault(product, []).append(data[key])
            cached_counts = json.loads(str(data["counts_json"].item())) if "counts_json" in data else {}
            counts["files"] += int(cached_counts.get("files", 0))
            counts["total_entries"] += int(cached_counts.get("total_entries", 0))
            counts["signal_entries"] += int(cached_counts.get("signal_entries", 0))
            counts["background_entries"] += int(cached_counts.get("background_entries", 0))
            counts["scored_entries"] += int(cached_counts.get("scored_entries", 0))
            counts["missing_tree_files"].extend(cached_counts.get("missing_tree_files", []))
            counts["missing_branches"].update(cached_counts.get("missing_branches", {}))

    frame = {}
    for name, parts in frame_parts.items():
        if name in MANDATORY_CACHE_COLUMNS or len(parts) == len(cache_paths):
            frame[name] = np.concatenate(parts) if parts else np.array([], dtype="float32")
    add_derived_features(frame)
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


def threshold_for_signal_efficiency(y, score, target):
    import numpy as np

    y = np.asarray(y)
    score = np.asarray(score)
    mask = np.isin(y, [0, 1]) & np.isfinite(score)
    sig = score[mask & (y == 1)]
    bkg = score[mask & (y == 0)]
    if sig.size <= 0 or bkg.size <= 0:
        return None
    # Tight is score > threshold. The lower (1-target) signal quantile gives
    # the requested retained signal fraction without assuming score calibration.
    threshold = float(np.quantile(sig, max(0.0, min(1.0, 1.0 - target))))
    sig_eff = float(np.mean(sig > threshold))
    bkg_fake = float(np.mean(bkg > threshold))
    return {
        "threshold": threshold,
        "signal_efficiency": sig_eff,
        "background_fake_rate": bkg_fake,
        "signal_entries": int(sig.size),
        "background_entries": int(bkg.size),
    }


def product_pt_edges(product_spec, fallback_edges=None):
    edges = set()
    for model in product_spec.get("models", []) if product_spec else []:
        pt_range = model.get("selector", {}).get("pt_range")
        if pt_range is not None and len(pt_range) == 2:
            edges.add(float(pt_range[0]))
            edges.add(float(pt_range[1]))
    if len(edges) >= 2:
        return sorted(edges)
    if fallback_edges:
        return list(fallback_edges)
    return [lo for lo, _ in PT_BINS] + [PT_BINS[-1][1]]


def derive_working_points(frame, scores, product_specs, target=0.80):
    import numpy as np

    y = frame["is_signal"].astype("int32")
    et = frame["cluster_Et"].astype("float64")
    products = {}
    for product, score in scores.items():
        spec = product_specs.get(product, {}) if product_specs else {}
        eligible = np.isin(y, [0, 1]) & np.isfinite(score)
        if spec:
            eligible &= product_eligible_mask(frame, spec)
        inclusive = threshold_for_signal_efficiency(y[eligible], score[eligible], target)
        if inclusive is None:
            products[product] = {
                "status": "skipped",
                "reason": "insufficient signal/background scores for target working point",
            }
            continue

        edges = product_pt_edges(spec)
        bins = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            mask = eligible & np.isfinite(et) & (et >= lo) & (et < hi)
            item = threshold_for_signal_efficiency(y[mask], score[mask], target)
            if item is None:
                continue
            item.update({"pt_min": float(lo), "pt_max": float(hi), "pt_center": float(0.5 * (lo + hi))})
            bins.append(item)

        mode = "linear"
        intercept = float(inclusive["threshold"])
        slope = 0.0
        fit_quality = {
            "max_abs_efficiency_error": math.nan,
            "max_abs_threshold_residual": math.nan,
            "n_fit_bins": len(bins),
        }
        if len(bins) >= 2:
            x = np.array([b["pt_center"] for b in bins], dtype="float64")
            t = np.array([b["threshold"] for b in bins], dtype="float64")
            w = np.sqrt(np.array([max(1, b["signal_entries"]) for b in bins], dtype="float64"))
            coeff = np.polyfit(x, t, deg=1, w=w)
            slope = float(coeff[0])
            intercept = float(coeff[1])
            residual = t - (intercept + slope * x)
            eff_errors = []
            for b in bins:
                mask = eligible & np.isfinite(et) & (et >= b["pt_min"]) & (et < b["pt_max"])
                sig_scores = score[mask & (y == 1)]
                if sig_scores.size:
                    eff_errors.append(float(np.mean(sig_scores > (intercept + slope * b["pt_center"])) - target))
            fit_quality = {
                "max_abs_efficiency_error": float(max((abs(v) for v in eff_errors), default=math.nan)),
                "max_abs_threshold_residual": float(np.max(np.abs(residual))) if residual.size else math.nan,
                "n_fit_bins": len(bins),
            }
            if (
                math.isfinite(fit_quality["max_abs_efficiency_error"])
                and fit_quality["max_abs_efficiency_error"] > 0.08
            ):
                mode = "binned"

        pt_min = float(edges[0]) if edges else math.nan
        pt_max = float(edges[-1]) if edges else math.nan
        if not math.isfinite(pt_min) or not math.isfinite(pt_max):
            finite_et = et[eligible & np.isfinite(et)]
            pt_min = float(np.nanmin(finite_et)) if finite_et.size else math.nan
            pt_max = float(np.nanmax(finite_et)) if finite_et.size else math.nan

        entry = {
            "status": "ready",
            "target_signal_efficiency": float(target),
            "mode": mode,
            "intercept": float(intercept),
            "slope": float(slope),
            "pt_min": pt_min,
            "pt_max": pt_max,
            "max_score": 1.0,
            "inclusive": inclusive,
            "bins": bins,
            "fit_quality": fit_quality,
        }
        if mode == "binned":
            entry["binned_edges"] = [float(bins[0]["pt_min"])] + [float(b["pt_max"]) for b in bins] if bins else [pt_min, pt_max]
            entry["binned_thresholds"] = [float(b["threshold"]) for b in bins] if bins else [float(inclusive["threshold"])]
        products[product] = entry
    return {
        "schema": "RECOILJETS_AUAU_BDT_WORKING_POINTS_V1",
        "target_signal_efficiency": float(target),
        "products": products,
    }


def runtime_entry(name, wp):
    if wp.get("status") != "ready":
        return None
    if wp.get("mode") == "grid2d":
        pt_edges = ";".join(f"{x:.8g}" for x in wp.get("pt_edges", []))
        cent_edges = ";".join(f"{x:.8g}" for x in wp.get("cent_edges", []))
        thresholds = []
        for row in wp.get("grid_thresholds", []):
            thresholds.extend(row)
        vals = ";".join(f"{float(x):.8g}" for x in thresholds)
        return (
            f"{name}|grid2d|{pt_edges}|{cent_edges}|{vals}|"
            f"{wp['pt_min']:.8g}|{wp['pt_max']:.8g}|{wp.get('max_score', 1.0):.8g}"
        )
    if wp.get("mode") == "binned":
        edges = ";".join(f"{x:.8g}" for x in wp.get("binned_edges", []))
        vals = ";".join(f"{x:.8g}" for x in wp.get("binned_thresholds", []))
        return f"{name}|binned|{edges}|{vals}|{wp['pt_min']:.8g}|{wp['pt_max']:.8g}|{wp.get('max_score', 1.0):.8g}"
    return (
        f"{name}|linear|{wp['intercept']:.8g}|{wp.get('slope', 0.0):.8g}|"
        f"{wp['pt_min']:.8g}|{wp['pt_max']:.8g}|{wp.get('max_score', 1.0):.8g}"
    )


def write_working_point_outputs(outdir: Path, manifest, target_label=None):
    import csv

    target = manifest.get("target_signal_efficiency", 0.80)
    label = target_label or f"target{int(round(100 * target)):02d}"
    json_path = outdir / f"bdt_working_points_{label}.json"
    yaml_path = outdir / f"bdt_working_points_{label}.yaml"
    fragment_path = outdir / f"bdt_working_points_{label}_runtime_fragment.yaml"
    csv_path = outdir / f"bdt_working_points_{label}.csv"
    json_path.write_text(json.dumps(json_ready(manifest), indent=2, sort_keys=True) + "\n")

    lines = [
        "schema: RECOILJETS_AUAU_BDT_WORKING_POINTS_V1",
        f"target_signal_efficiency: {target:.6g}",
        "products:",
    ]
    for product, wp in manifest.get("products", {}).items():
        lines.append(f"  {product}:")
        lines.append(f"    status: {wp.get('status', 'unknown')}")
        if wp.get("status") == "ready":
            lines.append(f"    mode: {wp.get('mode', 'linear')}")
            lines.append(f"    intercept: {wp.get('intercept', math.nan):.10g}")
            lines.append(f"    slope: {wp.get('slope', 0.0):.10g}")
            lines.append(f"    pt_min: {wp.get('pt_min', math.nan):.10g}")
            lines.append(f"    pt_max: {wp.get('pt_max', math.nan):.10g}")
            lines.append(f"    signal_efficiency: {wp.get('inclusive', {}).get('signal_efficiency', math.nan):.10g}")
            lines.append(f"    background_fake_rate: {wp.get('inclusive', {}).get('background_fake_rate', math.nan):.10g}")
            if wp.get("mode") == "grid2d":
                lines.append("    pt_edges: [" + ", ".join(f"{x:.10g}" for x in wp.get("pt_edges", [])) + "]")
                lines.append("    cent_edges: [" + ", ".join(f"{x:.10g}" for x in wp.get("cent_edges", [])) + "]")
                fq = wp.get("fit_quality", {})
                lines.append(f"    max_abs_cell_efficiency_error: {fq.get('max_abs_cell_efficiency_error', math.nan):.10g}")
        else:
            lines.append(f"    reason: {wp.get('reason', 'unknown')}")
    yaml_path.write_text("\n".join(lines) + "\n")

    entries = [runtime_entry(product, wp) for product, wp in manifest.get("products", {}).items()]
    entries = [x for x in entries if x]
    quoted = ", ".join(json.dumps(x) for x in entries)
    fragment_path.write_text(
        "# Product-keyed entries. Use auau_tight_bdt_wp_product_map or the pipeline\n"
        "# config generator to map products onto RecoilJets tight-mode names.\n"
        f"auau_tight_bdt_product_working_point_entries: [{quoted}]\n"
    )

    with csv_path.open("w", newline="") as f:
        fields = [
            "product", "status", "mode", "intercept", "slope", "pt_min", "pt_max",
            "signal_efficiency", "background_fake_rate",
            "max_abs_cell_efficiency_error", "min_cell_signal_entries", "min_cell_background_entries",
        ]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for product, wp in manifest.get("products", {}).items():
            inc = wp.get("inclusive", {})
            writer.writerow({
                "product": product,
                "status": wp.get("status", "unknown"),
                "mode": wp.get("mode", ""),
                "intercept": wp.get("intercept", ""),
                "slope": wp.get("slope", ""),
                "pt_min": wp.get("pt_min", ""),
                "pt_max": wp.get("pt_max", ""),
                "signal_efficiency": inc.get("signal_efficiency", ""),
                "background_fake_rate": inc.get("background_fake_rate", ""),
                "max_abs_cell_efficiency_error": wp.get("fit_quality", {}).get("max_abs_cell_efficiency_error", ""),
                "min_cell_signal_entries": wp.get("fit_quality", {}).get("min_cell_signal_entries", ""),
                "min_cell_background_entries": wp.get("fit_quality", {}).get("min_cell_background_entries", ""),
            })
    return {"json": json_path, "yaml": yaml_path, "runtime_fragment": fragment_path, "csv": csv_path}


def make_working_point_diagnostic(outdir: Path, manifest, target_label=None):
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception:
        return None

    target = manifest.get("target_signal_efficiency", 0.80)
    label = target_label or f"target{int(round(100 * target)):02d}"
    ready = [(p, wp) for p, wp in manifest.get("products", {}).items() if wp.get("status") == "ready"]
    if not ready:
        return None
    if any(wp.get("mode") == "grid2d" for _, wp in ready):
        n = min(len(ready), 8)
        fig, axes = plt.subplots(n, 2, figsize=(10.5, max(2.4 * n, 3.2)), dpi=160, squeeze=False)
        for row, (product, wp) in enumerate(ready[:n]):
            thresholds = np.array(wp.get("grid_thresholds", []), dtype=float)
            eff = np.array(wp.get("grid_signal_efficiency", []), dtype=float)
            pt_edges = np.array(wp.get("pt_edges", []), dtype=float)
            cent_edges = np.array(wp.get("cent_edges", []), dtype=float)
            if thresholds.size == 0 or pt_edges.size < 2 or cent_edges.size < 2:
                continue
            extent = [pt_edges[0], pt_edges[-1], cent_edges[-1], cent_edges[0]]
            ax = axes[row][0]
            im = ax.imshow(thresholds, aspect="auto", interpolation="nearest", extent=extent, vmin=np.nanmin(thresholds), vmax=np.nanmax(thresholds), cmap="viridis")
            ax.set_title(product_label(product), fontsize=10)
            ax.set_ylabel("Centrality [%]")
            if row == n - 1:
                ax.set_xlabel(r"Cluster $E_T$ [GeV]")
            cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.015)
            cb.set_label("BDT cut", fontsize=8)

            ax = axes[row][1]
            delta = eff - float(target)
            vmax = max(0.02, float(np.nanmax(np.abs(delta))) if np.isfinite(delta).any() else 0.02)
            im = ax.imshow(delta, aspect="auto", interpolation="nearest", extent=extent, vmin=-vmax, vmax=vmax, cmap="coolwarm")
            ax.set_title("Achieved signal efficiency minus target", fontsize=10)
            if row == n - 1:
                ax.set_xlabel(r"Cluster $E_T$ [GeV]")
            ax.set_yticklabels([])
            cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.015)
            cb.set_label(r"$\Delta\epsilon_{sig}$", fontsize=8)
        fig.suptitle(f"Validation-derived local BDT cuts for {100*target:.0f}% signal efficiency", y=0.995)
        fig.tight_layout(rect=[0, 0, 1, 0.985])
        path = outdir / f"bdt_working_point_{label}_diagnostics.png"
        fig.savefig(path)
        plt.close(fig)
        return path
    # Keep this diagnostic compact. It is a provenance plot, not a slide plot.
    ready = ready[:18]
    fig, ax = plt.subplots(figsize=(9.0, 5.5), dpi=160)
    for product, wp in ready:
        bins = wp.get("bins", [])
        if bins:
            x = np.array([b["pt_center"] for b in bins], dtype=float)
            y = np.array([b["threshold"] for b in bins], dtype=float)
            ax.plot(x, y, marker="o", linewidth=1.3, markersize=3.5, label=product_label(product))
        if wp.get("mode") == "linear" and math.isfinite(wp.get("pt_min", math.nan)) and math.isfinite(wp.get("pt_max", math.nan)):
            xs = np.linspace(wp["pt_min"], wp["pt_max"], 50)
            ax.plot(xs, wp["intercept"] + wp.get("slope", 0.0) * xs, linewidth=1.0, alpha=0.7)
    ax.set_xlabel(r"Cluster $E_T$ [GeV]")
    ax.set_ylabel(f"BDT score threshold for {100*target:.0f}% signal efficiency")
    ax.set_title("Validation-derived BDT working points")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=7, ncol=2)
    fig.tight_layout()
    path = outdir / f"bdt_working_point_{label}_diagnostics.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def _target_row(thresholds, target):
    rows = thresholds.get("by_signal_efficiency", []) if isinstance(thresholds, dict) else []
    if not rows:
        return None
    return min(rows, key=lambda r: abs(float(r.get("target_signal_efficiency", -999.0)) - target))


def _local_score_cache_manifest(report_dir: Path) -> Path | None:
    local_dir = report_dir / "score_caches"
    if local_dir.is_dir():
        paths = sorted(local_dir.glob("score_cache_*.npz"))
        if paths:
            manifest = report_dir / "score_caches.local.list"
            manifest.write_text("\n".join(str(p) for p in paths) + "\n")
            return manifest
    manifest = report_dir / "score_caches.list"
    if manifest.is_file():
        return manifest
    return None


def derive_centpt_working_points_from_caches(report_dir: Path, target=0.80, pt_edges=None, cent_edges=None):
    import numpy as np

    manifest_path = _local_score_cache_manifest(report_dir)
    if manifest_path is None:
        raise SystemExit(f"Cannot derive centrality x E_T working points; missing score caches under {report_dir}")
    frame, scores, _counts = load_score_caches(manifest_path)
    y = frame["is_signal"].astype("int32")
    et = frame["cluster_Et"].astype("float64")
    cent = frame["centrality"].astype("float64")
    pt_edges = list(pt_edges or ([lo for lo, _ in PT_BINS] + [PT_BINS[-1][1]]))
    cent_edges = list(cent_edges or ([lo for lo, _ in CENT_BINS] + [CENT_BINS[-1][1]]))
    products = {}

    for product, score in scores.items():
        eligible_base = np.isin(y, [0, 1]) & np.isfinite(score) & np.isfinite(et) & np.isfinite(cent)
        eligible_base &= (cent >= cent_edges[0]) & (cent < cent_edges[-1])
        valid_pt_indices = []
        for ip, (plo, phi) in enumerate(zip(pt_edges[:-1], pt_edges[1:])):
            pmask = eligible_base & (et >= plo) & (et < phi)
            if np.any(pmask & (y == 1)) and np.any(pmask & (y == 0)):
                valid_pt_indices.append(ip)
        if not valid_pt_indices:
            products[product] = {
                "status": "skipped",
                "reason": "no E_T bins with both signal and background scores for target working point",
            }
            continue
        first_pt = min(valid_pt_indices)
        last_pt = max(valid_pt_indices)
        product_pt_edges = list(pt_edges[first_pt:last_pt + 2])
        eligible = eligible_base & (et >= product_pt_edges[0]) & (et < product_pt_edges[-1])
        inclusive = threshold_for_signal_efficiency(y[eligible], score[eligible], target)
        if inclusive is None:
            products[product] = {
                "status": "skipped",
                "reason": "insufficient signal/background scores for centrality x E_T target working point",
            }
            continue

        pt_only = {}
        for ip, (plo, phi) in enumerate(zip(product_pt_edges[:-1], product_pt_edges[1:])):
            mask = eligible & (et >= plo) & (et < phi)
            pt_only[ip] = threshold_for_signal_efficiency(y[mask], score[mask], target)
        cent_only = {}
        for ic, (clo, chi) in enumerate(zip(cent_edges[:-1], cent_edges[1:])):
            mask = eligible & (cent >= clo) & (cent < chi)
            cent_only[ic] = threshold_for_signal_efficiency(y[mask], score[mask], target)

        thresholds = []
        eff_grid = []
        fake_grid = []
        cells = []
        min_sig = None
        min_bkg = None
        for ic, (clo, chi) in enumerate(zip(cent_edges[:-1], cent_edges[1:])):
            thr_row = []
            eff_row = []
            fake_row = []
            for ip, (plo, phi) in enumerate(zip(product_pt_edges[:-1], product_pt_edges[1:])):
                mask = eligible & (cent >= clo) & (cent < chi) & (et >= plo) & (et < phi)
                item = threshold_for_signal_efficiency(y[mask], score[mask], target)
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
                sig_mask = mask & (y == 1)
                bkg_mask = mask & (y == 0)
                sig_scores = score[sig_mask]
                bkg_scores = score[bkg_mask]
                sig_eff = float(np.mean(sig_scores > threshold)) if sig_scores.size else math.nan
                bkg_fake = float(np.mean(bkg_scores > threshold)) if bkg_scores.size else math.nan
                sig_n = int(sig_scores.size)
                bkg_n = int(bkg_scores.size)
                min_sig = sig_n if min_sig is None else min(min_sig, sig_n)
                min_bkg = bkg_n if min_bkg is None else min(min_bkg, bkg_n)
                thr_row.append(threshold)
                eff_row.append(sig_eff)
                fake_row.append(bkg_fake)
                cells.append({
                    "centrality_min": float(clo),
                    "centrality_max": float(chi),
                    "pt_min": float(plo),
                    "pt_max": float(phi),
                    "threshold": threshold,
                    "signal_efficiency": sig_eff,
                    "background_fake_rate": bkg_fake,
                    "signal_entries": sig_n,
                    "background_entries": bkg_n,
                    "source": source,
                })
            thresholds.append(thr_row)
            eff_grid.append(eff_row)
            fake_grid.append(fake_row)

        cell_eff = np.array([c["signal_efficiency"] for c in cells], dtype=float)
        finite_cell_eff = cell_eff[np.isfinite(cell_eff)]
        # A plane fit is recorded as a smoothness diagnostic. The runtime uses
        # the local grid thresholds so the requested efficiency is preserved
        # cell-by-cell instead of being washed out by a global fit.
        fit_points = []
        fit_values = []
        fit_weights = []
        for c in cells:
            if not math.isfinite(c["threshold"]) or c["signal_entries"] <= 0:
                continue
            fit_points.append([1.0, 0.5 * (c["pt_min"] + c["pt_max"]), 0.5 * (c["centrality_min"] + c["centrality_max"])])
            fit_values.append(c["threshold"])
            fit_weights.append(math.sqrt(max(1, c["signal_entries"])))
        plane = {"intercept": math.nan, "pt_slope": math.nan, "centrality_slope": math.nan, "max_abs_residual": math.nan}
        if len(fit_points) >= 3:
            x = np.array(fit_points, dtype=float)
            t = np.array(fit_values, dtype=float)
            w = np.array(fit_weights, dtype=float)
            coeff, *_ = np.linalg.lstsq(x * w[:, None], t * w, rcond=None)
            pred = x @ coeff
            plane = {
                "intercept": float(coeff[0]),
                "pt_slope": float(coeff[1]),
                "centrality_slope": float(coeff[2]),
                "max_abs_residual": float(np.max(np.abs(t - pred))) if t.size else math.nan,
            }

        products[product] = {
            "status": "ready",
            "target_signal_efficiency": float(target),
            "mode": "grid2d",
            "pt_min": float(product_pt_edges[0]),
            "pt_max": float(product_pt_edges[-1]),
            "max_score": 1.0,
            "pt_edges": [float(x) for x in product_pt_edges],
            "cent_edges": [float(x) for x in cent_edges],
            "grid_thresholds": thresholds,
            "grid_signal_efficiency": eff_grid,
            "grid_background_fake_rate": fake_grid,
            "cells": cells,
            "inclusive": inclusive,
            "fit_quality": {
                "max_abs_cell_efficiency_error": float(np.max(np.abs(finite_cell_eff - target))) if finite_cell_eff.size else math.nan,
                "mean_abs_cell_efficiency_error": float(np.mean(np.abs(finite_cell_eff - target))) if finite_cell_eff.size else math.nan,
                "min_cell_signal_entries": int(min_sig or 0),
                "min_cell_background_entries": int(min_bkg or 0),
                "plane_fit": plane,
            },
        }
    return {
        "schema": "RECOILJETS_AUAU_BDT_WORKING_POINTS_V2",
        "target_signal_efficiency": float(target),
        "working_point_mode": "centpt",
        "source": str(report_dir),
        "products": products,
    }


def derive_working_points_from_report(report_dir: Path, target=0.80, wp_mode="centpt", pt_edges=None, cent_edges=None):
    import numpy as np

    existing = report_dir / f"bdt_working_points_target{int(round(100 * target)):02d}.json"
    if wp_mode == "centpt":
        manifest = derive_centpt_working_points_from_caches(report_dir, target, pt_edges=pt_edges, cent_edges=cent_edges)
        paths = write_working_point_outputs(report_dir, manifest)
        make_working_point_diagnostic(report_dir, manifest)
        return paths
    if existing.is_file():
        manifest = json.loads(existing.read_text())
        paths = write_working_point_outputs(report_dir, manifest)
        make_working_point_diagnostic(report_dir, manifest)
        return paths

    diag_path = report_dir / "validation_deep_diagnostics.json"
    if not diag_path.is_file():
        raise SystemExit(f"Cannot derive working points; missing {diag_path}")
    diagnostics = json.loads(diag_path.read_text())
    manifest = {
        "schema": "RECOILJETS_AUAU_BDT_WORKING_POINTS_V1",
        "target_signal_efficiency": float(target),
        "source": str(report_dir),
        "products": {},
    }
    pt_bins = diagnostics.get("pt_bins", [])
    for product, product_diag in diagnostics.get("products", {}).items():
        inclusive = _target_row(product_diag.get("thresholds_inclusive", {}), target)
        if not inclusive:
            manifest["products"][product] = {"status": "skipped", "reason": "missing inclusive threshold"}
            continue
        bins = []
        for b in pt_bins:
            key = b.get("key")
            row = _target_row(product_diag.get("auc_by_pt", {}).get(key, {}).get("thresholds", {}), target)
            if not row:
                continue
            bins.append({
                "pt_min": float(b["low"]),
                "pt_max": float(b["high"]),
                "pt_center": 0.5 * (float(b["low"]) + float(b["high"])),
                "threshold": float(row["threshold"]),
                "signal_efficiency": float(row.get("signal_efficiency", math.nan)),
                "background_fake_rate": float(row.get("background_fake_rate", math.nan)),
                "signal_entries": int(product_diag.get("auc_by_pt", {}).get(key, {}).get("signal_entries", 0)),
                "background_entries": int(product_diag.get("auc_by_pt", {}).get(key, {}).get("background_entries", 0)),
            })
        mode = "linear"
        intercept = float(inclusive["threshold"])
        slope = 0.0
        if len(bins) >= 2:
            x = np.array([b["pt_center"] for b in bins], dtype=float)
            t = np.array([b["threshold"] for b in bins], dtype=float)
            w = np.sqrt(np.array([max(1, b["signal_entries"]) for b in bins], dtype=float))
            coeff = np.polyfit(x, t, 1, w=w)
            slope = float(coeff[0])
            intercept = float(coeff[1])
            residual = t - (intercept + slope * x)
            max_resid = float(np.max(np.abs(residual))) if residual.size else math.nan
            if math.isfinite(max_resid) and max_resid > 0.12:
                mode = "binned"
        pt_min = float(bins[0]["pt_min"]) if bins else math.nan
        pt_max = float(bins[-1]["pt_max"]) if bins else math.nan
        wp = {
            "status": "ready",
            "target_signal_efficiency": float(target),
            "mode": mode,
            "intercept": intercept,
            "slope": slope,
            "pt_min": pt_min,
            "pt_max": pt_max,
            "max_score": 1.0,
            "inclusive": inclusive,
            "bins": bins,
        }
        if mode == "binned":
            wp["binned_edges"] = [float(bins[0]["pt_min"])] + [float(b["pt_max"]) for b in bins]
            wp["binned_thresholds"] = [float(b["threshold"]) for b in bins]
        manifest["products"][product] = wp
    paths = write_working_point_outputs(report_dir, manifest)
    make_working_point_diagnostic(report_dir, manifest)
    return paths


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
        label = product_label(product)
        roc = metrics["products"].get(product, {}).get("_roc")
        if roc is not None:
            fpr, tpr, thresholds = roc
            rejection = 1.0 - fpr
            roc_payload["products"][product] = {
                "label": label,
                "color": product_color(product),
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
        product_hist = {"label": label, "color": product_color(product), "by_centrality": {}}
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
                {"label": label, "color": product_color(product)},
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
                {"label": label, "color": product_color(product)},
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

    all_products = list(scores)
    ranked_products = sorted(
        all_products,
        key=lambda p: metrics["products"].get(p, {}).get("auc_inclusive", float("-inf")),
        reverse=True,
    )
    max_detail = int(os.environ.get("RJ_AUAU_TIGHT_BDT_VALIDATE_MAX_DETAIL_PLOTS", "12"))
    detail_products = ranked_products[:max(1, max_detail)]

    fig, ax = plt.subplots(figsize=(7.0, 6.0))
    for product in detail_products:
        score = scores[product]
        auc = metrics["products"][product]["auc_inclusive"]
        roc = metrics["products"][product].get("_roc")
        if roc is None:
            continue
        fpr, tpr, _ = roc
        label = product_label(product)
        ax.plot(fpr, tpr, lw=2, color=product_color(product), label=f"{label} AUC={auc:.3f}")
    ax.plot([0, 1], [0, 1], "k--", lw=1)
    ax.set_xlabel("Background fake rate")
    ax.set_ylabel("Signal efficiency")
    ax.set_title("AuAu tight-BDT ROC on embedded simulation")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(outdir / "roc_inclusive.png", dpi=170)
    plt.close(fig)

    for product in detail_products:
        score = scores[product]
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
            ax.set_title(f"{product_label(product)} ROC by centrality")
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
            ax.set_title(f"{product_label(product)} ROC by photon E$_T$")
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
        ax.set_title(f"{product_label(product)}: signal/background score separation")
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
        fig.suptitle(f"{product_label(product)} score separation by centrality", y=1.02)
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
            ax.set_title(f"{product_label(product)}: {cent_label(lo, hi)} score separation")
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
            ax.set_title(f"{product_label(product)}: BDT score vs {xlabel}")
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
            ax.set_title(f"{product_label(product)} AUC by centrality and $E_T$")
            for ic in range(len(CENT_BINS)):
                for ip in range(len(PT_BINS)):
                    if np.isfinite(heat[ic, ip]):
                        ax.text(ip, ic, f"{heat[ic, ip]:.3f}", ha="center", va="center", color="white", fontsize=8)
            fig.colorbar(im, ax=ax, label="AUC")
            fig.tight_layout()
            fig.savefig(outdir / f"auc_heatmap_cent_pt_{product}.png", dpi=170)
            plt.close(fig)

    products = ranked_products[:25]
    x = np.arange(len(products))
    aucs = [metrics["products"][p]["auc_inclusive"] for p in products]
    fig, ax = plt.subplots(figsize=(max(7.0, 0.35 * len(products)), 4.8))
    ax.bar(x, aucs, color=[product_color(p) for p in products])
    ax.set_xticks(x, [product_label(p) for p in products], rotation=20, ha="right")
    ax.set_ylim(0.45, 1.0)
    ax.set_ylabel("Inclusive AUC")
    ax.set_title("AuAu tight-BDT inclusive AUC")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "auc_summary.png", dpi=170)
    plt.close(fig)


def build_metrics(frame, scores, args, product_specs=None):
    import numpy as np

    y = frame["is_signal"].astype("int32")
    et = frame["cluster_Et"].astype("float64")
    cent = frame["centrality"].astype("float64")
    eiso = frame["reco_eiso"].astype("float64")
    metrics = {"products": {}}
    for product, score in scores.items():
        eligible = np.isin(y, [0, 1])
        if product_specs is not None and product in product_specs:
            eligible &= product_eligible_mask(frame, product_specs[product])
        elif product == "centDepBDTs":
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
        by_pt_cent = {}
        for plo, phi in PT_BINS:
            for clo, chi in CENT_BINS:
                mask = eligible & (et >= plo) & (et < phi) & (cent >= clo) & (cent < chi)
                a, _ = auc_for(y[mask], score[mask])
                by_pt_cent[f"pt_{plo:g}_{phi:g}_cent_{clo:g}_{chi:g}"] = a
        corr_mask = finite & np.isfinite(eiso)
        corr = math.nan
        if corr_mask.sum() >= 3:
            corr = float(np.corrcoef(score[corr_mask], eiso[corr_mask])[0, 1])
        metrics["products"][product] = {
            "auc_inclusive": auc,
            "_roc": roc,
            "auc_by_pt": by_pt,
            "auc_by_centrality": by_cent,
            "auc_by_pt_centrality": by_pt_cent,
            "eligible_entries": int(eligible.sum()),
            "finite_scores": int(finite.sum()),
            "finite_score_fraction": float(finite.sum() / max(1, eligible.sum())),
            "score_eiso_pearson": corr,
            "signal_score_mean": float(np.nanmean(score[finite & (y == 1)])) if (finite & (y == 1)).any() else math.nan,
            "background_score_mean": float(np.nanmean(score[finite & (y == 0)])) if (finite & (y == 0)).any() else math.nan,
        }
    return metrics


def ensure_product_metrics(metrics, products):
    defaults = {
        "auc_inclusive": math.nan,
        "_roc": None,
        "auc_by_pt": {f"{lo:g}_{hi:g}": math.nan for lo, hi in PT_BINS},
        "auc_by_centrality": {f"{lo:g}_{hi:g}": math.nan for lo, hi in CENT_BINS},
        "auc_by_pt_centrality": {
            f"pt_{plo:g}_{phi:g}_cent_{clo:g}_{chi:g}": math.nan
            for plo, phi in PT_BINS
            for clo, chi in CENT_BINS
        },
        "eligible_entries": 0,
        "finite_scores": 0,
        "finite_score_fraction": 0.0,
        "score_eiso_pearson": math.nan,
        "signal_score_mean": math.nan,
        "background_score_mean": math.nan,
    }
    metrics.setdefault("products", {})
    for product in products:
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
        if args.model_registry is not None:
            f.write(f"model_registry={args.model_registry}\n")
        f.write(f"total_entries={counts['total_entries']}\n")
        f.write(f"signal_entries={counts['signal_entries']}\n")
        f.write(f"background_entries={counts['background_entries']}\n")
        f.write(f"scored_entries={counts['scored_entries']}\n")
        for product in products:
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


def write_model_rankings(outdir: Path, metrics) -> None:
    import csv
    import numpy as np

    rows = []
    for product, product_metrics in metrics.get("products", {}).items():
        by_pt = [v for v in product_metrics.get("auc_by_pt", {}).values() if isinstance(v, (int, float)) and math.isfinite(v)]
        by_cent = [v for v in product_metrics.get("auc_by_centrality", {}).values() if isinstance(v, (int, float)) and math.isfinite(v)]
        by_pt_cent = [v for v in product_metrics.get("auc_by_pt_centrality", {}).values() if isinstance(v, (int, float)) and math.isfinite(v)]
        rows.append(
            {
                "product": product,
                "label": product_label(product),
                "auc_inclusive": product_metrics.get("auc_inclusive", math.nan),
                "auc_pt_mean": float(np.mean(by_pt)) if by_pt else math.nan,
                "auc_cent_mean": float(np.mean(by_cent)) if by_cent else math.nan,
                "auc_pt_cent_mean": float(np.mean(by_pt_cent)) if by_pt_cent else math.nan,
                "finite_score_fraction": product_metrics.get("finite_score_fraction", math.nan),
                "eligible_entries": product_metrics.get("eligible_entries", 0),
                "signal_score_mean": product_metrics.get("signal_score_mean", math.nan),
                "background_score_mean": product_metrics.get("background_score_mean", math.nan),
            }
        )
    rows.sort(
        key=lambda r: (
            np.nan_to_num(r["auc_inclusive"], nan=-1.0),
            np.nan_to_num(r["auc_pt_cent_mean"], nan=-1.0),
            np.nan_to_num(r["auc_pt_mean"], nan=-1.0),
        ),
        reverse=True,
    )
    (outdir / "validation_model_rankings.json").write_text(json.dumps(json_ready(rows), indent=2, sort_keys=True) + "\n")
    with (outdir / "validation_model_rankings.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]) if rows else ["product"])
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()

    if args.derive_working_points_from_report is not None:
        pt_fallback = [lo for lo, _ in PT_BINS] + [PT_BINS[-1][1]]
        cent_fallback = [lo for lo, _ in CENT_BINS] + [CENT_BINS[-1][1]]
        pt_edges = parse_edge_list(args.working_point_pt_bins, pt_fallback)
        cent_edges = parse_edge_list(args.working_point_cent_bins, cent_fallback)
        paths = derive_working_points_from_report(
            args.derive_working_points_from_report,
            args.target_signal_efficiency,
            wp_mode=args.working_point_mode,
            pt_edges=pt_edges,
            cent_edges=cent_edges,
        )
        for name, path in paths.items():
            print(f"{name}={path}")
        return 0

    if args.source is None or args.model_dir is None:
        raise SystemExit("--source and --model-dir are required unless --derive-working-points-from-report is used")

    require_imports()
    import numpy as np

    outdir = args.outdir or (args.source / "reports" / f"model_validation_{os.environ.get('RJ_AUAU_TIGHT_BDT_VALIDATE_STAMP') or __import__('datetime').datetime.now().strftime('%Y%m%d_%H%M%S')}")
    outdir.mkdir(parents=True, exist_ok=True)
    product_specs = active_product_specs(args.model_dir, args.model_registry)
    product_names = list(product_specs)
    feature_names = sorted({name for spec in product_specs.values() for name in spec["features"]})

    if args.merge_score_caches is not None:
        frame, scores, counts = load_score_caches(args.merge_score_caches)
        product_names = list(scores)
    else:
        paths = read_manifest(args.manifest) if args.manifest is not None else discover_manifest(args.source)
        models = load_models(args.model_dir, product_specs)
        frame, counts = collect_rows(
            paths,
            args.tree,
            args.score_max_rows,
            args.progress_every,
            feature_names=feature_names,
        )
        scores = score_rows(frame, models, product_specs) if counts["scored_entries"] > 0 else {p: np.array([], dtype="float32") for p in product_names}
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

    metrics = build_metrics(frame, scores, args, product_specs) if counts["scored_entries"] > 0 else {"products": {}}
    metrics = ensure_product_metrics(metrics, product_names)

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
    if set(metrics["products"]) != set(product_names):
        status = "CHECK"
        notes.append("not_all_products_scored")

    if not args.no_plots:
        diagnostics = build_deep_diagnostics(frame, scores)
        (outdir / "validation_deep_diagnostics.json").write_text(
            json.dumps(json_ready(diagnostics), indent=2, sort_keys=True) + "\n"
        )
        write_diagnostic_tables(outdir, diagnostics)
        working_points = derive_working_points(frame, scores, product_specs, args.target_signal_efficiency)
        write_working_point_outputs(outdir, working_points)
        make_working_point_diagnostic(outdir, working_points)
        make_plots(outdir, frame, scores, metrics)
        write_curve_exports(outdir, frame, scores, metrics)
    metrics_json = json_ready({"counts": counts, **metrics, "status": status, "notes": notes})
    (outdir / "validation_metrics.json").write_text(json.dumps(metrics_json, indent=2, sort_keys=True) + "\n")
    write_model_rankings(outdir, metrics)
    write_parseable_summary(outdir / "validation_summary.txt", status, args, counts, metrics, notes)
    print((outdir / "validation_summary.txt").read_text(), end="")
    return 0 if status == "READY" else 3


if __name__ == "__main__":
    raise SystemExit(main())
