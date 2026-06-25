#!/usr/bin/env python3
"""Train AuAu photon BDT TMVA files from RecoilJets_AuAu training trees.

The expected input tree is ``AuAuPhotonIDTrainingTree``.  For tight photon-ID,
``is_signal`` labels truth-matched isolated prompt photons from embedded
photon+jet signal samples and non-signal clusters from embedded inclusive jet
background samples.  For NPB, this follows the PPG12 convention:
``npb_label=1`` means a physics-like cluster and ``npb_label=0`` means a
timing-tagged non-physics background cluster from data.  The resulting score is
therefore a physics-like score, so the analysis cut remains
``auau_npb_score > 0.5``.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import os
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable

PPG12_EXACT_WEIGHT_COLUMN = "__ppg12_exact_training_weight"
PPG12_EXACT_ETA_RANGE = (-0.7, 0.7)
PPG12_EXACT_N_BINS = 20
PPG12_EXACT_ET_WEIGHT_CAP = 800.0
PPG12_EXACT_EXPECTED_SAMPLES = (
    "run28_embeddedPhoton12",
    "run28_embeddedPhoton20",
    "run28_embeddedJet12",
    "run28_embeddedJet20",
    "run28_embeddedJet30",
)
PPG12_EXACT_SAMPLE_ALIASES = (
    ("run28_embeddedPhoton12", ("run28_embeddedPhoton12", "embeddedPhoton12", "Photon12")),
    ("run28_embeddedPhoton20", ("run28_embeddedPhoton20", "embeddedPhoton20", "Photon20")),
    ("run28_embeddedJet12", ("run28_embeddedJet12", "embeddedJet12", "Jet12")),
    ("run28_embeddedJet20", ("run28_embeddedJet20", "embeddedJet20", "Jet20")),
    ("run28_embeddedJet30", ("run28_embeddedJet30", "embeddedJet30", "Jet30")),
    ("run28_embeddedJet40", ("run28_embeddedJet40", "embeddedJet40", "Jet40")),
    ("run28_photonjet5", ("run28_photonjet5", "photonjet5", "PhotonJet5")),
    ("run28_photonjet10", ("run28_photonjet10", "photonjet10", "PhotonJet10")),
    ("run28_photonjet20", ("run28_photonjet20", "photonjet20", "PhotonJet20")),
    ("run28_jet8", ("run28_jet8",)),
    ("run28_jet12", ("run28_jet12",)),
    ("run28_jet20", ("run28_jet20",)),
    ("run28_jet30", ("run28_jet30",)),
    ("run28_jet40", ("run28_jet40",)),
)

GLOBAL_EVENT_KEY_COLUMNS = ["source_sample", "input_file_index", "run", "evt"]
GLOBAL_CANDIDATE_KEY_COLUMNS = GLOBAL_EVENT_KEY_COLUMNS + ["input_tree_entry"]
TRAINING_IDENTITY_OPTIONAL_COLUMNS = ["run", "evt", "global_event_key"]
EVENT_QUALITY_DIRECT_COLUMNS = ["centrality", "event_calo_log10_total_energy_plus1"]
EVENT_QUALITY_COMPONENT_COLUMNS = [
    "event_calo_cemc_energy",
    "event_calo_ihcal_energy",
    "event_calo_ohcal_energy",
    "event_calo_total_energy",
]


PPG12_TIGHT_FEATURES = [
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

PPG12_BASE_V1E_FEATURES = [
    "cluster_Et",
    "cluster_weta_cogx",
    "vertexz",
    "cluster_Eta",
    "e11_over_e33",
    "cluster_et1",
    "cluster_et2",
    "cluster_et3",
    "cluster_et4",
]

PPG12_TIGHT_FEATURES_3X3_WIDTHS = [
    "cluster_Et",
    "cluster_weta33_cogx",
    "cluster_wphi33_cogx",
    "vertexz",
    "cluster_Eta",
    "e11_over_e33",
    "cluster_et1",
    "cluster_et2",
    "cluster_et3",
    "cluster_et4",
    "e32_over_e35",
]

PPG12_TIGHT_FEATURES_BASE_AND_3X3_WIDTHS = [
    "cluster_Et",
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "cluster_weta33_cogx",
    "cluster_wphi33_cogx",
    "vertexz",
    "cluster_Eta",
    "e11_over_e33",
    "cluster_et1",
    "cluster_et2",
    "cluster_et3",
    "cluster_et4",
    "e32_over_e35",
]

WIDTH_RATIO_FEATURES = [
    "cluster_weta_over_wphi",
    "cluster_weta33_over_wphi33",
]

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

ISOLATION_DIAGNOSTIC_FEATURES = [
    "reco_eiso_clip30",
    "reco_eiso_over_cluster_Et",
    "reco_eiso_signed_log1p",
]
ISOLATION_R30_DIAGNOSTIC_FEATURES = [
    "reco_eiso_r30_clip30",
    "reco_eiso_r30_over_cluster_Et",
    "reco_eiso_r30_signed_log1p",
]
ISOLATION_R40_DIAGNOSTIC_FEATURES = [
    "reco_eiso_r40_clip30",
    "reco_eiso_r40_over_cluster_Et",
    "reco_eiso_r40_signed_log1p",
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
    "reco_eiso_r30_clip30": ["reco_eiso_r30"],
    "reco_eiso_r30_over_cluster_Et": ["reco_eiso_r30", "cluster_Et"],
    "reco_eiso_r30_signed_log1p": ["reco_eiso_r30"],
    "reco_eiso_r40_clip30": ["reco_eiso_r40"],
    "reco_eiso_r40_over_cluster_Et": ["reco_eiso_r40", "cluster_Et"],
    "reco_eiso_r40_signed_log1p": ["reco_eiso_r40"],
}

TIGHT_MODES = [
    "legacy",
    "ppg12BaseV1E",
    "centINDcontrol",
    "centAsFeat",
    "centAsFeatMinOpt",
    "centAsFeat3x3",
    "centAsFeatBase3x3",
    "centAsFeatWidthRatios",
    "isoDiagnosticFull",
    "centDepBDTs",
]
TMVA_EXPORT_LOCK = threading.Lock()


def global_sixpack_noiso_features() -> list[str]:
    features: list[str] = []
    for feature in (
        PPG12_TIGHT_FEATURES_BASE_AND_3X3_WIDTHS
        + EXTENDED_SHOWER_FEATURES
        + WIDTH_RATIO_FEATURES
        + ["centrality"]
    ):
        if feature not in features:
            features.append(feature)
    return features


def global_sixpack_iso_features() -> list[str]:
    features = global_sixpack_noiso_features()
    for feature in ISOLATION_DIAGNOSTIC_FEATURES:
        if feature not in features:
            features.append(feature)
    return features


def global_sixpack_eiso_r30_features() -> list[str]:
    features = global_sixpack_noiso_features()
    if "reco_eiso_r30" not in features:
        features.append("reco_eiso_r30")
    return features


def global_sixpack_eiso_r40_features() -> list[str]:
    features = global_sixpack_noiso_features()
    if "reco_eiso_r40" not in features:
        features.append("reco_eiso_r40")
    return features


def global_sixpack_eiso_r30r40_features() -> list[str]:
    features = global_sixpack_noiso_features()
    for feature in ("reco_eiso_r30", "reco_eiso_r40"):
        if feature not in features:
            features.append(feature)
    return features


def tight_mode_features(mode: str, override: list[str] | None) -> list[str]:
    if override is not None:
        features = list(override)
    elif mode == "ppg12BaseV1E":
        features = list(PPG12_BASE_V1E_FEATURES)
    elif mode == "centAsFeat3x3":
        features = list(PPG12_TIGHT_FEATURES_3X3_WIDTHS)
    elif mode == "centAsFeatBase3x3":
        features = list(PPG12_TIGHT_FEATURES_BASE_AND_3X3_WIDTHS)
    elif mode == "centAsFeatWidthRatios":
        features = list(PPG12_TIGHT_FEATURES) + [
            "cluster_weta33_cogx",
            "cluster_wphi33_cogx",
            *WIDTH_RATIO_FEATURES,
        ]
    elif mode == "isoDiagnosticFull":
        features = diagnostic_isolation_feature_family(include_centrality=True)
    else:
        features = list(PPG12_TIGHT_FEATURES)
    if mode in ("centAsFeat", "centAsFeatMinOpt", "centAsFeat3x3", "centAsFeatBase3x3", "centAsFeatWidthRatios", "isoDiagnosticFull") and "centrality" not in features and "cent" not in features:
        features.append("centrality")
    return features


def diagnostic_isolation_feature_family(include_centrality: bool = True) -> list[str]:
    features: list[str] = []
    for feature in (
        PPG12_TIGHT_FEATURES_BASE_AND_3X3_WIDTHS
        + EXTENDED_SHOWER_FEATURES
        + WIDTH_RATIO_FEATURES
        + ISOLATION_DIAGNOSTIC_FEATURES
        + (["centrality"] if include_centrality else [])
    ):
        if feature not in features:
            features.append(feature)
    return features

PPG12_NPB_FEATURES = [
    "cluster_Et",
    "cluster_Eta",
    "vertexz",
    "e11_over_e33",
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
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "cluster_et1",
    "cluster_et2",
    "cluster_et3",
    "cluster_et4",
    "cluster_w32",
    "cluster_w52",
    "cluster_w72",
]

EISO_CONE_ABLATION_PT_EDGES = [15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 30.0, 35.0]
EISO_CONE_ABLATION_COARSE_CENT_BINS = [(0.0, 20.0), (20.0, 50.0), (50.0, 80.0)]
EISO_CONE_ABLATION_FINE_CENT_BINS = [
    (0.0, 10.0),
    (10.0, 20.0),
    (20.0, 30.0),
    (30.0, 40.0),
    (40.0, 50.0),
    (50.0, 60.0),
    (60.0, 80.0),
]
EISO_CONE_ABLATION_EXPECTED_COUNTS = {
    "globalEtCent1535_bdt_eisoR30_ptCent3": 24,
    "globalEtCent1535_bdt_eisoR30_ptCent7": 56,
    "globalEtCent1535_bdt_eisoR40_ptCent3": 24,
    "globalEtCent1535_bdt_eisoR40_ptCent7": 56,
    "globalEtCent1535_bdt_eisoR30R40_ptCent3": 24,
    "globalEtCent1535_bdt_eisoR30R40_ptCent7": 56,
}


def parse_cent_bins(text: str) -> list[tuple[float, float]]:
    bins: list[tuple[float, float]] = []
    if not text:
        return bins
    for item in text.split(","):
        lo_s, hi_s = item.split(":", 1)
        bins.append((float(lo_s), float(hi_s)))
    return bins


def parse_range_list(text: str) -> list[tuple[float, float]]:
    ranges: list[tuple[float, float]] = []
    if not text:
        return ranges
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        if ":" not in item:
            raise SystemExit(f"Range items must be lo:hi, got {item!r}")
        lo_s, hi_s = item.split(":", 1)
        lo = float(lo_s)
        hi = float(hi_s)
        if hi <= lo:
            raise SystemExit(f"Range upper edge must exceed lower edge: {item!r}")
        ranges.append((lo, hi))
    return ranges


def pt_window_tag(lo: float, hi: float) -> str:
    return f"pt{int(round(lo)):g}to{int(round(hi)):g}".replace(".", "p")


def cent_tag(lo: float, hi: float) -> str:
    return f"cent_{int(round(lo)):03d}_{int(round(hi)):03d}"


def parse_float_edges(text: str) -> list[float]:
    vals = [float(item.strip()) for item in text.split(",") if item.strip()]
    if len(vals) < 2:
        raise SystemExit(f"Need at least two bin edges, got: {text}")
    for lo, hi in zip(vals, vals[1:]):
        if hi <= lo:
            raise SystemExit(f"Bin edges must increase strictly: {text}")
    return vals


def bins_from_edges(edges: list[float]) -> list[tuple[float, float]]:
    return list(zip(edges[:-1], edges[1:]))


def pt_tag(lo: float, hi: float) -> str:
    return f"pt_{int(round(lo)):03d}_{int(round(hi)):03d}"


def range_label(rng: tuple[float, float] | None) -> str:
    if rng is None:
        return "all"
    lo, hi = rng
    return f"{lo:g}_{hi:g}".replace(".", "p")


def apply_single_pt_range(frame, text: str):
    ranges = parse_range_list(text)
    if not ranges:
        return frame, None
    if len(ranges) != 1:
        raise SystemExit("--pt-range accepts exactly one lo:hi window")
    lo, hi = ranges[0]
    if "cluster_Et" not in frame.columns:
        raise SystemExit("--pt-range requires cluster_Et in the input training tree")
    before = len(frame)
    selected = frame[(frame["cluster_Et"] >= lo) & (frame["cluster_Et"] < hi)].copy()
    report = {"lo": lo, "hi": hi, "rows_before": before, "rows_after": len(selected)}
    return selected, report


def expand_required_columns(columns: Iterable[str]) -> list[str]:
    expanded: set[str] = set()
    for col in columns:
        deps = DERIVED_FEATURE_DEPS.get(col)
        if deps:
            expanded.update(deps)
        else:
            expanded.add(col)
    return sorted(expanded)


def add_derived_features(frame):
    import numpy as np

    eps = 1.0e-6

    def safe_ratio(num: str, den: str):
        n = frame[num].to_numpy(dtype="float64")
        d = frame[den].to_numpy(dtype="float64")
        out = np.full(len(frame), np.nan, dtype="float64")
        good = np.isfinite(n) & np.isfinite(d) & (np.abs(d) > 1.0e-9)
        out[good] = n[good] / d[good]
        return out

    def safe_log_ratio(num: str, den: str):
        n = frame[num].to_numpy(dtype="float64")
        d = frame[den].to_numpy(dtype="float64")
        out = np.full(len(frame), np.nan, dtype="float64")
        good = np.isfinite(n) & np.isfinite(d) & (n > eps) & (d > eps)
        out[good] = np.log(n[good] / d[good])
        return out

    def safe_logit(name: str):
        x = frame[name].to_numpy(dtype="float64")
        out = np.full(len(frame), np.nan, dtype="float64")
        good = np.isfinite(x)
        clipped = np.clip(x[good], eps, 1.0 - eps)
        out[good] = np.log(clipped / (1.0 - clipped))
        return out

    if "cluster_weta_over_wphi" not in frame.columns and {"cluster_weta_cogx", "cluster_wphi_cogx"}.issubset(frame.columns):
        frame["cluster_weta_over_wphi"] = safe_ratio("cluster_weta_cogx", "cluster_wphi_cogx")
    if "cluster_weta33_over_wphi33" not in frame.columns and {"cluster_weta33_cogx", "cluster_wphi33_cogx"}.issubset(frame.columns):
        frame["cluster_weta33_over_wphi33"] = safe_ratio("cluster_weta33_cogx", "cluster_wphi33_cogx")
    if "shape_tail_e37_e53_logratio" not in frame.columns and {"e22_over_e37", "e22_over_e53"}.issubset(frame.columns):
        frame["shape_tail_e37_e53_logratio"] = safe_log_ratio("e22_over_e37", "e22_over_e53")
    if "shape_long_tail_e17_e71_logratio" not in frame.columns and {"e11_over_e17", "e11_over_e71"}.issubset(frame.columns):
        frame["shape_long_tail_e17_e71_logratio"] = safe_log_ratio("e11_over_e17", "e11_over_e71")
    if (
        "shape_width_core_full_logratio" not in frame.columns
        and {"cluster_weta_cogx", "cluster_wphi_cogx", "cluster_weta33_cogx", "cluster_wphi33_cogx"}.issubset(frame.columns)
    ):
        core = safe_ratio("cluster_weta33_cogx", "cluster_wphi33_cogx")
        full = safe_ratio("cluster_weta_cogx", "cluster_wphi_cogx")
        out = np.full(len(frame), np.nan, dtype="float64")
        good = np.isfinite(core) & np.isfinite(full) & (core > eps) & (full > eps)
        out[good] = np.log(core[good] / full[good])
        frame["shape_width_core_full_logratio"] = out
    if "shape_core_tail_tension" not in frame.columns and {"e11_over_e33", "e22_over_e37", "e22_over_e53"}.issubset(frame.columns):
        frame["shape_core_tail_tension"] = safe_logit("e11_over_e33") - 0.5 * (
            safe_logit("e22_over_e37") + safe_logit("e22_over_e53")
        )
    if "shape_compactness_gradient" not in frame.columns and {"e11_over_e33", "e32_over_e35"}.issubset(frame.columns):
        frame["shape_compactness_gradient"] = safe_logit("e11_over_e33") - safe_logit("e32_over_e35")

    if (
        {"shape_template_diag_chi2", "shape_template_max_abs_z"}.difference(frame.columns)
        and {"is_signal", "cluster_Et", "centrality", *SHAPE_TEMPLATE_BASIS}.issubset(frame.columns)
    ):
        basis = SHAPE_TEMPLATE_BASIS
        values = frame[basis].to_numpy(dtype="float64")
        is_signal = frame["is_signal"].to_numpy(dtype="int32") == 1
        et = frame["cluster_Et"].to_numpy(dtype="float64")
        cent = frame["centrality"].to_numpy(dtype="float64")
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

        chi2 = np.full(len(frame), np.nan, dtype="float64")
        max_abs = np.full(len(frame), np.nan, dtype="float64")
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
        frame["shape_template_diag_chi2"] = chi2
        frame["shape_template_max_abs_z"] = max_abs

    cols = set(frame.columns)
    if "reco_eiso_clip30" not in cols and "reco_eiso" in cols:
        reco_eiso = frame["reco_eiso"].to_numpy(dtype="float64")
        frame["reco_eiso_clip30"] = np.where(
            np.isfinite(reco_eiso) & (np.abs(reco_eiso) < 1.0e8),
            np.clip(reco_eiso, -20.0, 30.0),
            np.nan,
        )
        cols.add("reco_eiso_clip30")
    if "reco_eiso_over_cluster_Et" not in cols and {"reco_eiso", "cluster_Et"}.issubset(cols):
        reco_eiso = frame["reco_eiso"].to_numpy(dtype="float64")
        cluster_et = frame["cluster_Et"].to_numpy(dtype="float64")
        out = np.full(len(reco_eiso), np.nan, dtype="float64")
        mask = (
            np.isfinite(reco_eiso)
            & np.isfinite(cluster_et)
            & (np.abs(reco_eiso) < 1.0e8)
            & (cluster_et > 1.0e-6)
        )
        out[mask] = np.clip(reco_eiso[mask] / cluster_et[mask], -2.0, 3.0)
        frame["reco_eiso_over_cluster_Et"] = out
        cols.add("reco_eiso_over_cluster_Et")
    if "reco_eiso_signed_log1p" not in cols and "reco_eiso" in cols:
        reco_eiso = frame["reco_eiso"].to_numpy(dtype="float64")
        clipped = np.where(
            np.isfinite(reco_eiso) & (np.abs(reco_eiso) < 1.0e8),
            np.clip(reco_eiso, -20.0, 60.0),
            np.nan,
        )
        frame["reco_eiso_signed_log1p"] = np.sign(clipped) * np.log1p(np.abs(clipped))
    return frame


def stable_seed(*items: object) -> int:
    text = "|".join(str(item) for item in items)
    digest = hashlib.sha256(text.encode("utf-8")).hexdigest()
    return int(digest[:8], 16)


def infer_source_sample(path: Path) -> str:
    text = str(path)
    for sample, aliases in PPG12_EXACT_SAMPLE_ALIASES:
        if any(alias in text for alias in aliases):
            return sample
    return "unknown"


def expand_input_paths(items: list[Path]) -> list[Path]:
    paths: list[Path] = []
    for item in items:
        text = str(item)
        if text.startswith("@"):
            manifest = Path(text[1:])
            if not manifest.is_file():
                raise SystemExit(f"Input manifest does not exist: {manifest}")
            for raw in manifest.read_text().splitlines():
                line = raw.strip()
                if line and not line.startswith("#"):
                    paths.append(Path(line))
        else:
            paths.append(item)
    if not paths:
        raise SystemExit("No input ROOT files supplied")
    if os.environ.get("RJ_AUAU_BDT_SKIP_INPUT_EXISTENCE_CHECK", "0") != "1":
        missing = [str(path) for path in paths if not path.is_file()]
        if missing:
            preview = "\n  ".join(missing[:20])
            extra = "" if len(missing) <= 20 else f"\n  ... {len(missing) - 20} more"
            raise SystemExit(f"Input ROOT files are missing:\n  {preview}{extra}")
    return paths


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_low_calo_cut(path: Path) -> dict:
    if not path.is_file():
        raise SystemExit(f"Low-calo cut JSON does not exist: {path}")
    payload = json.loads(path.read_text())
    envelope = payload.get("envelope")
    if not isinstance(envelope, list) or not envelope:
        raise SystemExit(f"Low-calo cut JSON has no envelope table: {path}")
    rows = []
    for item in envelope:
        try:
            row = {
                "cent_lo": float(item["cent_lo"]),
                "cent_hi": float(item["cent_hi"]),
                "threshold": float(item["threshold"]),
                "median": float(item.get("median", math.nan)),
                "mad_sigma": float(item.get("mad_sigma", math.nan)),
                "quantile_floor": float(item.get("quantile_floor", math.nan)),
                "n_events": int(item.get("n_events", 0)),
                "status": str(item.get("status", "unknown")),
            }
        except KeyError as exc:
            raise SystemExit(f"Low-calo cut envelope row missing required field {exc}: {path}") from exc
        rows.append(row)
    rows.sort(key=lambda row: (row["cent_lo"], row["cent_hi"]))
    return {
        "path": str(path),
        "sha256": file_sha256(path),
        "schema": payload.get("schema"),
        "source_blind": bool(payload.get("source_blind", False)),
        "truth_blind": bool(payload.get("truth_blind", False)),
        "bdt_score_blind": bool(payload.get("bdt_score_blind", False)),
        "centrality_source": payload.get("centrality_source"),
        "cut_variable": payload.get("cut_variable"),
        "mad_scale": payload.get("mad_scale"),
        "quantile_floor": payload.get("quantile_floor"),
        "slice_width_percent": payload.get("slice_width_percent"),
        "envelope": rows,
    }


def _series_to_str(frame, name: str):
    return frame[name].astype(str).to_numpy(dtype=object, copy=False)


def ensure_global_event_identity(frame):
    missing = [col for col in GLOBAL_EVENT_KEY_COLUMNS if col not in frame.columns]
    if missing:
        return frame
    return frame


def _factorized_source_codes(frame):
    import numpy as np
    import pandas as pd

    if "source_sample" not in frame.columns:
        return np.zeros(len(frame), dtype="int32"), ["unknown"]
    codes, labels = pd.factorize(frame["source_sample"].astype(str), sort=True)
    labels = [str(item) for item in labels.tolist()]
    return codes.astype("int32", copy=False), labels


def _numeric_event_records(frame, columns: list[str]):
    import numpy as np

    missing = [col for col in columns if col not in frame.columns]
    if missing:
        return None
    dtype = []
    payload = {}
    for col in columns:
        if col == "source_sample":
            values, _ = _factorized_source_codes(frame)
            dtype.append((col, "i4"))
            payload[col] = values
        else:
            values = frame[col].to_numpy()
            dtype.append((col, "i8"))
            payload[col] = values.astype("int64", copy=False)
    records = np.empty(len(frame), dtype=dtype)
    for col, values in payload.items():
        records[col] = values
    return records


def _numeric_event_identity(frame, *, require_file_qualified: bool = False):
    import numpy as np

    columns = list(GLOBAL_EVENT_KEY_COLUMNS)
    records = _numeric_event_records(frame, columns)
    if records is None and not require_file_qualified:
        columns = ["source_sample", "run", "evt"]
        records = _numeric_event_records(frame, columns)
    if records is not None:
        _, inverse, counts = np.unique(records, return_inverse=True, return_counts=True)
        return {
            "status": "ok",
            "columns": columns,
            "inverse": inverse.astype("int64", copy=False),
            "counts": counts.astype("int64", copy=False),
            "unique_events": int(len(counts)),
            "uses_numeric_identity": True,
        }
    if "global_event_key" in frame.columns:
        keys = frame["global_event_key"].astype(str).to_numpy(dtype=object, copy=False)
        _, inverse, counts = np.unique(keys, return_inverse=True, return_counts=True)
        return {
            "status": "ok",
            "columns": ["global_event_key"],
            "inverse": inverse.astype("int64", copy=False),
            "counts": counts.astype("int64", copy=False),
            "unique_events": int(len(counts)),
            "uses_numeric_identity": False,
        }
    missing = [col for col in GLOBAL_EVENT_KEY_COLUMNS if col not in frame.columns]
    if require_file_qualified:
        raise SystemExit("File-qualified global event key is required but missing: " + ", ".join(missing))
    fallback = ["source_sample", "run", "evt"]
    missing_fallback = [col for col in fallback if col not in frame.columns]
    if missing_fallback:
        raise SystemExit("Cannot build event key; missing: " + ", ".join(missing_fallback))
    keys, columns = global_event_keys(frame, require_file_qualified=False)
    _, inverse, counts = np.unique(keys, return_inverse=True, return_counts=True)
    return {
        "status": "ok",
        "columns": columns,
        "inverse": inverse.astype("int64", copy=False),
        "counts": counts.astype("int64", copy=False),
        "unique_events": int(len(counts)),
        "uses_numeric_identity": False,
    }


def _event_count_for_mask(event_inverse, n_events: int, mask) -> int:
    import numpy as np

    mask = np.asarray(mask, dtype=bool)
    if not mask.any():
        return 0
    counts = np.bincount(event_inverse[mask], minlength=n_events)
    return int(np.count_nonzero(counts))


def _event_key_audit_from_counts(frame, columns: list[str], counts):
    import numpy as np

    counts = np.asarray(counts, dtype="int64")
    counts = counts[counts > 0]
    event_key_columns_available = [col for col in GLOBAL_EVENT_KEY_COLUMNS if col in frame.columns]
    missing_event_key_columns = [col for col in GLOBAL_EVENT_KEY_COLUMNS if col not in frame.columns]
    report = {
        "present": all(col in frame.columns for col in GLOBAL_EVENT_KEY_COLUMNS),
        "available_event_key_columns": event_key_columns_available,
        "missing_event_key_columns": missing_event_key_columns,
        "event_key_columns": list(columns),
        "candidate_key_columns": [col for col in GLOBAL_CANDIDATE_KEY_COLUMNS if col in frame.columns],
        "candidate_rows": int(len(frame)),
        "unique_events": int(len(counts)),
        "candidate_multiplicity_mean": float(np.mean(counts)) if len(counts) else math.nan,
        "candidate_multiplicity_max": int(np.max(counts)) if len(counts) else 0,
        "event_dedup_relies_only_on_source_run_evt": list(columns) != GLOBAL_EVENT_KEY_COLUMNS,
        "source_run_evt_unique_keys": None,
        "source_run_evt_repeated_after_global_event_dedup": None,
        "audit_status": "ok",
    }
    if all(col in frame.columns for col in ["source_sample", "run", "evt"]):
        simple = _numeric_event_records(frame, ["source_sample", "run", "evt"])
        simple_unique = np.unique(simple)
        report["source_run_evt_unique_keys"] = int(len(simple_unique))
        report["source_run_evt_repeated_after_global_event_dedup"] = int(len(counts) - len(simple_unique))
    return report


def _centrality_thresholds(cent, envelope: list[dict]):
    import numpy as np

    thresholds = np.full(len(cent), math.nan, dtype="float64")
    for row in envelope:
        lo = float(row["cent_lo"])
        hi = float(row["cent_hi"])
        threshold = float(row["threshold"])
        thresholds[(cent >= lo) & (cent < hi)] = threshold
    return thresholds


def global_event_keys(frame, *, require_file_qualified: bool = False):
    import numpy as np

    if "global_event_key" in frame.columns:
        return frame["global_event_key"].astype(str).to_numpy(dtype=object, copy=False), list(GLOBAL_EVENT_KEY_COLUMNS)
    missing = [col for col in GLOBAL_EVENT_KEY_COLUMNS if col not in frame.columns]
    if not missing:
        source = _series_to_str(frame, "source_sample")
        file_index = frame["input_file_index"].to_numpy()
        run = frame["run"].to_numpy()
        evt = frame["evt"].to_numpy()
        return np.asarray(
            [f"{source[i]}|{int(file_index[i])}|{int(run[i])}|{int(evt[i])}" for i in range(len(frame))],
            dtype=object,
        ), list(GLOBAL_EVENT_KEY_COLUMNS)
    if require_file_qualified:
        raise SystemExit("File-qualified global event key is required but missing: " + ", ".join(missing))
    fallback = ["source_sample", "run", "evt"]
    missing_fallback = [col for col in fallback if col not in frame.columns]
    if missing_fallback:
        raise SystemExit("Cannot build event key; missing: " + ", ".join(missing_fallback))
    source = _series_to_str(frame, "source_sample")
    run = frame["run"].to_numpy()
    evt = frame["evt"].to_numpy()
    return np.asarray(
        [f"{source[i]}|{int(run[i])}|{int(evt[i])}" for i in range(len(frame))],
        dtype=object,
    ), fallback


def summarize_global_event_key(frame) -> dict:
    import numpy as np

    event_key_columns_available = [col for col in GLOBAL_EVENT_KEY_COLUMNS if col in frame.columns]
    missing_event_key_columns = [col for col in GLOBAL_EVENT_KEY_COLUMNS if col not in frame.columns]
    report = {
        "present": all(col in frame.columns for col in GLOBAL_EVENT_KEY_COLUMNS),
        "available_event_key_columns": event_key_columns_available,
        "missing_event_key_columns": missing_event_key_columns,
        "event_key_columns": [],
        "candidate_key_columns": [col for col in GLOBAL_CANDIDATE_KEY_COLUMNS if col in frame.columns],
        "candidate_rows": int(len(frame)),
        "unique_events": 0,
        "candidate_multiplicity_mean": math.nan,
        "candidate_multiplicity_max": 0,
        "event_dedup_relies_only_on_source_run_evt": True,
        "source_run_evt_unique_keys": None,
        "source_run_evt_repeated_after_global_event_dedup": None,
    }
    if len(frame) == 0:
        return report
    if not report["present"] and not all(col in frame.columns for col in ["source_sample", "run", "evt"]):
        report["audit_status"] = "missing_event_key_columns"
        return report
    identity = _numeric_event_identity(frame, require_file_qualified=False)
    report["audit_status"] = "ok"
    report["event_key_columns"] = identity["columns"]
    report["event_dedup_relies_only_on_source_run_evt"] = identity["columns"] != GLOBAL_EVENT_KEY_COLUMNS
    counts = identity["counts"]
    report["unique_events"] = int(identity["unique_events"])
    if len(counts):
        report["candidate_multiplicity_mean"] = float(np.mean(counts))
        report["candidate_multiplicity_max"] = int(np.max(counts))
    if all(col in frame.columns for col in ["source_sample", "run", "evt"]):
        simple = _numeric_event_records(frame, ["source_sample", "run", "evt"])
        simple_unique = np.unique(simple)
        report["source_run_evt_unique_keys"] = int(len(simple_unique))
        report["source_run_evt_repeated_after_global_event_dedup"] = int(identity["unique_events"] - len(simple_unique))
    return report


def centrality_slice_threshold(cent: float, envelope: list[dict]) -> float:
    for row in envelope:
        lo = float(row["cent_lo"])
        hi = float(row["cent_hi"])
        if cent >= lo and cent < hi:
            return float(row["threshold"])
    return math.nan


def _fraction_table(frame, reject_mask, by: str, event_inverse, n_events: int):
    import numpy as np

    rows = []
    if by == "centrality_bin":
        cent = frame["centrality"].to_numpy(dtype="float64")
        categories = [("0-20", (cent >= 0.0) & (cent < 20.0)), ("20-50", (cent >= 20.0) & (cent < 50.0)), ("50-80", (cent >= 50.0) & (cent < 80.0))]
    elif by == "source_sample":
        source = frame["source_sample"].astype(str).to_numpy() if "source_sample" in frame.columns else np.full(len(frame), "unknown", dtype=object)
        categories = [(name, source == name) for name in sorted(np.unique(source).tolist())]
    else:
        raise ValueError(by)
    for name, mask in categories:
        mask = np.asarray(mask, dtype=bool)
        if not mask.any():
            continue
        total_candidates = int(mask.sum())
        rejected_candidates = int(np.sum(mask & reject_mask))
        total_events = _event_count_for_mask(event_inverse, n_events, mask)
        rejected_events = _event_count_for_mask(event_inverse, n_events, mask & reject_mask)
        rows.append(
            {
                by: name,
                "total_events": total_events,
                "rejected_events": rejected_events,
                "event_rejection_fraction": rejected_events / total_events if total_events else math.nan,
                "total_candidates": total_candidates,
                "rejected_candidates": rejected_candidates,
                "candidate_rejection_fraction": rejected_candidates / total_candidates if total_candidates else math.nan,
                "candidate_minus_event_fraction": (
                    (rejected_candidates / total_candidates) - (rejected_events / total_events)
                    if total_candidates and total_events
                    else math.nan
                ),
            }
        )
    return rows


def _component_median_table(frame, reject_mask, event_inverse, n_events: int):
    import numpy as np

    components = [
        ("CEMC", "event_calo_cemc_energy"),
        ("IHCal", "event_calo_ihcal_energy"),
        ("OHCal", "event_calo_ohcal_energy"),
        ("total_calo", "event_calo_total_energy"),
    ]
    available = [(label, col) for label, col in components if col in frame.columns]
    if not available:
        return []
    cent = frame["centrality"].to_numpy(dtype="float64")
    component_arrays = {col: frame[col].to_numpy(dtype="float64") for _, col in available}
    event_reject = np.bincount(event_inverse[reject_mask], minlength=n_events) > 0
    rows = []
    for cent_label, cent_mask in [("0-20", (cent >= 0.0) & (cent < 20.0)), ("20-50", (cent >= 20.0) & (cent < 50.0)), ("50-80", (cent >= 50.0) & (cent < 80.0))]:
        if not cent_mask.any():
            continue
        idx = np.flatnonzero(cent_mask)
        first_row = np.full(n_events, len(frame), dtype="int64")
        np.minimum.at(first_row, event_inverse[idx], idx)
        valid = first_row < len(frame)
        first_idx = first_row[valid]
        reject_events = event_reject[valid]
        for label, col in available:
            values = component_arrays[col][first_idx]
            retained = values[~reject_events]
            rejected = values[reject_events]
            retained = retained[np.isfinite(retained)]
            rejected = rejected[np.isfinite(rejected)]
            retained_median = float(np.median(retained)) if len(retained) else math.nan
            rejected_median = float(np.median(rejected)) if len(rejected) else math.nan
            rows.append(
                {
                    "centrality_bin": cent_label,
                    "component": label,
                    "retained_median": retained_median,
                    "rejected_median": rejected_median,
                    "rejected_over_retained": (
                        rejected_median / retained_median
                        if math.isfinite(retained_median) and retained_median != 0.0
                        else math.nan
                    ),
                }
            )
    return rows


def apply_low_calo_event_quality_filter(frame, cut_json: Path | None, audit_output: Path | None = None, audit_only: bool = False):
    import numpy as np

    if cut_json is None:
        return frame, {"enabled": False}
    missing = [col for col in EVENT_QUALITY_DIRECT_COLUMNS if col not in frame.columns]
    if missing:
        raise SystemExit("Low-calo event-quality filter missing required column(s): " + ", ".join(missing))
    cut = load_low_calo_cut(cut_json)
    frame = ensure_global_event_identity(frame.copy())
    event_identity = _numeric_event_identity(frame, require_file_qualified=True)
    event_inverse = event_identity["inverse"]
    event_key_columns = event_identity["columns"]
    n_events = int(event_identity["unique_events"])
    cent = frame["centrality"].to_numpy(dtype="float64")
    log_calo = frame["event_calo_log10_total_energy_plus1"].to_numpy(dtype="float64")
    thresholds = _centrality_thresholds(cent, cut["envelope"])
    in_range = np.isfinite(cent) & np.isfinite(log_calo) & np.isfinite(thresholds)
    row_below = in_range & (log_calo < thresholds)
    event_reject = np.bincount(event_inverse[row_below], minlength=n_events) > 0
    reject = event_reject[event_inverse]
    keep = ~reject
    retained_below = keep & in_range & (log_calo < thresholds)
    kept_events = np.bincount(event_inverse[keep], minlength=n_events) > 0
    retained_below_events = np.bincount(event_inverse[retained_below], minlength=n_events) > 0
    report = {
        "schema": "AUAU_LOW_CALO_UPSTREAM_FILTER_AUDIT_V1",
        "enabled": True,
        "audit_only": bool(audit_only),
        "cut_json_path": str(cut_json),
        "cut_json_sha256": cut["sha256"],
        "cut_json_schema": cut.get("schema"),
        "source_blind": cut.get("source_blind"),
        "truth_blind": cut.get("truth_blind"),
        "bdt_score_blind": cut.get("bdt_score_blind"),
        "variables_used_for_cut": list(EVENT_QUALITY_DIRECT_COLUMNS),
        "variables_not_used_for_cut": [
            "source_sample",
            "is_signal",
            "truth photon label",
            "truth isolation label",
            "BDT score",
            "cluster_Et",
            "cluster_Eta",
            "candidate shower-shape variables",
            "train/test split",
            "sample weight",
            "WP80 behavior",
        ],
        "boundary_convention": "centrality slice uses cent_lo <= centrality < cent_hi; candidate/event is rejected only when log10(total calo + 1) < threshold; equality is retained",
        "threshold_table": cut["envelope"],
        "rows_before": int(len(frame)),
        "rows_after": int(keep.sum()),
        "rows_rejected": int(reject.sum()),
        "globally_unique_events_before": int(n_events),
        "globally_unique_events_after": int(np.count_nonzero(kept_events)),
        "globally_unique_events_rejected": int(np.count_nonzero(event_reject)),
        "retained_below_envelope_events": int(np.count_nonzero(retained_below_events)),
        "retained_below_envelope_candidates": int(retained_below.sum()),
        "out_of_envelope_range_candidates_retained": int((~in_range).sum()),
        "event_key_audit_before": _event_key_audit_from_counts(frame, event_key_columns, event_identity["counts"]),
        "event_key_columns_used": event_key_columns,
        "rejection_by_centrality": _fraction_table(frame, reject, "centrality_bin", event_inverse, n_events),
        "rejection_by_source": _fraction_table(frame, reject, "source_sample", event_inverse, n_events),
        "component_medians_by_centrality": _component_median_table(frame, reject, event_inverse, n_events),
    }
    retained = frame.loc[keep].copy()
    kept_event_counts = np.bincount(event_inverse[keep], minlength=n_events)
    report["event_key_audit_after"] = _event_key_audit_from_counts(retained, event_key_columns, kept_event_counts)
    if audit_output is not None:
        audit_output.parent.mkdir(parents=True, exist_ok=True)
        audit_output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        print(f"[OK] wrote low-calo event-quality filter audit: {audit_output}", flush=True)
    return retained, report


def assume_low_calo_event_quality_filter_applied(
    frame,
    cut_json: Path,
    audit_output: Path | None,
    outdir: Path,
):
    if audit_output is not None:
        candidates = [audit_output]
    else:
        candidates = [
            outdir / "event_quality_filter_training_audit.json",
            outdir / "event_quality_filter_audit.json",
        ]
    existing = next((path for path in candidates if path.is_file()), None)
    if existing is None:
        searched = ", ".join(str(path) for path in candidates)
        raise SystemExit(
            "Cannot assume low-calo event-quality filter is already applied: "
            f"no audit JSON found; searched {searched}"
        )
    cut = load_low_calo_cut(cut_json)
    report = json.loads(existing.read_text())
    if report.get("cut_json_sha256") != cut["sha256"]:
        raise SystemExit(
            "Cannot assume low-calo event-quality filter is already applied: "
            f"audit cut hash {report.get('cut_json_sha256')} does not match {cut['sha256']}"
        )
    if (
        report.get("retained_below_envelope_events", 0) != 0
        or report.get("retained_below_envelope_candidates", 0) != 0
    ):
        raise SystemExit(
            "Cannot assume low-calo event-quality filter is already applied: "
            "audit has retained below-envelope events/candidates"
        )
    rows_after = report.get("rows_after", report.get("candidate_rows_after"))
    if rows_after is None:
        raise SystemExit(
            "Cannot assume low-calo event-quality filter is already applied: "
            "audit JSON does not contain rows_after"
        )
    if int(rows_after) != len(frame):
        raise SystemExit(
            "Cannot assume low-calo event-quality filter is already applied: "
            f"audit rows_after={rows_after} but cache rows={len(frame)}"
        )
    assumed = dict(report)
    assumed["assumed_already_applied"] = True
    assumed["audit_output_used"] = str(existing)
    assumed["validated_filtered_rows"] = int(len(frame))
    print(f"[OK] using pre-applied low-calo event-quality filter audit: {existing}", flush=True)
    return frame, assumed


def finite_mask(frame, columns: Iterable[str]):
    import numpy as np

    mask = np.ones(len(frame), dtype=bool)
    for col in columns:
        mask &= np.isfinite(frame[col].to_numpy())
    return mask


def load_frame(
    paths: list[Path],
    tree_name: str,
    required_columns: list[str],
    optional_columns: list[str],
    missing_label_branch: str | None,
    missing_label_value: int | None,
    skip_missing_tree: bool = False,
    label_branch: str | None = None,
    max_load_rows_per_class: int = 0,
    max_load_rows: int = 0,
    load_sample_seed: int = 42,
):
    try:
        import pandas as pd
        import uproot
    except ImportError as exc:
        raise SystemExit(
            "Missing dependency. Use an environment with uproot, pandas, numpy, "
            "scikit-learn, xgboost, and PyROOT for TMVA export."
        ) from exc

    import numpy as np

    frames = []
    seen_optional: set[str] = set()
    skipped_missing_tree: list[str] = []
    load_class_counts = {0: 0, 1: 0}
    total_loaded_rows = 0
    load_cap_enabled = max_load_rows_per_class > 0 or max_load_rows > 0
    file_index_by_path = {str(path): idx for idx, path in enumerate(paths)}
    iter_paths = list(paths)
    rng = np.random.default_rng(load_sample_seed)
    if load_cap_enabled:
        rng.shuffle(iter_paths)
        print(
            "[INFO] load-time row cap enabled: "
            f"max_load_rows_per_class={max_load_rows_per_class} "
            f"max_load_rows={max_load_rows} seed={load_sample_seed}",
            flush=True,
        )
    for path in iter_paths:
        with uproot.open(path) as root_file:
            try:
                tree = root_file[tree_name]
            except Exception as exc:
                if skip_missing_tree:
                    skipped_missing_tree.append(str(path))
                    print(f"[WARN] skipping {path}: missing tree {tree_name}", flush=True)
                    continue
                raise SystemExit(f"{path} does not contain tree {tree_name}") from exc
            keys = set(tree.keys())
            missing = [col for col in required_columns if col not in keys]
            allow_missing_label = (
                missing_label_branch is not None
                and missing == [missing_label_branch]
                and missing_label_value is not None
            )
            if allow_missing_label:
                missing = []
            if missing:
                raise SystemExit(f"{path}:{tree_name} missing required branches: {', '.join(missing)}")
            present_optional = [col for col in optional_columns if col in keys]
            present_identity = [
                col
                for col in TRAINING_IDENTITY_OPTIONAL_COLUMNS
                if col in keys and col not in required_columns and col not in present_optional
            ]
            seen_optional.update(present_optional)
            read_columns = [col for col in required_columns if col in keys] + present_optional + present_identity
            frame = tree.arrays(read_columns, library="pd")
            if allow_missing_label:
                frame[missing_label_branch] = int(missing_label_value)
            if "source_sample" not in frame.columns:
                frame["source_sample"] = infer_source_sample(path)
            frame["input_file"] = str(path)
            frame["input_file_index"] = int(file_index_by_path.get(str(path), -1))
            frame["input_tree_entry"] = np.arange(len(frame), dtype="int64")
            if load_cap_enabled and label_branch and label_branch in frame.columns:
                keep_parts = []
                for cls in (0, 1):
                    idx = np.flatnonzero(frame[label_branch].to_numpy(dtype="int32", copy=False) == cls)
                    if len(idx) == 0:
                        continue
                    remaining_class = len(idx)
                    if max_load_rows_per_class > 0:
                        remaining_class = min(remaining_class, max_load_rows_per_class - load_class_counts[cls])
                    remaining_total = len(idx)
                    if max_load_rows > 0:
                        remaining_total = min(remaining_total, max_load_rows - total_loaded_rows)
                    n_take = min(len(idx), remaining_class, remaining_total)
                    if n_take <= 0:
                        continue
                    if n_take < len(idx):
                        idx = rng.choice(idx, size=n_take, replace=False)
                    keep_parts.append(frame.iloc[np.sort(idx)])
                    load_class_counts[cls] += int(n_take)
                    total_loaded_rows += int(n_take)
                if keep_parts:
                    frame = pd.concat(keep_parts, ignore_index=True, copy=False)
                else:
                    continue
            frames.append(frame)
        if load_cap_enabled:
            class_cap_done = (
                max_load_rows_per_class > 0
                and all(load_class_counts[cls] >= max_load_rows_per_class for cls in (0, 1))
            )
            total_cap_done = max_load_rows > 0 and total_loaded_rows >= max_load_rows
            if class_cap_done or total_cap_done:
                print(
                    "[INFO] stopping input read after load-time caps: "
                    f"counts={load_class_counts} total={total_loaded_rows}",
                    flush=True,
                )
                break

    if not frames:
        raise SystemExit(f"No input files contained tree {tree_name}")
    if skipped_missing_tree:
        print(
            f"[WARN] skipped {len(skipped_missing_tree)} input files with no {tree_name}; "
            f"first={skipped_missing_tree[0]}",
            flush=True,
        )
    frame = add_derived_features(pd.concat(frames, ignore_index=True))
    ensure_global_event_identity(frame)
    if load_cap_enabled:
        print(
            "[INFO] loaded capped frame: "
            f"rows={len(frame)} class_counts={load_class_counts} total_selected={total_loaded_rows}",
            flush=True,
        )
    for col in optional_columns:
        if col not in frame.columns:
            frame[col] = 1.0
    return frame, sorted(seen_optional)


def save_frame_cache(frame, path: Path, columns: list[str]) -> None:
    import numpy as np

    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {col: frame[col].to_numpy() for col in columns if col in frame.columns}
    payload["__columns__"] = np.asarray(sorted(payload), dtype=object)
    np.savez_compressed(path, **payload)


def load_frame_cache(path: Path):
    import pandas as pd
    import numpy as np

    if not path.is_file():
        raise SystemExit(f"Training cache does not exist: {path}")
    data = np.load(path, allow_pickle=True)
    columns = [str(x) for x in data["__columns__"].tolist()]
    return pd.DataFrame({col: data[col] for col in columns})


def load_or_build_frame(
    paths: list[Path],
    tree_name: str,
    required_columns: list[str],
    optional_columns: list[str],
    label_branch: str,
    missing_label_value: int | None,
    cache_file: Path | None,
    cache_only: bool = False,
    skip_missing_tree: bool = False,
    max_load_rows_per_class: int = 0,
    max_load_rows: int = 0,
    load_sample_seed: int = 42,
):
    if cache_file is not None and cache_file.is_file():
        frame = add_derived_features(load_frame_cache(cache_file))
        ensure_global_event_identity(frame)
        missing = [col for col in required_columns if col not in frame.columns]
        if missing:
            raise SystemExit(f"Training cache {cache_file} is missing required columns: {', '.join(missing)}")
        for col in optional_columns:
            if col not in frame.columns:
                frame[col] = 1.0
        return frame, ["cache"], True

    frame, optional_seen = load_frame(
        paths,
        tree_name,
        expand_required_columns(required_columns),
        optional_columns,
        label_branch,
        missing_label_value,
        skip_missing_tree=skip_missing_tree,
        label_branch=label_branch,
        max_load_rows_per_class=max_load_rows_per_class,
        max_load_rows=max_load_rows,
        load_sample_seed=load_sample_seed,
    )
    if cache_file is not None:
        cache_cols = sorted(
            set(
                expand_required_columns(required_columns)
                + required_columns
                + optional_columns
                + [
                    "source_sample",
                    "input_file_index",
                    "input_tree_entry",
                ]
                + TRAINING_IDENTITY_OPTIONAL_COLUMNS
            )
        )
        save_frame_cache(frame, cache_file, cache_cols)
        print(f"[OK] wrote training cache: {cache_file}")
    if cache_only:
        return frame, optional_seen, False
    return frame, optional_seen, False


def normalize_mean_one(values):
    import numpy as np

    values = np.asarray(values, dtype="float64")
    finite = np.isfinite(values) & (values > 0.0)
    if not finite.any():
        return np.ones_like(values, dtype="float64")
    mean = float(values[finite].mean())
    if mean <= 0.0 or not math.isfinite(mean):
        return np.ones_like(values, dtype="float64")
    values = np.where(finite, values / mean, 1.0)
    return values


def inverse_pdf_factors(values, labels, nbins: int, max_factor: float):
    import numpy as np

    values = np.asarray(values, dtype="float64")
    labels = np.asarray(labels, dtype="int32")
    factors = np.ones(len(values), dtype="float64")
    finite = np.isfinite(values)
    if finite.sum() < max(10, nbins):
        return factors

    lo, hi = np.nanpercentile(values[finite], [1.0, 99.0])
    if not math.isfinite(float(lo)) or not math.isfinite(float(hi)) or hi <= lo:
        return factors
    edges = np.linspace(lo, hi, nbins + 1)

    for cls in (0, 1):
        cls_mask = finite & (labels == cls) & (values >= lo) & (values <= hi)
        if cls_mask.sum() < nbins:
            continue
        counts, _ = np.histogram(values[cls_mask], bins=edges)
        counts = counts.astype("float64")
        nonzero = counts > 0.0
        if not nonzero.any():
            continue
        target = float(np.mean(counts[nonzero]))
        idx = np.clip(np.searchsorted(edges, values[cls_mask], side="right") - 1, 0, nbins - 1)
        local = np.ones(cls_mask.sum(), dtype="float64")
        valid = counts[idx] > 0.0
        local[valid] = target / counts[idx[valid]]
        local = np.clip(local, 1.0 / max_factor, max_factor)
        factors[cls_mask] *= normalize_mean_one(local)

    return factors


def ppg12_exact_inverse_pdf_weights(values, *, n_bins: int, fixed_range=None, weight_cap: float | None = None):
    import numpy as np
    from scipy.interpolate import UnivariateSpline

    values = np.asarray(values, dtype="float64")
    weights = np.ones(len(values), dtype="float64")
    finite = np.isfinite(values)
    report = {
        "n_entries": int(len(values)),
        "n_finite": int(finite.sum()),
        "n_bins": int(n_bins),
        "fixed_range": list(fixed_range) if fixed_range is not None else None,
        "weight_cap": float(weight_cap) if weight_cap is not None else None,
        "status": "ok",
    }
    if finite.sum() < max(10, n_bins):
        report["status"] = "insufficient_finite_values"
        return weights, report

    vals = values[finite]
    if fixed_range is None:
        lo, hi = float(np.min(vals)), float(np.max(vals))
    else:
        lo, hi = float(fixed_range[0]), float(fixed_range[1])
    report["range"] = [lo, hi]
    if not math.isfinite(lo) or not math.isfinite(hi) or hi <= lo:
        report["status"] = "invalid_range"
        return weights, report

    edges = np.linspace(lo, hi, n_bins + 1)
    try:
        hist, bin_edges = np.histogram(vals, bins=edges, density=True)
    except Exception as exc:  # noqa: BLE001
        report["status"] = f"histogram_failed: {exc}"
        return weights, report
    hist = hist.astype("float64") * float(n_bins)
    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    good_hist = np.isfinite(hist)
    if good_hist.sum() < 4:
        report["status"] = "insufficient_histogram_support"
        return weights, report
    try:
        spline = UnivariateSpline(centers[good_hist], hist[good_hist], s=0.0)
        pdf = spline(vals)
    except Exception as exc:  # noqa: BLE001
        report["status"] = f"spline_failed: {exc}"
        return weights, report

    report["pdf_min_before_clip"] = float(np.nanmin(pdf)) if len(pdf) else math.nan
    report["pdf_max_before_clip"] = float(np.nanmax(pdf)) if len(pdf) else math.nan
    report["pdf_negative_fraction_before_clip"] = float(np.mean(pdf < 0.0)) if len(pdf) else math.nan
    pdf = np.clip(pdf, a_min=1.0e-3, a_max=None)
    local = 1.0 / pdf
    if weight_cap is not None:
        local = np.clip(local, a_min=None, a_max=float(weight_cap))
    local = normalize_mean_one(local)
    weights[finite] = local
    report["min_weight"] = float(np.min(local)) if len(local) else math.nan
    report["max_weight"] = float(np.max(local)) if len(local) else math.nan
    report["mean_weight"] = float(np.mean(local)) if len(local) else math.nan
    return weights, report


def compute_ppg12_exact_global_weights(frame, label_branch: str):
    import numpy as np

    labels = frame[label_branch].to_numpy(dtype="int32")
    weights = np.ones(len(frame), dtype="float64")
    report: dict[str, object] = {
        "weight_mode": "ppg12-exact",
        "event_weight_used": False,
        "vertex_reweight": False,
        "centrality_event_weight": False,
        "cross_section_weight_used_for_training": False,
        "weights_computed_before_binning": True,
        "eta_range": list(PPG12_EXACT_ETA_RANGE),
        "eta_bins": PPG12_EXACT_N_BINS,
        "et_bins": PPG12_EXACT_N_BINS,
        "et_weight_cap": PPG12_EXACT_ET_WEIGHT_CAP,
    }

    class_counts: dict[str, int] = {}
    class_weight_factors: dict[str, float] = {}
    n_total = 0
    for cls in (0, 1):
        n_cls = int((labels == cls).sum())
        class_counts[str(cls)] = n_cls
        n_total += n_cls
    if class_counts["0"] <= 0 or class_counts["1"] <= 0:
        raise SystemExit(f"PPG12-exact weights need both classes; observed counts={class_counts}")
    for cls in (0, 1):
        factor = float(n_total) / (2.0 * float(class_counts[str(cls)]))
        class_weight_factors[str(cls)] = factor
        weights[labels == cls] *= factor

    eta_reports = {}
    et_reports = {}
    for cls in (0, 1):
        mask = labels == cls
        eta_w, eta_report = ppg12_exact_inverse_pdf_weights(
            frame.loc[mask, "cluster_Eta"].to_numpy(dtype="float64"),
            n_bins=PPG12_EXACT_N_BINS,
            fixed_range=PPG12_EXACT_ETA_RANGE,
            weight_cap=None,
        )
        weights[mask] *= eta_w
        eta_reports[str(cls)] = eta_report

        et_w, et_report = ppg12_exact_inverse_pdf_weights(
            frame.loc[mask, "cluster_Et"].to_numpy(dtype="float64"),
            n_bins=PPG12_EXACT_N_BINS,
            fixed_range=None,
            weight_cap=PPG12_EXACT_ET_WEIGHT_CAP,
        )
        weights[mask] *= et_w
        et_reports[str(cls)] = et_report

    finite_positive = np.isfinite(weights) & (weights > 0.0)
    if not finite_positive.all():
        bad = int((~finite_positive).sum())
        raise SystemExit(f"PPG12-exact weights produced {bad} non-finite/non-positive rows")

    report["class_counts"] = class_counts
    report["class_weight_factors"] = class_weight_factors
    report["eta_reweight"] = eta_reports
    report["et_reweight"] = et_reports
    report["sum_weight_class0"] = float(weights[labels == 0].sum())
    report["sum_weight_class1"] = float(weights[labels == 1].sum())
    report["min_weight"] = float(np.min(weights)) if len(weights) else math.nan
    report["max_weight"] = float(np.max(weights)) if len(weights) else math.nan
    report["mean_weight"] = float(np.mean(weights)) if len(weights) else math.nan
    return weights, report


def parse_expected_samples(text: str) -> tuple[str, ...]:
    samples = tuple(item.strip() for item in text.split(",") if item.strip())
    return samples or PPG12_EXACT_EXPECTED_SAMPLES


def validate_ppg12_exact_samples(frame, label_branch: str, expected_samples: tuple[str, ...]) -> dict:
    import numpy as np

    if "source_sample" not in frame.columns:
        raise SystemExit("PPG12-exact mode requires source_sample. Rebuild the training cache from ROOT inputs.")
    samples = sorted(str(x) for x in frame["source_sample"].dropna().unique())
    expected = sorted(expected_samples)
    missing = sorted(set(expected) - set(samples))
    unexpected = sorted(set(samples) - set(expected))
    if missing or unexpected:
        raise SystemExit(
            "PPG12-exact sample set mismatch: "
            f"missing={missing or []} unexpected={unexpected or []} observed={samples}"
        )

    labels = frame[label_branch].to_numpy(dtype="int32")
    sample_arr = frame["source_sample"].astype(str).to_numpy()
    inventory = []
    mixed_label_counts: dict[str, int] = {}
    for sample in expected_samples:
        mask = sample_arr == sample
        n_signal = int(np.sum(mask & (labels == 1)))
        n_background = int(np.sum(mask & (labels == 0)))
        if "Photon" in sample and n_background:
            mixed_label_counts[sample] = n_background
        elif "Jet" in sample and n_signal:
            mixed_label_counts[sample] = n_signal
        inventory.append(
            {
                "source_sample": sample,
                "n_rows": int(mask.sum()),
                "n_signal": n_signal,
                "n_background": n_background,
            }
        )
    if mixed_label_counts:
        print(
            "[WARN] PPG12-exact source_sample/truth-label mixture observed; "
            "treating source_sample as provenance and truth label as the BDT class: "
            f"{mixed_label_counts}",
            flush=True,
        )
    return {
        "expected_samples": list(expected_samples),
        "observed_samples": samples,
        "inventory": inventory,
        "mixed_label_counts": mixed_label_counts,
        "source_sample_semantics": "provenance",
        "truth_label_semantics": "per-candidate BDT class",
        "mixed_labels_are_fatal": False,
    }


def write_ppg12_exact_sample_inventory(frame, label_branch: str, weights, outdir: Path) -> dict:
    import numpy as np
    import pandas as pd

    outdir.mkdir(parents=True, exist_ok=True)
    labels = frame[label_branch].to_numpy(dtype="int32")
    sample_arr = frame["source_sample"].astype(str).to_numpy()
    et = frame["cluster_Et"].to_numpy(dtype="float64")
    eta = frame["cluster_Eta"].to_numpy(dtype="float64")
    cent = frame["centrality"].to_numpy(dtype="float64") if "centrality" in frame.columns else np.full(len(frame), np.nan)
    weights = np.asarray(weights, dtype="float64")

    rows = []
    for sample in sorted(set(sample_arr)):
        for cls in (0, 1):
            mask = (sample_arr == sample) & (labels == cls)
            rows.append(
                {
                    "source_sample": sample,
                    "class": cls,
                    "class_name": "signal" if cls == 1 else "background",
                    "n_rows": int(mask.sum()),
                    "sum_ppg12_exact_weight": float(weights[mask].sum()) if mask.any() else 0.0,
                    "mean_cluster_Et": float(np.nanmean(et[mask])) if mask.any() else math.nan,
                    "mean_cluster_Eta": float(np.nanmean(eta[mask])) if mask.any() else math.nan,
                }
            )
    inventory_csv = outdir / "ppg12_exact_sample_inventory.csv"
    pd.DataFrame(rows).to_csv(inventory_csv, index=False)

    et_edges = np.asarray([15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 30.0, 35.0])
    eta_edges = np.linspace(PPG12_EXACT_ETA_RANGE[0], PPG12_EXACT_ETA_RANGE[1], PPG12_EXACT_N_BINS + 1)
    cent_bins = [(0.0, 20.0), (20.0, 50.0), (50.0, 80.0)]
    binned_rows = []
    for sample in sorted(set(sample_arr)):
        sample_mask = sample_arr == sample
        for cls in (0, 1):
            class_mask = sample_mask & (labels == cls)
            for lo, hi in zip(et_edges[:-1], et_edges[1:]):
                mask = class_mask & (et >= lo) & (et < hi)
                binned_rows.append(
                    {
                        "source_sample": sample,
                        "class": cls,
                        "axis": "cluster_Et",
                        "bin_low": float(lo),
                        "bin_high": float(hi),
                        "n_rows": int(mask.sum()),
                        "sum_ppg12_exact_weight": float(weights[mask].sum()) if mask.any() else 0.0,
                    }
                )
            for lo, hi in zip(eta_edges[:-1], eta_edges[1:]):
                mask = class_mask & (eta >= lo) & (eta < hi)
                binned_rows.append(
                    {
                        "source_sample": sample,
                        "class": cls,
                        "axis": "cluster_Eta",
                        "bin_low": float(lo),
                        "bin_high": float(hi),
                        "n_rows": int(mask.sum()),
                        "sum_ppg12_exact_weight": float(weights[mask].sum()) if mask.any() else 0.0,
                    }
                )
            for lo, hi in cent_bins:
                mask = class_mask & (cent >= lo) & (cent < hi)
                binned_rows.append(
                    {
                        "source_sample": sample,
                        "class": cls,
                        "axis": "centrality",
                        "bin_low": float(lo),
                        "bin_high": float(hi),
                        "n_rows": int(mask.sum()),
                        "sum_ppg12_exact_weight": float(weights[mask].sum()) if mask.any() else 0.0,
                    }
                )
    binned_csv = outdir / "ppg12_exact_sample_inventory_binned.csv"
    pd.DataFrame(binned_rows).to_csv(binned_csv, index=False)
    return {"inventory_csv": str(inventory_csv), "binned_inventory_csv": str(binned_csv)}


def step_density(ax, values, bins, *, weights=None, label: str, color: str, linestyle: str = "-"):
    import numpy as np

    values = np.asarray(values, dtype="float64")
    finite = np.isfinite(values)
    if weights is not None:
        weights = np.asarray(weights, dtype="float64")
        finite &= np.isfinite(weights) & (weights > 0.0)
        weights = weights[finite]
    vals = values[finite]
    if len(vals) == 0:
        return np.zeros(len(bins) - 1, dtype="float64")
    hist, edges = np.histogram(vals, bins=bins, weights=weights)
    width = np.diff(edges)
    norm = float(np.sum(hist * width))
    density = hist / norm if norm > 0.0 else hist
    y = np.r_[density, density[-1] if len(density) else 0.0]
    ax.step(edges, y, where="post", label=label, color=color, linewidth=2.2, linestyle=linestyle)
    return density


def make_ppg12_exact_closure_plots(
    frame,
    label_branch: str,
    weights,
    report: dict,
    outdir: Path,
    expected_samples: tuple[str, ...] | None = None,
) -> dict:
    import numpy as np
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outdir.mkdir(parents=True, exist_ok=True)
    labels = frame[label_branch].to_numpy(dtype="int32")
    et = frame["cluster_Et"].to_numpy(dtype="float64")
    eta = frame["cluster_Eta"].to_numpy(dtype="float64")
    weights = np.asarray(weights, dtype="float64")
    sig = labels == 1
    bkg = labels == 0
    colors = {"signal": "#009E73", "background": "#6B7280"}

    procedure_png = outdir / "ppg12_exact_reweighting_procedure.png"
    fig, ax = plt.subplots(figsize=(14, 7.5))
    ax.axis("off")
    ax.text(0.02, 0.92, "PPG12-exact BDT training reweighting", fontsize=26, fontweight="bold")
    ax.text(
        0.02,
        0.80,
        "Global training weight is computed once before routed ET/centrality slicing.",
        fontsize=17,
        color="#374151",
    )
    lines = [
        "1. Load Photon12+20 signal and Jet12+20+30 background candidates.",
        "2. Class-balance signal and background total weight.",
        "3. For each class independently: flatten cluster eta with a spline inverse-PDF in -0.7 < eta < 0.7.",
        "4. For each class independently: flatten cluster ET with a spline inverse-PDF over the observed ET range.",
        "5. Cap only the ET inverse-PDF weight at 800, matching PPG12.",
        "6. Do not multiply Au+Au event, vertex, centrality, or cross-section weights into BDT training.",
    ]
    y = 0.66
    for line in lines:
        ax.text(0.06, y, line, fontsize=16, color="#111827")
        y -= 0.085
    ax.text(
        0.06,
        0.08,
        "The closure plots must show that ET/eta sample composition is controlled before model performance is interpreted.",
        fontsize=15,
        color="#4B5563",
    )
    fig.tight_layout()
    fig.savefig(procedure_png, dpi=180)
    plt.close(fig)

    et_finite = et[np.isfinite(et)]
    et_lo = float(np.nanmin(et_finite)) if len(et_finite) else 15.0
    et_hi = float(np.nanmax(et_finite)) if len(et_finite) else 35.0
    if et_hi <= et_lo:
        et_lo, et_hi = 15.0, 35.0
    et_bins = np.linspace(et_lo, et_hi, PPG12_EXACT_N_BINS + 1)
    eta_bins = np.linspace(PPG12_EXACT_ETA_RANGE[0], PPG12_EXACT_ETA_RANGE[1], PPG12_EXACT_N_BINS + 1)

    closure_png = outdir / "ppg12_exact_et_eta_weight_closure.png"
    fig, axs = plt.subplots(2, 3, figsize=(17, 9.5))
    axs = axs.ravel()
    step_density(axs[0], et[sig], et_bins, label="Signal raw", color=colors["signal"])
    step_density(axs[0], et[bkg], et_bins, label="Background raw", color=colors["background"])
    axs[0].set_title("Raw cluster ET")
    axs[0].set_xlabel("cluster ET [GeV]")
    axs[0].set_ylabel("Area-normalized density")

    sig_et_w = step_density(axs[1], et[sig], et_bins, weights=weights[sig], label="Signal weighted", color=colors["signal"])
    bkg_et_w = step_density(axs[1], et[bkg], et_bins, weights=weights[bkg], label="Background weighted", color=colors["background"])
    axs[1].set_title("PPG12-weighted cluster ET")
    axs[1].set_xlabel("cluster ET [GeV]")

    ratio = np.divide(sig_et_w, bkg_et_w, out=np.full_like(sig_et_w, np.nan), where=bkg_et_w > 0.0)
    centers = 0.5 * (et_bins[:-1] + et_bins[1:])
    axs[2].axhline(1.0, color="#111827", linewidth=1.4)
    axs[2].plot(centers, ratio, marker="o", color="#2563EB", linewidth=1.8)
    axs[2].set_title("Weighted signal/background ET ratio")
    axs[2].set_xlabel("cluster ET [GeV]")
    axs[2].set_ylabel("density ratio")
    axs[2].set_ylim(0.0, max(2.0, float(np.nanmax(ratio)) * 1.15 if np.isfinite(ratio).any() else 2.0))

    step_density(axs[3], eta[sig], eta_bins, label="Signal raw", color=colors["signal"])
    step_density(axs[3], eta[bkg], eta_bins, label="Background raw", color=colors["background"])
    axs[3].set_title("Raw cluster eta")
    axs[3].set_xlabel("cluster eta")
    axs[3].set_ylabel("Area-normalized density")

    step_density(axs[4], eta[sig], eta_bins, weights=weights[sig], label="Signal weighted", color=colors["signal"])
    step_density(axs[4], eta[bkg], eta_bins, weights=weights[bkg], label="Background weighted", color=colors["background"])
    axs[4].set_title("PPG12-weighted cluster eta")
    axs[4].set_xlabel("cluster eta")

    w_min = max(float(np.nanmin(weights[weights > 0.0])), 1.0e-6)
    w_max = max(float(np.nanmax(weights)), w_min * 1.01)
    weight_bins = np.geomspace(w_min, w_max, 50)
    axs[5].hist(weights[sig], bins=weight_bins, histtype="step", linewidth=2.2, density=True, label="Signal", color=colors["signal"])
    axs[5].hist(weights[bkg], bins=weight_bins, histtype="step", linewidth=2.2, density=True, label="Background", color=colors["background"])
    axs[5].set_xscale("log")
    axs[5].set_title("Final training-weight distribution")
    axs[5].set_xlabel("PPG12-exact training weight")
    axs[5].set_ylabel("Density")

    for ax in axs:
        ax.grid(True, color="#E5E7EB", linewidth=0.8)
        ax.tick_params(direction="in", top=True, right=True)
        handles, labels_local = ax.get_legend_handles_labels()
        if handles:
            ax.legend(frameon=False, fontsize=10)
    fig.suptitle("PPG12-exact ET/eta training-weight closure", fontsize=22, fontweight="bold", y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(closure_png, dpi=180)
    plt.close(fig)

    sample_png = outdir / "ppg12_exact_sample_mix_closure.png"
    sample_arr = frame["source_sample"].astype(str).to_numpy()
    sample_order = expected_samples or PPG12_EXACT_EXPECTED_SAMPLES
    samples = [s for s in sample_order if s in set(sample_arr)]
    extras = sorted(set(sample_arr) - set(samples))
    samples.extend(extras)
    raw_counts = np.asarray([float(np.sum(sample_arr == s)) for s in samples])
    weighted_counts = np.asarray([float(weights[sample_arr == s].sum()) for s in samples])
    raw_frac = raw_counts / raw_counts.sum() if raw_counts.sum() > 0.0 else raw_counts
    weighted_frac = weighted_counts / weighted_counts.sum() if weighted_counts.sum() > 0.0 else weighted_counts
    x = np.arange(len(samples))
    fig, ax = plt.subplots(figsize=(15, 7.5))
    ax.bar(x - 0.18, raw_frac, width=0.36, label="Raw candidate fraction", color="#9CA3AF")
    ax.bar(x + 0.18, weighted_frac, width=0.36, label="PPG12-weighted fraction", color="#14B8A6")
    ax.set_xticks(x)
    ax.set_xticklabels([s.replace("run28_embedded", "") for s in samples], rotation=0, fontsize=12)
    ax.set_ylabel("Fraction of total candidates / total training weight")
    ax.set_title("Source-sample composition before and after PPG12-exact weights", fontsize=18, fontweight="bold")
    ax.grid(True, axis="y", color="#E5E7EB")
    ax.tick_params(direction="in", top=True, right=True)
    ax.legend(frameon=False)
    for i, (r, w) in enumerate(zip(raw_frac, weighted_frac)):
        ax.text(i - 0.18, r + 0.01, f"{100*r:.1f}%", ha="center", va="bottom", fontsize=10)
        ax.text(i + 0.18, w + 0.01, f"{100*w:.1f}%", ha="center", va="bottom", fontsize=10)
    fig.tight_layout()
    fig.savefig(sample_png, dpi=180)
    plt.close(fig)

    return {
        "procedure_png": str(procedure_png),
        "closure_png": str(closure_png),
        "sample_mix_png": str(sample_png),
    }


def prepare_ppg12_exact_global_weights(frame, label_branch: str, args):
    outdir = Path(args.ppg12_exact_closure_dir) if args.ppg12_exact_closure_dir else (args.outdir / "slideReady" / "ppg12_exact_reweight_bdt")
    expected_samples = parse_expected_samples(args.ppg12_exact_expected_samples)
    artifact_mode = getattr(args, "ppg12_exact_closure_artifacts", "full")
    sample_report = validate_ppg12_exact_samples(frame, label_branch, expected_samples)
    weights, weight_report = compute_ppg12_exact_global_weights(frame, label_branch)
    frame = frame.copy()
    frame[PPG12_EXACT_WEIGHT_COLUMN] = weights

    outdir.mkdir(parents=True, exist_ok=True)
    inventory_paths = {}
    plot_paths = {}
    if artifact_mode == "full":
        inventory_paths = write_ppg12_exact_sample_inventory(frame, label_branch, weights, outdir)
        plot_paths = make_ppg12_exact_closure_plots(frame, label_branch, weights, weight_report, outdir, expected_samples)
    metadata = {
        "schema": "AUAU_BDT_PPG12_EXACT_WEIGHT_CLOSURE_V1",
        "artifact_mode": artifact_mode,
        "sample_validation": sample_report,
        "weighting": weight_report,
        "artifacts": {**inventory_paths, **plot_paths},
    }
    metadata_path = outdir / "ppg12_exact_reweighting_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    if artifact_mode == "full":
        print(f"[OK] PPG12-exact reweighting closure written: {outdir}", flush=True)
    else:
        print(f"[OK] PPG12-exact reweighting metadata written: {outdir} artifacts={artifact_mode}", flush=True)
    return frame, {**metadata, "metadata_json": str(metadata_path)}


def ppg12_exact_precomputed_weight_status(frame) -> dict:
    import numpy as np

    if PPG12_EXACT_WEIGHT_COLUMN not in frame.columns:
        return {"usable": False, "reason": "missing", "has_column": False}
    weights = frame[PPG12_EXACT_WEIGHT_COLUMN].to_numpy(dtype="float64")
    finite_positive = np.isfinite(weights) & (weights > 0.0)
    invalid_rows = int((~finite_positive).sum())
    if invalid_rows:
        return {
            "usable": False,
            "reason": "invalid_rows",
            "has_column": True,
            "invalid_rows": invalid_rows,
            "n_rows": int(len(weights)),
        }
    all_unit = bool(len(weights) and np.allclose(weights, 1.0, rtol=0.0, atol=1.0e-12))
    min_weight = float(np.min(weights)) if len(weights) else math.nan
    max_weight = float(np.max(weights)) if len(weights) else math.nan
    mean_weight = float(np.mean(weights)) if len(weights) else math.nan
    return {
        "usable": not all_unit,
        "reason": "ok" if not all_unit else "all_unit_placeholder",
        "has_column": True,
        "n_rows": int(len(weights)),
        "min_weight": min_weight,
        "max_weight": max_weight,
        "mean_weight": mean_weight,
        "all_unit_placeholder": all_unit,
    }


def require_usable_ppg12_exact_precomputed_weights(frame, context: str) -> dict:
    status = ppg12_exact_precomputed_weight_status(frame)
    if status.get("usable"):
        return status
    reason = status.get("reason", "unknown")
    if reason == "all_unit_placeholder":
        raise SystemExit(
            f"PPG12-exact precomputed weights for {context} are all 1.0; "
            "this is the unweighted placeholder, not a valid exact-weight cache."
        )
    if reason == "invalid_rows":
        raise SystemExit(
            f"PPG12-exact precomputed weights for {context} contain "
            f"{status.get('invalid_rows')} invalid rows"
        )
    raise SystemExit(f"PPG12-exact precomputed weights for {context} are not usable: {reason}")


def summarize_ppg12_exact_precomputed_weights(frame, label_branch: str, args) -> dict:
    import numpy as np

    expected_samples = parse_expected_samples(args.ppg12_exact_expected_samples)
    sample_report = validate_ppg12_exact_samples(frame, label_branch, expected_samples)
    status = require_usable_ppg12_exact_precomputed_weights(frame, "campaign setup")
    weights = frame[PPG12_EXACT_WEIGHT_COLUMN].to_numpy(dtype="float64")
    labels = frame[label_branch].to_numpy(dtype="int32")
    return {
        "schema": "AUAU_BDT_PPG12_EXACT_WEIGHT_CLOSURE_V1",
        "reused_precomputed_training_weight": True,
        "precomputed_weight_status": status,
        "sample_validation": sample_report,
        "weighting": {
            "weight_mode": "ppg12-exact",
            "event_weight_used": False,
            "vertex_reweight": False,
            "centrality_event_weight": False,
            "cross_section_weight_used_for_training": False,
            "weights_computed_before_binning": True,
            "source": PPG12_EXACT_WEIGHT_COLUMN,
            "sum_weight_class0": float(weights[labels == 0].sum()),
            "sum_weight_class1": float(weights[labels == 1].sum()),
            "min_weight": float(np.min(weights)) if len(weights) else math.nan,
            "max_weight": float(np.max(weights)) if len(weights) else math.nan,
            "mean_weight": float(np.mean(weights)) if len(weights) else math.nan,
        },
    }


def compute_weights(frame, label_branch: str, args) -> tuple[object, dict]:
    import numpy as np

    labels = frame[label_branch].to_numpy(dtype="int32")
    if getattr(args, "weight_mode", "legacy") == "ppg12-exact":
        if PPG12_EXACT_WEIGHT_COLUMN not in frame.columns:
            if getattr(args, "campaign", None):
                raise SystemExit(
                    "PPG12-exact campaign training requires globally precomputed "
                    f"{PPG12_EXACT_WEIGHT_COLUMN} before model slicing."
                )
            weights, diagnostics = compute_ppg12_exact_global_weights(frame, label_branch)
            return weights, diagnostics
        status = require_usable_ppg12_exact_precomputed_weights(frame, "model training")
        weights = frame[PPG12_EXACT_WEIGHT_COLUMN].to_numpy(dtype="float64")
        diagnostics: dict[str, object] = {
            "weight_mode": "ppg12-exact",
            "event_weight_used": False,
            "vertex_reweight": False,
            "centrality_event_weight": False,
            "cross_section_weight_used_for_training": False,
            "weights_computed_before_binning": True,
            "source": PPG12_EXACT_WEIGHT_COLUMN,
            "sum_weight_class0": float(weights[labels == 0].sum()),
            "sum_weight_class1": float(weights[labels == 1].sum()),
            "min_weight": float(np.min(weights)) if len(weights) else math.nan,
            "max_weight": float(np.max(weights)) if len(weights) else math.nan,
            "mean_weight": float(np.mean(weights)) if len(weights) else math.nan,
            "precomputed_weight_status": status,
        }
        return weights, diagnostics

    weights = np.ones(len(frame), dtype="float64")
    diagnostics: dict[str, object] = {
        "weight_mode": "legacy",
        "use_event_weight": bool(args.use_event_weight),
        "et_reweight": bool(args.et_reweight),
        "eta_reweight": bool(args.eta_reweight),
    }

    if args.use_event_weight and args.weight_branch in frame.columns:
        ew = frame[args.weight_branch].to_numpy(dtype="float64")
        weights *= np.where(np.isfinite(ew) & (ew > 0.0), ew, 1.0)

    class_sums: dict[str, float] = {}
    for cls in (0, 1):
        mask = labels == cls
        class_sums[str(cls)] = float(weights[mask].sum())
    target = 0.5 * sum(class_sums.values())
    if target > 0.0:
        for cls in (0, 1):
            mask = labels == cls
            denom = float(weights[mask].sum())
            if denom > 0.0:
                weights[mask] *= target / denom

    if args.et_reweight and "cluster_Et" in frame.columns:
        weights *= inverse_pdf_factors(frame["cluster_Et"].to_numpy(), labels, args.flatten_bins, args.max_flatten_factor)
    if args.eta_reweight and "cluster_Eta" in frame.columns:
        weights *= inverse_pdf_factors(frame["cluster_Eta"].to_numpy(), labels, args.flatten_bins, args.max_flatten_factor)

    finite_positive = np.isfinite(weights) & (weights > 0.0)
    if finite_positive.any():
        median = float(np.median(weights[finite_positive]))
        cap = max(args.max_total_weight_factor * median, 1.0e-12)
        weights = np.where(finite_positive, np.clip(weights, 1.0e-12, cap), 1.0)
    else:
        weights = np.ones(len(frame), dtype="float64")

    diagnostics["sum_weight_class0"] = float(weights[labels == 0].sum())
    diagnostics["sum_weight_class1"] = float(weights[labels == 1].sum())
    diagnostics["min_weight"] = float(np.min(weights)) if len(weights) else math.nan
    diagnostics["max_weight"] = float(np.max(weights)) if len(weights) else math.nan
    diagnostics["mean_weight"] = float(np.mean(weights)) if len(weights) else math.nan
    return weights, diagnostics


def ppg12_background_subsample(frame, label_branch: str, args):
    import numpy as np

    if args.background_subsample_fraction >= 1.0:
        return frame, {"enabled": False}
    if "cluster_Et" not in frame.columns:
        return frame, {"enabled": False, "reason": "missing cluster_Et"}

    bkg = frame[label_branch].to_numpy(dtype="int32") == 0
    et = frame["cluster_Et"].to_numpy(dtype="float64")
    low_bkg = bkg & np.isfinite(et) & (et < args.background_subsample_et_threshold)
    keep = np.ones(len(frame), dtype=bool)

    rng = np.random.default_rng(args.background_subsample_seed)
    if args.background_subsample_flatten and low_bkg.sum() > args.background_subsample_bins:
        edges = np.linspace(float(np.nanmin(et[low_bkg])), float(args.background_subsample_et_threshold), args.background_subsample_bins + 1)
        for ibin in range(args.background_subsample_bins):
            mask = low_bkg & (et >= edges[ibin]) & (et < edges[ibin + 1])
            idx = np.flatnonzero(mask)
            if len(idx) == 0:
                continue
            n_keep = max(1, int(round(len(idx) * args.background_subsample_fraction)))
            drop = np.setdiff1d(idx, rng.choice(idx, size=n_keep, replace=False), assume_unique=False)
            keep[drop] = False
    else:
        idx = np.flatnonzero(low_bkg)
        n_keep = max(1, int(round(len(idx) * args.background_subsample_fraction))) if len(idx) else 0
        if n_keep < len(idx):
            drop = np.setdiff1d(idx, rng.choice(idx, size=n_keep, replace=False), assume_unique=False)
            keep[drop] = False

    return frame.loc[keep].copy(), {
        "enabled": True,
        "fraction": float(args.background_subsample_fraction),
        "et_threshold": float(args.background_subsample_et_threshold),
        "rows_before": int(len(frame)),
        "rows_after": int(keep.sum()),
        "low_et_background_before": int(low_bkg.sum()),
        "low_et_background_after": int((low_bkg & keep).sum()),
    }


def adaptive_majority_cap(frame, label_branch: str, args, metadata: dict):
    import numpy as np

    cap_ratio = float(getattr(args, "majority_cap_ratio", 0.0) or 0.0)
    if cap_ratio <= 0.0 or len(frame) == 0:
        return frame, {"enabled": False}

    labels = frame[label_branch].to_numpy(dtype="int32")
    counts = {cls: int((labels == cls).sum()) for cls in (0, 1)}
    if counts[0] == 0 or counts[1] == 0:
        return frame, {"enabled": False, "reason": "missing class", "counts_before": counts}

    minority = 0 if counts[0] <= counts[1] else 1
    majority = 1 - minority
    target_majority = int(math.ceil(counts[minority] * cap_ratio))
    if counts[majority] <= target_majority:
        return frame, {
            "enabled": False,
            "reason": "ratio already within cap",
            "cap_ratio": cap_ratio,
            "counts_before": counts,
        }

    rng = np.random.default_rng(stable_seed(metadata.get("model_id", "model"), args.random_seed, "majority-cap"))
    minority_idx = np.flatnonzero(labels == minority)
    majority_idx = np.flatnonzero(labels == majority)
    et = frame["cluster_Et"].to_numpy(dtype="float64") if "cluster_Et" in frame.columns else np.full(len(frame), np.nan)
    cent = frame["centrality"].to_numpy(dtype="float64") if "centrality" in frame.columns else np.full(len(frame), np.nan)

    et_bins = max(1, int(getattr(args, "majority_cap_et_bins", 12) or 12))
    cent_bins = max(1, int(getattr(args, "majority_cap_cent_bins", 8) or 8))

    def bin_index(values, idx, nbins):
        finite = np.isfinite(values[idx])
        out = np.zeros(len(idx), dtype="int32")
        if finite.sum() < nbins:
            return out
        lo, hi = np.nanpercentile(values[idx][finite], [1.0, 99.0])
        if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
            return out
        edges = np.linspace(lo, hi, nbins + 1)
        out = np.clip(np.searchsorted(edges, values[idx], side="right") - 1, 0, nbins - 1).astype("int32")
        out[~np.isfinite(values[idx])] = 0
        return out

    et_i = bin_index(et, majority_idx, et_bins)
    cent_i = bin_index(cent, majority_idx, cent_bins)
    strata = et_i * 1000 + cent_i
    keep_majority = []
    remaining = target_majority
    unique, unique_counts = np.unique(strata, return_counts=True)
    allocations = {}
    for s, n in zip(unique, unique_counts):
        allocations[int(s)] = max(1, int(math.floor(target_majority * int(n) / counts[majority])))
    allocated = sum(allocations.values())
    if allocated > target_majority:
        for s in sorted(allocations, key=lambda x: allocations[x], reverse=True):
            if allocated <= target_majority:
                break
            if allocations[s] > 1:
                allocations[s] -= 1
                allocated -= 1
    elif allocated < target_majority:
        for s in sorted(allocations, key=lambda x: allocations[x]):
            if allocated >= target_majority:
                break
            allocations[s] += 1
            allocated += 1

    for s in unique:
        local = majority_idx[strata == s]
        n_keep = min(len(local), allocations.get(int(s), 0), remaining)
        if n_keep <= 0:
            continue
        keep_majority.extend(rng.choice(local, size=n_keep, replace=False).tolist())
        remaining -= n_keep

    if remaining > 0:
        chosen = set(keep_majority)
        rest = np.asarray([i for i in majority_idx if int(i) not in chosen], dtype="int64")
        if len(rest):
            extra = rng.choice(rest, size=min(remaining, len(rest)), replace=False)
            keep_majority.extend(extra.tolist())

    keep = np.asarray(sorted(set(minority_idx.tolist() + keep_majority)), dtype="int64")
    capped = frame.iloc[keep].copy()
    after_labels = capped[label_branch].to_numpy(dtype="int32")
    counts_after = {cls: int((after_labels == cls).sum()) for cls in (0, 1)}
    return capped, {
        "enabled": True,
        "cap_ratio": cap_ratio,
        "minority_class": minority,
        "majority_class": majority,
        "counts_before": counts,
        "counts_after": counts_after,
        "target_majority": target_majority,
        "strata_used": int(len(unique)),
    }


def event_level_train_test_split(frame, x, y, weights, args, metadata: dict):
    import numpy as np

    required = ["run", "evt"]
    missing = [col for col in required if col not in frame.columns]
    if missing:
        raise SystemExit(
            "--split-mode event50 requires event identifier branches in the training frame: "
            + ", ".join(missing)
        )
    if len(frame) != len(x) or len(frame) != len(y) or len(frame) != len(weights):
        raise SystemExit("--split-mode event50 internal length mismatch")

    if len(frame) == 0:
        raise SystemExit("--split-mode event50 received an empty frame")
    event_keys, event_key_columns = global_event_keys(
        frame,
        require_file_qualified=bool(getattr(args, "require_global_event_key", False)),
    )
    unique_keys = np.unique(event_keys)
    if len(unique_keys) < 4:
        raise SystemExit(
            "--split-mode event50 cannot build a stable 50/50 event split: "
            f"only {len(unique_keys)} unique event keys"
        )

    seed = int(getattr(args, "random_seed", 13))
    hashes = np.asarray([stable_seed("event50", seed, metadata.get("model_id", "model"), key) for key in unique_keys], dtype=np.uint64)
    order = np.argsort(hashes, kind="mergesort")
    test_fraction = float(getattr(args, "test_size", 0.5))
    if not (0.0 < test_fraction < 1.0):
        raise SystemExit(f"--split-mode event50 requires 0 < --test-size < 1, got {test_fraction}")
    n_test_events = int(round(len(unique_keys) * test_fraction))
    n_test_events = min(max(n_test_events, 1), len(unique_keys) - 1)
    test_keys = set(unique_keys[order[:n_test_events]].tolist())
    test_mask = np.asarray([key in test_keys for key in event_keys], dtype=bool)
    train_mask = ~test_mask

    def counts(mask):
        return {str(cls): int(np.sum(mask & (y == cls))) for cls in (0, 1)}

    train_counts = counts(train_mask)
    test_counts = counts(test_mask)
    if any(train_counts[str(cls)] <= 0 for cls in (0, 1)) or any(test_counts[str(cls)] <= 0 for cls in (0, 1)):
        raise SystemExit(
            "--split-mode event50 produced a split without both classes in train/test: "
            f"train={train_counts} test={test_counts}"
        )

    report = {
        "mode": "event50",
        "event_key_columns": event_key_columns,
        "file_qualified_event_key": event_key_columns == GLOBAL_EVENT_KEY_COLUMNS,
        "unique_events": int(len(unique_keys)),
        "test_fraction_requested": test_fraction,
        "train_events": int(len(unique_keys) - n_test_events),
        "test_events": int(n_test_events),
        "train_rows": int(train_mask.sum()),
        "test_rows": int(test_mask.sum()),
        "train_class_counts": train_counts,
        "test_class_counts": test_counts,
    }
    return (
        x[train_mask],
        x[test_mask],
        y[train_mask],
        y[test_mask],
        weights[train_mask],
        weights[test_mask],
        report,
    )


def train_one(frame, features: list[str], label_branch: str, output: Path, metadata: dict, args) -> dict | None:
    import csv
    import numpy as np
    from sklearn.metrics import brier_score_loss, log_loss, roc_auc_score
    from sklearn.model_selection import train_test_split
    from xgboost import XGBClassifier

    cols = features + [label_branch]
    frame = frame.loc[finite_mask(frame, cols)].copy()
    frame[label_branch] = frame[label_branch].astype(int)
    frame = frame[frame[label_branch].isin([0, 1])].copy()
    n_rows_after_finite = int(len(frame))
    if getattr(args, "weight_mode", "legacy") == "ppg12-exact":
        subsample_report = {"enabled": False, "reason": "ppg12-exact uses global pre-slice training weights"}
        majority_cap_report = {"enabled": False, "reason": "ppg12-exact disables local majority-cap mutation"}
    else:
        frame, subsample_report = ppg12_background_subsample(frame, label_branch, args)
        frame, majority_cap_report = adaptive_majority_cap(frame, label_branch, args, metadata)

    n_sig = int((frame[label_branch] == 1).sum())
    n_bkg = int((frame[label_branch] == 0).sum())
    enough = n_sig >= args.min_rows_per_class and n_bkg >= args.min_rows_per_class
    if not enough:
        msg = (
            f"Skipping {output.name}: class counts below minimum "
            f"(signal={n_sig}, background={n_bkg}, minimum={args.min_rows_per_class})"
        )
        if metadata.get("cent_range") == "all" and not args.campaign:
            raise SystemExit(msg)
        print(f"[WARN] {msg}")
        skip_report = {
            **metadata,
            "status": "skipped",
            "skip_reason": "class counts below minimum",
            "output_tmva": str(output),
            "output_xgb_json": str(output.with_suffix(".xgb.json")),
            "features": features,
            "label_branch": label_branch,
            "n_rows_after_finite": n_rows_after_finite,
            "n_rows": int(len(frame)),
            "n_signal": n_sig,
            "n_background": n_bkg,
            "min_rows_per_class": int(args.min_rows_per_class),
            "background_subsampling": subsample_report,
            "majority_class_optimization": majority_cap_report,
        }
        if args.campaign:
            output.parent.mkdir(parents=True, exist_ok=True)
            output.with_suffix(".metadata.json").write_text(json.dumps(skip_report, indent=2, sort_keys=True) + "\n")
            return skip_report
        return None

    x = frame[features].to_numpy(dtype="float32")
    y = frame[label_branch].to_numpy(dtype="int32")
    weights, weight_report = compute_weights(frame, label_branch, args)

    split_mode = str(getattr(args, "split_mode", "row"))
    if split_mode == "event50":
        x_train, x_test, y_train, y_test, w_train, w_test, split_report = event_level_train_test_split(
            frame, x, y, weights, args, metadata
        )
    else:
        stratify = y if min(n_sig, n_bkg) >= 2 else None
        x_train, x_test, y_train, y_test, w_train, w_test = train_test_split(
            x, y, weights, test_size=args.test_size, random_state=args.random_seed, stratify=stratify
        )
        split_report = {
            "mode": "row",
            "test_fraction_requested": float(args.test_size),
            "train_rows": int(len(y_train)),
            "test_rows": int(len(y_test)),
        }

    model = XGBClassifier(
        n_estimators=args.n_estimators,
        max_depth=args.max_depth,
        learning_rate=args.learning_rate,
        subsample=args.subsample,
        colsample_bytree=args.colsample_bytree,
        reg_alpha=args.reg_alpha,
        reg_lambda=args.reg_lambda,
        grow_policy=args.grow_policy,
        max_bin=args.max_bin,
        n_jobs=args.n_jobs,
        objective="binary:logistic",
        eval_metric=["auc", "logloss"],
        tree_method=args.tree_method,
        random_state=args.random_seed,
    )
    model.fit(
        x_train,
        y_train,
        sample_weight=w_train,
        eval_set=[(x_train, y_train), (x_test, y_test)],
        sample_weight_eval_set=[w_train, w_test],
        verbose=False,
    )

    train_pred = model.predict_proba(x_train)[:, 1]
    holdout_pred = model.predict_proba(x_test)[:, 1]

    def safe_auc(y_true, pred, weight) -> float:
        return float(roc_auc_score(y_true, pred, sample_weight=weight)) if len(np.unique(y_true)) == 2 else math.nan

    def safe_logloss(y_true, pred, weight) -> float:
        try:
            return float(log_loss(y_true, pred, sample_weight=weight, labels=[0, 1]))
        except ValueError:
            return math.nan

    def safe_brier(y_true, pred, weight) -> float:
        try:
            return float(brier_score_loss(y_true, pred, sample_weight=weight))
        except ValueError:
            return math.nan

    train_auc = safe_auc(y_train, train_pred, w_train)
    holdout_auc = safe_auc(y_test, holdout_pred, w_test)
    train_logloss = safe_logloss(y_train, train_pred, w_train)
    holdout_logloss = safe_logloss(y_test, holdout_pred, w_test)
    train_brier = safe_brier(y_train, train_pred, w_train)
    holdout_brier = safe_brier(y_test, holdout_pred, w_test)
    auc_gap = train_auc - holdout_auc if math.isfinite(train_auc) and math.isfinite(holdout_auc) else math.nan
    logloss_gap = (
        holdout_logloss - train_logloss
        if math.isfinite(holdout_logloss) and math.isfinite(train_logloss)
        else math.nan
    )
    auc = holdout_auc

    output.parent.mkdir(parents=True, exist_ok=True)
    json_model = output.with_suffix(".xgb.json")
    model.get_booster().save_model(json_model)

    eval_history = model.evals_result()
    history_path = output.with_suffix(".training_history.csv")
    train_hist = eval_history.get("validation_0", {})
    holdout_hist = eval_history.get("validation_1", {})
    n_history = max(
        len(train_hist.get("auc", [])),
        len(train_hist.get("logloss", [])),
        len(holdout_hist.get("auc", [])),
        len(holdout_hist.get("logloss", [])),
    )

    def history_value(history: dict, metric: str, index: int) -> float:
        values = history.get(metric, [])
        return values[index] if index < len(values) else math.nan

    with history_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["iteration", "train_auc", "train_logloss", "holdout_auc", "holdout_logloss"],
        )
        writer.writeheader()
        for idx in range(n_history):
            writer.writerow(
                {
                    "iteration": idx,
                    "train_auc": history_value(train_hist, "auc", idx),
                    "train_logloss": history_value(train_hist, "logloss", idx),
                    "holdout_auc": history_value(holdout_hist, "auc", idx),
                    "holdout_logloss": history_value(holdout_hist, "logloss", idx),
                }
            )

    export_status = "not_attempted"
    export_error = ""
    if args.skip_tmva_export:
        export_status = "skipped_by_request"
    else:
        try:
            with TMVA_EXPORT_LOCK:
                import ROOT

                booster = model.get_booster()
                booster.feature_names = [f"f{i}" for i in range(len(features))]
                # ROOT 6.32's SaveXGBoost cannot parse XGBoost >=3 vector-valued
                # base_score strings like "[5E-1]". Patch only the Python config view
                # handed to TMVA; the saved XGBoost JSON remains untouched for audit.
                original_save_config = booster.save_config

                def tmva_save_config():
                    import re

                    text = original_save_config()
                    return re.sub(r'"base_score":"\[([0-9eE+\-.]+)\]"', r'"base_score":"\1"', text)

                booster.save_config = tmva_save_config  # type: ignore[method-assign]
                ROOT.TMVA.Experimental.SaveXGBoost(model, "myBDT", str(output), num_inputs=len(features))
            export_status = "ok"
        except Exception as exc:  # noqa: BLE001
            export_status = f"failed: {exc}"
            export_error = str(exc)

    report = {
        **metadata,
        "status": "trained",
        "output_tmva": str(output),
        "output_xgb_json": str(json_model),
        "features": features,
        "label_branch": label_branch,
        "n_rows_after_finite": n_rows_after_finite,
        "n_rows": int(len(frame)),
        "n_signal": n_sig,
        "n_background": n_bkg,
        "auc": auc,
        "training_history_csv": str(history_path),
        "overfit_diagnostics": {
            "schema": "AUAU_BDT_OVERFIT_DIAGNOSTICS_V1",
            "train_auc": train_auc,
            "holdout_auc": holdout_auc,
            "train_logloss": train_logloss,
            "holdout_logloss": holdout_logloss,
            "train_brier": train_brier,
            "holdout_brier": holdout_brier,
            "auc_gap_train_minus_holdout": auc_gap,
            "logloss_gap_holdout_minus_train": logloss_gap,
            "train_rows": int(len(y_train)),
            "holdout_rows": int(len(y_test)),
            "history_csv": str(history_path),
            "eval_metric": ["auc", "logloss"],
            "eval_history_keys": sorted(eval_history.keys()),
        },
        "tmva_export": export_status,
        "weighting": weight_report,
        "split": split_report,
        "xgboost": {
            "n_estimators": args.n_estimators,
            "max_depth": args.max_depth,
            "learning_rate": args.learning_rate,
            "subsample": args.subsample,
            "colsample_bytree": args.colsample_bytree,
            "tree_method": args.tree_method,
            "reg_alpha": args.reg_alpha,
            "reg_lambda": args.reg_lambda,
            "grow_policy": args.grow_policy,
            "max_bin": args.max_bin,
            "n_jobs": args.n_jobs,
        },
        "background_subsampling": subsample_report,
        "majority_class_optimization": majority_cap_report,
    }
    output.with_suffix(".metadata.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if export_error and not args.allow_tmva_export_failure:
        raise SystemExit(
            f"TMVA export failed for {output}: {export_error}. "
            "The XGBoost JSON was written for diagnostics, but this pipeline needs "
            "the TMVA/RBDT ROOT file for analysis consumption."
        )
    return report


def campaign_specs(args, outdir: Path) -> list[dict]:
    pt_edges = parse_float_edges(args.pt_bins)
    pt_bins = bins_from_edges(pt_edges)
    coarse_cent_bins = parse_cent_bins(args.coarse_cent_bins)
    fine_cent_bins = parse_cent_bins(args.fine_cent_bins)
    specs: list[dict] = []

    def add(
        product: str,
        model_id: str,
        features: list[str],
        pt_range: tuple[float, float] | None = None,
        cent_range: tuple[float, float] | None = None,
        minority_optimized: bool = True,
        role: str = "single",
        majority_cap_ratio: float | None = None,
    ) -> None:
        safe = model_id.replace(".", "p")
        specs.append(
            {
                "model_id": safe,
                "product": product,
                "role": role,
                "features": list(features),
                "pt_range": list(pt_range) if pt_range is not None else None,
                "cent_range": list(cent_range) if cent_range is not None else None,
                "minority_optimized": bool(minority_optimized),
                "output_tmva": str(outdir / f"auau_tight_bdt_{safe}_tmva.root"),
                "output_xgb_json": str(outdir / f"auau_tight_bdt_{safe}_tmva.xgb.json"),
                "metadata": str(outdir / f"auau_tight_bdt_{safe}_tmva.metadata.json"),
                "majority_cap_ratio": majority_cap_ratio,
            }
        )

    base = list(PPG12_TIGHT_FEATURES)
    cent_feat = list(PPG12_TIGHT_FEATURES) + ["centrality"]
    cent_feat_3x3 = list(PPG12_TIGHT_FEATURES_3X3_WIDTHS) + ["centrality"]
    cent_feat_base3x3 = list(PPG12_TIGHT_FEATURES_BASE_AND_3X3_WIDTHS) + ["centrality"]
    cent_feat_width_ratios = list(PPG12_TIGHT_FEATURES) + [
        "cluster_weta33_cogx",
        "cluster_wphi33_cogx",
        *WIDTH_RATIO_FEATURES,
        "centrality",
    ]

    for suffix, pt_range in (("allRange", None), ("pt5to40", (5.0, 40.0))):
        add(f"centINDcontrol_{suffix}", f"centINDcontrol_{suffix}", base, pt_range, None, True)
        add(f"centAsFeat_{suffix}", f"centAsFeat_{suffix}", cent_feat, pt_range, None, True)
        if suffix == "pt5to40":
            add(
                "centAsFeat3x3_pt5to40",
                "centAsFeat3x3_pt5to40",
                cent_feat_3x3,
                pt_range,
                None,
                True,
                "single-3x3-width",
            )
        add(
            f"centAsFeatMinOpt_{suffix}",
            f"centAsFeatMinOpt_{suffix}",
            cent_feat,
            pt_range,
            None,
            True,
            majority_cap_ratio=args.minopt_majority_cap_ratio,
        )
        for clo, chi in coarse_cent_bins:
            add(
                f"centDepBDTs_{suffix}",
                f"centDepBDTs_{suffix}_{cent_tag(clo, chi)}",
                base,
                pt_range,
                (clo, chi),
                True,
                "cent-bin",
            )
        for clo, chi in fine_cent_bins:
            add(
                f"centDepFineBDTs_{suffix}",
                f"centDepFineBDTs_{suffix}_{cent_tag(clo, chi)}",
                base,
                pt_range,
                (clo, chi),
                True,
                "fine-cent-bin",
            )

    for lo, hi in parse_range_list(args.extra_cent_as_feat_pt_ranges):
        suffix = pt_window_tag(lo, hi)
        add(
            f"centAsFeat_{suffix}",
            f"centAsFeat_{suffix}",
            cent_feat,
            (lo, hi),
            None,
            True,
            "single-pt-window",
        )
    for lo, hi in parse_range_list(args.extra_cent_as_feat_3x3_pt_ranges):
        suffix = pt_window_tag(lo, hi)
        add(
            f"centAsFeat3x3_{suffix}",
            f"centAsFeat3x3_{suffix}",
            cent_feat_3x3,
            (lo, hi),
            None,
            True,
            "single-pt-window-3x3-width",
        )
    for lo, hi in parse_range_list(args.extra_cent_as_feat_base3x3_pt_ranges):
        suffix = pt_window_tag(lo, hi)
        add(
            f"centAsFeatBase3x3_{suffix}",
            f"centAsFeatBase3x3_{suffix}",
            cent_feat_base3x3,
            (lo, hi),
            None,
            True,
            "single-pt-window-base-and-3x3-width",
        )
    for lo, hi in parse_range_list(args.extra_cent_as_feat_width_ratio_pt_ranges):
        suffix = pt_window_tag(lo, hi)
        add(
            f"centAsFeatWidthRatios_{suffix}",
            f"centAsFeatWidthRatios_{suffix}",
            cent_feat_width_ratios,
            (lo, hi),
            None,
            True,
            "single-pt-window-width-ratios",
        )

    for plo, phi in pt_bins:
        add(
            "ptBinCentAsFeat",
            f"ptBinCentAsFeat_{pt_tag(plo, phi)}",
            cent_feat,
            (plo, phi),
            None,
            True,
            "pt-bin",
        )
        for clo, chi in coarse_cent_bins:
            add(
                "ptCentDep3",
                f"ptCentDep3_{pt_tag(plo, phi)}_{cent_tag(clo, chi)}",
                base,
                (plo, phi),
                (clo, chi),
                True,
                "pt-cent-bin",
            )
        for clo, chi in fine_cent_bins:
            add(
                "ptCentDepFine",
                f"ptCentDepFine_{pt_tag(plo, phi)}_{cent_tag(clo, chi)}",
                base,
                (plo, phi),
                (clo, chi),
                True,
                "pt-fine-cent-bin",
            )

    return specs


def etfine_centstudy_specs(args, outdir: Path) -> list[dict]:
    pt_edges = parse_float_edges(args.pt_bins)
    pt_bins = bins_from_edges(pt_edges)
    coarse_cent_bins = parse_cent_bins(args.coarse_cent_bins)
    fine_cent_bins = parse_cent_bins(args.fine_cent_bins)
    specs: list[dict] = []

    def add(
        product: str,
        model_id: str,
        features: list[str],
        pt_range: tuple[float, float] | None = None,
        cent_range: tuple[float, float] | None = None,
        role: str = "single",
    ) -> None:
        safe = model_id.replace(".", "p")
        specs.append(
            {
                "model_id": safe,
                "product": product,
                "role": role,
                "features": list(features),
                "pt_range": list(pt_range) if pt_range is not None else None,
                "cent_range": list(cent_range) if cent_range is not None else None,
                "minority_optimized": True,
                "output_tmva": str(outdir / f"auau_tight_bdt_{safe}_tmva.root"),
                "output_xgb_json": str(outdir / f"auau_tight_bdt_{safe}_tmva.xgb.json"),
                "metadata": str(outdir / f"auau_tight_bdt_{safe}_tmva.metadata.json"),
                "majority_cap_ratio": None,
            }
        )

    features_with_cent = list(PPG12_TIGHT_FEATURES_BASE_AND_3X3_WIDTHS) + ["centrality"]
    features_cent_binned = list(PPG12_TIGHT_FEATURES_BASE_AND_3X3_WIDTHS)
    features_no_cent = list(PPG12_TIGHT_FEATURES_BASE_AND_3X3_WIDTHS)
    if len(pt_edges) < 2:
        raise SystemExit("etfine-centstudy needs at least two pT edges")
    full_pt = (pt_edges[0], pt_edges[-1])

    add(
        "noCent_pt1535",
        "noCent_pt1535",
        features_no_cent,
        full_pt,
        None,
        "single-pt-window-base-and-3x3-width-no-cent",
    )
    add(
        "centInput_pt1535",
        "centInput_pt1535",
        features_with_cent,
        full_pt,
        None,
        "single-pt-window-base-and-3x3-width",
    )
    for plo, phi in pt_bins:
        add(
            "ptFine_noCent",
            f"ptFine_noCent_{pt_tag(plo, phi)}",
            features_no_cent,
            (plo, phi),
            None,
            "fine-pt-bin-no-cent",
        )
        add(
            "ptFine_centInput",
            f"ptFine_centInput_{pt_tag(plo, phi)}",
            features_with_cent,
            (plo, phi),
            None,
            "fine-pt-bin-cent-as-feature",
        )
        for clo, chi in coarse_cent_bins:
            add(
                "ptFine_cent3",
                f"ptFine_cent3_{pt_tag(plo, phi)}_{cent_tag(clo, chi)}",
                features_cent_binned,
                (plo, phi),
                (clo, chi),
                "fine-pt-coarse-cent-bin",
            )
        for clo, chi in fine_cent_bins:
            add(
                "ptFine_cent7",
                f"ptFine_cent7_{pt_tag(plo, phi)}_{cent_tag(clo, chi)}",
                features_cent_binned,
                (plo, phi),
                (clo, chi),
                "fine-pt-fine-cent-bin",
            )

    return specs


def isolation_diagnostic_specs(args, outdir: Path) -> list[dict]:
    pt_edges = parse_float_edges(args.pt_bins)
    pt_bins = bins_from_edges(pt_edges)
    fine_cent_bins = parse_cent_bins(args.fine_cent_bins)
    specs: list[dict] = []
    features_full = diagnostic_isolation_feature_family(include_centrality=True)
    warning = (
        "Uses reconstructed isolation-derived inputs. Diagnostic ceiling test only; "
        "not ABCD-safe photon ID and not for purity production without redesign."
    )

    def add(
        product: str,
        model_id: str,
        pt_range: tuple[float, float] | None,
        cent_range: tuple[float, float] | None,
        role: str,
    ) -> None:
        safe = model_id.replace(".", "p")
        specs.append(
            {
                "model_id": safe,
                "product": product,
                "role": role,
                "features": list(features_full),
                "pt_range": list(pt_range) if pt_range is not None else None,
                "cent_range": list(cent_range) if cent_range is not None else None,
                "minority_optimized": True,
                "output_tmva": str(outdir / f"auau_tight_bdt_{safe}_tmva.root"),
                "output_xgb_json": str(outdir / f"auau_tight_bdt_{safe}_tmva.xgb.json"),
                "metadata": str(outdir / f"auau_tight_bdt_{safe}_tmva.metadata.json"),
                "majority_cap_ratio": None,
                "diagnostic_only": True,
                "abcd_warning": warning,
            }
        )

    if len(pt_edges) < 2:
        raise SystemExit("iso-diagnostic needs at least two pT edges")
    full_pt = (pt_edges[0], pt_edges[-1])
    add(
        "isoBDT_global15to35_EtCent_full",
        "isoBDT_global15to35_EtCent_full",
        full_pt,
        None,
        "diagnostic-global-pt-window-et-cent-isolation-inputs",
    )
    for plo, phi in pt_bins:
        for clo, chi in fine_cent_bins:
            add(
                "isoBDT_ptFine15to35_cent7_full",
                f"isoBDT_ptFine15to35_cent7_full_{pt_tag(plo, phi)}_{cent_tag(clo, chi)}",
                (plo, phi),
                (clo, chi),
                "diagnostic-fine-pt-fine-cent-bin-isolation-inputs",
            )
    return specs


def global_sixpack_specs(args, outdir: Path) -> list[dict]:
    specs: list[dict] = []
    full_pt = (15.0, 35.0)

    def add(model_id: str, features: list[str], role: str, diagnostic_only: bool = False) -> None:
        specs.append(
            {
                "model_id": model_id,
                "product": model_id,
                "role": role,
                "features": list(features),
                "pt_range": list(full_pt),
                "cent_range": None,
                "minority_optimized": False,
                "output_tmva": str(outdir / f"auau_tight_bdt_{model_id}_tmva.root"),
                "output_xgb_json": str(outdir / f"auau_tight_bdt_{model_id}_tmva.xgb.json"),
                "metadata": str(outdir / f"auau_tight_bdt_{model_id}_tmva.metadata.json"),
                "majority_cap_ratio": None,
                "diagnostic_only": diagnostic_only,
                "abcd_warning": "Uses reconstructed isolation-derived inputs; diagnostic only." if diagnostic_only else None,
            }
        )

    add("globalEtCent1535_bdt_noIso", global_sixpack_noiso_features(), "global-15-35-et-cent-full-shower-no-iso")
    add(
        "globalEtCent1535_bdt_iso",
        global_sixpack_iso_features(),
        "global-15-35-et-cent-full-shower-isolation-inputs",
        diagnostic_only=True,
    )
    return specs


def basev3e_e22_ablation_specs(args, outdir: Path) -> list[dict]:
    specs: list[dict] = []
    full_pt = (15.0, 35.0)
    base = list(PPG12_TIGHT_FEATURES_BASE_AND_3X3_WIDTHS)
    if "centrality" not in base:
        base.append("centrality")

    def add(model_id: str, extra_features: list[str], role: str) -> None:
        features = list(base)
        for feature in extra_features:
            if feature not in features:
                features.append(feature)
        specs.append(
            {
                "model_id": model_id,
                "product": model_id,
                "role": role,
                "features": features,
                "pt_range": list(full_pt),
                "cent_range": None,
                "minority_optimized": False,
                "output_tmva": str(outdir / f"auau_tight_bdt_{model_id}_tmva.root"),
                "output_xgb_json": str(outdir / f"auau_tight_bdt_{model_id}_tmva.xgb.json"),
                "metadata": str(outdir / f"auau_tight_bdt_{model_id}_tmva.metadata.json"),
                "majority_cap_ratio": None,
                "diagnostic_only": False,
                "abcd_warning": None,
            }
        )

    add(
        "baseBDT_v3E_withCentrality_w33_E22E37",
        ["e22_over_e37"],
        "base-v3e-centrality-3x3-widths-plus-e22-over-e37",
    )
    add(
        "baseBDT_v3E_withCentrality_w33_E22E53",
        ["e22_over_e53"],
        "base-v3e-centrality-3x3-widths-plus-e22-over-e53",
    )
    add(
        "baseBDT_v3E_withCentrality_w33_E22E37_E22E53",
        ["e22_over_e37", "e22_over_e53"],
        "base-v3e-centrality-3x3-widths-plus-e22-over-e37-and-e22-over-e53",
    )
    return specs


def corrected_baseline_shower_ladder_specs(args, outdir: Path) -> list[dict]:
    specs: list[dict] = []
    full_pt = (15.0, 35.0)
    base11 = list(PPG12_TIGHT_FEATURES)
    base12 = list(PPG12_TIGHT_FEATURES) + ["centrality"]
    base14 = list(PPG12_TIGHT_FEATURES_BASE_AND_3X3_WIDTHS) + ["centrality"]
    full32 = global_sixpack_noiso_features()
    ladder_features = [feature for feature in full32 if feature not in base14]

    def safe_token(text: str) -> str:
        return "".join(ch if ch.isalnum() else "_" for ch in text).strip("_")

    def unique_features(features: list[str]) -> list[str]:
        out: list[str] = []
        for feature in features:
            if feature not in out:
                out.append(feature)
        return out

    def add(
        product: str,
        model_id: str,
        features: list[str],
        role: str,
        diagnostic_only: bool = False,
        warning: str | None = None,
    ) -> None:
        safe = model_id.replace(".", "p")
        specs.append(
            {
                "model_id": safe,
                "product": product,
                "role": role,
                "features": unique_features(features),
                "pt_range": list(full_pt),
                "cent_range": None,
                "minority_optimized": False,
                "output_tmva": str(outdir / f"auau_tight_bdt_{safe}_tmva.root"),
                "output_xgb_json": str(outdir / f"auau_tight_bdt_{safe}_tmva.xgb.json"),
                "metadata": str(outdir / f"auau_tight_bdt_{safe}_tmva.metadata.json"),
                "majority_cap_ratio": None,
                "diagnostic_only": diagnostic_only,
                "abcd_warning": warning,
            }
        )

    add(
        "baseV3E11_pt1535",
        "baseV3E11_pt1535",
        base11,
        "base-v3E-11-feature-control-no-centrality-no-3x3-widths",
    )
    add(
        "baseV3E11_cent12_pt1535",
        "baseV3E11_cent12_pt1535",
        base12,
        "base-v3E-12-feature-control-centrality-no-3x3-widths",
    )
    add(
        "centAsFeatBase3x3_pt15to35",
        "centAsFeatBase3x3_pt15to35",
        base14,
        "corrected-default-14-feature-reference",
    )
    for feature in ladder_features:
        token = safe_token(feature)
        add(
            f"base14_plus1_{token}_pt1535",
            f"base14_plus1_{token}_pt1535",
            base14 + [feature],
            f"corrected-default-14-plus-single-shower-feature:{feature}",
        )
    cumulative = list(base14)
    for index, feature in enumerate(ladder_features, start=1):
        cumulative.append(feature)
        token = safe_token(feature)
        model_id = (
            "globalEtCent1535_bdt_noIso"
            if index == len(ladder_features)
            else f"base14_ladder{index:02d}_{token}_pt1535"
        )
        product = (
            "globalEtCent1535_bdt_noIso"
            if index == len(ladder_features)
            else "base14_to32_cumulative_ladder"
        )
        add(
            product,
            model_id,
            cumulative,
            f"corrected-default-cumulative-shower-ladder-step-{index:02d}:{feature}",
        )
    return specs


def corrected_baseline_binned14_specs(args, outdir: Path) -> list[dict]:
    pt_edges = parse_float_edges(args.pt_bins)
    pt_bins = bins_from_edges(pt_edges)
    coarse_cent_bins = parse_cent_bins(args.coarse_cent_bins)
    fine_cent_bins = parse_cent_bins(args.fine_cent_bins)
    features = list(PPG12_TIGHT_FEATURES_BASE_AND_3X3_WIDTHS) + ["centrality"]
    specs: list[dict] = []

    def add(
        product: str,
        model_id: str,
        pt_range: tuple[float, float] | None,
        cent_range: tuple[float, float] | None,
        role: str,
    ) -> None:
        safe = model_id.replace(".", "p")
        specs.append(
            {
                "model_id": safe,
                "product": product,
                "role": role,
                "features": list(features),
                "pt_range": list(pt_range) if pt_range is not None else None,
                "cent_range": list(cent_range) if cent_range is not None else None,
                "minority_optimized": False,
                "output_tmva": str(outdir / f"auau_tight_bdt_{safe}_tmva.root"),
                "output_xgb_json": str(outdir / f"auau_tight_bdt_{safe}_tmva.xgb.json"),
                "metadata": str(outdir / f"auau_tight_bdt_{safe}_tmva.metadata.json"),
                "majority_cap_ratio": None,
                "diagnostic_only": False,
                "abcd_warning": None,
            }
        )

    if len(pt_edges) < 2:
        raise SystemExit("corrected-baseline-binned14 needs at least two pT edges")
    full_pt = (pt_edges[0], pt_edges[-1])
    for plo, phi in pt_bins:
        add(
            "base14_perEt",
            f"base14_perEt_{pt_tag(plo, phi)}",
            (plo, phi),
            None,
            "corrected-default-14-feature-fine-et-bin",
        )
    for clo, chi in coarse_cent_bins:
        add(
            "base14_perCent3",
            f"base14_perCent3_{cent_tag(clo, chi)}",
            full_pt,
            (clo, chi),
            "corrected-default-14-feature-coarse-centrality-bin",
        )
    for clo, chi in fine_cent_bins:
        add(
            "base14_perCent7",
            f"base14_perCent7_{cent_tag(clo, chi)}",
            full_pt,
            (clo, chi),
            "corrected-default-14-feature-fine-centrality-bin",
        )
    for plo, phi in pt_bins:
        for clo, chi in coarse_cent_bins:
            add(
                "base14_perEtCent3",
                f"base14_perEtCent3_{pt_tag(plo, phi)}_{cent_tag(clo, chi)}",
                (plo, phi),
                (clo, chi),
                "corrected-default-14-feature-fine-et-coarse-centrality-bin",
            )
        for clo, chi in fine_cent_bins:
            add(
                "base14_perEtCent7",
                f"base14_perEtCent7_{pt_tag(plo, phi)}_{cent_tag(clo, chi)}",
                (plo, phi),
                (clo, chi),
                "corrected-default-14-feature-fine-et-fine-centrality-bin",
            )
    return specs


def corrected_baseline_iso14_specs(args, outdir: Path) -> list[dict]:
    specs: list[dict] = []
    full_pt = (15.0, 35.0)
    base14 = list(PPG12_TIGHT_FEATURES_BASE_AND_3X3_WIDTHS) + ["centrality"]
    warning = (
        "Uses raw reconstructed isolation ET as a BDT input. Diagnostic only; "
        "not ABCD-safe photon ID without a separate purity redesign."
    )

    def add(model_id: str, extra_features: list[str], role: str) -> None:
        specs.append(
            {
                "model_id": model_id,
                "product": model_id,
                "role": role,
                "features": list(base14) + list(extra_features),
                "pt_range": list(full_pt),
                "cent_range": None,
                "minority_optimized": False,
                "output_tmva": str(outdir / f"auau_tight_bdt_{model_id}_tmva.root"),
                "output_xgb_json": str(outdir / f"auau_tight_bdt_{model_id}_tmva.xgb.json"),
                "metadata": str(outdir / f"auau_tight_bdt_{model_id}_tmva.metadata.json"),
                "majority_cap_ratio": None,
                "diagnostic_only": True,
                "abcd_warning": warning,
            }
        )

    add("base14_eisoR30_pt1535", ["reco_eiso_r30"], "corrected-default-14-plus-raw-etiso-r03")
    add("base14_eisoR40_pt1535", ["reco_eiso_r40"], "corrected-default-14-plus-raw-etiso-r04")
    add(
        "base14_eisoR30R40_pt1535",
        ["reco_eiso_r30", "reco_eiso_r40"],
        "corrected-default-14-plus-raw-etiso-r03-and-r04",
    )
    return specs


def etcent_binned_sixpack_specs(args, outdir: Path) -> list[dict]:
    pt_edges = parse_float_edges(args.pt_bins)
    pt_bins = bins_from_edges(pt_edges)
    coarse_cent_bins = parse_cent_bins(args.coarse_cent_bins)
    fine_cent_bins = parse_cent_bins(args.fine_cent_bins)
    specs: list[dict] = []
    noiso_features = global_sixpack_noiso_features()
    iso_features = global_sixpack_iso_features()

    def add(
        product: str,
        model_id: str,
        features: list[str],
        pt_range: tuple[float, float],
        cent_range: tuple[float, float],
        role: str,
        diagnostic_only: bool,
    ) -> None:
        safe = model_id.replace(".", "p")
        specs.append(
            {
                "model_id": safe,
                "product": product,
                "role": role,
                "features": list(features),
                "pt_range": list(pt_range),
                "cent_range": list(cent_range),
                "minority_optimized": False,
                "output_tmva": str(outdir / f"auau_tight_bdt_{safe}_tmva.root"),
                "output_xgb_json": str(outdir / f"auau_tight_bdt_{safe}_tmva.xgb.json"),
                "metadata": str(outdir / f"auau_tight_bdt_{safe}_tmva.metadata.json"),
                "majority_cap_ratio": None,
                "diagnostic_only": diagnostic_only,
                "abcd_warning": "Uses reconstructed isolation-derived inputs; diagnostic only." if diagnostic_only else None,
            }
        )

    if len(pt_edges) < 2:
        raise SystemExit("etcent-binned-sixpack needs at least two pT edges")

    groups = (
        (
            "globalEtCent1535_bdt_noIso_ptCent3",
            "globalEtCent1535_bdt_noIso_ptCent3",
            noiso_features,
            coarse_cent_bins,
            "fine-et-coarse-cent-bin-global-sixpack-no-iso-features",
            False,
        ),
        (
            "globalEtCent1535_bdt_noIso_ptCent7",
            "globalEtCent1535_bdt_noIso_ptCent7",
            noiso_features,
            fine_cent_bins,
            "fine-et-fine-cent-bin-global-sixpack-no-iso-features",
            False,
        ),
        (
            "globalEtCent1535_bdt_iso_ptCent3",
            "globalEtCent1535_bdt_iso_ptCent3",
            iso_features,
            coarse_cent_bins,
            "fine-et-coarse-cent-bin-global-sixpack-isolation-input-features",
            True,
        ),
        (
            "globalEtCent1535_bdt_iso_ptCent7",
            "globalEtCent1535_bdt_iso_ptCent7",
            iso_features,
            fine_cent_bins,
            "fine-et-fine-cent-bin-global-sixpack-isolation-input-features",
            True,
        ),
    )
    for product, prefix, features, cent_bins, role, diagnostic_only in groups:
        for plo, phi in pt_bins:
            for clo, chi in cent_bins:
                add(
                    product,
                    f"{prefix}_{pt_tag(plo, phi)}_{cent_tag(clo, chi)}",
                    features,
                    (plo, phi),
                    (clo, chi),
                    role,
                    diagnostic_only,
                )
    return specs


def etcent_binned_sixpack_noiso_ptcent7_specs(args, outdir: Path) -> list[dict]:
    specs = [
        spec
        for spec in etcent_binned_sixpack_specs(args, outdir)
        if spec.get("product") == "globalEtCent1535_bdt_noIso_ptCent7"
    ]
    expected = (len(parse_float_edges(args.pt_bins)) - 1) * len(parse_cent_bins(args.fine_cent_bins))
    if len(specs) != expected:
        raise SystemExit(
            "etcent-binned-sixpack-noiso-ptcent7 planned an unexpected model count: "
            f"{len(specs)} vs expected {expected}"
        )
    return specs


def global_and_etcent_binned_sixpack_noiso_specs(args, outdir: Path) -> list[dict]:
    global_specs = [
        spec
        for spec in global_sixpack_specs(args, outdir)
        if spec.get("product") == "globalEtCent1535_bdt_noIso"
    ]
    routed_specs = [
        spec
        for spec in etcent_binned_sixpack_specs(args, outdir)
        if spec.get("product") in {"globalEtCent1535_bdt_noIso_ptCent3", "globalEtCent1535_bdt_noIso_ptCent7"}
    ]
    expected = 1 + (len(parse_float_edges(args.pt_bins)) - 1) * (
        len(parse_cent_bins(args.coarse_cent_bins)) + len(parse_cent_bins(args.fine_cent_bins))
    )
    specs = global_specs + routed_specs
    if len(specs) != expected:
        raise SystemExit(
            "global-and-etcent-binned-sixpack-noiso planned an unexpected model count: "
            f"{len(specs)} vs expected {expected}"
        )
    return specs


def etcent_binned_eiso_cone_ablation_specs(args, outdir: Path) -> list[dict]:
    requested_pt_edges = parse_float_edges(args.pt_bins)
    if requested_pt_edges != EISO_CONE_ABLATION_PT_EDGES:
        print(
            "[WARN] etcent-binned-eiso-cone-ablation ignores --pt-bins="
            f"{args.pt_bins!r}; using fixed 15-35 GeV edges "
            f"{','.join(f'{edge:g}' for edge in EISO_CONE_ABLATION_PT_EDGES)}"
        )
    pt_edges = list(EISO_CONE_ABLATION_PT_EDGES)
    pt_bins = bins_from_edges(pt_edges)
    coarse_cent_bins = list(EISO_CONE_ABLATION_COARSE_CENT_BINS)
    fine_cent_bins = list(EISO_CONE_ABLATION_FINE_CENT_BINS)
    specs: list[dict] = []
    eiso_r30_features = global_sixpack_eiso_r30_features()
    eiso_r40_features = global_sixpack_eiso_r40_features()
    eiso_r30r40_features = global_sixpack_eiso_r30r40_features()

    def add(
        product: str,
        model_id: str,
        features: list[str],
        pt_range: tuple[float, float],
        cent_range: tuple[float, float],
        role: str,
    ) -> None:
        safe = model_id.replace(".", "p")
        specs.append(
            {
                "model_id": safe,
                "product": product,
                "role": role,
                "features": list(features),
                "pt_range": list(pt_range),
                "cent_range": list(cent_range),
                "minority_optimized": False,
                "output_tmva": str(outdir / f"auau_tight_bdt_{safe}_tmva.root"),
                "output_xgb_json": str(outdir / f"auau_tight_bdt_{safe}_tmva.xgb.json"),
                "metadata": str(outdir / f"auau_tight_bdt_{safe}_tmva.metadata.json"),
                "majority_cap_ratio": None,
                "diagnostic_only": True,
                "abcd_warning": (
                    "Uses raw reconstructed isolation ET as an input; diagnostic only. "
                    "Not ABCD-safe photon ID without a separate purity redesign."
                ),
            }
        )

    if len(pt_edges) < 2:
        raise SystemExit("etcent-binned-eiso-cone-ablation needs at least two pT edges")

    groups = (
        (
            "globalEtCent1535_bdt_eisoR30_ptCent3",
            "globalEtCent1535_bdt_eisoR30_ptCent3",
            eiso_r30_features,
            coarse_cent_bins,
            "fine-et-coarse-cent-bin-global-sixpack-plus-raw-r30-reco-eiso",
        ),
        (
            "globalEtCent1535_bdt_eisoR30_ptCent7",
            "globalEtCent1535_bdt_eisoR30_ptCent7",
            eiso_r30_features,
            fine_cent_bins,
            "fine-et-fine-cent-bin-global-sixpack-plus-raw-r30-reco-eiso",
        ),
        (
            "globalEtCent1535_bdt_eisoR40_ptCent3",
            "globalEtCent1535_bdt_eisoR40_ptCent3",
            eiso_r40_features,
            coarse_cent_bins,
            "fine-et-coarse-cent-bin-global-sixpack-plus-raw-r40-reco-eiso",
        ),
        (
            "globalEtCent1535_bdt_eisoR40_ptCent7",
            "globalEtCent1535_bdt_eisoR40_ptCent7",
            eiso_r40_features,
            fine_cent_bins,
            "fine-et-fine-cent-bin-global-sixpack-plus-raw-r40-reco-eiso",
        ),
        (
            "globalEtCent1535_bdt_eisoR30R40_ptCent3",
            "globalEtCent1535_bdt_eisoR30R40_ptCent3",
            eiso_r30r40_features,
            coarse_cent_bins,
            "fine-et-coarse-cent-bin-global-sixpack-plus-raw-r30-and-r40-reco-eiso",
        ),
        (
            "globalEtCent1535_bdt_eisoR30R40_ptCent7",
            "globalEtCent1535_bdt_eisoR30R40_ptCent7",
            eiso_r30r40_features,
            fine_cent_bins,
            "fine-et-fine-cent-bin-global-sixpack-plus-raw-r30-and-r40-reco-eiso",
        ),
    )
    for product, prefix, features, cent_bins, role in groups:
        for plo, phi in pt_bins:
            for clo, chi in cent_bins:
                add(
                    product,
                    f"{prefix}_{pt_tag(plo, phi)}_{cent_tag(clo, chi)}",
                    features,
                    (plo, phi),
                    (clo, chi),
                    role,
                )
    observed_counts: dict[str, int] = {}
    for spec in specs:
        observed_counts[spec["product"]] = observed_counts.get(spec["product"], 0) + 1
    expected_total = sum(EISO_CONE_ABLATION_EXPECTED_COUNTS.values())
    if len(specs) != expected_total or observed_counts != EISO_CONE_ABLATION_EXPECTED_COUNTS:
        raise SystemExit(
            "Internal raw-eiso campaign spec mismatch: "
            f"expected {EISO_CONE_ABLATION_EXPECTED_COUNTS} ({expected_total} total), "
            f"observed {observed_counts} ({len(specs)} total)"
        )
    return specs


def shape_residual_ptcent7_specs(args, outdir: Path) -> list[dict]:
    pt_edges = list(SHAPE_TEMPLATE_PT_EDGES)
    pt_bins = bins_from_edges(pt_edges)
    fine_cent_bins = parse_cent_bins(args.fine_cent_bins)
    specs: list[dict] = []
    full_features = global_sixpack_noiso_features()
    compact_features = list(PPG12_TIGHT_FEATURES_BASE_AND_3X3_WIDTHS)
    if "centrality" not in compact_features:
        compact_features.append("centrality")
    residual_features = list(SHAPE_RESIDUAL_FEATURES)
    template_features = list(SHAPE_TEMPLATE_FEATURES)

    variants = [
        (
            "globalEtCent1535_bdt_noIso_ptCent7_shapeResiduals",
            "globalEtCent1535_bdt_noIso_ptCent7_shapeResiduals",
            full_features + residual_features,
            "current-56-bdt-full-features-plus-tail-orientation-core-tension-residuals",
        ),
        (
            "globalEtCent1535_bdt_noIso_ptCent7_shapeTemplateDiag",
            "globalEtCent1535_bdt_noIso_ptCent7_shapeTemplateDiag",
            full_features + template_features,
            "current-56-bdt-full-features-plus-diagonal-photon-template-residual",
        ),
        (
            "globalEtCent1535_bdt_noIso_ptCent7_shapeTemplateAll",
            "globalEtCent1535_bdt_noIso_ptCent7_shapeTemplateAll",
            full_features + residual_features + template_features,
            "current-56-bdt-full-features-plus-all-abcd-safe-shower-residuals",
        ),
        (
            "baseV3E_w33_cent_ptCent7_shapeTemplateDiag",
            "baseV3E_w33_cent_ptCent7_shapeTemplateDiag",
            compact_features + template_features,
            "compact-base-v3e-plus-3x3-centrality-and-template-residual",
        ),
        (
            "baseV3E_w33_cent_ptCent7_shapeTemplateAll",
            "baseV3E_w33_cent_ptCent7_shapeTemplateAll",
            compact_features + residual_features + template_features,
            "compact-base-v3e-plus-3x3-centrality-and-all-shower-residuals",
        ),
    ]

    def unique_features(features: list[str]) -> list[str]:
        out: list[str] = []
        for feature in features:
            if feature not in out:
                out.append(feature)
        return out

    if len(pt_edges) < 2:
        raise SystemExit("shape-residual-ptcent7 needs at least two pT edges")
    for product, prefix, features, role in variants:
        deduped = unique_features(features)
        for plo, phi in pt_bins:
            for clo, chi in fine_cent_bins:
                model_id = f"{prefix}_{pt_tag(plo, phi)}_{cent_tag(clo, chi)}"
                specs.append(
                    {
                        "model_id": model_id.replace(".", "p"),
                        "product": product,
                        "role": role,
                        "features": list(deduped),
                        "pt_range": [plo, phi],
                        "cent_range": [clo, chi],
                        "minority_optimized": False,
                        "output_tmva": str(outdir / f"auau_tight_bdt_{model_id.replace('.', 'p')}_tmva.root"),
                        "output_xgb_json": str(outdir / f"auau_tight_bdt_{model_id.replace('.', 'p')}_tmva.xgb.json"),
                        "metadata": str(outdir / f"auau_tight_bdt_{model_id.replace('.', 'p')}_tmva.metadata.json"),
                        "majority_cap_ratio": None,
                        "diagnostic_only": True,
                        "abcd_warning": (
                            "ABCD-safe shower-shape-only sidecar. Template residuals are "
                            "computed from signal-template statistics in the Python validation path; "
                            "freeze constants and wire C++ runtime before production promotion."
                        ),
                    }
                )
    return specs


def ppg12_pp_sixpack_specs(args, outdir: Path) -> list[dict]:
    specs: list[dict] = []
    pt_edges = parse_float_edges(args.pt_bins)
    if len(pt_edges) < 2:
        raise SystemExit("ppg12-sixpack needs at least two --pt-bins edges")
    full_pt = (float(pt_edges[0]), float(pt_edges[-1]))
    pp_cent = (-1.0, 0.0)

    def add(model_id: str, features: list[str], role: str, diagnostic_only: bool = False) -> None:
        specs.append(
            {
                "model_id": model_id,
                "product": model_id,
                "role": role,
                "features": list(features),
                "pt_range": list(full_pt),
                "cent_range": list(pp_cent),
                "minority_optimized": False,
                "output_tmva": str(outdir / f"pp_tight_bdt_{model_id}_tmva.root"),
                "output_xgb_json": str(outdir / f"pp_tight_bdt_{model_id}_tmva.xgb.json"),
                "metadata": str(outdir / f"pp_tight_bdt_{model_id}_tmva.metadata.json"),
                "majority_cap_ratio": None,
                "diagnostic_only": diagnostic_only,
                "abcd_warning": "Uses reconstructed isolation-derived inputs; diagnostic only." if diagnostic_only else None,
            }
        )

    add(
        "ppg12_base_v1E_bdt_noIso",
        PPG12_BASE_V1E_FEATURES,
        "pp-ppg12-base-v1E-no-centrality-no-isolation",
    )
    add(
        "ppg12_base_v3E_bdt_noIso",
        PPG12_TIGHT_FEATURES,
        "pp-ppg12-base-v3E-no-centrality-no-isolation",
    )
    add(
        "ppg12_base_v1E_bdt_iso",
        PPG12_BASE_V1E_FEATURES + ISOLATION_DIAGNOSTIC_FEATURES,
        "pp-ppg12-base-v1E-no-centrality-with-isolation-diagnostics",
        diagnostic_only=True,
    )
    return specs


def registry_payload(specs: list[dict], reports: list[dict], args, status: str = "PLANNED") -> dict:
    report_by_id = {r.get("model_id"): r for r in reports}
    products: dict[str, list[str]] = {}
    for spec in specs:
        products.setdefault(spec["product"], []).append(spec["model_id"])
    if args.campaign == "etcent-binned-eiso-cone-ablation":
        pt_bins = list(EISO_CONE_ABLATION_PT_EDGES)
        coarse_cent_bins = [[lo, hi] for lo, hi in EISO_CONE_ABLATION_COARSE_CENT_BINS]
        fine_cent_bins = [[lo, hi] for lo, hi in EISO_CONE_ABLATION_FINE_CENT_BINS]
    else:
        pt_bins = parse_float_edges(args.pt_bins)
        coarse_cent_bins = [[lo, hi] for lo, hi in parse_cent_bins(args.coarse_cent_bins)]
        fine_cent_bins = [[lo, hi] for lo, hi in parse_cent_bins(args.fine_cent_bins)]
    return {
        "schema": "AUAU_TIGHT_BDT_EXPANDED_REGISTRY_V1",
        "status": status,
        "campaign": args.campaign,
        "expected_model_count": len(specs),
        "model_count": len(specs),
        "products": products,
        "pt_bins": pt_bins,
        "coarse_cent_bins": coarse_cent_bins,
        "fine_cent_bins": fine_cent_bins,
        "defaults": {
            "weight_mode": getattr(args, "weight_mode", "legacy"),
            "majority_cap_ratio": float(args.majority_cap_ratio),
            "minopt_majority_cap_ratio": float(args.minopt_majority_cap_ratio),
            "parallel_workers": int(args.parallel_workers),
            "xgboost_n_jobs": int(args.n_jobs),
            "no_cross_section_weights": getattr(args, "weight_mode", "legacy") == "ppg12-exact",
            "event_weight_used_for_training": (
                bool(getattr(args, "use_event_weight", False))
                and getattr(args, "weight_mode", "legacy") != "ppg12-exact"
            ),
            "ppg12_exact_expected_samples": (
                list(parse_expected_samples(getattr(args, "ppg12_exact_expected_samples", "")))
                if getattr(args, "weight_mode", "legacy") == "ppg12-exact"
                else None
            ),
            "ppg12_exact_closure_dir": (
                str(getattr(args, "ppg12_exact_closure_dir", "") or "")
                if getattr(args, "weight_mode", "legacy") == "ppg12-exact"
                else None
            ),
            "ppg12_exact_closure_artifacts": (
                str(getattr(args, "ppg12_exact_closure_artifacts", "full") or "full")
                if getattr(args, "weight_mode", "legacy") == "ppg12-exact"
                else None
            ),
            "event_quality_filter_enabled": bool(getattr(args, "event_quality_cut_json", None)),
            "event_quality_cut_json": (
                str(getattr(args, "event_quality_cut_json", "") or "")
                if getattr(args, "event_quality_cut_json", None)
                else None
            ),
            "require_global_event_key": bool(getattr(args, "require_global_event_key", False)),
            "random_seed": int(getattr(args, "random_seed", 13)),
            "n_estimators": int(getattr(args, "n_estimators", 450)),
            "max_depth": int(getattr(args, "max_depth", 4)),
            "learning_rate": float(getattr(args, "learning_rate", 0.035)),
            "subsample": float(getattr(args, "subsample", 0.85)),
            "colsample_bytree": float(getattr(args, "colsample_bytree", 0.85)),
            "tree_method": str(getattr(args, "tree_method", "hist")),
            "reg_alpha": float(getattr(args, "reg_alpha", 5.0)),
            "reg_lambda": float(getattr(args, "reg_lambda", 0.3)),
            "grow_policy": str(getattr(args, "grow_policy", "lossguide")),
            "max_bin": int(getattr(args, "max_bin", 256)),
        },
        "models": [{**spec, "report": report_by_id.get(spec["model_id"])} for spec in specs],
    }


def write_registry(path: Path, specs: list[dict], reports: list[dict], args, status: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(registry_payload(specs, reports, args, status), indent=2, sort_keys=True) + "\n")
    print(f"[OK] wrote {path}")


def filter_specs(specs: list[dict], args) -> list[dict]:
    ids = None
    if args.campaign_spec_list is not None:
        ids = {
            line.strip()
            for line in args.campaign_spec_list.read_text().splitlines()
            if line.strip() and not line.strip().startswith("#")
        }
    if args.campaign_spec_ids:
        ids = set(ids or set())
        ids.update(item.strip() for item in args.campaign_spec_ids.split(",") if item.strip())
    if ids is not None:
        specs = [spec for spec in specs if spec["model_id"] in ids]
        missing = sorted(ids - {spec["model_id"] for spec in specs})
        if missing:
            raise SystemExit("Unknown campaign spec ids:\n  " + "\n  ".join(missing[:40]))
    max_specs = int(args.campaign_max_specs or os.environ.get("RJ_AUAU_BDT_CAMPAIGN_MAX_SPECS", "0") or "0")
    if max_specs > 0:
        specs = specs[:max_specs]
    return specs


def frame_for_spec(frame, spec: dict):
    mask = None
    if spec.get("pt_range") is not None:
        lo, hi = spec["pt_range"]
        local = (frame["cluster_Et"] >= lo) & (frame["cluster_Et"] < hi)
        mask = local if mask is None else (mask & local)
    if spec.get("cent_range") is not None:
        lo, hi = spec["cent_range"]
        local = (frame["centrality"] >= lo) & (frame["centrality"] < hi)
        mask = local if mask is None else (mask & local)
    if mask is None:
        return frame
    return frame.loc[mask].copy()


def train_campaign_spec(frame, spec: dict, common_metadata: dict, args) -> dict | None:
    local = frame_for_spec(frame, spec)
    local_args = copy.copy(args)
    if spec.get("majority_cap_ratio") is not None:
        local_args.majority_cap_ratio = float(spec["majority_cap_ratio"])
    metadata = {
        **common_metadata,
        "model_id": spec["model_id"],
        "product": spec["product"],
        "role": spec["role"],
        "pt_range": spec.get("pt_range"),
        "cent_range": spec.get("cent_range"),
        "minority_optimized": spec.get("minority_optimized", False),
        "majority_cap_ratio": float(local_args.majority_cap_ratio),
        "campaign": args.campaign,
        "diagnostic_only": bool(spec.get("diagnostic_only", False)),
        "abcd_warning": spec.get("abcd_warning"),
    }
    output = Path(spec["output_tmva"])
    return train_one(local, spec["features"], common_metadata["label_branch"], output, metadata, local_args)


def run_campaign(args) -> int:
    label_branch = args.label_branch or "is_signal"
    if args.task != "tight":
        raise SystemExit("Expanded campaign currently supports --task tight only.")

    optional_columns = [args.weight_branch]
    if not args.input and not args.plan_only and not (args.cache_file and args.cache_file.is_file()):
        raise SystemExit("--input is required for campaign training unless an existing --cache-file is supplied")
    paths = expand_input_paths(args.input) if args.input else []
    if args.campaign == "etfine-centstudy":
        specs_all = etfine_centstudy_specs(args, args.outdir)
    elif args.campaign == "iso-diagnostic":
        specs_all = isolation_diagnostic_specs(args, args.outdir)
    elif args.campaign == "global-sixpack":
        specs_all = global_sixpack_specs(args, args.outdir)
    elif args.campaign == "basev3e-e22-ablation":
        specs_all = basev3e_e22_ablation_specs(args, args.outdir)
    elif args.campaign == "corrected-baseline-shower-ladder":
        specs_all = corrected_baseline_shower_ladder_specs(args, args.outdir)
    elif args.campaign == "corrected-baseline-binned14":
        specs_all = corrected_baseline_binned14_specs(args, args.outdir)
    elif args.campaign == "corrected-baseline-iso14":
        specs_all = corrected_baseline_iso14_specs(args, args.outdir)
    elif args.campaign == "etcent-binned-sixpack":
        specs_all = etcent_binned_sixpack_specs(args, args.outdir)
    elif args.campaign == "etcent-binned-sixpack-noiso-ptcent7":
        specs_all = etcent_binned_sixpack_noiso_ptcent7_specs(args, args.outdir)
    elif args.campaign == "global-and-etcent-binned-sixpack-noiso":
        specs_all = global_and_etcent_binned_sixpack_noiso_specs(args, args.outdir)
    elif args.campaign == "etcent-binned-eiso-cone-ablation":
        specs_all = etcent_binned_eiso_cone_ablation_specs(args, args.outdir)
    elif args.campaign == "shape-residual-ptcent7":
        specs_all = shape_residual_ptcent7_specs(args, args.outdir)
    elif args.campaign == "ppg12-sixpack":
        specs_all = ppg12_pp_sixpack_specs(args, args.outdir)
    else:
        specs_all = campaign_specs(args, args.outdir)
    specs = filter_specs(specs_all, args)
    split_columns = ["run", "evt"] if getattr(args, "split_mode", "row") == "event50" else []
    event_quality_enabled = args.event_quality_cut_json is not None
    event_quality_required = EVENT_QUALITY_DIRECT_COLUMNS + EVENT_QUALITY_COMPONENT_COLUMNS + ["run", "evt"] if event_quality_enabled else []
    all_features = sorted(
        set(
            expand_required_columns(
                [feature for spec in specs for feature in spec["features"]]
                + ["centrality", "cluster_Et", "cluster_Eta", label_branch]
                + split_columns
                + event_quality_required
            )
        )
    )
    if args.majority_cap_ratio <= 0.0 and args.campaign != "ppg12-sixpack" and args.weight_mode != "ppg12-exact":
        args.majority_cap_ratio = 4.0

    planned_path = args.registry_output or (args.outdir / "model_registry.planned.json")
    if args.plan_only:
        write_registry(planned_path, specs, [], args, "PLANNED")
        print(f"[OK] planned campaign specs={len(specs_all)} selected={len(specs)}")
        return 0

    args.outdir.mkdir(parents=True, exist_ok=True)
    if args.weight_mode == "ppg12-exact":
        args.majority_cap_ratio = 0.0
        if args.use_event_weight:
            print("[INFO] --weight-mode ppg12-exact ignores event_weight for BDT training.", flush=True)
    if args.parallel_workers > 1 and args.n_jobs > 1:
        args.n_jobs = 1

    optional_columns = [PPG12_EXACT_WEIGHT_COLUMN] if args.weight_mode == "ppg12-exact" else [args.weight_branch]
    frame, optional_seen, from_cache = load_or_build_frame(
        paths,
        args.tree,
        all_features,
        optional_columns,
        label_branch,
        args.missing_label_value,
        None if args.event_quality_audit_only else args.cache_file,
        args.cache_only,
        skip_missing_tree=args.skip_missing_tree,
        max_load_rows_per_class=int(args.max_load_rows_per_class or 0),
        max_load_rows=int(args.max_load_rows or 0),
        load_sample_seed=int(args.load_sample_seed or args.random_seed or 42),
    )
    event_quality_report = {"enabled": False}
    if event_quality_enabled:
        audit_output = args.event_quality_audit_output
        if audit_output is None and (args.campaign_spec_list is None or args.event_quality_audit_only):
            audit_output = args.outdir / "event_quality_filter_audit.json"
        if args.event_quality_assume_filtered:
            if args.event_quality_audit_only:
                raise SystemExit("--event-quality-assume-filtered cannot be combined with --event-quality-audit-only")
            frame, event_quality_report = assume_low_calo_event_quality_filter_applied(
                frame,
                args.event_quality_cut_json,
                audit_output,
                args.outdir,
            )
        else:
            frame, event_quality_report = apply_low_calo_event_quality_filter(
                frame,
                args.event_quality_cut_json,
                audit_output=audit_output,
                audit_only=bool(args.event_quality_audit_only),
            )
        if (
            event_quality_report.get("retained_below_envelope_events", 0) != 0
            or event_quality_report.get("retained_below_envelope_candidates", 0) != 0
        ):
            raise SystemExit(
                "Low-calo upstream filter failed closure target: retained below-envelope "
                f"events={event_quality_report.get('retained_below_envelope_events')} "
                f"candidates={event_quality_report.get('retained_below_envelope_candidates')}"
            )
        if args.cache_file is not None and not args.event_quality_audit_only:
            cache_optional_columns = [] if args.weight_mode == "ppg12-exact" else optional_columns
            cache_cols = sorted(
                set(
                    expand_required_columns(all_features)
                    + all_features
                    + cache_optional_columns
                    + [
                        "source_sample",
                        "input_file_index",
                        "input_tree_entry",
                    ]
                    + TRAINING_IDENTITY_OPTIONAL_COLUMNS
                )
            )
            save_frame_cache(frame, args.cache_file, cache_cols)
            print(f"[OK] wrote upstream-filtered training cache: {args.cache_file}", flush=True)
        if args.event_quality_audit_only:
            planned_path.parent.mkdir(parents=True, exist_ok=True)
            planned_path.write_text(
                json.dumps(
                    {
                        "schema": "AUAU_LOW_CALO_UPSTREAM_FILTER_DRYRUN_V1",
                        "status": "EVENT_QUALITY_FILTER_AUDIT_READY",
                        "campaign": args.campaign,
                        "selected_specs": len(specs),
                        "event_quality_filter": event_quality_report,
                    },
                    indent=2,
                    sort_keys=True,
                )
                + "\n"
            )
            print(f"[OK] dry-run event-quality audit passed: {audit_output}", flush=True)
            return 0
    if args.cache_only and args.weight_mode == "ppg12-exact" and args.cache_only_skip_ppg12_exact_weights:
        write_registry(planned_path, specs, [], args, "CACHE_READY_UNWEIGHTED")
        return 0
    ppg12_exact_closure = None
    if args.weight_mode == "ppg12-exact":
        precomputed_status = ppg12_exact_precomputed_weight_status(frame)
        if precomputed_status.get("usable"):
            ppg12_exact_closure = summarize_ppg12_exact_precomputed_weights(frame, label_branch, args)
        else:
            if precomputed_status.get("has_column"):
                print(
                    "[WARN] ignoring unusable PPG12-exact cached weight column: "
                    f"{precomputed_status.get('reason')}; recomputing before CACHE_READY",
                    flush=True,
                )
                frame = frame.drop(columns=[PPG12_EXACT_WEIGHT_COLUMN])
            frame, ppg12_exact_closure = prepare_ppg12_exact_global_weights(frame, label_branch, args)
            if args.cache_file is not None:
                cache_cols = sorted(
                    set(
                        expand_required_columns(all_features)
                        + all_features
                        + optional_columns
                        + [
                            "source_sample",
                            "input_file_index",
                            "input_tree_entry",
                            PPG12_EXACT_WEIGHT_COLUMN,
                        ]
                        + TRAINING_IDENTITY_OPTIONAL_COLUMNS
                        + split_columns
                    )
                )
                save_frame_cache(frame, args.cache_file, cache_cols)
                print(f"[OK] updated training cache with PPG12-exact weights: {args.cache_file}", flush=True)

    if args.cache_only:
        write_registry(planned_path, specs, [], args, "CACHE_READY")
        return 0

    common_metadata = {
        "task": args.task,
        "input_files": [str(path) for path in paths],
        "tree": args.tree,
        "optional_branches_seen": optional_seen,
        "loaded_from_cache": from_cache,
        "cache_file": str(args.cache_file) if args.cache_file else None,
        "python": sys.version,
        "label_branch": label_branch,
        "weight_mode": args.weight_mode,
        "ppg12_exact_closure": ppg12_exact_closure,
        "event_quality_filter": event_quality_report,
        "global_event_key_audit": summarize_global_event_key(frame),
    }

    reports: list[dict] = []
    workers = max(1, int(args.parallel_workers))
    if workers == 1 or len(specs) <= 1:
        for spec in specs:
            report = train_campaign_spec(frame, spec, common_metadata, args)
            if report is not None:
                reports.append(report)
    else:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            future_map = {pool.submit(train_campaign_spec, frame, spec, common_metadata, args): spec for spec in specs}
            for future in as_completed(future_map):
                spec = future_map[future]
                try:
                    report = future.result()
                except Exception as exc:  # noqa: BLE001
                    raise SystemExit(f"Campaign spec {spec['model_id']} failed: {exc}") from exc
                if report is not None:
                    reports.append(report)

    registry_path = args.registry_output or (args.outdir / "model_registry.json")
    write_registry(registry_path, specs, reports, args, "READY")
    summary = args.outdir / "expanded_campaign_summary.json"
    summary.write_text(json.dumps(reports, indent=2, sort_keys=True) + "\n")
    print(f"[OK] expanded campaign trained reports={len(reports)} selected_specs={len(specs)} registry={registry_path}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task", choices=["tight", "npb"], required=True)
    parser.add_argument("--input", nargs="+", type=Path, default=None, help="ROOT files or @manifest.list files")
    parser.add_argument("--tree", default="AuAuPhotonIDTrainingTree")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--prefix", default="")
    parser.add_argument("--cent-bins", default="0:10,10:20,20:40,40:60,60:80,80:100")
    parser.add_argument("--tight-mode", choices=TIGHT_MODES, default="legacy",
                        help="Tight-BDT training product. centAsFeat3x3 is the centrality-input model with 3x3 shower-width moments.")
    parser.add_argument("--label-branch", default=None)
    parser.add_argument(
        "--missing-label-value",
        type=int,
        choices=[0, 1],
        default=None,
        help="For files missing only the label branch, fill that label with this value. For PPG12-style NPB, use 1 for embedded physics-side sim.",
    )
    parser.add_argument(
        "--skip-missing-tree",
        action="store_true",
        help="Skip input ROOT files that do not contain the requested training tree instead of failing immediately.",
    )
    parser.add_argument("--features", default=None, help="Comma-separated override feature order")
    parser.add_argument("--weight-branch", default="event_weight")
    parser.add_argument(
        "--weight-mode",
        choices=["legacy", "ppg12-exact"],
        default="legacy",
        help="Training-weight convention. ppg12-exact computes global class/ET/eta weights before routed model slicing and ignores event_weight.",
    )
    parser.add_argument(
        "--ppg12-exact-expected-samples",
        default=",".join(PPG12_EXACT_EXPECTED_SAMPLES),
        help="Comma-separated required source_sample names for --weight-mode ppg12-exact.",
    )
    parser.add_argument(
        "--ppg12-exact-closure-dir",
        type=Path,
        default=None,
        help="Output directory for PPG12-exact sample-mix and ET/eta closure PNG/CSV/JSON artifacts.",
    )
    parser.add_argument(
        "--ppg12-exact-closure-artifacts",
        choices=["full", "metadata-only"],
        default="full",
        help="Control PPG12-exact closure artifacts. metadata-only preserves training weights while deferring heavy CSV/PNG diagnostics.",
    )
    parser.add_argument("--use-event-weight", dest="use_event_weight", action="store_true", default=True)
    parser.add_argument("--no-event-weight", dest="use_event_weight", action="store_false")
    parser.add_argument("--et-reweight", dest="et_reweight", action="store_true", default=True)
    parser.add_argument("--no-et-reweight", dest="et_reweight", action="store_false")
    parser.add_argument("--eta-reweight", dest="eta_reweight", action="store_true", default=True)
    parser.add_argument("--no-eta-reweight", dest="eta_reweight", action="store_false")
    parser.add_argument("--flatten-bins", type=int, default=20)
    parser.add_argument("--max-flatten-factor", type=float, default=8.0)
    parser.add_argument("--max-total-weight-factor", type=float, default=50.0)
    parser.add_argument("--min-rows-per-class", type=int, default=10)
    parser.add_argument("--test-size", type=float, default=0.10)
    parser.add_argument(
        "--split-mode",
        choices=["row", "event50"],
        default="row",
        help="Train/test split convention. event50 requires run/evt and splits whole events deterministically.",
    )
    parser.add_argument(
        "--require-global-event-key",
        action="store_true",
        default=os.environ.get("RJ_AUAU_BDT_REQUIRE_GLOBAL_EVENT_KEY", "0") == "1",
        help="Require source_sample + stable input_file_index + run + evt event keys where event-level deduplication is audited.",
    )
    parser.add_argument(
        "--event-quality-cut-json",
        type=Path,
        default=Path(os.environ["RJ_AUAU_BDT_EVENT_QUALITY_CUT_JSON"])
        if os.environ.get("RJ_AUAU_BDT_EVENT_QUALITY_CUT_JSON")
        else None,
        help="Explicitly enable the source/truth/score-blind event-quality filter using this stored low-calo cut JSON.",
    )
    parser.add_argument(
        "--event-quality-audit-output",
        type=Path,
        default=Path(os.environ["RJ_AUAU_BDT_EVENT_QUALITY_AUDIT_OUTPUT"])
        if os.environ.get("RJ_AUAU_BDT_EVENT_QUALITY_AUDIT_OUTPUT")
        else None,
        help="Write the upstream event-quality filter audit JSON to this path.",
    )
    parser.add_argument(
        "--event-quality-audit-only",
        action="store_true",
        default=os.environ.get("RJ_AUAU_BDT_EVENT_QUALITY_AUDIT_ONLY", "0") == "1",
        help="Run the upstream low-calo filter audit and exit before training.",
    )
    parser.add_argument(
        "--event-quality-assume-filtered",
        action="store_true",
        default=os.environ.get("RJ_AUAU_BDT_EVENT_QUALITY_ASSUME_FILTERED", "0") == "1",
        help="Resume from a cache that already has the upstream low-calo event-quality filter applied; requires a matching audit JSON.",
    )
    parser.add_argument("--random-seed", type=int, default=13)
    parser.add_argument("--n-estimators", type=int, default=450)
    parser.add_argument("--max-depth", type=int, default=4)
    parser.add_argument("--learning-rate", type=float, default=0.035)
    parser.add_argument("--subsample", type=float, default=0.85)
    parser.add_argument("--colsample-bytree", type=float, default=0.85)
    parser.add_argument("--tree-method", default="hist")
    parser.add_argument("--reg-alpha", type=float, default=5.0)
    parser.add_argument("--reg-lambda", type=float, default=0.3)
    parser.add_argument("--grow-policy", default="lossguide")
    parser.add_argument("--max-bin", type=int, default=256)
    parser.add_argument("--n-jobs", type=int, default=4)
    parser.add_argument(
        "--campaign",
        choices=[
            "expanded-tight",
            "etfine-centstudy",
            "iso-diagnostic",
            "global-sixpack",
            "basev3e-e22-ablation",
            "corrected-baseline-shower-ladder",
            "corrected-baseline-binned14",
            "corrected-baseline-iso14",
            "etcent-binned-sixpack",
            "etcent-binned-sixpack-noiso-ptcent7",
            "global-and-etcent-binned-sixpack-noiso",
            "etcent-binned-eiso-cone-ablation",
            "shape-residual-ptcent7",
            "ppg12-sixpack",
        ],
        default=None,
    )
    parser.add_argument("--plan-only", action="store_true")
    parser.add_argument("--cache-only", action="store_true")
    parser.add_argument(
        "--cache-only-skip-ppg12-exact-weights",
        action="store_true",
        default=os.environ.get("RJ_AUAU_BDT_CACHE_ONLY_SKIP_PPG12_EXACT_WEIGHTS", "0") == "1",
        help=(
            "For staged Condor cache builds only: write the post-selection frame cache and "
            "exit before computing global PPG12-exact weights. A later reduce stage must "
            "merge all shards and compute the global weights once before training."
        ),
    )
    parser.add_argument("--cache-file", type=Path, default=None)
    parser.add_argument("--registry-output", type=Path, default=None)
    parser.add_argument("--campaign-spec-list", type=Path, default=None)
    parser.add_argument("--campaign-spec-ids", default="")
    parser.add_argument("--campaign-max-specs", type=int, default=0)
    parser.add_argument("--pt-bins", default="5,8,10,12,14,16,18,20,22,24,26,35")
    parser.add_argument("--coarse-cent-bins", default="0:20,20:50,50:80")
    parser.add_argument("--fine-cent-bins", default="0:10,10:20,20:30,30:40,40:50,50:60,60:80")
    parser.add_argument("--extra-cent-as-feat-pt-ranges", default="")
    parser.add_argument("--extra-cent-as-feat-3x3-pt-ranges", default="")
    parser.add_argument("--extra-cent-as-feat-base3x3-pt-ranges", default="")
    parser.add_argument("--extra-cent-as-feat-width-ratio-pt-ranges", default="")
    parser.add_argument("--parallel-workers", type=int, default=int(os.environ.get("RJ_AUAU_BDT_TRAIN_PARALLEL", "4")))
    parser.add_argument("--majority-cap-ratio", type=float, default=0.0)
    parser.add_argument("--minopt-majority-cap-ratio", type=float, default=float(os.environ.get("RJ_AUAU_BDT_MINOPT_MAJORITY_CAP_RATIO", "2.0")))
    parser.add_argument("--majority-cap-et-bins", type=int, default=12)
    parser.add_argument("--majority-cap-cent-bins", type=int, default=8)
    parser.add_argument("--background-subsample-fraction", type=float, default=1.0)
    parser.add_argument("--background-subsample-et-threshold", type=float, default=15.0)
    parser.add_argument("--background-subsample-bins", type=int, default=20)
    parser.add_argument("--background-subsample-seed", type=int, default=42)
    parser.add_argument("--background-subsample-flatten", action="store_true")
    parser.add_argument("--max-load-rows", type=int, default=int(os.environ.get("RJ_AUAU_BDT_MAX_LOAD_ROWS", "0")),
                        help="Deterministically sample at most this many rows while reading ROOT files, before concatenating into memory.")
    parser.add_argument("--max-load-rows-per-class", type=int, default=int(os.environ.get("RJ_AUAU_BDT_MAX_LOAD_ROWS_PER_CLASS", "0")),
                        help="Deterministically sample at most this many rows for each binary label while reading ROOT files.")
    parser.add_argument("--load-sample-seed", type=int, default=int(os.environ.get("RJ_AUAU_BDT_LOAD_SAMPLE_SEED", "42")),
                        help="Seed used for load-time file shuffling and row sampling when max-load caps are enabled.")
    parser.add_argument("--allow-tmva-export-failure", action="store_true",
                        help="Write diagnostics but do not fail if TMVA export fails. Not recommended for production pipeline tests.")
    parser.add_argument("--skip-tmva-export", action="store_true",
                        help="Write XGBoost JSON/metadata only and do not import ROOT/TMVA. Intended for validation-only campaigns.")
    args = parser.parse_args()

    if args.campaign:
        return run_campaign(args)

    if not args.input:
        raise SystemExit("--input is required unless --campaign --plan-only is used")

    feature_override = [item.strip() for item in args.features.split(",")] if args.features else None
    features = tight_mode_features(args.tight_mode, feature_override) if args.task == "tight" else (
        list(feature_override) if feature_override is not None else list(PPG12_NPB_FEATURES)
    )
    label_branch = args.label_branch or ("is_signal" if args.task == "tight" else "npb_label")
    if args.task == "npb" and label_branch == "is_signal":
        raise SystemExit("NPB training needs a real NPB label branch, not is_signal.")
    if args.task == "npb" and label_branch == "is_npb":
        raise SystemExit(
            "Use npb_label for PPG12-style NPB training: 1=physics-like, 0=NPB. "
            "The is_npb branch is kept only as an audit inverse."
        )

    paths = expand_input_paths(args.input)
    split_columns = ["run", "evt"] if args.split_mode == "event50" else []
    event_quality_enabled = args.event_quality_cut_json is not None
    event_quality_required = EVENT_QUALITY_DIRECT_COLUMNS + EVENT_QUALITY_COMPONENT_COLUMNS + ["run", "evt"] if event_quality_enabled else []
    required_columns = sorted(set(features + [label_branch, "centrality"] + split_columns + event_quality_required))
    optional_columns = [PPG12_EXACT_WEIGHT_COLUMN] if args.weight_mode == "ppg12-exact" else [args.weight_branch]
    frame, optional_seen = load_frame(
        paths,
        args.tree,
        required_columns,
        optional_columns,
        label_branch,
        args.missing_label_value,
        skip_missing_tree=args.skip_missing_tree,
        label_branch=label_branch,
        max_load_rows_per_class=int(args.max_load_rows_per_class or 0),
        max_load_rows=int(args.max_load_rows or 0),
        load_sample_seed=int(args.load_sample_seed or args.random_seed),
    )

    args.outdir.mkdir(parents=True, exist_ok=True)
    event_quality_report = {"enabled": False}
    if event_quality_enabled:
        audit_output = args.event_quality_audit_output or (args.outdir / "event_quality_filter_audit.json")
        if args.event_quality_assume_filtered:
            if args.event_quality_audit_only:
                raise SystemExit("--event-quality-assume-filtered cannot be combined with --event-quality-audit-only")
            frame, event_quality_report = assume_low_calo_event_quality_filter_applied(
                frame,
                args.event_quality_cut_json,
                audit_output,
                args.outdir,
            )
        else:
            frame, event_quality_report = apply_low_calo_event_quality_filter(
                frame,
                args.event_quality_cut_json,
                audit_output=audit_output,
                audit_only=bool(args.event_quality_audit_only),
            )
        if (
            event_quality_report.get("retained_below_envelope_events", 0) != 0
            or event_quality_report.get("retained_below_envelope_candidates", 0) != 0
        ):
            raise SystemExit(
                "Low-calo upstream filter failed closure target: retained below-envelope "
                f"events={event_quality_report.get('retained_below_envelope_events')} "
                f"candidates={event_quality_report.get('retained_below_envelope_candidates')}"
            )
        if args.event_quality_audit_only:
            print(f"[OK] dry-run event-quality audit passed: {audit_output}", flush=True)
            return 0
    ppg12_exact_closure = None
    if args.weight_mode == "ppg12-exact":
        args.majority_cap_ratio = 0.0
        if args.use_event_weight:
            print("[INFO] --weight-mode ppg12-exact ignores event_weight for BDT training.", flush=True)
        frame, ppg12_exact_closure = prepare_ppg12_exact_global_weights(frame, label_branch, args)
    prefix = args.prefix or f"auau_{args.task}_bdt"
    reports = []
    common_metadata = {
        "task": args.task,
        "input_files": [str(path) for path in paths],
        "tree": args.tree,
        "optional_branches_seen": optional_seen,
        "python": sys.version,
        "tight_mode": args.tight_mode if args.task == "tight" else None,
        "weight_mode": args.weight_mode,
        "ppg12_exact_closure": ppg12_exact_closure,
        "event_quality_filter": event_quality_report,
        "global_event_key_audit": summarize_global_event_key(frame),
    }

    train_all_cent = args.task != "tight" or args.tight_mode in ("legacy", "ppg12BaseV1E", "centINDcontrol", "centAsFeat", "centAsFeatMinOpt", "centAsFeat3x3", "centAsFeatBase3x3", "centAsFeatWidthRatios")
    train_cent_bins = args.task != "tight" or args.tight_mode in ("legacy", "centDepBDTs")

    if train_all_cent:
        all_report = train_one(
            frame,
            features,
            label_branch,
            args.outdir / f"{prefix}_allCent_tmva.root",
            {**common_metadata, "cent_range": "all"},
            args,
        )
        if all_report is not None:
            reports.append(all_report)

    if train_cent_bins:
        for lo, hi in parse_cent_bins(args.cent_bins):
            sub = frame[(frame["centrality"] >= lo) & (frame["centrality"] < hi)]
            if len(sub) == 0:
                print(f"[WARN] Skipping {cent_tag(lo, hi)}: no rows")
                continue
            report = train_one(
                sub,
                features,
                label_branch,
                args.outdir / f"{prefix}_{cent_tag(lo, hi)}_tmva.root",
                {**common_metadata, "cent_range": [lo, hi]},
                args,
            )
            if report is not None:
                reports.append(report)

    summary = args.outdir / f"{prefix}_summary.json"
    summary.write_text(json.dumps(reports, indent=2, sort_keys=True) + "\n")
    print(f"[OK] wrote {summary}")
    for report in reports:
        print(
            f"{Path(report['output_tmva']).name}: rows={report['n_rows']} "
            f"sig={report['n_signal']} bkg={report['n_background']} "
            f"auc={report['auc']:.4f} tmva={report['tmva_export']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
