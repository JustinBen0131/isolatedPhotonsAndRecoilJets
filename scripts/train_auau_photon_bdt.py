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

DERIVED_FEATURE_DEPS = {
    "cluster_weta_over_wphi": ["cluster_weta_cogx", "cluster_wphi_cogx"],
    "cluster_weta33_over_wphi33": ["cluster_weta33_cogx", "cluster_wphi33_cogx"],
    "reco_eiso_clip30": ["reco_eiso"],
    "reco_eiso_over_cluster_Et": ["reco_eiso", "cluster_Et"],
    "reco_eiso_signed_log1p": ["reco_eiso"],
}

TIGHT_MODES = [
    "legacy",
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


def tight_mode_features(mode: str, override: list[str] | None) -> list[str]:
    if override is not None:
        features = list(override)
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

    def safe_ratio(num: str, den: str):
        n = frame[num].to_numpy(dtype="float64")
        d = frame[den].to_numpy(dtype="float64")
        out = np.full(len(frame), np.nan, dtype="float64")
        good = np.isfinite(n) & np.isfinite(d) & (np.abs(d) > 1.0e-9)
        out[good] = n[good] / d[good]
        return out

    if "cluster_weta_over_wphi" not in frame.columns and {"cluster_weta_cogx", "cluster_wphi_cogx"}.issubset(frame.columns):
        frame["cluster_weta_over_wphi"] = safe_ratio("cluster_weta_cogx", "cluster_wphi_cogx")
    if "cluster_weta33_over_wphi33" not in frame.columns and {"cluster_weta33_cogx", "cluster_wphi33_cogx"}.issubset(frame.columns):
        frame["cluster_weta33_over_wphi33"] = safe_ratio("cluster_weta33_cogx", "cluster_wphi33_cogx")
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
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        preview = "\n  ".join(missing[:20])
        extra = "" if len(missing) <= 20 else f"\n  ... {len(missing) - 20} more"
        raise SystemExit(f"Input ROOT files are missing:\n  {preview}{extra}")
    return paths


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
):
    try:
        import pandas as pd
        import uproot
    except ImportError as exc:
        raise SystemExit(
            "Missing dependency. Use an environment with uproot, pandas, numpy, "
            "scikit-learn, xgboost, and PyROOT for TMVA export."
        ) from exc

    frames = []
    seen_optional: set[str] = set()
    for path in paths:
        with uproot.open(path) as root_file:
            try:
                tree = root_file[tree_name]
            except Exception as exc:
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
            seen_optional.update(present_optional)
            read_columns = [col for col in required_columns if col in keys] + present_optional
            frame = tree.arrays(read_columns, library="pd")
            if allow_missing_label:
                frame[missing_label_branch] = int(missing_label_value)
            frames.append(frame)

    frame = add_derived_features(pd.concat(frames, ignore_index=True))
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
):
    if cache_file is not None and cache_file.is_file():
        frame = add_derived_features(load_frame_cache(cache_file))
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
    )
    if cache_file is not None:
        cache_cols = sorted(set(expand_required_columns(required_columns) + required_columns + optional_columns))
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


def compute_weights(frame, label_branch: str, args) -> tuple[object, dict]:
    import numpy as np

    labels = frame[label_branch].to_numpy(dtype="int32")
    weights = np.ones(len(frame), dtype="float64")
    diagnostics: dict[str, object] = {
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


def train_one(frame, features: list[str], label_branch: str, output: Path, metadata: dict, args) -> dict | None:
    import numpy as np
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import train_test_split
    from xgboost import XGBClassifier

    cols = features + [label_branch]
    frame = frame.loc[finite_mask(frame, cols)].copy()
    frame[label_branch] = frame[label_branch].astype(int)
    frame = frame[frame[label_branch].isin([0, 1])].copy()
    n_rows_after_finite = int(len(frame))
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

    stratify = y if min(n_sig, n_bkg) >= 2 else None
    x_train, x_test, y_train, y_test, w_train, w_test = train_test_split(
        x, y, weights, test_size=args.test_size, random_state=args.random_seed, stratify=stratify
    )

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
        eval_metric="auc",
        tree_method=args.tree_method,
        random_state=args.random_seed,
    )
    model.fit(x_train, y_train, sample_weight=w_train)

    pred = model.predict_proba(x_test)[:, 1]
    auc = float(roc_auc_score(y_test, pred, sample_weight=w_test)) if len(np.unique(y_test)) == 2 else math.nan

    output.parent.mkdir(parents=True, exist_ok=True)
    json_model = output.with_suffix(".xgb.json")
    model.get_booster().save_model(json_model)

    export_status = "not_attempted"
    export_error = ""
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
        "tmva_export": export_status,
        "weighting": weight_report,
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


def registry_payload(specs: list[dict], reports: list[dict], args, status: str = "PLANNED") -> dict:
    report_by_id = {r.get("model_id"): r for r in reports}
    products: dict[str, list[str]] = {}
    for spec in specs:
        products.setdefault(spec["product"], []).append(spec["model_id"])
    return {
        "schema": "AUAU_TIGHT_BDT_EXPANDED_REGISTRY_V1",
        "status": status,
        "campaign": args.campaign,
        "expected_model_count": len(specs),
        "model_count": len(specs),
        "products": products,
        "pt_bins": parse_float_edges(args.pt_bins),
        "coarse_cent_bins": [[lo, hi] for lo, hi in parse_cent_bins(args.coarse_cent_bins)],
        "fine_cent_bins": [[lo, hi] for lo, hi in parse_cent_bins(args.fine_cent_bins)],
        "defaults": {
            "majority_cap_ratio": float(args.majority_cap_ratio),
            "minopt_majority_cap_ratio": float(args.minopt_majority_cap_ratio),
            "parallel_workers": int(args.parallel_workers),
            "xgboost_n_jobs": int(args.n_jobs),
            "no_cross_section_weights": True,
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
    else:
        specs_all = campaign_specs(args, args.outdir)
    specs = filter_specs(specs_all, args)
    all_features = sorted(set(expand_required_columns([feature for spec in specs for feature in spec["features"]] + ["centrality", label_branch])))
    if args.majority_cap_ratio <= 0.0:
        args.majority_cap_ratio = 4.0

    planned_path = args.registry_output or (args.outdir / "model_registry.planned.json")
    if args.plan_only:
        write_registry(planned_path, specs, [], args, "PLANNED")
        print(f"[OK] planned campaign specs={len(specs_all)} selected={len(specs)}")
        return 0

    args.outdir.mkdir(parents=True, exist_ok=True)
    if args.parallel_workers > 1 and args.n_jobs > 1:
        args.n_jobs = 1

    frame, optional_seen, from_cache = load_or_build_frame(
        paths,
        args.tree,
        all_features,
        optional_columns,
        label_branch,
        args.missing_label_value,
        args.cache_file,
        args.cache_only,
    )
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
    parser.add_argument("--features", default=None, help="Comma-separated override feature order")
    parser.add_argument("--weight-branch", default="event_weight")
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
    parser.add_argument("--campaign", choices=["expanded-tight", "etfine-centstudy", "iso-diagnostic", "global-sixpack"], default=None)
    parser.add_argument("--plan-only", action="store_true")
    parser.add_argument("--cache-only", action="store_true")
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
    parser.add_argument("--allow-tmva-export-failure", action="store_true",
                        help="Write diagnostics but do not fail if TMVA export fails. Not recommended for production pipeline tests.")
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
    required_columns = sorted(set(features + [label_branch, "centrality"]))
    optional_columns = [args.weight_branch]
    frame, optional_seen = load_frame(
        paths,
        args.tree,
        required_columns,
        optional_columns,
        label_branch,
        args.missing_label_value,
    )

    args.outdir.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix or f"auau_{args.task}_bdt"
    reports = []
    common_metadata = {
        "task": args.task,
        "input_files": [str(path) for path in paths],
        "tree": args.tree,
        "optional_branches_seen": optional_seen,
        "python": sys.version,
        "tight_mode": args.tight_mode if args.task == "tight" else None,
    }

    train_all_cent = args.task != "tight" or args.tight_mode in ("legacy", "centINDcontrol", "centAsFeat", "centAsFeatMinOpt", "centAsFeat3x3", "centAsFeatBase3x3", "centAsFeatWidthRatios")
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
