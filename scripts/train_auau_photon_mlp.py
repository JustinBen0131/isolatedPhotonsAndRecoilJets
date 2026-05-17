#!/usr/bin/env python3
"""Train compact AuAu photon-ID MLP artifacts from RecoilJets training trees.

This is intentionally a sidecar to ``train_auau_photon_bdt.py``.  It uses the
same ``AuAuPhotonIDTrainingTree`` branches and weighting conventions, but it
exports a small deterministic JSON artifact that production C++ can evaluate
without PyTorch, TensorFlow, ONNX, or scikit-learn.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import os
import random
import sys
from pathlib import Path
from typing import Iterable


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

BASE_AND_3X3_FEATURES = [
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

MODEL_SPECS = {
    "globalEtCent1535_mlp_noIso": {
        "variant": "auauGlobalEtCent1535MLP",
        "label": "Global 15-35 GeV MLP with E_T, centrality, full shower-shape inputs",
        "features": BASE_AND_3X3_FEATURES + EXTENDED_SHOWER_FEATURES + WIDTH_RATIO_FEATURES + ["centrality"],
        "filename": "auau_globalEtCent1535_mlp_noIso.json",
    },
    "globalEtCent1535_mlp_iso": {
        "variant": "auauGlobalEtCent1535IsoDiagnosticMLP",
        "label": "Global 15-35 GeV MLP with E_T, centrality, full shower-shape and isolation inputs",
        "features": BASE_AND_3X3_FEATURES
        + EXTENDED_SHOWER_FEATURES
        + WIDTH_RATIO_FEATURES
        + ISOLATION_DIAGNOSTIC_FEATURES
        + ["centrality"],
        "filename": "auau_globalEtCent1535_mlp_iso.json",
        "diagnostic_only": True,
        "abcd_warning": "Uses reconstructed isolation-derived inputs; diagnostic only.",
    },
    "centInputBase3x3MLP_pt1535": {
        "variant": "auauCentInputBase3x3MLP",
        "label": "Centrality input, base + 3x3 widths MLP",
        "features": BASE_AND_3X3_FEATURES + ["centrality"],
        "filename": "auau_tight_mlp_centInputBase3x3_pt1535.json",
    },
    "centInputBase3x3WidthRatiosMLP_pt1535": {
        "variant": "auauCentInputBase3x3MLP",
        "label": "Centrality input, base + 3x3 widths + width ratios MLP",
        "features": BASE_AND_3X3_FEATURES + WIDTH_RATIO_FEATURES + ["centrality"],
        "filename": "auau_tight_mlp_centInputBase3x3_pt1535.json",
    },
    "centInputKitchenSinkMLP": {
        "variant": "auauCentInputKitchenSinkMLP",
        "label": "Centrality input, extended shower-shape and energy-ratio kitchen-sink MLP",
        "features": BASE_AND_3X3_FEATURES + EXTENDED_SHOWER_FEATURES + WIDTH_RATIO_FEATURES + ["centrality"],
        "filename": "auau_tight_mlp_centInputKitchenSink.json",
    },
    "highPtDistilledKitchenMLP_v2": {
        "variant": "auauHighPtDistilledKitchenMLP",
        "label": "ABCD-safe high-pT distilled kitchen-sink MLP v2",
        "features": BASE_AND_3X3_FEATURES + EXTENDED_SHOWER_FEATURES + WIDTH_RATIO_FEATURES + ["centrality"],
        "filename": "auau_tight_mlp_highPtDistilledKitchen_v2.json",
        "abcd_safe": True,
        "training_only_teacher": "auau_tight_bdt_score",
    },
    "centInputKitchenSinkIsoMLP": {
        "variant": "auauCentInputKitchenSinkIsoDiagnosticMLP",
        "label": "Diagnostic only: kitchen-sink MLP with reconstructed isolation inputs",
        "features": BASE_AND_3X3_FEATURES
        + EXTENDED_SHOWER_FEATURES
        + WIDTH_RATIO_FEATURES
        + ISOLATION_DIAGNOSTIC_FEATURES
        + ["centrality"],
        "filename": "auau_tight_mlp_centInputKitchenSinkIsoDiagnostic.json",
        "diagnostic_only": True,
        "abcd_warning": "Uses reconstructed isolation-derived inputs; do not promote into ABCD purity production without redesigning the purity method.",
    },
    "centInputMLP_pt1535": {
        "variant": "auauCentInputMLP",
        "label": "Centrality input, base widths MLP",
        "features": BASE_FEATURES + ["centrality"],
        "filename": "auau_tight_mlp_centInput_pt1535.json",
    },
    "noCentBase3x3MLP_pt1535": {
        "variant": "auauNoCentBase3x3MLP",
        "label": "No centrality input, base + 3x3 widths MLP",
        "features": BASE_AND_3X3_FEATURES,
        "filename": "auau_tight_mlp_noCentBase3x3_pt1535.json",
    },
}

DEFAULT_PRODUCTS = [
    "centInputBase3x3MLP_pt1535",
    "centInputMLP_pt1535",
    "noCentBase3x3MLP_pt1535",
]


def stable_seed(*items: object) -> int:
    text = "|".join(str(item) for item in items)
    digest = hashlib.sha256(text.encode("utf-8")).hexdigest()
    return int(digest[:8], 16)


def stable_unit_interval(*items: object) -> float:
    return stable_seed(*items) / float(0xFFFFFFFF)


def json_ready(obj):
    try:
        import numpy as np
    except Exception:
        np = None
    if np is not None:
        if isinstance(obj, np.generic):
            return obj.item()
        if isinstance(obj, np.ndarray):
            return [json_ready(x) for x in obj.tolist()]
    if isinstance(obj, dict):
        return {str(k): json_ready(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [json_ready(v) for v in obj]
    if isinstance(obj, float) and not math.isfinite(obj):
        return None
    return obj


def parse_range(text: str) -> tuple[float, float]:
    if ":" not in text:
        raise SystemExit(f"Range must be lo:hi, got {text!r}")
    lo_s, hi_s = text.split(":", 1)
    lo = float(lo_s)
    hi = float(hi_s)
    if hi <= lo:
        raise SystemExit(f"Range upper edge must exceed lower edge: {text!r}")
    return lo, hi


def parse_bin_spec(text: str | None) -> list[tuple[float, float]]:
    """Parse ``15,20,25`` edges or ``15:20,20:25`` bin specs."""
    if text is None or not str(text).strip():
        return []
    cleaned = str(text).replace(";", ",")
    parts = [item.strip() for item in cleaned.split(",") if item.strip()]
    if not parts:
        return []
    bins: list[tuple[float, float]] = []
    if all(":" in item for item in parts):
        for item in parts:
            bins.append(parse_range(item))
    else:
        edges = [float(item) for item in parts]
        if len(edges) < 2:
            raise SystemExit(f"Bin spec needs at least two edges: {text!r}")
        for lo, hi in zip(edges[:-1], edges[1:]):
            if hi <= lo:
                raise SystemExit(f"Bin edges must be strictly increasing: {text!r}")
            bins.append((float(lo), float(hi)))
    return bins


def bin_label(lo: float, hi: float) -> str:
    return f"{lo:g}-{hi:g}"


def train_pt_bins(args) -> list[tuple[float, float]]:
    bins = parse_bin_spec(getattr(args, "train_pt_bins", ""))
    if bins:
        return bins
    if getattr(args, "pt_range", ""):
        return [parse_range(args.pt_range)]
    return []


def parse_bin_weights(text: str | None, bins: list[tuple[float, float]], mode: str = "equal") -> list[float]:
    if not bins:
        return []
    if text and str(text).strip():
        parts = [item.strip() for item in str(text).replace(";", ",").split(",") if item.strip()]
        weights = [math.nan] * len(bins)
        if all(":" not in item for item in parts):
            raw = [float(item) for item in parts]
            if len(raw) != len(bins):
                raise SystemExit(f"Expected {len(bins)} bin weights, got {len(raw)} from {text!r}")
            weights = raw
        else:
            for item in parts:
                pieces = item.split(":")
                if len(pieces) != 3:
                    raise SystemExit(f"Bin weights must be lo:hi:weight, got {item!r}")
                lo, hi, weight = float(pieces[0]), float(pieces[1]), float(pieces[2])
                matched = False
                for idx, (blo, bhi) in enumerate(bins):
                    if abs(lo - blo) < 1.0e-6 and abs(hi - bhi) < 1.0e-6:
                        weights[idx] = weight
                        matched = True
                        break
                if not matched:
                    raise SystemExit(f"Weight bin {lo:g}:{hi:g} is not present in bin spec {bins}")
        if any(not math.isfinite(float(w)) or float(w) <= 0.0 for w in weights):
            raise SystemExit(f"Bin weights must be finite and positive: {text!r}")
        return [float(w) for w in weights]
    if mode == "highpt":
        centers = [0.5 * (lo + hi) for lo, hi in bins]
        lo_c = min(centers)
        span = max(max(centers) - lo_c, 1.0e-9)
        return [0.25 + 0.75 * ((c - lo_c) / span) ** 1.5 for c in centers]
    return [1.0 for _ in bins]


def parse_products(text: str) -> list[str]:
    if text.strip().lower() in ("all", "default"):
        return list(DEFAULT_PRODUCTS)
    if text.strip().lower() in ("primary", "primary-only"):
        return ["centInputBase3x3MLP_pt1535"]
    if text.strip().lower() in ("global-sixpack", "global-etcent-sixpack"):
        return ["globalEtCent1535_mlp_noIso", "globalEtCent1535_mlp_iso"]
    if text.strip().lower() in ("primary-ratios", "primary-width-ratios", "ratios", "width-ratios"):
        return ["centInputBase3x3WidthRatiosMLP_pt1535"]
    if text.strip().lower() in ("kitchen-sink", "kitchensink", "extended", "extended-shower"):
        return ["centInputKitchenSinkMLP"]
    if text.strip().lower() in (
        "highpt-distilled-kitchen-v2",
        "highpt-distilled-kitchen",
        "distilled-kitchen-v2",
        "bdt-distilled-kitchen",
        "v2",
    ):
        return ["highPtDistilledKitchenMLP_v2"]
    if text.strip().lower() in (
        "iso-kitchen-sink",
        "kitchen-sink-iso",
        "kitchensink-iso",
        "iso-aware",
        "iso-diagnostic",
        "oracle-iso",
    ):
        return ["centInputKitchenSinkIsoMLP"]
    products = [item.strip() for item in text.split(",") if item.strip()]
    unknown = [item for item in products if item not in MODEL_SPECS]
    if unknown:
        raise SystemExit(
            "Unknown MLP product(s): "
            + ", ".join(unknown)
            + "\nKnown: "
            + ", ".join(MODEL_SPECS)
        )
    return products


def expand_required_columns(features: Iterable[str]) -> list[str]:
    required: set[str] = set()
    for feature in features:
        deps = DERIVED_FEATURE_DEPS.get(feature)
        if deps:
            required.update(deps)
        else:
            required.add(feature)
    return sorted(required)


def frame_columns(frame) -> set[str]:
    if hasattr(frame, "columns"):
        return set(frame.columns)
    return set(frame.keys())


def add_derived_features(frame):
    import numpy as np

    cols = frame_columns(frame)

    def safe_ratio(num: str, den: str):
        numerator = np.asarray(frame[num], dtype="float64")
        denominator = np.asarray(frame[den], dtype="float64")
        out = np.full(len(numerator), np.nan, dtype="float64")
        mask = np.isfinite(numerator) & np.isfinite(denominator) & (np.abs(denominator) > 1.0e-12)
        out[mask] = numerator[mask] / denominator[mask]
        return out

    if "cluster_weta_over_wphi" not in cols and {"cluster_weta_cogx", "cluster_wphi_cogx"}.issubset(cols):
        frame["cluster_weta_over_wphi"] = safe_ratio("cluster_weta_cogx", "cluster_wphi_cogx").astype("float32")
        cols.add("cluster_weta_over_wphi")
    if "cluster_weta33_over_wphi33" not in cols and {"cluster_weta33_cogx", "cluster_wphi33_cogx"}.issubset(cols):
        frame["cluster_weta33_over_wphi33"] = safe_ratio("cluster_weta33_cogx", "cluster_wphi33_cogx").astype("float32")
        cols.add("cluster_weta33_over_wphi33")
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
        cols.add("reco_eiso_signed_log1p")
    return frame


def compact_numeric_frame(frame, label_branch: str = "is_signal"):
    """Keep ROOT-derived training frames compact after filtering.

    The full sixpack sample has O(10M) rows, so pandas' default float64
    columns are expensive. The training math promotes batch operations as
    needed, but storing the selected frame as float32/int32 avoids wasting
    memory before minibatch training starts.
    """
    try:
        import pandas as pd
    except Exception:
        pd = None
    for col in list(frame.columns):
        if col == label_branch:
            frame[col] = frame[col].astype("int8", copy=False)
        elif col in {"run", "evt"}:
            if pd is not None:
                frame[col] = pd.to_numeric(frame[col], downcast="integer")
            else:
                frame[col] = frame[col].astype("int64", copy=False)
        else:
            if pd is not None:
                frame[col] = pd.to_numeric(frame[col], downcast="float")
            else:
                frame[col] = frame[col].astype("float32", copy=False)
    return frame


def parse_hidden_layers(text: str) -> list[int]:
    hidden = [int(x) for x in text.split(",") if x.strip()]
    if not hidden:
        raise SystemExit("--hidden-layers must contain at least one layer width")
    if any(width <= 0 for width in hidden):
        raise SystemExit(f"Hidden-layer widths must be positive: {text!r}")
    return hidden


def parse_hidden_layer_grid(primary: str, grid: str | None) -> list[list[int]]:
    specs = [primary]
    if grid:
        specs.extend(item.strip() for item in grid.replace("|", ";").split(";") if item.strip())
    out: list[list[int]] = []
    seen: set[tuple[int, ...]] = set()
    for spec in specs:
        hidden = parse_hidden_layers(spec)
        key = tuple(hidden)
        if key in seen:
            continue
        seen.add(key)
        out.append(hidden)
    return out


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


def expand_input_paths(inputs: list[Path] | None, source: Path | None, check_exists: bool = False) -> list[Path]:
    paths: list[Path] = []
    if source is not None:
        paths.extend(discover_manifest(source))
    for item in inputs or []:
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
    paths = list(dict.fromkeys(paths))
    if not paths:
        raise SystemExit("No input ROOT files supplied. Use --source or --input.")
    print(
        "[trainAuAuPhotonMLP] input path expansion "
        f"source={source} explicit_inputs={len(inputs or [])} paths={len(paths)} "
        f"check_exists={check_exists}",
        flush=True,
    )
    if check_exists:
        missing = [str(path) for path in paths if not path.is_file()]
        if missing:
            preview = "\n  ".join(missing[:20])
            extra = "" if len(missing) <= 20 else f"\n  ... {len(missing) - 20} more"
            raise SystemExit(f"Input ROOT files are missing:\n  {preview}{extra}")
    return paths


def path_group_key(path: Path) -> str:
    parts = list(path.parts)
    for label in ("background", "signal"):
        if label in parts:
            idx = parts.index(label)
            sample = parts[idx + 1] if idx + 1 < len(parts) else "unknownSample"
            return f"{label}:{sample}"
    name = path.name
    if "isSimEmbeddedInclusive" in name or "embeddedJet" in name:
        return "background:unknownSample"
    if "isSimEmbedded" in name or "embeddedPhoton" in name:
        return "signal:unknownSample"
    return "unknown:unknownSample"


def summarize_paths(paths: list[Path]) -> dict:
    groups: dict[str, int] = {}
    for path in paths:
        key = path_group_key(path)
        groups[key] = groups.get(key, 0) + 1
    return {
        "total_files": len(paths),
        "groups": dict(sorted(groups.items())),
    }


def limit_paths_by_group(paths: list[Path], max_files_per_sample: int, random_seed: int) -> tuple[list[Path], dict]:
    before = summarize_paths(paths)
    if max_files_per_sample <= 0:
        return paths, {
            "input_files_before": before,
            "input_files_after": before,
            "max_files_per_sample": int(max_files_per_sample),
            "limited": False,
        }
    grouped: dict[str, list[tuple[int, Path]]] = {}
    for idx, path in enumerate(paths):
        grouped.setdefault(path_group_key(path), []).append((idx, path))
    rng = random.Random(random_seed)
    kept: list[tuple[int, Path]] = []
    chosen_counts: dict[str, int] = {}
    for key, entries in sorted(grouped.items()):
        if len(entries) <= max_files_per_sample:
            chosen = entries
        else:
            pick = sorted(rng.sample(range(len(entries)), max_files_per_sample))
            chosen = [entries[i] for i in pick]
        chosen_counts[key] = len(chosen)
        kept.extend(chosen)
    kept.sort(key=lambda item: item[0])
    limited = [path for _, path in kept]
    return limited, {
        "input_files_before": before,
        "input_files_after": summarize_paths(limited),
        "max_files_per_sample": int(max_files_per_sample),
        "chosen_counts": dict(sorted(chosen_counts.items())),
        "limited": len(limited) != len(paths),
    }

def class_counts(frame, label_branch: str) -> dict:
    if len(frame) == 0:
        return {"0": 0, "1": 0}
    labels = frame[label_branch].to_numpy(dtype="int32")
    return {str(cls): int((labels == cls).sum()) for cls in (0, 1)}


def balanced_row_cap(frame, label_branch: str, max_rows: int, max_rows_per_class: int, random_seed: int):
    import numpy as np

    if len(frame) == 0:
        return frame
    if max_rows_per_class <= 0 and max_rows <= 0:
        return frame
    labels = frame[label_branch].to_numpy(dtype="int32")
    rng = np.random.default_rng(random_seed)
    keep: list[int] = []
    if max_rows_per_class > 0:
        for cls in (0, 1):
            idx = np.flatnonzero(labels == cls)
            if len(idx) > max_rows_per_class:
                idx = rng.choice(idx, size=max_rows_per_class, replace=False)
            keep.extend(int(i) for i in idx)
    else:
        per_class = max(1, max_rows // 2)
        leftovers: list[int] = []
        for cls in (0, 1):
            idx = np.flatnonzero(labels == cls)
            if len(idx) > per_class:
                chosen = rng.choice(idx, size=per_class, replace=False)
                keep.extend(int(i) for i in chosen)
                chosen_set = set(int(i) for i in chosen)
                leftovers.extend(int(i) for i in idx if int(i) not in chosen_set)
            else:
                keep.extend(int(i) for i in idx)
        remaining = max_rows - len(keep)
        if remaining > 0 and leftovers:
            extra = rng.choice(np.asarray(leftovers, dtype="int64"), size=min(remaining, len(leftovers)), replace=False)
            keep.extend(int(i) for i in extra)
    if max_rows > 0 and len(keep) > max_rows:
        keep = [int(i) for i in rng.choice(np.asarray(keep, dtype="int64"), size=max_rows, replace=False)]
    keep = sorted(set(keep))
    if len(keep) == len(frame):
        return frame
    return frame.iloc[keep].copy()


def balanced_pt_bin_row_cap(
    frame,
    label_branch: str,
    bins: list[tuple[float, float]],
    max_rows_per_pt_bin_class: int,
    max_rows: int,
    random_seed: int,
):
    import numpy as np

    if len(frame) == 0 or max_rows_per_pt_bin_class <= 0 or not bins:
        return frame
    labels = frame[label_branch].to_numpy(dtype="int32")
    et = frame["cluster_Et"].to_numpy(dtype="float64")
    rng = np.random.default_rng(random_seed)
    keep: list[int] = []
    bin_reports = []
    for lo, hi in bins:
        for cls in (0, 1):
            idx = np.flatnonzero(np.isfinite(et) & (et >= lo) & (et < hi) & (labels == cls))
            before = int(len(idx))
            if before > max_rows_per_pt_bin_class:
                idx = rng.choice(idx, size=max_rows_per_pt_bin_class, replace=False)
            keep.extend(int(i) for i in idx)
            bin_reports.append(
                {
                    "pt_bin": bin_label(lo, hi),
                    "class": int(cls),
                    "before": before,
                    "after": int(len(idx)),
                }
            )
    keep = sorted(set(keep))
    if max_rows > 0 and len(keep) > max_rows:
        keep = [int(i) for i in rng.choice(np.asarray(keep, dtype="int64"), size=max_rows, replace=False)]
        keep = sorted(set(keep))
    if len(keep) == len(frame):
        return frame
    capped = frame.iloc[keep].copy()
    print(
        "[trainAuAuPhotonMLP] pT-bin row cap "
        f"bins={[bin_label(lo, hi) for lo, hi in bins]} "
        f"max_rows_per_pt_bin_class={max_rows_per_pt_bin_class} "
        f"before={len(frame)} after={len(capped)} reports={bin_reports}",
        flush=True,
    )
    return capped


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(json_ready(payload), indent=2, sort_keys=True) + "\n")


def load_frame(paths: list[Path], tree_name: str, required_columns: list[str], optional_columns: list[str] | None = None, args=None):
    try:
        import pandas as pd
        import uproot
    except ImportError as exc:
        raise SystemExit(
            "MLP training needs pandas and uproot. On SDCC use "
            "RJ_ML_PYTHON=/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python."
        ) from exc

    frames = []
    required = sorted(set(required_columns))
    optional_columns = sorted(set(optional_columns or []))
    seen_optional: set[str] = set()
    read_report = {
        "files_read": 0,
        "rows_read": 0,
        "rows_after_read_filter": 0,
        "read_filter": {
            "pt_range": getattr(args, "pt_range", "") if args is not None else "",
            "centrality_range": getattr(args, "centrality_range", "") if args is not None else "",
            "enabled": bool(args is not None and (getattr(args, "pt_range", "") or getattr(args, "centrality_range", ""))),
        },
        "groups": {},
    }
    pt_range = parse_range(args.pt_range) if args is not None and getattr(args, "pt_range", "") else None
    cent_range = parse_range(args.centrality_range) if args is not None and getattr(args, "centrality_range", "") else None
    for idx, path in enumerate(paths, 1):
        if idx == 1 or idx == len(paths) or idx % 100 == 0:
            print(f"[trainAuAuPhotonMLP] reading {idx}/{len(paths)}: {path}", flush=True)
        with uproot.open(path) as root_file:
            try:
                tree = root_file[tree_name]
            except Exception as exc:
                raise SystemExit(f"{path} does not contain tree {tree_name}") from exc
            keys = set(tree.keys())
            missing = [col for col in required if col not in keys]
            if missing:
                raise SystemExit(f"{path}:{tree_name} missing required branches: {', '.join(missing)}")
            present_optional = [col for col in optional_columns if col in keys]
            seen_optional.update(present_optional)
            chunk = tree.arrays(required + present_optional, library="pd")
            rows_before_filter = int(len(chunk))
            if pt_range is not None and "cluster_Et" in chunk.columns:
                lo, hi = pt_range
                et = chunk["cluster_Et"]
                chunk = chunk[et.notna() & (et >= lo) & (et < hi)]
            if cent_range is not None and "centrality" in chunk.columns:
                lo, hi = cent_range
                cent = chunk["centrality"]
                chunk = chunk[cent.notna() & (cent >= lo) & (cent < hi)]
            if len(chunk) != rows_before_filter:
                chunk = chunk.copy()
            chunk = compact_numeric_frame(chunk, getattr(args, "label_branch", "is_signal") if args is not None else "is_signal")
            frames.append(chunk)
            group = path_group_key(path)
            group_report = read_report["groups"].setdefault(group, {"files": 0, "rows": 0})
            group_report.setdefault("rows_after_read_filter", 0)
            group_report["files"] += 1
            group_report["rows"] += rows_before_filter
            group_report["rows_after_read_filter"] += int(len(chunk))
            read_report["files_read"] += 1
            read_report["rows_read"] += rows_before_filter
            read_report["rows_after_read_filter"] += int(len(chunk))
    frame = pd.concat(frames, ignore_index=True, copy=False)
    for col in optional_columns:
        if col not in frame.columns:
            frame[col] = 1.0
    frame = compact_numeric_frame(frame, getattr(args, "label_branch", "is_signal") if args is not None else "is_signal")
    read_report["groups"] = dict(sorted(read_report["groups"].items()))
    return frame, sorted(seen_optional), read_report


def finite_mask(frame, columns: Iterable[str]):
    import numpy as np

    mask = np.ones(len(frame), dtype=bool)
    for col in columns:
        mask &= np.isfinite(frame[col].to_numpy(dtype="float64"))
    return mask


def normalize_mean_one(values):
    import numpy as np

    values = np.asarray(values, dtype="float64")
    finite = np.isfinite(values) & (values > 0.0)
    if not finite.any():
        return np.ones_like(values, dtype="float64")
    mean = float(values[finite].mean())
    if mean <= 0.0 or not math.isfinite(mean):
        return np.ones_like(values, dtype="float64")
    return np.where(finite, values / mean, 1.0)


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


def pt_bin_weight_factors(weights, et, labels, bins: list[tuple[float, float]], args):
    import numpy as np

    mode = getattr(args, "pt_bin_weight_mode", "none")
    if mode == "none" or not bins:
        return np.ones(len(labels), dtype="float64"), {"mode": mode, "bins": []}
    weights = np.asarray(weights, dtype="float64")
    et = np.asarray(et, dtype="float64")
    labels = np.asarray(labels, dtype="int32")
    targets = parse_bin_weights(getattr(args, "pt_bin_weight_spec", ""), bins, mode)
    factors = np.ones(len(labels), dtype="float64")
    reports = []
    for cls in (0, 1):
        present = []
        for idx, (lo, hi) in enumerate(bins):
            mask = (labels == cls) & np.isfinite(et) & (et >= lo) & (et < hi)
            current = float(weights[mask].sum())
            if current > 0.0:
                present.append((idx, mask, current))
        if not present:
            continue
        total = sum(current for _, _, current in present)
        target_norm = sum(targets[idx] for idx, _, _ in present)
        if total <= 0.0 or target_norm <= 0.0:
            continue
        for idx, mask, current in present:
            desired = total * targets[idx] / target_norm
            factor = desired / current if current > 0.0 else 1.0
            factor = float(np.clip(factor, 1.0 / args.max_flatten_factor, args.max_flatten_factor))
            factors[mask] *= factor
            reports.append(
                {
                    "pt_bin": bin_label(*bins[idx]),
                    "class": int(cls),
                    "entries": int(mask.sum()),
                    "current_weight": current,
                    "target_fraction": float(targets[idx] / target_norm),
                    "factor": factor,
                }
            )
    return factors, {"mode": mode, "target_weights": targets, "bins": reports}


def hard_example_weight_factors(frame, labels, args):
    import numpy as np

    branch = getattr(args, "hard_example_branch", "")
    if not branch:
        return np.ones(len(labels), dtype="float64"), {"enabled": False}
    if branch not in frame.columns:
        return np.ones(len(labels), dtype="float64"), {"enabled": False, "missing_branch": branch}
    raw = frame[branch].to_numpy(dtype="float64")
    finite = np.isfinite(raw)
    if finite.sum() < 10:
        return np.ones(len(labels), dtype="float64"), {"enabled": False, "branch": branch, "reason": "too_few_finite_scores"}
    lo, hi = np.nanpercentile(raw[finite], [1.0, 99.0])
    if not math.isfinite(float(lo)) or not math.isfinite(float(hi)) or hi <= lo:
        norm = np.clip(raw, 0.0, 1.0)
    elif lo < -0.05 or hi > 1.05:
        norm = np.clip((raw - lo) / (hi - lo), 0.0, 1.0)
    else:
        norm = np.clip(raw, 0.0, 1.0)
    labels = np.asarray(labels, dtype="int32")
    power = max(float(args.hard_example_power), 0.25)
    bkg_factor = max(float(args.hard_background_factor), 0.0)
    sig_factor = max(float(args.hard_signal_factor), 0.0)
    factors = np.ones(len(labels), dtype="float64")
    bkg = finite & (labels == 0)
    sig = finite & (labels == 1)
    if bkg.any() and bkg_factor > 0.0:
        factors[bkg] *= 1.0 + bkg_factor * (norm[bkg] ** power)
    if sig.any() and sig_factor > 0.0:
        factors[sig] *= 1.0 + sig_factor * ((1.0 - norm[sig]) ** power)
    for cls in (0, 1):
        cls_mask = labels == cls
        mean = float(np.mean(factors[cls_mask])) if cls_mask.any() else 1.0
        if math.isfinite(mean) and mean > 0.0:
            factors[cls_mask] /= mean
    return factors, {
        "enabled": True,
        "branch": branch,
        "score_p01": float(lo),
        "score_p99": float(hi),
        "hard_background_factor": float(bkg_factor),
        "hard_signal_factor": float(sig_factor),
        "power": float(power),
        "finite_fraction": float(finite.mean()) if len(finite) else math.nan,
        "mean_factor_signal": float(np.mean(factors[labels == 1])) if (labels == 1).any() else math.nan,
        "mean_factor_background": float(np.mean(factors[labels == 0])) if (labels == 0).any() else math.nan,
        "max_factor": float(np.max(factors)) if len(factors) else math.nan,
    }


def normalized_teacher_targets(frame, labels, args):
    import numpy as np

    branch = getattr(args, "distillation_branch", "")
    strength = float(getattr(args, "distillation_strength", 0.0))
    if not branch or strength <= 0.0:
        return None, {"enabled": False}
    if branch not in frame.columns:
        if getattr(args, "require_distillation", False):
            raise SystemExit(f"Distillation branch is required but missing from selected frame: {branch}")
        return None, {"enabled": False, "missing_branch": branch}

    raw = frame[branch].to_numpy(dtype="float64")
    finite = np.isfinite(raw)
    if finite.sum() < 10:
        if getattr(args, "require_distillation", False):
            raise SystemExit(f"Distillation branch {branch} has too few finite values: {int(finite.sum())}")
        return None, {"enabled": False, "branch": branch, "reason": "too_few_finite_scores"}

    lo, hi = np.nanpercentile(raw[finite], [1.0, 99.0])
    if not math.isfinite(float(lo)) or not math.isfinite(float(hi)) or hi <= lo:
        target = np.clip(raw, 0.0, 1.0)
        score_p01, score_p99 = math.nan, math.nan
        normalization = "clip_0_1"
    elif lo < -0.05 or hi > 1.05:
        target = np.clip((raw - lo) / (hi - lo), 0.0, 1.0)
        score_p01, score_p99 = float(lo), float(hi)
        normalization = "percentile_linear_1_99"
    else:
        target = np.clip(raw, 0.0, 1.0)
        score_p01, score_p99 = float(lo), float(hi)
        normalization = "native_0_1"

    temp = float(getattr(args, "distillation_temperature", 1.0))
    if not math.isfinite(temp) or temp <= 0.0:
        raise SystemExit("--distillation-temperature must be finite and positive")
    if abs(temp - 1.0) > 1.0e-12:
        clipped = np.clip(target, 1.0e-6, 1.0 - 1.0e-6)
        logits = np.log(clipped / (1.0 - clipped))
        target = sigmoid(logits / temp)

    valid = finite & np.isfinite(target)
    if getattr(args, "require_distillation", False) and float(valid.mean()) < float(args.distillation_min_finite_fraction):
        raise SystemExit(
            f"Distillation branch {branch} finite target fraction "
            f"{float(valid.mean()):.6g} is below required {args.distillation_min_finite_fraction:.6g}"
        )

    labels = np.asarray(labels, dtype="int32")
    report = {
        "enabled": True,
        "branch": branch,
        "strength": float(strength),
        "temperature": float(temp),
        "normalization": normalization,
        "score_p01": score_p01,
        "score_p99": score_p99,
        "finite_fraction": float(valid.mean()) if len(valid) else math.nan,
        "target_mean_signal": float(np.nanmean(target[valid & (labels == 1)])) if (valid & (labels == 1)).any() else math.nan,
        "target_mean_background": float(np.nanmean(target[valid & (labels == 0)])) if (valid & (labels == 0)).any() else math.nan,
    }
    out = np.full(len(raw), np.nan, dtype="float64")
    out[valid] = target[valid]
    return out, report


def blended_targets(y, teacher, args):
    import numpy as np

    y = np.asarray(y, dtype="float64")
    strength = float(getattr(args, "distillation_strength", 0.0))
    if teacher is None or strength <= 0.0:
        return y
    strength = max(0.0, min(strength, 0.90))
    teacher = np.asarray(teacher, dtype="float64")
    valid = np.isfinite(teacher)
    target = np.array(y, copy=True, dtype="float64")
    target[valid] = (1.0 - strength) * y[valid] + strength * np.clip(teacher[valid], 0.0, 1.0)
    return target


def blended_bce_from_logits(logits, y, weights, teacher, args, params=None):
    return weighted_bce_from_logits(
        logits,
        blended_targets(y, teacher, args),
        weights,
        params,
        args.l2,
    )


def compute_weights(frame, label_branch: str, args):
    import numpy as np

    labels = frame[label_branch].to_numpy(dtype="int32")
    weights = np.ones(len(frame), dtype="float64")
    diagnostics = {
        "use_event_weight": bool(args.use_event_weight),
        "et_reweight": bool(args.et_reweight),
        "eta_reweight": bool(args.eta_reweight),
    }
    if args.use_event_weight and args.weight_branch in frame.columns:
        ew = frame[args.weight_branch].to_numpy(dtype="float64")
        weights *= np.where(np.isfinite(ew) & (ew > 0.0), ew, 1.0)

    class_sums = {str(cls): float(weights[labels == cls].sum()) for cls in (0, 1)}
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

    pt_factors, pt_report = pt_bin_weight_factors(weights, frame["cluster_Et"].to_numpy(), labels, train_pt_bins(args), args)
    weights *= pt_factors
    hard_factors, hard_report = hard_example_weight_factors(frame, labels, args)
    weights *= hard_factors

    finite_positive = np.isfinite(weights) & (weights > 0.0)
    if finite_positive.any():
        median = float(np.median(weights[finite_positive]))
        cap = max(args.max_total_weight_factor * median, 1.0e-12)
        weights = np.where(finite_positive, np.clip(weights, 1.0e-12, cap), 1.0)
    else:
        weights = np.ones(len(frame), dtype="float64")

    diagnostics.update(
        {
            "sum_weight_class0": float(weights[labels == 0].sum()),
            "sum_weight_class1": float(weights[labels == 1].sum()),
            "min_weight": float(np.min(weights)) if len(weights) else math.nan,
            "max_weight": float(np.max(weights)) if len(weights) else math.nan,
            "mean_weight": float(np.mean(weights)) if len(weights) else math.nan,
            "pt_bin_weighting": pt_report,
            "hard_example_weighting": hard_report,
        }
    )
    return weights, diagnostics


def select_training_rows(frame, features: list[str], args):
    import numpy as np

    label = args.label_branch
    needed = sorted(set(features + [label]))
    mask = finite_mask(frame, needed)
    labels = frame[label].to_numpy(dtype="int32")
    mask &= np.isin(labels, [0, 1])
    if args.pt_range:
        lo, hi = parse_range(args.pt_range)
        et = frame["cluster_Et"].to_numpy(dtype="float64")
        mask &= np.isfinite(et) & (et >= lo) & (et < hi)
    if args.centrality_range:
        lo, hi = parse_range(args.centrality_range)
        cent = frame["centrality"].to_numpy(dtype="float64")
        mask &= np.isfinite(cent) & (cent >= lo) & (cent < hi)
    out = frame.loc[mask].copy()
    before_cap = class_counts(out, label)
    pt_bins = train_pt_bins(args)
    if args.max_rows_per_pt_bin_class > 0 and pt_bins:
        out = balanced_pt_bin_row_cap(out, label, pt_bins, args.max_rows_per_pt_bin_class, args.max_rows, args.random_seed)
    else:
        out = balanced_row_cap(out, label, args.max_rows, args.max_rows_per_class, args.random_seed)
    after_cap = class_counts(out, label)
    if before_cap != after_cap:
        print(
            "[trainAuAuPhotonMLP] row cap "
            f"before={before_cap} after={after_cap} "
            f"max_rows={args.max_rows} max_rows_per_class={args.max_rows_per_class} "
            f"max_rows_per_pt_bin_class={args.max_rows_per_pt_bin_class}",
            flush=True,
        )
    return out


def split_masks(frame, label_branch: str, args):
    import numpy as np

    n = len(frame)
    val_frac = float(args.validation_fraction)
    test_frac = float(args.test_fraction)
    if val_frac <= 0.0 or test_frac < 0.0 or val_frac + test_frac >= 0.8:
        raise SystemExit("Use sensible --validation-fraction and --test-fraction values")
    if {"run", "evt"}.issubset(frame.columns):
        run = frame["run"].to_numpy()
        evt = frame["evt"].to_numpy()
        u = np.asarray([stable_unit_interval(args.random_seed, int(r), int(e)) for r, e in zip(run, evt)])
        test = u < test_frac
        val = (u >= test_frac) & (u < test_frac + val_frac)
        train = ~(test | val)
    else:
        y = frame[label_branch].to_numpy(dtype="int32")
        train = np.zeros(n, dtype=bool)
        val = np.zeros(n, dtype=bool)
        test = np.zeros(n, dtype=bool)
        rng = np.random.default_rng(args.random_seed)
        for cls in (0, 1):
            idx = np.flatnonzero(y == cls)
            rng.shuffle(idx)
            n_test = int(round(len(idx) * test_frac))
            n_val = int(round(len(idx) * val_frac))
            test[idx[:n_test]] = True
            val[idx[n_test:n_test + n_val]] = True
            train[idx[n_test + n_val:]] = True
    return train, val, test


def standardize_from_train(x, train_mask, clip: float):
    import numpy as np

    mean = np.nanmean(x[train_mask], axis=0)
    scale = np.nanstd(x[train_mask], axis=0)
    scale = np.where(np.isfinite(scale) & (scale > 1.0e-9), scale, 1.0)
    mean = np.where(np.isfinite(mean), mean, 0.0)
    z = (x - mean) / scale
    if clip > 0:
        z = np.clip(z, -clip, clip)
    return z.astype("float32"), mean.astype("float64"), scale.astype("float64")


def sigmoid(x):
    import numpy as np

    x = np.asarray(x, dtype="float64")
    out = np.empty_like(x, dtype="float64")
    pos = x >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-x[pos]))
    ex = np.exp(x[~pos])
    out[~pos] = ex / (1.0 + ex)
    return out


def weighted_bce_from_logits(logits, y, weights, params=None, l2=0.0):
    import numpy as np

    logits = np.asarray(logits, dtype="float64").reshape(-1)
    y = np.asarray(y, dtype="float64").reshape(-1)
    weights = np.asarray(weights, dtype="float64").reshape(-1)
    w = weights / max(float(weights.sum()), 1.0e-12)
    loss_terms = np.maximum(logits, 0.0) - logits * y + np.log1p(np.exp(-np.abs(logits)))
    loss = float(np.sum(w * loss_terms))
    if params and l2 > 0.0:
        loss += 0.5 * l2 * sum(float((p["W"] ** 2).sum()) for p in params)
    return loss


def init_params(n_in: int, hidden: list[int], seed: int):
    import numpy as np

    rng = np.random.default_rng(seed)
    sizes = [n_in] + hidden + [1]
    params = []
    for n0, n1 in zip(sizes[:-1], sizes[1:]):
        lim = math.sqrt(6.0 / max(1, n0 + n1))
        params.append(
            {
                "W": rng.uniform(-lim, lim, size=(n0, n1)).astype("float64"),
                "b": np.zeros(n1, dtype="float64"),
            }
        )
    return params


def forward(x, params):
    import numpy as np

    activations = [x]
    preacts = []
    a = x
    for idx, layer in enumerate(params):
        z = a @ layer["W"] + layer["b"]
        preacts.append(z)
        if idx + 1 == len(params):
            a = z
        else:
            a = np.tanh(z)
        activations.append(a)
    return activations[-1].reshape(-1), {"activations": activations, "preacts": preacts}


def forward_logits_batched(x, params, batch_size: int = 0):
    """Forward pass for evaluation without retaining full activation caches."""
    import numpy as np

    n = len(x)
    if batch_size <= 0 or n <= batch_size:
        logits, _ = forward(x, params)
        return logits
    out = np.empty(n, dtype="float64")
    for start in range(0, n, batch_size):
        stop = min(start + batch_size, n)
        logits, _ = forward(x[start:stop], params)
        out[start:stop] = logits
    return out


def fixed_eval_indices(n_rows: int, max_rows: int, seed: int):
    import numpy as np

    if max_rows <= 0 or n_rows <= max_rows:
        return np.arange(n_rows, dtype="int64")
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n_rows, size=max_rows, replace=False).astype("int64"))


def train_numpy_mlp(
    x_train,
    y_train,
    w_train,
    x_val,
    y_val,
    w_val,
    feature_names,
    args,
    seed: int,
    hidden: list[int],
    candidate_label: str,
    teacher_train=None,
    teacher_val=None,
):
    import numpy as np

    params = init_params(x_train.shape[1], hidden, seed)
    moments = [{"mW": np.zeros_like(p["W"]), "vW": np.zeros_like(p["W"]), "mb": np.zeros_like(p["b"]), "vb": np.zeros_like(p["b"])} for p in params]
    rng = np.random.default_rng(seed + 17)
    jitter_features = [feature_names.index(f) for f in ("cluster_Et", "centrality") if f in feature_names]
    train_eval_idx = fixed_eval_indices(len(x_train), int(args.train_eval_max_rows), seed + 23)
    x_train_eval = x_train[train_eval_idx]
    y_train_eval = y_train[train_eval_idx]
    w_train_eval = w_train[train_eval_idx]
    teacher_train_eval = teacher_train[train_eval_idx] if teacher_train is not None else None
    best = copy.deepcopy(params)
    best_loss = math.inf
    best_epoch = 0
    history = []
    step = 0
    beta1 = 0.9
    beta2 = 0.999
    eps = 1.0e-8

    for epoch in range(1, args.epochs + 1):
        order = rng.permutation(len(x_train))
        for start in range(0, len(order), args.batch_size):
            idx = order[start:start + args.batch_size]
            xb = np.array(x_train[idx], copy=True)
            if args.conditional_jitter > 0.0 and jitter_features:
                xb[:, jitter_features] += rng.normal(0.0, args.conditional_jitter, size=(len(idx), len(jitter_features)))
            yb = y_train[idx].astype("float64")
            wb = w_train[idx].astype("float64")
            wb_norm = wb / max(float(wb.sum()), 1.0e-12)
            tb = teacher_train[idx] if teacher_train is not None else None
            target_b = blended_targets(yb, tb, args)

            logits, cache = forward(xb, params)
            dz = ((sigmoid(logits) - target_b) * wb_norm).reshape(-1, 1)
            grads = [None] * len(params)
            for layer_idx in reversed(range(len(params))):
                a_prev = cache["activations"][layer_idx]
                grads[layer_idx] = {
                    "W": a_prev.T @ dz + args.l2 * params[layer_idx]["W"],
                    "b": dz.sum(axis=0),
                }
                if layer_idx > 0:
                    da = dz @ params[layer_idx]["W"].T
                    dz = da * (1.0 - np.tanh(cache["preacts"][layer_idx - 1]) ** 2)

            step += 1
            for p, g, m in zip(params, grads, moments):
                m["mW"] = beta1 * m["mW"] + (1.0 - beta1) * g["W"]
                m["vW"] = beta2 * m["vW"] + (1.0 - beta2) * (g["W"] ** 2)
                m["mb"] = beta1 * m["mb"] + (1.0 - beta1) * g["b"]
                m["vb"] = beta2 * m["vb"] + (1.0 - beta2) * (g["b"] ** 2)
                mW_hat = m["mW"] / (1.0 - beta1 ** step)
                vW_hat = m["vW"] / (1.0 - beta2 ** step)
                mb_hat = m["mb"] / (1.0 - beta1 ** step)
                vb_hat = m["vb"] / (1.0 - beta2 ** step)
                p["W"] -= args.learning_rate * mW_hat / (np.sqrt(vW_hat) + eps)
                p["b"] -= args.learning_rate * mb_hat / (np.sqrt(vb_hat) + eps)

        train_logits = forward_logits_batched(x_train_eval, params, int(args.eval_batch_size))
        val_logits = forward_logits_batched(x_val, params, int(args.eval_batch_size))
        train_loss = blended_bce_from_logits(train_logits, y_train_eval, w_train_eval, teacher_train_eval, args, params)
        val_loss = blended_bce_from_logits(val_logits, y_val, w_val, teacher_val, args, params)
        train_truth_loss = weighted_bce_from_logits(train_logits, y_train_eval, w_train_eval, params, args.l2)
        val_truth_loss = weighted_bce_from_logits(val_logits, y_val, w_val, params, args.l2)
        row = {
            "epoch": epoch,
            "train_loss": train_loss,
            "validation_loss": val_loss,
            "train_truth_loss": train_truth_loss,
            "validation_truth_loss": val_truth_loss,
            "validation_auc": auc_score(y_val, sigmoid(val_logits), w_val),
            "train_eval_rows": int(len(train_eval_idx)),
            "validation_eval_rows": int(len(x_val)),
        }
        history.append(row)
        if epoch == 1 or epoch % max(1, args.progress_every) == 0:
            print(
                "[trainAuAuPhotonMLP] "
                f"candidate={candidate_label} epoch={epoch} train_loss={train_loss:.6g} "
                f"val_loss={val_loss:.6g} truth_val_loss={val_truth_loss:.6g} "
                f"val_auc={row['validation_auc']:.5f}",
                flush=True,
            )
        if val_loss + args.min_delta < best_loss:
            best_loss = val_loss
            best_epoch = epoch
            best = copy.deepcopy(params)
        elif epoch - best_epoch >= args.patience:
            print(f"[trainAuAuPhotonMLP] candidate={candidate_label} early stopping at epoch {epoch}; best_epoch={best_epoch}", flush=True)
            break
    best_row = min(history, key=lambda row: row["validation_loss"]) if history else {}
    return best, history, {
        "seed": int(seed),
        "hidden_layers": [int(x) for x in hidden],
        "best_epoch": int(best_epoch),
        "best_validation_loss": float(best_loss),
        "best_validation_auc": float(best_row.get("validation_auc", math.nan)),
        "trained_epochs": int(history[-1]["epoch"]) if history else 0,
    }


def auc_score(y, score, weights=None) -> float:
    import numpy as np

    y = np.asarray(y, dtype="int32")
    score = np.asarray(score, dtype="float64")
    mask = np.isfinite(score) & np.isin(y, [0, 1])
    y = y[mask]
    score = score[mask]
    if weights is None:
        weights = np.ones_like(score, dtype="float64")
    else:
        weights = np.asarray(weights, dtype="float64")[mask]
    pos = y == 1
    neg = y == 0
    w_pos = float(weights[pos].sum())
    w_neg = float(weights[neg].sum())
    if w_pos <= 0.0 or w_neg <= 0.0:
        return math.nan
    order = np.argsort(score, kind="mergesort")
    y_s = y[order]
    w_s = weights[order]
    score_s = score[order]
    auc = 0.0
    neg_before = 0.0
    i = 0
    while i < len(score_s):
        j = i + 1
        while j < len(score_s) and score_s[j] == score_s[i]:
            j += 1
        w_pos_tie = float(w_s[i:j][y_s[i:j] == 1].sum())
        w_neg_tie = float(w_s[i:j][y_s[i:j] == 0].sum())
        auc += w_pos_tie * (neg_before + 0.5 * w_neg_tie)
        neg_before += w_neg_tie
        i = j
    return float(auc / (w_pos * w_neg))


def threshold_for_signal_efficiency(y, score, target: float):
    import numpy as np

    y = np.asarray(y, dtype="int32")
    score = np.asarray(score, dtype="float64")
    mask = np.isfinite(score) & np.isin(y, [0, 1])
    sig = score[mask & (y == 1)]
    bkg = score[mask & (y == 0)]
    if sig.size == 0 or bkg.size == 0:
        return None
    threshold = float(np.quantile(sig, max(0.0, min(1.0, 1.0 - target))))
    return {
        "threshold": threshold,
        "signal_efficiency": float(np.mean(sig > threshold)),
        "background_fake_rate": float(np.mean(bkg > threshold)),
        "signal_entries": int(sig.size),
        "background_entries": int(bkg.size),
    }


def calibration_report(y, prob, weights=None, nbins: int = 10):
    import numpy as np

    y = np.asarray(y, dtype="float64")
    prob = np.asarray(prob, dtype="float64")
    if weights is None:
        weights = np.ones_like(prob)
    else:
        weights = np.asarray(weights, dtype="float64")
    mask = np.isfinite(prob) & np.isfinite(y) & np.isfinite(weights) & (weights > 0.0)
    y = y[mask]
    prob = np.clip(prob[mask], 0.0, 1.0)
    weights = weights[mask]
    if len(prob) == 0:
        return {"ece": math.nan, "brier": math.nan, "bins": []}
    bins = np.linspace(0.0, 1.0, nbins + 1)
    rows = []
    ece = 0.0
    wtot = float(weights.sum())
    for lo, hi in zip(bins[:-1], bins[1:]):
        sel = (prob >= lo) & (prob < hi if hi < 1.0 else prob <= hi)
        if not sel.any():
            continue
        w = weights[sel]
        conf = float(np.average(prob[sel], weights=w))
        acc = float(np.average(y[sel], weights=w))
        frac = float(w.sum() / wtot)
        ece += frac * abs(acc - conf)
        rows.append({"lo": float(lo), "hi": float(hi), "entries": int(sel.sum()), "confidence": conf, "signal_fraction": acc})
    brier = float(np.average((prob - y) ** 2, weights=weights))
    return {"ece": float(ece), "brier": brier, "bins": rows}


def fit_temperature(logits, y, weights):
    import numpy as np

    candidates = np.unique(np.concatenate([np.asarray([1.0]), np.exp(np.linspace(math.log(0.4), math.log(5.0), 80))]))
    best_t = 1.0
    best_loss = math.inf
    for t in candidates:
        loss = weighted_bce_from_logits(np.asarray(logits) / float(t), y, weights)
        if loss < best_loss:
            best_loss = loss
            best_t = float(t)
    return best_t, best_loss


def write_history(path: Path, history: list[dict]) -> None:
    if not history:
        return
    keys = list(history[0])
    with path.open("w") as out:
        out.write(",".join(keys) + "\n")
        for row in history:
            out.write(",".join(str(row.get(k, "")) for k in keys) + "\n")


def artifact_from_params(product: str, spec: dict, features: list[str], mean, scale, params, temperature: float, args, report: dict, hidden_layers: list[int]):
    layers = []
    for layer in params:
        # C++ evaluator expects weights grouped by output node, then input node.
        layers.append(
            {
                "weights": layer["W"].T.tolist(),
                "bias": layer["b"].tolist(),
            }
        )
    return {
        "schema": "RJ_AUAU_TIGHT_MLP_V1",
        "product": product,
        "variant": spec["variant"],
        "label": spec["label"],
        "features": features,
        "preprocessing": {
            "mean": mean.tolist(),
            "scale": scale.tolist(),
            "clip": float(args.input_clip),
        },
        "activation": "tanh",
        "output_temperature": float(temperature),
        "score_definition": "sigmoid(logit / output_temperature)",
        "layers": layers,
        "training": {
            "pt_range": args.pt_range,
            "hidden_layers": [int(x) for x in hidden_layers],
            "epochs": int(args.epochs),
            "eval_batch_size": int(args.eval_batch_size),
            "train_eval_max_rows": int(args.train_eval_max_rows),
            "learning_rate": float(args.learning_rate),
            "l2": float(args.l2),
            "conditional_jitter": float(args.conditional_jitter),
            "max_files_per_sample": int(args.max_files_per_sample),
            "max_rows": int(args.max_rows),
            "max_rows_per_class": int(args.max_rows_per_class),
            "restarts": int(args.restarts),
            "selection_metric": args.selection_metric,
            "selection_target_signal_efficiency": float(args.selection_target_signal_efficiency),
            "random_seed": int(args.random_seed),
            "centrality_range": args.centrality_range,
            "train_pt_bins": args.train_pt_bins,
            "max_rows_per_pt_bin_class": int(args.max_rows_per_pt_bin_class),
            "pt_bin_weight_mode": args.pt_bin_weight_mode,
            "pt_bin_weight_spec": args.pt_bin_weight_spec,
            "highpt_selection_weights": args.highpt_selection_weights,
            "hard_example_branch": args.hard_example_branch,
            "hard_background_factor": float(args.hard_background_factor),
            "hard_signal_factor": float(args.hard_signal_factor),
            "hard_example_power": float(args.hard_example_power),
            "distillation_branch": args.distillation_branch,
            "distillation_strength": float(args.distillation_strength),
            "distillation_temperature": float(args.distillation_temperature),
            "distillation_min_finite_fraction": float(args.distillation_min_finite_fraction),
            "require_distillation": bool(args.require_distillation),
        },
        "report": report,
    }


def load_mlp_artifact(path: Path):
    data = json.loads(Path(path).read_text())
    if data.get("schema") != "RJ_AUAU_TIGHT_MLP_V1":
        raise ValueError(f"Unsupported MLP artifact schema in {path}: {data.get('schema')}")
    return data


def predict_mlp_array(model: dict, x_raw):
    import numpy as np

    features = model["features"]
    prep = model.get("preprocessing", {})
    mean = np.asarray(prep.get("mean", [0.0] * len(features)), dtype="float64")
    scale = np.asarray(prep.get("scale", [1.0] * len(features)), dtype="float64")
    clip = float(prep.get("clip", 8.0))
    x = (np.asarray(x_raw, dtype="float64") - mean) / np.where(np.abs(scale) > 1.0e-12, scale, 1.0)
    if clip > 0:
        x = np.clip(x, -clip, clip)
    a = x
    layers = model["layers"]
    for idx, layer in enumerate(layers):
        w = np.asarray(layer["weights"], dtype="float64").T
        b = np.asarray(layer["bias"], dtype="float64")
        z = a @ w + b
        a = z if idx + 1 == len(layers) else np.tanh(z)
    logits = a.reshape(-1)
    temperature = float(model.get("output_temperature", 1.0))
    if not math.isfinite(temperature) or temperature <= 0.0:
        temperature = 1.0
    return sigmoid(logits / temperature)


def make_report(y, prob, weights):
    return {
        "auc": auc_score(y, prob, weights),
        "thresholds": {
            f"wp{int(round(target * 100)):03d}": threshold_for_signal_efficiency(y, prob, target)
            for target in (0.50, 0.70, 0.80, 0.90, 0.95)
        },
        "calibration": calibration_report(y, prob, weights),
    }


def binned_score_metrics(et, y, prob, weights, bins: list[tuple[float, float]], target_eff: float):
    import numpy as np

    et = np.asarray(et, dtype="float64")
    rows = []
    for lo, hi in bins:
        mask = np.isfinite(et) & (et >= lo) & (et < hi)
        y_bin = y[mask]
        prob_bin = prob[mask]
        weights_bin = weights[mask]
        wp = threshold_for_signal_efficiency(y_bin, prob_bin, target_eff)
        rows.append(
            {
                "lo": float(lo),
                "hi": float(hi),
                "entries": int(mask.sum()),
                "signal_entries": int((y_bin == 1).sum()),
                "background_entries": int((y_bin == 0).sum()),
                "auc": auc_score(y_bin, prob_bin, weights_bin),
                f"wp{int(round(target_eff * 100)):03d}": wp,
            }
        )
    return rows


def highpt_selection_summary(bin_rows: list[dict], bins: list[tuple[float, float]], target_eff: float, args):
    weights = parse_bin_weights(getattr(args, "highpt_selection_weights", ""), bins, "highpt")
    weighted_fake = 0.0
    weighted_auc = 0.0
    weight_sum_fake = 0.0
    weight_sum_auc = 0.0
    worst_fake = -math.inf
    worst_bin = None
    wp_key = f"wp{int(round(target_eff * 100)):03d}"
    for row, weight in zip(bin_rows, weights):
        fake = (row.get(wp_key) or {}).get("background_fake_rate")
        auc = row.get("auc")
        if isinstance(fake, (int, float)) and math.isfinite(float(fake)):
            weighted_fake += float(weight) * float(fake)
            weight_sum_fake += float(weight)
            if float(fake) > worst_fake:
                worst_fake = float(fake)
                worst_bin = bin_label(row["lo"], row["hi"])
        if isinstance(auc, (int, float)) and math.isfinite(float(auc)):
            weighted_auc += float(weight) * float(auc)
            weight_sum_auc += float(weight)
    return {
        "target_signal_efficiency": float(target_eff),
        "bin_weights": weights,
        "weighted_fake_rate": float(weighted_fake / weight_sum_fake) if weight_sum_fake > 0 else math.inf,
        "worst_fake_rate": float(worst_fake) if math.isfinite(worst_fake) else math.inf,
        "worst_fake_bin": worst_bin,
        "weighted_auc": float(weighted_auc / weight_sum_auc) if weight_sum_auc > 0 else -math.inf,
    }


def train_product(product: str, frame, args, outdir: Path):
    import numpy as np

    spec = MODEL_SPECS[product]
    features = list(spec["features"])
    selected = select_training_rows(frame, features, args)
    labels = selected[args.label_branch].to_numpy(dtype="int32")
    counts = {str(cls): int((labels == cls).sum()) for cls in (0, 1)}
    if min(counts.values()) < args.min_rows_per_class:
        raise SystemExit(f"{product}: not enough rows per class after cuts: {counts}")
    weights, weight_report = compute_weights(selected, args.label_branch, args)
    teacher_targets, teacher_report = normalized_teacher_targets(selected, labels, args)
    train_mask, val_mask, test_mask = split_masks(selected, args.label_branch, args)
    if args.verbose_diagnostics:
        split_counts = {
            name: {
                "rows": int(mask.sum()),
                "class0": int(((labels == 0) & mask).sum()),
                "class1": int(((labels == 1) & mask).sum()),
            }
            for name, mask in [("train", train_mask), ("validation", val_mask), ("test", test_mask)]
        }
        print(
            "[trainAuAuPhotonMLP] "
            f"product={product} features={len(features)} selected_rows={len(selected)} "
            f"class_counts={counts} splits={split_counts}",
            flush=True,
        )
        if teacher_report.get("enabled"):
            print(
                "[trainAuAuPhotonMLP] "
                f"product={product} distillation=branch:{teacher_report.get('branch')} "
                f"strength:{teacher_report.get('strength')} temp:{teacher_report.get('temperature')} "
                f"finite:{teacher_report.get('finite_fraction')}",
                flush=True,
            )
        print(
            "[trainAuAuPhotonMLP] "
            f"product={product} weights=sum0:{weight_report['sum_weight_class0']:.6g} "
            f"sum1:{weight_report['sum_weight_class1']:.6g} "
            f"min:{weight_report['min_weight']:.6g} mean:{weight_report['mean_weight']:.6g} "
            f"max:{weight_report['max_weight']:.6g}",
            flush=True,
        )
    for name, mask in [("train", train_mask), ("validation", val_mask), ("test", test_mask)]:
        y = labels[mask]
        if len(np.unique(y)) < 2:
            raise SystemExit(f"{product}: {name} split is missing a class")

    x_raw = np.column_stack([selected[f].to_numpy(dtype="float32") for f in features])
    x, mean, scale = standardize_from_train(x_raw, train_mask, args.input_clip)
    del x_raw
    hidden_grid = parse_hidden_layer_grid(args.hidden_layers, args.hidden_layer_grid)
    target_eff = float(args.selection_target_signal_efficiency)
    metric_pt_bins = train_pt_bins(args)
    et_values = selected["cluster_Et"].to_numpy(dtype="float64")
    candidate_entries = []
    candidate_reports = []
    history_stem = outdir / Path(spec["filename"]).with_suffix("")
    for arch_idx, hidden in enumerate(hidden_grid, 1):
        for restart_idx in range(max(1, int(args.restarts))):
            seed = stable_seed(args.random_seed, product, tuple(hidden), restart_idx)
            candidate_label = f"{product}:a{arch_idx}r{restart_idx + 1}"
            params_i, history_i, train_summary = train_numpy_mlp(
                x[train_mask],
                labels[train_mask],
                weights[train_mask],
                x[val_mask],
                labels[val_mask],
                weights[val_mask],
                features,
                args,
                seed,
                hidden,
                candidate_label,
                teacher_targets[train_mask] if teacher_targets is not None else None,
                teacher_targets[val_mask] if teacher_targets is not None else None,
            )
            val_logits_i = forward_logits_batched(x[val_mask], params_i, int(args.eval_batch_size))
            temperature_i, temp_loss_i = fit_temperature(val_logits_i, labels[val_mask], weights[val_mask])
            val_prob_i = sigmoid(val_logits_i / temperature_i)
            val_report_i = make_report(labels[val_mask], val_prob_i, weights[val_mask])
            selection_wp = threshold_for_signal_efficiency(labels[val_mask], val_prob_i, target_eff) or {}
            pt_bin_report_i = binned_score_metrics(
                et_values[val_mask],
                labels[val_mask],
                val_prob_i,
                weights[val_mask],
                metric_pt_bins,
                target_eff,
            )
            highpt_summary_i = highpt_selection_summary(pt_bin_report_i, metric_pt_bins, target_eff, args)
            candidate_report = {
                "candidate": candidate_label,
                "architecture_index": int(arch_idx),
                "restart_index": int(restart_idx + 1),
                "hidden_layers": [int(width) for width in hidden],
                "seed": int(seed),
                "trained_epochs": int(train_summary["trained_epochs"]),
                "best_epoch": int(train_summary["best_epoch"]),
                "raw_best_validation_loss": float(train_summary["best_validation_loss"]),
                "raw_best_validation_auc": float(train_summary["best_validation_auc"]),
                "temperature": float(temperature_i),
                "temperature_validation_loss": float(temp_loss_i),
                "validation_auc": float(val_report_i["auc"]),
                "validation_ece": float(val_report_i["calibration"].get("ece", math.nan)),
                "selection_target_signal_efficiency": float(target_eff),
                "selection_threshold": selection_wp.get("threshold"),
                "selection_fake_rate": selection_wp.get("background_fake_rate"),
                "pt_bin_validation": pt_bin_report_i,
                "highpt_selection": highpt_summary_i,
            }
            write_history(
                history_stem.with_name(f"{history_stem.name}.candidate_a{arch_idx}_r{restart_idx + 1}.history.csv"),
                history_i,
            )
            candidate_entries.append(
                {
                    "params": params_i,
                    "history": history_i,
                    "temperature": temperature_i,
                    "temperature_loss": temp_loss_i,
                    "hidden_layers": hidden,
                    "report": candidate_report,
                }
            )
            candidate_reports.append(candidate_report)
            print(
                "[trainAuAuPhotonMLP] "
                f"candidate={candidate_label} hidden={hidden} "
                f"val_auc={candidate_report['validation_auc']:.5f} "
                f"wp{int(round(target_eff * 100)):03d}_fake={candidate_report['selection_fake_rate']} "
                f"highpt_fake={highpt_summary_i['weighted_fake_rate']} "
                f"worst_pt_fake={highpt_summary_i['worst_fake_rate']} "
                f"temp={temperature_i:.5g}",
                flush=True,
            )

    def selection_key(entry):
        report_i = entry["report"]
        fake = report_i.get("selection_fake_rate")
        auc = report_i.get("validation_auc")
        loss = report_i.get("temperature_validation_loss")
        ece = report_i.get("validation_ece")
        fake_key = fake if isinstance(fake, (int, float)) and math.isfinite(float(fake)) else math.inf
        auc_key = auc if isinstance(auc, (int, float)) and math.isfinite(float(auc)) else -math.inf
        loss_key = loss if isinstance(loss, (int, float)) and math.isfinite(float(loss)) else math.inf
        ece_key = ece if isinstance(ece, (int, float)) and math.isfinite(float(ece)) else math.inf
        highpt = report_i.get("highpt_selection") or {}
        highpt_fake = highpt.get("weighted_fake_rate")
        highpt_worst = highpt.get("worst_fake_rate")
        highpt_auc = highpt.get("weighted_auc")
        highpt_fake_key = highpt_fake if isinstance(highpt_fake, (int, float)) and math.isfinite(float(highpt_fake)) else math.inf
        highpt_worst_key = highpt_worst if isinstance(highpt_worst, (int, float)) and math.isfinite(float(highpt_worst)) else math.inf
        highpt_auc_key = highpt_auc if isinstance(highpt_auc, (int, float)) and math.isfinite(float(highpt_auc)) else -math.inf
        if args.selection_metric == "highpt_wp80":
            return (highpt_fake_key, highpt_worst_key, -highpt_auc_key, fake_key, loss_key, ece_key)
        if args.selection_metric == "worst_pt_wp80":
            return (highpt_worst_key, highpt_fake_key, -highpt_auc_key, fake_key, loss_key, ece_key)
        if args.selection_metric == "highpt_auc":
            return (-highpt_auc_key, highpt_fake_key, highpt_worst_key, fake_key, loss_key, ece_key)
        if args.selection_metric == "validation_auc":
            return (-auc_key, fake_key, loss_key, ece_key)
        if args.selection_metric == "validation_loss":
            return (loss_key, fake_key, -auc_key, ece_key)
        return (fake_key, -auc_key, loss_key, ece_key)

    best_entry = min(candidate_entries, key=selection_key)
    params = best_entry["params"]
    history = best_entry["history"]
    hidden_layers = [int(width) for width in best_entry["hidden_layers"]]
    temperature = float(best_entry["temperature"])
    temp_loss = float(best_entry["temperature_loss"])
    selected_candidate = best_entry["report"]["candidate"]
    print(
        "[trainAuAuPhotonMLP] "
        f"selected product={product} candidate={selected_candidate} "
        f"metric={args.selection_metric} hidden={hidden_layers} "
        f"val_auc={best_entry['report']['validation_auc']:.5f} "
        f"selection_fake_rate={best_entry['report']['selection_fake_rate']}",
        flush=True,
    )
    train_idx = np.flatnonzero(train_mask)
    train_eval_local = fixed_eval_indices(len(train_idx), int(args.train_eval_max_rows), stable_seed(args.random_seed, product, "final_train_report"))
    train_eval_idx = train_idx[train_eval_local]
    train_logits = forward_logits_batched(x[train_eval_idx], params, int(args.eval_batch_size))
    val_logits = forward_logits_batched(x[val_mask], params, int(args.eval_batch_size))
    test_logits = forward_logits_batched(x[test_mask], params, int(args.eval_batch_size))
    train_prob = sigmoid(train_logits / temperature)
    val_prob = sigmoid(val_logits / temperature)
    test_prob = sigmoid(test_logits / temperature)
    val_pt_bins = binned_score_metrics(et_values[val_mask], labels[val_mask], val_prob, weights[val_mask], metric_pt_bins, target_eff)
    test_pt_bins = binned_score_metrics(et_values[test_mask], labels[test_mask], test_prob, weights[test_mask], metric_pt_bins, target_eff)
    report = {
        "status": "trained",
        "product": product,
        "variant": spec["variant"],
        "features": features,
        "rows": {
            "total": int(len(selected)),
            "train": int(train_mask.sum()),
            "validation": int(val_mask.sum()),
            "test": int(test_mask.sum()),
            "class_counts": counts,
        },
        "weights": weight_report,
        "distillation": teacher_report,
        "temperature": {"value": float(temperature), "validation_loss": float(temp_loss)},
        "model_selection": {
            "metric": args.selection_metric,
            "target_signal_efficiency": float(target_eff),
            "hidden_layer_grid": [[int(x) for x in hidden] for hidden in hidden_grid],
            "restarts": int(args.restarts),
            "selected_candidate": selected_candidate,
            "candidates": candidate_reports,
            "pt_bins": [{"lo": float(lo), "hi": float(hi)} for lo, hi in metric_pt_bins],
        },
        "train": make_report(labels[train_eval_idx], train_prob, weights[train_eval_idx]),
        "validation": make_report(labels[val_mask], val_prob, weights[val_mask]),
        "test": make_report(labels[test_mask], test_prob, weights[test_mask]),
        "validation_pt_bins": val_pt_bins,
        "test_pt_bins": test_pt_bins,
        "train_report_rows": {
            "evaluated": int(len(train_eval_idx)),
            "total_train_rows": int(train_mask.sum()),
            "mode": "fixed_subsample" if len(train_eval_idx) != int(train_mask.sum()) else "full_train",
        },
    }
    artifact = artifact_from_params(product, spec, features, mean, scale, params, temperature, args, report, hidden_layers)
    out_json = outdir / spec["filename"]
    out_json.write_text(json.dumps(json_ready(artifact), indent=2, sort_keys=True) + "\n")
    meta = out_json.with_suffix(".metadata.json")
    meta.write_text(json.dumps(json_ready(report), indent=2, sort_keys=True) + "\n")
    write_json(out_json.with_suffix(".model_selection.json"), report["model_selection"])
    write_history(out_json.with_suffix(".history.csv"), history)
    print(
        "[trainAuAuPhotonMLP] "
        f"product={product} test_auc={report['test']['auc']:.5f} "
        f"val_auc={report['validation']['auc']:.5f} artifact={out_json}",
        flush=True,
    )
    return {
        "product": product,
        "variant": spec["variant"],
        "label": spec["label"],
        "features": features,
        "output_json": str(out_json),
        "metadata": str(meta),
        "report": report,
    }


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", type=Path, default=None, help="Extraction run root with manifests/training_roots.list.")
    ap.add_argument("--input", nargs="+", type=Path, default=None, help="ROOT files or @manifest.list files.")
    ap.add_argument(
        "--check-input-files",
        action="store_true",
        default=os.environ.get("RJ_AUAU_MLP_CHECK_INPUT_FILES", "0") == "1",
        help=(
            "Stat every input ROOT path before training. Off by default to avoid "
            "long GPFS stalls; the ROOT reader will still fail on the first "
            "missing file."
        ),
    )
    ap.add_argument("--tree", default="AuAuPhotonIDTrainingTree")
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--products", default="all", help="Comma list, 'primary', or 'all'.")
    ap.add_argument("--pt-range", default="15:35")
    ap.add_argument("--centrality-range", default=os.environ.get("RJ_AUAU_MLP_TRAIN_CENTRALITY_RANGE", ""))
    ap.add_argument(
        "--train-pt-bins",
        default=os.environ.get("RJ_AUAU_MLP_TRAIN_PT_BINS", ""),
        help="Optional pT edges or lo:hi bins for pT-aware caps/weights/selection diagnostics.",
    )
    ap.add_argument("--label-branch", default="is_signal")
    ap.add_argument("--weight-branch", default="event_weight")
    ap.add_argument("--use-event-weight", dest="use_event_weight", action="store_true", default=True)
    ap.add_argument("--no-event-weight", dest="use_event_weight", action="store_false")
    ap.add_argument("--et-reweight", dest="et_reweight", action="store_true", default=True)
    ap.add_argument("--no-et-reweight", dest="et_reweight", action="store_false")
    ap.add_argument("--eta-reweight", dest="eta_reweight", action="store_true", default=True)
    ap.add_argument("--no-eta-reweight", dest="eta_reweight", action="store_false")
    ap.add_argument("--flatten-bins", type=int, default=20)
    ap.add_argument("--max-flatten-factor", type=float, default=8.0)
    ap.add_argument("--max-total-weight-factor", type=float, default=50.0)
    ap.add_argument("--min-rows-per-class", type=int, default=100)
    ap.add_argument("--max-files-per-sample", type=int, default=int(os.environ.get("RJ_AUAU_MLP_TRAIN_MAX_FILES_PER_SAMPLE", "0")))
    ap.add_argument("--max-rows", type=int, default=int(os.environ.get("RJ_AUAU_MLP_TRAIN_MAX_ROWS", "0")))
    ap.add_argument("--max-rows-per-class", type=int, default=int(os.environ.get("RJ_AUAU_MLP_TRAIN_MAX_ROWS_PER_CLASS", "0")))
    ap.add_argument("--max-rows-per-pt-bin-class", type=int, default=int(os.environ.get("RJ_AUAU_MLP_TRAIN_MAX_ROWS_PER_PT_BIN_CLASS", "0")))
    ap.add_argument("--pt-bin-weight-mode", choices=["none", "equal", "highpt"], default=os.environ.get("RJ_AUAU_MLP_TRAIN_PT_BIN_WEIGHT_MODE", "none"))
    ap.add_argument("--pt-bin-weight-spec", default=os.environ.get("RJ_AUAU_MLP_TRAIN_PT_BIN_WEIGHT_SPEC", ""))
    ap.add_argument("--hard-example-branch", default=os.environ.get("RJ_AUAU_MLP_TRAIN_HARD_EXAMPLE_BRANCH", ""))
    ap.add_argument("--hard-background-factor", type=float, default=float(os.environ.get("RJ_AUAU_MLP_TRAIN_HARD_BACKGROUND_FACTOR", "0.0")))
    ap.add_argument("--hard-signal-factor", type=float, default=float(os.environ.get("RJ_AUAU_MLP_TRAIN_HARD_SIGNAL_FACTOR", "0.0")))
    ap.add_argument("--hard-example-power", type=float, default=float(os.environ.get("RJ_AUAU_MLP_TRAIN_HARD_EXAMPLE_POWER", "2.0")))
    ap.add_argument("--distillation-branch", default=os.environ.get("RJ_AUAU_MLP_TRAIN_DISTILLATION_BRANCH", ""))
    ap.add_argument("--distillation-strength", type=float, default=float(os.environ.get("RJ_AUAU_MLP_TRAIN_DISTILLATION_STRENGTH", "0.0")))
    ap.add_argument("--distillation-temperature", type=float, default=float(os.environ.get("RJ_AUAU_MLP_TRAIN_DISTILLATION_TEMPERATURE", "1.0")))
    ap.add_argument("--distillation-min-finite-fraction", type=float, default=float(os.environ.get("RJ_AUAU_MLP_TRAIN_DISTILLATION_MIN_FINITE_FRACTION", "0.95")))
    ap.add_argument(
        "--require-distillation",
        action="store_true",
        default=os.environ.get("RJ_AUAU_MLP_TRAIN_REQUIRE_DISTILLATION", "0") == "1",
        help="Fail if the requested distillation branch is missing or too sparse.",
    )
    ap.add_argument("--validation-fraction", type=float, default=0.20)
    ap.add_argument("--test-fraction", type=float, default=0.10)
    ap.add_argument("--random-seed", type=int, default=137)
    ap.add_argument("--hidden-layers", default="48,24")
    ap.add_argument("--hidden-layer-grid", default=os.environ.get("RJ_AUAU_MLP_TRAIN_HIDDEN_LAYER_GRID", ""))
    ap.add_argument("--restarts", type=int, default=int(os.environ.get("RJ_AUAU_MLP_TRAIN_RESTARTS", "1")))
    ap.add_argument(
        "--selection-metric",
        choices=["wp_fake_rate", "validation_auc", "validation_loss", "highpt_wp80", "worst_pt_wp80", "highpt_auc"],
        default=os.environ.get("RJ_AUAU_MLP_TRAIN_SELECTION_METRIC", "wp_fake_rate"),
    )
    ap.add_argument("--selection-target-signal-efficiency", type=float, default=float(os.environ.get("RJ_AUAU_MLP_TRAIN_SELECTION_TARGET_SIGNAL_EFF", "0.80")))
    ap.add_argument("--highpt-selection-weights", default=os.environ.get("RJ_AUAU_MLP_TRAIN_HIGHPT_SELECTION_WEIGHTS", ""))
    ap.add_argument("--epochs", type=int, default=80)
    ap.add_argument("--batch-size", type=int, default=4096)
    ap.add_argument("--eval-batch-size", type=int, default=int(os.environ.get("RJ_AUAU_MLP_EVAL_BATCH_SIZE", "131072")))
    ap.add_argument("--train-eval-max-rows", type=int, default=int(os.environ.get("RJ_AUAU_MLP_TRAIN_EVAL_MAX_ROWS", "0")))
    ap.add_argument("--learning-rate", type=float, default=1.0e-3)
    ap.add_argument("--l2", type=float, default=1.0e-4)
    ap.add_argument("--patience", type=int, default=12)
    ap.add_argument("--min-delta", type=float, default=1.0e-5)
    ap.add_argument("--input-clip", type=float, default=8.0)
    ap.add_argument("--conditional-jitter", type=float, default=0.03)
    ap.add_argument("--progress-every", type=int, default=5)
    ap.add_argument("--verbose-diagnostics", dest="verbose_diagnostics", action="store_true", default=True)
    ap.add_argument("--quiet-diagnostics", dest="verbose_diagnostics", action="store_false")
    ap.add_argument("--registry-output", type=Path, default=None)
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    if args.restarts <= 0:
        raise SystemExit("--restarts must be positive")
    if args.eval_batch_size < 0:
        raise SystemExit("--eval-batch-size must be non-negative")
    if args.train_eval_max_rows < 0:
        raise SystemExit("--train-eval-max-rows must be non-negative")
    if not (0.0 < args.selection_target_signal_efficiency < 1.0):
        raise SystemExit("--selection-target-signal-efficiency must be between 0 and 1")
    if args.max_rows_per_pt_bin_class < 0:
        raise SystemExit("--max-rows-per-pt-bin-class must be non-negative")
    if not (0.0 <= float(args.distillation_strength) <= 0.90):
        raise SystemExit("--distillation-strength must be in [0, 0.90]")
    if not (math.isfinite(float(args.distillation_temperature)) and float(args.distillation_temperature) > 0.0):
        raise SystemExit("--distillation-temperature must be finite and positive")
    if not (0.0 <= float(args.distillation_min_finite_fraction) <= 1.0):
        raise SystemExit("--distillation-min-finite-fraction must be in [0, 1]")
    pt_bins = train_pt_bins(args)
    if args.centrality_range:
        parse_range(args.centrality_range)
    if args.pt_bin_weight_mode != "none" and not pt_bins:
        raise SystemExit("--pt-bin-weight-mode requires --train-pt-bins or --pt-range")
    if args.pt_bin_weight_mode != "none":
        parse_bin_weights(args.pt_bin_weight_spec, pt_bins, args.pt_bin_weight_mode)
    if args.selection_metric in ("highpt_wp80", "worst_pt_wp80", "highpt_auc"):
        parse_bin_weights(args.highpt_selection_weights, pt_bins, "highpt")
    products = parse_products(args.products)
    required_columns = set([args.label_branch, "cluster_Et", "cluster_Eta", "centrality", "vertexz"])
    required_columns.update(["run", "evt"])
    if args.hard_example_branch:
        required_columns.add(args.hard_example_branch)
    if args.distillation_branch:
        required_columns.add(args.distillation_branch)
    for product in products:
        required_columns.update(expand_required_columns(MODEL_SPECS[product]["features"]))
    optional_columns = [args.weight_branch]
    args.outdir.mkdir(parents=True, exist_ok=True)
    all_paths = expand_input_paths(args.input, args.source, args.check_input_files)
    paths, manifest_report = limit_paths_by_group(all_paths, args.max_files_per_sample, args.random_seed)
    manifest_report.update(
        {
            "schema": "RJ_AUAU_TIGHT_MLP_TRAINING_MANIFEST_SUMMARY_V1",
            "source": str(args.source) if args.source else None,
            "tree": args.tree,
            "products": products,
            "product_specs": {product: MODEL_SPECS[product] for product in products},
            "pt_range": args.pt_range,
            "centrality_range": args.centrality_range,
            "train_pt_bins": [{"lo": lo, "hi": hi} for lo, hi in train_pt_bins(args)],
            "pt_bin_weight_mode": args.pt_bin_weight_mode,
            "distillation_branch": args.distillation_branch,
            "distillation_strength": float(args.distillation_strength),
            "distillation_temperature": float(args.distillation_temperature),
            "require_distillation": bool(args.require_distillation),
            "random_seed": int(args.random_seed),
        }
    )
    write_json(args.outdir / "training_manifest_summary.json", manifest_report)
    print(
        "[trainAuAuPhotonMLP] manifest "
        f"files_before={manifest_report['input_files_before']['total_files']} "
        f"files_after={manifest_report['input_files_after']['total_files']} "
        f"groups={manifest_report['input_files_after']['groups']} "
        f"summary={args.outdir / 'training_manifest_summary.json'}",
        flush=True,
    )
    frame, optional_seen, read_report = load_frame(paths, args.tree, sorted(required_columns), optional_columns, args=args)
    frame = add_derived_features(frame)
    frame = compact_numeric_frame(frame, args.label_branch)
    write_json(args.outdir / "training_read_summary.json", {
        "schema": "RJ_AUAU_TIGHT_MLP_TRAINING_READ_SUMMARY_V1",
        "read": read_report,
        "optional_branches_seen": optional_seen,
    })
    print(
        "[trainAuAuPhotonMLP] loaded "
        f"files={read_report['files_read']} rows={read_report['rows_read']} "
        f"rows_after_read_filter={read_report['rows_after_read_filter']} "
        f"optional={optional_seen}",
        flush=True,
    )
    reports = []
    for product in products:
        reports.append(train_product(product, frame, args, args.outdir))
    registry = {
        "schema": "RJ_AUAU_TIGHT_MLP_REGISTRY_V1",
        "status": "READY",
        "source": str(args.source) if args.source else None,
        "input_files": [str(path) for path in paths],
        "tree": args.tree,
        "pt_range": args.pt_range,
        "products": products,
        "models": reports,
        "optional_branches_seen": optional_seen,
        "python": sys.version,
    }
    registry_path = args.registry_output or (args.outdir / "model_registry.json")
    registry_path.write_text(json.dumps(json_ready(registry), indent=2, sort_keys=True) + "\n")
    print(f"[trainAuAuPhotonMLP] wrote registry: {registry_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
