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

DERIVED_FEATURE_DEPS = {
    "cluster_weta_over_wphi": ["cluster_weta_cogx", "cluster_wphi_cogx"],
    "cluster_weta33_over_wphi33": ["cluster_weta33_cogx", "cluster_wphi33_cogx"],
}

MODEL_SPECS = {
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


def parse_products(text: str) -> list[str]:
    if text.strip().lower() in ("all", "default"):
        return list(DEFAULT_PRODUCTS)
    if text.strip().lower() in ("primary", "primary-only"):
        return ["centInputBase3x3MLP_pt1535"]
    if text.strip().lower() in ("primary-ratios", "primary-width-ratios", "ratios", "width-ratios"):
        return ["centInputBase3x3WidthRatiosMLP_pt1535"]
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
        frame["cluster_weta_over_wphi"] = safe_ratio("cluster_weta_cogx", "cluster_wphi_cogx")
        cols.add("cluster_weta_over_wphi")
    if "cluster_weta33_over_wphi33" not in cols and {"cluster_weta33_cogx", "cluster_wphi33_cogx"}.issubset(cols):
        frame["cluster_weta33_over_wphi33"] = safe_ratio("cluster_weta33_cogx", "cluster_wphi33_cogx")
        cols.add("cluster_weta33_over_wphi33")
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


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(json_ready(payload), indent=2, sort_keys=True) + "\n")


def load_frame(paths: list[Path], tree_name: str, required_columns: list[str], optional_columns: list[str] | None = None):
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
        "groups": {},
    }
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
            frames.append(chunk)
            group = path_group_key(path)
            group_report = read_report["groups"].setdefault(group, {"files": 0, "rows": 0})
            group_report["files"] += 1
            group_report["rows"] += int(len(chunk))
            read_report["files_read"] += 1
            read_report["rows_read"] += int(len(chunk))
    frame = pd.concat(frames, ignore_index=True)
    for col in optional_columns:
        if col not in frame.columns:
            frame[col] = 1.0
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
    out = frame.loc[mask].copy()
    before_cap = class_counts(out, label)
    out = balanced_row_cap(out, label, args.max_rows, args.max_rows_per_class, args.random_seed)
    after_cap = class_counts(out, label)
    if before_cap != after_cap:
        print(
            "[trainAuAuPhotonMLP] row cap "
            f"before={before_cap} after={after_cap} "
            f"max_rows={args.max_rows} max_rows_per_class={args.max_rows_per_class}",
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
    return z.astype("float64"), mean.astype("float64"), scale.astype("float64")


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


def train_numpy_mlp(x_train, y_train, w_train, x_val, y_val, w_val, feature_names, args, seed: int, hidden: list[int], candidate_label: str):
    import numpy as np

    params = init_params(x_train.shape[1], hidden, seed)
    moments = [{"mW": np.zeros_like(p["W"]), "vW": np.zeros_like(p["W"]), "mb": np.zeros_like(p["b"]), "vb": np.zeros_like(p["b"])} for p in params]
    rng = np.random.default_rng(seed + 17)
    jitter_features = [feature_names.index(f) for f in ("cluster_Et", "centrality") if f in feature_names]
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

            logits, cache = forward(xb, params)
            dz = ((sigmoid(logits) - yb) * wb_norm).reshape(-1, 1)
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

        train_logits, _ = forward(x_train, params)
        val_logits, _ = forward(x_val, params)
        train_loss = weighted_bce_from_logits(train_logits, y_train, w_train, params, args.l2)
        val_loss = weighted_bce_from_logits(val_logits, y_val, w_val, params, args.l2)
        row = {
            "epoch": epoch,
            "train_loss": train_loss,
            "validation_loss": val_loss,
            "validation_auc": auc_score(y_val, sigmoid(val_logits), w_val),
        }
        history.append(row)
        if epoch == 1 or epoch % max(1, args.progress_every) == 0:
            print(
                "[trainAuAuPhotonMLP] "
                f"candidate={candidate_label} epoch={epoch} train_loss={train_loss:.6g} "
                f"val_loss={val_loss:.6g} val_auc={row['validation_auc']:.5f}",
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

    x_raw = np.column_stack([selected[f].to_numpy(dtype="float64") for f in features])
    x, mean, scale = standardize_from_train(x_raw, train_mask, args.input_clip)
    hidden_grid = parse_hidden_layer_grid(args.hidden_layers, args.hidden_layer_grid)
    target_eff = float(args.selection_target_signal_efficiency)
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
            )
            val_logits_i, _ = forward(x[val_mask], params_i)
            temperature_i, temp_loss_i = fit_temperature(val_logits_i, labels[val_mask], weights[val_mask])
            val_prob_i = sigmoid(val_logits_i / temperature_i)
            val_report_i = make_report(labels[val_mask], val_prob_i, weights[val_mask])
            selection_wp = threshold_for_signal_efficiency(labels[val_mask], val_prob_i, target_eff) or {}
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
    train_logits, _ = forward(x[train_mask], params)
    val_logits, _ = forward(x[val_mask], params)
    test_logits, _ = forward(x[test_mask], params)
    train_prob = sigmoid(train_logits / temperature)
    val_prob = sigmoid(val_logits / temperature)
    test_prob = sigmoid(test_logits / temperature)
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
        "temperature": {"value": float(temperature), "validation_loss": float(temp_loss)},
        "model_selection": {
            "metric": args.selection_metric,
            "target_signal_efficiency": float(target_eff),
            "hidden_layer_grid": [[int(x) for x in hidden] for hidden in hidden_grid],
            "restarts": int(args.restarts),
            "selected_candidate": selected_candidate,
            "candidates": candidate_reports,
        },
        "train": make_report(labels[train_mask], train_prob, weights[train_mask]),
        "validation": make_report(labels[val_mask], val_prob, weights[val_mask]),
        "test": make_report(labels[test_mask], test_prob, weights[test_mask]),
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
    ap.add_argument("--validation-fraction", type=float, default=0.20)
    ap.add_argument("--test-fraction", type=float, default=0.10)
    ap.add_argument("--random-seed", type=int, default=137)
    ap.add_argument("--hidden-layers", default="48,24")
    ap.add_argument("--hidden-layer-grid", default=os.environ.get("RJ_AUAU_MLP_TRAIN_HIDDEN_LAYER_GRID", ""))
    ap.add_argument("--restarts", type=int, default=int(os.environ.get("RJ_AUAU_MLP_TRAIN_RESTARTS", "1")))
    ap.add_argument("--selection-metric", choices=["wp_fake_rate", "validation_auc", "validation_loss"], default=os.environ.get("RJ_AUAU_MLP_TRAIN_SELECTION_METRIC", "wp_fake_rate"))
    ap.add_argument("--selection-target-signal-efficiency", type=float, default=float(os.environ.get("RJ_AUAU_MLP_TRAIN_SELECTION_TARGET_SIGNAL_EFF", "0.80")))
    ap.add_argument("--epochs", type=int, default=80)
    ap.add_argument("--batch-size", type=int, default=4096)
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
    if not (0.0 < args.selection_target_signal_efficiency < 1.0):
        raise SystemExit("--selection-target-signal-efficiency must be between 0 and 1")
    products = parse_products(args.products)
    required_columns = set([args.label_branch, "cluster_Et", "cluster_Eta", "centrality", "vertexz"])
    required_columns.update(["run", "evt"])
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
            "pt_range": args.pt_range,
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
    frame, optional_seen, read_report = load_frame(paths, args.tree, sorted(required_columns), optional_columns)
    frame = add_derived_features(frame)
    write_json(args.outdir / "training_read_summary.json", {
        "schema": "RJ_AUAU_TIGHT_MLP_TRAINING_READ_SUMMARY_V1",
        "read": read_report,
        "optional_branches_seen": optional_seen,
    })
    print(
        "[trainAuAuPhotonMLP] loaded "
        f"files={read_report['files_read']} rows={read_report['rows_read']} "
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
