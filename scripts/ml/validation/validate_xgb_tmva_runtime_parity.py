#!/usr/bin/env python3
"""Validate exported TMVA/RBDT scores against saved XGBoost holdout scores.

The trainer writes the exact held-out feature matrix and XGBoost score into an
NPZ artifact.  This validator evaluates the exported ROOT model on those same
rows, records deterministic numerical differences, and fails when the maximum
absolute difference exceeds the requested tolerance.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


EXPECTED_FEATURES = [
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
    "centrality",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def decode_features(values: np.ndarray) -> list[str]:
    return [
        item.decode() if isinstance(item, (bytes, np.bytes_)) else str(item)
        for item in values.tolist()
    ]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--holdout", type=Path, required=True)
    parser.add_argument("--tmva-root", type=Path, required=True)
    parser.add_argument("--model-metadata", type=Path, required=True)
    parser.add_argument("--json-out", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=100_000)
    parser.add_argument("--max-rows", type=int, default=0)
    parser.add_argument("--tolerance", type=float, default=2.0e-6)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    for path in (args.holdout, args.tmva_root, args.model_metadata):
        if not path.is_file():
            raise SystemExit(f"missing required input: {path}")
    if args.chunk_size <= 0:
        raise SystemExit("--chunk-size must be positive")
    if args.max_rows < 0:
        raise SystemExit("--max-rows must be non-negative")
    if args.tolerance <= 0.0:
        raise SystemExit("--tolerance must be positive")

    metadata = json.loads(args.model_metadata.read_text())
    if metadata.get("features") != EXPECTED_FEATURES:
        raise SystemExit("model metadata does not contain the frozen 14-feature order")

    # The trusted trainer artifact stores only its feature-name vector with
    # object dtype; the numerical matrices remain ordinary NumPy arrays.  The
    # artifact hash is recorded below, so enabling object decoding here is
    # bounded to this explicitly supplied, provenance-pinned holdout file.
    holdout = np.load(args.holdout, allow_pickle=True)
    required = {"x", "score_xgboost", "features"}
    missing = sorted(required.difference(holdout.files))
    if missing:
        raise SystemExit(f"holdout NPZ is missing required keys: {missing}")
    features = decode_features(np.asarray(holdout["features"]))
    if features != EXPECTED_FEATURES:
        raise SystemExit(
            "holdout feature order does not match the frozen contract: "
            + json.dumps(features)
        )

    x = np.asarray(holdout["x"], dtype=np.float32)
    reference = np.asarray(holdout["score_xgboost"], dtype=np.float64).reshape(-1)
    if x.ndim != 2 or x.shape[1] != len(EXPECTED_FEATURES):
        raise SystemExit(f"unexpected holdout matrix shape: {x.shape}")
    if x.shape[0] != reference.size:
        raise SystemExit(
            f"holdout row mismatch: x={x.shape[0]} scores={reference.size}"
        )
    n_available = int(reference.size)
    n_rows = min(n_available, args.max_rows) if args.max_rows else n_available
    x = np.ascontiguousarray(x[:n_rows], dtype=np.float32)
    reference = reference[:n_rows]

    import ROOT  # pylint: disable=import-outside-toplevel

    model = ROOT.TMVA.Experimental.RBDT("myBDT", str(args.tmva_root))
    observed = np.empty(n_rows, dtype=np.float64)
    for start in range(0, n_rows, args.chunk_size):
        stop = min(n_rows, start + args.chunk_size)
        values = np.asarray(model.Compute(x[start:stop]), dtype=np.float64).reshape(-1)
        if values.size != stop - start:
            raise SystemExit(
                f"TMVA output-size mismatch for rows {start}:{stop}: {values.shape}"
            )
        observed[start:stop] = values

    finite = np.isfinite(reference) & np.isfinite(observed)
    difference = observed - reference
    absolute = np.abs(difference)
    max_abs = float(np.max(absolute[finite])) if finite.any() else float("nan")
    mean_abs = float(np.mean(absolute[finite])) if finite.any() else float("nan")
    rms = (
        float(np.sqrt(np.mean(np.square(difference[finite]))))
        if finite.any()
        else float("nan")
    )
    p99_abs = (
        float(np.quantile(absolute[finite], 0.99)) if finite.any() else float("nan")
    )
    n_nonfinite = int((~finite).sum())
    n_above = int(np.count_nonzero(finite & (absolute > args.tolerance)))
    passed = n_nonfinite == 0 and np.isfinite(max_abs) and max_abs <= args.tolerance

    payload = {
        "schema": "XGBOOST_TMVA_RUNTIME_PARITY_V1",
        "status": "PASS" if passed else "FAIL",
        "promotion_status": "NOT_PROMOTED",
        "method": "TMVA::Experimental::RBDT evaluated on saved XGBoost holdout rows",
        "rows_available": n_available,
        "rows_evaluated": n_rows,
        "feature_order": features,
        "tolerance": args.tolerance,
        "metrics": {
            "max_abs_difference": max_abs,
            "mean_abs_difference": mean_abs,
            "rms_difference": rms,
            "p99_abs_difference": p99_abs,
            "nonfinite_rows": n_nonfinite,
            "rows_above_tolerance": n_above,
        },
        "provenance": {
            "holdout": str(args.holdout),
            "holdout_sha256": sha256(args.holdout),
            "tmva_root": str(args.tmva_root),
            "tmva_root_sha256": sha256(args.tmva_root),
            "model_metadata": str(args.model_metadata),
            "model_metadata_sha256": sha256(args.model_metadata),
        },
    }
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))
    if not passed:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
