#!/usr/bin/env python3
"""Validate typed Au+Au SUB1 isolation replay from THE-119 ROOT rows.

The runtime witness follows PhotonClusterBuilder's standard-SUB1 convention:
sum accepted SUB1 tower ET in the requested cone across EMCal, HCALIN, and
HCALOUT, then subtract the reconstructed photon ET once.  The normalized
constituent rows must reproduce that value independently for R=0.3 and R=0.4.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import uproot


TREE_NAMES = {
    "candidate": "RJPhotonCandidateV1",
    "constituent": "RJIsolationConstituentV1",
    "witness": "RJIsolationWitnessV1",
}
RADII = (0.3, 0.4)


def _id_key(hi: Any, lo: Any) -> tuple[int, int]:
    return int(hi), int(lo)


def _grouped_sum(
    ids_hi: np.ndarray,
    ids_lo: np.ndarray,
    values: np.ndarray,
) -> list[tuple[tuple[int, int], float, int]]:
    """Return stable-ID grouped sums and counts for one uproot chunk."""
    if len(values) == 0:
        return []
    keys = np.empty(len(values), dtype=[("hi", "<u8"), ("lo", "<u8")])
    keys["hi"] = ids_hi
    keys["lo"] = ids_lo
    unique, inverse = np.unique(keys, return_inverse=True)
    sums = np.bincount(inverse, weights=values.astype(np.float64, copy=False))
    counts = np.bincount(inverse)
    return [
        (_id_key(key["hi"], key["lo"]), float(sums[index]), int(counts[index]))
        for index, key in enumerate(unique)
    ]


def validate_file(
    path: Path, absolute_tolerance: float, minimum_bytes: int
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "path": str(path),
        "size_bytes": path.stat().st_size if path.exists() else 0,
        "absolute_tolerance_gev": absolute_tolerance,
        "status": "FAIL",
        "failures": [],
        "radii": {},
    }
    if not path.exists() or path.stat().st_size < minimum_bytes:
        result["failures"].append("missing_or_tiny_root")
        return result

    try:
        root = uproot.open(path)
    except Exception as exc:  # pragma: no cover - evidence path
        result["failures"].append(f"uproot_open_failed:{exc}")
        return result

    with root:
        replay = root["ReplayFoundationV1"] if "ReplayFoundationV1" in root else root
        key_names = {key.split(";")[0] for key in replay.keys()}
        missing_trees = sorted(set(TREE_NAMES.values()) - key_names)
        if missing_trees:
            result["failures"].append(f"missing_trees:{','.join(missing_trees)}")
            return result

        candidate_arrays = replay[TREE_NAMES["candidate"]].arrays(
            ["candidate_id_hi", "candidate_id_lo", "cluster_et"], library="np"
        )
        candidate_et: dict[tuple[int, int], float] = {}
        duplicate_candidates = 0
        for hi, lo, et in zip(
            candidate_arrays["candidate_id_hi"],
            candidate_arrays["candidate_id_lo"],
            candidate_arrays["cluster_et"],
            strict=True,
        ):
            key = _id_key(hi, lo)
            if key in candidate_et:
                duplicate_candidates += 1
            candidate_et[key] = float(et)
        result["candidate_count"] = len(candidate_et)
        result["duplicate_candidate_ids"] = duplicate_candidates
        if duplicate_candidates:
            result["failures"].append("duplicate_candidate_ids")

        cone_sums = {radius: {} for radius in RADII}
        cone_counts = {radius: {} for radius in RADII}
        invalid_quality_rows = 0
        accepted_rows = 0
        nonfinite_accepted_rows = 0
        branches = [
            "candidate_id_hi",
            "candidate_id_lo",
            "delta_r",
            "sub1_energy",
            "quality_state",
        ]
        for arrays in replay[TREE_NAMES["constituent"]].iterate(
            branches, step_size="128 MB", library="np"
        ):
            quality = arrays["quality_state"]
            invalid_quality_rows += int(np.count_nonzero(quality < 0))
            accepted = quality == 1
            finite = np.isfinite(arrays["delta_r"]) & np.isfinite(arrays["sub1_energy"])
            nonfinite_accepted_rows += int(np.count_nonzero(accepted & ~finite))
            base = accepted & finite
            accepted_rows += int(np.count_nonzero(base))
            for radius in RADII:
                mask = base & (arrays["delta_r"] < radius)
                grouped = _grouped_sum(
                    arrays["candidate_id_hi"][mask],
                    arrays["candidate_id_lo"][mask],
                    arrays["sub1_energy"][mask],
                )
                for key, value, count in grouped:
                    cone_sums[radius][key] = cone_sums[radius].get(key, 0.0) + value
                    cone_counts[radius][key] = cone_counts[radius].get(key, 0) + count
        result["accepted_constituent_rows"] = accepted_rows
        result["invalid_quality_rows"] = invalid_quality_rows
        result["nonfinite_accepted_rows"] = nonfinite_accepted_rows
        if nonfinite_accepted_rows:
            result["failures"].append("nonfinite_accepted_constituents")

        witness_arrays = replay[TREE_NAMES["witness"]].arrays(
            [
                "candidate_id_hi",
                "candidate_id_lo",
                "radius",
                "subtraction_method",
                "reconstructed_or_truth",
                "cone_sum",
            ],
            library="np",
        )
        witness_rows = 0
        duplicate_witnesses = 0
        seen_witnesses: set[tuple[tuple[int, int], float]] = set()
        diffs = {radius: [] for radius in RADII}
        missing_candidates = {radius: 0 for radius in RADII}
        for hi, lo, radius_value, method, reco_truth, witness in zip(
            witness_arrays["candidate_id_hi"],
            witness_arrays["candidate_id_lo"],
            witness_arrays["radius"],
            witness_arrays["subtraction_method"],
            witness_arrays["reconstructed_or_truth"],
            witness_arrays["cone_sum"],
            strict=True,
        ):
            if int(method) != 2 or int(reco_truth) != 0:
                continue
            radius = min(RADII, key=lambda value: abs(value - float(radius_value)))
            if abs(radius - float(radius_value)) > 1.0e-9:
                result["failures"].append(f"unexpected_radius:{float(radius_value)}")
                continue
            key = _id_key(hi, lo)
            witness_key = (key, radius)
            if witness_key in seen_witnesses:
                duplicate_witnesses += 1
                continue
            seen_witnesses.add(witness_key)
            witness_rows += 1
            if key not in candidate_et:
                missing_candidates[radius] += 1
                continue
            replay_value = cone_sums[radius].get(key, 0.0) - candidate_et[key]
            diffs[radius].append(replay_value - float(witness))

        result["typed_witness_rows"] = witness_rows
        result["duplicate_typed_witnesses"] = duplicate_witnesses
        if duplicate_witnesses:
            result["failures"].append("duplicate_typed_witnesses")

        for radius in RADII:
            values = np.asarray(diffs[radius], dtype=np.float64)
            abs_values = np.abs(values)
            radius_result = {
                "witness_count": int(len(values)),
                "missing_candidate_foreign_keys": missing_candidates[radius],
                "candidate_count_with_constituents": len(cone_sums[radius]),
                "accepted_constituent_count": int(sum(cone_counts[radius].values())),
                "max_abs_difference_gev": float(np.max(abs_values)) if len(values) else None,
                "mean_abs_difference_gev": float(np.mean(abs_values)) if len(values) else None,
                "rms_difference_gev": float(np.sqrt(np.mean(values * values))) if len(values) else None,
                "outside_tolerance": int(np.count_nonzero(abs_values > absolute_tolerance)),
            }
            result["radii"][f"R{int(round(radius * 10)):02d}"] = radius_result
            if missing_candidates[radius]:
                result["failures"].append(f"R{radius}:missing_candidate_foreign_keys")
            if len(values) != len(candidate_et):
                result["failures"].append(
                    f"R{radius}:witness_population_{len(values)}_expected_{len(candidate_et)}"
                )
            if radius_result["outside_tolerance"]:
                result["failures"].append(
                    f"R{radius}:isolation_replay_outside_tolerance"
                )

    result["failures"] = sorted(set(result["failures"]))
    result["status"] = "PASS" if not result["failures"] else "FAIL"
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root_files", nargs="+", type=Path)
    parser.add_argument("--absolute-tolerance", type=float, default=1.0e-4)
    parser.add_argument("--minimum-bytes", type=int, default=50_000)
    parser.add_argument("--output-json", type=Path)
    args = parser.parse_args()
    if not math.isfinite(args.absolute_tolerance) or args.absolute_tolerance < 0:
        parser.error("--absolute-tolerance must be finite and nonnegative")
    if args.minimum_bytes < 0:
        parser.error("--minimum-bytes must be nonnegative")

    report = {
        "schema": "THE119_AUAU_ISOLATION_REPLAY_VALIDATION_V1",
        "formula": "sum(good finite per-tower SUB1 ET within R across three calorimeters) - cluster_et",
        "files": [
            validate_file(path.resolve(), args.absolute_tolerance, args.minimum_bytes)
            for path in args.root_files
        ],
    }
    report["status"] = (
        "PASS" if report["files"] and all(item["status"] == "PASS" for item in report["files"]) else "FAIL"
    )
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output_json:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(payload)
    print(payload, end="")
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
