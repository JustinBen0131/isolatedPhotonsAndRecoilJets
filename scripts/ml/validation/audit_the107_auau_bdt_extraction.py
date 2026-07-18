#!/usr/bin/env python3
"""Audit the fresh THE-107 Au+Au photon-ID extraction before training.

The audit proves that truth class, truth isolation, nominal labels, and source
roles remain independently recoverable in every extraction file.  It reads
only the compact label/provenance branches and is safe to parallelize across
the six embedded source samples.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np


TREE_NAME = "AuAuPhotonIDTrainingTree"
SAMPLES = {
    "run28_embeddedPhoton12": (1, 12),
    "run28_embeddedPhoton20": (1, 20),
    "run28_embeddedJet12": (2, 12),
    "run28_embeddedJet20": (2, 20),
    "run28_embeddedJet30": (2, 30),
    "run28_embeddedJet40": (2, 40),
}
BRANCHES = [
    "is_signal",
    "truth_match_found",
    "truth_photon_class",
    "truth_is_prompt",
    "truth_iso_et",
    "truth_iso_pass",
    "cluster_truth_track_id",
    "cluster_truth_pid",
    "cluster_truth_barcode",
    "source_role",
    "source_sample_code",
    "ppg12_source_role_label",
    "minimum_bias_classifier_decision",
    "cluster_Et",
    "cluster_Eta",
    "centrality",
    "vertexz",
]


def sample_from_path(path: Path) -> str | None:
    text = str(path)
    matches = [sample for sample in SAMPLES if re.search(rf"(?:^|/){re.escape(sample)}(?:/|$)", text)]
    if len(matches) != 1:
        return None
    return matches[0]


def manifest_digest(paths: list[Path], base: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(paths):
        stat = path.stat()
        record = f"{path.relative_to(base)}\t{stat.st_size}\t{stat.st_mtime_ns}\n"
        digest.update(record.encode("utf-8"))
    return digest.hexdigest()


def audit_file(item: tuple[str, str, float]) -> dict:
    import uproot

    raw_path, sample, truth_iso_max = item
    path = Path(raw_path)
    role, sample_code = SAMPLES[sample]
    result = {
        "path": raw_path,
        "sample": sample,
        "rows": 0,
        "class_counts": {str(code): 0 for code in (0, 1, 2, 3)},
        "nominal_signal": 0,
        "nominal_background": 0,
        "ppg12_signal": 0,
        "ppg12_background": 0,
        "ppg12_discarded": 0,
        "mismatches": {},
        "error": None,
    }
    try:
        with uproot.open(path) as root_file:
            if TREE_NAME not in root_file:
                raise RuntimeError(f"missing tree {TREE_NAME}")
            tree = root_file[TREE_NAME]
            missing = [name for name in BRANCHES if name not in tree.keys()]
            if missing:
                raise RuntimeError("missing branches: " + ", ".join(missing))
            arrays = tree.arrays(BRANCHES, library="np", how=dict)
    except Exception as exc:  # noqa: BLE001
        result["error"] = str(exc)
        return result

    n_rows = len(arrays["is_signal"])
    result["rows"] = int(n_rows)
    if n_rows == 0:
        return result

    def values(name: str, dtype: str) -> np.ndarray:
        return np.asarray(arrays[name], dtype=dtype)

    label = values("is_signal", "int32")
    match = values("truth_match_found", "int32")
    photon_class = values("truth_photon_class", "int32")
    prompt = values("truth_is_prompt", "int32")
    iso_et = values("truth_iso_et", "float64")
    iso_pass = values("truth_iso_pass", "int32")
    source_role = values("source_role", "int32")
    source_code = values("source_sample_code", "int32")
    ppg12 = values("ppg12_source_role_label", "int32")
    mb_decision = values("minimum_bias_classifier_decision", "int32")

    prompt_expected = np.isin(photon_class, [1, 2])
    iso_expected = np.zeros(n_rows, dtype=bool)
    iso_expected[prompt_expected] = iso_et[prompt_expected] < float(truth_iso_max)
    nominal_expected = prompt_expected & iso_expected
    ppg12_expected = np.full(n_rows, -1, dtype="int32")
    ppg12_expected[(role == 1) & prompt_expected] = 1
    ppg12_expected[(role == 2) & ~prompt_expected] = 0

    mismatch_masks = {
        "invalid_is_signal": ~np.isin(label, [0, 1]),
        "invalid_truth_match_found": ~np.isin(match, [0, 1]),
        "invalid_truth_photon_class": ~np.isin(photon_class, [0, 1, 2, 3]),
        "truth_is_prompt_vs_class": prompt != prompt_expected.astype("int32"),
        "invalid_truth_iso_pass": ~np.isin(iso_pass, [0, 1]),
        "prompt_truth_iso_nonfinite": prompt_expected & ~np.isfinite(iso_et),
        "truth_iso_pass_vs_et": iso_pass != iso_expected.astype("int32"),
        "nominal_label_closure": label != nominal_expected.astype("int32"),
        "source_role": source_role != role,
        "source_sample_code": source_code != sample_code,
        "ppg12_label_closure": ppg12 != ppg12_expected,
        "minimum_bias_classifier_decision": mb_decision != 2,
    }
    result["mismatches"] = {
        name: int(mask.sum()) for name, mask in mismatch_masks.items() if int(mask.sum())
    }
    result["class_counts"] = {
        str(code): int(np.sum(photon_class == code)) for code in (0, 1, 2, 3)
    }
    result["nominal_signal"] = int(np.sum(label == 1))
    result["nominal_background"] = int(np.sum(label == 0))
    result["ppg12_signal"] = int(np.sum(ppg12 == 1))
    result["ppg12_background"] = int(np.sum(ppg12 == 0))
    result["ppg12_discarded"] = int(np.sum(ppg12 < 0))
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True, help="Extraction root containing the six samples")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--samples",
        default=",".join(SAMPLES),
        help=(
            "Comma-separated source samples to require and audit. The default "
            "retains the full six-sample production gate."
        ),
    )
    parser.add_argument("--expected-files-per-sample", type=int, default=715)
    parser.add_argument("--min-file-size", type=int, default=10_240)
    parser.add_argument("--truth-iso-max", type=float, default=4.0)
    parser.add_argument("--workers", type=int, default=12)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    selected_samples = [item.strip() for item in args.samples.split(",") if item.strip()]
    unknown_samples = sorted(set(selected_samples) - set(SAMPLES))
    if not selected_samples or unknown_samples or len(selected_samples) != len(set(selected_samples)):
        raise SystemExit(
            "invalid --samples selection: "
            f"selected={selected_samples} unknown={unknown_samples}; choices={list(SAMPLES)}"
        )
    root = args.root.resolve()
    if not root.is_dir():
        raise SystemExit(f"missing extraction root: {root}")
    all_roots = sorted(root.rglob("*.root"))
    too_small = [path for path in all_roots if path.stat().st_size < args.min_file_size]
    classified: dict[str, list[Path]] = {sample: [] for sample in selected_samples}
    unknown = []
    for path in all_roots:
        sample = sample_from_path(path)
        if sample is None or sample not in classified:
            unknown.append(path)
        else:
            classified[sample].append(path)

    file_count_failures = {
        sample: len(paths)
        for sample, paths in classified.items()
        if len(paths) != args.expected_files_per_sample
    }
    if unknown or too_small or file_count_failures:
        report = {
            "schema": "THE107_AUAU_BDT_EXTRACTION_AUDIT_V1",
            "status": "FAILED_FILE_INVENTORY",
            "root": str(root),
            "selected_samples": selected_samples,
            "expected_files_per_sample": args.expected_files_per_sample,
            "files_per_sample": {sample: len(paths) for sample, paths in classified.items()},
            "unknown_files": [str(path) for path in unknown],
            "files_below_minimum_size": [str(path) for path in too_small],
            "file_count_failures": file_count_failures,
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        raise SystemExit(f"extraction file inventory failed; see {args.output}")

    items = [
        (str(path), sample, float(args.truth_iso_max))
        for sample, paths in classified.items()
        for path in paths
    ]
    results = []
    with ProcessPoolExecutor(max_workers=max(1, int(args.workers))) as executor:
        futures = {executor.submit(audit_file, item): item[0] for item in items}
        for future in as_completed(futures):
            try:
                results.append(future.result())
            except Exception as exc:  # noqa: BLE001
                results.append(
                    {
                        "path": futures[future],
                        "sample": sample_from_path(Path(futures[future])) or "unknown",
                        "rows": 0,
                        "class_counts": {},
                        "nominal_signal": 0,
                        "nominal_background": 0,
                        "ppg12_signal": 0,
                        "ppg12_background": 0,
                        "ppg12_discarded": 0,
                        "mismatches": {},
                        "error": repr(exc),
                    }
                )

    by_sample = {}
    for sample in selected_samples:
        subset = [item for item in results if item["sample"] == sample]
        aggregate = {
            "files": len(subset),
            "rows": int(sum(item["rows"] for item in subset)),
            "class_counts": {
                str(code): int(sum(item.get("class_counts", {}).get(str(code), 0) for item in subset))
                for code in (0, 1, 2, 3)
            },
            "nominal_signal": int(sum(item["nominal_signal"] for item in subset)),
            "nominal_background": int(sum(item["nominal_background"] for item in subset)),
            "ppg12_signal": int(sum(item["ppg12_signal"] for item in subset)),
            "ppg12_background": int(sum(item["ppg12_background"] for item in subset)),
            "ppg12_discarded": int(sum(item["ppg12_discarded"] for item in subset)),
        }
        by_sample[sample] = aggregate

    bad_files = [item for item in results if item.get("error") or item.get("mismatches")]
    report = {
        "schema": "THE107_AUAU_BDT_EXTRACTION_AUDIT_V1",
        "status": "PASSED" if not bad_files else "FAILED_CONTENT_AUDIT",
        "root": str(root),
        "selected_samples": selected_samples,
        "tree": TREE_NAME,
        "truth_iso_max_gev": float(args.truth_iso_max),
        "expected_files_per_sample": int(args.expected_files_per_sample),
        "total_files": len(all_roots),
        "total_rows": int(sum(item["rows"] for item in results)),
        "fileset_manifest_sha256": manifest_digest(all_roots, root),
        "samples": by_sample,
        "bad_file_count": len(bad_files),
        "bad_files": bad_files,
        "checks": [
            "every file opens and contains the required tree/branches",
            "truth class and prompt flag agree",
            "truth-isolation flag equals E_T^iso < 4 GeV for prompt photons",
            "nominal is_signal equals prompt AND truth-isolated",
            "PPG12 source-role label equals photon-source prompt or jet-source non-prompt",
            "source role/sample code match path provenance",
            "MinimumBiasClassifier decision equals 2",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": report["status"], "files": len(all_roots), "rows": report["total_rows"], "output": str(args.output)}, sort_keys=True))
    return 0 if report["status"] == "PASSED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
