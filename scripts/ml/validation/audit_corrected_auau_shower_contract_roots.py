#!/usr/bin/env python3
"""Verify the corrected Au+Au shower-feature contract in every extraction ROOT."""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path
import re


TREE = "AuAuPhotonIDTrainingTree"
SAMPLES = (
    "run28_embeddedPhoton12",
    "run28_embeddedPhoton20",
    "run28_embeddedJet12",
    "run28_embeddedJet20",
    "run28_embeddedJet30",
    "run28_embeddedJet40",
)
FEATURES = (
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
    "is_signal",
)
EXPECTED_CONTRACT = {
    "cemc_shower_shape_energy_source": "towerinfo_full_good_grid",
    "cemc_shower_shape_raw_cluster_towermap": False,
    "cemc_shower_shape_tower_acceptance": "towerinfo_get_isgood",
    "cemc_shower_shape_tower_min_energy_gev": 0,
    "cemc_shower_shape_diagnostic_variant": "canonical",
}


def sample_from_path(path: Path) -> str | None:
    matches = [
        sample
        for sample in SAMPLES
        if re.search(rf"(?:^|/){re.escape(sample)}(?:/|$)", str(path))
    ]
    return matches[0] if len(matches) == 1 else None


def audit_one(raw_path: str) -> dict:
    import uproot
    import yaml

    path = Path(raw_path)
    result = {
        "path": raw_path,
        "sample": sample_from_path(path),
        "rows": 0,
        "contract": {},
        "error": None,
    }
    try:
        with uproot.open(path) as root_file:
            if TREE not in root_file:
                raise RuntimeError(f"missing tree {TREE}")
            tree = root_file[TREE]
            missing = [name for name in FEATURES if name not in tree.keys()]
            if missing:
                raise RuntimeError("missing feature/label branches: " + ", ".join(missing))
            if "analysis_config_yaml" not in root_file:
                raise RuntimeError("missing analysis_config_yaml")
            config = yaml.safe_load(str(root_file["analysis_config_yaml"]))
            if not isinstance(config, dict):
                raise RuntimeError("analysis_config_yaml did not decode to a mapping")
            observed = {key: config.get(key) for key in EXPECTED_CONTRACT}
            result["rows"] = int(tree.num_entries)
            result["contract"] = observed
            mismatches = {
                key: {"expected": expected, "observed": observed.get(key)}
                for key, expected in EXPECTED_CONTRACT.items()
                if observed.get(key) != expected
            }
            if mismatches:
                raise RuntimeError("shower-contract mismatch: " + json.dumps(mismatches))
    except Exception as exc:  # noqa: BLE001
        result["error"] = str(exc)
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-files-per-sample", type=int, default=715)
    parser.add_argument("--min-file-size", type=int, default=10_240)
    parser.add_argument("--workers", type=int, default=12)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.root.resolve()
    if not root.is_dir():
        raise SystemExit(f"missing extraction root: {root}")
    paths = sorted(root.rglob("*.root"))
    classified = Counter(sample_from_path(path) for path in paths)
    inventory_errors = {
        "unknown_sample_files": int(classified.get(None, 0)),
        "files_below_minimum_size": [
            str(path) for path in paths if path.stat().st_size < args.min_file_size
        ],
        "sample_file_count_mismatches": {
            sample: int(classified.get(sample, 0))
            for sample in SAMPLES
            if classified.get(sample, 0) != args.expected_files_per_sample
        },
    }
    inventory_ok = not inventory_errors["unknown_sample_files"] and not inventory_errors[
        "files_below_minimum_size"
    ] and not inventory_errors["sample_file_count_mismatches"]

    results: list[dict] = []
    if inventory_ok:
        with ProcessPoolExecutor(max_workers=max(1, int(args.workers))) as executor:
            futures = {executor.submit(audit_one, str(path)): path for path in paths}
            for future in as_completed(futures):
                try:
                    results.append(future.result())
                except Exception as exc:  # noqa: BLE001
                    results.append(
                        {
                            "path": str(futures[future]),
                            "sample": sample_from_path(futures[future]),
                            "rows": 0,
                            "contract": {},
                            "error": repr(exc),
                        }
                    )

    bad = [item for item in results if item.get("error")]
    rows_by_sample = {
        sample: int(
            sum(item["rows"] for item in results if item.get("sample") == sample)
        )
        for sample in SAMPLES
    }
    if not inventory_ok:
        status = "FAILED_FILE_INVENTORY"
    elif bad:
        status = "FAILED_CONTRACT_AUDIT"
    else:
        status = "PASSED"
    report = {
        "schema": "CORRECTED_AUAU_SHOWER_CONTRACT_ROOT_AUDIT_V1",
        "status": status,
        "root": str(root),
        "tree": TREE,
        "expected_contract": EXPECTED_CONTRACT,
        "required_feature_order": list(FEATURES),
        "expected_files_per_sample": int(args.expected_files_per_sample),
        "total_files": len(paths),
        "files_per_sample": {sample: int(classified.get(sample, 0)) for sample in SAMPLES},
        "rows_per_sample": rows_by_sample,
        "total_rows": int(sum(rows_by_sample.values())),
        "inventory_errors": inventory_errors,
        "bad_file_count": len(bad),
        "bad_files": bad,
        "checks": [
            "all six source samples contain the expected number of nontrivial ROOT files",
            "every ROOT contains the full ordered 14-feature vector and is_signal label",
            "every embedded ROOT uses calibrated full-grid TowerInfo shower inputs",
            "only TowerInfo::get_isGood() controls shower-cell acceptance",
            "RawCluster-owned cell routing is disabled for embedded shower shapes",
            "the Au+Au shower-cell floor is zero GeV and the diagnostic variant is canonical",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "status": status,
                "files": len(paths),
                "rows": report["total_rows"],
                "output": str(args.output),
            },
            sort_keys=True,
        )
    )
    return 0 if status == "PASSED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
