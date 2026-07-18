#!/usr/bin/env python3
"""Derive centrality-dependent Au+Au BDT working points with clean provenance.

This is a thin, identity-checking wrapper around the established weighted
working-point implementation.  It deliberately changes no numerical method;
it replaces stale campaign-facing labels and records the corrected shower
feature contract used by the retrained model.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
from types import ModuleType


EXPECTED_PRODUCT = "centAsFeatBase3x3_pt15to35"
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


def load_core() -> ModuleType:
    path = Path(__file__).with_name("derive_the57_full_weighted_wp80.py")
    spec = importlib.util.spec_from_file_location("weighted_wp_core", path)
    if spec is None or spec.loader is None:
        raise SystemExit(f"cannot load weighted-WP core: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_model_identity(metadata: dict) -> None:
    observed = {
        "product": metadata.get("product"),
        "features": metadata.get("features"),
        "pt_range": metadata.get("pt_range"),
        "label_branch": metadata.get("label_branch"),
        "weight_mode": metadata.get("weight_mode"),
    }
    expected = {
        "product": EXPECTED_PRODUCT,
        "features": EXPECTED_FEATURES,
        "pt_range": [15.0, 35.0],
        "label_branch": "is_signal",
        "weight_mode": "ppg12-exact",
    }
    if observed != expected:
        raise SystemExit(
            "model metadata does not match the frozen 14-feature contract:\n"
            + json.dumps({"expected": expected, "observed": observed}, indent=2)
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--model-metadata", type=Path, required=True)
    parser.add_argument("--extraction-audit", type=Path, required=True)
    parser.add_argument("--json-out", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=500_000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    for path in (args.matrix, args.model, args.model_metadata, args.extraction_audit):
        if not path.is_file():
            raise SystemExit(f"missing required input: {path}")

    metadata = json.loads(args.model_metadata.read_text())
    extraction_audit = json.loads(args.extraction_audit.read_text())
    require_model_identity(metadata)
    if extraction_audit.get("status") not in {"PASS", "PASSED"}:
        raise SystemExit(
            f"corrected extraction audit is not passing: {extraction_audit.get('status')}"
        )

    core = load_core()
    payload = core.derive_payload(args)
    payload.update(
        {
            "schema": "CORRECTED_AUAU_SHOWER_CONTRACT_WEIGHTED_WP_V1",
            "status": "DERIVED_NOT_PROMOTED",
            "model_identity": EXPECTED_PRODUCT,
            "source_label": "corrected Au+Au shower-contract full training matrix",
            "model_label": "14-feature Au+Au BDT retrained after shower-contract repair",
            "training_inputs": (
                "baseV3E plus weta33/wphi33 and centrality; calibrated full good-"
                "TowerInfo 7x7 shower grid; zero-GeV Au+Au shower-cell floor; "
                "PPG12-exact ET/eta weights"
            ),
            "shower_feature_contract": {
                "cemc_source": "TOWERINFO_CALIB_CEMC",
                "local_grid": "complete_7x7",
                "tower_acceptance": "TowerInfo::get_isGood()",
                "additional_local_mask": False,
                "rawcluster_owned_cells_for_embedding": False,
                "minimum_positive_cell_energy_gev": 0.0,
                "negative_cells_included": False,
            },
            "provenance": {
                "matrix": str(args.matrix),
                "model": str(args.model),
                "model_sha256": sha256(args.model),
                "model_metadata": str(args.model_metadata),
                "model_metadata_sha256": sha256(args.model_metadata),
                "extraction_audit": str(args.extraction_audit),
                "extraction_audit_sha256": sha256(args.extraction_audit),
            },
        }
    )
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(args.json_out)


if __name__ == "__main__":
    main()
