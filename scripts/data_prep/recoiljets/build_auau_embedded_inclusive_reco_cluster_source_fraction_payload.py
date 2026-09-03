#!/usr/bin/env python3
"""Build the hash-bound Au+Au 1 GeV reco-cluster source-fraction payload.

Input files are offline JSON row batches.  Each row must carry ``sample_id``,
``producer_event_weight``, ``maximum_truth_jet_pt_gev``,
``centrality_percent``, and ``reco_cluster_pt_gev``.  This command performs no
ROOT/DST access and refuses to overwrite an existing payload or receipt.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any


REPO = Path(__file__).resolve().parents[3]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.data_prep.recoiljets.auau_embedded_inclusive_schema10_weighting import (  # noqa: E402
    EmbeddedInclusiveAnalysisWeightProvider,
    build_reco_cluster_source_fraction_payload,
    sha256_file,
    write_analysis_weight_artifacts,
)


def load_rows(path: Path) -> list[dict[str, Any]]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(value, dict):
        if value.get("schema") != "CanonicalSchema10AuAuEmbeddedInclusiveWeightRowBatchV1":
            raise ValueError(f"unsupported row-batch schema: {path}")
        value = value.get("rows")
    if not isinstance(value, list) or any(not isinstance(row, dict) for row in value):
        raise ValueError(f"JSON array of row objects required: {path}")
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", type=Path, required=True)
    parser.add_argument("--source-stitch-assembly", type=Path, required=True)
    parser.add_argument("--source-stitch-receipt", type=Path, required=True)
    parser.add_argument("--centrality-receipt", type=Path, required=True)
    parser.add_argument("--centrality-dependency-fingerprint", required=True)
    parser.add_argument("--pt-min", type=int, default=15)
    parser.add_argument("--pt-max", type=int, default=35)
    parser.add_argument("--output-payload", type=Path, required=True)
    parser.add_argument("--output-receipt", type=Path, required=True)
    args = parser.parse_args()
    if args.pt_max <= args.pt_min:
        raise ValueError("pt-max must exceed pt-min")
    inputs = tuple(path.resolve() for path in args.input)
    if len(set(inputs)) != len(inputs):
        raise ValueError("duplicate input row-batch path")
    rows: list[dict[str, Any]] = []
    bindings: list[dict[str, Any]] = []
    for path in inputs:
        if path.is_symlink() or not path.is_file():
            raise ValueError(f"input row batch is not a regular file: {path}")
        rows.extend(load_rows(path))
        bindings.append({
            "path": str(path),
            "sha256": sha256_file(path),
            "size_bytes": path.stat().st_size,
        })
    provider = EmbeddedInclusiveAnalysisWeightProvider.load(
        args.source_stitch_assembly,
        args.source_stitch_receipt,
        args.centrality_receipt,
        expected_centrality_dependency_fingerprint=(
            args.centrality_dependency_fingerprint
        ),
        verify_dependency_files=True,
    )
    payload = build_reco_cluster_source_fraction_payload(
        rows,
        provider,
        reco_cluster_pt_edges_gev=tuple(
            float(value) for value in range(args.pt_min, args.pt_max + 2)
        ),
        input_bindings=bindings,
    )
    receipt = write_analysis_weight_artifacts(
        args.output_payload,
        args.output_receipt,
        payload,
        provider,
    )
    print(json.dumps({
        "status": "PASS",
        "payload": str(args.output_payload.resolve()),
        "payload_sha256": receipt["payload_sha256"],
        "receipt": str(args.output_receipt.resolve()),
        "dependency_fingerprint": receipt["dependency_fingerprint"],
        "accepted_row_count": payload["accepted_row_count"],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
