#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import uproot

TABLES = {
    "RJSourceOccurrenceV1": 1,
    "RJEventV1": 1,
    "RJPhotonCandidateV1": 1,
    "RJModelEvaluationV1": 1,
    "RJShowerCellV1": 1,
    "RJIsolationConstituentV1": 1,
    "RJIsolationWitnessV1": 1,
    "RJJetV1": 1,
    "RJJetConstituentV1": 1,
    "RJPhotonJetPairV1": 1,
    "RJTruthPhotonV1": 1,
    "RJTruthJetV1": 1,
    "RJRecoTruthLinkV1": 3,
    "RJWeightComponentV1": 1,
    "RJEventDisplaySnapshotV1": 1,
}

REQUIRED_BRANCHES = {
    "RJEventV1": {"event_id_hi", "event_id_lo", "source_occurrence_id_hi", "candidate_count", "recoil_count"},
    "RJPhotonCandidateV1": {"candidate_id_hi", "event_id_hi", "cluster_et", "ordered_features", "below15_retention_state"},
    "RJModelEvaluationV1": {"candidate_id_hi", "model_id_hi", "raw_score", "applicability_state", "wp80", "delta_wp80"},
    "RJIsolationWitnessV1": {"radius", "subtraction_method", "reconstructed_or_truth", "cone_sum", "threshold"},
    "RJRecoTruthLinkV1": {"reco_type", "reco_id_hi", "truth_type", "truth_id_hi", "link_class"},
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    failures: list[str] = []
    inventory: dict[str, dict[str, object]] = {}
    with uproot.open(args.root) as root:
        base = root["ReplayFoundationV1"]
        for table, expected_entries in TABLES.items():
            if table not in base:
                failures.append(f"missing_table:{table}")
                continue
            tree = base[table]
            branches = set(tree.keys())
            if tree.num_entries != expected_entries:
                failures.append(f"entry_count:{table}:{tree.num_entries}!={expected_entries}")
            missing = sorted(REQUIRED_BRANCHES.get(table, set()) - branches)
            if missing:
                failures.append(f"missing_branches:{table}:{','.join(missing)}")
            inventory[table] = {"entries": int(tree.num_entries), "branches": sorted(branches)}

        # Stable SHA-256-derived 128-bit identity oracle.
        expected = hashlib.sha256(b"source|pp_data|run1|seg2").hexdigest()[:32]
        src = base["RJSourceOccurrenceV1"].arrays(["source_occurrence_id_hi", "source_occurrence_id_lo"], library="np")
        observed = f"{int(src['source_occurrence_id_hi'][0]):016x}{int(src['source_occurrence_id_lo'][0]):016x}"
        if observed != expected:
            failures.append(f"identity_sha256:{observed}!={expected}")

        links = base["RJRecoTruthLinkV1"].arrays(
            ["reco_type", "truth_type", "link_class"], library="np"
        )
        classes = sorted(int(value) for value in links["link_class"])
        if classes != [0, 0, 5]:
            failures.append(f"link_classes:{classes}!=[0,0,5]")

        compression = str(root.file.compression)
        if "ZSTD" not in compression.upper():
            failures.append(f"compression:{compression}")

    payload = {
        "schema": "THE118_REPLAY_SCHEMA_SYNTHETIC_VALIDATION_V1",
        "status": "PASS" if not failures else "FAIL",
        "root": str(args.root.resolve()),
        "root_sha256": sha256(args.root),
        "compression": compression,
        "table_count": len(inventory),
        "inventory": inventory,
        "failures": failures,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": payload["status"], "table_count": len(inventory), "failures": failures}))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
