#!/usr/bin/env python3
"""Build the THE-236 schema-10 sparse-data PRODUCTION packet locally.

This is the N-row generalization of build_the236_schema10_sparse_data_canary_packet.py,
which hardcoded exactly four rows.

The important structural difference from the canary is that the 51,535 production
rows are NOT enumerated here. The packet declares the expansion contract (list
directory, grouping rule, exact row count per system) and
materialize_the236_schema10_data_production.py expands the rows on SDCC, where
the DST pair lists already live. That mirrors the simulation campaign, whose
production_plan.json declares 5 samples and lets the remote materializer produce
7,145 rows. It keeps ~500k source lines off the local machine entirely and keeps
the packet small and human-readable.

AuAu only for this submission. The pp run selection is already applied to the
paired lists (1,568 lists == 1,568 GRL runs), but pp additionally declares
21,494 sources against ~200,140 segments in those runs, i.e. a segment-level
selection under source_authority THE-110 that is not derivable from the run
list. AuAu has no such gap: its 2,776 paired lists are exactly the 2,776 GRL
runs. pp is 1,075 of 51,535 rows and its candidate path is unproven anyway.

Grouping, reproduced from build_the236_schema10_data_plan.py:
  pp    flat ordered groups of 20, may span runs   21,494 sources ->  1,075 rows
  auau  run-bounded ordered groups of 10           492,280 in 2,776 runs -> 50,460 rows

Production rows carry event_limit=0 so the bounded-canary path stays off, but
exact per-row expected_processed IS required: the worker regex rejects 0 and
verifies the count against the real processed events before publishing. Those
counts come from a read-only FileCatalog SELECT the route runs on SDCC.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from typing import Any


REPO = next(p for p in Path(__file__).resolve().parents if (p / "AGENTS.md").is_file())
PLAN = REPO / (
    "agent_context/local/plans/canonical_pp_auau_photonjet_foundation_20260720/"
    "the236_schema10_data_production_20260817"
)
CANARY_PLAN = REPO / (
    "agent_context/local/plans/canonical_pp_auau_photonjet_foundation_20260720/"
    "the236_schema10_sparse_data_canary_20260817"
)
RUNTIME = CANARY_PLAN / "SPARSE_DATA_CERTIFIED_RUNTIME_RECEIPT.json"
MATERIALIZER = REPO / "scripts/sdcc/workflows/submit/materialize_the236_schema10_data_production.py"
WORKER = REPO / "scripts/sdcc/runtime/condor/run_the236_schema10_data_row.sh"
OUTPUT = PLAN / "production_packet"

CAMPAIGN = "the236_schema10_data_prod_20260818_813a2538_user04"
REMOTE_BASE = "/sphenix/tg/tg01/bulk/jbennett/thesisAna/recoiljets"
REMOTE_PACKET = f"{REMOTE_BASE}/packets/{CAMPAIGN}"
REMOTE_OUTPUT = f"{REMOTE_BASE}/outputs/{CAMPAIGN}"
REMOTE_LOG = f"{REMOTE_BASE}/logs/{CAMPAIGN}"
REMOTE_EVIDENCE = f"{REMOTE_BASE}/evidence/{CAMPAIGN}"
REMOTE_EVENT = f"/tmp/patsfan753/{CAMPAIGN}"
CONTRACT_SHA256 = "1388f6e3c89703d3f7a6d9b7f66371403b2aa4d861376af18cc0d21bd3b54df4"

SCRATCH = "/sphenix/u/patsfan753/scratch/thesisAnalysis"
SYSTEMS: dict[str, dict[str, Any]] = {
    "pp": {
        "list_dir": f"{SCRATCH}/dst_lists_pp",
        "list_prefix": "dst_ppg12_pair-",
        "group_size": 20,
        "grouping": "FLAT_ORDERED_GROUPS_OF_20",
        "request_memory_mb": 3_072,
        "shower_definition": "H70",
        "source_authority": "THE-110",
    },
    "auau": {
        "list_dir": f"{SCRATCH}/dst_lists_auau",
        "list_prefix": "dst_auau_jet_pair-",
        "group_size": 10,
        "grouping": "RUN_BOUNDED_ORDERED_GROUPS_OF_10",
        "request_memory_mb": 2_560,
        "shower_definition": "H0",
        "source_authority": "THE-221",
    },
}
EXACT_TOTAL_ROWS = 50_460


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_exclusive(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o440)
    with os.fdopen(descriptor, "wb") as stream:
        stream.write(payload)


def build() -> dict[str, Any]:
    if OUTPUT.exists() or OUTPUT.is_symlink():
        raise SystemExit(f"fresh production packet namespace required: {OUTPUT}")
    for path in (RUNTIME, MATERIALIZER, WORKER):
        if path.is_symlink() or not path.is_file():
            raise SystemExit(f"regular frozen input required: {path}")

    runtime = json.loads(RUNTIME.read_text(encoding="utf-8"))

    # The compact plan consumed by the remote materializer.
    plan = {
        "schema": "THE236_SCHEMA10_DATA_PRODUCTION_PLAN_V1",
        "status": "APPROVED_FOR_MATERIALIZATION",
        "campaign_tag": CAMPAIGN,
        "science_schema_version": 10,
        "exact_total_rows": EXACT_TOTAL_ROWS,
        "expansion": "REMOTE_MATERIALIZATION_FROM_RESIDENT_DST_LISTS_V1",
        "event_limit": 0,
        "systems": [
            {
                "system": system,
                "list_dir": spec["list_dir"],
                "list_prefix": spec["list_prefix"],
                "group_size": spec["group_size"],
                "grouping": spec["grouping"],
                "shower_definition": spec["shower_definition"],
                "request_memory_mb": spec["request_memory_mb"],
            }
            for system, spec in SYSTEMS.items()
        ],
    }
    plan_bytes = (json.dumps(plan, indent=2, sort_keys=True) + "\n").encode()

    manifest = {
        # The worker hard-requires this exact schema string at
        # run_the236_schema10_data_row.sh:50. It is the canary manifest schema and
        # must not be renamed for production, or every row dies in preflight.
        "schema": "THE236Schema10SparseDataProductionManifestV1",
        "status": "PASS_FROZEN_UNSUBMITTED",
        "campaign_tag": CAMPAIGN,
        "science_schema_version": 10,
        "data_retention_profile": "sparse_photon_analysis_v1",
        "data_retention_contract_sha256": CONTRACT_SHA256,
        "runtime_receipt_path": str(RUNTIME),
        "runtime_receipt_sha256": sha256_file(RUNTIME),
        "runtime": {"pp": runtime["pp"], "auau": runtime["auau"]},
        "replay_environment": runtime["replay_environment"],
        "production_contract": {
            "exact_row_count": EXACT_TOTAL_ROWS,
            "rows_derived_remotely_from_lists": True,
            "rows_enumerated_locally": False,
            "expansion_controller": "materialize_production.py",
            "expansion_controller_sha256": sha256_file(MATERIALIZER),
            "production_plan_sha256": hashlib.sha256(plan_bytes).hexdigest(),
            "worker_sha256": sha256_file(WORKER),
            "bulk_submission_authorized": True,
            "event_limit": 0,
            "canary_precedent": {
                "campaign": "the236_schema10_sparse_data_canary_20260817_813a2538_showerdefclosed",
                "cluster": 670064,
                "rows_exit_zero": ["pp_typical", "auau_typical", "auau_photon_positive"],
                "auau_candidate_bearing_proven": True,
                "pp_candidate_bearing_proven": False,
                "note": (
                    "auau_photon_positive exited 0 in 4477 s at 1465 MB, proving the "
                    "shower-definition fix on the AuAu/H0 candidate path. The pp/H70 "
                    "candidate row was retired before completion, so pp candidate "
                    "handling enters production unproven; pp is 1,075 of 51,535 rows."
                ),
            },
        },
        "open_risks": [
            {
                "risk": "PRODUCTION_ROW_MEMORY_UNMEASURED",
                "detail": (
                    "Every canary row ran with an event_limit of at most 20,000. A "
                    "production row processes every event in its source group, so peak "
                    "memory at full row size has never been observed. The worker fails a "
                    "row whose peak exceeds its request, so an underestimate costs rows "
                    "rather than corrupting them."
                ),
                "measured_canary_peak_mb": {"pp": 1_954, "auau": 1_465},
                "requested_mb": {s: spec["request_memory_mb"] for s, spec in SYSTEMS.items()},
                "hard_stop_mb": 4_096,
            },
        ],
        "site_contract": {
            "getenv": False,
            "request_cpus": 1,
            "request_memory_mb": {s: spec["request_memory_mb"] for s, spec in SYSTEMS.items()},
            "memory_hard_stop_mb": 4_096,
            "condor_event_log": REMOTE_EVENT,
            "persistent_stdout_stderr": REMOTE_LOG,
            "persistent_output": REMOTE_OUTPUT,
            "candidate_publication": "ROOT_PART_VALIDATE_ATOMIC_RENAME",
            "shard_size": 1_000,
            "retry": False,
            "release": False,
            "maxjobs": False,
            "on_exit_hold": False,
        },
    }
    manifest_bytes = (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode()

    write_exclusive(OUTPUT / "production_plan.json", plan_bytes)
    write_exclusive(OUTPUT / "production_manifest.json", manifest_bytes)
    write_exclusive(OUTPUT / "materialize_production.py", MATERIALIZER.read_bytes())
    os.chmod(OUTPUT / "materialize_production.py", 0o555)
    write_exclusive(OUTPUT / "run_production_row.sh", WORKER.read_bytes())
    os.chmod(OUTPUT / "run_production_row.sh", 0o555)

    receipt = {
        "schema": "THE236Schema10DataProductionPacketReceiptV1",
        "status": "PASS_FROZEN_UNSUBMITTED",
        "campaign_tag": CAMPAIGN,
        "packet_root": REMOTE_PACKET,
        "output_root": REMOTE_OUTPUT,
        "log_root": REMOTE_LOG,
        "evidence_root": REMOTE_EVIDENCE,
        "event_root": REMOTE_EVENT,
        "exact_total_rows": EXACT_TOTAL_ROWS,
        "rows_derived_remotely_from_lists": True,
        "production_plan_sha256": hashlib.sha256(plan_bytes).hexdigest(),
        "production_manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
        "materializer_sha256": sha256_file(MATERIALIZER),
        "worker_sha256": sha256_file(WORKER),
        "runtime_receipt_sha256": sha256_file(RUNTIME),
        "systems": SYSTEMS,
    }
    write_exclusive(
        OUTPUT / "packet_receipt.json",
        (json.dumps(receipt, indent=2, sort_keys=True) + "\n").encode(),
    )
    return receipt


def main() -> int:
    print(json.dumps(build(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
