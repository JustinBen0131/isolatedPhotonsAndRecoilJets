#!/usr/bin/env python3
"""Build the RESUME packet for THE-236 schema-10 sparse data (pp AND Au+Au).

Why this campaign exists
------------------------
Campaign the236_schema10_data_prod_20260818_813a2538_user04 grouped pp sources
with FLAT_ORDERED_GROUPS_OF_20, which the plan documented as "groups MAY span
runs". Fun4All does not allow that:

  Fun4AllSyncManager.cc:266: Mixing run numbers (except runnumber=0 ...) is not
  supported ... Exiting now

A pp row whose 20 files straddle a run boundary therefore reads events until it
reaches the file where the run changes and then dies. That killed 1,605 of
10,007 pp rows, each having first burned up to eight hours of farm time, because
the boundary can fall anywhere inside the group. Au+Au was never exposed: its
grouping has always been run-bounded.

The materializer is now run-bounded for BOTH systems, so this packet reruns only
the sources that no successful row has already consumed.

Why this covers Au+Au too
-------------------------
Chris, 2026-08-18: submission of new Condor jobs is disabled ahead of the
19 Aug upgrade, and "the queues are not preserved ... Any jobs still in the
queue will be lost and will have to be resubmitted (same for jobs which still
run tomorrow)". At that moment 27,479 Au+Au rows were idle and 11,394 were
running, so most of the Au+Au campaign will not survive even though its
grouping was never broken. Published outputs DO survive, because the worker
renames its candidate into place on Lustre only after validation.

So the resume set is not "the pp defect" but "everything that has no
measurement receipt", derived remotely at submission time, which is the only
moment it can be correct.

Double counting is the hazard to respect
----------------------------------------
7,727 pp rows published clean output. Their sources are DONE. Regrouping all pp
sources and rerunning them would process those events a second time and inflate
every pp yield in the final dataset, which is a silent physics error rather than
a visible failure. The route therefore derives the surviving source set remotely
(every pp source, minus every source named in a records file whose row wrote a
measurement receipt) and passes it to the materializer as --restrict-sources.
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
OUTPUT = PLAN / "resume_packet"

CAMPAIGN = "the236_schema10_resume_20260819_runbounded"
PRIOR_CAMPAIGN = "the236_schema10_data_prod_20260818_813a2538_user04"
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
        # The whole point of this campaign.
        "grouping": "RUN_BOUNDED_ORDERED_GROUPS_OF_20",
        "request_memory_mb": 3_072,
        "shower_definition": "H70",
        "source_authority": "THE-110",
    },
    "auau": {
        "list_dir": f"{SCRATCH}/dst_lists_auau",
        "list_prefix": "dst_auau_jet_pair-",
        "group_size": 10,
        # Unchanged: Au+Au was always run-bounded, which is exactly why it never
        # hit the Fun4All run-mixing failure that killed 1,605 pp rows.
        "grouping": "RUN_BOUNDED_ORDERED_GROUPS_OF_10",
        "request_memory_mb": 2_560,
        "shower_definition": "H0",
        "source_authority": "THE-221",
    },
}


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
        raise SystemExit(f"fresh recovery packet namespace required: {OUTPUT}")
    for path in (RUNTIME, MATERIALIZER, WORKER):
        if path.is_symlink() or not path.is_file():
            raise SystemExit(f"regular frozen input required: {path}")

    runtime = json.loads(RUNTIME.read_text(encoding="utf-8"))

    plan = {
        "schema": "THE236_SCHEMA10_DATA_PRODUCTION_PLAN_V1",
        "status": "APPROVED_FOR_MATERIALIZATION",
        "campaign_tag": CAMPAIGN,
        "science_schema_version": 10,
        # Unknown until the surviving source set is derived on SDCC. The
        # materializer treats a declared count as advisory and reports a
        # mismatch rather than failing, because the lists are the truth.
        "exact_total_rows": 0,
        "expansion": "REMOTE_MATERIALIZATION_FROM_RESIDENT_DST_LISTS_V1",
        "event_limit": 0,
        "recovery_of_campaign": PRIOR_CAMPAIGN,
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
        # Renaming this kills every row in preflight: the worker hard-requires
        # this exact string at run_the236_schema10_data_row.sh:50.
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
            "exact_row_count": 0,
            "rows_derived_remotely_from_lists": True,
            "rows_enumerated_locally": False,
            "expansion_controller": "materialize_production.py",
            "expansion_controller_sha256": sha256_file(MATERIALIZER),
            "production_plan_sha256": hashlib.sha256(plan_bytes).hexdigest(),
            "worker_sha256": sha256_file(WORKER),
            "bulk_submission_authorized": True,
            "event_limit": 0,
            "recovery_contract": {
                "recovers": PRIOR_CAMPAIGN,
                "defect": "PP_GROUPS_SPANNED_RUN_BOUNDARIES_PLUS_CONDOR_QUEUE_LOSS",
                "evidence": (
                    "Fun4AllSyncManager.cc:266 'Mixing run numbers ... is not "
                    "supported'; 1,605 of 10,007 pp rows exited non-zero, wall "
                    "times spread from seconds to 29,810 s because the crash "
                    "lands wherever the boundary falls inside the group."
                ),
                "queue_loss": (
                    "Condor queues are not preserved across the 2026-08-19 "
                    "upgrade, so every row without a measurement receipt must be "
                    "resubmitted regardless of which defect, if any, it hit."
                ),
                "remedy": "RUN_BOUNDED_GROUPING_FOR_BOTH_SYSTEMS",
                "double_count_guard": "RESTRICT_SOURCES_TO_ROWS_WITHOUT_MEASUREMENT_RECEIPTS",
            },
        },
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
        "schema": "THE236Schema10ResumePacketReceiptV1",
        "status": "PASS_FROZEN_UNSUBMITTED",
        "campaign_tag": CAMPAIGN,
        "recovers_campaign": PRIOR_CAMPAIGN,
        "packet_root": REMOTE_PACKET,
        "output_root": REMOTE_OUTPUT,
        "log_root": REMOTE_LOG,
        "evidence_root": REMOTE_EVIDENCE,
        "event_root": REMOTE_EVENT,
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
