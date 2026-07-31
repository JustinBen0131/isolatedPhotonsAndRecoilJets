#!/usr/bin/env python3
"""Focused tests for the controller-bound THE-134 submitter admission."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import tempfile
import time
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[5]
SUBMITTER = REPO / "scripts/sdcc/runtime/condor/RecoilJets_Condor_submit.sh"


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FullExtractionAdmissionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.source_list = self.root / "source.list"
        self.source_list.write_text("/a\t/b\t/c\t/d\t/e\n", encoding="utf-8")
        self.config = self.root / "analysis.yaml"
        self.config.write_text("coneR: [0.3]\n", encoding="utf-8")
        self.authorization = self.root / "authorization.json"
        self.execution = self.root / "execution.json"
        self.frozen_sha = hashlib.sha256(b"frozen").hexdigest()
        self.row_id = "pp_signal_photon5"
        self.output_namespace = self.root / "output" / self.row_id
        self.submit_namespace = self.root / "submit" / self.row_id
        campaign = {
            "tag": "unit-the134-full",
            "output_root": str(self.root / "output"),
            "evidence_root": str(self.root / "evidence"),
            "submit_root": str(self.root / "submit"),
        }
        authorization = {
            "schema": "THE134_FULL_EXTRACTION_SUBMISSION_AUTHORIZATION_V1",
            "status": "PASS_EXACT_SUBMISSION_AUTHORIZED",
            "campaign": campaign,
            "bindings": {},
            "limits": {
                "row_count": 13,
                "job_count": 13,
                "group_size": 7,
                "request_memory_mb": 3000,
                "output_pair_count": 13,
            },
            "duplicate_guard": {
                "status": "PASS_NO_ACTIVE_DUPLICATE",
                "active_matching_jobs": 0,
                "fingerprint_sha256": self.frozen_sha,
            },
            "capacity_guard": {
                "status": "PASS",
                "group_size": 7,
                "request_memory_mb": 3000,
            },
            "quota_guard": {
                "status": "PASS",
                "authoritative": True,
                "certificate_sha256": self.frozen_sha,
            },
            "approval": {
                "explicit": True,
                "scope": "THE134_FULL_SOURCE_COMPLETE_GROUP7_EXTRACTION",
                "principal": "unit-test",
                "evidence_sha256": self.frozen_sha,
            },
            "provenance": {
                "codex_chat_name": "unit-test",
                "codex_thread_id": "unit-test",
            },
            "expires_at_unix": int(time.time()) + 3600,
        }
        self.authorization.write_text(
            json.dumps(authorization, sort_keys=True), encoding="utf-8"
        )
        authorization_link = {
            "path": str(self.authorization.resolve()),
            "sha256": file_sha256(self.authorization),
        }
        row = {
            "row_id": self.row_id,
            "system": "pp",
            "dataset": "isSim",
            "sample": "run28_photonjet5",
            "expected_job_count": 1,
            "row_fingerprint_sha256": self.frozen_sha,
            "full_source_manifest_sha256": self.frozen_sha,
            "config_sha256": file_sha256(self.config),
            "code_sha256": self.frozen_sha,
            "submitter": {
                "path": str(SUBMITTER.resolve()),
                "sha256": file_sha256(SUBMITTER),
            },
            "analysis_output_namespace": str(self.output_namespace),
            "submit_namespace": str(self.submit_namespace),
            "evidence_namespace": f"{campaign['evidence_root']}/{self.row_id}",
            "training_sidecar_template": (
                f"{self.output_namespace}/training_views/"
                "$(Cluster).$(Process).root"
            ),
            "materialization_environment": {
                "RJ_DAG_DRYRUN": "1",
                "RJ_AUTO_MERGE": "0",
                "RJ_AUTO_MEMORY_RETRY": "0",
                "RJ_AUTO_MEMORY_RETRY_MAX_RELEASES": "0",
                "RJ_REQUEST_MEMORY": "3000MB",
                "RJ_CONDOR_MAX_MATERIALIZE": "20",
                "RJ_CONDOR_MAX_IDLE": "5",
                "RJ_HOLD_FAILED_WORKERS": "1",
                "RJ_VALIDATE_SIM_INPUT_PATHS": "1",
                "RJ_VALIDATE_SIM_INPUT_MAX_LINES": "32",
                "RJ_LOGIN_NODE_MAX_PATH_VALIDATION_LINES": "256",
                "RJ_LOGIN_NODE_MAX_GROUP_FILES_PER_ROW": "2048",
                "RJ_DEST_BASE_OVERRIDE": str(self.output_namespace),
                "RJ_CONDOR_SUB_DIR": str(self.submit_namespace),
                "RJ_SUBMISSION_NAMESPACE": self.row_id,
                "RJ_THE134_MULTIVIEW_TRAINING_V1": "1",
                "RJ_THE134_MULTIVIEW_SIDECAR_ONLY_V1": "1",
                "RJ_THE134_EPHEMERAL_ANALYSIS_OUTPUT": "1",
                "RJ_PP_PHOTONID_EXTRACT_ONLY": "1",
                "RJ_PP_PHOTONID_TRAINING_TREE": "1",
                "RJ_PP_PHOTONID_TRAINING_TREE_MAX_ENTRIES": "0",
                "RJ_PPG12_PHOTON_YIELD": "1",
                "RJ_PPG12_PHOTON_YIELD_DOUBLE": "0",
                "RJ_REPLAY_SOURCE_MANIFEST_SHA256": self.frozen_sha,
                "RJ_REPLAY_CONFIG_SHA256": file_sha256(self.config),
                "RJ_REPLAY_CODE_SHA256": self.frozen_sha,
                "RJ_PROFILE_LABEL": f"{campaign['tag']}_{self.row_id}",
            },
        }
        rows = [row]
        rows.extend(
            {"row_id": f"dummy_{index}", "expected_job_count": 1}
            for index in range(1, 13)
        )
        execution = {
            "schema": "THE134_FULL_EXTRACTION_EXECUTION_MANIFEST_V1",
            "status": "READY_TO_SUBMIT_EXACT_GROUP7",
            "campaign": campaign,
            "counts": {
                "row_count": 13,
                "job_count": 13,
                "group_size": 7,
                "request_memory_mb": 3000,
                "output_pair_count": 13,
            },
            "bindings": {"authorization": authorization_link},
            "authorization": {
                "expires_at_unix": authorization["expires_at_unix"],
                "provenance": authorization["provenance"],
                "approval": authorization["approval"],
            },
            "authority": {
                "submission_authority": True,
                "broad_production_authority": False,
                "the121_authority": False,
                "the122_authority": False,
                "canonical_promotion": False,
            },
            "rows": rows,
        }
        self.execution.write_text(
            json.dumps(execution, sort_keys=True), encoding="utf-8"
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def environment(self) -> dict[str, str]:
        return {
            **os.environ,
            "ACTION": "condorDoAll",
            "DATASET": "isSim",
            "SIM_SAMPLE": "run28_photonjet5",
            "SIM_CLEAN_LIST": str(self.source_list),
            "GROUP_SIZE": "7",
            "GROUP_SIZE_EXPLICIT": "1",
            "RJ_AUTO_MERGE": "0",
            "RJ_CONFIG_YAML": str(self.config),
            "RJ_REQUEST_MEMORY": "3000MB",
            "RJ_SUBMISSION_NAMESPACE": self.row_id,
            "RJ_DEST_BASE_OVERRIDE": str(self.output_namespace),
            "RJ_CONDOR_SUB_DIR": str(self.submit_namespace),
            "RJ_THE134_MULTIVIEW_TRAINING_V1": "1",
            "RJ_THE134_MULTIVIEW_SIDECAR_ONLY_V1": "1",
            "RJ_THE134_EPHEMERAL_ANALYSIS_OUTPUT": "1",
            "RJ_PP_PHOTONID_EXTRACT_ONLY": "1",
            "RJ_PP_PHOTONID_TRAINING_TREE": "1",
            "RJ_PP_PHOTONID_TRAINING_TREE_MAX_ENTRIES": "0",
            "RJ_PPG12_PHOTON_YIELD": "1",
            "RJ_PPG12_PHOTON_YIELD_DOUBLE": "0",
            "RJ_THE134_EXTRACTION_EXECUTION_MANIFEST": str(self.execution),
            "RJ_THE134_EXTRACTION_EXECUTION_MANIFEST_SHA256": file_sha256(
                self.execution
            ),
            "RJ_THE134_EXTRACTION_AUTHORIZATION_RECEIPT": str(
                self.authorization
            ),
            "RJ_THE134_EXTRACTION_AUTHORIZATION_RECEIPT_SHA256": file_sha256(
                self.authorization
            ),
            "RJ_THE134_EXTRACTION_ROW_ID": self.row_id,
            "RJ_THE134_EXTRACTION_ROW_FINGERPRINT_SHA256": self.frozen_sha,
            "RJ_THE134_EXTRACTION_SUBMITTER_PATH": str(SUBMITTER.resolve()),
        }

    def run_gate(
        self, environment: dict[str, str] | None = None
    ) -> subprocess.CompletedProcess[str]:
        command = f"""
set -euo pipefail
err() {{ printf 'ERROR: %s\\n' "$*" >&2; }}
say() {{ printf '%s\\n' "$*"; }}
env_truthy() {{
  case "${{1:-}}" in 1|true|TRUE|yes|YES|on|ON) return 0 ;; *) return 1 ;; esac
}}
auto_merge_enabled() {{
  case "${{RJ_AUTO_MERGE:-1}}" in 0|false|FALSE|no|NO|off|OFF) return 1 ;; esac
  return 0
}}
eval "$(sed -n '/^ppg12_sha256_file()/,/^# Fail closed before any broad/p' '{SUBMITTER}' | sed '$d')"
validate_the134_full_extraction_admission
"""
        return subprocess.run(
            ["bash", "-c", command, str(SUBMITTER)],
            env=self.environment() if environment is None else environment,
            text=True,
            capture_output=True,
            check=False,
        )

    def test_exact_controller_bound_admission_passes(self) -> None:
        result = self.run_gate()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("THE134_FULL_EXTRACTION_ADMISSION_PASS", result.stdout)

    def test_exact_auau_admission_passes_without_pp_bindings(self) -> None:
        payload = json.loads(self.execution.read_text(encoding="utf-8"))
        row = payload["rows"][0]
        row["system"] = "auau"
        row["dataset"] = "isSimEmbedded"
        row["sample"] = "run28_embeddedPhoton12"
        for key in (
            "RJ_PP_PHOTONID_EXTRACT_ONLY",
            "RJ_PP_PHOTONID_TRAINING_TREE",
            "RJ_PP_PHOTONID_TRAINING_TREE_MAX_ENTRIES",
            "RJ_PPG12_PHOTON_YIELD",
            "RJ_PPG12_PHOTON_YIELD_DOUBLE",
        ):
            row["materialization_environment"].pop(key)
        self.execution.write_text(
            json.dumps(payload, sort_keys=True), encoding="utf-8"
        )
        environment = self.environment()
        environment["DATASET"] = "isSimEmbedded"
        environment["SIM_SAMPLE"] = "run28_embeddedPhoton12"
        for key in (
            "RJ_PP_PHOTONID_EXTRACT_ONLY",
            "RJ_PP_PHOTONID_TRAINING_TREE",
            "RJ_PP_PHOTONID_TRAINING_TREE_MAX_ENTRIES",
            "RJ_PPG12_PHOTON_YIELD",
            "RJ_PPG12_PHOTON_YIELD_DOUBLE",
        ):
            environment.pop(key)
        environment["RJ_THE134_EXTRACTION_EXECUTION_MANIFEST_SHA256"] = (
            file_sha256(self.execution)
        )
        result = self.run_gate(environment)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("THE134_FULL_EXTRACTION_ADMISSION_PASS", result.stdout)

    def test_full_extraction_routes_before_legacy_pp_sample_filter(self) -> None:
        source = SUBMITTER.read_text(encoding="utf-8")
        start = source.index("validate_ppg12_stitched_purity_admission() {")
        end = source.index("\n# Initializes paths for isSim mode", start)
        gate = source[start:end]
        full_extraction_route = gate.index(
            "if the134_full_extraction_requested; then"
        )
        pp_sample_filter = gate.index(
            '[[ "$sample" =~ ^run28_(photonjet(5|10|20)|jet(8|12|20|30|40))'
        )
        self.assertLess(full_extraction_route, pp_sample_filter)

    def test_missing_authorization_is_rejected(self) -> None:
        environment = self.environment()
        del environment["RJ_THE134_EXTRACTION_AUTHORIZATION_RECEIPT"]
        result = self.run_gate(environment)
        self.assertEqual(result.returncode, 99)
        self.assertIn("authorization receipt", result.stderr)

    def test_wrong_row_is_rejected(self) -> None:
        environment = self.environment()
        environment["RJ_THE134_EXTRACTION_ROW_ID"] = "pp_signal_photon10"
        result = self.run_gate(environment)
        self.assertEqual(result.returncode, 99)
        self.assertIn("missing or duplicated", result.stderr)

    def test_mutated_execution_row_is_rejected(self) -> None:
        payload = json.loads(self.execution.read_text(encoding="utf-8"))
        payload["rows"][0]["row_fingerprint_sha256"] = "f" * 64
        self.execution.write_text(
            json.dumps(payload, sort_keys=True), encoding="utf-8"
        )
        environment = self.environment()
        result = self.run_gate(environment)
        self.assertEqual(result.returncode, 99)
        self.assertIn("row fingerprint differs", result.stderr)


if __name__ == "__main__":
    unittest.main()
