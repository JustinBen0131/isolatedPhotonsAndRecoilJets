#!/usr/bin/env python3
"""Mutation tests for the SDCC safe-resume admission boundary."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
GATE_PATH = HERE.parent / "sdcc_safe_resume_gate.py"
SPEC = importlib.util.spec_from_file_location("sdcc_safe_resume_gate", GATE_PATH)
assert SPEC is not None and SPEC.loader is not None
gate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gate)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SafeResumeGateTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.now = 1_900_000_000
        self.records = {}
        for name, payload in {
            "science": {"status": "PASS"},
            "permissions": {"status": "PASS", "tree_walk_performed": False},
            "duplicate": {"status": "PASS", "active_matching_jobs": 0},
            "execution": {"status": "PASS_NON_SUBMITTING", "submission_performed": False},
            "capacity": {"status": "PASS", "group_size": 7},
            "measurement": {"status": "PASS", "max_memory_mb": 1221},
            "quota": {"status": "PASS", "authoritative": True},
        }.items():
            path = self.root / f"{name}.json"
            path.write_text(json.dumps(payload, sort_keys=True) + "\n", encoding="utf-8")
            self.records[name] = {
                "path": str(path),
                "sha256": digest(path),
                "assertions": [{"path": "status", "equals": payload["status"]}],
            }

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def contract(self, operation_class: str = "watched_capacity_canary") -> dict:
        full = operation_class == "full_extraction"
        operation = {
            "class": operation_class,
            "campaign_tag": "the134_safe_resume_fixture",
            "submit_host": "sphnxuser05",
            "row_count": 13 if full else 2,
            "logical_job_count": 18577 if full else 2,
            "cluster_count": 13 if full else 2,
            "group_size": 7,
            "request_memory_mb": 3000,
            "max_materialize_per_cluster": 20 if full else 1,
            "max_idle_per_cluster": 5 if full else 1,
            "attended_submission": True,
            "automatic_job_control": False,
            "auto_memory_retry": False,
            "hold_failed_workers": True,
        }
        return {
            "schema": gate.CONTRACT_SCHEMA,
            "status": "PREPARED_NOT_SUBMITTED",
            "operation": operation,
            "login_node": {
                "explicit_path_count": 48,
                "max_explicit_paths": 256,
                "manifest_validation_rows": 32,
                "hard_manifest_validation_limit": 256,
                "max_group_files_per_row": 2048,
                "tree_walk": False,
                "unbounded_glob": False,
                "per_manifest_row_subprocess": False,
                "max_concurrent_submitters": 1,
                "submitter_timeout_seconds": 600,
                "capture_limit_bytes": 262144,
                "high_cardinality_location": "condor",
            },
            "network": {
                "allowed_services": ["sdcc_ssh", "condor_schedd", "sphenix_cvmfs"],
                "scanning": False,
                "arbitrary_upload": False,
                "host_switch": False,
                "model_directed_exploration": False,
            },
            "resources": {
                "measured_peak_memory_mb": 1221,
                "measured_witness_count": 2,
                "request_memory_mb": 3000,
                "minimum_headroom_ratio": 1.5,
                "maximum_headroom_ratio": 3.0,
                "automatic_widening": False,
                "measurement_receipt": self.records["measurement"],
            },
            "evidence": {
                "scientific_certificate": self.records["science"],
                "permission_receipt": self.records["permissions"],
                "duplicate_receipt": self.records["duplicate"],
                "execution_binding": self.records["execution"],
                "capacity_certificate": self.records["capacity"],
                "quota_certificate": self.records["quota"] if full else None,
            },
            "authority": {
                "site_admin_acknowledged": True,
                "user_approved": True,
                "exact_scope": operation_class,
                "submission_authority": True,
                "full_extraction_authority": full,
                "the121_authority": False,
                "the122_authority": False,
                "canonical_promotion": False,
            },
            "lifecycle": {
                "LOCAL_CHECK": "PASS",
                "SCIENTIFIC_CERTIFICATE": "PASS",
                "SITE_ADMISSION": "PENDING",
                "SUBMISSION": "NOT_STARTED",
                "RUNNING": "NOT_STARTED",
                "TERMINAL": "NOT_STARTED",
                "PRODUCTION_AUTHORIZED": False,
            },
            "expires_at_unix": self.now + 1800,
        }

    def write_contract(self, payload: dict, name: str = "contract.json") -> Path:
        path = self.root / name
        path.write_text(json.dumps(payload, sort_keys=True) + "\n", encoding="utf-8")
        return path

    def assert_rejected(self, mutator) -> None:
        payload = self.contract()
        mutator(payload)
        with self.assertRaises(gate.SafeResumeError):
            gate.certify_contract(self.write_contract(payload), now_seconds=self.now)

    def test_capacity_certificate_passes_and_verifies(self) -> None:
        contract = self.write_contract(self.contract())
        certificate = gate.certify_contract(contract, now_seconds=self.now)
        self.assertEqual(certificate["status"], "SITE_ADMISSION_PASS")
        self.assertFalse(certificate["submission_performed"])
        self.assertEqual(certificate["lifecycle"]["SUBMISSION"], "NOT_STARTED")
        output = self.root / "certificate.json"
        gate.write_new_json(output, certificate)
        verified = gate.verify_certificate(
            output,
            digest(output),
            operation_class="watched_capacity_canary",
            campaign_tag="the134_safe_resume_fixture",
            submit_host="sphnxuser05",
            now_seconds=self.now,
        )
        self.assertEqual(verified, certificate)

    def test_full_extraction_requires_quota_and_passes_with_it(self) -> None:
        contract = self.write_contract(self.contract("full_extraction"))
        certificate = gate.certify_contract(contract, now_seconds=self.now)
        self.assertTrue(certificate["authority"]["full_extraction_authority"])
        payload = self.contract("full_extraction")
        payload["evidence"]["quota_certificate"] = None
        with self.assertRaises(gate.SafeResumeError):
            gate.certify_contract(self.write_contract(payload, "no_quota.json"), now_seconds=self.now)

    def test_rejects_high_cardinality_and_process_fanout(self) -> None:
        for field, value in (
            ("explicit_path_count", 257),
            ("manifest_validation_rows", 33),
            ("tree_walk", True),
            ("unbounded_glob", True),
            ("per_manifest_row_subprocess", True),
            ("max_concurrent_submitters", 2),
        ):
            self.assert_rejected(lambda payload, f=field, v=value: payload["login_node"].__setitem__(f, v))

    def test_rejects_network_or_host_escape(self) -> None:
        self.assert_rejected(lambda payload: payload["network"].__setitem__("scanning", True))
        self.assert_rejected(lambda payload: payload["network"].__setitem__("host_switch", True))
        self.assert_rejected(lambda payload: payload["network"]["allowed_services"].append("public_internet"))
        self.assert_rejected(lambda payload: payload["operation"].__setitem__("submit_host", "sphnxuser04"))

    def test_rejects_resource_and_scheduler_drift(self) -> None:
        self.assert_rejected(lambda payload: payload["resources"].__setitem__("measured_peak_memory_mb", 2500))
        self.assert_rejected(lambda payload: payload["operation"].__setitem__("request_memory_mb", 8000))
        self.assert_rejected(lambda payload: payload["operation"].__setitem__("max_materialize_per_cluster", 2))
        self.assert_rejected(lambda payload: payload["operation"].__setitem__("auto_memory_retry", True))

    def test_rejects_authority_and_status_conflation(self) -> None:
        self.assert_rejected(lambda payload: payload["operation"].__setitem__("attended_submission", False))
        self.assert_rejected(lambda payload: payload["authority"].__setitem__("site_admin_acknowledged", False))
        self.assert_rejected(lambda payload: payload["authority"].__setitem__("canonical_promotion", True))
        self.assert_rejected(lambda payload: payload["lifecycle"].__setitem__("SUBMISSION", "PASS"))

    def test_rejects_evidence_drift_and_weak_permissions(self) -> None:
        payload = self.contract()
        payload["evidence"]["scientific_certificate"]["sha256"] = "0" * 64
        with self.assertRaises(gate.SafeResumeError):
            gate.certify_contract(self.write_contract(payload), now_seconds=self.now)
        weak = json.loads((self.root / "permissions.json").read_text())
        weak["status"] = "FAIL"
        (self.root / "permissions.json").write_text(json.dumps(weak) + "\n")
        payload = self.contract()
        payload["evidence"]["permission_receipt"]["sha256"] = digest(self.root / "permissions.json")
        with self.assertRaises(gate.SafeResumeError):
            gate.certify_contract(self.write_contract(payload, "weak_permission.json"), now_seconds=self.now)

    def test_certificate_reopens_evidence_and_gate_identity(self) -> None:
        contract = self.write_contract(self.contract())
        certificate = gate.certify_contract(contract, now_seconds=self.now)
        output = self.root / "reopen_certificate.json"
        gate.write_new_json(output, certificate)
        science_path = self.root / "science.json"
        science_path.write_text('{"status":"DRIFT"}\n', encoding="utf-8")
        with self.assertRaisesRegex(gate.SafeResumeError, "artifact differs"):
            gate.verify_certificate(
                output,
                digest(output),
                operation_class="watched_capacity_canary",
                campaign_tag="the134_safe_resume_fixture",
                submit_host="sphnxuser05",
                now_seconds=self.now,
            )


if __name__ == "__main__":
    unittest.main()
