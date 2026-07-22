#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock


SCRIPT = Path(__file__).resolve().parents[1] / "manage_recoiljets_current_artifact.py"
SPEC = importlib.util.spec_from_file_location("manage_recoiljets_current_artifact", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class CurrentArtifactProductionGateTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.base = Path(self.tmp.name)
        self.registry = self.base / "registry.json"
        self.inclusive = self.base / "inclusive.root"
        self.photon = self.base / "photon.root"
        self.inclusive.write_bytes(b"inclusive-root")
        self.photon.write_bytes(b"photon-root")
        self.replay_patch = mock.patch.object(
            MODULE,
            "replay_ppg12_production_gate",
            side_effect=lambda *_args: json.loads(
                (self.base / "production_gate_report.json").read_text()
            ),
        )
        self.replay_patch.start()

    def tearDown(self) -> None:
        self.replay_patch.stop()
        self.tmp.cleanup()

    def args(
        self,
        *,
        sample_key: str = "pp_sim_inclusivejet_merged",
        status: str = "current",
        root: Path | None = None,
        production_gate_report: str = "",
    ) -> argparse.Namespace:
        return argparse.Namespace(
            registry=self.registry,
            sample_key=sample_key,
            sample_family="pp",
            artifact_kind="final_root",
            campaign_tag="unit_test",
            status=status,
            role="unit-test artifact",
            root_path=[str(root or self.inclusive)],
            remote_path=[],
            produced_at="2026-07-20T00:00:00+00:00",
            plot_policy="",
            sample_lane="",
            canonical_status="candidate_evidence",
            contract_report="",
            promotion_basis="unit test",
            production_gate_report=production_gate_report,
            waiver=[],
            notes="",
            allow_missing=False,
        )

    def write_report(self, *, corrupt_inclusive_hash: bool = False) -> Path:
        contract = {
            "schema": "ppg12-stitched-purity-closure-contract/v1",
            "unit_test": True,
        }
        contract_path = self.base / "closure_contract.json"
        contract_path.write_text(json.dumps(contract))
        contract_sha = MODULE.canonical_payload_sha256(contract)

        admission_dir = self.base / "admission"
        admission_dir.mkdir()
        admission_reference = admission_dir / "reference_manifest.json"
        admission_candidate = admission_dir / "candidate_manifest.json"
        admission_reference.write_text(json.dumps({"role": "reference"}))
        admission_candidate.write_text(json.dumps({"role": "candidate"}))
        admission_merge = admission_dir / "merge_audit.json"
        admission_merge.write_text(
            json.dumps(
                {
                    "schema": MODULE.PPG12_MERGE_AUDIT_SCHEMA,
                    "candidate_manifest": {
                        "path": str(admission_candidate.resolve()),
                        "sha256": MODULE.file_sha256(admission_candidate),
                    },
                    "audits": [{"status": "PASS", "failures": []}],
                }
            )
        )
        admission_gate = admission_dir / "gate_report.json"
        admission_gate_payload = {
            "schema": MODULE.PPG12_GATE_REPORT_SCHEMA,
            "mode": "admit",
            "status": "PASS",
            "failure_count": 0,
            "failures": [],
        }
        admission_gate.write_text(json.dumps(admission_gate_payload))
        admission = admission_dir / "admission_manifest.json"
        admission.write_text(
            json.dumps(
                {
                    "schema": MODULE.PPG12_ADMISSION_SCHEMA,
                    "status": "PASS",
                    "contract_sha256": contract_sha,
                    "reference_manifest": {
                        "path": str(admission_reference.resolve()),
                        "sha256": MODULE.file_sha256(admission_reference),
                    },
                    "candidate_manifest": {
                        "path": str(admission_candidate.resolve()),
                        "sha256": MODULE.file_sha256(admission_candidate),
                    },
                    "merge_audit": {
                        "path": str(admission_merge.resolve()),
                        "sha256": MODULE.file_sha256(admission_merge),
                    },
                    "gate_report_payload_sha256": MODULE.canonical_payload_sha256(
                        admission_gate_payload
                    ),
                }
            )
        )

        production_reference = self.base / "production_reference.json"
        production_candidate = self.base / "production_candidate.json"
        production_reference.write_text(json.dumps({"role": "reference"}))
        production_candidate.write_text(json.dumps({"role": "production"}))
        production_merge = self.base / "production_merge_audit.json"
        production_merge.write_text(
            json.dumps(
                {
                    "schema": MODULE.PPG12_MERGE_AUDIT_SCHEMA,
                    "candidate_manifest": {
                        "path": str(production_candidate.resolve()),
                        "sha256": MODULE.file_sha256(production_candidate),
                    },
                    "audits": [
                        {"family": "inclusive", "status": "PASS", "failures": []},
                        {"family": "photon", "status": "PASS", "failures": []},
                    ],
                }
            )
        )
        historical = self.base / "historical_comparison.json"
        historical.write_text(
            json.dumps(
                {"schema": "ppg12-stitched-purity-historical-comparison/v1"}
            )
        )

        artifacts = [
            {
                "family": "inclusive",
                "path": str(self.inclusive.resolve()),
                "sha256": (
                    "0" * 64
                    if corrupt_inclusive_hash
                    else MODULE.file_sha256(self.inclusive)
                ),
            },
            {
                "family": "photon",
                "path": str(self.photon.resolve()),
                "sha256": MODULE.file_sha256(self.photon),
            },
        ]
        artifacts.sort(key=lambda row: row["family"])
        wrapper = self.base / "production_wrapper.json"
        wrapper.write_text(
            json.dumps(
                {
                    "schema": MODULE.PPG12_PRODUCTION_WRAPPER_SCHEMA,
                    "admission_sha256": MODULE.file_sha256(admission),
                    "reference_manifest": {
                        "path": str(production_reference.resolve()),
                        "sha256": MODULE.file_sha256(production_reference),
                    },
                    "candidate_manifest": {
                        "path": str(production_candidate.resolve()),
                        "sha256": MODULE.file_sha256(production_candidate),
                    },
                    "merge_audit": {
                        "path": str(production_merge.resolve()),
                        "sha256": MODULE.file_sha256(production_merge),
                    },
                    "candidate_artifacts": artifacts,
                    "assembly": {
                        "contract_sha256": contract_sha,
                        "historical_comparison": {
                            "path": str(historical.resolve()),
                            "sha256": MODULE.file_sha256(historical),
                        },
                        "candidate_artifact_set_sha256": (
                            MODULE.canonical_payload_sha256(artifacts)
                        ),
                    },
                }
            )
        )
        gate_detail = self.base / "gate_report.json"
        gate_detail_payload = {
            "schema": MODULE.PPG12_GATE_REPORT_SCHEMA,
            "mode": "verify-production",
            "status": "PASS",
            "contract": {"path": str(contract_path.resolve()), "sha256": contract_sha},
            "reference_manifest": str(production_reference.resolve()),
            "candidate_manifest": str(production_candidate.resolve()),
            "failure_count": 0,
            "failures": [],
            "summary": {
                "admission_manifest": str(admission.resolve()),
                "production_wrapper": str(wrapper.resolve()),
                "candidate_artifacts": artifacts,
                "candidate_artifact_set_sha256": MODULE.canonical_payload_sha256(
                    artifacts
                ),
            },
        }
        gate_detail.write_text(json.dumps(gate_detail_payload))
        payload = {
            "schema": MODULE.PPG12_PRODUCTION_GATE_SCHEMA,
            "status": "PASS",
            "contract_sha256": contract_sha,
            "gate_report_payload_sha256": MODULE.canonical_payload_sha256(
                gate_detail_payload
            ),
            "admission_sha256": MODULE.file_sha256(admission),
            "production_manifest_sha256": MODULE.file_sha256(wrapper),
            "reference_manifest_sha256": MODULE.file_sha256(production_reference),
            "candidate_manifest_sha256": MODULE.file_sha256(production_candidate),
            "merge_audit_sha256": MODULE.file_sha256(production_merge),
            "candidate_artifacts": artifacts,
            "candidate_artifact_set_sha256": MODULE.canonical_payload_sha256(artifacts),
        }
        path = self.base / "production_gate_report.json"
        path.write_text(json.dumps(payload))
        return path

    def test_gated_current_registration_requires_report(self) -> None:
        self.assertEqual(MODULE.command_register(self.args()), 3)
        self.assertFalse(self.registry.exists())

    def test_sync_snapshots_cannot_recreate_ungated_current_pointer(self) -> None:
        entry_id = "legacy-pass-shaped-current"
        entry = {
            "id": entry_id,
            "sample_key": "pp_sim_inclusivejet_merged",
            "status": "current",
            "root_paths": [str(self.inclusive.resolve())],
        }
        self.registry.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "artifacts": [entry],
                    "current": {"pp_sim_inclusivejet_merged": entry_id},
                }
            )
        )
        code = MODULE.command_sync_snapshots(
            argparse.Namespace(registry=self.registry)
        )
        self.assertEqual(code, 3)
        self.assertFalse(
            (
                self.registry.parent
                / "current/pp_sim_inclusivejet_merged/current.json"
            ).exists()
        )
        self.assertFalse(
            (
                self.registry.parent
                / "entries/pp_sim_inclusivejet_merged"
                / f"{entry_id}.json"
            ).exists()
        )

    def test_candidate_registration_does_not_require_report(self) -> None:
        self.assertEqual(MODULE.command_register(self.args(status="candidate")), 0)
        registry = json.loads(self.registry.read_text())
        self.assertEqual(registry["artifacts"][0]["status"], "candidate")
        self.assertEqual(registry["current"], {})

    def test_report_must_bind_exact_root_hash(self) -> None:
        report = self.write_report(corrupt_inclusive_hash=True)
        self.assertEqual(
            MODULE.command_register(self.args(production_gate_report=str(report))),
            3,
        )
        self.assertFalse(self.registry.exists())

    def test_pass_shaped_receipt_with_arbitrary_hash_is_rejected(self) -> None:
        report = self.write_report()
        payload = json.loads(report.read_text())
        payload["admission_sha256"] = "a" * 64
        report.write_text(json.dumps(payload))
        self.assertEqual(
            MODULE.command_register(self.args(production_gate_report=str(report))),
            3,
        )
        self.assertFalse(self.registry.exists())

    def test_mutated_linked_gate_detail_is_rejected(self) -> None:
        report = self.write_report()
        detail = self.base / "gate_report.json"
        payload = json.loads(detail.read_text())
        payload["status"] = "FAIL"
        detail.write_text(json.dumps(payload))
        self.assertEqual(
            MODULE.command_register(self.args(production_gate_report=str(report))),
            3,
        )
        self.assertFalse(self.registry.exists())

    def test_mutated_linked_production_wrapper_is_rejected(self) -> None:
        report = self.write_report()
        wrapper = self.base / "production_wrapper.json"
        wrapper.write_text(wrapper.read_text() + "\n")
        self.assertEqual(
            MODULE.command_register(self.args(production_gate_report=str(report))),
            3,
        )
        self.assertFalse(self.registry.exists())

    def test_mutated_linked_candidate_manifest_is_rejected(self) -> None:
        report = self.write_report()
        candidate = self.base / "production_candidate.json"
        candidate.write_text(candidate.read_text() + "\n")
        self.assertEqual(
            MODULE.command_register(self.args(production_gate_report=str(report))),
            3,
        )
        self.assertFalse(self.registry.exists())

    def test_unresolvable_merge_audit_hash_is_rejected(self) -> None:
        report = self.write_report()
        payload = json.loads(report.read_text())
        payload["merge_audit_sha256"] = "b" * 64
        report.write_text(json.dumps(payload))
        self.assertEqual(
            MODULE.command_register(self.args(production_gate_report=str(report))),
            3,
        )
        self.assertFalse(self.registry.exists())

    def test_passing_report_allows_exact_current_pointer(self) -> None:
        report = self.write_report()
        self.assertEqual(
            MODULE.command_register(self.args(production_gate_report=str(report))),
            0,
        )
        pointer = json.loads(
            (
                self.registry.parent
                / "current"
                / "pp_sim_inclusivejet_merged"
                / "current.json"
            ).read_text()
        )
        self.assertEqual(pointer["root_paths"], [MODULE.normalize_path(str(self.inclusive))])
        self.assertEqual(
            pointer["production_gate_report"],
            MODULE.normalize_path(str(report)),
        )
        self.assertEqual(
            pointer["production_gate_report_sha256"],
            MODULE.file_sha256(report),
        )

    def test_self_consistent_receipt_that_canonical_replay_does_not_emit_is_rejected(self) -> None:
        report = self.write_report()
        replayed = json.loads(report.read_text())
        replayed["admission_sha256"] = "f" * 64
        with mock.patch.object(
            MODULE, "replay_ppg12_production_gate", return_value=replayed
        ):
            self.assertEqual(
                MODULE.command_register(
                    self.args(production_gate_report=str(report))
                ),
                3,
            )
        self.assertFalse(self.registry.exists())

    def test_unrelated_sample_key_remains_unchanged(self) -> None:
        self.assertEqual(
            MODULE.command_register(
                self.args(sample_key="pp_data_merged", production_gate_report="")
            ),
            0,
        )


if __name__ == "__main__":
    unittest.main()
