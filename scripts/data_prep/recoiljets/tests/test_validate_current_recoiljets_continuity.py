#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "validate_current_recoiljets_continuity.py"
)
SPEC = importlib.util.spec_from_file_location("continuity_validator", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class ContinuityValidatorTest(unittest.TestCase):
    def fixture(self, root: Path) -> tuple[Path, Path]:
        artifacts = []
        current = {}
        manifest_artifacts = []
        for index, sample_key in enumerate(sorted(MODULE.EXPECTED_SAMPLE_KEYS)):
            payload = root / f"{sample_key}.root"
            payload.write_bytes(f"payload-{index}".encode())
            artifact_id = f"artifact-{index}"
            current[sample_key] = artifact_id
            artifacts.append(
                {
                    "id": artifact_id,
                    "sample_key": sample_key,
                    "root_paths": [str(payload)],
                    "remote_paths": [f"/remote/{sample_key}.root"],
                }
            )
            manifest_artifacts.append(
                {
                    "sample_key": sample_key,
                    "artifact_id": artifact_id,
                    "classification": "QUALIFIED_CURRENT",
                    "local_path": str(payload),
                    "remote_path": f"/remote/{sample_key}.root",
                    "bytes": payload.stat().st_size,
                    "sha256": hashlib.sha256(payload.read_bytes()).hexdigest(),
                    "root_health": {
                        "status": "PASS_FRESH_TFILE",
                        "zombie": False,
                        "recovered": False,
                        "top_level_keys": 1,
                    },
                    "source_coverage": f"fixture coverage for {sample_key}",
                    "semantic_authority": "fixture-only current authority",
                    "known_omissions": ["not final H70"],
                    "consumers": ["regression"],
                }
            )
        registry = {
            "schema_version": 1,
            "current": current,
            "artifacts": artifacts,
        }
        registry_path = root / "registry.json"
        registry_path.write_text(json.dumps(registry))
        manifest = {
            "schema": "RecoilJetsCurrentContinuityFoundationV1",
            "registry": {
                "sha256": hashlib.sha256(registry_path.read_bytes()).hexdigest()
            },
            "rerun_performed": False,
            "copy_performed": False,
            "promotion_performed": False,
            "final_h70_authority": False,
            "artifacts": manifest_artifacts,
        }
        manifest_path = root / "manifest.json"
        manifest_path.write_text(json.dumps(manifest))
        return registry_path, manifest_path

    def test_metadata_only_passes_exact_six_lane_binding(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            registry, manifest = self.fixture(Path(temp))
            receipt = MODULE.validate(registry, manifest, metadata_only=True)
            self.assertEqual(receipt["status"], "PASS")
            self.assertEqual(receipt["validated_lanes"], 6)

    def test_current_pointer_drift_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            registry, manifest = self.fixture(Path(temp))
            payload = json.loads(registry.read_text())
            payload["current"]["pp_data_merged"] = "different-artifact"
            registry.write_text(json.dumps(payload))
            manifest_payload = json.loads(manifest.read_text())
            manifest_payload["registry"]["sha256"] = hashlib.sha256(
                registry.read_bytes()
            ).hexdigest()
            manifest.write_text(json.dumps(manifest_payload))
            receipt = MODULE.validate(registry, manifest, metadata_only=True)
            self.assertEqual(receipt["status"], "FAIL")
            self.assertTrue(
                any("artifact ID" in failure for failure in receipt["failures"])
            )

    def test_final_h70_mislabel_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            registry, manifest = self.fixture(Path(temp))
            payload = json.loads(manifest.read_text())
            payload["final_h70_authority"] = True
            manifest.write_text(json.dumps(payload))
            receipt = MODULE.validate(registry, manifest, metadata_only=True)
            self.assertEqual(receipt["status"], "FAIL")
            self.assertIn(
                "current continuity data is mislabeled as final H70",
                receipt["failures"],
            )

    def test_unhealthy_root_receipt_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            registry, manifest = self.fixture(Path(temp))
            payload = json.loads(manifest.read_text())
            payload["artifacts"][0]["root_health"]["zombie"] = True
            manifest.write_text(json.dumps(payload))
            receipt = MODULE.validate(registry, manifest, metadata_only=True)
            self.assertEqual(receipt["status"], "FAIL")
            self.assertTrue(
                any("zombie or unproven" in failure for failure in receipt["failures"])
            )


if __name__ == "__main__":
    unittest.main()
