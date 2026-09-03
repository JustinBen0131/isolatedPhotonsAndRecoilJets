from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts.data_prep.recoiljets.validate_the291_downstream_weighting import (
    COLLABORATOR_CONTRACT_PATHS,
    GOVERNANCE_CONTRACT_PATHS,
    build_policy,
    validate_repository,
)
from scripts.plotting.truth_purity.make_auau_current_reference_inclusive_candidate_yields_by_centrality import (
    read_candidate_block,
)


VALID_NOTE = (
    "Au+Au embedded-inclusive event_weight is not complete. "
    "Use ownership_effective_cross_section / Npass and then centrality. "
    "Otherwise report NOT_READY_FROM_PACKAGE_ALONE.\n"
)
VALID_WORKFLOW = """name: AuAu inclusive analysis-weight contract
jobs:
  fail-closed-weight-contract:
    steps:
      - run: python3 -m pip install numpy==1.26.4
      - run: python3 scripts/data_prep/recoiljets/validate_the291_downstream_weighting.py
      - run: python3 -m unittest test_auau_embedded_inclusive_schema10_weighting test_validate_the291_downstream_weighting test_plot_label_contract
"""


class THE291SurfaceValidatorTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.repo = Path(self.temporary.name)
        for relative in GOVERNANCE_CONTRACT_PATHS:
            path = self.repo / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(
                "!/.github/workflows/auau-inclusive-weight-contract.yml\n"
                if relative == ".gitignore"
                else VALID_WORKFLOW,
                encoding="utf-8",
            )
        for relative in COLLABORATOR_CONTRACT_PATHS:
            path = self.repo / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            if relative.endswith("README.md") or relative.endswith("BRANCH_REFERENCE.md"):
                path.write_text(VALID_NOTE, encoding="utf-8")
            else:
                path.write_text(
                    "// Au+Au embedded-inclusive package example\n"
                    "// raw/unweighted; event_weight is not applied\n",
                    encoding="utf-8",
                )
        self.legacy = self.repo / "scripts/plotting/legacy.py"
        self.legacy.parent.mkdir(parents=True, exist_ok=True)
        self.legacy.write_text(
            '"""Au+Au embedded-inclusive historical plot weight diagnostic."""\n',
            encoding="utf-8",
        )
        self.policy = self.repo / "surface.json"
        self._seal()

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _seal(self) -> None:
        self.policy.write_text(
            json.dumps(build_policy(self.repo), indent=2) + "\n", encoding="utf-8"
        )

    def _errors(self):
        return validate_repository(
            self.repo, self.policy, require_canonical_surface=False
        )[1]

    def test_exact_sealed_surface_passes(self) -> None:
        self.assertEqual(self._errors(), [])

    def test_added_relevant_consumer_fails_inventory(self) -> None:
        added = self.repo / "scripts/plotting/new_nominal.py"
        added.write_text(
            '"""Au+Au embedded inclusive plot using an analysis weight."""\n',
            encoding="utf-8",
        )
        self.assertTrue(any("inventory drifted" in error for error in self._errors()))

    def test_changed_sealed_file_fails_hash(self) -> None:
        self.legacy.write_text(
            '"""Au+Au embedded-inclusive historical plot weight diagnostic changed."""\n',
            encoding="utf-8",
        )
        self.assertTrue(any("hash drifted" in error for error in self._errors()))

    def test_resealed_workflow_drift_still_fails_semantics(self) -> None:
        workflow = self.repo / ".github/workflows/auau-inclusive-weight-contract.yml"
        workflow.write_text(
            VALID_WORKFLOW.replace(" test_plot_label_contract", ""),
            encoding="utf-8",
        )
        self._seal()
        self.assertTrue(any("test_plot_label_contract" in error for error in self._errors()))

    def test_policy_reseal_cannot_weaken_declared_prohibitions(self) -> None:
        payload = json.loads(self.policy.read_text(encoding="utf-8"))
        payload["prohibitions"] = []
        self.policy.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        self.assertTrue(any("policy prohibitions differ" in error for error in self._errors()))

    def test_resealed_collaborator_weight_application_still_fails_semantics(self) -> None:
        quick = self.repo / (
            "scripts/data_prep/recoiljets/collaboration/assets/examples/quick_start.C"
        )
        quick.write_text(
            '// Au+Au embedded-inclusive raw/unweighted\n'
            'pairs.SetBranchAddress("event_weight", &event_weight);\n'
            'hist->Fill(xjgamma, event_weight);\n',
            encoding="utf-8",
        )
        self._seal()
        self.assertTrue(any("quick_start.C still applies" in error for error in self._errors()))

    def test_resealed_generic_cross_section_over_ngen_consumer_still_fails_ast_gate(self) -> None:
        consumer = self.repo / "scripts/plotting/new_nominal.py"
        consumer.parent.mkdir(parents=True, exist_ok=True)
        consumer.write_text(
            '"""Au+Au embedded-inclusive plot weight builder."""\n'
            "bad_weight = cross_section / generated_events\n",
            encoding="utf-8",
        )
        self._seal()
        self.assertTrue(any("cross-section/Ngen" in error for error in self._errors()))

    def test_resealed_nominal_consumer_private_dependency_still_fails(self) -> None:
        consumer = self.repo / (
            "scripts/plotting/truth_purity/"
            "make_auau_current_reference_inclusive_candidate_yields_by_centrality.py"
        )
        consumer.parent.mkdir(parents=True, exist_ok=True)
        consumer.write_text(
            '"""Au+Au embedded inclusive plot using an analysis weight."""\n'
            "# agent_context private import must never be required here\n",
            encoding="utf-8",
        )
        self._seal()
        self.assertTrue(any("forbidden private dependency" in error for error in self._errors()))


class PublicNominalCandidateReaderTest(unittest.TestCase):
    def _write_block(self, path: Path, *, weight: float) -> None:
        metadata = {
            "schema": "THE259PurityFactorialCandidateBlockV1",
            "status": "COMPLETE_EXISTING_RICH_TREE_CANDIDATE_BLOCK",
            "candidate_count": 1,
            "feature_count": 2,
            "system": "auau",
            "sample_id": "auau_jet12",
            "source_origin": "inclusive",
            "source_occurrence_id": "occurrence-1",
        }
        encoded = np.frombuffer(
            (json.dumps(metadata, sort_keys=True, separators=(",", ":")) + "\n").encode(),
            dtype=np.uint8,
        )
        with path.open("wb") as stream:
            np.savez_compressed(
                stream,
                metadata_json_utf8=encoded,
                event_ids=np.asarray([[1, 2]], dtype=np.uint64),
                candidate_ids=np.asarray([[3, 4]], dtype=np.uint64),
                prompt_valid=np.asarray([0], dtype=np.uint8),
                prompt_ids=np.asarray([[0, 0]], dtype=np.uint64),
                ptgamma=np.asarray([20.0], dtype=np.float64),
                eta=np.asarray([0.1], dtype=np.float32),
                centrality=np.asarray([10.0], dtype=np.float32),
                vertex_z=np.asarray([0.0], dtype=np.float32),
                score=np.asarray([0.8], dtype=np.float32),
                isolation=np.asarray([1.0], dtype=np.float32),
                isolation_threshold=np.asarray([2.0], dtype=np.float32),
                isolation_sideband_threshold=np.asarray([2.0], dtype=np.float32),
                weight=np.asarray([weight], dtype=np.float64),
                active_preselection_state=np.asarray([1], dtype=np.int8),
                prompt_match_metric=np.asarray([np.inf], dtype=np.float32),
                ordered_features=np.asarray([[0.2, 0.3]], dtype=np.float32),
                xj_offsets=np.asarray([0, 1], dtype=np.uint64),
                xj_values=np.asarray([0.7], dtype=np.float32),
            )

    def test_public_reader_accepts_complete_block(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "block.npz"
            self._write_block(path, weight=2.5)
            metadata, records = read_candidate_block(path)
        self.assertEqual(metadata["source_occurrence_id"], "occurrence-1")
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].weight, 2.5)

    def test_public_reader_rejects_negative_weight(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad-block.npz"
            self._write_block(path, weight=-1.0)
            with self.assertRaisesRegex(ValueError, "nonnegative"):
                read_candidate_block(path)


if __name__ == "__main__":
    unittest.main()
