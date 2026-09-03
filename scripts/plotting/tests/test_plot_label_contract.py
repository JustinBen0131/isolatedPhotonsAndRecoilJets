import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

PLOTTING = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PLOTTING))

import plot_label_contract as plot_contract_module  # noqa: E402
from plot_label_contract import (  # noqa: E402
    FULL_ACCEPTED_SIMULATION_SCOPE,
    PlotLabelContract,
)
from data_prep.recoiljets.auau_centrality_weight_contract import (  # noqa: E402
    build_canonical_contract,
)
from data_prep.recoiljets.auau_embedded_inclusive_schema10_weighting import (  # noqa: E402
    sha256_file,
    write_analysis_weight_artifacts,
    write_downstream_artifact_receipt,
)


GUARD_SCRIPT = Path(os.environ.get(
    "THESIS_OS_GUARD_PATH",
    str(PLOTTING.parent / "os" / "safety" / "codex_os_guard.py"),
))
plot_guard = None
if GUARD_SCRIPT.is_file():
    GUARD_SPEC = importlib.util.spec_from_file_location("plot_contract_guard", GUARD_SCRIPT)
    assert GUARD_SPEC and GUARD_SPEC.loader
    plot_guard = importlib.util.module_from_spec(GUARD_SPEC)
    GUARD_SPEC.loader.exec_module(plot_guard)


DEPENDENCY_FINGERPRINT = "d" * 64


class FakeCentralityWeightContract:
    status = "READY_FOR_CANONICAL_AUAU_SIM_APPLICATION"
    contract_fingerprint = "c" * 64
    dependency_fingerprint = DEPENDENCY_FINGERPRINT

    def to_payload(self):
        return {
            "schema": "AuAuCentralityWeightContractV2",
            "status": self.status,
            "application_order": [
                "existing_event_weight",
                "canonical_source_stitch",
                "auau_centrality_reweight",
            ],
            "contract_fingerprint": self.contract_fingerprint,
            "dependency_fingerprint": self.dependency_fingerprint,
        }


def contract(tmp_path: Path, **overrides) -> PlotLabelContract:
    source = tmp_path / "source.npz"
    source.write_bytes(b"fixture")
    values = dict(
        system="auau",
        energy_label=r"$\sqrt{s_{NN}}=200$ GeV",
        centrality_label="50--80%",
        sample_label="Au+Au data",
        sample_kind="data",
        cut_lines=("tight: score > WP70",),
        input_paths=(str(source),),
    )
    values.update(overrides)
    if (
        values["system"] == "auau"
        and values["sample_kind"] in {"simulation", "mixed"}
        and "simulation_family" not in overrides
    ):
        values["simulation_family"] = "photon_signal"
    return PlotLabelContract(**values)


class PlotLabelContractTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_heavy_ion_contract_passes_with_centrality(self) -> None:
        self.assertEqual(contract(self.root).validate()["status"], "PASS")

    def test_heavy_ion_contract_refuses_missing_centrality(self) -> None:
        with self.assertRaisesRegex(ValueError, "centrality"):
            contract(self.root, centrality_label=None).validate()

    def test_heavy_ion_contract_refuses_pp_energy_symbol(self) -> None:
        with self.assertRaisesRegex(ValueError, "s_NN"):
            contract(self.root, energy_label=r"$\sqrt{s}=200$ GeV").validate()

    def test_bdt_qa_contract_requires_human_dataset_role(self) -> None:
        with self.assertRaisesRegex(ValueError, "training set"):
            contract(self.root, plot_kind="bdt_qa", dataset_role=None).validate()

    def test_bdt_qa_contract_requires_visible_dataset_role(self) -> None:
        with self.assertRaisesRegex(ValueError, "visibly label"):
            contract(
                self.root,
                plot_kind="bdt_qa",
                dataset_role="training set",
                sample_label="Embedded photon+jet simulation",
            ).validate()

    def test_bdt_qa_contract_accepts_visible_dataset_role(self) -> None:
        result = contract(
            self.root,
            plot_kind="bdt_qa",
            dataset_role="validation set",
            sample_label="Embedded photon+jet simulation — validation set",
            sample_kind="simulation",
            diagnostic_weight_exception="BDT validation-set performance QA",
        ).validate()
        self.assertEqual(result["plot_audit_fields"]["dataset_role"], "validation set")
        self.assertIsNone(result["plot_audit_fields"]["simulation_scope"])

    def test_downstream_science_refuses_dataset_role(self) -> None:
        with self.assertRaisesRegex(ValueError, "BDT-performance QA"):
            contract(
                self.root,
                plot_kind="abcd_closure",
                dataset_role="validation set",
            ).validate()

    def test_non_bdt_plot_refuses_partition_language(self) -> None:
        with self.assertRaisesRegex(ValueError, "partition language"):
            contract(
                self.root,
                plot_kind="diagnostic",
                sample_label="Embedded simulation, buckets 0--5",
            ).validate()

    def test_downstream_simulation_requires_full_accepted_scope(self) -> None:
        with self.assertRaisesRegex(ValueError, "full accepted simulation"):
            contract(
                self.root,
                plot_kind="purity_closure",
                sample_kind="simulation",
            ).validate()

    def test_downstream_simulation_refuses_split_scope(self) -> None:
        with self.assertRaisesRegex(ValueError, "full accepted simulation"):
            contract(
                self.root,
                plot_kind="prompt_leakage",
                sample_kind="simulation",
                simulation_scope="validation set only",
            ).validate()

    def test_downstream_science_cannot_hide_behind_diagnostic_name(self) -> None:
        with self.assertRaisesRegex(ValueError, "full accepted simulation"):
            contract(
                self.root,
                plot_kind="abcd_closure_diagnostic",
                sample_kind="simulation",
                diagnostic_weight_exception="method-development centrality diagnostic",
            ).validate()

    def test_contract_refuses_unknown_sample_kind(self) -> None:
        with self.assertRaisesRegex(ValueError, "sample kind"):
            contract(self.root, sample_kind="monte-carlo").validate()

    def test_auau_simulation_physics_requires_analysis_weight_receipt(self) -> None:
        with self.assertRaisesRegex(ValueError, "canonical analysis-weight receipt"):
            contract(
                self.root,
                sample_kind="simulation",
                simulation_scope=FULL_ACCEPTED_SIMULATION_SCOPE,
            ).validate()

    def test_auau_simulation_requires_explicit_family(self) -> None:
        with self.assertRaisesRegex(ValueError, "simulation_family"):
            contract(
                self.root,
                sample_kind="simulation",
                simulation_scope=FULL_ACCEPTED_SIMULATION_SCOPE,
                simulation_family=None,
            ).validate()

    def test_auau_mixed_physics_requires_analysis_weight_receipt(self) -> None:
        with self.assertRaisesRegex(ValueError, "canonical analysis-weight receipt"):
            contract(
                self.root,
                sample_kind="mixed",
                simulation_scope=FULL_ACCEPTED_SIMULATION_SCOPE,
            ).validate()

    def test_exact_current_analysis_weight_receipt_passes(self) -> None:
        receipt = self.root / "analysis_weight_receipt.json"
        receipt.write_text("{}\n", encoding="utf-8")
        with mock.patch.object(
            plot_contract_module,
            "load_contract_receipt",
            return_value=FakeCentralityWeightContract(),
        ) as loader:
            result = contract(
                self.root,
                sample_kind="simulation",
                simulation_scope=FULL_ACCEPTED_SIMULATION_SCOPE,
                analysis_weight_receipt_path=str(receipt),
                analysis_weight_dependency_fingerprint=DEPENDENCY_FINGERPRINT,
            ).validate()
        loader.assert_called_once_with(
            receipt,
            expected_dependency_fingerprint=DEPENDENCY_FINGERPRINT,
            require_ready=True,
            verify_dependency_files=False,
        )
        self.assertEqual(
            result["analysis_weight_provenance"]["receipt"]["contract_fingerprint"],
            "c" * 64,
        )
        self.assertEqual(
            result["plot_audit_fields"]["analysis_weight_provenance"],
            result["analysis_weight_provenance"],
        )
        self.assertEqual(
            result["plot_audit_fields"]["simulation_scope"],
            FULL_ACCEPTED_SIMULATION_SCOPE,
        )

    def test_stale_analysis_weight_receipt_is_rejected(self) -> None:
        receipt = self.root / "analysis_weight_receipt.json"
        receipt.write_text("{}\n", encoding="utf-8")
        with mock.patch.object(
            plot_contract_module,
            "load_contract_receipt",
            side_effect=ValueError("contract dependency fingerprint is stale"),
        ):
            with self.assertRaisesRegex(ValueError, "fingerprint is stale"):
                contract(
                    self.root,
                    sample_kind="simulation",
                    simulation_scope=FULL_ACCEPTED_SIMULATION_SCOPE,
                    analysis_weight_receipt_path=str(receipt),
                    analysis_weight_dependency_fingerprint=DEPENDENCY_FINGERPRINT,
                ).validate()

    def test_receipt_path_and_dependency_fingerprint_are_required_together(self) -> None:
        receipt = self.root / "analysis_weight_receipt.json"
        receipt.write_text("{}\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "must be provided together"):
            contract(
                self.root,
                sample_kind="simulation",
                simulation_scope=FULL_ACCEPTED_SIMULATION_SCOPE,
                analysis_weight_receipt_path=str(receipt),
            ).validate()

    def test_explicit_diagnostic_exception_passes_for_auau_simulation(self) -> None:
        result = contract(
            self.root,
            sample_kind="simulation",
            plot_kind="diagnostic",
            diagnostic_weight_exception="historical unweighted centrality-shape comparison",
        ).validate()
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(
            result["analysis_weight_provenance"],
            {
                "diagnostic_exception": {
                    "reason": "historical unweighted centrality-shape comparison"
                }
            },
        )

    def test_combined_simulation_family_is_diagnostic_only(self) -> None:
        result = contract(
            self.root,
            sample_kind="mixed",
            plot_kind="diagnostic",
            simulation_family="photon_signal_and_inclusive_background",
            diagnostic_weight_exception="unit-area combined-family shape diagnostic",
        ).validate()
        self.assertEqual(
            result["plot_audit_fields"]["simulation_family"],
            "photon_signal_and_inclusive_background",
        )

        receipt = self.root / "combined_receipt.json"
        receipt.write_text("{}\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "diagnostic-only"):
            contract(
                self.root,
                sample_kind="mixed",
                simulation_scope=FULL_ACCEPTED_SIMULATION_SCOPE,
                simulation_family="photon_signal_and_inclusive_background",
                analysis_weight_receipt_path=str(receipt),
                analysis_weight_dependency_fingerprint=DEPENDENCY_FINGERPRINT,
            ).validate()

    def test_apply_to_audit_stamps_guard_facing_fields(self) -> None:
        audit = {"x_label": "centrality", "y_label": "entries"}
        result = contract(
            self.root,
            sample_kind="simulation",
            plot_kind="diagnostic",
            diagnostic_weight_exception="historical unweighted comparison",
        ).apply_to_audit(audit)
        self.assertEqual(audit, {
            "x_label": "centrality",
            "y_label": "entries",
            **result["plot_audit_fields"],
        })

    @unittest.skipUnless(plot_guard is not None, "private render guard unavailable")
    def test_real_receipt_flows_from_label_contract_through_render_guard(self) -> None:
        roles = (
            "data_target",
            "sim_stitch",
            "centrality_calibration",
            "tree_selection",
            "derivation_code",
        )
        authorities = {
            role: {
                "authority_id": "integration-fixture-{}".format(role),
                "sha256": str(index + 1) * 64,
                "fixture": True,
            }
            for index, role in enumerate(roles)
        }
        weight_contract = build_canonical_contract(
            authorities=authorities,
            data_counts=tuple(float(index + 1) for index in range(16)),
            family_counts={
                "photonjet": tuple(float(16 - index) for index in range(16)),
                "inclusivejet": tuple(2.0 for _ in range(16)),
            },
        )
        receipt = self.root / "canonical_weight_receipt.json"
        weight_contract.write_receipt(receipt)
        audit = {"x_label": "centrality", "y_label": "weighted entries"}
        contract(
            self.root,
            sample_kind="simulation",
            sample_label="Embedded photon+jet simulation",
            simulation_scope=FULL_ACCEPTED_SIMULATION_SCOPE,
            analysis_weight_receipt_path=str(receipt),
            analysis_weight_dependency_fingerprint=weight_contract.dependency_fingerprint,
        ).apply_to_audit(audit)
        self.assertEqual(plot_guard.validate_plot_analysis_weight_provenance(audit), [])
        self.assertEqual(plot_guard.validate_plot_simulation_scope(audit), [])

    @unittest.skipUnless(plot_guard is not None, "private render guard unavailable")
    def test_composite_receipt_and_output_bytes_flow_through_render_guard(self) -> None:
        from scripts.data_prep.recoiljets.tests.test_auau_embedded_inclusive_schema10_weighting import (
            WeightFixture,
        )

        fixture_root = self.root / "composite"
        fixture_root.mkdir()
        fixture = WeightFixture(fixture_root)
        provider = fixture.provider()
        payload_path = fixture_root / "payload.json"
        analysis_receipt_path = fixture_root / "analysis_receipt.json"
        write_analysis_weight_artifacts(
            payload_path,
            analysis_receipt_path,
            fixture.analysis_payload(provider),
            provider,
        )
        candidate = fixture_root / "candidate.npz"
        png = fixture_root / "plot.png"
        audit_path = fixture_root / "audit.json"
        downstream_path = fixture_root / "downstream.json"
        candidate.write_bytes(b"exact candidate")
        png.write_bytes(b"exact plot bytes")
        audit = {
            "schema": "FixtureInclusivePlotAuditV1",
            "status": "PASS",
            "x_label": "x",
            "y_label": "weighted yield",
            "output_png": str(png.resolve()),
            "candidate_blocks": [str(candidate.resolve())],
            "analysis_weight_receipt_sha256": sha256_file(analysis_receipt_path),
            "analysis_weight_dependency_fingerprint": provider.dependency_fingerprint,
            "downstream_artifact_receipt_path": str(downstream_path.resolve()),
            "generator": {
                "path": str(Path(__file__).resolve()),
                "sha256": sha256_file(Path(__file__).resolve()),
                "size_bytes": Path(__file__).resolve().stat().st_size,
            },
        }
        PlotLabelContract(
            system="auau",
            energy_label=r"$\sqrt{s_{NN}}=200$ GeV",
            centrality_label="0--80%",
            sample_label="Embedded inclusive simulation",
            sample_kind="simulation",
            plot_kind="physics",
            simulation_scope=FULL_ACCEPTED_SIMULATION_SCOPE,
            simulation_family="inclusive_background",
            cut_lines=("inclusive candidates",),
            input_paths=(str(candidate),),
            analysis_weight_receipt_path=str(analysis_receipt_path),
            analysis_weight_dependency_fingerprint=provider.dependency_fingerprint,
        ).apply_to_audit(audit)
        audit_path.write_text(json.dumps(audit, sort_keys=True) + "\n", encoding="utf-8")
        write_downstream_artifact_receipt(
            downstream_path,
            analysis_weight_receipt_path=analysis_receipt_path,
            expected_dependency_fingerprint=provider.dependency_fingerprint,
            audit_path=audit_path,
            artifacts={"plot_png": png},
            input_artifacts=(candidate,),
        )
        self.assertEqual(
            plot_guard.validate_plot_analysis_weight_provenance(
                audit, rendered_path=png, audit_path=audit_path
            ),
            [],
        )
        candidate.write_bytes(b"drifted candidate")
        self.assertTrue(any(
            "downstream artifact receipt failed validation" in error
            for error in plot_guard.validate_plot_analysis_weight_provenance(
                audit, rendered_path=png, audit_path=audit_path
            )
        ))

    def test_diagnostic_exception_is_refused_for_physics_plot(self) -> None:
        with self.assertRaisesRegex(ValueError, "only for diagnostic"):
            contract(
                self.root,
                sample_kind="simulation",
                simulation_scope=FULL_ACCEPTED_SIMULATION_SCOPE,
                diagnostic_weight_exception="not nominal",
            ).validate()

    def test_receipt_and_exception_are_mutually_exclusive(self) -> None:
        receipt = self.root / "analysis_weight_receipt.json"
        receipt.write_text("{}\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "mutually exclusive"):
            contract(
                self.root,
                sample_kind="mixed",
                plot_kind="diagnostic",
                analysis_weight_receipt_path=str(receipt),
                analysis_weight_dependency_fingerprint=DEPENDENCY_FINGERPRINT,
                diagnostic_weight_exception="historical comparison",
            ).validate()

    def test_data_plot_cannot_claim_simulation_weight_exception(self) -> None:
        with self.assertRaisesRegex(ValueError, r"only to Au\+Au simulation or mixed"):
            contract(
                self.root,
                plot_kind="diagnostic",
                diagnostic_weight_exception="not applicable to data",
            ).validate()

    def test_contract_refuses_missing_input(self) -> None:
        with self.assertRaises(FileNotFoundError):
            contract(self.root, input_paths=(str(self.root / "missing.npz"),)).validate()


if __name__ == "__main__":
    unittest.main()
