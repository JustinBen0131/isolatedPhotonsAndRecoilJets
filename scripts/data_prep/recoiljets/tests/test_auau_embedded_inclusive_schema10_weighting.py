from __future__ import annotations

import json
import math
from pathlib import Path
import tempfile
import unittest

from scripts.data_prep.recoiljets.auau_centrality_weight_contract import (
    APPLICATION_ORDER,
    build_canonical_contract,
)
from scripts.data_prep.recoiljets.auau_embedded_inclusive_schema10_weighting import (
    ANALYSIS_STATUS,
    EmbeddedInclusiveAnalysisWeightProvider,
    EmbeddedInclusiveWeightError,
    SAMPLE_IDS,
    build_reco_cluster_source_fraction_payload,
    load_analysis_weight_receipt,
    load_downstream_artifact_receipt,
    load_source_stitch_artifact,
    sha256_file,
    write_analysis_weight_artifacts,
    write_downstream_artifact_receipt,
)
from scripts.data_prep.recoiljets.pp_schema10_weighting import (
    Schema10WeightingError,
    normalize_stitching_sample,
)
from scripts.plotting.plot_label_contract import (
    FULL_ACCEPTED_SIMULATION_SCOPE,
    PlotLabelContract,
)


XSECS = {
    "auau_jet12": 120.0,
    "auau_jet20": 21.0,
    "auau_jet30": 3.1,
    "auau_jet40": 0.41,
}
NPASS = {
    "auau_jet12": 100,
    "auau_jet20": 70,
    "auau_jet30": 31,
    "auau_jet40": 41,
}
WINDOWS = {
    "auau_jet12": (12.0, 21.0),
    "auau_jet20": (21.0, 31.0),
    "auau_jet30": (31.0, 41.0),
    "auau_jet40": (41.0, None),
}


class WeightFixture:
    def __init__(self, root: Path) -> None:
        self.root = root
        self.authority = root / "authority.json"
        self.assembly = root / "assembly.json"
        self.source_receipt = root / "source_receipt.json"
        self.centrality_receipt = root / "centrality.json"
        self.input_rows = root / "accepted_rows.json"
        self.authority.write_text('{"schema":"fixture"}\n', encoding="utf-8")
        rows = {}
        for sample_id in SAMPLE_IDS:
            low, high = WINDOWS[sample_id]
            denominator = NPASS[sample_id]
            weight = XSECS[sample_id] / denominator
            rows[sample_id] = {
                "sample_id": sample_id,
                "system": "auau",
                "source_class": "inclusive_background",
                "normalization_schema": "CanonicalSchema10EmbeddedInclusiveStitchingWeightingV1",
                "normalization_denominator_events": denominator,
                "normalization_denominator_kind": "current_schema10_events_in_same_truth_jet_ownership_window",
                "ownership_effective_cross_section_pb": XSECS[sample_id],
                "generated_events": denominator + 17,
                "ownership_window": {
                    "low_gev": low,
                    "high_gev": high,
                    "high_is_unbounded": high is None,
                    "lower_edge_inclusive": True,
                    "upper_edge_inclusive": False,
                },
                "stitching_weight_pb_per_owned_event": weight,
                "bin_edges_gev": [0.0, 1.0],
                "raw_counts": [denominator],
                "generator_stitching_density_pb_per_gev": [denominator * weight],
                "generator_stitching_density_sumw2_pb2_per_gev2": [
                    denominator * weight * weight
                ],
                "same_observable_denominator_proof": {
                    "event_level_ownership_audit_half_open_count": denominator,
                    "histogram_owned_raw_count": denominator,
                    "source_missing_raw_count": 0,
                    "source_underflow_raw_count": 0,
                    "source_overflow_raw_count": 0,
                    "source_sumw_equals_raw_counts": True,
                    "source_sumw2_equals_raw_counts": True,
                },
                "generator_stitching_channel": {
                    "centrality_reweighting_applied": False,
                    "density_scale_pb_per_gev_per_raw_event": weight,
                },
            }
        assembly = {
            "schema": "CanonicalSchema10EmbeddedInclusiveStitchingAssemblyV1",
            "status": "PASS",
            "authority_path": str(self.authority),
            "authority_sha256": sha256_file(self.authority),
            "sample_family": "inclusive_background",
            "sample_ids": list(SAMPLE_IDS),
            "normalization": {
                "schema": "CanonicalSchema10EmbeddedInclusiveStitchingWeightingV1",
                "centrality_reweighting_applied": False,
                "denominator_kind": "current_schema10_events_in_same_truth_jet_ownership_window",
            },
            "samples": rows,
        }
        self.assembly.write_text(json.dumps(assembly, sort_keys=True) + "\n", encoding="utf-8")
        source_receipt = {
            "schema": "CanonicalSchema10EmbeddedInclusiveStitchingValidationV1",
            "status": "pass",
            "output_schema": "CanonicalSchema10EmbeddedInclusiveStitchingAssemblyV1",
            "assembly_path": str(self.assembly),
            "assembly_sha256": sha256_file(self.assembly),
            "authority_path": str(self.authority),
            "authority_sha256": sha256_file(self.authority),
            "sample_ids": list(SAMPLE_IDS),
            "normalization_denominators": dict(NPASS),
        }
        self.source_receipt.write_text(
            json.dumps(source_receipt, sort_keys=True) + "\n", encoding="utf-8"
        )
        counts = tuple(float(index + 1) for index in range(16))
        contract = build_canonical_contract(
            authorities={
                role: {
                    "authority_id": f"fixture-{role}",
                    "sha256": str(index + 1) * 64,
                }
                for index, role in enumerate((
                    "data_target",
                    "sim_stitch",
                    "centrality_calibration",
                    "tree_selection",
                    "derivation_code",
                ))
            },
            data_counts=counts,
            family_counts={"photonjet": counts, "inclusivejet": tuple(1.0 for _ in counts)},
        )
        self.contract = contract
        self.centrality_receipt.write_text(
            json.dumps(contract.to_payload(), sort_keys=True) + "\n", encoding="utf-8"
        )
        self.input_rows.write_text("[]\n", encoding="utf-8")

    def provider(self) -> EmbeddedInclusiveAnalysisWeightProvider:
        return EmbeddedInclusiveAnalysisWeightProvider.load(
            self.assembly,
            self.source_receipt,
            self.centrality_receipt,
            expected_centrality_dependency_fingerprint=self.contract.dependency_fingerprint,
        )

    def input_bindings(self):
        return ({
            "path": str(self.input_rows),
            "sha256": sha256_file(self.input_rows),
            "size_bytes": self.input_rows.stat().st_size,
        },)

    def analysis_payload(self, provider=None):
        provider = provider or self.provider()
        return build_reco_cluster_source_fraction_payload(
            [{
                "sample_id": "auau_jet20",
                "producer_event_weight": 2.0,
                "maximum_truth_jet_pt_gev": 21.0,
                "centrality_percent": 5.0,
                "reco_cluster_pt_gev": 15.2,
            }],
            provider,
            reco_cluster_pt_edges_gev=(15.0, 16.0, 17.0),
            input_bindings=self.input_bindings(),
        )


class EmbeddedInclusiveWeightTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.fixture = WeightFixture(Path(self.temporary.name))

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_source_loader_accepts_sigma_eff_over_npass(self) -> None:
        source = load_source_stitch_artifact(
            self.fixture.assembly, self.fixture.source_receipt
        )
        self.assertAlmostEqual(
            source.sample("auau_jet20").stitching_weight_pb_per_owned_event,
            XSECS["auau_jet20"] / NPASS["auau_jet20"],
        )

    def test_source_loader_rejects_sigma_eff_over_ngen(self) -> None:
        payload = json.loads(self.fixture.assembly.read_text())
        row = payload["samples"]["auau_jet20"]
        row["stitching_weight_pb_per_owned_event"] = (
            row["ownership_effective_cross_section_pb"] / row["generated_events"]
        )
        self.fixture.assembly.write_text(json.dumps(payload, sort_keys=True) + "\n")
        receipt = json.loads(self.fixture.source_receipt.read_text())
        receipt["assembly_sha256"] = sha256_file(self.fixture.assembly)
        self.fixture.source_receipt.write_text(json.dumps(receipt, sort_keys=True) + "\n")
        with self.assertRaisesRegex(EmbeddedInclusiveWeightError, "sigma_eff/Npass"):
            load_source_stitch_artifact(self.fixture.assembly, self.fixture.source_receipt)

    def test_complete_weight_order_and_sumw2_payload(self) -> None:
        provider = self.fixture.provider()
        result = provider.apply(
            sample_id="auau_jet20",
            producer_event_weight=2.0,
            maximum_truth_jet_pt_gev=21.0,
            centrality_percent=5.0,
        )
        centrality = self.fixture.contract.weight(5.0, "inclusivejet")
        expected = 2.0 * XSECS["auau_jet20"] / NPASS["auau_jet20"] * centrality
        self.assertAlmostEqual(result.value, expected)
        self.assertEqual(result.components[0], APPLICATION_ORDER[0])
        self.assertEqual(result.components[1], APPLICATION_ORDER[1])
        self.assertTrue(result.components[2].startswith(APPLICATION_ORDER[2] + ":"))

        payload = build_reco_cluster_source_fraction_payload(
            [
                {
                    "sample_id": "auau_jet20",
                    "producer_event_weight": 2.0,
                    "maximum_truth_jet_pt_gev": 21.0,
                    "centrality_percent": 5.0,
                    "reco_cluster_pt_gev": 15.2,
                },
                {
                    "sample_id": "auau_jet20",
                    "producer_event_weight": 3.0,
                    "maximum_truth_jet_pt_gev": 30.999,
                    "centrality_percent": 5.0,
                    "reco_cluster_pt_gev": 15.8,
                },
            ],
            provider,
            reco_cluster_pt_edges_gev=(15.0, 16.0, 17.0),
            input_bindings=self.fixture.input_bindings(),
        )
        second = 3.0 * XSECS["auau_jet20"] / NPASS["auau_jet20"] * centrality
        self.assertEqual(payload["weight_state"], ANALYSIS_STATUS)
        self.assertAlmostEqual(payload["sumw"][0], expected + second)
        self.assertAlmostEqual(payload["sumw2"][0], expected * expected + second * second)

        payload_path = self.fixture.root / "payload.json"
        receipt_path = self.fixture.root / "receipt.json"
        write_analysis_weight_artifacts(payload_path, receipt_path, payload, provider)
        receipt = load_analysis_weight_receipt(
            receipt_path,
            expected_dependency_fingerprint=provider.dependency_fingerprint,
        )
        self.assertEqual(receipt["payload_sha256"], sha256_file(payload_path))
        label = PlotLabelContract(
            system="auau",
            energy_label="sqrt{s_{NN}} = 200 GeV",
            centrality_label="0-80%",
            sample_label="Embedded inclusive simulation",
            sample_kind="simulation",
            plot_kind="physics",
            simulation_scope=FULL_ACCEPTED_SIMULATION_SCOPE,
            simulation_family="inclusive_background",
            cut_lines=("15 <= reco cluster ET < 35 GeV",),
            input_paths=(str(payload_path),),
            analysis_weight_receipt_path=str(receipt_path),
            analysis_weight_dependency_fingerprint=provider.dependency_fingerprint,
        )
        validated = label.validate()
        self.assertEqual(
            validated["analysis_weight_provenance"]["receipt"]["schema"],
            receipt["schema"],
        )

    def test_half_open_ownership_and_nonunit_bins_fail_closed(self) -> None:
        provider = self.fixture.provider()
        with self.assertRaisesRegex(EmbeddedInclusiveWeightError, "outside"):
            provider.apply(
                sample_id="auau_jet12",
                producer_event_weight=1.0,
                maximum_truth_jet_pt_gev=21.0,
                centrality_percent=10.0,
            )
        with self.assertRaisesRegex(EmbeddedInclusiveWeightError, "exact 1 GeV"):
            build_reco_cluster_source_fraction_payload(
                [], provider, reco_cluster_pt_edges_gev=(15.0, 15.5, 16.0),
                input_bindings=self.fixture.input_bindings(),
            )

    def test_analysis_loader_rejects_numerical_payload_tamper_even_when_rehashed(self) -> None:
        provider = self.fixture.provider()
        payload_path = self.fixture.root / "tampered_payload.json"
        receipt_path = self.fixture.root / "tampered_receipt.json"
        write_analysis_weight_artifacts(
            payload_path, receipt_path, self.fixture.analysis_payload(provider), provider
        )
        payload = json.loads(payload_path.read_text())
        payload["sumw"][0] += 1.0
        payload_path.write_text(json.dumps(payload, sort_keys=True) + "\n")
        receipt = json.loads(receipt_path.read_text())
        receipt["payload_sha256"] = sha256_file(payload_path)
        receipt["payload_size_bytes"] = payload_path.stat().st_size
        receipt_path.write_text(json.dumps(receipt, sort_keys=True) + "\n")
        with self.assertRaisesRegex(EmbeddedInclusiveWeightError, "sumw totals differ"):
            load_analysis_weight_receipt(
                receipt_path,
                expected_dependency_fingerprint=provider.dependency_fingerprint,
            )

    def test_analysis_loader_rejects_transitive_input_drift(self) -> None:
        provider = self.fixture.provider()
        payload_path = self.fixture.root / "payload_with_input.json"
        receipt_path = self.fixture.root / "receipt_with_input.json"
        write_analysis_weight_artifacts(
            payload_path, receipt_path, self.fixture.analysis_payload(provider), provider
        )
        self.fixture.input_rows.write_text('[{"drift":true}]\n', encoding="utf-8")
        with self.assertRaisesRegex(EmbeddedInclusiveWeightError, "input binding 0"):
            load_analysis_weight_receipt(
                receipt_path,
                expected_dependency_fingerprint=provider.dependency_fingerprint,
            )

    def test_downstream_receipt_rejects_output_byte_drift(self) -> None:
        provider = self.fixture.provider()
        payload_path = self.fixture.root / "downstream_payload.json"
        analysis_receipt_path = self.fixture.root / "downstream_analysis_receipt.json"
        write_analysis_weight_artifacts(
            payload_path,
            analysis_receipt_path,
            self.fixture.analysis_payload(provider),
            provider,
        )
        plot_path = self.fixture.root / "plot.png"
        table_path = self.fixture.root / "table.csv"
        audit_path = self.fixture.root / "audit.json"
        downstream_path = self.fixture.root / "downstream_receipt.json"
        generator_path = self.fixture.root / "plot_generator.py"
        candidate_path = self.fixture.root / "candidate_block.npz"
        plot_path.write_bytes(b"fixture-png")
        table_path.write_text("x,y\n1,2\n", encoding="utf-8")
        generator_path.write_text("# exact generator\n", encoding="utf-8")
        candidate_path.write_bytes(b"exact-candidate-block")
        analysis_receipt = json.loads(analysis_receipt_path.read_text())
        audit_payload = {
            "schema": "FixturePlotAuditV1",
            "status": "PASS",
            "collision_system": "Au+Au",
            "sample_kind": "simulation",
            "plot_kind": "physics",
            "simulation_scope": "full accepted simulation",
            "simulation_family": "inclusive_background",
            "analysis_weight_receipt_sha256": sha256_file(analysis_receipt_path),
            "analysis_weight_dependency_fingerprint": provider.dependency_fingerprint,
            "analysis_weight_provenance": {
                "path": str(analysis_receipt_path.resolve()),
                "expected_dependency_fingerprint": provider.dependency_fingerprint,
                "receipt": analysis_receipt,
            },
            "downstream_artifact_receipt_path": str(downstream_path),
            "generator": {
                "path": str(generator_path.resolve()),
                "sha256": sha256_file(generator_path),
                "size_bytes": generator_path.stat().st_size,
            },
        }
        incomplete_audit = dict(audit_payload)
        incomplete_audit.pop("analysis_weight_provenance")
        audit_path.write_text(
            json.dumps(incomplete_audit, sort_keys=True) + "\n", encoding="utf-8"
        )
        with self.assertRaisesRegex(EmbeddedInclusiveWeightError, "PlotLabelContract"):
            write_downstream_artifact_receipt(
                downstream_path,
                analysis_weight_receipt_path=analysis_receipt_path,
                expected_dependency_fingerprint=provider.dependency_fingerprint,
                audit_path=audit_path,
                artifacts={"plot_png": plot_path, "table_csv": table_path},
                input_artifacts=(candidate_path,),
            )
        self.assertFalse(downstream_path.exists())
        audit_path.write_text(
            json.dumps(audit_payload, sort_keys=True) + "\n", encoding="utf-8"
        )
        write_downstream_artifact_receipt(
            downstream_path,
            analysis_weight_receipt_path=analysis_receipt_path,
            expected_dependency_fingerprint=provider.dependency_fingerprint,
            audit_path=audit_path,
            artifacts={"plot_png": plot_path, "table_csv": table_path},
            input_artifacts=(candidate_path,),
        )
        load_downstream_artifact_receipt(
            downstream_path,
            expected_dependency_fingerprint=provider.dependency_fingerprint,
            expected_audit_path=audit_path,
            expected_artifacts={"plot_png": plot_path, "table_csv": table_path},
            expected_inputs=(candidate_path,),
        )
        plot_path.write_bytes(b"fixture-png-drift")
        with self.assertRaisesRegex(EmbeddedInclusiveWeightError, "missing or stale"):
            load_downstream_artifact_receipt(
                downstream_path,
                expected_dependency_fingerprint=provider.dependency_fingerprint,
                expected_artifacts={"plot_png": plot_path, "table_csv": table_path},
            )
        plot_path.write_bytes(b"fixture-png")
        candidate_path.write_bytes(b"candidate-block-drift")
        with self.assertRaisesRegex(EmbeddedInclusiveWeightError, "input artifact"):
            load_downstream_artifact_receipt(
                downstream_path,
                expected_dependency_fingerprint=provider.dependency_fingerprint,
            )
        candidate_path.write_bytes(b"exact-candidate-block")
        generator_path.write_text("# generator drift\n", encoding="utf-8")
        with self.assertRaisesRegex(EmbeddedInclusiveWeightError, "generator"):
            load_downstream_artifact_receipt(
                downstream_path,
                expected_dependency_fingerprint=provider.dependency_fingerprint,
            )

    def test_inclusive_plot_rejects_centrality_only_receipt(self) -> None:
        label = PlotLabelContract(
            system="auau",
            energy_label="sqrt{s_{NN}} = 200 GeV",
            centrality_label="0-80%",
            sample_label="Embedded inclusive simulation",
            sample_kind="simulation",
            plot_kind="physics",
            simulation_scope=FULL_ACCEPTED_SIMULATION_SCOPE,
            simulation_family="inclusive_background",
            cut_lines=("inclusive background",),
            input_paths=(str(self.fixture.input_rows),),
            analysis_weight_receipt_path=str(self.fixture.centrality_receipt),
            analysis_weight_dependency_fingerprint=self.fixture.contract.dependency_fingerprint,
        )
        with self.assertRaisesRegex(ValueError, "centrality-only receipt is forbidden"):
            label.validate()

    def test_legacy_generic_weighting_rejects_auau_inclusive(self) -> None:
        merged = {
            "sample_id": "auau_jet20",
            "system": "auau",
            "source_class": "inclusive_background",
            "events_total": 87,
            "raw_counts": [1],
            "sumw": [1.0],
            "sumw2": [1.0],
        }
        catalog = {
            "sample_id": "auau_jet20",
            "system": "auau",
            "source_class": "inclusive_background",
            "generated_events": 87,
            "cross_section_pb": 21.0,
            "cross_section_weight_pb_per_event": 21.0 / 87.0,
            "production_campaign_tag": "fixture",
        }
        with self.assertRaisesRegex(Schema10WeightingError, "sigma_eff/Npass"):
            normalize_stitching_sample(merged, catalog, bin_width_gev=1.0)

    def test_pp_behavior_is_unchanged(self) -> None:
        merged = {
            "sample_id": "pp_jet20",
            "system": "pp",
            "source_class": "inclusive_background",
            "events_total": 10,
            "raw_counts": [2],
            "sumw": [4.0],
            "sumw2": [8.0],
        }
        catalog = {
            "sample_id": "pp_jet20",
            "system": "pp",
            "source_class": "inclusive_background",
            "generated_events": 10,
            "cross_section_pb": 5.0,
            "cross_section_weight_pb_per_event": 0.5,
            "production_campaign_tag": "fixture",
        }
        result = normalize_stitching_sample(merged, catalog, bin_width_gev=1.0)
        self.assertAlmostEqual(result["generator_stitching_density_pb_per_gev"][0], 1.0)
        self.assertAlmostEqual(
            result["analysis_weighted_density_pb_per_gev"][0], 4.0 * 7.3113 / 10.0
        )

    def test_plot_label_keeps_photon_centrality_receipt_route(self) -> None:
        label = PlotLabelContract(
            system="auau",
            energy_label="sqrt{s_{NN}} = 200 GeV",
            centrality_label="0-80%",
            sample_label="Embedded photon-signal simulation",
            sample_kind="simulation",
            plot_kind="response",
            simulation_scope=FULL_ACCEPTED_SIMULATION_SCOPE,
            simulation_family="photon_signal",
            cut_lines=("photon signal response",),
            input_paths=(str(self.fixture.input_rows),),
            analysis_weight_receipt_path=str(self.fixture.centrality_receipt),
            analysis_weight_dependency_fingerprint=self.fixture.contract.dependency_fingerprint,
        )
        result = label.validate()
        self.assertEqual(
            result["analysis_weight_provenance"]["receipt"]["status"],
            "READY_FOR_CANONICAL_AUAU_SIM_APPLICATION",
        )


if __name__ == "__main__":
    unittest.main()
