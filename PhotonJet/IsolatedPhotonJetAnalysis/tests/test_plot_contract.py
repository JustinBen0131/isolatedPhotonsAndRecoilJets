from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from photonjet.analysis.reduce import RecoilSelection, recoil_histogram
from photonjet.plotting.contract import (
    ROOT_EXPERIMENT_LABEL,
    compile_plot_contract,
    write_histogram_receipt,
)
from photonjet.plotting.render import render_histogram
from photonjet.plotting.root_bridge import (
    emit_root_annotation_header,
    verify_root_annotation_header,
)
from photonjet.provenance import canonical_json, write_json


ROOT = Path(__file__).resolve().parents[1]


class PlotContractTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.work = Path(self.temporary.name)
        self.trees = ROOT / "tests/fixtures/photonjet_trees_pp.root"
        self.selection = RecoilSelection(region="inclusive")
        self.dataset = ROOT / "config/datasets/pp_engineering_fixture.yaml"
        self.histogram = self.work / "histogram.json"
        self.histogram_receipt = self.work / "histogram.receipt.json"
        write_json(self.histogram, recoil_histogram([self.trees], self.selection))
        write_histogram_receipt(
            input_paths=[self.trees],
            histogram_path=self.histogram,
            selection=self.selection,
            dataset_manifest_path=self.dataset,
            receipt_path=self.histogram_receipt,
        )

    def _contract(self) -> dict:
        return compile_plot_contract(
            histogram_path=self.histogram,
            histogram_receipt_path=self.histogram_receipt,
            dataset_manifest_path=self.dataset,
        )

    def test_annotations_are_compiled_from_the_selection(self) -> None:
        contract = self._contract()
        root = contract["annotations"]["root_tlatex"]
        self.assertEqual(root["experiment"], ROOT_EXPERIMENT_LABEL)
        self.assertIn("15 #leq p_{T}^{#gamma} < 35 GeV", root["cuts"][0])
        self.assertIn("p_{T}^{jet} > 5 GeV", root["cuts"][2])
        self.assertEqual(root["cuts"][1], "Inclusive photon-jet pairs")
        histogram_receipt = json.loads(self.histogram_receipt.read_text(encoding="utf-8"))
        self.assertNotIn(self.temporary.name, json.dumps(histogram_receipt, sort_keys=True))
        self.assertEqual(histogram_receipt["histogram"]["name"], "histogram.json")

    def test_changed_cut_changes_label_without_a_label_argument(self) -> None:
        changed = RecoilSelection(region="inclusive", jet_pt_min=7.25)
        histogram = self.work / "changed.json"
        receipt = self.work / "changed.receipt.json"
        write_json(histogram, recoil_histogram([self.trees], changed))
        write_histogram_receipt(
            input_paths=[self.trees],
            histogram_path=histogram,
            selection=changed,
            dataset_manifest_path=self.dataset,
            receipt_path=receipt,
        )
        contract = compile_plot_contract(
            histogram_path=histogram,
            histogram_receipt_path=receipt,
            dataset_manifest_path=self.dataset,
        )
        self.assertIn("p_{T}^{jet} > 7.25 GeV", contract["annotations"]["root_tlatex"]["cuts"][2])

    def test_payload_or_receipt_tampering_fails_closed(self) -> None:
        payload = json.loads(self.histogram.read_text(encoding="utf-8"))
        payload["selection"]["jet_pt_min"] = 99.0
        write_json(self.histogram, payload)
        with self.assertRaisesRegex(ValueError, "histogram bytes"):
            self._contract()

    def test_dataset_substitution_fails_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "dataset manifest differs"):
            compile_plot_contract(
                histogram_path=self.histogram,
                histogram_receipt_path=self.histogram_receipt,
                dataset_manifest_path=ROOT / "config/datasets/auau_engineering_fixture.yaml",
            )

    def test_dataset_manifest_is_checked_against_tree_events(self) -> None:
        wrong = self.work / "wrong_dataset.json"
        write_json(
            wrong,
            {
                "schema": "PhotonJetDatasetManifestV1",
                "collision_system": "pp",
                "sample_kind": "engineering_fixture",
                "experiment_status": "Internal",
                "sqrt_s_gev": 200,
                "event_count": 999,
            },
        )
        with self.assertRaisesRegex(ValueError, "event_count=999"):
            write_histogram_receipt(
                input_paths=[self.trees],
                histogram_path=self.histogram,
                selection=self.selection,
                dataset_manifest_path=wrong,
                receipt_path=self.work / "wrong.receipt.json",
            )

    def test_collision_system_is_checked_against_tree_centrality(self) -> None:
        with self.assertRaisesRegex(ValueError, "Au\+Au dataset manifest"):
            write_histogram_receipt(
                input_paths=[self.trees],
                histogram_path=self.histogram,
                selection=self.selection,
                dataset_manifest_path=ROOT / "config/datasets/auau_engineering_fixture.yaml",
                receipt_path=self.work / "wrong_system.receipt.json",
            )

    def test_contract_annotation_tampering_fails_even_with_a_recomputed_hash(self) -> None:
        contract = deepcopy(self._contract())
        contract["annotations"]["matplotlib"]["dataset"] = "Au+Au, Data"
        body = {key: value for key, value in contract.items() if key != "contract_sha256"}
        contract["contract_sha256"] = hashlib.sha256(
            canonical_json(body).encode("utf-8")
        ).hexdigest()
        with self.assertRaisesRegex(ValueError, "annotations differ"):
            render_histogram(
                contract=contract,
                histogram_path=self.histogram,
                output_path=self.work / "tampered.png",
                receipt_path=self.work / "tampered.receipt.json",
            )

    def test_render_geometry_and_root_bridge_are_receipted(self) -> None:
        contract = self._contract()
        image = self.work / "plot.png"
        render_receipt = self.work / "plot.receipt.json"
        receipt = render_histogram(
            contract=contract,
            histogram_path=self.histogram,
            output_path=image,
            receipt_path=render_receipt,
        )
        self.assertTrue(image.is_file())
        self.assertEqual(receipt["image"]["name"], "plot.png")
        self.assertNotIn(self.temporary.name, json.dumps(receipt, sort_keys=True))
        self.assertFalse(receipt["qa"]["annotation_data_overlap"])
        self.assertFalse(receipt["qa"]["annotation_pair_overlap"])
        self.assertGreaterEqual(receipt["qa"]["minimum_annotation_font_points"], 11.0)
        header = emit_root_annotation_header(contract, self.work / "PlotAnnotation.generated.h")
        text = header.read_text(encoding="utf-8")
        self.assertIn(ROOT_EXPERIMENT_LABEL, text)
        self.assertIn(contract["contract_sha256"], text)
        self.assertIn("Do not edit physics labels", text)
        verification = verify_root_annotation_header(contract, header)
        self.assertTrue(verification["verified"])
        header.write_text(text + "// edited\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "edited or is stale"):
            verify_root_annotation_header(contract, header)

    def test_region_a_render_keeps_every_label_inside_the_plotting_frame(self) -> None:
        selection = RecoilSelection(region="A")
        histogram = self.work / "region_a.json"
        histogram_receipt = self.work / "region_a.receipt.json"
        write_json(histogram, recoil_histogram([self.trees], selection))
        write_histogram_receipt(
            input_paths=[self.trees],
            histogram_path=histogram,
            selection=selection,
            dataset_manifest_path=self.dataset,
            receipt_path=histogram_receipt,
        )
        contract = compile_plot_contract(
            histogram_path=histogram,
            histogram_receipt_path=histogram_receipt,
            dataset_manifest_path=self.dataset,
        )
        receipt = render_histogram(
            contract=contract,
            histogram_path=histogram,
            output_path=self.work / "region_a.png",
            receipt_path=self.work / "region_a.plot.receipt.json",
        )
        self.assertTrue(receipt["qa"]["all_annotations_inside_plotting_frame"])
        self.assertFalse(receipt["qa"]["external_header_band"])


if __name__ == "__main__":
    unittest.main()
