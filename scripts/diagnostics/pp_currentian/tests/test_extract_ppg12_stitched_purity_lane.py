#!/usr/bin/env python3
"""Regression tests for the fail-closed stitched-purity lane extractor."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import tempfile
import unittest
from array import array
from pathlib import Path


REPO = Path(__file__).resolve().parents[4]
MODULE_PATH = (
    REPO
    / "scripts/diagnostics/pp_currentian/extract_ppg12_stitched_purity_lane.py"
)
EXECUTION_PATH = (
    REPO
    / "scripts/diagnostics/pp_currentian/produce_ppg12_stitched_purity_execution_contract.py"
)
CONTRACT_PATH = (
    REPO / "agent_context/analysis_contracts/ppg12_stitched_purity_closure.yaml"
)
SPEC = importlib.util.spec_from_file_location(
    "extract_ppg12_stitched_purity_lane", MODULE_PATH
)
assert SPEC and SPEC.loader
EXTRACTOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EXTRACTOR)
EXECUTION_SPEC = importlib.util.spec_from_file_location(
    "produce_ppg12_stitched_purity_execution_contract", EXECUTION_PATH
)
assert EXECUTION_SPEC and EXECUTION_SPEC.loader
EXECUTION = importlib.util.module_from_spec(EXECUTION_SPEC)
EXECUTION_SPEC.loader.exec_module(EXECUTION)
CONTRACT = json.loads(CONTRACT_PATH.read_text())


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


class FakeReader:
    def __init__(self, histograms: dict[str, dict[str, object]]):
        self.histograms = histograms
        self.extra_paths: dict[str, list[str]] = {}
        self.cycles: dict[str, int] = {path: 1 for path in histograms}
        self.closed = False

    def paths_for_basename(self, basename: str) -> list[str]:
        result = [path for path in self.histograms if path.rsplit("/", 1)[-1] == basename]
        result.extend(self.extra_paths.get(basename, []))
        return sorted(set(result))

    def cycle_count(self, path: str) -> int:
        return self.cycles.get(path, 0)

    def histogram(self, path: str) -> dict[str, object]:
        return copy.deepcopy(self.histograms[path])

    def close(self) -> None:
        self.closed = True


class LaneExtractorTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.root_file = self.root / "lane.root"
        self.root_file.write_bytes(b"synthetic-root-bound-by-hash" + b"x" * 60_000)
        self.evidence_files: dict[str, Path] = {}
        self.event_rows = [
            f"NONE /x/js_pp200_signal/g4hits/file{index}.root "
            f"/x/js_pp200_signal/nopileup/jets/file{index}.root NONE NONE"
            for index in range(1, 6)
        ]
        for name in EXTRACTOR.EVIDENCE_KEYS:
            path = self.root / f"{name}.evidence"
            if name in {"source_list", "event_set"}:
                path.write_text("\n".join(self.event_rows) + "\n")
            else:
                path.write_text(f"exact {name} evidence\n")
            self.evidence_files[name] = path
        self.runtime_dir = self.root / "runtime"
        self.runtime_dir.mkdir()
        self.build_receipt = self.runtime_dir / "build_receipt.json"
        write_json(self.build_receipt, {"status": "source-locked-test-runtime"})
        runtime_roles = sorted(
            {
                role
                for roles in CONTRACT["lane_runtime_contract"]
                ["required_roles_by_lane_field"].values()
                for role in roles
                if role != "lane_config"
            }
        )
        runtime_files = []
        for role in runtime_roles:
            path = self.runtime_dir / role.replace("/", "_")
            path.write_text(f"executed bytes for {role}\n")
            runtime_files.append(
                {"role": role, "path": str(path), "sha256": file_sha256(path)}
            )
        self.runtime_manifest = self.runtime_dir / "runtime_manifest.json"
        write_json(
            self.runtime_manifest,
            {
                "schema_version": 1,
                **CONTRACT["lane_runtime_contract"]["required_manifest_values"],
                "build_receipt": str(self.build_receipt),
                "build_receipt_sha256": file_sha256(self.build_receipt),
                "files": runtime_files,
            },
        )
        self.parity_path = self.root / "candidate_rows.csv"
        self.parity_path.write_text("event,candidate,route,tag,abcd\n1,2,baseV3E,tight,A\n")
        self.metadata_path = self.root / "lane_input.json"
        self.fill_path = self.root / "fills.json"
        self.execution_path = self.root / "execution_contract.json"

    def tearDown(self) -> None:
        self.temp.cleanup()

    @staticmethod
    def required(family: str) -> list[str]:
        return list(CONTRACT["families"][family]["required_observables"])

    def make_reader(self, family: str) -> FakeReader:
        histograms: dict[str, dict[str, object]] = {}
        for index, name in enumerate(self.required(family)):
            path = f"SIM/{EXTRACTOR.REGION_OBJECTS[name]}"
            scale = float(index + 1)
            histograms[path] = {
                "bin_edges": [10.0, 12.0, 14.0],
                "sumw": [scale, 1.5 * scale],
                "sumw2": [0.5 * scale * scale, 0.75 * scale * scale],
                "entries": 5.0,
                "flow_sumw": [0.0, 0.0],
                "flow_sumw2": [0.0, 0.0],
                "classname": "TH1D",
            }
        return FakeReader(histograms)

    def make_inputs(
        self,
        family: str = "inclusive",
        sample: str | None = None,
        *,
        include_parity: bool | None = None,
    ) -> tuple[dict[str, object], dict[str, object], FakeReader]:
        sample = sample or ("jet8" if family == "inclusive" else "photon5")
        include_parity = family == "photon" if include_parity is None else include_parity
        lane_id = f"{family}:{sample}:0mrad:si"
        event_rows = list(self.event_rows)
        groups = EXTRACTOR._canonical_groups(lane_id, event_rows, 5)
        reader = self.make_reader(family)
        fills = {
            name: {
                "object": f"SIM/{EXTRACTOR.REGION_OBJECTS[name]}",
                "bin_edges": [10.0, 12.0, 14.0],
                "fills": [2, 3],
                "flow_fills": [0, 0],
            }
            for name in self.required(family)
        }
        fill_payload: dict[str, object] = {
            "schema": "ppg12-stitched-purity-fill-evidence/v1",
            "lane_id": lane_id,
            "root_sha256": file_sha256(self.root_file),
            "observables": fills,
        }
        write_json(self.fill_path, fill_payload)
        environment = EXECUTION._expected_environment(
            {
                "lane_id": lane_id,
                "family": family,
                "sample": sample,
                "period": "0mrad",
                "interaction": "si",
            },
            CONTRACT,
        )
        write_json(
            self.execution_path,
            EXECUTION._build_receipt(
                lane_id,
                self.evidence_files["source_list"],
                environment,
                CONTRACT_PATH,
            ),
        )
        metadata: dict[str, object] = {
            "schema": "ppg12-stitched-purity-lane-input/v1",
            "lane_id": lane_id,
            "family": family,
            "sample": sample,
            "period": "0mrad",
            "interaction": "si",
            "group_count": 1,
            "group_size": 5,
            "group_index_start": 0,
            "event_set_row_count": 5,
            "groups": groups,
            "group_set_sha256": EXTRACTOR._payload_sha256(groups),
            "external_scale": 1.0,
            "abcd_population": "unsuffixed",
            "object_prefix": "SIM",
            "root": {"path": str(self.root_file), "sha256": file_sha256(self.root_file)},
            "evidence": {
                name: {"path": str(path), "sha256": file_sha256(path)}
                for name, path in self.evidence_files.items()
            },
            "runtime_manifest": {
                "path": str(self.runtime_manifest),
                "sha256": file_sha256(self.runtime_manifest),
            },
            "execution_contract": {
                "path": str(self.execution_path),
                "sha256": file_sha256(self.execution_path),
            },
            "fill_evidence": {
                "path": str(self.fill_path),
                "sha256": file_sha256(self.fill_path),
            },
            "candidate_parity_evidence": (
                [
                    {
                        "role": "candidate_rows",
                        "path": str(self.parity_path),
                        "sha256": file_sha256(self.parity_path),
                    }
                ]
                if include_parity
                else []
            ),
        }
        write_json(self.metadata_path, metadata)
        return metadata, fill_payload, reader

    def run_extract(self, reader: FakeReader) -> dict[str, object]:
        return EXTRACTOR.extract_lane(
            self.metadata_path,
            CONTRACT_PATH,
            reader_factory=lambda _: reader,
        )

    def rewrite(self, payload: dict[str, object], path: Path) -> None:
        write_json(path, payload)

    def refresh_fill_link(self, metadata: dict[str, object]) -> None:
        metadata["fill_evidence"] = {
            "path": str(self.fill_path),
            "sha256": file_sha256(self.fill_path),
        }
        write_json(self.metadata_path, metadata)

    def test_inclusive_happy_path_is_assembler_compatible(self) -> None:
        _, _, reader = self.make_inputs()
        output = self.run_extract(reader)
        self.assertTrue(reader.closed)
        self.assertEqual(output["schema"], "ppg12-stitched-purity-lane/v1")
        self.assertEqual(output["lane_id"], "inclusive:jet8:0mrad:si")
        self.assertEqual(set(output["observables"]), set(self.required("inclusive")))
        self.assertEqual(output["observables"]["A"]["fills"], [2.0, 3.0])
        self.assertEqual(output["merge_input"]["sha256"], file_sha256(self.root_file))
        self.assertEqual(
            output["estimator_sha256"],
            output["extraction"]["runtime_source_sets"]["estimator_sha256"]
            ["portable_role_set_sha256"],
        )
        # Directly exercise the consuming assembler's lane validator.
        assembler_path = (
            REPO
            / "scripts/diagnostics/pp_currentian/assemble_ppg12_stitched_purity_manifest.py"
        )
        spec = importlib.util.spec_from_file_location("lane_assembler_for_test", assembler_path)
        assert spec and spec.loader
        assembler = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(assembler)
        assembler._reproduce_lane_extract = lambda *_args, **_kwargs: output
        validated, observed_edges = assembler._validate_lane(
            output, self.metadata_path, CONTRACT, require_merge_input=True
        )
        self.assertEqual(validated["lane_id"], output["lane_id"])
        self.assertEqual(observed_edges, [10.0, 12.0, 14.0])

    def test_photon_lane_links_candidate_parity(self) -> None:
        _, _, reader = self.make_inputs("photon")
        output = self.run_extract(reader)
        self.assertEqual(len(output["candidate_parity_evidence"]), 1)
        self.assertEqual(
            output["candidate_parity_evidence"][0]["role"], "candidate_rows"
        )

    def test_photon_lane_without_candidate_parity_fails(self) -> None:
        _, _, reader = self.make_inputs("photon", include_parity=False)
        with self.assertRaisesRegex(EXTRACTOR.ExtractionError, "roles must be exact"):
            self.run_extract(reader)

    def test_additional_complete_group_is_accepted_and_hash_bound(self) -> None:
        metadata, _, reader = self.make_inputs()
        event_rows = [
            f"NONE /x/js_pp200_signal/g4hits/file{index}.root "
            f"/x/js_pp200_signal/nopileup/jets/file{index}.root NONE NONE"
            for index in range(1, 11)
        ]
        self.evidence_files["source_list"].write_text("\n".join(event_rows) + "\n")
        metadata["evidence"]["source_list"] = {
            "path": str(self.evidence_files["source_list"]),
            "sha256": file_sha256(self.evidence_files["source_list"]),
        }
        self.evidence_files["event_set"].write_text("\n".join(event_rows) + "\n")
        metadata["evidence"]["event_set"] = {
            "path": str(self.evidence_files["event_set"]),
            "sha256": file_sha256(self.evidence_files["event_set"]),
        }
        environment = EXECUTION._expected_environment(
            {
                "lane_id": metadata["lane_id"],
                "family": metadata["family"],
                "sample": metadata["sample"],
                "period": metadata["period"],
                "interaction": metadata["interaction"],
            },
            CONTRACT,
        )
        write_json(
            self.execution_path,
            EXECUTION._build_receipt(
                metadata["lane_id"],
                self.evidence_files["source_list"],
                environment,
                CONTRACT_PATH,
            ),
        )
        metadata["execution_contract"] = {
            "path": str(self.execution_path),
            "sha256": file_sha256(self.execution_path),
        }
        groups = EXTRACTOR._canonical_groups(metadata["lane_id"], event_rows, 5)
        component_links = []
        for group in groups:
            component_path = self.root / f"component_{group['group_index']}.json"
            component_path.write_text("{}\n")
            component_links.append(
                {
                    "group_index": group["group_index"],
                    "group_id": group["group_id"],
                    "path": str(component_path),
                    "sha256": file_sha256(component_path),
                }
            )
        metadata.update(
            {
                "group_count": 2,
                "event_set_row_count": 10,
                "groups": groups,
                "group_set_sha256": EXTRACTOR._payload_sha256(groups),
                "group_input_sidecars": component_links,
            }
        )
        self.rewrite(metadata, self.metadata_path)
        output = self.run_extract(reader)
        self.assertEqual(output["group_count"], 2)
        self.assertEqual(output["event_set_row_count"], 10)
        self.assertEqual(output["groups"], groups)
        self.assertEqual(len(output["extraction"]["group_input_sidecars"]), 2)

    def test_partial_or_forged_group_evidence_fails(self) -> None:
        metadata, _, reader = self.make_inputs()
        metadata["groups"][0]["event_rows"].pop()
        self.rewrite(metadata, self.metadata_path)
        with self.assertRaisesRegex(EXTRACTOR.ExtractionError, "identities/hashes"):
            self.run_extract(reader)

    def test_nonunit_scale_fails(self) -> None:
        metadata, _, reader = self.make_inputs()
        metadata["external_scale"] = 2.3398
        self.rewrite(metadata, self.metadata_path)
        with self.assertRaisesRegex(EXTRACTOR.ExtractionError, "non-unit"):
            self.run_extract(reader)

    def test_lane_identity_mismatch_fails(self) -> None:
        metadata, _, reader = self.make_inputs()
        metadata["lane_id"] = "inclusive:jet12:0mrad:si"
        self.rewrite(metadata, self.metadata_path)
        with self.assertRaisesRegex(EXTRACTOR.ExtractionError, "lane identity mismatch"):
            self.run_extract(reader)

    def test_root_hash_mismatch_fails_before_read(self) -> None:
        metadata, _, reader = self.make_inputs()
        metadata["root"]["sha256"] = "0" * 64
        self.rewrite(metadata, self.metadata_path)
        with self.assertRaisesRegex(EXTRACTOR.ExtractionError, "root hash mismatch"):
            self.run_extract(reader)

    def test_empty_evidence_file_fails(self) -> None:
        metadata, _, reader = self.make_inputs()
        self.evidence_files["config"].write_bytes(b"")
        metadata["evidence"]["config"]["sha256"] = file_sha256(
            self.evidence_files["config"]
        )
        self.rewrite(metadata, self.metadata_path)
        with self.assertRaisesRegex(EXTRACTOR.ExtractionError, "evidence.config is empty"):
            self.run_extract(reader)

    def test_event_set_must_be_exact_ordered_source_list_slice(self) -> None:
        metadata, _, reader = self.make_inputs()
        reordered = [self.event_rows[1], self.event_rows[0], *self.event_rows[2:]]
        self.evidence_files["event_set"].write_text("\n".join(reordered) + "\n")
        metadata["evidence"]["event_set"] = {
            "path": str(self.evidence_files["event_set"]),
            "sha256": file_sha256(self.evidence_files["event_set"]),
        }
        event_rows = reordered
        groups = EXTRACTOR._canonical_groups(metadata["lane_id"], event_rows, 5)
        metadata["groups"] = groups
        metadata["group_set_sha256"] = EXTRACTOR._payload_sha256(groups)
        self.rewrite(metadata, self.metadata_path)
        with self.assertRaisesRegex(EXTRACTOR.ExtractionError, "production-list slice"):
            self.run_extract(reader)

    def test_runtime_hashes_are_recomputed_from_required_executed_roles(self) -> None:
        metadata, _, reader = self.make_inputs()
        runtime = json.loads(self.runtime_manifest.read_text())
        target = next(
            row for row in runtime["files"] if row["role"] == "ppg_apply_model_base_v3E"
        )
        Path(target["path"]).write_text("mutated model bytes\n")
        # A caller cannot rescue stale runtime evidence by changing only a lane hash.
        metadata["model_set_sha256"] = "0" * 64
        self.rewrite(metadata, self.metadata_path)
        with self.assertRaisesRegex(EXTRACTOR.ExtractionError, "hash mismatch"):
            self.run_extract(reader)

    def test_runtime_manifest_missing_required_role_fails_closed(self) -> None:
        metadata, _, reader = self.make_inputs()
        runtime = json.loads(self.runtime_manifest.read_text())
        runtime["files"] = [
            row for row in runtime["files"] if row["role"] != "ppg_apply_model_base_E"
        ]
        write_json(self.runtime_manifest, runtime)
        metadata["runtime_manifest"]["sha256"] = file_sha256(self.runtime_manifest)
        self.rewrite(metadata, self.metadata_path)
        with self.assertRaisesRegex(EXTRACTOR.ExtractionError, "lacks roles"):
            self.run_extract(reader)

    def test_missing_unsuffixed_object_fails_without_suffixed_fallback(self) -> None:
        _, _, reader = self.make_inputs()
        del reader.histograms["SIM/h_tight_iso_cluster_0"]
        with self.assertRaisesRegex(EXTRACTOR.ExtractionError, "class-suffixed substitutes"):
            self.run_extract(reader)

    def test_ambiguous_auto_object_resolution_fails(self) -> None:
        metadata, _, reader = self.make_inputs()
        metadata["object_prefix"] = "auto"
        self.rewrite(metadata, self.metadata_path)
        reader.extra_paths["h_tight_iso_cluster_0"] = ["ALT/h_tight_iso_cluster_0"]
        reader.cycles["ALT/h_tight_iso_cluster_0"] = 1
        with self.assertRaisesRegex(EXTRACTOR.ExtractionError, "ambiguous ROOT object"):
            self.run_extract(reader)

    def test_multiple_root_cycles_fail(self) -> None:
        _, _, reader = self.make_inputs()
        reader.cycles["SIM/h_tight_iso_cluster_0"] = 2
        with self.assertRaisesRegex(EXTRACTOR.ExtractionError, "ambiguous ROOT cycles"):
            self.run_extract(reader)

    def test_missing_sumw2_fails(self) -> None:
        _, _, reader = self.make_inputs()
        reader.histograms["SIM/h_tight_iso_cluster_0"]["sumw2"] = []
        with self.assertRaisesRegex(EXTRACTOR.ExtractionError, "sumw2"):
            self.run_extract(reader)

    def test_missing_fill_evidence_fails(self) -> None:
        metadata, _, reader = self.make_inputs()
        del metadata["fill_evidence"]
        self.rewrite(metadata, self.metadata_path)
        with self.assertRaisesRegex(EXTRACTOR.ExtractionError, "fill_evidence"):
            self.run_extract(reader)

    def test_fill_lane_mismatch_fails(self) -> None:
        metadata, fill_payload, reader = self.make_inputs()
        fill_payload["lane_id"] = "inclusive:jet12:0mrad:si"
        self.rewrite(fill_payload, self.fill_path)
        self.refresh_fill_link(metadata)
        with self.assertRaisesRegex(EXTRACTOR.ExtractionError, "fill-evidence lane mismatch"):
            self.run_extract(reader)

    def test_fill_object_mismatch_fails(self) -> None:
        metadata, fill_payload, reader = self.make_inputs()
        fill_payload["observables"]["A"]["object"] = "SIM/h_tight_iso_cluster_signal_0"
        self.rewrite(fill_payload, self.fill_path)
        self.refresh_fill_link(metadata)
        with self.assertRaisesRegex(EXTRACTOR.ExtractionError, "does not bind exact ROOT object"):
            self.run_extract(reader)

    def test_fill_counts_must_close_entries(self) -> None:
        metadata, fill_payload, reader = self.make_inputs()
        fill_payload["observables"]["A"]["fills"] = [1, 3]
        self.rewrite(fill_payload, self.fill_path)
        self.refresh_fill_link(metadata)
        with self.assertRaisesRegex(EXTRACTOR.ExtractionError, "does not close TH1 entries"):
            self.run_extract(reader)

    def test_bad_root_reader_fails_closed(self) -> None:
        self.make_inputs()

        def bad_reader(_: Path):
            raise OSError("corrupt ROOT header")

        with self.assertRaisesRegex(EXTRACTOR.ExtractionError, "bad ROOT input"):
            EXTRACTOR.extract_lane(
                self.metadata_path, CONTRACT_PATH, reader_factory=bad_reader
            )

    def test_pyroot_written_file_smoke(self) -> None:
        try:
            import ROOT  # type: ignore
        except ImportError:
            self.skipTest("PyROOT unavailable")
        metadata, fill_payload, _ = self.make_inputs("photon")
        root_file = ROOT.TFile(str(self.root_file), "RECREATE")
        directory = root_file.mkdir("SIM")
        directory.cd()
        for name in self.required("photon"):
            histogram = ROOT.TH1D(
                EXTRACTOR.REGION_OBJECTS[name], "", 2, array("d", [10.0, 12.0, 14.0])
            )
            histogram.Sumw2()
            for _ in range(2):
                histogram.Fill(11.0, 0.5)
            for _ in range(3):
                histogram.Fill(13.0, 0.5)
            histogram.Write()
        padding = "".join(hashlib.sha256(str(i).encode()).hexdigest() for i in range(3000))
        ROOT.TObjString(padding).Write("terminal_padding")
        root_file.Close()
        metadata["root"] = {
            "path": str(self.root_file),
            "sha256": file_sha256(self.root_file),
        }
        fill_payload["root_sha256"] = file_sha256(self.root_file)
        self.rewrite(fill_payload, self.fill_path)
        metadata["fill_evidence"] = {
            "path": str(self.fill_path),
            "sha256": file_sha256(self.fill_path),
        }
        self.rewrite(metadata, self.metadata_path)
        output = EXTRACTOR.extract_lane(self.metadata_path, CONTRACT_PATH)
        self.assertEqual(output["observables"]["A_signal"]["sumw"], [1.0, 1.5])
        self.assertEqual(output["observables"]["A_signal"]["sumw2"], [0.5, 0.75])


if __name__ == "__main__":
    unittest.main()
