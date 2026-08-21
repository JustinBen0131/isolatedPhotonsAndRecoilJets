#!/usr/bin/env python3
"""Pure local regression tests for THE-134 H70 ML-side contracts."""

from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from types import SimpleNamespace
from pathlib import Path
from unittest import mock

import numpy as np


ROOT = Path(__file__).resolve().parents[4]
CONTRACTS = ROOT / "scripts" / "ml" / "contracts"
sys.path.insert(0, str(CONTRACTS))
import the134_h70_contract as contract  # noqa: E402


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


prepare = load(
    "prepare_the134_h70_matrix",
    ROOT / "scripts" / "ml" / "training" / "prepare_the134_h70_matrix.py",
)
runner = load(
    "run_the134_h70_model",
    ROOT / "scripts" / "ml" / "training" / "run_the134_h70_model.py",
)
working_points = load(
    "derive_the134_h70_wp",
    ROOT / "scripts" / "ml" / "working_points" / "derive_the134_h70_wp.py",
)
holdout = load(
    "materialize_the134_h70_holdout",
    ROOT / "scripts" / "ml" / "validation" / "materialize_the134_h70_holdout.py",
)
model_validation = load(
    "validate_the134_h70_model",
    ROOT / "scripts" / "ml" / "validation" / "validate_the134_h70_model.py",
)
reuse_audit = load(
    "audit_the134_model_reuse",
    ROOT / "scripts" / "ml" / "validation" / "audit_the134_model_reuse.py",
)
registry = load(
    "build_the134_h70_model_registry",
    ROOT / "scripts" / "ml" / "validation" / "build_the134_h70_model_registry.py",
)


def synthetic_full_extraction_audit(
    system: str, matrix_sha256: str, *, view_name: str = "H70"
) -> dict:
    sources = {}
    for source in sorted(contract.expected_sources(system)):
        checks = {
            "authority_state": True,
            "full_source_manifest_sha256": True,
            "all_inputs_bound_to_full_manifest": True,
            "expected_input_count_positive": True,
            "input_count_exact": True,
            "expected_occurrence_count_positive": True,
            "occurrence_count_exact": True,
            "input_records_sha256": True,
            "input_record_set_exact": True,
            "record_system_exact": True,
            "record_source_exact": True,
            "training_view_root_hashes": True,
        }
        sources[source] = {
            "status": "PASS",
            "state": "SOURCE_COMPLETE",
            "full_source_manifest_sha256": "a" * 64,
            "input_records_sha256": "b" * 64,
            "expected_input_count": 1,
            "observed_input_count": 1,
            "expected_occurrence_count": 1,
            "observed_occurrence_count": 1,
            "selected_view_training_rows": 1,
            "checks": checks,
        }
    closure = {
        "schema": "THE134_EXACT_SOURCE_POPULATION_CLOSURE_V1",
        "status": "PASS",
        "sources": sources,
        "semantic_sha256": contract.canonical_json_sha256(sources),
    }
    return {
        "schema": "THE134_FACTORIAL_VIEW_TRAINING_MATRIX_AUDIT_V1",
        "status": "PASS",
        "system": system,
        "scope": "full",
        "full_training_authority": 1,
        "shower_definition": view_name,
        "shower_semantic_sha256": contract.shower_semantic_sha256(view_name),
        "feature_order": list(contract.FEATURES_BY_SYSTEM[system]),
        "required_sources": sorted(contract.expected_sources(system)),
        "observed_manifest_sources": sorted(contract.expected_sources(system)),
        "source_provenance_json_sha256": "c" * 64,
        "source_population_closure": closure,
        "matrix_sha256": matrix_sha256,
    }


def synthetic_arrays(
    system: str, *, auau_discard: bool = False, cluster_et: float = 20.0
) -> dict[str, np.ndarray]:
    features = contract.FEATURES_BY_SYSTEM[system]
    vector = np.arange(len(features), dtype=np.float32) / 10.0
    vector[features.index("cluster_Et")] = cluster_et
    vector[features.index("cluster_Eta")] = 0.2
    vector[features.index("vertexz")] = -3.0
    if system == "auau":
        vector[features.index("centrality")] = 25.0
    n = len(contract.ALL_SHOWER_VIEWS)
    labels = np.full(n, -1 if auau_discard else 1, dtype=np.int32)
    in_domain = 15.0 <= cluster_et < 35.0
    source = "run28_photonjet20" if system == "pp" else "run28_embeddedPhoton12"
    if system == "pp":
        weight_slice, weight_cross_section = 2.0, 2.0
        weight_vertex, weight_si_di = 3.0, 0.5
        weight_period, weight_exposure = 4.0, 4.0
    else:
        weight_slice, weight_cross_section = 2.0, 0.5
        weight_vertex, weight_si_di = 3.0, 0.25
        weight_period, weight_exposure = 4.0, 4.0
    weight_final = (
        weight_slice * weight_vertex * weight_si_di * weight_period
        if system == "pp"
        else weight_slice
        * weight_cross_section
        * weight_vertex
        * weight_si_di
        * weight_period
        * weight_exposure
    )
    source_id = contract.identity128_from_text(
        "|".join(
            [
                "TRAINING",
                "SYNTHETIC",
                source,
                "RUN28",
                "1",
                "0",
                "b" * 64,
                "c" * 64,
                "a" * 64,
            ]
        )
    )
    event_id = contract.identity128_from_text(f"TRAINING|{source}|1|0|1")
    candidate_id = contract.identity128_from_text(
        f"{contract.identity128_hex(event_id)}|candidate|0|7"
    )
    definition_ids = [
        contract.identity128_from_text(
            f"shower-definition|{view}|{contract.shower_semantic_text(view)}"
        )
        for view in contract.ALL_SHOWER_VIEWS
    ]
    training_view_ids = [
        contract.identity128_from_text(
            f"{contract.identity128_hex(candidate_id)}|training-view|"
            f"{contract.identity128_hex(definition_id)}"
        )
        for definition_id in definition_ids
    ]
    arrays = {
        "definition_name": np.asarray(contract.ALL_SHOWER_VIEWS, dtype=object),
        "shower_semantic_sha256": np.asarray(
            [contract.shower_semantic_sha256(view) for view in contract.ALL_SHOWER_VIEWS],
            dtype=object,
        ),
        "feature_contract_sha256": np.full(
            n, contract.feature_contract_sha256(system), dtype=object
        ),
        "ordered_features": np.asarray([vector.copy() for _ in range(n)], dtype=object),
        "candidate_id_hi": np.full(n, candidate_id[0], dtype=np.uint64),
        "candidate_id_lo": np.full(n, candidate_id[1], dtype=np.uint64),
        "training_view_id_hi": np.asarray([item[0] for item in training_view_ids], dtype=np.uint64),
        "training_view_id_lo": np.asarray([item[1] for item in training_view_ids], dtype=np.uint64),
        "source_occurrence_id_hi": np.full(n, source_id[0], dtype=np.uint64),
        "source_occurrence_id_lo": np.full(n, source_id[1], dtype=np.uint64),
        "event_id_hi": np.full(n, event_id[0], dtype=np.uint64),
        "event_id_lo": np.full(n, event_id[1], dtype=np.uint64),
        "definition_id_hi": np.asarray([item[0] for item in definition_ids], dtype=np.uint64),
        "definition_id_lo": np.asarray([item[1] for item in definition_ids], dtype=np.uint64),
        "feature_count": np.full(n, len(features), dtype=np.int32),
        "finite_feature_state": np.ones(n, dtype=np.int32),
        "system_code": np.full(n, contract.SYSTEM_CODES[system], dtype=np.int32),
        "source_sample": np.full(n, source, dtype=object),
        "source_lane": np.full(n, "TRAINING", dtype=object),
        "source_dataset": np.full(n, "SYNTHETIC", dtype=object),
        "source_period": np.full(n, "RUN28", dtype=object),
        "source_si_di_role": np.full(n, "SI", dtype=object),
        "source_ownership_state": np.full(n, "OWNED", dtype=object),
        "source_manifest_sha256": np.full(n, "a" * 64, dtype=object),
        "input_uri_sha256": np.full(n, "b" * 64, dtype=object),
        "input_file_sha256": np.full(n, "c" * 64, dtype=object),
        "source_run": np.ones(n, dtype=np.int32),
        "source_segment": np.zeros(n, dtype=np.int32),
        "training_label": labels,
        "is_signal": np.ones(n, dtype=np.int32),
        "source_role": np.ones(n, dtype=np.int32),
        "source_sample_code": np.full(
            n, 20 if system == "pp" else 12, dtype=np.int32
        ),
        "label_authority": np.full(n, "PPG12_SOURCE_ROLE", dtype=object),
        "model_domain_state": np.full(n, 0 if in_domain else 1, dtype=np.int32),
        "below15_retention_state": np.full(n, int(cluster_et < 15.0), dtype=np.int32),
        "nominal_training_eligible": ((labels >= 0) & in_domain).astype(np.int32),
        "working_point_state": np.full(n, -1, dtype=np.int32),
        "tag_state": np.full(n, -1, dtype=np.int32),
        "ppg12_source_role_label": labels.copy(),
        "truth_match_found": np.ones(n, dtype=np.int32),
        "truth_photon_class": np.ones(n, dtype=np.int32),
        "truth_is_prompt": np.ones(n, dtype=np.int32),
        "truth_iso_et": np.full(n, 1.25, dtype=np.float64),
        "truth_iso_pass": np.ones(n, dtype=np.int32),
        "cluster_truth_track_id": np.full(n, 101, dtype=np.int32),
        "cluster_truth_pid": np.full(n, 22, dtype=np.int32),
        "cluster_truth_barcode": np.full(n, 202, dtype=np.int32),
        "truth_energy_contribution": np.full(n, 0.75, dtype=np.float64),
        "npb_label": np.full(n, -1, dtype=np.int32),
        "is_npb": np.zeros(n, dtype=np.int32),
        "ppg12_analysis_window_pass": np.ones(n, dtype=np.int32),
        "ppg12_response_window_pass": np.ones(n, dtype=np.int32),
        "ppg12_truth_window_pass_r04": np.ones(n, dtype=np.int32),
        "ppg12_sample_bin": np.full(n, 20 if system == "pp" else 12, dtype=np.int32),
        "ppg12_xsec_pb": np.full(n, 2.5, dtype=np.float64),
        "ppg12_xsec_weight": np.full(n, 0.4, dtype=np.float64),
        "ppg12_window_low": np.full(n, 15.0, dtype=np.float64),
        "ppg12_window_high": np.full(n, 35.0, dtype=np.float64),
        "max_truth_jet_pt_r04": np.full(n, 8.0, dtype=np.float64),
        "minimum_bias_classifier_decision": np.full(
            n, 2 if system == "auau" else -1, dtype=np.int32
        ),
        "weight_application_count": np.ones(n, dtype=np.int32),
        "event_weight": np.full(n, weight_final, dtype=np.float64),
        "weight_slice": np.full(n, weight_slice, dtype=np.float64),
        "weight_cross_section": np.full(
            n, weight_cross_section, dtype=np.float64
        ),
        "weight_vertex": np.full(n, weight_vertex, dtype=np.float64),
        "weight_si_di": np.full(n, weight_si_di, dtype=np.float64),
        "weight_period": np.full(n, weight_period, dtype=np.float64),
        "weight_exposure": np.full(n, weight_exposure, dtype=np.float64),
        "weight_final": np.full(n, weight_final, dtype=np.float64),
        "cluster_Et": np.full(n, cluster_et, dtype=np.float32),
        "cluster_Eta": np.full(n, 0.2, dtype=np.float32),
        "cluster_Phi": np.full(n, 0.4, dtype=np.float32),
        "vertexz": np.full(n, -3.0, dtype=np.float32),
        "centrality": np.full(n, 25.0 if system == "auau" else -1.0, dtype=np.float32),
        "event_sequence": np.ones(n, dtype=np.int64),
        "run": np.ones(n, dtype=np.int32),
        "encounter_ordinal": np.zeros(n, dtype=np.int32),
        "cluster_index": np.zeros(n, dtype=np.int32),
        "cluster_map_key": np.full(n, 7, dtype=np.uint64),
    }
    return arrays


def synthetic_root_metadata(arrays: dict[str, np.ndarray]) -> dict[str, object]:
    return {
        "tree_num_entries": len(arrays["definition_name"]),
        "rj_photon_training_schema": contract.TRAINING_SCHEMA_NAME,
        "rj_photon_training_schema_version": contract.TRAINING_SCHEMA_VERSION,
        "schema_sha256": contract.training_schema_sha256(),
        "pp_feature_contract_sha256": contract.feature_contract_sha256("pp"),
        "auau_feature_contract_sha256": contract.feature_contract_sha256("auau"),
        "source_manifest_sha256": "a" * 64,
        "config_sha256": "d" * 64,
        "code_sha256": "e" * 64,
        "rj_photon_training_complete": "1",
        "rj_photon_training_entries": str(len(arrays["definition_name"])),
    }


def synthetic_source_provenance_record(
    path: Path,
    *,
    row_id: str = "pp_run28_photonjet20_000001",
    system: str = "pp",
    source: str = "run28_photonjet20",
) -> dict[str, str]:
    return {
        "path": str(path),
        "row_id": row_id,
        "system": system,
        "source_sample": source,
        "source_occurrence_id_hex": "0" * 31 + "1",
        "input_uri_sha256": "1" * 64,
        "input_file_sha256": "1" * 64,
        "source_manifest_sha256": "3" * 64,
        "config_sha256": "4" * 64,
        "code_sha256": "5" * 64,
        "training_view_root_sha256": "6" * 64,
    }


class ContractTests(unittest.TestCase):
    def test_factorial_matrix_reads_each_sidecar_once_and_emits_seven_audits(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_dir = root / "run28_photonjet20"
            source_dir.mkdir()
            inputs = [source_dir / "one.root", source_dir / "two.root"]
            for path in inputs:
                path.write_bytes(b"placeholder")

            signal = synthetic_arrays("pp")
            background = synthetic_arrays("pp")
            background["training_label"][:] = 0
            background["is_signal"][:] = 0
            background["ppg12_source_role_label"][:] = 0
            background["candidate_id_lo"][:] ^= np.uint64(1)

            def validation_side_effect(arrays, _metadata, **kwargs):
                view_name = kwargs["view_name"]
                selected = np.asarray(arrays["definition_name"] == view_name)
                vectors = [
                    list(map(float, arrays["ordered_features"][index]))
                    for index in np.flatnonzero(selected)
                ]
                return selected, vectors, {
                    "status": "PASS",
                    "population_state": "POPULATED",
                    "rows": len(arrays["definition_name"]),
                    "candidates": 1,
                    "selected_view": view_name,
                    "selected_view_rows": 1,
                    "selected_view_training_rows": 1,
                    "selected_view_input_invalid_rows": 0,
                    "discarded_selected_view_label_minus_one_rows": 0,
                    "below15_selected_view_rows": 0,
                    "below15_selected_view_eligible_rows": 0,
                    "out_of_domain_selected_view_rows": 0,
                    "out_of_domain_selected_view_eligible_rows": 0,
                    "definition_counts": {view: 1 for view in contract.ALL_SHOWER_VIEWS},
                    "source": kwargs["source"],
                    "source_occurrence_ids": ["0" * 31 + "1"],
                    "expected_label": 1,
                    "failures": [],
                }

            matrix = root / "matrix.npz"
            anchor_audit = root / "h70-audit.json"
            factorial_dir = root / "factorial-audits"
            args = SimpleNamespace(
                system="pp",
                view="H70",
                input=inputs,
                tree=contract.TREE_NAME,
                matrix_out=matrix,
                audit_out=anchor_audit,
                factorial_audit_directory=factorial_dir,
                scope="smoke",
                skip_input_hashes=True,
                source_provenance_json=None,
            )
            with mock.patch.object(prepare, "parse_args", return_value=args), mock.patch.object(
                prepare,
                "read_tree",
                side_effect=[
                    (signal, list(contract.CORE_BRANCHES), synthetic_root_metadata(signal)),
                    (
                        background,
                        list(contract.CORE_BRANCHES),
                        synthetic_root_metadata(background),
                    ),
                ],
            ) as read_tree, mock.patch.object(
                prepare, "validate_matrix_input", side_effect=validation_side_effect
            ) as validate:
                self.assertEqual(prepare.main(), 0)

            self.assertEqual(read_tree.call_count, 2)
            self.assertEqual(validate.call_count, 2 * len(contract.ALL_SHOWER_VIEWS))
            shared_sha = contract.sha256_file(matrix)
            anchor = json.loads(anchor_audit.read_text())
            self.assertEqual(anchor["matrix_sha256"], shared_sha)
            self.assertEqual(
                anchor["single_read_factorial_projection"]["root_reads"], 2
            )
            manifest = json.loads(
                (factorial_dir / "factorial_matrix_manifest.json").read_text()
            )
            self.assertEqual(set(manifest["audits"]), set(contract.ALL_SHOWER_VIEWS))
            for view_name, record in manifest["audits"].items():
                audit = json.loads(Path(record["path"]).read_text())
                self.assertEqual(audit["shower_definition"], view_name)
                self.assertEqual(audit["matrix_sha256"], shared_sha)

    def test_valid_empty_input_uses_external_occurrence_without_rows(self):
        arrays = {
            name: values[:0]
            for name, values in synthetic_arrays("pp").items()
        }
        metadata = synthetic_root_metadata(arrays)
        path = Path("/frozen/digest/valid-empty.root")
        record = synthetic_source_provenance_record(path)
        record.update(
            {
                "source_manifest_sha256": "a" * 64,
                "config_sha256": "d" * 64,
                "code_sha256": "e" * 64,
            }
        )
        selected, vectors, report = prepare.validate_matrix_input(
            arrays,
            metadata,
            system="pp",
            source="run28_photonjet20",
            view_name="H70",
            external=record,
        )
        self.assertEqual(report["status"], "PASS")
        self.assertEqual(report["population_state"], "VALID_EMPTY")
        self.assertEqual(
            report["source_occurrence_ids"],
            [record["source_occurrence_id_hex"]],
        )
        self.assertEqual(selected.size, 0)
        self.assertEqual(vectors, [])

    def test_nonempty_input_binds_rows_to_external_occurrence(self):
        arrays = synthetic_arrays("pp")
        arrays["input_file_sha256"][:] = "b" * 64
        source_id = contract.identity128_from_text(
            "|".join(
                [
                    "TRAINING",
                    "SYNTHETIC",
                    "run28_photonjet20",
                    "RUN28",
                    "1",
                    "0",
                    "b" * 64,
                    "b" * 64,
                    "a" * 64,
                ]
            )
        )
        arrays["source_occurrence_id_hi"][:] = source_id[0]
        arrays["source_occurrence_id_lo"][:] = source_id[1]
        metadata = synthetic_root_metadata(arrays)
        record = synthetic_source_provenance_record(
            Path("/frozen/digest/populated.root")
        )
        record.update(
            {
                "source_occurrence_id_hex": contract.identity128_hex(source_id),
                "input_uri_sha256": "b" * 64,
                "input_file_sha256": "b" * 64,
                "source_manifest_sha256": "a" * 64,
                "config_sha256": "d" * 64,
                "code_sha256": "e" * 64,
            }
        )
        _selected, _vectors, report = prepare.validate_matrix_input(
            arrays,
            metadata,
            system="pp",
            source="run28_photonjet20",
            view_name="H70",
            external=record,
        )
        self.assertEqual(report["status"], "PASS")
        bad_record = dict(record)
        bad_record["source_occurrence_id_hex"] = "f" * 32
        _selected, _vectors, bad = prepare.validate_matrix_input(
            arrays,
            metadata,
            system="pp",
            source="run28_photonjet20",
            view_name="H70",
            external=bad_record,
        )
        self.assertEqual(bad["status"], "FAIL")
        self.assertIn(
            "external/row source occurrence identity mismatch",
            bad["failures"],
        )

    def test_source_input_record_hash_binds_occurrence_identity(self):
        record = synthetic_source_provenance_record(
            Path("/frozen/digest/hash.root")
        )
        first = prepare.source_input_records_sha256([record])
        record["source_occurrence_id_hex"] = "f" * 32
        second = prepare.source_input_records_sha256([record])
        self.assertNotEqual(first, second)

    def test_valid_empty_occurrence_counts_toward_exact_source_closure(self):
        system = "pp"
        paths = []
        path_sources = []
        provenance = {}
        authority = {}
        reports = []
        selected_rows = {}
        for index, source in enumerate(contract.expected_sources(system)):
            path = Path(f"/frozen/{source}/{index}.root")
            occurrence = f"{index + 1:032x}"
            record = synthetic_source_provenance_record(
                path,
                row_id=f"{source}_{index:06d}",
                system=system,
                source=source,
            )
            record["source_occurrence_id_hex"] = occurrence
            record["source_manifest_sha256"] = f"{index + 1:064x}"
            paths.append(path)
            path_sources.append(source)
            provenance[str(path)] = record
            reports.append(
                {
                    "status": "PASS",
                    "population_state": (
                        "VALID_EMPTY" if index == 0 else "POPULATED"
                    ),
                    "source": source,
                    "source_occurrence_ids": [occurrence],
                }
            )
            selected_rows[source] = 0 if index == 0 else 1
            authority[source] = {
                "state": (
                    "SOURCE_COMPLETE_ZERO_IN_DOMAIN"
                    if index == 0
                    else "SOURCE_COMPLETE"
                ),
                "full_source_manifest_sha256": record[
                    "source_manifest_sha256"
                ],
                "expected_input_count": 1,
                "expected_occurrence_count": 1,
                "input_records_sha256": prepare.source_input_records_sha256(
                    [record]
                ),
            }
        closure = prepare.build_source_population_closure(
            system=system,
            paths=paths,
            path_sources=path_sources,
            external_provenance=provenance,
            source_coverage_authority=authority,
            file_reports=reports,
            selected_rows_by_source=selected_rows,
        )
        self.assertEqual(closure["status"], "PASS")
        first_source = contract.expected_sources(system)[0]
        self.assertEqual(
            closure["sources"][first_source]["observed_occurrence_count"],
            1,
        )

    def test_explicit_provenance_routes_opaque_path_without_token_inference(self):
        path = Path("/frozen/digest/0123456789abcdef.root")
        record = synthetic_source_provenance_record(path)
        with mock.patch.object(
            prepare,
            "source_from_path",
            side_effect=AssertionError("path inference must not run"),
        ):
            sources, failures = prepare.resolve_path_sources(
                [path], "pp", {str(path): record}
            )
        self.assertEqual(sources, ["run28_photonjet20"])
        self.assertEqual(failures, [])

    def test_explicit_provenance_rejects_wrong_system_without_fallback(self):
        path = Path("/frozen/digest/0123456789abcdef.root")
        record = synthetic_source_provenance_record(
            path,
            row_id="auau_run28_embeddedPhoton12_000001",
            system="auau",
            source="run28_embeddedPhoton12",
        )
        with mock.patch.object(
            prepare,
            "source_from_path",
            side_effect=AssertionError("path inference must not run"),
        ):
            with self.assertRaisesRegex(SystemExit, "system mismatch"):
                prepare.resolve_path_sources([path], "pp", {str(path): record})

    def test_source_provenance_manifest_rejects_unsafe_hash_and_duplicate_paths(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = Path("/frozen/digest/first.root")
            manifest = root / "provenance.json"
            payload = {
                "schema": "THE134_SOURCE_PROVENANCE_V1",
                "inputs": [synthetic_source_provenance_record(first)],
                "source_coverage_authority": {},
            }
            manifest.write_text(json.dumps(payload))
            indexed, coverage = prepare.load_source_provenance(manifest)
            self.assertEqual(
                indexed[str(first)]["source_sample"], "run28_photonjet20"
            )
            self.assertEqual(coverage, {})

            payload["inputs"][0]["path"] = "/frozen/../escape.root"
            manifest.write_text(json.dumps(payload))
            with self.assertRaisesRegex(SystemExit, "not canonical absolute"):
                prepare.load_source_provenance(manifest)

            payload["inputs"] = [synthetic_source_provenance_record(first)]
            payload["inputs"][0]["code_sha256"] = "not-a-hash"
            manifest.write_text(json.dumps(payload))
            with self.assertRaisesRegex(SystemExit, "code_sha256 is not SHA-256"):
                prepare.load_source_provenance(manifest)

            payload["inputs"] = [synthetic_source_provenance_record(first)]
            payload["inputs"][0]["source_occurrence_id_hex"] = "A" * 32
            manifest.write_text(json.dumps(payload))
            with self.assertRaisesRegex(
                SystemExit, "source_occurrence_id_hex is invalid"
            ):
                prepare.load_source_provenance(manifest)

            payload["inputs"] = [
                synthetic_source_provenance_record(first),
                synthetic_source_provenance_record(
                    first,
                    row_id="pp_run28_photonjet20_000002",
                ),
            ]
            manifest.write_text(json.dumps(payload))
            with self.assertRaisesRegex(
                SystemExit, "duplicate source provenance path"
            ):
                prepare.load_source_provenance(manifest)

    def test_legacy_smoke_without_manifest_retains_path_source_resolution(self):
        path = Path("/frozen/run28_photonjet20/input.root")
        sources, failures = prepare.resolve_path_sources([path], "pp", {})
        self.assertEqual(sources, ["run28_photonjet20"])
        self.assertEqual(failures, [])

    def test_h70_semantic_hash_matches_frozen_cpp_contract(self):
        self.assertEqual(
            contract.shower_semantic_sha256(),
            "fe6633650503f23cf3cfa662a2884d29efbcd7ca465b5b02cfe31ec220db3987",
        )

    def test_h70_definition_identity_matches_fixed_cpp_vector(self):
        identity = contract.identity128_from_text(
            "shower-definition|H70|" + contract.shower_semantic_text("H70")
        )
        self.assertEqual(identity, (1467324610346303874, 2533097968370047012))
        self.assertEqual(contract.identity128_hex(identity), "145cfbf579e1418223275f2e69f5d824")

    def test_feature_contract_hashes_are_system_specific(self):
        self.assertNotEqual(
            contract.feature_contract_sha256("pp"),
            contract.feature_contract_sha256("auau"),
        )
        self.assertEqual(len(contract.feature_contract_sha256("pp")), 64)

    def test_root_metadata_completion_and_schema_are_load_bearing(self):
        self.assertEqual(
            contract.training_schema_sha256(),
            "69842a307688422682067ea58424b88f260895b427db94f2ffc5872c1ad4a549",
        )
        arrays = synthetic_arrays("pp")
        metadata = synthetic_root_metadata(arrays)
        self.assertEqual(
            prepare.validate_artifact_metadata(metadata, arrays, external=None), []
        )
        metadata["rj_photon_training_complete"] = "0"
        failures = prepare.validate_artifact_metadata(metadata, arrays, external=None)
        self.assertTrue(any("training_complete" in item for item in failures))

    def test_authoritative_auau_order_rejects_observed_pp11_append_order(self):
        wrong = [
            *contract.PP_FEATURES,
            "centrality",
            "cluster_weta33_cogx",
            "cluster_wphi33_cogx",
        ]
        metadata = {
            "features": wrong,
            "pt_range": list(contract.MODEL_DOMAIN_GEV),
            "weight_mode": "ppg12-exact",
            "split": {"mode": "event50"},
            "the134_view_contract": {
                "shower_definition": "H70",
                "shower_semantic_sha256": contract.shower_semantic_sha256(),
            },
        }
        with self.assertRaises(SystemExit):
            working_points.require_model_metadata(metadata, "auau")
        self.assertNotEqual(wrong, list(contract.AUAU_FEATURES))

    def test_pp_factorial_and_order_pass(self):
        arrays = synthetic_arrays("pp")
        selected, vectors, report = prepare.validate_file_arrays(
            arrays, system="pp", source="run28_photonjet20"
        )
        self.assertEqual(report["status"], "PASS")
        self.assertEqual(int(np.sum(selected)), 1)
        self.assertEqual(len(vectors[0]), 11)

    def test_auau_minus_one_is_provenance_only(self):
        arrays = synthetic_arrays("auau", auau_discard=True)
        selected, vectors, report = prepare.validate_file_arrays(
            arrays, system="auau", source="run28_embeddedPhoton12"
        )
        self.assertEqual(report["status"], "PASS")
        self.assertEqual(int(np.sum(selected)), 0)
        self.assertEqual(vectors, [])
        self.assertEqual(report["discarded_selected_view_label_minus_one_rows"], 1)

    def test_missing_factorial_view_fails(self):
        arrays = synthetic_arrays("pp")
        arrays = {name: values[:-1] for name, values in arrays.items()}
        _selected, _vectors, report = prepare.validate_file_arrays(
            arrays, system="pp", source="run28_photonjet20"
        )
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(any("factorial" in failure for failure in report["failures"]))

    def test_below15_is_retained_but_never_training_eligible(self):
        arrays = synthetic_arrays("pp", cluster_et=12.0)
        selected, vectors, report = prepare.validate_file_arrays(
            arrays, system="pp", source="run28_photonjet20"
        )
        self.assertEqual(report["status"], "PASS")
        self.assertEqual(int(np.sum(selected)), 0)
        self.assertEqual(vectors, [])
        self.assertEqual(report["below15_selected_view_rows"], 1)
        self.assertEqual(report["below15_selected_view_eligible_rows"], 0)

        arrays["nominal_training_eligible"][:] = 1
        _selected, _vectors, bad = prepare.validate_file_arrays(
            arrays, system="pp", source="run28_photonjet20"
        )
        self.assertEqual(bad["status"], "FAIL")
        self.assertTrue(any("below-15" in failure for failure in bad["failures"]))

    def test_input_invalid_loose_capture_is_retained_but_not_selected(self):
        arrays = synthetic_arrays("pp")
        h70 = list(contract.ALL_SHOWER_VIEWS).index("H70")
        vector = np.asarray(arrays["ordered_features"][h70], dtype=np.float32).copy()
        vector[1] = np.nan
        arrays["ordered_features"][h70] = vector
        arrays["finite_feature_state"][h70] = 0
        arrays["model_domain_state"][h70] = 3
        arrays["nominal_training_eligible"][h70] = 0
        selected, vectors, report = prepare.validate_file_arrays(
            arrays, system="pp", source="run28_photonjet20"
        )
        self.assertEqual(report["status"], "PASS", report["failures"])
        self.assertEqual(int(np.sum(selected)), 0)
        self.assertEqual(vectors, [])
        self.assertEqual(report["selected_view_input_invalid_rows"], 1)

    def test_auau_requires_minimum_bias_classifier_pass(self):
        arrays = synthetic_arrays("auau")
        arrays["minimum_bias_classifier_decision"][:] = 1
        _selected, _vectors, report = prepare.validate_file_arrays(
            arrays, system="auau", source="run28_embeddedPhoton12"
        )
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(any("MinimumBias" in failure for failure in report["failures"]))

    def test_label_authority_and_weight_application_are_load_bearing(self):
        arrays = synthetic_arrays("pp")
        arrays["label_authority"][0] = "SILENT_RELABEL"
        arrays["weight_application_count"][1] = 2
        _selected, _vectors, report = prepare.validate_file_arrays(
            arrays, system="pp", source="run28_photonjet20"
        )
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(any("label authority" in failure for failure in report["failures"]))
        self.assertTrue(any("exactly one" in failure for failure in report["failures"]))

    def test_pp_nonunit_weight_ledger_and_source_authority_close(self):
        arrays = synthetic_arrays("pp")
        _selected, _vectors, report = prepare.validate_file_arrays(
            arrays, system="pp", source="run28_photonjet20"
        )
        self.assertEqual(report["status"], "PASS", report["failures"])
        self.assertTrue(np.all(arrays["weight_final"] == 12.0))

        arrays["weight_cross_section"][0] = 7.0
        arrays["weight_final"][1] = 11.0
        arrays["event_weight"][2] = 10.0
        arrays["source_role"][3] = 2
        arrays["source_sample_code"][4] = 10
        _selected, _vectors, bad = prepare.validate_file_arrays(
            arrays, system="pp", source="run28_photonjet20"
        )
        self.assertEqual(bad["status"], "FAIL")
        self.assertTrue(any("slice/cross-section" in failure for failure in bad["failures"]))
        self.assertTrue(any("product/final" in failure for failure in bad["failures"]))
        self.assertTrue(any("event_weight/weight_final" in failure for failure in bad["failures"]))
        self.assertTrue(any("source_role/source authority" in failure for failure in bad["failures"]))
        self.assertTrue(any("source_sample_code/source authority" in failure for failure in bad["failures"]))

    def test_auau_independent_nonunit_weight_ledger_is_load_bearing(self):
        arrays = synthetic_arrays("auau")
        _selected, _vectors, report = prepare.validate_file_arrays(
            arrays, system="auau", source="run28_embeddedPhoton12"
        )
        self.assertEqual(report["status"], "PASS", report["failures"])
        self.assertTrue(np.all(arrays["weight_final"] == 12.0))

        arrays["weight_exposure"][0] = 3.0
        _selected, _vectors, bad = prepare.validate_file_arrays(
            arrays, system="auau", source="run28_embeddedPhoton12"
        )
        self.assertEqual(bad["status"], "FAIL")
        self.assertTrue(any("product/final" in failure for failure in bad["failures"]))

        arrays = synthetic_arrays("auau")
        arrays["weight_slice"][0] = np.nan
        _selected, _vectors, bad = prepare.validate_file_arrays(
            arrays, system="auau", source="run28_embeddedPhoton12"
        )
        self.assertEqual(bad["status"], "FAIL")
        self.assertTrue(any("nonfinite/nonpositive" in failure for failure in bad["failures"]))

    def test_identity_foreign_key_mutations_are_rejected(self):
        arrays = synthetic_arrays("pp")
        arrays["source_occurrence_id_hi"][:] = 0
        arrays["source_occurrence_id_lo"][:] = 0
        arrays["definition_id_lo"][0] ^= np.uint64(1)
        _selected, _vectors, report = prepare.validate_file_arrays(
            arrays, system="pp", source="run28_photonjet20"
        )
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(any("source-occurrence" in failure for failure in report["failures"]))
        self.assertTrue(any("definition identity" in failure for failure in report["failures"]))

    def test_cluster_map_key_mutation_breaks_candidate_identity(self):
        arrays = synthetic_arrays("pp")
        arrays["cluster_map_key"][:] = np.uint64(8)
        _selected, _vectors, report = prepare.validate_file_arrays(
            arrays, system="pp", source="run28_photonjet20"
        )
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(
            any("cluster-map-key" in failure for failure in report["failures"])
        )

    def test_cluster_index_must_equal_encounter_ordinal(self):
        arrays = synthetic_arrays("pp")
        arrays["cluster_index"][0] = 1
        _selected, _vectors, report = prepare.validate_file_arrays(
            arrays, system="pp", source="run28_photonjet20"
        )
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(
            any("cluster_index/encounter_ordinal" in failure for failure in report["failures"])
        )

    def test_single_view_label_and_weight_ledger_drift_are_rejected(self):
        arrays = synthetic_arrays("auau")
        arrays["training_label"][0] = -1
        arrays["ppg12_source_role_label"][0] = -1
        arrays["nominal_training_eligible"][0] = 0
        _selected, _vectors, label_report = prepare.validate_file_arrays(
            arrays, system="auau", source="run28_embeddedPhoton12"
        )
        self.assertEqual(label_report["status"], "FAIL")
        self.assertTrue(
            any(
                "candidate-invariant training_label" in failure
                for failure in label_report["failures"]
            )
        )

        arrays = synthetic_arrays("auau")
        arrays["weight_slice"][0] = 3.0
        arrays["weight_final"][0] = 18.0
        arrays["event_weight"][0] = 18.0
        _selected, _vectors, weight_report = prepare.validate_file_arrays(
            arrays, system="auau", source="run28_embeddedPhoton12"
        )
        self.assertEqual(weight_report["status"], "FAIL")
        self.assertTrue(
            any(
                "candidate-invariant weight_slice" in failure
                for failure in weight_report["failures"]
            )
        )

    def test_valid_empty_auau_file_is_not_misclassified_as_file_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_dir = root / "run28_embeddedPhoton12"
            source_dir.mkdir()
            input_path = source_dir / "valid_empty.root"
            input_path.write_bytes(b"placeholder")
            audit = root / "audit.json"
            matrix = root / "matrix.npz"
            args = SimpleNamespace(
                system="auau",
                view="H70",
                input=[input_path],
                tree=contract.TREE_NAME,
                matrix_out=matrix,
                audit_out=audit,
                scope="smoke",
                skip_input_hashes=True,
                source_provenance_json=None,
            )
            arrays = synthetic_arrays("auau", auau_discard=True)
            with mock.patch.object(prepare, "parse_args", return_value=args), mock.patch.object(
                prepare,
                "read_tree",
                return_value=(
                    arrays,
                    list(contract.CORE_BRANCHES),
                    synthetic_root_metadata(arrays),
                ),
            ):
                self.assertEqual(prepare.main(), 2)
            payload = __import__("json").loads(audit.read_text())
            self.assertEqual(payload["file_reports"][0]["status"], "PASS")
            self.assertEqual(
                payload["file_reports"][0]["selected_view_training_rows"], 0
            )
            self.assertTrue(any("no valid H70 rows" in item for item in payload["failures"]))

    def test_auau_split_is_90_10_while_pp_remains_event50(self):
        self.assertEqual(runner.SYSTEM_DEFAULTS["pp"]["test_size"], 0.5)
        self.assertEqual(runner.SYSTEM_DEFAULTS["auau"]["test_size"], 0.1)
        self.assertEqual(contract.expected_model_split("auau", "H70")["mode"], "event50")

    def test_reuse_split_contracts_cannot_be_cross_used(self):
        pp = {"split": {"mode": "event50", "test_fraction_requested": 0.5, "random_seed": 42, "train_rows": 20, "test_rows": 20, "train_events": 10, "test_events": 10}}
        auau = {"split": {"mode": "row", "test_fraction_requested": 0.1, "random_seed": 13, "train_rows": 90, "test_rows": 10}}
        self.assertTrue(all(reuse_audit.reuse_split_checks(pp, "pp", "H70").values()))
        self.assertTrue(all(reuse_audit.reuse_split_checks(auau, "auau", "H0").values()))
        self.assertFalse(all(reuse_audit.reuse_split_checks(pp, "auau", "H0").values()))
        self.assertFalse(all(reuse_audit.reuse_split_checks(auau, "pp", "H70").values()))

    def test_the111_row_control_metadata_does_not_invent_event_counts(self):
        expected_split = contract.expected_model_split("auau", "H0")
        metadata = {
            "features": list(contract.AUAU_FEATURES),
            "pt_range": list(contract.MODEL_DOMAIN_GEV),
            "weight_mode": "ppg12-exact",
            "split": {
                "mode": "row",
                "test_fraction_requested": 0.1,
                "random_seed": 13,
                "train_rows": 90,
                "test_rows": 10,
            },
            "accepted_split_contract": expected_split,
            "model_origin": "REUSED_THE111",
            "reuse_pinned_hashes": {
                **contract.expected_reuse_pinned_hashes("auau", "H0")
            },
            "weighting": {
                "event_weight_used": False,
                "vertex_reweight": False,
                "centrality_event_weight": False,
                "cross_section_weight_used_for_training": False,
                "weights_computed_before_binning": True,
            },
            "event_quality_filter": {"enabled": False},
            "the134_view_contract": {
                "schema": "THE134_FACTORIAL_VIEW_MODEL_VIEW_CONTRACT_V1",
                "system": "auau",
                "shower_definition": "H0",
                "shower_semantic_sha256": contract.shower_semantic_sha256("H0"),
            },
        }
        checks = model_validation.validate_metadata(metadata, "auau", expected_view="H0")
        self.assertTrue(all(checks.values()), checks)
        metadata["split"]["test_rows"] = 0
        bad = model_validation.validate_metadata(metadata, "auau", expected_view="H0")
        self.assertFalse(bad["split_population"])

    def test_auau_h70_control_view_cross_use_is_rejected(self):
        certificate = {
            "schema": "THE134_CONTROL_VIEW_AUTHORITY_V1",
            "status": "PASS",
            "system": "auau",
            "shower_definition": "H70",
            "shower_semantic_sha256": contract.shower_semantic_sha256("H70"),
            "feature_order": list(contract.AUAU_FEATURES),
            "model_origin": "REUSED_THE116",
            "model_tmva_sha256": "a" * 64,
            "model_metadata_sha256": "b" * 64,
        }
        checks = model_validation.control_view_certificate_checks(
            certificate,
            "auau",
            expected_control_view="H0",
            model_tmva_sha256="a" * 64,
            model_metadata_sha256="b" * 64,
        )
        self.assertFalse(checks["shower_definition"])
        self.assertFalse(checks["shower_semantic_sha256"])
        self.assertFalse(all(checks.values()))

    def test_external_training_weight_modes_are_rejected(self):
        metadata = {
            "features": list(contract.PP_FEATURES),
            "pt_range": list(contract.MODEL_DOMAIN_GEV),
            "weight_mode": "ppg12-exact",
            "split": {
                "mode": "event50",
                "test_fraction_requested": 0.5,
                "random_seed": 42,
                "train_events": 10,
                "test_events": 10,
            },
            "weighting": {
                "event_weight_used": True,
                "vertex_reweight": False,
                "centrality_event_weight": False,
                "cross_section_weight_used_for_training": False,
                "weights_computed_before_binning": True,
            },
            "event_quality_filter": {"enabled": False},
            "model_origin": "REUSED_THE116",
            "reuse_pinned_hashes": {
                key: "a" * 64
                for key in (
                    "source_sha256",
                    "trainer_sha256",
                    "pipeline_sha256",
                    "config_sha256",
                    "working_points_sha256",
                )
            },
            "the134_view_contract": {
                "shower_definition": "H70",
                "shower_semantic_sha256": contract.shower_semantic_sha256("H70"),
            },
        }
        checks = model_validation.validate_metadata(metadata, "pp", expected_view="H70")
        self.assertFalse(checks["event_weight_unused"])
        self.assertFalse(all(checks.values()))


class AuthorityChainTests(unittest.TestCase):
    def test_smoke_audit_cannot_reach_training_runner(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            matrix = root / "matrix.npz"
            np.savez_compressed(matrix, value=np.asarray([1]))
            audit = synthetic_full_extraction_audit(
                "pp", contract.sha256_file(matrix)
            )
            audit.update({"scope": "smoke", "full_training_authority": 0})
            audit_path = root / "audit.json"
            audit_path.write_text(json.dumps(audit))
            trainer = root / "trainer.py"
            trainer.write_text("raise SystemExit('must not execute')\n")
            args = SimpleNamespace(
                system="pp",
                view="H70",
                matrix=matrix,
                extraction_audit=audit_path,
                outdir=root / "out",
                trainer=trainer,
                python=sys.executable,
                reuse_audit=None,
                execute=False,
            )
            with mock.patch.object(runner, "parse_args", return_value=args):
                with self.assertRaisesRegex(SystemExit, "full extraction authority"):
                    runner.main()
            self.assertFalse((root / "out").exists())

    def test_full_label_with_truncated_source_inventory_is_rejected(self):
        system = "pp"
        paths = []
        path_sources = []
        provenance = {}
        authority = {}
        reports = []
        rows = {}
        for index, source in enumerate(contract.expected_sources(system)):
            path = Path(f"/frozen/{source}/view_{index}.root")
            record = {
                "path": str(path),
                "system": system,
                "source_sample": source,
                "input_uri_sha256": "1" * 64,
                "input_file_sha256": "2" * 64,
                "source_manifest_sha256": "3" * 64,
                "config_sha256": "4" * 64,
                "code_sha256": "5" * 64,
                "training_view_root_sha256": "6" * 64,
            }
            paths.append(path)
            path_sources.append(source)
            provenance[str(path)] = record
            reports.append(
                {
                    "source": source,
                    "status": "PASS",
                    "source_occurrence_ids": [f"{index + 1:032x}"],
                }
            )
            rows[source] = 1
            authority[source] = {
                "state": "SOURCE_COMPLETE",
                "full_source_manifest_sha256": "3" * 64,
                # The frozen inventory declares two inputs.  Presenting one
                # tuple per source cannot be promoted merely by saying full.
                "expected_input_count": 2,
                "expected_occurrence_count": 2,
                "input_records_sha256": contract.source_input_records_sha256(
                    [record]
                ),
            }
        closure = prepare.build_source_population_closure(
            system=system,
            paths=paths,
            path_sources=path_sources,
            external_provenance=provenance,
            source_coverage_authority=authority,
            file_reports=reports,
            selected_rows_by_source=rows,
        )
        self.assertEqual(closure["status"], "FAIL")
        self.assertTrue(
            all(
                not source["checks"]["input_count_exact"]
                and not source["checks"]["occurrence_count_exact"]
                for source in closure["sources"].values()
            )
        )

    def test_auau_h0_is_reuse_only(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            matrix = root / "matrix.npz"
            np.savez_compressed(matrix, value=np.asarray([1]))
            audit = synthetic_full_extraction_audit(
                "auau", contract.sha256_file(matrix), view_name="H0"
            )
            audit_path = root / "audit.json"
            audit_path.write_text(json.dumps(audit))
            trainer = root / "trainer.py"
            trainer.write_text("raise SystemExit('must not execute')\n")
            args = SimpleNamespace(
                system="auau",
                view="H0",
                matrix=matrix,
                extraction_audit=audit_path,
                outdir=root / "out",
                trainer=trainer,
                python=sys.executable,
                reuse_audit=None,
                execute=False,
            )
            with mock.patch.object(runner, "parse_args", return_value=args):
                with self.assertRaisesRegex(SystemExit, "must reuse"):
                    runner.main()

    def test_reuse_weight_handoff_requires_finite_positive_weights(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            valid = root / "valid.npz"
            np.savez_compressed(
                valid,
                __ppg12_exact_training_weight=np.asarray([0.5, 1.5]),
            )
            report = runner.validate_reuse_weighted_matrix(valid)
            self.assertEqual(report["rows"], 2)
            invalid = root / "invalid.npz"
            np.savez_compressed(
                invalid,
                __ppg12_exact_training_weight=np.asarray([0.5, np.nan]),
            )
            with self.assertRaises(SystemExit):
                runner.validate_reuse_weighted_matrix(invalid)

    def test_trained_model_receipt_binds_extraction_authority_and_xgboost(self):
        binding = {
            "schema": "THE134_FULL_EXTRACTION_AUTHORITY_BINDING_V1",
            "scope": "full",
            "full_training_authority": 1,
            "source_provenance_json_sha256": "a" * 64,
            "source_population_closure_sha256": "b" * 64,
        }
        receipt = runner.build_training_completion_receipt(
            system="auau",
            view_name="H70",
            model_id="centAsFeatBase3x3_pt15to35",
            extraction_binding=binding,
            artifacts={
                "xgboost_sha256": "c" * 64,
                "metadata_sha256": "d" * 64,
                "weighted_matrix_sha256": "e" * 64,
            },
            enrichment={"status": "PASS"},
        )
        good = contract.model_artifact_receipt_checks(
            receipt,
            "auau",
            "H70",
            model_xgb_sha256="c" * 64,
            model_metadata_sha256="d" * 64,
            weighted_matrix_sha256="e" * 64,
        )
        self.assertTrue(all(good.values()), good)
        bad = contract.model_artifact_receipt_checks(
            receipt,
            "auau",
            "H70",
            model_xgb_sha256="f" * 64,
            model_metadata_sha256="d" * 64,
        )
        self.assertFalse(bad["model_xgb_sha256"])

        missing = dict(receipt)
        missing.pop("extraction_authority")
        missing_checks = contract.model_artifact_receipt_checks(
            missing,
            "auau",
            "H70",
            model_xgb_sha256="c" * 64,
            model_metadata_sha256="d" * 64,
        )
        self.assertFalse(missing_checks["extraction_authority"])

        mutated = dict(receipt)
        mutated["extraction_authority"] = {
            **binding,
            "full_training_authority": 0,
        }
        mutated_checks = contract.model_artifact_receipt_checks(
            mutated,
            "auau",
            "H70",
            model_xgb_sha256="c" * 64,
            model_metadata_sha256="d" * 64,
        )
        self.assertFalse(mutated_checks["extraction_authority"])

    def test_registry_binds_working_points_to_validated_xgboost(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            extraction = synthetic_full_extraction_audit("auau", "9" * 64)
            extraction_path = root / "extraction.json"
            extraction_path.write_text(json.dumps(extraction))
            binding = contract.extraction_authority_binding(extraction)
            extraction_sha = contract.sha256_file(extraction_path)
            artifact_paths = {}
            for name in (
                "model_xgb",
                "model_tmva",
                "model_metadata",
                "model_receipt",
                "holdout",
                "holdout_certificate",
                "score_sample",
                "population_certificate",
            ):
                path = root / name
                path.write_bytes(f"immutable-{name}".encode("utf-8"))
                artifact_paths[name] = path

            def artifact_sha(name: str) -> str:
                return contract.sha256_file(artifact_paths[name])

            validation = {
                "schema": "THE134_FACTORIAL_VIEW_MODEL_VALIDATION_V1",
                "status": "PASS",
                "system": "auau",
                "shower_definition": "H70",
                "shower_semantic_sha256": contract.shower_semantic_sha256(),
                "feature_order": list(contract.AUAU_FEATURES),
                "model_domain_gev": list(contract.MODEL_DOMAIN_GEV),
                "model_origin": "TRAINED_THE134",
                "reuse_pinned_hashes": None,
                "critical_gates": {"all": True},
                "extraction_authority": binding,
                "provenance": {
                    "model_xgb": str(artifact_paths["model_xgb"]),
                    "model_xgb_sha256": artifact_sha("model_xgb"),
                    "model_tmva": str(artifact_paths["model_tmva"]),
                    "model_tmva_sha256": artifact_sha("model_tmva"),
                    "model_metadata": str(artifact_paths["model_metadata"]),
                    "model_metadata_sha256": artifact_sha("model_metadata"),
                    "model_receipt": str(artifact_paths["model_receipt"]),
                    "model_receipt_sha256": artifact_sha("model_receipt"),
                    "holdout": str(artifact_paths["holdout"]),
                    "holdout_sha256": artifact_sha("holdout"),
                    "holdout_certificate": str(artifact_paths["holdout_certificate"]),
                    "holdout_certificate_sha256": artifact_sha("holdout_certificate"),
                    "extraction_audit_sha256": extraction_sha,
                },
            }
            wp = {
                "schema": "THE134_FACTORIAL_VIEW_WEIGHTED_WORKING_POINTS_V1",
                "status": "PASS",
                "system": "auau",
                "shower_definition": "H70",
                "shower_semantic_sha256": contract.shower_semantic_sha256(),
                "feature_order": list(contract.AUAU_FEATURES),
                "model_domain_gev": list(contract.MODEL_DOMAIN_GEV),
                "model_origin": "TRAINED_THE134",
                "reuse_pinned_hashes": None,
                "gates": {"all": True},
                "extraction_authority": binding,
                "provenance": {
                    "model_xgb_sha256": artifact_sha("model_xgb"),
                    "model_metadata_sha256": artifact_sha("model_metadata"),
                    "model_receipt_sha256": artifact_sha("model_receipt"),
                    "score_sample": str(artifact_paths["score_sample"]),
                    "score_sample_sha256": artifact_sha("score_sample"),
                    "population_certificate": str(
                        artifact_paths["population_certificate"]
                    ),
                    "population_certificate_sha256": artifact_sha(
                        "population_certificate"
                    ),
                },
            }
            validation_path = root / "validation.json"
            wp_path = root / "wp.json"
            validation_path.write_text(json.dumps(validation))
            wp_path.write_text(json.dumps(wp))
            lane = registry.validate_lane(
                "auau",
                "H70",
                "TRAINED_THE134",
                extraction_path,
                validation_path,
                wp_path,
            )
            self.assertEqual(lane["status"], "READY")
            artifact_paths["model_xgb"].write_bytes(b"swapped-after-validation")
            with self.assertRaisesRegex(SystemExit, "changed after certification"):
                registry.validate_lane(
                    "auau",
                    "H70",
                    "TRAINED_THE134",
                    extraction_path,
                    validation_path,
                    wp_path,
                )
            artifact_paths["model_xgb"].write_bytes(b"immutable-model_xgb")
            wp["provenance"]["model_xgb_sha256"] = "e" * 64
            wp_path.write_text(json.dumps(wp))
            with self.assertRaisesRegex(SystemExit, "XGBoost"):
                registry.validate_lane(
                    "auau",
                    "H70",
                    "TRAINED_THE134",
                    extraction_path,
                    validation_path,
                    wp_path,
                )


class WorkingPointTests(unittest.TestCase):
    def test_wp_gate_widening_is_rejected(self):
        args = SimpleNamespace(
            surface_max_rms=contract.MAX_SURFACE_RMS,
            surface_max_residual=contract.MAX_SURFACE_ABS_RESIDUAL,
            surface_max_efficiency_error=contract.MAX_SURFACE_EFFICIENCY_ERROR,
            exact_max_efficiency_error=contract.MAX_WP_EXACT_EFFICIENCY_ERROR + 1.0e-6,
        )
        with self.assertRaises(SystemExit):
            working_points.enforce_frozen_gate_maxima(args)

    def test_weighted_threshold_and_strict_efficiency(self):
        scores = np.asarray([0.1, 0.2, 0.8, 0.9])
        weights = np.ones(4)
        threshold = working_points.weighted_threshold(scores, weights, 0.5)
        self.assertEqual(threshold, 0.2)
        self.assertEqual(working_points.weighted_efficiency(scores > threshold, weights), 0.5)

    def test_surface_rejects_large_residual(self):
        class Args:
            surface_max_rms = 1.0e-6
            surface_max_residual = 1.0e-6
            surface_max_efficiency_error = 1.0e-6

        edges = np.asarray([0.0, 1.0, 2.0])
        labels = np.asarray([1, 0, 1, 0])
        scores = np.asarray([0.9, 0.1, 0.8, 0.2])
        weights = np.ones(4)
        axis = np.asarray([0.5, 0.5, 1.5, 1.5])
        rows = working_points.derive_rows(scores, labels, weights, axis, edges)
        surfaces = working_points.fit_surfaces(rows, scores, labels, weights, axis, edges, Args())
        self.assertTrue(all(item["status"] == "BINNED_THRESHOLDS_ONLY" for item in surfaces.values()))

    def test_exact_wp_tie_population_fails_efficiency_gate(self):
        edges = np.asarray([0.0, 1.0])
        labels = np.asarray([1, 1, 1, 0])
        scores = np.asarray([0.8, 0.8, 0.8, 0.1])
        weights = np.ones(4)
        axis = np.full(4, 0.5)
        rows = working_points.derive_rows(scores, labels, weights, axis, edges)
        self.assertTrue(all(row["signal_weight_fraction_tied_at_threshold"] == 1.0 for row in rows))
        self.assertTrue(all(row["abs_efficiency_error"] > 0.02 for row in rows))
        self.assertFalse(working_points.exact_efficiency_gate(rows, 0.02))


class SplitAndCapTests(unittest.TestCase):
    @staticmethod
    def _write_shared_factorial_matrix(path: Path) -> dict[str, np.ndarray]:
        rows = 3
        payload: dict[str, np.ndarray] = {
            "is_signal": np.asarray([1, 0, 1], dtype=np.int8),
            "training_label": np.asarray([1, 0, 1], dtype=np.int8),
            "candidate_id_hi": np.asarray([10, 11, 12], dtype=np.uint64),
            "candidate_id_lo": np.asarray([20, 21, 22], dtype=np.uint64),
            "input_file_index": np.asarray([0, 0, 1], dtype=np.int64),
            "input_tree_entry": np.asarray([0, 7, 14], dtype=np.int64),
            "finite_feature_state": np.ones(rows, dtype=np.int32),
            "model_domain_state": np.zeros(rows, dtype=np.int32),
            "definition_name": np.full(rows, "H70", dtype="U4"),
            "shower_semantic_sha256": np.full(
                rows, contract.shower_semantic_sha256("H70"), dtype="U64"
            ),
        }
        for feature_index, feature in enumerate(contract.PP_FEATURES):
            h70 = np.asarray(
                [feature_index + 0.1, feature_index + 0.2, feature_index + 0.3],
                dtype=np.float32,
            )
            payload[feature] = h70.copy()
            for view_index, view_name in enumerate(contract.ALL_SHOWER_VIEWS):
                payload[f"view_{view_name}__{feature}"] = h70 + np.float32(view_index)
        for view_index, view_name in enumerate(contract.ALL_SHOWER_VIEWS):
            payload[f"view_{view_name}__finite_feature_state"] = np.ones(
                rows, dtype=np.int32
            )
            payload[f"view_{view_name}__model_domain_state"] = np.zeros(
                rows, dtype=np.int32
            )
            payload[f"view_{view_name}__input_tree_entry"] = (
                np.asarray([0, 7, 14], dtype=np.int64) + view_index
            )
        payload["__columns__"] = np.asarray(sorted(payload), dtype=object)
        np.savez_compressed(path, **payload)
        return payload

    def test_view_projection_is_exact_and_preserves_row_authority(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "shared.npz"
            expected = self._write_shared_factorial_matrix(source)
            output = root / "h0.npz"
            report = runner.project_factorial_view(
                source, output, system="pp", view_name="H0"
            )
            projected = np.load(output, allow_pickle=True)
            self.assertEqual(report["status"], "PASS")
            for feature in contract.PP_FEATURES:
                np.testing.assert_array_equal(
                    projected[feature], expected[f"view_H0__{feature}"]
                )
            np.testing.assert_array_equal(
                projected["candidate_id_hi"], expected["candidate_id_hi"]
            )
            np.testing.assert_array_equal(
                projected["training_label"], expected["training_label"]
            )
            np.testing.assert_array_equal(
                projected["input_tree_entry"],
                expected["view_H0__input_tree_entry"],
            )
            self.assertEqual(set(projected["definition_name"].tolist()), {"H0"})

    def test_view_projection_rejects_missing_column_and_h70_anchor_drift(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "shared.npz"
            payload = self._write_shared_factorial_matrix(source)
            missing = dict(payload)
            missing.pop("view_H0__cluster_Et")
            missing.pop("__columns__")
            missing["__columns__"] = np.asarray(sorted(missing), dtype=object)
            missing_path = root / "missing.npz"
            np.savez_compressed(missing_path, **missing)
            with self.assertRaisesRegex(SystemExit, "missing projection column"):
                runner.project_factorial_view(
                    missing_path, root / "missing-out.npz", system="pp", view_name="H0"
                )

            drifted = dict(payload)
            drifted["cluster_Et"] = np.asarray(drifted["cluster_Et"]).copy()
            drifted["cluster_Et"][0] += np.float32(1.0)
            drifted_path = root / "drifted.npz"
            np.savez_compressed(drifted_path, **drifted)
            with self.assertRaisesRegex(SystemExit, "generic/H70 feature drift"):
                runner.project_factorial_view(
                    drifted_path, root / "drifted-out.npz", system="pp", view_name="H0"
                )

    def test_runtime_gate_widening_is_rejected(self):
        args = SimpleNamespace(
            runtime_tolerance=contract.MAX_RUNTIME_SCORE_ABS_DIFFERENCE * 2,
            max_auc_regression=contract.MAX_OVERALL_AUC_REGRESSION,
            max_binned_auc_regression=contract.MAX_BINNED_AUC_REGRESSION,
            max_train_holdout_auc_gap=contract.MAX_TRAIN_HOLDOUT_AUC_GAP,
        )
        with self.assertRaises(SystemExit):
            model_validation.enforce_frozen_gate_maxima(args)

    def test_holdout_reads_typed_control_view_columns(self):
        features = contract.PP_FEATURES
        audited = {
            feature: np.asarray([index], dtype=np.float32)
            for index, feature in enumerate(features)
        }
        audited.update(
            {
                f"view_H70__{feature}": np.asarray([100 + index], dtype=np.float32)
                for index, feature in enumerate(features)
            }
        )
        selected = holdout.aligned_view_matrix(
            audited, np.asarray([0]), features, view_name=None
        )
        control = holdout.aligned_view_matrix(
            audited, np.asarray([0]), features, view_name="H70"
        )
        self.assertEqual(selected[0, 0], 0.0)
        self.assertEqual(control[0, 0], 100.0)
    def test_event50_is_deterministic_and_disjoint(self):
        keys = np.asarray([f"event-{index // 2}" for index in range(40)])
        first, report1 = holdout.event50_mask(keys, seed=42, model_id="model", test_fraction=0.5)
        second, report2 = holdout.event50_mask(keys, seed=42, model_id="model", test_fraction=0.5)
        np.testing.assert_array_equal(first, second)
        self.assertEqual(report1, report2)
        self.assertEqual(report1["event_overlap"], 0)

    def test_pp_cap_is_deterministic_and_bounded(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.npz"
            columns = ["is_signal", "input_file_index", "value"]
            np.savez_compressed(
                source,
                __columns__=np.asarray(columns, dtype=object),
                is_signal=np.asarray([0, 0, 1, 1] * 5),
                input_file_index=np.repeat(np.arange(5), 4),
                value=np.arange(20),
            )
            out1 = root / "one.npz"
            out2 = root / "two.npz"
            report1 = runner.deterministic_pp_cap(source, out1, 4, 42)
            report2 = runner.deterministic_pp_cap(source, out2, 4, 42)
            self.assertEqual(report1["class_counts"], {"0": 4, "1": 4})
            self.assertEqual(report1["output_sha256"], report2["output_sha256"])
            one = np.load(out1, allow_pickle=True)
            self.assertEqual(len(one["is_signal"]), 8)
            self.assertEqual(one["input_file_index"].tolist(), [4, 4, 4, 4, 2, 2, 2, 2])

    def test_enriched_weighted_cache_restores_audited_provenance(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            weighted = root / "weighted.npz"
            snapshot = root / "snapshot.npz"
            output = root / "enriched.npz"
            keys = {
                "source_sample": np.asarray(["sample", "sample"]),
                "input_file_index": np.asarray([0, 0]),
                "input_tree_entry": np.asarray([1, 2]),
            }
            weighted_columns = [*keys, "__ppg12_exact_training_weight"]
            np.savez_compressed(
                weighted,
                __columns__=np.asarray(weighted_columns, dtype=object),
                **keys,
                __ppg12_exact_training_weight=np.asarray([0.5, 1.5]),
            )
            snapshot_columns = [*keys, "label_authority", "view_H0__cluster_Et"]
            np.savez_compressed(
                snapshot,
                __columns__=np.asarray(snapshot_columns, dtype=object),
                **keys,
                label_authority=np.asarray(["PPG12_SOURCE_ROLE"] * 2),
                view_H0__cluster_Et=np.asarray([20.0, 21.0]),
            )
            report = runner.enrich_weighted_cache(weighted, snapshot, output)
            self.assertIn("label_authority", report["reattached_columns"])
            enriched = np.load(output, allow_pickle=True)
            self.assertEqual(enriched["label_authority"].tolist(), ["PPG12_SOURCE_ROLE"] * 2)


if __name__ == "__main__":
    unittest.main()
