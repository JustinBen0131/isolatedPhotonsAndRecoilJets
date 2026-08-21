#!/usr/bin/env python3
"""Regression tests for THE-134 aggregate science and storage receipts."""

from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
VALIDATION = ROOT / "scripts" / "ml" / "validation"
CONTRACTS = ROOT / "scripts" / "ml" / "contracts"
sys.path.insert(0, str(VALIDATION))
sys.path.insert(0, str(CONTRACTS))


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


aggregate_contract = load_module(
    "the134_aggregate_contract",
    VALIDATION / "the134_aggregate_contract.py",
)
h70_contract = load_module(
    "the134_h70_contract",
    CONTRACTS / "the134_h70_contract.py",
)

BUILD = VALIDATION / "build_the134_science_freeze_certificate.py"
PROJECT = VALIDATION / "project_the134_storage_quota.py"
PUBLIC_COMMIT = "1" * 40
CODE_SHA256 = "2" * 64


class AggregateFixture:
    def __init__(self, root: Path):
        self.root = root
        self.registries = []
        self.certificates = []
        self.artifact_hashes = {}
        for system, source in aggregate_contract.REQUIRED_SOURCE_KEYS:
            self.artifact_hashes[(system, source)] = {
                "source_manifest_sha256": h70_contract.canonical_json_sha256(
                    ["manifest", system, source]
                ),
                "direct_artifact_sha256": h70_contract.canonical_json_sha256(
                    ["direct", system, source]
                ),
                "writer_artifact_sha256": h70_contract.canonical_json_sha256(
                    ["writer", system, source]
                ),
            }
        for view in aggregate_contract.REQUIRED_VIEW_KEYS:
            registry_path = root / f"registry_{view}.json"
            registry = {
                "schema": aggregate_contract.PAIRED_REGISTRY_SCHEMA,
                "status": "READY",
                "promotion_status": "CANDIDATE_CURRENT_NOT_CANONICAL",
                "model_count": 2,
                "systems": ["pp", "auau"],
                "shower_definition": view,
                "shower_semantic_sha256": h70_contract.shower_semantic_sha256(
                    view
                ),
                "public_commit": PUBLIC_COMMIT,
                "code_sha256": CODE_SHA256,
                "models": [
                    self.make_registry_model(system, view)
                    for system in aggregate_contract.SYSTEMS
                ],
                "boundaries": [
                    "registry readiness does not promote either model to CANONICAL",
                    "view-specific direct/writer/cache replay certification remains required",
                    "THE-121/THE-122 remain closed until final P5C certification",
                ],
            }
            registry["registry_semantic_sha256"] = (
                h70_contract.canonical_json_sha256(registry)
            )
            self.write_json(registry_path, registry)
            registry_sha = h70_contract.sha256_file(registry_path)
            self.registries.append(
                {"view": view, "path": str(registry_path), "sha256": registry_sha}
            )
            for system in aggregate_contract.SYSTEMS:
                certificate_path = root / f"replay_{system}_{view}.json"
                witnesses = []
                for source in aggregate_contract.REQUIRED_SOURCES_BY_SYSTEM[system]:
                    stable = self.artifact_hashes[(system, source)]
                    witnesses.append(
                        {
                            "system": system,
                            "shower_definition": view,
                            "source": source,
                            **stable,
                            "cache_receipt_sha256": (
                                h70_contract.canonical_json_sha256(
                                    ["cache", system, source, view]
                                )
                            ),
                            "gates": {
                                gate: True
                                for gate in aggregate_contract.SOURCE_WITNESS_GATES
                            },
                        }
                    )
                certificate = {
                    "schema": aggregate_contract.REPLAY_CERTIFICATE_SCHEMA,
                    "status": "PASS",
                    "promotion_status": "CANDIDATE_CURRENT_NOT_CANONICAL",
                    "system": system,
                    "shower_definition": view,
                    "shower_semantic_sha256": (
                        h70_contract.shower_semantic_sha256(view)
                    ),
                    "paired_registry_sha256": registry_sha,
                    "public_commit": PUBLIC_COMMIT,
                    "code_sha256": CODE_SHA256,
                    "gates": {
                        gate: True
                        for gate in aggregate_contract.REPLAY_CERTIFICATE_GATES
                    },
                    "source_witnesses": witnesses,
                }
                self.write_json(certificate_path, certificate)
                self.certificates.append(
                    {
                        "system": system,
                        "view": view,
                        "path": str(certificate_path),
                        "sha256": h70_contract.sha256_file(certificate_path),
                    }
                )
        self.manifest_path = root / "aggregate_manifest.json"
        self.write_manifest()

    def write_blob(self, name: str, content: str | None = None) -> Path:
        path = self.root / name
        path.write_text(content or f"fixture artifact {name}\n")
        return path

    def threshold_rows(self, system: str) -> list[dict]:
        edges = aggregate_contract.WORKING_POINT_EDGES_BY_SYSTEM[system]
        rows = []
        for target in h70_contract.TARGET_EFFICIENCIES:
            for index, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
                rows.append(
                    {
                        "wp": aggregate_contract.WORKING_POINT_LABEL_BY_TARGET[
                            target
                        ],
                        "target_signal_efficiency": target,
                        "bin_lo": lo,
                        "bin_hi": hi,
                        "bin_center": 0.5 * (lo + hi),
                        "threshold": 0.1 + target + 0.001 * index,
                        "achieved_signal_efficiency": target,
                        "achieved_minus_target": 0.0,
                        "abs_efficiency_error": 0.0,
                        "signal_weight_fraction_tied_at_threshold": 0.0,
                        "background_acceptance": 0.2,
                        "signal_rows": 10,
                        "background_rows": 20,
                    }
                )
        return rows

    def runtime_surfaces(self, system: str) -> dict:
        bin_count = (
            len(aggregate_contract.WORKING_POINT_EDGES_BY_SYSTEM[system]) - 1
        )
        return {
            aggregate_contract.WORKING_POINT_LABEL_BY_TARGET[target]: {
                "status": "ACCEPTED",
                "intercept": 0.1 + target,
                "slope_per_axis_unit": 0.001,
                "rms_residual": 0.001,
                "max_abs_residual": 0.002,
                "max_abs_efficiency_error": 0.001,
                "achieved_efficiency_by_bin": [target] * bin_count,
                "gates": {
                    "max_rms": h70_contract.MAX_SURFACE_RMS,
                    "max_abs_residual": (
                        h70_contract.MAX_SURFACE_ABS_RESIDUAL
                    ),
                    "max_abs_efficiency_error": (
                        h70_contract.MAX_SURFACE_EFFICIENCY_ERROR
                    ),
                },
            }
            for target in h70_contract.TARGET_EFFICIENCIES
        }

    def make_registry_model(self, system: str, view: str) -> dict:
        prefix = f"{system}_{view}"
        artifacts = {
            name: self.write_blob(f"{prefix}_{name}{suffix}")
            for name, suffix in (
                ("xgboost", ".json"),
                ("tmva", ".xml"),
                ("metadata", ".json"),
                ("extraction", ".json"),
                ("model_receipt", ".json"),
                ("holdout", ".npz"),
                ("holdout_certificate", ".json"),
                ("working_point_score_sample", ".npz"),
                ("working_point_population_certificate", ".json"),
            )
        }
        model_origin = h70_contract.expected_model_origin(system, view)
        reuse_pins = (
            dict(h70_contract.FROZEN_REUSE_PINNED_HASHES[(system, view)])
            if (system, view) in h70_contract.FROZEN_REUSE_PINNED_HASHES
            else None
        )
        model = {
            "xgboost": str(artifacts["xgboost"]),
            "xgboost_sha256": h70_contract.sha256_file(artifacts["xgboost"]),
            "tmva": str(artifacts["tmva"]),
            "tmva_sha256": h70_contract.sha256_file(artifacts["tmva"]),
            "metadata": str(artifacts["metadata"]),
            "metadata_sha256": h70_contract.sha256_file(artifacts["metadata"]),
        }
        required_sources = sorted(h70_contract.expected_sources(system))
        source_rows = {
            source: {
                "status": "PASS",
                "checks": {"fixture_exact_closure": True},
                "observed_input_count": 1,
                "observed_occurrence_count": 1,
                "full_source_manifest_sha256": (
                    h70_contract.canonical_json_sha256(
                        ["full source manifest", system, source]
                    )
                ),
                "input_records_sha256": h70_contract.canonical_json_sha256(
                    ["input records", system, source]
                ),
            }
            for source in required_sources
        }
        extraction = {
            "schema": "THE134_FACTORIAL_VIEW_TRAINING_MATRIX_AUDIT_V1",
            "status": "PASS",
            "system": system,
            "scope": "full",
            "full_training_authority": 1,
            "shower_definition": view,
            "shower_semantic_sha256": (
                h70_contract.shower_semantic_sha256(view)
            ),
            "feature_order": list(h70_contract.FEATURES_BY_SYSTEM[system]),
            "required_sources": required_sources,
            "observed_manifest_sources": required_sources,
            "source_provenance_json_sha256": (
                h70_contract.canonical_json_sha256(
                    ["source provenance", system, view]
                )
            ),
            "source_population_closure": {
                "status": "PASS",
                "sources": source_rows,
                "semantic_sha256": h70_contract.canonical_json_sha256(
                    source_rows
                ),
            },
        }
        self.write_json(artifacts["extraction"], extraction)
        extraction_binding = h70_contract.extraction_authority_binding(
            extraction
        )
        model_receipt = {
            "schema": (
                "THE134_FACTORIAL_VIEW_MODEL_REUSE_COMPLETE_V1"
                if model_origin.startswith("REUSED_")
                else "THE134_FACTORIAL_VIEW_MODEL_TRAINING_COMPLETE_V1"
            ),
            "status": "PASS",
            "system": system,
            "shower_definition": view,
            "shower_semantic_sha256": (
                h70_contract.shower_semantic_sha256(view)
            ),
            "model_origin": model_origin,
            "extraction_authority": extraction_binding,
            "artifacts": {
                "xgboost_sha256": model["xgboost_sha256"],
                "metadata_sha256": model["metadata_sha256"],
            },
        }
        self.write_json(artifacts["model_receipt"], model_receipt)
        validation_path = self.root / f"{prefix}_validation.json"
        validation = {
            "schema": aggregate_contract.MODEL_VALIDATION_SCHEMA,
            "status": "PASS",
            "promotion_status": "CANDIDATE_CURRENT_NOT_CANONICAL",
            "system": system,
            "shower_definition": view,
            "shower_semantic_sha256": (
                h70_contract.shower_semantic_sha256(view)
            ),
            "feature_order": list(h70_contract.FEATURES_BY_SYSTEM[system]),
            "model_domain_gev": list(h70_contract.MODEL_DOMAIN_GEV),
            "model_origin": model_origin,
            "reuse_pinned_hashes": reuse_pins,
            "extraction_authority": extraction_binding,
            "critical_gates": {
                gate: True for gate in aggregate_contract.MODEL_VALIDATION_GATES
            },
            "metadata_gates": {"fixture_metadata_identity": True},
            "control_metadata_gates": {"fixture_control_identity": True},
            "control_view_gates": {"fixture_control_view_authority": True},
            "extraction_gates": {"fixture_full_extraction_authority": True},
            "holdout_certificate_gates": {
                "fixture_exact_event_group_holdout": True
            },
            "runtime_parity": {
                "rows": 100,
                "max_python_tmva_abs_difference": 1.0e-8,
                "max_cached_python_abs_difference": 1.0e-12,
                "tolerance": h70_contract.MAX_RUNTIME_SCORE_ABS_DIFFERENCE,
            },
            "source_label_closure": {
                "observed_sources": required_sources,
                "wrong_source_rows": 0,
                "wrong_label_rows": 0,
            },
            "performance": {"fixture_weighted_holdout_auc": 0.9},
            "provenance": {
                "model_xgb": model["xgboost"],
                "model_xgb_sha256": model["xgboost_sha256"],
                "model_tmva": model["tmva"],
                "model_tmva_sha256": model["tmva_sha256"],
                "model_metadata": model["metadata"],
                "model_metadata_sha256": model["metadata_sha256"],
                "model_receipt": str(artifacts["model_receipt"]),
                "model_receipt_sha256": h70_contract.sha256_file(
                    artifacts["model_receipt"]
                ),
                "holdout": str(artifacts["holdout"]),
                "holdout_sha256": h70_contract.sha256_file(artifacts["holdout"]),
                "holdout_certificate": str(artifacts["holdout_certificate"]),
                "holdout_certificate_sha256": h70_contract.sha256_file(
                    artifacts["holdout_certificate"]
                ),
                "extraction_audit": str(artifacts["extraction"]),
                "extraction_audit_sha256": h70_contract.sha256_file(
                    artifacts["extraction"]
                ),
            },
        }
        self.write_json(validation_path, validation)

        exact_thresholds = self.threshold_rows(system)
        runtime_surfaces = self.runtime_surfaces(system)
        working_points_path = self.root / f"{prefix}_working_points.json"
        working_points = {
            "schema": aggregate_contract.WORKING_POINTS_SCHEMA,
            "status": "PASS",
            "promotion_status": "CANDIDATE_CURRENT_NOT_CANONICAL",
            "system": system,
            "shower_definition": view,
            "shower_semantic_sha256": (
                h70_contract.shower_semantic_sha256(view)
            ),
            "feature_order": list(h70_contract.FEATURES_BY_SYSTEM[system]),
            "model_domain_gev": list(h70_contract.MODEL_DOMAIN_GEV),
            "model_origin": model_origin,
            "reuse_pinned_hashes": reuse_pins,
            "extraction_authority": extraction_binding,
            "axis": aggregate_contract.WORKING_POINT_AXIS_BY_SYSTEM[system],
            "bin_edges": list(
                aggregate_contract.WORKING_POINT_EDGES_BY_SYSTEM[system]
            ),
            "targets": list(h70_contract.TARGET_EFFICIENCIES),
            "threshold_authority": (
                aggregate_contract.WORKING_POINT_THRESHOLD_AUTHORITY
            ),
            "comparison_operator": (
                aggregate_contract.WORKING_POINT_COMPARISON_OPERATOR
            ),
            "tie_policy": aggregate_contract.WORKING_POINT_TIE_POLICY,
            "exact_max_efficiency_error": (
                h70_contract.MAX_WP_EXACT_EFFICIENCY_ERROR
            ),
            "gates": {
                gate: True for gate in aggregate_contract.WORKING_POINT_GATES
            },
            "population_certificate_gates": {
                "fixture_population_identity": True
            },
            "exact_thresholds": exact_thresholds,
            "runtime_surfaces": runtime_surfaces,
            "provenance": {
                "score_sample": str(artifacts["working_point_score_sample"]),
                "score_sample_sha256": h70_contract.sha256_file(
                    artifacts["working_point_score_sample"]
                ),
                "population_certificate": str(
                    artifacts["working_point_population_certificate"]
                ),
                "population_certificate_sha256": h70_contract.sha256_file(
                    artifacts["working_point_population_certificate"]
                ),
                "model_metadata": model["metadata"],
                "model_metadata_sha256": model["metadata_sha256"],
                "model_xgb": model["xgboost"],
                "model_xgb_sha256": model["xgboost_sha256"],
                "model_receipt": str(artifacts["model_receipt"]),
                "model_receipt_sha256": h70_contract.sha256_file(
                    artifacts["model_receipt"]
                ),
            },
        }
        self.write_json(working_points_path, working_points)

        certificates = {
            "extraction": str(artifacts["extraction"]),
            "extraction_sha256": h70_contract.sha256_file(
                artifacts["extraction"]
            ),
            "validation": str(validation_path),
            "validation_sha256": h70_contract.sha256_file(validation_path),
            "working_points": str(working_points_path),
            "working_points_sha256": h70_contract.sha256_file(
                working_points_path
            ),
            "model_receipt": str(artifacts["model_receipt"]),
            "model_receipt_sha256": h70_contract.sha256_file(
                artifacts["model_receipt"]
            ),
            "holdout": str(artifacts["holdout"]),
            "holdout_sha256": h70_contract.sha256_file(artifacts["holdout"]),
            "holdout_certificate": str(artifacts["holdout_certificate"]),
            "holdout_certificate_sha256": h70_contract.sha256_file(
                artifacts["holdout_certificate"]
            ),
            "working_point_score_sample": str(
                artifacts["working_point_score_sample"]
            ),
            "working_point_score_sample_sha256": h70_contract.sha256_file(
                artifacts["working_point_score_sample"]
            ),
            "working_point_population_certificate": str(
                artifacts["working_point_population_certificate"]
            ),
            "working_point_population_certificate_sha256": (
                h70_contract.sha256_file(
                    artifacts["working_point_population_certificate"]
                )
            ),
        }
        return {
            "system": system,
            "status": "READY",
            "model_origin": model_origin,
            "reuse_pinned_hashes": reuse_pins,
            "shower_definition": view,
            "shower_semantic_sha256": (
                h70_contract.shower_semantic_sha256(view)
            ),
            "feature_order": list(h70_contract.FEATURES_BY_SYSTEM[system]),
            "model_domain_gev": list(h70_contract.MODEL_DOMAIN_GEV),
            "model": model,
            "exact_working_points": exact_thresholds,
            "working_point_axis": (
                aggregate_contract.WORKING_POINT_AXIS_BY_SYSTEM[system]
            ),
            "working_point_bin_edges": list(
                aggregate_contract.WORKING_POINT_EDGES_BY_SYSTEM[system]
            ),
            "threshold_authority": (
                aggregate_contract.WORKING_POINT_THRESHOLD_AUTHORITY
            ),
            "runtime_surfaces": runtime_surfaces,
            "certificates": certificates,
        }

    @staticmethod
    def write_json(path: Path, payload: dict) -> None:
        path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    def read_registry(self, index: int = 0) -> dict:
        return json.loads(Path(self.registries[index]["path"]).read_text())

    def rewrite_registry(self, registry: dict, index: int = 0) -> None:
        registry.pop("registry_semantic_sha256", None)
        registry["registry_semantic_sha256"] = (
            h70_contract.canonical_json_sha256(registry)
        )
        path = Path(self.registries[index]["path"])
        self.write_json(path, registry)
        self.registries[index]["sha256"] = h70_contract.sha256_file(path)
        self.write_manifest()

    def write_manifest(self) -> None:
        self.write_json(
            self.manifest_path,
            {
                "schema": aggregate_contract.SCIENCE_MANIFEST_SCHEMA,
                "public_commit": PUBLIC_COMMIT,
                "code_sha256": CODE_SHA256,
                "paired_registries": self.registries,
                "replay_certificates": self.certificates,
            },
        )


class The134AggregateCertificationTest(unittest.TestCase):
    def run_tool(self, script: Path, manifest: Path, output: Path):
        return subprocess.run(
            [
                sys.executable,
                str(script),
                "--manifest",
                str(manifest),
                "--json-out",
                str(output),
            ],
            text=True,
            capture_output=True,
            check=False,
        )

    def test_science_certificate_is_complete_and_deterministic(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = AggregateFixture(root)
            output_a = root / "science_a.json"
            output_b = root / "science_b.json"
            result_a = self.run_tool(BUILD, fixture.manifest_path, output_a)
            result_b = self.run_tool(BUILD, fixture.manifest_path, output_b)
            self.assertEqual(result_a.returncode, 0, result_a.stderr)
            self.assertEqual(result_b.returncode, 0, result_b.stderr)
            self.assertEqual(output_a.read_bytes(), output_b.read_bytes())
            payload = json.loads(output_a.read_text())
            self.assertEqual(payload["status"], "PASS")
            self.assertEqual(payload["paired_registry_count"], 7)
            self.assertEqual(payload["replay_certificate_count"], 14)
            self.assertEqual(payload["training_source_count"], 13)
            self.assertEqual(
                payload["training_source_view_witness_count"], 91
            )
            self.assertNotIn("source_count", payload)
            self.assertNotIn(
                "source_complete_direct_writer_cache_witnesses",
                payload["gates"],
            )
            self.assertTrue(
                payload["gates"][
                    "training_source_complete_direct_writer_cache_witnesses"
                ]
            )
            self.assertTrue(all(payload["gates"].values()))

    def test_science_certificate_rejects_missing_view_and_source(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = AggregateFixture(root)
            fixture.registries.pop()
            fixture.write_manifest()
            result = self.run_tool(
                BUILD, fixture.manifest_path, root / "missing_view.json"
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("inventory mismatch", result.stderr)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = AggregateFixture(root)
            certificate_record = fixture.certificates[0]
            certificate_path = Path(certificate_record["path"])
            certificate = json.loads(certificate_path.read_text())
            certificate["source_witnesses"].pop()
            fixture.write_json(certificate_path, certificate)
            certificate_record["sha256"] = h70_contract.sha256_file(certificate_path)
            fixture.write_manifest()
            result = self.run_tool(
                BUILD, fixture.manifest_path, root / "missing_source.json"
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("source witnesses inventory mismatch", result.stderr)

    def test_science_certificate_rejects_unpaired_registry(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = AggregateFixture(root)
            registry = fixture.read_registry()
            registry["models"][1]["system"] = "pp"
            fixture.rewrite_registry(registry)
            result = self.run_tool(
                BUILD, fixture.manifest_path, root / "unpaired.json"
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("exactly one pp and one auau", result.stderr)

    def test_science_certificate_rejects_legacy_or_invalid_wp_rows(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = AggregateFixture(root)
            registry = fixture.read_registry()
            registry["models"][0]["exact_working_points"] = {
                "WP70": 0.7,
                "WP80": 0.8,
                "WP90": 0.9,
            }
            fixture.rewrite_registry(registry)
            result = self.run_tool(
                BUILD, fixture.manifest_path, root / "legacy_wp.json"
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("must be a list", result.stderr)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = AggregateFixture(root)
            registry = fixture.read_registry()
            model = registry["models"][0]
            working_points_path = Path(
                model["certificates"]["working_points"]
            )
            working_points = json.loads(working_points_path.read_text())
            working_points["exact_thresholds"][0]["threshold"] = float("nan")
            model["exact_working_points"][0]["threshold"] = float("nan")
            fixture.write_json(working_points_path, working_points)
            model["certificates"]["working_points_sha256"] = (
                h70_contract.sha256_file(working_points_path)
            )
            fixture.rewrite_registry(registry)
            result = self.run_tool(
                BUILD, fixture.manifest_path, root / "nonfinite_wp.json"
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("must be a finite number", result.stderr)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = AggregateFixture(root)
            registry = fixture.read_registry()
            model = registry["models"][0]
            working_points_path = Path(
                model["certificates"]["working_points"]
            )
            working_points = json.loads(working_points_path.read_text())
            working_points["exact_thresholds"].pop()
            model["exact_working_points"].pop()
            fixture.write_json(working_points_path, working_points)
            model["certificates"]["working_points_sha256"] = (
                h70_contract.sha256_file(working_points_path)
            )
            fixture.rewrite_registry(registry)
            result = self.run_tool(
                BUILD, fixture.manifest_path, root / "missing_wp_bin.json"
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("row count", result.stderr)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = AggregateFixture(root)
            registry = fixture.read_registry()
            model = registry["models"][0]
            working_points_path = Path(
                model["certificates"]["working_points"]
            )
            working_points = json.loads(working_points_path.read_text())
            working_points["bin_edges"][1] += 0.5
            fixture.write_json(working_points_path, working_points)
            model["certificates"]["working_points_sha256"] = (
                h70_contract.sha256_file(working_points_path)
            )
            fixture.rewrite_registry(registry)
            result = self.run_tool(
                BUILD, fixture.manifest_path, root / "wrong_wp_edges.json"
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("bin_edges", result.stderr)

    def test_science_certificate_rejects_fabricated_receipt_gates(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = AggregateFixture(root)
            registry = fixture.read_registry()
            model = registry["models"][0]
            validation_path = Path(model["certificates"]["validation"])
            validation = json.loads(validation_path.read_text())
            validation["critical_gates"] = {"trust_me": True}
            fixture.write_json(validation_path, validation)
            model["certificates"]["validation_sha256"] = (
                h70_contract.sha256_file(validation_path)
            )
            fixture.rewrite_registry(registry)
            result = self.run_tool(
                BUILD, fixture.manifest_path, root / "fabricated_gates.json"
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("gate inventory mismatch", result.stderr)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = AggregateFixture(root)
            registry = fixture.read_registry()
            model = registry["models"][0]
            validation_path = Path(model["certificates"]["validation"])
            fixture.write_json(
                validation_path,
                {
                    "schema": aggregate_contract.MODEL_VALIDATION_SCHEMA,
                    "status": "PASS",
                    "critical_gates": {
                        gate: True
                        for gate in aggregate_contract.MODEL_VALIDATION_GATES
                    },
                },
            )
            model["certificates"]["validation_sha256"] = (
                h70_contract.sha256_file(validation_path)
            )
            fixture.rewrite_registry(registry)
            result = self.run_tool(
                BUILD, fixture.manifest_path, root / "gate_only_receipt.json"
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("key inventory mismatch", result.stderr)

    def test_science_certificate_rejects_wp_policy_and_model_hash_drift(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = AggregateFixture(root)
            registry = fixture.read_registry()
            model = registry["models"][0]
            working_points_path = Path(
                model["certificates"]["working_points"]
            )
            working_points = json.loads(working_points_path.read_text())
            working_points["comparison_operator"] = "score >= threshold"
            fixture.write_json(working_points_path, working_points)
            model["certificates"]["working_points_sha256"] = (
                h70_contract.sha256_file(working_points_path)
            )
            fixture.rewrite_registry(registry)
            result = self.run_tool(
                BUILD, fixture.manifest_path, root / "comparison_drift.json"
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("comparison_operator", result.stderr)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = AggregateFixture(root)
            registry = fixture.read_registry()
            model_path = Path(registry["models"][0]["model"]["xgboost"])
            model_path.write_text("mutated model\n")
            result = self.run_tool(
                BUILD, fixture.manifest_path, root / "model_drift.json"
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("model xgboost SHA-256 mismatch", result.stderr)

    def test_science_certificate_rejects_feature_and_runtime_identity_drift(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = AggregateFixture(root)
            registry = fixture.read_registry()
            model = registry["models"][0]
            validation_path = Path(model["certificates"]["validation"])
            validation = json.loads(validation_path.read_text())
            validation["feature_order"] = list(reversed(validation["feature_order"]))
            fixture.write_json(validation_path, validation)
            model["certificates"]["validation_sha256"] = (
                h70_contract.sha256_file(validation_path)
            )
            fixture.rewrite_registry(registry)
            result = self.run_tool(
                BUILD, fixture.manifest_path, root / "feature_drift.json"
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("feature_order mismatch", result.stderr)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = AggregateFixture(root)
            registry = fixture.read_registry()
            model = registry["models"][0]
            validation_path = Path(model["certificates"]["validation"])
            validation = json.loads(validation_path.read_text())
            validation["runtime_parity"]["tolerance"] = 4.0e-7
            fixture.write_json(validation_path, validation)
            model["certificates"]["validation_sha256"] = (
                h70_contract.sha256_file(validation_path)
            )
            fixture.rewrite_registry(registry)
            result = self.run_tool(
                BUILD, fixture.manifest_path, root / "runtime_drift.json"
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("runtime_parity.tolerance", result.stderr)

    def make_storage_manifest(
        self,
        path: Path,
        science_certificate: Path,
        *,
        quota_bytes: int = 1_000_000_000,
        used_bytes: int = 100_000_000,
        merge_count: int = 2,
    ) -> dict:
        projections = []
        for system, source in aggregate_contract.REQUIRED_SOURCE_KEYS:
            projections.append(
                {
                    "system": system,
                    "source": source,
                    "measurement_receipt_sha256": (
                        h70_contract.canonical_json_sha256(
                            ["measurement", system, source]
                        )
                    ),
                    "projected_direct_reference_bytes": 1_000_000,
                    "projected_replay_ttree_bytes": 2_000_000,
                    "projected_final_merged_bytes": 500_000,
                    "projected_evidence_bytes": 100_000,
                    "merge_workspace_bytes": (
                        4_000_000
                        if (system, source) == ("auau", "run28_embeddedJet40")
                        else 3_000_000
                    ),
                }
            )
        payload = {
            "schema": aggregate_contract.STORAGE_MANIFEST_SCHEMA,
            "science_freeze_certificate": {
                "path": str(science_certificate),
                "sha256": h70_contract.sha256_file(science_certificate),
            },
            "quota": {
                "quota_bytes": quota_bytes,
                "used_bytes": used_bytes,
                "usage_authority_sha256": "3" * 64,
                "minimum_free_fraction": 0.20,
            },
            "simultaneous_merge_workspace_count": merge_count,
            "source_projections": projections,
            "local_hot_tier": {
                "free_bytes": 100_000_000,
                "protected_free_bytes": 50_000_000,
                "planned_pull_bytes": 10_000_000,
            },
        }
        AggregateFixture.write_json(path, payload)
        return payload

    def test_storage_projection_reserves_two_workspaces_and_headroom(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = AggregateFixture(root)
            science = root / "science.json"
            self.assertEqual(
                self.run_tool(BUILD, fixture.manifest_path, science).returncode, 0
            )
            manifest = root / "storage_manifest.json"
            self.make_storage_manifest(manifest, science)
            output_a = root / "storage_a.json"
            output_b = root / "storage_b.json"
            result_a = self.run_tool(PROJECT, manifest, output_a)
            result_b = self.run_tool(PROJECT, manifest, output_b)
            self.assertEqual(result_a.returncode, 0, result_a.stderr)
            self.assertEqual(result_b.returncode, 0, result_b.stderr)
            self.assertEqual(output_a.read_bytes(), output_b.read_bytes())
            payload = json.loads(output_a.read_text())
            self.assertEqual(payload["status"], "PASS")
            self.assertEqual(
                payload["totals"]["simultaneous_merge_workspace_count"], 2
            )
            self.assertEqual(
                payload["totals"]["merge_workspace_reserve_bytes"], 8_000_000
            )
            self.assertGreaterEqual(
                payload["totals"]["projected_peak_free_fraction"], 0.20
            )

    def test_storage_projection_fails_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = AggregateFixture(root)
            science = root / "science.json"
            self.assertEqual(
                self.run_tool(BUILD, fixture.manifest_path, science).returncode, 0
            )
            manifest = root / "storage_manifest.json"
            payload = self.make_storage_manifest(
                manifest,
                science,
                quota_bytes=100_000_000,
                used_bytes=90_000_000,
            )
            result = self.run_tool(PROJECT, manifest, root / "headroom.json")
            self.assertEqual(result.returncode, 3)
            receipt = json.loads((root / "headroom.json").read_text())
            self.assertFalse(
                receipt["gates"]["minimum_free_headroom_at_least_20_percent"]
            )

            payload["source_projections"].pop()
            AggregateFixture.write_json(manifest, payload)
            result = self.run_tool(PROJECT, manifest, root / "missing_source.json")
            self.assertEqual(result.returncode, 2)
            self.assertIn("source projections inventory mismatch", result.stderr)

            self.make_storage_manifest(manifest, science, merge_count=1)
            result = self.run_tool(PROJECT, manifest, root / "merge_count.json")
            self.assertEqual(result.returncode, 2)
            self.assertIn("must equal exactly 2", result.stderr)

    def test_storage_projection_rejects_tampered_science_certificate(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = AggregateFixture(root)
            science = root / "science.json"
            self.assertEqual(
                self.run_tool(BUILD, fixture.manifest_path, science).returncode, 0
            )
            science_payload = json.loads(science.read_text())
            science_payload["boundaries"].append("unregistered mutation")
            AggregateFixture.write_json(science, science_payload)
            manifest = root / "storage_manifest.json"
            self.make_storage_manifest(manifest, science)
            result = self.run_tool(
                PROJECT, manifest, root / "semantic_tamper.json"
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("aggregate_semantic_sha256 mismatch", result.stderr)

            science_payload["boundaries"].pop()
            science_payload["gates"][
                "source_complete_direct_writer_cache_witnesses"
            ] = science_payload["gates"].pop(
                "training_source_complete_direct_writer_cache_witnesses"
            )
            science_payload.pop("aggregate_semantic_sha256")
            science_payload["aggregate_semantic_sha256"] = (
                h70_contract.canonical_json_sha256(science_payload)
            )
            AggregateFixture.write_json(science, science_payload)
            self.make_storage_manifest(manifest, science)
            result = self.run_tool(
                PROJECT, manifest, root / "retired_gate.json"
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("gate inventory mismatch", result.stderr)


if __name__ == "__main__":
    unittest.main()
