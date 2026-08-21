#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


HERE = Path(__file__).resolve()
SCRIPT = HERE.parents[1] / "build_the134_factorial_execution_plan.py"
SPEC = importlib.util.spec_from_file_location("the134_execution_plan", SCRIPT)
assert SPEC and SPEC.loader
planner = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(planner)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Fixture:
    def __init__(self, root: Path):
        self.root = root
        self.manifests: dict[str, Path] = {}
        self.reuse: dict[tuple[str, str], Path] = {}
        for system in planner.SYSTEMS:
            matrix = root / f"{system}.npz"
            matrix.write_bytes(f"matrix-{system}".encode())
            audits: dict[str, dict[str, str]] = {}
            for view in planner.ALL_SHOWER_VIEWS:
                audit_path = root / f"{system}-{view}.audit.json"
                audit = {
                    "schema": planner.MATRIX_AUDIT_SCHEMA,
                    "status": "PASS",
                    "system": system,
                    "shower_definition": view,
                    "shower_semantic_sha256": planner.shower_semantic_sha256(view),
                    "matrix": str(matrix),
                    "matrix_sha256": sha256(matrix),
                    "full_training_authority": 1,
                    "source_population_closure": {"status": "PASS"},
                    "single_read_factorial_projection": {"status": "PASS"},
                }
                audit_path.write_text(json.dumps(audit), encoding="utf-8")
                audits[view] = {"path": str(audit_path), "sha256": sha256(audit_path)}
            manifest_path = root / f"{system}.manifest.json"
            manifest = {
                "schema": planner.MATRIX_MANIFEST_SCHEMA,
                "status": "PASS",
                "system": system,
                "matrix": str(matrix),
                "matrix_sha256": sha256(matrix),
                "root_reads": 3,
                "input_count": 3,
                "audits": audits,
            }
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            self.manifests[system] = manifest_path

        for system, view in planner.REQUIRED_REUSE:
            path = root / f"reuse-{system}-{view}.json"
            path.write_text(
                json.dumps(
                    {
                        "schema": planner.REUSE_AUDIT_SCHEMA,
                        "status": "PASS",
                        "system": system,
                        "shower_definition": view,
                        "shower_semantic_sha256": planner.shower_semantic_sha256(view),
                        "model_origin": planner.expected_model_origin(system, view),
                    }
                ),
                encoding="utf-8",
            )
            self.reuse[(system, view)] = path

    def args(self, **overrides):
        values = {
            "pp_matrix_manifest": self.manifests["pp"],
            "auau_matrix_manifest": self.manifests["auau"],
            "reuse_audit": [
                f"pp:H70={self.reuse[('pp', 'H70')]}",
                f"auau:H0={self.reuse[('auau', 'H0')]}",
            ],
            "output_root": self.root / "models",
            "json_out": self.root / "plan.json",
            "python": "/usr/bin/python3",
            "runner": SCRIPT.with_name("run_the134_h70_model.py"),
            "max_parallel_views": 7,
        }
        values.update(overrides)
        return SimpleNamespace(**values)


class ExecutionPlanTests(unittest.TestCase):
    def test_exact_fourteen_lane_plan_and_reuse_routing(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = Fixture(Path(directory))
            plan = planner.build_plan(fixture.args())
            self.assertEqual(plan["schema"], planner.SCHEMA)
            self.assertEqual(plan["model_count"], 14)
            self.assertFalse(plan["execution_performed"])
            self.assertEqual(
                [lane["lane_id"] for lane in plan["lanes"]],
                [
                    f"{system}:{view}"
                    for system in planner.SYSTEMS
                    for view in planner.ALL_SHOWER_VIEWS
                ],
            )
            reused = {
                lane["lane_id"]
                for lane in plan["lanes"]
                if lane["reuse_audit"] is not None
            }
            self.assertEqual(reused, {"pp:H70", "auau:H0"})
            self.assertTrue(all(lane["command"][-1] == "--execute" for lane in plan["lanes"]))

    def test_plan_write_is_idempotent_and_drift_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = Fixture(Path(directory))
            args = fixture.args()
            plan = planner.build_plan(args)
            planner.write_json_idempotent(args.json_out, plan)
            planner.write_json_idempotent(args.json_out, plan)
            changed = dict(plan)
            changed["max_parallel_views"] = 1
            with self.assertRaisesRegex(planner.PlanError, "different content"):
                planner.write_json_idempotent(args.json_out, changed)

    def test_missing_view_hash_drift_and_early_authority_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = Fixture(Path(directory))
            manifest_path = fixture.manifests["pp"]
            manifest = json.loads(manifest_path.read_text())
            manifest["audits"].pop("R70")
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            with self.assertRaisesRegex(planner.PlanError, "audit inventory mismatch"):
                planner.build_plan(fixture.args())

        with tempfile.TemporaryDirectory() as directory:
            fixture = Fixture(Path(directory))
            Path(json.loads(fixture.manifests["auau"].read_text())["matrix"]).write_bytes(b"drift")
            with self.assertRaisesRegex(planner.PlanError, "SHA-256 drift"):
                planner.build_plan(fixture.args())

        with tempfile.TemporaryDirectory() as directory:
            fixture = Fixture(Path(directory))
            manifest = json.loads(fixture.manifests["pp"].read_text())
            audit_path = Path(manifest["audits"]["H70"]["path"])
            audit = json.loads(audit_path.read_text())
            audit["full_training_authority"] = 0
            audit_path.write_text(json.dumps(audit), encoding="utf-8")
            manifest["audits"]["H70"]["sha256"] = sha256(audit_path)
            fixture.manifests["pp"].write_text(json.dumps(manifest), encoding="utf-8")
            with self.assertRaisesRegex(planner.PlanError, "full_training_authority"):
                planner.build_plan(fixture.args())

    def test_reuse_inventory_parallel_cap_and_output_root_are_strict(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = Fixture(Path(directory))
            with self.assertRaisesRegex(planner.PlanError, "reuse audit inventory mismatch"):
                planner.build_plan(fixture.args(reuse_audit=[]))
            with self.assertRaisesRegex(planner.PlanError, "between 1 and 7"):
                planner.build_plan(fixture.args(max_parallel_views=8))
            with self.assertRaisesRegex(planner.PlanError, "absolute path"):
                planner.build_plan(fixture.args(output_root=Path("relative")))
            fixture.args().output_root.mkdir()
            with self.assertRaisesRegex(planner.PlanError, "already exists"):
                planner.build_plan(fixture.args())
            with self.assertRaisesRegex(planner.PlanError, "missing or unreadable"):
                planner.build_plan(
                    fixture.args(
                        python="/does/not/exist/python",
                        output_root=fixture.root / "different-models",
                    )
                )


if __name__ == "__main__":
    unittest.main()
