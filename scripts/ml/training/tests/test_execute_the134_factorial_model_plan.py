#!/usr/bin/env python3

from __future__ import annotations

import importlib.util
import hashlib
import json
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve()
TRAINING = HERE.parents[1]


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


planner = load("the134_plan_for_executor_test", TRAINING / "build_the134_factorial_execution_plan.py")
executor = load("the134_executor_test", TRAINING / "execute_the134_factorial_model_plan.py")
fixture_module = load(
    "the134_plan_fixture_for_executor_test",
    HERE.with_name("test_build_the134_factorial_execution_plan.py"),
)


class FactorialExecutorTests(unittest.TestCase):
    def make_plan(self, root: Path) -> tuple[Path, dict]:
        fixture = fixture_module.Fixture(root)
        plan = planner.build_plan(fixture.args(max_parallel_views=2))
        path = root / "execution-plan.json"
        path.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
        return path, plan

    def test_exact_plan_admits_one_system_without_execution(self):
        with tempfile.TemporaryDirectory() as directory:
            path, _ = self.make_plan(Path(directory))
            plan = executor.validate_plan(path, hashlib.sha256(path.read_bytes()).hexdigest())
            lanes, parallel = executor.admit_execution(
                plan,
                system="pp",
                requested_parallel=None,
                receipt_dir=None,
                execute=False,
            )
            self.assertEqual(len(lanes), 7)
            self.assertEqual(parallel, 2)
            self.assertTrue(all(lane["system"] == "pp" for lane in lanes))

    def test_hash_drift_and_parallelism_widening_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path, plan = self.make_plan(root)
            Path(plan["matrix_manifests"]["pp"]["matrix"]).write_bytes(b"drift")
            with self.assertRaisesRegex(executor.ExecutionError, "SHA-256 drift"):
                executor.validate_plan(path, hashlib.sha256(path.read_bytes()).hexdigest())

        with tempfile.TemporaryDirectory() as directory:
            path, _ = self.make_plan(Path(directory))
            plan = executor.validate_plan(path, hashlib.sha256(path.read_bytes()).hexdigest())
            with self.assertRaisesRegex(executor.ExecutionError, "cannot widen"):
                executor.admit_execution(
                    plan,
                    system="auau",
                    requested_parallel=3,
                    receipt_dir=None,
                    execute=False,
                )

    def test_login_host_existing_output_and_receipt_location_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path, _ = self.make_plan(root)
            plan = executor.validate_plan(path, hashlib.sha256(path.read_bytes()).hexdigest())
            with self.assertRaisesRegex(executor.ExecutionError, "login host"):
                executor.admit_execution(
                    plan,
                    system="pp",
                    requested_parallel=1,
                    receipt_dir=root / "receipts",
                    execute=True,
                    hostname="sphnxuser05.sdcc.bnl.gov",
                )
            with self.assertRaisesRegex(executor.ExecutionError, "outside"):
                executor.admit_execution(
                    plan,
                    system="pp",
                    requested_parallel=1,
                    receipt_dir=Path(plan["output_root"]) / "receipts",
                    execute=True,
                    hostname="worker01",
                )
            Path(plan["output_root"]).mkdir()
            with self.assertRaisesRegex(executor.ExecutionError, "output_root already exists"):
                executor.validate_plan(path, hashlib.sha256(path.read_bytes()).hexdigest())

    def test_plan_hash_is_required_before_internal_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            path, _ = self.make_plan(Path(directory))
            expected = hashlib.sha256(path.read_bytes()).hexdigest()
            payload = json.loads(path.read_text())
            payload["max_parallel_views"] = 1
            path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
            with self.assertRaisesRegex(executor.ExecutionError, "SHA-256 drift"):
                executor.validate_plan(path, expected)

    def test_first_bad_stops_new_lane_starts_and_never_retries(self):
        lanes = [
            {"lane_id": f"pp:V{index}", "command": ["unused"]}
            for index in range(5)
        ]
        calls: list[str] = []

        def run(lane):
            calls.append(lane["lane_id"])
            failed = lane["lane_id"] == "pp:V0"
            return {
                "lane_id": lane["lane_id"],
                "status": "FAIL" if failed else "PASS",
                "exit_code": 1 if failed else 0,
            }

        results, unstarted = executor.execute_lanes(lanes, 1, run)
        self.assertEqual(calls, ["pp:V0"])
        self.assertEqual([result["lane_id"] for result in results], ["pp:V0"])
        self.assertEqual(unstarted, ["pp:V1", "pp:V2", "pp:V3", "pp:V4"])


if __name__ == "__main__":
    unittest.main()
