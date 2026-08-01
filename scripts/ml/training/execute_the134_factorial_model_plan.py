#!/usr/bin/env python3
"""Execute one system from a certified THE-134 factorial model plan.

This is a bounded local/worker executor, never a submitter.  It reopens and
rehashes the complete plan surface, runs at most the plan-authorized number of
views, starts no new view after a first failure, and never retries.  Known
SDCC login hosts are rejected so model training cannot accidentally execute
as login-node work.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import socket
import subprocess
import sys
import time
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Callable

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import build_the134_factorial_execution_plan as planner  # noqa: E402


RECEIPT_SCHEMA = "THE134_FACTORIAL_MODEL_EXECUTION_RECEIPT_V1"
LOGIN_HOST_RE = re.compile(r"^(?:sphnxuser[0-9]+|login[0-9]*|ssh)$", re.IGNORECASE)


class ExecutionError(RuntimeError):
    """Fail-closed execution admission or runtime error."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.is_file():
        raise ExecutionError(f"{label} is missing: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ExecutionError(f"{label} is not readable JSON: {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ExecutionError(f"{label} must be a JSON object: {path}")
    return payload


def require_registered_file(path_value: object, sha_value: object, label: str) -> Path:
    try:
        path = Path(str(path_value)).expanduser().resolve(strict=True)
    except OSError as exc:
        raise ExecutionError(f"{label} is missing or unreadable: {path_value}") from exc
    if not path.is_file():
        raise ExecutionError(f"{label} is not a regular file: {path}")
    expected = planner.require_sha256(sha_value, f"{label} SHA-256")
    observed = sha256_file(path)
    if observed != expected:
        raise ExecutionError(
            f"{label} SHA-256 drift: expected={expected} observed={observed}"
        )
    return path


def command_option(command: list[str], option: str) -> str:
    if command.count(option) != 1:
        raise ExecutionError(f"command must contain {option} exactly once")
    index = command.index(option)
    if index + 1 >= len(command):
        raise ExecutionError(f"command option lacks a value: {option}")
    return command[index + 1]


def validate_plan(plan_path: Path, expected_sha256: str) -> dict[str, Any]:
    expected_plan_sha = planner.require_sha256(
        expected_sha256, "factorial execution plan SHA-256"
    )
    observed_plan_sha = sha256_file(plan_path)
    if observed_plan_sha != expected_plan_sha:
        raise ExecutionError(
            "factorial execution plan SHA-256 drift: "
            f"expected={expected_plan_sha} observed={observed_plan_sha}"
        )
    plan = read_json(plan_path, "factorial execution plan")
    checks = {
        "schema": plan.get("schema") == planner.SCHEMA,
        "status": plan.get("status") == "READY_NOT_EXECUTED",
        "execution_performed": plan.get("execution_performed") is False,
        "promotion": plan.get("promotion_status")
        == "CANDIDATE_CURRENT_NOT_CANONICAL",
        "systems": plan.get("systems") == list(planner.SYSTEMS),
        "views": plan.get("views") == list(planner.ALL_SHOWER_VIEWS),
        "model_count": plan.get("model_count") == 14,
        "parallelism_scope": plan.get("parallelism_scope")
        == "WITHIN_ONE_CERTIFIED_SYSTEM_AGGREGATE",
    }
    if not all(checks.values()):
        failed = sorted(name for name, passed in checks.items() if not passed)
        raise ExecutionError(f"factorial execution plan failed authority gates: {failed}")

    max_parallel = plan.get("max_parallel_views")
    if not isinstance(max_parallel, int) or not 1 <= max_parallel <= 7:
        raise ExecutionError("plan max_parallel_views must be an integer in [1, 7]")

    output_root = Path(str(plan.get("output_root", "")))
    if not output_root.is_absolute() or output_root == Path("/"):
        raise ExecutionError("plan output_root must be absolute and non-root")
    if output_root.exists():
        raise ExecutionError(
            f"plan output_root already exists; refusing possible duplicate work: {output_root}"
        )

    python = Path(str(plan.get("python", "")))
    if not python.is_file() or not os.access(python, os.X_OK):
        raise ExecutionError(f"plan Python runtime is not executable: {python}")
    runner_record = plan.get("runner")
    if not isinstance(runner_record, dict):
        raise ExecutionError("plan runner record is invalid")
    runner = require_registered_file(
        runner_record.get("path"), runner_record.get("sha256"), "model runner"
    )

    matrices = plan.get("matrix_manifests")
    if not isinstance(matrices, dict) or set(matrices) != set(planner.SYSTEMS):
        raise ExecutionError("plan matrix manifest inventory differs")
    for system in planner.SYSTEMS:
        record = matrices[system]
        if not isinstance(record, dict):
            raise ExecutionError(f"{system} matrix record is invalid")
        require_registered_file(
            record.get("manifest"),
            record.get("manifest_sha256"),
            f"{system} matrix manifest",
        )
        require_registered_file(
            record.get("matrix"), record.get("matrix_sha256"), f"{system} matrix"
        )
        audits = record.get("audits")
        if not isinstance(audits, dict) or set(audits) != set(
            planner.ALL_SHOWER_VIEWS
        ):
            raise ExecutionError(f"{system} extraction audit inventory differs")
        for view in planner.ALL_SHOWER_VIEWS:
            audit = audits[view]
            if not isinstance(audit, dict):
                raise ExecutionError(f"{system}/{view} extraction audit is invalid")
            require_registered_file(
                audit.get("path"),
                audit.get("sha256"),
                f"{system}/{view} extraction audit",
            )

    reuse = plan.get("reuse_authorities")
    expected_reuse = {f"{system}:{view}" for system, view in planner.REQUIRED_REUSE}
    if not isinstance(reuse, dict) or set(reuse) != expected_reuse:
        raise ExecutionError("plan reuse authority inventory differs")
    for lane_id, record in reuse.items():
        if not isinstance(record, dict):
            raise ExecutionError(f"{lane_id} reuse authority is invalid")
        require_registered_file(
            record.get("path"), record.get("sha256"), f"{lane_id} reuse authority"
        )

    lanes = plan.get("lanes")
    expected_lane_ids = [
        f"{system}:{view}"
        for system in planner.SYSTEMS
        for view in planner.ALL_SHOWER_VIEWS
    ]
    if not isinstance(lanes, list) or [lane.get("lane_id") for lane in lanes] != expected_lane_ids:
        raise ExecutionError("plan lane inventory/order differs")

    for lane in lanes:
        system = lane.get("system")
        view = lane.get("view")
        lane_id = lane.get("lane_id")
        if lane.get("status") != "READY_NOT_EXECUTED":
            raise ExecutionError(f"{lane_id} is not ready-not-executed")
        if lane.get("model_origin") != planner.expected_model_origin(system, view):
            raise ExecutionError(f"{lane_id} model origin differs")
        expected_output = output_root / str(system) / str(view).lower()
        if Path(str(lane.get("output_directory", ""))) != expected_output:
            raise ExecutionError(f"{lane_id} output directory differs")
        if expected_output.exists():
            raise ExecutionError(f"{lane_id} output already exists")
        command = lane.get("command")
        if not isinstance(command, list) or not all(
            isinstance(item, str) and item for item in command
        ):
            raise ExecutionError(f"{lane_id} command is invalid")
        command_checks = {
            "python": command[0] == str(python),
            "runner": len(command) > 1 and command[1] == str(runner),
            "system": command_option(command, "--system") == system,
            "view": command_option(command, "--view") == view,
            "matrix": command_option(command, "--matrix")
            == matrices[system]["matrix"],
            "audit": command_option(command, "--extraction-audit")
            == matrices[system]["audits"][view]["path"],
            "output": command_option(command, "--outdir") == str(expected_output),
            "execute": command[-1] == "--execute" and command.count("--execute") == 1,
        }
        reuse_record = lane.get("reuse_audit")
        if lane_id in expected_reuse:
            command_checks["reuse"] = (
                isinstance(reuse_record, dict)
                and command_option(command, "--reuse-audit") == reuse[lane_id]["path"]
            )
        else:
            command_checks["reuse"] = reuse_record is None and "--reuse-audit" not in command
        if not all(command_checks.values()):
            failed = sorted(name for name, passed in command_checks.items() if not passed)
            raise ExecutionError(f"{lane_id} command binding failed: {failed}")

    try:
        rebuilt = planner.build_plan(
            SimpleNamespace(
                pp_matrix_manifest=Path(matrices["pp"]["manifest"]),
                auau_matrix_manifest=Path(matrices["auau"]["manifest"]),
                reuse_audit=[
                    f"{lane_id}={reuse[lane_id]['path']}"
                    for lane_id in sorted(expected_reuse)
                ],
                output_root=output_root,
                json_out=plan_path,
                python=str(python),
                runner=runner,
                max_parallel_views=max_parallel,
            )
        )
    except (OSError, planner.PlanError) as exc:
        raise ExecutionError(f"factorial execution plan cannot be rebuilt: {exc}") from exc
    if rebuilt != plan:
        raise ExecutionError("factorial execution plan differs from deterministic rebuild")
    return plan


def admit_execution(
    plan: dict[str, Any],
    *,
    system: str,
    requested_parallel: int | None,
    receipt_dir: Path | None,
    execute: bool,
    hostname: str | None = None,
) -> tuple[list[dict[str, Any]], int]:
    if system not in planner.SYSTEMS:
        raise ExecutionError(f"unsupported system: {system}")
    plan_parallel = int(plan["max_parallel_views"])
    parallel = plan_parallel if requested_parallel is None else requested_parallel
    if not isinstance(parallel, int) or parallel < 1 or parallel > plan_parallel:
        raise ExecutionError(
            f"requested parallelism must be in [1, {plan_parallel}] and cannot widen the plan"
        )
    selected = [lane for lane in plan["lanes"] if lane["system"] == system]
    if len(selected) != 7:
        raise ExecutionError(f"{system} plan must contain exactly seven lanes")
    if not execute:
        return selected, parallel

    short_host = (hostname or socket.gethostname()).split(".", 1)[0]
    if LOGIN_HOST_RE.fullmatch(short_host):
        raise ExecutionError(f"execution on SDCC login host is forbidden: {short_host}")
    if receipt_dir is None or not receipt_dir.is_absolute():
        raise ExecutionError("execute mode requires an absolute receipt directory")
    receipt_dir = Path(os.path.normpath(str(receipt_dir)))
    output_root = Path(plan["output_root"])
    if receipt_dir == output_root or receipt_dir.is_relative_to(output_root):
        raise ExecutionError("receipt directory must be outside the model output root")
    if receipt_dir.exists():
        raise ExecutionError(f"receipt directory already exists: {receipt_dir}")
    return selected, parallel


def execute_lanes(
    lanes: list[dict[str, Any]],
    max_workers: int,
    run_lane: Callable[[dict[str, Any]], dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[str]]:
    pending = iter(lanes)
    active: dict[Future[dict[str, Any]], dict[str, Any]] = {}
    results: dict[str, dict[str, Any]] = {}
    failed = False

    def submit_next(pool: ThreadPoolExecutor) -> bool:
        try:
            lane = next(pending)
        except StopIteration:
            return False
        active[pool.submit(run_lane, lane)] = lane
        return True

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        for _ in range(min(max_workers, len(lanes))):
            submit_next(pool)
        while active:
            completed, _ = wait(active, return_when=FIRST_COMPLETED)
            for future in sorted(completed, key=lambda item: active[item]["lane_id"]):
                lane = active.pop(future)
                try:
                    result = future.result()
                except Exception as exc:  # preserve a bounded controller failure
                    result = {
                        "lane_id": lane["lane_id"],
                        "status": "CONTROLLER_EXCEPTION",
                        "exit_code": 125,
                        "error": f"{type(exc).__name__}: {exc}",
                    }
                results[lane["lane_id"]] = result
                if result.get("exit_code") != 0:
                    failed = True
            while not failed and len(active) < max_workers and submit_next(pool):
                pass

    ordered = [results[lane["lane_id"]] for lane in lanes if lane["lane_id"] in results]
    unstarted = [lane["lane_id"] for lane in lanes if lane["lane_id"] not in results]
    return ordered, unstarted


def run_lane_subprocess(lane: dict[str, Any], logs_root: Path) -> dict[str, Any]:
    lane_slug = lane["lane_id"].replace(":", "-").lower()
    stdout_path = logs_root / f"{lane_slug}.stdout.log"
    stderr_path = logs_root / f"{lane_slug}.stderr.log"
    started = datetime.now(timezone.utc).isoformat()
    start = time.monotonic()
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        completed = subprocess.run(
            lane["command"], stdout=stdout, stderr=stderr, check=False
        )
    return {
        "lane_id": lane["lane_id"],
        "status": "PASS" if completed.returncode == 0 else "FAIL",
        "exit_code": completed.returncode,
        "started_at": started,
        "elapsed_seconds": round(time.monotonic() - start, 6),
        "stdout": str(stdout_path),
        "stdout_sha256": sha256_file(stdout_path),
        "stderr": str(stderr_path),
        "stderr_sha256": sha256_file(stderr_path),
    }


def write_receipt(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True, type=Path)
    parser.add_argument("--plan-sha256", required=True)
    parser.add_argument("--system", required=True, choices=planner.SYSTEMS)
    parser.add_argument("--receipt-dir", type=Path)
    parser.add_argument("--max-parallel-views", type=int)
    parser.add_argument("--execute", action="store_true")
    return parser.parse_args()


def main() -> int:
    try:
        args = parse_args()
        plan_path = args.plan.expanduser().resolve(strict=True)
        plan = validate_plan(plan_path, args.plan_sha256)
        receipt_dir = args.receipt_dir.expanduser() if args.receipt_dir else None
        selected, parallel = admit_execution(
            plan,
            system=args.system,
            requested_parallel=args.max_parallel_views,
            receipt_dir=receipt_dir,
            execute=args.execute,
        )
        if not args.execute:
            print(
                json.dumps(
                    {
                        "schema": RECEIPT_SCHEMA,
                        "status": "READY_NOT_EXECUTED",
                        "system": args.system,
                        "lane_count": len(selected),
                        "max_parallel_views": parallel,
                        "plan_sha256": sha256_file(plan_path),
                    },
                    sort_keys=True,
                )
            )
            return 0

        assert receipt_dir is not None
        receipt_dir = Path(os.path.normpath(str(receipt_dir)))
        logs_root = receipt_dir / "logs"
        logs_root.mkdir(parents=True)
        results, unstarted = execute_lanes(
            selected,
            parallel,
            lambda lane: run_lane_subprocess(lane, logs_root),
        )
        passed = len(results) == 7 and not unstarted and all(
            result.get("exit_code") == 0 for result in results
        )
        receipt = {
            "schema": RECEIPT_SCHEMA,
            "status": "PASS" if passed else "FAIL_FIRST_BAD_PRESERVED",
            "execution_performed": True,
            "automatic_retry_performed": False,
            "system": args.system,
            "max_parallel_views": parallel,
            "plan": str(plan_path),
            "plan_sha256": sha256_file(plan_path),
            "results": results,
            "unstarted_after_first_bad": unstarted,
            "promotion_status": "CANDIDATE_CURRENT_NOT_CANONICAL",
        }
        receipt_path = receipt_dir / "execution_receipt.json"
        write_receipt(receipt_path, receipt)
        print(receipt_path)
        return 0 if passed else 2
    except (ExecutionError, OSError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    raise SystemExit(main())
