#!/usr/bin/env python3
"""Build a compact AuAu RecoilJets campaign timing/coverage checkpoint.

This helper is intentionally file-coverage-first. Condor state is useful for
monitoring, but final campaign acceptance is based on usable ROOT products.
Run it on an SDCC node, or anywhere the /sphenix paths are visible.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Any


ROOT_NAME_TEMPLATE = (
    "RecoilJets_isAuAu_"
    "preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_"
    "nonTightAuAuBDTComplement_baseVariant_run{run}_grp{grp}.root"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-id", required=True)
    parser.add_argument("--expected-list-root", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--cluster", action="append", default=[])
    parser.add_argument("--group-size", type=int)
    parser.add_argument("--request-memory-mb", type=int)
    parser.add_argument("--submitted-files", type=int)
    parser.add_argument("--submitted-runs", type=int)
    parser.add_argument("--non-tiny-threshold", type=int, default=50_000)
    parser.add_argument(
        "--recovery-target",
        action="append",
        default=[],
        help="Expected recovery target as RUN:GRP, for example 00076220:010.",
    )
    parser.add_argument("--previous-json")
    parser.add_argument("--write-json")
    parser.add_argument("--json", action="store_true", help="Print JSON only.")
    parser.add_argument("--run-condor", action="store_true")
    parser.add_argument("--include-history", action="store_true")
    parser.add_argument("--condor-timeout", type=int, default=20)
    return parser.parse_args()


def parse_expected_lists(list_root: Path) -> list[tuple[str, str, Path]]:
    expected: list[tuple[str, str, Path]] = []
    pattern = re.compile(r"run(\d+)_grp(\d+)\.list$")
    for path in sorted(list_root.glob("run*_grp*.list")):
        match = pattern.search(path.name)
        if not match:
            continue
        expected.append((match.group(1), match.group(2), path))
    return expected


def output_path(output_root: Path, run: str, grp: str) -> Path:
    filename = ROOT_NAME_TEMPLATE.format(run=run, grp=grp)
    canonical = output_root / run / filename
    if canonical.exists():
        return canonical
    legacy = output_root / f"run{run}" / filename
    if legacy.exists():
        return legacy
    return canonical


def classify_outputs(
    expected: list[tuple[str, str, Path]],
    output_root: Path,
    non_tiny_threshold: int,
    recovery_targets: set[tuple[str, str]],
) -> dict[str, Any]:
    by_run: dict[str, dict[str, int]] = {}
    recovery_status: list[dict[str, Any]] = []
    max_size = 0
    total_bytes = 0
    newest_mtime = 0.0
    counts = {"non_tiny": 0, "tiny": 0, "missing": 0}

    for run, grp, _ in expected:
        path = output_path(output_root, run, grp)
        try:
            stat = path.stat()
            size = stat.st_size
            newest_mtime = max(newest_mtime, stat.st_mtime)
        except OSError:
            size = -1
        status = "missing" if size < 0 else ("non_tiny" if size > non_tiny_threshold else "tiny")
        counts[status] += 1
        if size > 0:
            total_bytes += size
            max_size = max(max_size, size)

        run_counts = by_run.setdefault(run, {"total": 0, "non_tiny": 0, "tiny": 0, "missing": 0})
        run_counts["total"] += 1
        run_counts[status] += 1

        if (run, grp) in recovery_targets:
            recovery_status.append(
                {
                    "run": run,
                    "grp": grp,
                    "size_bytes": size,
                    "status": status,
                    "path": str(path),
                }
            )

    complete_runs = sum(1 for item in by_run.values() if item["non_tiny"] == item["total"])
    partial_runs = sum(1 for item in by_run.values() if 0 < item["non_tiny"] < item["total"])
    zero_runs = sum(1 for item in by_run.values() if item["non_tiny"] == 0)
    return {
        "expected_jobs": len(expected),
        "expected_runs": len(by_run),
        "jobs_non_tiny": counts["non_tiny"],
        "jobs_tiny": counts["tiny"],
        "jobs_missing": counts["missing"],
        "runs_complete_non_tiny": complete_runs,
        "runs_partial": partial_runs,
        "runs_zero_non_tiny": zero_runs,
        "total_bytes": total_bytes,
        "max_bytes": max_size,
        "newest_output_mtime": newest_mtime or None,
        "recovery_targets": sorted(recovery_status, key=lambda x: (x["run"], x["grp"])),
    }


def run_command(cmd: list[str], timeout: int) -> tuple[int, str]:
    try:
        proc = subprocess.run(
            cmd,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        output = exc.stdout or ""
        if isinstance(output, bytes):
            output = output.decode(errors="replace")
        return 124, output + f"\nTIMEOUT after {timeout}s"
    return proc.returncode, proc.stdout


def condor_queue(clusters: list[str], timeout: int) -> dict[str, Any]:
    if not clusters:
        return {}
    cluster_set = set(clusters)
    rc, output = run_command(
        [
            "condor_q",
            os.environ.get("USER", ""),
            "-af",
            "ClusterId",
            "ProcId",
            "JobStatus",
            "HoldReasonCode",
            "ExitCode",
            "RemoteWallClockTime",
            "ResidentSetSize_RAW",
            "Args",
        ],
        timeout,
    )
    summary: dict[str, Any] = {
        "returncode": rc,
        "raw_error": output if rc not in (0, 1) else "",
        "clusters": {},
        "held_rows": [],
    }
    for line in output.splitlines():
        parts = line.split(maxsplit=7)
        if len(parts) < 3 or parts[0] not in cluster_set:
            continue
        cluster, proc, status = parts[:3]
        cluster_summary = summary["clusters"].setdefault(
            cluster, {"total": 0, "idle": 0, "running": 0, "held": 0, "other": 0}
        )
        cluster_summary["total"] += 1
        if status == "1":
            cluster_summary["idle"] += 1
        elif status == "2":
            cluster_summary["running"] += 1
        elif status == "5":
            cluster_summary["held"] += 1
            summary["held_rows"].append({"cluster": cluster, "proc": proc, "line": line})
        else:
            cluster_summary["other"] += 1
    return summary


def condor_history(clusters: list[str], timeout: int) -> dict[str, Any]:
    if not clusters:
        return {}
    result: dict[str, Any] = {}
    for cluster in clusters:
        rc, output = run_command(
            [
                "condor_history",
                "-constraint",
                f"ClusterId=={cluster}",
                "-af",
                "ClusterId",
                "ProcId",
                "ExitCode",
                "RemoteWallClockTime",
                "RequestMemory",
                "MemoryUsage",
                "ResidentSetSize_RAW",
            ],
            timeout,
        )
        rows = 0
        exit0 = 0
        exit2 = 0
        other_exit = 0
        runtimes: list[float] = []
        memory_values: list[float] = []
        for line in output.splitlines():
            parts = line.split()
            if len(parts) < 4 or parts[0] != cluster:
                continue
            rows += 1
            try:
                exit_code = int(parts[2])
            except ValueError:
                exit_code = -999
            if exit_code == 0:
                exit0 += 1
            elif exit_code == 2:
                exit2 += 1
            else:
                other_exit += 1
            try:
                runtime = float(parts[3])
                if runtime > 0:
                    runtimes.append(runtime)
            except ValueError:
                pass
            for token in parts[4:]:
                try:
                    memory = float(token)
                except ValueError:
                    continue
                if memory > 0:
                    memory_values.append(memory)
                    break
        result[cluster] = {
            "returncode": rc,
            "history_rows": rows,
            "exit0": exit0,
            "exit2": exit2,
            "other_exit": other_exit,
            "runtime_mean_s": sum(runtimes) / len(runtimes) if runtimes else None,
            "runtime_max_s": max(runtimes) if runtimes else None,
            "memory_max_reported": max(memory_values) if memory_values else None,
            "raw_error": output if rc not in (0, 1) else "",
        }
    return result


def parse_recovery_targets(items: list[str]) -> set[tuple[str, str]]:
    targets: set[tuple[str, str]] = set()
    for item in items:
        if ":" not in item:
            raise SystemExit(f"invalid --recovery-target {item!r}; expected RUN:GRP")
        run, grp = item.split(":", 1)
        targets.add((run.strip(), grp.strip().zfill(3)))
    return targets


def load_previous(path: str | None) -> dict[str, Any] | None:
    if not path:
        return None
    with Path(path).open("r", encoding="utf-8") as handle:
        return json.load(handle)


def add_delta_and_eta(payload: dict[str, Any], previous: dict[str, Any] | None) -> None:
    if not previous:
        return
    current_time = payload["timestamp_epoch"]
    previous_time = previous.get("timestamp_epoch")
    if not isinstance(previous_time, (int, float)) or current_time <= previous_time:
        return
    coverage = payload["coverage"]
    previous_coverage = previous.get("coverage") if isinstance(previous.get("coverage"), dict) else {}
    elapsed = current_time - previous_time
    delta_non_tiny = coverage["jobs_non_tiny"] - int(previous_coverage.get("jobs_non_tiny", 0))
    delta_bytes = coverage["total_bytes"] - int(previous_coverage.get("total_bytes", 0))
    remaining = max(coverage["expected_jobs"] - coverage["jobs_non_tiny"], 0)
    jobs_per_hour = delta_non_tiny / elapsed * 3600.0 if elapsed > 0 else 0.0
    eta_hours = remaining / jobs_per_hour if jobs_per_hour > 0 else None
    payload["delta_from_previous"] = {
        "elapsed_seconds": elapsed,
        "delta_non_tiny_jobs": delta_non_tiny,
        "delta_total_bytes": delta_bytes,
        "jobs_per_hour": jobs_per_hour,
        "eta_hours_from_output_growth": eta_hours,
        "eta_caution": "Output-growth ETA is a lower-confidence estimate for long-tail AuAu jobs.",
    }


def render_text(payload: dict[str, Any]) -> str:
    coverage = payload["coverage"]
    lines = [
        f"campaign_id={payload['campaign_id']}",
        f"timestamp={payload['timestamp_iso']}",
        (
            "coverage "
            f"jobs={coverage['jobs_non_tiny']}/{coverage['expected_jobs']} non_tiny "
            f"tiny={coverage['jobs_tiny']} missing={coverage['jobs_missing']} "
            f"runs_complete={coverage['runs_complete_non_tiny']}/{coverage['expected_runs']} "
            f"runs_partial={coverage['runs_partial']} runs_zero={coverage['runs_zero_non_tiny']}"
        ),
        f"bytes total={coverage['total_bytes']} max={coverage['max_bytes']}",
    ]
    if payload.get("queue"):
        lines.append("queue " + json.dumps(payload["queue"].get("clusters", {}), sort_keys=True))
        if payload["queue"].get("held_rows"):
            lines.append(f"held_rows={len(payload['queue']['held_rows'])}")
    if payload.get("delta_from_previous"):
        delta = payload["delta_from_previous"]
        eta = delta["eta_hours_from_output_growth"]
        eta_text = "unavailable" if eta is None else f"{eta:.2f} h"
        lines.append(
            "delta "
            f"dt={delta['elapsed_seconds']:.0f}s "
            f"d_non_tiny={delta['delta_non_tiny_jobs']} "
            f"jobs_per_hour={delta['jobs_per_hour']:.2f} "
            f"eta={eta_text}"
        )
    if coverage["recovery_targets"]:
        ready = sum(1 for item in coverage["recovery_targets"] if item["status"] == "non_tiny")
        lines.append(f"recovery_targets_ready={ready}/{len(coverage['recovery_targets'])}")
    return "\n".join(lines)


def main() -> int:
    args = parse_args()
    list_root = Path(args.expected_list_root)
    output_root = Path(args.output_root)
    expected = parse_expected_lists(list_root)
    if not expected:
        print(f"ERROR: no run*_grp*.list files found under {list_root}", file=sys.stderr)
        return 2

    payload: dict[str, Any] = {
        "campaign_id": args.campaign_id,
        "timestamp_epoch": time.time(),
        "timestamp_iso": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "inputs": {
            "expected_list_root": str(list_root),
            "output_root": str(output_root),
            "clusters": args.cluster,
            "group_size": args.group_size,
            "request_memory_mb": args.request_memory_mb,
            "submitted_files": args.submitted_files,
            "submitted_runs": args.submitted_runs,
            "non_tiny_threshold": args.non_tiny_threshold,
        },
        "coverage": classify_outputs(
            expected,
            output_root,
            args.non_tiny_threshold,
            parse_recovery_targets(args.recovery_target),
        ),
    }
    if args.run_condor:
        payload["queue"] = condor_queue(args.cluster, args.condor_timeout)
    if args.include_history:
        payload["history"] = condor_history(args.cluster, args.condor_timeout)

    add_delta_and_eta(payload, load_previous(args.previous_json))

    if args.write_json:
        path = Path(args.write_json)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(render_text(payload))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
