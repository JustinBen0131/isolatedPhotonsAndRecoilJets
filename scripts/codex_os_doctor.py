#!/usr/bin/env python3
"""Run hardening checks for the ThesisAnalysis Codex operating system."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from codex_work_register_common import (
    DEFAULT_REGISTER,
    RegisterError,
    first_line,
    load_register,
    parse_when,
    validate_register,
)


LIVE_STATUSES = {"active", "running", "waiting", "blocked", "review"}
DONE_STATUSES = {"done_pending_review", "archived"}

BASE_REQUIRED_FILES = [
    Path("AGENTS.md"),
    Path("agent_context/CODEX_WORK_REGISTER.yaml"),
    Path("agent_context/ARTIFACT_REGISTRY.yaml"),
    Path("agent_context/THESIS_NARRATIVE_MAP.md"),
    Path("agent_context/policies/CODEX_OPERATING_SYSTEM.md"),
    Path("agent_context/policies/AGENTIC_OS_HARDENING.md"),
    Path("agent_context/policies/AGENTIC_OS_DREAMING.md"),
    Path("agent_context/policies/ASK_CHATGPT_DELEGATION.md"),
    Path("agent_context/policies/HARD_STOPS_AND_SAFETY.md"),
    Path("agent_context/policies/DUPLICATE_RUN_GUARD.md"),
    Path("agent_context/policies/MEMORY_AND_STATUS.md"),
    Path("agent_context/policies/LOAD_MAP.yaml"),
    Path("agent_context/templates/OS_POSTMORTEM_TEMPLATE.md"),
    Path("scripts/codex_artifact_registry.py"),
    Path("scripts/codex_context_pack.py"),
    Path("scripts/codex_os_guard.py"),
    Path("scripts/codex_os_dream.py"),
    Path("scripts/codex_chatgpt_research_pack.py"),
    Path("scripts/codex_os_snapshot.py"),
    Path("scripts/codex_thesis_radar.py"),
]

STRICT_REQUIRED_FILES = [
    Path("agent_context/local/os_events.jsonl"),
]


@dataclass
class Finding:
    severity: str
    message: str


def text_of(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except OSError:
        return ""


def add(findings: list[Finding], severity: str, message: str) -> None:
    findings.append(Finding(severity=severity, message=message))


def check_required_files(findings: list[Finding], profile: str) -> None:
    required_files = list(BASE_REQUIRED_FILES)
    if profile in {"strict", "release"}:
        required_files.extend(STRICT_REQUIRED_FILES)
    for path in required_files:
        if not path.exists():
            add(findings, "ERROR", f"missing required OS file: {path}")
        elif path.is_file() and not text_of(path).strip():
            add(findings, "ERROR", f"required OS file is empty: {path}")


def check_policy_routing(findings: list[Finding]) -> None:
    agents = text_of(Path("AGENTS.md"))
    load_map = text_of(Path("agent_context/policies/LOAD_MAP.yaml"))
    os_policy = text_of(Path("agent_context/policies/CODEX_OPERATING_SYSTEM.md"))

    required_agents_terms = [
        "CODEX_OPERATING_SYSTEM.md",
        "HARD_STOPS_AND_SAFETY.md",
        "CODEX_WORK_REGISTER.yaml",
    ]
    for term in required_agents_terms:
        if term not in agents:
            add(findings, "ERROR", f"AGENTS.md no longer references {term}")

    if "agentic_os_hardening" not in load_map:
        add(findings, "ERROR", "LOAD_MAP.yaml lacks agentic_os_hardening route")
    if "AGENTIC_OS_HARDENING.md" not in load_map:
        add(findings, "ERROR", "LOAD_MAP.yaml does not load AGENTIC_OS_HARDENING.md")
    if "agentic_os_dreaming" not in load_map:
        add(findings, "ERROR", "LOAD_MAP.yaml lacks agentic_os_dreaming route")
    if "AGENTIC_OS_DREAMING.md" not in load_map:
        add(findings, "ERROR", "LOAD_MAP.yaml does not load AGENTIC_OS_DREAMING.md")
    if "ask_chatgpt_delegation" not in load_map:
        add(findings, "ERROR", "LOAD_MAP.yaml lacks ask_chatgpt_delegation route")
    if "ASK_CHATGPT_DELEGATION.md" not in load_map:
        add(findings, "ERROR", "LOAD_MAP.yaml does not load ASK_CHATGPT_DELEGATION.md")
    if "codex_os_guard.py" not in load_map:
        add(findings, "WARN", "LOAD_MAP.yaml does not mention codex_os_guard.py")
    if "scripts/codex_os_doctor.py" not in os_policy:
        add(findings, "WARN", "CODEX_OPERATING_SYSTEM.md does not mention codex_os_doctor.py")
    if "codex_os_guard.py" not in os_policy:
        add(findings, "WARN", "CODEX_OPERATING_SYSTEM.md does not mention codex_os_guard.py")
    if "THESIS_NARRATIVE_MAP.md" not in os_policy:
        add(findings, "WARN", "CODEX_OPERATING_SYSTEM.md does not mention thesis spine")
    if "ARTIFACT_REGISTRY.yaml" not in os_policy:
        add(findings, "WARN", "CODEX_OPERATING_SYSTEM.md does not mention artifact registry")
    if "AGENTIC_OS_DREAMING.md" not in os_policy:
        add(findings, "WARN", "CODEX_OPERATING_SYSTEM.md does not mention dream policy")
    if "ASK_CHATGPT_DELEGATION.md" not in os_policy:
        add(findings, "WARN", "CODEX_OPERATING_SYSTEM.md does not mention ChatGPT delegation policy")


def check_register(data: dict[str, Any], findings: list[Finding], now: datetime) -> None:
    for error in validate_register(data):
        add(findings, "ERROR", error)

    daily = data.get("daily_cockpit")
    if not isinstance(daily, dict):
        add(findings, "ERROR", "daily_cockpit must exist in register")
    else:
        deck = daily.get("active_working_point_deck")
        if not isinstance(deck, dict) or not first_line(deck.get("url")):
            add(findings, "ERROR", "daily_cockpit.active_working_point_deck.url is missing")

    workstreams = data.get("workstreams")
    if not isinstance(workstreams, list):
        return

    session_to_workstreams: dict[str, list[str]] = {}
    known_ids = {
        str(item.get("workstream_id"))
        for item in workstreams
        if isinstance(item, dict) and item.get("workstream_id")
    }

    for item in workstreams:
        if not isinstance(item, dict):
            continue
        wid = str(item.get("workstream_id"))
        status = str(item.get("status"))
        next_action = first_line(item.get("current_next_action")).lower()
        evidence = item.get("evidence") if isinstance(item.get("evidence"), list) else []
        active_jobs = item.get("active_jobs") if isinstance(item.get("active_jobs"), list) else []
        active_session = first_line(item.get("active_codex_session"))

        if active_session and status in LIVE_STATUSES:
            session_to_workstreams.setdefault(active_session, []).append(wid)

        if status in LIVE_STATUSES:
            stale_after = parse_when(item.get("stale_after"))
            if stale_after is None:
                add(findings, "ERROR", f"{wid}: live workstream has unparseable stale_after")
            elif stale_after <= now:
                add(findings, "ERROR", f"{wid}: live workstream is stale after {item.get('stale_after')}")
            if not evidence:
                add(findings, "ERROR", f"{wid}: live workstream has no evidence")
            if item.get("linear_sync") == "synced" and not first_line(item.get("linear_issue")):
                add(findings, "ERROR", f"{wid}: linear_sync synced but linear_issue missing")
            completion_phrases = (
                "complete.",
                "completed.",
                "task complete",
                "no standing action",
            )
            if next_action.startswith(completion_phrases) and status not in {"review"}:
                add(findings, "WARN", f"{wid}: live status has completion-like next action")

        if status in DONE_STATUSES:
            if active_jobs:
                add(findings, "ERROR", f"{wid}: done/archived workstream still has active_jobs")
            if active_session:
                add(findings, "ERROR", f"{wid}: done/archived workstream still has active_codex_session")
            if any(word in next_action for word in ("wait", "monitor", "review the", "run ")):
                add(findings, "WARN", f"{wid}: done/archived next_action may still imply work")

        if item.get("linear_sync") == "synced" and not first_line(item.get("linear_issue")):
            add(findings, "ERROR", f"{wid}: synced Linear state but missing issue URL/key")

        deps = item.get("depends_on") if isinstance(item.get("depends_on"), list) else []
        for dep in deps:
            dep_text = str(dep)
            if dep_text.startswith("workstream:"):
                dep_id = dep_text.split(":", 1)[1]
                if dep_id not in known_ids:
                    add(findings, "ERROR", f"{wid}: missing dependency workstream {dep_id}")

    for session, owned in session_to_workstreams.items():
        if len(owned) > 2:
            add(
                findings,
                "WARN",
                f"session {session!r} owns many live workstreams: {', '.join(owned)}",
            )


def check_event_log(findings: list[Finding], profile: str) -> None:
    path = Path("agent_context/local/os_events.jsonl")
    if not path.exists():
        severity = "ERROR" if profile in {"strict", "release"} else "WARN"
        add(findings, severity, "private OS event log is missing: agent_context/local/os_events.jsonl")
        return

    warning_counts: dict[str, int] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            text = line.strip()
            if not text:
                continue
            try:
                event = json.loads(text)
            except json.JSONDecodeError:
                add(findings, "ERROR", f"os_events.jsonl line {line_number} is not valid JSON")
                continue
            if not isinstance(event, dict):
                add(findings, "ERROR", f"os_events.jsonl line {line_number} is not a JSON object")
                continue
            event_type = event.get("event_type")
            if not event_type:
                add(findings, "WARN", f"os_events.jsonl line {line_number} lacks event_type")
            if event_type == "doctor_warning":
                key = first_line(event.get("message")) or f"line {line_number}"
                warning_counts[key] = warning_counts.get(key, 0) + 1

    for message, count in warning_counts.items():
        if count >= 2:
            add(
                findings,
                "ERROR" if profile in {"strict", "release"} else "WARN",
                f"doctor warning repeated {count} times without encoded prevention: {message}",
            )


def run_helper(command: list[str]) -> tuple[int, str]:
    result = subprocess.run(command, text=True, capture_output=True)
    output = "\n".join(part for part in (result.stdout.strip(), result.stderr.strip()) if part)
    return result.returncode, output


def check_helper_scripts(findings: list[Finding], profile: str) -> None:
    checks = [
        [sys.executable, "scripts/codex_artifact_registry.py", "check"],
        [sys.executable, "scripts/codex_thesis_radar.py"],
        [sys.executable, "scripts/codex_os_dream.py", "validate", "--latest", "--allow-missing"],
    ]
    for command in checks:
        returncode, output = run_helper(command)
        if returncode != 0:
            add(findings, "ERROR", f"{' '.join(command)} failed: {first_line(output)}")
        elif profile in {"strict", "release"} and "UNMAPPED" in output:
            add(findings, "WARN", f"{' '.join(command)} reported unmapped non-P0 work")


def minimal_workstream(**overrides: Any) -> dict[str, Any]:
    item: dict[str, Any] = {
        "workstream_id": "selftest",
        "title": "Self Test",
        "status": "active",
        "priority": "P0",
        "goal": "exercise doctor failure detection",
        "current_next_action": "run self-test",
        "depends_on": [],
        "active_codex_session": None,
        "chat_label_or_thread": None,
        "linear_issue": None,
        "linear_sync": "not_needed",
        "today_doc_anchor": None,
        "workstream_spec": None,
        "active_jobs": [],
        "artifacts": [],
        "evidence": ["synthetic self-test evidence"],
        "last_verified": "2026-05-27T00:00:00+00:00",
        "next_check": "2026-05-28T00:00:00+00:00",
        "stale_after": "2026-05-28T00:00:00+00:00",
        "handoff_summary": "synthetic self-test row",
    }
    item.update(overrides)
    return item


def run_self_tests(findings: list[Finding], now: datetime) -> None:
    base = {
        "status_vocabulary": [
            "active",
            "running",
            "waiting",
            "blocked",
            "review",
            "backlog",
            "done_pending_review",
            "archived",
        ],
        "daily_cockpit": {"active_working_point_deck": {"url": "https://example.invalid/deck"}},
    }
    cases = [
        (
            "stale live workstream",
            [minimal_workstream(stale_after="2020-01-01T00:00:00+00:00")],
            "live workstream is stale",
        ),
        (
            "done workstream with active jobs",
            [minimal_workstream(status="done_pending_review", active_jobs=["cluster 1"], stale_after=None, next_check=None)],
            "done/archived workstream still has active_jobs",
        ),
        (
            "missing dependency",
            [minimal_workstream(depends_on=["workstream:missing"])],
            "missing dependency workstream missing",
        ),
    ]
    for label, workstreams, expected in cases:
        case_findings: list[Finding] = []
        check_register({**base, "workstreams": workstreams}, case_findings, now)
        if not any(expected in finding.message for finding in case_findings):
            add(findings, "ERROR", f"doctor self-test failed to catch {label}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("register", nargs="?", default=str(DEFAULT_REGISTER))
    parser.add_argument("--profile", choices=["daily", "strict", "release"], default="daily")
    parser.add_argument("--now", help="ISO timestamp override for testing")
    parser.add_argument("--warnings-as-errors", action="store_true")
    parser.add_argument("--self-test", action="store_true", help="run synthetic failure drills")
    args = parser.parse_args()

    now = parse_when(args.now) if args.now else datetime.now(timezone.utc)
    if now is None:
        print("ERROR: --now must be an ISO timestamp", file=sys.stderr)
        return 2

    findings: list[Finding] = []
    check_required_files(findings, args.profile)
    check_policy_routing(findings)
    check_event_log(findings, args.profile)

    try:
        data = load_register(Path(args.register))
    except (RegisterError, OSError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    check_register(data, findings, now)
    check_helper_scripts(findings, args.profile)
    if args.self_test or args.profile == "release":
        run_self_tests(findings, now)

    errors = [f for f in findings if f.severity == "ERROR"]
    warnings = [f for f in findings if f.severity == "WARN"]
    for finding in findings:
        print(f"{finding.severity}: {finding.message}")

    warnings_are_errors = args.warnings_as_errors or args.profile in {"strict", "release"}
    if errors or (warnings_are_errors and warnings):
        print(f"FAIL: {len(errors)} error(s), {len(warnings)} warning(s)")
        return 1
    print(f"OK: OS doctor profile={args.profile} passed with {len(warnings)} warning(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
