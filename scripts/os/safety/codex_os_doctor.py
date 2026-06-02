#!/usr/bin/env python3
"""Run hardening checks for the ThesisAnalysis Codex operating system."""

from __future__ import annotations
# Keep purpose-folder helpers runnable when invoked directly.
import sys as _codex_sys
from pathlib import Path as _CodexPath
_CODEX_THIS_FILE = _CodexPath(__file__).resolve()
_CODEX_SCRIPTS_DIR = next((p for p in _CODEX_THIS_FILE.parents if p.name == "scripts"), _CODEX_THIS_FILE.parent)
_CODEX_OS_DIR = _CODEX_SCRIPTS_DIR / "os"
_CODEX_IMPORT_DIRS = [_CODEX_SCRIPTS_DIR]
if _CODEX_OS_DIR.exists():
    _CODEX_IMPORT_DIRS.extend(p for p in _CODEX_OS_DIR.iterdir() if p.is_dir())
for _CODEX_IMPORT_DIR in _CODEX_IMPORT_DIRS:
    _CODEX_IMPORT_DIR_STR = str(_CODEX_IMPORT_DIR)
    if _CODEX_IMPORT_DIR_STR not in _codex_sys.path:
        _codex_sys.path.append(_CODEX_IMPORT_DIR_STR)
del _CODEX_THIS_FILE, _CODEX_SCRIPTS_DIR, _CODEX_OS_DIR, _CODEX_IMPORT_DIRS, _CODEX_IMPORT_DIR, _CODEX_IMPORT_DIR_STR

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
    Path("agent_context/memory/README.md"),
    Path("agent_context/memory/CONTEXT_RESONANCE_INDEX.yaml"),
    Path("agent_context/memory/SCHEMA_REGISTRY.yaml"),
    Path("agent_context/memory/NEGATIVE_MEMORY_MAP.yaml"),
    Path("agent_context/templates/OS_POSTMORTEM_TEMPLATE.md"),
    Path("scripts/os/artifacts/codex_artifact_registry.py"),
    Path("scripts/os/context/codex_context_pack.py"),
    Path("scripts/os/context/codex_context_resonance.py"),
    Path("scripts/codex_os_guard.py"),
    Path("scripts/codex_os_dream.py"),
    Path("scripts/codex_os_nightly_heartbeat.py"),
    Path("scripts/os/research/codex_chatgpt_research_pack.py"),
    Path("scripts/codex_os_snapshot.py"),
    Path("scripts/os/context/codex_thesis_radar.py"),
]

STRICT_REQUIRED_FILES = [
    Path("agent_context/local/os_events.jsonl"),
]

DREAM_ROOT = Path("agent_context/local/dreams")
MEMORY_REGISTRY_FILES = {
    Path("agent_context/memory/CONTEXT_RESONANCE_INDEX.yaml"): "memory_records",
    Path("agent_context/memory/SCHEMA_REGISTRY.yaml"): "schemas",
    Path("agent_context/memory/NEGATIVE_MEMORY_MAP.yaml"): "negative_memories",
}
MEMORY_FORBIDDEN_PATTERNS = (
    "BEGIN TRANSCRIPT",
    "chatgpt.com/c/",
    "sphnxuser",
    "ssh.sdcc",
    "/sphenix/u/",
    "RJ_CODEX_THREAD_ID=",
    "password",
    "private key",
)
EXPECTED_DREAM_LANE_IDS = (
    "status_provenance",
    "architecture_cohesion",
    "context_resonance",
    "cleanup_storage",
    "path_contract",
    "research_scout",
    "science_scout",
    "presentation_artifacts",
)


@dataclass
class Finding:
    severity: str
    message: str


def text_of(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except OSError:
        return ""


def json_of(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return payload if isinstance(payload, dict) else {}


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
    if "scripts/codex_os_dream.py lane --lane-id" not in os_policy:
        add(findings, "WARN", "CODEX_OPERATING_SYSTEM.md does not mention the lane dream entrypoint")
    if "scripts/codex_os_nightly_heartbeat.py nightly" in os_policy:
        add(findings, "WARN", "CODEX_OPERATING_SYSTEM.md still mentions the retired nightly heartbeat entrypoint")
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
    if "codex_context_resonance.py" not in os_policy:
        add(findings, "WARN", "CODEX_OPERATING_SYSTEM.md does not mention context resonance resolver")


def source_pointer_path(item: dict[str, Any]) -> str:
    pointer = item.get("source_pointer")
    if isinstance(pointer, dict):
        return first_line(pointer.get("path"))
    return first_line(item.get("source"))


def check_memory_architecture(findings: list[Finding]) -> None:
    for path, list_key in MEMORY_REGISTRY_FILES.items():
        text = text_of(path)
        upper_text = text.upper()
        for pattern in MEMORY_FORBIDDEN_PATTERNS:
            if pattern.upper() in upper_text:
                add(findings, "ERROR", f"{path} contains forbidden memory-registry text: {pattern}")
        try:
            data = load_register(path)
        except (RegisterError, OSError, RuntimeError) as exc:
            add(findings, "ERROR", f"{path} is not parseable: {exc}")
            continue
        if not isinstance(data, dict):
            add(findings, "ERROR", f"{path} did not parse as a mapping")
            continue
        rows = data.get(list_key)
        if not isinstance(rows, list) or not rows:
            add(findings, "ERROR", f"{path} lacks non-empty {list_key}")
            continue
        for index, item in enumerate(rows, start=1):
            if not isinstance(item, dict):
                add(findings, "ERROR", f"{path} {list_key}[{index}] is not a mapping")
                continue
            memory_id = first_line(item.get("memory_id") or item.get("schema_id") or item.get("id"))
            if not memory_id:
                add(findings, "ERROR", f"{path} {list_key}[{index}] lacks memory_id/schema_id")
            relation_type = first_line(item.get("relation_type"))
            if relation_type and relation_type not in {
                "same_artifact",
                "same_method",
                "same_failure_mode",
                "analogy_only",
                "warning_only",
                "visual_style",
                "path_contract",
            }:
                add(findings, "ERROR", f"{path} {memory_id}: invalid relation_type {relation_type}")
            evidence_class = first_line(item.get("evidence_class"))
            if evidence_class not in {"real_observed", "human_approved", "derived", "synthetic"}:
                add(findings, "ERROR", f"{path} {memory_id}: invalid evidence_class {evidence_class!r}")
            retrieval_policy = first_line(item.get("retrieval_policy"))
            if retrieval_policy and retrieval_policy not in {
                "conscious_context",
                "latent_nudge",
                "suppress",
                "quarantine_candidate",
            }:
                add(findings, "ERROR", f"{path} {memory_id}: invalid retrieval_policy {retrieval_policy}")
            if evidence_class == "synthetic" and retrieval_policy == "conscious_context":
                add(findings, "ERROR", f"{path} {memory_id}: synthetic material cannot be conscious_context")
            if not first_line(item.get("required_waking_check")):
                add(findings, "ERROR", f"{path} {memory_id}: missing required_waking_check")
            source = source_pointer_path(item)
            if not source:
                add(findings, "ERROR", f"{path} {memory_id}: missing source pointer")
            elif source.startswith(("agent_context/", "scripts/", "codex_notes/", "macros/")) and not Path(source).exists():
                add(findings, "ERROR", f"{path} {memory_id}: local source pointer does not exist: {source}")


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


def latest_lane_signals() -> dict[str, tuple[Path | None, dict[str, Any] | None]]:
    rows: dict[str, tuple[Path | None, dict[str, Any] | None]] = {}
    if not DREAM_ROOT.exists():
        return rows
    candidates = [path for path in DREAM_ROOT.iterdir() if path.is_dir()]
    for lane_id in EXPECTED_DREAM_LANE_IDS:
        suffix = f"-lane-{lane_id}"
        lane_dirs = [path for path in candidates if path.name.endswith(suffix)]
        if not lane_dirs:
            rows[lane_id] = (None, None)
            continue
        latest = max(lane_dirs, key=lambda path: path.stat().st_mtime)
        signal_path = latest / "lane_signal.json"
        if not signal_path.exists():
            rows[lane_id] = (latest, None)
            continue
        try:
            signal = json.loads(signal_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            rows[lane_id] = (latest, None)
            continue
        rows[lane_id] = (latest, signal if isinstance(signal, dict) else None)
    return rows


def aggregate_lane_signals(bundle: dict[str, tuple[Path | None, dict[str, Any] | None]]) -> tuple[Path | None, dict[str, Any] | None]:
    valid_rows = [(path, signal) for path, signal in bundle.values() if path is not None and signal is not None]
    if not valid_rows:
        return None, None
    latest_path, latest_signal = max(
        valid_rows,
        key=lambda item: parse_when(item[1].get("generated_at")) or datetime.fromtimestamp(item[0].stat().st_mtime, tz=timezone.utc),
    )
    aggregate = dict(latest_signal)
    summary = dict(aggregate.get("summary") or {})
    summary["missing_lane_count"] = sum(1 for path, signal in bundle.values() if path is None or signal is None)
    aggregate["summary"] = summary
    aggregate["lane_statuses"] = {
        lane_id: {
            "path": path.as_posix() if path else None,
            "valid": signal is not None,
        }
        for lane_id, (path, signal) in bundle.items()
    }
    top_findings: list[dict[str, Any]] = []
    for lane_id in EXPECTED_DREAM_LANE_IDS:
        _, signal = bundle.get(lane_id, (None, None))
        if not signal:
            continue
        for item in signal.get("top_findings") or []:
            if isinstance(item, dict):
                top_findings.append(item)
    aggregate["top_findings"] = top_findings[:14]
    return latest_path, aggregate


def check_dream_heartbeat(findings: list[Finding], profile: str, now: datetime) -> None:
    bundle = latest_lane_signals()
    latest_dir, signal = aggregate_lane_signals(bundle)
    if latest_dir is None:
        return
    if signal is None:
        severity = "ERROR" if profile in {"strict", "release"} else "WARN"
        add(findings, severity, "latest dream lanes are missing valid lane signals")
        return

    missing_lanes = [lane_id for lane_id, (path, lane_signal) in bundle.items() if path is None or lane_signal is None]
    if missing_lanes:
        add(
            findings,
            "ERROR" if profile in {"strict", "release"} else "WARN",
            f"dream lane signals missing or invalid for: {', '.join(missing_lanes)}",
        )

    for lane_id, (path, lane_signal) in bundle.items():
        if path is None or lane_signal is None:
            continue
        generated_at = parse_when(lane_signal.get("generated_at"))
        if generated_at is not None and (now - generated_at).total_seconds() > 48 * 3600:
            add(findings, "WARN", f"dream lane heartbeat is older than 48h: {lane_id} at {path.name}")

    summary = signal.get("summary") if isinstance(signal.get("summary"), dict) else {}
    automation_drift = int(summary.get("automation_drift_count") or 0)
    recurring = int(summary.get("recurring_hotspot_count") or 0)
    handled_recurring = int(summary.get("handled_recurring_hotspot_count") or 0)
    cleanup_candidates = int(summary.get("cleanup_candidate_count") or 0)
    schema_candidates = int(summary.get("schema_promotion_candidate_count") or 0)
    approval_ready_count = int(summary.get("approval_ready_count") or 0)
    cohesion_score = int(signal.get("cohesion_score") or 0)
    debt = signal.get("maintenance_debt") if isinstance(signal.get("maintenance_debt"), dict) else {}
    debt_score = int(debt.get("score") or 0)
    budget_remaining = int(debt.get("budget_remaining") or 0)
    debt_status = str(debt.get("status") or "")
    heartbeat_status = str(signal.get("status") or "")
    maintenance = signal.get("evolutionary_maintenance") if isinstance(signal.get("evolutionary_maintenance"), dict) else {}
    pilot = maintenance.get("shadow_pilot") if isinstance(maintenance.get("shadow_pilot"), dict) else {}
    validation = signal.get("validation") if isinstance(signal.get("validation"), dict) else {}
    doctor_error_count = int(summary.get("doctor_error_count") or 0)
    for lane_id, (path, lane_signal) in bundle.items():
        if not lane_signal:
            continue
        changed_actions = lane_signal.get("changed_actions") if isinstance(lane_signal.get("changed_actions"), dict) else {}
        if changed_actions.get("external_mutations_performed") not in {False, None}:
            add(findings, "ERROR", f"dream lane changed_actions reports external mutation: {lane_id}")
        if changed_actions.get("science_mutations_performed") not in {False, None}:
            add(findings, "ERROR", f"dream lane changed_actions reports science mutation: {lane_id}")
        if changed_actions.get("repo_tracked_mutations_performed") not in {False, None}:
            add(findings, "ERROR", f"dream lane changed_actions reports repo-tracked mutation: {lane_id}")
        if lane_id == "context_resonance":
            resonance = lane_signal.get("context_resonance") if isinstance(lane_signal.get("context_resonance"), dict) else {}
            nudges = resonance.get("latent_context_nudges") if isinstance(resonance.get("latent_context_nudges"), list) else []
            negative = resonance.get("negative_memories") if isinstance(resonance.get("negative_memories"), list) else []
            suppressed = resonance.get("suppressed_context") if isinstance(resonance.get("suppressed_context"), list) else []
            if not resonance:
                add(findings, "WARN", "context_resonance lane signal lacks context_resonance payload")
            if not nudges:
                add(findings, "WARN", "context_resonance lane produced no latent_context_nudges")
            if not negative:
                add(findings, "WARN", "context_resonance lane produced no negative_memory candidates")
            if not suppressed:
                add(findings, "WARN", "context_resonance lane produced no suppressed_context entries")

    if automation_drift:
        add(
            findings,
            "ERROR" if profile in {"strict", "release"} else "WARN",
            f"dream heartbeat reports automation drift count={automation_drift}",
        )
    if recurring:
        add(findings, "WARN", f"dream heartbeat reports unhandled recurring maintenance hotspots count={recurring}")
    if schema_candidates:
        add(findings, "WARN", f"dream heartbeat has schema-promotion candidates count={schema_candidates}")
    if approval_ready_count:
        add(
            findings,
            "ERROR" if profile in {"strict", "release"} else "WARN",
            f"dream heartbeat emitted approval-ready maintenance items during shadow pilot count={approval_ready_count}",
        )
    if pilot and pilot.get("approval_ready_allowed") is not False:
        add(
            findings,
            "ERROR" if profile in {"strict", "release"} else "WARN",
            "dream shadow pilot does not have approval_ready_allowed=false",
        )
    if cleanup_candidates >= 3:
        add(findings, "WARN", f"dream heartbeat reports local cleanup candidates count={cleanup_candidates}")
    if cohesion_score < 65:
        add(findings, "WARN", f"dream heartbeat cohesion score is low: {cohesion_score}/100")
    if debt_score >= 60:
        add(findings, "WARN", f"dream heartbeat maintenance debt is elevated: {debt_score}/100")
    if budget_remaining < 40:
        add(findings, "WARN", f"dream heartbeat reliability budget remaining is low: {budget_remaining}/100")
    if debt_status == "frozen_growth":
        add(findings, "WARN", "dream heartbeat requests frozen-growth mode until debt is reduced")
    if heartbeat_status == "error":
        add(
            findings,
            "ERROR" if profile in {"strict", "release"} else "WARN",
            f"latest lane heartbeat aggregate reports error status: {latest_dir.name}",
        )
    if validation and validation.get("dream_validate_ok") is False:
        add(
            findings,
            "ERROR" if profile in {"strict", "release"} else "WARN",
            f"latest lane heartbeat aggregate captured dream validation failure: {latest_dir.name}",
        )
    if doctor_error_count:
        add(findings, "WARN", f"latest lane heartbeat aggregate captured doctor errors count={doctor_error_count}")


def check_helper_scripts(findings: list[Finding], profile: str) -> None:
    checks = [
        [sys.executable, "scripts/os/artifacts/codex_artifact_registry.py", "check"],
        [sys.executable, "scripts/os/context/codex_thesis_radar.py"],
        [sys.executable, "scripts/os/context/codex_context_resonance.py", "canary", "--json"],
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


def result_payload(
    findings: list[Finding],
    profile: str,
    now: datetime,
    warnings_are_errors: bool,
) -> dict[str, Any]:
    errors = [finding for finding in findings if finding.severity == "ERROR"]
    warnings = [finding for finding in findings if finding.severity == "WARN"]
    status = "fail" if errors or (warnings_are_errors and warnings) else "ok"
    return {
        "generated_at": now.isoformat(timespec="seconds"),
        "profile": profile,
        "status": status,
        "warnings_as_errors": warnings_are_errors,
        "counts": {
            "errors": len(errors),
            "warnings": len(warnings),
            "findings": len(findings),
        },
        "errors": [finding.message for finding in errors],
        "warnings": [finding.message for finding in warnings],
        "findings": [{"severity": finding.severity, "message": finding.message} for finding in findings],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("register", nargs="?", default=str(DEFAULT_REGISTER))
    parser.add_argument("--profile", choices=["daily", "strict", "release"], default="daily")
    parser.add_argument("--now", help="ISO timestamp override for testing")
    parser.add_argument("--warnings-as-errors", action="store_true")
    parser.add_argument("--self-test", action="store_true", help="run synthetic failure drills")
    parser.add_argument("--json", action="store_true", help="emit structured JSON findings")
    args = parser.parse_args()

    now = parse_when(args.now) if args.now else datetime.now(timezone.utc)
    if now is None:
        print("ERROR: --now must be an ISO timestamp", file=sys.stderr)
        return 2

    findings: list[Finding] = []
    check_required_files(findings, args.profile)
    check_policy_routing(findings)
    check_memory_architecture(findings)
    check_event_log(findings, args.profile)

    try:
        data = load_register(Path(args.register))
    except (RegisterError, OSError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    check_register(data, findings, now)
    check_helper_scripts(findings, args.profile)
    check_dream_heartbeat(findings, args.profile, now)
    if args.self_test or args.profile == "release":
        run_self_tests(findings, now)

    errors = [f for f in findings if f.severity == "ERROR"]
    warnings = [f for f in findings if f.severity == "WARN"]
    warnings_are_errors = args.warnings_as_errors or args.profile in {"strict", "release"}
    payload = result_payload(findings, args.profile, now, warnings_are_errors)

    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 1 if payload["status"] == "fail" else 0

    for finding in findings:
        print(f"{finding.severity}: {finding.message}")

    if errors or (warnings_are_errors and warnings):
        print(f"FAIL: {len(errors)} error(s), {len(warnings)} warning(s)")
        return 1
    print(f"OK: OS doctor profile={args.profile} passed with {len(warnings)} warning(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
