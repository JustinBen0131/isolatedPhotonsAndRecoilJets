#!/usr/bin/env python3
"""Run the consolidated ThesisAnalysis nightly super-heartbeat."""

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
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from codex_work_register_common import DEFAULT_REGISTER, first_line, load_register, parse_when, sorted_workstreams


SYNTHETIC_HEADER = "SYNTHETIC DREAM OUTPUT - NOT USER APPROVAL - NOT REAL USER INTENT"
DREAM_ROOT = Path("agent_context/local/dreams")

LIVE_MORNING_STATUSES = {"active", "running", "review", "waiting", "blocked", "backlog"}
PRIORITY_ORDER = {"P0": 0, "P1": 1, "P2": 2, "P3": 3}
STATUS_ORDER = {"active": 0, "running": 1, "review": 2, "waiting": 3, "blocked": 4, "backlog": 5}


def now_utc() -> datetime:
    return datetime.now(timezone.utc)


def parse_now(value: str | None) -> datetime:
    if not value:
        return now_utc()
    parsed = parse_when(value)
    if parsed is None:
        raise ValueError("--now must be an ISO timestamp")
    return parsed.astimezone(timezone.utc)


def run_id_for(now: datetime) -> str:
    return f"{now.strftime('%Y%m%dT%H%M%SZ')}-nightly"


def safe_write(run_dir: Path, relative: str, content: str) -> Path:
    destination = (run_dir / relative).resolve()
    root = DREAM_ROOT.resolve()
    if root not in [destination, *destination.parents]:
        raise RuntimeError(f"refusing to write outside dream root: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(content, encoding="utf-8")
    return destination


def load_json(path: Path) -> dict[str, Any] | None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return payload if isinstance(payload, dict) else None


def run_command(command: list[str]) -> dict[str, Any]:
    result = subprocess.run(command, text=True, capture_output=True)
    stdout = result.stdout.strip()
    stderr = result.stderr.strip()
    output = "\n".join(part for part in (stdout, stderr) if part)
    return {
        "command": command,
        "returncode": result.returncode,
        "stdout": stdout,
        "stderr": stderr,
        "output": output,
    }


def nightly_commands(now: datetime, run_id: str, run_dir: Path) -> list[dict[str, Any]]:
    commands = [
        {
            "label": "dream_nightly",
            "command": [sys.executable, "scripts/codex_os_dream.py", "nightly", "--run-id", run_id],
        },
        {
            "label": "dream_validate",
            "command": [sys.executable, "scripts/codex_os_dream.py", "validate", str(run_dir)],
        },
        {
            "label": "doctor_daily",
            "command": [sys.executable, "scripts/codex_os_doctor.py", "--profile", "daily", "--json"],
        },
        {
            "label": "register_stale",
            "command": [sys.executable, "scripts/os/register/codex_work_register_stale.py", str(DEFAULT_REGISTER)],
        },
        {
            "label": "thesis_radar",
            "command": [sys.executable, "scripts/os/context/codex_thesis_radar.py"],
        },
        {
            "label": "artifact_registry_check",
            "command": [sys.executable, "scripts/os/artifacts/codex_artifact_registry.py", "check"],
        },
        {
            "label": "context_pack_daily",
            "command": [sys.executable, "scripts/os/context/codex_context_pack.py", "--cadence", "daily", "--max-top", "3"],
        },
    ]
    local_now = now.astimezone()
    if local_now.weekday() == 0:
        commands.append(
            {
                "label": "context_pack_weekly",
                "command": [sys.executable, "scripts/os/context/codex_context_pack.py", "--cadence", "weekly"],
            }
        )
    if local_now.day <= 3:
        commands.append(
            {
                "label": "context_pack_monthly",
                "command": [sys.executable, "scripts/os/context/codex_context_pack.py", "--cadence", "monthly"],
            }
        )
    return commands


def command_summary(result: dict[str, Any]) -> str:
    output = first_line(result.get("output"))
    if not output:
        return "OK"
    return output


def cleanup_candidates(trace: dict[str, Any] | None) -> list[dict[str, Any]]:
    if not isinstance(trace, dict):
        return []
    world = trace.get("world") if isinstance(trace.get("world"), dict) else {}
    local_state = world.get("local_state") if isinstance(world.get("local_state"), dict) else {}
    rows: list[dict[str, Any]] = []
    for label, summary in local_state.items():
        if not isinstance(summary, dict):
            continue
        for candidate in summary.get("cleanup_candidates") or []:
            if not isinstance(candidate, dict):
                continue
            rows.append(
                {
                    "surface": label,
                    "path": candidate.get("path"),
                    "age_days": candidate.get("age_days"),
                    "reason": candidate.get("reason"),
                }
            )
    return rows[:6]


def ranked_workstreams(register_data: dict[str, Any], signal: dict[str, Any] | None) -> list[dict[str, Any]]:
    risk_rank: dict[str, int] = {}
    for index, finding in enumerate((signal or {}).get("top_findings") or []):
        if isinstance(finding, dict):
            wid = first_line(finding.get("workstream_id"))
            if wid:
                risk_rank[wid] = index

    workstreams = []
    for item in sorted_workstreams(register_data):
        status = str(item.get("status"))
        if status not in LIVE_MORNING_STATUSES:
            continue
        workstreams.append(item)

    def sort_key(item: dict[str, Any]) -> tuple[int, int, int, str]:
        wid = first_line(item.get("workstream_id"))
        return (
            risk_rank.get(wid, 999),
            PRIORITY_ORDER.get(str(item.get("priority")), 9),
            STATUS_ORDER.get(str(item.get("status")), 9),
            first_line(item.get("title")) or wid,
        )

    return sorted(workstreams, key=sort_key)


def morning_appendix_payload(
    register_data: dict[str, Any],
    signal: dict[str, Any] | None,
    doctor: dict[str, Any],
    run_id: str,
) -> dict[str, Any]:
    ranked = ranked_workstreams(register_data, signal)
    do_first: list[dict[str, str]] = []
    if doctor.get("status") == "fail":
        do_first.append(
            {
                "label": "Resolve overnight heartbeat failure",
                "detail": f"Inspect `agent_context/local/dreams/{run_id}/nightly_heartbeat.md` before trusting the morning state.",
            }
        )
    for item in ranked:
        if len(do_first) >= 3:
            break
        do_first.append(
            {
                "label": f"{first_line(item.get('title'))} ({first_line(item.get('workstream_id'))})",
                "detail": first_line(item.get("current_next_action")) or "clarify the next action",
            }
        )

    waiting_blockers = []
    for item in ranked:
        status = str(item.get("status"))
        if status not in {"waiting", "blocked", "review"}:
            continue
        waiting_blockers.append(
            {
                "workstream_id": first_line(item.get("workstream_id")),
                "title": first_line(item.get("title")),
                "status": status,
                "decision": first_line(item.get("current_next_action")) or "clarify the pending decision",
            }
        )

    projection_candidates = []
    for item in ranked:
        if str(item.get("linear_sync")) != "synced" or item.get("today_doc_anchor") is None:
            projection_candidates.append(
                {
                    "workstream_id": first_line(item.get("workstream_id")),
                    "title": first_line(item.get("title")),
                    "status": str(item.get("status")),
                    "linear_sync": str(item.get("linear_sync")),
                }
            )

    active_deck = {}
    daily = register_data.get("daily_cockpit") if isinstance(register_data.get("daily_cockpit"), dict) else {}
    deck = daily.get("active_working_point_deck") if isinstance(daily.get("active_working_point_deck"), dict) else {}
    if deck:
        active_deck = {
            "label": first_line(deck.get("label")) or "WORKING POINT",
            "title": first_line(deck.get("title")) or "current working-point slides",
            "url": first_line(deck.get("url")),
        }

    return {
        "do_first": do_first[:3],
        "active_deck": active_deck,
        "waiting_blockers": waiting_blockers[:8],
        "projection_candidates": projection_candidates[:8],
    }


def render_morning_appendix(payload: dict[str, Any]) -> str:
    lines = [
        "# Morning Appendix",
        "",
        SYNTHETIC_HEADER,
        "",
        "## Do First Today",
    ]
    do_first = payload.get("do_first") or []
    if do_first:
        for index, item in enumerate(do_first, start=1):
            lines.append(f"{index}. {item['label']}: {item['detail']}")
    else:
        lines.append("1. No morning action proposals were generated.")

    lines.extend(["", "## Active Slide Deck"])
    active_deck = payload.get("active_deck") or {}
    if active_deck.get("url"):
        lines.append(f"- {active_deck.get('label')}: {active_deck.get('title')}")
        lines.append(f"- url: {active_deck.get('url')}")
    else:
        lines.append("- No active deck pointer is present in the register.")

    lines.extend(["", "## Waiting / Blocker Decisions"])
    waiting = payload.get("waiting_blockers") or []
    if waiting:
        for item in waiting:
            lines.append(
                f"- `{item['workstream_id']}` {item['status']}: {item['title']} -> {item['decision']}"
            )
    else:
        lines.append("- No waiting, blocked, or review workstreams need a morning decision.")

    lines.extend(["", "## Workstreams Likely To Need Linear / Today's Plan Projection"])
    candidates = payload.get("projection_candidates") or []
    if candidates:
        for item in candidates:
            lines.append(
                f"- `{item['workstream_id']}` {item['status']}: {item['title']} (linear_sync={item['linear_sync']})"
            )
    else:
        lines.append("- No projection candidates stand out from current register metadata.")

    lines.extend(
        [
            "",
            "## Boundary",
            "- Proposal only. Do not mutate Google Docs, Linear, SDCC, Condor, Gmail, Drive, Slides, or repo-tracked state from this appendix.",
        ]
    )
    return "\n".join(lines) + "\n"


def consecutive_clean_runs() -> int:
    if not DREAM_ROOT.exists():
        return 0
    candidates = [path for path in DREAM_ROOT.iterdir() if path.is_dir()]
    clean = 0
    for path in sorted(candidates, key=lambda item: item.stat().st_mtime, reverse=True):
        signal = load_json(path / "nightly_heartbeat_signal.json")
        if not signal:
            continue
        promotion = signal.get("promotion_gate") if isinstance(signal.get("promotion_gate"), dict) else {}
        current = promotion.get("current") if isinstance(promotion.get("current"), dict) else {}
        meets = (
            current.get("no_validation_failure") is True
            and current.get("no_automation_drift") is True
            and current.get("cohesion_score_ok") is True
            and current.get("maintenance_debt_ok") is True
            and current.get("no_doctor_errors") is True
        )
        if not meets:
            break
        clean += 1
    return clean


def promotion_gate(signal: dict[str, Any] | None, doctor: dict[str, Any], validation_ok: bool) -> dict[str, Any]:
    summary = signal.get("summary") if isinstance((signal or {}).get("summary"), dict) else {}
    debt = signal.get("maintenance_debt") if isinstance((signal or {}).get("maintenance_debt"), dict) else {}
    doctor_errors = doctor.get("errors") if isinstance(doctor.get("errors"), list) else []
    current = {
        "consecutive_clean_runs": consecutive_clean_runs(),
        "no_validation_failure": validation_ok,
        "no_automation_drift": int(summary.get("automation_drift_count") or 0) == 0,
        "cohesion_score_ok": int((signal or {}).get("cohesion_score") or 0) >= 70,
        "maintenance_debt_ok": int(debt.get("score") or 100) <= 45,
        "no_doctor_errors": len(doctor_errors) == 0,
        "explicit_user_approval": False,
    }
    eligible = (
        current["consecutive_clean_runs"] >= 3
        and current["no_validation_failure"]
        and current["no_automation_drift"]
        and current["cohesion_score_ok"]
        and current["maintenance_debt_ok"]
        and current["no_doctor_errors"]
        and current["explicit_user_approval"]
    )
    return {
        "enabled": False,
        "eligible": eligible,
        "criteria": {
            "consecutive_clean_runs": 3,
            "no_validation_failure": True,
            "no_automation_drift": True,
            "cohesion_score_min": 70,
            "maintenance_debt_max": 45,
            "no_doctor_errors": True,
            "explicit_user_approval": True,
        },
        "current": current,
    }


def heartbeat_status(
    results: list[dict[str, Any]],
    validation_ok: bool,
    doctor: dict[str, Any],
) -> str:
    if not validation_ok:
        return "error"
    if any(result["result"]["returncode"] != 0 for result in results):
        return "error"
    if doctor.get("status") == "fail":
        return "error"
    warnings = doctor.get("warnings") if isinstance(doctor.get("warnings"), list) else []
    if warnings:
        return "warn"
    return "ok"


def render_nightly_heartbeat(
    run_id: str,
    status: str,
    signal: dict[str, Any] | None,
    trace: dict[str, Any] | None,
    results: list[dict[str, Any]],
    doctor: dict[str, Any],
    appendix: dict[str, Any],
    promotion: dict[str, Any],
    validation_ok: bool,
    run_dir: Path,
) -> str:
    top_findings = (signal or {}).get("top_findings") or []
    summary = signal.get("summary") if isinstance((signal or {}).get("summary"), dict) else {}
    automation = signal.get("automation_state") if isinstance((signal or {}).get("automation_state"), dict) else {}
    debt = signal.get("maintenance_debt") if isinstance((signal or {}).get("maintenance_debt"), dict) else {}
    adaptation = (signal or {}).get("adaptation_cards") or []
    cleanup = cleanup_candidates(trace)
    doctor_errors = doctor.get("errors") if isinstance(doctor.get("errors"), list) else []
    doctor_warnings = doctor.get("warnings") if isinstance(doctor.get("warnings"), list) else []
    helper_failures = [item for item in results if item["label"] not in {"dream_nightly", "dream_validate", "doctor_daily"} and item["result"]["returncode"] != 0]

    lines = [
        f"# Nightly Heartbeat `{run_id}`",
        "",
        SYNTHETIC_HEADER,
        "",
        "## Overall Heartbeat Status",
        f"- status: `{status}`",
        f"- run directory: `{run_dir}`",
        f"- dream validation: {'passed' if validation_ok else 'failed'}",
        f"- doctor status: `{doctor.get('status', 'unknown')}`",
        f"- helper failures: {len(helper_failures)}",
        "",
        "## Top 3 Overnight Findings",
    ]
    if top_findings:
        for item in top_findings[:3]:
            target = item.get("workstream_id") or item.get("title") or item.get("evidence") or "global"
            lines.append(f"- `{item.get('kind')}` score={item.get('score')} target={target}")
    else:
        lines.append("- No overnight findings were captured.")

    lines.extend(["", "## Doctor Errors And Warnings"])
    if doctor_errors:
        for message in doctor_errors:
            lines.append(f"- ERROR: {message}")
    if doctor_warnings:
        for message in doctor_warnings:
            lines.append(f"- WARN: {message}")
    if not doctor_errors and not doctor_warnings:
        lines.append("- No doctor errors or warnings.")

    lines.extend(
        [
            "",
            "## Maintenance Debt And Cohesion",
            f"- cohesion score: {int((signal or {}).get('cohesion_score') or 0)}/100",
            f"- maintenance debt: {int(debt.get('score') or 0)}/100",
            f"- reliability budget remaining: {int(debt.get('budget_remaining') or 0)}/100",
            f"- debt status: `{debt.get('status') or 'unknown'}`",
            f"- adaptation cards: {len(adaptation)}",
        ]
    )

    lines.extend(
        [
            "",
            "## Automation Drift And Recurring Hotspots",
            f"- automation drift count: {int(summary.get('automation_drift_count') or 0)}",
            f"- recurring hotspot count: {int(summary.get('recurring_hotspot_count') or 0)}",
            f"- expected overnight ids: {', '.join(automation.get('expected_ids') or []) or 'none'}",
            f"- installed ids considered by dream: {', '.join(automation.get('installed_ids') or []) or 'none'}",
        ]
    )

    lines.extend(["", "## Adaptation Cards For Waking Review"])
    if adaptation:
        for item in adaptation[:3]:
            lines.append(
                f"- `{item.get('id')}` repeats={item.get('repeat_count')} status={item.get('promotion_status')} validator={item.get('validator_kind')}"
            )
    else:
        lines.append("- No adaptation cards crossed the recurrence threshold.")

    lines.extend(["", "## Cleanup Proposals Worth Review"])
    if cleanup:
        for item in cleanup[:3]:
            lines.append(
                f"- {item['surface']}: `{item['path']}` age={item['age_days']}d -> {item['reason']}"
            )
    else:
        lines.append("- No cleanup proposals stand out from the latest dream trace.")

    lines.extend(["", "## Top 3 Proposed Morning Actions"])
    for index, item in enumerate(appendix.get("do_first") or [], start=1):
        lines.append(f"{index}. {item['label']}: {item['detail']}")
    if not (appendix.get("do_first") or []):
        lines.append("1. No proposed morning actions were generated.")

    lines.extend(
        [
            "",
            "## Proposed Morning Appendix Summary",
            f"- do first count: {len(appendix.get('do_first') or [])}",
            f"- waiting/blocker decisions: {len(appendix.get('waiting_blockers') or [])}",
            f"- projection candidates: {len(appendix.get('projection_candidates') or [])}",
            f"- active deck present: {'yes' if (appendix.get('active_deck') or {}).get('url') else 'no'}",
        ]
    )

    lines.extend(
        [
            "",
            "## References",
            "- `dream_report.md`",
            "- `heartbeat_signal.json`",
            "- `maintenance_debt.md`",
            "- `adaptation_cards.md`",
            "- `morning_appendix.md`",
            "",
            "## No Live Mutations Performed",
            "- This heartbeat is proposal-only. It did not mutate SDCC, Condor, Gmail, Google Drive/Slides, Linear, or repo-tracked analysis state.",
            f"- Promotion gate enabled: {promotion.get('enabled')}",
        ]
    )
    return "\n".join(lines) + "\n"


def doctor_payload(result: dict[str, Any]) -> dict[str, Any]:
    parsed = None
    try:
        parsed = json.loads(result["stdout"]) if result["stdout"] else None
    except json.JSONDecodeError:
        parsed = None
    if isinstance(parsed, dict):
        return parsed
    return {
        "profile": "daily",
        "status": "fail" if result["returncode"] else "ok",
        "errors": [command_summary(result)] if result["returncode"] else [],
        "warnings": [],
    }


def orchestrate(args: argparse.Namespace) -> int:
    now = parse_now(args.now)
    run_id = args.run_id or run_id_for(now)
    run_dir = DREAM_ROOT / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for spec in nightly_commands(now, run_id, run_dir):
        result = run_command(spec["command"])
        results.append({"label": spec["label"], "result": result})

    dream_signal = load_json(run_dir / "heartbeat_signal.json")
    dream_trace = load_json(run_dir / "dream_trace.json")
    validation_result = next(item for item in results if item["label"] == "dream_validate")["result"]
    doctor_result = next(item for item in results if item["label"] == "doctor_daily")["result"]
    doctor = doctor_payload(doctor_result)
    validation_ok = validation_result["returncode"] == 0

    register_data = load_register(Path(args.register))
    appendix = morning_appendix_payload(register_data, dream_signal, doctor, run_id)
    appendix_path = safe_write(run_dir, "morning_appendix.md", render_morning_appendix(appendix))

    promotion = promotion_gate(dream_signal, doctor, validation_ok)
    status = heartbeat_status(results, validation_ok, doctor)
    nightly_signal = {
        "version": 1,
        "synthetic": True,
        "synthetic_header": SYNTHETIC_HEADER,
        "mutation_boundary": "proposal_only",
        "run_id": run_id,
        "mode": "nightly_heartbeat",
        "generated_at": (dream_signal or {}).get("generated_at") or now.isoformat(timespec="seconds"),
        "status": status,
        "cohesion_score": int((dream_signal or {}).get("cohesion_score") or 0),
        "maintenance_debt": (dream_signal or {}).get("maintenance_debt") or {},
        "summary": {
            **(((dream_signal or {}).get("summary") or {}) if isinstance((dream_signal or {}).get("summary"), dict) else {}),
            "doctor_error_count": len(doctor.get("errors") or []),
            "doctor_warning_count": len(doctor.get("warnings") or []),
            "helper_failure_count": sum(
                1
                for item in results
                if item["label"] not in {"dream_nightly", "dream_validate", "doctor_daily"}
                and item["result"]["returncode"] != 0
            ),
            "morning_action_count": len(appendix.get("do_first") or []),
        },
        "automation_state": (dream_signal or {}).get("automation_state") or {},
        "recurring_hotspots": (dream_signal or {}).get("recurring_hotspots") or [],
        "schema_promotion_candidates": (dream_signal or {}).get("schema_promotion_candidates") or [],
        "adaptation_cards": (dream_signal or {}).get("adaptation_cards") or [],
        "top_findings": ((dream_signal or {}).get("top_findings") or [])[:3],
        "morning_checks": (dream_signal or {}).get("morning_checks") or [],
        "references": {
            "dream_report": "dream_report.md",
            "dream_heartbeat": "heartbeat_signal.json",
            "maintenance_debt": "maintenance_debt.md",
            "adaptation_cards": "adaptation_cards.md",
            "morning_appendix": appendix_path.name,
            "nightly_heartbeat_report": "nightly_heartbeat.md",
        },
        "validation": {
            "dream_validate_ok": validation_ok,
            "dream_validate_output": validation_result.get("output") or "",
            "run_dir": run_dir.as_posix(),
        },
        "doctor": doctor,
        "helper_checks": [
            {
                "label": item["label"],
                "command": item["result"]["command"],
                "returncode": item["result"]["returncode"],
                "summary": command_summary(item["result"]),
            }
            for item in results
            if item["label"] not in {"doctor_daily"}
        ],
        "morning_appendix_summary": {
            "do_first_count": len(appendix.get("do_first") or []),
            "waiting_blocker_count": len(appendix.get("waiting_blockers") or []),
            "projection_candidate_count": len(appendix.get("projection_candidates") or []),
        },
        "promotion_gate": promotion,
    }

    report = render_nightly_heartbeat(
        run_id,
        status,
        dream_signal,
        dream_trace,
        results,
        doctor,
        appendix,
        promotion,
        validation_ok,
        run_dir,
    )
    report_path = safe_write(run_dir, "nightly_heartbeat.md", report)
    signal_path = safe_write(run_dir, "nightly_heartbeat_signal.json", json.dumps(nightly_signal, indent=2, sort_keys=True) + "\n")

    print(f"OK: nightly heartbeat wrote {run_dir}")
    print(f"report={report_path}")
    print(f"signal={signal_path}")
    print(f"appendix={appendix_path}")
    print(f"status={status}")
    return 0 if status != "error" else 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    nightly = subparsers.add_parser("nightly", help="run the consolidated nightly heartbeat")
    nightly.add_argument("--register", default=str(DEFAULT_REGISTER))
    nightly.add_argument("--run-id", help="override run id for testing")
    nightly.add_argument("--now", help="ISO timestamp override for testing")

    args = parser.parse_args()
    if args.command == "nightly":
        try:
            return orchestrate(args)
        except ValueError as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 2
    print(f"ERROR: unsupported command {args.command}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
