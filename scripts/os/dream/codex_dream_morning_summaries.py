#!/usr/bin/env python3
"""Aggregate eligible 03:30 dream lane summaries for Today's Plan."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, time, timedelta, timezone
from pathlib import Path
from typing import Any


DREAM_ROOT = Path("agent_context/local/dreams")
SYNTHETIC_HEADER = "SYNTHETIC DREAM OUTPUT - NOT USER APPROVAL - NOT REAL USER INTENT"
DREAM_SCHEDULED_SOURCE = "03:30_dream_automation"
EXPECTED_LANES = [
    "status_provenance",
    "architecture_cohesion",
    "context_resonance",
    "cleanup_storage",
    "path_contract",
    "research_scout",
    "science_scout",
    "presentation_artifacts",
]
STATUS_ORDER = {"changed": 0, "failed": 1, "deferred": 2, "no_safe_change": 3}


def parse_when(value: object) -> datetime | None:
    if not isinstance(value, str) or not value.strip():
        return None
    text = value.strip()
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    try:
        parsed = datetime.fromisoformat(text)
    except ValueError:
        return None
    if parsed.tzinfo is None:
        return parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


def now_utc() -> datetime:
    return datetime.now(timezone.utc)


def overnight_window(now: datetime) -> tuple[datetime, datetime]:
    local_now = now.astimezone()
    morning_cutoff = datetime.combine(local_now.date(), time(hour=12), tzinfo=local_now.tzinfo)
    if local_now < morning_cutoff:
        start_date = local_now.date() - timedelta(days=1)
    else:
        start_date = local_now.date()
    start = datetime.combine(start_date, time(hour=0), tzinfo=local_now.tzinfo)
    end = local_now + timedelta(hours=2)
    return start.astimezone(timezone.utc), end.astimezone(timezone.utc)


def read_json(path: Path) -> dict[str, Any] | None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return payload if isinstance(payload, dict) else None


def parse_txt_summary(path: Path, signal: dict[str, Any]) -> dict[str, Any] | None:
    try:
        text = path.read_text(encoding="utf-8")
    except OSError:
        return None
    if SYNTHETIC_HEADER not in text:
        return None
    fields: dict[str, str] = {}
    for line in text.splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        fields[key.strip().lower()] = value.strip()
    lane_id = fields.get("lane") or signal.get("lane_id")
    if not lane_id:
        return None
    scheduled = bool(signal.get("scheduled_morning_lane"))
    return {
        "synthetic": True,
        "synthetic_header": SYNTHETIC_HEADER,
        "lane_id": lane_id,
        "run_id": signal.get("run_id") or path.parent.name,
        "generated_at": signal.get("generated_at"),
        "status": fields.get("status") or "no_safe_change",
        "one_line_result": fields.get("one-line result") or "summary text fallback",
        "changed": fields.get("changed") or "none",
        "deferred": fields.get("deferred") or "none",
        "needs_justin": fields.get("needs justin", "").startswith("yes"),
        "needs_justin_reason": fields.get("needs justin", "no - none").split(" - ", 1)[-1],
        "evidence": fields.get("evidence") or str(path),
        "synthetic_boundary": fields.get("synthetic boundary") or "proposal-only, not user approval, not science evidence",
        "scheduled_morning_lane": scheduled,
        "scheduled_source": signal.get("scheduled_source"),
        "eligible_for_today_plan": bool(signal.get("eligible_for_today_plan")),
        "daily_plan_priority": 99,
    }


def load_summary(run_dir: Path) -> dict[str, Any] | None:
    payload = read_json(run_dir / "morning_lane_summary.json")
    if payload is not None:
        return payload
    signal = read_json(run_dir / "lane_signal.json") or {}
    if signal.get("eligible_for_today_plan") is not True:
        return None
    return parse_txt_summary(run_dir / "morning_lane_summary.txt", signal)


def within_window(payload: dict[str, Any], window: str, now: datetime) -> bool:
    if window == "all":
        return True
    generated = parse_when(payload.get("generated_at"))
    if generated is None:
        return False
    start, end = overnight_window(now)
    return start <= generated <= end


def eligible(payload: dict[str, Any], *, include_manual: bool) -> bool:
    if include_manual:
        return True
    return (
        payload.get("eligible_for_today_plan") is True
        and payload.get("scheduled_morning_lane") is True
        and payload.get("scheduled_source") == DREAM_SCHEDULED_SOURCE
    )


def select_latest_by_lane(
    dream_root: Path,
    *,
    window: str,
    include_manual: bool,
    now: datetime,
) -> tuple[dict[str, dict[str, Any]], dict[str, int]]:
    selected: dict[str, dict[str, Any]] = {}
    counts = {
        "scanned_count": 0,
        "skipped_validation_count": 0,
        "skipped_missing_summary_count": 0,
        "skipped_ineligible_count": 0,
        "skipped_window_count": 0,
    }
    if not dream_root.exists():
        return selected, counts
    for run_dir in sorted(path for path in dream_root.iterdir() if path.is_dir()):
        counts["scanned_count"] += 1
        if run_dir.name.startswith("validation-"):
            counts["skipped_validation_count"] += 1
            continue
        payload = load_summary(run_dir)
        if not payload:
            counts["skipped_missing_summary_count"] += 1
            continue
        if not eligible(payload, include_manual=include_manual):
            counts["skipped_ineligible_count"] += 1
            continue
        if not within_window(payload, window, now):
            counts["skipped_window_count"] += 1
            continue
        lane_id = str(payload.get("lane_id") or "")
        if lane_id not in EXPECTED_LANES:
            continue
        current = selected.get(lane_id)
        current_when = parse_when(current.get("generated_at")) if current else None
        payload_when = parse_when(payload.get("generated_at"))
        if current is None or (payload_when and current_when and payload_when > current_when) or (payload_when and current_when is None):
            payload = dict(payload)
            payload["run_dir"] = run_dir.as_posix()
            selected[lane_id] = payload
    return selected, counts


LANE_LABELS = {
    "status_provenance": "Status provenance",
    "architecture_cohesion": "Architecture cohesion",
    "context_resonance": "Context resonance",
    "cleanup_storage": "Cleanup/storage",
    "path_contract": "Path contract",
    "research_scout": "Research scout",
    "science_scout": "Science scout",
    "presentation_artifacts": "Presentation artifacts",
}

LANE_PURPOSES = {
    "status_provenance": "checks whether active work, jobs, and evidence look stale or need a read-only morning status pass",
    "architecture_cohesion": "checks whether the OS, doctor checks, automations, and register rules still agree with each other",
    "context_resonance": "checks the retrieval/memory feedback loop for drift, bad synthetic promotion, or missing context cues",
    "cleanup_storage": "checks local dream/research packs for safe retention, compaction, and cleanup proposals",
    "path_contract": "checks whether repo path, script, and helper changes are packaged for waking validation instead of silent self-application",
    "research_scout": "looks for research/process pressure points that may need a bounded waking research lane",
    "science_scout": "looks for science-analysis pressure points that may need waking validation or a new task",
    "presentation_artifacts": "checks slide/artifact hygiene and whether generated presentation assets need waking validation",
}


def lane_label(lane_id: object) -> str:
    lane = str(lane_id or "")
    return LANE_LABELS.get(lane, lane.replace("_", " ").title() or "Unknown lane")


def status_phrase(status: object) -> str:
    text = str(status or "no_safe_change")
    return {
        "changed": "Changed locally",
        "failed": "Failed",
        "deferred": "Deferred for waking Codex",
        "no_safe_change": "Quiet / no action",
    }.get(text, text.replace("_", " ").title())


def action_phrase(item: dict[str, Any]) -> str:
    status = str(item.get("status") or "no_safe_change")
    if status == "changed":
        return "No user action unless you want to inspect the local OS hygiene artifact."
    if status == "deferred":
        return "Treat this as a waking check, not as permission to mutate anything."
    if status == "failed":
        return "Inspect this before trusting the overnight dream bridge."
    return "No action needed."


def compact_text(value: object, *, limit: int = 240) -> str:
    text = " ".join(str(value or "").split())
    if len(text) <= limit:
        return text
    return text[: limit - 1].rstrip() + "..."


def markdown_bullets(payloads: list[dict[str, Any]], missing_lane_ids: list[str] | None = None) -> list[str]:
    lines: list[str] = []
    missing = missing_lane_ids or []
    if not payloads and not missing:
        return ["- No overnight dream update needs attention; internal checks stayed quiet."]
    lines.append("- Dream outputs are proposal-only: not user approval, not science evidence, not task completion.")
    visible = [item for item in payloads if item.get("status") in {"changed", "failed", "deferred", "no_safe_change"}]
    if not visible and not missing:
        lines.append("- No overnight dream update needs attention; internal checks stayed quiet.")
        return lines
    for item in visible:
        lane = str(item.get("lane_id") or "")
        result = compact_text(item.get("one_line_result") or "no summary")
        purpose = LANE_PURPOSES.get(lane, "checks one scheduled dream lane")
        lines.append(f"- **{lane_label(lane)} - {status_phrase(item.get('status'))}:** {purpose}.")
        lines.append(f"  - Result: {result}")
        if item.get("changed") and str(item.get("changed")) != "none":
            lines.append(f"  - Changed: {compact_text(item.get('changed'))}")
        if item.get("deferred") and str(item.get("deferred")) != "none":
            lines.append(f"  - Waking check: {compact_text(item.get('deferred'))}")
        evidence = item.get("evidence")
        if evidence:
            lines.append(f"  - Evidence: {compact_text(evidence)}")
        lines.append(f"  - Action: {action_phrase(item)}")
        if item.get("needs_justin"):
            lines.append(f"  - Needs Justin: {item.get('needs_justin_reason') or 'review requested'}")
    for lane in missing:
        lines.append(f"- **{lane_label(lane)} - No eligible 03:30 summary found:** {LANE_PURPOSES.get(lane, 'scheduled dream lane')}.")
        lines.append("  - Result: not shown in Today's Plan because the run did not produce an eligible scheduled-morning summary marker.")
        lines.append("  - Action: fix the dream-lane summary marker if this lane should be part of the morning digest.")
    return lines


def build_combined_summary(
    *,
    dream_root: Path = DREAM_ROOT,
    window: str = "overnight",
    include_manual: bool = False,
    now: datetime | None = None,
) -> dict[str, Any]:
    current = now or now_utc()
    latest, counts = select_latest_by_lane(dream_root, window=window, include_manual=include_manual, now=current)
    payloads = sorted(
        latest.values(),
        key=lambda item: (
            STATUS_ORDER.get(str(item.get("status")), 99),
            int(item.get("daily_plan_priority") or 99),
            str(item.get("lane_id")),
        ),
    )
    missing = [lane for lane in EXPECTED_LANES if lane not in latest]
    changed = [item for item in payloads if item.get("status") == "changed"]
    deferred = [item for item in payloads if item.get("status") == "deferred"]
    failed = [item for item in payloads if item.get("status") == "failed"]
    markdown = "\n".join(["## Overnight Dream Updates", *markdown_bullets(payloads, missing), ""])
    return {
        "generated_at": current.isoformat(timespec="seconds"),
        "window": window,
        "include_manual": include_manual,
        "scheduled_only": not include_manual,
        "summaries_found": len(payloads),
        "changed_count": len(changed),
        "deferred_count": len(deferred),
        "failed_count": len(failed),
        "missing_lane_ids": missing,
        "counts": counts,
        "top_updates": payloads[:6],
        "daily_doc_markdown": markdown,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dream-root", default=str(DREAM_ROOT))
    parser.add_argument("--window", choices=["overnight", "all"], default="overnight")
    parser.add_argument("--include-manual", action="store_true", help="debug: include manual or unapproved summaries")
    parser.add_argument("--scheduled-only", action="store_true", help="default; retained for explicit scripts")
    parser.add_argument("--now")
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--markdown", action="store_true")
    args = parser.parse_args()

    current = parse_when(args.now) if args.now else now_utc()
    summary = build_combined_summary(
        dream_root=Path(args.dream_root),
        window=args.window,
        include_manual=bool(args.include_manual),
        now=current or now_utc(),
    )
    if args.json:
        print(json.dumps(summary, indent=2, sort_keys=True))
    else:
        print(summary["daily_doc_markdown"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
