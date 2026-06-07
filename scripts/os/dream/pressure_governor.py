#!/usr/bin/env python3
"""Classify recurring dream pressure before it reaches daily surfaces."""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any


PRESSURE_STATUSES = {
    "actionable",
    "handled",
    "cooling",
    "appendix_only",
    "blocked_for_waking",
    "resolved",
}
REGISTER_PRESSURE_KINDS = {"stale_state", "active_job_status", "active_job_evidence"}
APPENDIX_ONLY_KINDS = {
    "latent_context_nudge",
    "negative_memory_candidate",
    "context_resonance_auto_maintenance",
    "research_scout_subtask",
    "dual_pro_research_plan",
    "cleanup_storage_deferred_improvement",
}
HANDLED_PROTOCOLS = {
    "stale_state": "register_workstream_refresh_contract",
    "active_job_status": "register_workstream_refresh_contract",
    "active_job_evidence": "register_workstream_refresh_contract",
    "wip_overload": "morning_cockpit_compression_contract",
    "latent_context_nudge": "context_resonance_homeostasis_contract",
    "negative_memory_candidate": "context_resonance_homeostasis_contract",
    "context_resonance_auto_maintenance": "context_resonance_homeostasis_contract",
    "research_scout_subtask": "research_scout_review_contract",
    "dual_pro_research_plan": "research_scout_review_contract",
    "cleanup_storage_deferred_improvement": "cleanup_retention_contract",
}
ACTIONABLE_STATUSES = {"actionable", "blocked_for_waking"}


def first_text(value: object) -> str:
    if isinstance(value, list):
        for item in value:
            text = first_text(item)
            if text:
                return text
        return ""
    if isinstance(value, dict):
        for key in ("workstream_id", "title", "label", "target", "evidence", "status"):
            text = first_text(value.get(key))
            if text:
                return text
        return ""
    text = " ".join(str(value or "").split())
    return text


def parse_datetime(value: object) -> datetime | None:
    text = first_text(value)
    if not text:
        return None
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    try:
        parsed = datetime.fromisoformat(text)
    except ValueError:
        return None
    if parsed.tzinfo is None:
        return parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


def is_umbrella_workstream(item: dict[str, Any]) -> bool:
    children = item.get("child_workstreams")
    return isinstance(children, list) and bool(children)


def active_job_is_live(job: object) -> bool:
    if not isinstance(job, dict):
        return bool(job)
    monitor = first_text(job.get("monitor")).lower()
    text = " ".join(
        first_text(job.get(key)).lower()
        for key in ("status", "monitor", "purpose", "label", "stamp")
    )
    if "not complete" in text or "no final" in text or "no artifact" in text:
        return True
    inactive_tokens = (
        "ready",
        "removed",
        "cleaned",
        "cleared",
        "exit status 0",
        "complete",
        "completed",
        "archived",
    )
    if any(token in text for token in inactive_tokens):
        return False
    active_tokens = ("submitted", "running", "idle", "held", "queued", "waiting", "retry")
    if any(token in text for token in active_tokens):
        return True
    return bool(monitor and monitor not in {"none", "null"})


def effective_active_jobs(item: dict[str, Any]) -> list[Any]:
    jobs = item.get("active_jobs") if isinstance(item.get("active_jobs"), list) else []
    if is_umbrella_workstream(item):
        return []
    return [job for job in jobs if active_job_is_live(job)]


def classify_register_workstream(item: dict[str, Any], now: datetime) -> dict[str, Any]:
    evidence = item.get("evidence") if isinstance(item.get("evidence"), list) else []
    raw_jobs = item.get("active_jobs") if isinstance(item.get("active_jobs"), list) else []
    live_jobs = effective_active_jobs(item)
    status = item.get("status")
    stale_after = parse_datetime(item.get("stale_after"))
    next_check = parse_datetime(item.get("next_check"))
    last_verified = parse_datetime(item.get("last_verified"))

    missing = []
    for field in ("evidence", "last_verified", "next_check", "stale_after"):
        value = item.get(field)
        if value in (None, "", []):
            missing.append(field)
    if stale_after is None:
        missing.append("stale_after_parseable")
    if next_check is None:
        missing.append("next_check_parseable")
    if last_verified is None:
        missing.append("last_verified_parseable")

    classification = "current"
    if status in {"done_pending_review", "archived"}:
        classification = "archive_review"
    elif missing or not evidence:
        classification = "needs_evidence"
    elif stale_after and stale_after <= now:
        classification = "stale"
    elif status in {"waiting", "blocked", "review"}:
        classification = "waiting"
    elif live_jobs:
        classification = "current"

    recommended_action = {
        "current": "keep hot; refresh next_check only after real evidence changes",
        "waiting": "preserve waiting/blocker reason and set next_check/stale_after from exact evidence",
        "stale": "inspect evidence before claiming status; update last_verified, next_check, and stale_after after waking validation",
        "needs_evidence": "add exact evidence or demote/archive; do not let this remain a live ambiguous workstream",
        "archive_review": "confirm whether this can remain out of live dream/status pressure",
    }[classification]

    return {
        "workstream_id": item.get("workstream_id"),
        "title": item.get("title"),
        "status": status,
        "classification": classification,
        "last_verified": item.get("last_verified"),
        "next_check": item.get("next_check"),
        "stale_after": item.get("stale_after"),
        "evidence_count": len(evidence),
        "active_job_count": len(live_jobs),
        "raw_active_job_count": len(raw_jobs),
        "umbrella_parent": is_umbrella_workstream(item),
        "missing_or_unparseable": sorted(set(missing)),
        "recommended_action": recommended_action,
        "mutation_boundary": "read_only_report",
    }


def source_lane_from_run_id(value: object) -> str:
    text = first_text(value)
    if "-lane-" not in text:
        return ""
    return text.rsplit("-lane-", 1)[-1]


def normalize_pressure(item: dict[str, Any], source_lane: str | None = None) -> dict[str, str]:
    raw_kind = first_text(item.get("kind")) or "unknown"
    target = (
        first_text(item.get("target"))
        or first_text(item.get("workstream_id"))
        or first_text(item.get("title"))
        or first_text(item.get("evidence"))
        or "global"
    )
    kind = raw_kind
    evidence = first_text(item.get("evidence"))
    if raw_kind == "recurring_hotspot" and "|" in evidence:
        split_kind, split_target = evidence.split("|", 1)
        kind = first_text(split_kind) or kind
        target = first_text(split_target) or target
    lane = first_text(source_lane) or first_text(item.get("source_lane")) or source_lane_from_run_id(item.get("latest_run_id"))
    return {
        "pressure_kind": kind,
        "pressure_target": target,
        "source_lane": lane or "global",
        "pressure_key": f"{kind}|{target}|{lane or 'global'}",
    }


def classify_pressure(
    item: dict[str, Any],
    *,
    register_rows_by_id: dict[str, dict[str, Any]] | None = None,
    source_lane: str | None = None,
) -> dict[str, Any]:
    normalized = normalize_pressure(item, source_lane=source_lane)
    kind = normalized["pressure_kind"]
    target = normalized["pressure_target"]
    row = (register_rows_by_id or {}).get(target)
    handler = HANDLED_PROTOCOLS.get(kind, "")
    status = "actionable"
    visibility = "daily"
    reason = "unhandled pressure remains waking-visible"

    if kind in REGISTER_PRESSURE_KINDS:
        handler = handler or "register_workstream_refresh_contract"
        if row and row.get("umbrella_parent"):
            status = "handled"
            visibility = "appendix"
            reason = "umbrella parent pressure is delegated to focused child workstreams"
        elif row and row.get("classification") in {"current", "waiting", "archive_review"}:
            status = "handled"
            visibility = "appendix"
            reason = f"register protocol classifies {target} as {row.get('classification')}"
        elif row and row.get("classification") in {"stale", "needs_evidence"}:
            status = "actionable"
            visibility = "daily"
            reason = f"register protocol classifies {target} as {row.get('classification')}"
        elif not row and target == "global":
            status = "cooling"
            visibility = "appendix"
            reason = "global register pressure lacks a concrete workstream target"
    elif kind in APPENDIX_ONLY_KINDS:
        status = "appendix_only"
        visibility = "appendix"
        reason = "proposal-only maintenance hint; keep raw artifact but suppress from daily pressure"
    elif kind == "wip_overload" and handler:
        status = "handled"
        visibility = "appendix"
        reason = "cockpit compression contract owns WIP pressure"
    else:
        score = item.get("score")
        try:
            score_int = int(score)
        except (TypeError, ValueError):
            score_int = 0
        if score_int and score_int < 75:
            status = "cooling"
            visibility = "appendix"
            reason = "below daily-action threshold"

    return {
        **normalized,
        "pressure_status": status,
        "handler": handler,
        "handled_by": handler if status == "handled" else "",
        "daily_visibility": visibility,
        "suppression_reason": reason if visibility != "daily" else "",
        "unhandled_actionable": status in ACTIONABLE_STATUSES and visibility == "daily",
    }


def decorate_pressure_item(
    item: dict[str, Any],
    *,
    register_rows_by_id: dict[str, dict[str, Any]] | None = None,
    source_lane: str | None = None,
) -> dict[str, Any]:
    return {
        **item,
        **classify_pressure(item, register_rows_by_id=register_rows_by_id, source_lane=source_lane),
    }


def pressure_summary(items: list[dict[str, Any]]) -> dict[str, Any]:
    status_counts = {status: 0 for status in sorted(PRESSURE_STATUSES)}
    visibility_counts = {"daily": 0, "appendix": 0}
    for item in items:
        status = first_text(item.get("pressure_status")) or "actionable"
        visibility = first_text(item.get("daily_visibility")) or "daily"
        status_counts[status] = status_counts.get(status, 0) + 1
        visibility_counts[visibility] = visibility_counts.get(visibility, 0) + 1
    actionable = [
        item
        for item in items
        if item.get("pressure_status") in ACTIONABLE_STATUSES and item.get("daily_visibility") == "daily"
    ]
    suppressed = [item for item in items if item.get("daily_visibility") != "daily"]
    return {
        "unhandled_actionable_count": len(actionable),
        "suppressed_count": len(suppressed),
        "daily_visible_count": visibility_counts.get("daily", 0),
        "appendix_only_count": visibility_counts.get("appendix", 0),
        "status_counts": status_counts,
        "visibility_counts": visibility_counts,
        "top_unhandled_keys": [first_text(item.get("pressure_key")) for item in actionable[:5]],
        "top_suppressed_keys": [first_text(item.get("pressure_key")) for item in suppressed[:5]],
    }
