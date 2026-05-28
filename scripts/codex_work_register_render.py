#!/usr/bin/env python3
"""Render the Codex work register into a Today's Plan cockpit summary."""

from __future__ import annotations

import argparse
from pathlib import Path

from codex_work_register_common import DEFAULT_REGISTER, first_line, load_register, sorted_workstreams


LIVE_STATUSES = {"active", "running", "waiting", "blocked", "review"}
ACTIVE_NOW_STATUSES = {"active", "running"}
TOP_ORDER = {
    "hp26_photon_id_talk": 0,
    "pp_exact_stitch_ppg12_baseline": 1,
    "jet40_embedded_inclusive_stitching": 2,
    "ml_final_model_ablation_map": 3,
    "et_eta_reweighting_explanation": 4,
}


def bullet(label: str, value: object) -> str:
    text = first_line(value)
    if not text:
        text = "none"
    return f"  - **{label}:** {text}"


def render_workstream(item: dict[str, object]) -> list[str]:
    title = item.get("title") or item.get("workstream_id")
    lines = [f"- **{title}** [{item.get('status')}, {item.get('priority')}]"]
    lines.append(bullet("Goal", item.get("goal")))
    lines.append(bullet("Next", item.get("current_next_action")))
    lines.append(bullet("Depends on", item.get("depends_on")))
    lines.append(bullet("Evidence", item.get("evidence")))
    lines.append(bullet("Linear", item.get("linear_issue") or item.get("linear_sync")))
    if item.get("workstream_spec"):
        lines.append(bullet("Spec", item.get("workstream_spec")))
    if item.get("active_jobs"):
        jobs = item.get("active_jobs") or []
        lines.append(f"  - **Active jobs:** {len(jobs)} registered")
    return lines


def linear_key(item: dict[str, object]) -> str:
    issue = first_line(item.get("linear_issue"))
    if not issue:
        return "Linear: pending"
    if "/issue/" in issue:
        slug = issue.rsplit("/issue/", 1)[-1]
        key = slug.split("/", 1)[0].upper()
        return f"Linear: {key}"
    return f"Linear: {issue}"


def daily_rank(item: dict[str, object]) -> tuple[int, int, str]:
    priority_order = {"P0": 0, "P1": 1, "P2": 2, "P3": 3}
    return (
        TOP_ORDER.get(str(item.get("workstream_id")), 99),
        priority_order.get(str(item.get("priority")), 99),
        str(item.get("title") or item.get("workstream_id")),
    )


def hp26_deck_lines(workstreams: list[dict[str, object]]) -> list[str]:
    hp = next((w for w in workstreams if w.get("workstream_id") == "hp26_photon_id_talk"), None)
    artifacts = hp.get("artifacts") if isinstance(hp, dict) else []
    if not isinstance(artifacts, list) or not artifacts:
        return [
            "- Current working-point deck: from Today's Plan or register context.",
            "- HP26 baseline deck: include when THE-5 is active.",
        ]
    lines: list[str] = []
    if len(artifacts) >= 1:
        lines.append(f"- HP truth-base deck: {artifacts[0]}")
    if len(artifacts) >= 2:
        lines.append(f"- HP reference/baseline deck: {artifacts[1]}")
    return lines


def working_point_line(data: dict[str, object]) -> str | None:
    daily = data.get("daily_cockpit")
    if not isinstance(daily, dict):
        return None
    deck = daily.get("active_working_point_deck")
    if not isinstance(deck, dict):
        return None
    label = first_line(deck.get("label")) or "WORKING POINT"
    url = first_line(deck.get("url"))
    if not url:
        return None
    return f"{label}: {url}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("register", nargs="?", default=str(DEFAULT_REGISTER))
    parser.add_argument("--max-top", type=int, default=3)
    args = parser.parse_args()

    data = load_register(Path(args.register))
    workstreams = sorted_workstreams(data)

    top = sorted([w for w in workstreams if w.get("status") in ACTIVE_NOW_STATUSES], key=daily_rank)
    active = sorted([w for w in workstreams if w.get("status") in ACTIVE_NOW_STATUSES], key=daily_rank)
    backlog = [w for w in workstreams if w.get("status") == "backlog"]
    live_for_checks = [w for w in workstreams if w.get("status") in {"active", "running", "waiting", "blocked", "review"}]

    print("# Today's Plan")
    print()
    print("Daily cockpit - Linear carries detail, repo register is truth")

    print()
    print("## Today's Work Contract")
    print("- Move at most 3 things forward today.")
    print("- Keep this doc short: daily execution only, not archive history.")
    print("- If Justin says \"add a task\", route it through repo register -> Linear -> Today's Plan if relevant.")
    print("- If a task needs detail, keep the detail in Linear and link it here.")
    print("- Tomorrow's plan should copy this structure, then update only changed facts.")
    active_working_point = working_point_line(data)
    if active_working_point:
        print(active_working_point)

    print()
    print("## Do First Today")
    for item in top[: args.max_top]:
        print(f"- **{item.get('title')}** - {first_line(item.get('current_next_action'))} {linear_key(item)}.")
    if not top:
        print("- No active registered work.")

    health = data.get("daily_cockpit", {}).get("morning_health_check", {})
    print()
    print("## Morning Health Check")
    if isinstance(health, dict) and health:
        print("Storage: report Justin-owned usage, not shared GPFS totals.")
        print("- `/sphenix/tg/tg01/bulk/jbennett`: used, quota, percent, file-count pressure.")
        print("- `/sphenix/u/patsfan753` and scratch: quota visibility and pressure.")
        print("Condor: read-only status only.")
        print("- Jobs: total, running, idle, held, and submit host.")
        print("- Tracked work: map live job groups to Linear/register workstreams.")
    else:
        print("- Morning health check is not configured.")

    print()
    print("## Slide Drivers")
    if active_working_point:
        print(f"- Default slide context: {active_working_point}")
    for line in hp26_deck_lines(workstreams):
        print(line)
    print("- Rule: if no deck link is specified, use the WORKING POINT deck above; accepted collaboration material first, new plots only when gaps remain.")

    print()
    print("## Active Now")
    for item in active:
        session = first_line(item.get("active_codex_session")) or "unclaimed"
        print(f"- **{str(item.get('status')).upper()} {item.get('title')}**. Next: {first_line(item.get('current_next_action'))} {linear_key(item)}. Session: {session}.")
    if not active:
        print("- No live workstreams.")

    print()
    print("## Waiting / Review")
    waiting_review = [w for w in workstreams if w.get("status") in {"waiting", "blocked", "review"}]
    for item in waiting_review:
        print(f"- **{str(item.get('status')).upper()} {item.get('title')}**: {first_line(item.get('current_next_action'))} {linear_key(item)}.")
    if not waiting_review:
        print("- No waiting, blocked, or review items.")

    print()
    print("## Active Jobs")
    active_jobs = [w for w in live_for_checks if w.get("active_jobs")]
    for item in active_jobs:
        jobs = item.get("active_jobs") or []
        print(f"- **{item.get('title')}**: {len(jobs)} registered job/watch item(s). Next: {first_line(item.get('current_next_action'))}")
    if not active_jobs:
        print("- No registered active jobs.")

    print()
    print("## Plot / Slide QA Queue")
    qa_items = []
    for item in live_for_checks:
        text = " ".join(
            first_line(value).lower()
            for value in (
                item.get("workstream_id"),
                item.get("title"),
                item.get("goal"),
                item.get("current_next_action"),
                item.get("artifacts"),
            )
        )
        if any(token in text for token in ("plot", "png", "slide", "deck", "overlay")):
            qa_items.append(item)
    for item in qa_items:
        print(f"- **{item.get('title')}** [{item.get('status')}]: {first_line(item.get('current_next_action'))} {linear_key(item)}.")
    if not qa_items:
        print("- No plot or slide QA items surfaced.")

    print()
    print("## Next Checks")
    for item in sorted(live_for_checks, key=lambda w: first_line(w.get("next_check"))):
        print(f"- **{item.get('title')}**: next_check={first_line(item.get('next_check'))}; stale_after={first_line(item.get('stale_after'))}.")
    if not live_for_checks:
        print("- No live next checks.")

    print()
    print("## Backlog - Do Not Let This Crowd Today")
    for item in backlog:
        print(f"- **{item.get('title')}**: {first_line(item.get('current_next_action'))} {linear_key(item)}.")
    if not backlog:
        print("- Backlog empty.")

    print()
    print("## Command Links")
    linear = data.get("linear", {}) if isinstance(data.get("linear"), dict) else {}
    if linear.get("project_url"):
        print(f"- Linear OS board: {linear.get('project_url')}")
    if linear.get("operating_rules_document"):
        print(f"- Operating rules: {linear.get('operating_rules_document')}")
    if linear.get("dashboard_doctrine_document"):
        print(f"- Linear doctrine: {linear.get('dashboard_doctrine_document')}")

    print()
    print("## Scratch / Inbox")
    print("Add raw notes here. Codex should route durable details into Linear, not let this page grow into an archive.")

    print()
    print("## Do Not Resurface As Active")
    print("GPFS outage status item; focused 21 GeV stitching comparison PNG; Blair fine reco-cluster ET artifacts.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
