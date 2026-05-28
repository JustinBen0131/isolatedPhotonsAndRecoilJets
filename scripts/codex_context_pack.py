#!/usr/bin/env python3
"""Generate compact boot, daily, weekly, and monthly Codex context packs."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from codex_work_register_common import DEFAULT_REGISTER, first_line, load_register, parse_when, sorted_workstreams
from codex_thesis_radar import analyze as analyze_thesis_radar


LIVE_STATUSES = {"active", "running", "waiting", "blocked", "review"}
ACTIVE_STATUSES = {"active", "running"}
CLAIM_ORDER = [
    "physics_target",
    "pp_baseline",
    "embedded_background",
    "ml_photon_id",
    "auau_response",
    "final_physics",
    "os_infrastructure",
]


def extract_section(path: Path, heading: str, max_lines: int = 8) -> list[str]:
    if not path.exists():
        return []
    lines = path.read_text(encoding="utf-8").splitlines()
    start = None
    for index, line in enumerate(lines):
        if line.strip() == heading:
            start = index + 1
            break
    if start is None:
        return []
    out: list[str] = []
    for line in lines[start:]:
        if line.startswith("## "):
            break
        if line.strip():
            out.append(line)
        if len(out) >= max_lines:
            break
    return out


def working_point(data: dict[str, Any]) -> str:
    daily = data.get("daily_cockpit") if isinstance(data.get("daily_cockpit"), dict) else {}
    deck = daily.get("active_working_point_deck") if isinstance(daily, dict) else {}
    if not isinstance(deck, dict):
        return "none"
    label = first_line(deck.get("label")) or "WORKING POINT"
    url = first_line(deck.get("url")) or "missing"
    return f"{label}: {url}"


def stale_rows(workstreams: list[dict[str, Any]], now: datetime) -> list[dict[str, Any]]:
    rows = []
    for item in workstreams:
        if item.get("status") not in LIVE_STATUSES:
            continue
        stale_after = parse_when(item.get("stale_after"))
        if stale_after and stale_after <= now:
            rows.append(item)
    return rows


def print_workstream_lines(items: list[dict[str, Any]], limit: int | None = None) -> None:
    subset = items[:limit] if limit is not None else items
    for item in subset:
        print(
            f"- **{item.get('title')}** [{item.get('status')}, {item.get('priority')}]: "
            f"{first_line(item.get('current_next_action'))}"
        )
        evidence = first_line(item.get("evidence"))
        if evidence:
            print(f"  Evidence: {evidence}")
        jobs = item.get("active_jobs") or []
        if jobs:
            print(f"  Active jobs: {len(jobs)} registered")
    if not subset:
        print("- none")


def artifact_claim_counts() -> dict[str, int]:
    try:
        data = load_register(Path("agent_context/ARTIFACT_REGISTRY.yaml"))
    except Exception:
        return {}
    counts = {claim: 0 for claim in CLAIM_ORDER}
    for item in data.get("artifacts") or []:
        if isinstance(item, dict):
            claim = str(item.get("thesis_claim"))
            counts[claim] = counts.get(claim, 0) + 1
    return counts


def render_daily(data: dict[str, Any], workstreams: list[dict[str, Any]], now: datetime, max_top: int) -> None:
    active = [w for w in workstreams if w.get("status") in ACTIVE_STATUSES]
    waiting = [w for w in workstreams if w.get("status") in {"waiting", "blocked", "review"}]
    stale = stale_rows(workstreams, now)
    radar_rows = analyze_thesis_radar(data)

    print("## Daily Cockpit")
    print(f"- Generated: {now.isoformat(timespec='seconds')}")
    print(f"- {working_point(data)}")
    print("- Guard before risky actions: `python3 scripts/codex_os_guard.py preflight ...`")
    print()
    print("### Top Active")
    print_workstream_lines(active, max_top)
    print()
    print("### Waiting / Review / Blocked")
    print_workstream_lines(waiting)
    print()
    print("### Stale Or Needs Check")
    print_workstream_lines(stale)
    print()
    print("### Thesis Radar")
    for row in radar_rows:
        claims = ", ".join(row["claims"]) if row["claims"] else "UNMAPPED"
        print(f"- {row['workstream_id']}: {claims}")


def render_weekly(data: dict[str, Any], workstreams: list[dict[str, Any]], now: datetime) -> None:
    live = [w for w in workstreams if w.get("status") in LIVE_STATUSES]
    done = [w for w in workstreams if w.get("status") == "done_pending_review"]
    backlog = [w for w in workstreams if w.get("status") == "backlog"]
    stale = stale_rows(workstreams, now)
    print("## Weekly Review")
    print("- Promote, demote, close, or split workstreams based on thesis value.")
    print("- Convert repeated friction into a policy, script, guard, or Linear follow-up.")
    print()
    print("### Live Work")
    print_workstream_lines(live)
    print()
    print("### Stale Work")
    print_workstream_lines(stale)
    print()
    print("### Done Pending Review")
    print_workstream_lines(done)
    print()
    print("### Backlog Pressure")
    print_workstream_lines(backlog)


def render_monthly(data: dict[str, Any], workstreams: list[dict[str, Any]], now: datetime) -> None:
    counts = artifact_claim_counts()
    radar_rows = analyze_thesis_radar(data)
    mapped_claims = {claim for row in radar_rows for claim in row["claims"]}
    print("## Monthly Thesis Audit")
    print("- Audit claim layers, evidence gaps, stale outputs, and approval-readiness.")
    print()
    print("### Claim-Layer Evidence")
    for claim in CLAIM_ORDER:
        count = counts.get(claim, 0)
        live = "live-work" if claim in mapped_claims else "no-live-work"
        print(f"- {claim}: artifacts={count}, {live}")
    print()
    print("### Evidence Gaps")
    for claim in CLAIM_ORDER:
        if counts.get(claim, 0) == 0 and claim != "os_infrastructure":
            print(f"- {claim}: no registered canonical artifact yet")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("register", nargs="?", default=str(DEFAULT_REGISTER))
    parser.add_argument("--cadence", choices=["boot", "daily", "weekly", "monthly"], default="boot")
    parser.add_argument("--max-top", type=int, default=3)
    parser.add_argument("--now")
    args = parser.parse_args()

    now = parse_when(args.now) if args.now else datetime.now(timezone.utc)
    if now is None:
        raise SystemExit("--now must be an ISO timestamp")

    data = load_register(Path(args.register))
    workstreams = sorted_workstreams(data)
    north_star = extract_section(Path("agent_context/THESIS_NARRATIVE_MAP.md"), "## North Star")

    print("# ThesisAnalysis Codex Context Pack")
    print()
    print("## Boot Order")
    print("- Read `AGENTS.md`.")
    print("- Load route policies from `agent_context/policies/LOAD_MAP.yaml`.")
    print("- Treat `agent_context/CODEX_WORK_REGISTER.yaml` as canonical state.")
    print("- Use the safety guard before risky mutations.")
    print()
    print("## North Star")
    for line in north_star:
        print(line)
    print()

    if args.cadence in {"boot", "daily"}:
        render_daily(data, workstreams, now, args.max_top)
    elif args.cadence == "weekly":
        render_weekly(data, workstreams, now)
    elif args.cadence == "monthly":
        render_monthly(data, workstreams, now)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
