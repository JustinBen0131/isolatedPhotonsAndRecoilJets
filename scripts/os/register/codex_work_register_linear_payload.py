#!/usr/bin/env python3
"""Emit Linear issue payloads from the Codex work register.

This script does not call Linear. It produces deterministic JSON that Codex can
use when the Linear connector is available.
"""

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
from pathlib import Path

from codex_work_register_common import DEFAULT_REGISTER, first_line, load_register, sorted_workstreams


STATUS_TO_LANE = {
    "active": "Active",
    "running": "Running",
    "waiting": "Waiting",
    "blocked": "Blocked",
    "review": "Review",
    "backlog": "Backlog",
    "done_pending_review": "Done Pending Archive",
    "archived": "Done Pending Archive",
}

CAMPAIGN_NAMES = {
    "ppg12_baseline_match": "P0 Campaign | Match pp Baseline Model to PPG12",
    "hp2026_photon_talk": "P0 Campaign | HP2026 Photon Talk",
    "final_model_choice": "P0 Campaign | Final Photon-ID Model Choice",
}

AREA_LABELS = {
    "approval": "Area: Approval",
    "hp26": "Area: Approval",
    "stitch": "Area: Stitching",
    "pp_exact": "Area: Stitching",
    "jet40": "Area: Stitching",
    "reweight": "Area: Stitching",
    "ppg12": "Area: Stitching",
    "shuhang": "Area: Stitching",
    "ml": "Area: ML",
    "model": "Area: ML",
    "bdt": "Area: ML",
    "feature": "Area: ML",
    "isolation": "Area: ML",
    "trigger": "Area: Trigger",
}


def linear_title_for(item: dict[str, object]) -> str:
    """Return a scan-first Linear title for a register workstream."""

    title = first_line(item.get("title")) or first_line(item.get("workstream_id"))
    workstream_id = str(item.get("workstream_id") or "").lower()
    status = str(item.get("status") or "")
    priority = str(item.get("priority") or "")
    text = f"{workstream_id} {title}".lower()

    if title.startswith(("P0 ", "P1 ", "RUN ", "WAIT ", "BKL ", "OS Baseline |")):
        return title
    campaign_id = str(item.get("campaign_id") or "")
    if campaign_id and workstream_id == campaign_id:
        return CAMPAIGN_NAMES.get(campaign_id, f"P0 Campaign | {title}")
    if "codex_os" in workstream_id or "operating system" in text:
        return f"OS Baseline | {title}"
    if status == "running" and ("stitch" in text or "jet40" in text or "pp_exact" in text):
        return f"RUN Stitching | {title}"
    if priority == "P0" and ("approval" in text or "hp26" in text or "talk" in text):
        return f"P0 Approval | {title}"
    if priority == "P0" and ("ml" in text or "model" in text or "bdt" in text):
        return f"P0 ML | {title}"
    if status == "waiting":
        return f"WAIT Validation | {title}"
    if priority == "P1":
        return f"P1 Closure | {title}"
    if status == "backlog":
        area = "ML Check" if any(token in text for token in ("ml", "bdt", "isolation", "feature")) else "Trigger" if "trigger" in text else "Backlog"
        return f"BKL {area} | {title}"
    return title


def labels_for(item: dict[str, object]) -> list[str]:
    """Return the baseline Linear label set implied by register state."""

    status = str(item.get("status") or "")
    priority = str(item.get("priority") or "")
    title = first_line(item.get("title"))
    workstream_id = str(item.get("workstream_id") or "")
    text = f"{workstream_id} {title}".lower()

    labels: list[str] = ["Workstream", "Codex OS"]
    if item.get("campaign_id"):
        labels.append("Campaign")
    if priority:
        labels.append(priority)

    if "operating system" in text or "codex_os" in text or status == "waiting":
        labels.append("Surface: Linear Detail")
    else:
        labels.append("Surface: Linear Detail" if status == "backlog" else "Surface: Today")

    for token, label in AREA_LABELS.items():
        if token in text and label not in labels:
            labels.append(label)

    if status == "running" or item.get("active_jobs"):
        labels.extend(["Running", "Ops: Active Jobs"])
    elif status in {"waiting", "blocked", "review", "backlog"}:
        labels.append(status.replace("_", " ").title())

    if "slide" in text or "talk" in text or "hp26" in text or "overlay" in text:
        labels.append("Slide-facing")
    if status == "review":
        labels.append("Needs: Justin Review")
    if "operating system" in text or "codex_os" in text:
        labels.append("Ops: Health Check")

    seen: set[str] = set()
    unique_labels: list[str] = []
    for label in labels:
        if label and label not in seen:
            unique_labels.append(label)
            seen.add(label)
    return unique_labels


def chat_map_for(item: dict[str, object]) -> list[str]:
    lines: list[str] = []
    campaign_id = first_line(item.get("campaign_id")) or "none"
    chat_label = first_line(item.get("chat_label_or_thread")) or "none"
    active_session = first_line(item.get("active_codex_session")) or "none"
    heartbeat = first_line(item.get("heartbeat_automation")) or "none"
    thread_id = first_line(item.get("codex_thread_id")) or "none"
    lines.extend(
        [
            f"- Campaign: `{campaign_id}`",
            f"- Chat label: `{chat_label}`",
            f"- Active Codex session: `{active_session}`",
            f"- Codex thread id: `{thread_id}`",
            f"- Heartbeat automation: `{heartbeat}`",
        ]
    )
    jobs = item.get("active_jobs") or []
    if isinstance(jobs, list) and jobs:
        for idx, job in enumerate(jobs[:4], start=1):
            if isinstance(job, dict):
                tag = first_line(job.get("stamp") or job.get("tag") or job.get("notes"))
                clusters = job.get("cluster_ids") or job.get("dag_cluster") or job.get("clusters")
                lines.append(f"- Active job {idx}: {tag or 'job'}; clusters={clusters}")
        if len(jobs) > 4:
            lines.append(f"- Additional active jobs: {len(jobs) - 4}")
    return lines


def body_for(item: dict[str, object]) -> str:
    lines = [
        "## At a glance",
        "",
        f"**Goal:** {first_line(item.get('goal'))}",
        "",
        f"**Current next action:** {first_line(item.get('current_next_action'))}",
        "",
        "## Codex chat map",
        "",
        *chat_map_for(item),
        "",
        f"**Depends on:** {first_line(item.get('depends_on')) or 'none'}",
        "",
        f"**Evidence:** {first_line(item.get('evidence')) or 'none'}",
        "",
        f"**Today's Plan:** {first_line(item.get('today_doc_anchor')) or 'pending'}",
        "",
        f"**Workstream spec:** {first_line(item.get('workstream_spec')) or 'none'}",
        "",
        f"**Active Codex session:** {first_line(item.get('active_codex_session')) or 'none'}",
        "",
        f"**Next check:** {first_line(item.get('next_check')) or 'none'}",
        "",
        f"**Stale after:** {first_line(item.get('stale_after')) or 'none'}",
        "",
        "## Latest handoff",
        "",
        first_line(item.get("handoff_summary")) or "none",
    ]
    jobs = item.get("active_jobs") or []
    if jobs:
        lines.extend(["", f"Active jobs registered: {len(jobs)}"])
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("register", nargs="?", default=str(DEFAULT_REGISTER))
    parser.add_argument("--include-archived", action="store_true")
    args = parser.parse_args()

    data = load_register(Path(args.register))
    payloads = []
    for item in sorted_workstreams(data):
        if item.get("status") == "archived" and not args.include_archived:
            continue
        payloads.append(
            {
                "workstream_id": item.get("workstream_id"),
                "title": linear_title_for(item),
                "status": item.get("status"),
                "lane": STATUS_TO_LANE.get(str(item.get("status")), "Backlog"),
                "priority": item.get("priority"),
                "labels": labels_for(item),
                "existing_linear_issue": item.get("linear_issue"),
                "linear_sync": item.get("linear_sync"),
                "body": body_for(item),
            }
        )

    print(json.dumps({"project": data.get("linear", {}).get("default_project"), "issues": payloads}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
