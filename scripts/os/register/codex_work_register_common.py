#!/usr/bin/env python3
"""Shared helpers for the ThesisAnalysis Codex work register."""

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

import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


DEFAULT_REGISTER = Path("agent_context/CODEX_WORK_REGISTER.yaml")

REQUIRED_WORKSTREAM_FIELDS = (
    "workstream_id",
    "title",
    "status",
    "priority",
    "goal",
    "current_next_action",
    "depends_on",
    "active_codex_session",
    "chat_label_or_thread",
    "linear_issue",
    "linear_sync",
    "today_doc_anchor",
    "workstream_spec",
    "active_jobs",
    "artifacts",
    "evidence",
    "last_verified",
    "next_check",
    "stale_after",
    "handoff_summary",
)

DEFAULT_STATUSES = {
    "active",
    "running",
    "waiting",
    "blocked",
    "review",
    "backlog",
    "done_pending_review",
    "archived",
}

STATUS_ORDER = {
    "active": 0,
    "running": 1,
    "blocked": 2,
    "waiting": 3,
    "review": 4,
    "backlog": 5,
    "done_pending_review": 6,
    "archived": 7,
}


class RegisterError(RuntimeError):
    """Raised when the work register cannot be loaded or validated."""


def load_register(path: Path = DEFAULT_REGISTER) -> dict[str, Any]:
    """Load YAML through PyYAML when available, falling back to Ruby Psych."""

    path = Path(path)
    if not path.exists():
        raise RegisterError(f"missing register: {path}")

    try:
        import yaml  # type: ignore

        with path.open("r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle)
    except ModuleNotFoundError:
        cmd = [
            "ruby",
            "-ryaml",
            "-rjson",
            "-rdate",
            "-e",
            "puts JSON.generate(YAML.load_file(ARGV[0]))",
            str(path),
        ]
        result = subprocess.run(cmd, check=True, text=True, capture_output=True)
        data = json.loads(result.stdout)

    if not isinstance(data, dict):
        raise RegisterError(f"register root must be a mapping: {path}")
    return data


def validate_register(data: dict[str, Any]) -> list[str]:
    """Return validation errors. Empty list means the register is usable."""

    errors: list[str] = []
    statuses = set(data.get("status_vocabulary") or DEFAULT_STATUSES)
    if not statuses:
        errors.append("status_vocabulary is empty")

    workstreams = data.get("workstreams")
    if not isinstance(workstreams, list):
        errors.append("workstreams must be a list")
        return errors

    seen_ids: set[str] = set()
    for idx, item in enumerate(workstreams, start=1):
        prefix = f"workstreams[{idx}]"
        if not isinstance(item, dict):
            errors.append(f"{prefix} must be a mapping")
            continue

        missing = [field for field in REQUIRED_WORKSTREAM_FIELDS if field not in item]
        for field in missing:
            errors.append(f"{prefix} missing required field: {field}")

        workstream_id = item.get("workstream_id")
        if not isinstance(workstream_id, str) or not workstream_id.strip():
            errors.append(f"{prefix}.workstream_id must be a non-empty string")
        elif workstream_id in seen_ids:
            errors.append(f"{prefix}.workstream_id duplicates {workstream_id}")
        else:
            seen_ids.add(workstream_id)

        status = item.get("status")
        if status not in statuses:
            errors.append(f"{prefix} has invalid status: {status!r}")

        for list_field in ("depends_on", "active_jobs", "artifacts", "evidence"):
            if list_field in item and not isinstance(item.get(list_field), list):
                errors.append(f"{prefix}.{list_field} must be a list")

        next_action = item.get("current_next_action")
        if not isinstance(next_action, str) or not next_action.strip():
            errors.append(f"{prefix}.current_next_action must be non-empty")

        if item.get("status") in {"active", "running", "waiting", "blocked", "review"}:
            if not item.get("stale_after"):
                errors.append(f"{prefix}.stale_after is required for live statuses")
            if not item.get("next_check"):
                errors.append(f"{prefix}.next_check is required for live statuses")

    return errors


def parse_when(value: Any) -> datetime | None:
    """Parse a register timestamp/date. Returns None for loose prose dates."""

    if value is None:
        return None
    if isinstance(value, datetime):
        return value
    text = str(value).strip()
    if not text:
        return None
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    try:
        parsed = datetime.fromisoformat(text)
    except ValueError:
        try:
            parsed = datetime.fromisoformat(text + "T23:59:59")
        except ValueError:
            return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed


def sorted_workstreams(data: dict[str, Any]) -> list[dict[str, Any]]:
    """Sort workstreams by status and priority for readable projections."""

    priority_order = {"P0": 0, "P1": 1, "P2": 2, "P3": 3}

    def key(item: dict[str, Any]) -> tuple[int, int, str]:
        return (
            STATUS_ORDER.get(str(item.get("status")), 99),
            priority_order.get(str(item.get("priority")), 99),
            str(item.get("title") or item.get("workstream_id")),
        )

    workstreams = data.get("workstreams") or []
    return sorted([w for w in workstreams if isinstance(w, dict)], key=key)


def first_line(value: Any) -> str:
    """Compact a scalar/list value for one-line status output."""

    if value is None:
        return ""
    if isinstance(value, list):
        if not value:
            return ""
        return first_line(value[0])
    text = " ".join(str(value).split())
    return text
