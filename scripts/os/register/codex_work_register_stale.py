#!/usr/bin/env python3
"""Report stale live workstreams from the Codex work register."""

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
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from codex_work_register_common import DEFAULT_REGISTER, first_line, load_register, parse_when, sorted_workstreams
from pressure_governor import classify_register_workstream


LIVE_STATUSES = {"active", "running", "waiting", "blocked", "review"}


def print_protocol(rows: list[dict[str, Any]], as_json: bool) -> None:
    if as_json:
        print(
            json.dumps(
                {
                    "protocol_id": "register_workstream_refresh_contract",
                    "mutation_boundary": "read_only_report",
                    "rows": rows,
                    "summary": {
                        classification: sum(1 for row in rows if row["classification"] == classification)
                        for classification in ("current", "waiting", "stale", "needs_evidence", "archive_review")
                    },
                },
                indent=2,
                sort_keys=True,
            )
        )
        return
    print("PROTOCOL register_workstream_refresh_contract mutation_boundary=read_only_report")
    for row in rows:
        missing = ",".join(row["missing_or_unparseable"]) or "none"
        print(
            f"{row['classification'].upper()} {row['workstream_id']}: {row['title']} | "
            f"status={row['status']} | evidence={row['evidence_count']} | "
            f"active_jobs={row['active_job_count']} raw_active_jobs={row.get('raw_active_job_count', row['active_job_count'])} | "
            f"last_verified={row['last_verified']} | next_check={row['next_check']} | stale_after={row['stale_after']} | "
            f"umbrella_parent={str(row.get('umbrella_parent', False)).lower()} | missing={missing} | action={row['recommended_action']}"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("register", nargs="?", default=str(DEFAULT_REGISTER))
    parser.add_argument("--now", help="ISO timestamp override for testing")
    parser.add_argument(
        "--protocol",
        action="store_true",
        help="Run the read-only register_workstream_refresh_contract classifier instead of the legacy stale-only report.",
    )
    parser.add_argument("--json", action="store_true", help="Emit JSON for --protocol output")
    args = parser.parse_args()

    now = parse_when(args.now) if args.now else datetime.now(timezone.utc)
    if now is None:
        raise SystemExit("--now must be an ISO timestamp")

    data = load_register(Path(args.register))
    if args.protocol:
        rows = [
            classify_register_workstream(item, now)
            for item in sorted_workstreams(data)
            if item.get("status") in LIVE_STATUSES or item.get("status") in {"done_pending_review", "archived"}
        ]
        print_protocol(rows, args.json)
        return 0

    stale = []
    unparseable = []
    for item in sorted_workstreams(data):
        if item.get("status") not in LIVE_STATUSES:
            continue
        stale_after = parse_when(item.get("stale_after"))
        if stale_after is None:
            unparseable.append(item)
        elif stale_after <= now:
            stale.append(item)

    for item in stale:
        print(
            f"STALE {item.get('workstream_id')}: {item.get('title')} | "
            f"status={item.get('status')} | stale_after={item.get('stale_after')} | "
            f"next={first_line(item.get('current_next_action'))}"
        )

    for item in unparseable:
        print(
            f"CHECK_DATE {item.get('workstream_id')}: stale_after is not parseable: "
            f"{item.get('stale_after')!r}"
        )

    if not stale and not unparseable:
        print("OK: no stale live workstreams")
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
