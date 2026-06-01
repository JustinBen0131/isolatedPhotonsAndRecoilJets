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
from datetime import datetime, timezone
from pathlib import Path

from codex_work_register_common import DEFAULT_REGISTER, first_line, load_register, parse_when, sorted_workstreams


LIVE_STATUSES = {"active", "running", "waiting", "blocked", "review"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("register", nargs="?", default=str(DEFAULT_REGISTER))
    parser.add_argument("--now", help="ISO timestamp override for testing")
    args = parser.parse_args()

    now = parse_when(args.now) if args.now else datetime.now(timezone.utc)
    if now is None:
        raise SystemExit("--now must be an ISO timestamp")

    data = load_register(Path(args.register))
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
