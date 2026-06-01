#!/usr/bin/env python3
"""Map active workstreams to the thesis narrative spine."""

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
import sys
from pathlib import Path
from typing import Any
import re

from codex_work_register_common import DEFAULT_REGISTER, first_line, load_register, sorted_workstreams


LIVE_STATUSES = {"active", "running", "waiting", "blocked", "review"}

CLAIM_KEYWORDS = {
    "physics_target": {
        "approval",
        "detector",
        "hp26",
        "photon",
        "prompt",
        "purity",
        "sideband",
        "talk",
    },
    "pp_baseline": {
        "basev3e",
        "overlay",
        "pp",
        "ppg12",
        "ppg18",
        "reweight",
        "shuhang",
        "stitch",
    },
    "embedded_background": {
        "background",
        "embedded",
        "inclusive",
        "jet12",
        "jet20",
        "jet30",
        "jet40",
        "stitch",
    },
    "ml_photon_id": {
        "ablation",
        "bdt",
        "closure",
        "fake",
        "feature",
        "isolation",
        "ml",
        "model",
        "wp80",
    },
    "auau_response": {
        "auau",
        "centrality",
        "response",
        "scaled",
        "trigger",
        "unfolding",
    },
    "final_physics": {
        "final",
        "nuclear",
        "observable",
        "xj",
        "xjgamma",
    },
    "os_infrastructure": {
        "codex",
        "doctor",
        "linear",
        "operating",
        "os",
        "register",
        "self-tuning",
        "symbiotic",
    },
}


def text_for(item: dict[str, Any]) -> str:
    parts = [
        item.get("workstream_id"),
        item.get("title"),
        item.get("goal"),
        item.get("current_next_action"),
        item.get("handoff_summary"),
    ]
    return " ".join(first_line(part).lower() for part in parts if part)


def claims_for(item: dict[str, Any]) -> list[str]:
    text = text_for(item)
    tokens = set(re.findall(r"[a-z0-9]+", text))
    claims = []
    for claim, keywords in CLAIM_KEYWORDS.items():
        if any((keyword in text if "-" in keyword or "_" in keyword else keyword in tokens) for keyword in keywords):
            claims.append(claim)
    return claims


def analyze(data: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for item in sorted_workstreams(data):
        if item.get("status") not in LIVE_STATUSES:
            continue
        claims = claims_for(item)
        rows.append(
            {
                "workstream_id": item.get("workstream_id"),
                "title": item.get("title"),
                "status": item.get("status"),
                "priority": item.get("priority"),
                "claims": claims,
                "mapped": bool(claims),
                "next": first_line(item.get("current_next_action")),
                "evidence": first_line(item.get("evidence")),
            }
        )
    return rows


def command_main(args: argparse.Namespace) -> int:
    data = load_register(Path(args.register))
    rows = analyze(data)
    unmapped = [row for row in rows if not row["mapped"]]
    unmapped_p0 = [row for row in unmapped if row.get("priority") == "P0"]

    if args.json:
        print(json.dumps({"rows": rows, "unmapped": unmapped}, indent=2, sort_keys=True))
    else:
        print("# Thesis Radar")
        print()
        for row in rows:
            marker = "OK" if row["mapped"] else "UNMAPPED"
            claims = ", ".join(row["claims"]) if row["claims"] else "none"
            print(
                f"{marker}: {row['workstream_id']} [{row['status']}, {row['priority']}] "
                f"claims={claims}"
            )
            if not row["mapped"]:
                print(f"  next={row['next']}")

    if unmapped_p0:
        for row in unmapped_p0:
            print(f"ERROR: active P0 does not map to thesis spine: {row['workstream_id']}", file=sys.stderr)
        return 1
    if args.strict and unmapped:
        for row in unmapped:
            print(f"ERROR: live workstream does not map to thesis spine: {row['workstream_id']}", file=sys.stderr)
        return 1
    print(f"OK: thesis radar mapped {len(rows) - len(unmapped)}/{len(rows)} live workstreams")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("register", nargs="?", default=str(DEFAULT_REGISTER))
    parser.add_argument("--strict", action="store_true", help="fail on any unmapped live workstream, not only P0")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    return command_main(args)


if __name__ == "__main__":
    raise SystemExit(main())
