#!/usr/bin/env python3
"""Approval-gated preflight guard for risky ThesisAnalysis OS actions."""

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
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from codex_work_register_common import DEFAULT_REGISTER, RegisterError, first_line, load_register


EVENT_LOG = Path("agent_context/local/os_events.jsonl")
SNAPSHOT_DIR = Path("agent_context/local/os_snapshots")

RISKY_ACTIONS = {
    "condor_submit": "level4",
    "condor_job_control": "level4",
    "merge_or_production": "level4",
    "remote_edit": "level4",
    "transfer_upload": "level4",
    "transfer_download": "level4",
    "cleanup": "level4",
    "google_slides_mutation": "level4",
    "gmail_mark_read": "level3",
    "linear_closeout": "level3",
    "broad_plot_campaign": "level2",
}

DUPLICATE_FINGERPRINT_REQUIRED = {
    "condor_submit",
    "merge_or_production",
    "remote_edit",
    "transfer_upload",
    "transfer_download",
    "broad_plot_campaign",
}

LIVE_STATUSES = {"active", "running", "waiting", "blocked", "review"}

WEAK_APPROVALS = {
    "yes",
    "y",
    "ok",
    "approved",
    "do it",
    "go",
    "proceed",
    "sure",
}

TOKEN_STOPWORDS = {
    "the",
    "and",
    "for",
    "with",
    "from",
    "this",
    "that",
    "run",
    "job",
    "all",
    "new",
    "test",
    "plot",
    "file",
}


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def compact_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def append_event(event: dict[str, Any]) -> None:
    EVENT_LOG.parent.mkdir(parents=True, exist_ok=True)
    event.setdefault("timestamp", now_iso())
    with EVENT_LOG.open("a", encoding="utf-8") as handle:
        handle.write(compact_json(event) + "\n")


def load_events() -> list[dict[str, Any]]:
    if not EVENT_LOG.exists():
        return []
    events: list[dict[str, Any]] = []
    with EVENT_LOG.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            text = line.strip()
            if not text:
                continue
            try:
                event = json.loads(text)
            except json.JSONDecodeError:
                events.append({"event_type": "invalid_jsonl", "line": line_number, "raw": text})
                continue
            if isinstance(event, dict):
                events.append(event)
    return events


def normalize_fingerprint(text: str) -> str:
    words = re.findall(r"[a-z0-9_.:/+-]+", text.lower())
    return " ".join(words)


def fingerprint_digest(text: str) -> str:
    normalized = normalize_fingerprint(text)
    return hashlib.sha256(normalized.encode("utf-8")).hexdigest()[:20]


def fingerprint_tokens(text: str) -> set[str]:
    words = set(re.findall(r"[a-z0-9_+-]{4,}", text.lower()))
    return {word for word in words if word not in TOKEN_STOPWORDS}


def snapshot_exists(snapshot_id: str) -> bool:
    if not snapshot_id:
        return False
    candidate = Path(snapshot_id)
    if candidate.exists():
        return True
    if (SNAPSHOT_DIR / f"{snapshot_id}.tar.gz").exists():
        return True
    if (SNAPSHOT_DIR / f"{snapshot_id}.manifest.json").exists():
        return True
    return False


def find_workstream(data: dict[str, Any], workstream_id: str) -> dict[str, Any] | None:
    for item in data.get("workstreams") or []:
        if isinstance(item, dict) and item.get("workstream_id") == workstream_id:
            return item
    return None


def approval_is_weak(approval: str) -> bool:
    text = " ".join(approval.lower().split())
    return len(text) < 12 or text in WEAK_APPROVALS


def scope_is_vague(scope: str) -> bool:
    text = scope.lower()
    if len(scope.strip()) < 16:
        return True
    wildcard_terms = ("*", "everything", "whatever", "all files", "all jobs", "whole repo")
    return any(term in text for term in wildcard_terms)


def existing_fingerprint_matches(digest: str) -> list[dict[str, Any]]:
    matches: list[dict[str, Any]] = []
    for event in load_events():
        if event.get("fingerprint_digest") == digest and event.get("event_type") == "guard_preflight":
            matches.append(event)
    return matches


def register_token_overlap(data: dict[str, Any], fingerprint: str) -> list[str]:
    tokens = fingerprint_tokens(fingerprint)
    if len(tokens) < 4:
        return []
    haystack = compact_json(data).lower()
    hits = sorted(token for token in tokens if token in haystack)
    if len(hits) >= min(6, max(4, len(tokens) // 2)):
        return hits
    return []


def read_fingerprint(args: argparse.Namespace) -> str:
    chunks: list[str] = []
    if args.fingerprint:
        chunks.extend(args.fingerprint)
    if args.fingerprint_file:
        for path_text in args.fingerprint_file:
            path = Path(path_text)
            chunks.append(path.read_text(encoding="utf-8"))
    if args.field:
        chunks.extend(args.field)
    return "\n".join(chunks).strip()


def print_findings(errors: list[str], warnings: list[str]) -> None:
    for warning in warnings:
        print(f"WARN: {warning}")
    for error in errors:
        print(f"ERROR: {error}", file=sys.stderr)


def preflight(args: argparse.Namespace) -> int:
    errors: list[str] = []
    warnings: list[str] = []

    try:
        data = load_register(Path(args.register))
    except (RegisterError, OSError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    workstream = find_workstream(data, args.workstream)
    if not workstream:
        errors.append(f"unknown workstream: {args.workstream}")
    elif workstream.get("status") not in LIVE_STATUSES and args.action != "linear_closeout":
        warnings.append(
            f"{args.workstream}: status={workstream.get('status')} is not live; confirm this action is really current"
        )

    if args.action not in RISKY_ACTIONS:
        errors.append(f"unsupported action: {args.action}")

    if scope_is_vague(args.scope) and not args.broad_scope_reason:
        errors.append("scope is too vague; provide exact dataset/files/deck/job IDs or --broad-scope-reason")

    if approval_is_weak(args.approval or ""):
        errors.append("approval is missing or too weak; pass the exact Justin approval/scope quote")

    if not args.snapshot_id and not args.no_snapshot_needed:
        errors.append("snapshot is required before risky mutation; pass --snapshot-id or --no-snapshot-needed with a reason")
    if args.snapshot_id and not snapshot_exists(args.snapshot_id):
        errors.append(f"snapshot does not exist: {args.snapshot_id}")

    fingerprint = read_fingerprint(args)
    digest = ""
    if args.action in DUPLICATE_FINGERPRINT_REQUIRED:
        if not fingerprint:
            errors.append("duplicate-sensitive action requires --fingerprint, --fingerprint-file, or --field entries")
        else:
            digest = fingerprint_digest(fingerprint)
            matches = existing_fingerprint_matches(digest)
            if matches and not args.allow_duplicate:
                prior = matches[-1]
                errors.append(
                    "duplicate fingerprint already passed guard: "
                    f"digest={digest} workstream={prior.get('workstream_id')} "
                    f"time={prior.get('timestamp')}"
                )
            overlap = register_token_overlap(data, fingerprint)
            if overlap and not args.allow_duplicate:
                warnings.append(
                    "fingerprint has substantial overlap with existing register text; "
                    f"tokens={', '.join(overlap[:10])}. Confirm this is not a duplicate."
                )

    if args.command and any(token in args.command for token in ("rm -rf", "git reset --hard", "git checkout --")):
        errors.append("command includes a destructive pattern that requires a more specific dedicated approval")

    print_findings(errors, warnings)
    if errors:
        if args.json:
            print(compact_json({"ok": False, "errors": errors, "warnings": warnings}))
        return 1

    event = {
        "event_type": "guard_preflight",
        "action": args.action,
        "risk_level": RISKY_ACTIONS[args.action],
        "workstream_id": args.workstream,
        "scope": args.scope,
        "approval": args.approval,
        "snapshot_id": args.snapshot_id,
        "no_snapshot_needed": args.no_snapshot_needed,
        "fingerprint_digest": digest or None,
        "allow_duplicate": args.allow_duplicate,
        "broad_scope_reason": args.broad_scope_reason,
        "evidence": args.evidence or [],
        "artifacts": args.artifact or [],
        "command": args.command,
    }

    if not args.dry_run:
        append_event(event)

    print(
        "OK: guard preflight passed "
        f"action={args.action} workstream={args.workstream} "
        f"risk={RISKY_ACTIONS[args.action]} dry_run={args.dry_run}"
    )
    if digest:
        print(f"fingerprint_digest={digest}")
    if args.json:
        print(compact_json({"ok": True, "warnings": warnings, "event": event}))
    return 0


def fingerprint(args: argparse.Namespace) -> int:
    text = read_fingerprint(args)
    if not text:
        print("ERROR: provide --fingerprint, --fingerprint-file, or --field entries", file=sys.stderr)
        return 1
    digest = fingerprint_digest(text)
    matches = existing_fingerprint_matches(digest)
    payload = {
        "fingerprint_digest": digest,
        "tokens": sorted(fingerprint_tokens(text)),
        "event_matches": matches,
    }
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(f"fingerprint_digest={digest}")
        if matches:
            for event in matches:
                print(f"MATCH {event.get('timestamp')} {event.get('workstream_id')} {event.get('scope')}")
        else:
            print("OK: no matching guard preflight in local event log")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    pre = subparsers.add_parser("preflight", help="validate a risky action before performing it")
    pre.add_argument("--register", default=str(DEFAULT_REGISTER))
    pre.add_argument("--action", required=True, choices=sorted(RISKY_ACTIONS))
    pre.add_argument("--workstream", required=True)
    pre.add_argument("--scope", required=True)
    pre.add_argument("--approval", default="")
    pre.add_argument("--snapshot-id")
    pre.add_argument("--no-snapshot-needed")
    pre.add_argument("--fingerprint", action="append")
    pre.add_argument("--fingerprint-file", action="append")
    pre.add_argument("--field", action="append", help="fingerprint field such as dataset=... or tag=...")
    pre.add_argument("--allow-duplicate")
    pre.add_argument("--broad-scope-reason")
    pre.add_argument("--evidence", action="append")
    pre.add_argument("--artifact", action="append")
    pre.add_argument("--command")
    pre.add_argument("--dry-run", action="store_true")
    pre.add_argument("--json", action="store_true")
    pre.set_defaults(func=preflight)

    fp = subparsers.add_parser("fingerprint", help="compute and search a duplicate-run fingerprint")
    fp.add_argument("--fingerprint", action="append")
    fp.add_argument("--fingerprint-file", action="append")
    fp.add_argument("--field", action="append")
    fp.add_argument("--json", action="store_true")
    fp.set_defaults(func=fingerprint)

    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
