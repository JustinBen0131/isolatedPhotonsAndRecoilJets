#!/usr/bin/env python3
"""Create a Codex thread-heartbeat automation TOML when the app handler is missing.

This is a fallback for phone-originated or connector-degraded Codex sessions
where the `codex_app.automation_update` tool is visible but returns
`No handler registered`. Prefer the native app tool when it works.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import time
from pathlib import Path


VALID_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


def die(message: str) -> "NoReturn":  # type: ignore[name-defined]
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(2)


def toml_string(value: str) -> str:
    if '"""' in value:
        die('triple quotes are not supported inside TOML strings; use a prompt file without """')
    return '"""' + value.rstrip() + '"""'


def inline_string(value: str) -> str:
    if "\n" in value or "\r" in value:
        die("inline TOML string contains a newline")
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'


def read_prompt(args: argparse.Namespace) -> str:
    if args.prompt_file and args.prompt:
        die("use either --prompt-file or --prompt, not both")
    if args.prompt_file:
        return Path(args.prompt_file).read_text()
    if args.prompt:
        return args.prompt
    die("missing --prompt-file or --prompt")


def existing_created_at(path: Path) -> int | None:
    if not path.exists():
        return None
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("created_at = "):
            raw = line.split("=", 1)[1].strip()
            if raw.isdigit():
                return int(raw)
    return None


def render(args: argparse.Namespace, prompt: str, created_at: int, updated_at: int) -> str:
    return "\n".join(
        [
            "version = 1",
            f"id = {inline_string(args.id)}",
            'kind = "heartbeat"',
            f"name = {inline_string(args.name)}",
            f"prompt = {toml_string(prompt)}",
            f"status = {inline_string(args.status)}",
            f"rrule = {inline_string(args.rrule)}",
            f"target_thread_id = {inline_string(args.target_thread_id)}",
            f"created_at = {created_at}",
            f"updated_at = {updated_at}",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Write a Codex heartbeat automation TOML fallback."
    )
    parser.add_argument("--id", required=True, help="Automation id / directory name.")
    parser.add_argument("--name", required=True, help="Human-readable automation name.")
    parser.add_argument("--target-thread-id", required=True, help="Codex thread id.")
    parser.add_argument("--rrule", default="FREQ=MINUTELY;INTERVAL=10")
    parser.add_argument("--status", default="ACTIVE", choices=["ACTIVE", "PAUSED"])
    parser.add_argument("--prompt-file", help="Path containing the heartbeat prompt.")
    parser.add_argument("--prompt", help="Heartbeat prompt text.")
    parser.add_argument(
        "--automations-dir",
        default=str(Path(os.environ.get("CODEX_HOME", Path.home() / ".codex")) / "automations"),
    )
    parser.add_argument(
        "--update",
        action="store_true",
        help="Overwrite an existing automation.toml while preserving created_at when present.",
    )
    args = parser.parse_args()

    if not VALID_ID.match(args.id):
        die("automation id may contain only letters, digits, dot, dash, and underscore")
    prompt = read_prompt(args)
    if not prompt.strip():
        die("prompt is empty")

    automation_dir = Path(args.automations_dir).expanduser() / args.id
    toml_path = automation_dir / "automation.toml"
    if toml_path.exists() and not args.update:
        die(f"{toml_path} exists; pass --update to replace it")

    now_ms = int(time.time() * 1000)
    created_at = existing_created_at(toml_path) or now_ms
    automation_dir.mkdir(parents=True, exist_ok=True)
    tmp_path = toml_path.with_suffix(".toml.tmp")
    tmp_path.write_text(render(args, prompt, created_at, now_ms))
    tmp_path.replace(toml_path)
    print(toml_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
