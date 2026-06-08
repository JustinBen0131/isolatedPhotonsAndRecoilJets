#!/usr/bin/env python3
"""Maintain local slide iteration ledgers and Claude distillation packs."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ROOT = REPO_ROOT / "agent_context" / "local" / "slide_iterations"
CLAUDE_WORKER = REPO_ROOT / "scripts" / "os" / "delegation" / "claude_slide_worker.py"


def now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def slugify(value: str, limit: int = 64) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "-", value.strip().lower()).strip("-")
    return (slug or "slide-iteration")[:limit].strip("-")


def repo_path(value: str | Path) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = REPO_ROOT / path
    return path.resolve()


def display_path(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def ledger_path(run_dir: Path) -> Path:
    return run_dir / "iteration_ledger.jsonl"


def state_path(run_dir: Path) -> Path:
    return run_dir / "iteration_state.json"


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def read_events(run_dir: Path) -> list[dict[str, Any]]:
    ledger = ledger_path(run_dir)
    if not ledger.exists():
        return []
    events: list[dict[str, Any]] = []
    for line in ledger.read_text(encoding="utf-8").splitlines():
        if line.strip():
            events.append(json.loads(line))
    return events


def append_event(run_dir: Path, event: dict[str, Any]) -> None:
    event = {"recorded_at": now(), **event}
    with ledger_path(run_dir).open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(event, sort_keys=True) + "\n")


def cmd_init(args: argparse.Namespace) -> int:
    root = repo_path(args.out_root or DEFAULT_ROOT)
    root.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().astimezone().strftime("%Y%m%dT%H%M%S%z")
    run_dir = root / f"{stamp}-{slugify(args.slide_key)}"
    run_dir.mkdir(parents=True, exist_ok=False)
    state = {
        "schema_version": 1,
        "created_at": now(),
        "updated_at": now(),
        "slide_key": args.slide_key,
        "objective": args.objective,
        "deck": args.deck,
        "audience": args.audience,
        "status": "open",
    }
    write_json(state_path(run_dir), state)
    append_event(
        run_dir,
        {
            "event_type": "conception",
            "source": "codex",
            "iteration": 0,
            "summary": args.objective,
            "artifacts": [],
            "categories": ["initial_conception"],
        },
    )
    print(display_path(run_dir))
    print(f"ledger={display_path(ledger_path(run_dir))}")
    return 0


def cmd_record(args: argparse.Namespace) -> int:
    run_dir = repo_path(args.run_dir)
    if not state_path(run_dir).exists():
        print(f"missing iteration_state.json: {display_path(run_dir)}", file=sys.stderr)
        return 2
    state = load_json(state_path(run_dir))
    artifacts = [display_path(repo_path(value)) for value in args.artifact or []]
    append_event(
        run_dir,
        {
            "event_type": args.event_type,
            "source": args.source,
            "iteration": args.iteration,
            "summary": args.summary,
            "artifacts": artifacts,
            "categories": args.category or [],
        },
    )
    if args.event_type in {"accepted", "rejected", "superseded"}:
        state["status"] = args.event_type
    state["updated_at"] = now()
    write_json(state_path(run_dir), state)
    print(f"recorded={args.event_type}")
    print(f"ledger={display_path(ledger_path(run_dir))}")
    return 0


def build_summary(run_dir: Path) -> Path:
    state = load_json(state_path(run_dir))
    events = read_events(run_dir)
    attempts = [event for event in events if event.get("event_type") == "attempt"]
    feedback = [event for event in events if event.get("event_type") == "feedback"]
    accepted = [event for event in events if event.get("event_type") == "accepted"]
    claude = [event for event in events if event.get("source") == "claude"]
    latest_artifacts: list[str] = []
    for event in events:
        for artifact in event.get("artifacts", []):
            if artifact not in latest_artifacts:
                latest_artifacts.append(artifact)
    lines = [
        "# Slide Iteration Summary",
        "",
        f"- Slide key: `{state.get('slide_key', 'unknown')}`",
        f"- Objective: {state.get('objective', '')}",
        f"- Deck: {state.get('deck') or 'unspecified'}",
        f"- Audience: {state.get('audience') or 'unspecified'}",
        f"- Status: `{state.get('status', 'unknown')}`",
        f"- Events logged: {len(events)}",
        f"- Attempts logged: {len(attempts)}",
        f"- Feedback events: {len(feedback)}",
        f"- Claude events: {len(claude)}",
        f"- Accepted events: {len(accepted)}",
        "",
        "## Latest Artifacts",
        "",
    ]
    lines.extend([f"- `{artifact}`" for artifact in latest_artifacts] or ["- None recorded."])
    lines.extend(["", "## Event Timeline", ""])
    for event in events:
        lines.append(
            "- "
            + f"{event.get('recorded_at')} | iteration {event.get('iteration')} | "
            + f"{event.get('event_type')} | {event.get('source')}: {event.get('summary')}"
        )
        for category in event.get("categories", []):
            lines.append(f"  - category: `{category}`")
        for artifact in event.get("artifacts", []):
            lines.append(f"  - artifact: `{artifact}`")
    summary_path = run_dir / "iteration_summary.md"
    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return summary_path


def cmd_summarize(args: argparse.Namespace) -> int:
    summary = build_summary(repo_path(args.run_dir))
    print(display_path(summary))
    return 0


def cmd_pack_claude(args: argparse.Namespace) -> int:
    run_dir = repo_path(args.run_dir)
    summary = build_summary(run_dir)
    state = load_json(state_path(run_dir))
    objective = args.objective or (
        "Distill slide iteration history into prevention rules and a compact "
        "Codex injection brief for future slide generation."
    )
    cmd = [
        sys.executable,
        str(CLAUDE_WORKER),
        "pack",
        "--contract",
        "iteration_digest",
        "--objective",
        objective,
        "--include-file",
        str(ledger_path(run_dir)),
        "--include-file",
        str(summary),
        "--policy",
        "agent_context/policies/DELEGATION_KERNEL.md",
        "--policy",
        "agent_context/policies/SLIDES_WORKFLOW.md",
        "--policy",
        "agent_context/SLIDE_STYLE_MAP.md",
        "--notes",
        (
            f"Slide key: {state.get('slide_key')}. Deck: {state.get('deck')}. "
            "Distill iteration-count and correction evidence into prevention "
            "rules. Do not overgeneralize from weak evidence."
        ),
    ]
    completed = subprocess.run(
        cmd,
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    print(completed.stdout, end="")
    sys.stdout.flush()
    if completed.stderr:
        print(completed.stderr, file=sys.stderr, end="")
    if completed.returncode != 0:
        return completed.returncode
    pack_dir = next((line.strip() for line in completed.stdout.splitlines() if line.strip()), "")
    if args.invoke and pack_dir:
        invoke_cmd = [
            sys.executable,
            str(CLAUDE_WORKER),
            "invoke",
            "--pack-dir",
            pack_dir,
            "--max-budget-usd",
            args.max_budget_usd,
            "--timeout",
            str(args.timeout),
        ]
        invoked = subprocess.run(invoke_cmd, cwd=REPO_ROOT, text=True, check=False)
        report_path = repo_path(pack_dir) / "worker_report.md"
        if invoked.returncode == 0 and report_path.exists():
            events = read_events(run_dir)
            max_iteration = max((int(event.get("iteration", 0)) for event in events), default=0)
            append_event(
                run_dir,
                {
                    "event_type": "claude_report",
                    "source": "claude",
                    "iteration": max_iteration,
                    "summary": "Claude produced a Slide Iteration Learning Digest.",
                    "artifacts": [display_path(report_path)],
                    "categories": ["iteration_digest"],
                },
            )
            build_summary(run_dir)
        return invoked.returncode
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    init = subparsers.add_parser("init", help="Create a local slide iteration run.")
    init.add_argument("--slide-key", required=True)
    init.add_argument("--objective", required=True)
    init.add_argument("--deck")
    init.add_argument("--audience")
    init.add_argument("--out-root")
    init.set_defaults(func=cmd_init)

    record = subparsers.add_parser("record", help="Append one iteration event.")
    record.add_argument("--run-dir", required=True)
    record.add_argument(
        "--event-type",
        required=True,
        choices=[
            "conception",
            "attempt",
            "feedback",
            "accepted",
            "rejected",
            "superseded",
            "claude_report",
            "codex_audit",
        ],
    )
    record.add_argument("--source", required=True, choices=["justin", "codex", "claude", "tool"])
    record.add_argument("--iteration", type=int, required=True)
    record.add_argument("--summary", required=True)
    record.add_argument("--artifact", action="append")
    record.add_argument("--category", action="append")
    record.set_defaults(func=cmd_record)

    summarize = subparsers.add_parser("summarize", help="Write iteration_summary.md.")
    summarize.add_argument("--run-dir", required=True)
    summarize.set_defaults(func=cmd_summarize)

    pack = subparsers.add_parser("pack-claude", help="Build an iteration digest pack for Claude.")
    pack.add_argument("--run-dir", required=True)
    pack.add_argument("--objective")
    pack.add_argument("--invoke", action="store_true")
    pack.add_argument("--max-budget-usd", default="0.35")
    pack.add_argument("--timeout", type=int, default=240)
    pack.set_defaults(func=cmd_pack_claude)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
