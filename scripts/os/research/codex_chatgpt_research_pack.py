#!/usr/bin/env python3
"""Create and validate private ChatGPT UI research transcript packs."""

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
import re
import sys
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path("agent_context/local/chatgpt_research")
SENSITIVE_PATTERNS = {
    "ssh_host": re.compile(r"\b(?:ssh|sphnxuser)\S*\.?(?:sdcc|bnl)\S*", re.IGNORECASE),
    "sphenix_path": re.compile(r"/sphenix/[^\s)]+"),
    "credential_assignment": re.compile(
        r"\b(password|token|secret|api[_ -]?key|otp|credential)\s*[:=]\s*\S+",
        re.IGNORECASE,
    ),
    "private_key_block": re.compile(r"BEGIN [A-Z ]*PRIVATE KEY"),
    "condor_cluster": re.compile(r"\bcluster\s+\d{5,}\b", re.IGNORECASE),
}
REQUIRED_FILES = [
    "prompt.md",
    "single_message_prompt.txt",
    "session_strategy.md",
    "chatgpt_response.md",
    "codex_synthesis.md",
    "source_leads.md",
    "proposed_changes.md",
]
MAX_SINGLE_MESSAGE_CHARS = 3200
MODE_ORDER = ["instant", "thinking", "heavy", "pro"]


def now_id() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def sanitize_slug(text: str) -> str:
    slug = re.sub(r"[^a-zA-Z0-9]+", "-", text.lower()).strip("-")
    return slug[:48] or "chatgpt-research"


def header() -> str:
    return (
        "CHATGPT UI DELEGATED RESEARCH - EXTERNAL MODEL OUTPUT - "
        "NOT USER APPROVAL - VERIFY BEFORE USE"
    )


def scan_sensitive(text: str) -> list[str]:
    hits = []
    for label, pattern in SENSITIVE_PATTERNS.items():
        if pattern.search(text):
            hits.append(label)
    return hits


def build_prompt(topic: str, decision: str) -> str:
    return f"""# {header()}

## Sanitized Research Prompt

Research topic: {topic}

Decision this should inform: {decision}

Context, intentionally sanitized: We are designing a local, approval-gated
agentic operating system for a scientific analysis repository. It has task
state, artifact provenance, safety guards, doctor checks, and a private
synthetic "dream" layer that rehearses likely future failures without touching
live external systems.

Please provide:

1. Primary-source leads and named concepts from neuroscience of dreaming,
   memory consolidation, replay, threat simulation, abstraction, and emotional
   salience.
2. Machine-learning and agent-system analogs such as experience replay,
   self-play, world models, curriculum generation, reflection, adversarial
   evaluation, and synthetic data quality control.
3. Operations/safety analogs such as chaos engineering, postmortems,
   monitoring, incident review, progressive delivery, and guardrails.
4. Concrete design recommendations for a nightly synthetic rehearsal system
   that improves next-day readiness without overfitting to a fake user.
5. Failure modes, anti-patterns, and evaluation metrics.

Do not assume access to private project facts. Mark uncertain claims. Prefer
primary sources, named methods, and testable design recommendations.
"""


def build_single_message_prompt(topic: str, decision: str) -> str:
    prompt = (
        "You are doing external design research for a local, approval-gated scientific-analysis "
        "agent operating system. It has task state, artifact provenance, safety guards, doctor "
        "checks, and a private synthetic dream layer that rehearses likely future failures without "
        "touching live systems.\n\n"
        f"Research topic: {topic}\n"
        f"Decision to inform: {decision}\n\n"
        "Deliver a concrete answer with:\n"
        "1. Named primary-source leads from neuroscience of sleep, replay, memory consolidation, "
        "forgetting, salience, and schema formation.\n"
        "2. ML and agent analogs: experience replay, world models, reflection, adversarial eval, "
        "synthetic-data QA, memory pruning, anomaly detection.\n"
        "3. Ops analogs: chaos engineering, postmortems, SRE alerting, progressive delivery, "
        "runbooks, reliability scoring, hygiene debt control.\n"
        "4. A reference architecture for a nightly dream plus waking doctor sharing one heartbeat "
        "while preserving a hard reality boundary.\n"
        "5. Concrete modules, signals, recurrence logic, retention and forgetting policies, "
        "cleanup proposal systems, anti-patterns, and evaluation metrics.\n"
        "6. A shortlist of the highest-leverage ideas beyond typical agent-memory systems.\n\n"
        "Do not assume access to private project facts. Mark uncertain claims. Prefer primary "
        "sources, named methods, and testable recommendations."
    )
    return prompt


def recommend_mode(topic: str, decision: str) -> tuple[str, list[str]]:
    text = f"{topic} {decision}".lower()
    heavy_triggers = [
        "architecture",
        "research",
        "deep",
        "policy",
        "doctor",
        "dream",
        "thesis",
        "multi",
        "system",
        "workflow",
        "compare",
        "critique",
        "brainstorm",
        "maintain",
        "infrastructure",
    ]
    instant_triggers = [
        "rewrite",
        "wording",
        "summarize",
        "summary",
        "typo",
        "naming",
        "title",
        "quick",
    ]
    if any(trigger in text for trigger in heavy_triggers):
        return (
            "thinking",
            ["thinking", "heavy", "pro"],
        )
    if any(trigger in text for trigger in instant_triggers):
        return (
            "instant",
            ["instant", "thinking", "heavy"],
        )
    return (
        "thinking",
        ["thinking", "heavy", "pro"],
    )


def build_session_strategy(topic: str, decision: str) -> str:
    starting_mode, escalation = recommend_mode(topic, decision)
    escalation_text = " -> ".join(escalation)
    return f"""# {header()}

## Session Objective

Topic: {topic}

Decision to inform: {decision}

## Recommended First Mode

Start with `{starting_mode}`.

Rationale: begin with the lightest mode likely to preserve answer quality for
this objective. Escalate only if the first response is too shallow, misses key
constraints, or leaves the decision materially underdetermined.

## Escalation Ladder

{escalation_text}

## First-Message Rule

- Send `single_message_prompt.txt` as one complete message.
- Do not split the first prompt across multiple sends.
- If a partial send happens, abandon that thread for this objective and start a
  fresh chat.

## Follow-Up Rule

- One bounded follow-up at a time.
- Each follow-up should have one purpose only: clarify, challenge, deepen,
  compare, or extract.
- De-escalate after the hard question is answered.
"""


def command_init(args: argparse.Namespace) -> int:
    run_dir = ROOT / f"{now_id()}-{sanitize_slug(args.topic)}"
    run_dir.mkdir(parents=True, exist_ok=False)
    prompt = build_prompt(args.topic, args.decision)
    single_message_prompt = build_single_message_prompt(args.topic, args.decision)
    session_strategy = build_session_strategy(args.topic, args.decision)
    sensitive = scan_sensitive(prompt)
    if sensitive:
        print(f"ERROR: generated prompt matched sensitive patterns: {', '.join(sensitive)}", file=sys.stderr)
        return 1
    if len(single_message_prompt) > MAX_SINGLE_MESSAGE_CHARS:
        print(
            f"ERROR: single-message prompt is too long ({len(single_message_prompt)} chars > {MAX_SINGLE_MESSAGE_CHARS})",
            file=sys.stderr,
        )
        return 1
    sensitive = scan_sensitive(single_message_prompt)
    if sensitive:
        print(
            f"ERROR: single-message prompt matched sensitive patterns: {', '.join(sensitive)}",
            file=sys.stderr,
        )
        return 1

    files = {
        "prompt.md": prompt,
        "single_message_prompt.txt": single_message_prompt + "\n",
        "session_strategy.md": session_strategy,
        "chatgpt_response.md": f"# {header()}\n\nPaste or save the ChatGPT UI response here.\n",
        "followups.md": (
            f"# {header()}\n\n"
            "First send `single_message_prompt.txt` as one clean message.\n"
            "Do not split the first prompt across multiple sends.\n"
            "Only after the first response lands, add bounded follow-up prompts here one at a time.\n"
        ),
        "codex_synthesis.md": f"# {header()}\n\n## Useful Claims\n\n## Rejected Or Unverified Claims\n\n## Local Changes Proposed\n",
        "source_leads.md": f"# {header()}\n\n## Primary Source Leads\n\n## Secondary Source Leads\n",
        "proposed_changes.md": f"# {header()}\n\n## Policy Patches\n\n## Script Patches\n\n## Task Proposals\n",
    }
    for name, content in files.items():
        (run_dir / name).write_text(content, encoding="utf-8")
    print(run_dir)
    return 0


def command_validate(args: argparse.Namespace) -> int:
    root = Path(args.path) if args.path else latest_run()
    if root is None or not root.exists():
        print("ERROR: no ChatGPT research pack found", file=sys.stderr)
        return 1
    errors = []
    resolved_root = root.resolve()
    try:
        resolved_root.relative_to(ROOT.resolve())
    except ValueError:
        errors.append(f"research pack is outside {ROOT}: {root}")
    for name in REQUIRED_FILES:
        path = root / name
        if not path.exists():
            errors.append(f"missing required file: {path}")
            continue
        text = path.read_text(encoding="utf-8", errors="ignore")
        if name != "single_message_prompt.txt" and header() not in text:
            errors.append(f"missing external-output marker: {path}")
        hits = scan_sensitive(text)
        if hits:
            errors.append(f"sensitive pattern(s) in {path}: {', '.join(hits)}")
    prompt_text = (root / "single_message_prompt.txt").read_text(encoding="utf-8", errors="ignore")
    if len(prompt_text.strip()) > MAX_SINGLE_MESSAGE_CHARS:
        errors.append(
            f"single_message_prompt.txt exceeds {MAX_SINGLE_MESSAGE_CHARS} chars"
        )
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(f"OK: {root} validates as a private ChatGPT research pack")
    return 0


def latest_run() -> Path | None:
    if not ROOT.exists():
        return None
    runs = sorted([path for path in ROOT.iterdir() if path.is_dir()])
    return runs[-1] if runs else None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    init = subparsers.add_parser("init", help="create a sanitized research pack")
    init.add_argument("--topic", required=True)
    init.add_argument("--decision", required=True)
    init.set_defaults(func=command_init)

    validate = subparsers.add_parser("validate", help="validate a research pack")
    validate.add_argument("path", nargs="?")
    validate.set_defaults(func=command_validate)

    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
