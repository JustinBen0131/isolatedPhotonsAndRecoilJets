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
import json
import re
import subprocess
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
    "clipboard_prompt.txt",
    "session_strategy.md",
    "mode_recommendation.json",
    "chatgpt_response.md",
    "codex_synthesis.md",
    "source_leads.md",
    "proposed_changes.md",
]
MAX_SINGLE_MESSAGE_CHARS = 3200
MODE_ORDER = [
    "instant",
    "thinking_light",
    "thinking_standard",
    "thinking_extended",
    "thinking_heavy",
    "pro_standard",
    "pro_extended",
    "deep_research",
]
USER_PASTE_HANDOFF_MODES = {"pro_standard", "pro_extended", "deep_research"}


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


def prompt_kind(topic: str, decision: str) -> str:
    text = f"{topic} {decision}".lower()
    if "chatgpt" in text and any(token in text for token in ("mode", "escalation", "routing", "submode")):
        return "chatgpt_mode_escalation"
    return "dream_architecture"


def build_prompt(topic: str, decision: str) -> str:
    if prompt_kind(topic, decision) == "chatgpt_mode_escalation":
        return f"""# {header()}

## Sanitized Research Prompt

Research topic: {topic}

Decision this should inform: {decision}

Context, intentionally sanitized: We are designing how a local coding agent
should use ChatGPT UI as a delegated research and critique worker. The agent
must route prompts to the fastest sufficient ChatGPT mode, prevent fragmented
prompt sends, avoid private-context leakage, and treat external-model output as
critique rather than truth.

Please provide:

1. A routing matrix for `instant`, `thinking-light`, `thinking-standard`,
   `thinking-extended`, `thinking-heavy`, `pro-standard`, `pro-extended`, and
   `deep-research`.
2. Escalation and de-escalation rules based on prompt length, number of
   constraints, citation/source burden, architecture depth, safety sensitivity,
   and prior-answer shallowness.
3. A wait-vs-human-pasteback policy for non-Pro vs Pro/extended Pro work.
4. A local scoring feature set a coding agent can implement to recommend mode
   and submode.
5. Anti-patterns and validation tests, especially fragmented prompting,
   overusing Pro, underusing reasoning, leaking private context, and treating
   ChatGPT output as project evidence.

Do not assume access to private project facts. Mark uncertain claims. Prefer
testable recommendations over generic advice.
"""
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
    if prompt_kind(topic, decision) == "chatgpt_mode_escalation":
        return (
            "Design a practical escalation policy for using ChatGPT UI modes as a delegated research "
            "and critique worker inside a local scientific-analysis coding agent. "
            f"Research topic: {topic}. Decision to inform: {decision}. "
            "Goal: route each prompt to the quickest mode likely to satisfy the accuracy and depth "
            "requirements, escalating only when measurable quality gates fail. Assume operating mode "
            "families and submodes are instant, thinking-light, thinking-standard, thinking-extended, "
            "thinking-heavy, pro-standard, pro-extended, and deep-research; treat labels as local UI "
            "contracts that can drift, not permanent product facts. Deliver: 1. a routing matrix from "
            "task complexity/risk to mode family and submode; 2. escalation and de-escalation rules "
            "based on prompt length, constraint count, citation/source-lead burden, architecture depth, "
            "safety sensitivity, and prior-answer shallowness; 3. when the coding agent should wait and "
            "collect vs submit once and ask the human to paste back Pro/extended-Pro output; 4. simple "
            "scoring features the coding agent can implement locally; 5. anti-patterns and validation "
            "tests, especially fragmented prompting, overusing Pro, underusing reasoning, leaking private "
            "context, and treating external-model output as truth. Do not assume private project facts. "
            "Mark uncertain claims. Prefer testable recommendations over generic advice."
        )
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


def build_clipboard_prompt(single_message_prompt: str) -> str:
    """Collapse layout line breaks so UI paste stays a single ChatGPT message."""
    paragraphs = re.split(r"\n\s*\n", single_message_prompt.strip())
    normalized = []
    for paragraph in paragraphs:
        normalized.append(re.sub(r"\s*\n\s*", " ", paragraph.strip()))
    return "\n\n".join(part for part in normalized if part)


def mode_family(mode: str) -> str:
    if mode.startswith("thinking_"):
        return "thinking"
    if mode.startswith("pro_"):
        return "pro"
    return mode


def recommend_mode(topic: str, decision: str) -> tuple[str, list[str]]:
    text = f"{topic} {decision}".lower()
    if "deep research" in text:
        return (
            "deep_research",
            ["deep_research"],
        )
    if any(trigger in text for trigger in ("extended pro", "maximum depth", "frontier", "exhaustive")):
        return (
            "pro_extended",
            ["pro_extended"],
        )
    if any(trigger in text for trigger in ("pro mode", "premium depth", "maximum rigor")):
        return (
            "pro_standard",
            ["pro_standard", "pro_extended"],
        )
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
        "escalation",
    ]
    extended_triggers = [
        "policy",
        "safety",
        "risk",
        "protocol",
        "validator",
        "runbook",
        "source",
        "citation",
        "matrix",
        "contract",
        "repo",
        "workflow",
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
    constraint_markers = sum(
        text.count(marker)
        for marker in (" and ", " or ", " but ", " while ", " unless ", " except ", "must", "should", "deliver")
    )
    if any(trigger in text for trigger in heavy_triggers) and constraint_markers >= 4:
        return (
            "thinking_heavy",
            ["thinking_heavy", "pro_standard", "pro_extended"],
        )
    if any(trigger in text for trigger in heavy_triggers + extended_triggers):
        return (
            "thinking_extended",
            ["thinking_extended", "thinking_heavy", "pro_standard"],
        )
    if any(trigger in text for trigger in instant_triggers):
        return (
            "instant",
            ["instant", "thinking_light", "thinking_standard"],
        )
    return (
        "thinking_standard",
        ["thinking_standard", "thinking_extended", "thinking_heavy"],
    )


def collection_policy_for_mode(mode: str) -> str:
    return "user_paste_handoff" if mode in USER_PASTE_HANDOFF_MODES else "codex_collects"


def score_mode_features(topic: str, decision: str) -> dict[str, int | bool]:
    text = f"{topic} {decision}".lower()
    constraint_markers = sum(
        text.count(marker)
        for marker in (" and ", " or ", " but ", " while ", " unless ", " except ", "must", "should", "deliver")
    )
    return {
        "prompt_chars_estimate": len(build_single_message_prompt(topic, decision)),
        "constraint_markers": constraint_markers,
        "source_burden": any(token in text for token in ("source", "citation", "literature", "research", "internet")),
        "architecture_depth": any(token in text for token in ("architecture", "system", "workflow", "os", "brain", "memory")),
        "safety_sensitivity": any(token in text for token in ("safety", "policy", "risk", "guard", "approval", "external")),
        "prior_answer_failure": any(token in text for token in ("shallow", "missed", "failed", "retry", "incomplete")),
        "pro_requested": any(token in text for token in ("pro mode", "extended pro", "premium depth", "maximum rigor")),
        "deep_research_requested": "deep research" in text,
        "quick_task": any(token in text for token in ("quick", "rewrite", "typo", "title", "wording")),
    }


def build_mode_recommendation_json(topic: str, decision: str) -> str:
    mode, escalation = recommend_mode(topic, decision)
    features = score_mode_features(topic, decision)
    payload = {
        "external_output_marker": header(),
        "recommended_mode": mode,
        "mode_family": mode_family(mode),
        "escalation_ladder": escalation,
        "response_collection_policy": collection_policy_for_mode(mode),
        "features": features,
        "rule": "Use the fastest mode that satisfies the evidence burden; escalate only on missed constraints, shallow reasoning, source burden, safety sensitivity, or prior failure.",
    }
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def render_collection_policy(mode: str, escalation: list[str]) -> str:
    policy = collection_policy_for_mode(mode)
    if policy == "user_paste_handoff":
        return (
            "Codex sends the one complete sanitized prompt, records the local pack and ChatGPT mode, "
            "then stops. Justin pastes the completed ChatGPT response back into Codex when it is done."
        )
    return (
        "Codex may send the one complete sanitized prompt, wait for the response, collect it, "
        "ask bounded follow-ups if needed, and synthesize the answer in the same Codex turn."
    )


def build_handoff_note(topic: str, decision: str) -> str:
    starting_mode, escalation = recommend_mode(topic, decision)
    policy = collection_policy_for_mode(starting_mode)
    return f"""# {header()}

## ChatGPT Response Collection Policy

- recommended_mode: `{starting_mode}`
- mode_family: `{mode_family(starting_mode)}`
- escalation_ladder: `{" -> ".join(escalation)}`
- response_collection_policy: `{policy}`

## Operating Rule

{render_collection_policy(starting_mode, escalation)}

## Pro / Deep Research Handoff Template

Use this message after submitting a Pro, extended Pro, or Deep Research prompt:

```text
I sent the sanitized prompt to ChatGPT in {starting_mode} mode.
Local research pack: <pack path>

This is a long-running external research job, so I am not going to spend Codex
tokens polling it. Paste the completed ChatGPT response back into this chat
when it finishes, and I will verify, synthesize, and turn only the useful parts
into local proposal artifacts.
```
"""


def build_thread_strategy(topic: str, decision: str) -> str:
    return f"""# {header()}

## Thread Context Decision

Default recommendation: `fresh`.

Topic: {topic}

Decision to inform: {decision}

## Decision Rule

Use `fresh` when the objective, artifact, audience, campaign, risk class, or
mode-escalation need has changed, or when an older thread has stale assumptions,
private details, fragmented prompts, or too many mixed goals.

Use `continue` only when the existing ChatGPT thread is about the same durable
object and the prior context materially lowers context cost:

- same paper, deck, slide family, code artifact, campaign, incident, or policy;
- prior attachments or wording conventions are still useful;
- thread has not drifted;
- no accidental partial send or stale input contaminated the thread;
- the next ask can be sent as one complete staged continuation message.

Use `fork_with_recap` when old context is useful but too long or drifted. Start
a new chat with only the compact context that should survive.

## Continuation Capsule Template

```text
Continuation objective:
Relevant prior context to preserve:
What changed since the prior answer:
Current evidence or artifact frame:
Hard constraints and exclusions:
Do not assume:
Deliverable:
Failure modes to check:
```

## Before UI Use

Record the chosen strategy here before sending anything:

```text
chosen_strategy: fresh | continue | fork_with_recap
existing_thread_title_or_url:
reason:
context_to_preserve:
context_to_drop:
```
"""


def build_session_strategy(topic: str, decision: str) -> str:
    starting_mode, escalation = recommend_mode(topic, decision)
    escalation_text = " -> ".join(escalation)
    collection_policy = collection_policy_for_mode(starting_mode)
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

## Response Collection Policy

- response_collection_policy: `{collection_policy}`
- policy: {render_collection_policy(starting_mode, escalation)}
- `instant` and `thinking_*` are Codex-managed collection modes.
- `pro_standard`, `pro_extended`, and `deep_research` are Justin-paste handoff modes
  after Codex submits one sanitized prompt.

## First-Message Rule

- Do not type the first ChatGPT prompt directly into the UI.
- Stage and inspect `single_message_prompt.txt`.
- Copy `clipboard_prompt.txt` to the clipboard and paste it into a fresh chat.
- Visually confirm the pasted input matches the staged prompt and contains no
  stale partial text.
- If a partial send happens, abandon that thread for this objective and start a
  fresh chat from the staged file.

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
    clipboard_prompt = build_clipboard_prompt(single_message_prompt)
    session_strategy = build_session_strategy(args.topic, args.decision)
    starting_mode, escalation = recommend_mode(args.topic, args.decision)
    collection_policy = collection_policy_for_mode(starting_mode)
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
    if len(clipboard_prompt) > MAX_SINGLE_MESSAGE_CHARS:
        print(
            f"ERROR: clipboard prompt is too long ({len(clipboard_prompt)} chars > {MAX_SINGLE_MESSAGE_CHARS})",
            file=sys.stderr,
        )
        return 1
    sensitive = scan_sensitive(clipboard_prompt)
    if sensitive:
        print(
            f"ERROR: clipboard prompt matched sensitive patterns: {', '.join(sensitive)}",
            file=sys.stderr,
        )
        return 1

    files = {
        "prompt.md": prompt,
        "single_message_prompt.txt": single_message_prompt + "\n",
        "clipboard_prompt.txt": clipboard_prompt + "\n",
        "session_strategy.md": session_strategy,
        "thread_strategy.md": build_thread_strategy(args.topic, args.decision),
        "mode_recommendation.json": build_mode_recommendation_json(args.topic, args.decision),
        "chatgpt_response.md": (
            f"# {header()}\n\n"
            f"response_collection_policy: `{collection_policy}`\n\n"
            "Paste or save the ChatGPT UI response here. For Pro, extended Pro, or Deep Research, "
            "Justin pastes the finished response back into Codex first; Codex then copies/synthesizes "
            "the useful content into this pack.\n"
        ),
        "response_collection_policy.md": build_handoff_note(args.topic, args.decision),
        "followups.md": (
            f"# {header()}\n\n"
            "First copy `clipboard_prompt.txt` to the clipboard and paste it into a fresh chat.\n"
            "Do not type or stream the first prompt directly into the ChatGPT UI.\n"
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
        if name not in {"single_message_prompt.txt", "clipboard_prompt.txt"} and header() not in text:
            errors.append(f"missing external-output marker: {path}")
        hits = scan_sensitive(text)
        if hits:
            errors.append(f"sensitive pattern(s) in {path}: {', '.join(hits)}")
    optional_policy = root / "response_collection_policy.md"
    if optional_policy.exists():
        text = optional_policy.read_text(encoding="utf-8", errors="ignore")
        if header() not in text:
            errors.append(f"missing external-output marker: {optional_policy}")
        if "response_collection_policy:" not in text:
            errors.append(f"missing response collection policy: {optional_policy}")
        hits = scan_sensitive(text)
        if hits:
            errors.append(f"sensitive pattern(s) in {optional_policy}: {', '.join(hits)}")
    optional_thread_strategy = root / "thread_strategy.md"
    if optional_thread_strategy.exists():
        text = optional_thread_strategy.read_text(encoding="utf-8", errors="ignore")
        if header() not in text:
            errors.append(f"missing external-output marker: {optional_thread_strategy}")
        if "chosen_strategy:" not in text:
            errors.append(f"missing chosen strategy field: {optional_thread_strategy}")
        hits = scan_sensitive(text)
        if hits:
            errors.append(f"sensitive pattern(s) in {optional_thread_strategy}: {', '.join(hits)}")
    prompt_text = (root / "single_message_prompt.txt").read_text(encoding="utf-8", errors="ignore")
    if len(prompt_text.strip()) > MAX_SINGLE_MESSAGE_CHARS:
        errors.append(
            f"single_message_prompt.txt exceeds {MAX_SINGLE_MESSAGE_CHARS} chars"
        )
    clipboard_path = root / "clipboard_prompt.txt"
    if clipboard_path.exists():
        clipboard_text = clipboard_path.read_text(encoding="utf-8", errors="ignore")
        if len(clipboard_text.strip()) > MAX_SINGLE_MESSAGE_CHARS:
            errors.append(
                f"clipboard_prompt.txt exceeds {MAX_SINGLE_MESSAGE_CHARS} chars"
            )
        if "\n\n\n" in clipboard_text:
            errors.append("clipboard_prompt.txt has excessive paragraph breaks")
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(f"OK: {root} validates as a private ChatGPT research pack")
    return 0


def command_clipboard(args: argparse.Namespace) -> int:
    root = Path(args.path) if args.path else latest_run()
    if root is None or not root.exists():
        print("ERROR: no ChatGPT research pack found", file=sys.stderr)
        return 1
    validation_args = argparse.Namespace(path=str(root))
    validation_status = command_validate(validation_args)
    if validation_status != 0:
        return validation_status
    prompt_path = root / "clipboard_prompt.txt"
    prompt = prompt_path.read_text(encoding="utf-8")
    try:
        subprocess.run(["pbcopy"], input=prompt, text=True, check=True)
    except (FileNotFoundError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: could not copy prompt to clipboard with pbcopy: {exc}", file=sys.stderr)
        return 1
    print(f"OK: copied {prompt_path} to clipboard")
    return 0


def command_save_response_from_clipboard(args: argparse.Namespace) -> int:
    root = Path(args.path) if args.path else latest_run()
    if root is None or not root.exists():
        print("ERROR: no ChatGPT research pack found", file=sys.stderr)
        return 1
    try:
        result = subprocess.run(["pbpaste"], text=True, capture_output=True, check=True)
    except (FileNotFoundError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: could not read clipboard with pbpaste: {exc}", file=sys.stderr)
        return 1
    response = result.stdout.strip()
    if not response:
        print("ERROR: clipboard is empty", file=sys.stderr)
        return 1
    hits = scan_sensitive(response)
    if hits:
        print(
            f"ERROR: ChatGPT response matched sensitive patterns: {', '.join(hits)}",
            file=sys.stderr,
        )
        return 1
    response_path = root / "chatgpt_response.md"
    response_path.write_text(
        f"# {header()}\n\n"
        "response_collection_policy: `codex_collects`\n\n"
        "## ChatGPT Response\n\n"
        f"{response}\n",
        encoding="utf-8",
    )
    print(f"OK: saved clipboard response to {response_path}")
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

    clipboard = subparsers.add_parser("clipboard", help="validate and copy clipboard_prompt.txt")
    clipboard.add_argument("path", nargs="?")
    clipboard.set_defaults(func=command_clipboard)

    save_response = subparsers.add_parser(
        "save-response-from-clipboard",
        help="save a copied ChatGPT response into chatgpt_response.md",
    )
    save_response.add_argument("path", nargs="?")
    save_response.set_defaults(func=command_save_response_from_clipboard)

    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
