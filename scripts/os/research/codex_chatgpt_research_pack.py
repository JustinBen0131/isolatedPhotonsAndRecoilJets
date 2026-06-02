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
from typing import Any


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
DUAL_PRO_RESPONSE_POLICY = "temporary_heartbeat_ingest"
DUAL_PRO_MODE = "pro_extended"
DUAL_PRO_SESSION_LABEL = "dual_pro_overnight_research"
DEFAULT_GITHUB_REPO = "https://github.com/JustinBen0131/isolatedPhotonsAndRecoilJets"
DUAL_PRO_PROMPT_KINDS = ("physics_thesis", "os_infrastructure")


def now_id() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def sanitize_slug(text: str) -> str:
    slug = re.sub(r"[^a-zA-Z0-9]+", "-", text.lower()).strip("-")
    return slug[:48] or "chatgpt-research"


def normalize_github_url(url: str) -> str:
    text = url.strip()
    match = re.match(r"git@github\.com:([^/]+)/(.+?)(?:\.git)?$", text)
    if match:
        return f"https://github.com/{match.group(1)}/{match.group(2)}"
    if text.endswith(".git") and text.startswith("https://github.com/"):
        return text[:-4]
    return text or DEFAULT_GITHUB_REPO


def git_text(*args: str) -> str:
    try:
        result = subprocess.run(["git", *args], text=True, capture_output=True, check=True)
    except (FileNotFoundError, subprocess.CalledProcessError):
        return ""
    return result.stdout.strip()


def repo_url_from_git() -> str:
    return normalize_github_url(git_text("remote", "get-url", "origin") or DEFAULT_GITHUB_REPO)


def git_preflight_payload(repo_url: str) -> dict[str, Any]:
    status_short = git_text("status", "--short")
    branch_status = git_text("status", "--short", "--branch")
    head = git_text("rev-parse", "HEAD")
    upstream = git_text("rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{u}")
    upstream_head = git_text("rev-parse", "@{u}") if upstream else ""
    ahead_behind = git_text("rev-list", "--left-right", "--count", "HEAD...@{u}") if upstream else ""
    return {
        "repo_url": normalize_github_url(repo_url),
        "required_before_submission": True,
        "rule": "Do not submit the OS prompt until the intended repo state has been committed and pushed.",
        "head": head,
        "upstream": upstream,
        "upstream_head": upstream_head,
        "ahead_behind": ahead_behind,
        "branch_status": branch_status,
        "working_tree_clean": status_short == "",
        "status_short": status_short,
        "suggested_pre_submit_checks": [
            "git status --short --branch",
            "git fetch origin",
            "git status --short --branch",
            "git push",
        ],
    }


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


def build_dual_physics_prompt() -> str:
    return (
        "Use ChatGPT Pro extended for a maximum-depth but compact research synthesis. "
        "Objective: find one fundamental, relevant, worthwhile physics insight for a thesis analysis "
        "centered on isolated photons and recoil jets in high-energy nuclear collisions. Context is "
        "sanitized: the analysis uses pp as a baseline, heavy-ion Au+Au as the medium-modified system, "
        "isolated direct-photon triggers, recoil-jet observables such as x_Jgamma, photon-ID/isolation/"
        "purity machinery, embedded backgrounds, response/unfolding, centrality, and systematic closure. "
        "Do not assume private project facts. Research the physics deeply enough to identify one insight "
        "that is more fundamental than a routine checklist item and would improve how the analysis is "
        "framed, checked, or remembered. Prefer concepts tied to color-charge-free photon calibration of "
        "parton energy loss, surface bias, medium response, recoil-jet energy redistribution, acoplanarity, "
        "fragmentation bias, centrality bias, isolation bias, or closure of pp-to-AuAu comparisons. "
        "Deliver exactly: 1. the single best insight in plain language; 2. why it matters for a gamma-jet "
        "thesis rather than generic jet quenching; 3. primary-source leads or canonical review/paper leads "
        "to verify; 4. one concrete null test or diagnostic that could be checked locally before treating "
        "the idea as evidence; 5. what should be recorded in daily memory/report if local verification "
        "supports it; 6. uncertainty and failure modes. Mark uncertain claims and do not present this "
        "response as project evidence."
    )


def build_dual_os_prompt(repo_url: str) -> str:
    repo = normalize_github_url(repo_url)
    return (
        "Use ChatGPT Pro extended for frontier agentic-OS infrastructure research and surgical engineering "
        f"critique. First inspect the pushed GitHub repository at {repo}. Start by understanding the current "
        "agentic OS files before recommending anything: AGENTS.md, agent_context/policies/, "
        "agent_context/CODEX_WORK_REGISTER.yaml, scripts/os/dream/codex_os_dream.py, "
        "scripts/os/research/codex_chatgpt_research_pack.py, scripts/os/context/codex_context_resonance.py, "
        "scripts/os/safety/codex_os_doctor.py, and scripts/os/README.md. If the repository is unavailable, "
        "say so and make only conditional recommendations. Objective: find exactly one fundamental "
        "infrastructure enhancement that advances the internal OS in a clean, reversible, no-regression way "
        "without changing the overall external functionality. Use premier research and best practices from "
        "agent systems, memory/retrieval, SRE, safety guardrails, workflow engines, postmortem loops, "
        "progressive delivery, and software maintainability. The change must be surgical: one small patch or "
        "runbook/validator refinement, OS-only, no science/task-status mutation, no SDCC/Condor/Drive/Slides/"
        "Gmail/Linear mutation, no secrets, no broad refactor, no behavioral churn. Deliver exactly: 1. what "
        "the current OS appears to do; 2. the one recommended enhancement and why it is fundamental; 3. the "
        "minimal target file(s) and patch shape; 4. validation commands and rollback; 5. risks or reasons not "
        "to apply it; 6. further optional improvements that should remain deferred. Mark uncertainty, cite "
        "source leads, and do not treat your answer as user approval or project evidence."
    )


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


def mode_recommendation_payload(
    topic: str,
    decision: str,
    *,
    mode_override: str | None = None,
    collection_policy_override: str | None = None,
    escalation_override: list[str] | None = None,
) -> dict[str, Any]:
    mode, escalation = recommend_mode(topic, decision)
    if mode_override:
        mode = mode_override
    if escalation_override is not None:
        escalation = escalation_override
    elif mode_override:
        escalation = [mode_override]
    features = score_mode_features(topic, decision)
    return {
        "external_output_marker": header(),
        "recommended_mode": mode,
        "mode_family": mode_family(mode),
        "escalation_ladder": escalation,
        "response_collection_policy": collection_policy_override or collection_policy_for_mode(mode),
        "features": features,
        "rule": "Use the fastest mode that satisfies the evidence burden; escalate only on missed constraints, shallow reasoning, source burden, safety sensitivity, or prior failure.",
    }


def build_mode_recommendation_json(
    topic: str,
    decision: str,
    *,
    mode_override: str | None = None,
    collection_policy_override: str | None = None,
    escalation_override: list[str] | None = None,
) -> str:
    payload = mode_recommendation_payload(
        topic,
        decision,
        mode_override=mode_override,
        collection_policy_override=collection_policy_override,
        escalation_override=escalation_override,
    )
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def render_collection_policy(
    mode: str,
    escalation: list[str],
    *,
    collection_policy_override: str | None = None,
) -> str:
    policy = collection_policy_override or collection_policy_for_mode(mode)
    if policy == DUAL_PRO_RESPONSE_POLICY:
        return (
            "Codex submits the one complete sanitized prompt in a fresh ChatGPT Pro extended conversation, "
            "records the local pack, then a single temporary local heartbeat repeats until the response is "
            "copied into chatgpt_response.md. The response remains external critique/source leads until "
            "local verification."
        )
    if policy == "user_paste_handoff":
        return (
            "Codex sends the one complete sanitized prompt, records the local pack and ChatGPT mode, "
            "then stops. Justin pastes the completed ChatGPT response back into Codex when it is done."
        )
    return (
        "Codex may send the one complete sanitized prompt, wait for the response, collect it, "
        "ask bounded follow-ups if needed, and synthesize the answer in the same Codex turn."
    )


def build_handoff_note(
    topic: str,
    decision: str,
    *,
    mode_override: str | None = None,
    collection_policy_override: str | None = None,
    escalation_override: list[str] | None = None,
) -> str:
    starting_mode, escalation = recommend_mode(topic, decision)
    if mode_override:
        starting_mode = mode_override
    if escalation_override is not None:
        escalation = escalation_override
    elif mode_override:
        escalation = [mode_override]
    policy = collection_policy_override or collection_policy_for_mode(starting_mode)
    if policy == DUAL_PRO_RESPONSE_POLICY:
        handoff_block = f"""## Temporary Heartbeat Template

Use this after submitting this explicit dual-Pro prompt:

```text
I sent the sanitized prompt to ChatGPT in {starting_mode} mode.
Local research pack: <pack path>

This pack is part of a bounded dual-Pro overnight research session. The local
temporary heartbeat will keep checking until this response is copied into
chatgpt_response.md. The response remains external critique/source leads until
Codex verifies it locally.
```
"""
    else:
        handoff_block = f"""## Pro / Deep Research Handoff Template

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
    return f"""# {header()}

## ChatGPT Response Collection Policy

- recommended_mode: `{starting_mode}`
- mode_family: `{mode_family(starting_mode)}`
- escalation_ladder: `{" -> ".join(escalation)}`
- response_collection_policy: `{policy}`

## Operating Rule

{render_collection_policy(starting_mode, escalation, collection_policy_override=collection_policy_override)}

{handoff_block}
"""


def build_thread_strategy(topic: str, decision: str, *, force_fresh: bool = False) -> str:
    default = "fresh_required" if force_fresh else "fresh"
    extra_rule = (
        "\nThis dual-Pro run requires a separate fresh conversation for this objective. "
        "Do not continue, fork, or reuse the other prompt's conversation.\n"
        if force_fresh
        else ""
    )
    return f"""# {header()}

## Thread Context Decision

Default recommendation: `{default}`.

Topic: {topic}

Decision to inform: {decision}
{extra_rule}

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


def build_session_strategy(
    topic: str,
    decision: str,
    *,
    mode_override: str | None = None,
    collection_policy_override: str | None = None,
    escalation_override: list[str] | None = None,
) -> str:
    starting_mode, escalation = recommend_mode(topic, decision)
    if mode_override:
        starting_mode = mode_override
    if escalation_override is not None:
        escalation = escalation_override
    elif mode_override:
        escalation = [mode_override]
    escalation_text = " -> ".join(escalation)
    collection_policy = collection_policy_override or collection_policy_for_mode(starting_mode)
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
- policy: {render_collection_policy(starting_mode, escalation, collection_policy_override=collection_policy_override)}
- `instant` and `thinking_*` are Codex-managed collection modes.
- `pro_standard`, `pro_extended`, and `deep_research` are Justin-paste handoff modes
  after Codex submits one sanitized prompt.
- `{DUAL_PRO_RESPONSE_POLICY}` is available only for an explicit bounded dual-Pro
  research automation: one local heartbeat repeats until both copied responses
  are present, and the responses remain unverified external critique.

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


def write_research_pack(
    run_dir: Path,
    *,
    topic: str,
    decision: str,
    prompt: str | None = None,
    single_message_prompt: str | None = None,
    mode_override: str | None = None,
    collection_policy_override: str | None = None,
    escalation_override: list[str] | None = None,
    force_fresh_thread: bool = False,
    metadata: dict[str, Any] | None = None,
) -> Path:
    run_dir.mkdir(parents=True, exist_ok=False)
    prompt = prompt or build_prompt(topic, decision)
    single_message_prompt = single_message_prompt or build_single_message_prompt(topic, decision)
    if header() not in prompt:
        prompt = f"# {header()}\n\n## Sanitized Research Prompt\n\n{prompt}\n"
    clipboard_prompt = build_clipboard_prompt(single_message_prompt)
    session_strategy = build_session_strategy(
        topic,
        decision,
        mode_override=mode_override,
        collection_policy_override=collection_policy_override,
        escalation_override=escalation_override,
    )
    starting_mode, escalation = recommend_mode(topic, decision)
    if mode_override:
        starting_mode = mode_override
    if escalation_override is not None:
        escalation = escalation_override
    elif mode_override:
        escalation = [mode_override]
    collection_policy = collection_policy_override or collection_policy_for_mode(starting_mode)
    sensitive = scan_sensitive(prompt)
    if sensitive:
        raise ValueError(f"generated prompt matched sensitive patterns: {', '.join(sensitive)}")
    if len(single_message_prompt) > MAX_SINGLE_MESSAGE_CHARS:
        raise ValueError(
            f"single-message prompt is too long ({len(single_message_prompt)} chars > {MAX_SINGLE_MESSAGE_CHARS})"
        )
    sensitive = scan_sensitive(single_message_prompt)
    if sensitive:
        raise ValueError(f"single-message prompt matched sensitive patterns: {', '.join(sensitive)}")
    if len(clipboard_prompt) > MAX_SINGLE_MESSAGE_CHARS:
        raise ValueError(f"clipboard prompt is too long ({len(clipboard_prompt)} chars > {MAX_SINGLE_MESSAGE_CHARS})")
    sensitive = scan_sensitive(clipboard_prompt)
    if sensitive:
        raise ValueError(f"clipboard prompt matched sensitive patterns: {', '.join(sensitive)}")

    files = {
        "prompt.md": prompt,
        "single_message_prompt.txt": single_message_prompt + "\n",
        "clipboard_prompt.txt": clipboard_prompt + "\n",
        "session_strategy.md": session_strategy,
        "thread_strategy.md": build_thread_strategy(topic, decision, force_fresh=force_fresh_thread),
        "mode_recommendation.json": build_mode_recommendation_json(
            topic,
            decision,
            mode_override=mode_override,
            collection_policy_override=collection_policy_override,
            escalation_override=escalation,
        ),
        "chatgpt_response.md": (
            f"# {header()}\n\n"
            f"response_collection_policy: `{collection_policy}`\n\n"
            "Paste or save the ChatGPT UI response here. For Pro, extended Pro, or Deep Research, "
            "Justin pastes the finished response back into Codex first; Codex then copies/synthesizes "
            "the useful content into this pack.\n"
        ),
        "response_collection_policy.md": build_handoff_note(
            topic,
            decision,
            mode_override=mode_override,
            collection_policy_override=collection_policy_override,
            escalation_override=escalation,
        ),
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
    if metadata:
        files["pack_metadata.json"] = json.dumps(metadata, indent=2, sort_keys=True) + "\n"
        files["pack_metadata.md"] = render_pack_metadata(metadata)
    for name, content in files.items():
        (run_dir / name).write_text(content, encoding="utf-8")
    return run_dir


def render_pack_metadata(metadata: dict[str, Any]) -> str:
    lines = [f"# {header()}", "", "## Pack Metadata", ""]
    for key, value in metadata.items():
        if isinstance(value, (dict, list)):
            lines.append(f"- {key}:")
            lines.append("")
            lines.append("```json")
            lines.append(json.dumps(value, indent=2, sort_keys=True))
            lines.append("```")
        else:
            lines.append(f"- {key}: `{value}`")
    return "\n".join(lines) + "\n"


def command_init(args: argparse.Namespace) -> int:
    run_dir = ROOT / f"{now_id()}-{sanitize_slug(args.topic)}"
    try:
        write_research_pack(run_dir, topic=args.topic, decision=args.decision)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
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


def save_response_text(root: Path, response: str, *, collection_policy: str | None = None) -> Path:
    hits = scan_sensitive(response)
    if hits:
        raise ValueError(f"ChatGPT response matched sensitive patterns: {', '.join(hits)}")
    mode_path = root / "mode_recommendation.json"
    inferred_policy = "codex_collects"
    if mode_path.exists():
        try:
            mode_payload = json.loads(mode_path.read_text(encoding="utf-8"))
            inferred_policy = str(mode_payload.get("response_collection_policy") or inferred_policy)
        except json.JSONDecodeError:
            inferred_policy = "codex_collects"
    response_path = root / "chatgpt_response.md"
    response_path.write_text(
        f"# {header()}\n\n"
        f"response_collection_policy: `{collection_policy or inferred_policy}`\n\n"
        "## ChatGPT Response\n\n"
        f"{response.strip()}\n",
        encoding="utf-8",
    )
    return response_path


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
    try:
        response_path = save_response_text(root, response)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(f"OK: saved clipboard response to {response_path}")
    return 0


def dual_session_id() -> str:
    return f"{now_id()}-{DUAL_PRO_SESSION_LABEL}"


def dual_manifest_path(session_id: str) -> Path:
    return ROOT / f"{session_id}.dual_manifest.json"


def dual_heartbeat_jsonl_path(session_id: str) -> Path:
    return ROOT / f"{session_id}.temporary_heartbeat.jsonl"


def dual_heartbeat_md_path(session_id: str) -> Path:
    return ROOT / f"{session_id}.temporary_heartbeat.md"


def dual_daily_note_path(session_id: str) -> Path:
    return ROOT / f"{session_id}.daily_report_note.md"


def dual_memory_candidate_path(session_id: str) -> Path:
    return ROOT / f"{session_id}.memory_candidate.md"


def dual_pack_specs(repo_url: str) -> list[dict[str, str]]:
    repo = normalize_github_url(repo_url)
    return [
        {
            "kind": "physics_thesis",
            "slug": "01-physics-thesis-fundamental",
            "topic": "Fundamental physics insight for isolated-photon recoil-jet thesis analysis",
            "decision": "Identify one locally verifiable physics idea worth reporting in memory and the daily update.",
            "prompt": build_dual_physics_prompt(),
            "single_message_prompt": build_dual_physics_prompt(),
        },
        {
            "kind": "os_infrastructure",
            "slug": "02-agentic-os-frontier-surgical",
            "topic": "Frontier agentic OS infrastructure research for a surgical no-regression enhancement",
            "decision": "Find exactly one reversible OS-only infrastructure improvement after inspecting the pushed GitHub repo.",
            "prompt": build_dual_os_prompt(repo),
            "single_message_prompt": build_dual_os_prompt(repo),
        },
    ]


def render_dual_manifest_markdown(manifest: dict[str, Any]) -> str:
    lines = [f"# {header()}", "", "## Dual Pro Extended Research Session", ""]
    lines.append(f"- session_id: `{manifest.get('session_id')}`")
    lines.append(f"- mode: `{manifest.get('mode')}`")
    lines.append(f"- response_collection_policy: `{manifest.get('response_collection_policy')}`")
    lines.append(f"- repo_url: `{manifest.get('repo_url')}`")
    git_preflight = manifest.get("git_preflight") if isinstance(manifest.get("git_preflight"), dict) else {}
    lines.append(f"- git working tree clean at pack creation: `{git_preflight.get('working_tree_clean')}`")
    lines.append(f"- git HEAD at pack creation: `{git_preflight.get('head')}`")
    lines.append("")
    lines.append("## Required Order")
    lines.append("- Submit `physics_thesis` first in a fresh ChatGPT Pro extended conversation.")
    lines.append("- Submit `os_infrastructure` second in a separate fresh ChatGPT Pro extended conversation.")
    lines.append("- Before submitting the OS prompt, confirm the repo has been committed and pushed.")
    lines.append("- Copy responses into the matching `chatgpt_response.md` files, then run the heartbeat command.")
    lines.append("")
    lines.append("## Packs")
    for item in manifest.get("packs") or []:
        lines.append(f"- `{item.get('kind')}`: {item.get('path')}")
    lines.append("")
    lines.append("## Heartbeat")
    heartbeat = manifest.get("heartbeat") if isinstance(manifest.get("heartbeat"), dict) else {}
    lines.append(f"- heartbeat_jsonl: `{heartbeat.get('jsonl')}`")
    lines.append(f"- heartbeat_markdown: `{heartbeat.get('markdown')}`")
    lines.append(f"- daily_report_note: `{heartbeat.get('daily_report_note')}`")
    lines.append("")
    lines.append("## Boundary")
    lines.append("- ChatGPT output is external critique/source leads only.")
    lines.append("- The physics response may become memory/report material only after local verification.")
    lines.append("- The OS response may drive at most one small reversible OS-only patch after local validation and rollback are explicit.")
    return "\n".join(lines) + "\n"


def command_init_dual_pro_overnight(args: argparse.Namespace) -> int:
    repo_url = normalize_github_url(args.repo_url or repo_url_from_git())
    session_id = args.session_id or dual_session_id()
    manifest_path = dual_manifest_path(session_id)
    if manifest_path.exists():
        print(f"ERROR: dual manifest already exists: {manifest_path}", file=sys.stderr)
        return 1
    ROOT.mkdir(parents=True, exist_ok=True)
    git_preflight = git_preflight_payload(repo_url)
    packs: list[dict[str, Any]] = []
    for spec in dual_pack_specs(repo_url):
        pack_dir = ROOT / f"{session_id}-{spec['slug']}"
        metadata = {
            "session_id": session_id,
            "prompt_kind": spec["kind"],
            "required_conversation": "fresh_separate_chat",
            "mode": DUAL_PRO_MODE,
            "response_collection_policy": DUAL_PRO_RESPONSE_POLICY,
            "submission_order": 1 if spec["kind"] == "physics_thesis" else 2,
            "repo_url": repo_url,
            "git_preflight": git_preflight,
            "boundary": "external model output is critique/source leads, not evidence or approval",
        }
        try:
            write_research_pack(
                pack_dir,
                topic=spec["topic"],
                decision=spec["decision"],
                prompt=spec["prompt"],
                single_message_prompt=spec["single_message_prompt"],
                mode_override=DUAL_PRO_MODE,
                collection_policy_override=DUAL_PRO_RESPONSE_POLICY,
                escalation_override=[DUAL_PRO_MODE],
                force_fresh_thread=True,
                metadata=metadata,
            )
        except ValueError as exc:
            print(f"ERROR: {spec['kind']} pack failed: {exc}", file=sys.stderr)
            return 1
        packs.append(
            {
                "kind": spec["kind"],
                "path": pack_dir.as_posix(),
                "clipboard_prompt": (pack_dir / "clipboard_prompt.txt").as_posix(),
                "chatgpt_response": (pack_dir / "chatgpt_response.md").as_posix(),
                "required_conversation": "fresh_separate_chat",
                "mode": DUAL_PRO_MODE,
                "response_collection_policy": DUAL_PRO_RESPONSE_POLICY,
            }
        )
    manifest = {
        "external_output_marker": header(),
        "schema_version": 1,
        "session_id": session_id,
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "mode": DUAL_PRO_MODE,
        "response_collection_policy": DUAL_PRO_RESPONSE_POLICY,
        "repo_url": repo_url,
        "git_preflight": git_preflight,
        "packs": packs,
        "heartbeat": {
            "repeat_until": "all_chatgpt_response_files_complete",
            "command": f"python3 scripts/os/research/codex_chatgpt_research_pack.py dual-heartbeat --manifest {manifest_path.as_posix()}",
            "jsonl": dual_heartbeat_jsonl_path(session_id).as_posix(),
            "markdown": dual_heartbeat_md_path(session_id).as_posix(),
            "daily_report_note": dual_daily_note_path(session_id).as_posix(),
            "memory_candidate": dual_memory_candidate_path(session_id).as_posix(),
        },
        "post_response_contract": {
            "physics_memory_rule": "record only locally verified physics insight; external response alone is not memory evidence",
            "os_change_rule": "apply at most one surgical reversible OS-only improvement after local validation, rollback, and doctor checks",
            "blocked_targets": [
                "science output or physics-status mutation",
                "SDCC or Condor mutation",
                "Google Drive/Slides/Gmail/Linear mutation",
                "secrets or credentials",
                "broad refactors or functionality changes",
            ],
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest_md = manifest_path.with_suffix(".dual_manifest.md")
    manifest_md.write_text(render_dual_manifest_markdown(manifest), encoding="utf-8")
    print(manifest_path)
    for item in packs:
        print(f"{item['kind']}={item['path']}")
    if not git_preflight.get("working_tree_clean"):
        print("WARN: working tree was dirty at pack creation; do not submit the OS prompt until commit/push is current.")
    return 0


def load_dual_manifest(path: str) -> tuple[Path, dict[str, Any]]:
    manifest_path = Path(path)
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("dual manifest is not a JSON object")
    if payload.get("response_collection_policy") != DUAL_PRO_RESPONSE_POLICY:
        raise ValueError("dual manifest has wrong response_collection_policy")
    packs = payload.get("packs")
    if not isinstance(packs, list) or len(packs) != 2:
        raise ValueError("dual manifest must contain exactly two packs")
    kinds = {str(item.get("kind")) for item in packs if isinstance(item, dict)}
    if kinds != set(DUAL_PRO_PROMPT_KINDS):
        raise ValueError(f"dual manifest prompt kinds are wrong: {sorted(kinds)}")
    return manifest_path, payload


def response_status_for_pack(pack_path: Path, *, min_chars: int = 500) -> dict[str, Any]:
    response_path = pack_path / "chatgpt_response.md"
    text = response_path.read_text(encoding="utf-8", errors="ignore") if response_path.exists() else ""
    marker = "## ChatGPT Response"
    response_text = text.split(marker, 1)[1].strip() if marker in text else ""
    placeholder = "Paste or save the ChatGPT UI response here" in text
    complete = bool(response_text and len(response_text) >= min_chars and not placeholder)
    return {
        "pack_path": pack_path.as_posix(),
        "response_path": response_path.as_posix(),
        "complete": complete,
        "response_chars": len(response_text),
        "min_response_chars": min_chars,
    }


def dual_status_payload(manifest: dict[str, Any], *, min_chars: int = 500) -> dict[str, Any]:
    rows = []
    for item in manifest.get("packs") or []:
        pack_path = Path(str(item.get("path") or ""))
        status = response_status_for_pack(pack_path, min_chars=min_chars)
        status["kind"] = item.get("kind")
        status["mode"] = item.get("mode")
        rows.append(status)
    complete = all(item.get("complete") for item in rows)
    return {
        "external_output_marker": header(),
        "session_id": manifest.get("session_id"),
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "status": "complete" if complete else "waiting",
        "complete": complete,
        "pending_kinds": [item.get("kind") for item in rows if not item.get("complete")],
        "packs": rows,
        "boundary": "external responses are critique/source leads until locally verified",
    }


def render_dual_heartbeat_markdown(status: dict[str, Any]) -> str:
    lines = [f"# {header()}", "", "## Temporary Dual-Pro Heartbeat", ""]
    lines.append(f"- session_id: `{status.get('session_id')}`")
    lines.append(f"- status: `{status.get('status')}`")
    lines.append(f"- complete: `{status.get('complete')}`")
    pending = ", ".join(str(item) for item in status.get("pending_kinds") or []) or "none"
    lines.append(f"- pending: {pending}")
    lines.append("")
    lines.append("## Pack Status")
    for item in status.get("packs") or []:
        lines.append(
            f"- `{item.get('kind')}` complete={item.get('complete')} "
            f"chars={item.get('response_chars')} path={item.get('response_path')}"
        )
    lines.append("")
    lines.append("## Boundary")
    lines.append("- Keep repeating this heartbeat until both response files are complete.")
    lines.append("- Do not treat either response as evidence, approval, or task completion without local verification.")
    return "\n".join(lines) + "\n"


def render_dual_daily_note(status: dict[str, Any], manifest: dict[str, Any]) -> str:
    lines = [f"# {header()}", "", "## Daily Report Note", ""]
    if status.get("complete"):
        lines.append("- Dual ChatGPT Pro extended research responses are copied into local packs.")
        lines.append("- Physics response is ready for local verification before memory/report promotion.")
        lines.append("- OS response is ready for one surgical OS-only patch candidate review and local validation.")
    else:
        lines.append("- Dual ChatGPT Pro extended research is still waiting on copied responses.")
        lines.append(f"- Pending prompt(s): {', '.join(status.get('pending_kinds') or []) or 'none'}")
    lines.append("")
    lines.append("## Local Packs")
    for item in manifest.get("packs") or []:
        lines.append(f"- `{item.get('kind')}`: {item.get('path')}")
    lines.append("")
    lines.append("## Allowed Promotion")
    lines.append("- Physics: source leads or insight can enter memory/daily report only after local verification.")
    lines.append("- OS: at most one reversible OS-only change may be made after validation and rollback are explicit.")
    lines.append("- Nothing here authorizes science output, task-status, SDCC, Condor, Drive/Slides/Gmail/Linear, or secret handling.")
    return "\n".join(lines) + "\n"


def render_dual_memory_candidate(status: dict[str, Any], manifest: dict[str, Any]) -> str:
    lines = [f"# {header()}", "", "## Memory Candidate", ""]
    if not status.get("complete"):
        lines.append("No memory candidate yet; both ChatGPT responses have not been copied.")
    else:
        lines.append("Candidate memory/update material exists only after Codex verifies the copied responses locally.")
        lines.append("")
        lines.append("- Physics candidate: extract one verified thesis-relevant insight plus source leads.")
        lines.append("- OS candidate: summarize the applied or deferred surgical OS improvement and validation evidence.")
    lines.append("")
    lines.append("## Source Packs")
    for item in manifest.get("packs") or []:
        lines.append(f"- `{item.get('kind')}`: {item.get('path')}")
    return "\n".join(lines) + "\n"


def command_dual_heartbeat(args: argparse.Namespace) -> int:
    manifest_path, manifest = load_dual_manifest(args.manifest)
    status = dual_status_payload(manifest, min_chars=args.min_response_chars)
    session_id = str(manifest.get("session_id"))
    heartbeat_jsonl = dual_heartbeat_jsonl_path(session_id)
    heartbeat_jsonl.parent.mkdir(parents=True, exist_ok=True)
    with heartbeat_jsonl.open("a", encoding="utf-8") as stream:
        stream.write(json.dumps(status, sort_keys=True) + "\n")
    dual_heartbeat_md_path(session_id).write_text(render_dual_heartbeat_markdown(status), encoding="utf-8")
    dual_daily_note_path(session_id).write_text(render_dual_daily_note(status, manifest), encoding="utf-8")
    dual_memory_candidate_path(session_id).write_text(render_dual_memory_candidate(status, manifest), encoding="utf-8")
    print(f"OK: dual heartbeat status={status['status']} manifest={manifest_path}")
    print(f"heartbeat={dual_heartbeat_md_path(session_id)}")
    print(f"daily_note={dual_daily_note_path(session_id)}")
    return 0 if status.get("complete") else 2


def command_validate_dual(args: argparse.Namespace) -> int:
    try:
        manifest_path, manifest = load_dual_manifest(args.manifest)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    errors: list[str] = []
    for item in manifest.get("packs") or []:
        pack_path = Path(str(item.get("path") or ""))
        status = command_validate(argparse.Namespace(path=str(pack_path)))
        if status != 0:
            errors.append(f"pack failed validation: {pack_path}")
        try:
            mode_payload = json.loads((pack_path / "mode_recommendation.json").read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            errors.append(f"mode_recommendation.json invalid in {pack_path}: {exc}")
            continue
        if mode_payload.get("recommended_mode") != DUAL_PRO_MODE:
            errors.append(f"{pack_path} is not {DUAL_PRO_MODE}")
        if mode_payload.get("response_collection_policy") != DUAL_PRO_RESPONSE_POLICY:
            errors.append(f"{pack_path} missing {DUAL_PRO_RESPONSE_POLICY}")
        thread_text = (pack_path / "thread_strategy.md").read_text(encoding="utf-8", errors="ignore")
        if "fresh_required" not in thread_text:
            errors.append(f"{pack_path} does not require a fresh separate conversation")
    os_pack = next((item for item in manifest.get("packs") or [] if item.get("kind") == "os_infrastructure"), None)
    if os_pack:
        prompt_text = (Path(os_pack["path"]) / "single_message_prompt.txt").read_text(encoding="utf-8", errors="ignore")
        for token in ("GitHub", "AGENTS.md", "scripts/os/dream/codex_os_dream.py"):
            if token not in prompt_text:
                errors.append(f"OS prompt missing required repo-context token: {token}")
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(f"OK: dual Pro research manifest validates at {manifest_path}")
    return 0


def command_save_dual_response_from_clipboard(args: argparse.Namespace) -> int:
    try:
        _, manifest = load_dual_manifest(args.manifest)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    matches = [item for item in manifest.get("packs") or [] if item.get("kind") == args.prompt_kind]
    if len(matches) != 1:
        print(f"ERROR: prompt kind not found in manifest: {args.prompt_kind}", file=sys.stderr)
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
    pack_path = Path(str(matches[0].get("path")))
    try:
        response_path = save_response_text(pack_path, response, collection_policy=DUAL_PRO_RESPONSE_POLICY)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(f"OK: saved {args.prompt_kind} response to {response_path}")
    heartbeat_status = command_dual_heartbeat(argparse.Namespace(manifest=args.manifest, min_response_chars=args.min_response_chars))
    return 0 if heartbeat_status in {0, 2} else heartbeat_status


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

    dual_init = subparsers.add_parser(
        "init-dual-pro-overnight",
        help="create two fresh-conversation ChatGPT Pro extended packs plus a temporary heartbeat manifest",
    )
    dual_init.add_argument("--repo-url", default="")
    dual_init.add_argument("--session-id", default="")
    dual_init.set_defaults(func=command_init_dual_pro_overnight)

    dual_validate = subparsers.add_parser("validate-dual", help="validate a dual Pro research manifest and its packs")
    dual_validate.add_argument("--manifest", required=True)
    dual_validate.set_defaults(func=command_validate_dual)

    dual_heartbeat = subparsers.add_parser("dual-heartbeat", help="write one temporary heartbeat status for a dual Pro session")
    dual_heartbeat.add_argument("--manifest", required=True)
    dual_heartbeat.add_argument("--min-response-chars", type=int, default=500)
    dual_heartbeat.set_defaults(func=command_dual_heartbeat)

    dual_save = subparsers.add_parser(
        "save-dual-response-from-clipboard",
        help="save clipboard text into one dual Pro pack and update the temporary heartbeat",
    )
    dual_save.add_argument("--manifest", required=True)
    dual_save.add_argument("--prompt-kind", required=True, choices=DUAL_PRO_PROMPT_KINDS)
    dual_save.add_argument("--min-response-chars", type=int, default=500)
    dual_save.set_defaults(func=command_save_dual_response_from_clipboard)

    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
