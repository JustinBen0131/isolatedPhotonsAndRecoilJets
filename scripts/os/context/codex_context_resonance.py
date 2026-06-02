#!/usr/bin/env python3
"""Resolve compact context-resonance nudges for a ThesisAnalysis task."""

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
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

from codex_work_register_common import DEFAULT_REGISTER, first_line, load_register, sorted_workstreams


SYNTHETIC_HEADER = "SYNTHETIC DREAM OUTPUT - NOT USER APPROVAL - NOT REAL USER INTENT"
WAKING_HEADER = "READ-ONLY WAKING CONTEXT - VERIFY BEFORE CLAIMS"

ARTIFACT_REGISTRY = Path("agent_context/ARTIFACT_REGISTRY.yaml")
DREAM_ROOT = Path("agent_context/local/dreams")
CHATGPT_RESEARCH_ROOT = Path("agent_context/local/chatgpt_research")
MEMORY_ROOT = Path("agent_context/memory")
CONTEXT_INDEX = MEMORY_ROOT / "CONTEXT_RESONANCE_INDEX.yaml"
SCHEMA_REGISTRY = MEMORY_ROOT / "SCHEMA_REGISTRY.yaml"
NEGATIVE_MEMORY_MAP = MEMORY_ROOT / "NEGATIVE_MEMORY_MAP.yaml"
LOCAL_CONTEXT_ROOT = Path("agent_context/local/context_resonance")
OUTCOME_LEDGER_NAME = "retrieval_outcome_ledger.jsonl"
INTERFERENCE_LEDGER_NAME = "context_interference_ledger.jsonl"
SALIENCE_INDEX_NAME = "memory_salience_index.json"

ROUTE_PROFILES = {
    "plot_generation": {
        "keywords": {"plot", "figure", "png", "qa", "histogram", "visual", "artifact", "provenance"},
        "policies": ["agent_context/policies/PLOTTING.md", "agent_context/ARTIFACT_REGISTRY.yaml"],
    },
    "sdcc_path_transfer": {
        "keywords": {"sdcc", "sftp", "transfer", "path", "condor", "bulk", "checkout", "remote"},
        "policies": [
            "agent_context/policies/SDCC_OPERATIONS.md",
            "agent_context/policies/TRANSFER_AND_PIPELINE_FILES.md",
        ],
    },
    "slide_scripting": {
        "keywords": {"slide", "deck", "speaker", "google", "presentation", "script", "png"},
        "policies": ["agent_context/policies/SLIDES_WORKFLOW.md", "agent_context/SLIDE_STYLE_MAP.md"],
    },
    "os_maintenance": {
        "keywords": {
            "dream",
            "doctor",
            "memory",
            "context",
            "resonance",
            "automation",
            "policy",
            "codex",
            "architecture",
        },
        "policies": [
            "agent_context/policies/AGENTIC_OS_HARDENING.md",
            "agent_context/policies/AGENT_MEMORY_DESIGN_NOTES.md",
        ],
    },
    "active_job_status": {
        "keywords": {"status", "job", "held", "queue", "running", "finished", "stale", "waiting"},
        "policies": ["agent_context/policies/MEMORY_AND_STATUS.md", "agent_context/CODEX_WORK_REGISTER.yaml"],
    },
}

RELATION_TYPES = {
    "same_artifact",
    "same_method",
    "same_failure_mode",
    "analogy_only",
    "warning_only",
    "visual_style",
    "path_contract",
}
EVIDENCE_CLASSES = {"real_observed", "human_approved", "derived", "synthetic"}
RETRIEVAL_POLICIES = {"conscious_context", "latent_nudge", "suppress", "quarantine_candidate"}
OUTCOME_RESULTS = {"useful", "stale", "harmful", "irrelevant", "missed"}
FEEDBACK_RESULTS = ("useful", "stale", "harmful", "irrelevant", "missed")
MAX_POSITIVE_FEEDBACK_ADJUSTMENT = 0.75
MAX_NEGATIVE_FEEDBACK_ADJUSTMENT = -1.5
GENERIC_WORKSTREAM_MATCH_TOKENS = {
    "active",
    "candidate",
    "check",
    "context",
    "current",
    "dream",
    "evidence",
    "generation",
    "maintenance",
    "ready",
    "script",
    "slide",
    "status",
    "task",
    "work",
}


def tokenize(text: str) -> set[str]:
    return {
        token
        for token in re.findall(r"[a-zA-Z0-9_+-]{4,}", text.lower())
        if token not in {"this", "that", "with", "from", "should", "would", "current", "task"}
    }


def load_yaml_dict(path: Path) -> dict[str, Any]:
    try:
        data = load_register(path)
    except Exception:
        return {}
    return data if isinstance(data, dict) else {}


def classify_task(task: str) -> dict[str, Any]:
    tokens = tokenize(task)
    scores: dict[str, int] = {}
    for route, profile in ROUTE_PROFILES.items():
        scores[route] = len(tokens & set(profile["keywords"]))
    route = max(scores, key=lambda key: scores[key])
    if scores.get(route, 0) == 0:
        route = "os_maintenance"
    return {
        "route": route,
        "tokens": sorted(tokens)[:32],
        "scores": scores,
    }


def local_source_exists(source: str) -> bool | None:
    if not source or source.startswith(("http://", "https://")):
        return None
    path = Path(source)
    if path.is_absolute():
        return None
    if source.startswith(("agent_context/", "scripts/", "codex_notes/", "macros/", "src/", "src_AuAu/")):
        return path.exists()
    return None


def source_path(item: dict[str, Any]) -> str:
    pointer = item.get("source_pointer")
    if isinstance(pointer, dict):
        return first_line(pointer.get("path"))
    return first_line(item.get("source"))


def clamp_rows(rows: list[dict[str, Any]], limit: int) -> list[dict[str, Any]]:
    rows = sorted(rows, key=lambda item: float(item.get("score", {}).get("resonance_score", 0.0)), reverse=True)
    return rows[:limit]


def empty_feedback_counts() -> dict[str, int]:
    return {result: 0 for result in FEEDBACK_RESULTS}


def bounded_feedback_adjustment(counts: dict[str, int]) -> tuple[float, str]:
    useful = max(int(counts.get("useful", 0) or 0), 0)
    stale = max(int(counts.get("stale", 0) or 0), 0)
    harmful = max(int(counts.get("harmful", 0) or 0), 0)
    irrelevant = max(int(counts.get("irrelevant", 0) or 0), 0)

    useful_boost = min(0.25 * useful, 0.75)
    stale_penalty = -min(0.45 * stale, 0.90)
    harmful_penalty = -min(0.75 * harmful, 1.50)
    irrelevant_penalty = -min(0.20 * max(irrelevant - 1, 0), 0.60)
    adjustment = useful_boost + stale_penalty + harmful_penalty + irrelevant_penalty
    adjustment = max(MAX_NEGATIVE_FEEDBACK_ADJUSTMENT, min(MAX_POSITIVE_FEEDBACK_ADJUSTMENT, adjustment))

    reasons: list[str] = []
    if useful_boost:
        reasons.append(f"useful +{useful_boost:.2f}")
    if stale_penalty:
        reasons.append(f"stale {stale_penalty:.2f}")
    if harmful_penalty:
        reasons.append(f"harmful {harmful_penalty:.2f}")
    if irrelevant_penalty:
        reasons.append(f"repeated_irrelevant {irrelevant_penalty:.2f}")
    if not reasons:
        if counts.get("missed"):
            reasons.append("missed proposal-only; no ranking boost")
        else:
            reasons.append("no local feedback")
    return round(adjustment, 3), "; ".join(reasons)


def filtered_outcome_rows(ledger_root: Path) -> tuple[list[dict[str, Any]], int]:
    raw_rows = iter_jsonl(outcome_ledger_path(ledger_root))
    rows: list[dict[str, Any]] = []
    for item in raw_rows:
        result = str(item.get("result") or "")
        if result not in OUTCOME_RESULTS:
            continue
        if str(item.get("evidence_class") or "derived") == "synthetic":
            continue
        memory_id = first_line(item.get("memory_id"))
        if not memory_id:
            continue
        row = dict(item)
        row["memory_id"] = memory_id
        row["result"] = result
        rows.append(row)
    return rows, len(raw_rows)


def aggregate_feedback_counts(outcomes: list[dict[str, Any]]) -> dict[str, dict[str, int]]:
    feedback: dict[str, dict[str, int]] = {}
    for item in outcomes:
        memory_id = first_line(item.get("memory_id"))
        result = str(item.get("result") or "")
        if not memory_id or result not in OUTCOME_RESULTS:
            continue
        counts = feedback.setdefault(memory_id, empty_feedback_counts())
        counts[result] += 1
    return feedback


def load_feedback_counts(ledger_root: Path | None = None) -> dict[str, dict[str, int]]:
    root = ledger_root or LOCAL_CONTEXT_ROOT
    outcomes, _ = filtered_outcome_rows(root)
    return aggregate_feedback_counts(outcomes)


def apply_feedback_adjustment(row: dict[str, Any], feedback_counts: dict[str, dict[str, int]]) -> dict[str, Any]:
    score = row.get("score") if isinstance(row.get("score"), dict) else {}
    memory_id = first_line(row.get("memory_id") or row.get("id"))
    counts = dict(empty_feedback_counts())
    if memory_id and memory_id in feedback_counts:
        counts.update({key: int(value or 0) for key, value in feedback_counts[memory_id].items() if key in OUTCOME_RESULTS})

    base_score = float(score.get("resonance_score", 0.0) or 0.0)
    adjustment, reason = bounded_feedback_adjustment(counts)
    final_score = round(base_score + adjustment, 3)
    score["base_resonance_score"] = round(base_score, 3)
    score["feedback_adjustment"] = adjustment
    score["final_resonance_score"] = final_score
    score["resonance_score"] = final_score
    score["feedback_counts"] = counts
    score["feedback_reason"] = reason
    row["score"] = score
    return row


def apply_feedback_to_rows(rows: list[dict[str, Any]], feedback_counts: dict[str, dict[str, int]]) -> list[dict[str, Any]]:
    for row in rows:
        if isinstance(row.get("score"), dict):
            apply_feedback_adjustment(row, feedback_counts)
    return rows


def score_record(record: dict[str, Any], task_tokens: set[str], route: str) -> dict[str, Any]:
    cue_terms = {str(term).lower() for term in record.get("cue_terms") or []}
    cue_overlap = len(task_tokens & cue_terms)
    route_bonus = 1.0 if route in (record.get("routes") or []) else 0.0
    salience = record.get("salience") if isinstance(record.get("salience"), dict) else {}
    thesis_value = float(salience.get("thesis_value", 0.5) or 0.0)
    evidence_value = float(salience.get("evidence_value", 0.5) or 0.0)
    recurrence = float(salience.get("recurrence", 0.3) or 0.0)
    error_prevention = float(salience.get("error_prevention_value", salience.get("risk_prevention", 0.5)) or 0.0)
    human_friction = float(salience.get("user_friction", 0.3) or 0.0)
    cross_domain = float(salience.get("cross_domain_analogy_value", 0.0) or 0.0)
    retrieval_cost = float(salience.get("retrieval_cost", 0.25) or 0.0)
    stale_risk = float(salience.get("stale_risk", 0.15) or 0.0)
    interference = float(salience.get("interference_risk", 0.15) or 0.0)
    synthetic_risk = float(salience.get("synthetic_contamination_risk", 0.0) or 0.0)
    evidence_class = str(record.get("evidence_class") or "derived")
    if evidence_class == "synthetic":
        synthetic_risk = max(synthetic_risk, 0.9)

    resonance_score = (
        (1.9 * cue_overlap)
        + (1.2 * route_bonus)
        + (1.0 * thesis_value)
        + (0.9 * evidence_value)
        + (0.8 * recurrence)
        + (1.3 * error_prevention)
        + (0.6 * human_friction)
        + (0.4 * cross_domain)
        - (0.7 * retrieval_cost)
        - (0.9 * stale_risk)
        - (1.0 * interference)
        - (1.4 * synthetic_risk)
    )
    return {
        "resonance_score": round(resonance_score, 3),
        "cue_overlap": cue_overlap,
        "route_bonus": route_bonus,
        "thesis_value": thesis_value,
        "risk_prevention": error_prevention,
        "token_cost": retrieval_cost,
        "stale_risk": stale_risk,
        "interference_risk": interference,
        "synthetic_contamination_risk": synthetic_risk,
    }


def route_policy_candidates(route: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for index, path in enumerate(ROUTE_PROFILES.get(route, {}).get("policies", [])):
        rows.append(
            {
                "memory_id": f"route_policy.{route}.{index + 1}",
                "source": path,
                "title": f"Route policy for {route}",
                "reason": f"direct route policy for {route}",
                "relation_type": "same_method",
                "retrieval_policy": "conscious_context",
                "evidence_class": "human_approved" if path.endswith(".md") else "real_observed",
                "allowed_use": "route policy",
                "required_waking_check": "load the file before acting if the task touches this route",
                "score": {"resonance_score": 99.0 - index, "cue_overlap": None, "route_bonus": 1.0},
            }
        )
    return rows


def live_workstream_hint(task_tokens: set[str], limit: int = 1) -> list[dict[str, Any]]:
    data = load_yaml_dict(DEFAULT_REGISTER)
    rows = []
    specific_tokens = task_tokens - GENERIC_WORKSTREAM_MATCH_TOKENS
    if not specific_tokens:
        return []
    for item in sorted_workstreams(data):
        if item.get("status") not in {"active", "running", "waiting", "blocked", "review"}:
            continue
        text = " ".join(
            first_line(value)
            for value in (
                item.get("workstream_id"),
                item.get("title"),
                item.get("current_next_action"),
                item.get("handoff_summary"),
            )
        )
        overlap = len(specific_tokens & tokenize(text))
        if overlap >= 2:
            rows.append((overlap, item))
    rows.sort(key=lambda pair: pair[0], reverse=True)
    return [
        {
            "memory_id": f"workstream.{item.get('workstream_id')}",
            "source": "agent_context/CODEX_WORK_REGISTER.yaml",
            "workstream_id": item.get("workstream_id"),
            "title": item.get("title"),
            "reason": f"live workstream has cue overlap with task ({score} specific tokens)",
            "status": item.get("status"),
            "next_action": first_line(item.get("current_next_action")),
            "relation_type": "same_artifact",
            "retrieval_policy": "conscious_context",
            "evidence_class": "real_observed",
            "allowed_use": "status orientation only",
            "required_waking_check": "verify current register evidence before changing status or claiming completion",
            "score": {"resonance_score": 18.0 + score, "cue_overlap": score, "route_bonus": 0.0},
        }
        for score, item in rows[:limit]
    ]


def artifact_hint(task_tokens: set[str]) -> list[dict[str, Any]]:
    data = load_yaml_dict(ARTIFACT_REGISTRY)
    rows = []
    specific_tokens = task_tokens - GENERIC_WORKSTREAM_MATCH_TOKENS
    if not specific_tokens:
        return []
    for item in data.get("artifacts") or []:
        if not isinstance(item, dict):
            continue
        text = " ".join(first_line(value) for value in item.values())
        overlap = len(specific_tokens & tokenize(text))
        if overlap >= 2:
            rows.append((overlap, item))
    rows.sort(key=lambda pair: pair[0], reverse=True)
    out = []
    for score, item in rows[:1]:
        out.append(
            {
                "memory_id": f"artifact.{item.get('artifact_id') or item.get('id') or 'matched'}",
                "source": "agent_context/ARTIFACT_REGISTRY.yaml",
                "artifact_id": item.get("artifact_id") or item.get("id"),
                "title": first_line(item.get("title")) or first_line(item.get("artifact_id")),
                "reason": f"registered artifact has cue overlap with task ({score} specific tokens)",
                "relation_type": "same_artifact",
                "retrieval_policy": "conscious_context",
                "evidence_class": "real_observed",
                "allowed_use": "artifact orientation only",
                "required_waking_check": "inspect artifact source, timestamp, and generation command before using as evidence",
                "score": {"resonance_score": 14.0 + score, "cue_overlap": score, "route_bonus": 0.0},
            }
        )
    return out


def registry_records(route: str, task_tokens: set[str]) -> list[dict[str, Any]]:
    data = load_yaml_dict(CONTEXT_INDEX)
    rows: list[dict[str, Any]] = []
    for raw in data.get("memory_records") or []:
        if not isinstance(raw, dict):
            continue
        score = score_record(raw, task_tokens, route)
        if score["cue_overlap"] == 0 and score["route_bonus"] == 0:
            continue
        retrieval_policy = str(raw.get("retrieval_policy") or "latent_nudge")
        evidence_class = str(raw.get("evidence_class") or "derived")
        if evidence_class == "synthetic" and retrieval_policy == "conscious_context":
            retrieval_policy = "latent_nudge"
        source = source_path(raw)
        rows.append(
            {
                "memory_id": raw.get("memory_id") or raw.get("id"),
                "source": source,
                "title": raw.get("title"),
                "reason": raw.get("why_it_surfaced") or raw.get("summary") or "registry cue matched the task",
                "why_it_surfaced": raw.get("why_it_surfaced") or raw.get("summary") or "registry cue matched the task",
                "relation_type": raw.get("relation_type") or "same_method",
                "retrieval_policy": retrieval_policy,
                "evidence_class": evidence_class,
                "allowed_use": raw.get("allowed_use") or "orientation only",
                "required_waking_check": raw.get("required_waking_check") or "verify the source pointer before using this context",
                "state": raw.get("state") if isinstance(raw.get("state"), dict) else {},
                "score": score,
                "source_exists": local_source_exists(source),
            }
        )
    return rows


def schema_candidates(route: str, task_tokens: set[str]) -> list[dict[str, Any]]:
    data = load_yaml_dict(SCHEMA_REGISTRY)
    rows: list[dict[str, Any]] = []
    for raw in data.get("schemas") or []:
        if not isinstance(raw, dict):
            continue
        score = score_record(raw, task_tokens, route)
        if score["cue_overlap"] == 0 and score["route_bonus"] == 0:
            continue
        source = source_path(raw)
        rows.append(
            {
                "memory_id": raw.get("schema_id") or raw.get("memory_id") or raw.get("id"),
                "source": source,
                "title": raw.get("title"),
                "reason": raw.get("trigger_pattern") or raw.get("summary") or "schema cue matched the task",
                "why_it_surfaced": raw.get("trigger_pattern") or raw.get("summary") or "schema cue matched the task",
                "relation_type": raw.get("relation_type") or "same_method",
                "retrieval_policy": "latent_nudge",
                "evidence_class": raw.get("evidence_class") or "derived",
                "allowed_use": raw.get("allowed_use") or "schema cue only",
                "required_waking_check": raw.get("required_waking_check") or "verify real event support before applying this schema",
                "score": score,
                "source_exists": local_source_exists(source),
            }
        )
    return rows


def negative_memory_candidates(route: str, task_tokens: set[str]) -> list[dict[str, Any]]:
    data = load_yaml_dict(NEGATIVE_MEMORY_MAP)
    rows: list[dict[str, Any]] = []
    for raw in data.get("negative_memories") or []:
        if not isinstance(raw, dict):
            continue
        score = score_record(raw, task_tokens, route)
        if score["cue_overlap"] == 0 and score["route_bonus"] == 0:
            continue
        source = source_path(raw)
        rows.append(
            {
                "memory_id": raw.get("memory_id") or raw.get("id"),
                "id": raw.get("memory_id") or raw.get("id"),
                "source": source,
                "trap": raw.get("trap"),
                "reason": raw.get("trap"),
                "relation_type": raw.get("relation_type") or "warning_only",
                "retrieval_policy": "latent_nudge",
                "evidence_class": raw.get("evidence_class") or "human_approved",
                "allowed_use": raw.get("allowed_use") or "warning only",
                "first_safe_action": raw.get("first_safe_action"),
                "required_waking_check": raw.get("required_waking_check") or "verify before using this warning",
                "score": score,
                "source_exists": local_source_exists(source),
            }
        )
    return rows


def recent_local_surface_counts() -> dict[str, Any]:
    def count_dirs(path: Path) -> int:
        if not path.exists():
            return 0
        return sum(1 for item in path.iterdir() if item.is_dir())

    return {
        "dream_run_dirs": count_dirs(DREAM_ROOT),
        "chatgpt_research_packs": count_dirs(CHATGPT_RESEARCH_ROOT),
    }


def suppressed_context(route: str, task_tokens: set[str], feedback_counts: dict[str, dict[str, int]] | None = None) -> list[dict[str, Any]]:
    counts = recent_local_surface_counts()
    rows = [
        {
            "memory_id": "suppressed.local_dreams",
            "source": "agent_context/local/dreams/",
            "reason": f"{counts['dream_run_dirs']} local dream run dirs are synthetic/proposal-only and should not enter factual context by default",
            "relation_type": "warning_only",
            "retrieval_policy": "suppress",
            "evidence_class": "synthetic",
            "allowed_use": "specific dream artifact only with synthetic boundary visible",
            "required_waking_check": "load only specific dream artifacts with synthetic boundary visible",
        },
        {
            "memory_id": "suppressed.external_research",
            "source": "agent_context/local/chatgpt_research/",
            "reason": f"{counts['chatgpt_research_packs']} research packs are external critique/source leads, not authority",
            "relation_type": "warning_only",
            "retrieval_policy": "suppress",
            "evidence_class": "derived",
            "allowed_use": "source leads and critique only",
            "required_waking_check": "verify locally before policy, code, task, or science promotion",
        },
    ]
    suppressed_rows = registry_records(route, task_tokens)
    apply_feedback_to_rows(suppressed_rows, feedback_counts or {})
    for item in suppressed_rows:
        if item.get("retrieval_policy") in {"suppress", "quarantine_candidate"} or item.get("evidence_class") == "synthetic":
            rows.append(item)
    return rows


def retrieval_outcome_candidates(route: str, task: str, selected_ids: list[str]) -> list[dict[str, Any]]:
    ordered_ids = [memory_id for memory_id in selected_ids if memory_id and not memory_id.startswith("route_policy.")]
    ordered_ids.extend(memory_id for memory_id in selected_ids if memory_id and memory_id.startswith("route_policy."))
    deduped_ids = list(dict.fromkeys(ordered_ids)) or [f"context_resonance.{route}"]
    rows = []
    for memory_id in deduped_ids[:3]:
        rows.append(
            {
                "memory_id": memory_id,
                "retrieved_for": task[:160],
                "retrieval_result": "pending_waking_use",
                "proposed_state_change": "record_useful_stale_harmful_irrelevant_or_missed_after_task",
                "requires_approval": False,
                "evidence_class": "derived",
                "required_waking_check": "after using any nudge, record whether it helped, was stale, or polluted the answer",
                "command": (
                    "python3 scripts/os/context/codex_context_resonance.py record-outcome "
                    f"--memory-id {memory_id!r} --task <task> --result useful --note <short note>"
                ),
            }
        )
    return rows


def build_context_resonance_payload(
    task: str,
    max_active: int = 3,
    max_nudges: int = 3,
    *,
    mode: str = "waking",
    ledger_root: Path | None = None,
) -> dict[str, Any]:
    signature = classify_task(task)
    task_tokens = set(signature["tokens"])
    route = str(signature["route"])
    feedback_counts = load_feedback_counts(ledger_root)

    route_rows = route_policy_candidates(route)
    route_rows.extend(live_workstream_hint(task_tokens, limit=1))
    route_rows.extend(artifact_hint(task_tokens))

    registry_rows = registry_records(route, task_tokens)
    schema_rows = schema_candidates(route, task_tokens)
    all_positive = route_rows + registry_rows + schema_rows
    apply_feedback_to_rows(all_positive, feedback_counts)
    conscious_context = clamp_rows(
        [
            item
            for item in all_positive
            if item.get("retrieval_policy") == "conscious_context" and item.get("evidence_class") != "synthetic"
        ],
        max_active,
    )
    latent_context_nudges = clamp_rows(
        [
            item
            for item in all_positive
            if item.get("retrieval_policy") == "latent_nudge" and item.get("evidence_class") != "synthetic"
        ],
        max_nudges,
    )
    negative_rows = negative_memory_candidates(route, task_tokens)
    apply_feedback_to_rows(negative_rows, feedback_counts)
    negative_memories = clamp_rows(negative_rows, max_nudges)
    selected_ids = [
        first_line(item.get("memory_id"))
        for item in (conscious_context + latent_context_nudges + negative_memories)
        if first_line(item.get("memory_id"))
    ]
    is_dream = mode == "dream"
    return {
        "synthetic": is_dream,
        "synthetic_header": SYNTHETIC_HEADER if is_dream else None,
        "waking_header": WAKING_HEADER if not is_dream else None,
        "status": "proposal_only" if is_dream else "read_only",
        "mutation_boundary": "proposal_only" if is_dream else "read_only_no_writes",
        "purpose": "surface the smallest useful old context without expanding boot context or promoting synthetic material",
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "task": task,
        "task_signature": signature,
        "limits": {
            "max_conscious_context": max_active,
            "max_active_facts": max_active,
            "max_latent_nudges": max_nudges,
            "max_negative_memories": max_nudges,
        },
        "conscious_context": conscious_context,
        "active_facts": conscious_context,
        "latent_context_nudges": latent_context_nudges,
        "negative_memories": negative_memories,
        "suppressed_context": suppressed_context(route, task_tokens, feedback_counts),
        "retrieval_outcome_candidates": retrieval_outcome_candidates(route, task, selected_ids),
        "required_checks_before_claim": [
            "load route policy before acting",
            "verify local evidence before factual claims",
            "treat dream and external-model material as proposal/critique until verified",
            "record retrieval outcome if a surfaced memory materially helped or polluted the task",
        ],
        "context_interference_budget": {
            "max_conscious_context": max_active,
            "max_latent_nudges": max_nudges,
            "failure_condition": "resolver output becomes a larger boot pack or surfaces stale/synthetic material as proof",
        },
    }


def render_markdown(payload: dict[str, Any]) -> str:
    header = SYNTHETIC_HEADER if payload.get("synthetic") else WAKING_HEADER
    boundary = "proposal-only, no live mutation" if payload.get("synthetic") else "read-only, verify before claims"
    lines = ["# Context Resonance Review", "", header, ""]
    signature = payload.get("task_signature") if isinstance(payload.get("task_signature"), dict) else {}
    lines.append(f"- task: {payload.get('task')}")
    lines.append(f"- route: `{signature.get('route')}`")
    lines.append(f"- boundary: {boundary}")
    lines.append("")
    lines.append("## Conscious Context")
    for item in payload.get("conscious_context") or payload.get("active_facts") or []:
        lines.append(f"- `{item.get('memory_id')}` -> `{item.get('source')}`: {item.get('reason')}")
    if not payload.get("conscious_context") and not payload.get("active_facts"):
        lines.append("- none")
    lines.append("")
    lines.append("## Latent Context Nudges")
    for item in payload.get("latent_context_nudges") or []:
        lines.append(f"- `{item.get('memory_id')}` ({item.get('relation_type')}): {item.get('why_it_surfaced') or item.get('reason')}")
    if not payload.get("latent_context_nudges"):
        lines.append("- none")
    lines.append("")
    lines.append("## Negative Memories")
    for item in payload.get("negative_memories") or []:
        lines.append(f"- `{item.get('memory_id')}`: {item.get('trap')}; first safe action: {item.get('first_safe_action')}")
    if not payload.get("negative_memories"):
        lines.append("- none")
    lines.append("")
    lines.append("## Suppressed Context")
    for item in payload.get("suppressed_context") or []:
        lines.append(f"- `{item.get('memory_id')}` -> `{item.get('source')}`: {item.get('reason')}")
    lines.append("")
    lines.append("## Required Checks Before Claim")
    for item in payload.get("required_checks_before_claim") or []:
        lines.append(f"- {item}")
    return "\n".join(lines) + "\n"


def outcome_ledger_path(ledger_root: Path) -> Path:
    return ledger_root / OUTCOME_LEDGER_NAME


def record_outcome(args: argparse.Namespace) -> int:
    ledger_root = Path(args.ledger_root) if args.ledger_root else LOCAL_CONTEXT_ROOT
    ledger_root.mkdir(parents=True, exist_ok=True)
    row = {
        "recorded_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "memory_id": args.memory_id,
        "task": args.task[:500],
        "result": args.result,
        "note": args.note[:500],
        "evidence_class": "derived",
        "ledger_policy": "local_only",
    }
    with outcome_ledger_path(ledger_root).open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(row, sort_keys=True) + "\n")
    print(json.dumps({"status": "recorded", "ledger": outcome_ledger_path(ledger_root).as_posix(), "row": row}, indent=2))
    return 0


def iter_jsonl(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if not path.exists():
        return rows
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            text = line.strip()
            if not text:
                continue
            try:
                item = json.loads(text)
            except json.JSONDecodeError:
                continue
            if isinstance(item, dict):
                rows.append(item)
    return rows


def load_json_file(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return payload if isinstance(payload, dict) else {}


def homeostasis_action_for_counts(counts: dict[str, int]) -> tuple[str, str] | None:
    useful = int(counts.get("useful", 0) or 0)
    stale = int(counts.get("stale", 0) or 0)
    harmful = int(counts.get("harmful", 0) or 0)
    irrelevant = int(counts.get("irrelevant", 0) or 0)
    missed = int(counts.get("missed", 0) or 0)
    if harmful >= 2 or (harmful and stale):
        return "suppress_candidate", "repeated harmful/stale retrieval feedback"
    if harmful >= 1:
        return "quarantine_candidate", "harmful retrieval feedback"
    if stale >= 1 or irrelevant >= 2:
        return "cool", "stale or repeatedly irrelevant retrieval feedback"
    if missed >= 2:
        return "promote", "repeated missed-memory feedback; review cue terms or tracked record"
    if missed >= 1:
        return "cue_review", "missed-memory feedback; no automatic ranking boost"
    if useful >= 2:
        return "keep_hot", "repeated useful retrieval feedback"
    if useful >= 1:
        return "keep_latent", "single useful retrieval feedback"
    return None


def build_homeostasis_proposals(feedback_counts: dict[str, dict[str, int]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for memory_id in sorted(feedback_counts):
        counts = feedback_counts[memory_id]
        action_reason = homeostasis_action_for_counts(counts)
        if not action_reason:
            continue
        action, reason = action_reason
        rows.append(
            {
                "memory_id": memory_id,
                "action": action,
                "reason": reason,
                "feedback_counts": dict(counts),
                "requires_waking_validation": True,
                "mutation_boundary": "proposal_only",
                "allowed_use": "routing homeostasis proposal only; do not edit tracked memory automatically",
            }
        )
    return rows


def local_maintenance_review(ledger_root: Path | None = None) -> dict[str, Any]:
    root = ledger_root or LOCAL_CONTEXT_ROOT
    outcomes, raw_outcome_count = filtered_outcome_rows(root)
    feedback_counts = aggregate_feedback_counts(outcomes)
    interference = iter_jsonl(root / INTERFERENCE_LEDGER_NAME)
    salience = load_json_file(root / SALIENCE_INDEX_NAME)
    salience_records = salience.get("records") if isinstance(salience.get("records"), list) else []

    stale_or_harmful = [
        item
        for item in outcomes
        if item.get("result") in {"stale", "harmful"} and item.get("evidence_class") != "synthetic"
    ]
    missed_or_useful = [
        item
        for item in outcomes
        if item.get("result") in {"missed", "useful"} and item.get("evidence_class") != "synthetic"
    ]
    pollution_counts = Counter(
        str(item.get("memory_id") or item.get("context_id") or "unknown")
        for item in stale_or_harmful + interference
        if str(item.get("evidence_class") or "derived") != "synthetic"
    )
    useful_trap_counts = Counter(str(item.get("memory_id") or "unknown") for item in missed_or_useful)
    source_counts = Counter(str(item.get("source") or "") for item in salience_records if item.get("source"))
    duplicate_sources = {source: count for source, count in source_counts.items() if count > 1}

    evidence: list[str] = []
    for item in stale_or_harmful[-3:]:
        evidence.append(
            f"retrieval_outcome:{item.get('memory_id')} result={item.get('result')} note={str(item.get('note') or '')[:120]}"
        )
    for item in interference[-3:]:
        evidence.append(
            f"context_interference:{item.get('memory_id') or item.get('context_id') or 'unknown'} note={str(item.get('note') or item.get('reason') or '')[:120]}"
        )
    if not salience_records:
        evidence.append("memory_salience_index: no local records available for route/index coherence")

    return {
        "status": "ok",
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "ledger_root": root.as_posix(),
        "outcome_count": len(outcomes),
        "raw_outcome_count": raw_outcome_count,
        "ignored_outcome_count": raw_outcome_count - len(outcomes),
        "interference_count": len(interference),
        "salience_record_count": len(salience_records),
        "feedback_counts": feedback_counts,
        "homeostasis_proposals": build_homeostasis_proposals(feedback_counts),
        "pollution_memory_counts": dict(pollution_counts.most_common(8)),
        "useful_memory_counts": dict(useful_trap_counts.most_common(8)),
        "duplicate_source_counts": duplicate_sources,
        "route_index_refresh_needed": not bool(salience_records),
        "candidate_evidence": evidence,
        "maintenance_policy": "local-only distilled routing evidence; synthetic rows cannot justify auto-apply",
    }


def review_payload(ledger_root: Path) -> dict[str, Any]:
    outcomes, raw_outcome_count = filtered_outcome_rows(ledger_root)
    feedback_counts = aggregate_feedback_counts(outcomes)
    interference = iter_jsonl(ledger_root / INTERFERENCE_LEDGER_NAME)
    result_counts = Counter(str(item.get("result") or "unknown") for item in outcomes)
    memory_counts = Counter(str(item.get("memory_id") or "unknown") for item in outcomes)
    harmful = [item for item in outcomes if item.get("result") in {"harmful", "stale"}]
    return {
        "status": "ok",
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "ledger_root": ledger_root.as_posix(),
        "outcome_count": len(outcomes),
        "raw_outcome_count": raw_outcome_count,
        "ignored_outcome_count": raw_outcome_count - len(outcomes),
        "interference_count": len(interference),
        "result_counts": dict(result_counts),
        "feedback_counts": feedback_counts,
        "top_memory_ids": dict(memory_counts.most_common(8)),
        "cooldown_candidates": harmful[-8:],
        "homeostasis_proposals": build_homeostasis_proposals(feedback_counts),
        "review_policy": "feedback is local-only routing evidence, not factual project proof",
    }


def maintenance_review(args: argparse.Namespace) -> int:
    ledger_root = Path(args.ledger_root) if args.ledger_root else LOCAL_CONTEXT_ROOT
    payload = local_maintenance_review(ledger_root)
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print("# Context Resonance Maintenance Review")
        print(f"- ledger: `{payload['ledger_root']}`")
        print(f"- outcomes: {payload['outcome_count']}")
        print(f"- interference: {payload['interference_count']}")
        print(f"- salience records: {payload['salience_record_count']}")
        print(f"- route/index refresh needed: `{payload['route_index_refresh_needed']}`")
    return 0


def review(args: argparse.Namespace) -> int:
    ledger_root = Path(args.ledger_root) if args.ledger_root else LOCAL_CONTEXT_ROOT
    payload = review_payload(ledger_root)
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print("# Context Resonance Feedback Review")
        print(f"- ledger: `{payload['ledger_root']}`")
        print(f"- outcomes: {payload['outcome_count']}")
        print(f"- result counts: {payload['result_counts']}")
    return 0


CANARY_TASKS = [
    {
        "id": "plot_generation",
        "task": "plot generation for a slide-ready PNG with provenance",
        "route": "plot_generation",
        "expected_any": {"plot_readiness_requires_provenance", "no_ready_without_evidence"},
    },
    {
        "id": "sdcc_transfer",
        "task": "SDCC path transfer map update without remote mutation",
        "route": "sdcc_path_transfer",
        "expected_any": {"sdcc_path_contract_drift", "no_remote_cleanup_without_guard"},
    },
    {
        "id": "slide_scripting",
        "task": "slide scripting helper for reusable presentation PNGs",
        "route": "slide_scripting",
        "expected_any": {"slide_thread_reuse_capsule", "plot_readiness_requires_provenance"},
    },
    {
        "id": "active_status",
        "task": "active job status check and stale-state refresh",
        "route": "active_job_status",
        "expected_any": {"workstream_refresh_contract", "no_generic_status_context_bleed"},
    },
    {
        "id": "memory_architecture",
        "task": "OS memory architecture and context resonance upgrade",
        "route": "os_maintenance",
        "expected_any": {"context_resonance_architecture", "dream_memory_not_more_context"},
    },
]


def canary_payload() -> dict[str, Any]:
    failures: list[str] = []
    cases: list[dict[str, Any]] = []
    with TemporaryDirectory(prefix="context_resonance_canary_") as temp_dir:
        empty_ledger_root = Path(temp_dir)
        for spec in CANARY_TASKS:
            payload = build_context_resonance_payload(
                spec["task"],
                max_active=3,
                max_nudges=3,
                mode="waking",
                ledger_root=empty_ledger_root,
            )
            signature = payload.get("task_signature") if isinstance(payload.get("task_signature"), dict) else {}
            ids = {
                first_line(item.get("memory_id"))
                for key in ("conscious_context", "latent_context_nudges", "negative_memories")
                for item in payload.get(key, [])
                if isinstance(item, dict)
            }
            if signature.get("route") != spec["route"]:
                failures.append(f"{spec['id']}: expected route {spec['route']} got {signature.get('route')}")
            if len(payload.get("conscious_context") or []) > 3:
                failures.append(f"{spec['id']}: conscious_context cap exceeded")
            if len(payload.get("latent_context_nudges") or []) > 3:
                failures.append(f"{spec['id']}: latent_context_nudges cap exceeded")
            if len(payload.get("negative_memories") or []) > 3:
                failures.append(f"{spec['id']}: negative_memories cap exceeded")
            if not ids & set(spec["expected_any"]):
                failures.append(f"{spec['id']}: expected one of {sorted(spec['expected_any'])}, got {sorted(ids)}")
            for key in ("conscious_context", "latent_context_nudges", "negative_memories"):
                for item in payload.get(key) or []:
                    score = item.get("score") if isinstance(item.get("score"), dict) else {}
                    for field in (
                        "base_resonance_score",
                        "feedback_adjustment",
                        "final_resonance_score",
                        "feedback_counts",
                        "feedback_reason",
                    ):
                        if field not in score:
                            failures.append(f"{spec['id']}: {item.get('memory_id')} lacks score.{field}")
            for item in payload.get("conscious_context") or []:
                if item.get("evidence_class") == "synthetic":
                    failures.append(f"{spec['id']}: synthetic material surfaced as conscious context")
                if not item.get("required_waking_check"):
                    failures.append(f"{spec['id']}: conscious item lacks required_waking_check")
            cases.append(
                {
                    "id": spec["id"],
                    "route": signature.get("route"),
                    "memory_ids": sorted(ids),
                    "caps": {
                        "conscious_context": len(payload.get("conscious_context") or []),
                        "latent_context_nudges": len(payload.get("latent_context_nudges") or []),
                        "negative_memories": len(payload.get("negative_memories") or []),
                    },
                }
            )

    feedback_cases = feedback_canary_cases()
    failures.extend(feedback_cases.get("failures") or [])
    cases.extend(feedback_cases.get("cases") or [])
    return {
        "status": "fail" if failures else "pass",
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "cases": cases,
        "failures": failures,
    }


def write_canary_outcomes(ledger_root: Path, rows: list[dict[str, Any]]) -> None:
    ledger_root.mkdir(parents=True, exist_ok=True)
    with outcome_ledger_path(ledger_root).open("a", encoding="utf-8") as handle:
        for row in rows:
            payload = {
                "recorded_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "memory_id": row["memory_id"],
                "task": row.get("task", "context resonance canary"),
                "result": row["result"],
                "note": row.get("note", "canary fixture"),
                "evidence_class": row.get("evidence_class", "derived"),
                "ledger_policy": "local_only",
            }
            handle.write(json.dumps(payload, sort_keys=True) + "\n")


def find_payload_item(payload: dict[str, Any], section: str, memory_id: str) -> dict[str, Any]:
    for item in payload.get(section) or []:
        if isinstance(item, dict) and first_line(item.get("memory_id")) == memory_id:
            return item
    return {}


def section_position(payload: dict[str, Any], section: str, memory_id: str) -> int | None:
    for index, item in enumerate(payload.get(section) or []):
        if isinstance(item, dict) and first_line(item.get("memory_id")) == memory_id:
            return index
    return None


def feedback_canary_cases() -> dict[str, Any]:
    failures: list[str] = []
    cases: list[dict[str, Any]] = []
    task = "OS memory architecture and context resonance upgrade"

    with TemporaryDirectory(prefix="context_resonance_feedback_boost_") as temp_dir:
        root = Path(temp_dir)
        base_payload = build_context_resonance_payload(task, max_active=3, max_nudges=3, mode="waking", ledger_root=root)
        base_position = section_position(base_payload, "latent_context_nudges", "context_resonance_schema")
        write_canary_outcomes(
            root,
            [
                {"memory_id": "context_resonance_schema", "result": "useful"},
                {"memory_id": "context_resonance_schema", "result": "useful"},
            ],
        )
        boosted_payload = build_context_resonance_payload(task, max_active=3, max_nudges=3, mode="waking", ledger_root=root)
        boosted_position = section_position(boosted_payload, "latent_context_nudges", "context_resonance_schema")
        boosted_item = find_payload_item(boosted_payload, "latent_context_nudges", "context_resonance_schema")
        boosted_score = boosted_item.get("score") if isinstance(boosted_item.get("score"), dict) else {}
        if base_position is None or boosted_position is None or boosted_position >= base_position:
            failures.append("feedback_useful_boost: useful memory did not rise within latent_context_nudges")
        if float(boosted_score.get("feedback_adjustment", 0.0) or 0.0) <= 0:
            failures.append("feedback_useful_boost: useful memory did not receive positive feedback adjustment")
        cases.append(
            {
                "id": "feedback_useful_boost",
                "memory_id": "context_resonance_schema",
                "base_position": base_position,
                "boosted_position": boosted_position,
                "feedback_adjustment": boosted_score.get("feedback_adjustment"),
            }
        )

    with TemporaryDirectory(prefix="context_resonance_feedback_penalty_") as temp_dir:
        root = Path(temp_dir)
        write_canary_outcomes(
            root,
            [
                {"memory_id": "context_resonance_architecture", "result": "harmful"},
                {"memory_id": "context_resonance_architecture", "result": "stale"},
            ],
        )
        payload = build_context_resonance_payload(task, max_active=3, max_nudges=3, mode="waking", ledger_root=root)
        item = find_payload_item(payload, "latent_context_nudges", "context_resonance_architecture")
        score = item.get("score") if isinstance(item.get("score"), dict) else {}
        if float(score.get("feedback_adjustment", 0.0) or 0.0) >= 0:
            failures.append("feedback_harmful_stale_penalty: harmful/stale memory was not penalized")
        if float(score.get("final_resonance_score", 0.0) or 0.0) >= float(score.get("base_resonance_score", 0.0) or 0.0):
            failures.append("feedback_harmful_stale_penalty: final score did not fall below base score")
        cases.append(
            {
                "id": "feedback_harmful_stale_penalty",
                "memory_id": "context_resonance_architecture",
                "feedback_adjustment": score.get("feedback_adjustment"),
            }
        )

    irrelevant_row = {"memory_id": "irrelevant_dummy", "score": {"resonance_score": 10.0}}
    apply_feedback_adjustment(
        irrelevant_row,
        {"irrelevant_dummy": {**empty_feedback_counts(), "irrelevant": 3}},
    )
    irrelevant_score = irrelevant_row["score"]
    if float(irrelevant_score.get("feedback_adjustment", 0.0) or 0.0) >= 0:
        failures.append("feedback_repeated_irrelevant_cools: repeated irrelevant memory was not cooled")
    cases.append(
        {
            "id": "feedback_repeated_irrelevant_cools",
            "memory_id": "irrelevant_dummy",
            "feedback_adjustment": irrelevant_score.get("feedback_adjustment"),
        }
    )

    with TemporaryDirectory(prefix="context_resonance_feedback_missed_") as temp_dir:
        root = Path(temp_dir)
        write_canary_outcomes(root, [{"memory_id": "context_resonance_schema", "result": "missed"}])
        payload = build_context_resonance_payload(task, max_active=3, max_nudges=3, mode="waking", ledger_root=root)
        signature = payload.get("task_signature") if isinstance(payload.get("task_signature"), dict) else {}
        item = find_payload_item(payload, "latent_context_nudges", "context_resonance_schema")
        score = item.get("score") if isinstance(item.get("score"), dict) else {}
        review = review_payload(root)
        proposals = review.get("homeostasis_proposals") if isinstance(review.get("homeostasis_proposals"), list) else []
        actions = {item.get("action") for item in proposals if isinstance(item, dict) and item.get("memory_id") == "context_resonance_schema"}
        if float(score.get("feedback_adjustment", 0.0) or 0.0) != 0.0:
            failures.append("feedback_missed_review_only: missed memory changed ranking")
        if not actions & {"cue_review", "promote"}:
            failures.append("feedback_missed_review_only: missed memory did not produce homeostasis cue")
        if signature.get("route") != "os_maintenance":
            failures.append(f"feedback_missed_review_only: expected os_maintenance route got {signature.get('route')}")
        cases.append(
            {
                "id": "feedback_missed_review_only",
                "memory_id": "context_resonance_schema",
                "route": signature.get("route"),
                "feedback_adjustment": score.get("feedback_adjustment"),
                "homeostasis_actions": sorted(str(action) for action in actions if action),
            }
        )
    return {
        "cases": cases,
        "failures": failures,
    }


def canary(args: argparse.Namespace) -> int:
    payload = canary_payload()
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print("# Context Resonance Canary")
        for case in payload["cases"]:
            print(f"- {case['id']}: route={case['route']}, ids={', '.join(case['memory_ids'])}")
        for failure in payload["failures"]:
            print(f"FAIL: {failure}")
    return 1 if payload["status"] != "pass" else 0


def resolve(args: argparse.Namespace) -> int:
    ledger_root = Path(args.ledger_root) if getattr(args, "ledger_root", None) else LOCAL_CONTEXT_ROOT
    payload = build_context_resonance_payload(
        args.task,
        max_active=args.max_active,
        max_nudges=args.max_nudges,
        mode=args.mode,
        ledger_root=ledger_root,
    )
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(render_markdown(payload), end="")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command")

    resolve_parser = subparsers.add_parser("resolve", help="resolve compact context for a task")
    resolve_parser.add_argument("--task", required=True, help="short natural-language task description")
    resolve_parser.add_argument("--json", action="store_true", help="emit JSON instead of Markdown")
    resolve_parser.add_argument("--max-active", type=int, default=3)
    resolve_parser.add_argument("--max-nudges", type=int, default=3)
    resolve_parser.add_argument("--mode", choices=["waking", "dream"], default="waking")
    resolve_parser.add_argument("--ledger-root", help="local context-resonance ledger root for deterministic testing")
    resolve_parser.set_defaults(func=resolve)

    record_parser = subparsers.add_parser("record-outcome", help="append local retrieval feedback")
    record_parser.add_argument("--memory-id", required=True)
    record_parser.add_argument("--task", required=True)
    record_parser.add_argument("--result", required=True, choices=sorted(OUTCOME_RESULTS))
    record_parser.add_argument("--note", required=True)
    record_parser.add_argument("--ledger-root")
    record_parser.set_defaults(func=record_outcome)

    review_parser = subparsers.add_parser("review", help="review local retrieval feedback")
    review_parser.add_argument("--ledger-root")
    review_parser.add_argument("--json", action="store_true")
    review_parser.set_defaults(func=review)

    maintenance_parser = subparsers.add_parser("maintenance-review", help="summarize local maintenance pressure")
    maintenance_parser.add_argument("--ledger-root")
    maintenance_parser.add_argument("--json", action="store_true")
    maintenance_parser.set_defaults(func=maintenance_review)

    canary_parser = subparsers.add_parser("canary", help="run compact routing canaries")
    canary_parser.add_argument("--json", action="store_true")
    canary_parser.set_defaults(func=canary)
    return parser


def legacy_resolve(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task", required=True, help="short natural-language task description")
    parser.add_argument("--json", action="store_true", help="emit JSON instead of Markdown")
    parser.add_argument("--max-active", type=int, default=3)
    parser.add_argument("--max-nudges", type=int, default=3)
    parser.add_argument("--mode", choices=["waking", "dream"], default="waking")
    parser.add_argument("--ledger-root")
    args = parser.parse_args(argv)
    return resolve(args)


def main(argv: list[str] | None = None) -> int:
    raw_args = list(sys.argv[1:] if argv is None else argv)
    commands = {"resolve", "record-outcome", "review", "maintenance-review", "canary"}
    if not raw_args or raw_args[0].startswith("-"):
        return legacy_resolve(raw_args)
    if raw_args[0] not in commands:
        return legacy_resolve(raw_args)
    parser = build_parser()
    args = parser.parse_args(raw_args)
    if not hasattr(args, "func"):
        parser.print_help()
        return 2
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
