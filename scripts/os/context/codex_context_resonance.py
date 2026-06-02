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
import hashlib
import json
import re
import sys
from collections import Counter
from datetime import datetime, timedelta, timezone
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
FEEDBACK_RECENT_DAYS = 30
FEEDBACK_AGING_DAYS = 90
FEEDBACK_OLD_DAYS = 180
CROSS_ROUTE_FEEDBACK_WEIGHT = 0.25
WEAK_EVIDENCE_POSITIVE_MULTIPLIER = 0.5
WEAK_EVIDENCE_NEGATIVE_MULTIPLIER = 0.2
PRIVATE_TEXT_PATTERNS = (
    "password",
    "passwd",
    "secret",
    "token",
    "api_key",
    "private key",
    "ssh-rsa",
    "aws_",
    "openai_api_key",
    "sphnxuser",
    "ssh.sdcc",
    "/sphenix/u/",
)
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


def empty_weighted_feedback_counts() -> dict[str, float]:
    return {result: 0.0 for result in FEEDBACK_RESULTS}


def bounded_feedback_adjustment(counts: dict[str, float | int]) -> tuple[float, str]:
    useful = max(float(counts.get("useful", 0.0) or 0.0), 0.0)
    stale = max(float(counts.get("stale", 0.0) or 0.0), 0.0)
    harmful = max(float(counts.get("harmful", 0.0) or 0.0), 0.0)
    irrelevant = max(float(counts.get("irrelevant", 0.0) or 0.0), 0.0)

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


def parse_recorded_at(value: Any) -> datetime | None:
    text = first_line(value)
    if not text:
        return None
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    try:
        parsed = datetime.fromisoformat(text)
    except ValueError:
        return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


def feedback_decay_for_recorded_at(value: Any, now: datetime | None = None) -> dict[str, Any]:
    parsed = parse_recorded_at(value)
    if parsed is None:
        return {"bucket": "malformed_recorded_at", "weight": 0.0, "age_days": None}
    reference = now or datetime.now(timezone.utc)
    age_days = max((reference - parsed).total_seconds() / 86400.0, 0.0)
    if age_days <= FEEDBACK_RECENT_DAYS:
        bucket, weight = "recent_full_effect", 1.0
    elif age_days <= FEEDBACK_AGING_DAYS:
        bucket, weight = "aging_partial_effect", 0.65
    elif age_days <= FEEDBACK_OLD_DAYS:
        bucket, weight = "old_weak_effect", 0.35
    else:
        bucket, weight = "very_old_review_only", 0.0
    return {"bucket": bucket, "weight": weight, "age_days": round(age_days, 1)}


def evidence_ref_for(item: dict[str, Any]) -> str:
    return first_line(item.get("evidence_ref") or item.get("evidence") or item.get("source") or "")


def row_contains_suspicious_text(item: dict[str, Any]) -> bool:
    text = " ".join(
        first_line(item.get(key))
        for key in ("memory_id", "task", "note", "evidence_ref", "evidence", "source")
    ).lower()
    return any(pattern in text for pattern in PRIVATE_TEXT_PATTERNS)


def evidence_quality_for(item: dict[str, Any]) -> str:
    evidence_ref = evidence_ref_for(item)
    if not evidence_ref:
        return "missing"
    lowered = evidence_ref.lower()
    if row_contains_suspicious_text(item):
        return "weak_private_text"
    if len(evidence_ref) < 8 or lowered in {"none", "n/a", "na", "unknown", "todo", "note"}:
        return "weak"
    if len(evidence_ref) < 12 and not any(marker in evidence_ref for marker in ("/", ":", "#", ".")):
        return "weak"
    return "strong"


def evidence_multiplier(result: str, quality: str) -> float:
    if quality == "strong":
        return 1.0
    if result in {"harmful", "stale"}:
        return WEAK_EVIDENCE_NEGATIVE_MULTIPLIER
    return WEAK_EVIDENCE_POSITIVE_MULTIPLIER


def task_fingerprint(task: Any) -> str:
    tokens = sorted(tokenize(first_line(task)))
    digest = hashlib.sha256(" ".join(tokens).encode("utf-8")).hexdigest()
    return digest[:12]


def feedback_route_for(item: dict[str, Any]) -> tuple[str, bool]:
    route = first_line(item.get("route"))
    if route in ROUTE_PROFILES:
        return route, True
    return "unknown", False


def add_audit_warning(warnings: list[dict[str, Any]], row_index: int, kind: str, item: dict[str, Any], detail: str) -> None:
    warnings.append(
        {
            "row_index": row_index,
            "kind": kind,
            "memory_id": first_line(item.get("memory_id")),
            "result": first_line(item.get("result")),
            "route": first_line(item.get("route")),
            "detail": detail,
        }
    )


def infer_feedback_trainable(row: dict[str, Any]) -> bool:
    memory_id = first_line(row.get("memory_id") or row.get("id"))
    source = first_line(row.get("source"))
    evidence_class = str(row.get("evidence_class") or "derived")
    retrieval_policy = str(row.get("retrieval_policy") or "latent_nudge")
    if evidence_class == "synthetic":
        return False
    if retrieval_policy in {"suppress", "quarantine_candidate"}:
        return False
    if memory_id.startswith(("route_policy.", "workstream.", "artifact.", "suppressed.")):
        return False
    if source.startswith(("agent_context/local/dreams", "agent_context/local/chatgpt_research")):
        return False
    if "external_model" in memory_id or "dream" in memory_id and evidence_class == "synthetic":
        return False
    return True


def feedback_detail_empty() -> dict[str, Any]:
    return {
        "global_counts": empty_feedback_counts(),
        "global_weighted_counts": empty_weighted_feedback_counts(),
        "route_counts": {},
        "route_weighted_counts": {},
        "decay_buckets": {},
        "evidence_quality_counts": {},
        "review_only_count": 0,
        "weak_evidence_negative_count": 0,
        "strong_evidence_refs": [],
        "task_fingerprints": [],
        "routes": [],
        "days": [],
    }


def round_weighted_counts(counts: dict[str, float | int]) -> dict[str, float]:
    return {result: round(float(counts.get(result, 0.0) or 0.0), 3) for result in FEEDBACK_RESULTS}


def normalize_feedback_detail(value: dict[str, Any] | None) -> dict[str, Any]:
    if not isinstance(value, dict):
        return feedback_detail_empty()
    if "global_counts" in value or "route_counts" in value:
        detail = feedback_detail_empty()
        detail.update(value)
        detail["global_counts"] = {**empty_feedback_counts(), **dict(detail.get("global_counts") or {})}
        detail["global_weighted_counts"] = {
            **empty_weighted_feedback_counts(),
            **dict(detail.get("global_weighted_counts") or {}),
        }
        return detail
    detail = feedback_detail_empty()
    for result in FEEDBACK_RESULTS:
        raw_value = value.get(result, 0) or 0
        detail["global_counts"][result] = int(raw_value)
        detail["global_weighted_counts"][result] = float(raw_value)
    return detail


def inspect_outcome_ledger(ledger_root: Path) -> dict[str, Any]:
    path = outcome_ledger_path(ledger_root)
    rows: list[dict[str, Any]] = []
    audit_warnings: list[dict[str, Any]] = []
    audit_warning_counts: Counter[str] = Counter()
    malformed_rows = 0
    unknown_outcome_rows = 0
    synthetic_rows = 0
    missing_memory_id_rows = 0
    unknown_route_rows = 0
    malformed_recorded_at_rows = 0
    weak_evidence_rows = 0
    suspicious_text_rows = 0
    raw_rows = 0
    if path.exists():
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                text = line.strip()
                if not text:
                    continue
                raw_rows += 1
                try:
                    item = json.loads(text)
                except json.JSONDecodeError:
                    malformed_rows += 1
                    audit_warning_counts["malformed_json"] += 1
                    audit_warnings.append(
                        {
                            "row_index": raw_rows,
                            "kind": "malformed_json",
                            "memory_id": "",
                            "result": "",
                            "route": "",
                            "detail": "row is not valid JSON",
                        }
                    )
                    continue
                if not isinstance(item, dict):
                    malformed_rows += 1
                    audit_warning_counts["malformed_json"] += 1
                    audit_warnings.append(
                        {
                            "row_index": raw_rows,
                            "kind": "malformed_json",
                            "memory_id": "",
                            "result": "",
                            "route": "",
                            "detail": "row is not a JSON object",
                        }
                    )
                    continue
                if row_contains_suspicious_text(item):
                    suspicious_text_rows += 1
                    audit_warning_counts["suspicious_private_text"] += 1
                    add_audit_warning(
                        audit_warnings,
                        raw_rows,
                        "suspicious_private_text",
                        item,
                        "row contains private-looking text; keep it local and do not promote blindly",
                    )
                result = str(item.get("result") or "")
                if result not in OUTCOME_RESULTS:
                    unknown_outcome_rows += 1
                    audit_warning_counts["unknown_result"] += 1
                    add_audit_warning(audit_warnings, raw_rows, "unknown_result", item, "result is not recognized")
                    continue
                if str(item.get("evidence_class") or "derived") == "synthetic":
                    synthetic_rows += 1
                    continue
                memory_id = first_line(item.get("memory_id"))
                if not memory_id:
                    missing_memory_id_rows += 1
                    audit_warning_counts["missing_memory_id"] += 1
                    add_audit_warning(audit_warnings, raw_rows, "missing_memory_id", item, "memory_id is required")
                    continue
                route, route_known = feedback_route_for(item)
                if not route_known:
                    unknown_route_rows += 1
                    audit_warning_counts["unknown_route"] += 1
                    add_audit_warning(
                        audit_warnings,
                        raw_rows,
                        "unknown_route",
                        item,
                        "route is missing or not one of the resolver route profiles",
                    )
                decay = feedback_decay_for_recorded_at(item.get("recorded_at"))
                if decay["bucket"] == "malformed_recorded_at":
                    malformed_recorded_at_rows += 1
                    audit_warning_counts["malformed_recorded_at"] += 1
                    add_audit_warning(
                        audit_warnings,
                        raw_rows,
                        "malformed_recorded_at",
                        item,
                        "recorded_at is missing or not ISO-8601 parseable",
                    )
                evidence_quality = evidence_quality_for(item)
                if evidence_quality != "strong":
                    weak_evidence_rows += 1
                    audit_warning_counts["missing_or_weak_evidence_ref"] += 1
                    add_audit_warning(
                        audit_warnings,
                        raw_rows,
                        "missing_or_weak_evidence_ref",
                        item,
                        "evidence_ref is missing, weak, or private-looking",
                    )
                row = dict(item)
                row["memory_id"] = memory_id
                row["result"] = result
                row["route"] = route
                row["route_known"] = route_known
                row["evidence_ref"] = evidence_ref_for(item)
                row["evidence_quality"] = evidence_quality
                row["decay"] = decay
                row["task_fingerprint"] = task_fingerprint(item.get("task"))
                rows.append(row)
    ignored_rows = raw_rows - len(rows)
    return {
        "path": path.as_posix(),
        "raw_rows": raw_rows,
        "accepted_rows": len(rows),
        "ignored_rows": ignored_rows,
        "synthetic_rows_ignored": synthetic_rows,
        "malformed_rows_ignored": malformed_rows,
        "unknown_outcome_rows_ignored": unknown_outcome_rows,
        "missing_memory_id_rows_ignored": missing_memory_id_rows,
        "unknown_route_rows": unknown_route_rows,
        "malformed_recorded_at_rows": malformed_recorded_at_rows,
        "missing_or_weak_evidence_ref_rows": weak_evidence_rows,
        "suspicious_private_text_rows": suspicious_text_rows,
        "audit_warning_count": len(audit_warnings),
        "audit_warning_counts": dict(audit_warning_counts),
        "audit_warnings": audit_warnings[:40],
        "rows": rows,
    }


def filtered_outcome_rows(ledger_root: Path) -> tuple[list[dict[str, Any]], int]:
    inspection = inspect_outcome_ledger(ledger_root)
    rows = inspection.get("rows") if isinstance(inspection.get("rows"), list) else []
    return rows, int(inspection.get("raw_rows", len(rows)) or 0)


def feedback_scoring_integrity() -> dict[str, Any]:
    useful_row = {"memory_id": "useful_cap", "score": {"resonance_score": 10.0}}
    stale_row = {"memory_id": "stale_cap", "score": {"resonance_score": 10.0}}
    harmful_row = {"memory_id": "harmful_cap", "score": {"resonance_score": 10.0}}
    irrelevant_once_row = {"memory_id": "irrelevant_once", "score": {"resonance_score": 10.0}}
    irrelevant_repeat_row = {"memory_id": "irrelevant_repeat", "score": {"resonance_score": 10.0}}
    missed_row = {"memory_id": "missed_review_only", "score": {"resonance_score": 10.0}}
    apply_feedback_adjustment(useful_row, {"useful_cap": {**empty_feedback_counts(), "useful": 10}})
    apply_feedback_adjustment(stale_row, {"stale_cap": {**empty_feedback_counts(), "stale": 10}})
    apply_feedback_adjustment(harmful_row, {"harmful_cap": {**empty_feedback_counts(), "harmful": 10}})
    apply_feedback_adjustment(irrelevant_once_row, {"irrelevant_once": {**empty_feedback_counts(), "irrelevant": 1}})
    apply_feedback_adjustment(irrelevant_repeat_row, {"irrelevant_repeat": {**empty_feedback_counts(), "irrelevant": 3}})
    apply_feedback_adjustment(missed_row, {"missed_review_only": {**empty_feedback_counts(), "missed": 4}})
    rows = [useful_row, stale_row, harmful_row, irrelevant_once_row, irrelevant_repeat_row, missed_row]
    metadata_fields = {
        "base_resonance_score",
        "feedback_adjustment",
        "final_resonance_score",
        "feedback_counts",
        "feedback_reason",
        "feedback_scope",
        "route_feedback_counts",
        "global_feedback_counts",
        "feedback_decay_summary",
        "evidence_quality_summary",
    }
    useful_adjustment = float(useful_row["score"].get("feedback_adjustment", 0.0) or 0.0)
    stale_adjustment = float(stale_row["score"].get("feedback_adjustment", 0.0) or 0.0)
    harmful_adjustment = float(harmful_row["score"].get("feedback_adjustment", 0.0) or 0.0)
    irrelevant_once_adjustment = float(irrelevant_once_row["score"].get("feedback_adjustment", 0.0) or 0.0)
    irrelevant_repeat_adjustment = float(irrelevant_repeat_row["score"].get("feedback_adjustment", 0.0) or 0.0)
    missed_adjustment = float(missed_row["score"].get("feedback_adjustment", 0.0) or 0.0)
    return {
        "caps_verified": useful_adjustment <= MAX_POSITIVE_FEEDBACK_ADJUSTMENT
        and stale_adjustment >= MAX_NEGATIVE_FEEDBACK_ADJUSTMENT
        and harmful_adjustment >= MAX_NEGATIVE_FEEDBACK_ADJUSTMENT,
        "useful_cap_adjustment": useful_adjustment,
        "stale_cap_adjustment": stale_adjustment,
        "harmful_cap_adjustment": harmful_adjustment,
        "harmful_cools_more_strongly_than_stale": harmful_adjustment < stale_adjustment,
        "one_irrelevant_not_over_penalized": irrelevant_once_adjustment == 0.0,
        "repeated_irrelevant_cools": irrelevant_repeat_adjustment < 0.0,
        "missed_review_only_verified": missed_adjustment == 0.0,
        "score_metadata_visible": all(metadata_fields <= set(row.get("score", {})) for row in rows),
    }


def feedback_loop_risks(ledger: dict[str, Any], feedback_counts: dict[str, dict[str, int]]) -> list[dict[str, str]]:
    risks: list[dict[str, str]] = []
    rows = ledger.get("rows") if isinstance(ledger.get("rows"), list) else []
    if int(ledger.get("accepted_rows", 0) or 0) < 3:
        risks.append(
            {
                "kind": "feedback_sparsity",
                "summary": "fewer than three accepted retrieval-outcome rows are available",
                "recommended_action": "record waking retrieval outcomes after real task use before promoting memory changes",
            }
        )
    for memory_id, counts in sorted(feedback_counts.items()):
        useful = int(counts.get("useful", 0) or 0)
        stale = int(counts.get("stale", 0) or 0)
        harmful = int(counts.get("harmful", 0) or 0)
        if useful and (stale or harmful):
            risks.append(
                {
                    "kind": "conflicting_feedback",
                    "summary": f"{memory_id} has both useful and stale/harmful feedback",
                    "recommended_action": "review route-specific use before suppressing or promoting this memory",
                }
            )
        if stale + harmful >= 2:
            risks.append(
                {
                    "kind": "stale_label_drift",
                    "summary": f"{memory_id} has repeated negative labels across local retrieval feedback",
                    "recommended_action": "verify whether this is route-specific drift or a global suppression candidate",
                }
            )
    if int(ledger.get("audit_warning_count", 0) or 0):
        risks.append(
            {
                "kind": "ledger_audit_warnings",
                "summary": f"{ledger.get('audit_warning_count')} retrieval-outcome ledger audit warning(s) need review",
                "recommended_action": "inspect review --json audit_warnings before using feedback for promotion or cooling",
            }
        )
    for item in rows:
        if item.get("result") in {"stale", "harmful"} and str(item.get("evidence_quality") or evidence_quality_for(item)) != "strong":
            risks.append(
                {
                    "kind": "negative_feedback_needs_evidence",
                    "summary": f"{item.get('memory_id')} has {item.get('result')} feedback without a strong evidence_ref",
                    "recommended_action": "add a short waking evidence_ref before using it as suppression pressure",
                }
            )
            break
    return risks[:8]


def feedback_loop_health_payload(
    *,
    run_id: str = "",
    ledger_root: Path | None = None,
) -> dict[str, Any]:
    root = ledger_root or LOCAL_CONTEXT_ROOT
    generated_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
    ledger = inspect_outcome_ledger(root)
    rows = ledger.get("rows") if isinstance(ledger.get("rows"), list) else []
    feedback_counts = aggregate_feedback_counts(rows)
    scoring = feedback_scoring_integrity()
    task = "SDCC path transfer map update without remote mutation"
    route_payload = build_context_resonance_payload(task, max_active=3, max_nudges=3, mode="waking", ledger_root=root)
    conscious_ids = [
        first_line(item.get("memory_id"))
        for item in route_payload.get("conscious_context") or []
        if isinstance(item, dict)
    ]
    route_policy_boundary_verified = any(memory_id.startswith("route_policy.sdcc_path_transfer") for memory_id in conscious_ids)
    synthetic_exclusion_verified = (
        int(ledger.get("synthetic_rows_ignored", 0) or 0) >= 0
        and all(item.get("evidence_class") != "synthetic" for item in route_payload.get("conscious_context") or [])
    )
    homeostasis = build_homeostasis_proposals(feedback_counts, rows)
    risks = feedback_loop_risks(ledger, feedback_counts)
    canary = canary_payload()
    scoring_integrity = {
        "caps_verified": bool(scoring.get("caps_verified")),
        "missed_review_only_verified": bool(scoring.get("missed_review_only_verified")),
        "synthetic_exclusion_verified": bool(synthetic_exclusion_verified),
        "route_policy_boundary_verified": bool(route_policy_boundary_verified),
        "feedback_metadata_visible": bool(scoring.get("score_metadata_visible")),
        "one_irrelevant_not_over_penalized": bool(scoring.get("one_irrelevant_not_over_penalized")),
        "repeated_irrelevant_cools": bool(scoring.get("repeated_irrelevant_cools")),
        "harmful_cools_more_strongly_than_stale": bool(scoring.get("harmful_cools_more_strongly_than_stale")),
        "canary_status": canary.get("status"),
    }
    failed_integrity = [key for key, value in scoring_integrity.items() if key != "canary_status" and value is not True]
    status = "healthy"
    if failed_integrity or canary.get("status") != "pass":
        status = "failed"
    elif risks:
        status = "warning"
    return {
        "lane_id": "context_resonance",
        "run_id": run_id,
        "generated_at": generated_at,
        "feedback_loop_status": status,
        "ledger": {key: value for key, value in ledger.items() if key != "rows"},
        "scoring_integrity": scoring_integrity,
        "homeostasis": {
            "proposal_count": len(homeostasis),
            "top_actions": homeostasis[:5],
            "tracked_registry_mutation": False,
        },
        "risks": risks,
        "auto_applied": {
            "action": "health_index_refresh",
            "changed_paths": [],
            "why_safe": [
                "local-only",
                "no tracked registry edit",
                "no science/task/external mutation",
                "synthetic rows excluded",
            ],
        },
        "deferred_for_waking": [
            {
                "action": "add_evidence_gate",
                "memory_id": first_line(risk.get("summary")).split(" ", 1)[0],
                "why": risk.get("summary"),
                "approval_needed": "waking Codex validation",
            }
            for risk in risks[:3]
        ],
        "canary_failures": canary.get("failures") or [],
    }


def render_feedback_loop_health_markdown(payload: dict[str, Any]) -> str:
    ledger = payload.get("ledger") if isinstance(payload.get("ledger"), dict) else {}
    scoring = payload.get("scoring_integrity") if isinstance(payload.get("scoring_integrity"), dict) else {}
    homeostasis = payload.get("homeostasis") if isinstance(payload.get("homeostasis"), dict) else {}
    lines = ["# Context Resonance Feedback Loop Health", "", WAKING_HEADER, ""]
    lines.append(f"- status: `{payload.get('feedback_loop_status')}`")
    lines.append(f"- ledger: `{ledger.get('path')}`")
    lines.append(
        "- rows: "
        f"raw={ledger.get('raw_rows', 0)}, accepted={ledger.get('accepted_rows', 0)}, "
        f"ignored={ledger.get('ignored_rows', 0)}, synthetic_ignored={ledger.get('synthetic_rows_ignored', 0)}, "
        f"malformed_ignored={ledger.get('malformed_rows_ignored', 0)}"
    )
    lines.append(
        "- integrity: "
        f"caps={scoring.get('caps_verified')}, missed_review_only={scoring.get('missed_review_only_verified')}, "
        f"synthetic_exclusion={scoring.get('synthetic_exclusion_verified')}, "
        f"route_policy_boundary={scoring.get('route_policy_boundary_verified')}"
    )
    lines.append(f"- homeostasis proposals: {homeostasis.get('proposal_count', 0)}")
    lines.append("- tracked registry mutation: false")
    lines.append("")
    lines.append("## Risks")
    risks = payload.get("risks") if isinstance(payload.get("risks"), list) else []
    if not risks:
        lines.append("- none")
    for risk in risks[:5]:
        if isinstance(risk, dict):
            lines.append(f"- `{risk.get('kind')}`: {risk.get('summary')} next: {risk.get('recommended_action')}")
    return "\n".join(lines) + "\n"


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


def aggregate_feedback_details(outcomes: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    feedback: dict[str, dict[str, Any]] = {}
    for item in outcomes:
        memory_id = first_line(item.get("memory_id"))
        result = str(item.get("result") or "")
        if not memory_id or result not in OUTCOME_RESULTS:
            continue
        detail = feedback.setdefault(memory_id, feedback_detail_empty())
        detail["global_counts"][result] += 1

        route = first_line(item.get("route")) or "unknown"
        route_known = bool(item.get("route_known")) and route in ROUTE_PROFILES
        route_counts = detail["route_counts"].setdefault(route, empty_feedback_counts())
        route_counts[result] += 1

        decay = item.get("decay") if isinstance(item.get("decay"), dict) else feedback_decay_for_recorded_at(item.get("recorded_at"))
        decay_bucket = str(decay.get("bucket") or "malformed_recorded_at")
        detail["decay_buckets"][decay_bucket] = int(detail["decay_buckets"].get(decay_bucket, 0) or 0) + 1

        evidence_quality = str(item.get("evidence_quality") or evidence_quality_for(item))
        detail["evidence_quality_counts"][evidence_quality] = int(detail["evidence_quality_counts"].get(evidence_quality, 0) or 0) + 1
        if result in {"harmful", "stale"} and evidence_quality != "strong":
            detail["weak_evidence_negative_count"] += 1

        evidence_ref = evidence_ref_for(item)
        if evidence_quality == "strong" and evidence_ref:
            detail["strong_evidence_refs"].append(evidence_ref)
        fingerprint = first_line(item.get("task_fingerprint")) or task_fingerprint(item.get("task"))
        if fingerprint:
            detail["task_fingerprints"].append(fingerprint)
        if route_known:
            detail["routes"].append(route)
        day = first_line(item.get("recorded_at"))[:10]
        if re.fullmatch(r"\d{4}-\d{2}-\d{2}", day):
            detail["days"].append(day)

        row_weight = float(decay.get("weight", 0.0) or 0.0) * evidence_multiplier(result, evidence_quality)
        if not route_known or row_weight <= 0:
            detail["review_only_count"] += 1
            continue
        detail["global_weighted_counts"][result] += row_weight
        route_weighted = detail["route_weighted_counts"].setdefault(route, empty_weighted_feedback_counts())
        route_weighted[result] += row_weight

    for detail in feedback.values():
        detail["global_weighted_counts"] = round_weighted_counts(detail["global_weighted_counts"])
        detail["strong_evidence_refs"] = sorted(set(detail["strong_evidence_refs"]))
        detail["task_fingerprints"] = sorted(set(detail["task_fingerprints"]))
        detail["routes"] = sorted(set(detail["routes"]))
        detail["days"] = sorted(set(detail["days"]))
        for route, counts in list(detail["route_weighted_counts"].items()):
            detail["route_weighted_counts"][route] = round_weighted_counts(counts)
    return feedback


def load_feedback_counts(ledger_root: Path | None = None) -> dict[str, dict[str, int]]:
    root = ledger_root or LOCAL_CONTEXT_ROOT
    outcomes, _ = filtered_outcome_rows(root)
    return aggregate_feedback_counts(outcomes)


def load_feedback_details(ledger_root: Path | None = None) -> dict[str, dict[str, Any]]:
    root = ledger_root or LOCAL_CONTEXT_ROOT
    outcomes, _ = filtered_outcome_rows(root)
    return aggregate_feedback_details(outcomes)


def feedback_counts_for_scope(detail: dict[str, Any], route: str | None) -> tuple[dict[str, float], dict[str, int], str]:
    if route is None:
        return round_weighted_counts(detail.get("global_weighted_counts") or {}), dict(detail.get("global_counts") or empty_feedback_counts()), "global_direct"
    route_weighted = dict(empty_weighted_feedback_counts())
    route_counts = dict(empty_feedback_counts())
    if route in detail.get("route_weighted_counts", {}):
        route_weighted.update(detail["route_weighted_counts"][route])
    if route in detail.get("route_counts", {}):
        route_counts.update(detail["route_counts"][route])
    cross_weighted = dict(empty_weighted_feedback_counts())
    for other_route, counts in (detail.get("route_weighted_counts") or {}).items():
        if other_route == route or other_route not in ROUTE_PROFILES:
            continue
        for result in FEEDBACK_RESULTS:
            cross_weighted[result] += float(counts.get(result, 0.0) or 0.0) * CROSS_ROUTE_FEEDBACK_WEIGHT
    effective = {
        result: round(float(route_weighted.get(result, 0.0) or 0.0) + float(cross_weighted.get(result, 0.0) or 0.0), 3)
        for result in FEEDBACK_RESULTS
    }
    if any(route_weighted.values()) and any(cross_weighted.values()):
        scope = "same_route_primary_with_cross_route_fallback"
    elif any(route_weighted.values()):
        scope = "same_route"
    elif any(cross_weighted.values()):
        scope = "cross_route_fallback"
    elif int(detail.get("review_only_count", 0) or 0):
        scope = "review_only_feedback"
    else:
        scope = "no_feedback"
    return effective, route_counts, scope


def apply_feedback_adjustment(
    row: dict[str, Any],
    feedback_counts: dict[str, dict[str, Any]],
    route: str | None = None,
) -> dict[str, Any]:
    score = row.get("score") if isinstance(row.get("score"), dict) else {}
    memory_id = first_line(row.get("memory_id") or row.get("id"))
    detail = normalize_feedback_detail(feedback_counts.get(memory_id) if memory_id else None)
    global_counts = dict(detail.get("global_counts") or empty_feedback_counts())
    effective_counts, route_counts, feedback_scope = feedback_counts_for_scope(detail, route)
    trainable = bool(row.get("feedback_trainable")) if "feedback_trainable" in row else infer_feedback_trainable(row)
    row["feedback_trainable"] = trainable

    base_score = float(score.get("resonance_score", 0.0) or 0.0)
    if trainable:
        adjustment, reason = bounded_feedback_adjustment(effective_counts)
    else:
        adjustment, reason = 0.0, "non-trainable candidate; local feedback ignored"
        feedback_scope = "not_trainable"
    final_score = round(base_score + adjustment, 3)
    score["base_resonance_score"] = round(base_score, 3)
    score["feedback_adjustment"] = adjustment
    score["final_resonance_score"] = final_score
    score["resonance_score"] = final_score
    score["feedback_counts"] = global_counts
    score["route_feedback_counts"] = route_counts
    score["global_feedback_counts"] = global_counts
    score["feedback_scope"] = feedback_scope
    score["feedback_decay_summary"] = {
        "decay_buckets": dict(detail.get("decay_buckets") or {}),
        "route": route,
        "effective_weighted_counts": effective_counts,
        "global_weighted_counts": round_weighted_counts(detail.get("global_weighted_counts") or {}),
        "review_only_count": int(detail.get("review_only_count", 0) or 0),
        "cross_route_fallback_weight": CROSS_ROUTE_FEEDBACK_WEIGHT,
    }
    score["evidence_quality_summary"] = {
        "quality_counts": dict(detail.get("evidence_quality_counts") or {}),
        "weak_evidence_negative_count": int(detail.get("weak_evidence_negative_count", 0) or 0),
        "strong_evidence_ref_count": len(detail.get("strong_evidence_refs") or []),
    }
    score["feedback_reason"] = reason
    row["score"] = score
    return row


def apply_feedback_to_rows(rows: list[dict[str, Any]], feedback_counts: dict[str, dict[str, Any]], route: str | None = None) -> list[dict[str, Any]]:
    for row in rows:
        if isinstance(row.get("score"), dict):
            apply_feedback_adjustment(row, feedback_counts, route)
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
                "feedback_trainable": False,
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
            "feedback_trainable": False,
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
                "feedback_trainable": False,
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
                "feedback_trainable": evidence_class != "synthetic" and retrieval_policy not in {"suppress", "quarantine_candidate"},
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
                "feedback_trainable": str(raw.get("evidence_class") or "derived") != "synthetic",
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
                "feedback_trainable": str(raw.get("evidence_class") or "human_approved") != "synthetic",
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


def suppressed_context(route: str, task_tokens: set[str], feedback_counts: dict[str, dict[str, Any]] | None = None) -> list[dict[str, Any]]:
    counts = recent_local_surface_counts()
    rows = [
        {
            "memory_id": "suppressed.local_dreams",
            "source": "agent_context/local/dreams/",
            "reason": f"{counts['dream_run_dirs']} local dream run dirs are synthetic/proposal-only and should not enter factual context by default",
            "relation_type": "warning_only",
            "retrieval_policy": "suppress",
            "evidence_class": "synthetic",
            "feedback_trainable": False,
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
            "feedback_trainable": False,
            "allowed_use": "source leads and critique only",
            "required_waking_check": "verify locally before policy, code, task, or science promotion",
        },
    ]
    suppressed_rows = registry_records(route, task_tokens)
    apply_feedback_to_rows(suppressed_rows, feedback_counts or {}, route)
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
                    f"--memory-id {memory_id!r} --task <task> --route {route} --result useful "
                    "--evidence-ref <local evidence pointer> --note <short note>"
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
    feedback_counts = load_feedback_details(ledger_root)

    route_rows = route_policy_candidates(route)
    route_rows.extend(live_workstream_hint(task_tokens, limit=1))
    route_rows.extend(artifact_hint(task_tokens))

    registry_rows = registry_records(route, task_tokens)
    schema_rows = schema_candidates(route, task_tokens)
    all_positive = route_rows + registry_rows + schema_rows
    apply_feedback_to_rows(all_positive, feedback_counts, route)
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
    apply_feedback_to_rows(negative_rows, feedback_counts, route)
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
    route = args.route or classify_task(args.task)["route"]
    row = {
        "recorded_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "memory_id": args.memory_id,
        "task": args.task[:500],
        "route": route,
        "result": args.result,
        "note": args.note[:500],
        "evidence_ref": (args.evidence_ref or "")[:500],
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


def homeostasis_support_summary(
    memory_id: str,
    outcomes: list[dict[str, Any]],
    relevant_results: set[str],
    counts: dict[str, int],
) -> dict[str, Any]:
    relevant_rows = [
        item
        for item in outcomes
        if first_line(item.get("memory_id")) == memory_id and str(item.get("result") or "") in relevant_results
    ]
    result_counts = {result: int(counts.get(result, 0) or 0) for result in relevant_results}
    task_fingerprints = {
        first_line(item.get("task_fingerprint")) or task_fingerprint(item.get("task"))
        for item in relevant_rows
        if first_line(item.get("task")) or first_line(item.get("task_fingerprint"))
    }
    routes = {
        first_line(item.get("route"))
        for item in relevant_rows
        if first_line(item.get("route")) in ROUTE_PROFILES
    }
    evidence_refs = {
        evidence_ref_for(item)
        for item in relevant_rows
        if str(item.get("evidence_quality") or evidence_quality_for(item)) == "strong" and evidence_ref_for(item)
    }
    days = {
        first_line(item.get("recorded_at"))[:10]
        for item in relevant_rows
        if re.fullmatch(r"\d{4}-\d{2}-\d{2}", first_line(item.get("recorded_at"))[:10])
    }
    total_relevant = sum(result_counts.values())
    max_relevant_result_count = max(result_counts.values()) if result_counts else 0
    criteria = {
        "repeated_outcomes": total_relevant >= 2 or max_relevant_result_count >= 2,
        "distinct_task_fingerprints": len(task_fingerprints) >= 2,
        "distinct_routes": len(routes) >= 2,
        "distinct_evidence_refs": len(evidence_refs) >= 2,
        "distinct_days": len(days) >= 2,
    }
    return {
        "relevant_results": sorted(relevant_results),
        "result_counts": result_counts,
        "criteria": criteria,
        "criteria_met_count": sum(1 for value in criteria.values() if value),
        "distinct_task_fingerprints": len(task_fingerprints),
        "distinct_routes": len(routes),
        "distinct_evidence_refs": len(evidence_refs),
        "distinct_days": len(days),
        "strong_evidence_ref_count": len(evidence_refs),
        "relevant_outcome_count": total_relevant,
    }


def homeostasis_action_for_counts(
    counts: dict[str, int],
    support: dict[str, Any] | None = None,
) -> tuple[str, str] | None:
    useful = int(counts.get("useful", 0) or 0)
    stale = int(counts.get("stale", 0) or 0)
    harmful = int(counts.get("harmful", 0) or 0)
    irrelevant = int(counts.get("irrelevant", 0) or 0)
    missed = int(counts.get("missed", 0) or 0)
    criteria_met = int((support or {}).get("criteria_met_count", 0) or 0)
    strong_evidence_refs = int((support or {}).get("strong_evidence_ref_count", 0) or 0)

    def strong_enough(*, negative: bool = False) -> bool:
        if criteria_met < 2:
            return False
        if negative and strong_evidence_refs < 1:
            return False
        return True

    if harmful >= 2 or (harmful and stale):
        if strong_enough(negative=True):
            return "suppress_candidate", "repeated harmful/stale retrieval feedback with independent support"
        return "cue_review", "negative feedback needs repeated, distinct, and evidence-backed support before suppression"
    if harmful >= 1:
        if strong_enough(negative=True):
            return "quarantine_candidate", "harmful retrieval feedback with independent support"
        return "cue_review", "harmful feedback needs review; no automatic quarantine pressure"
    if stale >= 1 or irrelevant >= 2:
        if strong_enough(negative=True):
            return "cool", "stale or repeatedly irrelevant retrieval feedback with independent support"
        return "cue_review", "cooling signal is not yet repeated or evidence-backed enough"
    if missed >= 2:
        if strong_enough():
            return "promote", "repeated missed-memory feedback with independent support; review cue terms or tracked record"
        return "cue_review", "missed-memory feedback needs independent support before promotion"
    if missed >= 1:
        return "cue_review", "missed-memory feedback; no automatic ranking boost"
    if useful >= 2:
        if strong_enough():
            return "keep_hot", "repeated useful retrieval feedback with independent support"
        return "keep_latent", "useful feedback exists but is not independently supported enough to keep hot"
    if useful >= 1:
        return "keep_latent", "single useful retrieval feedback"
    return None


def homeostasis_relevant_results(counts: dict[str, int]) -> set[str]:
    if int(counts.get("harmful", 0) or 0) or int(counts.get("stale", 0) or 0):
        return {"harmful", "stale"}
    if int(counts.get("irrelevant", 0) or 0):
        return {"irrelevant"}
    if int(counts.get("missed", 0) or 0):
        return {"missed"}
    if int(counts.get("useful", 0) or 0):
        return {"useful"}
    return set()


def build_homeostasis_proposals(
    feedback_counts: dict[str, dict[str, int]],
    outcomes: list[dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    outcome_rows = outcomes or []
    for memory_id in sorted(feedback_counts):
        counts = feedback_counts[memory_id]
        relevant_results = homeostasis_relevant_results(counts)
        support = homeostasis_support_summary(memory_id, outcome_rows, relevant_results, counts)
        action_reason = homeostasis_action_for_counts(counts, support)
        if not action_reason:
            continue
        action, reason = action_reason
        rows.append(
            {
                "memory_id": memory_id,
                "action": action,
                "reason": reason,
                "feedback_counts": dict(counts),
                "support_summary": support,
                "threshold_policy": (
                    "strong promote/cool/quarantine/suppress proposals require at least two independent support "
                    "signals; negative strong actions also need a strong evidence_ref"
                ),
                "requires_waking_validation": True,
                "mutation_boundary": "proposal_only",
                "allowed_use": "routing homeostasis proposal only; do not edit tracked memory automatically",
            }
        )
    return rows


def local_maintenance_review(ledger_root: Path | None = None) -> dict[str, Any]:
    root = ledger_root or LOCAL_CONTEXT_ROOT
    ledger_health = inspect_outcome_ledger(root)
    outcomes = ledger_health.get("rows") if isinstance(ledger_health.get("rows"), list) else []
    raw_outcome_count = int(ledger_health.get("raw_rows", len(outcomes)) or 0)
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
        "ledger_health": {key: value for key, value in ledger_health.items() if key != "rows"},
        "interference_count": len(interference),
        "salience_record_count": len(salience_records),
        "feedback_counts": feedback_counts,
        "homeostasis_proposals": build_homeostasis_proposals(feedback_counts, outcomes),
        "pollution_memory_counts": dict(pollution_counts.most_common(8)),
        "useful_memory_counts": dict(useful_trap_counts.most_common(8)),
        "duplicate_source_counts": duplicate_sources,
        "route_index_refresh_needed": not bool(salience_records),
        "candidate_evidence": evidence,
        "maintenance_policy": "local-only distilled routing evidence; synthetic rows cannot justify auto-apply",
    }


def review_payload(ledger_root: Path) -> dict[str, Any]:
    ledger_health = inspect_outcome_ledger(ledger_root)
    outcomes = ledger_health.get("rows") if isinstance(ledger_health.get("rows"), list) else []
    raw_outcome_count = int(ledger_health.get("raw_rows", len(outcomes)) or 0)
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
        "ledger_health": {key: value for key, value in ledger_health.items() if key != "rows"},
        "audit_warning_count": int(ledger_health.get("audit_warning_count", 0) or 0),
        "audit_warning_counts": dict(ledger_health.get("audit_warning_counts") or {}),
        "audit_warnings": ledger_health.get("audit_warnings") or [],
        "interference_count": len(interference),
        "result_counts": dict(result_counts),
        "feedback_counts": feedback_counts,
        "top_memory_ids": dict(memory_counts.most_common(8)),
        "cooldown_candidates": harmful[-8:],
        "homeostasis_proposals": build_homeostasis_proposals(feedback_counts, outcomes),
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
                        "feedback_scope",
                        "route_feedback_counts",
                        "global_feedback_counts",
                        "feedback_decay_summary",
                        "evidence_quality_summary",
                    ):
                        if field not in score:
                            failures.append(f"{spec['id']}: {item.get('memory_id')} lacks score.{field}")
                    if "feedback_trainable" not in item:
                        failures.append(f"{spec['id']}: {item.get('memory_id')} lacks feedback_trainable")
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
            task = row.get("task", "context resonance canary")
            payload = {
                "recorded_at": row.get("recorded_at", datetime.now(timezone.utc).isoformat(timespec="seconds")),
                "memory_id": row["memory_id"],
                "task": task,
                "route": row.get("route", classify_task(str(task))["route"]),
                "result": row["result"],
                "note": row.get("note", "canary fixture"),
                "evidence_ref": row.get("evidence_ref", f"canary/{row['memory_id']}/{row['result']}"),
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
    route = "os_maintenance"

    with TemporaryDirectory(prefix="context_resonance_feedback_boost_") as temp_dir:
        root = Path(temp_dir)
        base_payload = build_context_resonance_payload(task, max_active=3, max_nudges=3, mode="waking", ledger_root=root)
        base_position = section_position(base_payload, "latent_context_nudges", "context_resonance_schema")
        write_canary_outcomes(
            root,
            [
                {
                    "memory_id": "context_resonance_schema",
                    "result": "useful",
                    "route": route,
                    "task": task,
                    "evidence_ref": "canary/context-resonance-schema/useful-1",
                },
                {
                    "memory_id": "context_resonance_schema",
                    "result": "useful",
                    "route": route,
                    "task": "context resonance feedback stabilization for OS memory",
                    "evidence_ref": "canary/context-resonance-schema/useful-2",
                },
                {
                    "memory_id": "route_policy.os_maintenance.1",
                    "result": "useful",
                    "route": route,
                    "task": task,
                    "evidence_ref": "canary/route-policy/nontrainable",
                },
            ],
        )
        boosted_payload = build_context_resonance_payload(task, max_active=3, max_nudges=3, mode="waking", ledger_root=root)
        boosted_position = section_position(boosted_payload, "latent_context_nudges", "context_resonance_schema")
        boosted_item = find_payload_item(boosted_payload, "latent_context_nudges", "context_resonance_schema")
        boosted_score = boosted_item.get("score") if isinstance(boosted_item.get("score"), dict) else {}
        policy_item = find_payload_item(boosted_payload, "conscious_context", "route_policy.os_maintenance.1")
        policy_score = policy_item.get("score") if isinstance(policy_item.get("score"), dict) else {}
        if base_position is None or boosted_position is None or boosted_position >= base_position:
            failures.append("feedback_useful_boost: useful memory did not rise within latent_context_nudges")
        if float(boosted_score.get("feedback_adjustment", 0.0) or 0.0) <= 0:
            failures.append("feedback_useful_boost: useful memory did not receive positive feedback adjustment")
        if not boosted_item.get("feedback_trainable"):
            failures.append("feedback_useful_boost: durable schema row is not marked trainable")
        if policy_item.get("feedback_trainable"):
            failures.append("feedback_nontrainable_route_policy: route policy row is marked trainable")
        if float(policy_score.get("feedback_adjustment", 0.0) or 0.0) != 0.0:
            failures.append("feedback_nontrainable_route_policy: route policy row was feedback-adjusted")
        cases.append(
            {
                "id": "feedback_trainable_useful_boost",
                "memory_id": "context_resonance_schema",
                "base_position": base_position,
                "boosted_position": boosted_position,
                "feedback_adjustment": boosted_score.get("feedback_adjustment"),
                "feedback_scope": boosted_score.get("feedback_scope"),
                "route_policy_adjustment": policy_score.get("feedback_adjustment"),
                "route_policy_trainable": policy_item.get("feedback_trainable"),
            }
        )

    with TemporaryDirectory(prefix="context_resonance_route_scope_same_") as same_dir, TemporaryDirectory(
        prefix="context_resonance_route_scope_cross_"
    ) as cross_dir:
        same_root = Path(same_dir)
        cross_root = Path(cross_dir)
        write_canary_outcomes(
            same_root,
            [
                {
                    "memory_id": "context_resonance_schema",
                    "result": "useful",
                    "route": route,
                    "task": task,
                    "evidence_ref": "canary/route-scope/same",
                }
            ],
        )
        write_canary_outcomes(
            cross_root,
            [
                {
                    "memory_id": "context_resonance_schema",
                    "result": "useful",
                    "route": "plot_generation",
                    "task": "plot generation for a slide-ready PNG with provenance",
                    "evidence_ref": "canary/route-scope/cross",
                }
            ],
        )
        same_payload = build_context_resonance_payload(task, max_active=3, max_nudges=3, mode="waking", ledger_root=same_root)
        cross_payload = build_context_resonance_payload(task, max_active=3, max_nudges=3, mode="waking", ledger_root=cross_root)
        same_score = find_payload_item(same_payload, "latent_context_nudges", "context_resonance_schema").get("score", {})
        cross_score = find_payload_item(cross_payload, "latent_context_nudges", "context_resonance_schema").get("score", {})
        same_adjustment = float(same_score.get("feedback_adjustment", 0.0) or 0.0)
        cross_adjustment = float(cross_score.get("feedback_adjustment", 0.0) or 0.0)
        if not (same_adjustment > cross_adjustment > 0):
            failures.append("feedback_route_scope: same-route feedback was not stronger than cross-route fallback")
        if same_score.get("feedback_scope") != "same_route":
            failures.append("feedback_route_scope: same-route feedback scope was not explicit")
        if cross_score.get("feedback_scope") != "cross_route_fallback":
            failures.append("feedback_route_scope: cross-route fallback scope was not explicit")
        cases.append(
            {
                "id": "feedback_route_scope",
                "same_route_adjustment": same_adjustment,
                "cross_route_adjustment": cross_adjustment,
                "same_scope": same_score.get("feedback_scope"),
                "cross_scope": cross_score.get("feedback_scope"),
            }
        )

    with TemporaryDirectory(prefix="context_resonance_decay_recent_") as recent_dir, TemporaryDirectory(
        prefix="context_resonance_decay_old_"
    ) as old_dir, TemporaryDirectory(prefix="context_resonance_decay_very_old_") as very_old_dir:
        recent_root = Path(recent_dir)
        old_root = Path(old_dir)
        very_old_root = Path(very_old_dir)
        now = datetime.now(timezone.utc).replace(microsecond=0)
        old_recorded_at = (now - timedelta(days=120)).isoformat()
        very_old_recorded_at = (now - timedelta(days=240)).isoformat()
        common_row = {
            "memory_id": "context_resonance_schema",
            "result": "useful",
            "route": route,
            "task": task,
            "evidence_ref": "canary/decay/context-resonance-schema",
        }
        write_canary_outcomes(recent_root, [common_row])
        write_canary_outcomes(old_root, [{**common_row, "recorded_at": old_recorded_at}])
        write_canary_outcomes(very_old_root, [{**common_row, "recorded_at": very_old_recorded_at}])
        recent_score = find_payload_item(
            build_context_resonance_payload(task, max_active=3, max_nudges=3, mode="waking", ledger_root=recent_root),
            "latent_context_nudges",
            "context_resonance_schema",
        ).get("score", {})
        old_score = find_payload_item(
            build_context_resonance_payload(task, max_active=3, max_nudges=3, mode="waking", ledger_root=old_root),
            "latent_context_nudges",
            "context_resonance_schema",
        ).get("score", {})
        very_old_score = find_payload_item(
            build_context_resonance_payload(task, max_active=3, max_nudges=3, mode="waking", ledger_root=very_old_root),
            "latent_context_nudges",
            "context_resonance_schema",
        ).get("score", {})
        recent_adjustment = float(recent_score.get("feedback_adjustment", 0.0) or 0.0)
        old_adjustment = float(old_score.get("feedback_adjustment", 0.0) or 0.0)
        very_old_adjustment = float(very_old_score.get("feedback_adjustment", 0.0) or 0.0)
        if not (recent_adjustment > old_adjustment > very_old_adjustment):
            failures.append("feedback_recency_decay: old feedback was not weaker than recent feedback")
        if very_old_adjustment != 0.0:
            failures.append("feedback_recency_decay: very old feedback still affected scoring")
        cases.append(
            {
                "id": "feedback_recency_decay",
                "recent_adjustment": recent_adjustment,
                "old_adjustment": old_adjustment,
                "very_old_adjustment": very_old_adjustment,
                "old_decay": old_score.get("feedback_decay_summary"),
                "very_old_decay": very_old_score.get("feedback_decay_summary"),
            }
        )

    strong_negative_adjustment = 0.0
    with TemporaryDirectory(prefix="context_resonance_feedback_penalty_") as temp_dir:
        root = Path(temp_dir)
        write_canary_outcomes(
            root,
            [
                {
                    "memory_id": "context_resonance_architecture",
                    "result": "harmful",
                    "route": route,
                    "task": task,
                    "evidence_ref": "canary/context-resonance-architecture/harmful",
                },
                {
                    "memory_id": "context_resonance_architecture",
                    "result": "stale",
                    "route": route,
                    "task": "context resonance final stabilization against stale bias",
                    "evidence_ref": "canary/context-resonance-architecture/stale",
                },
            ],
        )
        payload = build_context_resonance_payload(task, max_active=3, max_nudges=3, mode="waking", ledger_root=root)
        item = find_payload_item(payload, "latent_context_nudges", "context_resonance_architecture")
        score = item.get("score") if isinstance(item.get("score"), dict) else {}
        strong_negative_adjustment = float(score.get("feedback_adjustment", 0.0) or 0.0)
        review = review_payload(root)
        actions = {
            item.get("action")
            for item in review.get("homeostasis_proposals", [])
            if isinstance(item, dict) and item.get("memory_id") == "context_resonance_architecture"
        }
        if strong_negative_adjustment >= 0:
            failures.append("feedback_harmful_stale_penalty: harmful/stale memory was not penalized")
        if float(score.get("final_resonance_score", 0.0) or 0.0) >= float(score.get("base_resonance_score", 0.0) or 0.0):
            failures.append("feedback_harmful_stale_penalty: final score did not fall below base score")
        if not actions & {"suppress_candidate", "quarantine_candidate", "cool"}:
            failures.append("feedback_harmful_stale_penalty: strong negative evidence did not create cooling proposal")
        cases.append(
            {
                "id": "feedback_harmful_stale_penalty",
                "memory_id": "context_resonance_architecture",
                "feedback_adjustment": strong_negative_adjustment,
                "homeostasis_actions": sorted(str(action) for action in actions if action),
            }
        )

    with TemporaryDirectory(prefix="context_resonance_weak_negative_") as temp_dir:
        root = Path(temp_dir)
        write_canary_outcomes(
            root,
            [
                {
                    "memory_id": "context_resonance_architecture",
                    "result": "harmful",
                    "route": route,
                    "task": task,
                    "evidence_ref": "",
                },
                {
                    "memory_id": "context_resonance_architecture",
                    "result": "stale",
                    "route": route,
                    "task": "context resonance stale feedback without evidence",
                    "evidence_ref": "",
                },
            ],
        )
        payload = build_context_resonance_payload(task, max_active=3, max_nudges=3, mode="waking", ledger_root=root)
        item = find_payload_item(payload, "latent_context_nudges", "context_resonance_architecture")
        score = item.get("score") if isinstance(item.get("score"), dict) else {}
        weak_adjustment = float(score.get("feedback_adjustment", 0.0) or 0.0)
        review = review_payload(root)
        actions = {
            item.get("action")
            for item in review.get("homeostasis_proposals", [])
            if isinstance(item, dict) and item.get("memory_id") == "context_resonance_architecture"
        }
        if weak_adjustment <= -0.5:
            failures.append("feedback_weak_negative_not_strong: missing evidence produced a strong suppression-like penalty")
        if strong_negative_adjustment and weak_adjustment <= strong_negative_adjustment:
            failures.append("feedback_weak_negative_not_strong: weak negative evidence cooled as strongly as strong evidence")
        if actions & {"suppress_candidate", "quarantine_candidate", "cool"}:
            failures.append("feedback_weak_negative_not_strong: weak negative evidence produced a strong homeostasis action")
        cases.append(
            {
                "id": "feedback_weak_negative_not_strong",
                "feedback_adjustment": weak_adjustment,
                "homeostasis_actions": sorted(str(action) for action in actions if action),
                "evidence_quality": score.get("evidence_quality_summary"),
            }
        )

    irrelevant_row = {"memory_id": "irrelevant_dummy", "score": {"resonance_score": 10.0}}
    irrelevant_once_row = {"memory_id": "irrelevant_once_dummy", "score": {"resonance_score": 10.0}}
    apply_feedback_adjustment(
        irrelevant_row,
        {"irrelevant_dummy": {**empty_feedback_counts(), "irrelevant": 3}},
    )
    apply_feedback_adjustment(
        irrelevant_once_row,
        {"irrelevant_once_dummy": {**empty_feedback_counts(), "irrelevant": 1}},
    )
    irrelevant_score = irrelevant_row["score"]
    irrelevant_once_score = irrelevant_once_row["score"]
    if float(irrelevant_once_score.get("feedback_adjustment", 0.0) or 0.0) != 0.0:
        failures.append("feedback_single_irrelevant_not_over_penalized: one irrelevant row changed ranking")
    if float(irrelevant_score.get("feedback_adjustment", 0.0) or 0.0) >= 0:
        failures.append("feedback_repeated_irrelevant_cools: repeated irrelevant memory was not cooled")
    cases.append(
        {
            "id": "feedback_repeated_irrelevant_cools",
            "memory_id": "irrelevant_dummy",
            "single_irrelevant_adjustment": irrelevant_once_score.get("feedback_adjustment"),
            "feedback_adjustment": irrelevant_score.get("feedback_adjustment"),
        }
    )

    scoring = feedback_scoring_integrity()
    if not scoring.get("caps_verified"):
        failures.append("feedback_caps_hold: useful/stale/harmful caps failed")
    if not scoring.get("harmful_cools_more_strongly_than_stale"):
        failures.append("feedback_harmful_stronger_than_stale: harmful did not cool more strongly than stale")
    if not scoring.get("missed_review_only_verified"):
        failures.append("feedback_missed_review_only_caps: missed feedback changed ranking")
    cases.append(
        {
            "id": "feedback_caps_and_metadata",
            "useful_cap_adjustment": scoring.get("useful_cap_adjustment"),
            "stale_cap_adjustment": scoring.get("stale_cap_adjustment"),
            "harmful_cap_adjustment": scoring.get("harmful_cap_adjustment"),
            "score_metadata_visible": scoring.get("score_metadata_visible"),
        }
    )

    with TemporaryDirectory(prefix="context_resonance_feedback_filter_") as temp_dir:
        root = Path(temp_dir)
        ledger = outcome_ledger_path(root)
        root.mkdir(parents=True, exist_ok=True)
        ledger.write_text(
            "\n".join(
                [
                    "{not json",
                    json.dumps({"memory_id": "context_resonance_schema", "result": "unknown", "evidence_class": "derived"}),
                    json.dumps({"memory_id": "context_resonance_schema", "result": "useful", "evidence_class": "synthetic"}),
                    json.dumps({"result": "useful", "evidence_class": "derived"}),
                    json.dumps(
                        {
                            "memory_id": "context_resonance_schema",
                            "result": "stale",
                            "route": "bogus_route",
                            "recorded_at": "not-a-date",
                            "evidence_ref": "n/a",
                            "note": "password token fixture",
                            "evidence_class": "derived",
                        }
                    ),
                    json.dumps(
                        {
                            "memory_id": "context_resonance_schema",
                            "result": "useful",
                            "route": route,
                            "recorded_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                            "evidence_ref": "canary/filtering/clean-useful",
                            "evidence_class": "derived",
                        }
                    ),
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        inspection = inspect_outcome_ledger(root)
        filtered_payload = build_context_resonance_payload(task, max_active=3, max_nudges=3, mode="waking", ledger_root=root)
        filtered_item = find_payload_item(filtered_payload, "latent_context_nudges", "context_resonance_schema")
        filtered_score = filtered_item.get("score") if isinstance(filtered_item.get("score"), dict) else {}
        warning_kinds = {item.get("kind") for item in inspection.get("audit_warnings") or [] if isinstance(item, dict)}
        expected_warning_kinds = {
            "malformed_json",
            "unknown_result",
            "missing_memory_id",
            "unknown_route",
            "malformed_recorded_at",
            "missing_or_weak_evidence_ref",
            "suspicious_private_text",
        }
        if inspection.get("accepted_rows") != 2:
            failures.append("feedback_ledger_filtering: expected exactly two accepted non-synthetic rows")
        if inspection.get("synthetic_rows_ignored") != 1 or inspection.get("malformed_rows_ignored") != 1:
            failures.append("feedback_ledger_filtering: malformed or synthetic row accounting failed")
        if not expected_warning_kinds <= warning_kinds:
            failures.append("feedback_ledger_filtering: audit warning coverage is incomplete")
        if float(filtered_score.get("feedback_adjustment", 0.0) or 0.0) != 0.25:
            failures.append("feedback_ledger_filtering: ignored rows affected scoring")
        cases.append(
            {
                "id": "feedback_ledger_filtering",
                "accepted_rows": inspection.get("accepted_rows"),
                "ignored_rows": inspection.get("ignored_rows"),
                "synthetic_rows_ignored": inspection.get("synthetic_rows_ignored"),
                "malformed_rows_ignored": inspection.get("malformed_rows_ignored"),
                "audit_warning_kinds": sorted(str(kind) for kind in warning_kinds if kind),
                "feedback_adjustment": filtered_score.get("feedback_adjustment"),
            }
        )

    with TemporaryDirectory(prefix="context_resonance_feedback_missed_") as temp_dir:
        root = Path(temp_dir)
        write_canary_outcomes(
            root,
            [
                {
                    "memory_id": "context_resonance_schema",
                    "result": "missed",
                    "route": route,
                    "task": task,
                    "evidence_ref": "canary/context-resonance-schema/missed",
                }
            ],
        )
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
            ids = case.get("memory_ids") or ([case.get("memory_id")] if case.get("memory_id") else [])
            print(f"- {case['id']}: route={case.get('route', '')}, ids={', '.join(str(item) for item in ids)}")
        for failure in payload["failures"]:
            print(f"FAIL: {failure}")
    return 1 if payload["status"] != "pass" else 0


def feedback_health(args: argparse.Namespace) -> int:
    ledger_root = Path(args.ledger_root) if args.ledger_root else LOCAL_CONTEXT_ROOT
    payload = feedback_loop_health_payload(run_id=args.run_id or "", ledger_root=ledger_root)
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(render_feedback_loop_health_markdown(payload), end="")
    return 1 if payload.get("feedback_loop_status") == "failed" else 0


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
    record_parser.add_argument("--route", choices=sorted(ROUTE_PROFILES))
    record_parser.add_argument("--result", required=True, choices=sorted(OUTCOME_RESULTS))
    record_parser.add_argument("--evidence-ref", default="")
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

    health_parser = subparsers.add_parser("feedback-health", help="check adaptive feedback-loop health")
    health_parser.add_argument("--ledger-root")
    health_parser.add_argument("--run-id", default="")
    health_parser.add_argument("--json", action="store_true")
    health_parser.set_defaults(func=feedback_health)
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
    commands = {"resolve", "record-outcome", "review", "maintenance-review", "canary", "feedback-health"}
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
