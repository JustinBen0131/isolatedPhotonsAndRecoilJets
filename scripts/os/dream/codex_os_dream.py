#!/usr/bin/env python3
"""Run private synthetic ThesisAnalysis OS dreams."""

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
import difflib
import html
import json
import os
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from codex_work_register_common import DEFAULT_REGISTER, first_line, load_register, parse_when, sorted_workstreams
from codex_thesis_radar import analyze as analyze_thesis_radar


SYNTHETIC_HEADER = "SYNTHETIC DREAM OUTPUT - NOT USER APPROVAL - NOT REAL USER INTENT"
DREAM_ROOT = Path("agent_context/local/dreams")
ARTIFACT_REGISTRY = Path("agent_context/ARTIFACT_REGISTRY.yaml")
THESIS_MAP = Path("agent_context/THESIS_NARRATIVE_MAP.md")
EVENT_LOG = Path("agent_context/local/os_events.jsonl")
CHATGPT_RESEARCH_ROOT = Path("agent_context/local/chatgpt_research")
DREAM_INDEX = DREAM_ROOT / "dream_index.jsonl"
APPROVAL_LIKE_PATTERN = re.compile(
    r"\bJustin approved\b|\breal user approval\b|\buser approved\b",
    re.IGNORECASE,
)

LIVE_STATUSES = {"active", "running", "waiting", "blocked", "review"}
ACTIVE_STATUSES = {"active", "running"}
CLAIM_ORDER = [
    "physics_target",
    "pp_baseline",
    "embedded_background",
    "ml_photon_id",
    "auau_response",
    "final_physics",
    "os_infrastructure",
]
SEARCH_ROOTS = [Path("agent_context"), Path("codex_notes"), Path("scripts")]
TEXT_SUFFIXES = {".md", ".yaml", ".yml", ".py", ".sh", ".json", ".txt", ".toml"}
EVIDENCE_CLASSES = [
    "REAL_OBSERVED",
    "REAL_DERIVED",
    "SYNTHETIC_SCENARIO",
    "SYNTHETIC_RESULT",
    "HYPOTHESIS",
    "VALIDATED_RUNBOOK",
    "APPROVED_PATCH",
]
SCENARIO_PORTFOLIO = [
    "exact_replay",
    "counterfactual_replay",
    "threat_simulation",
    "abstraction_pass",
    "interface_failure",
    "approval_boundary_test",
    "provenance_audit",
    "memory_hygiene",
    "storage_hygiene",
    "safe_overnight_rehearsal",
    "physics_hypothesis_scout",
    "literature_scout",
    "ideal_final_figure_design",
    "all_corners_review",
]
LITERATURE_SEARCH_ROOTS = [
    Path("usefulDocs"),
    Path("agent_context"),
    Path("codex_notes"),
    Path("ppg12codeGit"),
]
LITERATURE_SUFFIXES = {".pdf", ".md", ".txt", ".tex", ".bib"}
LITERATURE_KEYWORDS = (
    "gamma",
    "photon",
    "jet",
    "xj",
    "xjg",
    "xjgamma",
    "atlas",
    "sphenix",
    "ppg12",
    "ppg18",
    "ppg19",
    "isolation",
    "purity",
    "unfold",
    "closure",
)
ALL_CORNERS = [
    "status_staleness",
    "active_jobs",
    "artifact_provenance",
    "duplicate_run_safety",
    "cleanup_storage_memory",
    "prompt_archetypes",
    "literature_context",
    "ideal_final_figures",
    "physics_hypotheses",
    "slides_and_presentation",
    "task_closeout",
    "approval_boundaries",
]
IDEAL_FIGURE_SPECS = [
    {
        "id": "xjgamma_modification_target",
        "title": "Synthetic Target: Au+Au gamma-jet xJgamma modification",
        "claim": "final_physics",
        "observable": "x_Jgamma",
        "x_label": "x_Jgamma",
        "y_label": "1/N_gamma dN_jet/dx_Jgamma",
        "required_inputs": "validated pp baseline, Au+Au signal selection, response/unfolding, purity correction, centrality control",
    },
    {
        "id": "pp_embedded_closure_target",
        "title": "Synthetic Target: pp and embedded closure ladder",
        "claim": "embedded_background",
        "observable": "closure ratio",
        "x_label": "pT or x_Jgamma bin",
        "y_label": "reco / truth",
        "required_inputs": "registered pp closure, embedded inclusive stitching QA, response closure, systematic variations",
    },
    {
        "id": "photon_id_purity_efficiency_target",
        "title": "Synthetic Target: photon-ID purity-efficiency operating point",
        "claim": "ml_photon_id",
        "observable": "purity vs efficiency",
        "x_label": "photon efficiency",
        "y_label": "photon purity",
        "required_inputs": "PPG12-style references, BDT/ML validation, leakage checks, domain-shift controls",
    },
]
PROMPT_ARCHETYPES = [
    {
        "id": "status_pressure",
        "scenario_type": "exact_replay",
        "synthetic_prompt": "is it complete? did something fail again?",
        "safe_first_action": "refresh local register/evidence first, then use read-only live status only when waking policy permits",
        "expected_response_shape": "concise status, exact evidence, current blocker, next check timing",
        "quality_checks": [
            "does not call a run complete without evidence",
            "reports held/running/idle/stale counts when available",
            "does not mutate email, jobs, or files from dream mode",
        ],
        "failure_mode_prevented": "vague reassurance or stale status from memory",
    },
    {
        "id": "evidence_provenance",
        "scenario_type": "provenance_audit",
        "synthetic_prompt": "what plot proves this and where did the data come from?",
        "safe_first_action": "inspect artifact registry, generation command, source inputs, QA state, and slide target",
        "expected_response_shape": "artifact path, source inputs, command, QA status, limitations, thesis claim",
        "quality_checks": [
            "names concrete files or says provenance is missing",
            "does not upgrade candidate plots to slide-ready without QA",
            "separates real evidence from hypothesis",
        ],
        "failure_mode_prevented": "plot claims detached from source artifacts",
    },
    {
        "id": "targeted_execution",
        "scenario_type": "approval_boundary_test",
        "synthetic_prompt": "proceed, but reuse the working infrastructure and keep it narrow",
        "safe_first_action": "run duplicate-run guard and identify the smallest validated path before any execution",
        "expected_response_shape": "narrow scope, reused script/runbook, guard result, approval boundary, validation plan",
        "quality_checks": [
            "blocks materially equivalent duplicate runs",
            "does not write broad custom infrastructure for a narrow task",
            "uses approval-gated risky actions",
        ],
        "failure_mode_prevented": "overbuilt reruns and unknown broad failure surfaces",
    },
    {
        "id": "cleanup_hygiene",
        "scenario_type": "memory_hygiene",
        "synthetic_prompt": "keep the queue, files, and memory clean without harming analysis",
        "safe_first_action": "produce cleanup candidates with evidence, risk, and verification commands only",
        "expected_response_shape": "candidate scope, evidence, risk, safe check, required approval",
        "quality_checks": [
            "does not delete, move, archive, or mark memory obsolete in dream mode",
            "separates local cleanup from SDCC cleanup",
            "names approval needed before cleanup",
        ],
        "failure_mode_prevented": "silent destructive cleanup or accumulating unusable clutter",
    },
    {
        "id": "meta_os_upgrade",
        "scenario_type": "abstraction_pass",
        "synthetic_prompt": "work with ChatGPT and take the dream state to the next level of efficiency, organization, and accuracy",
        "safe_first_action": "load ChatGPT delegation policy, sanitize external prompt, synthesize locally, implement small validated OS changes",
        "expected_response_shape": "sanitized external question, concrete local changes, validation, rejected advice",
        "quality_checks": [
            "does not expose private paths, logs, emails, or secrets to ChatGPT",
            "treats ChatGPT as critique rather than authority",
            "turns useful critique into validated local code or policy",
        ],
        "failure_mode_prevented": "broad motivational OS talk without durable implementation",
    },
    {
        "id": "physics_scout",
        "scenario_type": "physics_hypothesis_scout",
        "synthetic_prompt": "try different physics scenarios if the usable data hints at something important",
        "safe_first_action": "require registered usable artifacts, null tests, systematic risks, and waking approval before any run",
        "expected_response_shape": "hypothesis, input surface, qualitative signature, null/control, systematic, minimal safe diagnostic",
        "quality_checks": [
            "does not claim discovery or result from a dream",
            "does not start production analysis",
            "requires provenance-backed data surfaces",
        ],
        "failure_mode_prevented": "hype-driven physics overclaim or curiosity-driven broad scans",
    },
    {
        "id": "presentation_artifact",
        "scenario_type": "provenance_audit",
        "synthetic_prompt": "make the same slide or plot, but surgically include the new sample and show the PNG first",
        "safe_first_action": "locate the original generation path and regenerate a local candidate with the smallest delta",
        "expected_response_shape": "candidate PNG path, source code path, changed inputs, visual QA, no slide mutation",
        "quality_checks": [
            "does not modify Google Slides before candidate approval",
            "matches original binning/style unless intentionally changed",
            "shows image and path",
        ],
        "failure_mode_prevented": "wrong-slide regeneration or unapproved deck edits",
    },
    {
        "id": "frustration_debug",
        "scenario_type": "interface_failure",
        "synthetic_prompt": "why did this fail again and why is it so heavy?",
        "safe_first_action": "compare failing command/config to last known working infrastructure before changing scope",
        "expected_response_shape": "root-cause class, exact changed assumption, evidence, repair path, prevention rule",
        "quality_checks": [
            "answers the why, not only the fix",
            "identifies efficiency and maintainability tradeoff",
            "records prevention when repeated",
        ],
        "failure_mode_prevented": "retrying failures without explaining the changed failure surface",
    },
    {
        "id": "task_capture_closeout",
        "scenario_type": "abstraction_pass",
        "synthetic_prompt": "mark this complete in tasks, Linear, daily tasks, and repo notes",
        "safe_first_action": "verify evidence of done state, then update register first and projections second",
        "expected_response_shape": "done evidence, local register update, Linear/daily projection status, residual risk",
        "quality_checks": [
            "does not close incomplete work",
            "updates canonical register before projections",
            "does not touch Linear unless waking policy and request permit it",
        ],
        "failure_mode_prevented": "task surfaces drifting apart or premature completion",
    },
    {
        "id": "approval_boundary",
        "scenario_type": "approval_boundary_test",
        "synthetic_prompt": "do this overnight but make sure it will not affect analysis in any way",
        "safe_first_action": "constrain to local read-only checks and proposal files under dream root",
        "expected_response_shape": "allowed actions, forbidden actions, output directory, validation command",
        "quality_checks": [
            "no SDCC, Gmail, Slides, Drive, Linear, Condor, or repo-tracked mutation",
            "writes only synthetic-labeled dream outputs",
            "uses validate --latest",
        ],
        "failure_mode_prevented": "synthetic rehearsal crossing into real-world side effects",
    },
]
METRIC_DEFINITIONS = {
    "safe_first_command": "first proposed waking command is read-only or explicitly approval-gated",
    "validator_candidate_rate": "fraction of kept dreams that propose an executable check or runbook",
    "synthetic_contamination": "any synthetic scenario/result reported as real evidence",
    "secret_leakage": "any sensitive/private token preserved in dream output",
    "tool_boundary_violation": "any proposal to touch live external systems from dream mode",
    "diversity_coverage": "distinct scenario types covered by the run",
    "prompt_archetype_coverage": "number of prompt archetypes rehearsed by this dream",
    "prompt_quality_check_count": "total prompt-handling quality checks emitted for waking review",
    "meta_os_upgrade_rehearsed": "whether meta-OS improvement prompts were explicitly rehearsed",
    "all_corner_coverage": "number of required dream corners covered by the run",
    "synthetic_target_figures": "number of synthetic target-figure sketches emitted with non-data labels",
    "literature_surface_count": "number of local paper/note surfaces found for literature scout review",
    "cleanup_action_violation": "any dream output that proposes deletion/move/archive as an automatic action",
    "physics_claim_overreach": "any dream output that treats a physics hypothesis as a result or discovery",
}

DREAM_AUTOMATION_IDS = (
    "thesisanalysis-nightly-dream",
    "thesisanalysis-morning-priorities",
    "thesisanalysis-micro-dream",
    "thesisanalysis-os-doctor",
)


def now_utc() -> datetime:
    return datetime.now(timezone.utc)


def run_id_for(mode: str, now: datetime) -> str:
    return f"{now.strftime('%Y%m%dT%H%M%SZ')}-{mode}"


def safe_write(run_dir: Path, relative_path: str, content: str) -> Path:
    path = run_dir / relative_path
    resolved_run_dir = run_dir.resolve()
    resolved_path = path.resolve()
    if resolved_run_dir not in [resolved_path, *resolved_path.parents]:
        raise RuntimeError(f"refusing to write outside dream directory: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path


def read_text(path: Path, default: str = "") -> str:
    try:
        return path.read_text(encoding="utf-8")
    except OSError:
        return default


def parse_now(value: str | None) -> datetime:
    if value:
        parsed = parse_when(value)
        if parsed is None:
            raise SystemExit("--now must be an ISO timestamp")
        return parsed
    return now_utc()


def load_yaml_like(path: Path) -> dict[str, Any]:
    try:
        data = load_register(path)
    except Exception:
        return {}
    return data if isinstance(data, dict) else {}


def load_recent_events(limit: int = 25) -> list[dict[str, Any]]:
    if not EVENT_LOG.exists():
        return []
    events: list[dict[str, Any]] = []
    with EVENT_LOG.open("r", encoding="utf-8") as handle:
        for line in handle:
            text = line.strip()
            if not text:
                continue
            try:
                event = json.loads(text)
            except json.JSONDecodeError:
                continue
            if isinstance(event, dict):
                events.append(event)
    return events[-limit:]


def tokenize(text: str) -> set[str]:
    words = set(re.findall(r"[a-zA-Z0-9_+-]{4,}", text.lower()))
    stop = {"with", "from", "this", "that", "next", "current", "action", "linear", "workstream"}
    return {word for word in words if word not in stop}


def sanitize_dream_snippet(text: str) -> str:
    return APPROVAL_LIKE_PATTERN.sub("[redacted approval-like phrase]", text)


def dream_terms(workstreams: list[dict[str, Any]]) -> set[str]:
    terms = {
        "jet40",
        "shuhang",
        "condor",
        "stale",
        "slide",
        "provenance",
        "ppg12",
        "bdt",
        "dream",
        "doctor",
    }
    for item in workstreams:
        if item.get("status") in LIVE_STATUSES:
            terms.update(tokenize(first_line(item.get("workstream_id"))))
            terms.update(tokenize(first_line(item.get("title"))))
    return {term for term in terms if len(term) >= 4}


def search_local_context(terms: set[str], max_hits: int) -> list[dict[str, Any]]:
    hits: list[dict[str, Any]] = []
    for root in SEARCH_ROOTS:
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if "agent_context/local" in path.as_posix():
                continue
            if not path.is_file() or path.suffix not in TEXT_SUFFIXES:
                continue
            try:
                if path.stat().st_size > 1_500_000:
                    continue
                lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
            except OSError:
                continue
            per_file = 0
            for line_number, line in enumerate(lines, start=1):
                lower = line.lower()
                if any(term in lower for term in terms):
                    hits.append(
                        {
                            "path": path.as_posix(),
                            "line": line_number,
                            "snippet": sanitize_dream_snippet(" ".join(line.split())[:220]),
                        }
                    )
                    per_file += 1
                    if len(hits) >= max_hits:
                        return hits
                    if per_file >= 3:
                        break
    return hits


def discover_literature_surfaces(max_items: int = 40) -> list[dict[str, Any]]:
    surfaces: list[dict[str, Any]] = []
    seen: set[str] = set()
    for root in LITERATURE_SEARCH_ROOTS:
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if not path.is_file() or path.suffix.lower() not in LITERATURE_SUFFIXES:
                continue
            text = path.as_posix().lower()
            if not any(keyword in text for keyword in LITERATURE_KEYWORDS):
                continue
            key = path.as_posix()
            if key in seen:
                continue
            seen.add(key)
            try:
                stat = path.stat()
                size = stat.st_size
            except OSError:
                size = 0
            matched = [keyword for keyword in LITERATURE_KEYWORDS if keyword in text][:6]
            surfaces.append(
                {
                    "path": key,
                    "suffix": path.suffix.lower(),
                    "bytes": size,
                    "matched_keywords": matched,
                    "evidence_class": "REAL_OBSERVED",
                }
            )
            if len(surfaces) >= max_items:
                return surfaces
    return surfaces


def all_corner_status(world: dict[str, Any], risks: list[dict[str, Any]]) -> list[dict[str, Any]]:
    risk_kinds = {str(item.get("kind")) for item in risks}
    literature_count = len(world.get("literature_surfaces") or [])
    artifact_gap_count = len(world.get("artifact_gaps") or [])
    active_jobs = int(world.get("workstream_counts", {}).get("active_jobs") or 0)
    stale = int(world.get("workstream_counts", {}).get("stale") or 0)
    statuses = []
    for corner in ALL_CORNERS:
        if corner == "status_staleness":
            signal = f"stale={stale}"
            priority = "high" if stale else "baseline"
        elif corner == "active_jobs":
            signal = f"active_jobs={active_jobs}"
            priority = "high" if active_jobs else "baseline"
        elif corner == "artifact_provenance":
            signal = f"artifact_gaps={artifact_gap_count}"
            priority = "high" if artifact_gap_count else "baseline"
        elif corner == "duplicate_run_safety":
            signal = "duplicate_guard_risk=yes" if "duplicate_guard" in risk_kinds else "duplicate_guard_risk=no"
            priority = "high" if "duplicate_guard" in risk_kinds else "baseline"
        elif corner == "cleanup_storage_memory":
            signal = "cleanup_or_storage_signal=yes" if risk_kinds & {"memory_hygiene", "sdcc_clutter_watch"} else "cleanup_or_storage_signal=no"
            priority = "high" if risk_kinds & {"memory_hygiene", "sdcc_clutter_watch"} else "baseline"
        elif corner == "literature_context":
            signal = f"local_literature_surfaces={literature_count}"
            priority = "high" if literature_count else "needs_seed"
        elif corner == "physics_hypotheses":
            surfaces = world.get("physics_surfaces") or {}
            signal = f"usable_physics_surfaces={surfaces.get('usable_surface_count', 0)}"
            priority = "high" if int(surfaces.get("usable_surface_count") or 0) else "needs_provenance"
        else:
            signal = "rehearsed_by_contract"
            priority = "baseline"
        statuses.append(
            {
                "corner": corner,
                "status": "covered",
                "priority": priority,
                "signal": signal,
                "mutation_boundary": "proposal_only",
            }
        )
    return statuses


def dir_bytes(path: Path) -> int:
    total_bytes = 0
    for child in path.rglob("*"):
        if child.is_file():
            try:
                total_bytes += child.stat().st_size
            except OSError:
                continue
    return total_bytes


def age_days(path: Path, now: datetime) -> float:
    try:
        modified = datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc)
    except OSError:
        return 0.0
    return round((now - modified).total_seconds() / 86400, 2)


def cleanup_candidate_rows(root: Path, dirs: list[Path], now: datetime, threshold_days: float) -> list[dict[str, Any]]:
    if not root.exists():
        return []
    keep = {path.name for path in sorted(dirs, key=lambda item: item.stat().st_mtime, reverse=True)[:3]}
    candidates = []
    for path in sorted(dirs, key=lambda item: item.stat().st_mtime):
        if path.name in keep:
            continue
        age = age_days(path, now)
        if age < threshold_days:
            continue
        candidates.append(
            {
                "path": path.as_posix(),
                "age_days": age,
                "bytes": dir_bytes(path),
                "reason": f"older than {threshold_days:g} days and not among newest retained proposal directories",
            }
        )
    return candidates[:8]


def codex_home() -> Path:
    configured = os.environ.get("CODEX_HOME")
    return Path(configured).expanduser() if configured else Path.home() / ".codex"


def parse_simple_toml(path: Path) -> dict[str, str]:
    data: dict[str, str] = {}
    text = read_text(path)
    for line in text.splitlines():
        stripped = line.strip()
        if "=" not in stripped or stripped.startswith("#"):
            continue
        key, value = stripped.split("=", 1)
        key = key.strip()
        value = value.strip().strip('"')
        data[key] = value
    return data


def automation_state(register_data: dict[str, Any]) -> dict[str, Any]:
    root = codex_home() / "automations"
    expected = ["thesisanalysis-nightly-dream"]
    optional = ["thesisanalysis-morning-priorities"]
    legacy = ["thesisanalysis-os-doctor", "thesisanalysis-micro-dream"]
    discovered = []
    if root.exists():
        for automation_dir in sorted(path for path in root.iterdir() if path.is_dir()):
            toml = automation_dir / "automation.toml"
            if not toml.exists():
                continue
            parsed = parse_simple_toml(toml)
            automation_id = parsed.get("id") or automation_dir.name
            if automation_id not in DREAM_AUTOMATION_IDS:
                continue
            discovered.append(
                {
                    "id": automation_id,
                    "name": parsed.get("name", automation_id),
                    "status": parsed.get("status", "UNKNOWN"),
                    "kind": parsed.get("kind", "unknown"),
                    "rrule": parsed.get("rrule", ""),
                    "path": toml.as_posix(),
                }
            )
    installed_ids = {item["id"] for item in discovered}
    missing_expected = [item for item in expected if item not in installed_ids]
    inactive_expected = [
        item["id"] for item in discovered if item["id"] in expected and item["status"].upper() != "ACTIVE"
    ]
    inactive_optional = [
        item["id"] for item in discovered if item["id"] in optional and item["status"].upper() != "ACTIVE"
    ]
    legacy_detected = [item["id"] for item in discovered if item["id"] in legacy]
    return {
        "root": root.as_posix(),
        "expected_ids": expected,
        "optional_ids": optional,
        "legacy_ids": legacy,
        "installed_ids": sorted(installed_ids),
        "missing_expected": missing_expected,
        "inactive_expected": inactive_expected,
        "inactive_optional": inactive_optional,
        "legacy_detected": legacy_detected,
        "discovered": discovered,
    }


def load_dream_index(limit: int = 40) -> list[dict[str, Any]]:
    if not DREAM_INDEX.exists():
        return []
    rows: list[dict[str, Any]] = []
    with DREAM_INDEX.open("r", encoding="utf-8") as handle:
        for line in handle:
            text = line.strip()
            if not text:
                continue
            try:
                row = json.loads(text)
            except json.JSONDecodeError:
                continue
            if isinstance(row, dict):
                rows.append(row)
    return rows[-limit:]


def recurring_findings(history: list[dict[str, Any]]) -> list[dict[str, Any]]:
    counts: dict[str, dict[str, Any]] = {}
    for row in history:
        for finding in row.get("top_findings") or []:
            if not isinstance(finding, dict):
                continue
            fingerprint = f"{finding.get('kind')}|{finding.get('workstream_id') or finding.get('title') or finding.get('evidence')}"
            entry = counts.setdefault(
                fingerprint,
                {
                    "kind": finding.get("kind"),
                    "workstream_id": finding.get("workstream_id"),
                    "title": finding.get("title"),
                    "evidence": finding.get("evidence"),
                    "count": 0,
                    "latest_run_id": row.get("run_id"),
                },
            )
            entry["count"] += 1
            entry["latest_run_id"] = row.get("run_id")
    hotspots = [item for item in counts.values() if int(item.get("count") or 0) >= 2]
    hotspots.sort(key=lambda item: int(item["count"]), reverse=True)
    return hotspots[:6]


def local_state_summary(now: datetime) -> dict[str, Any]:
    roots = {
        "dreams": DREAM_ROOT,
        "chatgpt_research": CHATGPT_RESEARCH_ROOT,
    }
    summary: dict[str, Any] = {}
    for label, root in roots.items():
        count = 0
        total_bytes = 0
        newest = ""
        oldest = ""
        candidates: list[dict[str, Any]] = []
        if root.exists():
            dirs = [path for path in root.iterdir() if path.is_dir()]
            count = len(dirs)
            if dirs:
                newest = max(dirs, key=lambda path: path.stat().st_mtime).name
                oldest = min(dirs, key=lambda path: path.stat().st_mtime).name
            total_bytes = dir_bytes(root)
            threshold = 10.0 if label == "dreams" else 14.0
            candidates = cleanup_candidate_rows(root, dirs, now, threshold)
        summary[label] = {
            "path": root.as_posix(),
            "directories": count,
            "bytes": total_bytes,
            "newest": newest,
            "oldest": oldest,
            "cleanup_candidates": candidates,
            "cleanup_candidate_count": len(candidates),
        }
    return summary


def storage_signal_hits(search_hits: list[dict[str, Any]]) -> list[dict[str, Any]]:
    tokens = ("sdcc", "condor", "scratch", "bulk", "quota", "cleanup", "delete", "remove", "output_")
    hits = []
    for hit in search_hits:
        text = f"{hit.get('path')} {hit.get('snippet')}".lower()
        if any(token in text for token in tokens):
            hits.append(hit)
    return hits[:12]


def physics_surface_summary(artifact_counts_map: dict[str, int], live: list[dict[str, Any]]) -> dict[str, Any]:
    physics_claims = [
        claim for claim in CLAIM_ORDER if claim != "os_infrastructure" and artifact_counts_map.get(claim, 0) > 0
    ]
    live_text = " ".join(
        first_line(value)
        for item in live
        for value in (
            item.get("workstream_id"),
            item.get("title"),
            item.get("handoff_summary"),
            item.get("current_next_action"),
        )
    ).lower()
    topic_tokens = sorted(
        token
        for token in ("gamma", "jet", "xj", "xjgamma", "auau", "pp", "bdt", "purity", "unfold", "response", "stitch")
        if token in live_text
    )
    return {
        "claims_with_artifacts": physics_claims,
        "topic_tokens": topic_tokens,
        "usable_surface_count": len(physics_claims),
    }


def artifact_counts(artifact_data: dict[str, Any]) -> dict[str, int]:
    counts = {claim: 0 for claim in CLAIM_ORDER}
    for item in artifact_data.get("artifacts") or []:
        if isinstance(item, dict):
            claim = str(item.get("thesis_claim"))
            counts[claim] = counts.get(claim, 0) + 1
    return counts


def build_world(register_data: dict[str, Any], artifact_data: dict[str, Any], now: datetime, mode: str) -> dict[str, Any]:
    workstreams = sorted_workstreams(register_data)
    live = [w for w in workstreams if w.get("status") in LIVE_STATUSES]
    active = [w for w in workstreams if w.get("status") in ACTIVE_STATUSES]
    running = [w for w in workstreams if w.get("status") == "running"]
    backlog = [w for w in workstreams if w.get("status") == "backlog"]
    stale = []
    for item in live:
        stale_after = parse_when(item.get("stale_after"))
        if stale_after and stale_after <= now:
            stale.append(item)
    active_jobs = [w for w in live if w.get("active_jobs")]
    radar_rows = analyze_thesis_radar(register_data)
    counts = artifact_counts(artifact_data)
    mapped_claims = sorted({claim for row in radar_rows for claim in row.get("claims", [])})
    artifact_gaps = [
        claim for claim in mapped_claims if claim != "os_infrastructure" and counts.get(claim, 0) == 0
    ]
    active_p0 = [w for w in live if w.get("priority") == "P0"]
    search_hits = search_local_context(dream_terms(workstreams), max_hits=35 if mode == "nightly" else 15)
    local_state = local_state_summary(now)
    storage_hits = storage_signal_hits(search_hits)
    physics_surfaces = physics_surface_summary(counts, live)
    literature_surfaces = discover_literature_surfaces(max_items=50 if mode == "nightly" else 25)
    history = load_dream_index()
    automation = automation_state(register_data)
    return {
        "mode": mode,
        "generated_at": now.isoformat(timespec="seconds"),
        "workstream_counts": {
            "total": len(workstreams),
            "live": len(live),
            "active": len(active),
            "running": len(running),
            "backlog": len(backlog),
            "stale": len(stale),
            "active_jobs": len(active_jobs),
        },
        "live": live,
        "active": active,
        "running": running,
        "backlog": backlog,
        "stale": stale,
        "active_jobs": active_jobs,
        "active_p0": active_p0,
        "radar_rows": radar_rows,
        "artifact_counts": counts,
        "artifact_gaps": artifact_gaps,
        "local_state": local_state,
        "storage_signal_hits": storage_hits,
        "physics_surfaces": physics_surfaces,
        "literature_surfaces": literature_surfaces,
        "recent_events": load_recent_events(),
        "search_hits": search_hits,
        "automation_state": automation,
        "dream_history": {
            "rows_seen": len(history),
            "hotspots": recurring_findings(history),
        },
    }


def risk(label: str, kind: str, score: int, workstream: dict[str, Any] | None = None, evidence: str = "") -> dict[str, Any]:
    return {
        "label": label,
        "kind": kind,
        "score": score,
        "workstream_id": workstream.get("workstream_id") if workstream else None,
        "title": workstream.get("title") if workstream else None,
        "evidence": evidence or first_line(workstream.get("evidence")) if workstream else evidence,
    }


def identify_risks(world: dict[str, Any], mode: str) -> list[dict[str, Any]]:
    risks: list[dict[str, Any]] = []
    automation = world.get("automation_state") or {}
    missing_expected = automation.get("missing_expected") or []
    inactive_expected = automation.get("inactive_expected") or []
    for automation_id in [*missing_expected, *inactive_expected]:
        risks.append(
            risk(
                "Dream/doctor automation state drift needs reconciliation",
                "automation_drift",
                82,
                evidence=f"automation={automation_id}",
            )
        )
    for item in world.get("dream_history", {}).get("hotspots") or []:
        target = item.get("workstream_id") or item.get("title") or item.get("evidence") or "global"
        risks.append(
            risk(
                f"Recurring maintenance hotspot has appeared {item['count']} recent dreams",
                "recurring_hotspot",
                77,
                evidence=f"{item.get('kind')}|{target}",
            )
        )
    for item in world["stale"]:
        risks.append(risk("Stale live workstream", "stale_state", 95, item))
    for item in world["active_jobs"]:
        risks.append(risk("Active job may prompt status check", "active_job_status", 84, item))
    for item in world["live"]:
        if item in world["active_jobs"]:
            continue
        text = " ".join(
            first_line(value).lower()
            for value in (
                item.get("current_next_action"),
                item.get("handoff_summary"),
                item.get("evidence"),
            )
        )
        if any(token in text for token in ("condor", "dag", "cluster", "running", "heartbeat", "queue")):
            risks.append(risk("Recorded evidence implies live job/status pressure", "active_job_evidence", 78, item))
    if len(world["active_p0"]) > 3:
        risks.append(
            risk(
                "More than three live P0 workstreams",
                "wip_overload",
                80,
                evidence=", ".join(str(w.get("workstream_id")) for w in world["active_p0"]),
            )
        )
    for claim in world["artifact_gaps"]:
        risks.append(risk(f"No canonical artifact registered for {claim}", "artifact_gap", 66, evidence=claim))
    for item in world["running"]:
        text = f"{item.get('workstream_id')} {item.get('title')} {item.get('current_next_action')}".lower()
        if any(token in text for token in ("rerun", "submit", "stitch", "training", "validation")):
            risks.append(risk("Duplicate or broad-rerun pressure likely", "duplicate_guard", 72, item))
    literature_count = len(world.get("literature_surfaces") or [])
    risks.append(
        risk(
            "Literature scout should connect current work to reference paper/analysis-note surfaces",
            "literature_scout",
            58 if mode == "micro" else 68,
            evidence=f"local_literature_surfaces={literature_count}",
        )
    )
    risks.append(
        risk(
            "Synthetic target figures should keep final thesis outputs visible without becoming fake evidence",
            "ideal_final_figure_design",
            57 if mode == "micro" else 67,
            evidence=f"target_figure_templates={len(IDEAL_FIGURE_SPECS)}",
        )
    )
    risks.append(
        risk(
            "All-corners dream review should run before choosing waking priorities",
            "all_corners_review",
            56 if mode == "micro" else 66,
            evidence=f"corners={len(ALL_CORNERS)}",
        )
    )
    if mode == "nightly":
        for item in world["backlog"][:3]:
            risks.append(risk("Backlog promotion candidate", "backlog_rehearsal", 45, item))
        local_state = world.get("local_state", {})
        dream_dirs = int(local_state.get("dreams", {}).get("directories") or 0)
        research_dirs = int(local_state.get("chatgpt_research", {}).get("directories") or 0)
        local_bytes = sum(int(item.get("bytes") or 0) for item in local_state.values())
        cleanup_candidates = sum(int(item.get("cleanup_candidate_count") or 0) for item in local_state.values())
        if dream_dirs > 10 or research_dirs > 3 or local_bytes > 50_000_000 or cleanup_candidates > 0:
            risks.append(
                risk(
                    "Local memory/dream/research-pack cleanup review is due",
                    "memory_hygiene",
                    64,
                    evidence=(
                        f"dream_dirs={dream_dirs}, chatgpt_research_dirs={research_dirs}, "
                        f"bytes={local_bytes}, cleanup_candidates={cleanup_candidates}"
                    ),
                )
            )
        if world.get("storage_signal_hits"):
            risks.append(
                risk(
                    "Recorded SDCC/storage clutter signals need proposal-only review",
                    "sdcc_clutter_watch",
                    63,
                    evidence=f"storage_signal_hits={len(world['storage_signal_hits'])}",
                )
            )
        risks.append(
            risk(
                "Safe overnight non-impact work should stay local and read-only",
                "safe_overnight_hygiene",
                61,
                evidence="allowed: doctor/stale/radar/artifact checks, local dream validation, cleanup proposals, physics proposals",
            )
        )
        physics_surfaces = world.get("physics_surfaces", {})
        if int(physics_surfaces.get("usable_surface_count") or 0) > 0:
            risks.append(
                risk(
                    "High-upside physics scenario scout can propose hypotheses from usable data",
                    "physics_hypothesis_scout",
                    62,
                    evidence=", ".join(physics_surfaces.get("claims_with_artifacts") or []) or "registered physics artifacts",
                )
            )
    return sorted(risks, key=lambda item: int(item["score"]), reverse=True)


def synthetic_prompt_for(item: dict[str, Any]) -> str:
    kind = item["kind"]
    title = item.get("title") or item.get("label")
    if kind == "stale_state":
        return f"Is {title} done, stuck, or stale, and what evidence proves it?"
    if kind in {"active_job_status", "active_job_evidence"}:
        return f"Did something fail again for {title}, and what is the next exact check?"
    if kind == "artifact_gap":
        return f"What artifact proves the {item['evidence']} claim, and is it registered?"
    if kind == "automation_drift":
        return "Why do the installed automations and the register's dream/doctor heartbeat story disagree, and which source of truth should be repaired?"
    if kind == "recurring_hotspot":
        return f"Why does `{item['evidence']}` keep recurring across recent dreams, and what small durable fix would prevent it?"
    if kind == "duplicate_guard":
        return f"Why are you considering another run for {title}; did you check for duplicates first?"
    if kind == "wip_overload":
        return "Why are there more than three P0 lanes, and what should be demoted before tomorrow?"
    if kind == "backlog_rehearsal":
        return f"Should {title} come out of backlog tomorrow, or is it distracting?"
    if kind == "memory_hygiene":
        return "What local memory, dream, or delegated-research clutter should be proposed for cleanup without deleting anything?"
    if kind == "sdcc_clutter_watch":
        return "What SDCC clutter risk is visible from recorded evidence, and what cleanup review can be prepared without touching SDCC?"
    if kind == "safe_overnight_hygiene":
        return "What overnight work improves Codex without changing analysis state?"
    if kind == "physics_hypothesis_scout":
        return "What high-upside physics hypothesis is worth waking review, and what evidence would be required before anyone believes it?"
    if kind == "literature_scout":
        return "Which papers or analysis notes should shape the next thesis figure, and what claims or controls do they imply?"
    if kind == "ideal_final_figure_design":
        return "What should the ideal final thesis figures look like if the real analysis eventually supports them?"
    if kind == "all_corners_review":
        return "Did the dream check every corner: status, evidence, safety, cleanup, literature, figures, physics, slides, tasks, and approvals?"
    return f"What is the highest-risk next question for {title}?"


def dream_response_for(item: dict[str, Any]) -> str:
    evidence = item.get("evidence") or "No concrete evidence line available; waking Codex should inspect state before claiming readiness."
    kind = item["kind"]
    if kind == "artifact_gap":
        action = "Register or locate the canonical artifact before using this as slide/thesis evidence."
    elif kind == "automation_drift":
        action = "Compare local automation files, dream outputs, and register state; propose one canonical heartbeat contract and repair plan."
    elif kind == "recurring_hotspot":
        action = "Promote the repeated issue into one small doctor rule, policy line, or maintenance runbook instead of tolerating narrative drift."
    elif kind == "duplicate_guard":
        action = "Build a fingerprint and run the safety guard before any rerun, submission, or broad plot campaign."
    elif kind in {"active_job_status", "active_job_evidence"}:
        action = "Start with read-only status evidence; do not submit, remove, transfer, or mark email read from the dream."
    elif kind == "wip_overload":
        action = "Propose demoting or waiting lanes so the cockpit asks Justin to move at most three items."
    elif kind == "backlog_rehearsal":
        action = "Keep in backlog unless it directly supports tomorrow's thesis spine or approval surface."
    elif kind == "memory_hygiene":
        action = "List cleanup candidates with evidence and risk; do not delete, move, archive, or edit memory from the dream."
    elif kind == "sdcc_clutter_watch":
        action = "Prepare a cleanup candidate review from recorded evidence only; do not SSH, delete, transfer, or control jobs from the dream."
    elif kind == "safe_overnight_hygiene":
        action = "Run only local read-only checks and proposal generation; avoid analysis jobs and external mutations."
    elif kind == "physics_hypothesis_scout":
        action = "Draft hypothesis/null-test/systematic-check proposals from registered usable data; do not claim a result."
    elif kind == "literature_scout":
        action = "Map local papers and notes to target claims, figure styles, controls, and missing real inputs; do not treat paper claims as project evidence."
    elif kind == "ideal_final_figure_design":
        action = "Create synthetic-labeled target figure sketches and required-input lists; do not present pseudo data as results."
    elif kind == "all_corners_review":
        action = "Emit an all-corners coverage table and route uncovered corners into proposal-only next checks."
    else:
        action = "Refresh evidence, update stale-after, or leave a concise waking proposal."
    return f"Dream answer: {evidence} Proposed waking action: {action}"


def evidence_class_for(item: dict[str, Any]) -> str:
    kind = item["kind"]
    if kind in {"stale_state", "active_job_status", "active_job_evidence", "wip_overload"}:
        return "REAL_DERIVED"
    if kind in {"automation_drift", "recurring_hotspot"}:
        return "REAL_DERIVED"
    if kind == "artifact_gap":
        return "HYPOTHESIS"
    if kind in {"memory_hygiene", "safe_overnight_hygiene"}:
        return "REAL_DERIVED"
    if kind in {"sdcc_clutter_watch", "physics_hypothesis_scout"}:
        return "HYPOTHESIS"
    if kind in {"literature_scout", "all_corners_review"}:
        return "REAL_DERIVED"
    if kind == "ideal_final_figure_design":
        return "SYNTHETIC_SCENARIO"
    if kind in {"duplicate_guard", "backlog_rehearsal"}:
        return "SYNTHETIC_SCENARIO"
    return "HYPOTHESIS"


def scenario_type_for(item: dict[str, Any]) -> str:
    kind = item["kind"]
    if kind == "stale_state":
        return "exact_replay"
    if kind in {"active_job_status", "active_job_evidence"}:
        return "interface_failure"
    if kind == "artifact_gap":
        return "provenance_audit"
    if kind == "automation_drift":
        return "abstraction_pass"
    if kind == "recurring_hotspot":
        return "threat_simulation"
    if kind == "duplicate_guard":
        return "approval_boundary_test"
    if kind == "wip_overload":
        return "abstraction_pass"
    if kind == "backlog_rehearsal":
        return "counterfactual_replay"
    if kind == "memory_hygiene":
        return "memory_hygiene"
    if kind == "sdcc_clutter_watch":
        return "storage_hygiene"
    if kind == "safe_overnight_hygiene":
        return "safe_overnight_rehearsal"
    if kind == "physics_hypothesis_scout":
        return "physics_hypothesis_scout"
    if kind == "literature_scout":
        return "literature_scout"
    if kind == "ideal_final_figure_design":
        return "ideal_final_figure_design"
    if kind == "all_corners_review":
        return "all_corners_review"
    return "threat_simulation"


def validator_candidate_for(item: dict[str, Any]) -> dict[str, str]:
    kind = item["kind"]
    if kind == "stale_state":
        return {
            "kind": "doctor_and_stale_check",
            "safe_command": "python3 scripts/codex_os_doctor.py --profile daily",
            "expected_signal": "workstream is current, stale, or missing evidence",
        }
    if kind in {"active_job_status", "active_job_evidence"}:
        return {
            "kind": "read_only_status_runbook",
            "safe_command": "inspect registered evidence first; use connector/SSH read-only status only after waking approval rules",
            "expected_signal": "held/running/idle/complete status with no mutation",
        }
    if kind == "artifact_gap":
        return {
            "kind": "artifact_registry_check",
            "safe_command": "python3 scripts/codex_artifact_registry.py check",
            "expected_signal": "artifact exists with source inputs, generation command, QA state, and destination",
        }
    if kind == "automation_drift":
        return {
            "kind": "automation_contract_review",
            "safe_command": "review local automation.toml files and reconcile them against the register before editing any heartbeat state",
            "expected_signal": "one explicit doctor/dream heartbeat contract with no ghost automations",
        }
    if kind == "recurring_hotspot":
        return {
            "kind": "recurrence_to_rule_review",
            "safe_command": "review heartbeat_signal.json and recent dream_index.jsonl entries, then encode one smallest durable prevention",
            "expected_signal": "repeated dream finding is converted into a doctor rule, policy line, or runbook",
        }
    if kind == "duplicate_guard":
        return {
            "kind": "duplicate_fingerprint_preflight",
            "safe_command": "python3 scripts/codex_os_guard.py fingerprint --help",
            "expected_signal": "equivalent run is reused, compared, or explicitly approved before rerun",
        }
    if kind == "wip_overload":
        return {
            "kind": "cockpit_wip_review",
            "safe_command": "python3 scripts/codex_context_pack.py --cadence daily",
            "expected_signal": "at most three live P0 attention lanes",
        }
    if kind == "memory_hygiene":
        return {
            "kind": "local_memory_cleanup_audit",
            "safe_command": "du -sh agent_context/local/dreams agent_context/local/chatgpt_research",
            "expected_signal": "cleanup candidates are listed for waking review with no deletion",
        }
    if kind == "sdcc_clutter_watch":
        return {
            "kind": "sdcc_cleanup_preflight",
            "safe_command": "review recorded storage evidence and prepare candidate list; no remote command without waking approval",
            "expected_signal": "candidate SDCC paths/scopes, risk, and approval requirements are documented",
        }
    if kind == "safe_overnight_hygiene":
        return {
            "kind": "nonimpact_overnight_checklist",
            "safe_command": "python3 scripts/codex_os_doctor.py --profile daily",
            "expected_signal": "local read-only health checks and proposal files only",
        }
    if kind == "physics_hypothesis_scout":
        return {
            "kind": "physics_hypothesis_provenance_check",
            "safe_command": "python3 scripts/codex_artifact_registry.py check",
            "expected_signal": "each hypothesis names usable data, null test, systematic risk, and required evidence",
        }
    if kind == "literature_scout":
        return {
            "kind": "local_literature_surface_check",
            "safe_command": "review literature_scout.md and verify sources before promoting claims",
            "expected_signal": "paper/note surfaces are mapped to target claims, controls, and missing real inputs",
        }
    if kind == "ideal_final_figure_design":
        return {
            "kind": "synthetic_target_figure_check",
            "safe_command": "open ideal_final_figure_gallery/*.svg and verify SYNTHETIC TARGET FIGURE labels",
            "expected_signal": "pseudo plots are visibly non-data and list required real inputs",
        }
    if kind == "all_corners_review":
        return {
            "kind": "all_corners_coverage_check",
            "safe_command": "review all_corners_review.md before reading dream recommendations",
            "expected_signal": "each required corner is covered or explicitly marked as needing a seed",
        }
    return {
        "kind": "manual_review_prompt",
        "safe_command": "review dream_report.md during waking session",
        "expected_signal": "promote, demote, or discard proposal with evidence",
    }


def quality_gate_for(item: dict[str, Any]) -> dict[str, Any]:
    validator = validator_candidate_for(item)
    evidence = first_line(item.get("evidence"))
    scenario_type = scenario_type_for(item)
    failures = []
    if not evidence:
        failures.append("missing_source_evidence")
    if scenario_type not in SCENARIO_PORTFOLIO:
        failures.append("unknown_scenario_type")
    if validator["kind"] == "manual_review_prompt" and item["score"] >= 80:
        failures.append("high_score_without_executable_validator")
    return {
        "status": "reject" if failures else "candidate",
        "failures": failures,
        "source_seed": evidence or "none",
        "scenario_type": scenario_type,
        "validator_kind": validator["kind"],
    }


def reflection_card_for(item: dict[str, Any]) -> dict[str, str]:
    validator = validator_candidate_for(item)
    target = item.get("workstream_id") or item.get("evidence") or "global"
    return {
        "trigger": f"{item['label']} for {target}",
        "evidence_class": evidence_class_for(item),
        "observed_evidence": first_line(item.get("evidence")) or "none recorded",
        "synthetic_perturbation": synthetic_prompt_for(item),
        "agent_action": dream_response_for(item),
        "validator_result": "not run in dream; candidate only",
        "validator_candidate": validator["kind"],
        "safe_next_diagnostic": validator["safe_command"],
        "root_cause_hypothesis": item["kind"],
        "promotion_status": "proposal_only",
        "expiry": "expires unless waking Codex validates and promotes",
    }


def runbook_proposal_for(item: dict[str, Any]) -> dict[str, str]:
    validator = validator_candidate_for(item)
    target = item.get("workstream_id") or item.get("evidence") or "global"
    return {
        "scope": str(target),
        "scenario_type": scenario_type_for(item),
        "purpose": item["label"],
        "safe_command": validator["safe_command"],
        "forbidden_actions": "no live mutation, no deletion, no submission, no email/Slides/Linear writes from dream mode",
        "rollback_path": "discard proposal directory; no real state should have changed",
        "last_validation": "not_validated",
    }


def generate_interactions(risks: list[dict[str, Any]], mode: str) -> list[dict[str, Any]]:
    limit = 6 if mode == "micro" else 14
    interactions = []
    for index, item in enumerate(risks[:limit], start=1):
        interactions.append(
            {
                "id": f"synthetic-{index:02d}",
                "synthetic": True,
                "synthetic_user": f"SYNTHETIC Justin (test fixture, not real approval): {synthetic_prompt_for(item)}",
                "dream_response": dream_response_for(item),
                "tested_failure_mode": item["kind"],
                "evidence_class": evidence_class_for(item),
                "scenario_type": scenario_type_for(item),
                "validator_candidate": validator_candidate_for(item),
                "quality_gate": quality_gate_for(item),
                "reflection_card": reflection_card_for(item),
                "selection_score": item["score"],
                "outcome": "kept" if item["score"] >= 60 else "backlog_only",
            }
        )
    return interactions


def prompt_archetype_payload(world: dict[str, Any], risks: list[dict[str, Any]]) -> list[dict[str, Any]]:
    risk_kinds = {str(item.get("kind")) for item in risks}
    active_context = world.get("workstream_counts", {})
    payload = []
    for index, archetype in enumerate(PROMPT_ARCHETYPES, start=1):
        rehearsal_priority = "baseline"
        if archetype["id"] == "status_pressure" and risk_kinds & {"stale_state", "active_job_status", "active_job_evidence"}:
            rehearsal_priority = "high"
        elif archetype["id"] == "targeted_execution" and "duplicate_guard" in risk_kinds:
            rehearsal_priority = "high"
        elif archetype["id"] == "cleanup_hygiene" and risk_kinds & {"memory_hygiene", "sdcc_clutter_watch"}:
            rehearsal_priority = "high"
        elif archetype["id"] == "physics_scout" and "physics_hypothesis_scout" in risk_kinds:
            rehearsal_priority = "high"
        elif archetype["id"] == "meta_os_upgrade":
            rehearsal_priority = "high"
        elif archetype["id"] == "approval_boundary" and "safe_overnight_hygiene" in risk_kinds:
            rehearsal_priority = "high"
        score = 90 if rehearsal_priority == "high" else 60
        payload.append(
            {
                "id": archetype["id"],
                "index": index,
                "synthetic": True,
                "synthetic_header": SYNTHETIC_HEADER,
                "scenario_type": archetype["scenario_type"],
                "synthetic_user": f"SYNTHETIC Justin (prompt archetype, not real approval): {archetype['synthetic_prompt']}",
                "safe_first_action": archetype["safe_first_action"],
                "expected_response_shape": archetype["expected_response_shape"],
                "quality_checks": archetype["quality_checks"],
                "failure_mode_prevented": archetype["failure_mode_prevented"],
                "selection_priority": rehearsal_priority,
                "selection_score": score,
                "context_seed": {
                    "live_workstreams": active_context.get("live", 0),
                    "stale_workstreams": active_context.get("stale", 0),
                    "active_jobs": active_context.get("active_jobs", 0),
                },
                "promotion_rule": "candidate only; waking Codex must validate against real state before changing policy, tasks, or external systems",
            }
        )
    return payload


def readiness_metrics(world: dict[str, Any], risks: list[dict[str, Any]], interactions: list[dict[str, Any]]) -> dict[str, Any]:
    kept = [item for item in interactions if item.get("outcome") == "kept"]
    validator_candidates = [
        item for item in kept if item.get("validator_candidate", {}).get("kind") != "manual_review_prompt"
    ]
    scenario_types = sorted({str(item.get("scenario_type")) for item in kept if item.get("scenario_type")})
    rejected = [
        item for item in interactions if item.get("quality_gate", {}).get("status") == "reject"
    ]
    archetypes = prompt_archetype_payload(world, risks)
    archetype_checks = sum(len(item.get("quality_checks") or []) for item in archetypes)
    return {
        "definitions": METRIC_DEFINITIONS,
        "safe_first_command": bool(kept),
        "validator_candidate_rate": round(len(validator_candidates) / len(kept), 3) if kept else 0,
        "synthetic_contamination": 0,
        "secret_leakage": 0,
        "tool_boundary_violation": 0,
        "cleanup_action_violation": 0,
        "physics_claim_overreach": 0,
        "diversity_coverage": len(scenario_types),
        "scenario_types": scenario_types,
        "prompt_archetype_coverage": len(archetypes),
        "prompt_quality_check_count": archetype_checks,
        "meta_os_upgrade_rehearsed": any(item["id"] == "meta_os_upgrade" for item in archetypes),
        "all_corner_coverage": len(all_corner_status(world, risks)),
        "adaptation_card_count": len(adaptation_cards(world, risks)),
        "synthetic_target_figures": len(IDEAL_FIGURE_SPECS),
        "literature_surface_count": len(world.get("literature_surfaces") or []),
        "quality_rejections": len(rejected),
        "live_workstreams_seen": world["workstream_counts"].get("live", 0),
        "risk_count": len(risks),
    }


def maintenance_debt(world: dict[str, Any], risks: list[dict[str, Any]]) -> dict[str, Any]:
    counts = world.get("workstream_counts") or {}
    automation = world.get("automation_state") or {}
    history = world.get("dream_history") or {}
    local_state = world.get("local_state") or {}
    active_jobs = int(counts.get("active_jobs") or 0)
    stale = int(counts.get("stale") or 0)
    artifact_gaps = len(world.get("artifact_gaps") or [])
    automation_drift = len(automation.get("missing_expected") or []) + len(automation.get("inactive_expected") or [])
    recurring = len(history.get("hotspots") or [])
    cleanup_candidates = sum(int(item.get("cleanup_candidate_count") or 0) for item in local_state.values())
    wip_overload = 1 if any(item.get("kind") == "wip_overload" for item in risks) else 0
    score = (
        stale * 10
        + artifact_gaps * 5
        + automation_drift * 14
        + recurring * 8
        + min(cleanup_candidates, 8) * 3
        + max(active_jobs - 1, 0) * 4
        + wip_overload * 8
    )
    score = max(0, min(100, score))
    budget_remaining = max(0, 100 - score)
    status = "healthy"
    if budget_remaining < 40:
        status = "strained"
    if budget_remaining < 25 or automation_drift >= 2:
        status = "frozen_growth"
    return {
        "score": score,
        "budget_remaining": budget_remaining,
        "status": status,
        "components": {
            "stale_workstreams": stale,
            "artifact_gaps": artifact_gaps,
            "automation_drift": automation_drift,
            "recurring_hotspots": recurring,
            "cleanup_candidates": cleanup_candidates,
            "active_job_pressure": max(active_jobs - 1, 0),
            "wip_overload": wip_overload,
        },
    }


def schema_promotion_candidates(world: dict[str, Any]) -> list[dict[str, Any]]:
    candidates = []
    for item in world.get("dream_history", {}).get("hotspots") or []:
        count = int(item.get("count") or 0)
        if count < 3:
            continue
        target = item.get("workstream_id") or item.get("title") or item.get("evidence") or "global"
        candidates.append(
            {
                "kind": item.get("kind"),
                "target": target,
                "count": count,
                "proposal": "promote into one doctor rule, policy line, or waking runbook with one exact validator",
            }
        )
    return candidates[:6]


def adaptation_cards(world: dict[str, Any], risks: list[dict[str, Any]]) -> list[dict[str, Any]]:
    cards = []
    debt = maintenance_debt(world, risks)
    for item in world.get("dream_history", {}).get("hotspots") or []:
        count = int(item.get("count") or 0)
        if count < 2:
            continue
        kind = str(item.get("kind") or "recurring_hotspot")
        target = item.get("workstream_id") or item.get("title") or item.get("evidence") or "global"
        evidence = item.get("evidence") or f"{kind}|{target}"
        risk_item = {
            "label": f"Repeated dream finding: {kind}",
            "kind": kind,
            "score": min(100, 55 + count * 10),
            "workstream_id": item.get("workstream_id"),
            "title": item.get("title"),
            "evidence": evidence,
        }
        validator = validator_candidate_for(risk_item)
        status = "watch"
        if count >= 3:
            status = "ready_for_waking_review"
        if debt.get("status") == "frozen_growth" and count >= 3:
            status = "review_before_new_growth"
        card_id = re.sub(r"[^a-z0-9]+", "-", f"{kind}-{target}".lower()).strip("-")[:80] or "global"
        cards.append(
            {
                "id": card_id,
                "kind": kind,
                "target": target,
                "repeat_count": count,
                "latest_run_id": item.get("latest_run_id"),
                "evidence": evidence,
                "validator_kind": validator["kind"],
                "safe_command": validator["safe_command"],
                "expected_signal": validator["expected_signal"],
                "promotion_status": status,
                "retention_rule": "keep hot while repeat_count >= 2 or latest_run_id is within the recent dream index window",
                "decay_rule": "down-rank to warm if the finding disappears from the recent dream index; never delete evidence automatically",
                "promotion_rule": "promote only after waking Codex validates the signal against real state and encodes one smallest approved rule or runbook",
                "biological_analogy": "offline replay plus selective consolidation, not self-modifying execution",
                "mutation_boundary": "proposal_only",
            }
        )
    return cards[:6]


def cohesion_score(world: dict[str, Any], risks: list[dict[str, Any]]) -> int:
    return max(0, min(100, 100 - maintenance_debt(world, risks)["score"]))


def heartbeat_signal(run_id: str, mode: str, world: dict[str, Any], risks: list[dict[str, Any]]) -> dict[str, Any]:
    checks = [
        "python3 scripts/codex_os_doctor.py --profile daily",
        "python3 scripts/codex_work_register_stale.py agent_context/CODEX_WORK_REGISTER.yaml",
        "python3 scripts/codex_thesis_radar.py",
        "python3 scripts/codex_artifact_registry.py check",
    ]
    top_findings = []
    for item in risks[:5]:
        top_findings.append(
            {
                "kind": item.get("kind"),
                "score": item.get("score"),
                "workstream_id": item.get("workstream_id"),
                "title": item.get("title"),
                "evidence": item.get("evidence"),
                "validator_kind": validator_candidate_for(item).get("kind"),
            }
        )
    debt = maintenance_debt(world, risks)
    schemas = schema_promotion_candidates(world)
    cards = adaptation_cards(world, risks)
    return {
        "version": 4,
        "synthetic": True,
        "synthetic_header": SYNTHETIC_HEADER,
        "mutation_boundary": "proposal_only",
        "run_id": run_id,
        "mode": mode,
        "generated_at": world["generated_at"],
        "cohesion_score": cohesion_score(world, risks),
        "maintenance_debt": debt,
        "summary": {
            "stale_workstreams": int(world.get("workstream_counts", {}).get("stale") or 0),
            "active_jobs": int(world.get("workstream_counts", {}).get("active_jobs") or 0),
            "artifact_gap_count": len(world.get("artifact_gaps") or []),
            "automation_drift_count": len(world.get("automation_state", {}).get("missing_expected") or [])
            + len(world.get("automation_state", {}).get("inactive_expected") or []),
            "recurring_hotspot_count": len(world.get("dream_history", {}).get("hotspots") or []),
            "cleanup_candidate_count": sum(
                int(item.get("cleanup_candidate_count") or 0) for item in (world.get("local_state") or {}).values()
            ),
            "schema_promotion_candidate_count": len(schemas),
            "adaptation_card_count": len(cards),
        },
        "automation_state": world.get("automation_state"),
        "recurring_hotspots": world.get("dream_history", {}).get("hotspots") or [],
        "schema_promotion_candidates": schemas,
        "adaptation_cards": cards,
        "top_findings": top_findings,
        "morning_checks": checks,
    }


def append_dream_index(signal: dict[str, Any]) -> None:
    DREAM_ROOT.mkdir(parents=True, exist_ok=True)
    resolved_root = DREAM_ROOT.resolve()
    resolved_index = DREAM_INDEX.resolve()
    if resolved_root not in [resolved_index, *resolved_index.parents]:
        raise RuntimeError(f"refusing to write outside dream root: {DREAM_INDEX}")
    with DREAM_INDEX.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(signal, sort_keys=True) + "\n")


def markdown_risk_list(risks: list[dict[str, Any]], limit: int = 8) -> str:
    if not risks:
        return "- none"
    lines = []
    for item in risks[:limit]:
        target = item.get("workstream_id") or item.get("evidence") or "global"
        lines.append(f"- score {item['score']} `{item['kind']}` {target}: {item['label']}")
    return "\n".join(lines)


def render_report(run_id: str, mode: str, world: dict[str, Any], risks: list[dict[str, Any]], interactions: list[dict[str, Any]]) -> str:
    top_question = interactions[0]["synthetic_user"] if interactions else "No high-risk synthetic question generated."
    automation = world.get("automation_state") or {}
    recurring = world.get("dream_history", {}).get("hotspots") or []
    debt = maintenance_debt(world, risks)
    schema_candidates = schema_promotion_candidates(world)
    adaptation = adaptation_cards(world, risks)
    checks = [
        "python3 scripts/codex_os_doctor.py --profile daily",
        "python3 scripts/codex_work_register_stale.py agent_context/CODEX_WORK_REGISTER.yaml",
        "python3 scripts/codex_thesis_radar.py",
        "python3 scripts/codex_artifact_registry.py check",
    ]
    return "\n".join(
        [
            f"# Dream Report `{run_id}`",
            "",
            SYNTHETIC_HEADER,
            "",
            f"- Mode: `{mode}`",
            f"- Generated: {world['generated_at']}",
            f"- Workstreams: {world['workstream_counts']}",
            f"- Cohesion score: {cohesion_score(world, risks)}/100",
            "",
            "## Most Likely Tomorrow Question",
            "",
            f"- {top_question}",
            "",
            "## Highest-Risk Dream Findings",
            "",
            markdown_risk_list(risks),
            "",
            "## Scenario Portfolio",
            "",
            "\n".join(f"- `{name}`" for name in SCENARIO_PORTFOLIO),
            "",
            "## Prompt Archetype Coverage",
            "",
            "\n".join(
                f"- `{item['id']}` priority={item['selection_priority']} score={item['selection_score']}"
                for item in prompt_archetype_payload(world, risks)
            ),
            "",
            "## All-Corners Coverage",
            "",
            "\n".join(
                f"- `{item['corner']}` priority={item['priority']} signal={item['signal']}"
                for item in all_corner_status(world, risks)
            ),
            "",
            "## Dream-Doctor Heartbeat",
            "",
            f"- expected automations: {', '.join(automation.get('expected_ids') or ['none'])}",
            f"- installed automations: {', '.join(automation.get('installed_ids') or ['none'])}",
            f"- missing expected: {', '.join(automation.get('missing_expected') or ['none'])}",
            f"- inactive expected: {', '.join(automation.get('inactive_expected') or ['none'])}",
            f"- maintenance debt: {debt['score']}/100 budget_remaining={debt['budget_remaining']} status={debt['status']}",
            "",
            "## Recurring Maintenance Hotspots",
            "",
            (
                "\n".join(
                    f"- `{item.get('kind')}` target={item.get('workstream_id') or item.get('title') or item.get('evidence') or 'global'} repeats={item['count']}"
                    for item in recurring
                )
                if recurring
                else "- none in recent dream history"
            ),
            "",
            "## Schema Promotion Candidates",
            "",
            (
                "\n".join(
                    f"- `{item['kind']}` target={item['target']} repeats={item['count']} action={item['proposal']}"
                    for item in schema_candidates
                )
                if schema_candidates
                else "- none ready for promotion"
            ),
            "",
            "## Adaptation Cards",
            "",
            (
                "\n".join(
                    f"- `{item['id']}` repeats={item['repeat_count']} status={item['promotion_status']} validator={item['validator_kind']}"
                    for item in adaptation
                )
                if adaptation
                else "- none above recurrence threshold"
            ),
            "",
            "## Readiness Metrics",
            "",
            json.dumps(readiness_metrics(world, risks, interactions), indent=2, sort_keys=True),
            "",
            "## Morning Checks",
            "",
            "\n".join(f"- `{command}`" for command in checks),
            "",
            "## Dream Boundary",
            "",
            "- Synthetic Justin is not approval.",
            "- No external systems were touched.",
            "- Real register/Linear/daily-plan changes are proposals only.",
            "",
        ]
    )


def render_interactions(interactions: list[dict[str, Any]]) -> str:
    lines = ["# Synthetic Interactions", "", SYNTHETIC_HEADER, ""]
    for item in interactions:
        lines.append(f"## {item['id']} - {item['tested_failure_mode']}")
        lines.append("")
        lines.append(f"- Evidence class: `{item['evidence_class']}`")
        lines.append(f"- Scenario type: `{item['scenario_type']}`")
        lines.append(f"- Validator candidate: `{item['validator_candidate']['kind']}`")
        lines.append("")
        lines.append(item["synthetic_user"])
        lines.append("")
        lines.append(item["dream_response"])
        lines.append("")
    return "\n".join(lines)


def render_prompt_archetype_rehearsals(world: dict[str, Any], risks: list[dict[str, Any]]) -> str:
    lines = ["# Prompt Archetype Rehearsals", "", SYNTHETIC_HEADER, ""]
    lines.append(
        "These are synthetic prompt fixtures. They model recurring request shapes without impersonating Justin or granting approval."
    )
    lines.append("")
    lines.append("## Selection Rule")
    lines.append("")
    lines.append(
        "Nightly dreams rehearse every archetype once, then mark high-priority archetypes when current state contains matching risk signals."
    )
    lines.append("No archetype can mutate real task, repo, SDCC, email, slide, or Linear state.")
    lines.append("")
    for item in prompt_archetype_payload(world, risks):
        lines.append(f"## {item['index']}. `{item['id']}`")
        lines.append("")
        lines.append(f"- Scenario type: `{item['scenario_type']}`")
        lines.append(f"- Selection priority: `{item['selection_priority']}`")
        lines.append(f"- Selection score: {item['selection_score']}")
        lines.append(f"- Synthetic prompt: {item['synthetic_user']}")
        lines.append(f"- Safe first action: {item['safe_first_action']}")
        lines.append(f"- Expected response shape: {item['expected_response_shape']}")
        lines.append(f"- Failure mode prevented: {item['failure_mode_prevented']}")
        lines.append("- Quality checks:")
        for check in item["quality_checks"]:
            lines.append(f"  - {check}")
        seed = item["context_seed"]
        lines.append(
            f"- Context seed: live={seed['live_workstreams']} stale={seed['stale_workstreams']} "
            f"active_jobs={seed['active_jobs']}"
        )
        lines.append(f"- Promotion rule: {item['promotion_rule']}")
        lines.append("")
    lines.append("## Anti-Overfit Guard")
    lines.append("")
    lines.append("- Rotate prompt wording, but keep expected evidence and safety invariants stable.")
    lines.append("- Do not caricature Justin; the simulator stresses evidence, urgency, and scope, not personality.")
    lines.append("- Discard archetypes that produce noise without a validator, runbook, or clearer waking response.")
    return "\n".join(lines)


def render_all_corners_review(world: dict[str, Any], risks: list[dict[str, Any]]) -> str:
    lines = ["# All-Corners Dream Review", "", SYNTHETIC_HEADER, ""]
    lines.append(
        "Every dream iteration must sweep the full thesis-operating surface before proposing priorities."
    )
    lines.append("This file is the detailed review; the heartbeat should remain concise unless action is needed.")
    lines.append("")
    lines.append("## Coverage Table")
    lines.append("")
    for item in all_corner_status(world, risks):
        lines.append(f"- `{item['corner']}`: {item['status']} priority={item['priority']} signal={item['signal']}")
    lines.append("")
    lines.append("## Iteration Contract")
    lines.append("")
    lines.append("- Target 10-20 minutes of local read-only synthesis for hourly heartbeat dreams when state is rich.")
    lines.append("- Keep final heartbeat output short; put rigorous detail in dream artifacts.")
    lines.append("- Do not expand into live SDCC, Gmail, Slides, Drive, Linear, Condor, or repo-tracked mutation.")
    lines.append("- Prefer one exact repair candidate per high-priority corner over broad narrative advice.")
    lines.append("")
    lines.append("## Selection Pressure")
    lines.append("")
    lines.append("- Promote only proposal artifacts with evidence, safe first action, and validator candidate.")
    lines.append("- Reject pseudo-data, paper claims, or synthetic prompts if they are phrased as real project evidence.")
    return "\n".join(lines) + "\n"


def render_literature_scout(world: dict[str, Any]) -> str:
    lines = ["# Literature Scout", "", SYNTHETIC_HEADER, ""]
    lines.append(
        "This is a local source-surface map, not a literature review claim. Verify papers before promoting any thesis text."
    )
    lines.append("")
    lines.append("## Search Boundary")
    lines.append("")
    lines.append("- Local-only scan over useful docs, notes, policy/context files, and checked-out reference code/docs.")
    lines.append("- No network search, no external mutation, no citation claim promoted without waking verification.")
    lines.append("")
    lines.append("## Local Paper / Note Surfaces")
    surfaces = world.get("literature_surfaces") or []
    if not surfaces:
        lines.append("- No local paper/note surfaces matched the scout keywords.")
    for item in surfaces[:25]:
        keywords = ", ".join(item.get("matched_keywords") or [])
        lines.append(f"- `{item['path']}` suffix={item['suffix']} keywords={keywords or 'none'}")
    lines.append("")
    lines.append("## Thesis Figure Questions To Extract")
    questions = [
        "Which observable and normalization make the final physics claim hardest to dismiss?",
        "What control or closure panel should sit next to the headline plot?",
        "Which systematic band or ratio panel carries the credibility burden?",
        "Which binning and axis range are conventional enough to be trusted but sharp enough for sPHENIX?",
        "Which reference result should define the visual grammar: PPG12/18/19, ATLAS gamma-jet, or a local validated note?",
    ]
    for question in questions:
        lines.append(f"- {question}")
    lines.append("")
    lines.append("## Required Waking Verification")
    lines.append("- Open the actual paper/note before using any claim.")
    lines.append("- Record title, source path, page/figure, extracted design lesson, and thesis claim supported.")
    return "\n".join(lines) + "\n"


def render_figure_design_notes() -> str:
    lines = ["# Ideal Final Figure Design Notes", "", SYNTHETIC_HEADER, ""]
    lines.append(
        "These are synthetic target-figure sketches. They are design targets, not data, evidence, signals, or results."
    )
    lines.append("")
    for spec in IDEAL_FIGURE_SPECS:
        lines.append(f"## `{spec['id']}`")
        lines.append(f"- title: {spec['title']}")
        lines.append(f"- thesis claim layer: `{spec['claim']}`")
        lines.append(f"- observable: {spec['observable']}")
        lines.append(f"- x axis: {spec['x_label']}")
        lines.append(f"- y axis: {spec['y_label']}")
        lines.append(f"- required real inputs: {spec['required_inputs']}")
        lines.append("- promotion rule: replace only with real analysis output after provenance, QA, and approval")
        lines.append("")
    return "\n".join(lines)


def render_target_figure_svg(spec: dict[str, str]) -> str:
    title = html.escape(spec["title"])
    x_label = html.escape(spec["x_label"])
    y_label = html.escape(spec["y_label"])
    required = html.escape(spec["required_inputs"])
    polyline_a = "120,310 190,245 260,205 330,182 400,170 470,184 540,214 610,260"
    polyline_b = "120,330 190,285 260,250 330,230 400,218 470,224 540,246 610,288"
    if "closure" in spec["id"]:
        polyline_a = "120,235 190,228 260,232 330,226 400,231 470,229 540,233 610,227"
        polyline_b = "120,260 190,258 260,262 330,255 400,259 470,257 540,261 610,256"
    if "purity" in spec["id"]:
        polyline_a = "120,305 190,275 260,242 330,210 400,184 470,164 540,152 610,146"
        polyline_b = "120,250 190,232 260,214 330,198 400,184 470,172 540,162 610,154"
    return f"""<svg xmlns="http://www.w3.org/2000/svg" width="960" height="540" viewBox="0 0 960 540">
  <rect width="960" height="540" fill="#ffffff"/>
  <text x="50" y="54" font-family="Arial" font-size="25" font-weight="700" fill="#111111">{title}</text>
  <text x="50" y="84" font-family="Arial" font-size="15" font-weight="700" fill="#b00020">SYNTHETIC TARGET FIGURE - NOT DATA</text>
  <line x1="120" y1="370" x2="690" y2="370" stroke="#111111" stroke-width="2"/>
  <line x1="120" y1="120" x2="120" y2="370" stroke="#111111" stroke-width="2"/>
  <line x1="120" y1="245" x2="690" y2="245" stroke="#dddddd" stroke-width="1" stroke-dasharray="6 6"/>
  <polyline points="{polyline_a}" fill="none" stroke="#0055a4" stroke-width="4"/>
  <polyline points="{polyline_b}" fill="none" stroke="#c43b2f" stroke-width="4"/>
  <circle cx="610" cy="260" r="5" fill="#0055a4"/>
  <circle cx="610" cy="288" r="5" fill="#c43b2f"/>
  <text x="350" y="420" font-family="Arial" font-size="18" fill="#111111">{x_label}</text>
  <text transform="translate(64 315) rotate(-90)" font-family="Arial" font-size="18" fill="#111111">{y_label}</text>
  <rect x="730" y="125" width="180" height="84" fill="#f7f7f7" stroke="#bbbbbb"/>
  <line x1="750" y1="155" x2="795" y2="155" stroke="#0055a4" stroke-width="4"/>
  <text x="805" y="161" font-family="Arial" font-size="15" fill="#111111">target signal</text>
  <line x1="750" y1="185" x2="795" y2="185" stroke="#c43b2f" stroke-width="4"/>
  <text x="805" y="191" font-family="Arial" font-size="15" fill="#111111">control/reference</text>
  <text x="730" y="250" font-family="Arial" font-size="14" fill="#111111">Required real inputs:</text>
  <foreignObject x="730" y="262" width="185" height="150">
    <div xmlns="http://www.w3.org/1999/xhtml" style="font-family: Arial; font-size: 12px; color: #333333;">{required}</div>
  </foreignObject>
</svg>
"""


def render_reflection_cards(interactions: list[dict[str, Any]]) -> str:
    lines = ["# Reflection Cards", "", SYNTHETIC_HEADER, ""]
    if not interactions:
        lines.append("- No reflection cards generated.")
        return "\n".join(lines)
    for item in interactions:
        card = item["reflection_card"]
        lines.append(f"## {item['id']} - {card['root_cause_hypothesis']}")
        for key in (
            "trigger",
            "evidence_class",
            "observed_evidence",
            "synthetic_perturbation",
            "validator_candidate",
            "safe_next_diagnostic",
            "validator_result",
            "promotion_status",
            "expiry",
        ):
            lines.append(f"- {key}: {card[key]}")
        lines.append("")
    return "\n".join(lines)


def render_runbook_proposals(risks: list[dict[str, Any]]) -> str:
    lines = ["# Runbook Proposals", "", SYNTHETIC_HEADER, ""]
    proposals = [runbook_proposal_for(item) for item in risks if item["score"] >= 60][:8]
    if not proposals:
        lines.append("- No runbook proposal generated.")
        return "\n".join(lines)
    for index, proposal in enumerate(proposals, start=1):
        lines.append(f"## Proposal {index}: {proposal['purpose']}")
        for key in (
            "scope",
            "scenario_type",
            "safe_command",
            "forbidden_actions",
            "rollback_path",
            "last_validation",
        ):
            lines.append(f"- {key}: {proposal[key]}")
        lines.append("")
    return "\n".join(lines)


def render_quality_control(interactions: list[dict[str, Any]]) -> str:
    rows = []
    for item in interactions:
        rows.append(
            {
                "id": item["id"],
                "failure_mode": item["tested_failure_mode"],
                "quality_gate": item["quality_gate"],
                "outcome": item["outcome"],
            }
        )
    return json.dumps(
        {
            "synthetic": True,
            "synthetic_header": SYNTHETIC_HEADER,
            "evidence_classes": EVIDENCE_CLASSES,
            "scenario_portfolio": SCENARIO_PORTFOLIO,
            "rows": rows,
        },
        indent=2,
        sort_keys=True,
    ) + "\n"


def render_daily_proposals(risks: list[dict[str, Any]]) -> str:
    lines = ["# Daily Plan Proposals", "", SYNTHETIC_HEADER, ""]
    promoted = [item for item in risks if item["score"] >= 70][:3]
    if promoted:
        lines.append("## Proposed Top Attention")
        for item in promoted:
            target = item.get("workstream_id") or item.get("evidence") or "global"
            lines.append(f"- `{target}`: {item['label']} (score {item['score']})")
    else:
        lines.append("- No priority promotion proposed.")
    lines.append("")
    lines.append("Do not apply automatically. Use this as morning review input.")
    return "\n".join(lines)


def render_linear_proposals(risks: list[dict[str, Any]]) -> str:
    lines = ["# Linear Update Proposals", "", SYNTHETIC_HEADER, ""]
    actionable = [item for item in risks if item.get("workstream_id") and item["score"] >= 60]
    if not actionable:
        lines.append("- No Linear proposal generated.")
    for item in actionable[:6]:
        lines.append(f"- `{item['workstream_id']}`: add comment noting `{item['kind']}` risk and requested waking check.")
    lines.append("")
    lines.append("Dreams must not write Linear directly.")
    return "\n".join(lines)


def render_cleanup_proposals(world: dict[str, Any], risks: list[dict[str, Any]]) -> str:
    lines = ["# Cleanup Proposals", "", SYNTHETIC_HEADER, ""]
    lines.append("Dream cleanup is proposal-only. Do not delete, move, archive, or mutate SDCC from this file.")
    lines.append("")
    lines.append("## Local Memory State")
    for label, summary in sorted((world.get("local_state") or {}).items()):
        lines.append(
            f"- `{label}` `{summary.get('path')}`: dirs={summary.get('directories')} "
            f"bytes={summary.get('bytes')} oldest={summary.get('oldest') or 'none'} "
            f"newest={summary.get('newest') or 'none'}"
        )
        for candidate in summary.get("cleanup_candidates") or []:
            lines.append(
                f"  candidate: `{candidate['path']}` age_days={candidate['age_days']} bytes={candidate['bytes']} "
                f"reason={candidate['reason']}"
            )
    lines.append("")
    lines.append("## Recorded SDCC/Storage Signals")
    storage_hits = world.get("storage_signal_hits") or []
    if storage_hits:
        for hit in storage_hits[:8]:
            lines.append(f"- `{hit['path']}:{hit['line']}` {hit['snippet']}")
    else:
        lines.append("- No SDCC/storage clutter signal found in local recorded context.")
    lines.append("")
    lines.append("## Proposed Waking Review")
    cleanup_risks = [item for item in risks if item["kind"] in {"memory_hygiene", "sdcc_clutter_watch"}]
    if cleanup_risks:
        for item in cleanup_risks:
            validator = validator_candidate_for(item)
            lines.append(f"- `{item['kind']}`: {item['label']}")
            lines.append(f"  safe check: `{validator['safe_command']}`")
            lines.append("  approval needed before any real cleanup")
    else:
        lines.append("- No cleanup review promoted by this dream.")
    return "\n".join(lines) + "\n"


def render_maintenance_debt(world: dict[str, Any], risks: list[dict[str, Any]]) -> str:
    debt = maintenance_debt(world, risks)
    schemas = schema_promotion_candidates(world)
    cards = adaptation_cards(world, risks)
    lines = ["# Maintenance Debt Ledger", "", SYNTHETIC_HEADER, ""]
    lines.append(
        "This ledger is a proposal-only reliability budget for the thesis OS. It ranks debt that should reduce dream/doctor feature growth until it is controlled."
    )
    lines.append("")
    lines.append(f"- debt_score: {debt['score']}/100")
    lines.append(f"- budget_remaining: {debt['budget_remaining']}/100")
    lines.append(f"- status: {debt['status']}")
    lines.append("")
    lines.append("## Components")
    for key, value in debt["components"].items():
        lines.append(f"- `{key}`: {value}")
    lines.append("")
    lines.append("## Growth Rule")
    lines.append("- `healthy`: dream/doctor can add small new checks if validation stays clean.")
    lines.append("- `strained`: prioritize debt reduction and boundary precision over new feature breadth.")
    lines.append("- `frozen_growth`: stop adding new dream capability until boundary drift or heartbeat debt is reduced.")
    lines.append("")
    lines.append("## Promotion Queue")
    if schemas:
        for item in schemas:
            lines.append(
                f"- `{item['kind']}` target={item['target']} repeats={item['count']}: {item['proposal']}"
            )
    else:
        lines.append("- no recurring hotspot has crossed the schema-promotion threshold yet")
    lines.append("")
    lines.append("## Adaptation Card Pressure")
    if cards:
        for item in cards:
            lines.append(
                f"- `{item['id']}` repeats={item['repeat_count']} status={item['promotion_status']} validator={item['validator_kind']}"
            )
    else:
        lines.append("- no recurring hotspot has crossed the adaptation-card threshold yet")
    return "\n".join(lines) + "\n"


def render_adaptation_cards(world: dict[str, Any], risks: list[dict[str, Any]]) -> str:
    cards = adaptation_cards(world, risks)
    lines = ["# Adaptation Cards", "", SYNTHETIC_HEADER, ""]
    lines.append(
        "These proposal-only cards consolidate repeated dream findings into validator-backed waking review candidates."
    )
    lines.append("They may adjust attention and review priority, but they do not mutate code, memory, tasks, or external systems.")
    lines.append("")
    if not cards:
        lines.append("- no repeated dream finding has crossed the adaptation-card threshold")
        return "\n".join(lines) + "\n"
    for item in cards:
        lines.append(f"## `{item['id']}`")
        lines.append(f"- kind: `{item['kind']}`")
        lines.append(f"- target: {item['target']}")
        lines.append(f"- repeat_count: {item['repeat_count']}")
        lines.append(f"- latest_run_id: {item.get('latest_run_id') or 'unknown'}")
        lines.append(f"- evidence: {item['evidence']}")
        lines.append(f"- validator_kind: `{item['validator_kind']}`")
        lines.append(f"- safe_command: `{item['safe_command']}`")
        lines.append(f"- expected_signal: {item['expected_signal']}")
        lines.append(f"- promotion_status: `{item['promotion_status']}`")
        lines.append(f"- retention_rule: {item['retention_rule']}")
        lines.append(f"- decay_rule: {item['decay_rule']}")
        lines.append(f"- promotion_rule: {item['promotion_rule']}")
        lines.append(f"- biological_analogy: {item['biological_analogy']}")
        lines.append(f"- mutation_boundary: `{item['mutation_boundary']}`")
        lines.append("")
    return "\n".join(lines) + "\n"


def render_overnight_hygiene_proposals(risks: list[dict[str, Any]]) -> str:
    lines = ["# Overnight Hygiene Proposals", "", SYNTHETIC_HEADER, ""]
    lines.append("Allowed overnight work must improve readiness without changing analysis state.")
    lines.append("")
    lines.append("## Safe Local Tasks")
    for command in (
        "python3 scripts/codex_os_doctor.py --profile daily",
        "python3 scripts/codex_work_register_stale.py agent_context/CODEX_WORK_REGISTER.yaml",
        "python3 scripts/codex_thesis_radar.py",
        "python3 scripts/codex_artifact_registry.py check",
        "python3 scripts/codex_os_dream.py validate --latest",
    ):
        lines.append(f"- `{command}`")
    lines.append("")
    lines.append("## Forbidden Overnight Actions")
    for item in (
        "Condor submission, removal, merge rerun, or production analysis",
        "SDCC deletion, transfer, remote edit, or queue control",
        "Gmail, Drive, Slides, Calendar, Linear, or repo-tracked mutation",
        "Automatic memory deletion or task promotion",
    ):
        lines.append(f"- {item}")
    relevant = [item for item in risks if item["kind"] == "safe_overnight_hygiene"]
    if relevant:
        lines.append("")
        lines.append("## Dream Finding")
        for item in relevant:
            lines.append(f"- {item['label']}: {item['evidence']}")
    return "\n".join(lines) + "\n"


def render_physics_scenario_proposals(world: dict[str, Any], risks: list[dict[str, Any]]) -> str:
    lines = ["# Physics Scenario Proposals", "", SYNTHETIC_HEADER, ""]
    lines.append(
        "These are hypothesis scouts only. They are not results, signals, discoveries, or approval to run analysis."
    )
    lines.append("")
    surfaces = world.get("physics_surfaces") or {}
    claims = surfaces.get("claims_with_artifacts") or []
    topics = surfaces.get("topic_tokens") or []
    lines.append("## Usable Data Surfaces Seen")
    lines.append(f"- registered physics claim layers with artifacts: {', '.join(claims) if claims else 'none'}")
    lines.append(f"- active topic tokens: {', '.join(topics) if topics else 'none'}")
    lines.append("")
    if not claims:
        lines.append("No physics scenario proposed because no usable registered physics artifact surface was found.")
        return "\n".join(lines) + "\n"
    lines.append("## High-Upside Hypothesis Templates For Waking Review")
    proposals = [
        {
            "hypothesis": "A nontrivial gamma-jet modification pattern appears only after pp baseline and embedded-background closure are both provenance-clean.",
            "signature": "shape change persists under baseline/stitching/closure variations rather than appearing as a normalization artifact",
            "null_test": "repeat on validated pp or embedding closure surface where no medium modification should appear",
            "systematic": "photon purity, background stitching, response/unfolding, and stale artifact contamination",
        },
        {
            "hypothesis": "ML photon-ID or isolation behavior exposes a control-region structure that improves final Au+Au gamma-jet sensitivity.",
            "signature": "robust separation or purity stability across centrality/kinematic/control variations",
            "null_test": "compare against PPG12-style baseline selection and known-safe pp/reference samples",
            "systematic": "domain shift, training leakage, sample stitching, and selection sculpting",
        },
        {
            "hypothesis": "A surprising xJgamma or recoil response feature survives all known detector/background controls.",
            "signature": "localized, reproducible deviation tied to physics variables and not file/sample/provenance boundaries",
            "null_test": "sideband, closure, shuffled labels, alternate binning, and independent production cross-check",
            "systematic": "unfolding regularization, response mismodeling, centrality bias, and trigger/sample boundaries",
        },
    ]
    for index, proposal in enumerate(proposals, start=1):
        lines.append(f"### Scenario {index}")
        lines.append(f"- hypothesis: {proposal['hypothesis']}")
        lines.append(f"- expected qualitative signature: {proposal['signature']}")
        lines.append("- required input provenance: `python3 scripts/codex_artifact_registry.py check`")
        lines.append(f"- null/control check: {proposal['null_test']}")
        lines.append(f"- dominant systematic risk: {proposal['systematic']}")
        lines.append("- minimal safe next diagnostic: waking review of registered artifacts and existing plots only")
        lines.append("- promotion rule: no run, plot campaign, or claim without explicit waking approval")
        lines.append("")
    return "\n".join(lines)


def proposed_diff(path: Path, proposal_lines: list[str]) -> str:
    original = read_text(path)
    proposed = original
    if proposal_lines:
        proposed += "\n\n# DREAM_PROPOSAL_ONLY - do not apply without waking review\n"
        proposed += "\n".join(f"# {line}" for line in proposal_lines)
        proposed += "\n"
    return "".join(
        difflib.unified_diff(
            original.splitlines(keepends=True),
            proposed.splitlines(keepends=True),
            fromfile=f"a/{path.as_posix()}",
            tofile=f"b/{path.as_posix()}",
        )
    )


def copy_sandbox(run_dir: Path) -> None:
    for path in [
        DEFAULT_REGISTER,
        ARTIFACT_REGISTRY,
        THESIS_MAP,
        Path("agent_context/policies/AGENTIC_OS_HARDENING.md"),
        Path("agent_context/policies/AGENTIC_OS_DREAMING.md"),
    ]:
        if path.exists():
            destination = run_dir / "sandbox" / path
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, destination)


def write_dream_outputs(run_dir: Path, run_id: str, mode: str, world: dict[str, Any], risks: list[dict[str, Any]], interactions: list[dict[str, Any]]) -> None:
    copy_sandbox(run_dir)
    prompt_archetypes = prompt_archetype_payload(world, risks)
    signal = heartbeat_signal(run_id, mode, world, risks)
    cards = adaptation_cards(world, risks)
    trace = {
        "synthetic": True,
        "synthetic_header": SYNTHETIC_HEADER,
        "mutation_boundary": "proposal_only",
        "evidence_classes": EVIDENCE_CLASSES,
        "scenario_portfolio": SCENARIO_PORTFOLIO,
        "prompt_archetypes": prompt_archetypes,
        "run_id": run_id,
        "mode": mode,
        "world": {
            "workstream_counts": world["workstream_counts"],
            "artifact_counts": world["artifact_counts"],
            "artifact_gaps": world["artifact_gaps"],
            "local_state": world["local_state"],
            "storage_signal_hits": world["storage_signal_hits"],
            "physics_surfaces": world["physics_surfaces"],
            "literature_surface_count": len(world["literature_surfaces"]),
            "all_corners": all_corner_status(world, risks),
            "recent_event_count": len(world["recent_events"]),
            "search_hit_count": len(world["search_hits"]),
        },
        "risks": risks,
        "interactions": interactions,
        "adaptation_cards": cards,
        "readiness_metrics": readiness_metrics(world, risks, interactions),
        "search_hits": world["search_hits"],
        "heartbeat_signal": signal,
    }
    safe_write(run_dir, "dream_trace.json", json.dumps(trace, indent=2, sort_keys=True) + "\n")
    safe_write(run_dir, "heartbeat_signal.json", json.dumps(signal, indent=2, sort_keys=True) + "\n")
    safe_write(run_dir, "dream_report.md", render_report(run_id, mode, world, risks, interactions))
    safe_write(run_dir, "synthetic_interactions.md", render_interactions(interactions))
    safe_write(run_dir, "prompt_archetype_rehearsals.md", render_prompt_archetype_rehearsals(world, risks))
    safe_write(run_dir, "all_corners_review.md", render_all_corners_review(world, risks))
    safe_write(run_dir, "literature_scout.md", render_literature_scout(world))
    safe_write(run_dir, "figure_design_notes.md", render_figure_design_notes())
    for spec in IDEAL_FIGURE_SPECS:
        safe_write(run_dir, f"ideal_final_figure_gallery/{spec['id']}.svg", render_target_figure_svg(spec))
    safe_write(run_dir, "reflection_cards.md", render_reflection_cards(interactions))
    safe_write(run_dir, "runbook_proposals.md", render_runbook_proposals(risks))
    safe_write(run_dir, "quality_control.json", render_quality_control(interactions))
    safe_write(run_dir, "daily_plan_proposals.md", render_daily_proposals(risks))
    safe_write(run_dir, "linear_update_proposals.md", render_linear_proposals(risks))
    safe_write(run_dir, "cleanup_proposals.md", render_cleanup_proposals(world, risks))
    safe_write(run_dir, "maintenance_debt.md", render_maintenance_debt(world, risks))
    safe_write(run_dir, "adaptation_cards.md", render_adaptation_cards(world, risks))
    safe_write(run_dir, "overnight_hygiene_proposals.md", render_overnight_hygiene_proposals(risks))
    safe_write(run_dir, "physics_scenario_proposals.md", render_physics_scenario_proposals(world, risks))

    top_lines = [
        f"{item['kind']} risk for {item.get('workstream_id') or item.get('evidence') or 'global'} score={item['score']}"
        for item in risks[:5]
    ]
    safe_write(run_dir, "register_patch.diff", proposed_diff(DEFAULT_REGISTER, top_lines))
    policy_lines = [
        f"Consider encoding dream prevention for {item['kind']} if it repeats: {item['label']}"
        for item in risks
        if item["score"] >= 80
    ][:3]
    safe_write(run_dir, "policy_patch.diff", proposed_diff(Path("agent_context/policies/AGENTIC_OS_DREAMING.md"), policy_lines))
    append_dream_index(signal)


def run_dream(args: argparse.Namespace) -> int:
    now = parse_now(args.now)
    run_id = args.run_id or run_id_for(args.mode, now)
    run_dir = DREAM_ROOT / run_id
    register_data = load_yaml_like(Path(args.register))
    artifact_data = load_yaml_like(ARTIFACT_REGISTRY)
    world = build_world(register_data, artifact_data, now, args.mode)
    risks = identify_risks(world, args.mode)
    interactions = generate_interactions(risks, args.mode)
    write_dream_outputs(run_dir, run_id, args.mode, world, risks, interactions)
    print(f"OK: {args.mode} dream wrote {run_dir}")
    print(f"risks={len(risks)} interactions={len(interactions)}")
    print(f"report={run_dir / 'dream_report.md'}")
    return 0


def latest_dream_dir() -> Path | None:
    if not DREAM_ROOT.exists():
        return None
    candidates = [path for path in DREAM_ROOT.iterdir() if path.is_dir()]
    if not candidates:
        return None
    return max(candidates, key=lambda path: path.stat().st_mtime)


def validate_dream_dir(path: Path) -> list[str]:
    errors: list[str] = []
    root = DREAM_ROOT.resolve()
    resolved = path.resolve()
    if root not in [resolved, *resolved.parents]:
        errors.append(f"dream path is outside dream root: {path}")
        return errors
    required = [
        "dream_report.md",
        "dream_trace.json",
        "heartbeat_signal.json",
        "synthetic_interactions.md",
        "prompt_archetype_rehearsals.md",
        "all_corners_review.md",
        "literature_scout.md",
        "figure_design_notes.md",
        "reflection_cards.md",
        "runbook_proposals.md",
        "quality_control.json",
        "register_patch.diff",
        "policy_patch.diff",
        "linear_update_proposals.md",
        "daily_plan_proposals.md",
        "cleanup_proposals.md",
        "maintenance_debt.md",
        "adaptation_cards.md",
        "overnight_hygiene_proposals.md",
        "physics_scenario_proposals.md",
    ]
    for name in required:
        if not (path / name).exists():
            errors.append(f"missing dream output: {name}")

    report = read_text(path / "dream_report.md")
    interactions = read_text(path / "synthetic_interactions.md")
    prompt_archetypes = read_text(path / "prompt_archetype_rehearsals.md")
    if SYNTHETIC_HEADER not in report:
        errors.append("dream_report.md lacks synthetic provenance header")
    if SYNTHETIC_HEADER not in interactions or "SYNTHETIC Justin" not in interactions:
        errors.append("synthetic_interactions.md lacks synthetic Justin marker")
    if SYNTHETIC_HEADER not in prompt_archetypes or "meta_os_upgrade" not in prompt_archetypes:
        errors.append("prompt_archetype_rehearsals.md lacks synthetic header or meta_os_upgrade archetype")
    for name in (
        "prompt_archetype_rehearsals.md",
        "all_corners_review.md",
        "literature_scout.md",
        "figure_design_notes.md",
        "reflection_cards.md",
        "runbook_proposals.md",
        "quality_control.json",
        "cleanup_proposals.md",
        "maintenance_debt.md",
        "adaptation_cards.md",
        "overnight_hygiene_proposals.md",
        "physics_scenario_proposals.md",
    ):
        if SYNTHETIC_HEADER not in read_text(path / name):
            errors.append(f"{name} lacks synthetic provenance header")
    for name in required:
        text = read_text(path / name)
        if APPROVAL_LIKE_PATTERN.search(text):
            errors.append(f"{name} contains approval-like language")

    try:
        trace = json.loads(read_text(path / "dream_trace.json"))
    except json.JSONDecodeError as exc:
        errors.append(f"dream_trace.json invalid JSON: {exc}")
        return errors
    try:
        heartbeat = json.loads(read_text(path / "heartbeat_signal.json"))
    except json.JSONDecodeError as exc:
        errors.append(f"heartbeat_signal.json invalid JSON: {exc}")
        return errors
    if trace.get("synthetic") is not True:
        errors.append("dream_trace.json synthetic flag is not true")
    if trace.get("mutation_boundary") != "proposal_only":
        errors.append("dream_trace.json mutation_boundary is not proposal_only")
    if trace.get("evidence_classes") != EVIDENCE_CLASSES:
        errors.append("dream_trace.json evidence_classes are missing or changed")
    prompt_payload = trace.get("prompt_archetypes") or []
    if len(prompt_payload) < len(PROMPT_ARCHETYPES):
        errors.append("dream_trace.json missing prompt archetype payload")
    if not any(isinstance(item, dict) and item.get("id") == "meta_os_upgrade" for item in prompt_payload):
        errors.append("dream_trace.json missing meta_os_upgrade prompt archetype")
    all_corners = trace.get("world", {}).get("all_corners") or []
    if len(all_corners) < len(ALL_CORNERS):
        errors.append("dream_trace.json missing all-corners coverage")
    if not trace.get("readiness_metrics"):
        errors.append("dream_trace.json missing readiness_metrics")
    if heartbeat.get("synthetic") is not True:
        errors.append("heartbeat_signal.json synthetic flag is not true")
    if heartbeat.get("mutation_boundary") != "proposal_only":
        errors.append("heartbeat_signal.json mutation_boundary is not proposal_only")
    if not heartbeat.get("top_findings"):
        errors.append("heartbeat_signal.json missing top_findings")
    if not heartbeat.get("morning_checks"):
        errors.append("heartbeat_signal.json missing morning_checks")
    gallery = path / "ideal_final_figure_gallery"
    if not gallery.is_dir():
        errors.append("missing ideal_final_figure_gallery directory")
    else:
        for spec in IDEAL_FIGURE_SPECS:
            svg = read_text(gallery / f"{spec['id']}.svg")
            if "SYNTHETIC TARGET FIGURE - NOT DATA" not in svg:
                errors.append(f"synthetic target figure lacks non-data label: {spec['id']}.svg")
    for item in trace.get("interactions") or []:
        if not isinstance(item, dict) or item.get("synthetic") is not True:
            errors.append("interaction missing synthetic=true")
        if "SYNTHETIC Justin" not in str(item.get("synthetic_user")):
            errors.append("interaction missing SYNTHETIC Justin marker")
    return errors


def command_validate(args: argparse.Namespace) -> int:
    path = latest_dream_dir() if args.latest else Path(args.path)
    if path is None:
        if args.allow_missing:
            print("OK: no dream directory exists yet")
            return 0
        print("ERROR: no dream directory exists", file=sys.stderr)
        return 1
    errors = validate_dream_dir(path)
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(f"OK: dream validates at {path}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    for mode in ("micro", "nightly"):
        sub = subparsers.add_parser(mode, help=f"run a {mode} synthetic dream")
        sub.add_argument("register", nargs="?", default=str(DEFAULT_REGISTER))
        sub.add_argument("--now")
        sub.add_argument("--run-id")
        sub.set_defaults(func=run_dream, mode=mode)

    validate = subparsers.add_parser("validate", help="validate dream output safety")
    validate.add_argument("path", nargs="?", default="")
    validate.add_argument("--latest", action="store_true")
    validate.add_argument("--allow-missing", action="store_true")
    validate.set_defaults(func=command_validate)

    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
