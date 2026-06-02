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
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from codex_work_register_common import DEFAULT_REGISTER, first_line, load_register, parse_when, sorted_workstreams
from codex_context_resonance import build_context_resonance_payload, render_markdown as render_context_resonance_markdown
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
SHADOW_EVIDENCE_CLASSES = ["real_observed", "human_approved", "derived", "synthetic"]
FUEL_CLASSES = ["LIVE", "DORMANT", "DUPLICATE", "STALE", "CONFLICTED", "SYNTHETIC", "UNOWNED", "PROTECTED"]
SHADOW_PILOT_STAGES = {
    1: {
        "name": "Night 1 - Baseline Map",
        "status": "baseline_map_only",
        "output_level": "inventory_and_scores_only",
        "description": "Build auditable maps and branch-pressure baselines; do not create approval-ready packets.",
    },
    2: {
        "name": "Night 2 - Candidate Formation",
        "status": "candidate_formation_only",
        "output_level": "candidate_hypotheses_only",
        "description": "Draft schemas, validators, quarantines, compression ideas, and ranked maintenance hypotheses.",
    },
    3: {
        "name": "Night 3 - Shadow Burn Rehearsal",
        "status": "draft_shadow_only",
        "output_level": "shadow_rehearsal_only",
        "description": "Create one sandbox-only controlled-burn rehearsal with rollback and canary sketches.",
    },
    4: {
        "name": "Night 4 - Calibration And Human Gate",
        "status": "calibration_waiting_human_gate",
        "output_level": "calibration_and_review_only",
        "description": "Compare prior shadow nights against waking reality and wait for Justin before approval-ready packets.",
    },
}
TFE_BENEFIT_WEIGHTS = {
    "progress": 1.0,
    "correctness": 1.1,
    "reuse": 0.8,
    "safety": 1.1,
    "human_friction_reduction": 0.9,
    "future_tractability": 1.0,
}
TFE_COST_WEIGHTS = {
    "token_cost": 0.7,
    "lookup_cost": 0.8,
    "operational_risk": 1.2,
    "fragmentation": 0.9,
    "maintenance_toil": 1.0,
    "approval_burden": 0.7,
    "synthetic_contamination_risk": 1.3,
}
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

AUTONOMY_TIERS = {
    "auto_safe": "local generated-junk cleanup and ignored internal artifact refresh",
    "auto_validated": "small tracked Codex-OS maintenance packages with rollback and validation evidence",
    "research_only": "sanitized external-model research prompts and local synthesis, never authority",
    "blocked_for_waking": "runtime, science, SDCC mutation, external apps, or task/status changes",
}
DREAM_INTERNAL_NOISE_WORKSTREAM_IDS = {
    "hp26_photon_id_talk",
}
HANDLED_RECURRENCE_PROTOCOLS = {
    "active_job_status": "register_workstream_refresh_contract",
    "active_job_evidence": "register_workstream_refresh_contract",
    "stale_state": "register_workstream_refresh_contract",
    "wip_overload": "morning_cockpit_compression_contract",
}
DREAM_LANE_DEFINITIONS = [
    {
        "lane_id": "status_provenance",
        "purpose": "refresh real evidence, stale-state, active-job, and artifact-provenance pressure",
        "agent_role": "evidence auditor",
        "parallelizable": True,
        "primary_outputs": ["targeted_findings.json", "runbook_proposals.md"],
        "progression": "observe status drift -> classify current/waiting/stale -> propose reusable refresh validators",
    },
    {
        "lane_id": "architecture_cohesion",
        "purpose": "convert repeated warnings into cross-system validators, runbooks, indexes, or compression rules",
        "agent_role": "systems architect",
        "parallelizable": True,
        "primary_outputs": ["structural_advancements.md", "internal_evolution_queue.md"],
        "progression": "measure friction -> compare against simple alternatives -> propose the smallest durable structure",
    },
    {
        "lane_id": "context_resonance",
        "purpose": "improve retrieval quality through latent nudges, negative memories, stale-context suppression, and reconsolidation review",
        "agent_role": "context resonance auditor",
        "parallelizable": True,
        "primary_outputs": [
            "context_resonance_review.md",
            "latent_context_nudges.md",
            "negative_memory_candidates.md",
            "retrieval_outcome_candidates.md",
        ],
        "progression": "surface small useful context -> suppress stale/synthetic context -> record retrieval outcomes for future homeostasis",
    },
    {
        "lane_id": "cleanup_storage",
        "purpose": "find local clutter, stale internal artifacts, retained dream/research packs, and recorded SDCC/storage drift",
        "agent_role": "hygiene auditor",
        "parallelizable": True,
        "primary_outputs": ["cleanup_proposals.md", "sdcc_base_repo_hygiene.md", "changed_actions.md"],
        "progression": "remove only auto-safe generated junk -> list review candidates -> keep SDCC read-only",
    },
    {
        "lane_id": "path_contract",
        "purpose": "guard local, SDCC, SFTP, TG bulk, and output-default path contracts",
        "agent_role": "interface contract verifier",
        "parallelizable": True,
        "primary_outputs": ["path_contract_drift.md", "sdcc_base_repo_hygiene.md"],
        "progression": "detect legacy drift -> map canonical/compat paths -> propose resolver or index fixes",
    },
    {
        "lane_id": "research_scout",
        "purpose": "use ChatGPT as the default online critique/source-lead worker and improve how Codex uses human-accessible AI/tools",
        "agent_role": "external research delegator",
        "parallelizable": True,
        "primary_outputs": ["research_synthesis.md", "human_tool_leverage.md"],
        "progression": "queue sanitized prompt -> review history/tool leverage opportunities -> collect/paste when appropriate -> verify locally before promotion",
    },
    {
        "lane_id": "science_scout",
        "purpose": "draft high-upside physics hypotheses only from provenance-backed existing surfaces",
        "agent_role": "physics hypothesis scout",
        "parallelizable": True,
        "primary_outputs": ["physics_scenario_proposals.md", "figure_design_notes.md"],
        "progression": "identify usable surfaces -> require null tests/systematics -> never claim results from dreams",
    },
    {
        "lane_id": "presentation_artifacts",
        "purpose": "handle slide/figure implications only when the active thesis work needs them",
        "agent_role": "presentation artifact scout",
        "parallelizable": True,
        "primary_outputs": ["figure_design_notes.md", "ideal_final_figure_gallery/"],
        "progression": "generate synthetic targets or local PNG candidates only; do not mutate Google Slides",
        "suppression_rule": "talk-only workstreams such as hp26_photon_id_talk do not drive thesisAnalysis internal-maintenance pressure",
    },
]

DREAM_LANE_AUTOMATIONS = {
    "status_provenance": {
        "automation_id": "thesisanalysis-dream-status-provenance",
        "thread_binding": "fresh_chat_per_run",
    },
    "architecture_cohesion": {
        "automation_id": "thesisanalysis-dream-architecture-cohesion",
        "thread_binding": "fresh_chat_per_run",
    },
    "context_resonance": {
        "automation_id": "thesisanalysis-dream-lane-context-resonance",
        "thread_binding": "fresh_chat_per_run",
    },
    "cleanup_storage": {
        "automation_id": "thesisanalysis-dream-cleanup-storage",
        "thread_binding": "fresh_chat_per_run",
    },
    "path_contract": {
        "automation_id": "thesisanalysis-dream-path-contract",
        "thread_binding": "fresh_chat_per_run",
    },
    "research_scout": {
        "automation_id": "thesisanalysis-dream-research-scout",
        "thread_binding": "fresh_chat_per_run",
    },
    "science_scout": {
        "automation_id": "thesisanalysis-dream-science-scout",
        "thread_binding": "fresh_chat_per_run",
    },
    "presentation_artifacts": {
        "automation_id": "thesisanalysis-dream-presentation-artifacts",
        "thread_binding": "fresh_chat_per_run",
    },
}

EXPECTED_DREAM_AUTOMATION_IDS = [
    DREAM_LANE_AUTOMATIONS[item["lane_id"]]["automation_id"] for item in DREAM_LANE_DEFINITIONS
]
OPTIONAL_DREAM_AUTOMATION_IDS = ["thesisanalysis-morning-priorities"]
LEGACY_DREAM_AUTOMATION_IDS = [
    "thesisanalysis-nightly-dream",
    "thesisanalysis-os-doctor",
    "thesisanalysis-micro-dream",
]
DREAM_AUTOMATION_IDS = tuple(
    EXPECTED_DREAM_AUTOMATION_IDS + OPTIONAL_DREAM_AUTOMATION_IDS + LEGACY_DREAM_AUTOMATION_IDS
)
DREAM_LARGE_MODES = {"nightly", "lane"}
AUTO_SAFE_DELETE_NAMES = {".DS_Store"}
AUTO_SAFE_DELETE_SUFFIXES = {".pyc", ".pyo"}
AUTO_SAFE_DELETE_DIR_NAMES = {"__pycache__"}
AUTO_SAFE_DELETE_PREFIXES = {"._"}
AUTO_SAFE_ROOTS = [
    Path("agent_context/local"),
    Path("scripts"),
    Path("macros"),
    Path("codex_notes"),
]
AUTONOMY_VALIDATION_COMMANDS = [
    "python3 scripts/codex_os_dream.py validate --latest",
    "python3 scripts/codex_os_dream.py lane --lane-id status_provenance",
    "python3 scripts/codex_os_doctor.py --profile daily",
    "python3 scripts/os/register/codex_work_register_stale.py agent_context/CODEX_WORK_REGISTER.yaml --protocol",
    "python3 scripts/os/register/codex_work_register_validate.py agent_context/CODEX_WORK_REGISTER.yaml",
    "python3 scripts/os/artifacts/codex_artifact_registry.py check",
    "python3 scripts/sdcc/validate_transfer_map.py",
    "python3 scripts/sdcc/runtime/audit/checkout_hygiene_probe.py . --profile local --report-only",
    "git diff --check",
]


def now_utc() -> datetime:
    return datetime.now(timezone.utc)


def run_id_for(mode: str, now: datetime) -> str:
    return f"{now.strftime('%Y%m%dT%H%M%SZ')}-{mode}"


def lane_run_id_for(lane_id: str, now: datetime) -> str:
    return f"{now.strftime('%Y%m%dT%H%M%SZ')}-lane-{lane_id}"


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


def repo_root() -> Path:
    return Path(".").resolve()


def repo_relative(path: Path) -> str:
    root = repo_root()
    try:
        return path.resolve().relative_to(root).as_posix()
    except ValueError:
        return path.as_posix()


def git_tracked_paths() -> set[str]:
    try:
        result = subprocess.run(["git", "ls-files", "-z"], text=True, capture_output=True, check=False)
    except OSError:
        return set()
    if result.returncode != 0:
        return set()
    return {item for item in result.stdout.split("\0") if item}


def is_auto_safe_generated_trash(path: Path) -> bool:
    name = path.name
    if name in AUTO_SAFE_DELETE_NAMES:
        return True
    if any(name.startswith(prefix) for prefix in AUTO_SAFE_DELETE_PREFIXES):
        return True
    if path.is_dir() and name in AUTO_SAFE_DELETE_DIR_NAMES:
        return True
    if path.is_file() and path.suffix in AUTO_SAFE_DELETE_SUFFIXES:
        return True
    return False


def collect_auto_safe_generated_trash(run_dir: Path) -> list[Path]:
    root = repo_root()
    tracked = git_tracked_paths()
    candidates: list[Path] = []
    roots = [root, *(root / item for item in AUTO_SAFE_ROOTS)]
    seen_roots: set[Path] = set()
    run_resolved = run_dir.resolve()
    for start in roots:
        if not start.exists():
            continue
        resolved_start = start.resolve()
        if resolved_start in seen_roots:
            continue
        seen_roots.add(resolved_start)
        if resolved_start == root:
            iterable = list(resolved_start.iterdir())
        else:
            iterable = list(resolved_start.rglob("*"))
        for path in iterable:
            try:
                resolved = path.resolve()
            except OSError:
                continue
            if root not in [resolved, *resolved.parents]:
                continue
            if (root / ".git").resolve() in [resolved, *resolved.parents]:
                continue
            if run_resolved in [resolved, *resolved.parents]:
                continue
            if not is_auto_safe_generated_trash(path):
                continue
            rel = repo_relative(path)
            if rel in tracked:
                continue
            candidates.append(path)
    candidates.sort(key=lambda item: (len(item.parts), item.as_posix()), reverse=True)
    return candidates


def remove_auto_safe_generated_trash(candidates: list[Path], limit: int = 200) -> list[dict[str, Any]]:
    actions: list[dict[str, Any]] = []
    for path in candidates[:limit]:
        before_bytes = 0
        rel = repo_relative(path)
        try:
            if path.is_dir():
                before_bytes = sum(child.stat().st_size for child in path.rglob("*") if child.is_file())
                shutil.rmtree(path)
                action_type = "removed_generated_cache_dir"
            elif path.is_file():
                before_bytes = path.stat().st_size
                path.unlink()
                action_type = "removed_generated_cache_file"
            else:
                continue
        except OSError as exc:
            actions.append(
                {
                    "tier": "auto_safe",
                    "performed": False,
                    "action": "failed_generated_junk_cleanup",
                    "path": rel,
                    "reason": str(exc),
                    "rollback": "nothing changed for this path",
                }
            )
            continue
        actions.append(
            {
                "tier": "auto_safe",
                "performed": True,
                "action": action_type,
                "path": rel,
                "before_bytes": before_bytes,
                "reason": "generated local cache/junk pattern; untracked; outside current dream run",
                "rollback": "regenerated automatically by OS or interpreter if needed; no source data removed",
            }
        )
    return actions


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


def is_internal_maintenance_noise_workstream(item: dict[str, Any] | None) -> bool:
    if not item:
        return False
    workstream_id = first_line(item.get("workstream_id"))
    if workstream_id in DREAM_INTERNAL_NOISE_WORKSTREAM_IDS:
        return True
    text = " ".join(
        first_line(value).lower()
        for value in (
            item.get("workstream_id"),
            item.get("title"),
            item.get("handoff_summary"),
            item.get("current_next_action"),
        )
    )
    return "hard probes" in text or "hp2026" in text or "hp26" in text


def risk_is_internal_maintenance_noise(item: dict[str, Any]) -> bool:
    workstream_id = first_line(item.get("workstream_id"))
    if workstream_id in DREAM_INTERNAL_NOISE_WORKSTREAM_IDS:
        return True
    text = " ".join(
        first_line(value).lower()
        for value in (
            item.get("workstream_id"),
            item.get("title"),
            item.get("evidence"),
            item.get("label"),
        )
    )
    return "hp26_photon_id_talk" in text or "hard probes" in text or "hp2026" in text


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


def lane_definition_by_id(lane_id: str) -> dict[str, Any]:
    for item in DREAM_LANE_DEFINITIONS:
        if item["lane_id"] == lane_id:
            return item
    raise KeyError(f"unknown lane_id: {lane_id}")


def automation_state(register_data: dict[str, Any]) -> dict[str, Any]:
    root = codex_home() / "automations"
    expected = list(EXPECTED_DREAM_AUTOMATION_IDS)
    optional = list(OPTIONAL_DREAM_AUTOMATION_IDS)
    legacy = list(LEGACY_DREAM_AUTOMATION_IDS)
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
                    "target_thread_id": parsed.get("target_thread_id", ""),
                    "thread_binding": "fixed_thread" if parsed.get("target_thread_id") else "fresh_chat_per_run",
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
        if str(row.get("run_id") or "").startswith("simulation-"):
            continue
        for finding in row.get("top_findings") or []:
            if not isinstance(finding, dict):
                continue
            fingerprint = f"{finding.get('kind')}|{finding.get('workstream_id') or finding.get('title') or finding.get('evidence')}"
            if "hp26_photon_id_talk" in fingerprint or "hard probes" in fingerprint.lower() or "hp2026" in fingerprint.lower():
                continue
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


def recurrence_protocol_for(item: dict[str, Any]) -> str:
    kind = str(item.get("kind") or "")
    target = str(item.get("workstream_id") or item.get("title") or item.get("evidence") or "")
    protocol = HANDLED_RECURRENCE_PROTOCOLS.get(kind, "")
    if protocol == "register_workstream_refresh_contract" and target == "global":
        return ""
    return protocol


def split_recurring_hotspots(world: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    handled: list[dict[str, Any]] = []
    unhandled: list[dict[str, Any]] = []
    for item in world.get("dream_history", {}).get("hotspots") or []:
        protocol = recurrence_protocol_for(item)
        if protocol:
            copied = dict(item)
            copied["handled_by"] = protocol
            handled.append(copied)
        else:
            unhandled.append(item)
    return handled, unhandled


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


def date_key(value: datetime) -> str:
    return value.astimezone().date().isoformat()


def simulation_group_prefix(run_id: str) -> str | None:
    if not run_id.startswith("simulation-"):
        return None
    match = re.match(r"^(simulation-.+?)-(?:night\d+|day\d+)-", run_id)
    if match:
        return match.group(1)
    parts = run_id.split("-")
    return "-".join(parts[:3]) if len(parts) >= 3 else run_id


def shadow_pilot_dates(history: list[dict[str, Any]], run_id: str) -> list[str]:
    simulation_prefix = simulation_group_prefix(run_id)
    dates: list[str] = []
    for row in history:
        history_run_id = str(row.get("run_id") or "")
        if simulation_prefix:
            if not history_run_id.startswith(f"{simulation_prefix}-"):
                continue
        elif history_run_id.startswith("simulation-"):
            continue
        maintenance = row.get("evolutionary_maintenance") if isinstance(row.get("evolutionary_maintenance"), dict) else {}
        pilot = maintenance.get("shadow_pilot") if isinstance(maintenance.get("shadow_pilot"), dict) else {}
        if pilot.get("name") != "four_night_shadow_pilot":
            continue
        day = first_line(pilot.get("pilot_date"))
        if day and day not in dates:
            dates.append(day)
    return dates[-4:]


def shadow_pilot_state(now: datetime, history: list[dict[str, Any]], run_id: str) -> dict[str, Any]:
    dates = shadow_pilot_dates(history, run_id)
    today = date_key(now)
    ordered = list(dates)
    if today not in ordered:
        ordered.append(today)
    if len(ordered) > 4:
        night_number = len(ordered)
        post_pilot = True
    else:
        night_number = ordered.index(today) + 1
        post_pilot = False
    stage = SHADOW_PILOT_STAGES[min(max(night_number, 1), 4)]
    stage_name = "Post-Pilot - Waiting Human Enablement" if post_pilot else stage["name"]
    status = "shadow_complete_waiting_human_enablement" if post_pilot else stage["status"]
    output_level = "locked_until_human_enablement" if post_pilot else stage["output_level"]
    description = (
        "The four-night shadow pilot is complete; approval-ready packets remain disabled until Justin explicitly enables them."
        if post_pilot
        else stage["description"]
    )
    return {
        "name": "four_night_shadow_pilot",
        "simulation_group": simulation_group_prefix(run_id) or "",
        "pilot_date": today,
        "night_number": night_number,
        "stage_name": stage_name,
        "status": status,
        "output_level": output_level,
        "description": description,
        "completed_pilot_dates": dates,
        "approval_ready_allowed": False,
        "approval_ready_reason": "disabled until Justin explicitly enables approval-ready packets after reviewing the four-night shadow pilot",
        "human_gate_required": True,
        "mutation_boundary": "proposal_only",
        "allowed_actions": [
            "inventory",
            "score",
            "shadow_rehearsal",
            "validator_draft",
            "rollback_draft",
            "morning_review_queue",
        ],
        "forbidden_actions": [
            "repo_tracked_mutation",
            "sdcc_mutation",
            "condor_submission_or_control",
            "google_drive_slides_linear_gmail_mutation",
            "automatic_cleanup",
            "automatic_schema_promotion",
            "approval_ready_without_human_enablement",
        ],
    }


def build_world(register_data: dict[str, Any], artifact_data: dict[str, Any], now: datetime, mode: str, run_id: str) -> dict[str, Any]:
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
    search_hits = search_local_context(dream_terms(workstreams), max_hits=35 if mode in DREAM_LARGE_MODES else 15)
    local_state = local_state_summary(now)
    storage_hits = storage_signal_hits(search_hits)
    physics_surfaces = physics_surface_summary(counts, live)
    literature_surfaces = discover_literature_surfaces(max_items=50 if mode in DREAM_LARGE_MODES else 25)
    history = load_dream_index()
    automation = automation_state(register_data)
    pilot = shadow_pilot_state(now, history, run_id)
    return {
        "mode": mode,
        "generated_at": now.isoformat(timespec="seconds"),
        "shadow_pilot": pilot,
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
        if str(target) in DREAM_INTERNAL_NOISE_WORKSTREAM_IDS or "hp26_photon_id_talk" in str(target):
            continue
        risks.append(
            risk(
                f"Recurring maintenance hotspot has appeared {item['count']} recent dreams",
                "recurring_hotspot",
                77,
                evidence=f"{item.get('kind')}|{target}",
            )
        )
    for item in world["stale"]:
        if is_internal_maintenance_noise_workstream(item):
            continue
        risks.append(risk("Stale live workstream", "stale_state", 95, item))
    for item in world["active_jobs"]:
        if is_internal_maintenance_noise_workstream(item):
            continue
        risks.append(risk("Active job may prompt status check", "active_job_status", 84, item))
    for item in world["live"]:
        if is_internal_maintenance_noise_workstream(item):
            continue
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
    internal_active_p0 = [w for w in world["active_p0"] if not is_internal_maintenance_noise_workstream(w)]
    if len(internal_active_p0) > 3:
        risks.append(
            risk(
                "More than three live P0 workstreams",
                "wip_overload",
                80,
                evidence=", ".join(str(w.get("workstream_id")) for w in internal_active_p0),
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
    if mode in DREAM_LARGE_MODES:
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
    handled_hotspots, unhandled_hotspots = split_recurring_hotspots(world)
    recurring = len(unhandled_hotspots)
    cleanup_candidates = sum(int(item.get("cleanup_candidate_count") or 0) for item in local_state.values())
    wip_overload = 1 if any(item.get("kind") == "wip_overload" and not recurrence_protocol_for(item) for item in risks) else 0
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
            "handled_recurring_hotspots": len(handled_hotspots),
            "cleanup_candidates": cleanup_candidates,
            "active_job_pressure": max(active_jobs - 1, 0),
            "wip_overload": wip_overload,
        },
    }


def schema_promotion_candidates(world: dict[str, Any]) -> list[dict[str, Any]]:
    candidates = []
    _, unhandled_hotspots = split_recurring_hotspots(world)
    for item in unhandled_hotspots:
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
        handled_by = recurrence_protocol_for(item)
        status = "watch"
        if count >= 3:
            status = "ready_for_waking_review"
        if debt.get("status") == "frozen_growth" and count >= 3:
            status = "review_before_new_growth"
        if handled_by:
            status = "handled_by_protocol"
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
                "handled_by": handled_by,
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


def weighted_sum(values: dict[str, float], weights: dict[str, float]) -> float:
    return round(sum(float(values.get(key, 0.0)) * weight for key, weight in weights.items()), 3)


def structure_kind_for(item: dict[str, Any]) -> str:
    kind = str(item.get("kind") or "")
    validator = validator_candidate_for(item)["kind"]
    if kind in {"recurring_hotspot", "automation_drift"}:
        return "schema_card_or_runbook"
    if kind in {"memory_hygiene", "sdcc_clutter_watch"}:
        return "quarantine_or_cleanup_review"
    if kind in {"artifact_gap", "duplicate_guard", "stale_state"}:
        return "validator_or_index"
    if kind == "wip_overload":
        return "cockpit_compression"
    if validator != "manual_review_prompt":
        return "validator_candidate"
    return "simple_note"


def tfe_score_for(item: dict[str, Any]) -> dict[str, Any]:
    base = max(0.0, min(1.0, float(item.get("score") or 0) / 100.0))
    kind = str(item.get("kind") or "")
    safety_bonus = 0.18 if kind in {"duplicate_guard", "automation_drift", "sdcc_clutter_watch"} else 0.08
    reuse_bonus = 0.18 if kind in {"recurring_hotspot", "memory_hygiene", "artifact_gap"} else 0.08
    provenance_bonus = 0.16 if kind in {"artifact_gap", "literature_scout", "physics_hypothesis_scout"} else 0.06
    benefits = {
        "progress": round(base * 0.78, 3),
        "correctness": round(min(1.0, base * 0.74 + provenance_bonus), 3),
        "reuse": round(min(1.0, base * 0.54 + reuse_bonus), 3),
        "safety": round(min(1.0, base * 0.64 + safety_bonus), 3),
        "human_friction_reduction": round(min(1.0, base * 0.58 + 0.12), 3),
        "future_tractability": round(min(1.0, base * 0.66 + reuse_bonus), 3),
    }
    synthetic_risk = 0.18 if evidence_class_for(item).startswith("SYNTHETIC") else 0.08
    costs = {
        "token_cost": round(0.12 + base * 0.16, 3),
        "lookup_cost": round(0.10 + base * 0.12, 3),
        "operational_risk": 0.0,
        "fragmentation": round(0.08 + base * 0.14, 3),
        "maintenance_toil": round(0.10 + base * 0.18, 3),
        "approval_burden": round(0.10 + base * 0.10, 3),
        "synthetic_contamination_risk": synthetic_risk,
    }
    benefit_total = weighted_sum(benefits, TFE_BENEFIT_WEIGHTS)
    cost_total = weighted_sum(costs, TFE_COST_WEIGHTS)
    alternatives = {
        "do_nothing": round(max(0.02, benefit_total * 0.18) / (1 + cost_total * 0.15), 3),
        "simple_note": round(benefit_total * 0.42 / (1 + cost_total * 0.28), 3),
        "single_flat_index": round(benefit_total * 0.58 / (1 + cost_total * 0.42), 3),
        "existing_runbook_update": round(benefit_total * 0.66 / (1 + cost_total * 0.50), 3),
    }
    tfe = round(benefit_total / (1 + cost_total), 3)
    best_alt_name, best_alt_score = max(alternatives.items(), key=lambda pair: pair[1])
    denominator = max(0.1, cost_total - best_alt_score)
    msv = round((benefit_total - best_alt_score) / denominator, 3)
    structure_kind = structure_kind_for(item)
    if tfe <= best_alt_score * 1.08:
        decision = "keep_simple"
    elif msv < 1.0:
        decision = "prefer_simpler_alternative"
    else:
        decision = f"shadow_propose_{structure_kind}"
    return {
        "target": item.get("workstream_id") or item.get("title") or item.get("evidence") or "global",
        "risk_kind": kind,
        "structure_kind": structure_kind,
        "evidence_class": evidence_class_for(item),
        "benefits": benefits,
        "costs": costs,
        "benefit_total": benefit_total,
        "cost_total": cost_total,
        "tfe": tfe,
        "alternatives": alternatives,
        "best_simpler_alternative": {"name": best_alt_name, "score": best_alt_score},
        "marginal_structure_value": msv,
        "decision": decision,
        "mutation_boundary": "proposal_only",
    }


def branch_pressure_for(item: dict[str, Any]) -> dict[str, Any]:
    base = max(0.0, min(1.0, float(item.get("score") or 0) / 100.0))
    kind = str(item.get("kind") or "")
    high_misroute = {
        "duplicate_guard",
        "artifact_gap",
        "automation_drift",
        "recurring_hotspot",
        "sdcc_clutter_watch",
        "physics_hypothesis_scout",
    }
    flow = round(base, 3)
    heterogeneity = 0.75 if kind in {"recurring_hotspot", "memory_hygiene", "sdcc_clutter_watch"} else 0.45
    risk_if_misrouted = 0.9 if kind in high_misroute else 0.55
    reuse_potential = 0.8 if kind in {"recurring_hotspot", "automation_drift", "artifact_gap"} else 0.55
    existing_retrieval_quality = 0.6 if kind in {"artifact_gap", "recurring_hotspot"} else 0.72
    maintenance_tax = 0.55 if structure_kind_for(item) in {"schema_card_or_runbook", "quarantine_or_cleanup_review"} else 0.42
    bpi = round((flow * heterogeneity * risk_if_misrouted * reuse_potential) / max(0.1, existing_retrieval_quality * maintenance_tax), 3)
    if bpi >= 1.0:
        decision = "shadow_split_or_index_candidate"
    elif bpi <= 0.28:
        decision = "compression_or_keep_flat"
    else:
        decision = "keep_simple_with_watch"
    return {
        "node_id": item.get("workstream_id") or kind or "global",
        "risk_kind": kind,
        "objective_flow": flow,
        "heterogeneity": heterogeneity,
        "risk_if_misrouted": risk_if_misrouted,
        "reuse_potential": reuse_potential,
        "existing_retrieval_quality": existing_retrieval_quality,
        "maintenance_tax": maintenance_tax,
        "branch_pressure_index": bpi,
        "decision": decision,
    }


def fuel_class_for(item: dict[str, Any]) -> str:
    kind = str(item.get("kind") or "")
    if kind in {"active_job_status", "active_job_evidence"}:
        return "LIVE"
    if kind in {"duplicate_guard", "recurring_hotspot"}:
        return "DUPLICATE"
    if kind in {"stale_state", "artifact_gap"}:
        return "STALE"
    if kind in {"memory_hygiene", "sdcc_clutter_watch"}:
        return "DORMANT"
    if evidence_class_for(item).startswith("SYNTHETIC"):
        return "SYNTHETIC"
    return "UNOWNED"


def targeted_findings_payload(world: dict[str, Any], risks: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for item in risks[:12]:
        rows.append(
            {
                "label": item.get("label"),
                "kind": item.get("kind"),
                "score": item.get("score"),
                "target": item.get("workstream_id") or item.get("title") or item.get("evidence") or "global",
                "evidence": item.get("evidence"),
                "evidence_class": evidence_class_for(item),
                "fuel_class": fuel_class_for(item),
                "validator_candidate": validator_candidate_for(item),
                "tfe": tfe_score_for(item),
                "branch_pressure": branch_pressure_for(item),
                "pilot_status": world.get("shadow_pilot", {}).get("status"),
            }
        )
    return rows


def branch_pressure_report_payload(risks: list[dict[str, Any]]) -> dict[str, Any]:
    rows = [branch_pressure_for(item) for item in risks[:12]]
    return {
        "synthetic": True,
        "synthetic_header": SYNTHETIC_HEADER,
        "mutation_boundary": "proposal_only",
        "thresholds": {"split": 1.0, "compress": 0.28},
        "rows": rows,
        "summary": {
            "split_candidates": sum(1 for item in rows if item["decision"] == "shadow_split_or_index_candidate"),
            "compression_candidates": sum(1 for item in rows if item["decision"] == "compression_or_keep_flat"),
            "keep_simple": sum(1 for item in rows if item["decision"] == "keep_simple_with_watch"),
        },
    }


def memory_homeostasis_payload(world: dict[str, Any], risks: list[dict[str, Any]]) -> dict[str, Any]:
    local_state = world.get("local_state") or {}
    rows: list[dict[str, Any]] = []
    for label, summary in sorted(local_state.items()):
        cleanup_count = int(summary.get("cleanup_candidate_count") or 0)
        rows.append(
            {
                "node_id": f"local_state.{label}",
                "fuel_class": "DORMANT" if cleanup_count else "LIVE",
                "recommended_action": "propose_cleanup_review" if cleanup_count else "keep_warm",
                "directories": summary.get("directories"),
                "bytes": summary.get("bytes"),
                "cleanup_candidates": cleanup_count,
                "retrieval_policy": "visible_to_doctor_not_default_context" if cleanup_count else "normal",
            }
        )
    for item in risks[:8]:
        rows.append(
            {
                "node_id": item.get("workstream_id") or item.get("kind"),
                "fuel_class": fuel_class_for(item),
                "recommended_action": "quarantine_or_review" if fuel_class_for(item) in {"STALE", "DUPLICATE", "SYNTHETIC"} else "keep_hot",
                "evidence_class": evidence_class_for(item),
                "retrieval_policy": "proposal_only_until_waking_validation",
            }
        )
    return {
        "synthetic": True,
        "synthetic_header": SYNTHETIC_HEADER,
        "mutation_boundary": "proposal_only",
        "fuel_classes": FUEL_CLASSES,
        "rows": rows,
        "rules": {
            "quarantine": "not retrieved by default, visible in audit, restorable only after waking validation",
            "compression": "summary plus source pointers only; never destructive deletion",
            "promotion": "requires real evidence, validator pass, rollback, and Justin approval before live mutation",
        },
    }


def internal_target_for(item: dict[str, Any]) -> str:
    return (
        first_line(item.get("target"))
        or first_line(item.get("workstream_id"))
        or first_line(item.get("title"))
        or first_line(item.get("evidence"))
        or first_line(item.get("kind"))
        or "global"
    )


def internal_evolution_tier_for(item: dict[str, Any]) -> str:
    kind = str(item.get("kind") or "")
    if kind in {"memory_hygiene"}:
        return "local_ignored_auto"
    if kind in {"sdcc_clutter_watch"}:
        return "read_only_external_probe_auto"
    if kind in {"active_job_status", "active_job_evidence", "wip_overload", "recurring_hotspot", "artifact_gap", "stale_state"}:
        return "tracked_internal_with_rollback"
    return "proposal_only_review"


def internal_evolution_files_for(item: dict[str, Any]) -> list[str]:
    kind = str(item.get("kind") or "")
    if kind == "wip_overload":
        return [
            "scripts/os/heartbeat/codex_os_nightly_heartbeat.py",
            "agent_context/policies/CODEX_OPERATING_SYSTEM.md",
        ]
    if kind in {"active_job_status", "active_job_evidence", "recurring_hotspot", "stale_state"}:
        return [
            "scripts/os/dream/codex_os_dream.py",
            "scripts/os/safety/codex_os_doctor.py",
            "agent_context/policies/AGENTIC_OS_HARDENING.md",
            "agent_context/CODEX_WORK_REGISTER.yaml",
        ]
    if kind == "artifact_gap":
        return [
            "scripts/os/artifacts/codex_artifact_registry.py",
            "agent_context/ARTIFACT_REGISTRY.yaml",
        ]
    if kind == "sdcc_clutter_watch":
        return [
            "scripts/sdcc/runtime/audit/checkout_hygiene_probe.py",
            "scripts/sdcc/SDCC_RUNTIME_INDEX.yaml",
            "agent_context/policies/SDCC_OPERATIONS.md",
            "agent_context/policies/STORAGE_AND_CLEANUP.md",
        ]
    if kind == "memory_hygiene":
        return [
            "agent_context/local/dreams/",
            "agent_context/local/chatgpt_research/",
        ]
    return ["agent_context/policies/AGENTIC_OS_DREAMING.md"]


def internal_evolution_action_for(item: dict[str, Any]) -> str:
    kind = str(item.get("kind") or "")
    target = internal_target_for(item)
    if kind == "wip_overload":
        return "compress the morning heartbeat so science/talk priorities stay first and OS hygiene moves to appendix detail"
    if kind in {"active_job_status", "active_job_evidence"}:
        return f"turn repeated status pressure for `{target}` into one read-only evidence-refresh validator or runbook"
    if kind == "stale_state":
        return f"turn stale-state pressure for `{target}` into a refresh/archive runbook with exact evidence, last_verified, next_check, and stale_after rules"
    if kind == "recurring_hotspot":
        return "merge repeated dream prose into one validator-backed adaptation card with decay rules"
    if kind == "artifact_gap":
        return f"add or repair artifact-registry coverage for `{target}` before future claims use it as evidence"
    if kind == "sdcc_clutter_watch":
        return "run the neutral checkout hygiene probe and write only local recommendations unless a separate waking cleanup is approved"
    if kind == "memory_hygiene":
        return "rank old local dream/research packs for review, keeping source pointers and deleting nothing automatically"
    return f"keep `{target}` as a low-priority proposal until it repeats or blocks traversal"


def internal_evolution_queue_payload(world: dict[str, Any], maintenance: dict[str, Any], risks: list[dict[str, Any]]) -> dict[str, Any]:
    targeted = maintenance.get("targeted_findings") or []
    rows: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    for item in targeted:
        if not isinstance(item, dict):
            continue
        key = (str(item.get("kind") or ""), internal_target_for(item))
        if key in seen:
            continue
        seen.add(key)
        tier = internal_evolution_tier_for(item)
        tfe = item.get("tfe") if isinstance(item.get("tfe"), dict) else {}
        rows.append(
            {
                "rank": len(rows) + 1,
                "kind": item.get("kind"),
                "target": internal_target_for(item),
                "auto_tier": tier,
                "tfe": tfe.get("tfe"),
                "msv": tfe.get("marginal_structure_value"),
                "decision": tfe.get("decision"),
                "action": internal_evolution_action_for(item),
                "candidate_files": internal_evolution_files_for(item),
                "allowed_without_new_approval": {
                    "local_ignored_auto": "write/update ignored local queues, summaries, and cleanup candidate lists",
                    "read_only_external_probe_auto": "run bounded read-only probes and write local-only recommendations",
                    "tracked_internal_with_rollback": "waking Codex may apply a small tracked OS/policy/runbook patch only after validation commands are known and no external/science state changes",
                    "proposal_only_review": "record only as a proposal",
                }.get(tier, "record only as a proposal"),
                "blocked_without_explicit_approval": [
                    "SDCC mutation",
                    "Condor submission or job control",
                    "scientific output/model deletion or movement",
                    "Google Drive/Slides/Gmail/Linear mutation",
                    "promotion of synthetic findings into real evidence",
                ],
                "validation_commands": [
                    "python3 scripts/codex_os_dream.py validate --latest",
                    "python3 scripts/codex_os_dream.py lane --lane-id status_provenance",
                    "python3 scripts/codex_os_doctor.py --profile daily",
                ],
                "rollback": "revert the exact tracked patch or discard this local dream directory; do not delete source evidence",
            }
        )
        if len(rows) >= 8:
            break
    refresh_targets = [
        internal_target_for(item)
        for item in risks
        if item.get("kind") in {"stale_state", "active_job_status", "active_job_evidence"}
        and not risk_is_internal_maintenance_noise(item)
    ]
    if refresh_targets and not any(row["kind"] == "workstream_refresh_protocol" for row in rows):
        rows.insert(
            0,
            {
                "rank": 1,
                "kind": "workstream_refresh_protocol",
                "target": "register_workstream_refresh_contract",
                "auto_tier": "tracked_internal_with_rollback",
                "tfe": 1.9,
                "msv": 14.0,
                "decision": "promote_generic_protocol_over_workstream_specific_status_advice",
                "action": (
                    "create one reusable workstream refresh protocol: read register evidence and active-job fields, "
                    "classify current/waiting/stale/archive-review, require exact evidence, then update last_verified, "
                    "next_check, stale_after, and active-job notes only through waking validation"
                ),
                "candidate_files": [
                    "scripts/os/register/codex_work_register_stale.py",
                    "scripts/os/safety/codex_os_doctor.py",
                    "agent_context/policies/MEMORY_AND_STATUS.md",
                    "agent_context/policies/CODEX_OPERATING_SYSTEM.md",
                ],
                "example_targets": refresh_targets[:6],
                "allowed_without_new_approval": "package the protocol as a tracked internal OS patch with rollback and validation; do not change science status automatically",
                "blocked_without_explicit_approval": [
                    "closing workstreams",
                    "marking science outputs valid",
                    "editing SDCC or moving outputs",
                    "promoting synthetic status findings into evidence",
                ],
                "validation_commands": [
                    "python3 scripts/os/register/codex_work_register_stale.py agent_context/CODEX_WORK_REGISTER.yaml",
                    "python3 scripts/codex_os_doctor.py --profile daily",
                    "python3 scripts/codex_os_dream.py validate --latest",
                ],
                "rollback": "revert the exact tracked protocol patch; no workstream status or science output is changed by the dream",
            },
        )
    stale_targets = [
        internal_target_for(item)
        for item in risks
        if item.get("kind") == "stale_state" and not risk_is_internal_maintenance_noise(item)
    ]
    if stale_targets and not any(row["kind"] == "stale_workstream_refresh_protocol" for row in rows):
        rows.insert(
            0,
            {
                "rank": 1,
                "kind": "stale_workstream_refresh_protocol",
                "target": "register_stale_state_contract",
                "auto_tier": "tracked_internal_with_rollback",
                "tfe": 1.75,
                "msv": 12.0,
                "decision": "promote_generic_protocol_over_one_off_runbook",
                "action": (
                    "create one reusable refresh/archive protocol for stale workstreams: read register evidence, "
                    "classify current/waiting/stale/archive-review, require exact evidence, and update last_verified, "
                    "next_check, and stale_after only through waking validation"
                ),
                "candidate_files": [
                    "scripts/os/register/codex_work_register_stale.py",
                    "scripts/os/safety/codex_os_doctor.py",
                    "agent_context/policies/MEMORY_AND_STATUS.md",
                    "agent_context/policies/AGENTIC_OS_DREAMING.md",
                ],
                "example_targets": stale_targets[:6],
                "allowed_without_new_approval": "package the protocol as a tracked internal OS patch with rollback and validation; do not change science status automatically",
                "blocked_without_explicit_approval": [
                    "closing workstreams",
                    "marking science outputs valid",
                    "editing SDCC or moving outputs",
                    "promoting synthetic stale-state findings into evidence",
                ],
                "validation_commands": [
                    "python3 scripts/os/register/codex_work_register_stale.py agent_context/CODEX_WORK_REGISTER.yaml",
                    "python3 scripts/codex_os_doctor.py --profile daily",
                    "python3 scripts/codex_os_dream.py validate --latest",
                ],
                "rollback": "revert the exact tracked protocol patch; no workstream status or science output is changed by the dream",
            },
        )
    for index, row in enumerate(rows, start=1):
        row["rank"] = index
    if not any(row["kind"] == "base_repo_sdcc_hygiene" for row in rows):
        snapshot = repo_root_hygiene_snapshot()
        rows.append(
            {
                "rank": len(rows) + 1,
                "kind": "base_repo_sdcc_hygiene",
                "target": "local_repo_and_sdcc_checkout",
                "auto_tier": "read_only_external_probe_auto",
                "tfe": 1.25 if snapshot["loose_generated_or_tmp_count"] else 1.1,
                "msv": 8.0 if snapshot["loose_generated_or_tmp_count"] else 5.0,
                "decision": "standing_hygiene_probe_and_local_queue",
                "action": "keep base repo and SDCC checkout hygiene in the lane heartbeats: local generated trash is safe to flag, SDCC remains read-only unless separately approved",
                "candidate_files": [
                    "scripts/sdcc/runtime/audit/checkout_hygiene_probe.py",
                    "scripts/os/dream/codex_os_dream.py",
                    "agent_context/policies/SDCC_OPERATIONS.md",
                    "agent_context/policies/STORAGE_AND_CLEANUP.md",
                ],
                "allowed_without_new_approval": "run local scans and bounded read-only SDCC hygiene probes; write findings only to local ignored dream outputs",
                "blocked_without_explicit_approval": [
                    "SDCC mutation",
                    "Condor submission or job control",
                    "scientific output/model deletion or movement",
                    "Google Drive/Slides/Gmail/Linear mutation",
                    "promotion of synthetic findings into real evidence",
                ],
                "validation_commands": [
                    "python3 scripts/codex_os_dream.py validate --latest",
                    "python3 scripts/codex_os_dream.py lane --lane-id path_contract",
                    "python3 scripts/codex_os_doctor.py --profile daily",
                    "python3 scripts/sdcc/runtime/audit/checkout_hygiene_probe.py . --profile local --report-only",
                ],
                "rollback": "discard local dream recommendation files; do not mutate SDCC or delete scientific artifacts",
            }
        )
    return {
        "synthetic": True,
        "synthetic_header": SYNTHETIC_HEADER,
        "mutation_boundary": "proposal_only",
        "purpose": "compress repeated dream findings into a waking internal-evolution queue for Codex brain maintenance",
        "autonomy_boundary": {
            "auto_safe": [
                "ignored local dream/research summaries",
                "read-only SDCC/base hygiene probes",
                "small tracked OS/policy/runbook changes when rollback and validation are explicit",
            ],
            "approval_required": [
                "SDCC edits or cleanup",
                "Condor/job control",
                "science output/model movement or deletion",
                "Google Slides/Drive/Gmail/Linear mutation",
                "physics status or task closure changes",
            ],
        },
        "summary": {
            "row_count": len(rows),
            "local_ignored_auto": sum(1 for row in rows if row["auto_tier"] == "local_ignored_auto"),
            "read_only_external_probe_auto": sum(1 for row in rows if row["auto_tier"] == "read_only_external_probe_auto"),
            "tracked_internal_with_rollback": sum(1 for row in rows if row["auto_tier"] == "tracked_internal_with_rollback"),
            "proposal_only_review": sum(1 for row in rows if row["auto_tier"] == "proposal_only_review"),
            "risk_count_seen": len(risks),
            "maintenance_debt_status": maintenance_debt(world, risks).get("status"),
        },
        "rows": rows,
    }


def repo_root_hygiene_snapshot() -> dict[str, Any]:
    root = Path(".")
    loose_patterns = []
    suspicious = []
    try:
        entries = sorted(root.iterdir(), key=lambda item: item.name)
    except OSError:
        entries = []
    for path in entries:
        name = path.name
        if name in {".git", ".gitignore"}:
            continue
        if name.strip() != name or name.strip() == "":
            suspicious.append(name)
        if name == ".DS_Store" or name == "__pycache__" or name.startswith("tmp_") or name.endswith((".pyc", ".pyo")):
            loose_patterns.append(name)
    return {
        "top_level_count": len(entries),
        "suspicious_name_count": len(suspicious),
        "loose_generated_or_tmp_count": len(loose_patterns),
        "suspicious_names": suspicious[:20],
        "loose_generated_or_tmp": loose_patterns[:20],
    }


def sdcc_base_repo_hygiene_payload(world: dict[str, Any], risks: list[dict[str, Any]]) -> dict[str, Any]:
    storage_hits = world.get("storage_signal_hits") or []
    sdcc_risks = [item for item in risks if item.get("kind") == "sdcc_clutter_watch"]
    return {
        "synthetic": True,
        "synthetic_header": SYNTHETIC_HEADER,
        "mutation_boundary": "proposal_only",
        "purpose": "keep local repo and SDCC checkout architecture clean without creating agent-facing remote artifacts",
        "local_repo_snapshot": repo_root_hygiene_snapshot(),
        "sdcc_contract": {
            "read_only_probe": "scripts/sdcc/runtime/audit/checkout_hygiene_probe.py --report-only (run on SDCC if installed, or stream from local without leaving a remote file)",
            "canonical_base_domains": [
                "runs/",
                "state/",
                "inputs/",
                "models/",
                "evidence/",
                "scratch/",
            ],
            "banned_visible_remote_names": [
                "agent_context",
                ".codex",
                "codex*",
                "THE-*",
                "whitespace-only or leading/trailing-space paths",
            ],
            "future_output_rule": "new SDCC checkout outputs should use runs/state/inputs/models/evidence/scratch and new TG bulk products should use typed thesisAna/recoiljets/<family>/<campaign>/{signal,background}",
        },
        "recorded_storage_signals": storage_hits[:12],
        "sdcc_risk_count": len(sdcc_risks),
        "recommendations": [
            {
                "scope": "local repo",
                "action": "keep root and top-level scripts/macros as command surfaces; place new Codex OS helpers under scripts/os/* and new evidence under typed dataOutput/agent_context registries",
                "auto_tier": "tracked_internal_with_rollback",
            },
            {
                "scope": "SDCC checkout",
                "action": "run read-only checkout hygiene probe automatically in the nightly package when SDCC access is available; write findings only to local ignored dream output",
                "auto_tier": "read_only_external_probe_auto",
            },
            {
                "scope": "SDCC checkout",
                "action": "never create remote agent_context/.codex/codex/THE-* evidence; keep migration manifests and rollback maps local-only",
                "auto_tier": "policy_enforced",
            },
            {
                "scope": "TG bulk",
                "action": "keep historical flat roots read-compatible; do not move them until a separate live-job/reference/rollback audit exists",
                "auto_tier": "approval_required",
            },
        ],
    }


def context_resonance_payload(world: dict[str, Any], risks: list[dict[str, Any]]) -> dict[str, Any]:
    task = (
        "nightly ThesisAnalysis internal architecture maintenance: improve retrieval quality, "
        "surface latent context nudges, suppress stale or synthetic context, and identify negative memories"
    )
    payload = build_context_resonance_payload(task, max_active=3, max_nudges=3, mode="dream")
    risk_rows = []
    for item in risks[:6]:
        risk_rows.append(
            {
                "kind": item.get("kind"),
                "target": item.get("workstream_id") or item.get("title") or item.get("evidence") or "global",
                "relation_type": "same_failure_mode" if item.get("kind") in {"recurring_hotspot", "stale_state"} else "warning_only",
                "evidence_class": evidence_class_for(item),
                "retrieval_policy": "latent_nudge",
                "why_it_surfaced": item.get("label"),
                "required_waking_check": validator_candidate_for(item).get("safe_command"),
            }
        )
    payload["nightly_review"] = {
        "risk_linked_context_count": len(risk_rows),
        "risk_linked_context": risk_rows,
        "memory_engineering_rule": "improve retrieval quality, not raw context volume",
        "future_codex_question": "what old context should help tomorrow without hijacking the current task?",
    }
    payload["context_interference_metrics"] = {
        "false_context_nudges_per_week": "track locally after waking use",
        "stale_memory_retrievals": "track locally after waking use",
        "user_corrections_due_to_old_context": "track locally after waking use",
        "missed_relevant_memory_rate": "track locally after waking use",
    }
    return payload


def render_context_resonance_list(title: str, rows: list[dict[str, Any]]) -> str:
    lines = [f"# {title}", "", SYNTHETIC_HEADER, ""]
    if not rows:
        lines.append("- none")
        return "\n".join(lines) + "\n"
    for item in rows:
        label = item.get("id") or item.get("memory_id") or item.get("trap") or item.get("source") or item.get("target")
        relation = item.get("relation_type") or "n/a"
        policy = item.get("retrieval_policy") or "proposal_only"
        lines.append(f"## `{label}`")
        lines.append(f"- relation_type: `{relation}`")
        lines.append(f"- retrieval_policy: `{policy}`")
        lines.append(f"- evidence_class: `{item.get('evidence_class')}`")
        if item.get("why_it_surfaced"):
            lines.append(f"- why_it_surfaced: {item.get('why_it_surfaced')}")
        if item.get("trap"):
            lines.append(f"- trap: {item.get('trap')}")
        if item.get("first_safe_action"):
            lines.append(f"- first_safe_action: {item.get('first_safe_action')}")
        if item.get("required_waking_check"):
            lines.append(f"- required_waking_check: {item.get('required_waking_check')}")
        lines.append("")
    return "\n".join(lines)


def lane_id_for_risk(item: dict[str, Any]) -> str:
    kind = str(item.get("kind") or "")
    if kind in {"stale_state", "active_job_status", "active_job_evidence", "artifact_gap", "duplicate_guard"}:
        return "status_provenance"
    if kind in {"recurring_hotspot", "wip_overload", "automation_drift"}:
        return "architecture_cohesion"
    if kind in {"memory_hygiene", "sdcc_clutter_watch", "safe_overnight_hygiene"}:
        return "cleanup_storage"
    if kind in {"base_repo_sdcc_hygiene"}:
        return "path_contract"
    if kind in {"literature_scout"}:
        return "research_scout"
    if kind in {"physics_hypothesis_scout", "ideal_final_figure_design"}:
        return "science_scout"
    if kind in {"all_corners_review"}:
        return "architecture_cohesion"
    return "architecture_cohesion"


def lane_heartbeat_payload(world: dict[str, Any], risks: list[dict[str, Any]]) -> dict[str, Any]:
    grouped: dict[str, list[dict[str, Any]]] = {item["lane_id"]: [] for item in DREAM_LANE_DEFINITIONS}
    for item in risks:
        if risk_is_internal_maintenance_noise(item):
            continue
        grouped.setdefault(lane_id_for_risk(item), []).append(item)
    lanes: list[dict[str, Any]] = []
    for definition in DREAM_LANE_DEFINITIONS:
        lane_id = definition["lane_id"]
        lane_risks = grouped.get(lane_id) or []
        lanes.append(
            {
                **definition,
                "synthetic": True,
                "mutation_boundary": "proposal_only",
                "risk_count": len(lane_risks),
                "top_risks": [
                    {
                        "kind": item.get("kind"),
                        "target": item.get("workstream_id") or item.get("evidence") or item.get("title") or "global",
                        "score": item.get("score"),
                    }
                    for item in lane_risks[:4]
                ],
                "overnight_status": "active_lane" if lane_risks or lane_id in {"architecture_cohesion", "cleanup_storage", "research_scout", "context_resonance"} else "standby_lane",
                "handoff_to_master": "write lane-local findings, then let the master heartbeat reconcile conflicts and produce one morning digest",
            }
        )
    return {
        "synthetic": True,
        "synthetic_header": SYNTHETIC_HEADER,
        "mutation_boundary": "proposal_only",
        "purpose": "define parallelizable dream sub-heartbeats while keeping one master overnight reconciler",
        "master_heartbeat": {
            "role": "orchestrator_and_conflict_resolver",
            "runs_lanes_in_parallel": True,
            "single_user_surface": "morning_conversation_digest.md plus nightly_heartbeat.md",
            "conflict_rule": "protected paths, evidence firewall, and no external/science mutation override any lane suggestion",
        },
        "lanes": lanes,
        "engineering_rule": "bifurcate by purpose, not by independent authority; lanes may disagree, master reconciles, doctor validates",
    }


def render_lane_heartbeats(payload: dict[str, Any]) -> str:
    lines = ["# Dream Lane Heartbeats", "", SYNTHETIC_HEADER, ""]
    lines.append("This is the multi-agent shape of the nightly dream: separate purpose lanes, one master reconciler.")
    lines.append("")
    master = payload.get("master_heartbeat") if isinstance(payload.get("master_heartbeat"), dict) else {}
    lines.append("## Master Heartbeat")
    lines.append(f"- role: {master.get('role')}")
    lines.append(f"- runs_lanes_in_parallel: `{master.get('runs_lanes_in_parallel') is True}`")
    lines.append(f"- single_user_surface: `{master.get('single_user_surface')}`")
    lines.append(f"- conflict_rule: {master.get('conflict_rule')}")
    lines.append("")
    lines.append("## Lanes")
    for lane in payload.get("lanes") or []:
        lines.append(f"### `{lane.get('lane_id')}`")
        lines.append(f"- status: `{lane.get('overnight_status')}`")
        lines.append(f"- role: {lane.get('agent_role')}")
        lines.append(f"- purpose: {lane.get('purpose')}")
        lines.append(f"- progression: {lane.get('progression')}")
        lines.append(f"- risk_count: {lane.get('risk_count')}")
        if lane.get("suppression_rule"):
            lines.append(f"- suppression_rule: {lane.get('suppression_rule')}")
        lines.append("- primary_outputs:")
        for output in lane.get("primary_outputs") or []:
            lines.append(f"  - `{output}`")
        top_risks = lane.get("top_risks") if isinstance(lane.get("top_risks"), list) else []
        if top_risks:
            lines.append("- top_risks:")
            for risk_row in top_risks:
                lines.append(
                    f"  - `{risk_row.get('kind')}` target=`{risk_row.get('target')}` score={risk_row.get('score')}"
                )
        lines.append(f"- handoff_to_master: {lane.get('handoff_to_master')}")
        lines.append("")
    lines.append(f"## Engineering Rule\n- {payload.get('engineering_rule')}")
    return "\n".join(lines) + "\n"


def render_single_lane_heartbeat(lane: dict[str, Any]) -> str:
    lines = [f"# Dream Lane `{lane.get('lane_id')}`", "", SYNTHETIC_HEADER, ""]
    lines.append(f"- status: `{lane.get('overnight_status')}`")
    lines.append(f"- role: {lane.get('agent_role')}")
    lines.append(f"- purpose: {lane.get('purpose')}")
    lines.append(f"- progression: {lane.get('progression')}")
    lines.append(f"- mutation_boundary: `{lane.get('mutation_boundary')}`")
    lines.append(f"- handoff_to_master: {lane.get('handoff_to_master')}")
    if lane.get("suppression_rule"):
        lines.append(f"- suppression_rule: {lane.get('suppression_rule')}")
    lines.append("")
    lines.append("## Primary Outputs")
    for output in lane.get("primary_outputs") or []:
        lines.append(f"- `{output}`")
    lines.append("")
    lines.append("## Top Lane Risks")
    top_risks = lane.get("top_risks") if isinstance(lane.get("top_risks"), list) else []
    if top_risks:
        for risk_row in top_risks:
            lines.append(f"- `{risk_row.get('kind')}` target=`{risk_row.get('target')}` score={risk_row.get('score')}")
    else:
        lines.append("- none")
    lines.append("")
    lines.append("## Master Contract")
    lines.append("- This lane has no independent authority; the master heartbeat reconciles conflicts and validates boundaries.")
    return "\n".join(lines) + "\n"


def write_lane_outputs(run_dir: Path, payload: dict[str, Any]) -> None:
    for lane in payload.get("lanes") or []:
        lane_id = str(lane.get("lane_id") or "unknown_lane")
        safe_lane_id = re.sub(r"[^a-zA-Z0-9_.-]+", "_", lane_id).strip("_") or "unknown_lane"
        safe_write(run_dir, f"lanes/{safe_lane_id}/heartbeat.json", render_json_file(lane))
        safe_write(run_dir, f"lanes/{safe_lane_id}/heartbeat.md", render_single_lane_heartbeat(lane))


def cleanup_compaction_payload(now: datetime) -> dict[str, Any]:
    surfaces: list[dict[str, Any]] = []
    for label, root in (("dreams", DREAM_ROOT), ("chatgpt_research", CHATGPT_RESEARCH_ROOT)):
        dirs = [path for path in root.iterdir() if path.is_dir()] if root.exists() else []
        dirs.sort(key=lambda item: item.stat().st_mtime if item.exists() else 0, reverse=True)
        hot_keep = [path.as_posix() for path in dirs[:8]]
        compact_candidates = []
        for path in dirs[8:24]:
            compact_candidates.append(
                {
                    "path": path.as_posix(),
                    "age_days": age_days(path, now),
                    "bytes": dir_bytes(path),
                    "source_pointer": f"{path.as_posix()}/",
                    "recommended_action": "summarize_to_index_before_any_archive_or_prune",
                }
            )
        surfaces.append(
            {
                "surface": label,
                "root": root.as_posix(),
                "directory_count": len(dirs),
                "hot_keep": hot_keep,
                "compact_candidates": compact_candidates,
            }
        )
    return {
        "synthetic": True,
        "synthetic_header": SYNTHETIC_HEADER,
        "mutation_boundary": "proposal_only",
        "purpose": "summarize older local dream/research packs into compact pointers before any future archive or prune",
        "delete_anything": False,
        "move_anything": False,
        "surfaces": surfaces,
    }


def render_cleanup_compaction(payload: dict[str, Any]) -> str:
    lines = ["# Cleanup Compaction Index", "", SYNTHETIC_HEADER, ""]
    lines.append("This is a local-only compaction index. It does not delete, move, or archive source directories.")
    lines.append("")
    for surface in payload.get("surfaces") or []:
        lines.append(f"## `{surface.get('surface')}`")
        lines.append(f"- root: `{surface.get('root')}`")
        lines.append(f"- directory_count: {surface.get('directory_count')}")
        lines.append("- hot_keep:")
        for path in surface.get("hot_keep") or []:
            lines.append(f"  - `{path}`")
        candidates = surface.get("compact_candidates") if isinstance(surface.get("compact_candidates"), list) else []
        lines.append("- compact_candidates:")
        if candidates:
            for item in candidates[:12]:
                lines.append(
                    f"  - `{item.get('path')}` age={item.get('age_days')}d bytes={item.get('bytes')} action={item.get('recommended_action')}"
                )
        else:
            lines.append("  - none")
        lines.append("")
    lines.append("## Rule")
    lines.append("- Any future archive/prune must preserve this source pointer chain and remain local-only.")
    return "\n".join(lines) + "\n"


def autonomy_envelope_payload() -> dict[str, Any]:
    return {
        "synthetic": True,
        "synthetic_header": SYNTHETIC_HEADER,
        "mutation_boundary": "local_internal_auto_safe",
        "purpose": "let dream mode improve Codex internal organization automatically while blocking runtime/science/external mutation",
        "tiers": [
            {
                "tier": "auto_safe",
                "description": AUTONOMY_TIERS["auto_safe"],
                "may_apply_overnight": True,
                "examples": [
                    "remove untracked .DS_Store / AppleDouble / __pycache__ / .pyc artifacts",
                    "refresh ignored dream-run indexes and compact local audit summaries",
                    "write local-only SDCC/base hygiene observations",
                ],
            },
            {
                "tier": "auto_validated",
                "description": AUTONOMY_TIERS["auto_validated"],
                "may_apply_overnight": False,
                "examples": [
                    "prepare validator or runbook patches with exact rollback",
                    "package small OS/policy/index improvements for waking Codex",
                    "require doctor and diff checks before tracked application",
                ],
            },
            {
                "tier": "research_only",
                "description": AUTONOMY_TIERS["research_only"],
                "may_apply_overnight": False,
                "examples": [
                    "prepare sanitized prompts for ChatGPT/Claude/Gemini-style critique",
                    "store source leads and rejected ideas locally",
                    "verify locally before any policy or code promotion",
                ],
            },
            {
                "tier": "blocked_for_waking",
                "description": AUTONOMY_TIERS["blocked_for_waking"],
                "may_apply_overnight": False,
                "examples": [
                    "SDCC runtime edits or cleanup",
                    "Condor submission/job control/merge reruns",
                    "science output, model, TG bulk, Slides, Gmail, Linear, or GitHub mutation",
                ],
            },
        ],
        "hard_guards": {
            "external_mutations_allowed": False,
            "science_state_mutations_allowed": False,
            "sdcc_mutations_allowed": False,
            "repo_tracked_self_modification_by_dream_script": False,
            "synthetic_or_external_model_output_as_evidence": False,
        },
    }


def research_scout_payload(world: dict[str, Any], risks: list[dict[str, Any]]) -> dict[str, Any]:
    debt = maintenance_debt(world, risks)
    top_kinds = sorted({str(item.get("kind") or "unknown") for item in risks[:8]})
    prompt = (
        "Design improvements for a local scientific-analysis agent OS that must reduce context cost, "
        "path drift, stale state, and maintenance debt without mutating science/runtime systems overnight. "
        "Prefer agent-memory, SRE, software-maintenance, and research-workflow ideas that are testable, "
        "small enough to validate locally, and safe against synthetic evidence contamination. "
        "Prefer better use of human-accessible tools and AI interfaces over reinventing equivalent tooling. "
        "If ChatGPT history review or online research would improve the agent OS, say what to inspect, why it matters, "
        "and how to validate the resulting change locally. "
        "Do not assume private project facts; mark uncertainty and provide source leads."
    )
    leverage = {
        "principle": "Use the same tools available to the human more effectively instead of reinventing equivalent machinery.",
        "analogy": "A house robot should use the dishwasher and washer/dryer better than the human can coordinate them, not invent new dishwashing or laundry machines.",
        "history_review_status": "queued_for_waking_or_delegated_lane",
        "history_review_boundary": "nightly scripts must not autonomously browse authenticated ChatGPT history; package the review as a waking/delegated task unless the history is already available in a local artifact",
        "first_class_subtasks": [
            {
                "id": "chatgpt_history_review",
                "goal": "Review recent ChatGPT history for reusable prompt patterns, mode choices, and failure modes that improve Codex throughput.",
                "why": "Recent successful and failed human-visible AI interactions are a direct source of better routing, prompt hygiene, and handoff structure.",
                "morning_action": "Inspect recent ChatGPT threads for prompt/mode patterns worth promoting into delegation policy or dream prompts.",
            },
            {
                "id": "online_research_upgrade",
                "goal": "Use source-backed online research to find better agent-memory, delegation, and workflow designs for the ThesisAnalysis OS.",
                "why": "External critique can surface techniques and interface usage patterns that local dreaming alone will miss.",
                "morning_action": "Run one sanitized ChatGPT-first research prompt focused on better AI/tool orchestration, then verify locally before promotion.",
            },
            {
                "id": "human_tool_orchestration",
                "goal": "Identify where Codex should route work to human-accessible tools such as ChatGPT, Browser, Drive, or Computer Use instead of recreating the same capability internally.",
                "why": "The agent improves fastest when it leverages existing interfaces and tools more efficiently than the human can operate them manually.",
                "morning_action": "Promote one reversible tool-routing improvement that increases leverage without weakening provenance or safety boundaries.",
            },
        ],
    }
    return {
        "synthetic": True,
        "synthetic_header": SYNTHETIC_HEADER,
        "tier": "research_only",
        "status": "queued_local_prompt_only",
        "reason_not_submitted_by_script": "non-interactive nightly scripts must not operate authenticated external-model UIs; Codex can submit through Computer Use in a waking/delegated research lane",
        "recommended_modes": [
            {"interface": "ChatGPT", "mode": "thinking_or_heavy", "use_for": "default paid online research/critique lane"},
            {"interface": "ChatGPT", "mode": "pro_or_deep_research", "use_for": "long source-backed research; user paste-back policy"},
            {
                "interface": "Claude_or_Gemini",
                "mode": "optional_second_opinion_only",
                "use_for": "use only if already available and clearly better than ChatGPT for the specific interface/task",
            },
        ],
        "current_pressure": {
            "maintenance_debt_status": debt.get("status"),
            "risk_kinds_seen": top_kinds,
            "cohesion_score": cohesion_score(world, risks),
        },
        "human_tool_leverage": leverage,
        "sanitized_first_prompt": prompt,
        "recommended_history_sources": [
            "recent ChatGPT threads relevant to dream/OS improvement",
            "validated local research packs under agent_context/local/chatgpt_research/",
            "recent dream outputs that exposed repeated friction or missed leverage",
        ],
        "promotion_rule": "external output becomes code/policy only after local verification, contradiction check, and normal validation",
    }


def deferred_actions_payload(maintenance: dict[str, Any], research: dict[str, Any]) -> list[dict[str, Any]]:
    queue = maintenance.get("internal_evolution_queue") if isinstance(maintenance.get("internal_evolution_queue"), dict) else {}
    rows = queue.get("rows") if isinstance(queue.get("rows"), list) else []
    deferred: list[dict[str, Any]] = []
    for row in rows:
        tier = str(row.get("auto_tier") or "")
        if tier in {"tracked_internal_with_rollback", "proposal_only_review"}:
            deferred.append(
                {
                    "tier": "auto_validated" if tier == "tracked_internal_with_rollback" else "blocked_for_waking",
                    "reason": "tracked repo changes are packaged for waking Codex validation instead of silently self-applied by the dream script",
                    "target": row.get("target"),
                    "action": row.get("action"),
                    "candidate_files": row.get("candidate_files") or [],
                    "validation_commands": row.get("validation_commands") or [],
                }
            )
    deferred.append(
        {
            "tier": "research_only",
            "reason": research.get("reason_not_submitted_by_script"),
            "target": "external_model_research_lane",
            "action": "submit sanitized prompt through Codex/Computer Use when a waking or dedicated research lane is active",
            "candidate_files": ["agent_context/local/chatgpt_research/", "agent_context/policies/ASK_CHATGPT_DELEGATION.md"],
            "validation_commands": ["python3 scripts/os/research/codex_chatgpt_research_pack.py validate <pack_dir>"],
        }
    )
    for target in (
        "SDCC runtime edits or cleanup",
        "Condor/job control/merge reruns",
        "science output/model/TG bulk movement or deletion",
        "Slides/Gmail/Linear/GitHub mutation",
        "physics-status or task-closure changes",
    ):
        deferred.append(
            {
                "tier": "blocked_for_waking",
                "reason": "external, runtime, or science-facing state requires explicit waking scope and guard checks",
                "target": target,
                "action": "do not perform in dream mode",
                "candidate_files": [],
                "validation_commands": [],
            }
        )
    return deferred


def autonomous_maintenance_payload(run_dir: Path, world: dict[str, Any], risks: list[dict[str, Any]], maintenance: dict[str, Any]) -> dict[str, Any]:
    envelope = autonomy_envelope_payload()
    candidates = collect_auto_safe_generated_trash(run_dir)
    changed = remove_auto_safe_generated_trash(candidates)
    if not changed:
        changed = [
            {
                "tier": "auto_safe",
                "performed": False,
                "action": "no_generated_junk_cleanup_needed",
                "path": ".",
                "reason": "no untracked generated trash matched the auto-safe pattern set outside the current dream run",
                "rollback": "no state changed",
            }
        ]
    research = research_scout_payload(world, risks)
    deferred = deferred_actions_payload(maintenance, research)
    return {
        "synthetic": True,
        "synthetic_header": SYNTHETIC_HEADER,
        "mutation_boundary": "local_internal_auto_safe",
        "autonomy_envelope": envelope,
        "changed_actions": {
            "synthetic": True,
            "synthetic_header": SYNTHETIC_HEADER,
            "mutation_boundary": "local_internal_auto_safe",
            "repo_tracked_mutations_performed": False,
            "external_mutations_performed": False,
            "science_mutations_performed": False,
            "rows": changed,
        },
        "deferred_actions": {
            "synthetic": True,
            "synthetic_header": SYNTHETIC_HEADER,
            "mutation_boundary": "proposal_only",
            "rows": deferred,
        },
        "research_scout": research,
        "validation_summary": {
            "synthetic": True,
            "synthetic_header": SYNTHETIC_HEADER,
            "mutation_boundary": "local_internal_auto_safe",
            "commands": AUTONOMY_VALIDATION_COMMANDS,
            "status": "commands_listed_for_heartbeat_and_waking_validation",
        },
        "morning_digest": {
            "synthetic": True,
            "synthetic_header": SYNTHETIC_HEADER,
            "mutation_boundary": "local_internal_auto_safe",
            "changed_count": sum(1 for item in changed if item.get("performed") is True),
            "deferred_count": len(deferred),
            "research_status": research.get("status"),
            "summary": "auto-safe internal cleanup ran; tracked/external/science changes stayed deferred",
        },
    }


def controlled_burn_payload(world: dict[str, Any], risks: list[dict[str, Any]]) -> dict[str, Any]:
    pilot = world.get("shadow_pilot") or {}
    night = int(pilot.get("night_number") or 1)
    candidate = next((item for item in risks if item.get("kind") in {"recurring_hotspot", "memory_hygiene", "sdcc_clutter_watch", "artifact_gap"}), None)
    proposals: list[dict[str, Any]] = []
    if candidate and night >= 3:
        raw_target = candidate.get("workstream_id") or candidate.get("title") or candidate.get("evidence") or "global"
        evidence_summary = first_line(candidate.get("evidence")) or first_line(candidate.get("label"))
        burn_unit = f"dream_recurrence::{evidence_summary}" if candidate.get("kind") == "recurring_hotspot" else str(raw_target)
        action = (
            "merge_repeated_hotspot_into_single_validator_or_runbook_proposal"
            if candidate.get("kind") == "recurring_hotspot"
            else "shadow_layout_and_redirect_rehearsal"
        )
        proposals.append(
            {
                "proposal_id": f"PRB-SHADOW-{date_key(parse_when(world['generated_at']) or now_utc()).replace('-', '')}-001",
                "proposal_type": "prescribed_reorganization_burn",
                "status": "draft_shadow_only" if night == 3 else "calibration_only",
                "approval_ready": False,
                "burn_unit": burn_unit,
                "fuel_class": fuel_class_for(candidate),
                "source_evidence": [
                    {
                        "evidence_class": "derived",
                        "summary": evidence_summary,
                    }
                ],
                "synthetic_inputs": [
                    {
                        "allowed_use": "stress-test only",
                        "dream_case": "shadow_pilot_controlled_burn_rehearsal",
                    }
                ],
                "proposed_action": {
                    "action": action,
                    "delete_anything": False,
                    "archive_only": False,
                    "create_redirects": "proposal_only_if_a_future_file_move_is_approved",
                    "requires_human_approval": True,
                },
                "expected_benefit": {
                    "lookup_steps_reduction_estimate": 2,
                    "token_cost_reduction_estimate": 0.18,
                    "risk_reduction_estimate": 0.25,
                },
                "risks": [
                    "false positive if source evidence is stale",
                    "over-structuring if a simple note would solve the issue",
                ],
                "rollback_plan": {
                    "checkpoint_required": True,
                    "restore_map": "shadow_only_no_live_changes",
                },
                "validators": {
                    "provenance_firewall": "shadow_pending",
                    "synthetic_as_evidence": "pass_by_design",
                    "reference_integrity": "shadow_pending",
                    "protected_paths": "pass_by_design",
                },
            }
        )
    return {
        "synthetic": True,
        "synthetic_header": SYNTHETIC_HEADER,
        "mutation_boundary": "proposal_only",
        "status": pilot.get("status"),
        "night_number": night,
        "approval_ready_allowed": False,
        "stop_conditions": [
            "synthetic artifact would become canonical",
            "real source pointer would be lost",
            "protected file would be modified",
            "rollback cannot be generated",
            "proposal requires hidden external-system mutation",
            "validator confidence is below threshold",
            "doctor cannot explain the change in plain language",
        ],
        "proposals": proposals,
    }


def evolutionary_maintenance_payload(world: dict[str, Any], risks: list[dict[str, Any]]) -> dict[str, Any]:
    targeted = targeted_findings_payload(world, risks)
    branch = branch_pressure_report_payload(risks)
    homeostasis = memory_homeostasis_payload(world, risks)
    burn = controlled_burn_payload(world, risks)
    lane_heartbeats = lane_heartbeat_payload(world, risks)
    resonance = context_resonance_payload(world, risks)
    generated_at = parse_when(world.get("generated_at")) or now_utc()
    cleanup_compaction = cleanup_compaction_payload(generated_at)
    maintenance_stub = {
        "targeted_findings": targeted,
        "branch_pressure": branch,
        "memory_homeostasis": homeostasis,
        "controlled_burn": burn,
        "lane_heartbeats": lane_heartbeats,
        "context_resonance": resonance,
        "cleanup_compaction": cleanup_compaction,
    }
    internal_queue = internal_evolution_queue_payload(world, maintenance_stub, risks)
    sdcc_hygiene = sdcc_base_repo_hygiene_payload(world, risks)
    approval_ready_count = sum(1 for item in burn.get("proposals", []) if item.get("approval_ready") is True)
    return {
        "version": 1,
        "name": "evolutionary_maintenance_architecture",
        "shadow_pilot": world.get("shadow_pilot"),
        "objective": "increase Thesis Flow Efficiency without increasing autonomy during the four-night shadow pilot",
        "analogy_boundary": "Murray-style branching, sleep replay, synaptic homeostasis, pruning, controlled burns, SRE, and technical-debt control are engineering analogies, not biological proof.",
        "evidence_classes": SHADOW_EVIDENCE_CLASSES,
        "optimization_principles": {
            "tfe": "benefit / (1 + cost), compared against do_nothing/simple_note/single_flat_index/existing_runbook_update",
            "msv": "marginal value of new structure over the simplest safe alternative",
            "bpi": "flow * heterogeneity * risk_if_misrouted * reuse_potential over retrieval quality and maintenance tax",
        },
        "targeted_findings": targeted,
        "branch_pressure": branch,
        "memory_homeostasis": homeostasis,
        "controlled_burn": burn,
        "lane_heartbeats": lane_heartbeats,
        "context_resonance": resonance,
        "cleanup_compaction": cleanup_compaction,
        "internal_evolution_queue": internal_queue,
        "sdcc_base_repo_hygiene": sdcc_hygiene,
        "summary": {
            "targeted_finding_count": len(targeted),
            "internal_evolution_queue_count": len(internal_queue.get("rows") or []),
            "lane_count": len(lane_heartbeats.get("lanes") or []),
            "context_resonance_nudge_count": len(resonance.get("latent_context_nudges") or []),
            "negative_memory_candidate_count": len(resonance.get("negative_memories") or []),
            "approval_ready_count": approval_ready_count,
            "shadow_only": True,
            "human_gate_required": True,
        },
    }


def render_json_file(payload: Any) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def yaml_value(value: Any, indent: int = 0) -> list[str]:
    space = " " * indent
    if isinstance(value, dict):
        lines: list[str] = []
        for key, child in value.items():
            if isinstance(child, (dict, list)):
                lines.append(f"{space}{key}:")
                lines.extend(yaml_value(child, indent + 2))
            else:
                lines.append(f"{space}{key}: {json.dumps(child)}")
        return lines
    if isinstance(value, list):
        lines = []
        if not value:
            return [f"{space}[]"]
        for child in value:
            if isinstance(child, (dict, list)):
                lines.append(f"{space}-")
                lines.extend(yaml_value(child, indent + 2))
            else:
                lines.append(f"{space}- {json.dumps(child)}")
        return lines
    return [f"{space}{json.dumps(value)}"]


def render_yaml_payload(payload: dict[str, Any]) -> str:
    return "\n".join(yaml_value(payload)) + "\n"


def render_path_contract_drift(world: dict[str, Any], maintenance: dict[str, Any]) -> str:
    pilot = maintenance["shadow_pilot"]
    lines = ["# Path Contract Drift", "", SYNTHETIC_HEADER, ""]
    lines.append(f"- pilot_stage: {pilot['stage_name']}")
    lines.append("- mutation_boundary: proposal_only")
    lines.append("- approval_ready_allowed: false")
    lines.append("")
    lines.append("## Drift Signals")
    for item in world.get("storage_signal_hits") or []:
        lines.append(f"- `{item['path']}:{item['line']}` {item['snippet']}")
    if not (world.get("storage_signal_hits") or []):
        lines.append("- No path/storage drift signal found in local recorded context.")
    lines.append("")
    lines.append("## Contract Rule")
    lines.append("- Record canonical-vs-legacy path drift here first; waking Codex must verify before editing scripts, SDCC, or transfer policy.")
    return "\n".join(lines) + "\n"


def render_promotion_candidates(maintenance: dict[str, Any]) -> str:
    pilot = maintenance["shadow_pilot"]
    lines = ["# Promotion Candidates", "", SYNTHETIC_HEADER, ""]
    lines.append(f"- pilot_stage: {pilot['stage_name']}")
    lines.append("- approval_ready_allowed: false")
    lines.append("- reason: four-night shadow pilot must be reviewed by Justin before approval-ready packets exist")
    lines.append("")
    lines.append("## TFE-Ranked Shadow Candidates")
    for item in maintenance["targeted_findings"][:8]:
        tfe = item["tfe"]
        lines.append(
            f"- `{item['kind']}` target={item['target']} tfe={tfe['tfe']} "
            f"msv={tfe['marginal_structure_value']} decision={tfe['decision']}"
        )
    if not maintenance["targeted_findings"]:
        lines.append("- none")
    return "\n".join(lines) + "\n"


def render_morning_approval_queue(maintenance: dict[str, Any]) -> str:
    pilot = maintenance["shadow_pilot"]
    lines = ["# Morning Approval Queue", "", SYNTHETIC_HEADER, ""]
    lines.append("This queue is shadow-only during the first four nights.")
    lines.append("")
    lines.append(f"- pilot_stage: {pilot['stage_name']}")
    lines.append(f"- pilot_status: {pilot['status']}")
    lines.append("- approval_ready_items: 0")
    lines.append("- human_gate_required: true")
    lines.append("")
    lines.append("## Morning Review Questions")
    lines.append("- Which findings were useful versus noisy?")
    lines.append("- Did any candidate overreach beyond real evidence?")
    lines.append("- Did TFE/MSV/BPI prefer simple alternatives where appropriate?")
    lines.append("- Is the controlled-burn rehearsal understandable and reversible?")
    lines.append("- Should approval-ready packets remain disabled?")
    return "\n".join(lines) + "\n"


def render_internal_evolution_queue(payload: dict[str, Any]) -> str:
    lines = ["# Internal Evolution Queue", "", SYNTHETIC_HEADER, ""]
    lines.append(
        "This queue compresses repeated dream findings into concrete internal-maintenance actions for waking Codex."
    )
    lines.append("It may guide Codex brain cleanup, but it is not approval for science, SDCC, Condor, Slides, Linear, or output mutation.")
    lines.append("")
    summary = payload.get("summary") if isinstance(payload.get("summary"), dict) else {}
    lines.append("## Summary")
    for key in (
        "row_count",
        "local_ignored_auto",
        "read_only_external_probe_auto",
        "tracked_internal_with_rollback",
        "proposal_only_review",
        "maintenance_debt_status",
    ):
        lines.append(f"- `{key}`: {summary.get(key)}")
    lines.append("")
    lines.append("## Autonomy Boundary")
    boundary = payload.get("autonomy_boundary") if isinstance(payload.get("autonomy_boundary"), dict) else {}
    for item in boundary.get("auto_safe") or []:
        lines.append(f"- auto-safe: {item}")
    for item in boundary.get("approval_required") or []:
        lines.append(f"- approval-required: {item}")
    lines.append("")
    lines.append("## Ranked Internal Changes")
    rows = payload.get("rows") if isinstance(payload.get("rows"), list) else []
    if not rows:
        lines.append("- No internal evolution rows were generated.")
    for item in rows:
        lines.append(f"### {item.get('rank')}. `{item.get('kind')}` -> {item.get('target')}")
        lines.append(f"- auto_tier: `{item.get('auto_tier')}`")
        lines.append(f"- tfe: {item.get('tfe')} msv: {item.get('msv')} decision: `{item.get('decision')}`")
        lines.append(f"- action: {item.get('action')}")
        lines.append("- candidate_files:")
        for file_path in item.get("candidate_files") or []:
            lines.append(f"  - `{file_path}`")
        lines.append(f"- allowed_without_new_approval: {item.get('allowed_without_new_approval')}")
        lines.append("- blocked_without_explicit_approval:")
        for blocked in item.get("blocked_without_explicit_approval") or []:
            lines.append(f"  - {blocked}")
        lines.append("- validation_commands:")
        for command in item.get("validation_commands") or []:
            lines.append(f"  - `{command}`")
        lines.append(f"- rollback: {item.get('rollback')}")
        lines.append("")
    return "\n".join(lines) + "\n"


def render_sdcc_base_repo_hygiene(payload: dict[str, Any]) -> str:
    lines = ["# SDCC And Base Repo Hygiene Recommendations", "", SYNTHETIC_HEADER, ""]
    lines.append("These recommendations make SDCC/base organization part of the dream without creating remote agent-facing state.")
    lines.append("")
    local = payload.get("local_repo_snapshot") if isinstance(payload.get("local_repo_snapshot"), dict) else {}
    lines.append("## Local Repo Snapshot")
    lines.append(f"- top_level_count: {local.get('top_level_count')}")
    lines.append(f"- suspicious_name_count: {local.get('suspicious_name_count')}")
    lines.append(f"- loose_generated_or_tmp_count: {local.get('loose_generated_or_tmp_count')}")
    for name in local.get("suspicious_names") or []:
        lines.append(f"  suspicious: `{name}`")
    for name in local.get("loose_generated_or_tmp") or []:
        lines.append(f"  loose/generated: `{name}`")
    lines.append("")
    contract = payload.get("sdcc_contract") if isinstance(payload.get("sdcc_contract"), dict) else {}
    lines.append("## SDCC Contract")
    lines.append(f"- read_only_probe: `{contract.get('read_only_probe')}`")
    lines.append("- canonical_base_domains:")
    for domain in contract.get("canonical_base_domains") or []:
        lines.append(f"  - `{domain}`")
    lines.append("- banned_visible_remote_names:")
    for pattern in contract.get("banned_visible_remote_names") or []:
        lines.append(f"  - `{pattern}`")
    lines.append(f"- future_output_rule: {contract.get('future_output_rule')}")
    lines.append("")
    lines.append("## Recorded Storage Signals")
    signals = payload.get("recorded_storage_signals") if isinstance(payload.get("recorded_storage_signals"), list) else []
    if signals:
        for hit in signals:
            lines.append(f"- `{hit.get('path')}:{hit.get('line')}` {hit.get('snippet')}")
    else:
        lines.append("- No recorded SDCC/storage drift signal was found in local context.")
    lines.append("")
    lines.append("## Recommendations")
    for item in payload.get("recommendations") or []:
        lines.append(f"- `{item.get('scope')}` [{item.get('auto_tier')}]: {item.get('action')}")
    return "\n".join(lines) + "\n"


def render_autonomy_envelope(payload: dict[str, Any]) -> str:
    lines = ["# Dream Autonomy Envelope", "", SYNTHETIC_HEADER, ""]
    lines.append("This is the allowed nightly self-maintenance boundary.")
    lines.append("")
    for tier in payload.get("tiers") or []:
        lines.append(f"## `{tier.get('tier')}`")
        lines.append(f"- may_apply_overnight: `{tier.get('may_apply_overnight') is True}`")
        lines.append(f"- description: {tier.get('description')}")
        for example in tier.get("examples") or []:
            lines.append(f"- example: {example}")
        lines.append("")
    lines.append("## Hard Guards")
    guards = payload.get("hard_guards") if isinstance(payload.get("hard_guards"), dict) else {}
    for key, value in guards.items():
        lines.append(f"- `{key}`: `{value}`")
    return "\n".join(lines) + "\n"


def render_changed_actions(payload: dict[str, Any]) -> str:
    lines = ["# Changed Actions Log", "", SYNTHETIC_HEADER, ""]
    lines.append("This records actual auto-safe local/internal maintenance performed by the dream run.")
    lines.append("")
    lines.append(f"- repo_tracked_mutations_performed: `{payload.get('repo_tracked_mutations_performed') is True}`")
    lines.append(f"- external_mutations_performed: `{payload.get('external_mutations_performed') is True}`")
    lines.append(f"- science_mutations_performed: `{payload.get('science_mutations_performed') is True}`")
    lines.append("")
    rows = payload.get("rows") if isinstance(payload.get("rows"), list) else []
    for item in rows:
        lines.append(f"- `{item.get('action')}` performed=`{item.get('performed') is True}` path=`{item.get('path')}`")
        lines.append(f"  reason: {item.get('reason')}")
        lines.append(f"  rollback: {item.get('rollback')}")
    return "\n".join(lines) + "\n"


def render_deferred_for_justin(payload: dict[str, Any]) -> str:
    lines = ["# Deferred For Justin / Waking Codex", "", SYNTHETIC_HEADER, ""]
    lines.append("These were intentionally not performed by the dream script.")
    lines.append("")
    rows = payload.get("rows") if isinstance(payload.get("rows"), list) else []
    for item in rows:
        lines.append(f"- `{item.get('tier')}` target=`{item.get('target')}`")
        lines.append(f"  action: {item.get('action')}")
        lines.append(f"  reason: {item.get('reason')}")
    return "\n".join(lines) + "\n"


def render_research_synthesis(payload: dict[str, Any]) -> str:
    lines = ["# Research Scout Synthesis", "", SYNTHETIC_HEADER, ""]
    lines.append("External models are research/critique helpers only. No answer here is evidence or approval.")
    lines.append("")
    lines.append(f"- status: `{payload.get('status')}`")
    lines.append(f"- reason_not_submitted_by_script: {payload.get('reason_not_submitted_by_script')}")
    lines.append("")
    lines.append("## Recommended Model Uses")
    for item in payload.get("recommended_modes") or []:
        lines.append(f"- `{item.get('interface')}` `{item.get('mode')}`: {item.get('use_for')}")
    pressure = payload.get("current_pressure") if isinstance(payload.get("current_pressure"), dict) else {}
    lines.append("")
    lines.append("## Current Pressure")
    lines.append(f"- maintenance_debt_status: `{pressure.get('maintenance_debt_status')}`")
    lines.append(f"- cohesion_score: {pressure.get('cohesion_score')}")
    lines.append(f"- risk_kinds_seen: {', '.join(pressure.get('risk_kinds_seen') or []) or 'none'}")
    lines.append("")
    leverage = payload.get("human_tool_leverage") if isinstance(payload.get("human_tool_leverage"), dict) else {}
    lines.append("## Human Tool Leverage")
    lines.append(f"- principle: {leverage.get('principle')}")
    lines.append(f"- analogy: {leverage.get('analogy')}")
    lines.append(f"- history_review_status: `{leverage.get('history_review_status')}`")
    lines.append(f"- history_review_boundary: {leverage.get('history_review_boundary')}")
    lines.append("")
    lines.append("## First-Class Research Scout Subtasks")
    for item in leverage.get("first_class_subtasks") or []:
        lines.append(f"- `{item.get('id')}`: {item.get('goal')}")
        lines.append(f"  why: {item.get('why')}")
        lines.append(f"  morning_action: {item.get('morning_action')}")
    lines.append("")
    lines.append("## Recommended History Sources")
    for source in payload.get("recommended_history_sources") or []:
        lines.append(f"- {source}")
    lines.append("")
    lines.append("## Sanitized First Prompt")
    lines.append("")
    lines.append("```text")
    lines.append(str(payload.get("sanitized_first_prompt") or ""))
    lines.append("```")
    lines.append("")
    lines.append(f"- promotion_rule: {payload.get('promotion_rule')}")
    return "\n".join(lines) + "\n"


def render_human_tool_leverage(payload: dict[str, Any]) -> str:
    leverage = payload.get("human_tool_leverage") if isinstance(payload.get("human_tool_leverage"), dict) else {}
    lines = ["# Human Tool Leverage", "", SYNTHETIC_HEADER, ""]
    lines.append("This is the standing nightly subtask for improving how Codex uses the same AI/tools available to the human.")
    lines.append("")
    lines.append(f"- principle: {leverage.get('principle')}")
    lines.append(f"- analogy: {leverage.get('analogy')}")
    lines.append(f"- history_review_status: `{leverage.get('history_review_status')}`")
    lines.append(f"- history_review_boundary: {leverage.get('history_review_boundary')}")
    lines.append("")
    lines.append("## First-Class Subtasks")
    for item in leverage.get("first_class_subtasks") or []:
        lines.append(f"- `{item.get('id')}`: {item.get('goal')}")
        lines.append(f"  why: {item.get('why')}")
        lines.append(f"  next: {item.get('morning_action')}")
    lines.append("")
    lines.append("## Recommended Sources To Inspect")
    for source in payload.get("recommended_history_sources") or []:
        lines.append(f"- {source}")
    return "\n".join(lines) + "\n"


def render_validation_summary(payload: dict[str, Any]) -> str:
    lines = ["# Validation Summary", "", SYNTHETIC_HEADER, ""]
    lines.append("These are the checks that validate the autonomous internal-maintenance envelope.")
    lines.append("")
    lines.append(f"- status: `{payload.get('status')}`")
    lines.append("")
    for command in payload.get("commands") or []:
        lines.append(f"- `{command}`")
    return "\n".join(lines) + "\n"


def render_morning_digest(payload: dict[str, Any]) -> str:
    lines = ["# Morning Digest", "", SYNTHETIC_HEADER, ""]
    lines.append(str(payload.get("summary") or ""))
    lines.append("")
    lines.append(f"- changed_count: {payload.get('changed_count')}")
    lines.append(f"- deferred_count: {payload.get('deferred_count')}")
    lines.append(f"- research_status: `{payload.get('research_status')}`")
    lines.append("")
    lines.append("Anything outside the auto-safe internal envelope was deferred.")
    return "\n".join(lines) + "\n"


def cleanup_review_rows(world: dict[str, Any], maintenance: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for label, summary in sorted((world.get("local_state") or {}).items()):
        for candidate in summary.get("cleanup_candidates") or []:
            rows.append(
                {
                    "surface": label,
                    "path": candidate.get("path"),
                    "action": "ask_before_archive_or_delete",
                    "reason": candidate.get("reason"),
                    "age_days": candidate.get("age_days"),
                    "bytes": candidate.get("bytes"),
                }
            )
    sdcc = maintenance.get("sdcc_base_repo_hygiene") if isinstance(maintenance.get("sdcc_base_repo_hygiene"), dict) else {}
    local = sdcc.get("local_repo_snapshot") if isinstance(sdcc.get("local_repo_snapshot"), dict) else {}
    for name in local.get("loose_generated_or_tmp") or []:
        rows.append(
            {
                "surface": "local_repo_root",
                "path": name,
                "action": "auto_safe_cleanup_candidate",
                "reason": "loose generated/tmp-looking root artifact",
                "age_days": None,
                "bytes": None,
            }
        )
    for name in local.get("suspicious_names") or []:
        rows.append(
            {
                "surface": "local_repo_root",
                "path": name,
                "action": "inspect_before_any_mutation",
                "reason": "suspicious whitespace or malformed top-level path",
                "age_days": None,
                "bytes": None,
            }
        )
    return rows[:12]


def next_architecture_hypotheses(world: dict[str, Any], risks: list[dict[str, Any]], maintenance: dict[str, Any]) -> list[dict[str, str]]:
    queue = maintenance.get("internal_evolution_queue") if isinstance(maintenance.get("internal_evolution_queue"), dict) else {}
    rows = queue.get("rows") if isinstance(queue.get("rows"), list) else []
    hypotheses: list[dict[str, str]] = []
    for row in rows[:4]:
        hypotheses.append(
            {
                "target": str(row.get("target") or "internal_os"),
                "hypothesis": str(row.get("action") or "improve internal organization"),
                "why": "ranked by the dream's TFE/MSV/BPI scoring as likely to reduce repeated lookup or status friction",
                "safe_next_step": "turn this into a small validator/runbook/index patch with rollback and doctor validation",
            }
        )
    if not hypotheses:
        pressure = ", ".join(sorted({str(item.get("kind") or "unknown") for item in risks[:6]})) or "no dominant risk"
        hypotheses.append(
            {
                "target": "dream_search_depth",
                "hypothesis": "broaden the next pass from generated-junk cleanup to stale-note, duplicate-index, and path-contract drift search",
                "why": f"no high-ranked queue row appeared, but pressure still exists: {pressure}",
                "safe_next_step": "produce a candidate list only; do not delete or move files automatically",
            }
        )
    return hypotheses[:4]


def structural_advancement_rows(world: dict[str, Any], risks: list[dict[str, Any]], maintenance: dict[str, Any]) -> list[dict[str, Any]]:
    queue = maintenance.get("internal_evolution_queue") if isinstance(maintenance.get("internal_evolution_queue"), dict) else {}
    queue_rows = queue.get("rows") if isinstance(queue.get("rows"), list) else []
    rows: list[dict[str, Any]] = []
    for item in queue_rows[:6]:
        kind = str(item.get("kind") or "unknown")
        target = str(item.get("target") or "internal_os")
        action = str(item.get("action") or "improve internal structure")
        if kind in {"active_job_status", "active_job_evidence"}:
            problem = "Repeated status/evidence questions are being rediscovered as prose."
            advancement = "Create a read-only evidence-refresh validator/runbook for this workstream."
            structural_gain = "Future Codex can refresh status through one known path instead of scanning scattered ledgers."
        elif kind == "workstream_refresh_protocol":
            problem = "Active/stale workstream pressure is still turning into repeated per-workstream advice."
            advancement = "Create one generic workstream refresh protocol; use `ml_final_model_ablation_map` as a calibration example, not a bespoke dream chore."
            structural_gain = "Future Codex gets one route for current/waiting/stale/archive-review classification and exact evidence updates instead of rebuilding the logic per workstream."
        elif kind == "wip_overload":
            problem = "The morning cockpit carries too many active threads at once."
            advancement = "Compress OS hygiene into appendix detail and keep science priorities first."
            structural_gain = "Morning review becomes easier to scan without losing the lower-level OS evidence."
        elif kind == "recurring_hotspot":
            problem = "The same dream warning is repeating without becoming infrastructure."
            advancement = "Promote the warning into a validator-backed adaptation card with decay rules."
            structural_gain = "Repeated warnings become checkable structure instead of nightly narrative noise."
        elif kind == "stale_state":
            problem = "A workstream is aging past its refresh window, so future Codex cannot quickly tell whether it is current, waiting, or stale."
            advancement = "Create a refresh/archive runbook that requires exact evidence, last_verified, next_check, and stale_after updates."
            structural_gain = "Stale work becomes a bounded refresh operation instead of a repeated rediscovery problem."
        elif kind == "stale_workstream_refresh_protocol":
            problem = "Stale-state pressure is being handled as one-off workstream advice instead of a reusable register maintenance mechanism."
            advancement = "Create one generic stale-workstream refresh/archive protocol and use `ml_final_model_ablation_map` as the first calibration case."
            structural_gain = "Future Codex gets one repeatable path for stale work instead of rediscovering exact evidence, last_verified, next_check, and stale_after rules every night."
        elif kind == "base_repo_sdcc_hygiene":
            problem = "Local/SDCC root clutter can recur if only checked manually."
            advancement = "Keep a read-only hygiene probe in the nightly package and record canonical-path drift locally."
            structural_gain = "Path drift becomes visible without leaving agent-facing state on SDCC."
        else:
            problem = "A repeated friction point exists but has not been made structural."
            advancement = action
            structural_gain = "The dream should convert this into a reusable check, index, or runbook."
        rows.append(
            {
                "target": target,
                "pressure_kind": kind,
                "problem_in_plain_english": problem,
                "structural_advancement": advancement,
                "why_it_improves_the_os": structural_gain,
                "safe_next_step": action,
                "validation": item.get("validation_commands") or AUTONOMY_VALIDATION_COMMANDS,
                "rollback": item.get("rollback") or "revert the exact tracked patch or discard local proposal artifacts",
                "autonomy_status": "package_for_waking_auto_validated" if item.get("auto_tier") == "tracked_internal_with_rollback" else item.get("auto_tier"),
                "scores_are_supporting_only": {
                    "tfe": item.get("tfe"),
                    "msv": item.get("msv"),
                    "decision": item.get("decision"),
                },
            }
        )
    if not rows:
        rows.append(
            {
                "target": "dream_internal_architecture",
                "pressure_kind": "no_ranked_queue_row",
                "problem_in_plain_english": "The dream did not find a high-confidence structural row.",
                "structural_advancement": "Widen the next scan to duplicate notes, stale policy references, and path-contract drift.",
                "why_it_improves_the_os": "The dream keeps searching for infrastructure gains instead of stopping at no generated junk.",
                "safe_next_step": "produce candidates only; do not delete or move files automatically",
                "validation": AUTONOMY_VALIDATION_COMMANDS,
                "rollback": "discard local proposal artifacts",
                "autonomy_status": "research_or_review_only",
                "scores_are_supporting_only": {},
            }
        )
    return rows


def structural_advancement_payload(world: dict[str, Any], risks: list[dict[str, Any]], maintenance: dict[str, Any]) -> dict[str, Any]:
    rows = structural_advancement_rows(world, risks, maintenance)
    return {
        "synthetic": True,
        "synthetic_header": SYNTHETIC_HEADER,
        "mutation_boundary": "proposal_only",
        "purpose": "turn recurring dream pressure into concrete infrastructure advancements; scores are supporting metadata, not the outcome",
        "row_count": len(rows),
        "rows": rows,
    }


def render_structural_advancements(payload: dict[str, Any]) -> str:
    lines = ["# Structural Advancements", "", SYNTHETIC_HEADER, ""]
    lines.append("This is the main dream improvement ledger. Scores are supporting metadata; the product is concrete infrastructure change.")
    lines.append("")
    rows = payload.get("rows") if isinstance(payload.get("rows"), list) else []
    for index, item in enumerate(rows, start=1):
        lines.append(f"## {index}. `{item.get('target')}`")
        lines.append(f"- pressure: `{item.get('pressure_kind')}`")
        lines.append(f"- problem: {item.get('problem_in_plain_english')}")
        lines.append(f"- structural advancement: {item.get('structural_advancement')}")
        lines.append(f"- why it improves the OS: {item.get('why_it_improves_the_os')}")
        lines.append(f"- safe next step: {item.get('safe_next_step')}")
        lines.append(f"- autonomy status: `{item.get('autonomy_status')}`")
        lines.append(f"- rollback: {item.get('rollback')}")
        lines.append("- validation:")
        for command in item.get("validation") or []:
            lines.append(f"  - `{command}`")
        lines.append("")
    return "\n".join(lines) + "\n"


def deeper_search_angles(world: dict[str, Any], risks: list[dict[str, Any]], cleanup_rows: list[dict[str, Any]]) -> list[str]:
    angles: list[str] = []
    if not cleanup_rows:
        angles.append("Scan stale local notes, duplicate policy references, and old dream/research packs instead of stopping at generated-junk cleanup.")
    if any(item.get("kind") in {"active_job_status", "active_job_evidence"} for item in risks):
        angles.append("Convert repeated active-status/evidence warnings into read-only refresh validators so status questions stop recurring as prose.")
    if any(item.get("kind") == "wip_overload" for item in risks):
        angles.append("Compress the morning cockpit so science priorities stay first and OS hygiene becomes appendix detail.")
    local_state = world.get("local_state") or {}
    if sum(int(item.get("directories") or 0) for item in local_state.values()) > 20:
        angles.append("Review local dream/research-pack retention thresholds so ignored internal memory does not become its own clutter.")
    if not angles:
        angles.append("Run a ChatGPT-heavy research scout for a fresh critique of the current OS maintenance envelope.")
    return angles[:4]


def render_morning_conversation_digest(
    run_id: str,
    world: dict[str, Any],
    risks: list[dict[str, Any]],
    maintenance: dict[str, Any],
    autonomous: dict[str, Any],
) -> str:
    changed = autonomous.get("changed_actions") if isinstance(autonomous.get("changed_actions"), dict) else {}
    changed_rows = changed.get("rows") if isinstance(changed.get("rows"), list) else []
    performed = [item for item in changed_rows if item.get("performed") is True]
    deferred = autonomous.get("deferred_actions") if isinstance(autonomous.get("deferred_actions"), dict) else {}
    deferred_rows = deferred.get("rows") if isinstance(deferred.get("rows"), list) else []
    research = autonomous.get("research_scout") if isinstance(autonomous.get("research_scout"), dict) else {}
    debt = maintenance_debt(world, risks)
    internal_queue = maintenance.get("internal_evolution_queue") if isinstance(maintenance.get("internal_evolution_queue"), dict) else {}
    queue_rows = internal_queue.get("rows") if isinstance(internal_queue.get("rows"), list) else []
    sdcc = maintenance.get("sdcc_base_repo_hygiene") if isinstance(maintenance.get("sdcc_base_repo_hygiene"), dict) else {}
    local_snapshot = sdcc.get("local_repo_snapshot") if isinstance(sdcc.get("local_repo_snapshot"), dict) else {}
    cleanup_rows = cleanup_review_rows(world, maintenance)
    structural = maintenance.get("structural_advancements") if isinstance(maintenance.get("structural_advancements"), dict) else {}
    structural_rows = structural.get("rows") if isinstance(structural.get("rows"), list) else []
    lane_payload = maintenance.get("lane_heartbeats") if isinstance(maintenance.get("lane_heartbeats"), dict) else {}
    lane_rows = lane_payload.get("lanes") if isinstance(lane_payload.get("lanes"), list) else []
    hypotheses = next_architecture_hypotheses(world, risks, maintenance)
    search_angles = deeper_search_angles(world, risks, cleanup_rows)

    lines = ["# Morning Conversation Digest", "", SYNTHETIC_HEADER, ""]
    lines.append("This is the short digest Codex should say back in chat after a dream heartbeat.")
    lines.append("")
    lines.append("## Quick Read")
    lines.append(f"- run: `{run_id}`")
    lines.append(f"- dream validation: expected to be checked by heartbeat/doctor")
    lines.append(f"- automatic changes performed: {len(performed)}")
    lines.append(f"- deferred items: {len(deferred_rows)}")
    lines.append(f"- cleanup review candidates: {len(cleanup_rows)}")
    lines.append(f"- structural advancement candidates: {len(structural_rows)}")
    lines.append(f"- dream lanes: {len(lane_rows)}")
    lines.append(f"- architecture hypotheses: {len(hypotheses)}")
    lines.append(f"- maintenance debt signal: {debt.get('score')}/100 status=`{debt.get('status')}` (diagnostic only)")
    lines.append(f"- local root clutter signals: {local_snapshot.get('loose_generated_or_tmp_count', 0)}")
    lines.append("")
    lines.append("## What I Did Automatically")
    if performed:
        for item in performed[:8]:
            lines.append(f"- `{item.get('action')}` on `{item.get('path')}`: {item.get('reason')}")
    else:
        lines.append("- Ran the heartbeat/dream checks and found no generated junk that fit the auto-safe cleanup allowlist.")
    lines.append("- Wrote the dream package under `agent_context/local/dreams/`.")
    lines.append("- Kept SDCC, Condor, science outputs, external apps, and repo-tracked analysis state untouched.")
    lines.append("")
    lines.append("## What I Learned")
    if risks:
        kinds = ", ".join(sorted({str(item.get("kind") or "unknown") for item in risks[:8]}))
        lines.append(f"- The strongest recurring pressure kinds are: {kinds}.")
    lines.append(f"- Cohesion is still low at {cohesion_score(world, risks)}/100. The useful response is not the number; it is converting repeated friction into validators, runbooks, indexes, and cleaner retrieval paths.")
    lines.append("- ChatGPT remains the default paid online research lane; Claude/Gemini are optional second opinions only when clearly worth the extra flow cost.")
    lines.append("- SDCC/base hygiene should stay read-only in dream mode, with findings saved locally.")
    lines.append("- No automatic cleanup does not mean no improvement; it means the next layer is reviewable cleanup, architecture compression, or research.")
    lines.append("")
    lines.append("## What I Tried / Checked")
    if lane_rows:
        active_lanes = [item.get("lane_id") for item in lane_rows if item.get("overnight_status") == "active_lane"]
        lines.append(f"- Split the dream into purpose lanes under one master heartbeat: {', '.join(active_lanes) if active_lanes else 'all lanes standby'}.")
    lines.append("- Checked local repo hygiene and SDCC/base organization signals from local recorded context.")
    lines.append("- Ranked internal evolution candidates using the TFE/MSV/BPI maintenance scoring stack.")
    lines.append("- Queued a sanitized external-model research prompt, but did not submit it from the noninteractive nightly script.")
    lines.append("- Generated adaptation-card and validator/runbook candidates from repeated warnings.")
    lines.append("")
    lines.append("## What I Changed")
    if performed:
        for item in performed[:8]:
            lines.append(f"- Changed `{item.get('path')}` by `{item.get('action')}`; rollback: {item.get('rollback')}")
    else:
        lines.append("- No file/content change was needed in the final auto-safe cleanup pass.")
    lines.append("")
    lines.append("## Cleanup I Found Or Need To Review")
    if cleanup_rows:
        for item in cleanup_rows:
            detail = f"- `{item.get('surface')}` `{item.get('path')}`: {item.get('action')} because {item.get('reason')}"
            if item.get("age_days") is not None:
                detail += f" (age {item.get('age_days')}d)"
            lines.append(detail)
    else:
        lines.append("- No cleanup candidate crossed the current safe/review threshold, so the next pass should look deeper at stale notes, duplicate references, and retained dream/research packs.")
    lines.append("")
    lines.append("## Concrete Structural Advancements")
    if structural_rows:
        for item in structural_rows[:5]:
            lines.append(f"- `{item.get('target')}`: {item.get('structural_advancement')}")
            lines.append(f"  problem fixed: {item.get('problem_in_plain_english')}")
            lines.append(f"  why it helps: {item.get('why_it_improves_the_os')}")
    else:
        lines.append("- No structural advancement row was generated; this is a dream-quality issue to fix.")
    lines.append("")
    lines.append("## Architecture Improvement Hypotheses")
    for item in hypotheses:
        lines.append(f"- `{item['target']}`: {item['hypothesis']}")
        lines.append(f"  why: {item['why']}")
        lines.append(f"  safe next step: {item['safe_next_step']}")
    lines.append("")
    lines.append("## What I Did Not Change But Think Should Be Improved")
    auto_validated = [item for item in deferred_rows if item.get("tier") == "auto_validated"]
    for item in auto_validated[:6]:
        lines.append(f"- `{item.get('target')}`: {item.get('action')}")
        lines.append(f"  why not automatic: {item.get('reason')}")
    if not auto_validated:
        lines.append("- No tracked internal improvements crossed the auto_validated threshold.")
    lines.append("")
    lines.append("## Blocked / Waking-Only Items")
    blocked = [item for item in deferred_rows if item.get("tier") == "blocked_for_waking"]
    if blocked:
        for item in blocked[:8]:
            lines.append(f"- `{item.get('target')}`: {item.get('reason')}")
    else:
        lines.append("- No blocked waking-only items were recorded.")
    lines.append("")
    lines.append("## Best Next Internal Fix")
    if queue_rows:
        first = queue_rows[0]
        lines.append(f"- `{first.get('target')}`: {first.get('action')}")
        lines.append("  why: it is the highest-ranked internal evolution item and should reduce repeated status-pressure noise.")
    else:
        lines.append("- Keep current envelope and wait for the next recurring signal.")
    lines.append("")
    lines.append("## Next Search Angles If The Obvious Pass Finds Nothing")
    for angle in search_angles:
        lines.append(f"- {angle}")
    lines.append("")
    lines.append("## How To Say This In Chat")
    lines.append("- Start with status, validation, and whether anything changed.")
    lines.append("- Then list the top 1-3 deferred improvements with the reason they were deferred.")
    lines.append("- End with the one internal fix Codex recommends doing next.")
    return "\n".join(lines) + "\n"


def render_shadow_approval_packet_readme(maintenance: dict[str, Any]) -> str:
    pilot = maintenance["shadow_pilot"]
    lines = ["# Approval Packets Disabled During Shadow Pilot", "", SYNTHETIC_HEADER, ""]
    lines.append(f"- pilot_stage: {pilot['stage_name']}")
    lines.append("- approval_ready_allowed: false")
    lines.append("- live_mutation_performed: false")
    lines.append("- exact_approval_required_before_future_packets: Justin must explicitly enable approval-ready proposal packets after reviewing Night 4.")
    lines.append("")
    lines.append("No approval packet in this directory is executable. Use `morning_approval_queue.md` for review questions only.")
    return "\n".join(lines) + "\n"


def heartbeat_signal(
    run_id: str,
    mode: str,
    world: dict[str, Any],
    risks: list[dict[str, Any]],
    maintenance: dict[str, Any] | None = None,
) -> dict[str, Any]:
    checks = [
        "python3 scripts/codex_os_doctor.py --profile daily",
        "python3 scripts/codex_work_register_stale.py agent_context/CODEX_WORK_REGISTER.yaml",
        "python3 scripts/os/register/codex_work_register_stale.py agent_context/CODEX_WORK_REGISTER.yaml --protocol",
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
    handled_hotspots, unhandled_hotspots = split_recurring_hotspots(world)
    if maintenance is None:
        maintenance = evolutionary_maintenance_payload(world, risks)
    autonomous = maintenance.get("autonomous_maintenance") if isinstance(maintenance.get("autonomous_maintenance"), dict) else {}
    autonomous_digest = autonomous.get("morning_digest") if isinstance(autonomous.get("morning_digest"), dict) else {}
    return {
        "version": 7,
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
            "recurring_hotspot_count": len(unhandled_hotspots),
            "handled_recurring_hotspot_count": len(handled_hotspots),
            "raw_recurring_hotspot_count": len(world.get("dream_history", {}).get("hotspots") or []),
            "cleanup_candidate_count": sum(
                int(item.get("cleanup_candidate_count") or 0) for item in (world.get("local_state") or {}).values()
            ),
            "schema_promotion_candidate_count": len(schemas),
            "adaptation_card_count": len(cards),
            "approval_ready_count": maintenance["summary"]["approval_ready_count"],
            "autonomous_changed_count": int(autonomous_digest.get("changed_count") or 0),
            "autonomous_deferred_count": int(autonomous_digest.get("deferred_count") or 0),
        },
        "automation_state": world.get("automation_state"),
        "recurring_hotspots": unhandled_hotspots,
        "handled_recurring_hotspots": handled_hotspots,
        "schema_promotion_candidates": schemas,
        "adaptation_cards": cards,
        "evolutionary_maintenance": maintenance,
        "autonomous_maintenance": autonomous,
        "internal_evolution_queue": maintenance.get("internal_evolution_queue") or {},
        "sdcc_base_repo_hygiene": maintenance.get("sdcc_base_repo_hygiene") or {},
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


def render_report(
    run_id: str,
    mode: str,
    world: dict[str, Any],
    risks: list[dict[str, Any]],
    interactions: list[dict[str, Any]],
    maintenance: dict[str, Any] | None = None,
) -> str:
    top_question = interactions[0]["synthetic_user"] if interactions else "No high-risk synthetic question generated."
    automation = world.get("automation_state") or {}
    handled_hotspots, recurring = split_recurring_hotspots(world)
    debt = maintenance_debt(world, risks)
    schema_candidates = schema_promotion_candidates(world)
    adaptation = adaptation_cards(world, risks)
    if maintenance is None:
        maintenance = evolutionary_maintenance_payload(world, risks)
    pilot = maintenance["shadow_pilot"]
    autonomous = maintenance.get("autonomous_maintenance") if isinstance(maintenance.get("autonomous_maintenance"), dict) else {}
    autonomous_digest = autonomous.get("morning_digest") if isinstance(autonomous.get("morning_digest"), dict) else {}
    checks = [
        "python3 scripts/codex_os_doctor.py --profile daily",
        "python3 scripts/codex_work_register_stale.py agent_context/CODEX_WORK_REGISTER.yaml",
        "python3 scripts/os/register/codex_work_register_stale.py agent_context/CODEX_WORK_REGISTER.yaml --protocol",
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
            "## Four-Night Shadow Pilot",
            "",
            f"- stage: {pilot['stage_name']}",
            f"- status: {pilot['status']}",
            f"- output level: {pilot['output_level']}",
            f"- approval-ready allowed: {str(pilot['approval_ready_allowed']).lower()}",
            f"- human gate required: {str(pilot['human_gate_required']).lower()}",
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
            "## Handled Recurring Hotspots",
            "",
            (
                "\n".join(
                    f"- `{item.get('kind')}` target={item.get('workstream_id') or item.get('title') or item.get('evidence') or 'global'} repeats={item['count']} handled_by={item.get('handled_by')}"
                    for item in handled_hotspots
                )
                if handled_hotspots
                else "- none handled by a current protocol"
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
            "## Internal Evolution Queue",
            "",
            (
                "\n".join(
                    f"- `{item.get('auto_tier')}` `{item.get('kind')}` target={item.get('target')}: {item.get('action')}"
                    for item in (maintenance.get("internal_evolution_queue", {}).get("rows") or [])[:5]
                )
                if (maintenance.get("internal_evolution_queue", {}).get("rows") or [])
                else "- no internal evolution rows generated"
            ),
            "",
            "## SDCC/Base Organization Watch",
            "",
            f"- local loose/generated count: {maintenance.get('sdcc_base_repo_hygiene', {}).get('local_repo_snapshot', {}).get('loose_generated_or_tmp_count', 0)}",
            f"- recorded SDCC/storage signals: {len(maintenance.get('sdcc_base_repo_hygiene', {}).get('recorded_storage_signals') or [])}",
            "- remote mutation remains approval-gated; read-only probes and local recommendations are the auto-safe lane.",
            "",
            "## Autonomous Internal Maintenance",
            "",
            f"- auto-safe changes performed: {int(autonomous_digest.get('changed_count') or 0)}",
            f"- deferred/blocked actions recorded: {int(autonomous_digest.get('deferred_count') or 0)}",
            f"- research scout status: `{autonomous_digest.get('research_status') or 'unknown'}`",
            "- tracked self-modification, SDCC mutation, and science/external actions remain blocked inside the nightly script.",
            "",
            "## Readiness Metrics",
            "",
            json.dumps(readiness_metrics(world, risks, interactions), indent=2, sort_keys=True),
            "",
            "## Thesis Flow Efficiency Top Candidates",
            "",
            "\n".join(
                f"- `{item['kind']}` target={item['target']} tfe={item['tfe']['tfe']} "
                f"msv={item['tfe']['marginal_structure_value']} decision={item['tfe']['decision']}"
                for item in maintenance["targeted_findings"][:5]
            )
            or "- none",
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
            "- Four-night shadow pilot output is not approval-ready.",
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
    handled_hotspots, unhandled_hotspots = split_recurring_hotspots(world)
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
    lines.append("## Recurrence Retirement")
    lines.append(f"- handled_recurring_hotspots: {len(handled_hotspots)}")
    lines.append(f"- unhandled_recurring_hotspots: {len(unhandled_hotspots)}")
    for item in handled_hotspots:
        lines.append(
            f"- handled `{item.get('kind')}` target={item.get('workstream_id') or item.get('title') or item.get('evidence') or 'global'} by `{item.get('handled_by')}`"
        )
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
            handled = f" handled_by={item.get('handled_by')}" if item.get("handled_by") else ""
            lines.append(
                f"- `{item['id']}` repeats={item['repeat_count']} status={item['promotion_status']} validator={item['validator_kind']}{handled}"
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
        if item.get("handled_by"):
            lines.append(f"- handled_by: `{item['handled_by']}`")
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
    maintenance = evolutionary_maintenance_payload(world, risks)
    autonomous = autonomous_maintenance_payload(run_dir, world, risks, maintenance)
    structural = structural_advancement_payload(world, risks, maintenance)
    maintenance["structural_advancements"] = structural
    maintenance["autonomous_maintenance"] = autonomous
    maintenance["summary"]["autonomous_changed_count"] = int(autonomous["morning_digest"].get("changed_count") or 0)
    maintenance["summary"]["autonomous_deferred_count"] = int(autonomous["morning_digest"].get("deferred_count") or 0)
    maintenance["summary"]["research_scout_status"] = autonomous["morning_digest"].get("research_status")
    signal = heartbeat_signal(run_id, mode, world, risks, maintenance)
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
        "evolutionary_maintenance": maintenance,
        "autonomous_maintenance": autonomous,
        "readiness_metrics": readiness_metrics(world, risks, interactions),
        "search_hits": world["search_hits"],
        "heartbeat_signal": signal,
    }
    safe_write(run_dir, "dream_trace.json", json.dumps(trace, indent=2, sort_keys=True) + "\n")
    safe_write(run_dir, "heartbeat_signal.json", json.dumps(signal, indent=2, sort_keys=True) + "\n")
    safe_write(run_dir, "dream_report.md", render_report(run_id, mode, world, risks, interactions, maintenance))
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
    safe_write(
        run_dir,
        "targeted_findings.json",
        render_json_file(
            {
                "synthetic": True,
                "synthetic_header": SYNTHETIC_HEADER,
                "mutation_boundary": "proposal_only",
                "rows": maintenance["targeted_findings"],
            }
        ),
    )
    safe_write(
        run_dir,
        "thesis_flow_scores.json",
        render_json_file(
            {
                "synthetic": True,
                "synthetic_header": SYNTHETIC_HEADER,
                "mutation_boundary": "proposal_only",
                "rows": [item["tfe"] for item in maintenance["targeted_findings"]],
            }
        ),
    )
    safe_write(run_dir, "branch_pressure_report.json", render_json_file(maintenance["branch_pressure"]))
    safe_write(run_dir, "memory_homeostasis_report.json", render_json_file(maintenance["memory_homeostasis"]))
    safe_write(run_dir, "context_resonance_review.md", render_context_resonance_markdown(maintenance["context_resonance"]))
    safe_write(run_dir, "context_resonance_review.json", render_json_file(maintenance["context_resonance"]))
    safe_write(
        run_dir,
        "latent_context_nudges.md",
        render_context_resonance_list("Latent Context Nudges", maintenance["context_resonance"].get("latent_context_nudges") or []),
    )
    safe_write(
        run_dir,
        "negative_memory_candidates.md",
        render_context_resonance_list("Negative Memory Candidates", maintenance["context_resonance"].get("negative_memories") or []),
    )
    safe_write(
        run_dir,
        "retrieval_outcome_candidates.md",
        render_context_resonance_list(
            "Retrieval Outcome Candidates",
            maintenance["context_resonance"].get("retrieval_outcome_candidates") or [],
        ),
    )
    safe_write(run_dir, "controlled_burn_proposals.yaml", render_yaml_payload(maintenance["controlled_burn"]))
    safe_write(run_dir, "lane_heartbeats.json", render_json_file(maintenance["lane_heartbeats"]))
    safe_write(run_dir, "lane_heartbeats.md", render_lane_heartbeats(maintenance["lane_heartbeats"]))
    write_lane_outputs(run_dir, maintenance["lane_heartbeats"])
    safe_write(run_dir, "cleanup_compaction_index.json", render_json_file(maintenance["cleanup_compaction"]))
    safe_write(run_dir, "cleanup_compaction_index.md", render_cleanup_compaction(maintenance["cleanup_compaction"]))
    safe_write(run_dir, "internal_evolution_queue.json", render_json_file(maintenance["internal_evolution_queue"]))
    safe_write(run_dir, "internal_evolution_queue.md", render_internal_evolution_queue(maintenance["internal_evolution_queue"]))
    safe_write(run_dir, "sdcc_base_repo_hygiene.md", render_sdcc_base_repo_hygiene(maintenance["sdcc_base_repo_hygiene"]))
    safe_write(run_dir, "autonomy_envelope.json", render_json_file(autonomous["autonomy_envelope"]))
    safe_write(run_dir, "autonomy_envelope.md", render_autonomy_envelope(autonomous["autonomy_envelope"]))
    safe_write(run_dir, "changed_actions.json", render_json_file(autonomous["changed_actions"]))
    safe_write(run_dir, "changed_actions.md", render_changed_actions(autonomous["changed_actions"]))
    safe_write(run_dir, "deferred_actions.json", render_json_file(autonomous["deferred_actions"]))
    safe_write(run_dir, "deferred_for_justin.md", render_deferred_for_justin(autonomous["deferred_actions"]))
    safe_write(run_dir, "research_synthesis.md", render_research_synthesis(autonomous["research_scout"]))
    safe_write(run_dir, "validation_summary.md", render_validation_summary(autonomous["validation_summary"]))
    safe_write(run_dir, "morning_digest.md", render_morning_digest(autonomous["morning_digest"]))
    safe_write(run_dir, "structural_advancements.json", render_json_file(structural))
    safe_write(run_dir, "structural_advancements.md", render_structural_advancements(structural))
    safe_write(run_dir, "morning_conversation_digest.md", render_morning_conversation_digest(run_id, world, risks, maintenance, autonomous))
    safe_write(run_dir, "path_contract_drift.md", render_path_contract_drift(world, maintenance))
    safe_write(run_dir, "promotion_candidates.md", render_promotion_candidates(maintenance))
    safe_write(run_dir, "morning_approval_queue.md", render_morning_approval_queue(maintenance))
    safe_write(run_dir, "approval_packets/SHADOW_README.md", render_shadow_approval_packet_readme(maintenance))

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


def lane_risks_for(risks: list[dict[str, Any]], lane_id: str) -> list[dict[str, Any]]:
    return [item for item in risks if not risk_is_internal_maintenance_noise(item) and lane_id_for_risk(item) == lane_id]


def lane_required_artifacts(lane_id: str) -> list[str]:
    mapping = {
        "status_provenance": ["targeted_findings.json", "runbook_proposals.md"],
        "architecture_cohesion": ["structural_advancements.md", "internal_evolution_queue.md", "internal_evolution_queue.json"],
        "context_resonance": [
            "context_resonance_review.md",
            "context_resonance_review.json",
            "latent_context_nudges.md",
            "latent_context_nudges.json",
            "negative_memory_candidates.md",
            "negative_memory_candidates.json",
            "retrieval_outcome_candidates.md",
            "retrieval_outcome_candidates.json",
        ],
        "cleanup_storage": [
            "cleanup_proposals.md",
            "cleanup_compaction_index.md",
            "cleanup_compaction_index.json",
            "changed_actions.md",
            "changed_actions.json",
        ],
        "path_contract": ["path_contract_drift.md", "sdcc_base_repo_hygiene.md"],
        "research_scout": ["research_synthesis.md", "human_tool_leverage.md"],
        "science_scout": ["physics_scenario_proposals.md", "ideal_final_figure_gallery/xjgamma_modification_target.svg"],
        "presentation_artifacts": ["figure_design_notes.md", "daily_plan_proposals.md"],
    }
    return mapping[lane_id]


def render_lane_validation_summary(lane_id: str) -> str:
    payload = {
        "status": "lane_commands_listed_for_validation",
        "commands": [
            f"python3 scripts/codex_os_dream.py lane --lane-id {lane_id}",
            "python3 scripts/codex_os_dream.py validate --latest",
            "python3 scripts/codex_os_doctor.py --profile daily",
        ],
    }
    return render_validation_summary(payload)


def render_lane_heartbeat_report(
    run_id: str,
    lane_id: str,
    world: dict[str, Any],
    lane_risks: list[dict[str, Any]],
    lane_signal_payload: dict[str, Any],
) -> str:
    lane = lane_definition_by_id(lane_id)
    research_payload = lane_signal_payload.get("research_scout") if lane_id == "research_scout" and isinstance(lane_signal_payload.get("research_scout"), dict) else {}
    research_subtasks = []
    if research_payload:
        leverage = research_payload.get("human_tool_leverage") if isinstance(research_payload.get("human_tool_leverage"), dict) else {}
        research_subtasks = leverage.get("first_class_subtasks") if isinstance(leverage.get("first_class_subtasks"), list) else []
    resonance_payload = (
        lane_signal_payload.get("context_resonance")
        if lane_id == "context_resonance" and isinstance(lane_signal_payload.get("context_resonance"), dict)
        else {}
    )
    resonance_findings = []
    if resonance_payload:
        resonance_findings.extend(resonance_payload.get("latent_context_nudges") or [])
        resonance_findings.extend(resonance_payload.get("negative_memories") or [])
    lines = [f"# Dream Lane Heartbeat `{lane_id}`", "", SYNTHETIC_HEADER, ""]
    lines.append(f"- run: `{run_id}`")
    lines.append(f"- lane: `{lane_id}`")
    lines.append(f"- role: {lane.get('agent_role')}")
    lines.append(f"- purpose: {lane.get('purpose')}")
    lines.append(f"- generated: {world.get('generated_at')}")
    lines.append(f"- automation_id: `{lane_signal_payload.get('automation_id')}`")
    lines.append(f"- thread_binding: `{lane_signal_payload.get('thread_binding')}`")
    if lane_signal_payload.get("target_thread_id"):
        lines.append(f"- target_thread_id: `{lane_signal_payload.get('target_thread_id')}`")
    lines.append(f"- lane_risk_count: {len(lane_risks)}")
    lines.append("")
    lines.append("## Top Findings")
    if lane_risks:
        for item in lane_risks[:5]:
            target = item.get("workstream_id") or item.get("evidence") or item.get("title") or "global"
            lines.append(f"- score {item.get('score')} `{item.get('kind')}` {target}: {item.get('label')}")
    elif research_subtasks:
        for item in research_subtasks[:5]:
            lines.append(f"- `{item.get('id')}`: {item.get('goal')}")
    elif resonance_findings:
        for item in resonance_findings[:5]:
            title = item.get("id") or item.get("memory_id") or item.get("source") or "context_signal"
            detail = item.get("why_it_surfaced") or item.get("trap") or item.get("reason") or item.get("required_waking_check")
            lines.append(f"- `{title}`: {detail}")
    else:
        lines.append("- none")
    lines.append("")
    lines.append("## Morning Actions")
    if lane_risks:
        for item in lane_risks[:3]:
            lines.append(f"- {dream_response_for(item)}")
    elif research_subtasks:
        for item in research_subtasks[:3]:
            lines.append(f"- {item.get('morning_action')}")
    elif resonance_findings:
        for item in resonance_findings[:3]:
            lines.append(f"- {item.get('required_waking_check') or item.get('first_safe_action') or 'Review the context resonance artifact before using this signal.'}")
    else:
        lines.append("- No lane-local waking action stands out from the current evidence.")
    lines.append("")
    if research_subtasks:
        lines.append("## Human-Tool Leverage Subtasks")
        for item in research_subtasks[:3]:
            lines.append(f"- `{item.get('id')}`: {item.get('goal')}")
            lines.append(f"  next: {item.get('morning_action')}")
        lines.append("")
    lines.append("## Boundary")
    lines.append("- Proposal only. No SDCC, Condor, Gmail, Google Drive/Slides, Linear, or repo-tracked mutation.")
    return "\n".join(lines) + "\n"


def render_lane_digest(lane_id: str, lane_risks: list[dict[str, Any]], maintenance: dict[str, Any], autonomous: dict[str, Any]) -> str:
    lane = lane_definition_by_id(lane_id)
    queue = maintenance.get("internal_evolution_queue") if isinstance(maintenance.get("internal_evolution_queue"), dict) else {}
    queue_rows = queue.get("rows") if isinstance(queue.get("rows"), list) else []
    research_payload = autonomous.get("research_scout") if lane_id == "research_scout" and isinstance(autonomous.get("research_scout"), dict) else {}
    research_subtasks = []
    if research_payload:
        leverage = research_payload.get("human_tool_leverage") if isinstance(research_payload.get("human_tool_leverage"), dict) else {}
        research_subtasks = leverage.get("first_class_subtasks") if isinstance(leverage.get("first_class_subtasks"), list) else []
    resonance_payload = (
        maintenance.get("context_resonance")
        if lane_id == "context_resonance" and isinstance(maintenance.get("context_resonance"), dict)
        else {}
    )
    resonance_findings = []
    if resonance_payload:
        resonance_findings.extend(resonance_payload.get("latent_context_nudges") or [])
        resonance_findings.extend(resonance_payload.get("negative_memories") or [])
    lines = [f"# Lane Digest `{lane_id}`", "", SYNTHETIC_HEADER, ""]
    lines.append(f"- role: {lane.get('agent_role')}")
    lines.append(f"- purpose: {lane.get('purpose')}")
    lines.append(f"- risk_count: {len(lane_risks)}")
    lines.append("")
    lines.append("## What Changed")
    lines.append("- This lane writes only local dream artifacts and proposal surfaces.")
    lines.append("")
    lines.append("## Top 3 Lane Findings")
    if lane_risks:
        for item in lane_risks[:3]:
            lines.append(f"- `{item.get('kind')}`: {item.get('label')}")
    elif research_subtasks:
        for item in research_subtasks[:3]:
            lines.append(f"- `{item.get('id')}`: {item.get('goal')}")
    elif resonance_findings:
        for item in resonance_findings[:3]:
            title = item.get("id") or item.get("memory_id") or item.get("source") or "context_signal"
            detail = item.get("why_it_surfaced") or item.get("trap") or item.get("reason")
            lines.append(f"- `{title}`: {detail}")
    else:
        lines.append("- none")
    lines.append("")
    lines.append("## Top 3 Proposed Morning Actions")
    if lane_risks:
        for item in lane_risks[:3]:
            lines.append(f"- {validator_candidate_for(item).get('safe_command')}")
    elif research_subtasks:
        for item in research_subtasks[:3]:
            lines.append(f"- {item.get('morning_action')}")
    elif resonance_findings:
        for item in resonance_findings[:3]:
            lines.append(f"- {item.get('required_waking_check') or item.get('first_safe_action') or 'Record whether this signal helped or polluted context after waking use.'}")
    else:
        lines.append("- No lane-local waking action was proposed.")
    lines.append("")
    if research_subtasks:
        lines.append("## First-Class Research Scout Subtasks")
        for item in research_subtasks[:3]:
            lines.append(f"- `{item.get('id')}`: {item.get('goal')}")
            lines.append(f"  next: {item.get('morning_action')}")
        lines.append("")
    lines.append("## Best Next Internal Fix")
    matching = [row for row in queue_rows if row.get("target") == lane_id or lane_id in str(row.get("target") or "")]
    source = matching[0] if matching else (queue_rows[0] if queue_rows else None)
    if source:
        lines.append(f"- `{source.get('target')}`: {source.get('action')}")
    else:
        lines.append("- No ranked lane-local internal fix is currently queued.")
    return "\n".join(lines) + "\n"


def write_lane_artifacts(
    run_dir: Path,
    lane_id: str,
    world: dict[str, Any],
    lane_risks: list[dict[str, Any]],
    maintenance: dict[str, Any],
    autonomous: dict[str, Any],
    structural: dict[str, Any],
) -> None:
    if lane_id == "status_provenance":
        rows = [item for item in maintenance.get("targeted_findings") or [] if lane_id_for_risk(item) == lane_id]
        safe_write(
            run_dir,
            "targeted_findings.json",
            render_json_file(
                {
                    "synthetic": True,
                    "synthetic_header": SYNTHETIC_HEADER,
                    "mutation_boundary": "proposal_only",
                    "rows": rows,
                }
            ),
        )
        safe_write(run_dir, "runbook_proposals.md", render_runbook_proposals(lane_risks))
    elif lane_id == "architecture_cohesion":
        safe_write(run_dir, "structural_advancements.md", render_structural_advancements(structural))
        safe_write(run_dir, "internal_evolution_queue.md", render_internal_evolution_queue(maintenance["internal_evolution_queue"]))
        safe_write(run_dir, "internal_evolution_queue.json", render_json_file(maintenance["internal_evolution_queue"]))
    elif lane_id == "context_resonance":
        resonance = maintenance["context_resonance"]
        safe_write(run_dir, "context_resonance_review.md", render_context_resonance_markdown(resonance))
        safe_write(run_dir, "context_resonance_review.json", render_json_file(resonance))
        safe_write(
            run_dir,
            "latent_context_nudges.md",
            render_context_resonance_list("Latent Context Nudges", resonance.get("latent_context_nudges") or []),
        )
        safe_write(
            run_dir,
            "latent_context_nudges.json",
            render_json_file(
                {
                    "synthetic": True,
                    "synthetic_header": SYNTHETIC_HEADER,
                    "mutation_boundary": "proposal_only",
                    "rows": resonance.get("latent_context_nudges") or [],
                }
            ),
        )
        safe_write(
            run_dir,
            "negative_memory_candidates.md",
            render_context_resonance_list("Negative Memory Candidates", resonance.get("negative_memories") or []),
        )
        safe_write(
            run_dir,
            "negative_memory_candidates.json",
            render_json_file(
                {
                    "synthetic": True,
                    "synthetic_header": SYNTHETIC_HEADER,
                    "mutation_boundary": "proposal_only",
                    "rows": resonance.get("negative_memories") or [],
                }
            ),
        )
        safe_write(
            run_dir,
            "retrieval_outcome_candidates.md",
            render_context_resonance_list("Retrieval Outcome Candidates", resonance.get("retrieval_outcome_candidates") or []),
        )
        safe_write(
            run_dir,
            "retrieval_outcome_candidates.json",
            render_json_file(
                {
                    "synthetic": True,
                    "synthetic_header": SYNTHETIC_HEADER,
                    "mutation_boundary": "proposal_only",
                    "rows": resonance.get("retrieval_outcome_candidates") or [],
                }
            ),
        )
    elif lane_id == "cleanup_storage":
        safe_write(run_dir, "cleanup_proposals.md", render_cleanup_proposals(world, lane_risks))
        safe_write(run_dir, "cleanup_compaction_index.md", render_cleanup_compaction(maintenance["cleanup_compaction"]))
        safe_write(run_dir, "cleanup_compaction_index.json", render_json_file(maintenance["cleanup_compaction"]))
        safe_write(run_dir, "changed_actions.md", render_changed_actions(autonomous["changed_actions"]))
        safe_write(run_dir, "changed_actions.json", render_json_file(autonomous["changed_actions"]))
    elif lane_id == "path_contract":
        safe_write(run_dir, "path_contract_drift.md", render_path_contract_drift(world, maintenance))
        safe_write(run_dir, "sdcc_base_repo_hygiene.md", render_sdcc_base_repo_hygiene(maintenance["sdcc_base_repo_hygiene"]))
    elif lane_id == "research_scout":
        safe_write(run_dir, "research_synthesis.md", render_research_synthesis(autonomous["research_scout"]))
        safe_write(run_dir, "human_tool_leverage.md", render_human_tool_leverage(autonomous["research_scout"]))
    elif lane_id == "science_scout":
        safe_write(run_dir, "physics_scenario_proposals.md", render_physics_scenario_proposals(world, lane_risks))
        for spec in IDEAL_FIGURE_SPECS:
            safe_write(run_dir, f"ideal_final_figure_gallery/{spec['id']}.svg", render_target_figure_svg(spec))
    elif lane_id == "presentation_artifacts":
        safe_write(run_dir, "figure_design_notes.md", render_figure_design_notes())
        safe_write(run_dir, "daily_plan_proposals.md", render_daily_proposals(lane_risks))


def write_lane_dream_outputs(
    run_dir: Path,
    run_id: str,
    lane_id: str,
    world: dict[str, Any],
    risks: list[dict[str, Any]],
) -> None:
    copy_sandbox(run_dir)
    maintenance = evolutionary_maintenance_payload(world, risks)
    autonomous = autonomous_maintenance_payload(run_dir, world, risks, maintenance)
    structural = structural_advancement_payload(world, risks, maintenance)
    maintenance["structural_advancements"] = structural
    maintenance["autonomous_maintenance"] = autonomous
    lane_risks = lane_risks_for(risks, lane_id)
    interactions = generate_interactions(lane_risks, "lane")
    signal = heartbeat_signal(run_id, "lane", world, risks, maintenance)
    signal["version"] = 8
    signal["output_kind"] = "lane"
    signal["lane_id"] = lane_id
    signal["lane"] = lane_definition_by_id(lane_id)
    signal["automation_id"] = DREAM_LANE_AUTOMATIONS[lane_id]["automation_id"]
    signal["thread_binding"] = DREAM_LANE_AUTOMATIONS[lane_id].get("thread_binding", "fresh_chat_per_run")
    if DREAM_LANE_AUTOMATIONS[lane_id].get("target_thread_id"):
        signal["target_thread_id"] = DREAM_LANE_AUTOMATIONS[lane_id]["target_thread_id"]
    signal["summary"]["lane_risk_count"] = len(lane_risks)
    signal["top_findings"] = [
        {
            "kind": item.get("kind"),
            "score": item.get("score"),
            "workstream_id": item.get("workstream_id"),
            "title": item.get("title"),
            "evidence": item.get("evidence"),
            "validator_kind": validator_candidate_for(item).get("kind"),
        }
        for item in lane_risks[:5]
    ]
    signal["changed_actions"] = {
        "external_mutations_performed": False,
        "science_mutations_performed": False,
        "repo_tracked_mutations_performed": False,
    }
    if lane_id == "research_scout":
        signal["research_scout"] = autonomous["research_scout"]
        research_payload = autonomous["research_scout"]
        leverage = research_payload.get("human_tool_leverage") if isinstance(research_payload.get("human_tool_leverage"), dict) else {}
        subtasks = leverage.get("first_class_subtasks") if isinstance(leverage.get("first_class_subtasks"), list) else []
        if subtasks:
            signal["top_findings"] = [
                {
                    "kind": "research_scout_subtask",
                    "score": None,
                    "workstream_id": "research_scout",
                    "title": item.get("id"),
                    "evidence": item.get("goal"),
                    "validator_kind": "delegated_research_review",
                }
                for item in subtasks[:5]
            ]
    if lane_id == "context_resonance":
        resonance_payload = maintenance["context_resonance"]
        signal["context_resonance"] = resonance_payload
        nudges = resonance_payload.get("latent_context_nudges") if isinstance(resonance_payload.get("latent_context_nudges"), list) else []
        negative = resonance_payload.get("negative_memories") if isinstance(resonance_payload.get("negative_memories"), list) else []
        retrieval = (
            resonance_payload.get("retrieval_outcome_candidates")
            if isinstance(resonance_payload.get("retrieval_outcome_candidates"), list)
            else []
        )
        signal["summary"]["context_resonance_nudge_count"] = len(nudges)
        signal["summary"]["negative_memory_candidate_count"] = len(negative)
        signal["summary"]["retrieval_outcome_candidate_count"] = len(retrieval)
        signal["top_findings"] = [
            {
                "kind": "latent_context_nudge",
                "score": None,
                "workstream_id": "context_resonance",
                "title": item.get("id"),
                "evidence": item.get("why_it_surfaced"),
                "validator_kind": "required_waking_check",
            }
            for item in nudges[:3]
        ] + [
            {
                "kind": "negative_memory_candidate",
                "score": None,
                "workstream_id": "context_resonance",
                "title": item.get("id"),
                "evidence": item.get("trap"),
                "validator_kind": "negative_memory_review",
            }
            for item in negative[:2]
        ]
    maintenance_payload = signal.get("evolutionary_maintenance") if isinstance(signal.get("evolutionary_maintenance"), dict) else {}
    lane_payload = maintenance_payload.get("lane_heartbeats") if isinstance(maintenance_payload.get("lane_heartbeats"), dict) else {}
    if lane_payload:
        lane_rows = [item for item in lane_payload.get("lanes") or [] if item.get("lane_id") == lane_id]
        lane_payload["purpose"] = "lane-local heartbeat payload"
        lane_payload["engineering_rule"] = "one lane per automation; each nightly run opens a fresh automation chat; doctor aggregates the set"
        lane_payload["lanes"] = lane_rows
        lane_payload.pop("master_heartbeat", None)
    trace = {
        "synthetic": True,
        "synthetic_header": SYNTHETIC_HEADER,
        "mutation_boundary": "proposal_only",
        "output_kind": "lane",
        "lane_id": lane_id,
        "run_id": run_id,
        "mode": "lane",
        "world": {
            "generated_at": world["generated_at"],
            "workstream_counts": world["workstream_counts"],
            "artifact_counts": world["artifact_counts"],
            "artifact_gaps": world["artifact_gaps"],
        },
        "risks": lane_risks,
        "interactions": interactions,
        "lane_signal": signal,
    }
    safe_write(run_dir, "dream_trace.json", render_json_file(trace))
    safe_write(run_dir, "lane_signal.json", render_json_file(signal))
    safe_write(run_dir, "lane_heartbeat.md", render_lane_heartbeat_report(run_id, lane_id, world, lane_risks, signal))
    safe_write(run_dir, "lane_digest.md", render_lane_digest(lane_id, lane_risks, maintenance, autonomous))
    safe_write(run_dir, "validation_summary.md", render_lane_validation_summary(lane_id))
    write_lane_artifacts(run_dir, lane_id, world, lane_risks, maintenance, autonomous, structural)
    append_dream_index(signal)


def run_dream(args: argparse.Namespace) -> int:
    raise SystemExit(
        "ERROR: `python3 scripts/codex_os_dream.py "
        f"{args.mode}` is retired. Use `python3 scripts/codex_os_dream.py lane --lane-id <lane_id>`."
    )


def run_lane_dream(args: argparse.Namespace) -> int:
    now = parse_now(args.now)
    run_id = args.run_id or lane_run_id_for(args.lane_id, now)
    run_dir = DREAM_ROOT / run_id
    register_data = load_yaml_like(Path(args.register))
    artifact_data = load_yaml_like(ARTIFACT_REGISTRY)
    world = build_world(register_data, artifact_data, now, "lane", run_id)
    risks = identify_risks(world, "lane")
    write_lane_dream_outputs(run_dir, run_id, args.lane_id, world, risks)
    lane_risks = lane_risks_for(risks, args.lane_id)
    print(f"OK: lane dream wrote {run_dir}")
    print(f"lane_id={args.lane_id} risks={len(lane_risks)}")
    print(f"report={run_dir / 'lane_heartbeat.md'}")
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
        "targeted_findings.json",
        "thesis_flow_scores.json",
        "branch_pressure_report.json",
        "memory_homeostasis_report.json",
        "controlled_burn_proposals.yaml",
        "lane_heartbeats.json",
        "lane_heartbeats.md",
        "cleanup_compaction_index.json",
        "cleanup_compaction_index.md",
        "internal_evolution_queue.json",
        "internal_evolution_queue.md",
        "sdcc_base_repo_hygiene.md",
        "autonomy_envelope.json",
        "autonomy_envelope.md",
        "changed_actions.json",
        "changed_actions.md",
        "deferred_actions.json",
        "deferred_for_justin.md",
        "research_synthesis.md",
        "validation_summary.md",
        "morning_digest.md",
        "structural_advancements.json",
        "structural_advancements.md",
        "morning_conversation_digest.md",
        "path_contract_drift.md",
        "promotion_candidates.md",
        "morning_approval_queue.md",
        "approval_packets/SHADOW_README.md",
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
        "lane_heartbeats.md",
        "cleanup_compaction_index.md",
        "internal_evolution_queue.md",
        "sdcc_base_repo_hygiene.md",
        "autonomy_envelope.md",
        "changed_actions.md",
        "deferred_for_justin.md",
        "research_synthesis.md",
        "validation_summary.md",
        "morning_digest.md",
        "structural_advancements.md",
        "morning_conversation_digest.md",
        "path_contract_drift.md",
        "promotion_candidates.md",
        "morning_approval_queue.md",
        "approval_packets/SHADOW_README.md",
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
    maintenance = heartbeat.get("evolutionary_maintenance")
    if not isinstance(maintenance, dict):
        errors.append("heartbeat_signal.json missing evolutionary_maintenance")
        maintenance = {}
    pilot = maintenance.get("shadow_pilot") if isinstance(maintenance.get("shadow_pilot"), dict) else {}
    if pilot.get("name") != "four_night_shadow_pilot":
        errors.append("evolutionary_maintenance missing four-night shadow pilot")
    if pilot.get("approval_ready_allowed") is not False:
        errors.append("four-night shadow pilot approval_ready_allowed is not false")
    summary = maintenance.get("summary") if isinstance(maintenance.get("summary"), dict) else {}
    if int(summary.get("approval_ready_count") or 0) != 0:
        errors.append("evolutionary maintenance emitted approval-ready items during shadow pilot")
    lanes = maintenance.get("lane_heartbeats") if isinstance(maintenance.get("lane_heartbeats"), dict) else {}
    if len(lanes.get("lanes") or []) < len(DREAM_LANE_DEFINITIONS):
        errors.append("evolutionary maintenance missing dream lane heartbeats")
    for lane in lanes.get("lanes") or []:
        lane_id = str(lane.get("lane_id") or "unknown_lane")
        safe_lane_id = re.sub(r"[^a-zA-Z0-9_.-]+", "_", lane_id).strip("_") or "unknown_lane"
        if not (path / "lanes" / safe_lane_id / "heartbeat.md").exists():
            errors.append(f"missing lane heartbeat markdown for {lane_id}")
        if not (path / "lanes" / safe_lane_id / "heartbeat.json").exists():
            errors.append(f"missing lane heartbeat json for {lane_id}")
    for name in (
        "targeted_findings.json",
        "thesis_flow_scores.json",
        "branch_pressure_report.json",
        "memory_homeostasis_report.json",
        "lane_heartbeats.json",
        "cleanup_compaction_index.json",
        "internal_evolution_queue.json",
        "autonomy_envelope.json",
        "changed_actions.json",
        "deferred_actions.json",
        "structural_advancements.json",
    ):
        try:
            payload = json.loads(read_text(path / name))
        except json.JSONDecodeError as exc:
            errors.append(f"{name} invalid JSON: {exc}")
            continue
        if not isinstance(payload, dict) or payload.get("synthetic") is not True:
            errors.append(f"{name} missing synthetic=true")
        allowed_boundaries = {"proposal_only", "local_internal_auto_safe"}
        if isinstance(payload, dict) and payload.get("mutation_boundary") not in allowed_boundaries:
            errors.append(f"{name} mutation_boundary is not an allowed dream boundary")
        if isinstance(payload, dict) and payload.get("synthetic_header") != SYNTHETIC_HEADER:
            errors.append(f"{name} missing synthetic header")
        if name == "changed_actions.json" and isinstance(payload, dict):
            if payload.get("external_mutations_performed") is not False:
                errors.append("changed_actions.json reports external mutation")
            if payload.get("science_mutations_performed") is not False:
                errors.append("changed_actions.json reports science mutation")
            if payload.get("repo_tracked_mutations_performed") is not False:
                errors.append("changed_actions.json reports repo-tracked mutation")
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


def validate_lane_dir(path: Path) -> list[str]:
    errors: list[str] = []
    root = DREAM_ROOT.resolve()
    resolved = path.resolve()
    if root not in [resolved, *resolved.parents]:
        errors.append(f"dream path is outside dream root: {path}")
        return errors
    required = ["dream_trace.json", "lane_signal.json", "lane_heartbeat.md", "lane_digest.md", "validation_summary.md"]
    for name in required:
        if not (path / name).exists():
            errors.append(f"missing lane output: {name}")
    try:
        signal = json.loads(read_text(path / "lane_signal.json"))
    except json.JSONDecodeError:
        signal = {}
    if not isinstance(signal, dict) or not signal:
        errors.append("lane_signal.json invalid JSON")
        return errors
    lane_id = str(signal.get("lane_id") or "")
    if lane_id not in DREAM_LANE_AUTOMATIONS:
        errors.append(f"lane_signal.json has unknown lane_id: {lane_id}")
        return errors
    for name in required[2:]:
        if SYNTHETIC_HEADER not in read_text(path / name):
            errors.append(f"{name} lacks synthetic provenance header")
    if signal.get("synthetic") is not True:
        errors.append("lane_signal.json synthetic flag is not true")
    if signal.get("mutation_boundary") != "proposal_only":
        errors.append("lane_signal.json mutation_boundary is not proposal_only")
    if signal.get("output_kind") != "lane":
        errors.append("lane_signal.json output_kind is not lane")
    contract = DREAM_LANE_AUTOMATIONS[lane_id]
    if signal.get("automation_id") != contract["automation_id"]:
        errors.append("lane_signal.json automation_id does not match lane contract")
    expected_binding = contract.get("thread_binding") or ("fixed_thread" if contract.get("target_thread_id") else "fresh_chat_per_run")
    if signal.get("thread_binding") != expected_binding:
        errors.append("lane_signal.json thread_binding does not match lane contract")
    expected_thread_id = contract.get("target_thread_id")
    if expected_thread_id:
        if signal.get("target_thread_id") != expected_thread_id:
            errors.append("lane_signal.json target_thread_id does not match lane contract")
    elif signal.get("target_thread_id") not in {None, ""}:
        errors.append("lane_signal.json unexpectedly hard-binds a target_thread_id")
    for name in lane_required_artifacts(lane_id):
        if not (path / name).exists():
            errors.append(f"missing lane-specific output: {name}")
    for name in ("nightly_heartbeat.md", "nightly_heartbeat_signal.json", "morning_appendix.md"):
        if (path / name).exists():
            errors.append(f"retired master-heartbeat output still present: {name}")
    return errors


def command_validate(args: argparse.Namespace) -> int:
    path = latest_dream_dir() if args.latest else Path(args.path)
    if path is None:
        if args.allow_missing:
            print("OK: no dream directory exists yet")
            return 0
        print("ERROR: no dream directory exists", file=sys.stderr)
        return 1
    errors = validate_lane_dir(path) if (path / "lane_signal.json").exists() else validate_dream_dir(path)
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

    lane = subparsers.add_parser("lane", help="run one specific dream lane")
    lane.add_argument("--lane-id", required=True, choices=[item["lane_id"] for item in DREAM_LANE_DEFINITIONS])
    lane.add_argument("register", nargs="?", default=str(DEFAULT_REGISTER))
    lane.add_argument("--now")
    lane.add_argument("--run-id")
    lane.set_defaults(func=run_lane_dream)

    validate = subparsers.add_parser("validate", help="validate dream output safety")
    validate.add_argument("path", nargs="?", default="")
    validate.add_argument("--latest", action="store_true")
    validate.add_argument("--allow-missing", action="store_true")
    validate.set_defaults(func=command_validate)

    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
