# Agent Policy Index

This directory turns the old monolithic `AGENTS.md` into a traversable memory
system for future Codex sessions.

## Design Principle

Root `AGENTS.md` is working memory: short, always loaded, full of hard stops
and pointers. This directory is long-term procedural memory: load only the
policy needed for the current task.

Useful memory concepts from a quick research pass:

- Working memory should be compact and task-facing.
- Long-term memory should be structured, searchable, and split by retrieval
  intent.
- Good agent memory is a write-manage-read loop: decide what is worth saving,
  keep it in the right store, retrieve it at the right moment, and forget or
  archive stale material.
- Entity/episodic memory matters for this project: campaign tags, ROOT paths,
  cluster IDs, Gmail stage messages, plot products, and slide IDs are entities;
  failed runs and successful recovery commands are episodes.
- Structured retrieval beats vague narrative. Prefer tables, YAML load maps,
  path ledgers, and short checklists.

## Pointer Types

- **Hot state**: live truth that changes often.
  `agent_context/TASK_BOARD.md`, `agent_context/STATUS_DASHBOARD.md`,
  `codex_notes/*.md`.
- **Policy**: stable operating rules. Files in this directory.
- **Evidence pointer**: exact path, tag, cluster ID, email subject, timestamp,
  ROOT object, plot PNG, or command output that proves a claim.
- **Run fingerprint**: fields that identify a training/production/merge so
  Codex can detect duplicates before spending compute.
- **Transfer pointer**: local path plus intended SDCC path and rebuild status.
- **Slide pointer**: deck ID, slide object ID, element object ID, thumbnail, and
  source plot path.

## Files

- `LOAD_MAP.yaml`: machine-readable-ish route table and preflight loads.
- `AGENT_MEMORY_DESIGN_NOTES.md`: research-grounded design rationale and
  pointer taxonomy.
- `LEGACY_AGENTS_COVERAGE.md`: old monolith heading-to-policy preservation map.
- `HARD_STOPS_AND_SAFETY.md`: non-negotiable stop rules.
- `DUPLICATE_RUN_GUARD.md`: duplicate detection and run fingerprints.
- `COLLABORATION_AND_HANDOFFS.md`: orchestration and delegated-lane handoffs.
- `SCIENCE_REFERENCE_HIERARCHY.md`: PPG/ATLAS science north star.
- `SDCC_OPERATIONS.md`: SSH, Condor, submissions, diagnostics.
- `EMAIL_AND_WATCHDOGS.md`: Gmail pipeline and automation rules.
- `MEMORY_AND_STATUS.md`: how to record durable state.
- `INPUTS_AND_ENVIRONMENT.md`: ROOT/env/input-file conventions.
- `PLOTTING.md`: plot generation and visual QA.
- `SLIDES_WORKFLOW.md`: Google Slides and deck style rules.
- `STORAGE_AND_CLEANUP.md`: cleanup, quotas, lifecycle.
- `TOOLS_AND_DEPENDENCIES.md`: missing-tool, runtime, and connector routing.
- `TRANSFER_AND_PIPELINE_FILES.md`: SFTP helpers, mapped files, upload rules.
- `CURRENT_PROJECT_STATE.md`: durable stale-output and current-analysis state.

## How Future Codex Should Use This

1. Read root `AGENTS.md`.
2. Open `LOAD_MAP.yaml`.
3. Load only the 1-4 policy files relevant to the user request.
4. Check hot state and evidence before acting.
5. After a consequential run, plot set, or slide update, write the new evidence
   back into the appropriate hot-state ledger.
