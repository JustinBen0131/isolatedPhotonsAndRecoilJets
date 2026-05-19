# Agent Memory Design Notes

## Research Patterns Folded In

These notes are for future Codex sessions trying to understand why the project
memory is structured this way.

Useful external patterns consulted on 2026-05-18:

- Microsoft AI Agents for Beginners, "Memory for AI Agents":
  `https://microsoft.github.io/ai-agents-for-beginners/13-agent-memory/`
- "Memory for Autonomous LLM Agents: Mechanisms, Evaluation, and Emerging
  Frontiers" arXiv overview:
  `https://arxiv.org/abs/2507.21046`

Key ideas adapted locally:

- working memory should stay small and task-local;
- long-term memory should be structured and explicitly retrievable;
- agent memory works best as a write/manage/read loop, not a pile of notes;
- memory types should be separated: policy, episodic run evidence, entity-like
  dataset/model facts, workflow habits, and cold-storage cleanup;
- forgetting and cleanup should be deliberate, evidence-backed actions.

## Local Implementation

- `AGENTS.md` is the short always-loaded working-memory control plane.
- `LOAD_MAP.yaml` is the router: classify request, load only relevant policy
  shards, then act.
- `TASK_BOARD.md` and `STATUS_DASHBOARD.md` are hot episodic state.
- `REFERENCE_MAP.md` and `SLIDE_STYLE_MAP.md` are stable entity/reference
  memory.
- `codex_notes/` remains compatibility memory for older durable notes.
- `DUPLICATE_RUN_GUARD.md` encodes run fingerprints so repeated compute can be
  detected before new submissions or analyses.

## Pointer Types

- Policy pointer: load this rule file before acting.
- Evidence pointer: path, timestamp, job ID, email subject, ROOT object, or
  command output that supports a claim.
- Run fingerprint: dataset, signal/background definition, model/features, cuts,
  tag, config, output root, source path, and objective.
- Slide pointer: deck ID, slide object ID, source PNG/artifact path, and visual
  approval state.
- Transfer pointer: local path, SDCC path, helper command, rebuild need, and
  last sync evidence.

## Maintenance Heuristic

Add a memory note only when it will change a future decision. If it is merely
interesting, put it in the chat; if it protects compute, prevents stale output,
or anchors a physics claim, give it a durable pointer.
