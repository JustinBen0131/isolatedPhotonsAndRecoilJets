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
- `agent_context/memory/` is the tracked neutral memory-routing layer:
  `CONTEXT_RESONANCE_INDEX.yaml` holds cue records, `SCHEMA_REGISTRY.yaml`
  holds repeated-method patterns, and `NEGATIVE_MEMORY_MAP.yaml` holds traps
  plus first safe actions.
- `scripts/os/context/codex_context_resonance.py` is the waking salience
  resolver: given a task, it returns a capped set of `conscious_context`
  entries, compatibility `active_facts`, latent nudges, negative-memory traps,
  suppressed context, and required evidence checks.
- `agent_context/local/context_resonance/` is the private home for memory
  salience, retrieval outcome, and context-interference ledgers. It should
  store distilled routing evidence, not raw transcripts.

## Context Resonance Layer

The project brain should not grow by stuffing more raw context into every boot.
It should improve the route between a task and the few old signals that matter.

Use these categories:

- `active_fact`: concise, directly relevant context that should be loaded or
  verified before acting.
- `latent_nudge`: a warning, analogy, method cue, visual-style cue, or path
  contract that may help but cannot support a claim by itself.
- `negative_memory`: a repeated trap to avoid, with a first safe action.
- `suppressed_context`: stale, synthetic, overly specific, or noisy material
  that should not surface by default.
- `retrieval_outcome`: after use, whether the surfaced item helped, polluted
  context, or should be cooled down.

Every resonance item should include `relation_type`, `evidence_class`,
`retrieval_policy`, and `required_waking_check`. Success is lower orientation
cost and fewer stale-context mistakes, not a larger memory bundle.

The primary waking command is:

```bash
python3 scripts/os/context/codex_context_resonance.py resolve --task "<task>" --json
```

After a retrieved cue materially affects a task, use `record-outcome` to append
local-only feedback as `useful`, `stale`, `harmful`, `irrelevant`, or `missed`.
This feedback tunes routing and suppression; it is not factual evidence.

Tracked registries must stay sanitized. They should contain pattern names,
source pointers, cue terms, salience values, evidence class, retrieval policy,
allowed use, and required waking checks. They should not contain raw ChatGPT
transcripts, private SDCC details, secrets, bulky chat logs, or unverified
external-model text.

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
