# Codex OS Helpers

Agent operating-system utilities organized into `safety/`, `register/`,
`dream/`, `heartbeat/`, `research/`, `context/`, and `artifacts/`. Top-level
`scripts/` entries are hard-contract compatibility aliases; edit the canonical
files in this tree first.

Lane dream output now includes `internal_evolution_queue.md`,
`internal_evolution_queue.json`, and `sdcc_base_repo_hygiene.md` under the
local ignored dream run directory when the relevant lane owns them. Use those
files as the first surface for Codex-brain cleanup recommendations, SDCC/base
read-only hygiene probes, and small reversible internal OS maintenance. They
are not approval to mutate SDCC, Condor, scientific outputs, Slides/Drive/
Gmail/Linear, or physics/task status.

The lane-first contract writes `lane_heartbeat.md`, `lane_signal.json`,
`lane_digest.md`, and `validation_summary.md` in every lane run directory.
Status/provenance, architecture cohesion, cleanup/storage, path contracts,
context resonance, research, science scouting, and presentation artifacts now
run as separate heartbeats that each open a fresh automation-created chat per
nightly run.

Use `python3 scripts/os/context/codex_context_resonance.py resolve --task "<task>" --json`
when a task needs compact retrieval guidance. It returns at most three active
facts, three latent nudges, negative-memory traps, suppressed stale/synthetic
context, and exact waking checks before any claim. The resolver is read-only;
the legacy flat `--task` interface still works for compatibility. Use
`record-outcome` to append local-only feedback after a memory helped or polluted
a real task, `review` to summarize that local feedback, and `canary --json` to
check routing caps and synthetic suppression. The `context_resonance` dream
lane writes proposal artifacts only under the local ignored dream run directory.
`cleanup_compaction_index.md/json` summarizes older local dream/research packs
with source pointers before any future archive/prune decision. Auto-safe dream
actions are limited to untracked generated junk and ignored local internal
artifacts; tracked OS changes are packaged for waking validation rather than
silently self-applied.
Use `morning_conversation_digest.md` as the first human-facing summary after a
dream run; it organizes what happened into did/learned/tried/changed/deferred/
next-fix sections. Even when no generated junk is removed, the digest must
surface cleanup review candidates, architecture hypotheses, and deeper search
angles instead of ending with "nothing to clean."
Use `structural_advancements.md` as the main improvement ledger: scores are
only diagnostics, while the ledger names the concrete validator, runbook,
index, retention, or path-contract improvement the dream thinks will make the
OS more traversable.

Use `python3 scripts/os/register/codex_work_register_stale.py
agent_context/CODEX_WORK_REGISTER.yaml --protocol` as the read-only
`register_workstream_refresh_contract`. It classifies live workstreams and
prevents repeated status/stale warnings from recurring as prose once the
protocol handles them.
