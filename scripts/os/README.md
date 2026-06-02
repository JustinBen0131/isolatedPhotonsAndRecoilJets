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
check routing caps and synthetic suppression. Use `feedback-health --json` to
verify the adaptive feedback loop itself: parseable ledger rows, ignored
malformed/unknown/synthetic rows, bounded score adjustments, route-policy
loading, missed-feedback review-only behavior, and local proposal-only
homeostasis. Use `maintenance-review --json` to summarize local retrieval
pressure for the context-resonance maintenance selector.

The default `python3 scripts/codex_os_dream.py lane --lane-id context_resonance`
run now selects at most one auto-safe local memory-maintenance action. Allowed
automatic writes are limited to ignored state under
`agent_context/local/context_resonance/` and the current ignored dream run
directory; tracked registries and policies remain waking-review proposals. Use
`--no-auto-maintain` for a debug run that emits the candidate and validators
without applying the local change. If a candidate is blocked, the lane records
the intended update, exact failure reason, and best next step, then tries the
next candidate in the priority ladder until one validates or the run emits a
single no-change/deferred summary. The lane writes `changed_actions.md/json`
and, when it applies a change, `rollback_manifest.json`. It also writes compact
`feedback_loop_health.md/json` artifacts every run so the morning summary can
say whether feedback is improving retrieval quality without synthetic
contamination, stale suppression drift, or uncontrolled memory mutation. Old
resonance report files are compatibility views derived from the same lane
signal.
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

For Justin's explicit dual-Pro overnight research lane, use
`python3 scripts/os/research/codex_chatgpt_research_pack.py init-dual-pro-overnight`
to stage two fresh-conversation ChatGPT Pro extended packs: one physics thesis
prompt and one OS-infrastructure prompt. The generated manifest also names the
single temporary heartbeat command, response-ingestion commands, daily-note
path, and git precondition. The shell dream script only emits the plan; a
Codex automation using Computer Use performs authenticated ChatGPT UI actions.

Use `python3 scripts/os/register/codex_work_register_stale.py
agent_context/CODEX_WORK_REGISTER.yaml --protocol` as the read-only
`register_workstream_refresh_contract`. It classifies live workstreams and
prevents repeated status/stale warnings from recurring as prose once the
protocol handles them.
