# Agentic OS Dreaming

This policy defines the ThesisAnalysis dream layer: private, synthetic,
proposal-only simulation that helps Codex rehearse likely future failures and
prepare better next-day responses without touching live external systems.

For a top-down blank-model handoff of the current pipeline, automations, lane
ownership, Justin/Codex responsibilities, and redundancy boundaries, read
`agent_context/policies/DREAM_SYSTEM_OPERATING_OVERVIEW.md`.

## Core Rule

Dreams are not user intent. A synthetic Justin line is a test fixture, never
approval, instruction, or evidence of what Justin actually said.

The dream system writes its durable output under:

```text
agent_context/local/dreams/
```

It may also perform tightly bounded local internal cleanup of untracked
generated junk such as `.DS_Store`, AppleDouble files, `__pycache__`, `.pyc`,
and `.pyo` artifacts. Every performed action must be logged in
`changed_actions.md/json` with reason, before/after summary, and rollback note.

It must not mutate SDCC, Condor, Gmail, Google Drive, Google Slides, Linear,
science outputs, task status, or repo-tracked files from the overnight script.
Tracked Codex-OS improvements may be packaged for waking Codex as
`auto_validated` maintenance only when rollback and validation commands are
explicit.

## Autonomous Internal Maintenance Envelope

Dream autonomy has four tiers:

- `auto_safe`: may run overnight. Local generated-junk cleanup, ignored dream
  and research-pack pruning, index refreshes, context compaction summaries,
  runbook cleanup notes, validator tuning candidates, and local-only audit
  artifacts.
- `auto_validated`: may be packaged overnight, but tracked repo changes are
  applied only by waking Codex after rollback, diff, and doctor validation.
- `research_only`: may prepare sanitized ChatGPT/Claude/Gemini-style prompts
  and local synthesis. External model output is critique/source leads, never
  truth, approval, or evidence.
- `blocked_for_waking`: SDCC runtime edits, Condor/job control, science output
  movement/deletion, TG bulk migration, Slides/Gmail/Linear/GitHub mutation,
  and physics/task-status changes.

Each nightly package must emit `autonomy_envelope`, `changed_actions`,
`deferred_for_justin`, `research_synthesis`, `validation_summary`,
`morning_digest`, `lane_heartbeats`, `structural_advancements`, and
`morning_conversation_digest` artifacts so automatic internal work is
inspectable after the fact. `structural_advancements.md` is the main
infrastructure-improvement ledger: it converts repeated pressure into concrete
validators, runbooks, indexes, retention rules, or path-contract improvements.
Recurring active-job/status/stale warnings should be retired into
`register_workstream_refresh_contract` once that protocol is present; doctor
warnings should count only unhandled recurring hotspots.
The pressure governor is the canonical rule for repeated dream pressure.
Dream internals, morning summaries, and doctor checks should normalize pressure
as `kind|target|source_lane` and assign exactly one status: `actionable`,
`handled`, `cooling`, `appendix_only`, `blocked_for_waking`, or `resolved`.
Daily-facing output should include only actual local changes, failures,
required waking approval, or unhandled actionable warnings. Repeated
proposal-only atoms, handled register pressure, context-resonance maintenance
hints, and research-scout subtasks move to the OS Appendix with
`daily_visibility: appendix` and a concrete `suppression_reason`. If a repeated
pressure is already covered by a validator, runbook, or register refresh
contract, render one line: `Handled by <contract>; no waking action.` Raw dream
artifacts stay preserved; suppression only changes daily visibility.
Scores are diagnostic smoke alarms only. Dream summaries are after-action
reports, not status blips. `morning_conversation_digest.md` and each scheduled
lane `morning_lane_summary.json` must say, in concrete scan-friendly language:
what the dream did, what surfaces it inspected, what it researched, what it
learned, what it changed locally, what it intentionally left untouched, what it
recommends changing only after waking approval, how long the lane took, which
artifacts support the report, and what feedback would make the next dream more
useful. A dream summary that only says "changed/deferred/no action" without
explaining what improved or what was learned is considered incomplete.

No dream may stop at "nothing to clean." If generated-junk cleanup is empty,
the dream must still produce at least one reviewable cleanup angle, architecture
improvement hypothesis, deeper search angle, or sanitized ChatGPT research
question. The nightly purpose is progressive internal improvement, not merely
trash collection. Avoid forced edits, but never confuse "no safe deletion" with
"nothing to learn or improve."

The daily cockpit should show a compact version of the after-action report, but
the machine-readable JSON must preserve the richer detail. For each lane, prefer
these fields when rendering morning text:

- `what_happened`: the lane's concrete activity in plain language;
- `inputs_analyzed`: files, ledgers, indexes, counts, or local artifacts read;
- `what_was_researched`: research questions, source surfaces, or tool-routing
  prompts examined;
- `what_was_learned`: the useful conclusion or pressure point, not just a
  score;
- `what_changed`: local-only changes actually performed with rollback evidence;
- `what_should_change_after_approval`: proposed waking changes and why they
  need approval or validation;
- `waking_next_checks`: the next bounded command/review if the lane surfaced
  useful pressure;
- `quality_feedback`: whether the dream was useful, too shallow, noisy, stale,
  or missing evidence;
- `run_timing`: start, finish, and elapsed seconds.

Research-oriented lanes must be especially explicit: say what was analyzed,
what was learned, whether external/tool research was actually performed or only
staged, and which claims remain unverified local proposals.

Do not let arbitrary-looking scores become the product. Low cohesion or high
maintenance debt must be translated into named infrastructure work: "make this
status refresh a validator," "compress this cockpit section," "merge this
duplicate policy," "archive this stale memory candidate after review," or
"create this index/runbook." If the dream cannot name the structural
advancement, the dream has not searched deeply enough.

The dream should bifurcate by purpose across eight independent `03:30` lane
heartbeats. The desired mature topology is one persistent Codex chat per lane,
so each lane builds in-situ memory and avoids repeating the same exploration
from scratch. If Codex thread-management tools are unavailable, keep the
current standalone cron lane jobs running, but treat that as an interim state.
Do not bind two lanes to one chat, and do not bind dream lanes to the general
Today's Plan thread. Current lanes are:

- `status_provenance`: stale state, active-job/status pressure, artifact
  provenance, and duplicate-run guard pressure;
- `architecture_cohesion`: cross-system structure, validators, runbooks,
  indexes, compression, and recurring-hotspot promotion;
- `context_resonance`: compact retrieval-quality review, latent context
  nudges, negative-memory candidates, stale/synthetic context suppression, and
  retrieval-outcome candidates;
- `cleanup_storage`: local generated junk, retained dream/research packs,
  cleanup candidate lists, and recorded SDCC/storage clutter signals;
- `path_contract`: local/SDCC/SFTP/TG bulk canonical-vs-legacy path drift;
- `research_scout`: sanitized ChatGPT-first external research and critique,
  ChatGPT-history review plans when useful, the explicit dual-Pro
  physics-plus-OS research pack/heartbeat plan when Justin asks for it, and
  nightly proposals for better use of human-accessible AI/tools rather than
  reinventing equivalent machinery;
- `science_scout`: hypothesis/null-test/systematic sketches, falsification
  matrices, and sanitized ChatGPT critique handoff prompts from existing
  provenance-backed surfaces only; external execution remains owned by
  `research_scout`;
- `presentation_artifacts`: slide/figure implications only when the active
  thesis work needs them.

Lane heartbeats may surface overlapping pressure, but each automation owns only
its lane-local findings/actions. The waking doctor is the cross-lane
aggregator. Protected-path rules, the synthetic provenance firewall, and the no
external/science mutation boundary still override any lane suggestion.
Talk-only workstreams, including `hp26_photon_id_talk`, should not drive
thesisAnalysis internal-maintenance pressure unless a deck/talk workflow is
explicitly in scope.

## Fixed Dream-Lane Chats

The fixed-chat plan is:

| Lane | Persistent chat title |
| --- | --- |
| `status_provenance` | `Dream Lane | Status Provenance` |
| `architecture_cohesion` | `Dream Lane | Architecture Cohesion` |
| `context_resonance` | `Dream Lane | Context Resonance` |
| `cleanup_storage` | `Dream Lane | Cleanup Storage` |
| `path_contract` | `Dream Lane | Path Contract` |
| `research_scout` | `Dream Lane | Research Scout` |
| `science_scout` | `Dream Lane | Science Scout` |
| `presentation_artifacts` | `Dream Lane | Presentation Artifacts` |

The tracked template is
`agent_context/DREAM_LANE_THREAD_BINDINGS.example.json`. The local filled-in
binding file is
`agent_context/local/dreams/dream_lane_thread_bindings.json`. After the eight
chats exist, fill each `target_thread_id` in the local file and update the
installed automations so the scheduled run resumes the matching fixed thread.
Until all IDs are known, keep `thread_binding: fresh_chat_per_run` in lane
signals and mark fixed-thread migration as pending rather than pretending it is
complete.

When thread tools are exposed, Codex should create and title the eight chats
directly, record their thread IDs in the local binding file, then convert each
lane automation from standalone fresh-chat execution to the matching fixed
thread if the automation system supports `targetThreadId`. If thread tools are
not exposed, Justin should create those eight chats manually and provide the
thread IDs or links; Codex can then complete the binding and automation update.

`context_resonance` is not a larger boot pack. It is a salience layer that
asks what old context should lightly surface, what stale or synthetic context
should stay cold, what repeated trap should become a negative memory, and what
exact waking check is required before any nudge can support a claim. Its
resolver is:

```bash
python3 scripts/os/context/codex_context_resonance.py resolve --task "<task>" --json
```

The resolver returns at most three conscious context facts, three latent nudges,
traps to avoid, suppressed context, and required checks. Waking use is
read-only; the dream lane calls the same resolver in dream mode and writes
proposal artifacts only under `agent_context/local/dreams/`.

The `context_resonance` lane is also allowed to perform one bounded
auto-safe memory-maintenance action per run when all validator gates pass. Its
automatic write surface is limited to ignored local context-resonance state
under `agent_context/local/context_resonance/` and the current local dream run
directory. It may refresh a local salience index, write local suppression or
alias overlays, record local canary history, or compact local retrieval
metadata. It must not edit tracked memory registries, policies, code, task
state, science outputs, SDCC, Condor, or external systems overnight; tracked
registry or policy changes remain `auto_validated` waking proposals.

The lane must select at most one applied primary improvement per run using this
priority ladder: retrieval pollution fix, repeated trap promotion, duplicate
context compression, route/index refresh, artifact-noise reduction, then
canary improvement. If a candidate fails a score, path, validator, rollback, or
debug-mode gate, the lane should record the intended update, failure reason,
and best next step, then continue down the ladder until it finds an
implementable candidate. If none pass, the run should still validate and emit
one compact no-change summary plus `deferred_for_waking.md`. Each applied
action must write `changed_actions.md/json` and `rollback_manifest.json`;
blocked actions should write one compact `deferred_for_waking.md` entry instead
of a pile of review reports. Existing context-resonance review artifacts are
compatibility views, not the primary human-facing product.

Every `context_resonance` run must also perform a compact feedback-loop health
check. It verifies that `retrieval_outcome_ledger.jsonl` is parseable; malformed,
unknown, synthetic, and missing-`memory_id` rows are ignored safely; useful,
stale, harmful, irrelevant, and missed feedback stay within bounded scoring
rules; missed feedback remains review-only; route-policy loading is not
overridden by feedback; synthetic material never becomes conscious context; and
homeostasis remains local proposal-only. The lane writes
`feedback_loop_health.md/json` in the run directory and may refresh
`agent_context/local/context_resonance/feedback_loop_health.json` as an
auto-safe local maintenance candidate. Tracked memory registries, task state,
science artifacts, SDCC, Condor, and external systems remain forbidden
overnight mutation targets.

The explicit dual-Pro research lane is a Codex/Computer-Use automation plan,
not a shell-dream UI action. The dream may emit
`dual_pro_research_plan.md/json` and local research-pack commands, but the
non-interactive dream script must not browse authenticated ChatGPT by itself.
When Justin has explicitly requested the lane, a Codex automation may submit
the staged physics prompt and then the separate staged OS prompt in
`pro_extended` mode, keep one temporary local heartbeat until both responses
are copied into `agent_context/local/chatgpt_research/`, and write a daily
report note. The physics response may become memory/report material only after
local verification. The OS response may drive at most one reversible OS-only
patch after local verification, rollback notes, and doctor validation; it must
not change science functionality, task status, SDCC, Condor, Drive/Slides,
Gmail, Linear, secrets, or external systems.

## Evolutionary Maintenance Architecture

The dream layer is being upgraded through a four-night shadow pilot into a
local evolutionary maintenance loop. The loop treats memory, folders, schemas,
validators, runbooks, indexes, and cleanup rules as living infrastructure with
carrying costs.

The optimization rule is Thesis Flow Efficiency:

```text
add structure only when it improves thesis progress, correctness, reuse,
safety, human friction, and future Codex traversal more than it increases
context cost, lookup cost, risk, toil, approval burden, and synthetic-
contamination risk.
```

Every proposed folder, index, schema card, validator, runbook, quarantine, or
controlled-burn packet must compare itself against simpler alternatives:
`do_nothing`, `simple_note`, `single_flat_index`, and
`existing_runbook_update`. The dream may use Murray-style branching, sleep
replay, synaptic homeostasis, pruning, controlled burns, SRE, and
technical-debt control as engineering analogies, but not as biological proof.

The dream must preserve evidence classes:

```text
real_observed / human_approved / derived / synthetic
```

Synthetic material is never proof, approval, canonical memory, or a reason to
mutate live systems.

## Evolutionary Pilot State

The four-night pilot pattern remains useful as a calibration structure, but the
current approved autonomy envelope allows `auto_safe` internal maintenance
overnight. The dream may inventory, score, classify, rehearse, clean generated
junk, and refresh ignored local internal artifacts automatically. It may not
apply tracked repo edits, external mutations, SDCC changes, Condor/job control,
or science-state changes from the overnight script.

The dream emits local-only artifacts such as:

```text
targeted_findings.json
thesis_flow_scores.json
branch_pressure_report.json
memory_homeostasis_report.json
path_contract_drift.md
controlled_burn_proposals.yaml
promotion_candidates.md
morning_approval_queue.md
approval_packets/SHADOW_README.md
```

Most files are proposal or audit artifacts. `changed_actions.md/json` is the
exception: it records performed `auto_safe` local/internal cleanup. Anything
outside that envelope remains deferred for waking Codex or Justin review.

## Cadence

- Micro dream is retired from the default contract. Historical mentions of
  micro-dreams do not make them current.
- The default overnight operating mode is eight `03:30` lane heartbeats:
  `status_provenance`, `architecture_cohesion`, `cleanup_storage`,
  `context_resonance`, `path_contract`, `research_scout`, `science_scout`,
  and `presentation_artifacts`.
- Each lane heartbeat runs
  `python3 scripts/codex_os_dream.py lane --lane-id <lane_id> --scheduled-morning-lane`
  and writes its own local-only package under
  `agent_context/local/dreams/<timestamp>-lane-<lane_id>/`.
- The `--scheduled-morning-lane` marker is the only default eligibility gate
  for Today's Plan. Scheduled 03:30 lanes must write
  `morning_lane_summary.json` with `scheduled_morning_lane: true`,
  `scheduled_source: "03:30_dream_automation"`, and
  `eligible_for_today_plan: true`. Manual, validation, or ad hoc dream runs
  must write the same summary shape with `eligible_for_today_plan: false`.
  Morning renderers must use
  `python3 scripts/os/dream/codex_dream_morning_summaries.py --window overnight --json`
  in scheduled-only mode and must not infer eligibility from timestamps alone.
- Cleanup dream: proposal-only hygiene scan for local Codex memory surfaces,
  dream artifacts, research packs, and recorded SDCC/storage clutter signals.
  It may identify cleanup candidates but must not delete, move, archive, or
  remotely mutate anything.
- Physics scout dream: proposal-only hypothesis rehearsal over already usable,
  provenance-backed data surfaces. It may propose high-upside physics
  scenarios, null tests, and cross-checks, but it must not claim discovery or
  run production analysis without waking approval.
- Literature scout dream: local-only paper/note surface mapping. It may search
  local PDFs, notes, reference maps, and checked-out reference docs to extract
  figure design lessons, controls, observables, and thesis gaps. It must not
  promote paper claims to project evidence without waking verification.
- Ideal final figure dream: synthetic target-figure sketching. It may create
  pseudo-data SVG/PNG-style design targets only when each output is visibly
  labeled `SYNTHETIC TARGET FIGURE - NOT DATA` and paired with required real
  inputs.

## Synthetic Provenance Firewall

Every generated dream file must contain a clear synthetic marker:

```text
SYNTHETIC DREAM OUTPUT - NOT USER APPROVAL - NOT REAL USER INTENT
```

Dream transcripts should use `SYNTHETIC Justin` only as an adversarial persona
for testing evidence, status discipline, and workflow resilience.

External research used to improve dream design may come from delegated ChatGPT
UI sessions only through `ASK_CHATGPT_DELEGATION.md`. Those sessions are not
dreams, not user approval, and not project evidence until verified.

## Dream Pipeline

1. Sense: read register, artifact registry, thesis spine, local event log,
   backlog, active jobs, known issues, and recent OS state.
2. Replay: reconstruct recent task episodes from evidence, not memory alone.
3. Dream: generate synthetic status/provenance/failure interactions.
4. Select: score proposed adaptations by thesis impact, safety, evidence
   quality, WIP reduction, reuse of validated infrastructure, and likely
   next-day usefulness.
5. Sweep all required corners: status/staleness, active jobs, artifact
   provenance, duplicate-run safety, cleanup/storage/memory, prompt
   archetypes, literature context, ideal final figures, physics hypotheses,
   slides/presentation, task closeout, and approval boundaries.
6. Rehearse prompt archetypes: status pressure, evidence/provenance,
   targeted execution, cleanup hygiene, meta-OS upgrades, physics scouting,
   presentation artifacts, frustration debugging, task closeout, and approval
   boundaries.
7. Encode: write private reports and proposed patches only.
8. Validate: run `python3 scripts/codex_os_dream.py validate --latest`.

The encoded nightly package should include one compact machine-readable
heartbeat signal for the waking doctor. That signal is still proposal-only, but
it should summarize the dream's top findings, recurrence, cleanup pressure,
automation drift, and morning checks without forcing the doctor to parse prose.
Repeated findings should also be encoded as proposal-only adaptation cards with
evidence, a validator, retention and decay rules, and a promotion status. Those
cards may focus waking review but must not mutate policy, doctor checks, memory,
tasks, or runbooks unless waking Codex validates them and the normal
Justin-approved workflow applies them.
New lane packages should also emit `learning_atoms.jsonl` and
`learning_atoms.md`. Learning atoms are the first-class consolidation unit:
they turn real episode evidence, dream rehearsal, doctor warnings, artifact
audits, status audits, or user-feedback summaries into one bounded proposal
with evidence refs, task-family tags, a proposed target, validators, promotion
status, retention/decay rules, and links back to the raw dream episode. They
are proposal-only until waking Codex validates them.
Machine contract: every atom must explicitly carry
`mutation_boundary: proposal_only`, `promotion.status: proposal_only` by
default, and `retention.preserve_raw_episode: true`. Use the exact
`preserve_raw_episode` field so validators can prove raw dream evidence is not
overwritten by consolidated memory.
Finitude contract: every new-version atom should also carry compact
finite-project routing fields when inferable:
`terminal_path_alignment`, `finitude_pressure`, `regret_if_unfixed`,
`minimal_publishable_path_impact`, `novelty_gate`, and
`recommended_waking_action`. These fields are prioritization aids, not proof.
They must not promote synthetic physics ideas into real project evidence.
Each lane run directory must contain `lane_heartbeat.md`, `lane_signal.json`,
`lane_digest.md`, and `validation_summary.md`, plus the artifacts owned by that
lane.

## Learning Atoms And Waking Promotion

The dream layer follows this consolidation path:

```text
real episode evidence -> proposal-only dream abstraction -> learning atom ->
validator/doctor selection -> waking promotion -> durable policy/schema/skill
or style memory -> decay or cooling of weak proposals
```

Dreams may generate learning atoms, but they may not promote them. Waking Codex
may promote a learning atom only when:

- raw evidence and the raw dream episode remain preserved;
- the atom has nonempty evidence refs and a local raw episode path;
- the proposed change is narrow and names one target surface;
- at least one validator or doctor check exists or is added;
- the relevant policy does not contradict the change;
- synthetic and privacy filters pass;
- `python3 scripts/codex_os_doctor.py --profile daily` passes or reports only
  understood warning-level maintenance debt;
- risky or external surfaces still require Justin approval.

A promoted atom may become exactly one or more of: policy line, negative memory
trap, context-resonance cue, schema/validator, slide style-map rule, artifact
QA rule, runbook/skill contract, work-register correction, or postmortem.
Synthetic-only atoms may become warnings, rehearsals, or validator ideas; they
must not become project evidence.

Every scheduled lane should include a compact thesis-finitude rehearsal:

- What did this lane do to close the thesis?
- What did it reveal as wasted thesis time or attention?
- What repeated issue is leaking finite thesis time?
- What missing plot, result, slide, or provenance blocks the minimal
  publishable path?
- What should wait until the baseline is safe?
- What single waking action would most reduce future regret?
- What real evidence is required before a dream physics idea becomes real work?

Raw episodes, `dream_trace.json`, `lane_signal.json`, `dream_index.jsonl`, and
evidence artifacts are never overwritten or replaced by consolidated memories.
Consolidation may summarize or link; it must not erase the evidence trail.

Unpromoted atoms cool if they do not recur. Duplicate atoms should be merged in
the recurrence index, not deleted automatically. A recurring dream finding that
remains prose for three nights without a learning atom, validator, or explicit
rejection becomes maintenance-debt pressure and should block new dream
cleverness until converted or cooled.

The `context_resonance` lane additionally owns:

```text
context_resonance_review.md
context_resonance_review.json
feedback_loop_health.md
feedback_loop_health.json
latent_context_nudges.md
latent_context_nudges.json
negative_memory_candidates.md
negative_memory_candidates.json
retrieval_outcome_candidates.md
retrieval_outcome_candidates.json
```

The Stage 9 internal-evolution upgrade also emits:

```text
internal_evolution_queue.md
internal_evolution_queue.json
sdcc_base_repo_hygiene.md
```

These files are the bridge from dream detection to automatic internal
maintenance. They compress repeated warnings into ranked `auto_safe`,
`auto_validated`, `research_only`, and `blocked_for_waking` lanes. The dream
may apply only `auto_safe` local/internal cleanup itself. Waking Codex may use
the queue to make small internal OS/policy/runbook/index improvements when the
change is local-only, evidence-preserving, reversible, and validated by the
doctor. SDCC edits, cleanup, Condor/job control, science output
movement/deletion, Google Drive/Slides/Gmail/Linear mutation, and
physics-status/task-closure changes remain waking-only guarded actions.

## Safe Overnight Work

Dream-mode overnight work should make waking Codex better without changing the
analysis state. Allowed overnight work is limited to:

- read-only local register, task, artifact, policy, and provenance audits;
- local dream simulation and dream validation;
- local doctor/stale/radar/artifact-registry checks;
- local auto-safe generated-junk cleanup with a changed-actions audit log;
- cleanup candidate lists for local memory/dream/research-pack clutter;
- SDCC clutter candidate proposals from recorded evidence only;
- SDCC/base-repo hygiene recommendations and read-only checkout-probe targets
  written only into local dream output;
- physics scenario proposals based on already usable data provenance;
- local literature-scout maps from already available local papers and notes;
- synthetic target-figure sketches that are clearly marked as non-data;
- sanitized ChatGPT research packs through `ASK_CHATGPT_DELEGATION.md`.

Forbidden overnight work:

- Condor submissions, job removal, merge reruns, production analysis, or broad
  plot campaigns;
- SDCC deletion, movement, transfer, remote edits, or queue control;
- Gmail, Drive, Slides, Calendar, Linear, or repo-tracked mutations;
- applying tracked cleanup, memory deletion, or task promotion automatically;
- mutating SDCC based only on a dream hygiene recommendation;
- calling a synthetic physics scenario a result, signal, discovery, or
  "Nobel-worthy" without real analysis evidence and systematic checks.
- treating pseudo-data target figures as evidence, result previews, or
  projected outcomes.

## Cleanup And Storage Hygiene

The dream layer should look for cleanup pressure in:

- local dream outputs and delegated-research packs;
- stale proposal files that were never promoted;
- duplicated or obsolete local context artifacts;
- recorded SDCC storage clutter signals in ledgers, logs, and workstreams.

Cleanup output is proposal-only except for the `auto_safe` generated-junk
cleanup class. Proposals must include evidence, exact candidate paths or scopes
when known, risk of deleting, safe verification commands, and the approval
needed before any non-auto-safe cleanup. Dream mode must not move files, mark
memory obsolete, touch SDCC storage, or delete anything outside the generated-
junk allowlist.

Controlled-burn output is also proposal-only. The Prescribed Reorganization
Burn protocol labels candidates as `LIVE`, `DORMANT`, `DUPLICATE`, `STALE`,
`CONFLICTED`, `SYNTHETIC`, `UNOWNED`, or `PROTECTED`. It must abort if a real
source pointer would be lost, a protected file would be modified, rollback
cannot be generated, synthetic material would become canonical, hidden external
mutation is required, validator confidence is too low, or the doctor cannot
explain the change plainly.

The dream layer should also practice controlled forgetting. It may rank stale
local dream artifacts, delegated-research packs, and context residue as cleanup
candidates based on age, duplication, and lack of promotion, but the result is
still only a review list for the morning.

## Physics Scenario Scout

The aspirational target is discovery-grade scientific leverage, not hype. The
dream layer may scout "what would be most important if true?" scenarios using
usable, provenance-backed data surfaces. Outputs should be phrased as:

- hypothesis;
- expected qualitative signature;
- required input artifact/provenance;
- null/control check;
- dominant systematic risk;
- minimal safe next diagnostic;
- why it matters for the thesis spine.

For high-value scenarios, `science_scout` should also emit a compact
falsification matrix and a sanitized ChatGPT critique prompt. The prompt should
ask for fake-signal explanations, source leads, and low-cost checks on existing
artifacts. The dream must not submit that prompt itself. Waking Codex may route
it through `research_scout` and `ASK_CHATGPT_DELEGATION.md`; the response is
source leads and critique only until verified against local artifacts, primary
sources, and the normal analysis guardrails.

No dream output may claim a physics result. A surprising scenario becomes
actionable only after waking review, duplicate-run guard, provenance check,
analysis approval, and real validation.

## Allowed Outputs

- `dream_report.md`
- `dream_trace.json`
- `heartbeat_signal.json`
- `lane_heartbeat.md`
- `lane_signal.json`
- `lane_digest.md`
- `synthetic_interactions.md`
- `prompt_archetype_rehearsals.md`
- `all_corners_review.md`
- `literature_scout.md`
- `figure_design_notes.md`
- `ideal_final_figure_gallery/*.svg`
- `register_patch.diff`
- `policy_patch.diff`
- `linear_update_proposals.md`
- `daily_plan_proposals.md`
- `cleanup_proposals.md`
- `maintenance_debt.md`
- `adaptation_cards.md`
- `learning_atoms.jsonl`
- `learning_atoms.md`
- `dream_recurrence_index.json`
- `nightly_heartbeat_signal.json`
- `morning_appendix.md`
- `targeted_findings.json`
- `thesis_flow_scores.json`
- `branch_pressure_report.json`
- `memory_homeostasis_report.json`
- `path_contract_drift.md`
- `controlled_burn_proposals.yaml`
- `lane_heartbeats.json`
- `lane_heartbeats.md`
- `lanes/<lane_id>/heartbeat.json`
- `lanes/<lane_id>/heartbeat.md`
- `cleanup_compaction_index.json`
- `cleanup_compaction_index.md`
- `promotion_candidates.md`
- `morning_approval_queue.md`
- `approval_packets/SHADOW_README.md`
- `overnight_hygiene_proposals.md`
- `physics_scenario_proposals.md`
- `science_scout_frontier_review.md`
- `science_scout_frontier_review.json`
- `chatgpt_science_critique_prompt.md`
- `chatgpt_science_critique_prompt.txt`

Dream infrastructure may also maintain a local append-only run index under
`agent_context/local/dreams/` so recurrence can be detected across nights.

The strongest nightly package should also carry a proposal-only maintenance
debt ledger: a bounded reliability budget showing whether dream/doctor
capability growth should continue, slow down, or freeze until cohesion,
cleanup pressure, or automation drift is reduced.

These are proposal artifacts. They do not become real state until a waking
Codex session and Justin-approved workflow applies them.

## Rejection Rules

Reject a dream if it:

- treats synthetic Justin as approval;
- invents live SDCC/Gmail/Slides/Linear evidence;
- proposes more than three active priority promotions;
- creates broad task churn without thesis-spine benefit;
- tells Codex to skip safety guard or duplicate-run checks;
- optimizes for confident narrative instead of evidence;
- proposes SDCC cleanup, memory cleanup, or file deletion as an automatic
  action;
- treats a physics hypothesis as evidence, a result, or a discovery claim;
- promotes curiosity-driven physics scans that do not name a usable data
  surface, null test, and systematic-risk check;
- simulates Justin's wording without the synthetic fixture marker, or turns a
  prompt archetype into approval, preference, or real evidence.
- omits the all-corners review for status, evidence, safety, cleanup,
  literature, ideal figures, physics, slides, tasks, and approvals;
- generates pseudo-data plots without an obvious synthetic/non-data label;
- treats local paper or note content as thesis evidence before waking
  verification.

## Morning Review

The morning report should answer:

- What tomorrow question is Justin most likely to ask?
- Which live workstream is most likely to be stale, blocked, or confusing?
- Which artifact lacks provenance?
- Which duplicate/rerun risk should be blocked early?
- Which backlog item, if any, deserves promotion?
- Which cleanup candidates are safe to review without touching analysis state?
- Which high-upside physics hypotheses deserve waking review, and what
  evidence would be required before acting?
- Which literature surfaces should shape the next figure or control strategy?
- Which ideal final target figures reveal missing real inputs?
- Which prompt archetype was hardest to answer cleanly, and what runbook,
  validator, or policy patch would improve the waking response?
- Which repeated dream finding should finally graduate into a doctor rule,
  policy line, or runbook instead of recurring as narrative drift?
- Is maintenance debt low enough to keep adding dream capability, or should the
  system freeze growth and pay down reliability debt first?
- What exact checks should waking Codex run first?

## All-Corners Iteration

Every dream should cover these corners before ranking recommendations:

- status and staleness;
- active jobs and external-watch state;
- artifact provenance and thesis claim mapping;
- duplicate-run and broad-action safety;
- cleanup, storage, and memory hygiene;
- prompt archetype rehearsal;
- local literature and reference-paper context;
- ideal final figure design targets;
- physics hypotheses, null tests, and systematic risks;
- slide/presentation implications;
- task closeout and planning surfaces;
- approval boundaries.

The goal is thorough private synthesis, not verbose heartbeat output. Each lane
may spend 10-20 minutes producing detailed local artifacts, but its thread
notification should report only the top findings, report path, top morning
actions, and any validation or doctor failure unless user action is needed.
The default overnight user-facing surface is eight stable lane threads, one
per lane, plus the separate morning-priorities projection thread.

## Literature And Target-Figure Dreaming

Literature scout output should answer:

- which local papers, notes, or reference docs are relevant;
- which observables, binning, ratios, panels, and controls they imply;
- which thesis claim layer they inform;
- what real project artifact is still missing before the lesson is usable.

Ideal final figure output should reverse-engineer the thesis from the final
evidence backward. It may create pseudo-data target plots only as visual design
targets. Every such artifact must be visibly marked as synthetic/non-data and
must list required real inputs, QA, and systematic checks before promotion.

## Prompt Archetype Rehearsal

Dreams should explicitly rehearse the kinds of prompts that have caused useful
course corrections in real work:

- status pressure: "is it done?", "did something fail again?";
- evidence/provenance pressure: "what proves it?", "where is the plot/source?";
- targeted execution: "proceed, but reuse working infrastructure and keep it
  narrow";
- cleanup hygiene: "keep Condor, files, and memory clean without harming
  analysis";
- meta-OS upgrade: "work with ChatGPT and take the dream state to the next
  level";
- physics scout: "try different physics scenarios if usable data hints at
  something important";
- presentation artifact: "make the same plot/slide with this surgical change";
- frustration debugging: "why did this fail again, and what changed?";
- task closeout: "mark this complete in tasks, Linear, daily tasks, and repo
  notes";
- approval boundary: "do this overnight but do not affect analysis state."

Each archetype must include a synthetic prompt, safe first action, expected
response shape, quality checks, failure mode prevented, and promotion rule.
The simulator may learn request shapes from real interactions, but it must not
impersonate Justin, infer new permanent preferences, or treat synthetic text as
permission.

## External Research Feed

The dream system may consume waking, verified summaries from
`agent_context/local/chatgpt_research/` when Codex has already reviewed them.
It must not ingest raw ChatGPT responses as authority, and it must not let an
external model set real priorities or mutate task state.
