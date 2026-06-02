# Agentic OS Dreaming

This policy defines the ThesisAnalysis dream layer: private, synthetic,
proposal-only simulation that helps Codex rehearse likely future failures and
prepare better next-day responses without touching live external systems.

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
Scores are diagnostic smoke alarms only. `morning_conversation_digest.md` is
the preferred user-facing summary: it must say what the dream did, learned,
tried, changed, deferred, why it deferred, which structural advancements it
found, and the best next internal fix.

No dream may stop at "nothing to clean." If generated-junk cleanup is empty,
the dream must still produce at least one reviewable cleanup angle, architecture
improvement hypothesis, deeper search angle, or sanitized ChatGPT research
question. The nightly purpose is progressive internal improvement, not merely
trash collection. Avoid forced edits, but never confuse "no safe deletion" with
"nothing to learn or improve."

Do not let arbitrary-looking scores become the product. Low cohesion or high
maintenance debt must be translated into named infrastructure work: "make this
status refresh a validator," "compress this cockpit section," "merge this
duplicate policy," "archive this stale memory candidate after review," or
"create this index/runbook." If the dream cannot name the structural
advancement, the dream has not searched deeply enough.

The dream should bifurcate by purpose across eight independent `03:30` lane
heartbeats, each opening a fresh automation-generated chat each night and each
writing its own local dream run directory. Current lanes are:

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
  ChatGPT-history review plans when useful, and nightly proposals for better
  use of human-accessible AI/tools rather than reinventing equivalent
  machinery;
- `science_scout`: hypothesis/null-test/systematic sketches from existing
  provenance-backed surfaces only;
- `presentation_artifacts`: slide/figure implications only when the active
  thesis work needs them.

Lane heartbeats may surface overlapping pressure, but each automation owns only
its lane-local findings/actions. The waking doctor is the cross-lane
aggregator. Protected-path rules, the synthetic provenance firewall, and the no
external/science mutation boundary still override any lane suggestion.
Talk-only workstreams, including `hp26_photon_id_talk`, should not drive
thesisAnalysis internal-maintenance pressure unless a deck/talk workflow is
explicitly in scope.

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
- Each lane heartbeat runs `python3 scripts/codex_os_dream.py lane --lane-id <lane_id>`
  and writes its own local-only package under
  `agent_context/local/dreams/<timestamp>-lane-<lane_id>/`.
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
Each lane run directory must contain `lane_heartbeat.md`, `lane_signal.json`,
`lane_digest.md`, and `validation_summary.md`, plus the artifacts owned by that
lane.

The `context_resonance` lane additionally owns:

```text
context_resonance_review.md
context_resonance_review.json
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
