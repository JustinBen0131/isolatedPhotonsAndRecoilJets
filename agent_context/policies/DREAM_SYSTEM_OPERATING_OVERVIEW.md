# ThesisAnalysis Dream System Operating Overview

This document is the top-down handoff for a blank Codex model. It explains
what the dream system is, which parts are automated, what each lane owns, what
Justin should and should not have to read, and how waking Codex should turn
nightly proposals into real improvements without corrupting project state.

For detailed rules, load `agent_context/policies/AGENTIC_OS_DREAMING.md`.
This file is the map; the dreaming policy is the law.

## One-Sentence Model

The dream system is a private, scheduled, proposal-only nightly rehearsal layer:
eight persistent Codex chats run at 03:30, each inspects one slice of the
ThesisAnalysis operating system, writes local synthetic artifacts under
`agent_context/local/dreams/`, and gives waking Codex compact evidence-gated
proposals. Dreams do not approve work, prove science, mutate external systems,
or replace Justin.

## The Three Layers

### 1. Night dream lanes

Eight heartbeat automations run every day at 03:30 local time. Each automation
is attached to one persistent Codex chat named `Dream Lane | ...`. The lane
prompt reads `AGENTS.md` and `AGENTIC_OS_DREAMING.md`, then runs:

```bash
python3 scripts/codex_os_dream.py lane --lane-id <lane_id> --scheduled-morning-lane
```

Each lane writes a local run package under:

```text
agent_context/local/dreams/<timestamp>-lane-<lane_id>/
```

Scheduled lane runs must mark themselves as eligible for the morning digest:

```text
scheduled_morning_lane: true
scheduled_source: "03:30_dream_automation"
eligible_for_today_plan: true
```

Manual, validation, and ad hoc runs must not appear in Today's Plan unless a
debug command explicitly asks for manual inclusion.

### 2. Morning projection

The daily planning automation runs at 10:00 local time in the main morning
thread. It uses:

```bash
python3 scripts/os/dream/codex_dream_morning_summaries.py --window overnight --json
```

Then the register renderer inserts `Overnight Dream Updates` into Today's Plan
after `Do First Today`. This section should be short, concrete, and
decision-oriented. It is not a raw dream log.

Justin should normally read only Today's Plan, not every dream artifact.

### 3. Waking promotion

Waking Codex, not the overnight dream, decides whether a proposal becomes a
real repo policy, memory trap, validator, task correction, runbook, style-map
rule, or postmortem. Promotion requires:

- real or preserved raw evidence;
- nonempty evidence refs;
- a narrow target surface;
- a validator or doctor check;
- privacy and synthetic filters passing;
- `python3 scripts/codex_os_doctor.py --profile daily` passing or producing
  only understood warning-level maintenance debt;
- explicit Justin approval for risky or external surfaces.

Dreams can say "this looks important." Waking Codex must prove "this is safe
and worth changing."

## Core Safety Boundaries

Dream output is synthetic. It is never:

- Justin approval;
- a science result;
- evidence that a plot is valid;
- evidence that a task is done;
- permission to submit jobs;
- permission to mutate SDCC, Condor, Drive, Slides, Gmail, Linear, Calendar,
  GitHub, or live analysis outputs.

Allowed overnight work is local and internal:

- read-only inspection of repo ledgers, register, artifact registry, policies,
  local notes, and local dream/research history;
- proposal-only failure rehearsal;
- local ignored index refreshes;
- generated-junk cleanup only when explicitly allowed by policy and logged in
  `changed_actions.md/json`;
- local dream package generation and validation.

Forbidden overnight work includes:

- Condor submissions or job control;
- SDCC edits, transfers, deletes, movement, or queue control;
- Google Slides, Drive, Gmail, Linear, Calendar, GitHub mutation;
- science output creation, deletion, or claim promotion;
- task-status closure or priority changes based only on dream output.

## Live Automations

These are the intended live dream and planning automations. The dream lanes are
heartbeat automations bound to persistent lane chats. The morning automation is
the daily user-facing projection.

| Automation id | Chat title | Time | Lane id | Purpose |
| --- | --- | --- | --- | --- |
| `thesisanalysis-dream-status-provenance` | `Dream Lane | Status Provenance` | 03:30 | `status_provenance` | Stale state, active jobs, artifact provenance, duplicate-run pressure. |
| `thesisanalysis-dream-architecture-cohesion` | `Dream Lane | Architecture Cohesion` | 03:30 | `architecture_cohesion` | Cross-system structure, validators, runbooks, indexes, compression, recurring-hotspot promotion. |
| `thesisanalysis-dream-lane-context-resonance` | `Dream Lane | Context Resonance` | 03:30 | `context_resonance` | Retrieval quality, latent nudges, negative-memory candidates, stale/synthetic suppression, feedback health. |
| `thesisanalysis-dream-cleanup-storage` | `Dream Lane | Cleanup Storage` | 03:30 | `cleanup_storage` | Generated junk, retained dream/research packs, cleanup candidates, storage clutter signals. |
| `thesisanalysis-dream-path-contract` | `Dream Lane | Path Contract` | 03:30 | `path_contract` | Local/SDCC/SFTP/TG bulk canonical-vs-legacy path drift. |
| `thesisanalysis-dream-research-scout` | `Dream Lane | Research Scout` | 03:30 | `research_scout` | Sanitized external-research plans, source-lead critique, ChatGPT/tool leverage proposals. |
| `thesisanalysis-dream-science-scout` | `Dream Lane | Science Scout` | 03:30 | `science_scout` | Physics hypotheses, null tests, systematic sketches, final-figure scaffolds from existing evidence surfaces. |
| `thesisanalysis-dream-presentation-artifacts` | `Dream Lane | Presentation Artifacts` | 03:30 | `presentation_artifacts` | Slide/figure implications, presentation artifact hygiene, style/regression pressure without Slides mutation. |
| `thesisanalysis-morning-priorities` | main morning priorities thread | 10:00 | none | Generates the daily user-facing Today's Plan from register, dream summaries, and active project surfaces. |

The local fixed-thread binding file is:

```text
agent_context/local/dreams/dream_lane_thread_bindings.json
```

The tracked example/template is:

```text
agent_context/DREAM_LANE_THREAD_BINDINGS.example.json
```

If a lane opens in the wrong chat, first check the automation `target_thread_id`
against the local binding file.

## What Each Lane Owns

### `status_provenance`

Goal: keep status claims honest.

This lane looks for stale active/running/review workstreams, active jobs without
fresh evidence, artifact readiness claims without registry support, and
duplicate-run risk. It should produce waking checks like "verify this job
state," "refresh this stale workstream," or "turn this repeated status gap into
a validator."

It must not mark tasks done, close workstreams, submit jobs, remove jobs, or
change production state.

### `architecture_cohesion`

Goal: keep the OS coherent without making it larger for its own sake.

This lane looks for repeated policy duplication, missing runbooks, weak
validators, context bloat, overgrown dream machinery, and places where a small
schema or check would prevent repeated work. It should prefer one small
validated improvement over another new subsystem.

It must not treat OS work as thesis progress unless it reduces future thesis
waste or protects evidence quality.

### `context_resonance`

Goal: improve what future Codex sessions retrieve, suppress, and remember.

This lane runs or packages `scripts/os/context/codex_context_resonance.py`.
It may refresh local ignored salience indexes under
`agent_context/local/context_resonance/`, inspect retrieval-outcome feedback,
surface negative-memory candidates, and propose routing fixes.

It must not edit tracked memory registries overnight. Tracked changes to
`CONTEXT_RESONANCE_INDEX.yaml`, `NEGATIVE_MEMORY_MAP.yaml`, or
`SCHEMA_REGISTRY.yaml` are waking-Codex changes only.

### `cleanup_storage`

Goal: reduce local clutter without risking evidence.

This lane reviews dream directories, delegated-research packs, generated junk,
local context artifacts, and recorded storage pressure. It may clean only
policy-allowed generated junk and must log any performed cleanup. Most cleanup
should remain a proposal.

It must not delete science outputs, remote files, SDCC files, raw evidence,
repo-tracked files, or anything needed for provenance.

### `path_contract`

Goal: keep local and remote path expectations sane.

This lane inspects local path conventions, command aliases, transfer maps,
canonical-vs-compat script locations, SDCC path-contract notes, and drift
between old and new organization. It should produce bounded path-contract
warnings and migration proposals.

It must not edit SDCC, upload files, move protected scripts, or change Fun4All,
Condor, transfer, or base pipeline entrypoints without waking approval.

### `research_scout`

Goal: make external research and ChatGPT-style delegation useful without
letting it become authority.

This lane creates sanitized research questions, identifies local reference
surfaces, proposes ChatGPT/deep-research prompts, and extracts source-lead
critique. It can suggest what to ask external tools, but the noninteractive
dream itself should not browse authenticated UI or treat external output as
truth.

It must not promote external model text into project evidence. Waking Codex
must verify sources locally before changing policies, tasks, plots, or science
claims.

### `science_scout`

Goal: search for high-value physics checks without hallucinating progress.

This lane sketches hypotheses, null tests, systematics, and target figures
using existing provenance-backed project surfaces. It may identify what would
be most thesis-closing if verified.

It must not call a hypothesis a result, start production, or add active science
work just because the idea is interesting. The Thesis Finitude & Care Kernel
keeps novelty behind minimal publishable baseline safety unless the check is
bounded and directly protects the baseline.

### `presentation_artifacts`

Goal: prevent repeated slide/plot presentation regressions.

This lane watches for slide-readiness gaps, missing speaker scripts, tiny text,
internal notes on the slide canvas, source-plot confusion, provenance footers,
plot readability, and reusable slide-story scaffolds. It should generate
proposal-only style-memory or slide-audit candidates.

It must not mutate Google Slides or claim slide readiness without waking
inspection and artifact evidence.

## Artifact Contract Per Lane Run

A modern lane package should include at least:

```text
lane_heartbeat.md
lane_signal.json
lane_digest.md
validation_summary.md
morning_lane_summary.txt
morning_lane_summary.json
morning_appendix.md
learning_atoms.jsonl
learning_atoms.md
dream_recurrence_index.json
nightly_heartbeat_signal.json
promotion_candidates.md
maintenance_debt.md
```

Some lanes also write specialized artifacts, for example:

```text
context_resonance_review.json
feedback_loop_health.json
latent_context_nudges.json
negative_memory_candidates.json
retrieval_outcome_candidates.json
changed_actions.json
rollback_manifest.json
structural_advancements.json
dual_pro_research_plan.json
human_tool_leverage.md
ideal_final_figure_gallery/
```

Do not assume every legacy run has every modern file. Validators are version
aware; older packages can remain valid if they predate new artifacts.

## Learning Atoms

Learning atoms are the bridge from night pressure to waking improvement.

Path:

```text
learning_atoms.jsonl
learning_atoms.md
```

Each atom is one bounded proposal. It should include:

- stable id and dream run id;
- source kind and evidence refs;
- task family;
- lesson type;
- problem symptom and recurrence;
- proposed change target;
- validator commands;
- promotion status, normally `proposal_only`;
- retention and decay rules;
- raw episode link;
- finitude fields such as `terminal_path_alignment`, `finitude_pressure`,
  `regret_if_unfixed`, `minimal_publishable_path_impact`, `novelty_gate`, and
  `recommended_waking_action`.

Waking Codex should not promote an atom because it sounds smart. It should ask:

1. Is there real evidence or only synthetic pressure?
2. Is the proposed change narrow?
3. Does it reduce repeated thesis-time waste?
4. Is there a validator?
5. Does doctor pass?
6. Does Justin need to approve because the target is risky, public, or
   external?

## Thesis Finitude Strategy

The dream system exists to help finish the thesis, not to become a hobby OS.
The current decision prior is:

```text
Choose the shortest safe evidence-backed action that most advances the
thesis-closing photon+jet artifact.
```

Before promoting dream output, waking Codex should classify the proposed work:

- `terminal_path`: directly closes a minimal publishable thesis artifact;
- `risk_reduction`: prevents wrong or fragile thesis-facing conclusions;
- `artifact_quality_multiplier`: improves plots/slides/scripts that are needed
  for thesis communication;
- `evidence_integrity`: strengthens provenance, validators, duplicate guards,
  or artifact registries;
- `workflow_compounding`: reduces future repeated work or orientation cost;
- `novelty_after_baseline`: interesting but should wait until baseline safety;
- `distraction_risk`: consumes attention without protecting the thesis path;
- `blocked_by_missing_evidence`: cannot progress until real evidence exists.

Novel ML, speculative physics, broad refactors, and big OS cleverness should
wait unless they directly unlock/protect the minimal publishable baseline or
are a small bounded test.

## Justin's Place In The System

Justin should not sift through the dream directories every morning.

Justin's normal surface is:

```text
Today's Plan
```

Justin should see:

- at most a short Overnight Dream Updates digest;
- only dream findings that change today, reveal a blocker, or need a decision;
- a clear `Needs Justin` line when approval or human judgment is actually
  needed.

Justin's role:

- set scientific direction and taste;
- approve risky changes;
- review plots/slides/results;
- correct Codex when the artifact, priority, or story is wrong;
- decide when a speculative idea becomes worth real work.

Justin should not be asked to:

- read all `learning_atoms.jsonl`;
- choose between internal validator designs;
- inspect every recurring-hotspot score;
- manually maintain Linear or register state;
- push a redundant approval button for reversible internal local-only
  maintenance that waking Codex can validate safely.

If a dream finding can be handled by waking Codex as a safe, local,
evidence-backed OS patch, Codex should handle it and report it. If it touches
science claims, production, public slides, external systems, priorities, or
anything irreversible, Justin must stay in the loop.

## Waking Codex's Place In The System

Waking Codex is the filter and implementer.

On a normal morning or after Justin asks what changed:

1. Read `AGENTS.md`.
2. Read `AGENTIC_OS_DREAMING.md` and this overview if dream behavior is in
   scope.
3. Run or inspect:

```bash
python3 scripts/os/dream/codex_dream_morning_summaries.py --window overnight --json
python3 scripts/codex_os_doctor.py --profile daily --json
```

4. If needed, inspect lane-local artifacts, not every run:

```text
morning_lane_summary.json
learning_atoms.jsonl
promotion_candidates.md
maintenance_debt.md
changed_actions.json
validation_summary.md
```

5. Decide whether there is:

- no user-facing action;
- a Today's Plan note;
- a register/Linear/task update;
- a small local OS patch;
- a Justin approval question;
- a real science/action blocker.

6. Promote only narrow, validated, evidence-backed changes.

## Manual Commands

Run one lane manually for debugging:

```bash
python3 scripts/codex_os_dream.py lane --lane-id <lane_id> --no-auto-maintain
```

Simulate a scheduled lane for validation only:

```bash
python3 scripts/codex_os_dream.py lane --lane-id <lane_id> --run-id validation-<name> --scheduled-morning-lane --no-auto-maintain
```

Validate a specific run:

```bash
python3 scripts/codex_os_dream.py validate agent_context/local/dreams/<run_dir>
```

Validate latest, but prefer exact paths during concurrent 03:30 lanes:

```bash
python3 scripts/codex_os_dream.py validate --latest --allow-missing
```

Aggregate scheduled-only morning summaries:

```bash
python3 scripts/os/dream/codex_dream_morning_summaries.py --window overnight --json
```

Run daily doctor:

```bash
python3 scripts/codex_os_doctor.py --profile daily --json
```

Resolve task-shaped context:

```bash
python3 scripts/os/context/codex_context_resonance.py resolve --task "<task>" --json
```

Check context-resonance health:

```bash
python3 scripts/os/context/codex_context_resonance.py canary --json
python3 scripts/os/context/codex_context_resonance.py review --json
python3 scripts/os/context/codex_context_resonance.py feedback-health --json
```

## What Should Reach Today's Plan

Today's Plan should receive only scheduled, eligible 03:30 lane summaries.

Good Today's Plan dream entry:

```text
Overnight Dream Updates
- Context Resonance: feedback loop stayed healthy; no tracked memory mutation.
- Science Scout: surfaced one proposal-only target-figure scaffold; no real
  science evidence or action until waking validation.
- Needs Justin: none.
```

Bad Today's Plan dream entry:

```text
Here are 80 lines of learning atoms, recurrence scores, and JSON fields...
```

The dream digest should help Justin decide what to do today, not force him to
debug the OS.

## What Should Stay Out Of Today's Plan

Keep these out of the daily surface unless something is broken:

- raw `learning_atoms.jsonl`;
- recurrence index details;
- homeostasis score mechanics;
- long maintenance-debt prose;
- per-lane validation command lists;
- local file-retention inventories;
- ChatGPT prompt-pack internals;
- duplicate versions of the same warning from multiple lanes.

Those belong in local dream artifacts and waking Codex summaries.

## Separation And Redundancy Rules

The system works only if each layer has a job.

Use this separation:

- Dream lanes: generate proposal-only pressure and local evidence.
- Doctor: aggregate, validate, and warn.
- Context resonance: retrieve/suppress task-shaped context.
- Register: canonical active work and ownership.
- Linear: durable workstream dashboard, not dream internals.
- Today's Plan: Justin-facing one-page cockpit.
- Waking Codex: promotion, implementation, and concise explanation.
- Justin: scientific approval, taste, priorities, and external/risky decisions.

Avoid these redundancies:

- Do not put every dream proposal into Linear.
- Do not put raw dream atoms into Today's Plan.
- Do not ask Justin to approve safe local OS hygiene that Codex can validate
  and report after the fact.
- Do not let multiple lanes open the same workstream unless the doctor or
  waking Codex merges the pressure into one action.
- Do not create a new dream lane when a memory trap, validator, or runbook
  would solve the problem.
- Do not promote speculative physics into active tasks until real inputs and
  null tests are named.

## When To Notify Justin

Notify Justin only when one of these is true:

- a dream found a stale active/running job or blocker that affects today's
  work;
- a validation or doctor failure needs action;
- a repeated issue is wasting thesis time and needs a decision;
- an external/risky mutation would be required;
- a top priority changed;
- a proposal is high-value but requires taste, science judgment, or approval.

Stay quiet or write a short digest when:

- dreams passed and only produced low-risk internal proposals;
- maintenance debt is known and unchanged;
- no scheduled eligible dream summary needs attention;
- all changes were local ignored indexes or proposal-only files.

## Debugging Checklist For A Blank Model

If the dream system looks broken:

1. Check the installed automation metadata under:

```text
/Users/patsfan753/.codex/automations/
```

2. Confirm the eight dream automations are `ACTIVE`, `kind = "heartbeat"`, and
   scheduled at `FREQ=DAILY;BYHOUR=3;BYMINUTE=30;BYSECOND=0`.
3. Confirm each has the expected `target_thread_id`.
4. Check the local binding file:

```text
agent_context/local/dreams/dream_lane_thread_bindings.json
```

5. Inspect the latest exact run directory for the affected lane.
6. Prefer exact validation paths over `--latest` during concurrent runs.
7. Run doctor.
8. If the failure is a policy/schema gap, patch the smallest tracked surface
   and rerun validation.
9. If the failure is external or risky, ask Justin instead of improvising.

## Current Design Judgment

The dream system is useful only if it reduces Justin's mental load. The best
shape is:

```text
many detailed local dream artifacts
-> one doctor/aggregator judgment
-> one short Today's Plan digest
-> waking Codex handles safe internals
-> Justin sees only decisions, risks, and thesis-facing next actions
```

The most important future improvement is not more dream volume. It is better
selection: fewer surfaced items, stronger evidence refs, cleaner promotion
rules, and more automatic consolidation of duplicate lane pressure into one
actionable waking step.

