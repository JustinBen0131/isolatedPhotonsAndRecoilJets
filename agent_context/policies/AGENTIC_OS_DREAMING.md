# Agentic OS Dreaming

This policy defines the ThesisAnalysis dream layer: private, synthetic,
proposal-only simulation that helps Codex rehearse likely future failures and
prepare better next-day responses without touching live external systems.

## Core Rule

Dreams are not user intent. A synthetic Justin line is a test fixture, never
approval, instruction, or evidence of what Justin actually said.

The dream system may write only under:

```text
agent_context/local/dreams/
```

It must not mutate SDCC, Condor, Gmail, Google Drive, Google Slides, Linear, or
real repo-tracked files. It may generate proposed patches and task updates in
the dream directory for waking review.

## Cadence

- Micro dream is retired from the default contract. Re-enable it only by an
  explicit policy change plus a real installed automation; historical mentions
  of micro-dreams do not make them current.
- Nightly dream: one deeper overnight pass. It replays recent work, simulates
  tomorrow's likely status/provenance pressure, and produces a morning report.
- Nightly super-heartbeat: the default overnight operating mode. It runs the
  nightly dream, validates it, runs the doctor and read-only support checks,
  and emits one consolidated nightly package plus one stable overnight chat.
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
When the nightly super-heartbeat orchestrator is used, the same dream run
directory should also contain `nightly_heartbeat.md`,
`nightly_heartbeat_signal.json`, and `morning_appendix.md` so one overnight
artifact package carries the whole read-only story.

## Safe Overnight Work

Dream-mode overnight work should make waking Codex better without changing the
analysis state. Allowed overnight work is limited to:

- read-only local register, task, artifact, policy, and provenance audits;
- local dream simulation and dream validation;
- local doctor/stale/radar/artifact-registry checks;
- cleanup candidate lists for local memory/dream/research-pack clutter;
- SDCC clutter candidate proposals from recorded evidence only;
- physics scenario proposals based on already usable data provenance;
- local literature-scout maps from already available local papers and notes;
- synthetic target-figure sketches that are clearly marked as non-data;
- sanitized ChatGPT research packs through `ASK_CHATGPT_DELEGATION.md`.

Forbidden overnight work:

- Condor submissions, job removal, merge reruns, production analysis, or broad
  plot campaigns;
- SDCC deletion, movement, transfer, remote edits, or queue control;
- Gmail, Drive, Slides, Calendar, Linear, or repo-tracked mutations;
- applying cleanup, memory deletion, or task promotion automatically;
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

Cleanup output is always a proposal. It must include evidence, exact candidate
paths or scopes when known, risk of deleting, safe verification commands, and
the approval needed before any real cleanup. Dream mode must not run `rm`, move
files, mark memory obsolete, or touch SDCC storage.

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
- `nightly_heartbeat.md`
- `nightly_heartbeat_signal.json`
- `morning_appendix.md`
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

The goal is thorough private synthesis, not verbose heartbeat output. The
dream may spend 10-20 minutes producing detailed local artifacts, but the
thread notification should report only the top findings, report path, top
morning actions, and any validation or doctor failure unless user action is
needed. The default overnight user-facing surface is one stable nightly
heartbeat thread, not multiple fragmented dream/doctor chats.

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
