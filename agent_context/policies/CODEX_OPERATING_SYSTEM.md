# Codex Operating System

This policy is the strict coordination layer for ThesisAnalysis. It exists so
separate Codex chats, heartbeats, Linear, Google Docs, and repo ledgers all
agree on what is active, running, waiting, blocked, or done.

## Core Rule

`agent_context/CODEX_WORK_REGISTER.yaml` is the canonical active-work register.
Google Docs and Linear are projections of that register. They are not the
source of truth.

For OS hardening, self-tuning workflow, symbiotic human-agent operation, or
thesis control-plane work, also load
`agent_context/policies/AGENTIC_OS_HARDENING.md` and
`agent_context/THESIS_NARRATIVE_MAP.md`.

For synthetic dream rehearsal, also load
`agent_context/policies/AGENTIC_OS_DREAMING.md`. Dream outputs are private
proposal artifacts under `agent_context/local/dreams/`; they are not user
approval and must not be promoted silently into real task state.
The canonical overnight coordinators are the eight lane commands
`python3 scripts/codex_os_dream.py lane --lane-id <lane_id>`, one for each
dream lane at `03:30`. Morning-priorities automation is a downstream waking
projection step, not part of the overnight mutation boundary.
For compact retrieval guidance before a task, use:

```bash
python3 scripts/os/context/codex_context_resonance.py resolve --task "<task>" --json
```

It returns a small set of `conscious_context` / compatibility `active_facts`,
latent nudges, negative-memory traps, suppressed context, and evidence checks.
It is a waking router, not proof; load and verify the cited files before
making factual claims. Tracked cue records live under `agent_context/memory/`;
dynamic retrieval feedback stays private under
`agent_context/local/context_resonance/`.
The Stage 9 evolutionary-maintenance layer is now a guarded internal
self-maintenance loop. Dream mode may automatically perform `auto_safe`
untracked generated-junk cleanup and refresh ignored local dream/research
artifacts, and it must record those actions in `changed_actions.md/json`.
`internal_evolution_queue.md` / `internal_evolution_queue.json` are the
control surface for `auto_validated` tracked OS/policy/runbook/index work:
waking Codex may apply those small internal improvements automatically when
the change is local-only, evidence-preserving, reversible, does not alter
science status or external systems, has exact validation commands, and passes
the doctor. SDCC edits/cleanup, Condor/job control, scientific output/model
movement or deletion, Google Drive/Slides/Gmail/Linear mutation, and task
closure or physics-status changes remain waking guarded actions.

For delegated ChatGPT UI research, also load
`agent_context/policies/ASK_CHATGPT_DELEGATION.md`. ChatGPT output is external
research/critique, not evidence or approval, and transcripts stay private under
`agent_context/local/chatgpt_research/`.
For `instant`, `thinking`, and `heavy`, Codex may collect and synthesize the
response. For `pro`, `extended pro`, or `deep research`, Codex submits one
sanitized prompt and stops; Justin pastes the completed response back before
Codex verifies or integrates it.
The explicit dual-Pro overnight research lane is the only current exception:
when Justin asks for that lane, Codex may stage two `pro_extended` packs, submit
them in separate fresh ChatGPT conversations through Computer Use, run one
temporary local heartbeat until both copied responses are present, and then
verify locally. The OS response may drive at most one reversible OS-only patch
after rollback and doctor validation; it must not change science state,
external systems, task status, or overall functionality.

For thesis-facing artifact provenance, use
`agent_context/ARTIFACT_REGISTRY.yaml`. For risky mutations, use the safety
kernel:

```bash
python3 scripts/codex_os_snapshot.py create --label "<scope>"
python3 scripts/codex_os_guard.py preflight ...
```

Before starting any meaningful multi-step work, Codex must either claim an
existing workstream in the register or create a new one.

Meaningful work includes:

- SDCC status checks, submissions, watchdogs, production, validation, merges,
  or transfers;
- new plot batches, slide candidates, Google Slides edits, or deck QA;
- training, ablation, model selection, stitching, or approval-facing analysis;
- heartbeat setup or follow-up that can leave work active after the chat ends;
- multi-chat or delegated work.

Tiny one-shot answers do not need a register claim unless they change project
state or create a follow-up obligation.

## Task Capture Trigger

When Justin says any close variant of:

- "add a task";
- "add this task";
- "add this to my task list";
- "put this on my todo list";
- "put this in the backlog";
- "track this as a task";
- "make this a Linear task";
- "remember this as a task";

Codex must treat the wording as an implementation command, not a suggestion.

Required task-capture pipeline:

1. Decide whether the request belongs to an existing workstream or needs a new
   durable workstream.
2. Update `agent_context/CODEX_WORK_REGISTER.yaml` first. The register is the
   canonical task source.
3. Sync Linear second: update the existing issue or create one issue for the
   durable workstream. Do not create noisy command-level issues.
4. Update Today's Plan only when the task is active today, waiting/review and
   worth surfacing, or changes the Top 3. Backlog-only tasks should usually
   stay in Linear detail and the register.
5. Update `TASK_BOARD.md`, `STATUS_DASHBOARD.md`, or other local status notes
   only when the task changes durable project state, evidence, jobs, or
   scientific decisions.
6. Leave a compact handoff in the register and, when Linear is available, a
   Linear comment if the task changed ownership, blockers, active jobs, or next
   action.

If Linear or Google Docs tools are unavailable, keep the register current and
set the affected `linear_sync` or `today_doc_anchor` state to pending/null.

Do not ask before applying this pipeline when the user uses explicit
task-capture wording. Ask only when the requested task is ambiguous enough that
placing it in the wrong workstream would create real confusion, or when carrying
out the task would require a risky mutation such as SDCC submission, deletion,
Slides mutation, credential handling, or job control.

Casual future-tense discussion is different. If Justin merely mentions that
something might be useful someday and does not use task-capture wording, ask
before adding it.

## Start-Work Gate

Before real work:

1. Read `AGENTS.md`.
2. Load this policy and any route-specific policy files.
3. Read `agent_context/CODEX_WORK_REGISTER.yaml`,
   `agent_context/TASK_BOARD.md`, and `agent_context/STATUS_DASHBOARD.md`.
4. Apply `DUPLICATE_RUN_GUARD.md` for any run, training, validation, merge,
   production-style analysis, or artifact that could already exist.
5. Claim or create exactly one workstream for the action.
6. Set one concrete `current_next_action`.
7. Record the chat/session label when available.
8. Record active jobs, artifacts, evidence, `next_check`, and `stale_after`
   whenever the work can continue after the current response.
9. Mirror state to Linear and Today's Plan when the tools are available. If
   Linear is unavailable, set `linear_sync: pending`.
10. For risky mutations, run the safety guard with exact scope, duplicate
    fingerprint when relevant, snapshot ID, and explicit approval before the
    mutation.

If a task is already claimed by another active Codex session, stop and report
the existing workstream, owner/session, next check, and evidence before doing
duplicate work.

## Doctor Gate

Before calling OS state healthy, before closing a major workstream, or after
changing register/Linear/daily-plan policy, run:

```bash
python3 scripts/codex_os_doctor.py --profile daily
```

The doctor is the executable backstop for the soft rules in this policy. A
passing doctor does not prove the physics is correct, but a failing doctor
means the operating system is not healthy enough to trust without repair.

Use `--profile strict` before claiming the OS is release-clean. Use
`--profile release` when testing the OS machinery itself because it includes
synthetic failure drills.

## Register Fields

Each workstream should keep these fields current when applicable:

- `workstream_id`: stable lowercase identifier.
- `title`: human-readable workstream name.
- `status`: `active`, `running`, `waiting`, `blocked`, `review`, `backlog`,
  `done_pending_review`, or `archived`.
- `priority`: `P0`, `P1`, `P2`, or `P3`.
- `goal`: what the work unlocks.
- `current_next_action`: one concrete next step.
- `depends_on`: prerequisite workstreams or external blockers.
- `active_codex_session`: chat, heartbeat, or agent label if active.
- `chat_label_or_thread`: user-facing chat/thread label when known.
- `linear_issue`: Linear issue key or URL when synced.
- `linear_sync`: `synced`, `pending`, or `not_needed`.
- `today_doc_anchor`: Today's Plan section/anchor when available.
- `workstream_spec`: durable Google Doc, Markdown, or Drive spec for complex
  studies.
- `active_jobs`: live DAGs, clusters, watches, or long-running checks.
- `artifacts`: outputs, roots, PNGs, CSVs, manifests, or decks.
- `evidence`: concise source lines supporting status.
- `last_verified`: date or timestamp of last evidence check.
- `next_check`: when Codex should look again.
- `stale_after`: timestamp/date when lack of update should be flagged.
- `handoff_summary`: compact state for the next Codex chat.

## Status Semantics

- `active`: ready to move now.
- `running`: external work is in flight; include job/watch evidence.
- `waiting`: waiting on approval, collaborator, connector, output, or user.
- `blocked`: cannot move until a named blocker is cleared.
- `review`: artifact exists and needs inspection or decision.
- `backlog`: important but not today's lane.
- `done_pending_review`: appears complete but Justin or evidence cleanup has
  not approved archival.
- `archived`: no longer shown as active or rolled into Today's Plan.

Do not mark work done without concrete evidence or a clear Justin statement.

## Done Protocol

Before moving a workstream to `done_pending_review`, Linear `Done`, or archive:

1. record concrete completion evidence;
2. clear `active_jobs` unless they are explicitly historical and moved into
   evidence;
3. clear `active_codex_session`;
4. set `current_next_action` to no standing action or an exact reopening
   condition;
5. remove it from Today's Plan active/waiting/review surfaces unless Justin
   asks for visibility;
6. update Linear with a completion note when the connector is available;
7. for artifact-producing work, register canonical outputs or explicitly state
   why no registry entry is needed;
8. for risky mutations, confirm the guard event and snapshot are recorded in
   `agent_context/local/os_events.jsonl`;
9. run `python3 scripts/codex_os_doctor.py --profile daily`.

If any of those steps fail, keep the workstream live or in review.

## Linear Projection

Codex operates Linear for Justin. Justin should not have to manually maintain
Linear.

Use one Linear issue per durable workstream. Do not create one issue per shell
command, plot attempt, or small chat action.

Linear's job is to be the detailed dashboard behind the one-page Today's Plan.
It should be useful to open, but Justin should not need to edit it.

### Campaign And Chat Hierarchy

For umbrella goals that contain multiple slides, jobs, validations, or chats,
Linear should have one campaign-level issue. A campaign issue is the user-facing
home base for the whole effort. Child/related workstream issues are allowed
only when they represent durable subgoals, not individual shell commands.

Every active or running workstream must identify:

- `campaign_id`: stable umbrella identifier, when the work belongs to a
  campaign.
- `chat_label_or_thread`: the short human-readable Codex chat name Justin can
  recognize.
- `active_codex_session`: the current chat/thread/heartbeat label when Codex
  is actively working.
- `RJ_CODEX_CHAT_NAME` and `RJ_CODEX_THREAD_ID` for SDCC submissions that can
  emit pipeline emails.

Linear issue descriptions and comments for active/running work must include a
`Codex chat map` block with the current chat label, thread id when available,
heartbeat automation id if one is monitoring the work, and active campaign
tags/DAGs. If Codex does work in a chat that changes a campaign, Codex must
update the matching Linear campaign/workstream issue before final response, or
mark `linear_sync: pending` in the register if the connector is unavailable.

Suggested chat names should be short and campaign-scoped:

- `CampaignShort | Action | YYYYMMDD` for ordinary chats.
- `CampaignShort | Watch <tag>` for heartbeat/watchdog chats.
- `CampaignShort | Slide <n>` for slide-production chats.

Example: `PPG12 Match | Watch npbet5env`.

### Linear UX Doctrine

Use scan-first titles. The first words should tell Justin why the issue matters
before he reads the rest:

- `P0 Campaign | ...`: umbrella campaign coordinating multiple workstreams,
  chats, jobs, and slide deliverables.
- `P0 Approval | ...`: deadline-governing approval work.
- `RUN Stitching | ...`: live stitching, Condor, or watchdog work.
- `P0 ML | ...`: active final-model-choice work.
- `P1 Closure | ...`: important explanation, validation, or closure work.
- `WAIT Validation | ...`: important validation that is waiting on named
  prerequisites.
- `BKL ... | ...`: preserved backlog that should not distract today.
- `OS Baseline | ...`: meta-infrastructure for the operating system.

Use labels as a visual grammar, not decoration:

- Surface labels: `Surface: Today`, `Surface: Linear Detail`.
- Area labels: `Area: Approval`, `Area: Stitching`, `Area: ML`,
  `Area: Trigger`.
- Ops labels: `Ops: Active Jobs`, `Ops: Health Check`.
- Decision labels: `Needs: Justin Review`, `Slide-facing`.
- Priority/status labels: `P0`, `P1`, `P2`, `Running`, `Review`,
  `Backlog`, `Blocked`, `Waiting`.

Every durable issue should include, in this order:

- `Goal`: what the work unlocks.
- `Current next action`: one concrete next step.
- `Why this matters`: the thesis or approval consequence.
- `Codex chat map`: active chat label/thread id, heartbeat id, campaign tag,
  and any active DAG clusters.
- `Depends on`: exact blocker, issue link, or `none`.
- `Evidence`: one compact source line.
- `Codex operating rule`: what future Codex sessions must preserve.
- `Daily surface`: whether the item belongs in Today's Plan today.

Linear lanes should map to register statuses:

- `active` -> In Progress with `Surface: Today` if it can move today.
- `running` -> In Progress with `Running` and `Ops: Active Jobs`.
- `waiting` -> waiting label/state when available, otherwise In Progress plus
  `Waiting`.
- `blocked` -> blocked label/state when available, otherwise In Progress plus
  `Blocked`.
- `review` -> In Progress plus `Review` and usually `Needs: Justin Review`.
- `backlog` -> Backlog plus `Surface: Linear Detail`.
- `done_pending_review` -> Done only after evidence and handoff are recorded.

Each Linear issue should include goal, current next action, blocker or
dependency, active jobs, Today's Plan link, workstream spec link, and latest
Codex handoff comment.

When a workstream is active in a Codex chat, add or update a comment with:

- what Codex is doing now;
- what files, jobs, docs, or slides it touched;
- what is still running or waiting;
- one next check or handoff action.

If the Linear connector is unavailable, keep the register current and set
`linear_sync: pending`. Sync later when tools are available.

## Today's Plan Projection

Today's Plan is Justin's human cockpit. It should be short, visual, and
generated from active register state.

Daily doc top structure:

- `Today's Work Contract`: 3-4 short rules that keep the page from becoming an
  archive.
- `WORKING POINT`: the current working-point slide deck/link. It must stay near
  the top of every Today's Plan. Rollover should copy it from the previous
  daily doc unchanged unless Justin edits it or explicitly provides a new link.
  If Justin references slides without specifying another deck link, Codex
  should treat this `WORKING POINT` deck as the default slide context.
  Before creating or regenerating a daily plan, read the previous/current
  Today's Plan and extract its `WORKING POINT:` line. If that line differs from
  the register, treat the doc as Justin's latest edit and update the register
  before rendering the next plan.
- `Do First Today`: at most three concrete actions, written as outcome plus
  next action.
- `Morning Health Check`: a compact read-only operational snapshot, including
  `/tg/tg01/bulk/jbennett` usage, Justin-visible `/sphenix/user` usage when
  available, and `./checkCondorQ` queue status for Justin's jobs.
- `Slide Drivers`: active working-point deck, HP26 baseline deck when
  relevant, and the current slide/narrative rule.
- `Active Now`: only the active/running workstreams that can affect today.
- `Waiting / Review`: blocked, waiting, and review items with one next
  condition each.
- `Backlog - Do Not Let This Crowd Today`: important preserved work that should
  not compete with the top three.
- `Command Links`: Linear board, spillover doc, operating rules, and key specs.
- `Scratch / Inbox`
- `Do Not Resurface As Active`

Visual contract:

- body text black, normal weight, 11 pt;
- links remain blue and underlined;
- section headers use restrained color bands and varied heading colors, not a
  monotone blue stack;
- status colors are green for Active, amber for Running/Review, blue-gray for
  Waiting, red for Blocked, gray for Backlog, muted green for Done;
- bold only task titles, `Next:` labels, status words, and urgent dates;
- the first screen should show the work contract, Top 3, and enough health/slide
  context to act without scrolling through the full history;
- do not roll completed history forward.

Useful visual mutations should propagate. If Codex tries an unrequested but
safe planning-surface improvement and Justin explicitly says he likes it, treat
that as approval to preserve the pattern in the OS. Record the pattern in the
register or policy, and use it in future daily cockpit renders unless Justin
later rejects it. As of 2026-06-03, the approved pattern is an HTML-first local
daily cockpit source under `agent_context/local/daily_plans/`, imported to a
native Google Doc and then patched with explicit Google Docs styling. This
pattern is preferred over raw markdown because it gives Justin a faster visual
scan surface and creates a reusable local artifact for future render/template
improvements.

If Google Drive import rejects markdown or strips formatting, do not leave a
plain-text dump as the daily plan. Either apply a native Google Docs styling
pass immediately from connector readback (title, colored section bands, real
bullets, neutral body text, status colors, surgical bolding) or stop and report
that the formatted cockpit could not be created. A raw markdown/plain-text
Today's Plan is not an acceptable completed morning surface.

Complex studies belong in linked workstream specs, not expanded inside the
daily cockpit.

The morning health check should stay short enough to scan in under 20 seconds:

- quota lines should show used, available, percent used, and any obvious
  pressure point;
- Condor lines should show running, idle, held, and total job counts;
- active job groups should include the registered workstream, chat/session
  label, intention, and next check;
- if live SDCC/Condor access is unavailable, write `Not checked` plus the
  reason instead of omitting the section;
- the health check is read-only. Do not clean, hold/release/remove jobs,
  submit recovery, or transfer files from this check without explicit approval.

## Handoff And Closeout

At the end of meaningful work, update the register with:

- final status;
- exact evidence paths, job IDs, or artifact links;
- whether Linear and Today's Plan are synced or pending;
- one next action;
- cleanup or duplicate-run implications;
- handoff summary for the next Codex chat.

If work remains active or running, make sure `next_check` and `stale_after` are
set.

## Staleness Rule

Any `active`, `running`, `waiting`, `blocked`, or `review` workstream with an
expired `stale_after` should be surfaced by the heartbeat. Notify only when a
stale item, blocker, completed job, changed priority, or new evidence requires
Justin's attention.
