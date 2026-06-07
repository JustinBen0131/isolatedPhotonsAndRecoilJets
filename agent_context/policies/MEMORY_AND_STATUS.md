# Memory And Status

## Memory Stores

- `agent_context/CODEX_WORK_REGISTER.yaml`: canonical active-work register for
  Codex ownership, Linear sync, Today's Plan anchors, active jobs, next checks,
  stale-after timing, and compact handoffs.
- `agent_context/ARTIFACT_REGISTRY.yaml`: canonical provenance registry for
  reusable plots, ROOT files, CSV/JSON summaries, slide candidates, scripts,
  and OS evidence bundles.
- `agent_context/TASK_BOARD.md`: clean milestone/blocker board.
- `agent_context/STATUS_DASHBOARD.md`: traffic-light dataset/job/output
  dashboard.
- `agent_context/REFERENCE_MAP.md`: compact reference decision map.
- `agent_context/SLIDE_STYLE_MAP.md`: durable slide design language.
- `agent_context/THESIS_NARRATIVE_MAP.md`: durable hierarchy from thesis goal
  to evidence, plots, jobs, and slide/story obligations.
- `agent_context/memory/`: tracked neutral memory-routing registries:
  `CONTEXT_RESONANCE_INDEX.yaml`, `SCHEMA_REGISTRY.yaml`, and
  `NEGATIVE_MEMORY_MAP.yaml`. These contain sanitized cue records and source
  pointers, not raw transcripts or private remote details.
- `agent_context/local/os_events.jsonl`: private append-only event ledger for
  guard preflights, approvals, snapshots, warnings, repairs, and accepted
  risks. Do not commit it.
- `agent_context/local/context_resonance/`: private salience, retrieval-outcome,
  and context-interference ledgers. Use it for distilled memory-routing
  evidence only, not raw chat transcripts.
- `codex_notes/PROJECT_BOARD.md`, `DATASET_STATUS.md`, `KNOWN_ISSUES.md`,
  `RUN_LOG.md`: older project-state ledgers, especially dataset validity.

## What To Record

Record durable decisions and evidence:

- active workstream ownership and the current Codex chat/session label;
- dataset choices and sample pairs;
- job IDs, DAG IDs, submit node, campaign tags;
- commands that worked;
- output validity/staleness;
- physics conclusions;
- recurring mistakes and hard stops;
- cleanup decisions and protected paths.
- agent/OS failures, near misses, and repeated friction points that should
  become a policy, script, preflight, Linear follow-up, or accepted risk.
- reusable plot, slide, script, speaker-note, or workflow corrections after
  Justin feedback, when the lesson can prevent the same mistake in future work.

Do not record speculative guesses, bulky conversational notes, or casual task
ideas unless Justin asks.

For context resonance, record whether a surfaced memory actually helped,
polluted the answer, or should be suppressed next time. The goal is better
retrieval quality: compact active facts, useful latent nudges, negative-memory
traps, and stale/synthetic context kept cold.

Use the waking resolver before a context-sensitive task:

```bash
python3 scripts/os/context/codex_context_resonance.py resolve --task "<task>" --json
```

After a retrieved memory materially affects the task, append local-only
feedback with `record-outcome`. Do not put raw ChatGPT transcripts, SDCC
private paths, passwords, or bulky chat text into tracked memory registries.

## Trigger Words

If Justin says "add a task", "add this task", "add this to my task list", "put
this on my todo list", "put this in the backlog", "track this as a task", "make
this a Linear task", "remember this as a task", or similar, run the full task
capture pipeline:

1. update `agent_context/CODEX_WORK_REGISTER.yaml` first;
2. update or create the matching Linear workstream issue;
3. update Today's Plan only if the task belongs on the daily cockpit;
4. update local status/memory notes only when there is durable evidence,
   project state, job state, or a science decision to preserve;
5. leave one next action and a compact handoff.

If Justin says "track", "remember", "record", "note this", or similar for a
fact that is not a task, update the relevant memory/status file. If he casually
mentions a possible task without task-capture wording, ask before adding it.

When a task appears finished, move it to `Done Pending Removal` and ask before
removing or archiving.

## Codex Work Register

For meaningful multi-step work, update
`agent_context/CODEX_WORK_REGISTER.yaml` before starting and again at closeout.

Record:

- the workstream being claimed or created;
- current next action;
- active Codex session or heartbeat label;
- Linear issue or `linear_sync: pending`;
- Today's Plan anchor or `today_doc_anchor: null`;
- active jobs, artifacts, evidence, `next_check`, and `stale_after`;
- compact handoff summary for the next Codex chat.

If the matching workstream is an umbrella campaign and the current deliverable
can finish without completing the campaign, record or create a focused child
workstream instead of overloading the parent. The child should own the active
Codex session, heartbeat automation, active jobs, artifact contract, QA
criteria, and done criteria. The parent should keep only the campaign goal,
child map, durable context, and campaign-level decision state.

Optional critical-path fields may be added when they improve routing or
decision quality. Do not add them mechanically to every task:

- `critical_path_class`: `terminal_path`, `risk_reduction`,
  `artifact_quality_multiplier`, `evidence_integrity`,
  `workflow_compounding`, `novelty_after_baseline`, `distraction_risk`, or
  `blocked_by_missing_evidence`.
- `terminal_artifact_supported`: the rung of the thesis artifact ladder this
  work supports.
- `minimal_publishable_path_impact`: low/medium/high/critical impact on the
  shortest safe publishable path.
- `regret_if_delayed`: what future Justin would lose if this waits.
- `opportunity_cost`: what thesis work or attention this displaces.
- `safe_next_action`: the smallest evidence-backed next step.
- `novelty_gate`: why novelty is allowed now or why it should wait.
- `defer_until_baseline_safe`: true when interesting work should wait for
  baseline safety.
- `next_thesis_closing_action`: the next action that most directly closes the
  thesis path.

If Justin corrects repeated drift that costs thesis time, classify the drift
and encode the smallest durable prevention mechanism: policy, negative memory,
schema, validator, style-map entry, artifact QA rule, runbook/skill,
work-register update, or postmortem.

If another active session already owns the same workstream, stop before doing
duplicate work and report the existing owner, next check, and evidence.

Google Docs and Linear are projections of the register. Do not rely on either
as the only source of active job truth.

## Campaign Recording

When Codex gives a submission command or sees Justin submit one, immediately
record:

- the matching workstream in `CODEX_WORK_REGISTER.yaml`;
- submit host/node;
- exact command;
- dataset/mode;
- timestamped run roots;
- DAG path and Condor cluster IDs;
- report/output roots;
- next expected checkpoint.

Use visible terminal prompt node as evidence when present.

### AuAu Campaign Timing And Efficiency

For AuAu RecoilJets campaigns, every heartbeat should treat the campaign as a
measured system, not just a Condor queue. Record submitted runs/files/jobs,
`groupSize`, requested memory, active/running/held counts, non-tiny ROOT
coverage, complete-run count, output byte growth since the previous checkpoint,
recovery-target status, and known failure modes. "Finished" means ROOT coverage
and sanity: non-tiny, non-zombie, expected keys/histograms, and merge-readiness,
not Condor history alone. Use
`scripts/sdcc/runtime/audit/auau_campaign_checkpoint.py` for compact
checkpoint JSON/text where possible, and store sidecar timing results
separately from mergeable campaign outputs before changing future broad
submission defaults.

## Evidence Discipline

Do not mark outputs valid/stale/ready/done without evidence. Prefer exact
paths, timestamps, ROOT object names, job IDs, Gmail subject lines, terminal
output, and compact command output.

## Cleanup Memory

For SDCC tests or sidecar workflows, record active run roots and cleanup
decisions while live. When a test is superseded, failed, or no longer needed,
remind Justin which timestamped artifacts can be cleaned after confirming no
live jobs reference them.

Prompt cleanup at the end of major tasks, production passes, or plot sets.
Never silently delete useful context.

## Failure-To-Policy Loop

When a mistake or near miss repeats, do not just apologize or remember it in
chat. Convert it into one durable mechanism:

- policy patch when the issue is behavioral;
- script or doctor check when it is mechanically detectable;
- task/register correction when ownership or next action was wrong;
- Linear follow-up when the fix needs later review;
- postmortem note when the failure could recur or caused user-visible waste.

Use `agent_context/templates/OS_POSTMORTEM_TEMPLATE.md` for postmortems and
run `python3 scripts/codex_os_doctor.py --profile daily` after OS-state
repairs. Repeated warnings should become encoded prevention, not permanent
daily noise.

## Feedback-To-Learning Loop

When Justin corrects a generated plot, full-slide candidate, speaker script,
task surface, or workflow behavior, classify the correction before moving on.
The categories are:

- style/narrative;
- physics/science;
- provenance/evidence;
- workflow/safety;
- retrieval failure;
- memory salience failure;
- task/status failure;
- tool/runtime failure;
- artifact-quality failure.

Choose exactly one primary durable mechanism when the lesson is reusable:
policy patch, style-map entry, negative-memory trap, schema/validator,
artifact QA rule, script/default change, register update, postmortem, or
Linear follow-up. Do not record every preference; record corrections that
prevent repeated user friction, scientific ambiguity, unsafe action, or
artifact-quality regressions.

Every durable reflection must cite concrete evidence: user correction, file
path, slide ID, artifact path, command output, validator result, or register
entry. Use:

- `agent_context/templates/OS_REFLECTION_TEMPLATE.md` for OS behavior lessons;
- `agent_context/templates/SLIDE_FEEDBACK_TEMPLATE.md` for slide/plot/script
  correction loops;
- `agent_context/templates/TRAJECTORY_LEARNING_TEMPLATE.md` for multi-step
  improvements across a campaign or dream lane.

Use this taxonomy when deciding where memory belongs:

- evidence: verified files, commands, logs, artifacts, user statements;
- experience: observed friction, repeated correction, failed workflow;
- belief: provisional interpretation that still needs verification;
- style: Justin-facing visual, wording, hierarchy, or script preferences;
- policy: hard behavioral rule or safety condition;
- skill/procedure: repeatable method, checklist, script default, or validator.

Dream learning atoms are not durable memory by themselves. Treat
`learning_atoms.jsonl` as a proposal queue generated from preserved raw dream
episodes and real evidence pointers. A waking session may promote one atom only
when it has a raw episode path, concrete evidence refs, a narrow target, a
validator or doctor check, and no synthetic-approval ambiguity. Promotion may
write exactly one or more of: policy line, negative-memory trap,
context-resonance cue, schema/validator, slide style-map rule, artifact QA
rule, runbook/skill contract, work-register correction, or postmortem.

Do not overwrite raw episodes or replace evidence with a consolidated atom.
Unpromoted atoms cool through the dream recurrence index; they are not deleted
automatically. If the same dream finding appears as prose for three scheduled
nights without an atom, validator, explicit rejection, or cooldown decision,
record it as maintenance debt rather than letting it remain recurring noise.

## Workstream Refresh Protocol

Use `register_workstream_refresh_contract` when live workstreams repeatedly
create status, active-job, stale-state, or evidence pressure. The read-only
check is:

```bash
python3 scripts/os/register/codex_work_register_stale.py agent_context/CODEX_WORK_REGISTER.yaml --protocol
```

It classifies live workstreams as `current`, `waiting`, `stale`,
`needs_evidence`, or `archive_review`. Any real register update must preserve
exact evidence and update `last_verified`, `next_check`, `stale_after`, and
active-job fields through waking validation; the dream must not close tasks or
change science status automatically.
