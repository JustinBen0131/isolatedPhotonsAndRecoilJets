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
