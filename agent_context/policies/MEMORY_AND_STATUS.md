# Memory And Status

## Memory Stores

- `agent_context/TASK_BOARD.md`: clean milestone/blocker board.
- `agent_context/STATUS_DASHBOARD.md`: traffic-light dataset/job/output
  dashboard.
- `agent_context/REFERENCE_MAP.md`: compact reference decision map.
- `agent_context/SLIDE_STYLE_MAP.md`: durable slide design language.
- `codex_notes/PROJECT_BOARD.md`, `DATASET_STATUS.md`, `KNOWN_ISSUES.md`,
  `RUN_LOG.md`: older project-state ledgers, especially dataset validity.

## What To Record

Record durable decisions and evidence:

- dataset choices and sample pairs;
- job IDs, DAG IDs, submit node, campaign tags;
- commands that worked;
- output validity/staleness;
- physics conclusions;
- recurring mistakes and hard stops;
- cleanup decisions and protected paths.

Do not record speculative guesses, bulky conversational notes, or casual task
ideas unless Justin asks.

## Trigger Words

If Justin says "track", "add", "remember", "record", "note this", or similar,
update the relevant memory file. If he casually mentions a possible task, ask
before adding it.

When a task appears finished, move it to `Done Pending Removal` and ask before
removing or archiving.

## Campaign Recording

When Codex gives a submission command or sees Justin submit one, immediately
record:

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
