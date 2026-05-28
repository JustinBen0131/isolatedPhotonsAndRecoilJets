# RecoilJets Gmail And Watchdogs

## Gmail Pipeline Label

Gmail label `RecoilJets Pipeline` is the live dashboard for
RecoilJets/SDCC/Condor/DAGMan notifications.

Unread hygiene is mandatory for this label. Codex should not leave old
pipeline emails unread after reading and digesting them into chat, Linear,
`STATUS_DASHBOARD.md`, `TASK_BOARD.md`, or `CODEX_WORK_REGISTER.yaml`.
Unread RecoilJets Pipeline messages must mean one of two things only:

- new/unconsumed information still needs triage; or
- Codex intentionally left the message unread and reported exactly why.

Stale consumed emails must not remain unread as future false alarms.

Treat messages matching `RECOILJETS_STAGE_EMAIL_V1`, `[RecoilJets]`, DAGMan,
Condor, HTCondor, `READY`, `CHECK`, `FAILED`, held, rescue, or removed as
pipeline status.

When reading a new pipeline email, extract:

- dataset, stage, status;
- `codex_chat_name` and `codex_thread_id` from the `CODEX SUBMISSION` block;
- manifest/output paths;
- DAG/log/error paths;
- rescue count;
- profile summary;
- next action;
- job/cluster IDs.

After consuming and reporting or recording the useful information, mark only
that consumed pipeline email read. Do not delete or archive fresh actionable
pipeline messages unless Justin asks.

If an email reports `FAILED` or `CHECK`, marking the email read does not mean
the stage is resolved.

After marking consumed pipeline messages read, verify the `RecoilJets Pipeline`
label unread count before closing the status check. If any unread pipeline
messages remain, report whether they are new/unconsumed or intentionally left
unread, so stale failure alerts do not survive as ambiguous future state.

For active watchdog heartbeats, the normal loop is: Gmail first, digest useful
pipeline evidence, update the relevant local/Linear state if the evidence
changes status, mark consumed messages read, then verify the label unread count
is clean. Do not rely only on Condor/ROOT checks while old consumed pipeline
emails remain unread.

If a newly consumed Codex-launched RecoilJets email has
`codex_chat_name=unknown` or `codex_thread_id=unknown`, treat that as a
provenance defect. It does not by itself mean the physics job failed, but it
means the originating submit command did not propagate the required
`RJ_CODEX_CHAT_NAME`/`RJ_CODEX_THREAD_ID` metadata into the remote DAG creation
environment. Report the affected email subjects/stages and fix the next submit
driver before any further Codex-launched rerun.

## Gmail Connector Auth Failures

If a Gmail tool call fails with `token_expired`, `401`, expired OAuth, or a
similar connector-auth error:

- Do not claim Gmail was inspected, searched, or marked read.
- Do not call the plugin broken or say the Gmail integration needs to be
  revamped. In long-running chats, this is often chat-session credential
  expiration rather than a real connector failure.
- Retry once using exact label IDs from `list_labels` if that call is
  available. If Gmail still fails, tell Justin: "Please fork/reopen this chat;
  Gmail auth may refresh cleanly in the new chat." Treat forking as the first
  recovery path before asking Justin to reconnect the Gmail plugin.
- State the exact auth error once, then continue with safe read-only fallback
  evidence such as Condor queue state, DAG logs, local ledgers, visible
  terminal output, and output ROOT counts.
- Record in `STATUS_DASHBOARD.md`, `TASK_BOARD.md`, or
  `CODEX_WORK_REGISTER.yaml` that Gmail was unavailable due chat/connector
  auth, so future chats know the missing email evidence is an auth blocker
  rather than a clean inbox.
- After Justin forks/reopens the chat or reconnects Gmail, resume the normal
  Gmail-first flow and mark only consumed pipeline messages read.

## Live Progress Checks

When Justin asks for progress, status, whether a job is okay/stuck, or "what is
going on" for an active RecoilJets/SDCC/Condor job:

1. Check Gmail pipeline messages first.
2. Summarize fresh READY/CHECK/FAILED/held/rescue/profiling signals.
3. Then summarize visible terminal or compact queue/log state.
4. Ask for exactly one compact SDCC diagnostic only if more evidence is needed.

Think like the user following production in real time: cluster/DAG IDs,
idle/running/held counts, worker stdout/err hints, DAGMan final/merge state,
profile rows, output ROOTs appearing, and whether the next stage started.

## Campaign-Scoped Watchdog

The `RecoilJets pipeline watchdog` should be paused by default. It is
campaign-scoped, not a permanent monitor.

Turn it on only with concrete evidence of an active production/training/merge:
user says jobs were submitted, terminal output shows cluster/DAG IDs, a
pipeline email reports running/check/ready/failure, or Justin explicitly asks
Codex to watch.

Before activating, record the campaign in `STATUS_DASHBOARD.md` or
`TASK_BOARD.md`: dataset, command/cluster ID, evidence source, and condition
that turns the watchdog off.

Do not turn it off merely because a stage failed. Diagnose, preserve compact
evidence, clean failed bulky artifacts when safe, patch/relaunch or provide the
exact next command, then stop only when the intended task is completed,
validated, explicitly stopped, or blocked on user approval.

Watchdog output should be action-oriented: failures, held/stale jobs, ready
outputs, or exact next commands. It must not submit, transfer, delete, type
passwords, or mutate SDCC.
