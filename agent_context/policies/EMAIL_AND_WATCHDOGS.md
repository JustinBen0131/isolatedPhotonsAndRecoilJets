# RecoilJets Gmail And Watchdogs

## Gmail Pipeline Label

Gmail label `RecoilJets Pipeline` is the live dashboard for
RecoilJets/SDCC/Condor/DAGMan notifications.

Treat messages matching `RECOILJETS_STAGE_EMAIL_V1`, `[RecoilJets]`, DAGMan,
Condor, HTCondor, `READY`, `CHECK`, `FAILED`, held, rescue, or removed as
pipeline status.

When reading a new pipeline email, extract:

- dataset, stage, status;
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
