# SDCC Operations

## Local-Only Access Details

Specific SDCC hostnames, usernames, remote checkout paths, SSH command
patterns, Condor cluster IDs, and active watchdog status belong in local-only
memory, not GitHub-tracked files.

Before any SDCC or Condor work in Justin's local checkout, load
`agent_context/local/SDCC_LOCAL_RULES.md` if it exists. If it is absent, ask
Justin for the current local access rule instead of inventing one or committing
site-specific access details.

Use read-only diagnostics by default. Submissions, merge reruns, cleanup, file
transfer, remote edits, and job control require explicit current-task approval.
Before a complex check, run a tiny probe on the target node. If that fails,
stop retrying quote-heavy SSH commands and fall back to the visible terminal or
a paste-ready handoff.

## Persistent Session Rule

For more than a tiny check, open one persistent SSH session and reuse it for
`condor_q`, `tmux`, `tail`, `find`, `python3`, and report inspection. Avoid
many quote-heavy one-shot SSH commands. Exit the remote shell when finished.

## Visible Terminal Handoff

If Justin is working in a visible SDCC terminal, initially treat it as
user-controlled. Codex may read visible output, summarize, and provide compact
commands. If Justin explicitly delegates the terminal for the current
diagnostic, Codex may run compact diagnostic commands there.

That delegation applies only to the visible SSH terminal and only while the
task remains related to the current question. If the visible terminal is not
clearly SSHed into the intended SDCC checkout, stop and fall back to a
paste-ready command.

Allowed diagnostics include `pwd`, `ls`, `find`, `rg`, `grep`, `sed`, `awk`,
`python3` one-offs, `jq`, `wc`, `head`, `tail`, `diff`, `sort`, `condor_q`,
`condor_history`, and project-specific read-only report commands.

Prefer compact summaries over huge raw dumps. If large output is unavoidable,
write it to a clear temporary path on SDCC and print only the key summary plus
that path.

Even under delegation, ask before commands that submit jobs, remove or
overwrite files, transfer files, edit remote files, change permissions,
kill/hold/release jobs, run `sftp`/`scp`/`rsync`, install software, or expose
secrets.

## Paste-Ready Diagnostic Handoff

When live SDCC evidence is needed but Codex does not have explicit delegation
to operate a terminal, give Justin one compact paste-ready command block rather
than several fragile one-liners.

The block should:

- print host, working directory, timestamp, and relevant campaign tag;
- read only small logs, queue summaries, reports, manifests, and file counts;
- avoid huge ROOT reads/transfers unless Justin asked for them;
- be safe to rerun and not mutate remote state.

If putting the command on the macOS clipboard would help, ask first and then
copy only the command text. Never copy secrets, never type passwords, and do
not use clipboard automation as a substitute for explaining what will run.

Only use `pbcopy` after Justin clearly approves, unless his current message
explicitly asks Codex to copy a specific command. Clipboard content must be the
command only: no prose, no Markdown fences, no expected-output notes. After
copying, tell Justin where to run it and ask for the terminal output verbatim.
If `pbcopy` is unavailable or blocked, provide the command in a fenced shell
block instead.

Keep the handoff continuous. After giving or running an SDCC diagnostic, inspect
the visible or pasted output before calling the task complete. Extract the
conclusion, next fix, affected campaign tag/output root/cluster IDs, and any
needed upload or retry command.

When checking a campaign, keep lanes separated: identify submit host, campaign
tag, YAML/config, dataset, output root, and cluster IDs before drawing
conclusions. Do not assign errors from one terminal/node to another campaign
without matching tag, output root, or cluster ID.

Track SDCC-side workflow state from visible or pasted output, including
dataset/mode, exact command, DAG/job IDs, manifest paths, output bases,
transfer commands/results, held/failed reasons, and parseable profile lines
such as `RECOILJETS_JOB_PROFILE_V1`.

## Submission Discipline

Before any new training/production/merge, run the duplicate guard. Then:

- explicitly pass Codex provenance into the remote submit environment for any
  RecoilJets/Condor workflow that can send pipeline emails. The remote command
  that creates the DAG must set both `RJ_CODEX_CHAT_NAME=<short chat label>` and
  `RJ_CODEX_THREAD_ID=<current Codex thread id>`; local `CODEX_*` environment
  variables do not automatically propagate through nested SSH or later Condor
  scheduler nodes. If either value is unavailable, stop before submission and
  choose a short deterministic `RJ_CODEX_CHAT_NAME` from the workstream/tag;
  never accept `codex_chat_name=unknown` or `codex_thread_id=unknown` as the
  intended state for new Codex-launched workflows.
- record dataset/mode, command, submit host, timestamped roots, DAG paths,
  cluster IDs, report/output roots, and next checkpoint in `STATUS_DASHBOARD.md`
  or `TASK_BOARD.md`;
- for targeted diagnostics or slide-driven validation, run the narrowest
  existing submit/merge path that can answer the question. Reuse the validated
  infrastructure, but constrain its matrix explicitly with row selectors,
  config matches, sample gates, `maxJobs`, or fanout caps so the job scope
  matches the requested deliverable. Before submitting, print or otherwise
  verify the intended cfg count, cfg tag(s), sample list, group size, memory
  request, and total job count. If the command would launch unrelated cfg rows,
  broad fanout, extra model families, or a production-size matrix, stop and
  narrow it instead of relying on the broad default;
- for paired or multi-dataset submissions, prefer one pasteable driver block
  with `set -euo pipefail`, preflight queue/stale-job checks, explicit memory
  knobs, sequential submission, timestamped log, and post-submit `condor_q`;
- for matrix-style submissions, use a queue-gated driver, submit bounded chunks,
  gate on matching campaign tag, count only active jobs, and stop on held jobs;
- use fresh campaign tags after canceled or partially submitted attempts.

For large RecoilJets target-WP or BDT/ML matrices, default to explicit queue
knobs like:

```bash
RJ_TARGETWP_QUEUE_GATE=1
RJ_TARGETWP_QUEUE_SCOPE=matching
RJ_TARGETWP_MAX_QUEUED=40000
RJ_TARGETWP_RESUME_BELOW=0
RJ_TARGETWP_POLL_SECONDS=300
```

Before telling Justin to rerun a gated campaign, verify when possible that the
old matching campaign has `0` idle, `0` running, and `0` held active jobs.
Condor `X`/removed and completed records are dead history, not active pressure.

For suggested submissions spanning multiple related commands, provide one
paste-ready driver block or temporary `/tmp/...` script plus one `bash`
invocation. The driver should fail fast, log output, print queue/report paths,
and leave enough evidence for future status checks without reconstructing from
scrollback.

The driver block should define provenance once near the top and reuse it for
all RecoilJets submit/merge commands, for example:

```bash
RJ_CODEX_CHAT_NAME="<short-workstream-label>"
RJ_CODEX_THREAD_ID="<codex-thread-id>"
export RJ_CODEX_CHAT_NAME RJ_CODEX_THREAD_ID
```

For nested SSH commands, place those assignments inside the quoted remote
command that actually runs `RecoilJets_Condor_submit.sh` or `mergeRecoilJets.sh`.
Setting them only in the local macOS shell is insufficient.

## Condor Practical Limits

- Input segment/list files are usually about 100k events per ROOT segment.
  Rough sizing: `groupSize=7` is about 700k events/job, `groupSize=20` is about
  2M events/job.
- For exact job counts, prefer submitter dry counters, for example
  `./RecoilJets_Condor_submit.sh isAuAu CHECKJOBS groupSize 7`.
- Justin has typically never seen more than about 15k of his own Condor jobs
  running simultaneously. Treat that as a practical concurrency ceiling.
- For large productions, stage the ramp: upload code, run model/workflow checks,
  submit small pilot/capped sample, inspect parseable stage email/logs, then
  scale once job count and failure behavior are understood.

## AuAu ML Environment

For AuAu ML training on SDCC, use the local-only `RJ_ML_PYTHON` path recorded
in `agent_context/local/SDCC_LOCAL_RULES.md` or current campaign status. The
intended env has `uproot`, `pandas`, `numpy`, `sklearn`, `xgboost`, and `ROOT`.
Default CVMFS Python may miss `uproot`/`xgboost`.

## Submit Host Helper

When Justin asks which SDCC node/submit host to use, Codex should run the
helper locally from the repository base as a read-only diagnostic and summarize
the rank-1 recommended host:

```bash
./checkCondorQ
```

For a known cluster, Codex should run:

```bash
./checkCondorQ --cluster <cluster_id>
```

This does not require separate confirmation when the current request is about
choosing or checking the submit host. If sandbox/network policy requires
escalation, request the narrow command approval and then run it.

The helper is still local-only: do not upload it to SDCC, do not add it to the
SFTP allowlist, and do not treat it as a production or mutation command. It
streams only small queue/log snippets and should be summarized back to Justin
with the recommended node and any obvious active-cluster context. The exact
node list and access details belong in local-only helper/config state, not in
tracked Codex policy text.

## ROOT Merge Strategy

For large ROOT merge/stitch workloads, prefer layered Condor merging over one
huge `hadd` or one huge recursive merge. Example: for about `10k` files, merge
groups of about `200`, then those intermediate outputs in groups of about
`10`, then a final small merge.

## Production Path

The production RecoilJets histogram path is direct RecoilJets fanout, not
pool/replay. Use `condor all` / `condorDoAll` for optimized production and
`allDirect` / `condorDoAllDirect` only for one-cfg-per-pass validation.
