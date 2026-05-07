# Analysis Environment

For any ROOT macro, `root`, `root-config`, ACLiC build, or ROOT-dependent script in this repo, use the analysis environment that matches the user's `use_root_analysis` fish function.

Canonical wrapper:

```bash
./scripts/root_in_analysis_env.sh /Users/patsfan753/Desktop/analysis/env/bin/root -l -q 'macros/MyMacro.C()'
```

Notes:

- This wrapper reproduces the `use_root_analysis` setup from `~/.config/fish/config.fish`.
- It forces the known-good RooUnfold install at `external/RooUnfold` to win over partial env copies.
- Do not rely on interactive shell startup files for ROOT in this repo; use the wrapper explicitly.

# Durable Analysis Context

- Treat this section as the unified workspace memory for stable, analysis-wide
  facts the user states in conversation. Add only high-value facts that will
  help future Codex runs make better technical decisions across chats: dataset
  semantics, event/file scale, SDCC workflow conventions, durable run policies,
  and analysis-output priorities. Do not add transient todo items, one-off run
  status, speculative guesses, or bulky conversational notes.
- RecoilJets analysis input segment/list files are typically about 100k events
  per segment/input ROOT file unless a specific list or production note says
  otherwise. Use this for rough Condor sizing: `groupSize=7` is about 700k
  events/job, and `groupSize=20` is about 2M events/job.
- For exact job counts, prefer the submitter's dry counters on the SDCC
  checkout, for example `./RecoilJets_Condor_submit.sh isAuAu CHECKJOBS
  groupSize 7` or the corresponding SIM dataset. Use those counts before
  recommending a large production submission.
- The user has typically never seen more than about 15k of their own Condor
  jobs running simultaneously, even when more jobs are submitted. Treat that as
  a practical concurrency ceiling for throughput planning: submitted DAG size
  can be larger, but memory requests should be low enough for many jobs to match
  slots quickly within that real simultaneous-running limit.
- For large productions, prefer a staged ramp: upload code, run model/workflow
  checks, submit a small pilot or capped segment/sample, inspect the parseable
  stage email/logs, then scale to the full DAG once the job count and failure
  behavior are understood.
- For AuAu ML training commands on SDCC, use the user's working ML Python
  environment explicitly instead of default `python3`:
  `RJ_ML_PYTHON=/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python`. This env
  has been verified to import `uproot`, `pandas`, `numpy`, `sklearn`,
  `xgboost`, and `ROOT`. The default sPHENIX CVMFS Python may have PyROOT but
  can be missing `uproot` and `xgboost`, causing `trainTightBDT`/`trainMLAll`
  preflight failures.
- The production RecoilJets histogram path is direct RecoilJets fanout, not
  pool/replay. Use `condor all` / `condorDoAll` for optimized production and
  `allDirect` / `condorDoAllDirect` only for one-cfg-per-pass validation.
- When RecoilJets trigger logic, GL1 bits, `ScaledVector`, live/scaled counts,
  or scaledown interpretation needs closer grounding, use the local manual
  `usefulDocs/Gl1-gtm_user_manual_v53.pdf` as the first reference before making
  trigger-semantics claims or code changes.
- When the user asks which SDCC node/submit host to use, or any variation of
  "which node should I submit to?", do not run the helper automatically.
  Instead, tell the user to run `./checkCondorQ` from the local repository base
  and submit from the rank-1 reachable host. The helper SSHes to `sphnxuser01`
  through `sphnxuser08`, checks Condor queue pressure, ranks primarily by
  all-user submit-node busyness rather than the user's own queue, summarizes the
  user's active clusters/configs when present, and does not submit or modify
  jobs. It is local-only and must not be uploaded to SDCC or added to the SFTP
  uploader allowlist. For a known Condor cluster, tell the user to run
  `./checkCondorQ --cluster <cluster_id>` to quickly check whether it is still
  active and which submit host/config/output it belongs to. The helper streams
  only small queue/log snippets from SDCC and stores temporary local copies
  under `/tmp/checkCondorQ.*` during the run.

# RecoilJets Gmail Pipeline Emails

- The Gmail label `RecoilJets Pipeline` is the live email dashboard for
  RecoilJets/SDCC/Condor/DAGMan stage notifications. Treat messages matching
  `RECOILJETS_STAGE_EMAIL_V1`, `[RecoilJets]`, DAGMan, Condor, HTCondor,
  `READY`, `CHECK`, `FAILED`, held, rescue, or removed as pipeline-status
  messages.
- When Codex reads a new RecoilJets pipeline email and uses its information in
  the chat, extract the actionable fields first: dataset, stage, status,
  manifest/output paths, DAG/log/error paths, rescue count, profile summary,
  next action, and any job/cluster IDs. Report the useful status compactly.
- After the message has been consumed and either reported to the user or folded
  into `codex_notes/` when explicitly requested, mark that RecoilJets pipeline
  email read by removing Gmail's `UNREAD` label. This keeps the label useful as
  a dashboard for genuinely new stage changes.
- Do not delete RecoilJets pipeline emails. Do not archive fresh actionable
  `READY`, `CHECK`, or `FAILED` emails unless the user explicitly asks. It is
  fine to label or mark stale already-consumed RecoilJets job emails read when
  the user asks for inbox/label cleanup.
- If an email reports a failure or `CHECK`, do not mark the stage itself
  resolved just because the email was marked read. Marking the email read only
  means Codex consumed the notification; the pipeline state still needs terminal
  logs, report pulls, or user confirmation before being called fixed.

# Input Files

- `InputFiles/` contains local analysis input ROOT files and should be treated as a first-class project directory like `macros/`, `scripts/`, `src/`, and `src_AuAu/`.
- When the user refers to input files or local ROOT inputs, first check `InputFiles/` from the repository base, especially `InputFiles/simEmbedded/`, `InputFiles/simPhotonJet/`, and data-specific subdirectories.
- The `@InputFiles` UI mention may not autocomplete in all clients; still use the literal path `InputFiles/...` when the user refers to those files.
- When the user asks what input variants, cuts, or local merged outputs are
  available for a dataset generated by `scripts/mergeRecoilJets.sh`, run the
  local-only inspector and summarize its output digestibly:

```bash
./checkAvailableVariant <dataset>
```

  Examples: `./checkAvailableVariant isSimEmbedded`,
  `./checkAvailableVariant isSimEmbeddedInclusive`,
  `./checkAvailableVariant isSim`, `./checkAvailableVariant isAuAu`,
  `./checkAvailableVariant isPP`, `./checkAvailableVariant isPPrun25`.
  The script parses local filenames in `InputFiles/...`, reports common cuts,
  available variant axes, and SIM sample coverage. It does not use ROOT, SSH,
  or SFTP.

- Canonical local SIM combinations are built by
  `scripts/sftp_get_recoiljets_outputs.sh` after SIM downloads, or on demand
  with:

```bash
./scripts/sftp_get_recoiljets_outputs.sh mergeLocalSim all
```

  Unless the user explicitly asks to analyze individual SIM slices, use the
  canonical merged products for the requested cfg tag:
  `dataOutput/combinedSimOnly/<cfg>/photonJet5and10and20merged_SIM/` for pp
  photon+jet 5+10+20, and
  `dataOutput/combinedSimOnlyEMBEDDED/<cfg>/photonJet12and20merged_SIM/` for
  embedded Au+Au photon+jet 12+20. Inclusive embedded jet 10+20 and inclusive
  jet5 canonical products are also materialized under the corresponding
  `combinedSimOnly*` cfg folders when those local inputs exist.

# Plot Output

- For generated analysis plots, write PNG files by default. Do not also write
  PDF copies unless the user explicitly asks for PDF output.
- Never draw ROOT statistical boxes on analysis plots unless the user explicitly
  asks for them. Always disable them with `gStyle->SetOptStat(0)` and call
  `SetStats(false)` on frame/data histograms used for plotting, including
  temporary frame histograms.
- Anytime Codex produces analysis plot images, it must view the generated PNGs
  before saying they are finalized or ready to view. Check legends, TLatex
  labels, annotations, and data markers/error bars for overlaps. If any text or
  legend blocks data, rerun the individual plot generation as needed until the
  panel is cleanly organized, data remain viewable, and text is large and
  readable. Do this efficiently and carefully for each plot or group of plots
  the user asks Codex to generate.
- Treat presentation readiness as an active visual QA pass, not a file-exists
  check. For every generated plot or plot family, inspect the PNGs and tune the
  macro coordinates, margins, axis ranges, legends, and TLatex positions until
  the plot clearly shows: what is being plotted, the relevant cuts/cfg/centrality
  or trigger context, the legend entries, and the `sPHENIX Internal` label block.
  Data points and error bars should remain as visible as possible; labels and
  legends must not cover the important distributions unless there is no better
  placement. If a template is reused across many plots, inspect representative
  outputs from each template and rerun after adjustments.
- Match the analysis label style used in `macros/AnalyzeRecoilJets_RunTriggerAna.cpp`:
  write the experiment label as `#it{#bf{sPHENIX}} Internal`, not
  `#bf{sPHENIX Internal}` or plain `sPHENIX Internal`. Use ROOT text font 42
  for the label block unless the existing macro helper already provides a
  stricter local style.
- When the user asks Codex to generate or regenerate plots for slides, Codex
  must first show the generated PNGs directly in the conversation and discuss
  any physics or formatting concerns before editing Google Slides. Iterate on
  the plot images in chat until the user explicitly says the plots are ready for
  slides. Only after that approval should Codex insert or replace the plots in
  the requested Google Slides deck.
- If a generated plot looks physically surprising or inconsistent with the
  user's expectation, pause the slide-update step and debug the plotted inputs,
  histogram normalization/scaling, numerator/denominator definitions, and
  relevant ROOT contents before presenting it as slide-ready.

# Slide Style Preferences

- For Google Slides and other presentation edits, prefer body text at least
  14 pt when possible. Use line spacing around 1.2-1.8 so text fills the slide
  space cleanly instead of looking cramped or stranded.
- For overview and conclusion slides, it is okay to use larger fonts and
  larger spacing to fill the space nicely, as long as the slide remains clear
  and easy to read.

# Change Control

- Before making any file edit, state exactly which file(s) will change and what the change is intended to do.
- Do not make the edit until the user explicitly approves it, unless the user has already directly requested that exact edit.

# Local Tool Dependencies

- If a missing local CLI tool or Python package would materially improve speed,
  accuracy, or reliability, Codex may ask to install it. The request must be
  concise and name the dependency plus the practical reason it helps.
- Do not silently install dependencies. Ask first, and use the approved
  escalation flow when the install needs network access or writes outside the
  workspace.
- Prefer small, standard developer tools when useful, for example `ripgrep`
  for fast code search.

# SSH Diagnostic CLI Handoff

- SDCC terminal default mode:
  When the user is working in a visible terminal in this chat and that terminal
  is SSHed into SDCC/BNL or otherwise operating in the remote SDCC analysis
  checkout, Codex should initially treat that terminal as user-controlled.
  Codex may inspect visible terminal output, scroll/read terminal state when
  available, summarize what happened, track submitted job IDs/DAG paths/
  manifests/output paths, and tell the user what command to run next.
- Open SSH terminal delegation:
  If the user explicitly authorizes Codex to drive the visible SSH terminal for
  the current diagnostic task or current chat, for example "you can drive the
  open SSH terminal for this diagnostic" or "run the needed SSH checks in the
  terminal below", Codex may type and run compact diagnostic commands directly
  in that already-open SSH terminal and then inspect the resulting output. Treat
  that authorization as applying only to the visible SSH terminal and only while
  the task remains related to the user's current question.
- Diagnostic commands that may be run under open SSH terminal delegation include
  read-only or low-risk inspection commands such as `pwd`, `ls`, `find`, `rg`,
  `grep`, `sed`, `awk`, `python3` one-off parsers, `jq`, `wc`, `head`, `tail`,
  `diff`, `comm`, `sort`, `condor_q`, `condor_history`, and project-specific
  status/report inspection commands. Prefer commands that print compact
  summaries over huge raw dumps. If a large output is unavoidable, write it to a
  clear temporary path and print only the key summary plus that path.
- Even under open SSH terminal delegation, Codex must still ask for explicit
  approval before running commands that submit jobs, remove or overwrite files,
  edit remote files, transfer files, change permissions, kill/hold/release
  jobs, run `sftp`/`scp`/`rsync`, install software, expose secrets, or otherwise
  mutate remote SDCC state. Do not ask for, read, store, repeat, or type SDCC
  passwords or other secrets.
- If the visible terminal is not clearly SSHed into the intended SDCC checkout,
  or if the needed command is risky or ambiguous, fall back to showing the
  paste-ready command and asking before running it.
- This SSH delegation policy applies only to visible SDCC/remote SSH terminals.
  For the local Mac checkout/workspace terminal, Codex may continue to run
  commands, edit files, inspect outputs, run tests, and operate normally within
  the active sandbox and approval rules.
- Track SDCC-side workflow state from visible terminal output or user-pasted
  output: which dataset/mode was submitted, exact command, job/DAG ids, manifest
  paths, output bases, transfer commands/results, held/failed reasons, and
  parseable profile lines such as `RECOILJETS_JOB_PROFILE_V1`.
- When Codex needs information that can only be obtained from the user's
  private Brookhaven/SDCC SSH terminal and open SSH terminal delegation has not
  been granted for the current task, prepare one paste-ready shell command. If
  the command is multi-line or otherwise awkward to copy, show it verbatim in
  the conversation, describe concisely what it checks or does, and ask the user
  to approve copying it to the clipboard.
- Only run `pbcopy` after the user confirms approval, for example by replying
  `y` or otherwise clearly approving. If the user's current message explicitly
  asks Codex to copy a specific command to the clipboard, that request itself is
  approval; do not ask for an additional `y`.
- The clipboard content must contain the command only: no prose, no Markdown
  fences, no expected-output notes.
- After copying, tell the user that the command is on their clipboard, where to
  run it, and ask them to paste the terminal output back verbatim.
- Keep the SSH diagnostic handoff continuous: after running or giving a
  diagnostic command, do not treat the task as complete until the visible
  terminal output or pasted output has been inspected and Codex has extracted
  the needed conclusion or next fix. If the command reveals a problem, respond
  with the precise diagnosis, the file/code change or operational next step,
  the relevant SDCC-safe resubmission/retry command for the affected jobs only,
  and any required local `scripts/sftp_push_recoiljets.sh ...` upload command.
  If the output confirms the expected state, say so clearly and record any
  important job IDs, manifests, output paths, or tuning evidence in the project
  notes when it affects future analysis continuity.
- For iterative SDCC debugging where Codex asks the user to run CLI diagnostics
  and paste or show the output, keep the assistant response logically open while
  inspecting the visible terminal output. Do not answer only the first fragment
  of a long pasted diagnostic; continue watching the terminal/readback in the
  same debugging turn until the command has finished or until the visible output
  clearly identifies the next needed command or code fix.
- If `pbcopy` is unavailable or blocked, say that clipboard copy failed and
  provide the command in a fenced shell block instead. Do not create a temp
  text file for this workflow unless the user explicitly asks for one.
- Prefer commands that print compact, diagnostic summaries over huge raw dumps.
  If a large output is unavoidable, have the SSH-side command write a log path
  and print only the key summary plus the path.

# sPHENIX Pipeline Paths

Changes to any of the following files are changes the user will copy into the
SSH analysis checkout and run from:

```bash
/sphenix/u/patsfan753/scratch/thesisAnalysis
```

Treat these as sPHENIX-side pipeline files, not local Mac-only files:

- `/Users/patsfan753/Desktop/ThesisAnalysis/scripts/mergeRecoilJets.sh`
- `/Users/patsfan753/Desktop/ThesisAnalysis/scripts/RecoilJets_Condor_AuAu.sh`
- `/Users/patsfan753/Desktop/ThesisAnalysis/scripts/RecoilJets_Condor_submit.sh`
- `/Users/patsfan753/Desktop/ThesisAnalysis/scripts/RecoilJets_Condor.sh`
- `/Users/patsfan753/Desktop/ThesisAnalysis/src_AuAu/RecoilJets_AuAu.h`
- `/Users/patsfan753/Desktop/ThesisAnalysis/src_AuAu/RecoilJets_AuAu.cc`
- `/Users/patsfan753/Desktop/ThesisAnalysis/src/PhotonClusterBuilder.cc`
- `/Users/patsfan753/Desktop/ThesisAnalysis/src/PhotonClusterBuilder.h`
- `/Users/patsfan753/Desktop/ThesisAnalysis/src/RecoilJets.cc`
- `/Users/patsfan753/Desktop/ThesisAnalysis/src/RecoilJets.h`
- `/Users/patsfan753/Desktop/ThesisAnalysis/macros/analysis_config.yaml`
- `/Users/patsfan753/Desktop/ThesisAnalysis/macros/Fun4All_recoilJets_AuAu.C`
- `/Users/patsfan753/Desktop/ThesisAnalysis/macros/Fun4All_recoilJets_unified_impl.C`
- `/Users/patsfan753/Desktop/ThesisAnalysis/macros/Fun4All_recoilJets.C`

Only offline macros that the user runs after the online pipeline completes
should be treated as local-only code.

# SDCC Output Storage Policy

- Large intermediate pipeline products belong on the SDCC bulk `/sphenix/tg/tg01/bulk/...`
  filesystem, matching the original production pipeline. This includes per-run
  folders, per-segment/per-chunk ROOT outputs, Condor production trees, and
  other high-volume artifacts created before final merge/analysis consumption.
- Final merged outputs, local analysis-ready ROOT files, pulled products, plots,
  slide inputs, and other user-facing final artifacts belong in the user's
  normal local/basis analysis directories, such as `InputFiles/...`,
  `dataOutput/...`, plot output folders, or the relevant Google Drive thesis
  workspace once explicitly requested.
- When changing pipeline paths, preserve this split: bulk for scalable
  production intermediates, user/local analysis areas for final products and
  presentation-ready outputs. Do not redirect large per-segment outputs
  into the user's basis/home directory unless the user explicitly asks.

# SDCC Transfer Workflow

Use `scripts/sftp_push_recoiljets.sh` as the only supported transfer helper for
copying known SSH pipeline files from the local Mac checkout to the SDCC
analysis checkout.

Use `scripts/sftp_get_recoiljets_outputs.sh` as the local-only helper for
pulling merged ROOT outputs from SDCC into `InputFiles/...`. It opens
interactive `sftp`, lists exact remote final ROOT files for the requested
dataset, prints the local overwrite preview, and then downloads only those
listed files after user confirmation. Do not SSH directly, do not run raw
`sftp` manually from Codex, and do not ask for or handle the SDCC password.

Codex may run these two local transfer helpers from the local Mac workspace
after explicitly asking the user for approval for the exact command:

```bash
./scripts/sftp_push_recoiljets.sh ...
./scripts/sftp_get_recoiljets_outputs.sh ...
```

If the helper opens an interactive password prompt, the user enters the SDCC
password manually. Codex must not ask for, read, store, repeat, or type the
password. If the user does not approve Codex running the helper, provide the
same command for the user to run manually.

Read-only SDCC/local transfer checks:

- `scripts/sftp_push_recoiljets.sh` also supports read-only `status` and `diff`
  modes. These fetch mapped SDCC files into a local temp directory with SFTP and
  compare them to the local checkout. They must not modify SDCC files.
- When the user asks whether SDCC is "up to date", "in sync", "the same as
  local", "all transferred", "all transferable files match", or similar, steer
  them to:

```bash
./scripts/sftp_push_recoiljets.sh status
```

- When the user asks whether the currently changed transferable files are on
  SDCC, or whether local work-in-progress has been uploaded, steer them to:

```bash
./scripts/sftp_push_recoiljets.sh status changed
```

- When the user asks whether the Condor/submission pipeline scripts match SDCC,
  steer them to:

```bash
./scripts/sftp_push_recoiljets.sh status condor
```

- When the user asks to see what differs for a specific mapped file, steer them
  to:

```bash
./scripts/sftp_push_recoiljets.sh diff <file-or-basename>
```

- When the user asks to see line-by-line differences for all changed
  transferable files, steer them to:

```bash
./scripts/sftp_push_recoiljets.sh diff changed
```

- Treat `MATCH` as local and SDCC byte-identical for that mapped file,
  `DIFFER` as local and SDCC both existing but not identical, and `MISSING` as
  absent on SDCC at the mapped remote path. If all entries are `MATCH` except a
  known file the user is actively editing, explain that the rest of the mapped
  transfer set is consistent and name the exception.
- These read-only check modes still open interactive SFTP and may prompt for the
  SDCC password. Codex may run them locally only after explicit user approval of
  the exact command; otherwise give the command for the user to run manually.

Local repo base:

```bash
/Users/patsfan753/Desktop/ThesisAnalysis
```

Remote SDCC repo base:

```bash
/sphenix/u/patsfan753/scratch/thesisAnalysis
```

Remote host:

```bash
patsfan753@sftp.sdcc.bnl.gov
```

Rules for future Codex runs:

- Treat the mapped files below as SSH pipeline files.
- When Codex edits any mapped file, its final response must include the exact
  `./scripts/sftp_push_recoiljets.sh ...` command the user should run.
  This applies even for tiny warning fixes, one-line changes, or follow-up
  edits to files changed earlier in the conversation.
- For substantive SDCC pipeline edits, assume the user usually wants the local
  GitHub repo updated after the successful SDCC transfer. Include the uploader's
  built-in commit/push form by default, with the exact selected files/groups and
  a descriptive message:

```bash
./scripts/sftp_push_recoiljets.sh <changed-files-or-groups> --commit-push -m "Concise change summary"
```

- If the user explicitly asks for upload-only, omit `--commit-push`.
- The commit message should describe the actual code/config change, not just say
  "changes".
- The uploader commits only the selected transferable files after successful
  SFTP. It must never use `git add .` for this workflow.
- Do not SSH directly.
- Do not run `sftp` directly.
- Do not ask for or handle the SDCC password.
- The local uploader/getter helpers may be run by Codex only after explicit
  user approval of the exact command. The user still enters any SDCC password
  prompt manually.
- When visible SDCC terminal output or parseable stage emails show that merged
  outputs are ready, Codex may offer to run
  `scripts/sftp_get_recoiljets_outputs.sh` locally to pull the relevant dataset
  into `InputFiles/...`; ask for approval before running it.
- Prefer the smallest upload command that covers the changed files.
- Prefer explicit changed file names for source files. Use `RecoilJets.cc`
  and `RecoilJets.h` directly for `src/`; use `RecoilJets_AuAu.cc` and
  `RecoilJets_AuAu.h` directly for `src_AuAu/`. Do not use `pp`, `auau`,
  or `both` source groups.
- If one known file changed, give a one-file command, for example:

```bash
./scripts/sftp_push_recoiljets.sh mergeRecoilJets.sh
```

- If several known files changed, give one command with all file names, for
  example:

```bash
./scripts/sftp_push_recoiljets.sh mergeRecoilJets.sh analysis_config.yaml RecoilJets_AuAu.cc
```

- If a whole non-source logical group changed, use the group, for example:

```bash
./scripts/sftp_push_recoiljets.sh condor
```

- If Codex creates or edits a new file that is needed on SDCC but is not in the
  uploader allowlist, explicitly say that the uploader allowlist needs to be
  extended. Ask before editing the uploader/AGENTS.md unless the user already
  requested that exact update.
- Do not use wildcard uploads for `macros/`, `scripts/`, `src/`, or `src_AuAu/`.
- New SSH-needed macro files must be added explicitly to the uploader allowlist
  with their intended remote path.
- `scripts/sftp_push_recoiljets.sh` is local-only. Do not include it in SDCC
  transfer commands.
- `scripts/sftp_get_recoiljets_outputs.sh` is local-only. Do not include it in
  SDCC transfer commands.
- Because the uploader/getter are local-only, do not rely on `--commit-push` to
  commit edits to those helper scripts or to `AGENTS.md`. Use normal local git
  commands for those files if needed.

Known transfer map:

```text
scripts/audit_auau_grl_projection.sh      -> scripts/audit_auau_grl_projection.sh
scripts/audit_auau_ml_training_smoke.py   -> scripts/audit_auau_ml_training_smoke.py
scripts/estimateEmbeddedPhotonXsec.sh     -> scripts/estimateEmbeddedPhotonXsec.sh
scripts/make_dstListsData.sh              -> scripts/make_dstListsData.sh
scripts/makeThesisSimLists.sh             -> scripts/makeThesisSimLists.sh
scripts/mergeRecoilJets.sh                -> scripts/mergeRecoilJets.sh
scripts/root_in_analysis_env.sh           -> scripts/root_in_analysis_env.sh
scripts/RecoilJets_Condor_AuAu.sh         -> RecoilJets_Condor_AuAu.sh
scripts/RecoilJets_Condor_submit.sh       -> RecoilJets_Condor_submit.sh
scripts/RecoilJets_Condor.sh              -> RecoilJets_Condor.sh
scripts/train_auau_jet_residual_bdt.py    -> scripts/train_auau_jet_residual_bdt.py
scripts/train_auau_photon_bdt.py          -> scripts/train_auau_photon_bdt.py
macros/analysis_config.yaml               -> macros/analysis_config.yaml
macros/Fun4All_recoilJets.C               -> macros/Fun4All_recoilJets.C
macros/Fun4All_recoilJets_AuAu.C          -> macros/Fun4All_recoilJets_AuAu.C
macros/Fun4All_recoilJets_unified_impl.C  -> macros/Fun4All_recoilJets_unified_impl.C
macros/PrintPPStitchDiagnostics.C         -> macros/PrintPPStitchDiagnostics.C
src/RecoilJets.cc                         -> src/RecoilJets.cc
src/RecoilJets.h                          -> src/RecoilJets.h
src_AuAu/RecoilJets_AuAu.cc               -> src_AuAu/RecoilJets_AuAu.cc
src_AuAu/RecoilJets_AuAu.h                -> src_AuAu/RecoilJets_AuAu.h
```

Uploader groups:

```text
condor    scripts/RecoilJets_Condor.sh + scripts/RecoilJets_Condor_AuAu.sh + scripts/RecoilJets_Condor_submit.sh
macros    known pipeline macros/configs
scripts   known pipeline helper scripts
pipeline  all known transferable pipeline files except local-only sftp helpers
```

Important mapping detail: most files preserve their relative path on the remote
side, but these three local files are under `scripts/` locally and live at the
remote repo base on SDCC:

```text
scripts/RecoilJets_Condor.sh       -> RecoilJets_Condor.sh
scripts/RecoilJets_Condor_AuAu.sh  -> RecoilJets_Condor_AuAu.sh
scripts/RecoilJets_Condor_submit.sh -> RecoilJets_Condor_submit.sh
```

# Durable Project State Tracking

`codex_notes/` is the local source of truth for running analysis state that
should survive long chats, context compaction, and future Codex sessions. Treat
it as local-only project state; it is not SDCC pipeline code.

Before analysis/output work, plot review, dataset validity claims, or Condor
status triage, read:

- `codex_notes/PROJECT_BOARD.md`
- `codex_notes/DATASET_STATUS.md`
- `codex_notes/KNOWN_ISSUES.md`

Use `codex_notes/RUN_LOG.md` for short dated updates about user-reported job
submissions, SDCC pasted summaries, local pulls, merges, plot-generation passes,
and important local inspection commands.

Recording rules:

- If the user says "track this", "record this", "note this", or otherwise
  explicitly asks Codex to remember a status update, update the appropriate
  `codex_notes/` files.
- If the user casually mentions a new task, issue, job submission, finished job,
  stale dataset, or interesting finding, ask before recording it.
- Keep notes concise, dated, and evidence-based. Do not mark outputs valid or
  stale without evidence such as local file timestamps, job IDs, SDCC pasted
  terminal output, local ROOT/input inspection, or a clear user statement.
- Do not silently add subjective "interesting findings"; ask first unless the
  user explicitly requested tracking.
- SDCC/Condor status must come from user-pasted terminal output or the approved
  clipboard handoff workflow. Do not SSH directly, run raw `sftp`, or ask for
  the SDCC password.

Current high-priority tracking context:

- The ABCD/purity tight-axis fix is present in local git commit `507eabd3`
  (`2026-05-05 21:17 EDT`) in `src/RecoilJets.cc` and
  `src_AuAu/RecoilJets_AuAu.cc`. It changes `fillIsoSSTagCounters(...)` to use
  the caller-computed `TightTag` from the active configured photon-ID variants.
- Do not consider affected output families production-valid until the source is
  transferred to SDCC, affected jobs rerun, fresh ROOT outputs pulled, and
  `codex_notes/DATASET_STATUS.md` updated with evidence.
- Existing `InputFiles/` ROOT products using a non-reference tight or
  non-tight photon-ID axis, including `tightVariantA`, `nonTightVariantA`,
  renamed/new `newPPG12`-style settings, or later AuAu BDT sideband settings,
  should be treated as stale only for the affected ABCD/purity-derived
  histogram families if produced before that local fix. The ROOT files are not
  globally unusable.
- Outputs that only changed preselection while keeping
  `tightReference/nonTightReference` are lower risk for this specific
  tight-axis mismatch; use normal provenance checks, but do not mark them stale
  for this bug without additional evidence.
- Affected histogram families for this bug are:
  - candidate ABCD counts: `h_isIsolated_isTight*`,
    `h_notIsolated_isTight*`, `h_isIsolated_notTight*`,
    `h_notIsolated_notTight*`;
  - candidate ABCD distributions: `h_Eiso_ABCD_[A-D]*` and
    `h_pTgamma_ABCD_[A-D]*`;
  - ABCD-tagged shower-shape templates:
    `h_ss_<var>_isIsolated_isTight*`,
    `h_ss_<var>_notIsolated_isTight*`,
    `h_ss_<var>_isIsolated_notTight*`,
    `h_ss_<var>_notIsolated_notTight*`;
  - truth-signal leakage: `h_sigABCD_MC*`;
  - event-leading xJ purity counts: `h_xJpurityLead_isIsolated_isTight*`,
    `h_xJpurityLead_notIsolated_isTight*`,
    `h_xJpurityLead_isIsolated_notTight*`,
    `h_xJpurityLead_notIsolated_notTight*`;
  - event-leading xJ leakage: `h_xJpurityLead_sigABCD_MC*`.
- Non-ABCD histograms, trees, event QA, jet QA, trigger QA, and other outputs
  can still be used if the requested plot/result does not read or derive from
  the affected families above. When discussing a plot, first check whether its
  macro or inputs consume those families before warning that the plot is stale.
- When fresh rerun outputs are pulled into `InputFiles/`, update
  `codex_notes/DATASET_STATUS.md` with timestamp/job evidence and retire
  obsolete stale warnings instead of leaving old global cautions in place.
- `tightReference/nonTightReference` outputs are lower risk for this specific
  bug because the old hard-coded ABCD tight axis matched reference behavior, but
  still require normal provenance checks before being used in final results.
