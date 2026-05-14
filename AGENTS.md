# Codex SDCC SSH Rule

- For Codex-run SDCC shell checks/submissions from the local Mac, do **not**
  rely on `ProxyJump` to `sphnxuserXX`: this environment can fail with
  `channel 0: open failed: administratively prohibited: open failed` /
  `stdio forwarding failed`.
- Preferred Codex SSH pattern is nested SSH through the login node, always with
  the macOS launchctl SSH agent:

```bash
SSH_AUTH_SOCK="$(launchctl getenv SSH_AUTH_SOCK)" \
ssh patsfan753@ssh.sdcc.bnl.gov \
  "ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null sphnxuser05.sdcc.bnl.gov 'cd /sphenix/u/patsfan753/scratch/thesisAnalysis && <command>'"
```

- Use this for read-only diagnostics by default. For submissions/merge reruns,
  only run it when the user explicitly authorizes the specific campaign action.
- For sustained SDCC work in one turn, do **not** repeatedly open fresh SSH
  logins for every small check. Open one persistent SSH session and reuse it,
  or use the already-visible SSH terminal when the user explicitly delegates
  that terminal for the current task. Repeated one-shot SSH commands are only
  acceptable for isolated, quick checks.
- File transfer still prefers the project SFTP helpers with the same
  `SSH_AUTH_SOCK="$(launchctl getenv SSH_AUTH_SOCK)"` prefix.

# Scientific Mission

- The ultimate analysis target is PPG19 Au+Au gamma-jet / `x_{J#gamma}`. The
  user is a co-spokesperson of PPG19, and agent decisions should keep that
  analysis as the north star.
- PPG18 pp gamma-jet is the validated reference baseline for the Au+Au work.
  Preserve pp comparisons whenever practical, and use pp behavior to anchor
  Au+Au extensions before treating new Au+Au machinery as final.
- PPG12 is the photon-ID, isolation, BDT/NPB, purity, and pp infrastructure
  reference. Use the existing local checkout at
  `/Users/patsfan753/Desktop/ThesisAnalysis/ppg12codeGit` before inventing new
  photon-ID, BDT, ABCD, stitching, or systematic logic.
- The final gamma-jet analysis should mimic the ATLAS gamma-jet analysis note
  as closely as possible while staying consistent with sPHENIX and PPG
  constraints.

# Reference Hierarchy

- Implementation reference for pp photon-ID infrastructure:
  `ppg12codeGit/`.
- PPG12 box-cut photon-ID reference:
  `usefulDocs/PPG12_analysis_note.pdf`.
- PPG12 BDT/NPB infrastructure reference:
  `usefulDocs/PPG12_analysis_note_withBDT.pdf`.
- ATLAS gamma-jet target-analysis reference:
  `usefulDocs/Gamma_Jet_Analysis_Note (1).pdf`.
- Au+Au xJ/dijet correction and presentation-style reference:
  `usefulDocs/PPG_08_dijet_xJ_in_Au_Au___draft_Conference_note (7).pdf`.
- Trigger semantics reference:
  `usefulDocs/Gl1-gtm_user_manual_v53.pdf`.
- Local agent reference map:
  `agent_context/REFERENCE_MAP.md`.

# Analysis Decision Rules

- Do not invent photon-ID, purity, unfolding, xJ, BDT/NPB, sample stitching, or
  plotting conventions when an existing PPG12/PPG18/ATLAS/PPG08 precedent can
  be reused or adapted.
- Default reasoning order for physics changes: validate pp baseline behavior,
  extend to Au+U/AuAu or embedded context, validate with embedded signal and
  embedded inclusive/background samples, then promote to final data products.
- For corrections, final plots, purity, ABCD, unfolding, normalization,
  stitching, response matrices, or surprising outputs, enter skeptic mode:
  actively check for mathematical mistakes, stale provenance, histogram-family
  mismatches, data/MC definition mismatches, and pipeline-stage mismatches.
- For routine implementation and operational fixes, stay surgical and follow
  existing architecture.
- Prefer evidence over memory. Mark outputs valid or stale only with concrete
  evidence such as ROOT inspection, file timestamps, job IDs, terminal output,
  emails, or a clear user statement.

# Plot And Architecture Reuse

- Reuse existing plotting infrastructure unless it is genuinely impossible.
  Quick plots should mutate/adapt existing macros, helpers, styles, axis-label
  conventions, and output organization instead of creating isolated one-off
  styles.
- Prefer `macros/AnalyzeRecoilJets*`, `macros/sPhenixStyle.*`, and existing
  plot helper patterns before writing new plotting code.
- Final or slide plots must use consistent sPHENIX styling, existing
  histogram naming/normalization conventions, and visual QA before being called
  ready.
- If a plot looks physically surprising, debug definitions, normalization,
  numerator/denominator choices, sample weights, and ROOT input provenance
  before presenting it as a result.

# Slide Integration Rules

- Google Slides Computer Use safety contract:
  - Do not use Computer Use on the user's browser or desktop for Google Slides
    unless it is genuinely useful for a visual/layout check or a narrowly
    targeted edit that cannot be done cleanly through the Google Drive/Slides
    API alone.
  - Before every Computer Use session, first tell the user exactly why Computer
    Use is needed, which Google Slides tab/deck will be touched, and what
    actions will be performed. Wait for explicit permission each time.
  - Computer Use is limited to Google Slides tabs for this project unless the
    user explicitly authorizes another app or tab in that moment.
  - Prefer read-only visual inspection when possible. If editing is needed,
    keep actions narrow, reversible, and targeted to the named slide/deck.
  - Do not switch to unrelated browser tabs, interact with email/calendar/chat,
    submit forms, present, share, delete, or change permissions through
    Computer Use unless the user explicitly asks for that specific action.
  - If the visible browser state is not the expected Google Slides deck/tab,
    stop and ask before navigating or clicking.
- Before substantial Google Slides edits, inspect relevant past Drive
  presentations when available. Use them to learn the user's existing slide
  style and organization for the specific topic.
- Gemini / Beautify Slides may be used as an optional design-critique and
  layout-ideation assistant for difficult slides, but do not let it become the
  source of truth for physics content or final layout. Use it to generate
  cleaner wording, grouping, hierarchy, and color ideas; then implement the
  accepted design with controlled Google Slides edits and verify with
  thumbnails. Never let it delete, overwrite, or reuse the Backup slide.
- Never remove, overwrite, or repurpose the deck's Backup slide. Backup is the
  divider for backup plots/material that go after it. Add new presentation
  slides before the Backup slide, and add backup material after it. Do not
  delete the Backup slide to insert, move, or replace any slide under any
  circumstance.
- Use `agent_context/SLIDE_STYLE_MAP.md` as the durable local guide to the
  user's preferred slide design language. It records examples from strong past
  decks, including teaching slides, colored comparison panels, subtle shadows,
  plot organization, font hierarchy, and BDT/ML translation rules.
- Match the user's preferred slide language: concise bullets above plots, clear
  hierarchy, bolded lead phrases, organized sub-bullets when helpful, body text
  at least 14 pt when possible, line spacing around 1.2-1.8, and a
  human-readable claim/evidence/implication flow.
- Preserve the essence of the user's original wording while making it easier
  to scan and explain aloud.
- Reference similar older slides when possible rather than inventing a new
  visual/text style from scratch.

# Collaboration Model

- Codex owns orchestration: SDCC status strategy, Gmail pipeline emails, local
  notes, dashboards, Google Slides, Drive presentation context, integration,
  transfer-command guidance, and final sanity checks.
- Claude owns bounded code/reference work: PPG12/reference inspection, local
  code changes, local static analysis, feature implementation, and compact
  handoff.
- Split work only across disjoint lanes by default. Avoid two agents editing
  the same files or solving the same unknown.
- Handoffs should be compact and structured: changed files, evidence/tests,
  physics assumptions, risks/open questions, and exact next command, including
  SDCC upload/check command when relevant.

# Local Agent Memory

- `agent_context/TASK_BOARD.md` is the clean milestone/blocker task board.
- `agent_context/STATUS_DASHBOARD.md` is the traffic-light dataset/job/output
  dashboard.
- `agent_context/REFERENCE_MAP.md` is the compact reference decision map.
- Track durable memory only for decisions and evidence: dataset choices, sample
  pairs, job IDs, commands that worked, output validity, physics conclusions,
  and recurring mistakes.
- If the user explicitly says "track", "add", "remember", "record", or similar,
  update local task/memory files.
- If the user casually mentions a possible task, ask before adding it.
- When a task appears finished, move it to `Done Pending Removal` and ask the
  user before removing or archiving it.
- Prompt cleanup at the end of major tasks, production passes, or plot sets.
  Never silently delete useful context.
- For SDCC tests or sidecar workflows, record active run roots and cleanup
  decisions while the work is live. When a test is superseded, failed, or no
  longer needed, remind the user which timestamped artifacts can be cleaned
  after confirming no live jobs still reference them.
- When Codex gives the user an SDCC submission command or sees the user submit
  one in the visible terminal, immediately note the submit host/node, exact
  command, dataset/mode, timestamped run roots, DAG path, Condor cluster IDs,
  report/output roots, and next expected checkpoint in
  `agent_context/STATUS_DASHBOARD.md` or `agent_context/TASK_BOARD.md`. If the
  terminal prompt shows a node such as `sphnxuser07`, use that as evidence.
  Do this for smoke tests, sidecar ML jobs, validation DAGs, and full
  production submissions so future chats can resume without reconstructing the
  campaign from terminal scrollback.
- For multi-dataset SDCC submissions that are meant to run together, especially
  paired `isSimEmbedded`/`isSimEmbeddedInclusive`, pp/SIM sample pairs, or
  signal/background validation campaigns, do not hand the user several fragile
  one-line commands by default. Provide one pasteable driver block or temporary
  shell script that:
  - sets `set -euo pipefail`;
  - prints the submit host, timestamp, config YAML, output roots, memory knobs,
    and dataset list before submitting;
  - runs a compact preflight queue/stale-job check;
  - uses explicit worker and merge memory variables instead of relying on
    ambiguous defaults;
  - uses isolated output roots and guarded cleanup flags only when the user has
    approved a clean resubmit;
  - submits the datasets sequentially so a first failure stops the second
    submission;
  - writes a timestamped local submit log with `tee`;
  - prints a post-submit `condor_q -nobatch` summary.
  After the user runs the block, inspect the terminal output, record the new
  DAG/cluster IDs immediately, and keep watching until the submit state is
  unambiguous.
- For large matrix-style RecoilJets submissions, such as many BDT/ML variants,
  many generated YAML configs, or paired signal/background MC campaigns, do not
  submit the whole matrix into Condor at once. Use a queue-gated driver by
  default:
  - compute the expected job count before submission and warn if the full matrix
    would exceed about `40000` queued jobs;
  - submit one bounded chunk at a time, usually one YAML config or one dataset
    pair, then wait before submitting the next chunk;
  - gate on the matching campaign tag, not the whole user queue, unless the user
    explicitly asks for a global user-queue gate;
  - count only active jobs when gating: ignore Condor removed/completed records
    such as `X`/removed and completed states, and stop immediately if any held
    jobs appear;
  - prefer explicit knobs like `RJ_TARGETWP_QUEUE_GATE=1`,
    `RJ_TARGETWP_QUEUE_SCOPE=matching`, `RJ_TARGETWP_MAX_QUEUED=40000`,
    `RJ_TARGETWP_RESUME_BELOW=0`, and `RJ_TARGETWP_POLL_SECONDS=300`;
  - use a fresh campaign tag after any canceled or partially submitted attempt
    so partial raw outputs cannot be confused with the clean gated campaign.
  Before telling the user to rerun a gated campaign, verify with an SSH-auth
  `condor_q` check when possible that the old matching campaign has `0` idle,
  `0` running, and `0` held jobs. If Condor still shows removed `X` rows, treat
  them as dead records, not active queue pressure.
- Keep an `Automation Backlog` of repeated manual checks that could become
  read-only local scripts, dry-run counters, QA summaries, status dashboards,
  or plot-inspection helpers. Ask before implementing new automation.

# Campaign-Scoped Watchdog Policy

- The `RecoilJets pipeline watchdog` heartbeat should be paused by default.
  Treat it as campaign-scoped automation, not a permanent always-on monitor.
- Turn the watchdog on only when there is concrete evidence of an active
  RecoilJets production campaign in the current conversation: the user says
  jobs were submitted, terminal output shows Condor/DAG cluster IDs, a
  RecoilJets pipeline email reports a running/check/ready/failure state, or the
  user explicitly asks Codex to watch a production.
- Before activating, record the campaign in `agent_context/STATUS_DASHBOARD.md`
  or `agent_context/TASK_BOARD.md` with dataset, command/cluster ID if known,
  evidence source, and what condition should turn the watchdog off.
- Turn the watchdog off when the campaign is done: outputs were pulled and
  validated, failure/held state was diagnosed and no longer needs monitoring,
  the user says to stop watching, or there has been no actionable campaign
  state left to track.
- Watchdog output must remain action-oriented: failures, held/stale jobs, ready
  outputs, or exact next commands. It must not submit jobs, transfer files,
  mutate SDCC, type passwords, or mark a pipeline stage resolved.

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
- When a RecoilJets/SDCC/Condor job or DAG is actively running and the user asks
  for progress, status, "what do you see", whether it is okay, whether it is
  stuck, or similar, always check relevant Gmail pipeline emails in addition to
  the visible terminal/queue state before answering. First report any fresh
  `READY`, `CHECK`, `FAILED`, held, rescue, or profiling email signal; then
  summarize the visible terminal state; then ask for exactly one compact SDCC
  diagnostic command only if more live evidence is needed.
- For active overnight or long-running RecoilJets campaigns recorded in
  `agent_context/STATUS_DASHBOARD.md`, the next progress/resume turn should be
  email-first even if the user only says "it finished", "what now", or "check
  this". Inspect fresh `RecoilJets Pipeline` Gmail messages, mark consumed
  pipeline emails read, update the dashboard with READY/CHECK/FAILED evidence,
  and only then request one compact SDCC diagnostic if email does not settle the
  state.
- In live-run progress checks, think like the user trying to follow the
  production in real time: look for current cluster/DAG ids, idle/running/held
  counts, most recent worker stdout/err/log hints, DAGMan final/merge stage
  state, profile rows, output roots appearing, and whether the next expected
  stage has started. Prefer one high-signal command that prints a compact
  summary over several piecemeal commands.

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
- Plot colors must be presentation-safe on projectors and large screens: avoid
  pale yellow or low-contrast hues for important lines, bars, annotations, or
  heatmap cells; use colorblind-aware, high-contrast palettes; and make heatmap
  annotation text switch between black and white based on the cell brightness.
  If a plot uses a very light color, verify the label/number remains readable
  after viewing the PNG at slide size.
- Plot labels must name the actual comparison in plain physics terms. Do not
  use shorthand such as `PPG12-like` in figure legends or titles when the plot
  can instead say what changed, for example `No centrality input`,
  `Centrality as input`, or `Separate centrality-bin BDTs`.
- Match the local ROOT style in `macros/sPhenixStyle.C`/`.h` for final or
  slide-candidate plots: call `SetsPhenixStyle()` or reproduce its essentials
  (font 42, no ROOT title/stat boxes, white background, borderless legends,
  readable margins, tick marks on both axes). The canonical experiment label for
  new plots is `#it{#bf{sPHENIX}} Internal`: the word `sPHENIX` is bold italic,
  `Internal` is regular upright. Do not use `#bf{sPHENIX Internal}` or plain
  `sPHENIX Internal`. If editing an existing legacy macro whose nearby plots
  consistently use `#bf{sPHENIX} #it{Internal}`, preserve that local convention
  rather than mixing styles within one plot family.
- For Python/matplotlib diagnostic plots intended for slides, emulate the same
  label visually: bold italic `sPHENIX` followed by upright `Internal`, with the
  collision-system/sample line below it. Do not let matplotlib degrade the label
  to unstyled plain `sPHENIX Internal`; if mathtext cannot reproduce it cleanly,
  split the label into separately styled text objects.
- When the user asks Codex to generate or regenerate plots for slides, Codex
  must first show the generated PNGs directly in the conversation and discuss
  any physics or formatting concerns before editing Google Slides. Iterate on
  the plot images in chat until the user explicitly says the plots are ready for
  slides. Only after that approval should Codex insert or replace the plots in
  the requested Google Slides deck.
- When the user asks to "show plots in chat" or asks for plot candidates, Codex
  should respond with the strongest slide-candidate PNGs directly embedded in
  the conversation, not only filenames. For each candidate, give a compact
  plain-language description of what is plotted, what is learned from it, why it
  might or might not belong on a slide, and any physics/formatting caveat. Work
  through candidates one by one when useful so the user can decide whether to
  add, revise, or skip each plot.
- After generating or selecting plot images, include clickable markdown links to
  each PNG and to the containing output folder when practical, so the user can
  open the exact files directly. Use absolute local paths in the links and keep
  the labels human-readable, for example `[ROC comparison](</abs/path/plot.png>)`
  and `[slideReady folder](</abs/path/slideReady>)`.
- For slide use, generated PNGs must contain only the scientific plot itself:
  axes, data, plot legend, and in-plot analysis labels that belong to the
  figure. Never bake slide titles, explanatory bullets, takeaway text, gray
  callout boxes, footers, or other slide layout material into the PNG. Those
  elements must be native editable Google Slides text/shapes so the user can
  revise them directly.
- Plot titles inside PNGs must describe the plotted observable or comparison
  directly, not the slide's narrative takeaway. Use publication-style figure
  titles such as `XGBoost split-gain fraction from shower-width inputs`, while
  keeping interpretive messages like "why this matters" in native slide text.
- If a generated plot looks physically surprising or inconsistent with the
  user's expectation, pause the slide-update step and debug the plotted inputs,
  histogram normalization/scaling, numerator/denominator definitions, and
  relevant ROOT contents before presenting it as slide-ready.

# Slide Style Preferences

- For Google Slides and presentation deck edits, always use Times New Roman
  for all native slide text: titles, body text, callouts, captions, footers,
  page numbers, and labels created as editable slide text. Preserve plot-internal
  fonts unless regenerating the plot itself.
- For Google Slides and other presentation edits, prefer body text at least
  14 pt when possible. Use line spacing around 1.2-1.8 so text fills the slide
  space cleanly instead of looking cramped or stranded.
- The user's preferred readable update-slide style is closer to their corrected
  slide 13: larger 14+ pt body text, well-filled white space rather than sparse
  small text, at least about 1.3 line spacing between body lines, and header or
  footer callouts placed in soft gray ovular/rounded boxes when that helps the
  slide read as a polished meeting slide.
- For gray rounded slide callouts/text boxes, avoid visible black borders.
  Match the outline color to the gray fill unless the user explicitly asks for
  a contrasting border.
- For colored slide callouts, use a two-object structure by default: create the
  rounded colored shape as a background object with no text, then place a
  separate transparent Times New Roman text box above it. This is especially
  important for bottom takeaway bands and RHS explanation boxes because it keeps
  text alignment, spacing, and later edits clean.
- For RHS gray explanation boxes, separate conceptual chunks with deliberate
  paragraph spacing after each short paragraph or heading block. The box should
  read as grouped thoughts with comfortable air between them, not as one dense
  paragraph stack.
- For bottom yellow takeaway bands, center the takeaway text in its own text box
  over the yellow rounded background. The text should look calm, horizontally
  centered, vertically balanced, and easy to read at presentation distance.
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
- Codex may also perform read-only SDCC diagnostics from the local Mac shell
  using the user's SSH agent when live campaign status is needed. The preferred
  order is:
  1. inspect the already-visible SSH terminal if it is on the relevant SDCC
     node;
  2. if the user has asked Codex to check directly and the task will involve
     multiple SDCC checks, open one persistent SSH session and reuse it for the
     burst of work;
  3. for a genuinely isolated quick check, run one compact read-only SSH
     diagnostic from the local Mac shell;
  4. if direct SSH/auth/host-key/proxy setup fails, stop retrying and fall back
     to the visible terminal or one paste-ready command for the user.
- Persistent SSH session rule:
  During an active SDCC workflow, Codex should avoid a pattern of many
  short-lived logins that each run one tiny command. Prefer a single SSH shell
  session for repeated `condor_q`, `tmux`, `tail`, `find`, `python3`, or report
  inspection commands. Keep the session scoped to the current task, exit it
  cleanly when the workflow is done, and do not leave unnecessary remote shells
  running. One-off SSH remains fine for small status probes, but it should not
  be the default for a sustained diagnostic or production-followup burst.
- For direct read-only SSH diagnostics, keep the command boring and robust.
  Start with the user's launchd SSH agent:
  `SSH_AUTH_SOCK="$(launchctl getenv SSH_AUTH_SOCK)"`. Prefer the established
  nested-login pattern over fragile `ProxyJump` if stdio forwarding is blocked:
  `ssh -o BatchMode=yes patsfan753@ssh.sdcc.bnl.gov "ssh sphnxuserNN.sdcc.bnl.gov '...read-only command...'"`.
  Before a complex check, run a tiny probe such as `hostname; pwd` on the target
  node. If that probe fails, do not keep improvising quote-heavy SSH commands.
- Keep direct SSH checks strictly read-only: queue/status inspection, tmux/log
  tails, `find`/`ls`/`grep`/`condor_q`/`condor_history`, and compact report
  summaries. Do not submit Condor jobs, run merge/production scripts,
  remove/hold/release jobs, transfer files, edit remote files, change
  permissions, or otherwise mutate SDCC through SSH unless the user explicitly
  asks for that specific action in the current task.
- For multi-line or quote-heavy SDCC diagnostics, avoid nested one-liner
  gymnastics. Either ask to drive the visible SSH terminal for that diagnostic,
  or give the user one paste-ready `bash` block to run on the target node. The
  command should print compact summaries and write any large raw output to an
  obvious temporary path.
- When checking a campaign by SSH, always keep lanes separated: identify the
  submit host, campaign tag, YAML/config, dataset, output root, and cluster ids
  before drawing conclusions. Do not treat errors visible in one terminal/node
  as belonging to another active campaign without matching the campaign tag,
  output root, or cluster id.
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
- When a suggested SDCC command is a submission command, tell the user that the
  run should be tracked after submission. After the user runs it or terminal
  output confirms it started, inspect the prompt/terminal for the submit node
  and record the active campaign details in the local tracking files before
  moving on.
- When a suggested SDCC submission consists of more than one related command,
  wrap it into a single paste-ready `cat > /tmp/... <<'EOF'` driver script plus
  `bash /tmp/...` invocation, or an equivalently robust one-block command.
  Avoid asking the user to manually run/edit multiple long environment-variable
  commands in sequence. The driver should fail fast, log everything, and print
  enough queue/report paths for Codex to watch the campaign without guessing.
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
- `/Users/patsfan753/Desktop/ThesisAnalysis/scripts/recoiljets_cleanup.sh`
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

`/Users/patsfan753/Desktop/ThesisAnalysis/macros/AnalyzeRecoilJets.h` is an
offline analysis header only. Do not make SDCC/Condor production, merge, or
final-stitch scripts include it or depend on it being present in the remote
analysis checkout. If remote stitching needs helper logic from that header,
copy the small required logic into the SDCC-side script or macro explicitly.

When a pushed SDCC-side change touches compiled source or headers, tell the
user to rebuild from the directory that owns the changed code, not vaguely from
an arbitrary shell location. Examples:

```bash
cd /sphenix/u/patsfan753/scratch/thesisAnalysis/src_AuAu
make clean
makeProject

cd /sphenix/u/patsfan753/scratch/thesisAnalysis/src
make clean
makeProject
```

For shell-script, YAML, or macro-only changes that do not require relinking,
say explicitly that no `makeProject` rebuild is needed.

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
- When a Condor pass, smoke test, scaled-trigger study, SIM production, merge
  stage, or local pull is finished or superseded, consider whether SDCC bulk
  artifacts are now disposable. If cleanup is applicable, proactively remind
  the user to run the guarded helper, starting with a dry-run:

```bash
./scripts/recoiljets_cleanup.sh dryrun dataset <dataset>
./scripts/recoiljets_cleanup.sh apply dataset <dataset>
./scripts/recoiljets_cleanup.sh dryrun smoke
./scripts/recoiljets_cleanup.sh apply smoke
```

  Use `dryrun all` / `apply all` only when the user is intentionally starting a
  broad fresh pass. Never imply cleanup is required while live Condor jobs may
  still reference the outputs; the helper should refuse bulk deletion unless
  the user's Condor queue is empty.
- AuAu tight-BDT sidecar artifacts should stay organized and explicitly
  lifecycle-managed:
  - extraction smoke/full roots under
    `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_<timestamp>/`
    are disposable once their row counts, model training, and any needed
    diagnostics have been consumed;
  - matching SDCC DAG/submit roots under
    `condor_sub/auauTightBDT_<timestamp>/` are disposable after the DAG has
    finished and its logs are no longer needed for debugging;
  - local SDCC checkout test products under `local_bdt_training_outputs/` are
    disposable after a better smoke/full extraction exists;
  - trained model directories under `bdt_models/tight_<timestamp>/` are not
    disposable unless a newer model is validated and the old one is explicitly
    superseded.
- Keep the working checkout tidy during diagnostics: avoid leaving dry-run
  `condor_sub/auauTightBDT_*` folders, temporary manifests, or fake-list test
  files in the repo after local checks. Prefer timestamped run roots for real
  SDCC artifacts, and remove only generated test clutter you created.

# ROOT Merge Strategy

- For large ROOT merge/stitch workloads, prefer layered Condor merging over one
  huge `hadd` or one huge recursive ROOT merge. The user's experience is that
  staged merging is usually fastest and most Condor-friendly: for example,
  with about `10k` files, first merge groups of roughly `200`, then merge those
  intermediate outputs in groups of roughly `10`, then do one final small
  merge. Apply the same principle to RecoilJets SIM stitching and large data
  merges: split memory-heavy work into many moderate Condor jobs, then assemble
  deterministic partials into the exact canonical final ROOT output. Avoid
  raising memory to very large values when a clean layered merge can keep jobs
  smaller and start faster.

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
when the current task clearly calls for SDCC transfer, status comparison, or
output retrieval. Before every helper run, Codex must explicitly say what it is
about to do, show the exact command, classify the interaction as read-only
status/diff, upload, download, local-only merge, or commit/push after upload,
and name the expected local/SDCC paths or files when practical:

```bash
./scripts/sftp_push_recoiljets.sh ...
./scripts/sftp_get_recoiljets_outputs.sh ...
```

If the helper opens an interactive password prompt, the user enters the SDCC
password manually. Codex must not ask for, read, store, repeat, or type the
password. If the user objects, pauses, or asks to run the helper manually,
Codex must stop and provide the same command for the user to run manually.

Routine low-risk helper runs do not require an additional yes/no confirmation
after Codex has given that explicit preflight notice. Routine low-risk runs
include read-only `status`/`diff`, pulling requested ready outputs, local-only
`mergeLocalSim`, uploading the exact mapped files Codex just edited for a
trivial or requested pipeline fix, and using the uploader's selected-file
`--commit-push` form after a successful upload.

Codex must still ask for explicit user confirmation before broad or risky
transfer actions: uploading broad groups such as `pipeline`, `all`, or
`changed`; downloading large datasets when the requested product is ambiguous
or may overwrite valuable local outputs; using `--commit-push` for anything
substantive or surprising; extending the uploader allowlist; or any action
that could overwrite, delete, submit, hold, release, or otherwise mutate
remote SDCC state beyond copying the selected mapped files.

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
  SDCC password. Codex may run them locally after the explicit preflight notice
  above; if the user objects or wants to run the check manually, give the
  command instead.

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
- When Codex edits any mapped file and has not already run the upload helper,
  its final response must include the exact `./scripts/sftp_push_recoiljets.sh
  ...` command the user should run. This applies even for tiny warning fixes,
  one-line changes, or follow-up edits to files changed earlier in the
  conversation.
- For substantive SDCC pipeline edits, assume the user usually wants the local
  GitHub repo updated after the successful SDCC transfer. Include the uploader's
  built-in commit/push form by default, with the exact selected files/groups and
  a descriptive message:

```bash
./scripts/sftp_push_recoiljets.sh <changed-files-or-groups> --commit-push -m "Concise change summary"
```

- For clean requested edits to mapped SDCC pipeline files, Codex should offer to
  run the exact selected-file `scripts/sftp_push_recoiljets.sh ... --commit-push`
  command itself after tests pass, not merely leave the command for the user.
  Ask first and explain what changed, which local paths will be uploaded, which
  SDCC paths will be overwritten, and whether `makeProject` or another rebuild
  is required afterward.
- Be especially explicit before uploads that add or remove mapped files. Name
  each new/deleted file, why it is part of the requested change, and wait for
  user approval before running the uploader.
- If the user explicitly asks for upload-only, omit `--commit-push`.
- The commit message should describe the actual code/config change, not just say
  "changes".
- The uploader commits only the selected transferable files after successful
  SFTP. It must never use `git add .` for this workflow.
- Do not SSH directly.
- Do not run `sftp` directly.
- Do not ask for or handle the SDCC password.
- The local uploader/getter helpers may be run by Codex after the explicit
  preflight notice above. The user still enters any SDCC password prompt
  manually.
- When visible SDCC terminal output or parseable stage emails show that merged
  outputs are ready, Codex may run `scripts/sftp_get_recoiljets_outputs.sh`
  locally after the explicit preflight notice above to pull the relevant
  dataset into `InputFiles/...`, unless the product or overwrite risk is
  ambiguous enough to require confirmation.
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
scripts/recoiljets_cleanup.sh             -> scripts/recoiljets_cleanup.sh
scripts/root_in_analysis_env.sh           -> scripts/root_in_analysis_env.sh
scripts/submit_auau_bdt_widthstudy_pt1530_wp080.sh -> scripts/submit_auau_bdt_widthstudy_pt1530_wp080.sh
scripts/RecoilJets_Condor_AuAu.sh         -> RecoilJets_Condor_AuAu.sh
scripts/RecoilJets_Condor_submit.sh       -> RecoilJets_Condor_submit.sh
scripts/RecoilJets_Condor.sh              -> RecoilJets_Condor.sh
scripts/train_auau_jet_residual_bdt.py    -> scripts/train_auau_jet_residual_bdt.py
scripts/train_auau_photon_bdt.py          -> scripts/train_auau_photon_bdt.py
scripts/validate_auau_tight_bdt_on_sim.py -> scripts/validate_auau_tight_bdt_on_sim.py
macros/analysis_config.yaml               -> macros/analysis_config.yaml
macros/analysis_config_auau_bdt_widthstudy_pt1530_wp080.yaml -> macros/analysis_config_auau_bdt_widthstudy_pt1530_wp080.yaml
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
