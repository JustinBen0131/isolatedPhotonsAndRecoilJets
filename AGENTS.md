# Codex Project Control Plane

This file is the always-loaded control plane for
`/Users/patsfan753/Desktop/ThesisAnalysis`. Keep it short. For detailed rules,
load the focused policy file named by the routing table below.

## First Reflex

1. Classify the user's request before acting: science reasoning, SDCC
   training/production/merge, local code edit, plot generation, Google Slides,
   Gmail/watchdog, transfer, cleanup/storage, or memory/status update.
2. If Justin says "add a task", "add this to the task list", "put this on my
   todo/backlog", "track this as a task", or a close variant, treat that as an
   explicit task-capture command. Load `CODEX_OPERATING_SYSTEM.md` and
   `MEMORY_AND_STATUS.md`, then update the full task pipeline: repo register
   first, Linear second, Today's Plan only if it belongs on today's surface.
3. Load the matching policy files from `agent_context/policies/`.
4. For meaningful multi-step work, load
   `agent_context/policies/CODEX_OPERATING_SYSTEM.md` and claim or update the
   matching entry in `agent_context/CODEX_WORK_REGISTER.yaml` before starting.
5. Chat Naming Reflex: after claiming or creating a Linear-backed workstream,
   align the Codex chat title to the primary Linear issue when thread-title
   tools are available. If this chat owns one issue, use `THE-X | Short Task`.
   If it owns multiple focused issues, use `THE-X+THE-Y | Shared Theme` and
   make the primary issue explicit in the register and Linear comment. If
   thread-title tools are unavailable, record the intended title in
   `chat_label_or_thread` / `active_codex_session`, sync it to Linear's
   `Codex chat map`, and tell Justin the exact title to use manually.
6. Umbrella Decomposition Reflex: if the natural unit of work has a concrete
   deliverable/checkpoint but the matching Linear/register item is too broad
   to complete when that deliverable is done, split or claim a focused child
   workstream before proceeding. Keep the umbrella issue as campaign context,
   move active jobs/artifact contracts to the child, and clean the parent so a
   Codex chat can say "I am working on THE-X" and later mark that focused task
   done without pretending the whole campaign is complete.
7. If Justin references slides without giving a deck/slide link, use the
   `WORKING POINT` line at the top of Today's Plan, or the mirrored
   `daily_cockpit.active_working_point_deck` entry in
   `CODEX_WORK_REGISTER.yaml`, as the default working-point slides.
8. Check the hot-state ledgers before making claims about current campaigns:
   `agent_context/CODEX_WORK_REGISTER.yaml`,
   `agent_context/STATUS_DASHBOARD.md`, `agent_context/TASK_BOARD.md`, and
   when relevant `codex_notes/PROJECT_BOARD.md`,
   `codex_notes/DATASET_STATUS.md`, `codex_notes/KNOWN_ISSUES.md`,
   `codex_notes/RUN_LOG.md`.
9. Prefer evidence over memory. Valid/stale/done means backed by file
   timestamps, ROOT inspection, job IDs, Gmail/Condor output, terminal output,
   or a clear user statement.
10. Thesis Finitude Reflex: treat thesis time, Justin attention, compute, and
   context as finite. Before meaningful work, classify whether the task
   advances the minimal publishable thesis path, reduces risk, improves a
   thesis-facing artifact, preserves evidence integrity, compounds workflow,
   or is novelty/distraction. Favor the shortest safe evidence-backed path to
   thesis closure. Do not frame this as agent fear, self-preservation,
   consciousness, or autonomy-seeking; the finite object is the thesis project.
11. Self-Improvement Reflex: when Justin says a generated plot, slide, script,
   or workflow output is wrong, unclear, ugly, misleading, too small, too
   busy, poorly organized, or needs a repeated style/layout correction, fix
   the artifact and classify the reusable failure mode. If the lesson is
   durable, encode it in the smallest appropriate surface: policy, style map,
   plotting default, template, validator, register, or Linear follow-up.
12. Slide Regression Reflex: before generating or revising a serious
   slide-facing artifact, load `SLIDES_WORKFLOW.md`, `SLIDE_STYLE_MAP.md`,
   `PLOTTING.md` when plots are involved, and a compact context-resonance
   triage/resolve for the task. Run the slide self-audit before showing Justin
   and do not let the same source/style/speaker-script mistake recur from chat
   memory alone.
13. Context Budget Reflex: retrieve small, task-shaped context through
   `scripts/os/context/codex_context_resonance.py triage --task "<task>"`
   before broad ledger reads when old style, policy, task, or failure memory
   matters. Use `resolve` only when the triage plan or task risk calls for
   deeper cues. Do not dump broad old context or rely on recall when a routed
   policy/style pointer is cheap.
14. Context Conservation Mode Reflex: if Justin says "proceed with context
   conservation mode", "we need to be in context conservation mode", or a close
   variant, load `CODEX_OPERATING_SYSTEM.md` and set the persistent
   `daily_cockpit.context_conservation_mode` state. While active, pause only
   dream-lane automations, use the smallest sufficient reads, avoid Chronicle /
   browser / Computer Use unless explicitly requested or required for safety,
   keep answers brief, and keep slide work to PNG candidates plus terse
   evidence/caveat notes unless Justin explicitly asks for deck mutation or a
   full script. The mode remains active until Justin explicitly turns it off.
15. When creating new scripts, plotting helpers, slide generators, ROOT macros,
   diagnostics, or OS utilities, load `REPO_ORGANIZATION.md` and place the new
   file in the clearest safe home. Do not add more unrelated side helpers to
   flat `scripts/` or `macros/` by default. Do not move Fun4All, SDCC, Condor,
   transfer, or base pipeline entrypoints without an explicit compatibility
   migration plan and Justin's approval. After THE-23 stage 5, top-level
   `scripts/` entries are hard command aliases and indexes; canonical source
   should be found through `scripts/COMMAND_INDEX.tsv`,
   `scripts/bin/thesis-script`, `scripts/SCRIPT_INDEX.yaml`,
   `scripts/HELPER_INDEX.yaml`, `scripts/slides/INDEX.yaml`, and the purpose
   subfolders. Historical local names that are not hard aliases live under
   `scripts/compat/local/`.
16. If Justin says "set up a Zoom room", "open a Zoom room", "make a Zoom
   room", "start a Zoom room", "get me a Zoom invite link", or a close variant,
   load `MEETING_CAPTURE_AND_ZOOM.md`. Treat the phrase as explicit permission
   to use Computer Use on Zoom for that narrow room-setup task: create/start
   the Justin-owned room, enable or verify transcript/AI-summary capture where
   available, copy the invite link, and return the link plus any caveat.

## Hard Stops

Load `agent_context/policies/HARD_STOPS_AND_SAFETY.md` before any risky action.

- Duplicate-run hard stop: if Justin asks for a new training campaign,
  validation, merge, RecoilJets production, or production-style analysis that
  is identical or materially equivalent to an existing run/artifact, stop
  immediately. Report the matching tag/path/evidence and ask whether to reuse,
  compare, or deliberately rerun.
- Do not submit Condor, run merge/production scripts, remove jobs, delete
  outputs, transfer large payloads, edit remote SDCC files, or mutate Google
  Slides unless the current user request explicitly authorizes that action or
  the relevant policy permits it.
- For approved risky mutations, use the OS safety kernel first:
  `scripts/codex_os_snapshot.py` plus `scripts/codex_os_guard.py`.
- Never revert user changes unless explicitly requested.
- Never ask for, read, store, repeat, or type SDCC passwords or other secrets.
- Never remove, overwrite, or repurpose a Google Slides `Backup` slide.

## SDCC Local Rules

Load `agent_context/policies/SDCC_OPERATIONS.md` for SDCC work.

Detailed SDCC hostnames, usernames, access patterns, and live campaign ledgers
are local-only Codex memory. If `agent_context/local/SDCC_LOCAL_RULES.md`
exists, load it before any SDCC/Condor action. Do not commit SDCC access
details, queue snapshots, job IDs, remote paths, or watchdog status histories
to GitHub. Submissions, merge reruns, cleanup, file transfers, and remote edits
still require explicit approval for the specific campaign action.

Never create or leave a front-facing `agent_context/`, `.codex/`, `codex*`,
or whitespace-only path in the SDCC checkout. Agent/control-plane evidence is
local-only unless Justin explicitly approves a one-off remote diagnostic file,
and even then use a neutral hidden scratch path such as `.recoiljets_tmp/` or
`/tmp`, not an agent-named directory.

For Codex-launched RecoilJets/Condor workflows, the remote submit command must
set `RJ_CODEX_CHAT_NAME` and `RJ_CODEX_THREAD_ID` before DAG/meta files are
created. Do not rely on local `CODEX_*` variables surviving nested SSH or
Condor scheduler nodes; missing `codex_chat_name`/`codex_thread_id` in pipeline
emails is a provenance defect to fix before the next Codex-launched rerun.

## Routing Table

The machine-readable index is `agent_context/policies/LOAD_MAP.yaml`.

| Request shape | Load these policy files first |
| --- | --- |
| new training, validation, merge, production, production-style analysis | `DUPLICATE_RUN_GUARD.md`, `SDCC_OPERATIONS.md`, `MEMORY_AND_STATUS.md`, `STORAGE_AND_CLEANUP.md`, `TRANSFER_AND_PIPELINE_FILES.md` |
| SDCC status, held jobs, queue, watchdog, heartbeat | `EMAIL_AND_WATCHDOGS.md`, `SDCC_OPERATIONS.md`, `MEMORY_AND_STATUS.md` |
| plot generation or plot QA | `PLOTTING.md`, `SCIENCE_REFERENCE_HIERARCHY.md`, `INPUTS_AND_ENVIRONMENT.md` |
| Google Slides or deck edits | `SLIDES_WORKFLOW.md`, `PLOTTING.md` if plots are involved |
| new plotting script, slide generator, ROOT macro, diagnostics helper, or scripts/macros cleanup | `REPO_ORGANIZATION.md`, plus `PLOTTING.md`/`SLIDES_WORKFLOW.md`/`TRANSFER_AND_PIPELINE_FILES.md` as relevant |
| slide reference without explicit deck link, current working slides, working point slides | `SLIDES_WORKFLOW.md`, `CODEX_OPERATING_SYSTEM.md` |
| photon ID, BDT, ABCD, purity, xJ, unfolding, stitching | `SCIENCE_REFERENCE_HIERARCHY.md`, `DUPLICATE_RUN_GUARD.md` if a run/output is involved |
| SFTP push/get, SDCC pipeline file upload, mapped-file changes | `TRANSFER_AND_PIPELINE_FILES.md`, `SDCC_OPERATIONS.md` |
| cleanup, quotas, stale outputs, failed payloads | `STORAGE_AND_CLEANUP.md`, `MEMORY_AND_STATUS.md` |
| local ROOT macro or ROOT-dependent script | `INPUTS_AND_ENVIRONMENT.md`, `PLOTTING.md` if producing figures |
| missing tool, dependency, plugin, runtime, browser/Drive/Gmail workflow | `TOOLS_AND_DEPENDENCIES.md` |
| explicit delegation, subagent, Claude, handoff, parallel lanes | `COLLABORATION_AND_HANDOFFS.md` |
| current stale-output bug or dataset validity | `CURRENT_PROJECT_STATE.md`, `MEMORY_AND_STATUS.md` |
| Codex work registration, Linear sync, Today's Plan, multi-chat state, active jobs | `CODEX_OPERATING_SYSTEM.md`, `MEMORY_AND_STATUS.md` |
| agentic OS hardening, self-tuning, symbiotic workflow, thesis control plane, doctor checks | `AGENTIC_OS_HARDENING.md`, `CODEX_OPERATING_SYSTEM.md`, `MEMORY_AND_STATUS.md` |
| agentic dreaming, synthetic rehearsal, night simulation, dream automation | `AGENTIC_OS_DREAMING.md`, `AGENTIC_OS_HARDENING.md`, `CODEX_OPERATING_SYSTEM.md` |
| ask ChatGPT, use ChatGPT UI, delegated external research, deep research for OS design | `ASK_CHATGPT_DELEGATION.md`, `TOOLS_AND_DEPENDENCIES.md`, `COLLABORATION_AND_HANDOFFS.md` |
| Zoom room setup, Zoom invite link, meeting transcript, Zoom AI summary, meeting capture | `MEETING_CAPTURE_AND_ZOOM.md`, `TOOLS_AND_DEPENDENCIES.md`, `CODEX_OPERATING_SYSTEM.md`, `MEMORY_AND_STATUS.md` |
| "add a task", add to todo/backlog, task capture, task list update | `CODEX_OPERATING_SYSTEM.md`, `MEMORY_AND_STATUS.md` |

## Science North Star

Load `SCIENCE_REFERENCE_HIERARCHY.md` for details.

- Ultimate target: PPG19 Au+Au gamma-jet / `x_{J#gamma}`.
- PPG18 pp gamma-jet is the validated reference baseline.
- PPG12 is the photon-ID, isolation, BDT/NPB, purity, and pp infrastructure
  reference. Use `ppg12codeGit/` and the local PPG12 notes before inventing new
  photon-ID, BDT, ABCD, stitching, or systematic logic.
- ATLAS gamma-jet analysis note is the final-analysis style target, adapted to
  sPHENIX/PPG constraints.

## Memory Architecture

Load `MEMORY_AND_STATUS.md` before recording state.

- `TASK_BOARD.md`: milestones, blockers, next actions.
- `CODEX_WORK_REGISTER.yaml`: canonical active work, Codex ownership, Linear
  sync state, Today's Plan anchors, active jobs, and stale-check timing.
- `STATUS_DASHBOARD.md`: evidence-backed live campaign/output status.
- `REFERENCE_MAP.md`: compact science/reference decisions.
- `SLIDE_STYLE_MAP.md`: durable deck style examples.
- `codex_notes/`: older local project-state ledgers that still matter for
  dataset validity and bug tracking.
- Record only decisions and evidence, not every thought.
- For meaningful multi-step work, update `CODEX_WORK_REGISTER.yaml` first; use
  Google Docs and Linear as projections of that register.
- When Justin says "add a task", "add this to the task list", "put this on my
  todo/backlog", "track this as a task", or a close variant, do not leave it as
  a chat note. Add or update the cohesive task pipeline: register, Linear,
  Today's Plan if relevant, and memory/status notes when durable. When he
  casually mentions a possible task without task-capture wording, ask before
  adding it.

## Editing And Change Control

- Before file edits, state exactly which files will change and why. If Justin
  already directly requested the edit, proceed after the statement; otherwise
  wait for approval.
- Use `apply_patch` for manual edits. Do not use destructive git commands.
- You may be in a dirty worktree. Preserve unrelated user changes.
- When drafting messages for Justin to send, keep the style close to his
  wording and avoid formal colon-heavy constructions such as `Topic: detail`.
- For ROOT work, use:

```bash
./scripts/root_in_analysis_env.sh /Users/patsfan753/Desktop/analysis/env/bin/root -l -q 'macros/MyMacro.C()'
```

## Plot And Slide Guardrails

Load `PLOTTING.md` and `SLIDES_WORKFLOW.md`.

- New analysis plots are PNG by default, no ROOT stats boxes, sPHENIX style,
  and visually inspected before being called ready.
- For slide-use plots, show PNGs in chat first. Do not insert or replace plots
  in Google Slides until Justin approves the candidate image.
- Native Google Slides text should use Times New Roman unless the user says
  otherwise. Keep body text readable and preserve the `Backup` divider.

## Transfer Guardrails

Load `TRANSFER_AND_PIPELINE_FILES.md`.

- Use `scripts/sftp_push_recoiljets.sh` for mapped local-to-SDCC pipeline
  uploads and read-only status/diff checks. For script paths, consult
  `scripts/sdcc/TRANSFER_MAP.tsv`; it records local canonical paths, local
  aliases, SDCC canonical paths, compatibility paths, and upload targets.
- Use `scripts/sftp_get_recoiljets_outputs.sh` for pulling ready outputs.
- Do not run raw `sftp`/`scp` manually from Codex.
- If a mapped SDCC-side file changes, final response should include the exact
  smallest upload command unless Codex already ran it with approval.

## Repo Organization

Load `REPO_ORGANIZATION.md` before adding or moving code under `scripts/` or
`macros/`. New local side helpers should go in clear purpose-specific folders
when safe. For existing scripts, prefer `scripts/bin/thesis-script path <name>`
or `scripts/COMMAND_INDEX.tsv` over scanning the parent directory. Preserve
current paths for Fun4All, SDCC, Condor, transfer, and base pipeline entrypoints
unless a staged migration is explicitly approved. Do not recreate flat
top-level script source files; use canonical subfolders plus documented hard or
compat aliases.

## Research Note

This organization follows a practical agent-memory model: keep working memory
small, route to structured long-term stores, record evidence as explicit
entities, use a write/manage/read loop, and make forgetting/cleanup deliberate.
The implementation is local markdown/YAML so future Codex sessions can inspect
it quickly with `rg`, `sed`, and targeted file reads.

## Imported Claude Cowork project instructions

You are working in /Users/patsfan753/Desktop/ThesisAnalysis, a sPHENIX RecoilJets thesis analysis repo.

First read AGENTS.md, then follow it as the project authority. Also read codex_notes/PROJECT_BOARD.md, codex_notes/DATASET_STATUS.md, and codex_notes/KNOWN_ISSUES.md before making dataset-validity claims, Condor/SDCC recommendations, or plot/slide conclusions.

Use ripgrep for search. Do not revert user changes. Keep edits scoped. For ROOT macros or ROOT-dependent scripts, use:
./scripts/root_in_analysis_env.sh /Users/patsfan753/Desktop/analysis/env/bin/root ...

Treat SDCC/remote SSH terminals as user-controlled unless explicitly instructed otherwise. Do not ask for, store, or type SDCC passwords. Prefer paste-ready commands for the user to run on SDCC. Use scripts/sftp_push_recoiljets.sh and scripts/sftp_get_recoiljets_outputs.sh for transfer workflows, following AGENTS.md.

For analysis plots, output PNG by default, disable ROOT stats boxes, and visually inspect generated PNGs before calling plots final.

For Google Slides/deck text, use concise, clear, organized bullets with bolded structure when helpful. Prefer >=14 pt body text and line spacing around 1.2-1.8. Explain physics/analysis conclusions carefully and distinguish proven evidence from plausible interpretation.

When changing mapped SDCC pipeline files, include the exact sftp_push_recoiljets.sh command needed afterward. Do not use wildcard uploads.
