# Codex Project Control Plane

This file is the always-loaded control plane for
`/Users/patsfan753/Desktop/ThesisAnalysis`. Keep it short. For detailed rules,
load the focused policy file named by the routing table below.

## First Reflex

1. Classify the user's request before acting: science reasoning, SDCC
   training/production/merge, local code edit, plot generation, Google Slides,
   Gmail/watchdog, transfer, cleanup/storage, or memory/status update.
2. Load the matching policy files from `agent_context/policies/`.
3. Check the hot-state ledgers before making claims about current campaigns:
   `agent_context/STATUS_DASHBOARD.md`, `agent_context/TASK_BOARD.md`, and
   when relevant `codex_notes/PROJECT_BOARD.md`,
   `codex_notes/DATASET_STATUS.md`, `codex_notes/KNOWN_ISSUES.md`,
   `codex_notes/RUN_LOG.md`.
4. Prefer evidence over memory. Valid/stale/done means backed by file
   timestamps, ROOT inspection, job IDs, Gmail/Condor output, terminal output,
   or a clear user statement.

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
- Never revert user changes unless explicitly requested.
- Never ask for, read, store, repeat, or type SDCC passwords or other secrets.
- Never remove, overwrite, or repurpose a Google Slides `Backup` slide.

## SDCC SSH Rule

Load `agent_context/policies/SDCC_OPERATIONS.md` for SDCC work.

Do not rely on `ProxyJump` to `sphnxuserXX` from Codex. Preferred read-only
diagnostic pattern from the local Mac:

```bash
SSH_AUTH_SOCK="$(launchctl getenv SSH_AUTH_SOCK)" \
ssh patsfan753@ssh.sdcc.bnl.gov \
  "ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null sphnxuser05.sdcc.bnl.gov 'cd /sphenix/u/patsfan753/scratch/thesisAnalysis && <command>'"
```

Use nested SSH for compact read-only diagnostics. For sustained SDCC work, open
one persistent SSH session instead of many tiny logins. Submissions, merge
reruns, cleanup, file transfers, and remote edits require explicit approval for
the specific campaign action.

## Routing Table

The machine-readable index is `agent_context/policies/LOAD_MAP.yaml`.

| Request shape | Load these policy files first |
| --- | --- |
| new training, validation, merge, production, production-style analysis | `DUPLICATE_RUN_GUARD.md`, `SDCC_OPERATIONS.md`, `MEMORY_AND_STATUS.md`, `STORAGE_AND_CLEANUP.md`, `TRANSFER_AND_PIPELINE_FILES.md` |
| SDCC status, held jobs, queue, watchdog, heartbeat | `EMAIL_AND_WATCHDOGS.md`, `SDCC_OPERATIONS.md`, `MEMORY_AND_STATUS.md` |
| plot generation or plot QA | `PLOTTING.md`, `SCIENCE_REFERENCE_HIERARCHY.md`, `INPUTS_AND_ENVIRONMENT.md` |
| Google Slides or deck edits | `SLIDES_WORKFLOW.md`, `PLOTTING.md` if plots are involved |
| photon ID, BDT, ABCD, purity, xJ, unfolding, stitching | `SCIENCE_REFERENCE_HIERARCHY.md`, `DUPLICATE_RUN_GUARD.md` if a run/output is involved |
| SFTP push/get, SDCC pipeline file upload, mapped-file changes | `TRANSFER_AND_PIPELINE_FILES.md`, `SDCC_OPERATIONS.md` |
| cleanup, quotas, stale outputs, failed payloads | `STORAGE_AND_CLEANUP.md`, `MEMORY_AND_STATUS.md` |
| local ROOT macro or ROOT-dependent script | `INPUTS_AND_ENVIRONMENT.md`, `PLOTTING.md` if producing figures |
| missing tool, dependency, plugin, runtime, browser/Drive/Gmail workflow | `TOOLS_AND_DEPENDENCIES.md` |
| explicit delegation, subagent, Claude, handoff, parallel lanes | `COLLABORATION_AND_HANDOFFS.md` |
| current stale-output bug or dataset validity | `CURRENT_PROJECT_STATE.md`, `MEMORY_AND_STATUS.md` |

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
- `STATUS_DASHBOARD.md`: evidence-backed live campaign/output status.
- `REFERENCE_MAP.md`: compact science/reference decisions.
- `SLIDE_STYLE_MAP.md`: durable deck style examples.
- `codex_notes/`: older local project-state ledgers that still matter for
  dataset validity and bug tracking.
- Record only decisions and evidence, not every thought.
- When Justin says "track", "add", "remember", "record", update the relevant
  memory file. When he casually mentions a task, ask before adding it.

## Editing And Change Control

- Before file edits, state exactly which files will change and why. If Justin
  already directly requested the edit, proceed after the statement; otherwise
  wait for approval.
- Use `apply_patch` for manual edits. Do not use destructive git commands.
- You may be in a dirty worktree. Preserve unrelated user changes.
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
  uploads and read-only status/diff checks.
- Use `scripts/sftp_get_recoiljets_outputs.sh` for pulling ready outputs.
- Do not run raw `sftp`/`scp` manually from Codex.
- If a mapped SDCC-side file changes, final response should include the exact
  smallest upload command unless Codex already ran it with approval.

## Research Note

This organization follows a practical agent-memory model: keep working memory
small, route to structured long-term stores, record evidence as explicit
entities, use a write/manage/read loop, and make forgetting/cleanup deliberate.
The implementation is local markdown/YAML so future Codex sessions can inspect
it quickly with `rg`, `sed`, and targeted file reads.
