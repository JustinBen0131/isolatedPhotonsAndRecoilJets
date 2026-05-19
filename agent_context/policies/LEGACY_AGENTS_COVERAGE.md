# Legacy AGENTS Coverage

This file maps the original monolithic `AGENTS.md` headings to the organized
policy shards. It exists to prevent memory loss during future refactors.

## Coverage Map

| Original heading | New home |
| --- | --- |
| `Codex SDCC SSH Rule` | `AGENTS.md`, `SDCC_OPERATIONS.md` |
| `Scientific Mission` | `AGENTS.md`, `SCIENCE_REFERENCE_HIERARCHY.md` |
| `Reference Hierarchy` | `SCIENCE_REFERENCE_HIERARCHY.md` |
| `Analysis Decision Rules` | `SCIENCE_REFERENCE_HIERARCHY.md`, `CURRENT_PROJECT_STATE.md` |
| `Plot And Architecture Reuse` | `PLOTTING.md`, `SCIENCE_REFERENCE_HIERARCHY.md` |
| `Slide Integration Rules` | `SLIDES_WORKFLOW.md`, `PLOTTING.md` |
| `Collaboration Model` | `COLLABORATION_AND_HANDOFFS.md` |
| `Local Agent Memory` | `MEMORY_AND_STATUS.md`, `SDCC_OPERATIONS.md`, `EMAIL_AND_WATCHDOGS.md` |
| `Campaign-Scoped Watchdog Policy` | `EMAIL_AND_WATCHDOGS.md`, `STORAGE_AND_CLEANUP.md` |
| `Analysis Environment` | `INPUTS_AND_ENVIRONMENT.md`, `TOOLS_AND_DEPENDENCIES.md` |
| `Durable Analysis Context` | `SDCC_OPERATIONS.md`, `SCIENCE_REFERENCE_HIERARCHY.md`, `INPUTS_AND_ENVIRONMENT.md` |
| `RecoilJets Gmail Pipeline Emails` | `EMAIL_AND_WATCHDOGS.md` |
| `Input Files` | `INPUTS_AND_ENVIRONMENT.md` |
| `Plot Output` | `PLOTTING.md` |
| `Slide Style Preferences` | `SLIDES_WORKFLOW.md` |
| `Change Control` | `AGENTS.md`, `HARD_STOPS_AND_SAFETY.md` |
| `Local Tool Dependencies` | `TOOLS_AND_DEPENDENCIES.md` |
| `SSH Diagnostic CLI Handoff` | `SDCC_OPERATIONS.md`, `MEMORY_AND_STATUS.md` |
| `sPHENIX Pipeline Paths` | `TRANSFER_AND_PIPELINE_FILES.md` |
| `SDCC Output Storage Policy` | `STORAGE_AND_CLEANUP.md` |
| `ROOT Merge Strategy` | `SDCC_OPERATIONS.md` |
| `SDCC Transfer Workflow` | `TRANSFER_AND_PIPELINE_FILES.md`, `SDCC_OPERATIONS.md` |
| `Durable Project State Tracking` | `CURRENT_PROJECT_STATE.md`, `MEMORY_AND_STATUS.md` |

## Intentional Update

The original submit-host helper rule said Codex should tell Justin to run
`./checkCondorQ`. Justin corrected this on 2026-05-18: Codex should run
`./checkCondorQ` itself as a read-only local diagnostic when asked which SDCC
node/submit host to use. That correction lives in `SDCC_OPERATIONS.md`.

## Audit Rule

If `AGENTS.md` is refactored again, first compare the old committed headings
with this coverage map and preserve each heading's durable rule in either the
root control plane or a named policy shard.
