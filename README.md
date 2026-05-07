# ThesisAnalysis

This repository contains the source, macros, and scripts for Justin Bennett's
sPHENIX RecoilJets thesis analysis pipeline. It is intentionally kept small:
GitHub tracks only the code needed to understand, build, submit, and analyze the
pipeline.

## Repository Layout

- `src/`: pp-style RecoilJets analysis module and photon-cluster helpers.
- `src_AuAu/`: Au+Au / embedded RecoilJets analysis module.
- `macros/`: Fun4All steering macros, analysis configuration, ROOT plotting,
  QA, unfolding, stitching, and presentation-plot macros.
- `scripts/`: Condor submission wrappers, SDCC transfer helpers, list builders,
  merge helpers, ROOT environment wrapper, and training/diagnostic scripts.

Local analysis inputs, ROOT outputs, presentations, notes, scratch directories,
external dependencies, and copied reference repositories are intentionally not
tracked here.

## ROOT Environment

Run ROOT-dependent commands through the repository wrapper so the analysis uses
the same environment as the local `use_root_analysis` setup:

```bash
./scripts/root_in_analysis_env.sh /Users/patsfan753/Desktop/analysis/env/bin/root -l -q 'macros/MyMacro.C()'
```

Do not rely on interactive shell startup files for ROOT, RooUnfold, or analysis
library paths.

## Pipeline Overview

1. Build or update dataset/file lists with scripts such as
   `scripts/make_dstListsData.sh` and `scripts/makeThesisSimLists.sh`.
2. Configure analysis choices in `macros/analysis_config.yaml`.
3. Run the appropriate Fun4All steering macro:
   `macros/Fun4All_recoilJets.C` for pp-style workflows or
   `macros/Fun4All_recoilJets_AuAu.C` for Au+Au/embedded workflows.
4. Submit staged SDCC/Condor production with
   `scripts/RecoilJets_Condor_submit.sh`.
5. Pull merged outputs with `scripts/sftp_get_recoiljets_outputs.sh` into the
   local `InputFiles/` area.
6. Run plotting, QA, purity, unfolding, stitching, and comparison macros from
   `macros/`.

## SDCC Transfer Helpers

Use the local SFTP helper to compare or upload known pipeline files:

```bash
./scripts/sftp_push_recoiljets.sh status
./scripts/sftp_push_recoiljets.sh status changed
./scripts/sftp_push_recoiljets.sh diff changed
./scripts/sftp_push_recoiljets.sh <file-or-group>
```

Use the getter helper for merged ROOT outputs:

```bash
./scripts/sftp_get_recoiljets_outputs.sh <dataset-or-mode>
```

These helpers are the supported path for moving pipeline code and outputs
between the local checkout and SDCC.

## What Is Not Tracked

The GitHub repository intentionally excludes:

- `InputFiles/` and `dataOutput/`
- `codex_notes/`, `AGENTS.md`, and local assistant/project state
- `Presentations/`, `usefulDocs/`, and local reference PDFs
- `coresoftware_local/`, `ppg12codeGit/`, `external/`, and copied dependencies
- `.tmp*`, `.pycache/`, `.codex/`, `.claude/`, and other local scratch state

Those files can remain in the local thesis-analysis directory, but they should
not be part of the GitHub transfer surface.
