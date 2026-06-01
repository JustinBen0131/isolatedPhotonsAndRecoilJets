# Repo Organization Migration Spec

Updated: 2026-06-01
Workstream: `repo_scripts_macros_organization` / Linear THE-23
Snapshots:

- Stage 1: `20260601T184437Z-the23-scripts-macros-slide-organization`
- Stage 2: `20260601T185641Z-the23-stage2-artifact-helper-cleanup`
- Stage 3: `20260601T191159Z-the23-stage3-scripts-local-sdcc-map`
- Stage 4: `20260601T192755Z-the23-stage4-sdcc-scripts-organization`
- Stage 5: `20260601T195915Z-the23-stage5-scripts-systems-architecture`

## Purpose

Make local helper code easier for Codex to traverse without breaking protected
SDCC, Fun4All, Condor, SFTP, ROOT, or operator command contracts. The current
state keeps fixed-path command contracts available, moves implementation code
into purpose folders, and leaves compatibility wrappers or symlinks at
historical command paths when those paths are still active contracts.
After stage 5, the local and SDCC `scripts/` top-level surfaces are minimal
command layers: regular source files live in deeper domain folders, hard active
command paths stay as symlinks, and retired historical aliases live under
`scripts/compat/local/` locally or `scripts/compat/sdcc/` remotely. Use
`scripts/bin/thesis-script`, `scripts/COMMAND_INDEX.tsv`, and the generated
indexes before scanning or creating new script files.

## Inventory

The full classification is in `agent_context/REPO_ORGANIZATION_INVENTORY.tsv`.
Current classification counts after stage 5:

- `compatibility_symlink`: 89
- `data_prep`: 17
- `diagnostic`: 7
- `diagnostic_macro`: 8
- `environment_helper`: 1
- `hard_alias`: 20
- `index_or_readme`: 29
- `ml`: 16
- `os_utility`: 14
- `plotting_macro`: 48
- `plotting_utility`: 70
- `protected`: 15
- `sdcc_runtime_or_transfer`: 43
- `slide_builder`: 67
- `command_gateway`: 1
- `unknown`: 70 mostly remaining macro-side files outside this scripts-focused stage

## Protected Anchors

Do not move these without a separate staged migration and explicit approval:

- `scripts/RecoilJets_Condor*.sh`
- `scripts/mergeRecoilJets.sh`
- `scripts/recoiljets_cleanup.sh`
- `scripts/make_dstListsData.sh`
- `scripts/makeThesisSimLists.sh`
- `scripts/root_in_analysis_env.sh`
- `scripts/sftp_push_recoiljets.sh`
- `scripts/sftp_get_recoiljets_outputs.sh`
- `scripts/MergeDownloadedRecoilJetsSim.C`, which is called by the local SFTP get helper after SIM downloads
- known transfer-mapped helpers such as `scripts/audit_auau_ml_training_smoke.py`, `scripts/train_auau_photon_bdt.py`, `scripts/train_auau_jet_residual_bdt.py`, and `scripts/validate_auau_tight_bdt_on_sim.py`
- `macros/Fun4All_*`
- `macros/analysis_config*.yaml`
- `macros/AnalyzeRecoilJets*`, `macros/AnalyzeTriggerGroupings.h`, `macros/sPhenixStyle.*`, `macros/Calo_Calib.C`, and `macros/HIJetReco.C`
- `src/` and `src_AuAu/`

## Stage 1 Implemented

- Created `scripts/slides/` as the canonical slide-generator subsystem.
- Moved 61 local slide-generation or slide-prep Python sources into
  `scripts/slides/{auau_bdt,pp_currentian,split_studies,working_point}/`.
- Left top-level compatibility wrappers for every moved slide generator.
- Added `scripts/slides/INDEX.yaml` and `scripts/slides/common/slide_defaults.py`.

## Stage 2 Implemented

- Moved 105 more local non-slide Python helper sources into
  `scripts/{os,diagnostics,data_prep,plotting}/` with top-level wrappers.
- Moved 56 local/offline ROOT macro sources into
  `macros/{plotting,diagnostics}/` with top-level include wrappers.
- Added `scripts/HELPER_INDEX.yaml` and `macros/MACRO_INDEX.yaml`.
- Removed 42 generated front-facing build/cache/platform artifacts
  listed in `agent_context/REPO_ORGANIZATION_CLEANUP_MANIFEST.tsv`.
- Kept protected SDCC/Fun4All/transfer anchors at their historical paths.

## Stage 3 Implemented

- Focused on `scripts/` only.
- Moved the remaining 64 top-level real script sources into canonical homes:
  `scripts/ml/`, `scripts/sdcc/runtime/`, `scripts/sdcc/pipelines/`,
  `scripts/sdcc/workflows/`, `scripts/sdcc/transfer/`,
  `scripts/sdcc/local_merge/`, `scripts/env/`, `scripts/diagnostics/`,
  `scripts/plotting/`, and `scripts/slides/auau_bdt/`.
- Converted 166 prior top-level Python compatibility wrapper files into
  symlinks to canonical source, improving import compatibility and eliminating
  duplicate wrapper code.
- The top-level `scripts/` surface now contains README/index files plus 230
  compatibility symlinks; no real helper source should be added there.
- Added `scripts/SCRIPT_INDEX.yaml` as the exhaustive scripts surface map and
  `scripts/sdcc/SDCC_RUNTIME_INDEX.yaml` as the local SDCC runtime/transfer
  map.
- Updated `scripts/sdcc/transfer/sftp_push_recoiljets.sh` so selected symlink
  paths upload the canonical target contents and `--commit-push` stages both
  compatibility paths and canonical targets.
- No remote SDCC mutation, SFTP upload, deletion, Condor submission, or Google
  Slides mutation was run in stage 3.

## Stage 4 Implemented

- Focused on local and SDCC `scripts/` organization after Justin explicitly
  approved SDCC scripts-folder updates.
- Reduced local top-level `scripts/` aliases from 230 to 88 active
  compatibility symlinks. The remaining top-level regular files are
  `scripts/README.md`, `scripts/SCRIPT_INDEX.yaml`, and
  `scripts/HELPER_INDEX.yaml`.
- Updated `scripts/sdcc/transfer/sftp_push_recoiljets.sh` so script uploads now
  target canonical SDCC paths such as `scripts/sdcc/runtime/`,
  `scripts/sdcc/pipelines/`, `scripts/sdcc/workflows/`, `scripts/ml/`,
  `scripts/diagnostics/`, `scripts/data_prep/`, `scripts/plotting/`,
  `scripts/slides/`, and `scripts/env/`. `RecoilJets_Condor*.sh` remain mapped
  to the SDCC repo base as protected Condor anchors.
- Removed stale static transfer-map entry for generated
  `macros/analysis_config_auau_bdt_mlp_stack_wp080.yaml`; the stack driver
  generates that config from `analysis_config_auau_bdt_mlp_stack_template.yaml`

## Stage 5 Implemented

- Focused on professional systems architecture for local and SDCC `scripts/`.
- Local canonical sources were deepened into domain folders:
  `scripts/sdcc/runtime/{condor,merge,lists,cleanup,xsec,audit}/`,
  `scripts/sdcc/workflows/{submit,target_wp,width_study,stacking,scan,diagnostics}/`,
  `scripts/sdcc/pipelines/{auau,pp}/`,
  `scripts/ml/{training,validation,working_points,stacking,audits}/`,
  `scripts/plotting/{auau_bdt,pp_currentian,stitching,efficiency,trigger,truth_purity}/`,
  `scripts/diagnostics/{pp_shuhang,auau_split,ml_validation}/`,
  `scripts/data_prep/{manifests,ppg12,recoiljets,stitching}/`, and
  `scripts/os/{safety,register,dream,research,context,artifacts}/`.
- Local top-level `scripts/` is now four regular files
  (`README.md`, `SCRIPT_INDEX.yaml`, `HELPER_INDEX.yaml`,
  `COMMAND_INDEX.tsv`), ten directories, and 21 hard aliases. Retired local
  aliases live under `scripts/compat/local/`.
- Added `scripts/bin/thesis-script` and `scripts/COMMAND_INDEX.tsv` so Codex
  can resolve old names, basenames, IDs, and canonical paths without keeping a
  crowded parent directory.
- Added `scripts/sdcc/TRANSFER_MAP.tsv` and
  `scripts/sdcc/validate_transfer_map.py`; regenerated the SFTP push helper
  arrays from canonical paths while preserving retired `scripts/foo` argument
  resolution through `scripts/compat/local/`.
- Mirrored the Stage 5 organization on SDCC for the remote `scripts/` tree
  only. Remote `scripts/` now has `README.md`,
  `SCRIPT_INDEX.remote.tsv`, seven domain directories, and nine hard parent
  aliases; retired aliases live under `scripts/compat/sdcc/`.
- Remote verification passed on `sphnxuser05`: no broken symlinks, no
  top-level trash, `bash -n` passed for representative hard and canonical shell
  scripts, and `python3 -m py_compile` passed for representative ML,
  diagnostic, and slide helpers with external pycache.
  instead of requiring a tracked local upload file.
- Mutated only the SDCC checkout `scripts/` tree plus a local-only SDCC archive
  directory. Remote generated/trash entries archived under
  `agent_context/local/sdcc_scripts_archived/20260601T193900Z/`.
- SDCC verification after migration: remote top-level `scripts/` contains two
  regular files (`README.md`, `SCRIPT_INDEX.remote.tsv`), 77 compatibility
  symlinks, eight top-level directories, no front-facing `__pycache__`/
  AppleDouble/backup trash, and no broken `scripts/` symlinks.
- Representative remote compatibility checks passed: `bash -n` on old and
  canonical SDCC shell paths and `python3 -m py_compile` on representative ML,
  diagnostic, and pp CurrentIAN slide helpers using external pycache.
- No Condor submission, production run, merge rerun, SFTP broad upload, Google
  Slides mutation, or non-scripts SDCC mutation was performed.

## Compatibility Rule

Old paths that appear in policies, register entries, status logs, task-board
notes, active SFTP maps, or likely SDCC/operator command history stay
executable as compatibility symlinks or wrappers. Historical helper paths that
are no longer active contracts are recorded in `scripts/SCRIPT_INDEX.yaml`,
`scripts/HELPER_INDEX.yaml`, and `scripts/slides/INDEX.yaml` instead of being
kept as parent-directory clutter. New local helper source should be created
directly in the appropriate purpose folder.
