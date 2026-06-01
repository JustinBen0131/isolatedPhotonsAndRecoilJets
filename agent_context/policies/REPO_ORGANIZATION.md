# Repo Organization Policy

This policy controls where Codex should put new helper code in
ThesisAnalysis. The goal is an obvious, navigable repo where a new collaborator
can tell what is SDCC pipeline code, ROOT/Fun4All entrypoint code, plotting
code, slide-candidate code, ML/validation code, and Codex OS machinery.

## Core Rule

Codex must be organization-conscious by default. Before creating a new script,
macro, notebook-like helper, slide-generation helper, plotting utility, or
diagnostic one-off, choose the narrowest existing or planned home that makes
the file's purpose obvious.

Do not keep adding unrelated side helpers to the top level of `scripts/` or
`macros/` unless the file is a protected entrypoint that must remain there for
SDCC, ROOT, transfer, or backward-compatibility reasons.

## Protected Compatibility Anchors

Do not move, rename, or repurpose SDCC/Fun4All/base pipeline entrypoints unless
Justin explicitly asks for a coordinated refactor and a compatibility plan has
been written first.

Protected compatibility anchors include:

- `scripts/RecoilJets_Condor*.sh`
- `scripts/mergeRecoilJets.sh`
- `scripts/recoiljets_cleanup.sh`
- `scripts/make_dstListsData.sh`
- `scripts/makeThesisSimLists.sh`
- `scripts/root_in_analysis_env.sh`
- `scripts/sftp_push_recoiljets.sh`
- `scripts/sftp_get_recoiljets_outputs.sh`
- `macros/Fun4All_recoilJets*.C`
- `macros/Fun4All_auauTightBDTTraining.C`
- `macros/analysis_config*.yaml`
- `src/`, `src_AuAu/`, and any file loaded by the SDCC build/runtime
  environment at a fixed path.

These command paths may still be edited when the requested work requires it,
but the paths are part of the operating contract. After THE-23 stage 5, local
and SDCC `scripts/` parents are minimal command surfaces: canonical source
lives under deeper purpose folders, hard aliases stay front-facing only for
active runtime/operator contracts, and historical local aliases live under
`scripts/compat/local/`. Edit the canonical source, keep active compatibility
paths available, and use the SFTP transfer map for any SDCC upload.

After THE-23 stage 8, local `macros/` follows the same contract. Protected
Fun4All, analysis-config, style, and `AnalyzeRecoilJets*` runtime anchors stay
at the top level. Offline plotting and diagnostic implementation source lives
under deeper domain folders such as `macros/plotting/auau_bdt/`,
`macros/plotting/stitching/`, `macros/plotting/width_study/`,
`macros/diagnostics/stitching/`, and `macros/diagnostics/preselection/`.
Only actively referenced macro hard aliases remain at top level; lower-risk
old wrappers live under `macros/compat/local/`. Use
`macros/bin/thesis-macro path <id-or-basename>` and `macros/MACRO_INDEX.yaml`
before editing or adding offline ROOT macros.

Only active OS/SDCC/operator command contracts should remain as top-level
`scripts/` symlinks. Historical helper aliases that are not active contracts
belong in `scripts/compat/local/`, `scripts/SCRIPT_INDEX.yaml`,
`scripts/HELPER_INDEX.yaml`, or `scripts/slides/INDEX.yaml`, not as extra
parent-directory entries. Use `scripts/bin/thesis-script path <id-or-name>` to
resolve the canonical home.

## Side-Code Placement

For new local-side code, prefer a purpose-specific subfolder once that folder
exists. If it does not exist, create it when the new file is clearly part of
that category and the move cannot affect protected SDCC entrypoints.

Recommended `scripts/` categories:

- `scripts/os/{safety,register,dream,research,context,artifacts}/`: Codex OS,
  register, doctor, dream, Linear/daily-cockpit, context-pack, and bookkeeping
  utilities.
- `scripts/sdcc/runtime/{condor,merge,lists,cleanup,xsec,audit}/`: fixed-path
  SDCC runtime scripts.
- `scripts/sdcc/workflows/{submit,target_wp,width_study,stacking,scan,diagnostics}/`:
  submit, merge-staging, queue-gated, and campaign drivers.
- `scripts/sdcc/pipelines/{auau,pp}/`: AuAu/pp pipeline entrypoints.
- `scripts/plotting/{auau_bdt,pp_currentian,stitching,efficiency,trigger,truth_purity}/`:
  local plotting utilities.
- `scripts/slides/<campaign>/<family>/`: full-slide PNG builders,
  deck-specific slide candidate generators, and slide-story assembly helpers.
- `scripts/ml/{training,validation,working_points,stacking,audits}/`:
  BDT/MLP/logreg/stack training, validation, scoring, and comparison helpers.
- `scripts/diagnostics/{pp_shuhang,auau_split,row_contracts,sdcc_audit,ml_validation}/`:
  focused audits, equivalence checks, row-contract checks, and short-lived
  investigative tools.
- `scripts/data_prep/{manifests,ppg12,recoiljets,stitching}/`: manifest
  builders, CSV/JSON extractors, compact-output converters, and data shaping.

Recommended `macros/` categories:

- `macros/plotting/{auau_bdt,target_wp,width_study,stitching,pp_currentian,ssqa}/`:
  ROOT plotting macros that are local/offline only.
- `macros/diagnostics/{auau_bdt,stitching,preselection,ssqa,model_checks,feasibility}/`:
  ROOT inspection, validation, audit, and exploratory macros.
- `macros/compat/local/`: lower-risk old macro wrappers retired from the
  parent command surface.
- `macros/bin/`: macro resolver/gateway helpers.
- `macros/fun4all/`: future Fun4All entrypoints only after a compatibility
  migration plan; current Fun4All anchors stay at existing paths.
- `macros/config/`: future analysis YAML templates only after confirming no
  SDCC or helper path expects the current top-level location.
- `macros/styles/`: shared ROOT style helpers if/when split from current
  top-level style anchors.

Top-level `scripts/` should not receive new source files. Use it only for
active compatibility symlinks, README/index files, and explicitly approved
command contracts.
Top-level `macros/` should not receive new offline plotting or diagnostic
source files. Use it only for protected runtime/config/style anchors, README
and macro indexes, and documented hard aliases.

## Before Creating A New File

1. Search first with `rg --files scripts macros` and reuse or extend an
   existing helper if that is the obvious path.
2. Decide whether the file is a protected pipeline/Fun4All/SDCC entrypoint or
   a local-side helper.
3. If protected, keep the expected path and update transfer maps or SDCC
   upload instructions only when needed.
4. If local-side helper, place it in the most specific folder available, or
   create the folder when safe.
5. Name the file by purpose and campaign. Avoid vague names like
   `make_plot.py`, `new_slide.py`, or `test.C`.
6. Write outputs under a matching campaign/dataOutput folder, not beside the
   source code.
7. Record any new reusable convention in the relevant policy or README when it
   prevents future confusion.

## Before Moving Existing Files

Do not move existing files just to clean up while active SDCC jobs, slide
deadlines, or collaborator-facing outputs depend on them.

Before any move:

1. Inventory references with `rg` across policies, scripts, macros, submit
   wrappers, docs, and notebooks.
2. Classify each target as protected anchor, reusable side helper, stale
   one-off, or unknown.
3. Write a migration map from old path to new path.
4. Add compatibility wrappers or symlinks when old paths are likely referenced
   by SDCC, docs, commands, or muscle memory.
5. For shell and Python paths that are imported or invoked both locally and on
   SDCC, prefer symlinks over Python wrappers so imports and command execution
   still load the canonical source.
6. Run syntax checks and any direct smoke tests.
7. For SDCC-mapped files, follow `TRANSFER_AND_PIPELINE_FILES.md`. Since
   THE-23 stage 5, `scripts/sdcc/TRANSFER_MAP.tsv` records local canonical
   path, local alias, SDCC canonical path, SDCC compatibility path, and upload
   target. Historical remote `scripts/foo` command paths remain only where
   they are hard contracts.

## Slide And Plot Helper Rule

New slide/plot generators should make their audience and scope obvious from
the path:

- reusable plotting primitive: `scripts/plotting/...`
- full-slide PNG builder: `scripts/slides/...`
- one-campaign diagnostic plot: `scripts/diagnostics/...` or
  `scripts/slides/<campaign>...` if it is slide-facing
- ROOT offline plot macro: `macros/plotting/...` unless it must stay top-level
  for an existing wrapper

When in doubt, prefer a new clearly named subfolder plus a short README over
another top-level helper file.

## Refactor Safety

Organization is valuable only if it does not break the analysis. Preserve
current placement for anything SDCC, Fun4All, Condor, transfer, or ROOT wrapper
infrastructure expects at a fixed path. Make the repo more obvious by default,
but never at the cost of silent path breakage.
