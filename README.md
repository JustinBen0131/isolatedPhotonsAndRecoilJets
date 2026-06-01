# Isolated Photons And Recoil Jets

Analysis code for sPHENIX isolated-photon and recoil-jet studies. The
repository contains the C++ analysis modules, Fun4All steering macros, ROOT
analysis macros, Condor/workflow helpers, plotting utilities, and ML helpers
used to produce and study pp, Au+Au, photon+jet simulation, embedded
photon+jet simulation, and embedded inclusive-jet background samples.

The code is organized so production entrypoints stay stable while local
analysis helpers live in purpose-specific folders.

## Repository Layout

- `src/`: pp-style `RecoilJets` module and photon-cluster helper code.
- `src_AuAu/`: Au+Au and embedded-sample `RecoilJets_AuAu` module.
- `macros/`: ROOT and Fun4All macro surface.
  - `Fun4All_*.C`, `analysis_config*.yaml`, `Calo_Calib.C`,
    `HIJetReco.C`, `AnalyzeRecoilJets*`, and `sPhenixStyle.*` are stable
    top-level analysis/runtime anchors.
  - `macros/plotting/`: offline ROOT plotting macros, grouped by analysis
    area such as `auau_bdt`, `target_wp`, `width_study`, `stitching`,
    `pp_currentian`, and `ssqa`.
  - `macros/diagnostics/`: ROOT inspection and validation helpers.
  - `macros/MACRO_INDEX.yaml` and `macros/bin/thesis-macro` resolve macro
    names to canonical implementation paths.
- `scripts/`: command entrypoints plus organized helper subsystems.
  - `scripts/sdcc/`: production/runtime/workflow helpers and transfer maps.
  - `scripts/ml/`: training, validation, working-point, and stacking helpers.
  - `scripts/plotting/`: Python plotting utilities.
  - `scripts/slides/`: full-slide PNG builders and slide-specific assets.
  - `scripts/diagnostics/`: row-contract, split-study, ML, and workflow audits.
  - `scripts/data_prep/`: manifest, table, stitching, and compact-data helpers.
  - `scripts/os/`: project operating-system, register, guard, and maintenance
    utilities.
  - `scripts/bin/thesis-script` resolves script ids/basenames to canonical
    paths.
- `agent_context/`: project-local coordination, policies, indexes, and status
  ledgers used to keep multi-step analysis work reproducible.

## Main Analysis Flow

1. Build DST or simulation input lists.
   - Data list helpers live under `scripts/sdcc/runtime/lists/`.
   - The compatibility commands `scripts/make_dstListsData.sh` and
     `scripts/makeThesisSimLists.sh` remain available.

2. Configure photon, isolation, jet, and dataset variants.
   - Main configuration starts from `macros/analysis_config.yaml`.
   - Campaign-specific configs live beside it as `analysis_config*.yaml`.

3. Run Fun4All/RecoilJets.
   - pp-style workflows use `macros/Fun4All_recoilJets.C`.
   - Au+Au and embedded workflows use `macros/Fun4All_recoilJets_AuAu.C`.
   - Shared steering logic lives in `macros/Fun4All_recoilJets_unified_impl.C`.

4. Submit or stage production workflows.
   - Condor wrappers remain available as top-level compatibility commands:
     `scripts/RecoilJets_Condor_submit.sh`,
     `scripts/RecoilJets_Condor.sh`, and
     `scripts/RecoilJets_Condor_AuAu.sh`.
   - Canonical workflow implementations live under `scripts/sdcc/`.

5. Merge and pull analysis-ready outputs.
   - `scripts/mergeRecoilJets.sh` handles production-side merging.
   - `scripts/sftp_get_recoiljets_outputs.sh` pulls merged ROOT outputs into
     local analysis input areas.

6. Produce downstream studies.
   - ROOT analysis and QA use the `AnalyzeRecoilJets*` macro family and the
     organized `macros/plotting/` / `macros/diagnostics/` helpers.
   - ML studies use `scripts/ml/`.
   - Plot and slide candidates use `scripts/plotting/` and `scripts/slides/`.

## Dataset Modes

The production and analysis helpers are built around these dataset modes:

- `isPP`: pp data-style RecoilJets production.
- `isAuAu`: Au+Au data-style RecoilJets production.
- `isSim`: pp photon+jet simulation.
- `isSimInclusive`: pp inclusive-jet simulation/background.
- `isSimEmbedded`: Au+Au embedded photon+jet simulation.
- `isSimEmbeddedInclusive`: Au+Au embedded inclusive-jet background.

## Finding Code

Prefer the indexes and resolver commands before adding a new helper:

```bash
scripts/bin/thesis-script path <script-id-or-basename>
macros/bin/thesis-macro path <macro-id-or-basename>
```

Common indexes:

- `scripts/SCRIPT_INDEX.yaml`
- `scripts/HELPER_INDEX.yaml`
- `scripts/slides/INDEX.yaml`
- `scripts/sdcc/SDCC_RUNTIME_INDEX.yaml`
- `scripts/sdcc/IO_RUNTIME_INDEX.yaml`
- `macros/MACRO_INDEX.yaml`

New code should go directly into the most specific existing folder rather than
adding new flat files to the top of `scripts/` or `macros/`.

## ROOT And Environment

ROOT-dependent commands should use the repository wrapper so the analysis
environment is configured consistently:

```bash
./scripts/root_in_analysis_env.sh root -l -q 'macros/MyMacro.C()'
```

Use external cache/build locations for syntax checks and ACLiC products. Do
not commit generated ROOT dictionaries, `.so` files, `.pcm` files, Python
`__pycache__`, or local output products.

## Outputs And Tracked Scope

This repository tracks source code, macros, scripts, indexes, and lightweight
project documentation. It does not track large ROOT inputs, generated analysis
outputs, plot batches, slide decks, copied external repositories, local
environment state, or private machine-specific configuration.

Analysis-ready local ROOT files should be kept outside the tracked source tree
or in ignored data/input areas. Compact plots, QA tables, and slide assets
should be regenerated from the tracked scripts and macros when possible.
