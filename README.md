# ThesisAnalysis RecoilJets Pipeline

This repository contains the core source, macros, and scripts for the sPHENIX
RecoilJets analysis pipeline. The tracked code is organized around producing
RecoilJets ROOT outputs from pp, Au+Au, photon+jet simulation, embedded
photon+jet simulation, and embedded inclusive-jet simulation samples, then
building the downstream QA, stitching, purity, response, and unfolding products.

## Code Layout

- `src/`: pp-style RecoilJets module and photon-cluster helper code.
- `src_AuAu/`: Au+Au and embedded RecoilJets module.
- `macros/`: Fun4All steering macros, analysis configuration, ROOT QA,
  stitching, plotting, photon-purity, response, and unfolding workflows.
- `scripts/`: dataset list builders, Condor submission wrappers, merge helpers,
  transfer helpers, environment wrappers, diagnostics, and ML training scripts.

## Main Production Flow

1. Build or update DST/input lists.
   - Data lists are handled by `scripts/make_dstListsData.sh`.
   - Thesis simulation lists are handled by `scripts/makeThesisSimLists.sh`.

2. Configure analysis settings.
   - The main analysis configuration lives in `macros/analysis_config.yaml`.
   - This controls photon-ID working points, isolation/preselection variants,
     jet/recoil settings, Au+Au cluster-UE modes, and related analysis axes.

3. Run the Fun4All analysis module.
   - pp-style workflows use `macros/Fun4All_recoilJets.C`.
   - Au+Au and embedded workflows use `macros/Fun4All_recoilJets_AuAu.C`.
   - Shared implementation details live in
     `macros/Fun4All_recoilJets_unified_impl.C`.

4. Submit or stage production on Condor.
   - The main entry point is `scripts/RecoilJets_Condor_submit.sh`.
   - Dataset-specific execution is handled through
     `scripts/RecoilJets_Condor.sh` and
     `scripts/RecoilJets_Condor_AuAu.sh`.

5. Merge and collect outputs.
   - `scripts/mergeRecoilJets.sh` handles production-side merging.
   - `scripts/sftp_get_recoiljets_outputs.sh` pulls merged outputs into the
     local analysis input area.
   - `scripts/MergeDownloadedRecoilJetsSim.C` builds canonical local SIM
     combinations after downloads.

6. Run downstream analysis products.
   - QA, stitching, photon-purity, response-matrix, unfolding, trigger, and
     comparison studies are implemented in the `macros/AnalyzeRecoilJets*`
     macro family and focused helper macros.

## ROOT Environment

ROOT-dependent commands should be run through the repository environment
wrapper:

```bash
./scripts/root_in_analysis_env.sh /Users/patsfan753/Desktop/analysis/env/bin/root -l -q 'macros/MyMacro.C()'
```

This keeps ROOT, RooUnfold, and analysis-library paths consistent across macro
runs and ACLiC builds.

## Dataset Modes

The submission and macro paths are structured around these main dataset modes:

- `isPP`: pp data-style RecoilJets production.
- `isAuAu`: Au+Au data-style RecoilJets production.
- `isSim`: pp photon+jet simulation.
- `isSimInclusive`: pp inclusive-jet simulation/background.
- `isSimEmbedded`: Au+Au embedded photon+jet simulation.
- `isSimEmbeddedInclusive`: Au+Au embedded inclusive-jet
  simulation/background.

Embedded-inclusive analysis is intended to use the current Jet12 and Jet20
sample pair once the corresponding cross-section weights and stitched products
are validated.

## Important Analysis Products

- Photon-ID and isolation QA are produced from the RecoilJets module outputs
  and the `AnalyzeRecoilJets` macro family.
- Photon-purity and Region-C subtraction use ABCD-style isolated/tight sideband
  histograms, with separate handling for event-leading xJ purity counters.
- xJ unfolding uses photon-yield unfolding for the `N_gamma` normalization and
  two-dimensional response matrices in `(pT_gamma, xJgamma)`.
- Embedded photon and embedded inclusive samples feed the Au+Au photon-BDT and
  ML-response studies.
- z-vertex and centrality reweighting should be derived and compared for
  embedded signal/background BDT training.

## GitHub Scope

The GitHub repository intentionally tracks only source, macros, scripts, and
this README. Local ROOT inputs, generated outputs, presentations, reference
PDFs, scratch files, external dependencies, copied reference repositories, and
machine-local state are not part of the tracked repository.
