# Project Board

Last updated: 2026-05-06

## Now

- Pipeline infrastructure/speed update is active.
- AuAu scaled-trigger efficiency study is running on SDCC from user-reported
  submission cluster `1127133` on `sphnxuser04`. This is a special
  `scaledTriggerStudy` production intended to make scaled trigger turn-on
  overlays from max-cluster-energy information, with scaledowns applied during
  merge.
- ABCD/purity tight-axis fix is present in local git commit `507eabd3` in `src/RecoilJets.cc` and `src_AuAu/RecoilJets_AuAu.cc`, but affected ABCD/purity-derived histogram families are not valid until the fix is transferred to SDCC and rerun.
- Treat current `InputFiles/` ROOT products using non-reference tight/non-tight photon-ID axes such as `tightVariantA`, `nonTightVariantA`, or `newPPG12`-style settings as stale only for affected ABCD/purity-derived histograms and downstream corrections unless later evidence says otherwise. The files are still usable for unrelated histograms and QA.

## Next

- Standing production target: regenerate the local final analysis inputs under
  `InputFiles/auau25`, `InputFiles/pp24`, `InputFiles/simPhotonJet`,
  `InputFiles/simEmbedded`, `InputFiles/InclusiveJetSIM`, and
  `InputFiles/InclusiveJetSIM_EMBEDDED` with the new pool/replay YAML/view
  architecture.
- Expected regenerated output shape:
  - pp/data-style final outputs: 15 ROOT files per dataset output, one per
    active `photon_id_sets` working point.
  - pp photon-jet SIM and inclusive SIM: 15 ROOT files per SIM sample before
    canonical sample combination, then 15 combined/canonical products when the
    local merge helper materializes them.
  - AuAu-like outputs with both active `clusterUEpipeline` modes: 30 ROOT files
    per dataset/sample family, equal to 15 photon-ID working points times
    `noSub`/`baseVariant`, unless a production intentionally restricts the UE
    axis.
- After these regenerated inputs are available, all new plot/regeneration work
  should treat the 15/30 ROOT files plus `AnalysisViewCatalog`/internal view
  directories as the primary organization. Plot code should mutate/adapt from
  the existing `AnalyzeRecoilJets` pipeline's plotting structure, mathematical
  corrections, and functionality so user requests can resolve semantic choices
  like photon-ID working point, coneR, jet pT cut, dphi, vz, isolation mode, and
  UE pipeline without relying on the old one-file-per-cut-variant layout.
- Transfer the source fix to SDCC when ready:
  `./scripts/sftp_push_recoiljets.sh RecoilJets.cc RecoilJets_AuAu.cc`
- Rerun affected pp, embedded SIM, and AuAu jobs that used non-reference tight/non-tight photon-ID variants.
- Pull fresh outputs into `InputFiles/` and update `codex_notes/DATASET_STATUS.md` with file timestamps and job evidence.

## Waiting On SDCC

- `scaledTriggerStudy` AuAu production cluster `1127133` should be monitored
  until it leaves the queue. The prompt line in the pasted terminal showed
  `condor_q 1127132`, but the scaled-trigger submission itself reported
  `6404 job(s) submitted to cluster 1127133`; use `1127133` for this task.
- When cluster `1127133` is done: run the AuAu merge stages for the
  `_scaledTriggerStudy` output, pull the single final AuAu ROOT file locally,
  then make scaled-trigger efficiency turn-on overlays using the logic in
  `macros/AnalyzeRecoilJets_RunTriggerAna.cpp`.
- Future Condor status should be recorded from user-pasted SDCC terminal output or approved clipboard handoff commands.

## Needs User Decision

- Decide which affected variant set gets rerun first after the pipeline infrastructure/speed update is ready.
- Decide when the local ABCD/purity tight-axis fix is ready to transfer to SDCC and rerun.

## Done Recently

- 2026-05-05: Identified that the ABCD/purity tight-axis fix is in local commit `507eabd3` and later than current affected `InputFiles/` products.
- 2026-05-05: Created durable Codex tracking notes for project state, dataset status, known issues, and run logs.
