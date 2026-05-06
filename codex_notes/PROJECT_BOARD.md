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
- HIGH PRIORITY after the regenerated local final inputs are present in every
  target folder: rerun the pp isolation-efficiency fit workflow and use it to
  tune/validate the sliding-window configuration for pp and AuAu. The concrete
  reference output to reproduce/update first is
  `dataOutput/pp/ppg12Style_isoCutEfficiencyFits_ppMerged_fixedIso2GeV_reference.png`.
  The follow-up should cover the relevant pp `isoR=0.30`/`isoR=0.40`
  variants and any working-point variant whose photon selection changes the
  derived sliding-window behavior. After the correct pp/AuAu sliding-window
  configuration is established, rerun all affected `isSliding` variants from
  the earliest required pipeline stage, not from scratch unless the stored pool
  variables are insufficient.
- HIGH PRIORITY after the regenerated local final inputs are present: reproduce
  the pp photon pT reweighting / Figure-30-like comparison and implement the
  corresponding simulation response-matrix weighting policy. The concrete
  reference output to regenerate first is
  `dataOutput/combinedSimOnly/pp_reference_vs_variantA_unfolding_overlay/figure30LikePhotonReweighting/reference_leadingPhoton_figure30_like_purityCorrectedData_vs_pythiaSignal.png`.
  Match the PPG12-style pp pT reweighting behavior, apply it independently to
  the relevant pp simulation response matrices, then decide and implement the
  analogous AuAu treatment so AuAu response matrices receive the correct
  independently derived weighting where needed.
- Fix and validate the non-embedded pp inclusive-jet SIM path
  (`isSimInclusive` / `InputFiles/InclusiveJetSIM`) before relying on it for
  final comparisons. The intended task is to stitch/combine the inclusive jet
  samples in the same spirit as the pp PhotonJet5/10/20 stitching path, using
  the proper PPG12-style sample boundaries/weights where applicable, instead
  of treating only one `jet5` file as the canonical inclusive product.
  Preserve the user's ability to choose the final configuration later, but
  make the offline products organized so the inclusive-jet stitching choice can
  be compared cleanly.
- Embedded inclusive-jet list generation for the current desired samples is
  complete on SDCC: terminal evidence from 2026-05-06 shows
  `run28_embeddedJet12` and `run28_embeddedJet20` each produced 10,000 raw,
  matched, pair, and triplet list rows, with `DST_MBD_EPD` correctly represented
  by a 10,000-line `NONE` placeholder. All `isSimEmbeddedInclusive` analysis
  purposes should now use `run28_embeddedJet12` + `run28_embeddedJet20`, not the
  old `run28_embeddedJet10` + `run28_embeddedJet20` pair.
- HIGH PRIORITY now: run the embedded-inclusive cross-section estimator for
  Jet12 and Jet20, copy the resulting `sigma_eff_pb` values into
  `macros/AnalyzeRecoilJets.h`, then generate/pull/merge the
  `isSimEmbeddedInclusive` RecoilJets outputs and validate the
  `embeddedJet12and20merged_SIM` stitched products before using them in final
  analysis. After that, proceed to the AuAu/embedded ML training path. Train the
  available ML modes, including photon tight-BDT variants, JetML/residual
  correction variants, and the NPB variant. For the NPB variant, first run
  cluster-timing QA and confirm the timing behavior is plausibly analogous to
  the pp NPB usage before treating NPB training labels or outputs as
  analysis-ready.
- Before the next broad pool production, extend the photon pool schema to
  intentionally over-store cheap scalar cluster information from
  `PhotonClusterBuilder` and `RawClusterBuilderTemplate`. The goal is to make
  future photon-ID cuts, BDT/NPB feature studies, cluster-timing QA, HCAL-leakage
  checks, and shower-shape scans replay-only whenever the needed information is
  already computed in the first DST pass. Store scalar values, not raw tower
  maps: all PhotonClusterBuilder shower-shape scalars, isolation pieces,
  timing scalars, HCAL nearest-tower/leakage scalars, raw/corrected CoG style
  variables, cluster probability/chi2/ecore/mean-time/incidence-angle metadata,
  and enough cluster-builder configuration metadata to know how the clusters
  were made. Avoid storing full tower maps or DST-like constituent arrays by
  default unless a later task explicitly needs them.
- Transfer the source fix to SDCC when ready:
  `./scripts/sftp_push_recoiljets.sh RecoilJets.cc RecoilJets_AuAu.cc`
- Rerun affected pp, embedded SIM, and AuAu jobs that used non-reference tight/non-tight photon-ID variants.
- Pull fresh outputs into `InputFiles/` and update `codex_notes/DATASET_STATUS.md` with file timestamps and job evidence.

## Waiting On SDCC

- Overnight smoke-test DAGs were submitted on `sphnxuser03` on 2026-05-06.
  Check them first on the next session before any full production decision:
  `isPP` cluster `3016331`, `isAuAu` cluster `3016333`, `isSim` cluster
  `3016336`, and `isSimEmbedded` cluster `3016340`. The visible `condor_q`
  showed the DAGMan jobs plus first child nodes in the queue with no held jobs
  for the user's jobs at submission time. When these finish, pull smoke reports
  with `scripts/sftp_get_recoiljets_outputs.sh smokeTestLatest <dataset>` for
  `isPP`, `isAuAu`, `isSim`, and `isSimEmbedded`; use the profiling summaries
  and `tuning_inputs.json` to choose full-production `groupSize` and
  `request_memory` for each dataset. Then inspect the smoke ROOT outputs enough
  to confirm all expected working-point/view information exists, compare SIM
  output against existing local `InputFiles/` variants for matching
  configurations where possible, and only then decide whether to promote the
  full-SIM smoke outputs as final products and move pp/AuAu to full submission.
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
- Decide/confirm the desired `isSimInclusive` pp inclusive-jet stitching
  prescription after inspecting current PPG12 treatment and the available
  jet5/jet10/jet20 samples, then implement it in the local merge/plot path.
- Do not promote embedded inclusive `jet12`/`jet20` products to final-use
  inputs until the new cross-section estimator has been run, the resulting
  constants have been propagated, and the ROOT production/merge/stitching
  prescription has been implemented/validated.
- Do not promote embedded-inclusive ML models trained from the new `jet12`/`jet20`
  sample set until file-count/list sanity checks, cluster-timing QA for NPB,
  and basic training/application QA have been recorded.

## Future Analysis Development

- After a general analysis cut set is pinned down, analyze pi0/diphoton
  background structure in more depth: generate pi0 mass and related diphoton
  distributions, fit them, and implement/use diphoton cluster tagging so the
  background sample composition can be understood beyond the current photon-ID
  sideband machinery. This should be a focused rerun/scan around the chosen
  baseline cuts, not another broad all-cuts production pass.
- After a general analysis cut set is pinned down and the main production
  pipeline is stable, explore additional UE-subtraction variations. This
  includes trying cluster-level UE-subtracted inputs, exercising all currently
  available `clusterUEpipeline` variants, and adding more finely spaced
  UE-subtraction variations with smaller controlled changes between adjacent
  configurations. The point is to make these reruns quick and interpretable by
  anchoring them to the chosen baseline cuts.

## Done Recently

- 2026-05-05: Identified that the ABCD/purity tight-axis fix is in local commit `507eabd3` and later than current affected `InputFiles/` products.
- 2026-05-05: Created durable Codex tracking notes for project state, dataset status, known issues, and run logs.
