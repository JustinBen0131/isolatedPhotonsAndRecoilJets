# Dataset Status

Last updated: 2026-05-06

## Status Labels

- `Valid`: output is trusted for the named use case and has dated evidence.
- `Stale`: output exists but should not be trusted for the named use case.
- `Needs Rerun`: output must be regenerated before use.
- `Running`: rerun has been submitted but fresh output is not yet pulled/validated.
- `Unknown`: not enough evidence has been recorded.

## ABCD/Purity Tight-Axis Rule

- Source fix evidence: local git commit `507eabd3` at `2026-05-05 21:17:40 -0400`.
- Fix location: `src/RecoilJets.cc` and `src_AuAu/RecoilJets_AuAu.cc`.
- The fix changes `fillIsoSSTagCounters(...)` to use the caller-computed `TightTag` from the active configured photon-ID variants.
- SDCC transfer, affected reruns, and fresh `InputFiles/` pulls have not been recorded yet.
- `tightReference/nonTightReference` outputs are lower risk for this specific bug because the old hard-coded ABCD tight axis matched reference behavior.
- Outputs with non-reference tight or non-tight photon-ID axes, including `tightVariantA`, `nonTightVariantA`, renamed/new `newPPG12`-style settings, or later AuAu BDT sideband settings, are `Stale` only for the affected ABCD/purity-derived histogram families if produced before the local fix. The ROOT files are not globally invalid.
- Outputs that only changed preselection while keeping `tightReference/nonTightReference` are lower risk for this specific tight-axis mismatch; do not mark them stale for this bug without additional evidence.

## Affected Histogram Families

These are the histogram families directly affected by the old hard-coded reference tight axis in `fillIsoSSTagCounters(...)`, plus downstream products that read them:

- Candidate ABCD counts: `h_isIsolated_isTight*`, `h_notIsolated_isTight*`, `h_isIsolated_notTight*`, `h_notIsolated_notTight*`.
- Candidate ABCD distributions: `h_Eiso_ABCD_[A-D]*`, `h_pTgamma_ABCD_[A-D]*`.
- ABCD-tagged shower-shape templates: `h_ss_<var>_isIsolated_isTight*`, `h_ss_<var>_notIsolated_isTight*`, `h_ss_<var>_isIsolated_notTight*`, `h_ss_<var>_notIsolated_notTight*`.
- Truth-signal ABCD leakage: `h_sigABCD_MC*`.
- Event-leading xJ purity counts: `h_xJpurityLead_isIsolated_isTight*`, `h_xJpurityLead_notIsolated_isTight*`, `h_xJpurityLead_isIsolated_notTight*`, `h_xJpurityLead_notIsolated_notTight*`.
- Event-leading xJ leakage: `h_xJpurityLead_sigABCD_MC*`.
- Derived outputs: raw/corrected purity, ABCD count tables, leakage factors, purity- or leakage-corrected xJ, and RooUnfold purity-corrected outputs that read the affected families.

## Not Globally Invalid

- Existing pre-fix ROOT files can still be used for histograms, trees, QA, and plots that do not read or derive from the affected families above.
- Before warning that a plot is stale, check the macro/input histogram names. If it does not consume the affected families, this specific bug is not a reason to reject it.
- Generic tight/non-tight QA histograms are not automatically invalidated by this issue unless the plot or correction derives from the affected ABCD/purity families.

## Current Local Evidence

Evidence source: local file timestamps from `InputFiles/` inspected on 2026-05-05.

- `InputFiles/pp24/RecoilJets_pp_ALL_jetMinPt5_7pi_8_vz60_isoR40_fixedIso2GeV_preselectionVariantA_tightVariantA_nonTightVariantA.root`
  - Timestamp: `2026-05-01 10:09:22 -0400`
  - Status for affected ABCD/purity families: `Stale`
- `InputFiles/simPhotonJet/*preselectionVariantA_tightVariantA_nonTightVariantA.root`
  - Timestamps: around `2026-05-04 18:22 -0400`
  - Status for affected ABCD/purity families: `Stale`
- `InputFiles/simEmbedded/*tightVariantA*` and `InputFiles/simEmbedded/*nonTightVariantA*`
  - Timestamps: mostly `2026-05-04 11:18-11:19 -0400`
  - Status for affected ABCD/purity families: `Stale`
- `InputFiles/simEmbedded/merged/*tightVariantA*` and `InputFiles/simEmbedded/merged/*nonTightVariantA*`
  - Timestamps: mostly `2026-05-04 19:09-19:12 -0400`
  - Status for affected ABCD/purity families: `Stale`
- `InputFiles/auau25/RecoilJets_auau_ALL_jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionVariantA_tightReference_nonTightReference.root`
  - Timestamp: `2026-05-05 11:30:07 -0400`
  - Status for this tight-axis mismatch: not marked `Stale` by this bug because the tight/non-tight axes are reference; keep normal provenance checks before final use.

## Regeneration Target For New YAML/View Architecture

Standing goal: reproduce the final local analysis products in these folders
using the current 15 active `photon_id_sets` and the new pool/replay
organization:

- `InputFiles/auau25`
- `InputFiles/pp24`
- `InputFiles/simPhotonJet`
- `InputFiles/simEmbedded`
- `InputFiles/InclusiveJetSIM`
- `InputFiles/InclusiveJetSIM_EMBEDDED`

Expected file-count shape after regeneration:

- pp/data-style products: 15 ROOT files per final dataset output, one per
  active photon-ID working point.
- pp photon-jet SIM and inclusive SIM: 15 ROOT files per SIM sample before
  sample combination; canonical combined products should also collapse to 15
  working-point files when local merge helpers combine samples.
- AuAu-like data/embedded products: 30 ROOT files per dataset/sample family if
  both active `clusterUEpipeline` values remain enabled, equal to 15
  photon-ID working points times `noSub`/`baseVariant`.

After these regenerated products exist, plot-generation work should not assume
the old one-cfg-file-per-slider layout. Plot macros and helper code should
adapt the existing `AnalyzeRecoilJets` plotting/correction structure to resolve
requests through the new ROOT organization: working-point files plus
`AnalysisViewCatalog`/internal view directories for tunable scan axes such as
`coneR`, `jet_pt_min`, `back_to_back_dphi_min_pi_fraction`, `vz_cut_cm`, and
isolation mode.

## Active SDCC Productions

- Embedded-inclusive SIM list generation
  - Status: `Valid` for list-building/matched-list input preparation only.
  - User-visible SDCC terminal evidence from 2026-05-06 showed both
    `run28_embeddedJet12` and `run28_embeddedJet20` completed
    `scripts/makeThesisSimLists.sh` with 10,000 raw rows, 10,000 matched rows,
    and 10,000 pair/triplet rows. `DST_MBD_EPD` is an expected `NONE`
    placeholder for these embedded packs.
  - Current desired `isSimEmbeddedInclusive` sample pair:
    `run28_embeddedJet12` + `run28_embeddedJet20`.
  - Not yet `Valid` for final analysis merged ROOT use: cross-section estimates
    for Jet12/Jet20, propagated constants, RecoilJets ROOT production, local
    pulls, and stitched `embeddedJet12and20merged_SIM` validation still need
    recorded evidence.
- Overnight pool/replay smoke tests
  - Status: `Needs Rerun`.
  - Clusters: `isPP` `3016331`, `isAuAu` `3016333`, `isSim` `3016336`,
    `isSimEmbedded` `3016340`.
  - Gmail evidence showed the `isPP` `20260506_022824` DAG reached its final
    node, but only one pool ROOT was reported from hundreds of planned capture
    jobs, so the smoke result is not a valid production-tuning measurement.
    Local submitter fix has been made so rerun smoke DAGs use args-file queue
    syntax and final emails validate expected vs produced/profiled jobs.
  - Output/report bases are the timestamped `thesisAnaSmoke`,
    `thesisAnaPoolsSmoke`, and `condor_sub/pool_workflow_*_20260506_*`
    directories printed by the submitter.
  - Next evidence needed: queue completion, final parseable stage emails or
    terminal/log evidence, pulled smoke reports via
    `scripts/sftp_get_recoiljets_outputs.sh smokeTestLatest <dataset>`, tuning
    summaries reviewed, ROOT output sanity checked, and SIM comparisons against
    existing local `InputFiles/` variants where configurations overlap.
- AuAu `scaledTriggerStudy`
  - Status: `Valid` for the current presentation plot-generation task based on
    the user's 2026-05-06 report that the scaled-trigger plots are now on slide
    9 of the trigger-analysis deck.
  - Cluster: `1127133` on `sphnxuser04`.
  - Submitted command:
    `RJ_SCALED_TRIGGER_REQUEST_MEMORY=2000MB ./RecoilJets_Condor_submit.sh isAuAu scaledTriggerStudy condorDoAll`.
  - Jobs/runs: `6404` jobs over `620` runs, `groupSize=20`.
  - Output base:
    `/sphenix/tg/tg01/bulk/jbennett/thesisAna/auau/jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionReference_tightReference_nonTightReference_scaledTriggerStudy`.
  - Intended use: scaled trigger efficiency turn-ons from max-cluster-energy
    distributions and trigger-relevant information, with scaledowns applied in
    the merge stage; offline overlays should follow
    `macros/AnalyzeRecoilJets_RunTriggerAna.cpp`.
  - Current follow-up: no further plot-generation task is open for this study.
    Clean and explain slide 9 before the 2026-05-07 group meeting.

## Required Evidence To Mark Affected Outputs Valid

- Source fix transferred to SDCC from local commit `507eabd3` or a later commit containing the same fix.
- `RecoilJets.cc` and `RecoilJets_AuAu.cc` transferred to SDCC.
- Affected Condor jobs rerun after transfer.
- Fresh ROOT outputs pulled into `InputFiles/`.
- This file updated with timestamps, job IDs or pasted SDCC evidence, and any relevant cfg tags.
- Once a variant has fresh post-fix evidence, move it out of `Stale` for the affected families and remove any obsolete broad warning that no longer applies.
