# Run Log

Use this file for short dated notes about local inspections, user-reported SDCC/Condor activity, pulls, merges, plot-generation passes, and status changes.

## 2026-05-06

- Added standing regeneration target: reproduce final local analysis inputs in
  `InputFiles/auau25`, `InputFiles/pp24`, `InputFiles/simPhotonJet`,
  `InputFiles/simEmbedded`, `InputFiles/InclusiveJetSIM`, and
  `InputFiles/InclusiveJetSIM_EMBEDDED` using the current 15 active
  `photon_id_sets` and new pool/replay internal-view architecture.
- Recorded expected output organization: 15 ROOT working-point files per
  pp-style final dataset/sample output, 15 times SIM samples before canonical
  combination, and 30 ROOT files for AuAu-like outputs when both
  `clusterUEpipeline` modes remain active.
- Recorded plotting implication: once these regenerated inputs are available,
  future plot generation should adapt from the existing `AnalyzeRecoilJets`
  pipeline to resolve semantic requests through working-point files plus
  `AnalysisViewCatalog`/internal view directories, while preserving the
  mathematical correction and plotting functionality the user relies on.
- Added running task for the non-embedded pp inclusive-jet SIM path:
  `isSimInclusive` / `InputFiles/InclusiveJetSIM` should not remain a
  one-sample `jet5` placeholder if jet5/jet10/jet20 are available. Need to
  inspect the current PPG12 inclusive-jet stitching/weighting prescription and
  implement a canonical stitched product analogous in spirit to the pp
  PhotonJet5/10/20 combination, so later comparisons can choose among
  inclusive-jet stitching configurations cleanly.
- Added prerequisite task for embedded inclusive-jet samples: after embedded
  inclusive `jet12` production finishes, generate/download the `jet12` and
  `jet20` outputs, run the cross-section estimator over that sample set, then
  implement and validate proper stitching before using the products in any
  final analysis code. This should mirror the care used for `isSimEmbedded`
  stitching rather than treating the files as immediately analysis-ready.
- Added two high-priority post-regeneration tasks. First, after all new final
  ROOT inputs are available in the target `InputFiles/` folders, rerun the pp
  isolation-efficiency fit workflow beginning with
  `dataOutput/pp/ppg12Style_isoCutEfficiencyFits_ppMerged_fixedIso2GeV_reference.png`
  and use it to tune/validate the pp and AuAu sliding-window configuration for
  the relevant `isoR=0.30`, `isoR=0.40`, and photon-ID variants before rerunning
  affected `isSliding` outputs from the correct pipeline stage. Second,
  regenerate the Figure-30-like pp photon pT reweighting comparison beginning
  with
  `dataOutput/combinedSimOnly/pp_reference_vs_variantA_unfolding_overlay/figure30LikePhotonReweighting/reference_leadingPhoton_figure30_like_purityCorrectedData_vs_pythiaSignal.png`,
  then implement PPG12-style pT reweighting for pp response matrices and decide
  the analogous independently derived AuAu treatment before final response use.
- Added future analysis-development tasks to pursue after a general analysis
  cut set is pinned down: pi0/diphoton distribution studies with fitting and
  diphoton cluster tagging, plus a later UE-subtraction scan that exercises
  existing `clusterUEpipeline` options and adds more finely spaced
  cluster-level UE-subtraction variants. These should be focused quick reruns
  around the chosen baseline, not another broad all-cuts production pass.
- Added high-priority embedded-inclusive ML follow-up: once embedded inclusive
  `jet12` production exists, run `scripts/makeThesisSimLists.sh`, sanity-check
  the `jet12` file count against embedded inclusive `jet20`, then train the
  available ML modes including photon BDT variants, JetML/residual variants,
  and NPB. NPB must be gated by cluster-timing QA first, checking whether the
  timing behavior is similar enough to the pp NPB use case before trusting NPB
  labels or model outputs.
- Recorded pool-schema direction from the user's workflow: for the next broad
  DST pass, it is worth over-storing cheap scalar photon/cluster information so
  future feature/cut/model scans do not require reopening DSTs. Inspection of
  `PhotonClusterBuilder.cc` and `RawClusterBuilderTemplate.cc` shows useful
  scalar candidates beyond the current pool subset: raw energy-window scalars
  (`e11`, `e22`, `e33`, `e55`, `e77`, `e13`, `e15`, `e17`, `e31`, `e51`,
  `e71`, `e35`, `e37`, `e53`, `e73`, `e57`, `e75`, `e32`, `e52`, `e72`),
  widths/CoG/position helpers (`weta`, `wphi`, `weta_cog`, `wphi_cog`,
  `detamax`, `dphimax`, `detacog`, `dphicog`, `drad`, raw/corrected tower CoG
  where accessible), HCAL leakage/nearest-tower scalars (`ihcal_et`,
  `ohcal_et`, `ihcal_et22`, `ohcal_et22`, `ihcal_et33`, `ohcal_et33`,
  `ihcal_ieta`, `ihcal_iphi`, `ohcal_ieta`, `ohcal_iphi`), extra isolation
  scalars (`iso_02_emcal`, `iso_01_emcal`, `iso_005_emcal`,
  `ppg12_iso_axis_eta`, `ppg12_iso_axis_phi`), and RawClusterBuilderTemplate
  cluster metadata (`energy`, `ecore`, `prob`, `chi2/ndf`, mean time,
  incidence alpha phi/eta, tower thresholds, peak threshold, profile/noise
  settings, input/output node names, tower-info mode, tower selection flag,
  detailed-geometry flag, alternate-vertex mode, and subcluster splitting).
  Do not store full tower maps by default; the intended payload is scalar-rich,
  not DST-like.

## 2026-05-05

- Created durable Codex tracking notes in `codex_notes/`.
- Recorded ABCD/purity tight-axis mismatch status: fix exists locally in `src/RecoilJets.cc` and `src_AuAu/RecoilJets_AuAu.cc`, but affected outputs need SDCC transfer, rerun, pull, and validation before being trusted.
- Local timestamp inspection found current affected `InputFiles/` products were produced before the local fix time of approximately `2026-05-05 20:08 EDT`.
- No affected rerun Condor submissions have been recorded yet.
- User reported an active AuAu `scaledTriggerStudy` Condor production submitted
  from `sphnxuser04` with:
  `RJ_SCALED_TRIGGER_REQUEST_MEMORY=2000MB ./RecoilJets_Condor_submit.sh isAuAu scaledTriggerStudy condorDoAll`.
  Submission details from pasted terminal output:
  cluster `1127133`, `6404` jobs, `620` runs, `groupSize=20`,
  request memory `2000MB`, snapshot
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_snapshots/auau_20260504_131105`,
  submit file `RecoilJets_auau_20260504_131117.sub`, output base
  `/sphenix/tg/tg01/bulk/jbennett/thesisAna/auau/jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionReference_tightReference_nonTightReference_scaledTriggerStudy`.
  Note: the terminal prompt later showed `condor_q 1127132`, but this
  specific scaled-trigger production is cluster `1127133`.
- Purpose of cluster `1127133`: produce scaled-trigger efficiency turn-on
  inputs using only trigger-relevant information, especially max cluster energy
  distributions, then apply scaledowns during merging and make offline trigger
  overlay plots following `macros/AnalyzeRecoilJets_RunTriggerAna.cpp`.
  Follow-up when done: merge per-run AuAu outputs, perform sliceRuns/final
  addChunks as needed, pull the single final AuAu file locally, then generate
  scaled-trigger turn-on overlays.
