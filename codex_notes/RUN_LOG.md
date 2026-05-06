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
