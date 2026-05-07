# Run Log

Use this file for short dated notes about local inspections, user-reported SDCC/Condor activity, pulls, merges, plot-generation passes, and status changes.

## 2026-05-07

- Priority decision: stop broad redesign and focus on getting the pipeline to
  production quickly. The code is believed to be close, but not yet proven; the
  remaining risk is practical SDCC/pipeline bugs such as wrapper environment
  mismatches, empty pool outputs, missing replay branches/features, strict final
  summary checks, or merge/pull mapping issues. The next useful information
  should come from smoke outputs, not more speculative architecture work.
- Immediate execution priority is embedded SIM first, because it unlocks the
  AuAu ML tight-ID path the user cares about most:
  `isSimEmbedded` for embedded photon signal samples and
  `isSimEmbeddedInclusive` for embedded inclusive-jet/background samples. Use
  quick 1k-event smoke DAGs first to test list resolution, pool capture,
  replay fanout, output ROOT production, emails, and reports:
  `RJ_NOTIFY_EMAILS=just0131@gmail.com RJ_SMOKE_SIM_NEVENTS=1000
  RJ_SMOKE_SIM_MAX_JOBS_PER_SAMPLE=1 ./RecoilJets_Condor_submit.sh
  isSimEmbedded condorDoAllSmoke`, then the same command for
  `isSimEmbeddedInclusive`.
- If the embedded smoke DAGs finish cleanly, proceed aggressively to real
  embedded SIM production and only do the minimum additional smoke testing
  needed for pp/AuAu DATA before full submissions. The intended near-term order
  is: finish `isSimEmbedded` smoke, finish `isSimEmbeddedInclusive` smoke, fix
  only concrete failures, submit real embedded SIM production, run quick pp/AuAu
  smoke checks, then submit full pp/AuAu passes. Do not spend many more days
  tuning before getting full passes over data moving.
- Clarified current `sphnxuser02` status from visible terminal output:
  `isSimEmbedded condorDoAllSmoke` was building the DAG, not yet running
  worker jobs. The output showed matched 10,000-line lists for
  `run28_embeddedPhoton12` and `run28_embeddedPhoton20`, split into
  `groupSize=2` chunks and capped to one capture job/sample for the quick smoke.
  It then built replay fanout over 2,160 view rows, 15 working-point output
  roots, and 15 replay shards per capture axis/sample. This is the intended
  mini full-pipeline shape: capture small DST chunks to pools, replay into
  cfg-tagged smoke ROOT outputs, then final summary/email.

## 2026-05-06

- User reported that the AuAu scaled-trigger efficiency plots are now placed
  on slide 9 of the trigger-analysis Google Slides deck
  (`1yrtM1Xxyb-uSrSyQqAbLOwUFlg9gTE41F4i0Uyr7Sbk`, slide
  `id.g3dd90ada7d5_0_161`). Scaled-trigger plot generation is therefore no
  longer an open todo for the current presentation pass. Remaining priority:
  clean the slide formatting and write/explain the interpretation by
  2026-05-07 11:00 EDT if possible, and definitely before the 14:00 EDT group
  meeting.
- User-visible SDCC terminal output showed `scripts/makeThesisSimLists.sh`
  completed the embedded-inclusive list packs for both current target samples:
  `run28_embeddedJet12` and `run28_embeddedJet20`. Each sample reported 10,000
  raw rows for `DST_CALO_CLUSTER`, `G4Hits`, `DST_JETS`, `DST_GLOBAL`, and the
  placeholder `DST_MBD_EPD`; 10,000 matched rows for the same lists; and 10,000
  rows in all pair/triplet convenience lists. The terminal also reported
  embedded completeness as `COMPLETE (10000/10000 OutDirs have
  calo+G4+truth-jet+global)`. Treat `run28_embeddedJet12` +
  `run28_embeddedJet20` as the active `isSimEmbeddedInclusive` sample pair.
  Remaining work before final-use promotion: run the inclusive Jet12/Jet20
  cross-section estimator, propagate constants, run/pull RecoilJets outputs,
  and validate `embeddedJet12and20merged_SIM` stitching.
- Gmail inspection of RecoilJets stage emails showed the `isPP`
  `dataPoolSmoke_pp_20260506_022824` DAG reached its final node, but the
  payload was not a valid 10-run smoke validation: it reported only
  `pool_root_files=1` despite `588` planned capture jobs, and the profile
  summary was contaminated by an older failed `20260505` profile row because
  profile globs did not include the workflow timestamp. Local submitter fix:
  `scripts/RecoilJets_Condor_submit.sh` now writes Condor worker arguments to
  explicit `.args` files and uses `queue arguments from ...`, avoiding repeated
  `arguments/queue` blocks inside DAG submit files; pool/DAG profile prefixes
  now include the workflow timestamp; final stage emails now include
  `expected_capture_jobs`, `expected_replay_jobs`, `profile_rows`, and
  `profile_failures`, and mark missing pool/profile evidence as `FAILED` or
  `CHECK` instead of reporting a misleading `READY`.
- User submitted the overnight smoke-test DAGs from `sphnxuser03` after the
  pool fanout arithmetic-return fix was pushed. Terminal evidence showed all
  four dataset commands reached DAG submission:
  `isPP` DAGMan cluster `3016331`,
  `isAuAu` DAGMan cluster `3016333`,
  `isSim` DAGMan cluster `3016336`,
  and `isSimEmbedded` DAGMan cluster `3016340`.
  The submitted shapes were: pp `588` capture jobs and `90` replay jobs; AuAu
  two capture/replay axes (`noSub`, `baseVariant`) with `216+216` capture jobs
  and `20+20` replay jobs; pp SIM three photon-jet samples with `1429` capture
  jobs and `72` replay jobs per sample; embedded SIM two embedded photon
  samples times two UE axes with `1429` capture jobs and `72` replay jobs per
  sample/axis. The visible `condor_q` after submission showed the DAGMan
  clusters and initial child nodes queued with no user-held jobs. Next morning:
  check `sphnxuser03`, pull smoke reports, tune `groupSize`/memory from
  `tuning_inputs.json`, inspect smoke ROOT outputs, compare matching SIM
  configurations to existing local `InputFiles/` variants, and decide whether
  full-SIM smoke products can be renamed/promoted rather than rerun.
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
