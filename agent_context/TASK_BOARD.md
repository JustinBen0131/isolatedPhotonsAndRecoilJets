# Task Board

Last updated: 2026-05-11

Task tracking is milestone/blocker level. Add tasks when Justin explicitly asks
to track/add/remember them. If a task only appears casually in conversation,
ask before adding it. When a task appears complete, move it to `Done Pending
Removal` and ask Justin before removing or archiving.

## Now

- FUTURE FLAGSHIP R&D: next-generation pp+AuAu photon-finder ML program.
  Justin wants this captured as the long-horizon ambitious goal: after the
  current BDT validation and JSTG slides are stable, build a unified photon-ID
  learning program trained and validated across pp and AuAu so it can become
  the smartest, most physically trustworthy photon finder we can make. Keep it
  separate from the immediate WP0.80 BDT production validation. Proposed
  staged scope:
  1. Establish the current expanded BDT as the transparent baseline with
     efficiency, fake/leakage, isolation, shower-shape, xJ, and centrality/pT
     validation.
  2. Build a pp+AuAu common training table with matched production features,
     explicit domain labels, truth labels, event/candidate weights, and strict
     sample ownership.
  3. Train stronger baselines first: tuned XGBoost/BDT scans, calibrated
     working points, monotonic/regularized variants where physically useful,
     and robust class/phase-space balancing.
  4. Then explore neural-network candidates such as compact MLPs, mixture-of-
     experts by collision system/centrality/pT, and uncertainty/calibration
     heads, always compared against the BDT at fixed signal efficiency and
     fixed background rejection.
  5. Require interpretability before promotion: feature/permutation studies,
     score vs isolation/ABCD safety, domain-shift checks, pp closure, AuAu
     embedded closure, and constrained RecoilJets production validation.
  6. Promote nothing to production unless it beats the BDT in a way that is
     stable, explainable, and useful for PPG19 xJgamma physics.
- ACTIVE ADD-ON: 3x3 shower-width centrality-input AuAu tight-BDT comparison.
  Justin asked to train one model identical in concept to the centrality-input
  BDT but replacing `cluster_weta_cogx/cluster_wphi_cogx` with
  `cluster_weta33_cogx/cluster_wphi33_cogx`, then run it through the same
  `isSimEmbedded` and `isSimEmbeddedInclusive` RecoilJets MC validation
  pipeline. Local implementation target: model id
  `centAsFeat3x3_pt5to40`, output file
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_expanded_20260509_152604/auau_tight_bdt_centAsFeat3x3_pt5to40_tmva.root`,
  cfg tag `preselectionNewPPG12_tightAuAuCentInput3x3BDT_nonTightAuAuBDTComplement_baseVariant`.
  Ordered checklist:
  1. DONE locally 2026-05-10: add trainer spec/features and runtime cfg row.
  2. DONE on SDCC 2026-05-11: uploaded changes, rebuilt `src_AuAu` with
     `make clean; makeProject`, and trained only `centAsFeat3x3_pt5to40` from
     fixed extraction source `auauTightBDT_20260508_233049`.
     Evidence: `sphnxuser01` terminal and Gmail
     `[RecoilJets][auauTightBDT_trainCentInput3x3FromExtraction][READY]`;
     output TMVA:
     `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_expanded_20260509_152604/auau_tight_bdt_centAsFeat3x3_pt5to40_tmva.root`.
  3. DONE on SDCC 2026-05-10/11: isolated one-row runtime/full MC
     validation jobs were submitted from `sphnxuser01` with
     `RJ_PHOTON_ID_ROW_MATCH=AuAuCentInput3x3BDT`, fanout1, isolated WP0.50
     and WP0.80 output roots, and all four final auto-workflow emails reached
     `READY` with `rescue_file_count=0`.
  4. NEXT: pull/pair the four 3x3 one-row outputs from
     `output_bdt3x3_20260510_223529_wp050` and
     `output_bdt3x3_20260510_223529_wp080`, then compare the 3x3-width BDT
     against the standard centrality-input and pT/centrality-binned BDTs.
- ACTIVE RUNBOOK: MC-first RecoilJets validation of expanded AuAu tight-BDT
  families. Goal is to validate the 8 expanded BDT families plus 2 baselines
  in normal RecoilJets embedded-MC outputs, while leaving the main production
  YAML untouched. Campaign config:
  `macros/analysis_config_auau_bdt_validation.yaml`; model source:
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_expanded_20260509_152604`;
  working point: shared AuAu BDT score `> 0.50`; common preselection:
  `newPPG12`; non-tight: BDT complement. Rebuild evidence: on 2026-05-09 from
  `sphnxuser07`, `src_AuAu` finished `make clean; makeProject` successfully
  after the compile fix, installing to
  `/sphenix/u/patsfan753/thesisAnalysis_auau/install`. Ordered checklist to
  track:
  1. DONE 2026-05-09 on `sphnxuser07`: Local signal-side test
     `RJ_CONFIG_YAML=macros/analysis_config_auau_bdt_validation.yaml
     ./RecoilJets_Condor_submit.sh isSimEmbedded local 100
     SAMPLE=run28_embeddedPhoton20 VERBOSE=0` completed all 10 cfg rows,
     including `tightAuAuPtCent3BDT` and `tightAuAuPtCent7BDT`, with
     `exit_code=0` and normal MB-scale ROOT outputs.
  2. SKIPPED BY USER CHOICE 2026-05-09: Local inclusive-background test
     `isSimEmbeddedInclusive local 100 SAMPLE=run28_embeddedJet20`. Residual
     risk is accepted for the overnight full run.
  3. SUPERSEDED BY FULL SUBMISSION: Dry counts:
     `RJ_CONFIG_YAML=macros/analysis_config_auau_bdt_validation.yaml
     ./RecoilJets_Condor_submit.sh isSimEmbedded CHECKJOBS groupSize
     7` and the same for `isSimEmbeddedInclusive`.
  4. SKIPPED BY USER CHOICE 2026-05-09: Smoke `isSimEmbedded`:
     `RJ_NOTIFY_EMAILS=just0131@gmail.com
     RJ_CONFIG_YAML=macros/analysis_config_auau_bdt_validation.yaml
     ./RecoilJets_Condor_submit.sh isSimEmbedded condorDoAllSmoke
     groupSize 7`.
  5. SKIPPED BY USER CHOICE 2026-05-09: Smoke `isSimEmbeddedInclusive` with the same `RJ_CONFIG_YAML` and
     `condorDoAllSmoke groupSize 7`.
  6. DONE / RUNNING 2026-05-09: Full `isSimEmbedded condorDoAll groupSize 7`
     submitted from `sphnxuser07`, cluster `432408`.
  7. DONE / RUNNING 2026-05-09: Full `isSimEmbeddedInclusive condorDoAll
     groupSize 7` submitted from `sphnxuser07`, cluster `432413`.
- OVERNIGHT ACTIVE 2026-05-09: Justin chose to skip smoke after local testing
  and submitted full MC validation campaigns from `sphnxuser07`.
  `isSimEmbedded` command used `RJ_NOTIFY_EMAILS=just0131@gmail.com`,
  `RJ_PROFILE_JOB=1`, `RJ_ID_FANOUT_MAX_ROWS=1`,
  `RJ_REQUEST_MEMORY=4000MB`, `RJ_SIM_FIRSTROUND_REQUEST_MEMORY=6000MB`,
  `RJ_CONFIG_YAML=macros/analysis_config_auau_bdt_validation.yaml`,
  `./RecoilJets_Condor_submit.sh isSimEmbedded condorDoAll groupSize 7`;
  timestamp `simembedded_20260509_234130`, DAG cluster `432408`, DAG path
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembedded_20260509_234130/RecoilJets_auto_simembedded_20260509_234130.dag`.
  `isSimEmbeddedInclusive` was submitted with the same knobs and timestamp
  `simembeddedinclusive_20260509_234610`; terminal evidence showed all 20
  cfg/sample analysis nodes were added, with Jet12 `9998` rows, Jet20 `10000`
  rows, and `1429` groups per sample/cfg; final DAG cluster `432413`; DAG path
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembeddedinclusive_20260509_234610/RecoilJets_auto_simembeddedinclusive_20260509_234610.dag`.
  Gmail check on 2026-05-10: `isSimEmbedded` reached final READY with
  `rescue_file_count=0`, final output base
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/output/simembedded`, and next
  action `./scripts/sftp_get_recoiljets_outputs.sh isSimEmbedded`. The
  inclusive-background campaign has only a current-campaign firstRound STARTED
  email so far; no current `isSimEmbeddedInclusive` READY/CHECK/FAILED email
  was found. Next: request one compact SDCC diagnostic for cluster `432413` or
  the `simembeddedinclusive_20260509_234610` DAG. Once both are READY, pull
  outputs and make the full comparison suite: reco efficiency, isolation
  discrimination/ABCD behavior, tight/complement shower-shape distributions,
  truth leakage/fake rates, and model-family comparisons across centrality and
  photon `E_T`.
- PIVOT FOR JSTG NEXT WEDNESDAY: full AuAu data jobs were removed by Justin on
  2026-05-08. Do not spend Condor space on AuAu data until the simulation/BDT
  story is ready. The priority is a clean 15-minute JSTG presentation draft by
  Monday showing simulation validation and a defensible AuAu photon tight-BDT
  path.
- JSTG priority 1: estimate and validate embedded photon/inclusive cross
  sections. Use `scripts/estimateEmbeddedPhotonXsec.sh` and reproduce the
  existing reference stitching diagnostic
  `dataOutput/combinedSimOnlyEMBEDDED/jetMinPt5_7pi_8_vz60_isoR40_fixedIso2GeV_baseVariant_preselectionReference_tightReference_nonTightReference/photonJet12and20merged_SIM/embeddedPhoton_stitchedTruthFilterPtSpectrum.png`
  in the new compact `preselectionReference_tightReference_nonTightReference_baseVariant`
  pipeline. This must validate the current infrastructure for both
  `isSimEmbedded` and `isSimEmbeddedInclusive`, including Jet12/Jet20 or
  Photon12/Photon20 weights, event counts, and pT spectrum continuity. Active
  xsec estimator pass: on 2026-05-08 from `sphnxuser08`, Justin submitted
  `./scripts/estimateEmbeddedPhotonXsec.sh firstPass --family inclusive
  --xsec-shards 50 --xsec-events 1000000 --mode interpreted` from
  `/sphenix/u/patsfan753/scratch/thesisAnalysis`. Workdir:
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_snapshots/pythia_xsec_20260508_154701`;
  manifest:
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/pythia_xsec_firstPass_20260508_154701.txt`;
  submit file:
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/pythia_xsec_20260508_154701.sub`;
  Condor cluster `631685`; 100 jobs total, 50 shards each for `EmbeddedJet12`
  and `EmbeddedJet20`, 1,000,000 raw events/shard, interpreted mode, 900 MB.
  The 2026-05-08 16:04 secondPass aggregation found all 100 shard CSVs.
  Derived values: `EmbeddedJet12 sigma_eff_pb=1.21692467e+06`
  (`n_pass=482371`, `rel_stat=0.001440`) and
  `EmbeddedJet20 sigma_eff_pb=5.56198698e+04` (`n_pass=532317`,
  `rel_stat=0.001371`). Combined CSV:
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_snapshots/pythia_xsec_20260508_154701/combined_results.csv`;
  summary:
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_snapshots/pythia_xsec_20260508_154701/combined_summary.txt`.
- JSTG priority 2: enhance AuAu tight-BDT training outputs and QA. The training
  package should expose feature distributions, scores, feature importance/model
  internals, centrality-binned validation with finer centrality binning, and
  both pT-dependent and pT-independent checks analogous to the centrality
  treatment. Try minority-class optimization or explicit class-weight scans
  where useful, and record the tradeoff in signal efficiency/background
  rejection.
- JSTG priority 3: validate the BDT clearly in simulation. Produce PPG12-style
  QA that can convince reviewers: ROC curves including centrality bins,
  signal/background score distributions, reco efficiency curves versus photon
  pT and centrality, and comparison to the reference tight ID working point.
  The desired story is a reasonable photon ID across pT/centrality with a
  clear MC efficiency or rejection improvement from BDT tight ID.
- JSTG priority 4: prepare clean, concise 15-minute slides by Monday for draft
  circulation. Slides should cover the three evidence blocks above with enough
  detail to explain the stitching/cross sections, why the BDT is trustworthy,
  and what remains before data promotion.
- Simulation pipeline priority: finish Condor final-stitch support for
  `isSimEmbedded` and `isSimEmbeddedInclusive`, then pull only canonical
  stitched products with `./scripts/sftp_get_recoiljets_outputs.sh
  isSimEmbedded` / `isSimEmbeddedInclusive`. Do not rely on slow local serial
  stitching for the compact internalized files. After finalStitch is `READY`,
  run the strict legacy-vs-compact histogram comparison before treating the new
  files as analysis inputs: old scalar histograms such as
  `h_Eiso_pT_10_12_cent_0_20` from the matching old cfg should compare 1:1 to
  the corresponding new internal-view histogram such as
  `h_Eiso_isoR40_fixedIso4GeV_pT_10_12_cent_0_20` in
  `preselectionReference_tightReference_nonTightReference_baseVariant`, with
  all cuts matched. Pay special attention to `vz`: old scalar files may carry
  `vz60`, while new AuAu-like embedded defaults may use `vz10` unless
  overridden; strict bin-content parity requires the same `vz` acceptance.
- TOP PRIORITY COMPLETE / NEXT PLOTTING: AuAu tight-BDT
  `validateOnSimCondor` succeeded from `sphnxuser07` on 2026-05-08 after the
  sPHENIX/PyROOT worker environment patch.
  Visible terminal evidence shows the command:
  `RJ_NOTIFY_EMAILS=just0131@gmail.com
  RJ_ML_PYTHON=/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python
  ./RecoilJets_Condor_submit.sh isSimEmbeddedAndInclusive trainTightBDT
  validateOnSimCondor
  SOURCE=/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_015049
  MODEL_DIR=/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_20260508_104530
  groupSize 100`. Current run root:
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/condor_sub/auauTightBDTValidate_20260508_154538`;
  report root:
  `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_015049/reports/model_validation_condor_20260508_154538`;
  DAG path
  Gmail reported `[RecoilJets][auauTightBDT_validateOnSimCondor][READY]`;
  email was consumed and marked read. Validation metrics:
  `total_entries=12427173`, `signal_entries=8751401`,
  `background_entries=3675772`, `scored_entries=400000`,
  `finite_score_fraction=1`, AUCs `centINDcontrol=0.881426`,
  `centAsFeat=0.883994`, `centDepBDTs=0.88403`, notes `none`. The report
  includes `validation_metrics.json`, `validation_deep_diagnostics.json`,
  `validation_curves.root`, `validation_auc_table.csv`,
  `validation_threshold_table.csv`, and `validation_feature_summary.csv`.
  Next action: pull/copy these report artifacts and make final centrality/pT
  ROC, AUC, threshold, and feature-summary plots; then run the constrained
  data+MC BDT-variant validation pass before adding the models to broad
  production.
- AuAu data validation remains important, but it is deliberately deprioritized
  until the simulation/BDT JSTG package is ready. Do not restart full AuAu data
  production unless the user explicitly decides to resume it.
- Track the `isSimEmbedded` direct-fanout smoke submitted from `sphnxuser03`
  with `RJ_NOTIFY_EMAILS=just0131@gmail.com RJ_PROFILE_JOB=1
  RJ_SMOKE_SIM_NEVENTS=10000 RJ_SMOKE_SIM_SAMPLE_LIMIT=2
  RJ_SMOKE_REQUEST_MEMORY_EMBEDDED=3000MB ./RecoilJets_Condor_submit.sh
  isSimEmbedded condorDoAllSmoke groupSize 4`. This smoke held on memory:
  clusters `3016478`/`3016479` exceeded the `3072 MB` cgroup limit with
  `RequestMemory=3000`. Remove the held DAG, then rerun with reduced fanout
  (`RJ_ID_FANOUT_MAX_ROWS=5`) and/or higher memory before choosing
  full-production `groupSize` and memory. Follow-up reduced-fanout smoke at
  `RJ_ID_FANOUT_MAX_ROWS=5` also held all `72` jobs at `3000MB`, so the next
  step was a tiny single-sample higher-memory ladder. The `RJ_ID_FANOUT_MAX_ROWS=1`,
  `6000MB`, one-sample smoke reached SIM firstRound and final auto workflow
  `READY`. A follow-up 2026-05-08 full-chain smoke with automatic secondRound
  reached `[RecoilJets][auto_simembedded_final_ready][READY]`; 15 final
  `embeddedPhoton12` smoke ROOT files were pulled under
  `InputFiles/pipelineSmoke/isSimEmbedded/simembedded_smokeTest_20260508_005005_final`
  with nonzero sizes and representative ROOT catalog checks.
- Track the full `isSimEmbedded` production submitted from `sphnxuser03` on
  2026-05-08 after the successful smoke validation. Intended command/settings:
  `RJ_NOTIFY_EMAILS=just0131@gmail.com RJ_PROFILE_JOB=1
  RJ_ID_FANOUT_MAX_ROWS=1 RJ_REQUEST_MEMORY=4000MB RJ_CLEAN_OUTPUT_BASE=1
  ./RecoilJets_Condor_submit.sh isSimEmbedded condorDoAll groupSize 7`.
  Visible terminal evidence confirms `groupSize=7`, `nEvents=0`, both
  `run28_embeddedPhoton12` and `run28_embeddedPhoton20`, timestamp namespace
  `simembedded_condorDoAll_20260508_012157`, one cfg output per fanout pass,
  and `1429` jobs per cfg/sample. Tomorrow before pulling final outputs, move
  or copy the current `InputFiles/simEmbedded` contents into a temporary
  comparison folder such as `InputFiles/crossCheck/simEmbedded_before_20260508`
  so the old files are not overwritten. After pulling, run a strict 1:1
  correctness comparison between the old and new
  `preselectionReference_tightReference_nonTightReference_baseVariant` outputs
  for the same sample/cut toggles: compare ROOT object paths, classes, binning,
  duplicate/unique histogram names, representative bin contents/integrals, and
  content uniqueness across suffixed internal-view histograms for every
  relevant histogram family. This must explicitly catch over-tagging: if a
  quantity should be canonical because iso cone, iso WP, jet-pT threshold, or
  dphi cut does not physically affect it, the output should not contain
  multiple differently suffixed histograms with identical contents. Conversely,
  histograms whose physics meaning really does vary by an internal axis should
  have distinct, correctly suffixed names and non-accidentally duplicated
  contents unless low statistics explain equality. Expected result, if the new direct
  fanout/internal-view implementation is correct: the matching reference/
  reference/reference legacy-toggle histograms map exactly to the previous
  files, except for clearly intentional differences from internalized
  iso/cone, jet-pT, dphi, or metadata naming. Any difference must be written
  down with the histogram name, old/new location, and whether it is intended or
  a bug.
- Track the fixed full `isSimEmbeddedInclusive` production submitted from
  `sphnxuser08` on 2026-05-08 after the embedded-inclusive stitching fix was
  uploaded and `src_AuAu` was rebuilt. The old `sphnxuser01` campaign
  `simembeddedinclusive_condorDoAll_20260508_112828` / DAGMan cluster
  `2670772` is superseded and later produced a stale `[FAILED]` email with
  `rescue_file_count=1`; do not use its outputs. Current fixed command from
  `/sphenix/u/patsfan753/scratch/thesisAnalysis`:
  `RJ_NOTIFY_EMAILS=just0131@gmail.com RJ_PROFILE_JOB=1
  RJ_ID_FANOUT_MAX_ROWS=1 RJ_REQUEST_MEMORY=4000MB
  RJ_SIM_FIRSTROUND_REQUEST_MEMORY=6000MB RJ_CLEAN_OUTPUT_BASE=1
  ./RecoilJets_Condor_submit.sh isSimEmbeddedInclusive condorDoAll groupSize 7`.
  Current timestamp namespace:
  `simembeddedinclusive_condorDoAll_20260508_163858`; auto DAG:
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembeddedinclusive_20260508_163858/RecoilJets_auto_simembeddedinclusive_20260508_163858.dag`;
  snapshot dir:
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_snapshots/simembeddedinclusive_20260508_163858`;
  YAML/fanout dir:
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_yaml_overrides/simembeddedinclusive_condorDoAll_20260508_163858`;
  bulk output base:
  `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembeddedinclusive`.
  Terminal evidence shows production `vz=[10]`, `run28_embeddedJet12` has
  `9998` matched rows, `run28_embeddedJet20` has `10000` matched rows, and
  each visible sample/cfg fanout plans `1429` analysis jobs at `groupSize=7`.
  Heartbeat check at 2026-05-08 17:31 EDT found visible Condor worker cluster
  `631694` with `4714` jobs total, `3117` running, `1597` idle, and `0` held.
  Gmail on 2026-05-08 22:20 EDT reported
  `[RecoilJets][auto_simembeddedinclusive_final_ready][READY]` from
  `sphnxuser08` for the `163858` auto DAG, with `rescue_file_count=0`, final
  output base `/sphenix/u/patsfan753/scratch/thesisAnalysis/output/simembeddedinclusive`,
  and next action `./scripts/sftp_get_recoiljets_outputs.sh
  isSimEmbeddedInclusive`; the email was consumed and marked read. The
  workflow completed `analysis -> SIM_FIRSTROUND -> SIM_SECONDROUND ->
  SIM_FINALSTITCH`; the `finalStitch` Condor stage is where the weighted
  Jet12/Jet20 stitch is applied with `EmbeddedJet12
  sigma_eff_pb=1.21692467e+06` and `EmbeddedJet20 sigma_eff_pb=5.56198698e+04`.
  Slide 5 in
  `https://docs.google.com/presentation/d/1-bLijV9dHONak_7WV9aFAC_eZXf7LLytaqs74bo2NKQ/edit?slide=id.g3dfaa525a20_5_2#slide=id.g3dfaa525a20_5_2`
  is waiting on the reference/reference/reference final stitched output:
  `preselectionReference_tightReference_nonTightReference_baseVariant/embeddedJet12and20merged_SIM/RecoilJets_embeddedJet12plus20_MERGED.root`.
  Next: pull or validate that output with
  `./scripts/sftp_get_recoiljets_outputs.sh isSimEmbeddedInclusive`,
  regenerate the inclusive stitched truth-filter pT QA plot, show the PNG in
  chat for approval, then replace the old photon plot on slide 5.
- Validate the fixed SIM automatic merge chain. On 2026-05-08,
  `scripts/RecoilJets_Condor_submit.sh` was updated and uploaded so SIM-family
  auto workflows run `analysis -> SIM_FIRSTROUND -> SIM_SECONDROUND -> final
  READY/CHECK email`, matching the data-side expectation. Run an SDCC dry-run
  and then a tiny real smoke before full `isSim`, `isSimInclusive`,
  `isSimEmbedded`, or `isSimEmbeddedInclusive` production.
- AuAu tight-BDT sidecar extraction is smoke-validated and ready for full
  production. The 2026-05-08 `sphnxuser02` smoke produced 24 ROOT files with
  nonzero signal/background (`entries=37320 signal=26316 background=11004`) and
  trained/read back all three model families. The notifier environment bug was
  patched/uploaded in `scripts/auau_tight_bdt_pipeline.sh`. Next action is full
  extraction over embeddedPhoton12/20 and embeddedJet12/20, then
  `trainFromExtraction` and `applyCheck` on the full model output before adding
  model paths to `analysis_config.yaml`. After the full models exist, read and
  sanity-check the TMVA outputs/summaries, add the three AuAu tight-BDT model
  variants to `macros/analysis_config.yaml`, then run a fixed quick data+MC
  validation configuration only over those new BDT variants. The purpose of
  that quick pass is to compare output behavior against the reference tight ID
  before promoting the new model variants into broader production.
- Next production-validation pass should be simulation-first: `isSimEmbedded`,
  `isSimEmbeddedInclusive`, then `isSim`/`isSimInclusive` only as needed for
  pp-like reference/background cross-checks relevant to the BDT and stitching
  presentation.
- After the next embedded signal/background passes, derive z-vertex and
  centrality reweighting and compare two AuAu photon-BDT model families:
  reweighted and non-reweighted.
- QA the Region-C purity-correction path in depth. Confirm candidate-level
  `h_isIsolated_*` purity summaries are not confused with event-leading
  `h_xJpurityLead_*` counters used for xJ purity correction and per-photon
  normalization.
- Analyze centrality-dependent efficiencies, not only inclusive/all-centrality
  trigger and photon-efficiency summaries.

## Next

- Regenerate final local analysis inputs under `InputFiles/auau25`,
  `InputFiles/pp24`, `InputFiles/simPhotonJet`, `InputFiles/simEmbedded`,
  `InputFiles/InclusiveJetSIM`, and `InputFiles/InclusiveJetSIM_EMBEDDED`.
- Rerun pp isolation-efficiency fits and retune sliding isolation windows for
  pp/AuAu and cone-specific `R=0.3`/`R=0.4` behavior.
- Reproduce pp photon-pT reweighting / Figure-30-like comparison and propagate
  the response-matrix weighting policy.
- Fix and validate the pp inclusive-jet SIM stitching path.
- Embedded-inclusive Jet12/Jet20 cross-section estimation from Condor cluster
  `631685` on `sphnxuser08` aggregated cleanly with all 100 shard CSVs. Derived
  constants are `EmbeddedJet12 sigma_eff_pb=1.21692467e+06` and
  `EmbeddedJet20 sigma_eff_pb=5.56198698e+04`. Propagate constants only after
  sanity checks, produce/pull RecoilJets outputs, and validate
  `embeddedJet12and20merged_SIM`.
- Before locking the AuAu photon-BDT feature set, study `wEta`/`wPhi` variants
  in `src/PhotonClusterBuilder.cc`, compare pp vs AuAu behavior/definitions,
  and decide which should be stored/used as BDT features.
- Ensure future ML/tree products store cheap scalar cluster information needed
  for photon-ID, BDT/NPB, timing, HCAL leakage, and shower-shape studies.

## Waiting

- Full `isAuAu` production from `sphnxuser04`: monitor DAGMan cluster `1134652`,
  the downstream analysis worker cluster once it appears, READY/CHECK/FAILED
  emails, and held-job state. The older `sphnxuser01` 3500MB campaign and
  `sphnxuser07` two-run smoke are superseded by this full submission unless
  later evidence shows they still need cleanup.
- `isSimEmbedded` smoke recovery from `sphnxuser03`: remove held clusters
  `3016477 3016478 3016479` and `3016487 3016488 3016489 3016490 3016491
  3016492 3016493`; pull the successful `RJ_ID_FANOUT_MAX_ROWS=1`, `6000MB`
  smoke reports and use them for final memory/groupSize/fanout decision.
- Full `isSimEmbedded` production from `sphnxuser03`: wait for READY/CHECK
  emails for analysis/firstRound/secondRound/final, then preserve old
  `InputFiles/simEmbedded` before pulling and do the strict old-vs-new ROOT
  comparison before using the outputs as final analysis inputs.
- Fixed full `isSimEmbeddedInclusive` production from `sphnxuser08`: remote
  `163858` auto DAG was READY with `rescue_file_count=0`, and the old
  `112828`/`2670772` run is superseded. On 2026-05-08 Codex pulled 15 final
  merged files with `./scripts/sftp_get_recoiljets_outputs.sh
  isSimEmbeddedInclusive`; the reference/reference/reference target exists at
  `dataOutput/combinedSimOnlyEMBEDDED/preselectionReference_tightReference_nonTightReference_baseVariant/embeddedJet12and20merged_SIM/RecoilJets_embeddedJet12plus20_MERGED.root`.
  ROOT inspection found `SIM/MERGE_INFO` with Jet12/Jet20 derived constants and
  `SIM/h_embedInclusiveStitch_filterJetPt_kept`, then generated
  `embeddedInclusiveJet_finalMergedStitchedTruthJetPtSpectrum.png` next to the
  ROOT file. Codex then pulled the raw per-sample secondRound files with
  `SFTP_GET_SIM_RAW=1 ./scripts/sftp_get_recoiljets_outputs.sh
  isSimEmbeddedInclusive` and generated the slide-style split component plot
  `embeddedInclusiveJet_stitchedTruthFilterPtSpectrum.png` from the Jet12 and
  Jet20 raw inputs. The helper's slow automatic all-config local merge was
  stopped after the raw files and plot existed. Next: approve the split PNG in
  chat, then replace the plot on slide 5.
- Fresh SDCC smoke or production evidence for direct-fanout pp and SIM jobs.
- User decision on when to transfer the ABCD/purity tight-axis fix and rerun
  affected outputs.
- User decision on final `isSimInclusive` pp inclusive-jet stitching
  prescription.

## Done Pending Removal

- Scaled-trigger efficiency plot generation for the current presentation pass
  was reported complete; remaining value is slide explanation/polish and any
  future follow-up plots. Ask Justin before archiving.

## Archived

- 2026-05-05: Identified ABCD/purity tight-axis fix in local commit `507eabd3`.
- 2026-05-05: Created durable Codex tracking notes for project state, dataset
  status, known issues, and run logs.
- 2026-05-06: Embedded-inclusive list generation for `run28_embeddedJet12` and
  `run28_embeddedJet20` reported complete with 10,000 matched rows each.

## Automation Backlog

- Campaign-scoped email-first production watchdog for `RecoilJets Pipeline`
  Gmail label every 30 minutes during active production. Keep it paused by
  default; activate only after submitted/running job evidence or an explicit
  user request to watch a campaign; pause again after outputs are pulled/
  validated, failure is diagnosed, or Justin says to stop watching. Summarize
  failures/ready outputs and provide exact next commands without mutating SDCC.
- Lightweight Markdown status dashboard updates from pipeline emails, visible
  terminal output, pulled files, and local ROOT inspection.
- Local read-only QA summaries for purity/xJ normalization, stitched SIM
  spectra, and output validity/provenance.
- Plot-generation helpers that reuse existing `AnalyzeRecoilJets` style and
  inspect representative PNGs before slide integration.
