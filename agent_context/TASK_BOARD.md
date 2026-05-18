# Task Board

Last updated: 2026-05-18

Task tracking is milestone/blocker level. Add tasks when Justin explicitly asks
to track/add/remember them. If a task only appears casually in conversation,
ask before adding it. When a task appears complete, move it to `Done Pending
Removal` and ask Justin before removing or archiving.

## Now

- TOP PRIORITY / JSTG SLIDES BY MONDAY EOD: build the 15-minute Jet Structure
  Topical Group story around one simple, defensible model:
  `globalEtCent1535_bdt_noIso`, the no-isolation global BDT trained on
  Photon12+20 signal and Jet12+20+30 inclusive background over
  `15 < E_T < 35 GeV`. Main message: the BDT uses shower-shape /
  energy-sharing evidence plus cluster `E_T` and centrality context; isolation
  stays out of the tight-ID model so the ABCD purity logic remains defensible.
  Required main-flow slides: motivation for AuAu-trained photon ID, model
  definition, sample/label sanity with Jet30 background necessity, score
  separation/AUC, WP80 operating point / fake-rate and purity behavior, and a
  production-style RecoilJets check from the completed noIso WP80 merged signal
  and background outputs. Variants such as isolation-input BDT, MLP, stackers,
  routed/fine-bin BDTs, and EtFine studies should appear only as one compact
  end/backup slide unless Justin explicitly pivots the story.
  Evidence already in hand: BDT validation is READY for stamp
  `20260516_135439`; noIso AUC `0.840908`, finite score fraction `0.999994`,
  signal mean score `0.670902`, background mean score `0.302054`, and WP80
  threshold `0.5619116425514221` for the original flat-threshold cross-check.
  Corrected production output now exists for the preferred `8 E_T x 3
  centrality` grid2d WP80 cuts plus cone-specific sliding-isolation fix: tag
  `global_etcent_inclusive3_bdt_noIso_cent3grid_isoFix_20260517_2030` on
  `sphnxuser02` has final MERGED ROOTs for signal Photon12+20 and
  inclusive-background Jet12+20+30. Optional backup tag
  `global_etcent_inclusive3_bdt_iso_wp80_20260517_001911` is also merged, but
  should not headline the JSTG talk.
- TOP PRIORITY / CORRECT WP80 RUNTIME CUTS: derive the no-isolation BDT WP80
  thresholds as `E_T`-dependent cuts within centrality bins from the existing
  validation score caches. Show Justin two PNGs before rerunning production:
  one with `3` centrality bins (`0-20`, `20-50`, `50-80`) and one with `7`
  centrality bins (`0-10`, `10-20`, `20-30`, `30-40`, `40-50`, `50-60`,
  `60-80`). Moving forward, default AuAu BDT production should use these
  centrality-binned `grid2d` WP80 cuts, not the prior single flat threshold
  `0.5619116425514221`, unless Justin explicitly asks for a flat-threshold
  cross-check.
- ACTIVE SIDECAR / ETCENT-BINNED SIXPACK-FEATURE BDT COMPARISON: this
  user-requested backup/model comparison lane is READY after a finalize-only
  recovery. Four groups were trained from the current sixpack Photon12+20 plus
  Jet12+20+30 staging sample: noIso `8x3`, noIso `8x7`, iso `8x3`, and iso
  `8x7`. Output root:
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/global_etcent_inclusive3_sixpack_20260516_135439/bdt_binned_sidecars/etcent_binned_sixpack_features_20260517_174151`.
  Important invariant: the only intended difference from the global BDTs is
  the independent `E_T x centrality` binning; every binned model still keeps
  `cluster_Et` and `centrality` as input features and uses the same noIso/iso
  feature families as the validated global sixpack BDTs. Evidence from
  `sphnxuser02`: original DAG cluster `2122802` and cache child `2122803` left
  the queue with `0` held jobs; all `160` models were present. The DAG
  `FINALIZE` node failed with status `127` from a missing ML-venv
  `libpython3.13.so.1.0`, so Codex did not retrain and instead merged the
  registries directly. Final summary now reports `status=READY`,
  `xgb_json_count=160`, `metadata_json_count=160`, `tmva_root_count=160`,
  `status_by_model={'trained': 160}`, and `missing_paths=0`; merged registry is
  `bdt_models/model_registry.json`. Compact cleanup evidence was preserved at
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/cleanup_preserved_summaries/etcent_binned_sixpack_features_20260517_174151`;
  the official `312M` merged registry was not duplicated. Justin then asked to
  validate the sidecar and compare it to the global BDT. Current baseline-only
  validation is READY from `sphnxuser02`: `8x3` DAGMan cluster `2123919` wrote
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/global_etcent_inclusive3_sixpack_20260516_135439/validation/bdt_binned_noIso_cent3`,
  and `8x7` DAGMan cluster `2123938` wrote
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/global_etcent_inclusive3_sixpack_20260516_135439/validation/bdt_binned_noIso_cent7`.
  Coarse-centrality AUCs are global BDT `0.811/0.838/0.854`, binned `8x3`
  `0.826/0.846/0.863`, and binned `8x7` `0.838/0.860/0.871`.
  Compact PNG/CSV artifacts were pulled locally and visually checked:
  `/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439/slideReady/binned_bdt_comparison/cent7/global_noiso_bdt_vs_binned_noiso_bdt_ptCent7_score_separation_by_centrality.png`
  and backup
  `/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439/slideReady/binned_bdt_comparison/cent3/global_noiso_bdt_vs_binned_noiso_bdt_ptCent3_score_separation_by_centrality.png`.
  Next checkpoint: decide whether the binned-BDT comparison belongs in backup
  slides or should remain an internal model-selection note. Do not delete it as
  obsolete.
- TOP PRIORITY / PP CROSS-CHECK: run the no-centrality-as-feature BDT
  validation in pp, using the combined SIM photon-jet `5+10+20` signal sample
  and the corresponding inclusive-jet background counterpart. Goal: check the
  same BDT variant logic outside AuAu/centrality so the pp baseline can anchor
  the AuAu no-centrality comparison cleanly. Required outputs: validation AUC /
  ROC, WP80-style fake-rate and threshold diagnostics, score-separation plots,
  and a concise comparison against the normal pp BDT/reference photon-ID
  baseline.
- TOP PRIORITY / FINISH TODAY: validate truth signal/background tagging behind
  the BDT/MLP score-separation plots. Blair explicitly asked whether the signal
  not peaking at 1 could indicate truth-tagging issues. Check event-level
  signal-tag multiplicity in signal embedded MC and verify the inclusive
  embedded sample has only a rare truth-signal tag rate, roughly at the
  accidental/prompt level rather than order-one per event. Required outputs:
  per-event tagged-photon multiplicity histograms/tables for signal and
  inclusive samples, duplicate matched-truth checks for any signal event with
  more than one signal-tagged reco photon, and a short conclusion on whether
  the original BDT score-separation plot is stable after this sanity check.
- TODO / BDT CLOSURE QA: show signal/background BDT-score separation on the
  training sample for the current best BDT variants, then compare it directly
  with the independent embedded-validation score separation. Goal: closure-test
  whether the BDT learned a stable photon/background ordering or whether the
  training sample looks artificially sharper than validation. Required outputs:
  training-sample score overlays, validation-sample score overlays using the
  same binning/style, AUC or separation summary table, and a short overtraining
  conclusion.
- DONE PENDING REMOVAL / OFFLINE: target80 BDT all-available campaign
  `bdt_target80_gated_20260512_001012` on `sphnxuser05`. Codex fixed the
  remaining EtFine merge issue: `scripts/mergeRecoilJets.sh` now recognizes
  `auauEtFineCentInputBDT`, `auauEtFineCent3BDT`, and
  `auauEtFineCent7BDT`, and
  `scripts/merge_auau_bdt_target80_ready.sh` passes `MERGE_CFG_MATCH` through
  to the merge script. Local `bash -n` and `git diff --check` passed; both
  scripts were uploaded to SDCC. Remote dry scan with `MERGE_CFG_MATCH=EtFine`
  reports `analysis_config_etfine_15to35_target80` as `READY` with
  `finals=6/12`, while all other target80 configs are already complete
  (`expanded_5to40=22/22`, `widthstudy_pt10to35=8/8`,
  `widthstudy_pt1530=8/8`, `widthstudy_pt15to35=8/8`,
  `widthstudy_pt5to35=8/8`). Active repair tmux is
  `target80_merge_repair_etfine_20260513_1233` with
  `MERGE_CFG_MATCH=EtFine`, `RJ_TARGET80_MERGE_DO_RUN=1`,
  `RJ_TARGET80_MERGE_LOOP=0`, and `RJ_TARGET80_MERGE_MAX_CONFIGS=1`. Launch
  evidence: EtFine cfg tags were found, first-round merge submitted worker
  cluster `5348890.*` and DAGMan cluster `5348891`, with `80` active/idle jobs
  and `0` held. 2026-05-13 12:55 EDT update: the local SFTP pull script was
  also patched to normalize `auauCentInputBase3x3BDT`, EtFine BDT modes, and
  current MLP/stack mode names consistently with the merge script. All
  already-complete non-EtFine target80 finals were pulled offline to
  `dataOutput/target80_all_available/bdt_target80_gated_20260512_001012/`;
  local counts are now EtFine `6/12`, expanded `22/22`, and all four
  width-study configs `8/8` for a total of `60` nonzero final ROOT files
  (`~11 GB`). Latest remote check showed EtFine `isSimEmbedded firstRound`
  still running with `114` merge jobs and `0` held. 13:36 EDT update: EtFine
  advanced to `9/12` remote final stitched files and is now in
  `isSimEmbeddedInclusive firstRound` with `57` running merge jobs and `0`
  held; Gmail had no unread pipeline messages. Next checkpoint: wait for
  EtFine `12/12` final stitched files, then pull the remaining EtFine finals
  offline.
  13:40 EDT gym-mode automation update: heartbeat
  `watch-gated-target80-campaign` now checks every 15 minutes while Justin is
  away. It is authorized to do targeted helper/script fixes, upload those
  fixes, continue/restart only the already-authorized EtFine merge repair if a
  targeted fix requires it, and pull final ROOT outputs offline. It is not
  authorized to remove jobs, delete remote outputs, or submit unrelated new
  production campaigns. Completion means all remote final counts are present
  (`EtFine=12/12`, expanded `22/22`, all width-study configs `8/8`), local
  nonzero final ROOT files are pulled under
  `dataOutput/target80_all_available/bdt_target80_gated_20260512_001012/`, and
  an inventory table is written/reported with config, expected count, remote
  count, local count, and local path.
  17:14 EDT final check: SDCC `sphnxuser05` has `0` Justin jobs and `0` held
  jobs for this lane. Remote final counts are complete:
  `EtFine=12/12`, `expanded_5to40=22/22`, `widthstudy_pt10to35=8/8`,
  `widthstudy_pt1530=8/8`, `widthstudy_pt15to35=8/8`, and
  `widthstudy_pt5to35=8/8`. Local target80 folder contains all `66` expected
  nonzero final ROOT files (`~12 GB`) under
  `dataOutput/target80_all_available/bdt_target80_gated_20260512_001012/`.
  Heartbeat `watch-gated-target80-campaign` was deleted because the watch task
  is complete. Next: generate target80 ID efficiency, reco efficiency,
  background tight-rate/fake-rate, shower-shape, isolation/ABCD, and xJ
  comparison plots from this local matrix.
- READY FOR FIRST TARGET80 PLOTS: first completed target80 BDT
  signal/background MC pair is offline locally. Folder:
  `/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/target80_first_offline/bdt_target80_gated_20260512_001012/analysis_config_etfine_15to35_target80`.
  It contains 6 nonzero finalStitch ROOT files for config
  `analysis_config_etfine_15to35_target80`: signal `isSimEmbedded` and
  inclusive-background `isSimEmbeddedInclusive` for `reference`, `newPPG12`,
  and `auauCentInputBase3x3BDT`. Next plotting priority: target80 ID efficiency
  and reco efficiency for this first available pair.
- OFFLINE / READY FOR WP80 BDT-vs-MLP QA: MLP WP80 production campaign
  `mlp_deep_primary_ratios_wp080_20260512_152538` finished cleanly on
  `sphnxuser06` after strict recovery of the three BDT-comparator Jet12 groups.
  Clean watcher
  `mlp_wp80_release_merge_watch_clean_20260512_195445` ended with `0` Justin
  held jobs, `320` chunkMerge ROOT files, `16` ALL ROOT files, and `8` final
  MERGED ROOT files. Codex pulled all 8 merged files locally to
  `/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/combinedSimOnlyEMBEDDED/mlp_wp80_release_20260512_152538`
  and verified each ROOT file opens with one top-level key. This folder now has
  paired `isSimEmbedded` photon12+20 and `isSimEmbeddedInclusive` jet12+20
  outputs for `reference`, `newPPG12`, `auauCentInputBase3x3BDT`, and
  `auauCentInputBase3x3MLP`. Next plotting priority: BDT-vs-MLP WP80
  signal/background comparison, then the matched tight-photon diagnostic
  pipeline Justin asked to track.
- RUNNING IN TMUX / DIAGNOSTIC ONLY: iso-aware kitchen-sink AuAu tight-MLP
  side test on `sphnxuser02`. Started 2026-05-12 19:03 EDT in tmux session
  `mlp_iso_kitchensink_tmux_20260512_190357`, model dir
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/tight_mlp_iso_kitchensink_tmux_20260512_190357`,
  and live log
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/train_iso_kitchensink_tmux_20260512_190357.log`.
  This is the deliberately non-ABCD-safe "how good can we do if isolation is
  visible?" model: kitchen-sink shower/ratio features plus
  `reco_eiso_clip30`, `reco_eiso_over_cluster_Et`, and
  `reco_eiso_signed_log1p`. Evidence at 20:04 EDT: still training, near epoch
  `20`, with validation AUC about `0.846`; no registry yet. Next checkpoint:
  watch for artifact writeout and let the validation watchdog submit the
  diagnostic validation once `applyCheck` passes.
- RUNNING IN TMUX / AUTO-VALIDATION WATCHDOG: SDCC validation handoff for
  standalone kitchen-sink MLPs on `sphnxuser02`. Session
  `mlp_validation_watchdog_20260512_190453`, log
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/mlp_validation_watchdog_20260512_190453.log`,
  script `/tmp/mlp_validation_watchdog_20260512_190453.sh`. It checks the
  non-iso and iso standalone model dirs every 15 minutes; once a registry
  exists and `applyCheck` passes, it submits `validateOnSimCondor` once with
  `groupSize=100`, `SCORE_MAX_ROWS=600000`, and `4000MB` validation workers.
  Codex heartbeat automation `auau-mlp-watchdog` checks the active MLP runs
  every 30 minutes and is authorized to restart the validation handoff if the
  SDCC watchdog dies before validation submission. Evidence at 20:04 EDT: the
  session is alive and correctly waiting because neither standalone registry is
  ready yet.
- RUNNING IN TMUX: standalone kitchen-sink AuAu tight-MLP side test on
  `sphnxuser02`. The previous 48 GB Condor worker `2066012` was removed at
  Justin's request; a follow-up query showed `0` jobs in that cluster. Active
  tmux session is `mlp_kitchensink_tmux_20260512_185247`, model dir
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/tight_mlp_kitchensink_tmux_20260512_185247`,
  and live log
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/train_kitchensink_tmux_20260512_185247.log`.
  This side test uses extended shower-shape / energy-ratio features, high-pT
  WP80 selection, pT-balanced caps/weights, and BDT-guided hard-example
  weighting during training only. Evidence at launch: tmux process tree reached
  `/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python ... train_auau_photon_mlp.py`,
  and the log saw `paths=8000` with balanced signal/background manifest groups.
  Evidence at 20:04 EDT: still training near epoch `108`, with validation AUC
  about `0.746`; no registry yet. Next checkpoint: tail the tmux log for
  artifact writeout; after training finishes, run full embedded validation
  because the existing primary-ratios validation cache lacks the new extended
  input features.
- SUPERSEDED BY ETFINE-ONLY REPAIR: earlier target80 guarded-merge note.
  The raw-output guard worked, but the remaining issue was cfg-tag
  normalization for EtFine BDT modes. Use active tmux
  `target80_merge_repair_etfine_20260513_1233` as the live source of truth.
- FAILED / FIX BEFORE RESUBMIT: high-pT AuAu tight-MLP Condor sweep on
  `sphnxuser02`. Submitted
  2026-05-12 18:23 EDT with DAGMan cluster `2066007`, submit dir
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/condor_sub/auauTightMLPHighPtSweep_20260512_182344`,
  and model root
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/tight_mlp_highpt_sweep_20260512_182344`.
  It trains `global15_highPtBalanced`, `global5_broadReach`, routed
  `pt3_15to35`, and routed `cent3_15to35` with `MAXJOBS=4` and `24000MB`
  requested per worker, then should finalize a sweep manifest and rescored
  validation report from the existing 80-shard MLP validation cache. Evidence
  at 20:04 EDT: no active `patsfan753` Condor jobs; DAG `2066007` aborted with
  `DAG_STATUS_NODE_FAILED`; failed nodes were `GLOBAL5` and `PT3_015_020`; no
  `sweep_manifest.json` and no `validation_rank_table.csv`. `GLOBAL5` failed
  from a finite-positive bin-weight validation error for the broad `5-35`
  model. `PT3_015_020` trained a promising artifact anyway
  (`val_auc=0.78278`, `test_auc=0.78421`, selection fake rate `0.3589`) but
  the wrapper exited nonzero with `unexpected EOF while looking for matching
  '"'`. This does not look memory-related; worker memory use was only a few GB
  against `24000MB`. Next: fix the broad-model pT-bin weighting config/parser
  and the generated worker quoting bug, then resubmit only a corrected sweep
  or salvage/rescore completed artifacts intentionally.
- ACTIVE SLIDE TARGET UNTIL FRIDAY PPG MEETING: working-point Slides deck is
  `AuAuWP_5_12_2026`, presentation id
  `1n4N74TMiB-HiWvZvW0LYW0Emg9AL25yuUSJAAfJzaWg`, URL
  `https://docs.google.com/presentation/d/1n4N74TMiB-HiWvZvW0LYW0Emg9AL25yuUSJAAfJzaWg/edit`.
  Justin explicitly set this as the working point slides deck "from now on"
  until the PPG meeting Friday. Before future slide edits, re-confirm this deck
  through the Google Drive/Slides connector and use live slide object IDs.
- FIXED / RERUN PREP NEEDED: all-available target80 BDT config staging.
  The expanded 5-40 validation rerun `model_validation_condor_20260511_221510`
  is READY and wrote target80 WP files. The first `sphnxuser03` prep attempt
  with tag `bdt_target80_all_available_20260511_224158` stopped only because
  its expanded product map requested optional product `centAsFeat3x3_pt5to40`,
  which is not in that expanded validation report. Codex patched/uploaded
  `scripts/prepare_auau_bdt_target80_available_campaigns.sh` so each product
  map is filtered against the actual WP JSON products; missing optional entries
  now print `[WARN]` and are skipped. The rerun on `sphnxuser03` completed and
  staged configs under
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_generated_configs/bdt_target80_all_available_20260511_224158`.
  No MC jobs were submitted because `RJ_TARGET80_DO_SUBMIT=0`. Next: inspect or
  pull target80 diagnostics, run the config-directory local smoke/preflight,
  then submit paired MC only if that check is clean.
- READY TOOLING: target80 config-directory submit helper. Codex added and
  uploaded `scripts/submit_auau_bdt_target80_config_dir.sh` on 2026-05-11. It
  is the one-command driver for the final target80 MC matrix once the missing
  expanded `5-40` target80 YAML exists: pass
  `RJ_TARGET80_CONFIG_DIR=<dir>`, leave `RJ_TARGET80_DO_SUBMIT=0` for a dry
  smoke, then rerun with `RJ_TARGET80_DO_SUBMIT=1` to submit every
  `analysis_config_*_target80.yaml` sequentially through
  `submit_auau_bdt_targetwp_pair.sh`.
- READY TO RESUBMIT SMARTER: all-available target80 BDT paired MC submission on
  `sphnxuser03`. Justin launched the first all-at-once attempt with
  `RJ_TARGET80_RUN_LOCAL_SMOKE=0 RJ_TARGET80_DO_SUBMIT=1 bash
  ./scripts/submit_auau_bdt_target80_config_dir.sh` for config dir
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_generated_configs/bdt_target80_all_available_20260511_224158`.
  First config `analysis_config_etfine_15to35_target80.yaml` is queued for both
  embedded signal and inclusive-background samples with `34296` raw worker jobs
  and `0` held in the visible queue. The second config
  `analysis_config_expanded_5to40_target80.yaml` has started; observed expanded
  clusters so far are `3032999`-`3033004`, and the latest visible row was
  submitting `AuAuCentInputBDT` for embedded Photon12. The job count was too
  high for one queue, so Justin interrupted it and chose to remove the target80
  jobs. Codex patched/uploaded the submit drivers so the redo submits one
  dataset/config chunk at a time and waits for the matching campaign queue to
  drain before continuing. Next: run the gated driver with a fresh campaign tag
  and leave it alive, preferably in `tmux` or `nohup`, so it can advance
  automatically.
- RUNNING IN TMUX: clean queue-gated target80 BDT campaign. New campaign tag
  `bdt_target80_gated_20260512_001012` is running on `sphnxuser05` in tmux
  session `target80_bdt_target80_gated_20260512_001012`. SSH-auth check showed
  `4287` active matching jobs, clusters `5340030`, `5340031`, `5340032`,
  `0` held. This is the intended bounded behavior. Next status check should
  inspect tmux/logs and Gmail before any new submission. Follow-up terminal
  evidence showed the gate paused correctly after the `isSimEmbedded` half of
  the first YAML: `queued=17148 held=0` before `isSimEmbeddedInclusive`. To
  watch without entering tmux, use
  `tail -f /sphenix/u/patsfan753/scratch/thesisAnalysis/condor_generated_configs/bdt_target80_gated_20260512_001012/tmux_bdt_target80_gated_20260512_001012.log`;
  `Ctrl-C` from that tail is safe and does not stop the tmux submitter.
- MLP FULL-FILE VALIDATION READY / SEPARATE LANE: Gmail reported
  `[RecoilJets][auauTightMLP_validateOnSimCondor][READY]` from
  `sphnxuser06` at 2026-05-11 21:59 for the environment-fixed full-file MLP
  validation (`primary_full_envfix_20260511_215801`); Codex consumed and marked
  the email read. This does not unblock the expanded BDT target80 YAML; it is
  a separate neural-network diagnostic lane. Next: pull and inspect MLP
  validation outputs only when we intentionally return to that lane.
- WAITING ON READY EMAIL / DO NOT SUBMIT MC YET: expanded 5-40 AuAu tight-BDT
  validation rerun for target-80 working points. On 2026-05-11 at 21:53 EDT,
  Justin submitted from `sphnxuser06` after target-80 prep failed on old
  expanded report `model_validation_condor_20260509_192942` because its score
  caches lacked the updated diagnostic column `cluster_weta33_cogx`. Active
  rerun command used
  `RJ_NOTIFY_EMAILS=just0131@gmail.com`,
  `RJ_ML_PYTHON=/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python`,
  `./scripts/auau_tight_bdt_pipeline.sh validateOnSimCondor`,
  source
  `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049`,
  model dir
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_expanded_20260509_152604`,
  registry
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_expanded_20260509_152604/model_registry.json`,
  and `groupSize 100`. Terminal evidence: report
  `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/model_validation_condor_20260511_215329`,
  submit dir
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/condor_sub/auauTightBDTValidate_20260511_215329`,
  DAG
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/condor_sub/auauTightBDTValidate_20260511_215329/auau_tight_bdt_validateOnSimCondor.dag`,
  score-cache dir
  `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/model_validation_condor_20260511_215329/score_caches`,
  `8000` ROOT files, `80` shards, request memory `2500MB`, and DAGMan cluster
  `1274333`. Next: when this report is READY, derive target-80 working points
  from `20260511_215329`, regenerate the missing expanded frozen YAML
  `analysis_config_expanded_5to40_target80.yaml`, inspect the target-80 fit
  plots in chat, then submit the paired MC matrix only after the expanded YAML
  joins the five already staged target80 configs. Update at 22:03 EDT:
  `condor_q` on `sphnxuser06` drained to zero, but
  `cat "$expval/validation_summary.txt"` returned `No such file or directory`;
  this means the expanded rerun is not confirmed READY. Diagnose the DAG/report
  before trying to derive target80 cuts or submit MC from the expanded family.
- VALIDATION READY / HOLD BEFORE MC: fine-`E_T` AuAu BDT target-80 staging.
  On 2026-05-11, `sphnxuser06` training finished READY for
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_etfine_centstudy_current`
  with products `centInput_pt1535`, `ptFine_centInput`, `ptFine_cent3`, and
  `ptFine_cent7`. Evidence: training-tree validation passed with
  `entries=12000316`, `signal=8526032`, `background=3474284`, `files=8000`;
  cache `training_matrix_etfine_centstudy.npz` and registry
  `model_registry.json` were written; `expanded campaign trained reports=89
  selected_specs=89`; `expanded applyCheck opened 89 TMVA ROOT files`; final
  summary printed `RECOILJETS_AUAU_TIGHT_BDT_ETFINE_CENTSTUDY_TRAINING_V1
  status=READY`. Justin then submitted `validateOnSimCondor` from `sphnxuser06`
  with DAG cluster `1273669` and report
  `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/model_validation_condor_20260511_194832`.
  Gmail at 2026-05-11 19:51 reported validation `status=READY`, finite score
  fraction `1`, `scored_entries=400000`, and AUCs: `centInput_pt1535=0.773122`,
  `ptFine_centInput=0.776523`, `ptFine_cent3=0.78874`, `ptFine_cent7=0.808249`;
  Codex consumed and marked the READY email read. Justin then ran
  `deriveWorkingPointsFromValidation` with `TARGET=0.80`, which wrote
  `bdt_working_points_target80.{json,yaml,csv}`,
  `bdt_working_points_target80_runtime_fragment.yaml`, and
  `bdt_working_point_target80_diagnostics.png`. Codex pulled the full report
  locally to
  `/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auauTightBDTValidation/model_validation_condor_20260511_194832`.
  Hard stop now: inspect the threshold/fit plots in chat before any
  `isSimEmbedded`/`isSimEmbeddedInclusive` full submission.
- WATCH / MERGE LATER: fixed width-window WP0.50 paired RecoilJets MC campaign.
  Justin clean-resubmitted from `sphnxuser04` on 2026-05-11 after removing old
  jobs on that submit host. Campaign tag:
  `widthstudy_windows_wp050_fixed_20260511_180500`; submit log:
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/submit_widthstudy_windows_wp050_fixed_20260511_180500.log`.
  This is the fixed analysis-only run: `RJ_WIDTHSTUDY_AUTO_MERGE=0`,
  `RJ_HOLD_FAILED_WORKERS=1`, `RJ_REQUEST_MEMORY=10000MB`,
  `RJ_AUTO_MEMORY_RETRY_CAP_MB=16000`, `RJ_SIM_FIRSTROUND_REQUEST_MEMORY=8000MB`.
  It submitted all three pT windows (`pt5to35`, `pt10to35`, `pt15to35`) for
  both signal `isSimEmbedded` and inclusive-background `isSimEmbeddedInclusive`.
  Worker clusters span `1140092`-`1140103`, with `1429` jobs each; submit order
  per window is Photon12, Photon20, Jet12, Jet20. Bulk roots are under
  `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded_widthstudy_windows_wp050_fixed_20260511_180500_pt*to35`
  and
  `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembeddedinclusive_widthstudy_windows_wp050_fixed_20260511_180500_pt*to35`.
  Merge/output roots are
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_widthstudy_windows_wp050_fixed_20260511_180500_pt5to35`,
  `_pt10to35`, and `_pt15to35`. Final submission snapshot showed `11674`
  idle, `5474` running, and `0` held for the user query; `7210 removed` were
  the deliberate preflight cleanup of stale jobs. Gmail check right after
  submission found two fresh unread
  `[RecoilJets][auto_simembeddedinclusive_final_ready][FAILED]` messages from
  the removed pre-fix auto workflow; Codex marked them read after confirming
  they do not belong to the new analysis-only cluster run. Next: watch
  Gmail/queue until clusters `1140092`-`1140103` drain; if clean, run the
  wrapper-printed manual `mergeRecoilJets.sh firstRound/secondRound/finalStitch`
  commands window by window. If any held jobs appear, inspect the hold reason
  before removing or releasing.
- VALIDATION READY / READY FOR MC SUBMISSION: 9-model width-window AuAu
  tight-BDT campaign.
  Justin trained on `sphnxuser04` from the fixed extraction
  source `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049`
  into model dir
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_centinput_widthstudy_windows_current`.
  Initial command observed in terminal:
  `RJ_NOTIFY_EMAILS=just0131@gmail.com RJ_ML_PYTHON=/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python RJ_AUAU_BDT_TRAIN_PARALLEL=3 ./scripts/auau_tight_bdt_pipeline.sh trainWidthStudyWindowsFromExtraction SOURCE=/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049 MODEL_DIR=/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_centinput_widthstudy_windows_current PT_WINDOWS=5:35,10:35,15:35`.
  That first attempt reached `8000/8000` and passed with
  `entries=12000316`, `signal=8526032`, `background=3474284`, `files=8000`; the
  shared cache `training_matrix_widthstudy_windows.npz` was written, then the
  Python training process was killed by the system on line 887. Treat as a
  recoverable memory-pressure failure. Current follow-up evidence on
  2026-05-11: Justin reran serially with
  `RJ_AUAU_BDT_TRAIN_PARALLEL=1 RJ_AUAU_BDT_XGB_N_JOBS=1` in the same
  `MODEL_DIR`; a `ps` check from another `sphnxuser04` terminal showed wrapper
  PID `3591301` and active Python trainer PID `3593547` using about `5.4 GB`
  RSS and `99.8%` CPU with `--parallel-workers 1 --n-jobs 1`. Expected model
  matrix is three centrality-input
  width variants (`base widths`, `3x3 widths`, `base+3x3 widths`) for three
  training/applicability windows (`5-35`, `10-35`, `15-35` GeV). Final
  terminal evidence on 2026-05-11: `model_registry.json` was written,
  `expanded campaign trained reports=9 selected_specs=9`, `expanded applyCheck
  opened 9 TMVA ROOT files`, and
  `RECOILJETS_AUAU_TIGHT_BDT_WIDTHSTUDY_WINDOWS_TRAINING_V1 status=READY`.
  Validation submission evidence on 2026-05-11 from `sphnxuser04`:
  Justin submitted `validateOnSimCondor` for this registry, report
  `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/model_validation_condor_20260511_171943`,
  submit dir
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/condor_sub/auauTightBDTValidate_20260511_171943`,
  and DAGMan cluster `1139989`. Queue evidence later drained to zero. Gmail
  on 2026-05-11 reported
  `[RecoilJets][auauTightBDT_validateOnSimCondor][READY]` from `sphnxuser04`
  for this report; Codex consumed the email and marked it read. Next steps:
  pull and inspect the validation report, then submit paired `isSimEmbedded` and
  `isSimEmbeddedInclusive` campaigns using
  `scripts/submit_auau_bdt_widthstudy_windows_wp050.sh` so production uses BDT
  score `> 0.50`.
- JSTG POSTPONED ONE WEEK / BDT DEVELOPMENT SPRINT:
  On 2026-05-11 Justin reported JSTG was postponed by one week. Treat this as
  valuable runway to push the AuAu tight-BDT program from promising validation
  into a working production-style candidate by next week, not as permission to
  rush uninspected MC. Priority lanes for the extra week: finish validation-only
  passes for fine-`E_T`, width-window, and available expanded registries; derive
  target-80 working points and inspect threshold/fit diagnostics in chat before
  any full MC submissions; run approved `isSimEmbedded` and
  `isSimEmbeddedInclusive` campaigns for the best variants; pull outputs and
  build the comparison suite covering ID/reco efficiency, background tight
  rate/fake leakage, tight/complement shower shapes, isolation/ABCD behavior,
  centrality/`E_T` stability, and `x_{J#gamma}` stability. Use the slide time
  to distill a coherent next-week JSTG story, with deeper backup detail saved
  for the following PPG19/PPG-style discussion.
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
  4. DONE locally 2026-05-11: pulled/pair-checked the four 3x3 one-row
     outputs from `output_bdt3x3_20260510_223529_wp050` and
     `output_bdt3x3_20260510_223529_wp080` into
     `dataOutput/auau_bdt3x3_mc_validation/{wp050,wp080}`. Added
     `macros/CompareAuAu3x3IDEfficiency.C` and ran a ROOT batch-mode
     ID-efficiency check against the standard full-cluster centrality-input
     WP0.50 BDT. Outputs:
     `dataOutput/auau_bdt3x3_mc_validation/id_efficiency_compare/bdt3x3_vs_fullcluster_id_efficiency.csv`,
     `bdt3x3_vs_fullcluster_signal_tight_fraction.png`, and
     `bdt3x3_vs_fullcluster_background_tight_fraction.png`.
  5. NEXT: compare the 3x3-width BDT against the standard centrality-input and
     pT/centrality-binned BDTs in reco efficiency, isolation, ABCD behavior,
     shower-shape templates, leakage/fake rates, and xJ-sensitive inputs.
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

- LATER PRIORITY / PPG12 CLUSTERING SETTINGS: understand and play with the
  PPG12 splitting and topo-clustering configurations before making any
  long-term photon-ID clustering choices. Reference `ppg12codeGit/` and
  distinguish the two axes cleanly: split vs no-split EMCal photon-cluster
  containers (`CLUSTERINFO_CEMC` vs `CLUSTERINFO_CEMC_NO_SPLIT`) and
  topo-cluster isolation choices (`use_topo_iso`). Goal for later: decide
  whether Au+Au should add a controlled split/no-split BDT comparison and/or a
  separate topo-isolation ABCD/systematics comparison.
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
- Build a future BDT-vs-MLP matched tight-photon validation pipeline. After
  the full-stat MLP and selected BDT working points are available in matching
  `isSimEmbedded` / `isSimEmbeddedInclusive` outputs, compare candidate-level
  overlap and disagreement: tight in both, BDT-only tight, MLP-only tight, and
  rejected by both. Summarize truth matching, photon `E_T`, centrality,
  isolation, shower-shape distributions, WP80 fake rate, and RecoilJets-level
  recovery/leakage for each category so the MLP/BDT score difference becomes a
  physics diagnostic.
- Run the BDT-beating AuAu tight-MLP v2 decision lane. Local code now has a
  guarded SDCC driver at `scripts/auau_mlp_bdt_beating_driver.sh`; mutating
  modes require `RJ_DO_RUN=1`, Condor validation additionally requires
  `RJ_ALLOW_CONDOR=1`. Next operational sequence: recover `global5_broadReach`
  if still missing, write/rescore the high-pT manifest, train
  `auauHighPtDistilledKitchenMLP_v2`, validate it against current MLP /
  high-pT sweep / best BDT references, then only promote to production if the
  `20-35 GeV` WP80 behavior improves cleanly.
- DONE / USE AS CEILING DIAGNOSTIC: stacked BDT+MLP calibrator sidecar.
  Codex added `scripts/train_auau_stacked_bdt_mlp_calibrator.py` and
  `scripts/submit_auau_stacked_bdt_mlp_calibrator.sh`, verified local
  synthetic-cache tests and SDCC two-shard smoke, then ran the full stacker on
  `sphnxuser02` in tmux session
  `mlp_stacked_bdt_mlp_calib_20260512_225152`. Output dir:
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/stacked_bdt_mlp_calibrator_20260512_225152`.
  Best test model was `stack_interactions_tiny_gbm` with AUC `0.97004` and
  WP80 fake `0.0405`; strongest pT-bin test AUCs were around `0.75-0.81` in
  `20-35 GeV`. This is not an ABCD-safe production ID because the BDT score is
  an input, but it strongly supports using the BDT as a teacher/ensemble
  diagnostic for the next clean MLP distillation step.
- TRAINED / VALIDATION FALLBACK NEEDED: fine-pT distilled kitchen MLP recovery
  sweep. The first `sphnxuser02` DAGMan cluster `2066023` failed fast because
  route-local experts inherited global high-pT selection weights
  (`15:20,20:25,25:35`) that did not match single-route bins such as `15:18`.
  The first recovery DAG `2066036` still failed because the remote v2 pipeline
  defaults overrode the route-local CLI argument. Codex hardened
  `scripts/submit_auau_mlp_finept_distilled_sweep.sh` so generated workers
  export route-local pT/ranking env vars before invoking the pipeline, and
  updated `scripts/auau_tight_mlp_pipeline.sh` to print effective high-pT
  selection weights in the training banner. Local and remote checks passed:
  `bash -n`, Python compile, route-weight parser smoke, local DAG dry-run,
  generated-worker syntax, and remote `/tmp` DAG dry-run from the uploaded SDCC
  copy. SDCC upload/status returned `MATCH` for both files. Fresh recovery was
  submitted from `sphnxuser02` with stamp `20260512_215558`: DAGMan cluster
  `2066049`, sweep dir
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/tight_mlp_finept_distilled_kitchen_v2_recover2_20260512_215558`,
  submit root
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/condor_sub/auauTightMLPFinePtDistilled_recover2_20260512_215558`.
  Early queue check showed first route workers `2066050` and `2066051`
  running with `0` held; their `.out` files printed
  `route_highpt_selection_weights=15:18:1.0` / `18:20:1.0`, the pipeline
  banner printed matching `high-pT selection weights`, `.err` files were empty,
  and both workers were reading ROOT files. This resolves the repeated
  `Weight bin 15:20 is not present` failure for the first routes. 22:55 EDT
  update: all six route experts trained, all six `model_registry.json` files
  exist, and `finept_distilled_kitchen_v2_sweep_manifest.json` exists. The
  automatic validation-cache rescore failed because the old primary-MLP cache
  does not contain kitchen-sink feature columns needed by these artifacts
  (`cluster_weta35_cogx`, `cluster_wphi53_cogx`, `cluster_w32`, and related
  ratios). Next: run a fresh ROOT-backed routed validation using the sweep
  manifest, or generate a kitchen-feature score cache; do not use training AUC
  as the promotion evidence.

## Waiting

- Width-window WP0.50 RecoilJets MC validation:
  - Current failed auto-DAG evidence says sampled workers exited cleanly and
    were removed only because DAGMan hit its max-held-job abort limit.
  - Local fix staged 2026-05-11: `RecoilJets_Condor_submit.sh` adds
    `on_exit_hold` for nonzero/signal worker exits; width-window WP0.50 submit
    wrapper defaults to `RJ_WIDTHSTUDY_AUTO_MERGE=0` and prints explicit
    post-worker merge commands.
  - Uploaded fix to SDCC on 2026-05-11:
    `./scripts/sftp_push_recoiljets.sh RecoilJets_Condor_submit.sh
    submit_auau_bdt_widthstudy_windows_wp050.sh`. No rebuild needed.
  - Next: cleanly remove stale/superseded failed campaign jobs with Justin's
    approval, then resubmit the width-window WP0.50 campaign from the chosen
    submit host.
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
- Three-slice `isSimEmbeddedInclusive` inclusive3 stitching QA from
  `sphnxuser02`: DAG `2116548` completed, the final merged ROOT was opened, and
  the slide-10-equivalent PNG was regenerated/pulled at
  `dataOutput/embeddedXsec/inclusive3_20260515_172346/embeddedInclusiveJet_finalMergedStitchedTruthJetPtSpectrum.png`.
  The visible sample colors are now Jet12 blue, Jet20-to-30 orange, and Jet30
  pink; the plot was visually inspected in chat on 2026-05-16. Next: use this
  PNG for the three-sample stitching slide or ask before archiving the task.
- Fresh SDCC smoke or production evidence for direct-fanout pp and SIM jobs.
- User decision on when to transfer the ABCD/purity tight-axis fix and rerun
  affected outputs.
- User decision on final `isSimInclusive` pp inclusive-jet stitching
  prescription.
- BDT-beating AuAu tight-MLP automation lane: the new driver and v2
  runtime/trainer support were uploaded to SDCC on 2026-05-12 and verified
  `MATCH`. Use read-only `status` first, then gated `recoverGlobal5Tmux`,
  `writeHighPtManifest`, `rescoreHighPt`, and `trainV2Tmux` as evidence
  becomes available. The `auau-mlp-watchdog` heartbeat is active and read-only;
  it must not submit jobs, remove jobs, transfer files, or edit SDCC without
  Justin explicitly asking in the live thread.

## Done Pending Removal

- Fixed width-window WP0.50 RecoilJets MC campaign
  `widthstudy_windows_wp050_fixed_20260511_180500` completed finalStitch and
  was pulled offline on 2026-05-12. Local folder:
  `/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auau_widthstudy_windows_wp050_fixed_20260511_180500_complete`.
  Verification: `18` nonzero ROOT files, `3.9G` total. Contents cover three
  pT training/applicability windows (`pt5to35`, `pt10to35`, `pt15to35`), both
  embedded MC datasets (`photonJet12and20merged_SIM`,
  `embeddedJet12and20merged_SIM`), and three rows: reference box-cuts,
  centrality-input BDT with base widths, and centrality-input BDT with 3x3
  widths. The base+3x3 row was not produced by this remote WP0.5 campaign.
  Next analysis task is offline plot comparison; ask Justin before archiving
  this completed operational task.
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
