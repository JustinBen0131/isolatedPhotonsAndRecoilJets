# Status Dashboard

Last updated: 2026-05-12

Use traffic-light status for datasets/jobs/outputs. Keep rows evidence-based:
email, terminal output, job ID, pulled file timestamp, ROOT inspection, or clear
user statement.

Latest active update:
- 2026-05-12 16:42 EDT MLP WP80 production held-job diagnosis: live
  `sphnxuser06` Condor query found exactly `3` held jobs, all in cluster
  `1274730` and all belonging to the BDT comparator row, not the MLP row.
  Held jobs: `1274730.218` on `sphnx0252.sdcc.bnl.gov` group/list `219`,
  `1274730.219` on `sphnx1117.sdcc.bnl.gov` group/list `220`, and
  `1274730.221` on `sphnx1285.sdcc.bnl.gov` group/list `222`. Dataset/sample:
  `isSimEmbeddedInclusive` / `run28_embeddedJet12`; cfg row:
  `preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant`.
  Hold reason for all three:
  `RecoilJets worker exited nonzero or by signal; inspect stdout/err before release/remove`.
  Profiles show `exit_code=1`, `elapsed_seconds=25-34`, `max_rss_kb` about
  `775852-786680`, and `request_memory_mb=12000`, so this is not a memory
  failure. Stderr only showed ROOT class-table warnings; stdout ended with
  `[ERROR] Fun4All macro failed (rc=1)`. A tab-aware stat check of the held
  group input lists found all `84` ROOT input files exist and are nonzero
  (`28` files per held group); a deeper PyROOT open probe was stopped because
  it was too slow/noisy. Separate MLP-only cluster check at 16:39 showed
  clusters `1274724`, `1274725`, `1274732`, `1274733` had `0` held jobs,
  `1561` running and `2858` idle; MLP signal outputs were already complete
  with `2858` ROOT files (`1429` photon12 + `1429` photon20), while MLP
  inclusive jobs were still queued. Next: do not touch the clean MLP lane; if
  the held count stays isolated, debug/resubmit only these BDT-comparator
  `embeddedJet12` groups after the queue drains or decide whether the
  comparator row is dispensable for the MLP production comparison.
- 2026-05-12 16:33 EDT Codex submitted the authorized guarded `finalStitch` stage on `sphnxuser04` for flat BDT >0.5 width-window campaign `widthstudy_windows_wp050_fixed_20260511_180500`. Preflight required `0` active/held `patsfan753` jobs and verified each window/dataset (`pt5to35`, `pt10to35`, `pt15to35` x `simembedded`, `simembeddedinclusive`) had `6` secondRound flat files. The staged helper submitted 18 finalStitch DAGs covering 3 windows x 2 datasets x 3 cfg rows; observed DAG/worker clusters span `1140257` through `1140292`. Stage log: `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_generated_configs/widthstudy_windows_wp050_fixed_20260511_180500/merge_finalStitch_20260512_163253.log`. At 16:39 EDT queue showed `144` idle, `0` running, `0` held, and final merged count still `0/18`; this is waiting for slots, not failed. The active heartbeat now watches this lane plus the target80 lane every 15 minutes and is authorized to pull offline automatically once final merged count reaches `18/18` with queue `0` and no held jobs.
- 2026-05-12 15:38 EDT MLP deep primary-ratios WP80 paired RecoilJets MC
  production is submitted from `sphnxuser06`. Driver header:
  `RECOILJETS_AUAU_MLP_DEEP_PRIMARY_RATIOS_WP80_PAIR_DRIVER_V1`; campaign tag
  `mlp_deep_primary_ratios_wp080_20260512_152538`; generated YAML:
  `condor_generated_configs/mlp_deep_primary_ratios_nostat_fullval_20260512_145449_wp080/analysis_config_mlp_deep_primary_ratios_nostat_fullval_20260512_145449_wp080.yaml`;
  submit log:
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_generated_configs/mlp_deep_primary_ratios_wp080_20260512_152538/submit_mlp_deep_primary_ratios_wp080_20260512_152538.log`.
  Knobs: `groupSize=7`, worker memory `12000MB`, retry cap `16000MB`,
  merge memory hint `8000MB`, merge group `75`, `auto_merge=0`, notify
  `just0131@gmail.com`, analysis-only workers. Signal bulk root:
  `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded_mlp_deep_primary_ratios_wp080_20260512_152538`;
  inclusive bulk root:
  `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembeddedinclusive_mlp_deep_primary_ratios_wp080_20260512_152538`;
  shared merge/output root:
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_mlp_deep_primary_ratios_wp080_20260512_152538`.
  Signal snapshot/fanout: `simembedded_20260512_152539`,
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_yaml_overrides/simembedded_condorDoAll_20260512_152539`;
  inclusive snapshot/fanout: `simembeddedinclusive_20260512_153305`,
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_yaml_overrides/simembeddedinclusive_condorDoAll_20260512_153305`.
  Worker clusters are `1274718`-`1274725` for `isSimEmbedded` and
  `1274726`-`1274733` for `isSimEmbeddedInclusive`, covering the four YAML
  rows (reference, newPPG12, best current BDT reference, and MLP) times the
  two samples in each dataset. Terminal post-submit queue summary:
  `22864` jobs total, `7213` running, `15651` idle, `0` held. Printed tracker:
  `TRACK_THIS_CAMPAIGN=mlp_deep_primary_ratios_wp080_20260512_152538`.
  Next checkpoint: watch Gmail `RecoilJets Pipeline` and/or compact SDCC queue
  status for held jobs or READY/CHECK/FAILED emails; after both halves finish,
  run the printed/manual merge step because `auto_merge=0`, then pull the
  isolated MLP output root for BDT-vs-MLP comparisons.
- 2026-05-12 EDT deep primary-ratios MLP WP80 RecoilJets runtime check passed
  on `sphnxuser06`. Justin first rebuilt `src_AuAu` successfully with
  `makeProject`, installing to
  `/sphenix/u/patsfan753/thesisAnalysis_auau/install`. The first local check
  failed only because Codex gave the wrong submitter path
  `./scripts/RecoilJets_Condor_submit.sh`; corrected command used repo-base
  `./RecoilJets_Condor_submit.sh`. Runtime check command used YAML
  `condor_generated_configs/mlp_deep_primary_ratios_nostat_fullval_20260512_145449_wp080/analysis_config_mlp_deep_primary_ratios_nostat_fullval_20260512_145449_wp080.yaml`,
  `RJ_PHOTON_ID_ROW_MATCH=auauCentInputBase3x3MLP`,
  `RJ_ID_FANOUT_MAX_ROWS=1`, dataset `isSimEmbedded`, sample
  `run28_embeddedPhoton20`, local `200` events and `NFILES=1`. Terminal
  evidence: row filter matched `1/4` photon-ID rows; YAML parsed OK; C++
  printed `[AuAuTightMLP] loaded auauCentInputBase3x3MLP model ... with 16
  features and 4 layers`; job profile ended `exit_code=0 elapsed_seconds=54`;
  output ROOT file written under
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/local_sim_outputs/simembedded/preselectionNewPPG12_tightAuAuCentInputBase3x3MLP_nonTightAuAuMLPComplement_baseVariant/run28_embeddedPhoton20/`
  with size `7917112` bytes. Next: submit paired `isSimEmbedded` and
  `isSimEmbeddedInclusive` MLP WP80 production with isolated output roots,
  `auto_merge=0`, and track clusters/output roots immediately.
- 2026-05-12 EDT WP80 YAML generation succeeded for the deep primary-ratios
  MLP. Generated YAML:
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_generated_configs/mlp_deep_primary_ratios_nostat_fullval_20260512_145449_wp080/analysis_config_mlp_deep_primary_ratios_nostat_fullval_20260512_145449_wp080.yaml`.
  Terminal grep confirmed photon-ID row `[newPPG12,
  auauCentInputBase3x3MLP, auauMLPComplement]`, model dir
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/tight_mlp_deep_primary_ratios_nostat_20260512_005242`,
  and runtime WP entry
  `auauCentInputBase3x3MLP|binned|15;20;25;35|0.528583622;0.5030803084;0.5226714969|15|35|1`.
  Next: rebuild/check RecoilJets runtime with the uploaded MLP width-ratio C++
  support, then submit paired `isSimEmbedded`/`isSimEmbeddedInclusive` with
  isolated MLP output roots if the local/runtime check passes.
- 2026-05-12 EDT score-cache label audit for the deep primary-ratios MLP
  validation is clean. Visible terminal evidence from `sphnxuser06`: report
  has `80` score-cache `.npz` files plus `validation_metrics.json` and
  `mlp_working_points_target80.json`; corrected cache audit found
  `cache_label_entries=400000`, `cache_signal_entries=177133`, and
  `cache_background_entries=222867`. First cache keys include `is_signal`,
  diagnostic features, width-ratio inputs, `centrality`, `reco_eiso`, and
  `score_centInputBase3x3WidthRatiosMLP_pt1535`. Metrics include AUC
  `0.8707484856214256`, finite fraction `0.9998975`, ECE
  `0.07544835409327035`, score-vs-centrality `-0.0671`, score-vs-`E_T`
  `0.4838`, and score-vs-Eiso `-0.3870`. Centrality-bin AUCs improve from
  central to peripheral (`0-20`: `0.8234`, `20-50`: `0.8604`, `50-80`:
  `0.8909`); pT-bin AUCs are strong at `15-20` (`0.8281`) but degrade at
  `20-25` (`0.7121`) and `25-35` (`0.6942`). Next: generate WP80 runtime YAML
  and decide whether to run paired MC production as a candidate despite the
  high-pT degradation.
- 2026-05-12 14:57 EDT Gmail `RecoilJets Pipeline` reported
  `[RecoilJets][auauTightMLP_validateOnSimCondor][READY]` for the deep
  primary-ratios MLP validation submitted from `sphnxuser06`; Codex consumed
  and marked the email read. Report dir:
  `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/mlp_model_validation_condor_deep_primary_ratios_nostat_fullval_20260512_145449`;
  model dir:
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/tight_mlp_deep_primary_ratios_nostat_20260512_005242`.
  Email metrics: `scored_entries=400000`, `signal_entries=8526032`,
  `background_entries=3474284`,
  `product.centInputBase3x3WidthRatiosMLP_pt1535.auc=0.8707484856214256`,
  finite fraction `0.9998975`, WP80 threshold `0.5027055621147155`, WP80 fake
  rate `0.18650068883194745`, score-vs-Eiso correlation
  `-0.38702614706090505`. The summary still reports
  `scored_signal_entries=0` and `scored_background_entries=0`, matching the
  earlier known summary quirk; next checkpoint is a score-cache label audit
  before generating/promoting WP YAML.
- 2026-05-12 14:55 EDT AuAu tight-MLP deep primary-ratios full embedded
  validation was submitted from `sphnxuser06`. Driver header:
  `RECOILJETS_AUAU_MLP_DEEP_PRIMARY_RATIOS_FULL_VALIDATE_V1`; stamp:
  `deep_primary_ratios_nostat_fullval_20260512_145449`; source:
  `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049`;
  model dir:
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/tight_mlp_deep_primary_ratios_nostat_20260512_005242`;
  report:
  `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/mlp_model_validation_condor_deep_primary_ratios_nostat_fullval_20260512_145449`;
  submit dir:
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/condor_sub/auauTightMLPValidate_deep_primary_ratios_nostat_fullval_20260512_145449`;
  DAG:
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/condor_sub/auauTightMLPValidate_deep_primary_ratios_nostat_fullval_20260512_145449/auau_tight_mlp_validateOnSimCondor.dag`;
  DAGMan cluster `1274599`. Preflight queue on `sphnxuser06` showed `0`
  Justin jobs; submitter planned `8000` ROOT files, `80` shards, `groupSize=100`,
  `scoreMaxPerShard=5000`, `RJ_AUAU_TIGHT_MLP_VALIDATE_TOTAL_SCORE_MAX_ROWS=400000`,
  and request memory `3000MB`. Post-submit queue showed the DAGMan parent
  running. Next checkpoint: watch for worker jobs, then a READY/CHECK Gmail
  stage email or inspect report score caches for `80/80` shards.
- 2026-05-12 EDT AuAu tight-MLP deep primary-ratios no-stat training finished
  cleanly on `sphnxuser06`. Session/log/model from
  `mlp_deep_primary_ratios_nostat_20260512_005242`. Terminal evidence:
  `early stopping at epoch 156; best_epoch=101`,
  selected product `centInputBase3x3WidthRatiosMLP_pt1535` with hidden
  `[128,64,32]`, validation AUC `0.70904`, test AUC `0.70986`, WP80 fake rate
  `0.4451753747634988`, artifact
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/tight_mlp_deep_primary_ratios_nostat_20260512_005242/auau_tight_mlp_centInputBase3x3_pt1535.json`,
  registry
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/tight_mlp_deep_primary_ratios_nostat_20260512_005242/model_registry.json`,
  `[OK] applyCheck loaded and scored 1 MLP JSON artifacts`, and
  `DONE_MLP_MODEL_DIR=...tight_mlp_deep_primary_ratios_nostat_20260512_005242`.
  Next: run full embedded validation against the current BDT extraction source
  before any production promotion. The train-side AUC is lower than hoped, so
  treat this as a serious candidate to validate, not as automatically better
  than the BDT.
- 2026-05-12 00:52 EDT AuAu tight-MLP deep primary-ratios no-stat run is live
  in `tmux` on `sphnxuser06`. Session:
  `mlp_deep_primary_ratios_nostat_20260512_005242`; source:
  `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049`;
  model dir:
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/tight_mlp_deep_primary_ratios_nostat_20260512_005242`;
  log:
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/train_deep_primary_ratios_nostat_20260512_005242.log`.
  Terminal evidence: ML env check printed `ML env ok`; `tmux ls` showed one
  window; trainer printed `paths=8000 check_exists=False`, a balanced manifest
  with `2000` files each for embedded jet12, jet20, photon12, and photon20, and
  reached `reading 1/8000`. This confirms the no-stat patch is active and the
  run is no longer stuck in the old pre-read GPFS wait. Next checkpoint: watch
  for read progress beyond early files, then the first epoch line with
  validation AUC.
- 2026-05-12 EDT AuAu tight-MLP deep single-candidate restart attempt
  `mlp_deep_primary_ratios_20260512_003834` on SDCC was intentionally killed
  after terminal diagnostics showed it alive but stuck before first training
  read (`cl_sync_io_wait`, low RSS, log still at the pipeline header). Likely
  culprit: the trainer was stat-checking all 8000 GPFS input files before
  printing any read progress. Codex patched
  `scripts/train_auau_photon_mlp.py` so the all-file existence check is opt-in
  via `--check-input-files` / `RJ_AUAU_MLP_CHECK_INPUT_FILES=1`; default now
  prints input expansion immediately and lets the ROOT reader fail on the first
  truly missing file. Local `py_compile` passed and SFTP upload completed for
  `scripts/train_auau_photon_mlp.py`,
  `scripts/validate_auau_tight_mlp_on_sim.py`, and
  `scripts/auau_tight_mlp_pipeline.sh`. Next: relaunch with a fresh `nostat`
  stamp and watch for `[trainAuAuPhotonMLP] input path expansion ...`, then
  `manifest`, then `reading 1/8000`.
- 2026-05-12 EDT AuAu tight-MLP deep single-candidate restart is staged.
  Codex locally patched and SFTP-uploaded
  `scripts/train_auau_photon_mlp.py`,
  `scripts/validate_auau_tight_mlp_on_sim.py`,
  `scripts/auau_tight_mlp_pipeline.sh`, and
  `src_AuAu/RecoilJets_AuAu.cc`. New mode:
  `trainPrimaryDeepFromExtraction`, defaulting to product `primary-ratios`
  (`centInputBase3x3WidthRatiosMLP_pt1535`) exported as runtime variant
  `auauCentInputBase3x3MLP` and filename
  `auau_tight_mlp_centInputBase3x3_pt1535.json`. Default deep knobs:
  pT `15:35`, `max_rows_per_class=350000`, `epochs=260`, `patience=55`,
  hidden layers `128,64,32`, `restarts=1`, selection metric
  `validation_auc`, batch `4096`, learning rate `7e-4`, L2 `5e-5`,
  conditional jitter `0.015`, input clip `7.0`. Local checks passed:
  `bash -n`, Python `py_compile`, and product-alias/feature-order smoke.
  Next: run the paste-ready tmux block from SDCC, watch log
  `train_deep_primary_ratios_<STAMP>.log`, and only promote if the full
  registry writes cleanly and validation beats the current smoke/full-file
  baseline.
- 2026-05-12 00:11 EDT durable queue-gated target80 campaign is running in
  detached `tmux` on `sphnxuser05`. Session:
  `target80_bdt_target80_gated_20260512_001012`; campaign tag:
  `bdt_target80_gated_20260512_001012`. Codex checked via SSH auth:
  `tmux ls` shows the session alive, `tmux` exists at `/usr/bin/tmux`, active
  matching target80 jobs = `4287`, held target80 jobs = `0`. Current clusters:
  `5340030`, `5340031`, `5340032`, all under the first YAML family
  `analysis_config_etfine_15to35_target80`; visible queue tail showed
  `isSimEmbedded`/`run28_embeddedPhoton12` `newPPG12` row in cluster `5340032`.
  Bulk root in visible jobs:
  `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded_bdt_target80_gated_20260512_001012_analysis_config_etfine_15to35_target80`.
  Latest log directory observed:
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_generated_configs/bdt_target80_gated_20260512_001012_analysis_config_etfine_15to35_target80`;
  helper log:
  `submit_bdt_target80_gated_20260512_001012_analysis_config_etfine_15to35_target80.log`.
  The queue gate is doing what we want: it submitted a small active chunk, not
  the full ~183k-job matrix. Tomorrow/resume: SSH to `sphnxuser05`, inspect
  `tmux ls`, tail the matching `condor_generated_configs/bdt_target80_gated_*`
  logs, and check `condor_q -nobatch patsfan753`; do not launch another
  target80 submit unless this tmux session exits cleanly or is intentionally
  stopped.
- 2026-05-12 00:21 EDT visible `sphnxuser05` terminal confirmed the queue gate
  worked for the first target80 chunk. After submitting the `isSimEmbedded`
  half of `analysis_config_etfine_15to35_target80`, it stopped at:
  `[queue_gate] before isSimEmbeddedInclusive ... queued=17148 held=0
  resume_below=0 max_queued=40000`. This means it will not submit the
  inclusive-background half until the current `isSimEmbedded` jobs drain to
  zero. A later visible `condor_q` showed `17065` active user jobs
  (`11317` idle, `5748` running, `0` held), still below the 40k cap. Safe
  monitoring command that does not affect tmux:
  `tail -f /sphenix/u/patsfan753/scratch/thesisAnalysis/condor_generated_configs/bdt_target80_gated_20260512_001012/tmux_bdt_target80_gated_20260512_001012.log`;
  pressing `Ctrl-C` from this `tail -f` stops only the tail, not the tmux
  submitter. Inner submitter log:
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_generated_configs/bdt_target80_gated_20260512_001012/submit_bdt_target80_gated_20260512_001012.log`.
  The heartbeat watchdog was adjusted to check hourly until this first queue
  drains and the next chunk submits; after that, switch to 5-hour cadence if
  all remains healthy.
- 2026-05-12 00:0x EDT target80 all-available campaign was intentionally
  stopped before it could submit the full matrix at once. Justin hit `Ctrl-C`
  on `sphnxuser03`; last visible fully submitted cluster was `3033032`
  (`isSimEmbeddedInclusive`, `run28_embeddedJet20`,
  `AuAuCentInputMinOptBDT` target80 row). Approximate visible submission before
  interruption was `58` raw-worker clusters, about `82,882` jobs, from
  `etfine_15to35` plus part of `expanded_5to40`. Justin then chose to remove
  the target80 jobs and redo the campaign in bounded queue partitions rather
  than leave ~183k jobs queued. Codex patched and uploaded
  `scripts/submit_auau_bdt_target80_config_dir.sh` and
  `scripts/submit_auau_bdt_targetwp_pair.sh`: both now default to
  `RJ_TARGETWP_QUEUE_GATE=1`, `RJ_TARGETWP_MAX_QUEUED=40000`,
  `RJ_TARGETWP_RESUME_BELOW=0`, and `RJ_TARGETWP_POLL_SECONDS=300`. The pair
  submitter waits between `isSimEmbedded` and `isSimEmbeddedInclusive`, and the
  config-dir driver waits between YAML configs, stopping immediately if held
  jobs appear. Local checks passed: `bash -n` for both scripts and
  `git diff --check`; SFTP upload completed successfully. Next: resubmit with
  a new campaign tag, not the interrupted tag, so partial raw outputs cannot
  be confused with the clean gated campaign. Follow-up SSH-auth queue check on
  `sphnxuser03` confirmed the removed campaign had `0 idle`, `0 running`, and
  `0 held` target80 jobs; Condor still showed removed `X` records temporarily.
  Codex patched both queue gates again so active-queue counts ignore removed
  and completed jobs (`JobStatus` 3/4), preventing the gated driver from
  waiting on already-dead records. Re-upload of both submit helpers completed
  successfully.
- 2026-05-11 23:4x EDT target80 full paired MC submission is still actively
  queuing from `sphnxuser03`. No fresh unread `RecoilJets Pipeline` emails were
  found for this live submit check. Visible terminal evidence shows the first
  target80 config
  `analysis_config_etfine_15to35_target80.yaml` has been submitted for both
  `isSimEmbedded` and `isSimEmbeddedInclusive`; the queue summary printed
  `34296` jobs for that config family, with `34196 idle`, `100 running`,
  `0 held`, and `0 removed`. The driver then advanced to the second config
  `analysis_config_expanded_5to40_target80.yaml`, campaign
  `bdt_target80_all_available_20260511_224158_analysis_config_expanded_5to40_target80`,
  with bulk root
  `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded_bdt_target80_all_available_20260511_224158_analysis_config_expanded_5to40_target80`,
  merge root
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_bdt_target80_all_available_20260511_224158_analysis_config_expanded_5to40_target80`,
  snapshot
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_snapshots/simembedded_20260511_233928`,
  YAML/fanout dir
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_yaml_overrides/simembedded_condorDoAll_20260511_233928`,
  `groupSize=7`, worker memory `12000MB`, merge memory hint `8000MB`, and
  `auto_merge=0`. Observed expanded clusters so far: `3032999` =
  `isSimEmbedded` `run28_embeddedPhoton12` reference row and `3033000` =
  `isSimEmbedded` `run28_embeddedPhoton20` reference row. Follow-up terminal
  evidence added `3033001` = `isSimEmbedded` `run28_embeddedPhoton12`
  `newPPG12` tight/reference non-tight row, `3033002` = `isSimEmbedded`
  `run28_embeddedPhoton20` `newPPG12` tight/reference non-tight row,
  `3033003` = `isSimEmbedded` `run28_embeddedPhoton12`
  `AuAuNoCentBDT` target80 row, and `3033004` = `isSimEmbedded`
  `run28_embeddedPhoton20` `AuAuNoCentBDT` target80 row. Latest visible line
  was submitting `AuAuCentInputBDT` for `run28_embeddedPhoton12`. Submission
  is not complete yet; keep the terminal running and do not start manual
  merges until the driver finishes all 6 YAMLs and the raw workers later drain.
- 2026-05-11 22:19 EDT `sphnxuser04` merge recovery: after the fixed
  width-window WP0.50 analysis-only workers had finished, the manual merge was
  accidentally started with fragile chained commands and the wrong SIM input
  override name. Justin intentionally ran `condor_rm patsfan753` on the
  `sphnxuser04` schedd; visible terminal evidence immediately afterward showed
  `0` remaining `patsfan753` jobs. This removed active Condor merge DAGs/jobs,
  not the already written raw worker ROOT files under the campaign bulk roots.
  Local recovery patch now staged: `scripts/submit_auau_bdt_widthstudy_windows_wp050.sh`
  prints `MERGE_SIM_INPUT_BASE_OVERRIDE` for SIM merges, and new helper
  `scripts/merge_auau_bdt_widthstudy_windows_wp050_staged.sh` runs exactly one
  stage at a time (`firstRound`, `secondRound`, or `finalStitch`) so later
  stages cannot race ahead of missing first-round outputs. Local syntax and
  `git diff --check` passed. Uploaded both files to SDCC and verified `MATCH`
  with `./scripts/sftp_push_recoiljets.sh status ...`. Next checkpoint: run
  staged `firstRound` only on `sphnxuser04`, wait for it to drain, then run
  staged `secondRound`, then staged `finalStitch`.
- 2026-05-11 22:25 EDT `sphnxuser04` staged firstRound merge is running for
  campaign `widthstudy_windows_wp050_fixed_20260511_180500`. Justin launched
  `RJ_WIDTHSTUDY_CAMPAIGN_TAG=widthstudy_windows_wp050_fixed_20260511_180500
  RJ_WIDTHSTUDY_WINDOWS=5:35,10:35,15:35
  RJ_SIM_FIRSTROUND_REQUEST_MEMORY=8000MB RJ_SIM_MERGE_GROUP_SIZE=75 bash
  scripts/merge_auau_bdt_widthstudy_windows_wp050_staged.sh firstRound`.
  Visible terminal evidence shows the helper using the corrected raw input
  bases, e.g. `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded_widthstudy_windows_wp050_fixed_20260511_180500_pt5to35`
  and `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembeddedinclusive_widthstudy_windows_wp050_fixed_20260511_180500_pt5to35`,
  writing chunk outputs under
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_widthstudy_windows_wp050_fixed_20260511_180500_pt5to35/...`.
  It grouped each 1429-file sample into 20 firstRound jobs and submitted quiet
  DAG wrappers; observed clusters include `1140128` through `1140137` while
  scanning/queuing the `pt5to35` signal side and starting the inclusive side.
  Follow-up terminal evidence shows the `pt5to35` inclusive-background
  `embeddedJet20` side submitted quiet firstRound DAG cluster `1140146`, then
  the helper advanced to the `pt10to35` signal merge with the correct
  `simembedded_widthstudy_windows_wp050_fixed_20260511_180500_pt10to35` input
  and output root. Later visible terminal evidence shows `pt10to35` signal
  firstRound submitted cluster `1140161`, and `pt10to35` inclusive-background
  submitted at least cluster `1140167` while continuing through the reference
  cfg. Do not run secondRound until this firstRound queue drains cleanly.
- 2026-05-11 22:30 EDT `sphnxuser04` staged firstRound merge submission
  completed. Visible terminal queue after the helper returned showed `720`
  user jobs, all `hadd_condor.sh`, `720 idle`, `0 running`, `0 held`,
  `0 removed`; all-user query showed the same 720 idle plus 117 unrelated held
  jobs from other users. This is the expected firstRound merge layer waiting
  for slots, not a failure. Do not submit `secondRound` until these 720 jobs
  drain and no held jobs appear for `patsfan753`.
- 2026-05-11 22:15 EDT / visible terminal later switched to `sphnxuser06`:
  Justin resubmitted the expanded 5-40 AuAu tight-BDT validation after the
  diagnostic-column cache fix. Command used `RJ_NOTIFY_EMAILS=just0131@gmail.com`,
  `RJ_ML_PYTHON=/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python`,
  `./scripts/auau_tight_bdt_pipeline.sh validateOnSimCondor`, source
  `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049`,
  model dir
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_expanded_20260509_152604`,
  registry
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_expanded_20260509_152604/model_registry.json`,
  and `groupSize 100`. New report root:
  `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/model_validation_condor_20260511_221510`;
  submit dir:
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/condor_sub/auauTightBDTValidate_20260511_221510`;
  DAG cluster `1274500`; expected `8000` ROOT files, `80` shards, request
  memory `2500MB`. This is independent of the `sphnxuser04` merge recovery.
- 2026-05-11 22:25 EDT Gmail confirmed the above expanded BDT validation
  rerun reached `[RecoilJets][auauTightBDT_validateOnSimCondor][READY]` for
  `model_validation_condor_20260511_221510`. This unblocks deriving/regenerating
  target80 configs from the expanded model family. However, the visible
  `sphnxuser06` terminal then showed an already-live large RecoilJets queue:
  `22861` user jobs, `13236 idle`, `9625 running`, `0 held`, with
  `RecoilJets_Condor_AuAu.sh run28_embeddedJet20` entries such as cluster
  `1274595`. Before any further "full campaign" submit, identify what those
  live clusters belong to; do not duplicate-submit while this queue is active.
- 2026-05-11 22:41 EDT new BDT target80 prep started on `sphnxuser03` to keep
  it separate from the large MLP queue on `sphnxuser06`. Justin ran the staged
  prep command with tag `bdt_target80_all_available_20260511_224158`, output
  config dir
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_generated_configs/bdt_target80_all_available_20260511_224158`,
  `TARGET_SIGNAL_EFF=0.80`, `RJ_TARGET80_DO_SUBMIT=0`, expanded validation
  `model_validation_condor_20260511_221510`, width 15-30 validation
  `model_validation_condor_20260511_203438`, width-window validation
  `model_validation_condor_20260511_203440`, and et-fine validation
  `model_validation_condor_20260511_203441`. Visible terminal evidence showed
  `RECOILJETS_AUAU_BDT_TARGET80_AVAILABLE_PREP_V1`,
  `host=sphnxuser03.sdcc.bnl.gov`, and the first step deriving expanded
  5-40 target80 working points in `centpt` mode with centrality bins
  `0,20,50,80`. This is prep only; no full MC jobs should submit until the
  generated diagnostics/configs are inspected and the separate submit driver is
  run.
- 2026-05-11 22:5x EDT visible `sphnxuser03` terminal showed the above prep
  successfully wrote expanded target80 artifacts
  `bdt_working_points_target80.{json,yaml,csv}` and
  `bdt_working_points_target80_runtime_fragment.yaml` under validation report
  `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/model_validation_condor_20260511_221510`,
  then stopped during YAML generation because the expanded product map requested
  optional product `centAsFeat3x3_pt5to40`, which is not present in that
  expanded validation report. Local fix uploaded to SDCC:
  `scripts/prepare_auau_bdt_target80_available_campaigns.sh` now filters each
  `variant=product` map against the products actually present in the generated
  WP JSON, prints `[WARN] dropping unavailable product map entry...`, and
  continues preparing configs for available products. The helper also now reuses
  an existing `bdt_working_points_target80.json` unless
  `RJ_TARGET80_FORCE_DERIVE=1`, so rerunning after a late YAML-generation
  failure does not waste time re-merging all score caches. Checks passed locally:
  `bash -n scripts/prepare_auau_bdt_target80_available_campaigns.sh` and
  `git diff --check -- scripts/prepare_auau_bdt_target80_available_campaigns.sh`.
  Next: rerun the same prep command on `sphnxuser03` with the same tag/config
  directory to finish staging YAMLs, then inspect target80 diagnostics before
  any MC submission.
- 2026-05-11 later visible `sphnxuser03` terminal confirmed the rerun finished
  cleanly. Prepared config directory:
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_generated_configs/bdt_target80_all_available_20260511_224158`.
  The terminal showed generated target80 YAMLs for `etfine_15to35`,
  `widthstudy_pt5to35`, `widthstudy_pt10to35`, and `widthstudy_pt15to35`, with
  runtime entries `4`, `3`, `3`, and `3` respectively, and staged WP artifacts
  under `working_point_artifacts/<name>`. It also prepared width-window
  templates in
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_generated_configs/bdt_target80_all_available_20260511_224158/widthstudy_window_templates`.
  The script was run with `RJ_TARGET80_DO_SUBMIT=0`, so no MC production jobs
  were submitted by this prep step. Next checkpoint: inspect/pull the
  `bdt_working_point_target80_diagnostics.png` artifacts and run the target80
  config-directory local smoke/preflight before full paired MC submission.
- 2026-05-11 follow-up on `sphnxuser03`: first local-smoke command failed with
  `Permission denied` because `scripts/submit_auau_bdt_target80_config_dir.sh`
  was not executable on SDCC. Running it as `bash ./scripts/...` then reached
  the helper but failed with
  `./scripts/RecoilJets_Condor_submit.sh: No such file or directory`; the
  RecoilJets submitter is root-level on SDCC, not under `scripts/`. Codex
  patched and uploaded `scripts/submit_auau_bdt_target80_config_dir.sh` so the
  smoke path calls `${repo_root}/RecoilJets_Condor_submit.sh` and emits a clear
  error if the root submitter is missing. Local checks passed: `bash -n` and
  `git diff --check` for the helper. Next: rerun the same smoke command with
  `bash ./scripts/submit_auau_bdt_target80_config_dir.sh`.
- 2026-05-11 23:22 EDT target80 full paired MC submission started from
  `sphnxuser03` after Justin interrupted the slow local 1000-event smoke. Driver:
  `RJ_TARGET80_CONFIG_DIR=/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_generated_configs/bdt_target80_all_available_20260511_224158`,
  `RJ_TARGET80_SUBMIT_GROUP=bdt_target80_all_available_20260511_224158`,
  `RJ_TARGET80_RUN_LOCAL_SMOKE=0`, `RJ_TARGET80_DO_SUBMIT=1`,
  `bash ./scripts/submit_auau_bdt_target80_config_dir.sh`. The driver found
  6 YAMLs: `analysis_config_etfine_15to35_target80.yaml`,
  `analysis_config_expanded_5to40_target80.yaml`,
  `analysis_config_widthstudy_pt10to35_target80.yaml`,
  `analysis_config_widthstudy_pt1530_target80.yaml`,
  `analysis_config_widthstudy_pt15to35_target80.yaml`, and
  `analysis_config_widthstudy_pt5to35_target80.yaml`. Log file:
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_generated_configs/bdt_target80_all_available_20260511_224158/submit_bdt_target80_all_available_20260511_224158.log`.
  First pair in progress:
  `bdt_target80_all_available_20260511_224158_analysis_config_etfine_15to35_target80`,
  YAML `analysis_config_etfine_15to35_target80.yaml`, `auto_merge=0`,
  `groupSize=7`, worker memory `12000MB`, retry cap `16000`, merge memory
  `8000MB`. Bulk roots:
  `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded_bdt_target80_all_available_20260511_224158_analysis_config_etfine_15to35_target80`
  and expected inclusive partner under `simembeddedinclusive_...`; merge root:
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_bdt_target80_all_available_20260511_224158_analysis_config_etfine_15to35_target80`.
  Observed clusters so far for `isSimEmbedded`: `3032975` =
  `run28_embeddedPhoton12` reference row, `3032976` =
  `run28_embeddedPhoton20` reference row, `3032977` =
  `run28_embeddedPhoton12` `newPPG12` tight row, `3032978` =
  `run28_embeddedPhoton20` `newPPG12` tight row, `3032979` =
  `run28_embeddedPhoton12` `AuAuCentInputBase3x3BDT` row, and `3032980` =
  `run28_embeddedPhoton20` `AuAuCentInputBase3x3BDT` row. Submitter is still
  running through the remaining et-fine rows.
- 2026-05-11 23:34 EDT visible `sphnxuser04` terminal confirmed the
  width-window WP0.50 manual `firstRound` merge recovery drained cleanly:
  `condor_q` showed `Total for query: 0 jobs; 0 completed, 0 removed, 0 idle,
  0 running, 0 held, 0 suspended` for `patsfan753`. Just before draining, only
  two DAGMan parents and a few `hadd_condor.sh`/notify jobs remained, with
  `0` held. Gmail check found no fresh unread `RecoilJets Pipeline` email for
  this 04 merge. Next safe step is staged `secondRound` for campaign
  `widthstudy_windows_wp050_fixed_20260511_180500`; do not run `finalStitch`
  until secondRound drains to `0` jobs and `0` held.
- 2026-05-12 afternoon EDT `sphnxuser04` width-window WP0.50 merge
  inventory was rechecked through nested SSH after ProxyJump forwarding failed.
  Queue was clean for `patsfan753`: `0` idle/running/held. For campaign
  `widthstudy_windows_wp050_fixed_20260511_180500`, every window
  (`pt5to35`, `pt10to35`, `pt15to35`) and both datasets (`simembedded`,
  `simembeddedinclusive`) has exactly `6` flat secondRound sample-level ROOT
  files named `RecoilJets_<sample>_ALL_<cfg>.root`; all final stitched counts
  are still `0`. So worker output and secondRound are complete; only
  `finalStitch` is missing before offline pull. Earlier manual secondRound
  attempts failed because they discovered the wrong input base or because
  finalStitch could not discover flat secondRound files.
- 2026-05-12 follow-up local fix: `scripts/mergeRecoilJets.sh` was patched so
  SIM `secondRound` discovery accepts firstRound partial files
  `chunkMerge_<sample>_grp*.root` under the merge output tree instead of
  requiring raw sample subdirectories; SIM `finalStitch` discovery now accepts
  sample-level secondRound files `RecoilJets_<sample>_ALL_<cfg>.root`.
  `scripts/merge_auau_bdt_widthstudy_windows_wp050_staged.sh` was patched so
  `secondRound`/`finalStitch` use the stage output tree as input and pass
  `condor` automatically. Local checks passed: `bash -n` for both files and
  `git diff --check`. Upload via `./scripts/sftp_push_recoiljets.sh
  mergeRecoilJets.sh scripts/merge_auau_bdt_widthstudy_windows_wp050_staged.sh`
  failed twice because local DNS could not resolve `sftp.sdcc.bnl.gov`; retry
  upload when SDCC DNS is reachable before rerunning the staged merge helper.
- 2026-05-11 fixed the expanded BDT validation merge failure that blocked
  target-80 generation. Diagnosis from `sphnxuser06` report
  `model_validation_condor_20260511_215329`: all 80 worker shards completed,
  but MERGE failed with
  `Missing diagnostic column cluster_weta33_cogx in score cache ... score_cache_00001.npz`.
  Local patch: `scripts/validate_auau_tight_bdt_on_sim.py` now always includes
  `DIAGNOSTIC_FEATURES` in the worker-side `collect_rows` branch list, so score
  caches contain `cluster_weta33_cogx` and `cluster_wphi33_cogx` even when a
  product family's feature list would not otherwise force them. Checks passed:
  `PYTHONPYCACHEPREFIX=/private/tmp/pycache-thesis python3 -m py_compile
  scripts/validate_auau_tight_bdt_on_sim.py`, `bash -n
  scripts/auau_tight_bdt_pipeline.sh scripts/sftp_push_recoiljets.sh`, and
  `git diff --check`. Uploaded the validator to SDCC with
  `./scripts/sftp_push_recoiljets.sh scripts/validate_auau_tight_bdt_on_sim.py`.
  Next checkpoint: rerun only the expanded 5-40 `validateOnSimCondor` on
  `sphnxuser06`, then derive target80 WPs from the new READY report.
- 2026-05-11 target-80 full-matrix submit helper added locally and uploaded to
  SDCC: `scripts/submit_auau_bdt_target80_config_dir.sh`. Local checks passed:
  `bash -n` for the helper, target-WP pair submitter, target80 prep helper, and
  SFTP uploader; `git diff --check` passed. The helper takes a generated
  target80 config directory, verifies each YAML has
  `auau_tight_bdt_working_point_entries`, optionally runs one local
  `isSimEmbedded` smoke, then sequentially submits every
  `analysis_config_*_target80.yaml` using
  `scripts/submit_auau_bdt_targetwp_pair.sh` with isolated campaign tags. This
  is the intended one-command submission path after the expanded target80 YAML
  is generated and inspected.
- 2026-05-11 21:59 EDT Gmail reported
  `[RecoilJets][auauTightMLP_validateOnSimCondor][READY]` from
  `sphnxuser06` for the environment-fixed full-file MLP validation submitted
  with stamp `primary_full_envfix_20260511_215801`; Codex consumed the email
  and marked it read. This is separate from the expanded BDT target-80
  prerequisite validation. Next checkpoint for the MLP lane: pull/inspect the
  MLP validation report before drawing physics conclusions or promoting it to
  any RecoilJets MC production test.
- 2026-05-11 21:51 EDT terminal evidence on `sphnxuser04` showed the fixed
  width-window WP0.50 RecoilJets analysis-only worker campaign drained to
  `0` user jobs, `0` idle, `0` running, `0` held, and `0` removed in the live
  query. This was campaign
  `widthstudy_windows_wp050_fixed_20260511_180500`, worker clusters
  `1140092`-`1140103`. Because the fixed wrapper intentionally used
  `RJ_WIDTHSTUDY_AUTO_MERGE=0`, the drained queue means the raw analysis worker
  layer is complete; it does not mean final merged SIM ROOTs exist yet. The
  next required step is to run the wrapper-printed manual
  `mergeRecoilJets.sh firstRound/secondRound/finalStitch` commands for each
  pT-window/dataset pair, then pull the final merged outputs. The associated
  9-model width-window validation report is already available locally at
  `/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auauTightBDTValidation/model_validation_condor_20260511_171943`.
- 2026-05-11 21:58 EDT MLP full-file Condor validation was resubmitted from
  `sphnxuser06` after the MLP worker `LD_LIBRARY_PATH` fix. Driver header:
  `RECOILJETS_AUAU_MLP_PRIMARY_FULL_VALIDATE_ENVFIX_SUBMIT_V1`; stamp:
  `primary_full_envfix_20260511_215801`; source:
  `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049`;
  model dir:
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/tight_mlp_primary_competition_smoke_20260511_210945`;
  report:
  `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/mlp_model_validation_condor_primary_full_envfix_20260511_215801`;
  DAG:
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/condor_sub/auauTightMLPValidate_primary_full_envfix_20260511_215801/auau_tight_mlp_validateOnSimCondor.dag`;
  DAGMan cluster `1274415`; `8000` ROOT files, `80` shards,
  `groupSize=100`, `scoreMaxPerShard=5000`, request memory `3000MB`. Gmail
  check immediately after submission found no new unread MLP READY/CHECK/FAILED
  email yet. At 22:02 EDT, visible terminal `condor_q` on `sphnxuser06` showed
  `0` `patsfan753` jobs and `0` held; Gmail searches over unread and recent
  RecoilJets pipeline messages still found no fresh MLP READY/CHECK/FAILED
  message. Next checkpoint: inspect the report root and merge/worker logs on
  SDCC to determine READY/CHECK status.
- 2026-05-11 22:02 EDT MLP full-file Condor validation env-fix report was
  inspected in the visible SDCC terminal and is READY. Report:
  `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/mlp_model_validation_condor_primary_full_envfix_20260511_215801`.
  The report has `score_caches=80`, `status=READY`, `scored_entries=400000`,
  `signal_entries=8526032`, `background_entries=3474284`,
  `finite_fraction=1.0`, AUC `0.8510415464617735`, WP80 threshold
  `0.48582918047904966`, WP80 fake rate `0.20916062045973607`, and
  score-vs-Eiso correlation `-0.3566281184596356`. The summary printed
  `scored_signal_entries=0` and `scored_background_entries=0` despite finite
  AUC/WP metrics, so this should be treated as a merge-summary bookkeeping
  oddity to inspect before freezing production YAML. Next checkpoint: inspect
  `validation_metrics.json`/working-point JSON and derive WP80 YAML if the
  class counts inside the metrics/caches are sane.
- 2026-05-11 22:10 EDT score-cache label audit for the READY MLP full-file
  validation passed on SDCC. For report
  `mlp_model_validation_condor_primary_full_envfix_20260511_215801`, cache
  audit printed `cache_files=80`, `cache_label_entries=400000`,
  `cache_signal_entries=177133`, `cache_background_entries=222867`, AUC
  `0.8510415464617735`, finite fraction `1.0`, WP80
  `signal_entries=177133`, `background_entries=222867`,
  `signal_efficiency=0.7999977418098265`, `background_fake_rate=0.20916062045973607`,
  threshold `0.48582918047904966`, and runtime entry
  `auauCentInputBase3x3MLP|binned|15;20;25;35|0.4994844258;0.4940122843;0.5239367485|15|35|1`.
  The subsequent WP-derivation command failed because `RJ_ML_PYTHON` was set
  but not exported, so the child pipeline script used default `python3` and
  failed to import `uproot`; with `set -e` active, the SSH shell exited back to
  the local Mac. Next checkpoint: reconnect to SDCC and run WP80 generation via
  a temporary script that exports `RJ_ML_PYTHON`.
- 2026-05-11 MLP WP80 YAML generation succeeded on `sphnxuser06`. Generated
  YAML:
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_generated_configs/mlp_primary_full_envfix_20260511_215801_wp080/analysis_config_mlp_primary_full_envfix_20260511_215801_wp080.yaml`.
  The YAML contains model dir
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/tight_mlp_primary_competition_smoke_20260511_210945`
  and runtime WP entry
  `auauCentInputBase3x3MLP|binned|15;20;25;35|0.4994844258;0.4940122843;0.5239367485|15|35|1`.
  Next checkpoint: submit paired analysis-only `isSimEmbedded` and
  `isSimEmbeddedInclusive` MLP MC with isolated output roots and explicit memory.
- 2026-05-11 22:13 EDT paired MLP WP80 MC submission started from
  `sphnxuser06` with driver `/tmp/submit_mlp_wp80_pair.sh`, header
  `RECOILJETS_AUAU_MLP_TARGETWP_PAIR_SUBMIT_V1`. Campaign tag:
  `mlp_primary_full_envfix_wp080_20260511_221343`; YAML:
  `condor_generated_configs/mlp_primary_full_envfix_20260511_215801_wp080/analysis_config_mlp_primary_full_envfix_20260511_215801_wp080.yaml`;
  `groupSize=7`, `RJ_REQUEST_MEMORY=12000MB`,
  `RJ_AUTO_MEMORY_RETRY_CAP_MB=16000`, `RJ_SIM_FIRSTROUND_REQUEST_MEMORY=8000MB`,
  `RJ_MLP_TARGETWP_AUTO_MERGE=0` (analysis-only), notify
  `just0131@gmail.com`. Signal bulk root:
  `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded_mlp_primary_full_envfix_wp080_20260511_221343`;
  merge root:
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_mlp_primary_full_envfix_wp080_20260511_221343`;
  signal namespace `simembedded_20260511_221344`;
  snapshot dir
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_snapshots/simembedded_20260511_221344`;
  YAML/fanout dir
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_yaml_overrides/simembedded_condorDoAll_20260511_221344`.
  Visible clusters so far: `isSimEmbedded`
  `preselectionNewPPG12_tightReference_nonTightReference_baseVariant`
  `run28_embeddedPhoton12` `1429` jobs -> cluster `1274497` and
  `run28_embeddedPhoton20` `1429` jobs -> cluster `1274498`;
  `preselectionNewPPG12_tightNewPPG12_nonTightReference_baseVariant`
  `run28_embeddedPhoton12` `1429` jobs -> cluster `1274499` and
  `run28_embeddedPhoton20` `1429` jobs -> cluster `1274581`;
  `preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant`
  `run28_embeddedPhoton12` `1429` jobs -> cluster `1274582`.
  Follow-up terminal check at 22:18 EDT showed the matching
  `run28_embeddedPhoton20` BDT-reference pass submitted as cluster `1274584`.
  Follow-up terminal check at 22:20 EDT confirmed the MLP WP80 row:
  `preselectionNewPPG12_tightAuAuCentInputBase3x3MLP_nonTightAuAuMLPComplement_baseVariant`
  `run28_embeddedPhoton12` `1429` jobs -> cluster `1274585` and
  `run28_embeddedPhoton20` `1429` jobs -> cluster `1274586`. The driver then
  completed the signal-side submission and moved to
  `isSimEmbeddedInclusive`. Inclusive bulk root:
  `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembeddedinclusive_mlp_primary_full_envfix_wp080_20260511_221343`;
  inclusive namespace `simembeddedinclusive_20260511_222048`; snapshot dir
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_snapshots/simembeddedinclusive_20260511_222048`;
  YAML/fanout dir
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_yaml_overrides/simembeddedinclusive_condorDoAll_20260511_222048`.
  First inclusive cluster visible:
  `preselectionNewPPG12_tightReference_nonTightReference_baseVariant`
  `run28_embeddedJet12` `1429` jobs -> cluster `1274587` and
  `run28_embeddedJet20` `1429` jobs -> cluster `1274588`.
  Follow-up terminal check showed the second inclusive reference row
  `preselectionNewPPG12_tightNewPPG12_nonTightReference_baseVariant`
  completed: `run28_embeddedJet12` `1429` jobs -> cluster `1274589` and
  `run28_embeddedJet20` `1429` jobs -> cluster `1274590`. The submitter then
  completed the inclusive BDT-reference row
  `preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant`:
  `run28_embeddedJet12` `1429` jobs -> cluster `1274592` and
  `run28_embeddedJet20` `1429` jobs -> cluster `1274593`. The submitter then
  completed the inclusive MLP WP80 row; terminal tail showed
  `preselectionNewPPG12_tightAuAuCentInputBase3x3MLP_nonTightAuAuMLPComplement_baseVariant`
  `run28_embeddedJet20` entries under cluster `1274595`. The final footer
  printed `TRACK_THIS_CAMPAIGN=mlp_primary_full_envfix_wp080_20260511_221343`,
  `GENERATED_YAML=condor_generated_configs/mlp_primary_full_envfix_20260511_215801_wp080/analysis_config_mlp_primary_full_envfix_20260511_215801_wp080.yaml`,
  and submit log
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_generated_configs/mlp_primary_full_envfix_wp080_20260511_221343/submit_mlp_primary_full_envfix_wp080_20260511_221343.log`.
  Final queue summary at submission end: `22864` jobs matching the campaign,
  `16808` idle, `6056` running, `0` held, `0` removed. Gmail label
  `RecoilJets Pipeline` had `0` unread messages immediately after submission.
  Follow-up user-run SDCC snapshot at about 22:25 EDT showed the campaign still
  healthy: `22863` matching jobs, `16010` idle, `6853` running, `0` held,
  `0` removed. Follow-up terminal status at about 22:30 EDT showed continued
  healthy scaling: `22863` jobs, `13238` idle, `9625` running, `0` held,
  `0` removed; Gmail `RecoilJets Pipeline` had `0` unread messages. Next
  checkpoint: parse the submit log for the exact inclusive
  MLP `run28_embeddedJet12` cluster, then watch queue/emails until worker layer
  drains; auto-merge is off, so manual merge commands will be needed after
  workers finish.
- 2026-05-11 about 22:35 EDT the partial first-pass MLP WP80 MC campaign
  `mlp_primary_full_envfix_wp080_20260511_221343` was intentionally removed
  from `sphnxuser06` before merge/output promotion because the model had only
  been trained on the competition-smoke 300-files-per-sample subset. Visible
  terminal evidence after removal showed jobs in Condor state `X` and summary:
  `19034` removed, `0` idle, `0` running, `0` held for `patsfan753`; Gmail
  `RecoilJets Pipeline` had `0` unread messages. Treat this campaign as
  canceled/runtime-smoke-only. Next checkpoint: start full-stat primary MLP
  training in a fresh model directory, then validate full-file before launching
  a new paired MC campaign.
- 2026-05-11 22:44 EDT full-stat primary MLP training started in `tmux` on
  `sphnxuser06`. Session: `mlp_full_primary`; script:
  `/tmp/train_mlp_full_primary.sh`; source:
  `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049`;
  model dir:
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/tight_mlp_full_primary_20260511_224401`;
  log:
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/train_full_primary_20260511_224401.log`.
  Training config: `PRODUCTS=primary`, `PT_RANGE=15:35`,
  `MAX_FILES_PER_SAMPLE=0`, `MAX_ROWS=0`, `MAX_ROWS_PER_CLASS=0`,
  `EPOCHS=120`, `PATIENCE=18`, `RESTARTS=3`, hidden grid
  `64,32;96,48;64,32,16;128,64,32`, selection metric `wp_fake_rate` at signal
  efficiency `0.80`, thread knobs set to `4`. Terminal evidence showed
  `files_before=8000 files_after=8000` with all four groups at `2000` files and
  the trainer reading `1/8000`, confirming this is full-stat, not the previous
  300-file-per-sample smoke. Next checkpoint: monitor the log for row-read
  summary, candidate training progress, selected candidate, `applyCheck`, and
  final `DONE_FULL_MLP_MODEL_DIR=...`.
- 2026-05-11 22:25 EDT Gmail reported
  `[RecoilJets][auauTightBDT_validateOnSimCondor][READY]` from `sphnxuser06`
  for the expanded BDT validation lane. Codex consumed the email and marked it
  read. This is separate from the active MLP WP80 MC campaign, but it likely
  unblocks the parallel expanded-BDT target-80 working-point lane.
- 2026-05-11 21:53 EDT expanded AuAu tight-BDT validation was rerun from
  `sphnxuser06` to regenerate modern score caches for target-80 working-point
  derivation. Reason: target-80 prep against the older expanded report
  `model_validation_condor_20260509_192942` failed because
  `score_cache_00001.npz` was missing the updated diagnostic column
  `cluster_weta33_cogx`. Submit command used
  `RJ_NOTIFY_EMAILS=just0131@gmail.com`,
  `RJ_ML_PYTHON=/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python`,
  `./scripts/auau_tight_bdt_pipeline.sh validateOnSimCondor`,
  source
  `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049`,
  model dir
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_expanded_20260509_152604`,
  registry
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_expanded_20260509_152604/model_registry.json`,
  and `groupSize 100`. New report root:
  `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/model_validation_condor_20260511_215329`;
  submit dir:
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/condor_sub/auauTightBDTValidate_20260511_215329`;
  DAG:
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/condor_sub/auauTightBDTValidate_20260511_215329/auau_tight_bdt_validateOnSimCondor.dag`;
  score caches:
  `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/model_validation_condor_20260511_215329/score_caches`;
  DAGMan cluster `1274333`; `8000` ROOT files, `80` shards,
  request memory `2500MB`. Next checkpoint: wait for READY/CHECK/FAILED email
  or queue evidence, then derive target-80 WPs from this new report, generate
  `analysis_config_expanded_5to40_target80.yaml` in
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_generated_configs/target80_all_available_20260511_212451`,
  inspect the target-80 fit diagnostics in chat, and only then submit paired
  `isSimEmbedded`/`isSimEmbeddedInclusive` MC. Local pull attempt at 22:00 EDT
  found the remote report path not visible through `sftp.sdcc.bnl.gov`; if the
  SDCC terminal shows the report exists, copy or tar the report to a GPFS/SFTP-
  visible transfer directory before pulling. Follow-up terminal evidence at
  22:03 EDT showed `condor_q` drained to zero on `sphnxuser06`, but
  `cat "$expval/validation_summary.txt"` returned `No such file or directory`.
  Treat `model_validation_condor_20260511_215329` as `CHECK / no summary` until
  DAGMan and worker logs explain whether shards failed, merge failed, or the
  report directory was never finalized.
- 2026-05-11 MLP full-file Condor validation of the current competition-smoke
  artifact was submitted from `sphnxuser06` and immediately returned
  `[CHECK]`. Submit stamp `primary_full_20260511_214302`; DAGMan cluster
  `1274171`; DAG
  `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/condor_sub/auauTightMLPValidate_primary_full_20260511_214302/auau_tight_mlp_validateOnSimCondor.dag`;
  report
  `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/mlp_model_validation_condor_primary_full_20260511_214302`.
  Pipeline email at 21:44:38 reported `status=CHECK`,
  `expected_score_caches=80`, `found_score_caches=0`,
  `notes=missing score cache shards`. Next checkpoint: inspect the DAGMan log
  and worker submit/log/error tails on SDCC before resubmitting or modifying the
  validator.
- 2026-05-11 MLP Condor validation root cause fixed locally and uploaded to
  SDCC: all `validate_*.err` files from
  `primary_full_20260511_214302` showed
  `/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python: error while loading shared libraries: libpython3.13.so.1.0`.
  `scripts/auau_tight_mlp_pipeline.sh` now writes generated
  `validate_worker.sh` and `validate_merge.sh` scripts that set up the
  sPHENIX stack, prepend the ML venv `lib`/`lib64` directories to
  `LD_LIBRARY_PATH`, unset `PYTHONHOME`, and call `"$ml_python"`.
  Verification: local `bash -n`, dry-run Condor script generation, generated
  worker/merge `bash -n`, and SFTP remote status `MATCH` for
  `scripts/auau_tight_mlp_pipeline.sh`. Next checkpoint: resubmit the full-file
  MLP validation with a fresh stamp and watch for score caches/READY email.
- 2026-05-11 target-80 rerun automation was hardened and uploaded to SDCC.
  Updated scripts:
  `scripts/prepare_auau_bdt_target80_available_campaigns.sh` and
  `scripts/submit_auau_bdt_widthstudy_windows_wp050.sh`. Static checks passed:
  `bash -n` for both scripts, `git diff --check`, and Python compile of the
  working-point/config helpers with a writable pycache. The prep helper now
  supports campaign-specific target-WP pT bins and includes the width-window
  validation family (`pt5to35`, `pt10to35`, `pt15to35`) by generating
  prep-only YAML templates from the existing width-window wrapper. It can stage
  frozen target-80 YAMLs for validated available families: 15-30 width study,
  5/10/15-35 width-window study, and fine-`E_T` 15-35 study. The expanded
  5-40 validation report `203437` remains excluded unless it is diagnosed or
  rerun because it had no summary. Next checkpoint: on SDCC rebuild `src_AuAu`,
  run the target-80 prep block with validation paths `203438`, `203440`, and
  `203441`, run one local runtime check per generated YAML, then submit paired
  `isSimEmbedded`/`isSimEmbeddedInclusive` target-80 campaigns only after the
  generated configs and diagnostics are confirmed.
- 2026-05-11 target-80 config staging completed on `sphnxuser06` from
  `/sphenix/u/patsfan753/scratch/thesisAnalysis`. Staging tag:
  `target80_all_available_20260511_212451`; config directory:
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_generated_configs/target80_all_available_20260511_212451`.
  Generated frozen YAMLs:
  `analysis_config_widthstudy_pt1530_target80.yaml`,
  `analysis_config_widthstudy_pt5to35_target80.yaml`,
  `analysis_config_widthstudy_pt10to35_target80.yaml`,
  `analysis_config_widthstudy_pt15to35_target80.yaml`, and
  `analysis_config_etfine_15to35_target80.yaml`. Each YAML has copied
  working-point artifacts under
  `condor_generated_configs/target80_all_available_20260511_212451/working_point_artifacts/<campaign>`.
  After model-matrix audit, this set is known to be incomplete for the original
  expanded comparison axes: no-centrality-input, global centrality-input,
  minority-balanced centrality-input, pT-binned centrality-input, and
  pT×centrality BDTs are in the older good expanded validation report
  `model_validation_condor_20260509_192942`, not in the currently staged five
  YAMLs. Add an expanded target80 YAML from `20260509_192942` before declaring
  the MC target80 campaign matrix complete. Running trace:
  `agent_context/AUAU_BDT_MODEL_VARIANT_TRACE.md`. Next checkpoint: rebuild
  `src_AuAu` if not already rebuilt after target80 runtime changes, generate
  the expanded target80 config, run local runtime smoke over all YAMLs, then
  submit paired MC reruns only if the smoke confirms finite scores and target80
  working-point lookup.
- 2026-05-11 JSTG was postponed by one week. The near-term strategy is now a
  focused AuAu tight-BDT development sprint: use the extra runway to validate
  the fine-`E_T`, width-window, and available expanded model registries; derive
  and visually inspect target-80 working-point fits before any full
  `isSimEmbedded`/`isSimEmbeddedInclusive` submissions; then run approved MC
  campaigns and build the offline comparison suite for ID/reco efficiency,
  background leakage, isolation/ABCD behavior, shower-shape templates,
  centrality/`E_T` stability, and `x_{J#gamma}` stability. Do not rush full MC
  from a new working-point config until the diagnostics have been pulled and
  inspected in chat.
- 2026-05-11 fixed width-window WP0.50 RecoilJets MC campaign
  `widthstudy_windows_wp050_fixed_20260511_180500` was clean-resubmitted from
  `sphnxuser04` after removing all old `patsfan753` Condor jobs on that submit
  host. The fixed wrapper ran in analysis-only mode
  (`RJ_WIDTHSTUDY_AUTO_MERGE=0`), requested `RJ_REQUEST_MEMORY=10000MB`,
  `RJ_AUTO_MEMORY_RETRY_CAP_MB=16000`, and printed explicit merge commands for
  each pT window instead of launching the fragile parent auto-merge DAG. It
  submitted all three windows (`pt5to35`, `pt10to35`, `pt15to35`) for both
  `isSimEmbedded` and `isSimEmbeddedInclusive`, with four worker clusters per
  window. Terminal evidence shows the final cluster span `1140092`-`1140103`
  with `1429` jobs per cluster. Final queue snapshot: `24358` jobs in the
  query, `7210 removed` from the prior cleanup, `11674 idle`, `5474 running`,
  and `0 held` for the current user query. Track log:
  `/sphenix/u/patsfan753/scratch/thesisAnalysis/submit_widthstudy_windows_wp050_fixed_20260511_180500.log`.
  Gmail check immediately after submission found two fresh unread
  `[RecoilJets][auto_simembeddedinclusive_final_ready][FAILED]` messages at
  18:04 and 18:09, but those are stale auto-DAG notifications from the
  removed pre-fix workflow because the new fixed run is analysis-only and does
  not launch `auto_*_final_ready` parent DAGs. Codex marked those stale
  consumed emails read. Next checkpoint: watch for fresh RecoilJets emails and,
  once all worker clusters leave the queue with no held jobs, run the printed
  manual merge commands per window.

| Status | Item | Evidence | Blocker / Risk | Next Command Or Action |
|---|---|---|---|---|
| 🟢 READY / Full-File Validation Passed | AuAu tight-MLP full-file Condor validation of competition-smoke artifact on `sphnxuser06` | Visible terminal evidence on 2026-05-11 from `/sphenix/u/patsfan753/scratch/thesisAnalysis` on `sphnxuser06`: Justin ran `/tmp/submit_mlp_primary_full_validation_envfix.sh`, header `RECOILJETS_AUAU_MLP_PRIMARY_FULL_VALIDATE_ENVFIX_SUBMIT_V1`; `host=sphnxuser06.sdcc.bnl.gov`; stamp `primary_full_envfix_20260511_215801`; source `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049`; model dir `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/tight_mlp_primary_competition_smoke_20260511_210945`; `score_max_rows=400000`; `request_memory=3000MB`; `group_size=100`; report `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/mlp_model_validation_condor_primary_full_envfix_20260511_215801`; DAG `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/condor_sub/auauTightMLPValidate_primary_full_envfix_20260511_215801/auau_tight_mlp_validateOnSimCondor.dag`; DAGMan cluster `1274415`; `8000` ROOT files, `80` shards, `scoreMaxPerShard=5000`. Post-submit queue showed DAGMan parent running; later `condor_q` showed `0` `patsfan753` jobs and `0` held. Report inspection showed `score_caches=80`, `status=READY`, `scored_entries=400000`, `signal_entries=8526032`, `background_entries=3474284`, product `centInputBase3x3MLP_pt1535` AUC `0.8510415464617735`, finite fraction `1.0`, WP80 threshold `0.48582918047904966`, WP80 fake rate `0.20916062045973607`, and score-vs-Eiso correlation `-0.3566281184596356`. | Full-file validation is green and close to the competition smoke metrics. The summary reports `scored_signal_entries=0` and `scored_background_entries=0` even though AUC/WP metrics are finite, so inspect metrics JSON/cache label counts before freezing the production YAML. | Inspect `validation_metrics.json` and `mlp_working_points*.json`; if class counts are sane, derive/generate the WP80 YAML from this report and run the runtime parity/apply check before paired `isSimEmbedded`/`isSimEmbeddedInclusive` MLP MC submission. |
| 🟢 Cache Labels Verified / WP80 YAML Pending | AuAu tight-MLP WP80 config generation from full-file validation | Visible terminal evidence on 2026-05-11 from `sphnxuser06`: cache audit over `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/mlp_model_validation_condor_primary_full_envfix_20260511_215801/score_caches.list` printed `cache_files=80`, `cache_label_entries=400000`, `cache_signal_entries=177133`, `cache_background_entries=222867`, AUC `0.8510415464617735`, finite fraction `1.0`, WP80 threshold `0.48582918047904966`, WP80 signal efficiency `0.7999977418098265`, WP80 fake rate `0.20916062045973607`, and runtime entry `auauCentInputBase3x3MLP|binned|15;20;25;35|0.4994844258;0.4940122843;0.5239367485|15|35|1`. | The class-label audit resolves the summary `0/0` concern. The WP-generation block failed only because `RJ_ML_PYTHON` was not exported before calling the pipeline script, causing default `python3` to lack `uproot`; the SSH shell exited because `set -e` was active. | Reconnect to SDCC and run a temporary WP80 generation script with `export RJ_ML_PYTHON=/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python`, then grep the generated YAML for the MLP model dir and working-point entries. |
| 🟢 WP80 YAML Ready | AuAu tight-MLP WP80 config from full-file validation | Visible terminal evidence on 2026-05-11 from `sphnxuser06`: `deriveWorkingPointsFromValidation` merged all 80 score caches and wrote `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/mlp_model_validation_condor_primary_full_envfix_20260511_215801/mlp_working_points_target80.{json,yaml}`. `generateWorkingPointConfig` then wrote `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_generated_configs/mlp_primary_full_envfix_20260511_215801_wp080/analysis_config_mlp_primary_full_envfix_20260511_215801_wp080.yaml`. Grep confirmed row `[newPPG12, auauCentInputBase3x3MLP, auauMLPComplement]`, `auau_tight_mlp_model_dir: /gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/tight_mlp_primary_competition_smoke_20260511_210945`, and `auau_tight_mlp_working_point_entries: ["auauCentInputBase3x3MLP|binned|15;20;25;35|0.4994844258;0.4940122843;0.5239367485|15|35|1"]`. | Ready for production-like MC submission. Keep output roots isolated from BDT products; submit analysis-only first so workers remain inspectable before merge. | Submit paired `isSimEmbedded` and `isSimEmbeddedInclusive` with `RJ_MLP_TARGETWP_CONFIG_YAML=condor_generated_configs/mlp_primary_full_envfix_20260511_215801_wp080/analysis_config_mlp_primary_full_envfix_20260511_215801_wp080.yaml`, `RJ_MLP_TARGETWP_AUTO_MERGE=0`, `RJ_REQUEST_MEMORY=12000MB`, and explicit campaign tag. |
| 🟡 Submitting / Signal Clusters Started | MLP WP80 paired MC campaign `mlp_primary_full_envfix_wp080_20260511_221343` | Visible terminal evidence on 2026-05-11 from `sphnxuser06`: driver `/tmp/submit_mlp_wp80_pair.sh` launched `scripts/submit_auau_mlp_targetwp_pair.sh` with YAML `condor_generated_configs/mlp_primary_full_envfix_20260511_215801_wp080/analysis_config_mlp_primary_full_envfix_20260511_215801_wp080.yaml`, `groupSize=7`, worker memory `12000MB`, retry cap `16000`, merge memory `8000MB`, `auto_merge=0`, notify `just0131@gmail.com`. Signal bulk root `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded_mlp_primary_full_envfix_wp080_20260511_221343`; merge root `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_mlp_primary_full_envfix_wp080_20260511_221343`; `isSimEmbedded` namespace `simembedded_20260511_221344`; snapshot dir `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_snapshots/simembedded_20260511_221344`; fanout dir `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_yaml_overrides/simembedded_condorDoAll_20260511_221344`. Submitted `isSimEmbedded` `run28_embeddedPhoton12` `1429` jobs to cluster `1274497` and `run28_embeddedPhoton20` `1429` jobs to cluster `1274498`. | Submitter was still running when last checked and appeared to begin another `run28_embeddedPhoton12` pass. Need the rest of the terminal output before calling the pair fully submitted. | Continue watching terminal. Record any additional signal clusters, then the `isSimEmbeddedInclusive` bulk root/clusters, final post-submit queue, and printed manual merge commands. |
| 🟡 Fixed / Ready To Resubmit | AuAu tight-MLP full-file Condor validation of competition-smoke artifact on `sphnxuser06` | Visible terminal evidence on 2026-05-11 from `/sphenix/u/patsfan753/scratch/thesisAnalysis` on `sphnxuser06`: Justin submitted `/tmp/submit_mlp_primary_full_validation.sh`, header `RECOILJETS_AUAU_MLP_PRIMARY_FULL_VALIDATE_SUBMIT_V1`; source `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049`; model dir `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/tight_mlp_primary_competition_smoke_20260511_210945`; `score_max_rows=400000`; `request_memory=3000MB`; `group_size=100`; submit stamp `primary_full_20260511_214302`; report `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/mlp_model_validation_condor_primary_full_20260511_214302`; DAG `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/condor_sub/auauTightMLPValidate_primary_full_20260511_214302/auau_tight_mlp_validateOnSimCondor.dag`; DAGMan cluster `1274171`. The queue showed DAGMan running at 21:43, then no `patsfan753` jobs by 21:46. Gmail `RecoilJets Pipeline` email at 21:44:38 reported `[RecoilJets][auauTightMLP_validateOnSimCondor][CHECK]` with `expected_score_caches=80`, `found_score_caches=0`, and `notes=missing score cache shards`; Codex consumed and marked that email read. Justin ran the focused diagnostic; all worker stderr files showed missing `libpython3.13.so.1.0`. Codex patched `scripts/auau_tight_mlp_pipeline.sh` so generated validation worker/merge scripts set the ML venv `LD_LIBRARY_PATH`, verified local syntax/dry-run generation, uploaded to SDCC with `SSH_AUTH_SOCK="$(launchctl getenv SSH_AUTH_SOCK)" ./scripts/sftp_push_recoiljets.sh auau_tight_mlp_pipeline.sh`, and confirmed remote `MATCH`. | The failed report is stale and should not be reused. The issue was a Condor worker environment failure, not MLP score quality. | Resubmit full-file validation with a fresh stamp using the fixed uploaded script; then watch for `expected_score_caches=80 found_score_caches=80` and a READY email/report. |
| 🟢 READY / Competition Smoke Passed | AuAu tight-MLP competition smoke on `sphnxuser06` | Visible terminal evidence on 2026-05-11 from `/sphenix/u/patsfan753/scratch/thesisAnalysis` on `sphnxuser06`: Justin launched `/tmp/run_mlp_competition_smoke.sh`, logging with `tee` to `/tmp/run_mlp_competition_smoke_*.log`. Driver header `RECOILJETS_AUAU_MLP_COMPETITION_SMOKE_V1`; source `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049`; model dir `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/tight_mlp_primary_competition_smoke_20260511_210945`; validation report `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/mlp_validation_competition_smoke_20260511_210945`; Python `/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python`; products `primary`; pT range `15:35`; max files/sample `300`; max rows `200000`; max rows/class `100000`; restarts `3`; hidden grid `64,32;96,48;64,32,16`; selection metric `wp_fake_rate` at signal efficiency `0.80`. Trainer printed balanced manifest `files_before=8000 files_after=1200` with groups `background:run28_embeddedJet12=300`, `background:run28_embeddedJet20=300`, `signal:run28_embeddedPhoton12=300`, `signal:run28_embeddedPhoton20=300`, then selected `centInputBase3x3MLP_pt1535:a4r2` with hidden layers `[64,32,16]`, `val_auc=0.69704`, `test_auc=0.69467`, and selection WP80 fake rate `0.4581393138838734`. Artifact `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/tight_mlp_primary_competition_smoke_20260511_210945/auau_tight_mlp_centInputBase3x3_pt1535.json`; registry `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/tight_mlp_primary_competition_smoke_20260511_210945/model_registry.json`; `applyCheck` passed (`loaded and scored 1 MLP JSON artifacts`). Validation used the same 1200-file balanced manifest, scored `239997` entries (`106410` signal, `133587` background), wrote score cache/metrics/WP files, and ended `status=READY` with finite fraction `1.0`, validation AUC `0.8498507602114719`, WP80 threshold `0.4843840599060058`, WP80 fake rate `0.211397815655715`, and score-vs-Eiso correlation `-0.3535302621966963`. | This is a much cleaner and more honest smoke than the first sparse pass, but it is still a capped validation subset. Do not promote directly to MC production until the model-selection JSON, training/validation manifest summaries, and read summaries are inspected. | Inspect `auau_tight_mlp_centInputBase3x3_pt1535.model_selection.json`, `training_manifest_summary.json`, `training_read_summary.json`, and `validation_manifest_summary.json`. If sane, run full embedded validation: `RJ_NOTIFY_EMAILS=just0131@gmail.com RJ_ML_PYTHON=/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python ./scripts/auau_tight_mlp_pipeline.sh validateOnSimCondor SOURCE=/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049 MODEL_DIR=/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/tight_mlp_primary_competition_smoke_20260511_210945 groupSize 100`. |
| 🟢 READY / First Smoke Passed | AuAu tight-MLP primary smoke on `sphnxuser06` | Visible terminal evidence on 2026-05-11 from `/sphenix/u/patsfan753/scratch/thesisAnalysis` on `sphnxuser06`: after rebuilding `src_AuAu` cleanly with `make clean; makeProject`, Justin launched `/tmp/run_mlp_smoke_*.sh`. Driver header `RECOILJETS_AUAU_MLP_SMOKE_DRIVER_V1`; source `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049`; model dir `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/tight_mlp_primary_smoke_20260511_204704`; validation report `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/mlp_validation_smoke_20260511_204704`; Python `/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python`; products `primary`; pT range `15:35`. The run read all 8000 extraction ROOT files, trained `centInputBase3x3MLP_pt1535`, early-stopped at epoch 70 with `best_epoch=58`, wrote artifact `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/tight_mlp_primary_smoke_20260511_204704/auau_tight_mlp_centInputBase3x3_pt1535.json`, wrote `model_registry.json`, and passed `applyCheck` (`loaded and scored 1 MLP JSON artifacts`). Validation wrote score cache, metrics, WP files, and summary under the report dir with `status=READY`, `scored_entries=96000`, `finite_fraction=1.0`, `product.centInputBase3x3MLP_pt1535.auc=0.9300795397885159`, `wp080_threshold=0.4507956326007843`, `wp080_fake_rate=0.1010262086561732`, and `score_vs_eiso_corr=-0.42998208333272886`. Training summary printed `test_auc=0.73553` and `val_auc=0.72890`. Runtime WP fragment: `auauCentInputBase3x3MLP|binned|15;20;25;35|0.4891304076;0.4713326871;0.6083169103|15|35|1`. After this pass, Codex uploaded improved verbose/balanced smoke scripts to SDCC and verified `MATCH` for `scripts/train_auau_photon_mlp.py`, `scripts/validate_auau_tight_mlp_on_sim.py`, and `scripts/auau_tight_mlp_pipeline.sh`. Codex then enhanced the trainer again to run an internal architecture/restart tournament and select the exported JSON artifact by WP80 fake rate, with AUC/loss/calibration tie-breaks, while preserving the same C++ runtime JSON schema; SDCC upload/status check showed `MATCH` for updated `scripts/train_auau_photon_mlp.py` and `scripts/auau_tight_mlp_pipeline.sh`. | First-pass mechanics are green, but do not promote this smoke model directly to production: the training-vs-validation AUC difference and the old validation smoke's sparse per-file scoring need the new verbose composition summaries before physics interpretation. | Run the new competition smoke, inspect `*.model_selection.json`, `training_manifest_summary.json`, `training_read_summary.json`, and `validation_manifest_summary.json`, then launch `validateOnSimCondor` for the selected MLP artifact if the cleaner smoke remains READY. |
| 🟡 Target-80 Derived / 203440 Pull Workaround Needed | All-variant AuAu BDT registry validation sweep for target-80 inputs | Visible terminal evidence on 2026-05-11 from `/sphenix/u/patsfan753/scratch/thesisAnalysis` on `sphnxuser06`: Justin ran the one-block driver `/tmp/submit_all_auau_bdt_validations_target80_inputs.sh`, which printed `RECOILJETS_AUAU_BDT_ALL_VARIANT_VALIDATION_SUBMIT_V1`, `source=/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049`, `group_size=100`, and `TRACK_THIS_ALL_VARIANT_VALIDATION=20260511_203443`. The driver submitted registry-driven `validateOnSimCondor` jobs for available model dirs: expanded 5-40 (`tight_expanded_20260509_152604`), width-study 15-30 (`tight_centinput_widthstudy_pt1530_current`), width-study windows (`tight_centinput_widthstudy_windows_current`), and fine-`E_T` centrality study (`tight_etfine_centstudy_current`). Terminal queue snapshot after submission showed report/submit timestamp `model_validation_condor_20260511_203437` actively queued with validation workers `1273814`-`1273841` and DAGMan parents for later timestamps `20260511_203440` and `20260511_203441`; query summary was `80 jobs`, `80 idle`, `0 running`, `0 held`, `0 removed`, and all-user summary `110 jobs`, `30 running`. Later terminal snapshots at 20:37 and 20:42 showed `condor_q` on `sphnxuser06` drained to `0` user jobs and `0` held. Gmail search found three new unread `[RecoilJets][auauTightBDT_validateOnSimCondor][READY]` emails from `sphnxuser06` at 20:37:31, 20:37:34, and 20:38:26; Codex consumed and marked those three READY emails read. Justin ran the report-accounting diagnostic. It showed `model_validation_condor_20260511_203438` READY for `tight_centinput_widthstudy_pt1530_current`, `model_validation_condor_20260511_203440` READY for `tight_centinput_widthstudy_windows_current`, and `model_validation_condor_20260511_203441` READY for `tight_etfine_centstudy_current`, all with `finite_score_fraction=1` and `notes=none`. The first report, `model_validation_condor_20260511_203437`, printed `NO SUMMARY`, so the expanded 5-40 validation report did not finalize a summary in that report dir. Justin then ran target-80 derivation for `203438`, `203440`, and `203441` with `MODE=centpt` and `CENT_BINS=0,20,50,80`; terminal evidence showed final `json=`, `yaml=`, `runtime_fragment=`, and `csv=` lines through `203441`. Codex pulled `203438` and `203441` locally to `dataOutput/auauTightBDTValidation/model_validation_condor_20260511_203438` and `_203441`. Both contain `bdt_working_points_target80.{json,yaml,csv}`, `bdt_working_points_target80_runtime_fragment.yaml`, and `bdt_working_point_target80_diagnostics.png`; visual inspection shows the achieved-efficiency-minus-target panels are flat near zero, so the local cell-by-cell 80% calibration is working. Pulling `203440` twice via the direct `/sphenix/tg/.../reports/model_validation_condor_20260511_203440` path failed with SFTP `File not found`, despite the SDCC terminal seeing the READY summary; use a small SDCC-side copy/tar to a GPFS transfer directory for that one if needed. | The three newer focused validations have derived target-80 WPs. Local inspection is available for `203438` and `203441`. `203440` is validated and derived remotely but not pulled locally due SFTP namespace/path issue. The expanded 5-40 report `203437` remains unusable until diagnosed or rerun. | If using width-window `203440` plots, copy its key files to a GPFS transfer dir from SDCC and pull that lite dir. Before target-80 MC submission, generate frozen YAMLs from approved WPs, rebuild `src_AuAu` once for the `grid2d` runtime, and run one local RecoilJets runtime check. |
| 🟢 Pulled / Inspect Target-80 Fits | Fine-`E_T` 89-model AuAu BDT `validateOnSimCondor` on `sphnxuser06` | Visible terminal evidence on 2026-05-11 from `/sphenix/u/patsfan753/scratch/thesisAnalysis` on `sphnxuser06`: Justin submitted `RJ_NOTIFY_EMAILS=just0131@gmail.com RJ_ML_PYTHON=/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python ./scripts/auau_tight_bdt_pipeline.sh validateOnSimCondor SOURCE=/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049 MODEL_DIR=/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_etfine_centstudy_current MODEL_REGISTRY=/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_etfine_centstudy_current/model_registry.json groupSize 100`. Path plan: report `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/model_validation_condor_20260511_194832`; submit dir `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/condor_sub/auauTightBDTValidate_20260511_194832`; DAG `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/condor_sub/auauTightBDTValidate_20260511_194832/auau_tight_bdt_validateOnSimCondor.dag`; cache dir `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/model_validation_condor_20260511_194832/score_caches`; `root files=8000`, `shards=80`, `groupSize=100`, `scoreMaxPerShard=5000`, `request_memory=2500MB`; DAGMan cluster `1273669`. Gmail evidence at 2026-05-11 19:51 from `sphnxuser06` reported `[RecoilJets][auauTightBDT_validateOnSimCondor][READY]` with `total_entries=12000316`, `signal_entries=8526032`, `background_entries=3474284`, `scored_entries=400000`, finite score fraction `1`, report dir `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/model_validation_condor_20260511_194832`, and AUCs: `centInput_pt1535=0.773122`, `ptFine_centInput=0.776523`, `ptFine_cent3=0.78874`, `ptFine_cent7=0.808249`. Codex consumed and marked the email read; terminal queue also showed `0` user jobs on `sphnxuser06`. Justin then ran `deriveWorkingPointsFromValidation TARGET=0.80`, which wrote `bdt_working_points_target80.{json,yaml,csv}`, `bdt_working_points_target80_runtime_fragment.yaml`, and `bdt_working_point_target80_diagnostics.png`. Codex pulled the full report locally to `/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auauTightBDTValidation/model_validation_condor_20260511_194832`. Current target-80 entries are all linear: `centInput_pt1535` intercept `0.503146`, slope `0`; `ptFine_centInput` intercept `0.528476`, slope `-0.000767`; `ptFine_cent3` intercept `0.531206`, slope `-0.000960`; `ptFine_cent7` intercept `0.543716`, slope `-0.001513`. | This is still a hard stop before MC. Need visually inspect/approve target-80 threshold diagnostics and, ideally, compare actual per-bin target-efficiency behavior before generating/submitting any RecoilJets MC YAML. | Review `/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auauTightBDTValidation/model_validation_condor_20260511_194832/bdt_working_point_target80_diagnostics.png` in chat. If approved, generate frozen target-80 YAML, do one local runtime test, then submit paired `isSimEmbedded`/`isSimEmbeddedInclusive`. |
| 🟢 Implemented / Staged | Per-model target-80 AuAu BDT working-point contract | 2026-05-11 local implementation completed and SFTP upload succeeded for `scripts/validate_auau_tight_bdt_on_sim.py`, `scripts/make_auau_bdt_target_wp_config.py`, `scripts/auau_tight_bdt_pipeline.sh`, `scripts/submit_auau_bdt_etfine_centstudy_target80.sh`, `src_AuAu/RecoilJets_AuAu.cc`, `src_AuAu/RecoilJets_AuAu.h`, and `macros/Fun4All_recoilJets_unified_impl.C`. A second SFTP upload added staged multi-campaign helpers `scripts/submit_auau_bdt_targetwp_pair.sh` and `scripts/prepare_auau_bdt_target80_available_campaigns.sh`. Static checks passed locally: Python compile with `PYTHONPYCACHEPREFIX=/tmp/recoiljets_pycache`, `bash -n` for pipeline/target80 wrappers/uploader, and `git diff --check`. New validation outputs are `bdt_working_points_target80.json`, `.yaml`, `_runtime_fragment.yaml`, `.csv`, and `bdt_working_point_target80_diagnostics.png`. Runtime now supports per-tight-mode `auau_tight_bdt_working_point_entries` with linear or binned thresholds, falling back to the old global cut only when no matching entry is configured. The available-campaign prep script currently stages target-80 YAMLs for already-trained/validated registries: expanded 5-40, width-study 15-30, and fine-`E_T` 15-35, if their validation report paths are provided. | Needs SDCC rebuild because `src_AuAu/RecoilJets_AuAu.{cc,h}` and `macros/Fun4All_recoilJets_unified_impl.C` changed. Do not submit target-80 MC until the validation report is READY and the target-80 diagnostics/linear fits have been inspected. The exact full matrix with 3x3-only/base+3x3 centrality-binned and `E_T x centrality` variants for both 5-40 and 15-30 is not fully trained yet; the prep helper covers currently available products and skips missing validation reports. | On the new node: rebuild `src_AuAu`, run/finish `validateOnSimCondor` for the desired model dir(s), run `prepare_auau_bdt_target80_available_campaigns.sh` with validation report env vars, pull/inspect the target-80 WP diagnostics, do one local runtime test per generated YAML, then only submit paired MC with `submit_auau_bdt_targetwp_pair.sh` if the cuts look sane. |
| 🟢 Trained / Ready For Validation | Fine-`E_T` centrality BDT sidecar campaign for `15-35 GeV` | 2026-05-11 local implementation added campaign `etfine-centstudy` with 89 planned models (`centInput_pt1535`, `ptFine_centInput`, `ptFine_cent3`, `ptFine_cent7`), runtime modes `auauEtFineCentInputBDT`, `auauEtFineCent3BDT`, `auauEtFineCent7BDT`, YAML `macros/analysis_config_auau_bdt_etfine_centstudy_wp050.yaml`, and paired MC wrapper `scripts/submit_auau_bdt_etfine_centstudy_wp050.sh`. Visible terminal evidence on `sphnxuser06`: Justin ran `RJ_NOTIFY_EMAILS=just0131@gmail.com RJ_ML_PYTHON=/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python RJ_AUAU_BDT_TRAIN_PARALLEL=1 RJ_AUAU_BDT_XGB_N_JOBS=1 ./scripts/auau_tight_bdt_pipeline.sh trainEtFineCentStudyFromExtraction SOURCE=/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049 MODEL_DIR=/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_etfine_centstudy_current`. Training-tree validation passed with `entries=12000316 signal=8526032 background=3474284 files=8000`; training wrote cache `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_etfine_centstudy_current/training_matrix_etfine_centstudy.npz`, registry `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_etfine_centstudy_current/model_registry.json`, trained `reports=89 selected_specs=89`, and `expanded applyCheck opened 89 TMVA ROOT files`. Summary printed `RECOILJETS_AUAU_TIGHT_BDT_ETFINE_CENTSTUDY_TRAINING_V1 status=READY`, `pt_bins=15,17,19,21,23,25,27,30,35`, and model `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_etfine_centstudy_current/auau_tight_bdt_centInput_pt1535_tmva.root`. The later `unexpected EOF while looking for matching '"'` was caused by the pipeline script being overwritten during the running process after the READY summary; follow-up terminal check showed `bash -n scripts/auau_tight_bdt_pipeline.sh` clean, `registry OK`, and exactly `89` `*_tmva.root` files. | No blocker for validation. Do not run MC submission until `validateOnSimCondor` is READY, target-80 working-point diagnostics are pulled and inspected in chat, and a local runtime test passes. | Run registry validation on this model dir from `sphnxuser06` or another chosen node: `RJ_NOTIFY_EMAILS=just0131@gmail.com RJ_ML_PYTHON=/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python ./scripts/auau_tight_bdt_pipeline.sh validateOnSimCondor SOURCE=/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049 MODEL_DIR=/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_etfine_centstudy_current MODEL_REGISTRY=/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_etfine_centstudy_current/model_registry.json groupSize 100`. |
| 🟡 Running / Analysis-Only | Fixed width-window WP0.50 RecoilJets MC validation on `sphnxuser04` | Visible terminal evidence on 2026-05-11 from `/sphenix/u/patsfan753/scratch/thesisAnalysis` on `sphnxuser04`: after diagnosing the earlier DAGMan mass-removal, Codex uploaded fixed `RecoilJets_Condor_submit.sh` and `submit_auau_bdt_widthstudy_windows_wp050.sh`. Justin then ran `/tmp/resubmit_widthstudy_windows_wp050_from_scratch.sh`, typed `YES`, and removed all old `patsfan753` jobs on that submit host. The fixed campaign tag is `widthstudy_windows_wp050_fixed_20260511_180500`; submit log `/sphenix/u/patsfan753/scratch/thesisAnalysis/submit_widthstudy_windows_wp050_fixed_20260511_180500.log`; config dir `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_generated_configs/widthstudy_windows_wp050_fixed_20260511_180500`; model dir `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_centinput_widthstudy_windows_current`; working point `0.50`; windows `5-35`, `10-35`, `15-35`; `groupSize=7`; worker memory `10000MB`; memory retry cap `16000MB`; merge memory hint `8000MB`; `RJ_WIDTHSTUDY_AUTO_MERGE=0`; `RJ_HOLD_FAILED_WORKERS=1`. Worker clusters were submitted as one continuous span: `pt5to35` clusters `1140092`-`1140095`, `pt10to35` clusters `1140096`-`1140099`, `pt15to35` clusters `1140100`-`1140103`, each with `1429` jobs, in submit order signal Photon12/Photon20 then inclusive-background Jet12/Jet20. Bulk roots are `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded_widthstudy_windows_wp050_fixed_20260511_180500_pt5to35`, `/simembeddedinclusive_widthstudy_windows_wp050_fixed_20260511_180500_pt5to35`, and corresponding `_pt10to35` and `_pt15to35` roots. Merge/output roots are `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_widthstudy_windows_wp050_fixed_20260511_180500_pt5to35`, `_pt10to35`, and `_pt15to35`. The wrapper printed manual `mergeRecoilJets.sh firstRound/secondRound/finalStitch` commands for each dataset/window instead of launching parent auto-DAGs. Final queue snapshot after submission showed `11674 idle`, `5474 running`, and `0 held` for the user query; `7210 removed` are from the deliberate preflight cleanup of old jobs. Gmail check found two unread stale `auto_simembeddedinclusive_final_ready FAILED` emails from the removed pre-fix auto workflow; Codex marked them read after identifying them as not belonging to the new analysis-only clusters. | Analysis workers are running, but merge is intentionally not automatic. Need emails/queue evidence that worker clusters complete cleanly before running manual merges. If any worker fails now it should stay held for diagnosis because the fixed submit file sets worker failure hold behavior. | Watch Gmail and queue for clusters `1140092`-`1140103`. If no held jobs and workers drain cleanly, run the printed merge commands window-by-window. If held jobs appear, inspect with `condor_q -hold patsfan753 -af ClusterId ProcId HoldReason HoldReasonCode RequestMemory MemoryUsage NumJobStarts Args | head -80` before removing or releasing anything. |
| 🔴 Partial Failure / Diagnose | Width-window WP0.50 RecoilJets MC validation on `sphnxuser04` | Visible terminal evidence on 2026-05-11 from `/sphenix/u/patsfan753/scratch/thesisAnalysis` on `sphnxuser04`: Justin ran `RJ_NOTIFY_EMAILS=just0131@gmail.com RJ_WIDTHSTUDY_MODEL_DIR=/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_centinput_widthstudy_windows_current ./scripts/submit_auau_bdt_widthstudy_windows_wp050.sh`. Parseable header: `RECOILJETS_AUAU_BDT_WIDTHSTUDY_WINDOWS_SUBMIT_V1`, campaign `widthstudy_windows_wp050_20260511_172931`, config dir `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_generated_configs/widthstudy_windows_wp050_20260511_172931`, windows `5:35,10:35,15:35`, working point `0.50`, `group_size=7`, worker memory `6000MB`, merge memory `6000MB`, merge group size `75`, notify `just0131@gmail.com`. Submitted pairs: `pt5to35` signal DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembedded_20260511_172932/RecoilJets_auto_simembedded_20260511_172932.dag` cluster `1140074`, inclusive DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembeddedinclusive_20260511_172956/RecoilJets_auto_simembeddedinclusive_20260511_172956.dag` cluster `1140077`; `pt10to35` signal DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembedded_20260511_173021/RecoilJets_auto_simembedded_20260511_173021.dag` cluster `1140080`, inclusive DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembeddedinclusive_20260511_173045/RecoilJets_auto_simembeddedinclusive_20260511_173045.dag` cluster `1140083`; `pt15to35` signal DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembedded_20260511_173111/RecoilJets_auto_simembedded_20260511_173111.dag` cluster `1140086`, inclusive DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembeddedinclusive_20260511_173138/RecoilJets_auto_simembeddedinclusive_20260511_173138.dag` cluster `1140089`. Output roots are `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_widthstudy_windows_wp050_20260511_172931_pt5to35`, `_pt10to35`, and `_pt15to35`; bulk roots are `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded[_inclusive]_widthstudy_windows_wp050_20260511_172931_pt*to35`. Wrapper printed `TRACK_THIS_CAMPAIGN=widthstudy_windows_wp050_20260511_172931` and submit log `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_generated_configs/widthstudy_windows_wp050_20260511_172931/submit_widthstudy_windows_wp050_20260511_172931.log`. Post-submit queue snapshot showed `14290` jobs for query: `11806 idle`, `2470 running`, `14 held`. Later terminal showed `11159` jobs: `10001 running`, `0 held`, `1158 removed`. Gmail reported two consumed `FAILED` final emails for `isSimEmbedded` signal DAGs: `pt5to35` DAG `auto_workflow_simembedded_20260511_172932` and `pt10to35` DAG `auto_workflow_simembedded_20260511_173021`, both with `rescue_file_count=1`. | Signal `pt5to35` and `pt10to35` failed and some jobs were removed by DAG abort/rescue behavior. Need inspect DAGMan/node logs to identify whether this is memory, missing model/config, output collision, or fast worker failure. Other DAGs may still be running; do not pull/use outputs until final email evidence is settled. | Run one compact SDCC diagnostic for the two failed signal DAGs, then decide whether to patch/rerun only failed windows or stop all six and relaunch with a fixed wrapper/config. |
| 🟢 READY / Pull + Submit MC | 9-model width-window AuAu tight-BDT `validateOnSimCondor` on `sphnxuser04` | Visible terminal evidence on 2026-05-11 from `/sphenix/u/patsfan753/scratch/thesisAnalysis` on `sphnxuser04`: after the 9-model width-window training reached `READY`, Justin submitted `RJ_NOTIFY_EMAILS=just0131@gmail.com RJ_ML_PYTHON=/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python ./scripts/auau_tight_bdt_pipeline.sh validateOnSimCondor SOURCE=/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049 MODEL_DIR=/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_centinput_widthstudy_windows_current MODEL_REGISTRY=/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_centinput_widthstudy_windows_current/model_registry.json groupSize 100`. Path plan: report `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/model_validation_condor_20260511_171943`; submit dir `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/condor_sub/auauTightBDTValidate_20260511_171943`; DAG `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/condor_sub/auauTightBDTValidate_20260511_171943/auau_tight_bdt_validateOnSimCondor.dag`; score caches `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/model_validation_condor_20260511_171943/score_caches`; `root files=8000`, `shards=80`, `groupSize=100`, `scoreMaxPerShard=5000`, `request_memory=2500MB`; DAGMan cluster `1139989`. Queue later drained to zero in the visible terminal. Gmail on 2026-05-11 reported `[RecoilJets][auauTightBDT_validateOnSimCondor][READY]` from `sphnxuser04` for this report; Codex consumed the email and marked it read. | None for validation. Need local pull and then production-style paired MC outputs for `isSimEmbedded` and `isSimEmbeddedInclusive`. | Pull with `SSH_AUTH_SOCK="$(launchctl getenv SSH_AUTH_SOCK)" ./scripts/sftp_get_recoiljets_outputs.sh auauTightBDTValidation /sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/model_validation_condor_20260511_171943`, inspect rankings/AUCs, then submit paired MC with `RJ_NOTIFY_EMAILS=just0131@gmail.com RJ_WIDTHSTUDY_MODEL_DIR=/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_centinput_widthstudy_windows_current ./scripts/submit_auau_bdt_widthstudy_windows_wp050.sh`. |
| 🟢 Trained / Ready For Validation | 9-model width-window AuAu tight-BDT training on `sphnxuser04` | Terminal evidence on 2026-05-11 from `/sphenix/u/patsfan753/scratch/thesisAnalysis` on `sphnxuser04`: first parallel attempt with `RJ_AUAU_BDT_TRAIN_PARALLEL=3` passed 8000-file validation and wrote cache but was killed by memory pressure. Justin reran serially with `RJ_AUAU_BDT_TRAIN_PARALLEL=1 RJ_AUAU_BDT_XGB_N_JOBS=1` in the same `MODEL_DIR`. The serial run passed training-tree validation with `entries=12000316`, `signal=8526032`, `background=3474284`, `files=8000`, wrote `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_centinput_widthstudy_windows_current/model_registry.json`, trained `reports=9 selected_specs=9`, and `expanded applyCheck opened 9 TMVA ROOT files`. Final parseable summary: `RECOILJETS_AUAU_TIGHT_BDT_WIDTHSTUDY_WINDOWS_TRAINING_V1 status=READY`, source `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049`, model dir `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_centinput_widthstudy_windows_current`, registry `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_centinput_widthstudy_windows_current/model_registry.json`, pT windows `5:35,10:35,15:35`. Trained specs are base widths, 3x3 widths only, and base+3x3 widths for each pT window. | None for training. Validation has not yet been run for this 9-model window set. | Run registry validation: `RJ_NOTIFY_EMAILS=just0131@gmail.com RJ_ML_PYTHON=/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python ./scripts/auau_tight_bdt_pipeline.sh validateOnSimCondor SOURCE=/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049 MODEL_DIR=/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_centinput_widthstudy_windows_current MODEL_REGISTRY=/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_centinput_widthstudy_windows_current/model_registry.json groupSize 100`. If validation is READY, submit both SIM samples for all three windows at BDT `>0.50` with `RJ_NOTIFY_EMAILS=just0131@gmail.com RJ_WIDTHSTUDY_MODEL_DIR=/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_centinput_widthstudy_windows_current ./scripts/submit_auau_bdt_widthstudy_windows_wp050.sh`. |
| 🟢 Pulled / Ready For QA | WP0.50 width-study signal+background pair local | Gmail final READY email `19e1890f5e7c8728` for the resubmitted signal campaign `widthstudy_pt1530_wp050_resub_mem10_20260511_132314` reported DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembedded_20260511_132315/RecoilJets_auto_simembedded_20260511_132315.dag`, `rescue_file_count=0`, and final output `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_widthstudy_pt1530_wp050_resub_mem10_20260511_132314/simembedded`. FinalStitch READY email `19e1890bea7abc46` reported clusters `433787 433788 433789 433790`. On 2026-05-11 Codex pulled the four merged photon12+20 signal ROOT products into `/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auau_widthstudy_pt1530_wp050/combinedSimOnlyEMBEDDED`, the same local base already holding the four original clean WP0.50 embeddedJet12+20 inclusive-background files. Local check found 8 ROOT files total under that folder, size `1.3G`. Consumed READY/STARTED pipeline email IDs `19e1890f5e7c8728`, `19e1890bea7abc46`, `19e188cf79a28f0f`, and `19e188cc886d107d` were marked read. | Signal/background are paired locally, but remote provenance differs: signal came from the mem10 resubmission after the original signal failed at 6 GB memory; inclusive background came from the original clean WP0.50 run. | Use `/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auau_widthstudy_pt1530_wp050/combinedSimOnlyEMBEDDED` for WP0.50 QA and plots; keep the mixed remote provenance noted in plot notes if needed. |
| 🟡 Resubmitted / Watch Email | WP0.50 `isSimEmbedded` signal rerun with higher worker memory on `sphnxuser07` | Visible terminal evidence on 2026-05-11 13:23 EDT from `/sphenix/u/patsfan753/scratch/thesisAnalysis` on `sphnxuser07`: after the first paste attempt from `~` failed with `env: './RecoilJets_Condor_submit.sh': No such file or directory`, Justin ran the resubmit block from the repo checkout. Campaign tag `widthstudy_pt1530_wp050_resub_mem10_20260511_132314`; YAML `macros/analysis_config_auau_bdt_widthstudy_pt1530_wp050.yaml`; dataset `isSimEmbedded`; worker memory `RJ_REQUEST_MEMORY=10000MB`; memory retry cap `RJ_AUTO_MEMORY_RETRY_CAP_MB=16000`; firstRound merge memory `RJ_SIM_FIRSTROUND_REQUEST_MEMORY=8000MB`; `groupSize=7`; bulk base `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded_widthstudy_pt1530_wp050_resub_mem10_20260511_132314`; merge base `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_widthstudy_pt1530_wp050_resub_mem10_20260511_132314`; snapshot `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_snapshots/simembedded_20260511_132315`; DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembedded_20260511_132315/RecoilJets_auto_simembedded_20260511_132315.dag`; DAGMan cluster `433725`. Submitter printed two analysis nodes (`run28_embeddedPhoton12` and `run28_embeddedPhoton20`), `1429` jobs per sample, fanout `4 cfg tags from one DST pass`, `vz=[10]`, and `TRACK_THIS_CAMPAIGN=widthstudy_pt1530_wp050_resub_mem10_20260511_132314`. | This replaces the failed WP0.50 signal DAG `auto_workflow_simembedded_20260511_114620`, which failed because embeddedPhoton20 workers exceeded the 6 GB cgroup memory limit. Need READY/CHECK/FAILED email before pulling signal. | Watch Gmail for `auto_workflow_simembedded_20260511_132315` and campaign `widthstudy_pt1530_wp050_resub_mem10_20260511_132314`. When READY with `rescue_file_count=0`, pull with `SFTP_GET_CONFIG_YAML=macros/analysis_config_auau_bdt_widthstudy_pt1530_wp050.yaml SFTP_GET_REMOTE_DIR_OVERRIDE=/sphenix/u/patsfan753/scratch/thesisAnalysis/output_widthstudy_pt1530_wp050_resub_mem10_20260511_132314/simembedded SFTP_GET_LOCAL_COMBINED_BASE=/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auau_widthstudy_pt1530_wp050_resub_mem10/combinedSimOnlyEMBEDDED ./scripts/sftp_get_recoiljets_outputs.sh isSimEmbedded`. |
| 🟢 Pulled / Waiting For Pair | WP0.50 `isSimEmbeddedInclusive` completed inclusive-background half | Gmail final READY email `19e17ff2cf1bacc6` reported `isSimEmbeddedInclusive` output `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_widthstudy_pt1530_wp050_20260511_114619/simembeddedinclusive`, DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembeddedinclusive_20260511_114707/RecoilJets_auto_simembeddedinclusive_20260511_114707.dag`, `rescue_file_count=0`. On 2026-05-11 Codex pulled the four merged embeddedJet12+20 ROOT products with `SFTP_GET_CONFIG_YAML=macros/analysis_config_auau_bdt_widthstudy_pt1530_wp050.yaml`, `SFTP_GET_REMOTE_DIR_OVERRIDE=/sphenix/u/patsfan753/scratch/thesisAnalysis/output_widthstudy_pt1530_wp050_20260511_114619/simembeddedinclusive`, and `SFTP_GET_LOCAL_COMBINED_BASE=/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auau_widthstudy_pt1530_wp050/combinedSimOnlyEMBEDDED`; helper ended `[OK] Download complete. Downloaded 4 file(s).` | Signal half is not yet usable; original WP0.50 `isSimEmbedded` failed and is being rerun as `widthstudy_pt1530_wp050_resub_mem10_20260511_132314`. Keep the original inclusive-background provenance separate from the resubmitted signal provenance in plots/notes. | After WP0.50 signal rerun is READY and pulled, decide whether to pair this inclusive half with the resubmitted signal outputs for WP0.50 ROC/reco-efficiency/ID-efficiency plots. |
| 🟢 Pulled / Ready For QA | WP0.80 signal and inclusive-background width-study outputs | Gmail READY emails reported WP0.80 campaign `widthstudy_pt1530_wp080_20260511_111743` clean for both halves: `isSimEmbedded` DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembedded_20260511_111744/RecoilJets_auto_simembedded_20260511_111744.dag`, output `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_widthstudy_pt1530_wp080_20260511_111743/simembedded`; `isSimEmbeddedInclusive` DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembeddedinclusive_20260511_111810/RecoilJets_auto_simembeddedinclusive_20260511_111810.dag`, output `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_widthstudy_pt1530_wp080_20260511_111743/simembeddedinclusive`; both `rescue_file_count=0`. On 2026-05-11 Codex pulled four merged photon12+20 signal files and four merged embeddedJet12+20 inclusive-background files under `/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auau_widthstudy_pt1530_wp080/combinedSimOnlyEMBEDDED`; both helper calls ended `[OK] Download complete. Downloaded 4 file(s).` | Needs ROOT/plot QA before slide use. | Use the WP0.80 local folder for first reco-efficiency, ID-efficiency, isolation, and ROC-facing QA while WP0.50 signal rerun completes. |
| 🟢/🟡 Email checkpoint | 2026-05-11 `sphnxuser07` 15-30 GeV width-study RecoilJets MC validation | Gmail check on 2026-05-11 13:04 EDT consumed and marked read the latest `RecoilJets Pipeline` final emails. WP0.80 campaign `widthstudy_pt1530_wp080_20260511_111743` is READY for both halves: `isSimEmbedded` DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembedded_20260511_111744/RecoilJets_auto_simembedded_20260511_111744.dag`, `rescue_file_count=0`, final output `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_widthstudy_pt1530_wp080_20260511_111743/simembedded`; and `isSimEmbeddedInclusive` DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembeddedinclusive_20260511_111810/RecoilJets_auto_simembeddedinclusive_20260511_111810.dag`, `rescue_file_count=0`, final output `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_widthstudy_pt1530_wp080_20260511_111743/simembeddedinclusive`. WP0.50 campaign `widthstudy_pt1530_wp050_20260511_114619` has READY evidence for `isSimEmbeddedInclusive`: DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembeddedinclusive_20260511_114707/RecoilJets_auto_simembeddedinclusive_20260511_114707.dag`, `rescue_file_count=0`, final output `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_widthstudy_pt1530_wp050_20260511_114619/simembeddedinclusive`. | No matching final email was found yet for WP0.50 `isSimEmbedded` / signal DAG `auto_workflow_simembedded_20260511_114620` in the Gmail searches. It may still be running, delayed, or require a compact SDCC queue/DAG check. | Pull WP0.80 signal+inclusive when ready to QA. For WP0.50, either wait for the signal final email or inspect `auto_workflow_simembedded_20260511_114620` / cluster `433538` before pulling the pair. |
| 🟡 Submitted / Watch Email | 15-30 GeV width-study WP0.50 RecoilJets MC validation on `sphnxuser07` | 2026-05-11 local change created `macros/analysis_config_auau_bdt_widthstudy_pt1530_wp050.yaml`, identical to WP0.80 except `auau_tight_bdt_min_intercept: 0.50`, and uploaded it to SDCC with `./scripts/sftp_push_recoiljets.sh analysis_config_auau_bdt_widthstudy_pt1530_wp050.yaml`. Visible terminal evidence at 2026-05-11 11:46-11:48 EDT from `/sphenix/u/patsfan753/scratch/thesisAnalysis` on `sphnxuser07`: Justin ran `RJ_NOTIFY_EMAILS=just0131@gmail.com RJ_WIDTHSTUDY_CONFIG_YAML=macros/analysis_config_auau_bdt_widthstudy_pt1530_wp050.yaml RJ_WIDTHSTUDY_CAMPAIGN_TAG=widthstudy_pt1530_wp050_$(date +%Y%m%d_%H%M%S) ./scripts/submit_auau_bdt_widthstudy_pt1530_wp080.sh`. A harmless pasted-shell artifact printed `-bash: can: command not found`, then the actual wrapper completed. Campaign tag printed `widthstudy_pt1530_wp050_20260511_114619`. `isSimEmbedded`: bulk `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded_widthstudy_pt1530_wp050_20260511_114619`, merge `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_widthstudy_pt1530_wp050_20260511_114619`, snapshot `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_snapshots/simembedded_20260511_114620`, DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembedded_20260511_114620/RecoilJets_auto_simembedded_20260511_114620.dag`, DAGMan cluster `433538`, two analysis nodes (`run28_embeddedPhoton12`, `run28_embeddedPhoton20`), `1429` jobs per sample, fanout `4 cfg tags from one DST pass`. `isSimEmbeddedInclusive`: bulk `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembeddedinclusive_widthstudy_pt1530_wp050_20260511_114619`, same merge root `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_widthstudy_pt1530_wp050_20260511_114619`, snapshot `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_snapshots/simembeddedinclusive_20260511_114707`, DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembeddedinclusive_20260511_114707/RecoilJets_auto_simembeddedinclusive_20260511_114707.dag`, DAGMan cluster `433542`, two analysis nodes (`run28_embeddedJet12`, `run28_embeddedJet20`), `1429` jobs per sample, fanout `4 cfg tags from one DST pass`. Worker and firstRound merge memory were both `6000MB`; production `vz=[10]`; wrapper printed `TRACK_THIS_CAMPAIGN=widthstudy_pt1530_wp050_20260511_114619`. | Earlier WP0.80 campaign is also live/running on the same node. Keep output roots separate when pulling. Need final READY/CHECK/FAILED emails and queue state. | Watch Gmail for `widthstudy_pt1530_wp050_20260511_114619` final emails. If email is incomplete after the queue drains, inspect DAGs `auto_workflow_simembedded_20260511_114620` and `auto_workflow_simembeddedinclusive_20260511_114707`. |
| 🟡 Running / Watch Email | 15-30 GeV width-study WP0.80 RecoilJets MC validation on `sphnxuser07` | Visible terminal evidence on 2026-05-11 11:17 EDT from `/sphenix/u/patsfan753/scratch/thesisAnalysis` on `sphnxuser07`: after `src_AuAu` `make clean && makeProject` finished and installed `libRecoilJetsAuAu` into `/sphenix/u/patsfan753/thesisAnalysis_auau/install`, Justin ran `cd /sphenix/u/patsfan753/scratch/thesisAnalysis; RJ_NOTIFY_EMAILS=just0131@gmail.com ./scripts/submit_auau_bdt_widthstudy_pt1530_wp080.sh`. Wrapper campaign tag printed `widthstudy_pt1530_wp080_20260511_111743`. The visible first submission is `isSimEmbedded` with YAML `macros/analysis_config_auau_bdt_widthstudy_pt1530_wp080.yaml`, bulk base `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded_widthstudy_pt1530_wp080_20260511_111743`, merge base `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_widthstudy_pt1530_wp080_20260511_111743`, `groupSize=7`, worker memory `6000MB`, merge memory hint `6000MB`, production `vz=[10]`, snapshot `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_snapshots/simembedded_20260511_111744`, DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembedded_20260511_111744/RecoilJets_auto_simembedded_20260511_111744.dag`, and fanout `4 cfg tags from one DST pass` for the width-study matrix. Follow-up terminal queue evidence at 2026-05-11 ~11:31 EDT shows active Condor jobs on `sphnxuser07`: `Total for patsfan753: 5715 jobs; 0 idle, 5715 running, 0 held`, with visible workers from cluster `433535` running `RecoilJets_Condor_AuAu.sh run28_embeddedJet20`, consistent with the inclusive half also being submitted and running. The YAML name and config keys define strict WP0.80: `auau_tight_bdt_min_intercept: 0.80`, `auau_tight_bdt_apply_pt_min: 15.0`, `auau_tight_bdt_apply_pt_max: 30.0`. Gmail search at that time found no READY/CHECK/FAILED emails yet for this campaign. | Need final DAGMan cluster IDs and READY/CHECK/FAILED emails for both final auto workflows. Current live evidence is healthy: many running, zero held. | Watch Gmail for `widthstudy_pt1530_wp080_20260511_111743` final emails. If email is incomplete after the queue drains, inspect DAG dirs `auto_workflow_simembedded_20260511_111744` and the matching `auto_workflow_simembeddedinclusive_*` directory. |
| 🟢 READY / Submit MC | 15-30 GeV width-study BDT validation on `sphnxuser07` | 2026-05-11 terminal/email evidence from `sphnxuser07`: training summary `RECOILJETS_AUAU_TIGHT_BDT_WIDTHSTUDY_PT1530_TRAINING_V1 status=READY` wrote model dir `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_centinput_widthstudy_pt1530_current` with three TMVA models: base widths, 3x3 widths, and base+3x3 widths, all for `15 <= cluster_Et < 30` and centrality as input. Validation command submitted DAGMan cluster `433444`, report `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/model_validation_condor_20260511_110255`, then Gmail reported `[RecoilJets][auauTightBDT_validateOnSimCondor][READY]` with `total_entries=12000316`, `signal_entries=8526032`, `background_entries=3474284`, `scored_entries=400000`, finite score fraction `1`. AUCs: base widths `0.769382`, 3x3-only `0.768251`, base+3x3 `0.773786`. Codex pulled the full validation report locally to `dataOutput/auauTightBDTValidation/model_validation_condor_20260511_110255`; key files include `validation_curves.root`, `validation_auc_table.csv`, `validation_threshold_table.csv`, `validation_feature_summary.csv`, and PNG ROC/score/heatmap diagnostics. The consumed READY email was marked read. | Validation AUC is lower than the broad 5-40 campaign because this is a strict `15-30 GeV` focused study, but the apples-to-apples result favors base+3x3 widths. MC RecoilJets validation still needs to confirm the strict runtime pT guard and WP0.80 behavior in final histograms. | From `/sphenix/u/patsfan753/scratch/thesisAnalysis` on the chosen submit node, run `RJ_NOTIFY_EMAILS=just0131@gmail.com ./scripts/submit_auau_bdt_widthstudy_pt1530_wp080.sh`. Track the printed `TRACK_THIS_CAMPAIGN`, DAGs, and output roots. If `src_AuAu` was not rebuilt after the pT-guard/model-selector patch, first run `cd src_AuAu && make clean && makeProject`. |
| 🟢 Pulled / Ready For QA | `bdt3x3_20260510_223529_wp080` one-row MC outputs | 2026-05-11 Codex pulled the READY sphnxuser01 WP0.80 `AuAuCentInput3x3BDT` campaign outputs using `SFTP_GET_REMOTE_DIR_OVERRIDE` and `SFTP_GET_CFG_MATCH=AuAuCentInput3x3BDT`. Local signal file: `dataOutput/auau_bdt_mc_validation/bdt3x3_20260510_223529_wp080/preselectionNewPPG12_tightAuAuCentInput3x3BDT_nonTightAuAuBDTComplement_baseVariant/photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root` (`196M`). Local inclusive-background file: `dataOutput/auau_bdt_mc_validation/bdt3x3_20260510_223529_wp080/preselectionNewPPG12_tightAuAuCentInput3x3BDT_nonTightAuAuBDTComplement_baseVariant/embeddedJet12and20merged_SIM/RecoilJets_embeddedJet12plus20_MERGED.root` (`81M`). Remote roots were `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_bdt3x3_20260510_223529_wp080/simembedded` and `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_bdt3x3_20260510_223529_wp080/simembeddedinclusive`. | This pull contains only the one 3x3 centrality-input BDT cfg row at WP0.80, not the broader no-3x3 WP0.80 comparison set. | Use these files for focused WP0.80 3x3-vs-prior photon efficiency, ID-efficiency, isolation, purity/leakage, and xJ-sensitive QA. Keep plots labeled as `15?`/`5-40` according to the actual model used by this campaign. |
| 🟢 Pulled / ID Checked | AuAuCentInput3x3BDT WP0.50/WP0.80 one-row MC validation on `sphnxuser01` | Visible terminal evidence on 2026-05-10 22:35-22:37 EDT from `/sphenix/u/patsfan753/scratch/thesisAnalysis` on `sphnxuser01`: Justin submitted the separate 3x3-only campaign with `TRACK_THIS_CAMPAIGN=bdt3x3_20260510_223529`. All four commands used `RJ_NOTIFY_EMAILS=just0131@gmail.com`, `RJ_PROFILE_JOB=1`, `RJ_PHOTON_ID_ROW_MATCH=AuAuCentInput3x3BDT`, `RJ_ID_FANOUT_MAX_ROWS=1`, `RJ_REQUEST_MEMORY=6000MB`, `RJ_AUTO_MEMORY_RETRY_CAP_MB=10000`, `RJ_SIM_FIRSTROUND_REQUEST_MEMORY=6000MB`, `RJ_SIM_MERGE_GROUP_SIZE=75`, `groupSize=7`, and isolated bulk/merge roots. The submitter printed `RJ_PHOTON_ID_ROW_MATCH=AuAuCentInput3x3BDT -> 1/11 photon-ID row(s)` and `Fanout outputs: 1 cfg tags from one DST pass` for every sample. WP0.50 `isSimEmbedded`: YAML `macros/analysis_config_auau_bdt_validation.yaml`, bulk `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded_bdt3x3_20260510_223529_wp050`, merge `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_bdt3x3_20260510_223529_wp050`, DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembedded_20260510_223530/RecoilJets_auto_simembedded_20260510_223530.dag`, DAGMan cluster `2674754`. WP0.50 `isSimEmbeddedInclusive`: bulk `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembeddedinclusive_bdt3x3_20260510_223529_wp050`, merge `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_bdt3x3_20260510_223529_wp050`, DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembeddedinclusive_20260510_223553/RecoilJets_auto_simembeddedinclusive_20260510_223553.dag`, DAGMan cluster `2674757`. WP0.80 `isSimEmbedded`: YAML `macros/analysis_config_auau_bdt_validation_wp080.yaml`, bulk `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded_bdt3x3_20260510_223529_wp080`, merge `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_bdt3x3_20260510_223529_wp080`, DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembedded_20260510_223618/RecoilJets_auto_simembedded_20260510_223618.dag`, DAGMan cluster `2674760`. WP0.80 `isSimEmbeddedInclusive`: bulk `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembeddedinclusive_bdt3x3_20260510_223529_wp080`, merge `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_bdt3x3_20260510_223529_wp080`, DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembeddedinclusive_20260510_223645/RecoilJets_auto_simembeddedinclusive_20260510_223645.dag`, DAGMan cluster `2674763`. Each dataset has two analysis sample nodes (`Photon12/20` or `Jet12/20`) with about `1429` worker jobs per sample. Gmail READY evidence consumed on 2026-05-11: all four final auto-workflow emails reported `status=READY` with `rescue_file_count=0`: WP0.50 `isSimEmbedded` final output `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_bdt3x3_20260510_223529_wp050/simembedded`; WP0.50 `isSimEmbeddedInclusive` final output `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_bdt3x3_20260510_223529_wp050/simembeddedinclusive`; WP0.80 `isSimEmbedded` final output `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_bdt3x3_20260510_223529_wp080/simembedded`; WP0.80 `isSimEmbeddedInclusive` final output `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_bdt3x3_20260510_223529_wp080/simembeddedinclusive`. Consumed Gmail message IDs were marked read. Codex pulled the four one-row outputs locally under `dataOutput/auau_bdt3x3_mc_validation/wp050` and `dataOutput/auau_bdt3x3_mc_validation/wp080`, checked the merged ROOT files by reading the expected `SIM/h_Eiso_tight` and `SIM/h_Eiso_nonTight` histogram families, and wrote the first ID-efficiency comparison to `dataOutput/auau_bdt3x3_mc_validation/id_efficiency_compare/`. | No current blocker from email: all four isolated 3x3 one-row workflows are READY. Still keep these outputs separate from the larger no-3x3 WP0.80 work on `sphnxuser07`. | Continue QA from the pulled local outputs: compare 3x3 vs full-cluster centrality-input BDT and pT/centrality-binned BDTs in reco efficiency, isolation, ABCD behavior, shower-shape templates, leakage/fake rates, and xJ-sensitive inputs. |
| 🟢 Pulled / Plot QA Started | Clean WP0.80 no-3x3 fanout1 MC validation reset on `sphnxuser07` | Visible terminal evidence on 2026-05-10 22:27-22:43 EDT from `/sphenix/u/patsfan753/scratch/thesisAnalysis` on `sphnxuser07`: Justin ran the clean reset driver `/tmp/recoiljets_submit_wp080_no3x3_fanout1_clean_both.sh` after removing stale WP0.80 clusters `433083-433100`. The driver set `RJ_CONFIG_YAML=macros/analysis_config_auau_bdt_validation_wp080_no3x3.yaml`, `RJ_ID_FANOUT_MAX_ROWS=1`, `RJ_REQUEST_MEMORY=4000MB`, `RJ_SIM_FIRSTROUND_REQUEST_MEMORY=6000MB`, `RJ_PROFILE_JOB=1`, `groupSize=7`, `RJ_SIMEMBED_DEST_BASE=/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded_bdtwp080_no3x3_fanout1_clean`, `RJ_SIMEMBEDINCLUSIVE_DEST_BASE=/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembeddedinclusive_bdtwp080_no3x3_fanout1_clean`, and `RJ_MERGE_OUT_BASE_OVERRIDE=/sphenix/u/patsfan753/scratch/thesisAnalysis/output_bdtwp080_no3x3_fanout1_clean`. Clean DAGs: `isSimEmbedded` namespace `simembedded_20260510_222739`, DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembedded_20260510_222739/RecoilJets_auto_simembedded_20260510_222739.dag`, active DAGMan cluster `433101`; `isSimEmbeddedInclusive` namespace `simembeddedinclusive_20260510_223238`, DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembeddedinclusive_20260510_223238/RecoilJets_auto_simembeddedinclusive_20260510_223238.dag`, active DAGMan cluster `433105`. Crucial evidence: each cfg/sample node printed `Fanout outputs: 1 cfg tags from one DST pass`, matching the successful WP0.50 execution shape. A queue mapper showed older duplicate no-3x3 DAGs `433089` (`simembedded_20260510_221722`, output `output_bdtwp080_no3x3_fanout1`) and `433094` (`simembeddedinclusive_20260510_222122`, same old output root) were duplicate older attempts; Justin ran `condor_rm 433089 433094`, and the post-settle check showed `433089/433094` status `3` removed with child jobs status `X`. Remaining clean live work: `433101` children `433102/433103/433104` and `433105` children `433106/433107/433108/433109`, all `mem=4000`; no held jobs were visible. Gmail check at 22:29 found only stale failures from old bad attempts (`214818`, `214844`, `220121`) and no new clean-run failure. | The submitter still prints a stale 12GB floor warning even though this clean mode intentionally uses fanout1 with `4000MB`; patch later so the warning only applies to multi-row fanout. Gmail READY evidence was found for both clean final outputs; Codex pulled 20 merged ROOT files locally under `dataOutput/auau_bdt_mc_validation/wp080_no3x3_fanout1_clean/combinedSimOnlyEMBEDDED` and generated WP0.80 slide-31/32-style efficiency PNGs/CSVs under `dataOutput/auau_bdt_mc_validation/wp080_no3x3_fanout1_clean/photon_efficiency_qa`. | Use the pulled WP0.80 no-3x3 files for comparison plots. Current finding: literal `score > 0.80` is very tight and strongly suppresses BDT photon efficiency relative to box-cuts, so label these as strict-threshold diagnostics, not 80% signal-efficiency working points. |
| 🟡 Split State / Memory Cliff | WP0.80 no-3x3 `isSimEmbedded` and `isSimEmbeddedInclusive` MC validation on `sphnxuser07` | Visible terminal evidence on 2026-05-10 from `/sphenix/u/patsfan753/scratch/thesisAnalysis` on `sphnxuser07`: Justin ran `/tmp/recoiljets_submit_wp080_no3x3_both.sh`. Config is `macros/analysis_config_auau_bdt_validation_wp080_no3x3.yaml`, explicitly removing the not-yet-validated `AuAuCentInput3x3BDT` row. Worker memory was `RJ_REQUEST_MEMORY=9000MB`, retry cap `RJ_AUTO_MEMORY_RETRY_CAP_MB=12000`, merge memory `RJ_SIM_FIRSTROUND_REQUEST_MEMORY=6000MB`, `groupSize=7`. `isSimEmbedded` timestamp namespace `simembedded_20260510_214818`; auto DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembedded_20260510_214818/RecoilJets_auto_simembedded_20260510_214818.dag`; DAGMan cluster `433080`; output roots `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded_bdtwp080_no3x3` and `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_bdtwp080_no3x3`; terminal confirmed `Fanout outputs: 10 cfg tags from one DST pass` for both embedded photon samples. Gmail on 2026-05-10 21:52 reported `[RecoilJets][auto_simembedded_final_ready][FAILED]` for this no-3x3 signal DAG with `rescue_file_count=1`; email was consumed and marked read. Follow-up terminal diagnostics showed the actual root cause: many `433081.*` signal workers were held/evicted with Condor code `34/102`, `Job has gone over cgroup memory limit of 9088 megabytes`, last measured usage about `9063-9067 MB`. Codex patched and uploaded `RecoilJets_Condor_submit.sh` on 2026-05-10 so future `analysis_config_auau_bdt_validation*` embedded workers default to a `12000MB` floor and `RJ_AUTO_MEMORY_RETRY_CAP_MB=16000` when the caller does not explicitly override memory. `isSimEmbeddedInclusive` timestamp namespace `simembeddedinclusive_20260510_214844`; auto DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembeddedinclusive_20260510_214844/RecoilJets_auto_simembeddedinclusive_20260510_214844.dag`; DAGMan cluster visible as `433083`; worker cluster visible as `433085` for `run28_embeddedJet20`; terminal showed the inclusive queue progressing with `0` held and idle count reaching `0` in the DAG tail. | Signal-side no-3x3 needs more than 9GB; it failed at the cgroup boundary, not because of the removed 3x3 row. Inclusive side looked healthy in the latest terminal view but still needs final READY/CHECK evidence. | Let inclusive finish or fail independently. Resubmit only `isSimEmbedded` with the uploaded submitter and either explicit `RJ_REQUEST_MEMORY=12000MB` or no override; use a distinct output root if inclusive is still live. |
| 🟡 Submitting / Watch Holds | WP0.80 no-3x3 `isSimEmbedded` 12GB rerun on `sphnxuser07` | Visible terminal evidence on 2026-05-10 from `/sphenix/u/patsfan753/scratch/thesisAnalysis` on `sphnxuser07`: Justin reran only the failed signal-side campaign with `RJ_CONFIG_YAML=macros/analysis_config_auau_bdt_validation_wp080_no3x3.yaml`, `RJ_REQUEST_MEMORY=12000MB`, `RJ_AUTO_MEMORY_RETRY_CAP_MB=16000`, `RJ_SIM_FIRSTROUND_REQUEST_MEMORY=6000MB`, `RJ_SIMEMBED_DEST_BASE=/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded_bdtwp080_no3x3_12gb`, `RJ_MERGE_OUT_BASE_OVERRIDE=/sphenix/u/patsfan753/scratch/thesisAnalysis/output_bdtwp080_no3x3_12gb`, `./RecoilJets_Condor_submit.sh isSimEmbedded condorDoAll groupSize 7`. Submitter printed `SIM worker memory request: 12000MB (12000 MB)`, production `vz=[10]`, artifact namespace `simembedded_condorDoAll_20260510_220121`, snapshot dir `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_snapshots/simembedded_20260510_220121`, and auto DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembedded_20260510_220121/RecoilJets_auto_simembedded_20260510_220121.dag`. | Need final DAGMan/worker cluster IDs and first queue summary after submission finishes. | After submit output completes, record cluster IDs, run `condor_q`, and watch for any held jobs or final READY/CHECK email. |
| 🟢 Submitter Default Updated | AuAu BDT validation fanout safety default | 2026-05-10 Codex patched and uploaded `/sphenix/u/patsfan753/scratch/thesisAnalysis/RecoilJets_Condor_submit.sh`: for `isSimEmbedded`/`isSimEmbeddedInclusive condorDoAll` using any `analysis_config_auau_bdt_validation*` YAML, the submitter now defaults `RJ_ID_FANOUT_MAX_ROWS=1` when the user has not set it explicitly. This matches the successful WP0.50 campaign execution shape and prevents accidental 10-row TMVA fanout memory blowups. Syntax checks passed locally before upload. | None for the submitter default. Existing failed/removed WP0.80 jobs still need cleanup/resubmission. | Remove the active failed WP0.80 clusters, then rerun both samples with the fanout1 default or explicit `RJ_ID_FANOUT_MAX_ROWS=1`, using dedicated WP0.80 no3x3 output roots. |
| 🟡 Running / Watch Holds | WP0.80 `isSimEmbedded` full MC validation on `sphnxuser07` | Visible terminal evidence on 2026-05-10 from `/sphenix/u/patsfan753/scratch/thesisAnalysis` on `sphnxuser07`: after the memory patch/upload, Justin ran the one-block WP0.80 submit driver. It submitted `isSimEmbedded` with config `macros/analysis_config_auau_bdt_validation_wp080.yaml`, worker memory `RJ_REQUEST_MEMORY=9000MB`, retry cap `RJ_AUTO_MEMORY_RETRY_CAP_MB=12000`, merge memory `RJ_SIM_FIRSTROUND_REQUEST_MEMORY=6000MB`, output roots `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded_bdtwp080` and `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_bdtwp080`. Timestamp namespace `simembedded_20260510_213932`; auto DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembedded_20260510_213932/RecoilJets_auto_simembedded_20260510_213932.dag`; DAGMan cluster `433077`; analysis worker clusters visible as `433078`/`433079` with `2858` jobs total. Terminal at 21:40 showed `826 running`, `2032 idle`, `0 held`. | The paired `isSimEmbeddedInclusive` did not submit from the same driver because `RJ_CLEAN_OUTPUT_BASE=1` cleanup guard refused to run while the first DAG was already in the queue. This is a driver sequencing/cleanup issue, not a memory/physics failure. | Watch `433077`/worker clusters for held jobs and final email. Do not submit inclusive with cleanup while this DAG is live. After signal sample finishes or if we intentionally skip cleanup for an isolated inclusive root, submit `isSimEmbeddedInclusive` with the same 9000MB settings. |
| 🔴 Failed / Ready To Resubmit With Higher Memory | WP0.80 `isSimEmbedded` resubmit on `sphnxuser07` | Gmail on 2026-05-10 21:29 EDT reported `[RecoilJets][auto_simembedded_final_ready][FAILED]` for the new WP0.80 auto workflow. DAG path: `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembedded_20260510_212618/RecoilJets_auto_simembedded_20260510_212618.dag`; DAGMan out: `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembedded_20260510_212618/RecoilJets_auto_simembedded_20260510_212618.dag.dagman.out`; nodes log: `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembedded_20260510_212618/RecoilJets_auto_simembedded_20260510_212618.dag.nodes.log`; `rescue_file_count=1`; final output base `/sphenix/u/patsfan753/scratch/thesisAnalysis/output_bdtwp080/simembedded`. Visible terminal diagnostic showed worker cluster `433073` queued `2858` jobs, then many workers held/aborted. Worker log examples reported `Job has gone over cgroup memory limit of 4608 megabytes. Last measured usage: 4592 megabytes` with `RequestMemory=4500`, before useful event processing. This caused `DAG Limit: Max number of held jobs was reached`. On 2026-05-10 Codex patched and uploaded `RecoilJets_Condor_submit.sh` so embedded SIM with `analysis_config_auau_bdt_validation*` defaults to a `9000MB` worker floor and `RJ_AUTO_MEMORY_RETRY_CAP_MB=12000` when the user does not explicitly override memory. | Failed partial WP0.80 output roots from the 21:26 attempt may exist; clean them only through the guarded WP0.80 resubmission driver. | Resubmit both WP0.80 MC datasets with one driver block using explicit `RJ_REQUEST_MEMORY=9000MB RJ_AUTO_MEMORY_RETRY_CAP_MB=12000 RJ_CLEAN_OUTPUT_BASE=1`; watch terminal and Gmail for new DAG IDs and held-job state. |
| 🔴 Failed / Ready To Resubmit | WP0.80 `isSimEmbedded` MC validation on `sphnxuser07` | Visible terminal evidence on 2026-05-10 from `/sphenix/u/patsfan753/scratch/thesisAnalysis` on `sphnxuser07`: user submitted `RJ_NOTIFY_EMAILS=just0131@gmail.com RJ_CONFIG_YAML=macros/analysis_config_auau_bdt_validation_wp080.yaml RJ_SIMEMBED_DEST_BASE=/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded_bdtwp080 RJ_MERGE_OUT_BASE_OVERRIDE=/sphenix/u/patsfan753/scratch/thesisAnalysis/output_bdtwp080 RJ_SIM_FIRSTROUND_REQUEST_MEMORY=6000MB ./RecoilJets_Condor_submit.sh isSimEmbedded condorDoAll groupSize 7`. YAML/fanout namespace `simembedded_condorDoAll_20260510_211424`; DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembedded_20260510_211424/RecoilJets_auto_simembedded_20260510_211424.dag`; DAGMan cluster `433065`; two analysis nodes for `run28_embeddedPhoton12` and `run28_embeddedPhoton20`, `1429` jobs each, `11` cfg fanout outputs per DST pass. Terminal DAG tail showed `STATUS_ERROR`, `Node return val=-1002`, `ULOG_JOB_ABORTED` on worker cluster `433066`, then DAGMan aborted and removed remaining jobs. Root cause was operational: command set `RJ_SIM_FIRSTROUND_REQUEST_MEMORY=6000MB` for merge jobs but omitted analysis-worker `RJ_REQUEST_MEMORY`, so workers fell back to generic low SIM memory. Local submitter patched and uploaded on 2026-05-10 so embedded `condorDoAll` defaults to a 4000MB worker floor and prints worker vs merge memory separately. | Failed partial WP0.80 outputs may exist under the dedicated WP0.80 bulk/output roots; clean them on resubmit. | After confirming no old WP0.80 jobs remain, resubmit with explicit `RJ_REQUEST_MEMORY=4500MB RJ_CLEAN_OUTPUT_BASE=1` and the same WP0.80 YAML/output roots. Then track the new DAG cluster. |
| 🔴 Failed / Ready To Resubmit | WP0.80 `isSimEmbeddedInclusive` MC validation on `sphnxuser07` | Visible terminal evidence on 2026-05-10 from `/sphenix/u/patsfan753/scratch/thesisAnalysis` on `sphnxuser07`: user submitted `RJ_NOTIFY_EMAILS=just0131@gmail.com RJ_CONFIG_YAML=macros/analysis_config_auau_bdt_validation_wp080.yaml RJ_SIMEMBEDINCLUSIVE_DEST_BASE=/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembeddedinclusive_bdtwp080 RJ_MERGE_OUT_BASE_OVERRIDE=/sphenix/u/patsfan753/scratch/thesisAnalysis/output_bdtwp080 RJ_SIM_FIRSTROUND_REQUEST_MEMORY=6000MB ./RecoilJets_Condor_submit.sh isSimEmbeddedInclusive condorDoAll groupSize 7`. YAML/fanout namespace `simembeddedinclusive_condorDoAll_20260510_211550`; DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembeddedinclusive_20260510_211550/RecoilJets_auto_simembeddedinclusive_20260510_211550.dag`; DAGMan cluster `433068`; two analysis nodes for `run28_embeddedJet12` and `run28_embeddedJet20`, `1429` jobs each, `11` cfg fanout outputs per DST pass. Terminal DAG tail showed `STATUS_ERROR`, `Node return val=-1002`, `ULOG_JOB_ABORTED` on worker cluster `433070`, then DAGMan aborted and removed remaining jobs. Root cause was operational: command set `RJ_SIM_FIRSTROUND_REQUEST_MEMORY=6000MB` for merge jobs but omitted analysis-worker `RJ_REQUEST_MEMORY`, so workers fell back to generic low SIM memory. Local submitter patched and uploaded on 2026-05-10 so embedded `condorDoAll` defaults to a 4000MB worker floor and prints worker vs merge memory separately. | Failed partial WP0.80 outputs may exist under the dedicated WP0.80 bulk/output roots; clean them on resubmit. | After confirming no old WP0.80 jobs remain, resubmit with explicit `RJ_REQUEST_MEMORY=4500MB RJ_CLEAN_OUTPUT_BASE=1` and the same WP0.80 YAML/output roots. Then track the new DAG cluster. |
| 🟢 Trained / Ready For Runtime Test | 3x3 shower-width AuAu tight-BDT add-on on `sphnxuser01` | 2026-05-11 terminal and Gmail evidence show `trainCentInput3x3FromExtraction` completed from `/sphenix/u/patsfan753/scratch/thesisAnalysis` on `sphnxuser01` after `src_AuAu` rebuild. Training-tree validation passed with `entries=12000316 signal=8526032 background=3474284 files=8000`. The one-model campaign wrote `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_expanded_20260509_152604/training_matrix_3x3.npz`, `model_registry_3x3.json`, and TMVA file `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_expanded_20260509_152604/auau_tight_bdt_centAsFeat3x3_pt5to40_tmva.root`; applyCheck opened 1 TMVA ROOT file. Gmail subject `[RecoilJets][auauTightBDT_trainCentInput3x3FromExtraction][READY]` was consumed and marked read. | Needs one tiny RecoilJets runtime check that the new cfg row loads and scores inside `RecoilJets_AuAu`; after that, submit only this row for `isSimEmbedded` and `isSimEmbeddedInclusive`. | Run `RJ_CONFIG_YAML=macros/analysis_config_auau_bdt_validation.yaml RJ_PHOTON_ID_ROW_MATCH=AuAuCentInput3x3BDT ./scripts/RecoilJets_Condor_submit.sh isSimEmbedded local 200 NFILES=1 SAMPLE=run28_embeddedPhoton20 VERBOSE=1`; if clean, submit both full one-row MC validations with `RJ_PHOTON_ID_ROW_MATCH=AuAuCentInput3x3BDT`. |
| 🟢 READY / 3x3-Only Validated | 3x3 `validateOnSimCondor` registry targeting on `sphnxuser01` | The earlier 2026-05-10 run at report root `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/model_validation_condor_20260510_215850`, submit root `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/condor_sub/auauTightBDTValidate_20260510_215850`, DAGMan cluster `2674590`, completed READY but validated the full expanded registry because `MODEL_REGISTRY` was ignored. After local patches and upload, Justin reran from `/sphenix/u/patsfan753/scratch/thesisAnalysis` on `sphnxuser01` with `RJ_NOTIFY_EMAILS=just0131@gmail.com RJ_ML_PYTHON=/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python RJ_AUAU_TIGHT_BDT_VALIDATE_STAMP=cent3x3_$(date +%Y%m%d_%H%M%S) ./scripts/auau_tight_bdt_pipeline.sh validateOnSimCondor SOURCE=/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049 MODEL_DIR=/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_expanded_20260509_152604 MODEL_REGISTRY=/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_expanded_20260509_152604/model_registry_3x3.json groupSize 100`. Terminal explicitly printed `registry : .../model_registry_3x3.json` and `model registry: .../model_registry_3x3.json`. Fresh report root: `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049/reports/model_validation_condor_cent3x3_20260510_221419`; submit root: `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/condor_sub/auauTightBDTValidate_cent3x3_20260510_221419`; DAGMan cluster `2674672`; 8000 ROOT files split into 80 shards, scoreMaxPerShard 5000, request memory `2500MB`. Gmail on 2026-05-10 22:15 EDT reported `[RecoilJets][auauTightBDT_validateOnSimCondor][READY]` with `model_registry=/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_expanded_20260509_152604/model_registry_3x3.json`, `centAsFeat3x3_pt5to40_auc=0.883141`, finite score fraction `1`, signal score mean `0.611924`, background score mean `0.255937`, `score_eiso_pearson=-0.373558`, `scored_entries=400000`, and `notes=none`; email was consumed and marked read. | The 3x3 model is now ML-validation READY. It still needs isolated RecoilJets production-style validation through `isSimEmbedded` and `isSimEmbeddedInclusive` before being mixed back into larger comparison matrices. A separate fresh email also reported `[RecoilJets][auto_simembedded_final_ready][FAILED]` for the WP0.80 no-3x3 12GB signal rerun; that is separate from this 3x3 validation and remains a production-memory/debug item. | Pull the 3x3 validation report if plots/tables are needed. Next production-style step: run only the `AuAuCentInput3x3BDT` row for `isSimEmbedded` and `isSimEmbeddedInclusive` using isolated output roots and `RJ_PHOTON_ID_ROW_MATCH=AuAuCentInput3x3BDT`. |
| 🟢 Pulled / Initial ROOT Sanity Passed | Full `isSimEmbedded` MC-first BDT validation production on `sphnxuser07` | User submitted on 2026-05-09 from `/sphenix/u/patsfan753/scratch/thesisAnalysis` with `RJ_NOTIFY_EMAILS=just0131@gmail.com RJ_PROFILE_JOB=1 RJ_ID_FANOUT_MAX_ROWS=1 RJ_REQUEST_MEMORY=4000MB RJ_SIM_FIRSTROUND_REQUEST_MEMORY=6000MB RJ_CONFIG_YAML=macros/analysis_config_auau_bdt_validation.yaml ./RecoilJets_Condor_submit.sh isSimEmbedded condorDoAll groupSize 7`. Timestamp namespace `simembedded_20260509_234130`; auto DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembedded_20260509_234130/RecoilJets_auto_simembedded_20260509_234130.dag`; DAGMan cluster `432408`. Gmail on 2026-05-10 reported old auto merge `READY`, but terminal inspection showed only default-YAML tags were promoted. Codex patched/uploaded merge scripts. User reran manual firstRound from `sphnxuser07` with explicit `MERGE_CONFIG_YAML`, then secondRound/finalStitch completed. Local pull on 2026-05-10 used `SFTP_GET_CONFIG_YAML=macros/analysis_config_auau_bdt_validation.yaml ./scripts/sftp_get_recoiljets_outputs.sh isSimEmbedded`, preserved old colliding files under `previous/20260510_190502_isSimEmbedded`, and downloaded the 10 campaign files. Local ROOT inspection found the 8 expanded AuAu BDT cfg products with photonJet12+20 final ROOT files, each paired with embedded-inclusive outputs after the separate pull. Representative BDT file contains `h_xJ` (5104), recoil response (768), truth/miss (20/16), ABCD Eiso (704), isolation/tight counters, and tight-isolated photon candidates. | Initial numeric scan shows global BDT families are promising; centrality-dependent production application has very low truth-reco efficiency in this first pass and needs diagnostic plots/inspection before being promoted. | Generate the MC comparison suite: truth-photon reco efficiency, embedded-jet leakage, isolation/ABCD behavior, tight/complement shower shapes, recoil/xJ stability, and centrality/photon-E_T family comparisons. |
| 🟢 Pulled / Initial ROOT Sanity Passed | Full `isSimEmbeddedInclusive` MC-first BDT validation production on `sphnxuser07` | User submitted on 2026-05-09 from `/sphenix/u/patsfan753/scratch/thesisAnalysis` with `RJ_NOTIFY_EMAILS=just0131@gmail.com RJ_PROFILE_JOB=1 RJ_ID_FANOUT_MAX_ROWS=1 RJ_REQUEST_MEMORY=4000MB RJ_SIM_FIRSTROUND_REQUEST_MEMORY=6000MB RJ_CONFIG_YAML=macros/analysis_config_auau_bdt_validation.yaml ./RecoilJets_Condor_submit.sh isSimEmbeddedInclusive condorDoAll groupSize 7`. Timestamp namespace `simembeddedinclusive_20260509_234610`; auto DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembeddedinclusive_20260509_234610/RecoilJets_auto_simembeddedinclusive_20260509_234610.dag`; DAGMan cluster `432413`. Gmail on 2026-05-10 reported old auto merge `[RecoilJets][auto_simembeddedinclusive_final_ready][FAILED]`; diagnosis was default-YAML merge discovery despite complete BDT worker outputs. Codex patched/uploaded merge scripts. User reran manual firstRound from `sphnxuser07` with explicit `MERGE_CONFIG_YAML`, then secondRound/finalStitch completed. Local pull on 2026-05-10 used `SFTP_GET_CONFIG_YAML=macros/analysis_config_auau_bdt_validation.yaml ./scripts/sftp_get_recoiljets_outputs.sh isSimEmbeddedInclusive`, preserved old colliding files under `previous/20260510_191148_isSimEmbeddedInclusive`, and downloaded the 10 campaign files. Representative embeddedJet12+20 BDT ROOT contains the same key analysis families as the signal file: `h_xJ` (5104), recoil response, truth/miss, ABCD Eiso, isolation/tight counters, and tight-isolated photon candidates. | First-pass background leakage is nonzero as expected and must be compared with signal efficiency and physics weights before calling any working point final. | Use these files alongside `isSimEmbedded` to rank the BDT families by signal efficiency vs embedded-inclusive fake/leakage, then select 2-3 finalists for constrained production/data validation. |
| 🟢 Local Gate Passed / Full Submitted | MC-first RecoilJets validation of expanded AuAu tight-BDT families | Local/SDCC code path is prepared for the campaign-specific config `macros/analysis_config_auau_bdt_validation.yaml`. It compares 2 baselines plus 8 expanded AuAu BDT families using `newPPG12` preselection, shared BDT working point `score > 0.50`, and BDT-complement non-tight. Model directory is `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_expanded_20260509_152604`. Visible terminal evidence on 2026-05-09 from `sphnxuser07` shows `src_AuAu` rebuild completed successfully after the compile fix, then `isSimEmbedded local 100 SAMPLE=run28_embeddedPhoton20 VERBOSE=0` completed all 10 cfg rows. Representative final rows `tightAuAuPtCent3BDT` and `tightAuAuPtCent7BDT` reported `exit_code=0`, `elapsed_seconds=33`, and normal ROOT outputs around `5.7 MB`; earlier baseline/no-cent/cent-input/minority rows also reported `exit_code=0` with normal MB-scale outputs. Justin intentionally skipped smoke/local inclusive-background testing and submitted both full MC campaigns for overnight running; those live runs are tracked in the two rows above. | Residual risk is the intentionally skipped smoke/background-local gate. The live full campaigns will settle this through the normal email/DAG stages. | Tomorrow, do not rerun local/smoke by default. Check `RecoilJets Pipeline` Gmail first for clusters `432408` and `432413`; if READY, pull outputs and build the comparison suite. If CHECK/FAILED, inspect the matching DAG logs from the paths above. |
| ⚪ Removed / Paused | Full `isAuAu` data production on `sphnxuser04` | User reported on 2026-05-08 that all AuAu jobs were already removed after the strategic pivot. The paused campaign was DAGMan cluster `1134652`, timestamp `auau_condorAll_20260508_111239`, command `RJ_NOTIFY_EMAILS=just0131@gmail.com RJ_PROFILE_JOB=1 RJ_REQUEST_MEMORY=4500MB RJ_CLEAN_OUTPUT_BASE=1 ./RecoilJets_Condor_submit.sh isAuAu condor all groupSize 4`. | Do not restart AuAu data production until the simulation/BDT JSTG package is ready and explicitly re-approved. | Treat AuAu data as paused. If needed later, verify no stale matching jobs remain with `condor_q -nobatch -wide | grep -E '1134652|auau_condorAll_20260508_111239|auto_workflow_auau_20260508_111239|RecoilJets_Condor_AuAu.sh' || true`. |
| 🟡 JSTG Priority | Simulation + AuAu tight-BDT presentation package | User explicitly shifted focus on 2026-05-08 toward fast simulation iteration: cross-section estimation with `scripts/estimateEmbeddedPhotonXsec.sh`, reproducing `embeddedPhoton_stitchedTruthFilterPtSpectrum.png` in the new reference/reference/reference compact pipeline for both `isSimEmbedded` and `isSimEmbeddedInclusive`, richer BDT training QA, simulation validation, and a Monday draft for a 15-minute JSTG presentation next Wednesday. A local compact `isSimEmbedded` reference plot was produced with `macros/MakeEmbeddedPhotonXsecNormQA.C("preselectionReference_tightReference_nonTightReference_baseVariant")`; the truth-filter stitching plot itself is valid, while side iso-efficiency warnings are expected until that offline macro learns internal-iso suffixes. | Need canonical final-stitch outputs, cross-section/weight validation, legacy-vs-compact histogram parity, and BDT plots. Avoid slow local serial stitching and avoid new full AuAu data submissions. | Finish/pull stitched `isSimEmbedded` and `isSimEmbeddedInclusive` products; compare old scalar histograms to new internal-view suffixed histograms with matched cuts, e.g. old `h_Eiso_pT_10_12_cent_0_20` vs new `h_Eiso_isoR40_fixedIso4GeV_pT_10_12_cent_0_20`, guarding especially against `vz60` vs `vz10` mismatches; run the cross-section estimator; reproduce the reference stitched truth-filter pT spectrum; build BDT ROC/score/efficiency plots binned in centrality and pT; assemble draft slides by Monday. |
| 🟢 Aggregated | Embedded-inclusive Jet12/Jet20 xsec estimator | Visible terminal evidence on 2026-05-08 from `sphnxuser08` shows the patched `scripts/estimateEmbeddedPhotonXsec.sh` was run from `/sphenix/u/patsfan753/scratch/thesisAnalysis` with `./scripts/estimateEmbeddedPhotonXsec.sh firstPass --family inclusive --xsec-shards 50 --xsec-events 1000000 --mode interpreted`. Workdir: `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_snapshots/pythia_xsec_20260508_154701`; manifest: `/sphenix/u/patsfan753/scratch/thesisAnalysis/pythia_xsec_firstPass_20260508_154701.txt`; submit file: `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/pythia_xsec_20260508_154701.sub`; Condor cluster `631685`; `100` jobs total, `50` shards each for `EmbeddedJet12` and `EmbeddedJet20`, `1,000,000` raw events per shard, `900MB` request memory. The secondPass aggregation at 2026-05-08 16:04 found `Expected=100` and `Found=100` shard CSVs. Results: `EmbeddedJet12 n_ok=49,999,989 n_pass=482,371 eff_weight=9.64742212e-03 sigma_gen_mb=1.26139880e-01 sigma_eff_pb=1.21692467e+06 rel_stat=0.001440`; `EmbeddedJet20 n_ok=49,999,987 n_pass=532,317 eff_weight=1.06463428e-02 sigma_gen_mb=5.22431703e-03 sigma_eff_pb=5.56198698e+04 rel_stat=0.001371`. Combined CSV: `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_snapshots/pythia_xsec_20260508_154701/combined_results.csv`; summary: `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_snapshots/pythia_xsec_20260508_154701/combined_summary.txt`. | Numbers are generator-derived and high-stat, but should be sanity-checked before replacing constants. Compared with old PPG12-style constants, derived Jet12 is about `18.3%` lower and derived Jet20 about `11.2%` lower. The estimator follows Brian `embed_2025` production mapping (`Jet12 -> phpythia8_10GeV_JS_MDC2.cfg + 12 <= pT < 20`, `Jet20 -> phpythia8_20GeV_JS_MDC2.cfg + pT >= 20`). | After approval, update embedded-inclusive constants in `macros/AnalyzeRecoilJets.h`, `scripts/mergeRecoilJets.sh`, and `macros/MakeEmbeddedPhotonXsecNormQA.C`; then rerun/rebuild only the affected merge/QA products as needed. |
| 🟢 READY | AuAu tight-BDT `validateOnSimCondor` on `sphnxuser07` | Gmail on 2026-05-08 15:47 reported `[RecoilJets][auauTightBDT_validateOnSimCondor][READY]` from `sphnxuser07`; email was consumed and marked read. Source `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_015049`; model dir `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_20260508_104530`; report dir `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_015049/reports/model_validation_condor_20260508_154538`. Metrics: `total_entries=12427173`, `signal_entries=8751401`, `background_entries=3675772`, `scored_entries=400000`, finite score fraction `1`, AUCs `centINDcontrol=0.881426`, `centAsFeat=0.883994`, `centDepBDTs=0.88403`, score/eiso Pearson about `-0.380` to `-0.390`, notes `none`. Rich outputs exist remotely: `validation_metrics.json`, `validation_deep_diagnostics.json`, `validation_curves.root`, `validation_auc_table.csv`, `validation_threshold_table.csv`, `validation_feature_summary.csv`. | Need pull/inspect the richer diagnostics and generate final slide/analysis-note-quality validation plots. This validates the Condor-sharded validation machinery, not yet the constrained data+MC BDT-variant production pass. | Pull or copy the report artifacts from `model_validation_condor_20260508_154538`; then make centrality/pT ROC, AUC, threshold, and feature-summary plots for slide/update use; next physics step is constrained data+MC BDT-variant validation. |
| 🟢 READY Email | `isAuAu` auto workflow from `sphnxuser01` | Gmail heartbeat check on 2026-05-08 consumed `[RecoilJets][auto_auau_final_ready][READY]`, `[RecoilJets][data_final_auau_all][READY]`, and `[RecoilJets][data_sliceRuns_auau_all][READY]` from `sphnxuser01`. Reported DAG: `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_auau_20260507_221758/RecoilJets_auto_auau_20260507_221758.dag`; `rescue_file_count=0`; final stage clusters included `2658331`-`2658346`; next action from email is `./scripts/sftp_get_recoiljets_outputs.sh isAuAu`. | Need pull and ROOT-level old-vs-new validation before treating local `InputFiles/auau25` as replaced. Separate later `sphnxuser01` 3500MB attempt (`auto_workflow_auau_20260507_232024`, clusters `2659455/2659456`) previously showed a memory hold and should not be confused with this READY workflow. | Preserve current `InputFiles/auau25` into a cross-check folder, then pull with `./scripts/sftp_get_recoiljets_outputs.sh isAuAu` and compare representative legacy/root catalog content before full analysis use. |
| 🔴 Held Smoke | `isSimEmbedded` direct-fanout smoke | User submitted from `sphnxuser03` with `RJ_NOTIFY_EMAILS=just0131@gmail.com RJ_PROFILE_JOB=1 RJ_SMOKE_SIM_NEVENTS=10000 RJ_SMOKE_SIM_SAMPLE_LIMIT=2 RJ_SMOKE_REQUEST_MEMORY_EMBEDDED=3000MB ./RecoilJets_Condor_submit.sh isSimEmbedded condorDoAllSmoke groupSize 4`. DAG cluster `3016477`; worker clusters `3016478`/`3016479`. `condor_q ... -hold -af` showed cgroup memory holds: request `3000`, limit `3072 MB`, last measured usage about `3059-3062 MB`, image size `3250000`. | This was a memory-limit hold before useful profiling/merge evidence; do not tune full production from this pass. | Remove `3016477 3016478 3016479`, then rerun a smaller fanout smoke using `RJ_ID_FANOUT_MAX_ROWS=5` and/or higher memory to learn whether module fanout dominates RSS. |
| 🔴 Held Smoke | `isSimEmbedded` reduced-fanout smoke | User reran from `sphnxuser03` with reduced ID fanout (`RJ_ID_FANOUT_MAX_ROWS=5`) and `RJ_SMOKE_REQUEST_MEMORY_EMBEDDED=3000MB`. Terminal showed `72` jobs all held across worker clusters `3016488`-`3016493`; Gmail final email `[RecoilJets][auto_simembedded_firstRound_ready][FAILED]` for DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembedded_20260507_233440/RecoilJets_auto_simembedded_20260507_233440.dag` reported a rescue file. | Reducing fanout to 5 did not solve 3 GB memory holds, so baseline embedded/AuAu worker footprint is already near/above the 3 GB cgroup limit. | Remove `3016487 3016488 3016489 3016490 3016491 3016492 3016493`; next test should be a tiny single-sample memory ladder with higher request memory, not another broad 72-job smoke. |
| 🟢 Smoke Validated | `isSimEmbedded` single-fanout smoke | User submitted from `sphnxuser03` with `RJ_ID_FANOUT_MAX_ROWS=1`, `RJ_SMOKE_SIM_NEVENTS=10000`, `RJ_SMOKE_SIM_SAMPLE_LIMIT=1`, `RJ_SMOKE_SIM_MAX_JOBS_PER_SAMPLE_EMBEDDED=1`, `RJ_SMOKE_REQUEST_MEMORY_EMBEDDED=6000MB`, `groupSize=4`, `SAMPLE=run28_embeddedPhoton12`. Terminal showed the analysis jobs ran instead of holding, then entered SIM firstRound merge. Gmail reported `[RecoilJets][sim_firstRound_simembedded_all][READY]` and `[RecoilJets][auto_simembedded_firstRound_ready][READY]` for DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembedded_20260507_235449/RecoilJets_auto_simembedded_20260507_235449.dag`, `rescue_file_count=0`, final output base `/sphenix/u/patsfan753/scratch/thesisAnalysis/outputSmoke/simembedded_smokeTest_20260507_235449/simembedded`. | Smoke only covered one sample and one job per fanout shard; need profile reports before choosing full production settings. | Pull smoke reports with `./scripts/sftp_get_recoiljets_outputs.sh smokeTestLatest isSimEmbedded`, inspect max RSS/elapsed/output sanity, then decide whether full embedded uses `RJ_ID_FANOUT_MAX_ROWS=1` at `6000MB` or tests `fanout=2`. |
| 🟢 Validated | SIM auto-DAG secondRound integration | 2026-05-08 `isSimEmbedded` smoke from `sphnxuser03` reached `[RecoilJets][sim_firstRound_simembedded_all][READY]`, then `[RecoilJets][sim_secondRound_simembedded_all][READY]`, then `[RecoilJets][auto_simembedded_final_ready][READY]`. Final output base reported: `/sphenix/u/patsfan753/scratch/thesisAnalysis/outputSmoke/simembedded_smokeTest_20260508_005005/simembedded`; DAG: `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembedded_20260508_005006/RecoilJets_auto_simembedded_20260508_005006.dag`; rescue count `0`. | Need pull/inspect final smoke ROOTs and profile rows before full embedded production settings are final. | Pull with `./scripts/sftp_get_recoiljets_outputs.sh smokeFinal isSimEmbedded simembedded_smokeTest_20260508_005005`, inspect output count/catalog, then tune full `isSimEmbedded`. |
| 🟢 Smoke Validated | Direct-fanout AuAu smoke path | 2026-05-07 smoke tests produced 15 final ROOT files under pulled `InputFiles/pipelineSmoke/isAuAu/...`; catalog checks found expected iso views, jet pT/dphi internal keys, and representative `h_Eiso`, `h_ss`, ABCD/isolation, xJ, and jet QA families. Final notify/post-script email chain was fixed and validated. | Smoke cannot prove every rare histogram branch has statistics, but production mechanics and representative catalog were coherent. | Use full production evidence for final validation. |
| ⚪ Superseded | Two-run full-event `isAuAu` smoke on `sphnxuser07` | Earlier 2026-05-08 terminal evidence showed `sphnxuser07` worker cluster `431679` running 222-223 AuAu jobs with `0` held and memory sizes mostly below the later 4500 MB request. User then chose to remove it and submit full production instead. | Smoke outputs are not intended to be preserved as final validation unless later evidence says they completed before removal. | Do not chase this smoke unless Gmail/terminal later shows it completed and produced useful outputs. Current AuAu focus is the `sphnxuser06` full production row. |
| 🟡 Priority | `isSimEmbedded` | Needed for AuAu photon-ID ML signal. Smoke is now submitted from `sphnxuser03`. | Must validate nonempty outputs/reports before full production. | Use active smoke row above to tune memory/groupSize, then scale if clean. |
| 🔴 Merge Held | Full `isSimEmbedded` production on `sphnxuser03` | User submitted after the successful full-chain smoke on 2026-05-08 with `RJ_NOTIFY_EMAILS=just0131@gmail.com RJ_PROFILE_JOB=1 RJ_ID_FANOUT_MAX_ROWS=1 RJ_REQUEST_MEMORY=4000MB RJ_CLEAN_OUTPUT_BASE=1 ./RecoilJets_Condor_submit.sh isSimEmbedded condorDoAll groupSize 7`. Visible terminal confirmed `simembedded_condorDoAll_20260508_012157`, `1429` jobs per cfg/sample, both embeddedPhoton12/20, and later showed the workflow reached SIM firstRound merging. Terminal on 2026-05-08 then showed `150 held` `hadd_condor.sh` jobs across clusters around `3016723`-`3016755`, each with `SIZE 4151.0`, plus running merge DAGMan parents; Gmail had no fresh failure email yet. Local fix uploaded `mergeRecoilJets.sh` + `RecoilJets_Condor_submit.sh` so future embedded firstRound defaults to `groupSize=75`, supports `RJ_SIM_FIRSTROUND_REQUEST_MEMORY`, and emits general Condor memory-hold auto-release logic capped by `RJ_AUTO_MEMORY_RETRY_CAP_MB` default `8000` with two releases max. | This is a merge-stage memory/chunking problem, not a reason to redo DST analysis. Need either bump/release held firstRound hadd jobs or remove only the merge layer and rerun firstRound with smaller groups from existing raw ROOT outputs. | First inspect hold reason. If memory, try `condor_qedit` held `hadd_condor.sh` jobs to `RequestMemory 6000` and `condor_release`; if they hold again, remove only merge DAG/hadd clusters and rerun `RJ_SIM_FIRSTROUND_REQUEST_MEMORY=6000MB ./scripts/mergeRecoilJets.sh isSimEmbedded firstRound groupSize 75`, then `./scripts/mergeRecoilJets.sh isSimEmbedded secondRound condor`. Preserve old local `InputFiles/simEmbedded` before any final pull. |
| 🔴 Superseded / Failed | Old full `isSimEmbeddedInclusive` production on `sphnxuser01` | Terminal evidence on 2026-05-08 showed full `isSimEmbeddedInclusive` submission from `sphnxuser01`, timestamp namespace `simembeddedinclusive_condorDoAll_20260508_112828`, separated output base `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembeddedinclusive`, DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembeddedinclusive_20260508_112828/RecoilJets_auto_simembeddedinclusive_20260508_112828.dag`, top-level DAGMan cluster `2670772`, samples `run28_embeddedJet12`/`run28_embeddedJet20`, and `groupSize=7`. On 2026-05-08 the code was found stale for this campaign: `src_AuAu/RecoilJets_AuAu.cc` still used old PPG12-style embedded-inclusive ownership (`Jet12: 14<=pT<21`, `Jet20: 21<=pT<32`), and finalStitch weights still used `1.4903e6`/`6.2623e4` pb instead of the derived production-matched values. Gmail later reported `[RecoilJets][auto_simembeddedinclusive_final_ready][FAILED]` for the same `112828` DAG with `rescue_file_count=1`; that email is stale relative to the fixed `163858` resubmission. | Existing outputs from this campaign must not be used for final embedded-inclusive background normalization or stitched QA. Full worker rerun is required because the event rejection happened inside the RecoilJets worker before merge. | Ignore this campaign except for cleanup/debug archaeology. Use the active fixed-stitching `163858` campaign below. |
| 🟢 Pulled / Plot Ready | Fixed full `isSimEmbeddedInclusive` production on `sphnxuser08` | Visible terminal evidence on 2026-05-08 shows `src_AuAu` rebuild completed on `sphnxuser08` with `makeProject finished` into `/sphenix/u/patsfan753/thesisAnalysis_auau/install`. Justin then submitted from `/sphenix/u/patsfan753/scratch/thesisAnalysis`: `RJ_NOTIFY_EMAILS=just0131@gmail.com RJ_PROFILE_JOB=1 RJ_ID_FANOUT_MAX_ROWS=1 RJ_REQUEST_MEMORY=4000MB RJ_SIM_FIRSTROUND_REQUEST_MEMORY=6000MB RJ_CLEAN_OUTPUT_BASE=1 ./RecoilJets_Condor_submit.sh isSimEmbeddedInclusive condorDoAll groupSize 7`. Timestamp namespace: `simembeddedinclusive_condorDoAll_20260508_163858`; YAML/fanout dir: `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_yaml_overrides/simembeddedinclusive_condorDoAll_20260508_163858`; snapshot dir: `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_snapshots/simembeddedinclusive_20260508_163858`; auto DAG: `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembeddedinclusive_20260508_163858/RecoilJets_auto_simembeddedinclusive_20260508_163858.dag`; bulk dest base cleaned/reused: `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembeddedinclusive`; production `vz` selected `[10]`. Terminal showed `run28_embeddedJet12` has `9998` matched rows and `run28_embeddedJet20` has `10000` matched rows, with `1429` analysis jobs per sample/cfg at `groupSize=7`. Heartbeat check at 2026-05-08 17:31 EDT found active Condor worker cluster `631694` with `4714` jobs total, `3117` running, `1597` idle, and `0` held. Gmail on 2026-05-08 22:20 EDT reported `[RecoilJets][auto_simembeddedinclusive_final_ready][READY]` from `sphnxuser08` for DAG `/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/auto_workflow_simembeddedinclusive_20260508_163858/RecoilJets_auto_simembeddedinclusive_20260508_163858.dag`, `rescue_file_count=0`, final output base `/sphenix/u/patsfan753/scratch/thesisAnalysis/output/simembeddedinclusive`, next action `./scripts/sftp_get_recoiljets_outputs.sh isSimEmbeddedInclusive`; email was consumed and marked read. On 2026-05-08 Codex pulled 15 final merged inclusive-embedded SIM files with `./scripts/sftp_get_recoiljets_outputs.sh isSimEmbeddedInclusive`. The target file exists locally at `dataOutput/combinedSimOnlyEMBEDDED/preselectionReference_tightReference_nonTightReference_baseVariant/embeddedJet12and20merged_SIM/RecoilJets_embeddedJet12plus20_MERGED.root` (`143M`). ROOT inspection of `SIM/MERGE_INFO` showed `Nslices=2 shards=8 [embeddedJet12 Nraw=7116378 sigma_pb=1216924.67 w=23.3580] [embeddedJet20 Nraw=7597337 sigma_pb=55619.8698 w=1]`; `SIM/h_embedInclusiveStitch_filterJetPt_kept` exists and was plotted to `embeddedInclusiveJet_finalMergedStitchedTruthJetPtSpectrum.png`. Codex then pulled the raw per-sample secondRound files with `SFTP_GET_SIM_RAW=1 ./scripts/sftp_get_recoiljets_outputs.sh isSimEmbeddedInclusive`; the reference raw files exist at `InputFiles/InclusiveJetSIM_EMBEDDED/RecoilJets_embeddedJet12_ALL_preselectionReference_tightReference_nonTightReference_baseVariant.root` (`81M`) and `InputFiles/InclusiveJetSIM_EMBEDDED/RecoilJets_embeddedJet20_ALL_preselectionReference_tightReference_nonTightReference_baseVariant.root` (`120M`). ROOT inspection found `Jet12 kept=7.11638e6 raw=7.11638e6`, `Jet20 kept=7.5973e6 raw=7.59734e6`. The split component plot was generated at `dataOutput/combinedSimOnlyEMBEDDED/preselectionReference_tightReference_nonTightReference_baseVariant/embeddedJet12and20merged_SIM/embeddedInclusiveJet_stitchedTruthFilterPtSpectrum.png`. | The final merged file cannot be unmixed into components after the fact; the blue/red component plot comes from the per-sample secondRound inputs weighted offline with the same derived finalStitch constants. The helper's automatic all-config local merge was intentionally stopped after raw files and plot were produced because it was extra and slow. | Show/approve the split Jet12/Jet20 PNG in chat, then replace the plot on slide 5 of `JSTG_5_13_26`. The `recoiljets-pipeline-watchdog` remains paused because its off-condition is met for this campaign. |
| 🟢 Hardened Locally | `isSimEmbedded` vs `isSimEmbeddedInclusive` workflow separation | 2026-05-08 local code update separates `isSimEmbeddedInclusive` bulk and merge input base to `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembeddedinclusive`, keeps `isSimEmbedded` under `/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded`, and stops local smoke cleanup from scanning/removing `tmp_recoil_merge_*` directories that may belong to active merge DAGs on another submit host. | Needs upload to SDCC before rerunning inclusive local/smoke tests. | Upload `RecoilJets_Condor_submit.sh`, `mergeRecoilJets.sh`, and `recoiljets_cleanup.sh`; then rerun the `isSimEmbeddedInclusive` local test from the SDCC checkout. |
| 🟡 Priority | `isSim` production readiness | User wants the next testing pass after AuAu to validate SIM direct-fanout, memory, merge, and output availability before full submission. | Needs same quick chain validation as AuAu: analysis -> SIM firstRound -> READY/CHECK email. | After AuAu submit is stable, run `isSim` smoke/checkjobs and tune memory/grouping. |
| 🟡 Priority | `isSimEmbeddedInclusive` | Jet12/Jet20 lists reported complete with 10,000 matched rows each. | Needs cross-section constants, RecoilJets outputs, local pull, stitched validation. | Run Jet12/Jet20 xsec estimator, propagate constants, produce/pull/validate. |
| 🟡 Planned | z-vertex and centrality reweighting | User explicitly marked as top priority after next embedded passes. | Needs fresh embedded signal/background outputs. | Derive weights after next passes; train reweighted and non-reweighted BDT families. |
| 🟡 QA Needed | Region-C xJ purity correction | Code distinguishes `h_isIsolated_*` candidate counters from `h_xJpurityLead_*` event-leading counters. | Need confirm plotted/quoted purity and `N_gamma` normalization use intended pipeline. | Inspect RooUnfold purity-corrected mode and generated QA outputs. |
| 🔴 Stale For Affected Families | Non-reference ABCD/purity outputs produced before `507eabd3` | `codex_notes/DATASET_STATUS.md` records source fix after current affected products. | ABCD/purity-derived families may be invalid until SDCC transfer/rerun/pull. | Transfer `RecoilJets.cc RecoilJets_AuAu.cc`, rerun affected jobs, pull outputs, update evidence. |
| 🟢 Valid For Current Slide Task | AuAu scaled-trigger study | User reported plots are on slide 9; cluster `1127133`; 6404 jobs over 620 runs. | Presentation explanation/polish may still be needed. | Ask before archiving current slide task. |
| 🟢 Ready For Full Extraction | AuAu tight-BDT sidecar + expanded campaign infrastructure | Validated smoke from `sphnxuser02`: `RJ_AUAU_TIGHT_BDT_NEVENTS=5000 ./scripts/auau_tight_bdt_pipeline.sh smokeTest groupSize 5 maxJobs 6`, DAGMan cluster `2028148`, 24 extraction jobs, run root `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_012029`. DAG completed cleanly. Direct validation passed with `entries=37320 signal=26316 background=11004 files=24`; training wrote `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_20260508_014325`, all 5 TMVA files had `tmva=ok`, and `applyCheck` opened all 5 files. The false `CHECK` email/summary was traced to scheduler Python missing `libpython3.13.so.1.0`; `scripts/auau_tight_bdt_pipeline.sh` was patched locally and uploaded on 2026-05-08 to export the ML venv lib path in notify jobs and capture Python startup errors. On 2026-05-08, Justin ran the final local infrastructure test from `sphnxuser02` after the extract-only embedded ownership fix and expanded-campaign code: `RJ_ML_PYTHON=/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python RJ_AUAU_TIGHT_BDT_LOCAL_NEVENTS=1000 ./scripts/auau_tight_bdt_pipeline.sh localTest 1000 NFILES=1`, followed by `RJ_AUAU_BDT_CAMPAIGN_MAX_SPECS=6 RJ_AUAU_BDT_TRAIN_PARALLEL=2 RJ_ML_PYTHON=/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python ./scripts/auau_tight_bdt_pipeline.sh trainExpandedFromExtraction SOURCE=local_bdt_training_outputs/tight_20260508_231909`. Terminal evidence: local source `local_bdt_training_outputs/tight_20260508_231909`, validation `PASS` with `entries=1174 signal=825 background=349 files=4`; normal three-model training wrote `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_20260508_232602` with all TMVA exports `ok`; expanded subset training wrote `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_expanded_20260508_232623`, `model_registry.json`, trained `6` selected specs, and `expanded applyCheck opened 6 TMVA ROOT files`; final status `RECOILJETS_AUAU_TIGHT_BDT_EXPANDED_TRAINING_V1 status=READY`. The trailing `-bash: : command not found` was from a malformed shell diagnostic after the pipeline completed, not from the BDT workflow. | The old full extraction/model pair `auauTightBDT_20260508_015049` / `tight_20260508_104530` predates the extract-only embedded ownership fix and remains workflow-validation only for final model-choice evidence. The expanded campaign has been locally validated but not yet submitted at full scale. | Submit the fresh post-fix extraction from `sphnxuser02` or the chosen submit node with `RJ_NOTIFY_EMAILS=just0131@gmail.com RJ_ML_PYTHON=/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python RJ_AUAU_TIGHT_BDT_REQUEST_MEMORY=3000MB ./scripts/auau_tight_bdt_pipeline.sh condorExtract groupSize 5`. After the READY email gives a new `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_<timestamp>` source, run `RJ_NOTIFY_EMAILS=just0131@gmail.com RJ_ML_PYTHON=/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python ./scripts/auau_tight_bdt_pipeline.sh trainExpandedFromExtractionCondor SOURCE=<new post-fix source> groupSize 4`. |
| 🟢 Ready For Expanded Validation | AuAu tight-BDT expanded model campaign | Visible terminal evidence on 2026-05-09 from `sphnxuser02`: the expanded campaign DAGMan cluster `2028345` finished with cache + all 40 training shards complete. Submit root `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/condor_sub/auauTightBDTExpanded_20260509_152604`; model dir `/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_expanded_20260509_152604`; source `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049`. Finalized summary now reports `status=READY_WITH_SKIPS`, `expected_model_count=158`, `trained_model_count=147`, `skipped_model_count=11`, `missing_model_count=0`, `bad_model_count=0`. The skipped models are exactly the dedicated `35-40 GeV` pT-binned products, each with zero selected signal/background rows. On 2026-05-09 local code was patched/uploaded so future default expanded campaigns use pT-binned edges only through `26-35 GeV`, while broad `5-40 GeV` products remain. The validator now ignores registry reports with `status=skipped`. | The 147 trained models are usable for model-study validation. Do not promote dedicated `35-40 GeV` pT-bin models; high-pT tail is covered by broad `5-40 GeV` models. | Run registry-driven `validateOnSimCondor` on the expanded model dir, pull the resulting report with `./scripts/sftp_get_recoiljets_outputs.sh auauTightBDTValidation <remote-report-dir>`, then make ranked pT/centrality plots from the 147 trained models. |

## Status Legend

- 🟢 Valid / usable for stated purpose.
- 🟡 Needs evidence, QA, decision, or follow-up.
- 🔴 Stale, blocked, failed, or unsafe for stated use.

## Update Rules

- Do not mark an item green without dated evidence.
- If an item is done, move the corresponding task to
  `agent_context/TASK_BOARD.md` `Done Pending Removal` and ask Justin before
  archiving.
- Keep next commands exact and paste-ready when possible.
- The `RecoilJets pipeline watchdog` should be active only for rows that
  represent submitted/running jobs, ready outputs, failures, or explicit user
  monitoring requests. Add the off-condition when activating it, and pause it
  when that condition is met.
