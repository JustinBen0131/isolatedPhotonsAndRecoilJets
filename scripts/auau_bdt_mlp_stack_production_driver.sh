#!/usr/bin/env bash
set -euo pipefail

RJ_REPO_BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
readonly RJ_REPO_BASE

ML_PYTHON="${RJ_ML_PYTHON:-/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python}"
MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/matplotlib-${USER:-patsfan753}}"
SOURCE="${SOURCE:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049}"
STAMP="${STACK_STAMP:-$(date +%Y%m%d_%H%M%S)}"
RUN_ROOT="${RUN_ROOT:-/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/bdt_mlp_stack_production_${STAMP}}"
TEMPLATE="${TEMPLATE:-macros/analysis_config_auau_bdt_mlp_stack_template.yaml}"
CONFIG_OUT="${CONFIG_OUT:-macros/analysis_config_auau_bdt_mlp_stack_wp080.yaml}"
MLP_CACHE="${MLP_CACHE:-${SOURCE}/reports/mlp_model_validation_condor_stack_full_features_20260512_2338/score_caches.list}"
BDT_CACHE="${BDT_CACHE:-${SOURCE}/reports/model_validation_condor_stack_full_features_20260512_2338/score_caches.list}"
BDT_MODE="${BDT_MODE:-auauEtFineCent7BDT}"
MLP_MODE="${MLP_MODE:-auauCentInputBase3x3MLP}"
VARIANT="${VARIANT:-ptFine5to35_cent7_full}"
ALGORITHM="${ALGORITHM:-gbm}"
PT_EDGES="${PT_EDGES:-15,17,19,21,23,25,27,30,35}"
CENT_EDGES="${CENT_EDGES:-0,10,20,30,40,50,60,80}"
TARGET="${TARGET:-0.80}"
FIT_SPLIT="${FIT_SPLIT:-trainval}"
WP_SPLIT="${WP_SPLIT:-test}"
MAX_ROWS="${MAX_ROWS:-0}"
MAX_SHARDS="${MAX_SHARDS:-0}"
GBM_ESTIMATORS="${GBM_ESTIMATORS:-90}"
GBM_LEARNING_RATE="${GBM_LEARNING_RATE:-0.045}"
GBM_MAX_DEPTH="${GBM_MAX_DEPTH:-3}"
GBM_MAX_LEAF_NODES="${GBM_MAX_LEAF_NODES:-8}"
LINEAR_BACKEND="${LINEAR_BACKEND:-numpy}"
MAX_LINEAR_STEPS="${MAX_LINEAR_STEPS:-1800}"
NN_HIDDEN="${NN_HIDDEN:-64,32}"
NN_EPOCHS="${NN_EPOCHS:-180}"
NN_PATIENCE="${NN_PATIENCE:-24}"
NN_BATCH_SIZE="${NN_BATCH_SIZE:-8192}"
NN_LEARNING_RATE="${NN_LEARNING_RATE:-0.0015}"
NN_L2="${NN_L2:-0.001}"
REQUEST_MEMORY="${REQUEST_MEMORY:-4000MB}"
SIM_FIRSTROUND_REQUEST_MEMORY="${SIM_FIRSTROUND_REQUEST_MEMORY:-6000MB}"
GROUP_SIZE="${GROUP_SIZE:-7}"
OUTPUT_BASE="${OUTPUT_BASE:-/sphenix/u/patsfan753/scratch/thesisAnalysis/output_bdt_mlp_stack_wp080_${STAMP}}"
DEST_BASE="${DEST_BASE:-/sphenix/tg/tg01/bulk/jbennett/thesisAna/bdt_mlp_stack_wp080_${STAMP}}"

say() { printf '\033[1;36m[auauBDTMLPStack]\033[0m %s\n' "$*"; }
die() { printf '\033[1;31m[auauBDTMLPStack][ERR]\033[0m %s\n' "$*" >&2; exit 2; }
repo_path_or_absolute() {
  case "$1" in
    /*) printf '%s\n' "$1" ;;
    *) printf '%s/%s\n' "$RJ_REPO_BASE" "$1" ;;
  esac
}
mkdir -p "$MPLCONFIGDIR"
export MPLCONFIGDIR
CONFIG_PATH="$(repo_path_or_absolute "$CONFIG_OUT")"

usage() {
  cat <<'EOF'
Usage:
  ./scripts/auau_bdt_mlp_stack_production_driver.sh preflight
  ./scripts/auau_bdt_mlp_stack_production_driver.sh fullCache
  ./scripts/auau_bdt_mlp_stack_production_driver.sh train
  ./scripts/auau_bdt_mlp_stack_production_driver.sh trainingCurves
  ./scripts/auau_bdt_mlp_stack_production_driver.sh deriveWP80
  ./scripts/auau_bdt_mlp_stack_production_driver.sh makeConfig
  ./scripts/auau_bdt_mlp_stack_production_driver.sh parity
  ./scripts/auau_bdt_mlp_stack_production_driver.sh smoke
  ./scripts/auau_bdt_mlp_stack_production_driver.sh status
  RJ_DO_RUN=1 RJ_ALLOW_SUBMIT=1 RJ_WP_DIAGNOSTICS_APPROVED=1 ./scripts/auau_bdt_mlp_stack_production_driver.sh submitPair
  ./scripts/auau_bdt_mlp_stack_production_driver.sh pullQA

This driver promotes exactly one isolation-input-free BDT+MLP stacker.  It does
not submit production unless both RJ_DO_RUN=1 and RJ_ALLOW_SUBMIT=1 are set.
Inspect RUN_ROOT/diagnostics before submitPair.
EOF
}

require_wp_outputs() {
  local missing=0
  local path
  for path in \
    "${RUN_ROOT}/promotion_summary.json" \
    "${RUN_ROOT}/stack_working_points_target80.json" \
    "${RUN_ROOT}/stack_working_points_target80_cells.csv" \
    "${RUN_ROOT}/stack_training_history.csv" \
    "${RUN_ROOT}/diagnostics/stack_wp80_threshold_grid.png" \
    "${RUN_ROOT}/diagnostics/stack_wp80_signal_efficiency_grid.png" \
    "${RUN_ROOT}/diagnostics/stack_wp80_fake_rate_grid.png" \
    "${RUN_ROOT}/diagnostics/stack_wp80_threshold_vs_et_by_centrality.png"; do
    if [[ ! -s "$path" ]]; then
      say "MISSING required WP output: $path"
      missing=1
    fi
  done
  (( missing == 0 )) || die "WP80 diagnostics are incomplete; do not generate/submit production config yet"
}

mode="${1:-status}"

print_plan() {
  say "RECOILJETS_AUAU_BDT_MLP_STACK_PRODUCTION_DRIVER_V1"
  say "host=$(hostname -f 2>/dev/null || hostname)"
  say "repo=${RJ_REPO_BASE}"
  say "stamp=${STAMP}"
  say "run_root=${RUN_ROOT}"
  say "source=${SOURCE}"
  say "mlp_cache=${MLP_CACHE}"
  say "bdt_cache=${BDT_CACHE}"
  say "variant=${VARIANT} algorithm=${ALGORITHM}"
  say "bdt_mode=${BDT_MODE} mlp_mode=${MLP_MODE}"
  say "pt_edges=${PT_EDGES} cent_edges=${CENT_EDGES} target=${TARGET}"
  say "fit_split=${FIT_SPLIT} wp_split=${WP_SPLIT} max_rows=${MAX_ROWS} max_shards=${MAX_SHARDS}"
  say "gbm_estimators=${GBM_ESTIMATORS} gbm_learning_rate=${GBM_LEARNING_RATE} gbm_max_depth=${GBM_MAX_DEPTH} gbm_max_leaf_nodes=${GBM_MAX_LEAF_NODES}"
  say "linear_backend=${LINEAR_BACKEND} max_linear_steps=${MAX_LINEAR_STEPS}"
  say "nn_hidden=${NN_HIDDEN} nn_epochs=${NN_EPOCHS} nn_patience=${NN_PATIENCE} nn_batch_size=${NN_BATCH_SIZE} nn_lr=${NN_LEARNING_RATE} nn_l2=${NN_L2}"
  say "template=${TEMPLATE}"
  say "config_out=${CONFIG_OUT}"
  say "output_base=${OUTPUT_BASE}"
  say "dest_base=${DEST_BASE}"
}

case "$mode" in
  preflight)
    print_plan
    [[ -x "$ML_PYTHON" ]] || die "ML_PYTHON is not executable: $ML_PYTHON"
    [[ -s "$MLP_CACHE" ]] || die "Missing MLP_CACHE=$MLP_CACHE"
    [[ -s "$BDT_CACHE" ]] || die "Missing BDT_CACHE=$BDT_CACHE"
    "$ML_PYTHON" -m py_compile \
      "${RJ_REPO_BASE}/scripts/promote_auau_stacked_bdt_mlp.py" \
      "${RJ_REPO_BASE}/scripts/train_auau_stacked_bdt_mlp_sweep.py" \
      "${RJ_REPO_BASE}/scripts/make_auau_stacked_training_curves.py"
    say "preflight OK"
    ;;
  train)
    print_plan
    mkdir -p "$RUN_ROOT"
    "$ML_PYTHON" "${RJ_REPO_BASE}/scripts/promote_auau_stacked_bdt_mlp.py" train-wp \
      --mlp-cache "$MLP_CACHE" \
      --bdt-cache "$BDT_CACHE" \
      --outdir "$RUN_ROOT" \
      --variant "$VARIANT" \
      --algorithm "$ALGORITHM" \
      --target-signal-efficiency "$TARGET" \
      --fit-split "$FIT_SPLIT" \
      --wp-split "$WP_SPLIT" \
      --wp-pt-edges "$PT_EDGES" \
      --wp-cent-edges "$CENT_EDGES" \
      --max-rows "$MAX_ROWS" \
      --max-shards "$MAX_SHARDS" \
      --linear-backend "$LINEAR_BACKEND" \
      --max-linear-steps "$MAX_LINEAR_STEPS" \
      --gbm-estimators "$GBM_ESTIMATORS" \
      --gbm-learning-rate "$GBM_LEARNING_RATE" \
      --gbm-max-depth "$GBM_MAX_DEPTH" \
      --gbm-max-leaf-nodes "$GBM_MAX_LEAF_NODES" \
      --nn-hidden "$NN_HIDDEN" \
      --nn-epochs "$NN_EPOCHS" \
      --nn-patience "$NN_PATIENCE" \
      --nn-batch-size "$NN_BATCH_SIZE" \
      --nn-learning-rate "$NN_LEARNING_RATE" \
      --nn-l2 "$NN_L2"
    say "inspect diagnostics before makeConfig/submitPair: ${RUN_ROOT}/diagnostics"
    ;;
  trainingCurves|training-curves)
    print_plan
    [[ -s "${RUN_ROOT}/stack_training_history.csv" ]] || die "Missing ${RUN_ROOT}/stack_training_history.csv"
    "$ML_PYTHON" "${RJ_REPO_BASE}/scripts/make_auau_stacked_training_curves.py" \
      --history-csv "${RUN_ROOT}/stack_training_history.csv" \
      --outdir "${RUN_ROOT}/diagnostics/training_curves" \
      --models "${VARIANT}_${ALGORITHM}" \
      --top-n 1
    ;;
  fullCache)
    print_plan
    [[ -s "$MLP_CACHE" ]] || die "Missing MLP_CACHE=$MLP_CACHE"
    [[ -s "$BDT_CACHE" ]] || die "Missing BDT_CACHE=$BDT_CACHE"
    say "using existing aligned score-cache manifests"
    wc -l "$MLP_CACHE" "$BDT_CACHE"
    ;;
  deriveWP80)
    print_plan
    "$0" train
    ;;
  makeConfig)
    print_plan
    require_wp_outputs
    [[ -s "${RUN_ROOT}/stack_working_points_target80.json" ]] || die "Missing ${RUN_ROOT}/stack_working_points_target80.json"
    "$ML_PYTHON" "${RJ_REPO_BASE}/scripts/promote_auau_stacked_bdt_mlp.py" make-config \
      --template "${RJ_REPO_BASE}/${TEMPLATE}" \
      --working-points "${RUN_ROOT}/stack_working_points_target80.json" \
      --out "$CONFIG_PATH" \
      --bdt-mode "$BDT_MODE" \
      --mlp-mode "$MLP_MODE"
    say "config ready: ${CONFIG_OUT}"
    ;;
  parity)
    print_plan
    [[ -s "${RUN_ROOT}/artifacts/${VARIANT}_${ALGORITHM}.json" ]] || die "Missing stack artifact; run train first"
    say "parity gate is C++ runtime smoke on SDCC after upload/build"
    say "required acceptance: auau_tight_bdt_mlp_score finite and Python/C++ agreement within 1e-5 on representative rows"
    ;;
  smoke)
    print_plan
    [[ "${RJ_DO_RUN:-0}" == "1" ]] || die "Refusing smoke without RJ_DO_RUN=1"
    die "Smoke execution is intentionally site-specific; run after SDCC upload/build with the generated config and one representative input list."
    ;;
  status)
    print_plan
    for path in \
      "${RUN_ROOT}/promotion_summary.json" \
      "${RUN_ROOT}/stack_working_points_target80.json" \
      "${RUN_ROOT}/stack_working_points_target80.yaml" \
      "${RUN_ROOT}/stack_working_points_target80_cells.csv" \
      "${RUN_ROOT}/stack_training_history.csv" \
      "${RUN_ROOT}/diagnostics/stack_wp80_threshold_grid.png" \
      "${RUN_ROOT}/diagnostics/training_curves/${VARIANT}_${ALGORITHM}_training_curves.png" \
      "${RUN_ROOT}/diagnostics/stack_wp80_signal_efficiency_grid.png" \
      "${RUN_ROOT}/diagnostics/stack_wp80_fake_rate_grid.png" \
      "${RUN_ROOT}/diagnostics/stack_wp80_threshold_vs_et_by_centrality.png" \
      "$CONFIG_PATH"; do
      if [[ -s "$path" ]]; then say "FOUND $path"; else say "MISSING $path"; fi
    done
    ;;
  submitPair)
    print_plan
    [[ "${RJ_DO_RUN:-0}" == "1" ]] || die "Refusing submitPair without RJ_DO_RUN=1"
    [[ "${RJ_ALLOW_SUBMIT:-0}" == "1" ]] || die "Refusing submitPair without RJ_ALLOW_SUBMIT=1"
    [[ "${RJ_WP_DIAGNOSTICS_APPROVED:-0}" == "1" ]] || die "Refusing submitPair until WP plots/fits are inspected and RJ_WP_DIAGNOSTICS_APPROVED=1 is set"
    require_wp_outputs
    [[ -s "$CONFIG_PATH" ]] || die "Missing config ${CONFIG_OUT}"
    grep -q "auau_tight_bdt_mlp_stack_working_point_entries" "$CONFIG_PATH" || die "Config lacks stack working-point entries"
    log="/tmp/submit_bdt_mlp_stack_${STAMP}.log"
    {
      print_plan
      cd "$RJ_REPO_BASE"
      export RJ_CONFIG_YAML="$CONFIG_OUT"
      export RJ_NOTIFY_EMAILS="${RJ_NOTIFY_EMAILS:-just0131@gmail.com}"
      export RJ_PROFILE_JOB=1
      export RJ_ID_FANOUT_MAX_ROWS=1
      export RJ_REQUEST_MEMORY="$REQUEST_MEMORY"
      export RJ_SIM_FIRSTROUND_REQUEST_MEMORY="$SIM_FIRSTROUND_REQUEST_MEMORY"
      export RJ_MERGE_OUT_BASE_OVERRIDE="$OUTPUT_BASE"
      export RJ_SIMEMBED_DEST_BASE="${DEST_BASE}/simembedded"
      export RJ_SIMEMBEDINCLUSIVE_DEST_BASE="${DEST_BASE}/simembeddedinclusive"
      ./RecoilJets_Condor_submit.sh isSimEmbedded condorDoAll groupSize "$GROUP_SIZE"
      ./RecoilJets_Condor_submit.sh isSimEmbeddedInclusive condorDoAll groupSize "$GROUP_SIZE"
      condor_q -nobatch "${USER:-patsfan753}" | tail -80
    } | tee "$log"
    say "submit_log=${log}"
    ;;
  pullQA)
    print_plan
    say "pullQA should run only after Gmail/Condor reports READY for both paired campaigns."
    say "Use the standard sftp_get helper with RJ_CONFIG_YAML=${CONFIG_OUT}, then ROOT-open every merged file."
    ;;
  -h|--help|help)
    usage
    ;;
  *)
    usage
    die "Unknown mode: $mode"
    ;;
esac
