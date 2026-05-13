#!/usr/bin/env bash
set -euo pipefail

RJ_REPO_BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
readonly RJ_REPO_BASE
TRAIN_BASE="${RJ_AUAU_TIGHT_MLP_TRAIN_BASE:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining}"
MODEL_BASE="${RJ_AUAU_MLP_MODEL_BASE:-${RJ_REPO_BASE}/mlp_models}"
TRAIN_SCRIPT="${RJ_AUAU_TIGHT_MLP_TRAIN_SCRIPT:-${RJ_REPO_BASE}/scripts/train_auau_photon_mlp.py}"
VALIDATE_SCRIPT="${RJ_AUAU_TIGHT_MLP_VALIDATE_SCRIPT:-${RJ_REPO_BASE}/scripts/validate_auau_tight_mlp_on_sim.py}"
ML_PYTHON="${RJ_ML_PYTHON:-python3}"
NOTIFY_EMAILS="${RJ_NOTIFY_EMAILS:-just0131@gmail.com}"

ts() { date +%Y%m%d_%H%M%S; }
say() { printf '\033[1;36m[auauTightMLP]\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[auauTightMLP][WARN]\033[0m %s\n' "$*" >&2; }
err() { printf '\033[1;31m[auauTightMLP][ERR]\033[0m %s\n' "$*" >&2; }
die() { err "$*"; exit 2; }

setup_sphenix_stack_env() {
  export USER="${USER:-$(id -u -n)}"
  export LOGNAME="${LOGNAME:-$USER}"
  export HOME="/sphenix/u/${LOGNAME}"
  local myinstall="/sphenix/u/${USER}/thesisAnalysis/install"
  local myinstall_auau="/sphenix/u/${USER}/thesisAnalysis_auau/install"
  set +u
  # shellcheck disable=SC1091
  source /opt/sphenix/core/bin/sphenix_setup.sh -n
  [[ -d "$myinstall" ]] && source /opt/sphenix/core/bin/setup_local.sh "$myinstall" || true
  [[ -d "$myinstall_auau" ]] && source /opt/sphenix/core/bin/setup_local.sh "$myinstall_auau" || true
  set -u
}

setup_ml_python_env() {
  setup_sphenix_stack_env
  local ml_python_prefix="" ml_python_real="" ml_python_real_prefix=""
  if [[ "$ML_PYTHON" == */* ]]; then
    ml_python_prefix="$(cd "$(dirname "$ML_PYTHON")/.." && pwd -P 2>/dev/null || true)"
    ml_python_real="$(readlink -f "$ML_PYTHON" 2>/dev/null || true)"
    [[ -n "$ml_python_real" ]] && ml_python_real_prefix="$(cd "$(dirname "$ml_python_real")/.." && pwd -P 2>/dev/null || true)"
  fi
  local ld_joined="" d
  for d in "$ml_python_prefix/lib" "$ml_python_prefix/lib64" "$ml_python_real_prefix/lib" "$ml_python_real_prefix/lib64"; do
    [[ -n "$d" && -d "$d" ]] || continue
    case ":$ld_joined:" in *":$d:"*) ;; *) ld_joined="${ld_joined:+$ld_joined:}$d" ;; esac
  done
  [[ -n "$ld_joined" ]] && export LD_LIBRARY_PATH="$ld_joined:${LD_LIBRARY_PATH:-}"
  unset PYTHONHOME
}

send_summary_email() {
  local subject="$1"
  local summary="$2"
  if command -v mail >/dev/null 2>&1 && [[ -s "$summary" ]]; then
    mail -s "$subject" "$NOTIFY_EMAILS" < "$summary" || true
  fi
}

guard_generated_path() {
  local label="$1" path="$2"
  [[ -n "$path" ]] || die "$label resolved to an empty path"
  case "$path" in
    /cvmfs/*|/opt/sphenix/*|/usr/*|/bin/*|/lib/*|/lib64/*|/etc/*)
      die "$label resolved to protected path: $path"
      ;;
  esac
}

make_root_manifest() {
  local search_root="$1" out="$2"
  mkdir -p "$(dirname "$out")"
  find "$search_root" -type f -name '*.root' | sort -V > "$out"
}

usage() {
  cat <<'EOF'
Usage:
  ./scripts/auau_tight_mlp_pipeline.sh trainFromExtraction SOURCE=/path [MODEL_DIR=/path] [PRODUCTS=primary|all]
  ./scripts/auau_tight_mlp_pipeline.sh trainPrimaryDeepFromExtraction SOURCE=/path [MODEL_DIR=/path]
  ./scripts/auau_tight_mlp_pipeline.sh trainHighPtBalancedFromExtraction SOURCE=/path [MODEL_DIR=/path]
  ./scripts/auau_tight_mlp_pipeline.sh trainKitchenSinkFromExtraction SOURCE=/path [MODEL_DIR=/path]
  ./scripts/auau_tight_mlp_pipeline.sh trainIsoKitchenSinkFromExtraction SOURCE=/path [MODEL_DIR=/path]
  ./scripts/auau_tight_mlp_pipeline.sh trainHighPtDistilledKitchenV2FromExtraction SOURCE=/path [MODEL_DIR=/path]
  ./scripts/auau_tight_mlp_pipeline.sh smokeTrainFromExtraction SOURCE=/path [MODEL_DIR=/path]
  ./scripts/auau_tight_mlp_pipeline.sh applyCheck MODEL_DIR=/path
  ./scripts/auau_tight_mlp_pipeline.sh validateOnSim SOURCE=/path MODEL_DIR=/path [OUTDIR=/path]
  ./scripts/auau_tight_mlp_pipeline.sh validateOnSimCondor SOURCE=/path MODEL_DIR=/path [groupSize N]
  ./scripts/auau_tight_mlp_pipeline.sh rescoreValidationCache CACHE=/path/score_caches.list OUTDIR=/path [MODEL_DIR=/path|SWEEP_MANIFEST=/path]
  ./scripts/auau_tight_mlp_pipeline.sh deriveWorkingPointsFromValidation VALIDATION=/path [TARGET=0.80]
  ./scripts/auau_tight_mlp_pipeline.sh generateWorkingPointConfig TEMPLATE=/path WORKING_POINTS=/path OUT=/path [MODEL_DIR=/path]

Sidecar AuAu tight-MLP workflow. It consumes the same AuAuPhotonIDTrainingTree
extraction products as the BDT pipeline and writes runtime JSON artifacts plus
MLP-specific working-point YAML fragments.
EOF
}

train_from_extraction() {
  local source="" model_dir="" products="${RJ_AUAU_TIGHT_MLP_PRODUCTS:-all}" pt_range="${RJ_AUAU_TIGHT_MLP_PT_RANGE:-15:35}"
  local centrality_range="${RJ_AUAU_MLP_TRAIN_CENTRALITY_RANGE:-}"
  local train_pt_bins="${RJ_AUAU_MLP_TRAIN_PT_BINS:-}"
  local max_rows="${RJ_AUAU_MLP_TRAIN_MAX_ROWS:-0}"
  local max_rows_per_class="${RJ_AUAU_MLP_TRAIN_MAX_ROWS_PER_CLASS:-0}"
  local max_rows_per_pt_bin_class="${RJ_AUAU_MLP_TRAIN_MAX_ROWS_PER_PT_BIN_CLASS:-0}"
  local max_files_per_sample="${RJ_AUAU_MLP_TRAIN_MAX_FILES_PER_SAMPLE:-0}"
  local pt_bin_weight_mode="${RJ_AUAU_MLP_TRAIN_PT_BIN_WEIGHT_MODE:-none}"
  local pt_bin_weight_spec="${RJ_AUAU_MLP_TRAIN_PT_BIN_WEIGHT_SPEC:-}"
  local highpt_selection_weights="${RJ_AUAU_MLP_TRAIN_HIGHPT_SELECTION_WEIGHTS:-}"
  local hard_example_branch="${RJ_AUAU_MLP_TRAIN_HARD_EXAMPLE_BRANCH:-}"
  local hard_background_factor="${RJ_AUAU_MLP_TRAIN_HARD_BACKGROUND_FACTOR:-0}"
  local hard_signal_factor="${RJ_AUAU_MLP_TRAIN_HARD_SIGNAL_FACTOR:-0}"
  local hard_example_power="${RJ_AUAU_MLP_TRAIN_HARD_EXAMPLE_POWER:-2.0}"
  local distillation_branch="${RJ_AUAU_MLP_TRAIN_DISTILLATION_BRANCH:-}"
  local distillation_strength="${RJ_AUAU_MLP_TRAIN_DISTILLATION_STRENGTH:-0.0}"
  local distillation_temperature="${RJ_AUAU_MLP_TRAIN_DISTILLATION_TEMPERATURE:-1.0}"
  local distillation_min_finite_fraction="${RJ_AUAU_MLP_TRAIN_DISTILLATION_MIN_FINITE_FRACTION:-0.95}"
  local require_distillation="${RJ_AUAU_MLP_TRAIN_REQUIRE_DISTILLATION:-0}"
  local epochs="${RJ_AUAU_MLP_TRAIN_EPOCHS:-80}"
  local patience="${RJ_AUAU_MLP_TRAIN_PATIENCE:-12}"
  local progress_every="${RJ_AUAU_MLP_TRAIN_PROGRESS_EVERY:-5}"
  local restarts="${RJ_AUAU_MLP_TRAIN_RESTARTS:-1}"
  local hidden_layers="${RJ_AUAU_MLP_TRAIN_HIDDEN_LAYERS:-48,24}"
  local hidden_grid="${RJ_AUAU_MLP_TRAIN_HIDDEN_LAYER_GRID:-}"
  local selection_metric="${RJ_AUAU_MLP_TRAIN_SELECTION_METRIC:-wp_fake_rate}"
  local selection_target="${RJ_AUAU_MLP_TRAIN_SELECTION_TARGET_SIGNAL_EFF:-0.80}"
  local batch_size="${RJ_AUAU_MLP_TRAIN_BATCH_SIZE:-4096}"
  local learning_rate="${RJ_AUAU_MLP_TRAIN_LEARNING_RATE:-1.0e-3}"
  local l2="${RJ_AUAU_MLP_TRAIN_L2:-1.0e-4}"
  local min_delta="${RJ_AUAU_MLP_TRAIN_MIN_DELTA:-1.0e-5}"
  local conditional_jitter="${RJ_AUAU_MLP_TRAIN_CONDITIONAL_JITTER:-0.03}"
  local input_clip="${RJ_AUAU_MLP_TRAIN_INPUT_CLIP:-8.0}"
  local tok
  for tok in "$@"; do
    case "$tok" in
      SOURCE=*) source="${tok#SOURCE=}" ;;
      MODEL_DIR=*|MODELDIR=*|model_dir=*) model_dir="${tok#*=}" ;;
      PRODUCTS=*|products=*) products="${tok#*=}" ;;
      PT_RANGE=*|ptRange=*|pt_range=*) pt_range="${tok#*=}" ;;
      CENTRALITY_RANGE=*|centralityRange=*|centrality_range=*) centrality_range="${tok#*=}" ;;
      TRAIN_PT_BINS=*|trainPtBins=*|train_pt_bins=*) train_pt_bins="${tok#*=}" ;;
      MAX_ROWS=*|maxRows=*) max_rows="${tok#*=}" ;;
      MAX_ROWS_PER_CLASS=*|maxRowsPerClass=*) max_rows_per_class="${tok#*=}" ;;
      MAX_ROWS_PER_PT_BIN_CLASS=*|maxRowsPerPtBinClass=*) max_rows_per_pt_bin_class="${tok#*=}" ;;
      MAX_FILES_PER_SAMPLE=*|maxFilesPerSample=*) max_files_per_sample="${tok#*=}" ;;
      PT_BIN_WEIGHT_MODE=*|ptBinWeightMode=*) pt_bin_weight_mode="${tok#*=}" ;;
      PT_BIN_WEIGHT_SPEC=*|ptBinWeightSpec=*) pt_bin_weight_spec="${tok#*=}" ;;
      HIGHPT_SELECTION_WEIGHTS=*|highptSelectionWeights=*) highpt_selection_weights="${tok#*=}" ;;
      HARD_EXAMPLE_BRANCH=*|hardExampleBranch=*) hard_example_branch="${tok#*=}" ;;
      HARD_BACKGROUND_FACTOR=*|hardBackgroundFactor=*) hard_background_factor="${tok#*=}" ;;
      HARD_SIGNAL_FACTOR=*|hardSignalFactor=*) hard_signal_factor="${tok#*=}" ;;
      HARD_EXAMPLE_POWER=*|hardExamplePower=*) hard_example_power="${tok#*=}" ;;
      DISTILLATION_BRANCH=*|distillationBranch=*) distillation_branch="${tok#*=}" ;;
      DISTILLATION_STRENGTH=*|distillationStrength=*) distillation_strength="${tok#*=}" ;;
      DISTILLATION_TEMPERATURE=*|distillationTemperature=*) distillation_temperature="${tok#*=}" ;;
      DISTILLATION_MIN_FINITE_FRACTION=*|distillationMinFiniteFraction=*) distillation_min_finite_fraction="${tok#*=}" ;;
      REQUIRE_DISTILLATION=*|requireDistillation=*) require_distillation="${tok#*=}" ;;
      EPOCHS=*|epochs=*) epochs="${tok#*=}" ;;
      PATIENCE=*|patience=*) patience="${tok#*=}" ;;
      PROGRESS_EVERY=*|progressEvery=*) progress_every="${tok#*=}" ;;
      RESTARTS=*|restarts=*) restarts="${tok#*=}" ;;
      HIDDEN_LAYERS=*|hiddenLayers=*) hidden_layers="${tok#*=}" ;;
      HIDDEN_LAYER_GRID=*|hiddenLayerGrid=*) hidden_grid="${tok#*=}" ;;
      SELECTION_METRIC=*|selectionMetric=*) selection_metric="${tok#*=}" ;;
      SELECTION_TARGET_SIGNAL_EFF=*|selectionTargetSignalEff=*) selection_target="${tok#*=}" ;;
      BATCH_SIZE=*|batchSize=*) batch_size="${tok#*=}" ;;
      LEARNING_RATE=*|learningRate=*) learning_rate="${tok#*=}" ;;
      L2=*|l2=*) l2="${tok#*=}" ;;
      MIN_DELTA=*|minDelta=*) min_delta="${tok#*=}" ;;
      CONDITIONAL_JITTER=*|conditionalJitter=*) conditional_jitter="${tok#*=}" ;;
      INPUT_CLIP=*|inputClip=*) input_clip="${tok#*=}" ;;
    esac
  done
  [[ -n "$source" && -d "$source" ]] || die "trainFromExtraction requires SOURCE=/path/to/extraction"
  [[ -s "$TRAIN_SCRIPT" ]] || die "Missing trainer: $TRAIN_SCRIPT"
  local stamp="${RJ_AUAU_TIGHT_MLP_STAMP:-$(ts)}"
  model_dir="${model_dir:-${MODEL_BASE}/tight_mlp_${stamp}}"
  guard_generated_path "model_dir" "$model_dir"
  mkdir -p "$model_dir"
  setup_ml_python_env
  say "Training AuAu tight-MLP products"
  say "  source   : $source"
  say "  model dir: $model_dir"
  say "  products : $products"
  say "  pt range : $pt_range"
  [[ -n "$centrality_range" ]] && say "  centrality range: $centrality_range"
  [[ -n "$train_pt_bins" ]] && say "  train pT bins: $train_pt_bins"
  say "  max files/sample: $max_files_per_sample"
  say "  max rows: $max_rows  max rows/class: $max_rows_per_class  max rows/pT-bin/class: $max_rows_per_pt_bin_class"
  say "  pT-bin weighting: $pt_bin_weight_mode ${pt_bin_weight_spec:+($pt_bin_weight_spec)}"
  say "  high-pT selection weights: ${highpt_selection_weights:-none}"
  [[ -n "$hard_example_branch" ]] && say "  hard-example weighting: branch=$hard_example_branch bkg=$hard_background_factor sig=$hard_signal_factor power=$hard_example_power"
  [[ -n "$distillation_branch" ]] && say "  distillation: branch=$distillation_branch strength=$distillation_strength temp=$distillation_temperature required=$require_distillation minFinite=$distillation_min_finite_fraction"
  say "  epochs/patience: $epochs/$patience"
  say "  hidden layers: $hidden_layers"
  say "  restarts: $restarts  selection: $selection_metric @ signal eff $selection_target"
  say "  batch/lr/l2: $batch_size / $learning_rate / $l2"
  say "  min_delta/jitter/input_clip: $min_delta / $conditional_jitter / $input_clip"
  [[ -n "$hidden_grid" ]] && say "  hidden grid: $hidden_grid"
  local -a train_args=(
    "$TRAIN_SCRIPT"
    --source "$source" \
    --outdir "$model_dir" \
    --products "$products" \
    --pt-range "$pt_range" \
    --epochs "$epochs" \
    --patience "$patience" \
    --progress-every "$progress_every" \
    --restarts "$restarts" \
    --hidden-layers "$hidden_layers" \
    --selection-metric "$selection_metric" \
    --selection-target-signal-efficiency "$selection_target" \
    --batch-size "$batch_size" \
    --learning-rate "$learning_rate" \
    --l2 "$l2" \
    --min-delta "$min_delta" \
    --conditional-jitter "$conditional_jitter" \
    --input-clip "$input_clip"
  )
  [[ -n "$centrality_range" ]] && train_args+=( --centrality-range "$centrality_range" )
  [[ -n "$train_pt_bins" ]] && train_args+=( --train-pt-bins "$train_pt_bins" )
  [[ "$max_rows_per_pt_bin_class" =~ ^[0-9]+$ && "$max_rows_per_pt_bin_class" -gt 0 ]] && train_args+=( --max-rows-per-pt-bin-class "$max_rows_per_pt_bin_class" )
  [[ "$pt_bin_weight_mode" != "none" ]] && train_args+=( --pt-bin-weight-mode "$pt_bin_weight_mode" )
  [[ -n "$pt_bin_weight_spec" ]] && train_args+=( --pt-bin-weight-spec "$pt_bin_weight_spec" )
  [[ -n "$highpt_selection_weights" ]] && train_args+=( --highpt-selection-weights "$highpt_selection_weights" )
  [[ -n "$hard_example_branch" ]] && train_args+=( --hard-example-branch "$hard_example_branch" --hard-background-factor "$hard_background_factor" --hard-signal-factor "$hard_signal_factor" --hard-example-power "$hard_example_power" )
  [[ -n "$distillation_branch" ]] && train_args+=( --distillation-branch "$distillation_branch" --distillation-strength "$distillation_strength" --distillation-temperature "$distillation_temperature" --distillation-min-finite-fraction "$distillation_min_finite_fraction" )
  [[ "$require_distillation" == "1" || "$require_distillation" == "true" || "$require_distillation" == "TRUE" || "$require_distillation" == "yes" ]] && train_args+=( --require-distillation )
  [[ -n "$hidden_grid" ]] && train_args+=( --hidden-layer-grid "$hidden_grid" )
  [[ "$max_rows" =~ ^[0-9]+$ && "$max_rows" -gt 0 ]] && train_args+=( --max-rows "$max_rows" )
  [[ "$max_rows_per_class" =~ ^[0-9]+$ && "$max_rows_per_class" -gt 0 ]] && train_args+=( --max-rows-per-class "$max_rows_per_class" )
  [[ "$max_files_per_sample" =~ ^[0-9]+$ && "$max_files_per_sample" -gt 0 ]] && train_args+=( --max-files-per-sample "$max_files_per_sample" )
  "$ML_PYTHON" "${train_args[@]}"
  say "MLP training registry: ${model_dir}/model_registry.json"
}

train_primary_deep_from_extraction() {
  export RJ_AUAU_TIGHT_MLP_PRODUCTS="${RJ_AUAU_TIGHT_MLP_PRODUCTS:-primary-ratios}"
  export RJ_AUAU_MLP_TRAIN_MAX_FILES_PER_SAMPLE="${RJ_AUAU_MLP_TRAIN_MAX_FILES_PER_SAMPLE:-0}"
  export RJ_AUAU_MLP_TRAIN_MAX_ROWS="${RJ_AUAU_MLP_TRAIN_MAX_ROWS:-0}"
  export RJ_AUAU_MLP_TRAIN_MAX_ROWS_PER_CLASS="${RJ_AUAU_MLP_TRAIN_MAX_ROWS_PER_CLASS:-350000}"
  export RJ_AUAU_MLP_TRAIN_EPOCHS="${RJ_AUAU_MLP_TRAIN_EPOCHS:-260}"
  export RJ_AUAU_MLP_TRAIN_PATIENCE="${RJ_AUAU_MLP_TRAIN_PATIENCE:-55}"
  export RJ_AUAU_MLP_TRAIN_PROGRESS_EVERY="${RJ_AUAU_MLP_TRAIN_PROGRESS_EVERY:-2}"
  export RJ_AUAU_MLP_TRAIN_RESTARTS="${RJ_AUAU_MLP_TRAIN_RESTARTS:-1}"
  export RJ_AUAU_MLP_TRAIN_HIDDEN_LAYERS="${RJ_AUAU_MLP_TRAIN_HIDDEN_LAYERS:-128,64,32}"
  export RJ_AUAU_MLP_TRAIN_HIDDEN_LAYER_GRID="${RJ_AUAU_MLP_TRAIN_HIDDEN_LAYER_GRID:-}"
  export RJ_AUAU_MLP_TRAIN_SELECTION_METRIC="${RJ_AUAU_MLP_TRAIN_SELECTION_METRIC:-validation_auc}"
  export RJ_AUAU_MLP_TRAIN_SELECTION_TARGET_SIGNAL_EFF="${RJ_AUAU_MLP_TRAIN_SELECTION_TARGET_SIGNAL_EFF:-0.80}"
  export RJ_AUAU_MLP_TRAIN_BATCH_SIZE="${RJ_AUAU_MLP_TRAIN_BATCH_SIZE:-4096}"
  export RJ_AUAU_MLP_TRAIN_LEARNING_RATE="${RJ_AUAU_MLP_TRAIN_LEARNING_RATE:-7.0e-4}"
  export RJ_AUAU_MLP_TRAIN_L2="${RJ_AUAU_MLP_TRAIN_L2:-5.0e-5}"
  export RJ_AUAU_MLP_TRAIN_MIN_DELTA="${RJ_AUAU_MLP_TRAIN_MIN_DELTA:-2.0e-6}"
  export RJ_AUAU_MLP_TRAIN_CONDITIONAL_JITTER="${RJ_AUAU_MLP_TRAIN_CONDITIONAL_JITTER:-0.015}"
  export RJ_AUAU_MLP_TRAIN_INPUT_CLIP="${RJ_AUAU_MLP_TRAIN_INPUT_CLIP:-7.0}"
  train_from_extraction "$@"
}

train_highpt_balanced_from_extraction() {
  export RJ_AUAU_TIGHT_MLP_PRODUCTS="${RJ_AUAU_TIGHT_MLP_PRODUCTS:-primary-ratios}"
  export RJ_AUAU_MLP_TRAIN_MAX_FILES_PER_SAMPLE="${RJ_AUAU_MLP_TRAIN_MAX_FILES_PER_SAMPLE:-0}"
  export RJ_AUAU_MLP_TRAIN_MAX_ROWS="${RJ_AUAU_MLP_TRAIN_MAX_ROWS:-0}"
  export RJ_AUAU_MLP_TRAIN_MAX_ROWS_PER_CLASS="${RJ_AUAU_MLP_TRAIN_MAX_ROWS_PER_CLASS:-0}"
  export RJ_AUAU_MLP_TRAIN_MAX_ROWS_PER_PT_BIN_CLASS="${RJ_AUAU_MLP_TRAIN_MAX_ROWS_PER_PT_BIN_CLASS:-120000}"
  export RJ_AUAU_MLP_TRAIN_PT_BIN_WEIGHT_MODE="${RJ_AUAU_MLP_TRAIN_PT_BIN_WEIGHT_MODE:-highpt}"
  export RJ_AUAU_MLP_TRAIN_SELECTION_METRIC="${RJ_AUAU_MLP_TRAIN_SELECTION_METRIC:-highpt_wp80}"
  export RJ_AUAU_MLP_TRAIN_SELECTION_TARGET_SIGNAL_EFF="${RJ_AUAU_MLP_TRAIN_SELECTION_TARGET_SIGNAL_EFF:-0.80}"
  export RJ_AUAU_MLP_TRAIN_EPOCHS="${RJ_AUAU_MLP_TRAIN_EPOCHS:-220}"
  export RJ_AUAU_MLP_TRAIN_PATIENCE="${RJ_AUAU_MLP_TRAIN_PATIENCE:-45}"
  export RJ_AUAU_MLP_TRAIN_PROGRESS_EVERY="${RJ_AUAU_MLP_TRAIN_PROGRESS_EVERY:-2}"
  export RJ_AUAU_MLP_TRAIN_RESTARTS="${RJ_AUAU_MLP_TRAIN_RESTARTS:-1}"
  export RJ_AUAU_MLP_TRAIN_HIDDEN_LAYERS="${RJ_AUAU_MLP_TRAIN_HIDDEN_LAYERS:-128,64,32}"
  export RJ_AUAU_MLP_TRAIN_HIDDEN_LAYER_GRID="${RJ_AUAU_MLP_TRAIN_HIDDEN_LAYER_GRID:-160,80,40}"
  export RJ_AUAU_MLP_TRAIN_BATCH_SIZE="${RJ_AUAU_MLP_TRAIN_BATCH_SIZE:-4096}"
  export RJ_AUAU_MLP_TRAIN_LEARNING_RATE="${RJ_AUAU_MLP_TRAIN_LEARNING_RATE:-7.0e-4}"
  export RJ_AUAU_MLP_TRAIN_L2="${RJ_AUAU_MLP_TRAIN_L2:-5.0e-5}"
  export RJ_AUAU_MLP_TRAIN_MIN_DELTA="${RJ_AUAU_MLP_TRAIN_MIN_DELTA:-2.0e-6}"
  export RJ_AUAU_MLP_TRAIN_CONDITIONAL_JITTER="${RJ_AUAU_MLP_TRAIN_CONDITIONAL_JITTER:-0.015}"
  export RJ_AUAU_MLP_TRAIN_INPUT_CLIP="${RJ_AUAU_MLP_TRAIN_INPUT_CLIP:-7.0}"
  train_from_extraction "$@"
}

train_kitchen_sink_from_extraction() {
  export RJ_AUAU_TIGHT_MLP_PRODUCTS="${RJ_AUAU_TIGHT_MLP_PRODUCTS:-kitchen-sink}"
  export RJ_AUAU_TIGHT_MLP_PT_RANGE="${RJ_AUAU_TIGHT_MLP_PT_RANGE:-5:35}"
  export RJ_AUAU_MLP_TRAIN_PT_BINS="${RJ_AUAU_MLP_TRAIN_PT_BINS:-5,10,15,20,25,35}"
  export RJ_AUAU_MLP_TRAIN_MAX_FILES_PER_SAMPLE="${RJ_AUAU_MLP_TRAIN_MAX_FILES_PER_SAMPLE:-0}"
  export RJ_AUAU_MLP_TRAIN_MAX_ROWS="${RJ_AUAU_MLP_TRAIN_MAX_ROWS:-0}"
  export RJ_AUAU_MLP_TRAIN_MAX_ROWS_PER_CLASS="${RJ_AUAU_MLP_TRAIN_MAX_ROWS_PER_CLASS:-0}"
  export RJ_AUAU_MLP_TRAIN_MAX_ROWS_PER_PT_BIN_CLASS="${RJ_AUAU_MLP_TRAIN_MAX_ROWS_PER_PT_BIN_CLASS:-80000}"
  export RJ_AUAU_MLP_TRAIN_PT_BIN_WEIGHT_MODE="${RJ_AUAU_MLP_TRAIN_PT_BIN_WEIGHT_MODE:-highpt}"
  export RJ_AUAU_MLP_TRAIN_PT_BIN_WEIGHT_SPEC="${RJ_AUAU_MLP_TRAIN_PT_BIN_WEIGHT_SPEC:-5:10:0.05,10:15:0.10,15:20:0.20,20:25:0.30,25:35:0.35}"
  export RJ_AUAU_MLP_TRAIN_HIGHPT_SELECTION_WEIGHTS="${RJ_AUAU_MLP_TRAIN_HIGHPT_SELECTION_WEIGHTS:-5:10:0.05,10:15:0.10,15:20:0.20,20:25:0.30,25:35:0.35}"
  export RJ_AUAU_MLP_TRAIN_SELECTION_METRIC="${RJ_AUAU_MLP_TRAIN_SELECTION_METRIC:-highpt_wp80}"
  export RJ_AUAU_MLP_TRAIN_SELECTION_TARGET_SIGNAL_EFF="${RJ_AUAU_MLP_TRAIN_SELECTION_TARGET_SIGNAL_EFF:-0.80}"
  export RJ_AUAU_MLP_TRAIN_HARD_EXAMPLE_BRANCH="${RJ_AUAU_MLP_TRAIN_HARD_EXAMPLE_BRANCH:-auau_tight_bdt_score}"
  export RJ_AUAU_MLP_TRAIN_HARD_BACKGROUND_FACTOR="${RJ_AUAU_MLP_TRAIN_HARD_BACKGROUND_FACTOR:-3.0}"
  export RJ_AUAU_MLP_TRAIN_HARD_SIGNAL_FACTOR="${RJ_AUAU_MLP_TRAIN_HARD_SIGNAL_FACTOR:-1.0}"
  export RJ_AUAU_MLP_TRAIN_HARD_EXAMPLE_POWER="${RJ_AUAU_MLP_TRAIN_HARD_EXAMPLE_POWER:-2.0}"
  export RJ_AUAU_MLP_TRAIN_EPOCHS="${RJ_AUAU_MLP_TRAIN_EPOCHS:-260}"
  export RJ_AUAU_MLP_TRAIN_PATIENCE="${RJ_AUAU_MLP_TRAIN_PATIENCE:-55}"
  export RJ_AUAU_MLP_TRAIN_PROGRESS_EVERY="${RJ_AUAU_MLP_TRAIN_PROGRESS_EVERY:-2}"
  export RJ_AUAU_MLP_TRAIN_RESTARTS="${RJ_AUAU_MLP_TRAIN_RESTARTS:-1}"
  export RJ_AUAU_MLP_TRAIN_HIDDEN_LAYERS="${RJ_AUAU_MLP_TRAIN_HIDDEN_LAYERS:-192,96,48}"
  export RJ_AUAU_MLP_TRAIN_HIDDEN_LAYER_GRID="${RJ_AUAU_MLP_TRAIN_HIDDEN_LAYER_GRID:-256,128,64}"
  export RJ_AUAU_MLP_TRAIN_BATCH_SIZE="${RJ_AUAU_MLP_TRAIN_BATCH_SIZE:-4096}"
  export RJ_AUAU_MLP_TRAIN_LEARNING_RATE="${RJ_AUAU_MLP_TRAIN_LEARNING_RATE:-5.0e-4}"
  export RJ_AUAU_MLP_TRAIN_L2="${RJ_AUAU_MLP_TRAIN_L2:-7.0e-5}"
  export RJ_AUAU_MLP_TRAIN_MIN_DELTA="${RJ_AUAU_MLP_TRAIN_MIN_DELTA:-2.0e-6}"
  export RJ_AUAU_MLP_TRAIN_CONDITIONAL_JITTER="${RJ_AUAU_MLP_TRAIN_CONDITIONAL_JITTER:-0.012}"
  export RJ_AUAU_MLP_TRAIN_INPUT_CLIP="${RJ_AUAU_MLP_TRAIN_INPUT_CLIP:-7.0}"
  train_from_extraction "$@"
}

train_iso_kitchen_sink_from_extraction() {
  export RJ_AUAU_TIGHT_MLP_PRODUCTS="${RJ_AUAU_TIGHT_MLP_PRODUCTS:-iso-kitchen-sink}"
  export RJ_AUAU_TIGHT_MLP_PT_RANGE="${RJ_AUAU_TIGHT_MLP_PT_RANGE:-5:35}"
  export RJ_AUAU_MLP_TRAIN_PT_BINS="${RJ_AUAU_MLP_TRAIN_PT_BINS:-5,10,15,20,25,35}"
  export RJ_AUAU_MLP_TRAIN_MAX_FILES_PER_SAMPLE="${RJ_AUAU_MLP_TRAIN_MAX_FILES_PER_SAMPLE:-0}"
  export RJ_AUAU_MLP_TRAIN_MAX_ROWS="${RJ_AUAU_MLP_TRAIN_MAX_ROWS:-0}"
  export RJ_AUAU_MLP_TRAIN_MAX_ROWS_PER_CLASS="${RJ_AUAU_MLP_TRAIN_MAX_ROWS_PER_CLASS:-0}"
  export RJ_AUAU_MLP_TRAIN_MAX_ROWS_PER_PT_BIN_CLASS="${RJ_AUAU_MLP_TRAIN_MAX_ROWS_PER_PT_BIN_CLASS:-80000}"
  export RJ_AUAU_MLP_TRAIN_PT_BIN_WEIGHT_MODE="${RJ_AUAU_MLP_TRAIN_PT_BIN_WEIGHT_MODE:-highpt}"
  export RJ_AUAU_MLP_TRAIN_PT_BIN_WEIGHT_SPEC="${RJ_AUAU_MLP_TRAIN_PT_BIN_WEIGHT_SPEC:-5:10:0.05,10:15:0.10,15:20:0.20,20:25:0.30,25:35:0.35}"
  export RJ_AUAU_MLP_TRAIN_HIGHPT_SELECTION_WEIGHTS="${RJ_AUAU_MLP_TRAIN_HIGHPT_SELECTION_WEIGHTS:-5:10:0.05,10:15:0.10,15:20:0.20,20:25:0.30,25:35:0.35}"
  export RJ_AUAU_MLP_TRAIN_SELECTION_METRIC="${RJ_AUAU_MLP_TRAIN_SELECTION_METRIC:-highpt_wp80}"
  export RJ_AUAU_MLP_TRAIN_SELECTION_TARGET_SIGNAL_EFF="${RJ_AUAU_MLP_TRAIN_SELECTION_TARGET_SIGNAL_EFF:-0.80}"
  export RJ_AUAU_MLP_TRAIN_HARD_EXAMPLE_BRANCH="${RJ_AUAU_MLP_TRAIN_HARD_EXAMPLE_BRANCH:-auau_tight_bdt_score}"
  export RJ_AUAU_MLP_TRAIN_HARD_BACKGROUND_FACTOR="${RJ_AUAU_MLP_TRAIN_HARD_BACKGROUND_FACTOR:-3.0}"
  export RJ_AUAU_MLP_TRAIN_HARD_SIGNAL_FACTOR="${RJ_AUAU_MLP_TRAIN_HARD_SIGNAL_FACTOR:-1.0}"
  export RJ_AUAU_MLP_TRAIN_HARD_EXAMPLE_POWER="${RJ_AUAU_MLP_TRAIN_HARD_EXAMPLE_POWER:-2.0}"
  export RJ_AUAU_MLP_TRAIN_EPOCHS="${RJ_AUAU_MLP_TRAIN_EPOCHS:-260}"
  export RJ_AUAU_MLP_TRAIN_PATIENCE="${RJ_AUAU_MLP_TRAIN_PATIENCE:-55}"
  export RJ_AUAU_MLP_TRAIN_PROGRESS_EVERY="${RJ_AUAU_MLP_TRAIN_PROGRESS_EVERY:-2}"
  export RJ_AUAU_MLP_TRAIN_RESTARTS="${RJ_AUAU_MLP_TRAIN_RESTARTS:-1}"
  export RJ_AUAU_MLP_TRAIN_HIDDEN_LAYERS="${RJ_AUAU_MLP_TRAIN_HIDDEN_LAYERS:-192,96,48}"
  export RJ_AUAU_MLP_TRAIN_HIDDEN_LAYER_GRID="${RJ_AUAU_MLP_TRAIN_HIDDEN_LAYER_GRID:-256,128,64}"
  export RJ_AUAU_MLP_TRAIN_BATCH_SIZE="${RJ_AUAU_MLP_TRAIN_BATCH_SIZE:-4096}"
  export RJ_AUAU_MLP_TRAIN_LEARNING_RATE="${RJ_AUAU_MLP_TRAIN_LEARNING_RATE:-5.0e-4}"
  export RJ_AUAU_MLP_TRAIN_L2="${RJ_AUAU_MLP_TRAIN_L2:-7.0e-5}"
  export RJ_AUAU_MLP_TRAIN_MIN_DELTA="${RJ_AUAU_MLP_TRAIN_MIN_DELTA:-2.0e-6}"
  export RJ_AUAU_MLP_TRAIN_CONDITIONAL_JITTER="${RJ_AUAU_MLP_TRAIN_CONDITIONAL_JITTER:-0.012}"
  export RJ_AUAU_MLP_TRAIN_INPUT_CLIP="${RJ_AUAU_MLP_TRAIN_INPUT_CLIP:-7.0}"
  say "NOTE: iso-aware kitchen-sink MLP is diagnostic-only; it uses reco isolation-derived inputs and is not ABCD-purity safe."
  train_from_extraction "$@"
}

train_highpt_distilled_kitchen_v2_from_extraction() {
  export RJ_AUAU_TIGHT_MLP_PRODUCTS="${RJ_AUAU_TIGHT_MLP_PRODUCTS:-highpt-distilled-kitchen-v2}"
  export RJ_AUAU_TIGHT_MLP_PT_RANGE="${RJ_AUAU_TIGHT_MLP_PT_RANGE:-15:35}"
  export RJ_AUAU_MLP_TRAIN_PT_BINS="${RJ_AUAU_MLP_TRAIN_PT_BINS:-15,20,25,35}"
  export RJ_AUAU_MLP_TRAIN_MAX_FILES_PER_SAMPLE="${RJ_AUAU_MLP_TRAIN_MAX_FILES_PER_SAMPLE:-0}"
  export RJ_AUAU_MLP_TRAIN_MAX_ROWS="${RJ_AUAU_MLP_TRAIN_MAX_ROWS:-0}"
  export RJ_AUAU_MLP_TRAIN_MAX_ROWS_PER_CLASS="${RJ_AUAU_MLP_TRAIN_MAX_ROWS_PER_CLASS:-0}"
  export RJ_AUAU_MLP_TRAIN_MAX_ROWS_PER_PT_BIN_CLASS="${RJ_AUAU_MLP_TRAIN_MAX_ROWS_PER_PT_BIN_CLASS:-120000}"
  export RJ_AUAU_MLP_TRAIN_PT_BIN_WEIGHT_MODE="${RJ_AUAU_MLP_TRAIN_PT_BIN_WEIGHT_MODE:-highpt}"
  export RJ_AUAU_MLP_TRAIN_PT_BIN_WEIGHT_SPEC="${RJ_AUAU_MLP_TRAIN_PT_BIN_WEIGHT_SPEC:-15:20:0.15,20:25:0.35,25:35:0.50}"
  export RJ_AUAU_MLP_TRAIN_HIGHPT_SELECTION_WEIGHTS="${RJ_AUAU_MLP_TRAIN_HIGHPT_SELECTION_WEIGHTS:-15:20:0.15,20:25:0.35,25:35:0.50}"
  export RJ_AUAU_MLP_TRAIN_SELECTION_METRIC="${RJ_AUAU_MLP_TRAIN_SELECTION_METRIC:-highpt_wp80}"
  export RJ_AUAU_MLP_TRAIN_SELECTION_TARGET_SIGNAL_EFF="${RJ_AUAU_MLP_TRAIN_SELECTION_TARGET_SIGNAL_EFF:-0.80}"
  export RJ_AUAU_MLP_TRAIN_HARD_EXAMPLE_BRANCH="${RJ_AUAU_MLP_TRAIN_HARD_EXAMPLE_BRANCH:-auau_tight_bdt_score}"
  export RJ_AUAU_MLP_TRAIN_HARD_BACKGROUND_FACTOR="${RJ_AUAU_MLP_TRAIN_HARD_BACKGROUND_FACTOR:-2.0}"
  export RJ_AUAU_MLP_TRAIN_HARD_SIGNAL_FACTOR="${RJ_AUAU_MLP_TRAIN_HARD_SIGNAL_FACTOR:-0.8}"
  export RJ_AUAU_MLP_TRAIN_HARD_EXAMPLE_POWER="${RJ_AUAU_MLP_TRAIN_HARD_EXAMPLE_POWER:-2.0}"
  export RJ_AUAU_MLP_TRAIN_DISTILLATION_BRANCH="${RJ_AUAU_MLP_TRAIN_DISTILLATION_BRANCH:-auau_tight_bdt_score}"
  export RJ_AUAU_MLP_TRAIN_DISTILLATION_STRENGTH="${RJ_AUAU_MLP_TRAIN_DISTILLATION_STRENGTH:-0.15}"
  export RJ_AUAU_MLP_TRAIN_DISTILLATION_TEMPERATURE="${RJ_AUAU_MLP_TRAIN_DISTILLATION_TEMPERATURE:-1.0}"
  export RJ_AUAU_MLP_TRAIN_DISTILLATION_MIN_FINITE_FRACTION="${RJ_AUAU_MLP_TRAIN_DISTILLATION_MIN_FINITE_FRACTION:-0.95}"
  export RJ_AUAU_MLP_TRAIN_REQUIRE_DISTILLATION="${RJ_AUAU_MLP_TRAIN_REQUIRE_DISTILLATION:-1}"
  export RJ_AUAU_MLP_TRAIN_EPOCHS="${RJ_AUAU_MLP_TRAIN_EPOCHS:-280}"
  export RJ_AUAU_MLP_TRAIN_PATIENCE="${RJ_AUAU_MLP_TRAIN_PATIENCE:-60}"
  export RJ_AUAU_MLP_TRAIN_PROGRESS_EVERY="${RJ_AUAU_MLP_TRAIN_PROGRESS_EVERY:-2}"
  export RJ_AUAU_MLP_TRAIN_RESTARTS="${RJ_AUAU_MLP_TRAIN_RESTARTS:-2}"
  export RJ_AUAU_MLP_TRAIN_HIDDEN_LAYERS="${RJ_AUAU_MLP_TRAIN_HIDDEN_LAYERS:-192,96,48}"
  export RJ_AUAU_MLP_TRAIN_HIDDEN_LAYER_GRID="${RJ_AUAU_MLP_TRAIN_HIDDEN_LAYER_GRID:-160,80,40}"
  export RJ_AUAU_MLP_TRAIN_BATCH_SIZE="${RJ_AUAU_MLP_TRAIN_BATCH_SIZE:-4096}"
  export RJ_AUAU_MLP_TRAIN_LEARNING_RATE="${RJ_AUAU_MLP_TRAIN_LEARNING_RATE:-4.0e-4}"
  export RJ_AUAU_MLP_TRAIN_L2="${RJ_AUAU_MLP_TRAIN_L2:-8.0e-5}"
  export RJ_AUAU_MLP_TRAIN_MIN_DELTA="${RJ_AUAU_MLP_TRAIN_MIN_DELTA:-2.0e-6}"
  export RJ_AUAU_MLP_TRAIN_CONDITIONAL_JITTER="${RJ_AUAU_MLP_TRAIN_CONDITIONAL_JITTER:-0.012}"
  export RJ_AUAU_MLP_TRAIN_INPUT_CLIP="${RJ_AUAU_MLP_TRAIN_INPUT_CLIP:-7.0}"
  say "NOTE: v2 uses the BDT score only as a training teacher/hard-example guide; the exported artifact has no BDT-score runtime input."
  train_from_extraction "$@"
}

smoke_train_from_extraction() {
  export RJ_AUAU_TIGHT_MLP_PRODUCTS="${RJ_AUAU_TIGHT_MLP_PRODUCTS:-primary}"
  export RJ_AUAU_MLP_TRAIN_MAX_FILES_PER_SAMPLE="${RJ_AUAU_MLP_TRAIN_MAX_FILES_PER_SAMPLE:-300}"
  export RJ_AUAU_MLP_TRAIN_MAX_ROWS="${RJ_AUAU_MLP_TRAIN_MAX_ROWS:-200000}"
  export RJ_AUAU_MLP_TRAIN_MAX_ROWS_PER_CLASS="${RJ_AUAU_MLP_TRAIN_MAX_ROWS_PER_CLASS:-100000}"
  export RJ_AUAU_MLP_TRAIN_EPOCHS="${RJ_AUAU_MLP_TRAIN_EPOCHS:-80}"
  export RJ_AUAU_MLP_TRAIN_PATIENCE="${RJ_AUAU_MLP_TRAIN_PATIENCE:-12}"
  export RJ_AUAU_MLP_TRAIN_RESTARTS="${RJ_AUAU_MLP_TRAIN_RESTARTS:-3}"
  export RJ_AUAU_MLP_TRAIN_HIDDEN_LAYER_GRID="${RJ_AUAU_MLP_TRAIN_HIDDEN_LAYER_GRID:-64,32;96,48;64,32,16}"
  export RJ_AUAU_MLP_TRAIN_SELECTION_METRIC="${RJ_AUAU_MLP_TRAIN_SELECTION_METRIC:-wp_fake_rate}"
  export RJ_AUAU_MLP_TRAIN_SELECTION_TARGET_SIGNAL_EFF="${RJ_AUAU_MLP_TRAIN_SELECTION_TARGET_SIGNAL_EFF:-0.80}"
  train_from_extraction "$@"
}

apply_check() {
  local model_dir=""
  local tok
  for tok in "$@"; do
    case "$tok" in MODEL_DIR=*|MODELDIR=*|model_dir=*) model_dir="${tok#*=}" ;; esac
  done
  [[ -n "$model_dir" && -d "$model_dir" ]] || die "applyCheck requires MODEL_DIR=/path"
  [[ -s "${model_dir}/model_registry.json" ]] || die "Missing model registry: ${model_dir}/model_registry.json"
  setup_ml_python_env
  "$ML_PYTHON" - "$model_dir" <<'PY'
import json, math, sys
from pathlib import Path
sys.path.insert(0, str(Path("scripts").resolve()))
from train_auau_photon_mlp import load_mlp_artifact, predict_mlp_array
model_dir = Path(sys.argv[1])
reg = json.loads((model_dir / "model_registry.json").read_text())
bad = []
for spec in reg.get("models", []):
    path = model_dir / Path(spec["output_json"]).name
    try:
        model = load_mlp_artifact(path)
        x = [[0.0 for _ in model["features"]]]
        y = float(predict_mlp_array(model, x)[0])
        if not math.isfinite(y) or y < 0.0 or y > 1.0:
            bad.append(f"{path}: non-finite/out-of-range smoke score {y}")
    except Exception as exc:
        bad.append(f"{path}: {exc}")
if bad:
    raise SystemExit("MLP applyCheck failed:\n  " + "\n  ".join(bad))
print(f"[OK] applyCheck loaded and scored {len(reg.get('models', []))} MLP JSON artifacts")
PY
}

validate_on_sim() {
  local source="" model_dir="" model_registry="" outdir=""
  local score_max="${RJ_AUAU_TIGHT_MLP_VALIDATE_SCORE_MAX_ROWS:-300000}"
  local max_files_per_sample="${RJ_AUAU_TIGHT_MLP_VALIDATE_MAX_FILES_PER_SAMPLE:-0}"
  local pt_bins="${RJ_AUAU_TIGHT_MLP_VALIDATE_PT_BINS:-}"
  local centrality_bins="${RJ_AUAU_TIGHT_MLP_VALIDATE_CENTRALITY_BINS:-}"
  local tok
  for tok in "$@"; do
    case "$tok" in
      SOURCE=*) source="${tok#SOURCE=}" ;;
      MODEL_DIR=*|MODELDIR=*|model_dir=*) model_dir="${tok#*=}" ;;
      MODEL_REGISTRY=*|REGISTRY=*|model_registry=*) model_registry="${tok#*=}" ;;
      OUTDIR=*|outdir=*) outdir="${tok#*=}" ;;
      SCORE_MAX_ROWS=*|scoreMaxRows=*) score_max="${tok#*=}" ;;
      MAX_FILES_PER_SAMPLE=*|maxFilesPerSample=*) max_files_per_sample="${tok#*=}" ;;
      PT_BINS=*|ptBins=*) pt_bins="${tok#*=}" ;;
      CENTRALITY_BINS=*|centralityBins=*) centrality_bins="${tok#*=}" ;;
    esac
  done
  [[ -n "$source" && -d "$source" ]] || die "validateOnSim requires SOURCE=/path"
  [[ -n "$model_dir" && -d "$model_dir" ]] || die "validateOnSim requires MODEL_DIR=/path"
  setup_ml_python_env
  local -a args=( "$VALIDATE_SCRIPT" --source "$source" --model-dir "$model_dir" )
  [[ -n "$model_registry" ]] && args+=( --model-registry "$model_registry" )
  [[ -n "$outdir" ]] && args+=( --outdir "$outdir" )
  [[ -n "$pt_bins" ]] && args+=( --pt-bins "$pt_bins" )
  [[ -n "$centrality_bins" ]] && args+=( --centrality-bins "$centrality_bins" )
  [[ "$score_max" =~ ^[0-9]+$ && "$score_max" -gt 0 ]] && args+=( --score-max-rows "$score_max" )
  [[ "$max_files_per_sample" =~ ^[0-9]+$ && "$max_files_per_sample" -gt 0 ]] && args+=( --max-files-per-sample "$max_files_per_sample" )
  say "Validation score rows: $score_max  max files/sample: $max_files_per_sample"
  "$ML_PYTHON" "${args[@]}"
  local summary="${outdir:-${source}/reports/mlp_model_validation}/validation_summary.txt"
  if [[ -s "$summary" ]]; then
    local status
    status="$(awk -F= '/^status=/ {print $2; exit}' "$summary")"
    send_summary_email "[RecoilJets][auauTightMLP_validateOnSim][${status:-CHECK}]" "$summary"
    say "validation summary: $summary"
  fi
}

validate_on_sim_condor() {
  local source="" model_dir="" model_registry="" outdir="" group_size="${RJ_AUAU_TIGHT_MLP_VALIDATE_GROUP_SIZE:-100}"
  local total_score_max="${RJ_AUAU_TIGHT_MLP_VALIDATE_TOTAL_SCORE_MAX_ROWS:-400000}"
  local reqmem="${RJ_AUAU_TIGHT_MLP_VALIDATE_REQUEST_MEMORY:-2500MB}"
  local tok
  for tok in "$@"; do
    case "$tok" in
      SOURCE=*) source="${tok#SOURCE=}" ;;
      MODEL_DIR=*|MODELDIR=*|model_dir=*) model_dir="${tok#*=}" ;;
      MODEL_REGISTRY=*|REGISTRY=*|model_registry=*) model_registry="${tok#*=}" ;;
      OUTDIR=*|outdir=*) outdir="${tok#*=}" ;;
      groupSize=*) group_size="${tok#groupSize=}" ;;
      SCORE_MAX_ROWS=*|scoreMaxRows=*) total_score_max="${tok#*=}" ;;
    esac
  done
  while (($#)); do
    case "$1" in
      groupSize) group_size="${2:?missing value after groupSize}"; shift 2 ;;
      scoreMaxRows|totalScoreMaxRows) total_score_max="${2:?missing value after $1}"; shift 2 ;;
      *) shift ;;
    esac
  done
  [[ -n "$source" && -d "$source" ]] || die "validateOnSimCondor requires SOURCE=/path"
  [[ -n "$model_dir" && -d "$model_dir" ]] || die "validateOnSimCondor requires MODEL_DIR=/path"
  [[ "$group_size" =~ ^[0-9]+$ && "$group_size" -gt 0 ]] || die "groupSize must be positive"
  local stamp="${RJ_AUAU_TIGHT_MLP_VALIDATE_STAMP:-$(ts)}"
  local report_root="${outdir:-${source}/reports/mlp_model_validation_condor_${stamp}}"
  local sub_root="${RJ_REPO_BASE}/condor_sub/auauTightMLPValidate_${stamp}"
  local shard_dir="${sub_root}/shards" cache_dir="${report_root}/score_caches"
  guard_generated_path "validation report root" "$report_root"
  guard_generated_path "validation submit root" "$sub_root"
  mkdir -p "$report_root" "$sub_root" "$shard_dir" "$cache_dir"
  local root_manifest="${source}/manifests/training_roots.list"
  local search_root="$source"
  [[ -d "${source}/extraction" ]] && search_root="${source}/extraction"
  make_root_manifest "$search_root" "$root_manifest"
  local nroots
  nroots="$(wc -l < "$root_manifest" | tr -d ' ')"
  [[ "$nroots" != "0" ]] || die "No ROOT files available for validation"
  split -l "$group_size" -d -a 5 "$root_manifest" "${shard_dir}/roots_"
  local shard_count
  shard_count="$(find "$shard_dir" -maxdepth 1 -type f -name 'roots_*' | wc -l | tr -d ' ')"
  local score_max_per_shard=0
  if [[ "$total_score_max" =~ ^[0-9]+$ && "$total_score_max" -gt 0 ]]; then
    score_max_per_shard=$(( (total_score_max + shard_count - 1) / shard_count ))
  fi
  local worker="${sub_root}/validate_worker.sh"
  cat > "$worker" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export RJ_ML_PYTHON="${ML_PYTHON}"
echo "[auauTightMLP] worker env setup start" >&2
export USER="\${USER:-\$(id -u -n)}"
export LOGNAME="\${LOGNAME:-\$USER}"
export HOME="/sphenix/u/\${LOGNAME}"
MYINSTALL="/sphenix/u/\${USER}/thesisAnalysis/install"
MYINSTALL_AUAU="/sphenix/u/\${USER}/thesisAnalysis_auau/install"
set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n
echo "[auauTightMLP] worker sphenix_setup rc=\$?" >&2
if [[ -d "\$MYINSTALL" ]]; then
  source /opt/sphenix/core/bin/setup_local.sh "\$MYINSTALL" || true
  echo "[auauTightMLP] worker setup_local rc=\$?" >&2
fi
if [[ -d "\$MYINSTALL_AUAU" ]]; then
  source /opt/sphenix/core/bin/setup_local.sh "\$MYINSTALL_AUAU" || true
  echo "[auauTightMLP] worker setup_local_auau rc=\$?" >&2
fi
set -u
ml_python="${ML_PYTHON}"
ml_python_prefix="\$(cd "\$(dirname "\$ml_python")/.." && pwd -P 2>/dev/null || true)"
ml_python_real="\$(readlink -f "\$ml_python" 2>/dev/null || true)"
ml_python_real_prefix=""
if [[ -n "\$ml_python_real" ]]; then
  ml_python_real_prefix="\$(cd "\$(dirname "\$ml_python_real")/.." && pwd -P 2>/dev/null || true)"
fi
ld_joined=""
for d in "\$ml_python_prefix/lib" "\$ml_python_prefix/lib64" "\$ml_python_real_prefix/lib" "\$ml_python_real_prefix/lib64"; do
  [[ -n "\$d" && -d "\$d" ]] || continue
  case ":\$ld_joined:" in *":\$d:"*) ;; *) ld_joined="\${ld_joined:+\$ld_joined:}\$d" ;; esac
done
[[ -n "\$ld_joined" ]] && export LD_LIBRARY_PATH="\$ld_joined:\${LD_LIBRARY_PATH:-}"
unset PYTHONHOME
echo "[auauTightMLP] worker env setup done" >&2
cd "${RJ_REPO_BASE}"
manifest="\${1:?manifest}"
outdir="\${2:?outdir}"
cache="\${3:?cache}"
score_max="\${4:?score_max}"
model_registry="${model_registry}"
model_registry_arg=()
[[ -n "\$model_registry" ]] && model_registry_arg=(--model-registry "\$model_registry")
"\$ml_python" "${VALIDATE_SCRIPT}" --source "${source}" --model-dir "${model_dir}" "\${model_registry_arg[@]}" --manifest "\$manifest" --outdir "\$outdir" --write-score-cache "\$cache" --score-max-rows "\$score_max" --no-plots
EOF
  chmod +x "$worker"
  local args_file="${sub_root}/validate_args.txt"
  : > "$args_file"
  local idx=0 shard
  for shard in "${shard_dir}/roots_"*; do
    [[ -s "$shard" ]] || continue
    idx=$((idx + 1))
    printf '%s %s %s %s\n' "$shard" "${report_root}/shards/shard_$(printf '%05d' "$idx")" "${cache_dir}/score_cache_$(printf '%05d' "$idx").npz" "$score_max_per_shard" >> "$args_file"
  done
  local dag="${sub_root}/auau_tight_mlp_validateOnSimCondor.dag"
  : > "$dag"
  local shard_index=0 shard_manifest shard_out shard_cache shard_scoremax
  while read -r shard_manifest shard_out shard_cache shard_scoremax; do
    shard_index=$((shard_index + 1))
    mkdir -p "$shard_out"
    local label sub
    label="$(printf '%05d' "$shard_index")"
    sub="${sub_root}/validate_${label}.sub"
    cat > "$sub" <<EOF
universe = vanilla
executable = ${worker}
arguments = ${shard_manifest} ${shard_out} ${shard_cache} ${shard_scoremax}
output = ${sub_root}/validate_${label}.out
error = ${sub_root}/validate_${label}.err
log = ${sub_root}/validate_${label}.log
request_memory = ${reqmem}
notification = Never
queue
EOF
    echo "JOB VALIDATE_${label} ${sub}" >> "$dag"
    echo "RETRY VALIDATE_${label} 1" >> "$dag"
  done < "$args_file"
  local merge="${sub_root}/validate_merge.sh"
  cat > "$merge" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export RJ_ML_PYTHON="${ML_PYTHON}"
echo "[auauTightMLP] merge env setup start" >&2
export USER="\${USER:-\$(id -u -n)}"
export LOGNAME="\${LOGNAME:-\$USER}"
export HOME="/sphenix/u/\${LOGNAME}"
MYINSTALL="/sphenix/u/\${USER}/thesisAnalysis/install"
MYINSTALL_AUAU="/sphenix/u/\${USER}/thesisAnalysis_auau/install"
set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n
echo "[auauTightMLP] merge sphenix_setup rc=\$?" >&2
if [[ -d "\$MYINSTALL" ]]; then
  source /opt/sphenix/core/bin/setup_local.sh "\$MYINSTALL" || true
  echo "[auauTightMLP] merge setup_local rc=\$?" >&2
fi
if [[ -d "\$MYINSTALL_AUAU" ]]; then
  source /opt/sphenix/core/bin/setup_local.sh "\$MYINSTALL_AUAU" || true
  echo "[auauTightMLP] merge setup_local_auau rc=\$?" >&2
fi
set -u
ml_python="${ML_PYTHON}"
ml_python_prefix="\$(cd "\$(dirname "\$ml_python")/.." && pwd -P 2>/dev/null || true)"
ml_python_real="\$(readlink -f "\$ml_python" 2>/dev/null || true)"
ml_python_real_prefix=""
if [[ -n "\$ml_python_real" ]]; then
  ml_python_real_prefix="\$(cd "\$(dirname "\$ml_python_real")/.." && pwd -P 2>/dev/null || true)"
fi
ld_joined=""
for d in "\$ml_python_prefix/lib" "\$ml_python_prefix/lib64" "\$ml_python_real_prefix/lib" "\$ml_python_real_prefix/lib64"; do
  [[ -n "\$d" && -d "\$d" ]] || continue
  case ":\$ld_joined:" in *":\$d:"*) ;; *) ld_joined="\${ld_joined:+\$ld_joined:}\$d" ;; esac
done
[[ -n "\$ld_joined" ]] && export LD_LIBRARY_PATH="\$ld_joined:\${LD_LIBRARY_PATH:-}"
unset PYTHONHOME
echo "[auauTightMLP] merge env setup done" >&2
cd "${RJ_REPO_BASE}"
cache_manifest="${report_root}/score_caches.list"
find "${cache_dir}" -type f -name 'score_cache_*.npz' | sort -V > "\$cache_manifest" || true
expected=${idx}
found=\$(wc -l < "\$cache_manifest" | tr -d ' ')
summary="${report_root}/validation_summary.txt"
if [[ "\$found" != "\$expected" || "\$found" == "0" ]]; then
  {
    echo "RECOILJETS_AUAU_TIGHT_MLP_SIM_VALIDATION_V1"
    echo "status=CHECK"
    echo "source=${source}"
    echo "model_dir=${model_dir}"
    echo "report_dir=${report_root}"
    echo "expected_score_caches=\$expected"
    echo "found_score_caches=\$found"
    echo "notes=missing score cache shards"
  } > "\$summary"
else
  model_registry="${model_registry}"
  model_registry_arg=()
  [[ -n "\$model_registry" ]] && model_registry_arg=(--model-registry "\$model_registry")
  "\$ml_python" "${VALIDATE_SCRIPT}" --source "${source}" --model-dir "${model_dir}" "\${model_registry_arg[@]}" --merge-score-caches "\$cache_manifest" --outdir "${report_root}"
fi
status=\$(awk -F= '/^status=/ {print \$2; exit}' "\$summary" 2>/dev/null || echo CHECK)
if command -v mail >/dev/null 2>&1 && [[ -s "\$summary" ]]; then
  mail -s "[RecoilJets][auauTightMLP_validateOnSimCondor][\${status:-CHECK}]" "${NOTIFY_EMAILS}" < "\$summary" || true
fi
cat "\$summary"
[[ "\${status:-CHECK}" == "READY" ]] || exit 3
EOF
  chmod +x "$merge"
  local merge_sub="${sub_root}/validate_merge.sub"
  cat > "$merge_sub" <<EOF
universe = scheduler
executable = ${merge}
output = ${sub_root}/validate_merge.out
error = ${sub_root}/validate_merge.err
log = ${sub_root}/validate_merge.log
notification = Never
queue
EOF
  echo "FINAL MERGE ${merge_sub}" >> "$dag"
  say "Condor validation DAG: $dag"
  say "source=$source  model_dir=$model_dir  report=$report_root"
  say "root files=${nroots} shards=${idx} groupSize=${group_size} scoreMaxPerShard=${score_max_per_shard} request_memory=${reqmem}"
  if [[ "${RJ_DAG_DRYRUN:-0}" == "1" ]]; then
    echo "RECOILJETS_AUAU_TIGHT_MLP_VALIDATE_DRYRUN_V1"
    echo "dag=${dag}"
    echo "report_root=${report_root}"
    return 0
  fi
  condor_submit_dag "$dag"
}

rescore_validation_cache() {
  local cache_manifest="" model_dir="" model_registry="" sweep_manifest="" outdir=""
  local pt_bins="${RJ_AUAU_TIGHT_MLP_VALIDATE_PT_BINS:-}"
  local centrality_bins="${RJ_AUAU_TIGHT_MLP_VALIDATE_CENTRALITY_BINS:-}"
  local tok
  for tok in "$@"; do
    case "$tok" in
      CACHE=*|CACHE_MANIFEST=*|cache=*) cache_manifest="${tok#*=}" ;;
      MODEL_DIR=*|MODELDIR=*|model_dir=*) model_dir="${tok#*=}" ;;
      MODEL_REGISTRY=*|REGISTRY=*|model_registry=*) model_registry="${tok#*=}" ;;
      SWEEP_MANIFEST=*|sweepManifest=*) sweep_manifest="${tok#*=}" ;;
      OUTDIR=*|outdir=*) outdir="${tok#*=}" ;;
      PT_BINS=*|ptBins=*) pt_bins="${tok#*=}" ;;
      CENTRALITY_BINS=*|centralityBins=*) centrality_bins="${tok#*=}" ;;
    esac
  done
  [[ -n "$cache_manifest" && -s "$cache_manifest" ]] || die "rescoreValidationCache requires CACHE=/path/to/score_caches.list"
  [[ -n "$outdir" ]] || die "rescoreValidationCache requires OUTDIR=/path"
  if [[ -z "$sweep_manifest" ]]; then
    [[ -n "$model_dir" && -d "$model_dir" ]] || die "rescoreValidationCache requires MODEL_DIR=/path or SWEEP_MANIFEST=/path"
  fi
  setup_ml_python_env
  local -a args=( "$VALIDATE_SCRIPT" --rescore-score-caches "$cache_manifest" --outdir "$outdir" )
  [[ -n "$model_dir" ]] && args+=( --model-dir "$model_dir" )
  [[ -n "$model_registry" ]] && args+=( --model-registry "$model_registry" )
  [[ -n "$sweep_manifest" ]] && args+=( --sweep-manifest "$sweep_manifest" )
  [[ -n "$pt_bins" ]] && args+=( --pt-bins "$pt_bins" )
  [[ -n "$centrality_bins" ]] && args+=( --centrality-bins "$centrality_bins" )
  say "Rescoring validation cache"
  say "  cache : $cache_manifest"
  [[ -n "$model_dir" ]] && say "  model : $model_dir"
  [[ -n "$sweep_manifest" ]] && say "  sweep : $sweep_manifest"
  say "  outdir: $outdir"
  "$ML_PYTHON" "${args[@]}"
}

derive_working_points_from_validation() {
  local validation="" target="0.80" tok
  for tok in "$@"; do
    case "$tok" in VALIDATION=*) validation="${tok#VALIDATION=}" ;; TARGET=*) target="${tok#TARGET=}" ;; esac
  done
  [[ -n "$validation" && -d "$validation" ]] || die "deriveWorkingPointsFromValidation requires VALIDATION=/path"
  setup_ml_python_env
  "$ML_PYTHON" "$VALIDATE_SCRIPT" --derive-working-points-from-report "$validation" --target-signal-efficiency "$target"
}

generate_working_point_config() {
  local template="" wp="" out="" model_dir="" tok
  for tok in "$@"; do
    case "$tok" in
      TEMPLATE=*) template="${tok#TEMPLATE=}" ;;
      WORKING_POINTS=*) wp="${tok#WORKING_POINTS=}" ;;
      OUT=*) out="${tok#OUT=}" ;;
      MODEL_DIR=*) model_dir="${tok#MODEL_DIR=}" ;;
    esac
  done
  [[ -n "$template" ]] || die "generateWorkingPointConfig requires TEMPLATE=/path"
  [[ -n "$wp" ]] || die "generateWorkingPointConfig requires WORKING_POINTS=/path"
  [[ -n "$out" ]] || die "generateWorkingPointConfig requires OUT=/path"
  local cmd=(python3 "${RJ_REPO_BASE}/scripts/make_auau_mlp_target_wp_config.py" --template "$template" --working-points "$wp" --out "$out")
  [[ -n "$model_dir" ]] && cmd+=(--model-dir "$model_dir")
  "${cmd[@]}"
}

main() {
  local mode="${1:-}"
  [[ -n "$mode" ]] || { usage; exit 2; }
  shift || true
  case "$mode" in
    -h|--help|help) usage ;;
    trainFromExtraction) train_from_extraction "$@" ;;
    trainPrimaryDeepFromExtraction|trainPrimaryDeep|trainPrimaryRatiosDeepFromExtraction|trainPrimaryRatiosDeep) train_primary_deep_from_extraction "$@" ;;
    trainHighPtBalancedFromExtraction|trainHighPtBalanced|trainHighPt) train_highpt_balanced_from_extraction "$@" ;;
    trainKitchenSinkFromExtraction|trainKitchenSink|kitchenSink) train_kitchen_sink_from_extraction "$@" ;;
    trainIsoKitchenSinkFromExtraction|trainIsoKitchenSink|isoKitchenSink|isoDiagnostic) train_iso_kitchen_sink_from_extraction "$@" ;;
    trainHighPtDistilledKitchenV2FromExtraction|trainHighPtDistilledKitchenV2|trainDistilledKitchenV2|trainBDTBeatingV2) train_highpt_distilled_kitchen_v2_from_extraction "$@" ;;
    smokeTrainFromExtraction|smokeTrain) smoke_train_from_extraction "$@" ;;
    applyCheck|smokeTestApplyExisting) apply_check "$@" ;;
    validateOnSim|validateSim|simValidation) validate_on_sim "$@" ;;
    validateOnSimCondor|condorValidateOnSim|validateSimCondor) validate_on_sim_condor "$@" ;;
    rescoreValidationCache|rescoreCache|validateFromCache) rescore_validation_cache "$@" ;;
    deriveWorkingPointsFromValidation|deriveWPFromValidation) derive_working_points_from_validation "$@" ;;
    generateWorkingPointConfig|generateTargetWPConfig) generate_working_point_config "$@" ;;
    *) usage; die "Unknown mode: $mode" ;;
  esac
}

main "$@"
