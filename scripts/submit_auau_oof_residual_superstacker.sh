#!/usr/bin/env bash
set -euo pipefail

RJ_REPO_BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
ML_PYTHON="${RJ_ML_PYTHON:-/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python}"
SOURCE="${RJ_AUAU_STACK_SOURCE:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049}"
MODEL_BASE="${RJ_AUAU_MLP_MODEL_BASE:-/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models}"
STAMP="${RJ_AUAU_SUPERSTACK_STAMP:-$(date +%Y%m%d_%H%M%S)}"
OUTDIR="${OUTDIR:-${MODEL_BASE}/oof_residual_superstacker_pt1535_${STAMP}}"
SUB_ROOT="${SUB_ROOT:-${RJ_REPO_BASE}/condor_sub/auauOOFResidualSuperStacker_${STAMP}}"
LOG="${LOG:-${OUTDIR}/oof_residual_superstacker_${STAMP}.log}"
MLP_CACHE="${MLP_CACHE:-${SOURCE}/reports/mlp_model_validation_condor_stack_full_features_uncapped_20260513_123934/score_caches.list}"
BDT_CACHE="${BDT_CACHE:-${SOURCE}/reports/model_validation_condor_stack_full_features_uncapped_20260513_123934/score_caches.list}"
MLP_SCORE="${MLP_SCORE:-score_centInputBase3x3WidthRatiosMLP_pt1535}"
BDT_SCORE="${BDT_SCORE:-score_ptFine_cent7}"
REQUEST_MEMORY="${REQUEST_MEMORY:-24000MB}"
REQUEST_CPUS="${REQUEST_CPUS:-1}"
MAX_ROWS="${MAX_ROWS:-0}"
MAX_SHARDS="${MAX_SHARDS:-0}"
REQUIRE_FULL_STAT="${REQUIRE_FULL_STAT:-1}"
EXPECTED_SHARDS="${EXPECTED_SHARDS:-80}"
FOLDS="${FOLDS:-5}"
EPOCHS="${EPOCHS:-180}"
PATIENCE="${PATIENCE:-28}"
BATCH_SIZE="${BATCH_SIZE:-32768}"
HIDDEN="${HIDDEN:-16}"
LOWER_ALGORITHMS="${LOWER_ALGORITHMS:-logistic,gbm,nn}"
BASE_SCORE="${BASE_SCORE:-nn}"
LOWER_FEATURE_MODE="${LOWER_FEATURE_MODE:-score_context}"
SUPER_FEATURE_MODE="${SUPER_FEATURE_MODE:-scores_context}"
FINAL_MODE="${FINAL_MODE:-residual}"
MODEL_NAME="${MODEL_NAME:-}"
TRAIN_FRACTION="${TRAIN_FRACTION:-0.80}"
VAL_FRACTION="${VAL_FRACTION:-0.10}"
INCLUDE_ISOLATION_CONTEXT="${INCLUDE_ISOLATION_CONTEXT:-0}"

say() { printf '\033[1;36m[auauSuperStack]\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[auauSuperStack][WARN]\033[0m %s\n' "$*" >&2; }
die() { printf '\033[1;31m[auauSuperStack][ERR]\033[0m %s\n' "$*" >&2; exit 2; }

usage() {
  cat <<'EOF'
Usage:
  ./scripts/submit_auau_oof_residual_superstacker.sh selfTest
  ./scripts/submit_auau_oof_residual_superstacker.sh train [OUTDIR=/path] [MAX_ROWS=N]
  ./scripts/submit_auau_oof_residual_superstacker.sh status OUTDIR=/path [SUB_ROOT=/path]

Runs the out-of-fold residual NN super-stacker diagnostic.  This is a ceiling
test using BDT/MLP score inputs.  Production is not submitted by this driver.
For the stricter heavy superlearner use:
  FINAL_MODE=direct LOWER_FEATURE_MODE=full_features SUPER_FEATURE_MODE=scores_plus_full_features HIDDEN=256,128,64
EOF
}

mode="${1:-train}"
if [[ $# -gt 0 ]]; then shift; fi
for tok in "$@"; do
  case "$tok" in
    OUTDIR=*) OUTDIR="${tok#OUTDIR=}" ;;
    SUB_ROOT=*) SUB_ROOT="${tok#SUB_ROOT=}" ;;
    LOG=*) LOG="${tok#LOG=}" ;;
    MLP_CACHE=*) MLP_CACHE="${tok#MLP_CACHE=}" ;;
    BDT_CACHE=*) BDT_CACHE="${tok#BDT_CACHE=}" ;;
    MLP_SCORE=*) MLP_SCORE="${tok#MLP_SCORE=}" ;;
    BDT_SCORE=*) BDT_SCORE="${tok#BDT_SCORE=}" ;;
    REQUEST_MEMORY=*) REQUEST_MEMORY="${tok#REQUEST_MEMORY=}" ;;
    REQUEST_CPUS=*) REQUEST_CPUS="${tok#REQUEST_CPUS=}" ;;
    MAX_ROWS=*) MAX_ROWS="${tok#MAX_ROWS=}" ;;
    MAX_SHARDS=*) MAX_SHARDS="${tok#MAX_SHARDS=}" ;;
    REQUIRE_FULL_STAT=*) REQUIRE_FULL_STAT="${tok#REQUIRE_FULL_STAT=}" ;;
    EXPECTED_SHARDS=*) EXPECTED_SHARDS="${tok#EXPECTED_SHARDS=}" ;;
    FOLDS=*) FOLDS="${tok#FOLDS=}" ;;
    EPOCHS=*) EPOCHS="${tok#EPOCHS=}" ;;
    PATIENCE=*) PATIENCE="${tok#PATIENCE=}" ;;
    BATCH_SIZE=*) BATCH_SIZE="${tok#BATCH_SIZE=}" ;;
    HIDDEN=*) HIDDEN="${tok#HIDDEN=}" ;;
    LOWER_ALGORITHMS=*) LOWER_ALGORITHMS="${tok#LOWER_ALGORITHMS=}" ;;
    BASE_SCORE=*) BASE_SCORE="${tok#BASE_SCORE=}" ;;
    LOWER_FEATURE_MODE=*) LOWER_FEATURE_MODE="${tok#LOWER_FEATURE_MODE=}" ;;
    SUPER_FEATURE_MODE=*) SUPER_FEATURE_MODE="${tok#SUPER_FEATURE_MODE=}" ;;
    FINAL_MODE=*) FINAL_MODE="${tok#FINAL_MODE=}" ;;
    MODEL_NAME=*) MODEL_NAME="${tok#MODEL_NAME=}" ;;
    TRAIN_FRACTION=*) TRAIN_FRACTION="${tok#TRAIN_FRACTION=}" ;;
    VAL_FRACTION=*) VAL_FRACTION="${tok#VAL_FRACTION=}" ;;
    INCLUDE_ISOLATION_CONTEXT=*) INCLUDE_ISOLATION_CONTEXT="${tok#INCLUDE_ISOLATION_CONTEXT=}" ;;
    -h|--help|help) usage; exit 0 ;;
    *) die "Unknown argument: $tok" ;;
  esac
done

print_plan() {
  say "RECOILJETS_AUAU_OOF_RESIDUAL_SUPERSTACKER_SUBMIT_V1"
  say "host=$(hostname -f 2>/dev/null || hostname)"
  say "timestamp=$STAMP"
  say "repo=$RJ_REPO_BASE"
  say "outdir=$OUTDIR"
  say "sub_root=$SUB_ROOT"
  say "log=$LOG"
  say "mlp_cache=$MLP_CACHE"
  say "bdt_cache=$BDT_CACHE"
  say "mlp_score=$MLP_SCORE bdt_score=$BDT_SCORE"
  say "request_memory=$REQUEST_MEMORY request_cpus=$REQUEST_CPUS"
  say "max_rows=$MAX_ROWS max_shards=$MAX_SHARDS require_full_stat=$REQUIRE_FULL_STAT expected_shards=$EXPECTED_SHARDS folds=$FOLDS lower_algorithms=$LOWER_ALGORITHMS base_score=$BASE_SCORE"
  say "lower_feature_mode=$LOWER_FEATURE_MODE super_feature_mode=$SUPER_FEATURE_MODE final_mode=$FINAL_MODE model_name=${MODEL_NAME:-<auto>}"
  say "hidden=$HIDDEN epochs=$EPOCHS patience=$PATIENCE batch_size=$BATCH_SIZE"
  say "train_fraction=$TRAIN_FRACTION val_fraction=$VAL_FRACTION include_isolation_context=$INCLUDE_ISOLATION_CONTEXT"
}

python_run="$ML_PYTHON"
if [[ ! -x "$python_run" ]]; then
  python_run="${RJ_LOCAL_PYTHON:-python3}"
fi

case "$mode" in
  selfTest|self-test)
    print_plan
    mkdir -p "$OUTDIR/self_test"
    "$python_run" "${RJ_REPO_BASE}/scripts/train_auau_oof_residual_superstacker.py" \
      --outdir "$OUTDIR/self_test" \
      --self-test \
      --folds 3 \
      --epochs 12 \
      --patience 4 \
      --hidden "$HIDDEN" \
      --lower-algorithms "$LOWER_ALGORITHMS" \
      --base-score "$BASE_SCORE" \
      --lower-feature-mode "$LOWER_FEATURE_MODE" \
      --super-feature-mode "$SUPER_FEATURE_MODE" \
      --final-mode "$FINAL_MODE" \
      --model-name "$MODEL_NAME" \
      --train-fraction "$TRAIN_FRACTION" \
      --val-fraction "$VAL_FRACTION" \
      $(if [[ "$INCLUDE_ISOLATION_CONTEXT" == "1" ]]; then printf '%s' "--include-isolation-context"; fi)
    ;;
  train)
    print_plan
    [[ -s "$MLP_CACHE" ]] || die "Missing MLP_CACHE=$MLP_CACHE"
    [[ -s "$BDT_CACHE" ]] || die "Missing BDT_CACHE=$BDT_CACHE"
    mkdir -p "$OUTDIR" "$SUB_ROOT" "$(dirname "$LOG")"
    worker="${SUB_ROOT}/run_oof_residual_superstacker.sh"
    cat > "$worker" <<EOF
#!/usr/bin/env bash
set -euo pipefail
cd "$RJ_REPO_BASE"
export RJ_ML_PYTHON="$ML_PYTHON"
WORKER_USER="\${USER:-patsfan753}"
MYINSTALL="/sphenix/u/\${WORKER_USER}/thesisAnalysis/install"
MYINSTALL_AUAU="/sphenix/u/\${WORKER_USER}/thesisAnalysis_auau/install"
set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n
[[ -d "\$MYINSTALL" ]] && source /opt/sphenix/core/bin/setup_local.sh "\$MYINSTALL" || true
[[ -d "\$MYINSTALL_AUAU" ]] && source /opt/sphenix/core/bin/setup_local.sh "\$MYINSTALL_AUAU" || true
set -u
ml_python="$ML_PYTHON"
ml_prefix="\$(cd "\$(dirname "\$ml_python")/.." && pwd -P 2>/dev/null || true)"
[[ -d "\$ml_prefix/lib" ]] && export LD_LIBRARY_PATH="\$ml_prefix/lib:\${LD_LIBRARY_PATH:-}"
echo "RECOILJETS_AUAU_OOF_RESIDUAL_SUPERSTACKER_WORKER_V1"
echo "host=\$(hostname -f 2>/dev/null || hostname)"
echo "start=\$(date)"
full_stat_args=()
if [[ "$REQUIRE_FULL_STAT" == "1" ]]; then
  full_stat_args+=(--require-full-stat --expected-shards "$EXPECTED_SHARDS")
fi
isolation_args=()
if [[ "$INCLUDE_ISOLATION_CONTEXT" == "1" ]]; then
  isolation_args+=(--include-isolation-context)
fi
"\$ml_python" scripts/train_auau_oof_residual_superstacker.py \\
  --mlp-cache "$MLP_CACHE" \\
  --bdt-cache "$BDT_CACHE" \\
  --outdir "$OUTDIR" \\
  --mlp-score "$MLP_SCORE" \\
  --bdt-score "$BDT_SCORE" \\
  --max-rows "$MAX_ROWS" \\
  --max-shards "$MAX_SHARDS" \\
  "\${full_stat_args[@]}" \\
  "\${isolation_args[@]}" \\
  --train-fraction "$TRAIN_FRACTION" \\
  --val-fraction "$VAL_FRACTION" \\
  --folds "$FOLDS" \\
  --lower-algorithms "$LOWER_ALGORITHMS" \\
  --base-score "$BASE_SCORE" \\
  --lower-feature-mode "$LOWER_FEATURE_MODE" \\
  --super-feature-mode "$SUPER_FEATURE_MODE" \\
  --final-mode "$FINAL_MODE" \\
  --model-name "$MODEL_NAME" \\
  --hidden "$HIDDEN" \\
  --epochs "$EPOCHS" \\
  --patience "$PATIENCE" \\
  --batch-size "$BATCH_SIZE"
echo "finish=\$(date)"
EOF
    chmod +x "$worker"
    sub="${SUB_ROOT}/oof_residual_superstacker.sub"
    cat > "$sub" <<EOF
universe = vanilla
executable = ${worker}
output = ${SUB_ROOT}/oof_residual_superstacker.out
error = ${SUB_ROOT}/oof_residual_superstacker.err
log = ${SUB_ROOT}/oof_residual_superstacker.log
request_memory = ${REQUEST_MEMORY}
request_cpus = ${REQUEST_CPUS}
+JobBatchName = "auau_oof_residual_superstacker_${STAMP}"
notification = Never
queue
EOF
    say "worker=$worker"
    say "submit=$sub"
    if [[ "${RJ_DAG_DRYRUN:-0}" == "1" ]]; then
      say "RJ_DAG_DRYRUN=1; not submitting"
      exit 0
    fi
    command -v condor_submit >/dev/null 2>&1 || die "condor_submit not found"
    condor_submit "$sub" | tee "$LOG"
    say "rank_table=${OUTDIR}/oof_residual_supernn_rank_table.csv"
    say "metrics=${OUTDIR}/oof_residual_supernn_metrics.json"
    ;;
  status)
    print_plan
    for path in \
      "$OUTDIR/oof_residual_supernn_rank_table.csv" \
      "$OUTDIR/oof_direct_supernn_rank_table.csv" \
      "$OUTDIR/oof_residual_supernn_metrics.json" \
      "$OUTDIR/oof_direct_supernn_metrics.json" \
      "$OUTDIR/oof_residual_supernn_training_history.csv" \
      "$OUTDIR/oof_direct_supernn_training_history.csv" \
      "$OUTDIR/oof_residual_supernn_manifest.json"; do
      if [[ -s "$path" ]]; then
        say "exists: $path"
        [[ "$path" == *.csv ]] && head -8 "$path" || true
      else
        warn "missing: $path"
      fi
    done
    [[ -d "$SUB_ROOT" ]] && find "$SUB_ROOT" -maxdepth 1 -type f -printf '%TY-%Tm-%Td %TH:%TM %p\n' | sort | tail -20 || true
    ;;
  -h|--help|help)
    usage
    ;;
  *)
    die "Unknown mode: $mode"
    ;;
esac
