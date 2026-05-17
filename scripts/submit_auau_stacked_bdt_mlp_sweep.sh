#!/usr/bin/env bash
set -euo pipefail

RJ_REPO_BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
readonly RJ_REPO_BASE

ML_PYTHON="${RJ_ML_PYTHON:-/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python}"
SOURCE="${RJ_AUAU_STACK_SOURCE:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049}"
MODEL_BASE="${RJ_AUAU_MLP_MODEL_BASE:-/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models}"
STAMP="${RJ_AUAU_STACK_SWEEP_STAMP:-$(date +%Y%m%d_%H%M%S)}"
OUTDIR="${RJ_AUAU_STACK_SWEEP_OUTDIR:-${MODEL_BASE}/stacked_bdt_mlp_full_feature_sweep_${STAMP}}"
SUB_ROOT="${RJ_AUAU_STACK_SWEEP_SUB_ROOT:-${RJ_REPO_BASE}/condor_sub/auauStackedBDTMLPFullFeature_${STAMP}}"
LOG="${RJ_AUAU_STACK_SWEEP_LOG:-${OUTDIR}/stacked_bdt_mlp_full_feature_sweep_${STAMP}.log}"
LOG_EXPLICIT=0
[[ -n "${RJ_AUAU_STACK_SWEEP_LOG+x}" ]] && LOG_EXPLICIT=1
MLP_CACHE="${RJ_AUAU_STACK_MLP_CACHE:-${SOURCE}/reports/mlp_model_validation_condor_deep_primary_ratios_nostat_fullval_20260512_145449/score_caches.list}"
BDT_CACHE="${RJ_AUAU_STACK_BDT_CACHE:-${SOURCE}/reports/model_validation_condor_20260511_203441/score_caches.list}"
LOGREG_CACHE="${RJ_AUAU_STACK_LOGREG_CACHE:-}"
MLP_SCORE="${RJ_AUAU_STACK_MLP_SCORE:-score_centInputBase3x3WidthRatiosMLP_pt1535}"
BDT_SCORE="${RJ_AUAU_STACK_BDT_SCORE:-score_ptFine_cent7}"
LOGREG_SCORE="${RJ_AUAU_STACK_LOGREG_SCORE:-score_centEtFullLogReg_pt1535}"
REQUEST_MEMORY="${RJ_AUAU_STACK_SWEEP_REQUEST_MEMORY:-12000MB}"
REQUEST_CPUS="${RJ_AUAU_STACK_SWEEP_REQUEST_CPUS:-1}"
MAX_ROWS="${RJ_AUAU_STACK_SWEEP_MAX_ROWS:-0}"
MAX_SHARDS="${RJ_AUAU_STACK_SWEEP_MAX_SHARDS:-0}"
ALGORITHMS="${RJ_AUAU_STACK_SWEEP_ALGORITHMS:-logistic,gbm}"
VARIANTS="${RJ_AUAU_STACK_SWEEP_VARIANTS:-}"
SWEEP="${RJ_AUAU_STACK_SWEEP_PRESET:-full}"
INCLUDE_ISOLATION_CONTEXT="${RJ_AUAU_STACK_INCLUDE_ISOLATION_CONTEXT:-0}"
TOP_N="${RJ_AUAU_STACK_SWEEP_TOP_N:-4}"
ALLOW_MISSING_FULL_FEATURES="${RJ_AUAU_STACK_ALLOW_MISSING_FULL_FEATURES:-0}"
LINEAR_BACKEND="${RJ_AUAU_STACK_LINEAR_BACKEND:-numpy}"
NN_HIDDEN="${RJ_AUAU_STACK_NN_HIDDEN:-64,32}"
NN_EPOCHS="${RJ_AUAU_STACK_NN_EPOCHS:-180}"
NN_PATIENCE="${RJ_AUAU_STACK_NN_PATIENCE:-24}"
NN_BATCH_SIZE="${RJ_AUAU_STACK_NN_BATCH_SIZE:-8192}"
NN_LEARNING_RATE="${RJ_AUAU_STACK_NN_LEARNING_RATE:-0.0015}"
NN_L2="${RJ_AUAU_STACK_NN_L2:-0.001}"

say() { printf '\033[1;36m[auauStackSweep]\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[auauStackSweep][WARN]\033[0m %s\n' "$*" >&2; }
err() { printf '\033[1;31m[auauStackSweep][ERR]\033[0m %s\n' "$*" >&2; }
die() { err "$*"; exit 2; }

usage() {
  cat <<'EOF'
Usage:
  ./scripts/submit_auau_stacked_bdt_mlp_sweep.sh preflight [MLP_CACHE=/path] [BDT_CACHE=/path]
  ./scripts/submit_auau_stacked_bdt_mlp_sweep.sh selfTest
  ./scripts/submit_auau_stacked_bdt_mlp_sweep.sh train [MLP_CACHE=/path] [BDT_CACHE=/path]
  ./scripts/submit_auau_stacked_bdt_mlp_sweep.sh status OUTDIR=/path

Runs the diagnostic full-feature BDT+MLP stacker sweep.  This is a ceiling test:
the stackers use the runtime BDT score and are not ABCD-safe production photon
ID candidates unless separately reviewed.

Useful knobs:
  RJ_ML_PYTHON=/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python
  OUTDIR=/gpfs/.../stacked_bdt_mlp_full_feature_sweep_<stamp>
  MLP_CACHE=/path/score_caches.list
  BDT_CACHE=/path/score_caches.list
  LOGREG_CACHE=/path/score_caches.list
  MAX_ROWS=0 MAX_SHARDS=0
  ALGORITHMS=logistic,gbm
  VARIANTS=ptFine5to35_cent7_full,ptFine15to35_cent7_full
  For neural side tests: ALGORITHMS=nn or ALGORITHMS=logistic,gbm,nn
  LINEAR_BACKEND=numpy|sklearn
  NN_HIDDEN=64,32 NN_EPOCHS=180 NN_PATIENCE=24
  SWEEP=full|compact
  ALLOW_MISSING_FULL_FEATURES=1
  INCLUDE_ISOLATION_CONTEXT=1
  RJ_DAG_DRYRUN=1
EOF
}

mode="${1:-train}"
if [[ $# -gt 0 ]]; then shift; fi

for tok in "$@"; do
  case "$tok" in
    SOURCE=*) SOURCE="${tok#SOURCE=}" ;;
    OUTDIR=*) OUTDIR="${tok#OUTDIR=}" ;;
    SUB_ROOT=*) SUB_ROOT="${tok#SUB_ROOT=}" ;;
    LOG=*) LOG="${tok#LOG=}"; LOG_EXPLICIT=1 ;;
    MLP_CACHE=*) MLP_CACHE="${tok#MLP_CACHE=}" ;;
    BDT_CACHE=*) BDT_CACHE="${tok#BDT_CACHE=}" ;;
    LOGREG_CACHE=*) LOGREG_CACHE="${tok#LOGREG_CACHE=}" ;;
    MLP_SCORE=*) MLP_SCORE="${tok#MLP_SCORE=}" ;;
    BDT_SCORE=*) BDT_SCORE="${tok#BDT_SCORE=}" ;;
    LOGREG_SCORE=*) LOGREG_SCORE="${tok#LOGREG_SCORE=}" ;;
    REQUEST_MEMORY=*) REQUEST_MEMORY="${tok#REQUEST_MEMORY=}" ;;
    REQUEST_CPUS=*) REQUEST_CPUS="${tok#REQUEST_CPUS=}" ;;
    MAX_ROWS=*) MAX_ROWS="${tok#MAX_ROWS=}" ;;
    MAX_SHARDS=*) MAX_SHARDS="${tok#MAX_SHARDS=}" ;;
    ALGORITHMS=*) ALGORITHMS="${tok#ALGORITHMS=}" ;;
    VARIANTS=*) VARIANTS="${tok#VARIANTS=}" ;;
    SWEEP=*) SWEEP="${tok#SWEEP=}" ;;
    INCLUDE_ISOLATION_CONTEXT=*) INCLUDE_ISOLATION_CONTEXT="${tok#INCLUDE_ISOLATION_CONTEXT=}" ;;
    TOP_N=*) TOP_N="${tok#TOP_N=}" ;;
    LINEAR_BACKEND=*) LINEAR_BACKEND="${tok#LINEAR_BACKEND=}" ;;
    NN_HIDDEN=*) NN_HIDDEN="${tok#NN_HIDDEN=}" ;;
    NN_EPOCHS=*) NN_EPOCHS="${tok#NN_EPOCHS=}" ;;
    NN_PATIENCE=*) NN_PATIENCE="${tok#NN_PATIENCE=}" ;;
    NN_BATCH_SIZE=*) NN_BATCH_SIZE="${tok#NN_BATCH_SIZE=}" ;;
    NN_LEARNING_RATE=*) NN_LEARNING_RATE="${tok#NN_LEARNING_RATE=}" ;;
    NN_L2=*) NN_L2="${tok#NN_L2=}" ;;
    ALLOW_MISSING_FULL_FEATURES=*) ALLOW_MISSING_FULL_FEATURES="${tok#ALLOW_MISSING_FULL_FEATURES=}" ;;
    -h|--help|help) usage; exit 0 ;;
    *) die "Unknown argument: $tok" ;;
  esac
done

if [[ "$LOG_EXPLICIT" == "0" ]]; then
  LOG="${OUTDIR}/stacked_bdt_mlp_full_feature_sweep_${STAMP}.log"
fi

common_args=(
  --mlp-cache "$MLP_CACHE"
  --bdt-cache "$BDT_CACHE"
  --outdir "$OUTDIR"
  --mlp-score "$MLP_SCORE"
  --bdt-score "$BDT_SCORE"
  --max-rows "$MAX_ROWS"
  --max-shards "$MAX_SHARDS"
  --algorithms "$ALGORITHMS"
  --variants "$VARIANTS"
  --sweep "$SWEEP"
  --top-n "$TOP_N"
  --linear-backend "$LINEAR_BACKEND"
  --nn-hidden "$NN_HIDDEN"
  --nn-epochs "$NN_EPOCHS"
  --nn-patience "$NN_PATIENCE"
  --nn-batch-size "$NN_BATCH_SIZE"
  --nn-learning-rate "$NN_LEARNING_RATE"
  --nn-l2 "$NN_L2"
  --include-controls
  --note "Diagnostic full-feature BDT+MLP stacker sweep; runtime BDT score is an input."
)
if [[ -n "$LOGREG_CACHE" ]]; then
  common_args+=( --logreg-cache "$LOGREG_CACHE" --logreg-score "$LOGREG_SCORE" )
fi
if [[ "$ALLOW_MISSING_FULL_FEATURES" == "1" ]]; then
  common_args+=( --allow-missing-full-features )
fi
if [[ "$INCLUDE_ISOLATION_CONTEXT" == "1" ]]; then
  common_args+=( --include-isolation-context )
fi
PYTHON_RUN="$ML_PYTHON"
if [[ ! -x "$PYTHON_RUN" ]]; then
  PYTHON_RUN="${RJ_LOCAL_PYTHON:-python3}"
fi

print_plan() {
  say "RECOILJETS_AUAU_STACKED_BDT_MLP_FULL_FEATURE_SWEEP_SUBMIT_V1"
  say "host=$(hostname -f 2>/dev/null || hostname)"
  say "timestamp=$STAMP"
  say "repo=$RJ_REPO_BASE"
  say "source=$SOURCE"
  say "outdir=$OUTDIR"
  say "sub_root=$SUB_ROOT"
  say "log=$LOG"
  say "mlp_cache=$MLP_CACHE"
  say "bdt_cache=$BDT_CACHE"
  say "logreg_cache=${LOGREG_CACHE:-<none>}"
  say "mlp_score=$MLP_SCORE bdt_score=$BDT_SCORE logreg_score=$LOGREG_SCORE"
  say "request_memory=$REQUEST_MEMORY request_cpus=$REQUEST_CPUS"
  say "max_rows=$MAX_ROWS max_shards=$MAX_SHARDS algorithms=$ALGORITHMS variants=${VARIANTS:-<sweep-default>} sweep=$SWEEP top_n=$TOP_N"
  say "include_isolation_context=$INCLUDE_ISOLATION_CONTEXT"
  say "linear_backend=$LINEAR_BACKEND nn_hidden=$NN_HIDDEN nn_epochs=$NN_EPOCHS nn_patience=$NN_PATIENCE nn_batch_size=$NN_BATCH_SIZE nn_lr=$NN_LEARNING_RATE nn_l2=$NN_L2"
  say "allow_missing_full_features=$ALLOW_MISSING_FULL_FEATURES"
}

case "$mode" in
  selfTest|self-test)
    print_plan
    mkdir -p "$OUTDIR"
    "$PYTHON_RUN" "${RJ_REPO_BASE}/scripts/train_auau_stacked_bdt_mlp_sweep.py" \
      --outdir "$OUTDIR/self_test" \
      --self-test \
      --algorithms "$ALGORITHMS" \
      --variants "$VARIANTS" \
      --sweep "$SWEEP" \
      --linear-backend "$LINEAR_BACKEND" \
      --nn-hidden "$NN_HIDDEN" \
      --nn-epochs "$NN_EPOCHS" \
      --nn-patience "$NN_PATIENCE" \
      --nn-batch-size "$NN_BATCH_SIZE" \
      --nn-learning-rate "$NN_LEARNING_RATE" \
      --nn-l2 "$NN_L2" \
      --top-n "$TOP_N" \
      $([[ "$INCLUDE_ISOLATION_CONTEXT" == "1" ]] && printf '%s' "--include-isolation-context")
    ;;
  preflight)
    print_plan
    [[ -s "$MLP_CACHE" ]] || die "Missing MLP_CACHE=$MLP_CACHE"
    [[ -s "$BDT_CACHE" ]] || die "Missing BDT_CACHE=$BDT_CACHE"
    [[ -z "$LOGREG_CACHE" || -s "$LOGREG_CACHE" ]] || die "Missing LOGREG_CACHE=$LOGREG_CACHE"
    mkdir -p "$OUTDIR"
    "$PYTHON_RUN" "${RJ_REPO_BASE}/scripts/train_auau_stacked_bdt_mlp_sweep.py" \
      "${common_args[@]}" \
      --preflight-only
    ;;
  train)
    print_plan
    [[ -s "$MLP_CACHE" ]] || die "Missing MLP_CACHE=$MLP_CACHE"
    [[ -s "$BDT_CACHE" ]] || die "Missing BDT_CACHE=$BDT_CACHE"
    [[ -z "$LOGREG_CACHE" || -s "$LOGREG_CACHE" ]] || die "Missing LOGREG_CACHE=$LOGREG_CACHE"
    mkdir -p "$OUTDIR" "$SUB_ROOT" "$(dirname "$LOG")"
    worker="${SUB_ROOT}/run_stacked_bdt_mlp_full_feature_sweep.sh"
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
if [[ -d "\$MYINSTALL" ]]; then
  source /opt/sphenix/core/bin/setup_local.sh "\$MYINSTALL" || true
  echo "[auauStackSweep] worker setup_local rc=\$?" >&2
fi
if [[ -d "\$MYINSTALL_AUAU" ]]; then
  source /opt/sphenix/core/bin/setup_local.sh "\$MYINSTALL_AUAU" || true
  echo "[auauStackSweep] worker setup_local_auau rc=\$?" >&2
fi
set -u

ml_python="${ML_PYTHON}"
ml_python_prefix="\$(cd "\$(dirname "\$ml_python")/.." && pwd -P 2>/dev/null || true)"
ml_python_real="\$(readlink -f "\$ml_python" 2>/dev/null || true)"
ml_python_real_prefix=""
if [[ -n "\$ml_python_real" ]]; then
  ml_python_real_prefix="\$(cd "\$(dirname "\$ml_python_real")/.." && pwd -P 2>/dev/null || true)"
fi
ld_candidates=()
for d in "\$ml_python_prefix/lib" "\$ml_python_prefix/lib64" "\$ml_python_real_prefix/lib" "\$ml_python_real_prefix/lib64"; do
  [[ -d "\$d" ]] && ld_candidates+=("\$d")
done
if [[ \${#ld_candidates[@]} -gt 0 ]]; then
  ld_joined="\$(IFS=:; echo "\${ld_candidates[*]}")"
  export LD_LIBRARY_PATH="\$ld_joined:\${LD_LIBRARY_PATH:-}"
fi

extra_args=()
if [[ "$ALLOW_MISSING_FULL_FEATURES" == "1" ]]; then
  extra_args+=(--allow-missing-full-features)
fi
if [[ "$INCLUDE_ISOLATION_CONTEXT" == "1" ]]; then
  extra_args+=(--include-isolation-context)
fi
if [[ -n "$LOGREG_CACHE" ]]; then
  extra_args+=(--logreg-cache "$LOGREG_CACHE" --logreg-score "$LOGREG_SCORE")
fi
echo "RECOILJETS_AUAU_STACKED_BDT_MLP_FULL_FEATURE_SWEEP_WORKER_V1"
echo "host=\$(hostname -f 2>/dev/null || hostname)"
echo "start=\$(date)"
"\$ml_python" scripts/train_auau_stacked_bdt_mlp_sweep.py \\
  --mlp-cache "$MLP_CACHE" \\
  --bdt-cache "$BDT_CACHE" \\
  --outdir "$OUTDIR" \\
  --mlp-score "$MLP_SCORE" \\
  --bdt-score "$BDT_SCORE" \\
  --max-rows "$MAX_ROWS" \\
  --max-shards "$MAX_SHARDS" \\
  --algorithms "$ALGORITHMS" \\
  --variants "$VARIANTS" \\
  --sweep "$SWEEP" \\
  --top-n "$TOP_N" \\
  --linear-backend "$LINEAR_BACKEND" \\
  --nn-hidden "$NN_HIDDEN" \\
  --nn-epochs "$NN_EPOCHS" \\
  --nn-patience "$NN_PATIENCE" \\
  --nn-batch-size "$NN_BATCH_SIZE" \\
  --nn-learning-rate "$NN_LEARNING_RATE" \\
  --nn-l2 "$NN_L2" \\
  --include-controls \\
  "\${extra_args[@]}" \\
  --note "Diagnostic full-feature BDT+MLP stacker sweep; runtime BDT score is an input."
echo "finish=\$(date)"
EOF
    chmod +x "$worker"
    sub="${SUB_ROOT}/stacked_bdt_mlp_full_feature_sweep.sub"
    cat > "$sub" <<EOF
universe = vanilla
executable = ${worker}
output = ${SUB_ROOT}/stacked_bdt_mlp_full_feature_sweep.out
error = ${SUB_ROOT}/stacked_bdt_mlp_full_feature_sweep.err
log = ${SUB_ROOT}/stacked_bdt_mlp_full_feature_sweep.log
request_memory = ${REQUEST_MEMORY}
request_cpus = ${REQUEST_CPUS}
+JobBatchName = "auau_stack_bdt_mlp_sweep_${STAMP}"
notification = Never
queue
EOF
    say "worker=$worker"
    say "submit=$sub"
    if [[ "${RJ_DAG_DRYRUN:-0}" == "1" ]]; then
      say "RJ_DAG_DRYRUN=1; not submitting"
      exit 0
    fi
    command -v condor_submit >/dev/null 2>&1 || die "condor_submit not found in PATH"
    condor_submit "$sub" | tee "$LOG"
    say "rank_table=${OUTDIR}/stacked_sweep_rank_table.csv"
    say "top4=${OUTDIR}/stacked_sweep_top4.json"
    ;;
  status)
    print_plan
    for path in \
      "$OUTDIR/stacked_sweep_preflight.json" \
      "$OUTDIR/stacked_sweep_rank_table.csv" \
      "$OUTDIR/stacked_sweep_top4.json" \
      "$OUTDIR/stacked_sweep_metrics.json" \
      "$OUTDIR/stacked_sweep_training_history.csv"; do
      if [[ -s "$path" ]]; then
        say "exists: $path"
        [[ "$path" == *.csv ]] && head -5 "$path" || true
      else
        warn "missing: $path"
      fi
    done
    if [[ -d "$SUB_ROOT" ]]; then
      say "submit-root tail:"
      find "$SUB_ROOT" -maxdepth 1 -type f \( -name '*.out' -o -name '*.err' -o -name '*.log' \) -printf '%TY-%Tm-%Td %TH:%TM %p\n' | sort | tail -20 || true
    fi
    ;;
  -h|--help|help)
    usage
    ;;
  *)
    die "Unknown mode: $mode"
    ;;
esac
