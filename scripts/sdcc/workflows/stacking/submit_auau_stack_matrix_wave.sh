#!/usr/bin/env bash
set -euo pipefail

RJ_REPO_BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
readonly RJ_REPO_BASE

ML_PYTHON="${RJ_ML_PYTHON:-/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python}"
SOURCE="${RJ_AUAU_STACK_SOURCE:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049}"
MODEL_BASE="${RJ_AUAU_MLP_MODEL_BASE:-/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models}"
STAMP="${RJ_AUAU_MATRIX_WAVE_STAMP:-$(date +%Y%m%d_%H%M%S)}"
WAVE="${RJ_AUAU_MATRIX_WAVE:-A}"
REQUEST_CPUS="${RJ_AUAU_MATRIX_WAVE_REQUEST_CPUS:-1}"
MAXJOBS="${RJ_AUAU_MATRIX_WAVE_MAXJOBS:-4}"
EXPECTED_SHARDS="${RJ_AUAU_MATRIX_WAVE_EXPECTED_SHARDS:-80}"
TOP_N="${RJ_AUAU_MATRIX_WAVE_TOP_N:-4}"
ALGORITHMS="${RJ_AUAU_MATRIX_WAVE_ALGORITHMS:-logistic,gbm,nn}"
LINEAR_BACKEND="${RJ_AUAU_MATRIX_WAVE_LINEAR_BACKEND:-numpy}"
NN_HIDDEN="${RJ_AUAU_MATRIX_WAVE_NN_HIDDEN:-96,48}"
NN_EPOCHS="${RJ_AUAU_MATRIX_WAVE_NN_EPOCHS:-160}"
NN_PATIENCE="${RJ_AUAU_MATRIX_WAVE_NN_PATIENCE:-24}"
NN_BATCH_SIZE="${RJ_AUAU_MATRIX_WAVE_NN_BATCH_SIZE:-8192}"
NN_LEARNING_RATE="${RJ_AUAU_MATRIX_WAVE_NN_LR:-0.0015}"
NN_L2="${RJ_AUAU_MATRIX_WAVE_NN_L2:-0.001}"
MATRIX_ROUTINGS="${RJ_AUAU_MATRIX_WAVE_ROUTINGS:-all}"

MLP_CACHE="${RJ_AUAU_STACK_MLP_CACHE:-${SOURCE}/reports/mlp_model_validation_condor_stack_full_features_uncapped_20260513_123934/score_caches.list}"
BDT_CACHE="${RJ_AUAU_STACK_BDT_CACHE:-${SOURCE}/reports/model_validation_condor_stack_full_features_uncapped_20260513_123934/score_caches.list}"
LOGREG_CACHE="${RJ_AUAU_STACK_LOGREG_CACHE:-${SOURCE}/reports/logreg_model_validation_condor_20260513_213044/score_caches.list}"
MLP_SCORE="${RJ_AUAU_STACK_MLP_SCORE:-score_centInputBase3x3WidthRatiosMLP_pt1535}"
BDT_SCORE="${RJ_AUAU_STACK_BDT_SCORE:-score_ptFine_cent7}"
LOGREG_SCORE="${RJ_AUAU_STACK_LOGREG_SCORE:-score_centEtFullLogReg_pt1535}"

OUTDIR="${RJ_AUAU_MATRIX_WAVE_OUTDIR:-${MODEL_BASE}/stack_matrix_wave_${WAVE}_${STAMP}}"
SUB_ROOT="${RJ_AUAU_MATRIX_WAVE_SUB_ROOT:-${RJ_REPO_BASE}/condor_sub/auauStackMatrixWave_${WAVE}_${STAMP}}"
LOG="${RJ_AUAU_MATRIX_WAVE_LOG:-${OUTDIR}/stack_matrix_wave_${WAVE}_${STAMP}.log}"
OUTDIR_EXPLICIT=0
SUB_ROOT_EXPLICIT=0
LOG_EXPLICIT=0
[[ -n "${RJ_AUAU_MATRIX_WAVE_OUTDIR+x}" ]] && OUTDIR_EXPLICIT=1
[[ -n "${RJ_AUAU_MATRIX_WAVE_SUB_ROOT+x}" ]] && SUB_ROOT_EXPLICIT=1
[[ -n "${RJ_AUAU_MATRIX_WAVE_LOG+x}" ]] && LOG_EXPLICIT=1

say() { printf '\033[1;36m[auauStackMatrixWave]\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[auauStackMatrixWave][WARN]\033[0m %s\n' "$*" >&2; }
err() { printf '\033[1;31m[auauStackMatrixWave][ERR]\033[0m %s\n' "$*" >&2; }
die() { err "$*"; exit 2; }

usage() {
  cat <<'EOF'
Usage:
  ./scripts/submit_auau_stack_matrix_wave.sh selfTest
  ./scripts/submit_auau_stack_matrix_wave.sh preflight WAVE=A
  ./scripts/submit_auau_stack_matrix_wave.sh train WAVE=A
  ./scripts/submit_auau_stack_matrix_wave.sh status OUTDIR=/gpfs/...

Waves:
  A       15-35 pair score-only
  B       15-35 tri-score-only
  C       15-35 pair/tri context
  D       15-35 pair/tri full and one-score rescue controls
  E/F/G/H Same as A/B/C/D for 5-35
  auto15  15-35 A/B/C/D cohort workers with queue-gated concurrency
  auto35  5-35 E/F/G/H cohort workers with queue-gated concurrency

All matrix waves are validation/ranking only. They never submit production.
EOF
}

mode="${1:-train}"
if [[ $# -gt 0 ]]; then shift; fi

for tok in "$@"; do
  case "$tok" in
    WAVE=*) WAVE="${tok#WAVE=}" ;;
    OUTDIR=*) OUTDIR="${tok#OUTDIR=}"; OUTDIR_EXPLICIT=1 ;;
    SUB_ROOT=*) SUB_ROOT="${tok#SUB_ROOT=}"; SUB_ROOT_EXPLICIT=1 ;;
    LOG=*) LOG="${tok#LOG=}"; LOG_EXPLICIT=1 ;;
    MLP_CACHE=*) MLP_CACHE="${tok#MLP_CACHE=}" ;;
    BDT_CACHE=*) BDT_CACHE="${tok#BDT_CACHE=}" ;;
    LOGREG_CACHE=*) LOGREG_CACHE="${tok#LOGREG_CACHE=}" ;;
    MLP_SCORE=*) MLP_SCORE="${tok#MLP_SCORE=}" ;;
    BDT_SCORE=*) BDT_SCORE="${tok#BDT_SCORE=}" ;;
    LOGREG_SCORE=*) LOGREG_SCORE="${tok#LOGREG_SCORE=}" ;;
    ALGORITHMS=*) ALGORITHMS="${tok#ALGORITHMS=}" ;;
    MATRIX_ROUTINGS=*) MATRIX_ROUTINGS="${tok#MATRIX_ROUTINGS=}" ;;
    MAXJOBS=*) MAXJOBS="${tok#MAXJOBS=}" ;;
    REQUEST_CPUS=*) REQUEST_CPUS="${tok#REQUEST_CPUS=}" ;;
    EXPECTED_SHARDS=*) EXPECTED_SHARDS="${tok#EXPECTED_SHARDS=}" ;;
    -h|--help|help) usage; exit 0 ;;
    *) die "Unknown argument: $tok" ;;
  esac
done

if [[ "$OUTDIR_EXPLICIT" == "0" ]]; then
  OUTDIR="${MODEL_BASE}/stack_matrix_wave_${WAVE}_${STAMP}"
fi
if [[ "$SUB_ROOT_EXPLICIT" == "0" ]]; then
  SUB_ROOT="${RJ_REPO_BASE}/condor_sub/auauStackMatrixWave_${WAVE}_${STAMP}"
fi
if [[ "$LOG_EXPLICIT" == "0" ]]; then
  LOG="${OUTDIR}/stack_matrix_wave_${WAVE}_${STAMP}.log"
fi

wave_specs() {
  case "$WAVE" in
    A) printf '%s\n' "15-35 stack_pair_score_only 24000MB" ;;
    B) printf '%s\n' "15-35 stack_tri_score_only 24000MB" ;;
    C) printf '%s\n' "15-35 stack_pair_context 24000MB" "15-35 stack_tri_context 24000MB" ;;
    D) printf '%s\n' "15-35 stack_pair_full_no_iso 32000MB" "15-35 stack_tri_full_no_iso 32000MB" "15-35 stack_bdt_score_full_no_iso 32000MB" "15-35 stack_mlp_score_full_no_iso 32000MB" ;;
    E) printf '%s\n' "5-35 stack_pair_score_only 24000MB" ;;
    F) printf '%s\n' "5-35 stack_tri_score_only 24000MB" ;;
    G) printf '%s\n' "5-35 stack_pair_context 24000MB" "5-35 stack_tri_context 24000MB" ;;
    H) printf '%s\n' "5-35 stack_pair_full_no_iso 32000MB" "5-35 stack_tri_full_no_iso 32000MB" "5-35 stack_bdt_score_full_no_iso 32000MB" "5-35 stack_mlp_score_full_no_iso 32000MB" ;;
    auto15) printf '%s\n' "15-35 stack_pair_score_only 24000MB" "15-35 stack_tri_score_only 24000MB" "15-35 stack_pair_context 24000MB" "15-35 stack_tri_context 24000MB" "15-35 stack_pair_full_no_iso 32000MB" "15-35 stack_tri_full_no_iso 32000MB" "15-35 stack_bdt_score_full_no_iso 32000MB" "15-35 stack_mlp_score_full_no_iso 32000MB" ;;
    auto35) printf '%s\n' "5-35 stack_pair_score_only 24000MB" "5-35 stack_tri_score_only 24000MB" "5-35 stack_pair_context 24000MB" "5-35 stack_tri_context 24000MB" "5-35 stack_pair_full_no_iso 32000MB" "5-35 stack_tri_full_no_iso 32000MB" "5-35 stack_bdt_score_full_no_iso 32000MB" "5-35 stack_mlp_score_full_no_iso 32000MB" ;;
    *) die "Unknown WAVE=$WAVE" ;;
  esac
}

print_plan() {
  say "RECOILJETS_AUAU_STACK_MATRIX_WAVE_V1"
  say "host=$(hostname -f 2>/dev/null || hostname)"
  say "stamp=$STAMP"
  say "wave=$WAVE"
  say "repo=$RJ_REPO_BASE"
  say "source=$SOURCE"
  say "outdir=$OUTDIR"
  say "sub_root=$SUB_ROOT"
  say "log=$LOG"
  say "mlp_cache=$MLP_CACHE"
  say "bdt_cache=$BDT_CACHE"
  say "logreg_cache=$LOGREG_CACHE"
  say "algorithms=$ALGORITHMS routings=$MATRIX_ROUTINGS maxjobs=$MAXJOBS expected_shards=$EXPECTED_SHARDS"
  wave_specs | sed 's/^/[auauStackMatrixWave] worker_spec: /'
}

common_python_args() {
  local range="$1"
  local cohort="$2"
  local out="$3"
  printf '%q ' \
    scripts/train_auau_stacked_bdt_mlp_sweep.py \
    --mlp-cache "$MLP_CACHE" \
    --bdt-cache "$BDT_CACHE" \
    --logreg-cache "$LOGREG_CACHE" \
    --outdir "$out" \
    --mlp-score "$MLP_SCORE" \
    --bdt-score "$BDT_SCORE" \
    --logreg-score "$LOGREG_SCORE" \
    --algorithms "$ALGORITHMS" \
    --matrix-cohort "$cohort" \
    --matrix-wave-name "$WAVE" \
    --training-range "$range" \
    --matrix-routings "$MATRIX_ROUTINGS" \
    --require-full-stat \
    --expected-shards "$EXPECTED_SHARDS" \
    --stack-training-safety disjoint_base_scores_heldout_test \
    --linear-backend "$LINEAR_BACKEND" \
    --nn-hidden "$NN_HIDDEN" \
    --nn-epochs "$NN_EPOCHS" \
    --nn-patience "$NN_PATIENCE" \
    --nn-batch-size "$NN_BATCH_SIZE" \
    --nn-learning-rate "$NN_LEARNING_RATE" \
    --nn-l2 "$NN_L2" \
    --top-n "$TOP_N" \
    --note "Autonomous full-stat matrix wave ${WAVE}; stack training uses first-stage scores from disjoint validation-cache provenance plus a held-out test split."
}

write_worker() {
  local idx="$1"
  local range="$2"
  local cohort="$3"
  local mem="$4"
  local node_out="$OUTDIR/${range//-/_}/${cohort}"
  local worker="$SUB_ROOT/run_${idx}_${range//-/_}_${cohort}.sh"
  local args
  args="$(common_python_args "$range" "$cohort" "$node_out")"
  cat > "$worker" <<EOF
#!/usr/bin/env bash
set -euo pipefail
cd "$RJ_REPO_BASE"
export RJ_ML_PYTHON="$ML_PYTHON"
set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n
set -u
ml_python="$ML_PYTHON"
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
  export LD_LIBRARY_PATH="\$(IFS=:; echo "\${ld_candidates[*]}"):\${LD_LIBRARY_PATH:-}"
fi
echo "RECOILJETS_AUAU_STACK_MATRIX_WAVE_WORKER_V1"
echo "host=\$(hostname -f 2>/dev/null || hostname)"
echo "start=\$(date)"
echo "range=$range cohort=$cohort outdir=$node_out"
"\$ml_python" $args
echo "finish=\$(date)"
EOF
  chmod +x "$worker"
  printf '%s|%s|%s|%s|%s\n' "$idx" "$range" "$cohort" "$mem" "$worker"
}

case "$mode" in
  selfTest|self-test)
    print_plan
    tmp="${OUTDIR}/self_test"
    mkdir -p "$tmp"
    read -r test_range test_cohort _ < <(wave_specs | head -1)
    "$ML_PYTHON" "$RJ_REPO_BASE/scripts/train_auau_stacked_bdt_mlp_sweep.py" \
      --outdir "$tmp" \
      --self-test \
      --matrix-cohort "$test_cohort" \
      --training-range "$test_range" \
      --matrix-routings global_centEtInput,etFine_cent7_routed \
      --algorithms logistic,nn \
      --expected-shards 1 \
      --stack-training-safety disjoint_base_scores_heldout_test
    ;;
  preflight)
    print_plan
    [[ -s "$MLP_CACHE" ]] || die "Missing MLP_CACHE=$MLP_CACHE"
    [[ -s "$BDT_CACHE" ]] || die "Missing BDT_CACHE=$BDT_CACHE"
    [[ -s "$LOGREG_CACHE" ]] || die "Missing LOGREG_CACHE=$LOGREG_CACHE"
    ;;
  train)
    print_plan
    [[ -s "$MLP_CACHE" ]] || die "Missing MLP_CACHE=$MLP_CACHE"
    [[ -s "$BDT_CACHE" ]] || die "Missing BDT_CACHE=$BDT_CACHE"
    [[ -s "$LOGREG_CACHE" ]] || die "Missing LOGREG_CACHE=$LOGREG_CACHE"
    mkdir -p "$OUTDIR" "$SUB_ROOT" "$(dirname "$LOG")"
    worker_manifest="$SUB_ROOT/matrix_wave_workers.list"
    : > "$worker_manifest"
    idx=0
    while read -r range cohort mem; do
      [[ -n "$range" ]] || continue
      idx=$((idx + 1))
      write_worker "$idx" "$range" "$cohort" "$mem" >> "$worker_manifest"
    done < <(wave_specs)
    dag="$SUB_ROOT/stack_matrix_wave.dag"
    : > "$dag"
    while IFS='|' read -r node_id range cohort mem worker; do
      node="N${node_id}_${range//-/_}_${cohort}"
      sub="$SUB_ROOT/${node}.sub"
      cat > "$sub" <<EOF
universe = vanilla
executable = $worker
output = $SUB_ROOT/${node}.out
error = $SUB_ROOT/${node}.err
log = $SUB_ROOT/${node}.condor.log
request_memory = $mem
request_cpus = $REQUEST_CPUS
+JobFlavour = "tomorrow"
queue
EOF
      {
        echo "JOB $node $sub"
        echo "CATEGORY $node matrix"
      } >> "$dag"
    done < "$worker_manifest"
    echo "MAXJOBS matrix $MAXJOBS" >> "$dag"
    {
      print_plan
      say "worker_manifest=$worker_manifest"
      say "dag=$dag"
      condor_submit_dag "$dag"
    } | tee "$LOG"
    ;;
  status)
    print_plan
    find "$OUTDIR" -maxdepth 4 \( -name 'stacked_sweep_rank_table.csv' -o -name 'stacked_sweep_metrics.json' -o -name 'matrix_wave_manifest.json' -o -name 'stacked_sweep_top4.json' \) -print | sort
    ;;
  *)
    usage
    die "Unknown mode=$mode"
    ;;
esac
