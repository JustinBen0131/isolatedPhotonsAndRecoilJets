#!/usr/bin/env bash
set -euo pipefail

RJ_REPO_BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
readonly RJ_REPO_BASE

ML_PYTHON="${RJ_ML_PYTHON:-/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python}"
SOURCE="${RJ_AUAU_STACK_SOURCE:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049}"
MODEL_BASE="${RJ_AUAU_MLP_MODEL_BASE:-/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models}"
STAMP="${RJ_AUAU_STACK_CHAIN_STAMP:-$(date +%Y%m%d_%H%M%S)}"
CHAIN_ROOT="${RJ_AUAU_STACK_CHAIN_ROOT:-${RJ_REPO_BASE}/condor_sub/auauStackedBDTMLPFullFeatureChain_${STAMP}}"
SESSION="${RJ_AUAU_STACK_CHAIN_TMUX:-stack_full_feature_chain_${STAMP}}"
LOG="${RJ_AUAU_STACK_CHAIN_LOG:-${CHAIN_ROOT}/stack_full_feature_chain_${STAMP}.log}"
MLP_MODEL_DIR="${RJ_AUAU_STACK_MLP_MODEL_DIR:-/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/tight_mlp_deep_primary_ratios_nostat_20260512_005242}"
BDT_MODEL_DIR="${RJ_AUAU_STACK_BDT_MODEL_DIR:-/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_etfine_centstudy_current}"
BDT_MODEL_REGISTRY="${RJ_AUAU_STACK_BDT_MODEL_REGISTRY:-${BDT_MODEL_DIR}/model_registry.json}"
MLP_REPORT="${RJ_AUAU_STACK_MLP_REPORT:-${SOURCE}/reports/mlp_model_validation_condor_stack_full_features_${STAMP}}"
BDT_REPORT="${RJ_AUAU_STACK_BDT_REPORT:-${SOURCE}/reports/model_validation_condor_stack_full_features_${STAMP}}"
STACK_OUTDIR="${RJ_AUAU_STACK_SWEEP_OUTDIR:-${MODEL_BASE}/stacked_bdt_mlp_full_feature_sweep_${STAMP}}"
STACK_SUB_ROOT="${RJ_AUAU_STACK_SWEEP_SUB_ROOT:-${RJ_REPO_BASE}/condor_sub/auauStackedBDTMLPFullFeature_${STAMP}}"
GROUP_SIZE="${RJ_AUAU_STACK_VALIDATE_GROUP_SIZE:-100}"
SCORE_MAX_ROWS="${RJ_AUAU_STACK_VALIDATE_SCORE_MAX_ROWS:-400000}"
VALIDATE_MEMORY="${RJ_AUAU_STACK_VALIDATE_MEMORY:-3000MB}"
STACK_MEMORY="${RJ_AUAU_STACK_SWEEP_REQUEST_MEMORY:-12000MB}"
POLL_SECONDS="${RJ_AUAU_STACK_CHAIN_POLL_SECONDS:-300}"

say() { printf '\033[1;36m[auauStackChain]\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[auauStackChain][WARN]\033[0m %s\n' "$*" >&2; }
die() { printf '\033[1;31m[auauStackChain][ERR]\033[0m %s\n' "$*" >&2; exit 2; }

usage() {
  cat <<'EOF'
Usage:
  ./scripts/submit_auau_stacked_bdt_mlp_full_feature_chain.sh start
  ./scripts/submit_auau_stacked_bdt_mlp_full_feature_chain.sh status [CHAIN_ROOT=/path]

Starts one tmux watcher that:
  1. submits fresh BDT and MLP validation-cache DAGs with full-feature caches;
  2. waits until both reports are READY and score_caches.list exists;
  3. submits the full-feature stacked BDT+MLP sweep Condor job.

Mutating start mode requires RJ_DO_RUN=1.  The stackers are diagnostic MC
ceiling tests because the runtime BDT score is an input.
EOF
}

mode="${1:-start}"
if [[ $# -gt 0 ]]; then shift; fi
for tok in "$@"; do
  case "$tok" in
    SOURCE=*) SOURCE="${tok#SOURCE=}" ;;
    CHAIN_ROOT=*) CHAIN_ROOT="${tok#CHAIN_ROOT=}" ;;
    SESSION=*|TMUX=*) SESSION="${tok#*=}" ;;
    LOG=*) LOG="${tok#LOG=}" ;;
    MLP_MODEL_DIR=*) MLP_MODEL_DIR="${tok#MLP_MODEL_DIR=}" ;;
    BDT_MODEL_DIR=*) BDT_MODEL_DIR="${tok#BDT_MODEL_DIR=}" ;;
    BDT_MODEL_REGISTRY=*) BDT_MODEL_REGISTRY="${tok#BDT_MODEL_REGISTRY=}" ;;
    MLP_REPORT=*) MLP_REPORT="${tok#MLP_REPORT=}" ;;
    BDT_REPORT=*) BDT_REPORT="${tok#BDT_REPORT=}" ;;
    STACK_OUTDIR=*) STACK_OUTDIR="${tok#STACK_OUTDIR=}" ;;
    STACK_SUB_ROOT=*) STACK_SUB_ROOT="${tok#STACK_SUB_ROOT=}" ;;
    GROUP_SIZE=*) GROUP_SIZE="${tok#GROUP_SIZE=}" ;;
    SCORE_MAX_ROWS=*) SCORE_MAX_ROWS="${tok#SCORE_MAX_ROWS=}" ;;
    VALIDATE_MEMORY=*) VALIDATE_MEMORY="${tok#VALIDATE_MEMORY=}" ;;
    STACK_MEMORY=*) STACK_MEMORY="${tok#STACK_MEMORY=}" ;;
    POLL_SECONDS=*) POLL_SECONDS="${tok#POLL_SECONDS=}" ;;
    -h|--help|help) usage; exit 0 ;;
    *) die "Unknown argument: $tok" ;;
  esac
done

print_plan() {
  say "RECOILJETS_AUAU_STACKED_BDT_MLP_FULL_FEATURE_CHAIN_V1"
  say "host=$(hostname -f 2>/dev/null || hostname)"
  say "stamp=$STAMP"
  say "repo=$RJ_REPO_BASE"
  say "source=$SOURCE"
  say "chain_root=$CHAIN_ROOT"
  say "session=$SESSION"
  say "log=$LOG"
  say "mlp_model_dir=$MLP_MODEL_DIR"
  say "bdt_model_dir=$BDT_MODEL_DIR"
  say "bdt_model_registry=$BDT_MODEL_REGISTRY"
  say "mlp_report=$MLP_REPORT"
  say "bdt_report=$BDT_REPORT"
  say "stack_outdir=$STACK_OUTDIR"
  say "stack_sub_root=$STACK_SUB_ROOT"
  say "group_size=$GROUP_SIZE score_max_rows=$SCORE_MAX_ROWS validate_memory=$VALIDATE_MEMORY stack_memory=$STACK_MEMORY poll_seconds=$POLL_SECONDS"
}

case "$mode" in
  start)
    print_plan
    [[ "${RJ_DO_RUN:-0}" == "1" ]] || die "Set RJ_DO_RUN=1 to start the chain tmux and submit Condor validation jobs."
    [[ -d "$SOURCE" ]] || die "Missing SOURCE=$SOURCE"
    [[ -d "$MLP_MODEL_DIR" ]] || die "Missing MLP_MODEL_DIR=$MLP_MODEL_DIR"
    [[ -d "$BDT_MODEL_DIR" ]] || die "Missing BDT_MODEL_DIR=$BDT_MODEL_DIR"
    [[ -s "$BDT_MODEL_REGISTRY" ]] || die "Missing BDT_MODEL_REGISTRY=$BDT_MODEL_REGISTRY"
    command -v tmux >/dev/null 2>&1 || die "tmux not found"
    mkdir -p "$CHAIN_ROOT" "$(dirname "$LOG")"
    watcher="${CHAIN_ROOT}/run_stack_full_feature_chain.sh"
    cat > "$watcher" <<EOF
#!/usr/bin/env bash
set -euo pipefail
cd "$RJ_REPO_BASE"
export RJ_ML_PYTHON="$ML_PYTHON"
export RJ_NOTIFY_EMAILS="\${RJ_NOTIFY_EMAILS:-just0131@gmail.com}"
echo "RECOILJETS_AUAU_STACKED_BDT_MLP_FULL_FEATURE_CHAIN_WATCHER_V1"
echo "host=\$(hostname -f 2>/dev/null || hostname)"
echo "start=\$(date)"
echo "mlp_report=$MLP_REPORT"
echo "bdt_report=$BDT_REPORT"
echo "stack_outdir=$STACK_OUTDIR"

ready_summary_file() {
  local report="\$1"
  local candidate
  for candidate in "\$report/summary.txt" "\$report/validation_summary.txt"; do
    [[ -s "\$candidate" ]] && { printf '%s\n' "\$candidate"; return 0; }
  done
  return 1
}

ready_report() {
  local report="\$1"
  local summary
  [[ -s "\$report/score_caches.list" ]] || return 1
  summary="\$(ready_summary_file "\$report")" || return 1
  grep -q '^status=READY' "\$summary"
}

if ! ready_report "$MLP_REPORT"; then
  echo "[stackChain] submitting MLP full-feature validation cache DAG"
  RJ_AUAU_TIGHT_MLP_VALIDATE_STAMP="mlp_stack_full_features_${STAMP}" \\
  RJ_AUAU_TIGHT_MLP_VALIDATE_REQUEST_MEMORY="$VALIDATE_MEMORY" \\
  ./scripts/auau_tight_mlp_pipeline.sh validateOnSimCondor \\
    SOURCE="$SOURCE" \\
    MODEL_DIR="$MLP_MODEL_DIR" \\
    OUTDIR="$MLP_REPORT" \\
    groupSize "$GROUP_SIZE" \\
    scoreMaxRows "$SCORE_MAX_ROWS"
else
  echo "[stackChain] MLP report already READY"
fi

if ! ready_report "$BDT_REPORT"; then
  echo "[stackChain] submitting BDT full-feature validation cache DAG"
  RJ_AUAU_TIGHT_BDT_VALIDATE_STAMP="bdt_stack_full_features_${STAMP}" \\
  RJ_AUAU_TIGHT_BDT_VALIDATE_REQUEST_MEMORY="$VALIDATE_MEMORY" \\
  ./scripts/auau_tight_bdt_pipeline.sh validateOnSimCondor \\
    SOURCE="$SOURCE" \\
    MODEL_DIR="$BDT_MODEL_DIR" \\
    MODEL_REGISTRY="$BDT_MODEL_REGISTRY" \\
    OUTDIR="$BDT_REPORT" \\
    groupSize "$GROUP_SIZE" \\
    scoreMaxRows "$SCORE_MAX_ROWS"
else
  echo "[stackChain] BDT report already READY"
fi

while true; do
  echo "[stackChain] poll=\$(date)"
  if condor_q -hold "\${USER:-patsfan753}" 2>/dev/null | grep -q 'patsfan753'; then
    echo "[stackChain][WARN] held jobs detected for \${USER:-patsfan753}; inspect condor_q -hold" >&2
  fi
  if ready_report "$MLP_REPORT" && ready_report "$BDT_REPORT"; then
    echo "[stackChain] both validation reports READY"
    break
  fi
  echo "[stackChain] waiting for reports: mlp_ready=\$(ready_report "$MLP_REPORT" && echo 1 || echo 0) bdt_ready=\$(ready_report "$BDT_REPORT" && echo 1 || echo 0)"
  sleep "$POLL_SECONDS"
done

./scripts/submit_auau_stacked_bdt_mlp_sweep.sh preflight \\
  MLP_CACHE="$MLP_REPORT/score_caches.list" \\
  BDT_CACHE="$BDT_REPORT/score_caches.list" \\
  OUTDIR="$STACK_OUTDIR" \\
  MAX_SHARDS=1 \\
  MAX_ROWS=2000 \\
  ALGORITHMS=logistic \\
  SWEEP=compact

./scripts/submit_auau_stacked_bdt_mlp_sweep.sh train \\
  MLP_CACHE="$MLP_REPORT/score_caches.list" \\
  BDT_CACHE="$BDT_REPORT/score_caches.list" \\
  OUTDIR="$STACK_OUTDIR" \\
  SUB_ROOT="$STACK_SUB_ROOT" \\
  REQUEST_MEMORY="$STACK_MEMORY" \\
  ALGORITHMS=logistic,gbm \\
  SWEEP=full \\
  TOP_N=4

echo "DONE_STACK_FULL_FEATURE_CHAIN_ROOT=$CHAIN_ROOT"
echo "DONE_STACK_FULL_FEATURE_MLP_CACHE=$MLP_REPORT/score_caches.list"
echo "DONE_STACK_FULL_FEATURE_BDT_CACHE=$BDT_REPORT/score_caches.list"
echo "DONE_STACK_FULL_FEATURE_STACK_OUTDIR=$STACK_OUTDIR"
echo "finish=\$(date)"
EOF
    chmod +x "$watcher"
    if tmux has-session -t "$SESSION" 2>/dev/null; then
      die "tmux session already exists: $SESSION"
    fi
    tmux new-session -d -s "$SESSION" "bash '$watcher' 2>&1 | tee '$LOG'"
    say "started tmux session=$SESSION"
    say "tail with: tail -f $LOG"
    ;;
  status)
    print_plan
    [[ -s "$LOG" ]] && { say "log_tail=$LOG"; tail -80 "$LOG"; } || warn "missing log: $LOG"
    for path in \
      "$MLP_REPORT/summary.txt" "$MLP_REPORT/validation_summary.txt" \
      "$BDT_REPORT/summary.txt" "$BDT_REPORT/validation_summary.txt" \
      "$STACK_OUTDIR/stacked_sweep_rank_table.csv" "$STACK_OUTDIR/stacked_sweep_top4.json"; do
      [[ -s "$path" ]] && say "exists: $path" || warn "missing: $path"
    done
    ;;
  -h|--help|help)
    usage
    ;;
  *)
    die "Unknown mode: $mode"
    ;;
esac
