#!/usr/bin/env bash
set -euo pipefail

RJ_REPO_BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
readonly RJ_REPO_BASE

ML_PYTHON="${RJ_ML_PYTHON:-/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python}"
SOURCE="${RJ_AUAU_STACK_SOURCE:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049}"
MODEL_BASE="${RJ_AUAU_MLP_MODEL_BASE:-/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models}"
STAMP="${RJ_AUAU_STACK_STAMP:-$(date +%Y%m%d_%H%M%S)}"
OUTDIR="${RJ_AUAU_STACK_OUTDIR:-${MODEL_BASE}/stacked_bdt_mlp_calibrator_${STAMP}}"
SESSION="${RJ_AUAU_STACK_TMUX:-mlp_stacked_bdt_mlp_calib_${STAMP}}"
LOG="${RJ_AUAU_STACK_LOG:-${OUTDIR}/train_stacked_bdt_mlp_calibrator_${STAMP}.log}"
MLP_CACHE="${RJ_AUAU_STACK_MLP_CACHE:-${SOURCE}/reports/mlp_model_validation_condor_deep_primary_ratios_nostat_fullval_20260512_145449/score_caches.list}"
BDT_CACHE="${RJ_AUAU_STACK_BDT_CACHE:-${SOURCE}/reports/model_validation_condor_20260511_203441/score_caches.list}"
MLP_SCORE="${RJ_AUAU_STACK_MLP_SCORE:-score_centInputBase3x3WidthRatiosMLP_pt1535}"
BDT_SCORE="${RJ_AUAU_STACK_BDT_SCORE:-score_ptFine_cent7}"
PT_BINS="${RJ_AUAU_STACK_PT_BINS:-15,18,20,22.5,25,30,35}"
MAX_ROWS="${RJ_AUAU_STACK_MAX_ROWS:-0}"
MAX_SHARDS="${RJ_AUAU_STACK_MAX_SHARDS:-0}"
MODELS="${RJ_AUAU_STACK_MODELS:-logistic,gbm}"

say() { printf '\033[1;36m[auauStackBDTMLP]\033[0m %s\n' "$*"; }
err() { printf '\033[1;31m[auauStackBDTMLP][ERR]\033[0m %s\n' "$*" >&2; }
die() { err "$*"; exit 2; }

usage() {
  cat <<'EOF'
Usage:
  ./scripts/submit_auau_stacked_bdt_mlp_calibrator.sh [OUTDIR=/path] [MAX_ROWS=N]

Launches a tmux sidecar that trains diagnostic BDT+MLP stackers from existing
row-compatible validation score caches. This is a ceiling/complementarity test,
not a production runtime mode.

Useful knobs:
  MLP_CACHE=/path/score_caches.list
  BDT_CACHE=/path/score_caches.list
  MLP_SCORE=score_centInputBase3x3WidthRatiosMLP_pt1535
  BDT_SCORE=score_ptFine_cent7
  PT_BINS=15,18,20,22.5,25,30,35
  MAX_ROWS=0
  MAX_SHARDS=0
  MODELS=logistic,gbm
  RJ_DRYRUN=1
EOF
}

for tok in "$@"; do
  case "$tok" in
    SOURCE=*) SOURCE="${tok#SOURCE=}" ;;
    OUTDIR=*) OUTDIR="${tok#OUTDIR=}" ;;
    SESSION=*|TMUX=*) SESSION="${tok#*=}" ;;
    LOG=*) LOG="${tok#LOG=}" ;;
    MLP_CACHE=*) MLP_CACHE="${tok#MLP_CACHE=}" ;;
    BDT_CACHE=*) BDT_CACHE="${tok#BDT_CACHE=}" ;;
    MLP_SCORE=*) MLP_SCORE="${tok#MLP_SCORE=}" ;;
    BDT_SCORE=*) BDT_SCORE="${tok#BDT_SCORE=}" ;;
    PT_BINS=*) PT_BINS="${tok#PT_BINS=}" ;;
    MAX_ROWS=*) MAX_ROWS="${tok#MAX_ROWS=}" ;;
    MAX_SHARDS=*) MAX_SHARDS="${tok#MAX_SHARDS=}" ;;
    MODELS=*) MODELS="${tok#MODELS=}" ;;
    -h|--help|help) usage; exit 0 ;;
    *) die "Unknown argument: $tok" ;;
  esac
done

[[ -s "$MLP_CACHE" ]] || die "Missing MLP_CACHE=$MLP_CACHE"
[[ -s "$BDT_CACHE" ]] || die "Missing BDT_CACHE=$BDT_CACHE"
[[ "$MAX_ROWS" =~ ^[0-9]+$ ]] || die "MAX_ROWS must be a nonnegative integer"
[[ "$MAX_SHARDS" =~ ^[0-9]+$ ]] || die "MAX_SHARDS must be a nonnegative integer"

mkdir -p "$OUTDIR" "$(dirname "$LOG")"

say "RECOILJETS_AUAU_STACKED_BDT_MLP_SUBMIT_V1"
say "host=$(hostname -f 2>/dev/null || hostname)"
say "stamp=$STAMP"
say "session=$SESSION"
say "outdir=$OUTDIR"
say "log=$LOG"
say "mlp_cache=$MLP_CACHE"
say "bdt_cache=$BDT_CACHE"
say "mlp_score=$MLP_SCORE bdt_score=$BDT_SCORE"
say "pt_bins=$PT_BINS max_rows=$MAX_ROWS max_shards=$MAX_SHARDS models=$MODELS"

run_script="/tmp/${SESSION}.sh"
cat > "$run_script" <<EOF
#!/usr/bin/env bash
set -euo pipefail
cd "${RJ_REPO_BASE}"
export RJ_ML_PYTHON="${ML_PYTHON}"
echo "RECOILJETS_AUAU_STACKED_BDT_MLP_TMUX_V1"
echo "host=\$(hostname -f 2>/dev/null || hostname)"
echo "start=\$(date)"
echo "outdir=${OUTDIR}"
echo "log=${LOG}"
"${ML_PYTHON}" scripts/train_auau_stacked_bdt_mlp_calibrator.py \\
  --mlp-cache "${MLP_CACHE}" \\
  --bdt-cache "${BDT_CACHE}" \\
  --outdir "${OUTDIR}" \\
  --mlp-score "${MLP_SCORE}" \\
  --bdt-score "${BDT_SCORE}" \\
  --pt-bins "${PT_BINS}" \\
  --max-rows "${MAX_ROWS}" \\
  --max-shards "${MAX_SHARDS}" \\
  --models "${MODELS}" \\
  --note "Diagnostic BDT+primary-MLP stacker; not a production runtime mode."
echo "finish=\$(date)"
EOF
chmod +x "$run_script"

if [[ "${RJ_DRYRUN:-0}" == "1" ]]; then
  say "RJ_DRYRUN=1; not starting tmux"
  say "run_script=$run_script"
  exit 0
fi

if tmux has-session -t "$SESSION" 2>/dev/null; then
  die "tmux session already exists: $SESSION"
fi

tmux new-session -d -s "$SESSION" "bash '$run_script' 2>&1 | tee '$LOG'"
say "started tmux session=$SESSION"
say "tail with: tail -f $LOG"
say "attach with: tmux attach -t $SESSION"
