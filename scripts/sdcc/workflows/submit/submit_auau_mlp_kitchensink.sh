#!/usr/bin/env bash
set -euo pipefail

RJ_REPO_BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
readonly RJ_REPO_BASE

ML_PYTHON="${RJ_ML_PYTHON:-/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python}"
SOURCE="${RJ_AUAU_MLP_KITCHENSINK_SOURCE:-}"
MODEL_BASE="${RJ_AUAU_MLP_MODEL_BASE:-${RJ_REPO_BASE}/mlp_models}"
STAMP="${RJ_AUAU_MLP_KITCHENSINK_STAMP:-$(date +%Y%m%d_%H%M%S)}"
MODEL_DIR="${RJ_AUAU_MLP_KITCHENSINK_MODEL_DIR:-${MODEL_BASE}/tight_mlp_kitchensink_${STAMP}}"
SUB_ROOT="${RJ_AUAU_MLP_KITCHENSINK_SUB_ROOT:-${RJ_REPO_BASE}/condor_sub/auauTightMLPKitchenSink_${STAMP}}"
REQUEST_MEMORY="${RJ_AUAU_MLP_KITCHENSINK_REQUEST_MEMORY:-48000MB}"

say() { printf '\033[1;36m[auauMLPKitchenSink]\033[0m %s\n' "$*"; }
err() { printf '\033[1;31m[auauMLPKitchenSink][ERR]\033[0m %s\n' "$*" >&2; }
die() { err "$*"; exit 2; }

usage() {
  cat <<'EOF'
Usage:
  ./scripts/submit_auau_mlp_kitchensink.sh SOURCE=/path/to/extraction

Submits one isolated Condor job for a high-capacity AuAu tight-MLP kitchen-sink
side test. It uses extended shower-shape / energy-ratio inputs, high-pT WP80
model selection, and BDT-guided hard-example weighting during training only.
EOF
}

for tok in "$@"; do
  case "$tok" in
    SOURCE=*) SOURCE="${tok#SOURCE=}" ;;
    MODEL_DIR=*) MODEL_DIR="${tok#MODEL_DIR=}" ;;
    SUB_ROOT=*) SUB_ROOT="${tok#SUB_ROOT=}" ;;
    REQUEST_MEMORY=*) REQUEST_MEMORY="${tok#REQUEST_MEMORY=}" ;;
    -h|--help|help) usage; exit 0 ;;
  esac
done

[[ -n "$SOURCE" && -d "$SOURCE" ]] || die "SOURCE=/path/to/extraction is required"
mkdir -p "$MODEL_DIR" "$SUB_ROOT"

say "RECOILJETS_AUAU_MLP_KITCHENSINK_SUBMIT_V1"
say "submit_host=$(hostname -f 2>/dev/null || hostname)"
say "timestamp=$STAMP"
say "source=$SOURCE"
say "model_dir=$MODEL_DIR"
say "sub_root=$SUB_ROOT"
say "request_memory=$REQUEST_MEMORY"
say "features=extended shower-shape + energy-ratio + width-ratio + centrality"
say "training=5:35 with pT-bin caps/weights and highpt_wp80 selection"
say "hard_examples=auau_tight_bdt_score is used only for training weights, not as an MLP input"

worker="${SUB_ROOT}/train_kitchensink_worker.sh"
cat > "$worker" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export RJ_ML_PYTHON="${ML_PYTHON}"
cd "${RJ_REPO_BASE}"
echo "[auauMLPKitchenSink] worker start host=\$(hostname -f 2>/dev/null || hostname) date=\$(date)"
./scripts/auau_tight_mlp_pipeline.sh trainKitchenSinkFromExtraction SOURCE="${SOURCE}" MODEL_DIR="${MODEL_DIR}"
./scripts/auau_tight_mlp_pipeline.sh applyCheck MODEL_DIR="${MODEL_DIR}"
echo "DONE_KITCHENSINK_MLP_MODEL_DIR=${MODEL_DIR}"
echo "[auauMLPKitchenSink] worker done date=\$(date)"
EOF
chmod +x "$worker"

sub="${SUB_ROOT}/kitchensink_train.sub"
cat > "$sub" <<EOF
universe = vanilla
executable = ${worker}
output = ${SUB_ROOT}/kitchensink_train.out
error = ${SUB_ROOT}/kitchensink_train.err
log = ${SUB_ROOT}/kitchensink_train.log
request_memory = ${REQUEST_MEMORY}
request_cpus = 1
+JobBatchName = "auau_mlp_kitchensink_${STAMP}"
notification = Never
queue
EOF

say "submit_file=$sub"
if [[ "${RJ_DAG_DRYRUN:-0}" == "1" ]]; then
  say "RJ_DAG_DRYRUN=1; not submitting"
  exit 0
fi

condor_submit "$sub"
