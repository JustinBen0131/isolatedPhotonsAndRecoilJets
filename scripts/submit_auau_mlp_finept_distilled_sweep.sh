#!/usr/bin/env bash
set -euo pipefail

RJ_REPO_BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
readonly RJ_REPO_BASE

ML_PYTHON="${RJ_ML_PYTHON:-/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python}"
SOURCE="${RJ_AUAU_MLP_FINEPT_SOURCE:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049}"
MODEL_BASE="${RJ_AUAU_MLP_MODEL_BASE:-/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models}"
STAMP="${RJ_AUAU_MLP_FINEPT_STAMP:-$(date +%Y%m%d_%H%M%S)}"
SWEEP_DIR="${RJ_AUAU_MLP_FINEPT_SWEEP_DIR:-${MODEL_BASE}/tight_mlp_finept_distilled_kitchen_v2_${STAMP}}"
SUB_ROOT="${RJ_AUAU_MLP_FINEPT_SUB_ROOT:-${RJ_REPO_BASE}/condor_sub/auauTightMLPFinePtDistilled_${STAMP}}"
REQUEST_MEMORY="${RJ_AUAU_MLP_FINEPT_REQUEST_MEMORY:-24000MB}"
MAXJOBS="${RJ_AUAU_MLP_FINEPT_MAXJOBS:-2}"
TRAIN_EPOCHS="${RJ_AUAU_MLP_FINEPT_EPOCHS:-220}"
TRAIN_PATIENCE="${RJ_AUAU_MLP_FINEPT_PATIENCE:-45}"
TRAIN_RESTARTS="${RJ_AUAU_MLP_FINEPT_RESTARTS:-1}"
TRAIN_HIDDEN="${RJ_AUAU_MLP_FINEPT_HIDDEN_LAYERS:-128,64,32}"
TRAIN_GRID="${RJ_AUAU_MLP_FINEPT_HIDDEN_LAYER_GRID:-160,80,40}"
TRAIN_BATCH="${RJ_AUAU_MLP_FINEPT_BATCH_SIZE:-4096}"
TRAIN_LR="${RJ_AUAU_MLP_FINEPT_LEARNING_RATE:-4.0e-4}"
TRAIN_L2="${RJ_AUAU_MLP_FINEPT_L2:-8.0e-5}"
TRAIN_MIN_DELTA="${RJ_AUAU_MLP_FINEPT_MIN_DELTA:-2.0e-6}"
ROUTE_CAP="${RJ_AUAU_MLP_FINEPT_MAX_ROWS_PER_PT_BIN_CLASS:-160000}"
DISTILLATION_STRENGTH="${RJ_AUAU_MLP_FINEPT_DISTILLATION_STRENGTH:-0.18}"
CONDITIONAL_JITTER="${RJ_AUAU_MLP_FINEPT_CONDITIONAL_JITTER:-0.008}"
VALIDATION_CACHE="${RJ_AUAU_MLP_FINEPT_VALIDATION_CACHE:-${SOURCE}/reports/mlp_model_validation_condor_deep_primary_ratios_nostat_fullval_20260512_145449/score_caches.list}"
VALIDATION_OUTDIR="${RJ_AUAU_MLP_FINEPT_VALIDATION_OUTDIR:-}"
VALIDATION_PT_BINS="${RJ_AUAU_MLP_FINEPT_VALIDATION_PT_BINS:-15,18,20,22.5,25,30,35}"

say() { printf '\033[1;36m[auauMLPFinePtDistilled]\033[0m %s\n' "$*"; }
err() { printf '\033[1;31m[auauMLPFinePtDistilled][ERR]\033[0m %s\n' "$*" >&2; }
die() { err "$*"; exit 2; }

usage() {
  cat <<'EOF'
Usage:
  ./scripts/submit_auau_mlp_finept_distilled_sweep.sh [SOURCE=/path] [VALIDATION_CACHE=/path/score_caches.list]

Submits an isolated Condor DAG for a fine-pT routed ABCD-safe distilled kitchen
AuAu tight-MLP validation candidate:
  15-18, 18-20, 20-22.5, 22.5-25, 25-30, 30-35 GeV.

Useful knobs:
  SWEEP_DIR=/path
  REQUEST_MEMORY=24000MB
  MAXJOBS=2
  EPOCHS=220
  PATIENCE=45
  RESTARTS=1
  DISTILLATION_STRENGTH=0.18
  VALIDATION_OUTDIR=/path
  RJ_DAG_DRYRUN=1
EOF
}

for tok in "$@"; do
  case "$tok" in
    SOURCE=*) SOURCE="${tok#SOURCE=}" ;;
    SWEEP_DIR=*) SWEEP_DIR="${tok#SWEEP_DIR=}" ;;
    SUB_ROOT=*) SUB_ROOT="${tok#SUB_ROOT=}" ;;
    REQUEST_MEMORY=*) REQUEST_MEMORY="${tok#REQUEST_MEMORY=}" ;;
    MAXJOBS=*) MAXJOBS="${tok#MAXJOBS=}" ;;
    EPOCHS=*) TRAIN_EPOCHS="${tok#EPOCHS=}" ;;
    PATIENCE=*) TRAIN_PATIENCE="${tok#PATIENCE=}" ;;
    RESTARTS=*) TRAIN_RESTARTS="${tok#RESTARTS=}" ;;
    HIDDEN_LAYERS=*) TRAIN_HIDDEN="${tok#HIDDEN_LAYERS=}" ;;
    HIDDEN_LAYER_GRID=*) TRAIN_GRID="${tok#HIDDEN_LAYER_GRID=}" ;;
    BATCH_SIZE=*) TRAIN_BATCH="${tok#BATCH_SIZE=}" ;;
    LEARNING_RATE=*) TRAIN_LR="${tok#LEARNING_RATE=}" ;;
    L2=*) TRAIN_L2="${tok#L2=}" ;;
    MIN_DELTA=*) TRAIN_MIN_DELTA="${tok#MIN_DELTA=}" ;;
    MAX_ROWS_PER_PT_BIN_CLASS=*) ROUTE_CAP="${tok#MAX_ROWS_PER_PT_BIN_CLASS=}" ;;
    DISTILLATION_STRENGTH=*) DISTILLATION_STRENGTH="${tok#DISTILLATION_STRENGTH=}" ;;
    CONDITIONAL_JITTER=*) CONDITIONAL_JITTER="${tok#CONDITIONAL_JITTER=}" ;;
    VALIDATION_CACHE=*|CACHE=*) VALIDATION_CACHE="${tok#*=}" ;;
    VALIDATION_OUTDIR=*) VALIDATION_OUTDIR="${tok#VALIDATION_OUTDIR=}" ;;
    VALIDATION_PT_BINS=*) VALIDATION_PT_BINS="${tok#VALIDATION_PT_BINS=}" ;;
    -h|--help|help) usage; exit 0 ;;
    *) die "Unknown argument: $tok" ;;
  esac
done

VALIDATION_OUTDIR="${VALIDATION_OUTDIR:-${SWEEP_DIR}/validation_rescore}"

[[ -n "$SOURCE" && -d "$SOURCE" ]] || die "SOURCE=/path/to/extraction is required"
[[ -n "$VALIDATION_CACHE" && -s "$VALIDATION_CACHE" ]] || die "VALIDATION_CACHE=/path/score_caches.list is required"
[[ "$MAXJOBS" =~ ^[0-9]+$ && "$MAXJOBS" -gt 0 ]] || die "MAXJOBS must be positive"
[[ "$ROUTE_CAP" =~ ^[0-9]+$ && "$ROUTE_CAP" -gt 0 ]] || die "MAX_ROWS_PER_PT_BIN_CLASS must be positive"
if [[ "${RJ_DAG_DRYRUN:-0}" != "1" ]]; then
  command -v condor_submit_dag >/dev/null 2>&1 || die "condor_submit_dag not found in PATH"
fi

mkdir -p "$SWEEP_DIR" "$SUB_ROOT"

say "RECOILJETS_AUAU_MLP_FINEPT_DISTILLED_SUBMIT_V1"
say "submit_host=$(hostname -f 2>/dev/null || hostname)"
say "timestamp=$STAMP"
say "source=$SOURCE"
say "sweep_dir=$SWEEP_DIR"
say "sub_root=$SUB_ROOT"
say "request_memory=$REQUEST_MEMORY maxjobs=$MAXJOBS"
say "epochs=$TRAIN_EPOCHS patience=$TRAIN_PATIENCE restarts=$TRAIN_RESTARTS hidden=$TRAIN_HIDDEN grid=$TRAIN_GRID"
say "route_cap=$ROUTE_CAP distillation_strength=$DISTILLATION_STRENGTH conditional_jitter=$CONDITIONAL_JITTER"
say "validation_cache=$VALIDATION_CACHE"
say "validation_outdir=$VALIDATION_OUTDIR validation_pt_bins=$VALIDATION_PT_BINS"

worker="${SUB_ROOT}/train_finept_worker.sh"
cat > "$worker" <<EOF
#!/usr/bin/env bash
set -euo pipefail
label="\${1:?label}"
model_dir="\${2:?model_dir}"
pt_range="\${3:?pt_range}"
train_pt_bins="\${4:?train_pt_bins}"
export RJ_ML_PYTHON="${ML_PYTHON}"
route_weight="\${pt_range}:1.0"
export RJ_AUAU_TIGHT_MLP_PT_RANGE="\$pt_range"
export RJ_AUAU_MLP_TRAIN_PT_BINS="\$train_pt_bins"
export RJ_AUAU_MLP_TRAIN_PT_BIN_WEIGHT_MODE=none
export RJ_AUAU_MLP_TRAIN_HIGHPT_SELECTION_WEIGHTS="\$route_weight"
export RJ_AUAU_MLP_TRAIN_SELECTION_METRIC=wp_fake_rate
cd "${RJ_REPO_BASE}"
echo "[auauMLPFinePtDistilled] start label=\$label model_dir=\$model_dir pt_range=\$pt_range train_pt_bins=\$train_pt_bins"
echo "[auauMLPFinePtDistilled] route_highpt_selection_weights=\$RJ_AUAU_MLP_TRAIN_HIGHPT_SELECTION_WEIGHTS"
./scripts/auau_tight_mlp_pipeline.sh trainHighPtDistilledKitchenV2FromExtraction \\
  SOURCE="${SOURCE}" \\
  MODEL_DIR="\$model_dir" \\
  PT_RANGE="\$pt_range" \\
  TRAIN_PT_BINS="\$train_pt_bins" \\
  PT_BIN_WEIGHT_MODE=none \\
  HIGHPT_SELECTION_WEIGHTS="\$route_weight" \\
  SELECTION_METRIC=wp_fake_rate \\
  MAX_ROWS_PER_PT_BIN_CLASS="${ROUTE_CAP}" \\
  EPOCHS="${TRAIN_EPOCHS}" \\
  PATIENCE="${TRAIN_PATIENCE}" \\
  RESTARTS="${TRAIN_RESTARTS}" \\
  HIDDEN_LAYERS="${TRAIN_HIDDEN}" \\
  HIDDEN_LAYER_GRID="${TRAIN_GRID}" \\
  BATCH_SIZE="${TRAIN_BATCH}" \\
  LEARNING_RATE="${TRAIN_LR}" \\
  L2="${TRAIN_L2}" \\
  MIN_DELTA="${TRAIN_MIN_DELTA}" \\
  DISTILLATION_STRENGTH="${DISTILLATION_STRENGTH}" \\
  CONDITIONAL_JITTER="${CONDITIONAL_JITTER}"
./scripts/auau_tight_mlp_pipeline.sh applyCheck MODEL_DIR="\$model_dir"
echo "[auauMLPFinePtDistilled] done label=\$label model_dir=\$model_dir"
EOF
chmod +x "$worker"

dag="${SUB_ROOT}/auau_mlp_finept_distilled_sweep.dag"
: > "$dag"
echo "MAXJOBS train ${MAXJOBS}" >> "$dag"
: > "${SUB_ROOT}/_parents.tmp"

add_job() {
  local label="$1" model_dir="$2" pt_range="$3" train_pt_bins="$4"
  mkdir -p "$model_dir"
  local sub="${SUB_ROOT}/${label}.sub"
  cat > "$sub" <<EOF
universe = vanilla
executable = ${worker}
arguments = ${label} ${model_dir} ${pt_range} ${train_pt_bins}
output = ${SUB_ROOT}/${label}.out
error = ${SUB_ROOT}/${label}.err
log = ${SUB_ROOT}/${label}.log
request_memory = ${REQUEST_MEMORY}
request_cpus = 1
+JobBatchName = "auau_mlp_finept_distilled_${STAMP}"
notification = Never
queue
EOF
  echo "JOB ${label} ${sub}" >> "$dag"
  echo "CATEGORY ${label} train" >> "$dag"
  echo "RETRY ${label} 1" >> "$dag"
  echo "PARENT ${label} CHILD FINALIZE" >> "${SUB_ROOT}/_parents.tmp"
}

add_job PT_015_018 "${SWEEP_DIR}/pt_015_018" "15:18" "15,18"
add_job PT_018_020 "${SWEEP_DIR}/pt_018_020" "18:20" "18,20"
add_job PT_020_022p5 "${SWEEP_DIR}/pt_020_022p5" "20:22.5" "20,22.5"
add_job PT_022p5_025 "${SWEEP_DIR}/pt_022p5_025" "22.5:25" "22.5,25"
add_job PT_025_030 "${SWEEP_DIR}/pt_025_030" "25:30" "25,30"
add_job PT_030_035 "${SWEEP_DIR}/pt_030_035" "30:35" "30,35"

finalize="${SUB_ROOT}/finalize_finept_distilled.sh"
cat > "$finalize" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export RJ_ML_PYTHON="${ML_PYTHON}"
cd "${RJ_REPO_BASE}"
manifest="${SWEEP_DIR}/finept_distilled_kitchen_v2_sweep_manifest.json"
python3 - <<'PY'
import json
from pathlib import Path

sweep_dir = Path("${SWEEP_DIR}")
artifact_name = "auau_tight_mlp_highPtDistilledKitchen_v2.json"
routes = [
    (15.0, 18.0, "pt_015_018", sweep_dir / "pt_015_018" / artifact_name),
    (18.0, 20.0, "pt_018_020", sweep_dir / "pt_018_020" / artifact_name),
    (20.0, 22.5, "pt_020_022p5", sweep_dir / "pt_020_022p5" / artifact_name),
    (22.5, 25.0, "pt_022p5_025", sweep_dir / "pt_022p5_025" / artifact_name),
    (25.0, 30.0, "pt_025_030", sweep_dir / "pt_025_030" / artifact_name),
    (30.0, 35.0, "pt_030_035", sweep_dir / "pt_030_035" / artifact_name),
]
missing = [str(path) for _, _, _, path in routes if not path.is_file()]
if missing:
    raise SystemExit("Missing trained fine-pT MLP artifacts:\\n  " + "\\n  ".join(missing))
payload = {
    "schema": "RJ_AUAU_TIGHT_MLP_SWEEP_V1",
    "source": "${SOURCE}",
    "sweep_dir": str(sweep_dir),
    "variants": [
        {
            "name": "finePtDistilledKitchenMLP_v2",
            "kind": "routed",
            "variant": "auauHighPtDistilledKitchenMLP_finePtValidationOnly",
            "axis": "cluster_Et",
            "routes": [
                {"lo": lo, "hi": hi, "label": label, "artifact": str(path)}
                for lo, hi, label, path in routes
            ],
        }
    ],
}
sweep_dir.mkdir(parents=True, exist_ok=True)
Path("${SWEEP_DIR}/finept_distilled_kitchen_v2_sweep_manifest.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\\n")
print("${SWEEP_DIR}/finept_distilled_kitchen_v2_sweep_manifest.json")
PY
./scripts/auau_tight_mlp_pipeline.sh rescoreValidationCache \\
  CACHE="${VALIDATION_CACHE}" \\
  SWEEP_MANIFEST="\$manifest" \\
  OUTDIR="${VALIDATION_OUTDIR}" \\
  PT_BINS="${VALIDATION_PT_BINS}"
echo "DONE_FINEPT_DISTILLED_SWEEP_DIR=${SWEEP_DIR}"
echo "DONE_FINEPT_DISTILLED_MANIFEST=\$manifest"
echo "DONE_FINEPT_DISTILLED_VALIDATION=${VALIDATION_OUTDIR}"
EOF
chmod +x "$finalize"

final_sub="${SUB_ROOT}/FINALIZE.sub"
cat > "$final_sub" <<EOF
universe = scheduler
executable = ${finalize}
output = ${SUB_ROOT}/FINALIZE.out
error = ${SUB_ROOT}/FINALIZE.err
log = ${SUB_ROOT}/FINALIZE.log
notification = Never
queue
EOF
echo "JOB FINALIZE ${final_sub}" >> "$dag"
cat "${SUB_ROOT}/_parents.tmp" >> "$dag"
rm -f "${SUB_ROOT}/_parents.tmp"

say "dag=$dag"
say "sweep_manifest=${SWEEP_DIR}/finept_distilled_kitchen_v2_sweep_manifest.json"
if [[ "${RJ_DAG_DRYRUN:-0}" == "1" ]]; then
  say "RJ_DAG_DRYRUN=1; not submitting"
  exit 0
fi

condor_submit_dag "$dag"
