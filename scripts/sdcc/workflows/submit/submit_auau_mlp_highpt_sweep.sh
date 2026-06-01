#!/usr/bin/env bash
set -euo pipefail

RJ_REPO_BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
readonly RJ_REPO_BASE

ML_PYTHON="${RJ_ML_PYTHON:-/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python}"
MODEL_BASE="${RJ_AUAU_MLP_MODEL_BASE:-${RJ_REPO_BASE}/mlp_models}"
SOURCE="${RJ_AUAU_MLP_SWEEP_SOURCE:-}"
STAMP="${RJ_AUAU_MLP_SWEEP_STAMP:-$(date +%Y%m%d_%H%M%S)}"
SWEEP_DIR="${RJ_AUAU_MLP_SWEEP_DIR:-${MODEL_BASE}/tight_mlp_highpt_sweep_${STAMP}}"
SUB_ROOT="${RJ_AUAU_MLP_SWEEP_SUB_ROOT:-${RJ_REPO_BASE}/condor_sub/auauTightMLPHighPtSweep_${STAMP}}"
REQUEST_MEMORY="${RJ_AUAU_MLP_SWEEP_REQUEST_MEMORY:-24000MB}"
MAXJOBS="${RJ_AUAU_MLP_SWEEP_MAXJOBS:-4}"
TRAIN_EPOCHS="${RJ_AUAU_MLP_SWEEP_EPOCHS:-220}"
TRAIN_PATIENCE="${RJ_AUAU_MLP_SWEEP_PATIENCE:-45}"
TRAIN_RESTARTS="${RJ_AUAU_MLP_SWEEP_RESTARTS:-1}"
TRAIN_HIDDEN="${RJ_AUAU_MLP_SWEEP_HIDDEN_LAYERS:-128,64,32}"
TRAIN_GRID="${RJ_AUAU_MLP_SWEEP_HIDDEN_LAYER_GRID:-160,80,40}"
TRAIN_BATCH="${RJ_AUAU_MLP_SWEEP_BATCH_SIZE:-4096}"
VALIDATION_CACHE="${RJ_AUAU_MLP_SWEEP_VALIDATION_CACHE:-}"
VALIDATION_OUTDIR="${RJ_AUAU_MLP_SWEEP_VALIDATION_OUTDIR:-${SWEEP_DIR}/validation_rescore}"
VALIDATION_OUTDIR_SET=0
VALIDATION_PT_BINS="${RJ_AUAU_MLP_SWEEP_VALIDATION_PT_BINS:-15,20,25,35}"

say() { printf '\033[1;36m[auauMLPHighPtSweep]\033[0m %s\n' "$*"; }
err() { printf '\033[1;31m[auauMLPHighPtSweep][ERR]\033[0m %s\n' "$*" >&2; }
die() { err "$*"; exit 2; }

usage() {
  cat <<'EOF'
Usage:
  ./scripts/submit_auau_mlp_highpt_sweep.sh SOURCE=/path/to/extraction [VALIDATION_CACHE=/path/score_caches.list]

Submits an isolated Condor DAG that trains the high-pT AuAu tight-MLP sweep:
  global15_highPtBalanced, global5_broadReach, pt3_15to35, cent3_15to35.

Useful knobs:
  SWEEP_DIR=/path
  REQUEST_MEMORY=24000MB
  MAXJOBS=4
  EPOCHS=220
  PATIENCE=45
  VALIDATION_OUTDIR=/path
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
    VALIDATION_CACHE=*|CACHE=*) VALIDATION_CACHE="${tok#*=}" ;;
    VALIDATION_OUTDIR=*) VALIDATION_OUTDIR="${tok#VALIDATION_OUTDIR=}"; VALIDATION_OUTDIR_SET=1 ;;
    VALIDATION_PT_BINS=*) VALIDATION_PT_BINS="${tok#VALIDATION_PT_BINS=}" ;;
    -h|--help|help) usage; exit 0 ;;
  esac
done

if [[ "$VALIDATION_OUTDIR_SET" == "0" && -z "${RJ_AUAU_MLP_SWEEP_VALIDATION_OUTDIR:-}" ]]; then
  VALIDATION_OUTDIR="${SWEEP_DIR}/validation_rescore"
fi

[[ -n "$SOURCE" && -d "$SOURCE" ]] || die "SOURCE=/path/to/extraction is required"
[[ "$MAXJOBS" =~ ^[0-9]+$ && "$MAXJOBS" -gt 0 ]] || die "MAXJOBS must be positive"
mkdir -p "$SWEEP_DIR" "$SUB_ROOT"

say "RECOILJETS_AUAU_MLP_HIGHPT_SWEEP_SUBMIT_V1"
say "submit_host=$(hostname -f 2>/dev/null || hostname)"
say "timestamp=$STAMP"
say "source=$SOURCE"
say "sweep_dir=$SWEEP_DIR"
say "sub_root=$SUB_ROOT"
say "request_memory=$REQUEST_MEMORY maxjobs=$MAXJOBS"
say "epochs=$TRAIN_EPOCHS patience=$TRAIN_PATIENCE restarts=$TRAIN_RESTARTS hidden=$TRAIN_HIDDEN grid=$TRAIN_GRID"
[[ -n "$VALIDATION_CACHE" ]] && say "validation_cache=$VALIDATION_CACHE"
say "validation_outdir=$VALIDATION_OUTDIR validation_pt_bins=$VALIDATION_PT_BINS"

worker="${SUB_ROOT}/train_worker.sh"
cat > "$worker" <<EOF
#!/usr/bin/env bash
set -euo pipefail
variant="\${1:?variant}"
model_dir="\${2:?model_dir}"
pt_range="\${3:?pt_range}"
centrality_range="\${4:?centrality_range}"
train_pt_bins="\${5:?train_pt_bins}"
pt_weight_spec="\${6:?pt_weight_spec}"
highpt_weights="\${7:?highpt_weights}"
pt_bin_cap="\${8:?pt_bin_cap}"
export RJ_ML_PYTHON="${ML_PYTHON}"
export RJ_AUAU_TIGHT_MLP_PRODUCTS="primary-ratios"
export RJ_AUAU_MLP_TRAIN_EPOCHS="${TRAIN_EPOCHS}"
export RJ_AUAU_MLP_TRAIN_PATIENCE="${TRAIN_PATIENCE}"
export RJ_AUAU_MLP_TRAIN_RESTARTS="${TRAIN_RESTARTS}"
export RJ_AUAU_MLP_TRAIN_HIDDEN_LAYERS="${TRAIN_HIDDEN}"
export RJ_AUAU_MLP_TRAIN_HIDDEN_LAYER_GRID="${TRAIN_GRID}"
export RJ_AUAU_MLP_TRAIN_BATCH_SIZE="${TRAIN_BATCH}"
export RJ_AUAU_MLP_TRAIN_PROGRESS_EVERY="\${RJ_AUAU_MLP_TRAIN_PROGRESS_EVERY:-2}"
export RJ_AUAU_MLP_TRAIN_SELECTION_METRIC="highpt_wp80"
export RJ_AUAU_MLP_TRAIN_SELECTION_TARGET_SIGNAL_EFF="0.80"
export RJ_AUAU_MLP_TRAIN_PT_BIN_WEIGHT_MODE="highpt"
export RJ_AUAU_MLP_TRAIN_MAX_ROWS_PER_CLASS="0"
export RJ_AUAU_MLP_TRAIN_MAX_ROWS_PER_PT_BIN_CLASS="\$pt_bin_cap"
cd "${RJ_REPO_BASE}"
args=(trainHighPtBalancedFromExtraction SOURCE="${SOURCE}" MODEL_DIR="\$model_dir" PT_RANGE="\$pt_range" TRAIN_PT_BINS="\$train_pt_bins")
[[ "\$centrality_range" != "-" ]] && args+=(CENTRALITY_RANGE="\$centrality_range")
[[ "\$pt_weight_spec" != "-" ]] && args+=(PT_BIN_WEIGHT_SPEC="\$pt_weight_spec")
[[ "\$highpt_weights" != "-" ]] && args+=(HIGHPT_SELECTION_WEIGHTS="\$highpt_weights")
echo "[auauMLPHighPtSweep] start variant=\$variant model_dir=\$model_dir pt_range=\$pt_range centrality_range=\$centrality_range train_pt_bins=\$train_pt_bins"
./scripts/auau_tight_mlp_pipeline.sh "\${args[@]}"
./scripts/auau_tight_mlp_pipeline.sh applyCheck MODEL_DIR="\$model_dir"
echo "[auauMLPHighPtSweep] done variant=\$variant model_dir=\$model_dir"
EOF
chmod +x "$worker"

dag="${SUB_ROOT}/auau_mlp_highpt_sweep.dag"
: > "$dag"
echo "MAXJOBS train ${MAXJOBS}" >> "$dag"

add_job() {
  local label="$1" model_dir="$2" pt_range="$3" cent_range="$4" train_pt_bins="$5" pt_weight_spec="$6" highpt_weights="$7" cap="$8"
  mkdir -p "$model_dir"
  local sub="${SUB_ROOT}/${label}.sub"
  cat > "$sub" <<EOF
universe = vanilla
executable = ${worker}
arguments = ${label} ${model_dir} ${pt_range} ${cent_range} ${train_pt_bins} ${pt_weight_spec} ${highpt_weights} ${cap}
output = ${SUB_ROOT}/${label}.out
error = ${SUB_ROOT}/${label}.err
log = ${SUB_ROOT}/${label}.log
request_memory = ${REQUEST_MEMORY}
request_cpus = 1
+JobBatchName = "auau_mlp_highpt_sweep_${STAMP}"
notification = Never
queue
EOF
  echo "JOB ${label} ${sub}" >> "$dag"
  echo "CATEGORY ${label} train" >> "$dag"
  echo "RETRY ${label} 1" >> "$dag"
  echo "PARENT ${label} CHILD FINALIZE" >> "${SUB_ROOT}/_parents.tmp"
}

: > "${SUB_ROOT}/_parents.tmp"
add_job GLOBAL15 "${SWEEP_DIR}/global15_highPtBalanced" "15:35" "-" "15,20,25,35" "15:20:0.20,20:25:0.35,25:35:0.45" "15:20:0.20,20:25:0.35,25:35:0.45" "120000"
add_job GLOBAL5 "${SWEEP_DIR}/global5_broadReach" "5:35" "-" "5,10,15,20,25,35" "5:10:0.10,10:15:0.15,15:20:0.20,20:25:0.25,25:35:0.30" "5:10:0.03,10:15:0.07,15:20:0.20,20:25:0.30,25:35:0.40" "90000"
add_job PT3_015_020 "${SWEEP_DIR}/pt3_15to35/pt_015_020" "15:20" "-" "15,20" "-" "-" "180000"
add_job PT3_020_025 "${SWEEP_DIR}/pt3_15to35/pt_020_025" "20:25" "-" "20,25" "-" "-" "180000"
add_job PT3_025_035 "${SWEEP_DIR}/pt3_15to35/pt_025_035" "25:35" "-" "25,35" "-" "-" "180000"
add_job CENT3_000_020 "${SWEEP_DIR}/cent3_15to35/cent_000_020" "15:35" "0:20" "15,20,25,35" "15:20:0.20,20:25:0.35,25:35:0.45" "15:20:0.20,20:25:0.35,25:35:0.45" "90000"
add_job CENT3_020_050 "${SWEEP_DIR}/cent3_15to35/cent_020_050" "15:35" "20:50" "15,20,25,35" "15:20:0.20,20:25:0.35,25:35:0.45" "15:20:0.20,20:25:0.35,25:35:0.45" "90000"
add_job CENT3_050_080 "${SWEEP_DIR}/cent3_15to35/cent_050_080" "15:35" "50:80" "15,20,25,35" "15:20:0.20,20:25:0.35,25:35:0.45" "15:20:0.20,20:25:0.35,25:35:0.45" "90000"

finalize="${SUB_ROOT}/finalize_sweep.sh"
cat > "$finalize" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export RJ_ML_PYTHON="${ML_PYTHON}"
cd "${RJ_REPO_BASE}"
manifest="${SWEEP_DIR}/sweep_manifest.json"
python3 - <<'PY'
import json
from pathlib import Path
sweep_dir = Path("${SWEEP_DIR}")
artifact_name = "auau_tight_mlp_centInputBase3x3_pt1535.json"
variants = [
    {
        "name": "global15_highPtBalanced",
        "kind": "single",
        "variant": "auauCentInputBase3x3MLP",
        "artifact": str(sweep_dir / "global15_highPtBalanced" / artifact_name),
        "training": {"pt_range": "15:35", "train_pt_bins": "15,20,25,35"},
    },
    {
        "name": "global5_broadReach",
        "kind": "single",
        "variant": "auauCentInputBase3x3MLP",
        "artifact": str(sweep_dir / "global5_broadReach" / artifact_name),
        "training": {"pt_range": "5:35", "train_pt_bins": "5,10,15,20,25,35"},
    },
    {
        "name": "pt3_15to35",
        "kind": "routed",
        "variant": "auauCentInputBase3x3MLP_pt3ValidationOnly",
        "axis": "cluster_Et",
        "routes": [
            {"lo": 15.0, "hi": 20.0, "label": "pt_015_020", "artifact": str(sweep_dir / "pt3_15to35" / "pt_015_020" / artifact_name)},
            {"lo": 20.0, "hi": 25.0, "label": "pt_020_025", "artifact": str(sweep_dir / "pt3_15to35" / "pt_020_025" / artifact_name)},
            {"lo": 25.0, "hi": 35.0, "label": "pt_025_035", "artifact": str(sweep_dir / "pt3_15to35" / "pt_025_035" / artifact_name)},
        ],
    },
    {
        "name": "cent3_15to35",
        "kind": "routed",
        "variant": "auauCentInputBase3x3MLP_cent3ValidationOnly",
        "axis": "centrality",
        "routes": [
            {"lo": 0.0, "hi": 20.0, "label": "cent_000_020", "artifact": str(sweep_dir / "cent3_15to35" / "cent_000_020" / artifact_name)},
            {"lo": 20.0, "hi": 50.0, "label": "cent_020_050", "artifact": str(sweep_dir / "cent3_15to35" / "cent_020_050" / artifact_name)},
            {"lo": 50.0, "hi": 80.0, "label": "cent_050_080", "artifact": str(sweep_dir / "cent3_15to35" / "cent_050_080" / artifact_name)},
        ],
    },
]
missing = []
for variant in variants:
    if variant["kind"] == "single":
        if not Path(variant["artifact"]).is_file():
            missing.append(variant["artifact"])
    else:
        for route in variant["routes"]:
            if not Path(route["artifact"]).is_file():
                missing.append(route["artifact"])
if missing:
    raise SystemExit("Missing trained MLP artifacts:\\n  " + "\\n  ".join(missing))
payload = {
    "schema": "RJ_AUAU_TIGHT_MLP_SWEEP_V1",
    "source": "${SOURCE}",
    "sweep_dir": str(sweep_dir),
    "variants": variants,
}
sweep_dir.mkdir(parents=True, exist_ok=True)
(sweep_dir / "sweep_manifest.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\\n")
print(sweep_dir / "sweep_manifest.json")
PY
if [[ -n "${VALIDATION_CACHE}" && -s "${VALIDATION_CACHE}" ]]; then
  ./scripts/auau_tight_mlp_pipeline.sh rescoreValidationCache CACHE="${VALIDATION_CACHE}" SWEEP_MANIFEST="\$manifest" OUTDIR="${VALIDATION_OUTDIR}" PT_BINS="${VALIDATION_PT_BINS}"
fi
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
say "sweep_manifest=${SWEEP_DIR}/sweep_manifest.json"
if [[ "${RJ_DAG_DRYRUN:-0}" == "1" ]]; then
  say "RJ_DAG_DRYRUN=1; not submitting"
  exit 0
fi

condor_submit_dag "$dag"
