#!/usr/bin/env bash
set -euo pipefail

MODE="${1:-plan}"
if [[ $# -gt 0 ]]; then
  shift
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
RJ_REPO_BASE="${RJ_REPO_BASE:-$(cd "${SCRIPT_DIR}/../../../.." && pwd -P)}"
PIPELINE="${RJ_THE38_BDT_PIPELINE:-${RJ_REPO_BASE}/scripts/auau_tight_bdt_pipeline.sh}"
ML_PYTHON="${RJ_ML_PYTHON:-/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python}"
TRAIN_SCRIPT="${RJ_THE38_TRAIN_SCRIPT:-${RJ_REPO_BASE}/scripts/ml/training/train_auau_photon_bdt.py}"
VALIDATE_SCRIPT="${RJ_THE38_VALIDATE_SCRIPT:-${RJ_REPO_BASE}/scripts/ml/validation/validate_auau_tight_bdt_on_sim.py}"

STAMP="${RJ_THE38_STAMP:-$(date +%Y%m%d_%H%M%S)}"
TAG="${RJ_THE38_TAG:-THE38_tree_depth_capacity_${STAMP}}"
DEPTHS_CSV="${RJ_THE38_DEPTHS:-2,3,4,5,6}"
SOURCE="${RJ_THE38_SOURCE:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/THE8_branchA_ladder_jet12_20_30_40_20260527}"
MODEL_BASE="${RJ_THE38_MODEL_BASE:-/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models}"
LOCAL_OUTPUT_ROOT="${RJ_THE38_LOCAL_OUTPUT_ROOT:-dataOutput/auauMLDiagnosticRuns/${TAG}}"
LOCAL_REGISTRY_ROOT="${RJ_THE38_LOCAL_REGISTRY_ROOT:-${LOCAL_OUTPUT_ROOT}/registries}"
SLIDE_OUTDIR="${RJ_THE38_SLIDE_OUTDIR:-${LOCAL_OUTPUT_ROOT}/slideReady}"
SLIDE_GENERATOR="${RJ_THE38_SLIDE_GENERATOR:-${RJ_REPO_BASE}/scripts/slides/split_studies/tree_depth/make_the38_tree_depth_capacity_slide.py}"

SPEC_ID="${RJ_THE38_SPEC_ID:-centInput_pt1535}"
CAMPAIGN="${RJ_THE38_CAMPAIGN:-etfine-centstudy}"
WEIGHT_MODE="${RJ_THE38_WEIGHT_MODE:-ppg12-exact}"
TEST_SIZE="${RJ_THE38_TEST_SIZE:-0.10}"
SPLIT_MODE="${RJ_THE38_SPLIT_MODE:-row}"
EXPECTED_SAMPLES="${RJ_THE38_EXPECTED_SAMPLES:-run28_embeddedPhoton12,run28_embeddedPhoton20,run28_embeddedJet12,run28_embeddedJet20,run28_embeddedJet30,run28_embeddedJet40}"
TRAIN_GROUP_SIZE="${RJ_THE38_TRAIN_GROUP_SIZE:-1}"
TRAIN_REQUEST_MEMORY="${RJ_THE38_TRAIN_REQUEST_MEMORY:-16000MB}"
CACHE_REQUEST_MEMORY="${RJ_THE38_CACHE_REQUEST_MEMORY:-24000MB}"
VALIDATE_GROUP_SIZE="${RJ_THE38_VALIDATE_GROUP_SIZE:-50}"
VALIDATE_REQUEST_MEMORY="${RJ_THE38_VALIDATE_REQUEST_MEMORY:-2500MB}"
VALIDATE_TOTAL_SCORE_MAX="${RJ_THE38_VALIDATE_TOTAL_SCORE_MAX:-0}"
VALIDATE_STAMP="${RJ_THE38_VALIDATE_STAMP:-${STAMP}}"
ALLOW_EXISTING="${RJ_THE38_ALLOW_EXISTING:-0}"
REUSE_EXISTING_MANIFEST="${RJ_THE38_REUSE_EXISTING_MANIFEST:-1}"
SKIP_TRAINING_TREE_VALIDATE="${RJ_THE38_SKIP_TRAINING_TREE_VALIDATE:-1}"

say() { printf '\033[1;36m[the38Depth]\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[the38Depth][WARN]\033[0m %s\n' "$*" >&2; }
die() { printf '\033[1;31m[the38Depth][ERR]\033[0m %s\n' "$*" >&2; exit 2; }

usage() {
  cat <<'EOF'
Usage:
  submit_the38_tree_depth_capacity_campaign.sh plan [KEY=VALUE ...]
  submit_the38_tree_depth_capacity_campaign.sh submit-training [KEY=VALUE ...]
  submit_the38_tree_depth_capacity_campaign.sh status [KEY=VALUE ...]
  submit_the38_tree_depth_capacity_campaign.sh render-primary-slide [KEY=VALUE ...]
  submit_the38_tree_depth_capacity_campaign.sh submit-fullscore-check [KEY=VALUE ...]

Modes:
  plan
      Path guard plus RJ_DAG_DRYRUN=1 training dry-runs for each depth.
  submit-training
      Submit the five depth-specific training DAGs. Requires
      RJ_THE38_APPROVED=1, RJ_THE38_DUPLICATE_GUARD_PASSED=1,
      RJ_CODEX_CHAT_NAME, and RJ_CODEX_THREAD_ID.
  status
      Read-only local/SDCC status summary for depth registries and queue rows.
  render-primary-slide
      Render the local train/holdout explanatory PNG from pulled/read
      registry JSON files. This does not submit validation.
  submit-fullscore-check
      Optional later full score-cache validation sanity check. Requires the
      same approval/provenance gates as submit-training.

Common overrides:
  SOURCE=/path TAG=name STAMP=stamp DEPTHS=2,3,4,5,6
  LOCAL_REGISTRY_ROOT=/local/path SLIDE_OUTDIR=/local/path
EOF
}

tag_override=0
local_output_override=0
local_registry_override=0
slide_outdir_override=0

for tok in "$@"; do
  case "$tok" in
    SOURCE=*) SOURCE="${tok#SOURCE=}" ;;
    TAG=*) TAG="${tok#TAG=}"; tag_override=1 ;;
    STAMP=*) STAMP="${tok#STAMP=}" ;;
    DEPTHS=*|DEPTHS_CSV=*) DEPTHS_CSV="${tok#*=}" ;;
    MODEL_BASE=*) MODEL_BASE="${tok#MODEL_BASE=}" ;;
    LOCAL_OUTPUT_ROOT=*) LOCAL_OUTPUT_ROOT="${tok#LOCAL_OUTPUT_ROOT=}"; local_output_override=1 ;;
    LOCAL_REGISTRY_ROOT=*) LOCAL_REGISTRY_ROOT="${tok#LOCAL_REGISTRY_ROOT=}"; local_registry_override=1 ;;
    SLIDE_OUTDIR=*) SLIDE_OUTDIR="${tok#SLIDE_OUTDIR=}"; slide_outdir_override=1 ;;
    VALIDATE_TOTAL_SCORE_MAX=*|scoreMaxRows=*) VALIDATE_TOTAL_SCORE_MAX="${tok#*=}" ;;
    -h|--help|help) usage; exit 0 ;;
    *) die "Unknown argument: $tok" ;;
  esac
done

if [[ "$tag_override" == "0" && -z "${RJ_THE38_TAG:-}" ]]; then
  TAG="THE38_tree_depth_capacity_${STAMP}"
fi
if [[ "$local_output_override" == "0" && -z "${RJ_THE38_LOCAL_OUTPUT_ROOT:-}" ]]; then
  LOCAL_OUTPUT_ROOT="dataOutput/auauMLDiagnosticRuns/${TAG}"
fi
if [[ "$local_registry_override" == "0" && -z "${RJ_THE38_LOCAL_REGISTRY_ROOT:-}" ]]; then
  LOCAL_REGISTRY_ROOT="${LOCAL_OUTPUT_ROOT}/registries"
fi
if [[ "$slide_outdir_override" == "0" && -z "${RJ_THE38_SLIDE_OUTDIR:-}" ]]; then
  SLIDE_OUTDIR="${LOCAL_OUTPUT_ROOT}/slideReady"
fi

IFS=',' read -r -a DEPTHS <<< "$DEPTHS_CSV"

model_dir_for_depth() {
  local depth="$1"
  printf '%s/%s_d%s' "$MODEL_BASE" "$TAG" "$depth"
}

train_stamp_for_depth() {
  local depth="$1"
  printf '%s_d%s' "$TAG" "$depth"
}

registry_for_depth() {
  local depth="$1"
  printf '%s/%s_d%s/model_registry.json' "$MODEL_BASE" "$TAG" "$depth"
}

local_registry_for_depth() {
  local depth="$1"
  local override_var="RJ_THE38_REGISTRY_D${depth}"
  if [[ -n "${!override_var:-}" ]]; then
    printf '%s' "${!override_var}"
  else
    printf '%s/d%s/model_registry.json' "$LOCAL_REGISTRY_ROOT" "$depth"
  fi
}

report_dir_for_depth() {
  local depth="$1"
  printf '%s/reports/model_validation_condor_%s_d%s_fullscore_%s' "$SOURCE" "$TAG" "$depth" "$VALIDATE_STAMP"
}

print_plan() {
  say "schema=THE38_TREE_DEPTH_CAPACITY_CAMPAIGN_V1"
  say "mode=$MODE"
  say "stamp=$STAMP"
  say "tag=$TAG"
  say "repo=$RJ_REPO_BASE"
  say "pipeline=$PIPELINE"
  say "train_script=$TRAIN_SCRIPT"
  say "validate_script=$VALIDATE_SCRIPT"
  say "source=$SOURCE"
  say "depths=$DEPTHS_CSV"
  say "campaign=$CAMPAIGN spec_id=$SPEC_ID weight_mode=$WEIGHT_MODE"
  say "split=${SPLIT_MODE} test_size=${TEST_SIZE}"
  say "expected_samples=$EXPECTED_SAMPLES"
  say "model_base=$MODEL_BASE"
  say "reuse_existing_manifest=$REUSE_EXISTING_MANIFEST"
  say "skip_training_tree_validate=$SKIP_TRAINING_TREE_VALIDATE"
  say "local_registry_root=$LOCAL_REGISTRY_ROOT"
  say "slide_outdir=$SLIDE_OUTDIR"
  for depth in "${DEPTHS[@]}"; do
    say "depth_${depth}_model=$(model_dir_for_depth "$depth")"
  done
}

need_file() {
  [[ -s "$1" ]] || die "missing required file: $1"
}

need_dir() {
  [[ -d "$1" ]] || die "missing required directory: $1"
}

source_manifest_guard() {
  [[ "$REUSE_EXISTING_MANIFEST" == "1" ]] || return 0
  local manifest="${SOURCE}/manifests/training_roots.list"
  need_file "$manifest"
  local nroots
  nroots="$(wc -l < "$manifest" | tr -d ' ')"
  [[ "$nroots" =~ ^[0-9]+$ && "$nroots" -gt 100 ]] || die "source manifest has too few roots: $manifest lines=$nroots"
  if grep -F "/reports/" "$manifest" >/dev/null 2>&1; then
    die "source manifest includes report products; refusing to train from $manifest"
  fi
  local sample
  IFS=',' read -r -a expected_arr <<< "$EXPECTED_SAMPLES"
  for sample in "${expected_arr[@]}"; do
    grep -F "$sample" "$manifest" >/dev/null 2>&1 || die "source manifest missing expected sample ${sample}: $manifest"
  done
}

require_submit_approval() {
  [[ "${RJ_THE38_APPROVED:-0}" == "1" ]] || die "$MODE requires RJ_THE38_APPROVED=1"
  [[ "${RJ_THE38_DUPLICATE_GUARD_PASSED:-0}" == "1" ]] || die "$MODE requires RJ_THE38_DUPLICATE_GUARD_PASSED=1"
  [[ -n "${RJ_CODEX_CHAT_NAME:-}" ]] || die "$MODE requires RJ_CODEX_CHAT_NAME"
  [[ -n "${RJ_CODEX_THREAD_ID:-}" ]] || die "$MODE requires RJ_CODEX_THREAD_ID"
}

fingerprint() {
  cat <<EOF
purpose=THE-38 tree-depth capacity control
dataset=Jet40-inclusive Branch A AuAu photon-ID training source
source=${SOURCE}
campaign=${CAMPAIGN}
spec=${SPEC_ID}
model=BDT centInput_pt1535
depths=${DEPTHS_CSV}
split=${SPLIT_MODE} test_size=${TEST_SIZE}
weight_mode=${WEIGHT_MODE}
expected_samples=${EXPECTED_SAMPLES}
tag=${TAG}
EOF
}

path_guard() {
  need_file "$PIPELINE"
  need_file "$TRAIN_SCRIPT"
  need_file "$VALIDATE_SCRIPT"
  need_dir "$SOURCE"
  source_manifest_guard
  local depth model_dir registry tmva
  for depth in "${DEPTHS[@]}"; do
    model_dir="$(model_dir_for_depth "$depth")"
    registry="${model_dir}/model_registry.json"
    if [[ "$ALLOW_EXISTING" == "1" ]]; then
      continue
    fi
    if [[ -s "$registry" ]]; then
      die "target already has final registry: $registry"
    fi
    shopt -s nullglob
    for tmva in "${model_dir}"/*_tmva.root; do
      shopt -u nullglob
      die "target already has trained TMVA file: $tmva"
    done
    shopt -u nullglob
  done
}

run_training_depth() {
  local depth="$1"
  local dryrun="${2:-0}"
  local model_dir train_stamp
  model_dir="$(model_dir_for_depth "$depth")"
  train_stamp="$(train_stamp_for_depth "$depth")"
  say "training_depth=${depth} dryrun=${dryrun} model_dir=${model_dir}"
  RJ_ML_PYTHON="$ML_PYTHON" \
  RJ_AUAU_TIGHT_BDT_TRAIN_SCRIPT="$TRAIN_SCRIPT" \
  RJ_AUAU_TIGHT_BDT_VALIDATE_SCRIPT="$VALIDATE_SCRIPT" \
  RJ_AUAU_TIGHT_BDT_TRAIN_STAMP="$train_stamp" \
  RJ_AUAU_TIGHT_BDT_MODEL_DIR="$model_dir" \
  RJ_AUAU_BDT_CACHE_FILE="${model_dir}/training_matrix.npz" \
  RJ_AUAU_BDT_CAMPAIGN="$CAMPAIGN" \
  RJ_AUAU_BDT_CAMPAIGN_SPEC_IDS="$SPEC_ID" \
  RJ_AUAU_BDT_WEIGHT_MODE="$WEIGHT_MODE" \
  RJ_AUAU_BDT_PPG12_EXACT_EXPECTED_SAMPLES="$EXPECTED_SAMPLES" \
  RJ_AUAU_BDT_TEST_SIZE="$TEST_SIZE" \
  RJ_AUAU_BDT_SPLIT_MODE="$SPLIT_MODE" \
  RJ_AUAU_BDT_MAX_DEPTH="$depth" \
  RJ_AUAU_BDT_REUSE_EXISTING_MANIFEST="$REUSE_EXISTING_MANIFEST" \
  RJ_AUAU_BDT_SKIP_TRAINING_TREE_VALIDATE="$SKIP_TRAINING_TREE_VALIDATE" \
  RJ_AUAU_BDT_EXPANDED_GROUP_SIZE="$TRAIN_GROUP_SIZE" \
  RJ_AUAU_BDT_EXPANDED_REQUEST_MEMORY="$TRAIN_REQUEST_MEMORY" \
  RJ_AUAU_BDT_CACHE_REQUEST_MEMORY="$CACHE_REQUEST_MEMORY" \
  RJ_AUAU_BDT_TRAIN_PARALLEL_PER_JOB=1 \
  RJ_AUAU_BDT_XGB_N_JOBS=1 \
  RJ_DAG_DRYRUN="$dryrun" \
  "$PIPELINE" trainExpandedFromExtractionCondor SOURCE="$SOURCE" groupSize "$TRAIN_GROUP_SIZE"
}

plan_mode() {
  print_plan
  say "duplicate_fingerprint_begin"
  fingerprint
  say "duplicate_fingerprint_end"
  path_guard
  for depth in "${DEPTHS[@]}"; do
    run_training_depth "$depth" 1
  done
  say "plan_complete=no Condor submission performed"
}

submit_training() {
  print_plan
  require_submit_approval
  path_guard
  local depth
  for depth in "${DEPTHS[@]}"; do
    run_training_depth "$depth" 0
  done
  status_mode
}

status_mode() {
  print_plan
  if command -v condor_q >/dev/null 2>&1; then
    say "condor_q matching tag rows:"
    condor_q "${USER:-patsfan753}" -nobatch 2>/dev/null | grep -F "$TAG" || true
  else
    warn "condor_q not available"
  fi
  local registry_dirs=()
  local depth
  for depth in "${DEPTHS[@]}"; do
    registry_dirs+=("$(model_dir_for_depth "$depth")")
  done
  "$ML_PYTHON" - "$TAG" "$EXPECTED_SAMPLES" "${registry_dirs[@]}" <<'PY'
import json
import sys
from pathlib import Path

tag = sys.argv[1]
expected_samples = set(filter(None, sys.argv[2].split(",")))
registry_paths = [Path(item) / "model_registry.json" for item in sys.argv[3:]]
print(f"THE38_TREE_DEPTH_STATUS tag={tag}")
for registry_path in registry_paths:
    depth = registry_path.parent.name.rsplit("_d", 1)[-1]
    if not registry_path.exists():
        print(f"depth={depth} status=MISSING registry={registry_path}")
        continue
    payload = json.loads(registry_path.read_text())
    models = payload.get("models", [])
    model = next((m for m in models if m.get("model_id") == "centInput_pt1535"), models[0] if models else None)
    report = (model or {}).get("report") or {}
    overfit = report.get("overfit_diagnostics") or {}
    xgb = report.get("xgboost") or {}
    sample_validation_locations = [
        (report.get("weighting") or {}).get("sample_validation") or {},
        (report.get("ppg12_exact_closure") or {}).get("sample_validation") or {},
        report.get("sample_validation") or {},
    ]
    observed = set()
    for sample_validation in sample_validation_locations:
        observed = set(sample_validation.get("observed_samples") or [])
        if observed:
            break
    missing_expected = sorted(expected_samples - observed) if observed else []
    print(
        "depth={depth} status={status} max_depth={max_depth} "
        "train_auc={train_auc} holdout_auc={holdout_auc} "
        "auc_gap={auc_gap} logloss_gap={logloss_gap} "
        "missing_expected_samples={missing}".format(
            depth=depth,
            status=payload.get("status"),
            max_depth=xgb.get("max_depth"),
            train_auc=overfit.get("train_auc"),
            holdout_auc=overfit.get("holdout_auc"),
            auc_gap=overfit.get("auc_gap_train_minus_holdout"),
            logloss_gap=overfit.get("logloss_gap_holdout_minus_train"),
            missing=missing_expected,
        )
    )
PY
}

render_primary_slide() {
  print_plan
  need_file "$SLIDE_GENERATOR"
  local args=()
  local depth local_registry
  for depth in "${DEPTHS[@]}"; do
    local_registry="$(local_registry_for_depth "$depth")"
    args+=(--registry "${depth}=${local_registry}")
  done
  "$ML_PYTHON" "$SLIDE_GENERATOR" \
    "${args[@]}" \
    --outdir "$SLIDE_OUTDIR" \
    --tag "$TAG" \
    --baseline-depth 4 \
    --require-sample run28_embeddedJet40
}

submit_fullscore_check() {
  print_plan
  require_submit_approval
  local depth model_dir registry report_dir
  for depth in "${DEPTHS[@]}"; do
    model_dir="$(model_dir_for_depth "$depth")"
    registry="$(registry_for_depth "$depth")"
    report_dir="$(report_dir_for_depth "$depth")"
    need_file "$registry"
    say "fullscore_depth=${depth} model_dir=${model_dir} report_dir=${report_dir}"
    RJ_ML_PYTHON="$ML_PYTHON" \
    RJ_AUAU_TIGHT_BDT_TRAIN_SCRIPT="$TRAIN_SCRIPT" \
    RJ_AUAU_TIGHT_BDT_VALIDATE_SCRIPT="$VALIDATE_SCRIPT" \
    RJ_AUAU_TIGHT_BDT_VALIDATE_STAMP="${TAG}_d${depth}_fullscore_${VALIDATE_STAMP}" \
    RJ_AUAU_TIGHT_BDT_VALIDATE_GROUP_SIZE="$VALIDATE_GROUP_SIZE" \
    RJ_AUAU_TIGHT_BDT_VALIDATE_TOTAL_SCORE_MAX_ROWS="$VALIDATE_TOTAL_SCORE_MAX" \
    RJ_AUAU_TIGHT_BDT_VALIDATE_REQUEST_MEMORY="$VALIDATE_REQUEST_MEMORY" \
    RJ_AUAU_BDT_REUSE_EXISTING_MANIFEST="$REUSE_EXISTING_MANIFEST" \
    "$PIPELINE" validateOnSimCondor \
      SOURCE="$SOURCE" \
      MODEL_DIR="$model_dir" \
      MODEL_REGISTRY="$registry" \
      OUTDIR="$report_dir" \
      groupSize "$VALIDATE_GROUP_SIZE" \
      scoreMaxRows "$VALIDATE_TOTAL_SCORE_MAX"
  done
}

case "$MODE" in
  plan) plan_mode ;;
  submit-training) submit_training ;;
  status) status_mode ;;
  render-primary-slide) render_primary_slide ;;
  submit-fullscore-check) submit_fullscore_check ;;
  fingerprint) fingerprint ;;
  -h|--help|help) usage ;;
  *) die "unknown mode: $MODE" ;;
esac
