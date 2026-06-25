#!/usr/bin/env bash
set -euo pipefail

MODE="${1:-manifest}"
if [[ $# -gt 0 ]]; then
  shift
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
RJ_REPO_BASE="${RJ_REPO_BASE:-$(cd "${SCRIPT_DIR}/../../../.." && pwd -P)}"
PIPELINE="${RJ_THE8_BDT_PIPELINE:-${RJ_REPO_BASE}/scripts/auau_tight_bdt_pipeline.sh}"
TRAIN_SCRIPT="${RJ_THE8_TRAIN_SCRIPT:-${RJ_REPO_BASE}/scripts/ml/training/train_auau_photon_bdt.py}"
ML_PYTHON="${RJ_ML_PYTHON:-/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python}"

STAMP="${RJ_THE8_EXPANSION_STAMP:-$(date +%Y%m%d_%H%M%S)}"
DEFAULT_TAG="THE8_corrected_truthdefault_20260615"
DEFAULT_SOURCE="/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260615_002939"
DEFAULT_ROOT_MANIFEST="/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/the57_combined_training_roots.list"
DEFAULT_EVENT_QUALITY_CUT_JSON="/sphenix/u/patsfan753/scratch/thesisAnalysis/dataOutput/auauTightBDTValidation/THE58_lowCaloSmoothCut_20260614/the58_low_calo_cut_v1.json"
TAG="${RJ_THE8_EXPANSION_TAG:-${DEFAULT_TAG}}"
SOURCE="${RJ_THE8_EXPANSION_SOURCE:-${DEFAULT_SOURCE}}"
ROOT_MANIFEST="${RJ_THE8_EXPANSION_ROOT_MANIFEST:-${DEFAULT_ROOT_MANIFEST}}"
MODEL_BASE="${RJ_THE8_MODEL_BASE:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/the8_models}"
BASELINE_MODEL_DIR="${RJ_THE8_BASELINE_MODEL_DIR:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/the57_models/a1_withcut_20260615_stagedcache_streamreduce}"
EVENT_QUALITY_CUT_JSON="${RJ_THE8_EVENT_QUALITY_CUT_JSON:-${DEFAULT_EVENT_QUALITY_CUT_JSON}}"
ALLOW_NO_EVENT_QUALITY_CUT="${RJ_THE8_ALLOW_NO_EVENT_QUALITY_CUT:-0}"
LOCAL_MANIFEST_ROOT="${RJ_THE8_LOCAL_MANIFEST_ROOT:-${RJ_REPO_BASE}/condor_generated_configs/${TAG}}"

DEFAULT_SIGNAL_SAMPLES="run28_embeddedPhoton12 run28_embeddedPhoton20"
DEFAULT_BACKGROUND_SAMPLES="run28_embeddedJet12 run28_embeddedJet20 run28_embeddedJet30 run28_embeddedJet40"
DEFAULT_EXPECTED_SAMPLES="run28_embeddedPhoton12,run28_embeddedPhoton20,run28_embeddedJet12,run28_embeddedJet20,run28_embeddedJet30,run28_embeddedJet40"

WEIGHT_MODE="${RJ_THE8_WEIGHT_MODE:-ppg12-exact}"
SPLIT_MODE="${RJ_THE8_SPLIT_MODE:-row}"
TRAIN_GROUP_SIZE="${RJ_THE8_TRAIN_GROUP_SIZE:-1}"
TRAIN_REQUEST_MEMORY="${RJ_THE8_TRAIN_REQUEST_MEMORY:-24000MB}"
CACHE_REQUEST_MEMORY="${RJ_THE8_CACHE_REQUEST_MEMORY:-16000MB}"
CACHE_PART_REQUEST_MEMORY="${RJ_THE8_CACHE_PART_REQUEST_MEMORY:-6000MB}"
CACHE_REDUCE_REQUEST_MEMORY="${RJ_THE8_CACHE_REDUCE_REQUEST_MEMORY:-16000MB}"
CACHE_SHARDS="${RJ_THE8_CACHE_SHARDS:-12}"
CACHE_PART_MAXJOBS="${RJ_THE8_CACHE_PART_MAXJOBS:-4}"
TRAIN_PARALLEL_PER_JOB="${RJ_THE8_TRAIN_PARALLEL_PER_JOB:-1}"
XGB_N_JOBS="${RJ_THE8_XGB_N_JOBS:-1}"
PPG12_CLOSURE_ARTIFACTS="${RJ_THE8_PPG12_CLOSURE_ARTIFACTS:-metadata-only}"

ET_BINS="${RJ_THE8_ET_BINS:-15,17,19,21,23,25,27,30,35}"
COARSE_CENT_BINS="${RJ_THE8_COARSE_CENT_BINS:-0:20,20:50,50:80}"
FINE_CENT_BINS="${RJ_THE8_FINE_CENT_BINS:-0:10,10:20,20:30,30:40,40:50,50:60,60:80}"
DEPTHS_CSV="${RJ_THE8_DEPTHS:-7,8,9,10}"
SEEDS_CSV="${RJ_THE8_SEEDS:-101,202,303}"
LANES_CSV="${RJ_THE8_LANES:-}"

say() { printf '\033[1;36m[the8Diag]\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[the8Diag][WARN]\033[0m %s\n' "$*" >&2; }
die() { printf '\033[1;31m[the8Diag][ERR]\033[0m %s\n' "$*" >&2; exit 2; }

usage() {
  cat <<'EOF'
Usage:
  submit_the8_corrected_baseline_diagnostic_expansion.sh manifest [KEY=VALUE ...]
  submit_the8_corrected_baseline_diagnostic_expansion.sh plan [KEY=VALUE ...]
  submit_the8_corrected_baseline_diagnostic_expansion.sh submit-training [KEY=VALUE ...]
  submit_the8_corrected_baseline_diagnostic_expansion.sh status [KEY=VALUE ...]
  submit_the8_corrected_baseline_diagnostic_expansion.sh fingerprint [KEY=VALUE ...]
  submit_the8_corrected_baseline_diagnostic_expansion.sh heartbeat-prompt [KEY=VALUE ...]

Modes:
  manifest
      Write and print the suite manifest. Does not require SDCC source paths.
  plan
      Build dry-run DAG handles for every lane. Defaults to the THE-57 A1
      corrected truth baseline source, 20,003-root manifest, and THE-58 low-calo
      cut JSON; override SOURCE/ROOT_MANIFEST/EVENT_QUALITY_CUT_JSON only when
      deliberately regenerating from a newer baseline.
  submit-training
      Submit all training DAGs. Requires RJ_THE8_EXPANSION_APPROVED=1,
      RJ_THE8_DUPLICATE_GUARD_PASSED=1, RJ_CODEX_CHAT_NAME, and
      RJ_CODEX_THREAD_ID.
  status
      Read-only status summary for expected registries and matching queue rows.
  fingerprint
      Print the campaign identity and hard gates.
  heartbeat-prompt
      Write a paste-ready monitoring prompt under the manifest directory.

Common overrides:
  SOURCE=/path ROOT_MANIFEST=/path TAG=name STAMP=YYYYmmdd_HHMMSS
  EVENT_QUALITY_CUT_JSON=/path MODEL_BASE=/path
  LANES=sample_jet12_20,sample_jet12_20_30
EOF
}

tag_override=0
manifest_root_override=0
for tok in "$@"; do
  case "$tok" in
    SOURCE=*) SOURCE="${tok#SOURCE=}" ;;
    ROOT_MANIFEST=*) ROOT_MANIFEST="${tok#ROOT_MANIFEST=}" ;;
    EVENT_QUALITY_CUT_JSON=*) EVENT_QUALITY_CUT_JSON="${tok#EVENT_QUALITY_CUT_JSON=}" ;;
    TAG=*) TAG="${tok#TAG=}"; tag_override=1 ;;
    STAMP=*) STAMP="${tok#STAMP=}" ;;
    MODEL_BASE=*) MODEL_BASE="${tok#MODEL_BASE=}" ;;
    BASELINE_MODEL_DIR=*) BASELINE_MODEL_DIR="${tok#BASELINE_MODEL_DIR=}" ;;
    LOCAL_MANIFEST_ROOT=*) LOCAL_MANIFEST_ROOT="${tok#LOCAL_MANIFEST_ROOT=}"; manifest_root_override=1 ;;
    DEPTHS=*) DEPTHS_CSV="${tok#DEPTHS=}" ;;
    SEEDS=*) SEEDS_CSV="${tok#SEEDS=}" ;;
    LANES=*) LANES_CSV="${tok#LANES=}" ;;
    -h|--help|help) usage; exit 0 ;;
    *) die "Unknown argument: $tok" ;;
  esac
done

if [[ "$manifest_root_override" == "0" && -z "${RJ_THE8_LOCAL_MANIFEST_ROOT:-}" ]]; then
  LOCAL_MANIFEST_ROOT="${RJ_REPO_BASE}/condor_generated_configs/${TAG}"
fi

IFS=',' read -r -a DEPTHS <<< "$DEPTHS_CSV"
IFS=',' read -r -a SEEDS <<< "$SEEDS_CSV"

EXTRA_SHOWER_FEATURES=(
  cluster_weta35_cogx
  cluster_wphi53_cogx
  cluster_w32
  cluster_w52
  cluster_w72
  e11_over_e22
  e11_over_e13
  e11_over_e15
  e11_over_e17
  e11_over_e31
  e11_over_e51
  e11_over_e71
  e22_over_e33
  e22_over_e35
  e22_over_e37
  e22_over_e53
  cluster_weta_over_wphi
  cluster_weta33_over_wphi33
)

safe_token() {
  local text="$1"
  printf '%s' "$text" | tr -c '[:alnum:]' '_'
}

join_csv() {
  local IFS=,
  printf '%s' "$*"
}

build_feature_ladder_spec_ids() {
  local ids=(baseV3E11_pt1535 baseV3E11_cent12_pt1535)
  local feature token idx=1
  for feature in "${EXTRA_SHOWER_FEATURES[@]}"; do
    token="$(safe_token "$feature")"
    ids+=("base14_plus1_${token}_pt1535")
  done
  for feature in "${EXTRA_SHOWER_FEATURES[@]}"; do
    token="$(safe_token "$feature")"
    if [[ "$idx" -eq "${#EXTRA_SHOWER_FEATURES[@]}" ]]; then
      ids+=("globalEtCent1535_bdt_noIso")
    else
      ids+=("base14_ladder$(printf '%02d' "$idx")_${token}_pt1535")
    fi
    idx=$((idx + 1))
  done
  join_csv "${ids[@]}"
}

LANES=(
  etwindow_5to40
  feature_ladder
  binned14
  sample_jet12_20
  sample_jet12_20_30
  split_50_50
  split_10_90
  iso14
)
for depth in "${DEPTHS[@]}"; do
  LANES+=("depth_${depth}")
done
for seed in "${SEEDS[@]}"; do
  LANES+=("seed_${seed}")
done
LANES+=(hyper_lr002_n900 hyper_regstrong hyper_subsample070)
if [[ -n "${LANES_CSV//[[:space:]]/}" ]]; then
  IFS=',' read -r -a LANES <<< "$LANES_CSV"
fi

lane_config() {
  local lane="$1"
  LANE_FAMILY=""
  LANE_DESCRIPTION=""
  LANE_CAMPAIGN="expanded-tight"
  LANE_SPEC_IDS="centAsFeatBase3x3_pt15to35"
  LANE_EXTRA_BASE3X3_RANGES="15:35"
  LANE_SIGNAL_SAMPLES="$DEFAULT_SIGNAL_SAMPLES"
  LANE_BACKGROUND_SAMPLES="$DEFAULT_BACKGROUND_SAMPLES"
  LANE_EXPECTED_SAMPLES="$DEFAULT_EXPECTED_SAMPLES"
  LANE_TEST_SIZE="0.10"
  LANE_MAX_DEPTH=""
  LANE_RANDOM_SEED=""
  LANE_N_ESTIMATORS=""
  LANE_LEARNING_RATE=""
  LANE_SUBSAMPLE=""
  LANE_COLSAMPLE_BYTREE=""
  LANE_REG_ALPHA=""
  LANE_REG_LAMBDA=""
  LANE_MAX_BIN=""

  case "$lane" in
    etwindow_5to40)
      LANE_FAMILY="training-window"
      LANE_DESCRIPTION="5-40 GeV counterpart of the corrected 14-feature default"
      LANE_SPEC_IDS="centAsFeatBase3x3_pt5to40"
      LANE_EXTRA_BASE3X3_RANGES="5:40"
      ;;
    feature_ladder)
      LANE_FAMILY="feature-ladder"
      LANE_DESCRIPTION="baseV3E 11/12 controls, base14 plus-one/cumulative shower ladder, and full 32-feature endpoint"
      LANE_CAMPAIGN="corrected-baseline-shower-ladder"
      LANE_SPEC_IDS="$(build_feature_ladder_spec_ids)"
      LANE_EXTRA_BASE3X3_RANGES=""
      ;;
    binned14)
      LANE_FAMILY="routing"
      LANE_DESCRIPTION="default 14-feature model trained per ET, per centrality, and ET x centrality"
      LANE_CAMPAIGN="corrected-baseline-binned14"
      LANE_SPEC_IDS=""
      LANE_EXTRA_BASE3X3_RANGES=""
      ;;
    sample_jet12_20)
      LANE_FAMILY="sample-composition"
      LANE_DESCRIPTION="default 14-feature model with Jet12+20 background only"
      LANE_BACKGROUND_SAMPLES="run28_embeddedJet12 run28_embeddedJet20"
      LANE_EXPECTED_SAMPLES="run28_embeddedPhoton12,run28_embeddedPhoton20,run28_embeddedJet12,run28_embeddedJet20"
      ;;
    sample_jet12_20_30)
      LANE_FAMILY="sample-composition"
      LANE_DESCRIPTION="default 14-feature model with Jet12+20+30 background only"
      LANE_BACKGROUND_SAMPLES="run28_embeddedJet12 run28_embeddedJet20 run28_embeddedJet30"
      LANE_EXPECTED_SAMPLES="run28_embeddedPhoton12,run28_embeddedPhoton20,run28_embeddedJet12,run28_embeddedJet20,run28_embeddedJet30"
      ;;
    sig12_bkg12_20)
      LANE_FAMILY="sample-composition-signal"
      LANE_DESCRIPTION="Photon12-only signal with Jet12+20 background"
      LANE_SIGNAL_SAMPLES="run28_embeddedPhoton12"
      LANE_BACKGROUND_SAMPLES="run28_embeddedJet12 run28_embeddedJet20"
      LANE_EXPECTED_SAMPLES="run28_embeddedPhoton12,run28_embeddedJet12,run28_embeddedJet20"
      ;;
    sig12_bkg12_20_30)
      LANE_FAMILY="sample-composition-signal"
      LANE_DESCRIPTION="Photon12-only signal with Jet12+20+30 background"
      LANE_SIGNAL_SAMPLES="run28_embeddedPhoton12"
      LANE_BACKGROUND_SAMPLES="run28_embeddedJet12 run28_embeddedJet20 run28_embeddedJet30"
      LANE_EXPECTED_SAMPLES="run28_embeddedPhoton12,run28_embeddedJet12,run28_embeddedJet20,run28_embeddedJet30"
      ;;
    sig12_bkg12_20_30_40)
      LANE_FAMILY="sample-composition-signal"
      LANE_DESCRIPTION="Photon12-only signal with Jet12+20+30+40 background"
      LANE_SIGNAL_SAMPLES="run28_embeddedPhoton12"
      LANE_EXPECTED_SAMPLES="run28_embeddedPhoton12,run28_embeddedJet12,run28_embeddedJet20,run28_embeddedJet30,run28_embeddedJet40"
      ;;
    sig20_bkg12_20)
      LANE_FAMILY="sample-composition-signal"
      LANE_DESCRIPTION="Photon20-only signal with Jet12+20 background"
      LANE_SIGNAL_SAMPLES="run28_embeddedPhoton20"
      LANE_BACKGROUND_SAMPLES="run28_embeddedJet12 run28_embeddedJet20"
      LANE_EXPECTED_SAMPLES="run28_embeddedPhoton20,run28_embeddedJet12,run28_embeddedJet20"
      ;;
    sig20_bkg12_20_30)
      LANE_FAMILY="sample-composition-signal"
      LANE_DESCRIPTION="Photon20-only signal with Jet12+20+30 background"
      LANE_SIGNAL_SAMPLES="run28_embeddedPhoton20"
      LANE_BACKGROUND_SAMPLES="run28_embeddedJet12 run28_embeddedJet20 run28_embeddedJet30"
      LANE_EXPECTED_SAMPLES="run28_embeddedPhoton20,run28_embeddedJet12,run28_embeddedJet20,run28_embeddedJet30"
      ;;
    sig20_bkg12_20_30_40)
      LANE_FAMILY="sample-composition-signal"
      LANE_DESCRIPTION="Photon20-only signal with Jet12+20+30+40 background"
      LANE_SIGNAL_SAMPLES="run28_embeddedPhoton20"
      LANE_EXPECTED_SAMPLES="run28_embeddedPhoton20,run28_embeddedJet12,run28_embeddedJet20,run28_embeddedJet30,run28_embeddedJet40"
      ;;
    split_50_50)
      LANE_FAMILY="split-control"
      LANE_DESCRIPTION="default 14-feature model with 50/50 train/holdout split"
      LANE_TEST_SIZE="0.50"
      ;;
    split_10_90)
      LANE_FAMILY="split-control"
      LANE_DESCRIPTION="default 14-feature model with 10/90 train/holdout split"
      LANE_TEST_SIZE="0.90"
      ;;
    iso14)
      LANE_FAMILY="isolation-diagnostic"
      LANE_DESCRIPTION="base14 plus raw ETiso R=0.3, R=0.4, and both; diagnostic only"
      LANE_CAMPAIGN="corrected-baseline-iso14"
      LANE_SPEC_IDS=""
      LANE_EXTRA_BASE3X3_RANGES=""
      ;;
    depth_*)
      local depth="${lane#depth_}"
      LANE_FAMILY="tree-depth"
      LANE_DESCRIPTION="default 14-feature model with max_depth=${depth}"
      LANE_MAX_DEPTH="$depth"
      ;;
    seed_*)
      local seed="${lane#seed_}"
      LANE_FAMILY="seed-replica"
      LANE_DESCRIPTION="default 14-feature model with random_seed=${seed}"
      LANE_RANDOM_SEED="$seed"
      ;;
    hyper_lr002_n900)
      LANE_FAMILY="hyperparameter"
      LANE_DESCRIPTION="lower learning rate with more estimators: eta=0.02, n_estimators=900"
      LANE_LEARNING_RATE="0.02"
      LANE_N_ESTIMATORS="900"
      ;;
    hyper_regstrong)
      LANE_FAMILY="hyperparameter"
      LANE_DESCRIPTION="stronger regularization: reg_alpha=10, reg_lambda=1"
      LANE_REG_ALPHA="10.0"
      LANE_REG_LAMBDA="1.0"
      ;;
    hyper_subsample070)
      LANE_FAMILY="hyperparameter"
      LANE_DESCRIPTION="stronger row/column subsampling: subsample=0.70, colsample_bytree=0.70"
      LANE_SUBSAMPLE="0.70"
      LANE_COLSAMPLE_BYTREE="0.70"
      ;;
    *)
      die "unknown lane: $lane"
      ;;
  esac
}

model_dir_for_lane() {
  printf '%s/%s_%s' "$MODEL_BASE" "$TAG" "$1"
}

report_dir_for_lane() {
  if [[ -n "${RJ_THE8_REPORT_BASE:-}" ]]; then
    printf '%s/%s_%s' "$RJ_THE8_REPORT_BASE" "$TAG" "$1"
  else
    printf '%s/reports/%s_%s' "$SOURCE" "$TAG" "$1"
  fi
}

need_file() {
  [[ -s "$1" ]] || die "missing required file: $1"
}

need_dir() {
  [[ -d "$1" ]] || die "missing required directory: $1"
}

path_guard() {
  need_file "$PIPELINE"
  need_file "$TRAIN_SCRIPT"
  [[ -n "$SOURCE" ]] || die "$MODE requires SOURCE=/path or RJ_THE8_EXPANSION_SOURCE"
  need_dir "$SOURCE"
  if [[ -n "$ROOT_MANIFEST" ]]; then
    need_file "$ROOT_MANIFEST"
    local nroots
    nroots="$(wc -l < "$ROOT_MANIFEST" | tr -d ' ')"
    [[ "$nroots" =~ ^[0-9]+$ && "$nroots" -gt 100 ]] || die "ROOT_MANIFEST has too few roots: $ROOT_MANIFEST lines=$nroots"
  fi
  if [[ "$ALLOW_NO_EVENT_QUALITY_CUT" != "1" ]]; then
    [[ -n "$EVENT_QUALITY_CUT_JSON" ]] || die "$MODE requires EVENT_QUALITY_CUT_JSON=/path or RJ_THE8_EVENT_QUALITY_CUT_JSON"
    need_file "$EVENT_QUALITY_CUT_JSON"
  fi
}

write_lane_root_manifest() {
  local lane="$1"
  local expected_samples="$2"
  local out="${LOCAL_MANIFEST_ROOT}/root_manifests/${TAG}_${lane}_training_roots.list"
  [[ -n "$ROOT_MANIFEST" ]] || die "write_lane_root_manifest requires ROOT_MANIFEST"
  mkdir -p "$(dirname "$out")"
  python3 - "$ROOT_MANIFEST" "$out" "$expected_samples" <<'PY'
import re
import sys
from pathlib import Path

source = Path(sys.argv[1])
out = Path(sys.argv[2])
expected = tuple(item.strip() for item in sys.argv[3].split(",") if item.strip())
if not expected:
    raise SystemExit("no expected samples supplied for lane ROOT manifest")

patterns = {
    sample: re.compile(r"(^|/)" + re.escape(sample) + r"([/_]|$)")
    for sample in expected
}
selected = []
counts = {sample: 0 for sample in expected}
for raw in source.read_text().splitlines():
    line = raw.strip()
    if not line:
        continue
    for sample, pattern in patterns.items():
        if pattern.search(line):
            selected.append(line)
            counts[sample] += 1
            break

missing = [sample for sample, count in counts.items() if count <= 0]
if missing or not selected:
    raise SystemExit(
        "lane ROOT manifest filter failed: "
        f"missing={missing} selected={len(selected)} source={source}"
    )
out.write_text("\n".join(selected) + "\n")
print(f"lane_root_manifest={out} roots={len(selected)} sample_counts={counts}", file=sys.stderr)
PY
  [[ -s "$out" ]] || die "built empty lane ROOT manifest: $out"
  printf '%s\n' "$out"
}

require_submit_approval() {
  [[ "${RJ_THE8_EXPANSION_APPROVED:-0}" == "1" ]] || die "$MODE requires RJ_THE8_EXPANSION_APPROVED=1"
  [[ "${RJ_THE8_DUPLICATE_GUARD_PASSED:-0}" == "1" ]] || die "$MODE requires RJ_THE8_DUPLICATE_GUARD_PASSED=1"
  [[ -n "${RJ_CODEX_CHAT_NAME:-}" ]] || die "$MODE requires RJ_CODEX_CHAT_NAME"
  [[ -n "${RJ_CODEX_THREAD_ID:-}" ]] || die "$MODE requires RJ_CODEX_THREAD_ID"
}

write_manifest() {
  mkdir -p "$LOCAL_MANIFEST_ROOT"
  local manifest="${LOCAL_MANIFEST_ROOT}/the8_diagnostic_training_suite.tsv"
  {
    printf 'lane\tfamily\tcampaign\tspec_ids\ttest_size\tmax_depth\trandom_seed\tsamples\tdescription\tmodel_dir\n'
    printf 'reference_existing\tbaseline-reference\tno-retrain\tcentAsFeatBase3x3_pt15to35\t0.10\t\t\t%s\tExisting corrected THE-57 default baseline; reused unless deliberately rerun\t%s\n' "$DEFAULT_EXPECTED_SAMPLES" "$BASELINE_MODEL_DIR"
    local lane sample_text
    for lane in "${LANES[@]}"; do
      lane_config "$lane"
      sample_text="$LANE_EXPECTED_SAMPLES"
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$lane" "$LANE_FAMILY" "$LANE_CAMPAIGN" "${LANE_SPEC_IDS:-<all>}" \
        "$LANE_TEST_SIZE" "$LANE_MAX_DEPTH" "$LANE_RANDOM_SEED" "$sample_text" \
        "$LANE_DESCRIPTION" "$(model_dir_for_lane "$lane")"
    done
  } > "$manifest"
  say "manifest=$manifest"
  column -t -s $'\t' "$manifest" || cat "$manifest"
}

print_plan_header() {
  say "schema=THE8_CORRECTED_BASELINE_DIAGNOSTIC_EXPANSION_V1"
  say "mode=$MODE"
  say "tag=$TAG"
  say "source=${SOURCE:-<required for plan/submit>}"
  say "root_manifest=${ROOT_MANIFEST:-<pipeline-generated from source>}"
  say "baseline_model_dir=$BASELINE_MODEL_DIR"
  say "model_base=$MODEL_BASE"
  say "event_quality_cut_json=${EVENT_QUALITY_CUT_JSON:-<required unless RJ_THE8_ALLOW_NO_EVENT_QUALITY_CUT=1>}"
  say "lanes=${#LANES[@]} plus reference_existing no-retrain row"
  say "staged_cache=1 cache_shards=$CACHE_SHARDS cache_part_maxjobs=$CACHE_PART_MAXJOBS"
  say "memory train=$TRAIN_REQUEST_MEMORY cache=$CACHE_REQUEST_MEMORY cache_part=$CACHE_PART_REQUEST_MEMORY cache_reduce=$CACHE_REDUCE_REQUEST_MEMORY"
  say "depths=$DEPTHS_CSV seeds=$SEEDS_CSV"
}

run_lane() {
  local lane="$1"
  local dryrun="$2"
  lane_config "$lane"
  local lane_model_dir lane_report_dir lane_root_manifest
  lane_model_dir="$(model_dir_for_lane "$lane")"
  lane_report_dir="$(report_dir_for_lane "$lane")"
  lane_root_manifest=""
  if [[ -n "$ROOT_MANIFEST" ]]; then
    lane_root_manifest="$(write_lane_root_manifest "$lane" "$LANE_EXPECTED_SAMPLES")"
  fi
  say "lane=$lane family=$LANE_FAMILY campaign=$LANE_CAMPAIGN model_dir=$lane_model_dir dryrun=$dryrun"
  (
    export ML_PYTHON="$ML_PYTHON"
    export RJ_AUAU_TIGHT_BDT_TRAIN_SCRIPT="$TRAIN_SCRIPT"
    export RJ_AUAU_TIGHT_BDT_TRAIN_STAMP="${TAG}_${lane}"
    export RJ_AUAU_TIGHT_BDT_MODEL_DIR="$lane_model_dir"
    export RJ_AUAU_BDT_REPORT_DIR="$lane_report_dir"
    export RJ_AUAU_BDT_CAMPAIGN="$LANE_CAMPAIGN"
    export RJ_AUAU_BDT_CAMPAIGN_SPEC_IDS="$LANE_SPEC_IDS"
    export RJ_AUAU_BDT_WEIGHT_MODE="$WEIGHT_MODE"
    export RJ_AUAU_BDT_TEST_SIZE="$LANE_TEST_SIZE"
    export RJ_AUAU_BDT_SPLIT_MODE="$SPLIT_MODE"
    export RJ_AUAU_BDT_EXTRA_CENT_AS_FEAT_BASE3X3_PT_RANGES="$LANE_EXTRA_BASE3X3_RANGES"
    export RJ_AUAU_BDT_ETFINE_PT_BINS="$ET_BINS"
    export RJ_AUAU_BDT_ETFINE_COARSE_CENT_BINS="$COARSE_CENT_BINS"
    export RJ_AUAU_BDT_ETFINE_FINE_CENT_BINS="$FINE_CENT_BINS"
    export RJ_AUAU_BDT_PPG12_EXACT_EXPECTED_SAMPLES="$LANE_EXPECTED_SAMPLES"
    export RJ_AUAU_BDT_PPG12_EXACT_CLOSURE_ARTIFACTS="$PPG12_CLOSURE_ARTIFACTS"
    export RJ_AUAU_TIGHT_BDT_SIGNAL_SAMPLES="$LANE_SIGNAL_SAMPLES"
    export RJ_AUAU_TIGHT_BDT_BACKGROUND_SAMPLES="$LANE_BACKGROUND_SAMPLES"
    export RJ_AUAU_BDT_STAGED_CACHE=1
    export RJ_AUAU_BDT_EXPANDED_GROUP_SIZE="$TRAIN_GROUP_SIZE"
    export RJ_AUAU_BDT_EXPANDED_REQUEST_MEMORY="$TRAIN_REQUEST_MEMORY"
    export RJ_AUAU_BDT_CACHE_REQUEST_MEMORY="$CACHE_REQUEST_MEMORY"
    export RJ_AUAU_BDT_CACHE_PART_REQUEST_MEMORY="$CACHE_PART_REQUEST_MEMORY"
    export RJ_AUAU_BDT_CACHE_REDUCE_REQUEST_MEMORY="$CACHE_REDUCE_REQUEST_MEMORY"
    export RJ_AUAU_BDT_CACHE_SHARDS="$CACHE_SHARDS"
    export RJ_AUAU_BDT_CACHE_PART_MAXJOBS="$CACHE_PART_MAXJOBS"
    export RJ_AUAU_BDT_TRAIN_PARALLEL_PER_JOB="$TRAIN_PARALLEL_PER_JOB"
    export RJ_AUAU_BDT_XGB_N_JOBS="$XGB_N_JOBS"
    export RJ_DAG_DRYRUN="$dryrun"
    [[ -z "$lane_root_manifest" ]] || export RJ_AUAU_BDT_ROOT_MANIFEST="$lane_root_manifest"
    [[ -z "$EVENT_QUALITY_CUT_JSON" ]] || export RJ_AUAU_BDT_EVENT_QUALITY_CUT_JSON="$EVENT_QUALITY_CUT_JSON"
    [[ -z "$LANE_MAX_DEPTH" ]] || export RJ_AUAU_BDT_MAX_DEPTH="$LANE_MAX_DEPTH"
    [[ -z "$LANE_RANDOM_SEED" ]] || export RJ_AUAU_BDT_RANDOM_SEED="$LANE_RANDOM_SEED"
    [[ -z "$LANE_N_ESTIMATORS" ]] || export RJ_AUAU_BDT_N_ESTIMATORS="$LANE_N_ESTIMATORS"
    [[ -z "$LANE_LEARNING_RATE" ]] || export RJ_AUAU_BDT_LEARNING_RATE="$LANE_LEARNING_RATE"
    [[ -z "$LANE_SUBSAMPLE" ]] || export RJ_AUAU_BDT_SUBSAMPLE="$LANE_SUBSAMPLE"
    [[ -z "$LANE_COLSAMPLE_BYTREE" ]] || export RJ_AUAU_BDT_COLSAMPLE_BYTREE="$LANE_COLSAMPLE_BYTREE"
    [[ -z "$LANE_REG_ALPHA" ]] || export RJ_AUAU_BDT_REG_ALPHA="$LANE_REG_ALPHA"
    [[ -z "$LANE_REG_LAMBDA" ]] || export RJ_AUAU_BDT_REG_LAMBDA="$LANE_REG_LAMBDA"
    [[ -z "$LANE_MAX_BIN" ]] || export RJ_AUAU_BDT_MAX_BIN="$LANE_MAX_BIN"
    "$PIPELINE" trainExpandedFromExtractionCondor "SOURCE=$SOURCE" groupSize "$TRAIN_GROUP_SIZE"
  )
}

run_all_lanes() {
  local dryrun="$1"
  local lane
  for lane in "${LANES[@]}"; do
    run_lane "$lane" "$dryrun"
  done
}

status_summary() {
  print_plan_header
  local lane registry status
  for lane in "${LANES[@]}"; do
    registry="$(model_dir_for_lane "$lane")/model_registry.json"
    if [[ -s "$registry" ]]; then
      status="$("$ML_PYTHON" - "$registry" <<'PY'
import json, sys
data = json.load(open(sys.argv[1]))
print(f"{data.get('status')} trained={data.get('trained_model_count', data.get('model_count'))} expected={data.get('expected_model_count')}")
PY
)"
    else
      status="MISSING"
    fi
    say "lane=$lane registry=$registry status=$status"
  done
  if command -v condor_q >/dev/null 2>&1; then
    say "matching_condor_rows:"
    condor_q "${USER:-$LOGNAME}" -wide 2>/dev/null | grep -F "$TAG" || true
  else
    warn "condor_q not available on this host; registry status only."
  fi
}

fingerprint() {
  cat <<EOF
purpose=THE-8 corrected-baseline diagnostic training expansion
linear_issue=THE-75
parent=THE-8
tag=${TAG}
source=${SOURCE:-<required>}
root_manifest=${ROOT_MANIFEST:-<pipeline-generated from source>}
baseline_model_dir=${BASELINE_MODEL_DIR}
baseline_definition=AuAu baseV3E+centrality+weta33/wphi33 14 features, 15-35 GeV, PPG12-exact ET/eta reweight, THE-58 cleaned
default_expected_samples=${DEFAULT_EXPECTED_SAMPLES}
event_quality_cut_json=${EVENT_QUALITY_CUT_JSON:-<required unless explicitly bypassed>}
lanes=${#LANES[@]}
depths=${DEPTHS_CSV}
seeds=${SEEDS_CSV}
submit_gates=RJ_THE8_EXPANSION_APPROVED,RJ_THE8_DUPLICATE_GUARD_PASSED,RJ_CODEX_CHAT_NAME,RJ_CODEX_THREAD_ID
EOF
}

heartbeat_prompt() {
  mkdir -p "$LOCAL_MANIFEST_ROOT"
  local prompt="${LOCAL_MANIFEST_ROOT}/the8_diagnostic_training_heartbeat_prompt.txt"
  cat > "$prompt" <<EOF
Read-only heartbeat for THE-75 / THE-8Q corrected-baseline diagnostic training expansion.

1. Run the wrapper status mode:
   ${SCRIPT_DIR}/submit_the8_corrected_baseline_diagnostic_expansion.sh status TAG=${TAG} SOURCE=${SOURCE:-<SOURCE>} ROOT_MANIFEST=${ROOT_MANIFEST:-<ROOT_MANIFEST>} EVENT_QUALITY_CUT_JSON=${EVENT_QUALITY_CUT_JSON:-<CUT_JSON>}
2. For each lane, report: family, model_dir, registry status, trained/skipped/missing counts, DAG/sub_root if present, held/running/idle rows, and next action.
3. Do not call the campaign valid from Condor history alone. Require non-tiny ROOT/TMVA outputs, model_registry status, expected feature counts, and finite-score smoke/validation before promotion.
4. Do not submit, remove, rerun, merge, upload, or edit remote files without a fresh Justin approval gate.
EOF
  say "heartbeat_prompt=$prompt"
}

case "$MODE" in
  manifest)
    print_plan_header
    write_manifest
    ;;
  plan)
    path_guard
    print_plan_header
    write_manifest
    run_all_lanes 1
    ;;
  submit-training)
    path_guard
    require_submit_approval
    print_plan_header
    write_manifest
    run_all_lanes 0
    ;;
  status)
    status_summary
    ;;
  fingerprint)
    fingerprint
    ;;
  heartbeat-prompt)
    heartbeat_prompt
    ;;
  *)
    usage
    die "unknown mode: $MODE"
    ;;
esac
