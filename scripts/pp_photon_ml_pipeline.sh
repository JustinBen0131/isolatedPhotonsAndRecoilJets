#!/usr/bin/env bash
set -euo pipefail

# pp Photon-ID ML pipeline: PPG12-matched extraction first, then the same
# trainer/validation family used by the Au+Au photon-ID work.

REPO_BASE="${RJ_REPO_BASE:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
STAMP="${RJ_PP_PHOTON_ML_STAMP:-$(date +%Y%m%d_%H%M%S)}"
RUN_ROOT="${RJ_PP_PHOTON_ML_RUN_ROOT:-${REPO_BASE}/dataOutput/ppPhotonMLPipeline/ppg12_matched_${STAMP}}"
REMOTE_DEST_ROOT="${RJ_PP_PHOTON_ML_DEST_ROOT:-/sphenix/tg/tg01/bulk/jbennett/thesisAna/ppPhotonMLPipeline/${STAMP}}"
ML_PYTHON="${RJ_ML_PYTHON:-/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python}"
TREE_NAME="${RJ_PP_PHOTON_ML_TREE:-AuAuPhotonIDTrainingTree}"
MANIFEST="${MANIFEST:-${RUN_ROOT}/training_roots.list}"
BDT_OUTDIR="${BDT_OUTDIR:-${RUN_ROOT}/models/bdt_ppg12_sixpack}"
MLP_OUTDIR="${MLP_OUTDIR:-${RUN_ROOT}/models/mlp_ppg12_sixpack}"
STACK_NOISO_OUTDIR="${STACK_NOISO_OUTDIR:-${RUN_ROOT}/models/oof_stack_ppg12_noIso}"
STACK_ISO_OUTDIR="${STACK_ISO_OUTDIR:-${RUN_ROOT}/models/oof_stack_ppg12_iso}"
VALIDATION_OUTDIR="${VALIDATION_OUTDIR:-${RUN_ROOT}/validation}"
BDT_VALIDATION_OUTDIR="${BDT_VALIDATION_OUTDIR:-${VALIDATION_OUTDIR}/bdt}"
MLP_VALIDATION_OUTDIR="${MLP_VALIDATION_OUTDIR:-${VALIDATION_OUTDIR}/mlp}"
GROUP_SIZE="${GROUP_SIZE:-10}"
LOCAL_EVENTS="${LOCAL_EVENTS:-5000}"
VERBOSE="${VERBOSE:-1}"
PHOTON_ID_ROW_MATCH="${RJ_PHOTON_ID_ROW_MATCH:-preselectionReference_tightReference}"

PPG12_BASE_V1E_FEATURES="cluster_Et,cluster_weta_cogx,vertexz,cluster_Eta,e11_over_e33,cluster_et1,cluster_et2,cluster_et3,cluster_et4"
PPG12_PT_BINS="6,10,15,20,25,35"
PP_BDT_MAX_LOAD_ROWS_PER_CLASS="${PP_BDT_MAX_LOAD_ROWS_PER_CLASS:-2000000}"
PP_MLP_MAX_LOAD_ROWS_PER_CLASS="${PP_MLP_MAX_LOAD_ROWS_PER_CLASS:-1500000}"
PP_VALIDATION_MAX_LOAD_ROWS_PER_CLASS="${PP_VALIDATION_MAX_LOAD_ROWS_PER_CLASS:-1000000}"
PP_LOAD_SAMPLE_SEED="${PP_LOAD_SAMPLE_SEED:-42}"
PP_BDT_N_JOBS="${PP_BDT_N_JOBS:-2}"
PP_SKIP_TMVA_EXPORT="${PP_SKIP_TMVA_EXPORT:-1}"

PHOTON_SAMPLES=(run28_photonjet5 run28_photonjet10 run28_photonjet20)
JET_SAMPLES=(run28_jet5 run28_jet12 run28_jet20 run28_jet30 run28_jet40)

log() { printf '[ppPhotonML] %s\n' "$*"; }
die() { printf '[ppPhotonML][ERROR] %s\n' "$*" >&2; exit 2; }
need_file() { [[ -s "$1" ]] || die "missing required file: $1"; }
mkdir_run() { mkdir -p "$RUN_ROOT" "$BDT_OUTDIR" "$MLP_OUTDIR" "$STACK_NOISO_OUTDIR" "$STACK_ISO_OUTDIR" "$BDT_VALIDATION_OUTDIR" "$MLP_VALIDATION_OUTDIR"; }

usage() {
  cat <<EOF
Usage: $0 <mode>

Modes
  localTest                 Small photon5 + jet5 pp extraction smoke test.
  condorExtract             Submit full pp table extraction for photon5/10/20 and jet5/12/20/30/40.
  buildManifest ROOT_DIR=... Build ${TREE_NAME} ROOT manifest from extracted files.
  trainBDTFromExtraction    Train PPG12-like no-iso and iso BDTs from MANIFEST.
  trainMLPFromExtraction    Train matching no-iso and iso MLPs from MANIFEST.
  validateTables            Score held pp tables and write stack-compatible score caches.
  trainStackFromExtraction  Train no-iso and iso BDT+MLP NN stacks after score caches exist.
  validateOnSim             Alias for validateTables.
  runAll                    buildManifest, train BDT/MLP, validate tables, train stacks.
  status                    Print resolved paths and expected inputs.

Important variables
  RUN_ROOT=$RUN_ROOT
  RJ_PP_PHOTON_ML_DEST_ROOT=$REMOTE_DEST_ROOT
  MANIFEST=$MANIFEST
  ROOT_DIR=<directory with extracted ROOT files> for buildManifest
  RJ_DO_RUN=1 is required for condorExtract submissions.
EOF
}

write_metadata() {
  mkdir_run
  cat >"${RUN_ROOT}/pp_photon_ml_manifest.json" <<EOF
{
  "schema": "RJ_PP_PHOTON_ML_PIPELINE_V1",
  "run_root": "${RUN_ROOT}",
  "remote_dest_root": "${REMOTE_DEST_ROOT}",
  "tree": "${TREE_NAME}",
  "signal_samples": ["${PHOTON_SAMPLES[0]}", "${PHOTON_SAMPLES[1]}", "${PHOTON_SAMPLES[2]}"],
  "background_samples": ["${JET_SAMPLES[0]}", "${JET_SAMPLES[1]}", "${JET_SAMPLES[2]}", "${JET_SAMPLES[3]}", "${JET_SAMPLES[4]}"],
  "features_ppg12_base_v1E_noIso": "${PPG12_BASE_V1E_FEATURES}",
  "features_ppg12_base_v1E_iso_additions": "reco_eiso_clip30,reco_eiso_over_cluster_Et,reco_eiso_signed_log1p",
  "pt_bins": "${PPG12_PT_BINS}",
  "bdt_max_load_rows_per_class": ${PP_BDT_MAX_LOAD_ROWS_PER_CLASS},
  "mlp_max_load_rows_per_class": ${PP_MLP_MAX_LOAD_ROWS_PER_CLASS},
  "validation_max_load_rows_per_class": ${PP_VALIDATION_MAX_LOAD_ROWS_PER_CLASS},
  "load_sample_seed": ${PP_LOAD_SAMPLE_SEED},
  "skip_tmva_export": ${PP_SKIP_TMVA_EXPORT},
  "notes": [
    "Compatibility tree is named AuAuPhotonIDTrainingTree so existing trainers can be reused.",
    "pp centrality branch is filled with -1 and is not included in pp feature presets.",
    "PPG12 parity mode uses row-level split seed 42 and raw event/sample weights off.",
    "photon-ID row match: ${PHOTON_ID_ROW_MATCH}"
  ]
}
EOF
}

extract_one() {
  local dataset="$1"
  local sample="$2"
  local role="$3"
  local mode="$4"
  shift 4

  log "extract ${dataset} sample=${sample} role=${role} mode=${mode}"
  (
    cd "$REPO_BASE"
    local submitter="${RECOILJETS_SUBMIT:-./RecoilJets_Condor_submit.sh}"
    if [[ ! -x "$submitter" && -x "./scripts/RecoilJets_Condor_submit.sh" ]]; then
      submitter="./scripts/RecoilJets_Condor_submit.sh"
    fi
    [[ -x "$submitter" ]] || die "RecoilJets submitter not executable: $submitter"
    local dest_base="${RJ_PP_PHOTON_ML_DEST_ROOT:-$REMOTE_DEST_ROOT}/${dataset}"
    RJ_PP_PHOTONID_EXTRACT_ONLY=1 \
    RJ_PP_PHOTONID_TRAINING_TREE=1 \
    RJ_PP_PHOTONID_TRAINING_TREE_MAX_ENTRIES="${RJ_PP_PHOTONID_TRAINING_TREE_MAX_ENTRIES:-0}" \
    RJ_PP_PHOTONID_SOURCE_ROLE="$role" \
    RJ_PP_PHOTONID_PPG12_FILTER=1 \
    RJ_PHOTON_ID_ROW_MATCH="$PHOTON_ID_ROW_MATCH" \
    RJ_DEST_BASE_OVERRIDE="$dest_base" \
    "$submitter" "$dataset" "$mode" "$@" SAMPLE="$sample" VERBOSE="$VERBOSE"
  )
}

local_test() {
  mkdir_run
  write_metadata
  extract_one isSim run28_photonjet5 signal local "$LOCAL_EVENTS"
  extract_one isSimInclusive run28_jet5 background local "$LOCAL_EVENTS"
  log "local test submitted/ran; inspect normal RecoilJets local outputs, then build a manifest with ROOT_DIR=..."
}

condor_extract() {
  [[ "${RJ_DO_RUN:-0}" == "1" ]] || die "condorExtract is mutating; rerun with RJ_DO_RUN=1 after checking queue pressure."
  mkdir_run
  write_metadata
  log "submitting photon signal samples: ${PHOTON_SAMPLES[*]}"
  for sample in "${PHOTON_SAMPLES[@]}"; do
    extract_one isSim "$sample" signal condorDoAll groupSize "$GROUP_SIZE"
  done
  log "submitting inclusive-jet background samples: ${JET_SAMPLES[*]}"
  for sample in "${JET_SAMPLES[@]}"; do
    extract_one isSimInclusive "$sample" background condorDoAll groupSize "$GROUP_SIZE"
  done
  log "submitted extraction jobs. After outputs are pulled/merged, run: ROOT_DIR=/path/to/roots $0 buildManifest"
}

build_manifest() {
  mkdir_run
  local root_dir="${ROOT_DIR:-}"
  [[ -n "$root_dir" ]] || die "buildManifest requires ROOT_DIR=/path/to/extracted/root/files"
  [[ -d "$root_dir" ]] || die "ROOT_DIR is not a directory: $root_dir"
  local all_roots="${RUN_ROOT}/all_candidate_roots.list"
  local qa_json="${RUN_ROOT}/manifest_tree_qa.json"
  find "$root_dir" -type f -name '*.root' | sort > "$all_roots"
  local n_all
  n_all="$(wc -l < "$all_roots" | tr -d ' ')"
  [[ "$n_all" -gt 0 ]] || die "no ROOT files found under $root_dir"
  "$ML_PYTHON" - "$all_roots" "$MANIFEST" "$qa_json" "$TREE_NAME" <<'PY'
import json
import sys
from pathlib import Path

import uproot

all_roots = Path(sys.argv[1])
manifest = Path(sys.argv[2])
qa_json = Path(sys.argv[3])
tree_name = sys.argv[4]

roots = [line.strip() for line in all_roots.read_text().splitlines() if line.strip()]
good = []
missing = []
bad = []
entries_total = 0
signal_entries = 0
background_entries = 0

for path in roots:
    try:
        with uproot.open(path) as f:
            if tree_name not in f:
                missing.append(path)
                continue
            tree = f[tree_name]
            entries = int(tree.num_entries)
            entries_total += entries
            if entries:
                labels = tree["is_signal"].array(library="np")
                signal_entries += int((labels == 1).sum())
                background_entries += int((labels == 0).sum())
            good.append(path)
    except Exception as exc:
        bad.append({"path": path, "error": str(exc)})

manifest.parent.mkdir(parents=True, exist_ok=True)
manifest.write_text("\n".join(good) + ("\n" if good else ""))
qa = {
    "tree_name": tree_name,
    "candidate_root_files": len(roots),
    "files_with_tree": len(good),
    "files_missing_tree": len(missing),
    "unreadable_files": len(bad),
    "tree_entries_total": entries_total,
    "signal_entries": signal_entries,
    "background_entries": background_entries,
    "missing_tree_examples": missing[:20],
    "unreadable_examples": bad[:20],
}
qa_json.parent.mkdir(parents=True, exist_ok=True)
qa_json.write_text(json.dumps(qa, indent=2, sort_keys=True) + "\n")
print(json.dumps(qa, sort_keys=True))
if not good:
    raise SystemExit(f"no ROOT files contained tree {tree_name}; see {qa_json}")
if entries_total <= 0:
    raise SystemExit(f"tree {tree_name} exists but has zero total entries; see {qa_json}")
if signal_entries <= 0 or background_entries <= 0:
    raise SystemExit(
        f"tree {tree_name} does not contain both labels "
        f"(signal={signal_entries}, background={background_entries}); see {qa_json}"
    )
PY
  local n
  n="$(wc -l < "$MANIFEST" | tr -d ' ')"
  write_metadata
  log "manifest=$MANIFEST files_with_tree=$n candidate_files=$n_all qa=$qa_json"
}

train_bdt() {
  mkdir_run
  need_file "$MANIFEST"
  local tmva_args=()
  if [[ "$PP_SKIP_TMVA_EXPORT" == "1" ]]; then
    tmva_args+=(--skip-tmva-export)
  fi
  "$ML_PYTHON" "${REPO_BASE}/scripts/train_auau_photon_bdt.py" \
    --task tight \
    --input "@${MANIFEST}" \
    --tree "$TREE_NAME" \
    --outdir "$BDT_OUTDIR" \
    --campaign ppg12-sixpack \
    --pt-bins "$PPG12_PT_BINS" \
    --cent-bins=-1:0 \
    --test-size 0.20 \
    --random-seed 42 \
    --n-estimators 750 \
    --max-depth 5 \
    --learning-rate 0.1 \
    --subsample 0.5 \
    --colsample-bytree 0.6 \
    --tree-method hist \
    --reg-alpha 5.0 \
    --reg-lambda 0.3 \
    --grow-policy lossguide \
    --max-bin 256 \
    --n-jobs "$PP_BDT_N_JOBS" \
    --no-event-weight \
    --et-reweight \
    --eta-reweight \
    --max-load-rows-per-class "$PP_BDT_MAX_LOAD_ROWS_PER_CLASS" \
    --load-sample-seed "$PP_LOAD_SAMPLE_SEED" \
    --majority-cap-ratio 0 \
    --skip-missing-tree \
    "${tmva_args[@]}" \
    --registry-output "${BDT_OUTDIR}/model_registry.json"
  log "BDT done: ${BDT_OUTDIR}"
}

train_mlp() {
  mkdir_run
  need_file "$MANIFEST"
  "$ML_PYTHON" "${REPO_BASE}/scripts/train_auau_photon_mlp.py" \
    --input "@${MANIFEST}" \
    --tree "$TREE_NAME" \
    --outdir "$MLP_OUTDIR" \
    --products ppg12-sixpack \
    --pt-range 6:35 \
    --centrality-range=-1:0 \
    --train-pt-bins "$PPG12_PT_BINS" \
    --validation-fraction 0.10 \
    --test-fraction 0.10 \
    --random-seed 137 \
    --hidden-layer-grid "256,128,64;384,192,96,48" \
    --epochs 280 \
    --patience 55 \
    --batch-size 8192 \
    --learning-rate 5.0e-4 \
    --l2 7.0e-5 \
    --no-event-weight \
    --et-reweight \
    --eta-reweight \
    --max-load-rows-per-class "$PP_MLP_MAX_LOAD_ROWS_PER_CLASS" \
    --load-sample-seed "$PP_LOAD_SAMPLE_SEED" \
    --selection-metric validation_auc \
    --skip-missing-tree \
    --registry-output "${MLP_OUTDIR}/model_registry.json"
  log "MLP done: ${MLP_OUTDIR}"
}

validate_tables() {
  mkdir_run
  need_file "$MANIFEST"
  need_file "${BDT_OUTDIR}/model_registry.json"
  need_file "${MLP_OUTDIR}/model_registry.json"
  "$ML_PYTHON" "${REPO_BASE}/scripts/validate_pp_photon_ml_tables.py" \
    --input "@${MANIFEST}" \
    --tree "$TREE_NAME" \
    --outdir "$BDT_VALIDATION_OUTDIR" \
    --kind bdt \
    --bdt-registry "${BDT_OUTDIR}/model_registry.json" \
    --pt-range 6:35 \
    --centrality-range=-1:0 \
    --pt-bins "$PPG12_PT_BINS" \
    --max-load-rows-per-class "$PP_VALIDATION_MAX_LOAD_ROWS_PER_CLASS" \
    --random-seed "$PP_LOAD_SAMPLE_SEED" \
    --skip-missing-tree
  "$ML_PYTHON" "${REPO_BASE}/scripts/validate_pp_photon_ml_tables.py" \
    --input "@${MANIFEST}" \
    --tree "$TREE_NAME" \
    --outdir "$MLP_VALIDATION_OUTDIR" \
    --kind mlp \
    --mlp-registry "${MLP_OUTDIR}/model_registry.json" \
    --pt-range 6:35 \
    --centrality-range=-1:0 \
    --pt-bins "$PPG12_PT_BINS" \
    --max-load-rows-per-class "$PP_VALIDATION_MAX_LOAD_ROWS_PER_CLASS" \
    --random-seed "$PP_LOAD_SAMPLE_SEED" \
    --skip-missing-tree
  log "validation done: ${VALIDATION_OUTDIR}"
}

train_stack_one() {
  local tag="$1"
  local outdir="$2"
  local mlp_score="$3"
  local bdt_score="$4"
  local extra_args=()
  if [[ "$tag" == "iso" ]]; then
    extra_args+=(--include-isolation-context)
  fi
  mkdir_run
  local mlp_cache="${MLP_CACHE:-${MLP_VALIDATION_OUTDIR}/score_caches.list}"
  local bdt_cache="${BDT_CACHE:-${BDT_VALIDATION_OUTDIR}/score_caches.list}"
  need_file "$mlp_cache"
  need_file "$bdt_cache"
  "$ML_PYTHON" "${REPO_BASE}/scripts/train_auau_oof_residual_superstacker.py" \
    --mlp-cache "$mlp_cache" \
    --bdt-cache "$bdt_cache" \
    --outdir "$outdir" \
    --mlp-score "$mlp_score" \
    --bdt-score "$bdt_score" \
    --require-full-stat \
    --expected-shards 1 \
    --pt-min 6 \
    --pt-max 35 \
    --cent-min -1 \
    --cent-max 0 \
    --report-pt-bins "$PPG12_PT_BINS" \
    --report-cent-bins=-1,0 \
    --folds 5 \
    --lower-algorithms logistic,gbm,nn \
    --lower-feature-mode full_features \
    --super-feature-mode scores_plus_full_features \
    --final-mode direct \
    --model-name "ppg12_base_v1E_oofSuperNN_${tag}" \
    "${extra_args[@]}"
  log "stack ${tag} done: ${outdir}"
}

train_stack() {
  train_stack_one noIso "$STACK_NOISO_OUTDIR" score_ppg12_base_v1E_mlp_noIso score_ppg12_base_v1E_bdt_noIso
  train_stack_one iso "$STACK_ISO_OUTDIR" score_ppg12_base_v1E_mlp_iso score_ppg12_base_v1E_bdt_iso
}

status() {
  write_metadata
  cat <<EOF
RUN_ROOT=$RUN_ROOT
REMOTE_DEST_ROOT=$REMOTE_DEST_ROOT
MANIFEST=$MANIFEST
BDT_OUTDIR=$BDT_OUTDIR
MLP_OUTDIR=$MLP_OUTDIR
STACK_NOISO_OUTDIR=$STACK_NOISO_OUTDIR
STACK_ISO_OUTDIR=$STACK_ISO_OUTDIR
VALIDATION_OUTDIR=$VALIDATION_OUTDIR
BDT_VALIDATION_OUTDIR=$BDT_VALIDATION_OUTDIR
MLP_VALIDATION_OUTDIR=$MLP_VALIDATION_OUTDIR
TREE_NAME=$TREE_NAME
FEATURES=$PPG12_BASE_V1E_FEATURES
PP_BDT_MAX_LOAD_ROWS_PER_CLASS=$PP_BDT_MAX_LOAD_ROWS_PER_CLASS
PP_MLP_MAX_LOAD_ROWS_PER_CLASS=$PP_MLP_MAX_LOAD_ROWS_PER_CLASS
PP_VALIDATION_MAX_LOAD_ROWS_PER_CLASS=$PP_VALIDATION_MAX_LOAD_ROWS_PER_CLASS
PP_LOAD_SAMPLE_SEED=$PP_LOAD_SAMPLE_SEED
PP_BDT_N_JOBS=$PP_BDT_N_JOBS
PP_SKIP_TMVA_EXPORT=$PP_SKIP_TMVA_EXPORT
EOF
}

mode="${1:-}"
case "$mode" in
  localTest) local_test ;;
  condorExtract) condor_extract ;;
  buildManifest) build_manifest ;;
  trainBDTFromExtraction) train_bdt ;;
  trainMLPFromExtraction) train_mlp ;;
  validateTables) validate_tables ;;
  trainStackFromExtraction) train_stack ;;
  validateOnSim) validate_tables ;;
  runAll) build_manifest; train_bdt; train_mlp; validate_tables; train_stack ;;
  status) status ;;
  ""|-h|--help|help) usage ;;
  *) usage; die "unknown mode: $mode" ;;
esac
