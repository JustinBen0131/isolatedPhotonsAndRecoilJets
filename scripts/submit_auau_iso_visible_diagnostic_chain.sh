#!/usr/bin/env bash
set -euo pipefail

RJ_REPO_BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
readonly RJ_REPO_BASE

ML_PYTHON="${RJ_ML_PYTHON:-/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python}"
SOURCE="${RJ_AUAU_ISO_DIAG_SOURCE:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049}"
MODEL_BASE="${RJ_AUAU_MLP_MODEL_BASE:-/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models}"
BDT_MODEL_BASE="${RJ_AUAU_BDT_MODEL_BASE:-/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models}"
STAMP="${RJ_AUAU_ISO_DIAG_STAMP:-$(date +%Y%m%d_%H%M%S)}"
CHAIN_ROOT="${RJ_AUAU_ISO_DIAG_CHAIN_ROOT:-${MODEL_BASE}/isolation_visible_diagnostic_${STAMP}}"
BDT_MODEL_DIR="${RJ_AUAU_ISO_DIAG_BDT_MODEL_DIR:-${BDT_MODEL_BASE}/tight_iso_visible_diagnostic_${STAMP}}"
MLP_MODEL_DIR="${RJ_AUAU_ISO_DIAG_MLP_MODEL_DIR:-${MODEL_BASE}/tight_mlp_iso_visible_diagnostic_${STAMP}}"
BDT_REPORT="${RJ_AUAU_ISO_DIAG_BDT_REPORT:-${SOURCE}/reports/bdt_iso_visible_validation_${STAMP}}"
MLP_REPORT="${RJ_AUAU_ISO_DIAG_MLP_REPORT:-${SOURCE}/reports/mlp_iso_visible_validation_${STAMP}}"
STACK_ISO_OUTDIR="${RJ_AUAU_ISO_DIAG_STACK_ISO_OUTDIR:-${MODEL_BASE}/nnstack_iso_scores_iso_context_${STAMP}}"
STACK_CLEAN_OUTDIR="${RJ_AUAU_ISO_DIAG_STACK_CLEAN_OUTDIR:-${MODEL_BASE}/nnstack_clean_scores_iso_context_${STAMP}}"
SUMMARY_DIR="${RJ_AUAU_ISO_DIAG_SUMMARY_DIR:-${CHAIN_ROOT}/summary}"
LOG="${RJ_AUAU_ISO_DIAG_LOG:-${CHAIN_ROOT}/isolation_visible_diagnostic_chain_${STAMP}.log}"

PT_BINS="${RJ_AUAU_ISO_DIAG_PT_BINS:-15,18,20,22.5,25,30,35}"
PT_RANGE="${RJ_AUAU_ISO_DIAG_PT_RANGE:-15:35}"
FINE_CENT_BINS="${RJ_AUAU_ISO_DIAG_FINE_CENT_BINS:-0:10,10:20,20:30,30:40,40:50,50:60,60:80}"
VALIDATE_GROUP_SIZE="${RJ_AUAU_ISO_DIAG_VALIDATE_GROUP_SIZE:-100}"
VALIDATE_TOTAL_SCORE_MAX="${RJ_AUAU_ISO_DIAG_VALIDATE_TOTAL_SCORE_MAX:-600000}"
VALIDATE_REQUEST_MEMORY="${RJ_AUAU_ISO_DIAG_VALIDATE_REQUEST_MEMORY:-3000MB}"
POLL_SECONDS="${RJ_AUAU_ISO_DIAG_POLL_SECONDS:-300}"
STACK_MAX_ROWS="${RJ_AUAU_ISO_DIAG_STACK_MAX_ROWS:-0}"
STACK_ALGORITHMS="${RJ_AUAU_ISO_DIAG_STACK_ALGORITHMS:-nn}"
STACK_VARIANTS="${RJ_AUAU_ISO_DIAG_STACK_VARIANTS:-ptFine15to35_cent7_full}"
NN_HIDDEN="${RJ_AUAU_ISO_DIAG_NN_HIDDEN:-64,32}"
NN_EPOCHS="${RJ_AUAU_ISO_DIAG_NN_EPOCHS:-180}"
NN_PATIENCE="${RJ_AUAU_ISO_DIAG_NN_PATIENCE:-24}"

CLEAN_MLP_CACHE="${RJ_AUAU_ISO_DIAG_CLEAN_MLP_CACHE:-${SOURCE}/reports/mlp_model_validation_condor_stack_full_features_20260512_2338/score_caches.list}"
CLEAN_BDT_CACHE="${RJ_AUAU_ISO_DIAG_CLEAN_BDT_CACHE:-${SOURCE}/reports/model_validation_condor_stack_full_features_20260512_2338/score_caches.list}"
CLEAN_MLP_SCORE="${RJ_AUAU_ISO_DIAG_CLEAN_MLP_SCORE:-score_centInputBase3x3WidthRatiosMLP_pt1535}"
CLEAN_BDT_SCORE="${RJ_AUAU_ISO_DIAG_CLEAN_BDT_SCORE:-score_ptFine_cent7}"

say() { printf '\033[1;36m[auauIsoVisibleDiag]\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[auauIsoVisibleDiag][WARN]\033[0m %s\n' "$*" >&2; }
err() { printf '\033[1;31m[auauIsoVisibleDiag][ERR]\033[0m %s\n' "$*" >&2; }
die() { err "$*"; exit 2; }

usage() {
  cat <<'EOF'
Usage:
  ./scripts/submit_auau_iso_visible_diagnostic_chain.sh selfTest
  ./scripts/submit_auau_iso_visible_diagnostic_chain.sh run [SOURCE=/path]
  ./scripts/submit_auau_iso_visible_diagnostic_chain.sh status [CHAIN_ROOT=/path]

This is a diagnostic-only isolation-visible photon-ID lane. It trains and
validates BDT/MLP/NN-stack variants that use reconstructed isolation-derived
inputs. It never submits isSimEmbedded/isSimEmbeddedInclusive production.
EOF
}

mode="${1:-run}"
if [[ $# -gt 0 ]]; then shift; fi

for tok in "$@"; do
  case "$tok" in
    SOURCE=*) SOURCE="${tok#SOURCE=}" ;;
    CHAIN_ROOT=*) CHAIN_ROOT="${tok#CHAIN_ROOT=}" ;;
    BDT_MODEL_DIR=*) BDT_MODEL_DIR="${tok#BDT_MODEL_DIR=}" ;;
    MLP_MODEL_DIR=*) MLP_MODEL_DIR="${tok#MLP_MODEL_DIR=}" ;;
    BDT_REPORT=*) BDT_REPORT="${tok#BDT_REPORT=}" ;;
    MLP_REPORT=*) MLP_REPORT="${tok#MLP_REPORT=}" ;;
    STACK_ISO_OUTDIR=*) STACK_ISO_OUTDIR="${tok#STACK_ISO_OUTDIR=}" ;;
    STACK_CLEAN_OUTDIR=*) STACK_CLEAN_OUTDIR="${tok#STACK_CLEAN_OUTDIR=}" ;;
    CLEAN_MLP_CACHE=*) CLEAN_MLP_CACHE="${tok#CLEAN_MLP_CACHE=}" ;;
    CLEAN_BDT_CACHE=*) CLEAN_BDT_CACHE="${tok#CLEAN_BDT_CACHE=}" ;;
    CLEAN_MLP_SCORE=*) CLEAN_MLP_SCORE="${tok#CLEAN_MLP_SCORE=}" ;;
    CLEAN_BDT_SCORE=*) CLEAN_BDT_SCORE="${tok#CLEAN_BDT_SCORE=}" ;;
    VALIDATE_TOTAL_SCORE_MAX=*) VALIDATE_TOTAL_SCORE_MAX="${tok#VALIDATE_TOTAL_SCORE_MAX=}" ;;
    STACK_MAX_ROWS=*) STACK_MAX_ROWS="${tok#STACK_MAX_ROWS=}" ;;
    STACK_ALGORITHMS=*) STACK_ALGORITHMS="${tok#STACK_ALGORITHMS=}" ;;
    STACK_VARIANTS=*) STACK_VARIANTS="${tok#STACK_VARIANTS=}" ;;
    -h|--help|help) usage; exit 0 ;;
    *) die "Unknown argument: $tok" ;;
  esac
done

print_plan() {
  say "RECOILJETS_AUAU_ISOLATION_VISIBLE_DIAGNOSTIC_CHAIN_V1"
  say "host=$(hostname -f 2>/dev/null || hostname)"
  say "stamp=$STAMP"
  say "repo=$RJ_REPO_BASE"
  say "source=$SOURCE"
  say "chain_root=$CHAIN_ROOT"
  say "bdt_model_dir=$BDT_MODEL_DIR"
  say "mlp_model_dir=$MLP_MODEL_DIR"
  say "bdt_report=$BDT_REPORT"
  say "mlp_report=$MLP_REPORT"
  say "stack_iso_outdir=$STACK_ISO_OUTDIR"
  say "stack_clean_outdir=$STACK_CLEAN_OUTDIR"
  say "summary_dir=$SUMMARY_DIR"
  say "pt_range=$PT_RANGE pt_bins=$PT_BINS fine_cent_bins=$FINE_CENT_BINS"
  say "diagnostic_warning=uses isolation-derived inputs; diagnostic ceiling test only; not ABCD-safe photon ID"
}

setup_ml_env() {
  export RJ_ML_PYTHON="$ML_PYTHON"
  local ml_python_prefix ml_python_real ml_python_real_prefix ld_joined d
  ml_python_prefix="$(cd "$(dirname "$ML_PYTHON")/.." && pwd -P 2>/dev/null || true)"
  ml_python_real="$(readlink -f "$ML_PYTHON" 2>/dev/null || true)"
  ml_python_real_prefix=""
  if [[ -n "$ml_python_real" ]]; then
    ml_python_real_prefix="$(cd "$(dirname "$ml_python_real")/.." && pwd -P 2>/dev/null || true)"
  fi
  ld_joined=""
  for d in "$ml_python_prefix/lib" "$ml_python_prefix/lib64" "$ml_python_real_prefix/lib" "$ml_python_real_prefix/lib64"; do
    [[ -n "$d" && -d "$d" ]] || continue
    case ":$ld_joined:" in *":$d:"*) ;; *) ld_joined="${ld_joined:+$ld_joined:}$d" ;; esac
  done
  [[ -n "$ld_joined" ]] && export LD_LIBRARY_PATH="$ld_joined:${LD_LIBRARY_PATH:-}"
  unset PYTHONHOME
}

make_training_manifest() {
  local search_root="$SOURCE"
  [[ -d "${SOURCE}/extraction" ]] && search_root="${SOURCE}/extraction"
  mkdir -p "$CHAIN_ROOT"
  find "$search_root" -type f -name '*.root' | sort -V > "${CHAIN_ROOT}/training_roots.list"
  [[ -s "${CHAIN_ROOT}/training_roots.list" ]] || die "No ROOT files found under $search_root"
  say "training_manifest=${CHAIN_ROOT}/training_roots.list entries=$(wc -l < "${CHAIN_ROOT}/training_roots.list" | tr -d ' ')"
}

apply_check_bdt_registry() {
  "$ML_PYTHON" - "$BDT_MODEL_DIR/model_registry.json" <<'PY'
import json
import sys
from pathlib import Path
try:
    import ROOT
except Exception as exc:
    raise SystemExit(f"PyROOT import failed: {exc}")
registry = Path(sys.argv[1])
payload = json.loads(registry.read_text())
bad = []
checked = 0
for model in payload.get("models", []):
    path = Path(model.get("output_tmva", ""))
    if not path.is_file() or path.stat().st_size <= 0:
        bad.append(f"missing {path}")
        continue
    handle = ROOT.TFile.Open(str(path))
    if not handle or handle.IsZombie():
        bad.append(f"unreadable {path}")
    else:
        checked += 1
    if handle:
        handle.Close()
if bad:
    raise SystemExit("BDT applyCheck failed:\n  " + "\n  ".join(bad[:50]))
print(f"[OK] BDT applyCheck opened {checked} TMVA ROOT files from {registry}")
PY
}

wait_for_ready() {
  local label="${1:?label}"
  local summary="${2:?summary}"
  say "waiting for ${label}: ${summary}"
  while true; do
    if [[ -s "$summary" ]]; then
      local status
      status="$(awk -F= '/^status=/ {print $2; exit}' "$summary")"
      status="${status:-CHECK}"
      cat "$summary"
      if [[ "$status" == "READY" ]]; then
        say "${label} READY"
        return 0
      fi
      die "${label} finished with status=${status}; inspect ${summary}"
    fi
    if command -v condor_q >/dev/null 2>&1; then
      condor_q -nobatch "${USER:-patsfan753}" | tail -25 || true
    fi
    sleep "$POLL_SECONDS"
  done
}

train_bdt() {
  say "Training isolation-input BDT diagnostics"
  mkdir -p "$BDT_MODEL_DIR"
  "$ML_PYTHON" scripts/train_auau_photon_bdt.py \
    --task tight \
    --campaign iso-diagnostic \
    --input "@${CHAIN_ROOT}/training_roots.list" \
    --outdir "$BDT_MODEL_DIR" \
    --cache-file "${BDT_MODEL_DIR}/training_matrix_iso_diagnostic.npz" \
    --registry-output "${BDT_MODEL_DIR}/model_registry.json" \
    --pt-bins "$PT_BINS" \
    --fine-cent-bins "$FINE_CENT_BINS" \
    --parallel-workers "${RJ_AUAU_ISO_DIAG_BDT_TRAIN_PARALLEL:-2}" \
    --n-jobs "${RJ_AUAU_ISO_DIAG_BDT_XGB_N_JOBS:-1}" \
    --majority-cap-ratio "${RJ_AUAU_ISO_DIAG_BDT_MAJORITY_CAP_RATIO:-4.0}"
  apply_check_bdt_registry
}

validate_bdt() {
  say "Submitting isolation-input BDT Condor validation"
  RJ_AUAU_TIGHT_BDT_VALIDATE_STAMP="iso_visible_${STAMP}" \
  RJ_AUAU_TIGHT_BDT_VALIDATE_REQUEST_MEMORY="$VALIDATE_REQUEST_MEMORY" \
  ./scripts/auau_tight_bdt_pipeline.sh validateOnSimCondor \
    SOURCE="$SOURCE" \
    MODEL_DIR="$BDT_MODEL_DIR" \
    MODEL_REGISTRY="${BDT_MODEL_DIR}/model_registry.json" \
    OUTDIR="$BDT_REPORT" \
    groupSize "$VALIDATE_GROUP_SIZE" \
    SCORE_MAX_ROWS="$VALIDATE_TOTAL_SCORE_MAX"
  wait_for_ready "BDT validation" "${BDT_REPORT}/validation_summary.txt"
}

train_mlp() {
  say "Training isolation-input MLP diagnostic"
  RJ_AUAU_TIGHT_MLP_PRODUCTS="iso-kitchen-sink" \
  RJ_AUAU_TIGHT_MLP_PT_RANGE="$PT_RANGE" \
  RJ_AUAU_MLP_TRAIN_PT_BINS="$PT_BINS" \
  RJ_AUAU_MLP_TRAIN_PT_BIN_WEIGHT_SPEC="15:18:0.12,18:20:0.16,20:22.5:0.20,22.5:25:0.22,25:30:0.16,30:35:0.14" \
  RJ_AUAU_MLP_TRAIN_HIGHPT_SELECTION_WEIGHTS="15:18:0.12,18:20:0.16,20:22.5:0.20,22.5:25:0.22,25:30:0.16,30:35:0.14" \
  ./scripts/auau_tight_mlp_pipeline.sh trainIsoKitchenSinkFromExtraction SOURCE="$SOURCE" MODEL_DIR="$MLP_MODEL_DIR"
  ./scripts/auau_tight_mlp_pipeline.sh applyCheck MODEL_DIR="$MLP_MODEL_DIR"
}

validate_mlp() {
  say "Submitting isolation-input MLP Condor validation"
  RJ_AUAU_TIGHT_MLP_VALIDATE_STAMP="iso_visible_${STAMP}" \
  RJ_AUAU_TIGHT_MLP_VALIDATE_REQUEST_MEMORY="$VALIDATE_REQUEST_MEMORY" \
  RJ_AUAU_TIGHT_MLP_VALIDATE_PT_BINS="$PT_BINS" \
  ./scripts/auau_tight_mlp_pipeline.sh validateOnSimCondor \
    SOURCE="$SOURCE" \
    MODEL_DIR="$MLP_MODEL_DIR" \
    OUTDIR="$MLP_REPORT" \
    groupSize "$VALIDATE_GROUP_SIZE" \
    SCORE_MAX_ROWS="$VALIDATE_TOTAL_SCORE_MAX"
  wait_for_ready "MLP validation" "${MLP_REPORT}/validation_summary.txt"
}

train_stack() {
  say "Training NN stack: iso-aware BDT score + iso-aware MLP score + isolation context"
  "$ML_PYTHON" scripts/train_auau_stacked_bdt_mlp_sweep.py \
    --mlp-cache "${MLP_REPORT}/score_caches.list" \
    --bdt-cache "${BDT_REPORT}/score_caches.list" \
    --outdir "$STACK_ISO_OUTDIR" \
    --mlp-score score_centInputKitchenSinkIsoMLP \
    --bdt-score score_isoBDT_ptFine15to35_cent7_full \
    --algorithms "$STACK_ALGORITHMS" \
    --variants "$STACK_VARIANTS" \
    --sweep full \
    --include-isolation-context \
    --include-controls \
    --max-rows "$STACK_MAX_ROWS" \
    --nn-hidden "$NN_HIDDEN" \
    --nn-epochs "$NN_EPOCHS" \
    --nn-patience "$NN_PATIENCE" \
    --note "Diagnostic-only NN stack with iso-aware BDT/MLP scores and isolation context; not ABCD-safe."

  say "Training NN stack ablation: clean BDT/MLP scores + isolation context"
  [[ -s "$CLEAN_MLP_CACHE" ]] || die "Missing CLEAN_MLP_CACHE=$CLEAN_MLP_CACHE"
  [[ -s "$CLEAN_BDT_CACHE" ]] || die "Missing CLEAN_BDT_CACHE=$CLEAN_BDT_CACHE"
  "$ML_PYTHON" scripts/train_auau_stacked_bdt_mlp_sweep.py \
    --mlp-cache "$CLEAN_MLP_CACHE" \
    --bdt-cache "$CLEAN_BDT_CACHE" \
    --outdir "$STACK_CLEAN_OUTDIR" \
    --mlp-score "$CLEAN_MLP_SCORE" \
    --bdt-score "$CLEAN_BDT_SCORE" \
    --algorithms "$STACK_ALGORITHMS" \
    --variants "$STACK_VARIANTS" \
    --sweep full \
    --include-isolation-context \
    --include-controls \
    --max-rows "$STACK_MAX_ROWS" \
    --nn-hidden "$NN_HIDDEN" \
    --nn-epochs "$NN_EPOCHS" \
    --nn-patience "$NN_PATIENCE" \
    --note "Diagnostic-only NN stack with clean scores plus isolation context; not ABCD-safe."
}

summarize() {
  "$ML_PYTHON" scripts/make_auau_iso_visible_diagnostic_summary.py \
    --bdt-report "$BDT_REPORT" \
    --mlp-report "$MLP_REPORT" \
    --iso-stack-dir "$STACK_ISO_OUTDIR" \
    --clean-context-stack-dir "$STACK_CLEAN_OUTDIR" \
    --outdir "$SUMMARY_DIR"
  cat > "${CHAIN_ROOT}/isolation_visible_diagnostic_manifest.json" <<EOF
{
  "schema": "RECOILJETS_AUAU_ISOLATION_VISIBLE_DIAGNOSTIC_CHAIN_MANIFEST_V1",
  "stamp": "${STAMP}",
  "source": "${SOURCE}",
  "diagnostic_only": true,
  "abcd_warning": "uses isolation-derived inputs; diagnostic ceiling test only; not ABCD-safe photon ID",
  "bdt_model_dir": "${BDT_MODEL_DIR}",
  "mlp_model_dir": "${MLP_MODEL_DIR}",
  "bdt_report": "${BDT_REPORT}",
  "mlp_report": "${MLP_REPORT}",
  "stack_iso_outdir": "${STACK_ISO_OUTDIR}",
  "stack_clean_outdir": "${STACK_CLEAN_OUTDIR}",
  "summary_dir": "${SUMMARY_DIR}",
  "next_action": "Inspect summary/rank plots. Do not submit isSimEmbedded/isSimEmbeddedInclusive from this diagnostic lane."
}
EOF
  say "DONE_ISO_VISIBLE_DIAGNOSTIC_CHAIN=${CHAIN_ROOT}"
}

status_mode() {
  print_plan
  for path in \
    "$BDT_MODEL_DIR/model_registry.json" \
    "$MLP_MODEL_DIR/model_registry.json" \
    "$BDT_REPORT/validation_summary.txt" \
    "$MLP_REPORT/validation_summary.txt" \
    "$STACK_ISO_OUTDIR/stacked_sweep_rank_table.csv" \
    "$STACK_CLEAN_OUTDIR/stacked_sweep_rank_table.csv" \
    "$SUMMARY_DIR/iso_visible_diagnostic_summary.txt"
  do
    if [[ -s "$path" ]]; then
      say "READY_FILE=$path"
    else
      say "missing=$path"
    fi
  done
}

run_mode() {
  print_plan
  mkdir -p "$CHAIN_ROOT" "$SUMMARY_DIR" "$(dirname "$LOG")"
  exec > >(tee -a "$LOG") 2>&1
  print_plan
  setup_ml_env
  make_training_manifest
  train_bdt
  validate_bdt
  train_mlp
  validate_mlp
  train_stack
  summarize
}

self_test_mode() {
  print_plan
  setup_ml_env
  "$ML_PYTHON" -m py_compile \
    scripts/train_auau_photon_bdt.py \
    scripts/validate_auau_tight_bdt_on_sim.py \
    scripts/train_auau_photon_mlp.py \
    scripts/validate_auau_tight_mlp_on_sim.py \
    scripts/train_auau_stacked_bdt_mlp_sweep.py \
    scripts/make_auau_iso_visible_diagnostic_summary.py
  bash -n scripts/submit_auau_iso_visible_diagnostic_chain.sh
  bash -n scripts/submit_auau_stacked_bdt_mlp_sweep.sh
  local out="${CHAIN_ROOT}/self_test"
  "$ML_PYTHON" scripts/train_auau_stacked_bdt_mlp_sweep.py \
    --outdir "$out/stack_self_test" \
    --self-test \
    --algorithms nn \
    --include-isolation-context \
    --nn-epochs 2 \
    --nn-patience 1 \
    --top-n 2
  "$ML_PYTHON" scripts/train_auau_photon_bdt.py \
    --task tight \
    --campaign iso-diagnostic \
    --outdir "$out/bdt_plan" \
    --plan-only \
    --pt-bins "$PT_BINS" \
    --fine-cent-bins "$FINE_CENT_BINS"
  say "selfTest PASS"
}

case "$mode" in
  run) run_mode ;;
  selfTest|self-test) self_test_mode ;;
  status) status_mode ;;
  -h|--help|help) usage ;;
  *) usage; die "Unknown mode: $mode" ;;
esac
