#!/usr/bin/env bash
set -euo pipefail

RJ_REPO_BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
ML_PYTHON="${RJ_ML_PYTHON:-/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python}"
STAMP="${RJ_AUAU_SIXPACK_STAMP:-$(date +%Y%m%d_%H%M%S)}"

BASE_SOURCE="${RJ_AUAU_SIXPACK_BASE_SOURCE:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049}"
TRAIN_BASE="${RJ_AUAU_SIXPACK_TRAIN_BASE:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining}"
JET30_SOURCE="${RJ_AUAU_SIXPACK_JET30_SOURCE:-${TRAIN_BASE}/auauTightBDT_${STAMP}_jet30}"
MODEL_BASE="${RJ_AUAU_SIXPACK_MODEL_BASE:-/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models}"
OUTDIR="${OUTDIR:-${MODEL_BASE}/global_etcent_inclusive3_sixpack_${STAMP}}"
STAGE_SOURCE="${STAGE_SOURCE:-${TRAIN_BASE}/globalEtCentInclusive3Sixpack_${STAMP}}"
SUB_ROOT="${SUB_ROOT:-${RJ_REPO_BASE}/condor_sub/auauGlobalSixpack_${STAMP}}"
LOG="${LOG:-${OUTDIR}/global_sixpack_${STAMP}.log}"

VALIDATE_GROUP_SIZE="${VALIDATE_GROUP_SIZE:-100}"
DIRECT_REQUEST_MEMORY="${DIRECT_REQUEST_MEMORY:-32000MB}"
BDT_REQUEST_MEMORY="${BDT_REQUEST_MEMORY:-$DIRECT_REQUEST_MEMORY}"
MLP_REQUEST_MEMORY="${MLP_REQUEST_MEMORY:-64000MB}"
MLP_EVAL_BATCH_SIZE="${MLP_EVAL_BATCH_SIZE:-131072}"
MLP_TRAIN_EVAL_MAX_ROWS="${MLP_TRAIN_EVAL_MAX_ROWS:-500000}"
STACK_REQUEST_MEMORY="${STACK_REQUEST_MEMORY:-64000MB}"
STACK_HIDDEN="${STACK_HIDDEN:-384,192,96,48}"
STACK_EPOCHS="${STACK_EPOCHS:-260}"
STACK_PATIENCE="${STACK_PATIENCE:-45}"
STACK_BATCH_SIZE="${STACK_BATCH_SIZE:-16384}"
STACK_FOLDS="${STACK_FOLDS:-5}"

JET30_RESCUE_GROUP_SIZE="${JET30_RESCUE_GROUP_SIZE:-7}"
JET30_RESCUE_MEMORY="${JET30_RESCUE_MEMORY:-4000MB}"

say() { printf '\033[1;36m[auauSixpack]\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[auauSixpack][WARN]\033[0m %s\n' "$*" >&2; }
die() { printf '\033[1;31m[auauSixpack][ERR]\033[0m %s\n' "$*" >&2; exit 2; }

setup_env() {
  export USER="${USER:-$(id -u -n)}"
  export LOGNAME="${LOGNAME:-$USER}"
  export HOME="/sphenix/u/${LOGNAME}"
  set +u
  source /opt/sphenix/core/bin/sphenix_setup.sh -n
  [[ -d "/sphenix/u/${USER}/thesisAnalysis/install" ]] && source /opt/sphenix/core/bin/setup_local.sh "/sphenix/u/${USER}/thesisAnalysis/install" || true
  [[ -d "/sphenix/u/${USER}/thesisAnalysis_auau/install" ]] && source /opt/sphenix/core/bin/setup_local.sh "/sphenix/u/${USER}/thesisAnalysis_auau/install" || true
  set -u
  local prefix real real_prefix ld=""
  prefix="$(cd "$(dirname "$ML_PYTHON")/.." && pwd -P 2>/dev/null || true)"
  real="$(readlink -f "$ML_PYTHON" 2>/dev/null || true)"
  real_prefix=""
  [[ -n "$real" ]] && real_prefix="$(cd "$(dirname "$real")/.." && pwd -P 2>/dev/null || true)"
  for d in "$prefix/lib" "$prefix/lib64" "$real_prefix/lib" "$real_prefix/lib64"; do
    [[ -d "$d" ]] || continue
    case ":$ld:" in *":$d:"*) ;; *) ld="${ld:+$ld:}$d" ;; esac
  done
  [[ -n "$ld" ]] && export LD_LIBRARY_PATH="$ld:${LD_LIBRARY_PATH:-}"
  unset PYTHONHOME
}

print_plan() {
  say "RECOILJETS_AUAU_GLOBAL_ETCENT_INCLUSIVE3_SIXPACK_V1"
  say "host=$(hostname -f 2>/dev/null || hostname)"
  say "stamp=$STAMP"
  say "repo=$RJ_REPO_BASE"
  say "base_source=$BASE_SOURCE"
  say "jet30_source=$JET30_SOURCE"
  say "stage_source=$STAGE_SOURCE"
  say "outdir=$OUTDIR"
  say "sub_root=$SUB_ROOT"
  say "bdt_memory=$BDT_REQUEST_MEMORY mlp_memory=$MLP_REQUEST_MEMORY mlp_eval_batch=$MLP_EVAL_BATCH_SIZE mlp_train_eval_rows=$MLP_TRAIN_EVAL_MAX_ROWS stack_memory=$STACK_REQUEST_MEMORY stack_hidden=$STACK_HIDDEN stack_batch=$STACK_BATCH_SIZE"
}

jet30_ready() {
  local summary="${JET30_SOURCE}/reports/final_summary.txt"
  [[ -s "$summary" ]] || return 1
  grep -q '^status=READY$' "$summary" || return 1
  grep -q '^expected_background_samples=1$' "$summary" || return 1
}

submit_jet30_rescue() {
  say "Jet30 extraction is not READY; submitting strict Jet30-only rescue extraction."
  (
    cd "$RJ_REPO_BASE"
    export RJ_ML_PYTHON="$ML_PYTHON"
    export RJ_AUAU_TIGHT_BDT_STAMP="${STAMP}_jet30"
    export RJ_AUAU_TIGHT_BDT_SIGNAL_SAMPLES=""
    export RJ_AUAU_TIGHT_BDT_BACKGROUND_SAMPLES="run28_embeddedJet30"
    export RJ_AUAU_TIGHT_BDT_REQUEST_MEMORY="$JET30_RESCUE_MEMORY"
    unset RJ_AUAU_TIGHT_BDT_TOLERATE_ROOT_ABORT_WITH_OUTPUT
    ./scripts/auau_tight_bdt_pipeline.sh condorExtract groupSize "$JET30_RESCUE_GROUP_SIZE"
  )
}

require_jet30_ready() {
  if jet30_ready; then
    say "Jet30 training-tree extraction READY: ${JET30_SOURCE}"
    return 0
  fi
  submit_jet30_rescue
  wait_summary_ready "${JET30_SOURCE}/reports/final_summary.txt" jet30_rescue
  jet30_ready || die "Jet30 rescue summary is READY but provenance checks failed: ${JET30_SOURCE}/reports/final_summary.txt"
  say "Jet30 rescue READY; continuing sixpack training."
}

link_roots() {
  local src="$1" dest="$2"
  mkdir -p "$dest"
  while IFS= read -r path; do
    [[ -n "$path" ]] || continue
    ln -sfn "$path" "${dest}/$(basename "$path")"
  done < <(find "$src" -type f -name '*.root' | sort -V)
}

prepare_source() {
  require_jet30_ready
  mkdir -p "${STAGE_SOURCE}/extraction/signal" "${STAGE_SOURCE}/extraction/background" "${STAGE_SOURCE}/manifests" "${STAGE_SOURCE}/reports"
  link_roots "${BASE_SOURCE}/extraction/signal/run28_embeddedPhoton12" "${STAGE_SOURCE}/extraction/signal/run28_embeddedPhoton12"
  link_roots "${BASE_SOURCE}/extraction/signal/run28_embeddedPhoton20" "${STAGE_SOURCE}/extraction/signal/run28_embeddedPhoton20"
  link_roots "${BASE_SOURCE}/extraction/background/run28_embeddedJet12" "${STAGE_SOURCE}/extraction/background/run28_embeddedJet12"
  link_roots "${BASE_SOURCE}/extraction/background/run28_embeddedJet20" "${STAGE_SOURCE}/extraction/background/run28_embeddedJet20"
  link_roots "${JET30_SOURCE}/extraction/background/run28_embeddedJet30" "${STAGE_SOURCE}/extraction/background/run28_embeddedJet30"
  find "${STAGE_SOURCE}/extraction" -type l -name '*.root' | sort -V > "${STAGE_SOURCE}/manifests/training_roots.list"
  local nroots
  nroots="$(wc -l < "${STAGE_SOURCE}/manifests/training_roots.list" | tr -d ' ')"
  [[ "$nroots" != "0" ]] || die "staged source has no ROOT files"
  cat > "${STAGE_SOURCE}/reports/source_provenance.json" <<EOF
{
  "schema": "AUAU_GLOBAL_ETCENT_INCLUSIVE3_SIXPACK_SOURCE_V1",
  "signal_samples": ["run28_embeddedPhoton12", "run28_embeddedPhoton20"],
  "background_samples": ["run28_embeddedJet12", "run28_embeddedJet20", "run28_embeddedJet30"],
  "base_source": "${BASE_SOURCE}",
  "jet30_source": "${JET30_SOURCE}",
  "staged_source": "${STAGE_SOURCE}",
  "root_count": ${nroots},
  "inclusive3_weights_pb": {
    "run28_embeddedJet12": 1.21692467e6,
    "run28_embeddedJet20": 5.44464934e4,
    "run28_embeddedJet30": 2.40291630e3
  }
}
EOF
  "$ML_PYTHON" "$RJ_REPO_BASE/scripts/train_auau_photon_bdt.py" --task tight --input "@${STAGE_SOURCE}/manifests/training_roots.list" --outdir "${OUTDIR}/preflight_bdt_plan" --campaign global-sixpack --plan-only >/dev/null
  say "Staged source ready: roots=${nroots}"
}

submit_job() {
  local name="$1" mem="$2" script="$3"
  local sub="${SUB_ROOT}/${name}.sub"
  cat > "$sub" <<EOF
universe = vanilla
executable = ${script}
output = ${SUB_ROOT}/${name}.out
error = ${SUB_ROOT}/${name}.err
log = ${SUB_ROOT}/${name}.log
request_memory = ${mem}
notification = Never
+JobBatchName = "auau_global_sixpack_${STAMP}_${name}"
queue
EOF
  condor_submit "$sub" | tee "${SUB_ROOT}/${name}.submit"
}

cluster_from_submit() {
  awk '/submitted to cluster/ {gsub("\\.","",$NF); print $NF; exit}' "$1"
}

wait_cluster_done() {
  local cluster="$1" label="$2"
  while true; do
    local held active
    held="$(condor_q "$cluster" -nobatch 2>/dev/null | awk '/Total for query:/ {for(i=1;i<=NF;i++) if($i=="held,") print $(i-1)}' | tail -1)"
    active="$(condor_q "$cluster" -nobatch 2>/dev/null | awk '/Total for query:/ {print $4; exit}')"
    held="${held:-0}"
    active="${active:-0}"
    if [[ "$held" != "0" ]]; then
      die "${label} cluster ${cluster} has held jobs"
    fi
    [[ "$active" == "0" ]] && break
    say "${label} cluster ${cluster} still active=${active}; sleeping"
    sleep 300
  done
}

wait_summary_ready() {
  local summary="$1" label="$2"
  while true; do
    if [[ -s "$summary" ]]; then
      local status
      status="$(awk -F= '/^status=/ {print $2; exit}' "$summary")"
      [[ "$status" == "READY" ]] && break
      [[ "$status" == "CHECK" || "$status" == "FAILED" ]] && die "${label} summary is ${status}: ${summary}"
    fi
    say "waiting for ${label}: ${summary}"
    sleep 300
  done
}

write_workers() {
  mkdir -p "$OUTDIR" "$SUB_ROOT"
  cat > "${SUB_ROOT}/train_bdt.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
cd "$RJ_REPO_BASE"
export RJ_ML_PYTHON="$ML_PYTHON"
set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n
set -u
"$ML_PYTHON" scripts/train_auau_photon_bdt.py --task tight --input "@${STAGE_SOURCE}/manifests/training_roots.list" --outdir "${OUTDIR}/bdt_models" --campaign global-sixpack --test-size 0.10 --random-seed 13 --parallel-workers 2
EOF
  cat > "${SUB_ROOT}/train_mlp_noiso.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
cd "$RJ_REPO_BASE"
export RJ_ML_PYTHON="$ML_PYTHON"
set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n
set -u
"$ML_PYTHON" scripts/train_auau_photon_mlp.py --source "$STAGE_SOURCE" --outdir "${OUTDIR}/mlp_models" --products globalEtCent1535_mlp_noIso --registry-output "${OUTDIR}/mlp_models/model_registry_noIso.json" --pt-range 15:35 --validation-fraction 0.10 --test-fraction 0.10 --random-seed 137 --max-rows 0 --max-rows-per-class 0 --max-rows-per-pt-bin-class 0 --hidden-layers 256,128,64 --hidden-layer-grid "256,128,64;384,192,96,48" --epochs 280 --patience 55 --batch-size 8192 --eval-batch-size "$MLP_EVAL_BATCH_SIZE" --train-eval-max-rows "$MLP_TRAIN_EVAL_MAX_ROWS" --learning-rate 5.0e-4 --l2 7.0e-5 --selection-metric validation_auc
cp -f "${OUTDIR}/mlp_models/training_manifest_summary.json" "${OUTDIR}/mlp_models/training_manifest_summary_noIso.json"
cp -f "${OUTDIR}/mlp_models/training_read_summary.json" "${OUTDIR}/mlp_models/training_read_summary_noIso.json"
EOF
  cat > "${SUB_ROOT}/train_mlp_iso.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
cd "$RJ_REPO_BASE"
export RJ_ML_PYTHON="$ML_PYTHON"
set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n
set -u
"$ML_PYTHON" scripts/train_auau_photon_mlp.py --source "$STAGE_SOURCE" --outdir "${OUTDIR}/mlp_models" --products globalEtCent1535_mlp_iso --registry-output "${OUTDIR}/mlp_models/model_registry_iso.json" --pt-range 15:35 --validation-fraction 0.10 --test-fraction 0.10 --random-seed 137 --max-rows 0 --max-rows-per-class 0 --max-rows-per-pt-bin-class 0 --hidden-layers 256,128,64 --hidden-layer-grid "256,128,64;384,192,96,48" --epochs 280 --patience 55 --batch-size 8192 --eval-batch-size "$MLP_EVAL_BATCH_SIZE" --train-eval-max-rows "$MLP_TRAIN_EVAL_MAX_ROWS" --learning-rate 5.0e-4 --l2 7.0e-5 --selection-metric validation_auc
cp -f "${OUTDIR}/mlp_models/training_manifest_summary.json" "${OUTDIR}/mlp_models/training_manifest_summary_iso.json"
cp -f "${OUTDIR}/mlp_models/training_read_summary.json" "${OUTDIR}/mlp_models/training_read_summary_iso.json"
EOF
  cat > "${SUB_ROOT}/stack_noiso.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
cd "$RJ_REPO_BASE"
export RJ_ML_PYTHON="$ML_PYTHON"
set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n
set -u
"$ML_PYTHON" scripts/train_auau_oof_residual_superstacker.py --mlp-cache "${OUTDIR}/validation/mlp/score_caches.list" --bdt-cache "${OUTDIR}/validation/bdt/score_caches.list" --outdir "${OUTDIR}/stack_noIso" --mlp-score score_globalEtCent1535_mlp_noIso --bdt-score score_globalEtCent1535_bdt_noIso --require-full-stat --expected-shards 1 --train-fraction 0.80 --val-fraction 0.10 --folds "$STACK_FOLDS" --lower-algorithms logistic,gbm,nn --lower-feature-mode full_features --super-feature-mode scores_plus_full_features --final-mode direct --model-name globalEtCent1535_oofSuperNN_noIso --hidden "$STACK_HIDDEN" --epochs "$STACK_EPOCHS" --patience "$STACK_PATIENCE" --batch-size "$STACK_BATCH_SIZE"
EOF
  cat > "${SUB_ROOT}/stack_iso.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
cd "$RJ_REPO_BASE"
export RJ_ML_PYTHON="$ML_PYTHON"
set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n
set -u
"$ML_PYTHON" scripts/train_auau_oof_residual_superstacker.py --mlp-cache "${OUTDIR}/validation/mlp/score_caches.list" --bdt-cache "${OUTDIR}/validation/bdt/score_caches.list" --outdir "${OUTDIR}/stack_iso" --mlp-score score_globalEtCent1535_mlp_iso --bdt-score score_globalEtCent1535_bdt_iso --require-full-stat --expected-shards 1 --train-fraction 0.80 --val-fraction 0.10 --folds "$STACK_FOLDS" --lower-algorithms logistic,gbm,nn --lower-feature-mode full_features --super-feature-mode scores_plus_full_features --include-isolation-context --final-mode direct --model-name globalEtCent1535_oofSuperNN_iso --hidden "$STACK_HIDDEN" --epochs "$STACK_EPOCHS" --patience "$STACK_PATIENCE" --batch-size "$STACK_BATCH_SIZE"
EOF
  chmod +x "${SUB_ROOT}"/*.sh
}

merge_mlp_registries() {
  "$ML_PYTHON" - "${OUTDIR}/mlp_models/model_registry_noIso.json" "${OUTDIR}/mlp_models/model_registry_iso.json" "${OUTDIR}/mlp_models/model_registry.json" <<'PY'
import json
import sys
from pathlib import Path

noiso_path, iso_path, out_path = map(Path, sys.argv[1:4])
registries = [json.loads(noiso_path.read_text()), json.loads(iso_path.read_text())]
merged = dict(registries[0])
merged["schema"] = "RJ_AUAU_TIGHT_MLP_REGISTRY_V1"
merged["status"] = "READY"
merged["products"] = []
merged["models"] = []
merged["merged_from"] = [str(noiso_path), str(iso_path)]
seen = set()
for registry in registries:
    for product in registry.get("products", []):
        if product not in seen:
            merged["products"].append(product)
            seen.add(product)
    merged["models"].extend(registry.get("models", []))
optional = set()
for registry in registries:
    optional.update(registry.get("optional_branches_seen", []))
merged["optional_branches_seen"] = sorted(optional)
out_path.write_text(json.dumps(merged, indent=2, sort_keys=True) + "\n")
print(f"[auauSixpack] merged MLP registry: {out_path}")
PY
}

validate_bdt() {
  RJ_AUAU_TIGHT_BDT_VALIDATE_TOTAL_SCORE_MAX_ROWS=0 RJ_AUAU_TIGHT_BDT_VALIDATE_REQUEST_MEMORY=4000MB \
    ./scripts/auau_tight_bdt_pipeline.sh validateOnSimCondor SOURCE="$STAGE_SOURCE" MODEL_DIR="${OUTDIR}/bdt_models" MODEL_REGISTRY="${OUTDIR}/bdt_models/model_registry.json" OUTDIR="${OUTDIR}/validation/bdt" groupSize "$VALIDATE_GROUP_SIZE"
}

validate_mlp() {
  RJ_AUAU_TIGHT_MLP_VALIDATE_TOTAL_SCORE_MAX_ROWS=0 RJ_AUAU_TIGHT_MLP_VALIDATE_REQUEST_MEMORY=4000MB \
    ./scripts/auau_tight_mlp_pipeline.sh validateOnSimCondor SOURCE="$STAGE_SOURCE" MODEL_DIR="${OUTDIR}/mlp_models" MODEL_REGISTRY="${OUTDIR}/mlp_models/model_registry.json" OUTDIR="${OUTDIR}/validation/mlp" groupSize "$VALIDATE_GROUP_SIZE"
}

write_final_summary() {
  cat > "${OUTDIR}/sixpack_summary.txt" <<EOF
RECOILJETS_AUAU_GLOBAL_ETCENT_INCLUSIVE3_SIXPACK_SUMMARY_V1
status=READY
outdir=${OUTDIR}
stage_source=${STAGE_SOURCE}
bdt_validation=${OUTDIR}/validation/bdt
mlp_validation=${OUTDIR}/validation/mlp
stack_noiso=${OUTDIR}/stack_noIso
stack_iso=${OUTDIR}/stack_iso
EOF
  cat "${OUTDIR}/sixpack_summary.txt"
}

orchestrate() {
  setup_env
  print_plan | tee -a "$LOG"
  mkdir -p "$OUTDIR" "$SUB_ROOT"
  prepare_source | tee -a "$LOG"
  write_workers
  submit_job train_bdt "$BDT_REQUEST_MEMORY" "${SUB_ROOT}/train_bdt.sh"
  wait_cluster_done "$(cluster_from_submit "${SUB_ROOT}/train_bdt.submit")" train_bdt
  run_remaining_after_bdt
}

run_remaining_after_bdt() {
  setup_env
  print_plan | tee -a "$LOG"
  mkdir -p "$OUTDIR" "$SUB_ROOT"
  [[ -s "${OUTDIR}/bdt_models/model_registry.json" ]] || die "missing BDT registry: ${OUTDIR}/bdt_models/model_registry.json"
  [[ -s "${STAGE_SOURCE}/manifests/training_roots.list" ]] || die "missing staged training roots: ${STAGE_SOURCE}/manifests/training_roots.list"
  write_workers
  submit_job train_mlp_noiso "$MLP_REQUEST_MEMORY" "${SUB_ROOT}/train_mlp_noiso.sh"
  wait_cluster_done "$(cluster_from_submit "${SUB_ROOT}/train_mlp_noiso.submit")" train_mlp_noiso
  submit_job train_mlp_iso "$MLP_REQUEST_MEMORY" "${SUB_ROOT}/train_mlp_iso.sh"
  wait_cluster_done "$(cluster_from_submit "${SUB_ROOT}/train_mlp_iso.submit")" train_mlp_iso
  merge_mlp_registries
  validate_bdt | tee -a "$LOG"
  validate_mlp | tee -a "$LOG"
  wait_summary_ready "${OUTDIR}/validation/bdt/validation_summary.txt" bdt_validation
  wait_summary_ready "${OUTDIR}/validation/mlp/validation_summary.txt" mlp_validation
  cp "${OUTDIR}/validation/bdt/score_caches.list" "${OUTDIR}/validation/bdt/score_caches.fullstat.list"
  cp "${OUTDIR}/validation/mlp/score_caches.list" "${OUTDIR}/validation/mlp/score_caches.fullstat.list"
  submit_job stack_noiso "$STACK_REQUEST_MEMORY" "${SUB_ROOT}/stack_noiso.sh"
  wait_cluster_done "$(cluster_from_submit "${SUB_ROOT}/stack_noiso.submit")" stack_noiso
  submit_job stack_iso "$STACK_REQUEST_MEMORY" "${SUB_ROOT}/stack_iso.sh"
  wait_cluster_done "$(cluster_from_submit "${SUB_ROOT}/stack_iso.submit")" stack_iso
  [[ -s "${OUTDIR}/stack_noIso/oof_direct_supernn_metrics.json" ]] || die "missing no-iso stack metrics"
  [[ -s "${OUTDIR}/stack_iso/oof_direct_supernn_metrics.json" ]] || die "missing iso stack metrics"
  write_final_summary
}

status() {
  print_plan
  for path in \
    "${JET30_SOURCE}/reports/final_summary.txt" \
    "${OUTDIR}/validation/bdt/validation_summary.txt" \
    "${OUTDIR}/validation/mlp/validation_summary.txt" \
    "${OUTDIR}/stack_noIso/oof_direct_supernn_metrics.json" \
    "${OUTDIR}/stack_iso/oof_direct_supernn_metrics.json" \
    "${OUTDIR}/sixpack_summary.txt"; do
    if [[ -s "$path" ]]; then
      say "exists: $path"
      head -20 "$path" || true
    else
      warn "missing: $path"
    fi
  done
}

mode="${1:-train}"
case "$mode" in
  train)
    if [[ "${RJ_AUAU_SIXPACK_IN_TMUX:-0}" == "1" ]]; then
      orchestrate
    else
      mkdir -p "$OUTDIR" "$SUB_ROOT"
      tmux new-session -d -s "auau_sixpack_${STAMP}" "cd '$RJ_REPO_BASE' && RJ_AUAU_SIXPACK_IN_TMUX=1 RJ_AUAU_SIXPACK_STAMP='$STAMP' OUTDIR='$OUTDIR' STAGE_SOURCE='$STAGE_SOURCE' SUB_ROOT='$SUB_ROOT' LOG='$LOG' BDT_REQUEST_MEMORY='$BDT_REQUEST_MEMORY' MLP_REQUEST_MEMORY='$MLP_REQUEST_MEMORY' MLP_EVAL_BATCH_SIZE='$MLP_EVAL_BATCH_SIZE' MLP_TRAIN_EVAL_MAX_ROWS='$MLP_TRAIN_EVAL_MAX_ROWS' STACK_REQUEST_MEMORY='$STACK_REQUEST_MEMORY' STACK_BATCH_SIZE='$STACK_BATCH_SIZE' '$RJ_REPO_BASE/scripts/submit_auau_global_sixpack_oof_stack.sh' train"
      say "tmux=auau_sixpack_${STAMP}"
      say "log=$LOG"
    fi
    ;;
  orchestrate) orchestrate ;;
  resumeAfterBdt)
    if [[ "${RJ_AUAU_SIXPACK_IN_TMUX:-0}" == "1" ]]; then
      run_remaining_after_bdt
    else
      mkdir -p "$OUTDIR" "$SUB_ROOT"
      tmux new-session -d -s "auau_sixpack_resume_${STAMP}" "cd '$RJ_REPO_BASE' && RJ_AUAU_SIXPACK_IN_TMUX=1 RJ_AUAU_SIXPACK_STAMP='$STAMP' OUTDIR='$OUTDIR' STAGE_SOURCE='$STAGE_SOURCE' SUB_ROOT='$SUB_ROOT' LOG='$LOG' BDT_REQUEST_MEMORY='$BDT_REQUEST_MEMORY' MLP_REQUEST_MEMORY='$MLP_REQUEST_MEMORY' MLP_EVAL_BATCH_SIZE='$MLP_EVAL_BATCH_SIZE' MLP_TRAIN_EVAL_MAX_ROWS='$MLP_TRAIN_EVAL_MAX_ROWS' STACK_REQUEST_MEMORY='$STACK_REQUEST_MEMORY' STACK_BATCH_SIZE='$STACK_BATCH_SIZE' '$RJ_REPO_BASE/scripts/submit_auau_global_sixpack_oof_stack.sh' resumeAfterBdt"
      say "tmux=auau_sixpack_resume_${STAMP}"
      say "log=$LOG"
    fi
    ;;
  preflight) setup_env; print_plan; require_jet30_ready; prepare_source ;;
  rescueJet30) setup_env; print_plan; submit_jet30_rescue ;;
  status) status ;;
  *) die "Unknown mode: $mode" ;;
esac
