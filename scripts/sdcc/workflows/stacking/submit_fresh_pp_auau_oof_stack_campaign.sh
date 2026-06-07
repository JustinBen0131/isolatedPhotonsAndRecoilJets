#!/usr/bin/env bash
set -euo pipefail

MODE="${1:-plan}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
RJ_REPO_BASE="${RJ_REPO_BASE:-$(cd "${SCRIPT_DIR}/../../../.." && pwd -P)}"
ML_PYTHON="${RJ_ML_PYTHON:-/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python}"
STAMP="${RJ_FRESH_OOF_STACK_STAMP:-$(date +%Y%m%d_%H%M%S)}"
STACK_MODE="${RJ_FRESH_OOF_STACK_MODE:-oof5}"
SIMPLE_STACK_FOLD="${RJ_FRESH_OOF_SIMPLE_STACK_FOLD:-0}"

MODEL_BASE="${RJ_FRESH_OOF_MODEL_BASE:-/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models}"
RUN_KIND_PREFIX="${RJ_FRESH_OOF_RUN_KIND_PREFIX:-fresh_pp_auau_bdt_mlp_oof_stack}"
SUB_KIND_PREFIX="${RJ_FRESH_OOF_SUB_KIND_PREFIX:-freshOOFStack}"
if [[ "$STACK_MODE" == "simple_holdout" && -z "${RJ_FRESH_OOF_RUN_KIND_PREFIX:-}" ]]; then
  RUN_KIND_PREFIX="fresh_pp_auau_bdt_mlp_simple_holdout_stack"
fi
if [[ "$STACK_MODE" == "simple_holdout" && -z "${RJ_FRESH_OOF_SUB_KIND_PREFIX:-}" ]]; then
  SUB_KIND_PREFIX="freshSimpleHoldoutStack"
fi
RUN_ROOT="${RJ_FRESH_OOF_RUN_ROOT:-${MODEL_BASE}/${RUN_KIND_PREFIX}_${STAMP}}"
SUB_ROOT="${RJ_FRESH_OOF_SUB_ROOT:-${RJ_REPO_BASE}/condor_sub/${SUB_KIND_PREFIX}_${STAMP}}"
LOG_ROOT="${RJ_FRESH_OOF_LOG_ROOT:-${RUN_ROOT}/logs}"

AUAU_SOURCE="${RJ_FRESH_OOF_AUAU_SOURCE:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/THE8_branchA_ladder_jet12_20_30_40_20260527}"
AUAU_MANIFEST="${RJ_FRESH_OOF_AUAU_MANIFEST:-${AUAU_SOURCE}/manifests/training_roots.list}"

PP_RUN_ROOT="${RJ_FRESH_OOF_PP_RUN_ROOT:-/sphenix/u/patsfan753/scratch/thesisAnalysis/dataOutput/ppPhotonMLPipeline/ppg12_matched_20260521_230310}"
PP_TRAIN_MANIFEST="${RJ_FRESH_OOF_PP_TRAIN_MANIFEST:-${PP_RUN_ROOT}/training_roots_currentIAN.list}"
PP_SIGNAL_MANIFEST="${RJ_FRESH_OOF_PP_SIGNAL_MANIFEST:-${PP_RUN_ROOT}/signal_roots_currentIAN.list}"
PP_INCLUSIVE_MANIFEST="${RJ_FRESH_OOF_PP_INCLUSIVE_MANIFEST:-${PP_RUN_ROOT}/inclusive_roots_currentIAN.list}"

PP_TREE="${RJ_FRESH_OOF_PP_TREE:-AuAuPhotonIDTrainingTree}"
AUAU_TREE="${RJ_FRESH_OOF_AUAU_TREE:-AuAuPhotonIDTrainingTree}"
REQUEST_CPUS="${RJ_FRESH_OOF_REQUEST_CPUS:-2}"
REQUEST_MEMORY="${RJ_FRESH_OOF_REQUEST_MEMORY:-128000MB}"
MAX_LOAD_ROWS_PER_CLASS="${RJ_FRESH_OOF_MAX_LOAD_ROWS_PER_CLASS:-0}"
MAX_LOAD_ROWS="${RJ_FRESH_OOF_MAX_LOAD_ROWS:-0}"
BDT_ESTIMATORS="${RJ_FRESH_OOF_BDT_ESTIMATORS:-750}"
PREDICT_CHUNK_ROWS="${RJ_FRESH_OOF_PREDICT_CHUNK_ROWS:-750000}"
STACK_MLP_EPOCHS="${RJ_FRESH_OOF_STACK_MLP_EPOCHS:-180}"
STACK_MLP_HIDDEN="${RJ_FRESH_OOF_STACK_MLP_HIDDEN:-256,128,64}"
STACK_MLP_BATCH_SIZE="${RJ_FRESH_OOF_STACK_MLP_BATCH_SIZE:-8192}"

say() { printf '\033[1;36m[freshOOFStack]\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[freshOOFStack][WARN]\033[0m %s\n' "$*" >&2; }
die() { printf '\033[1;31m[freshOOFStack][ERR]\033[0m %s\n' "$*" >&2; exit 2; }

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
  export OMP_NUM_THREADS="$REQUEST_CPUS"
  export OPENBLAS_NUM_THREADS="$REQUEST_CPUS"
  export MKL_NUM_THREADS="$REQUEST_CPUS"
  export NUMEXPR_NUM_THREADS="$REQUEST_CPUS"
  export XGBOOST_NUM_THREADS="$REQUEST_CPUS"
  export MALLOC_ARENA_MAX="${RJ_FRESH_OOF_MALLOC_ARENA_MAX:-2}"
  unset PYTHONHOME
}

print_plan() {
  say "schema=RJ_FRESH_PP_AUAU_BDT_MLP_STACK_CAMPAIGN_V2"
  say "mode=$MODE"
  say "stack_training_mode=$STACK_MODE"
  say "simple_stack_fold=$SIMPLE_STACK_FOLD"
  say "stamp=$STAMP"
  say "repo=$RJ_REPO_BASE"
  say "run_root=$RUN_ROOT"
  say "sub_root=$SUB_ROOT"
  say "ml_python=$ML_PYTHON"
  say "auau_manifest=$AUAU_MANIFEST"
  say "pp_train_manifest=$PP_TRAIN_MANIFEST"
  say "pp_overlay_manifests=$PP_SIGNAL_MANIFEST $PP_INCLUSIVE_MANIFEST"
  say "trees pp=$PP_TREE auau=$AUAU_TREE"
  say "request_cpus=$REQUEST_CPUS request_memory=$REQUEST_MEMORY"
  say "max_load_rows_per_class=$MAX_LOAD_ROWS_PER_CLASS max_load_rows=$MAX_LOAD_ROWS"
  say "predict_chunk_rows=$PREDICT_CHUNK_ROWS"
}

need_file() {
  [[ -s "$1" ]] || die "missing required file: $1"
}

require_driver() {
  need_file "${RJ_REPO_BASE}/scripts/ml/stacking/train_photon_bdt_mlp_oof_stack.py"
}

require_stack_mode() {
  case "$STACK_MODE" in
    oof5|simple_holdout) ;;
    *) die "RJ_FRESH_OOF_STACK_MODE must be oof5 or simple_holdout, got: $STACK_MODE" ;;
  esac
  [[ "$SIMPLE_STACK_FOLD" =~ ^[0-9]+$ ]] || die "RJ_FRESH_OOF_SIMPLE_STACK_FOLD must be an integer"
}

require_manifests() {
  need_file "$AUAU_MANIFEST"
  need_file "$PP_TRAIN_MANIFEST"
  need_file "$PP_SIGNAL_MANIFEST"
  need_file "$PP_INCLUSIVE_MANIFEST"
}

require_submit_provenance() {
  [[ "${RJ_FRESH_OOF_APPROVED:-0}" == "1" ]] || die "submit requires RJ_FRESH_OOF_APPROVED=1"
  [[ -n "${RJ_CODEX_CHAT_NAME:-}" ]] || die "submit requires RJ_CODEX_CHAT_NAME"
  [[ -n "${RJ_CODEX_THREAD_ID:-}" ]] || die "submit requires RJ_CODEX_THREAD_ID"
}

queue_snapshot() {
  say "condor queue snapshot:"
  condor_q "$USER" -nobatch 2>/dev/null | tail -20 || true
}

preflight() {
  setup_env
  print_plan
  require_driver
  require_stack_mode
  require_manifests
  queue_snapshot
  mkdir -p "$RUN_ROOT" "$LOG_ROOT"
  "$ML_PYTHON" "${RJ_REPO_BASE}/scripts/ml/stacking/train_photon_bdt_mlp_oof_stack.py" \
    --domain pp \
    --self-test \
    --self-test-rows 1200 \
    --outdir "${RUN_ROOT}/selftest_pp" \
    --feature-preset pp_basev3e_noiso \
    --pt-range 5:35 \
    --centrality-range=-1:0 \
    --weight-mode none \
    --folds 3 \
    --stack-training-mode "$STACK_MODE" \
    --simple-stack-fold "$SIMPLE_STACK_FOLD" \
    --locked-test-fraction 0.20 \
    --bdt-estimators 20 \
    --stack-gbm-estimators 10 \
    --stack-mlp-epochs 3 \
    --stack-mlp-patience 2 \
    --stack-mlp-batch-size 256
  "$ML_PYTHON" "${RJ_REPO_BASE}/scripts/ml/stacking/train_photon_bdt_mlp_oof_stack.py" \
    --domain auau \
    --self-test \
    --self-test-rows 1200 \
    --outdir "${RUN_ROOT}/selftest_auau" \
    --feature-preset auau_global_noiso \
    --pt-range 15:35 \
    --centrality-range 0:80 \
    --weight-mode none \
    --folds 3 \
    --stack-training-mode "$STACK_MODE" \
    --simple-stack-fold "$SIMPLE_STACK_FOLD" \
    --locked-test-fraction 0.20 \
    --bdt-estimators 20 \
    --stack-gbm-estimators 10 \
    --stack-mlp-epochs 3 \
    --stack-mlp-patience 2 \
    --stack-mlp-batch-size 256
  write_campaign_metadata "PREFLIGHT_PASS"
}

write_campaign_metadata() {
  local status="$1"
  mkdir -p "$RUN_ROOT"
  cat > "${RUN_ROOT}/campaign_submit_manifest.json" <<EOF
{
  "schema": "RJ_FRESH_PP_AUAU_BDT_MLP_STACK_SUBMIT_V2",
  "status": "${status}",
  "stack_training_mode": "${STACK_MODE}",
  "simple_stack_fold": ${SIMPLE_STACK_FOLD},
  "stamp": "${STAMP}",
  "run_root": "${RUN_ROOT}",
  "sub_root": "${SUB_ROOT}",
  "repo": "${RJ_REPO_BASE}",
  "ml_python": "${ML_PYTHON}",
  "request_cpus": "${REQUEST_CPUS}",
  "request_memory": "${REQUEST_MEMORY}",
  "predict_chunk_rows": ${PREDICT_CHUNK_ROWS},
  "max_load_rows_per_class": "${MAX_LOAD_ROWS_PER_CLASS}",
  "max_load_rows": "${MAX_LOAD_ROWS}",
  "codex_chat_name": "${RJ_CODEX_CHAT_NAME:-}",
  "codex_thread_id": "${RJ_CODEX_THREAD_ID:-}",
  "domains": {
    "pp": {
      "train_manifest": "${PP_TRAIN_MANIFEST}",
      "overlay_manifests": ["${PP_SIGNAL_MANIFEST}", "${PP_INCLUSIVE_MANIFEST}"],
      "tree": "${PP_TREE}",
      "outdir": "${RUN_ROOT}/pp",
      "feature_preset": "pp_basev3e_noiso",
      "training_label_definition": "is_signal == 1 vs is_signal == 0",
      "training_sources": {
        "Photon+Jet": ["run28_photonjet5", "run28_photonjet10", "run28_photonjet20"],
        "inclusive_training": ["run28_jet8", "run28_jet12", "run28_jet20", "run28_jet30"],
        "inclusive_overlay": ["run28_jet8", "run28_jet12", "run28_jet20", "run28_jet30", "run28_jet40"]
      }
    },
    "auau": {
      "train_manifest": "${AUAU_MANIFEST}",
      "overlay_manifests": ["${AUAU_MANIFEST}"],
      "tree": "${AUAU_TREE}",
      "outdir": "${RUN_ROOT}/auau",
      "feature_preset": "auau_global_noiso",
      "training_label_definition": "is_signal == 1 vs is_signal == 0",
      "training_sources": {
        "Photon+Jet": ["run28_embeddedPhoton12", "run28_embeddedPhoton20"],
        "inclusive_training": ["run28_embeddedJet12", "run28_embeddedJet20", "run28_embeddedJet30", "run28_embeddedJet40"]
      }
    }
  }
}
EOF
}

write_pp_worker() {
  mkdir -p "$SUB_ROOT" "$RUN_ROOT" "$LOG_ROOT"
  cat > "${SUB_ROOT}/run_pp.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
cd "$RJ_REPO_BASE"
ML_PYTHON="$ML_PYTHON"
REQUEST_CPUS="$REQUEST_CPUS"
PREDICT_CHUNK_ROWS="$PREDICT_CHUNK_ROWS"
export RJ_CODEX_CHAT_NAME="${RJ_CODEX_CHAT_NAME:-}"
export RJ_CODEX_THREAD_ID="${RJ_CODEX_THREAD_ID:-}"
$(declare -f setup_env)
setup_env
"$ML_PYTHON" scripts/ml/stacking/train_photon_bdt_mlp_oof_stack.py \\
  --domain pp \\
  --input "@${PP_TRAIN_MANIFEST}" \\
  --overlay-input "@${PP_SIGNAL_MANIFEST}" "@${PP_INCLUSIVE_MANIFEST}" \\
  --tree "$PP_TREE" \\
  --outdir "${RUN_ROOT}/pp" \\
  --feature-preset pp_basev3e_noiso \\
  --pt-range 5:35 \\
  --centrality-range=-1:0 \\
  --weight-mode ppg12-exact \\
  --ppg12-exact-expected-samples run28_photonjet5,run28_photonjet10,run28_photonjet20,run28_jet8,run28_jet12,run28_jet20,run28_jet30 \\
  --training-inclusive-sources run28_jet8,run28_jet12,run28_jet20,run28_jet30 \\
  --overlay-signal-sources run28_photonjet5,run28_photonjet10,run28_photonjet20 \\
  --overlay-inclusive-sources run28_jet8,run28_jet12,run28_jet20,run28_jet30,run28_jet40 \\
  --locked-test-fraction 0.20 \\
  --folds 5 \\
  --stack-training-mode "$STACK_MODE" \\
  --simple-stack-fold "$SIMPLE_STACK_FOLD" \\
  --random-seed 260604 \\
  --report-et-bins 5,10,14,18,22,35 \\
  --max-load-rows-per-class "$MAX_LOAD_ROWS_PER_CLASS" \\
  --max-load-rows "$MAX_LOAD_ROWS" \\
  --skip-missing-tree \\
  --n-jobs "$REQUEST_CPUS" \\
  --predict-chunk-rows "$PREDICT_CHUNK_ROWS" \\
  --bdt-estimators "$BDT_ESTIMATORS" \\
  --stack-mlp-hidden "$STACK_MLP_HIDDEN" \\
  --stack-mlp-epochs "$STACK_MLP_EPOCHS" \\
  --stack-mlp-batch-size "$STACK_MLP_BATCH_SIZE"
EOF
  chmod +x "${SUB_ROOT}/run_pp.sh"
}

write_auau_worker() {
  mkdir -p "$SUB_ROOT" "$RUN_ROOT" "$LOG_ROOT"
  cat > "${SUB_ROOT}/run_auau.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
cd "$RJ_REPO_BASE"
ML_PYTHON="$ML_PYTHON"
REQUEST_CPUS="$REQUEST_CPUS"
PREDICT_CHUNK_ROWS="$PREDICT_CHUNK_ROWS"
export RJ_CODEX_CHAT_NAME="${RJ_CODEX_CHAT_NAME:-}"
export RJ_CODEX_THREAD_ID="${RJ_CODEX_THREAD_ID:-}"
$(declare -f setup_env)
setup_env
"$ML_PYTHON" scripts/ml/stacking/train_photon_bdt_mlp_oof_stack.py \\
  --domain auau \\
  --input "@${AUAU_MANIFEST}" \\
  --overlay-input "@${AUAU_MANIFEST}" \\
  --tree "$AUAU_TREE" \\
  --outdir "${RUN_ROOT}/auau" \\
  --feature-preset auau_global_noiso \\
  --pt-range 15:35 \\
  --centrality-range 0:80 \\
  --weight-mode ppg12-exact \\
  --ppg12-exact-expected-samples run28_embeddedPhoton12,run28_embeddedPhoton20,run28_embeddedJet12,run28_embeddedJet20,run28_embeddedJet30,run28_embeddedJet40 \\
  --training-inclusive-sources run28_embeddedJet12,run28_embeddedJet20,run28_embeddedJet30,run28_embeddedJet40 \\
  --overlay-signal-sources run28_embeddedPhoton12,run28_embeddedPhoton20 \\
  --overlay-inclusive-sources run28_embeddedJet12,run28_embeddedJet20,run28_embeddedJet30,run28_embeddedJet40 \\
  --locked-test-fraction 0.20 \\
  --folds 5 \\
  --stack-training-mode "$STACK_MODE" \\
  --simple-stack-fold "$SIMPLE_STACK_FOLD" \\
  --random-seed 260604 \\
  --report-et-bins 15,17,19,21,23,25,27,30,35 \\
  --report-cent-bins 0,20,40,60,80 \\
  --max-load-rows-per-class "$MAX_LOAD_ROWS_PER_CLASS" \\
  --max-load-rows "$MAX_LOAD_ROWS" \\
  --skip-missing-tree \\
  --n-jobs "$REQUEST_CPUS" \\
  --predict-chunk-rows "$PREDICT_CHUNK_ROWS" \\
  --bdt-estimators "$BDT_ESTIMATORS" \\
  --stack-mlp-hidden "$STACK_MLP_HIDDEN" \\
  --stack-mlp-epochs "$STACK_MLP_EPOCHS" \\
  --stack-mlp-batch-size "$STACK_MLP_BATCH_SIZE"
EOF
  chmod +x "${SUB_ROOT}/run_auau.sh"
}

write_workers() {
  write_pp_worker
  write_auau_worker
}

submit_one() {
  local name="$1"
  local executable="$2"
  local sub="${SUB_ROOT}/${name}.sub"
  cat > "$sub" <<EOF
universe = vanilla
executable = ${executable}
output = ${SUB_ROOT}/${name}.out
error = ${SUB_ROOT}/${name}.err
log = ${SUB_ROOT}/${name}.log
request_cpus = ${REQUEST_CPUS}
request_memory = ${REQUEST_MEMORY}
notification = Never
getenv = True
+JobBatchName = "fresh_${STACK_MODE}_stack_${STAMP}_${name}"
queue
EOF
  condor_submit "$sub" | tee "${SUB_ROOT}/${name}.submit"
}

submit_campaign() {
  setup_env
  print_plan
  require_submit_provenance
  require_driver
  require_stack_mode
  require_manifests
  mkdir -p "$RUN_ROOT" "$SUB_ROOT" "$LOG_ROOT"
  write_campaign_metadata "SUBMITTING"
  write_workers
  submit_one pp "${SUB_ROOT}/run_pp.sh"
  submit_one auau "${SUB_ROOT}/run_auau.sh"
  write_campaign_metadata "SUBMITTED"
  status
}

next_retry_name() {
  local domain="$1" n=1
  while [[ -e "${SUB_ROOT}/${domain}_retry${n}.sub" || -e "${SUB_ROOT}/${domain}_retry${n}.submit" ]]; do
    n=$((n + 1))
  done
  printf '%s_retry%d' "$domain" "$n"
}

submit_single_domain() {
  local domain="$1" worker=""
  setup_env
  print_plan
  require_submit_provenance
  require_driver
  require_stack_mode
  require_manifests
  mkdir -p "$RUN_ROOT" "$SUB_ROOT" "$LOG_ROOT"
  case "$domain" in
    pp)
      write_pp_worker
      worker="${SUB_ROOT}/run_pp.sh"
      ;;
    auau)
      write_auau_worker
      worker="${SUB_ROOT}/run_auau.sh"
      ;;
    *) die "submit_single_domain requires pp or auau, got: $domain" ;;
  esac
  local submit_name="$domain"
  if [[ -e "${SUB_ROOT}/${domain}.submit" || -e "${SUB_ROOT}/${domain}.sub" ]]; then
    submit_name="$(next_retry_name "$domain")"
    warn "existing ${domain} submit file found; submitting scoped retry as ${submit_name}"
  fi
  write_campaign_metadata "SUBMITTING_${submit_name}"
  submit_one "$submit_name" "$worker"
  write_campaign_metadata "SUBMITTED_${submit_name}"
  status
}

cluster_from_submit() {
  awk '/submitted to cluster/ {gsub("\\.","",$NF); print $NF; exit}' "$1" 2>/dev/null || true
}

latest_submit_file() {
  local domain="$1"
  local files=()
  shopt -s nullglob
  files=("${SUB_ROOT}/${domain}.submit" "${SUB_ROOT}/${domain}"_retry*.submit)
  shopt -u nullglob
  (( ${#files[@]} > 0 )) || return 0
  ls -t "${files[@]}" 2>/dev/null | head -1
}

status() {
  print_plan
  for domain in pp auau; do
    local submit_file=""
    local cluster=""
    submit_file="$(latest_submit_file "$domain")"
    [[ -s "$submit_file" ]] && cluster="$(cluster_from_submit "$submit_file")"
    if [[ -n "$cluster" ]]; then
      say "${domain}_cluster=${cluster} submit_file=${submit_file}"
      condor_q "$cluster" -nobatch 2>/dev/null || true
    else
      warn "no submit cluster recorded yet for ${domain}"
    fi
    for path in \
      "${RUN_ROOT}/${domain}/campaign_manifest.json" \
      "${RUN_ROOT}/${domain}/model_metrics.csv" \
      "${RUN_ROOT}/${domain}/stratified_metrics.csv" \
      "${RUN_ROOT}/${domain}/overlay_histograms.csv" \
      "${RUN_ROOT}/${domain}/score_correlations.json"; do
      if [[ -s "$path" ]]; then
        say "exists: $path"
      else
        warn "missing: $path"
      fi
    done
  done
}

wait_ready() {
  setup_env
  print_plan
  local interval="${RJ_FRESH_OOF_WAIT_INTERVAL:-900}"
  while true; do
    local pp_ready=0 auau_ready=0
    [[ -s "${RUN_ROOT}/pp/campaign_manifest.json" ]] && grep -q '"status": "READY"' "${RUN_ROOT}/pp/campaign_manifest.json" && pp_ready=1
    [[ -s "${RUN_ROOT}/auau/campaign_manifest.json" ]] && grep -q '"status": "READY"' "${RUN_ROOT}/auau/campaign_manifest.json" && auau_ready=1
    status
    if [[ "$pp_ready" == "1" && "$auau_ready" == "1" ]]; then
      write_campaign_metadata "READY"
      say "READY: pp and auau manifests are complete"
      break
    fi
    sleep "$interval"
  done
}

case "$MODE" in
  plan) print_plan ;;
  preflight) preflight ;;
  submit) submit_campaign ;;
  submit-pp) submit_single_domain pp ;;
  submit-auau) submit_single_domain auau ;;
  status) setup_env; status ;;
  wait) wait_ready ;;
  *)
    die "unknown mode: $MODE (expected plan|preflight|submit|submit-pp|submit-auau|status|wait)"
    ;;
esac
