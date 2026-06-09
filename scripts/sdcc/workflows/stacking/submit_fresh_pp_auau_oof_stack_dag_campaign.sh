#!/usr/bin/env bash
set -euo pipefail

MODE="${1:-plan}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
RJ_REPO_BASE="${RJ_REPO_BASE:-$(cd "${SCRIPT_DIR}/../../../.." && pwd -P)}"
ML_PYTHON="${RJ_ML_PYTHON:-/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python}"
STAMP="${RJ_FRESH_STACK_DAG_STAMP:-fresh_pp_auau_stack_dag_$(date +%Y%m%d_%H%M%S)}"
MODEL_BASE="${RJ_FRESH_STACK_DAG_MODEL_BASE:-/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models}"
RUN_ROOT="${RJ_FRESH_STACK_DAG_RUN_ROOT:-${MODEL_BASE}/fresh_pp_auau_bdt_mlp_oof_stack_dag_${STAMP}}"
SUB_ROOT="${RJ_FRESH_STACK_DAG_SUB_ROOT:-${RJ_REPO_BASE}/condor_sub/freshOOFStackDAG_${STAMP}}"
LOG_ROOT="${RJ_FRESH_STACK_DAG_LOG_ROOT:-${RUN_ROOT}/logs}"

AUAU_SOURCE="${RJ_FRESH_OOF_AUAU_SOURCE:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/THE8_branchA_ladder_jet12_20_30_40_20260527}"
AUAU_MANIFEST="${RJ_FRESH_OOF_AUAU_MANIFEST:-${AUAU_SOURCE}/manifests/training_roots.list}"

PP_RUN_ROOT="${RJ_FRESH_OOF_PP_RUN_ROOT:-/sphenix/u/patsfan753/scratch/thesisAnalysis/dataOutput/ppPhotonMLPipeline/ppg12_matched_20260521_230310}"
PP_TRAIN_MANIFEST="${RJ_FRESH_OOF_PP_TRAIN_MANIFEST:-${PP_RUN_ROOT}/training_roots_currentIAN.list}"
PP_SIGNAL_MANIFEST="${RJ_FRESH_OOF_PP_SIGNAL_MANIFEST:-${PP_RUN_ROOT}/signal_roots_currentIAN.list}"
PP_INCLUSIVE_MANIFEST="${RJ_FRESH_OOF_PP_INCLUSIVE_MANIFEST:-${PP_RUN_ROOT}/inclusive_roots_currentIAN.list}"

PP_TREE="${RJ_FRESH_OOF_PP_TREE:-AuAuPhotonIDTrainingTree}"
AUAU_TREE="${RJ_FRESH_OOF_AUAU_TREE:-AuAuPhotonIDTrainingTree}"
REQUEST_CPUS_TRAIN="${RJ_FRESH_STACK_DAG_REQUEST_CPUS_TRAIN:-2}"
REQUEST_CPUS_LIGHT="${RJ_FRESH_STACK_DAG_REQUEST_CPUS_LIGHT:-1}"
PP_MATRIX_MEMORY="${RJ_FRESH_STACK_DAG_PP_MATRIX_MEMORY:-64000MB}"
AUAU_MATRIX_MEMORY="${RJ_FRESH_STACK_DAG_AUAU_MATRIX_MEMORY:-48000MB}"
PP_BDT_MEMORY="${RJ_FRESH_STACK_DAG_PP_BDT_MEMORY:-64000MB}"
AUAU_BDT_MEMORY="${RJ_FRESH_STACK_DAG_AUAU_BDT_MEMORY:-48000MB}"
PP_MLP_MEMORY="${RJ_FRESH_STACK_DAG_PP_MLP_MEMORY:-32000MB}"
AUAU_MLP_MEMORY="${RJ_FRESH_STACK_DAG_AUAU_MLP_MEMORY:-24000MB}"
SCORE_MEMORY="${RJ_FRESH_STACK_DAG_SCORE_MEMORY:-12000MB}"
STACK_MEMORY="${RJ_FRESH_STACK_DAG_STACK_MEMORY:-16000MB}"
REDUCE_MEMORY="${RJ_FRESH_STACK_DAG_REDUCE_MEMORY:-24000MB}"
PREDICT_CHUNK_ROWS="${RJ_FRESH_STACK_DAG_PREDICT_CHUNK_ROWS:-100000}"
BDT_ESTIMATORS="${RJ_FRESH_OOF_BDT_ESTIMATORS:-750}"
STACK_MLP_EPOCHS="${RJ_FRESH_OOF_STACK_MLP_EPOCHS:-180}"
STACK_MLP_HIDDEN="${RJ_FRESH_OOF_STACK_MLP_HIDDEN:-256,128,64}"
STACK_MLP_BATCH_SIZE="${RJ_FRESH_OOF_STACK_MLP_BATCH_SIZE:-8192}"

say() { printf '\033[1;36m[freshOOFStackDAG]\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[freshOOFStackDAG][WARN]\033[0m %s\n' "$*" >&2; }
die() { printf '\033[1;31m[freshOOFStackDAG][ERR]\033[0m %s\n' "$*" >&2; exit 2; }

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
  export OMP_NUM_THREADS="${REQUEST_CPUS_TRAIN}"
  export OPENBLAS_NUM_THREADS="${REQUEST_CPUS_TRAIN}"
  export MKL_NUM_THREADS="${REQUEST_CPUS_TRAIN}"
  export NUMEXPR_NUM_THREADS="${REQUEST_CPUS_TRAIN}"
  export XGBOOST_NUM_THREADS="${REQUEST_CPUS_TRAIN}"
  export MALLOC_ARENA_MAX="${RJ_FRESH_OOF_MALLOC_ARENA_MAX:-2}"
  unset PYTHONHOME
}

print_plan() {
  say "schema=RJ_FRESH_PP_AUAU_BDT_MLP_STACK_DAG_CAMPAIGN_V1"
  say "mode=$MODE"
  say "stamp=$STAMP"
  say "repo=$RJ_REPO_BASE"
  say "run_root=$RUN_ROOT"
  say "sub_root=$SUB_ROOT"
  say "ml_python=$ML_PYTHON"
  say "stage_driver=${RJ_REPO_BASE}/scripts/ml/stacking/train_photon_bdt_mlp_oof_stack_staged.py"
  say "pp_manifest=$PP_TRAIN_MANIFEST"
  say "pp_overlay_manifests=$PP_SIGNAL_MANIFEST $PP_INCLUSIVE_MANIFEST"
  say "auau_manifest=$AUAU_MANIFEST"
  say "memory matrix pp/auau=$PP_MATRIX_MEMORY/$AUAU_MATRIX_MEMORY bdt pp/auau=$PP_BDT_MEMORY/$AUAU_BDT_MEMORY mlp pp/auau=$PP_MLP_MEMORY/$AUAU_MLP_MEMORY score=$SCORE_MEMORY stack=$STACK_MEMORY reduce=$REDUCE_MEMORY"
  say "throttle categories: MAXJOBS pp_bdt=2 pp_mlp=3 auau_bdt=2"
}

need_file() {
  [[ -s "$1" ]] || die "missing required file: $1"
}

require_inputs() {
  need_file "${RJ_REPO_BASE}/scripts/ml/stacking/train_photon_bdt_mlp_oof_stack.py"
  need_file "${RJ_REPO_BASE}/scripts/ml/stacking/train_photon_bdt_mlp_oof_stack_staged.py"
  need_file "$AUAU_MANIFEST"
  need_file "$PP_TRAIN_MANIFEST"
  need_file "$PP_SIGNAL_MANIFEST"
  need_file "$PP_INCLUSIVE_MANIFEST"
}

require_submit_provenance() {
  [[ "${RJ_FRESH_STACK_DAG_APPROVED:-0}" == "1" ]] || die "submit requires RJ_FRESH_STACK_DAG_APPROVED=1"
  [[ -n "${RJ_CODEX_CHAT_NAME:-}" ]] || die "submit requires RJ_CODEX_CHAT_NAME"
  [[ -n "${RJ_CODEX_THREAD_ID:-}" ]] || die "submit requires RJ_CODEX_THREAD_ID"
}

domain_base_args() {
  local domain="$1" outdir="$2" matrix_dir="$3"
  case "$domain" in
    pp)
      printf '%q ' \
        --domain pp \
        --input "@${PP_TRAIN_MANIFEST}" \
        --overlay-input "@${PP_SIGNAL_MANIFEST}" "@${PP_INCLUSIVE_MANIFEST}" \
        --tree "$PP_TREE" \
        --outdir "$outdir" \
        --matrix-dir "$matrix_dir" \
        --feature-preset pp_basev3e_noiso \
        --pt-range 5:35 \
        --centrality-range=-1:0 \
        --weight-mode ppg12-exact \
        --ppg12-exact-expected-samples run28_photonjet5,run28_photonjet10,run28_photonjet20,run28_jet8,run28_jet12,run28_jet20,run28_jet30 \
        --training-inclusive-sources run28_jet8,run28_jet12,run28_jet20,run28_jet30 \
        --overlay-signal-sources run28_photonjet5,run28_photonjet10,run28_photonjet20 \
        --overlay-inclusive-sources run28_jet8,run28_jet12,run28_jet20,run28_jet30,run28_jet40 \
        --locked-test-fraction 0.20 \
        --folds 5 \
        --random-seed 260604 \
        --report-et-bins 5,10,14,18,22,35 \
        --skip-missing-tree
      ;;
    auau)
      printf '%q ' \
        --domain auau \
        --input "@${AUAU_MANIFEST}" \
        --overlay-input "@${AUAU_MANIFEST}" \
        --tree "$AUAU_TREE" \
        --outdir "$outdir" \
        --matrix-dir "$matrix_dir" \
        --feature-preset auau_global_noiso \
        --pt-range 15:35 \
        --centrality-range 0:80 \
        --weight-mode ppg12-exact \
        --ppg12-exact-expected-samples run28_embeddedPhoton12,run28_embeddedPhoton20,run28_embeddedJet12,run28_embeddedJet20,run28_embeddedJet30,run28_embeddedJet40 \
        --training-inclusive-sources run28_embeddedJet12,run28_embeddedJet20,run28_embeddedJet30,run28_embeddedJet40 \
        --overlay-signal-sources run28_embeddedPhoton12,run28_embeddedPhoton20 \
        --overlay-inclusive-sources run28_embeddedJet12,run28_embeddedJet20,run28_embeddedJet30,run28_embeddedJet40 \
        --locked-test-fraction 0.20 \
        --folds 5 \
        --random-seed 260604 \
        --report-et-bins 15,17,19,21,23,25,27,30,35 \
        --report-cent-bins 0,20,40,60,80 \
        --skip-missing-tree
      ;;
    *) die "unknown domain: $domain" ;;
  esac
}

lane_arg() {
  local lane="$1"
  case "$lane" in
    current_oof) printf '%q ' --stack-training-mode oof5 --stack-validation-fold 0 ;;
    simple_holdout) printf '%q ' --stack-training-mode simple_holdout --simple-stack-fold 0 --stack-validation-fold 0 ;;
    *) die "unknown lane: $lane" ;;
  esac
}

node_memory() {
  local domain="$1" kind="$2"
  case "${domain}:${kind}" in
    pp:matrix) echo "$PP_MATRIX_MEMORY" ;;
    auau:matrix) echo "$AUAU_MATRIX_MEMORY" ;;
    pp:bdt) echo "$PP_BDT_MEMORY" ;;
    auau:bdt) echo "$AUAU_BDT_MEMORY" ;;
    pp:mlp) echo "$PP_MLP_MEMORY" ;;
    auau:mlp) echo "$AUAU_MLP_MEMORY" ;;
    *:score) echo "$SCORE_MEMORY" ;;
    *:stack) echo "$STACK_MEMORY" ;;
    *:reduce) echo "$REDUCE_MEMORY" ;;
    *) die "unknown memory kind ${domain}:${kind}" ;;
  esac
}

node_cpus() {
  local kind="$1"
  case "$kind" in
    matrix|bdt|mlp|stack) echo "$REQUEST_CPUS_TRAIN" ;;
    score|reduce) echo "$REQUEST_CPUS_LIGHT" ;;
    *) die "unknown cpu kind: $kind" ;;
  esac
}

write_node() {
  local node="$1" domain="$2" kind="$3" command="$4"
  local mem cpus worker sub
  mem="$(node_memory "$domain" "$kind")"
  cpus="$(node_cpus "$kind")"
  worker="${SUB_ROOT}/nodes/${node}.sh"
  sub="${SUB_ROOT}/nodes/${node}.sub"
  mkdir -p "${SUB_ROOT}/nodes" "${SUB_ROOT}/node_logs"
  cat > "$worker" <<EOF
#!/usr/bin/env bash
set -euo pipefail
cd "$RJ_REPO_BASE"
ML_PYTHON="$ML_PYTHON"
REQUEST_CPUS_TRAIN="$cpus"
export RJ_CODEX_CHAT_NAME="${RJ_CODEX_CHAT_NAME:-}"
export RJ_CODEX_THREAD_ID="${RJ_CODEX_THREAD_ID:-}"
$(declare -f setup_env)
setup_env
export OMP_NUM_THREADS="$cpus"
export OPENBLAS_NUM_THREADS="$cpus"
export MKL_NUM_THREADS="$cpus"
export NUMEXPR_NUM_THREADS="$cpus"
export XGBOOST_NUM_THREADS="$cpus"
$command
EOF
  chmod +x "$worker"
  cat > "$sub" <<EOF
universe = vanilla
executable = ${worker}
output = ${SUB_ROOT}/node_logs/${node}.out
error = ${SUB_ROOT}/node_logs/${node}.err
log = ${SUB_ROOT}/node_logs/${node}.log
request_cpus = ${cpus}
request_memory = ${mem}
notification = Never
getenv = True
+JobBatchName = "fresh_stack_dag_${STAMP}_${node}"
queue
EOF
}

append_node() {
  local node="$1" sub="$2"
  printf 'JOB %s %s\n' "$node" "$sub" >> "${SUB_ROOT}/fresh_stack.dag"
}

append_parent() {
  local parent="$1" child="$2"
  printf 'PARENT %s CHILD %s\n' "$parent" "$child" >> "${SUB_ROOT}/fresh_stack.dag"
}

append_category() {
  local node="$1" category="$2"
  printf 'CATEGORY %s %s\n' "$node" "$category" >> "${SUB_ROOT}/fresh_stack.dag"
}

stage_cmd() {
  local stage="$1" args="$2"
  printf '%q ' "$ML_PYTHON" "${RJ_REPO_BASE}/scripts/ml/stacking/train_photon_bdt_mlp_oof_stack_staged.py" --stage "$stage"
  printf '%s ' "$args"
  printf '%q ' --predict-chunk-rows "$PREDICT_CHUNK_ROWS" --bdt-estimators "$BDT_ESTIMATORS" --stack-mlp-hidden "$STACK_MLP_HIDDEN" --stack-mlp-epochs "$STACK_MLP_EPOCHS" --stack-mlp-batch-size "$STACK_MLP_BATCH_SIZE" --n-jobs "$REQUEST_CPUS_TRAIN"
}

write_campaign_manifest() {
  local status="$1"
  mkdir -p "$RUN_ROOT"
  cat > "${RUN_ROOT}/dag_submit_manifest.json" <<EOF
{
  "schema": "RJ_FRESH_PP_AUAU_BDT_MLP_STACK_DAG_SUBMIT_V1",
  "status": "${status}",
  "stamp": "${STAMP}",
  "run_root": "${RUN_ROOT}",
  "sub_root": "${SUB_ROOT}",
  "repo": "${RJ_REPO_BASE}",
  "ml_python": "${ML_PYTHON}",
  "execution_model": "DAGMan staged matrix/base-score/stack/reduce",
  "codex_chat_name": "${RJ_CODEX_CHAT_NAME:-}",
  "codex_thread_id": "${RJ_CODEX_THREAD_ID:-}",
  "lanes": ["current_oof", "simple_holdout"],
  "domains": ["pp", "auau"],
  "memory_contract": {
    "pp_matrix": "${PP_MATRIX_MEMORY}",
    "auau_matrix": "${AUAU_MATRIX_MEMORY}",
    "pp_bdt": "${PP_BDT_MEMORY}",
    "auau_bdt": "${AUAU_BDT_MEMORY}",
    "pp_mlp": "${PP_MLP_MEMORY}",
    "auau_mlp": "${AUAU_MLP_MEMORY}",
    "score": "${SCORE_MEMORY}",
    "stack": "${STACK_MEMORY}",
    "reduce": "${REDUCE_MEMORY}"
  }
}
EOF
}

build_dag() {
  setup_env
  print_plan
  require_inputs
  mkdir -p "$RUN_ROOT" "$SUB_ROOT" "$LOG_ROOT"
  : > "${SUB_ROOT}/fresh_stack.dag"
  printf 'MAXJOBS pp_bdt 2\nMAXJOBS pp_mlp 3\nMAXJOBS auau_bdt 2\nMAXJOBS matrix 2\n\n' >> "${SUB_ROOT}/fresh_stack.dag"

  local domain lane domain_out matrix_out matrix_dir base_args lane_args matrix_node cmd
  for domain in pp auau; do
    matrix_out="${RUN_ROOT}/matrix_store/${domain}"
    matrix_dir="${matrix_out}/matrix"
    base_args="$(domain_base_args "$domain" "$matrix_out" "$matrix_dir")"
    matrix_node="matrix_${domain}"
    cmd="$(stage_cmd build-matrix "${base_args} $(lane_arg current_oof)")"
    write_node "$matrix_node" "$domain" matrix "$cmd"
    append_node "$matrix_node" "${SUB_ROOT}/nodes/${matrix_node}.sub"
    append_category "$matrix_node" matrix
  done

  for lane in current_oof simple_holdout; do
    for domain in pp auau; do
      domain_out="${RUN_ROOT}/${lane}/${domain}"
      matrix_dir="${RUN_ROOT}/matrix_store/${domain}/matrix"
      base_args="$(domain_base_args "$domain" "$domain_out" "$matrix_dir")"
      lane_args="$(lane_arg "$lane")"
      local score_nodes=()
      if [[ "$lane" == "current_oof" ]]; then
        for fold in 0 1 2 3 4; do
          for kind in bdt mlp; do
            local train_node="train_${lane}_${domain}_fold${fold}_${kind}"
            local score_node="score_${lane}_${domain}_fold${fold}_${kind}"
            cmd="$(stage_cmd train-base "${base_args} ${lane_args} --model-role fold --fold ${fold} --model-kind ${kind}")"
            write_node "$train_node" "$domain" "$kind" "$cmd"
            append_node "$train_node" "${SUB_ROOT}/nodes/${train_node}.sub"
            append_parent "matrix_${domain}" "$train_node"
            append_category "$train_node" "${domain}_${kind}"
            cmd="$(stage_cmd score-base "${base_args} ${lane_args} --model-role fold --fold ${fold} --model-kind ${kind} --score-region fold")"
            write_node "$score_node" "$domain" score "$cmd"
            append_node "$score_node" "${SUB_ROOT}/nodes/${score_node}.sub"
            append_parent "$train_node" "$score_node"
            score_nodes+=("$score_node")
          done
        done
      else
        for kind in bdt mlp; do
          local train_node="train_${lane}_${domain}_simple_${kind}"
          local score_node="score_${lane}_${domain}_simple_${kind}"
          cmd="$(stage_cmd train-base "${base_args} ${lane_args} --model-role simple --model-kind ${kind}")"
          write_node "$train_node" "$domain" "$kind" "$cmd"
          append_node "$train_node" "${SUB_ROOT}/nodes/${train_node}.sub"
          append_parent "matrix_${domain}" "$train_node"
          append_category "$train_node" "${domain}_${kind}"
          cmd="$(stage_cmd score-base "${base_args} ${lane_args} --model-role simple --model-kind ${kind} --score-region simple_holdout")"
          write_node "$score_node" "$domain" score "$cmd"
          append_node "$score_node" "${SUB_ROOT}/nodes/${score_node}.sub"
          append_parent "$train_node" "$score_node"
          score_nodes+=("$score_node")
        done
      fi
      for kind in bdt mlp; do
        local train_node="train_${lane}_${domain}_final_${kind}"
        local score_node="score_${lane}_${domain}_final_${kind}"
        cmd="$(stage_cmd train-base "${base_args} ${lane_args} --model-role final --model-kind ${kind}")"
        write_node "$train_node" "$domain" "$kind" "$cmd"
        append_node "$train_node" "${SUB_ROOT}/nodes/${train_node}.sub"
        append_parent "matrix_${domain}" "$train_node"
        append_category "$train_node" "${domain}_${kind}"
        cmd="$(stage_cmd score-base "${base_args} ${lane_args} --model-role final --model-kind ${kind} --score-region locked_test")"
        write_node "$score_node" "$domain" score "$cmd"
        append_node "$score_node" "${SUB_ROOT}/nodes/${score_node}.sub"
        append_parent "$train_node" "$score_node"
        score_nodes+=("$score_node")
      done
      local stack_node="stack_${lane}_${domain}"
      local reduce_node="reduce_${lane}_${domain}"
      cmd="$(stage_cmd train-stack "${base_args} ${lane_args}")"
      write_node "$stack_node" "$domain" stack "$cmd"
      append_node "$stack_node" "${SUB_ROOT}/nodes/${stack_node}.sub"
      for parent_node in "${score_nodes[@]}"; do
        append_parent "$parent_node" "$stack_node"
      done
      cmd="$(stage_cmd reduce-domain "${base_args} ${lane_args}")"
      write_node "$reduce_node" "$domain" reduce "$cmd"
      append_node "$reduce_node" "${SUB_ROOT}/nodes/${reduce_node}.sub"
      append_parent "$stack_node" "$reduce_node"
    done
  done
  write_campaign_manifest "DAG_WRITTEN"
  say "dag=${SUB_ROOT}/fresh_stack.dag"
  say "node_count=$(grep -c '^JOB ' "${SUB_ROOT}/fresh_stack.dag")"
}

preflight() {
  build_dag
  bash -n "${SUB_ROOT}"/nodes/*.sh
  "$ML_PYTHON" "${RJ_REPO_BASE}/scripts/ml/stacking/train_photon_bdt_mlp_oof_stack_staged.py" \
    --stage build-matrix \
    --domain pp \
    --self-test \
    --self-test-rows 1000 \
    --outdir "${RUN_ROOT}/selftest_matrix/pp" \
    --feature-preset pp_basev3e_noiso \
    --pt-range 5:35 \
    --centrality-range=-1:0 \
    --weight-mode none \
    --folds 3 \
    --locked-test-fraction 0.20
  write_campaign_manifest "PREFLIGHT_PASS"
}

submit_dag() {
  setup_env
  require_submit_provenance
  build_dag
  write_campaign_manifest "SUBMITTING"
  (cd "$SUB_ROOT" && condor_submit_dag -force fresh_stack.dag) | tee "${SUB_ROOT}/dag.submit"
  write_campaign_manifest "SUBMITTED"
  status
}

cluster_from_submit() {
  awk '/submitted to cluster/ {gsub("\\.","",$NF); print $NF; exit}' "${SUB_ROOT}/dag.submit" 2>/dev/null || true
}

status() {
  print_plan
  if [[ -s "${SUB_ROOT}/dag.submit" ]]; then
    local cluster
    cluster="$(cluster_from_submit)"
    say "dag_cluster=${cluster:-unknown}"
    [[ -n "$cluster" ]] && condor_q "$cluster" -nobatch 2>/dev/null || true
  else
    warn "no DAG submit record yet: ${SUB_ROOT}/dag.submit"
  fi
  for lane in current_oof simple_holdout; do
    for domain in pp auau; do
      local root="${RUN_ROOT}/${lane}/${domain}"
      if [[ -s "${root}/campaign_manifest.json" ]]; then
        say "READY? ${root}/campaign_manifest.json"
      else
        warn "missing final manifest: ${root}/campaign_manifest.json"
      fi
    done
  done
}

case "$MODE" in
  plan) print_plan ;;
  preflight) preflight ;;
  write-dag) build_dag ;;
  submit) submit_dag ;;
  status) setup_env; status ;;
  *)
    die "unknown mode: $MODE (expected plan|preflight|write-dag|submit|status)"
    ;;
esac
