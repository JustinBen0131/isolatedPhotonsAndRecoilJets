#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

group="${RJ_TARGET80_GROUP:-bdt_target80_gated_20260512_001012}"
config_dir="${RJ_TARGET80_CONFIG_DIR:-${repo_root}/condor_generated_configs/bdt_target80_all_available_20260511_224158}"
merge_group="${RJ_SIM_MERGE_GROUP_SIZE:-75}"
merge_memory="${RJ_SIM_FIRSTROUND_REQUEST_MEMORY:-8000MB}"
retry_cap="${RJ_AUTO_MEMORY_RETRY_CAP_MB:-16000}"
poll_seconds="${RJ_TARGET80_MERGE_POLL_SECONDS:-300}"
do_run="${RJ_TARGET80_MERGE_DO_RUN:-0}"
loop="${RJ_TARGET80_MERGE_LOOP:-0}"
max_configs="${RJ_TARGET80_MERGE_MAX_CONFIGS:-999}"
merge_cfg_match="${MERGE_CFG_MATCH:-}"

log_dir="${RJ_TARGET80_MERGE_LOG_DIR:-${repo_root}/condor_generated_configs/${group}}"
mkdir -p "$log_dir"
log="${log_dir}/target80_merge_ready_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$log") 2>&1

echo "RECOILJETS_TARGET80_READY_MERGE_DRIVER_V1"
echo "host=$(hostname -f 2>/dev/null || hostname)"
echo "repo_root=${repo_root}"
echo "group=${group}"
echo "config_dir=${config_dir}"
echo "merge_group=${merge_group}"
echo "merge_memory=${merge_memory}"
echo "retry_cap=${retry_cap}"
echo "do_run=${do_run}"
echo "loop=${loop}"
echo "max_configs=${max_configs}"
echo "poll_seconds=${poll_seconds}"
echo "merge_cfg_match=${merge_cfg_match}"
echo "log=${log}"
echo

[[ -d "$config_dir" ]] || { echo "[ERROR] Missing config_dir: ${config_dir}" >&2; exit 2; }

queue_cache="$(mktemp "${TMPDIR:-/tmp}/target80_merge_q.XXXXXX")"
trap 'rm -f "$queue_cache"' EXIT

refresh_queue_cache() {
  condor_q "${USER:-patsfan753}" -af JobStatus Args > "$queue_cache" 2>/dev/null || : > "$queue_cache"
}

count_active_matching() {
  local pattern="$1"
  awk -v pat="$pattern" '$1 != 3 && $1 != 4 && index($0, pat) {n++} END {print n+0}' "$queue_cache"
}

count_held_matching() {
  local pattern="$1"
  awk -v pat="$pattern" '$1 == 5 && index($0, pat) {n++} END {print n+0}' "$queue_cache"
}

count_roots() {
  local dir="$1"
  find "$dir" -type f -name "*.root" 2>/dev/null | wc -l | awk '{print $1+0}'
}

count_finals() {
  local dir="$1"
  find "$dir" -type f \( \
    -path "*/photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root" -o \
    -path "*/embeddedJet12and20merged_SIM/RecoilJets_embeddedJet12plus20_MERGED.root" -o \
    -path "*/embeddedJet12and20and30merged_SIM/RecoilJets_embeddedJet12plus20plus30_MERGED.root" \
  \) 2>/dev/null | wc -l | awk '{print $1+0}'
}

count_cfg_rows_from_yaml() {
  local yaml="$1"
  awk '
    /^[[:space:]]*photon_id_sets:[[:space:]]*$/ { in_sets=1; next }
    in_sets && /^[^[:space:]-]/ { in_sets=0 }
    in_sets && /^[[:space:]]*-[[:space:]]*\[/ { n++ }
    END { print n+0 }
  ' "$yaml"
}

count_cfg_dirs() {
  local dir="$1"
  find "$dir" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l | awk '{print $1+0}'
}

wait_for_pattern_clear() {
  local label="$1"
  local pattern="$2"
  while true; do
    refresh_queue_cache
    local active held
    active="$(count_active_matching "$pattern")"
    held="$(count_held_matching "$pattern")"
    echo "[wait] ${label}: active=${active} held=${held} pattern=${pattern}"
    if (( held > 0 )); then
      echo "[ERROR] Held jobs appeared while waiting for ${label}; stopping." >&2
      condor_q -nobatch "${USER:-patsfan753}" 2>/dev/null | tail -80 || true
      exit 11
    fi
    (( active == 0 )) && break
    sleep "$poll_seconds"
  done
}

run_merge_stage() {
  local tag="$1"
  local yaml="$2"
  local dataset="$3"
  local stage="$4"
  local input_base="$5"
  local merge_base="$6"
  local -a stage_args=( "$stage" )

  if [[ "$stage" == "secondRound" || "$stage" == "finalStitch" ]]; then
    stage_args+=( "condor" )
  fi

  echo
  echo "====================================================================="
  echo "MERGE tag=${tag} dataset=${dataset} stage=${stage}"
  echo "  yaml       : ${yaml}"
  echo "  input base : ${input_base}"
  echo "  output base: ${merge_base}"
  echo "====================================================================="

  [[ -f "$yaml" ]] || { echo "[ERROR] Missing YAML: ${yaml}" >&2; exit 20; }
  [[ -d "$input_base" ]] || { echo "[ERROR] Missing input base: ${input_base}" >&2; exit 21; }

  env \
    RJ_SIMEMBEDDEDINCLUSIVE_THREE_SAMPLES="${RJ_SIMEMBEDDEDINCLUSIVE_THREE_SAMPLES:-0}" \
    RJ_SIMEMBEDDEDINCLUSIVE_INCLUDE_JET30="${RJ_SIMEMBEDDEDINCLUSIVE_INCLUDE_JET30:-0}" \
    MERGE_CONFIG_YAML="${yaml}" \
    MERGE_CFG_MATCH="${MERGE_CFG_MATCH:-}" \
    MERGE_SIM_INPUT_BASE_OVERRIDE="${input_base}" \
    MERGE_OUT_BASE_OVERRIDE="${merge_base}" \
    RJ_SIM_FIRSTROUND_REQUEST_MEMORY="${merge_memory}" \
    RJ_AUTO_MEMORY_RETRY_CAP_MB="${retry_cap}" \
    ./scripts/mergeRecoilJets.sh "${dataset}" "${stage_args[@]}" groupSize "${merge_group}"

  wait_for_pattern_clear "${tag} ${dataset} ${stage}" "${merge_base}"
}

process_ready_config() {
  local yaml="$1"
  local name tag signal_bulk background_bulk merge_base active held sigraw bkgraw finals
  local cfg_rows expected_finals sigcfgs bkgcfgs
  name="$(basename "$yaml" .yaml)"
  tag="${group}_${name}"
  signal_bulk="/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded_${tag}"
  background_bulk="/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembeddedinclusive_${tag}"
  merge_base="/sphenix/u/patsfan753/scratch/thesisAnalysis/output_${tag}"

  cfg_rows="$(count_cfg_rows_from_yaml "$yaml")"
  expected_finals=$((2 * cfg_rows))
  active="$(count_active_matching "$tag")"
  held="$(count_held_matching "$tag")"
  sigraw="$(count_roots "$signal_bulk")"
  bkgraw="$(count_roots "$background_bulk")"
  sigcfgs="$(count_cfg_dirs "$signal_bulk")"
  bkgcfgs="$(count_cfg_dirs "$background_bulk")"
  finals="$(count_finals "$merge_base")"

  printf "%-48s active=%6s held=%3s sigcfg=%3s/%-3s bkgcfg=%3s/%-3s sigraw=%7s bkgraw=%7s finals=%4s/%-4s" \
    "$name" "$active" "$held" "$sigcfgs" "$cfg_rows" "$bkgcfgs" "$cfg_rows" "$sigraw" "$bkgraw" "$finals" "$expected_finals"

  if (( held > 0 )); then
    echo "  -> HELD: inspect, do not merge"
    return 2
  fi
  if (( active > 0 )); then
    echo "  -> WAIT: matching jobs still active"
    return 1
  fi
  if (( cfg_rows == 0 )); then
    echo "  -> BLOCKED: no photon_id_sets rows parsed from YAML"
    return 2
  fi
  if (( finals >= expected_finals )); then
    echo "  -> DONE/SKIP: all expected merged outputs exist"
    return 1
  fi
  if (( sigcfgs < cfg_rows || bkgcfgs < cfg_rows || sigraw == 0 || bkgraw == 0 )); then
    echo "  -> WAIT: both raw datasets are not complete yet"
    return 1
  fi

  if (( finals > 0 )); then
    echo "  -> READY: partial final outputs exist; rerunning merge stages to fill missing finals"
  else
    echo "  -> READY"
  fi
  if [[ "$do_run" != "1" ]]; then
    return 0
  fi

  run_merge_stage "$tag" "$yaml" isSimEmbedded firstRound "$signal_bulk" "$merge_base"
  run_merge_stage "$tag" "$yaml" isSimEmbedded secondRound "${merge_base}/simembedded" "$merge_base"
  run_merge_stage "$tag" "$yaml" isSimEmbedded finalStitch "${merge_base}/simembedded" "$merge_base"

  run_merge_stage "$tag" "$yaml" isSimEmbeddedInclusive firstRound "$background_bulk" "$merge_base"
  run_merge_stage "$tag" "$yaml" isSimEmbeddedInclusive secondRound "${merge_base}/simembeddedinclusive" "$merge_base"
  run_merge_stage "$tag" "$yaml" isSimEmbeddedInclusive finalStitch "${merge_base}/simembeddedinclusive" "$merge_base"

  local final_after
  final_after="$(count_finals "$merge_base")"
  echo "[done] ${name}: final_count=${final_after} output=${merge_base}"
}

scan_once() {
  local processed=0 blockers=0 ready=0
  refresh_queue_cache
  echo
  echo "== target80 merge readiness scan $(date) =="
  shopt -s nullglob
  for yaml in "${config_dir}"/analysis_config_*_target80.yaml; do
    set +e
    process_ready_config "$yaml"
    rc=$?
    set -e
    case "$rc" in
      0) (( ready+=1 )); if [[ "$do_run" == "1" ]]; then (( processed+=1 )); fi ;;
      1) : ;;
      2) (( blockers+=1 )) ;;
      *) echo "[ERROR] Unexpected status ${rc} for ${yaml}" >&2; exit "$rc" ;;
    esac
    if [[ "$do_run" == "1" ]] && (( processed >= max_configs )); then
      echo "[limit] processed=${processed}; stopping due RJ_TARGET80_MERGE_MAX_CONFIGS=${max_configs}"
      break
    fi
  done
  shopt -u nullglob
  echo "summary_ready=${ready} summary_processed=${processed} summary_blockers=${blockers}"
  (( blockers == 0 ))
}

if [[ "$loop" == "1" ]]; then
  while true; do
    scan_once
    echo "[loop] sleeping ${poll_seconds}s"
    sleep "$poll_seconds"
  done
else
  scan_once
fi

echo
echo "TARGET80_MERGE_LOG=${log}"
