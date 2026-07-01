#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
if [[ -f "${script_dir}/../../runtime/io/recoiljets_io_paths.sh" ]]; then
  source "${script_dir}/../../runtime/io/recoiljets_io_paths.sh"
elif [[ -f scripts/sdcc/runtime/io/recoiljets_io_paths.sh ]]; then
  source scripts/sdcc/runtime/io/recoiljets_io_paths.sh
else
  echo "[ERROR] Cannot locate scripts/sdcc/runtime/io/recoiljets_io_paths.sh" >&2
  exit 2
fi
repo_root="$(rj_find_repo_root "$script_dir")"
cd "$repo_root"

campaign_tag="${RJ_THE79_RAA_CAMPAIGN_TAG:-the79_phenix_like_raa_bdt98_wp80_$(date +%Y%m%d_%H%M)}"
rj_validate_campaign_tag "$campaign_tag" || exit 2

family="${RJ_THE79_RAA_FAMILY:-phenix_like_raa}"
pp_yaml="${RJ_THE79_RAA_PP_YAML:-macros/analysis_config_the79_phenix_like_raa_pp_ppg12.yaml}"
auau_yaml="${RJ_THE79_RAA_AUAU_YAML:-macros/analysis_config_the79_phenix_like_raa_auau_bdt98_wp80.yaml}"
preflight_py="${RJ_THE79_RAA_PREFLIGHT:-scripts/diagnostics/ml_validation/check_the79_phenix_like_raa_inputs.py}"
model_root="${RJ_THE79_AUAU_MODEL_ROOT:-/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/the79_fullsrc_binned14_ppg12pt5to40_20260627_2052}"
taa_table="${RJ_THE79_TAA_TABLE:-}"

[[ -f "$pp_yaml" ]] || { echo "[ERROR] pp YAML not found: $pp_yaml" >&2; exit 2; }
[[ -f "$auau_yaml" ]] || { echo "[ERROR] AuAu YAML not found: $auau_yaml" >&2; exit 2; }
[[ -f "$preflight_py" ]] || { echo "[ERROR] preflight helper not found: $preflight_py" >&2; exit 2; }

bulk_root="$(rj_recoiljets_bulk_root)/${family}/${campaign_tag}"
pp_bulk="${RJ_THE79_PP_BULK:-${bulk_root}/pp}"
auau_bulk="${RJ_THE79_AUAU_BULK:-${bulk_root}/auau}"
signal_bulk="${RJ_THE79_SIGNAL_BULK:-$(rj_bulk_signal_root_for_campaign "$family" "$campaign_tag")}"
background_bulk="${RJ_THE79_BACKGROUND_BULK:-$(rj_bulk_background_root_for_campaign "$family" "$campaign_tag")}"
merge_base="${RJ_THE79_MERGE_BASE:-$(rj_merge_root_for_campaign "$campaign_tag")}"

log_dir="${RJ_THE79_RAA_LOG_DIR:-${repo_root}/condor_generated_configs/${campaign_tag}}"
mkdir -p "$log_dir"
submit_log="${log_dir}/submit_${campaign_tag}.log"
exec > >(tee -a "$submit_log") 2>&1

notify="${RJ_NOTIFY_EMAILS:-just0131@gmail.com}"
pp_group="${RJ_THE79_PP_GROUP_SIZE:-500}"
auau_group="${RJ_THE79_AUAU_GROUP_SIZE:-7}"
sim_group="${RJ_THE79_SIM_GROUP_SIZE:-7}"
pp_memory="${RJ_THE79_PP_MEMORY:-3000MB}"
auau_memory="${RJ_THE79_AUAU_MEMORY:-12000MB}"
sim_memory="${RJ_THE79_SIM_MEMORY:-12000MB}"
merge_memory="${RJ_THE79_MERGE_MEMORY:-16000MB}"
merge_group="${RJ_THE79_MERGE_GROUP_SIZE:-75}"
allow_existing="${RJ_THE79_ALLOW_EXISTING:-0}"
do_run="${RJ_DO_RUN:-0}"

echo "RECOILJETS_THE79_PHENIX_LIKE_RAA_SUBMIT_V1"
echo "submit_host=$(hostname -f 2>/dev/null || hostname)"
echo "repo_root=${repo_root}"
echo "campaign_tag=${campaign_tag}"
echo "family=${family}"
echo "pp_yaml=${pp_yaml}"
echo "auau_yaml=${auau_yaml}"
echo "model_root=${model_root}"
echo "taa_table=${taa_table:-TAA_PENDING}"
echo "pp_bulk=${pp_bulk}"
echo "auau_bulk=${auau_bulk}"
echo "signal_bulk=${signal_bulk}"
echo "background_bulk=${background_bulk}"
echo "merge_base=${merge_base}"
echo "pp_group_size=${pp_group}"
echo "auau_group_size=${auau_group}"
echo "sim_group_size=${sim_group}"
echo "notify=${notify}"
echo "RJ_CODEX_CHAT_NAME=${RJ_CODEX_CHAT_NAME:-}"
echo "RJ_CODEX_THREAD_ID=${RJ_CODEX_THREAD_ID:-}"
echo "RJ_DAG_DRYRUN=${RJ_DAG_DRYRUN:-0}"
echo "RJ_DO_RUN=${do_run}"
echo

if [[ -z "${RJ_CODEX_CHAT_NAME:-}" || -z "${RJ_CODEX_THREAD_ID:-}" ]]; then
  echo "[ERROR] RJ_CODEX_CHAT_NAME and RJ_CODEX_THREAD_ID must be set before Codex-launched submit." >&2
  exit 3
fi

if condor_q "${USER:-patsfan753}" -af ClusterId JobStatus Args 2>/dev/null | grep -F "$campaign_tag" >/dev/null; then
  echo "[ERROR] Active Condor jobs already mention campaign tag ${campaign_tag}; refusing duplicate submit." >&2
  condor_q "${USER:-patsfan753}" -af ClusterId JobStatus Args 2>/dev/null | grep -F "$campaign_tag" || true
  exit 4
fi

check_empty_or_allowed() {
  local path="$1"
  local label="$2"
  if [[ -e "$path" ]] && [[ "$allow_existing" != "1" ]]; then
    echo "[ERROR] ${label} already exists: ${path}" >&2
    echo "[ERROR] Set RJ_THE79_ALLOW_EXISTING=1 only for an intentional resume/dry-run inspection." >&2
    exit 5
  fi
}
check_empty_or_allowed "$pp_bulk" "pp bulk output"
check_empty_or_allowed "$auau_bulk" "AuAu bulk output"
check_empty_or_allowed "$signal_bulk" "embedded signal bulk output"
check_empty_or_allowed "$background_bulk" "embedded background bulk output"
check_empty_or_allowed "$merge_base" "merge output"

preflight_out="${log_dir}/the79_phenix_like_raa_preflight.json"
preflight_args=(python3 "$preflight_py" --pp-yaml "$pp_yaml" --auau-yaml "$auau_yaml" --model-root "$model_root" --allow-taa-pending --out "$preflight_out")
if [[ -n "$taa_table" ]]; then
  preflight_args+=(--taa-table "$taa_table")
fi
"${preflight_args[@]}"

cat <<EOF

== submit contract ==
pp denominator:
  dataset      : isPP
  photon path  : PPG12 newPPG12 tight/sideband baseV3E
  env          : RJ_PPG12_PHOTON_YIELD=1, RJ_PHOTON_ID_ROW_MATCH=preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12
AuAu numerator:
  dataset      : isAuAu
  photon path  : newPPG12 + auauEtFineCent7BDT + auauBDTSideband
  AuAu BDT     : base14_perEtCent7, 98 BDTs, 8 <= E_T < 40 WP application
  non-tight    : relativeToTight, T80(c)-0.20 < score < T80(c)-0.03
embedded QA:
  signal       : isSimEmbedded Photon12+Photon20
  background   : isSimEmbeddedInclusive Jet12+Jet20+Jet30+Jet40 via RJ_SIMEMBEDDEDINCLUSIVE_FOUR_SAMPLES=1
normalization:
  TAA status   : ${taa_table:+READY CANDIDATE}${taa_table:-TAA_PENDING; final normalized R_AA blocked}
EOF

submit_data_pp() {
  echo
  echo "====================================================================="
  echo "Submitting pp denominator"
  echo "====================================================================="
  env \
    RJ_NOTIFY_EMAILS="$notify" \
    RJ_PROFILE_JOB=1 \
    RJ_CONFIG_YAML="$pp_yaml" \
    RJ_DEST_BASE_OVERRIDE="$pp_bulk" \
    RJ_MERGE_OUT_BASE_OVERRIDE="$merge_base" \
    RJ_REQUEST_MEMORY="$pp_memory" \
    RJ_AUTO_MERGE=1 \
    RJ_SIM_FIRSTROUND_REQUEST_MEMORY="$merge_memory" \
    RJ_SIM_MERGE_GROUP_SIZE="$merge_group" \
    RJ_ID_FANOUT_MAX_ROWS=1 \
    RJ_PPG12_PHOTON_YIELD=1 \
    RJ_PHOTON_ID_ROW_MATCH=preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12 \
    RJ_CODEX_CHAT_NAME="$RJ_CODEX_CHAT_NAME" \
    RJ_CODEX_THREAD_ID="$RJ_CODEX_THREAD_ID" \
    ./RecoilJets_Condor_submit.sh isPP condor all groupSize "$pp_group"
}

submit_data_auau() {
  echo
  echo "====================================================================="
  echo "Submitting AuAu numerator"
  echo "====================================================================="
  env \
    RJ_NOTIFY_EMAILS="$notify" \
    RJ_PROFILE_JOB=1 \
    RJ_CONFIG_YAML="$auau_yaml" \
    RJ_DEST_BASE_OVERRIDE="$auau_bulk" \
    RJ_MERGE_OUT_BASE_OVERRIDE="$merge_base" \
    RJ_REQUEST_MEMORY="$auau_memory" \
    RJ_AUTO_MERGE=1 \
    RJ_SIM_FIRSTROUND_REQUEST_MEMORY="$merge_memory" \
    RJ_SIM_MERGE_GROUP_SIZE="$merge_group" \
    RJ_ID_FANOUT_MAX_ROWS=1 \
    RJ_CODEX_CHAT_NAME="$RJ_CODEX_CHAT_NAME" \
    RJ_CODEX_THREAD_ID="$RJ_CODEX_THREAD_ID" \
    ./RecoilJets_Condor_submit.sh isAuAu condor all groupSize "$auau_group"
}

submit_embedded_signal() {
  echo
  echo "====================================================================="
  echo "Submitting embedded signal QA"
  echo "====================================================================="
  env \
    RJ_NOTIFY_EMAILS="$notify" \
    RJ_PROFILE_JOB=1 \
    RJ_CONFIG_YAML="$auau_yaml" \
    RJ_SIMEMBED_DEST_BASE="$signal_bulk" \
    RJ_MERGE_OUT_BASE_OVERRIDE="$merge_base" \
    RJ_REQUEST_MEMORY="$sim_memory" \
    RJ_AUTO_MERGE=1 \
    RJ_SIM_FIRSTROUND_REQUEST_MEMORY="$merge_memory" \
    RJ_SIM_MERGE_GROUP_SIZE="$merge_group" \
    RJ_ID_FANOUT_MAX_ROWS=1 \
    RJ_CODEX_CHAT_NAME="$RJ_CODEX_CHAT_NAME" \
    RJ_CODEX_THREAD_ID="$RJ_CODEX_THREAD_ID" \
    ./RecoilJets_Condor_submit.sh isSimEmbedded condorDoAll groupSize "$sim_group"
}

submit_embedded_background() {
  echo
  echo "====================================================================="
  echo "Submitting embedded background QA"
  echo "====================================================================="
  env \
    RJ_NOTIFY_EMAILS="$notify" \
    RJ_PROFILE_JOB=1 \
    RJ_CONFIG_YAML="$auau_yaml" \
    RJ_SIMEMBEDINCLUSIVE_DEST_BASE="$background_bulk" \
    RJ_MERGE_OUT_BASE_OVERRIDE="$merge_base" \
    RJ_REQUEST_MEMORY="$sim_memory" \
    RJ_AUTO_MERGE=1 \
    RJ_SIM_FIRSTROUND_REQUEST_MEMORY="$merge_memory" \
    RJ_SIM_MERGE_GROUP_SIZE="$merge_group" \
    RJ_ID_FANOUT_MAX_ROWS=1 \
    RJ_SIMEMBEDDEDINCLUSIVE_FOUR_SAMPLES=1 \
    RJ_CODEX_CHAT_NAME="$RJ_CODEX_CHAT_NAME" \
    RJ_CODEX_THREAD_ID="$RJ_CODEX_THREAD_ID" \
    ./RecoilJets_Condor_submit.sh isSimEmbeddedInclusive condorDoAll groupSize "$sim_group"
}

if [[ "$do_run" != "1" ]]; then
  cat <<EOF

DRY_INTENT_ONLY=1
Set RJ_DO_RUN=1 to execute the four submits. Set RJ_DAG_DRYRUN=1 with RJ_DO_RUN=1 to build/print DAGs without condor_submit_dag.
TRACK_THIS_CAMPAIGN=${campaign_tag}
SUBMIT_LOG=${submit_log}
PREFLIGHT_JSON=${preflight_out}
EOF
  exit 0
fi

echo "== pre-submit queue snapshot =="
condor_q -nobatch "${USER:-patsfan753}" 2>/dev/null | tail -50 || true

submit_data_pp
submit_data_auau
submit_embedded_signal
submit_embedded_background

echo
echo "== post-submit queue snapshot =="
condor_q -nobatch "${USER:-patsfan753}" 2>/dev/null | tail -80 || true
echo
echo "TRACK_THIS_CAMPAIGN=${campaign_tag}"
echo "SUBMIT_LOG=${submit_log}"
echo "PREFLIGHT_JSON=${preflight_out}"
echo "PP_BULK=${pp_bulk}"
echo "AUAU_BULK=${auau_bulk}"
echo "SIGNAL_BULK=${signal_bulk}"
echo "BACKGROUND_BULK=${background_bulk}"
echo "MERGE_BASE=${merge_base}"
