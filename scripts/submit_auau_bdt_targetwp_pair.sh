#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

yaml="${RJ_TARGETWP_CONFIG_YAML:-${1:-}}"
[[ -n "$yaml" ]] || { echo "[ERROR] Set RJ_TARGETWP_CONFIG_YAML=/path/to/generated_targetwp.yaml or pass it as argv[1]" >&2; exit 2; }
[[ -f "$yaml" ]] || { echo "[ERROR] Config YAML not found: $yaml" >&2; exit 2; }

base_name="$(basename "$yaml" .yaml)"
campaign_tag="${RJ_TARGETWP_CAMPAIGN_TAG:-${base_name}_$(date +%Y%m%d_%H%M%S)}"
group_size="${RJ_TARGETWP_GROUP_SIZE:-7}"
notify="${RJ_NOTIFY_EMAILS:-just0131@gmail.com}"
request_memory="${RJ_REQUEST_MEMORY:-12000MB}"
retry_cap="${RJ_AUTO_MEMORY_RETRY_CAP_MB:-16000}"
merge_memory="${RJ_SIM_FIRSTROUND_REQUEST_MEMORY:-8000MB}"
merge_group="${RJ_SIM_MERGE_GROUP_SIZE:-75}"
auto_merge="${RJ_TARGETWP_AUTO_MERGE:-0}"
queue_gate="${RJ_TARGETWP_QUEUE_GATE:-1}"
queue_scope="${RJ_TARGETWP_QUEUE_SCOPE:-matching}"
max_queued="${RJ_TARGETWP_MAX_QUEUED:-40000}"
resume_below="${RJ_TARGETWP_RESUME_BELOW:-0}"
poll_seconds="${RJ_TARGETWP_POLL_SECONDS:-300}"

log_dir="${RJ_TARGETWP_LOG_DIR:-${repo_root}/condor_generated_configs/${campaign_tag}}"
mkdir -p "$log_dir"
submit_log="${log_dir}/submit_${campaign_tag}.log"
exec > >(tee -a "$submit_log") 2>&1

echo "RECOILJETS_AUAU_BDT_TARGETWP_PAIR_SUBMIT_V1"
echo "submit_host=$(hostname -f 2>/dev/null || hostname)"
echo "repo_root=${repo_root}"
echo "campaign_tag=${campaign_tag}"
echo "yaml=${yaml}"
echo "group_size=${group_size}"
echo "request_memory=${request_memory}"
echo "retry_cap=${retry_cap}"
echo "merge_memory=${merge_memory}"
echo "merge_group_size=${merge_group}"
echo "auto_merge=${auto_merge}"
echo "queue_gate=${queue_gate}"
echo "queue_scope=${queue_scope}"
echo "max_queued=${max_queued}"
echo "resume_below=${resume_below}"
echo "poll_seconds=${poll_seconds}"
echo "notify=${notify}"
if [[ "$auto_merge" != "1" ]]; then
  echo "analysis_only=1"
  echo "note=Manual merge is intentional so failed workers remain inspectable."
fi
echo

if ! grep -q "auau_tight_bdt_working_point_entries" "$yaml"; then
  echo "[ERROR] YAML has no auau_tight_bdt_working_point_entries; refusing target-WP submission." >&2
  exit 3
fi

count_jobs() {
  local scope="$1"
  local pattern="$2"
  if [[ "$scope" == "user" ]]; then
    condor_q "${USER:-patsfan753}" -af JobStatus 2>/dev/null \
      | awk '$1 != 3 && $1 != 4 {n++} END {print n+0}'
  else
    condor_q "${USER:-patsfan753}" -af JobStatus Args 2>/dev/null \
      | awk -v pat="$pattern" '$1 != 3 && $1 != 4 && index($0, pat) {n++} END {print n+0}'
  fi
}

count_held_jobs() {
  local scope="$1"
  local pattern="$2"
  if [[ "$scope" == "user" ]]; then
    condor_q "${USER:-patsfan753}" -af JobStatus 2>/dev/null \
      | awk '$1 == 5 {n++} END {print n+0}'
  else
    condor_q "${USER:-patsfan753}" -af JobStatus Args 2>/dev/null \
      | awk -v pat="$pattern" '$1 == 5 && index($0, pat) {n++} END {print n+0}'
  fi
}

wait_for_queue_gate() {
  local label="$1"
  local pattern="$2"
  [[ "$queue_gate" == "1" ]] || return 0

  while true; do
    local queued held
    queued="$(count_jobs "$queue_scope" "$pattern")"
    held="$(count_held_jobs "$queue_scope" "$pattern")"
    echo "[queue_gate] ${label}: scope=${queue_scope} pattern=${pattern} queued=${queued} held=${held} resume_below=${resume_below} max_queued=${max_queued}"
    if (( held > 0 )); then
      echo "[ERROR] Queue gate sees held jobs for ${label}; stopping before submitting more." >&2
      condor_q -nobatch "${USER:-patsfan753}" 2>/dev/null | tail -80 || true
      exit 11
    fi
    if (( queued <= resume_below )); then
      break
    fi
    if (( queued > max_queued )); then
      echo "[WARN] ${label}: queued=${queued} is above max_queued=${max_queued}; waiting and not submitting more."
    fi
    sleep "$poll_seconds"
  done
}

submit_one() {
  local dataset="$1"
  local bulk_var="$2"
  local bulk_base="$3"
  local merge_base="$4"
  echo
  echo "====================================================================="
  echo "Submitting ${dataset}"
  echo "  campaign   : ${campaign_tag}"
  echo "  yaml       : ${yaml}"
  echo "  bulk base  : ${bulk_base}"
  echo "  merge base : ${merge_base}"
  echo "====================================================================="
  env \
    RJ_NOTIFY_EMAILS="${notify}" \
    RJ_PROFILE_JOB=1 \
    RJ_CONFIG_YAML="${yaml}" \
    RJ_ID_FANOUT_MAX_ROWS=1 \
    RJ_REQUEST_MEMORY="${request_memory}" \
    RJ_AUTO_MEMORY_RETRY_CAP_MB="${retry_cap}" \
    RJ_HOLD_FAILED_WORKERS=1 \
    RJ_AUTO_MERGE="${auto_merge}" \
    RJ_SIM_FIRSTROUND_REQUEST_MEMORY="${merge_memory}" \
    RJ_SIM_MERGE_GROUP_SIZE="${merge_group}" \
    RJ_MERGE_OUT_BASE_OVERRIDE="${merge_base}" \
    "${bulk_var}=${bulk_base}" \
    ./RecoilJets_Condor_submit.sh "${dataset}" condorDoAll groupSize "${group_size}"
}

echo "== preflight queue snapshot =="
condor_q -nobatch "${USER:-patsfan753}" 2>/dev/null | tail -30 || true

merge_base="/sphenix/u/patsfan753/scratch/thesisAnalysis/output_${campaign_tag}"
signal_bulk="/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded_${campaign_tag}"
background_bulk="/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembeddedinclusive_${campaign_tag}"

wait_for_queue_gate "before isSimEmbedded" "$campaign_tag"
submit_one isSimEmbedded RJ_SIMEMBED_DEST_BASE "$signal_bulk" "$merge_base"
wait_for_queue_gate "before isSimEmbeddedInclusive" "$campaign_tag"
submit_one isSimEmbeddedInclusive RJ_SIMEMBEDINCLUSIVE_DEST_BASE "$background_bulk" "$merge_base"
wait_for_queue_gate "after isSimEmbeddedInclusive" "$campaign_tag"

if [[ "$auto_merge" != "1" ]]; then
  cat <<EOF

== merge commands after all matching worker jobs are clean ==
env MERGE_CONFIG_YAML=${yaml} MERGE_RUN_BASE_OVERRIDE=${signal_bulk} MERGE_OUT_BASE_OVERRIDE=${merge_base} RJ_SIM_FIRSTROUND_REQUEST_MEMORY=${merge_memory} ./scripts/mergeRecoilJets.sh isSimEmbedded firstRound groupSize ${merge_group}
env MERGE_CONFIG_YAML=${yaml} MERGE_RUN_BASE_OVERRIDE=${signal_bulk} MERGE_OUT_BASE_OVERRIDE=${merge_base} RJ_SIM_FIRSTROUND_REQUEST_MEMORY=${merge_memory} ./scripts/mergeRecoilJets.sh isSimEmbedded secondRound groupSize ${merge_group}
env MERGE_CONFIG_YAML=${yaml} MERGE_RUN_BASE_OVERRIDE=${signal_bulk} MERGE_OUT_BASE_OVERRIDE=${merge_base} RJ_SIM_FIRSTROUND_REQUEST_MEMORY=${merge_memory} ./scripts/mergeRecoilJets.sh isSimEmbedded finalStitch groupSize ${merge_group}
env MERGE_CONFIG_YAML=${yaml} MERGE_RUN_BASE_OVERRIDE=${background_bulk} MERGE_OUT_BASE_OVERRIDE=${merge_base} RJ_SIM_FIRSTROUND_REQUEST_MEMORY=${merge_memory} ./scripts/mergeRecoilJets.sh isSimEmbeddedInclusive firstRound groupSize ${merge_group}
env MERGE_CONFIG_YAML=${yaml} MERGE_RUN_BASE_OVERRIDE=${background_bulk} MERGE_OUT_BASE_OVERRIDE=${merge_base} RJ_SIM_FIRSTROUND_REQUEST_MEMORY=${merge_memory} ./scripts/mergeRecoilJets.sh isSimEmbeddedInclusive secondRound groupSize ${merge_group}
env MERGE_CONFIG_YAML=${yaml} MERGE_RUN_BASE_OVERRIDE=${background_bulk} MERGE_OUT_BASE_OVERRIDE=${merge_base} RJ_SIM_FIRSTROUND_REQUEST_MEMORY=${merge_memory} ./scripts/mergeRecoilJets.sh isSimEmbeddedInclusive finalStitch groupSize ${merge_group}
EOF
fi

echo
echo "== post-submit queue snapshot =="
condor_q -nobatch "${USER:-patsfan753}" 2>/dev/null | tail -60 || true
echo
echo "TRACK_THIS_CAMPAIGN=${campaign_tag}"
echo "GENERATED_YAML=${yaml}"
echo "SUBMIT_LOG=${submit_log}"
