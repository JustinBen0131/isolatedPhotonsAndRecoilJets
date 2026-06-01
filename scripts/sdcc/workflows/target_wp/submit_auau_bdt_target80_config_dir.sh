#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

config_dir="${RJ_TARGET80_CONFIG_DIR:-${1:-}}"
[[ -n "$config_dir" ]] || { echo "[ERROR] Set RJ_TARGET80_CONFIG_DIR=/path/to/target80/configs or pass it as argv[1]" >&2; exit 2; }
[[ -d "$config_dir" ]] || { echo "[ERROR] Config directory not found: $config_dir" >&2; exit 2; }

campaign_group="${RJ_TARGET80_SUBMIT_GROUP:-target80_all_$(date +%Y%m%d_%H%M%S)}"
notify="${RJ_NOTIFY_EMAILS:-just0131@gmail.com}"
smoke="${RJ_TARGET80_RUN_LOCAL_SMOKE:-1}"
submit="${RJ_TARGET80_DO_SUBMIT:-0}"
pattern="${RJ_TARGET80_YAML_GLOB:-analysis_config_*_target80.yaml}"
queue_gate="${RJ_TARGETWP_QUEUE_GATE:-1}"
queue_scope="${RJ_TARGETWP_QUEUE_SCOPE:-matching}"
max_queued="${RJ_TARGETWP_MAX_QUEUED:-40000}"
resume_below="${RJ_TARGETWP_RESUME_BELOW:-0}"
poll_seconds="${RJ_TARGETWP_POLL_SECONDS:-300}"
log_dir="${RJ_TARGET80_LOG_DIR:-${repo_root}/condor_generated_configs/${campaign_group}}"
mkdir -p "$log_dir"
log_file="${log_dir}/submit_${campaign_group}.log"

exec > >(tee -a "$log_file") 2>&1

echo "RECOILJETS_AUAU_BDT_TARGET80_CONFIG_DIR_DRIVER_V1"
echo "submit_host=$(hostname -f 2>/dev/null || hostname)"
echo "repo_root=${repo_root}"
echo "config_dir=${config_dir}"
echo "campaign_group=${campaign_group}"
echo "yaml_glob=${pattern}"
echo "run_local_smoke=${smoke}"
echo "do_submit=${submit}"
echo "queue_gate=${queue_gate}"
echo "queue_scope=${queue_scope}"
echo "max_queued=${max_queued}"
echo "resume_below=${resume_below}"
echo "poll_seconds=${poll_seconds}"
echo "notify=${notify}"
echo "log_file=${log_file}"
echo

submitter="${repo_root}/RecoilJets_Condor_submit.sh"
if [[ ! -f "$submitter" ]]; then
  echo "[ERROR] Missing root-level submitter: $submitter" >&2
  echo "[ERROR] Run this from the RecoilJets SDCC checkout that contains RecoilJets_Condor_submit.sh." >&2
  exit 2
fi

mapfile -t yamls < <(find "$config_dir" -maxdepth 1 -type f -name "$pattern" | sort)
if (( ${#yamls[@]} == 0 )); then
  echo "[ERROR] No target80 YAMLs found under $config_dir matching $pattern" >&2
  exit 3
fi

echo "== target80 YAMLs =="
printf '  %s\n' "${yamls[@]}"
echo

for yaml in "${yamls[@]}"; do
  if ! grep -q "auau_tight_bdt_working_point_entries" "$yaml"; then
    echo "[ERROR] Missing auau_tight_bdt_working_point_entries in $yaml" >&2
    exit 4
  fi
done

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

if [[ "$smoke" == "1" ]]; then
  smoke_yaml=""
  for candidate in "${yamls[@]}"; do
    case "$(basename "$candidate")" in
      analysis_config_expanded_5to40_target80.yaml)
        smoke_yaml="$candidate"
        break
        ;;
    esac
  done
  [[ -n "$smoke_yaml" ]] || smoke_yaml="${yamls[0]}"
  echo "== local smoke =="
  echo "yaml=${smoke_yaml}"
  RJ_CONFIG_YAML="$smoke_yaml" \
    "$submitter" isSimEmbedded local 1000 NFILES=1 SAMPLE=run28_embeddedPhoton20 VERBOSE=1
  echo "[OK] local smoke completed for ${smoke_yaml}"
  echo
else
  echo "== local smoke skipped by RJ_TARGET80_RUN_LOCAL_SMOKE=${smoke} =="
fi

if [[ "$submit" != "1" ]]; then
  cat <<EOF
== DRY RUN COMPLETE ==
Set RJ_TARGET80_DO_SUBMIT=1 to submit all YAMLs sequentially.
Example:
  RJ_NOTIFY_EMAILS=${notify} \\
  RJ_TARGET80_CONFIG_DIR=${config_dir} \\
  RJ_TARGET80_SUBMIT_GROUP=${campaign_group} \\
  RJ_TARGET80_RUN_LOCAL_SMOKE=0 \\
  RJ_TARGET80_DO_SUBMIT=1 \\
  ./scripts/submit_auau_bdt_target80_config_dir.sh
EOF
  echo "TRACK_THIS_TARGET80_GROUP=${campaign_group}"
  echo "DRY_RUN_YAML_COUNT=${#yamls[@]}"
  exit 0
fi

echo "== preflight queue snapshot =="
condor_q -nobatch "${USER:-patsfan753}" 2>/dev/null | tail -40 || true
echo

for yaml in "${yamls[@]}"; do
  name="$(basename "$yaml" .yaml)"
  tag="${campaign_group}_${name}"
  wait_for_queue_gate "before ${name}" "$campaign_group"
  echo
  echo "====================================================================="
  echo "Submitting target80 pair"
  echo "  yaml : ${yaml}"
  echo "  tag  : ${tag}"
  echo "====================================================================="
  RJ_NOTIFY_EMAILS="$notify" \
  RJ_TARGETWP_CONFIG_YAML="$yaml" \
  RJ_TARGETWP_CAMPAIGN_TAG="$tag" \
  RJ_TARGETWP_QUEUE_GATE="$queue_gate" \
  RJ_TARGETWP_QUEUE_SCOPE="$queue_scope" \
  RJ_TARGETWP_MAX_QUEUED="$max_queued" \
  RJ_TARGETWP_RESUME_BELOW="$resume_below" \
  RJ_TARGETWP_POLL_SECONDS="$poll_seconds" \
    ./scripts/submit_auau_bdt_targetwp_pair.sh
  wait_for_queue_gate "after ${name}" "$campaign_group"
done

echo
echo "== post-submit queue snapshot =="
condor_q -nobatch "${USER:-patsfan753}" 2>/dev/null | tail -80 || true
echo
echo "TRACK_THIS_TARGET80_GROUP=${campaign_group}"
echo "SUBMITTED_YAML_COUNT=${#yamls[@]}"
echo "SUBMIT_LOG=${log_file}"
