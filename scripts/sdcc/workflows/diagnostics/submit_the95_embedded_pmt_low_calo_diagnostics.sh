#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
if [[ -f "${script_dir}/../../runtime/io/recoiljets_io_paths.sh" ]]; then
  source "${script_dir}/../../runtime/io/recoiljets_io_paths.sh"
elif [[ -f scripts/sdcc/runtime/io/recoiljets_io_paths.sh ]]; then
  source scripts/sdcc/runtime/io/recoiljets_io_paths.sh
else
  echo "[ERROR] cannot locate recoiljets_io_paths.sh" >&2
  exit 2
fi

repo_root="$(rj_find_repo_root "$script_dir")"
cd "$repo_root"

mode="${1:-print}"
campaign_tag="${RJ_THE95_PMT_CAMPAIGN_TAG:-the95_embedded_pmt_low_calo_$(date +%Y%m%d_%H%M)}"
family="${RJ_THE95_PMT_FAMILY:-the95_pmt_low_calo}"
yaml="${RJ_THE95_PMT_YAML:-${repo_root}/macros/analysis_config_the95_embedded_pmt_low_calo_diagnostics.yaml}"
group_size="${RJ_THE95_PMT_GROUP_SIZE:-10}"
signal_memory="${RJ_THE95_PMT_SIGNAL_MEMORY:-5000MB}"
inclusive_memory="${RJ_THE95_PMT_INCLUSIVE_MEMORY:-6000MB}"
do_run="${RJ_DO_RUN:-0}"
allow_existing="${RJ_THE95_PMT_ALLOW_EXISTING:-0}"
require_embedded_minbias="${RJ_THE95_REQUIRE_EMBEDDED_MINBIAS:-0}"

rj_validate_campaign_tag "$campaign_tag" || exit 2
rj_validate_campaign_tag "$family" || exit 2

signal_base="${RJ_THE95_PMT_SIGNAL_BASE:-$(rj_bulk_signal_root_for_campaign "$family" "$campaign_tag")}"
inclusive_base="${RJ_THE95_PMT_INCLUSIVE_BASE:-$(rj_bulk_background_root_for_campaign "$family" "$campaign_tag")}"
merge_base="${RJ_THE95_PMT_MERGE_BASE:-$(rj_merge_root_for_campaign "$campaign_tag")}"
state_dir="${RJ_THE95_PMT_STATE_DIR:-${repo_root}/condor_generated_configs/${campaign_tag}}"
manifest="${state_dir}/campaign_manifest.tsv"

export RJ_CODEX_CHAT_NAME="${RJ_CODEX_CHAT_NAME:-THE-95_PMT_Low-Calo_Diagnostic}"
export RJ_CODEX_THREAD_ID="${RJ_CODEX_THREAD_ID:-${CODEX_THREAD_ID:-THE-95-local-codex-thread}}"

die() { printf '[THE95-PMT][ERROR] %s\n' "$*" >&2; exit 2; }
say() { printf '[THE95-PMT] %s\n' "$*"; }

write_manifest() {
  mkdir -p "$state_dir"
  cat > "$manifest" <<EOF
dataset	samples	output_base	yaml
isSimEmbedded	run28_embeddedPhoton12,run28_embeddedPhoton20	${signal_base}	${yaml}
isSimEmbeddedInclusive	run28_embeddedJet12,run28_embeddedJet20,run28_embeddedJet30,run28_embeddedJet40	${inclusive_base}	${yaml}
EOF
}

preflight() {
  [[ -f "$yaml" ]] || die "missing YAML: $yaml"
  grep -Fq 'analysis_config_the95_embedded_pmt_low_calo_diagnostics' "$yaml" || true
  if [[ "$allow_existing" != "1" ]]; then
    [[ ! -e "$signal_base" ]] || die "signal output exists: $signal_base"
    [[ ! -e "$inclusive_base" ]] || die "inclusive output exists: $inclusive_base"
    [[ ! -e "$merge_base" ]] || die "merge output exists: $merge_base"
  fi
  if condor_q "${USER:-patsfan753}" -af ClusterId Args 2>/dev/null | grep -F "$campaign_tag" >/dev/null; then
    condor_q "${USER:-patsfan753}" -af ClusterId JobStatus Args 2>/dev/null | grep -F "$campaign_tag" || true
    die "active Condor rows already mention $campaign_tag"
  fi
}

common_env() {
  printf '%s\n' \
    "RJ_CONFIG_YAML=${yaml}" \
    "RJ_MBD_PMT_LOW_CALO_DIAGNOSTICS=1" \
    "RJ_REQUIRE_EMBEDDED_MINBIAS_CLASSIFIER=${require_embedded_minbias}" \
    "RJ_FORCE_CALO_TOWER_STATUS_FOR_EMBEDDED=1" \
    "RJ_EVENT_CALO_REQUIRE_ISGOOD=1" \
    "RJ_DISABLE_ID_FANOUT=1" \
    "RJ_ID_FANOUT_MAX_ROWS=1" \
    "RJ_PHOTON_ID_ROW_MATCH=preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTSideband" \
    "RJ_DISABLE_ISO_CONE_INTERNALIZATION=1" \
    "RJ_DISABLE_JET_PT_INTERNALIZATION=1" \
    "RJ_DISABLE_DPHI_INTERNALIZATION=1" \
    "RJ_AUTO_MERGE=0" \
    "RJ_PROFILE_JOB=1" \
    "RJ_CODEX_CHAT_NAME=${RJ_CODEX_CHAT_NAME}" \
    "RJ_CODEX_THREAD_ID=${RJ_CODEX_THREAD_ID}"
}

run_dataset() {
  local dataset="$1"
  local output_base="$2"
  local memory="$3"
  local sample="${4:-}"
  local submit_mode="condorDoAll"
  local -a args=(groupSize "$group_size")
  if [[ "$mode" == "smoke" ]]; then
    submit_mode="condorDoAllSmoke"
    args=(groupSize 1 maxJobs 1 "SAMPLE=${sample}")
  fi
  mapfile -t env_args < <(common_env)
  env_args+=("RJ_REQUEST_MEMORY=${memory}")
  if [[ "$dataset" == "isSimEmbedded" ]]; then
    env_args+=("RJ_SIMEMBED_DEST_BASE=${output_base}")
  else
    env_args+=("RJ_SIMEMBEDINCLUSIVE_DEST_BASE=${output_base}" "RJ_SIMEMBEDDEDINCLUSIVE_FOUR_SAMPLES=1")
  fi
  say "submit ${dataset} output=${output_base} mode=${submit_mode} sample=${sample:-all}"
  env "${env_args[@]}" ./RecoilJets_Condor_submit.sh "$dataset" "$submit_mode" "${args[@]}"
}

print_contract() {
  cat <<EOF
RECOILJETS_THE95_EMBEDDED_PMT_LOW_CALO_V1
campaign_tag=${campaign_tag}
family=${family}
yaml=${yaml}
signal_base=${signal_base}
inclusive_base=${inclusive_base}
merge_base=${merge_base}
manifest=${manifest}
group_size=${group_size}
signal_memory=${signal_memory}
inclusive_memory=${inclusive_memory}
samples_signal=run28_embeddedPhoton12,run28_embeddedPhoton20
samples_inclusive=run28_embeddedJet12,run28_embeddedJet20,run28_embeddedJet30,run28_embeddedJet40
automatic_merge=disabled
diagnostic_only=1
require_embedded_minbias_classifier=${require_embedded_minbias}
id_fanout=disabled
RJ_DO_RUN=${do_run}
EOF
}

case "$mode" in
  print)
    write_manifest
    preflight
    print_contract
    ;;
  smoke)
    write_manifest
    preflight
    print_contract
    [[ "$do_run" == "1" ]] || die "smoke requires RJ_DO_RUN=1"
    run_dataset isSimEmbedded "${signal_base}_smoke_photon12" "$signal_memory" run28_embeddedPhoton12
    run_dataset isSimEmbedded "${signal_base}_smoke_photon20" "$signal_memory" run28_embeddedPhoton20
    run_dataset isSimEmbeddedInclusive "${inclusive_base}_smoke_jet12" "$inclusive_memory" run28_embeddedJet12
    run_dataset isSimEmbeddedInclusive "${inclusive_base}_smoke_jet40" "$inclusive_memory" run28_embeddedJet40
    ;;
  submit)
    write_manifest
    preflight
    print_contract
    [[ "$do_run" == "1" ]] || die "submit requires RJ_DO_RUN=1"
    run_dataset isSimEmbedded "$signal_base" "$signal_memory"
    run_dataset isSimEmbeddedInclusive "$inclusive_base" "$inclusive_memory"
    ;;
  *)
    die "usage: $0 {print|smoke|submit}"
    ;;
esac
