#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
if [[ -f "${script_dir}/../../runtime/io/recoiljets_io_paths.sh" ]]; then
  source "${script_dir}/../../runtime/io/recoiljets_io_paths.sh"
elif [[ -f scripts/sdcc/runtime/io/recoiljets_io_paths.sh ]]; then
  source scripts/sdcc/runtime/io/recoiljets_io_paths.sh
else
  echo "[THE100][ERROR] Cannot locate recoiljets_io_paths.sh" >&2
  exit 2
fi

repo_root="$(rj_find_repo_root "$script_dir")"
cd "$repo_root"

say() { printf '[THE100] %s\n' "$*"; }
die() { printf '[THE100][ERROR] %s\n' "$1" >&2; exit "${2:-2}"; }

mode="${1:-print}"
campaign_tag="${RJ_THE100_CAMPAIGN_TAG:-the100_auau_dualview_20260714}"
family="${RJ_THE100_FAMILY:-auau_dualview}"
yaml="${RJ_THE100_CONFIG_YAML:-${repo_root}/macros/analysis_config_the100_auau_dualview_comparison.yaml}"
build_root="${RJ_THE100_BUILD_ROOT:-${repo_root}/.recoiljets_tmp/the100_dualview_build_20260714}"
build_dir="${build_root}/build"
install_dir="${build_root}/install"
auau_library="${RJ_THE100_AUAU_LIBRARY:-${install_dir}/lib/libRecoilJetsAuAu.so}"
evidence_dir="${RJ_THE100_EVIDENCE_DIR:-${repo_root}/evidence/qa/the100_dualview_20260714}"
manifest="${evidence_dir}/campaign_manifest.tsv"
submit_log="${evidence_dir}/submission.log"

bulk_root="${RJ_THE100_BULK_ROOT:-$(rj_recoiljets_bulk_root)/${family}/${campaign_tag}}"
data_bulk="${bulk_root}/data"
signal_bulk="${bulk_root}/signal"
inclusive_bulk="${bulk_root}/inclusive"
canary_root="${RJ_THE100_CANARY_ROOT:-$(rj_recoiljets_bulk_root)/smoke/${family}/${campaign_tag}}"
canary_data="${canary_root}/data"
canary_signal="${canary_root}/signal"
canary_inclusive="${canary_root}/inclusive"
merge_base="${RJ_THE100_MERGE_BASE:-$(rj_merge_root_for_campaign "$campaign_tag")}"

data_group="${RJ_THE100_DATA_GROUP_SIZE:-7}"
sim_group="${RJ_THE100_SIM_GROUP_SIZE:-1}"
data_memory="${RJ_THE100_DATA_MEMORY:-5000MB}"
sim_memory="${RJ_THE100_SIM_MEMORY:-6500MB}"
canary_events="${RJ_THE100_CANARY_EVENTS:-10000}"
allow_existing="${RJ_THE100_ALLOW_EXISTING:-0}"

export RJ_CODEX_CHAT_NAME="${RJ_CODEX_CHAT_NAME:-THE-100 | AuAu Dual-View ABCD Comparison}"
export RJ_CODEX_THREAD_ID="${RJ_CODEX_THREAD_ID:-019f4346-8021-7e73-af6d-16ddfb1e438d}"

usage() {
  cat <<'EOF'
Usage: submit_the100_auau_dualview_campaign.sh MODE

Modes:
  print                Print the immutable campaign contract.
  build-isolated       Build libRecoilJetsAuAu.so in the campaign-only prefix.
  canary-dryrun        Disabled: the generic smoke dry-run path submits jobs.
  canary-submit        Submit data, Photon12/20, and Jet12/20/30/40 canaries.
  production-dryrun    Print exact full-production job counts without submitting.
  production-submit    Submit full data, signal, and four-sample inclusive lanes.
  production-submit-inclusive-sample SAMPLE
                       Resume one canonical inclusive sample after a queue-cap gate.
  status               Report matching queue rows and output coverage.

Production submission is analysis-only. Canonical merge is separately guarded.
EOF
}

require_inputs() {
  rj_validate_campaign_tag "$campaign_tag" || exit 2
  [[ -s "$yaml" ]] || die "Missing THE-100 config: ${yaml}"
  [[ -s "$auau_library" ]] || die "Missing isolated AuAu library: ${auau_library}. Run build-isolated first."
  [[ -x ./RecoilJets_Condor_submit.sh ]] || die "Missing executable RecoilJets_Condor_submit.sh in ${repo_root}"
}

assert_fresh_path() {
  local path="$1" label="$2"
  if [[ -e "$path" && "$allow_existing" != "1" ]]; then
    die "${label} already exists: ${path}. Use a fresh campaign tag; set RJ_THE100_ALLOW_EXISTING=1 only for an intentional resume."
  fi
}

assert_no_active_duplicate() {
  local output_root="${1:-$bulk_root}"
  if condor_q "${USER:-patsfan753}" -af Args 2>/dev/null | grep -F "$output_root" >/dev/null; then
    condor_q "${USER:-patsfan753}" -af ClusterId ProcId JobStatus Args 2>/dev/null | grep -F "$output_root" | head -20 || true
    die "Active Condor jobs already target ${output_root}; refusing duplicate submission." 4
  fi
}

write_manifest() {
  mkdir -p "$evidence_dir"
  {
    printf 'lane\tdataset\tsamples\toutput_base\n'
    printf 'data\tisAuAu\tTanner_GRL_2776_paired_DST_JET_DST_JETCALO\t%s\n' "$data_bulk"
    printf 'signal\tisSimEmbedded\trun28_embeddedPhoton12+run28_embeddedPhoton20\t%s\n' "$signal_bulk"
    printf 'inclusive\tisSimEmbeddedInclusive\trun28_embeddedJet12+run28_embeddedJet20+run28_embeddedJet30+run28_embeddedJet40\t%s\n' "$inclusive_bulk"
  } > "$manifest"
}

print_contract() {
  write_manifest
  cat <<EOF
RECOILJETS_THE100_AUAU_DUALVIEW_V1
campaign_tag=${campaign_tag}
config=${yaml}
isolated_library=${auau_library}
data_input=Tanner 2776-run GRL, paired DST_JET + DST_JETCALO
signal_samples=run28_embeddedPhoton12,run28_embeddedPhoton20
inclusive_samples=run28_embeddedJet12,run28_embeddedJet20,run28_embeddedJet30,run28_embeddedJet40
photon_id_rows=auauBDTSideband,auauBDTComplement
bounded_definition=T80(centrality)-0.20 < score < T80(centrality)-0.03
unrestricted_definition=score < T80(centrality)
diagnostic_surface=score-T80(centrality) versus Eiso versus photon-pT by centrality
data_bulk=${data_bulk}
signal_bulk=${signal_bulk}
inclusive_bulk=${inclusive_bulk}
canary_root=${canary_root}
merge_base=${merge_base}
automatic_merge=disabled
RJ_CODEX_CHAT_NAME=${RJ_CODEX_CHAT_NAME}
RJ_CODEX_THREAD_ID=${RJ_CODEX_THREAD_ID}
manifest=${manifest}
EOF
}

build_isolated() {
  [[ -s "${repo_root}/src_AuAu/autogen.sh" ]] || die "Missing src_AuAu/autogen.sh"
  case "$build_root" in
    "${repo_root}/.recoiljets_tmp/the100_dualview_build_"*) ;;
    *) die "Refusing to rebuild outside the THE-100-owned hidden path: ${build_root}" ;;
  esac
  rm -rf "$build_dir" "$install_dir"
  mkdir -p "$build_dir" "$install_dir" "$evidence_dir"
  cp -a "${repo_root}/src_AuAu/." "$build_dir/"
  rm -rf "${build_dir}/autom4te.cache" "${build_dir}/.deps"
  rm -f "${build_dir}/config.status" "${build_dir}/config.log" \
        "${build_dir}/Makefile" "${build_dir}/libtool" "${build_dir}/stamp-h1"
  set +u
  source /opt/sphenix/core/bin/sphenix_setup.sh -n
  source /opt/sphenix/core/bin/setup_local.sh "${RJ_THE100_DEP_INSTALL:-/sphenix/u/patsfan753/thesisAnalysis/install}"
  set -u
  (
    cd "$build_dir"
    bash ./autogen.sh --prefix="$install_dir"
    make -j"${RJ_THE100_BUILD_JOBS:-4}"
    make install
  ) 2>&1 | tee "${evidence_dir}/isolated_build.log"
  [[ -s "$auau_library" ]] || die "Build finished without ${auau_library}"
  sha256sum "$auau_library" | tee "${evidence_dir}/isolated_library.sha256"
}

common_env() {
  printf '%s\n' \
    "RJ_CONFIG_YAML=${yaml}" \
    "RJ_AUAU_LIBRARY_OVERRIDE=${auau_library}" \
    "RJ_SUBMIT_EXTRA_ENV=RJ_AUAU_DUALVIEW_DIAGNOSTICS=1" \
    "RJ_ID_FANOUT_MAX_ROWS=2" \
    "RJ_AUTO_MERGE=0" \
    "RJ_REQUIRE_NON_TINY_OUTPUT=1" \
    "RJ_MIN_OUTPUT_BYTES=50000" \
    "RJ_PROFILE_JOB=1" \
    "RJ_JOB_HEARTBEAT_SECONDS=120" \
    "RJ_INTERNAL_FIXED_ISO_GEV_AUAU=4.0" \
    "RJ_AUAU_BUILD_TOPOCLUSTER_ISOLATION=0" \
    "RJ_AUAU_USE_TOPOCLUSTER_ISOLATION=0" \
    "RJ_CODEX_CHAT_NAME=${RJ_CODEX_CHAT_NAME}" \
    "RJ_CODEX_THREAD_ID=${RJ_CODEX_THREAD_ID}"
}

submit_canary_lane() {
  local dryrun="$1" dataset="$2" sample="$3" output_base="$4" memory="$5"
  local -a env_args extra_args
  mapfile -t env_args < <(common_env)
  extra_args=()
  if [[ "$dataset" == "isSimEmbeddedInclusive" ]]; then
    extra_args+=("RJ_SIMEMBEDDEDINCLUSIVE_FOUR_SAMPLES=1")
  fi
  env "${env_args[@]}" "${extra_args[@]}" \
    "RJ_DAG_DRYRUN=${dryrun}" \
    "RJ_REQUEST_MEMORY=${memory}" \
    "RJ_SMOKE_OUTPUT_BASE=${output_base}" \
    "RJ_SMOKE_SIM_NEVENTS=${canary_events}" \
    ./RecoilJets_Condor_submit.sh "$dataset" condorDoAllSmoke groupSize 1 maxJobs 1 "SAMPLE=${sample}"
}

run_canaries() {
  local dryrun="$1"
  require_inputs
  if [[ "$dryrun" == "1" ]]; then
    canary_root="${canary_root}_dryrun"
    canary_data="${canary_root}/data"
    canary_signal="${canary_root}/signal"
    canary_inclusive="${canary_root}/inclusive"
  fi
  assert_no_active_duplicate "$canary_root"
  if [[ "$dryrun" == "0" ]]; then
    assert_fresh_path "$canary_root" "canary output root"
  fi
  mkdir -p "$evidence_dir"
  exec > >(tee -a "${evidence_dir}/canary_submission.log") 2>&1
  local -a env_args
  mapfile -t env_args < <(common_env)
  env "${env_args[@]}" \
    "RJ_DAG_DRYRUN=${dryrun}" \
    "RJ_REQUEST_MEMORY=${data_memory}" \
    "RJ_SMOKE_OUTPUT_BASE=${canary_data}" \
    "RJ_SMOKE_DATA_RUNS=1" \
    "RJ_SMOKE_DATA_NEVENTS=${canary_events}" \
    ./RecoilJets_Condor_submit.sh isAuAu condor smokeTest groupSize 1

  submit_canary_lane "$dryrun" isSimEmbedded run28_embeddedPhoton12 "$canary_signal" "$sim_memory"
  submit_canary_lane "$dryrun" isSimEmbedded run28_embeddedPhoton20 "$canary_signal" "$sim_memory"
  submit_canary_lane "$dryrun" isSimEmbeddedInclusive run28_embeddedJet12 "$canary_inclusive" "$sim_memory"
  submit_canary_lane "$dryrun" isSimEmbeddedInclusive run28_embeddedJet20 "$canary_inclusive" "$sim_memory"
  submit_canary_lane "$dryrun" isSimEmbeddedInclusive run28_embeddedJet30 "$canary_inclusive" "$sim_memory"
  submit_canary_lane "$dryrun" isSimEmbeddedInclusive run28_embeddedJet40 "$canary_inclusive" "$sim_memory"
}

production_counts() {
  require_inputs
  local -a env_args
  mapfile -t env_args < <(common_env)
  say "Full data count"
  env "${env_args[@]}" "RJ_REQUEST_MEMORY=${data_memory}" \
    ./RecoilJets_Condor_submit.sh isAuAu CHECKJOBS groupSize "$data_group"
  say "Full Photon12+20 count"
  env "${env_args[@]}" "RJ_REQUEST_MEMORY=${sim_memory}" \
    ./RecoilJets_Condor_submit.sh isSimEmbedded CHECKJOBS groupSize "$sim_group"
  say "Full Jet12+20+30+40 count"
  env "${env_args[@]}" "RJ_REQUEST_MEMORY=${sim_memory}" \
    "RJ_SIMEMBEDDEDINCLUSIVE_FOUR_SAMPLES=1" \
    ./RecoilJets_Condor_submit.sh isSimEmbeddedInclusive CHECKJOBS groupSize "$sim_group"
}

submit_production() {
  require_inputs
  assert_no_active_duplicate "$bulk_root"
  assert_fresh_path "$bulk_root" "production bulk root"
  assert_fresh_path "$merge_base" "production merge root"
  write_manifest
  mkdir -p "$evidence_dir"
  exec > >(tee -a "$submit_log") 2>&1
  local -a env_args
  mapfile -t env_args < <(common_env)

  say "Submitting full Tanner-GRL paired AuAu data lane"
  env "${env_args[@]}" \
    "RJ_REQUEST_MEMORY=${data_memory}" \
    "RJ_DEST_BASE_OVERRIDE=${data_bulk}" \
    "RJ_MERGE_OUT_BASE_OVERRIDE=${merge_base}/data" \
    ./RecoilJets_Condor_submit.sh isAuAu condor all groupSize "$data_group"

  say "Submitting full embedded Photon12+20 signal lane"
  env "${env_args[@]}" \
    "RJ_REQUEST_MEMORY=${sim_memory}" \
    "RJ_SIMEMBED_DEST_BASE=${signal_bulk}" \
    "RJ_MERGE_OUT_BASE_OVERRIDE=${merge_base}/signal" \
    ./RecoilJets_Condor_submit.sh isSimEmbedded condorDoAll groupSize "$sim_group"

  say "Submitting full embedded Jet12+20+30+40 inclusive lane"
  env "${env_args[@]}" \
    "RJ_REQUEST_MEMORY=${sim_memory}" \
    "RJ_SIMEMBEDINCLUSIVE_DEST_BASE=${inclusive_bulk}" \
    "RJ_MERGE_OUT_BASE_OVERRIDE=${merge_base}/inclusive" \
    "RJ_SIMEMBEDDEDINCLUSIVE_FOUR_SAMPLES=1" \
    ./RecoilJets_Condor_submit.sh isSimEmbeddedInclusive condorDoAll groupSize "$sim_group"

  condor_q "${USER:-patsfan753}" -af ClusterId JobStatus Args 2>/dev/null | grep -F "$campaign_tag" || true
}

submit_inclusive_sample() {
  local sample="${1:-}"
  case "$sample" in
    run28_embeddedJet12|run28_embeddedJet20|run28_embeddedJet30|run28_embeddedJet40) ;;
    *) die "Inclusive resume requires one canonical sample: run28_embeddedJet12, run28_embeddedJet20, run28_embeddedJet30, or run28_embeddedJet40." ;;
  esac
  require_inputs
  if condor_q "${USER:-patsfan753}" -af Args 2>/dev/null | \
      awk -v root="$inclusive_bulk" -v sample="$sample" \
        'index($0, root) && index($0, sample) { found=1 } END { exit(found ? 0 : 1) }'; then
    die "Active Condor jobs already target ${inclusive_bulk} for ${sample}; refusing duplicate submission." 4
  fi
  mkdir -p "$evidence_dir"
  exec > >(tee -a "${evidence_dir}/inclusive_resume_${sample}.log") 2>&1
  local -a env_args
  mapfile -t env_args < <(common_env)
  say "Submitting queue-cap recovery for ${sample}"
  env "${env_args[@]}" \
    "RJ_REQUEST_MEMORY=${sim_memory}" \
    "RJ_SIMEMBEDINCLUSIVE_DEST_BASE=${inclusive_bulk}" \
    "RJ_MERGE_OUT_BASE_OVERRIDE=${merge_base}/inclusive" \
    "RJ_SIMEMBEDDEDINCLUSIVE_FOUR_SAMPLES=1" \
    ./RecoilJets_Condor_submit.sh isSimEmbeddedInclusive condorDoAll groupSize "$sim_group" "SAMPLE=${sample}"
}

status_report() {
  print_contract
  printf '\nQUEUE\n'
  condor_q "${USER:-patsfan753}" -af ClusterId ProcId JobStatus HoldReason Args 2>/dev/null | grep -F "$campaign_tag" || true
  printf '\nROOT_COUNTS\n'
  for path in "$canary_data" "$canary_signal" "$canary_inclusive" "$data_bulk" "$signal_bulk" "$inclusive_bulk" "$merge_base"; do
    count="$(find "$path" -type f -name '*.root' -size +50k 2>/dev/null | wc -l | tr -d ' ')"
    printf '%s\t%s\n' "$count" "$path"
  done
}

case "$mode" in
  print) print_contract ;;
  build-isolated) build_isolated ;;
  canary-dryrun) die "Disabled: the generic smoke dry-run path still submits jobs. Use production-dryrun for read-only counts and canary-submit only for an intentional fresh canary tag." 63 ;;
  canary-submit) run_canaries 0 ;;
  production-dryrun) production_counts ;;
  production-submit) submit_production ;;
  production-submit-inclusive-sample) submit_inclusive_sample "${2:-}" ;;
  status) status_report ;;
  help|-h|--help) usage ;;
  *) usage >&2; die "Unknown mode: ${mode}" ;;
esac
