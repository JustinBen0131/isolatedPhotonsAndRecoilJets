#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
if [[ -f "${script_dir}/../../runtime/io/recoiljets_io_paths.sh" ]]; then
  source "${script_dir}/../../runtime/io/recoiljets_io_paths.sh"
elif [[ -f scripts/sdcc/runtime/io/recoiljets_io_paths.sh ]]; then
  source scripts/sdcc/runtime/io/recoiljets_io_paths.sh
else
  echo "[THE105][ERROR] Cannot locate recoiljets_io_paths.sh" >&2
  exit 2
fi

repo_root="$(rj_find_repo_root "$script_dir")"
cd "$repo_root"

say() { printf '[THE105] %s\n' "$*"; }
die() { printf '[THE105][ERROR] %s\n' "$1" >&2; exit "${2:-2}"; }

mode="${1:-print}"
campaign_tag="${RJ_THE105_CAMPAIGN_TAG:-the105_auau_shower_contract_factorial_canary_20260717_2155}"
yaml="${RJ_THE105_CONFIG_YAML:-${repo_root}/macros/analysis_config_the88_bounded_sideband_default_bdt.yaml}"
build_root="${RJ_THE105_BUILD_ROOT:-${repo_root}/.recoiljets_tmp/the105_shower_contract_build_20260717}"
calo_build="${build_root}/caloreco_build"
calo_install="${build_root}/caloreco_install"
auau_build="${build_root}/auau_build"
auau_install="${build_root}/auau_install"
calo_library="${RJ_THE105_CALO_LIBRARY:-${calo_install}/lib/libcalo_reco.so}"
calo_header="${RJ_THE105_CALO_HEADER:-${calo_install}/include/caloreco/PhotonClusterBuilder.h}"
auau_library="${RJ_THE105_AUAU_LIBRARY:-${auau_install}/lib/libRecoilJetsAuAu.so}"
evidence_dir="${RJ_THE105_EVIDENCE_DIR:-${repo_root}/evidence/qa/the105_auau_shower_contract_factorial_20260717}"
canary_root="${RJ_THE105_CANARY_ROOT:-$(rj_recoiljets_bulk_root)/smoke/auau_shower_contract/${campaign_tag}}"

data_group="${RJ_THE105_DATA_GROUP_SIZE:-5}"
data_runs="${RJ_THE105_DATA_RUNS:-3}"
data_events="${RJ_THE105_DATA_EVENTS:-20000}"
sim_group="${RJ_THE105_SIM_GROUP_SIZE:-5}"
sim_events="${RJ_THE105_SIM_EVENTS:-20000}"
data_memory="${RJ_THE105_DATA_MEMORY:-6500MB}"
sim_memory="${RJ_THE105_SIM_MEMORY:-7500MB}"
allow_existing="${RJ_THE105_ALLOW_EXISTING:-0}"

variants=(historical towerinfo70 canonical)
signal_samples=(run28_embeddedPhoton12 run28_embeddedPhoton20)
inclusive_samples=(run28_embeddedJet12 run28_embeddedJet20 run28_embeddedJet30 run28_embeddedJet40)

export RJ_CODEX_CHAT_NAME="${RJ_CODEX_CHAT_NAME:-THE-105 | AuAu Shower-Contract Factorial Canary}"
export RJ_CODEX_THREAD_ID="${RJ_CODEX_THREAD_ID:-019f4c73-f704-78f0-9731-19e9cc9560c3}"

usage() {
  cat <<'EOF'
Usage: submit_the105_auau_shower_contract_factorial.sh MODE

Modes:
  print            Print the frozen three-arm diagnostic contract.
  build-isolated   Build campaign-local libcalo_reco and libRecoilJetsAuAu.
  canary-dryrun    Expand the exact bounded submission without submitting.
  canary-submit    Submit the three variants on matched bounded inputs.
  status           Report matching queue rows and output ROOT coverage.

This is a diagnostic canary only. It does not merge, promote, or replace any
canonical AuAu product. The tight-ID column evaluates the currently frozen
classifier and is not a retrained-model claim.
EOF
}

require_inputs() {
  rj_validate_campaign_tag "$campaign_tag" || exit 2
  [[ -s "$yaml" ]] || die "Missing analysis config: ${yaml}"
  [[ -s "$calo_library" ]] || die "Missing isolated CaloReco library: ${calo_library}. Run build-isolated first."
  [[ -s "$calo_header" ]] || die "Missing isolated PhotonClusterBuilder header: ${calo_header}. Run build-isolated first."
  [[ -s "$auau_library" ]] || die "Missing isolated AuAu library: ${auau_library}. Run build-isolated first."
  [[ -x ./RecoilJets_Condor_submit.sh ]] || die "Missing RecoilJets_Condor_submit.sh"
}

assert_fresh() {
  if [[ -e "$canary_root" && "$allow_existing" != "1" ]]; then
    die "Canary output root already exists: ${canary_root}. Use a fresh tag; set RJ_THE105_ALLOW_EXISTING=1 only for a proven resume."
  fi
  if condor_q "${USER:-patsfan753}" -af Args 2>/dev/null | grep -F "$canary_root" >/dev/null; then
    die "Active Condor rows already target ${canary_root}" 4
  fi
}

write_manifest() {
  mkdir -p "$evidence_dir"
  {
    printf 'variant\tpopulation\treconstruction_source\ttower_floor_gev\tselection_stages\n'
    printf 'historical\tdata\tTowerInfo full good grid\t0.070\tbefore,after_NCB_preselection,after_tight\n'
    printf 'historical\tembedded\tRawCluster towermap membership\t0.070\tbefore,after_NCB_preselection,after_tight\n'
    printf 'towerinfo70\tall\tTowerInfo full good grid\t0.070\tbefore,after_NCB_preselection,after_tight\n'
    printf 'canonical\tall\tTowerInfo full good grid\t0.000\tbefore,after_NCB_preselection,after_tight\n'
  } > "${evidence_dir}/campaign_contract.tsv"
  {
    printf 'campaign_tag=%s\n' "$campaign_tag"
    printf 'source_commit=%s\n' "$(git rev-parse HEAD 2>/dev/null || printf unknown)"
    printf 'config=%s\n' "$yaml"
    printf 'data_contract=first_%s_resolved_GRL_runs_groupSize_%s_up_to_%s_events_each\n' "$data_runs" "$data_group" "$data_events"
    printf 'sim_contract=first_group_of_%s_paired_rows_per_sample_up_to_%s_events\n' "$sim_group" "$sim_events"
    printf 'candidate_contract=15<=ET<35,abs_eta<0.7,one_skim_row_per_candidate\n'
    printf 'normalization_contract=full_finite_candidate_denominator_with_zero_underflow_overflow_reported_separately\n'
  } > "${evidence_dir}/campaign_manifest.txt"
}

print_contract() {
  write_manifest
  cat <<EOF
RECOILJETS_THE105_AUAU_SHOWER_CONTRACT_FACTORIAL_V1
campaign_tag=${campaign_tag}
config=${yaml}
variants=${variants[*]}
data=first ${data_runs} resolved GRL runs, groupSize ${data_group}, <=${data_events} events/job
signal_samples=${signal_samples[*]}
inclusive_samples=${inclusive_samples[*]}
sim=first ${sim_group} paired rows/sample, <=${sim_events} events/job
stages=before preselection; after complete NCB preselection; after frozen tight ID
candidate_skim=enabled, one row per candidate
maximum_jobs=$(( ${#variants[@]} * (data_runs + ${#signal_samples[@]} + ${#inclusive_samples[@]}) ))
automatic_merge=disabled
canonical_replacement=forbidden
canary_root=${canary_root}
isolated_calo_library=${calo_library}
isolated_auau_library=${auau_library}
EOF
}

reset_build_dir() {
  local path="$1"
  case "$path" in
    "${build_root}"/*) rm -rf "$path" ;;
    *) die "Refusing build cleanup outside ${build_root}: ${path}" ;;
  esac
  mkdir -p "$path"
}

build_component() {
  local source_dir="$1" build_dir="$2" install_dir="$3" label="$4"
  [[ -s "${source_dir}/autogen.sh" ]] || die "Missing ${label} autogen.sh: ${source_dir}/autogen.sh"
  reset_build_dir "$build_dir"
  reset_build_dir "$install_dir"
  cp -a "${source_dir}/." "$build_dir/"
  rm -rf "${build_dir}/autom4te.cache" "${build_dir}/.deps"
  rm -f "${build_dir}/config.status" "${build_dir}/config.log" \
        "${build_dir}/Makefile" "${build_dir}/libtool" "${build_dir}/stamp-h1"
  (
    cd "$build_dir"
    bash ./autogen.sh --prefix="$install_dir"
    make -j"${RJ_THE105_BUILD_JOBS:-4}"
    make install
  ) 2>&1 | tee "${evidence_dir}/${label}_build.log"
}

build_isolated() {
  case "$build_root" in
    "${repo_root}/.recoiljets_tmp/the105_shower_contract_build_"*) ;;
    *) die "Refusing isolated build outside the THE-105 hidden build root: ${build_root}" ;;
  esac
  mkdir -p "$build_root" "$evidence_dir"
  set +u
  source /opt/sphenix/core/bin/sphenix_setup.sh -n
  source /opt/sphenix/core/bin/setup_local.sh /sphenix/u/patsfan753/thesisAnalysis/install
  set -u
  build_component \
    "${repo_root}/coresoftware_local/offline/packages/CaloReco" \
    "$calo_build" "$calo_install" caloreco
  set +u
  source /opt/sphenix/core/bin/setup_local.sh "$calo_install"
  set -u
  build_component "${repo_root}/src_AuAu" "$auau_build" "$auau_install" recoiljets_auau
  [[ -s "$calo_library" ]] || die "CaloReco build did not produce ${calo_library}"
  [[ -s "$calo_header" ]] || die "CaloReco build did not install ${calo_header}"
  [[ -s "$auau_library" ]] || die "AuAu build did not produce ${auau_library}"
  sha256sum "$calo_library" "$calo_header" "$auau_library" | tee "${evidence_dir}/isolated_runtime.sha256"
}

common_env() {
  local variant="$1"
  printf '%s\n' \
    "RJ_CONFIG_YAML=${yaml}" \
    "RJ_CALO_RECO_LIBRARY_OVERRIDE=${calo_library}" \
    "RJ_PHOTON_CLUSTER_BUILDER_HEADER_OVERRIDE=${calo_header}" \
    "RJ_AUAU_LIBRARY_OVERRIDE=${auau_library}" \
    "RJ_SUBMIT_EXTRA_ENV=RJ_REQUIRE_EMBEDDED_MINBIAS_CLASSIFIER=1;RJ_AUAU_SHOWER_SHAPE_DIAGNOSTIC_VARIANT=${variant};RJ_AUAU_PHOTON_CANDIDATE_SKIM=1;RJ_AUAU_PHOTON_CANDIDATE_SKIM_MAX_ENTRIES=0" \
    "RJ_ID_FANOUT_MAX_ROWS=1" \
    "RJ_PHOTON_ID_ROW_MATCH=preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTSideband" \
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

submit_variant() {
  local variant="$1" dryrun="$2"
  local variant_root="${canary_root}/${variant}"
  local -a env_args
  mapfile -t env_args < <(common_env "$variant")

  say "Expanding ${variant}/data"
  env "${env_args[@]}" \
    "RJ_DAG_DRYRUN=${dryrun}" \
    "RJ_REQUEST_MEMORY=${data_memory}" \
    "RJ_SMOKE_OUTPUT_BASE=${variant_root}/data" \
    "RJ_SMOKE_DATA_RUNS=${data_runs}" \
    "RJ_SMOKE_DATA_NEVENTS=${data_events}" \
    ./RecoilJets_Condor_submit.sh isAuAu condor smokeTest groupSize "$data_group"

  local sample
  for sample in "${signal_samples[@]}"; do
    say "Expanding ${variant}/signal/${sample}"
    env "${env_args[@]}" \
      "RJ_DAG_DRYRUN=${dryrun}" \
      "RJ_REQUEST_MEMORY=${sim_memory}" \
      "RJ_SMOKE_OUTPUT_BASE=${variant_root}/signal" \
      "RJ_SMOKE_SIM_NEVENTS=${sim_events}" \
      ./RecoilJets_Condor_submit.sh isSimEmbedded condorDoAllSmoke \
        groupSize "$sim_group" maxJobs 1 "SAMPLE=${sample}"
  done

  for sample in "${inclusive_samples[@]}"; do
    say "Expanding ${variant}/inclusive/${sample}"
    env "${env_args[@]}" \
      "RJ_DAG_DRYRUN=${dryrun}" \
      "RJ_REQUEST_MEMORY=${sim_memory}" \
      "RJ_SMOKE_OUTPUT_BASE=${variant_root}/inclusive" \
      "RJ_SMOKE_SIM_NEVENTS=${sim_events}" \
      "RJ_SIMEMBEDDEDINCLUSIVE_FOUR_SAMPLES=1" \
      ./RecoilJets_Condor_submit.sh isSimEmbeddedInclusive condorDoAllSmoke \
        groupSize "$sim_group" maxJobs 1 "SAMPLE=${sample}"
  done
}

run_canary() {
  local dryrun="$1"
  require_inputs
  assert_fresh
  write_manifest
  mkdir -p "$evidence_dir"
  exec > >(tee -a "${evidence_dir}/canary_${dryrun}.log") 2>&1
  local variant
  for variant in "${variants[@]}"; do
    submit_variant "$variant" "$dryrun"
  done
}

status() {
  printf '%s\n' 'cluster proc status starts variant args'
  condor_q "${USER:-patsfan753}" -af ClusterId ProcId JobStatus NumJobStarts Args 2>/dev/null \
    | grep -F "$canary_root" || true
  printf '\nROOT coverage under %s\n' "$canary_root"
  find "$canary_root" -type f -name '*.root' -size +50000c 2>/dev/null \
    | awk -F/ '{print $(NF-4) "/" $(NF-3) "/" $(NF-2)}' \
    | sort | uniq -c || true
}

case "$mode" in
  print) print_contract ;;
  build-isolated) build_isolated ;;
  canary-dryrun) run_canary 1 ;;
  canary-submit) run_canary 0 ;;
  status) status ;;
  -h|--help|help) usage ;;
  *) usage >&2; die "Unknown mode: ${mode}" ;;
esac
