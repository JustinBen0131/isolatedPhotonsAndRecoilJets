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
photon_builder_build="${build_root}/photon_builder_build"
photon_builder_install="${build_root}/photon_builder_install"
auau_build="${build_root}/auau_build"
auau_install="${build_root}/auau_install"
calo_library="${RJ_THE105_CALO_LIBRARY:-/sphenix/u/${USER:-patsfan753}/thesisAnalysis/install/lib/libcalo_reco.so}"
photon_builder_library="${RJ_THE105_PHOTON_BUILDER_LIBRARY:-${photon_builder_install}/lib/libphoton_cluster_builder_override.so}"
calo_header="${RJ_THE105_CALO_HEADER:-${photon_builder_install}/include/caloreco/PhotonClusterBuilder.h}"
auau_library="${RJ_THE105_AUAU_LIBRARY:-${auau_install}/lib/libRecoilJetsAuAu.so}"
evidence_dir="${RJ_THE105_EVIDENCE_DIR:-${repo_root}/evidence/qa/the105_auau_shower_contract_factorial_20260717}"
canary_root="${RJ_THE105_CANARY_ROOT:-$(rj_recoiljets_bulk_root)/smoke/auau_shower_contract/${campaign_tag}}"

data_group="${RJ_THE105_DATA_GROUP_SIZE:-5}"
data_runs="${RJ_THE105_DATA_RUNS:-3}"
data_max_jobs="${RJ_THE105_DATA_MAX_JOBS:-${data_runs}}"
data_events="${RJ_THE105_DATA_EVENTS:-20000}"
sim_sample_rows="${RJ_THE105_SIM_SAMPLE_ROWS:-100}"
sim_group="${RJ_THE105_SIM_GROUP_SIZE:-${sim_sample_rows}}"
sim_events="${RJ_THE105_SIM_EVENTS:-20000}"
sampled_list_root="${RJ_THE105_SAMPLED_LIST_ROOT:-${repo_root}/.recoiljets_tmp/the105_shower_contract_sampled_lists/${campaign_tag}}"
data_memory="${RJ_THE105_DATA_MEMORY:-12000MB}"
sim_memory="${RJ_THE105_SIM_MEMORY:-12000MB}"
memory_retry_cap_mb="${RJ_THE105_MEMORY_RETRY_CAP_MB:-16000}"
allow_existing="${RJ_THE105_ALLOW_EXISTING:-0}"
runtime_trace="${RJ_THE105_RUNTIME_TRACE:-0}"
[[ "$runtime_trace" == "0" || "$runtime_trace" == "1" ]] || \
  die "RJ_THE105_RUNTIME_TRACE must be 0 or 1, got: ${runtime_trace}"
runtime_trace_verbosity="${RJ_THE105_RUNTIME_TRACE_VERBOSITY:-1}"
[[ "$runtime_trace_verbosity" =~ ^[0-9]+$ ]] || \
  die "RJ_THE105_RUNTIME_TRACE_VERBOSITY must be a non-negative integer, got: ${runtime_trace_verbosity}"
[[ "$sim_sample_rows" =~ ^[0-9]+$ && "$sim_sample_rows" -gt 0 ]] || \
  die "RJ_THE105_SIM_SAMPLE_ROWS must be a positive integer, got: ${sim_sample_rows}"
[[ "$sim_group" =~ ^[0-9]+$ && "$sim_group" -gt 0 ]] || \
  die "RJ_THE105_SIM_GROUP_SIZE must be a positive integer, got: ${sim_group}"
[[ "$sim_group" == "$sim_sample_rows" ]] || \
  die "Diagnostic stratification requires RJ_THE105_SIM_GROUP_SIZE (${sim_group}) to equal RJ_THE105_SIM_SAMPLE_ROWS (${sim_sample_rows})"
[[ "$data_runs" =~ ^[0-9]+$ && "$data_runs" -gt 0 ]] || \
  die "RJ_THE105_DATA_RUNS must be a positive integer, got: ${data_runs}"
[[ "$data_max_jobs" =~ ^[0-9]+$ && "$data_max_jobs" -gt 0 ]] || \
  die "RJ_THE105_DATA_MAX_JOBS must be a positive integer, got: ${data_max_jobs}"

variants=(historical towerinfo70 canonical)
signal_samples=(run28_embeddedPhoton12 run28_embeddedPhoton20)
inclusive_samples=(run28_embeddedJet12 run28_embeddedJet20 run28_embeddedJet30 run28_embeddedJet40)
if [[ "${RJ_THE105_VARIANTS+x}" == "x" ]]; then
  variants=()
  [[ -n "$RJ_THE105_VARIANTS" ]] && read -r -a variants <<< "$RJ_THE105_VARIANTS"
fi
if [[ "${RJ_THE105_SIGNAL_SAMPLES+x}" == "x" ]]; then
  signal_samples=()
  [[ -n "$RJ_THE105_SIGNAL_SAMPLES" ]] && read -r -a signal_samples <<< "$RJ_THE105_SIGNAL_SAMPLES"
fi
if [[ "${RJ_THE105_INCLUSIVE_SAMPLES+x}" == "x" ]]; then
  inclusive_samples=()
  [[ -n "$RJ_THE105_INCLUSIVE_SAMPLES" ]] && read -r -a inclusive_samples <<< "$RJ_THE105_INCLUSIVE_SAMPLES"
fi
include_data="${RJ_THE105_INCLUDE_DATA:-1}"
[[ "$include_data" == "0" || "$include_data" == "1" ]] || \
  die "RJ_THE105_INCLUDE_DATA must be 0 or 1, got: ${include_data}"
(( ${#variants[@]} > 0 )) || die "RJ_THE105_VARIANTS selected zero variants"

export RJ_CODEX_CHAT_NAME="${RJ_CODEX_CHAT_NAME:-THE-105 | AuAu Shower-Contract Factorial Canary}"
export RJ_CODEX_THREAD_ID="${RJ_CODEX_THREAD_ID:-019f4c73-f704-78f0-9731-19e9cc9560c3}"

usage() {
  cat <<'EOF'
Usage: submit_the105_auau_shower_contract_factorial.sh MODE

Modes:
  print            Print the frozen three-arm diagnostic contract.
  build-isolated   Build a PhotonClusterBuilder-only override and libRecoilJetsAuAu.
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
  [[ -s "$calo_library" ]] || die "Missing known-good CaloReco library: ${calo_library}."
  [[ -s "$photon_builder_library" ]] || die "Missing PhotonClusterBuilder override: ${photon_builder_library}. Run build-isolated first."
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

prepare_stratified_sim_lists() {
  local -a samples=("${signal_samples[@]}" "${inclusive_samples[@]}")
  (( ${#samples[@]} > 0 )) || return 0

  case "$sampled_list_root" in
    "${repo_root}/.recoiljets_tmp/"*) ;;
    *) die "Refusing sampled-list output outside ${repo_root}/.recoiljets_tmp: ${sampled_list_root}" ;;
  esac

  local sample source_dir target_dir reference_file total_rows file file_rows
  local -a matched_files=(
    DST_CALO_CLUSTER.matched.list
    G4Hits.matched.list
    DST_JETS.matched.list
    DST_GLOBAL.matched.list
    DST_MBD_EPD.matched.list
  )

  mkdir -p "$sampled_list_root"
  : > "${sampled_list_root}/sampled_row_indices.tsv"
  printf 'sample\tsource_rows\tsampled_rows\trow_index\n' \
    >> "${sampled_list_root}/sampled_row_indices.tsv"

  for sample in "${samples[@]}"; do
    source_dir="${repo_root}/simListFiles/${sample}"
    target_dir="${sampled_list_root}/${sample}"
    reference_file="${source_dir}/G4Hits.matched.list"
    [[ -s "$reference_file" ]] || die "Missing canonical matched list: ${reference_file}"
    total_rows="$(wc -l < "$reference_file" | tr -d ' ')"
    (( sim_sample_rows <= total_rows )) || \
      die "Requested ${sim_sample_rows} sampled rows from ${sample}, but only ${total_rows} exist"
    mkdir -p "$target_dir"

    awk -v sample="$sample" -v n="$total_rows" -v k="$sim_sample_rows" '
      BEGIN {
        for (i = 0; i < k; ++i) {
          row = int((i + 0.5) * n / k) + 1;
          printf "%s\t%d\t%d\t%d\n", sample, n, k, row;
        }
      }
    ' >> "${sampled_list_root}/sampled_row_indices.tsv"

    for file in "${matched_files[@]}"; do
      [[ -s "${source_dir}/${file}" ]] || die "Missing canonical matched list: ${source_dir}/${file}"
      file_rows="$(wc -l < "${source_dir}/${file}" | tr -d ' ')"
      [[ "$file_rows" == "$total_rows" ]] || \
        die "Matched-list row mismatch for ${sample}/${file}: ${file_rows} vs ${total_rows}"
      awk -v n="$total_rows" -v k="$sim_sample_rows" '
        BEGIN {
          for (i = 0; i < k; ++i) {
            row = int((i + 0.5) * n / k) + 1;
            keep[row] = 1;
          }
        }
        keep[FNR]
      ' "${source_dir}/${file}" > "${target_dir}/${file}"
      [[ "$(wc -l < "${target_dir}/${file}" | tr -d ' ')" == "$sim_sample_rows" ]] || \
        die "Failed to materialize ${sim_sample_rows} rows for ${sample}/${file}"
    done
  done

  say "Prepared ${sim_sample_rows} evenly spaced matched rows/sample under ${sampled_list_root}"
}

write_manifest() {
  mkdir -p "$evidence_dir"
  local photon_builder_cc="src/PhotonClusterBuilder.cc"
  local photon_builder_h="src/PhotonClusterBuilder.h"
  if [[ ! -s "$photon_builder_cc" || ! -s "$photon_builder_h" ]]; then
    photon_builder_cc="coresoftware_local/offline/packages/CaloReco/PhotonClusterBuilder.cc"
    photon_builder_h="coresoftware_local/offline/packages/CaloReco/PhotonClusterBuilder.h"
  fi
  [[ -s "$photon_builder_cc" ]] || die "Missing mapped PhotonClusterBuilder source: ${photon_builder_cc}"
  [[ -s "$photon_builder_h" ]] || die "Missing mapped PhotonClusterBuilder header: ${photon_builder_h}"
  {
    printf 'variant\tpopulation\treconstruction_source\ttower_floor_gev\tselection_stages\n'
    printf 'historical\tdata\tTowerInfo grid plus historical local chi2/CDB mask\t0.070\tbefore,after_NCB_preselection,after_tight\n'
    printf 'historical\tembedded\tRawCluster towermap membership\t0.070\tbefore,after_NCB_preselection,after_tight\n'
    printf 'towerinfo70\tall\tTowerInfo full grid with get_isGood acceptance\t0.070\tbefore,after_NCB_preselection,after_tight\n'
    printf 'canonical\tall\tTowerInfo full grid with get_isGood acceptance\t0.000\tbefore,after_NCB_preselection,after_tight\n'
  } > "${evidence_dir}/campaign_contract.tsv"
  {
    printf 'campaign_tag=%s\n' "$campaign_tag"
    printf 'source_checkout_head=%s\n' "$(git rev-parse HEAD 2>/dev/null || printf unknown)"
    printf 'config=%s\n' "$yaml"
    printf 'variants=%s\n' "${variants[*]}"
    printf 'include_data=%s\n' "$include_data"
    printf 'signal_samples=%s\n' "${signal_samples[*]:-none}"
    printf 'inclusive_samples=%s\n' "${inclusive_samples[*]:-none}"
    if [[ "$include_data" == "1" ]]; then
      printf 'data_contract=%s_largest_stat_GRL_run_pool_groupSize_%s_global_job_cap_%s_up_to_%s_events_each\n' \
        "$data_runs" "$data_group" "$data_max_jobs" "$data_events"
      printf 'data_sampling_note=job_cap_is_applied_across_groups_in_run_list_order_and_does_not_guarantee_one_job_per_selected_run\n'
    else
      printf 'data_contract=disabled\n'
    fi
    printf 'sim_contract=%s_evenly_spaced_paired_rows_per_sample_one_group_up_to_%s_events\n' "$sim_sample_rows" "$sim_events"
    printf 'sim_sampled_list_root=%s\n' "$sampled_list_root"
    printf 'sim_sampled_row_index_manifest=%s\n' "${sampled_list_root}/sampled_row_indices.tsv"
    printf 'analysis_mode=candidate_skim_only_no_truth_matching_no_training_tree_no_production_histogram_booking\n'
    printf 'candidate_contract=15<=ET<35,abs_eta<0.7,one_skim_row_per_candidate\n'
    printf 'normalization_contract=full_finite_candidate_denominator_with_zero_underflow_overflow_reported_separately\n'
  } > "${evidence_dir}/campaign_manifest.txt"
  local -a source_files=(
    macros/Fun4All_recoilJets_unified_impl.C
    "$photon_builder_cc"
    "$photon_builder_h"
    src_AuAu/RecoilJets_AuAu.cc
    src_AuAu/RecoilJets_AuAu.h
    scripts/sdcc/runtime/condor/RecoilJets_Condor_submit.sh
    scripts/sdcc/workflows/diagnostics/submit_the105_auau_shower_contract_factorial.sh
    scripts/sdcc/workflows/diagnostics/the105_preserve_invalid_shower_shapes.patch
  )
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "${source_files[@]}" > "${evidence_dir}/source_files.sha256"
  else
    shasum -a 256 "${source_files[@]}" > "${evidence_dir}/source_files.sha256"
  fi
}

print_contract() {
  write_manifest
  local data_jobs=0
  local data_description="disabled"
  if [[ "$include_data" == "1" ]]; then
    data_jobs="$data_max_jobs"
    data_description="${data_runs}-run largest-statistics pool, groupSize ${data_group}, global cap ${data_max_jobs} jobs in run-list order, <=${data_events} events/job"
  fi
  cat <<EOF
RECOILJETS_THE105_AUAU_SHOWER_CONTRACT_FACTORIAL_V1
campaign_tag=${campaign_tag}
config=${yaml}
variants=${variants[*]}
data=${data_description}
signal_samples=${signal_samples[*]-}
inclusive_samples=${inclusive_samples[*]-}
sim=${sim_sample_rows} evenly spaced paired rows/sample, one ${sim_group}-row group, <=${sim_events} events/job
stages=before preselection; after complete NCB preselection; after frozen tight ID
analysis_mode=candidate-skim-only; truth matching, training tree, and production histogram suite disabled
candidate_skim=enabled, one row per candidate
maximum_jobs=$(( ${#variants[@]} * (data_jobs + ${#signal_samples[@]} + ${#inclusive_samples[@]}) ))
automatic_merge=disabled
canonical_replacement=forbidden
canary_root=${canary_root}
isolated_calo_library=${calo_library}
isolated_photon_builder_override=${photon_builder_library}
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

build_photon_builder_override() {
  local source_dir="${repo_root}/coresoftware_local/offline/packages/CaloReco"
  local source_cc="${source_dir}/PhotonClusterBuilder.cc"
  local source_h="${source_dir}/PhotonClusterBuilder.h"
  [[ -s "$source_cc" && -s "$source_h" ]] || \
    die "Missing PhotonClusterBuilder source under ${source_dir}"

  reset_build_dir "$photon_builder_build"
  reset_build_dir "$photon_builder_install"
  mkdir -p "${photon_builder_install}/lib" \
           "${photon_builder_install}/include/caloreco"
  cp -f "$source_cc" "$source_h" "$photon_builder_build/"

  local -a link_libs=(
    -lmbd_io -lcalo_io -lcdbobjects -lCLHEP -lffamodules
    -lffarawobjects -lgsl -lgslcblas -lglobalvertex_io -lsph_onnx
    -lphg4hit -lphparameter_io -lphool -lSubsysReco -lTMVA -lTMVAUtils
  )
  (
    cd "$photon_builder_build"
    g++ -std=c++20 -O2 -g -fPIC -shared -Wl,-z,defs \
      -I. \
      -I/sphenix/u/${USER:-patsfan753}/thesisAnalysis/install/include \
      -isystem "${OFFLINE_MAIN}/include" \
      -isystem "${ROOTSYS}/include" \
      PhotonClusterBuilder.cc \
      -L/sphenix/u/${USER:-patsfan753}/thesisAnalysis/install/lib \
      -L"${OFFLINE_MAIN}/lib64" -L"${OFFLINE_MAIN}/lib" \
      "${link_libs[@]}" $(root-config --libs) \
      -Wl,-soname,libphoton_cluster_builder_override.so \
      -o "${photon_builder_install}/lib/libphoton_cluster_builder_override.so"
    cp -f PhotonClusterBuilder.h \
      "${photon_builder_install}/include/caloreco/PhotonClusterBuilder.h"
  ) 2>&1 | tee "${evidence_dir}/photon_builder_override_build.log"
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
  build_photon_builder_override
  build_component "${repo_root}/src_AuAu" "$auau_build" "$auau_install" recoiljets_auau
  [[ -s "$calo_library" ]] || die "Known-good CaloReco library is missing: ${calo_library}"
  [[ -s "$photon_builder_library" ]] || die "PhotonClusterBuilder override build did not produce ${photon_builder_library}"
  [[ -s "$calo_header" ]] || die "PhotonClusterBuilder override did not install ${calo_header}"
  [[ -s "$auau_library" ]] || die "AuAu build did not produce ${auau_library}"
  sha256sum "$calo_library" "$photon_builder_library" "$calo_header" "$auau_library" | tee "${evidence_dir}/isolated_runtime.sha256"
}

common_env() {
  local variant="$1"
  local trace_suffix=""
  if [[ "$runtime_trace" == "1" ]]; then
    trace_suffix=";RJ_VERBOSITY=${runtime_trace_verbosity};RJ_F4A_VERBOSE=1;RJ_STEP_EVENTS=1"
  fi
  printf '%s\n' \
    "RJ_CONFIG_YAML=${yaml}" \
    "RJ_CALO_RECO_LIBRARY_OVERRIDE=${calo_library}" \
    "RJ_PHOTON_CLUSTER_BUILDER_LIBRARY_OVERRIDE=${photon_builder_library}" \
    "RJ_PHOTON_CLUSTER_BUILDER_HEADER_OVERRIDE=${calo_header}" \
    "RJ_AUAU_LIBRARY_OVERRIDE=${auau_library}" \
    "RJ_SUBMIT_EXTRA_ENV=RJ_REQUIRE_EMBEDDED_MINBIAS_CLASSIFIER=1;RJ_AUAU_SHOWER_SHAPE_DIAGNOSTIC_VARIANT=${variant};RJ_AUAU_CANDIDATE_SKIM_ONLY=1;RJ_AUAU_BDT_EXTRACT_ONLY=1;RJ_AUAU_PHOTON_CANDIDATE_SKIM=1;RJ_AUAU_PHOTON_CANDIDATE_SKIM_MAX_ENTRIES=0${trace_suffix}" \
    "RJ_ID_FANOUT_MAX_ROWS=1" \
    "RJ_SIM_ROOT_OVERRIDE=${sampled_list_root}" \
    "RJ_PHOTON_ID_ROW_MATCH=preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTSideband" \
    "RJ_AUTO_MERGE=0" \
    "RJ_AUTO_MEMORY_RETRY_CAP_MB=${memory_retry_cap_mb}" \
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

  if [[ "$include_data" == "1" ]]; then
    say "Expanding ${variant}/data"
    env "${env_args[@]}" \
      "RJ_DAG_DRYRUN=${dryrun}" \
      "RJ_REQUEST_MEMORY=${data_memory}" \
      "RJ_SMOKE_OUTPUT_BASE=${variant_root}/data" \
      "RJ_SMOKE_DATA_RUNS=${data_runs}" \
      "RJ_SMOKE_DATA_MAX_JOBS=${data_max_jobs}" \
      "RJ_SMOKE_DATA_NEVENTS=${data_events}" \
      ./RecoilJets_Condor_submit.sh isAuAu condor smokeTest groupSize "$data_group"
  fi

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
  prepare_stratified_sim_lists
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
