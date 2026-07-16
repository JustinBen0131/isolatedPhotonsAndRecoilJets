#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
if [[ -f "${script_dir}/../../runtime/io/recoiljets_io_paths.sh" ]]; then
  source "${script_dir}/../../runtime/io/recoiljets_io_paths.sh"
elif [[ -f scripts/sdcc/runtime/io/recoiljets_io_paths.sh ]]; then
  source scripts/sdcc/runtime/io/recoiljets_io_paths.sh
else
  echo "[THE104][ERROR] Cannot locate recoiljets_io_paths.sh" >&2
  exit 2
fi

repo_root="$(rj_find_repo_root "$script_dir")"
cd "$repo_root"

say() { printf '[THE104] %s\n' "$*"; }
die() { printf '[THE104][ERROR] %s\n' "$1" >&2; exit "${2:-2}"; }

mode="${1:-print}"
campaign_tag="${RJ_THE104_CAMPAIGN_TAG:-the104_canonical_minbias_isolation_20260715}"
base_yaml="${RJ_THE104_BASE_YAML:-${repo_root}/macros/analysis_config_the104_canonical_minbias_isolation_base.yaml}"
phosub_yaml="${RJ_THE104_PHOSUB_YAML:-${repo_root}/macros/analysis_config_the104_canonical_minbias_isolation_phosub.yaml}"
build_root="${RJ_THE104_BUILD_ROOT:-${repo_root}/.recoiljets_tmp/the104_minbias_isolation_build_20260715}"
build_dir="${build_root}/build"
install_dir="${build_root}/install"
auau_library="${RJ_THE104_AUAU_LIBRARY:-${install_dir}/lib/libRecoilJetsAuAu.so}"
evidence_dir="${RJ_THE104_EVIDENCE_DIR:-${repo_root}/evidence/qa/the104_minbias_isolation_20260715}"
bulk_root="${RJ_THE104_BULK_ROOT:-$(rj_recoiljets_bulk_root)/auau_isolation/${campaign_tag}}"
canary_root="${RJ_THE104_CANARY_ROOT:-$(rj_recoiljets_bulk_root)/smoke/auau_isolation/${campaign_tag}}"
merge_root="${RJ_THE104_MERGE_ROOT:-$(rj_merge_root_for_campaign "$campaign_tag")}"
group_size="${RJ_THE104_GROUP_SIZE:-7}"
base_memory="${RJ_THE104_BASE_MEMORY:-6500MB}"
phosub_memory="${RJ_THE104_PHOSUB_MEMORY:-12000MB}"
canary_events="${RJ_THE104_CANARY_EVENTS:-10000}"
allow_existing="${RJ_THE104_ALLOW_EXISTING:-0}"
samples=(run28_embeddedPhoton12 run28_embeddedPhoton20)
variants=(baseVariant variantB)

export RJ_CODEX_CHAT_NAME="${RJ_CODEX_CHAT_NAME:-THE-104 | Canonical MinBias Isolation Slides}"
export RJ_CODEX_THREAD_ID="${RJ_CODEX_THREAD_ID:-019f663e-8454-7b10-a1f0-575ae832e84f}"

usage() {
  cat <<'EOF'
Usage: submit_the104_canonical_minbias_isolation_slides.sh MODE [VARIANT] [STAGE]

Modes:
  print                  Print the immutable campaign contract.
  build-isolated         Build a campaign-only libRecoilJetsAuAu.so.
  canary-submit          Submit Photon12/20 canaries for baseVariant and variantB.
  production-count       Print exact job counts for all four production lanes.
  production-submit      Submit all four production lanes.
  production-submit-lane VARIANT SAMPLE
                         Submit one missing lane after a queue-cap gate.
  merge-submit VARIANT STAGE
                         Submit a guarded merge stage: firstRound, secondRound,
                         or finalStitch. The caller must first prove drained raw
                         coverage and the preceding merge stage.
  status                 Report matching queue rows and non-tiny ROOT coverage.

Automatic merge is disabled. Every worker requires the embedded
MinimumBiasClassifier before physics histogram filling.
EOF
}

yaml_for_variant() {
  case "$1" in
    baseVariant) printf '%s\n' "$base_yaml" ;;
    variantB) printf '%s\n' "$phosub_yaml" ;;
    *) die "Unknown variant: $1" ;;
  esac
}

memory_for_variant() {
  case "$1" in
    baseVariant) printf '%s\n' "$base_memory" ;;
    variantB) printf '%s\n' "$phosub_memory" ;;
    *) die "Unknown variant: $1" ;;
  esac
}

require_inputs() {
  rj_validate_campaign_tag "$campaign_tag" || exit 2
  [[ -s "$base_yaml" ]] || die "Missing baseline config: ${base_yaml}"
  [[ -s "$phosub_yaml" ]] || die "Missing PHOSUB config: ${phosub_yaml}"
  [[ -s "$auau_library" ]] || die "Missing isolated AuAu library: ${auau_library}. Run build-isolated first."
  [[ -x ./RecoilJets_Condor_submit.sh ]] || die "Missing RecoilJets_Condor_submit.sh"
}

assert_no_active_duplicate() {
  local output_root="$1"
  if condor_q "${USER:-patsfan753}" -af Args 2>/dev/null | grep -F "$output_root" >/dev/null; then
    die "Active Condor jobs already target ${output_root}" 4
  fi
}

assert_fresh_path() {
  local path="$1" label="$2"
  if [[ -e "$path" && "$allow_existing" != "1" ]]; then
    die "${label} exists: ${path}. Use a fresh tag or explicitly set RJ_THE104_ALLOW_EXISTING=1 for a proven resume."
  fi
}

common_env() {
  local yaml="$1"
  printf '%s\n' \
    "RJ_CONFIG_YAML=${yaml}" \
    "RJ_AUAU_LIBRARY_OVERRIDE=${auau_library}" \
    "RJ_SUBMIT_EXTRA_ENV=RJ_REQUIRE_EMBEDDED_MINBIAS_CLASSIFIER=1" \
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

print_contract() {
  cat <<EOF
RECOILJETS_THE104_CANONICAL_MINBIAS_ISOLATION_V1
campaign_tag=${campaign_tag}
base_config=${base_yaml}
phosub_config=${phosub_yaml}
isolated_library=${auau_library}
samples=${samples[*]}
variants=${variants[*]}
canonical_bdt=/sphenix/user/patsfan753/forBlair/auauBDTdefault/auau_tight_bdt_centAsFeatBase3x3_pt15to35_tmva.root
canonical_wp80=T80(c)=0.5378806890+0.0011122896*c
embedded_minbias_classifier_filter=ON
photon_pt_edges=15,17,19,21,23,26,30,35
centrality_edges=0:5:80
isolation_cones=R0.3,R0.4
bulk_root=${bulk_root}
canary_root=${canary_root}
merge_root=${merge_root}
automatic_merge=disabled
RJ_CODEX_CHAT_NAME=${RJ_CODEX_CHAT_NAME}
RJ_CODEX_THREAD_ID=${RJ_CODEX_THREAD_ID}
EOF
}

build_isolated() {
  [[ -s "${repo_root}/src_AuAu/autogen.sh" ]] || die "Missing src_AuAu/autogen.sh"
  case "$build_root" in
    "${repo_root}/.recoiljets_tmp/the104_minbias_isolation_build_"*) ;;
    *) die "Refusing build outside THE-104 hidden path: ${build_root}" ;;
  esac
  rm -rf "$build_dir" "$install_dir"
  mkdir -p "$build_dir" "$install_dir" "$evidence_dir"
  cp -a "${repo_root}/src_AuAu/." "$build_dir/"
  rm -rf "${build_dir}/autom4te.cache" "${build_dir}/.deps"
  rm -f "${build_dir}/config.status" "${build_dir}/config.log" \
        "${build_dir}/Makefile" "${build_dir}/libtool" "${build_dir}/stamp-h1"
  set +u
  source /opt/sphenix/core/bin/sphenix_setup.sh -n
  source /opt/sphenix/core/bin/setup_local.sh "${RJ_THE104_DEP_INSTALL:-/sphenix/u/patsfan753/thesisAnalysis/install}"
  set -u
  (
    cd "$build_dir"
    bash ./autogen.sh --prefix="$install_dir"
    make -j"${RJ_THE104_BUILD_JOBS:-4}"
    make install
  ) 2>&1 | tee "${evidence_dir}/isolated_build.log"
  [[ -s "$auau_library" ]] || die "Build did not produce ${auau_library}"
  sha256sum "$auau_library" | tee "${evidence_dir}/isolated_library.sha256"
}

submit_lane() {
  local variant="$1" sample="$2" smoke="$3"
  local yaml memory output_root
  local -a env_args
  yaml="$(yaml_for_variant "$variant")"
  memory="$(memory_for_variant "$variant")"
  if [[ "$smoke" == "1" ]]; then
    output_root="${canary_root}/${variant}"
  else
    output_root="${bulk_root}/${variant}"
  fi
  mapfile -t env_args < <(common_env "$yaml")
  if [[ "$smoke" == "1" ]]; then
    env "${env_args[@]}" \
      "RJ_REQUEST_MEMORY=${memory}" \
      "RJ_SMOKE_OUTPUT_BASE=${output_root}" \
      "RJ_SMOKE_SIM_NEVENTS=${canary_events}" \
      ./RecoilJets_Condor_submit.sh isSimEmbedded condorDoAllSmoke groupSize 1 maxJobs 1 "SAMPLE=${sample}"
  else
    env "${env_args[@]}" \
      "RJ_REQUEST_MEMORY=${memory}" \
      "RJ_SIMEMBED_DEST_BASE=${output_root}" \
      "RJ_MERGE_OUT_BASE_OVERRIDE=${merge_root}/${variant}" \
      ./RecoilJets_Condor_submit.sh isSimEmbedded condorDoAll groupSize "$group_size" "SAMPLE=${sample}"
  fi
}

submit_canaries() {
  require_inputs
  assert_no_active_duplicate "$canary_root"
  assert_fresh_path "$canary_root" "canary output root"
  mkdir -p "$evidence_dir"
  exec > >(tee -a "${evidence_dir}/canary_submission.log") 2>&1
  for variant in "${variants[@]}"; do
    for sample in "${samples[@]}"; do
      say "Submitting ${variant}/${sample} canary"
      submit_lane "$variant" "$sample" 1
    done
  done
}

production_count() {
  require_inputs
  local -a env_args
  local yaml memory
  for variant in "${variants[@]}"; do
    yaml="$(yaml_for_variant "$variant")"
    memory="$(memory_for_variant "$variant")"
    mapfile -t env_args < <(common_env "$yaml")
    for sample in "${samples[@]}"; do
      say "Counting ${variant}/${sample}"
      env "${env_args[@]}" "RJ_REQUEST_MEMORY=${memory}" \
        ./RecoilJets_Condor_submit.sh isSimEmbedded CHECKJOBS groupSize "$group_size" "SAMPLE=${sample}"
    done
  done
}

submit_one_production_lane() {
  local variant="$1" sample="$2" output_root="${bulk_root}/${variant}"
  case "$sample" in run28_embeddedPhoton12|run28_embeddedPhoton20) ;; *) die "Unknown signal sample: ${sample}" ;; esac
  yaml_for_variant "$variant" >/dev/null
  require_inputs
  if condor_q "${USER:-patsfan753}" -af Args 2>/dev/null | grep -F "$output_root" | grep -F "$sample" >/dev/null; then
    die "Active jobs already target ${variant}/${sample}" 4
  fi
  mkdir -p "$evidence_dir"
  submit_lane "$variant" "$sample" 0
}

submit_production() {
  require_inputs
  assert_no_active_duplicate "$bulk_root"
  assert_fresh_path "$bulk_root" "production bulk root"
  assert_fresh_path "$merge_root" "production merge root"
  mkdir -p "$evidence_dir"
  exec > >(tee -a "${evidence_dir}/production_submission.log") 2>&1
  for variant in "${variants[@]}"; do
    for sample in "${samples[@]}"; do
      say "Submitting ${variant}/${sample} production"
      submit_lane "$variant" "$sample" 0
    done
  done
  condor_q "${USER:-patsfan753}" -af ClusterId JobStatus Args 2>/dev/null | grep -F "$campaign_tag" || true
}

submit_merge_stage() {
  local variant="${2:-}" stage="${3:-}"
  local yaml input_root out_root
  yaml="$(yaml_for_variant "$variant")"
  case "$stage" in firstRound|secondRound|finalStitch) ;; *) die "Invalid merge stage: ${stage}" ;; esac
  input_root="${bulk_root}/${variant}"
  out_root="${merge_root}/${variant}"
  [[ -d "$input_root" ]] || die "Missing raw input root: ${input_root}"
  env MERGE_CONFIG_YAML="$yaml" \
      MERGE_SIM_INPUT_BASE_OVERRIDE="$input_root" \
      MERGE_OUT_BASE_OVERRIDE="$out_root" \
      RJ_SIM_FIRSTROUND_REQUEST_MEMORY="${RJ_THE104_MERGE_MEMORY:-16000MB}" \
      ./scripts/mergeRecoilJets.sh isSimEmbedded "$stage" groupSize "${RJ_THE104_MERGE_GROUP_SIZE:-75}"
}

status_report() {
  print_contract
  printf '\nQUEUE\n'
  condor_q "${USER:-patsfan753}" -af ClusterId ProcId JobStatus HoldReason Args 2>/dev/null | grep -F "$campaign_tag" || true
  printf '\nROOT_COUNTS\n'
  for variant in "${variants[@]}"; do
    for sample in "${samples[@]}"; do
      count="$(find "${bulk_root}/${variant}/${sample}" -type f -name '*.root' -size +50k 2>/dev/null | wc -l | tr -d ' ')"
      printf '%s\t%s\t%s\n' "$variant" "$sample" "$count"
    done
  done
}

case "$mode" in
  print) print_contract ;;
  build-isolated) build_isolated ;;
  canary-submit) submit_canaries ;;
  production-count) production_count ;;
  production-submit) submit_production ;;
  production-submit-lane) submit_one_production_lane "${2:-}" "${3:-}" ;;
  merge-submit) submit_merge_stage "$@" ;;
  status) status_report ;;
  help|-h|--help) usage ;;
  *) usage >&2; die "Unknown mode: ${mode}" ;;
esac
