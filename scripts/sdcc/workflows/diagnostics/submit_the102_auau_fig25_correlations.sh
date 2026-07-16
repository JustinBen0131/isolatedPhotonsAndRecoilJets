#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
if [[ -f "${script_dir}/../../runtime/io/recoiljets_io_paths.sh" ]]; then
  source "${script_dir}/../../runtime/io/recoiljets_io_paths.sh"
elif [[ -f scripts/sdcc/runtime/io/recoiljets_io_paths.sh ]]; then
  source scripts/sdcc/runtime/io/recoiljets_io_paths.sh
else
  echo "[THE102][ERROR] Cannot locate recoiljets_io_paths.sh" >&2
  exit 2
fi

repo_root="$(rj_find_repo_root "$script_dir")"
cd "$repo_root"

say() { printf '[THE102] %s\n' "$*"; }
die() { printf '[THE102][ERROR] %s\n' "$1" >&2; exit "${2:-2}"; }

mode="${1:-print}"
campaign_tag="${RJ_THE102_CAMPAIGN_TAG:-the102_auau_fig25_correlations_20260714}"
yaml="${RJ_THE102_CONFIG_YAML:-${repo_root}/macros/analysis_config_the102_auau_fig25_correlations.yaml}"
build_root="${RJ_THE102_BUILD_ROOT:-${repo_root}/.recoiljets_tmp/the102_fig25_build_20260714}"
build_dir="${build_root}/build"
install_dir="${build_root}/install"
auau_library="${RJ_THE102_AUAU_LIBRARY:-${install_dir}/lib/libRecoilJetsAuAu.so}"
evidence_dir="${RJ_THE102_EVIDENCE_DIR:-${repo_root}/evidence/qa/the102_fig25_20260714}"
bulk_root="${RJ_THE102_BULK_ROOT:-$(rj_recoiljets_bulk_root)/auau_fig25/${campaign_tag}}"
canary_root="${RJ_THE102_CANARY_ROOT:-$(rj_recoiljets_bulk_root)/smoke/auau_fig25/${campaign_tag}}"
merge_base="${RJ_THE102_MERGE_BASE:-$(rj_merge_root_for_campaign "$campaign_tag")/inclusive}"
group_size="${RJ_THE102_GROUP_SIZE:-10}"
memory="${RJ_THE102_MEMORY:-6500MB}"
canary_events="${RJ_THE102_CANARY_EVENTS:-10000}"
allow_existing="${RJ_THE102_ALLOW_EXISTING:-0}"
samples=(run28_embeddedJet12 run28_embeddedJet20 run28_embeddedJet30 run28_embeddedJet40)

export RJ_CODEX_CHAT_NAME="${RJ_CODEX_CHAT_NAME:-THE-102 | AuAu Fig25 Correlations}"
export RJ_CODEX_THREAD_ID="${RJ_CODEX_THREAD_ID:-019f61f3-e7bc-7cb1-ae5d-48ed6d88eb95}"

usage() {
  cat <<'EOF'
Usage: submit_the102_auau_fig25_correlations.sh MODE

Modes:
  print             Print the immutable campaign contract.
  build-isolated    Build a campaign-only libRecoilJetsAuAu.so.
  canary-submit     Submit one Jet12/20/30/40 canary each.
  production-count Print the focused four-sample production counts.
  production-submit Submit Jet12/20/30/40 with 10 input groups per job.
  status            Report matching queue rows and ROOT coverage.

Automatic merge is disabled. The embedded MinimumBiasClassifier filter is
forced off even though the YAML retains the data-production classifier flag.
EOF
}

require_inputs() {
  rj_validate_campaign_tag "$campaign_tag" || exit 2
  [[ -s "$yaml" ]] || die "Missing config: ${yaml}"
  [[ -s "$auau_library" ]] || die "Missing isolated AuAu library: ${auau_library}"
  [[ -x ./RecoilJets_Condor_submit.sh ]] || die "Missing RecoilJets_Condor_submit.sh"
}

assert_no_active_duplicate() {
  local output_root="$1"
  if condor_q "${USER:-patsfan753}" -af Args 2>/dev/null | grep -F "$output_root" >/dev/null; then
    die "Active Condor jobs already target ${output_root}" 4
  fi
}

common_env() {
  printf '%s\n' \
    "RJ_CONFIG_YAML=${yaml}" \
    "RJ_AUAU_LIBRARY_OVERRIDE=${auau_library}" \
    "RJ_SUBMIT_EXTRA_ENV=RJ_AUAU_FIG25_CORRELATION_DIAGNOSTICS=1;RJ_REQUIRE_EMBEDDED_MINBIAS_CLASSIFIER=0" \
    "RJ_ID_FANOUT_MAX_ROWS=1" \
    "RJ_AUTO_MERGE=0" \
    "RJ_REQUIRE_NON_TINY_OUTPUT=1" \
    "RJ_MIN_OUTPUT_BYTES=50000" \
    "RJ_PROFILE_JOB=1" \
    "RJ_JOB_HEARTBEAT_SECONDS=120" \
    "RJ_INTERNAL_FIXED_ISO_GEV_AUAU=4.0" \
    "RJ_AUAU_BUILD_TOPOCLUSTER_ISOLATION=0" \
    "RJ_AUAU_USE_TOPOCLUSTER_ISOLATION=0" \
    "RJ_SIMEMBEDDEDINCLUSIVE_FOUR_SAMPLES=1" \
    "RJ_CODEX_CHAT_NAME=${RJ_CODEX_CHAT_NAME}" \
    "RJ_CODEX_THREAD_ID=${RJ_CODEX_THREAD_ID}"
}

print_contract() {
  cat <<EOF
RECOILJETS_THE102_AUAU_FIG25_V1
campaign_tag=${campaign_tag}
config=${yaml}
isolated_library=${auau_library}
samples=${samples[*]}
model=THE-95 no-veto baseline14 pT15to35
embedded_minbias_classifier_filter=OFF
photon_pt=15 <= pT < 35 GeV
centrality=0-20,20-50,50-80
surfaces=e11/e33_vs_Eiso_background,raw_BDT_score_vs_Eiso_background
group_size=${group_size}
bulk_root=${bulk_root}
canary_root=${canary_root}
merge_base=${merge_base}
automatic_merge=disabled
RJ_CODEX_CHAT_NAME=${RJ_CODEX_CHAT_NAME}
RJ_CODEX_THREAD_ID=${RJ_CODEX_THREAD_ID}
EOF
}

build_isolated() {
  [[ -s "${repo_root}/src_AuAu/autogen.sh" ]] || die "Missing src_AuAu/autogen.sh"
  case "$build_root" in
    "${repo_root}/.recoiljets_tmp/the102_fig25_build_"*) ;;
    *) die "Refusing build outside THE-102 hidden path: ${build_root}" ;;
  esac
  rm -rf "$build_dir" "$install_dir"
  mkdir -p "$build_dir" "$install_dir" "$evidence_dir"
  cp -a "${repo_root}/src_AuAu/." "$build_dir/"
  rm -rf "${build_dir}/autom4te.cache" "${build_dir}/.deps"
  rm -f "${build_dir}/config.status" "${build_dir}/config.log" \
        "${build_dir}/Makefile" "${build_dir}/libtool" "${build_dir}/stamp-h1"
  set +u
  source /opt/sphenix/core/bin/sphenix_setup.sh -n
  source /opt/sphenix/core/bin/setup_local.sh "${RJ_THE102_DEP_INSTALL:-/sphenix/u/patsfan753/thesisAnalysis/install}"
  set -u
  (
    cd "$build_dir"
    bash ./autogen.sh --prefix="$install_dir"
    make -j"${RJ_THE102_BUILD_JOBS:-4}"
    make install
  ) 2>&1 | tee "${evidence_dir}/isolated_build.log"
  [[ -s "$auau_library" ]] || die "Build did not produce ${auau_library}"
  sha256sum "$auau_library" | tee "${evidence_dir}/isolated_library.sha256"
}

submit_sample() {
  local sample="$1" output_root="$2" smoke="$3"
  local -a env_args smoke_args
  mapfile -t env_args < <(common_env)
  smoke_args=()
  if [[ "$smoke" == "1" ]]; then
    smoke_args+=("RJ_SMOKE_OUTPUT_BASE=${output_root}" "RJ_SMOKE_SIM_NEVENTS=${canary_events}")
    env "${env_args[@]}" "${smoke_args[@]}" "RJ_REQUEST_MEMORY=${memory}" \
      ./RecoilJets_Condor_submit.sh isSimEmbeddedInclusive condorDoAllSmoke groupSize 1 maxJobs 1 "SAMPLE=${sample}"
  else
    env "${env_args[@]}" "RJ_REQUEST_MEMORY=${memory}" \
      "RJ_SIMEMBEDINCLUSIVE_DEST_BASE=${output_root}" \
      "RJ_MERGE_OUT_BASE_OVERRIDE=${merge_base}" \
      ./RecoilJets_Condor_submit.sh isSimEmbeddedInclusive condorDoAll groupSize "$group_size" "SAMPLE=${sample}"
  fi
}

submit_canaries() {
  require_inputs
  assert_no_active_duplicate "$canary_root"
  [[ ! -e "$canary_root" || "$allow_existing" == "1" ]] || die "Canary root exists: ${canary_root}"
  mkdir -p "$evidence_dir"
  exec > >(tee -a "${evidence_dir}/canary_submission.log") 2>&1
  for sample in "${samples[@]}"; do submit_sample "$sample" "$canary_root" 1; done
}

production_count() {
  require_inputs
  local -a env_args
  mapfile -t env_args < <(common_env)
  for sample in "${samples[@]}"; do
    say "Counting ${sample}"
    env "${env_args[@]}" "RJ_REQUEST_MEMORY=${memory}" \
      ./RecoilJets_Condor_submit.sh isSimEmbeddedInclusive CHECKJOBS groupSize "$group_size" "SAMPLE=${sample}"
  done
}

submit_production() {
  require_inputs
  assert_no_active_duplicate "$bulk_root"
  [[ ! -e "$bulk_root" || "$allow_existing" == "1" ]] || die "Production root exists: ${bulk_root}"
  mkdir -p "$evidence_dir"
  exec > >(tee -a "${evidence_dir}/production_submission.log") 2>&1
  for sample in "${samples[@]}"; do submit_sample "$sample" "$bulk_root" 0; done
  condor_q "${USER:-patsfan753}" -af ClusterId JobStatus Args 2>/dev/null | grep -F "$campaign_tag" || true
}

status() {
  print_contract
  condor_q "${USER:-patsfan753}" -af ClusterId ProcId JobStatus Args 2>/dev/null | grep -F "$campaign_tag" || true
  if [[ -d "$bulk_root" ]]; then
    find "$bulk_root" -type f -name '*.root' -size +50k -printf '%p\n' 2>/dev/null |
      awk -F/ '{for(i=1;i<=NF;i++) if($i ~ /embeddedJet(12|20|30|40)/) c[$i]++} END{for(k in c) print k,c[k]}' | sort
  fi
}

case "$mode" in
  print) print_contract ;;
  build-isolated) build_isolated ;;
  canary-submit) submit_canaries ;;
  production-count) production_count ;;
  production-submit) submit_production ;;
  status) status ;;
  -h|--help|help) usage ;;
  *) usage; die "Unknown mode: ${mode}" ;;
esac
