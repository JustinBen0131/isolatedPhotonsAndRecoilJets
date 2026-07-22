#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd -P)"
cd "$repo_root"

mode="${1:-preflight}"
tag="${RJ_THE119_TAG:-the119_six_lane_replay_canary_20260720_v2}"
base="${RJ_THE119_OUTPUT_ROOT:-/sphenix/tg/tg01/bulk/jbennett/thesisAna/recoiljets/smoke/replay_foundation/${tag}}"
evidence="${RJ_THE119_EVIDENCE_ROOT:-/sphenix/u/patsfan753/scratch/thesisAnalysis/evidence/qa/${tag}}"
pp_cfg="${RJ_THE119_PP_CONFIG:-${repo_root}/macros/analysis_config_the119_pp_replay_foundation.yaml}"
auau_cfg="${RJ_THE119_AUAU_CONFIG:-${repo_root}/macros/analysis_config_the112_auau_combined_bdt_triplet.yaml}"
pp_lib="${RJ_THE119_PP_LIBRARY:-/sphenix/u/patsfan753/scratch/thesisAnalysis/.recoiljets_tmp/the118_replay_runtime_build_20260720/install_pp/lib/libRecoilJets.so}"
auau_lib="${RJ_THE119_AUAU_LIBRARY:-/sphenix/u/patsfan753/scratch/thesisAnalysis/.recoiljets_tmp/the118_replay_runtime_build_20260720/install_auau/lib/libRecoilJetsAuAu.so}"
pp_lib_sha="${RJ_THE119_PP_LIBRARY_SHA256:-58ac460753344844309cc3ff6bf408d049e49c944a4be28d03a9797e3edcb364}"
auau_lib_sha="${RJ_THE119_AUAU_LIBRARY_SHA256:-42f1a860b4f76ace51baa79e0c1413ac3a8fa42664e868c831d5044079e53019}"
pp_model="/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/the116_models/the116_pp_matched_basev3e_15to35_20260720/models/bdt_ppg12_basev3e_15to35/pp_tight_bdt_ppg12_base_v3E_bdt_noIso_tmva.root"
pp_ref="/sphenix/user/shuhangli/ppg12/FunWithxgboost/binned_models/model_base_v3E_split_single_tmva.root"
auau_model="/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/the111_models/the111_combined_corrected_shower_ppg12_labels_20260719_1618/combined/auau_tight_bdt_centAsFeatBase3x3_pt15to35_tmva.root"
schema_sha="$(sha256sum src/RJReplayFoundationV1.h | awk '{print $1}')"
semantic_sha="97402b8d1e51e11082015ffbd920a346d19fc3c43ae189bf6367df49939f5c6a"
photon_capture_et_min="${RJ_THE119_PHOTON_CAPTURE_ET_MIN:-5.0}"
canary_nevents="${RJ_THE119_NEVENTS:-3000}"
replay_trace="${RJ_THE119_REPLAY_TRACE:-0}"
code_commit="${RJ_THE119_CODE_COMMIT:-$(git rev-parse HEAD)}"
code_sha="${RJ_THE119_CODE_SHA256:-$(
  sha256sum \
    src/RecoilJets.cc \
    src/RecoilJets.h \
    src_AuAu/RecoilJets_AuAu.cc \
    src_AuAu/RecoilJets_AuAu.h \
    src/RJReplayFoundationV1.h \
    src/RJReplayRuntimeV1.h \
    macros/Fun4All_recoilJets_unified_impl.C \
  | sha256sum | awk '{print $1}'
)}"
only_keys="${RJ_THE119_ONLY_KEYS:-}"

export RJ_CODEX_CHAT_NAME="THE-114+THE-119 | pp/AuAu Replay Foundation"
export RJ_CODEX_THREAD_ID="019f80b5-dc56-7330-9ee7-56ef417547dc"

say(){ printf '[THE119] %s\n' "$*"; }
die(){ printf '[THE119][ERROR] %s\n' "$*" >&2; exit 2; }
sha(){ printf '%s' "$1" | sha256sum | awk '{print $1}'; }

want_row(){
  local arm="$1" lane="$2" sample="$3"
  [[ -z "$only_keys" ]] && return 0
  [[ ",${only_keys}," == *",${arm}:${lane}:${sample},"* ]]
}

require_inputs(){
  [[ -x ./RecoilJets_Condor_submit.sh ]] || die "missing submitter"
  for f in "$pp_cfg" "$auau_cfg" "$pp_lib" "$auau_lib" "$pp_model" "$pp_ref" "$auau_model"; do
    [[ -s "$f" ]] || die "missing required input: $f"
  done
  [[ "$pp_lib_sha" =~ ^[0-9a-f]{64}$ ]] || die "expected p+p library identity must be a 64-character SHA-256"
  [[ "$auau_lib_sha" =~ ^[0-9a-f]{64}$ ]] || die "expected Au+Au library identity must be a 64-character SHA-256"
  [[ "$(sha256sum "$pp_lib" | awk '{print $1}')" == "$pp_lib_sha" ]] || die "p+p library hash drift"
  [[ "$(sha256sum "$auau_lib" | awk '{print $1}')" == "$auau_lib_sha" ]] || die "Au+Au library hash drift"
  [[ "$(sha256sum "$pp_model" | awk '{print $1}')" == "228d4cb73f7dc945a613c5a604add71a372b7540c2dc8c630b533d215bb17b30" ]] || die "THE-116 model hash drift"
  [[ "$(sha256sum "$pp_ref" | awk '{print $1}')" == "7679e634260402fb3815b2733767182690eec7587f9e09bffc307a05d00d59df" ]] || die "PPG12 model hash drift"
  [[ "$(sha256sum "$auau_model" | awk '{print $1}')" == "d50c69ec98558cb80730ab45fe6801d4accbcf2221c482af91e8899cede1c925" ]] || die "THE-111 model hash drift"
  [[ "$code_commit" =~ ^[0-9a-f]{40}$ ]] || die "campaign code commit must be a 40-character Git object id"
  [[ "$code_sha" =~ ^[0-9a-f]{64}$ ]] || die "replay code identity must be a 64-character SHA-256"
  [[ "$photon_capture_et_min" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]] || die "photon capture threshold must be numeric"
  awk -v x="$photon_capture_et_min" 'BEGIN { exit !(x >= 0.0 && x < 15.0) }' || die "photon capture threshold must satisfy 0 <= ET < 15 GeV"
  [[ "$canary_nevents" =~ ^[1-9][0-9]*$ ]] || die "canary event count must be a positive integer"
  [[ "$replay_trace" == 0 || "$replay_trace" == 1 ]] || die "RJ_THE119_REPLAY_TRACE must be 0 or 1"
}

common_extra(){
  local lane="$1" dataset="$2" sample="$3" arm="$4" cfg="$5" model_sha="$6"
  local enabled=0
  [[ "$arm" == writer ]] && enabled=1
  local source_sha
  source_sha="$(sha "${tag}|${lane}|${dataset}|${sample}|accepted-canary-source-v1")"
  printf '%s' "RJ_REPLAY_FOUNDATION_V1=${enabled};RJ_REPLAY_FOUNDATION_CANARY=1;RJ_REPLAY_TRACE=${replay_trace};RJ_REPLAY_PHOTON_CAPTURE_ET_MIN=${photon_capture_et_min};RJ_REPLAY_LANE=${lane};RJ_REPLAY_DATASET=${dataset};RJ_REPLAY_SAMPLE=${sample};RJ_REPLAY_SOURCE_MANIFEST_SHA256=${source_sha};RJ_REPLAY_SCHEMA_SHA256=${schema_sha};RJ_REPLAY_SEMANTIC_SHA256=${semantic_sha};RJ_REPLAY_SOURCE_SHA256=${source_sha};RJ_REPLAY_MODEL_SHA256=${model_sha};RJ_REPLAY_CONFIG_SHA256=$(sha256sum "$cfg" | awk '{print $1}');RJ_REPLAY_CODE_SHA256=${code_sha}"
}

pp_extra(){
  local lane="$1" dataset="$2" sample="$3" arm="$4"
  printf '%s;%s' "$(common_extra "$lane" "$dataset" "$sample" "$arm" "$pp_cfg" 228d4cb73f7dc945a613c5a604add71a372b7540c2dc8c630b533d215bb17b30)" \
    "RJ_REPLAY_MODEL_SCORE_NAME=tight_bdt_score;RJ_REPLAY_REFERENCE_MODEL_FILE=${pp_ref};RJ_REPLAY_REFERENCE_MODEL_SHA256=7679e634260402fb3815b2733767182690eec7587f9e09bffc307a05d00d59df;RJ_REPLAY_REFERENCE_SCORE_NAME=ppg12_reference_bdt_score;RJ_REPLAY_WP70_BINS=0.79682856798172,0.766527533531189,0.764809787273407,0.7529897093772888,0.7708977460861206,0.8068315982818604,0.8892104029655457,0.972591757774353;RJ_REPLAY_WP80_BINS=0.7195994257926941,0.682415783405304,0.6793394684791565,0.6720289587974548,0.6960929036140442,0.7344872951507568,0.8287723064422607,0.9486955404281616;RJ_REPLAY_WP90_BINS=0.5593066215515137,0.5007686018943787,0.5124438405036926,0.5229008793830872,0.5534335374832153,0.5945547223091125,0.6987603902816772,0.8794801831245422"
}

auau_extra(){
  local lane="$1" dataset="$2" sample="$3" arm="$4"
  printf '%s;%s' "$(common_extra "$lane" "$dataset" "$sample" "$arm" "$auau_cfg" d50c69ec98558cb80730ab45fe6801d4accbcf2221c482af91e8899cede1c925)" \
    "RJ_REPLAY_MODEL_SCORE_NAME=auau_tight_bdt_score;RJ_REPLAY_WP70_INTERCEPT=0.6529177794;RJ_REPLAY_WP70_SLOPE=0.0013378442;RJ_REPLAY_WP80_INTERCEPT=0.5544148693;RJ_REPLAY_WP80_SLOPE=0.0015499421;RJ_REPLAY_WP90_INTERCEPT=0.4046618113;RJ_REPLAY_WP90_SLOPE=0.0014896756"
}

assert_fresh(){
  if condor_q "${USER:-patsfan753}" -af Args 2>/dev/null | grep -F "$base" >/dev/null; then die "active duplicate targets $base"; fi
  [[ ! -e "$base" ]] || die "output already exists: $base"
}

submit_pp(){
  local lane="$1" dataset="$2" sample="$3" arm="$4"
  if ! want_row "$arm" "$lane" "$sample"; then
    say "SKIP ${arm}:${lane}:${sample} (not in RJ_THE119_ONLY_KEYS)"
    return 0
  fi
  local out="$base/$arm/$lane/$sample"
  env RJ_CONFIG_YAML="$pp_cfg" RJ_PP_LIBRARY_OVERRIDE="$pp_lib" RJ_AUTO_MERGE=0 \
    RJ_REQUEST_MEMORY=8000MB \
    RJ_REPLAY_FOUNDATION_CANARY=1 RJ_REPLAY_LANE="$lane" RJ_REPLAY_SCHEMA_SHA256="$schema_sha" \
    RJ_REQUIRE_NON_TINY_OUTPUT=1 RJ_MIN_OUTPUT_BYTES=50000 RJ_PROFILE_JOB=1 \
    RJ_JOB_HEARTBEAT_SECONDS=120 RJ_SMOKE_OUTPUT_BASE="$out" RJ_SMOKE_SIM_NEVENTS="$canary_nevents" \
    RJ_SMOKE_DATA_RUNS=1 RJ_SMOKE_DATA_MAX_JOBS=1 RJ_SMOKE_DATA_NEVENTS="$canary_nevents" RJ_SUBMIT_EXTRA_ENV="$(pp_extra "$lane" "$dataset" "$sample" "$arm")" \
    ./RecoilJets_Condor_submit.sh "$dataset" $([[ "$dataset" == isPP ]] && printf 'condor smokeTest groupSize 1' || printf 'condorDoAllSmoke groupSize 1 maxJobs 1 SAMPLE=%s' "$sample")
}

submit_auau(){
  local lane="$1" dataset="$2" sample="$3" arm="$4"
  if ! want_row "$arm" "$lane" "$sample"; then
    say "SKIP ${arm}:${lane}:${sample} (not in RJ_THE119_ONLY_KEYS)"
    return 0
  fi
  local out="$base/$arm/$lane/$sample"
  env RJ_CONFIG_YAML="$auau_cfg" RJ_AUAU_LIBRARY_OVERRIDE="$auau_lib" RJ_AUTO_MERGE=0 \
    RJ_REQUEST_MEMORY=8000MB \
    RJ_REPLAY_FOUNDATION_CANARY=1 RJ_REPLAY_LANE="$lane" RJ_REPLAY_SCHEMA_SHA256="$schema_sha" \
    RJ_REQUIRE_NON_TINY_OUTPUT=1 RJ_MIN_OUTPUT_BYTES=50000 RJ_PROFILE_JOB=1 \
    RJ_JOB_HEARTBEAT_SECONDS=120 RJ_SMOKE_OUTPUT_BASE="$out" RJ_SMOKE_SIM_NEVENTS="$canary_nevents" \
    RJ_SMOKE_DATA_RUNS=1 RJ_SMOKE_DATA_MAX_JOBS=1 RJ_SMOKE_DATA_NEVENTS="$canary_nevents" RJ_INTERNAL_FIXED_ISO_GEV_AUAU=4.0 \
    RJ_AUAU_BUILD_TOPOCLUSTER_ISOLATION=0 RJ_AUAU_USE_TOPOCLUSTER_ISOLATION=0 \
    RJ_SUBMIT_EXTRA_ENV="$(auau_extra "$lane" "$dataset" "$sample" "$arm")" \
    ./RecoilJets_Condor_submit.sh "$dataset" $([[ "$dataset" == isAuAu ]] && printf 'condor smokeTest groupSize 1' || printf 'condorDoAllSmoke groupSize 1 maxJobs 1 SAMPLE=%s' "$sample")
}

preflight(){
  require_inputs
  bash -n "$0" scripts/sdcc/runtime/condor/RecoilJets_Condor.sh scripts/sdcc/runtime/condor/RecoilJets_Condor_AuAu.sh
  mkdir -p "$evidence"
  {
    printf 'tag=%s\nbase=%s\ncode_commit=%s\ncode_sha256=%s\nschema_sha=%s\nsemantic_sha=%s\nphoton_capture_et_min_gev=%s\ncanary_nevents=%s\nreplay_trace=%s\nonly_keys=%s\n' \
      "$tag" "$base" "$code_commit" "$code_sha" "$schema_sha" "$semantic_sha" "$photon_capture_et_min" "$canary_nevents" "$replay_trace" "$only_keys"
    sha256sum "$pp_cfg" "$auau_cfg" "$pp_lib" "$auau_lib" "$pp_model" "$pp_ref" "$auau_model"
  } > "$evidence/preflight_receipt.txt"
  say "PREFLIGHT_PASS evidence=$evidence/preflight_receipt.txt"
}

submit(){
  preflight
  assert_fresh
  mkdir -p "$evidence"
  exec > >(tee -a "$evidence/submission.log") 2>&1
  for arm in direct writer; do
    submit_pp pp_data isPP pp_data "$arm"
    submit_pp pp_photon_sim isSim run28_photonjet5 "$arm"
    submit_pp pp_photon_sim isSim run28_photonjet20 "$arm"
    submit_pp pp_inclusive_sim isSimInclusive run28_jet8 "$arm"
    submit_pp pp_inclusive_sim isSimInclusive run28_jet40 "$arm"
    submit_auau auau_data isAuAu auau_data "$arm"
    submit_auau auau_photon_embedded isSimEmbedded run28_embeddedPhoton12 "$arm"
    submit_auau auau_photon_embedded isSimEmbedded run28_embeddedPhoton20 "$arm"
    submit_auau auau_inclusive_embedded isSimEmbeddedInclusive run28_embeddedJet12 "$arm"
    submit_auau auau_inclusive_embedded isSimEmbeddedInclusive run28_embeddedJet40 "$arm"
  done
  condor_q "${USER:-patsfan753}" -af ClusterId ProcId JobStatus Args | grep -F "$base" | tee "$evidence/initial_queue.tsv" || true
}

status(){
  condor_q "${USER:-patsfan753}" -af ClusterId ProcId JobStatus HoldReason Args 2>/dev/null | grep -F "$base" || true
  find "$base" -type f -name '*.root' -printf '%s\t%p\n' 2>/dev/null | sort -k2 || true
}

case "$mode" in
  preflight) preflight ;;
  submit) submit ;;
  status) status ;;
  *) die "usage: $0 preflight|submit|status" ;;
esac
