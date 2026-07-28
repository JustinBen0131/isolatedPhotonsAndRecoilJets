#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd -P)"
cd "$repo_root"

requested_mode="${1:-preflight}"
sidecar_only_mode=0
mode="$requested_mode"
case "$requested_mode" in
  sidecar-preflight|sidecar-submit|sidecar-status|sidecar-validate)
    sidecar_only_mode=1
    mode="${requested_mode#sidecar-}"
    ;;
esac
if (( sidecar_only_mode )); then
  tag="${RJ_THE134_TAG:-the134_h70_sidecar_differential_canary_20260728_v1}"
else
  tag="${RJ_THE134_TAG:-the134_h70_shower_factorial_replay_canary_20260722_v1}"
fi
build_root="${RJ_THE134_BUILD_ROOT:-/sphenix/u/patsfan753/scratch/thesisAnalysis/.recoiljets_tmp/the134_h70_schema_v2_20260722_r3}"
pp_library="${RJ_THE134_PP_LIBRARY:-${build_root}/install/pp/lib/libRecoilJets.so}"
auau_library="${RJ_THE134_AUAU_LIBRARY:-${build_root}/install/auau/lib/libRecoilJetsAuAu.so}"
: "${RJ_THE134_PP_LIBRARY_SHA256:?set the frozen p+p library SHA-256}"
: "${RJ_THE134_AUAU_LIBRARY_SHA256:?set the frozen Au+Au library SHA-256}"

factorial_semantic_sha="$({
  sha256sum \
    src/RJReplayFoundationV1.h \
    src/RJReplayRuntimeV1.h \
    src/RJShowerFactorialV1.h
  printf '%s\n' 'THE134|H70,H0,G70,G0,O70,O0,R70|nominal-domain=15<=ET<35|schema=2'
} | sha256sum | awk '{print $1}')"
if [[ ! "$factorial_semantic_sha" =~ ^[0-9a-f]{64}$ ]]; then
  echo "[THE134][ERROR] could not derive the factorial semantic identity" >&2
  exit 2
fi

# One source witness per detector lane, each with a direct and writer arm.  The
# completed THE-119 matrix remains the base replay-infrastructure certificate;
# this bounded matrix certifies only the new multi-view shower payload.
if (( sidecar_only_mode )); then
  default_only_keys="$({
    for arm in direct writer; do
      printf '%s\n' \
        "${arm}:pp_inclusive_sim:run28_jet8" \
        "${arm}:auau_inclusive_embedded:run28_embeddedJet12"
    done
  } | paste -sd, -)"
else
  default_only_keys="$({
    for arm in direct writer; do
      printf '%s\n' \
        "${arm}:pp_data:pp_data" \
        "${arm}:pp_photon_sim:run28_photonjet20" \
        "${arm}:pp_inclusive_sim:run28_jet8" \
        "${arm}:auau_data:auau_data" \
        "${arm}:auau_photon_embedded:run28_embeddedPhoton12" \
        "${arm}:auau_inclusive_embedded:run28_embeddedJet12"
    done
  } | paste -sd, -)"
fi
only_keys="${RJ_THE134_ONLY_KEYS:-$default_only_keys}"
if (( sidecar_only_mode )) && [[ "$only_keys" != "$default_only_keys" ]]; then
  echo "[THE134][ERROR] sidecar differential mode owns exactly the frozen four rows" >&2
  exit 2
fi

export RJ_CODEX_CHAT_NAME="THE-114+THE-134 | pp/AuAu H70 Production Gate"
export RJ_CODEX_THREAD_ID="019f80b5-dc56-7330-9ee7-56ef417547dc"
export RJ_THE119_TAG="$tag"
export RJ_THE119_OUTPUT_ROOT="${RJ_THE134_OUTPUT_ROOT:-/sphenix/tg/tg01/bulk/jbennett/thesisAna/recoiljets/smoke/replay_foundation/${tag}}"
export RJ_THE119_EVIDENCE_ROOT="${RJ_THE134_EVIDENCE_ROOT:-/sphenix/u/patsfan753/scratch/thesisAnalysis/evidence/qa/${tag}}"
export RJ_THE119_PP_LIBRARY="$pp_library"
export RJ_THE119_AUAU_LIBRARY="$auau_library"
export RJ_THE119_PP_LIBRARY_SHA256="$RJ_THE134_PP_LIBRARY_SHA256"
export RJ_THE119_AUAU_LIBRARY_SHA256="$RJ_THE134_AUAU_LIBRARY_SHA256"
export RJ_THE119_PP_MODEL_SHOWER_DEFINITION=H70
export RJ_THE119_PP_REFERENCE_MODEL_SHOWER_DEFINITION=H70
export RJ_THE119_AUAU_MODEL_SHOWER_DEFINITION=H0
export RJ_THE119_SEMANTIC_SHA256="$factorial_semantic_sha"
export RJ_THE119_ONLY_KEYS="$only_keys"
export RJ_THE119_NEVENTS="${RJ_THE134_NEVENTS:-3000}"
export RJ_THE119_PHOTON_CAPTURE_ET_MIN=5.0
export RJ_THE119_JET_CONSTITUENT_PT_MIN=5.0
export RJ_THE119_PP_CAPTURE_WITNESS_QA=1
export RJ_THE119_CODE_COMMIT="${RJ_THE134_CODE_COMMIT:-$(git rev-parse HEAD)}"
if (( sidecar_only_mode )); then
  export RJ_THE119_EXTRA_ENV_TEMPLATE=''
  export RJ_THE119_PP_WRITER_EXTRA_ENV_TEMPLATE=''
  export RJ_THE119_AUAU_WRITER_EXTRA_ENV_TEMPLATE=''
  export RJ_THE119_WRITER_EXTRA_ENV_TEMPLATE='RJ_THE134_MULTIVIEW_TRAINING_V1=1;RJ_THE134_MULTIVIEW_TRAINING_FILE=__OUTPUT__/RJPhotonTrainingViewV1.root;RJ_THE134_MULTIVIEW_SIDECAR_ONLY_V1=1;RJ_THE134_EXPECTED_SOURCE_ROLE=background'
  export RJ_THE119_PP_EXTRA_ENV_TEMPLATE='RJ_REPLAY_PERIOD=0mrad;RJ_REPLAY_SI_DI_ROLE=SI;RJ_REPLAY_OWNERSHIP_STATE=source_role_frozen;RJ_SIM_ALLOW_NONE_LISTS=0;RJ_PPG12_PHOTON_YIELD=1;RJ_PPG12_PHOTON_YIELD_DOUBLE=0;RJ_PPG12_PERIOD=0mrad;RJ_PPG12_PERIOD_USE_LUMI_WEIGHT=1;RJ_PPG12_PERIOD_STRICT_DI=0;RJ_PPG12_PERIOD_ALLOW_ALL_SIM=0;RJ_PPG12_PERIOD_ALLOW_MIX_OVERRIDE=0;RJ_PPG12_PERIOD_ALLOW_VERTEX_FILE_OVERRIDE=0;RJ_PP_PHOTONID_EXTRACT_ONLY=1;RJ_PP_PHOTONID_TRAINING_TREE=1;RJ_PP_PHOTONID_TRAINING_TREE_MAX_ENTRIES=0;RJ_PP_PHOTONID_SOURCE_ROLE=background;RJ_PP_PHOTONID_PPG12_FILTER=1;RJ_PP_PHOTONID_REQUIRE_PRESELECTION=0'
  export RJ_THE119_AUAU_EXTRA_ENV_TEMPLATE='RJ_REPLAY_PERIOD=AUAU_RUN24;RJ_REPLAY_SI_DI_ROLE=EMBEDDED;RJ_REPLAY_OWNERSHIP_STATE=source_role_frozen;RJ_SIM_ALLOW_NONE_LISTS=1;RJ_AUAU_BDT_EXTRACT_ONLY=1;RJ_AUAU_BDT_TRAINING_TREE=1;RJ_AUAU_BDT_TRAINING_TREE_MAX_ENTRIES=0;RJ_AUAU_BDT_NPB_DATA_TAGGING=0;RJ_REQUIRE_EMBEDDED_MINBIAS_CLASSIFIER=1;RJ_AUAU_BUILD_TOPOCLUSTER_ISOLATION=0;RJ_AUAU_USE_TOPOCLUSTER_ISOLATION=0'
fi

validate_terminal_rows() {
  local initial_queue="${RJ_THE119_EVIDENCE_ROOT}/initial_queue.tsv"
  local cluster proc cluster_proc queue_state history_state job_status exit_code
  [[ -s "$initial_queue" ]] ||
    { echo "[THE134][ERROR] missing exact initial queue receipt: $initial_queue" >&2; return 2; }
  [[ "$(wc -l < "$initial_queue" | tr -d ' ')" == 4 ]] ||
    { echo "[THE134][ERROR] sidecar differential canary must own exactly four rows" >&2; return 2; }
  while read -r cluster proc _status _args; do
    [[ "$cluster" =~ ^[0-9]+$ && "$proc" =~ ^[0-9]+$ ]] ||
      { echo "[THE134][ERROR] malformed exact Condor identity in $initial_queue" >&2; return 2; }
    cluster_proc="${cluster}.${proc}"
    queue_state="$(condor_q "$cluster_proc" -af JobStatus HoldReason NumJobStarts 2>/dev/null || true)"
    [[ -z "$queue_state" ]] ||
      { echo "[THE134][ERROR] row is not terminal: $cluster_proc $queue_state" >&2; return 2; }
    history_state="$(condor_history "$cluster_proc" -limit 1 -af JobStatus ExitCode 2>/dev/null || true)"
    read -r job_status exit_code <<< "$history_state"
    [[ "$job_status" == 4 && "$exit_code" == 0 ]] ||
      { echo "[THE134][ERROR] terminal gate failed: $cluster_proc ${history_state:-NOT_FOUND}" >&2; return 2; }
  done < "$initial_queue"
}

case "$mode" in
  preflight|submit|status)
    exec scripts/sdcc/workflows/diagnostics/submit_the119_replay_foundation_canaries.sh "$mode"
    ;;
  validate)
    (( sidecar_only_mode )) ||
      { echo "[THE134][ERROR] validate is reserved for sidecar-* modes" >&2; exit 2; }
    validate_terminal_rows
    exec python3 scripts/diagnostics/replay_foundation/validate_the134_sidecar_differential.py \
      --output-root "$RJ_THE119_OUTPUT_ROOT" \
      --preflight-receipt "$RJ_THE119_EVIDENCE_ROOT/preflight_receipt.txt" \
      --output-json "$RJ_THE119_EVIDENCE_ROOT/sidecar_differential_certificate.json"
    ;;
  *)
    echo "usage: $0 preflight|submit|status|sidecar-preflight|sidecar-submit|sidecar-status|sidecar-validate" >&2
    exit 2
    ;;
esac
