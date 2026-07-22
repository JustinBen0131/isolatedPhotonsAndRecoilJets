#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd -P)"
cd "$repo_root"

mode="${1:-preflight}"
tag="${RJ_THE134_TAG:-the134_h70_shower_factorial_replay_canary_20260722_v1}"
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
only_keys="${RJ_THE134_ONLY_KEYS:-$default_only_keys}"

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

case "$mode" in
  preflight|submit|status)
    exec scripts/sdcc/workflows/diagnostics/submit_the119_replay_foundation_canaries.sh "$mode"
    ;;
  *)
    echo "usage: $0 preflight|submit|status" >&2
    exit 2
    ;;
esac
