#!/usr/bin/env bash
# Submit PhotonJet12/20 generator xsec scans for stitching-boundary tests.
#
# This is generator-only. It runs estimateEmbeddedPhotonXsec.sh with matching
# PhotonJet stitching windows so downstream RecoilJets stitching can use a
# self-consistent ownership window and sigma_eff pair.
#
# Variants:
#   exclusive21:
#     PhotonJet12 owns 12 <= pT_gamma < 21, PhotonJet20 owns pT_gamma >= 21.
#   exclusive22:
#     PhotonJet12 owns 12 <= pT_gamma < 22, PhotonJet20 owns pT_gamma >= 22.
#   lowerOnly21:
#     PhotonJet12 uses pT_gamma >= 12, PhotonJet20 uses pT_gamma >= 21.
#   lowerOnly22:
#     PhotonJet12 uses pT_gamma >= 12, PhotonJet20 uses pT_gamma >= 22.
#     lowerOnly variants are diagnostics only; threshold samples overlap unless
#     a separate ownership, averaging, or subtraction prescription is defined.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

BOUNDARIES="${BOUNDARIES:-21 22}"
VARIANTS="${VARIANTS:-}"
XSEC_SHARDS="${XSEC_SHARDS:-50}"
XSEC_EVENTS="${XSEC_EVENTS:-1000000}"
XSEC_MEMORY="${XSEC_MEMORY:-900MB}"
XSEC_MODE="${XSEC_MODE:-interpreted}"
XSEC_BASE_SEED="${XSEC_BASE_SEED:-753210}"
MANIFEST_LIST="${BASE_DIR}/photon_stitch_xsec_manifests_$(date '+%Y%m%d_%H%M%S').txt"

if ! command -v condor_submit >/dev/null 2>&1; then
  echo "[ERROR] condor_submit not found. Run this on an SDCC submit host." >&2
  exit 1
fi

cd "${BASE_DIR}"

if [[ -z "${VARIANTS}" ]]; then
  VARIANTS=""
  for boundary in ${BOUNDARIES}; do
    VARIANTS="${VARIANTS} exclusive${boundary}"
  done
  VARIANTS="${VARIANTS# }"
fi

echo "PhotonJet12/20 shifted-boundary xsec scan"
echo "Base       : ${BASE_DIR}"
echo "Variants   : ${VARIANTS}"
echo "Shards     : ${XSEC_SHARDS} per sample"
echo "Events     : ${XSEC_EVENTS} raw attempts per shard"
echo "Memory     : ${XSEC_MEMORY}"
echo "Mode       : ${XSEC_MODE}"
echo "Manifest list: ${MANIFEST_LIST}"
echo

: > "${MANIFEST_LIST}"

idx=0
for variant in ${VARIANTS}; do
  unset RJ_XSEC_PHOTONJET12_MIN RJ_XSEC_PHOTONJET12_MAX
  unset RJ_XSEC_PHOTONJET20_MIN RJ_XSEC_PHOTONJET20_MAX

  case "${variant}" in
    exclusive*)
      boundary="${variant#exclusive}"
      if ! [[ "${boundary}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
        echo "[ERROR] Bad boundary in variant '${variant}'" >&2
        exit 2
      fi
      export RJ_XSEC_PHOTONJET12_MIN=12
      export RJ_XSEC_PHOTONJET12_MAX="${boundary}"
      export RJ_XSEC_PHOTONJET20_MIN="${boundary}"
      export RJ_XSEC_PHOTONJET20_MAX=-1
      description="exclusive: PhotonJet12 12-${boundary}, PhotonJet20 >=${boundary}"
      ;;
    lowerOnly*)
      boundary="${variant#lowerOnly}"
      if ! [[ "${boundary}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
        echo "[ERROR] Bad boundary in variant '${variant}'" >&2
        exit 2
      fi
      export RJ_XSEC_PHOTONJET12_MIN=12
      export RJ_XSEC_PHOTONJET12_MAX=-1
      export RJ_XSEC_PHOTONJET20_MIN="${boundary}"
      export RJ_XSEC_PHOTONJET20_MAX=-1
      description="diagnostic threshold-only: PhotonJet12 >=12, PhotonJet20 >=${boundary}"
      ;;
    *)
      echo "[ERROR] Unknown variant '${variant}'. Supported: exclusive<N> lowerOnly<N>" >&2
      exit 2
      ;;
  esac

  if ! [[ "${boundary}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    echo "[ERROR] Bad boundary '${boundary}'" >&2
    exit 2
  fi

  seed=$(( XSEC_BASE_SEED + 1000 * idx ))
  before="$(ls -t "${BASE_DIR}"/pythia_xsec_firstPass_*.txt 2>/dev/null | head -n 1 || true)"

  echo "====================================================================="
  echo "[SUBMIT] ${variant}"
  echo "  ${description}"
  echo "  xsec seed: ${seed}"
  echo "====================================================================="

  ./scripts/estimateEmbeddedPhotonXsec.sh firstPass \
    --family photon \
    --xsec-shards "${XSEC_SHARDS}" \
    --xsec-events "${XSEC_EVENTS}" \
    --xsec-memory "${XSEC_MEMORY}" \
    --xsec-base-seed "${seed}" \
    --mode "${XSEC_MODE}"

  after="$(ls -t "${BASE_DIR}"/pythia_xsec_firstPass_*.txt 2>/dev/null | head -n 1 || true)"
  if [[ -z "${after}" || "${after}" == "${before}" ]]; then
    echo "[ERROR] Could not identify new manifest after variant ${variant}" >&2
    exit 3
  fi

  printf 'variant=%s manifest=%s\n' "${variant}" "${after}" | tee -a "${MANIFEST_LIST}"
  idx=$(( idx + 1 ))
done

echo
echo "[DONE] Submitted PhotonJet xsec scans."
echo "[NEXT] Check Condor, then aggregate each manifest with:"
awk '{
  manifest = $2
  sub(/^manifest=/, "", manifest)
  print "  ./scripts/estimateEmbeddedPhotonXsec.sh secondPass --manifest " manifest
}' "${MANIFEST_LIST}"
