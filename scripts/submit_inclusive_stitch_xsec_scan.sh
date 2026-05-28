#!/usr/bin/env bash
# Submit embedded-inclusive Jet12/20/30 generator xsec scans for stitching tests.
#
# Variants:
#   exclusive31:
#     Jet12 owns 12 <= pT < 21, Jet20 owns 21 <= pT < 31, Jet30 owns pT >= 31.
#   lowerOnly31:
#     Jet12 uses pT >= 12, Jet20 uses pT >= 21, Jet30 uses pT >= 31.
#     This is a diagnostic only; summing threshold samples double-counts
#     overlapping phase space unless a separate averaging/ownership rule is used.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

VARIANTS="${VARIANTS:-exclusive31 lowerOnly31}"
XSEC_SHARDS="${XSEC_SHARDS:-50}"
XSEC_EVENTS="${XSEC_EVENTS:-1000000}"
XSEC_MEMORY="${XSEC_MEMORY:-900MB}"
XSEC_MODE="${XSEC_MODE:-interpreted}"
XSEC_BASE_SEED="${XSEC_BASE_SEED:-753310}"
MANIFEST_LIST="${BASE_DIR}/inclusive_stitch_xsec_manifests_$(date '+%Y%m%d_%H%M%S').txt"

if ! command -v condor_submit >/dev/null 2>&1; then
  echo "[ERROR] condor_submit not found. Run this on an SDCC submit host." >&2
  exit 1
fi

cd "${BASE_DIR}"

echo "Embedded-inclusive Jet12/20/30 xsec scan"
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
  unset RJ_XSEC_EMBEDDED_JET12_MIN RJ_XSEC_EMBEDDED_JET12_MAX
  unset RJ_XSEC_EMBEDDED_JET20TO30_MIN RJ_XSEC_EMBEDDED_JET20TO30_MAX
  unset RJ_XSEC_EMBEDDED_JET30_MIN RJ_XSEC_EMBEDDED_JET30_MAX

  case "${variant}" in
    exclusive31)
      export RJ_XSEC_EMBEDDED_JET12_MIN=12
      export RJ_XSEC_EMBEDDED_JET12_MAX=21
      export RJ_XSEC_EMBEDDED_JET20TO30_MIN=21
      export RJ_XSEC_EMBEDDED_JET20TO30_MAX=31
      export RJ_XSEC_EMBEDDED_JET30_MIN=31
      unset RJ_XSEC_EMBEDDED_JET30_MAX
      description="exclusive: Jet12 12-21, Jet20 21-31, Jet30 >=31"
      ;;
    lowerOnly31)
      export RJ_XSEC_EMBEDDED_JET12_MIN=12
      export RJ_XSEC_EMBEDDED_JET12_MAX=-1
      export RJ_XSEC_EMBEDDED_JET20TO30_MIN=21
      export RJ_XSEC_EMBEDDED_JET20TO30_MAX=-1
      export RJ_XSEC_EMBEDDED_JET30_MIN=31
      export RJ_XSEC_EMBEDDED_JET30_MAX=-1
      description="diagnostic threshold-only: Jet12 >=12, Jet20 >=21, Jet30 >=31"
      ;;
    *)
      echo "[ERROR] Unknown variant '${variant}'. Supported: exclusive31 lowerOnly31" >&2
      exit 2
      ;;
  esac

  seed=$(( XSEC_BASE_SEED + 1000 * idx ))
  before="$(ls -t "${BASE_DIR}"/pythia_xsec_firstPass_*.txt 2>/dev/null | head -n 1 || true)"

  echo "====================================================================="
  echo "[SUBMIT] ${variant}"
  echo "  ${description}"
  echo "  xsec seed: ${seed}"
  echo "====================================================================="

  ./scripts/estimateEmbeddedPhotonXsec.sh firstPass \
    --family inclusive3 \
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
echo "[DONE] Submitted inclusive-stitch xsec scans."
echo "[NEXT] Check Condor, then aggregate each manifest with:"
awk '{
  manifest = $2
  sub(/^manifest=/, "", manifest)
  print "  ./scripts/estimateEmbeddedPhotonXsec.sh secondPass --manifest " manifest
}' "${MANIFEST_LIST}"
