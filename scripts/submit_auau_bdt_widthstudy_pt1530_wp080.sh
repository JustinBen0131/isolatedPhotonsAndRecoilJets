#!/usr/bin/env bash
set -euo pipefail

campaign_tag="${RJ_WIDTHSTUDY_CAMPAIGN_TAG:-widthstudy_pt1530_wp080_$(date +%Y%m%d_%H%M%S)}"
yaml="${RJ_WIDTHSTUDY_CONFIG_YAML:-macros/analysis_config_auau_bdt_widthstudy_pt1530_wp080.yaml}"
group_size="${RJ_WIDTHSTUDY_GROUP_SIZE:-7}"
notify="${RJ_NOTIFY_EMAILS:-just0131@gmail.com}"

submit_one() {
  local dataset="$1"
  local bulk_var="$2"
  local bulk_base="$3"
  local merge_base="$4"

  echo
  echo "====================================================================="
  echo "Submitting ${dataset} width-study validation"
  echo "  tag        : ${campaign_tag}"
  echo "  yaml       : ${yaml}"
  echo "  bulk base  : ${bulk_base}"
  echo "  merge base : ${merge_base}"
  echo "  groupSize  : ${group_size}"
  echo "====================================================================="

  env \
    RJ_NOTIFY_EMAILS="${notify}" \
    RJ_PROFILE_JOB=1 \
    RJ_CONFIG_YAML="${yaml}" \
    RJ_REQUEST_MEMORY="${RJ_REQUEST_MEMORY:-6000MB}" \
    RJ_AUTO_MEMORY_RETRY_CAP_MB="${RJ_AUTO_MEMORY_RETRY_CAP_MB:-10000}" \
    RJ_SIM_FIRSTROUND_REQUEST_MEMORY="${RJ_SIM_FIRSTROUND_REQUEST_MEMORY:-6000MB}" \
    RJ_SIM_MERGE_GROUP_SIZE="${RJ_SIM_MERGE_GROUP_SIZE:-75}" \
    RJ_MERGE_OUT_BASE_OVERRIDE="${merge_base}" \
    "${bulk_var}=${bulk_base}" \
    ./RecoilJets_Condor_submit.sh "${dataset}" condorDoAll groupSize "${group_size}"
}

submit_one \
  isSimEmbedded \
  RJ_SIMEMBED_DEST_BASE \
  "/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded_${campaign_tag}" \
  "/sphenix/u/patsfan753/scratch/thesisAnalysis/output_${campaign_tag}"

submit_one \
  isSimEmbeddedInclusive \
  RJ_SIMEMBEDINCLUSIVE_DEST_BASE \
  "/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembeddedinclusive_${campaign_tag}" \
  "/sphenix/u/patsfan753/scratch/thesisAnalysis/output_${campaign_tag}"

echo
echo "TRACK_THIS_CAMPAIGN=${campaign_tag}"
