#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
if [[ -f "${script_dir}/../../runtime/io/recoiljets_io_paths.sh" ]]; then
  source "${script_dir}/../../runtime/io/recoiljets_io_paths.sh"
elif [[ -f scripts/sdcc/runtime/io/recoiljets_io_paths.sh ]]; then
  source scripts/sdcc/runtime/io/recoiljets_io_paths.sh
else
  echo "[THE98][ERROR] Cannot locate recoiljets_io_paths.sh" >&2
  exit 2
fi

repo_root="$(rj_find_repo_root "$script_dir")"
cd "$repo_root"

say() { printf '[THE98] %s\n' "$*"; }
die() { printf '[THE98][ERROR] %s\n' "$1" >&2; exit "${2:-2}"; }

mode="${1:-print}"
campaign_tag="${RJ_THE98_CAMPAIGN_TAG:-the98_auau_emcal_phosub_variantB_$(date +%Y%m%d_%H%M)}"
family="${RJ_THE98_FAMILY:-auau_emcal_phosub_iso}"
config_dir="${RJ_THE98_CONFIG_DIR:-${repo_root}/condor_generated_configs/${campaign_tag}}"
yaml="${config_dir}/analysis_config_${campaign_tag}_nominal_15to35.yaml"
manifest="${config_dir}/the98_campaign_manifest.tsv"
submit_log="${config_dir}/submit_${campaign_tag}.log"

sim_group="${RJ_THE98_SIM_GROUP_SIZE:-7}"
sim_memory="${RJ_THE98_SIM_MEMORY:-12000MB}"
merge_group="${RJ_THE98_MERGE_GROUP_SIZE:-75}"
merge_memory="${RJ_THE98_MERGE_MEMORY:-16000MB}"
notify="${RJ_NOTIFY_EMAILS:-just0131@gmail.com}"
allow_existing="${RJ_THE98_ALLOW_EXISTING:-0}"

export RJ_CODEX_CHAT_NAME="${RJ_CODEX_CHAT_NAME:-THE-98_AuAu_EMCal_UE_Sub_Iso}"
export RJ_CODEX_THREAD_ID="${RJ_CODEX_THREAD_ID:-}"
[[ -n "$RJ_CODEX_THREAD_ID" ]] || die "RJ_CODEX_THREAD_ID is required for Codex-launched submission provenance."

rj_validate_campaign_tag "$campaign_tag" || exit 2
rj_validate_campaign_tag "$family" || exit 2

bulk_tag="${campaign_tag}_signal"
canary_tag="${campaign_tag}_canary"
signal_bulk="$(rj_bulk_signal_root_for_campaign "$family" "$bulk_tag")"
merge_base="$(rj_merge_root_for_campaign "$bulk_tag")"
canary_bulk="$(rj_bulk_signal_root_for_campaign "$family" "$canary_tag")"
canary_merge="$(rj_merge_root_for_campaign "$canary_tag")"

usage() {
  cat <<EOF
Usage: $0 <print|canary-dryrun|canary-submit|production-submit|status|heartbeat-prompt>

Required for submission:
  RJ_CODEX_THREAD_ID=<current Codex thread id>

Optional:
  RJ_THE98_CAMPAIGN_TAG=${campaign_tag}
  RJ_THE98_ALLOW_EXISTING=1   allow intentional resume inspection

Physics scope:
  isSimEmbedded Photon12+20 only; clusterUEpipeline=variantB; 15-35 GeV;
  5% centrality bins; R=0.3 and R=0.4 truth-matched isolation histograms.
EOF
}

write_yaml() {
  mkdir -p "$config_dir"
  cat > "$yaml" <<'EOF'
# THE-98: EMCal tower-level UE subtraction followed by PHOSUB reclustering.
photon_eta_abs_max: 0.7
jet_pt_min: [5.0]
jet_radii: [0.2, 0.3, 0.4]
back_to_back_dphi_min_pi_fraction: [0.875]
use_vz_cut: true
vz_cut_cm: [60]
setMinBiasClassifer: true
centrality_edges: [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80]
clusterUEpipeline: [variantB]

isSlidingIso: false
isSlidingAndFixed: true
fixedGeV: [4.0]
coneR: [0.40]
isolation_wp: {aGeV: 0.081, bPerGeV: 0.019, sideGapGeV: 0.0, truthIsoGeV: 4.0, towerMin: 0.0}
auau_cent_iso_wp:
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
auau_cent_iso_wp_r30:
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
auau_cent_iso_wp_r40:
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}
  - {aGeV: 4.0, bPerGeV: 0.0, sideGapGeV: 0.0}

photon_id_pre:   {e11e33_max: 0.98, et1_min: 0.60, et1_max: 1.00, e32e35_min: 0.80, e32e35_max: 1.00, weta_max: 0.60}
photon_id_tight: {w_lo: 0.0, w_hi_intercept: 0.15, w_hi_slope: 0.006, e11e33_min: 0.40, e11e33_max: 0.98, et1_min: 0.90, et1_max: 1.00, e32e35_min: 0.92, e32e35_max: 1.00}
matching: {pho_dr_max: 0.05, jet_dr_max: 0.3}

jes3_photon_pt_bins: [15, 17, 19, 21, 23, 26, 30, 35]
unfold_reco_photon_pt_bins:  [15, 17, 19, 21, 23, 26, 30, 35]
unfold_truth_photon_pt_bins: [15, 17, 19, 21, 23, 26, 30, 35]
unfold_jet_pt_binning: {start: 0.0, stop: 60.0, step: 0.5}
unfold_xj_bins: [0.0, 0.20, 0.24, 0.29, 0.35, 0.41, 0.50, 0.60, 0.72, 0.86, 1.03, 1.24, 1.49, 1.78, 2.14, 3.0]

photon_id_sets:
  - [newPPG12, reference, reference]

notify_emails: [just0131@gmail.com]
npb_model_file: /sphenix/user/shuhangli/ppg12/FunWithxgboost/npb_models/npb_score_split_tmva.root
npb_cut: 0.5
npb_features: [cluster_Et, cluster_Eta, vertexz, e11_over_e33, e32_over_e35, e11_over_e22, e11_over_e13, e11_over_e15, e11_over_e17, e11_over_e31, e11_over_e51, e11_over_e71, e22_over_e33, e22_over_e35, e22_over_e37, e22_over_e53, cluster_weta_cogx, cluster_wphi_cogx, cluster_et1, cluster_et2, cluster_et3, cluster_et4, cluster_w32, cluster_w52, cluster_w72]
tight_bdt_model_file: /sphenix/user/shuhangli/ppg12/FunWithxgboost/binned_models/model_base_v3E_split_single_tmva.root
tight_bdt_min_intercept: 0.815625
tight_bdt_min_slope: -0.0015625
tight_bdt_max: 1.0
nontight_bdt_min_intercept: 0.7333333333333333
nontight_bdt_min_slope: -0.01333333333333333
nontight_bdt_max_intercept: 0.684375
nontight_bdt_max_slope: 0.0015625
tight_bdt_features: [cluster_Et, cluster_weta_cogx, cluster_wphi_cogx, vertexz, cluster_Eta, e11_over_e33, cluster_et1, cluster_et2, cluster_et3, cluster_et4, e32_over_e35]
auau_bdt_training_tree: false
auau_bdt_training_tree_max_entries: 0
auau_bdt_npb_data_tagging: false
vertex_reweight_on_auau: false
centrality_reweight_on: false
event_display_tree: false
EOF

  printf 'dataset\tyaml\tclusterUEpipeline\tsamples\tbulk_base\tmerge_base\n' > "$manifest"
  printf 'isSimEmbedded\t%s\tvariantB\trun28_embeddedPhoton12+run28_embeddedPhoton20\t%s\t%s\n' \
    "$yaml" "$signal_bulk" "$merge_base" >> "$manifest"
}

check_unused_path() {
  local path="$1"
  local label="$2"
  if [[ -e "$path" && "$allow_existing" != "1" ]]; then
    die "${label} already exists: ${path}. Use a fresh tag or set RJ_THE98_ALLOW_EXISTING=1 for intentional resume inspection." 5
  fi
}

preflight() {
  [[ -s "$yaml" ]] || die "missing generated YAML: $yaml"
  if condor_q "${USER:-patsfan753}" -af Args 2>/dev/null | grep -F "$campaign_tag" >/dev/null; then
    die "active Condor jobs already mention ${campaign_tag}" 4
  fi
}

common_env=(
  "RJ_NOTIFY_EMAILS=${notify}"
  "RJ_PROFILE_JOB=1"
  "RJ_CONFIG_YAML=${yaml}"
  "RJ_REQUEST_MEMORY=${sim_memory}"
  "RJ_AUTO_MERGE=1"
  "RJ_SIM_FIRSTROUND_REQUEST_MEMORY=${merge_memory}"
  "RJ_SIM_MERGE_GROUP_SIZE=${merge_group}"
  "RJ_ID_FANOUT_MAX_ROWS=1"
  "RJ_PHOTON_ID_ROW_MATCH=preselectionNewPPG12_tightReference_nonTightReference"
  "RJ_INTERNAL_FIXED_ISO_GEV_AUAU=4.0"
  "RJ_AUAU_BUILD_TOPOCLUSTER_ISOLATION=0"
  "RJ_AUAU_USE_TOPOCLUSTER_ISOLATION=0"
  "RJ_HIUE_VERBOSITY=${RJ_HIUE_VERBOSITY:-1}"
  "RJ_CODEX_CHAT_NAME=${RJ_CODEX_CHAT_NAME}"
  "RJ_CODEX_THREAD_ID=${RJ_CODEX_THREAD_ID}"
)

print_contract() {
  cat <<EOF
RECOILJETS_THE98_AUAU_EMCAL_PHOSUB_ISO_V1
campaign_tag=${campaign_tag}
config_dir=${config_dir}
yaml=${yaml}
manifest=${manifest}
clusterUEpipeline=variantB
dataset=isSimEmbedded
samples=run28_embeddedPhoton12,run28_embeddedPhoton20
pt_bins=15,17,19,21,23,26,30,35
centrality_bins=0:5:80
signal_bulk=${signal_bulk}
merge_base=${merge_base}
canary_bulk=${canary_bulk}
canary_merge=${canary_merge}
RJ_CODEX_CHAT_NAME=${RJ_CODEX_CHAT_NAME}
RJ_CODEX_THREAD_ID=${RJ_CODEX_THREAD_ID}
EOF
}

write_yaml
print_contract

case "$mode" in
  print)
    preflight
    column -t -s $'\t' "$manifest" || sed -n '1,5p' "$manifest"
    ;;
  canary-dryrun|canary-submit)
    preflight
    check_unused_path "$canary_bulk" "canary bulk output"
    check_unused_path "$canary_merge" "canary merge output"
    dryrun=0
    [[ "$mode" == "canary-dryrun" ]] && dryrun=1
    env "${common_env[@]}" \
      "RJ_DAG_DRYRUN=${dryrun}" \
      "RJ_SMOKE_OUTPUT_BASE=${canary_bulk}" \
      "RJ_MERGE_OUT_BASE_OVERRIDE=${canary_merge}" \
      "RJ_SMOKE_SIM_NEVENTS=${RJ_THE98_CANARY_NEVENTS:-1000}" \
      "RJ_SMOKE_GROUPSIZE_EMBEDDED=1" \
      "RJ_SMOKE_SIM_MAX_JOBS_PER_SAMPLE_EMBEDDED=1" \
      ./RecoilJets_Condor_submit.sh isSimEmbedded condorDoAllSmoke groupSize 1 maxJobs 1 SAMPLE=run28_embeddedPhoton12
    ;;
  production-submit)
    preflight
    check_unused_path "$signal_bulk" "production bulk output"
    check_unused_path "$merge_base" "production merge output"
    exec > >(tee -a "$submit_log") 2>&1
    condor_q -nobatch "${USER:-patsfan753}" 2>/dev/null | tail -40 || true
    env "${common_env[@]}" \
      "RJ_SIMEMBED_DEST_BASE=${signal_bulk}" \
      "RJ_MERGE_OUT_BASE_OVERRIDE=${merge_base}" \
      ./RecoilJets_Condor_submit.sh isSimEmbedded condorDoAll groupSize "$sim_group"
    condor_q -nobatch "${USER:-patsfan753}" 2>/dev/null | tail -80 || true
    ;;
  status)
    printf 'QUEUE_ROWS\n'
    condor_q "${USER:-patsfan753}" -af ClusterId ProcId JobStatus Args 2>/dev/null | grep -F "$campaign_tag" || true
    printf 'OUTPUTS\n'
    find "$canary_merge" "$merge_base" -type f -name '*.root' -size +100k 2>/dev/null | sort || true
    ;;
  heartbeat-prompt)
    cat <<EOF
Read-only heartbeat for THE-98 AuAu EMCal PHOSUB isolation campaign.

Campaign tag: ${campaign_tag}
Config: ${yaml}
Manifest: ${manifest}
Canary merge root: ${canary_merge}
Production merge root: ${merge_base}

Objective: finish the variantB canary gate, then monitor the explicitly submitted
Photon12+20 production through merge, pull, ROOT QA, and the two local R=0.3/R=0.4 slides.

Allowed: read-only condor_q, condor_history, log, find, and ROOT inspection.
Forbidden without fresh approval: resubmission, deletion, Condor job control,
data/inclusive-jet submissions, or Google Slides mutation.
EOF
    ;;
  help|-h|--help)
    usage
    ;;
  *)
    usage >&2
    die "unknown mode: $mode"
    ;;
esac
