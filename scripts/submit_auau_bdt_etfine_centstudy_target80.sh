#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

campaign_tag="${RJ_ETFINE_CAMPAIGN_TAG:-etfine_centstudy_target80_$(date +%Y%m%d_%H%M%S)}"
target_eff="${TARGET_SIGNAL_EFF:-0.80}"
wp_pt_bins="${RJ_AUAU_BDT_WP_PT_BINS:-15,17,19,21,23,25,27,30,35}"
wp_cent_bins="${RJ_AUAU_BDT_WP_CENT_BINS:-0,20,50,80}"
target_label="target$(python3 - <<PY
print(f"{int(round(100*float('${target_eff}'))):02d}")
PY
)"
model_dir="${RJ_ETFINE_MODEL_DIR:-/gpfs/mnt/gpfs02/sphenix/user/${USER:-patsfan753}/thesisAnalysis/bdt_models/tight_etfine_centstudy_current}"
validation_report="${RJ_ETFINE_VALIDATION_REPORT:-}"
template_yaml="${RJ_ETFINE_TEMPLATE_YAML:-macros/analysis_config_auau_bdt_etfine_centstudy_wp050.yaml}"
config_dir="${RJ_ETFINE_CONFIG_DIR:-${repo_root}/condor_generated_configs/${campaign_tag}}"
group_size="${RJ_ETFINE_GROUP_SIZE:-7}"
notify="${RJ_NOTIFY_EMAILS:-just0131@gmail.com}"
request_memory="${RJ_REQUEST_MEMORY:-12000MB}"
retry_cap="${RJ_AUTO_MEMORY_RETRY_CAP_MB:-16000}"
merge_memory="${RJ_SIM_FIRSTROUND_REQUEST_MEMORY:-8000MB}"
merge_group="${RJ_SIM_MERGE_GROUP_SIZE:-75}"
auto_merge="${RJ_ETFINE_AUTO_MERGE:-0}"

mkdir -p "$config_dir"
submit_log="${config_dir}/submit_${campaign_tag}.log"
exec > >(tee -a "$submit_log") 2>&1

echo "RECOILJETS_AUAU_BDT_ETFINE_CENTSTUDY_TARGET80_SUBMIT_V1"
echo "submit_host=$(hostname -f 2>/dev/null || hostname)"
echo "repo_root=${repo_root}"
echo "campaign_tag=${campaign_tag}"
echo "model_dir=${model_dir}"
echo "validation_report=${validation_report}"
echo "template_yaml=${template_yaml}"
echo "config_dir=${config_dir}"
echo "target_signal_efficiency=${target_eff}"
echo "working_point_mode=centpt"
echo "working_point_pt_bins=${wp_pt_bins}"
echo "working_point_cent_bins=${wp_cent_bins}"
echo "valid_pt_window=15:35"
echo "group_size=${group_size}"
echo "request_memory=${request_memory}"
echo "retry_cap=${retry_cap}"
echo "merge_memory=${merge_memory}"
echo "merge_group_size=${merge_group}"
echo "auto_merge=${auto_merge}"
echo "notify=${notify}"
if [[ "$auto_merge" != "1" ]]; then
  echo "analysis_only=1"
  echo "note=RJ_ETFINE_AUTO_MERGE defaults to 0 so failed workers stay inspectable; run printed merge commands after worker clusters are clean."
fi
echo

[[ -f "$template_yaml" ]] || { echo "[ERROR] Template YAML not found: $template_yaml" >&2; exit 2; }
[[ -n "$validation_report" ]] || { echo "[ERROR] Set RJ_ETFINE_VALIDATION_REPORT=/path/to/model_validation_report" >&2; exit 2; }
[[ -d "$validation_report" ]] || { echo "[ERROR] Validation report directory not found: $validation_report" >&2; exit 2; }

echo "== deriving working points =="
RJ_ML_PYTHON="${RJ_ML_PYTHON:-/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python}" \
./scripts/auau_tight_bdt_pipeline.sh deriveWorkingPointsFromValidation \
  VALIDATION="$validation_report" \
  TARGET="$target_eff" \
  MODE=centpt \
  PT_BINS="$wp_pt_bins" \
  CENT_BINS="$wp_cent_bins"

wp_json="${validation_report}/bdt_working_points_${target_label}.json"
[[ -f "$wp_json" ]] || { echo "[ERROR] Expected working-point JSON not produced: $wp_json" >&2; exit 3; }

product_map="auauCentInputBase3x3BDT=centInput_pt1535,auauEtFineCentInputBDT=ptFine_centInput,auauEtFineCent3BDT=ptFine_cent3,auauEtFineCent7BDT=ptFine_cent7"
yaml="${config_dir}/analysis_config_auau_bdt_etfine_centstudy_${target_label}.yaml"
./scripts/auau_tight_bdt_pipeline.sh generateWorkingPointConfig \
  TEMPLATE="$template_yaml" \
  WORKING_POINTS="$wp_json" \
  PRODUCT_MAP="$product_map" \
  OUT="$yaml" \
  MODEL_DIR="$model_dir"
echo "[OK] wrote frozen target-WP YAML: ${yaml}"

required_models=(
  "${model_dir}/auau_tight_bdt_centInput_pt1535_tmva.root"
)
pt_edges=(15 17 19 21 23 25 27 30 35)
cent3_edges=(0 20 50 80)
cent7_edges=(0 10 20 30 40 50 60 80)
fmt_pt() { printf 'pt_%03d_%03d' "$1" "$2"; }
fmt_cent() { printf 'cent_%03d_%03d' "$1" "$2"; }
for ((ip=0; ip+1<${#pt_edges[@]}; ip++)); do
  pt_tag="$(fmt_pt "${pt_edges[$ip]}" "${pt_edges[$((ip+1))]}")"
  required_models+=( "${model_dir}/auau_tight_bdt_ptFine_centInput_${pt_tag}_tmva.root" )
  for ((ic=0; ic+1<${#cent3_edges[@]}; ic++)); do
    cent_tag="$(fmt_cent "${cent3_edges[$ic]}" "${cent3_edges[$((ic+1))]}")"
    required_models+=( "${model_dir}/auau_tight_bdt_ptFine_cent3_${pt_tag}_${cent_tag}_tmva.root" )
  done
  for ((ic=0; ic+1<${#cent7_edges[@]}; ic++)); do
    cent_tag="$(fmt_cent "${cent7_edges[$ic]}" "${cent7_edges[$((ic+1))]}")"
    required_models+=( "${model_dir}/auau_tight_bdt_ptFine_cent7_${pt_tag}_${cent_tag}_tmva.root" )
  done
done

missing=0
echo "== model preflight =="
for model in "${required_models[@]}"; do
  if [[ -f "$model" ]]; then
    echo "  [OK] $model"
  else
    echo "  [MISSING] $model"
    missing=1
  fi
done
if (( missing )); then
  echo "[ERROR] Missing one or more fine-E_T BDT model ROOT files; train/validate before submitting MC." >&2
  exit 3
fi

submit_one() {
  local dataset="$1"
  local bulk_var="$2"
  local bulk_base="$3"
  local merge_base="$4"
  echo
  echo "====================================================================="
  echo "Submitting ${dataset}"
  echo "  campaign   : ${campaign_tag}"
  echo "  yaml       : ${yaml}"
  echo "  bulk base  : ${bulk_base}"
  echo "  merge base : ${merge_base}"
  echo "====================================================================="
  env \
    RJ_NOTIFY_EMAILS="${notify}" \
    RJ_PROFILE_JOB=1 \
    RJ_CONFIG_YAML="${yaml}" \
    RJ_ID_FANOUT_MAX_ROWS=1 \
    RJ_REQUEST_MEMORY="${request_memory}" \
    RJ_AUTO_MEMORY_RETRY_CAP_MB="${retry_cap}" \
    RJ_HOLD_FAILED_WORKERS=1 \
    RJ_AUTO_MERGE="${auto_merge}" \
    RJ_SIM_FIRSTROUND_REQUEST_MEMORY="${merge_memory}" \
    RJ_SIM_MERGE_GROUP_SIZE="${merge_group}" \
    RJ_MERGE_OUT_BASE_OVERRIDE="${merge_base}" \
    "${bulk_var}=${bulk_base}" \
    ./RecoilJets_Condor_submit.sh "${dataset}" condorDoAll groupSize "${group_size}"
}

echo
echo "== preflight queue snapshot =="
condor_q -nobatch "${USER:-patsfan753}" 2>/dev/null | tail -30 || true

merge_base="/sphenix/u/patsfan753/scratch/thesisAnalysis/output_${campaign_tag}"
signal_bulk="/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded_${campaign_tag}"
background_bulk="/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembeddedinclusive_${campaign_tag}"

submit_one isSimEmbedded RJ_SIMEMBED_DEST_BASE "$signal_bulk" "$merge_base"
submit_one isSimEmbeddedInclusive RJ_SIMEMBEDINCLUSIVE_DEST_BASE "$background_bulk" "$merge_base"

if [[ "$auto_merge" != "1" ]]; then
  cat <<EOF

== merge commands after all matching worker jobs are clean ==
env MERGE_CONFIG_YAML=${yaml} MERGE_RUN_BASE_OVERRIDE=${signal_bulk} MERGE_OUT_BASE_OVERRIDE=${merge_base} RJ_SIM_FIRSTROUND_REQUEST_MEMORY=${merge_memory} ./scripts/mergeRecoilJets.sh isSimEmbedded firstRound groupSize ${merge_group}
env MERGE_CONFIG_YAML=${yaml} MERGE_RUN_BASE_OVERRIDE=${signal_bulk} MERGE_OUT_BASE_OVERRIDE=${merge_base} RJ_SIM_FIRSTROUND_REQUEST_MEMORY=${merge_memory} ./scripts/mergeRecoilJets.sh isSimEmbedded secondRound groupSize ${merge_group}
env MERGE_CONFIG_YAML=${yaml} MERGE_RUN_BASE_OVERRIDE=${signal_bulk} MERGE_OUT_BASE_OVERRIDE=${merge_base} RJ_SIM_FIRSTROUND_REQUEST_MEMORY=${merge_memory} ./scripts/mergeRecoilJets.sh isSimEmbedded finalStitch groupSize ${merge_group}
env MERGE_CONFIG_YAML=${yaml} MERGE_RUN_BASE_OVERRIDE=${background_bulk} MERGE_OUT_BASE_OVERRIDE=${merge_base} RJ_SIM_FIRSTROUND_REQUEST_MEMORY=${merge_memory} ./scripts/mergeRecoilJets.sh isSimEmbeddedInclusive firstRound groupSize ${merge_group}
env MERGE_CONFIG_YAML=${yaml} MERGE_RUN_BASE_OVERRIDE=${background_bulk} MERGE_OUT_BASE_OVERRIDE=${merge_base} RJ_SIM_FIRSTROUND_REQUEST_MEMORY=${merge_memory} ./scripts/mergeRecoilJets.sh isSimEmbeddedInclusive secondRound groupSize ${merge_group}
env MERGE_CONFIG_YAML=${yaml} MERGE_RUN_BASE_OVERRIDE=${background_bulk} MERGE_OUT_BASE_OVERRIDE=${merge_base} RJ_SIM_FIRSTROUND_REQUEST_MEMORY=${merge_memory} ./scripts/mergeRecoilJets.sh isSimEmbeddedInclusive finalStitch groupSize ${merge_group}
EOF
fi

echo
echo "== post-submit queue snapshot =="
condor_q -nobatch "${USER:-patsfan753}" 2>/dev/null | tail -60 || true
echo
echo "TRACK_THIS_CAMPAIGN=${campaign_tag}"
echo "TARGET_WP_JSON=${wp_json}"
echo "GENERATED_YAML=${yaml}"
echo "SUBMIT_LOG=${submit_log}"
