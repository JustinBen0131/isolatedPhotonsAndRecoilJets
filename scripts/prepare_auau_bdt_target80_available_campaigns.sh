#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

tag="${RJ_TARGET80_PREP_TAG:-target80_available_$(date +%Y%m%d_%H%M%S)}"
out_dir="${RJ_TARGET80_CONFIG_DIR:-${repo_root}/condor_generated_configs/${tag}}"
target="${TARGET_SIGNAL_EFF:-0.80}"
wp_mode="${RJ_AUAU_BDT_WP_MODE:-centpt}"
default_wp_pt_bins="${RJ_AUAU_BDT_WP_PT_BINS:-}"
default_wp_cent_bins="${RJ_AUAU_BDT_WP_CENT_BINS:-0,20,50,80}"
do_submit="${RJ_TARGET80_DO_SUBMIT:-0}"
mkdir -p "$out_dir"

echo "RECOILJETS_AUAU_BDT_TARGET80_AVAILABLE_PREP_V1"
echo "host=$(hostname -f 2>/dev/null || hostname)"
echo "repo_root=${repo_root}"
echo "tag=${tag}"
echo "out_dir=${out_dir}"
echo "target_signal_efficiency=${target}"
echo "working_point_mode=${wp_mode}"
echo "default_working_point_pt_bins=${default_wp_pt_bins:-validator default}"
echo "default_working_point_cent_bins=${default_wp_cent_bins}"
echo "do_submit=${do_submit}"
echo

make_one() {
  local name="$1"
  local validation="$2"
  local template="$3"
  local model_dir="$4"
  local product_map="$5"
  local pt_bins="${6:-$default_wp_pt_bins}"
  local cent_bins="${7:-$default_wp_cent_bins}"

  if [[ -z "$validation" ]]; then
    echo "[SKIP] ${name}: validation report variable not set"
    return 0
  fi
  if [[ ! -d "$validation" ]]; then
    echo "[SKIP] ${name}: validation report directory missing: $validation"
    return 0
  fi
  if [[ ! -f "$template" ]]; then
    echo "[SKIP] ${name}: template missing: $template"
    return 0
  fi

  echo
  echo "====================================================================="
  echo "Preparing ${name}"
  echo "  validation : ${validation}"
  echo "  template   : ${template}"
  echo "  model dir  : ${model_dir}"
  echo "  wp pt bins : ${pt_bins:-validator default}"
  echo "  wp cent bins: ${cent_bins}"
  echo "====================================================================="

  local target_label
  target_label="$(python3 - <<PY
print(f"target{int(round(100*float('${target}'))):02d}")
PY
)"
  local wp_json="${validation}/bdt_working_points_${target_label}.json"
  if [[ -f "$wp_json" && "${RJ_TARGET80_FORCE_DERIVE:-0}" != "1" ]]; then
    echo "[REUSE] ${name}: existing working-point JSON: ${wp_json}"
    echo "[REUSE] set RJ_TARGET80_FORCE_DERIVE=1 to regenerate it"
  else
    RJ_ML_PYTHON="${RJ_ML_PYTHON:-/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python}" \
    ./scripts/auau_tight_bdt_pipeline.sh deriveWorkingPointsFromValidation \
      VALIDATION="$validation" \
      TARGET="$target" \
      MODE="$wp_mode" \
      PT_BINS="$pt_bins" \
      CENT_BINS="$cent_bins"
  fi
  local yaml="${out_dir}/analysis_config_${name}_${target_label}.yaml"
  [[ -f "$wp_json" ]] || { echo "[ERROR] missing expected working-point JSON: $wp_json" >&2; exit 3; }

  local artifact_dir="${out_dir}/working_point_artifacts/${name}"
  mkdir -p "$artifact_dir"
  for artifact in \
    "${validation}/bdt_working_points_${target_label}.json" \
    "${validation}/bdt_working_points_${target_label}.yaml" \
    "${validation}/bdt_working_points_${target_label}.csv" \
    "${validation}/bdt_working_points_${target_label}_runtime_fragment.yaml" \
    "${validation}/bdt_working_point_${target_label}_diagnostics.png"
  do
    if [[ -f "$artifact" ]]; then
      cp -f "$artifact" "${artifact_dir}/"
    fi
  done

  local filtered_product_map
  filtered_product_map="$(python3 - "$wp_json" "$product_map" <<'PY'
import json
import sys

wp_json, raw_map = sys.argv[1], sys.argv[2]
with open(wp_json, "r", encoding="utf-8") as f:
    payload = json.load(f)

products = set()
if isinstance(payload, dict):
    for key in ("products", "working_points", "entries"):
        block = payload.get(key)
        if isinstance(block, dict):
            products.update(str(k) for k in block.keys())
        elif isinstance(block, list):
            for item in block:
                if isinstance(item, dict):
                    for pkey in ("product", "product_id", "id", "name"):
                        if item.get(pkey):
                            products.add(str(item[pkey]))

kept = []
missing = []
for part in raw_map.split(","):
    part = part.strip()
    if not part:
        continue
    if "=" not in part:
        missing.append(f"{part} (malformed)")
        continue
    variant, product = part.split("=", 1)
    variant = variant.strip()
    product = product.strip()
    if product in products:
        kept.append(f"{variant}={product}")
    else:
        missing.append(f"{variant}={product}")

for item in missing:
    print(f"[WARN] dropping unavailable product map entry for this validation: {item}", file=sys.stderr)
print(",".join(kept))
PY
)"
  if [[ -z "$filtered_product_map" ]]; then
    echo "[SKIP] ${name}: no requested product-map entries are present in ${wp_json}"
    return 0
  fi

  ./scripts/auau_tight_bdt_pipeline.sh generateWorkingPointConfig \
    TEMPLATE="$template" \
    WORKING_POINTS="$wp_json" \
    PRODUCT_MAP="$filtered_product_map" \
    OUT="$yaml" \
    MODEL_DIR="$model_dir"

  echo "[READY] ${name}_yaml=${yaml}"
  echo "[READY] ${name}_wp_json=${wp_json}"
  echo "[READY] ${name}_wp_artifacts=${artifact_dir}"
  echo "[NEXT] inspect ${validation}/bdt_working_point_${target_label}_diagnostics.png"
  echo "[NEXT] local test:"
  echo "RJ_CONFIG_YAML=${yaml} ./scripts/RecoilJets_Condor_submit.sh isSimEmbedded local 1000 NFILES=1 SAMPLE=run28_embeddedPhoton20 VERBOSE=1"
  echo "[NEXT] paired MC submit after inspection:"
  echo "RJ_NOTIFY_EMAILS=just0131@gmail.com RJ_TARGETWP_CONFIG_YAML=${yaml} ./scripts/submit_auau_bdt_targetwp_pair.sh"

  if [[ "$do_submit" == "1" ]]; then
    RJ_NOTIFY_EMAILS="${RJ_NOTIFY_EMAILS:-just0131@gmail.com}" \
    RJ_TARGETWP_CONFIG_YAML="$yaml" \
    RJ_TARGETWP_CAMPAIGN_TAG="${name}_${target_label}_$(date +%Y%m%d_%H%M%S)" \
    ./scripts/submit_auau_bdt_targetwp_pair.sh
  fi
}

expanded_map="auauNoCentBDT=centINDcontrol_pt5to40,auauCentInputBDT=centAsFeat_pt5to40,auauCentInput3x3BDT=centAsFeat3x3_pt5to40,auauCentInputMinOptBDT=centAsFeatMinOpt_pt5to40,auauCent3BDT=centDepBDTs_pt5to40,auauCent7BDT=centDepFineBDTs_pt5to40,auauPtBinCentInputBDT=ptBinCentAsFeat,auauPtCent3BDT=ptCentDep3,auauPtCent7BDT=ptCentDepFine"
width1530_map="auauCentInputBDT=centAsFeat_pt15to30,auauCentInput3x3BDT=centAsFeat3x3_pt15to30,auauCentInputBase3x3BDT=centAsFeatBase3x3_pt15to30"
etfine_map="auauCentInputBase3x3BDT=centInput_pt1535,auauEtFineCentInputBDT=ptFine_centInput,auauEtFineCent3BDT=ptFine_cent3,auauEtFineCent7BDT=ptFine_cent7"
width_windows_csv="${RJ_TARGET80_WIDTHWINDOWS_WINDOWS:-5:35,10:35,15:35}"
width_windows_model_dir="${RJ_TARGET80_WIDTHWINDOWS_MODEL_DIR:-/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_centinput_widthstudy_windows_current}"
width_windows_validation="${RJ_TARGET80_WIDTHWINDOWS_VALIDATION:-}"
width_windows_template_dir="${out_dir}/widthstudy_window_templates"

pt_tag() {
  local lo="$1"
  local hi="$2"
  printf 'pt%sto%s' "${lo%.*}" "${hi%.*}"
}

bins_for_window_csv() {
  local tag="$1"
  case "$tag" in
    pt5to35) printf '5,8,10,12,14,16,18,20,22,24,26,30,35' ;;
    pt10to35) printf '10,12,14,16,18,20,22,24,26,30,35' ;;
    pt15to35) printf '15,16,18,20,22,24,26,30,35' ;;
    *) echo "No target-WP pT binning preset for ${tag}; add it before preparing." >&2; return 2 ;;
  esac
}

make_one \
  "expanded_5to40" \
  "${RJ_TARGET80_EXPANDED_VALIDATION:-}" \
  "macros/analysis_config_auau_bdt_validation.yaml" \
  "${RJ_TARGET80_EXPANDED_MODEL_DIR:-/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_expanded_20260509_152604}" \
  "$expanded_map" \
  "${RJ_TARGET80_EXPANDED_PT_BINS:-5,8,10,12,14,16,18,20,22,24,26,35,40}" \
  "${RJ_TARGET80_EXPANDED_CENT_BINS:-$default_wp_cent_bins}"

make_one \
  "widthstudy_pt1530" \
  "${RJ_TARGET80_WIDTH1530_VALIDATION:-}" \
  "macros/analysis_config_auau_bdt_widthstudy_pt1530_wp050.yaml" \
  "${RJ_TARGET80_WIDTH1530_MODEL_DIR:-/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_centinput_widthstudy_pt1530_current}" \
  "$width1530_map" \
  "${RJ_TARGET80_WIDTH1530_PT_BINS:-15,17,19,21,23,25,27,30}" \
  "${RJ_TARGET80_WIDTH1530_CENT_BINS:-$default_wp_cent_bins}"

make_one \
  "etfine_15to35" \
  "${RJ_TARGET80_ETFINE_VALIDATION:-}" \
  "macros/analysis_config_auau_bdt_etfine_centstudy_wp050.yaml" \
  "${RJ_TARGET80_ETFINE_MODEL_DIR:-/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/tight_etfine_centstudy_current}" \
  "$etfine_map" \
  "${RJ_TARGET80_ETFINE_PT_BINS:-15,17,19,21,23,25,27,30,35}" \
  "${RJ_TARGET80_ETFINE_CENT_BINS:-$default_wp_cent_bins}"

if [[ -n "$width_windows_validation" ]]; then
  echo
  echo "====================================================================="
  echo "Preparing width-window templates"
  echo "  validation : ${width_windows_validation}"
  echo "  model dir  : ${width_windows_model_dir}"
  echo "  windows    : ${width_windows_csv}"
  echo "  template dir: ${width_windows_template_dir}"
  echo "====================================================================="
  RJ_WIDTHSTUDY_PREP_ONLY=1 \
  RJ_WIDTHSTUDY_MODEL_DIR="$width_windows_model_dir" \
  RJ_WIDTHSTUDY_CONFIG_DIR="$width_windows_template_dir" \
  RJ_WIDTHSTUDY_WINDOWS="$width_windows_csv" \
  RJ_WIDTHSTUDY_CAMPAIGN_TAG="${tag}_widthstudy_templates" \
  ./scripts/submit_auau_bdt_widthstudy_windows_wp050.sh

  IFS=',' read -r -a window_items <<< "$width_windows_csv"
  for raw in "${window_items[@]}"; do
    item="${raw//[[:space:]]/}"
    [[ -n "$item" ]] || continue
    lo="${item%%:*}"
    hi="${item##*:}"
    wtag="$(pt_tag "$lo" "$hi")"
    template="${width_windows_template_dir}/analysis_config_auau_bdt_widthstudy_${wtag}_wp050.yaml"
    product_map="auauCentInputBDT=centAsFeat_${wtag},auauCentInput3x3BDT=centAsFeat3x3_${wtag},auauCentInputBase3x3BDT=centAsFeatBase3x3_${wtag}"
    make_one \
      "widthstudy_${wtag}" \
      "$width_windows_validation" \
      "$template" \
      "$width_windows_model_dir" \
      "$product_map" \
      "$(bins_for_window_csv "$wtag")" \
      "${RJ_TARGET80_WIDTHWINDOWS_CENT_BINS:-$default_wp_cent_bins}"
  done
else
  echo "[SKIP] widthstudy_windows: RJ_TARGET80_WIDTHWINDOWS_VALIDATION not set"
fi

cat <<EOF

Summary:
  Prepared configs live in: ${out_dir}
  Default mode only stages configs and prints local-test/submit commands.
  To submit automatically after preparation, rerun with RJ_TARGET80_DO_SUBMIT=1.
EOF
