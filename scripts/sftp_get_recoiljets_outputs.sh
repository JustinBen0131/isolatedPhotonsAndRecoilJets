#!/usr/bin/env bash
set -euo pipefail

LOCAL_BASE="/Users/patsfan753/Desktop/ThesisAnalysis"
REMOTE_BASE="/sphenix/u/patsfan753/scratch/thesisAnalysis"
REMOTE_HOST="patsfan753@sftp.sdcc.bnl.gov"

usage() {
  cat <<'EOF'
Usage:
  ./scripts/sftp_get_recoiljets_outputs.sh <dataset>
  ./scripts/sftp_get_recoiljets_outputs.sh trainingLatest <tight|npb|jetResidual>
  ./scripts/sftp_get_recoiljets_outputs.sh tightBDTSmokeLatest
  ./scripts/sftp_get_recoiljets_outputs.sh tightBDTSmoke <remote-path-or-dir-name>
  ./scripts/sftp_get_recoiljets_outputs.sh mlIntegrationLatest
  ./scripts/sftp_get_recoiljets_outputs.sh mlIntegration <remote-path-or-dir-name>
  ./scripts/sftp_get_recoiljets_outputs.sh smokeTestLatest <dataset> [--roots]
  ./scripts/sftp_get_recoiljets_outputs.sh smokeTest <dataset> <remote-path-or-dir-name> [--roots]
  ./scripts/sftp_get_recoiljets_outputs.sh scaledTriggerStudy

Datasets:
  isAuAu                    -> InputFiles/auau25
  isPP                      -> InputFiles/pp24
  isPPrun25                 -> InputFiles/pp25
  isSim                     -> InputFiles/simPhotonJet
  isSimEmbedded             -> InputFiles/simEmbedded
  isSimEmbeddedInclusive    -> InputFiles/InclusiveJetSIM_EMBEDDED
  isSimInclusive            -> InputFiles/InclusiveJetSIM
  isSimJet5                 -> InputFiles/InclusiveJetSIM
  mergeLocalSim [dataset|all] -> build canonical local merged SIM products
                                from already-present InputFiles without sftp

The script builds the expected final ROOT filenames from the current
macros/analysis_config.yaml matrix, prints the local overwrite preview, then
opens interactive sftp once to fetch those files. No password is stored.
For SIM datasets it then materializes canonical merged outputs under
dataOutput/combinedSimOnly or dataOutput/combinedSimOnlyEMBEDDED.

trainingLatest pulls the newest SDCC local smoke-test training output directory
from:
  /sphenix/u/patsfan753/scratch/thesisAnalysis/local_bdt_training_outputs
and the matching model directory from:
  /sphenix/u/patsfan753/scratch/thesisAnalysis/bdt_models

tightBDTSmokeLatest/tightBDTSmoke are focused aliases for the tight-BDT-only
local/smoke workflow. They download into:
  InputFiles/trainingSmoke/tightBDT

mlIntegrationLatest pulls the newest full local ML integration test directory
from:
  /sphenix/u/patsfan753/scratch/thesisAnalysis/local_ml_pipeline_tests

mlIntegration with an explicit path or directory name bypasses latest-directory
discovery and downloads that exact run.

smokeTestLatest/smokeTest do the same for the overnight per-dataset smokeTest
workflow under thesisAnaSmoke. Pass --roots only when you intentionally want the
disposable smoke ROOT files too.

scaledTriggerStudy pulls the one-off AuAu scaled-trigger final ROOT file into:
  InputFiles/auau25
EOF
}

trim_ws() {
  local s="$1"
  s="${s#"${s%%[![:space:]]*}"}"
  s="${s%"${s##*[![:space:]]}"}"
  printf "%s\n" "$s"
}

yaml_path() {
  printf "%s\n" "${LOCAL_BASE}/macros/analysis_config.yaml"
}

make_tmp_file() {
  local prefix="$1"
  local tmp_dir="${LOCAL_BASE}/.tmp"
  mkdir -p "$tmp_dir"
  mktemp "${tmp_dir}/${prefix}.XXXXXX"
}

yaml_get_scalar_bool() {
  local yaml="$1" key="$2" default="$3"
  local line val
  line="$(grep -E "^[[:space:]]*${key}:" "$yaml" | head -n 1 || true)"
  if [[ -z "$line" ]]; then
    printf "%s\n" "$default"
    return 0
  fi
  val="${line#*:}"
  val="${val%%#*}"
  val="$(trim_ws "$val")"
  [[ -n "$val" ]] || val="$default"
  printf "%s\n" "$val"
}

yaml_get_inline_list() {
  local yaml="$1" key="$2"
  local line inner
  line="$(grep -E "^[[:space:]]*${key}:" "$yaml" | head -n 1 || true)"
  [[ -n "$line" ]] || return 0
  inner="${line#*:}"
  inner="${inner%%#*}"
  inner="$(trim_ws "$inner")"
  inner="${inner#[}"
  inner="${inner%]}"
  awk -v s="$inner" '
    BEGIN{
      n=split(s,a,",");
      for(i=1;i<=n;i++){
        gsub(/^[[:space:]]+|[[:space:]]+$/,"",a[i]);
        gsub(/^["'\'']|["'\'']$/,"",a[i]);
        if(a[i]!="") print a[i];
      }
    }'
}

sim_is_close() {
  awk -v a="$1" -v b="$2" 'BEGIN{d=a-b; if (d<0) d=-d; exit(d<1e-9 ? 0 : 1)}'
}

sim_pt_tag() {
  local pt="$1"
  if [[ "$pt" =~ ^([0-9]+)\.0+$ ]]; then
    echo "${BASH_REMATCH[1]}"
  else
    local s="$pt"
    s="${s//./p}"
    s="${s//-/m}"
    echo "$s"
  fi
}

sim_b2b_dir_tag() {
  local frac="$1"
  if sim_is_close "$frac" "0.5"; then
    echo "pi_2"
  elif sim_is_close "$frac" "0.875"; then
    echo "7pi_8"
  else
    local s="$frac"
    s="${s//./p}"
    s="${s//-/m}"
    echo "piFrac${s}"
  fi
}

dphi_internal_enabled() {
  case "${RJ_DISABLE_DPHI_INTERNALIZATION:-0}" in
    1|true|TRUE|yes|YES|on|ON) return 1 ;;
  esac
  case "${RJ_INTERNALIZE_DPHI:-1}" in
    0|false|FALSE|no|NO|off|OFF) return 1 ;;
  esac
  return 0
}

dphi_submit_values() {
  if dphi_internal_enabled && (( "$#" > 0 )); then
    local v
    for v in "$@"; do
      if sim_is_close "$v" "0.875"; then
        printf '%s\n' "$v"
        return 0
      fi
    done
    printf '%s\n' "$1"
  else
    printf '%s\n' "$@"
  fi
}

dphi_dir_tag_component() {
  local frac="$1"
  if dphi_internal_enabled; then
    echo "dphiScan"
  else
    sim_b2b_dir_tag "$frac"
  fi
}

sim_vz_tag() {
  local vz="$1"
  if [[ "$vz" =~ ^([0-9]+)\.0+$ ]]; then
    echo "vz${BASH_REMATCH[1]}"
  else
    local s="$vz"
    s="${s//./p}"
    s="${s//-/m}"
    echo "vz${s}"
  fi
}

sim_cone_tag() {
  local r="$1"
  local r100
  r100="$(awk -v r="$r" 'BEGIN{v=int(r*100+0.5); printf "%d", v}')"
  echo "isoR${r100}"
}

sim_iso_tag() {
  local sliding="$1" fixed="$2"
  if [[ "$sliding" == "true" ]]; then
    echo "isSliding"
  elif [[ "$fixed" =~ ^([0-9]+)\.0+$ ]]; then
    echo "fixedIso${BASH_REMATCH[1]}GeV"
  else
    local s="$fixed"
    s="${s//./p}"
    echo "fixedIso${s}GeV"
  fi
}

is_merge_dataset() {
  case "$1" in
    isSim|isSimEmbedded|isSimEmbeddedInclusive|isSimInclusive)
      return 0
      ;;
    *)
      return 1
      ;;
  esac
}

local_input_dir_for_merge_dataset() {
  case "$1" in
    isSim) echo "${LOCAL_BASE}/InputFiles/simPhotonJet" ;;
    isSimEmbedded) echo "${LOCAL_BASE}/InputFiles/simEmbedded" ;;
    isSimEmbeddedInclusive) echo "${LOCAL_BASE}/InputFiles/InclusiveJetSIM_EMBEDDED" ;;
    isSimInclusive) echo "${LOCAL_BASE}/InputFiles/InclusiveJetSIM" ;;
    *) return 1 ;;
  esac
}

primary_sample_for_merge_dataset() {
  case "$1" in
    isSim) echo "photonjet5" ;;
    isSimEmbedded) echo "embeddedPhoton12" ;;
    isSimEmbeddedInclusive) echo "embeddedJet12" ;;
    isSimInclusive) echo "jet5" ;;
    *) return 1 ;;
  esac
}

required_samples_for_merge_dataset() {
  case "$1" in
    isSim) printf "%s\n" photonjet5 photonjet10 photonjet20 ;;
    isSimEmbedded) printf "%s\n" embeddedPhoton12 embeddedPhoton20 ;;
    isSimEmbeddedInclusive) printf "%s\n" embeddedJet12 embeddedJet20 ;;
    isSimInclusive) printf "%s\n" jet5 ;;
    *) return 1 ;;
  esac
}

discover_complete_local_sim_cfg_tags() {
  local label="$1" dir primary f base cfg sample ok
  is_merge_dataset "$label" || return 0
  dir="$(local_input_dir_for_merge_dataset "$label")"
  primary="$(primary_sample_for_merge_dataset "$label")"
  [[ -d "$dir" ]] || return 0

  shopt -s nullglob
  for f in "${dir}/RecoilJets_${primary}_ALL_"*.root; do
    base="${f##*/}"
    cfg="${base#RecoilJets_${primary}_ALL_}"
    cfg="${cfg%.root}"
    ok=1
    while IFS= read -r sample; do
      [[ -n "$sample" ]] || continue
      if [[ ! -f "${dir}/RecoilJets_${sample}_ALL_${cfg}.root" ]]; then
        ok=0
        break
      fi
    done < <(required_samples_for_merge_dataset "$label")
    if (( ok )); then
      printf "%s\n" "$cfg"
    fi
  done
  shopt -u nullglob
}

merge_recoiljets_sim_outputs() {
  local label="$1"
  shift || true
  local cfgs=( "$@" )
  local cfg_list root_cmd

  is_merge_dataset "$label" || return 0
  if (( ${#cfgs[@]} == 0 )); then
    echo "[MERGE] No complete local SIM cfg tags found for ${label}; skipping."
    return 0
  fi

  cfg_list="$(make_tmp_file "recoiljets_merge_cfgs")"
  printf "%s\n" "${cfgs[@]}" | sort -u > "$cfg_list"
  root_cmd="scripts/MergeDownloadedRecoilJetsSim.C(\"${label}\",\"${cfg_list}\")"

  echo
  echo "[MERGE] Building canonical local merged SIM output for ${label}"
  echo "[MERGE] cfg count: $(wc -l < "$cfg_list" | tr -d ' ')"
  echo "[MERGE] ROOT helper: ${root_cmd}"

  if ( cd "$LOCAL_BASE" && TMPDIR="${LOCAL_BASE}/.tmp" ./scripts/root_in_analysis_env.sh /Users/patsfan753/Desktop/analysis/env/bin/root -l -b -q "$root_cmd" ); then
    rm -f "$cfg_list"
    echo "[MERGE][OK] Canonical merged SIM outputs are current for ${label}."
  else
    local status=$?
    rm -f "$cfg_list"
    echo "[MERGE][ERROR] Canonical SIM merge failed for ${label}." >&2
    exit "$status"
  fi
}

merge_local_existing_sim_outputs() {
  local requested="${1:-all}"
  local labels=()
  local label cfgs=()

  case "$requested" in
    all|"")
      labels=( isSim isSimEmbedded isSimEmbeddedInclusive isSimInclusive )
      ;;
    isSim|isSimEmbedded|isSimEmbeddedInclusive|isSimInclusive)
      labels=( "$requested" )
      ;;
    *)
      echo "[ERROR] mergeLocalSim target must be one of: all, isSim, isSimEmbedded, isSimEmbeddedInclusive, isSimInclusive" >&2
      exit 2
      ;;
  esac

  for label in "${labels[@]}"; do
    cfgs=()
    while IFS= read -r cfg; do
      [[ -n "$cfg" ]] && cfgs+=( "$cfg" )
    done < <(discover_complete_local_sim_cfg_tags "$label" | sort -u)
    merge_recoiljets_sim_outputs "$label" "${cfgs[@]}"
  done
}

selection_mode_normalize() {
  selection_mode_normalize_for_key "" "$1"
}

selection_mode_normalize_for_key() {
  local key="$1"
  local mode
  mode="$(trim_ws "$2")"
  case "$key:$mode" in
    preselection:variantA|preselection:VariantA|preselection:varianta|preselection:newPPG12|preselection:NewPPG12|preselection:newppg12) echo "newPPG12"; return 0 ;;
    preselection:variantB|preselection:VariantB|preselection:variantb|preselection:noPreCriteria|preselection:NoPreCriteria|preselection:noprecriteria) echo "noPreCriteria"; return 0 ;;
    preselection:variantC|preselection:VariantC|preselection:variantc|preselection:onlyNPB|preselection:OnlyNPB|preselection:onlynpb) echo "onlyNPB"; return 0 ;;
    preselection:variantD|preselection:VariantD|preselection:variantd|preselection:refPlusNPB|preselection:RefPlusNPB|preselection:refplusnpb) echo "refPlusNPB"; return 0 ;;
    preselection:variantE|preselection:VariantE|preselection:variante|preselection:auauOnlyNPB|preselection:AuAuOnlyNPB|preselection:auauonlynpb) echo "auauOnlyNPB"; return 0 ;;
    tight:variantA|tight:VariantA|tight:varianta|tight:newPPG12|tight:NewPPG12|tight:newppg12) echo "newPPG12"; return 0 ;;
    tight:variantB|tight:VariantB|tight:variantb|tight:auauEmbeddedBDT|tight:AuAuEmbeddedBDT|tight:auauembeddedbdt) echo "auauEmbeddedBDT"; return 0 ;;
    nonTight:variantA|nonTight:VariantA|nonTight:varianta|nonTight:bdtSideband|nonTight:BDTSideband|nonTight:bdtsideband|nonTight:newPPG12|nonTight:NewPPG12|nonTight:newppg12) echo "newPPG12"; return 0 ;;
    nonTight:variantB|nonTight:VariantB|nonTight:variantb|nonTight:auauBDTSideband|nonTight:AuAuBDTSideband|nonTight:auaubdtsideband) echo "auauBDTSideband"; return 0 ;;
    nonTight:variantC|nonTight:VariantC|nonTight:variantc|nonTight:auauBDTComplement|nonTight:AuAuBDTComplement|nonTight:auaubdtcomplement) echo "auauBDTComplement"; return 0 ;;
  esac
  case "$mode" in
    ""|reference|Reference) echo "reference" ;;
    variantA|VariantA|varianta|newPPG12|NewPPG12|newppg12) echo "newPPG12" ;;
    auauEmbeddedBDT|AuAuEmbeddedBDT|auauembeddedbdt) echo "auauEmbeddedBDT" ;;
    auauBDTSideband|AuAuBDTSideband|auaubdtsideband) echo "auauBDTSideband" ;;
    auauBDTComplement|AuAuBDTComplement|auaubdtcomplement) echo "auauBDTComplement" ;;
    variantB|VariantB|variantb) echo "variantB" ;;
    *) echo "$mode" ;;
  esac
}

selection_mode_tag() {
  local key="$1"
  local mode
  mode="$(selection_mode_normalize_for_key "$key" "$2")"
  case "$mode" in
    reference) echo "${key}Reference" ;;
    newPPG12) echo "${key}NewPPG12" ;;
    noPreCriteria) echo "${key}NoPreCriteria" ;;
    onlyNPB) echo "${key}OnlyNPB" ;;
    refPlusNPB) echo "${key}RefPlusNPB" ;;
    auauOnlyNPB) echo "${key}AuAuOnlyNPB" ;;
    auauEmbeddedBDT) echo "${key}AuAuEmbeddedBDT" ;;
    auauBDTSideband) echo "${key}AuAuBDTSideband" ;;
    auauBDTComplement) echo "${key}AuAuBDTComplement" ;;
    variantB) echo "${key}VariantB" ;;
    *)
      awk -v key="$key" -v mode="$mode" '
        BEGIN{
          first=toupper(substr(mode,1,1));
          rest=substr(mode,2);
          print key first rest;
        }'
      ;;
  esac
}

yaml_get_photon_id_sets() {
  local yaml="$1"
  awk '
    function trim(s){ gsub(/^[[:space:]]+|[[:space:]]+$/, "", s); return s }
    BEGIN{inset=0}
    {
      line=$0
      sub(/#.*/, "", line)
      if (line ~ /^[[:space:]]*photon_id_sets[[:space:]]*:/) { inset=1; next }
      if (inset && line ~ /^[[:alnum:]_][[:alnum:]_[:space:]-]*:/) { exit }
      if (!inset || line !~ /\[/) next
      sub(/^.*\[/, "", line)
      sub(/\].*$/, "", line)
      n=split(line, a, ",")
      if (n >= 3) print trim(a[1]) "|" trim(a[2]) "|" trim(a[3])
    }' "$yaml"
}

build_iso_mode_tags_from_yaml() {
  local yaml="$1"
  local isSlidingIso isSlidingAndFixed
  isSlidingIso="$(yaml_get_scalar_bool "$yaml" "isSlidingIso" "false")"
  isSlidingAndFixed="$(yaml_get_scalar_bool "$yaml" "isSlidingAndFixed" "false")"

  local _fixeds=()
  while IFS= read -r line; do
    _fixeds+=( "$line" )
  done < <(yaml_get_inline_list "$yaml" "fixedGeV")
  if (( ${#_fixeds[@]} == 0 )); then
    _fixeds=( "2.0" )
  fi

  if [[ "$isSlidingAndFixed" == "true" ]]; then
    echo "isSliding"
    local _f
    for _f in "${_fixeds[@]}"; do
      sim_iso_tag "false" "$_f"
    done
    return 0
  fi

  if [[ "$isSlidingIso" == "true" ]]; then
    echo "isSliding"
    return 0
  fi

  local _f
  for _f in "${_fixeds[@]}"; do
    sim_iso_tag "false" "$_f"
  done
}

dataset_includes_uepipe_in_tag() {
  case "$1" in
    isAuAu|auau|AuAu|AUAU|isSimEmbedded|simembedded|SIMEMBEDDED|isSimEmbeddedInclusive|simembeddedinclusive|SIMEMBEDDEDINCLUSIVE)
      return 0
      ;;
    *)
      return 1
      ;;
  esac
}

build_cfg_tags_from_yaml() {
  local dataset_token="$1"
  local yaml
  yaml="$(yaml_path)"
  [[ -f "$yaml" ]] || { echo "[ERROR] YAML not found: $yaml" >&2; exit 40; }

  local -a jet_pts b2bs b2bs_submit vzs cones iso_base_tags uepipes
  jet_pts=()
  b2bs=()
  vzs=()
  cones=()
  iso_base_tags=()
  uepipes=()

  local line
  while IFS= read -r line; do jet_pts+=( "$line" ); done < <(yaml_get_inline_list "$yaml" "jet_pt_min")
  while IFS= read -r line; do b2bs+=( "$line" ); done < <(yaml_get_inline_list "$yaml" "back_to_back_dphi_min_pi_fraction")
  while IFS= read -r line; do vzs+=( "$line" ); done < <(yaml_get_inline_list "$yaml" "vz_cut_cm")
  while IFS= read -r line; do cones+=( "$line" ); done < <(yaml_get_inline_list "$yaml" "coneR")
  while IFS= read -r line; do iso_base_tags+=( "$line" ); done < <(build_iso_mode_tags_from_yaml "$yaml")
  local -a photon_id_rows=()
  while IFS= read -r line; do photon_id_rows+=( "$line" ); done < <(yaml_get_photon_id_sets "$yaml")
  (( ${#photon_id_rows[@]} > 0 )) || { echo "[ERROR] YAML must define photon_id_sets for cfg-tag generation: $yaml" >&2; exit 41; }

  (( ${#jet_pts[@]} )) || jet_pts=( "5.0" )
  (( ${#b2bs[@]} )) || b2bs=( "0.875" )
  while IFS= read -r line; do b2bs_submit+=( "$line" ); done < <(dphi_submit_values "${b2bs[@]}")
  (( ${#vzs[@]} )) || vzs=( "30.0" )
  (( ${#cones[@]} )) || cones=( "0.30" )
  (( ${#iso_base_tags[@]} )) || iso_base_tags=( "fixedIso2GeV" )

  if dataset_includes_uepipe_in_tag "$dataset_token"; then
    while IFS= read -r line; do uepipes+=( "$line" ); done < <(yaml_get_inline_list "$yaml" "clusterUEpipeline")
    (( ${#uepipes[@]} )) || uepipes=( "noSub" )
  else
    uepipes=( "noSub" )
  fi

  local pt frac vz cone iso uep pre tight nonTight tag selection_tag full_tag
  local cfg_suffix="${SFTP_GET_CFG_SUFFIX:-}"
  local tight_norm nonTight_norm pre_norm
  for pt in "${jet_pts[@]}"; do
    for frac in "${b2bs_submit[@]}"; do
      for vz in "${vzs[@]}"; do
        for cone in "${cones[@]}"; do
          for iso in "${iso_base_tags[@]}"; do
            local row
            for row in "${photon_id_rows[@]}"; do
              IFS='|' read -r pre tight nonTight <<< "$row"
              pre_norm="$(selection_mode_normalize_for_key "preselection" "$pre")"
              tight_norm="$(selection_mode_normalize_for_key "tight" "$tight")"
              nonTight_norm="$(selection_mode_normalize_for_key "nonTight" "$nonTight")"
              selection_tag="$(selection_mode_tag "preselection" "$pre_norm")_$(selection_mode_tag "tight" "$tight_norm")_$(selection_mode_tag "nonTight" "$nonTight_norm")"
              tag="jetMinPt$(sim_pt_tag "$pt")_$(dphi_dir_tag_component "$frac")_$(sim_vz_tag "$vz")_$(sim_cone_tag "$cone")_${iso}"
              for uep in "${uepipes[@]}"; do
                if dataset_includes_uepipe_in_tag "$dataset_token"; then
                  full_tag="${tag}_${uep}_${selection_tag}"
                else
                  full_tag="${tag}_${selection_tag}"
                fi
                echo "${full_tag}${cfg_suffix}"
              done
            done
          done
        done
      done
    done
  done | sort -u
}

training_prefix_for_task() {
  case "$1" in
    tight|trainTightBDT) echo "tight_" ;;
    npb|trainNPB) echo "npb_" ;;
    jetResidual|jetML|trainJetMLResidual) echo "jetResidual_" ;;
    *)
      echo "[ERROR] Unknown training task: $1" >&2
      echo "[ERROR] Use one of: tight, npb, jetResidual" >&2
      exit 2
      ;;
  esac
}

training_latest_remote_dir() {
  local remote_parent="$1"
  local prefix="$2"
  local optional="${3:-false}"
  local ls_batch ls_out latest
  ls_batch="$(make_tmp_file "sftp_get_recoiljets_ls")"
  {
    printf 'ls %s/%s*\n' "$remote_parent" "$prefix"
  } > "$ls_batch"

  set +e
  ls_out="$(sftp \
      -oBatchMode=no \
      -oPreferredAuthentications=publickey,password,keyboard-interactive \
      -b "$ls_batch" \
      "$REMOTE_HOST" 2>&1)"
  local status=$?
  set -e
  rm -f "$ls_batch"

  if (( status != 0 )); then
    if [[ "$optional" == "true" ]]; then
      return 1
    fi
    echo "[ERROR] Could not list ${remote_parent}/${prefix}* on SDCC." >&2
    echo "$ls_out" >&2
    exit "$status"
  fi

  latest="$(printf "%s\n" "$ls_out" \
    | awk -v p="$prefix" '
        {
          for (i=1; i<=NF; i++) {
            v=$i;
            gsub(/:$/, "", v);
            gsub(/\/$/, "", v);
            n=split(v,a,"/");
            v=a[n];
            gsub(/^[[:space:]]+|[[:space:]]+$/,"",v);
            if (v ~ "^" p) print v;
          }
        }' \
    | sort \
    | tail -n 1)"

  if [[ -z "$latest" ]]; then
    if [[ "$optional" == "true" ]]; then
      return 1
    fi
    echo "[ERROR] No remote training directories found under ${remote_parent}/${prefix}*" >&2
    echo "[ERROR] Raw sftp listing output:" >&2
    printf "%s\n" "$ls_out" >&2
    exit 4
  fi
  printf "%s\n" "$latest"
}

download_training_latest() {
  local task="$1"
  local prefix train_parent model_parent latest_train latest_model local_dir batch
  prefix="$(training_prefix_for_task "$task")"
  train_parent="${REMOTE_BASE}/local_bdt_training_outputs"
  model_parent="${REMOTE_BASE}/bdt_models"
  latest_train="$(training_latest_remote_dir "$train_parent" "$prefix")"

  latest_model=""
  if latest_model="$(training_latest_remote_dir "$model_parent" "$prefix" true 2>/dev/null)"; then
    :
  else
    latest_model=""
  fi

  local_dir="${LOCAL_BASE}/InputFiles/trainingSmoke/${prefix%_}"
  mkdir -p "$local_dir"
  batch="$(make_tmp_file "sftp_get_recoiljets_training")"
  cleanup_training() { rm -f "$batch"; }
  trap cleanup_training EXIT

  {
    printf 'lcd %s\n' "$local_dir"
    printf 'get -r %s/%s %s\n' "$train_parent" "$latest_train" "$latest_train"
    if [[ -n "$latest_model" ]]; then
      printf 'get -r %s/%s %s\n' "$model_parent" "$latest_model" "$latest_model"
    fi
  } > "$batch"

  echo
  echo "Remote host       : ${REMOTE_HOST}"
  echo "Training remote   : ${train_parent}/${latest_train}"
  if [[ -n "$latest_model" ]]; then
    echo "Model remote      : ${model_parent}/${latest_model}"
  else
    echo "Model remote      : not found yet"
  fi
  echo "Local dir         : ${local_dir}"
  echo
  echo "This will overwrite matching local files/directories."
  read -r -p "Continue? [y/N]: " confirm
  case "$confirm" in
    y|Y|yes|YES|Yes) ;;
    *) echo "Aborted."; exit 0 ;;
  esac

  echo
  echo "sftp batch commands:"
  sed 's/^/  /' "$batch"
  echo
  echo "Opening interactive sftp to download selected training outputs."
  if sftp \
      -oBatchMode=no \
      -oPreferredAuthentications=publickey,password,keyboard-interactive \
      -b "$batch" \
      "$REMOTE_HOST"; then
    echo
    echo "[OK] Training output download complete."
    echo "Downloaded into: ${local_dir}"
    trap - EXIT
    rm -f "$batch"
  else
    status=$?
    echo
    echo "[ERROR] sftp download failed with exit code ${status}." >&2
    exit "$status"
  fi
}

download_tight_bdt_smoke() {
  local requested="${1:-}"
  local parent prefix latest remote_dir local_dir batch
  parent="${REMOTE_BASE}/local_bdt_training_outputs"
  prefix="tight"
  if [[ -n "$requested" ]]; then
    if [[ "$requested" == /* ]]; then
      remote_dir="${requested%/}"
      latest="${remote_dir##*/}"
    else
      latest="${requested%/}"
      remote_dir="${parent}/${latest}"
    fi
  else
    latest="$(training_latest_remote_dir "$parent" "$prefix")"
    remote_dir="${parent}/${latest}"
  fi

  local_dir="${LOCAL_BASE}/InputFiles/trainingSmoke/tightBDT"
  mkdir -p "$local_dir"
  batch="$(make_tmp_file "sftp_get_recoiljets_tightbdt")"
  cleanup_tightbdt() { rm -f "$batch"; }
  trap cleanup_tightbdt EXIT

  {
    printf 'lcd %s\n' "$local_dir"
    printf 'get -r %s %s\n' "$remote_dir" "$latest"
    printf 'get -r %s/bdt_models/%s %s_models\n' "$REMOTE_BASE" "$latest" "$latest"
  } > "$batch"

  echo
  echo "Remote host     : ${REMOTE_HOST}"
  echo "Tight-BDT remote: ${remote_dir}"
  echo "Model remote    : ${REMOTE_BASE}/bdt_models/${latest}"
  echo "Local dir       : ${local_dir}"
  echo
  echo "This will overwrite matching local files/directories."
  read -r -p "Continue? [y/N]: " confirm
  case "$confirm" in
    y|Y|yes|YES|Yes) ;;
    *) echo "Aborted."; exit 0 ;;
  esac

  echo
  echo "sftp batch commands:"
  sed 's/^/  /' "$batch"
  echo
  echo "Opening interactive sftp to download selected tight-BDT outputs."
  if sftp \
      -oBatchMode=no \
      -oPreferredAuthentications=publickey,password,keyboard-interactive \
      -b "$batch" \
      "$REMOTE_HOST"; then
    echo
    echo "[OK] Tight-BDT smoke download complete."
    echo "Downloaded into: ${local_dir}"
    trap - EXIT
    rm -f "$batch"
  else
    status=$?
    echo
    echo "[ERROR] sftp download failed with exit code ${status}." >&2
    exit "$status"
  fi
}

download_ml_integration_latest() {
  local requested="${1:-}"
  local parent prefix latest remote_dir local_dir batch
  parent="${REMOTE_BASE}/local_ml_pipeline_tests"
  prefix="mlIntegration_"
  if [[ -n "$requested" ]]; then
    if [[ "$requested" == /* ]]; then
      remote_dir="${requested%/}"
      latest="${remote_dir##*/}"
    else
      latest="${requested%/}"
      remote_dir="${parent}/${latest}"
    fi
    if [[ "$latest" != ${prefix}* ]]; then
      echo "[ERROR] ML integration directory must start with ${prefix}: ${latest}" >&2
      exit 2
    fi
  else
    latest="$(training_latest_remote_dir "$parent" "$prefix")"
    remote_dir="${parent}/${latest}"
  fi
  local_dir="${LOCAL_BASE}/InputFiles/trainingSmoke/mlIntegration"
  mkdir -p "$local_dir"
  batch="$(make_tmp_file "sftp_get_recoiljets_mlint")"
  cleanup_mlint() { rm -f "$batch"; }
  trap cleanup_mlint EXIT

  {
    printf 'lcd %s\n' "$local_dir"
    printf 'get -r %s %s\n' "$remote_dir" "$latest"
  } > "$batch"

  echo
  echo "Remote host       : ${REMOTE_HOST}"
  echo "Integration remote: ${remote_dir}"
  echo "Local dir         : ${local_dir}"
  echo
  echo "This will overwrite matching local files/directories."
  read -r -p "Continue? [y/N]: " confirm
  case "$confirm" in
    y|Y|yes|YES|Yes) ;;
    *) echo "Aborted."; exit 0 ;;
  esac

  echo
  echo "sftp batch commands:"
  sed 's/^/  /' "$batch"
  echo
  echo "Opening interactive sftp to download selected ML integration outputs."
  if sftp \
      -oBatchMode=no \
      -oPreferredAuthentications=publickey,password,keyboard-interactive \
      -b "$batch" \
      "$REMOTE_HOST"; then
    echo
    echo "[OK] ML integration download complete."
    echo "Downloaded into: ${local_dir}"
    trap - EXIT
    rm -f "$batch"
  else
    status=$?
    echo
    echo "[ERROR] sftp download failed with exit code ${status}." >&2
    exit "$status"
  fi
}

resolve_smoke_dataset() {
  case "$1" in
    isAuAu|auau|AuAu|AUAU)
      SMOKE_LABEL="isAuAu"; SMOKE_REMOTE_TAG="auau" ;;
    isPP|pp|PP)
      SMOKE_LABEL="isPP"; SMOKE_REMOTE_TAG="pp" ;;
    isPPrun25|pprun25|pp25|PP25)
      SMOKE_LABEL="isPPrun25"; SMOKE_REMOTE_TAG="pp25" ;;
    isSim|sim|SIM)
      SMOKE_LABEL="isSim"; SMOKE_REMOTE_TAG="sim" ;;
    isSimEmbedded|simembedded|SIMEMBEDDED)
      SMOKE_LABEL="isSimEmbedded"; SMOKE_REMOTE_TAG="simembedded" ;;
    isSimEmbeddedInclusive|simembeddedinclusive|SIMEMBEDDEDINCLUSIVE)
      SMOKE_LABEL="isSimEmbeddedInclusive"; SMOKE_REMOTE_TAG="simembeddedinclusive" ;;
    isSimInclusive|siminclusive|SIMINCLUSIVE|isSimJet5|simjet5|SIMJET5)
      SMOKE_LABEL="isSimInclusive"; SMOKE_REMOTE_TAG="simjet5" ;;
    *)
      echo "[ERROR] Unknown smokeTest dataset: $1" >&2
      echo "[ERROR] Use one of: isPP, isAuAu, isSim, isSimEmbedded, isSimEmbeddedInclusive, isSimInclusive" >&2
      exit 2
      ;;
  esac
}

download_smoke_test_report() {
  local requested_dataset="${1:-}"
  local requested="${2:-}"
  local smoke_kind="${3:-smokeTest}"
  local include_roots="${4:-0}"
  local parent prefix latest remote_dir local_dir batch
  [[ -n "$requested_dataset" ]] || { echo "[ERROR] ${smoke_kind}Latest requires a dataset, e.g. isPP" >&2; exit 2; }
  resolve_smoke_dataset "$requested_dataset"
  parent="/sphenix/tg/tg01/bulk/jbennett/thesisAnaSmoke"
  prefix="${SMOKE_REMOTE_TAG}_smokeTest_"

  if [[ -n "$requested" ]]; then
    if [[ "$requested" == /* ]]; then
      remote_dir="${requested%/}"
      latest="${remote_dir##*/}"
    else
      latest="${requested%/}"
      remote_dir="${parent}/${latest}"
    fi
  else
    latest="$(training_latest_remote_dir "$parent" "$prefix")"
    remote_dir="${parent}/${latest}"
  fi

  if [[ "$latest" != ${prefix}* ]]; then
    echo "[ERROR] ${smoke_kind} directory must start with ${prefix}: ${latest}" >&2
    exit 2
  fi

  local_dir="${LOCAL_BASE}/InputFiles/pipelineSmoke/${SMOKE_LABEL}"
  mkdir -p "$local_dir"
  batch="$(make_tmp_file "sftp_get_recoiljets_${smoke_kind}")"
  cleanup_smoketest() { rm -f "$batch"; }
  trap cleanup_smoketest EXIT

  {
    printf 'lcd %s\n' "$local_dir"
    printf 'get -r %s/_pipeline_reports %s_reports\n' "$remote_dir" "$latest"
    if [[ "$include_roots" == "1" ]]; then
      printf 'get -r %s %s_roots\n' "$remote_dir" "$latest"
    fi
  } > "$batch"

  echo
  echo "Remote host       : ${REMOTE_HOST}"
  echo "Smoke remote      : ${remote_dir}"
  echo "Report remote     : ${remote_dir}/_pipeline_reports"
  echo "Local dir         : ${local_dir}"
  echo
  if [[ "$include_roots" == "1" ]]; then
    echo "This downloads reports plus disposable smoke ROOT outputs."
  else
    echo "This downloads only small smoke tuning/report files, not ROOT outputs."
  fi
  read -r -p "Continue? [y/N]: " confirm
  case "$confirm" in
    y|Y|yes|YES|Yes) ;;
    *) echo "Aborted."; exit 0 ;;
  esac

  echo
  echo "sftp batch commands:"
  sed 's/^/  /' "$batch"
  echo
  echo "Opening interactive sftp to download selected ${smoke_kind} reports."
  if sftp \
      -oBatchMode=no \
      -oPreferredAuthentications=publickey,password,keyboard-interactive \
      -b "$batch" \
      "$REMOTE_HOST"; then
    echo
    echo "[OK] ${smoke_kind} report download complete."
    echo "Downloaded into: ${local_dir}/${latest}_reports"
    [[ "$include_roots" == "1" ]] && echo "Smoke ROOT copy: ${local_dir}/${latest}_roots"
    trap - EXIT
    rm -f "$batch"
  else
    status=$?
    echo
    echo "[ERROR] sftp download failed with exit code ${status}." >&2
    exit "$status"
  fi
}

download_scaled_trigger_study() {
  local cfg file remote_dir local_dir batch
  cfg="jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionReference_tightReference_nonTightReference_scaledTriggerStudy"
  file="RecoilJets_auau_ALL_${cfg}.root"
  remote_dir="${REMOTE_BASE}/output/auau"
  local_dir="${LOCAL_BASE}/InputFiles/auau25"

  mkdir -p "$local_dir"
  batch="$(make_tmp_file "sftp_get_scaled_trigger_study")"
  cleanup_scaled_trigger() { rm -f "$batch"; }
  trap cleanup_scaled_trigger EXIT

  echo
  echo "Remote host : ${REMOTE_HOST}"
  echo "Remote file : ${remote_dir}/${file}"
  echo "Local file  : ${local_dir}/${file}"
  echo "Study       : scaledTriggerStudy"
  echo

  if [[ -e "${local_dir}/${file}" ]]; then
    echo "Existing local file found:"
    echo "  ${local_dir}/${file}"
    echo
    echo "Choose how to handle it before download:"
    echo "  o = overwrite in place"
    echo "  p = move existing file to ${local_dir}/previous/<timestamp>/ first"
    echo "  a = abort"
    read -r -p "Action? [o/p/a]: " conflict_action
    case "$conflict_action" in
      o|O|overwrite|OVERWRITE)
        echo "Proceeding with overwrite in place."
        ;;
      p|P|previous|PREVIOUS)
        previous_dir="${local_dir}/previous/$(date +%Y%m%d_%H%M%S)_scaledTriggerStudy"
        mkdir -p "$previous_dir"
        echo "Moving existing file to: ${previous_dir}"
        mv "${local_dir}/${file}" "${previous_dir}/"
        ;;
      a|A|abort|ABORT|"")
        echo "Aborted."
        exit 0
        ;;
      *)
        echo "[ERROR] Unknown action: ${conflict_action}" >&2
        exit 2
        ;;
    esac
  fi

  echo "Proceeding without an extra y/N prompt."

  {
    printf 'lcd %s\n' "$local_dir"
    printf 'cd %s\n' "$remote_dir"
    printf 'get %s %s\n' "$file" "$file"
  } > "$batch"

  echo
  echo "sftp batch commands:"
  sed 's/^/  /' "$batch"
  echo
  echo "Opening interactive sftp to download the scaled-trigger study ROOT file."
  if sftp \
      -oBatchMode=no \
      -oPreferredAuthentications=publickey,password,keyboard-interactive \
      -b "$batch" \
      "$REMOTE_HOST"; then
    echo
    echo "[OK] scaledTriggerStudy download complete."
    echo "Downloaded into: ${local_dir}/${file}"
    trap - EXIT
    rm -f "$batch"
  else
    status=$?
    echo
    echo "[ERROR] sftp download failed with exit code ${status}." >&2
    exit "$status"
  fi
}

dataset="${1:-}"
case "$dataset" in
  -h|--help|help|"")
    usage
    [[ -n "$dataset" ]] && exit 0
    exit 2
    ;;
esac

if [[ "$dataset" == "trainingLatest" || "$dataset" == "trainingSmoke" ]]; then
  download_training_latest "${2:-}"
  exit 0
fi

if [[ "$dataset" == "mergeLocalSim" ]]; then
  merge_local_existing_sim_outputs "${2:-all}"
  exit 0
fi

if [[ "$dataset" == "tightBDTSmokeLatest" ]]; then
  download_tight_bdt_smoke
  exit 0
fi

if [[ "$dataset" == "tightBDTSmoke" ]]; then
  download_tight_bdt_smoke "${2:-}"
  exit 0
fi

if [[ "$dataset" == "mlIntegrationLatest" ]]; then
  download_ml_integration_latest
  exit 0
fi

if [[ "$dataset" == "mlIntegration" ]]; then
  download_ml_integration_latest "${2:-}"
  exit 0
fi

if [[ "$dataset" == "smokeTestLatest" ]]; then
  include_roots=0
  [[ "${3:-}" == "--roots" ]] && include_roots=1
  download_smoke_test_report "${2:-}" "" "smokeTest" "$include_roots"
  exit 0
fi

if [[ "$dataset" == "smokeTest" ]]; then
  include_roots=0
  [[ "${4:-}" == "--roots" ]] && include_roots=1
  download_smoke_test_report "${2:-}" "${3:-}" "smokeTest" "$include_roots"
  exit 0
fi

if [[ "$dataset" == "scaledTriggerStudy" ]]; then
  download_scaled_trigger_study
  exit 0
fi

remote_tag=""
local_subdir=""
sample_tags=()
label=""

case "$dataset" in
  isAuAu|auau|AuAu|AUAU)
    label="isAuAu"
    remote_tag="auau"
    local_subdir="InputFiles/auau25"
    sample_tags=( "auau" )
    ;;
  isPP|pp|PP)
    label="isPP"
    remote_tag="pp"
    local_subdir="InputFiles/pp24"
    sample_tags=( "pp" )
    ;;
  isPPrun25|pprun25|pp25|PP25)
    label="isPPrun25"
    remote_tag="pp25"
    local_subdir="InputFiles/pp25"
    sample_tags=( "pp25" )
    ;;
  isSim|sim|SIM)
    label="isSim"
    remote_tag="sim"
    local_subdir="InputFiles/simPhotonJet"
    sample_tags=( "photonjet5" "photonjet10" "photonjet20" )
    ;;
  isSimEmbedded|simembedded|SIMEMBEDDED)
    label="isSimEmbedded"
    remote_tag="simembedded"
    local_subdir="InputFiles/simEmbedded"
    sample_tags=( "embeddedPhoton12" "embeddedPhoton20" )
    ;;
  isSimEmbeddedInclusive|simembeddedinclusive|SIMEMBEDDEDINCLUSIVE)
    label="isSimEmbeddedInclusive"
    remote_tag="simembeddedinclusive"
    local_subdir="InputFiles/InclusiveJetSIM_EMBEDDED"
    sample_tags=( "embeddedJet12" "embeddedJet20" )
    ;;
  isSimInclusive|siminclusive|SIMINCLUSIVE|isSimJet5|simjet5|SIMJET5)
    label="isSimInclusive"
    remote_tag="simjet5"
    local_subdir="InputFiles/InclusiveJetSIM"
    sample_tags=( "jet5" )
    ;;
  *)
    echo "[ERROR] Unknown dataset: $dataset" >&2
    echo >&2
    usage >&2
    exit 2
    ;;
esac

remote_dir="${REMOTE_BASE}/output/${remote_tag}"
local_dir="${LOCAL_BASE}/${local_subdir}"

if [[ ! -d "$LOCAL_BASE" ]]; then
  echo "[ERROR] Local base does not exist: $LOCAL_BASE" >&2
  exit 3
fi

mkdir -p "$local_dir"

get_batch="$(make_tmp_file "sftp_get_recoiljets_get")"
cleanup() {
  rm -f "$get_batch"
}
trap cleanup EXIT

cfg_tags=()
while IFS= read -r cfg; do
  cfg_tags+=( "$cfg" )
done < <(build_cfg_tags_from_yaml "$dataset")

if (( ${#cfg_tags[@]} == 0 )); then
  echo "[ERROR] No cfg tags could be built from $(yaml_path)" >&2
  exit 5
fi

remote_files=()
local_files=()
existing_files=()
for cfg in "${cfg_tags[@]}"; do
  for sample in "${sample_tags[@]}"; do
    file="RecoilJets_${sample}_ALL_${cfg}.root"
    remote_files+=( "$file" )
    local_files+=( "${local_dir}/${file}" )
    if [[ -e "${local_dir}/${file}" ]]; then
      existing_files+=( "${local_dir}/${file}" )
    fi
  done
done

echo
echo "Remote host : ${REMOTE_HOST}"
echo "Remote dir  : ${remote_dir}"
echo "Local dir   : ${local_dir}"
echo "Variant     : ${label}"
echo "Config YAML : $(yaml_path)"
echo
echo "Samples:"
for sample in "${sample_tags[@]}"; do
  echo "  ${sample}"
done
echo
echo "Files to download (${#remote_files[@]}):"
for f in "${remote_files[@]}"; do
  echo "  ${remote_dir}/${f}"
  echo "    -> ${local_dir}/${f}"
done
echo
if (( ${#existing_files[@]} > 0 )); then
  echo "Existing local files with the same requested tag were found (${#existing_files[@]}):"
  for f in "${existing_files[@]}"; do
    echo "  ${f}"
  done
  echo
  echo "Choose how to handle existing files before download:"
  echo "  o = overwrite in place"
  echo "  p = move existing files to ${local_dir}/previous/<timestamp>/ first"
  echo "  a = abort"
  read -r -p "Action? [o/p/a]: " conflict_action
  case "$conflict_action" in
    o|O|overwrite|OVERWRITE)
      echo "Proceeding with overwrite in place."
      ;;
    p|P|previous|PREVIOUS)
      previous_dir="${local_dir}/previous/$(date +%Y%m%d_%H%M%S)_${label}"
      mkdir -p "$previous_dir"
      echo "Moving existing files to: ${previous_dir}"
      for f in "${existing_files[@]}"; do
        mv "$f" "${previous_dir}/"
      done
      ;;
    a|A|abort|ABORT|"")
      echo "Aborted."
      exit 0
      ;;
    *)
      echo "[ERROR] Unknown action: ${conflict_action}" >&2
      exit 2
      ;;
  esac
else
  echo "No matching local files exist for this requested tag set."
fi

echo
echo "Ready to download ${#remote_files[@]} file(s)."
read -r -p "Continue? [y/N]: " confirm
case "$confirm" in
  y|Y|yes|YES|Yes)
    ;;
  *)
    echo "Aborted."
    exit 0
    ;;
esac

{
  printf 'lcd %s\n' "$local_dir"
  printf 'cd %s\n' "$remote_dir"
  for f in "${remote_files[@]}"; do
    printf 'get %s %s\n' "$f" "$f"
  done
} > "$get_batch"

echo
echo "sftp batch commands:"
sed 's/^/  /' "$get_batch"
echo
echo "Opening interactive sftp to download selected files."
echo "Public-key auth is tried first; no password is needed if your SDCC key is installed."
if sftp \
    -oBatchMode=no \
    -oPreferredAuthentications=publickey,password,keyboard-interactive \
    -b "$get_batch" \
    "$REMOTE_HOST"; then
  echo
  echo "[OK] Download complete."
else
  status=$?
  echo
  echo "[ERROR] sftp download failed with exit code ${status}." >&2
  echo "[ERROR] No success confirmation was received from sftp." >&2
  exit "$status"
fi

echo
echo "Downloaded ${#remote_files[@]} file(s) into: ${local_dir}"

if is_merge_dataset "$label"; then
  merge_recoiljets_sim_outputs "$label" "${cfg_tags[@]}"
fi
