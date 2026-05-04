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

Datasets:
  isAuAu                    -> InputFiles/auau25
  isPP                      -> InputFiles/pp24
  isPPrun25                 -> InputFiles/pp25
  isSim                     -> InputFiles/simPhotonJet
  isSimEmbedded             -> InputFiles/simEmbedded
  isSimEmbeddedInclusive    -> InputFiles/InclusiveJetSIM_EMBEDDED
  isSimInclusive            -> InputFiles/InclusiveJetSIM
  isSimJet5                 -> InputFiles/InclusiveJetSIM

The script builds the expected final ROOT filenames from the current
macros/analysis_config.yaml matrix, prints the local overwrite preview, then
opens interactive sftp once to fetch those files. No password is stored.

trainingLatest pulls the newest SDCC local smoke-test training output directory
from:
  /sphenix/u/patsfan753/scratch/thesisAnalysis/local_bdt_training_outputs
and the matching model directory from:
  /sphenix/u/patsfan753/scratch/thesisAnalysis/bdt_models
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

selection_mode_normalize() {
  local mode
  mode="$(trim_ws "$1")"
  case "$mode" in
    ""|reference|Reference) echo "reference" ;;
    variantA|VariantA|varianta) echo "variantA" ;;
    variantB|VariantB|variantb) echo "variantB" ;;
    *) echo "$mode" ;;
  esac
}

selection_mode_tag() {
  local key="$1"
  local mode
  mode="$(selection_mode_normalize "$2")"
  case "$mode" in
    reference) echo "${key}Reference" ;;
    variantA) echo "${key}VariantA" ;;
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

  local -a jet_pts b2bs vzs cones iso_base_tags uepipes preselection_modes tight_modes nonTight_modes
  jet_pts=()
  b2bs=()
  vzs=()
  cones=()
  iso_base_tags=()
  uepipes=()
  preselection_modes=()
  tight_modes=()
  nonTight_modes=()

  local line
  while IFS= read -r line; do jet_pts+=( "$line" ); done < <(yaml_get_inline_list "$yaml" "jet_pt_min")
  while IFS= read -r line; do b2bs+=( "$line" ); done < <(yaml_get_inline_list "$yaml" "back_to_back_dphi_min_pi_fraction")
  while IFS= read -r line; do vzs+=( "$line" ); done < <(yaml_get_inline_list "$yaml" "vz_cut_cm")
  while IFS= read -r line; do cones+=( "$line" ); done < <(yaml_get_inline_list "$yaml" "coneR")
  while IFS= read -r line; do iso_base_tags+=( "$line" ); done < <(build_iso_mode_tags_from_yaml "$yaml")
  while IFS= read -r line; do preselection_modes+=( "$line" ); done < <(yaml_get_inline_list "$yaml" "preselection")
  while IFS= read -r line; do tight_modes+=( "$line" ); done < <(yaml_get_inline_list "$yaml" "tight")
  while IFS= read -r line; do nonTight_modes+=( "$line" ); done < <(yaml_get_inline_list "$yaml" "nonTight")

  (( ${#jet_pts[@]} )) || jet_pts=( "5.0" )
  (( ${#b2bs[@]} )) || b2bs=( "0.875" )
  (( ${#vzs[@]} )) || vzs=( "30.0" )
  (( ${#cones[@]} )) || cones=( "0.30" )
  (( ${#iso_base_tags[@]} )) || iso_base_tags=( "fixedIso2GeV" )
  (( ${#preselection_modes[@]} )) || preselection_modes=( "reference" )
  (( ${#tight_modes[@]} )) || tight_modes=( "reference" )
  (( ${#nonTight_modes[@]} )) || nonTight_modes=( "reference" )

  if dataset_includes_uepipe_in_tag "$dataset_token"; then
    while IFS= read -r line; do uepipes+=( "$line" ); done < <(yaml_get_inline_list "$yaml" "clusterUEpipeline")
    (( ${#uepipes[@]} )) || uepipes=( "noSub" )
  else
    uepipes=( "noSub" )
  fi

  local pt frac vz cone iso uep pre tight nonTight tag selection_tag
  for pt in "${jet_pts[@]}"; do
    for frac in "${b2bs[@]}"; do
      for vz in "${vzs[@]}"; do
        for cone in "${cones[@]}"; do
          for iso in "${iso_base_tags[@]}"; do
            for pre in "${preselection_modes[@]}"; do
              for tight in "${tight_modes[@]}"; do
                for nonTight in "${nonTight_modes[@]}"; do
                  selection_tag="$(selection_mode_tag "preselection" "$pre")_$(selection_mode_tag "tight" "$tight")_$(selection_mode_tag "nonTight" "$nonTight")"
                  tag="jetMinPt$(sim_pt_tag "$pt")_$(sim_b2b_dir_tag "$frac")_$(sim_vz_tag "$vz")_$(sim_cone_tag "$cone")_${iso}"
                  for uep in "${uepipes[@]}"; do
                    if dataset_includes_uepipe_in_tag "$dataset_token"; then
                      echo "${tag}_${uep}_${selection_tag}"
                    else
                      echo "${tag}_${selection_tag}"
                    fi
                  done
                done
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
  ls_batch="$(mktemp "${TMPDIR:-/tmp}/sftp_get_recoiljets_ls.XXXXXX")"
  {
    printf 'cd %s\n' "$remote_parent"
    printf 'ls -1 %s*\n' "$prefix"
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
        index($0,p)>0 {
          n=split($0,a,"/");
          v=a[n];
          gsub(/^[[:space:]]+|[[:space:]]+$/,"",v);
          if (v ~ "^" p) print v;
        }' \
    | sort \
    | tail -n 1)"

  if [[ -z "$latest" ]]; then
    if [[ "$optional" == "true" ]]; then
      return 1
    fi
    echo "[ERROR] No remote training directories found under ${remote_parent}/${prefix}*" >&2
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
  batch="$(mktemp "${TMPDIR:-/tmp}/sftp_get_recoiljets_training.XXXXXX")"
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
    sample_tags=( "embeddedJet10" "embeddedJet20" )
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

get_batch="$(mktemp "${TMPDIR:-/tmp}/sftp_get_recoiljets_get.XXXXXX")"
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
for cfg in "${cfg_tags[@]}"; do
  for sample in "${sample_tags[@]}"; do
    remote_files+=( "RecoilJets_${sample}_ALL_${cfg}.root" )
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
echo "This will overwrite matching local files."
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
