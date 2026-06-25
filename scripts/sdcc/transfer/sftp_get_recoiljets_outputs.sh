#!/usr/bin/env bash
set -euo pipefail

LOCAL_BASE="/Users/patsfan753/Desktop/ThesisAnalysis"
REMOTE_BASE="${SFTP_GET_REMOTE_BASE:-/sphenix/u/patsfan753/scratch/thesisAnalysis}"
REMOTE_HOST="patsfan753@sftp.sdcc.bnl.gov"

usage() {
  cat <<'EOF'
Usage:
  ./scripts/sftp_get_recoiljets_outputs.sh <dataset>
  ./scripts/sftp_get_recoiljets_outputs.sh trainingLatest <tight|npb|jetResidual>
  ./scripts/sftp_get_recoiljets_outputs.sh tightBDTSmokeLatest
  ./scripts/sftp_get_recoiljets_outputs.sh tightBDTSmoke <remote-path-or-dir-name>
  ./scripts/sftp_get_recoiljets_outputs.sh auauTightBDTValidation <remote-report-dir>
  ./scripts/sftp_get_recoiljets_outputs.sh auauTightBDTValidationScores <remote-report-dir> <local-dir>
  ./scripts/sftp_get_recoiljets_outputs.sh auauMLDiagnosticCompact <remote-dir> <local-dir> <file...>
  ./scripts/sftp_get_recoiljets_outputs.sh stitchDiagnosticsCompact <remote-dir> <local-dir> <file...>
  ./scripts/sftp_get_recoiljets_outputs.sh auauBDTMLPStackPromotion <remote-run-dir>
  ./scripts/sftp_get_recoiljets_outputs.sh freshOOFStackCompact <remote-run-dir> [local-dir]
  ./scripts/sftp_get_recoiljets_outputs.sh freshOOFStackDAGCompact <remote-run-dir> [local-dir]
  ./scripts/sftp_get_recoiljets_outputs.sh ppPhotonMLCompact <remote-dir> <local-dir> <file...>
  ./scripts/sftp_get_recoiljets_outputs.sh ppPhotonMLValidation <remote-validation-dir> [local-dir] [file ...]
  ./scripts/sftp_get_recoiljets_outputs.sh selectedRootFiles <remote-dir> <local-dir> <relative-root-file...>
  ./scripts/sftp_get_recoiljets_outputs.sh mlIntegrationLatest
  ./scripts/sftp_get_recoiljets_outputs.sh mlIntegration <remote-path-or-dir-name>
  ./scripts/sftp_get_recoiljets_outputs.sh smokeTestLatest <dataset> [--roots]
  ./scripts/sftp_get_recoiljets_outputs.sh smokeTest <dataset> <remote-path-or-dir-name> [--roots]
  ./scripts/sftp_get_recoiljets_outputs.sh smokeFinalLatest <dataset>
  ./scripts/sftp_get_recoiljets_outputs.sh smokeFinal <dataset> <remote-path-or-dir-name>
  ./scripts/sftp_get_recoiljets_outputs.sh scaledTriggerStudy
  ./scripts/sftp_get_recoiljets_outputs.sh scaledTriggerCentStudy
  ./scripts/sftp_get_recoiljets_outputs.sh scaledTriggerRunByRunQA <remote-dir> [local-dir]

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
Set SFTP_GET_CAMPAIGN_TAG to pull from the Stage-7 canonical
runs/recoiljets/current/<campaign>/<dataset> path. Set
SFTP_GET_REMOTE_DIR_OVERRIDE for an explicit legacy or non-default remote
output directory, and SFTP_GET_LOCAL_COMBINED_BASE to keep combined SIM pulls
in a campaign-specific local folder.

trainingLatest pulls the newest SDCC local smoke-test training output directory
from:
  /sphenix/u/patsfan753/scratch/thesisAnalysis/local_bdt_training_outputs
and the matching model directory from:
  /sphenix/u/patsfan753/scratch/thesisAnalysis/bdt_models

tightBDTSmokeLatest/tightBDTSmoke are focused aliases for the tight-BDT-only
local/smoke workflow. They download into:
  InputFiles/trainingSmoke/tightBDT

auauTightBDTValidation pulls a single finished simulation-validation report
directory into:
  dataOutput/auauTightBDTValidation

auauTightBDTValidationScores pulls only score_caches.list plus score_caches/
from an explicit AuAu tight-BDT validation report directory. It is intended for
regenerating working-point plots from event-level scores without pulling the
entire report.

auauMLDiagnosticCompact pulls selected compact PNG/JSON/CSV/TXT artifacts from
an explicit SDCC dataOutput/auauMLDiagnosticRuns directory.

stitchDiagnosticsCompact pulls selected compact PNG/JSON/CSV/TXT artifacts from
an explicit SDCC dataOutput/stitchDiagnostics directory.

auauBDTMLPStackPromotion pulls a single stack-promotion / WP-diagnostic run
directory into:
  dataOutput/auauBDTMLPStackPromotion

freshOOFStackCompact pulls only compact CSV/JSON artifacts from a fresh
pp/AuAu BDT+MLP OOF stack run under the SDCC mlp_models area into:
  dataOutput/fresh_pp_auau_oof_stack/<run-name>

freshOOFStackDAGCompact pulls the compact CSV/JSON plus locked-test NPZ
artifacts from a clean-slate staged DAG run under the SDCC mlp_models area into:
  dataOutput/fresh_pp_auau_oof_stack/<run-name>

ppPhotonMLValidation pulls selected compact validation artifacts from an
explicit ppPhotonMLPipeline validation directory. If file names are not passed,
it defaults to the corrected Fig. 19 NPB overlay PNG/JSON/LOG set.

selectedRootFiles pulls an explicit list of ROOT files from a known RecoilJets
output directory. It is intended for narrow diagnostics of partial or staged
outputs; it refuses wildcards, absolute file arguments, traversal, and
non-RecoilJets output roots.

ppPhotonMLCompact pulls selected compact PNG/JSON/CSV/TXT/LOG artifacts from
an explicit ppPhotonMLPipeline directory such as validation/insitu_stitching or
slide_assets.

mlIntegrationLatest pulls the newest full local ML integration test directory
from:
  /sphenix/u/patsfan753/scratch/thesisAnalysis/local_ml_pipeline_tests

mlIntegration with an explicit path or directory name bypasses latest-directory
discovery and downloads that exact run.

smokeTestLatest/smokeTest do the same for the overnight per-dataset smokeTest
workflow under thesisAnaSmoke. Pass --roots only when you intentionally want the
disposable smoke ROOT files too.

smokeFinalLatest/smokeFinal pull only the final merged smoke ROOT files from
outputSmoke/<dataset>_smokeTest_<timestamp>/<tag>, which is the right mode for
quick post-DAG output validation.

scaledTriggerStudy pulls the one-off AuAu scaled-trigger final ROOT file into:
  InputFiles/auau25

scaledTriggerCentStudy pulls the centrality-sliced AuAu scaled-trigger final
ROOT file into:
  InputFiles/auau25

scaledTriggerRunByRunQA pulls a curated scaled-trigger run-by-run QA artifact
directory produced on SDCC. It downloads summary CSV/Markdown files plus the
selected near-unity and problematic PNG subsets, not the full 620-image set.
EOF
}

trim_ws() {
  local s="$1"
  s="${s#"${s%%[![:space:]]*}"}"
  s="${s%"${s##*[![:space:]]}"}"
  printf "%s\n" "$s"
}

validate_remote_path() {
  local path="$1"
  local label="${2:-remote path}"
  case "$path" in
    ""|*[[:space:]]*|*agent_context*|*.codex*|*codex*|*THE-*)
      echo "[ERROR] Unsafe ${label}: ${path}" >&2
      echo "[ERROR] Refusing empty, whitespace, agent-context, codex, or THE-* SDCC paths." >&2
      exit 2
      ;;
  esac
  case "$path" in
    *"/../"*|*"//"*|*/..|../*)
      echo "[ERROR] Unsafe ${label}: ${path}" >&2
      echo "[ERROR] Refusing traversal-like SDCC paths." >&2
      exit 2
      ;;
  esac
}

validate_remote_relative_file() {
  local path="$1"
  local label="${2:-remote relative file}"
  case "$path" in
    ""|/*|*[[:space:]]*|*agent_context*|*.codex*|*codex*|*THE-*|*"*"*|*"?"*|*"["*|*"]"*)
      echo "[ERROR] Unsafe ${label}: ${path}" >&2
      echo "[ERROR] Refusing empty, absolute, wildcard, whitespace, agent-context, codex, or THE-* paths." >&2
      exit 2
      ;;
  esac
  case "$path" in
    *"/../"*|*"//"*|*/..|../*|..)
      echo "[ERROR] Unsafe ${label}: ${path}" >&2
      echo "[ERROR] Refusing traversal-like relative paths." >&2
      exit 2
      ;;
  esac
}

yaml_path() {
  local yaml="${SFTP_GET_CONFIG_YAML:-${RJ_CONFIG_YAML:-${LOCAL_BASE}/macros/analysis_config.yaml}}"
  case "$yaml" in
    /*) printf "%s\n" "$yaml" ;;
    *) printf "%s\n" "${LOCAL_BASE}/${yaml}" ;;
  esac
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
    isSimInclusive) printf "%s\n" jet5 jet8 jet12 jet20 jet30 jet40 ;;
    *) return 1 ;;
  esac
}

sim_combined_remote_file() {
  local label="$1" cfg="$2"
  local embedded_inclusive_dir="${SFTP_GET_EMBEDDED_INCLUSIVE_MERGED_DIR:-embeddedJet12and20merged_SIM}"
  local embedded_inclusive_file="${SFTP_GET_EMBEDDED_INCLUSIVE_MERGED_FILE:-RecoilJets_embeddedJet12plus20_MERGED.root}"
  case "$label" in
    isSim)
      echo "${cfg}/photonJet5and10and20merged_SIM/RecoilJets_photonjet5plus10plus20_MERGED.root" ;;
    isSimEmbedded)
      echo "${cfg}/photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root" ;;
    isSimEmbeddedInclusive)
      echo "${cfg}/${embedded_inclusive_dir}/${embedded_inclusive_file}" ;;
    isSimInclusive)
      echo "${cfg}/inclusiveJet5to40_SIM/RecoilJets_jet5plus8plus12plus20plus30plus40_MERGED.root" ;;
    *)
      return 1 ;;
  esac
}

sim_combined_local_file() {
  local label="$1" cfg="$2"
  local sim_base="${SFTP_GET_LOCAL_COMBINED_BASE:-${LOCAL_BASE}/dataOutput/combinedSimOnly}"
  local embed_base="${SFTP_GET_LOCAL_COMBINED_BASE:-${LOCAL_BASE}/dataOutput/combinedSimOnlyEMBEDDED}"
  local embedded_inclusive_dir="${SFTP_GET_EMBEDDED_INCLUSIVE_MERGED_DIR:-embeddedJet12and20merged_SIM}"
  local embedded_inclusive_file="${SFTP_GET_EMBEDDED_INCLUSIVE_MERGED_FILE:-RecoilJets_embeddedJet12plus20_MERGED.root}"
  case "$label" in
    isSim)
      echo "${sim_base}/${cfg}/photonJet5and10and20merged_SIM/RecoilJets_photonjet5plus10plus20_MERGED.root" ;;
    isSimEmbedded)
      echo "${embed_base}/${cfg}/photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root" ;;
    isSimEmbeddedInclusive)
      echo "${embed_base}/${cfg}/${embedded_inclusive_dir}/${embedded_inclusive_file}" ;;
    isSimInclusive)
      echo "${sim_base}/${cfg}/inclusiveJet5to40_SIM/RecoilJets_jet5plus8plus12plus20plus30plus40_MERGED.root" ;;
    *)
      return 1 ;;
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
    tight:auauNoCentBDT|tight:AuAuNoCentBDT|tight:auaunocentbdt) echo "auauNoCentBDT"; return 0 ;;
    tight:auauCentInputBDT|tight:AuAuCentInputBDT|tight:auaucentinputbdt) echo "auauCentInputBDT"; return 0 ;;
    tight:auauCentInput3x3BDT|tight:AuAuCentInput3x3BDT|tight:auaucentinput3x3bdt) echo "auauCentInput3x3BDT"; return 0 ;;
    tight:auauCentInputBase3x3BDT|tight:AuAuCentInputBase3x3BDT|tight:auaucentinputbase3x3bdt) echo "auauCentInputBase3x3BDT"; return 0 ;;
    tight:auauCentInputMinOptBDT|tight:AuAuCentInputMinOptBDT|tight:auaucentinputminoptbdt) echo "auauCentInputMinOptBDT"; return 0 ;;
    tight:auauCent3BDT|tight:AuAuCent3BDT|tight:auaucent3bdt) echo "auauCent3BDT"; return 0 ;;
    tight:auauCent7BDT|tight:AuAuCent7BDT|tight:auaucent7bdt) echo "auauCent7BDT"; return 0 ;;
    tight:auauPtBinCentInputBDT|tight:AuAuPtBinCentInputBDT|tight:auauptbincentinputbdt) echo "auauPtBinCentInputBDT"; return 0 ;;
    tight:auauPtCent3BDT|tight:AuAuPtCent3BDT|tight:auauptcent3bdt) echo "auauPtCent3BDT"; return 0 ;;
    tight:auauPtCent7BDT|tight:AuAuPtCent7BDT|tight:auauptcent7bdt) echo "auauPtCent7BDT"; return 0 ;;
    tight:auauEtFineCentInputBDT|tight:AuAuEtFineCentInputBDT|tight:auauetfinecentinputbdt) echo "auauEtFineCentInputBDT"; return 0 ;;
    tight:auauEtFineCent3BDT|tight:AuAuEtFineCent3BDT|tight:auauetfinecent3bdt) echo "auauEtFineCent3BDT"; return 0 ;;
    tight:auauEtFineCent7BDT|tight:AuAuEtFineCent7BDT|tight:auauetfinecent7bdt) echo "auauEtFineCent7BDT"; return 0 ;;
    tight:auauCentInputMLP|tight:AuAuCentInputMLP|tight:auaucentinputmlp) echo "auauCentInputMLP"; return 0 ;;
    tight:auauNoCentBase3x3MLP|tight:AuAuNoCentBase3x3MLP|tight:auaunocentbase3x3mlp) echo "auauNoCentBase3x3MLP"; return 0 ;;
    tight:auauCentInputBase3x3MLP|tight:AuAuCentInputBase3x3MLP|tight:auaucentinputbase3x3mlp) echo "auauCentInputBase3x3MLP"; return 0 ;;
    tight:auauBDTMLPStack|tight:AuAuBDTMLPStack|tight:auaubdtmlpstack|tight:bdtmlpstack) echo "auauBDTMLPStack"; return 0 ;;
    nonTight:variantA|nonTight:VariantA|nonTight:varianta|nonTight:bdtSideband|nonTight:BDTSideband|nonTight:bdtsideband|nonTight:newPPG12|nonTight:NewPPG12|nonTight:newppg12) echo "newPPG12"; return 0 ;;
    nonTight:variantB|nonTight:VariantB|nonTight:variantb|nonTight:auauBDTSideband|nonTight:AuAuBDTSideband|nonTight:auaubdtsideband) echo "auauBDTSideband"; return 0 ;;
    nonTight:variantC|nonTight:VariantC|nonTight:variantc|nonTight:auauBDTComplement|nonTight:AuAuBDTComplement|nonTight:auaubdtcomplement) echo "auauBDTComplement"; return 0 ;;
    nonTight:auauMLPSideband|nonTight:AuAuMLPSideband|nonTight:auaumlpsideband) echo "auauMLPSideband"; return 0 ;;
    nonTight:auauMLPComplement|nonTight:AuAuMLPComplement|nonTight:auaumlpcomplement) echo "auauMLPComplement"; return 0 ;;
    nonTight:auauBDTMLPStackSideband|nonTight:AuAuBDTMLPStackSideband|nonTight:auaubdtmlpstacksideband) echo "auauBDTMLPStackSideband"; return 0 ;;
    nonTight:auauBDTMLPStackComplement|nonTight:AuAuBDTMLPStackComplement|nonTight:auaubdtmlpstackcomplement) echo "auauBDTMLPStackComplement"; return 0 ;;
  esac
  case "$mode" in
    ""|reference|Reference) echo "reference" ;;
    variantA|VariantA|varianta|newPPG12|NewPPG12|newppg12) echo "newPPG12" ;;
    auauEmbeddedBDT|AuAuEmbeddedBDT|auauembeddedbdt) echo "auauEmbeddedBDT" ;;
    auauNoCentBDT|AuAuNoCentBDT|auaunocentbdt) echo "auauNoCentBDT" ;;
    auauCentInputBDT|AuAuCentInputBDT|auaucentinputbdt) echo "auauCentInputBDT" ;;
    auauCentInput3x3BDT|AuAuCentInput3x3BDT|auaucentinput3x3bdt) echo "auauCentInput3x3BDT" ;;
    auauCentInputBase3x3BDT|AuAuCentInputBase3x3BDT|auaucentinputbase3x3bdt) echo "auauCentInputBase3x3BDT" ;;
    auauCentInputMinOptBDT|AuAuCentInputMinOptBDT|auaucentinputminoptbdt) echo "auauCentInputMinOptBDT" ;;
    auauCent3BDT|AuAuCent3BDT|auaucent3bdt) echo "auauCent3BDT" ;;
    auauCent7BDT|AuAuCent7BDT|auaucent7bdt) echo "auauCent7BDT" ;;
    auauPtBinCentInputBDT|AuAuPtBinCentInputBDT|auauptbincentinputbdt) echo "auauPtBinCentInputBDT" ;;
    auauPtCent3BDT|AuAuPtCent3BDT|auauptcent3bdt) echo "auauPtCent3BDT" ;;
    auauPtCent7BDT|AuAuPtCent7BDT|auauptcent7bdt) echo "auauPtCent7BDT" ;;
    auauEtFineCentInputBDT|AuAuEtFineCentInputBDT|auauetfinecentinputbdt) echo "auauEtFineCentInputBDT" ;;
    auauEtFineCent3BDT|AuAuEtFineCent3BDT|auauetfinecent3bdt) echo "auauEtFineCent3BDT" ;;
    auauEtFineCent7BDT|AuAuEtFineCent7BDT|auauetfinecent7bdt) echo "auauEtFineCent7BDT" ;;
    auauBDTMLPStack|AuAuBDTMLPStack|auaubdtmlpstack|bdtmlpstack) echo "auauBDTMLPStack" ;;
    auauBDTSideband|AuAuBDTSideband|auaubdtsideband) echo "auauBDTSideband" ;;
    auauBDTComplement|AuAuBDTComplement|auaubdtcomplement) echo "auauBDTComplement" ;;
    auauMLPSideband|AuAuMLPSideband|auaumlpsideband) echo "auauMLPSideband" ;;
    auauMLPComplement|AuAuMLPComplement|auaumlpcomplement) echo "auauMLPComplement" ;;
    auauBDTMLPStackSideband|AuAuBDTMLPStackSideband|auaubdtmlpstacksideband) echo "auauBDTMLPStackSideband" ;;
    auauBDTMLPStackComplement|AuAuBDTMLPStackComplement|auaubdtmlpstackcomplement) echo "auauBDTMLPStackComplement" ;;
    auauCentInputMLP|AuAuCentInputMLP|auaucentinputmlp) echo "auauCentInputMLP" ;;
    auauNoCentBase3x3MLP|AuAuNoCentBase3x3MLP|auaunocentbase3x3mlp) echo "auauNoCentBase3x3MLP" ;;
    auauCentInputBase3x3MLP|AuAuCentInputBase3x3MLP|auaucentinputbase3x3mlp) echo "auauCentInputBase3x3MLP" ;;
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
    auauNoCentBDT) echo "${key}AuAuNoCentBDT" ;;
    auauCentInputBDT) echo "${key}AuAuCentInputBDT" ;;
    auauCentInput3x3BDT) echo "${key}AuAuCentInput3x3BDT" ;;
    auauCentInputBase3x3BDT) echo "${key}AuAuCentInputBase3x3BDT" ;;
    auauCentInputMinOptBDT) echo "${key}AuAuCentInputMinOptBDT" ;;
    auauCent3BDT) echo "${key}AuAuCent3BDT" ;;
    auauCent7BDT) echo "${key}AuAuCent7BDT" ;;
    auauPtBinCentInputBDT) echo "${key}AuAuPtBinCentInputBDT" ;;
    auauPtCent3BDT) echo "${key}AuAuPtCent3BDT" ;;
    auauPtCent7BDT) echo "${key}AuAuPtCent7BDT" ;;
    auauEtFineCentInputBDT) echo "${key}AuAuEtFineCentInputBDT" ;;
    auauEtFineCent3BDT) echo "${key}AuAuEtFineCent3BDT" ;;
    auauEtFineCent7BDT) echo "${key}AuAuEtFineCent7BDT" ;;
    auauBDTMLPStack) echo "${key}AuAuBDTMLPStack" ;;
    auauBDTSideband) echo "${key}AuAuBDTSideband" ;;
    auauBDTComplement) echo "${key}AuAuBDTComplement" ;;
    auauMLPSideband) echo "${key}AuAuMLPSideband" ;;
    auauMLPComplement) echo "${key}AuAuMLPComplement" ;;
    auauBDTMLPStackSideband) echo "${key}AuAuBDTMLPStackSideband" ;;
    auauBDTMLPStackComplement) echo "${key}AuAuBDTMLPStackComplement" ;;
    auauCentInputMLP) echo "${key}AuAuCentInputMLP" ;;
    auauNoCentBase3x3MLP) echo "${key}AuAuNoCentBase3x3MLP" ;;
    auauCentInputBase3x3MLP) echo "${key}AuAuCentInputBase3x3MLP" ;;
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

  if [[ -n "${SFTP_GET_CFG_MATCH:-}" ]]; then
    local match_lc
    match_lc="$(printf "%s" "$SFTP_GET_CFG_MATCH" | tr '[:upper:]' '[:lower:]')"
    local -a filtered_rows=()
    local row pre tight nonTight pre_norm tight_norm nonTight_norm selection_tag row_norm
    for row in "${photon_id_rows[@]}"; do
      IFS='|' read -r pre tight nonTight <<< "$row"
      pre_norm="$(selection_mode_normalize_for_key "preselection" "$pre")"
      tight_norm="$(selection_mode_normalize_for_key "tight" "$tight")"
      nonTight_norm="$(selection_mode_normalize_for_key "nonTight" "$nonTight")"
      selection_tag="$(selection_mode_tag "preselection" "$pre_norm")_$(selection_mode_tag "tight" "$tight_norm")_$(selection_mode_tag "nonTight" "$nonTight_norm")"
      row_norm="${pre_norm}|${tight_norm}|${nonTight_norm}|${selection_tag}"
      [[ "$(printf "%s" "$row_norm" | tr '[:upper:]' '[:lower:]')" == *"${match_lc}"* ]] && filtered_rows+=( "$row" )
    done
    (( ${#filtered_rows[@]} > 0 )) || { echo "[ERROR] SFTP_GET_CFG_MATCH='${SFTP_GET_CFG_MATCH}' matched no cfg rows in $yaml" >&2; exit 42; }
    photon_id_rows=( "${filtered_rows[@]}" )
  fi

  # Production outputs now internalize jet pT, dphi, vz, and iso/cone views.
  # The analysis-facing final ROOT files are keyed only by photon-ID selection
  # plus the real tagged UE axis for AuAu-like datasets. Keep an explicit
  # legacy mode for archived scalar outputs from older productions.
  if [[ "${SFTP_GET_LEGACY_SCALAR_TAGS:-0}" != "1" ]]; then
    if dataset_includes_uepipe_in_tag "$dataset_token"; then
      while IFS= read -r line; do uepipes+=( "$line" ); done < <(yaml_get_inline_list "$yaml" "clusterUEpipeline")
      (( ${#uepipes[@]} )) || uepipes=( "baseVariant" )
    else
      uepipes=( "" )
    fi

    local row pre tight nonTight selection_tag full_tag
    local cfg_suffix="${SFTP_GET_CFG_SUFFIX:-}"
    local tight_norm nonTight_norm pre_norm uep
    for row in "${photon_id_rows[@]}"; do
      IFS='|' read -r pre tight nonTight <<< "$row"
      pre_norm="$(selection_mode_normalize_for_key "preselection" "$pre")"
      tight_norm="$(selection_mode_normalize_for_key "tight" "$tight")"
      nonTight_norm="$(selection_mode_normalize_for_key "nonTight" "$nonTight")"
      selection_tag="$(selection_mode_tag "preselection" "$pre_norm")_$(selection_mode_tag "tight" "$tight_norm")_$(selection_mode_tag "nonTight" "$nonTight_norm")"
      for uep in "${uepipes[@]}"; do
        if [[ -n "$uep" ]]; then
          full_tag="${selection_tag}_${uep}"
        else
          full_tag="${selection_tag}"
        fi
        echo "${full_tag}${cfg_suffix}"
      done
    done | sort -u
    return 0
  fi

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
  validate_remote_path "$remote_parent" "training/latest remote parent"
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
  validate_remote_path "$train_parent" "training output parent"
  validate_remote_path "$model_parent" "training model parent"
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
    validate_remote_path "${train_parent}/${latest_train}" "training output directory"
    printf 'get -r %s/%s %s\n' "$train_parent" "$latest_train" "$latest_train"
    if [[ -n "$latest_model" ]]; then
      validate_remote_path "${model_parent}/${latest_model}" "training model directory"
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
  validate_remote_path "$remote_dir" "tight-BDT smoke remote directory"
  validate_remote_path "${REMOTE_BASE}/bdt_models/${latest}" "tight-BDT model remote directory"

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

download_auau_tight_bdt_validation() {
  local remote_dir="${1:-}"
  local report_name local_dir batch
  if [[ -z "$remote_dir" ]]; then
    echo "[ERROR] auauTightBDTValidation requires the remote model_validation_* report directory." >&2
    echo "[ERROR] Example:" >&2
    echo "  ./scripts/sftp_get_recoiljets_outputs.sh auauTightBDTValidation /sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_YYYYMMDD_HHMMSS/reports/model_validation_YYYYMMDD_HHMMSS" >&2
    exit 2
  fi
  remote_dir="${remote_dir%/}"
  report_name="${remote_dir##*/}"

  case "$remote_dir" in
    /sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_*/reports/model_validation_*) ;;
    *)
      echo "[ERROR] Refusing to pull non-validation path:" >&2
      echo "  ${remote_dir}" >&2
      echo "[ERROR] Expected /sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_*/reports/model_validation_*" >&2
      exit 2
      ;;
  esac
  validate_remote_path "$remote_dir" "AuAu tight-BDT validation remote directory"

  local_dir="${LOCAL_BASE}/dataOutput/auauTightBDTValidation"
  mkdir -p "$local_dir"
  batch="$(make_tmp_file "sftp_get_recoiljets_auau_bdt_validation")"
  cleanup_auau_bdt_validation() { rm -f "$batch"; }
  trap cleanup_auau_bdt_validation EXIT

  {
    printf 'lcd %s\n' "$local_dir"
    printf 'get -r %s %s\n' "$remote_dir" "$report_name"
  } > "$batch"

  echo
  echo "Remote host       : ${REMOTE_HOST}"
  echo "Validation remote : ${remote_dir}"
  echo "Local dir         : ${local_dir}/${report_name}"
  echo
  echo "sftp batch commands:"
  sed 's/^/  /' "$batch"
  echo
  echo "Opening interactive sftp to download the AuAu tight-BDT validation report."
  if sftp \
      -oBatchMode=no \
      -oPreferredAuthentications=publickey,password,keyboard-interactive \
      -b "$batch" \
      "$REMOTE_HOST"; then
    echo
    echo "[OK] AuAu tight-BDT validation report download complete."
    echo "Downloaded into: ${local_dir}/${report_name}"
    trap - EXIT
    rm -f "$batch"
  else
    status=$?
    echo
    echo "[ERROR] sftp download failed with exit code ${status}." >&2
    exit "$status"
  fi
}

download_auau_tight_bdt_validation_scores() {
  local remote_dir="${1:-}"
  local local_dir_arg="${2:-}"
  local local_dir batch
  if [[ -z "$remote_dir" || -z "$local_dir_arg" ]]; then
    echo "[ERROR] auauTightBDTValidationScores requires remote-report-dir and local-dir." >&2
    exit 2
  fi
  remote_dir="${remote_dir%/}"
  local_dir="${LOCAL_BASE}/${local_dir_arg#/}"

  case "$remote_dir" in
    /sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/*/validation_*|\
    /sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_*/reports/model_validation_*) ;;
    *)
      echo "[ERROR] Refusing to pull score caches from non-validation path:" >&2
      echo "  ${remote_dir}" >&2
      exit 2
      ;;
  esac
  validate_remote_path "$remote_dir" "AuAu tight-BDT validation score-cache remote directory"

  mkdir -p "$local_dir"
  batch="$(make_tmp_file "sftp_get_recoiljets_auau_bdt_validation_scores")"
  cleanup_auau_bdt_validation_scores() { rm -f "$batch"; }
  trap cleanup_auau_bdt_validation_scores EXIT

  {
    printf 'lcd %s\n' "$local_dir"
    printf 'get %s/score_caches.list score_caches.list\n' "$remote_dir"
    printf 'get -r %s/score_caches score_caches\n' "$remote_dir"
  } > "$batch"

  echo
  echo "Remote host       : ${REMOTE_HOST}"
  echo "Validation remote : ${remote_dir}"
  echo "Local dir         : ${local_dir}"
  echo
  echo "sftp batch commands:"
  sed 's/^/  /' "$batch"
  echo
  echo "Opening interactive sftp to download validation score caches only."
  if sftp \
      -oBatchMode=no \
      -oPreferredAuthentications=publickey,password,keyboard-interactive \
      -b "$batch" \
      "$REMOTE_HOST"; then
    echo
    echo "[OK] AuAu tight-BDT validation score-cache download complete."
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

download_auau_ml_diagnostic_compact() {
  local remote_dir="${1:-}"
  local local_dir="${2:-}"
  shift 2 || true
  local files=("$@")
  local batch

  if [[ -z "$remote_dir" || -z "$local_dir" || ${#files[@]} -eq 0 ]]; then
    echo "[ERROR] auauMLDiagnosticCompact requires remote-dir local-dir file..." >&2
    exit 2
  fi

  remote_dir="${remote_dir%/}"
  case "$remote_dir" in
    dataOutput/auauMLDiagnosticRuns/*)
      remote_dir="${REMOTE_BASE}/${remote_dir}"
      ;;
    "${REMOTE_BASE}"/dataOutput/auauMLDiagnosticRuns/*)
      ;;
    /sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_eiso_cone_raw_*/reports/model_validation_condor_*)
      ;;
    /sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/THE8_branchA_ladder_jet12_20*_20260527/reports/model_validation_condor_THE8_branchA_jet12_20*_scorecache_fullstat_20260527)
      ;;
    /sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/the57_models/*/validation_*)
      ;;
    /gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/THE38_tree_depth_capacity_*_d[0-9])
      ;;
    /sphenix/user/patsfan753/thesisAnalysis/bdt_models/THE38_tree_depth_capacity_*_d[0-9])
      ;;
    /sphenix/u/patsfan753/thesisAnalysis/bdt_models/THE38_tree_depth_capacity_*_d[0-9])
      ;;
    /gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/THE32_lowcalo_upstreamcut_*)
      ;;
    /sphenix/user/patsfan753/thesisAnalysis/bdt_models/THE32_lowcalo_upstreamcut_*)
      ;;
    /sphenix/u/patsfan753/scratch/thesisAnalysis/bdt_models/THE32_lowcalo_upstreamcut_*)
      ;;
    *)
      echo "[ERROR] Refusing non-compact AuAu ML diagnostic path:" >&2
      echo "  ${remote_dir}" >&2
      echo "[ERROR] Expected dataOutput/auauMLDiagnosticRuns/*, an auauTightBDT_eiso_cone_raw_* model_validation_condor_* report, a THE8_branchA_ladder compact validation report, or a THE38_tree_depth_capacity compact model registry directory." >&2
      exit 2
      ;;
  esac
  validate_remote_path "$remote_dir" "AuAu ML diagnostic remote directory"

  if [[ "$local_dir" != /* ]]; then
    local_dir="${LOCAL_BASE}/${local_dir}"
  fi

  local f
  for f in "${files[@]}"; do
    case "$f" in
      */*|*.root|*.ROOT|*.npz|*.tgz|*.tar|*.gz)
        echo "[ERROR] Refusing non-compact or path-like file argument: ${f}" >&2
        exit 2
        ;;
      *.png|*.json|*.log|*.txt|*.csv|*.md|*.yaml|*.yml) ;;
      *)
        echo "[ERROR] Refusing unsupported compact artifact extension: ${f}" >&2
        exit 2
        ;;
    esac
  done

  mkdir -p "$local_dir"
  batch="$(make_tmp_file "sftp_get_recoiljets_auau_ml_diagnostic_compact")"
  cleanup_auau_ml_diag() { rm -f "$batch"; }
  trap cleanup_auau_ml_diag EXIT

  {
    printf 'lcd %s\n' "$local_dir"
    for f in "${files[@]}"; do
      printf 'get %s/%s\n' "$remote_dir" "$f"
    done
  } > "$batch"

  echo
  echo "Remote host          : ${REMOTE_HOST}"
  echo "Diagnostic remote dir: ${remote_dir}"
  echo "Local dir            : ${local_dir}"
  echo
  echo "This downloads only compact PNG/JSON/CSV/TXT artifacts."
  echo
  echo "sftp batch commands:"
  sed 's/^/  /' "$batch"
  echo
  if sftp \
      -oBatchMode=no \
      -oPreferredAuthentications=publickey,password,keyboard-interactive \
      -b "$batch" \
      "$REMOTE_HOST"; then
    echo
    echo "[OK] AuAu ML diagnostic compact download complete."
    trap - EXIT
    rm -f "$batch"
  else
    status=$?
    echo
    echo "[ERROR] sftp download failed with exit code ${status}." >&2
    exit "$status"
  fi
}

download_stitch_diagnostics_compact() {
  local remote_dir="${1:-}"
  local local_dir="${2:-}"
  shift 2 || true
  local files=("$@")
  local batch

  if [[ -z "$remote_dir" || -z "$local_dir" || ${#files[@]} -eq 0 ]]; then
    echo "[ERROR] stitchDiagnosticsCompact requires remote-dir local-dir file..." >&2
    exit 2
  fi

  remote_dir="${remote_dir%/}"
  case "$remote_dir" in
    dataOutput/stitchDiagnostics/*)
      remote_dir="${REMOTE_BASE}/${remote_dir}"
      ;;
    "${REMOTE_BASE}"/dataOutput/stitchDiagnostics/*)
      ;;
    *)
      echo "[ERROR] Refusing non-stitchDiagnostics path:" >&2
      echo "  ${remote_dir}" >&2
      exit 2
      ;;
  esac
  validate_remote_path "$remote_dir" "stitch diagnostics remote directory"

  if [[ "$local_dir" != /* ]]; then
    local_dir="${LOCAL_BASE}/${local_dir}"
  fi

  local f
  for f in "${files[@]}"; do
    case "$f" in
      */*|*.root|*.ROOT|*.npz|*.tgz|*.tar|*.gz)
        echo "[ERROR] Refusing non-compact or path-like file argument: ${f}" >&2
        exit 2
        ;;
      *.png|*.json|*.log|*.txt|*.csv) ;;
      *)
        echo "[ERROR] Refusing unsupported compact artifact extension: ${f}" >&2
        exit 2
        ;;
    esac
  done

  mkdir -p "$local_dir"
  batch="$(make_tmp_file "sftp_get_recoiljets_stitch_diag_compact")"
  cleanup_stitch_diag() { rm -f "$batch"; }
  trap cleanup_stitch_diag EXIT

  {
    printf 'lcd %s\n' "$local_dir"
    for f in "${files[@]}"; do
      printf 'get %s/%s\n' "$remote_dir" "$f"
    done
  } > "$batch"

  echo
  echo "Remote host          : ${REMOTE_HOST}"
  echo "Stitch diagnostic dir: ${remote_dir}"
  echo "Local dir            : ${local_dir}"
  echo
  echo "This downloads only compact PNG/JSON/CSV/TXT artifacts."
  echo
  echo "sftp batch commands:"
  sed 's/^/  /' "$batch"
  echo
  if sftp \
      -oBatchMode=no \
      -oPreferredAuthentications=publickey,password,keyboard-interactive \
      -b "$batch" \
      "$REMOTE_HOST"; then
    echo
    echo "[OK] stitch diagnostics compact download complete."
    trap - EXIT
    rm -f "$batch"
  else
    status=$?
    echo
    echo "[ERROR] sftp download failed with exit code ${status}." >&2
    exit "$status"
  fi
}

download_pp_photon_ml_compact() {
  local remote_dir="${1:-}"
  local local_dir="${2:-}"
  shift 2 || true
  local files=("$@")
  local batch

  if [[ -z "$remote_dir" || -z "$local_dir" || ${#files[@]} -eq 0 ]]; then
    echo "[ERROR] ppPhotonMLCompact requires remote-dir local-dir file..." >&2
    exit 2
  fi

  remote_dir="${remote_dir%/}"
  case "$remote_dir" in
    dataOutput/ppPhotonMLPipeline/*)
      remote_dir="${REMOTE_BASE}/${remote_dir}"
      ;;
    "${REMOTE_BASE}"/dataOutput/ppPhotonMLPipeline/*)
      ;;
    *)
      echo "[ERROR] Refusing non-ppPhotonMLPipeline path:" >&2
      echo "  ${remote_dir}" >&2
      exit 2
      ;;
  esac
  validate_remote_path "$remote_dir" "pp photon-ML remote directory"

  if [[ "$local_dir" != /* ]]; then
    local_dir="${LOCAL_BASE}/${local_dir}"
  fi

  local f
  for f in "${files[@]}"; do
    case "$f" in
      */*|*.root|*.ROOT|*.npz|*.tgz|*.tar|*.gz)
        echo "[ERROR] Refusing non-compact or path-like file argument: ${f}" >&2
        exit 2
        ;;
      *.png|*.json|*.log|*.txt|*.csv) ;;
      *)
        echo "[ERROR] Refusing unsupported compact artifact extension: ${f}" >&2
        exit 2
        ;;
    esac
  done

  mkdir -p "$local_dir"
  batch="$(make_tmp_file "sftp_get_recoiljets_pp_photon_ml_compact")"
  cleanup_pp_ml_compact() { rm -f "$batch"; }
  trap cleanup_pp_ml_compact EXIT

  {
    printf 'lcd %s\n' "$local_dir"
    for f in "${files[@]}"; do
      printf 'get %s/%s\n' "$remote_dir" "$f"
    done
  } > "$batch"

  echo
  echo "Remote host          : ${REMOTE_HOST}"
  echo "ppPhotonML remote dir: ${remote_dir}"
  echo "Local dir            : ${local_dir}"
  echo
  echo "This downloads only compact PNG/JSON/CSV/TXT/LOG artifacts."
  echo
  echo "sftp batch commands:"
  sed 's/^/  /' "$batch"
  echo
  if sftp \
      -oBatchMode=no \
      -oPreferredAuthentications=publickey,password,keyboard-interactive \
      -b "$batch" \
      "$REMOTE_HOST"; then
    echo
    echo "[OK] pp photon-ML compact download complete."
    trap - EXIT
    rm -f "$batch"
  else
    status=$?
    echo
    echo "[ERROR] sftp download failed with exit code ${status}." >&2
    exit "$status"
  fi
}

download_selected_root_files() {
  local remote_dir="${1:-}"
  local local_dir="${2:-}"
  shift 2 || true
  local files=("$@")
  local batch

  if [[ -z "$remote_dir" || -z "$local_dir" || ${#files[@]} -eq 0 ]]; then
    echo "[ERROR] selectedRootFiles requires remote-dir local-dir relative-root-file..." >&2
    exit 2
  fi

  remote_dir="${remote_dir%/}"
  case "$remote_dir" in
    /sphenix/tg/tg01/bulk/jbennett/thesisAna/auau/*|\
    /sphenix/tg/tg01/bulk/jbennett/thesisAna/recoiljets/*|\
    /sphenix/tg/tg01/bulk/jbennett/thesisAnaSmoke/*|\
    /sphenix/u/patsfan753/scratch/thesisAnalysis/runs/recoiljets/current/*)
      ;;
    /sphenix/user/patsfan753/thesisAnalysis/bdt_models/THE32_lowcalo_upstreamcut_*|\
    /sphenix/u/patsfan753/scratch/thesisAnalysis/bdt_models/THE32_lowcalo_upstreamcut_*)
      ;;
    *)
      echo "[ERROR] Refusing selected ROOT pull from non-RecoilJets output path:" >&2
      echo "  ${remote_dir}" >&2
      exit 2
      ;;
  esac
  validate_remote_path "$remote_dir" "selected ROOT remote directory"

  if [[ "$local_dir" != /* ]]; then
    local_dir="${LOCAL_BASE}/${local_dir}"
  fi

  local f
  for f in "${files[@]}"; do
    validate_remote_relative_file "$f" "selected ROOT file"
    case "$f" in
      *.root|*.ROOT) ;;
      *)
        echo "[ERROR] selectedRootFiles accepts only .root file arguments: ${f}" >&2
        exit 2
        ;;
    esac
  done

  mkdir -p "$local_dir"
  batch="$(make_tmp_file "sftp_get_recoiljets_selected_roots")"
  local listfile encoded archive
  listfile="$(make_tmp_file "sftp_get_recoiljets_selected_roots_files")"
  encoded="$(make_tmp_file "sftp_get_recoiljets_selected_roots_payload.b64")"
  archive="$(make_tmp_file "sftp_get_recoiljets_selected_roots_payload.tgz")"
  cleanup_selected_roots() { rm -f "$batch" "$listfile" "$encoded" "$archive"; }
  trap cleanup_selected_roots EXIT
  printf '%s\n' "${files[@]}" > "$listfile"

  {
    printf 'lcd %s\n' "$local_dir"
    printf 'cd %s\n' "$remote_dir"
    for f in "${files[@]}"; do
      printf 'get %s %s\n' "$f" "${f##*/}"
    done
  } > "$batch"

  if [[ "${SFTP_GET_NESTED_SSH_TAR:-0}" == "1" ]]; then
    local gateway="${SFTP_GET_NESTED_SSH_GATEWAY:-patsfan753@ssh.sdcc.bnl.gov}"
    local target="${SFTP_GET_NESTED_SSH_TARGET:-sphnxuser05.sdcc.bnl.gov}"
    local remote_tar_cmd
    remote_tar_cmd="cd $(printf '%q' "$remote_dir") && printf '__RJ_TAR_BEGIN__\\n' && tar -czf - -T - | base64 && printf '\\n__RJ_TAR_END__\\n'"

    echo
    echo "Nested SSH gateway   : ${gateway}"
    echo "Nested SSH target    : ${target}"
    echo "Selected ROOT dir    : ${remote_dir}"
    echo "Local dir            : ${local_dir}"
    echo
    echo "This downloads only the explicitly listed ROOT files via nested SSH tar."
    echo
    echo "Relative ROOT files:"
    sed 's/^/  /' "$listfile"
    echo
    if SSH_AUTH_SOCK="${SSH_AUTH_SOCK:-$(launchctl getenv SSH_AUTH_SOCK 2>/dev/null || true)}" \
        ssh -q -o BatchMode=yes "$gateway" \
        "ssh -q -T -o LogLevel=ERROR -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null ${target} $(printf '%q' "$remote_tar_cmd")" \
        < "$listfile" \
        | awk '
            BEGIN { seen = 0; done = 0 }
            /^__RJ_TAR_BEGIN__$/ { seen = 1; next }
            /^__RJ_TAR_END__$/ { done = 1; exit 0 }
            seen { print }
            END { if (!seen || !done) exit 42 }
          ' > "$encoded" \
        && base64 -D -i "$encoded" -o "$archive" \
        && tar -xzf "$archive" -C "$local_dir"; then
      echo
      echo "[OK] Download complete."
      echo "Downloaded into: ${local_dir}"
      trap - EXIT
      rm -f "$batch" "$listfile" "$encoded" "$archive"
    else
      status=$?
      echo
      echo "[ERROR] nested SSH tar download failed with exit code ${status}." >&2
      echo "File list was: ${listfile}" >&2
      exit "$status"
    fi

    echo
    echo "Downloaded ${#files[@]} selected ROOT file(s)."
    return 0
  fi

  echo
  echo "Remote host          : ${REMOTE_HOST}"
  echo "Selected ROOT dir    : ${remote_dir}"
  echo "Local dir            : ${local_dir}"
  echo
  echo "This downloads only the explicitly listed ROOT files."
  echo
  echo "sftp batch commands:"
  sed 's/^/  /' "$batch"
  echo
  echo "Opening interactive sftp to download selected ROOT files."
  echo "Public-key auth is tried first; no password is needed if your SDCC key is installed."
  if sftp \
      -oBatchMode=no \
      -oPreferredAuthentications=publickey,password,keyboard-interactive \
      -b "$batch" \
      "$REMOTE_HOST"; then
    echo
    echo "[OK] Download complete."
    trap - EXIT
    rm -f "$batch"
  else
    status=$?
    echo
    echo "[ERROR] sftp download failed with exit code ${status}." >&2
    echo "[ERROR] No success confirmation was received from sftp." >&2
    exit "$status"
  fi

  echo
  echo "Downloaded ${#files[@]} selected ROOT file(s)."
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
  validate_remote_path "$remote_dir" "ML integration remote directory"
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

download_pp_photon_ml_validation() {
  local remote_dir="${1:-}"
  local local_dir="${2:-}"
  shift 2 || true
  local files=("$@")
  local batch

  if [[ -z "$remote_dir" ]]; then
    echo "[ERROR] ppPhotonMLValidation requires the remote validation directory." >&2
    echo "[ERROR] Example:" >&2
    echo "  ./scripts/sftp_get_recoiljets_outputs.sh ppPhotonMLValidation ${REMOTE_BASE}/dataOutput/ppPhotonMLPipeline/<run>/validation/<report> dataOutput/ppPhotonMLPipeline/<local-report>" >&2
    exit 2
  fi

  remote_dir="${remote_dir%/}"
  case "$remote_dir" in
    dataOutput/ppPhotonMLPipeline/*/validation/*)
      remote_dir="${REMOTE_BASE}/${remote_dir}"
      ;;
    "${REMOTE_BASE}"/dataOutput/ppPhotonMLPipeline/*/validation/*)
      ;;
    *)
      echo "[ERROR] Refusing to pull non-ppPhotonML validation path:" >&2
      echo "  ${remote_dir}" >&2
      echo "[ERROR] Expected ${REMOTE_BASE}/dataOutput/ppPhotonMLPipeline/*/validation/*" >&2
      exit 2
      ;;
  esac
  validate_remote_path "$remote_dir" "pp photon-ML validation remote directory"

  if [[ -z "$local_dir" ]]; then
    local_dir="${LOCAL_BASE}/dataOutput/ppPhotonMLPipeline/${remote_dir##*/}"
  elif [[ "$local_dir" != /* ]]; then
    local_dir="${LOCAL_BASE}/${local_dir}"
  fi

  if (( ${#files[@]} == 0 )); then
    files=(
      ppg12_fig19_equiv_pp_noCent_bdt_18_22_fullInclusive_unitnorm.png
      ppg12_fig19_equiv_pp_noCent_bdt_18_22_summary.json
      fig19_npbaudit_fix5_20260521_0920.log
    )
  fi

  local f
  for f in "${files[@]}"; do
    case "$f" in
      */*|*.root|*.ROOT)
        echo "[ERROR] Refusing non-compact or path-like file argument: ${f}" >&2
        exit 2
        ;;
      *.png|*.json|*.log|*.txt|*.csv) ;;
      *)
        echo "[ERROR] Refusing unsupported compact artifact extension: ${f}" >&2
        exit 2
        ;;
    esac
  done

  mkdir -p "$local_dir"
  batch="$(make_tmp_file "sftp_get_recoiljets_pp_photon_ml_validation")"
  cleanup_pp_photon_ml_validation() { rm -f "$batch"; }
  trap cleanup_pp_photon_ml_validation EXIT

  {
    printf 'lcd %s\n' "$local_dir"
    for f in "${files[@]}"; do
      printf 'get %s/%s\n' "$remote_dir" "$f"
    done
  } > "$batch"

  echo
  echo "Remote host          : ${REMOTE_HOST}"
  echo "Validation remote dir: ${remote_dir}"
  echo "Local dir            : ${local_dir}"
  echo
  echo "This downloads only compact PNG/JSON/LOG validation artifacts."
  echo "This will overwrite matching local files."
  read -r -p "Continue? [y/N]: " confirm
  case "$confirm" in
    y|Y|yes|YES|Yes) ;;
    *) echo "Aborted."; exit 0 ;;
  esac

  echo
  echo "sftp batch commands:"
  sed 's/^/  /' "$batch"
  echo
  echo "Opening interactive sftp to download compact pp photon-ML validation outputs."
  if sftp \
      -oBatchMode=no \
      -oPreferredAuthentications=publickey,password,keyboard-interactive \
      -b "$batch" \
      "$REMOTE_HOST"; then
    echo
    echo "[OK] pp photon-ML validation download complete."
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
      SMOKE_LABEL="isSimInclusive"; SMOKE_REMOTE_TAG="siminclusive" ;;
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
  validate_remote_path "$remote_dir" "${smoke_kind} remote directory"

  local_dir="${LOCAL_BASE}/InputFiles/pipelineSmoke/${SMOKE_LABEL}"
  mkdir -p "$local_dir"
  batch="$(make_tmp_file "sftp_get_recoiljets_${smoke_kind}")"
  cleanup_smoketest() { rm -f "$batch"; }
  trap cleanup_smoketest EXIT

  {
    printf 'lcd %s\n' "$local_dir"
    printf '%s\n' "-get -r ${remote_dir}/_pipeline_reports ${latest}_reports"
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

download_smoke_final_outputs() {
  local requested_dataset="${1:-}"
  local requested="${2:-}"
  local smoke_kind="${3:-smokeFinal}"
  local parent prefix latest remote_dir remote_final_dir local_dir local_final_dir batch
  [[ -n "$requested_dataset" ]] || { echo "[ERROR] ${smoke_kind}Latest requires a dataset, e.g. isAuAu" >&2; exit 2; }
  resolve_smoke_dataset "$requested_dataset"
  parent="${REMOTE_BASE}/outputSmoke"
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

  remote_final_dir="${remote_dir}/${SMOKE_REMOTE_TAG}"
  validate_remote_path "$remote_final_dir" "${smoke_kind} final remote directory"
  local_dir="${LOCAL_BASE}/InputFiles/pipelineSmoke/${SMOKE_LABEL}"
  local_final_dir="${local_dir}/${latest}_final"
  mkdir -p "$local_final_dir"
  batch="$(make_tmp_file "sftp_get_recoiljets_${smoke_kind}")"
  cleanup_smokefinal() { rm -f "$batch"; }
  trap cleanup_smokefinal EXIT

  {
    printf 'lcd %s\n' "$local_final_dir"
    printf 'mget %s/RecoilJets_*_ALL_*.root\n' "$remote_final_dir"
  } > "$batch"

  echo
  echo "Remote host       : ${REMOTE_HOST}"
  echo "Smoke final remote: ${remote_final_dir}"
  echo "Local final dir   : ${local_final_dir}"
  echo
  echo "This downloads only final merged smoke ROOT outputs, not per-run/per-segment trees."
  echo "Proceeding without an extra y/N prompt."

  echo
  echo "sftp batch commands:"
  sed 's/^/  /' "$batch"
  echo
  echo "Opening interactive sftp to download selected ${smoke_kind} ROOT outputs."
  if sftp \
      -oBatchMode=no \
      -oPreferredAuthentications=publickey,password,keyboard-interactive \
      -b "$batch" \
      "$REMOTE_HOST"; then
    echo
    echo "[OK] ${smoke_kind} final ROOT download complete."
    echo "Downloaded into: ${local_final_dir}"
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
  local study="${1:-scaledTriggerStudy}"
  local cfg file remote_dir local_dir batch previous_label
  case "$study" in
    scaledTriggerStudy)
      cfg="jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionReference_tightReference_nonTightReference_scaledTriggerStudy"
      previous_label="scaledTriggerStudy"
      ;;
    scaledTriggerCentStudy)
      cfg="jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionReference_tightReference_nonTightReference_scaledTriggerCentStudy_cent0_20_50_80"
      previous_label="scaledTriggerCentStudy"
      ;;
    *)
      echo "[ERROR] Unknown scaled-trigger study: ${study}" >&2
      exit 2
      ;;
  esac
  file="RecoilJets_auau_ALL_${cfg}.root"
  remote_dir="${REMOTE_BASE}/output/auau"
  validate_remote_path "$remote_dir" "scaled-trigger remote directory"
  local_dir="${LOCAL_BASE}/InputFiles/auau25"

  mkdir -p "$local_dir"
  batch="$(make_tmp_file "sftp_get_scaled_trigger_study")"
  cleanup_scaled_trigger() { rm -f "$batch"; }
  trap cleanup_scaled_trigger EXIT

  echo
  echo "Remote host : ${REMOTE_HOST}"
  echo "Remote file : ${remote_dir}/${file}"
  echo "Local file  : ${local_dir}/${file}"
  echo "Study       : ${study}"
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
        previous_dir="${local_dir}/previous/$(date +%Y%m%d_%H%M%S)_${previous_label}"
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
  echo "Opening interactive sftp to download the ${study} ROOT file."
  if sftp \
      -oBatchMode=no \
      -oPreferredAuthentications=publickey,password,keyboard-interactive \
      -b "$batch" \
      "$REMOTE_HOST"; then
    echo
    echo "[OK] ${study} download complete."
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

download_scaled_trigger_run_by_run_qa() {
  local remote_dir="${1:-}"
  local local_dir="${2:-}"
  local batch

  if [[ -z "$remote_dir" ]]; then
    echo "[ERROR] scaledTriggerRunByRunQA requires a remote artifact directory." >&2
    echo "Usage: $0 scaledTriggerRunByRunQA <remote-dir> [local-dir]" >&2
    exit 2
  fi

  if [[ -z "$local_dir" ]]; then
    local_dir="${LOCAL_BASE}/dataOutput/auau/scaledTriggerRunByRunQA/$(basename "$remote_dir")"
  fi
  validate_remote_path "$remote_dir" "scaled-trigger run-by-run QA remote directory"

  mkdir -p "$local_dir/png_good_near_unity" "$local_dir/png_problematic" "$local_dir/png_tail_only_rejected"
  if [[ -n "${SCALED_TRIGGER_GET_ALL_PNGS:-}" ]]; then
    mkdir -p "$local_dir/png_all"
  fi
  batch="$(make_tmp_file "sftp_get_scaled_trigger_run_by_run_qa")"
  cleanup_scaled_trigger_run_by_run_qa() { rm -f "$batch"; }
  trap cleanup_scaled_trigger_run_by_run_qa EXIT

  echo
  echo "Remote host : ${REMOTE_HOST}"
  echo "Remote dir  : ${remote_dir}"
  echo "Local dir   : ${local_dir}"
  echo "Study       : scaledTriggerRunByRunQA"
  echo "Scope       : summaries + curated near-unity/problematic PNG subsets"
  if [[ -n "${SCALED_TRIGGER_EXTRA_PNGS:-}" ]]; then
    echo "Extra PNGs  : ${SCALED_TRIGGER_EXTRA_PNGS}"
  fi
  if [[ -n "${SCALED_TRIGGER_GET_ALL_PNGS:-}" ]]; then
    echo "All PNGs    : enabled"
  fi
  echo

  {
    printf 'lcd %s\n' "$local_dir"
    printf 'cd %s\n' "$remote_dir"
    printf 'get run_metrics.csv run_metrics.csv\n'
    printf 'get classification_summary.txt classification_summary.txt\n'
    printf 'get near_unity_plot_table.md near_unity_plot_table.md\n'
    printf 'get problematic_plot_table.md problematic_plot_table.md\n'
    printf 'lcd %s/png_good_near_unity\n' "$local_dir"
    printf 'cd %s/png_good_near_unity\n' "$remote_dir"
    printf 'mget *.png\n'
    printf 'lcd %s/png_problematic\n' "$local_dir"
    printf 'cd %s/png_problematic\n' "$remote_dir"
    printf 'mget *.png\n'
    if [[ -n "${SCALED_TRIGGER_GET_ALL_PNGS:-}" ]]; then
      printf 'lcd %s/png_all\n' "$local_dir"
      printf 'cd %s/png_all\n' "$remote_dir"
      printf 'mget *.png\n'
    fi
    if [[ -n "${SCALED_TRIGGER_EXTRA_PNGS:-}" ]]; then
      printf 'lcd %s/png_tail_only_rejected\n' "$local_dir"
      printf 'cd %s/png_all\n' "$remote_dir"
      for png in ${SCALED_TRIGGER_EXTRA_PNGS}; do
        printf 'get %s %s\n' "$png" "$png"
      done
    fi
  } > "$batch"

  echo "sftp batch commands:"
  sed 's/^/  /' "$batch"
  echo
  echo "Opening interactive sftp to download curated scaled-trigger run-by-run QA artifacts."
  if sftp \
      -oBatchMode=no \
      -oPreferredAuthentications=publickey,password,keyboard-interactive \
      -b "$batch" \
      "$REMOTE_HOST"; then
    echo
    echo "[OK] scaledTriggerRunByRunQA download complete."
    echo "Downloaded into: ${local_dir}"
    trap - EXIT
    rm -f "$batch"
  else
    status=$?
    echo
    echo "[ERROR] scaledTriggerRunByRunQA download failed with status ${status}." >&2
    echo "Batch file was: ${batch}" >&2
    exit "$status"
  fi
}

download_auau_bdt_mlp_stack_promotion() {
  local remote_dir="${1:-}"
  local run_name local_dir batch
  if [[ -z "$remote_dir" ]]; then
    echo "[ERROR] auauBDTMLPStackPromotion requires the remote stack-promotion run directory." >&2
    echo "[ERROR] Example:" >&2
    echo "  ./scripts/sftp_get_recoiljets_outputs.sh auauBDTMLPStackPromotion /gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/bdt_mlp_stack_nn_wp80_diagnostic_YYYYMMDD_HHMMSS" >&2
    exit 2
  fi
  remote_dir="${remote_dir%/}"
  run_name="${remote_dir##*/}"

  case "$remote_dir" in
    /gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/bdt_mlp_stack_*|\
    /gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/stacked_bdt_mlp_*|\
    /sphenix/u/patsfan753/scratch/thesisAnalysis/stack_promotion_pulls/bdt_mlp_stack_*|\
    /sphenix/u/patsfan753/scratch/thesisAnalysis/stack_promotion_pulls/stacked_bdt_mlp_*) ;;
    *)
      echo "[ERROR] Refusing to pull non-stack-promotion path:" >&2
      echo "  ${remote_dir}" >&2
      echo "[ERROR] Expected a stack run under /gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/" >&2
      echo "[ERROR] or a staged stack run under /sphenix/u/patsfan753/scratch/thesisAnalysis/stack_promotion_pulls/" >&2
      exit 2
      ;;
  esac
  validate_remote_path "$remote_dir" "AuAu BDT+MLP stack promotion remote directory"

  local_dir="${LOCAL_BASE}/dataOutput/auauBDTMLPStackPromotion"
  mkdir -p "$local_dir"
  batch="$(make_tmp_file "sftp_get_recoiljets_stack_promotion")"
  cleanup_stack_promotion() { rm -f "$batch"; }
  trap cleanup_stack_promotion EXIT

  {
    printf 'lcd %s\n' "$local_dir"
    printf 'get -r %s %s\n' "$remote_dir" "$run_name"
  } > "$batch"

  echo
  echo "Remote host      : ${REMOTE_HOST}"
  echo "Stack run remote : ${remote_dir}"
  echo "Local dir        : ${local_dir}/${run_name}"
  echo
  echo "sftp batch commands:"
  sed 's/^/  /' "$batch"
  echo
  echo "Opening interactive sftp to download the AuAu BDT+MLP stack promotion run."
  if sftp \
      -oBatchMode=no \
      -oPreferredAuthentications=publickey,password,keyboard-interactive \
      -b "$batch" \
      "$REMOTE_HOST"; then
    echo
    echo "[OK] AuAu BDT+MLP stack promotion download complete."
    echo "Downloaded into: ${local_dir}/${run_name}"
    trap - EXIT
    rm -f "$batch"
  else
    status=$?
    echo
    echo "[ERROR] sftp download failed with exit code ${status}." >&2
    exit "$status"
  fi
}

download_fresh_oof_stack_compact() {
  local remote_dir="${1:-}"
  local local_dir="${2:-}"
  local run_name batch
  local files=(
    campaign_manifest.json
    partition_qa.json
    feature_contract.json
    base_model_artifacts.json
    stack_model_artifacts.json
    model_metrics.csv
    model_metrics.json
    stratified_metrics.csv
    overlay_histograms.csv
    score_correlations.json
    leakage_qa.json
  )
  if [[ -z "$remote_dir" ]]; then
    echo "[ERROR] freshOOFStackCompact requires the remote fresh OOF stack run directory." >&2
    echo "[ERROR] Example:" >&2
    echo "  ./scripts/sftp_get_recoiljets_outputs.sh freshOOFStackCompact /gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/fresh_pp_auau_bdt_mlp_oof_stack_YYYYMMDD_HHMMSS" >&2
    exit 2
  fi
  remote_dir="${remote_dir%/}"
  run_name="${remote_dir##*/}"

  case "$remote_dir" in
    /gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/fresh_pp_auau_bdt_mlp_oof_stack_*) ;;
    /gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/fresh_pp_auau_bdt_mlp_simple_holdout_stack_*) ;;
    *)
      echo "[ERROR] Refusing to pull non-fresh-stack path:" >&2
      echo "  ${remote_dir}" >&2
      echo "[ERROR] Expected a fresh_pp_auau_bdt_mlp_oof_stack_* or fresh_pp_auau_bdt_mlp_simple_holdout_stack_* run under the SDCC mlp_models area." >&2
      exit 2
      ;;
  esac
  validate_remote_path "$remote_dir" "fresh OOF stack remote directory"

  if [[ -z "$local_dir" ]]; then
    local_dir="${LOCAL_BASE}/dataOutput/fresh_pp_auau_oof_stack/${run_name}"
  elif [[ "$local_dir" != /* ]]; then
    local_dir="${LOCAL_BASE}/${local_dir}"
  fi

  mkdir -p "$local_dir/pp" "$local_dir/auau"
  batch="$(make_tmp_file "sftp_get_recoiljets_fresh_oof_stack_compact")"
  cleanup_fresh_oof_stack_compact() { rm -f "$batch"; }
  trap cleanup_fresh_oof_stack_compact EXIT

  {
    printf 'lcd %s\n' "$local_dir"
    printf 'get %s/campaign_submit_manifest.json campaign_submit_manifest.json\n' "$remote_dir"
    local domain f
    for domain in pp auau; do
      printf 'lcd %s/%s\n' "$local_dir" "$domain"
      for f in "${files[@]}"; do
        printf 'get %s/%s/%s\n' "$remote_dir" "$domain" "$f"
      done
    done
  } > "$batch"

  echo
  echo "Remote host          : ${REMOTE_HOST}"
  echo "Fresh OOF remote dir : ${remote_dir}"
  echo "Local dir            : ${local_dir}"
  echo
  echo "This downloads only compact CSV/JSON artifacts; it excludes ROOT, NPZ, and model binaries."
  echo
  echo "sftp batch commands:"
  sed 's/^/  /' "$batch"
  echo
  if sftp \
      -oBatchMode=no \
      -oPreferredAuthentications=publickey,password,keyboard-interactive \
      -b "$batch" \
      "$REMOTE_HOST"; then
    echo
    echo "[OK] fresh OOF stack compact download complete."
    echo "Downloaded into: ${local_dir}"
    trap - EXIT
    rm -f "$batch"
  else
    status=$?
    echo
    echo "[ERROR] sftp download failed with exit code ${status}." >&2
    echo "Batch file was: ${batch}" >&2
    exit "$status"
  fi
}

download_fresh_oof_stack_dag_compact() {
  local remote_dir="${1:-}"
  local local_dir="${2:-}"
  local run_name batch listfile
  local files=(
    campaign_manifest.json
    partition_qa.json
    feature_contract.json
    base_model_artifacts.json
    stack_model_artifacts.json
    model_metrics.csv
    model_metrics.json
    stratified_metrics.csv
    overlay_histograms.csv
    score_correlations.json
    leakage_qa.json
    locked_test_score_table.npz
  )
  local rel_paths=()
  if [[ -z "$remote_dir" ]]; then
    echo "[ERROR] freshOOFStackDAGCompact requires the remote clean-slate DAG run directory." >&2
    echo "[ERROR] Example:" >&2
    echo "  ./scripts/sftp_get_recoiljets_outputs.sh freshOOFStackDAGCompact /gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/fresh_pp_auau_bdt_mlp_oof_stack_dag_YYYYMMDD_HHMM_clean_dag" >&2
    exit 2
  fi
  remote_dir="${remote_dir%/}"
  run_name="${remote_dir##*/}"

  case "$remote_dir" in
    /gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/fresh_pp_auau_bdt_mlp_oof_stack_dag_*) ;;
    *)
      echo "[ERROR] Refusing to pull non-fresh-stack-DAG path:" >&2
      echo "  ${remote_dir}" >&2
      echo "[ERROR] Expected a fresh_pp_auau_bdt_mlp_oof_stack_dag_* run under the SDCC mlp_models area." >&2
      exit 2
      ;;
  esac
  validate_remote_path "$remote_dir" "fresh OOF stack DAG remote directory"

  if [[ -z "$local_dir" ]]; then
    local_dir="${LOCAL_BASE}/dataOutput/fresh_pp_auau_oof_stack/${run_name}"
  elif [[ "$local_dir" != /* ]]; then
    local_dir="${LOCAL_BASE}/${local_dir}"
  fi

  local lane domain f
  for lane in current_oof simple_holdout; do
    for domain in pp auau; do
      mkdir -p "$local_dir/$lane/$domain"
      for f in "${files[@]}"; do
        rel_paths+=("$lane/$domain/$f")
      done
    done
  done

  if [[ "${SFTP_GET_NESTED_SSH_TAR:-0}" == "1" ]]; then
    local gateway="${SFTP_GET_NESTED_SSH_GATEWAY:-patsfan753@ssh.sdcc.bnl.gov}"
    local target="${SFTP_GET_NESTED_SSH_TARGET:-sphnxuser05.sdcc.bnl.gov}"
    local remote_tar_cmd encoded archive
    listfile="$(make_tmp_file "sftp_get_recoiljets_fresh_oof_stack_dag_files")"
    encoded="$(make_tmp_file "sftp_get_recoiljets_fresh_oof_stack_dag_payload.b64")"
    archive="$(make_tmp_file "sftp_get_recoiljets_fresh_oof_stack_dag_payload.tgz")"
    cleanup_fresh_oof_stack_dag_tar() { rm -f "$listfile" "$encoded" "$archive"; }
    trap cleanup_fresh_oof_stack_dag_tar EXIT
    printf '%s\n' "${rel_paths[@]}" > "$listfile"
    remote_tar_cmd="cd $(printf '%q' "$remote_dir") && printf '__RJ_TAR_BEGIN__\\n' && tar -czf - -T - | base64 && printf '\\n__RJ_TAR_END__\\n'"

    echo
    echo "Nested SSH gateway      : ${gateway}"
    echo "Nested SSH target       : ${target}"
    echo "Fresh OOF DAG remote dir: ${remote_dir}"
    echo "Local dir               : ${local_dir}"
    echo
    echo "This downloads compact CSV/JSON artifacts plus locked_test_score_table.npz via tar over nested SSH; it excludes ROOT and model binaries."
    echo
    echo "Relative files:"
    sed 's/^/  /' "$listfile"
    echo
    if ssh -q -o BatchMode=yes "$gateway" \
        "ssh -q -T -o LogLevel=ERROR -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null ${target} $(printf '%q' "$remote_tar_cmd")" \
        < "$listfile" \
        | awk '
            BEGIN { seen = 0; done = 0 }
            /^__RJ_TAR_BEGIN__$/ { seen = 1; next }
            /^__RJ_TAR_END__$/ { done = 1; exit 0 }
            seen { print }
            END { if (!seen || !done) exit 42 }
          ' > "$encoded" \
        && base64 -D -i "$encoded" -o "$archive" \
        && tar -xzf "$archive" -C "$local_dir"; then
      echo
      echo "[OK] fresh OOF stack DAG compact download complete."
      echo "Downloaded into: ${local_dir}"
      trap - EXIT
      rm -f "$listfile" "$encoded" "$archive"
      return 0
    else
      status=$?
      echo
      echo "[ERROR] nested SSH tar download failed with exit code ${status}." >&2
      echo "File list was: ${listfile}" >&2
      exit "$status"
    fi
  fi

  batch="$(make_tmp_file "sftp_get_recoiljets_fresh_oof_stack_dag_compact")"
  cleanup_fresh_oof_stack_dag_compact() { rm -f "$batch"; }
  trap cleanup_fresh_oof_stack_dag_compact EXIT

  {
    printf 'lcd %s\n' "$local_dir"
    for lane in current_oof simple_holdout; do
      for domain in pp auau; do
        printf 'lcd %s/%s/%s\n' "$local_dir" "$lane" "$domain"
        for f in "${files[@]}"; do
          printf 'get %s/%s/%s/%s\n' "$remote_dir" "$lane" "$domain" "$f"
        done
      done
    done
  } > "$batch"

  echo
  echo "Remote host              : ${REMOTE_HOST}"
  echo "Fresh OOF DAG remote dir : ${remote_dir}"
  echo "Local dir                : ${local_dir}"
  echo
  echo "This downloads compact CSV/JSON artifacts plus locked_test_score_table.npz; it excludes ROOT and model binaries."
  echo
  echo "sftp batch commands:"
  sed 's/^/  /' "$batch"
  echo
  if sftp \
      -oBatchMode=no \
      -oPreferredAuthentications=publickey,password,keyboard-interactive \
      -b "$batch" \
      "$REMOTE_HOST"; then
    echo
    echo "[OK] fresh OOF stack DAG compact download complete."
    echo "Downloaded into: ${local_dir}"
    trap - EXIT
    rm -f "$batch"
  else
    status=$?
    echo
    echo "[ERROR] sftp download failed with exit code ${status}." >&2
    echo "Batch file was: ${batch}" >&2
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

if [[ "$dataset" == "auauTightBDTValidation" ]]; then
  download_auau_tight_bdt_validation "${2:-}"
  exit 0
fi

if [[ "$dataset" == "auauTightBDTValidationScores" ]]; then
  download_auau_tight_bdt_validation_scores "${2:-}" "${3:-}"
  exit 0
fi

if [[ "$dataset" == "auauMLDiagnosticCompact" ]]; then
  download_auau_ml_diagnostic_compact "${2:-}" "${3:-}" "${@:4}"
  exit 0
fi

if [[ "$dataset" == "stitchDiagnosticsCompact" ]]; then
  download_stitch_diagnostics_compact "${2:-}" "${3:-}" "${@:4}"
  exit 0
fi

if [[ "$dataset" == "auauBDTMLPStackPromotion" ]]; then
  download_auau_bdt_mlp_stack_promotion "${2:-}"
  exit 0
fi

if [[ "$dataset" == "freshOOFStackCompact" ]]; then
  download_fresh_oof_stack_compact "${2:-}" "${3:-}"
  exit 0
fi

if [[ "$dataset" == "freshOOFStackDAGCompact" ]]; then
  download_fresh_oof_stack_dag_compact "${2:-}" "${3:-}"
  exit 0
fi

if [[ "$dataset" == "ppPhotonMLCompact" ]]; then
  download_pp_photon_ml_compact "${2:-}" "${3:-}" "${@:4}"
  exit 0
fi

if [[ "$dataset" == "ppPhotonMLValidation" ]]; then
  download_pp_photon_ml_validation "${2:-}" "${3:-}" "${@:4}"
  exit 0
fi

if [[ "$dataset" == "selectedRootFiles" ]]; then
  download_selected_root_files "${2:-}" "${3:-}" "${@:4}"
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

if [[ "$dataset" == "smokeFinalLatest" ]]; then
  download_smoke_final_outputs "${2:-}" "" "smokeFinal"
  exit 0
fi

if [[ "$dataset" == "smokeFinal" ]]; then
  download_smoke_final_outputs "${2:-}" "${3:-}" "smokeFinal"
  exit 0
fi

if [[ "$dataset" == "scaledTriggerStudy" ]]; then
  download_scaled_trigger_study "scaledTriggerStudy"
  exit 0
fi

if [[ "$dataset" == "scaledTriggerCentStudy" ]]; then
  download_scaled_trigger_study "scaledTriggerCentStudy"
  exit 0
fi

if [[ "$dataset" == "scaledTriggerRunByRunQA" ]]; then
  download_scaled_trigger_run_by_run_qa "${2:-}" "${3:-}"
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
    remote_tag="siminclusive"
    local_subdir="InputFiles/InclusiveJetSIM"
    sample_tags=( "jet5" "jet8" "jet12" "jet20" "jet30" "jet40" )
    ;;
  *)
    echo "[ERROR] Unknown dataset: $dataset" >&2
    echo >&2
    usage >&2
    exit 2
    ;;
esac

if [[ -n "${SFTP_GET_REMOTE_DIR_OVERRIDE:-}" ]]; then
  remote_dir="${SFTP_GET_REMOTE_DIR_OVERRIDE}"
elif [[ -n "${SFTP_GET_CAMPAIGN_TAG:-}" ]]; then
  validate_remote_path "${SFTP_GET_CAMPAIGN_TAG}" "campaign tag"
  remote_dir="${REMOTE_BASE}/runs/recoiljets/current/${SFTP_GET_CAMPAIGN_TAG}/${remote_tag}"
else
  remote_dir="${REMOTE_BASE}/output/${remote_tag}"
fi
validate_remote_path "$remote_dir" "dataset output remote directory"
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
sim_combined_pull=0
if is_merge_dataset "$label" && [[ "${SFTP_GET_SIM_RAW:-0}" != "1" ]]; then
  sim_combined_pull=1
  local_dir="${LOCAL_BASE}"
fi

for cfg in "${cfg_tags[@]}"; do
  if (( sim_combined_pull )); then
    file="$(sim_combined_remote_file "$label" "$cfg")"
    local_file="$(sim_combined_local_file "$label" "$cfg")"
    remote_files+=( "$file" )
    local_files+=( "$local_file" )
    mkdir -p "$(dirname "$local_file")"
    if [[ -e "$local_file" ]]; then
      existing_files+=( "$local_file" )
    fi
  else
    for sample in "${sample_tags[@]}"; do
      file="RecoilJets_${sample}_ALL_${cfg}.root"
      remote_files+=( "$file" )
      local_files+=( "${local_dir}/${file}" )
      if [[ -e "${local_dir}/${file}" ]]; then
        existing_files+=( "${local_dir}/${file}" )
      fi
    done
  fi
done

echo
echo "Remote host : ${REMOTE_HOST}"
echo "Remote dir  : ${remote_dir}"
echo "Local dir   : ${local_dir}"
echo "Variant     : ${label}"
echo "Config YAML : $(yaml_path)"
echo
if (( sim_combined_pull )); then
  echo "Pull mode   : combined SIM products built by remote Condor finalStitch"
  echo "             set SFTP_GET_SIM_RAW=1 to pull raw per-sample secondRound files instead"
else
  echo "Pull mode   : raw/final dataset ROOT files"
fi
echo
echo "Samples:"
for sample in "${sample_tags[@]}"; do
  echo "  ${sample}"
done
echo
echo "Files to download (${#remote_files[@]}):"
for i in "${!remote_files[@]}"; do
  f="${remote_files[$i]}"
  echo "  ${remote_dir}/${f}"
  echo "    -> ${local_files[$i]}"
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
        rel="${f#${local_dir}/}"
        mkdir -p "${previous_dir}/$(dirname "$rel")"
        mv "$f" "${previous_dir}/${rel}"
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

if [[ "${SFTP_GET_NESTED_SSH_TAR:-0}" == "1" ]]; then
  nested_gateway="${SFTP_GET_NESTED_SSH_GATEWAY:-patsfan753@ssh.sdcc.bnl.gov}"
  nested_target="${SFTP_GET_NESTED_SSH_TARGET:-sphnxuser05.sdcc.bnl.gov}"
  nested_listfile="$(make_tmp_file "sftp_get_recoiljets_standard_files")"
  nested_encoded="$(make_tmp_file "sftp_get_recoiljets_standard_payload.b64")"
  nested_archive="$(make_tmp_file "sftp_get_recoiljets_standard_payload.tgz")"
  nested_extract_dir="$(mktemp -d "${LOCAL_BASE}/.tmp/sftp_get_recoiljets_standard_extract.XXXXXX")"
  cleanup_standard_nested_tar() {
    rm -f "$get_batch" "$nested_listfile" "$nested_encoded" "$nested_archive"
    rm -rf "$nested_extract_dir"
  }
  trap cleanup_standard_nested_tar EXIT

  for f in "${remote_files[@]}"; do
    validate_remote_relative_file "$f" "nested SSH ROOT pull file"
  done
  printf '%s\n' "${remote_files[@]}" > "$nested_listfile"
  nested_remote_tar_cmd="cd $(printf '%q' "$remote_dir") && printf '__RJ_TAR_BEGIN__\\n' && tar -czf - -T - | base64 && printf '\\n__RJ_TAR_END__\\n'"

  echo
  echo "Nested SSH gateway: ${nested_gateway}"
  echo "Nested SSH target : ${nested_target}"
  echo "Remote dir        : ${remote_dir}"
  echo "Local target base : ${local_dir}"
  echo
  echo "This downloads the requested dataset files through the RecoilJets transfer helper using nested SSH tar."
  echo
  echo "Relative files:"
  sed 's/^/  /' "$nested_listfile"
  echo

  if SSH_AUTH_SOCK="${SSH_AUTH_SOCK:-$(launchctl getenv SSH_AUTH_SOCK 2>/dev/null || true)}" \
      ssh -q -o BatchMode=yes "$nested_gateway" \
      "ssh -q -T -o LogLevel=ERROR -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null ${nested_target} $(printf '%q' "$nested_remote_tar_cmd")" \
      < "$nested_listfile" \
      | awk '
          BEGIN { seen = 0; done = 0 }
          /^__RJ_TAR_BEGIN__$/ { seen = 1; next }
          /^__RJ_TAR_END__$/ { done = 1; exit 0 }
          seen { print }
          END { if (!seen || !done) exit 42 }
        ' > "$nested_encoded" \
      && base64 -D -i "$nested_encoded" -o "$nested_archive" \
      && tar -xzf "$nested_archive" -C "$nested_extract_dir"; then
    for i in "${!remote_files[@]}"; do
      src="${nested_extract_dir}/${remote_files[$i]}"
      dst="${local_files[$i]}"
      if [[ ! -f "$src" ]]; then
        echo "[ERROR] Nested SSH tar did not produce expected file: ${remote_files[$i]}" >&2
        exit 43
      fi
      mkdir -p "$(dirname "$dst")"
      mv "$src" "$dst"
    done
    echo
    echo "[OK] Nested SSH tar download complete."
    trap - EXIT
    rm -f "$get_batch" "$nested_listfile" "$nested_encoded" "$nested_archive"
    rm -rf "$nested_extract_dir"
  else
    status=$?
    echo
    echo "[ERROR] nested SSH tar download failed with exit code ${status}." >&2
    exit "$status"
  fi

  echo
  echo "Downloaded ${#remote_files[@]} file(s)."

  if is_merge_dataset "$label" && (( sim_combined_pull == 0 )); then
    merge_recoiljets_sim_outputs "$label" "${cfg_tags[@]}"
  fi
  exit 0
fi

{
  printf 'lcd %s\n' "$local_dir"
  printf 'cd %s\n' "$remote_dir"
  for i in "${!remote_files[@]}"; do
    printf 'get %s %s\n' "${remote_files[$i]}" "${local_files[$i]}"
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
echo "Downloaded ${#remote_files[@]} file(s)."

if is_merge_dataset "$label" && (( sim_combined_pull == 0 )); then
  merge_recoiljets_sim_outputs "$label" "${cfg_tags[@]}"
fi
