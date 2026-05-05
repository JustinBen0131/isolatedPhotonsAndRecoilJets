#!/usr/bin/env bash
###############################################################################
# mergeRecoilJets.sh
#
# Merge RecoilJets ROOT outputs produced by RecoilJets_Condor_submit.sh.
# Supports DATA (pp, pp25, auau, oo) and SIM (isSim, isSimJet5, etc.).
#
# DISCOVERY-BASED — this script discovers cfg_tag subdirectories from the
# filesystem and NEVER reads analysis_config.yaml at merge time.
#
# ═══════════════════════════════════════════════════════════════════════════════
# DATA MERGE — 3-STEP FLOW (recommended for large datasets like AuAu)
# ═══════════════════════════════════════════════════════════════════════════════
#
# STEP 1 — Per-run partials (one Condor hadd job per run number, per cfg_tag)
#   Produces: perRun/<cfg_tag>/chunkMerge_run_<run8>.root  (e.g. 1086 files)
#   Cleans:   rm -rf on the entire perRun/<cfg_tag>/ directory (nukes old
#             chunkMerge_run_* AND sliceMerge_grp* files from prior runs).
#
#     ./mergeRecoilJets.sh condor auau
#     ./mergeRecoilJets.sh condor pp
#     ./mergeRecoilJets.sh condor pp test          # single-run smoke test
#     ./mergeRecoilJets.sh condor auau firstHalf   # first half of runs only
#
# STEP 2 (optional) — sliceRuns intermediate merge
#   Groups the per-run partials into batches of 200 and submits Condor jobs
#   that produce sliceMerge_grpNNN.root files inside perRun/<cfg_tag>/.
#   For 1086 partials → ceil(1086/200) = 6 Condor jobs per cfg_tag.
#   This step does NOT produce a final file — it only creates the slices.
#
#     ./mergeRecoilJets.sh addChunks auau condor sliceRuns
#
# STEP 3 — Final merge (one Condor job per cfg_tag)
#   Auto-detects: if sliceMerge_grp*.root files exist in perRun/<cfg_tag>/,
#   merges those (e.g. 6 files). Otherwise falls back to merging all 1086
#   chunkMerge_run_* files directly (original behavior).
#   Produces: output/<tag>/RecoilJets_<tag>_ALL_<cfg_tag>.root
#
#     ./mergeRecoilJets.sh addChunks auau condor    # Condor job
#     ./mergeRecoilJets.sh addChunks auau           # local hadd
#
# SAFETY: Step 1 nukes the entire perRun dir, so sliceMerge files can only
# exist after an explicit Step 2. The detection in Step 3 is automatic.
#
# ─── QUICK COPY/PASTE: Full AuAu 3-step flow ────────────────────────────────
#
#   ./mergeRecoilJets.sh condor auau                      # step 1
#   # (wait for all Condor jobs to finish)
#   ./mergeRecoilJets.sh addChunks auau condor sliceRuns  # step 2
#   # (wait for slice jobs to finish)
#   ./mergeRecoilJets.sh addChunks auau condor            # step 3
#
# ─── QUICK COPY/PASTE: Skip step 2 (small datasets / pp) ────────────────────
#
#   ./mergeRecoilJets.sh condor pp                        # step 1
#   ./mergeRecoilJets.sh addChunks pp condor              # step 3 directly
#
# ═══════════════════════════════════════════════════════════════════════════════
# DATA MERGE — CFG_FILTER (merge a subset of cfg_tag variants)
# ═══════════════════════════════════════════════════════════════════════════════
#
# Use CFG_FILTER env var to select which cfg_tags to process.
# Recognised values:
#   allButVariantB  — keep every tag EXCEPT those ending in _variantB
#   variantB        — keep ONLY tags ending in _variantB
#
#   CFG_FILTER=allButVariantB ./mergeRecoilJets.sh condor auau
#   CFG_FILTER=allButVariantB ./mergeRecoilJets.sh addChunks auau condor sliceRuns
#   CFG_FILTER=allButVariantB ./mergeRecoilJets.sh addChunks auau condor
#
#   CFG_FILTER=variantB ./mergeRecoilJets.sh condor auau
#   CFG_FILTER=variantB ./mergeRecoilJets.sh addChunks auau condor
#
# Omit CFG_FILTER to process all discovered cfg_tags (default).
#
# ═══════════════════════════════════════════════════════════════════════════════
# DATA MERGE — INVENTORY / DRY RUN
# ═══════════════════════════════════════════════════════════════════════════════
#
#   ./mergeRecoilJets.sh checkFileOutput pp       # report only, no side effects
#   ./mergeRecoilJets.sh checkFileOutput auau
#
#   DRYRUN=1 ./mergeRecoilJets.sh condor auau     # plan only, no submissions
#   DRYRUN=1 SKIP_TRACE=1 ./mergeRecoilJets.sh condor pp test
#
# ═══════════════════════════════════════════════════════════════════════════════
# SIM MERGE — 2-STEP FLOW
# ═══════════════════════════════════════════════════════════════════════════════
#
# STEP 1 — firstRound (group raw segment outputs into chunk-partials)
#   ./mergeRecoilJets.sh isSim firstRound
#   ./mergeRecoilJets.sh isSim firstRound groupSize 300 SAMPLE=run28_photonjet10
#   ./mergeRecoilJets.sh isSimJet5 firstRound
#   ./mergeRecoilJets.sh isSimEmbedded firstRound
#
# STEP 2 — secondRound (merge chunk-partials into one final file per cfg_tag)
#   ./mergeRecoilJets.sh isSim secondRound
#   ./mergeRecoilJets.sh isSim secondRound condor
#   ./mergeRecoilJets.sh isSim secondRound condor SAMPLE=run28_photonjet10
#
# ═══════════════════════════════════════════════════════════════════════════════
# OUTPUT LOCATIONS
# ═══════════════════════════════════════════════════════════════════════════════
#
# DATA stage-1 per-run partials:
#   output/<tag>/perRun/<cfg_tag>/chunkMerge_run_<run8>.root
#
# DATA sliceRuns intermediates (optional):
#   output/<tag>/perRun/<cfg_tag>/sliceMerge_grpNNN.root
#
# DATA final merged files:
#   output/<tag>/RecoilJets_<tag>_ALL_<cfg_tag>.root
#
# SIM firstRound chunk-partials:
#   output/<simTag>/<cfg_tag>/chunkMerge_<sampleTag>_grpNNN.root
#
# SIM secondRound final files:
#   output/<simTag>/RecoilJets_<sampleTag>_ALL_<cfg_tag>.root
#
# ═══════════════════════════════════════════════════════════════════════════════
# SAFETY / RESUME BEHAVIOR
# ═══════════════════════════════════════════════════════════════════════════════
#
#   • This script NEVER runs condor_rm (it will not kill running jobs).
#   • Stage-1 (condor mode) EXCLUDES any output ROOT file whose producing
#     RecoilJets_Condor.sh job is still in condor_q (IDLE/RUNNING/HELD/etc).
#   • DRYRUN=1   → NO deletions, NO hadd, NO condor_submit (prints the plan).
#   • SKIP_TRACE=1 → prints per-run total/busy/eligible counts.
#
# NOTES
#   • Only top-level *.root files are considered (maxdepth=1).
#   • Sorting uses: sort -V (natural ordering of chunk indices).
#   • SLICE_BATCH_SIZE=200 (tunable constant inside the addChunks section).
###############################################################################
set -euo pipefail

# DRYRUN=1  -> NO deletions, NO hadd, NO condor_submit. Only prints what would happen.
DRYRUN="${DRYRUN:-0}"

# SKIP_TRACE=1 -> print per-run (total/busy/eligible) summaries during planning/submission.
SKIP_TRACE="${SKIP_TRACE:-0}"

# Condor merge jobs should jump ahead of default-priority analysis jobs.
# Override with MERGE_CONDOR_PRIORITY=N if a different local priority is needed.
MERGE_CONDOR_PRIORITY="${MERGE_CONDOR_PRIORITY:-5}"
if [[ ! "$MERGE_CONDOR_PRIORITY" =~ ^-?[0-9]+$ ]]; then
  echo "mergeRecoilJets.sh: MERGE_CONDOR_PRIORITY must be an integer, got '${MERGE_CONDOR_PRIORITY}'" >&2
  exit 2
fi

# ---------- Pretty printing ----------
BOLD=$'\e[1m'; RED=$'\e[31m'; YEL=$'\e[33m'; GRN=$'\e[32m'; BLU=$'\e[34m'; RST=$'\e[0m'
say()  { printf "${BLU}➜${RST} %s\n" "$*"; }
warn() { printf "${YEL}⚠ %s${RST}\n" "$*" >&2; }
err()  { printf "${RED}✘ %s${RST}\n" "$*" >&2; }

# ---------- Fixed dataset roots ----------
RUN_BASE_PP="/sphenix/tg/tg01/bulk/jbennett/thesisAna/pp"
RUN_BASE_PP25="/sphenix/tg/tg01/bulk/jbennett/thesisAna/pp25"
RUN_BASE_AA="/sphenix/tg/tg01/bulk/jbennett/thesisAna/auau"
RUN_BASE_OO="/sphenix/tg/tg01/bulk/jbennett/thesisAna/oo"

# ---------- Output base (required by you) ----------
OUT_BASE="/sphenix/u/patsfan753/scratch/thesisAnalysis/output"

# Per-dataset output dirs end up as:
#   /sphenix/u/patsfan753/scratch/thesisAnalysis/output/pp
#   /sphenix/u/patsfan753/scratch/thesisAnalysis/output/auau

# ---------- Logs / temp / round files ----------
BASE="/sphenix/u/patsfan753/scratch/thesisAnalysis"
HOSTTAG="$(hostname -s 2>/dev/null || echo unknownhost)"

LOG_DIR="${BASE}/log"
OUT_DIR="${BASE}/stdout"
ERR_DIR="${BASE}/error"

# NOTE: TMP_DIR must be invocation-specific so one SIM firstRound submission
# cannot delete another submission's queued listfiles/wrapper before Condor runs.
TMP_DIR="${BASE}/tmp_recoil_merge_${HOSTTAG}_$$"

# Where splitGoldenRunList round files live:
#   ${ROUND_BASE}/<pp|auau>/goldenRuns_<pp|auau>_segmentK.txt
ROUND_BASE="${BASE}/condor_segments"

SCALED_TRIG_RUNLIST="${BASE}/dst_lists_auau/scaledEffRuns_MBD_NS_geq_2_vtx_lt_150__Pho10_12.list"
SCALED_TRIG_CONFIG_TXT="${BASE}/dst_lists_auau/scaledEffConfig_MBD_NS_geq_2_vtx_lt_150__Pho10_12.txt"

CONDOR_EXEC="${TMP_DIR}/hadd_condor.sh"   # small wrapper emitted on-the-fly

# ---------- Naming ----------
PARTIAL_PREFIX="chunkMerge_run"          # per-run partial
FINAL_PREFIX="RecoilJets"                # final combined file prefix

# ---------- Helpers ----------
usage() {
  cat <<USAGE
${BOLD}Usage:${RST}
  $0 condor <pp|pp25|auau|oo> [test|firstHalf]
  $0 addChunks <pp|pp25|auau|oo> [condor]
  $0 checkFileOutput <pp|pp25|auau|oo>

Examples:
  $0 condor pp
  $0 condor auau firstHalf
  $0 condor oo
  $0 addChunks pp
  $0 addChunks auau condor
  $0 addChunks oo
USAGE
  exit 2
}

need_cmd(){ command -v "$1" >/dev/null 2>&1 || { err "Missing command: $1"; exit 3; }; }

final_root_is_good() {
  local f="$1"
  [[ -s "$f" ]]
}

cleanup_current_tmp_dir_after_local_final() {
  [[ -n "${TMP_DIR:-}" ]] || return 0

  case "$TMP_DIR" in
    "${BASE}"/tmp_recoil_merge*)
      if [[ -d "$TMP_DIR" ]]; then
        say "Cleaning current merge tmp dir: ${TMP_DIR}"
        rm -rf -- "$TMP_DIR"
      fi
      ;;
    *)
      warn "Refusing to clean suspicious TMP_DIR: ${TMP_DIR}"
      ;;
  esac
}

scale_per_run_corrected_histograms_in_file() {
  local rootfile="$1"
  SCALED_TRIG_LAST_STATUS="not_checked"
  [[ "$TAG" == "auau" ]] || return 0
  [[ -s "$rootfile" ]] || return 0

  local base run run_num
  base="$(basename "$rootfile")"
  if [[ ! "$base" =~ ^${PARTIAL_PREFIX}_([0-9]+)\.root$ ]]; then
    return 0
  fi

  run="${BASH_REMATCH[1]}"
  run_num=$((10#$run))

  SCALED_TRIG_LAST_STATUS="not_in_runlist"
  grep -Eq "^0*${run_num}$" "$SCALED_TRIG_RUNLIST" || return 0

  local scales
  scales="$(awk -v run="$run_num" '
    $1=="CONFIG" {
      r=""; bs=""; s10=""; s12="";
      for (i=1; i<=NF; ++i) {
        split($i, a, "=");
        if (a[1] == "run") r = a[2] + 0;
        else if (a[1] == "baselineScale") bs = a[2];
        else if (a[1] == "pho10Scale") s10 = a[2];
        else if (a[1] == "pho12Scale") s12 = a[2];
      }
      if (r == run) {
        print bs, s10, s12;
        exit;
      }
    }' "$SCALED_TRIG_CONFIG_TXT")"

  if [[ -z "$scales" ]]; then
    warn "scaledTrigQA: no config entry found for run ${run_num}; leaving ${rootfile} unchanged"
    SCALED_TRIG_LAST_STATUS="no_config"
    return 0
  fi

  local baselineScale pho10Scale pho12Scale
  read -r baselineScale pho10Scale pho12Scale <<< "$scales"
  local numeric_re='^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$'
  if [[ ! "$baselineScale" =~ $numeric_re ]] || [[ ! "$pho10Scale" =~ $numeric_re ]] || [[ ! "$pho12Scale" =~ $numeric_re ]]; then
    warn "scaledTrigQA: invalid scale config for run ${run_num}; leaving ${rootfile} unchanged"
    SCALED_TRIG_LAST_STATUS="no_config"
    return 0
  fi

  SCALED_TRIG_LAST_STATUS="unknown"

  local macro="${TMP_DIR}/scaledTrigQA_scale_run_${run_num}.C"
  cat > "$macro" <<EOF
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TNamed.h>
#include <TObject.h>
#include <iostream>

{
  TFile f("${rootfile}", "UPDATE");
  if (!f.IsOpen() || f.IsZombie())
  {
    std::cerr << "[scaledTrigQA] Could not open ${rootfile} for UPDATE\\n";
    std::cout << "__SCALEDTRIGQA_STATUS__ open_failed\\n";
  }
  else if (!f.Get("scaledTrigQA_perRunCorrected_applied"))
  {
    const char* objPaths[3] = {
      "MBD_NS_geq_2_vtx_lt_150/h_maxEnergyClus_NewTriggerFilling_perRunCorrected_MBD_NS_geq_2_vtx_lt_150",
      "Photon_10/h_maxEnergyClus_NewTriggerFilling_perRunCorrected_Photon_10",
      "Photon_12/h_maxEnergyClus_NewTriggerFilling_perRunCorrected_Photon_12"
    };
    const char* dirNames[3] = {
      "MBD_NS_geq_2_vtx_lt_150",
      "Photon_10",
      "Photon_12"
    };
    const char* histNames[3] = {
      "h_maxEnergyClus_NewTriggerFilling_perRunCorrected_MBD_NS_geq_2_vtx_lt_150",
      "h_maxEnergyClus_NewTriggerFilling_perRunCorrected_Photon_10",
      "h_maxEnergyClus_NewTriggerFilling_perRunCorrected_Photon_12"
    };
    double factors[3] = {
      ${baselineScale},
      ${pho10Scale},
      ${pho12Scale}
    };

    int foundTargets = 0;
    for (int i = 0; i < 3; ++i)
    {
      if (factors[i] <= 0.0) continue;
      TH1* h = dynamic_cast<TH1*>(f.Get(objPaths[i]));
      if (h) ++foundTargets;
    }

    if (foundTargets == 0)
    {
      std::cout << "__SCALEDTRIGQA_STATUS__ skip_no_targets\\n";
    }
    else
    {
      int scaledTargets = 0;
      for (int i = 0; i < 3; ++i)
      {
        if (factors[i] <= 0.0) continue;
        TH1* h = dynamic_cast<TH1*>(f.Get(objPaths[i]));
        if (!h) continue;
        TDirectory* dir = f.GetDirectory(dirNames[i]);
        if (!dir) continue;

        h->Scale(factors[i]);
        dir->cd();
        h->Write(histNames[i], TObject::kOverwrite);
        f.cd();
        ++scaledTargets;
      }

      if (scaledTargets > 0)
      {
        TNamed marker("scaledTrigQA_perRunCorrected_applied", "1");
        marker.Write("scaledTrigQA_perRunCorrected_applied", TObject::kOverwrite);
        std::cout << "__SCALEDTRIGQA_STATUS__ scaled " << scaledTargets << "\\n";
      }
      else
      {
        std::cout << "__SCALEDTRIGQA_STATUS__ skip_no_writable_targets\\n";
      }
    }
  }
  else
  {
    std::cout << "__SCALEDTRIGQA_STATUS__ already_applied\\n";
  }
  f.Close();
}
EOF

  local root_output root_rc
  set +e
  root_output="$(root -l -b -q "$macro" 2>&1)"
  root_rc=$?
  set -e
  rm -f "$macro"

  if (( root_rc != 0 )); then
    if (( root_rc == 130 || root_rc == 143 )); then
      warn "scaledTrigQA: interrupted while scaling ${rootfile}; stopping"
      return "$root_rc"
    fi
    warn "scaledTrigQA: ROOT scaling failed for ${rootfile}"
    SCALED_TRIG_LAST_STATUS="failed"
    return 0
  fi

  case "$root_output" in
    *"__SCALEDTRIGQA_STATUS__ scaled "*)
      SCALED_TRIG_LAST_STATUS="scaled"
      ;;
    *"__SCALEDTRIGQA_STATUS__ skip_no_targets"*)
      SCALED_TRIG_LAST_STATUS="skip_no_targets"
      ;;
    *"__SCALEDTRIGQA_STATUS__ skip_no_writable_targets"*)
      SCALED_TRIG_LAST_STATUS="skip_no_writable_targets"
      ;;
    *"__SCALEDTRIGQA_STATUS__ already_applied"*)
      SCALED_TRIG_LAST_STATUS="already_applied"
      ;;
    *"__SCALEDTRIGQA_STATUS__ open_failed"*)
      SCALED_TRIG_LAST_STATUS="failed"
      ;;
    *)
      SCALED_TRIG_LAST_STATUS="unknown"
      ;;
  esac
}

scale_scaled_trig_qa_partials() {
  local partial_dir="$1"
  [[ "$TAG" == "auau" ]] || return 0
  [[ -d "$partial_dir" ]] || return 0

  if [[ ! -f "$SCALED_TRIG_RUNLIST" ]]; then
    warn "scaledTrigQA: missing run list ${SCALED_TRIG_RUNLIST}; skipping perRunCorrected scaling"
    return 0
  fi

  if [[ ! -f "$SCALED_TRIG_CONFIG_TXT" ]]; then
    warn "scaledTrigQA: missing config ${SCALED_TRIG_CONFIG_TXT}; skipping perRunCorrected scaling"
    return 0
  fi

  need_cmd root

  local -a partials=()
  mapfile -t partials < <(ls -1 "${partial_dir}/${PARTIAL_PREFIX}_"*.root 2>/dev/null | sort -V || true)
  if (( ${#partials[@]} == 0 )); then
    return 0
  fi

  say "scaledTrigQA: checking ${#partials[@]} per-run partials in ${partial_dir}"

  local f
  local checked=0 scaled=0 skipped_no_targets=0 skipped_no_writable_targets=0 already_applied=0 no_config=0 failed=0 unknown=0
  for f in "${partials[@]}"; do
    scale_per_run_corrected_histograms_in_file "$f"
    case "${SCALED_TRIG_LAST_STATUS:-unknown}" in
      not_checked|not_in_runlist)
        ;;
      scaled)
        (( checked += 1 ))
        (( scaled += 1 ))
        ;;
      skip_no_targets)
        (( checked += 1 ))
        (( skipped_no_targets += 1 ))
        ;;
      skip_no_writable_targets)
        (( checked += 1 ))
        (( skipped_no_writable_targets += 1 ))
        ;;
      already_applied)
        (( checked += 1 ))
        (( already_applied += 1 ))
        ;;
      no_config)
        (( checked += 1 ))
        (( no_config += 1 ))
        ;;
      failed)
        (( checked += 1 ))
        (( failed += 1 ))
        ;;
      *)
        (( checked += 1 ))
        (( unknown += 1 ))
        ;;
    esac
  done

  if (( checked > 0 )); then
    say "scaledTrigQA: inspected ${checked} run-listed partials: scaled=${scaled}, skipped_no_targets=${skipped_no_targets}, already_applied=${already_applied}, missing_config=${no_config}, failed=${failed}"
    if (( skipped_no_writable_targets > 0 || unknown > 0 )); then
      warn "scaledTrigQA: additional skipped/unknown statuses: skip_no_writable_targets=${skipped_no_writable_targets}, unknown=${unknown}"
    fi
  else
    say "scaledTrigQA: no per-run partials matched the scaled-trigger run list"
  fi
}

# ------------------------ Tag helpers (QA cross-check only) ----------------
trim_ws() {
  local s="$1"
  s="${s#"${s%%[![:space:]]*}"}"
  s="${s%"${s##*[![:space:]]}"}"
  printf "%s" "$s"
}

sim_is_close() { awk -v x="$1" -v y="$2" 'BEGIN{d=x-y; if(d<0)d=-d; exit !(d<1e-6)}'; }

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
  r100=$(awk -v r="$r" 'BEGIN{v=int(r*100+0.5); printf "%d", v}')
  echo "isoR${r100}"
}

sim_iso_tag() {
  local sliding="$1" fixed="$2"
  if [[ "$sliding" == "true" ]]; then
    echo "isSliding"
  else
    if [[ "$fixed" =~ ^([0-9]+)\.0+$ ]]; then
      echo "fixedIso${BASH_REMATCH[1]}GeV"
    else
      local s="$fixed"
      s="${s//./p}"
      echo "fixedIso${s}GeV"
    fi
  fi
}

# ---------------------------------------------------------------------------
# YAML-driven cfg_tag generation
#
# Goal:
#   Restrict mergeRecoilJets.sh to ONLY the cfg_tag combinations that are
#   currently present in the master analysis_config.yaml, instead of merging
#   every historical cfg_tag directory found on disk.
#
# Default YAML path:
#   ../macros/analysis_config.yaml  (relative to this script)
# Override:
#   export MERGE_CONFIG_YAML=/path/to/analysis_config.yaml
# ---------------------------------------------------------------------------
merge_yaml_path() {
  if [[ -n "${MERGE_CONFIG_YAML:-}" ]]; then
    printf "%s\n" "$MERGE_CONFIG_YAML"
    return 0
  fi
  local here
  here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  printf "%s\n" "${here}/../macros/analysis_config.yaml"
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
        if(a[i]!="") print a[i];
      }
    }'
}

selection_mode_normalize() {
  local mode
  mode="$(trim_ws "$1")"
  case "$mode" in
    ""|reference|Reference) echo "reference" ;;
    variantA|VariantA|varianta) echo "variantA" ;;
    variantB|VariantB|variantb) echo "variantB" ;;
    variantC|VariantC|variantc) echo "variantC" ;;
    variantD|VariantD|variantd) echo "variantD" ;;
    variantE|VariantE|variante) echo "variantE" ;;
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
    variantC) echo "${key}VariantC" ;;
    variantD) echo "${key}VariantD" ;;
    variantE) echo "${key}VariantE" ;;
    *)
      local first="${mode:0:1}"
      local rest="${mode:1}"
      echo "${key}${first^^}${rest}"
      ;;
  esac
}

build_iso_mode_tags_from_yaml() {
  local yaml="$1"
  local isSlidingIso isSlidingAndFixed
  isSlidingIso="$(yaml_get_scalar_bool "$yaml" "isSlidingIso" "false")"
  isSlidingAndFixed="$(yaml_get_scalar_bool "$yaml" "isSlidingAndFixed" "false")"

  mapfile -t _fixeds < <(yaml_get_inline_list "$yaml" "fixedGeV")
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

build_cfg_tags_from_yaml() {
  local dataset_token="$1"
  local yaml
  yaml="$(merge_yaml_path)"
  [[ -f "$yaml" ]] || { err "YAML not found for cfg-tag generation: $yaml"; exit 40; }

  local -a jet_pts b2bs vzs cones iso_base_tags uepipes preselection_modes tight_modes nonTight_modes
  local include_uepipe_in_tag=0
  mapfile -t jet_pts       < <(yaml_get_inline_list "$yaml" "jet_pt_min")
  mapfile -t b2bs          < <(yaml_get_inline_list "$yaml" "back_to_back_dphi_min_pi_fraction")
  mapfile -t vzs           < <(yaml_get_inline_list "$yaml" "vz_cut_cm")
  mapfile -t cones         < <(yaml_get_inline_list "$yaml" "coneR")
  mapfile -t iso_base_tags < <(build_iso_mode_tags_from_yaml "$yaml")
  mapfile -t preselection_modes < <(yaml_get_inline_list "$yaml" "preselection")
  mapfile -t tight_modes       < <(yaml_get_inline_list "$yaml" "tight")
  mapfile -t nonTight_modes    < <(yaml_get_inline_list "$yaml" "nonTight")

  if (( ${#jet_pts[@]} == 0 )); then jet_pts=( "5.0" ); fi
  if (( ${#b2bs[@]} == 0 )); then b2bs=( "0.875" ); fi
  if (( ${#vzs[@]} == 0 )); then vzs=( "30.0" ); fi
  if (( ${#cones[@]} == 0 )); then cones=( "0.30" ); fi
  if (( ${#iso_base_tags[@]} == 0 )); then iso_base_tags=( "fixedIso2GeV" ); fi
  if (( ${#preselection_modes[@]} == 0 )); then preselection_modes=( "reference" ); fi
  if (( ${#tight_modes[@]} == 0 )); then tight_modes=( "reference" ); fi
  if (( ${#nonTight_modes[@]} == 0 )); then nonTight_modes=( "reference" ); fi

  case "$dataset_token" in
    auau|oo|isSimEmbedded|issimembedded|simembedded|SIMEMBEDDED|isSimEmbeddedInclusive|issimembeddedinclusive|simembeddedinclusive|SIMEMBEDDEDINCLUSIVE)
      mapfile -t uepipes < <(yaml_get_inline_list "$yaml" "clusterUEpipeline")
      (( ${#uepipes[@]} > 0 )) || uepipes=( "noSub" )
      include_uepipe_in_tag=1
      ;;
    *)
      uepipes=( "noSub" )
      include_uepipe_in_tag=0
      ;;
  esac

  local pt frac vz cone iso uep pre tight nonTight tag selection_tag full_tag
  local cfg_suffix="${MERGE_CFG_SUFFIX:-}"
  local tight_norm nonTight_norm
  for pt in "${jet_pts[@]}"; do
    for frac in "${b2bs[@]}"; do
      for vz in "${vzs[@]}"; do
        for cone in "${cones[@]}"; do
          for iso in "${iso_base_tags[@]}"; do
            for pre in "${preselection_modes[@]}"; do
              for tight in "${tight_modes[@]}"; do
                for nonTight in "${nonTight_modes[@]}"; do
                  tight_norm="$(selection_mode_normalize "$tight")"
                  nonTight_norm="$(selection_mode_normalize "$nonTight")"
                  if [[ "$tight_norm" == "reference" && "$nonTight_norm" != "reference" ]]; then
                    continue
                  fi
                  selection_tag="$(selection_mode_tag "preselection" "$pre")_$(selection_mode_tag "tight" "$tight")_$(selection_mode_tag "nonTight" "$nonTight")"
                  tag="jetMinPt$(sim_pt_tag "$pt")_$(sim_b2b_dir_tag "$frac")_$(sim_vz_tag "$vz")_$(sim_cone_tag "$cone")_${iso}"
                  for uep in "${uepipes[@]}"; do
                    if (( include_uepipe_in_tag )); then
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
      done
    done
  done | sort -u
}

# ---------------------------------------------------------------------------
# Extract the analysis cfg tag directly from a ROOT file's stamped YAML.
# Used ONLY for QA cross-checking (directory name vs stamped tag).
# Every segment ROOT file contains a TObjString "analysis_config_yaml" with
# the per-job override YAML (all scalar values). This function reads it,
# parses the scalar keys, and prints the same cfg tag string the submission
# scripts would produce. Returns empty string on any failure.
# ---------------------------------------------------------------------------
extract_cfg_tag_from_root() {
  local rootfile="$1"
  [[ -f "$rootfile" ]] || return 1

  # Write a temp ROOT macro (quoted heredoc prevents bash expansion of C++ $variables)
  local _macro="${TMP_DIR}/_extract_cfg_$$.C"
  cat > "$_macro" <<'ENDMACRO'
#include <TFile.h>
#include <TObjString.h>
#include <iostream>
void _extract_cfg(const char* fname) {
  TFile* f = TFile::Open(fname, "READ");
  if (!f || f->IsZombie()) { std::cout << "__FAIL__" << std::endl; return; }
  TObjString* obj = dynamic_cast<TObjString*>(f->Get("analysis_config_yaml"));
  if (!obj) { std::cout << "__FAIL__" << std::endl; f->Close(); return; }
  std::cout << "__YAMLBEGIN__" << std::endl;
  std::cout << obj->GetString().Data() << std::endl;
  std::cout << "__YAMLEND__" << std::endl;
  f->Close();
}
ENDMACRO

  local _raw
  _raw="$(root -b -l -q "${_macro}(\"${rootfile}\")" 2>/dev/null || true)"
  rm -f "$_macro"

  local _yaml
  _yaml="$(echo "$_raw" | sed -n '/__YAMLBEGIN__/,/__YAMLEND__/{ /__YAML/d; p; }')"
  [[ -n "$_yaml" ]] || return 1

  local _tmpyaml="${TMP_DIR}/_extract_yaml_$$.yaml"
  echo "$_yaml" > "$_tmpyaml"

  # Parse scalar keys from the per-job override YAML
  local _pt _frac _vz _cone _sliding _fixed _preselection _tight _nonTight

  _pt="$(grep -E '^[[:space:]]*jet_pt_min:' "$_tmpyaml" | head -1 || true)"
  _pt="${_pt#*:}"; _pt="${_pt%%#*}"; _pt="$(trim_ws "$_pt")"
  [[ -n "$_pt" ]] || _pt="5.0"

  _frac="$(grep -E '^[[:space:]]*back_to_back_dphi_min_pi_fraction:' "$_tmpyaml" | head -1 || true)"
  _frac="${_frac#*:}"; _frac="${_frac%%#*}"; _frac="$(trim_ws "$_frac")"
  [[ -n "$_frac" ]] || _frac="0.875"

  _vz="$(grep -E '^[[:space:]]*vz_cut_cm:' "$_tmpyaml" | head -1 || true)"
  _vz="${_vz#*:}"; _vz="${_vz%%#*}"; _vz="$(trim_ws "$_vz")"
  [[ -n "$_vz" ]] || _vz="30.0"

  _cone="$(grep -E '^[[:space:]]*coneR:' "$_tmpyaml" | head -1 || true)"
  _cone="${_cone#*:}"; _cone="${_cone%%#*}"; _cone="$(trim_ws "$_cone")"
  [[ -n "$_cone" ]] || _cone="0.30"

  _sliding="$(grep -E '^[[:space:]]*isSlidingIso:' "$_tmpyaml" | head -1 || true)"
  _sliding="${_sliding#*:}"; _sliding="${_sliding%%#*}"; _sliding="$(trim_ws "$_sliding")"
  [[ -n "$_sliding" ]] || _sliding="false"

  _fixed="$(grep -E '^[[:space:]]*fixedGeV:' "$_tmpyaml" | head -1 || true)"
  _fixed="${_fixed#*:}"; _fixed="${_fixed%%#*}"; _fixed="$(trim_ws "$_fixed")"
  [[ -n "$_fixed" ]] || _fixed="2.0"

  _preselection="$(grep -E '^[[:space:]]*preselection:' "$_tmpyaml" | head -1 || true)"
  _preselection="${_preselection#*:}"; _preselection="${_preselection%%#*}"; _preselection="$(trim_ws "$_preselection")"
  [[ -n "$_preselection" ]] || _preselection="reference"

  _tight="$(grep -E '^[[:space:]]*tight:' "$_tmpyaml" | head -1 || true)"
  _tight="${_tight#*:}"; _tight="${_tight%%#*}"; _tight="$(trim_ws "$_tight")"
  [[ -n "$_tight" ]] || _tight="reference"

  _nonTight="$(grep -E '^[[:space:]]*nonTight:' "$_tmpyaml" | head -1 || true)"
  _nonTight="${_nonTight#*:}"; _nonTight="${_nonTight%%#*}"; _nonTight="$(trim_ws "$_nonTight")"
  [[ -n "$_nonTight" ]] || _nonTight="reference"

  local _uepipe
  _uepipe="$(grep -E '^[[:space:]]*clusterUEpipeline:' "$_tmpyaml" | head -1 || true)"
  _uepipe="${_uepipe#*:}"; _uepipe="${_uepipe%%#*}"; _uepipe="$(trim_ws "$_uepipe")"
  case "$_uepipe" in
    true|1) _uepipe="variantA" ;;
    false|0|"") _uepipe="noSub" ;;
  esac

  rm -f "$_tmpyaml"

  local _selection_tag
  _selection_tag="$(selection_mode_tag "preselection" "$_preselection")_$(selection_mode_tag "tight" "$_tight")_$(selection_mode_tag "nonTight" "$_nonTight")"

  local _tag="jetMinPt$(sim_pt_tag "$_pt")_$(sim_b2b_dir_tag "$_frac")_$(sim_vz_tag "$_vz")_$(sim_cone_tag "$_cone")_$(sim_iso_tag "$_sliding" "$_fixed")"
  if [[ "$_uepipe" != "noSub" ]]; then
    _tag="${_tag}_${_uepipe}_${_selection_tag}"
  else
    _tag="${_tag}_${_selection_tag}"
  fi
  echo "$_tag"
}

to_tag() {
  case "${1:-}" in
    pp|PP|isPP|PP_DATA|pp_data)   echo "pp" ;;
    pp25|PP25|isPPrun25|pprun25|PP_RUN25|pp_run25) echo "pp25" ;;
    auau|AA|isAuAu|AuAu|aa|AA_DATA|auau_data) echo "auau" ;;
    oo|OO|isOO|OO_DATA|oo_data) echo "oo" ;;
    *) err "Dataset must be 'pp', 'pp25', 'auau', or 'oo'"; exit 4 ;;
  esac
}

resolve_dataset() {
  TAG="$(to_tag "${1:-}")"
  case "$TAG" in
    pp)   RUN_BASE="$RUN_BASE_PP" ;;
    pp25) RUN_BASE="$RUN_BASE_PP25" ;;
    auau) RUN_BASE="$RUN_BASE_AA" ;;
    oo)   RUN_BASE="$RUN_BASE_OO" ;;
  esac

  DEST_DIR="${OUT_BASE}/${TAG}"      # where partials and final live
  PERRUN_DIR="${OUT_BASE}/${TAG}/perRun"  # per-run partials (stage-1 merge output)
  ROUND_DIR="${ROUND_BASE}/${TAG}"   # where splitGoldenRunList round files live

  mkdir -p "$DEST_DIR" "$PERRUN_DIR" "$ROUND_DIR" "$LOG_DIR" "$OUT_DIR" "$ERR_DIR" "$TMP_DIR"
}

emit_hadd_wrapper() {
  local exe="$1"
  cat > "$exe" <<'EOS'
#!/usr/bin/env bash
set -eo pipefail
set +u
export USER="$(id -un)"; export LOGNAME="$USER"; export HOME="/sphenix/u/$USER"
MYINSTALL="/sphenix/u/$USER/thesisAnalysis/install"
source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh "$MYINSTALL" || true
set -u
if [[ $# -ne 2 ]]; then
  echo "[ERROR] Usage: $0 <listfile> <outroot>" >&2
  exit 1
fi
LIST="$1"; OUT="$2"
echo "[hadd_condor] inputs=$(wc -l < "$LIST")  ->  $OUT"
hadd -v 3 -f "$OUT" @"$LIST"
EOS
  chmod +x "$exe"
}

# Collect run directories that actually contain at least one *.root file
discover_runs() {
  local base="$1"
  mapfile -t RUNS < <(find "$base" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort -V)
  # Filter to only runs with at least one .root at top-level
  local valid=()
  for r in "${RUNS[@]}"; do
    if compgen -G "${base}/${r}/*.root" > /dev/null; then
      valid+=( "$r" )
    fi
  done
  RUNS=( "${valid[@]}" )
}

# Normalize a run token to 8 digits (handles "47289" or "00047289")
pad_run8() {
  local x="${1:-}"
  x="${x%%#*}"
  x="${x//[[:space:]]/}"
  [[ -n "$x" ]] || return 1
  [[ "$x" =~ ^[0-9]+$ ]] || return 1
  printf "%08d" "$((10#$x))"
}

# Optional verbose preview of a run list (smoke-test friendly)
print_runs_preview() {
  local runs=( "$@" )
  local n="${#runs[@]}"
  (( n == 0 )) && return 0

  local head_n=5
  local tail_n=5

  if (( ${PRINT_ALL_RUNS:-0} )); then
    say "  All runs (${n}):"
    printf "%s\n" "${runs[@]}" | sed 's/^/    /'
    return 0
  fi

  if (( n <= head_n + tail_n )); then
    say "  Runs (${n}): ${runs[*]}"
    return 0
  fi

  say "  First ${head_n} runs: ${runs[*]:0:${head_n}}"
  say "  Last  ${tail_n} runs: ${runs[*]:$((n-tail_n)):${tail_n}}"
  say "  (Set PRINT_ALL_RUNS=1 to print every run)"
}

# Load RUNS from splitGoldenRunList segment files:
#   ${ROUND_DIR}/goldenRuns_${TAG}_segmentK.txt
# and restrict to those runs (so a given schedd only merges “its” segments).
load_runs_from_segments() {
  local segs=( "$@" )
  (( ${#segs[@]} > 0 )) || { err "No segment indices provided to load_runs_from_segments()"; exit 30; }

  [[ -d "${ROUND_DIR:-}" ]] || { err "ROUND_DIR not found: ${ROUND_DIR:-<unset>}"; exit 31; }

  local tmp="${TMP_DIR}/runs_${TAG}_segments_${HOSTTAG}.txt"
  : > "$tmp"

  say "Run selection mode: ${BOLD}fromSplitRunList${RST} (host=${HOSTTAG})"
  say "  Round dir: ${ROUND_DIR}"

  for s in "${segs[@]}"; do
    [[ "$s" =~ ^[0-9]+$ ]] || { err "Bad segment index '$s' (must be integer)"; exit 32; }
    local f="${ROUND_DIR}/goldenRuns_${TAG}_segment${s}.txt"
    [[ -s "$f" ]] || { err "Segment file missing/empty: $f"; exit 33; }

    local nr
    nr=$(grep -E -v '^[[:space:]]*($|#)' "$f" | wc -l | awk '{print $1}')
    say "  Segment ${s}: ${nr} runs  ->  $(basename "$f")"

    while IFS= read -r line; do
      [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]] && continue
      r8="$(pad_run8 "$line" || true)"
      [[ -n "${r8:-}" ]] || continue
      printf "%s\n" "$r8" >> "$tmp"
    done < "$f"
  done

  # Unique + sort
  mapfile -t RUNS < <(sort -u -V "$tmp")

  # Filter to runs that actually exist under RUN_BASE (directory exists),
  # to avoid trying to merge nonsense.
  local filtered=()
  for r in "${RUNS[@]}"; do
    if [[ -d "${RUN_BASE}/${r}" ]]; then
      filtered+=( "$r" )
    else
      warn "[segments] run ${r} listed in segments but no directory under RUN_BASE: ${RUN_BASE}/${r}"
    fi
  done
  RUNS=( "${filtered[@]}" )

  say "  Runs loaded (after filtering to existing dirs): ${#RUNS[@]}"
  print_runs_preview "${RUNS[@]}"
}

# -----------------------------------------------------------------------------
# Active-job skiplist
#   Build ONCE per script run: list of output ROOT files that correspond to
#   RecoilJets_Condor.sh jobs currently still in condor_q (IDLE/RUNNING/HELD/etc).
#   Those files are excluded from per-run hadd inputs.
# -----------------------------------------------------------------------------
SKIP_BUILT=0
SKIP_FILE=""

build_active_skiplist() {
  (( SKIP_BUILT )) && return 0

  SKIP_FILE="${TMP_DIR}/skip_active_RecoilJets_${TAG}.txt"
  : > "$SKIP_FILE"

  local want="isPP"
  local outds="isPP"
  if [[ "$TAG" == "auau" ]]; then
    want="isAuAu"
    outds="isAuAu"
  fi
  if [[ "$TAG" == "pp25" ]]; then
    want="isPPrun25"
    outds="isPP"
  fi
  if [[ "$TAG" == "oo" ]]; then
    want="isOO"
    outds="isOO"
  fi

  if command -v condor_q >/dev/null 2>&1; then
    (
      set +e +o pipefail
      condor_q "${USER:-$(id -un)}" \
        -constraint '(regexp("RecoilJets_Condor.sh",Cmd) || regexp("RecoilJets_Condor_AuAu.sh",Cmd)) && (JobStatus==1 || JobStatus==2 || JobStatus==5 || JobStatus==6 || JobStatus==7)' \
        -af Args 2>/dev/null |
      awk -v want="$want" -v outds="$outds" '
        # Args format (from your submit):
        #   run8  chunkList  isPP|isPPrun25|isAuAu|isOO  Cluster  0  grpIdx  NONE  destBase
        ($3 == want) {
          run  = $1
          lst  = $2
          dest = $NF
          if (run=="" || lst=="" || dest=="") next

          # chunk tag = basename(list) without .list
          n = lst
          sub(/^.*\//,"",n)
          sub(/\.list$/,"",n)

          # output ROOT path matches RecoilJets_Condor.sh / RecoilJets_Condor_AuAu.sh naming:
          #   destBase/run8/RecoilJets_<outds>_<chunkTag>.root
          printf "%s/%s/RecoilJets_%s_%s.root\n", dest, run, outds, n
        }
      ' | sort -u > "$SKIP_FILE"
    ) || true
  fi

  SKIP_BUILT=1

  local nskip
  nskip=$(wc -l < "$SKIP_FILE" | awk '{print $1}')
  local want2="$want"

  if (( nskip > 0 )); then
    printf "${YEL}⚠${RST} [skiplist] %d active %s outputs will be excluded (from condor_q)\n" "$nskip" "$want2" >&2
    if (( SKIP_TRACE )); then
      printf "${DIM:-}${YEL}⚠${RST} [skiplist] first 5:\n" >&2
      head -n 5 "$SKIP_FILE" >&2
    fi
  else
    printf "${GRN}✔${RST} [skiplist] No active %s RecoilJets jobs found in condor_q\n" "$want2" >&2
  fi

  return 0
}

# Build a per-run list of files (sorted) to be merged
# NEW behavior:
#   - Excludes any outputs whose producing RecoilJets_Condor.sh job is still in condor_q.
#   - Supports DRYRUN=1 / SKIP_TRACE=1 summaries without polluting stdout (list path only).
make_run_list() {
  local run8="$1"
  local srcdir="${RUN_BASE}/${run8}"
  local list="${TMP_DIR}/run_${_cfg:-flat}_${run8}.txt"
  local all="${list}.00_all"

  build_active_skiplist

  # 1) All ROOTs on disk for this run (exclude LOCAL test outputs)
  find "$srcdir" -maxdepth 1 -type f -name "*.root" -not -name "*_LOCAL_*" | sort -V > "$all"
  [[ -s "$all" ]] || { rm -f "$all" 2>/dev/null || true; return 1; }

  local total busy_present eligible busy_inq
  total=$(wc -l < "$all" | awk '{print $1}')

  # 2) Count how many ACTIVE-job outputs are already present on disk (intersection)
  busy_present=0
  if [[ -s "$SKIP_FILE" ]]; then
    busy_present=$(grep -F -x -f "$SKIP_FILE" "$all" 2>/dev/null | wc -l | awk '{print $1}')
  fi

  # 3) Remove active-job outputs from the input list
  if [[ -s "$SKIP_FILE" ]]; then
    grep -F -v -f "$SKIP_FILE" "$all" > "$list" || true
    rm -f "$all" 2>/dev/null || true
  else
    mv "$all" "$list"
  fi

  [[ -s "$list" ]] || { rm -f "$list" 2>/dev/null || true; return 1; }

  eligible=$(wc -l < "$list" | awk '{print $1}')

  # Optional: how many outputs for this run are still in condor_q (may not exist yet)
  busy_inq=0
  if [[ -s "$SKIP_FILE" ]]; then
    busy_inq=$(awk -v pre="${RUN_BASE}/${run8}/" 'index($0,pre)==1{c++} END{print c+0}' "$SKIP_FILE")
  fi

  if (( DRYRUN )) || (( SKIP_TRACE )); then
    printf "${BLU}➜${RST} [run ${run8}] total=%d  busyInQ=%d  busyPresent=%d  eligible=%d\n" \
      "$total" "$busy_inq" "$busy_present" "$eligible" >&2
  fi

  echo "$list"
}

# ---------- Parse CLI ----------
[[ $# -ge 2 ]] || usage

# ============================================================
# SIM mode (DISCOVERY-BASED — does NOT read YAML):
#   ./mergeRecoilJets.sh isSim firstRound  [groupSize N] [SAMPLE=run28_photonjet10]
#   ./mergeRecoilJets.sh isSim secondRound [condor]      [SAMPLE=run28_photonjet10]
#
# cfg_tag directories are DISCOVERED from the output tree written by
# RecoilJets_Condor_submit.sh. The YAML is only used at submission time;
# editing it between submission and merging has NO effect on the merge.
#
# Output (flat, one dir per dataset type):
#   /sphenix/u/patsfan753/scratch/thesisAnalysis/output/<simTag>/
#     RecoilJets_<sampleTag>_ALL_<cfg_tag>.root
# ============================================================
if [[ "${1}" =~ ^(isSim|sim|SIM|isSimJet5|isSimjet5|simjet5|SIMJET5|isSimMB|simmb|SIMMB|isSimEmbedded|issimembedded|simembedded|SIMEMBEDDED|isSimEmbeddedInclusive|issimembeddedinclusive|simembeddedinclusive|SIMEMBEDDEDINCLUSIVE)$ ]]; then
  SIM_DATASET_TOKEN="${1}"
  SIM_ACTION="${2:-}"
  shift 2

  # Defaults (overrideable); variant-specific below
  SIM_SAMPLE="run28_photonjet10"
  SIM_SAMPLE_EXPLICIT=0
  SIM_GROUP_SIZE="300"          # number of ROOT files per firstRound hadd (default = 300)
  SIM_PREFER_CONDOR=false       # only used for secondRound

  # Defaults for optional behavior flags
  SIM_FIRSTROUND_LOCAL=false

  # Parse optional tokens
  while [[ $# -gt 0 ]]; do
    case "$1" in
      groupSize)
        [[ $# -ge 2 ]] || { err "groupSize requires a value"; exit 2; }
        SIM_GROUP_SIZE="$2"
        shift 2
        ;;
      SAMPLE=*)
        SIM_SAMPLE="${1#SAMPLE=}"
        SIM_SAMPLE_EXPLICIT=1
        shift
        ;;
      condor)
        SIM_PREFER_CONDOR=true
        shift
        ;;
      local)
        SIM_FIRSTROUND_LOCAL=true
        shift
        ;;
      *)
        # ignore unknown tokens to stay backward-compatible
        shift
        ;;
    esac
  done

  [[ "$SIM_GROUP_SIZE" =~ ^[0-9]+$ ]] || { err "groupSize must be an integer"; exit 2; }
  (( SIM_GROUP_SIZE > 0 )) || { err "groupSize must be > 0"; exit 2; }

  # Input sim outputs live here (matches Condor output layout per variant)
  case "$SIM_DATASET_TOKEN" in
    isSimJet5|isSimjet5|simjet5|SIMJET5)
      SIM_INPUT_BASE="/sphenix/tg/tg01/bulk/jbennett/thesisAna/simjet5"
      SIM_OUTPUT_TAG="simjet5"
      [[ "$SIM_SAMPLE_EXPLICIT" -eq 0 ]] && SIM_SAMPLE="run28_jet5"
      ;;
    isSimMB|simmb|SIMMB)
      SIM_INPUT_BASE="/sphenix/tg/tg01/bulk/jbennett/thesisAna/simmb"
      SIM_OUTPUT_TAG="simmb"
      [[ "$SIM_SAMPLE_EXPLICIT" -eq 0 ]] && SIM_SAMPLE="run28_detroit"
      ;;
    isSimEmbedded|issimembedded|simembedded|SIMEMBEDDED)
      SIM_INPUT_BASE="/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded"
      SIM_OUTPUT_TAG="simembedded"
      [[ "$SIM_SAMPLE_EXPLICIT" -eq 0 ]] && SIM_SAMPLE="run28_embeddedPhoton20"
      ;;
    isSimEmbeddedInclusive|issimembeddedinclusive|simembeddedinclusive|SIMEMBEDDEDINCLUSIVE)
      SIM_INPUT_BASE="/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded"
      SIM_OUTPUT_TAG="simembeddedinclusive"
      [[ "$SIM_SAMPLE_EXPLICIT" -eq 0 ]] && SIM_SAMPLE="run28_embeddedJet20"
      ;;
    *)
      SIM_INPUT_BASE="/sphenix/tg/tg01/bulk/jbennett/thesisAna/sim"
      SIM_OUTPUT_TAG="sim"
      ;;
  esac

  # Embedded firstRound hadd jobs can exceed the generic 2 GB request when
  # merging larger Au+Au-style ROOT outputs. Keep pp/MB/jet SIM unchanged.
  SIM_FIRSTROUND_REQUEST_MEMORY="2GB"
  case "$SIM_DATASET_TOKEN" in
    isSimEmbedded|issimembedded|simembedded|SIMEMBEDDED|isSimEmbeddedInclusive|issimembeddedinclusive|simembeddedinclusive|SIMEMBEDDEDINCLUSIVE)
      SIM_FIRSTROUND_REQUEST_MEMORY="4GB"
      ;;
  esac

  # Build sample list (same defaults as submit script, or explicit SAMPLE=)
  samples=()
  if [[ "${SIM_SAMPLE_EXPLICIT:-0}" -eq 0 ]]; then
    case "$SIM_DATASET_TOKEN" in
      isSimJet5|isSimjet5|simjet5|SIMJET5) samples=( "run28_jet5" ) ;;
      isSimMB|simmb|SIMMB)       samples=( "run28_detroit" ) ;;
      isSimEmbedded|issimembedded|simembedded|SIMEMBEDDED) samples=( "run28_embeddedPhoton12" "run28_embeddedPhoton20" ) ;;
      isSimEmbeddedInclusive|issimembeddedinclusive|simembeddedinclusive|SIMEMBEDDEDINCLUSIVE) samples=( "run28_embeddedJet10" "run28_embeddedJet20" ) ;;
      *)                         samples=( "run28_photonjet5" "run28_photonjet10" "run28_photonjet20" ) ;;
    esac
  else
    samples=( "${SIM_SAMPLE}" )
  fi

  # ---------------------------------------------------------------
  # DISCOVER cfg_tag directories from the filesystem.
  #
  # IMPORTANT:
  #   isSimEmbedded and isSimEmbeddedInclusive currently share the same TG
  #   input base:
  #     /sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded
  #
  #   So after the raw cfg_tag discovery, we must FILTER to only cfg_tag
  #   directories that actually contain one of the requested sample
  #   subdirectories for THIS dataset token. This prevents
  #   isSimEmbeddedInclusive from trying to merge the old photon-embedded
  #   output tree, and vice versa.
  # ---------------------------------------------------------------
  FLAT_OUT_DIR="${OUT_BASE}/${SIM_OUTPUT_TAG}"
  mkdir -p "$FLAT_OUT_DIR" "$LOG_DIR" "$OUT_DIR" "$ERR_DIR" "$TMP_DIR"

  if [[ "$SIM_ACTION" == "firstRound" ]]; then
    _discover_base="$SIM_INPUT_BASE"
  elif [[ "$SIM_ACTION" == "secondRound" ]]; then
    # Keep discovery anchored to the TG input base, but filter cfg_tags below
    # by the requested sample directories for this dataset token.
    _discover_base="$SIM_INPUT_BASE"
  else
    err "Unknown isSim action '${SIM_ACTION}'. Allowed: firstRound | secondRound"
    exit 2
  fi

  mapfile -t _ALL_SIM_CFG_TAGS < <(
    find "$_discover_base" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort -V
  )

  mapfile -t _YAML_SIM_CFG_TAGS < <(build_cfg_tags_from_yaml "$SIM_DATASET_TOKEN")

  SIM_CFG_TAGS=()
  for _cfg in "${_ALL_SIM_CFG_TAGS[@]}"; do
    _keep_cfg=0

    # First require that the cfg_tag is part of the CURRENT master YAML matrix
    for _y in "${_YAML_SIM_CFG_TAGS[@]}"; do
      if [[ "$_cfg" == "$_y" ]]; then
        _keep_cfg=1
        break
      fi
    done
    (( _keep_cfg )) || continue

    # Then require that this cfg_tag actually contains one of the requested samples
    _keep_cfg=0
    for _samp in "${samples[@]}"; do
      if [[ -d "${_discover_base}/${_cfg}/${_samp}" ]]; then
        _keep_cfg=1
        break
      fi
    done
    (( _keep_cfg )) && SIM_CFG_TAGS+=( "$_cfg" )
  done

  if (( ${#SIM_CFG_TAGS[@]} == 0 )); then
    err "No cfg_tag subdirectories found under ${_discover_base} that match the CURRENT YAML matrix for requested samples: ${samples[*]}"
    err "YAML used: $(merge_yaml_path)"
    err "Run RecoilJets_Condor_submit.sh ${SIM_DATASET_TOKEN} condorDoAll first to create the output tree."
    exit 20
  fi

  say "Dataset : ${BOLD}${SIM_DATASET_TOKEN}${RST}"
  say "Action  : ${SIM_ACTION}"
  say "Input base  : ${SIM_INPUT_BASE}"
  say "Output base : ${FLAT_OUT_DIR}"
  say "groupSize   : ${SIM_GROUP_SIZE}"
  say "Samples     : ${samples[*]}"
  say "YAML cfg_tags allowed: ${#_YAML_SIM_CFG_TAGS[@]}"
  say "Discovered cfg_tags (filtered by CURRENT YAML + requested samples): ${#SIM_CFG_TAGS[@]}"
  echo

  # Clean TMP_DIR at start of firstRound to avoid stale listfiles/submits
  if [[ "$SIM_ACTION" == "firstRound" ]]; then
    rm -rf "$TMP_DIR" && mkdir -p "$TMP_DIR"
  fi

  final_paths=()
  sim_cleanup_partial_dirs=()

  for cfg_tag in "${SIM_CFG_TAGS[@]}"; do

    COMBO_INPUT_BASE="${SIM_INPUT_BASE}/${cfg_tag}"
    DEST_DIR="${FLAT_OUT_DIR}/${cfg_tag}"
    mkdir -p "$DEST_DIR"

    say "-----------------------------"
    say "cfg_tag : ${BOLD}${cfg_tag}${RST}"
    say "  Input : ${COMBO_INPUT_BASE}"
    say "  Work  : ${DEST_DIR}"
    echo

    for samp in "${samples[@]}"; do
      SIM_SAMPLE="$samp"
      SIM_TAG="${SIM_SAMPLE##*_}"

      if [[ "$SIM_ACTION" == "firstRound" ]]; then
        SIM_INPUT_DIR="${COMBO_INPUT_BASE}/${SIM_SAMPLE}"

        if [[ ! -d "$SIM_INPUT_DIR" ]]; then
          warn "SIM input directory not found (skipping): $SIM_INPUT_DIR"
          continue
        fi

        say "  [scan] sample=${SIM_SAMPLE}"
        say "  [scan] input dir: ${SIM_INPUT_DIR}"
        say "  [scan] discovering candidate ROOT files..."
        mapfile -t SIM_INPUTS_RAW < <(find "$SIM_INPUT_DIR" -maxdepth 1 -type f -name "*.root" -not -name "*_LOCAL_*" -not -name "*_condorTest_*" | sort -V || true)
        say "  [scan] raw ROOT files found: ${#SIM_INPUTS_RAW[@]}"

        if (( ${#SIM_INPUTS_RAW[@]} == 0 )); then
          warn "No raw *.root files found in: $SIM_INPUT_DIR (skipping)"
          continue
        fi

        say "  [fast-filter] using cfg_tag directory + filename match instead of opening ROOT per file"
        say "  [fast-filter] required filename token: ${cfg_tag}"

        mapfile -t SIM_INPUTS < <(find "$SIM_INPUT_DIR" -maxdepth 1 -type f -name "*${cfg_tag}*.root" -not -name "*_LOCAL_*" -not -name "*_condorTest_*" | sort -V || true)

        SIM_CFG_FILENAME_MISMATCH_COUNT=$(( ${#SIM_INPUTS_RAW[@]} - ${#SIM_INPUTS[@]} ))
        say "  [fast-filter] kept=${#SIM_INPUTS[@]} rejected_by_filename=${SIM_CFG_FILENAME_MISMATCH_COUNT}"

        if (( SIM_CFG_FILENAME_MISMATCH_COUNT > 0 )); then
          warn "Some ROOT files in ${SIM_INPUT_DIR} do not contain cfg_tag in the filename and were excluded."
          if (( ${MERGE_SIM_FAST_FILTER_VERBOSE:-0} )); then
            say "  [fast-filter] first rejected filenames:"
            _nrej=0
            for _raw in "${SIM_INPUTS_RAW[@]}"; do
              _matched=0
              for _keep in "${SIM_INPUTS[@]}"; do
                if [[ "$_raw" == "$_keep" ]]; then
                  _matched=1
                  break
                fi
              done
              if (( !_matched )); then
                printf "    %s\n" "$(basename "$_raw")"
                (( _nrej+=1 ))
              fi
              (( _nrej >= 10 )) && break
            done
          fi
        fi

        if (( ${MERGE_SIM_ROOT_STAMP_AUDIT:-0} && ${#SIM_INPUTS[@]} > 0 )); then
          SIM_STAMP_AUDIT_LIMIT="${MERGE_SIM_ROOT_STAMP_AUDIT_LIMIT:-3}"
          [[ "$SIM_STAMP_AUDIT_LIMIT" =~ ^[0-9]+$ ]] || SIM_STAMP_AUDIT_LIMIT=3
          (( SIM_STAMP_AUDIT_LIMIT > 0 )) || SIM_STAMP_AUDIT_LIMIT=3

          say "  [stamp-audit] optional QA enabled: checking up to ${SIM_STAMP_AUDIT_LIMIT} kept files"
          SIM_STAMP_AUDIT_CHECKED=0
          SIM_STAMP_AUDIT_EMPTY=0
          SIM_STAMP_AUDIT_MISMATCH=0

          for _cand in "${SIM_INPUTS[@]}"; do
            (( SIM_STAMP_AUDIT_CHECKED+=1 ))
            _cand_tag="$(extract_cfg_tag_from_root "$_cand" 2>/dev/null || true)"

            if [[ -z "$_cand_tag" ]]; then
              (( SIM_STAMP_AUDIT_EMPTY+=1 ))
            elif [[ "$_cand_tag" != "$cfg_tag" ]]; then
              (( SIM_STAMP_AUDIT_MISMATCH+=1 ))
              warn "[stamp-audit] filename-selected file has mismatched ROOT-stamped tag: $(basename "$_cand")"
              warn "[stamp-audit] expected=${cfg_tag}"
              warn "[stamp-audit] stamped =${_cand_tag}"
            fi

            (( SIM_STAMP_AUDIT_CHECKED >= SIM_STAMP_AUDIT_LIMIT )) && break
          done

          say "  [stamp-audit] checked=${SIM_STAMP_AUDIT_CHECKED} empty=${SIM_STAMP_AUDIT_EMPTY} mismatched=${SIM_STAMP_AUDIT_MISMATCH}"
        fi

        if (( ${#SIM_INPUTS[@]} == 0 )); then
          warn "No cfg-matching *.root files found by filename in: $SIM_INPUT_DIR (skipping)"
          continue
        fi

        SIM_PARTIAL_PREFIX="chunkMerge_${SIM_TAG}_grp"

        if $SIM_FIRSTROUND_LOCAL; then
          say "SIM firstRound (LOCAL): cfg=${cfg_tag} sample=${SIM_SAMPLE} inputs=${#SIM_INPUTS[@]} -> grouped hadd on this node"
        else
          need_cmd condor_submit
          say "SIM firstRound (CONDOR): cfg=${cfg_tag} sample=${SIM_SAMPLE} inputs=${#SIM_INPUTS[@]} -> grouped hadd jobs on Condor"
        fi

        # Clean previous partials for this sample
        find "$DEST_DIR" -maxdepth 1 -type f -name "${SIM_PARTIAL_PREFIX}*.root" -delete || true
        # Clean stale secondRound final so it can't coexist with fresh partials
        rm -f "${FLAT_OUT_DIR}/${FINAL_PREFIX}_${SIM_TAG}_ALL_${cfg_tag}.root" 2>/dev/null || true

        total="${#SIM_INPUTS[@]}"
        grp=0

        if $SIM_FIRSTROUND_LOCAL; then
          need_cmd hadd

          for ((i=0; i<total; i+=SIM_GROUP_SIZE)); do
            (( grp+=1 ))
            grpTag="$(printf "%03d" "$grp")"

            listfile="${TMP_DIR}/sim_${cfg_tag}_${SIM_TAG}_grp${grpTag}.txt"
            : > "$listfile"

            for ((j=i; j<i+SIM_GROUP_SIZE && j<total; j++)); do
              printf "%s\n" "${SIM_INPUTS[$j]}" >> "$listfile"
            done

            out="${DEST_DIR}/${SIM_PARTIAL_PREFIX}${grpTag}.root"
            say "[LOCAL firstRound] cfg=${cfg_tag} sample=${SIM_TAG} grp=${grpTag} inputs=$(wc -l < "$listfile") -> $(basename "$out")"
            hadd -v 3 -f "$out" @"$listfile"
          done

          say "LOCAL firstRound complete for cfg=${cfg_tag} sample=${SIM_SAMPLE}. Partials are under: ${DEST_DIR}"
        else
          emit_hadd_wrapper "$CONDOR_EXEC"

          SUB="${TMP_DIR}/recoil_sim_${cfg_tag}_${SIM_TAG}_firstRound.sub"
          rm -f "$SUB"
          cat > "$SUB" <<EOT
universe   = vanilla
executable = $CONDOR_EXEC
output     = $OUT_DIR/recoil.sim.${cfg_tag}.${SIM_TAG}.\$(Cluster).\$(Process).out
error      = $ERR_DIR/recoil.sim.${cfg_tag}.${SIM_TAG}.\$(Cluster).\$(Process).err
log        = $LOG_DIR/recoil.sim.${cfg_tag}.${SIM_TAG}.\$(Cluster).\$(Process).log
request_memory = $SIM_FIRSTROUND_REQUEST_MEMORY
priority = $MERGE_CONDOR_PRIORITY
getenv = True
should_transfer_files = NO
stream_output = True
stream_error  = True
EOT

          ngroups=$(( (total + SIM_GROUP_SIZE - 1) / SIM_GROUP_SIZE ))
          say "firstRound plan (cfg=${cfg_tag} sample=${SIM_SAMPLE} tag=${SIM_TAG})"
          say "  inputs total     : ${total}"
          say "  groupSize        : ${SIM_GROUP_SIZE}"
          say "  request_memory   : ${SIM_FIRSTROUND_REQUEST_MEMORY}"
          say "  expected groups  : ${ngroups}"
          say "  submit file      : ${SUB}"
          say "  output partials  : ${DEST_DIR}/${SIM_PARTIAL_PREFIX}NNN.root"
          echo

          for ((i=0; i<total; i+=SIM_GROUP_SIZE)); do
            (( grp+=1 ))
            grpTag="$(printf "%03d" "$grp")"

            listfile="${TMP_DIR}/sim_${cfg_tag}_${SIM_TAG}_grp${grpTag}.txt"
            : > "$listfile"

            for ((j=i; j<i+SIM_GROUP_SIZE && j<total; j++)); do
              printf "%s\n" "${SIM_INPUTS[$j]}" >> "$listfile"
            done

            if [[ ! -s "$listfile" ]]; then
              err "firstRound: produced empty listfile (cfg=${cfg_tag} sample=${SIM_SAMPLE} grp=${grpTag}) → ${listfile}"
              exit 26
            fi

            out="${DEST_DIR}/${SIM_PARTIAL_PREFIX}${grpTag}.root"
            n_in=$(wc -l < "$listfile" | awk '{print $1}')
            first_in=$(head -n 1 "$listfile" 2>/dev/null || true)
            last_in=$(tail -n 1 "$listfile" 2>/dev/null || true)

            say "[firstRound queue] cfg=${cfg_tag}  sample=${SIM_TAG}  grp=${grpTag}/${ngroups}  inputs=${n_in}"
            say "  list : ${listfile}"
            say "  out  : ${out}"
            (( n_in > 0 )) && { say "  first: ${first_in}"; say "  last : ${last_in}"; }
            echo

            printf 'arguments = %s %s\nqueue\n\n' "$listfile" "$out" >> "$SUB"
          done

          if [[ ! -s "$SUB" ]]; then
            err "firstRound: submit file is empty/unwritten → ${SUB}"
            exit 27
          fi

          say "Submitting ${BOLD}${grp}${RST} firstRound Condor merge jobs → $(basename "$SUB")"
          condor_submit "$SUB"
          say "FirstRound submitted. Partials will appear under: ${DEST_DIR}"
        fi

      elif [[ "$SIM_ACTION" == "secondRound" ]]; then
        SIM_PARTIAL_PREFIX="chunkMerge_${SIM_TAG}_grp"

        # Merge the firstRound partials into ONE final file (per cfg_tag × sample)
        mapfile -t partials < <(ls -1 "${DEST_DIR}/${SIM_PARTIAL_PREFIX}"*.root 2>/dev/null | sort -V || true)
        if (( ${#partials[@]} == 0 )); then
          warn "No firstRound partials found for cfg=${cfg_tag} sample=${SIM_SAMPLE} in ${DEST_DIR} (expected ${SIM_PARTIAL_PREFIX}*.root). Skipping."
          continue
        fi

        LIST="${TMP_DIR}/sim_${cfg_tag}_${SIM_TAG}_partialList.txt"
        printf "%s\n" "${partials[@]}" > "$LIST"

        # Output uses the directory name as the cfg tag (YAML-immune)
        SIM_FINAL="${FLAT_OUT_DIR}/${FINAL_PREFIX}_${SIM_TAG}_ALL_${cfg_tag}.root"

        # Optional QA cross-check: verify ROOT-stamped tag matches directory name
        _root_tag="$(extract_cfg_tag_from_root "${partials[0]}" 2>/dev/null || true)"
        if [[ -n "$_root_tag" ]]; then
          if [[ "$_root_tag" != "$cfg_tag" ]]; then
            warn "[QA] Directory tag '${cfg_tag}' != ROOT-stamped tag '${_root_tag}' (using directory tag for naming)"
          else
            say "  [QA] ROOT-stamped tag matches directory: ${_root_tag}"
          fi
        fi

        say "SIM secondRound: cfg=${cfg_tag} sample=${SIM_SAMPLE} partials=${#partials[@]} -> ${SIM_FINAL}"

        if $SIM_PREFER_CONDOR; then
          need_cmd condor_submit
          emit_hadd_wrapper "$CONDOR_EXEC"

          SUB="${TMP_DIR}/recoil_sim_${cfg_tag}_${SIM_TAG}_secondRound.sub"
          rm -f "$SUB"
          cat > "$SUB" <<EOT
universe   = vanilla
executable = $CONDOR_EXEC
output     = $OUT_DIR/recoil.sim.${cfg_tag}.${SIM_TAG}.final.\$(Cluster).\$(Process).out
error      = $ERR_DIR/recoil.sim.${cfg_tag}.${SIM_TAG}.final.\$(Cluster).\$(Process).err
log        = $LOG_DIR/recoil.sim.${cfg_tag}.${SIM_TAG}.final.\$(Cluster).\$(Process).log
request_memory = 6GB
priority = $MERGE_CONDOR_PRIORITY
getenv = True
should_transfer_files = NO
stream_output = True
stream_error  = True
arguments = $LIST $SIM_FINAL
queue
EOT
          say "Submitting secondRound final merge on Condor → $(basename "$SUB")"
          condor_submit "$SUB"
        else
          need_cmd hadd
          say "Running secondRound final merge locally (ROOT hadd)…"
          hadd -v 3 -f "$SIM_FINAL" @"$LIST"
          if final_root_is_good "$SIM_FINAL"; then
            say "Created ${SIM_FINAL}"
            final_paths+=( "$SIM_FINAL" )
            sim_cleanup_partial_dirs+=( "$DEST_DIR" )
          else
            err "Local secondRound did not produce a non-empty final ROOT file: ${SIM_FINAL}"
            exit 31
          fi
        fi
      fi

      echo
    done
  done

  if [[ "$SIM_ACTION" == "secondRound" ]] && (( ${#final_paths[@]} > 0 )); then
    say "====================================================================="
    say "SecondRound outputs (final files), grouped by cfg tag:"

    declare -A group_to_paths
    declare -a group_keys=()

    for p in "${final_paths[@]}"; do
      b="$(basename "$p")"
      cfg="${b#*_ALL_}"
      cfg="${cfg%.root}"
      if [[ -z "$cfg" || "$cfg" == "$b" ]]; then
        cfg="UNKNOWN_CFG"
      fi
      if [[ -z "${group_to_paths[$cfg]+x}" ]]; then
        group_keys+=( "$cfg" )
        group_to_paths["$cfg"]="$p"
      else
        group_to_paths["$cfg"]+=$'\n'"$p"
      fi
    done

    IFS=$'\n' sorted_keys=($(printf "%s\n" "${group_keys[@]}" | sort -V))
    unset IFS

    for k in "${sorted_keys[@]}"; do
      say "---------------------------------------------------------------------"
      say "CFG: ${BOLD}${k}${RST}"
      printf "%s\n" "${group_to_paths[$k]}" | sort -V | sed 's/^/  /'
    done

    say "====================================================================="
  fi

  # Clean up firstRound partial directories only after LOCAL secondRound success.
  # Do not clean after Condor final submission; those final jobs may still be pending.
  if [[ "$SIM_ACTION" == "secondRound" ]] && (( ${#sim_cleanup_partial_dirs[@]} > 0 )); then
    say "Cleaning firstRound partial directories from ${FLAT_OUT_DIR} after verified local final merge…"

    mapfile -t _sim_cleanup_dirs_unique < <(printf "%s\n" "${sim_cleanup_partial_dirs[@]}" | sort -u)
    for _partials_dir in "${_sim_cleanup_dirs_unique[@]}"; do
      case "$_partials_dir" in
        "${FLAT_OUT_DIR}"/*)
          if [[ -d "$_partials_dir" ]]; then
            rm -rf -- "$_partials_dir"
            say "  removed: ${_partials_dir}"
          fi
          ;;
        *)
          warn "Refusing to clean suspicious SIM partial dir: ${_partials_dir}"
          ;;
      esac
    done

    say "Done. Only final merged SIM files remain under: ${FLAT_OUT_DIR}"
    cleanup_current_tmp_dir_after_local_final
  fi

  exit 0
fi

# ============================================================
# DATA behavior (DISCOVERY-BASED — cfg_tag aware):
#   condor / addChunks / checkFileOutput
#
# The submit script creates:
#   <RUN_BASE>/<cfg_tag>/<run8>/*.root
# This script discovers cfg_tag subdirectories from the filesystem,
# so it is immune to YAML edits between submission and merging.
#
# Final outputs (flat, one file per cfg_tag):
#   output/<TAG>/RecoilJets_<TAG>_ALL_<cfg_tag>.root
# ============================================================
if [[ "${1:-}" =~ ^(pp|PP|isPP|PP_DATA|pp_data|pp25|PP25|isPPrun25|pprun25|PP_RUN25|pp_run25|auau|AA|isAuAu|AuAu|aa|AA_DATA|auau_data|oo|OO|isOO|OO_DATA|oo_data)$ && "${2:-}" =~ ^(condor|addChunks|checkFileOutput)$ ]]; then
  MODE="$2"
  DATASET_REQ="$1"
else
  MODE="$1"
  DATASET_REQ="$2"
fi
SUBMODE="${3:-}"

[[ "$MODE" =~ ^(condor|addChunks|checkFileOutput)$ ]] || usage
resolve_dataset "$DATASET_REQ"

say "Dataset: ${BOLD}${TAG}${RST}"
say "Run base: ${RUN_BASE}"
say "Output dir: ${DEST_DIR}"

# ---------------------------------------------------------------
# Discover cfg_tag subdirectories under RUN_BASE.
# Skip bare 8-digit run dirs (legacy flat layout).
# Then FILTER those discovered cfg_tags against the CURRENT
# master YAML matrix so we only merge the cfg tags that belong to
# the current analysis configuration.
# ---------------------------------------------------------------
mapfile -t _ALL_DATA_CFG_TAGS < <(
  find "$RUN_BASE" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' |
  grep -v -E '^[0-9]{8}$' | sort -V
)

if (( ${#_ALL_DATA_CFG_TAGS[@]} == 0 )); then
  # Fallback: old flat layout (no cfg_tag subdirs, runs directly under RUN_BASE)
  DATA_CFG_TAGS=( "" )
  say "No cfg_tag subdirectories found; using flat layout (legacy)"
else
  mapfile -t _YAML_DATA_CFG_TAGS < <(build_cfg_tags_from_yaml "$TAG")
  DATA_CFG_TAGS=()
  for _ct in "${_ALL_DATA_CFG_TAGS[@]}"; do
    for _y in "${_YAML_DATA_CFG_TAGS[@]}"; do
      if [[ "$_ct" == "$_y" ]]; then
        DATA_CFG_TAGS+=( "$_ct" )
        break
      fi
    done
  done

  if (( ${#DATA_CFG_TAGS[@]} == 0 )); then
    err "No discovered cfg_tag directories under ${RUN_BASE} match the CURRENT YAML matrix."
    err "YAML used: $(merge_yaml_path)"
    exit 21
  fi

  say "YAML cfg_tags allowed: ${#_YAML_DATA_CFG_TAGS[@]}"
  say "Discovered cfg_tags (filtered by CURRENT YAML): ${#DATA_CFG_TAGS[@]}"
  for _ct in "${DATA_CFG_TAGS[@]}"; do say "  ${_ct}"; done
fi

# ── Optional CFG_FILTER: select a subset of discovered cfg_tags ──
# Recognised values:
#   allButVariantB  – keep every tag EXCEPT those ending in _variantB
#   variantB        – keep ONLY tags ending in _variantB
if [[ -n "${CFG_FILTER:-}" ]]; then
  _filtered=()
  for _ct in "${DATA_CFG_TAGS[@]}"; do
    case "$CFG_FILTER" in
      allButVariantB)
        [[ "$_ct" == *_variantB ]] && continue ;;
      variantB)
        [[ "$_ct" != *_variantB ]] && continue ;;
      *)
        die "Unknown CFG_FILTER='${CFG_FILTER}'. Recognised: allButVariantB, variantB" ;;
    esac
    _filtered+=( "$_ct" )
  done
  DATA_CFG_TAGS=( "${_filtered[@]}" )
  say "CFG_FILTER=${CFG_FILTER} → kept ${#DATA_CFG_TAGS[@]} cfg_tag(s):"
  for _ct in "${DATA_CFG_TAGS[@]}"; do say "  ${_ct}"; done
  if (( ${#DATA_CFG_TAGS[@]} == 0 )); then
    die "CFG_FILTER=${CFG_FILTER} removed all cfg_tags – nothing to do."
  fi
fi
echo

# ---------- Modes ----------
if [[ "$MODE" == "checkFileOutput" ]]; then
  # Build the active-job skiplist once (read-only condor_q)
  build_active_skiplist

  for _cfg in "${DATA_CFG_TAGS[@]}"; do
    if [[ -n "$_cfg" ]]; then
      _run_base="${RUN_BASE}/${_cfg}"
      say "═══════════════════════════════════════════════════"
      say "cfg_tag: ${BOLD}${_cfg}${RST}"
    else
      _run_base="$RUN_BASE"
    fi

    say "Scanning available run directories under ${_run_base}"
    discover_runs "$_run_base"
    printf "Runs with *.root present: %s\n" "${#RUNS[@]}"
    if ((${#RUNS[@]}==0)); then
      warn "No runs found with *.root files."
      continue
    fi

    say "Per-run summary (TOTAL vs BUSY in condor_q vs BUSY present on disk vs ELIGIBLE)"
    printf "  %-10s  %8s  %8s  %12s  %10s\n" "run" "total" "busyInQ" "busyPresent" "eligible"
    printf "  %-10s  %8s  %8s  %12s  %10s\n" "----------" "--------" "--------" "------------" "----------"

    eligibleRuns=0
    for r in "${RUNS[@]}"; do
      srcdir="${_run_base}/${r}"

      all="${TMP_DIR}/check_${_cfg:-flat}_${r}.all.txt"
      find "$srcdir" -maxdepth 1 -type f -name "*.root" -not -name "*_LOCAL_*" | sort -V > "$all"
      total=$(wc -l < "$all" | awk '{print $1}')

      busyInQ=0
      busyPresent=0
      if [[ -s "$SKIP_FILE" ]]; then
        busyInQ=$(awk -v pre="${_run_base}/${r}/" 'index($0,pre)==1{c++} END{print c+0}' "$SKIP_FILE")
        busyPresent=$(grep -F -x -f "$SKIP_FILE" "$all" 2>/dev/null | wc -l | awk '{print $1}')
      fi

      eligible=$(( total - busyPresent ))
      (( eligible > 0 )) && (( eligibleRuns += 1 ))

      printf "  %-10s  %8d  %8d  %12d  %10d\n" "$r" "$total" "$busyInQ" "$busyPresent" "$eligible"
      rm -f "$all" 2>/dev/null || true
    done

    say "Eligible runs (eligible>0): ${eligibleRuns} / ${#RUNS[@]}"
    echo
  done
  say "NOTE: This mode only REPORTS. It does not delete files, submit jobs, or run hadd."
  exit 0
fi

if [[ "$MODE" == "condor" ]]; then
  if (( DRYRUN )); then
    warn "DRYRUN=1 → planning only (NO deletes, NO condor_submit, NO hadd). Reading condor_q is safe."
  else
    need_cmd condor_submit
  fi

  # Clean TMP_DIR at start of stage-1 merge to avoid stale listfiles/submits
  rm -rf "$TMP_DIR" && mkdir -p "$TMP_DIR"

  RUN_BASE_ORIG="$RUN_BASE"
  DEST_DIR_ORIG="$DEST_DIR"
  PERRUN_DIR_ORIG="$PERRUN_DIR"

  # Build the skiplist once up-front (covers all cfg_tags; paths are absolute)
  build_active_skiplist

  for _cfg in "${DATA_CFG_TAGS[@]}"; do
    if [[ -n "$_cfg" ]]; then
      RUN_BASE="${RUN_BASE_ORIG}/${_cfg}"
      DEST_DIR="${DEST_DIR_ORIG}/${_cfg}"
      PERRUN_CFG_DIR="${PERRUN_DIR_ORIG}/${_cfg}"
      say "═══════════════════════════════════════════════════"
      say "cfg_tag: ${BOLD}${_cfg}${RST}"
      say "  RUN_BASE: ${RUN_BASE}"
      say "  PERRUN_CFG_DIR: ${PERRUN_CFG_DIR}"
    else
      RUN_BASE="$RUN_BASE_ORIG"
      DEST_DIR="$DEST_DIR_ORIG"
      PERRUN_CFG_DIR="$PERRUN_DIR_ORIG"
    fi

    mkdir -p "$PERRUN_CFG_DIR"

    say "Preparing per-run Condor merges (one job per run)"

    SEGMENT_MODE=0
    SEGMENTS=()

    if [[ "$SUBMODE" =~ ^(fromSplitRunList|fromSplit|rounds|round|segments|seg)$ ]]; then
      SEGMENT_MODE=1
      for tok in "${@:4}"; do
        [[ "$tok" =~ ^[0-9]+$ ]] && SEGMENTS+=( "$tok" )
      done
      if (( ${#SEGMENTS[@]} == 0 )); then
        err "Segment-restricted mode '${SUBMODE}' requires one or more segment indices."
        err "Example: $0 condor ${TAG} fromSplitRunList round 1 round 3 currentNode"
        exit 6
      fi

      load_runs_from_segments "${SEGMENTS[@]}"

      if ((${#RUNS[@]}==0)); then
        err "No runs loaded from segments (${SEGMENTS[*]}). Check: ${ROUND_DIR}/goldenRuns_${TAG}_segment*.txt"
        exit 5
      fi
    else
      discover_runs "$RUN_BASE"
      if ((${#RUNS[@]}==0)); then
        warn "No runs found with *.root files under ${RUN_BASE}; skipping cfg_tag '${_cfg:-flat}'"
        continue
      fi

      # Apply legacy submodes
      case "$SUBMODE" in
        test)        RUNS=( "${RUNS[0]}" ) ;;
        firstHalf)   RUNS=( "${RUNS[@]:0:$(( (${#RUNS[@]}+1)/2 ))}" ) ;;
        "" )         ;;
        * )          err "Unknown submode '$SUBMODE' (allowed: test, firstHalf, fromSplitRunList/rounds/segments)"; exit 6 ;;
      esac
    fi

    if (( DRYRUN )); then
      say "DRYRUN plan: which per-run partial merges WOULD be submitted"
      planned=0
      for r in "${RUNS[@]}"; do
        listfile="$(make_run_list "$r" || true)"
        if [[ -z "${listfile:-}" || ! -s "$listfile" ]]; then
          warn "Run $r: no eligible ROOT files after skipping active jobs (or list build failed)"
          continue
        fi
        nfiles=$(wc -l < "$listfile" | awk '{print $1}')
        out="${PERRUN_CFG_DIR}/${PARTIAL_PREFIX}_${r}.root"
        say "  would merge run ${r}: inputs=${nfiles} -> ${out}"
        ((planned+=1))
      done
      say "DRYRUN summary (cfg=${_cfg:-flat}): planned runs=${planned}. (No files deleted; no jobs submitted.)"
      continue
    fi

    # Normal mode:
    if (( SEGMENT_MODE )); then
      say "Segment-restricted stage-1: ${BOLD}NOT deleting${RST} existing partials in ${PERRUN_CFG_DIR}"
      say "  (This avoids clobbering partials produced from another schedd/node.)"
      say "  Any run we submit will overwrite its own partial: ${PARTIAL_PREFIX}_<run8>.root"
      rm -f "${DEST_DIR_ORIG}/${FINAL_PREFIX}_${TAG}_ALL_${_cfg}"*.root 2>/dev/null || true
    else
      say "Cleaning old perRun partials in ${PERRUN_CFG_DIR}"
      rm -rf "${PERRUN_CFG_DIR}" 2>/dev/null || true
      mkdir -p "${PERRUN_CFG_DIR}"
      say "Cleaning stale finals in ${DEST_DIR_ORIG}"
      rm -f "${DEST_DIR_ORIG}/${FINAL_PREFIX}_${TAG}_ALL_${_cfg}"*.root 2>/dev/null || true
    fi

    emit_hadd_wrapper "$CONDOR_EXEC"

    SUB="${TMP_DIR}/recoil_partials_${TAG}_${_cfg:-flat}.sub"
    rm -f "$SUB"
    cat > "$SUB" <<EOT
universe   = vanilla
executable = $CONDOR_EXEC
output     = $OUT_DIR/recoil.\$(Cluster).\$(Process).out
error      = $ERR_DIR/recoil.\$(Cluster).\$(Process).err
log        = $LOG_DIR/recoil.\$(Cluster).\$(Process).log
request_memory = 2GB
priority = $MERGE_CONDOR_PRIORITY
getenv = True
should_transfer_files = NO
stream_output = True
stream_error  = True
EOT

    queued=0
    skipped=0
    for r in "${RUNS[@]}"; do
      listfile="$(make_run_list "$r" || true)"
      if [[ -z "${listfile:-}" || ! -s "$listfile" ]]; then
        warn "Run $r: no eligible files after skipping active jobs; skipping merge submit"
        if (( SEGMENT_MODE )); then
          [[ -f "${DEST_DIR}/${PARTIAL_PREFIX}_${r}.root" ]] && \
            warn "  (segment-mode) existing partial left as-is: ${DEST_DIR}/${PARTIAL_PREFIX}_${r}.root"
        fi
        ((skipped+=1))
        continue
      fi

      out="${PERRUN_CFG_DIR}/${PARTIAL_PREFIX}_${r}.root"
      nfiles=$(wc -l < "$listfile" | awk '{print $1}')
      (( SKIP_TRACE )) && say "Queueing run ${r}: inputs=${nfiles} -> $(basename "$out")"

      printf 'arguments = %s %s\nqueue\n\n' "$listfile" "$out" >> "$SUB"
      ((queued+=1))
    done

    (( SKIP_TRACE )) && say "Stage-1 queue summary (cfg=${_cfg:-flat}): queued=${queued}, skipped=${skipped}"

    if (( queued == 0 )); then
      warn "No Condor jobs to submit for cfg_tag '${_cfg:-flat}' (no non-empty eligible run lists)."
      continue
    fi

    say "Submitting ${BOLD}${queued}${RST} Condor merge jobs → $(basename "$SUB")"
    condor_submit "$SUB"
    say "Stage-1 submitted. Partial outputs will appear under: ${PERRUN_CFG_DIR}"
    echo
  done

  if (( DRYRUN )); then
    say "DRYRUN complete."
  fi
  exit 0
fi

if [[ "$MODE" == "addChunks" ]]; then
  # Final merge, prefer local unless "condor" explicitly provided
  PREFER_CONDOR=false
  [[ "${SUBMODE:-}" == "condor" ]] && PREFER_CONDOR=true

  # sliceRuns: intermediate step that groups per-run partials into batches
  # before the final merge, reducing single-job input counts dramatically.
  SLICE_RUNS=false
  [[ "${4:-}" == "sliceRuns" ]] && SLICE_RUNS=true
  SLICE_BATCH_SIZE=200   # number of per-run partials per slice hadd job

  DEST_DIR_ORIG="$DEST_DIR"
  _PERRUN_BASE="${OUT_BASE}/${TAG}/perRun"

  # Discover cfg_tag subdirs under perRun (where stage-1 partials were written)
  mapfile -t _ALL_MERGE_CFG_TAGS < <(
    find "$_PERRUN_BASE" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' 2>/dev/null | sort -V
  )

  if (( ${#_ALL_MERGE_CFG_TAGS[@]} == 0 )); then
    # Fallback: flat layout (partials directly in perRun base)
    MERGE_CFG_TAGS=( "" )
  else
    # Filter discovered perRun cfg_tags against the CURRENT YAML matrix
    mapfile -t _YAML_MERGE_CFG_TAGS < <(build_cfg_tags_from_yaml "$TAG")
    MERGE_CFG_TAGS=()
    for _ct in "${_ALL_MERGE_CFG_TAGS[@]}"; do
      for _y in "${_YAML_MERGE_CFG_TAGS[@]}"; do
        [[ "$_ct" == "$_y" ]] && { MERGE_CFG_TAGS+=( "$_ct" ); break; }
      done
    done
    say "YAML cfg_tags allowed: ${#_YAML_MERGE_CFG_TAGS[@]}"
    say "Discovered cfg_tags (filtered by CURRENT YAML): ${#MERGE_CFG_TAGS[@]}"
    for _ct in "${MERGE_CFG_TAGS[@]}"; do say "  ${_ct}"; done
    if (( ${#MERGE_CFG_TAGS[@]} == 0 )); then
      die "No perRun cfg_tag directories under ${_PERRUN_BASE} match the CURRENT YAML matrix."
    fi
  fi

  # ── Optional CFG_FILTER (same env-var as stage-1) ──
  if [[ -n "${CFG_FILTER:-}" ]]; then
    _filtered=()
    for _ct in "${MERGE_CFG_TAGS[@]}"; do
      case "$CFG_FILTER" in
        allButVariantB)
          [[ "$_ct" == *_variantB ]] && continue ;;
        variantB)
          [[ "$_ct" != *_variantB ]] && continue ;;
        *)
          die "Unknown CFG_FILTER='${CFG_FILTER}'. Recognised: allButVariantB, variantB" ;;
      esac
      _filtered+=( "$_ct" )
    done
    MERGE_CFG_TAGS=( "${_filtered[@]}" )
    say "CFG_FILTER=${CFG_FILTER} → kept ${#MERGE_CFG_TAGS[@]} merge cfg_tag(s):"
    for _ct in "${MERGE_CFG_TAGS[@]}"; do say "  ${_ct}"; done
    if (( ${#MERGE_CFG_TAGS[@]} == 0 )); then
      die "CFG_FILTER=${CFG_FILTER} removed all merge cfg_tags – nothing to do."
    fi
  fi

  data_cleanup_slice_dirs=()
  data_local_final_success=0

  for _cfg in "${MERGE_CFG_TAGS[@]}"; do
    if [[ -n "$_cfg" ]]; then
      _partial_dir="${_PERRUN_BASE}/${_cfg}"
      say "═══════════════════════════════════════════════════"
      say "cfg_tag: ${BOLD}${_cfg}${RST}"
    else
      _partial_dir="$_PERRUN_BASE"
    fi

    scale_scaled_trig_qa_partials "$_partial_dir"

    # ── sliceRuns: intermediate merge of per-run partials into batches ──
    if $SLICE_RUNS; then
      mapfile -t _run_partials < <(ls -1 "${_partial_dir}/${PARTIAL_PREFIX}_"*.root 2>/dev/null | sort -V || true)
      if (( ${#_run_partials[@]} == 0 )); then
        warn "No per-run partials in ${_partial_dir} for sliceRuns. Skipping."
        continue
      fi

      _total=${#_run_partials[@]}
      _ngroups=$(( (_total + SLICE_BATCH_SIZE - 1) / SLICE_BATCH_SIZE ))

      say "sliceRuns: ${_total} per-run partials → ${_ngroups} slice groups (batch size ${SLICE_BATCH_SIZE})"

      # Clean old sliceMerge files before writing new ones
      rm -f "${_partial_dir}/sliceMerge_grp"*.root 2>/dev/null || true

      need_cmd condor_submit
      emit_hadd_wrapper "$CONDOR_EXEC"

      SUB="${TMP_DIR}/recoil_slice_${TAG}_${_cfg:-flat}.sub"
      rm -f "$SUB"
      cat > "$SUB" <<EOT
universe   = vanilla
executable = $CONDOR_EXEC
output     = $OUT_DIR/recoil.slice.\$(Cluster).\$(Process).out
error      = $ERR_DIR/recoil.slice.\$(Cluster).\$(Process).err
log        = $LOG_DIR/recoil.slice.\$(Cluster).\$(Process).log
request_memory = 4GB
priority = $MERGE_CONDOR_PRIORITY
getenv = True
should_transfer_files = NO
stream_output = True
stream_error  = True
EOT

      _offset=0
      for (( _g=1; _g<=_ngroups; _g++ )); do
        _grpTag="$(printf "%03d" "$_g")"
        _end=$(( _offset + SLICE_BATCH_SIZE ))
        (( _end > _total )) && _end=$_total
        _sz=$(( _end - _offset ))

        _list="${TMP_DIR}/sliceList_${TAG}_${_cfg:-flat}_grp${_grpTag}.txt"
        printf "%s\n" "${_run_partials[@]:${_offset}:${_sz}}" > "$_list"
        _out="${_partial_dir}/sliceMerge_grp${_grpTag}.root"

        say "  slice grp${_grpTag}: ${_sz} partials → $(basename "$_out")"
        printf 'arguments = %s %s\nqueue\n\n' "$_list" "$_out" >> "$SUB"

        (( _offset = _end ))
      done

      say "Submitting ${_ngroups} slice merge jobs → $(basename "$SUB")"
      condor_submit "$SUB"
      echo
      continue
    fi

    # Collect partials: prefer sliceMerge files if present, else per-run partials
    _using_slice_partials=0
    mapfile -t partials < <(ls -1 "${_partial_dir}/sliceMerge_grp"*.root 2>/dev/null | sort -V || true)
    if (( ${#partials[@]} > 0 )); then
      _using_slice_partials=1
      say "  Using ${#partials[@]} sliceMerge partials (from prior sliceRuns step)"
    else
      mapfile -t partials < <(ls -1 "${_partial_dir}/${PARTIAL_PREFIX}_"*.root 2>/dev/null | sort -V || true)
    fi
    if (( ${#partials[@]} == 0 )); then
      warn "No partials found in ${_partial_dir} (expected sliceMerge_grp*.root or ${PARTIAL_PREFIX}_*.root). Skipping."
      continue
    fi

    LIST="${TMP_DIR}/partialList_${TAG}_${_cfg:-flat}.txt"
    printf "%s\n" "${partials[@]}" > "$LIST"

    # Output filename uses the directory name as the cfg tag (YAML-immune)
    if [[ -n "$_cfg" ]]; then
      FINAL="${DEST_DIR_ORIG}/${FINAL_PREFIX}_${TAG}_ALL_${_cfg}.root"

      # Optional QA cross-check: verify ROOT-stamped tag matches directory name
      _root_tag="$(extract_cfg_tag_from_root "${partials[0]}" 2>/dev/null || true)"
      if [[ -n "$_root_tag" ]]; then
        if [[ "$_root_tag" != "$_cfg" ]]; then
          warn "[QA] Directory tag '${_cfg}' != ROOT-stamped tag '${_root_tag}' (using directory tag for naming)"
        else
          say "  [QA] ROOT-stamped tag matches directory: ${_root_tag}"
        fi
      fi
    else
      FINAL="${DEST_DIR_ORIG}/${FINAL_PREFIX}_${TAG}_ALL.root"
    fi

    say "Final merge target: ${FINAL}"
    say "Inputs: ${#partials[@]} partials"

    if $PREFER_CONDOR; then
      need_cmd condor_submit
      emit_hadd_wrapper "$CONDOR_EXEC"

      SUB="${TMP_DIR}/recoil_final_${TAG}_${_cfg:-flat}.sub"
      rm -f "$SUB"
      cat > "$SUB" <<EOT
universe   = vanilla
executable = $CONDOR_EXEC
output     = $OUT_DIR/recoil.final.\$(Cluster).\$(Process).out
error      = $ERR_DIR/recoil.final.\$(Cluster).\$(Process).err
log        = $LOG_DIR/recoil.final.\$(Cluster).\$(Process).log
request_memory = 3GB
priority = $MERGE_CONDOR_PRIORITY
getenv = True
should_transfer_files = NO
stream_output = True
stream_error  = True
arguments = $LIST $FINAL
queue
EOT
      say "Submitting final merge on Condor → $(basename "$SUB")"
      condor_submit "$SUB"
    else
      say "Running final merge locally (ROOT hadd)…"
      hadd -v 3 -f "$FINAL" @"$LIST"
      if final_root_is_good "$FINAL"; then
        say "Created ${FINAL}"
        data_local_final_success=1
        if (( _using_slice_partials )); then
          data_cleanup_slice_dirs+=( "$_partial_dir" )
        fi
      else
        err "Local addChunks did not produce a non-empty final ROOT file: ${FINAL}"
        exit 32
      fi
    fi
    echo
  done

  if ! $PREFER_CONDOR && (( data_local_final_success )); then
    if (( ${#data_cleanup_slice_dirs[@]} > 0 )); then
      say "Cleaning DATA sliceMerge intermediates after verified local final merge…"

      mapfile -t _data_cleanup_slice_dirs_unique < <(printf "%s\n" "${data_cleanup_slice_dirs[@]}" | sort -u)
      for _slice_dir in "${_data_cleanup_slice_dirs_unique[@]}"; do
        case "$_slice_dir" in
          "${_PERRUN_BASE}"|"${_PERRUN_BASE}"/*)
            if compgen -G "${_slice_dir}/sliceMerge_grp*.root" > /dev/null; then
              find "$_slice_dir" -maxdepth 1 -type f -name "sliceMerge_grp*.root" -delete
              say "  removed sliceMerge_grp*.root from: ${_slice_dir}"
            fi
            ;;
          *)
            warn "Refusing to clean suspicious DATA slice dir: ${_slice_dir}"
            ;;
        esac
      done

      say "Done. DATA per-run chunkMerge_run_<run8>.root files were kept."
    fi

    cleanup_current_tmp_dir_after_local_final
  fi

  exit 0
fi

# Fallback (should not reach here)
usage
