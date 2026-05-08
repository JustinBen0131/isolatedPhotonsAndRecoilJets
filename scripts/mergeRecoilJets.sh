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
#   ./mergeRecoilJets.sh isSimInclusive firstRound
#   ./mergeRecoilJets.sh isSimEmbedded firstRound
#
# STEP 2 — secondRound (merge chunk-partials into one final file per cfg_tag)
#   ./mergeRecoilJets.sh isSim secondRound
#   ./mergeRecoilJets.sh isSim secondRound condor
#   ./mergeRecoilJets.sh isSim secondRound condor SAMPLE=run28_photonjet10
#
# STEP 3 — finalStitch (Condor weighted stitch across SIM samples per cfg_tag)
#   ./mergeRecoilJets.sh isSim finalStitch condor
#   ./mergeRecoilJets.sh isSimEmbedded finalStitch condor
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
# SIM finalStitch combined files:
#   output/<simTag>/<cfg_tag>/<comboDir>/RecoilJets_*_MERGED.root
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
#   • RJ_STAGE_EMAIL_MODE=per_cfg (default) sends one READY/CHECK email per
#     Condor merge stage when notify emails are configured.
#   • RJ_STAGE_EMAIL_MODE=none suppresses those emails while still wrapping the
#     merge in a tracking DAG when RJ_STAGE_EMAIL_STRICT=1 is set. This is used
#     by the top-level production DAG to avoid inbox spam while preserving
#     failure propagation.
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

# scaledTrigQA correction diagnostics. These are intentionally verbose by
# default because this step opens many ROOT files one-by-one before sliceRuns.
SCALED_TRIG_VERBOSE="${SCALED_TRIG_VERBOSE:-1}"
SCALED_TRIG_HEARTBEAT_SECONDS="${SCALED_TRIG_HEARTBEAT_SECONDS:-30}"
SCALED_TRIG_ADDCHUNKS_APPLY="${SCALED_TRIG_ADDCHUNKS_APPLY:-0}"

# ---------- Pretty printing ----------
BOLD=$'\e[1m'; RED=$'\e[31m'; YEL=$'\e[33m'; GRN=$'\e[32m'; BLU=$'\e[34m'; RST=$'\e[0m'
say()  { printf "${BLU}➜${RST} %s\n" "$*"; }
warn() { printf "${YEL}⚠ %s${RST}\n" "$*" >&2; }
err()  { printf "${RED}✘ %s${RST}\n" "$*" >&2; }

memory_spec_to_mb() {
  local spec="${1:-}"
  spec="${spec//[[:space:]]/}"
  local upper="${spec^^}"
  if [[ "$upper" =~ ^([0-9]+)GB$ ]]; then
    printf '%d\n' "$(( ${BASH_REMATCH[1]} * 1024 ))"
  elif [[ "$upper" =~ ^([0-9]+)MB$ ]]; then
    printf '%d\n' "${BASH_REMATCH[1]}"
  elif [[ "$upper" =~ ^[0-9]+$ ]]; then
    printf '%d\n' "$upper"
  else
    return 1
  fi
}

condor_auto_memory_retry_block() {
  local base_spec="${1:?base memory required}"
  local base_mb
  base_mb="$(memory_spec_to_mb "$base_spec")" || {
    printf 'request_memory = %s\n' "$base_spec"
    return 0
  }
  local cap_mb="${RJ_AUTO_MEMORY_RETRY_CAP_MB:-8000}"
  local retries="${RJ_AUTO_MEMORY_RETRY_MAX_RELEASES:-2}"
  [[ "$cap_mb" =~ ^[0-9]+$ ]] || cap_mb=8000
  [[ "$retries" =~ ^[0-9]+$ ]] || retries=2
  if [[ "${RJ_AUTO_MEMORY_RETRY:-1}" == "0" || "$base_mb" -ge "$cap_mb" || "$retries" -le 0 ]]; then
    printf 'request_memory = %s\n' "$base_mb"
    return 0
  fi
  cat <<EOT
+RJBaseRequestMemoryMb = ${base_mb}
+RJMemoryRetryCapMb = ${cap_mb}
+RJMemoryRetryMaxReleases = ${retries}
request_memory = ifThenElse(MemoryUsage =!= undefined, ifThenElse(int(MemoryUsage * 1.35 + 512) > ${cap_mb}, ${cap_mb}, ifThenElse(int(MemoryUsage * 1.35 + 512) > ${base_mb}, int(MemoryUsage * 1.35 + 512), ${base_mb})), ${base_mb})
periodic_release = (JobStatus == 5) && (NumJobStarts <= ${retries}) && (RequestMemory < ${cap_mb}) && ((HoldReasonCode == 34) || ((HoldReason =!= undefined) && regexp("memory|Memory|cgroup|request_memory|request memory", HoldReason)))
EOT
}

# ---------- Fixed dataset roots ----------
RUN_BASE_PP="/sphenix/tg/tg01/bulk/jbennett/thesisAna/pp"
RUN_BASE_PP25="/sphenix/tg/tg01/bulk/jbennett/thesisAna/pp25"
RUN_BASE_AA="/sphenix/tg/tg01/bulk/jbennett/thesisAna/auau"
RUN_BASE_OO="/sphenix/tg/tg01/bulk/jbennett/thesisAna/oo"

# ---------- Output base (required by you) ----------
OUT_BASE="${MERGE_OUT_BASE_OVERRIDE:-/sphenix/u/patsfan753/scratch/thesisAnalysis/output}"

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
  std::cout << "[scaledTrigQA] run=${run_num} file=${rootfile}\\n";
  std::cout << "[scaledTrigQA] factors baseline=${baselineScale} pho10=${pho10Scale} pho12=${pho12Scale}\\n";
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
      if (h)
      {
        ++foundTargets;
        std::cout << "[scaledTrigQA] target found: " << objPaths[i] << " factor=" << factors[i] << "\\n";
      }
      else
      {
        std::cout << "[scaledTrigQA] target missing: " << objPaths[i] << "\\n";
      }
    }

    if (foundTargets == 0)
    {
      std::cout << "[scaledTrigQA] no target histograms exist in this file; leaving it unchanged\\n";
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
        std::cout << "[scaledTrigQA] scaled: " << objPaths[i] << " factor=" << factors[i] << "\\n";
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
        std::cout << "[scaledTrigQA] targets existed but no writable directory/positive factor combination was available\\n";
        std::cout << "__SCALEDTRIGQA_STATUS__ skip_no_writable_targets\\n";
      }
    }
  }
  else
  {
    std::cout << "[scaledTrigQA] marker already present; leaving file unchanged\\n";
    std::cout << "__SCALEDTRIGQA_STATUS__ already_applied\\n";
  }
  f.Close();
}
EOF

  local root_log="${TMP_DIR}/scaledTrigQA_scale_run_${run_num}.log"
  local start_ts elapsed root_rc root_pid watcher_pid
  start_ts="$(date +%s)"
  if [[ "$SCALED_TRIG_VERBOSE" != "0" ]]; then
    say "scaledTrigQA: run ${run} ROOT check/update start: ${rootfile}"
    say "scaledTrigQA: run ${run} scale factors: baseline=${baselineScale}, pho10=${pho10Scale}, pho12=${pho12Scale}"
  fi

  local old_int_trap
  old_int_trap="$(trap -p INT || true)"
  SCALED_TRIG_ROOT_PID=""
  trap 'warn "scaledTrigQA: interrupt received; terminating ROOT child"; if [[ -n "${SCALED_TRIG_ROOT_PID:-}" ]]; then kill "$SCALED_TRIG_ROOT_PID" 2>/dev/null || true; fi; exit 130' INT

  set +e
  root -l -b -q "$macro" > "$root_log" 2>&1 &
  root_pid=$!
  SCALED_TRIG_ROOT_PID="$root_pid"
  (
    sleep "$SCALED_TRIG_HEARTBEAT_SECONDS"
    while kill -0 "$root_pid" 2>/dev/null; do
      elapsed=$(( $(date +%s) - start_ts ))
      warn "scaledTrigQA: run ${run} still inside ROOT after ${elapsed}s: ${rootfile}"
      sleep "$SCALED_TRIG_HEARTBEAT_SECONDS"
    done
  ) &
  watcher_pid=$!
  wait "$root_pid"
  root_rc=$?
  kill "$watcher_pid" 2>/dev/null || true
  wait "$watcher_pid" 2>/dev/null || true
  set -e
  SCALED_TRIG_ROOT_PID=""
  if [[ -n "$old_int_trap" ]]; then
    eval "$old_int_trap"
  else
    trap - INT
  fi

  elapsed=$(( $(date +%s) - start_ts ))
  rm -f "$macro"

  if (( root_rc != 0 )); then
    if (( root_rc == 130 || root_rc == 143 )); then
      warn "scaledTrigQA: interrupted while scaling ${rootfile}; stopping"
      return "$root_rc"
    fi
    warn "scaledTrigQA: ROOT scaling failed for ${rootfile} after ${elapsed}s"
    if [[ -s "$root_log" ]]; then
      warn "scaledTrigQA: ROOT log for failed run ${run}:"
      sed 's/^/[scaledTrigQA ROOT] /' "$root_log" >&2 || true
    fi
    SCALED_TRIG_LAST_STATUS="failed"
    rm -f "$root_log"
    return 0
  fi

  local root_output
  root_output="$(cat "$root_log" 2>/dev/null || true)"
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

  if [[ "$SCALED_TRIG_VERBOSE" != "0" ]]; then
    case "$SCALED_TRIG_LAST_STATUS" in
      scaled|skip_no_writable_targets|already_applied|failed|unknown)
        sed 's/^/[scaledTrigQA ROOT] /' "$root_log" >&2 || true
        ;;
      skip_no_targets)
        grep -E 'target missing|no target histograms|__SCALEDTRIGQA_STATUS__' "$root_log" | sed 's/^/[scaledTrigQA ROOT] /' >&2 || true
        ;;
    esac
    say "scaledTrigQA: run ${run} status=${SCALED_TRIG_LAST_STATUS} elapsed=${elapsed}s"
  fi

  rm -f "$root_log"
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

merge_notify_emails_csv() {
  if [[ -n "${RJ_NOTIFY_EMAILS:-}" ]]; then
    printf "%s\n" "$RJ_NOTIFY_EMAILS"
    return 0
  fi
  local yaml
  yaml="$(merge_yaml_path)"
  [[ -f "$yaml" ]] || return 0
  local -a emails=()
  mapfile -t emails < <(yaml_get_inline_list "$yaml" "notify_emails")
  local out=""
  local e
  for e in "${emails[@]}"; do
    e="$(trim_ws "$e")"
    [[ -z "$e" ]] && continue
    out="${out:+${out},}${e}"
  done
  printf "%s\n" "$out"
}

merge_stage_email_mode() {
  local mode="${RJ_STAGE_EMAIL_MODE:-per_cfg}"
  case "$mode" in
    none|quiet|off|0|false|FALSE|no|NO)
      printf '%s\n' "none"
      ;;
    per_cfg|per-cfg|cfg|default|"")
      printf '%s\n' "per_cfg"
      ;;
    *)
      warn "Unknown RJ_STAGE_EMAIL_MODE='${mode}', using per_cfg"
      printf '%s\n' "per_cfg"
      ;;
  esac
}

merge_stage_strict_enabled() {
  case "${RJ_STAGE_EMAIL_STRICT:-0}" in
    1|true|TRUE|yes|YES|on|ON) return 0 ;;
  esac
  return 1
}

merge_stage_tracking_description() {
  local mode
  mode="$(merge_stage_email_mode)"
  if [[ "$mode" == "none" ]]; then
    if merge_stage_strict_enabled; then
      printf '%s\n' "quiet strict DAG validation (no per-cfg email)"
    else
      printf '%s\n' "quiet DAG tracking (no per-cfg email)"
    fi
  else
    printf '%s\n' "one-shot per-cfg READY/CHECK email"
  fi
}

submit_condor_stage_with_ready_email() {
  local sub="$1"
  local stage_key="$2"
  local legacy_subject="$3"
  local body="$4"
  local validation_file="${5:-}"
  : "$legacy_subject"

  local emails
  emails="$(merge_notify_emails_csv)"
  local email_mode
  email_mode="$(merge_stage_email_mode)"
  local strict="0"
  merge_stage_strict_enabled && strict="1"
  if [[ "$email_mode" == "none" ]]; then
    emails="__none__"
  fi

  if [[ -z "$emails" && "$strict" != "1" ]]; then
    condor_submit "$sub"
    return $?
  fi
  [[ -n "$emails" ]] || emails="__none__"

  if ! command -v condor_submit_dag >/dev/null 2>&1; then
    warn "condor_submit_dag not found; submitting merge stage directly without DAG-stage tracking"
    condor_submit "$sub"
    return $?
  fi

  local notify_exec="${TMP_DIR}/notify_${stage_key}_$$.sh"
  local notify_sub="${TMP_DIR}/notify_${stage_key}_$$.sub"
  local notify_body="${TMP_DIR}/notify_${stage_key}_$$.txt"
  local dag="${TMP_DIR}/notify_${stage_key}_$$.dag"
  local validation_arg="__none__"
  [[ -n "$validation_file" ]] && validation_arg="$validation_file"

  printf "%s\n" "$body" > "$notify_body"
  cat > "$notify_exec" <<'EOS'
#!/usr/bin/env bash
set -euo pipefail
emails="${1:?emails required}"
stage_key="${2:?stage key required}"
body_file="${3:?body file required}"
dag_file="${4:?dag file required}"
submit_file="${5:?submit file required}"
out_dir="${6:?out dir required}"
err_dir="${7:?err dir required}"
log_dir="${8:?log dir required}"
validation_file="${9:-__none__}"
email_mode="${10:-per_cfg}"
strict="${11:-0}"

status="READY"
status_note="Condor stage DAG reached its final notification node."
dagman_out="${dag_file}.dagman.out"
nodes_log="${dag_file}.nodes.log"
rescue_count=0
expected_count=0
present_count=0
missing_count=0
empty_count=0
first_missing=""
first_empty=""
if compgen -G "${dag_file}.rescue*" >/dev/null 2>&1; then
  rescue_count=$(compgen -G "${dag_file}.rescue*" | wc -l | awk '{print $1}')
  status="FAILED"
  status_note="DAGMan rescue file(s) were produced; inspect the DAG logs before continuing."
elif [[ -s "$dagman_out" ]] && grep -Eiq 'failed with|DAG abort|aborted|Job was held|held job|POST Script failed|Node Status:[[:space:]]*STATUS_ERROR|Error:[[:space:]].*failed|return value [1-9][0-9]*' "$dagman_out"; then
  status="CHECK"
  status_note="DAGMan log contains error/hold-like text; inspect logs before treating outputs as ready."
fi

if [[ "$validation_file" != "__none__" && -s "$validation_file" ]]; then
  while IFS= read -r expected_path; do
    [[ -n "$expected_path" ]] || continue
    [[ "$expected_path" == \#* ]] && continue
    expected_count=$((expected_count + 1))
    if [[ ! -e "$expected_path" ]]; then
      missing_count=$((missing_count + 1))
      [[ -z "$first_missing" ]] && first_missing="$expected_path"
    elif [[ ! -s "$expected_path" ]]; then
      empty_count=$((empty_count + 1))
      [[ -z "$first_empty" ]] && first_empty="$expected_path"
    else
      present_count=$((present_count + 1))
    fi
  done < "$validation_file"
  if (( missing_count > 0 || empty_count > 0 )); then
    status="CHECK"
    status_note="Stage DAG completed, but expected output validation failed; inspect merge outputs before continuing."
  fi
fi

subject="[RecoilJets][${stage_key}][${status}]"
message_file="$(mktemp "${TMPDIR:-/tmp}/recoiljets_stage_notify.XXXXXX")"
{
  echo "RECOILJETS_STAGE_EMAIL_V1"
  echo "status=${status}"
  echo "status_note=${status_note}"
  echo "stage=${stage_key}"
  echo "stage_type=merge_dag"
  echo "email_mode=${email_mode}"
  echo "strict_validation=${strict}"
  echo "submit_file=${submit_file}"
  echo "dag_file=${dag_file}"
  echo "dagman_out=${dagman_out}"
  echo "nodes_log=${nodes_log}"
  echo "rescue_file_count=${rescue_count}"
  echo "job_output_dir=${out_dir}"
  echo "job_error_dir=${err_dir}"
  echo "job_log_dir=${log_dir}"
  echo "validation_file=${validation_file}"
  echo "expected_outputs=${expected_count}"
  echo "present_outputs=${present_count}"
  echo "missing_outputs=${missing_count}"
  echo "empty_outputs=${empty_count}"
  echo "first_missing_output=${first_missing}"
  echo "first_empty_output=${first_empty}"
  echo "next_action=If status is READY, continue with the next merge/pull/plot step named below. If status is CHECK or FAILED, inspect dagman_out, nodes_log, and the job .err files first."
  echo
  echo "message:"
  sed 's/^/  /' "$body_file"
} > "$message_file"

if [[ "$email_mode" == "none" || "$emails" == "__none__" ]]; then
  echo "[notify] ${stage_key}: email suppressed by RJ_STAGE_EMAIL_MODE=none; status=${status}" >&2
  cat "$message_file" >&2
elif command -v mail >/dev/null 2>&1; then
  mail -s "$subject" "$emails" < "$message_file" || echo "[notify] mail command failed for ${emails}; status=${status}" >&2
elif command -v mailx >/dev/null 2>&1; then
  mailx -s "$subject" "$emails" < "$message_file" || echo "[notify] mailx command failed for ${emails}; status=${status}" >&2
else
  echo "[notify] mail/mailx not found; notification skipped for ${emails}" >&2
fi
rm -f "$message_file"

if [[ "$strict" == "1" && "$status" != "READY" ]]; then
  echo "[notify] ${stage_key}: strict validation status=${status}; failing merge stage DAG" >&2
  exit 42
fi
EOS
  chmod +x "$notify_exec"

  cat > "$notify_sub" <<EOT
universe   = vanilla
executable = $notify_exec
output     = $OUT_DIR/notify.${stage_key}.\$(Cluster).\$(Process).out
error      = $ERR_DIR/notify.${stage_key}.\$(Cluster).\$(Process).err
log        = $LOG_DIR/notify.${stage_key}.\$(Cluster).\$(Process).log
getenv = True
should_transfer_files = NO
stream_output = True
stream_error  = True
notification = Never
arguments = $emails $stage_key $notify_body $dag $sub $OUT_DIR $ERR_DIR $LOG_DIR $validation_arg $email_mode $strict
queue
EOT

  cat > "$dag" <<EOT
JOB MERGE $sub
FINAL NOTIFY $notify_sub
EOT

  say "Stage tracking: $(merge_stage_tracking_description)"
  condor_submit_dag -notification Never "$dag"
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
    variantC|VariantC|variantc) echo "variantC" ;;
    variantD|VariantD|variantd) echo "variantD" ;;
    variantE|VariantE|variante) echo "variantE" ;;
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

  local -a jet_pts b2bs b2bs_submit vzs cones iso_base_tags uepipes
  local include_uepipe_in_tag=0
  mapfile -t jet_pts       < <(yaml_get_inline_list "$yaml" "jet_pt_min")
  mapfile -t b2bs          < <(yaml_get_inline_list "$yaml" "back_to_back_dphi_min_pi_fraction")
  mapfile -t vzs           < <(yaml_get_inline_list "$yaml" "vz_cut_cm")
  mapfile -t cones         < <(yaml_get_inline_list "$yaml" "coneR")
  mapfile -t iso_base_tags < <(build_iso_mode_tags_from_yaml "$yaml")
  local -a photon_id_rows
  mapfile -t photon_id_rows < <(yaml_get_photon_id_sets "$yaml")
  (( ${#photon_id_rows[@]} > 0 )) || { err "YAML must define photon_id_sets for cfg-tag generation: $yaml"; exit 41; }

  if (( ${#jet_pts[@]} == 0 )); then jet_pts=( "5.0" ); fi
  if (( ${#b2bs[@]} == 0 )); then b2bs=( "0.875" ); fi
  mapfile -t b2bs_submit < <(dphi_submit_values "${b2bs[@]}")
  if (( ${#vzs[@]} == 0 )); then vzs=( "30.0" ); fi
  if (( ${#cones[@]} == 0 )); then cones=( "0.30" ); fi
  if (( ${#iso_base_tags[@]} == 0 )); then iso_base_tags=( "fixedIso2GeV" ); fi

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
  local tight_norm nonTight_norm pre_norm

  # Optimized direct-fanout production writes one public cfg directory per
  # photon-ID triplet, with jet pT/dphi/iso/cone held as internal histogram
  # views. Emit those tags first so merge discovery accepts the new layout.
  local row
  for row in "${photon_id_rows[@]}"; do
    IFS='|' read -r pre tight nonTight <<< "$row"
    pre_norm="$(selection_mode_normalize_for_key "preselection" "$pre")"
    tight_norm="$(selection_mode_normalize_for_key "tight" "$tight")"
    nonTight_norm="$(selection_mode_normalize_for_key "nonTight" "$nonTight")"
    selection_tag="$(selection_mode_tag "preselection" "$pre_norm")_$(selection_mode_tag "tight" "$tight_norm")_$(selection_mode_tag "nonTight" "$nonTight_norm")"
    for uep in "${uepipes[@]}"; do
      if (( include_uepipe_in_tag )); then
        echo "${selection_tag}_${uep}${cfg_suffix}"
      else
        echo "${selection_tag}${cfg_suffix}"
      fi
    done
  done

  # Legacy scalar cfg tags remain accepted for allDirect / old productions.
  for pt in "${jet_pts[@]}"; do
    for frac in "${b2bs_submit[@]}"; do
      for vz in "${vzs[@]}"; do
        for cone in "${cones[@]}"; do
          for iso in "${iso_base_tags[@]}"; do
            for row in "${photon_id_rows[@]}"; do
              IFS='|' read -r pre tight nonTight <<< "$row"
              pre_norm="$(selection_mode_normalize_for_key "preselection" "$pre")"
              tight_norm="$(selection_mode_normalize_for_key "tight" "$tight")"
              nonTight_norm="$(selection_mode_normalize_for_key "nonTight" "$nonTight")"
              selection_tag="$(selection_mode_tag "preselection" "$pre_norm")_$(selection_mode_tag "tight" "$tight_norm")_$(selection_mode_tag "nonTight" "$nonTight_norm")"
              tag="jetMinPt$(sim_pt_tag "$pt")_$(dphi_dir_tag_component "$frac")_$(sim_vz_tag "$vz")_$(sim_cone_tag "$cone")_${iso}"
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
  done | sort -u
}

build_cfg_tags_from_manifest() {
  local root_base="$1"
  local manifest=""
  # Avoid letting an old .recoiljets_latest_manifest.txt silently constrain a
  # fresh production. Use a manifest only when the caller explicitly provides
  # one, or opts into the latest marker for a known matching production.
  if [[ -n "${MERGE_PRODUCTION_MANIFEST:-}" ]]; then
    manifest="$MERGE_PRODUCTION_MANIFEST"
  elif [[ "${RJ_MERGE_USE_LATEST_MANIFEST:-0}" == "1" ]]; then
    manifest="${root_base}/.recoiljets_latest_manifest.txt"
  else
    return 1
  fi
  [[ -s "$manifest" ]] || return 1
  awk -F= '$1 == "cfg_tag" && $2 != "" { print $2 }' "$manifest" | sort -u
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

  local _internal_dphi
  _internal_dphi="$(grep -E '^[[:space:]]*internal_back_to_back_dphi_min_pi_fraction:' "$_tmpyaml" | head -1 || true)"

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

  local _dphi_tag
  if [[ -n "$_internal_dphi" ]] && dphi_internal_enabled; then
    _dphi_tag="dphiScan"
  else
    _dphi_tag="$(sim_b2b_dir_tag "$_frac")"
  fi

  local _tag="jetMinPt$(sim_pt_tag "$_pt")_${_dphi_tag}_$(sim_vz_tag "$_vz")_$(sim_cone_tag "$_cone")_$(sim_iso_tag "$_sliding" "$_fixed")"
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
  if [[ -n "${MERGE_RUN_BASE_OVERRIDE:-}" ]]; then
    RUN_BASE="$MERGE_RUN_BASE_OVERRIDE"
  fi

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

scale_scaled_trig_after_hadd() {
  local rootfile="$1"
  [[ "${RJ_SCALED_TRIG_AFTER_HADD:-0}" == "1" ]] || return 0
  [[ -s "$rootfile" ]] || { echo "[scaledTrigQA stage1] skip: output is missing or empty: $rootfile"; return 0; }

  local base run run_num
  base="$(basename "$rootfile")"
  if [[ ! "$base" =~ ^chunkMerge_run_([0-9]+)\.root$ ]]; then
    echo "[scaledTrigQA stage1] skip: output is not a per-run chunkMerge file: $base"
    return 0
  fi

  run="${BASH_REMATCH[1]}"
  run_num=$((10#$run))

  local analysis_base="${RJ_ANALYSIS_BASE:-/sphenix/u/${USER}/scratch/thesisAnalysis}"
  local runlist="${RJ_SCALED_TRIG_RUNLIST:-${analysis_base}/dst_lists_auau/scaledEffRuns_MBD_NS_geq_2_vtx_lt_150__Pho10_12.list}"
  local config="${RJ_SCALED_TRIG_CONFIG_TXT:-${analysis_base}/dst_lists_auau/scaledEffConfig_MBD_NS_geq_2_vtx_lt_150__Pho10_12.txt}"

  if [[ ! -s "$runlist" ]]; then
    echo "[scaledTrigQA stage1] skip: missing run list: $runlist"
    return 0
  fi
  if [[ ! -s "$config" ]]; then
    echo "[scaledTrigQA stage1] skip: missing scale config: $config"
    return 0
  fi
  if ! grep -Eq "^0*${run_num}$" "$runlist"; then
    echo "[scaledTrigQA stage1] skip: run ${run} is not in scaled-trigger run list"
    return 0
  fi

  local scales baselineScale pho10Scale pho12Scale
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
    }' "$config")"

  if [[ -z "$scales" ]]; then
    echo "[scaledTrigQA stage1] skip: no scale config entry for run ${run_num}"
    return 0
  fi

  read -r baselineScale pho10Scale pho12Scale <<< "$scales"
  local numeric_re='^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$'
  if [[ ! "$baselineScale" =~ $numeric_re ]] || [[ ! "$pho10Scale" =~ $numeric_re ]] || [[ ! "$pho12Scale" =~ $numeric_re ]]; then
    echo "[scaledTrigQA stage1] skip: invalid scale config for run ${run_num}: $scales"
    return 0
  fi

  local macro
  macro="$(mktemp "${TMPDIR:-/tmp}/scaledTrigQA_stage1_${run_num}.XXXXXX.C")"
  cat > "$macro" <<EOF2
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TNamed.h>
#include <TObject.h>
#include <iostream>

{
  TFile f("${rootfile}", "UPDATE");
  std::cout << "[scaledTrigQA stage1] run=${run_num} file=${rootfile}\\n";
  std::cout << "[scaledTrigQA stage1] factors baseline=${baselineScale} pho10=${pho10Scale} pho12=${pho12Scale}\\n";
  if (!f.IsOpen() || f.IsZombie())
  {
    std::cerr << "[scaledTrigQA stage1] could not open file for UPDATE\\n";
  }
  else if (f.Get("scaledTrigQA_perRunCorrected_applied"))
  {
    std::cout << "[scaledTrigQA stage1] marker already present; leaving output unchanged\\n";
  }
  else
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

    int scaledTargets = 0;
    int foundTargets = 0;
    for (int i = 0; i < 3; ++i)
    {
      if (factors[i] <= 0.0) continue;
      TH1* h = dynamic_cast<TH1*>(f.Get(objPaths[i]));
      if (!h)
      {
        std::cout << "[scaledTrigQA stage1] target missing: " << objPaths[i] << "\\n";
        continue;
      }
      ++foundTargets;
      TDirectory* dir = f.GetDirectory(dirNames[i]);
      if (!dir)
      {
        std::cout << "[scaledTrigQA stage1] directory missing for target: " << dirNames[i] << "\\n";
        continue;
      }
      h->Scale(factors[i]);
      dir->cd();
      h->Write(histNames[i], TObject::kOverwrite);
      f.cd();
      ++scaledTargets;
      std::cout << "[scaledTrigQA stage1] scaled: " << objPaths[i] << " factor=" << factors[i] << "\\n";
    }

    if (scaledTargets > 0)
    {
      TNamed marker("scaledTrigQA_perRunCorrected_applied", "1");
      marker.Write("scaledTrigQA_perRunCorrected_applied", TObject::kOverwrite);
      std::cout << "[scaledTrigQA stage1] wrote marker scaledTrigQA_perRunCorrected_applied; scaledTargets=" << scaledTargets << "\\n";
    }
    else if (foundTargets == 0)
    {
      std::cout << "[scaledTrigQA stage1] no target histograms found; leaving output unchanged\\n";
    }
    else
    {
      std::cout << "[scaledTrigQA stage1] targets existed, but none were writable/scalable; leaving marker absent\\n";
    }
  }
  f.Close();
}
EOF2

  echo "[scaledTrigQA stage1] post-hadd check/update start: run=${run} out=${rootfile}"
  if ! root -l -b -q "$macro"; then
    echo "[scaledTrigQA stage1] WARNING: ROOT scaling command failed for ${rootfile}; keeping hadd output" >&2
  fi
  rm -f "$macro"
}

send_recoiljets_merge_notification() {
  [[ -n "${RJ_NOTIFY_EMAILS:-}" ]] || return 0
  [[ -n "${RJ_NOTIFY_MESSAGE:-}" ]] || return 0
  local subject="${RJ_NOTIFY_SUBJECT:-RecoilJets merge complete}"
  local body="${RJ_NOTIFY_MESSAGE}"
  if command -v mail >/dev/null 2>&1; then
    printf "%s\n" "$body" | mail -s "$subject" "$RJ_NOTIFY_EMAILS" || true
  elif command -v mailx >/dev/null 2>&1; then
    printf "%s\n" "$body" | mailx -s "$subject" "$RJ_NOTIFY_EMAILS" || true
  else
    echo "[notify] mail/mailx not found; notification skipped for ${RJ_NOTIFY_EMAILS}" >&2
  fi
}

echo "[hadd_condor] inputs=$(wc -l < "$LIST")  ->  $OUT"
hadd -v 3 -f "$OUT" @"$LIST"
scale_scaled_trig_after_hadd "$OUT"
send_recoiljets_merge_notification
EOS
  chmod +x "$exe"
}

emit_sim_stitch_wrapper() {
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
if [[ $# -ne 1 ]]; then
  echo "[ERROR] Usage: $0 <stitch_spec>" >&2
  exit 1
fi
SPEC="$1"
[[ -s "$SPEC" ]] || { echo "[ERROR] Missing/empty stitch spec: $SPEC" >&2; exit 2; }

cxx_quote() {
  local s="$1"
  s="${s//\\/\\\\}"
  s="${s//\"/\\\"}"
  printf '"%s"' "$s"
}

mapfile -t rows < "$SPEC"
(( ${#rows[@]} >= 5 )) || { echo "[ERROR] Stitch spec has too few rows: $SPEC" >&2; exit 3; }
dataset="${rows[0]}"
cfg="${rows[1]}"
out="${rows[2]}"
topdir="${rows[3]}"
nslices="${rows[4]}"
[[ "$nslices" =~ ^[0-9]+$ ]] || { echo "[ERROR] Invalid nslices in stitch spec: $nslices" >&2; exit 4; }
(( ${#rows[@]} == 5 + nslices )) || { echo "[ERROR] Stitch spec row count mismatch: rows=${#rows[@]} nslices=${nslices}" >&2; exit 5; }

for ((i=0; i<nslices; ++i)); do
  IFS='|' read -r label sigma path <<< "${rows[$((5+i))]}"
  [[ -n "$label" && -n "$sigma" && -n "$path" ]] || { echo "[ERROR] Bad slice row: ${rows[$((5+i))]}" >&2; exit 6; }
  [[ -s "$path" ]] || { echo "[ERROR] Missing/empty input slice: $path" >&2; exit 7; }
done

mkdir -p "$(dirname "$out")"
rm -f "$out"

macro_dir="$(mktemp -d "${TMPDIR:-/tmp}/recoiljets_sim_stitch_${dataset}_${cfg}.XXXXXX")"
macro="${macro_dir}/recoiljets_sim_stitch.C"
cat > "$macro" <<'ROOTMACRO'
#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TObject.h>
#include <TNamed.h>
#include <TH1.h>
#include <TTree.h>
#include <TSystem.h>

#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace RJSimStitch
{
using std::cout;
using std::endl;
using std::set;
using std::string;
using std::vector;

void EnsureParentDirForFile(const string& filepath)
{
  const size_t pos = filepath.find_last_of('/');
  if (pos == string::npos) return;
  const string dir = filepath.substr(0, pos);
  if (!dir.empty()) gSystem->mkdir(dir.c_str(), true);
}

double ReadEventCountFromFile(TFile* f, const string& topDirName)
{
  if (!f) return 0.0;
  TDirectory* d = f->GetDirectory(topDirName.c_str());
  if (!d) return 0.0;
  TH1* cnt = dynamic_cast<TH1*>(d->Get(("cnt_" + topDirName).c_str()));
  if (!cnt) return 0.0;
  return cnt->GetBinContent(1);
}

void AddScaledRecursive(TDirectory* outDir, TDirectory* inDir, double w)
{
  if (!outDir || !inDir) return;

  auto IsDirClass = [](const std::string& cls)->bool
  {
    return (cls == "TDirectoryFile" || cls == "TDirectory");
  };

  TIter next(inDir->GetListOfKeys());
  while (TKey* key = static_cast<TKey*>(next()))
  {
    const std::string name = key->GetName();
    const std::string cls = key->GetClassName();

    if (IsDirClass(cls))
    {
      TDirectory* subIn = dynamic_cast<TDirectory*>(inDir->Get(name.c_str()));
      if (!subIn) continue;

      outDir->cd();
      TDirectory* subOut = outDir->GetDirectory(name.c_str());
      if (!subOut) subOut = outDir->mkdir(name.c_str());
      if (!subOut) continue;

      AddScaledRecursive(subOut, subIn, w);
      continue;
    }

    TObject* objIn = inDir->Get(name.c_str());
    if (!objIn) continue;

    if (auto* hIn = dynamic_cast<TH1*>(objIn))
    {
      outDir->cd();
      TH1* hOut = dynamic_cast<TH1*>(outDir->Get(name.c_str()));

      if (!hOut)
      {
        TH1* hNew = dynamic_cast<TH1*>(hIn->Clone(name.c_str()));
        if (!hNew) continue;
        hNew->SetDirectory(outDir);
        if (hNew->GetSumw2N() == 0) hNew->Sumw2();
        if (w != 0.0) hNew->Scale(w);
        hNew->Write(name.c_str(), TObject::kOverwrite);
      }
      else
      {
        TH1* tmp = dynamic_cast<TH1*>(hIn->Clone((name + "_tmpAdd").c_str()));
        if (!tmp) continue;
        tmp->SetDirectory(nullptr);
        if (tmp->GetSumw2N() == 0) tmp->Sumw2();
        if (w != 0.0) tmp->Scale(w);
        hOut->Add(tmp);
        hOut->Write(name.c_str(), TObject::kOverwrite);
        delete tmp;
      }

      continue;
    }

    outDir->cd();
    if (!outDir->Get(name.c_str()))
    {
      TObject* c = objIn->Clone(name.c_str());
      if (c) c->Write(name.c_str(), TObject::kOverwrite);
    }
  }
}

bool BuildMergedSIMFile_PhotonSlices(const vector<string>& inFiles,
                                     const vector<double>& sigmas_pb,
                                     const string& outMerged,
                                     const string& topDirName,
                                     const vector<string>& sliceLabels)
{
  cout << "\n[MERGE SIM] Building merged SIM file with cross-section weights\n"
       << "  out    = " << outMerged << "\n"
       << "  topDir = " << topDirName << "\n";

  const size_t n = inFiles.size();
  if (n < 2 || sigmas_pb.size() != n)
  {
    cout << "[MERGE SIM][ERROR] Need >=2 inputs and a matching sigma list.\n";
    return false;
  }

  vector<TFile*> fin(n, nullptr);
  vector<TDirectory*> din(n, nullptr);
  vector<double> Nraw(n, 0.0);
  vector<double> w(n, 0.0);

  auto CloseAll = [&](){
    for (auto* f : fin) { if (f) f->Close(); }
  };

  for (size_t i = 0; i < n; ++i)
  {
    cout << "  in[" << i << "] = " << inFiles[i]
         << "   sigma_pb=" << std::setprecision(12) << sigmas_pb[i] << "\n";

    if (sigmas_pb[i] <= 0.0)
    {
      cout << "[MERGE SIM][ERROR] sigma_pb <= 0 for input: " << inFiles[i] << "\n";
      CloseAll();
      return false;
    }

    fin[i] = TFile::Open(inFiles[i].c_str(), "READ");
    if (!fin[i] || fin[i]->IsZombie())
    {
      cout << "[MERGE SIM][ERROR] Cannot open input SIM file: " << inFiles[i] << "\n";
      CloseAll();
      return false;
    }

    din[i] = fin[i]->GetDirectory(topDirName.c_str());
    if (!din[i])
    {
      cout << "[MERGE SIM][ERROR] Missing topDir '" << topDirName
           << "' in input SIM file: " << inFiles[i] << "\n";
      CloseAll();
      return false;
    }

    Nraw[i] = ReadEventCountFromFile(fin[i], topDirName);
    if (Nraw[i] <= 0.0)
    {
      cout << "[MERGE SIM][ERROR] Nraw <= 0 for input: " << inFiles[i] << "\n";
      CloseAll();
      return false;
    }
    w[i] = sigmas_pb[i] / Nraw[i];
  }

  double wRef = w.back();
  if (wRef <= 0.0)
  {
    for (size_t i = n; i > 0; --i)
    {
      if (w[i - 1] > 0.0)
      {
        wRef = w[i - 1];
        break;
      }
    }
  }
  if (wRef <= 0.0)
  {
    cout << "[MERGE SIM][ERROR] wRef <= 0; cannot normalize slice weights.\n";
    CloseAll();
    return false;
  }
  for (size_t i = 0; i < n; ++i) w[i] /= wRef;

  cout << "[MERGE SIM] Slice weights: w = (sigma_eff/Nraw)/(sigma_eff_ref/Nraw_ref)\n";
  for (size_t i = 0; i < n; ++i)
  {
    const string lab = (!sliceLabels.empty() && sliceLabels.size() == n) ? sliceLabels[i] : std::to_string(i);
    cout << "  [" << lab << "]  Nraw=" << std::fixed << std::setprecision(0) << Nraw[i]
         << "   sigma_pb=" << std::setprecision(12) << sigmas_pb[i]
         << "   w=" << std::setprecision(12) << w[i] << "\n";
  }

  EnsureParentDirForFile(outMerged);
  TFile* fout = TFile::Open(outMerged.c_str(), "RECREATE");
  if (!fout || fout->IsZombie())
  {
    cout << "[MERGE SIM][ERROR] Cannot create merged output file: " << outMerged << "\n";
    CloseAll();
    return false;
  }

  fout->cd();
  TDirectory* outTop = fout->mkdir(topDirName.c_str());
  if (!outTop)
  {
    cout << "[MERGE SIM][ERROR] Cannot create output topDir: " << topDirName << "\n";
    fout->Close();
    CloseAll();
    return false;
  }

  for (size_t i = 0; i < n; ++i) AddScaledRecursive(outTop, din[i], w[i]);

  TTree* tOutED = nullptr;
  for (size_t i = 0; i < n; ++i)
  {
    TTree* tInED = dynamic_cast<TTree*>(fin[i] ? fin[i]->Get("EventDisplayTree") : nullptr);
    if (!tInED && din[i]) tInED = dynamic_cast<TTree*>(din[i]->Get("EventDisplayTree"));
    if (!tInED) continue;

    if (!tOutED)
    {
      outTop->cd();
      tOutED = tInED->CloneTree(0);
      if (tOutED) tOutED->SetDirectory(outTop);
    }
    if (tOutED) tOutED->CopyEntries(tInED);
  }
  if (tOutED)
  {
    outTop->cd();
    tOutED->Write("EventDisplayTree", TObject::kOverwrite);
  }

  std::ostringstream oss;
  oss << "Merged photonJet slices. Nslices=" << n << " ";
  for (size_t i = 0; i < n; ++i)
  {
    const string lab = (!sliceLabels.empty() && sliceLabels.size() == n) ? sliceLabels[i] : std::to_string(i);
    oss << "[" << lab
        << " Nraw=" << std::fixed << std::setprecision(0) << Nraw[i]
        << " sigma_pb=" << std::setprecision(12) << sigmas_pb[i]
        << " w=" << std::setprecision(12) << w[i]
        << "] ";
  }

  outTop->cd();
  TNamed meta("MERGE_INFO", oss.str().c_str());
  meta.Write("MERGE_INFO", TObject::kOverwrite);

  fout->Write();
  fout->Close();
  CloseAll();

  cout << "[MERGE SIM] Done. Merged file written: " << outMerged << "\n";
  return true;
}
} // namespace RJSimStitch
ROOTMACRO
{
  printf '\n'
  printf 'void %s() {\n' "recoiljets_sim_stitch"
  printf '  std::vector<std::string> inputs = {'
  for ((i=0; i<nslices; ++i)); do
    IFS='|' read -r label sigma path <<< "${rows[$((5+i))]}"
    (( i > 0 )) && printf ', '
    cxx_quote "$path"
  done
  printf '};\n'
  printf '  std::vector<double> sigmas = {'
  for ((i=0; i<nslices; ++i)); do
    IFS='|' read -r label sigma path <<< "${rows[$((5+i))]}"
    (( i > 0 )) && printf ', '
    printf '%s' "$sigma"
  done
  printf '};\n'
  printf '  std::vector<std::string> labels = {'
  for ((i=0; i<nslices; ++i)); do
    IFS='|' read -r label sigma path <<< "${rows[$((5+i))]}"
    (( i > 0 )) && printf ', '
    cxx_quote "$label"
  done
  printf '};\n'
  printf '  const std::string out = '
  cxx_quote "$out"
  printf ';\n'
  printf '  const std::string topdir = '
  cxx_quote "$topdir"
  printf ';\n'
  printf '  const bool ok = RJSimStitch::BuildMergedSIMFile_PhotonSlices(inputs, sigmas, out, topdir, labels);\n'
  printf '  if (!ok) gSystem->Exit(10);\n'
  printf '}\n'
} >> "$macro"

echo "[sim_stitch] dataset=${dataset} cfg=${cfg} nslices=${nslices}"
echo "[sim_stitch] out=${out}"
root -l -b -q "$macro"
rc=$?
rm -rf "$macro_dir"
exit "$rc"
EOS
  chmod +x "$exe"
}

sim_stitch_plan_for_dataset() {
  local token="$1"
  case "$token" in
    isSim|sim|SIM)
      SIM_STITCH_COMBO_DIR="photonJet5and10and20merged_SIM"
      SIM_STITCH_OUTPUT_FILE="RecoilJets_photonjet5plus10plus20_MERGED.root"
      SIM_STITCH_TOPDIR="SIM"
      SIM_STITCH_ROWS=(
        "photonJet5|146359.3|photonjet5"
        "photonJet10|6944.675|photonjet10"
        "photonJet20|130.4461|photonjet20"
      )
      ;;
    isSimInclusive|issiminclusive|siminclusive|SIMINCLUSIVE|isSimJet5|isSimjet5|simjet5|SIMJET5)
      SIM_STITCH_COMBO_DIR="inclusiveJet5to40_SIM"
      SIM_STITCH_OUTPUT_FILE="RecoilJets_jet5plus8plus12plus20plus30plus40_MERGED.root"
      SIM_STITCH_TOPDIR="SIM"
      SIM_STITCH_ROWS=(
        "jet5|1.3878e8|jet5"
        "jet8|1.15e7|jet8"
        "jet12|1.4903e6|jet12"
        "jet20|6.2623e4|jet20"
        "jet30|2.5298e3|jet30"
        "jet40|1.3553e2|jet40"
      )
      ;;
    isSimEmbedded|issimembedded|simembedded|SIMEMBEDDED)
      SIM_STITCH_COMBO_DIR="photonJet12and20merged_SIM"
      SIM_STITCH_OUTPUT_FILE="RecoilJets_embeddedPhoton12plus20_MERGED.root"
      SIM_STITCH_TOPDIR="SIM"
      SIM_STITCH_ROWS=(
        "embeddedPhoton12|2598.12425|embeddedPhoton12"
        "embeddedPhoton20|133.317866|embeddedPhoton20"
      )
      ;;
    isSimEmbeddedInclusive|issimembeddedinclusive|simembeddedinclusive|SIMEMBEDDEDINCLUSIVE)
      SIM_STITCH_COMBO_DIR="embeddedJet12and20merged_SIM"
      SIM_STITCH_OUTPUT_FILE="RecoilJets_embeddedJet12plus20_MERGED.root"
      SIM_STITCH_TOPDIR="SIM"
      SIM_STITCH_ROWS=(
        "embeddedJet12|1.21692467e6|embeddedJet12"
        "embeddedJet20|5.56198698e4|embeddedJet20"
      )
      ;;
    *)
      return 1
      ;;
  esac
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
if [[ "${1}" =~ ^(isSim|sim|SIM|isSimJet5|isSimjet5|isSimInclusive|issiminclusive|simjet5|siminclusive|SIMJET5|SIMINCLUSIVE|isSimMB|simmb|SIMMB|isSimEmbedded|issimembedded|simembedded|SIMEMBEDDED|isSimEmbeddedInclusive|issimembeddedinclusive|simembeddedinclusive|SIMEMBEDDEDINCLUSIVE)$ ]]; then
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
      SIM_INPUT_BASE="/sphenix/tg/tg01/bulk/jbennett/thesisAna/siminclusive"
      SIM_OUTPUT_TAG="siminclusive"
      [[ "$SIM_SAMPLE_EXPLICIT" -eq 0 ]] && SIM_SAMPLE="run28_jet5"
      ;;
    isSimInclusive|issiminclusive|siminclusive|SIMINCLUSIVE)
      SIM_INPUT_BASE="/sphenix/tg/tg01/bulk/jbennett/thesisAna/siminclusive"
      SIM_OUTPUT_TAG="siminclusive"
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
      SIM_INPUT_BASE="/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembeddedinclusive"
      SIM_OUTPUT_TAG="simembeddedinclusive"
      [[ "$SIM_SAMPLE_EXPLICIT" -eq 0 ]] && SIM_SAMPLE="run28_embeddedJet20"
      ;;
    *)
      SIM_INPUT_BASE="/sphenix/tg/tg01/bulk/jbennett/thesisAna/sim"
      SIM_OUTPUT_TAG="sim"
      ;;
  esac
  if [[ -n "${MERGE_SIM_INPUT_BASE_OVERRIDE:-}" ]]; then
    SIM_INPUT_BASE="$MERGE_SIM_INPUT_BASE_OVERRIDE"
  fi

  # Embedded firstRound hadd jobs can exceed the generic 2 GB request when
  # merging larger Au+Au-style ROOT outputs. Keep pp/MB/jet SIM unchanged.
  # RJ_SIM_FIRSTROUND_REQUEST_MEMORY is intentionally an operational override:
  # it lets a held merge stage be retried without rerunning the DST analysis.
  SIM_FIRSTROUND_REQUEST_MEMORY="2GB"
  case "$SIM_DATASET_TOKEN" in
    isSimEmbedded|issimembedded|simembedded|SIMEMBEDDED|isSimEmbeddedInclusive|issimembeddedinclusive|simembeddedinclusive|SIMEMBEDDEDINCLUSIVE)
      SIM_FIRSTROUND_REQUEST_MEMORY="4GB"
      ;;
  esac
  if [[ -n "${RJ_SIM_FIRSTROUND_REQUEST_MEMORY:-}" ]]; then
    SIM_FIRSTROUND_REQUEST_MEMORY="${RJ_SIM_FIRSTROUND_REQUEST_MEMORY}"
  fi

  # Build sample list (same defaults as submit script, or explicit SAMPLE=)
  samples=()
  if [[ "${SIM_SAMPLE_EXPLICIT:-0}" -eq 0 ]]; then
    case "$SIM_DATASET_TOKEN" in
      isSimJet5|isSimjet5|simjet5|SIMJET5|isSimInclusive|issiminclusive|siminclusive|SIMINCLUSIVE)
        samples=( "run28_jet5" "run28_jet8" "run28_jet12" "run28_jet20" "run28_jet30" "run28_jet40" )
        ;;
      isSimMB|simmb|SIMMB)       samples=( "run28_detroit" ) ;;
      isSimEmbedded|issimembedded|simembedded|SIMEMBEDDED) samples=( "run28_embeddedPhoton12" "run28_embeddedPhoton20" ) ;;
      isSimEmbeddedInclusive|issimembeddedinclusive|simembeddedinclusive|SIMEMBEDDEDINCLUSIVE) samples=( "run28_embeddedJet12" "run28_embeddedJet20" ) ;;
      *)                         samples=( "run28_photonjet5" "run28_photonjet10" "run28_photonjet20" ) ;;
    esac
  else
    samples=( "${SIM_SAMPLE}" )
  fi

  # ---------------------------------------------------------------
  # DISCOVER cfg_tag directories from the filesystem.
  #
  # IMPORTANT:
  #   isSimEmbedded and isSimEmbeddedInclusive use separate TG input bases.
  #   We still filter cfg_tag directories by requested sample subdirectories
  #   so stale or manually copied trees cannot be merged into the wrong dataset.
  # ---------------------------------------------------------------
  FLAT_OUT_DIR="${OUT_BASE}/${SIM_OUTPUT_TAG}"
  mkdir -p "$FLAT_OUT_DIR" "$LOG_DIR" "$OUT_DIR" "$ERR_DIR" "$TMP_DIR"

  if [[ "$SIM_ACTION" == "firstRound" ]]; then
    _discover_base="$SIM_INPUT_BASE"
  elif [[ "$SIM_ACTION" == "secondRound" || "$SIM_ACTION" == "finalStitch" ]]; then
    # Keep discovery anchored to the TG input base, but filter cfg_tags below
    # by the requested sample directories for this dataset token.
    _discover_base="$SIM_INPUT_BASE"
  else
    err "Unknown isSim action '${SIM_ACTION}'. Allowed: firstRound | secondRound | finalStitch"
    exit 2
  fi

  mapfile -t _ALL_SIM_CFG_TAGS < <(
    find "$_discover_base" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort -V
  )

  _CFG_SOURCE="manifest"
  mapfile -t _YAML_SIM_CFG_TAGS < <(build_cfg_tags_from_manifest "$_discover_base" || true)
  if (( ${#_YAML_SIM_CFG_TAGS[@]} == 0 )); then
    _CFG_SOURCE="YAML"
    mapfile -t _YAML_SIM_CFG_TAGS < <(build_cfg_tags_from_yaml "$SIM_DATASET_TOKEN")
  fi

  SIM_CFG_TAGS=()
  for _cfg in "${_ALL_SIM_CFG_TAGS[@]}"; do
    _keep_cfg=0

    # First require that the cfg_tag is part of the production manifest, or
    # the current YAML matrix for older productions without a manifest.
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
    err "No cfg_tag subdirectories found under ${_discover_base} that match the ${_CFG_SOURCE} cfg-tag list for requested samples: ${samples[*]}"
    [[ "$_CFG_SOURCE" == "YAML" ]] && err "YAML used: $(merge_yaml_path)"
    [[ "$_CFG_SOURCE" == "manifest" ]] && err "Manifest used: ${MERGE_PRODUCTION_MANIFEST:-${_discover_base}/.recoiljets_latest_manifest.txt}"
    err "Run RecoilJets_Condor_submit.sh ${SIM_DATASET_TOKEN} condorDoAll first to create the output tree."
    exit 20
  fi

  say "Dataset : ${BOLD}${SIM_DATASET_TOKEN}${RST}"
  say "Action  : ${SIM_ACTION}"
  say "Input base  : ${SIM_INPUT_BASE}"
  say "Output base : ${FLAT_OUT_DIR}"
  say "groupSize   : ${SIM_GROUP_SIZE}"
  say "Samples     : ${samples[*]}"
  say "${_CFG_SOURCE} cfg_tags allowed: ${#_YAML_SIM_CFG_TAGS[@]}"
  say "Discovered cfg_tags (filtered by ${_CFG_SOURCE} + requested samples): ${#SIM_CFG_TAGS[@]}"
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

    if [[ "$SIM_ACTION" == "finalStitch" ]]; then
      sim_stitch_plan_for_dataset "$SIM_DATASET_TOKEN" || {
        err "No SIM final-stitch plan is configured for dataset token: ${SIM_DATASET_TOKEN}"
        exit 41
      }

      need_cmd condor_submit
      emit_sim_stitch_wrapper "$CONDOR_EXEC"

      STITCH_SUB="${TMP_DIR}/recoil_sim_${cfg_tag}_finalStitch.sub"
      STITCH_ARGS="${TMP_DIR}/recoil_sim_${cfg_tag}_finalStitch.args"
      STITCH_EXPECTED="${TMP_DIR}/recoil_sim_${cfg_tag}_finalStitch.expected"
      STITCH_SPEC="${TMP_DIR}/recoil_sim_${cfg_tag}_finalStitch.spec"
      STITCH_OUT="${DEST_DIR}/${SIM_STITCH_COMBO_DIR}/${SIM_STITCH_OUTPUT_FILE}"

      : > "$STITCH_SPEC"
      printf '%s\n' "$SIM_DATASET_TOKEN" >> "$STITCH_SPEC"
      printf '%s\n' "$cfg_tag" >> "$STITCH_SPEC"
      printf '%s\n' "$STITCH_OUT" >> "$STITCH_SPEC"
      printf '%s\n' "$SIM_STITCH_TOPDIR" >> "$STITCH_SPEC"
      printf '%s\n' "${#SIM_STITCH_ROWS[@]}" >> "$STITCH_SPEC"

      _stitch_missing=0
      for _row in "${SIM_STITCH_ROWS[@]}"; do
        IFS='|' read -r _label _sigma _sample_tag <<< "$_row"
        _sample_final="${FLAT_OUT_DIR}/${FINAL_PREFIX}_${_sample_tag}_ALL_${cfg_tag}.root"
        if [[ ! -s "$_sample_final" ]]; then
          warn "Missing sample-level secondRound input for finalStitch: ${_sample_final}"
          _stitch_missing=1
        fi
        printf '%s|%s|%s\n' "$_label" "$_sigma" "$_sample_final" >> "$STITCH_SPEC"
      done
      (( _stitch_missing == 0 )) || {
        err "Cannot submit finalStitch for cfg=${cfg_tag}; one or more sample-level secondRound files are missing."
        exit 42
      }

      printf '%s\n' "$STITCH_SPEC" > "$STITCH_ARGS"
      printf '%s\n' "$STITCH_OUT" > "$STITCH_EXPECTED"

      cat > "$STITCH_SUB" <<EOT
universe   = vanilla
executable = $CONDOR_EXEC
output     = $OUT_DIR/recoil.sim.${cfg_tag}.finalStitch.\$(Cluster).\$(Process).out
error      = $ERR_DIR/recoil.sim.${cfg_tag}.finalStitch.\$(Cluster).\$(Process).err
log        = $LOG_DIR/recoil.sim.${cfg_tag}.finalStitch.\$(Cluster).\$(Process).log
$(condor_auto_memory_retry_block "${RJ_SIM_FINAL_STITCH_REQUEST_MEMORY:-6GB}")
priority = $MERGE_CONDOR_PRIORITY
getenv = True
should_transfer_files = NO
stream_output = True
stream_error  = True
notification = Never
queue arguments from ${STITCH_ARGS}
EOT

      say "SIM finalStitch: cfg=${cfg_tag} samples=${#SIM_STITCH_ROWS[@]} -> ${STITCH_OUT}"
      submit_condor_stage_with_ready_email \
        "$STITCH_SUB" \
        "sim_finalStitch_${cfg_tag}" \
        "RecoilJets_SIM_finalStitch_ready" \
        "SIM final weighted stitch validation finished for cfg=${cfg_tag}. If status is READY, canonical combined ROOT output: ${STITCH_OUT}." \
        "$STITCH_EXPECTED"
      echo
      continue
    fi

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
          ARGS="${TMP_DIR}/recoil_sim_${cfg_tag}_${SIM_TAG}_firstRound.args"
          EXPECTED="${TMP_DIR}/recoil_sim_${cfg_tag}_${SIM_TAG}_firstRound.expected"
          rm -f "$SUB"
          : > "$ARGS"
          : > "$EXPECTED"
          cat > "$SUB" <<EOT
universe   = vanilla
executable = $CONDOR_EXEC
output     = $OUT_DIR/recoil.sim.${cfg_tag}.${SIM_TAG}.\$(Cluster).\$(Process).out
error      = $ERR_DIR/recoil.sim.${cfg_tag}.${SIM_TAG}.\$(Cluster).\$(Process).err
log        = $LOG_DIR/recoil.sim.${cfg_tag}.${SIM_TAG}.\$(Cluster).\$(Process).log
$(condor_auto_memory_retry_block "$SIM_FIRSTROUND_REQUEST_MEMORY")
priority = $MERGE_CONDOR_PRIORITY
getenv = True
should_transfer_files = NO
stream_output = True
stream_error  = True
notification = Never
queue arguments from ${ARGS}
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

            printf '%s %s\n' "$listfile" "$out" >> "$ARGS"
            printf '%s\n' "$out" >> "$EXPECTED"
          done

          if [[ ! -s "$ARGS" ]]; then
            err "firstRound: args file is empty/unwritten → ${ARGS}"
            exit 27
          fi

          say "Submitting ${BOLD}${grp}${RST} firstRound Condor merge jobs → $(basename "$SUB")"
          submit_condor_stage_with_ready_email \
            "$SUB" \
            "sim_firstRound_${cfg_tag}_${SIM_TAG}" \
            "RecoilJets_SIM_firstRound_ready" \
            "SIM firstRound merge is complete for cfg=${cfg_tag}, sample=${SIM_TAG}. You can now run the local secondRound merge for this sample/cfg set." \
            "$EXPECTED"
          say "FirstRound submitted. Stage tracking: $(merge_stage_tracking_description). Partials will appear under: ${DEST_DIR}"
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
          EXPECTED="${TMP_DIR}/recoil_sim_${cfg_tag}_${SIM_TAG}_secondRound.expected"
          rm -f "$SUB"
          printf '%s\n' "$SIM_FINAL" > "$EXPECTED"
          cat > "$SUB" <<EOT
universe   = vanilla
executable = $CONDOR_EXEC
output     = $OUT_DIR/recoil.sim.${cfg_tag}.${SIM_TAG}.final.\$(Cluster).\$(Process).out
error      = $ERR_DIR/recoil.sim.${cfg_tag}.${SIM_TAG}.final.\$(Cluster).\$(Process).err
log        = $LOG_DIR/recoil.sim.${cfg_tag}.${SIM_TAG}.final.\$(Cluster).\$(Process).log
$(condor_auto_memory_retry_block "6GB")
priority = $MERGE_CONDOR_PRIORITY
getenv = True
should_transfer_files = NO
stream_output = True
stream_error  = True
notification = Never
arguments = $LIST $SIM_FINAL
queue
EOT
          say "Submitting secondRound final merge on Condor → $(basename "$SUB")"
          submit_condor_stage_with_ready_email \
            "$SUB" \
            "sim_secondRound_${cfg_tag}_${SIM_TAG}" \
            "RecoilJets_SIM_secondRound_ready" \
            "SIM secondRound final merge is complete for cfg=${cfg_tag}, sample=${SIM_TAG}. Final ROOT output: ${SIM_FINAL}. You can now pull/plot this sample or continue the combined local SIM merger." \
            "$EXPECTED"
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
  _CFG_SOURCE="manifest"
  mapfile -t _YAML_DATA_CFG_TAGS < <(build_cfg_tags_from_manifest "$RUN_BASE" || true)
  if (( ${#_YAML_DATA_CFG_TAGS[@]} == 0 )); then
    _CFG_SOURCE="YAML"
    mapfile -t _YAML_DATA_CFG_TAGS < <(build_cfg_tags_from_yaml "$TAG")
  fi
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
    err "No discovered cfg_tag directories under ${RUN_BASE} match the ${_CFG_SOURCE} cfg-tag list."
    [[ "$_CFG_SOURCE" == "YAML" ]] && err "YAML used: $(merge_yaml_path)"
    [[ "$_CFG_SOURCE" == "manifest" ]] && err "Manifest used: ${MERGE_PRODUCTION_MANIFEST:-${RUN_BASE}/.recoiljets_latest_manifest.txt}"
    exit 21
  fi

  say "${_CFG_SOURCE} cfg_tags allowed: ${#_YAML_DATA_CFG_TAGS[@]}"
  say "Discovered cfg_tags (filtered by ${_CFG_SOURCE}): ${#DATA_CFG_TAGS[@]}"
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
    ARGS="${TMP_DIR}/recoil_partials_${TAG}_${_cfg:-flat}.args"
    EXPECTED="${TMP_DIR}/recoil_partials_${TAG}_${_cfg:-flat}.expected"
    rm -f "$SUB"
    : > "$ARGS"
    : > "$EXPECTED"
    cat > "$SUB" <<EOT
universe   = vanilla
executable = $CONDOR_EXEC
output     = $OUT_DIR/recoil.\$(Cluster).\$(Process).out
error      = $ERR_DIR/recoil.\$(Cluster).\$(Process).err
log        = $LOG_DIR/recoil.\$(Cluster).\$(Process).log
$(condor_auto_memory_retry_block "2GB")
priority = $MERGE_CONDOR_PRIORITY
getenv = True
should_transfer_files = NO
stream_output = True
stream_error  = True
notification = Never
EOT
    if [[ "$TAG" == "auau" ]]; then
      cat >> "$SUB" <<EOT
environment = "RJ_SCALED_TRIG_AFTER_HADD=1 RJ_ANALYSIS_BASE=${BASE} RJ_SCALED_TRIG_RUNLIST=${SCALED_TRIG_RUNLIST} RJ_SCALED_TRIG_CONFIG_TXT=${SCALED_TRIG_CONFIG_TXT}"
EOT
    fi
    cat >> "$SUB" <<EOT
queue arguments from ${ARGS}
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

      printf '%s %s\n' "$listfile" "$out" >> "$ARGS"
      printf '%s\n' "$out" >> "$EXPECTED"
      ((queued+=1))
    done

    (( SKIP_TRACE )) && say "Stage-1 queue summary (cfg=${_cfg:-flat}): queued=${queued}, skipped=${skipped}"

    if (( queued == 0 )); then
      warn "No Condor jobs to submit for cfg_tag '${_cfg:-flat}' (no non-empty eligible run lists)."
      continue
    fi

    say "Submitting ${BOLD}${queued}${RST} Condor merge jobs → $(basename "$SUB")"
    submit_condor_stage_with_ready_email \
      "$SUB" \
      "data_perRun_${TAG}_${_cfg:-flat}" \
      "RecoilJets_${TAG}_perRun_ready" \
      "Data per-run merge is complete for dataset=${TAG}, cfg=${_cfg:-flat}. AuAu scaled-trigger corrections, when configured, have been applied inside the per-run hadd job. You can now run the sliceRuns merge stage." \
      "$EXPECTED"
    say "Stage-1 submitted. Stage tracking: $(merge_stage_tracking_description). Partial outputs will appear under: ${PERRUN_CFG_DIR}"
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
    # Filter discovered perRun cfg_tags against the submit-time manifest when
    # available; fall back to the current YAML for older productions.
    _CFG_SOURCE="manifest"
    mapfile -t _YAML_MERGE_CFG_TAGS < <(build_cfg_tags_from_manifest "$RUN_BASE" || true)
    if (( ${#_YAML_MERGE_CFG_TAGS[@]} == 0 )); then
      _CFG_SOURCE="YAML"
      mapfile -t _YAML_MERGE_CFG_TAGS < <(build_cfg_tags_from_yaml "$TAG")
    fi
    MERGE_CFG_TAGS=()
    for _ct in "${_ALL_MERGE_CFG_TAGS[@]}"; do
      for _y in "${_YAML_MERGE_CFG_TAGS[@]}"; do
        [[ "$_ct" == "$_y" ]] && { MERGE_CFG_TAGS+=( "$_ct" ); break; }
      done
    done
    say "${_CFG_SOURCE} cfg_tags allowed: ${#_YAML_MERGE_CFG_TAGS[@]}"
    say "Discovered cfg_tags (filtered by ${_CFG_SOURCE}): ${#MERGE_CFG_TAGS[@]}"
    for _ct in "${MERGE_CFG_TAGS[@]}"; do say "  ${_ct}"; done
    if (( ${#MERGE_CFG_TAGS[@]} == 0 )); then
      die "No perRun cfg_tag directories under ${_PERRUN_BASE} match the ${_CFG_SOURCE} cfg-tag list."
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

    if [[ "$TAG" == "auau" && "$SCALED_TRIG_ADDCHUNKS_APPLY" == "1" ]]; then
      warn "scaledTrigQA: SCALED_TRIG_ADDCHUNKS_APPLY=1; applying legacy addChunks pre-pass to existing partials"
      scale_scaled_trig_qa_partials "$_partial_dir"
    elif [[ "$TAG" == "auau" ]]; then
      say "scaledTrigQA: skipping addChunks pre-pass; future per-run Condor hadd jobs apply this after stage-1 hadd"
      say "scaledTrigQA: for old partials that still need correction, rerun with SCALED_TRIG_ADDCHUNKS_APPLY=1"
    fi

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
      ARGS="${TMP_DIR}/recoil_slice_${TAG}_${_cfg:-flat}.args"
      EXPECTED="${TMP_DIR}/recoil_slice_${TAG}_${_cfg:-flat}.expected"
      rm -f "$SUB"
      : > "$ARGS"
      : > "$EXPECTED"
      cat > "$SUB" <<EOT
universe   = vanilla
executable = $CONDOR_EXEC
output     = $OUT_DIR/recoil.slice.\$(Cluster).\$(Process).out
error      = $ERR_DIR/recoil.slice.\$(Cluster).\$(Process).err
log        = $LOG_DIR/recoil.slice.\$(Cluster).\$(Process).log
$(condor_auto_memory_retry_block "4GB")
priority = $MERGE_CONDOR_PRIORITY
getenv = True
should_transfer_files = NO
stream_output = True
stream_error  = True
notification = Never
queue arguments from ${ARGS}
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
        printf '%s %s\n' "$_list" "$_out" >> "$ARGS"
        printf '%s\n' "$_out" >> "$EXPECTED"

        (( _offset = _end ))
      done

      say "Submitting ${_ngroups} slice merge jobs → $(basename "$SUB")"
      submit_condor_stage_with_ready_email \
        "$SUB" \
        "sliceRuns_${TAG}_${_cfg:-flat}" \
        "RecoilJets_${TAG}_sliceRuns_ready" \
        "Data sliceRuns merge is complete for dataset=${TAG}, cfg=${_cfg:-flat}. You can now run the final local addChunks merge." \
        "$EXPECTED"
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
      EXPECTED="${TMP_DIR}/recoil_final_${TAG}_${_cfg:-flat}.expected"
      rm -f "$SUB"
      printf '%s\n' "$FINAL" > "$EXPECTED"
      cat > "$SUB" <<EOT
universe   = vanilla
executable = $CONDOR_EXEC
output     = $OUT_DIR/recoil.final.\$(Cluster).\$(Process).out
error      = $ERR_DIR/recoil.final.\$(Cluster).\$(Process).err
log        = $LOG_DIR/recoil.final.\$(Cluster).\$(Process).log
$(condor_auto_memory_retry_block "3GB")
priority = $MERGE_CONDOR_PRIORITY
getenv = True
should_transfer_files = NO
stream_output = True
stream_error  = True
notification = Never
arguments = $LIST $FINAL
queue
EOT
      say "Submitting final merge on Condor → $(basename "$SUB")"
      submit_condor_stage_with_ready_email \
        "$SUB" \
        "data_final_${TAG}_${_cfg:-flat}" \
        "RecoilJets_${TAG}_final_ready" \
        "Data final addChunks merge is complete for dataset=${TAG}, cfg=${_cfg:-flat}. Final ROOT output: ${FINAL}. You can now pull the output locally and make plots." \
        "$EXPECTED"
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
