#!/usr/bin/env bash
###############################################################################
# RecoilJets_Condor.sh  —  Condor executable (per-job wrapper)
#   • Runs the Fun4All macro Fun4All_recoilJets.C on a listfile chunk.
#   • Automatically sets the data type (isPP / isAuAu) via env RJ_DATASET.
#   • Writes outputs under the dataset-specific base path, in a subdir
#     named by the run number, with a filename derived from the chunk list.
#
# Usage (from Condor submit file or locally):
#   RecoilJets_Condor.sh <run8> <chunkList> <isPP|isAuAu> <Cluster|LOCAL> <nEvents> <chunkIdx> NONE <destBase>
#
# Example:
#   RecoilJets_Condor.sh 00048721 /.../condor_lists/pp/run00048721_grp001.list isPP 12345 0 1 NONE /sphenix/tg/tg01/bulk/jbennett/thesisAna/pp
###############################################################################
set -euo pipefail

# ------------------------ Arguments ------------------------
run8="${1:?run8 (DATA) or SIM_SAMPLE (SIM) required}"
chunk_list="${2:?chunk list (.list) required}"
dataset_raw="${3:-isPP}"            # isPP | isAuAu | isSim
cluster_id="${4:-LOCAL}"            # informational
nevents="${5:-0}"                   # 0 → all events in the list
chunk_idx="${6:-0}"                 # informational; not used for naming
_ignored2="${7:-}"                  # keep slot for compatibility (NONE)
dest_base="${8:-}"                  # If empty, we derive from dataset

# ------------------------ Fixed paths ----------------------
BASE="/sphenix/u/patsfan753/scratch/thesisAnalysis"
MACRO="${BASE}/macros/Fun4All_recoilJets.C"

# Condor logging (handled by submit file), but we echo as well
echo "====================================================================="
echo "[INFO] RecoilJets_Condor.sh starting"
echo "       Host: $(hostname -f)"
echo "       CWD : $(pwd)"
echo "       Args: run=$run8  chunk=$(basename "$chunk_list")  dataset=$dataset_raw  nevents=$nevents"
echo "====================================================================="

# ------------------------ Environment ----------------------
export USER="${USER:-$(id -u -n)}"
export LOGNAME="${LOGNAME:-$USER}"
export HOME="/sphenix/u/${LOGNAME}"

# sPHENIX offline setup (system + local area for custom libs)
MYINSTALL="/sphenix/u/${USER}/thesisAnalysis/install"

# Disable 'nounset' while sourcing env scripts; they may read unset vars (e.g. PGHOST)
set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n
if [[ -d "$MYINSTALL" ]]; then
  # do not fail if local area is not present; macro has R__LOAD_LIBRARY with absolute path
  source /opt/sphenix/core/bin/setup_local.sh "$MYINSTALL" || true
fi
set -u

# ------------------------ Dataset routing ------------------
# Normalize dataset and set defaults:
#  - isSim must remain isSim end-to-end so the analysis module can detect it.
#  - Fun4All macro will treat isSim as pp-style reconstruction internally.
analysis_tag="isPP"
case "$dataset_raw" in
  isPP|pp|PP)
    dataset="isPP"
    analysis_tag="isPP"
    export RJ_DATASET="isPP"
    export RJ_IS_SIM=0
    ;;
  isPPrun25|pprun25|pp25|PP25)
    dataset="isPPrun25"
    analysis_tag="isPP"
    export RJ_DATASET="isPPrun25"
    export RJ_IS_SIM=0
    ;;
  isSim|sim|SIM)
    dataset="isSim"
    analysis_tag="isSim"
    export RJ_DATASET="isSim"
    export RJ_IS_SIM=1
    ;;
  isAuAu|auau|AA)
    echo "[FATAL] RecoilJets_Condor.sh is pp-style only. Use RecoilJets_Condor_AuAu.sh for isAuAu."
    exit 50
    ;;
  *)
    echo "[WARN] Unknown dataset '$dataset_raw' → defaulting to 'isPP'"
    dataset="isPP"
    analysis_tag="isPP"
    export RJ_DATASET="isPP"
    export RJ_IS_SIM=0
    ;;
esac

# Destination base (if not supplied as arg 8)
if [[ -z "$dest_base" ]]; then
  if [[ "$analysis_tag" == "isSim" ]]; then
    dest_base="/sphenix/tg/tg01/bulk/jbennett/thesisAna/sim"
  elif [[ "$dataset" == "isPPrun25" ]]; then
    dest_base="/sphenix/tg/tg01/bulk/jbennett/thesisAna/pp25"
  elif [[ "$dataset" == "isPP" ]]; then
    dest_base="/sphenix/tg/tg01/bulk/jbennett/thesisAna/pp"
  else
    dest_base="/sphenix/tg/tg01/bulk/jbennett/thesisAna/auau"
  fi
fi

# ------------------------ Paths & naming -------------------
# Output directory (one folder per run)
out_dir="${dest_base}/${run8}"
mkdir -p "$out_dir"

# The output file name follows the chunk list name (group name) for consistency
chunk_base="$(basename "$chunk_list")"               # e.g. run00048721_grp001.list
chunk_tag="${chunk_base%.list}"                      # e.g. run00048721_grp001
out_root="${out_dir}/RecoilJets_${analysis_tag}_${chunk_tag}.root"

echo "[INFO] Output path = $out_root"

# The macro expects to find the ROOT files listed inside chunk_list. If those
# entries are relative file names, make sure the CWD is the directory where
# the ROOT files live.
if [[ "$dataset" == "isPPrun25" ]]; then
  run_dec=$((10#$run8))
  range_begin=$(( (run_dec/100) * 100 ))
  range_end=$(( (run_dec/100 + 1) * 100 ))
  prod_dir="$(printf "/sphenix/lustre01/sphnxpro/production2/run3pp/physics/calofitting/new_newcdbtag_v008/run_%08d_%08d" "$range_begin" "$range_end")"
  cd "$prod_dir" || { echo "[FATAL] Cannot cd to CALOFITTING production dir: $prod_dir"; exit 5; }
else
  list_dir="$(dirname "$chunk_list")"
  cd "$list_dir"
fi

# ------------------------ Sanity checks --------------------
[[ -f "$MACRO" ]] || { echo "[FATAL] Macro not found: $MACRO"; exit 2; }
[[ -f "$chunk_list" ]] || { echo "[FATAL] Chunk list not found: $chunk_list"; exit 3; }
if [[ ! -s "$chunk_list" ]]; then
  echo "[FATAL] Chunk list is empty: $chunk_list"
  exit 4
fi

# ------------------------ Run ROOT macro -------------------
set +e
echo "[INFO] Running ROOT:"
echo "root -b -q -l \"${MACRO}(${nevents}, \\\"${chunk_list}\\\", \\\"${out_root}\\\", false)\""
root -b -q -l "${MACRO}(${nevents}, \"${chunk_list}\", \"${out_root}\", false)"
rc=$?
set -e

echo "---------------------------------------------------------------------"
if (( rc != 0 )); then
  echo "[ERROR] Fun4All macro failed (rc=$rc)"
  exit $rc
fi
echo "[OK]   Finished successfully → $(ls -l "$out_root" 2>/dev/null || echo '(file not found!)')"
exit 0
