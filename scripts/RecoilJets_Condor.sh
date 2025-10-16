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
run8="${1:?run8 required (8-digit)}"
chunk_list="${2:?chunk list (.list) required}"
dataset_raw="${3:-isAuAu}"          # isPP | isAuAu
cluster_id="${4:-LOCAL}"            # informational
nevents="${5:-0}"                   # 0 → all events in the list
chunk_idx="${6:-0}"                 # informational; not used for naming
_ignored2="${7:-}"                  # keep slot for compatibility (NONE)
dest_base="${8:-}"                  # If empty, we derive from dataset

# ------------------------ Fixed paths ----------------------
BASE="/sphenix/u/patsfan753/scratch/thesisAnalysis"
MACRO="${BASE}/Fun4All_recoilJets.C"

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
MYINSTALL="/sphenix/user/${USER}/install"
source /opt/sphenix/core/bin/sphenix_setup.sh -n
if [[ -d "$MYINSTALL" ]]; then
  # do not fail if local area is not present; macro has R__LOAD_LIBRARY with absolute path
  source /opt/sphenix/core/bin/setup_local.sh "$MYINSTALL" || true
fi

# ------------------------ Dataset routing ------------------
# Normalize dataset and set defaults:
case "$dataset_raw" in
  isPP|pp|PP)     dataset="isPP" ;;
  isAuAu|auau|AA) dataset="isAuAu" ;;
  *)              echo "[WARN] Unknown dataset '$dataset_raw' → defaulting to 'isAuAu'"; dataset="isAuAu" ;;
esac
export RJ_DATASET="$dataset"   # Fun4All macro will read this and call setDataType()

# Destination base (if not supplied as arg 8)
if [[ -z "$dest_base" ]]; then
  if [[ "$dataset" == "isPP" ]]; then
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
out_root="${out_dir}/RecoilJets_${dataset}_${chunk_tag}.root"

echo "[INFO] Output path = $out_root"

# The macro expects to find the ROOT files listed inside chunk_list. If those
# entries are relative file names, make sure the CWD is the directory where
# the ROOT files live (we use the list’s directory as a safe default).
list_dir="$(dirname "$chunk_list")"
cd "$list_dir"

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
