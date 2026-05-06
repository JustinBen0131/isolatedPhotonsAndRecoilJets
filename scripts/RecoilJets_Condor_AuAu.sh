#!/usr/bin/env bash
###############################################################################
# RecoilJets_Condor_AuAu.sh  —  Condor executable (per-job wrapper)
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
dataset_raw="${3:-isAuAu}"          # isPP | isAuAu | isOO | isSim
cluster_id="${4:-LOCAL}"            # informational
nevents="${5:-0}"                   # 0 → all events in the list
chunk_idx="${6:-0}"                 # informational; not used for naming
_ignored2="${7:-}"                  # keep slot for compatibility (NONE)
dest_base="${8:-}"                  # If empty, we derive from dataset

# ------------------------ Fixed paths ----------------------
BASE="/sphenix/u/patsfan753/scratch/thesisAnalysis"
MACRO="${RJ_MACRO_PATH:-${BASE}/macros/Fun4All_recoilJets_AuAu.C}"

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
MYINSTALL_AUAU="/sphenix/u/${USER}/thesisAnalysis_auau/install"

# Disable 'nounset' while sourcing env scripts; they may read unset vars (e.g. PGHOST)
set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n
if [[ -d "$MYINSTALL" ]]; then
  # do not fail if local area is not present; macro has R__LOAD_LIBRARY with absolute path
  source /opt/sphenix/core/bin/setup_local.sh "$MYINSTALL" || true
fi
if [[ -d "$MYINSTALL_AUAU" ]]; then
  # AuAu local area (for libRecoilJetsAuAu.so and any future AuAu-only libs)
  source /opt/sphenix/core/bin/setup_local.sh "$MYINSTALL_AUAU" || true
fi
set -u

# ------------------------ Dataset routing ------------------
# Normalize dataset and set defaults:
#  - isSim must remain isSim end-to-end so the analysis module can detect it.
#  - Fun4All macro will treat isSim as pp-style reconstruction internally.
analysis_tag="isAuAu"
case "$dataset_raw" in
  isPP|pp|PP)
    dataset="isPP"
    analysis_tag="isPP"
    export RJ_DATASET="isPP"
    export RJ_IS_SIM=0
    ;;
  isAuAu|auau|AA)
    dataset="isAuAu"
    analysis_tag="isAuAu"
    export RJ_DATASET="isAuAu"
    export RJ_IS_SIM=0
    ;;
  isOO|oo|OO)
    dataset="isOO"
    analysis_tag="isOO"
    export RJ_DATASET="isAuAu"
    export RJ_IS_SIM=0
    ;;
  isSim|sim|SIM)
    dataset="isSim"
    analysis_tag="isSim"
    export RJ_DATASET="isSim"
    export RJ_IS_SIM=1
    ;;
  isSimEmbedded|issimembedded|simembedded|SIMEMBEDDED)
    dataset="isSimEmbedded"
    analysis_tag="isSimEmbedded"
    export RJ_DATASET="isSimEmbedded"
    export RJ_IS_SIM=1
    ;;
  isSimEmbeddedInclusive|issimembeddedinclusive|simembeddedinclusive|SIMEMBEDDEDINCLUSIVE)
    dataset="isSimEmbeddedInclusive"
    analysis_tag="isSimEmbeddedInclusive"
    export RJ_DATASET="isSimEmbedded"
    export RJ_IS_SIM=1
    ;;
  *)
    echo "[WARN] Unknown dataset '$dataset_raw' → defaulting to 'isAuAu'"
    dataset="isAuAu"
    analysis_tag="isAuAu"
    export RJ_DATASET="isAuAu"
    export RJ_IS_SIM=0
    ;;
esac
export RJ_SIM_SAMPLE="$run8"

# Destination base (if not supplied as arg 8)
if [[ -z "$dest_base" ]]; then
  if [[ "$analysis_tag" == "isSimEmbedded" || "$analysis_tag" == "isSimEmbeddedInclusive" ]]; then
    dest_base="/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded"
  elif [[ "$analysis_tag" == "isSim" ]]; then
    dest_base="/sphenix/tg/tg01/bulk/jbennett/thesisAna/sim"
  elif [[ "$analysis_tag" == "isOO" ]]; then
    dest_base="/sphenix/tg/tg01/bulk/jbennett/thesisAna/oo"
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
echo "[INFO] Debug env  : RJ_VERBOSITY=${RJ_VERBOSITY:-unset}  RJ_F4A_VERBOSE=${RJ_F4A_VERBOSE:-unset}  RJ_STEP_EVENTS=${RJ_STEP_EVENTS:-unset}  RJ_CRASH_BACKTRACE=${RJ_CRASH_BACKTRACE:-unset}"

fanout_outputs=()
if [[ -n "${RJ_ID_FANOUT_DIRS_FILE:-}" ]]; then
  [[ -f "$RJ_ID_FANOUT_DIRS_FILE" ]] || { echo "[FATAL] RJ_ID_FANOUT_DIRS_FILE not found: $RJ_ID_FANOUT_DIRS_FILE"; exit 6; }
  fanout_file="${TMPDIR:-/tmp}/rj_id_fanout_$$_${chunk_tag}.txt"
  : > "$fanout_file"
  while IFS= read -r fan_line; do
    [[ -z "${fan_line:-}" || "${fan_line:0:1}" == "#" ]] && continue
    IFS='|' read -r -a fan_cols <<< "$fan_line"
    fan_dest="${fan_cols[0]:-}"
    fan_cfg="${fan_cols[1]:-}"
    fan_yaml=""
    if (( ${#fan_cols[@]} >= 6 )); then
      fan_yaml="${fan_cols[2]:-}"
      fan_pre="${fan_cols[3]:-}"
      fan_tight="${fan_cols[4]:-}"
      fan_nonTight="${fan_cols[5]:-}"
      fan_view="${fan_cols[6]:-}"
      fan_materialize="${fan_cols[7]:-}"
    else
      fan_pre="${fan_cols[2]:-}"
      fan_tight="${fan_cols[3]:-}"
      fan_nonTight="${fan_cols[4]:-}"
      fan_view=""
      fan_materialize=""
    fi
    [[ -z "${fan_dest:-}" || "${fan_dest:0:1}" == "#" ]] && continue
    fan_out_dir="${fan_dest}/${run8}"
    mkdir -p "$fan_out_dir"
    fan_out_root="${fan_out_dir}/RecoilJets_${analysis_tag}_${chunk_tag}.root"
    if [[ -n "$fan_yaml" ]]; then
      printf '%s|%s|%s|%s|%s|%s|%s|%s\n' "$fan_out_root" "$fan_cfg" "$fan_yaml" "$fan_pre" "$fan_tight" "$fan_nonTight" "$fan_view" "$fan_materialize" >> "$fanout_file"
    else
      printf '%s|%s|%s|%s|%s\n' "$fan_out_root" "$fan_cfg" "$fan_pre" "$fan_tight" "$fan_nonTight" >> "$fanout_file"
    fi
    fanout_outputs+=( "$fan_out_root" )
  done < "$RJ_ID_FANOUT_DIRS_FILE"
  [[ -s "$fanout_file" ]] || { echo "[FATAL] fanout dirs file produced no output rows: $RJ_ID_FANOUT_DIRS_FILE"; exit 7; }
  export RJ_ID_FANOUT_FILE="$fanout_file"
  export RJ_REPLAY_FANOUT_FILE="$fanout_file"
  out_root="${fanout_outputs[0]}"
  echo "[INFO] ID fanout enabled: $(wc -l < "$fanout_file") outputs from one Fun4All pass"
  echo "[INFO] ID fanout file   : $fanout_file"
  echo "[INFO] Primary output   : $out_root"
fi

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
if (( ${#fanout_outputs[@]} > 0 )); then
  missing=0
  for f in "${fanout_outputs[@]}"; do
    if [[ ! -s "$f" ]]; then
      echo "[ERROR] Missing fanout output: $f"
      missing=1
    fi
  done
  (( missing == 0 )) || exit 8
  echo "[OK]   Finished successfully → ${#fanout_outputs[@]} fanout ROOT files"
else
  echo "[OK]   Finished successfully → $(ls -l "$out_root" 2>/dev/null || echo '(file not found!)')"
fi
exit 0
