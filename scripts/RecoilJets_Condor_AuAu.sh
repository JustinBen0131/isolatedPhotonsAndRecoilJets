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
    export RJ_DATASET="isSimEmbeddedInclusive"
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
  declare -A fanout_dir_seen=()
  while IFS= read -r fan_line; do
    [[ -z "${fan_line:-}" || "${fan_line:0:1}" == "#" ]] && continue
    IFS='|' read -r -a fan_cols <<< "$fan_line"
    fan_dest="${fan_cols[0]:-}"
    fan_cfg="${fan_cols[1]:-}"
    fan_pre="${fan_cols[2]:-}"
    fan_tight="${fan_cols[3]:-}"
    fan_nonTight="${fan_cols[4]:-}"
    [[ -z "${fan_dest:-}" || "${fan_dest:0:1}" == "#" ]] && continue
    fan_out_dir="${fan_dest}/${run8}"
    if [[ -z "${fanout_dir_seen[$fan_out_dir]:-}" ]]; then
      mkdir -p "$fan_out_dir"
      fanout_dir_seen["$fan_out_dir"]=1
    fi
    fan_out_root="${fan_out_dir}/RecoilJets_${analysis_tag}_${fan_cfg}_${chunk_tag}.root"
    if (( ${#fan_cols[@]} >= 8 )); then
      printf '%s|%s|%s|%s|%s|%s|%s|%s\n' \
        "$fan_out_root" "$fan_cfg" "$fan_pre" "$fan_tight" "$fan_nonTight" \
        "${fan_cols[5]}" "${fan_cols[6]}" "${fan_cols[7]}" >> "$fanout_file"
    else
      printf '%s|%s|%s|%s|%s\n' "$fan_out_root" "$fan_cfg" "$fan_pre" "$fan_tight" "$fan_nonTight" >> "$fanout_file"
    fi
    fanout_outputs+=( "$fan_out_root" )
  done < "$RJ_ID_FANOUT_DIRS_FILE"
  [[ -s "$fanout_file" ]] || { echo "[FATAL] fanout dirs file produced no output rows: $RJ_ID_FANOUT_DIRS_FILE"; exit 7; }
  export RJ_ID_FANOUT_FILE="$fanout_file"
  out_root="${fanout_outputs[0]}"
  echo "[INFO] ID fanout enabled: $(wc -l < "$fanout_file") cfg rows, ${#fanout_dir_seen[@]} output directories from one Fun4All pass"
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
profile_start_epoch="$(date +%s)"
profile_enabled="${RJ_PROFILE_JOB:-0}"
profile_stage="${RJ_PROFILE_STAGE:-analysis}"
profile_label="${RJ_PROFILE_LABEL:-${chunk_tag}}"
profile_file="${TMPDIR:-/tmp}/recoiljets_time_$$_${chunk_tag}.txt"
input_files="$(grep -Ev '^[[:space:]]*($|#)' "$chunk_list" | wc -l | awk '{print $1}')"

file_size_bytes() {
  local f="$1"
  if [[ ! -f "$f" ]]; then
    echo 0
  elif stat -c '%s' "$f" >/dev/null 2>&1; then
    stat -c '%s' "$f"
  else
    wc -c < "$f" | awk '{print $1}'
  fi
}

emit_profile_summary() {
  local exit_code="$1"
  local end_epoch elapsed max_rss_kb user_cpu_s system_cpu_s cpu_percent major_faults minor_faults voluntary_cs involuntary_cs fs_inputs fs_outputs output_files output_bytes f sz fanout_view_count fanout_output_roots
  end_epoch="$(date +%s)"
  elapsed=$(( end_epoch - profile_start_epoch ))
  max_rss_kb="unknown"
  user_cpu_s="unknown"
  system_cpu_s="unknown"
  cpu_percent="unknown"
  major_faults="unknown"
  minor_faults="unknown"
  voluntary_cs="unknown"
  involuntary_cs="unknown"
  fs_inputs="unknown"
  fs_outputs="unknown"
  if [[ -s "$profile_file" ]]; then
    max_rss_kb="$(awk -F: '/Maximum resident set size/ {gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2); print $2; exit}' "$profile_file")"
    [[ -n "$max_rss_kb" ]] || max_rss_kb="unknown"
    user_cpu_s="$(awk -F: '/User time \(seconds\)/ {gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2); print $2; exit}' "$profile_file")"
    system_cpu_s="$(awk -F: '/System time \(seconds\)/ {gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2); print $2; exit}' "$profile_file")"
    cpu_percent="$(awk -F: '/Percent of CPU this job got/ {gsub(/^[[:space:]]+|[[:space:]]+|%$/, "", $2); print $2; exit}' "$profile_file")"
    major_faults="$(awk -F: '/Major \(requiring I\/O\) page faults/ {gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2); print $2; exit}' "$profile_file")"
    minor_faults="$(awk -F: '/Minor \(reclaiming a frame\) page faults/ {gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2); print $2; exit}' "$profile_file")"
    voluntary_cs="$(awk -F: '/Voluntary context switches/ {gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2); print $2; exit}' "$profile_file")"
    involuntary_cs="$(awk -F: '/Involuntary context switches/ {gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2); print $2; exit}' "$profile_file")"
    fs_inputs="$(awk -F: '/File system inputs/ {gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2); print $2; exit}' "$profile_file")"
    fs_outputs="$(awk -F: '/File system outputs/ {gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2); print $2; exit}' "$profile_file")"
    user_cpu_s="${user_cpu_s:-unknown}"
    system_cpu_s="${system_cpu_s:-unknown}"
    cpu_percent="${cpu_percent:-unknown}"
    major_faults="${major_faults:-unknown}"
    minor_faults="${minor_faults:-unknown}"
    voluntary_cs="${voluntary_cs:-unknown}"
    involuntary_cs="${involuntary_cs:-unknown}"
    fs_inputs="${fs_inputs:-unknown}"
    fs_outputs="${fs_outputs:-unknown}"
  fi
  output_files=0
  output_bytes=0
  if (( ${#fanout_outputs[@]} > 0 )); then
    for f in "${fanout_outputs[@]}"; do
      [[ -f "$f" ]] || continue
      sz="$(file_size_bytes "$f")"
      output_bytes=$(( output_bytes + sz ))
      output_files=$(( output_files + 1 ))
    done
  elif [[ -f "$out_root" ]]; then
    output_bytes="$(file_size_bytes "$out_root")"
    output_files=1
  fi
  fanout_view_count=0
  fanout_output_roots=0
  if [[ -n "${RJ_ID_FANOUT_FILE:-}" && -s "$RJ_ID_FANOUT_FILE" ]]; then
    fanout_view_count="$(awk 'NF && $0 !~ /^[[:space:]]*#/ {c++} END{print c+0}' "$RJ_ID_FANOUT_FILE" 2>/dev/null)"
    fanout_output_roots="$(awk -F'|' 'NF && $1 !~ /^#/ && $1 != "" {seen[$1]=1} END{for(k in seen)c++; print c+0}' "$RJ_ID_FANOUT_FILE" 2>/dev/null || echo 0)"
  fi

  echo "RECOILJETS_JOB_PROFILE_V1 stage=${profile_stage} label=${profile_label} dataset=${dataset} analysis_tag=${analysis_tag} run=${run8} chunk=${chunk_tag} input_files=${input_files} nevents=${nevents} cluster_id=${cluster_id} exit_code=${exit_code} elapsed_seconds=${elapsed} max_rss_kb=${max_rss_kb} user_cpu_s=${user_cpu_s} system_cpu_s=${system_cpu_s} cpu_percent=${cpu_percent} major_page_faults=${major_faults} minor_page_faults=${minor_faults} voluntary_context_switches=${voluntary_cs} involuntary_context_switches=${involuntary_cs} fs_inputs=${fs_inputs} fs_outputs=${fs_outputs} output_files=${output_files} output_bytes=${output_bytes} request_memory_mb=${RJ_REQUEST_MEMORY_MB:-unknown} fanout_view_count=${fanout_view_count} fanout_output_roots=${fanout_output_roots} macro=${MACRO} config=${RJ_CONFIG_YAML:-unset}"
  if [[ -s "$profile_file" ]]; then
    sed 's/^/[time-v] /' "$profile_file"
  fi
}

heartbeat_pid=""
start_heartbeat() {
  local hb="${RJ_JOB_HEARTBEAT_SECONDS:-0}"
  [[ "$hb" =~ ^[0-9]+$ && "$hb" -gt 0 ]] || return 0
  (
    while true; do
      sleep "$hb" || exit 0
      local now elapsed out_bytes
      now="$(date +%s)"
      elapsed=$(( now - profile_start_epoch ))
      out_bytes=0
      if [[ -f "$out_root" ]]; then
        out_bytes="$(file_size_bytes "$out_root")"
      fi
      echo "RECOILJETS_JOB_HEARTBEAT_V1 stage=${profile_stage} label=${profile_label} dataset=${dataset} run=${run8} chunk=${chunk_tag} elapsed_seconds=${elapsed} input_files=${input_files} nevents=${nevents} output_bytes=${out_bytes}"
    done
  ) &
  heartbeat_pid="$!"
}

stop_heartbeat() {
  if [[ -n "${heartbeat_pid:-}" ]]; then
    kill "$heartbeat_pid" >/dev/null 2>&1 || true
    wait "$heartbeat_pid" >/dev/null 2>&1 || true
    heartbeat_pid=""
  fi
}

set +e
echo "[INFO] Running ROOT:"
echo "root -b -q -l \"${MACRO}(${nevents}, \\\"${chunk_list}\\\", \\\"${out_root}\\\", false)\""
start_heartbeat
if [[ "$profile_enabled" == "1" || "$profile_enabled" == "true" || "$profile_enabled" == "TRUE" ]] && command -v /usr/bin/time >/dev/null 2>&1; then
  /usr/bin/time -v -o "$profile_file" root -b -q -l "${MACRO}(${nevents}, \"${chunk_list}\", \"${out_root}\", false)"
else
  root -b -q -l "${MACRO}(${nevents}, \"${chunk_list}\", \"${out_root}\", false)"
fi
rc=$?
stop_heartbeat
set -e

echo "---------------------------------------------------------------------"
emit_profile_summary "$rc"
rm -f "$profile_file"
if (( rc != 0 )); then
  if [[ "${RJ_ALLOW_NONZERO_WITH_ROOT_OUTPUT:-0}" == "1" && -s "$out_root" ]]; then
    echo "[WARN] Fun4All macro returned rc=$rc, but RJ_ALLOW_NONZERO_WITH_ROOT_OUTPUT=1 and primary ROOT output exists."
    echo "[WARN] Treating this worker as successful; downstream validation must still open and verify the ROOT file."
    rc=0
  fi
fi
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
