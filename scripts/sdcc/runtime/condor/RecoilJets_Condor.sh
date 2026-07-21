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

validate_archived_ppg12_di_chunk_list() {
  local list="$1"
  local sample_slice="$2"
  awk -F '\t' \
    -v g4="/js_pp200_signal_dual/g4hits/run0028/${sample_slice}/" \
    -v jets="/js_pp200_signal_dual/nopileup/jets/run0028/${sample_slice}/" '
    BEGIN { bad=0 }
    NF != 5 || $1 != "NONE" || index($2, g4) == 0 ||
    index($3, jets) == 0 || $4 != "NONE" || $5 != "NONE" {
      if (bad < 5) {
        printf "Archived PPG12 DI row %d mismatch: CALO=%s G4Hits=%s DST_JETS=%s DST_GLOBAL=%s MBD=%s\n", NR, $1, $2, $3, $4, $5 > "/dev/stderr"
      }
      bad++
    }
    END { exit bad == 0 ? 0 : 1 }
  ' "$list"
}

# ------------------------ Fixed paths ----------------------
BASE="/sphenix/u/patsfan753/scratch/thesisAnalysis"
MACRO="${RJ_MACRO_PATH:-${BASE}/macros/Fun4All_recoilJets.C}"

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
PPG12_ARCHIVED_OFFLINE_MAIN="/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_ana/ana.541"
ppg12_archived_di_lane=0
case "$run8" in
  run28_jet8_double|run28_jet12_double|run28_jet20_double|run28_jet30_double|run28_jet40_double)
    [[ "$dataset_raw" == "isSimInclusive" ]] && ppg12_archived_di_lane=1
    ;;
  run28_jet5_double)
    if [[ "$dataset_raw" == "isSimInclusive" ]]; then
      echo "[FATAL] run28_jet5_double is outside the canonical archived-DI inclusive replacement contract."
      exit 98
    fi
    ;;
esac

if (( ppg12_archived_di_lane )); then
  archived_runtime_manifest="${RJ_PPG12_DI_RUNTIME_MANIFEST:-}"
  archived_runtime_manifest_sha256="${RJ_PPG12_DI_RUNTIME_MANIFEST_SHA256:-}"
  archived_period="${RJ_PPG12_PERIOD:-}"
  [[ "$archived_period" =~ ^(0mrad|1p5mrad)$ ]] || {
    echo "[FATAL] Archived PPG12 DI runtime requires RJ_PPG12_PERIOD=0mrad or 1p5mrad; got '${archived_period:-<unset>}'."
    exit 98
  }

  # Condor uses getenv=True. Remove inherited values that could widen or
  # redirect this branch, then pin the accepted archived production contract.
  unset RJ_PPG12_DI_ARCHIVED_RECO_CHAIN RJ_TRUTH_JETS_MODE
  unset RJ_PPG12_DI_ARCHIVED_RELEASE RJ_PPG12_DI_RUNTIME_MANIFEST
  unset RJ_PPG12_DI_RUNTIME_MANIFEST_SHA256
  unset RJ_PPG12_DI_ARCHIVED_EXPECT_PEDESTAL RJ_PPG12_PEDESTAL_OVERRIDE
  unset RJ_PPG12_PERIOD_ALLOW_ALL_SIM RJ_PPG12_PHOTON_YIELD_DOUBLE
  unset RJ_PPG12_PERIOD_ALLOW_MIX_OVERRIDE RJ_PPG12_PERIOD_ALLOW_VERTEX_FILE_OVERRIDE
  unset RJ_PPG12_PERIOD_USE_LUMI_WEIGHT RJ_PPG12_PHOTON_YIELD_MIX_WEIGHT
  unset RJ_PPG12_CROSSING_PERIOD RJ_PP_VERTEX_REWEIGHT_FILE RJ_PP_VERTEX_REWEIGHT_HIST
  unset RJ_PPG12_PERIOD_STRICT_DI RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4
  unset RJ_PPG12_PPSIM_G4_ONLY RJ_SIM_ALLOW_NONE_LISTS
  unset RJ_FORCE_RELEASE_CORE_LIBS RJ_RELEASE_CORE_LIB_DIR
  unset RJ_RELEASE_CORE_LIB64_DIR RJ_FORCE_RELEASE_CALO_IO
  unset RJ_RELEASE_CALO_IO_PATH RJ_SNAPSHOT_LIB_DIR
  unset RJ_SIM_SAMPLE RJ_EMBEDDED_INCLUSIVE_JET_SAMPLE
  unset LD_PRELOAD LD_LIBRARY_PATH ROOT_INCLUDE_PATH
  export RJ_PPG12_DI_ARCHIVED_RECO_CHAIN=1
  export RJ_TRUTH_JETS_MODE=DST
  export RJ_PPG12_DI_ARCHIVED_RELEASE="$PPG12_ARCHIVED_OFFLINE_MAIN"
  export RJ_PPG12_PERIOD_ALLOW_ALL_SIM=0
  export RJ_PPG12_PERIOD_ALLOW_MIX_OVERRIDE=0
  export RJ_PPG12_PERIOD_ALLOW_VERTEX_FILE_OVERRIDE=0
  export RJ_PPG12_PERIOD_USE_LUMI_WEIGHT=1
  export RJ_PPG12_PERIOD="$archived_period"
  export RJ_PPG12_PHOTON_YIELD_DOUBLE=1
  export RJ_PPG12_PERIOD_STRICT_DI=1
  export RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4=1
  export RJ_PPG12_PPSIM_G4_ONLY=1
  export RJ_SIM_ALLOW_NONE_LISTS=1
  export RJ_SIM_SAMPLE="$run8"
  export RJ_EMBEDDED_INCLUSIVE_JET_SAMPLE="$run8"
  export RJ_PPG12_DI_RUNTIME_MANIFEST="$archived_runtime_manifest"
  export RJ_PPG12_DI_RUNTIME_MANIFEST_SHA256="$archived_runtime_manifest_sha256"
  archived_wrapper_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
  MACRO="${archived_wrapper_dir}/Fun4All_recoilJets.C"
  export RJ_MACRO_PATH="$MACRO"

  [[ "$archived_runtime_manifest" == /* && -s "$archived_runtime_manifest" ]] || {
    echo "[FATAL] Archived PPG12 DI runtime manifest is missing or not absolute: ${archived_runtime_manifest:-<unset>}"
    exit 96
  }
  [[ "$archived_runtime_manifest_sha256" =~ ^[0-9a-f]{64}$ ]] || {
    echo "[FATAL] Archived PPG12 DI runtime-manifest digest is invalid."
    exit 96
  }
  command -v sha256sum >/dev/null 2>&1 || {
    echo "[FATAL] sha256sum is required to verify the archived PPG12 DI runtime."
    exit 96
  }
  actual_runtime_manifest_sha256="$(sha256sum "$archived_runtime_manifest" | awk '{print $1}')"
  [[ "$actual_runtime_manifest_sha256" == "$archived_runtime_manifest_sha256" ]] || {
    echo "[FATAL] Archived PPG12 DI runtime-manifest digest mismatch."
    exit 96
  }
  (
    cd "$(dirname "$archived_runtime_manifest")"
    sha256sum -c "$(basename "$archived_runtime_manifest")"
  ) || {
    echo "[FATAL] Archived PPG12 DI runtime snapshot failed hash verification."
    exit 96
  }

  # Source only the fixed historical release in this branch.
  set +u
  source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.541
  set -u
  [[ "${OFFLINE_MAIN:-}" == "$PPG12_ARCHIVED_OFFLINE_MAIN" ]] || {
    echo "[FATAL] Archived PPG12 DI runtime resolved OFFLINE_MAIN='${OFFLINE_MAIN:-<unset>}', expected '${PPG12_ARCHIVED_OFFLINE_MAIN}'."
    exit 96
  }
else
  # Existing runtime for every non-archived lane remains unchanged.
  set +u
  source /opt/sphenix/core/bin/sphenix_setup.sh -n
  if [[ -d "$MYINSTALL" ]]; then
    source /opt/sphenix/core/bin/setup_local.sh "$MYINSTALL" || true
  fi
  if [[ -d "${MYINSTALL}/lib" ]]; then
    export LD_LIBRARY_PATH="${MYINSTALL}/lib:${LD_LIBRARY_PATH:-}"
  fi
  if [[ -d "${MYINSTALL}/include" ]]; then
    export ROOT_INCLUDE_PATH="${MYINSTALL}/include:${ROOT_INCLUDE_PATH:-}"
  fi
  set -u
fi

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
    export RJ_PPG12_PP_DATA_PAIRED="${RJ_PPG12_PP_DATA_PAIRED:-1}"
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
  isSimInclusive|issiminclusive|siminclusive|SIMINCLUSIVE)
    dataset="isSimInclusive"
    analysis_tag="isSimInclusive"
    export RJ_DATASET="isSimInclusive"
    export RJ_IS_SIM=1
    ;;
  isSimJet5|simjet5|SIMJET5)
    dataset="isSimJet5"
    analysis_tag="isSimInclusive"
    export RJ_DATASET="isSimInclusive"
    export RJ_IS_SIM=1
    ;;
  isSimMB|simmb|SIMMB)
    dataset="isSimMB"
    analysis_tag="isSimMB"
    export RJ_DATASET="isSimMB"
    export RJ_IS_SIM=1
    ;;
  isAuAu|auau|AA)
    echo "[FATAL] RecoilJets_Condor.sh is pp-style only. Use RecoilJets_Condor_AuAu.sh for isAuAu."
    exit 50
    ;;
  isSimEmbedded|issimembedded|simembedded|SIMEMBEDDED)
    echo "[FATAL] RecoilJets_Condor.sh is pp-style only. Use RecoilJets_Condor_AuAu.sh for isSimEmbedded."
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
export RJ_SIM_SAMPLE="$run8"

if [[ "$dataset" == "isSim" ]]; then
  # For pp photon+jet production run8 is the slice label passed by the submitter
  # (PhotonJet5, PhotonJet10, PhotonJet20). RecoilJets uses this to apply the
  # same PPG12 stitching window before any weighted SIM histogram filling.
  export RJ_PPG12_PHOTON_SAMPLE="$run8"
fi

# Destination base (if not supplied as arg 8)
if [[ -z "$dest_base" ]]; then
  if [[ "$analysis_tag" == "isSimInclusive" ]]; then
    dest_base="/sphenix/tg/tg01/bulk/jbennett/thesisAna/siminclusive"
  elif [[ "$analysis_tag" == "isSimMB" ]]; then
    dest_base="/sphenix/tg/tg01/bulk/jbennett/thesisAna/simmb"
  elif [[ "$analysis_tag" == "isSim" ]]; then
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
if (( ppg12_archived_di_lane )); then
  sample_slice="${run8#run28_}"
  sample_slice="${sample_slice%_double}"
  if ! validate_archived_ppg12_di_chunk_list "$chunk_list" "$sample_slice"; then
    echo "[FATAL] Archived PPG12 DI worker requires the exact G4-only graph: CALO/GLOBAL/MBD=NONE plus run-28 ${sample_slice} dual G4Hits and truth-jet sources."
    exit 98
  fi
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

require_non_tiny_output() {
  case "${RJ_REQUIRE_NON_TINY_OUTPUT:-0}" in
    1|true|TRUE|yes|YES|on|ON) return 0 ;;
    *) return 1 ;;
  esac
}

min_output_bytes() {
  local min="${RJ_MIN_OUTPUT_BYTES:-50000}"
  if [[ ! "$min" =~ ^[0-9]+$ ]]; then
    echo "[WARN] Invalid RJ_MIN_OUTPUT_BYTES='$min'; using 50000" >&2
    min=50000
  fi
  # Sparse pp data chunks can be valid around 30-45 KB after event/candidate
  # cuts. Keep the stub guard, but do not hold those successful chunks.
  if [[ "${dataset:-${RJ_DATASET:-}}" == "isPP" && "$min" -gt 25000 ]]; then
    min=25000
  fi
  echo "$min"
}

root_output_structurally_valid() {
  local f="$1"
  [[ -f "$f" ]] || return 1
  python3 - "$f" <<'PY'
import sys

path = sys.argv[1]
try:
    import ROOT
except Exception as exc:
    print(f"[WARN] Could not import ROOT for small-output validation: {exc}", file=sys.stderr)
    sys.exit(2)

ROOT.gROOT.SetBatch(True)
tf = ROOT.TFile.Open(path)
if not tf or tf.IsZombie():
    sys.exit(1)

keys = list(tf.GetListOfKeys())
has_config = any(key.GetName() == "analysis_config_yaml" for key in keys)
has_directory = any(key.ReadObj().InheritsFrom("TDirectory") for key in keys)
tf.Close()
sys.exit(0 if (has_config and has_directory) else 1)
PY
}

rj_truthy() {
  case "${1:-0}" in
    1|true|TRUE|yes|YES|on|ON) return 0 ;;
    *) return 1 ;;
  esac
}

root_invoke_prefix=()
if rj_truthy "${RJ_FORCE_RELEASE_CALO_IO:-0}"; then
  release_calo_io="${RJ_RELEASE_CALO_IO_PATH:-/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_ana/ana.558/lib/libcalo_io.so.0}"
  if [[ ! -r "$release_calo_io" ]]; then
    echo "[FATAL] RJ_FORCE_RELEASE_CALO_IO requested, but libcalo_io is not readable: $release_calo_io"
    exit 9
  fi
  root_invoke_prefix=( env "LD_PRELOAD=${release_calo_io}${LD_PRELOAD:+:${LD_PRELOAD}}" )
  echo "[INFO] RJ_FORCE_RELEASE_CALO_IO=1: preloading $release_calo_io for ROOT invocation only"
fi

check_required_output_file() {
  local f="$1"
  local label="$2"
  local min bytes
  min="$(min_output_bytes)"
  if [[ ! -f "$f" ]]; then
    echo "[ERROR] Missing required ${label} output: $f"
    return 1
  fi
  bytes="$(file_size_bytes "$f")"
  if (( bytes < min )); then
    if [[ "${dataset:-${RJ_DATASET:-}}" == "isAuAu" ]] && root_output_structurally_valid "$f"; then
      echo "[INFO] Required ${label} output passed small AuAu ROOT structural validation: $f (${bytes} bytes < ${min})"
      return 0
    fi
    echo "[ERROR] Required ${label} output is too small: $f (${bytes} bytes < ${min})"
    return 1
  fi
  return 0
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

  echo "RECOILJETS_JOB_PROFILE_V1 stage=${profile_stage} label=${profile_label} dataset=${dataset} analysis_tag=${analysis_tag} run=${run8} chunk=${chunk_tag} input_files=${input_files} nevents=${nevents} cluster_id=${cluster_id} exit_code=${exit_code} elapsed_seconds=${elapsed} max_rss_kb=${max_rss_kb} user_cpu_s=${user_cpu_s} system_cpu_s=${system_cpu_s} cpu_percent=${cpu_percent} major_page_faults=${major_faults} minor_page_faults=${minor_faults} voluntary_context_switches=${voluntary_cs} involuntary_context_switches=${involuntary_cs} fs_inputs=${fs_inputs} fs_outputs=${fs_outputs} output_files=${output_files} output_bytes=${output_bytes} request_memory_mb=${RJ_REQUEST_MEMORY_MB:-unknown} fanout_view_count=${fanout_view_count} fanout_output_roots=${fanout_output_roots} archived_chain=${RJ_PPG12_DI_ARCHIVED_RECO_CHAIN:-0} truth_jets_mode=${RJ_TRUTH_JETS_MODE:-unset} release=${RJ_PPG12_DI_ARCHIVED_RELEASE:-${OFFLINE_MAIN:-unset}} sample=${RJ_SIM_SAMPLE:-$run8} runtime_manifest_sha256=${RJ_PPG12_DI_RUNTIME_MANIFEST_SHA256:-unset} macro=${MACRO} config=${RJ_CONFIG_YAML:-unset}"
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

# Defensive default: with set -u enabled, never let wrapper bookkeeping turn
# into an unbound rc hold. ROOT success/failure overwrites this immediately.
rc=125
set +e
echo "[INFO] Running ROOT:"
echo "root -b -q -l \"${MACRO}(${nevents}, \\\"${chunk_list}\\\", \\\"${out_root}\\\", false)\""
start_heartbeat
if [[ "$profile_enabled" == "1" || "$profile_enabled" == "true" || "$profile_enabled" == "TRUE" ]] && command -v /usr/bin/time >/dev/null 2>&1; then
  /usr/bin/time -v -o "$profile_file" "${root_invoke_prefix[@]}" root -b -q -l "${MACRO}(${nevents}, \"${chunk_list}\", \"${out_root}\", false)"
else
  "${root_invoke_prefix[@]}" root -b -q -l "${MACRO}(${nevents}, \"${chunk_list}\", \"${out_root}\", false)"
fi
rc=$?
stop_heartbeat
set -e

echo "---------------------------------------------------------------------"
emit_profile_summary "$rc"
rm -f "$profile_file"
if (( rc != 0 )); then
  echo "[ERROR] Fun4All macro failed (rc=$rc)"
  exit $rc
fi
if (( ${#fanout_outputs[@]} > 0 )); then
  missing=0
  primary_fanout_output="${fanout_outputs[0]}"
  allow_missing_secondary_fanout=0
  if rj_truthy "${RJ_ALLOW_MISSING_SECONDARY_FANOUT:-0}"; then
    allow_missing_secondary_fanout=1
  fi
  for f in "${fanout_outputs[@]}"; do
    if [[ ! -s "$f" ]]; then
      if (( allow_missing_secondary_fanout )) && [[ "$f" != "$primary_fanout_output" ]]; then
        echo "[WARN] Missing secondary fanout output allowed by RJ_ALLOW_MISSING_SECONDARY_FANOUT=1: $f"
        continue
      fi
      echo "[ERROR] Missing fanout output: $f"
      missing=1
    elif require_non_tiny_output && ! check_required_output_file "$f" "fanout"; then
      missing=1
    fi
  done
  (( missing == 0 )) || exit 8
  echo "[OK]   Finished successfully → ${#fanout_outputs[@]} fanout ROOT files"
else
  if require_non_tiny_output; then
    check_required_output_file "$out_root" "primary" || exit 8
  elif [[ ! -f "$out_root" ]]; then
    echo "[WARN] ROOT exited with rc=0 but primary output is missing and RJ_REQUIRE_NON_TINY_OUTPUT=0: $out_root"
  fi
  echo "[OK]   Finished successfully → $(ls -l "$out_root" 2>/dev/null || echo '(file not found!)')"
fi
exit 0
