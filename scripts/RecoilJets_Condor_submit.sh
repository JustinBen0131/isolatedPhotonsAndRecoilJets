#!/usr/bin/env bash
###############################################################################
# RecoilJets_Condor_submit.sh  —  one-stop driver for local & Condor running
#
# WHAT THIS SCRIPT DOES (DATA ONLY; NO SIM):
#   • Supports two datasets: isPP and isAuAu
#   • Reads golden run lists and per-run .list files you provided
#   • Splits each run’s .list into fixed-size groups (no cross-run mixing)
#   • Builds Condor submit files and queues one job per group
#   • Optionally: create “rounds” (segments) of runs with a job cap
#   • Local test mode: run first file of the first run for N events
#
# DIRECTORY ASSUMPTIONS (from your message):
#   Base:  /sphenix/u/patsfan753/scratch/thesisAnalysis
#     • Fun4All macro:     ${BASE}/Fun4All_recoilJets.C
#     • Golden lists:
#         - PP:   ${BASE}/Full_ppGoldenRunList_Version3.list
#         - AuAu: ${BASE}/Full_AuAuGoldenRunList_Version1_personal.list
#     • Per-run input lists:
#         - PP:   ${BASE}/dst_lists_pp/dst_calofitting-<run8>.list
#         - AuAu: ${BASE}/dst_lists_auau/dst_calofitting-<run8>.list
#     • Condor staging (this script creates):
#         - ${BASE}/condor_lists/<pp|auau>/*.list      (grouped lists)
#         - ${BASE}/condor_segments/<pp|auau>/*.txt    (rounds of runs)
#     • Condor logs (requested):
#         - ${BASE}/stdout, ${BASE}/error, ${BASE}/log
#
# OUTPUT DESTINATION:
#   • PP  → /sphenix/tg/tg01/bulk/jbennett/thesisAna/pp/<run8>/
#   • AuAu→ /sphenix/tg/tg01/bulk/jbennett/thesisAna/auau/<run8>/
#
# USAGE EXAMPLES:
#   # Local quick check on PP (first run, first file), 5000 events:
#   ./RecoilJets_Condor_submit.sh isPP  local 5000
#
#   # Create round files (segments) for AuAu with group size 3 and 15000 job cap:
#   ./RecoilJets_Condor_submit.sh isAuAu splitGoldenRunList groupSize 3 maxJobs 15000
#
#   # Submit AuAu “round 1” with the defaults (groupSize=3, maxJobs=15000):
#   ./RecoilJets_Condor_submit.sh isAuAu condor round 1
#
#   # Submit PP “round 2” but only the first group job of each run (smoke test):
#   ./RecoilJets_Condor_submit.sh isPP condor round 2 firstChunk
#
#   # Submit PP without pre-splitting (all runs currently in the golden list):
#   ./RecoilJets_Condor_submit.sh isPP condor all
###############################################################################
set -euo pipefail

# ------------------------ Pretty printing ------------------
BOLD=$'\e[1m'; DIM=$'\e[2m'; RED=$'\e[31m'; YEL=$'\e[33m'; GRN=$'\e[32m'; BLU=$'\e[34m'; RST=$'\e[0m'
say()  { printf "${BLU}➜${RST} %s\n" "$*"; }
warn() { printf "${YEL}⚠ %s${RST}\n" "$*" >&2; }
err()  { printf "${RED}✘ %s${RST}\n" "$*" >&2; }

# ------------------------ Fixed paths ----------------------
BASE="/sphenix/u/patsfan753/scratch/thesisAnalysis"
MACRO="${BASE}/Fun4All_recoilJets.C"
EXE="${BASE}/RecoilJets_Condor.sh"

# Logs (as requested)
LOG_DIR="${BASE}/log"
OUT_DIR="${BASE}/stdout"
ERR_DIR="${BASE}/error"

# Golden lists provided by you
PP_GOLDEN="${BASE}/Full_ppGoldenRunList_Version3.list"
AA_GOLDEN="${BASE}/Full_AuAuGoldenRunList_Version1_personal.list"

# Per-run input list directories
PP_LIST_DIR="${BASE}/dst_lists_pp"
AA_LIST_DIR="${BASE}/dst_lists_auau"

# Where we stage grouped chunk .list files and the round files
STAGE_ROOT="${BASE}/condor_lists"
ROUND_ROOT="${BASE}/condor_segments"

# Dataset-specific output trees on TG storage
PP_DEST_BASE="/sphenix/tg/tg01/bulk/jbennett/thesisAna/pp"
AA_DEST_BASE="/sphenix/tg/tg01/bulk/jbennett/thesisAna/auau"

# ------------------------ Defaults -------------------------
GROUP_SIZE=3         # files per Condor job (never mixes runs)
MAX_JOBS=15000      # job budget per round file
LOCAL_EVENTS=5000   # default N for "local" if not given

# ------------------------ Helpers --------------------------
usage() {
  cat <<USAGE
${BOLD}Usage:${RST}
  ${BOLD}$0 <isPP|isAuAu> local [Nevents]${RST}
  ${BOLD}$0 <isPP|isAuAu> splitGoldenRunList [groupSize N] [maxJobs M]${RST}
  ${BOLD}$0 <isPP|isAuAu> condor round <N> [groupSize N] [maxJobs M] [firstChunk]${RST}
  ${BOLD}$0 <isPP|isAuAu> condor all [groupSize N]${RST}

Examples:
  $0 isPP local 5000
  $0 isAuAu splitGoldenRunList groupSize 3 maxJobs 12000
  $0 isAuAu condor round 1
  $0 isPP  condor round 2 firstChunk
  $0 isPP  condor all groupSize 4
USAGE
  exit 2
}

need_cmd() { command -v "$1" >/dev/null 2>&1 || { err "Missing required command: $1"; exit 3; }; }

# Dataset resolver (sets globals: DATASET, GOLDEN, LIST_DIR, DEST_BASE, TAG, STAGE_DIR, ROUND_DIR)
resolve_dataset() {
  case "${1:-}" in
    isPP|pp|PP)
      DATASET="isPP"
      GOLDEN="$PP_GOLDEN"
      LIST_DIR="$PP_LIST_DIR"
      DEST_BASE="$PP_DEST_BASE"
      TAG="pp"
      ;;
    isAuAu|auau|AA)
      DATASET="isAuAu"
      GOLDEN="$AA_GOLDEN"
      LIST_DIR="$AA_LIST_DIR"
      DEST_BASE="$AA_DEST_BASE"
      TAG="auau"
      ;;
    *)
      usage
      ;;
  esac
  STAGE_DIR="${STAGE_ROOT}/${TAG}"
  ROUND_DIR="${ROUND_ROOT}/${TAG}"
  mkdir -p "$STAGE_DIR" "$ROUND_DIR" "$LOG_DIR" "$OUT_DIR" "$ERR_DIR"
}

# Count ceil(n/d)
ceil_div() { local n="$1" d="$2"; echo $(( (n + d - 1) / d )); }

# Format run to 8 digits
run8() { printf "%08d" "$((10#$1))"; }

# Build per-run grouped list files; prints absolute paths to grouped lists, one per line
#   make_groups <run8> <groupSize>  → writes STAGE_DIR/run<run8>_grpXXX.list files
make_groups() {
  local r8="$1" gs="$2"
  local src="${LIST_DIR}/dst_calofitting-${r8}.list"
  [[ -s "$src" ]] || { warn "List is missing/empty for run ${r8} → $src"; return 1; }

  # Clean old groups for this run
  rm -f "${STAGE_DIR}/run${r8}_grp"*.list 2>/dev/null || true

  local nfiles; nfiles=$(wc -l < "$src" | awk '{print $1}')
  local ngroups; ngroups=$(ceil_div "$nfiles" "$gs")

  local start=1
  local g=1
  while (( g <= ngroups )); do
    local end=$(( start + gs - 1 ))
    local out="${STAGE_DIR}/run${r8}_grp$(printf "%03d" "$g").list"
    # sed is 1-indexed on lines; clamp 'end' to nfiles
    if (( end > nfiles )); then end="$nfiles"; fi
    sed -n "${start},${end}p" "$src" > "$out"
    echo "$out"
    start=$(( end + 1 ))
    g=$(( g + 1 ))
  done
}

# Split golden run list into “round” files so that each round’s sum of jobs
# (sum over runs of ceil(nFiles/GROUP_SIZE)) ≤ MAX_JOBS
split_golden() {
  local gs="$1" cap="$2"
  say "Splitting golden run list → groups of ${BOLD}${gs}${RST} files/job, cap ${BOLD}${cap}${RST} jobs/round"
  [[ -s "$GOLDEN" ]] || { err "Golden run list not found or empty: $GOLDEN"; exit 4; }

  # Clean previous rounds for this dataset
  rm -f "${ROUND_DIR}/goldenRuns_${TAG}_segment"*.txt 2>/dev/null || true

  local seg=1 jobs_in_seg=0 runs_in_seg=0
  local segfile="${ROUND_DIR}/goldenRuns_${TAG}_segment${seg}.txt"
  : > "$segfile"

  local total_jobs=0 total_runs=0
  while IFS= read -r rn; do
    [[ -z "$rn" || "$rn" =~ ^# ]] && continue
    local r8; r8="$(run8 "$rn")"
    local lf="${LIST_DIR}/dst_calofitting-${r8}.list"
    if [[ ! -s "$lf" ]]; then
      warn "No per-run list for ${r8}; skipping"
      continue
    fi
    local nfiles; nfiles=$(wc -l < "$lf" | awk '{print $1}')
    local nj; nj=$(ceil_div "$nfiles" "$gs")

    # If adding this run would exceed the cap, flush segment
    if (( jobs_in_seg + nj > cap )) && (( runs_in_seg > 0 )); then
      say "Round ${BOLD}${seg}${RST}: runs=${runs_in_seg}, jobs=${jobs_in_seg} → ${segfile}"
      seg=$(( seg + 1 ))
      segfile="${ROUND_DIR}/goldenRuns_${TAG}_segment${seg}.txt"
      : > "$segfile"
      jobs_in_seg=0; runs_in_seg=0
    fi

    echo "$r8" >> "$segfile"
    jobs_in_seg=$(( jobs_in_seg + nj ))
    runs_in_seg=$(( runs_in_seg + 1 ))
    total_jobs=$(( total_jobs + nj ))
    total_runs=$(( total_runs + 1 ))
  done < "$GOLDEN"

  say "Round ${BOLD}${seg}${RST}: runs=${runs_in_seg}, jobs=${jobs_in_seg} → ${segfile}"
  say "Total runs=${BOLD}${total_runs}${RST}, total condensed jobs≈${BOLD}${total_jobs}${RST}"
}

# Submit a set of runs (from a round file or the whole golden list) to Condor
#   submit_condor <runs_source> [firstChunk]
submit_condor() {
  local source="$1" first_chunk="${2:-}"
  [[ -s "$source" ]] || { err "Run source not found or empty: $source"; exit 5; }

  local stamp; stamp="$(date +%Y%m%d_%H%M%S)"
  local sub="${BASE}/RecoilJets_${TAG}_${stamp}.sub"

  cat > "$sub" <<SUB
universe      = vanilla
executable    = ${EXE}
initialdir    = ${BASE}
getenv        = True
log           = ${LOG_DIR}/job.\$(Cluster).\$(Process).log
output        = ${OUT_DIR}/job.\$(Cluster).\$(Process).out
error         = ${ERR_DIR}/job.\$(Cluster).\$(Process).err
request_memory= 2000MB
should_transfer_files = NO
stream_output = True
stream_error  = True
SUB

  local queued=0
  while IFS= read -r r8; do
    [[ -z "$r8" || "$r8" =~ ^# ]] && continue
    # Build group lists for this run
    mapfile -t groups < <( make_groups "$r8" "$GROUP_SIZE" )
    (( ${#groups[@]} )) || { warn "No groups produced for run $r8; skipping"; continue; }

    # If "firstChunk" is requested, only submit the first group (smoke test)
    if [[ "$first_chunk" == "firstChunk" ]]; then
      groups=( "${groups[0]}" )
    fi

    local gidx=0
    for glist in "${groups[@]}"; do
      (( gidx+=1 ))
      # args: <run8> <chunkList> <isPP|isAuAu> <Cluster> <nEvents=0> <chunkIdx> NONE <destBase>
      printf 'arguments = %s %s %s $(Cluster) 0 %d NONE %s\nqueue\n\n' \
             "$r8" "$glist" "$DATASET" "$gidx" "$DEST_BASE" >> "$sub"
      (( queued+=1 ))
    done
  done < "$source"

  say "Submitting ${BOLD}${queued}${RST} jobs → $(basename "$sub")"
  need_cmd condor_submit
  condor_submit "$sub"
}

# ------------------------ Parse CLI ------------------------
[[ $# -ge 1 ]] || usage
resolve_dataset "$1"
ACTION="${2:-condor}"

# Allow optional overrides: groupSize N   maxJobs M
shift_ct=2
while [[ $# -gt $shift_ct ]]; do
  case "${!shift_ct}" in
    groupSize) GROUP_SIZE="${!((shift_ct+1))}"; shift_ct=$((shift_ct+2)) ;;
    maxJobs)   MAX_JOBS="${!((shift_ct+1))}";   shift_ct=$((shift_ct+2)) ;;
    *) break ;;
  esac
done

say "Dataset=${BOLD}${DATASET}${RST}  Tag=${TAG}"
say "Macro=${MACRO}"
say "Input lists dir=${LIST_DIR}"
say "Golden runs=${GOLDEN}"
say "Stage dir=${STAGE_DIR}"
say "Rounds dir=${ROUND_DIR}"
say "Dest base=${DEST_BASE}"
say "Group size=${GROUP_SIZE}, Max jobs/round=${MAX_JOBS}"
echo

# ------------------------ Actions --------------------------
case "$ACTION" in
  local)
    # Run locally on the FIRST file of the FIRST run, for N events (default LOCAL_EVENTS)
    nevt="${3:-$LOCAL_EVENTS}"
    say "Local test on ${DATASET}  (events=${nevt})"

    [[ -s "$GOLDEN" ]] || { err "Golden list is empty: $GOLDEN"; exit 6; }
    first_run="$(grep -m1 -E '^[0-9]+' "$GOLDEN" | head -n1 || true)"
    [[ -n "$first_run" ]] || { err "No run found in $GOLDEN"; exit 6; }
    r8="$(run8 "$first_run")"

    src="${LIST_DIR}/dst_calofitting-${r8}.list"
    [[ -s "$src" ]] || { err "Per-run list missing or empty: $src"; exit 7; }

    # Build a 1-file list (first file only)
    tmp="${STAGE_DIR}/run${r8}_LOCAL_firstfile.list"
    head -n 1 "$src" > "$tmp"
    say "Temp list → $tmp"
    say "Invoking wrapper locally…"
    RJ_DATASET="$DATASET" bash "$EXE" "$r8" "$tmp" "$DATASET" LOCAL "$nevt" 1 NONE "$DEST_BASE"
    ;;

  splitGoldenRunList)
    split_golden "$GROUP_SIZE" "$MAX_JOBS"
    ;;

  condor)
    # condor round <N> [firstChunk]    or    condor all
    submode="${3:-}"
    case "$submode" in
      round)
        seg="${4:?round number required}"
        firstChunk="${5:-}"
        round_file="${ROUND_DIR}/goldenRuns_${TAG}_segment${seg}.txt"
        [[ -s "$round_file" ]] || { err "Round file not found: $round_file. Run 'splitGoldenRunList' first."; exit 8; }
        submit_condor "$round_file" "$firstChunk"
        ;;
      all|"")
        # Submit all runs currently in GOLDEN (no round cap)
        # Write a temp "all runs" source file
        tmp_src="${ROUND_DIR}/ALL_${TAG}_$(date +%s).txt"
        grep -E '^[0-9]+' "$GOLDEN" > "$tmp_src"
        submit_condor "$tmp_src"
        ;;
      *)
        usage
        ;;
    esac
    ;;

  *)
    usage
    ;;
esac

say "Done."
