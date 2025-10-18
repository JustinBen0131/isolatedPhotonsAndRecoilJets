#!/usr/bin/env bash
###############################################################################
# RecoilJets_Condor_submit.sh — one-stop driver for LOCAL testing and CONDOR
# batch submission of the sPHENIX RecoilJets analysis (DATA ONLY; no MC).
#
# OVERVIEW
#   • Two datasets are supported: isPP and isAuAu.
#   • Input comes from per‑run *.list files (one ROOT path per line).
#   • Jobs never mix files across runs. Each job processes a fixed-size group
#     of files from a single run (grouping is configurable).
#   • Optional “rounds” (segments) let you cap the total number of jobs
#     per submit wave; you can submit a specific round later.
#
# DATASET TOKENS (case-insensitive)
#   isPP  | pp  | PP    → proton-proton
#   isAuAu| auau| AA    → Au+Au
#
# DIRECTORY LAYOUT (fixed)
#   BASE          = /sphenix/u/patsfan753/scratch/thesisAnalysis
#   MACRO         = ${BASE}/macros/Fun4All_recoilJets.C   # for display only here
#   EXE (wrapper) = ${BASE}/RecoilJets_Condor.sh          # actually runs the macro
#
# INPUTS YOU PROVIDE
#   Golden run lists:
#     • PP  : ${BASE}/Full_ppGoldenRunList_Version3.list
#     • AuAu: ${BASE}/Full_AuAuGoldenRunList_Version1_personal.list
#   Per-run file lists (one per run8):
#     • PP  : ${BASE}/dst_lists_pp/dst_calofitting-<run8>.list
#     • AuAu: ${BASE}/dst_lists_auau/dst_calofitting-<run8>.list
#
# ARTIFACTS THIS SCRIPT CREATES (and where)
#   Grouped chunk lists (per run) — used by Condor jobs:
#     ${BASE}/condor_lists/<pp|auau>/run<run8>_grpNNN.list
#     • CLEANING: When grouping a given run, all existing run<run8>_grp*.list
#       in that dataset’s condor_lists directory are removed first, then
#       regenerated according to the current groupSize.
#
#   Round (segment) files — lists of run8 identifiers that define a “round”:
#     ${BASE}/condor_segments/<pp|auau>/goldenRuns_<tag>_segment<k>.txt
#     • CLEANING: When you run `splitGoldenRunList`, all existing
#       goldenRuns_<tag>_segment*.txt for that dataset are deleted and
#       re-created using the current groupSize and maxJobs cap.
#
#   Condor submit files (auto-named):
#     ${BASE}/RecoilJets_<tag>_YYYYMMDD_HHMMSS.sub
#     • Settings baked in: request_memory=2000MB, should_transfer_files=NO,
#       stream_output/error=True, getenv=True.
#     • One queued job per chunk list line (i.e., one job per group).
#     • Environment exported to each job:
#         RJ_DATASET=<isPP|isAuAu>   (used by the macro to set data type)
#         RJ_VERBOSITY=0             (quiet macro in batch by default)
#
# LOG/STDOUT/STDERR DIRECTORIES (auto-created)
#   ${BASE}/log     — Condor *.log
#   ${BASE}/stdout  — Condor *.out
#   ${BASE}/error   — Condor *.err
#
# OUTPUT ROOT FILE DESTINATION TREE (used by the wrapper EXE)
#   • PP  → /sphenix/tg/tg01/bulk/jbennett/thesisAna/pp/<run8>/
#   • AuAu→ /sphenix/tg/tg01/bulk/jbennett/thesisAna/auau/<run8>/
#
# VERBOSITY POLICY (propagates to the Fun4All macro via RJ_VERBOSITY)
#   • LOCAL mode: default RJ_VERBOSITY=10 (chatty). You can override with:
#       - trailing token  VERBOSE=N   (preferred), or
#       - environment     VERBOSE=N   (fallback)
#     The script exports RJ_VERBOSITY accordingly to the wrapper.
#   • CONDOR mode: RJ_VERBOSITY is forced to 0 in the submit file.
#
# COMMANDS (argument grammar and exact behavior)
#
#   1) LOCAL quick test on the FIRST RUN (first file only)
#      $ ./RecoilJets_Condor_submit.sh <isPP|isAuAu> local [Nevents] [VERBOSE=N]
#        • Picks the first uncommented run from the dataset’s golden list.
#        • Uses only the first file of that run’s per‑run list.
#        • Creates a temporary chunk list:
#            ${BASE}/condor_lists/<tag>/run<run8>_LOCAL_firstfile.list
#        • Invokes the wrapper with:
#            RJ_DATASET=<dataset>, RJ_VERBOSITY=<10 or override>,
#            nevents=<Nevents or 5000 default>.
#
#   2) PREPARE ROUND FILES (segment the golden list by job budget)
#      $ ./RecoilJets_Condor_submit.sh <isPP|isAuAu> splitGoldenRunList \
#            [groupSize N] [maxJobs M]
#        • Computes number of jobs per run as ceil(nFiles / groupSize).
#        • Sequentially fills round<k> until adding next run would exceed cap M,
#          then starts a new round file.
#        • Writes: ${BASE}/condor_segments/<tag>/goldenRuns_<tag>_segment<k>.txt
#        • CLEANING: deletes all previous segment*.txt for this dataset first.
#
#   3) SUBMIT A SPECIFIC ROUND TO CONDOR
#      $ ./RecoilJets_Condor_submit.sh <isPP|isAuAu> condor round <K> \
#            [groupSize N] [firstChunk]
#        • Reads runs from the prebuilt round file <K>.
#        • Regroups each run’s per‑run list using the (possibly overridden)
#          groupSize and queues one job per group.
#        • If ‘firstChunk’ is provided, only the first group of each run is
#          submitted (smoke test).
#        • RJ_DATASET is exported; RJ_VERBOSITY is forced to 0.
#        • NOTE: maxJobs is IGNORED here — it’s only used by splitGoldenRunList.
#
#   4) SUBMIT ALL RUNS CURRENTLY IN THE GOLDEN LIST
#      $ ./RecoilJets_Condor_submit.sh <isPP|isAuAu> condor all [groupSize N]
#        • Builds a temporary ALL‑runs source from the golden list and submits.
#        • Regroups with the given groupSize if provided.
#        • RJ_DATASET is exported; RJ_VERBOSITY is forced to 0.
#        • NOTE: ‘firstChunk’ is NOT supported in ‘condor all’.
#        • NOTE: maxJobs is IGNORED here.
#
# OPTIONS (placement matters)
#   • Place [groupSize N] and [maxJobs M] AFTER the subcommand, exactly as in
#     the examples below (the parser expects this order).
#
# EXAMPLES
#   • Local quick check (PP, first run’s first file), 5000 events, verbose=10:
#       ./RecoilJets_Condor_submit.sh isPP  local 5000
#   • Local with explicit verbosity 4 (events default to 5000 when VERBOSE=N
#     is given as the 3rd token):
#       ./RecoilJets_Condor_submit.sh isAuAu local VERBOSE=4
#   • Make round files for AuAu with 3 files/job and cap of 15000 jobs/round:
#       ./RecoilJets_Condor_submit.sh isAuAu splitGoldenRunList groupSize 3 maxJobs 15000
#   • Submit AuAu round 1 (uses existing segment file):
#       ./RecoilJets_Condor_submit.sh isAuAu condor round 1
#   • Submit PP round 2 but only the first chunk per run:
#       ./RecoilJets_Condor_submit.sh isPP  condor round 2 firstChunk
#   • Submit ALL PP runs with 4 files per job:
#       ./RecoilJets_Condor_submit.sh isPP  condor all groupSize 4
#
# REQUIREMENTS
#   • condor_submit must be on PATH for Condor submissions.
#   • Golden lists and per‑run list files must exist as described above.
###############################################################################

set -euo pipefail

# ------------------------ Pretty printing ------------------
BOLD=$'\e[1m'; DIM=$'\e[2m'; RED=$'\e[31m'; YEL=$'\e[33m'; GRN=$'\e[32m'; BLU=$'\e[34m'; RST=$'\e[0m'
say()  { printf "${BLU}➜${RST} %s\n" "$*"; }
warn() { printf "${YEL}⚠ %s${RST}\n" "$*" >&2; }
err()  { printf "${RED}✘ %s${RST}\n" "$*" >&2; }

# ------------------------ Fixed paths ----------------------
BASE="/sphenix/u/patsfan753/scratch/thesisAnalysis"
MACRO="${BASE}/macros/Fun4All_recoilJets.C"
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
# Force dataset & quiet macro on Condor:
environment   = RJ_DATASET=${DATASET};RJ_VERBOSITY=0
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
    # Events and local verbosity policy:
    #   - default events: LOCAL_EVENTS
    #   - default verbosity: 10
    #   - allow override via a trailing "VERBOSE=N" token or environment VERBOSE
    nevt="${3:-$LOCAL_EVENTS}"
    local_v_override=""
    if [[ "${3:-}" =~ ^VERBOSE=([0-9]+)$ ]]; then
      local_v_override="${BASH_REMATCH[1]}"
      nevt="$LOCAL_EVENTS"
    elif [[ "${4:-}" =~ ^VERBOSE=([0-9]+)$ ]]; then
      local_v_override="${BASH_REMATCH[1]}"
    fi
    RJV="${local_v_override:-${VERBOSE:-10}}"
    say "Local test on ${DATASET}  (events=${nevt}, RJ_VERBOSITY=${RJV})"

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
    RJ_DATASET="$DATASET" RJ_VERBOSITY="$RJV" bash "$EXE" "$r8" "$tmp" "$DATASET" LOCAL "$nevt" 1 NONE "$DEST_BASE"
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
