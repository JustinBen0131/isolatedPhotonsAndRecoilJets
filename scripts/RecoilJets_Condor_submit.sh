#!/usr/bin/env bash
###############################################################################
# RecoilJets_Condor_submit.sh — one-stop driver for LOCAL testing and CONDOR
#
# BASELINE SIM COMMAND TO SUBMIT ALL JOBS I CARE ABOUT IN SIM FOR 10 AND 20 samples --> ./RecoilJets_Condor_submit.sh isSim condorDoAll
#
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
MAX_JOBS=50000      # job budget per round file
LOCAL_EVENTS=5000   # default N for "local" if not given
TRIGGER_BIT=""      # optional: filter runs by GL1 scaledown bit (e.g., TRIGGER=26)

# ------------------------ Simulation (MC) defaults --------
# isSim expects a master list with ONE line per segment, whitespace-separated:
#   col1 = DST_CALO_CLUSTER
#   col2 = G4Hits
#   col3 = DST_JETS        (truth jets DST)
#   col4 = DST_GLOBAL      (GlobalVertexMap lives here)
#   col5 = DST_MBD_EPD     (MBD inputs; needed for reco MBD vertex)
SIM_ROOT="${BASE}/simListFiles"
SIM_SAMPLE_DEFAULT="run28_photonjet10"
SIM_SAMPLE="${SIM_SAMPLE_DEFAULT}"

# Name for the staged 5-column master list (built by sim_init from matched lists)
SIM_MASTER5_NAME="DST_SIM_MASTER_5COL.list"

# Output dir root for sim is: ${SIM_DEST_BASE}/<cfgTag>/<SIM_SAMPLE>
SIM_DEST_BASE="/sphenix/tg/tg01/bulk/jbennett/thesisAna/sim"

# YAML config (Fun4All default is next to the macro unless RJ_CONFIG_YAML is set)
SIM_YAML_DEFAULT="${BASE}/macros/analysis_config.yaml"
SIM_YAML_OVERRIDE_DIR="${BASE}/condor_yaml_overrides"
SIM_CFG_TAG=""

# ------------------------ Helpers --------------------------
usage() {
  cat <<USAGE
${BOLD}Usage:${RST}

${BOLD}DATA modes:${RST}
  ${BOLD}$0 <isPP|isAuAu> local [Nevents] [VERBOSE=N]${RST}
  ${BOLD}$0 <isPP|isAuAu> CHECKJOBS [groupSize N]${RST}
  ${BOLD}$0 <isPP|isAuAu> splitGoldenRunList [groupSize N] [maxJobs M]${RST}
  ${BOLD}$0 <isPP|isAuAu> condor testJob${RST}
  ${BOLD}$0 <isPP|isAuAu> condor round <N> [groupSize N] [firstChunk]${RST}
  ${BOLD}$0 <isPP|isAuAu> condor all [groupSize N]${RST}

${BOLD}SIM mode:${RST}
  ${BOLD}$0 isSim local [Nevents] [VERBOSE=N] [SAMPLE=run28_photonjet10]${RST}
  ${BOLD}$0 isSim CHECKJOBS [groupSize N] [SAMPLE=run28_photonjet10]${RST}
  ${BOLD}$0 isSim condorTest [SAMPLE=run28_photonjet10]${RST}
  ${BOLD}$0 isSim condorDoAll [groupSize N] [SAMPLE=run28_photonjet10]${RST}

Examples:
  $0 isPP  local 5000
  $0 isAuAu splitGoldenRunList groupSize 3 maxJobs 12000
  $0 isAuAu condor round 1
  $0 isPP  condor all groupSize 4
  $0 isPP  CHECKJOBS groupSize 4

  $0 isSim CHECKJOBS groupSize 5
  $0 isSim local 5000
  $0 isSim condorTest
  $0 isSim condorDoAll groupSize 5
USAGE
  exit 2
}

need_cmd() { command -v "$1" >/dev/null 2>&1 || { err "Missing required command: $1"; exit 3; }; }

# ------------------------ SIM YAML sweep helpers --------------------------
trim_ws() {
  local s="$1"
  s="${s#"${s%%[![:space:]]*}"}"
  s="${s%"${s##*[![:space:]]}"}"
  printf "%s" "$s"
}

sim_yaml_master_path() {
  if [[ -n "${RJ_CONFIG_YAML:-}" ]]; then
    local p; p="$(trim_ws "${RJ_CONFIG_YAML}")"
    [[ -n "$p" ]] && { echo "$p"; return; }
  fi
  echo "$SIM_YAML_DEFAULT"
}

yaml_get_values() {
  local key="$1" file="$2"
  local line rhs inside
  line=$(grep -E "^[[:space:]]*${key}:" "$file" | head -n 1 || true)
  [[ -n "$line" ]] || { err "Missing YAML key '${key}' in ${file}"; exit 70; }
  line="${line%%#*}"
  rhs="${line#*:}"
  rhs="$(trim_ws "$rhs")"

  if [[ "$rhs" == \[*\] ]]; then
    inside="${rhs#\[}"
    inside="${inside%\]}"
    IFS=',' read -ra parts <<< "$inside"
    for p in "${parts[@]}"; do
      p="$(trim_ws "$p")"
      [[ -n "$p" ]] && echo "$p"
    done
  else
    [[ -n "$rhs" ]] && echo "$rhs"
  fi
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

sim_b2b_tag() {
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

sim_make_yaml_override() {
  local master="$1" pt="$2" frac="$3" tag="$4"
  mkdir -p "$SIM_YAML_OVERRIDE_DIR"
  local out="${SIM_YAML_OVERRIDE_DIR}/analysis_config_${tag}.yaml"

  [[ -s "$master" ]] || { err "Master YAML not found or empty: $master"; exit 71; }
  grep -Eq '^[[:space:]]*jet_pt_min:' "$master" || { err "YAML missing key 'jet_pt_min' in $master"; exit 71; }
  grep -Eq '^[[:space:]]*back_to_back_dphi_min_pi_fraction:' "$master" || { err "YAML missing key 'back_to_back_dphi_min_pi_fraction' in $master"; exit 71; }

  sed -E \
    -e "s|^([[:space:]]*jet_pt_min:).*|\\1 ${pt}|" \
    -e "s|^([[:space:]]*back_to_back_dphi_min_pi_fraction:).*|\\1 ${frac}|" \
    "$master" > "$out"

  echo "$out"
}

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
    isSim|sim|SIM)
      DATASET="isSim"
      GOLDEN=""        # not used in sim mode
      LIST_DIR=""      # not used in sim mode
      DEST_BASE="$SIM_DEST_BASE"   # output dir will be DEST_BASE/SIM_SAMPLE
      TAG="sim"
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

# Check if a GL1 trigger bit is active for a run (scaledown != -1).
# Usage: is_trigger_active <run8> <bit> ; returns 0 if active, 1 otherwise.
is_trigger_active() {
  local run8="$1" bit="$2"
  [[ -z "$bit" ]] && return 0
  local run_dec=$((10#$run8))
  local val
  val=$(psql -h sphnxdaqdbreplica -d daq -At -q -F $'\t' \
        -c "SELECT scaledown${bit} FROM gl1_scaledown WHERE runnumber=${run_dec};" \
        | tr -d '[:space:]')
  [[ -n "$val" && "$val" != "-1" ]]
}

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

# ------------------------ Simulation helpers --------------------------
# Initializes paths for isSim mode and prepares a cleaned master list.
sim_init() {
  SIM_DIR="${SIM_ROOT}/${SIM_SAMPLE}"
  [[ -d "$SIM_DIR" ]] || { err "Sim sample directory not found: $SIM_DIR"; exit 20; }

  SIM_STAGE_DIR="${STAGE_DIR}/${SIM_SAMPLE}"
  mkdir -p "$SIM_STAGE_DIR"

  # Build a staged 5-column master list from the matched lists created by makeThesisSimLists.sh
  local calo="${SIM_DIR}/DST_CALO_CLUSTER.matched.list"
  local g4="${SIM_DIR}/G4Hits.matched.list"
  local jets="${SIM_DIR}/DST_JETS.matched.list"
  local glob="${SIM_DIR}/DST_GLOBAL.matched.list"
  local mbd="${SIM_DIR}/DST_MBD_EPD.matched.list"

  [[ -s "$calo" ]] || { err "Missing: $calo"; err "Run makeThesisSimLists.sh for ${SIM_SAMPLE}"; exit 23; }
  [[ -s "$g4"   ]] || { err "Missing: $g4";   err "Run makeThesisSimLists.sh for ${SIM_SAMPLE}"; exit 23; }
  [[ -s "$jets" ]] || { err "Missing: $jets"; err "Run makeThesisSimLists.sh for ${SIM_SAMPLE}"; exit 23; }
  [[ -s "$glob" ]] || { err "Missing: $glob"; err "Run makeThesisSimLists.sh for ${SIM_SAMPLE}"; exit 23; }
  [[ -s "$mbd"  ]] || { err "Missing: $mbd";  err "Run makeThesisSimLists.sh for ${SIM_SAMPLE}"; exit 23; }

  SIM_MASTER_LIST="${SIM_STAGE_DIR}/sim_${SIM_SAMPLE}_${SIM_MASTER5_NAME}"
  paste "$calo" "$g4" "$jets" "$glob" "$mbd" > "$SIM_MASTER_LIST"

  # Clean master list: strip blank lines + comment lines (keeps ALL columns)
  SIM_CLEAN_LIST="${SIM_STAGE_DIR}/sim_${SIM_SAMPLE}_PAIR_MASTER_CLEAN.list"
  grep -E -v '^[[:space:]]*($|#)' "$SIM_MASTER_LIST" > "$SIM_CLEAN_LIST" || true
  [[ -s "$SIM_CLEAN_LIST" ]] || { err "Sim master list empty after cleaning: $SIM_MASTER_LIST"; exit 22; }

  SIM_OUT_DIR="${DEST_BASE}/${SIM_SAMPLE}"
  mkdir -p "$SIM_OUT_DIR"

  if [[ -n "${SIM_CFG_TAG:-}" ]]; then
    SIM_JOB_PREFIX="sim_${SIM_SAMPLE}_${SIM_CFG_TAG}"
  else
    SIM_JOB_PREFIX="sim_${SIM_SAMPLE}"
  fi
}

# Build grouped chunk lists for isSim (one job per chunk list)
make_sim_groups() {
  local gs="$1"
  sim_init

  rm -f "${SIM_STAGE_DIR}/${SIM_JOB_PREFIX}_grp"*.list 2>/dev/null || true

  local nfiles; nfiles=$(wc -l < "$SIM_CLEAN_LIST" | awk '{print $1}')
  local ngroups; ngroups=$(ceil_div "$nfiles" "$gs")

  local start=1
  local g=1
  while (( g <= ngroups )); do
    local end=$(( start + gs - 1 ))
    local out="${SIM_STAGE_DIR}/${SIM_JOB_PREFIX}_grp$(printf "%03d" "$g").list"
    if (( end > nfiles )); then end="$nfiles"; fi
    sed -n "${start},${end}p" "$SIM_CLEAN_LIST" > "$out"
    echo "$out"
    start=$(( end + 1 ))
    g=$(( g + 1 ))
  done
}

# Dry-run job count for isSim
check_jobs_sim() {
  local gs="$GROUP_SIZE"
  sim_init

  local nfiles; nfiles=$(wc -l < "$SIM_CLEAN_LIST" | awk '{print $1}')
  local njobs; njobs=$(ceil_div "$nfiles" "$gs")

  say "CHECKJOBS (dataset=isSim, sample=${SIM_SAMPLE})"
  say "  groupSize         : ${gs}"
  say "  total input files : ${nfiles}"
  say "  total jobs (all)  : ${njobs}"
  say "  master list       : ${SIM_MASTER_LIST}"
}

# Wipe previous sim outputs + sim-tagged logs/out/err (ONLY used before condorDoAll)
wipe_sim_artifacts() {
  sim_init

  say "WIPING previous isSim artifacts for sample ${BOLD}${SIM_SAMPLE}${RST}"
  say "  Output ROOT dir : ${SIM_OUT_DIR}"
  say "  Logs prefix     : ${SIM_JOB_PREFIX}"

  rm -f "${SIM_OUT_DIR}/RecoilJets_"*.root 2>/dev/null || true
  rm -f "${SIM_STAGE_DIR}/${SIM_JOB_PREFIX}_grp"*.list 2>/dev/null || true
  rm -f "${SIM_STAGE_DIR}/${SIM_JOB_PREFIX}_LOCAL_"*.list 2>/dev/null || true
  rm -f "${SIM_STAGE_DIR}/${SIM_JOB_PREFIX}_condorTest_"*.list 2>/dev/null || true

  rm -f "${LOG_DIR}/${SIM_JOB_PREFIX}.job."*.log 2>/dev/null || true
  rm -f "${OUT_DIR}/${SIM_JOB_PREFIX}.job."*.out 2>/dev/null || true
  rm -f "${ERR_DIR}/${SIM_JOB_PREFIX}.job."*.err 2>/dev/null || true
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

# Dry-run job count for a hypothetical "condor all" at current GROUP_SIZE
#   check_jobs_all
# Reads the golden list and per-run lists; prints totals only. No side effects.
check_jobs_all() {
  local gs="$GROUP_SIZE"
  [[ -s "$GOLDEN" ]] || { err "Golden run list not found or empty: $GOLDEN"; exit 4; }

  local total_runs=0 total_files=0 total_jobs=0 missing=0

  while IFS= read -r rn; do
    [[ -z "$rn" || "$rn" =~ ^# ]] && continue
    local r8; r8="$(run8 "$rn")"
    local lf="${LIST_DIR}/dst_calofitting-${r8}.list"

    if [[ ! -s "$lf" ]]; then
      warn "No per-run list for ${r8}; skipping"
      ((missing++))
      continue
    fi

    local nfiles; nfiles=$(wc -l < "$lf" | awk '{print $1}')
    local nj;     nj=$(ceil_div "$nfiles" "$gs")

    (( total_runs  += 1    ))
    (( total_files += nfiles ))
    (( total_jobs  += nj   ))
  done < "$GOLDEN"

  say "CHECKJOBS (dataset=${DATASET})"
  say "  groupSize         : ${gs}"
  say "  runs (with lists) : ${total_runs}"
  say "  total input files : ${total_files}"
  say "  total jobs (all)  : ${total_jobs}"
  (( missing > 0 )) && warn "runs skipped due to missing per-run lists: ${missing}"
}

# Submit a set of runs (from a round file or the whole golden list) to Condor
#   submit_condor <runs_source> [firstChunk]
submit_condor() {
  local source="$1"
  local first_chunk="${2:-}"
  local wipe_mode="${3:-never}"   # never | once | always

  [[ -s "$source" ]] || { err "Run source not found or empty: $source"; exit 5; }

  # -------------------------------------------------------------------
  # DATA output wipe policy (pp/auau only):
  #   wipe_mode=always -> wipe DEST_BASE every time before submit
  #   wipe_mode=once   -> wipe DEST_BASE only if stamp is not present
  #   wipe_mode=never  -> never wipe
  #
  # Stamp file lives in BASE and persists across rounds:
  #   ${BASE}/.RJ_OUTPUTS_WIPED_<tag>.stamp
  # -------------------------------------------------------------------
  if [[ "$DATASET" != "isSim" ]]; then
    local stamp_file="${BASE}/.RJ_OUTPUTS_WIPED_${TAG}.stamp"
    local do_wipe=0

    if [[ "$wipe_mode" == "always" ]]; then
      do_wipe=1
    elif [[ "$wipe_mode" == "once" ]]; then
      [[ -f "$stamp_file" ]] || do_wipe=1
    fi

    if (( do_wipe )); then
      # Safety checks to avoid catastrophic wipes
      [[ -n "${DEST_BASE:-}" ]] || { err "DEST_BASE is empty; refusing to wipe."; exit 60; }
      [[ "$DEST_BASE" != "/" ]] || { err "DEST_BASE is '/' ; refusing to wipe."; exit 60; }

      case "$DEST_BASE" in
        */thesisAna/pp|*/thesisAna/auau) ;;
        *)
          err "Refusing to wipe DEST_BASE='$DEST_BASE' (not an expected thesisAna/{pp|auau} path)"
          exit 61
          ;;
      esac

      say "WIPING OUTPUT TREE (dataset=${DATASET}, mode=${wipe_mode}) → ${DEST_BASE}"
      mkdir -p "$DEST_BASE"

      # Remove everything *inside* DEST_BASE (files + subdirs), keep DEST_BASE itself
      find "$DEST_BASE" -mindepth 1 -maxdepth 1 -exec rm -rf {} + 2>/dev/null || true

      date +"%Y-%m-%d %H:%M:%S" > "$stamp_file"
      say "Wipe complete. Stamp written: $stamp_file"
      say "To force another wipe later, delete the stamp: rm -f $stamp_file"
    else
      if [[ "$wipe_mode" != "never" ]]; then
        say "Not wiping outputs (mode=${wipe_mode}); stamp exists: ${BASE}/.RJ_OUTPUTS_WIPED_${TAG}.stamp"
      fi
    fi
  fi

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
  local run_counter=0
  local t0; t0="$(date +%s)"

  # Controls:
  #   RJ_SUBMIT_TRACE=1               -> print one line per run
  #   RJ_SUBMIT_PROGRESS_EVERY=25     -> print progress every N runs (default 25)
  local submit_trace="${RJ_SUBMIT_TRACE:-0}"
  local progress_every="${RJ_SUBMIT_PROGRESS_EVERY:-25}"
  [[ "$progress_every" =~ ^[0-9]+$ ]] || progress_every=25
  (( progress_every > 0 )) || progress_every=25

  say "Building Condor submit file: $(basename "$sub")"
  say "Progress: export RJ_SUBMIT_TRACE=1 for per-run lines (or set RJ_SUBMIT_PROGRESS_EVERY=N)"

  while IFS= read -r r8_raw; do
    [[ -z "$r8_raw" || "$r8_raw" =~ ^# ]] && continue
    r8="$(run8 "$r8_raw")"   # ← zero-pad to 8 digits
    (( run_counter+=1 ))

    if (( submit_trace )); then
      say "[submit] run ${r8} (run #${run_counter})"
    elif (( run_counter % progress_every == 0 )); then
      local now elapsed
      now="$(date +%s)"
      elapsed=$(( now - t0 ))
      say "[submit] processed ${run_counter} runs, queued ${queued} jobs so far (elapsed ${elapsed}s)"
    fi

    # Optional trigger filter: require GL1 scaledown bit to be active
    if [[ -n "${TRIGGER_BIT}" ]]; then
      if ! is_trigger_active "$r8" "$TRIGGER_BIT"; then
        say "Skipping run ${r8} (trigger bit ${TRIGGER_BIT} not active)"
        continue
      fi
    fi

    # Build group lists for this run
    mapfile -t groups < <( make_groups "$r8" "$GROUP_SIZE" )

    (( ${#groups[@]} )) || { warn "No groups produced for run $r8; skipping"; continue; }

    if (( submit_trace )); then
      say "    files per job=${GROUP_SIZE} -> groups produced=${#groups[@]}"
    fi

    # If "firstChunk" is requested, only submit the first group (smoke test)
    if [[ "$first_chunk" == "firstChunk" ]]; then
      groups=( "${groups[0]}" )
      (( submit_trace )) && say "    firstChunk enabled -> submitting 1 group for run ${r8}"
    fi

    local gidx=0
    for glist in "${groups[@]}"; do
      (( gidx+=1 ))
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

# Parse remaining tokens (order-agnostic):
ACTION=""
SIM_SAMPLE_EXPLICIT=0
GROUP_SIZE_EXPLICIT=0
tokens=( "${@:2}" )
for (( idx=0; idx<${#tokens[@]}; idx++ )); do
  tok="${tokens[$idx]}"
  case "$tok" in
    local|condor|splitGoldenRunList|condorTest|condorDoAll)
      ACTION="$tok"
      ;;
    groupSize)
      next=$((idx+1))
      [[ $next -lt ${#tokens[@]} ]] || { err "groupSize requires a value"; exit 2; }
      GROUP_SIZE="${tokens[$next]}"
      GROUP_SIZE_EXPLICIT=1
      idx=$next
      ;;
    maxJobs)
      next=$((idx+1))
      [[ $next -lt ${#tokens[@]} ]] || { err "maxJobs requires a value"; exit 2; }
      MAX_JOBS="${tokens[$next]}"
      idx=$next
      ;;
    TRIGGER=*)
      TRIGGER_BIT="${tok#TRIGGER=}"
      [[ "$TRIGGER_BIT" =~ ^[0-9]+$ ]] || { err "TRIGGER must be an integer bit index (e.g., TRIGGER=26)"; exit 2; }
      ;;
    SAMPLE=*)
      SIM_SAMPLE="${tok#SAMPLE=}"
      SIM_SAMPLE_EXPLICIT=1
      ;;
    CHECKJOBS)
      ACTION="CHECKJOBS"
      ;;
    *)
      :  # ignore unrecognized tokens
      ;;
  esac
done

# Default action if none provided
if [[ -z "$ACTION" ]]; then
  if [[ "$DATASET" == "isSim" ]]; then
    ACTION="CHECKJOBS"
  else
    ACTION="condor"
  fi
fi

# If a trigger filter is requested, require psql
if [[ -n "${TRIGGER_BIT}" ]]; then
  need_cmd psql
fi

# Only print the full banner when we're actually running work (not CHECKJOBS)
if [[ "$ACTION" != "CHECKJOBS" ]]; then
  say "Dataset=${BOLD}${DATASET}${RST}  Tag=${TAG}"
  say "Macro=${MACRO}"
  say "Input lists dir=${LIST_DIR}"
  say "Golden runs=${GOLDEN}"
  say "Stage dir=${STAGE_DIR}"
  say "Rounds dir=${ROUND_DIR}"
  say "Dest base=${DEST_BASE}"
  if [[ -n "${TRIGGER_BIT}" ]]; then
    say "Trigger filter      : bit=${TRIGGER_BIT} (scaledown != -1 required)"
  fi
  say "Group size=${GROUP_SIZE}, Max jobs/round=${MAX_JOBS}"
  echo
fi

# ------------------------ Actions --------------------------
case "$ACTION" in
  CHECKJOBS)
    # Dry-run only
    if [[ "$DATASET" == "isSim" ]]; then
      check_jobs_sim
    else
      check_jobs_all
    fi
    exit 0
    ;;

  local)
    # Parse optional [Nevents] and VERBOSE=N from tokens after 'local'
    nevt="$LOCAL_EVENTS"
    RJV="10"                     # default for local runs
    rest=( "${@:3}" )
    for t in "${rest[@]}"; do
      if [[ "$t" =~ ^VERBOSE=([0-9]+)$ ]]; then
        RJV="${BASH_REMATCH[1]}"
      elif [[ "$t" =~ ^[0-9]+$ ]]; then
        nevt="$t"
      fi
    done

    if [[ "$DATASET" == "isSim" ]]; then
      master_yaml="$(sim_yaml_master_path)"
      [[ -s "$master_yaml" ]] || { err "Master YAML not found or empty: $master_yaml"; exit 72; }

      mapfile -t sim_pts   < <( yaml_get_values "jet_pt_min" "$master_yaml" )
      mapfile -t sim_fracs < <( yaml_get_values "back_to_back_dphi_min_pi_fraction" "$master_yaml" )
      (( ${#sim_pts[@]} ))   || { err "No values found for jet_pt_min in $master_yaml"; exit 72; }
      (( ${#sim_fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $master_yaml"; exit 72; }

      pt0="${sim_pts[0]}"
      frac0="${sim_fracs[0]}"
      SIM_CFG_TAG="jetMinPt$(sim_pt_tag "$pt0")_$(sim_b2b_tag "$frac0")"
      DEST_BASE="${SIM_DEST_BASE}/${SIM_CFG_TAG}"
      yaml_override="$(sim_make_yaml_override "$master_yaml" "$pt0" "$frac0" "$SIM_CFG_TAG")"

      sim_init

      say "Local test on isSim sample=${SIM_SAMPLE} (tag=${SIM_CFG_TAG}, events=${nevt}, RJ_VERBOSITY=${RJV})"
      say "YAML override: ${yaml_override}"

      tmp="${SIM_STAGE_DIR}/${SIM_JOB_PREFIX}_LOCAL_firstfile.list"
      head -n 1 "$SIM_CLEAN_LIST" > "$tmp"
      say "Temp list → $tmp"
      say "Invoking wrapper locally…"

      RJ_VERBOSITY="$RJV" RJ_CONFIG_YAML="$yaml_override" bash "$EXE" "$SIM_SAMPLE" "$tmp" "isSim" LOCAL "$nevt" 1 NONE "$DEST_BASE"
    else
      say "Local test on ${DATASET}  (events=${nevt}, RJ_VERBOSITY=${RJV})"

      [[ -s "$GOLDEN" ]] || { err "Golden list is empty: $GOLDEN"; exit 6; }

      # Pick first golden run; if TRIGGER_BIT set, pick the first run with that bit active
      r8=""
      while IFS= read -r rn; do
        [[ -z "$rn" || "$rn" =~ ^# ]] && continue
        cand="$(run8 "$rn")"
        if [[ -n "${TRIGGER_BIT}" ]]; then
          if is_trigger_active "$cand" "$TRIGGER_BIT"; then r8="$cand"; break; fi
        else
          r8="$cand"; break
        fi
      done < "$GOLDEN"

      [[ -n "$r8" ]] || { err "No run in $GOLDEN satisfies the requested trigger filter (TRIGGER=${TRIGGER_BIT:-none})."; exit 6; }

      src="${LIST_DIR}/dst_calofitting-${r8}.list"
      [[ -s "$src" ]] || { err "Per-run list missing or empty: $src"; exit 7; }

      # Build a 1-file list (first file only)
      tmp="${STAGE_DIR}/run${r8}_LOCAL_firstfile.list"
      head -n 1 "$src" > "$tmp"
      say "Temp list → $tmp"
      say "Invoking wrapper locally…"
      RJ_DATASET="$DATASET" RJ_VERBOSITY="$RJV" bash "$EXE" "$r8" "$tmp" "$DATASET" LOCAL "$nevt" 1 NONE "$DEST_BASE"
    fi
    ;;

  condorTest)
    # One verbose Condor job on first sim file
    [[ "$DATASET" == "isSim" ]] || { err "condorTest is only valid for isSim"; exit 2; }
    need_cmd condor_submit

    master_yaml="$(sim_yaml_master_path)"
    [[ -s "$master_yaml" ]] || { err "Master YAML not found or empty: $master_yaml"; exit 72; }

    mapfile -t sim_pts   < <( yaml_get_values "jet_pt_min" "$master_yaml" )
    mapfile -t sim_fracs < <( yaml_get_values "back_to_back_dphi_min_pi_fraction" "$master_yaml" )
    (( ${#sim_pts[@]} ))   || { err "No values found for jet_pt_min in $master_yaml"; exit 72; }
    (( ${#sim_fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $master_yaml"; exit 72; }

    pt0="${sim_pts[0]}"
    frac0="${sim_fracs[0]}"
    SIM_CFG_TAG="jetMinPt$(sim_pt_tag "$pt0")_$(sim_b2b_tag "$frac0")"
    DEST_BASE="${SIM_DEST_BASE}/${SIM_CFG_TAG}"
    yaml_override="$(sim_make_yaml_override "$master_yaml" "$pt0" "$frac0" "$SIM_CFG_TAG")"

    sim_init

    tmp="${SIM_STAGE_DIR}/${SIM_JOB_PREFIX}_condorTest_firstfile.list"
    head -n 1 "$SIM_CLEAN_LIST" > "$tmp"

    stamp="$(date +%Y%m%d_%H%M%S)"
    sub="${BASE}/RecoilJets_sim_${SIM_CFG_TAG}_${SIM_SAMPLE}_${stamp}_TEST.sub"

    cat > "$sub" <<SUB
universe      = vanilla
executable    = ${EXE}
initialdir    = ${BASE}
getenv        = True
log           = ${LOG_DIR}/${SIM_JOB_PREFIX}.job.\$(Cluster).\$(Process).log
output        = ${OUT_DIR}/${SIM_JOB_PREFIX}.job.\$(Cluster).\$(Process).out
error         = ${ERR_DIR}/${SIM_JOB_PREFIX}.job.\$(Cluster).\$(Process).err
request_memory= 2000MB
should_transfer_files = NO
stream_output = True
stream_error  = True
environment   = RJ_VERBOSITY=10;RJ_CONFIG_YAML=${yaml_override}
arguments     = ${SIM_SAMPLE} ${tmp} isSim \$(Cluster) 0 1 NONE ${DEST_BASE}
queue
SUB

    say "Submitting 1 isSim condorTest job (sample=${SIM_SAMPLE}, tag=${SIM_CFG_TAG}) → $(basename "$sub")"
    say "Output ROOT dir: ${DEST_BASE}/${SIM_SAMPLE}"
    say "YAML override: ${yaml_override}"
    condor_submit "$sub"
    ;;

  condorDoAll)
    # Submit all sim files in grouped chunks. MUST wipe outputs/logs first.
    [[ "$DATASET" == "isSim" ]] || { err "condorDoAll is only valid for isSim"; exit 2; }
    need_cmd condor_submit

    gs_doall="$GROUP_SIZE"
    if [[ "${GROUP_SIZE_EXPLICIT:-0}" -eq 0 ]]; then
      gs_doall="5"
    fi

    master_yaml="$(sim_yaml_master_path)"
    [[ -s "$master_yaml" ]] || { err "Master YAML not found or empty: $master_yaml"; exit 72; }

    mapfile -t sim_pts   < <( yaml_get_values "jet_pt_min" "$master_yaml" )
    mapfile -t sim_fracs < <( yaml_get_values "back_to_back_dphi_min_pi_fraction" "$master_yaml" )
    (( ${#sim_pts[@]} ))   || { err "No values found for jet_pt_min in $master_yaml"; exit 72; }
    (( ${#sim_fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $master_yaml"; exit 72; }

    samples=()
    if [[ "${SIM_SAMPLE_EXPLICIT:-0}" -eq 0 ]]; then
      samples=( "run28_photonjet10" "run28_photonjet20" )
    else
      samples=( "${SIM_SAMPLE}" )
    fi

    for pt in "${sim_pts[@]}"; do
      for frac in "${sim_fracs[@]}"; do
        SIM_CFG_TAG="jetMinPt$(sim_pt_tag "$pt")_$(sim_b2b_tag "$frac")"
        DEST_BASE="${SIM_DEST_BASE}/${SIM_CFG_TAG}"
        yaml_override="$(sim_make_yaml_override "$master_yaml" "$pt" "$frac" "$SIM_CFG_TAG")"

        for samp in "${samples[@]}"; do
          SIM_SAMPLE="$samp"
          GROUP_SIZE="$gs_doall"

          wipe_sim_artifacts

          mapfile -t groups < <( make_sim_groups "$GROUP_SIZE" )
          (( ${#groups[@]} )) || { err "No sim groups produced (sample=${SIM_SAMPLE}, tag=${SIM_CFG_TAG})"; exit 30; }

          stamp="$(date +%Y%m%d_%H%M%S)"
          sub="${BASE}/RecoilJets_sim_${SIM_CFG_TAG}_${SIM_SAMPLE}_${stamp}.sub"

          cat > "$sub" <<SUB
universe      = vanilla
executable    = ${EXE}
initialdir    = ${BASE}
getenv        = True
log           = ${LOG_DIR}/${SIM_JOB_PREFIX}.job.\$(Cluster).\$(Process).log
output        = ${OUT_DIR}/${SIM_JOB_PREFIX}.job.\$(Cluster).\$(Process).out
error         = ${ERR_DIR}/${SIM_JOB_PREFIX}.job.\$(Cluster).\$(Process).err
request_memory= 2000MB
should_transfer_files = NO
stream_output = True
stream_error  = True
environment   = RJ_VERBOSITY=0;RJ_CONFIG_YAML=${yaml_override}
SUB

          gidx=0
          for glist in "${groups[@]}"; do
            (( gidx+=1 ))
            printf 'arguments = %s %s %s $(Cluster) 0 %d NONE %s\nqueue\n\n' \
                   "$SIM_SAMPLE" "$glist" "isSim" "$gidx" "$DEST_BASE" >> "$sub"
          done

          say "Submitting isSim condorDoAll (tag=${SIM_CFG_TAG}, sample=${SIM_SAMPLE}, groupSize=${GROUP_SIZE}) → jobs=${BOLD}${#groups[@]}${RST}"
          say "Output ROOT dir: ${DEST_BASE}/${SIM_SAMPLE}"
          say "YAML override: ${yaml_override}"
          condor_submit "$sub"
        done
      done
    done
    ;;

  splitGoldenRunList)
    [[ "$DATASET" == "isSim" ]] && { err "splitGoldenRunList is not used in isSim mode"; exit 2; }
    split_golden "$GROUP_SIZE" "$MAX_JOBS"
    ;;

  condor)
    [[ "$DATASET" == "isSim" ]] && { err "Use: isSim condorTest | isSim condorDoAll (not 'condor')"; exit 2; }

    # DATA ONLY: condor testJob | condor round <N> [firstChunk] | condor all
    submode="${3:-}"
    case "$submode" in
      testJob)
        [[ -s "$GOLDEN" ]] || { err "Golden list is empty: $GOLDEN"; exit 6; }
        r8=""
        while IFS= read -r rn; do
          [[ -z "$rn" || "$rn" =~ ^# ]] && continue
          cand="$(run8 "$rn")"
          if [[ -n "${TRIGGER_BIT}" ]]; then
            if is_trigger_active "$cand" "$TRIGGER_BIT"; then r8="$cand"; break; fi
          else
            r8="$cand"; break
          fi
        done < "$GOLDEN"
        [[ -n "$r8" ]] || { err "No run in $GOLDEN satisfies the requested trigger filter (TRIGGER=${TRIGGER_BIT:-none})."; exit 6; }

        mapfile -t groups < <( make_groups "$r8" 1 )
        (( ${#groups[@]} )) || { err "No groups produced for run $r8"; exit 9; }
        glist="${groups[0]}"

        stamp="$(date +%Y%m%d_%H%M%S)"
        sub="${BASE}/RecoilJets_${TAG}_${stamp}_TEST.sub"

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
environment   = RJ_DATASET=${DATASET};RJ_VERBOSITY=10
arguments     = ${r8} ${glist} ${DATASET} \$(Cluster) 0 1 NONE ${DEST_BASE}
queue
SUB
        say "Submitting 1 test job on run ${BOLD}${r8}${RST} (first chunk, groupSize=1) → $(basename "$sub")"
        condor_submit "$sub"
        ;;
      round)
        seg="${4:?round number required}"
        firstChunk="${5:-}"
        round_file="${ROUND_DIR}/goldenRuns_${TAG}_segment${seg}.txt"
        [[ -s "$round_file" ]] || { err "Round file not found: $round_file. Run 'splitGoldenRunList' first."; exit 8; }

        # Wipe outputs only the FIRST time you submit any round (stamp-based)
        submit_condor "$round_file" "$firstChunk" "once"
        ;;
      all|"")
        say "Preparing CONDOR ALL submission (dataset=${DATASET}, groupSize=${GROUP_SIZE})"
        say "This step generates per-run grouped chunk lists and a large submit file before calling condor_submit."
        tmp_src="${ROUND_DIR}/ALL_${TAG}_$(date +%s).txt"
        grep -E '^[0-9]+' "$GOLDEN" > "$tmp_src"

        # Always wipe outputs before a full 'condor all' submit
        submit_condor "$tmp_src" "" "always"
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
