#!/usr/bin/env bash
###############################################################################
# mergeRecoilJets.sh  —  FAST REFERENCE
#
#./mergeRecoilJets.sh isSim firstRound groupSize 300 SAMPLE=run28_photonjet10
#./mergeRecoilJets.sh isSim secondRound SAMPLE=run28_photonjet10
#
#
# PURPOSE
#   Merge RecoilJets outputs produced by RecoilJets_Condor_submit.sh /
#   RecoilJets_Condor.sh. Supports:
#     • DATA (pp, auau): per-run merge → then grand-total merge
#     • SIM  (isSim):   two-step chunk merge (firstRound → secondRound)
#
# INPUT ROOT LOCATIONS (what you are merging)
#   DATA:
#     pp   : /sphenix/tg/tg01/bulk/jbennett/thesisAna/pp/<run8>/*.root
#     auau : /sphenix/tg/tg01/bulk/jbennett/thesisAna/auau/<run8>/*.root
#   SIM:
#     SAMPLE dir: /sphenix/tg/tg01/bulk/jbennett/thesisAna/<SAMPLE>/*.root
#       e.g.      /sphenix/tg/tg01/bulk/jbennett/thesisAna/run28_photonjet10/*.root
#
# OUTPUT MERGE LOCATION (where merged artifacts are written)
#   /sphenix/u/patsfan753/scratch/thesisAnalysis/output/<tag>/
#     DATA tags: pp | auau
#     SIM  tag : simTag = suffix after last "_" in SAMPLE
#               run28_photonjet10 → photonjet10
#
# ─────────────────────────────────────────────────────────────────────────────
# QUICK COMMANDS (copy/paste)
#
# DATA: inventory only (no merges)
#   ./mergeRecoilJets.sh checkFileOutput pp
#   ./mergeRecoilJets.sh checkFileOutput auau
#
# DATA: stage-1 per-run partials (Condor; one hadd job per run)
#   ./mergeRecoilJets.sh condor pp
#   ./mergeRecoilJets.sh condor pp test
#   ./mergeRecoilJets.sh condor auau firstHalf
#
# DATA: stage-2 grand total (merge all per-run partials)
#   ./mergeRecoilJets.sh addChunks pp
#   ./mergeRecoilJets.sh addChunks pp condor
#   ./mergeRecoilJets.sh addChunks auau
#   ./mergeRecoilJets.sh addChunks auau condor
#
# SIM: firstRound (Condor chunk-partials; groupSize = files per hadd)
#   ./mergeRecoilJets.sh isSim firstRound groupSize 300 SAMPLE=run28_photonjet10
#
# SIM: secondRound (final merge of chunk-partials)
#   ./mergeRecoilJets.sh isSim secondRound SAMPLE=run28_photonjet10
#   ./mergeRecoilJets.sh isSim secondRound condor SAMPLE=run28_photonjet10
#
# ─────────────────────────────────────────────────────────────────────────────
# WHAT EACH MODE PRODUCES
#
# DATA stage-1 outputs (per-run partials):
#   /.../output/<pp|auau>/chunkMerge_run_<run8>.root
#
# DATA stage-2 output (grand total):
#   /.../output/<pp|auau>/RecoilJets_<pp|auau>_ALL.root
#
# SIM firstRound outputs (chunk-partials):
#   /.../output/<simTag>/chunkMerge_<simTag>_grpNNN.root
#
# SIM secondRound output (final):
#   /.../output/<simTag>/RecoilJets_<simTag>_ALL.root
#
# ─────────────────────────────────────────────────────────────────────────────
# SAFETY / RESUME BEHAVIOR (DATA stage-1)
#   • This script NEVER runs condor_rm (it will not kill jobs).
#   • When building per-run hadd lists, it EXCLUDES any output ROOT file whose
#     producing RecoilJets_Condor.sh job is still present in condor_q
#     (IDLE/RUNNING/HELD/TRANSFERRING/SUSPENDED). This prevents merging
#     partially-produced outputs.
#
# TESTING (NO side effects)
#   DRYRUN=1      → NO deletions, NO hadd, NO condor_submit (prints the plan)
#   SKIP_TRACE=1  → prints per-run totals/busy/eligible counts
#   Example:
#     DRYRUN=1 SKIP_TRACE=1 ./mergeRecoilJets.sh condor pp test
#
# NOTES
#   • Only top-level *.root files are considered (maxdepth=1).
#   • Sorting uses: sort -V (natural ordering of chunk indices).
###############################################################################
set -euo pipefail

# DRYRUN=1  -> NO deletions, NO hadd, NO condor_submit. Only prints what would happen.
DRYRUN="${DRYRUN:-0}"

# SKIP_TRACE=1 -> print per-run (total/busy/eligible) summaries during planning/submission.
SKIP_TRACE="${SKIP_TRACE:-0}"

# ---------- Pretty printing ----------
BOLD=$'\e[1m'; RED=$'\e[31m'; YEL=$'\e[33m'; GRN=$'\e[32m'; BLU=$'\e[34m'; RST=$'\e[0m'
say()  { printf "${BLU}➜${RST} %s\n" "$*"; }
warn() { printf "${YEL}⚠ %s${RST}\n" "$*" >&2; }
err()  { printf "${RED}✘ %s${RST}\n" "$*" >&2; }

# ---------- Fixed dataset roots ----------
RUN_BASE_PP="/sphenix/tg/tg01/bulk/jbennett/thesisAna/pp"
RUN_BASE_AA="/sphenix/tg/tg01/bulk/jbennett/thesisAna/auau"

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

# NOTE: TMP_DIR is host-specific so merges run from sphnxuser02 and sphnxuser03
# cannot clobber each other's skiplists/listfiles when used close together.
TMP_DIR="${BASE}/tmp_recoil_merge_${HOSTTAG}"

# Where splitGoldenRunList round files live:
#   ${ROUND_BASE}/<pp|auau>/goldenRuns_<pp|auau>_segmentK.txt
ROUND_BASE="${BASE}/condor_segments"

CONDOR_EXEC="hadd_condor.sh"   # small wrapper emitted on-the-fly

# ---------- Naming ----------
PARTIAL_PREFIX="chunkMerge_run"          # per-run partial
FINAL_PREFIX="RecoilJets"                # final combined file prefix

# ---------- Helpers ----------
usage() {
  cat <<USAGE
${BOLD}Usage:${RST}
  $0 condor <pp|auau> [test|firstHalf]
  $0 addChunks <pp|auau> [condor]
  $0 checkFileOutput <pp|auau>

Examples:
  $0 condor pp
  $0 condor auau firstHalf
  $0 addChunks pp
  $0 addChunks auau condor
USAGE
  exit 2
}

need_cmd(){ command -v "$1" >/dev/null 2>&1 || { err "Missing command: $1"; exit 3; }; }

to_tag() {
  case "${1:-}" in
    pp|PP|isPP|PP_DATA|pp_data)   echo "pp" ;;
    auau|AA|isAuAu|AuAu|aa|AA_DATA|auau_data) echo "auau" ;;
    *) err "Dataset must be 'pp' or 'auau'"; exit 4 ;;
  esac
}

resolve_dataset() {
  TAG="$(to_tag "${1:-}")"
  case "$TAG" in
    pp)   RUN_BASE="$RUN_BASE_PP" ;;
    auau) RUN_BASE="$RUN_BASE_AA" ;;
  esac

  DEST_DIR="${OUT_BASE}/${TAG}"      # where partials and final live
  ROUND_DIR="${ROUND_BASE}/${TAG}"   # where splitGoldenRunList round files live

  mkdir -p "$DEST_DIR" "$ROUND_DIR" "$LOG_DIR" "$OUT_DIR" "$ERR_DIR" "$TMP_DIR"
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
  [[ "$TAG" == "auau" ]] && want="isAuAu"

  if command -v condor_q >/dev/null 2>&1; then
    (
      set +e +o pipefail
      condor_q "${USER:-$(id -un)}" \
        -constraint 'regexp("RecoilJets_Condor.sh",Cmd) && (JobStatus==1 || JobStatus==2 || JobStatus==5 || JobStatus==6 || JobStatus==7)' \
        -af Args 2>/dev/null |
      awk -v want="$want" '
        # Args format (from your submit):
        #   run8  chunkList  isPP|isAuAu  Cluster  0  grpIdx  NONE  destBase
        ($3 == want) {
          run  = $1
          lst  = $2
          ds   = $3
          dest = $NF
          if (run=="" || lst=="" || dest=="") next

          # chunk tag = basename(list) without .list
          n = lst
          sub(/^.*\//,"",n)
          sub(/\.list$/,"",n)

          # output ROOT path matches RecoilJets_Condor.sh naming:
          #   destBase/run8/RecoilJets_<ds>_<chunkTag>.root
          printf "%s/%s/RecoilJets_%s_%s.root\n", dest, run, ds, n
        }
      ' | sort -u > "$SKIP_FILE"
    ) || true
  fi

  SKIP_BUILT=1

  local nskip
  nskip=$(wc -l < "$SKIP_FILE" | awk '{print $1}')
  local want2="isPP"
  [[ "$TAG" == "auau" ]] && want2="isAuAu"

  if (( nskip > 0 )); then
    printf "${YEL}⚠${RST} [skiplist] %d active %s outputs will be excluded (from condor_q)\n" "$nskip" "$want2" >&2
    (( SKIP_TRACE )) && { printf "${DIM:-}${YEL}⚠${RST} [skiplist] first 5:\n" >&2; head -n 5 "$SKIP_FILE" >&2; }
  else
    printf "${GRN}✔${RST} [skiplist] No active %s RecoilJets jobs found in condor_q\n" "$want2" >&2
  fi
}

# Build a per-run list of files (sorted) to be merged
# NEW behavior:
#   - Excludes any outputs whose producing RecoilJets_Condor.sh job is still in condor_q.
#   - Supports DRYRUN=1 / SKIP_TRACE=1 summaries without polluting stdout (list path only).
make_run_list() {
  local run8="$1"
  local srcdir="${RUN_BASE}/${run8}"
  local list="${TMP_DIR}/run_${run8}.txt"
  local all="${list}.00_all"

  build_active_skiplist

  # 1) All ROOTs on disk for this run
  find "$srcdir" -maxdepth 1 -type f -name "*.root" | sort -V > "$all"
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
# SIM mode:
#   ./mergeRecoilJets.sh isSim firstRound  [groupSize N] [SAMPLE=run28_photonjet10]
#   ./mergeRecoilJets.sh isSim secondRound [condor]      [SAMPLE=run28_photonjet10]
#   - firstRound  : Condor hadd of *existing sim outputs* in groups of groupSize
#   - secondRound : merge firstRound partials into ONE final file
# Output directory (as requested):
#   /sphenix/u/patsfan753/scratch/thesisAnalysis/output/<simTag>
# where simTag defaults to suffix after last "_" in SAMPLE:
#   run28_photonjet10 -> photonjet10
# ============================================================
if [[ "${1}" =~ ^(isSim|sim|SIM)$ ]]; then
  SIM_ACTION="${2:-}"
  shift 2

  # Defaults (overrideable)
  SIM_SAMPLE="run28_photonjet10"
  SIM_GROUP_SIZE="200"          # number of ROOT files per firstRound hadd
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

  # Input sim outputs live here (your RecoilJets sim job output directory)
  SIM_INPUT_BASE="/sphenix/tg/tg01/bulk/jbennett/thesisAna"
  SIM_INPUT_DIR="${SIM_INPUT_BASE}/${SIM_SAMPLE}"

  [[ -d "$SIM_INPUT_DIR" ]] || { err "SIM input directory not found: $SIM_INPUT_DIR"; exit 20; }

  # Output tag + dir (your requested: output/photonjet10)
  SIM_TAG="${SIM_SAMPLE##*_}"
  DEST_DIR="${OUT_BASE}/${SIM_TAG}"
  mkdir -p "$DEST_DIR" "$LOG_DIR" "$OUT_DIR" "$ERR_DIR" "$TMP_DIR"

  say "Dataset : ${BOLD}isSim${RST}"
  say "Sample  : ${BOLD}${SIM_SAMPLE}${RST}"
  say "Input   : ${SIM_INPUT_DIR}"
  say "Output  : ${DEST_DIR}"
  say "groupSize (merge) : ${SIM_GROUP_SIZE}"
  echo

  # Collect the existing sim ROOT outputs produced by your RecoilJets jobs
  mapfile -t SIM_INPUTS < <(find "$SIM_INPUT_DIR" -maxdepth 1 -type f -name "*.root" | sort -V || true)
  if (( ${#SIM_INPUTS[@]} == 0 )); then
    err "No *.root files found in: $SIM_INPUT_DIR"
    exit 21
  fi

  SIM_PARTIAL_PREFIX="chunkMerge_${SIM_TAG}_grp"
  SIM_FINAL="${DEST_DIR}/${FINAL_PREFIX}_${SIM_TAG}_ALL.root"

  if [[ "$SIM_ACTION" == "firstRound" ]]; then

    if $SIM_FIRSTROUND_LOCAL; then
      say "SIM firstRound (LOCAL): ${#SIM_INPUTS[@]} inputs -> grouped hadd on this node"
    else
      need_cmd condor_submit
      say "SIM firstRound: ${#SIM_INPUTS[@]} inputs -> grouped hadd jobs on Condor"
    fi

    # Clean previous sim partials + final (ONLY sim-tagged files)
    find "$DEST_DIR" -maxdepth 1 -type f \
      \( -name "${SIM_PARTIAL_PREFIX}*.root" -o -name "${FINAL_PREFIX}_${SIM_TAG}_ALL.root" \) -delete || true

    total="${#SIM_INPUTS[@]}"
    grp=0

    if $SIM_FIRSTROUND_LOCAL; then
      # LOCAL: build listfiles and run hadd sequentially (same outputs as Condor mode)
      need_cmd hadd

      for ((i=0; i<total; i+=SIM_GROUP_SIZE)); do
        (( grp+=1 ))
        grpTag="$(printf "%03d" "$grp")"

        listfile="${TMP_DIR}/sim_${SIM_TAG}_grp${grpTag}.txt"
        : > "$listfile"

        for ((j=i; j<i+SIM_GROUP_SIZE && j<total; j++)); do
          printf "%s\n" "${SIM_INPUTS[$j]}" >> "$listfile"
        done

        out="${DEST_DIR}/${SIM_PARTIAL_PREFIX}${grpTag}.root"
        say "[LOCAL firstRound] grp=${grpTag} inputs=$(wc -l < "$listfile") -> $(basename "$out")"
        hadd -v 3 -f "$out" @"$listfile"
      done

      say "LOCAL firstRound complete. Partials are under: ${DEST_DIR}"
      exit 0
    fi

    # CONDOR: original behavior (submit one hadd job per group)
    emit_hadd_wrapper "$CONDOR_EXEC"

    SUB="${TMP_DIR}/recoil_sim_${SIM_TAG}_firstRound.sub"
    rm -f "$SUB"
    cat > "$SUB" <<EOT
universe   = vanilla
executable = $CONDOR_EXEC
output     = $OUT_DIR/recoil.sim.${SIM_TAG}.\$(Cluster).\$(Process).out
error      = $ERR_DIR/recoil.sim.${SIM_TAG}.\$(Cluster).\$(Process).err
log        = $LOG_DIR/recoil.sim.${SIM_TAG}.\$(Cluster).\$(Process).log
request_memory = 4GB
getenv = True
should_transfer_files = NO
stream_output = True
stream_error  = True
EOT

    # Build listfiles + one Condor job per group
    for ((i=0; i<total; i+=SIM_GROUP_SIZE)); do
      (( grp+=1 ))
      grpTag="$(printf "%03d" "$grp")"

      listfile="${TMP_DIR}/sim_${SIM_TAG}_grp${grpTag}.txt"
      : > "$listfile"

      for ((j=i; j<i+SIM_GROUP_SIZE && j<total; j++)); do
        printf "%s\n" "${SIM_INPUTS[$j]}" >> "$listfile"
      done

      out="${DEST_DIR}/${SIM_PARTIAL_PREFIX}${grpTag}.root"
      printf 'arguments = %s %s\nqueue\n\n' "$listfile" "$out" >> "$SUB"
    done

    say "Submitting ${BOLD}${grp}${RST} firstRound Condor merge jobs → $(basename "$SUB")"
    condor_submit "$SUB"
    say "FirstRound submitted. Partials will appear under: ${DEST_DIR}"
    exit 0

  elif [[ "$SIM_ACTION" == "secondRound" ]]; then
    # Merge the firstRound partials into ONE final file
    mapfile -t partials < <(ls -1 "${DEST_DIR}/${SIM_PARTIAL_PREFIX}"*.root 2>/dev/null | sort -V || true)
    if (( ${#partials[@]} == 0 )); then
      err "No firstRound partials found in ${DEST_DIR} (expected ${SIM_PARTIAL_PREFIX}*.root). Run: ./mergeRecoilJets.sh isSim firstRound"
      exit 22
    fi

    LIST="${TMP_DIR}/sim_${SIM_TAG}_partialList.txt"
    printf "%s\n" "${partials[@]}" > "$LIST"

    say "SIM secondRound: ${#partials[@]} partials -> ${SIM_FINAL}"

    if $SIM_PREFER_CONDOR; then
      need_cmd condor_submit
      emit_hadd_wrapper "$CONDOR_EXEC"

      SUB="${TMP_DIR}/recoil_sim_${SIM_TAG}_secondRound.sub"
      rm -f "$SUB"
      cat > "$SUB" <<EOT
universe   = vanilla
executable = $CONDOR_EXEC
output     = $OUT_DIR/recoil.sim.${SIM_TAG}.final.\$(Cluster).\$(Process).out
error      = $ERR_DIR/recoil.sim.${SIM_TAG}.final.\$(Cluster).\$(Process).err
log        = $LOG_DIR/recoil.sim.${SIM_TAG}.final.\$(Cluster).\$(Process).log
request_memory = 6GB
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
      say "Running secondRound final merge locally (ROOT hadd)…"
      hadd -v 3 -f "$SIM_FINAL" @"$LIST"
      say "Created ${SIM_FINAL}"
    fi

    exit 0
  else
    err "Unknown isSim action '${SIM_ACTION}'. Allowed: firstRound | secondRound"
    exit 2
  fi
fi

# ============================================================
# Existing DATA behavior (unchanged): condor/addChunks/checkFileOutput
# ============================================================
MODE="$1"
DATASET_REQ="$2"
SUBMODE="${3:-}"

[[ "$MODE" =~ ^(condor|addChunks|checkFileOutput)$ ]] || usage
resolve_dataset "$DATASET_REQ"

say "Dataset: ${BOLD}${TAG}${RST}"
say "Run base: ${RUN_BASE}"
say "Output dir: ${DEST_DIR}"
echo

# ---------- Modes ----------
if [[ "$MODE" == "checkFileOutput" ]]; then
  say "Scanning available run directories under ${RUN_BASE}"
  discover_runs "$RUN_BASE"
  printf "Runs with *.root present: %s\n" "${#RUNS[@]}"
  if ((${#RUNS[@]}==0)); then
    warn "No runs found with *.root files."
    exit 0
  fi

  # Build the active-job skiplist once (read-only condor_q)
  build_active_skiplist

  say "Per-run summary (TOTAL vs BUSY in condor_q vs BUSY present on disk vs ELIGIBLE)"
  printf "  %-10s  %8s  %8s  %12s  %10s\n" "run" "total" "busyInQ" "busyPresent" "eligible"
  printf "  %-10s  %8s  %8s  %12s  %10s\n" "----------" "--------" "--------" "------------" "----------"

  eligibleRuns=0
  for r in "${RUNS[@]}"; do
    srcdir="${RUN_BASE}/${r}"

    # list all present ROOTs (absolute paths)
    all="${TMP_DIR}/check_${r}.all.txt"
    find "$srcdir" -maxdepth 1 -type f -name "*.root" | sort -V > "$all"
    total=$(wc -l < "$all" | awk '{print $1}')

    busyInQ=0
    busyPresent=0
    if [[ -s "$SKIP_FILE" ]]; then
      busyInQ=$(awk -v pre="${RUN_BASE}/${r}/" 'index($0,pre)==1{c++} END{print c+0}' "$SKIP_FILE")
      busyPresent=$(grep -F -x -f "$SKIP_FILE" "$all" 2>/dev/null | wc -l | awk '{print $1}')
    fi

    eligible=$(( total - busyPresent ))
    (( eligible > 0 )) && (( eligibleRuns += 1 ))

    printf "  %-10s  %8d  %8d  %12d  %10d\n" "$r" "$total" "$busyInQ" "$busyPresent" "$eligible"
    rm -f "$all" 2>/dev/null || true
  done

  say "Eligible runs (eligible>0): ${eligibleRuns} / ${#RUNS[@]}"
  say "NOTE: This mode only REPORTS. It does not delete files, submit jobs, or run hadd."
  exit 0
fi

if [[ "$MODE" == "condor" ]]; then
  if (( DRYRUN )); then
    warn "DRYRUN=1 → planning only (NO deletes, NO condor_submit, NO hadd). Reading condor_q is safe."
  else
    need_cmd condor_submit
  fi

  say "Preparing per-run Condor merges (one job per run)"

  SEGMENT_MODE=0
  SEGMENTS=()

  # ------------------------------------------------------------
  # RUN SELECTION
  #
  # Default (backward-compatible):
  #   - discover_runs(): scan RUN_BASE/<run8> dirs that already have *.root
  #
  # New (splitGoldenRunList aware):
  #   - If SUBMODE is one of:
  #       fromSplitRunList | fromSplit | rounds | round | segments | seg
  #     then we read run8 values from:
  #       ${ROUND_DIR}/goldenRuns_${TAG}_segmentK.txt
  #     using the segment indices provided after SUBMODE.
  #
  # This supports your “round … currentNode” token style because we just
  # extract integers from the remaining args.
  #
  # Example:
  #   DRYRUN=1 SKIP_TRACE=1 ./mergeRecoilJets.sh condor pp fromSplitRunList round 1 round 3 currentNode
  # ------------------------------------------------------------
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
      err "No runs found with *.root files under ${RUN_BASE}"
      exit 5
    fi

    # Apply legacy submodes
    case "$SUBMODE" in
      test)        RUNS=( "${RUNS[0]}" ) ;;
      firstHalf)   RUNS=( "${RUNS[@]:0:$(( (${#RUNS[@]}+1)/2 ))}" ) ;;
      "" )         ;;
      * )          err "Unknown submode '$SUBMODE' (allowed: test, firstHalf, fromSplitRunList/rounds/segments)"; exit 6 ;;
    esac
  fi

  # Build the skiplist once up-front (read-only condor_q)
  build_active_skiplist

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
      out="${DEST_DIR}/${PARTIAL_PREFIX}_${r}.root"
      say "  would merge run ${r}: inputs=${nfiles} -> ${out}"
      ((planned+=1))
    done
    say "DRYRUN summary: planned runs=${planned}. (No files deleted; no jobs submitted.)"
    exit 0
  fi

  # Normal mode:
  #   - Legacy behavior (auto-discovery): wipe ALL partials in DEST_DIR (original behavior)
  #   - Segment-restricted behavior     : DO NOT wipe all partials (safe for multi-schedd workflow)
  if (( SEGMENT_MODE )); then
    say "Segment-restricted stage-1: ${BOLD}NOT deleting${RST} existing partials in ${DEST_DIR}"
    say "  (This avoids clobbering partials produced from another schedd/node.)"
    say "  Any run we submit will overwrite its own partial: ${PARTIAL_PREFIX}_<run8>.root"
    # Optional: remove stale ALL file so it's obvious it must be rebuilt later
    rm -f "${DEST_DIR}/${FINAL_PREFIX}_${TAG}_ALL.root" 2>/dev/null || true
  else
    say "Cleaning old partials/final in ${DEST_DIR}"
    find "$DEST_DIR" -maxdepth 1 -type f \
         \( -name "${PARTIAL_PREFIX}_*.root" -o -name "${FINAL_PREFIX}_${TAG}_ALL.root" \) -delete || true
  fi

  emit_hadd_wrapper "$CONDOR_EXEC"

  SUB="${TMP_DIR}/recoil_partials_${TAG}.sub"
  rm -f "$SUB"
  cat > "$SUB" <<EOT
universe   = vanilla
executable = $CONDOR_EXEC
output     = $OUT_DIR/recoil.\$(Cluster).\$(Process).out
error      = $ERR_DIR/recoil.\$(Cluster).\$(Process).err
log        = $LOG_DIR/recoil.\$(Cluster).\$(Process).log
request_memory = 2GB
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

    out="${DEST_DIR}/${PARTIAL_PREFIX}_${r}.root"
    nfiles=$(wc -l < "$listfile" | awk '{print $1}')
    (( SKIP_TRACE )) && say "Queueing run ${r}: inputs=${nfiles} -> $(basename "$out")"

    printf 'arguments = %s %s\nqueue\n\n' "$listfile" "$out" >> "$SUB"
    ((queued+=1))
  done

  (( SKIP_TRACE )) && say "Stage-1 queue summary: queued=${queued}, skipped=${skipped}"

  if (( queued == 0 )); then
    err "No Condor jobs to submit (no non-empty eligible run lists)."
    exit 7
  fi

  say "Submitting ${BOLD}${queued}${RST} Condor merge jobs → $(basename "$SUB")"
  condor_submit "$SUB"
  say "Stage-1 submitted. Partial outputs will appear under: ${DEST_DIR}"
  exit 0
fi

if [[ "$MODE" == "addChunks" ]]; then
  # Final merge, prefer local unless "condor" explicitly provided
  PREFER_CONDOR=false
  [[ "${SUBMODE:-}" == "condor" ]] && PREFER_CONDOR=true

  # Collect partials
  mapfile -t partials < <(ls -1 "${DEST_DIR}/${PARTIAL_PREFIX}_"*.root 2>/dev/null | sort -V || true)
  if (( ${#partials[@]} == 0 )); then
    err "No partials found in ${DEST_DIR} (expected ${PARTIAL_PREFIX}_*.root)."
    exit 8
  fi

  LIST="${TMP_DIR}/partialList_${TAG}.txt"
  printf "%s\n" "${partials[@]}" > "$LIST"
  FINAL="${DEST_DIR}/${FINAL_PREFIX}_${TAG}_ALL.root"

  say "Final merge target: ${FINAL}"
  say "Inputs: ${#partials[@]} partials"

  if $PREFER_CONDOR; then
    need_cmd condor_submit
    emit_hadd_wrapper "$CONDOR_EXEC"

    SUB="${TMP_DIR}/recoil_final_${TAG}.sub"
    rm -f "$SUB"
    cat > "$SUB" <<EOT
universe   = vanilla
executable = $CONDOR_EXEC
output     = $OUT_DIR/recoil.final.\$(Cluster).\$(Process).out
error      = $ERR_DIR/recoil.final.\$(Cluster).\$(Process).err
log        = $LOG_DIR/recoil.final.\$(Cluster).\$(Process).log
request_memory = 3GB
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
    say "Created ${FINAL}"
  fi
  exit 0
fi

# Fallback (should not reach here)
usage

