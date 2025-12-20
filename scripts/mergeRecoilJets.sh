#!/usr/bin/env bash
###############################################################################
# mergeRecoilJets.sh
#
# Data-only hierarchical merge, tailored to RecoilJets_Condor_submit.sh outputs.
#
# STAGE 1 (CONDOR, one job per run):
#   • Find per-run data outputs produced by RecoilJets:
#       PP   : /sphenix/tg/tg01/bulk/jbennett/thesisAna/pp/<run8>/*.root
#       AuAu : /sphenix/tg/tg01/bulk/jbennett/thesisAna/auau/<run8>/*.root
#   • Build a list file per run (sorted).
#   • Submit exactly one hadd job per run → partial:
#       /sphenix/u/patsfan753/scratch/thesisAnalysis/output/<tag>/chunkMerge_run_<run8>.root
#
# STAGE 2 (ADDCHUNKS, local or Condor):
#   • Merge all partials chunkMerge_run_*.root in that <tag> output dir into:
#       /sphenix/u/patsfan753/scratch/thesisAnalysis/output/<tag>/RecoilJets_<tag>_ALL.root
#
# USAGE
#   ./mergeRecoilJets.sh condor <pp|auau> [test|firstHalf]
#   ./mergeRecoilJets.sh addChunks <pp|auau> [condor]
#   ./mergeRecoilJets.sh checkFileOutput <pp|auau>
#
# EXAMPLES
#   ./mergeRecoilJets.sh condor pp
#   ./mergeRecoilJets.sh condor auau firstHalf
#   ./mergeRecoilJets.sh addChunks pp
#   ./mergeRecoilJets.sh addChunks auau condor
#
# NOTES
#   • We intentionally use one Condor job per run (simple, robust).
#   • Inputs per run are all *.root files in that run’s directory (maxdepth=1).
#   • Sorting uses -V (version sort) to keep natural order of chunk indices.
###############################################################################
set -euo pipefail

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

# ---------- Logs and temp ----------
LOG_DIR="/sphenix/u/patsfan753/scratch/thesisAnalysis/log"
OUT_DIR="/sphenix/u/patsfan753/scratch/thesisAnalysis/stdout"
ERR_DIR="/sphenix/u/patsfan753/scratch/thesisAnalysis/error"
TMP_DIR="/sphenix/u/patsfan753/scratch/thesisAnalysis/tmp_recoil_merge"
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
  DEST_DIR="${OUT_BASE}/${TAG}"         # where partials and final live
  mkdir -p "$DEST_DIR" "$LOG_DIR" "$OUT_DIR" "$ERR_DIR" "$TMP_DIR"
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

# Build a per-run list of files (sorted) to be merged
make_run_list() {
  local run8="$1"
  local srcdir="${RUN_BASE}/${run8}"
  local list="${TMP_DIR}/run_${run8}.txt"
  # Only top-level ROOT files; version-sorted
  find "$srcdir" -maxdepth 1 -type f -name "*.root" | sort -V > "$list"
  [[ -s "$list" ]] || return 1
  echo "$list"
}

# ---------- Parse CLI ----------
[[ $# -ge 2 ]] || usage
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
  # Show a succinct summary per run
  for r in "${RUNS[@]}"; do
    n=$(find "${RUN_BASE}/${r}" -maxdepth 1 -type f -name "*.root" | wc -l | awk '{print $1}')
    printf "  run %-8s : %5d files\n" "$r" "$n"
  done
  exit 0
fi

if [[ "$MODE" == "condor" ]]; then
  need_cmd condor_submit

  say "Preparing per-run Condor merges (one job per run)"
  discover_runs "$RUN_BASE"
  if ((${#RUNS[@]}==0)); then
    err "No runs found with *.root files under ${RUN_BASE}"
    exit 5
  fi

  # Apply submodes
  case "$SUBMODE" in
    test)        RUNS=( "${RUNS[0]}" ) ;;
    firstHalf)   RUNS=( "${RUNS[@]:0:$(( (${#RUNS[@]}+1)/2 ))}" ) ;;
    "" )         ;;
    * )          err "Unknown submode '$SUBMODE' (allowed: test, firstHalf)"; exit 6 ;;
  esac

  # Clean slate policy for this dataset's output dir
  say "Cleaning old partials/final in ${DEST_DIR}"
  find "$DEST_DIR" -maxdepth 1 -type f \
       \( -name "${PARTIAL_PREFIX}_*.root" -o -name "${FINAL_PREFIX}_${TAG}_ALL.root" \) -delete || true

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
  for r in "${RUNS[@]}"; do
    listfile="$(make_run_list "$r" || true)"
    if [[ -z "${listfile:-}" || ! -s "$listfile" ]]; then
      warn "Run $r has no files (or list build failed); skipping"
      continue
    fi
    out="${DEST_DIR}/${PARTIAL_PREFIX}_${r}.root"
    printf 'arguments = %s %s\nqueue\n\n' "$listfile" "$out" >> "$SUB"
    ((queued+=1))
  done

  if (( queued == 0 )); then
    err "No Condor jobs to submit (no non-empty run lists)."
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

