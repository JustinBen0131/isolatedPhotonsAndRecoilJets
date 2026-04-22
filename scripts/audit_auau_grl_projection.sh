#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'
shopt -s nullglob

###############################################################################
# audit_auau_grl_projection.sh
#
# Purpose
#   Two-pass Au+Au GRL projection audit for DST_CALOFITTING.
#
# Core idea
#   This script separates the expensive per-file event counting from the final
#   audit/summary logic:
#
#     firstPass
#       • reads the canonical Au+Au GRL
#       • builds the canonical expected per-run DST_CALOFITTING lists with
#         CreateDstList.pl
#       • compares those against the current project lists in:
#           /sphenix/u/patsfan753/scratch/thesisAnalysis/dst_lists_auau
#       • collects the UNIQUE referenced DST ROOT filenames
#       • splits that unique-file inventory into chunk list files
#       • prepares a ROOT counting macro + wrapper
#       • either:
#           - submits Condor counting jobs for all chunks, or
#           - runs one LOCAL validation test on the first available chunk/file
#       • writes a persistent first-pass manifest in the CURRENT directory
#
#     secondPass
#       • loads the first-pass manifest
#       • reads all completed chunk result TSVs from the persisted work area
#       • aggregates file-level event counts
#       • performs the full run-by-run audit:
#           expected lists vs current lists vs counted events
#       • prints the final summary tables
#
# Supported commands
#
#   1) Full Condor first pass
#        ./audit_auau_grl_projection.sh firstPass
#
#      What it does
#        • creates a persisted work area under:
#            /sphenix/u/patsfan753/scratch/thesisAnalysis/condor_snapshots/
#              audit_auau_<timestamp>/
#        • writes chunk lists under:
#            .../chunks/
#        • writes per-chunk result TSV targets under:
#            .../results/
#        • writes Condor submit files under:
#            /sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/
#        • writes Condor logs/stdout/stderr under:
#            /sphenix/u/patsfan753/scratch/thesisAnalysis/log/
#            /sphenix/u/patsfan753/scratch/thesisAnalysis/stdout/
#            /sphenix/u/patsfan753/scratch/thesisAnalysis/error/
#        • writes a manifest in the CURRENT directory:
#            ./audit_auau_firstPass_<timestamp>.txt
#
#   2) Local first-pass validation
#        ./audit_auau_grl_projection.sh firstPass LOCAL
#
#      What it does
#        • performs all first-pass planning exactly as above
#        • does NOT submit Condor jobs
#        • runs the counting wrapper locally on the first available chunk
#        • prints terminal diagnostics so you can verify the parallel logic
#        • still writes a manifest in the CURRENT directory
#
#   3) Final aggregation / audit
#        ./audit_auau_grl_projection.sh secondPass
#
#      What it does
#        • finds the most recent manifest in the CURRENT directory
#        • loads the persisted work area from that manifest
#        • reads completed chunk TSV results
#        • performs the final audit and prints the summary
#
#   4) Final aggregation with an explicit manifest
#        ./audit_auau_grl_projection.sh secondPass --manifest ./audit_auau_firstPass_<timestamp>.txt
#
# Important behavior
#   • This script does NOT modify your real dst_lists_auau directory.
#   • Expected lists produced by CreateDstList.pl are written only inside the
#     persisted first-pass work area.
#   • The expensive per-file event counting is done chunk-by-chunk, suitable
#     for Condor parallelization.
#   • The default cap is:
#         MAX_CONDOR_JOBS=15000
#     and chunk size is chosen automatically from the unique-file count.
#   • secondPass is fast because it only reads completed TSV results rather
#     than reopening all DST files itself.
#
# Useful options
#   --manifest <path>          Use an explicit first-pass manifest in secondPass
#   --grl <path>               Override the GRL
#   --lists-dir <path>         Override the current dst list directory
#   --dataset <token>          Override dataset token
#   --tag <token>              Override CDB tag token
#   --prefix <token>           Override DST prefix
#   --max-condor-jobs <N>      Maximum number of chunk jobs to create
#   --run-progress-every <N>   Progress cadence for run-level loops
#   --root-progress-every <N>  Progress cadence inside the ROOT counter macro
#   --summary-only             Suppress the full per-run table in secondPass
#   -q / --quiet               Reduce chatter
#   -v / --verbose             Increase chatter
#
# Typical workflow
#
#   Pass 1: submit chunk-count jobs
#     ./audit_auau_grl_projection.sh firstPass
#
#   Optional validation instead of submit
#     ./audit_auau_grl_projection.sh firstPass LOCAL
#
#   Pass 2: after jobs finish, aggregate + audit
#     ./audit_auau_grl_projection.sh secondPass
#
###############################################################################

########################################
# Defaults
########################################
BASE="/sphenix/u/patsfan753/scratch/thesisAnalysis"
GRL="${BASE}/GRLs_tanner/run3auau_new_newcdbtag_v008_dst_calofitting_grl.list"
LIST_DIR="${BASE}/dst_lists_auau"

DATASET="run3auau"
TAG="new_newcdbtag_v008"
PREFIX="DST_CALOFITTING"

VERBOSE=1
RUN_PROGRESS_EVERY=100
ROOT_PROGRESS_EVERY=250
SHOW_ALL_RUN_TABLES=1

MAX_CONDOR_JOBS=15000
ACTION=""
LOCAL_MODE=0
MANIFEST_PATH=""
WORK_STAMP=""
TMP_IS_EPHEMERAL=0

CONDOR_SUB_DIR="${BASE}/condor_sub"
CONDOR_LOG_DIR="${BASE}/log"
CONDOR_OUT_DIR="${BASE}/stdout"
CONDOR_ERR_DIR="${BASE}/error"
CONDOR_SNAPSHOT_ROOT="${BASE}/condor_snapshots"

########################################
# Globals set later
########################################
PREFIX_STEM=""
TMP_ROOT=""
TMP_EXPECTED_DIR=""
TMP_EXPECTED_LOG=""
TMP_UNIQUE_FILES=""
TMP_ROOT_MACRO=""
TMP_ROOT_RESULTS=""
TMP_CHUNK_DIR=""
TMP_RESULTS_DIR=""
TMP_WRAPPER=""
CHUNK_SIZE=0
CHUNK_COUNT=0

########################################
# ANSI styles
########################################
if [[ -t 1 ]]; then
  BOLD=$'\033[1m'
  DIM=$'\033[2m'
  RED=$'\033[31m'
  GREEN=$'\033[32m'
  YELLOW=$'\033[33m'
  BLUE=$'\033[34m'
  MAGENTA=$'\033[35m'
  CYAN=$'\033[36m'
  RESET=$'\033[0m'
else
  BOLD=""
  DIM=""
  RED=""
  GREEN=""
  YELLOW=""
  BLUE=""
  MAGENTA=""
  CYAN=""
  RESET=""
fi

########################################
# Cleanup / error handling
########################################
cleanup() {
  if [[ "${TMP_IS_EPHEMERAL:-0}" -eq 1 && -n "${TMP_ROOT:-}" && -d "${TMP_ROOT:-}" ]]; then
    rm -rf "${TMP_ROOT}"
  fi
}

on_err() {
  local line="$1"
  local cmd="$2"
  local rc="$3"
  printf "\n%s[FATAL]%s line %s | rc=%s | command: %s\n" "${RED}${BOLD}" "${RESET}" "$line" "$rc" "$cmd" >&2
  exit "$rc"
}

trap 'on_err "$LINENO" "$BASH_COMMAND" "$?"' ERR
trap cleanup EXIT INT TERM

########################################
# Messaging helpers
########################################
usage() {
  cat <<EOF
Usage:
  $(basename "$0") firstPass [LOCAL] [options]
  $(basename "$0") secondPass [options]

Options:
  --manifest <path>            Override the first-pass manifest used by secondPass
  --grl <path>                 Override GRL path
  --lists-dir <path>           Override current list directory
  --dataset <token>            Default: ${DATASET}
  --tag <token>                Default: ${TAG}
  --prefix <token>             Default: ${PREFIX}
  --max-condor-jobs <N>        Default: ${MAX_CONDOR_JOBS}
  --run-progress-every <N>     Default: ${RUN_PROGRESS_EVERY}
  --root-progress-every <N>    Default: ${ROOT_PROGRESS_EVERY}
  --summary-only               Suppress full run tables in secondPass
  -q, --quiet                  Less chatter
  -v, --verbose                More chatter
  -h, --help                   Show this help

Examples:
  $(basename "$0") firstPass
  $(basename "$0") firstPass LOCAL
  $(basename "$0") secondPass
  $(basename "$0") secondPass --manifest ./audit_auau_firstPass_20260101_120000.txt
EOF
  exit 0
}
say()  { printf "%s➜%s %s\n" "${BLUE}"  "${RESET}" "$*"; }
good() { printf "%s✔%s %s\n" "${GREEN}" "${RESET}" "$*"; }
warn() { printf "%s⚠%s %s\n" "${YELLOW}" "${RESET}" "$*" >&2; }
note() { printf "%s•%s %s\n" "${CYAN}" "${RESET}" "$*"; }

section() {
  local title="$1"
  printf "\n%s%s%s\n" "${BOLD}${MAGENTA}" "$title" "${RESET}"
  printf "%s\n" "$(printf '%*s' 96 '' | tr ' ' '=')"
}

subsection() {
  local title="$1"
  printf "\n%s%s%s\n" "${BOLD}${CYAN}" "$title" "${RESET}"
  printf "%s\n" "$(printf '%*s' 96 '' | tr ' ' '-')"
}

########################################
# Small utilities
########################################
need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    printf "%s[FATAL]%s missing required command: %s\n" "${RED}${BOLD}" "${RESET}" "$1" >&2
    exit 1
  }
}

run8() {
  printf "%08d" "$((10#$1))"
}

num_or_zero() {
  local x="${1:-}"
  if [[ "$x" =~ ^-?[0-9]+([.][0-9]+)?$ ]]; then
    printf '%s' "$x"
  else
    printf '0'
  fi
}

fmt_num() {
  awk -v n="${1:-0}" '
    BEGIN {
      s = sprintf("%.0f", n + 0);
      neg = "";
      if (substr(s,1,1) == "-") {
        neg = "-";
        s = substr(s,2);
      }
      out = "";
      while (length(s) > 3) {
        out = "," substr(s, length(s)-2) out;
        s = substr(s, 1, length(s)-3);
      }
      printf "%s%s%s", neg, s, out;
    }'
}

fmt_pct() {
  awk -v n="${1:-0}" -v d="${2:-0}" '
    BEGIN {
      if (d <= 0) printf "n/a";
      else printf "%.2f%%", 100.0 * n / d;
    }'
}

fmt_m() {
  awk -v n="${1:-0}" 'BEGIN { printf "%.3f", (n / 1.0e6); }'
}

fmt_h() {
  awk -v s="${1:-0}" 'BEGIN { printf "%.2f", (s / 3600.0); }'
}

yn() {
  if [[ "${1:-0}" -ne 0 ]]; then
    printf 'Y'
  else
    printf 'N'
  fi
}

truncate_text() {
  local s="${1:-}"
  local max="${2:-24}"
  if (( ${#s} <= max )); then
    printf '%s' "$s"
  else
    printf '%s...' "${s:0:max-3}"
  fi
}

normalize_name() {
  printf '%s' "${1:-}" | tr '[:upper:]' '[:lower:]' | awk '{$1=$1; print}'
}

append_note() {
  local -n _ref="$1"
  local _note="$2"
  if [[ -z "${_ref}" ]]; then
    _ref="${_note}"
  else
    _ref="${_ref},${_note}"
  fi
}

########################################
# Run-list helpers
########################################
read_runs_from_file() {
  local f="$1"
  mapfile -t RUN_ORDER < <(
    awk '
      {
        sub(/#.*/, "", $0);
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", $0);
        if ($0 ~ /^[0-9]+$/) print $0;
      }' "$f" \
    | sort -n \
    | awk '{printf "%08d\n", $1}' \
    | uniq
  )
}
########################################
# List-file helpers
########################################
recompute_prefix_stem() {
  PREFIX_STEM="$(printf '%s' "${PREFIX#DST_}" | tr '[:upper:]' '[:lower:]')"
}

list_candidates_for_run() {
  local dir="$1"
  local r8="$2"
  printf '%s\n' \
    "${dir}/dst_${PREFIX_STEM}-${r8}.list" \
    "${dir}/${PREFIX}-${r8}.list" \
    "${dir}/${PREFIX}_${DATASET}_${TAG}-${r8}.list"
}

find_best_list_file() {
  local dir="$1"
  local r8="$2"
  local -a existing=()
  local -a nonempty=()
  local c

  while IFS= read -r c; do
    [[ -e "$c" ]] && existing+=("$c")
    [[ -s "$c" ]] && nonempty+=("$c")
  done < <(list_candidates_for_run "$dir" "$r8")

  if (( ${#nonempty[@]} > 1 )) && (( VERBOSE > 0 )); then
    warn "multiple non-empty list candidates for run ${r8} in ${dir}; using ${nonempty[0]}"
  fi

  if (( ${#nonempty[@]} > 0 )); then
    printf '%s\n' "${nonempty[0]}"
    return 0
  fi

  if (( ${#existing[@]} > 1 )) && (( VERBOSE > 0 )); then
    warn "multiple empty/existing list candidates for run ${r8} in ${dir}; using ${existing[0]}"
  fi

  if (( ${#existing[@]} > 0 )); then
    printf '%s\n' "${existing[0]}"
    return 0
  fi

  return 1
}

normalize_list() {
  local f="$1"
  [[ -n "$f" && -f "$f" ]] || return 0
  awk '
    {
      sub(/#.*/, "", $0);
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", $0);
      if ($0 != "") print $0;
    }' "$f" | sort -u
}

read_list_into_array() {
  local f="$1"
  local -n out="$2"
  out=()
  [[ -n "$f" && -f "$f" ]] || return 0
  mapfile -t out < <(normalize_list "$f")
}

########################################
# Expected-list builder
########################################
build_expected_lists() {
  section "Step 1/5 — Build canonical DB-backed expected DST_CALOFITTING lists"

  TMP_EXPECTED_DIR="${TMP_ROOT}/expected_lists"
  TMP_EXPECTED_LOG="${TMP_ROOT}/CreateDstList.log"
  mkdir -p "$TMP_EXPECTED_DIR"

  say "GRL            : ${GRL}"
  say "Temp output dir: ${TMP_EXPECTED_DIR}"
  say "CreateDstList  : tag=${TAG} dataset=${DATASET} prefix=${PREFIX}"

  (
    cd "$TMP_EXPECTED_DIR"
    CreateDstList.pl --tag "$TAG" --dataset "$DATASET" --list "$GRL" "$PREFIX"
  ) >"${TMP_EXPECTED_LOG}" 2>&1

  local n_lists
  n_lists=$(
    {
      compgen -G "${TMP_EXPECTED_DIR}/dst_${PREFIX_STEM}-*.list" || true
      compgen -G "${TMP_EXPECTED_DIR}/${PREFIX}-*.list" || true
      compgen -G "${TMP_EXPECTED_DIR}/${PREFIX}_${DATASET}_${TAG}-*.list" || true
    } | sort -u | awk 'NF{n++} END{print n+0}'
  )

  good "CreateDstList.pl completed"
  note "Temporary expected per-run lists generated: $(fmt_num "$n_lists")"

  if (( VERBOSE > 1 )); then
    subsection "CreateDstList.pl log tail"
    tail -n 30 "${TMP_EXPECTED_LOG}" || true
  fi
}

########################################
# ROOT event-entry counter (no full analysis)
########################################
write_root_counter_macro() {
  TMP_ROOT_MACRO="${TMP_ROOT}/count_dst_entries.C"
  cat > "${TMP_ROOT_MACRO}" <<'EOF'
#include <TFile.h>
#include <TKey.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TROOT.h>
#include <TError.h>
#include <TSystem.h>
#include <frog/FROG.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

namespace dstcount
{
  void trim(std::string& s)
  {
    const char* ws = " \t\r\n";
    const auto b = s.find_first_not_of(ws);
    if (b == std::string::npos) { s.clear(); return; }
    const auto e = s.find_last_not_of(ws);
    s = s.substr(b, e - b + 1);
  }

  void scan_dir(TDirectory* dir, const std::string& base, Long64_t& bestEntries, std::string& bestPath)
  {
    if (!dir) return;
    TIter next(dir->GetListOfKeys());
    while (TKey* key = static_cast<TKey*>(next()))
    {
      TObject* obj = key->ReadObj();
      if (!obj) continue;

      if (obj->InheritsFrom(TTree::Class()))
      {
        TTree* t = dynamic_cast<TTree*>(obj);
        if (t)
        {
          const Long64_t n = t->GetEntriesFast();
          const std::string path = base.empty() ? std::string(t->GetName()) : base + "/" + t->GetName();
          if (n > bestEntries)
          {
            bestEntries = n;
            bestPath = path;
          }
        }
      }
      else if (obj->InheritsFrom(TDirectory::Class()))
      {
        TDirectory* sub = dynamic_cast<TDirectory*>(obj);
        if (sub)
        {
          const std::string subbase = base.empty() ? std::string(sub->GetName()) : base + "/" + sub->GetName();
          scan_dir(sub, subbase, bestEntries, bestPath);
        }
      }

      delete obj;
    }
  }

  Long64_t best_tree_entries(TFile* f, std::string& treeName)
  {
    treeName.clear();
    if (!f) return -1;

    if (TTree* t = dynamic_cast<TTree*>(f->Get("T")))
    {
      treeName = "T";
      return t->GetEntriesFast();
    }

    Long64_t bestEntries = -1;
    std::string bestPath;
    scan_dir(f, "", bestEntries, bestPath);
    treeName = bestPath;
    return bestEntries;
  }
}

void count_dst_entries(const char* input_list,
                       const char* output_file,
                       int progress_every = 250,
                       int verbose = 1)
{
  gROOT->SetBatch(kTRUE);
  gErrorIgnoreLevel = kWarning;

  std::ifstream in(input_list);
  if (!in.is_open())
  {
    std::cerr << "[count_dst_entries] FATAL: cannot open input list: " << input_list << std::endl;
    gSystem->Exit(1);
    return;
  }

  std::ofstream out(output_file);
  if (!out.is_open())
  {
    std::cerr << "[count_dst_entries] FATAL: cannot open output file: " << output_file << std::endl;
    gSystem->Exit(1);
    return;
  }

  std::vector<std::string> files;
  for (std::string line; std::getline(in, line); )
  {
    dstcount::trim(line);
    if (line.empty()) continue;
    if (!line.empty() && line[0] == '#') continue;
    files.push_back(line);
  }

  Long64_t ok_files = 0;
  Long64_t bad_open = 0;
  Long64_t zombies = 0;
  Long64_t no_tree = 0;
  Long64_t total_entries = 0;

  for (std::size_t i = 0; i < files.size(); ++i)
  {
    const std::string& path = files[i];
    if (verbose > 0 && progress_every > 0)
    {
      if ((i + 1) == 1 || ((i + 1) % static_cast<std::size_t>(progress_every) == 0) || (i + 1) == files.size())
      {
        std::cerr << "[count_dst_entries] " << (i + 1) << " / " << files.size() << " : " << path << std::endl;
      }
    }

    std::string openPath = path;
    if (!openPath.empty() && openPath[0] != '/' && gSystem->AccessPathName(openPath.c_str()))
    {
      FROG frog;
      const char* resolved = frog.location(openPath);
      if (resolved && std::string(resolved).size())
      {
        openPath = resolved;
      }
    }

    TFile* f = TFile::Open(openPath.c_str(), "READ");
    if (!f)
    {
      out << "OPENFAIL\t0\t-\t" << path << "\n";
      ++bad_open;
      continue;
    }
    if (f->IsZombie())
    {
      out << "ZOMBIE\t0\t-\t" << path << "\n";
      ++zombies;
      f->Close();
      delete f;
      continue;
    }

    std::string treeName;
    const Long64_t entries = dstcount::best_tree_entries(f, treeName);
    if (entries < 0)
    {
      out << "NOTREE\t0\t-\t" << path << "\n";
      ++no_tree;
    }
    else
    {
      std::replace(treeName.begin(), treeName.end(), '\t', ' ');
      out << "OK\t" << entries << "\t" << treeName << "\t" << path << "\n";
      ++ok_files;
      total_entries += entries;
    }

    f->Close();
    delete f;
  }

  std::cerr << "[count_dst_entries] SUMMARY"
            << " files=" << files.size()
            << " ok=" << ok_files
            << " openfail=" << bad_open
            << " zombie=" << zombies
            << " notree=" << no_tree
            << " totalEntries=" << total_entries
            << std::endl;
}
EOF
}

declare -A FILE_STATUS=()
declare -A FILE_ENTRIES=()
declare -A FILE_TREE=()
declare -A TREE_COUNT=()

ROOT_OK_FILES=0
ROOT_OPENFAIL_FILES=0
ROOT_ZOMBIE_FILES=0
ROOT_NOTREE_FILES=0
ROOT_TOTAL_COUNTED_ENTRIES=0
UNIQUE_FILE_COUNT=0

collect_unique_files() {
  section "Step 2/5 — Collect unique expected/current ROOT files for exact entry counting"

  TMP_UNIQUE_FILES="${TMP_ROOT}/unique_files.txt"
  : > "${TMP_UNIQUE_FILES}"

  local idx=0
  local r8 expected_list current_list
  for r8 in "${RUN_ORDER[@]}"; do
    ((idx+=1))
    expected_list="$(find_best_list_file "$TMP_EXPECTED_DIR" "$r8" || true)"
    current_list="$(find_best_list_file "$LIST_DIR" "$r8" || true)"

    [[ -n "$expected_list" ]] && normalize_list "$expected_list" >> "${TMP_UNIQUE_FILES}"
    [[ -n "$current_list" ]] && normalize_list "$current_list" >> "${TMP_UNIQUE_FILES}"

    if (( VERBOSE > 0 )) && { (( idx == 1 )) || (( idx % RUN_PROGRESS_EVERY == 0 )) || (( idx == ${#RUN_ORDER[@]} )); }; then
      say "unique-file collection progress: ${idx} / ${#RUN_ORDER[@]}"
    fi
  done

  sort -u -o "${TMP_UNIQUE_FILES}" "${TMP_UNIQUE_FILES}"
  UNIQUE_FILE_COUNT="$(awk 'NF{n++} END{print n+0}' "${TMP_UNIQUE_FILES}")"

  good "Unique ROOT files to scan: $(fmt_num "${UNIQUE_FILE_COUNT}")"
}

choose_chunk_size() {
  local total="$1"
  local max_jobs="$2"

  if (( total <= 0 )); then
    echo 1
    return 0
  fi

  if (( max_jobs <= 0 )); then
    echo "$total"
    return 0
  fi

  echo $(( (total + max_jobs - 1) / max_jobs ))
}

write_counter_wrapper() {
  TMP_WRAPPER="${TMP_ROOT}/audit_auau_count_wrapper.sh"
  cat > "${TMP_WRAPPER}" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

chunk_list="${1:?chunk list required}"
output_tsv="${2:?output tsv required}"
macro_path="${3:?root macro path required}"
progress_every="${4:-250}"
verbose="${5:-1}"

export USER="${USER:-$(id -un)}"
export LOGNAME="${LOGNAME:-$USER}"
export HOME="/sphenix/u/${LOGNAME}"

set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n
set -u

root -l -b -q "${macro_path}(\"${chunk_list}\",\"${output_tsv}\",${progress_every},${verbose})"
EOF
  chmod +x "${TMP_WRAPPER}"
}

prepare_persistent_workdir() {
  mkdir -p "${CONDOR_SUB_DIR}" "${CONDOR_LOG_DIR}" "${CONDOR_OUT_DIR}" "${CONDOR_ERR_DIR}" "${CONDOR_SNAPSHOT_ROOT}"

  WORK_STAMP="$(date +%Y%m%d_%H%M%S)"
  TMP_ROOT="${CONDOR_SNAPSHOT_ROOT}/audit_auau_${WORK_STAMP}"
  TMP_IS_EPHEMERAL=0

  mkdir -p "${TMP_ROOT}"
  TMP_EXPECTED_DIR="${TMP_ROOT}/expected_lists"
  TMP_EXPECTED_LOG="${TMP_ROOT}/CreateDstList.log"
  TMP_UNIQUE_FILES="${TMP_ROOT}/unique_files.txt"
  TMP_ROOT_MACRO="${TMP_ROOT}/count_dst_entries.C"
  TMP_ROOT_RESULTS="${TMP_ROOT}/root_results.tsv"
  TMP_CHUNK_DIR="${TMP_ROOT}/chunks"
  TMP_RESULTS_DIR="${TMP_ROOT}/results"

  mkdir -p "${TMP_EXPECTED_DIR}" "${TMP_CHUNK_DIR}" "${TMP_RESULTS_DIR}"
}

build_count_chunks() {
  TMP_CHUNK_DIR="${TMP_ROOT}/chunks"
  TMP_RESULTS_DIR="${TMP_ROOT}/results"
  mkdir -p "${TMP_CHUNK_DIR}" "${TMP_RESULTS_DIR}"

  rm -f "${TMP_CHUNK_DIR}/audit_chunk_"*.list 2>/dev/null || true
  rm -f "${TMP_RESULTS_DIR}/audit_chunk_"*.tsv 2>/dev/null || true

  CHUNK_SIZE="$(choose_chunk_size "${UNIQUE_FILE_COUNT}" "${MAX_CONDOR_JOBS}")"
  (( CHUNK_SIZE < 1 )) && CHUNK_SIZE=1

  split -l "${CHUNK_SIZE}" -d -a 5 "${TMP_UNIQUE_FILES}" "${TMP_CHUNK_DIR}/audit_chunk_"

  CHUNK_COUNT=0
  local raw out
  for raw in "${TMP_CHUNK_DIR}/audit_chunk_"*; do
    [[ -f "$raw" ]] || continue
    out="${raw}.list"
    mv "$raw" "$out"
    ((CHUNK_COUNT+=1))
  done
}

write_firstpass_manifest() {
  local manifest_path="$1"
  local condor_sub_path="$2"

  {
    printf 'MANIFEST_VERSION=%q\n' "1"
    printf 'ACTION=%q\n' "firstPass"
    printf 'BASE=%q\n' "$BASE"
    printf 'GRL=%q\n' "$GRL"
    printf 'LIST_DIR=%q\n' "$LIST_DIR"
    printf 'DATASET=%q\n' "$DATASET"
    printf 'TAG=%q\n' "$TAG"
    printf 'PREFIX=%q\n' "$PREFIX"
    printf 'PREFIX_STEM=%q\n' "$PREFIX_STEM"
    printf 'VERBOSE=%q\n' "$VERBOSE"
    printf 'RUN_PROGRESS_EVERY=%q\n' "$RUN_PROGRESS_EVERY"
    printf 'ROOT_PROGRESS_EVERY=%q\n' "$ROOT_PROGRESS_EVERY"
    printf 'SHOW_ALL_RUN_TABLES=%q\n' "$SHOW_ALL_RUN_TABLES"
    printf 'MAX_CONDOR_JOBS=%q\n' "$MAX_CONDOR_JOBS"
    printf 'WORK_STAMP=%q\n' "$WORK_STAMP"
    printf 'TMP_ROOT=%q\n' "$TMP_ROOT"
    printf 'TMP_EXPECTED_DIR=%q\n' "$TMP_EXPECTED_DIR"
    printf 'TMP_EXPECTED_LOG=%q\n' "$TMP_EXPECTED_LOG"
    printf 'TMP_UNIQUE_FILES=%q\n' "$TMP_UNIQUE_FILES"
    printf 'TMP_ROOT_MACRO=%q\n' "$TMP_ROOT_MACRO"
    printf 'TMP_ROOT_RESULTS=%q\n' "$TMP_ROOT_RESULTS"
    printf 'TMP_CHUNK_DIR=%q\n' "$TMP_CHUNK_DIR"
    printf 'TMP_RESULTS_DIR=%q\n' "$TMP_RESULTS_DIR"
    printf 'TMP_WRAPPER=%q\n' "$TMP_WRAPPER"
    printf 'CHUNK_SIZE=%q\n' "$CHUNK_SIZE"
    printf 'CHUNK_COUNT=%q\n' "$CHUNK_COUNT"
    printf 'UNIQUE_FILE_COUNT=%q\n' "$UNIQUE_FILE_COUNT"
    printf 'RUN_COUNT=%q\n' "${#RUN_ORDER[@]}"
    printf 'LOCAL_MODE=%q\n' "$LOCAL_MODE"
    printf 'CONDOR_SUB=%q\n' "$condor_sub_path"
  } > "${manifest_path}"
}

latest_manifest_in_pwd() {
  ls -1t ./audit_auau_firstPass_*.txt 2>/dev/null | head -n 1
}

load_firstpass_manifest() {
  local manifest_path="$1"
  [[ -f "${manifest_path}" ]] || fatal "manifest not found: ${manifest_path}"
  # shellcheck disable=SC1090
  source "${manifest_path}"

  TMP_IS_EPHEMERAL=0
  TMP_EXPECTED_DIR="${TMP_ROOT}/expected_lists"
  TMP_EXPECTED_LOG="${TMP_ROOT}/CreateDstList.log"
  TMP_UNIQUE_FILES="${TMP_ROOT}/unique_files.txt"
  TMP_ROOT_MACRO="${TMP_ROOT}/count_dst_entries.C"
  TMP_ROOT_RESULTS="${TMP_ROOT}/root_results.tsv"
  TMP_CHUNK_DIR="${TMP_ROOT}/chunks"
  TMP_RESULTS_DIR="${TMP_ROOT}/results"
  TMP_WRAPPER="${TMP_ROOT}/audit_auau_count_wrapper.sh"
}

submit_condor_counter_jobs() {
  need_cmd condor_submit

  local sub="${CONDOR_SUB_DIR}/audit_auau_${WORK_STAMP}.sub"
  cat > "${sub}" <<SUB
universe      = vanilla
executable    = ${TMP_WRAPPER}
initialdir    = ${BASE}
getenv        = True
log           = ${CONDOR_LOG_DIR}/audit_auau.${WORK_STAMP}.\$(Cluster).\$(Process).log
output        = ${CONDOR_OUT_DIR}/audit_auau.${WORK_STAMP}.\$(Cluster).\$(Process).out
error         = ${CONDOR_ERR_DIR}/audit_auau.${WORK_STAMP}.\$(Cluster).\$(Process).err
request_memory= 2000MB
should_transfer_files = NO
stream_output = True
stream_error  = True
SUB

  local chunk_list out_tsv queued=0
  for chunk_list in "${TMP_CHUNK_DIR}/audit_chunk_"*.list; do
    [[ -f "${chunk_list}" ]] || continue
    out_tsv="${TMP_RESULTS_DIR}/$(basename "${chunk_list%.list}").tsv"
    printf 'arguments = %s %s %s %s %s\nqueue\n\n' \
      "${chunk_list}" "${out_tsv}" "${TMP_ROOT_MACRO}" "${ROOT_PROGRESS_EVERY}" "${VERBOSE}" >> "${sub}"
    ((queued+=1))
  done

  (( queued > 0 )) || fatal "no condor jobs were prepared"

  say "Submitting ${queued} Condor counting jobs"
  say "Condor submit file : ${sub}"
  say "Chunk dir          : ${TMP_CHUNK_DIR}"
  say "Results dir        : ${TMP_RESULTS_DIR}"

  condor_submit "${sub}"
  MANIFEST_PATH="${PWD}/audit_auau_firstPass_${WORK_STAMP}.txt"
  write_firstpass_manifest "${MANIFEST_PATH}" "${sub}"
  good "First-pass manifest written: ${MANIFEST_PATH}"
  note "Run second pass later with: ./$(basename "$0") secondPass"
}

run_local_counter_test() {
  local first_chunk
  first_chunk="$(ls -1 "${TMP_CHUNK_DIR}"/audit_chunk_*.list 2>/dev/null | head -n 1 || true)"
  [[ -n "${first_chunk}" ]] || fatal "no chunk list available for LOCAL test"

  local first_file
  first_file="$(head -n 1 "${first_chunk}" 2>/dev/null || true)"
  local local_out="${TMP_RESULTS_DIR}/$(basename "${first_chunk%.list}")_LOCAL.tsv"

  MANIFEST_PATH="${PWD}/audit_auau_firstPass_${WORK_STAMP}.txt"
  write_firstpass_manifest "${MANIFEST_PATH}" "LOCAL_TEST"

  subsection "LOCAL first-pass validation"
  say "Manifest path      : ${MANIFEST_PATH}"
  say "Test chunk list    : ${first_chunk}"
  say "First file in test : ${first_file}"
  say "Output TSV         : ${local_out}"
  note "This runs exactly the same counting wrapper locally on the first available file."

  bash "${TMP_WRAPPER}" "${first_chunk}" "${local_out}" "${TMP_ROOT_MACRO}" "${ROOT_PROGRESS_EVERY}" "${VERBOSE}"

  good "LOCAL first-pass validation completed"
  note "For the full submission, run: ./$(basename "$0") firstPass"
}

collect_count_results() {
  ROOT_OK_FILES=0
  ROOT_OPENFAIL_FILES=0
  ROOT_ZOMBIE_FILES=0
  ROOT_NOTREE_FILES=0
  ROOT_TOTAL_COUNTED_ENTRIES=0
  TREE_COUNT=()
  FILE_STATUS=()
  FILE_ENTRIES=()
  FILE_TREE=()

  local results_found=0
  local tsv status entries tree path tree_key

  for tsv in "${TMP_RESULTS_DIR}"/audit_chunk_*.tsv "${TMP_RESULTS_DIR}"/audit_chunk_*_LOCAL.tsv; do
    [[ -s "${tsv}" ]] || continue
    ((results_found+=1))
    while IFS=$'\t' read -r status entries tree path; do
      [[ -n "${path:-}" ]] || continue
      FILE_STATUS["$path"]="$status"
      FILE_ENTRIES["$path"]="$(num_or_zero "${entries:-0}")"
      FILE_TREE["$path"]="$tree"
    done < "${tsv}"
  done

  local expected_results="${CHUNK_COUNT:-0}"
  if (( LOCAL_MODE == 1 )); then
    note "Collected LOCAL validation result files: $(fmt_num "${results_found}")"
  else
    note "Collected job result files: $(fmt_num "${results_found}") / $(fmt_num "${expected_results}")"
  fi

  local unique_path
  while IFS= read -r unique_path; do
    [[ -n "${unique_path:-}" ]] || continue
    case "${FILE_STATUS[$unique_path]:-UNSCANNED}" in
      OK)
        ((ROOT_OK_FILES+=1))
        ROOT_TOTAL_COUNTED_ENTRIES=$(( ROOT_TOTAL_COUNTED_ENTRIES + ${FILE_ENTRIES[$unique_path]:-0} ))
        tree_key="${FILE_TREE[$unique_path]:-UNKNOWN}"
        TREE_COUNT["$tree_key"]=$(( ${TREE_COUNT["$tree_key"]:-0} + 1 ))
        ;;
      OPENFAIL) ((ROOT_OPENFAIL_FILES+=1)) ;;
      ZOMBIE)   ((ROOT_ZOMBIE_FILES+=1))   ;;
      NOTREE)   ((ROOT_NOTREE_FILES+=1))   ;;
      *)        ;;
    esac
  done < "${TMP_UNIQUE_FILES}"

  good "Collected event-count results"
  note "Files with counts available     : $(fmt_num "${ROOT_OK_FILES}")"
  note "Open failures                  : $(fmt_num "${ROOT_OPENFAIL_FILES}")"
  note "Zombie files                   : $(fmt_num "${ROOT_ZOMBIE_FILES}")"
  note "Files with no TTree            : $(fmt_num "${ROOT_NOTREE_FILES}")"
  note "Total counted entries available: $(fmt_num "${ROOT_TOTAL_COUNTED_ENTRIES}")"
}
########################################
# Per-run audit storage
########################################
declare -A RUN_EXPECTED_LIST=()
declare -A RUN_CURRENT_LIST=()

declare -A RUN_EXPECTED_SEG=()
declare -A RUN_CURRENT_SEG=()
declare -A RUN_PRESENT_SEG=()
declare -A RUN_MISSING_SEG=()
declare -A RUN_EXTRA_SEG=()

declare -A RUN_EXPECTED_EVT=()
declare -A RUN_CURRENT_EVT=()
declare -A RUN_PRESENT_EVT=()
declare -A RUN_MISSING_EVT=()
declare -A RUN_EXTRA_EVT=()

declare -A RUN_STATUS=()
declare -A RUN_NOTES=()

TOT_EXPECTED_SEG=0
TOT_CURRENT_SEG=0
TOT_PRESENT_SEG=0
TOT_MISSING_SEG=0
TOT_EXTRA_SEG=0

TOT_EXPECTED_EVT=0
TOT_CURRENT_EVT=0
TOT_PRESENT_EVT=0
TOT_MISSING_EVT=0
TOT_EXTRA_EVT=0

RUNS_WITH_EXPECTED=0
RUNS_WITH_CURRENT=0
RUNS_WITH_PRESENT=0
RUNS_COMPLETE=0
RUNS_PARTIAL=0
RUNS_ZERO_PRESENT=0
RUNS_WITH_EXTRA=0
RUNS_NO_EXPECTED=0

audit_runs() {
  section "Step 4/5 — Run-by-run audit (DB-backed expectation vs current lists vs collected count results)"

  local idx=0
  local r8 expected_list current_list
  for r8 in "${RUN_ORDER[@]}"; do
    ((idx+=1))

    expected_list="$(find_best_list_file "$TMP_EXPECTED_DIR" "$r8" || true)"
    current_list="$(find_best_list_file "$LIST_DIR" "$r8" || true)"

    RUN_EXPECTED_LIST["$r8"]="$expected_list"
    RUN_CURRENT_LIST["$r8"]="$current_list"

    local -a expected_files=()
    local -a current_files=()
    read_list_into_array "$expected_list" expected_files
    read_list_into_array "$current_list" current_files

    local expected_seg="${#expected_files[@]}"
    local current_seg="${#current_files[@]}"

    RUN_EXPECTED_SEG["$r8"]="$expected_seg"
    RUN_CURRENT_SEG["$r8"]="$current_seg"

    (( TOT_EXPECTED_SEG += expected_seg ))
    (( TOT_CURRENT_SEG  += current_seg  ))

    if (( expected_seg > 0 )); then (( RUNS_WITH_EXPECTED += 1 )); else (( RUNS_NO_EXPECTED += 1 )); fi
    if (( current_seg  > 0 )); then (( RUNS_WITH_CURRENT  += 1 )); fi

    local -A expected_set=()
    local -A current_set=()
    local f
    for f in "${expected_files[@]}"; do expected_set["$f"]=1; done
    for f in "${current_files[@]}"; do current_set["$f"]=1; done

    local present_seg=0
    local missing_seg=0
    local extra_seg=0

    local expected_evt=0
    local current_evt=0
    local present_evt=0
    local extra_evt=0

    local status ent
    for f in "${expected_files[@]}"; do
      status="${FILE_STATUS[$f]:-UNSCANNED}"
      ent="${FILE_ENTRIES[$f]:-0}"
      if [[ "$status" == "OK" ]]; then
        (( expected_evt += ent ))
      fi

      if [[ -n "${current_set[$f]:-}" ]]; then
        (( present_seg += 1 ))
        if [[ "$status" == "OK" ]]; then
          (( present_evt += ent ))
        fi
      else
        (( missing_seg += 1 ))
      fi
    done

    for f in "${current_files[@]}"; do
      status="${FILE_STATUS[$f]:-UNSCANNED}"
      ent="${FILE_ENTRIES[$f]:-0}"
      if [[ "$status" == "OK" ]]; then
        (( current_evt += ent ))
      fi

      if [[ -z "${expected_set[$f]:-}" ]]; then
        (( extra_seg += 1 ))
        if [[ "$status" == "OK" ]]; then
          (( extra_evt += ent ))
        fi
      fi
    done

    local missing_evt=$(( expected_evt - present_evt ))

    RUN_PRESENT_SEG["$r8"]="$present_seg"
    RUN_MISSING_SEG["$r8"]="$missing_seg"
    RUN_EXTRA_SEG["$r8"]="$extra_seg"

    RUN_EXPECTED_EVT["$r8"]="$expected_evt"
    RUN_CURRENT_EVT["$r8"]="$current_evt"
    RUN_PRESENT_EVT["$r8"]="$present_evt"
    RUN_MISSING_EVT["$r8"]="$missing_evt"
    RUN_EXTRA_EVT["$r8"]="$extra_evt"

    (( TOT_PRESENT_SEG += present_seg ))
    (( TOT_MISSING_SEG += missing_seg ))
    (( TOT_EXTRA_SEG   += extra_seg   ))

    (( TOT_EXPECTED_EVT += expected_evt ))
    (( TOT_CURRENT_EVT  += current_evt  ))
    (( TOT_PRESENT_EVT  += present_evt  ))
    (( TOT_MISSING_EVT  += missing_evt  ))
    (( TOT_EXTRA_EVT    += extra_evt    ))

    if (( present_seg > 0 )); then (( RUNS_WITH_PRESENT += 1 )); fi
    if (( expected_seg > 0 && missing_seg == 0 && extra_seg == 0 )); then (( RUNS_COMPLETE += 1 )); fi
    if (( expected_seg > 0 && missing_seg > 0 )); then (( RUNS_PARTIAL += 1 )); fi
    if (( expected_seg > 0 && present_seg == 0 )); then (( RUNS_ZERO_PRESENT += 1 )); fi
    if (( extra_seg > 0 )); then (( RUNS_WITH_EXTRA += 1 )); fi

    local notes=""
    if [[ -z "$expected_list" ]]; then append_note notes "NO_DB_EXPECTED_LIST"; fi
    if [[ -n "$expected_list" && ! -s "$expected_list" ]]; then append_note notes "EMPTY_DB_EXPECTED_LIST"; fi
    if [[ -z "$current_list" ]]; then append_note notes "NO_CURRENT_LIST"; fi
    if [[ -n "$current_list" && ! -s "$current_list" ]]; then append_note notes "EMPTY_CURRENT_LIST"; fi

    local run_status
    if (( expected_seg == 0 && current_seg == 0 )); then
      run_status="NOEXP"
    elif (( expected_seg == 0 && current_seg > 0 )); then
      run_status="ORPHAN"
    elif (( expected_seg > 0 && current_seg == 0 )); then
      run_status="EMPTY"
    elif (( expected_seg > 0 && current_seg > 0 && present_seg == 0 )); then
      run_status="DISJOINT"
    elif (( missing_seg > 0 && extra_seg > 0 )); then
      run_status="MIXED"
    elif (( missing_seg > 0 )); then
      run_status="MISS"
    elif (( extra_seg > 0 )); then
      run_status="EXTRA"
    else
      run_status="OK"
    fi

    RUN_STATUS["$r8"]="$run_status"
    RUN_NOTES["$r8"]="$notes"

    if (( VERBOSE > 0 )) && { (( idx == 1 )) || (( idx % RUN_PROGRESS_EVERY == 0 )) || (( idx == ${#RUN_ORDER[@]} )); }; then
      say "run audit progress: ${idx} / ${#RUN_ORDER[@]}"
    fi
  done

  good "Per-run audit completed"
}
########################################
# Printing helpers / tables
########################################
print_core_summary() {
  section "Step 5/5 — Final terminal summary"

  note "Expected segments are the current DB-backed canonical expectation (CreateDstList.pl in the persisted first-pass work dir)."
  note "Event totals come from collected Condor/local counting job results."
  note "No project DST list files were modified."

  subsection "A) GRL / run coverage"
  printf "  %-42s : %12s\n" "GRL runs" "$(fmt_num "${#RUN_ORDER[@]}")"
  printf "  %-42s : %12s\n" "Runs with DB-backed expected DST segments" "$(fmt_num "${RUNS_WITH_EXPECTED}")"
  printf "  %-42s : %12s\n" "Runs with current non-empty project lists" "$(fmt_num "${RUNS_WITH_CURRENT}")"
  printf "  %-42s : %12s\n" "Runs with any expected segment currently present" "$(fmt_num "${RUNS_WITH_PRESENT}")"
  printf "  %-42s : %12s\n" "Runs with full expected coverage" "$(fmt_num "${RUNS_COMPLETE}")"
  printf "  %-42s : %12s\n" "Runs with partial missing expected coverage" "$(fmt_num "${RUNS_PARTIAL}")"
  printf "  %-42s : %12s\n" "Runs with zero expected segments present" "$(fmt_num "${RUNS_ZERO_PRESENT}")"
  printf "  %-42s : %12s\n" "Runs with extra current-listed segments" "$(fmt_num "${RUNS_WITH_EXTRA}")"
  printf "  %-42s : %12s\n" "Runs with no DB-backed expected DST segments" "$(fmt_num "${RUNS_NO_EXPECTED}")"

  subsection "B) Segment / file availability"
  printf "  %-42s : %12s\n" "Expected segments total (DB-backed)" "$(fmt_num "${TOT_EXPECTED_SEG}")"
  printf "  %-42s : %12s\n" "Present expected segments in current lists" "$(fmt_num "${TOT_PRESENT_SEG}")"
  printf "  %-42s : %12s  (%s)\n" "Present/expected completeness" "$(fmt_num "${TOT_PRESENT_SEG}")" "$(fmt_pct "${TOT_PRESENT_SEG}" "${TOT_EXPECTED_SEG}")"
  printf "  %-42s : %12s\n" "Missing expected segments" "$(fmt_num "${TOT_MISSING_SEG}")"
  printf "  %-42s : %12s\n" "Current listed segments total" "$(fmt_num "${TOT_CURRENT_SEG}")"
  printf "  %-42s : %12s\n" "Extra current-listed segments (not in expected)" "$(fmt_num "${TOT_EXTRA_SEG}")"

  subsection "C) Exact DST event availability"
  printf "  %-42s : %12s\n" "Expected exact DST events" "$(fmt_num "${TOT_EXPECTED_EVT}")"
  printf "  %-42s : %12s\n" "Present expected exact DST events" "$(fmt_num "${TOT_PRESENT_EVT}")"
  printf "  %-42s : %12s  (%s)\n" "Present/expected event completeness" "$(fmt_num "${TOT_PRESENT_EVT}")" "$(fmt_pct "${TOT_PRESENT_EVT}" "${TOT_EXPECTED_EVT}")"
  printf "  %-42s : %12s\n" "Missing expected exact DST events" "$(fmt_num "${TOT_MISSING_EVT}")"
  printf "  %-42s : %12s\n" "Current listed exact DST events" "$(fmt_num "${TOT_CURRENT_EVT}")"
  printf "  %-42s : %12s\n" "Extra current-listed exact DST events" "$(fmt_num "${TOT_EXTRA_EVT}")"

  subsection "D) Entry-count health"
  printf "  %-42s : %12s\n" "Unique referenced ROOT files" "$(fmt_num "${UNIQUE_FILE_COUNT}")"
  printf "  %-42s : %12s\n" "Files with counts available" "$(fmt_num "${ROOT_OK_FILES}")"
  printf "  %-42s : %12s\n" "Open failures" "$(fmt_num "${ROOT_OPENFAIL_FILES}")"
  printf "  %-42s : %12s\n" "Zombie files" "$(fmt_num "${ROOT_ZOMBIE_FILES}")"
  printf "  %-42s : %12s\n" "Files with no tree" "$(fmt_num "${ROOT_NOTREE_FILES}")"
  printf "  %-42s : %12s\n" "Total counted entries available" "$(fmt_num "${ROOT_TOTAL_COUNTED_ENTRIES}")"
}

print_ranked_tables() {
  subsection "E) Ranked high-impact runs"

  printf "  %-8s | %-8s | %-8s | %-8s | %-12s | %-12s | %-10s\n" \
    "Run" "ExpSeg" "PresSeg" "MissSeg" "PresEvt(M)" "MissEvt(M)" "Status"
  printf "  %s-+-%s-+-%s-+-%s-+-%s-+-%s-+-%s\n" \
    "$(printf '%*s' 8 '' | tr ' ' '-')" \
    "$(printf '%*s' 8 '' | tr ' ' '-')" \
    "$(printf '%*s' 8 '' | tr ' ' '-')" \
    "$(printf '%*s' 8 '' | tr ' ' '-')" \
    "$(printf '%*s' 12 '' | tr ' ' '-')" \
    "$(printf '%*s' 12 '' | tr ' ' '-')" \
    "$(printf '%*s' 10 '' | tr ' ' '-')"

  while IFS=$'\t' read -r miss_evt miss_seg r8; do
    [[ -n "${r8:-}" ]] || continue
    printf "  %-8s | %8s | %8s | %8s | %12s | %12s | %-10s\n" \
      "$r8" \
      "$(fmt_num "${RUN_EXPECTED_SEG[$r8]:-0}")" \
      "$(fmt_num "${RUN_PRESENT_SEG[$r8]:-0}")" \
      "$(fmt_num "${RUN_MISSING_SEG[$r8]:-0}")" \
      "$(fmt_m "${RUN_PRESENT_EVT[$r8]:-0}")" \
      "$(fmt_m "${RUN_MISSING_EVT[$r8]:-0}")" \
      "${RUN_STATUS[$r8]:-UNK}"
  done < <(
    for r8 in "${RUN_ORDER[@]}"; do
      printf '%015d\t%015d\t%s\n' "${RUN_MISSING_EVT[$r8]:-0}" "${RUN_MISSING_SEG[$r8]:-0}" "$r8"
    done | sort -t $'\t' -k1,1nr -k2,2nr -k3,3 | head -n 15
  )

  subsection "F) Top runs by currently available exact DST events"
  printf "  %-8s | %-8s | %-12s | %-10s\n" \
    "Run" "CurSeg" "CurEvt(M)" "Status"
  printf "  %s-+-%s-+-%s-+-%s\n" \
    "$(printf '%*s' 8 '' | tr ' ' '-')" \
    "$(printf '%*s' 8 '' | tr ' ' '-')" \
    "$(printf '%*s' 12 '' | tr ' ' '-')" \
    "$(printf '%*s' 10 '' | tr ' ' '-')"

  while IFS=$'\t' read -r cur_evt r8; do
    [[ -n "${r8:-}" ]] || continue
    printf "  %-8s | %8s | %12s | %-10s\n" \
      "$r8" \
      "$(fmt_num "${RUN_CURRENT_SEG[$r8]:-0}")" \
      "$(fmt_m "${RUN_CURRENT_EVT[$r8]:-0}")" \
      "${RUN_STATUS[$r8]:-UNK}"
  done < <(
    for r8 in "${RUN_ORDER[@]}"; do
      printf '%015d\t%s\n' "${RUN_CURRENT_EVT[$r8]:-0}" "$r8"
    done | sort -t $'\t' -k1,1nr -k2,2 | head -n 15
  )
}

print_run_status_table() {
  if (( SHOW_ALL_RUN_TABLES == 0 )); then
    return 0
  fi

  subsection "G) Full per-run audit table"

  printf "  %-8s | %-6s | %-6s | %-6s | %-6s | %-6s | %-10s | %-10s | %-28s\n" \
    "Run" "Exp" "Cur" "Pres" "Miss" "Extra" "PresEvt(M)" "MissEvt(M)" "Notes"
  printf "  %s-+-%s-+-%s-+-%s-+-%s-+-%s-+-%s-+-%s-+-%s\n" \
    "$(printf '%*s' 8 '' | tr ' ' '-')" \
    "$(printf '%*s' 6 '' | tr ' ' '-')" \
    "$(printf '%*s' 6 '' | tr ' ' '-')" \
    "$(printf '%*s' 6 '' | tr ' ' '-')" \
    "$(printf '%*s' 6 '' | tr ' ' '-')" \
    "$(printf '%*s' 6 '' | tr ' ' '-')" \
    "$(printf '%*s' 10 '' | tr ' ' '-')" \
    "$(printf '%*s' 10 '' | tr ' ' '-')" \
    "$(printf '%*s' 28 '' | tr ' ' '-')"

  local r8
  for r8 in "${RUN_ORDER[@]}"; do
    printf "  %-8s | %6s | %6s | %6s | %6s | %6s | %10s | %10s | %-28s\n" \
      "$r8" \
      "$(fmt_num "${RUN_EXPECTED_SEG[$r8]:-0}")" \
      "$(fmt_num "${RUN_CURRENT_SEG[$r8]:-0}")" \
      "$(fmt_num "${RUN_PRESENT_SEG[$r8]:-0}")" \
      "$(fmt_num "${RUN_MISSING_SEG[$r8]:-0}")" \
      "$(fmt_num "${RUN_EXTRA_SEG[$r8]:-0}")" \
      "$(fmt_m "${RUN_PRESENT_EVT[$r8]:-0}")" \
      "$(fmt_m "${RUN_MISSING_EVT[$r8]:-0}")" \
      "$(truncate_text "${RUN_STATUS[$r8]:-UNK}${RUN_NOTES[$r8]:+,${RUN_NOTES[$r8]}}" 28)"
  done
}
fatal() {
  printf "%s[FATAL]%s %s\n" "${RED}${BOLD}" "${RESET}" "$*" >&2
  exit 1
}

parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      firstPass|secondPass)
        ACTION="$1"
        shift
        ;;
      LOCAL)
        LOCAL_MODE=1
        shift
        ;;
      --manifest)
        shift
        [[ $# -gt 0 ]] || fatal "missing value after --manifest"
        MANIFEST_PATH="$1"
        shift
        ;;
      --grl)
        shift
        [[ $# -gt 0 ]] || fatal "missing value after --grl"
        GRL="$1"
        shift
        ;;
      --lists-dir)
        shift
        [[ $# -gt 0 ]] || fatal "missing value after --lists-dir"
        LIST_DIR="$1"
        shift
        ;;
      --dataset)
        shift
        [[ $# -gt 0 ]] || fatal "missing value after --dataset"
        DATASET="$1"
        shift
        ;;
      --tag)
        shift
        [[ $# -gt 0 ]] || fatal "missing value after --tag"
        TAG="$1"
        shift
        ;;
      --prefix)
        shift
        [[ $# -gt 0 ]] || fatal "missing value after --prefix"
        PREFIX="$1"
        shift
        ;;
      --max-condor-jobs)
        shift
        [[ $# -gt 0 ]] || fatal "missing value after --max-condor-jobs"
        MAX_CONDOR_JOBS="$1"
        shift
        ;;
      --run-progress-every)
        shift
        [[ $# -gt 0 ]] || fatal "missing value after --run-progress-every"
        RUN_PROGRESS_EVERY="$1"
        shift
        ;;
      --root-progress-every)
        shift
        [[ $# -gt 0 ]] || fatal "missing value after --root-progress-every"
        ROOT_PROGRESS_EVERY="$1"
        shift
        ;;
      --summary-only)
        SHOW_ALL_RUN_TABLES=0
        shift
        ;;
      -q|--quiet)
        VERBOSE=0
        shift
        ;;
      -v|--verbose)
        ((VERBOSE+=1))
        shift
        ;;
      -h|--help)
        usage
        ;;
      *)
        fatal "unrecognized argument: $1"
        ;;
    esac
  done

  [[ -n "${ACTION}" ]] || ACTION="firstPass"
  [[ "$RUN_PROGRESS_EVERY" =~ ^[0-9]+$ ]] || fatal "--run-progress-every must be an integer"
  [[ "$ROOT_PROGRESS_EVERY" =~ ^[0-9]+$ ]] || fatal "--root-progress-every must be an integer"
  [[ "$MAX_CONDOR_JOBS" =~ ^[0-9]+$ ]] || fatal "--max-condor-jobs must be an integer"
  (( MAX_CONDOR_JOBS > 0 )) || fatal "--max-condor-jobs must be > 0"
}

print_startup_context() {
  section "Au+Au GRL projection audit"

  printf "  %-24s : %s\n" "Script" "$(basename "$0")"
  printf "  %-24s : %s\n" "Action" "${ACTION}"
  printf "  %-24s : %s\n" "Mode" "$( [[ "${LOCAL_MODE}" -eq 1 ]] && printf 'LOCAL test' || printf 'CONDOR / aggregate' )"
  printf "  %-24s : %s\n" "BASE" "${BASE}"
  printf "  %-24s : %s\n" "GRL" "${GRL}"
  printf "  %-24s : %s\n" "Current list dir" "${LIST_DIR}"
  printf "  %-24s : %s\n" "Dataset" "${DATASET}"
  printf "  %-24s : %s\n" "Tag" "${TAG}"
  printf "  %-24s : %s\n" "Prefix" "${PREFIX}"
  printf "  %-24s : %s\n" "Verbosity" "${VERBOSE}"
  printf "  %-24s : %s\n" "Show full run tables" "$(yn "${SHOW_ALL_RUN_TABLES}")"
  printf "  %-24s : %s\n" "Max Condor jobs" "${MAX_CONDOR_JOBS}"
  if [[ -n "${MANIFEST_PATH}" ]]; then
    printf "  %-24s : %s\n" "Manifest" "${MANIFEST_PATH}"
  fi
}

main() {
  parse_args "$@"

  need_cmd awk
  need_cmd sed
  need_cmd sort
  need_cmd uniq
  need_cmd head
  need_cmd tail
  need_cmd tr
  need_cmd xargs
  need_cmd split
  need_cmd root
  need_cmd CreateDstList.pl

  [[ -f "${GRL}" ]] || fatal "GRL not found: ${GRL}"
  [[ -s "${GRL}" ]] || fatal "GRL is empty: ${GRL}"
  [[ -d "${LIST_DIR}" ]] || fatal "current list directory not found: ${LIST_DIR}"

  recompute_prefix_stem

  if [[ "${ACTION}" == "firstPass" ]]; then
    prepare_persistent_workdir
    read_runs_from_file "${GRL}"
    (( ${#RUN_ORDER[@]} > 0 )) || fatal "no valid runs found in GRL: ${GRL}"

    print_startup_context
    note "GRL runs loaded: $(fmt_num "${#RUN_ORDER[@]}")"

    build_expected_lists
    collect_unique_files
    write_root_counter_macro
    write_counter_wrapper
    build_count_chunks

    subsection "First-pass plan"
    note "Persisted work dir : ${TMP_ROOT}"
    note "Chunk dir          : ${TMP_CHUNK_DIR}"
    note "Results dir        : ${TMP_RESULTS_DIR}"
    note "Unique files       : $(fmt_num "${UNIQUE_FILE_COUNT}")"
    note "Chunk size         : $(fmt_num "${CHUNK_SIZE}")"
    note "Chunk count        : $(fmt_num "${CHUNK_COUNT}")"

    if (( LOCAL_MODE == 1 )); then
      run_local_counter_test
    else
      submit_condor_counter_jobs
    fi
    exit 0
  fi

  if [[ -z "${MANIFEST_PATH}" ]]; then
    MANIFEST_PATH="$(latest_manifest_in_pwd)"
  fi
  [[ -n "${MANIFEST_PATH}" ]] || fatal "no first-pass manifest found in the current directory; run firstPass first or use --manifest"

  load_firstpass_manifest "${MANIFEST_PATH}"
  read_runs_from_file "${GRL}"
  (( ${#RUN_ORDER[@]} > 0 )) || fatal "no valid runs found in GRL: ${GRL}"

  print_startup_context
  note "GRL runs loaded: $(fmt_num "${#RUN_ORDER[@]}")"

  collect_count_results
  audit_runs
  print_core_summary
  print_ranked_tables
  print_run_status_table

  good "Second-pass audit complete"
}

main "$@"
