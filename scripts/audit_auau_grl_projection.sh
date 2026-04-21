#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'
shopt -s nullglob

###############################################################################
# audit_auau_grl_projection.sh
#
# Purpose
#   Top-down, one-shot Au+Au GRL availability / projection audit.
#
# What it does
#   1) Reads the canonical Au+Au GRL.
#   2) Builds the canonical DB-backed expected DST_CALOFITTING per-run lists
#      on the fly via CreateDstList.pl in a TEMP dir (no permanent artifacts).
#   3) Compares those expected lists against the CURRENT project lists in:
#        /sphenix/u/patsfan753/scratch/thesisAnalysis/dst_lists_auau
#   4) Counts exact available/expected DST events by opening the ROOT files and
#      reading ONLY the event-bearing TTree entry counts (no full analysis run).
#   5) Pulls per-run runtime / GL1 / trigger bookkeeping from psql.
#   6) Prints meticulous terminal summaries + run-by-run tables.
#
# No permanent CSVs / artifacts are written.
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

TRIG_BIT=14
TRIG_EXPECT_NAME="MBD N&S >= 2, vtx < 150 cm"

PSQL_HOST="sphnxdaqdbreplica"
PSQL_DB="daq"

VERBOSE=1
RUN_PROGRESS_EVERY=100
ROOT_PROGRESS_EVERY=250
SHOW_ALL_RUN_TABLES=1

# Optional approximate proxy block ONLY if supplied by user
MBD_XSEC_MB=""

# Event-count method:
#   auto = try SQL/filelist first, then ROOT fallback for unresolved files
#   sql  = try SQL/filelist only
#   root = use the original ROOT-per-file scan only
COUNT_METHOD="auto"

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
  if [[ -n "${TMP_ROOT:-}" && -d "${TMP_ROOT:-}" ]]; then
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
  $(basename "$0") [options]

Options:
  --grl <path>                 Override GRL path
  --lists-dir <path>           Override current list directory
  --dataset <token>            Default: ${DATASET}
  --tag <token>                Default: ${TAG}
  --prefix <token>             Default: ${PREFIX}
  --trigger-bit <int>          Default: ${TRIG_BIT}
  --trigger-name <string>      Default: "${TRIG_EXPECT_NAME}"
  --mbd-xsec-mb <value>        Optional approximate proxy block
  --run-progress-every <N>     Default: ${RUN_PROGRESS_EVERY}
  --root-progress-every <N>    Default: ${ROOT_PROGRESS_EVERY}
  --summary-only               Suppress full run tables
  -q, --quiet                  Less chatter
  -v, --verbose                More chatter
  -h, --help                   Show this help

Examples:
  $(basename "$0")
  $(basename "$0") --summary-only
  $(basename "$0") --mbd-xsec-mb 25.0
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
# psql helpers
########################################
PSQL=()

sql() {
  "${PSQL[@]}" -c "$1" 2>/dev/null || true
}

check_psql_connectivity() {
  local probe
  probe="$(sql "SELECT 1;" | head -n 1 | tr -d '[:space:]')"
  [[ "$probe" == "1" ]] || {
    printf "%s[FATAL]%s unable to query psql at host=%s db=%s\n" \
      "${RED}${BOLD}" "${RESET}" "$PSQL_HOST" "$PSQL_DB" >&2
    exit 1
  }
}

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

read_run_runtime_gl1() {
  local run="$1"
  local rt gl1
  rt="$(sql "SELECT FLOOR(EXTRACT(EPOCH FROM (ertimestamp-brtimestamp)))::INT FROM run WHERE runnumber=${run};" \
        | head -n 1 | xargs)"
  gl1="$(sql "SELECT COALESCE(SUM(lastevent-firstevent+1),0)::BIGINT
              FROM filelist
              WHERE runnumber=${run}
                AND filename LIKE '%GL1_physics_gl1daq%.evt';" \
        | head -n 1 | xargs)"
  printf '%s\t%s\n' "$(num_or_zero "$rt")" "$(num_or_zero "$gl1")"
}

read_trigger_stats() {
  local run="$1"
  local name row menu_present active raw live scaled

  name="$(sql "SELECT COALESCE(triggername,'')
               FROM gl1_triggernames
               WHERE index=${TRIG_BIT}
                 AND ${run} BETWEEN runnumber AND runnumber_last
               ORDER BY runnumber DESC
               LIMIT 1;" | head -n 1 | sed 's/[[:space:]]*$//')"

  if [[ -n "$name" ]]; then
    menu_present=1
  else
    menu_present=0
  fi

  row="$(sql "SELECT
                COALESCE(MAX(CASE WHEN scaled != -1 THEN 1 ELSE 0 END),0),
                COALESCE(SUM(CASE WHEN raw    > 0 THEN raw    ELSE 0 END),0),
                COALESCE(SUM(CASE WHEN live   > 0 THEN live   ELSE 0 END),0),
                COALESCE(SUM(CASE WHEN scaled > 0 THEN scaled ELSE 0 END),0)
              FROM gl1_scalers
              WHERE runnumber=${run}
                AND index=${TRIG_BIT};" | head -n 1)"

  IFS=$'\t' read -r active raw live scaled <<< "${row:-0    0    0    0}"
  active="$(num_or_zero "${active:-0}")"
  raw="$(num_or_zero "${raw:-0}")"
  live="$(num_or_zero "${live:-0}")"
  scaled="$(num_or_zero "${scaled:-0}")"

  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${menu_present}" "${name}" "${active}" "${raw}" "${live}" "${scaled}"
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

run_sql_entry_scan() {
  local sql_input="${TMP_ROOT}/unique_files_sql.txt"
  local sql_results="${TMP_ROOT}/sql_results.tsv"
  local sql_file="${TMP_ROOT}/filelist_entry_lookup.sql"

  awk '
    {
      s = $0;
      sub(/^.*\//, "", s);
      if (s != "") print s;
    }' "${TMP_UNIQUE_FILES}" | sort -u > "${sql_input}"

  if [[ ! -s "${sql_input}" ]]; then
    warn "No filenames available for SQL filelist lookup."
    return 1
  fi

  cat > "${sql_file}" <<EOF
CREATE TEMP TABLE audit_files(filename text) ON COMMIT DROP;
\copy audit_files(filename) FROM '${sql_input}'
COPY (
  SELECT a.filename,
         COALESCE(SUM(COALESCE(f.lastevent - f.firstevent + 1, 0)),0)::BIGINT AS entries,
         COUNT(f.filename)::BIGINT AS matches
  FROM audit_files a
  LEFT JOIN filelist f
    ON f.filename = a.filename
  GROUP BY a.filename
  ORDER BY a.filename
) TO STDOUT WITH (FORMAT csv, DELIMITER E'\t', HEADER false);
EOF

  "${PSQL[@]}" -f "${sql_file}" > "${sql_results}" 2>/dev/null || return 1

  local -A sql_entries=()
  local -A sql_matches=()

  local fname entries matches
  while IFS=$'\t' read -r fname entries matches; do
    [[ -n "${fname:-}" ]] || continue
    sql_entries["$fname"]="$(num_or_zero "${entries:-0}")"
    sql_matches["$fname"]="$(num_or_zero "${matches:-0}")"
  done < "${sql_results}"

  local matched=0
  local path base match_count
  while IFS= read -r path; do
    [[ -n "${path:-}" ]] || continue
    base="${path##*/}"
    match_count="$(num_or_zero "${sql_matches[$base]:-0}")"
    if (( match_count > 0 )); then
      FILE_STATUS["$path"]="OK"
      FILE_ENTRIES["$path"]="${sql_entries[$base]:-0}"
      FILE_TREE["$path"]="SQL:filelist"
      ((matched+=1))
    fi
  done < "${TMP_UNIQUE_FILES}"

  if (( matched > 0 )); then
    note "SQL filelist matches found: $(fmt_num "${matched}")"
    return 0
  fi

  warn "SQL filelist lookup returned zero matches for the referenced DST filenames."
  return 1
}

run_root_entry_scan() {
  section "Step 3/5 — Event-entry counting (fast SQL/filelist first, ROOT fallback only if needed)"

  if (( UNIQUE_FILE_COUNT == 0 )); then
    warn "No unique files found to scan; event totals will be zero."
    TMP_ROOT_RESULTS="${TMP_ROOT}/root_results.tsv"
    : > "${TMP_ROOT_RESULTS}"
    return 0
  fi

  local unresolved_list="${TMP_ROOT}/unique_files_unresolved.txt"
  : > "${unresolved_list}"

  if [[ "${COUNT_METHOD}" != "root" ]]; then
    subsection "Step 3a — Fast SQL/filelist lookup"
    if run_sql_entry_scan; then
      :
    else
      if [[ "${COUNT_METHOD}" == "sql" ]]; then
        warn "COUNT_METHOD=sql requested, but SQL lookup did not resolve the referenced DST files."
      else
        warn "SQL lookup unavailable or incomplete; falling back to ROOT for unresolved files."
      fi
    fi
  fi

  local path
  while IFS= read -r path; do
    [[ -n "${path:-}" ]] || continue
    if [[ -z "${FILE_STATUS[$path]:-}" ]]; then
      printf '%s\n' "$path" >> "${unresolved_list}"
    fi
  done < "${TMP_UNIQUE_FILES}"

  local unresolved_count
  unresolved_count="$(awk 'NF{n++} END{print n+0}' "${unresolved_list}")"

  TMP_ROOT_RESULTS="${TMP_ROOT}/root_results.tsv"
  : > "${TMP_ROOT_RESULTS}"

  if [[ "${COUNT_METHOD}" == "sql" ]] && (( unresolved_count > 0 )); then
    warn "COUNT_METHOD=sql left $(fmt_num "${unresolved_count}") files unresolved; those files will remain uncounted."
  elif (( unresolved_count > 0 )); then
    subsection "Step 3b — ROOT fallback for unresolved files"
    write_root_counter_macro

    say "ROOT macro        : ${TMP_ROOT_MACRO}"
    say "Input file list   : ${unresolved_list}"
    say "Output result file: ${TMP_ROOT_RESULTS}"
    note "This does a fast SQL/filelist pass first, then opens only unresolved DST files in ROOT. It does NOT run your recoil analysis."

    local cmd
    cmd="${TMP_ROOT_MACRO}(\"${unresolved_list}\",\"${TMP_ROOT_RESULTS}\",${ROOT_PROGRESS_EVERY},${VERBOSE})"

    root -l -b -q "${cmd}"

    local status entries tree rootpath
    while IFS=$'\t' read -r status entries tree rootpath; do
      [[ -n "${rootpath:-}" ]] || continue
      FILE_STATUS["$rootpath"]="$status"
      FILE_ENTRIES["$rootpath"]="$(num_or_zero "${entries:-0}")"
      FILE_TREE["$rootpath"]="$tree"
    done < "${TMP_ROOT_RESULTS}"
  else
    good "No ROOT fallback needed; SQL resolved all referenced files."
  fi

  ROOT_OK_FILES=0
  ROOT_OPENFAIL_FILES=0
  ROOT_ZOMBIE_FILES=0
  ROOT_NOTREE_FILES=0
  ROOT_TOTAL_COUNTED_ENTRIES=0
  TREE_COUNT=()

  local tree_key
  while IFS= read -r path; do
    [[ -n "${path:-}" ]] || continue
    case "${FILE_STATUS[$path]:-UNSCANNED}" in
      OK)
        ((ROOT_OK_FILES+=1))
        ROOT_TOTAL_COUNTED_ENTRIES=$(( ROOT_TOTAL_COUNTED_ENTRIES + ${FILE_ENTRIES[$path]:-0} ))
        tree_key="${FILE_TREE[$path]:-UNKNOWN}"
        TREE_COUNT["$tree_key"]=$(( ${TREE_COUNT["$tree_key"]:-0} + 1 ))
        ;;
      OPENFAIL) ((ROOT_OPENFAIL_FILES+=1)) ;;
      ZOMBIE)   ((ROOT_ZOMBIE_FILES+=1))   ;;
      NOTREE)   ((ROOT_NOTREE_FILES+=1))   ;;
      *)        ;;
    esac
  done < "${TMP_UNIQUE_FILES}"

  good "Entry-count scan parsed"
  note "Files with counts available     : $(fmt_num "${ROOT_OK_FILES}")"
  note "Open failures                  : $(fmt_num "${ROOT_OPENFAIL_FILES}")"
  note "Zombie files                   : $(fmt_num "${ROOT_ZOMBIE_FILES}")"
  note "Files with no TTree            : $(fmt_num "${ROOT_NOTREE_FILES}")"
  note "Total counted entries available: $(fmt_num "${ROOT_TOTAL_COUNTED_ENTRIES}")"

  if (( ${#TREE_COUNT[@]} > 0 )); then
    subsection "Top detected count sources / event-tree paths"
    printf "  %-40s | %12s\n" "SourceOrTreePath" "Files"
    printf "  %s-+-%s\n" "$(printf '%*s' 40 '' | tr ' ' '-')" "$(printf '%*s' 12 '' | tr ' ' '-')"
    while IFS=$'\t' read -r count tree; do
      printf "  %-40s | %12s\n" "$(truncate_text "$tree" 40)" "$(fmt_num "$count")"
    done < <(
      for tree in "${!TREE_COUNT[@]}"; do
        printf '%s\t%s\n' "${TREE_COUNT[$tree]}" "$tree"
      done | sort -t $'\t' -k1,1nr -k2,2 | head -n 10
    )
  fi

  if (( ROOT_OPENFAIL_FILES > 0 || ROOT_ZOMBIE_FILES > 0 || ROOT_NOTREE_FILES > 0 )); then
    warn "Some files could not contribute event counts; totals below reflect only files with counts available."
  fi
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

declare -A RUN_RUNTIME_S=()
declare -A RUN_GL1_EVT=()

declare -A RUN_TRIG_MENU_PRESENT=()
declare -A RUN_TRIG_NAME=()
declare -A RUN_TRIG_ACTIVE=()
declare -A RUN_TRIG_RAW=()
declare -A RUN_TRIG_LIVE=()
declare -A RUN_TRIG_SCALED=()

declare -A RUN_STATUS=()
declare -A RUN_NOTES=()

# Totals / counters
TOT_RUNTIME_ALL=0
TOT_GL1_ALL=0
TOT_TRIG_RAW_ALL=0
TOT_TRIG_LIVE_ALL=0
TOT_TRIG_SCALED_ALL=0
TOT_TRIG_ACTIVE_RUNS_ALL=0

TOT_RUNTIME_CURRENT=0
TOT_GL1_CURRENT=0
TOT_TRIG_RAW_CURRENT=0
TOT_TRIG_LIVE_CURRENT=0
TOT_TRIG_SCALED_CURRENT=0
TOT_TRIG_ACTIVE_RUNS_CURRENT=0

TOT_RUNTIME_COMPLETE=0
TOT_GL1_COMPLETE=0
TOT_TRIG_RAW_COMPLETE=0
TOT_TRIG_LIVE_COMPLETE=0
TOT_TRIG_SCALED_COMPLETE=0
TOT_TRIG_ACTIVE_RUNS_COMPLETE=0

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
  section "Step 4/5 — Run-by-run audit (DB-backed expectation vs current lists vs exact ROOT entries vs psql)"

  local idx=0
  local r8 run expected_list current_list
  for r8 in "${RUN_ORDER[@]}"; do
    ((idx+=1))
    run=$((10#$r8))

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

    # Build sets
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

    # psql run bookkeeping
    local rt gl1
    IFS=$'\t' read -r rt gl1 < <(read_run_runtime_gl1 "$run")
    rt="$(num_or_zero "$rt")"
    gl1="$(num_or_zero "$gl1")"

    RUN_RUNTIME_S["$r8"]="$rt"
    RUN_GL1_EVT["$r8"]="$gl1"

    (( TOT_RUNTIME_ALL += rt ))
    (( TOT_GL1_ALL     += gl1 ))

    # psql trigger stats
    local menu_present trig_name trig_active trig_raw trig_live trig_scaled
    IFS=$'\t' read -r menu_present trig_name trig_active trig_raw trig_live trig_scaled < <(read_trigger_stats "$run")

    RUN_TRIG_MENU_PRESENT["$r8"]="$(num_or_zero "$menu_present")"
    RUN_TRIG_NAME["$r8"]="$trig_name"
    RUN_TRIG_ACTIVE["$r8"]="$(num_or_zero "$trig_active")"
    RUN_TRIG_RAW["$r8"]="$(num_or_zero "$trig_raw")"
    RUN_TRIG_LIVE["$r8"]="$(num_or_zero "$trig_live")"
    RUN_TRIG_SCALED["$r8"]="$(num_or_zero "$trig_scaled")"

    (( TOT_TRIG_RAW_ALL    += RUN_TRIG_RAW["$r8"]    ))
    (( TOT_TRIG_LIVE_ALL   += RUN_TRIG_LIVE["$r8"]   ))
    (( TOT_TRIG_SCALED_ALL += RUN_TRIG_SCALED["$r8"] ))
    if (( RUN_TRIG_ACTIVE["$r8"] != 0 )); then (( TOT_TRIG_ACTIVE_RUNS_ALL += 1 )); fi

    if (( current_seg > 0 )); then
      (( TOT_RUNTIME_CURRENT    += rt ))
      (( TOT_GL1_CURRENT        += gl1 ))
      (( TOT_TRIG_RAW_CURRENT   += RUN_TRIG_RAW["$r8"]    ))
      (( TOT_TRIG_LIVE_CURRENT  += RUN_TRIG_LIVE["$r8"]   ))
      (( TOT_TRIG_SCALED_CURRENT+= RUN_TRIG_SCALED["$r8"] ))
      if (( RUN_TRIG_ACTIVE["$r8"] != 0 )); then (( TOT_TRIG_ACTIVE_RUNS_CURRENT += 1 )); fi
    fi

    if (( expected_seg > 0 && missing_seg == 0 && extra_seg == 0 )); then
      (( TOT_RUNTIME_COMPLETE     += rt ))
      (( TOT_GL1_COMPLETE         += gl1 ))
      (( TOT_TRIG_RAW_COMPLETE    += RUN_TRIG_RAW["$r8"]    ))
      (( TOT_TRIG_LIVE_COMPLETE   += RUN_TRIG_LIVE["$r8"]   ))
      (( TOT_TRIG_SCALED_COMPLETE += RUN_TRIG_SCALED["$r8"] ))
      if (( RUN_TRIG_ACTIVE["$r8"] != 0 )); then (( TOT_TRIG_ACTIVE_RUNS_COMPLETE += 1 )); fi
    fi

    # Status / notes
    local notes=""
    if [[ -z "$expected_list" ]]; then append_note notes "NO_DB_EXPECTED_LIST"; fi
    if [[ -n "$expected_list" && ! -s "$expected_list" ]]; then append_note notes "EMPTY_DB_EXPECTED_LIST"; fi
    if [[ -z "$current_list" ]]; then append_note notes "NO_CURRENT_LIST"; fi
    if [[ -n "$current_list" && ! -s "$current_list" ]]; then append_note notes "EMPTY_CURRENT_LIST"; fi

    if [[ -n "$trig_name" ]]; then
      if [[ "$(normalize_name "$trig_name")" != "$(normalize_name "$TRIG_EXPECT_NAME")" ]]; then
        append_note notes "BIT${TRIG_BIT}_NAME_MISMATCH"
      fi
    else
      append_note notes "BIT${TRIG_BIT}_NO_MENU_ENTRY"
    fi

    if (( RUN_TRIG_ACTIVE["$r8"] == 0 )); then
      append_note notes "BIT${TRIG_BIT}_OFF"
    fi

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

  note "Expected segments are the current DB-backed canonical expectation (CreateDstList.pl in temp dir)."
  note "Event totals come from a fast SQL/filelist lookup first, with ROOT fallback only for unresolved files when COUNT_METHOD=auto."
  note "No permanent CSVs or output artifacts were created."

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

  subsection "C) Exact DST event availability (readable ROOT files only)"
  printf "  %-42s : %12s\n" "Expected exact DST events" "$(fmt_num "${TOT_EXPECTED_EVT}")"
  printf "  %-42s : %12s\n" "Present expected exact DST events" "$(fmt_num "${TOT_PRESENT_EVT}")"
  printf "  %-42s : %12s  (%s)\n" "Present/expected event completeness" "$(fmt_num "${TOT_PRESENT_EVT}")" "$(fmt_pct "${TOT_PRESENT_EVT}" "${TOT_EXPECTED_EVT}")"
  printf "  %-42s : %12s\n" "Missing expected exact DST events" "$(fmt_num "${TOT_MISSING_EVT}")"
  printf "  %-42s : %12s\n" "Current listed exact DST events" "$(fmt_num "${TOT_CURRENT_EVT}")"
  printf "  %-42s : %12s\n" "Extra current-listed exact DST events" "$(fmt_num "${TOT_EXTRA_EVT}")"

  subsection "D) Entry-count health"
  printf "  %-42s : %12s\n" "Unique referenced ROOT files scanned" "$(fmt_num "${UNIQUE_FILE_COUNT}")"
  printf "  %-42s : %12s\n" "Readable files with event tree" "$(fmt_num "${ROOT_OK_FILES}")"
  printf "  %-42s : %12s\n" "Open failures" "$(fmt_num "${ROOT_OPENFAIL_FILES}")"
  printf "  %-42s : %12s\n" "Zombie files" "$(fmt_num "${ROOT_ZOMBIE_FILES}")"
  printf "  %-42s : %12s\n" "Files with no tree" "$(fmt_num "${ROOT_NOTREE_FILES}")"

  subsection "E) psql exposure bookkeeping — ALL GRL runs"
  printf "  %-42s : %12s h\n" "Runtime over all GRL runs" "$(fmt_h "${TOT_RUNTIME_ALL}")"
  printf "  %-42s : %12s\n"   "GL1 .evt over all GRL runs" "$(fmt_num "${TOT_GL1_ALL}")"
  printf "  %-42s : %12s / %s\n" "Bit ${TRIG_BIT} active runs / GRL runs" "$(fmt_num "${TOT_TRIG_ACTIVE_RUNS_ALL}")" "$(fmt_num "${#RUN_ORDER[@]}")"
  printf "  %-42s : %12s\n"   "Bit ${TRIG_BIT} RAW total" "$(fmt_num "${TOT_TRIG_RAW_ALL}")"
  printf "  %-42s : %12s\n"   "Bit ${TRIG_BIT} LIVE total" "$(fmt_num "${TOT_TRIG_LIVE_ALL}")"
  printf "  %-42s : %12s\n"   "Bit ${TRIG_BIT} SCALED total" "$(fmt_num "${TOT_TRIG_SCALED_ALL}")"

  subsection "F) psql exposure bookkeeping — CURRENT analyzable subset (non-empty current lists)"
  printf "  %-42s : %12s h\n" "Runtime over current-list subset" "$(fmt_h "${TOT_RUNTIME_CURRENT}")"
  printf "  %-42s : %12s\n"   "GL1 .evt over current-list subset" "$(fmt_num "${TOT_GL1_CURRENT}")"
  printf "  %-42s : %12s / %s\n" "Bit ${TRIG_BIT} active runs / current-list runs" "$(fmt_num "${TOT_TRIG_ACTIVE_RUNS_CURRENT}")" "$(fmt_num "${RUNS_WITH_CURRENT}")"
  printf "  %-42s : %12s\n"   "Bit ${TRIG_BIT} RAW total" "$(fmt_num "${TOT_TRIG_RAW_CURRENT}")"
  printf "  %-42s : %12s\n"   "Bit ${TRIG_BIT} LIVE total" "$(fmt_num "${TOT_TRIG_LIVE_CURRENT}")"
  printf "  %-42s : %12s\n"   "Bit ${TRIG_BIT} SCALED total" "$(fmt_num "${TOT_TRIG_SCALED_CURRENT}")"

  subsection "G) psql exposure bookkeeping — FULLY COMPLETE runs only"
  printf "  %-42s : %12s h\n" "Runtime over fully complete runs" "$(fmt_h "${TOT_RUNTIME_COMPLETE}")"
  printf "  %-42s : %12s\n"   "GL1 .evt over fully complete runs" "$(fmt_num "${TOT_GL1_COMPLETE}")"
  printf "  %-42s : %12s / %s\n" "Bit ${TRIG_BIT} active runs / complete runs" "$(fmt_num "${TOT_TRIG_ACTIVE_RUNS_COMPLETE}")" "$(fmt_num "${RUNS_COMPLETE}")"
  printf "  %-42s : %12s\n"   "Bit ${TRIG_BIT} RAW total" "$(fmt_num "${TOT_TRIG_RAW_COMPLETE}")"
  printf "  %-42s : %12s\n"   "Bit ${TRIG_BIT} LIVE total" "$(fmt_num "${TOT_TRIG_LIVE_COMPLETE}")"
  printf "  %-42s : %12s\n"   "Bit ${TRIG_BIT} SCALED total" "$(fmt_num "${TOT_TRIG_SCALED_COMPLETE}")"

  if [[ -n "${MBD_XSEC_MB}" ]]; then
    subsection "H) Optional approximate luminosity proxy block"
    note "These are simple trigger-count / cross-section proxies, printed only because --mbd-xsec-mb was supplied."
    note "Interpret with caution; this script's main mission is availability/completeness, not authoritative luminosity calibration."

    local all_raw_l all_live_l all_scaled_l
    local cur_raw_l cur_live_l cur_scaled_l
    local comp_raw_l comp_live_l comp_scaled_l

    all_raw_l="$(awk -v n="${TOT_TRIG_RAW_ALL}"    -v x="${MBD_XSEC_MB}" 'BEGIN{ if (x > 0) printf "%.6f", n / x; else print "0"; }')"
    all_live_l="$(awk -v n="${TOT_TRIG_LIVE_ALL}"   -v x="${MBD_XSEC_MB}" 'BEGIN{ if (x > 0) printf "%.6f", n / x; else print "0"; }')"
    all_scaled_l="$(awk -v n="${TOT_TRIG_SCALED_ALL}" -v x="${MBD_XSEC_MB}" 'BEGIN{ if (x > 0) printf "%.6f", n / x; else print "0"; }')"

    cur_raw_l="$(awk -v n="${TOT_TRIG_RAW_CURRENT}"    -v x="${MBD_XSEC_MB}" 'BEGIN{ if (x > 0) printf "%.6f", n / x; else print "0"; }')"
    cur_live_l="$(awk -v n="${TOT_TRIG_LIVE_CURRENT}"   -v x="${MBD_XSEC_MB}" 'BEGIN{ if (x > 0) printf "%.6f", n / x; else print "0"; }')"
    cur_scaled_l="$(awk -v n="${TOT_TRIG_SCALED_CURRENT}" -v x="${MBD_XSEC_MB}" 'BEGIN{ if (x > 0) printf "%.6f", n / x; else print "0"; }')"

    comp_raw_l="$(awk -v n="${TOT_TRIG_RAW_COMPLETE}"    -v x="${MBD_XSEC_MB}" 'BEGIN{ if (x > 0) printf "%.6f", n / x; else print "0"; }')"
    comp_live_l="$(awk -v n="${TOT_TRIG_LIVE_COMPLETE}"   -v x="${MBD_XSEC_MB}" 'BEGIN{ if (x > 0) printf "%.6f", n / x; else print "0"; }')"
    comp_scaled_l="$(awk -v n="${TOT_TRIG_SCALED_COMPLETE}" -v x="${MBD_XSEC_MB}" 'BEGIN{ if (x > 0) printf "%.6f", n / x; else print "0"; }')"

    printf "  %-42s : %12s mb^-1\n" "Approx proxy L from ALL-RUN RAW / sigma" "${all_raw_l}"
    printf "  %-42s : %12s mb^-1\n" "Approx proxy L from ALL-RUN LIVE / sigma" "${all_live_l}"
    printf "  %-42s : %12s mb^-1\n" "Approx proxy L from ALL-RUN SCALED / sigma" "${all_scaled_l}"

    printf "  %-42s : %12s mb^-1\n" "Approx proxy L from CURRENT-LIST RAW / sigma" "${cur_raw_l}"
    printf "  %-42s : %12s mb^-1\n" "Approx proxy L from CURRENT-LIST LIVE / sigma" "${cur_live_l}"
    printf "  %-42s : %12s mb^-1\n" "Approx proxy L from CURRENT-LIST SCALED / sigma" "${cur_scaled_l}"

    printf "  %-42s : %12s mb^-1\n" "Approx proxy L from COMPLETE-RUN RAW / sigma" "${comp_raw_l}"
    printf "  %-42s : %12s mb^-1\n" "Approx proxy L from COMPLETE-RUN LIVE / sigma" "${comp_live_l}"
    printf "  %-42s : %12s mb^-1\n" "Approx proxy L from COMPLETE-RUN SCALED / sigma" "${comp_scaled_l}"
  fi
}

print_ranked_tables() {
  subsection "I) Ranked high-impact runs"

  printf "  %-8s | %-8s | %-8s | %-8s | %-12s | %-12s | %-8s | %-10s\n" \
    "Run" "ExpSeg" "PresSeg" "MissSeg" "PresEvt(M)" "MissEvt(M)" "TrigOn" "Status"
  printf "  %s-+-%s-+-%s-+-%s-+-%s-+-%s-+-%s-+-%s\n" \
    "$(printf '%*s' 8 '' | tr ' ' '-')" \
    "$(printf '%*s' 8 '' | tr ' ' '-')" \
    "$(printf '%*s' 8 '' | tr ' ' '-')" \
    "$(printf '%*s' 8 '' | tr ' ' '-')" \
    "$(printf '%*s' 12 '' | tr ' ' '-')" \
    "$(printf '%*s' 12 '' | tr ' ' '-')" \
    "$(printf '%*s' 8 '' | tr ' ' '-')" \
    "$(printf '%*s' 10 '' | tr ' ' '-')"

  while IFS=$'\t' read -r miss_evt miss_seg r8; do
    [[ -n "${r8:-}" ]] || continue
    printf "  %-8s | %8s | %8s | %8s | %12s | %12s | %8s | %-10s\n" \
      "$r8" \
      "$(fmt_num "${RUN_EXPECTED_SEG[$r8]:-0}")" \
      "$(fmt_num "${RUN_PRESENT_SEG[$r8]:-0}")" \
      "$(fmt_num "${RUN_MISSING_SEG[$r8]:-0}")" \
      "$(fmt_m "${RUN_PRESENT_EVT[$r8]:-0}")" \
      "$(fmt_m "${RUN_MISSING_EVT[$r8]:-0}")" \
      "$(yn "${RUN_TRIG_ACTIVE[$r8]:-0}")" \
      "${RUN_STATUS[$r8]:-UNK}"
  done < <(
    for r8 in "${RUN_ORDER[@]}"; do
      printf '%015d\t%015d\t%s\n' "${RUN_MISSING_EVT[$r8]:-0}" "${RUN_MISSING_SEG[$r8]:-0}" "$r8"
    done | sort -t $'\t' -k1,1nr -k2,2nr -k3,3 | head -n 15
  )

  subsection "J) Top runs by currently available exact DST events"
  printf "  %-8s | %-8s | %-12s | %-12s | %-8s | %-8s | %-10s\n" \
    "Run" "CurSeg" "CurEvt(M)" "GL1evt(M)" "TrigOn" "Live(M)" "Status"
  printf "  %s-+-%s-+-%s-+-%s-+-%s-+-%s-+-%s\n" \
    "$(printf '%*s' 8 '' | tr ' ' '-')" \
    "$(printf '%*s' 8 '' | tr ' ' '-')" \
    "$(printf '%*s' 12 '' | tr ' ' '-')" \
    "$(printf '%*s' 12 '' | tr ' ' '-')" \
    "$(printf '%*s' 8 '' | tr ' ' '-')" \
    "$(printf '%*s' 8 '' | tr ' ' '-')" \
    "$(printf '%*s' 10 '' | tr ' ' '-')"

  while IFS=$'\t' read -r cur_evt r8; do
    [[ -n "${r8:-}" ]] || continue
    printf "  %-8s | %8s | %12s | %12s | %8s | %8s | %-10s\n" \
      "$r8" \
      "$(fmt_num "${RUN_CURRENT_SEG[$r8]:-0}")" \
      "$(fmt_m "${RUN_CURRENT_EVT[$r8]:-0}")" \
      "$(fmt_m "${RUN_GL1_EVT[$r8]:-0}")" \
      "$(yn "${RUN_TRIG_ACTIVE[$r8]:-0}")" \
      "$(fmt_m "${RUN_TRIG_LIVE[$r8]:-0}")" \
      "${RUN_STATUS[$r8]:-UNK}"
  done < <(
    for r8 in "${RUN_ORDER[@]}"; do
      printf '%015d\t%s\n' "${RUN_CURRENT_EVT[$r8]:-0}" "$r8"
    done | sort -t $'\t' -k1,1nr -k2,2 | head -n 15
  )

  subsection "K) Runs with current lists but trigger bit ${TRIG_BIT} OFF"
  local printed=0
  printf "  %-8s | %-8s | %-12s | %-12s | %-24s\n" \
    "Run" "CurSeg" "CurEvt(M)" "GL1evt(M)" "Notes"
  printf "  %s-+-%s-+-%s-+-%s-+-%s\n" \
    "$(printf '%*s' 8 '' | tr ' ' '-')" \
    "$(printf '%*s' 8 '' | tr ' ' '-')" \
    "$(printf '%*s' 12 '' | tr ' ' '-')" \
    "$(printf '%*s' 12 '' | tr ' ' '-')" \
    "$(printf '%*s' 24 '' | tr ' ' '-')"

  while IFS=$'\t' read -r r8; do
    [[ -n "${r8:-}" ]] || continue
    printf "  %-8s | %8s | %12s | %12s | %-24s\n" \
      "$r8" \
      "$(fmt_num "${RUN_CURRENT_SEG[$r8]:-0}")" \
      "$(fmt_m "${RUN_CURRENT_EVT[$r8]:-0}")" \
      "$(fmt_m "${RUN_GL1_EVT[$r8]:-0}")" \
      "$(truncate_text "${RUN_NOTES[$r8]:-}" 24)"
    ((printed+=1))
  done < <(
    for r8 in "${RUN_ORDER[@]}"; do
      if (( ${RUN_CURRENT_SEG[$r8]:-0} > 0 )) && (( ${RUN_TRIG_ACTIVE[$r8]:-0} == 0 )); then
        printf '%s\n' "$r8"
      fi
    done | sort
  )

  if (( printed == 0 )); then
    printf "  %-8s | %-8s | %-12s | %-12s | %-24s\n" "-" "-" "-" "-" "none"
  fi
}

print_run_status_table() {
  if (( SHOW_ALL_RUN_TABLES == 0 )); then
    return 0
  fi

  subsection "L) Full per-run audit table"

  printf "  %-8s | %-6s | %-6s | %-6s | %-6s | %-6s | %-10s | %-10s | %-5s | %-9s | %-8s | %-28s\n" \
    "Run" "Exp" "Cur" "Pres" "Miss" "Extra" "PresEvt(M)" "MissEvt(M)" "Bit14" "Live(M)" "GL1(M)" "Notes"
  printf "  %s-+-%s-+-%s-+-%s-+-%s-+-%s-+-%s-+-%s-+-%s-+-%s-+-%s-+-%s\n" \
    "$(printf '%*s' 8 '' | tr ' ' '-')" \
    "$(printf '%*s' 6 '' | tr ' ' '-')" \
    "$(printf '%*s' 6 '' | tr ' ' '-')" \
    "$(printf '%*s' 6 '' | tr ' ' '-')" \
    "$(printf '%*s' 6 '' | tr ' ' '-')" \
    "$(printf '%*s' 6 '' | tr ' ' '-')" \
    "$(printf '%*s' 10 '' | tr ' ' '-')" \
    "$(printf '%*s' 10 '' | tr ' ' '-')" \
    "$(printf '%*s' 5 '' | tr ' ' '-')" \
    "$(printf '%*s' 9 '' | tr ' ' '-')" \
    "$(printf '%*s' 8 '' | tr ' ' '-')" \
    "$(printf '%*s' 28 '' | tr ' ' '-')"

  local r8
  for r8 in "${RUN_ORDER[@]}"; do
    printf "  %-8s | %6s | %6s | %6s | %6s | %6s | %10s | %10s | %5s | %9s | %8s | %-28s\n" \
      "$r8" \
      "$(fmt_num "${RUN_EXPECTED_SEG[$r8]:-0}")" \
      "$(fmt_num "${RUN_CURRENT_SEG[$r8]:-0}")" \
      "$(fmt_num "${RUN_PRESENT_SEG[$r8]:-0}")" \
      "$(fmt_num "${RUN_MISSING_SEG[$r8]:-0}")" \
      "$(fmt_num "${RUN_EXTRA_SEG[$r8]:-0}")" \
      "$(fmt_m "${RUN_PRESENT_EVT[$r8]:-0}")" \
      "$(fmt_m "${RUN_MISSING_EVT[$r8]:-0}")" \
      "$(yn "${RUN_TRIG_ACTIVE[$r8]:-0}")" \
      "$(fmt_m "${RUN_TRIG_LIVE[$r8]:-0}")" \
      "$(fmt_m "${RUN_GL1_EVT[$r8]:-0}")" \
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
      --trigger-bit)
        shift
        [[ $# -gt 0 ]] || fatal "missing value after --trigger-bit"
        TRIG_BIT="$1"
        shift
        ;;
      --trigger-name)
        shift
        [[ $# -gt 0 ]] || fatal "missing value after --trigger-name"
        TRIG_EXPECT_NAME="$1"
        shift
        ;;
      --mbd-xsec-mb)
        shift
        [[ $# -gt 0 ]] || fatal "missing value after --mbd-xsec-mb"
        MBD_XSEC_MB="$1"
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

  [[ "$TRIG_BIT" =~ ^[0-9]+$ ]] || fatal "--trigger-bit must be an integer"
  [[ "$RUN_PROGRESS_EVERY" =~ ^[0-9]+$ ]] || fatal "--run-progress-every must be an integer"
  [[ "$ROOT_PROGRESS_EVERY" =~ ^[0-9]+$ ]] || fatal "--root-progress-every must be an integer"
  if [[ -n "${MBD_XSEC_MB}" ]]; then
    [[ "${MBD_XSEC_MB}" =~ ^[0-9]+([.][0-9]+)?$ ]] || fatal "--mbd-xsec-mb must be numeric"
  fi
}

print_startup_context() {
  section "Au+Au GRL projection audit"

  printf "  %-24s : %s\n" "Script" "$(basename "$0")"
  printf "  %-24s : %s\n" "BASE" "${BASE}"
  printf "  %-24s : %s\n" "GRL" "${GRL}"
  printf "  %-24s : %s\n" "Current list dir" "${LIST_DIR}"
  printf "  %-24s : %s\n" "Dataset" "${DATASET}"
  printf "  %-24s : %s\n" "Tag" "${TAG}"
  printf "  %-24s : %s\n" "Prefix" "${PREFIX}"
  printf "  %-24s : %s\n" "Trigger bit" "${TRIG_BIT}"
  printf "  %-24s : %s\n" "Expected trigger name" "${TRIG_EXPECT_NAME}"
  printf "  %-24s : %s\n" "Verbosity" "${VERBOSE}"
  printf "  %-24s : %s\n" "Show full run tables" "$(yn "${SHOW_ALL_RUN_TABLES}")"
  if [[ -n "${MBD_XSEC_MB}" ]]; then
    printf "  %-24s : %s mb\n" "Optional MBD xsec" "${MBD_XSEC_MB}"
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
  need_cmd mktemp
  need_cmd psql
  need_cmd root
  need_cmd CreateDstList.pl

  [[ -f "${GRL}" ]] || fatal "GRL not found: ${GRL}"
  [[ -s "${GRL}" ]] || fatal "GRL is empty: ${GRL}"
  [[ -d "${LIST_DIR}" ]] || fatal "current list directory not found: ${LIST_DIR}"

  recompute_prefix_stem

  TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/audit_auau_grl_projection.XXXXXX")"
  PSQL=(psql -h "${PSQL_HOST}" -d "${PSQL_DB}" -At -F $'\t' -q)

  check_psql_connectivity
  read_runs_from_file "${GRL}"
  (( ${#RUN_ORDER[@]} > 0 )) || fatal "no valid runs found in GRL: ${GRL}"

  print_startup_context
  note "GRL runs loaded: $(fmt_num "${#RUN_ORDER[@]}")"

  build_expected_lists
  collect_unique_files
  run_root_entry_scan
  audit_runs
  print_core_summary
  print_ranked_tables
  print_run_status_table

  good "Audit complete"
}

main "$@"
