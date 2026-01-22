#!/usr/bin/env bash
# ==============================================================================
# makeThesisSimLists.sh
# ==============================================================================
# PURPOSE
#   Generate analysis-ready *.list files for Run-23 Pythia8 jets embedded into
#   Au+Au (HIJING) at two energies:
#     • type 11 → 30 GeV jets
#     • type 12 → 10 GeV jets
#
#   For each type, we request the following filetypes from the catalog:
#     DST_CALO_CLUSTER  DST_CALO_G4HIT  DST_GLOBAL  DST_MBD_EPD
#     DST_TRUTH_G4HIT   DST_TRUTH_JET   G4Hits
#
#   Output (created if missing):
#     /sphenix/u/patsfan753/scratch/thesisAnalysis/simFiles/
#       ├─ run23_type11_auau/   (30 GeV jets)
#       └─ run23_type12_auau/   (10 GeV jets)
#
#   Diagnostics per pack:
#     • Line counts per list
#     • Pairing check: segment-key overlap across lists (anchor = DST_CALO_CLUSTER)
#       - Reports coverage and a small sample of missing keys, if any
#
# USAGE
#   ./makeThesisSimLists.sh
#   Optional flags:
#     --outroot DIR       (default: /sphenix/u/patsfan753/scratch/thesisAnalysis/simFiles)
#     --head N            (preview first N lines of each list; default: 8; 0 = no preview)
#     --quiet             (suppress INFO logs)
#     -h | --help
#
# REQUIREMENTS
#   • sPHENIX environment sourced; CreateFileList.pl must be in PATH
#   • The catalog must know about: -type {11,12} -run 23 -embed auau
#
# NOTES
#   • We normalize list filenames to UPPERCASE (e.g., DST_TRUTH_JET.list, G4Hits.list)
#   • Segment-pairing keys are inferred from the trailing "-<number>-<number>.root"
#     (best effort; robust against most standard sPHENIX filenames).
# ==============================================================================

set -Eeuo pipefail

# --------------------------- Defaults -----------------------------------------
# Base output directory you asked for:
OUTROOT="/sphenix/u/patsfan753/scratch/thesisAnalysis/simListFiles"
HEAD_TAIL_LINES="8"
VERBOSE="true"

# Fixed parameters for THIS script: Run-28 Pythia8 photon+jet (filesystem truth)
# We will NOT use the catalog. We will scan lustre directories directly.
RUNNUM="28"
#PHOTONJET_SAMPLES=( "photonjet10" "photonjet20" )
PHOTONJET_SAMPLES=( "photonjet10" "photonjet20" )

# Source base (what you just proved exists)
MDC2_BASE="/sphenix/lustre01/sphnxpro/mdc2/js_pp200_signal"

# These are the exact directories you found:
#   g4hits/run0028/photonjetX
#   nopileup/calocluster/run0028/photonjetX
#   nopileup/global/run0028/photonjetX
#   nopileup/jets/run0028/photonjetX
#   nopileup/mbdepd/run0028/photonjetX
#   nopileup/trkrhit/run0028/photonjetX


# --------------------------- CLI ----------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --outroot) OUTROOT="${2:-}"; shift 2 ;;
    --head)    HEAD_TAIL_LINES="${2:-8}"; shift 2 ;;
    --quiet)   VERBOSE="false"; shift 1 ;;
    -h|--help)
      sed -n '1,200p' "$0" | sed 's/^# \{0,1\}//g'
      exit 0
      ;;
    *) echo "[WARN $(date '+%H:%M:%S')] Unknown arg: $1"; shift 1 ;;
  esac
done

# --------------------------- Logging helpers ----------------------------------
ts() { date '+%Y-%m-%d %H:%M:%S'; }

# Optional deep debug: export DEBUG=1 or DEBUG=true
DEBUG="${DEBUG:-false}"
[[ "$DEBUG" == "1" ]] && DEBUG="true"

# ANSI color (stderr only). Auto-disables if stderr isn't a TTY or NO_COLOR is set.
USE_COLOR="true"
if [[ ! -t 2 || -n "${NO_COLOR:-}" ]]; then USE_COLOR="false"; fi

if [[ "$USE_COLOR" == "true" ]]; then
  C_RESET=$'\033[0m'
  C_BOLD=$'\033[1m'
  C_DIM=$'\033[2m'
  C_RED=$'\033[31m'
  C_GRN=$'\033[32m'
  C_YLW=$'\033[33m'
  C_BLU=$'\033[34m'
  C_MAG=$'\033[35m'
  C_CYN=$'\033[36m'
  C_WHT=$'\033[37m'
else
  C_RESET=""; C_BOLD=""; C_DIM=""
  C_RED=""; C_GRN=""; C_YLW=""
  C_BLU=""; C_MAG=""; C_CYN=""; C_WHT=""
fi

# Pretty separator + step counter (all output goes to STDERR)
RULE_W=79
hr() { printf "%b%s%b\n" "$C_DIM" "$(printf '%*s' "$RULE_W" '' | tr ' ' '-')" "$C_RESET" >&2; }
rule() { hr; }  # backward-compatible with your existing calls

_log() {
  local tag="$1"; shift
  local color="$1"; shift
  printf "%b[%s %s]%b %s\n" "${color}${C_BOLD}" "$tag" "$(ts)" "$C_RESET" "$*" >&2
}

say() {
  if [[ "$VERBOSE" == "true" ]]; then
    _log "INFO" "$C_CYN" "$*"
  fi
  return 0
}

ok() {
  if [[ "$VERBOSE" == "true" ]]; then
    _log "OK  " "$C_GRN" "$*"
  fi
  return 0
}

dbg() {
  if [[ "$DEBUG" == "true" ]]; then
    _log "DBG " "$C_MAG" "$*"
  fi
  return 0
}
warn() { _log "WARN" "$C_YLW" "$*"; }
die()  { _log "ERR " "$C_RED" "$*"; exit 1; }

STEP=0
step() {
  STEP=$((STEP + 1))
  printf "%b==>%b %b[%02d]%b %s\n" "$C_BLU" "$C_RESET" "$C_BOLD" "$STEP" "$C_RESET" "$*" >&2
}

# Safe line counter (never trips set -e)
count_lines() {
  local f="$1"
  if [[ -f "$f" ]]; then
    local n
    n="$(wc -l "$f" 2>/dev/null | awk '{print $1}' || true)"
    [[ -n "$n" ]] && echo "$n" || echo "0"
  else
    echo "MISSING"
  fi
}

# Global run summary accumulators (filled per PACK, printed at end)
START_EPOCH="$(date +%s)"
START_HUMAN="$(ts)"

declare -a PACK_ORDER=()
declare -A PACK_OUTDIR=()
declare -A PACK_REMOVED=()
declare -A PACK_RAW_COUNTS=()
declare -A PACK_MATCH_COUNTS=()
declare -A PACK_PAIR_COUNTS=()

record_pack() {
  local tag="$1"
  local outdir="$2"

  PACK_ORDER+=( "$tag" )
  PACK_OUTDIR["$tag"]="$outdir"

  PACK_RAW_COUNTS["$tag"]="DST_CALO_CLUSTER=$(count_lines "${outdir}/DST_CALO_CLUSTER.list")  G4Hits=$(count_lines "${outdir}/G4Hits.list")  DST_JETS=$(count_lines "${outdir}/DST_JETS.list")  DST_GLOBAL=$(count_lines "${outdir}/DST_GLOBAL.list")  DST_MBD_EPD=$(count_lines "${outdir}/DST_MBD_EPD.list")"
  PACK_MATCH_COUNTS["$tag"]="DST_CALO_CLUSTER=$(count_lines "${outdir}/DST_CALO_CLUSTER.matched.list")  G4Hits=$(count_lines "${outdir}/G4Hits.matched.list")  DST_JETS=$(count_lines "${outdir}/DST_JETS.matched.list")  DST_GLOBAL=$(count_lines "${outdir}/DST_GLOBAL.matched.list")  DST_MBD_EPD=$(count_lines "${outdir}/DST_MBD_EPD.matched.list")"
  PACK_PAIR_COUNTS["$tag"]="Calo__G4Hits=$(count_lines "${outdir}/DST_CALO_CLUSTER__G4Hits.pairs.list")  Calo__DST_JETS=$(count_lines "${outdir}/DST_CALO_CLUSTER__DST_JETS.pairs.list")  Calo__DST_GLOBAL=$(count_lines "${outdir}/DST_CALO_CLUSTER__DST_GLOBAL.pairs.list")  Calo__DST_MBD_EPD=$(count_lines "${outdir}/DST_CALO_CLUSTER__DST_MBD_EPD.pairs.list")  TripA=$(count_lines "${outdir}/DST_CALO_CLUSTER__G4Hits__DST_JETS.triplets.list")  TripB=$(count_lines "${outdir}/DST_CALO_CLUSTER__DST_GLOBAL__DST_MBD_EPD.triplets.list")"
}

final_summary() {
  local end_epoch end_human dt
  end_epoch="$(date +%s)"
  end_human="$(ts)"
  dt=$(( end_epoch - START_EPOCH ))

  rule
  echo "[FINAL SUMMARY] makeThesisSimLists.sh" >&2
  echo "Started : ${START_HUMAN}" >&2
  echo "Finished: ${end_human}" >&2
  echo "Runtime : ${dt}s" >&2
  echo "OUTROOT : ${OUTROOT}" >&2
  echo "Packs   : ${#PACK_ORDER[@]}" >&2
  echo "" >&2

  local i=0
  local tag
  for tag in "${PACK_ORDER[@]}"; do
    i=$((i + 1))
    local outdir="${PACK_OUTDIR[$tag]}"

    echo "[$i] ${tag}" >&2
    echo "    outdir  : ${outdir}" >&2
    echo "    cleanup : ${PACK_REMOVED[$tag]:-(none)}" >&2
    echo "    raw     : ${PACK_RAW_COUNTS[$tag]:-N/A}" >&2
    echo "    matched : ${PACK_MATCH_COUNTS[$tag]:-N/A}" >&2
    echo "    pairs   : ${PACK_PAIR_COUNTS[$tag]:-N/A}" >&2
    echo "    report  : ${outdir}/pair_report.txt" >&2
    echo "    summary : ${outdir}/summary.txt" >&2
    echo "    pipeline input (paired list):" >&2
    echo "      ${outdir}/DST_CALO_CLUSTER__G4Hits.pairs.list" >&2

    if [[ -d "$outdir" ]]; then
      echo "    files created in outdir:" >&2
      ls -1 "$outdir" 2>/dev/null | sed 's/^/      - /' >&2 || true
    fi
    echo "" >&2
  done

  rule
}

# Print a hard, immediate failure context on any error
on_err() {
  local ec=$?
  local line="${BASH_LINENO[0]:-unknown}"
  local src="${BASH_SOURCE[1]:-${BASH_SOURCE[0]:-unknown}}"
  local cmd="${BASH_COMMAND:-unknown}"

  echo >&2
  echo "[ERR  $(ts)] Script FAILED (exit=${ec})" >&2
  echo "[ERR  $(ts)] Location: ${src}:${line}" >&2
  echo "[ERR  $(ts)] Command : ${cmd}" >&2

  if [[ ${#FUNCNAME[@]} -gt 1 ]]; then
    echo "[ERR  $(ts)] Call stack:" >&2
    local i
    for (( i=1; i<${#FUNCNAME[@]}; ++i )); do
      echo "  at ${FUNCNAME[$i]}()  ${BASH_SOURCE[$i]}:${BASH_LINENO[$((i-1))]}" >&2
    done
  fi
  exit "$ec"
}
trap on_err ERR

# Optional xtrace with line numbers
if [[ "$DEBUG" == "true" ]]; then
  export PS4='+ [TRACE ${BASH_SOURCE##*/}:${LINENO}:${FUNCNAME[0]}] '
  set -x
fi

# --------------------------- Utilities -----------------------------------------
resolve_catalog() {
  for c in CreateFileList.pl CreateFileLists.pl CreateDstList.pl; do
    if command -v "$c" >/dev/null 2>&1; then echo "$c"; return 0; fi
  done
  die "No catalog tool found. Source env (e.g., 'source /opt/sphenix/core/bin/sphenix_setup.sh -n')."
}

# Map API filetype → expected lowercase list filename produced by the catalog
type_to_lower_list() {
  local t="$1"
  if [[ "$t" == "G4Hits" ]]; then
    echo "g4hits.list"
  else
    # Lowercase the API token; typical catalog behavior
    echo "$(tr 'A-Z' 'a-z' <<<"$t").list"
  fi
}

# Extract a robust "segment key" from a ROOT filepath
# Returns "<NUM>-<NUM>" if the last two dash-separated fields are numeric; otherwise a best-effort token
extract_keys() {
  # shellcheck disable=SC2013
  for p in $(cat "$1"); do
    local b; b="$(basename "$p" .root)"
    # Split on '-' and try to use the last two numeric fields as the key:
    # e.g., foo-bar-0000001234-0000  -> 0000001234-0000
    local IFS='-'; read -r -a parts <<< "$b"
    local n="${#parts[@]}"
    if (( n >= 2 )); then
      local last="${parts[$((n-1))]}"
      local prev="${parts[$((n-2))]}"
      if [[ "$last" =~ ^[0-9]+$ && "$prev" =~ ^[0-9]+$ ]]; then
        echo "${prev}-${last}"
        continue
      fi
    fi
    # Fallback: emit the basename if no standard pattern detected
    echo "$b"
  done | sort -u
}

pair_check() {
  # Args: <outdir> <anchor_list_name> <other_list_name>...
  local outdir="$1"; shift
  local anchor="$1"; shift
  local others=( "$@" )

  local report="${outdir}/pair_report.txt"
  : > "$report"

  say "Pairing: anchor=${anchor}  report=$(basename "$report")"
  dbg "Pair report path: $report"

  local anchor_list="${outdir}/${anchor}"
  [[ -s "$anchor_list" ]] || die "Anchor list missing/empty: ${anchor_list}"

  # Temp dir (function-local cleanup)
  local tmpdir; tmpdir="$(mktemp -d)"
  trap '[[ -n "${tmpdir:-}" ]] && rm -rf "${tmpdir}"' RETURN

  local anchor_base="${anchor%.list}"

  # ---------------- Build anchor keys + map ----------------
  awk '
    function key_of(path,    b, n, a) {
      b=path; sub(/^.*\//,"",b); sub(/\.root$/,"",b);
      n=split(b,a,"-");
      if(n>=2 && a[n] ~ /^[0-9]+$/ && a[n-1] ~ /^[0-9]+$/) return a[n-1] "-" a[n];
      return b;
    }
    { print key_of($0) }' "$anchor_list" | sort -u > "${tmpdir}/${anchor_base}.keys"

  awk '
    function key_of(path,    b, n, a) {
      b=path; sub(/^.*\//,"",b); sub(/\.root$/,"",b);
      n=split(b,a,"-");
      if(n>=2 && a[n] ~ /^[0-9]+$/ && a[n-1] ~ /^[0-9]+$/) return a[n-1] "-" a[n];
      return b;
    }
    { k=key_of($0); print k "\t" $0 }' "$anchor_list" > "${tmpdir}/${anchor_base}.map"

  local a_count; a_count=$(wc -l < "${tmpdir}/${anchor_base}.keys" | tr -d ' ')
  printf "ANCHOR: %s  keys=%d\n" "$anchor" "$a_count" >> "$report"

  # Terminal header table
  if [[ "$VERBOSE" == "true" ]]; then
    printf "%b%-9s%b %-22s  %s\n" "$C_DIM" "STATUS" "$C_RESET" "LIST" "common/anchor (pct)   miss  extra" >&2
    printf "%b%-9s%b %-22s  %s\n" "$C_DIM" "--------" "$C_RESET" "----------------------" "------------------------------" >&2
  fi

  # common_all starts as anchor.keys, then intersect with each other.keys
  cp -f "${tmpdir}/${anchor_base}.keys" "${tmpdir}/common_all.keys"

  # Store maps so we can write matched/pairs AFTER the final intersection is known
  declare -A MAP
  MAP["${anchor_base}"]="${tmpdir}/${anchor_base}.map"

  # ---------------- Process each required "other" list ----------------
  local oname
  for oname in "${others[@]}"; do
    local other_list="${outdir}/${oname}"
    [[ -s "$other_list" ]] || die "Required list missing/empty: ${other_list}"

    local obase="${oname%.list}"

    # other keys + map
    awk '
      function key_of(path,    b, n, a) {
        b=path; sub(/^.*\//,"",b); sub(/\.root$/,"",b);
        n=split(b,a,"-");
        if(n>=2 && a[n] ~ /^[0-9]+$/ && a[n-1] ~ /^[0-9]+$/) return a[n-1] "-" a[n];
        return b;
      }
      { print key_of($0) }' "$other_list" | sort -u > "${tmpdir}/${obase}.keys"

    awk '
      function key_of(path,    b, n, a) {
        b=path; sub(/^.*\//,"",b); sub(/\.root$/,"",b);
        n=split(b,a,"-");
        if(n>=2 && a[n] ~ /^[0-9]+$/ && a[n-1] ~ /^[0-9]+$/) return a[n-1] "-" a[n];
        return b;
      }
      { k=key_of($0); print k "\t" $0 }' "$other_list" > "${tmpdir}/${obase}.map"

    MAP["${obase}"]="${tmpdir}/${obase}.map"

    # reporting vs ANCHOR
    comm -12 "${tmpdir}/${anchor_base}.keys" "${tmpdir}/${obase}.keys" > "${tmpdir}/common_${obase}.keys"
    comm -23 "${tmpdir}/${anchor_base}.keys" "${tmpdir}/${obase}.keys" > "${tmpdir}/missing_in_${obase}.keys"
    comm -13 "${tmpdir}/${anchor_base}.keys" "${tmpdir}/${obase}.keys" > "${tmpdir}/extra_in_${obase}.keys"

    local c_count m_count e_count
    c_count=$(wc -l < "${tmpdir}/common_${obase}.keys" | tr -d ' ')
    m_count=$(wc -l < "${tmpdir}/missing_in_${obase}.keys" | tr -d ' ')
    e_count=$(wc -l < "${tmpdir}/extra_in_${obase}.keys" | tr -d ' ')

    local pct
    pct="$(awk -v c="$c_count" -v a="$a_count" 'BEGIN{if(a>0) printf("%.1f",100*c/a); else print "0.0"}')"

    printf "LIST  : %s  common=%d / %d (%s%%)  missing=%d  extra=%d\n" \
      "$oname" "$c_count" "$a_count" "$pct" "$m_count" "$e_count" >> "$report"

    # Terminal: compact status row
    if [[ "$VERBOSE" == "true" ]]; then
      local badge="[OK]"
      local bcol="$C_GRN"
      if (( m_count > 0 )); then
        badge="[MISS]"
        bcol="$C_RED"
      elif (( e_count > 0 )); then
        badge="[EXTRA]"
        bcol="$C_YLW"
      fi
      printf "%b%-9s%b %-22s  %3d/%3d (%5s%%)   %4d  %5d\n" \
        "$bcol" "$badge" "$C_RESET" "$oname" "$c_count" "$a_count" "$pct" "$m_count" "$e_count" >&2
    fi

    # Persist missing keys (if any)
    if (( m_count > 0 )); then
      local miss_keys_out="${outdir}/missing_in_${obase}.keys"
      cp -f "${tmpdir}/missing_in_${obase}.keys" "$miss_keys_out"
      printf "  missing-keys file -> %s\n" "$miss_keys_out" >> "$report"
      printf "%b  missing-keys%b file -> %s\n" "$C_RED" "$C_RESET" "$miss_keys_out" >&2

      if [[ "$DEBUG" == "true" ]]; then
        printf "%b  sample missing:%b %s\n" "$C_RED" "$C_RESET" "$(head -n 8 "$miss_keys_out" | paste -sd',' -)" >&2
      fi
    fi

    # Persist extra keys (if any)
    if (( e_count > 0 )); then
      local extra_keys_out="${outdir}/extra_in_${obase}.keys"
      cp -f "${tmpdir}/extra_in_${obase}.keys" "$extra_keys_out"
      printf "  extra-keys file -> %s\n" "$extra_keys_out" >> "$report"
      printf "%b  extra-keys%b  file -> %s\n" "$C_YLW" "$C_RESET" "$extra_keys_out" >&2

      if [[ "$DEBUG" == "true" ]]; then
        printf "%b  sample extra:%b %s\n" "$C_YLW" "$C_RESET" "$(head -n 8 "$extra_keys_out" | paste -sd',' -)" >&2
      fi
    fi

    # dropped-from-anchor paths (anchor keys missing in this other list)
    local dropped_anchor="${outdir}/dropped_from_${anchor_base}_vs_${obase}.list"
    awk 'NR==FNR{drop[$1]=1; next} ($1 in drop){print $2}' \
      "${tmpdir}/missing_in_${obase}.keys" "${tmpdir}/${anchor_base}.map" | sort -V > "$dropped_anchor"

    # update common_all = intersection(common_all, other.keys)
    comm -12 "${tmpdir}/common_all.keys" "${tmpdir}/${obase}.keys" > "${tmpdir}/common_all.next"
    mv -f "${tmpdir}/common_all.next" "${tmpdir}/common_all.keys"
  done

  local common_all_count
  common_all_count=$(wc -l < "${tmpdir}/common_all.keys" | tr -d ' ')
  printf "COMMON(all lists): keys=%d\n" "$common_all_count" >> "$report"

  (( common_all_count > 0 )) || die "No common segment keys across ALL lists in ${outdir}. See: ${report}"

  ok "Pairing OK: common keys across ALL required lists = ${common_all_count}"
  dbg "Report: ${report}"

  # ---------------- Write matched lists (aligned to common_all) ----------------
  local anchor_matched="${outdir}/${anchor_base}.matched.list"
  awk 'NR==FNR{keys[++n]=$1; next} {path[$1]=$2}
       END{for(i=1;i<=n;i++){k=keys[i]; if(k in path) print path[k]}}' \
      "${tmpdir}/common_all.keys" "${MAP[$anchor_base]}" > "$anchor_matched"

  for oname in "${others[@]}"; do
    local obase="${oname%.list}"
    local other_matched="${outdir}/${obase}.matched.list"
    awk 'NR==FNR{keys[++n]=$1; next} {path[$1]=$2}
         END{for(i=1;i<=n;i++){k=keys[i]; if(k in path) print path[k]}}' \
        "${tmpdir}/common_all.keys" "${MAP[$obase]}" > "$other_matched"
  done

  # ---------------- Write 2-column pair lists (aligned to common_all) ----------------
  for oname in "${others[@]}"; do
    local obase="${oname%.list}"
    local pairfile="${outdir}/${anchor_base}__${obase}.pairs.list"

    awk '
      FNR==NR { keys[++n]=$1; next }
      FILENAME==ARGV[2] { a[$1]=$2; next }
      FILENAME==ARGV[3] { b[$1]=$2; next }
      END {
        for(i=1;i<=n;i++){
          k=keys[i];
          if((k in a) && (k in b)) print a[k], b[k];
        }
      }' "${tmpdir}/common_all.keys" "${MAP[$anchor_base]}" "${MAP[$obase]}" > "$pairfile"
  done

  # ---------------- Optional: 3-column triplets lists (aligned to common_all) ----------------
  if [[ -n "${MAP[G4Hits]:-}" && -n "${MAP[DST_JETS]:-}" ]]; then
    local tripA="${outdir}/${anchor_base}__G4Hits__DST_JETS.triplets.list"
    awk '
      FNR==NR { keys[++n]=$1; next }
      FILENAME==ARGV[2] { a[$1]=$2; next }
      FILENAME==ARGV[3] { b[$1]=$2; next }
      FILENAME==ARGV[4] { c[$1]=$2; next }
      END {
        for(i=1;i<=n;i++){
          k=keys[i];
          if((k in a) && (k in b) && (k in c)) print a[k], b[k], c[k];
        }
      }' "${tmpdir}/common_all.keys" "${MAP[$anchor_base]}" "${MAP[G4Hits]}" "${MAP[DST_JETS]}" > "$tripA"
  fi

  if [[ -n "${MAP[DST_GLOBAL]:-}" && -n "${MAP[DST_MBD_EPD]:-}" ]]; then
    local tripB="${outdir}/${anchor_base}__DST_GLOBAL__DST_MBD_EPD.triplets.list"
    awk '
      FNR==NR { keys[++n]=$1; next }
      FILENAME==ARGV[2] { a[$1]=$2; next }
      FILENAME==ARGV[3] { b[$1]=$2; next }
      FILENAME==ARGV[4] { c[$1]=$2; next }
      END {
        for(i=1;i<=n;i++){
          k=keys[i];
          if((k in a) && (k in b) && (k in c)) print a[k], b[k], c[k];
        }
      }' "${tmpdir}/common_all.keys" "${MAP[$anchor_base]}" "${MAP[DST_GLOBAL]}" "${MAP[DST_MBD_EPD]}" > "$tripB"
  fi

  # ---------------- Report what we produced (FILE ONLY; keep it clean) ----------------
  {
    echo ""
    echo "MATCHED OUTPUTS (common across ALL lists):"
    printf "  %s (%d)\n" "$(basename "$anchor_matched")" "$(wc -l < "$anchor_matched" | tr -d ' ')"
    for oname in "${others[@]}"; do
      local obase="${oname%.list}"
      local other_matched="${outdir}/${obase}.matched.list"
      printf "  %s (%d)\n" "$(basename "$other_matched")" "$(wc -l < "$other_matched" | tr -d ' ')"
    done

    echo "PAIR/TRIPLET LISTS:"
    for oname in "${others[@]}"; do
      local obase="${oname%.list}"
      local pairfile="${outdir}/${anchor_base}__${obase}.pairs.list"
      printf "  %s (%d)\n" "$(basename "$pairfile")" "$(wc -l < "$pairfile" | tr -d ' ')"
    done

    if [[ -f "${outdir}/${anchor_base}__G4Hits__DST_JETS.triplets.list" ]]; then
      printf "  %s (%d)\n" "$(basename "${outdir}/${anchor_base}__G4Hits__DST_JETS.triplets.list")" \
        "$(wc -l < "${outdir}/${anchor_base}__G4Hits__DST_JETS.triplets.list" | tr -d ' ')"
    fi
    if [[ -f "${outdir}/${anchor_base}__DST_GLOBAL__DST_MBD_EPD.triplets.list" ]]; then
      printf "  %s (%d)\n" "$(basename "${outdir}/${anchor_base}__DST_GLOBAL__DST_MBD_EPD.triplets.list")" \
        "$(wc -l < "${outdir}/${anchor_base}__DST_GLOBAL__DST_MBD_EPD.triplets.list" | tr -d ' ')"
    fi
  } >> "$report"
}


preview_head() {
  local f="$1" k="$2"
  (( k > 0 )) || return 0
  [[ -s "$f" ]] || return 0

  if [[ "$DEBUG" == "true" ]]; then
    printf "%b  preview (first %d):%b\n" "$C_DIM" "$k" "$C_RESET" >&2
    head -n "$k" "$f" 2>/dev/null | sed 's/^/    /' >&2 || true
  else
    printf "%b  preview (first %d basenames):%b\n" "$C_DIM" "$k" "$C_RESET" >&2
    head -n "$k" "$f" 2>/dev/null | awk -F/ '{print "    " $NF}' >&2 || true
  fi
}

build_pack() {
  # Args: <sample_dirname>  (photonjet5 | photonjet10 | photonjet20)
  local sample="$1"

  local tag="run${RUNNUM}_${sample}"
  local outdir="${OUTROOT}/${tag}"

  step "BEGIN PACK: ${tag}"

  # Start fresh every time (and remember what we removed for the final summary)
  if [[ -d "$outdir" ]]; then
    local n_old
    n_old=$(find "$outdir" -maxdepth 1 -type f 2>/dev/null | wc -l | tr -d ' ')
    PACK_REMOVED["$tag"]="removed existing outdir (files=${n_old})"
    step "Cleanup: removing existing outdir → ${outdir} (files=${n_old})"
    rm -rf "$outdir"
  else
    PACK_REMOVED["$tag"]="none (outdir did not exist)"
    step "Cleanup: no prior outdir to remove → ${outdir}"
  fi

  step "Create: output directory → ${outdir}"
  mkdir -p "$outdir"

  # Source directories (exactly as on lustre)
  local g4dir="${MDC2_BASE}/g4hits/run00${RUNNUM}/${sample}"
  local calodir="${MDC2_BASE}/nopileup/calocluster/run00${RUNNUM}/${sample}"
  local gldir="${MDC2_BASE}/nopileup/global/run00${RUNNUM}/${sample}"
  local jetsdir="${MDC2_BASE}/nopileup/jets/run00${RUNNUM}/${sample}"
  local mbddir="${MDC2_BASE}/nopileup/mbdepd/run00${RUNNUM}/${sample}"
  local trkdir="${MDC2_BASE}/nopileup/trkrhit/run00${RUNNUM}/${sample}"

  rule
  say "PACK: ${tag}"
  say "OUT : ${outdir}"
  say "SRC :"
  say "  CaloClus   = ${calodir}   (anchor)"
  say "  G4Hits     = ${g4dir}"
  say "  Jets       = ${jetsdir}"
  rule

  # Helper: write a list from a directory (root files only)
  write_list() {
    local src="$1"
    local out="$2"
    local label="$3"

    if [[ ! -d "$src" ]]; then
      die "Missing directory for ${label}: ${src}"
    fi

    find "$src" \( -type f -o -type l \) -iname '*.root' -print 2>/dev/null | sort -V > "$out"
    local n; n=$(wc -l < "$out" | tr -d ' ')
    if (( n == 0 )); then
      echo "[ERR  $(ts)] 0 ROOT files found for ${label} in: ${src}" >&2
      echo "[ERR  $(ts)] Debug: listing directory (first 30 entries):" >&2
      ls -la "$src" 2>/dev/null | head -n 30 >&2 || true
      echo "[ERR  $(ts)] Debug: any symlinks here? (first 30):" >&2
      find "$src" -maxdepth 1 -type l 2>/dev/null | head -n 30 >&2 || true
      die "Stopping because list would be empty: ${out}"
    fi


    ok "Wrote $(printf '%6d' "$n") -> $(basename "$out")"
    dbg "Path: $out"
    preview_head "$out" "$HEAD_TAIL_LINES"
  }

  step "Write raw lists (filesystem scan → *.list)"
  write_list "$calodir" "${outdir}/DST_CALO_CLUSTER.list" "DST_CALO_CLUSTER (calocluster) [ANCHOR]"
  write_list "$jetsdir" "${outdir}/DST_JETS.list"        "DST_JETS (nopileup/jets)"
  write_list "$gldir"   "${outdir}/DST_GLOBAL.list"      "DST_GLOBAL (nopileup/global) [VERTEX]"
  write_list "$mbddir"  "${outdir}/DST_MBD_EPD.list"     "DST_MBD_EPD (nopileup/mbdepd) [VERTEX]"

  # --------------------------- G4Hits (optional) ---------------------------
  # If missing/empty, warn in RED and create a placeholder list of "NONE"
  # with the same number of lines as the anchor list so downstream stays aligned.
  {
    if [[ -d "$g4dir" ]]; then
      find "$g4dir" \( -type f -o -type l \) -iname '*.root' -print 2>/dev/null | sort -V > "${outdir}/G4Hits.list"
    else
      : > "${outdir}/G4Hits.list"
    fi
  } || true

  n_g4="$(wc -l < "${outdir}/G4Hits.list" 2>/dev/null | tr -d ' ' || echo 0)"
  n_calo="$(wc -l < "${outdir}/DST_CALO_CLUSTER.list" 2>/dev/null | tr -d ' ' || echo 0)"

  G4_OK="true"
  if [[ "${n_g4}" == "0" ]]; then
    G4_OK="false"
    # RED warning
    printf "\033[31m[WARN %s] G4Hits missing/empty for %s (dir=%s). Continuing WITHOUT real G4Hits.\033[0m\n" \
      "$(ts)" "${tag}" "${g4dir}" >&2
    printf "\033[31m[WARN %s] Writing placeholder G4Hits.list filled with 'NONE' (lines=%d) so 5-col lists can still be built.\033[0m\n" \
      "$(ts)" "${n_calo}" >&2

    awk '{print "NONE"}' "${outdir}/DST_CALO_CLUSTER.list" > "${outdir}/G4Hits.list"
    n_g4="$(wc -l < "${outdir}/G4Hits.list" | tr -d ' ')"
  else
    ok "Wrote $(printf '%6d' "$n_g4") -> G4Hits.list"
    dbg "Path: ${outdir}/G4Hits.list"
    preview_head "${outdir}/G4Hits.list" "$HEAD_TAIL_LINES"
  fi

  # --------------------------- Pairing (include G4 if present) ---------------------------
  # Rule:
  #   - If we found >=1 real G4Hits file, G4Hits becomes REQUIRED for pairing and must produce nonzero matches.
  #   - Only run "NO G4" mode if truly NONE exist (then we will generate placeholder 'NONE' aligned lists later).
  if [[ "$G4_OK" == "true" ]]; then
    step "Pairing check + generate matched/pairs/triplets (key-aligned) [WITH G4]"
    pair_check "$outdir" "DST_CALO_CLUSTER.list" \
      "G4Hits.list" "DST_JETS.list" "DST_GLOBAL.list" "DST_MBD_EPD.list"
  else
    step "Pairing check + generate matched/pairs/triplets (key-aligned) [NO G4]"
    pair_check "$outdir" "DST_CALO_CLUSTER.list" \
      "DST_JETS.list" "DST_GLOBAL.list" "DST_MBD_EPD.list"
  fi

  # --------------------------- Placeholder matched/pairs for G4 (if missing) ---------------------------
  # Create key-aligned placeholder matched/pairs/triplets so downstream tooling that expects these files
  # doesn't explode, while still clearly encoding "NONE".
  if [[ "$G4_OK" == "false" ]]; then
    awk '{print "NONE"}' "${outdir}/DST_CALO_CLUSTER.matched.list" > "${outdir}/G4Hits.matched.list"
    paste "${outdir}/DST_CALO_CLUSTER.matched.list" "${outdir}/G4Hits.matched.list" > "${outdir}/DST_CALO_CLUSTER__G4Hits.pairs.list"
    paste "${outdir}/DST_CALO_CLUSTER.matched.list" "${outdir}/G4Hits.matched.list" "${outdir}/DST_JETS.matched.list" > "${outdir}/DST_CALO_CLUSTER__G4Hits__DST_JETS.triplets.list"

    # Also append a note to the pairing report
    {
      echo ""
      echo "NOTE: G4Hits missing/empty for this pack; wrote placeholder 'NONE' lists."
      echo "  - G4Hits.list and G4Hits.matched.list are placeholders"
      echo "  - DST_CALO_CLUSTER__G4Hits.pairs.list and ...triplets.list are placeholders"
    } >> "${outdir}/pair_report.txt" 2>/dev/null || true
  fi

  step "Write summary.txt (counts + pointers)"
  {
    echo "================ PACK SUMMARY ================"
    echo "Tag   : ${tag}"
    echo "Outdir: ${outdir}"
    echo ""

    echo "-- Raw list counts (filesystem scan) --"
    for f in DST_CALO_CLUSTER G4Hits DST_JETS; do
      printf "%-26s : %6d  (%s)\n" \
        "${f}.list" \
        "$(wc -l < "${outdir}/${f}.list" | tr -d ' ')" \
        "${outdir}/${f}.list"
    done
    echo ""

    echo "-- Matched counts (intersection across ALL lists) --"
    for f in DST_CALO_CLUSTER G4Hits DST_JETS; do
      printf "%-26s : %6d  (%s)\n" \
        "${f}.matched.list" \
        "$(wc -l < "${outdir}/${f}.matched.list" | tr -d ' ')" \
        "${outdir}/${f}.matched.list"
    done
    echo ""

    echo "-- Pair / triplet lists (key-aligned) --"
    printf "%-26s : %6d  (%s)\n" \
      "DST_CALO_CLUSTER__G4Hits.pairs.list" \
      "$(wc -l < "${outdir}/DST_CALO_CLUSTER__G4Hits.pairs.list" | tr -d ' ')" \
      "${outdir}/DST_CALO_CLUSTER__G4Hits.pairs.list"

    printf "%-26s : %6d  (%s)\n" \
      "DST_CALO_CLUSTER__DST_JETS.pairs.list" \
      "$(wc -l < "${outdir}/DST_CALO_CLUSTER__DST_JETS.pairs.list" | tr -d ' ')" \
      "${outdir}/DST_CALO_CLUSTER__DST_JETS.pairs.list"

    if [[ -f "${outdir}/DST_CALO_CLUSTER__G4Hits__DST_JETS.triplets.list" ]]; then
      printf "%-26s : %6d  (%s)\n" \
        "DST_CALO_CLUSTER__G4Hits__DST_JETS.triplets.list" \
        "$(wc -l < "${outdir}/DST_CALO_CLUSTER__G4Hits__DST_JETS.triplets.list" | tr -d ' ')" \
        "${outdir}/DST_CALO_CLUSTER__G4Hits__DST_JETS.triplets.list"
    fi

    echo ""
    echo "-- Diagnostics --"
    echo "pair_report.txt : ${outdir}/pair_report.txt"
    echo "summary.txt     : ${outdir}/summary.txt"
    echo "dropped lists   : ${outdir}/dropped_from_*.list"
    echo "missing keys    : ${outdir}/missing_in_*.keys (only if any missing)"
    echo "============================================="
  } | tee "${outdir}/summary.txt"

  say "Final recommended inputs:"
  if [[ "${G4_OK:-true}" == "true" ]]; then
    say "  (1) Calo + G4Hits paired list:"
    say "      ${outdir}/DST_CALO_CLUSTER__G4Hits.pairs.list"
  else
    warn "$(printf "\033[31mG4Hits missing for this pack; G4 pairs/triplets are placeholders ('NONE').\033[0m")"
    say "  (1) Calo + Jets paired list (real):"
    say "      ${outdir}/DST_CALO_CLUSTER__DST_JETS.pairs.list"
  fi
  say "  (2) Matched single-column lists (same key-order):"
  say "      ${outdir}/DST_CALO_CLUSTER.matched.list"
  say "      ${outdir}/G4Hits.matched.list"
  say "      ${outdir}/DST_JETS.matched.list"
  say "  (3) Optional convenience:"
  say "      ${outdir}/DST_CALO_CLUSTER__G4Hits__DST_JETS.triplets.list"

  record_pack "$tag" "$outdir"
  step "END PACK: ${tag}"

  rule
  echo "[DONE $(ts)] Lists & diagnostics ready → ${outdir}" >&2
}

# --------------------------- Photon+Jet type discovery ------------------------

# Probe one (type, run) for a single file and return the first path line (or empty)
probe_first_path() {
  # Args: <type> <probe_log_file>
  local typ="$1"
  local plog="${2:-}"
  local CF; CF="$(resolve_catalog)"
  local tmp; tmp="$(mktemp -d)" || return 1

  local -a cmd=( "$CF" -type "$typ" -run "$RUNNUM" "${CATALOG_EXTRA_ARGS[@]}" -n 1 "$PROBE_FILETYPE" )

  dbg "PROBE type=${typ}  cmd=${cmd[*]}"

  # Run catalog probe and capture stderr/stdout for transparency (but don't kill the scan)
  local out ec
  out="$( (cd "$tmp" && "${cmd[@]}") 2>&1 )"
  ec=$?

  if [[ -n "$plog" ]]; then
    {
      echo "[$(ts)] type=${typ} exit=${ec} cmd=${cmd[*]}"
      if (( ec != 0 )); then
        echo "[$(ts)] type=${typ} catalog output (first 40 lines):"
        echo "$out" | head -n 40
        echo "----"
      fi
    } >> "$plog"
  fi

  # Locate the produced list (some tools write uppercase, some lowercase)
  local upper="${PROBE_FILETYPE}.list"
  local lower; lower="$(type_to_lower_list "$PROBE_FILETYPE")"

  local f=""
  if [[ -s "$tmp/$upper" ]]; then f="$tmp/$upper"; fi
  if [[ -z "$f" && -s "$tmp/$lower" ]]; then f="$tmp/$lower"; fi

  local line=""
  if [[ -n "$f" ]]; then
    line="$(head -n 1 "$f" || true)"
  else
    dbg "PROBE type=${typ} produced no list (${upper} or ${lower})"
  fi

  rm -rf "$tmp"
  echo "$line"
}

# Case-insensitive match: does a path look like Photon5 / Photon10 / Photon20?
path_matches_photon_sample() {
  local path="$1"
  local sample="$2"

  local p; p="$(tr 'A-Z' 'a-z' <<<"$path")"
  local e=""
  case "$sample" in
    Photon5)  e="5"  ;;
    Photon10) e="10" ;;
    Photon20) e="20" ;;
    *) return 1 ;;
  esac

  # Match "photon" then any non-digit run of chars, then the energy, then a non-digit boundary
  # Examples that will match: photon5, photon_5, photonjet_10, photon-20gev, etc.
  [[ "$p" =~ photon[^0-9]*${e}([^0-9]|$) ]]
}

discover_photon_types() {
  # Outputs lines (to STDOUT only): "<Sample> <Type> <ExamplePath>"
  # All logs go to STDERR via say/warn/dbg.

  local candidates="${OUTROOT}/run${RUNNUM}_candidate_types_${PROBE_FILETYPE}.txt"
  local probe_log="${OUTROOT}/run${RUNNUM}_probe_${PROBE_FILETYPE}.log"
  : > "$candidates"
  : > "$probe_log"

  say "Scanning catalog types ${TYPE_SCAN_MIN}..${TYPE_SCAN_MAX} for Run=${RUNNUM} (probe=${PROBE_FILETYPE})"
  say "Catalog tool: $(command -v "$(resolve_catalog)" || true)"
  say "Probe log     → ${probe_log}"
  say "Candidate list→ ${candidates}"

  declare -A first_by_type=()
  local t line found=0
  local every=5   # progress cadence; reduce to 1 if you want every type printed

  for t in $(seq "$TYPE_SCAN_MIN" "$TYPE_SCAN_MAX"); do
    line="$(probe_first_path "$t" "$probe_log")"

    if [[ -n "$line" ]]; then
      first_by_type["$t"]="$line"
      printf "type=%s  %s\n" "$t" "$line" >> "$candidates"
      found=$((found+1))
      dbg "FOUND: type=${t} -> ${line}"
    else
      dbg "MISS : type=${t}"
    fi

    if (( t % every == 0 )); then
      say "Progress: scanned type=${t}/${TYPE_SCAN_MAX}  candidates_found=${found}"
    fi
  done

  if [[ ! -s "$candidates" ]]; then
    die "No candidates found for Run=${RUNNUM} with probe=${PROBE_FILETYPE}. Try PROBE_FILETYPE=DST_CALO_CLUSTER or expand TYPE_SCAN_MAX. Probe log: ${probe_log}"
  fi

  say "Scan complete: candidates_found=${found}"
  say "Quick grep photon (first 20):"
  grep -i "photon" "$candidates" | head -n 20 >&2 || true

  local found_any="false"
  local s
  for s in "${PHOTON_SAMPLES[@]}"; do
    local found_type=""
    local best_path=""

    # deterministically scan types in ascending order
    local t_sorted
    for t_sorted in $(printf '%s\n' "${!first_by_type[@]}" | sort -n); do
      line="${first_by_type[$t_sorted]}"
      if path_matches_photon_sample "$line" "$s"; then
        found_type="$t_sorted"
        best_path="$line"
        break
      fi
    done

    if [[ -n "$found_type" ]]; then
      say "MATCH: ${s} -> type=${found_type}  example=$(basename "$best_path")"
      found_any="true"
      echo "$s $found_type $best_path"
    else
      warn "NO MATCH for ${s}. Inspect ${candidates} and ${probe_log}."
    fi
  done

  [[ "$found_any" == "true" ]] || return 1
  return 0
}

# --------------------------- Main ---------------------------------------------
mkdir -p "$OUTROOT"

# Hard reset: always start fresh for this Run's photonjet packs (+ any stale run-level logs)
clean_glob="${OUTROOT}/run${RUNNUM}_photonjet"*
if compgen -G "$clean_glob" >/dev/null 2>&1; then
  step "Clean: removing previous outputs -> ${clean_glob}"
  rm -rf $clean_glob 2>/dev/null || true
else
  step "Clean: no prior outputs found for run${RUNNUM}_photonjet*"
fi

rm -f "${OUTROOT}/run${RUNNUM}_candidate_types_"*.txt "${OUTROOT}/run${RUNNUM}_probe_"*.log 2>/dev/null || true

rule
say "Building Run-${RUNNUM} Pythia photon+jet lists under: ${OUTROOT}"
say "Mode: filesystem scan (no catalog)"
say "Base: ${MDC2_BASE}"
say "Samples: ${PHOTONJET_SAMPLES[*]}"
rule

say "Building Run-${RUNNUM} photon+jet lists from filesystem (no catalog)…"
say "Base: ${MDC2_BASE}"
say "Samples: ${PHOTONJET_SAMPLES[*]}"

for sample in "${PHOTONJET_SAMPLES[@]}"; do
  build_pack "$sample"
done

final_summary
say "All packs completed. Root: ${OUTROOT}"

