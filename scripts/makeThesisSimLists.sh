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

set -euo pipefail

# --------------------------- Defaults -----------------------------------------
OUTROOT="/sphenix/u/patsfan753/scratch/thesisAnalysis/simFiles"
HEAD_TAIL_LINES="8"
VERBOSE="true"

# Fixed parameters for this script
RUNNUM="23"
TYPES=( "11" "12" )   # 11 = 30 GeV, 12 = 10 GeV
LABEL_11="30GeV"
LABEL_12="10GeV"

# Requested filetypes (order matters only for presentation)
# Note: G4Hits is not supported with -embed; omit it here.
FILETYPES=(DST_CALO_CLUSTER DST_CALO_G4HIT DST_GLOBAL DST_MBD_EPD DST_TRUTH_G4HIT DST_TRUTH_JET)


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
ts()   { date '+%Y-%m-%d %H:%M:%S'; }
say()  { if [[ "$VERBOSE" == "true" ]]; then echo "[INFO $(ts)] $*"; fi; }
warn() { echo "[WARN $(ts)] $*" >&2; }
die()  { echo "[ERR  $(ts)] $*"  >&2; exit 1; }
rule() { echo "-------------------------------------------------------------------------------"; }

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

  say "Pairing check (anchor=${anchor}) → ${report}"

  local ak="${outdir}/${anchor}"
  [[ -s "$ak" ]] || { warn "Anchor list missing or empty: ${anchor}"; return 0; }

  local tmpdir; tmpdir="$(mktemp -d)"; trap 'rm -rf "$tmpdir"' EXIT
  extract_keys "$ak" > "${tmpdir}/anchor.keys"

  local a_count; a_count=$(wc -l < "${tmpdir}/anchor.keys" | tr -d ' ')
  printf "ANCHOR: %s  keys=%d\n" "$anchor" "$a_count" | tee -a "$report"

  for oname in "${others[@]}"; do
    local ok="${outdir}/${oname}"
    if [[ ! -s "$ok" ]]; then
      printf "LIST  : %s  (missing/empty)\n" "$oname" | tee -a "$report"
      continue
    fi
    extract_keys "$ok" > "${tmpdir}/other.keys"

    # Compute intersections and differences
    comm -12 "${tmpdir}/anchor.keys" "${tmpdir}/other.keys" > "${tmpdir}/common.keys"
    comm -23 "${tmpdir}/anchor.keys" "${tmpdir}/other.keys" > "${tmpdir}/missing_in_other.keys"
    comm -13 "${tmpdir}/anchor.keys" "${tmpdir}/other.keys" > "${tmpdir}/extra_in_other.keys"

    local c_count m_count e_count
    c_count=$(wc -l < "${tmpdir}/common.keys" | tr -d ' ')
    m_count=$(wc -l < "${tmpdir}/missing_in_other.keys" | tr -d ' ')
    e_count=$(wc -l < "${tmpdir}/extra_in_other.keys" | tr -d ' ')

    printf "LIST  : %s  common=%d / %d (%.1f%%)  missing=%d  extra=%d\n" \
      "$oname" "$c_count" "$a_count" "$(awk -v c="$c_count" -v a="$a_count" 'BEGIN{if(a>0) printf("%.1f",100*c/a); else print "0.0"}')" \
      "$m_count" "$e_count" | tee -a "$report"

    if (( m_count > 0 )); then
      printf "  sample missing → " | tee -a "$report"
      head -n 5 "${tmpdir}/missing_in_other.keys" | paste -sd',' - | tee -a "$report"
      echo | tee -a "$report"
    fi
  done
}

preview_head() {
  local f="$1" k="$2"
  if (( k > 0 )); then
    head -n "$k" "$f" || true
  fi
}

build_pack() {
  # Args: <type> <energy_label>  (type=11→30GeV, 12→10GeV)
  local typ="$1" elabel="$2"
  local CF; CF="$(resolve_catalog)"

  local tag="run${RUNNUM}_type${typ}_auau"
  local outdir="${OUTROOT}/${tag}"
  mkdir -p "$outdir"

  rule
  say "PACK: ${tag}  (embed=auau, Run=${RUNNUM}, Type=${typ} ~ ${elabel})"
  say "TOOL: $(command -v "$CF")"
  say "OUT : ${outdir}"
  rule

  (
    cd "$outdir"
    # Run the catalog once, requesting all desired filetypes in one shot
    # This writes one *.list per filetype into $outdir
    "$CF" -type "$typ" -run "$RUNNUM" -embed auau "${FILETYPES[@]}" 2>&1 | tee catalog.log
  )

  # Normalize list names to UPPERCASE and validate
  local produced_total=0
  for ft in "${FILETYPES[@]}"; do
    local lower upper nlines
    lower="$(type_to_lower_list "$ft")"
    upper="${ft}.list"
    if [[ -f "${outdir}/${lower}" ]]; then
      mv -f "${outdir}/${lower}" "${outdir}/${upper}"
    fi
    if [[ ! -s "${outdir}/${upper}" ]]; then
      die "Missing or empty list: ${outdir}/${upper} (see ${outdir}/catalog.log)"
    fi
    nlines=$(wc -l < "${outdir}/${upper}" | tr -d ' ')
    produced_total=$(( produced_total + nlines ))
    printf "[INFO %s] Wrote %6d → %s\n" "$(ts)" "$nlines" "${outdir}/${upper}"
    preview_head "${outdir}/${upper}" "$HEAD_TAIL_LINES"
  done

  # Summary counts
  {
    echo "================ PACK SUMMARY ================"
    echo "Tag: ${tag}  (Run ${RUNNUM}, Type ${typ} ~ ${elabel}, embed=auau)"
    for ft in "${FILETYPES[@]}"; do
      printf "%-16s : %6d\n" "${ft}" "$(wc -l < "${outdir}/${ft}.list" | tr -d ' ')"
    done
    echo "Total (sum of lines across lists): ${produced_total}"
    echo "============================================="
  } | tee "${outdir}/summary.txt"

  pair_check "$outdir" "DST_CALO_CLUSTER.list" \
      "DST_CALO_G4HIT.list" "DST_GLOBAL.list" "DST_MBD_EPD.list" \
      "DST_TRUTH_G4HIT.list" "DST_TRUTH_JET.list"

  rule
  echo "[DONE $(ts)] Lists & diagnostics ready → ${outdir}"
}

# --------------------------- Main ---------------------------------------------
mkdir -p "$OUTROOT"

build_pack "11" "$LABEL_11"   # 30 GeV jets
build_pack "12" "$LABEL_12"   # 10 GeV jets

rule
say "All packs completed. Root: ${OUTROOT}"

