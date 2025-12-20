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
PHOTONJET_SAMPLES=( "photonjet5" "photonjet10" "photonjet20" )

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

# IMPORTANT: send INFO to STDERR so it doesn't get captured by mapfile/process-substitution.
say()  { if [[ "$VERBOSE" == "true" ]]; then echo "[INFO $(ts)] $*" >&2; fi; }
dbg()  { if [[ "$DEBUG" == "true" ]]; then echo "[DBG  $(ts)] $*" >&2; fi; }
warn() { echo "[WARN $(ts)] $*" >&2; }
die()  { echo "[ERR  $(ts)] $*" >&2; exit 1; }
rule() { echo "-------------------------------------------------------------------------------" >&2; }

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

  say "Pairing check (anchor=${anchor}) → ${report}"

  local ak="${outdir}/${anchor}"
  [[ -s "$ak" ]] || { warn "Anchor list missing or empty: ${anchor}"; return 0; }

  # Temp dir (function-local cleanup; no global EXIT trap)
  local tmpdir; tmpdir="$(mktemp -d)"
  trap '[[ -n "${tmpdir:-}" ]] && rm -rf "${tmpdir}"' RETURN

  # Build anchor keys and a key->path map
  awk '
    function key_of(path,    b, n, a, last, prev) {
      b=path; sub(/^.*\//,"",b); sub(/\.root$/,"",b);
      n=split(b,a,"-");
      if(n>=2 && a[n] ~ /^[0-9]+$/ && a[n-1] ~ /^[0-9]+$/) return a[n-1] "-" a[n];
      return b;
    }
    { k=key_of($0); print k }' "$ak" | sort -u > "${tmpdir}/anchor.keys"

  awk '
    function key_of(path,    b, n, a, last, prev) {
      b=path; sub(/^.*\//,"",b); sub(/\.root$/,"",b);
      n=split(b,a,"-");
      if(n>=2 && a[n] ~ /^[0-9]+$/ && a[n-1] ~ /^[0-9]+$/) return a[n-1] "-" a[n];
      return b;
    }
    { k=key_of($0); print k "\t" $0 }' "$ak" > "${tmpdir}/anchor.map"

  local a_count; a_count=$(wc -l < "${tmpdir}/anchor.keys" | tr -d ' ')
  printf "ANCHOR: %s  keys=%d\n" "$anchor" "$a_count" | tee -a "$report"

  # For each other list, compute intersection and (optionally) write matched lists
  for oname in "${others[@]}"; do
    local ok="${outdir}/${oname}"
    if [[ ! -s "$ok" ]]; then
      printf "LIST  : %s  (missing/empty)\n" "$oname" | tee -a "$report"
      continue
    fi

    # other keys and map
    awk '
      function key_of(path,    b, n, a, last, prev) {
        b=path; sub(/^.*\//,"",b); sub(/\.root$/,"",b);
        n=split(b,a,"-");
        if(n>=2 && a[n] ~ /^[0-9]+$/ && a[n-1] ~ /^[0-9]+$/) return a[n-1] "-" a[n];
        return b;
      }
      { k=key_of($0); print k }' "$ok" | sort -u > "${tmpdir}/other.keys"

    awk '
      function key_of(path,    b, n, a, last, prev) {
        b=path; sub(/^.*\//,"",b); sub(/\.root$/,"",b);
        n=split(b,a,"-");
        if(n>=2 && a[n] ~ /^[0-9]+$/ && a[n-1] ~ /^[0-9]+$/) return a[n-1] "-" a[n];
        return b;
      }
      { k=key_of($0); print k "\t" $0 }' "$ok" > "${tmpdir}/other.map"

    # Intersections and diffs
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

    # Write detailed “what’s missing” key list + sample
    if (( m_count > 0 )); then
      local miss_keys="${outdir}/missing_in_${oname%.list}.keys"
      cp -f "${tmpdir}/missing_in_other.keys" "$miss_keys"
      printf "  missing-keys file → %s\n" "$miss_keys" | tee -a "$report"
      printf "  sample missing → " | tee -a "$report"
      head -n 8 "${tmpdir}/missing_in_other.keys" | paste -sd',' - | tee -a "$report"
      echo | tee -a "$report"
    fi

    # --- New behavior: produce matched-only lists for this pair ---
    # Anchor matched list
    local anchor_matched="${outdir}/${anchor%.list}.matched.list"
    awk 'NR==FNR{keep[$1]=1; next} ($1 in keep){print $2}' \
      "${tmpdir}/common.keys" "${tmpdir}/anchor.map" | sort -V > "$anchor_matched"

    # Other matched list
    local other_matched="${outdir}/${oname%.list}.matched.list"
    awk 'NR==FNR{keep[$1]=1; next} ($1 in keep){print $2}' \
      "${tmpdir}/common.keys" "${tmpdir}/other.map" | sort -V > "$other_matched"

    # Dropped-from-anchor paths (those anchor keys missing in other)
    local dropped_anchor="${outdir}/dropped_from_${anchor%.list}_vs_${oname%.list}.list"
    awk 'NR==FNR{drop[$1]=1; next} ($1 in drop){print $2}' \
      "${tmpdir}/missing_in_other.keys" "${tmpdir}/anchor.map" | sort -V > "$dropped_anchor"

    local n_anchor_matched n_other_matched n_dropped
    n_anchor_matched=$(wc -l < "$anchor_matched" | tr -d ' ')
    n_other_matched=$(wc -l < "$other_matched" | tr -d ' ')
    n_dropped=$(wc -l < "$dropped_anchor" | tr -d ' ')

    printf "  MATCHED OUTPUTS:\n" | tee -a "$report"
    printf "    %s  (%d)\n" "$(basename "$anchor_matched")" "$n_anchor_matched" | tee -a "$report"
    printf "    %s  (%d)\n" "$(basename "$other_matched")" "$n_other_matched" | tee -a "$report"
    printf "  DROPPED FROM ANCHOR (because missing in %s): %s  (%d)\n" \
      "$oname" "$(basename "$dropped_anchor")" "$n_dropped" | tee -a "$report"
  done
}

preview_head() {
  local f="$1" k="$2"
  if (( k > 0 )); then
    head -n "$k" "$f" || true
  fi
}

build_pack() {
  # Args: <sample_dirname>  (photonjet5 | photonjet10 | photonjet20)
  local sample="$1"

  local tag="run${RUNNUM}_${sample}"
  local outdir="${OUTROOT}/${tag}"

  # Start fresh every time
  if [[ -d "$outdir" ]]; then
    rm -rf "$outdir"
  fi
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

    find "$src" -type f -name '*.root' | sort -V > "$out"
    local n; n=$(wc -l < "$out" | tr -d ' ')
    if (( n == 0 )); then
      die "0 ROOT files found for ${label} in: ${src}"
    fi

    printf "[INFO %s] Wrote %6d → %s\n" "$(ts)" "$n" "$out" >&2
    preview_head "$out" "$HEAD_TAIL_LINES"
  }

  # Raw lists (everything that exists on lustre)
  write_list "$calodir" "${outdir}/DST_CALO_CLUSTER.list" "DST_CALO_CLUSTER (calocluster) [ANCHOR]"
  write_list "$g4dir"   "${outdir}/G4Hits.list"          "G4Hits (g4hits)"
  write_list "$jetsdir" "${outdir}/DST_JETS.list"        "DST_JETS (nopileup/jets)"

  # Match everything to the anchor (calocluster) and write *.matched.list outputs
  pair_check "$outdir" "DST_CALO_CLUSTER.list" \
    "G4Hits.list" "DST_JETS.list"

  # Summary counts: raw + matched
  {
    echo "================ PACK SUMMARY ================"
    echo "Tag: ${tag}"
    echo ""
    echo "-- Raw list counts --"
    for f in DST_CALO_CLUSTER G4Hits DST_JETS; do
      printf "%-16s : %6d\n" "${f}" "$(wc -l < "${outdir}/${f}.list" | tr -d ' ')"
    done
    echo ""
    echo "-- Matched-to-anchor counts (recommended) --"
    for f in DST_CALO_CLUSTER G4Hits DST_JETS; do
      printf "%-24s : %6d\n" "${f}.matched" "$(wc -l < "${outdir}/${f}.matched.list" | tr -d ' ')"
    done
    echo "============================================="
  } | tee "${outdir}/summary.txt"

  say "Final recommended inputs (matched to calocluster):"
  say "  ${outdir}/DST_CALO_CLUSTER.matched.list"
  say "  ${outdir}/G4Hits.matched.list"
  say "  ${outdir}/DST_JETS.matched.list"

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

rule
say "All packs completed. Root: ${OUTROOT}"

