#!/bin/bash
# make_dstListsData.sh — unified data DST list builder
#
# Usage:
#   ./make_dstListsData.sh pp24
#   ./make_dstListsData.sh pp25
#   ./make_dstListsData.sh auau
#   ./make_dstListsData.sh oo25
#   ./make_dstListsData.sh run2auau
#   ./make_dstListsData.sh ALL
#
# Optional trigger QA only:
#   ./make_dstListsData.sh pp24 QA
#   ./make_dstListsData.sh pp25 QA
#   ./make_dstListsData.sh auau QA
#   ./make_dstListsData.sh oo25 QA
#   ./make_dstListsData.sh run2auau QA
#   ./make_dstListsData.sh ALL QA

set -uo pipefail  # (intentionally NOT using -e)

IN_BASE="/sphenix/u/patsfan753/scratch/thesisAnalysis"
GRL_BASE="$IN_BASE/GRLs_tanner"
MODE="${1:-}"
ACTION="${2:-}"
EXTRA_ACTION="${3:-}"

usage() {
  cat <<'USAGE'
Usage:
  ./make_dstListsData.sh pp24 [QA]
  ./make_dstListsData.sh pp25 [QA]
  ./make_dstListsData.sh auau [QA]
  ./make_dstListsData.sh auau QA doQUICKcheck
  ./make_dstListsData.sh auau QA scaledTriggerAna
  ./make_dstListsData.sh oo25 [QA]
  ./make_dstListsData.sh run2auau [QA]
  ./make_dstListsData.sh ALL [QA]
USAGE
}

read_runs_from_file() {
  local f="$1"
  local -n out="$2"
  out=()
  while IFS= read -r raw; do
    local rn rn_dec
    rn="${raw%%#*}"
    rn="$(echo -n "$rn" | tr -d '[:space:]')"
    [[ -z "$rn" ]] && continue
    rn_dec=$(printf "%d" "$rn" 2>/dev/null || echo "$rn")
    [[ "$rn_dec" =~ ^[0-9]+$ ]] || continue
    out+=( "$rn_dec" )
  done < "$f"
}

mode_requested_type() {
  if [[ "${BUILD_MODE:-}" == "per_run" ]]; then
    printf '%s\n' "${TYPE:-}"
  else
    printf '%s\n' "${PREFIX:-}"
  fi
}

mode_output_files() {
  if [[ "${BUILD_MODE:-}" == "per_run" ]]; then
    {
      compgen -G "${OUT_DIR}/dst_jetcalo-*.list" || true
      compgen -G "${OUT_DIR}/DST_JETCALO-*.list" || true
    } | sort -u
  else
    local stem
    stem="${PREFIX#DST_}"
    stem="$(echo "$stem" | tr '[:upper:]' '[:lower:]')"
    {
      compgen -G "${OUT_DIR}/dst_${stem}-*.list" || true
      compgen -G "${OUT_DIR}/${PREFIX}-*.list" || true
      compgen -G "${OUT_DIR}/${PREFIX}_${DATASET}_${TAG}-*.list" || true
    } | sort -u
  fi
}

print_naming_scheme_summary() {
  local stem
  echo "Naming scheme summary:"
  if [[ "${BUILD_MODE:-}" == "per_run" ]]; then
    echo "  Primary emitted pattern : dst_jetcalo-<RUN8>.list"
    echo "  Also accepted           : DST_JETCALO-<RUN8>.list"
  else
    stem="${PREFIX#DST_}"
    stem="$(echo "$stem" | tr '[:upper:]' '[:lower:]')"
    echo "  Primary emitted pattern : dst_${stem}-<RUN8>.list"
    echo "  Also accepted           : ${PREFIX}-<RUN8>.list"
    echo "  Also accepted           : ${PREFIX}_${DATASET}_${TAG}-<RUN8>.list"
  fi
}

probe_type_counts() {
  local probe_type="$1"
  local tmpdir rc files entries f probe_stem

  command -v CreateDstList.pl >/dev/null 2>&1 || { printf '0|0|CreateDstList.pl unavailable\n'; return; }
  [[ -n "${DATASET:-}" ]] || { printf '0|0|dataset unavailable\n'; return; }
  [[ -f "$LIST_FILE" ]] || { printf '0|0|missing input list\n'; return; }
  [[ -s "$LIST_FILE" ]] || { printf '0|0|empty input list\n'; return; }

  tmpdir="$(mktemp -d)"
  (
    cd "$tmpdir" || exit 1
    CreateDstList.pl --tag "$TAG" --dataset "$DATASET" --list "$LIST_FILE" "$probe_type" >/dev/null 2>&1
  )
  rc=$?

  if (( rc != 0 )); then
    rm -rf "$tmpdir"
    printf '0|0|rc=%d\n' "$rc"
    return
  fi

  probe_stem="${probe_type#DST_}"
  probe_stem="$(echo "$probe_stem" | tr '[:upper:]' '[:lower:]')"

  mapfile -t _probe_files < <(
    {
      compgen -G "${tmpdir}/dst_${probe_stem}-*.list" || true
      compgen -G "${tmpdir}/${probe_type}-*.list" || true
      compgen -G "${tmpdir}/${probe_type}_${DATASET}_${TAG}-*.list" || true
    } | sort -u
  )
  files=${#_probe_files[@]}
  entries=0
  for f in "${_probe_files[@]}"; do
    entries=$(( entries + $(awk 'NF{n++} END{print n+0}' "$f") ))
  done

  rm -rf "$tmpdir"
  printf '%d|%d|ok\n' "$files" "$entries"
}

probe_other_types_summary() {
  local requested probe_type files entries status
  local -a notes=()

  requested="$(mode_requested_type)"

  for probe_type in DST_CALOFITTING DST_JETCALO; do
    [[ "$probe_type" == "$requested" ]] && continue
    IFS='|' read -r files entries status < <(probe_type_counts "$probe_type")
    if [[ "$status" == "ok" ]]; then
      if (( files > 0 || entries > 0 )); then
        notes+=( "${probe_type}: ${files} lists / ${entries} entries" )
      else
        notes+=( "${probe_type}: none" )
      fi
    else
      notes+=( "${probe_type}: probe failed (${status})" )
    fi
  done

  if ((${#notes[@]} == 0)); then
    printf 'none\n'
    return
  fi

  local IFS='; '
  printf '%s\n' "${notes[*]}"
}

find_run_list_for_type() {
  local list_type="$1"
  local run8="$2"
  local stem
  stem="${list_type#DST_}"
  stem="$(echo "$stem" | tr '[:upper:]' '[:lower:]')"

  local candidates=(
    "${OUT_DIR}/dst_${stem}-${run8}.list"
    "${OUT_DIR}/${list_type}-${run8}.list"
    "${OUT_DIR}/${list_type}_${DATASET}_${TAG}-${run8}.list"
  )

  local f
  for f in "${candidates[@]}"; do
    if [[ -s "$f" ]]; then
      printf '%s\n' "$f"
      return 0
    fi
  done

  return 1
}

dst_runseg_key() {
  local path="$1"
  local base
  base="${path##*/}"

  if [[ "$base" =~ -([0-9]{8})-([0-9]{5})([-.]|\.root$) ]]; then
    printf '%s_%s\n' "${BASH_REMATCH[1]}" "${BASH_REMATCH[2]}"
    return 0
  fi

  if [[ "$base" =~ -([0-9]{8})-([0-9]+)\.root$ ]]; then
    printf '%s_%s\n' "${BASH_REMATCH[1]}" "${BASH_REMATCH[2]}"
    return 0
  fi

  return 1
}

pair_calo_zdc_for_run() {
  local run8="$1"
  local calo_file="$2"
  local zdc_file="$3"
  local out_file="$4"
  local tmp_file="${out_file}.tmp_pair_$$"

  PAIR_MATCHED_SEG=0
  PAIR_CALO_ONLY_SEG=0
  PAIR_ZDC_ONLY_SEG=0
  PAIR_MODE="none"

  local -a calo_lines=()
  local -a zdc_lines=()
  mapfile -t calo_lines < <(awk 'NF{print $1}' "$calo_file")
  mapfile -t zdc_lines  < <(awk 'NF{print $1}' "$zdc_file")

  if (( ${#calo_lines[@]} == 0 )); then
    echo "[WARN] Empty primary list for run ${run8}: ${calo_file}"
    return 1
  fi
  if (( ${#zdc_lines[@]} == 0 )); then
    echo "[WARN] Empty DST_ZDC_RAW list for run ${run8}: ${zdc_file}"
    return 1
  fi

  local use_keyed_pairing=1
  declare -A zdc_by_key=()
  declare -A zdc_seen=()
  declare -A zdc_used=()
  local zdc key calo

  for zdc in "${zdc_lines[@]}"; do
    key="$(dst_runseg_key "$zdc" || true)"
    if [[ -z "$key" || -n "${zdc_by_key[$key]:-}" ]]; then
      use_keyed_pairing=0
      break
    fi
    zdc_by_key["$key"]="$zdc"
    zdc_seen["$key"]=1
  done

  : > "$tmp_file"

  if (( use_keyed_pairing )); then
    for calo in "${calo_lines[@]}"; do
      key="$(dst_runseg_key "$calo" || true)"
      if [[ -z "$key" ]]; then
        use_keyed_pairing=0
        break
      fi

      if [[ -n "${zdc_by_key[$key]:-}" ]]; then
        printf '%s %s\n' "$calo" "${zdc_by_key[$key]}" >> "$tmp_file"
        zdc_used["$key"]=1
        ((PAIR_MATCHED_SEG+=1))
      else
        ((PAIR_CALO_ONLY_SEG+=1))
      fi
    done

    if (( use_keyed_pairing )); then
      for key in "${!zdc_seen[@]}"; do
        if [[ -z "${zdc_used[$key]:-}" ]]; then
          ((PAIR_ZDC_ONLY_SEG+=1))
        fi
      done

      if (( PAIR_MATCHED_SEG == 0 )); then
        echo "[WARN] No matching run/segment pairs for run ${run8}: primary entries=${#calo_lines[@]} DST_ZDC_RAW entries=${#zdc_lines[@]}"
        rm -f "$tmp_file"
        return 1
      fi

      PAIR_MODE="segment_intersection"
      mv -f "$tmp_file" "$out_file"
      return 0
    fi
  fi

  : > "$tmp_file"
  PAIR_MATCHED_SEG=0
  PAIR_CALO_ONLY_SEG=0
  PAIR_ZDC_ONLY_SEG=0

  if (( ${#calo_lines[@]} != ${#zdc_lines[@]} )); then
    echo "[WARN] Cannot pair run ${run8}: primary entries=${#calo_lines[@]} DST_ZDC_RAW entries=${#zdc_lines[@]}"
    rm -f "$tmp_file"
    return 1
  fi

  echo "[WARN] Falling back to line-by-line DST_ZDC_RAW pairing for run ${run8}; could not parse run/segment from every filename."
  local i
  for (( i=0; i<${#calo_lines[@]}; i++ )); do
    printf '%s %s\n' "${calo_lines[$i]}" "${zdc_lines[$i]}" >> "$tmp_file"
    ((PAIR_MATCHED_SEG+=1))
  done

  PAIR_MODE="line_by_line"
  mv -f "$tmp_file" "$out_file"
  return 0
}

runlist_unique_runs() {
  local f="$1"
  [[ -s "$f" ]] || return 0

  awk '
    function emit(x) {
      gsub(/^0+/, "", x);
      if (x == "") x = "0";
      printf "%08d\n", x + 0;
    }
    {
      line = $0;
      sub(/#.*/, "", line);
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", line);
      if (line == "") next;

      split(line, fields, /[[:space:]]+/);
      x = fields[1];

      if (x ~ /^[0-9]+$/) {
        emit(x);
        next;
      }

      if (match(x, /-([0-9]{8})-[0-9]+/, m)) {
        print m[1];
        next;
      }

      if (match(line, /-([0-9]{8})-[0-9]+/, m)) {
        print m[1];
        next;
      }

      if (match(x, /(^|[^0-9])([0-9]{8})([^0-9]|$)/, m)) {
        print m[2];
        next;
      }

      if (match(line, /(^|[^0-9])([0-9]{8})([^0-9]|$)/, m)) {
        print m[2];
        next;
      }

      if (x ~ /^[0-9]{5,7}$/) {
        emit(x);
        next;
      }
    }
  ' "$f" | sort -u
}

run_list_overlap_row() {
  local label_a="$1"
  local file_a="$2"
  local label_b="$3"
  local file_b="$4"

  if [[ ! -s "$file_a" || ! -s "$file_b" ]]; then
    printf "  %-32s | %-32s | %8s | %8s | %8s | %8s | %8s | %10s | %10s\n" \
      "$label_a" "$label_b" "MISSING" "MISSING" "-" "-" "-" "-" "-"
    return 0
  fi

  local tmp_a tmp_b
  tmp_a="$(mktemp)"
  tmp_b="$(mktemp)"

  runlist_unique_runs "$file_a" > "$tmp_a"
  runlist_unique_runs "$file_b" > "$tmp_b"

  awk -v label_a="$label_a" -v label_b="$label_b" '
    FNR == NR {
      if ($1 != "" && !($1 in A)) {
        A[$1] = 1;
        na++;
      }
      next;
    }
    {
      if ($1 != "" && !($1 in B)) {
        B[$1] = 1;
        nb++;
      }
    }
    END {
      overlap = 0;
      for (r in A) {
        if (r in B) overlap++;
      }
      only_a = na - overlap;
      only_b = nb - overlap;
      pct_a = (na > 0 ? 100.0 * overlap / na : 0.0);
      pct_b = (nb > 0 ? 100.0 * overlap / nb : 0.0);
      printf "  %-32s | %-32s | %8d | %8d | %8d | %8d | %8d | %9.2f%% | %9.2f%%\n",
             label_a, label_b, na, nb, overlap, only_a, only_b, pct_a, pct_b;
    }
  ' "$tmp_a" "$tmp_b"

  rm -f "$tmp_a" "$tmp_b"
}

tmp_find_run_list_for_type() {
  local tmpdir="$1"
  local list_type="$2"
  local dataset="$3"
  local tag="$4"
  local run8="$5"
  local stem
  stem="${list_type#DST_}"
  stem="$(echo "$stem" | tr '[:upper:]' '[:lower:]')"

  local candidates=(
    "${tmpdir}/dst_${stem}-${run8}.list"
    "${tmpdir}/${list_type}-${run8}.list"
    "${tmpdir}/${list_type}_${dataset}_${tag}-${run8}.list"
  )

  local f
  for f in "${candidates[@]}"; do
    if [[ -s "$f" ]]; then
      printf '%s\n' "$f"
      return 0
    fi
  done

  return 1
}

segment_key_file() {
  local f="$1"
  awk '
    NF {
      base = $1;
      sub(/^.*\//, "", base);

      if (match(base, /-([0-9]{8})-([0-9]{5})([-.]|\.root$)/, m)) {
        print m[1] "_" m[2];
      } else if (match(base, /-([0-9]{8})-([0-9]+)\.root$/, m)) {
        printf "%s_%05d\n", m[1], m[2] + 0;
      }
    }
  ' "$f" | sort -u
}

run_list_segment_inventory_row() {
  local label="$1"
  local run_file="$2"
  local calo_type="$3"
  local zdc_type="$4"
  local dataset="$5"
  local tag="$6"

  if [[ ! -s "$run_file" ]]; then
    printf "  %-32s | %8s | %8s | %8s | %8s | %10s | %10s | %12s | %10s | %10s\n" \
      "$label" "MISSING" "-" "-" "-" "-" "-" "-" "-" "-"
    return 0
  fi

  local tmpdir norm_runs
  tmpdir="$(mktemp -d)"
  norm_runs="${tmpdir}/runs.list"
  runlist_unique_runs "$run_file" > "$norm_runs"

  local requested_runs
  requested_runs=$(awk 'NF{n++} END{print n+0}' "$norm_runs")

  (
    cd "$tmpdir" || exit 1
    CreateDstList.pl --tag "$tag" --dataset "$dataset" --list "$norm_runs" "$calo_type" >/dev/null 2>&1
    CreateDstList.pl --tag "$tag" --dataset "$dataset" --list "$norm_runs" "$zdc_type"  >/dev/null 2>&1
  )

  local calo_runs=0
  local zdc_runs=0
  local both_runs=0
  local matched_runs=0
  local calo_segments=0
  local zdc_segments=0
  local matched_segments=0
  local calo_only_segments=0
  local zdc_only_segments=0
  local r8 calo_list zdc_list calo_seg_file zdc_seg_file n_calo n_zdc n_match n_calo_only n_zdc_only

  while IFS= read -r r8; do
    [[ -n "$r8" ]] || continue

    calo_list="$(tmp_find_run_list_for_type "$tmpdir" "$calo_type" "$dataset" "$tag" "$r8" || true)"
    zdc_list="$(tmp_find_run_list_for_type "$tmpdir" "$zdc_type" "$dataset" "$tag" "$r8" || true)"

    if [[ -n "$calo_list" ]]; then
      ((calo_runs+=1))
      calo_seg_file="${tmpdir}/calo_${r8}.seg"
      segment_key_file "$calo_list" > "$calo_seg_file"
      n_calo=$(awk 'NF{n++} END{print n+0}' "$calo_seg_file")
      calo_segments=$(( calo_segments + n_calo ))
    else
      calo_seg_file=""
      n_calo=0
    fi

    if [[ -n "$zdc_list" ]]; then
      ((zdc_runs+=1))
      zdc_seg_file="${tmpdir}/zdc_${r8}.seg"
      segment_key_file "$zdc_list" > "$zdc_seg_file"
      n_zdc=$(awk 'NF{n++} END{print n+0}' "$zdc_seg_file")
      zdc_segments=$(( zdc_segments + n_zdc ))
    else
      zdc_seg_file=""
      n_zdc=0
    fi

    if [[ -n "$calo_seg_file" && -n "$zdc_seg_file" ]]; then
      ((both_runs+=1))

      n_match=$(comm -12 "$calo_seg_file" "$zdc_seg_file" | awk 'NF{n++} END{print n+0}')
      n_calo_only=$(comm -23 "$calo_seg_file" "$zdc_seg_file" | awk 'NF{n++} END{print n+0}')
      n_zdc_only=$(comm -13 "$calo_seg_file" "$zdc_seg_file" | awk 'NF{n++} END{print n+0}')

      if (( n_match > 0 )); then
        ((matched_runs+=1))
      fi

      matched_segments=$(( matched_segments + n_match ))
      calo_only_segments=$(( calo_only_segments + n_calo_only ))
      zdc_only_segments=$(( zdc_only_segments + n_zdc_only ))
    fi
  done < "$norm_runs"

  printf "  %-32s | %8d | %8d | %8d | %8d | %10d | %10d | %12d | %10d | %10d\n" \
    "$label" "$requested_runs" "$calo_runs" "$zdc_runs" "$matched_runs" "$calo_segments" "$zdc_segments" "$matched_segments" "$calo_only_segments" "$zdc_only_segments"

  rm -rf "$tmpdir"
}

pair_zdc_raw_for_auau_lists() {
  case "$LABEL" in
    auau|run2auau) ;;
    *) return 0 ;;
  esac

  local zdc_type="${ZDC_TYPE:-DST_ZDC_RAW}"
  local zdc_tag="${ZDC_TAG:-$TAG}"
  local zdc_dataset="${ZDC_DATASET:-$DATASET}"
  local zdc_rc=0
  local calo_present=0
  local zdc_present=0
  local both_present=0
  local paired=0
  local missing_calo=0
  local missing_zdc=0
  local pair_fail=0
  local calo_total_segments=0
  local zdc_total_segments=0
  local matched_segment_pairs=0
  local calo_only_segments=0
  local zdc_only_segments=0
  local partial_pair_runs=0
  local r r8 calo_file zdc_file primary_stem primary_out
  local report_file="${OUT_DIR}/dst_pairing_report_${LABEL}.txt"
  local unpaired_dir="${OUT_DIR}/unpaired_no_complete_pair"
  local -a calo_present_runs=()
  local -a zdc_present_runs=()
  local -a both_present_runs=()
  local -a paired_runs=()
  local -a partial_pair_runs_list=()
  local -a missing_calo_runs=()
  local -a missing_zdc_runs=()
  local -a missing_both_runs=()
  local -a missing_calo_only_runs=()
  local -a missing_zdc_only_runs=()
  local -a pair_fail_runs=()

  primary_stem="${PREFIX#DST_}"
  primary_stem="$(echo "$primary_stem" | tr '[:upper:]' '[:lower:]')"

  mkdir -p "$unpaired_dir"

  echo "[INFO] Building paired AuAu DATA lists: ${PREFIX} + ${zdc_type}"
  echo "[INFO] Running CreateDstList.pl for dataset=${zdc_dataset} prefix=${zdc_type} tag=${zdc_tag}"
  CreateDstList.pl --tag "$zdc_tag" --dataset "$zdc_dataset" --list "$LIST_FILE" "$zdc_type"
  zdc_rc=$?
  echo "[INFO] CreateDstList.pl for ${zdc_type} finished with exit code: $zdc_rc"

  if (( zdc_rc != 0 )); then
    {
      echo "DST pairing report for $LABEL"
      echo "Input golden run list: $LIST_FILE"
      echo "Output dir: $OUT_DIR"
      echo "Primary type: $PREFIX"
      echo "Primary dataset: $DATASET"
      echo "Primary tag: $TAG"
      echo "ZDC type: $zdc_type"
      echo "ZDC dataset: $zdc_dataset"
      echo "ZDC tag: $zdc_tag"
      echo
      echo "CreateDstList.pl for ZDC failed with exit code: $zdc_rc"
      echo "Pairing status: ERROR_CREATE_ZDC_LISTS"
    } > "$report_file"
    echo "[WARN] DST pairing report saved to: $report_file"
    echo "[WARN] Pairing status: ERROR_CREATE_ZDC_LISTS"
    return 1
  fi

  for r in "${RUNS_ALL[@]}"; do
    r8=$(printf "%08d" "$((10#$r))")
    calo_file="$(find_run_list_for_type "$PREFIX" "$r8" || true)"
    zdc_file="$(find_run_list_for_type "$zdc_type" "$r8" || true)"
    primary_out="${OUT_DIR}/dst_${primary_stem}-${r8}.list"

    if [[ -n "$calo_file" ]]; then
      ((calo_present+=1))
      calo_present_runs+=( "$r8" )
      calo_total_segments=$(( calo_total_segments + $(awk 'NF{n++} END{print n+0}' "$calo_file") ))
    fi

    if [[ -n "$zdc_file" ]]; then
      ((zdc_present+=1))
      zdc_present_runs+=( "$r8" )
      zdc_total_segments=$(( zdc_total_segments + $(awk 'NF{n++} END{print n+0}' "$zdc_file") ))
    fi

    if [[ -z "$calo_file" && -z "$zdc_file" ]]; then
      ((missing_calo+=1))
      ((missing_zdc+=1))
      missing_calo_runs+=( "$r8" )
      missing_zdc_runs+=( "$r8" )
      missing_both_runs+=( "$r8" )
      continue
    fi

    if [[ -z "$calo_file" ]]; then
      ((missing_calo+=1))
      missing_calo_runs+=( "$r8" )
      missing_calo_only_runs+=( "$r8" )
      continue
    fi

    if [[ -z "$zdc_file" ]]; then
      ((missing_zdc+=1))
      missing_zdc_runs+=( "$r8" )
      missing_zdc_only_runs+=( "$r8" )

      mv -f "$calo_file" "${unpaired_dir}/$(basename "$calo_file")" 2>/dev/null || true
      if [[ "$calo_file" != "$primary_out" ]]; then
        rm -f "$primary_out"
      fi

      continue
    fi

    ((both_present+=1))
    both_present_runs+=( "$r8" )

    if pair_calo_zdc_for_run "$r8" "$calo_file" "$zdc_file" "$primary_out"; then
      ((paired+=1))
      paired_runs+=( "$r8" )

      matched_segment_pairs=$(( matched_segment_pairs + PAIR_MATCHED_SEG ))
      calo_only_segments=$(( calo_only_segments + PAIR_CALO_ONLY_SEG ))
      zdc_only_segments=$(( zdc_only_segments + PAIR_ZDC_ONLY_SEG ))

      if (( PAIR_CALO_ONLY_SEG > 0 || PAIR_ZDC_ONLY_SEG > 0 )); then
        ((partial_pair_runs+=1))
        partial_pair_runs_list+=( "${r8} matched=${PAIR_MATCHED_SEG} caloOnly=${PAIR_CALO_ONLY_SEG} zdcOnly=${PAIR_ZDC_ONLY_SEG} mode=${PAIR_MODE}" )
        echo "[WARN] Partial segment overlap for run ${r8}: matched=${PAIR_MATCHED_SEG} caloOnly=${PAIR_CALO_ONLY_SEG} zdcOnly=${PAIR_ZDC_ONLY_SEG}"
      fi

      if [[ "$calo_file" != "$primary_out" ]]; then
        rm -f "$calo_file"
      fi
    else
      ((pair_fail+=1))
      pair_fail_runs+=( "$r8" )

      mv -f "$calo_file" "${unpaired_dir}/$(basename "$calo_file")" 2>/dev/null || true
      if [[ "$calo_file" != "$primary_out" ]]; then
        rm -f "$primary_out"
      fi
    fi
  done

  local overlap_pair_status="NO_OVERLAP"
  if (( both_present > 0 && pair_fail == 0 && paired == both_present )); then
    overlap_pair_status="YES"
  elif (( both_present > 0 && pair_fail > 0 )); then
    overlap_pair_status="NO"
  fi

  local missing_sets_identical="NO"
  if (( missing_calo == missing_zdc && ${#missing_calo_only_runs[@]} == 0 && ${#missing_zdc_only_runs[@]} == 0 )); then
    missing_sets_identical="YES"
  fi

  local coverage_status="COMPLETE"
  if (( missing_calo > 0 || missing_zdc > 0 )); then
    coverage_status="PARTIAL"
  fi

  local pairing_status="OK"
  if (( paired == 0 )); then
    pairing_status="ERROR_NO_PAIRED_RUNS"
  elif (( pair_fail > 0 )); then
    pairing_status="ERROR_SEGMENT_PAIRING"
  elif (( partial_pair_runs > 0 )); then
    pairing_status="OK_PARTIAL_SEGMENTS_EXCLUDED"
  elif (( missing_calo > 0 || missing_zdc > 0 )); then
    pairing_status="OK_PARTIAL_RUN_COVERAGE"
  fi

  local original_run_count="${#RUNS_ALL[@]}"
  local run_loss=$(( original_run_count - paired ))
  local run_retained_pct
  local run_lost_pct
  local segment_retained_pct
  local segment_lost_pct
  run_retained_pct="$(awk -v n="$paired" -v d="$original_run_count" 'BEGIN{ if (d>0) printf "%.2f", 100*n/d; else printf "0.00"; }')"
  run_lost_pct="$(awk -v n="$run_loss" -v d="$original_run_count" 'BEGIN{ if (d>0) printf "%.2f", 100*n/d; else printf "0.00"; }')"
  segment_retained_pct="$(awk -v n="$matched_segment_pairs" -v d="$calo_total_segments" 'BEGIN{ if (d>0) printf "%.2f", 100*n/d; else printf "0.00"; }')"
  segment_lost_pct="$(awk -v n="$calo_only_segments" -v d="$calo_total_segments" 'BEGIN{ if (d>0) printf "%.2f", 100*n/d; else printf "0.00"; }')"

  {
    echo "DST pairing report for $LABEL"
    echo "Input golden run list: $LIST_FILE"
    echo "Output dir: $OUT_DIR"
    echo "Primary type: $PREFIX"
    echo "Primary dataset: $DATASET"
    echo "Primary tag: $TAG"
    echo "ZDC type: $zdc_type"
    echo "ZDC dataset: $zdc_dataset"
    echo "ZDC tag: $zdc_tag"
    echo
    echo "Run-level summary:"
    echo "Total requested runs: ${#RUNS_ALL[@]}"
    echo "Runs with ${PREFIX} list: $calo_present"
    echo "Runs with ${zdc_type} list: $zdc_present"
    echo "Runs with both lists available: $both_present"
    echo "Runs with at least one matched CALOFITTING/ZDC segment pair: $paired"
    echo "Runs with partial segment overlap kept: $partial_pair_runs"
    echo "Runs with zero segment-pairing overlap/failure: $pair_fail"
    echo "Coverage status: $coverage_status"
    echo "Pairing status: $pairing_status"
    echo
    echo "Segment-level summary:"
    echo "${PREFIX} segments available in covered runs: $calo_total_segments"
    echo "${zdc_type} segments available in covered runs: $zdc_total_segments"
    echo "Matched CALOFITTING/ZDC segment pairs written: $matched_segment_pairs"
    echo "CALOFITTING-only segments excluded: $calo_only_segments"
    echo "ZDC-only segments excluded: $zdc_only_segments"
    echo "Segment retained fraction vs available ${PREFIX}: ${segment_retained_pct}%"
    echo "Segment loss fraction vs available ${PREFIX}: ${segment_lost_pct}%"
    echo
    echo "Missing-run summary:"
    echo "Missing ${PREFIX} runs: $missing_calo"
    echo "Missing ${zdc_type} runs: $missing_zdc"
    echo "Missing both ${PREFIX} and ${zdc_type}: ${#missing_both_runs[@]}"
    echo "Missing ${PREFIX} only: ${#missing_calo_only_runs[@]}"
    echo "Missing ${zdc_type} only: ${#missing_zdc_only_runs[@]}"
    echo "Missing-run sets identical between ${PREFIX} and ${zdc_type}: $missing_sets_identical"
    echo
    echo "Retention relative to input GRL:"
    echo "Input GRL runs: $original_run_count"
    echo "Final analyzable paired runs: $paired"
    echo "Dropped runs: $run_loss"
    echo "Run retained fraction: ${run_retained_pct}%"
    echo "Run loss fraction: ${run_lost_pct}%"
    echo "Input available ${PREFIX} segments: $calo_total_segments"
    echo "Final matched segment pairs: $matched_segment_pairs"
    echo "Excluded ${PREFIX} segments: $calo_only_segments"
    echo "Segment retained fraction: ${segment_retained_pct}%"
    echo "Segment loss fraction: ${segment_lost_pct}%"
    echo
    echo "Run-list policy:"
    echo "Using the configured input GRL directly. No old/new AuAu production overlap"
    echo "comparison is produced by this script."
    echo "Configured input GRL: $LIST_FILE"
    echo
    echo "Runs with ${PREFIX} list:"
    if (( ${#calo_present_runs[@]} > 0 )); then
      printf '%s\n' "${calo_present_runs[@]}"
    else
      echo "(none)"
    fi
    echo
    echo "Runs with ${zdc_type} list:"
    if (( ${#zdc_present_runs[@]} > 0 )); then
      printf '%s\n' "${zdc_present_runs[@]}"
    else
      echo "(none)"
    fi
    echo
    echo "Runs with both lists available:"
    if (( ${#both_present_runs[@]} > 0 )); then
      printf '%s\n' "${both_present_runs[@]}"
    else
      echo "(none)"
    fi
    echo
    echo "Runs with at least one matched CALOFITTING/ZDC segment pair:"
    if (( ${#paired_runs[@]} > 0 )); then
      printf '%s\n' "${paired_runs[@]}"
    else
      echo "(none)"
    fi
    echo
    echo "Runs with partial segment overlap kept:"
    if (( ${#partial_pair_runs_list[@]} > 0 )); then
      printf '%s\n' "${partial_pair_runs_list[@]}"
    else
      echo "(none)"
    fi
    echo
    echo "Missing ${PREFIX} runs:"
    if (( ${#missing_calo_runs[@]} > 0 )); then
      printf '%s\n' "${missing_calo_runs[@]}"
    else
      echo "(none)"
    fi
    echo
    echo "Missing ${zdc_type} runs:"
    if (( ${#missing_zdc_runs[@]} > 0 )); then
      printf '%s\n' "${missing_zdc_runs[@]}"
    else
      echo "(none)"
    fi
    echo
    echo "Missing both:"
    if (( ${#missing_both_runs[@]} > 0 )); then
      printf '%s\n' "${missing_both_runs[@]}"
    else
      echo "(none)"
    fi
    echo
    echo "Missing ${PREFIX} only:"
    if (( ${#missing_calo_only_runs[@]} > 0 )); then
      printf '%s\n' "${missing_calo_only_runs[@]}"
    else
      echo "(none)"
    fi
    echo
    echo "Missing ${zdc_type} only:"
    if (( ${#missing_zdc_only_runs[@]} > 0 )); then
      printf '%s\n' "${missing_zdc_only_runs[@]}"
    else
      echo "(none)"
    fi
    echo
    echo "Zero-overlap / segment-pairing failure runs:"
    if (( ${#pair_fail_runs[@]} > 0 )); then
      printf '%s\n' "${pair_fail_runs[@]}"
    else
      echo "(none)"
    fi
    echo
    echo "Unpaired one-column ${PREFIX} lists moved to:"
    echo "$unpaired_dir"
  } > "$report_file"

  echo "----------------------------------------"
  echo "DST pairing diagnostic for $LABEL"
  echo "  total requested runs                       : ${#RUNS_ALL[@]}"
  echo "  runs with ${PREFIX} list                    : $calo_present"
  echo "  runs with ${zdc_type} list                   : $zdc_present"
  echo "  runs with both lists available              : $both_present"
  echo "  runs with matched segment pairs             : $paired"
  echo "  partial-overlap runs kept                   : $partial_pair_runs"
  echo "  zero-overlap / segment-pairing failures     : $pair_fail"
  echo "  ${PREFIX} segments available                : $calo_total_segments"
  echo "  ${zdc_type} segments available               : $zdc_total_segments"
  echo "  matched CALOFITTING/ZDC segment pairs       : $matched_segment_pairs"
  echo "  CALOFITTING-only segments excluded          : $calo_only_segments"
  echo "  ZDC-only segments excluded                  : $zdc_only_segments"
  echo "  run retained fraction                       : ${run_retained_pct}%"
  echo "  segment retained fraction                   : ${segment_retained_pct}%"
  echo "  missing ${PREFIX} runs                      : $missing_calo"
  echo "  missing ${zdc_type} runs                     : $missing_zdc"
  echo "  missing both                               : ${#missing_both_runs[@]}"
  echo "  missing ${PREFIX} only                      : ${#missing_calo_only_runs[@]}"
  echo "  missing ${zdc_type} only                     : ${#missing_zdc_only_runs[@]}"
  echo "  missing-run sets identical                  : $missing_sets_identical"
  echo "  coverage status                            : $coverage_status"
  echo "  pairing status                             : $pairing_status"
  echo "  report                                     : $report_file"

  if (( paired == 0 || pair_fail > 0 )); then
    return 1
  fi

  return 0
}

collect_mode_build_metrics() {
  SUMMARY_REQUESTED_TYPE="$(mode_requested_type)"

  if [[ -f "$LIST_FILE" ]]; then
    read_runs_from_file "$LIST_FILE" _summary_runs
    SUMMARY_RUNS=${#_summary_runs[@]}
  else
    SUMMARY_RUNS=0
  fi

  mapfile -t _summary_files < <(mode_output_files)
  SUMMARY_LISTS=${#_summary_files[@]}
  SUMMARY_ENTRIES=0

  for _summary_file in "${_summary_files[@]}"; do
    SUMMARY_ENTRIES=$(( SUMMARY_ENTRIES + $(awk 'NF{n++} END{print n+0}' "$_summary_file") ))
  done

  SUMMARY_OTHER="$(probe_other_types_summary)"
}

print_all_build_summary_block() {
  local status="$1"

  echo "----------------------------------------"
  echo "ALL summary for $LABEL"
  echo "Requested type: $SUMMARY_REQUESTED_TYPE"
  echo "Requested runs: $SUMMARY_RUNS"
  echo "List files present: $SUMMARY_LISTS"
  echo "Total DST entries: $SUMMARY_ENTRIES"
  echo "Other available types: $SUMMARY_OTHER"
  echo "Status: $status"
}

print_all_final_build_summary() {
  local -n rows_ref="$1"
  local total_runs=0
  local total_lists=0
  local total_entries=0
  local errors=0

  echo
  echo "=================================================="
  echo "FINAL ALL BUILD SUMMARY"
  echo "=================================================="

  for row in "${rows_ref[@]}"; do
    IFS='|' read -r label requested runs lists entries other status <<< "$row"
    echo "[$label] $status"
    echo "  requested type: $requested"
    echo "  requested runs: $runs"
    echo "  list files present: $lists"
    echo "  total DST entries: $entries"
    echo "  other available types: $other"
    total_runs=$(( total_runs + runs ))
    total_lists=$(( total_lists + lists ))
    total_entries=$(( total_entries + entries ))
    [[ "$status" == OK ]] || ((errors+=1))
  done

  echo "--------------------------------------------------"
  echo "Total requested runs across ALL: $total_runs"
  echo "Total list files present across ALL: $total_lists"
  echo "Total DST entries across ALL: $total_entries"
  echo "Modes with errors: $errors / ${#rows_ref[@]}"
}

print_all_final_qa_summary() {
  local -n rows_ref="$1"
  local total_runs=0
  local errors=0

  echo
  echo "=================================================="
  echo "FINAL ALL QA SUMMARY"
  echo "=================================================="

  for row in "${rows_ref[@]}"; do
    IFS='|' read -r label runs status <<< "$row"
    echo "[$label] $status"
    echo "  input runs: $runs"
    total_runs=$(( total_runs + runs ))
    [[ "$status" == OK ]] || ((errors+=1))
  done

  echo "--------------------------------------------------"
  echo "Total input runs across ALL QA scans: $total_runs"
  echo "Modes with errors: $errors / ${#rows_ref[@]}"
}

run_all_modes() {
  local -a modes=(auau pp24 pp25 oo25 run2auau)
  local -a summary_rows=()
  local mode rc status

  echo
  echo "=================================================="
  echo "BEGIN ALL${ACTION:+ $ACTION}"
  echo "=================================================="
  echo "[INFO] Modes to process: ${modes[*]}"

  for mode in "${modes[@]}"; do
    MODE="$mode"
    echo
    echo "[INFO] Setting up mode: $MODE"
    setup_mode

    echo
    echo "=================================================="
    echo "BEGIN $LABEL${ACTION:+ $ACTION}"
    echo "=================================================="
    echo "[INFO] Input list: $LIST_FILE"
    echo "[INFO] Output dir: $OUT_DIR"
    echo "[INFO] Requested type: $(mode_requested_type)"
    echo "[INFO] Build mode: $BUILD_MODE"

    if [[ "$ACTION" == "QA" ]]; then
      echo "[INFO] QA-only mode for $LABEL: using only the input golden run list."
      echo "[INFO] QA-only mode for $LABEL: no DST lists will be created, removed, modified, or scanned."
      echo "[INFO] Starting trigger QA for $LABEL"
      if run_trigger_qa; then
        rc=0
      else
        rc=$?
      fi

      if [[ -f "$LIST_FILE" ]]; then
        read_runs_from_file "$LIST_FILE" _qa_runs
        qa_runs=${#_qa_runs[@]}
      else
        qa_runs=0
      fi

      if (( rc == 0 )); then
        status="OK"
      else
        status="ERROR($rc)"
      fi

      echo "[INFO] Finished $LABEL${ACTION:+ $ACTION} with status: $status"
      summary_rows+=( "$LABEL|$qa_runs|$status" )
      continue
    fi

    echo "[INFO] Starting build for $LABEL"
    if [[ "$BUILD_MODE" == "per_run" ]]; then
      if build_pp24_lists; then
        rc=0
      else
        rc=$?
      fi
    else
      if build_dataset_lists; then
        rc=0
      else
        rc=$?
      fi
    fi

    echo "[INFO] Collecting post-build summary for $LABEL"
    collect_mode_build_metrics

    if (( rc == 0 )); then
      status="OK"
    else
      status="ERROR($rc)"
    fi

    print_all_build_summary_block "$status"
    echo "[INFO] Finished $LABEL${ACTION:+ $ACTION} with status: $status"
    summary_rows+=( "$LABEL|$SUMMARY_REQUESTED_TYPE|$SUMMARY_RUNS|$SUMMARY_LISTS|$SUMMARY_ENTRIES|$SUMMARY_OTHER|$status" )
  done

  if [[ "$ACTION" == "QA" ]]; then
    print_all_final_qa_summary summary_rows
  else
    print_all_final_build_summary summary_rows
  fi

  for row in "${summary_rows[@]}"; do
    [[ "$row" == *"|ERROR("* ]] && return 1
  done

  return 0
}

setup_mode() {
  case "$MODE" in
    pp24)
      LABEL="pp24"
      LIST_FILE="$GRL_BASE/run2pp_ana509_2024p022_v001_dst_calofitting_grl.list"
      OUT_DIR="$IN_BASE/dst_lists_pp"
      TAG="ana521_2025p007_v001"
      TYPE="DST_JETCALO"
      DATASET="run2pp"
      PREFIX="DST_JETCALO"
      BUILD_MODE="per_run"
      ;;
    pp25)
      LABEL="pp25"
      LIST_FILE="$GRL_BASE/run3pp_new_newcdbtag_v008_dst_calofitting_grl.list"
      OUT_DIR="$IN_BASE/dst_lists_pp_run25"
      TAG="new_newcdbtag_v008"
      DATASET="run3pp"
      PREFIX="DST_CALOFITTING"
      TYPE="DST_CALOFITTING"
      BUILD_MODE="dataset_list"
      ;;
    auau)
      LABEL="auau"
      LIST_FILE="${AUAU_GRL_FILE:-$GRL_BASE/run3auau_pro001_pcdb001_v001_dst_calofitting_grl.list}"
      OUT_DIR="$IN_BASE/dst_lists_auau"
      TAG="${AUAU_TAG:-pro001_pcdb001_v001}"
      DATASET="run3auau"
      PREFIX="DST_CALOFITTING"
      TYPE="DST_CALOFITTING"
      BUILD_MODE="dataset_list"
      ;;
    oo25)
      LABEL="oo25"
      LIST_FILE="$GRL_BASE/run3oo_ana536_2025p010_v001_dst_calofitting_grl.list"
      OUT_DIR="$IN_BASE/dst_lists_oo"
      TAG="ana536_2025p010_v001"
      DATASET="run3oo"
      PREFIX="DST_CALOFITTING"
      TYPE="DST_CALOFITTING"
      BUILD_MODE="dataset_list"
      ;;
    run2auau)
      LABEL="run2auau"
      LIST_FILE="$GRL_BASE/run2auau_ana509_2024p022_v001_dst_calofitting_grl.list"
      OUT_DIR="$IN_BASE/dst_lists_auau_run2"
      TAG="ana509_2024p022_v001"
      DATASET="run2auau"
      PREFIX="DST_CALOFITTING"
      TYPE="DST_CALOFITTING"
      BUILD_MODE="dataset_list"
      ;;
    *)
      usage
      exit 1
      ;;
  esac
}

run_trigger_qa() {
  command -v psql >/dev/null 2>&1 || { echo "[ERROR] psql not found in PATH"; exit 1; }
  [[ -f "$LIST_FILE" ]] || { echo "[ERROR] Run list not found: $LIST_FILE"; exit 1; }
  [[ -s "$LIST_FILE" ]] || { echo "[ERROR] Run list is empty: $LIST_FILE"; exit 1; }

  read_runs_from_file "$LIST_FILE" RUNS_ALL
  (( ${#RUNS_ALL[@]} > 0 )) || { echo "[ERROR] No valid run numbers found in: $LIST_FILE"; exit 1; }

  PSQL=(psql -h sphnxdaqdbreplica -d daq -At -F $'\t' -q)
  sql()         { "${PSQL[@]}" -c "$1" 2>/dev/null || true; }
  num_or_zero() { [[ $1 =~ ^-?[0-9]+([.][0-9]+)?$ ]] && printf '%s' "$1" || printf 0; }

  local ANSI_BOLD ANSI_DIM ANSI_CYAN ANSI_GREEN ANSI_YELLOW ANSI_RED ANSI_RESET
  ANSI_BOLD=$'\033[1m'
  ANSI_DIM=$'\033[2m'
  ANSI_CYAN=$'\033[36m'
  ANSI_GREEN=$'\033[32m'
  ANSI_YELLOW=$'\033[33m'
  ANSI_RED=$'\033[31m'
  ANSI_RESET=$'\033[0m'

  echo
  printf "%sTrigger QA summary for: %s%s\n" "$ANSI_BOLD$ANSI_CYAN" "$LABEL" "$ANSI_RESET"
  echo "Input golden run list: $LIST_FILE"
  echo "Runs in input list: ${#RUNS_ALL[@]}"
  echo "Requested QA metrics per trigger:"
  echo "  - trigger-menu bit -> name stability across all runs"
  echo "  - first adjacent run where the trigger bit -> name map changes"
  echo "  - menu presence by run"
  echo "  - raw-positive run coverage and total raw counts"
  echo "  - live-positive run coverage and total live counts"
  echo "  - scaled-positive run coverage and total scaled counts (reference only)"
  echo "  - scaledown averages (reference only)"
  echo "QA mode behavior:"
  echo "  - reads only the input golden run list"
  echo "  - queries trigger metadata from daq tables"
  echo "  - does not create, clean, scan, or modify any DST list files"
  echo "doNotScale interpretation guidance:"
  echo "  - raw/live are the primary quantities to inspect for raw-trigger efficiency studies"
  echo "  - scaled/scaledown are kept only as secondary reference QA"
  echo

  echo "[INFO] Verifying psql connectivity expectation: host=sphnxdaqdbreplica db=daq"
  echo "[INFO] Beginning trigger-menu scan over ${#RUNS_ALL[@]} golden-list runs"

  declare -A TRIG_NAME
  declare -A TRIG_MENU_RUNS
  declare -A TRIG_MENU_FIRSTRUN
  declare -A TRIG_RAWON_RUNS
  declare -A TRIG_RAW_FIRSTRUN
  declare -A TRIG_SUMRAW
  declare -A TRIG_LIVEON_RUNS
  declare -A TRIG_LIVE_FIRSTRUN
  declare -A TRIG_SUMLIVE
  declare -A TRIG_SCALEDON_RUNS
  declare -A TRIG_SCALED_FIRSTRUN
  declare -A TRIG_SUMSCALED
  declare -A TRIG_SUMSDFACTOR
  declare -A TRIG_SDCOUNT

  local total_menu_rows=0
  local runs_with_no_menu_rows=0
  local total_scaler_rows=0
  local runs_with_no_scaler_rows=0
  local rows_bad_idx=0
  local rows_raw_positive=0
  local rows_live_positive=0
  local rows_scaled_minus_one=0
  local rows_scaled_zero=0
  local rows_scaled_positive=0
  local rows_used_for_avg=0
  local first_problem_run=""

  local first_menu_run=""
  local prev_menu_run=""
  local first_menu_change_run=""
  local first_menu_change_prev_run=""
  local menu_change_count=0
  local -a first_menu_change_rows=()

  declare -A PREV_MENU=()

  midx=0
  for run in "${RUNS_ALL[@]}"; do
    ((midx+=1))

    unset CUR_MENU
    declare -A CUR_MENU=()

    local menu_row_count=0

    while IFS=$'\t' read -r idx trg; do
      [[ -z "$idx" ]] && continue

      idx=$(num_or_zero "$idx")

      if [[ ! "$idx" =~ ^[0-9]+$ ]]; then
        ((rows_bad_idx+=1))
        [[ -z "$first_problem_run" ]] && first_problem_run="$run"
        continue
      fi

      [[ -z "$trg" ]] && trg="(missing-name)"

      CUR_MENU["$idx"]="$trg"
      [[ -z "${TRIG_NAME["$idx"]:-}" ]] && TRIG_NAME["$idx"]="$trg"

      TRIG_MENU_RUNS["$idx"]=$(( ${TRIG_MENU_RUNS["$idx"]:-0} + 1 ))
      [[ -z "${TRIG_MENU_FIRSTRUN["$idx"]:-}" ]] && TRIG_MENU_FIRSTRUN["$idx"]="$run"

      ((menu_row_count+=1))
      ((total_menu_rows+=1))
    done < <(sql "SELECT index, triggername
                  FROM gl1_triggernames
                  WHERE $run BETWEEN runnumber AND runnumber_last
                  ORDER BY index;")

    if (( menu_row_count == 0 )); then
      ((runs_with_no_menu_rows+=1))
      [[ -z "$first_problem_run" ]] && first_problem_run="$run"
      printf "%s[WARN]%s Run %s returned zero rows from gl1_triggernames menu query\n" "$ANSI_YELLOW" "$ANSI_RESET" "$run"
    fi

    if [[ -z "$first_menu_run" ]]; then
      first_menu_run="$run"
    else
      unset MENU_KEYS
      declare -A MENU_KEYS=()

      for idx in "${!PREV_MENU[@]}"; do
        MENU_KEYS["$idx"]=1
      done

      for idx in "${!CUR_MENU[@]}"; do
        MENU_KEYS["$idx"]=1
      done

      local menu_changed=0
      local -a diff_rows=()

      for idx in "${!MENU_KEYS[@]}"; do
        prev_name="${PREV_MENU["$idx"]:-"(absent)"}"
        cur_name="${CUR_MENU["$idx"]:-"(absent)"}"

        if [[ "$prev_name" != "$cur_name" ]]; then
          menu_changed=1
          diff_rows+=( "$idx|$prev_name|$cur_name" )
        fi
      done

      if (( menu_changed )); then
        ((menu_change_count+=1))

        if [[ -z "$first_menu_change_run" ]]; then
          first_menu_change_run="$run"
          first_menu_change_prev_run="$prev_menu_run"
          first_menu_change_rows=( "${diff_rows[@]}" )
        fi
      fi
    fi

    unset PREV_MENU
    declare -A PREV_MENU=()

    for idx in "${!CUR_MENU[@]}"; do
      PREV_MENU["$idx"]="${CUR_MENU["$idx"]}"
    done

    prev_menu_run="$run"

    if (( midx <= 5 || midx % 200 == 0 || midx == ${#RUNS_ALL[@]} )); then
      printf "[INFO] Menu QA progress: %d / %d | run=%s | menuBits=%d\n" \
        "$midx" "${#RUNS_ALL[@]}" "$run" "$menu_row_count"
    fi
  done

  echo "[INFO] Beginning raw/live/scaled scaler scan over ${#RUNS_ALL[@]} golden-list runs"

  declare -A BIT_ACTIVE=()

  gidx=0
  for run in "${RUNS_ALL[@]}"; do
    ((gidx+=1))

    unset SEEN_RAW
    unset SEEN_LIVE
    unset SEEN_SCALED
    declare -A SEEN_RAW=()
    declare -A SEEN_LIVE=()
    declare -A SEEN_SCALED=()

    local run_scaler_row_count=0
    local run_raw_positive_count=0
    local run_live_positive_count=0
    local run_scaled_positive_count=0

    while IFS=$'\t' read -r trg idx scaled raw live; do
      [[ -z "$idx" ]] && continue

      ((run_scaler_row_count+=1))
      ((total_scaler_rows+=1))

      idx=$(num_or_zero "$idx")
      scaled=$(num_or_zero "$scaled")
      raw=$(num_or_zero "$raw")
      live=$(num_or_zero "$live")

      if [[ ! "$idx" =~ ^[0-9]+$ ]]; then
        ((rows_bad_idx+=1))
        [[ -z "$first_problem_run" ]] && first_problem_run="$run"
        continue
      fi

      [[ -z "$trg" ]] && trg="${TRIG_NAME["$idx"]:-"(missing-name)"}"
      [[ -z "${TRIG_NAME["$idx"]:-}" ]] && TRIG_NAME["$idx"]="$trg"

      if (( scaled == -1 )); then
        ((rows_scaled_minus_one+=1))
      elif (( scaled == 0 )); then
        ((rows_scaled_zero+=1))
      elif (( scaled > 0 )); then
        ((rows_scaled_positive+=1))
      fi

      if (( scaled != -1 )); then
        BIT_ACTIVE["${run}_${idx}"]=1
      fi

      if (( raw > 0 )); then
        TRIG_SUMRAW["$idx"]=$(( ${TRIG_SUMRAW["$idx"]:-0} + raw ))
        ((run_raw_positive_count+=1))
        ((rows_raw_positive+=1))

        if [[ -z "${SEEN_RAW["$idx"]:-}" ]]; then
          TRIG_RAWON_RUNS["$idx"]=$(( ${TRIG_RAWON_RUNS["$idx"]:-0} + 1 ))
          SEEN_RAW["$idx"]=1
        fi

        [[ -z "${TRIG_RAW_FIRSTRUN["$idx"]:-}" ]] && TRIG_RAW_FIRSTRUN["$idx"]="$run"
      fi

      if (( live > 0 )); then
        TRIG_SUMLIVE["$idx"]=$(( ${TRIG_SUMLIVE["$idx"]:-0} + live ))
        ((run_live_positive_count+=1))
        ((rows_live_positive+=1))

        if [[ -z "${SEEN_LIVE["$idx"]:-}" ]]; then
          TRIG_LIVEON_RUNS["$idx"]=$(( ${TRIG_LIVEON_RUNS["$idx"]:-0} + 1 ))
          SEEN_LIVE["$idx"]=1
        fi

        [[ -z "${TRIG_LIVE_FIRSTRUN["$idx"]:-}" ]] && TRIG_LIVE_FIRSTRUN["$idx"]="$run"
      fi

      if (( scaled > 0 )); then
        TRIG_SUMSCALED["$idx"]=$(( ${TRIG_SUMSCALED["$idx"]:-0} + scaled ))
        ((run_scaled_positive_count+=1))

        if [[ -z "${SEEN_SCALED["$idx"]:-}" ]]; then
          TRIG_SCALEDON_RUNS["$idx"]=$(( ${TRIG_SCALEDON_RUNS["$idx"]:-0} + 1 ))
          SEEN_SCALED["$idx"]=1
        fi

        [[ -z "${TRIG_SCALED_FIRSTRUN["$idx"]:-}" ]] && TRIG_SCALED_FIRSTRUN["$idx"]="$run"
      fi

      if (( scaled > 0 && live >= 0 )); then
        sdf=$(awk -v l="$live" -v s="$scaled" 'BEGIN{ if (s>0) printf "%.6f", l/s; else printf "0"; }')
        TRIG_SUMSDFACTOR["$idx"]="$(awk -v a="${TRIG_SUMSDFACTOR["$idx"]:-0}" -v b="$sdf" 'BEGIN{ printf "%.6f", a+b }')"
        TRIG_SDCOUNT["$idx"]=$(( ${TRIG_SDCOUNT["$idx"]:-0} + 1 ))
        ((rows_used_for_avg+=1))
      fi
    done < <(sql "SELECT COALESCE(t.triggername,''), s.index, s.scaled, s.raw, s.live
                   FROM gl1_scalers s
                   LEFT JOIN gl1_triggernames t
                     ON s.index = t.index
                    AND s.runnumber BETWEEN t.runnumber AND t.runnumber_last
                  WHERE s.runnumber = $run
                  ORDER BY s.index;")

    if (( run_scaler_row_count == 0 )); then
      ((runs_with_no_scaler_rows+=1))
      [[ -z "$first_problem_run" ]] && first_problem_run="$run"
      printf "%s[WARN]%s Run %s returned zero rows from gl1_scalers query\n" "$ANSI_YELLOW" "$ANSI_RESET" "$run"
    fi

    if (( gidx <= 5 || gidx % 200 == 0 || gidx == ${#RUNS_ALL[@]} )); then
      printf "[INFO] Scaler QA progress: %d / %d | run=%s | rows=%d | raw>0 rows=%d | live>0 rows=%d | scaled>0 rows=%d\n" \
        "$gidx" "${#RUNS_ALL[@]}" "$run" "$run_scaler_row_count" "$run_raw_positive_count" "$run_live_positive_count" "$run_scaled_positive_count"
    fi
  done

  tmp_rows=()
  maxlen=11
  for idx in "${!TRIG_NAME[@]}"; do
    trg="${TRIG_NAME["$idx"]}"
    [[ -z "$trg" ]] && trg="(missing-name)"
    (( ${#trg} > maxlen )) && maxlen=${#trg}

    menu_runs=${TRIG_MENU_RUNS["$idx"]:-0}
    raw_runs=${TRIG_RAWON_RUNS["$idx"]:-0}
    live_runs=${TRIG_LIVEON_RUNS["$idx"]:-0}
    scaled_runs=${TRIG_SCALEDON_RUNS["$idx"]:-0}
    sumraw=${TRIG_SUMRAW["$idx"]:-0}
    sumlive=${TRIG_SUMLIVE["$idx"]:-0}
    sumscaled=${TRIG_SUMSCALED["$idx"]:-0}

    if (( sumraw > 0 )); then
      live_over_raw=$(awk -v l="$sumlive" -v r="$sumraw" 'BEGIN{ if (r>0) printf "%.3f", l/r; else printf "0.000"; }')
    else
      live_over_raw="0.000"
    fi

    if (( ${TRIG_SDCOUNT["$idx"]:-0} > 0 )); then
      avg_sd=$(awk -v a="${TRIG_SUMSDFACTOR["$idx"]:-0}" -v n="${TRIG_SDCOUNT["$idx"]:-0}" 'BEGIN{ if (n>0) printf "%.3f", a/n; else printf "0.000"; }')
    else
      avg_sd="0.000"
    fi

    tmp_rows+=( "$idx|$trg|$menu_runs|$raw_runs|$live_runs|$scaled_runs|$sumraw|$sumlive|$sumscaled|$live_over_raw|$avg_sd" )
  done
  (( maxlen+=2 ))

  echo
  printf "%sTrigger-menu QA%s\n" "$ANSI_BOLD$ANSI_CYAN" "$ANSI_RESET"
  printf "  Total menu rows read: %d\n" "$total_menu_rows"
  printf "  Runs with zero menu rows: %d\n" "$runs_with_no_menu_rows"
  printf "  Distinct trigger-menu epochs across adjacent GRL runs: %d\n" "$(( menu_change_count + 1 ))"
  printf "  Adjacent-run trigger-menu changes observed: %d\n" "$menu_change_count"
  printf "  First GRL run used as menu baseline: %s\n" "$first_menu_run"

  if [[ -n "$first_menu_change_run" ]]; then
    echo "  First adjacent-run menu change detected:"
    printf "    previous run: %s\n" "$first_menu_change_prev_run"
    printf "    current  run: %s\n" "$first_menu_change_run"
    echo "    Changed bits:"
    printf "%s\n" "${first_menu_change_rows[@]}" \
    | awk -F'|' '{printf "%04d|%s|%s\n",$1,$2,$3}' \
    | sort -n \
    | while IFS='|' read -r padIdx oldname newname; do
        idx=$((10#$padIdx))
        printf "      bit %d : %s  -->  %s\n" "$idx" "$oldname" "$newname"
      done
  else
    printf "  No adjacent-run trigger-menu changes were detected across the provided GRL.\n"
  fi

  echo
  printf "%sRaw/Live/Scaled bookkeeping diagnostics%s\n" "$ANSI_BOLD$ANSI_CYAN" "$ANSI_RESET"
  printf "  Total scaler rows read: %d\n" "$total_scaler_rows"
  printf "  Runs with zero scaler rows: %d\n" "$runs_with_no_scaler_rows"
  printf "  Row counts by content class:\n"
  printf "    raw   > 0 : %d\n" "$rows_raw_positive"
  printf "    live  > 0 : %d\n" "$rows_live_positive"
  printf "    scaled = -1 : %d\n" "$rows_scaled_minus_one"
  printf "    scaled = 0  : %d\n" "$rows_scaled_zero"
  printf "    scaled > 0  : %d\n" "$rows_scaled_positive"
  printf "  Rows used for average prescale calculation, live/scaled (reference only): %d\n" "$rows_used_for_avg"
  printf "  Rows with unusable trigger bit values: %d\n" "$rows_bad_idx"
  if [[ -n "$first_problem_run" ]]; then
    printf "%s  First run showing a warning condition: %s%s\n" "$ANSI_YELLOW" "$first_problem_run" "$ANSI_RESET"
  else
    printf "%s  No per-run structural warning conditions were encountered during scan%s\n" "$ANSI_GREEN" "$ANSI_RESET"
  fi

  if ((${#tmp_rows[@]})); then
    echo
    printf "%sDataset: %s%s\n" "$ANSI_BOLD$ANSI_GREEN" "$LABEL" "$ANSI_RESET"
    printf "%sTotal trigger bits seen in menu/scaler scans: %d%s\n" "$ANSI_DIM" "${#tmp_rows[@]}" "$ANSI_RESET"
    printf "%sTable notes: MenuRuns counts runs where the bit->name exists in the trigger menu; RawRuns and LiveRuns are the primary doNotScale-relevant QA; scaled columns are printed only for reference.%s\n" "$ANSI_DIM" "$ANSI_RESET"
    echo

    printf "%s  %4s | %-*s | %8s | %8s | %8s | %10s | %14s | %14s | %8s | %11s%s\n" \
      "$ANSI_BOLD" "Bit" "$maxlen" "TriggerName" "MenuRuns" "RawRuns" "LiveRuns" "ScaledRuns" "SumRaw" "SumLive" "L/R" "AvgPrescale" "$ANSI_RESET"
    printf "%s  %-4s-+-%-*s-+-%-8s-+-%-8s-+-%-8s-+-%-10s-+-%-14s-+-%-14s-+-%-8s-+-%-11s%s\n" \
      "$ANSI_DIM" \
      "$(printf '─%.0s' $(seq 1 4))" \
      "$maxlen" "$(printf '─%.0s' $(seq 1 $maxlen))" \
      "$(printf '─%.0s' $(seq 1 8))" \
      "$(printf '─%.0s' $(seq 1 8))" \
      "$(printf '─%.0s' $(seq 1 8))" \
      "$(printf '─%.0s' $(seq 1 10))" \
      "$(printf '─%.0s' $(seq 1 14))" \
      "$(printf '─%.0s' $(seq 1 14))" \
      "$(printf '─%.0s' $(seq 1 8))" \
      "$(printf '─%.0s' $(seq 1 11))" \
      "$ANSI_RESET"

    printf "%s\n" "${tmp_rows[@]}" \
    | awk -F'|' '{printf "%016d|%04d|%08d|%08d|%08d|%010d|%016d|%016d|%s|%s|%s\n",$7,$1,$3,$4,$5,$6,$8,$9,$10,$11,$2}' \
    | sort -r -n \
    | while IFS='|' read -r padSumRaw padIdx padMenu padRawRuns padLiveRuns padScaledRuns padSumLive padSumScaled liveOverRaw avgScale trg; do
        sumraw=$((10#$padSumRaw))
        idx=$((10#$padIdx))
        menu_runs=$((10#$padMenu))
        raw_runs=$((10#$padRawRuns))
        live_runs=$((10#$padLiveRuns))
        scaled_runs=$((10#$padScaledRuns))
        sumlive=$((10#$padSumLive))
        printf "  %4d | %-*s | %8d | %8d | %8d | %10d | %14d | %14d | %8s | %11s\n" \
          "$idx" "$maxlen" "$trg" "$menu_runs" "$raw_runs" "$live_runs" "$scaled_runs" "$sumraw" "$sumlive" "$liveOverRaw" "$avgScale"
      done
    echo

    print_focus_table() {
      local title="$1"
      shift
      local -a bits=( "$@" )
      local focus_maxlen=11
      local bit name

      for bit in "${bits[@]}"; do
        name="${TRIG_NAME["$bit"]:-"(missing-name)"}"
        (( ${#name} > focus_maxlen )) && focus_maxlen=${#name}
      done
      (( focus_maxlen+=2 ))

      printf "%s%s%s\n" "$ANSI_BOLD$ANSI_GREEN" "$title" "$ANSI_RESET"
      printf "%s  %4s | %-*s | %8s | %8s | %8s | %8s | %8s | %8s | %12s | %12s%s\n" \
        "$ANSI_BOLD" "Bit" "$focus_maxlen" "TriggerName" "1stMenu" "MenuRuns" "1stRaw" "RawRuns" "1stLive" "LiveRuns" "SumRaw" "SumLive" "$ANSI_RESET"
      printf "%s  %-4s-+-%-*s-+-%-8s-+-%-8s-+-%-8s-+-%-8s-+-%-8s-+-%-8s-+-%-12s-+-%-12s%s\n" \
        "$ANSI_DIM" \
        "$(printf '─%.0s' $(seq 1 4))" \
        "$focus_maxlen" "$(printf '─%.0s' $(seq 1 $focus_maxlen))" \
        "$(printf '─%.0s' $(seq 1 8))" \
        "$(printf '─%.0s' $(seq 1 8))" \
        "$(printf '─%.0s' $(seq 1 8))" \
        "$(printf '─%.0s' $(seq 1 8))" \
        "$(printf '─%.0s' $(seq 1 8))" \
        "$(printf '─%.0s' $(seq 1 8))" \
        "$(printf '─%.0s' $(seq 1 12))" \
        "$(printf '─%.0s' $(seq 1 12))" \
        "$ANSI_RESET"

      for bit in "${bits[@]}"; do
        name="${TRIG_NAME["$bit"]:-"(missing-name)"}"
        first_menu="${TRIG_MENU_FIRSTRUN["$bit"]:--}"
        menu_runs="${TRIG_MENU_RUNS["$bit"]:-0}"
        first_raw="${TRIG_RAW_FIRSTRUN["$bit"]:--}"
        raw_runs="${TRIG_RAWON_RUNS["$bit"]:-0}"
        first_live="${TRIG_LIVE_FIRSTRUN["$bit"]:--}"
        live_runs="${TRIG_LIVEON_RUNS["$bit"]:-0}"
        sumraw="${TRIG_SUMRAW["$bit"]:-0}"
        sumlive="${TRIG_SUMLIVE["$bit"]:-0}"

        printf "  %4d | %-*s | %8s | %8d | %8s | %8d | %8s | %8d | %12d | %12d\n" \
          "$bit" "$focus_maxlen" "$name" "$first_menu" "$menu_runs" "$first_raw" "$raw_runs" "$first_live" "$live_runs" "$sumraw" "$sumlive"
      done
      echo
    }

    print_trigger_groupings() {
      local outfile="$OUT_DIR/trigger_groupings_${LABEL}.txt"
      local epoch_idx=0
      local epoch_start=""
      local prev_run=""
      local run idx trg

      mkdir -p "$OUT_DIR"

      emit_epoch_groupings() {
        local run_start="$1"
        local run_end="$2"
        local map_name="$3"
        local -n epoch_menu_ref="$map_name"
        local mbd_bit mbd_name mbd_thr mbd_vtx pho_bit pho_name pho_thr pho_vtx pho_e
        local -a novtx_rows=()
        local -a vtx_rows=()

        ((epoch_idx+=1))

        for mbd_bit in "${!epoch_menu_ref[@]}"; do
          mbd_name="${epoch_menu_ref["$mbd_bit"]:-}"

          if [[ "$mbd_name" =~ ^MBD[[:space:]]+N\&S[[:space:]]+\>\=[[:space:]]+([0-9]+)$ ]]; then
            mbd_thr="${BASH_REMATCH[1]}"

            for pho_bit in "${!epoch_menu_ref[@]}"; do
              pho_name="${epoch_menu_ref["$pho_bit"]:-}"

              if [[ "$pho_name" =~ ^Photon[[:space:]]+([0-9]+([.][0-9]+)?)[[:space:]]+GeV[[:space:]]*\+[[:space:]]+MBD[[:space:]]+N\&?S[[:space:]]+\>\=[[:space:]]+([0-9]+)$ ]]; then
                pho_e="${BASH_REMATCH[1]}"
                pho_thr="${BASH_REMATCH[3]}"

                if [[ "$pho_thr" == "$mbd_thr" ]]; then
                  novtx_rows+=( "$mbd_thr|$pho_e|$mbd_bit|$mbd_name|$pho_bit|$pho_name" )
                fi
              fi
            done
          elif [[ "$mbd_name" =~ ^MBD[[:space:]]+N\&S[[:space:]]+\>\=[[:space:]]+([0-9]+),[[:space:]]+vtx[[:space:]]+\<[[:space:]]+([0-9]+([.][0-9]+)?)[[:space:]]+cm$ ]]; then
            mbd_thr="${BASH_REMATCH[1]}"
            mbd_vtx="${BASH_REMATCH[2]}"

            for pho_bit in "${!epoch_menu_ref[@]}"; do
              pho_name="${epoch_menu_ref["$pho_bit"]:-}"

              if [[ "$pho_name" =~ ^Photon[[:space:]]+([0-9]+([.][0-9]+)?)[[:space:]]+GeV([[:space:]]*\+|,)[[:space:]]+MBD[[:space:]]+N\&?S[[:space:]]+\>\=[[:space:]]+([0-9]+),[[:space:]]+vtx[[:space:]]+\<[[:space:]]+([0-9]+([.][0-9]+)?)[[:space:]]+cm$ ]]; then
                pho_e="${BASH_REMATCH[1]}"
                pho_thr="${BASH_REMATCH[4]}"
                pho_vtx="${BASH_REMATCH[5]}"

                if [[ "$pho_thr" == "$mbd_thr" && "$pho_vtx" == "$mbd_vtx" ]]; then
                  vtx_rows+=( "$mbd_thr|$mbd_vtx|$pho_e|$mbd_bit|$mbd_name|$pho_bit|$pho_name" )
                fi
              fi
            done
          fi
        done

        echo "Epoch $epoch_idx"
        echo "  Runs: $run_start -> $run_end"

        echo "  No vertex cut combinations"
        if ((${#novtx_rows[@]})); then
          printf "%s\n" "${novtx_rows[@]}" | sort -t'|' -k1,1n -k2,2n | while IFS='|' read -r _thr _e mbd_bit mbd_name pho_bit pho_name; do
            printf "    MBD bit %-3s %-35s || Photon bit %-3s %s\n" "$mbd_bit" "$mbd_name" "$pho_bit" "$pho_name"
            printf "    CONFIG runMin=%s runMax=%s baselineBit=%s probeBit=%s baselineName=\"%s\" probeName=\"%s\"\n" \
              "$run_start" "$run_end" "$mbd_bit" "$pho_bit" "$mbd_name" "$pho_name"
          done
        else
          echo "    (none)"
        fi

        echo "  Vertex cut combinations"
        if ((${#vtx_rows[@]})); then
          printf "%s\n" "${vtx_rows[@]}" | sort -t'|' -k1,1n -k2,2n -k3,3n | while IFS='|' read -r _thr _vtx _e mbd_bit mbd_name pho_bit pho_name; do
            printf "    MBD bit %-3s %-35s || Photon bit %-3s %s\n" "$mbd_bit" "$mbd_name" "$pho_bit" "$pho_name"
            printf "    CONFIG runMin=%s runMax=%s baselineBit=%s probeBit=%s baselineName=\"%s\" probeName=\"%s\"\n" \
              "$run_start" "$run_end" "$mbd_bit" "$pho_bit" "$mbd_name" "$pho_name"
          done
        else
          echo "    (none)"
        fi

        echo
      }

      {
        echo "Trigger groupings for $LABEL"
        echo "Input golden run list: $LIST_FILE"
        echo "Runs in input list: ${#RUNS_ALL[@]}"
        echo "Trigger-menu epochs observed: $(( menu_change_count + 1 ))"
        if [[ -n "$first_menu_change_run" ]]; then
          echo "First adjacent-run menu change: ${first_menu_change_prev_run} -> ${first_menu_change_run}"
        else
          echo "First adjacent-run menu change: none"
        fi
        echo

        declare -A EPOCH_MENU=()

        for run in "${RUNS_ALL[@]}"; do
          unset CUR_MENU
          declare -A CUR_MENU=()

          while IFS=$'\t' read -r idx trg; do
            [[ -z "$idx" ]] && continue
            idx=$(num_or_zero "$idx")
            [[ "$idx" =~ ^[0-9]+$ ]] || continue
            [[ -z "$trg" ]] && trg="(missing-name)"
            CUR_MENU["$idx"]="$trg"
          done < <(sql "SELECT index, triggername
                        FROM gl1_triggernames
                        WHERE $run BETWEEN runnumber AND runnumber_last
                        ORDER BY index;")

          if [[ -z "$epoch_start" ]]; then
            epoch_start="$run"
            unset EPOCH_MENU
            declare -A EPOCH_MENU=()
            for idx in "${!CUR_MENU[@]}"; do
              EPOCH_MENU["$idx"]="${CUR_MENU["$idx"]}"
            done
          else
            unset MENU_KEYS
            declare -A MENU_KEYS=()

            for idx in "${!EPOCH_MENU[@]}"; do
              MENU_KEYS["$idx"]=1
            done

            for idx in "${!CUR_MENU[@]}"; do
              MENU_KEYS["$idx"]=1
            done

            local epoch_changed=0
            for idx in "${!MENU_KEYS[@]}"; do
              if [[ "${EPOCH_MENU["$idx"]:-"(absent)"}" != "${CUR_MENU["$idx"]:-"(absent)"}" ]]; then
                epoch_changed=1
                break
              fi
            done

            if (( epoch_changed )); then
              emit_epoch_groupings "$epoch_start" "$prev_run" EPOCH_MENU
              epoch_start="$run"
              unset EPOCH_MENU
              declare -A EPOCH_MENU=()
              for idx in "${!CUR_MENU[@]}"; do
                EPOCH_MENU["$idx"]="${CUR_MENU["$idx"]}"
              done
            fi
          fi

          prev_run="$run"
        done

        if [[ -n "$epoch_start" ]]; then
          emit_epoch_groupings "$epoch_start" "$prev_run" EPOCH_MENU
        fi
      } | tee "$outfile"

      echo "[INFO] Trigger grouping summary saved to: $outfile"
      echo
    }
    printf "%sDetailed bit-level audit helpers for %s (reference only)%s\n" "$ANSI_BOLD$ANSI_CYAN" "$LABEL" "$ANSI_RESET"
    echo
    print_focus_table "No vertex cut focus table" 10 24 25 26 27
    print_focus_table "Vertex cut focus table" 12 36 37 38
    print_trigger_groupings

print_compact_activation_summary() {
      local csv_file="$1"

      local -A FAMILY_BASE_ACTIVE=()
      local -A FAMILY_COMBO_COUNTS=()
      local -A FAMILY_PHOTON_PRESENT=()
      local -A FAMILY_FIRST_SEEN=()
      local -A GLOBAL_PHOTON_SEEN=()
      local -A GLOBAL_COMBO_SEEN=()
      local -a FAMILY_ORDER=()

      local run trg idx scaled raw live family base_thr base_vtx pho_e pho_label

      compact_family_label() {
        local s="$1"
        s="${s/MBD N&S >= /MBD>=}"
        s="${s/, vtx < /,vz<}"
        s="${s/ cm/}"
        printf "%s" "$s"
      }

      compact_combo_label() {
        local s="$1"
        if [[ -z "$s" || "$s" == "(none)" ]]; then
          printf -- "-"
          return
        fi
        if [[ "$s" == "(baseline only)" ]]; then
          printf "base"
          return
        fi

        local out=""
        local item num
        IFS=',' read -ra _items <<< "${s//, /,}"
        for item in "${_items[@]}"; do
          item="${item#"${item%%[![:space:]]*}"}"
          item="${item%"${item##*[![:space:]]}"}"
          [[ -z "$item" ]] && continue
          num="${item#Pho}"
          if [[ -n "$out" ]]; then
            out="${out}+${num}"
          else
            out="${num}"
          fi
        done
        [[ -z "$out" ]] && out="-"
        printf "%s" "$out"
      }

      compact_pho_set_label() {
        local s="$1"
        if [[ -z "$s" || "$s" == "(none)" ]]; then
          printf -- "-"
          return
        fi

        local out=""
        local item num
        IFS=',' read -ra _items <<< "${s//, /,}"
        for item in "${_items[@]}"; do
          item="${item#"${item%%[![:space:]]*}"}"
          item="${item%"${item##*[![:space:]]}"}"
          [[ -z "$item" ]] && continue
          num="${item#Pho}"
          if [[ -n "$out" ]]; then
            out="${out},${num}"
          else
            out="${num}"
          fi
        done
        [[ -z "$out" ]] && out="-"
        printf "%s" "$out"
      }

      hrule() {
        local n="$1"
        printf "%*s" "$n" "" | tr ' ' '-'
      }

      for run in "${RUNS_ALL[@]}"; do
        unset RUN_BASE_ACTIVE RUN_FAMILY_PHOTONS RUN_FAMILY_PHOTON_SEEN
        declare -A RUN_BASE_ACTIVE=()
        declare -A RUN_FAMILY_PHOTONS=()
        declare -A RUN_FAMILY_PHOTON_SEEN=()

        while IFS=$'\t' read -r trg idx scaled raw live; do
          [[ -z "$idx" ]] && continue

          idx=$(num_or_zero "$idx")
          scaled=$(num_or_zero "$scaled")

          [[ "$idx" =~ ^[0-9]+$ ]] || continue
          (( scaled != -1 )) || continue

          [[ -z "$trg" ]] && trg="${TRIG_NAME["$idx"]:-"(missing-name)"}"

          if [[ "$trg" =~ ^MBD[[:space:]]+N\&S[[:space:]]+\>\=[[:space:]]+([0-9]+)(,[[:space:]]+vtx[[:space:]]+\<[[:space:]]+([0-9]+([.][0-9]+)?)[[:space:]]+cm)?$ ]]; then
            base_thr="${BASH_REMATCH[1]}"
            base_vtx="${BASH_REMATCH[3]:-}"

            family="MBD N&S >= ${base_thr}"
            [[ -n "$base_vtx" ]] && family+=", vtx < ${base_vtx} cm"

            RUN_BASE_ACTIVE["$family"]=1

            if [[ -z "${FAMILY_FIRST_SEEN["$family"]:-}" ]]; then
              FAMILY_FIRST_SEEN["$family"]="$run"
              FAMILY_ORDER+=( "$family" )
            fi

          elif [[ "$trg" =~ ^Photon[[:space:]]+([0-9]+([.][0-9]+)?)[[:space:]]+GeV([[:space:]]*\+|,)[[:space:]]+MBD[[:space:]]+N\&?S[[:space:]]+\>\=[[:space:]]+([0-9]+)(,[[:space:]]+vtx[[:space:]]+\<[[:space:]]+([0-9]+([.][0-9]+)?)[[:space:]]+cm)?$ ]]; then
            pho_e="${BASH_REMATCH[1]}"
            base_thr="${BASH_REMATCH[4]}"
            base_vtx="${BASH_REMATCH[6]:-}"

            family="MBD N&S >= ${base_thr}"
            [[ -n "$base_vtx" ]] && family+=", vtx < ${base_vtx} cm"

            pho_label="Pho${pho_e}"

            if [[ -z "${FAMILY_FIRST_SEEN["$family"]:-}" ]]; then
              FAMILY_FIRST_SEEN["$family"]="$run"
              FAMILY_ORDER+=( "$family" )
            fi

            GLOBAL_PHOTON_SEEN["$pho_label"]=1

            if [[ -z "${RUN_FAMILY_PHOTON_SEEN["$family|$pho_label"]:-}" ]]; then
              if [[ -n "${RUN_FAMILY_PHOTONS["$family"]:-}" ]]; then
                RUN_FAMILY_PHOTONS["$family"]+=", ${pho_label}"
              else
                RUN_FAMILY_PHOTONS["$family"]="${pho_label}"
              fi
              RUN_FAMILY_PHOTON_SEEN["$family|$pho_label"]=1
            fi
          fi
        done < <(sql "SELECT COALESCE(t.triggername,''), s.index, s.scaled, s.raw, s.live
                       FROM gl1_scalers s
                       LEFT JOIN gl1_triggernames t
                         ON s.index = t.index
                        AND s.runnumber BETWEEN t.runnumber AND t.runnumber_last
                      WHERE s.runnumber = $run
                      ORDER BY s.index;")

        local combo npho key
        for family in "${!RUN_BASE_ACTIVE[@]}"; do
          FAMILY_BASE_ACTIVE["$family"]=$(( ${FAMILY_BASE_ACTIVE["$family"]:-0} + 1 ))

          combo="${RUN_FAMILY_PHOTONS["$family"]:-"(baseline only)"}"
          if [[ "$combo" == "(baseline only)" ]]; then
            npho=0
          else
            npho=$(awk -v s="$combo" 'BEGIN{ n=split(s,a,/, */); print n+0 }')
            IFS=',' read -ra _labels <<< "${combo//, /,}"
            for pho_label in "${_labels[@]}"; do
              pho_label="${pho_label#"${pho_label%%[![:space:]]*}"}"
              pho_label="${pho_label%"${pho_label##*[![:space:]]}"}"
              [[ -z "$pho_label" ]] && continue
              FAMILY_PHOTON_PRESENT["$family|$pho_label"]=1
            done
            GLOBAL_COMBO_SEEN["$combo"]=1
          fi

          key="${family}|${npho}|${combo}"
          FAMILY_COMBO_COUNTS["$key"]=$(( ${FAMILY_COMBO_COUNTS["$key"]:-0} + 1 ))
        done
      done

      local -a GLOBAL_PHOTON_LABELS=()
      if ((${#GLOBAL_PHOTON_SEEN[@]})); then
        mapfile -t GLOBAL_PHOTON_LABELS < <(
          for pho_label in "${!GLOBAL_PHOTON_SEEN[@]}"; do
            printf "%s\n" "$pho_label"
          done | awk '{x=$0; gsub(/^Pho/,"",x); printf "%08.3f|%s\n", x+0, $0}' | sort -n | cut -d'|' -f2-
        )
      fi

      local -a GLOBAL_COMBO_KEYS=()
      local -a GLOBAL_COMBO_LABELS=()
      if ((${#GLOBAL_COMBO_SEEN[@]})); then
        mapfile -t _combo_rows < <(
          for combo in "${!GLOBAL_COMBO_SEEN[@]}"; do
            npho=$(awk -v s="$combo" 'BEGIN{ n=split(s,a,/, */); print n+0 }')
            printf "%04d|%s|%s\n" "$npho" "$(compact_combo_label "$combo")" "$combo"
          done | sort -t'|' -k1,1n -k2,2
        )

        local _row _n _label _combo
        for _row in "${_combo_rows[@]}"; do
          IFS='|' read -r _n _label _combo <<< "$_row"
          GLOBAL_COMBO_LABELS+=( "$_label" )
          GLOBAL_COMBO_KEYS+=( "$_combo" )
        done
      fi

      printf "BaselineFamily,DisplayLabel,ActiveRuns,TotalGRLRuns,InactiveRuns,PhoSet,BaseOnlyRuns" > "$csv_file"
      local pho_label combo_label
      for pho_label in "${GLOBAL_PHOTON_LABELS[@]}"; do
        printf ",Pho%s" "${pho_label#Pho}" >> "$csv_file"
      done
      for combo_label in "${GLOBAL_COMBO_LABELS[@]}"; do
        printf ",Combo_%s" "${combo_label//+/_}" >> "$csv_file"
      done
      printf "\n" >> "$csv_file"

      local family_w=28 active_w=10 inactive_w=8 phoset_w=16 baseonly_w=8
      local -a combo_widths=()
      for i in "${!GLOBAL_COMBO_LABELS[@]}"; do
        local cw=${#GLOBAL_COMBO_LABELS[$i]}
        (( cw < 6 )) && cw=6
        combo_widths+=( "$cw" )
      done

      printf "\n%s═══════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════%s\n" "$ANSI_BOLD$ANSI_CYAN" "$ANSI_RESET"
      printf "%s  CANONICAL BASELINE + PHOTON ACTIVATION SNAPSHOT — %s (%d GRL runs)%s\n" "$ANSI_BOLD$ANSI_CYAN" "$LABEL" "${#RUNS_ALL[@]}" "$ANSI_RESET"
      printf "%s═══════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════%s\n" "$ANSI_BOLD$ANSI_CYAN" "$ANSI_RESET"

      printf "\n%s  %-*s  %*s  %*s  %-*s  %*s" \
        "$ANSI_BOLD" \
        "$family_w" "Baseline" \
        "$active_w" "Active/GRL" \
        "$inactive_w" "Inactive" \
        "$phoset_w" "PhoSet" \
        "$baseonly_w" "BaseOnly"
      for i in "${!GLOBAL_COMBO_LABELS[@]}"; do
        printf "  %*s" "${combo_widths[$i]}" "${GLOBAL_COMBO_LABELS[$i]}"
      done
      printf "%s\n" "$ANSI_RESET"

      printf "  "
      hrule "$family_w"
      printf "  "
      hrule "$active_w"
      printf "  "
      hrule "$inactive_w"
      printf "  "
      hrule "$phoset_w"
      printf "  "
      hrule "$baseonly_w"
      for i in "${!GLOBAL_COMBO_LABELS[@]}"; do
        printf "  "
        hrule "${combo_widths[$i]}"
      done
      printf "\n"

      local base_active inactive_runs base_only_count combo_count combo_key
      for family in "${FAMILY_ORDER[@]}"; do
        display_family="$(compact_family_label "$family")"
        base_active="${FAMILY_BASE_ACTIVE["$family"]:-0}"
        inactive_runs=$(( ${#RUNS_ALL[@]} - base_active ))
        base_only_count="${FAMILY_COMBO_COUNTS["$family|0|(baseline only)"]:-0}"

        declare -A _family_pho_seen=()
        local -a family_pho_labels=()
        local k label
        for k in "${!FAMILY_PHOTON_PRESENT[@]}"; do
          [[ "$k" == "$family|"* ]] || continue
          label="${k#"$family|"}"
          _family_pho_seen["$label"]=1
        done
        if ((${#_family_pho_seen[@]})); then
          mapfile -t family_pho_labels < <(
            for label in "${!_family_pho_seen[@]}"; do
              printf "%s\n" "$label"
            done | awk '{x=$0; gsub(/^Pho/,"",x); printf "%08.3f|%s\n", x+0, $0}' | sort -n | cut -d'|' -f2-
          )
          phoset_display="$(compact_pho_set_label "$(IFS=', '; echo "${family_pho_labels[*]}")")"
        else
          phoset_display="-"
        fi

        local _active_str="${base_active}/${#RUNS_ALL[@]}"
        printf "  %-*s  %*s  %*d  %-*s  %*d" \
          "$family_w" "$display_family" \
          "$active_w" "$_active_str" \
          "$inactive_w" "$inactive_runs" \
          "$phoset_w" "$phoset_display" \
          "$baseonly_w" "$base_only_count"

        for i in "${!GLOBAL_COMBO_KEYS[@]}"; do
          combo_key="${GLOBAL_COMBO_KEYS[$i]}"
          combo_count="${FAMILY_COMBO_COUNTS["$family|$(awk -v s="$combo_key" 'BEGIN{ n=split(s,a,/, */); print n+0 }')|$combo_key"]:-0}"
          if (( combo_count > 0 )); then
            printf "  %*d" "${combo_widths[$i]}" "$combo_count"
          else
            printf "  %*s" "${combo_widths[$i]}" "-"
          fi
        done
        printf "\n"

        printf "\"%s\",\"%s\",%d,%d,%d,\"%s\",%d" \
          "$family" \
          "$display_family" \
          "$base_active" \
          "${#RUNS_ALL[@]}" \
          "$inactive_runs" \
          "$phoset_display" \
          "$base_only_count" >> "$csv_file"

        for i in "${!GLOBAL_COMBO_KEYS[@]}"; do
          combo_key="${GLOBAL_COMBO_KEYS[$i]}"
          combo_count="${FAMILY_COMBO_COUNTS["$family|$(awk -v s="$combo_key" 'BEGIN{ n=split(s,a,/, */); print n+0 }')|$combo_key"]:-0}"
          printf ",%d" "$combo_count" >> "$csv_file"
        done
        printf "\n" >> "$csv_file"
      done

      printf "\n%s  CSV saved to: %s%s\n\n" "$ANSI_DIM" "$csv_file" "$ANSI_RESET"
    }
    # ── end canonical baseline + photon activation summary ──────────────────────

    generate_scaled_efficiency_study_tables() {
      local txt_file="$1"
      local csv_file="$2"

      if [[ "${LABEL:-}" != "auau" ]]; then
        return 0
      fi

      local -a study_order=(
        "MBD_NS_geq_2_vtx_lt_150__Pho6_8_10_12"
        "MBD_NS_geq_2_vtx_lt_10__Pho6_8_10_12"
      )

      local -A study_family=(
        ["MBD_NS_geq_2_vtx_lt_150__Pho6_8_10_12"]="MBD N&S >= 2, vtx < 150 cm"
        ["MBD_NS_geq_2_vtx_lt_10__Pho6_8_10_12"]="MBD N&S >= 2, vtx < 10 cm"
      )

      local -A study_baseline_key=(
        ["MBD_NS_geq_2_vtx_lt_150__Pho6_8_10_12"]="MBD_NS_geq_2_vtx_lt_150"
        ["MBD_NS_geq_2_vtx_lt_10__Pho6_8_10_12"]="MBD_NS_geq_2_vtx_lt_10"
      )

      local -A study_list_file=(
        ["MBD_NS_geq_2_vtx_lt_150__Pho6_8_10_12"]="$OUT_DIR/scaledEffRuns_MBD_NS_geq_2_vtx_lt_150__Pho6_8_10_12.list"
        ["MBD_NS_geq_2_vtx_lt_10__Pho6_8_10_12"]="$OUT_DIR/scaledEffRuns_MBD_NS_geq_2_vtx_lt_10__Pho6_8_10_12.list"
      )

      csv_quote() {
        local s="$1"
        s="${s//\"/\"\"}"
        printf '"%s"' "$s"
      }

      pho_key_for_study() {
        local study="$1"
        local pho="$2"
        case "$study" in
          MBD_NS_geq_2_vtx_lt_150__Pho6_8_10_12)
            printf "photon_%s_plus_MBD_NS_geq_2_vtx_lt_150" "$pho"
            ;;
          MBD_NS_geq_2_vtx_lt_10__Pho6_8_10_12)
            printf "photon_%s_plus_MBD_NS_geq_2_vtx_lt_10" "$pho"
            ;;
          *)
            printf "photon_%s_UNKNOWN" "$pho"
            ;;
        esac
      }

      mkdir -p "$OUT_DIR"

      {
        echo "# Scaled trigger-efficiency run table for $LABEL"
        echo "# Input golden run list: $LIST_FILE"
        echo "# Produced by: ./make_dstListsData.sh $MODE QA"
        echo "# Rule: each CONFIG row requires baseline + Photon 6/8/10/12 to have scaled > 0 in the same run."
        echo "# These rows are intended only for scaled max-cluster-energy trigger-efficiency studies."
        echo "# Format:"
        echo "# CONFIG study=<studyKey> run=<run> baselineBit=<bit> baselineKey=<key> baselineScale=<live/scaled> pho6Bit=<bit> pho6Key=<key> pho6Scale=<live/scaled> ..."
      } > "$txt_file"

      printf "StudyGroup,Run,BaselineFamily,BaselineBit,BaselineKey,BaselineName,BaselineLiveOverScaled" > "$csv_file"
      printf ",Pho6Bit,Pho6Key,Pho6Name,Pho6LiveOverScaled" >> "$csv_file"
      printf ",Pho8Bit,Pho8Key,Pho8Name,Pho8LiveOverScaled" >> "$csv_file"
      printf ",Pho10Bit,Pho10Key,Pho10Name,Pho10LiveOverScaled" >> "$csv_file"
      printf ",Pho12Bit,Pho12Key,Pho12Name,Pho12LiveOverScaled\n" >> "$csv_file"

      local study
      for study in "${study_order[@]}"; do
        : > "${study_list_file["$study"]}"
      done

      local run trg idx scaled raw live live_scaled
      local base_thr base_vtx family pho_e pho_int
      local study_count

      for run in "${RUNS_ALL[@]}"; do
        unset RUN_BASE_BIT RUN_BASE_NAME RUN_BASE_SCALE
        unset RUN_PHO_BIT RUN_PHO_NAME RUN_PHO_SCALE
        declare -A RUN_BASE_BIT=()
        declare -A RUN_BASE_NAME=()
        declare -A RUN_BASE_SCALE=()
        declare -A RUN_PHO_BIT=()
        declare -A RUN_PHO_NAME=()
        declare -A RUN_PHO_SCALE=()

        while IFS=$'\t' read -r trg idx scaled raw live; do
          [[ -z "$idx" ]] && continue

          idx=$(num_or_zero "$idx")
          scaled=$(num_or_zero "$scaled")
          raw=$(num_or_zero "$raw")
          live=$(num_or_zero "$live")

          [[ "$idx" =~ ^[0-9]+$ ]] || continue
          (( scaled > 0 )) || continue

          [[ -z "$trg" ]] && trg="${TRIG_NAME["$idx"]:-"(missing-name)"}"

          live_scaled=$(awk -v l="$live" -v s="$scaled" 'BEGIN{ if (s>0) printf "%.6f", l/s; else printf "0.000000"; }')

          if [[ "$trg" =~ ^MBD[[:space:]]+N\&S[[:space:]]+\>\=[[:space:]]+([0-9]+)(,[[:space:]]+vtx[[:space:]]+\<[[:space:]]+([0-9]+([.][0-9]+)?)[[:space:]]+cm)?$ ]]; then
            base_thr="${BASH_REMATCH[1]}"
            base_vtx="${BASH_REMATCH[3]:-}"

            family="MBD N&S >= ${base_thr}"
            [[ -n "$base_vtx" ]] && family+=", vtx < ${base_vtx} cm"

            RUN_BASE_BIT["$family"]="$idx"
            RUN_BASE_NAME["$family"]="$trg"
            RUN_BASE_SCALE["$family"]="$live_scaled"

          elif [[ "$trg" =~ ^Photon[[:space:]]+([0-9]+([.][0-9]+)?)[[:space:]]+GeV([[:space:]]*\+|,)[[:space:]]+MBD[[:space:]]+N\&?S[[:space:]]+\>\=[[:space:]]+([0-9]+)(,[[:space:]]+vtx[[:space:]]+\<[[:space:]]+([0-9]+([.][0-9]+)?)[[:space:]]+cm)?$ ]]; then
            pho_e="${BASH_REMATCH[1]}"
            base_thr="${BASH_REMATCH[4]}"
            base_vtx="${BASH_REMATCH[6]:-}"

            pho_int=$(awk -v x="$pho_e" 'BEGIN{ printf "%g", x+0 }')

            family="MBD N&S >= ${base_thr}"
            [[ -n "$base_vtx" ]] && family+=", vtx < ${base_vtx} cm"

            RUN_PHO_BIT["$family|$pho_int"]="$idx"
            RUN_PHO_NAME["$family|$pho_int"]="$trg"
            RUN_PHO_SCALE["$family|$pho_int"]="$live_scaled"
          fi
        done < <(sql "SELECT COALESCE(t.triggername,''), s.index, s.scaled, s.raw, s.live
                       FROM gl1_scalers s
                       LEFT JOIN gl1_triggernames t
                         ON s.index = t.index
                        AND s.runnumber BETWEEN t.runnumber AND t.runnumber_last
                      WHERE s.runnumber = $run
                      ORDER BY s.index;")

        for study in "${study_order[@]}"; do
          family="${study_family["$study"]}"

          [[ -n "${RUN_BASE_BIT["$family"]:-}" ]] || continue
          [[ -n "${RUN_PHO_BIT["$family|6"]:-}" ]] || continue
          [[ -n "${RUN_PHO_BIT["$family|8"]:-}" ]] || continue
          [[ -n "${RUN_PHO_BIT["$family|10"]:-}" ]] || continue
          [[ -n "${RUN_PHO_BIT["$family|12"]:-}" ]] || continue

          local baseline_key="${study_baseline_key["$study"]}"
          local pho6_key pho8_key pho10_key pho12_key
          pho6_key="$(pho_key_for_study "$study" 6)"
          pho8_key="$(pho_key_for_study "$study" 8)"
          pho10_key="$(pho_key_for_study "$study" 10)"
          pho12_key="$(pho_key_for_study "$study" 12)"

          printf "CONFIG study=%s run=%s baselineBit=%s baselineKey=%s baselineScale=%s" \
            "$study" \
            "$run" \
            "${RUN_BASE_BIT["$family"]}" \
            "$baseline_key" \
            "${RUN_BASE_SCALE["$family"]}" >> "$txt_file"
          printf " pho6Bit=%s pho6Key=%s pho6Scale=%s" \
            "${RUN_PHO_BIT["$family|6"]}" \
            "$pho6_key" \
            "${RUN_PHO_SCALE["$family|6"]}" >> "$txt_file"
          printf " pho8Bit=%s pho8Key=%s pho8Scale=%s" \
            "${RUN_PHO_BIT["$family|8"]}" \
            "$pho8_key" \
            "${RUN_PHO_SCALE["$family|8"]}" >> "$txt_file"
          printf " pho10Bit=%s pho10Key=%s pho10Scale=%s" \
            "${RUN_PHO_BIT["$family|10"]}" \
            "$pho10_key" \
            "${RUN_PHO_SCALE["$family|10"]}" >> "$txt_file"
          printf " pho12Bit=%s pho12Key=%s pho12Scale=%s\n" \
            "${RUN_PHO_BIT["$family|12"]}" \
            "$pho12_key" \
            "${RUN_PHO_SCALE["$family|12"]}" >> "$txt_file"

          csv_quote "$study" >> "$csv_file"
          printf ",%s," "$run" >> "$csv_file"
          csv_quote "$family" >> "$csv_file"
          printf ",%s," "${RUN_BASE_BIT["$family"]}" >> "$csv_file"
          csv_quote "$baseline_key" >> "$csv_file"
          printf "," >> "$csv_file"
          csv_quote "${RUN_BASE_NAME["$family"]}" >> "$csv_file"
          printf ",%s" "${RUN_BASE_SCALE["$family"]}" >> "$csv_file"

          printf ",%s," "${RUN_PHO_BIT["$family|6"]}" >> "$csv_file"
          csv_quote "$pho6_key" >> "$csv_file"
          printf "," >> "$csv_file"
          csv_quote "${RUN_PHO_NAME["$family|6"]}" >> "$csv_file"
          printf ",%s" "${RUN_PHO_SCALE["$family|6"]}" >> "$csv_file"

          printf ",%s," "${RUN_PHO_BIT["$family|8"]}" >> "$csv_file"
          csv_quote "$pho8_key" >> "$csv_file"
          printf "," >> "$csv_file"
          csv_quote "${RUN_PHO_NAME["$family|8"]}" >> "$csv_file"
          printf ",%s" "${RUN_PHO_SCALE["$family|8"]}" >> "$csv_file"

          printf ",%s," "${RUN_PHO_BIT["$family|10"]}" >> "$csv_file"
          csv_quote "$pho10_key" >> "$csv_file"
          printf "," >> "$csv_file"
          csv_quote "${RUN_PHO_NAME["$family|10"]}" >> "$csv_file"
          printf ",%s" "${RUN_PHO_SCALE["$family|10"]}" >> "$csv_file"

          printf ",%s," "${RUN_PHO_BIT["$family|12"]}" >> "$csv_file"
          csv_quote "$pho12_key" >> "$csv_file"
          printf "," >> "$csv_file"
          csv_quote "${RUN_PHO_NAME["$family|12"]}" >> "$csv_file"
          printf ",%s\n" "${RUN_PHO_SCALE["$family|12"]}" >> "$csv_file"

          printf "%s\n" "$run" >> "${study_list_file["$study"]}"
        done
      done

      printf "\n%sScaled trigger-efficiency study files for %s%s\n" "$ANSI_BOLD$ANSI_CYAN" "$LABEL" "$ANSI_RESET"
      printf "  CONFIG table: %s\n" "$txt_file"
      printf "  CSV table   : %s\n" "$csv_file"

      for study in "${study_order[@]}"; do
        study_count=$(awk 'NF{n++} END{print n+0}' "${study_list_file["$study"]}")
        printf "  Run list    : %s  (%d runs)\n" "${study_list_file["$study"]}" "$study_count"
      done
      printf "\n"
    }

    local activation_csv="$OUT_DIR/trigger_activation_summary_${LABEL}.csv"
    print_compact_activation_summary "$activation_csv"

    local scaled_eff_txt="$OUT_DIR/trigger_scaled_efficiency_studies_${LABEL}.txt"
    local scaled_eff_csv="$OUT_DIR/trigger_scaled_efficiency_studies_${LABEL}.csv"
    generate_scaled_efficiency_study_tables "$scaled_eff_txt" "$scaled_eff_csv"

  else
    printf "%s[WARN] No trigger rows were found for the provided golden run list.%s\n" "$ANSI_YELLOW" "$ANSI_RESET"
    printf "%s[WARN] This usually means the queries returned no menu/scaler rows for the supplied run list.%s\n" "$ANSI_YELLOW" "$ANSI_RESET"
  fi

  if (( runs_with_no_menu_rows > 0 || runs_with_no_scaler_rows > 0 || rows_bad_idx > 0 )); then
    printf "%s[WARN] QA completed with transparency warnings. Review the trigger-menu and raw/live diagnostics above before trusting the tables blindly.%s\n" "$ANSI_YELLOW" "$ANSI_RESET"
  else
    printf "%s[INFO] QA completed without structural warning conditions.%s\n" "$ANSI_GREEN" "$ANSI_RESET"
  fi
}



run_scaled_trigger_ana() {
  command -v psql >/dev/null 2>&1 || { echo "[ERROR] psql not found in PATH"; exit 1; }

  if [[ "${LABEL:-}" != "auau" ]]; then
    echo "[ERROR] scaledTriggerAna is currently defined only for auau mode."
    echo "        Run: ./make_dstListsData.sh auau QA scaledTriggerAna"
    return 1
  fi

  [[ -f "$LIST_FILE" ]] || { echo "[ERROR] Run list not found: $LIST_FILE"; return 1; }
  [[ -s "$LIST_FILE" ]] || { echo "[ERROR] Run list is empty: $LIST_FILE"; return 1; }

  read_runs_from_file "$LIST_FILE" RUNS_ALL
  (( ${#RUNS_ALL[@]} > 0 )) || { echo "[ERROR] No valid run numbers found in: $LIST_FILE"; return 1; }

  PSQL=(psql -h sphnxdaqdbreplica -d daq -At -F $'\t' -q)
  sql()         { "${PSQL[@]}" -c "$1" 2>/dev/null || true; }
  num_or_zero() { [[ $1 =~ ^-?[0-9]+([.][0-9]+)?$ ]] && printf '%s' "$1" || printf 0; }

  local ANSI_BOLD ANSI_DIM ANSI_CYAN ANSI_GREEN ANSI_YELLOW ANSI_RED ANSI_RESET
  ANSI_BOLD=$'\033[1m'
  ANSI_DIM=$'\033[2m'
  ANSI_CYAN=$'\033[36m'
  ANSI_GREEN=$'\033[32m'
  ANSI_YELLOW=$'\033[33m'
  ANSI_RED=$'\033[31m'
  ANSI_RESET=$'\033[0m'

  mkdir -p "$OUT_DIR"

  local selected_study="MBD_NS_geq_2_vtx_lt_150__Pho10_12"
  local selected_family="MBD N&S >= 2, vtx < 150 cm"
  local selected_baseline_key="MBD_NS_geq_2_vtx_lt_150"
  local selected_runlist="$OUT_DIR/scaledEffRuns_${selected_study}.list"

  local txt_file="$OUT_DIR/trigger_scaled_efficiency_studies_${LABEL}.txt"
  local run_csv_file="$OUT_DIR/trigger_scaled_efficiency_run_diagnostics_${LABEL}.csv"
  local summary_csv_file="$OUT_DIR/trigger_scaled_efficiency_summary_${LABEL}.csv"
  local reject_file="$OUT_DIR/trigger_scaled_efficiency_rejections_${LABEL}.txt"

  rm -f "$OUT_DIR/scaledEffRuns_MBD_NS_geq_2_vtx_lt_150__Pho6_8_10_12.list"
  rm -f "$OUT_DIR/scaledEffRuns_MBD_NS_geq_2_vtx_lt_10__Pho6_8_10_12.list"
  : > "$selected_runlist"

  csv_quote() {
    local s="$1"
    s="${s//\"/\"\"}"
    printf '"%s"' "$s"
  }

  pho_key_for_family() {
    local family="$1"
    local pho="$2"

    case "$family" in
      "MBD N&S >= 2, vtx < 150 cm")
        printf "photon_%s_plus_MBD_NS_geq_2_vtx_lt_150" "$pho"
        ;;
      "MBD N&S >= 2, vtx < 10 cm")
        printf "photon_%s_plus_MBD_NS_geq_2_vtx_lt_10" "$pho"
        ;;
      *)
        printf "photon_%s_UNKNOWN" "$pho"
        ;;
    esac
  }

  join_reason_parts() {
    local out=""
    local part
    for part in "$@"; do
      if [[ -n "$out" ]]; then
        out="${out}; ${part}"
      else
        out="$part"
      fi
    done
    printf "%s" "$out"
  }

  {
    echo "# Scaled trigger-efficiency run table for $LABEL"
    echo "# Input golden run list: $LIST_FILE"
    echo "# Produced by: ./make_dstListsData.sh $MODE QA scaledTriggerAna"
    echo "# Accepted CONFIG rows require:"
    echo "#   - MBD N&S >= 2, vtx < 150 cm scaled > 0"
    echo "#   - Photon 10 GeV + MBD NS >= 2, vtx < 150 cm scaled > 0"
    echo "#   - Photon 12 GeV + MBD NS >= 2, vtx < 150 cm scaled > 0"
    echo "# This is the recommended common-run scaled study because Photon 6/8 have no scaled-positive sample."
    echo "# Run list for source code:"
    echo "#   $selected_runlist"
    echo "# Format:"
    echo "# CONFIG study=<studyKey> run=<run> baselineBit=<bit> baselineKey=<key> baselineScale=<live/scaled> pho10Bit=<bit> pho10Key=<key> pho10Scale=<live/scaled> pho12Bit=<bit> pho12Key=<key> pho12Scale=<live/scaled>"
  } > "$txt_file"

  {
    echo "# Rejection diagnostics for selected scaled trigger-efficiency study"
    echo "# Input golden run list: $LIST_FILE"
    echo "# Produced by: ./make_dstListsData.sh $MODE QA scaledTriggerAna"
    echo "# Selected study: $selected_study"
    echo "# A run is rejected from the selected run list if MBD vtx<150, Photon 10 vtx<150, or Photon 12 vtx<150 has scaled <= 0 or is missing."
    echo
  } > "$reject_file"

  printf "Run,SelectedStudy,Status,Reason" > "$run_csv_file"
  printf ",BaselineBit,BaselineKey,BaselineName,BaselineRaw,BaselineLive,BaselineScaled,BaselineLiveOverScaled" >> "$run_csv_file"
  printf ",Pho10Bit,Pho10Key,Pho10Name,Pho10Raw,Pho10Live,Pho10Scaled,Pho10LiveOverScaled" >> "$run_csv_file"
  printf ",Pho12Bit,Pho12Key,Pho12Name,Pho12Raw,Pho12Live,Pho12Scaled,Pho12LiveOverScaled\n" >> "$run_csv_file"

  printf "BaselineFamily,Photon,PhotonRows,RawPositiveRuns,LivePositiveRuns,ScaledPositiveRuns,RowButScaledZeroRuns,CommonWithBaselineScaledRuns,Decision\n" > "$summary_csv_file"

  declare -A FAMILY_BASE_ROWS=()
  declare -A FAMILY_BASE_RAW_POS=()
  declare -A FAMILY_BASE_LIVE_POS=()
  declare -A FAMILY_BASE_SCALED_POS=()

  declare -A PHO_ROW_RUNS=()
  declare -A PHO_RAW_POS=()
  declare -A PHO_LIVE_POS=()
  declare -A PHO_SCALED_POS=()
  declare -A PHO_COMMON_WITH_BASE_SCALED=()

  local SELECT_ACCEPT=0
  local SELECT_REJECT=0
  local REJ_BASE_MISSING=0
  local REJ_BASE_SCALED_ZERO=0
  local REJ_PHO10_MISSING=0
  local REJ_PHO10_SCALED_ZERO=0
  local REJ_PHO12_MISSING=0
  local REJ_PHO12_SCALED_ZERO=0
  local FIRST_SELECT_REJECTION=""

  printf "\n%sScaled trigger-efficiency run-list diagnostic%s\n" "$ANSI_BOLD$ANSI_CYAN" "$ANSI_RESET"
  echo "Input golden run list: $LIST_FILE"
  echo "Runs in input list: ${#RUNS_ALL[@]}"
  echo "Output dir: $OUT_DIR"
  echo
  echo "Goal:"
  echo "  1. Summarize all AuAu photon triggers attached to MBD N&S >= 2, vtx < 150 cm."
  echo "  2. Show how many runs have scaled > 0 for each photon trigger."
  echo "  3. Build the simplest valid scaled-efficiency source-code run list:"
  echo "       $selected_study"
  echo "     requiring MBD vtx<150 + Photon 10 + Photon 12 all scaled-positive in the same run."
  echo "  4. Also summarize why the MBD N&S >= 2, vtx < 10 cm subset is not useful for the same scaled-trigger study."
  echo

  local run trg idx scaled raw live live_scaled
  local base_thr base_vtx family pho_e pho_int key
  local processed=0

  for run in "${RUNS_ALL[@]}"; do
    ((processed+=1))

    unset RUN_BASE_BIT RUN_BASE_NAME RUN_BASE_RAW RUN_BASE_LIVE RUN_BASE_SCALED RUN_BASE_SCALE
    unset RUN_PHO_BIT RUN_PHO_NAME RUN_PHO_RAW RUN_PHO_LIVE RUN_PHO_SCALED RUN_PHO_SCALE
    declare -A RUN_BASE_BIT=()
    declare -A RUN_BASE_NAME=()
    declare -A RUN_BASE_RAW=()
    declare -A RUN_BASE_LIVE=()
    declare -A RUN_BASE_SCALED=()
    declare -A RUN_BASE_SCALE=()
    declare -A RUN_PHO_BIT=()
    declare -A RUN_PHO_NAME=()
    declare -A RUN_PHO_RAW=()
    declare -A RUN_PHO_LIVE=()
    declare -A RUN_PHO_SCALED=()
    declare -A RUN_PHO_SCALE=()

    while IFS=$'\t' read -r idx trg scaled raw live; do
      [[ -z "$idx" ]] && continue

      idx=$(num_or_zero "$idx")
      scaled=$(num_or_zero "$scaled")
      raw=$(num_or_zero "$raw")
      live=$(num_or_zero "$live")

      [[ "$idx" =~ ^[0-9]+$ ]] || continue

      [[ -z "$trg" ]] && trg="(missing-name)"
      live_scaled=$(awk -v l="$live" -v s="$scaled" 'BEGIN{ if (s>0) printf "%.6f", l/s; else printf "0.000000"; }')

      if [[ "$trg" =~ ^MBD[[:space:]]+N\&S[[:space:]]+\>\=[[:space:]]+([0-9]+)(,[[:space:]]+vtx[[:space:]]+\<[[:space:]]+([0-9]+([.][0-9]+)?)[[:space:]]+cm)?$ ]]; then
        base_thr="${BASH_REMATCH[1]}"
        base_vtx="${BASH_REMATCH[3]:-}"

        family="MBD N&S >= ${base_thr}"
        [[ -n "$base_vtx" ]] && family+=", vtx < ${base_vtx} cm"

        RUN_BASE_BIT["$family"]="$idx"
        RUN_BASE_NAME["$family"]="$trg"
        RUN_BASE_RAW["$family"]="$raw"
        RUN_BASE_LIVE["$family"]="$live"
        RUN_BASE_SCALED["$family"]="$scaled"
        RUN_BASE_SCALE["$family"]="$live_scaled"

      elif [[ "$trg" =~ ^Photon[[:space:]]+([0-9]+([.][0-9]+)?)[[:space:]]+GeV([[:space:]]*\+|,)[[:space:]]+MBD[[:space:]]+N\&?S[[:space:]]+\>\=[[:space:]]+([0-9]+)(,[[:space:]]+vtx[[:space:]]+\<[[:space:]]+([0-9]+([.][0-9]+)?)[[:space:]]+cm)?$ ]]; then
        pho_e="${BASH_REMATCH[1]}"
        base_thr="${BASH_REMATCH[4]}"
        base_vtx="${BASH_REMATCH[6]:-}"

        pho_int=$(awk -v x="$pho_e" 'BEGIN{ printf "%g", x+0 }')

        family="MBD N&S >= ${base_thr}"
        [[ -n "$base_vtx" ]] && family+=", vtx < ${base_vtx} cm"

        key="$family|$pho_int"
        RUN_PHO_BIT["$key"]="$idx"
        RUN_PHO_NAME["$key"]="$trg"
        RUN_PHO_RAW["$key"]="$raw"
        RUN_PHO_LIVE["$key"]="$live"
        RUN_PHO_SCALED["$key"]="$scaled"
        RUN_PHO_SCALE["$key"]="$live_scaled"
      fi
    done < <(sql "SELECT s.index, COALESCE(t.triggername,''), s.scaled, s.raw, s.live
                   FROM gl1_scalers s
                   LEFT JOIN gl1_triggernames t
                     ON s.index = t.index
                    AND s.runnumber BETWEEN t.runnumber AND t.runnumber_last
                  WHERE s.runnumber = $run
                  ORDER BY s.index;")

    for family in "${!RUN_BASE_BIT[@]}"; do
      FAMILY_BASE_ROWS["$family"]=$(( ${FAMILY_BASE_ROWS["$family"]:-0} + 1 ))
      (( ${RUN_BASE_RAW["$family"]:-0} > 0 )) && FAMILY_BASE_RAW_POS["$family"]=$(( ${FAMILY_BASE_RAW_POS["$family"]:-0} + 1 ))
      (( ${RUN_BASE_LIVE["$family"]:-0} > 0 )) && FAMILY_BASE_LIVE_POS["$family"]=$(( ${FAMILY_BASE_LIVE_POS["$family"]:-0} + 1 ))
      (( ${RUN_BASE_SCALED["$family"]:-0} > 0 )) && FAMILY_BASE_SCALED_POS["$family"]=$(( ${FAMILY_BASE_SCALED_POS["$family"]:-0} + 1 ))
    done

    for key in "${!RUN_PHO_BIT[@]}"; do
      family="${key%|*}"
      PHO_ROW_RUNS["$key"]=$(( ${PHO_ROW_RUNS["$key"]:-0} + 1 ))
      (( ${RUN_PHO_RAW["$key"]:-0} > 0 )) && PHO_RAW_POS["$key"]=$(( ${PHO_RAW_POS["$key"]:-0} + 1 ))
      (( ${RUN_PHO_LIVE["$key"]:-0} > 0 )) && PHO_LIVE_POS["$key"]=$(( ${PHO_LIVE_POS["$key"]:-0} + 1 ))
      (( ${RUN_PHO_SCALED["$key"]:-0} > 0 )) && PHO_SCALED_POS["$key"]=$(( ${PHO_SCALED_POS["$key"]:-0} + 1 ))

      if (( ${RUN_BASE_SCALED["$family"]:-0} > 0 && ${RUN_PHO_SCALED["$key"]:-0} > 0 )); then
        PHO_COMMON_WITH_BASE_SCALED["$key"]=$(( ${PHO_COMMON_WITH_BASE_SCALED["$key"]:-0} + 1 ))
      fi
    done

    local -a reject_reasons=()
    local status="ACCEPT"
    local reason="OK"

    if [[ -z "${RUN_BASE_BIT["$selected_family"]:-}" ]]; then
      reject_reasons+=( "missing baseline row: $selected_family" )
      ((REJ_BASE_MISSING+=1))
    elif (( ${RUN_BASE_SCALED["$selected_family"]:-0} <= 0 )); then
      reject_reasons+=( "baseline scaled<=0: bit=${RUN_BASE_BIT["$selected_family"]} raw=${RUN_BASE_RAW["$selected_family"]:-0} live=${RUN_BASE_LIVE["$selected_family"]:-0} scaled=${RUN_BASE_SCALED["$selected_family"]:-0}" )
      ((REJ_BASE_SCALED_ZERO+=1))
    fi

    if [[ -z "${RUN_PHO_BIT["$selected_family|10"]:-}" ]]; then
      reject_reasons+=( "missing Photon 10 row for $selected_family" )
      ((REJ_PHO10_MISSING+=1))
    elif (( ${RUN_PHO_SCALED["$selected_family|10"]:-0} <= 0 )); then
      reject_reasons+=( "Photon 10 scaled<=0: bit=${RUN_PHO_BIT["$selected_family|10"]} raw=${RUN_PHO_RAW["$selected_family|10"]:-0} live=${RUN_PHO_LIVE["$selected_family|10"]:-0} scaled=${RUN_PHO_SCALED["$selected_family|10"]:-0}" )
      ((REJ_PHO10_SCALED_ZERO+=1))
    fi

    if [[ -z "${RUN_PHO_BIT["$selected_family|12"]:-}" ]]; then
      reject_reasons+=( "missing Photon 12 row for $selected_family" )
      ((REJ_PHO12_MISSING+=1))
    elif (( ${RUN_PHO_SCALED["$selected_family|12"]:-0} <= 0 )); then
      reject_reasons+=( "Photon 12 scaled<=0: bit=${RUN_PHO_BIT["$selected_family|12"]} raw=${RUN_PHO_RAW["$selected_family|12"]:-0} live=${RUN_PHO_LIVE["$selected_family|12"]:-0} scaled=${RUN_PHO_SCALED["$selected_family|12"]:-0}" )
      ((REJ_PHO12_SCALED_ZERO+=1))
    fi

    if ((${#reject_reasons[@]})); then
      status="REJECT"
      reason="$(join_reason_parts "${reject_reasons[@]}")"
      ((SELECT_REJECT+=1))
      [[ -z "$FIRST_SELECT_REJECTION" ]] && FIRST_SELECT_REJECTION="run $run: $reason"
      printf "REJECT study=%s run=%s reason=\"%s\"\n" "$selected_study" "$run" "$reason" >> "$reject_file"
    else
      ((SELECT_ACCEPT+=1))
      printf "%s\n" "$run" >> "$selected_runlist"

      printf "CONFIG study=%s run=%s baselineBit=%s baselineKey=%s baselineScale=%s" \
        "$selected_study" \
        "$run" \
        "${RUN_BASE_BIT["$selected_family"]}" \
        "$selected_baseline_key" \
        "${RUN_BASE_SCALE["$selected_family"]}" >> "$txt_file"
      printf " pho10Bit=%s pho10Key=%s pho10Scale=%s" \
        "${RUN_PHO_BIT["$selected_family|10"]}" \
        "$(pho_key_for_family "$selected_family" 10)" \
        "${RUN_PHO_SCALE["$selected_family|10"]}" >> "$txt_file"
      printf " pho12Bit=%s pho12Key=%s pho12Scale=%s\n" \
        "${RUN_PHO_BIT["$selected_family|12"]}" \
        "$(pho_key_for_family "$selected_family" 12)" \
        "${RUN_PHO_SCALE["$selected_family|12"]}" >> "$txt_file"
    fi

    printf "%s," "$run" >> "$run_csv_file"
    csv_quote "$selected_study" >> "$run_csv_file"
    printf "," >> "$run_csv_file"
    csv_quote "$status" >> "$run_csv_file"
    printf "," >> "$run_csv_file"
    csv_quote "$reason" >> "$run_csv_file"

    printf ",%s," "${RUN_BASE_BIT["$selected_family"]:-}" >> "$run_csv_file"
    csv_quote "$selected_baseline_key" >> "$run_csv_file"
    printf "," >> "$run_csv_file"
    csv_quote "${RUN_BASE_NAME["$selected_family"]:-}" >> "$run_csv_file"
    printf ",%s,%s,%s,%s" "${RUN_BASE_RAW["$selected_family"]:-0}" "${RUN_BASE_LIVE["$selected_family"]:-0}" "${RUN_BASE_SCALED["$selected_family"]:-0}" "${RUN_BASE_SCALE["$selected_family"]:-0.000000}" >> "$run_csv_file"

    key="$selected_family|10"
    printf ",%s," "${RUN_PHO_BIT["$key"]:-}" >> "$run_csv_file"
    csv_quote "$(pho_key_for_family "$selected_family" 10)" >> "$run_csv_file"
    printf "," >> "$run_csv_file"
    csv_quote "${RUN_PHO_NAME["$key"]:-}" >> "$run_csv_file"
    printf ",%s,%s,%s,%s" "${RUN_PHO_RAW["$key"]:-0}" "${RUN_PHO_LIVE["$key"]:-0}" "${RUN_PHO_SCALED["$key"]:-0}" "${RUN_PHO_SCALE["$key"]:-0.000000}" >> "$run_csv_file"

    key="$selected_family|12"
    printf ",%s," "${RUN_PHO_BIT["$key"]:-}" >> "$run_csv_file"
    csv_quote "$(pho_key_for_family "$selected_family" 12)" >> "$run_csv_file"
    printf "," >> "$run_csv_file"
    csv_quote "${RUN_PHO_NAME["$key"]:-}" >> "$run_csv_file"
    printf ",%s,%s,%s,%s\n" "${RUN_PHO_RAW["$key"]:-0}" "${RUN_PHO_LIVE["$key"]:-0}" "${RUN_PHO_SCALED["$key"]:-0}" "${RUN_PHO_SCALE["$key"]:-0.000000}" >> "$run_csv_file"

    if (( processed <= 5 || processed % 100 == 0 || processed == ${#RUNS_ALL[@]} )); then
      printf "[INFO] scaledTriggerAna progress: %d / %d | run=%s\n" "$processed" "${#RUNS_ALL[@]}" "$run"
    fi
  done

  print_family_scaled_summary() {
    local family="$1"
    local title="$2"

    local base_rows="${FAMILY_BASE_ROWS["$family"]:-0}"
    local base_raw="${FAMILY_BASE_RAW_POS["$family"]:-0}"
    local base_live="${FAMILY_BASE_LIVE_POS["$family"]:-0}"
    local base_scaled="${FAMILY_BASE_SCALED_POS["$family"]:-0}"
    local any_common=0

    printf "\n%s%s%s\n" "$ANSI_BOLD$ANSI_CYAN" "$title" "$ANSI_RESET"
    printf "  Baseline family         : %s\n" "$family"
    printf "  Baseline row runs       : %d / %d\n" "$base_rows" "${#RUNS_ALL[@]}"
    printf "  Baseline raw>0 runs     : %d\n" "$base_raw"
    printf "  Baseline live>0 runs    : %d\n" "$base_live"
    printf "  Baseline scaled>0 runs  : %d\n" "$base_scaled"
    printf "\n"
    printf "%s  %8s | %8s | %8s | %8s | %9s | %11s | %-38s%s\n" \
      "$ANSI_BOLD" "Photon" "Rows" "Raw>0" "Live>0" "Scaled>0" "RowsSc0" "Decision" "$ANSI_RESET"
    printf "%s  %-8s-+-%-8s-+-%-8s-+-%-8s-+-%-9s-+-%-11s-+-%-38s%s\n" \
      "$ANSI_DIM" "--------" "--------" "--------" "--------" "---------" "-----------" "--------------------------------------" "$ANSI_RESET"

    local -a phos=()
    local k pho
    mapfile -t phos < <(
      for k in "${!PHO_ROW_RUNS[@]}"; do
        [[ "$k" == "$family|"* ]] || continue
        printf "%s\n" "${k#"$family|"}"
      done | sort -n
    )

    if ((${#phos[@]} == 0)); then
      printf "  %-8s | %8s | %8s | %8s | %9s | %11s | %-38s\n" "-" "-" "-" "-" "-" "-" "no photon rows found"
      return
    fi

    for pho in "${phos[@]}"; do
      k="$family|$pho"

      local rows="${PHO_ROW_RUNS["$k"]:-0}"
      local rawpos="${PHO_RAW_POS["$k"]:-0}"
      local livepos="${PHO_LIVE_POS["$k"]:-0}"
      local scaledpos="${PHO_SCALED_POS["$k"]:-0}"
      local common="${PHO_COMMON_WITH_BASE_SCALED["$k"]:-0}"
      local rows_scaled_zero=$(( rows - scaledpos ))
      local decision=""

      if (( common > 0 )); then
        any_common=1
      fi

      if [[ "$family" == "$selected_family" && ( "$pho" == "10" || "$pho" == "12" ) ]]; then
        if (( common > 0 )); then
          decision="USE in ${selected_study}"
        else
          decision="wanted, but no common scaled sample"
        fi
      elif (( scaledpos == 0 )); then
        decision="do not use: scaled=0 in all row runs"
      elif [[ "$family" == "$selected_family" ]]; then
        decision="available but not selected for first pass"
      else
        decision="diagnostic only for this pass"
      fi

      printf "  %8s | %8d | %8d | %8d | %9d | %11d | %-38s\n" \
        "$pho" "$rows" "$rawpos" "$livepos" "$scaledpos" "$rows_scaled_zero" "$decision"

      csv_quote "$family" >> "$summary_csv_file"
      printf "," >> "$summary_csv_file"
      csv_quote "$pho" >> "$summary_csv_file"
      printf ",%d,%d,%d,%d,%d,%d," "$rows" "$rawpos" "$livepos" "$scaledpos" "$rows_scaled_zero" "$common" >> "$summary_csv_file"
      csv_quote "$decision" >> "$summary_csv_file"
      printf "\n" >> "$summary_csv_file"
    done

    printf "\n"
    printf "  Common scaled-positive runs with baseline:\n"
    for pho in "${phos[@]}"; do
      k="$family|$pho"
      printf "    Photon %-2s : %d\n" "$pho" "${PHO_COMMON_WITH_BASE_SCALED["$k"]:-0}"
    done

    if [[ "$family" == "$selected_family" ]]; then
      printf "  Common scaled-positive subset for first-pass overlay:\n"
      printf "    %s : %d runs\n" "$selected_study" "$SELECT_ACCEPT"
      printf "    requirement = baseline scaled>0 + Photon 10 scaled>0 + Photon 12 scaled>0 in the same run\n"
    elif (( any_common == 0 )); then
      printf "  Conclusion:\n"
      printf "    no photon trigger has any scaled-positive overlap with this baseline, so this subset is not used for the first scaled-trigger study\n"
    fi
  }

  print_family_scaled_summary "MBD N&S >= 2, vtx < 150 cm" "Available photon triggers for MBD N&S >= 2, vtx < 150 cm"
  print_family_scaled_summary "MBD N&S >= 2, vtx < 10 cm" "Diagnostic photon-trigger summary for MBD N&S >= 2, vtx < 10 cm"

  printf "\n%sSelected scaled-efficiency run list%s\n" "$ANSI_BOLD$ANSI_GREEN" "$ANSI_RESET"
  printf "  Selected study : %s\n" "$selected_study"
  printf "  Rule           : MBD vtx<150 scaled>0 AND Photon10 vtx<150 scaled>0 AND Photon12 vtx<150 scaled>0\n"
  printf "  Accepted runs  : %d\n" "$SELECT_ACCEPT"
  printf "  Rejected runs  : %d\n" "$SELECT_REJECT"
  printf "  Run list       : %s\n" "$selected_runlist"
  printf "  CONFIG table   : %s\n" "$txt_file"
  printf "  Run diagnostics: %s\n" "$run_csv_file"
  printf "  Summary CSV    : %s\n" "$summary_csv_file"
  printf "  Rejections     : %s\n" "$reject_file"

  printf "\n%sSelected-study rejection counters%s\n" "$ANSI_BOLD" "$ANSI_RESET"
  printf "  missing baseline      : %d\n" "$REJ_BASE_MISSING"
  printf "  baseline scaled<=0    : %d\n" "$REJ_BASE_SCALED_ZERO"
  printf "  missing Photon 10     : %d\n" "$REJ_PHO10_MISSING"
  printf "  Photon 10 scaled<=0   : %d\n" "$REJ_PHO10_SCALED_ZERO"
  printf "  missing Photon 12     : %d\n" "$REJ_PHO12_MISSING"
  printf "  Photon 12 scaled<=0   : %d\n" "$REJ_PHO12_SCALED_ZERO"

  if (( SELECT_REJECT > 0 )); then
    printf "%s  First rejection       : %s%s\n" "$ANSI_YELLOW" "${FIRST_SELECT_REJECTION:-none recorded}" "$ANSI_RESET"
  fi

  printf "\n%sInterpretation%s\n" "$ANSI_BOLD$ANSI_CYAN" "$ANSI_RESET"
  echo "  - Photon 6/8 can appear in menu/raw/live bookkeeping, but they are not valid for a scaled-bit correction if scaled>0 is zero."
  echo "  - The first formal scaled-efficiency check should therefore use the common MBD vtx<150 + Photon 10 + Photon 12 scaled-positive run set."
  echo "  - The vtx<10 table is printed as a diagnostic to show why that subset is not chosen for the first scaled-trigger study."
  echo

  return 0
}

run_trigger_quickcheck() {
  command -v psql >/dev/null 2>&1 || { echo "[ERROR] psql not found in PATH"; exit 1; }

  if [[ "${LABEL:-}" != "auau" ]]; then
    echo "[ERROR] doQUICKcheck is currently defined only for auau mode."
    echo "        Run: ./make_dstListsData.sh auau QA doQUICKcheck"
    return 1
  fi

  PSQL=(psql -h sphnxdaqdbreplica -d daq -At -F $'\t' -q)
  sql()         { "${PSQL[@]}" -c "$1" 2>/dev/null || true; }
  num_or_zero() { [[ $1 =~ ^-?[0-9]+([.][0-9]+)?$ ]] && printf '%s' "$1" || printf 0; }

  local GOOD_RUN=78280
  local -a BAD_RUNS=(68414 68491 71578 75299)

  local ANSI_BOLD ANSI_DIM ANSI_CYAN ANSI_GREEN ANSI_YELLOW ANSI_RESET
  ANSI_BOLD=$'\033[1m'
  ANSI_DIM=$'\033[2m'
  ANSI_CYAN=$'\033[36m'
  ANSI_GREEN=$'\033[32m'
  ANSI_YELLOW=$'\033[33m'
  ANSI_RESET=$'\033[0m'

  echo
  printf "%sAuAu trigger live/raw quick check%s\n" "$ANSI_BOLD$ANSI_CYAN" "$ANSI_RESET"
  echo "Good reference run: $GOOD_RUN"
  echo "Bad-reference runs: ${BAD_RUNS[*]}"
  echo "Focus: MBD/photon trigger rows from gl1_scalers joined to gl1_triggernames"
  echo "Columns:"
  echo "  L/R       = live/raw for that run and trigger"
  echo "  Ref L/R   = live/raw for the same trigger name in reference run $GOOD_RUN"
  echo "  Δ(L/R)    = L/R - Ref L/R"
  echo "  Live/Scal = live/scaled, effective prescale reference when scaled > 0"
  echo

  declare -A REF_LR
  declare -A REF_RAW
  declare -A REF_LIVE
  declare -A REF_SCALED
  declare -A REF_BIT

  while IFS=$'\t' read -r idx trg raw live scaled; do
    [[ -z "$idx" || -z "$trg" ]] && continue

    idx=$(num_or_zero "$idx")
    raw=$(num_or_zero "$raw")
    live=$(num_or_zero "$live")
    scaled=$(num_or_zero "$scaled")

    [[ "$idx" =~ ^[0-9]+$ ]] || continue
    [[ -n "$trg" ]] || continue

    lr=$(awk -v l="$live" -v r="$raw" 'BEGIN{ if (r>0) printf "%.6f", l/r; else printf "nan"; }')

    REF_LR["$trg"]="$lr"
    REF_RAW["$trg"]="$raw"
    REF_LIVE["$trg"]="$live"
    REF_SCALED["$trg"]="$scaled"
    REF_BIT["$trg"]="$idx"
  done < <(sql "SELECT s.index, COALESCE(t.triggername,''), s.raw, s.live, s.scaled
                 FROM gl1_scalers s
                 LEFT JOIN gl1_triggernames t
                   ON s.index = t.index
                  AND s.runnumber BETWEEN t.runnumber AND t.runnumber_last
                WHERE s.runnumber = ${GOOD_RUN}
                  AND (t.triggername ILIKE '%MBD%' OR t.triggername ILIKE '%Photon%')
                ORDER BY s.index;")

  if ((${#REF_LR[@]} == 0)); then
    echo "[ERROR] Reference run $GOOD_RUN returned no MBD/photon trigger rows."
    return 1
  fi

  printf "%s  %8s | %-5s | %-46s | %14s | %14s | %9s | %9s | %9s | %9s%s\n" \
    "$ANSI_BOLD" "Run" "Bit" "TriggerName" "Raw" "Live" "L/R" "Ref L/R" "Δ(L/R)" "Live/Scal" "$ANSI_RESET"
  printf "%s  %-8s-+-%-5s-+-%-46s-+-%-14s-+-%-14s-+-%-9s-+-%-9s-+-%-9s-+-%-9s%s\n" \
    "$ANSI_DIM" \
    "--------" "-----" "----------------------------------------------" "--------------" "--------------" "---------" "---------" "---------" "---------" \
    "$ANSI_RESET"

  print_run_rows() {
    local run="$1"
    local tag="$2"

    while IFS=$'\t' read -r idx trg raw live scaled; do
      [[ -z "$idx" || -z "$trg" ]] && continue

      idx=$(num_or_zero "$idx")
      raw=$(num_or_zero "$raw")
      live=$(num_or_zero "$live")
      scaled=$(num_or_zero "$scaled")

      [[ "$idx" =~ ^[0-9]+$ ]] || continue

      lr=$(awk -v l="$live" -v r="$raw" 'BEGIN{ if (r>0) printf "%.6f", l/r; else printf "nan"; }')
      live_scaled=$(awk -v l="$live" -v s="$scaled" 'BEGIN{ if (s>0) printf "%.3f", l/s; else printf "0.000"; }')

      ref_lr="${REF_LR["$trg"]:-nan}"
      if [[ "$ref_lr" == "nan" || "$lr" == "nan" ]]; then
        delta="nan"
      else
        delta=$(awk -v a="$lr" -v b="$ref_lr" 'BEGIN{ printf "%.6f", a-b; }')
      fi

      printf "  %8s | %5d | %-46.46s | %14d | %14d | %9s | %9s | %9s | %9s\n" \
        "$run" "$idx" "$trg" "$raw" "$live" "$lr" "$ref_lr" "$delta" "$live_scaled"
    done < <(sql "SELECT s.index, COALESCE(t.triggername,''), s.raw, s.live, s.scaled
                   FROM gl1_scalers s
                   LEFT JOIN gl1_triggernames t
                     ON s.index = t.index
                    AND s.runnumber BETWEEN t.runnumber AND t.runnumber_last
                  WHERE s.runnumber = ${run}
                    AND (t.triggername ILIKE '%MBD%' OR t.triggername ILIKE '%Photon%')
                  ORDER BY s.index;")

    if [[ "$tag" == "good" ]]; then
      printf "%s  %-8s-+-%-5s-+-%-46s-+-%-14s-+-%-14s-+-%-9s-+-%-9s-+-%-9s-+-%-9s%s\n" \
        "$ANSI_DIM" \
        "--------" "-----" "----------------------------------------------" "--------------" "--------------" "---------" "---------" "---------" "---------" \
        "$ANSI_RESET"
    fi
  }

  print_run_rows "$GOOD_RUN" "good"
  for run in "${BAD_RUNS[@]}"; do
    print_run_rows "$run" "bad"
  done

  echo
  printf "%sInterpretation guide:%s\n" "$ANSI_BOLD$ANSI_GREEN" "$ANSI_RESET"
  echo "  - Compare bad-run L/R against the same trigger-name Ref L/R from run $GOOD_RUN."
  echo "  - Large negative Δ(L/R) means that trigger has lower live/raw fraction than the good reference."
  echo "  - Live/Scal is the effective run-level prescale correction for scaled-trigger histograms."
  echo
}

build_pp24_lists() {
  command -v CreateDstList.pl >/dev/null 2>&1 || { echo "[ERROR] CreateDstList.pl not in PATH"; exit 1; }
  [[ -f "$LIST_FILE" ]] || { echo "[ERROR] Run list not found: $LIST_FILE"; exit 1; }
  [[ -s "$LIST_FILE" ]] || { echo "[ERROR] Run list is empty: $LIST_FILE"; exit 1; }

  CREATE_DST_TOOL="$(command -v CreateDstList.pl || true)"
  [[ -n "$CREATE_DST_TOOL" ]] || { echo "[ERROR] CreateDstList.pl not in PATH"; exit 1; }

  mkdir -p "$OUT_DIR"
  cd "$OUT_DIR" || exit 1

  echo "Mode: $LABEL"
  echo "Tag: $TAG"
  echo "Type: $TYPE"
  echo "Input: $LIST_FILE"
  echo "Output dir: $OUT_DIR"
  echo "Tool: $CREATE_DST_TOOL"
  echo "----------------------------------------"
  echo "[INFO] Reading run list for $LABEL"

  read_runs_from_file "$LIST_FILE" RUNS_ALL
  echo "[INFO] Valid runs found: ${#RUNS_ALL[@]}"
  echo "[INFO] Removing old dst_*.list and DST_*.list files from $OUT_DIR"

  total=0
  made=0

  rm -f dst_*.list DST_*.list

  for runnum in "${RUNS_ALL[@]}"; do
    ((total++))
    pad=$(printf "%08d" "$runnum")

    echo "[INFO] ($total / ${#RUNS_ALL[@]}) Creating $TYPE list for run $runnum"
    perl "$CREATE_DST_TOOL" --tag "$TAG" --run "$runnum" "$TYPE"
    rc=$?

    if [[ -s "dst_jetcalo-$pad.list" || -s "DST_JETCALO-$pad.list" ]]; then
      ((made++))
      echo "[INFO] Created list for run $runnum"
    else
      echo "[WARN] No list file found (or file empty) for run $runnum (exit=$rc)"
    fi
  done

  echo "----------------------------------------"
  echo "Requested runs: $total"
  echo "List files made: $made"
  echo "Done. Files are in: $OUT_DIR"
  print_naming_scheme_summary
}

build_dataset_lists() {
  command -v CreateDstList.pl >/dev/null 2>&1 || { echo "[ERROR] CreateDstList.pl not in PATH"; exit 1; }
  [[ -f "$LIST_FILE" ]] || { echo "[ERROR] Run list not found: $LIST_FILE"; exit 1; }
  [[ -s "$LIST_FILE" ]] || { echo "[ERROR] Run list is empty: $LIST_FILE"; exit 1; }

  echo "[INFO] Reading run list for $LABEL from $LIST_FILE"
  read_runs_from_file "$LIST_FILE" RUNS_ALL
  (( ${#RUNS_ALL[@]} > 0 )) || { echo "[ERROR] No valid run numbers found in: $LIST_FILE"; exit 1; }
  echo "[INFO] Valid runs found: ${#RUNS_ALL[@]}"

  echo "[INFO] Ensuring output directory exists: $OUT_DIR"
  mkdir -p "$OUT_DIR"
  case "$OUT_DIR" in
    /sphenix/u/patsfan753/scratch/thesisAnalysis/dst_lists_pp_run25|\
    /sphenix/u/patsfan753/scratch/thesisAnalysis/dst_lists_auau|\
    /sphenix/u/patsfan753/scratch/thesisAnalysis/dst_lists_oo|\
    /sphenix/u/patsfan753/scratch/thesisAnalysis/dst_lists_auau_run2)
      echo "[INFO] Clearing existing list files in $OUT_DIR"
      ( shopt -s dotglob nullglob; rm -rf "${OUT_DIR}/"* )
      ;;
    *)
      echo "[ERROR] Refusing to purge unknown OUT_DIR: $OUT_DIR"
      exit 1
      ;;
  esac

  cd "$OUT_DIR" || exit 1

  echo "Mode: $LABEL"
  echo "Dataset: $DATASET"
  echo "Tag: $TAG"
  echo "Prefix: $PREFIX"
  echo "Input: $LIST_FILE"
  echo "Output dir: $OUT_DIR"
  echo "----------------------------------------"
  echo "[INFO] Running CreateDstList.pl for dataset=$DATASET prefix=$PREFIX tag=$TAG"

  CreateDstList.pl --tag "$TAG" --dataset "$DATASET" --list "$LIST_FILE" "$PREFIX"
  rc=$?
  echo "[INFO] CreateDstList.pl finished with exit code: $rc"

  pair_rc=0
  if ! pair_zdc_raw_for_auau_lists; then
    pair_rc=1
  fi

  echo "[INFO] Scanning output files for expected per-run list coverage"

  local stem
  stem="${PREFIX#DST_}"
  stem="$(echo "$stem" | tr '[:upper:]' '[:lower:]')"

  made=0
  for r in "${RUNS_ALL[@]}"; do
    r8=$(printf "%08d" "$((10#$r))")
    if [[ -s "${OUT_DIR}/dst_${stem}-${r8}.list" || -s "${OUT_DIR}/${PREFIX}-${r8}.list" || -s "${OUT_DIR}/${PREFIX}_${DATASET}_${TAG}-${r8}.list" ]]; then
      ((made+=1))
    else
      echo "[WARN] Missing expected list for run $r"
    fi
  done
  miss_count=$(( ${#RUNS_ALL[@]} - made ))

  echo "----------------------------------------"
  echo "Build output coverage summary for $LABEL"
  echo "Requested runs: ${#RUNS_ALL[@]}"
  echo "Submit-ready canonical list files made: $made"
  echo "Runs without submit-ready canonical list: $miss_count"
  echo "CreateDstList exit code: $rc"

  if [[ "$LABEL" == "auau" || "$LABEL" == "run2auau" ]]; then
    local report_file="${OUT_DIR}/dst_pairing_report_${LABEL}.txt"
    if [[ -s "$report_file" ]]; then
      echo "----------------------------------------"
      echo "DST pairing diagnostic summary from: $report_file"
      awk -v prefix="${PREFIX}" '
        /^Run-level summary:/ {keep=1}
        keep && $0 == "Runs with " prefix " list:" {exit}
        keep {print}
      ' "$report_file" || true
    else
      echo "[WARN] No DST pairing report found at: $report_file"
    fi

    if (( pair_rc != 0 )); then
      echo "DST pairing status: ERROR_SEGMENT_PAIRING_OR_NO_PAIRS"
    else
      echo "DST pairing status: OK_FOR_AVAILABLE_OVERLAP"
    fi
  elif (( pair_rc != 0 )); then
    echo "DST pairing status: ERROR"
  fi

  echo "Done. Files are in: $OUT_DIR"
  print_naming_scheme_summary

  if (( rc != 0 || pair_rc != 0 )); then
    return 1
  fi

  if [[ "$LABEL" != "auau" && "$LABEL" != "run2auau" && "$miss_count" -gt 0 ]]; then
    return 1
  fi

  return 0
}

if [[ -n "$ACTION" && "$ACTION" != "QA" ]]; then
  usage
  exit 1
fi

if [[ -n "$EXTRA_ACTION" && "$EXTRA_ACTION" != "doQUICKcheck" && "$EXTRA_ACTION" != "scaledTriggerAna" ]]; then
  usage
  exit 1
fi

if [[ ( "$EXTRA_ACTION" == "doQUICKcheck" || "$EXTRA_ACTION" == "scaledTriggerAna" ) && ( "$MODE" != "auau" || "$ACTION" != "QA" ) ]]; then
  usage
  exit 1
fi

if [[ "$MODE" == "ALL" ]]; then
  run_all_modes
  exit $?
fi

setup_mode

if [[ "$ACTION" == "QA" ]]; then
  if [[ "$EXTRA_ACTION" == "doQUICKcheck" ]]; then
    run_trigger_quickcheck
  elif [[ "$EXTRA_ACTION" == "scaledTriggerAna" ]]; then
    run_scaled_trigger_ana
  else
    run_trigger_qa
  fi
  exit 0
fi

if [[ "$BUILD_MODE" == "per_run" ]]; then
  build_pp24_lists
  exit 0
fi

build_dataset_lists
