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

usage() {
  cat <<'USAGE'
Usage:
  ./make_dstListsData.sh pp24 [QA]
  ./make_dstListsData.sh pp25 [QA]
  ./make_dstListsData.sh auau [QA]
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
      LIST_FILE="$GRL_BASE/run3auau_new_newcdbtag_v008_dst_calofitting_grl.list"
      OUT_DIR="$IN_BASE/dst_lists_auau"
      TAG="new_newcdbtag_v008"
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

      if (( scaled > 0 && raw >= 0 )); then
        sdf=$(awk -v r="$raw" -v s="$scaled" 'BEGIN{ if (s>0) printf "%.6f", r/s; else printf "0"; }')
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
  printf "  Rows used for average scaledown calculation (reference only): %d\n" "$rows_used_for_avg"
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
      "$ANSI_BOLD" "Bit" "$maxlen" "TriggerName" "MenuRuns" "RawRuns" "LiveRuns" "ScaledRuns" "SumRaw" "SumLive" "L/R" "AvgScale" "$ANSI_RESET"
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
      local mbd_bit mbd_name mbd_thr mbd_vtx pho_bit pho_name pho_thr pho_vtx pho_e
      local -a novtx_rows=()
      local -a vtx_rows=()

      mkdir -p "$OUT_DIR"

      for mbd_bit in "${!TRIG_NAME[@]}"; do
        mbd_name="${TRIG_NAME["$mbd_bit"]:-}"

        if [[ "$mbd_name" =~ ^MBD[[:space:]]+N\&S[[:space:]]+\>\=[[:space:]]+([0-9]+)(,[[:space:]]+vtx[[:space:]]+\<[[:space:]]+([0-9]+)[[:space:]]+cm)?$ ]]; then
          mbd_thr="${BASH_REMATCH[1]}"
          mbd_vtx="${BASH_REMATCH[3]:-}"

          for pho_bit in "${!TRIG_NAME[@]}"; do
            pho_name="${TRIG_NAME["$pho_bit"]:-}"

            if [[ "$pho_name" =~ ^Photon[[:space:]]+([0-9]+)[[:space:]]+GeV([[:space:]]+\+[[:space:]]+MBD[[:space:]]+NS|,[[:space:]]+MBD[[:space:]]+N\&S)[[:space:]]+\>\=[[:space:]]+([0-9]+)(,[[:space:]]+vtx[[:space:]]+\<[[:space:]]+([0-9]+)[[:space:]]+cm)?$ ]]; then
              pho_e="${BASH_REMATCH[1]}"
              pho_thr="${BASH_REMATCH[3]}"
              pho_vtx="${BASH_REMATCH[5]:-}"

              if [[ "$pho_thr" == "$mbd_thr" ]]; then
                if [[ -z "$mbd_vtx" && -z "$pho_vtx" ]]; then
                  novtx_rows+=( "$mbd_thr|$pho_e|$mbd_bit|$mbd_name|$pho_bit|$pho_name" )
                elif [[ -n "$mbd_vtx" && "$pho_vtx" == "$mbd_vtx" ]]; then
                  vtx_rows+=( "$mbd_thr|$mbd_vtx|$pho_e|$mbd_bit|$mbd_name|$pho_bit|$pho_name" )
                fi
              fi
            fi
          done
        fi
      done

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
        echo "No vertex cut combinations"
        if ((${#novtx_rows[@]})); then
          printf "%s\n" "${novtx_rows[@]}" | sort -t'|' -k1,1n -k2,2n | while IFS='|' read -r _thr _e mbd_bit mbd_name pho_bit pho_name; do
            printf "  MBD bit %-3s %-35s || Photon bit %-3s %s\n" "$mbd_bit" "$mbd_name" "$pho_bit" "$pho_name"
          done
        else
          echo "  (none)"
        fi
        echo
        echo "Vertex cut combinations"
        if ((${#vtx_rows[@]})); then
          printf "%s\n" "${vtx_rows[@]}" | sort -t'|' -k1,1n -k2,2n -k3,3n | while IFS='|' read -r _thr _vtx _e mbd_bit mbd_name pho_bit pho_name; do
            printf "  MBD bit %-3s %-35s || Photon bit %-3s %s\n" "$mbd_bit" "$mbd_name" "$pho_bit" "$pho_name"
          done
        else
          echo "  (none)"
        fi
      } | tee "$outfile"

      echo "[INFO] Trigger grouping summary saved to: $outfile"
      echo
    }

    printf "%sTrigger-family focus for %s%s\n" "$ANSI_BOLD$ANSI_CYAN" "$LABEL" "$ANSI_RESET"
    echo
    print_focus_table "No vertex cut focus table" 10 24 25 26 27
    print_focus_table "Vertex cut focus table" 12 36 37 38
    print_trigger_groupings
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
  echo "Requested runs: ${#RUNS_ALL[@]}"
  echo "List files made: $made"
  if (( miss_count > 0 )); then
    echo "Runs without a list: $miss_count"
    echo "CreateDstList exit code: $rc"
  fi
  echo "Done. Files are in: $OUT_DIR"
  print_naming_scheme_summary

  if (( rc != 0 || miss_count > 0 )); then
    return 1
  fi

  return 0
}

if [[ -n "$ACTION" && "$ACTION" != "QA" ]]; then
  usage
  exit 1
fi

if [[ "$MODE" == "ALL" ]]; then
  run_all_modes
  exit $?
fi

setup_mode

if [[ "$ACTION" == "QA" ]]; then
  run_trigger_qa
  exit 0
fi

if [[ "$BUILD_MODE" == "per_run" ]]; then
  build_pp24_lists
  exit 0
fi

build_dataset_lists
