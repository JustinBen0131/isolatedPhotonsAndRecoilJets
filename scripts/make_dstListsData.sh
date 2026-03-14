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
  local -a modes=(pp24 pp25 auau oo25 run2auau)
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
  echo "  - trigger name"
  echo "  - trigger bit"
  echo "  - active runs (scaledown != -1)"
  echo "  - average scaledown factor over active runs"
  echo "  - total scaled events"
  echo "Ranking: most scaled events to least scaled events"
  echo "QA mode behavior:"
  echo "  - reads only the input golden run list"
  echo "  - queries trigger metadata from daq tables"
  echo "  - does not create, clean, scan, or modify any DST list files"
  echo "Active-run rule:"
  echo "  - a trigger is counted active for a run only when scaledown != -1"
  echo "Average scaledown rule:"
  echo "  - computed only from rows with scaled > 0 and raw >= 0"
  echo "Scaled-event rule:"
  echo "  - accumulated only from rows with scaled > 0"
  echo

  echo "[INFO] Verifying psql connectivity expectation: host=sphnxdaqdbreplica db=daq"
  echo "[INFO] Beginning trigger scan over ${#RUNS_ALL[@]} golden-list runs"

  declare -A TRIG_NAME
  declare -A TRIG_ACTIVE_RUNS
  declare -A TRIG_SUMSCALED
  declare -A TRIG_SUMSDFACTOR
  declare -A TRIG_SDCOUNT

  local total_query_rows=0
  local runs_with_no_rows=0
  local runs_with_only_inactive=0
  local runs_with_active_rows=0
  local rows_scaled_minus_one=0
  local rows_scaled_zero=0
  local rows_scaled_positive=0
  local rows_bad_idx=0
  local rows_missing_name=0
  local rows_used_for_avg=0
  local first_problem_run=""

  gidx=0
  for run in "${RUNS_ALL[@]}"; do
    ((gidx+=1))
    declare -A SEEN_ACTIVE=()

    local run_row_count=0
    local run_active_row_count=0
    local run_positive_scaled_count=0
    local run_avg_usable_count=0

    while IFS=$'\t' read -r trg idx scaled raw; do
      [[ -z "$idx" ]] && continue
      ((run_row_count+=1))
      ((total_query_rows+=1))

      idx=$(num_or_zero "$idx")
      scaled=$(num_or_zero "$scaled")
      raw=$(num_or_zero "$raw")

      if [[ ! "$idx" =~ ^[0-9]+$ ]]; then
        ((rows_bad_idx+=1))
        [[ -z "$first_problem_run" ]] && first_problem_run="$run"
        continue
      fi

      [[ -z "$trg" ]] && ((rows_missing_name+=1))

      TRIG_NAME["$idx"]="$trg"

      if (( scaled == -1 )); then
        ((rows_scaled_minus_one+=1))
      elif (( scaled == 0 )); then
        ((rows_scaled_zero+=1))
      elif (( scaled > 0 )); then
        ((rows_scaled_positive+=1))
      fi

      if (( scaled != -1 )); then
        ((run_active_row_count+=1))

        if [[ -z "${SEEN_ACTIVE[$idx]:-}" ]]; then
          TRIG_ACTIVE_RUNS["$idx"]=$(( ${TRIG_ACTIVE_RUNS["$idx"]:-0} + 1 ))
          SEEN_ACTIVE["$idx"]=1
        fi

        if (( scaled > 0 )); then
          TRIG_SUMSCALED["$idx"]=$(( ${TRIG_SUMSCALED["$idx"]:-0} + scaled ))
          ((run_positive_scaled_count+=1))
        fi

        if (( scaled > 0 && raw >= 0 )); then
          sdf=$(awk -v r="$raw" -v s="$scaled" 'BEGIN{ if (s>0) printf "%.6f", r/s; else printf "0"; }')
          TRIG_SUMSDFACTOR["$idx"]="$(awk -v a="${TRIG_SUMSDFACTOR["$idx"]:-0}" -v b="$sdf" 'BEGIN{ printf "%.6f", a+b }')"
          TRIG_SDCOUNT["$idx"]=$(( ${TRIG_SDCOUNT["$idx"]:-0} + 1 ))
          ((run_avg_usable_count+=1))
          ((rows_used_for_avg+=1))
        fi
      fi
    done < <(sql "SELECT t.triggername, s.index, s.scaled, s.raw
                   FROM gl1_scalers s
                   JOIN gl1_triggernames t
                     ON s.index = t.index
                    AND s.runnumber BETWEEN t.runnumber AND t.runnumber_last
                  WHERE s.runnumber = $run;")

    if (( run_row_count == 0 )); then
      ((runs_with_no_rows+=1))
      [[ -z "$first_problem_run" ]] && first_problem_run="$run"
      printf "%s[WARN]%s Run %s returned zero rows from gl1_scalers/gl1_triggernames join\n" "$ANSI_YELLOW" "$ANSI_RESET" "$run"
    elif (( run_active_row_count == 0 )); then
      ((runs_with_only_inactive+=1))
      printf "%s[WARN]%s Run %s returned %d trigger rows, but all had scaledown = -1\n" \
        "$ANSI_YELLOW" "$ANSI_RESET" "$run" "$run_row_count"
    else
      ((runs_with_active_rows+=1))
    fi

    if (( gidx <= 5 || gidx % 100 == 0 || gidx == ${#RUNS_ALL[@]} )); then
      printf "[INFO] Trigger QA progress: %d / %d | run=%s | rows=%d | activeRows=%d | scaled>0 rows=%d | avgUsableRows=%d\n" \
        "$gidx" "${#RUNS_ALL[@]}" "$run" "$run_row_count" "$run_active_row_count" "$run_positive_scaled_count" "$run_avg_usable_count"
    fi
  done

  tmp_rows=()
  maxlen=11
  for idx in "${!TRIG_ACTIVE_RUNS[@]}"; do
    trg="${TRIG_NAME["$idx"]}"
    [[ -z "$trg" ]] && trg="(missing-name)"
    (( ${#trg} > maxlen )) && maxlen=${#trg}
    active_runs=${TRIG_ACTIVE_RUNS["$idx"]:-0}
    sumscaled=${TRIG_SUMSCALED["$idx"]:-0}
    sdcount=${TRIG_SDCOUNT["$idx"]:-0}
    if (( sdcount > 0 )); then
      avg_sd=$(awk -v a="${TRIG_SUMSDFACTOR["$idx"]:-0}" -v n="$sdcount" 'BEGIN{ if (n>0) printf "%.3f", a/n; else printf "0.000"; }')
    else
      avg_sd="0.000"
    fi
    tmp_rows+=( "$idx|$trg|$active_runs|$avg_sd|$sumscaled|$sdcount" )
  done
  (( maxlen+=2 ))

  echo
  printf "%sQA bookkeeping diagnostics%s\n" "$ANSI_BOLD$ANSI_CYAN" "$ANSI_RESET"
  printf "  Total query rows read: %d\n" "$total_query_rows"
  printf "  Runs with active trigger rows: %d / %d\n" "$runs_with_active_rows" "${#RUNS_ALL[@]}"
  printf "  Runs with zero returned rows: %d\n" "$runs_with_no_rows"
  printf "  Runs with only scaledown=-1 rows: %d\n" "$runs_with_only_inactive"
  printf "  Row counts by scaled value class:\n"
  printf "    scaled = -1 : %d\n" "$rows_scaled_minus_one"
  printf "    scaled = 0  : %d\n" "$rows_scaled_zero"
  printf "    scaled > 0  : %d\n" "$rows_scaled_positive"
  printf "  Rows used for average scaledown calculation: %d\n" "$rows_used_for_avg"
  printf "  Rows with missing trigger names: %d\n" "$rows_missing_name"
  printf "  Rows with unusable trigger bit values: %d\n" "$rows_bad_idx"
  if [[ -n "$first_problem_run" ]]; then
    printf "%s  First run showing a warning condition: %s%s\n" "$ANSI_YELLOW" "$first_problem_run" "$ANSI_RESET"
  else
    printf "%s  No per-run warning conditions were encountered during scan%s\n" "$ANSI_GREEN" "$ANSI_RESET"
  fi

  if ((${#tmp_rows[@]})); then
    echo
    printf "%sDataset: %s%s\n" "$ANSI_BOLD$ANSI_GREEN" "$LABEL" "$ANSI_RESET"
    printf "%sTotal available active triggers in golden-run scan: %d%s\n" "$ANSI_DIM" "${#tmp_rows[@]}" "$ANSI_RESET"
    printf "%sTable notes: ActiveRuns uses scaledown != -1, AvgScale uses only scaled > 0 rows, ScaledEvents is the total of scaled > 0 rows.%s\n" "$ANSI_DIM" "$ANSI_RESET"
    echo

    printf "%s  %4s | %-*s | %11s | %11s | %14s | %10s%s\n" \
      "$ANSI_BOLD" "Bit" "$maxlen" "TriggerName" "ActiveRuns" "AvgScale" "ScaledEvents" "AvgNRows" "$ANSI_RESET"
    printf "%s  %-4s-+-%-*s-+-%-11s-+-%-11s-+-%-14s-+-%-10s%s\n" \
      "$ANSI_DIM" \
      "$(printf '─%.0s' $(seq 1 4))" \
      "$maxlen" "$(printf '─%.0s' $(seq 1 $maxlen))" \
      "$(printf '─%.0s' $(seq 1 11))" \
      "$(printf '─%.0s' $(seq 1 11))" \
      "$(printf '─%.0s' $(seq 1 14))" \
      "$(printf '─%.0s' $(seq 1 10))" \
      "$ANSI_RESET"

    printf "%s\n" "${tmp_rows[@]}" \
    | awk -F'|' '{printf "%015d|%04d|%011d|%s|%s|%010d\n",$5,$1,$3,$4,$2,$6}' \
    | sort -r -n \
    | while IFS='|' read -r padScaled padIdx padActive avgsd trg padAvgN; do
        sumscaled=$((10#$padScaled))
        idx=$((10#$padIdx))
        active_runs=$((10#$padActive))
        avg_n=$((10#$padAvgN))
        printf "  %4d | %-*s | %11d | %11s | %14d | %10d\n" \
          "$idx" "$maxlen" "$trg" "$active_runs" "$avgsd" "$sumscaled" "$avg_n"
      done
    echo

    local special_rows=()
    local special_maxlen=11
    local special_bit special_name special_active special_avg special_scaled special_avgn

    add_special_row() {
      local want_idx="$1"
      local have_row=""
      local row
      for row in "${tmp_rows[@]}"; do
        IFS='|' read -r special_bit special_name special_active special_avg special_scaled special_avgn <<< "$row"
        if [[ "$special_bit" == "$want_idx" ]]; then
          have_row="$row"
          break
        fi
      done

      if [[ -n "$have_row" ]]; then
        IFS='|' read -r special_bit special_name special_active special_avg special_scaled special_avgn <<< "$have_row"
      else
        special_bit="$want_idx"
        case "$want_idx" in
          10) special_name="MBD N&S >= 1" ;;
          12) special_name="MBD N&S >= 1, vtx < 10 cm" ;;
          24) special_name="Photon 2 GeV+ MBD NS >= 1" ;;
          25) special_name="Photon 3 GeV + MBD NS >= 1" ;;
          26) special_name="Photon 4 GeV + MBD NS >= 1" ;;
          27) special_name="Photon 5 GeV + MBD NS >= 1" ;;
          36) special_name="Photon 3 GeV, MBD N&S >= 1, vtx < 10 cm" ;;
          37) special_name="Photon 4 GeV, MBD N&S >= 1, vtx < 10 cm" ;;
          38) special_name="Photon 5 GeV, MBD N&S >= 1, vtx < 10 cm" ;;
          *)  special_name="(missing-name)" ;;
        esac
        special_active=0
        special_avg="0.000"
        special_scaled=0
        special_avgn=0
      fi

      (( ${#special_name} > special_maxlen )) && special_maxlen=${#special_name}
      special_rows+=( "$special_bit|$special_name|$special_active|$special_avg|$special_scaled" )
    }

    add_special_row 10
    add_special_row 24
    add_special_row 25
    add_special_row 26
    add_special_row 27
    add_special_row 12
    add_special_row 36
    add_special_row 37
    add_special_row 38
    (( special_maxlen+=2 ))

    echo
    printf "%sRare-trigger organization for %s%s\n" "$ANSI_BOLD$ANSI_CYAN" "$LABEL" "$ANSI_RESET"
    echo

    printf "%sNo vertex cut table%s\n" "$ANSI_BOLD$ANSI_GREEN" "$ANSI_RESET"
    printf "%s  %4s | %-*s | %11s | %11s | %14s%s\n" \
      "$ANSI_BOLD" "Bit" "$special_maxlen" "TriggerName" "ActiveRuns" "AvgScale" "ScaledEventSum" "$ANSI_RESET"
    printf "%s  %-4s-+-%-*s-+-%-11s-+-%-11s-+-%-14s%s\n" \
      "$ANSI_DIM" \
      "$(printf '─%.0s' $(seq 1 4))" \
      "$special_maxlen" "$(printf '─%.0s' $(seq 1 $special_maxlen))" \
      "$(printf '─%.0s' $(seq 1 11))" \
      "$(printf '─%.0s' $(seq 1 11))" \
      "$(printf '─%.0s' $(seq 1 14))" \
      "$ANSI_RESET"

    for special_bit in 10 24 25 26 27; do
      for row in "${special_rows[@]}"; do
        IFS='|' read -r sb sn sa savg sscaled <<< "$row"
        if [[ "$sb" == "$special_bit" ]]; then
          printf "  %4d | %-*s | %11d | %11s | %14d\n" \
            "$sb" "$special_maxlen" "$sn" "$sa" "$savg" "$sscaled"
          break
        fi
      done
    done

    echo
    printf "%sVertex cut table%s\n" "$ANSI_BOLD$ANSI_GREEN" "$ANSI_RESET"
    printf "%s  %4s | %-*s | %11s | %11s | %14s%s\n" \
      "$ANSI_BOLD" "Bit" "$special_maxlen" "TriggerName" "ActiveRuns" "AvgScale" "ScaledEventSum" "$ANSI_RESET"
    printf "%s  %-4s-+-%-*s-+-%-11s-+-%-11s-+-%-14s%s\n" \
      "$ANSI_DIM" \
      "$(printf '─%.0s' $(seq 1 4))" \
      "$special_maxlen" "$(printf '─%.0s' $(seq 1 $special_maxlen))" \
      "$(printf '─%.0s' $(seq 1 11))" \
      "$(printf '─%.0s' $(seq 1 11))" \
      "$(printf '─%.0s' $(seq 1 14))" \
      "$ANSI_RESET"

    for special_bit in 12 36 37 38; do
      for row in "${special_rows[@]}"; do
        IFS='|' read -r sb sn sa savg sscaled <<< "$row"
        if [[ "$sb" == "$special_bit" ]]; then
          printf "  %4d | %-*s | %11d | %11s | %14d\n" \
            "$sb" "$special_maxlen" "$sn" "$sa" "$savg" "$sscaled"
          break
        fi
      done
    done
    echo
  else
    printf "%s[WARN] No active trigger rows found for the provided golden run list.%s\n" "$ANSI_YELLOW" "$ANSI_RESET"
    printf "%s[WARN] This usually means the query returned no rows, all rows had scaledown = -1, or the run list does not match the expected DB content.%s\n" "$ANSI_YELLOW" "$ANSI_RESET"
  fi

  if (( runs_with_no_rows > 0 || rows_bad_idx > 0 )); then
    printf "%s[WARN] QA completed with transparency warnings. Review the diagnostics block above before trusting the table blindly.%s\n" "$ANSI_YELLOW" "$ANSI_RESET"
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
