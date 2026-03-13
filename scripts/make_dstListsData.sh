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
    compgen -G "${OUT_DIR}/dst_jetcalo-*.list" || true
    compgen -G "${OUT_DIR}/DST_JETCALO-*.list" || true
  else
    compgen -G "${OUT_DIR}/${PREFIX}_${DATASET}_${TAG}-*.list" || true
  fi
}

probe_type_counts() {
  local probe_type="$1"
  local tmpdir rc files entries f

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

  mapfile -t _probe_files < <(compgen -G "${tmpdir}/${probe_type}_${DATASET}_${TAG}-*.list" || true)
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
  num_or_zero() { [[ $1 =~ ^-?[0-9]+$ ]] && printf '%s' "$1" || printf 0; }

  ROWS=()
  tot_rt_all=0
  tot_ev_all=0
  idx=0
  for run in "${RUNS_ALL[@]}"; do
    ((idx+=1))
    rt=$(sql "SELECT FLOOR(EXTRACT(EPOCH FROM (ertimestamp-brtimestamp)))::INT FROM run WHERE runnumber=$run;" | tr -d ' ')
    rt=$(num_or_zero "$rt")
    ev=$(sql "SELECT COALESCE(SUM(lastevent-firstevent+1),0)::BIGINT
              FROM filelist WHERE runnumber=$run AND filename LIKE '%GL1_physics_gl1daq%.evt';" | tr -d ' ')
    ev=$(num_or_zero "$ev")
    ROWS+=("$run"$'\t'"$rt"$'\t'"$ev")
    tot_rt_all=$(( tot_rt_all + rt ))
    tot_ev_all=$(( tot_ev_all + ev ))
    (( idx % 200 == 0 )) && echo "[INFO] Runtime/GL1 bookkeeping progress: $idx / ${#RUNS_ALL[@]}"
  done

  echo
  echo "Trigger QA summary for: $LABEL"
  echo "Input list: $LIST_FILE"
  echo "Runs in input list: ${#RUNS_ALL[@]}"
  printf "Total GL1 .evt counts: %s\n" "$tot_ev_all"
  printf "Total runtime (hours): %.2f\n" "$(bc -l <<< "$tot_rt_all/3600")"

  declare -A TRIG_NAME
  declare -A TRIG_ONRUNS
  declare -A TRIG_SUMSCALED

  gidx=0
  for run in "${RUNS_ALL[@]}"; do
    ((gidx+=1))
    declare -A SEEN_RUN=()
    while IFS=$'\t' read -r trg idx scaled; do
      [[ -z "$idx" ]] && continue
      idx=$(num_or_zero "$idx")
      scaled=$(num_or_zero "$scaled")
      if (( scaled > 0 )); then
        TRIG_NAME["$idx"]="$trg"
        TRIG_SUMSCALED["$idx"]=$(( ${TRIG_SUMSCALED["$idx"]:-0} + scaled ))
        if [[ -z "${SEEN_RUN[$idx]:-}" ]]; then
          TRIG_ONRUNS["$idx"]=$(( ${TRIG_ONRUNS["$idx"]:-0} + 1 ))
          SEEN_RUN["$idx"]=1
        fi
      fi
    done < <(sql "SELECT t.triggername, s.index, s.scaled
                   FROM gl1_scalers s
                   JOIN gl1_triggernames t
                     ON s.index = t.index
                    AND s.runnumber BETWEEN t.runnumber AND t.runnumber_last
                  WHERE s.runnumber = $run;")
    (( gidx % 200 == 0 )) && echo "[INFO] Trigger scan progress: $gidx / ${#RUNS_ALL[@]}"
  done

  tmp_rows=()
  maxlen=7
  for idx in "${!TRIG_ONRUNS[@]}"; do
    trg="${TRIG_NAME["$idx"]}"
    (( ${#trg} > maxlen )) && maxlen=${#trg}
    onr=${TRIG_ONRUNS["$idx"]}
    sum=${TRIG_SUMSCALED["$idx"]}
    tmp_rows+=( "$idx|$trg|$onr|$sum" )
  done
  (( maxlen+=2 ))

  if ((${#tmp_rows[@]})); then
    echo
    printf "  %4s | %-*s | %8s | %14s\n" "Bit" "$maxlen" "Trigger" "ONruns" "ΣScaled"
    printf "  %-4s-+-%-*s-+-%-8s-+-%-14s\n" \
      "$(printf '─%.0s' $(seq 1 4))" \
      "$maxlen" "$(printf '─%.0s' $(seq 1 $maxlen))" \
      "$(printf '─%.0s' $(seq 1 8))" \
      "$(printf '─%.0s' $(seq 1 14))"

    printf "%s\n" "${tmp_rows[@]}" \
    | awk -F'|' '{printf "%015d|%04d|%08d|%s\n",$4,$1,$3,$2}' \
    | sort -n \
    | while IFS='|' read -r padSum padIdx padOn trg; do
        sum=$((10#$padSum))
        idx=$((10#$padIdx))
        on=$((10#$padOn))
        printf "  %4d | %-*s | %8d | %14d\n" "$idx" "$maxlen" "$trg" "$on" "$sum"
      done
    echo
  else
    echo "[WARN] No trigger activity found for the provided run list."
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

  made=0
  for r in "${RUNS_ALL[@]}"; do
    r8=$(printf "%08d" "$((10#$r))")
    if ls -1 "${OUT_DIR}/${PREFIX}_${DATASET}_${TAG}-${r8}.list" >/dev/null 2>&1; then
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
