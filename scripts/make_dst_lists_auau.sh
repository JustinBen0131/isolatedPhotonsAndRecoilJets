#!/bin/bash
###############################################################################
# make_dst_lists_auau.sh — Run-3 Au+Au CALOFITTING list builder (new_newcdbtag_v008)
# - Collect run list via CreateDstList.pl
# - QA: keep runs with runtime ≥ 300 s AND GL1 events ≥ 1e5
# - Generate per-run DST_CALOFITTING lists for the golden runs
###############################################################################
set -euo pipefail
IFS=$'\n\t' ; shopt -s nullglob

########################  COLOUR / LOG HELPERS  ###############################
esc=$'\e['
clr_red=${esc}0\;31m ; clr_grn=${esc}0\;32m
clr_yel=${esc}1\;33m ; clr_blu=${esc}1\;34m
clr_bld=${esc}1m     ; clr_rst=${esc}0m

say()   { printf "${clr_blu}➜${clr_rst} %s\n" "$*"; }
good()  { printf "${clr_grn}%s${clr_rst}\n"   "$*"; }
warn()  { printf "${clr_yel}⚠ %s${clr_rst}\n" "$*" >&2; }
fatal() { printf "${clr_red}✘ %s${clr_rst}\n" "$*" >&2; exit 1; }

trap 'fatal "Script aborted (line $LINENO) – $BASH_COMMAND"' ERR

############################  CONSTANTS  #####################################
# Output directory for all .list files:
list_dir="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/dst_list"

# Fixed production & dataset (as requested):
dataset="run3auau"
tag="new_newcdbtag_v008"
prefix="DST_CALOFITTING"

# QA thresholds:
min_runtime=300          # seconds (≥ 5 minutes)
min_gl1_evt=100000       # GL1 events (≥ 1e5)

###########################  PRE-FLIGHT CHECKS  ###############################
command -v CreateDstList.pl >/dev/null 2>&1 || fatal "CreateDstList.pl not found in PATH (source the sPHENIX env)."
command -v psql            >/dev/null 2>&1 || warn  "psql not found in PATH (will fail at DB query)."

rm -rf "${list_dir:?}"/* 2>/dev/null || true
mkdir -p "$list_dir"
good "Output directory prepared: $list_dir"

say "Building CALOFITTING lists for ${clr_bld}${tag}${clr_rst} (${dataset})"

##############################  HELPERS  ######################################
PSQL=(psql -h sphnxdaqdbreplica -d daq -At -F $'\t' -q)
sql()          { "${PSQL[@]}" -c "$1" 2>/dev/null || true; }
num_or_zero()  { [[ $1 =~ ^-?[0-9]+$ ]] && printf '%s' "$1" || printf 0; }

###########################  1) COLLECT RUNS  #################################
runs_file="${list_dir}/${dataset}-${tag}.runs"
say "Collecting run numbers via CreateDstList.pl --printruns"
CreateDstList.pl --tag "$tag" --dataset "$dataset" --printruns "$prefix" >"$runs_file"
total_runs=$(wc -l <"$runs_file" || echo 0)
(( total_runs > 0 )) || fatal "No runs returned by CreateDstList.pl"
good "Got $total_runs runs → $runs_file"

mapfile -t runs <"$runs_file"

###########################  2) QA FILTERING  #################################
say "Running QA (runtime ≥ ${min_runtime}s AND GL1 ≥ ${min_gl1_evt}) …"

rows=()
golden=()
tot_rt=0
tot_ev=0
fail_cnt=0

for run in "${runs[@]}"; do
  # Runtime (seconds)
  rt=$(sql "SELECT FLOOR(EXTRACT(EPOCH FROM (ertimestamp-brtimestamp)))::INT
            FROM run WHERE runnumber=$run;" | tr -d ' ')
  rt=$(num_or_zero "$rt")

  # GL1 events (from .evt filelist)
  ev=$(sql "SELECT COALESCE(SUM(lastevent-firstevent+1),0)::BIGINT
            FROM filelist WHERE runnumber=$run
             AND filename LIKE '%GL1_physics_gl1daq%.evt';" | tr -d ' ')
  ev=$(num_or_zero "$ev")

  note="GOOD"
  (( rt < min_runtime )) && note="SHORT"
  (( ev < min_gl1_evt )) && note=$([[ $note == GOOD ]] && echo "LOWEVT" || echo "${note}+LOWEVT")

  rows+=("$run"$'\t'"$rt"$'\t'"$ev"$'\t'"$note")
  tot_rt=$(( tot_rt + rt ))
  tot_ev=$(( tot_ev + ev ))

  if [[ $note == GOOD ]]; then
    golden+=("$run")
  else
    ((++fail_cnt))
  fi
done

# Pretty print per-run table
echo -e "\n${clr_grn}================ Per-run summary ================${clr_rst}"
printf '%s\n' "Run    rt[s]    GL1ev    Note"
printf '%s\n' "${rows[@]}" | column -t -s $'\t'

# Cut-flow
all_runs=${#runs[@]}
keep_runs=${#golden[@]}

echo -e "\n${clr_grn}================ Cut-flow summary ===============${clr_rst}"
printf '%-10s %6s %12s %12s\n' "Stage" "Runs" "GL1evt" "Runtime[h]"
printf '%-10s %6d %12d %12.2f\n' "All"   "$all_runs" "$tot_ev" "$(bc -l <<<"$tot_rt/3600")"

# After runtime cut
timeRuns=0; timeRt=0; timeEv=0
while IFS=$'\t' read -r _ rt ev note; do
  [[ $note == SHORT* ]] && continue
  ((timeRuns++)); timeRt=$((timeRt+rt)); timeEv=$((timeEv+ev))
done <<<"$(printf '%s\n' "${rows[@]}")"
printf '%-10s %6d %12d %12.2f\n' "Time≥5m" "$timeRuns" "$timeEv" "$(bc -l <<<"$timeRt/3600")"

# After both cuts (golden)
eventRuns=$keep_runs; eventEv=0; eventRt=0
for run in "${golden[@]}"; do
  # Lookup values from rows (small N; fine to loop)
  while IFS=$'\t' read -r rr rt ev note; do
    [[ $rr == "$run" && $note == GOOD ]] || continue
    eventRt=$((eventRt+rt)); eventEv=$((eventEv+ev))
  done <<<"$(printf '%s\n' "${rows[@]}")"
done
printf '%-10s %6d %12d %12.2f\n' "Evt≥1e5" "$eventRuns" "$eventEv" "$(bc -l <<<"$eventRt/3600")"

# Write golden manifest
golden_txt="${list_dir}/run3goldenruns.txt"
printf '%s\n' "${golden[@]}" >"$golden_txt"
good "Golden run list → $golden_txt   (${#golden[@]} runs)"

###########################  3) BUILD LISTS  ##################################
say "Creating ${prefix} lists via CreateDstList.pl (tag=${tag}, dataset=${dataset})"
(
  cd "$list_dir" || fatal "Cannot cd into $list_dir"
  CreateDstList.pl --tag "$tag" \
                   --dataset "$dataset" \
                   --list "$golden_txt" \
                   "$prefix"
)

# Quick count of produced lists
made=$(ls -1 "$list_dir"/${prefix}_*_"$tag"-*.list 2>/dev/null | wc -l || true)
say "Produced ${made} ${prefix} list file(s) in ${clr_bld}$list_dir${clr_rst}"

echo
good "All done."
