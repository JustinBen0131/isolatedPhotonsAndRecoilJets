#!/bin/bash
###############################################################################
# make_dst_lists_auau.sh — Run-3 Au+Au CALOFITTING list builder (v008, strict, verbose)
#
# Hard gates (ALL must pass):
#  • MAGNET ON (magnet_info.magnet_on == 't')
#  • Calo QA GOLDEN simultaneously for EMCal, IHCal, and OHCal (Production_write.goodruns)
#  • Runtime ≥ 300 s
#  • GL1 .evt counts ≥ 1×10^5
#
# Files written (and ONLY these):
#   /sphenix/u/patsfan753/scratch/thesisAnalysis/dst_lists_auau/run3goldenruns.txt
#   /sphenix/u/patsfan753/scratch/thesisAnalysis/dst_lists_auau/DST_CALOFITTING_run3auau_new_newcdbtag_v008-<run8>.list
###############################################################################
set -euo pipefail
IFS=$'\n\t' ; shopt -s nullglob

# ---------- Styles ----------
BOLD=$'\e[1m' ; RESET=$'\e[0m'
GREEN=$'\e[32m' ; CYAN=$'\e[36m' ; MAGENTA=$'\e[35m' ; YELLOW=$'\e[33m' ; RED=$'\e[31m' ; BLUE=$'\e[34m'
say()   { printf "${BLUE}➜${RESET} %s\n" "$*"; }
good()  { printf "${GREEN}%s${RESET}\n"   "$*"; }
warn()  { printf "${YELLOW}⚠ %s${RESET}\n" "$*" >&2; }
fatal() { printf "${RED}✘ %s${RESET}\n"   "$*" >&2; exit 1; }
trap 'fatal "Script aborted (line $LINENO) – $BASH_COMMAND"' ERR

# ---------- Fixed settings ----------
DATASET="run3auau"
TAG="new_newcdbtag_v008"
PREFIX="DST_CALOFITTING"
MIN_RUNTIME=300
MIN_GL1_EVT=100000
OUT_DIR="/sphenix/u/patsfan753/scratch/thesisAnalysis/dst_lists_auau"

# ---------- Env checks ----------
command -v CreateDstList.pl >/dev/null 2>&1 || fatal "CreateDstList.pl not found in PATH."
command -v psql            >/dev/null 2>&1 || fatal "psql not found in PATH."
command -v python3         >/dev/null 2>&1 || fatal "python3 not found in PATH."
python3 - <<'PY' >/dev/null 2>&1 || fatal "Python 'pyodbc' is required for Production_write Calo QA."
import pyodbc
PY

mkdir -p "$OUT_DIR"
# Purge OUT_DIR safely
case "$OUT_DIR" in
  /sphenix/u/patsfan753/scratch/thesisAnalysis/dst_lists_auau)
    ( shopt -s dotglob nullglob; rm -rf "${OUT_DIR}/"* )
    ;;
  *) fatal "Refusing to purge unknown OUT_DIR: $OUT_DIR" ;;
esac

say "Starting ${BOLD}Run-3 Au+Au CALOFITTING${RESET} list build"
say "Dataset: ${BOLD}${DATASET}${RESET}   Tag: ${BOLD}${TAG}${RESET}   Prefix: ${BOLD}${PREFIX}${RESET}"
good "Output directory: $OUT_DIR"

# ---------- DAQ helpers ----------
PSQL=(psql -h sphnxdaqdbreplica -d daq -At -F $'\t' -q)
sql()         { "${PSQL[@]}" -c "$1" 2>/dev/null || true; }
num_or_zero() { [[ $1 =~ ^-?[0-9]+$ ]] && printf '%s' "$1" || printf 0; }

# ---------- Stage 1: All printruns ----------
say "Collecting run numbers  (CreateDstList.pl --printruns ${PREFIX})"
mapfile -t RUNS_ALL < <(CreateDstList.pl --tag "$TAG" --dataset "$DATASET" --printruns "$PREFIX")
TOTAL_STAGE1_RUNS=${#RUNS_ALL[@]}
say "  Stage 1 count: ${TOTAL_STAGE1_RUNS} run(s) from printruns"
(( TOTAL_STAGE1_RUNS > 0 )) || fatal "No runs returned by CreateDstList.pl"

# ---------- Stage 2: Magnet ON ----------
say "Checking MAGNET status (magnet_info.magnet_on=='t') …"
RUNS_MAGNET_ON=()
rej_magnet=0
declare -A MAGNET_T
idx=0
for run in "${RUNS_ALL[@]}"; do
  ((idx+=1))
  mg=$(sql "SELECT magnet_on FROM magnet_info WHERE runnumber=$run;" | tr -d '[:space:]')
  if [[ "$mg" == "t" ]]; then
    MAGNET_T["$run"]=1
    RUNS_MAGNET_ON+=("$run")
  else
    ((rej_magnet+=1))
  fi
  (( idx % 200 == 0 )) && say "  MAGNET progress: $idx / $TOTAL_STAGE1_RUNS (kept: ${#RUNS_MAGNET_ON[@]})"
done
say "  Stage 2 count (MAGNET ON): ${#RUNS_MAGNET_ON[@]}  | dropped: ${rej_magnet}"

# ---------- Stage 3: Calo QA GOLDEN (EMC∩IHC∩OHC) ----------
say "Querying Calo QA (Production_write.goodruns) → require GOLDEN for EMCal, IHCal, and OHCal …"

CALO_QA_PY=$(cat <<'PY'
import sys, pyodbc
def log(msg): print(f"[CaloQA] {msg}", file=sys.stderr, flush=True)
runs = [int(x.strip()) for x in sys.stdin if x.strip()]
log(f"Received {len(runs)} run(s) from printruns.")
if not runs: sys.exit(0)
try:
    cn = pyodbc.connect("DSN=Production_write"); cur = cn.cursor(); log("Connected to DSN=Production_write.")
except Exception as e:
    log(f"ERROR connecting to Production_write: {e!r}"); sys.exit(0)

def has_column(name):
    try:
        cur.execute("""SELECT 1 FROM information_schema.columns
                       WHERE table_name='goodruns' AND column_name=? LIMIT 1""",(name,))
        return cur.fetchone() is not None
    except Exception as e:
        log(f"Schema probe failed for '{name}': {e!r}"); return False

style_struct = has_column('emcal_auto')
style_flat   = has_column('emcal_auto_runclass')
log(f"Schema detect → struct={style_struct}, flat={style_flat}")
if not style_struct and not style_flat:
    log("Neither style detected; will try BOTH syntaxes."); style_struct=True; style_flat=True

gold=None
for det in ('emcal','ihcal','ohcal'):
    S=set(); tried=False
    if style_struct:
        tried=True
        try:
            cur.execute(f"SELECT runnumber FROM goodruns WHERE ({det}_auto).runclass='GOLDEN'")
            rows=cur.fetchall(); log(f"{det.upper()}: struct GOLDEN rows = {len(rows)}")
            S |= {int(r.runnumber) for r in rows}
        except Exception as e:
            log(f"{det.upper()}: struct query failed: {e!r}")
    if style_flat:
        tried=True
        try:
            cur.execute(f"SELECT runnumber FROM goodruns WHERE {det}_auto_runclass='GOLDEN'")
            rows=cur.fetchall(); log(f"{det.upper()}: flat   GOLDEN rows = {len(rows)}")
            S |= {int(r.runnumber) for r in rows}
        except Exception as e:
            log(f"{det.upper()}: flat query failed: {e!r}")
    if not tried:
        log(f"{det.upper()}: no query attempted.")
    gold = S if gold is None else (gold & S)

gold_ct = 0 if gold is None else len(gold)
log(f"Intersection GOLDEN(EMC∩IHC∩OHC) size = {gold_ct}")
keep = sorted(set(runs) & (gold if gold is not None else set()))
log(f"Post-intersection with printruns: kept {len(keep)} run(s).")
for r in keep: print(r)
PY
)

mapfile -t RUNS_CALO_GOLDEN < <( printf '%s\n' "${RUNS_ALL[@]}" | python3 -c "$CALO_QA_PY" )
declare -A CALO_G; for r in "${RUNS_CALO_GOLDEN[@]}"; do CALO_G["$r"]=1; done
rej_calo=$(( TOTAL_STAGE1_RUNS - ${#RUNS_CALO_GOLDEN[@]} ))
say "  Stage 3 count (Calo QA GOLDEN EMC∩IHC∩OHC): ${#RUNS_CALO_GOLDEN[@]}  | dropped vs. Stage1: ${rej_calo}"

# ---------- Stage 4 & 5: Runtime / GL1 ----------
say "Fetching per-run runtime and GL1 .evt counts, applying thresholds …"
ROWS=() RUNS_TIME_OK=() RUNS_EVENT_OK=()
tot_rt_all=0 tot_ev_all=0
rej_short=0 rej_lowevt=0 rej_both=0
idx=0
for run in "${RUNS_ALL[@]}"; do
  ((idx+=1))
  rt=$(sql "SELECT FLOOR(EXTRACT(EPOCH FROM (ertimestamp-brtimestamp)))::INT FROM run WHERE runnumber=$run;" | tr -d ' ')
  rt=$(num_or_zero "$rt")
  ev=$(sql "SELECT COALESCE(SUM(lastevent-firstevent+1),0)::BIGINT
            FROM filelist WHERE runnumber=$run AND filename LIKE '%GL1_physics_gl1daq%.evt';" | tr -d ' ')
  ev=$(num_or_zero "$ev")
  short=$(( rt < MIN_RUNTIME ? 1 : 0 ))
  lowev=$(( ev < MIN_GL1_EVT ? 1 : 0 ))
  flag="OK"
  if (( short==1 && lowev==0 )); then flag="SHORT";        ((rej_short+=1));  fi
  if (( short==0 && lowev==1 )); then flag="LOWEVT";       ((rej_lowevt+=1)); fi
  if (( short==1 && lowev==1 )); then flag="SHORT+LOWEVT"; ((rej_both+=1));   fi
  ROWS+=("$run"$'\t'"$rt"$'\t'"$ev"$'\t'"$flag")
  tot_rt_all=$(( tot_rt_all + rt ))
  tot_ev_all=$(( tot_ev_all + ev ))
  (( rt >= MIN_RUNTIME )) && RUNS_TIME_OK+=("$run")
  (( ev >= MIN_GL1_EVT )) && RUNS_EVENT_OK+=("$run")
  (( idx % 200 == 0 )) && say "  RUNTIME/GL1 progress: $idx / $TOTAL_STAGE1_RUNS (time-ok: ${#RUNS_TIME_OK[@]}, evt-ok: ${#RUNS_EVENT_OK[@]})"
done
say "  Stage 4 count (Runtime ≥ ${MIN_RUNTIME}s): ${#RUNS_TIME_OK[@]}"
say "  Stage 5 count (GL1 .evt ≥ ${MIN_GL1_EVT}): ${#RUNS_EVENT_OK[@]}"

# ---------- Utility: sum events/runtime over a run list ----------
_sum_ev_rt_over() {
  local -n _runs=$1
  local sum_ev=0 sum_rt=0
  local line rr rt ev flag
  (( ${#_runs[@]}==0 )) && { echo -e "0\t0"; return; }
  declare -A _set; for rr in "${_runs[@]}"; do _set["$rr"]=1; done
  for line in "${ROWS[@]}"; do
    IFS=$'\t' read -r rr rt ev flag <<< "$line"
    [[ -n "${_set[$rr]:-}" ]] || continue
    sum_rt=$((sum_rt + rt))
    sum_ev=$((sum_ev + ev))
  done
  printf '%s\t%s\n' "$sum_ev" "$sum_rt"
}

# ---------- Build per-stage vectors ----------
STAGE1_LIST=("${RUNS_ALL[@]}")
STAGE2_LIST=("${RUNS_MAGNET_ON[@]}")
STAGE3_LIST=("${RUNS_CALO_GOLDEN[@]}")
STAGE4_LIST=("${RUNS_TIME_OK[@]}")
STAGE5_LIST=("${RUNS_EVENT_OK[@]}")

# Final GOLDEN = 2 ∩ 3 ∩ 4 ∩ 5
declare -A Mset Eset Tset Cset
for r in "${RUNS_MAGNET_ON[@]}";   do Mset["$r"]=1; done
for r in "${RUNS_EVENT_OK[@]}";    do Eset["$r"]=1; done
for r in "${RUNS_TIME_OK[@]}";     do Tset["$r"]=1; done
for r in "${RUNS_CALO_GOLDEN[@]}"; do Cset["$r"]=1; done
RUNS_GOLDEN_FINAL=()
for r in "${RUNS_ALL[@]}"; do
  if [[ -n "${Mset[$r]:-}" && -n "${Cset[$r]:-}" && -n "${Tset[$r]:-}" && -n "${Eset[$r]:-}" ]]; then
    RUNS_GOLDEN_FINAL+=("$r")
  fi
done

# ---------- Per-stage scalars ----------
s1_runs=${#STAGE1_LIST[@]}; read s1_ev s1_rt < <(printf '%s\t%s\n' "$tot_ev_all" "$tot_rt_all")
s2_runs=${#STAGE2_LIST[@]}; read s2_ev s2_rt < <(_sum_ev_rt_over STAGE2_LIST)
s3_runs=${#STAGE3_LIST[@]}; read s3_ev s3_rt < <(_sum_ev_rt_over STAGE3_LIST)
s4_runs=${#STAGE4_LIST[@]}; read s4_ev s4_rt < <(_sum_ev_rt_over STAGE4_LIST)
s5_runs=${#STAGE5_LIST[@]}; read s5_ev s5_rt < <(_sum_ev_rt_over STAGE5_LIST)
gold_runs=${#RUNS_GOLDEN_FINAL[@]}; read gold_ev gold_rt < <(_sum_ev_rt_over RUNS_GOLDEN_FINAL)

pct() { local p=${1:-0} t=${2:-1}; [[ "$t" -eq 0 ]] && printf "0.00" || printf "%0.2f" "$(bc -l <<< "100.0*$p/$t")"; }

# ---------- Tabulation (clean, linearized, auto-sorted by drop size) ----------
# cuts[] rows: "key|label|runs|events|drop"
declare -a CUTS
CUTS+=("magnet|2) MAGNET ON (magnet_info=='t')|${s2_runs}|${s2_ev}|$((s1_runs - s2_runs))")
CUTS+=("caloqa|3) Calo QA GOLDEN (EMC∩IHC∩OHC)|${s3_runs}|${s3_ev}|$((s1_runs - s3_runs))")
CUTS+=("runtime|4) Runtime ≥ ${MIN_RUNTIME}s|${s4_runs}|${s4_ev}|$((s1_runs - s4_runs))")
CUTS+=("gl1|5) GL1 .evt ≥ ${MIN_GL1_EVT}|${s5_runs}|${s5_ev}|$((s1_runs - s5_runs))")

# Sort by drop ascending (least → most), stable on ties by label
SORTED=$(printf '%s\n' "${CUTS[@]}" | awk -F'|' '{printf "%09d|%s|%s|%s|%s\n",$5,$1,$2,$3,$4}' | sort -n)

# ---- Pretty, fixed-width table (one set of widths everywhere) ----
LBLW=38   # stage label column width
EVW=15    # events column width
RUNW=8    # runs column width
RETW=9    # retained column width (including %)

# helper to draw a horizontal rule matching the table width
draw_hr() {
  # 2 spaces + label + " | " + events + " | " + runs + " | " + retained
  local total=$((2 + LBLW + 3 + EVW + 3 + RUNW + 3 + RETW))
  printf " %s\n" "$(printf '%*s' "$total" '' | tr ' ' '-')"
}

echo
printf " ${BOLD}${MAGENTA}CALOFITTING QA Summary  (${TAG})${RESET}\n"
draw_hr
printf "  %-*s | %*s | %*s | %*s\n" \
  $LBLW "Stage" $EVW ".evt Events" $RUNW "Runs" $RETW "(retained)"
# header underline
printf "  %-${LBLW}s-+-%-${EVW}s-+-%-${RUNW}s-+-%-${RETW}s\n" \
  "$(printf '─%.0s' $(seq $LBLW))" \
  "$(printf '─%.0s' $(seq $EVW))"  \
  "$(printf '─%.0s' $(seq $RUNW))" \
  "$(printf '─%.0s' $(seq $RETW))"

# Stage 1
printf "  %-*s | %*s | %*s | %*.2f%%\n" \
  $LBLW "1) All (CreateDstList printruns)" \
  $EVW "$s1_ev" $RUNW "$s1_runs" $RETW 100.00

# Sorted cuts (least → most drop)
while IFS='|' read -r drop key label runs ev ; do
  rpct=$(pct "$runs" "$s1_runs")
  printf "  %-*s | %*s | %*s | %*.2f%%\n" \
    $LBLW "$label" \
    $EVW "$ev" \
    $RUNW "$runs" \
    $RETW "$rpct"
done <<< "$SORTED"

# Final line
printf "  %-*s | %*s | %*s | %*.2f%%\n" \
  $LBLW "6) FINAL GOLDEN (2∩3∩4∩5)" \
  $EVW "$gold_ev" \
  $RUNW "$gold_runs" \
  $RETW "$(pct "$gold_runs" "$s1_runs")"

draw_hr
printf "  %-24s %10.2f h   |  %-24s %10.2f h\n" \
  "Total runtime (All)"   "$(bc -l <<< "$s1_rt/3600")" \
  "Total runtime (Final)" "$(bc -l <<< "$gold_rt/3600")"


echo "  Rejection tallies (not mutually exclusive):"
printf "    • %-22s : %5d\n" "MAGNET OFF/unknown" "$rej_magnet"
printf "    • %-22s : %5d\n" "Calo QA not GOLDEN" "$(( s1_runs - s3_runs ))"
printf "    • %-22s : %5d\n" "Runtime < ${MIN_RUNTIME}s" "$rej_short"
printf "    • %-22s : %5d\n" "GL1 .evt < ${MIN_GL1_EVT}" "$rej_lowevt"
printf "    • %-22s : %5d\n" "Runtime+GL1 both fail" "$rej_both"
echo

# ==================== Trigger summary (FINAL GOLDEN ONLY; no scaledown) ====================
# Aggregates over RUNS_GOLDEN_FINAL using gl1_scalers only:
#  - ONruns : number of golden runs where a trigger recorded scaled>0
#  - ΣScaled: total 'scaled' counts over golden runs
# Sorted by ΣScaled (least → most).

if (( gold_runs > 0 )); then
  say "Building trigger summary for FINAL GOLDEN runs (${gold_runs})…"

  # Accumulators
  declare -A TRIG_ONRUNS     # triggername -> runs with scaled>0
  declare -A TRIG_SUMSCALED  # triggername -> sum of 'scaled' across runs

  gidx=0
  for run in "${RUNS_GOLDEN_FINAL[@]}"; do
    ((gidx+=1))
    # Map index -> triggername and read scaled
    while IFS=$'\t' read -r trg idx scaled; do
      [[ -z "$trg" ]] && continue
      scaled=$(num_or_zero "$scaled")
      if (( scaled > 0 )); then
        TRIG_SUMSCALED["$trg"]=$(( ${TRIG_SUMSCALED["$trg"]:-0} + scaled ))
        TRIG_ONRUNS["$trg"]=$(( ${TRIG_ONRUNS["$trg"]:-0} + 1 ))
      fi
    done < <(sql "SELECT t.triggername, s.index, s.scaled
                   FROM gl1_scalers s
                   JOIN gl1_triggernames t
                     ON s.index = t.index
                    AND s.runnumber BETWEEN t.runnumber AND t.runnumber_last
                  WHERE s.runnumber = $run;")
    (( gidx % 200 == 0 )) && say "  Trigger scan progress: $gidx / $gold_runs"
  done

  # Prepare rows: trigger | ONruns | ΣScaled ; sort by ΣScaled asc
  tmp_rows=()
  maxlen=7
  for trg in "${!TRIG_ONRUNS[@]}"; do
    (( ${#trg} > maxlen )) && maxlen=${#trg}
    onr=${TRIG_ONRUNS["$trg"]}
    sum=${TRIG_SUMSCALED["$trg"]}
    tmp_rows+=( "$trg|$onr|$sum" )
  done
  (( maxlen+=2 ))

  if ((${#tmp_rows[@]})); then
    echo
    echo " ${BOLD}${GREEN}Trigger totals over FINAL GOLDEN runs (least → most events)${RESET}"
    hdr_fmt="  %-$(printf %d "$maxlen")s | %8s | %14s\n"
    printf "$hdr_fmt" "Trigger" "ONruns" "ΣScaled"
    printf "  %-${maxlen}s-+-%-8s-+-%-14s\n" \
      "$(printf '─%.0s' $(seq 1 $maxlen))" \
      "$(printf '─%.0s' $(seq 1 8))" \
      "$(printf '─%.0s' $(seq 1 14))"

    printf "%s\n" "${tmp_rows[@]}" \
    | awk -F'|' '{printf "%015d|%08d|%s\n",$3,$2,$1}' \
    | sort -n \
    | while IFS='|' read -r padSum padOn trg; do
        sum=$((10#$padSum))
        on=$((10#$padOn))
        printf "  %-$(printf %d "$maxlen")s | %8d | %14d\n" "$trg" "$on" "$sum"
      done
    echo
  else
    warn "No trigger activity found within FINAL GOLDEN runs."
  fi
fi
# ================== End trigger summary (FINAL GOLDEN ONLY; no scaledown) ==================



# ---------- Write manifest ----------
GOLDEN_TXT="${OUT_DIR}/run3goldenruns.txt"
PERSONAL_LIST_PATH="/sphenix/u/patsfan753/scratch/thesisAnalysis/Full_AuAuGoldenRunList_Version1_personal.list"

# Write the final golden runs to the standard manifest
printf '%s\n' "${RUNS_GOLDEN_FINAL[@]}" > "$GOLDEN_TXT"

# Also write a personal .list copy at the requested path
printf '%s\n' "${RUNS_GOLDEN_FINAL[@]}" > "$PERSONAL_LIST_PATH"

good "Golden run manifest → $GOLDEN_TXT  (${gold_runs} runs)"
good "Personal .list      → $PERSONAL_LIST_PATH"



# ---------- Build lists (guard empty) ----------
say "Creating ${PREFIX} per-run lists  (tag=${TAG}, dataset=${DATASET})"
(
  cd "$OUT_DIR" || fatal "Cannot cd into $OUT_DIR"
  if [[ -s "$GOLDEN_TXT" && "$gold_runs" -gt 0 ]]; then
    CreateDstList.pl --tag "$TAG" --dataset "$DATASET" --list "$GOLDEN_TXT" "$PREFIX"
  else
    warn "Golden manifest is empty — skipping CreateDstList.pl"
  fi
)

# ---------- Verify coverage ----------
made_count=0
for r in "${RUNS_GOLDEN_FINAL[@]}"; do
  r8=$(printf "%08d" "$((10#$r))")
  if ls -1 "${OUT_DIR}/${PREFIX}_${DATASET}_${TAG}-${r8}.list" >/dev/null 2>&1; then
    ((made_count+=1))
  fi
done
miss_count=$((gold_runs - made_count))

echo
echo " ${BOLD}${MAGENTA}CreateDstList Outcome${RESET}"
printf "  Golden runs requested : %6d\n" "$gold_runs"
printf "  Per-run lists created : %6d\n" "$made_count"
if (( miss_count > 0 )); then
  printf "  Runs without a list   : %6d\n" "$miss_count"
else
  good "All GOLDEN runs produced list files."
fi

echo
echo " ${BOLD}${CYAN}Final Outputs${RESET}"
echo "   • GOLDEN run manifest : $GOLDEN_TXT"
echo "   • Per-run lists       : $OUT_DIR/${PREFIX}_*_${DATASET}_${TAG}-<run>.list"
echo
good "Done."
