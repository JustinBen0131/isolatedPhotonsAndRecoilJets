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
# Optional gates (OFF by default — enable via CLI arguments):
#  • removeBadTowerMaps
#      Require a CEMC bad-tower map to exist in the CDB for the run.
#      When enabled:
#        – Stage 6 (“Bad-tower map available (CEMC)”) appears in the summary.
#        – Runs dropped by this gate are printed and saved to:
#          /sphenix/u/patsfan753/scratch/thesisAnalysis/dst_lists_auau/runs_missing_bad_tower_map.txt
#
#  • MBD_NS_GEQ2_VTX_150
#      Require GL1 trigger “MBD N&S >= 2, vtx < 150 cm” to be enabled (scaledown != -1).
#      When enabled:
#        – Stage 7 (“MBD_NS_geq_2_vtx_lt_150 active (scaledown != -1)”) appears in the summary.
#        – Runs dropped by this gate are printed and saved to:
#          /sphenix/u/patsfan753/scratch/thesisAnalysis/dst_lists_auau/runs_missing_MBD_NS_geq_2_vtx_lt_150.txt
#
# Notes on the summary table:
#  • Stages 2–7 are always shown. If an optional gate is not used, its stage mirrors the
#    previous stage so the table shape remains fixed.
#  • The final “8) FINAL GOLDEN (2∩3∩4∩5∩maps∩trig)” row always reflects the Stage-7 totals.
#
# Usage examples:
#   ./make_dst_lists_auau.sh
#       → Standard GOLDEN selection (no optional gates)
#
#   ./make_dst_lists_auau.sh removeBadTowerMaps
#       → GOLDEN selection + drop runs missing CEMC bad-tower maps
#
#   ./make_dst_lists_auau.sh MBD_NS_GEQ2_VTX_150
#       → GOLDEN selection + require the MBD_NS_geq_2_vtx_lt_150 trigger to be enabled
#
#   ./make_dst_lists_auau.sh removeBadTowerMaps MBD_NS_GEQ2_VTX_150
#       → GOLDEN selection + both optional gates applied (maps then trigger)
#
# Files written (and ONLY these):
#   /sphenix/u/patsfan753/scratch/thesisAnalysis/dst_lists_auau/run3goldenruns.txt
#   /sphenix/u/patsfan753/scratch/thesisAnalysis/dst_lists_auau/DST_CALOFITTING_run3auau_new_newcdbtag_v008-<run8>.list
#   (optional) runs_missing_bad_tower_map.txt
#   (optional) runs_missing_MBD_NS_geq_2_vtx_lt_150.txt
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

# Optional CLI flags:
#  • 'removeBadTowerMaps' → drop runs without a CEMC bad-tower map
#  • 'MBD_NS_GEQ2_VTX_150' → keep only runs where this trigger is enabled (scaledown != -1)
REMOVE_BAD_TOWER_MAPS=false
FILTER_MBD_NS_GEQ2_VTX_150=false
TRIG_BIT=14
TRIG_NAME_DB="MBD N&S >= 2, vtx < 150 cm"
TRIG_KEY="MBD_NS_geq_2_vtx_lt_150"

for arg in "$@"; do
  [[ "$arg" == "removeBadTowerMaps" ]] && REMOVE_BAD_TOWER_MAPS=true
  [[ "$arg" == "MBD_NS_GEQ2_VTX_150" ]] && FILTER_MBD_NS_GEQ2_VTX_150=true
done



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

# ---------- Build per-stage vectors (PROGRESSIVE gates) ----------
STAGE1_LIST=("${RUNS_ALL[@]}")
STAGE2_LIST=("${RUNS_MAGNET_ON[@]}")

# Stage 3 should be magnet∩caloQA (not caloQA alone)
declare -A Mset Cset
for r in "${RUNS_MAGNET_ON[@]}";   do Mset["$r"]=1; done
for r in "${RUNS_CALO_GOLDEN[@]}"; do Cset["$r"]=1; done

STAGE3_LIST=()
for r in "${RUNS_ALL[@]}"; do
  if [[ -n "${Mset[$r]:-}" && -n "${Cset[$r]:-}" ]]; then
    STAGE3_LIST+=("$r")
  fi
done

# Stage 4/5 should be applied after Stage 3 (not globally)
declare -A Tset Eset
for r in "${RUNS_TIME_OK[@]}";   do Tset["$r"]=1; done
for r in "${RUNS_EVENT_OK[@]}";  do Eset["$r"]=1; done

STAGE4_LIST=()
for r in "${STAGE3_LIST[@]}"; do
  if [[ -n "${Tset[$r]:-}" ]]; then
    STAGE4_LIST+=("$r")
  fi
done

STAGE5_LIST=()
for r in "${STAGE4_LIST[@]}"; do
  if [[ -n "${Eset[$r]:-}" ]]; then
    STAGE5_LIST+=("$r")
  fi
done

# Final GOLDEN starts as Stage 5, then optional gates refine RUNS_GOLDEN_FINAL
RUNS_GOLDEN_FINAL=("${STAGE5_LIST[@]}")

# ---------- Stage 6 (optional): Bad-tower map availability (CEMC) ----------
RUNS_MISSING_MAPS=()
if $REMOVE_BAD_TOWER_MAPS; then
  CDB_DIR="/cvmfs/sphenix.sdcc.bnl.gov/calibrations/sphnxpro/cdb/CEMC_BadTowerMap"
  say "Checking availability of CEMC bad-tower maps (CDB)…"
  say "CDB directory: ${BOLD}${CDB_DIR}${RESET}"

  MAP_FILE_COUNT=$(find "$CDB_DIR" -type f 2>/dev/null | wc -l | tr -d ' ')
  say "Total files under CDB dir: ${BOLD}${MAP_FILE_COUNT}${RESET}"

  say "Sampling first 30 CDB files and extracted run ids:"
  mapfile -t _MAP_FILES_SAMPLE < <(find "$CDB_DIR" -type f 2>/dev/null | head -n 30)
  idx=0
  for f in "${_MAP_FILES_SAMPLE[@]}"; do
    ((idx+=1))
    rid=$(echo "$f" | sed -E 's#.*[^0-9]([0-9]{5,8})[^0-9]*$#\1#')
    printf "    [%02d] %s -> runid=%s\n" "$idx" "$f" "${rid:-N/A}"
  done

  mapfile -t _MAP_RUNS < <(
    find "$CDB_DIR" -type f 2>/dev/null \
    | sed -E 's#.*[^0-9]([0-9]{5,8})[^0-9]*$#\1#' \
    | awk '{ printf "%d\n", $1 }' \
    | sort -u
  )

  say "Unique run-ids extracted from CDB: ${BOLD}${#_MAP_RUNS[@]}${RESET}"

  printf '%s\n' "${_MAP_RUNS[@]}"                 > "${OUT_DIR}/debug_map_runids.txt"
  printf '%s\n' "${RUNS_GOLDEN_FINAL[@]}"         > "${OUT_DIR}/debug_golden_before_mapgate.txt"

  if ((${#_MAP_RUNS[@]})); then
    say "CDB run-ids (first 20): $(printf '%s ' "${_MAP_RUNS[@]:0:20}")"
    tail_start=$(( ${#_MAP_RUNS[@]} > 20 ? ${#_MAP_RUNS[@]}-20 : 0 ))
    say "CDB run-ids (last 20):  $(printf '%s ' "${_MAP_RUNS[@]:$tail_start}")"
  fi

  declare -A MAPSET
  for r in "${_MAP_RUNS[@]}"; do MAPSET["$r"]=1; done

  RUNS_WITH_MAP=()
  RUNS_GOLDEN_FINAL_BEFORE_MAP=("${RUNS_GOLDEN_FINAL[@]}")
  for r in "${RUNS_GOLDEN_FINAL[@]}"; do
    if [[ -n "${MAPSET[$r]:-}" ]]; then
      RUNS_WITH_MAP+=("$r")
    else
      RUNS_MISSING_MAPS+=("$r")
    fi
  done
  RUNS_GOLDEN_FINAL=("${RUNS_WITH_MAP[@]}")

  printf '%s\n' "${RUNS_GOLDEN_FINAL_BEFORE_MAP[@]}" > "${OUT_DIR}/debug_golden_before_mapgate_dup.txt"
  printf '%s\n' "${RUNS_GOLDEN_FINAL[@]}"            > "${OUT_DIR}/debug_golden_after_mapgate.txt"
  printf '%s\n' "${RUNS_MISSING_MAPS[@]}"            > "${OUT_DIR}/debug_missing_map_runs.txt"

  say "Map-gate summary:"
  say "  Golden runs before map gate : ${BOLD}${#RUNS_GOLDEN_FINAL_BEFORE_MAP[@]}${RESET}"
  say "  Golden runs WITH maps       : ${BOLD}${#RUNS_GOLDEN_FINAL[@]}${RESET}"
  say "  Golden runs MISSING maps    : ${BOLD}${#RUNS_MISSING_MAPS[@]}${RESET}"

  if ((${#RUNS_GOLDEN_FINAL[@]})); then
    say "  Sample WITH maps (first 20): $(printf '%08d ' "${RUNS_GOLDEN_FINAL[@]:0:20}")"
  else
    warn "  No runs passed the map gate (WITH maps = 0)."
  fi
  if ((${#RUNS_MISSING_MAPS[@]})); then
    say "  Sample MISSING maps (first 20): $(printf '%08d ' "${RUNS_MISSING_MAPS[@]:0:20}")"
  fi

  say "Directory listing sample from CDB (depth 2, first 40 entries):"
  find "$CDB_DIR" -maxdepth 2 -type f 2>/dev/null | head -n 40 | sed 's/^/    /'

  # Note about normalization
  say "Note: run-id normalization uses integer conversion (drops leading zeros)."

  # Record Stage 6 scalars (post-map gate)
  s6_runs=${#RUNS_GOLDEN_FINAL[@]}
  read s6_ev s6_rt < <(_sum_ev_rt_over RUNS_GOLDEN_FINAL)
fi

# ---------- Stage 7 (optional): Require trigger active (scaledown != -1) ----------
RUNS_MISSING_TRIG=()
if $FILTER_MBD_NS_GEQ2_VTX_150; then
  say "Applying trigger gate: ${BOLD}${TRIG_NAME_DB}${RESET}  (bit=${TRIG_BIT}, scaledown != -1)"
  RUNS_GOLDEN_FINAL_BEFORE_TRIG=("${RUNS_GOLDEN_FINAL[@]}")
  RUNS_WITH_TRIG=()
  for r in "${RUNS_GOLDEN_FINAL[@]}"; do
    val=$(sql "SELECT scaledown${TRIG_BIT} FROM gl1_scaledown WHERE runnumber=${r};" | tr -d '[:space:]')
    # treat empty as disabled
    if [[ -n "$val" && "$val" != "-1" ]]; then
      RUNS_WITH_TRIG+=("$r")
    else
      RUNS_MISSING_TRIG+=("$r")
    fi
  done
  RUNS_GOLDEN_FINAL=("${RUNS_WITH_TRIG[@]}")

  # Record Stage 7 scalars (post-trigger gate)
  s7_runs=${#RUNS_GOLDEN_FINAL[@]}
  read s7_ev s7_rt < <(_sum_ev_rt_over RUNS_GOLDEN_FINAL)

  printf '%s\n' "${RUNS_MISSING_TRIG[@]}" > "${OUT_DIR}/runs_missing_${TRIG_KEY}.txt"
  say "Trigger gate summary:"
  say "  Golden runs WITH ${TRIG_KEY} : ${BOLD}${#RUNS_GOLDEN_FINAL[@]}${RESET}"
  say "  Golden runs MISSING trigger  : ${BOLD}${#RUNS_MISSING_TRIG[@]}${RESET}"
fi



# ---------- Per-stage scalars ----------
s1_runs=${#STAGE1_LIST[@]}; read s1_ev s1_rt < <(printf '%s\t%s\n' "$tot_ev_all" "$tot_rt_all")
s2_runs=${#STAGE2_LIST[@]}; read s2_ev s2_rt < <(_sum_ev_rt_over STAGE2_LIST)
s3_runs=${#STAGE3_LIST[@]}; read s3_ev s3_rt < <(_sum_ev_rt_over STAGE3_LIST)
s4_runs=${#STAGE4_LIST[@]}; read s4_ev s4_rt < <(_sum_ev_rt_over STAGE4_LIST)
s5_runs=${#STAGE5_LIST[@]}; read s5_ev s5_rt < <(_sum_ev_rt_over STAGE5_LIST)

# Ensure Stage 6 exists: if maps gate was not applied, Stage 6 = Stage 5
if ! $REMOVE_BAD_TOWER_MAPS; then
  s6_runs=$s5_runs
  s6_ev=$s5_ev
  s6_rt=$s5_rt
fi

# Ensure Stage 7 exists: if trigger gate was not applied, Stage 7 = Stage 6
if ! $FILTER_MBD_NS_GEQ2_VTX_150; then
  s7_runs=$s6_runs
  s7_ev=$s6_ev
  s7_rt=$s6_rt
fi

# Final “gold” must reflect the actual RUNS_GOLDEN_FINAL list (after optional gates)
gold_runs=${#RUNS_GOLDEN_FINAL[@]}
read gold_ev gold_rt < <(_sum_ev_rt_over RUNS_GOLDEN_FINAL)





pct() { local p=${1:-0} t=${2:-1}; [[ "$t" -eq 0 ]] && printf "0.00" || printf "%0.2f" "$(bc -l <<< "100.0*$p/$t")"; }

# ---------- Tabulation (clean, linearized, auto-sorted by drop size) ----------
# cuts[] rows: "key|label|runs|events|drop"
declare -a CUTS
CUTS+=("magnet|2) MAGNET ON (magnet_info=='t')|${s2_runs}|${s2_ev}|$((s1_runs - s2_runs))")
CUTS+=("caloqa|3) Calo QA GOLDEN (EMC+IHC+OHC)|${s3_runs}|${s3_ev}|$((s1_runs - s3_runs))")
CUTS+=("runtime|4) Runtime >= ${MIN_RUNTIME}s|${s4_runs}|${s4_ev}|$((s1_runs - s4_runs))")
CUTS+=("gl1|5) GL1 .evt >= ${MIN_GL1_EVT}|${s5_runs}|${s5_ev}|$((s1_runs - s5_runs))")
CUTS+=("maps|6) Bad-tower map available (CEMC)|${s6_runs}|${s6_ev}|$((s1_runs - s6_runs))")
CUTS+=("trig|7) ${TRIG_KEY} active|${s7_runs}|${s7_ev}|$((s1_runs - s7_runs))")





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
  "$(printf '%*s' "$LBLW" '' | tr ' ' '-')" \
  "$(printf '%*s' "$EVW"  '' | tr ' ' '-')" \
  "$(printf '%*s' "$RUNW" '' | tr ' ' '-')" \
  "$(printf '%*s' "$RETW" '' | tr ' ' '-')"

# Stage 1
printf "  %-*.*s | %*s | %*s | %*.2f%%\n" \
  $LBLW $LBLW "1) All (CreateDstList printruns)" \
  $EVW "$s1_ev" $RUNW "$s1_runs" $RETW 100.00

# Sorted cuts (least → most drop)
while IFS='|' read -r drop key label runs ev ; do
  rpct=$(pct "$runs" "$s1_runs")
  printf "  %-*.*s | %*s | %*s | %*.2f%%\n" \
    $LBLW $LBLW "$label" \
    $EVW "$ev" \
    $RUNW "$runs" \
    $RETW "$rpct"
done <<< "$SORTED"

# Final line (fixed shape)
_finalStage="8) FINAL GOLDEN (2+3+4+5+maps+trig)"
printf "  %-*.*s | %*s | %*s | %*.2f%%\n" \
  $LBLW $LBLW "$_finalStage" \
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
if $REMOVE_BAD_TOWER_MAPS; then
  printf "    • %-22s : %5d\n" "Missing bad-tower map" "$(( ${#RUNS_MISSING_MAPS[@]} ))"
fi
echo

# Only list missing-map runs if the gate was applied
if $REMOVE_BAD_TOWER_MAPS && ((${#RUNS_MISSING_MAPS[@]})); then
  echo "  Runs removed due to missing CEMC bad-tower maps:"
  for r in "${RUNS_MISSING_MAPS[@]}"; do
    printf "    %08d\n" "$((10#$r))"
  done
  printf '%s\n' "${RUNS_MISSING_MAPS[@]}" > "${OUT_DIR}/runs_missing_bad_tower_map.txt"
  echo "  Saved list → ${OUT_DIR}/runs_missing_bad_tower_map.txt"
  echo
fi

# Only list runs that failed the trigger gate if the gate was applied
if $FILTER_MBD_NS_GEQ2_VTX_150 && ((${#RUNS_MISSING_TRIG[@]})); then
  echo "  Runs removed due to disabled trigger (${TRIG_KEY}):"
  for r in "${RUNS_MISSING_TRIG[@]}"; do
    printf "    %08d\n" "$((10#$r))"
  done
  printf '%s\n' "${RUNS_MISSING_TRIG[@]}" > "${OUT_DIR}/runs_missing_${TRIG_KEY}.txt"
  echo "  Saved list → ${OUT_DIR}/runs_missing_${TRIG_KEY}.txt"
  echo
fi



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
