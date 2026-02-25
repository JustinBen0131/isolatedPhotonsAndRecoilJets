#!/bin/bash
# make_dst_lists_pp.sh — generate per-run DST lists (no side files)
# Production: ana509_2024p022_v001
# DST type   : DST_CALOFITTING

set -uo pipefail  # (intentionally NOT using -e)

IN_BASE="/sphenix/u/patsfan753/scratch/thesisAnalysis"
LIST_FILE="$IN_BASE/Full_ppGoldenRunList_Version3.list"
OUT_DIR="$IN_BASE/dst_lists_pp"

TAG="ana521_2025p007_v001"
TYPE="DST_JETCALO"

# Optional mode: ./make_dst_lists_pp.sh run25pp
#   - writes to dst_lists_pp_run25
#   - generates ONLY CALOFITTING list files (no tracking DSTs)
#   - uses the CALOFITTING production tag used by the run25pp macro snippet
MODE="${1:-}"

# ----------------------------------------------------------------------------
# QUERY MODE: checkTriggers
#   - NO cleaning
#   - NO writing any dst lists
#   - Compare run24 vs run25 after HARD GRL gates:
#       • MAGNET ON
#       • Calo QA GOLDEN (EMC∩IHC∩OHC)
#       • Runtime ≥ 300 s
#       • GL1 .evt counts ≥ 1×10^5
#     then compare how many of those runs have:
#       • Photon 4 GeV + MBD N&S >= 1 active  (scaledown != -1)
# ----------------------------------------------------------------------------
if [[ "$MODE" == "checkTriggers" ]]; then
  set -euo pipefail
  IFS=$'\n\t'

  # ---- Inputs (fixed) ----
  RUN24_LIST="$LIST_FILE"

  # Run25 baseline universe should come from printruns for the CALOFITTING production tag/dataset.
  RUN25_DATASET="run3pp"
  RUN25_TAG="new_newcdbtag_v008"
  RUN25_PREFIX="DST_CALOFITTING"

  command -v CreateDstList.pl >/dev/null 2>&1 || { echo "[ERROR] CreateDstList.pl not in PATH"; exit 1; }
  command -v psql            >/dev/null 2>&1 || { echo "[ERROR] psql not in PATH"; exit 1; }
  command -v python3         >/dev/null 2>&1 || { echo "[ERROR] python3 not in PATH"; exit 1; }
  python3 - <<'PY' >/dev/null 2>&1 || { echo "[ERROR] Python 'pyodbc' is required for Production_write Calo QA."; exit 1; }
import pyodbc
PY

  [[ -s "$RUN24_LIST" ]] || { echo "[ERROR] Run24 list missing/empty: $RUN24_LIST"; exit 2; }

  # ---- DAQ helpers ----
  PSQL=(psql -h sphnxdaqdbreplica -d daq -At -F $'\t' -q)
  sql()         { "${PSQL[@]}" -c "$1" 2>/dev/null || true; }
  num_or_zero() { [[ $1 =~ ^-?[0-9]+$ ]] && printf '%s' "$1" || printf 0; }

  MIN_RUNTIME=300
  MIN_GL1_EVT=100000

  # ---- Resolve the trigger bit for "Photon 4 GeV + MBD N&S >= 1" ----
  TRIG4_BIT=$(sql "SELECT index FROM gl1_triggernames WHERE triggername ILIKE '%Photon 4%' AND triggername ILIKE '%MBD%' AND triggername ILIKE '%>= 1%' ORDER BY index LIMIT 1;" | tr -d '[:space:]')
  [[ -n "$TRIG4_BIT" ]] || { echo "[ERROR] Could not resolve trigger bit for Photon 4 GeV + MBD N&S >= 1 from gl1_triggernames."; exit 3; }
  TRIG4_NAME=$(sql "SELECT triggername FROM gl1_triggernames WHERE index=${TRIG4_BIT} ORDER BY runnumber DESC LIMIT 1;" | tr -d '\r')
  [[ -n "$TRIG4_NAME" ]] || TRIG4_NAME="Photon 4 GeV + MBD N&S >= 1"

  # ---- Calo QA GOLDEN query blob (shared) ----
  CALO_QA_PY=$(cat <<'PY'
import sys, pyodbc
runs = [int(x.strip()) for x in sys.stdin if x.strip()]
if not runs: sys.exit(0)
try:
    cn = pyodbc.connect("DSN=Production_write"); cur = cn.cursor()
except Exception:
    sys.exit(0)

def has_column(name):
    try:
        cur.execute("""SELECT 1 FROM information_schema.columns
                       WHERE table_name='goodruns' AND column_name=? LIMIT 1""",(name,))
        return cur.fetchone() is not None
    except Exception:
        return False

style_struct = has_column('emcal_auto')
style_flat   = has_column('emcal_auto_runclass')
if not style_struct and not style_flat:
    style_struct=True; style_flat=True

gold=None
for det in ('emcal','ihcal','ohcal'):
    S=set()
    if style_struct:
        try:
            cur.execute(f"SELECT runnumber FROM goodruns WHERE ({det}_auto).runclass='GOLDEN'")
            S |= {int(r.runnumber) for r in cur.fetchall()}
        except Exception:
            pass
    if style_flat:
        try:
            cur.execute(f"SELECT runnumber FROM goodruns WHERE {det}_auto_runclass='GOLDEN'")
            S |= {int(r.runnumber) for r in cur.fetchall()}
        except Exception:
            pass
    gold = S if gold is None else (gold & S)

keep = sorted(set(runs) & (gold if gold is not None else set()))
for r in keep: print(r)
PY
  )

  # ---- read a run list file into an array of ints (strips comments) ----
  read_runs_from_file() {
    local f="$1"
    local -n out="$2"
    out=()
    while IFS= read -r raw; do
      rn="${raw%%#*}"
      rn="$(echo -n "$rn" | tr -d '[:space:]')"
      [[ -z "$rn" ]] && continue
      rn_dec=$(printf "%d" "$rn" 2>/dev/null || echo "$rn")
      [[ "$rn_dec" =~ ^[0-9]+$ ]] || continue
      out+=( "$rn_dec" )
    done < "$f"
  }

  # ---- compute pipeline for a given run list ----
  run_pipeline() {
    local label="$1"
    local -n input_runs="$2"

    local RUNS_ALL=("${input_runs[@]}")

    # Stage 2: Magnet ON
    local RUNS_MAGNET_ON=()
    local rej_magnet=0
    for r in "${RUNS_ALL[@]}"; do
      mg=$(sql "SELECT magnet_on FROM magnet_info WHERE runnumber=$r;" | tr -d '[:space:]')
      if [[ "$mg" == "t" ]]; then
        RUNS_MAGNET_ON+=( "$r" )
      else
        ((rej_magnet+=1))
      fi
    done

    # Stage 3: Calo QA GOLDEN (EMC∩IHC∩OHC)
    mapfile -t RUNS_CALO_GOLDEN < <( printf '%s\n' "${RUNS_ALL[@]}" | python3 -c "$CALO_QA_PY" )

    # Build Stage 3 = magnet ∩ caloQA (progressive)
    declare -A Mset Cset
    for r in "${RUNS_MAGNET_ON[@]}";   do Mset["$r"]=1; done
    for r in "${RUNS_CALO_GOLDEN[@]}"; do Cset["$r"]=1; done

    local STAGE3=()
    for r in "${RUNS_ALL[@]}"; do
      if [[ -n "${Mset[$r]:-}" && -n "${Cset[$r]:-}" ]]; then
        STAGE3+=( "$r" )
      fi
    done

    # Stage 4/5: runtime + GL1 .evt thresholds (computed per-run)
    local ROWS=()
    local RUNS_TIME_OK=()
    local RUNS_EVT_OK=()
    local tot_ev_all=0
    local tot_rt_all=0
    local rej_short=0
    local rej_lowevt=0
    local rej_both=0

    for r in "${RUNS_ALL[@]}"; do
      rt=$(sql "SELECT FLOOR(EXTRACT(EPOCH FROM (ertimestamp-brtimestamp)))::INT FROM run WHERE runnumber=$r;" | tr -d ' ')
      rt=$(num_or_zero "$rt")
      ev=$(sql "SELECT COALESCE(SUM(lastevent-firstevent+1),0)::BIGINT
                FROM filelist WHERE runnumber=$r AND filename LIKE '%GL1_physics_gl1daq%.evt';" | tr -d ' ')
      ev=$(num_or_zero "$ev")

      short=$(( rt < MIN_RUNTIME ? 1 : 0 ))
      lowev=$(( ev < MIN_GL1_EVT ? 1 : 0 ))
      if (( short==1 && lowev==0 )); then ((rej_short+=1));  fi
      if (( short==0 && lowev==1 )); then ((rej_lowevt+=1)); fi
      if (( short==1 && lowev==1 )); then ((rej_both+=1));   fi

      ROWS+=("$r"$'\t'"$rt"$'\t'"$ev")
      tot_rt_all=$(( tot_rt_all + rt ))
      tot_ev_all=$(( tot_ev_all + ev ))
      (( rt >= MIN_RUNTIME )) && RUNS_TIME_OK+=( "$r" )
      (( ev >= MIN_GL1_EVT )) && RUNS_EVT_OK+=( "$r" )
    done

    _sum_ev_over() {
      local -n _runs=$1
      local sum_ev=0
      local rr rt ev
      (( ${#_runs[@]}==0 )) && { echo 0; return; }
      declare -A _set; for rr in "${_runs[@]}"; do _set["$rr"]=1; done
      for line in "${ROWS[@]}"; do
        IFS=$'\t' read -r rr rt ev <<< "$line"
        [[ -n "${_set[$rr]:-}" ]] || continue
        sum_ev=$((sum_ev + ev))
      done
      echo "$sum_ev"
    }

    # Stage 4 = Stage3 ∩ runtimeOK
    declare -A Tset Eset
    for r in "${RUNS_TIME_OK[@]}"; do Tset["$r"]=1; done
    for r in "${RUNS_EVT_OK[@]}";  do Eset["$r"]=1; done

    local STAGE4=()
    for r in "${STAGE3[@]}"; do
      [[ -n "${Tset[$r]:-}" ]] && STAGE4+=( "$r" )
    done

    # Stage 5 = Stage4 ∩ evtOK
    local STAGE5=()
    for r in "${STAGE4[@]}"; do
      [[ -n "${Eset[$r]:-}" ]] && STAGE5+=( "$r" )
    done

    # Stage 6 = Stage5 ∩ (Photon4 trigger active: scaledown != -1)
    local STAGE6=()
    for r in "${STAGE5[@]}"; do
      v=$(sql "SELECT scaledown${TRIG4_BIT} FROM gl1_scaledown WHERE runnumber=${r};" | tr -d '[:space:]')
      if [[ -n "$v" && "$v" != "-1" ]]; then
        STAGE6+=( "$r" )
      fi
    done

    # Scalars (runs + events at key stages)
    local s1_runs=${#RUNS_ALL[@]}
    local s2_runs=${#RUNS_MAGNET_ON[@]}
    local s3_runs=${#STAGE3[@]}
    local s4_runs=${#STAGE4[@]}
    local s5_runs=${#STAGE5[@]}
    local s6_runs=${#STAGE6[@]}

    local s1_ev="$tot_ev_all"
    local s5_ev; s5_ev="$(_sum_ev_over STAGE5)"
    local s6_ev; s6_ev="$(_sum_ev_over STAGE6)"

    printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n" \
      "$label" "$s1_runs" "$s2_runs" "$s3_runs" "$s4_runs" "$s5_runs" "$s6_runs" "$s1_ev" "$s5_ev" "$s6_ev"
  }

  # ---- Run both datasets ----
  read_runs_from_file "$RUN24_LIST" RUN24_RUNS

  mapfile -t RUN25_RUNS < <(CreateDstList.pl --tag "$RUN25_TAG" --dataset "$RUN25_DATASET" --printruns "$RUN25_PREFIX")
  (( ${#RUN25_RUNS[@]} > 0 )) || { echo "[ERROR] Run25 printruns returned 0 runs (tag=${RUN25_TAG}, dataset=${RUN25_DATASET}, prefix=${RUN25_PREFIX})"; exit 2; }

  # Header
  echo
  echo "checkTriggers: Run24 vs Run25 (NO cleaning, NO list writes)"
  echo "  Trigger gate: ${TRIG4_NAME}  (bit=${TRIG4_BIT}, scaledown != -1)"
  echo

  # Compute stats rows
  row24="$(run_pipeline "run24pp" RUN24_RUNS)"
  row25="$(run_pipeline "run25pp" RUN25_RUNS)"

  # Pretty table
  printf "  %-8s | %8s | %8s | %8s | %8s | %8s | %8s || %14s | %14s | %14s\n" \
    "dataset" "S1runs" "S2mag" "S3qa" "S4rt" "S5gl1" "S6trig" "S1ev(all)" "S5ev(gl1)" "S6ev(trig)"
  printf "  %-8s-+-%8s-+-%8s-+-%8s-+-%8s-+-%8s-+-%8s-++-%14s-+-%14s-+-%14s\n" \
    "--------" "--------" "--------" "--------" "--------" "--------" "--------" "--------------" "--------------" "--------------"

  IFS=$'\t' read -r lab s1 s2 s3 s4 s5 s6 e1 e5 e6 <<< "$row24"
  printf "  %-8s | %8d | %8d | %8d | %8d | %8d | %8d || %14s | %14s | %14s\n" \
    "$lab" "$s1" "$s2" "$s3" "$s4" "$s5" "$s6" "$e1" "$e5" "$e6"

  IFS=$'\t' read -r lab s1 s2 s3 s4 s5 s6 e1 e5 e6 <<< "$row25"
  printf "  %-8s | %8d | %8d | %8d | %8d | %8d | %8d || %14s | %14s | %14s\n" \
    "$lab" "$s1" "$s2" "$s3" "$s4" "$s5" "$s6" "$e1" "$e5" "$e6"

  echo
  echo "Stages:"
  echo "  S1runs  = runs from the input list file"
  echo "  S2mag   = magnet_info.magnet_on == 't'"
  echo "  S3qa    = (S2) ∩ Calo QA GOLDEN (EMCal ∩ IHCal ∩ OHCal)"
  echo "  S4rt    = (S3) ∩ runtime ≥ ${MIN_RUNTIME}s"
  echo "  S5gl1   = (S4) ∩ GL1 .evt ≥ ${MIN_GL1_EVT}"
  echo "  S6trig  = (S5) ∩ trigger active (scaledown != -1)"
  echo

  exit 0
fi

if [[ "$MODE" == "run25pp" ]]; then
  OUT_DIR="$IN_BASE/dst_lists_pp_run25"
  TAG="new_newcdbtag_v008"
  TYPE="DST_CALOFITTING"
fi

CREATE_DST_TOOL="$(command -v CreateDstList.pl || true)"
[[ -n "$CREATE_DST_TOOL" ]] || { echo "[ERROR] CreateDstList.pl not in PATH"; exit 1; }

# ============================================================================
# run25pp mode: AuAu-style GOLDEN checks + tabulation + manifest + CreateDstList
#   - ignores LIST_FILE (uses printruns for this dataset/tag/prefix)
#   - applies the same HARD gates as make_dst_lists_auau.sh:
#       • MAGNET ON
#       • Calo QA GOLDEN (EMC∩IHC∩OHC) via Production_write.goodruns (pyodbc)
#       • Runtime ≥ 300 s
#       • GL1 .evt counts ≥ 1×10^5
#   - writes:
#       • ${OUT_DIR}/run3goldenruns.txt
#       • ${OUT_DIR}/${PREFIX}_${DATASET}_${TAG}-<run8>.list
# ============================================================================
if [[ "$MODE" == "run25pp" ]]; then
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

  # ---------- Fixed settings (run25pp) ----------
  DATASET="run3pp"
  PREFIX="DST_CALOFITTING"
  MIN_RUNTIME=300
  MIN_GL1_EVT=100000

  # ---------- Env checks ----------
  command -v CreateDstList.pl >/dev/null 2>&1 || fatal "CreateDstList.pl not found in PATH."
  command -v psql            >/dev/null 2>&1 || fatal "psql not found in PATH."
  command -v python3         >/dev/null 2>&1 || fatal "python3 not found in PATH."
  python3 - <<'PY' >/dev/null 2>&1 || fatal "Python 'pyodbc' is required for Production_write Calo QA."
import pyodbc
PY

  mkdir -p "$OUT_DIR"

  # Safety: NEVER touch baseline dst_lists_pp when running run25pp
  BASELINE_OUT_DIR="${IN_BASE}/dst_lists_pp"
  if [[ "$OUT_DIR" == "$BASELINE_OUT_DIR" ]]; then
    fatal "run25pp mode refuses to operate on baseline OUT_DIR: $OUT_DIR"
  fi

  # Purge OUT_DIR safely (run25pp only)
  RUN25_OUT_DIR="${IN_BASE}/dst_lists_pp_run25"
  if [[ "$OUT_DIR" == "$RUN25_OUT_DIR" ]]; then
    say "Purging output directory: $OUT_DIR"
    ( shopt -s dotglob nullglob; rm -rf "${OUT_DIR}/"* )
  else
    fatal "Refusing to purge unknown OUT_DIR: $OUT_DIR"
  fi

  say "Starting ${BOLD}Run-3 pp CALOFITTING${RESET} list build (run25pp mode)"
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
  idx=0
  for run in "${RUNS_ALL[@]}"; do
    ((idx+=1))
    mg=$(sql "SELECT magnet_on FROM magnet_info WHERE runnumber=$run;" | tr -d '[:space:]')
    if [[ "$mg" == "t" ]]; then
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

keep = sorted(set(runs) & (gold if gold is not None else set()))
log(f"Post-intersection with printruns: kept {len(keep)} run(s).")
for r in keep: print(r)
PY
  )

  mapfile -t RUNS_CALO_GOLDEN < <( printf '%s\n' "${RUNS_ALL[@]}" | python3 -c "$CALO_QA_PY" )
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

  declare -A Mset Cset
  for r in "${RUNS_MAGNET_ON[@]}";   do Mset["$r"]=1; done
  for r in "${RUNS_CALO_GOLDEN[@]}"; do Cset["$r"]=1; done

  STAGE3_LIST=()
  for r in "${RUNS_ALL[@]}"; do
    if [[ -n "${Mset[$r]:-}" && -n "${Cset[$r]:-}" ]]; then
      STAGE3_LIST+=("$r")
    fi
  done

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

  # ---------- Trigger gates (required): bit25 AND bit10 must both be active ----------
  TRIG25_BIT=25
  TRIG25_LABEL="Photon 3 GeV + MBD NS >= 1"
  TRIG10_BIT=10
  TRIG10_LABEL="MBD N&S >= 2"

  say "Applying trigger gate: ${BOLD}${TRIG25_LABEL}${RESET}  (bit=${TRIG25_BIT}, scaledown != -1)"
  STAGE6_LIST=()
  RUNS_MISSING_BIT25=()
  idx=0
  for r in "${STAGE5_LIST[@]}"; do
    ((idx+=1))
    val=$(sql "SELECT scaledown${TRIG25_BIT} FROM gl1_scaledown WHERE runnumber=${r};" | tr -d '[:space:]')
    if [[ -n "$val" && "$val" != "-1" ]]; then
      STAGE6_LIST+=("$r")
    else
      RUNS_MISSING_BIT25+=("$r")
    fi
    (( idx % 200 == 0 )) && say "  Bit25 progress: $idx / ${#STAGE5_LIST[@]} (kept: ${#STAGE6_LIST[@]})"
  done

  say "Applying trigger gate: ${BOLD}${TRIG10_LABEL}${RESET}  (bit=${TRIG10_BIT}, scaledown != -1)"
  STAGE7_LIST=()
  RUNS_MISSING_BIT10=()
  idx=0
  for r in "${STAGE6_LIST[@]}"; do
    ((idx+=1))
    val=$(sql "SELECT scaledown${TRIG10_BIT} FROM gl1_scaledown WHERE runnumber=${r};" | tr -d '[:space:]')
    if [[ -n "$val" && "$val" != "-1" ]]; then
      STAGE7_LIST+=("$r")
    else
      RUNS_MISSING_BIT10+=("$r")
    fi
    (( idx % 200 == 0 )) && say "  Bit10 progress: $idx / ${#STAGE6_LIST[@]} (kept: ${#STAGE7_LIST[@]})"
  done

  rej_trig25=${#RUNS_MISSING_BIT25[@]}
  rej_trig10=${#RUNS_MISSING_BIT10[@]}

  RUNS_GOLDEN_FINAL=("${STAGE7_LIST[@]}")

  # ---------- Per-stage scalars ----------
  s1_runs=${#STAGE1_LIST[@]}; read s1_ev s1_rt < <(printf '%s\t%s\n' "$tot_ev_all" "$tot_rt_all")
  s2_runs=${#STAGE2_LIST[@]}; read s2_ev s2_rt < <(_sum_ev_rt_over STAGE2_LIST)
  s3_runs=${#STAGE3_LIST[@]}; read s3_ev s3_rt < <(_sum_ev_rt_over STAGE3_LIST)
  s4_runs=${#STAGE4_LIST[@]}; read s4_ev s4_rt < <(_sum_ev_rt_over STAGE4_LIST)
  s5_runs=${#STAGE5_LIST[@]}; read s5_ev s5_rt < <(_sum_ev_rt_over STAGE5_LIST)
  s6_runs=${#STAGE6_LIST[@]}; read s6_ev s6_rt < <(_sum_ev_rt_over STAGE6_LIST)
  s7_runs=${#STAGE7_LIST[@]}; read s7_ev s7_rt < <(_sum_ev_rt_over STAGE7_LIST)

  gold_runs=${#RUNS_GOLDEN_FINAL[@]}
  read gold_ev gold_rt < <(_sum_ev_rt_over RUNS_GOLDEN_FINAL)

  pct() { local p=${1:-0} t=${2:-1}; [[ "$t" -eq 0 ]] && printf "0.00" || printf "%0.2f" "$(bc -l <<< "100.0*$p/$t")"; }

  # ---------- Tabulation ----------
  LBLW=38
  EVW=15
  RUNW=8
  RETW=9
  draw_hr() {
    local total=$((2 + LBLW + 3 + EVW + 3 + RUNW + 3 + RETW))
    printf " %s\n" "$(printf '%*s' "$total" '' | tr ' ' '-')"
  }

  echo
  printf " ${BOLD}${MAGENTA}CALOFITTING QA Summary  (${TAG})${RESET}\n"
  draw_hr
  printf "  %-*s | %*s | %*s | %*s\n" \
    $LBLW "Stage" $EVW ".evt Events" $RUNW "Runs" $RETW "(retained)"
  printf "  %-${LBLW}s-+-%-${EVW}s-+-%-${RUNW}s-+-%-${RETW}s\n" \
    "$(printf '%*s' "$LBLW" '' | tr ' ' '-')" \
    "$(printf '%*s' "$EVW"  '' | tr ' ' '-')" \
    "$(printf '%*s' "$RUNW" '' | tr ' ' '-')" \
    "$(printf '%*s' "$RETW" '' | tr ' ' '-')"

  printf "  %-*.*s | %*s | %*s | %*.2f%%\n" \
    $LBLW $LBLW "1) All (CreateDstList printruns)" \
    $EVW "$s1_ev" $RUNW "$s1_runs" $RETW 100.00

  printf "  %-*.*s | %*s | %*s | %*.2f%%\n" \
    $LBLW $LBLW "2) MAGNET ON (magnet_info=='t')" \
    $EVW "$s2_ev" $RUNW "$s2_runs" $RETW "$(pct "$s2_runs" "$s1_runs")"

  printf "  %-*.*s | %*s | %*s | %*.2f%%\n" \
    $LBLW $LBLW "3) Calo QA GOLDEN (EMC+IHC+OHC)" \
    $EVW "$s3_ev" $RUNW "$s3_runs" $RETW "$(pct "$s3_runs" "$s1_runs")"

  printf "  %-*.*s | %*s | %*s | %*.2f%%\n" \
    $LBLW $LBLW "4) Runtime >= ${MIN_RUNTIME}s" \
    $EVW "$s4_ev" $RUNW "$s4_runs" $RETW "$(pct "$s4_runs" "$s1_runs")"

  printf "  %-*.*s | %*s | %*s | %*.2f%%\n" \
    $LBLW $LBLW "5) GL1 .evt >= ${MIN_GL1_EVT}" \
    $EVW "$s5_ev" $RUNW "$s5_runs" $RETW "$(pct "$s5_runs" "$s1_runs")"

  printf "  %-*.*s | %*s | %*s | %*.2f%%\n" \
    $LBLW $LBLW "6) Photon3+MBD_NS>=1 active (bit 25)" \
    $EVW "$s6_ev" $RUNW "$s6_runs" $RETW "$(pct "$s6_runs" "$s1_runs")"

  printf "  %-*.*s | %*s | %*s | %*.2f%%\n" \
    $LBLW $LBLW "7) MBD_NS>=1 active (bit10)" \
    $EVW "$s7_ev" $RUNW "$s7_runs" $RETW "$(pct "$s7_runs" "$s1_runs")"

  _finalStage="8) FINAL GOLDEN (2+3+4+5+6+7)"
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
  printf "    • %-22s : %5d\n" "Bit25 inactive" "$rej_trig25"
  printf "    • %-22s : %5d\n" "Bit10 inactive" "$rej_trig10"
  echo

  # ==================== Trigger summary (FINAL GOLDEN ONLY; no scaledown) ====================
  # Aggregates over RUNS_GOLDEN_FINAL using gl1_scalers only:
  #  - ONruns : number of golden runs where a trigger recorded scaled>0
  #  - ΣScaled: total 'scaled' counts over golden runs
  # Sorted by ΣScaled (least → most).

  if (( gold_runs > 0 )); then
    say "Building trigger summary for FINAL GOLDEN runs (${gold_runs})…"

    # Accumulators (keyed by trigger bit/index)
    declare -A TRIG_NAME
    declare -A TRIG_ONRUNS
    declare -A TRIG_SUMSCALED

    gidx=0
    for run in "${RUNS_GOLDEN_FINAL[@]}"; do
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
      (( gidx % 200 == 0 )) && say "  Trigger scan progress: $gidx / $gold_runs"
    done

    # Prepare rows: bit | trigger | ONruns | ΣScaled ; sort by ΣScaled asc
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
      echo " ${BOLD}${GREEN}Trigger totals over FINAL GOLDEN runs (least → most events)${RESET}"
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
      warn "No trigger activity found within FINAL GOLDEN runs."
    fi
  fi
  # ================== End trigger summary (FINAL GOLDEN ONLY; no scaledown) ==================

  # ---------- Write manifest ----------
  GOLDEN_TXT="${OUT_DIR}/run3goldenruns.txt"
  printf '%s\n' "${RUNS_GOLDEN_FINAL[@]}" > "$GOLDEN_TXT"
  good "Golden run manifest → $GOLDEN_TXT  (${gold_runs} runs)"

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

  exit 0
fi

# ============================================================================
# Baseline behavior (unchanged): build per-run lists from LIST_FILE by calling
# CreateDstList.pl per-run and checking the resulting output file.
# ============================================================================
[[ -f "$LIST_FILE" ]]      || { echo "[ERROR] Run list not found: $LIST_FILE"; exit 1; }
[[ -s "$LIST_FILE" ]]      || { echo "[ERROR] Run list is empty: $LIST_FILE"; exit 1; }

mkdir -p "$OUT_DIR"
cd "$OUT_DIR" || exit 1

echo "Tag: $TAG"
echo "Type: $TYPE"
echo "Input: $LIST_FILE"
echo "Output dir: $OUT_DIR"
echo "Tool: $CREATE_DST_TOOL"
echo "----------------------------------------"

total=0
made=0

# Always start fresh: wipe previous output list files in this directory
rm -f dst_*.list DST_*.list

while IFS= read -r raw; do
  # strip comments and whitespace
  run="${raw%%#*}"
  run="$(echo -n "$run" | tr -d '[:space:]')"
  [[ -z "$run" ]] && continue

  ((total++))
  # normalize number, keep a zero-padded version for filename check
  runnum=$(printf "%d" "$run" 2>/dev/null || echo "$run")
  pad=$(printf "%08d" "$runnum")

  # Call the tool exactly as when it worked, just changing the DST type
  perl "$CREATE_DST_TOOL" --tag "$TAG" --run "$runnum" "$TYPE"
  rc=$?

  # Accept success even if tool is noisy; verify file existence (both case variants)
  if [[ -s "dst_jetcalo-$pad.list" || -s "DST_JETCALO-$pad.list" ]]; then
    ((made++))
  else
    echo "[WARN] No list file found (or file empty) for run $runnum (exit=$rc)"
  fi
done < "$LIST_FILE"

echo "----------------------------------------"
echo "Requested runs: $total"
echo "List files made: $made"
echo "Done. Files are in: $OUT_DIR"
