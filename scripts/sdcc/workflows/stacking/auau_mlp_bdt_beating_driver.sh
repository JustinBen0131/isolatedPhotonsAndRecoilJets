#!/usr/bin/env bash
set -euo pipefail

REPO_BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
readonly REPO_BASE

ML_PYTHON="${RJ_ML_PYTHON:-/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python}"
SOURCE="${RJ_AUAU_MLP_SOURCE:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049}"
MODEL_BASE="${RJ_AUAU_MLP_MODEL_BASE:-/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models}"
SWEEP_DIR="${RJ_AUAU_MLP_SWEEP_DIR:-${MODEL_BASE}/tight_mlp_highpt_sweep_20260512_182344}"
VALIDATION_CACHE="${RJ_AUAU_MLP_VALIDATION_CACHE:-${SOURCE}/reports/mlp_model_validation_condor_deep_primary_ratios_nostat_fullval_20260512_145449/score_caches.list}"
V2_TEMPLATE="${RJ_AUAU_MLP_V2_TEMPLATE:-macros/analysis_config_auau_mlp_v2_validation.yaml}"

say() { printf '\033[1;36m[auauMLPBDTBeat]\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[auauMLPBDTBeat][WARN]\033[0m %s\n' "$*" >&2; }
die() { printf '\033[1;31m[auauMLPBDTBeat][ERR]\033[0m %s\n' "$*" >&2; exit 2; }

stamp() { date +%Y%m%d_%H%M%S; }

require_run() {
  if [[ "${RJ_DO_RUN:-0}" != "1" ]]; then
    warn "Dry-run guard: set RJ_DO_RUN=1 to start tmux jobs, write manifests, rescore, or submit validation."
    return 1
  fi
}

require_condor() {
  if [[ "${RJ_ALLOW_CONDOR:-0}" != "1" ]]; then
    warn "Condor guard: set RJ_ALLOW_CONDOR=1 as well as RJ_DO_RUN=1 before submitting validation."
    return 1
  fi
}

usage() {
  cat <<'EOF'
Usage:
  ./scripts/auau_mlp_bdt_beating_driver.sh status
  ./scripts/auau_mlp_bdt_beating_driver.sh recoverGlobal5Tmux
  ./scripts/auau_mlp_bdt_beating_driver.sh writeHighPtManifest
  ./scripts/auau_mlp_bdt_beating_driver.sh rescoreHighPt
  ./scripts/auau_mlp_bdt_beating_driver.sh trainV2Tmux
  ./scripts/auau_mlp_bdt_beating_driver.sh validateV2Condor MODEL_DIR=/path
  ./scripts/auau_mlp_bdt_beating_driver.sh watchAndValidateV2 MODEL_DIR=/path

Safe-by-default orchestration for the BDT-beating AuAu tight-MLP lane.
Mutating modes require RJ_DO_RUN=1. Condor submission additionally requires
RJ_ALLOW_CONDOR=1. Read-only status is always allowed.
EOF
}

latest_dir() {
  local pattern="$1"
  compgen -G "$pattern" | sort -V | tail -1 || true
}

artifact_exists() {
  [[ -s "$1" ]]
}

print_common() {
  say "host=$(hostname -f 2>/dev/null || hostname)"
  say "repo=${REPO_BASE}"
  say "source=${SOURCE}"
  say "model_base=${MODEL_BASE}"
  say "sweep_dir=${SWEEP_DIR}"
  say "validation_cache=${VALIDATION_CACHE}"
}

status() {
  cd "$REPO_BASE"
  print_common
  say "tmux sessions matching MLP work:"
  tmux ls 2>/dev/null | grep -E 'mlp|MLP|kitchen|highpt|wp80' || true
  say "recent model registries:"
  find "$MODEL_BASE" -maxdepth 4 -name model_registry.json -printf '%TY-%Tm-%Td %TH:%TM %p\n' 2>/dev/null | sort | tail -20 || true
  say "high-pT rescore rank tables:"
  find "$SWEEP_DIR" -maxdepth 3 -name validation_rank_table.csv -printf '%TY-%Tm-%Td %TH:%TM %p\n' 2>/dev/null | sort | tail -10 || true
  say "kitchen-sink log tails:"
  for log in \
    "$MODEL_BASE"/train_kitchensink_tmux_*.log \
    "$MODEL_BASE"/train_iso_kitchensink_tmux_*.log \
    "$MODEL_BASE"/train_highpt_distilled_kitchen_v2_*.log \
    "$MODEL_BASE"/train_global5_broadReach_recover_*.log; do
    [[ -s "$log" ]] || continue
    echo "---- $log"
    tail -12 "$log" || true
  done
  say "Condor queue summary:"
  condor_q -nobatch "${USER:-patsfan753}" 2>/dev/null | tail -50 || true
  say "held jobs:"
  condor_q -hold "${USER:-patsfan753}" 2>/dev/null | tail -40 || true
}

recover_global5_tmux() {
  cd "$REPO_BASE"
  local now model_dir log runner session
  now="$(stamp)"
  model_dir="${GLOBAL5_RECOVERY_DIR:-${SWEEP_DIR}/global5_broadReach_recover_${now}}"
  log="${GLOBAL5_LOG:-${MODEL_BASE}/train_global5_broadReach_recover_${now}.log}"
  runner="/tmp/train_global5_broadReach_recover_${now}.sh"
  session="${GLOBAL5_SESSION:-mlp_global5_recover_${now}}"
  print_common
  say "global5 recovery model_dir=${model_dir}"
  say "global5 recovery log=${log}"
  say "global5 recovery tmux=${session}"
  require_run || return 0
  mkdir -p "$model_dir" "$(dirname "$log")"
  cat > "$runner" <<EOF
#!/usr/bin/env bash
set -euo pipefail
cd "$REPO_BASE"
export RJ_ML_PYTHON="$ML_PYTHON"
./scripts/auau_tight_mlp_pipeline.sh trainHighPtBalancedFromExtraction \\
  SOURCE="$SOURCE" \\
  MODEL_DIR="$model_dir" \\
  PRODUCTS=primary-ratios \\
  PT_RANGE=5:35 \\
  TRAIN_PT_BINS=5,10,15,20,25,35 \\
  PT_BIN_WEIGHT_SPEC=5:10:0.10,10:15:0.15,15:20:0.20,20:25:0.25,25:35:0.30 \\
  HIGHPT_SELECTION_WEIGHTS=5:10:0.03,10:15:0.07,15:20:0.20,20:25:0.30,25:35:0.40 \\
  MAX_ROWS_PER_PT_BIN_CLASS="${GLOBAL5_MAX_ROWS_PER_PT_BIN_CLASS:-90000}" \\
  EPOCHS="${GLOBAL5_EPOCHS:-220}" \\
  PATIENCE="${GLOBAL5_PATIENCE:-45}" \\
  RESTARTS="${GLOBAL5_RESTARTS:-1}" \\
  HIDDEN_LAYERS="${GLOBAL5_HIDDEN_LAYERS:-128,64,32}" \\
  HIDDEN_LAYER_GRID="${GLOBAL5_HIDDEN_LAYER_GRID:-160,80,40}"
./scripts/auau_tight_mlp_pipeline.sh applyCheck MODEL_DIR="$model_dir"
echo "DONE_GLOBAL5_BROADREACH_MODEL_DIR=$model_dir"
EOF
  chmod +x "$runner"
  tmux new-session -d -s "$session" "bash '$runner' 2>&1 | tee '$log'"
  say "started tmux session: $session"
  say "tail with: tmux attach -t $session"
  say "detach with: Ctrl-b then d"
}

write_highpt_manifest() {
  cd "$REPO_BASE"
  local manifest global5_dir now
  now="$(stamp)"
  manifest="${HIGHPT_MANIFEST:-${SWEEP_DIR}/sweep_manifest_current_${now}.json}"
  global5_dir="${GLOBAL5_DIR:-}"
  if [[ -z "$global5_dir" ]]; then
    global5_dir="$(latest_dir "${SWEEP_DIR}/global5_broadReach_recover_*")"
  fi
  [[ -n "$global5_dir" ]] || global5_dir="${SWEEP_DIR}/global5_broadReach"
  say "manifest=${manifest}"
  say "global5_dir=${global5_dir}"
  require_run || return 0
  mkdir -p "$(dirname "$manifest")"
  python3 - "$manifest" "$SWEEP_DIR" "$global5_dir" "${CURRENT_PRIMARY_MODEL_DIR:-${MODEL_BASE}/tight_mlp_current}" <<'PY'
import json
import sys
from pathlib import Path

manifest = Path(sys.argv[1])
sweep_dir = Path(sys.argv[2])
global5_dir = Path(sys.argv[3])
current_dir = Path(sys.argv[4])
artifact = "auau_tight_mlp_centInputBase3x3_pt1535.json"
variants = []

def add_single(name, model_dir, variant="auauCentInputBase3x3MLP", filename=artifact):
    path = Path(model_dir) / filename
    if path.is_file():
        variants.append({"name": name, "kind": "single", "variant": variant, "artifact": str(path)})
    else:
        print(f"[auauMLPBDTBeat][WARN] missing single artifact for {name}: {path}", file=sys.stderr)

def add_routed(name, axis, route_specs):
    routes = []
    missing = []
    for lo, hi, label, route_dir in route_specs:
        path = Path(route_dir) / artifact
        if path.is_file():
            routes.append({"lo": float(lo), "hi": float(hi), "label": label, "artifact": str(path)})
        else:
            missing.append(str(path))
    if missing:
        print(f"[auauMLPBDTBeat][WARN] skipping routed {name}; missing " + "; ".join(missing), file=sys.stderr)
        return
    variants.append({"name": name, "kind": "routed", "variant": name, "axis": axis, "routes": routes})

add_single("current_primary_ratios_mlp", current_dir)
add_single("global15_highPtBalanced", sweep_dir / "global15_highPtBalanced")
add_single("global5_broadReach", global5_dir)
add_routed(
    "pt3_15to35",
    "cluster_Et",
    [
        (15.0, 20.0, "pt_015_020", sweep_dir / "pt3_15to35" / "pt_015_020"),
        (20.0, 25.0, "pt_020_025", sweep_dir / "pt3_15to35" / "pt_020_025"),
        (25.0, 35.0, "pt_025_035", sweep_dir / "pt3_15to35" / "pt_025_035"),
    ],
)
add_routed(
    "cent3_15to35",
    "centrality",
    [
        (0.0, 20.0, "cent_000_020", sweep_dir / "cent3_15to35" / "cent_000_020"),
        (20.0, 50.0, "cent_020_050", sweep_dir / "cent3_15to35" / "cent_020_050"),
        (50.0, 80.0, "cent_050_080", sweep_dir / "cent3_15to35" / "cent_050_080"),
    ],
)
if not variants:
    raise SystemExit("No sweep artifacts found; refusing to write an empty manifest")
payload = {
    "schema": "RJ_AUAU_TIGHT_MLP_SWEEP_V1",
    "sweep_dir": str(sweep_dir),
    "variants": variants,
}
manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(manifest)
PY
  say "wrote manifest: ${manifest}"
}

rescore_highpt() {
  cd "$REPO_BASE"
  local manifest outdir now
  now="$(stamp)"
  manifest="${HIGHPT_MANIFEST:-}"
  if [[ -z "$manifest" ]]; then
    manifest="$(latest_dir "${SWEEP_DIR}/sweep_manifest_current_*.json")"
  fi
  [[ -n "$manifest" && -s "$manifest" ]] || die "Set HIGHPT_MANIFEST=/path or run writeHighPtManifest first"
  outdir="${HIGHPT_RESCORE_OUTDIR:-${SWEEP_DIR}/validation_rescore_current_${now}}"
  say "cache=${VALIDATION_CACHE}"
  say "manifest=${manifest}"
  say "outdir=${outdir}"
  require_run || return 0
  ./scripts/auau_tight_mlp_pipeline.sh rescoreValidationCache \
    CACHE="$VALIDATION_CACHE" \
    SWEEP_MANIFEST="$manifest" \
    OUTDIR="$outdir" \
    PT_BINS="${HIGHPT_RESCORE_PT_BINS:-15,20,25,35}"
}

train_v2_tmux() {
  cd "$REPO_BASE"
  local now model_dir log runner session report_dir
  now="$(stamp)"
  model_dir="${V2_MODEL_DIR:-${MODEL_BASE}/tight_mlp_highpt_distilled_kitchen_v2_${now}}"
  report_dir="${V2_VALIDATION_REPORT_DIR:-${model_dir}/validation_condor_${now}}"
  log="${V2_LOG:-${MODEL_BASE}/train_highpt_distilled_kitchen_v2_${now}.log}"
  runner="/tmp/train_highpt_distilled_kitchen_v2_${now}.sh"
  session="${V2_SESSION:-mlp_highpt_distilled_v2_${now}}"
  print_common
  say "v2 model_dir=${model_dir}"
  say "v2 validation_report=${report_dir}"
  say "v2 log=${log}"
  say "v2 tmux=${session}"
  require_run || return 0
  mkdir -p "$model_dir" "$(dirname "$log")"
  cat > "$runner" <<EOF
#!/usr/bin/env bash
set -euo pipefail
cd "$REPO_BASE"
export RJ_ML_PYTHON="$ML_PYTHON"
./scripts/auau_tight_mlp_pipeline.sh trainHighPtDistilledKitchenV2FromExtraction \\
  SOURCE="$SOURCE" \\
  MODEL_DIR="$model_dir"
./scripts/auau_tight_mlp_pipeline.sh applyCheck MODEL_DIR="$model_dir"
echo "DONE_HIGHPT_DISTILLED_KITCHEN_V2_MODEL_DIR=$model_dir"
if [[ "\${RJ_MLP_V2_AUTO_VALIDATE:-0}" == "1" && "\${RJ_ALLOW_CONDOR:-0}" == "1" ]]; then
  ./scripts/auau_tight_mlp_pipeline.sh validateOnSimCondor SOURCE="$SOURCE" MODEL_DIR="$model_dir" OUTDIR="$report_dir"
else
  echo "V2 validation not auto-submitted; set RJ_MLP_V2_AUTO_VALIDATE=1 RJ_ALLOW_CONDOR=1 before starting this tmux to allow it."
fi
EOF
  chmod +x "$runner"
  tmux new-session -d -s "$session" "bash '$runner' 2>&1 | tee '$log'"
  say "started tmux session: $session"
  say "tail with: tmux attach -t $session"
  say "detach with: Ctrl-b then d"
}

validate_v2_condor() {
  cd "$REPO_BASE"
  local model_dir="" outdir="" tok now
  now="$(stamp)"
  for tok in "$@"; do
    case "$tok" in
      MODEL_DIR=*|model_dir=*) model_dir="${tok#*=}" ;;
      OUTDIR=*|outdir=*) outdir="${tok#*=}" ;;
    esac
  done
  [[ -n "$model_dir" && -d "$model_dir" ]] || die "validateV2Condor requires MODEL_DIR=/path"
  outdir="${outdir:-${model_dir}/validation_condor_${now}}"
  say "v2 validation model_dir=${model_dir}"
  say "v2 validation outdir=${outdir}"
  require_run || return 0
  require_condor || return 0
  ./scripts/auau_tight_mlp_pipeline.sh validateOnSimCondor SOURCE="$SOURCE" MODEL_DIR="$model_dir" OUTDIR="$outdir"
}

watch_and_validate_v2() {
  cd "$REPO_BASE"
  local model_dir="" outdir="" sleep_s tok now
  now="$(stamp)"
  sleep_s="${RJ_MLP_DRIVER_SLEEP_SECONDS:-900}"
  for tok in "$@"; do
    case "$tok" in
      MODEL_DIR=*|model_dir=*) model_dir="${tok#*=}" ;;
      OUTDIR=*|outdir=*) outdir="${tok#*=}" ;;
      SLEEP=*|sleep=*) sleep_s="${tok#*=}" ;;
    esac
  done
  [[ -n "$model_dir" ]] || die "watchAndValidateV2 requires MODEL_DIR=/path"
  outdir="${outdir:-${model_dir}/validation_condor_watch_${now}}"
  say "watching model_dir=${model_dir}"
  say "will submit validation to ${outdir} after model_registry.json appears"
  require_run || return 0
  require_condor || return 0
  while true; do
    if [[ -s "${model_dir}/model_registry.json" ]]; then
      say "registry found; starting validation"
      ./scripts/auau_tight_mlp_pipeline.sh validateOnSimCondor SOURCE="$SOURCE" MODEL_DIR="$model_dir" OUTDIR="$outdir"
      return 0
    fi
    say "registry not ready yet; sleeping ${sleep_s}s"
    sleep "$sleep_s"
  done
}

main() {
  local mode="${1:-}"
  [[ -n "$mode" ]] || { usage; exit 2; }
  shift || true
  case "$mode" in
    -h|--help|help) usage ;;
    status) status "$@" ;;
    recoverGlobal5Tmux|recoverGlobal5) recover_global5_tmux "$@" ;;
    writeHighPtManifest|writeManifest) write_highpt_manifest "$@" ;;
    rescoreHighPt|rescore) rescore_highpt "$@" ;;
    trainV2Tmux|trainHighPtDistilledKitchenV2) train_v2_tmux "$@" ;;
    validateV2Condor|validateV2) validate_v2_condor "$@" ;;
    watchAndValidateV2|watchV2) watch_and_validate_v2 "$@" ;;
    *) usage; die "Unknown mode: $mode" ;;
  esac
}

main "$@"
