#!/usr/bin/env bash
set -euo pipefail

RJ_REPO_BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
readonly RJ_REPO_BASE
TRAIN_BASE="${RJ_AUAU_TIGHT_LOGREG_TRAIN_BASE:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining}"
MODEL_BASE="${RJ_AUAU_LOGREG_MODEL_BASE:-${RJ_REPO_BASE}/logreg_models}"
TRAIN_SCRIPT="${RJ_AUAU_TIGHT_LOGREG_TRAIN_SCRIPT:-${RJ_REPO_BASE}/scripts/train_auau_photon_logreg.py}"
VALIDATE_SCRIPT="${RJ_AUAU_TIGHT_LOGREG_VALIDATE_SCRIPT:-${RJ_REPO_BASE}/scripts/validate_auau_tight_logreg_on_sim.py}"
ML_PYTHON="${RJ_ML_PYTHON:-python3}"
NOTIFY_EMAILS="${RJ_NOTIFY_EMAILS:-just0131@gmail.com}"

ts() { date +%Y%m%d_%H%M%S; }
say() { printf '\033[1;36m[auauTightLogReg]\033[0m %s\n' "$*"; }
err() { printf '\033[1;31m[auauTightLogReg][ERR]\033[0m %s\n' "$*" >&2; }
die() { err "$*"; exit 2; }

setup_sphenix_stack_env() {
  export USER="${USER:-$(id -u -n)}"
  export LOGNAME="${LOGNAME:-$USER}"
  export HOME="/sphenix/u/${LOGNAME}"
  local myinstall="/sphenix/u/${USER}/thesisAnalysis/install"
  local myinstall_auau="/sphenix/u/${USER}/thesisAnalysis_auau/install"
  set +u
  # shellcheck disable=SC1091
  source /opt/sphenix/core/bin/sphenix_setup.sh -n
  [[ -d "$myinstall" ]] && source /opt/sphenix/core/bin/setup_local.sh "$myinstall" || true
  [[ -d "$myinstall_auau" ]] && source /opt/sphenix/core/bin/setup_local.sh "$myinstall_auau" || true
  set -u
}

setup_ml_python_env() {
  setup_sphenix_stack_env
  local ml_python_prefix="" ml_python_real="" ml_python_real_prefix="" d ld_joined=""
  if [[ "$ML_PYTHON" == */* ]]; then
    ml_python_prefix="$(cd "$(dirname "$ML_PYTHON")/.." && pwd -P 2>/dev/null || true)"
    ml_python_real="$(readlink -f "$ML_PYTHON" 2>/dev/null || true)"
    [[ -n "$ml_python_real" ]] && ml_python_real_prefix="$(cd "$(dirname "$ml_python_real")/.." && pwd -P 2>/dev/null || true)"
  fi
  for d in "$ml_python_prefix/lib" "$ml_python_prefix/lib64" "$ml_python_real_prefix/lib" "$ml_python_real_prefix/lib64"; do
    [[ -n "$d" && -d "$d" ]] || continue
    case ":$ld_joined:" in *":$d:"*) ;; *) ld_joined="${ld_joined:+$ld_joined:}$d" ;; esac
  done
  [[ -n "$ld_joined" ]] && export LD_LIBRARY_PATH="$ld_joined:${LD_LIBRARY_PATH:-}"
  unset PYTHONHOME
}

guard_generated_path() {
  local label="$1" path="$2"
  [[ -n "$path" ]] || die "$label resolved to an empty path"
  case "$path" in /cvmfs/*|/opt/sphenix/*|/usr/*|/bin/*|/lib/*|/lib64/*|/etc/*) die "$label resolved to protected path: $path" ;; esac
}

make_root_manifest() {
  local search_root="$1" out="$2"
  mkdir -p "$(dirname "$out")"
  find "$search_root" -type f -name '*.root' | sort -V > "$out"
}

send_summary_email() {
  local subject="$1" summary="$2"
  if command -v mail >/dev/null 2>&1 && [[ -s "$summary" ]]; then
    mail -s "$subject" "$NOTIFY_EMAILS" < "$summary" || true
  fi
}

usage() {
  cat <<'EOF'
Usage:
  ./scripts/auau_tight_logreg_pipeline.sh trainFromExtraction SOURCE=/path [MODEL_DIR=/path]
  ./scripts/auau_tight_logreg_pipeline.sh smokeTrainFromExtraction SOURCE=/path [MODEL_DIR=/path]
  ./scripts/auau_tight_logreg_pipeline.sh applyCheck MODEL_DIR=/path
  ./scripts/auau_tight_logreg_pipeline.sh validateOnSim SOURCE=/path MODEL_DIR=/path [OUTDIR=/path]
  ./scripts/auau_tight_logreg_pipeline.sh validateOnSimCondor SOURCE=/path MODEL_DIR=/path [groupSize N]
  ./scripts/auau_tight_logreg_pipeline.sh deriveWorkingPointsFromValidation VALIDATION=/path [TARGET=0.80]
  ./scripts/auau_tight_logreg_pipeline.sh generateWorkingPointConfig TEMPLATE=/path WORKING_POINTS=/path ARTIFACT=/path OUT=/path
EOF
}

train_from_extraction() {
  local source="" model_dir="" products="${RJ_AUAU_TIGHT_LOGREG_PRODUCTS:-all}"
  local pt_range="${RJ_AUAU_LOGREG_TRAIN_PT_RANGE:-15:35}"
  local max_files="${RJ_AUAU_LOGREG_TRAIN_MAX_FILES_PER_SAMPLE:-0}"
  local max_rows="${RJ_AUAU_LOGREG_TRAIN_MAX_ROWS:-0}"
  local max_rows_per_class="${RJ_AUAU_LOGREG_TRAIN_MAX_ROWS_PER_CLASS:-0}"
  local max_rows_per_pt_bin_class="${RJ_AUAU_LOGREG_TRAIN_MAX_ROWS_PER_PT_BIN_CLASS:-160000}"
  local epochs="${RJ_AUAU_LOGREG_TRAIN_EPOCHS:-60}"
  local patience="${RJ_AUAU_LOGREG_TRAIN_PATIENCE:-14}"
  local batch_size="${RJ_AUAU_LOGREG_TRAIN_BATCH_SIZE:-262144}"
  local lr="${RJ_AUAU_LOGREG_TRAIN_LEARNING_RATE:-0.018}"
  local l2="${RJ_AUAU_LOGREG_TRAIN_L2:-0.002}"
  local progress="${RJ_AUAU_LOGREG_TRAIN_PROGRESS_EVERY:-2}"
  local tok
  for tok in "$@"; do
    case "$tok" in
      SOURCE=*) source="${tok#SOURCE=}" ;;
      MODEL_DIR=*|MODELDIR=*|model_dir=*) model_dir="${tok#*=}" ;;
      PRODUCTS=*|products=*) products="${tok#*=}" ;;
      PT_RANGE=*|ptRange=*|pt_range=*) pt_range="${tok#*=}" ;;
      MAX_FILES_PER_SAMPLE=*|maxFilesPerSample=*) max_files="${tok#*=}" ;;
      MAX_ROWS=*|maxRows=*) max_rows="${tok#*=}" ;;
      MAX_ROWS_PER_CLASS=*|maxRowsPerClass=*) max_rows_per_class="${tok#*=}" ;;
      MAX_ROWS_PER_PT_BIN_CLASS=*|maxRowsPerPtBinClass=*) max_rows_per_pt_bin_class="${tok#*=}" ;;
      EPOCHS=*|epochs=*) epochs="${tok#*=}" ;;
      PATIENCE=*|patience=*) patience="${tok#*=}" ;;
      BATCH_SIZE=*|batchSize=*) batch_size="${tok#*=}" ;;
      LEARNING_RATE=*|learningRate=*) lr="${tok#*=}" ;;
      L2=*|l2=*) l2="${tok#*=}" ;;
    esac
  done
  [[ -n "$source" && -d "$source" ]] || die "trainFromExtraction requires SOURCE=/path"
  [[ -s "$TRAIN_SCRIPT" ]] || die "Missing trainer: $TRAIN_SCRIPT"
  local stamp="${RJ_AUAU_TIGHT_LOGREG_STAMP:-$(ts)}"
  model_dir="${model_dir:-${MODEL_BASE}/tight_logreg_${stamp}}"
  guard_generated_path "model_dir" "$model_dir"
  mkdir -p "$model_dir"
  setup_ml_python_env
  say "Training AuAu tight-logreg products"
  say "  source=$source"
  say "  model_dir=$model_dir"
  say "  products=$products pt_range=$pt_range"
  "$ML_PYTHON" "$TRAIN_SCRIPT" \
    --source "$source" \
    --outdir "$model_dir" \
    --products "$products" \
    --pt-range "$pt_range" \
    --max-files-per-sample "$max_files" \
    --max-rows "$max_rows" \
    --max-rows-per-class "$max_rows_per_class" \
    --max-rows-per-pt-bin-class "$max_rows_per_pt_bin_class" \
    --epochs "$epochs" \
    --patience "$patience" \
    --batch-size "$batch_size" \
    --learning-rate "$lr" \
    --l2 "$l2" \
    --progress-every "$progress"
}

smoke_train_from_extraction() {
  export RJ_AUAU_TIGHT_LOGREG_PRODUCTS="${RJ_AUAU_TIGHT_LOGREG_PRODUCTS:-all}"
  export RJ_AUAU_LOGREG_TRAIN_MAX_FILES_PER_SAMPLE="${RJ_AUAU_LOGREG_TRAIN_MAX_FILES_PER_SAMPLE:-60}"
  export RJ_AUAU_LOGREG_TRAIN_MAX_ROWS="${RJ_AUAU_LOGREG_TRAIN_MAX_ROWS:-120000}"
  export RJ_AUAU_LOGREG_TRAIN_MAX_ROWS_PER_PT_BIN_CLASS="${RJ_AUAU_LOGREG_TRAIN_MAX_ROWS_PER_PT_BIN_CLASS:-8000}"
  export RJ_AUAU_LOGREG_TRAIN_EPOCHS="${RJ_AUAU_LOGREG_TRAIN_EPOCHS:-18}"
  export RJ_AUAU_LOGREG_TRAIN_PATIENCE="${RJ_AUAU_LOGREG_TRAIN_PATIENCE:-8}"
  export RJ_AUAU_LOGREG_TRAIN_BATCH_SIZE="${RJ_AUAU_LOGREG_TRAIN_BATCH_SIZE:-32768}"
  train_from_extraction "$@"
}

apply_check() {
  local model_dir="" tok
  for tok in "$@"; do
    case "$tok" in MODEL_DIR=*|MODELDIR=*|model_dir=*) model_dir="${tok#*=}" ;; esac
  done
  [[ -n "$model_dir" && -d "$model_dir" ]] || die "applyCheck requires MODEL_DIR=/path"
  [[ -s "${model_dir}/model_registry.json" ]] || die "Missing model registry: ${model_dir}/model_registry.json"
  setup_ml_python_env
  "$ML_PYTHON" - "$model_dir" <<'PY'
import json, math, sys
from pathlib import Path
sys.path.insert(0, str(Path("scripts").resolve()))
from train_auau_photon_logreg import load_artifact, payload_to_model, score_model
import pandas as pd
model_dir = Path(sys.argv[1])
reg = json.loads((model_dir / "model_registry.json").read_text())
bad = []
n = 0
for row in reg.get("models", []):
    path = Path(row["output_json"])
    if not path.is_absolute():
        path = model_dir / path.name
    art = load_artifact(path)
    features = art.get("features", [])
    frame = pd.DataFrame({f: [0.1] for f in features})
    if "cluster_Et" in frame: frame["cluster_Et"] = [20.0]
    if "centrality" in frame: frame["centrality"] = [30.0]
    if art.get("routes"):
        y = 0.5
    else:
        y = float(score_model(payload_to_model(art["model"]), frame)[0])
    if not math.isfinite(y) or y < 0.0 or y > 1.0:
        bad.append(f"{path}: bad smoke score {y}")
    n += 1
if bad:
    raise SystemExit("LogReg applyCheck failed:\n  " + "\n  ".join(bad))
print(f"[OK] applyCheck loaded {n} logreg JSON artifacts")
PY
}

validate_on_sim() {
  local source="" model_dir="" outdir="" score_max="${RJ_AUAU_TIGHT_LOGREG_VALIDATE_SCORE_MAX_ROWS:-300000}" tok
  for tok in "$@"; do
    case "$tok" in
      SOURCE=*) source="${tok#SOURCE=}" ;;
      MODEL_DIR=*|MODELDIR=*|model_dir=*) model_dir="${tok#*=}" ;;
      OUTDIR=*|outdir=*) outdir="${tok#*=}" ;;
      SCORE_MAX_ROWS=*|scoreMaxRows=*) score_max="${tok#*=}" ;;
    esac
  done
  [[ -n "$source" && -d "$source" ]] || die "validateOnSim requires SOURCE=/path"
  [[ -n "$model_dir" && -d "$model_dir" ]] || die "validateOnSim requires MODEL_DIR=/path"
  setup_ml_python_env
  local args=( "$VALIDATE_SCRIPT" --source "$source" --model-dir "$model_dir" --score-max-rows "$score_max" )
  [[ -n "$outdir" ]] && args+=( --outdir "$outdir" )
  "$ML_PYTHON" "${args[@]}"
}

validate_on_sim_condor() {
  local source="" model_dir="" outdir="" group_size="${RJ_AUAU_TIGHT_LOGREG_VALIDATE_GROUP_SIZE:-100}"
  local total_score_max="${RJ_AUAU_TIGHT_LOGREG_VALIDATE_TOTAL_SCORE_MAX_ROWS:-600000}"
  local reqmem="${RJ_AUAU_TIGHT_LOGREG_VALIDATE_REQUEST_MEMORY:-3500MB}" tok
  for tok in "$@"; do
    case "$tok" in
      SOURCE=*) source="${tok#SOURCE=}" ;;
      MODEL_DIR=*|MODELDIR=*|model_dir=*) model_dir="${tok#*=}" ;;
      OUTDIR=*|outdir=*) outdir="${tok#*=}" ;;
      groupSize=*) group_size="${tok#groupSize=}" ;;
      SCORE_MAX_ROWS=*|scoreMaxRows=*) total_score_max="${tok#*=}" ;;
    esac
  done
  while (($#)); do
    case "$1" in
      groupSize) group_size="${2:?missing value after groupSize}"; shift 2 ;;
      scoreMaxRows|totalScoreMaxRows) total_score_max="${2:?missing value after $1}"; shift 2 ;;
      *) shift ;;
    esac
  done
  [[ -n "$source" && -d "$source" ]] || die "validateOnSimCondor requires SOURCE=/path"
  [[ -n "$model_dir" && -d "$model_dir" ]] || die "validateOnSimCondor requires MODEL_DIR=/path"
  local stamp="${RJ_AUAU_TIGHT_LOGREG_VALIDATE_STAMP:-$(ts)}"
  local report_root="${outdir:-${source}/reports/logreg_model_validation_condor_${stamp}}"
  local sub_root="${RJ_REPO_BASE}/condor_sub/auauTightLogRegValidate_${stamp}"
  local shard_dir="${sub_root}/shards" cache_dir="${report_root}/score_caches"
  guard_generated_path "validation report root" "$report_root"
  guard_generated_path "validation submit root" "$sub_root"
  mkdir -p "$report_root" "$sub_root" "$shard_dir" "$cache_dir"
  local root_manifest="${source}/manifests/training_roots.list"
  local search_root="$source"
  [[ -d "${source}/extraction" ]] && search_root="${source}/extraction"
  make_root_manifest "$search_root" "$root_manifest"
  local nroots
  nroots="$(wc -l < "$root_manifest" | tr -d ' ')"
  [[ "$nroots" != "0" ]] || die "No ROOT files available for validation"
  split -l "$group_size" -d -a 5 "$root_manifest" "${shard_dir}/roots_"
  local shard_count
  shard_count="$(find "$shard_dir" -maxdepth 1 -type f -name 'roots_*' | wc -l | tr -d ' ')"
  local score_max_per_shard=0
  if [[ "$total_score_max" =~ ^[0-9]+$ && "$total_score_max" -gt 0 ]]; then
    score_max_per_shard=$(( (total_score_max + shard_count - 1) / shard_count ))
  fi
  local worker="${sub_root}/validate_worker.sh"
  cat > "$worker" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export USER="\${USER:-\$(id -u -n)}"
export LOGNAME="\${LOGNAME:-\$USER}"
export HOME="/sphenix/u/\${LOGNAME}"
MYINSTALL="/sphenix/u/\${USER}/thesisAnalysis/install"
MYINSTALL_AUAU="/sphenix/u/\${USER}/thesisAnalysis_auau/install"
set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n
[[ -d "\$MYINSTALL" ]] && source /opt/sphenix/core/bin/setup_local.sh "\$MYINSTALL" || true
[[ -d "\$MYINSTALL_AUAU" ]] && source /opt/sphenix/core/bin/setup_local.sh "\$MYINSTALL_AUAU" || true
set -u
ml_python="${ML_PYTHON}"
ml_python_prefix="\$(cd "\$(dirname "\$ml_python")/.." && pwd -P 2>/dev/null || true)"
ml_python_real="\$(readlink -f "\$ml_python" 2>/dev/null || true)"
ml_python_real_prefix=""
[[ -n "\$ml_python_real" ]] && ml_python_real_prefix="\$(cd "\$(dirname "\$ml_python_real")/.." && pwd -P 2>/dev/null || true)"
ld_joined=""
for d in "\$ml_python_prefix/lib" "\$ml_python_prefix/lib64" "\$ml_python_real_prefix/lib" "\$ml_python_real_prefix/lib64"; do
  [[ -n "\$d" && -d "\$d" ]] || continue
  case ":\$ld_joined:" in *":\$d:"*) ;; *) ld_joined="\${ld_joined:+\$ld_joined:}\$d" ;; esac
done
[[ -n "\$ld_joined" ]] && export LD_LIBRARY_PATH="\$ld_joined:\${LD_LIBRARY_PATH:-}"
unset PYTHONHOME
cd "${RJ_REPO_BASE}"
"\$ml_python" "${VALIDATE_SCRIPT}" --source "${source}" --model-dir "${model_dir}" --manifest "\${1:?manifest}" --outdir "\${2:?outdir}" --write-score-cache "\${3:?cache}" --score-max-rows "\${4:?score_max}" --no-plots
EOF
  chmod +x "$worker"
  local dag="${sub_root}/auau_tight_logreg_validateOnSimCondor.dag"
  : > "$dag"
  local idx=0 shard
  for shard in "${shard_dir}/roots_"*; do
    [[ -s "$shard" ]] || continue
    idx=$((idx + 1))
    local label sub out cache
    label="$(printf '%05d' "$idx")"
    out="${report_root}/shards/shard_${label}"
    cache="${cache_dir}/score_cache_${label}.npz"
    mkdir -p "$out"
    sub="${sub_root}/validate_${label}.sub"
    cat > "$sub" <<EOF
universe = vanilla
executable = ${worker}
arguments = ${shard} ${out} ${cache} ${score_max_per_shard}
output = ${sub_root}/validate_${label}.out
error = ${sub_root}/validate_${label}.err
log = ${sub_root}/validate_${label}.log
request_memory = ${reqmem}
notification = Never
queue
EOF
    echo "JOB VALIDATE_${label} ${sub}" >> "$dag"
    echo "RETRY VALIDATE_${label} 1" >> "$dag"
  done
  local merge="${sub_root}/validate_merge.sh"
  cat > "$merge" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export USER="\${USER:-\$(id -u -n)}"
export LOGNAME="\${LOGNAME:-\$USER}"
export HOME="/sphenix/u/\${LOGNAME}"
set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n
set -u
ml_python="${ML_PYTHON}"
ml_python_prefix="\$(cd "\$(dirname "\$ml_python")/.." && pwd -P 2>/dev/null || true)"
ml_python_real="\$(readlink -f "\$ml_python" 2>/dev/null || true)"
ml_python_real_prefix=""
[[ -n "\$ml_python_real" ]] && ml_python_real_prefix="\$(cd "\$(dirname "\$ml_python_real")/.." && pwd -P 2>/dev/null || true)"
ld_joined=""
for d in "\$ml_python_prefix/lib" "\$ml_python_prefix/lib64" "\$ml_python_real_prefix/lib" "\$ml_python_real_prefix/lib64"; do
  [[ -n "\$d" && -d "\$d" ]] || continue
  case ":\$ld_joined:" in *":\$d:"*) ;; *) ld_joined="\${ld_joined:+\$ld_joined:}\$d" ;; esac
done
[[ -n "\$ld_joined" ]] && export LD_LIBRARY_PATH="\$ld_joined:\${LD_LIBRARY_PATH:-}"
unset PYTHONHOME
cd "${RJ_REPO_BASE}"
cache_manifest="${report_root}/score_caches.list"
find "${cache_dir}" -type f -name 'score_cache_*.npz' | sort -V > "\$cache_manifest" || true
expected=${idx}
found=\$(wc -l < "\$cache_manifest" | tr -d ' ')
summary="${report_root}/validation_summary.txt"
if [[ "\$found" != "\$expected" || "\$found" == "0" ]]; then
  {
    echo "RECOILJETS_AUAU_TIGHT_LOGREG_SIM_VALIDATION_V1"
    echo "status=CHECK"
    echo "source=${source}"
    echo "model_dir=${model_dir}"
    echo "report_dir=${report_root}"
    echo "expected_score_caches=\$expected"
    echo "found_score_caches=\$found"
    echo "notes=missing score cache shards"
  } > "\$summary"
else
  "\$ml_python" "${VALIDATE_SCRIPT}" --source "${source}" --model-dir "${model_dir}" --merge-score-caches "\$cache_manifest" --outdir "${report_root}"
fi
status=\$(awk -F= '/^status=/ {print \$2; exit}' "\$summary" 2>/dev/null || echo CHECK)
if command -v mail >/dev/null 2>&1 && [[ -s "\$summary" ]]; then
  mail -s "[RecoilJets][auauTightLogReg_validateOnSimCondor][\${status:-CHECK}]" "${NOTIFY_EMAILS}" < "\$summary" || true
fi
cat "\$summary"
[[ "\${status:-CHECK}" == "READY" ]] || exit 3
EOF
  chmod +x "$merge"
  local merge_sub="${sub_root}/validate_merge.sub"
  cat > "$merge_sub" <<EOF
universe = scheduler
executable = ${merge}
output = ${sub_root}/validate_merge.out
error = ${sub_root}/validate_merge.err
log = ${sub_root}/validate_merge.log
notification = Never
queue
EOF
  echo "FINAL MERGE ${merge_sub}" >> "$dag"
  say "Condor validation DAG: $dag"
  say "source=$source model_dir=$model_dir report=$report_root"
  say "root files=$nroots shards=$idx groupSize=$group_size scoreMaxPerShard=$score_max_per_shard request_memory=$reqmem"
  if [[ "${RJ_DAG_DRYRUN:-0}" == "1" ]]; then
    echo "RECOILJETS_AUAU_TIGHT_LOGREG_VALIDATE_DRYRUN_V1"
    echo "dag=${dag}"
    echo "report_root=${report_root}"
    return 0
  fi
  condor_submit_dag "$dag"
}

derive_working_points_from_validation() {
  local validation="" target="0.80" tok
  for tok in "$@"; do
    case "$tok" in VALIDATION=*) validation="${tok#VALIDATION=}" ;; TARGET=*) target="${tok#TARGET=}" ;; esac
  done
  [[ -n "$validation" && -d "$validation" ]] || die "deriveWorkingPointsFromValidation requires VALIDATION=/path"
  setup_ml_python_env
  "$ML_PYTHON" "$VALIDATE_SCRIPT" --derive-working-points-from-report "$validation" --target-signal-efficiency "$target"
}

generate_working_point_config() {
  local template="" wp="" artifact="" out="" product="" tok
  for tok in "$@"; do
    case "$tok" in
      TEMPLATE=*) template="${tok#TEMPLATE=}" ;;
      WORKING_POINTS=*) wp="${tok#WORKING_POINTS=}" ;;
      ARTIFACT=*) artifact="${tok#ARTIFACT=}" ;;
      OUT=*) out="${tok#OUT=}" ;;
      PRODUCT=*) product="${tok#PRODUCT=}" ;;
    esac
  done
  [[ -n "$template" ]] || die "generateWorkingPointConfig requires TEMPLATE=/path"
  [[ -n "$wp" ]] || die "generateWorkingPointConfig requires WORKING_POINTS=/path"
  [[ -n "$artifact" ]] || die "generateWorkingPointConfig requires ARTIFACT=/path"
  [[ -n "$out" ]] || die "generateWorkingPointConfig requires OUT=/path"
  local args=( "${RJ_REPO_BASE}/scripts/make_auau_logreg_target_wp_config.py" --template "$template" --working-points "$wp" --artifact "$artifact" --out "$out" )
  [[ -n "$product" ]] && args+=( --product "$product" )
  python3 "${args[@]}"
}

main() {
  local mode="${1:-}"
  [[ -n "$mode" ]] || { usage; exit 2; }
  shift || true
  case "$mode" in
    -h|--help|help) usage ;;
    trainFromExtraction) train_from_extraction "$@" ;;
    smokeTrainFromExtraction) smoke_train_from_extraction "$@" ;;
    applyCheck) apply_check "$@" ;;
    validateOnSim) validate_on_sim "$@" ;;
    validateOnSimCondor|condorValidateOnSim|validateSimCondor) validate_on_sim_condor "$@" ;;
    deriveWorkingPointsFromValidation) derive_working_points_from_validation "$@" ;;
    generateWorkingPointConfig) generate_working_point_config "$@" ;;
    *) die "Unknown mode: $mode" ;;
  esac
}

main "$@"
