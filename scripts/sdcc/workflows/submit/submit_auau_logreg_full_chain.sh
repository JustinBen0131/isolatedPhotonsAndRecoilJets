#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$repo_root"

stamp="${RJ_AUAU_LOGREG_CHAIN_STAMP:-$(date +%Y%m%d_%H%M%S)}"
source_dir="${SOURCE:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049}"
export RJ_ML_PYTHON="${RJ_ML_PYTHON:-/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python}"

model_base="${MODEL_BASE:-/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/logreg_models}"
model_dir="${MODEL_DIR:-${model_base}/tight_logreg_full_chain_${stamp}}"
smoke_dir="${SMOKE_DIR:-${model_base}/tight_logreg_smoke_${stamp}}"
report_dir="${REPORT_DIR:-${source_dir}/reports/logreg_model_validation_condor_${stamp}}"
config_dir="${CONFIG_DIR:-${repo_root}/condor_generated_configs/logreg_targetwp_${stamp}}"
log_dir="${LOG_DIR:-${repo_root}/condor_sub/auauTightLogRegFullChain_${stamp}}"
mkdir -p "$log_dir" "$config_dir"
log_file="${log_dir}/logreg_full_chain_${stamp}.log"
exec > >(tee -a "$log_file") 2>&1

echo "RECOILJETS_AUAU_TIGHT_LOGREG_FULL_CHAIN_V1"
echo "host=$(hostname -f 2>/dev/null || hostname)"
echo "stamp=${stamp}"
echo "source=${source_dir}"
echo "ml_python=${RJ_ML_PYTHON}"
echo "smoke_dir=${smoke_dir}"
echo "model_dir=${model_dir}"
echo "report_dir=${report_dir}"
echo "config_dir=${config_dir}"
echo "log_file=${log_file}"
echo

if [[ ! -d "$source_dir" ]]; then
  echo "[ERROR] missing source dir: ${source_dir}" >&2
  exit 2
fi

echo "== static checks =="
bash -n scripts/auau_tight_logreg_pipeline.sh scripts/submit_auau_logreg_targetwp_pair.sh
"${RJ_ML_PYTHON}" -m py_compile \
  scripts/train_auau_photon_logreg.py \
  scripts/validate_auau_tight_logreg_on_sim.py \
  scripts/make_auau_logreg_target_wp_config.py
"${RJ_ML_PYTHON}" scripts/train_auau_photon_logreg.py --self-test

echo
echo "== src_AuAu rebuild =="
(
  cd "${repo_root}/src_AuAu"
  if type makee >/dev/null 2>&1; then
    makee clean
  else
    make clean
  fi
  makeProject
)

echo
echo "== capped smoke train =="
RJ_AUAU_TIGHT_LOGREG_PRODUCTS="${RJ_AUAU_TIGHT_LOGREG_PRODUCTS:-all}" \
RJ_AUAU_LOGREG_TRAIN_MAX_FILES_PER_SAMPLE="${RJ_AUAU_LOGREG_SMOKE_MAX_FILES_PER_SAMPLE:-50}" \
RJ_AUAU_LOGREG_TRAIN_MAX_ROWS="${RJ_AUAU_LOGREG_SMOKE_MAX_ROWS:-120000}" \
RJ_AUAU_LOGREG_TRAIN_MAX_ROWS_PER_PT_BIN_CLASS="${RJ_AUAU_LOGREG_SMOKE_MAX_ROWS_PER_PT_BIN_CLASS:-6000}" \
RJ_AUAU_LOGREG_TRAIN_EPOCHS="${RJ_AUAU_LOGREG_SMOKE_EPOCHS:-14}" \
RJ_AUAU_LOGREG_TRAIN_PATIENCE="${RJ_AUAU_LOGREG_SMOKE_PATIENCE:-6}" \
RJ_AUAU_LOGREG_TRAIN_BATCH_SIZE="${RJ_AUAU_LOGREG_SMOKE_BATCH_SIZE:-32768}" \
  bash ./scripts/auau_tight_logreg_pipeline.sh smokeTrainFromExtraction SOURCE="${source_dir}" MODEL_DIR="${smoke_dir}"
bash ./scripts/auau_tight_logreg_pipeline.sh applyCheck MODEL_DIR="${smoke_dir}"

echo
echo "== full logistic training =="
RJ_AUAU_TIGHT_LOGREG_PRODUCTS="${RJ_AUAU_TIGHT_LOGREG_PRODUCTS:-all}" \
RJ_AUAU_LOGREG_TRAIN_MAX_FILES_PER_SAMPLE="${RJ_AUAU_LOGREG_TRAIN_MAX_FILES_PER_SAMPLE:-0}" \
RJ_AUAU_LOGREG_TRAIN_MAX_ROWS="${RJ_AUAU_LOGREG_TRAIN_MAX_ROWS:-0}" \
RJ_AUAU_LOGREG_TRAIN_MAX_ROWS_PER_PT_BIN_CLASS="${RJ_AUAU_LOGREG_TRAIN_MAX_ROWS_PER_PT_BIN_CLASS:-160000}" \
RJ_AUAU_LOGREG_TRAIN_EPOCHS="${RJ_AUAU_LOGREG_TRAIN_EPOCHS:-70}" \
RJ_AUAU_LOGREG_TRAIN_PATIENCE="${RJ_AUAU_LOGREG_TRAIN_PATIENCE:-16}" \
RJ_AUAU_LOGREG_TRAIN_BATCH_SIZE="${RJ_AUAU_LOGREG_TRAIN_BATCH_SIZE:-262144}" \
RJ_AUAU_LOGREG_TRAIN_LEARNING_RATE="${RJ_AUAU_LOGREG_TRAIN_LEARNING_RATE:-0.018}" \
RJ_AUAU_LOGREG_TRAIN_L2="${RJ_AUAU_LOGREG_TRAIN_L2:-0.002}" \
  bash ./scripts/auau_tight_logreg_pipeline.sh trainFromExtraction SOURCE="${source_dir}" MODEL_DIR="${model_dir}"
bash ./scripts/auau_tight_logreg_pipeline.sh applyCheck MODEL_DIR="${model_dir}"

echo
echo "== validation DAG submit =="
RJ_AUAU_TIGHT_LOGREG_VALIDATE_STAMP="${stamp}" \
RJ_AUAU_TIGHT_LOGREG_VALIDATE_TOTAL_SCORE_MAX_ROWS="${RJ_AUAU_TIGHT_LOGREG_VALIDATE_TOTAL_SCORE_MAX_ROWS:-0}" \
RJ_AUAU_TIGHT_LOGREG_VALIDATE_GROUP_SIZE="${RJ_AUAU_TIGHT_LOGREG_VALIDATE_GROUP_SIZE:-100}" \
RJ_AUAU_TIGHT_LOGREG_VALIDATE_REQUEST_MEMORY="${RJ_AUAU_TIGHT_LOGREG_VALIDATE_REQUEST_MEMORY:-4000MB}" \
  bash ./scripts/auau_tight_logreg_pipeline.sh validateOnSimCondor SOURCE="${source_dir}" MODEL_DIR="${model_dir}" OUTDIR="${report_dir}"

echo
echo "== validation polling =="
summary="${report_dir}/validation_summary.txt"
poll_seconds="${RJ_AUAU_LOGREG_POLL_SECONDS:-300}"
while true; do
  if [[ -s "$summary" ]]; then
    cat "$summary"
    status="$(awk -F= '/^status=/ {print $2; exit}' "$summary")"
    if [[ "$status" == "READY" ]]; then
      break
    fi
    echo "[ERROR] validation ended with status=${status:-UNKNOWN}" >&2
    exit 3
  fi
  held="$(condor_q -hold -nobatch "${USER:-patsfan753}" 2>/dev/null | awk 'END{print NR+0}')"
  echo "waiting_for_validation_summary report=${report_dir} held_lines=${held} next_poll_seconds=${poll_seconds}"
  sleep "$poll_seconds"
done

echo
echo "== derive/rank WP80 =="
bash ./scripts/auau_tight_logreg_pipeline.sh deriveWorkingPointsFromValidation VALIDATION="${report_dir}" TARGET=0.80
best_product="$(awk -F= '/^best_product=/ {print $2; exit}' "$summary")"
if [[ -z "$best_product" ]]; then
  best_product="$(awk -F, 'NR==2 {print $1; exit}' "${report_dir}/validation_rank_table.csv")"
fi
artifact="$("${RJ_ML_PYTHON}" - "$model_dir" "$best_product" <<'PY'
import json, sys
from pathlib import Path
model_dir = Path(sys.argv[1])
best = sys.argv[2]
reg = json.loads((model_dir / "model_registry.json").read_text())
for row in reg.get("models", []):
    if row.get("product") == best:
        p = Path(row["output_json"])
        print(p if p.is_absolute() else model_dir / p.name)
        break
else:
    raise SystemExit(f"best product not found in registry: {best}")
PY
)"

template="${TEMPLATE:-macros/analysis_config_auau_bdt_mlp_stack_template.yaml}"
if [[ ! -f "$template" ]]; then
  template="macros/analysis_config.yaml"
fi
yaml="${config_dir}/analysis_config_auau_logreg_${best_product}_wp080_${stamp}.yaml"
bash ./scripts/auau_tight_logreg_pipeline.sh generateWorkingPointConfig \
  TEMPLATE="${template}" \
  WORKING_POINTS="${report_dir}/logreg_working_points_target80.json" \
  ARTIFACT="${artifact}" \
  PRODUCT="${best_product}" \
  OUT="${yaml}"

production_block="${config_dir}/run_logreg_targetwp_pair_${stamp}.sh"
cat > "$production_block" <<EOF
#!/usr/bin/env bash
set -euo pipefail
cd "${repo_root}"
export RJ_LOGREG_TARGETWP_CONFIG_YAML="${yaml}"
export RJ_LOGREG_TARGETWP_CAMPAIGN_TAG="logreg_${best_product}_wp080_${stamp}"
export RJ_NOTIFY_EMAILS="${RJ_NOTIFY_EMAILS:-just0131@gmail.com}"
export RJ_REQUEST_MEMORY="${RJ_REQUEST_MEMORY:-12000MB}"
export RJ_AUTO_MEMORY_RETRY_CAP_MB="${RJ_AUTO_MEMORY_RETRY_CAP_MB:-16000}"
export RJ_SIM_FIRSTROUND_REQUEST_MEMORY="${RJ_SIM_FIRSTROUND_REQUEST_MEMORY:-8000MB}"
export RJ_SIM_MERGE_GROUP_SIZE="${RJ_SIM_MERGE_GROUP_SIZE:-75}"
bash ./scripts/submit_auau_logreg_targetwp_pair.sh "\$RJ_LOGREG_TARGETWP_CONFIG_YAML"
EOF
chmod +x "$production_block"

cat <<EOF

RECOILJETS_AUAU_TIGHT_LOGREG_FULL_CHAIN_READY_V1
status=READY_FOR_USER_REVIEW
best_product=${best_product}
artifact=${artifact}
validation_report=${report_dir}
rank_table=${report_dir}/validation_rank_table.csv
working_points=${report_dir}/logreg_working_points_target80.json
generated_yaml=${yaml}
production_block=${production_block}
note=Production pair is prepared but not submitted. Review validation/WP80 diagnostics before running production_block.
EOF
