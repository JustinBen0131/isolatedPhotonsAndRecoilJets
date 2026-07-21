#!/usr/bin/env bash
set -euo pipefail

# Submit exactly one validation-only THE-116 job.  It reconstructs the frozen
# deterministic event50 holdout, trains a non-promoted 5--35 GeV compatibility
# control from the identical weighted cache, certifies XGBoost/TMVA parity, and
# derives exact binned WP70/WP80/WP90 thresholds for the accepted candidate.

REPO_BASE="${RJ_REPO_BASE:-/sphenix/u/patsfan753/scratch/thesisAnalysis}"
TRAIN_TAG="${RJ_THE116_TRAIN_TAG:-the116_pp_matched_basev3e_15to35_20260720}"
TAG="${RJ_THE116_VALIDATION_TAG:-the116_pp_exact_holdout_validation_20260720}"
RUN_ROOT="${RJ_THE116_RUN_ROOT:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/the116_models/${TRAIN_TAG}}"
VALIDATION_ROOT="${RJ_THE116_VALIDATION_ROOT:-${RUN_ROOT}/validation/exact_holdout}"
SUBMIT_ROOT="${RJ_THE116_VALIDATION_SUBMIT_ROOT:-/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/${TAG}}"
ML_PYTHON="${RJ_ML_PYTHON:-/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python}"
MANIFEST="${RJ_THE116_MANIFEST:-/sphenix/u/patsfan753/scratch/thesisAnalysis/dataOutput/ppPhotonMLPipeline/ppg12_matched_ppg12_baseV3E_currentIAN_npbet5env_20260527_1356/training_roots_currentIAN.list}"
TRAINER="${REPO_BASE}/scripts/ml/training/train_auau_photon_bdt.py"
VALIDATOR="${REPO_BASE}/scripts/ml/validation/validate_the116_pp_bdt.py"
CANDIDATE_REGISTRY="${RUN_ROOT}/models/bdt_ppg12_basev3e_15to35/model_registry.json"

EXPECTED_MANIFEST_SHA="e9852be2112c3d91d9ab248731ee9d199e03213b6dd5c93e6f48bfdaaa206258"
EXPECTED_TRAINER_SHA="06de85c527e1f7d41a48b66046051bad70f238acc4af13bd5520a263a514b2c3"
EXPECTED_VALIDATOR_SHA="d8a879b304ef64861d27e5be710d6ab59fa14434a244ada0b83c7c6cb158740a"
EXPECTED_CANDIDATE_REGISTRY_SHA="ba29c46745aa79bcc70d0b4458043a411802045443edc129418cf62c788b3d96"
EXPECTED_CANDIDATE_XGB_SHA="6e1ccb0d2c76e10dbdb4ce0dd5aaa7e56bae81eb5b20f3a77ab6bcf22c749727"
EXPECTED_CANDIDATE_TMVA_SHA="228d4cb73f7dc945a613c5a604add71a372b7540c2dc8c630b533d215bb17b30"

die() { printf '[THE-116 validation][ERROR] %s\n' "$*" >&2; exit 2; }
sha256_file() { sha256sum "$1" | awk '{print $1}'; }
require_hash() {
  local path="$1" expected="$2" actual
  [[ -s "$path" ]] || die "missing required file: $path"
  actual="$(sha256_file "$path")"
  [[ "$actual" == "$expected" ]] || die "SHA mismatch for $path: expected=$expected actual=$actual"
}

[[ -n "${RJ_CODEX_CHAT_NAME:-}" ]] || die "RJ_CODEX_CHAT_NAME is required"
[[ -n "${RJ_CODEX_THREAD_ID:-}" ]] || die "RJ_CODEX_THREAD_ID is required"
[[ -x "$ML_PYTHON" ]] || die "ML Python is not executable: $ML_PYTHON"
require_hash "$MANIFEST" "$EXPECTED_MANIFEST_SHA"
require_hash "$TRAINER" "$EXPECTED_TRAINER_SHA"
require_hash "$VALIDATOR" "$EXPECTED_VALIDATOR_SHA"
require_hash "$CANDIDATE_REGISTRY" "$EXPECTED_CANDIDATE_REGISTRY_SHA"
[[ "$(wc -l < "$MANIFEST")" -eq 7000 ]] || die "training manifest does not contain exactly 7000 paths"

CANDIDATE_XGB="$($ML_PYTHON - "$CANDIDATE_REGISTRY" <<'PY'
import json, sys
print(json.load(open(sys.argv[1]))["models"][0]["output_xgb_json"])
PY
)"
CANDIDATE_TMVA="$($ML_PYTHON - "$CANDIDATE_REGISTRY" <<'PY'
import json, sys
print(json.load(open(sys.argv[1]))["models"][0]["output_tmva"])
PY
)"
require_hash "$CANDIDATE_XGB" "$EXPECTED_CANDIDATE_XGB_SHA"
require_hash "$CANDIDATE_TMVA" "$EXPECTED_CANDIDATE_TMVA_SHA"

mkdir -p "$SUBMIT_ROOT" "$VALIDATION_ROOT"
[[ ! -e "$VALIDATION_ROOT/validation_job_complete.json" ]] || die "validation completion marker already exists; refusing duplicate"
[[ ! -e "$SUBMIT_ROOT/submission_receipt.txt" ]] || die "submission receipt already exists; refusing duplicate"

CONFIG_ENV="$SUBMIT_ROOT/config.env"
{
  printf 'export %s=%q\n' TAG "$TAG"
  printf 'export %s=%q\n' REPO_BASE "$REPO_BASE"
  printf 'export %s=%q\n' RUN_ROOT "$RUN_ROOT"
  printf 'export %s=%q\n' VALIDATION_ROOT "$VALIDATION_ROOT"
  printf 'export %s=%q\n' ML_PYTHON "$ML_PYTHON"
  printf 'export %s=%q\n' MANIFEST "$MANIFEST"
  printf 'export %s=%q\n' TRAINER "$TRAINER"
  printf 'export %s=%q\n' VALIDATOR "$VALIDATOR"
  printf 'export %s=%q\n' CANDIDATE_REGISTRY "$CANDIDATE_REGISTRY"
  printf 'export %s=%q\n' EXPECTED_MANIFEST_SHA "$EXPECTED_MANIFEST_SHA"
  printf 'export %s=%q\n' EXPECTED_TRAINER_SHA "$EXPECTED_TRAINER_SHA"
  printf 'export %s=%q\n' EXPECTED_VALIDATOR_SHA "$EXPECTED_VALIDATOR_SHA"
  printf 'export %s=%q\n' EXPECTED_CANDIDATE_REGISTRY_SHA "$EXPECTED_CANDIDATE_REGISTRY_SHA"
  printf 'export %s=%q\n' EXPECTED_CANDIDATE_XGB_SHA "$EXPECTED_CANDIDATE_XGB_SHA"
  printf 'export %s=%q\n' EXPECTED_CANDIDATE_TMVA_SHA "$EXPECTED_CANDIDATE_TMVA_SHA"
  printf 'export %s=%q\n' RJ_CODEX_CHAT_NAME "$RJ_CODEX_CHAT_NAME"
  printf 'export %s=%q\n' RJ_CODEX_THREAD_ID "$RJ_CODEX_THREAD_ID"
} > "$CONFIG_ENV"
chmod 0600 "$CONFIG_ENV"

cat > "$VALIDATION_ROOT/submission_contract.json" <<EOF
{
  "schema": "THE116_PP_BDT_EXACT_HOLDOUT_SUBMISSION_V1",
  "tag": "$TAG",
  "training_tag": "$TRAIN_TAG",
  "manifest": "$MANIFEST",
  "manifest_sha256": "$EXPECTED_MANIFEST_SHA",
  "manifest_rows": 7000,
  "trainer": "$TRAINER",
  "trainer_sha256": "$EXPECTED_TRAINER_SHA",
  "validator": "$VALIDATOR",
  "validator_sha256": "$EXPECTED_VALIDATOR_SHA",
  "candidate_registry": "$CANDIDATE_REGISTRY",
  "candidate_registry_sha256": "$EXPECTED_CANDIDATE_REGISTRY_SHA",
  "candidate_xgb_sha256": "$EXPECTED_CANDIDATE_XGB_SHA",
  "candidate_tmva_sha256": "$EXPECTED_CANDIDATE_TMVA_SHA",
  "purpose": "validation-only deterministic event50 parity, 5-35 compatibility control, and exact weighted working points",
  "promotion_status": "NOT_PROMOTED",
  "codex_chat_name": "$RJ_CODEX_CHAT_NAME",
  "codex_thread_id": "$RJ_CODEX_THREAD_ID"
}
EOF

cat > "$SUBMIT_ROOT/worker.sh" <<'WORKER'
#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/config.env"

die() { printf '[THE-116 validation worker][ERROR] %s\n' "$*" >&2; exit 2; }
sha256_file() { sha256sum "$1" | awk '{print $1}'; }
require_hash() {
  local path="$1" expected="$2" actual
  [[ -s "$path" ]] || die "missing required file: $path"
  actual="$(sha256_file "$path")"
  [[ "$actual" == "$expected" ]] || die "SHA mismatch for $path: expected=$expected actual=$actual"
}

require_hash "$MANIFEST" "$EXPECTED_MANIFEST_SHA"
require_hash "$TRAINER" "$EXPECTED_TRAINER_SHA"
require_hash "$VALIDATOR" "$EXPECTED_VALIDATOR_SHA"
require_hash "$CANDIDATE_REGISTRY" "$EXPECTED_CANDIDATE_REGISTRY_SHA"

CANDIDATE_XGB="$($ML_PYTHON - "$CANDIDATE_REGISTRY" <<'PY'
import json, sys
print(json.load(open(sys.argv[1]))["models"][0]["output_xgb_json"])
PY
)"
CANDIDATE_TMVA="$($ML_PYTHON - "$CANDIDATE_REGISTRY" <<'PY'
import json, sys
print(json.load(open(sys.argv[1]))["models"][0]["output_tmva"])
PY
)"
require_hash "$CANDIDATE_XGB" "$EXPECTED_CANDIDATE_XGB_SHA"
require_hash "$CANDIDATE_TMVA" "$EXPECTED_CANDIDATE_TMVA_SHA"

CACHE="$VALIDATION_ROOT/ppg12_weighted_training_cache.npz"
CACHE_BUILDER="$VALIDATION_ROOT/cache_builder"
CONTROL_OUT="$VALIDATION_ROOT/control_5to35"
FINAL_OUT="$VALIDATION_ROOT/final"
LOG_DIR="$VALIDATION_ROOT/logs"
mkdir -p "$CACHE_BUILDER" "$CONTROL_OUT" "$FINAL_OUT" "$LOG_DIR"
exec > >(tee "$LOG_DIR/worker.stdout.log") 2> >(tee "$LOG_DIR/worker.stderr.log" >&2)
date -u +%FT%TZ > "$VALIDATION_ROOT/validation_job_started_utc.txt"

COMMON=(
  --task tight
  --tree AuAuPhotonIDTrainingTree
  --campaign ppg12-sixpack
  --campaign-spec-ids ppg12_base_v3E_bdt_noIso
  --cent-bins=-1:0
  --test-size 0.50
  --split-mode event50
  --random-seed 42
  --n-estimators 750
  --max-depth 5
  --learning-rate 0.1
  --subsample 0.5
  --colsample-bytree 0.6
  --tree-method hist
  --reg-alpha 5.0
  --reg-lambda 0.3
  --grow-policy lossguide
  --max-bin 256
  --n-jobs 2
  --no-event-weight
  --weight-mode ppg12-exact
  --ppg12-exact-expected-samples run28_photonjet5,run28_photonjet10,run28_photonjet20,run28_jet8,run28_jet12,run28_jet20,run28_jet30
  --ppg12-exact-closure-artifacts metadata-only
  --max-load-rows-per-class 2000000
  --load-sample-seed 42
  --majority-cap-ratio 0
)

if [[ ! -s "$CACHE" ]]; then
  "$ML_PYTHON" "$TRAINER" \
    "${COMMON[@]}" \
    --input "@$MANIFEST" \
    --outdir "$CACHE_BUILDER" \
    --pt-bins 5,10,14,18,22,35 \
    --skip-missing-tree \
    --cache-file "$CACHE" \
    --cache-only \
    --registry-output "$CACHE_BUILDER/model_registry.json"
fi
[[ -s "$CACHE" ]] || die "weighted training cache was not created"

if [[ ! -s "$CONTROL_OUT/model_registry.json" ]]; then
  "$ML_PYTHON" "$TRAINER" \
    "${COMMON[@]}" \
    --outdir "$CONTROL_OUT" \
    --pt-bins 5,10,14,18,22,35 \
    --cache-file "$CACHE" \
    --skip-tmva-export \
    --registry-output "$CONTROL_OUT/model_registry.json"
fi

"$ML_PYTHON" "$VALIDATOR" \
  --cache "$CACHE" \
  --candidate-registry "$CANDIDATE_REGISTRY" \
  --control-registry "$CONTROL_OUT/model_registry.json" \
  --outdir "$FINAL_OUT"

"$ML_PYTHON" - "$VALIDATION_ROOT" <<'PY'
import hashlib, json, pathlib, sys, time
root = pathlib.Path(sys.argv[1])
validation = json.loads((root / "final/the116_pp_bdt_validation.json").read_text())
if validation.get("status") != "PASS":
    raise SystemExit("THE-116 exact-holdout validation is not PASS")
files = {}
for relative in (
    "ppg12_weighted_training_cache.npz",
    "control_5to35/model_registry.json",
    "final/the116_pp_bdt_validation.json",
    "final/the116_pp_bdt_working_points.csv",
    "logs/worker.stdout.log",
):
    path = root / relative
    if not path.is_file() or path.stat().st_size <= 0:
        raise SystemExit(f"missing validation artifact: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    files[relative] = {"bytes": path.stat().st_size, "sha256": digest.hexdigest()}
(root / "validation_job_complete.json").write_text(json.dumps({
    "schema": "THE116_PP_BDT_EXACT_HOLDOUT_COMPLETE_V1",
    "completed_unix": time.time(),
    "status": validation["status"],
    "promotion_status": "NOT_PROMOTED",
    "files": files,
}, indent=2, sort_keys=True) + "\n")
PY
date -u +%FT%TZ > "$VALIDATION_ROOT/validation_job_finished_utc.txt"
WORKER
chmod 0755 "$SUBMIT_ROOT/worker.sh"

cat > "$SUBMIT_ROOT/the116_exact_holdout.submit" <<EOF
universe = vanilla
executable = $SUBMIT_ROOT/worker.sh
arguments =
initialdir = $SUBMIT_ROOT
output = $SUBMIT_ROOT/condor.\$(ClusterId).\$(ProcId).out
error = $SUBMIT_ROOT/condor.\$(ClusterId).\$(ProcId).err
log = $SUBMIT_ROOT/condor.\$(ClusterId).log
request_cpus = 2
request_memory = 48000MB
request_disk = 30000MB
should_transfer_files = NO
getenv = True
+JobBatchName = "$TAG"
+RJCampaign = "$TAG"
+RJCodexChatName = "${RJ_CODEX_CHAT_NAME}"
+RJCodexThreadId = "${RJ_CODEX_THREAD_ID}"
queue 1
EOF

if [[ "${RJ_DO_RUN:-0}" != "1" ]]; then
  printf '[THE-116 validation] dry-run ready\nsubmit=%s\nrun=%s\n' \
    "$SUBMIT_ROOT/the116_exact_holdout.submit" "$VALIDATION_ROOT"
  exit 0
fi

condor_submit -terse "$SUBMIT_ROOT/the116_exact_holdout.submit" | tee "$SUBMIT_ROOT/submission_receipt.txt"
