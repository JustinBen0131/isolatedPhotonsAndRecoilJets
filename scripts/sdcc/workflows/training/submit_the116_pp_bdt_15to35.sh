#!/usr/bin/env bash
set -euo pipefail

# Submit exactly one frozen THE-116 p+p baseV3E 15--35 GeV model job.
# The accepted 7,000-file extraction is reused; this script never submits DST
# extraction or a hyperparameter scan.

REPO_BASE="${RJ_REPO_BASE:-/sphenix/u/patsfan753/scratch/thesisAnalysis}"
TAG="${RJ_THE116_TAG:-the116_pp_matched_basev3e_15to35_20260720}"
RUN_ROOT="${RJ_THE116_RUN_ROOT:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/the116_models/${TAG}}"
SUBMIT_ROOT="${RJ_THE116_SUBMIT_ROOT:-/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/${TAG}}"
ML_PYTHON="${RJ_ML_PYTHON:-/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python}"
MANIFEST="${RJ_THE116_MANIFEST:-/sphenix/u/patsfan753/scratch/thesisAnalysis/dataOutput/ppPhotonMLPipeline/ppg12_matched_ppg12_basev3E_currentIAN_npbet5env_20260527_1356/training_roots_currentIAN.list}"
PREDECESSOR_REGISTRY="${RJ_THE116_PREDECESSOR_REGISTRY:-/sphenix/u/patsfan753/scratch/thesisAnalysis/dataOutput/ppPhotonMLPipeline/ppg12_matched_ppg12_basev3E_currentIAN_npbet5env_20260527_1356/models/bdt_ppg12_currentIAN_basev3E/model_registry.json}"
PIPELINE="${REPO_BASE}/scripts/sdcc/pipelines/pp/pp_photon_ml_pipeline.sh"
TRAINER="${REPO_BASE}/scripts/ml/training/train_auau_photon_bdt.py"
VALIDATOR="${REPO_BASE}/scripts/ml/validation/validate_pp_photon_ml_tables.py"

EXPECTED_MANIFEST_SHA="e9852be2112c3d91d9ab248731ee9d199e03213b6dd5c93e6f48bfdaaa206258"
EXPECTED_PREDECESSOR_SHA="34d90c17e3ba28b30231acdee6251dc5b8e10a87b3dc5fff16b95fbda4cde12d"
EXPECTED_PIPELINE_SHA="7395397e70fc834f7534fc44a993c9d55829bdf465de5da8e2fb9c58fe851233"
EXPECTED_TRAINER_SHA="06de85c527e1f7d41a48b66046051bad70f238acc4af13bd5520a263a514b2c3"
EXPECTED_PUBLIC_COMMIT="0cbec15a7aacd4a30beb1962af4ad009f63af501"

die() { printf '[THE-116][ERROR] %s\n' "$*" >&2; exit 2; }
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
[[ -x "$PIPELINE" ]] || die "pipeline is not executable: $PIPELINE"
[[ -s "$VALIDATOR" ]] || die "validator is missing: $VALIDATOR"
require_hash "$MANIFEST" "$EXPECTED_MANIFEST_SHA"
require_hash "$PREDECESSOR_REGISTRY" "$EXPECTED_PREDECESSOR_SHA"
require_hash "$PIPELINE" "$EXPECTED_PIPELINE_SHA"
require_hash "$TRAINER" "$EXPECTED_TRAINER_SHA"
[[ "$(wc -l < "$MANIFEST")" -eq 7000 ]] || die "training manifest does not contain exactly 7000 paths"

mkdir -p "$SUBMIT_ROOT" "$RUN_ROOT" "$RUN_ROOT/logs"
[[ ! -e "$RUN_ROOT/model_job_complete.json" ]] || die "THE-116 completion marker already exists; refusing duplicate training"

cat > "$RUN_ROOT/submission_contract.json" <<EOF
{
  "schema": "THE116_PP_BDT_SUBMISSION_CONTRACT_V1",
  "tag": "$TAG",
  "public_commit": "$EXPECTED_PUBLIC_COMMIT",
  "manifest": "$MANIFEST",
  "manifest_sha256": "$EXPECTED_MANIFEST_SHA",
  "manifest_rows": 7000,
  "predecessor_registry": "$PREDECESSOR_REGISTRY",
  "predecessor_registry_sha256": "$EXPECTED_PREDECESSOR_SHA",
  "pipeline_sha256": "$EXPECTED_PIPELINE_SHA",
  "trainer_sha256": "$EXPECTED_TRAINER_SHA",
  "pt_edges_gev": [15,17,19,21,23,25,28,32,35],
  "model_range_gev": [15,35],
  "signal_sources": ["run28_photonjet5","run28_photonjet10","run28_photonjet20"],
  "background_sources": ["run28_jet8","run28_jet12","run28_jet20","run28_jet30"],
  "weight_mode": "ppg12-exact",
  "split_mode": "event50",
  "split_seed": 42,
  "row_cap_per_class": 2000000,
  "load_sample_seed": 42,
  "chat_name": "${RJ_CODEX_CHAT_NAME}",
  "thread_id": "${RJ_CODEX_THREAD_ID}"
}
EOF

cat > "$SUBMIT_ROOT/worker.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export RJ_CODEX_CHAT_NAME='${RJ_CODEX_CHAT_NAME}'
export RJ_CODEX_THREAD_ID='${RJ_CODEX_THREAD_ID}'
export RJ_REPO_BASE='$REPO_BASE'
export RJ_ML_PYTHON='$ML_PYTHON'
export RJ_PP_PHOTON_ML_RUN_ROOT='$RUN_ROOT'
export CURRENT_IAN_TRAIN_MANIFEST='$MANIFEST'
export CURRENT_IAN_BDT_OUTDIR='$RUN_ROOT/models/bdt_ppg12_basev3e_15to35'
export CURRENT_IAN_CLOSURE_DIR='$RUN_ROOT/validation/ppg12_exact_reweight_closure'
export CURRENT_IAN_VALIDATION_OUTDIR='$RUN_ROOT/validation/currentian_bdt_15to35'
export RJ_PP_CURRENT_IAN_TRAIN_PT_BINS='15,17,19,21,23,25,28,32,35'
export RJ_PP_CURRENT_IAN_VALIDATION_PT_RANGE='15:35'
export PP_BDT_MAX_LOAD_ROWS_PER_CLASS=2000000
export PP_VALIDATION_MAX_LOAD_ROWS_PER_CLASS=1000000
export PP_LOAD_SAMPLE_SEED=42
export PP_BDT_N_JOBS=2
export PP_SKIP_TMVA_EXPORT=0
RUN_ROOT='$RUN_ROOT'
REPO_BASE='$REPO_BASE'
ML_PYTHON='$ML_PYTHON'
mkdir -p "\$RUN_ROOT/logs"
exec > >(tee "\$RUN_ROOT/logs/worker.stdout.log") 2> >(tee "\$RUN_ROOT/logs/worker.stderr.log" >&2)
date -u +%FT%TZ > "\$RUN_ROOT/model_job_started_utc.txt"
cd "\$REPO_BASE"
"\$REPO_BASE/scripts/sdcc/pipelines/pp/pp_photon_ml_pipeline.sh" trainCurrentIANBDT
"\$REPO_BASE/scripts/sdcc/pipelines/pp/pp_photon_ml_pipeline.sh" validateCurrentIANBDT
"\$ML_PYTHON" - "\$RUN_ROOT" <<'PY'
import hashlib, json, pathlib, sys, time
root = pathlib.Path(sys.argv[1])
registry = root / "models/bdt_ppg12_basev3e_15to35/model_registry.json"
payload = json.loads(registry.read_text())
if payload.get("status") != "READY" or payload.get("model_count") != 1:
    raise SystemExit("THE-116 registry is not READY with exactly one model")
files = {}
for path in sorted(root.rglob("*")):
    if path.is_file() and path.stat().st_size:
        files[str(path.relative_to(root))] = {
            "bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
(root / "model_job_complete.json").write_text(json.dumps({
    "schema": "THE116_MODEL_JOB_COMPLETE_V1",
    "completed_unix": time.time(),
    "registry_status": payload.get("status"),
    "model_count": payload.get("model_count"),
    "files": files,
}, indent=2, sort_keys=True) + "\n")
PY
date -u +%FT%TZ > "\$RUN_ROOT/model_job_finished_utc.txt"
EOF
chmod 0755 "$SUBMIT_ROOT/worker.sh"

cat > "$SUBMIT_ROOT/the116.submit" <<EOF
universe = vanilla
executable = $SUBMIT_ROOT/worker.sh
arguments =
initialdir = $SUBMIT_ROOT
output = $SUBMIT_ROOT/condor.\$(ClusterId).\$(ProcId).out
error = $SUBMIT_ROOT/condor.\$(ClusterId).\$(ProcId).err
log = $SUBMIT_ROOT/condor.\$(ClusterId).log
request_cpus = 2
request_memory = 48000MB
request_disk = 20000MB
should_transfer_files = NO
getenv = True
+JobBatchName = "$TAG"
+RJCampaign = "$TAG"
+RJCodeCommit = "$EXPECTED_PUBLIC_COMMIT"
+RJCodexChatName = "${RJ_CODEX_CHAT_NAME}"
+RJCodexThreadId = "${RJ_CODEX_THREAD_ID}"
queue 1
EOF

if [[ "${RJ_DO_RUN:-0}" != "1" ]]; then
  printf '[THE-116] dry-run ready\nsubmit=%s\nrun=%s\n' "$SUBMIT_ROOT/the116.submit" "$RUN_ROOT"
  exit 0
fi

condor_submit -terse "$SUBMIT_ROOT/the116.submit" | tee "$SUBMIT_ROOT/submission_receipt.txt"
