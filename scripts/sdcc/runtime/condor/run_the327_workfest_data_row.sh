#!/usr/bin/env bash
# DATA adapter candidate: final runtime binding and native offset/batch validation required.
set -euo pipefail
umask 022

# Immutable one-pair schema-14 DATA worker for the THE-327 Workfest successor.
# The helper modes below are the same validators used by the real worker and
# exist so packet builders can exercise the fail-closed contracts without ROOT
# or SDCC. They do not execute analysis, publish output, or contact a scheduler.

sha_file(){ /usr/bin/sha256sum "$1" | /usr/bin/awk '{print $1}'; }

if [[ "${1:-}" == "__validate_row" ]]; then
  shift
  exec python3 - "$@" <<'PY'
import hashlib, json, math, re, sys
from pathlib import Path

if len(sys.argv) != 19:
    raise SystemExit("row validator argument count differs")
(manifest_raw, manifest_sha, row_id, sample_id, records_raw, records_sha,
 offset_raw, size_raw, chunk_sha, planned_raw, plan_raw, plan_sha,
 original_raw, source_ordinal_raw, source_entry_raw, global_entry_raw,
 worker_self_raw, plan_locator) = sys.argv[1:]
SHA = re.compile(r"[0-9a-f]{64}")

def regular(raw, label):
    path = Path(raw)
    if not path.is_absolute() or not path.is_file() or path.is_symlink():
        raise SystemExit(f"{label} is absent or unsafe")
    return path

def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def positive(raw, label):
    if not re.fullmatch(r"[1-9][0-9]*", raw):
        raise SystemExit(f"{label} is not positive")
    return int(raw)

def nonnegative(raw, label):
    if not re.fullmatch(r"[0-9]+", raw):
        raise SystemExit(f"{label} is invalid")
    return int(raw)

manifest = regular(manifest_raw, "production manifest")
if not SHA.fullmatch(manifest_sha) or digest(manifest) != manifest_sha:
    raise SystemExit("production manifest differs")
payload = json.loads(manifest.read_text(encoding="utf-8"))
if (payload.get("schema") != "THE327Schema14DATAProductionManifestV1" or
        payload.get("status") != "PASS_FROZEN_UNSUBMITTED" or
        payload.get("submission_ready") is not True):
    raise SystemExit("production manifest state differs")
if payload.get("source_schema_version") != 14:
    raise SystemExit("source schema version differs")
source_schema = regular(payload.get("source_schema_path", ""), "source schema")
if (not SHA.fullmatch(str(payload.get("source_schema_sha256", ""))) or
        digest(source_schema) != payload["source_schema_sha256"]):
    raise SystemExit("source schema hash differs")
if not re.fullmatch(r"[A-Za-z0-9_.:-]+", row_id):
    raise SystemExit("row identity differs")

samples = payload.get("samples")
if not isinstance(samples, dict) or sample_id not in samples:
    raise SystemExit("sample is not admitted")
sample = samples[sample_id]
for key in ("system", "lane", "dataset", "sample", "si_di_role", "period"):
    if not isinstance(sample.get(key), str) or not sample[key]:
        raise SystemExit(f"sample field absent: {key}")
if (sample["system"], sample["si_di_role"]) not in {("pp", "DATA"), ("auau", "DATA")}:
    raise SystemExit("only paired pp/AuAu DATA is admitted")
normalization = sample.get("normalization")
if not isinstance(normalization, dict) or not normalization:
    raise SystemExit("normalization contract is absent")
roles = sample.get("source_roles")
expected_roles = ["jets", "jetcalo"]
if roles != expected_roles:
    raise SystemExit("source role order differs")

runtime = payload.get("runtime")
if not isinstance(runtime, dict):
    raise SystemExit("runtime binding is absent")
system = sample["system"]
prefix = "pp" if system == "pp" else "auau"
common_files = (
    "worker", "unified_macro", "calo_calib", "progress_header",
    "capture_contract", "dependency_manifest", "runtime_seal_receipt",
    "source_entry_audit", "source_entry_contract",
    "pp_library", "auau_library",
)
system_files = (f"{prefix}_wrapper", f"{prefix}_macro", f"{prefix}_config")
if system == "auau":
    system_files += ("centrality_contract", "centrality_manifest")
for stem in common_files + system_files:
    path = regular(runtime.get(f"{stem}_path", ""), f"runtime file {stem}")
    wanted = runtime.get(f"{stem}_sha256")
    if not isinstance(wanted, str) or not SHA.fullmatch(wanted) or digest(path) != wanted:
        raise SystemExit(f"runtime file differs: {stem}")
dependencies = json.loads(Path(runtime["dependency_manifest_path"]).read_text(encoding="utf-8"))
entries = dependencies.get("entries")
if (dependencies.get("schema") != "WorkfestExternalDependenciesV1" or
        dependencies.get("source_schema_version") != 14 or
        not isinstance(entries, list) or not 1 <= len(entries) <= 64):
    raise SystemExit("external dependency manifest shape differs")
seen_dependencies = set()
selected_dependencies = 0
for entry in entries:
    if not isinstance(entry, dict):
        raise SystemExit("external dependency entry differs")
    raw = entry.get("path")
    systems = entry.get("systems")
    size = entry.get("size_bytes")
    wanted = entry.get("sha256")
    physical = entry.get("physical_path")
    if (not isinstance(raw, str) or raw in seen_dependencies or
            not isinstance(systems, list) or not systems or
            len(systems) != len(set(systems)) or not set(systems) <= {"pp", "auau"} or
            isinstance(size, bool) or not isinstance(size, int) or not 0 < size <= 536870912 or
            not isinstance(wanted, str) or not SHA.fullmatch(wanted) or
            not isinstance(physical, str) or not Path(physical).is_absolute()):
        raise SystemExit("external dependency identity differs")
    seen_dependencies.add(raw)
    if system not in systems:
        continue
    path = regular(raw, "external dependency")
    if (str(path.resolve(strict=True)) != physical or path.stat().st_size != size or
            digest(path) != wanted):
        raise SystemExit(f"external dependency differs: {raw}")
    selected_dependencies += 1
if selected_dependencies == 0:
    raise SystemExit("system has no bound external dependencies")
worker_self = regular(worker_self_raw, "executing worker")
if worker_self.resolve(strict=True) != Path(runtime["worker_path"]).resolve(strict=True):
    raise SystemExit("executing worker path differs from immutable runtime binding")
for key in ("release", "offline_main", "durable_runtime_root", "dependency_runtime_root",
            "pp_install_prefix", "auau_install_prefix"):
    if not isinstance(runtime.get(key), str) or not runtime[key]:
        raise SystemExit(f"runtime field absent: {key}")
for key in ("durable_runtime_root", "dependency_runtime_root",
            "pp_install_prefix", "auau_install_prefix"):
    path = Path(runtime[key])
    if not path.is_absolute() or not path.is_dir() or path.is_symlink():
        raise SystemExit(f"runtime directory differs: {key}")
if not Path(runtime["offline_main"]).is_absolute():
    raise SystemExit("offline runtime root differs")

environment = runtime.get(f"{prefix}_replay_environment")
if not isinstance(environment, dict):
    raise SystemExit("replay environment is absent")
allowed = set(runtime.get("allowed_environment_keys", []))
if not allowed:
    raise SystemExit("allowed environment keys are absent")
protected = {
    "RJ_REPLAY_DATASET", "RJ_REPLAY_LANE", "RJ_REPLAY_SAMPLE",
    "RJ_REPLAY_SI_DI_ROLE", "RJ_REPLAY_OWNERSHIP_STATE",
    "RJ_REPLAY_INPUT_URI_SHA256", "RJ_REPLAY_INPUT_FILE_SHA256",
    "RJ_REPLAY_SOURCE_MANIFEST_SHA256", "RJ_REPLAY_SOURCE_SHA256",
    "RJ_REPLAY_PROVENANCE_MANIFEST_SHA256", "RJ_REPLAY_SEGMENT",
    "RJ_REPLAY_SEMANTIC_SHA256",
    "RJ_REPLAY_PERIOD", "RJ_REPLAY_SOURCE_FILE_ORDINAL",
    "RJ_REPLAY_SOURCE_ENTRY_BEGIN", "RJ_REPLAY_SOURCE_GLOBAL_ENTRY_BEGIN",
    "RJ_PROGRESS_PATH", "RJ_PROGRESS_TASK_ID", "RJ_PROGRESS_WORKSTREAM_ID",
    "RJ_PROGRESS_CAMPAIGN_TAG", "RJ_PROGRESS_ROW_ID", "RJ_PROGRESS_EXECUTOR",
    "RJ_PROGRESS_SKIP_TOTAL", "RJ_PROGRESS_FD", "RJ_OUTPUT_ROOT_CANDIDATE",
    "RJ_EVENT_OFFSET", "RJ_REPLAY_RUN", "RJ_RUNNUMBER", "RJ_TIMESTAMP", "RJ_CDB_TIMESTAMP",
}
for key, value in environment.items():
    if (key not in allowed or key in protected or
            not re.fullmatch(r"RJ_[A-Z0-9_]+", key) or
            not isinstance(value, str) or any(char in value for char in "\t\r\n")):
        raise SystemExit(f"unsafe replay environment: {key}")
snapshot_raw = environment.get("RJ_SNAPSHOT_LIB_DIR", "")
snapshot = Path(snapshot_raw)
if not snapshot.is_absolute() or not snapshot.is_dir() or snapshot.is_symlink():
    raise SystemExit("native-tested snapshot library directory is absent")
# The worker never substitutes a historical constituent cutoff. The frozen
# runtime must state the capture value explicitly and the worker exports it.
capture_value = environment.get("RJ_REPLAY_JET_CONSTITUENT_PT_MIN")
try:
    capture_number = float(capture_value) if capture_value is not None else float("nan")
    if not math.isfinite(capture_number) or capture_number < 0:
        raise ValueError
except ValueError:
    raise SystemExit("explicit capture threshold is absent or invalid")

if environment.get("RJ_SCHEMA10_DATA_RETENTION_PROFILE") != "workfest_replay_capture_v1":
    raise SystemExit("DATA capture must be Workfest, not sparse")
if environment.get("RJ_REPLAY_SAVE_JET_CONSTITUENTS") != "0":
    raise SystemExit("jet constituent saving must be disabled")
if system == "auau" and any(k.startswith("RJ_AUAU_CENTRALITY_") for k in environment):
    raise SystemExit("AuAu centrality must be derived per run, not inherited from a static environment")
records = regular(records_raw, "tuple records")
if (not SHA.fullmatch(records_sha) or sample.get("tuple_records_sha256") != records_sha
        or sample.get("tuple_records_path") != str(records)):
    raise SystemExit("tuple records binding differs")
offset = nonnegative(offset_raw, "tuple offset")
size = positive(size_raw, "tuple byte count")
if offset + size > records.stat().st_size or size > 16384:
    raise SystemExit("tuple slice exceeds records")
with records.open("rb") as stream:
    stream.seek(offset)
    chunk = stream.read(size)
if not SHA.fullmatch(chunk_sha) or hashlib.sha256(chunk).hexdigest() != chunk_sha:
    raise SystemExit("tuple hash differs")
if not chunk.endswith(b"\n") or chunk.count(b"\n") != 1:
    raise SystemExit("row is not one atomic tuple")
values = chunk[:-1].decode("utf-8").split("\t")
if len(values) != 2:
    raise SystemExit("tuple shape differs")
tuple_row = dict(zip(expected_roles, values))
for role, value in tuple_row.items():
    allowed_none = False
    if value == "NONE":
        if not allowed_none:
            raise SystemExit(f"required stream is NONE: {role}")
    elif not value.startswith("/") or any(char.isspace() for char in value):
        raise SystemExit(f"source path differs: {role}")

planned = positive(planned_raw, "planned GLOBAL event count")
if planned > (35000 if system == "auau" else 100000):
    raise SystemExit("planned GLOBAL event count exceeds row ceiling")
original = nonnegative(original_raw, "original ordinal")
source_ordinal = positive(source_ordinal_raw, "source file ordinal")
source_entry = nonnegative(source_entry_raw, "source entry begin")
global_entry = nonnegative(global_entry_raw, "source global entry begin")
if original + 1 != source_ordinal:
    raise SystemExit("source ordinal/entry coordinate differs")

plan = regular(plan_raw, "tuple event plan")
if (not SHA.fullmatch(plan_sha) or sample.get("tuple_plan_sha256") != plan_sha
        or sample.get("tuple_plan_path") != str(plan)):
    raise SystemExit("tuple event plan binding differs")
locator = re.fullmatch(r"([0-9]+):([1-9][0-9]*):([0-9a-f]{64})", plan_locator)
if not locator:
    raise SystemExit("tuple event plan locator differs")
plan_offset, plan_size = map(int, locator.groups()[:2])
if plan_size > 16384 or plan_offset + plan_size > plan.stat().st_size:
    raise SystemExit("tuple event plan slice exceeds bounds")
# The create-once materializer checks the full plan, unique occurrences and
# prefix sums once. Each job checks only its hash-bound row, not 59,998 rows.
with plan.open("rb") as stream:
    stream.seek(plan_offset)
    plan_bytes = stream.read(plan_size)
if hashlib.sha256(plan_bytes).hexdigest() != locator.group(3):
    raise SystemExit("tuple event plan row hash differs")
if not plan_bytes.endswith(b"\n") or plan_bytes.count(b"\n") != 1:
    raise SystemExit("tuple event plan slice is not one row")
plan_row = json.loads(plan_bytes)
if plan_row.get("sample_id") != sample_id or plan_row.get("source_file_ordinal") != source_ordinal:
    raise SystemExit("DATA plan occurrence differs")
if plan_row.get("row_id") != row_id:
    raise SystemExit("DATA row identity differs")
for key in ("run", "segment", "source_total_events", "event_offset", "event_count", "source_global_entry_begin"):
    value = plan_row.get(key)
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise SystemExit(f"invalid DATA plan field: {key}")
run = plan_row["run"]; segment = plan_row["segment"]
if run < 1000 or plan_row["source_total_events"] <= 0:
    raise SystemExit("invalid DATA run/source count")
if (plan_row["event_offset"] != source_entry or plan_row["event_count"] != planned
        or source_entry + planned > plan_row["source_total_events"]
        or plan_row["source_global_entry_begin"] != global_entry):
    raise SystemExit("DATA contiguous range differs")
if plan_row.get("source_paths") != tuple_row:
    raise SystemExit("DATA exact pair differs")
for role, path in tuple_row.items():
    match = re.search(r"-(\d{8})-(\d{5})\.root$", path)
    if not match or (int(match.group(1)), int(match.group(2))) != (run, segment):
        raise SystemExit("DATA pair run/segment identity differs")
if not Path(tuple_row["jets"]).name.startswith(("DST_Jet_", "DST_JET_")):
    raise SystemExit("JET source role differs")
if not Path(tuple_row["jetcalo"]).name.startswith("DST_JETCALO_"):
    raise SystemExit("JETCALO source role differs")

if system == "auau":
    import importlib.util
    module_spec=importlib.util.spec_from_file_location("canonical_centrality_contract",runtime["centrality_contract_path"])
    contract=importlib.util.module_from_spec(module_spec)
    sys.modules[module_spec.name]=contract
    module_spec.loader.exec_module(contract)
    approvals=()
    fallback=plan_row.get("centrality_fallback_approval")
    if fallback is not None:
        approval_path=regular(fallback.get("path", ""),"centrality fallback qualification")
        wanted=fallback.get("sha256", "")
        if not SHA.fullmatch(wanted) or digest(approval_path)!=wanted:
            raise SystemExit("centrality fallback evidence hash differs")
        approval=json.loads(approval_path.read_text())
        if (approval.get("schema")!="FrozenCentralityFallbackApprovalV1" or
                approval.get("status")!="PASS_COMPLETE_FALLBACK_QUALIFICATION" or
                approval.get("requested_run")!=run or approval.get("payload_run")!=68144 or
                approval.get("own_triplet_absent") is not True or
                approval.get("calibration_manifest_sha256")!=contract.NOMINAL_MANIFEST_SHA256):
            raise SystemExit("centrality fallback qualification does not cover this run/basis")
        approvals=(contract.FallbackApproval(run,wanted,True),)
    frozen=contract.manifest_from_frozen_bytes(Path(runtime["centrality_manifest_path"]).read_bytes(),
                                              [run],fallbacks=approvals)
    environment=contract.centrality_environment_for_run(run,frozen,inherited=environment)
    for field in ("DIVS","SCALE","VERTEX_SCALE"):
        key="RJ_AUAU_CENTRALITY_"+field
        if digest(regular(environment[key],"frozen centrality payload"))!=environment[key+"_SHA256"]:
            raise SystemExit("frozen centrality payload bytes differ")
    # Preserve the established AuAu finish-packet recovery binding. The C++
    # guard preserves saved DST status and restores HOT flags only when the
    # original source has no map and its input HOT flags are all absent.
    binding = runtime.get("cemc_recovery_by_run", {}).get(str(run))
    if not isinstance(binding, dict):
        raise SystemExit("CEMC recovery binding missing for run")
    map_path = regular(binding.get("path", ""), "CEMC recovery payload")
    wanted = binding.get("sha256", "")
    if (not SHA.fullmatch(str(wanted)) or digest(map_path) != wanted or
            "/cdb/CEMC_BadTowerMap/" not in str(map_path) or
            not str(map_path).endswith(f"_{run}cdb.root")):
        raise SystemExit("CEMC recovery payload run/hash differs")
    environment["RJ_AUAU_CEMC_RESTORE_MISSING_STATUS"] = "1"
    environment["RJ_AUAU_CEMC_RECOVERY_MAP"] = str(map_path)

canonical = runtime.get("canonical_assembly")
canonical_projection = None
if canonical is not None:
    base = Path(environment.get('RJ_RUNTIME_BASE', ''))
    if (not base.is_absolute() or not base.is_dir() or base.is_symlink()
            or environment.get('RJ_SNAPSHOT_LIB_PRECEDENCE') != 'prepend'):
        raise SystemExit('canonical capture requires explicit runtime base and prepend library binding')
    if not isinstance(canonical, dict) or canonical.get("schema") != "CanonicalBatchAssemblyConfigV1":
        raise SystemExit("canonical assembly configuration differs")
    entry = regular(canonical.get("entrypoint", ""), "canonical assembler")
    if entry.name != "finalize_canonical_capture.py":
        raise SystemExit("canonical assembler entrypoint differs")
    modules = canonical.get("module_sha256")
    if not isinstance(modules, dict) or not modules or len(modules) > 256:
        raise SystemExit("canonical module manifest differs")
    if set(modules) != {p.name for p in entry.parent.glob("*.py")}:
        raise SystemExit("canonical module directory is not completely sealed")
    for name, wanted in modules.items():
        if Path(name).name != name or not name.endswith(".py") or not SHA.fullmatch(str(wanted)):
            raise SystemExit("canonical module identity differs")
        if digest(regular(str(entry.parent/name), "canonical module")) != wanted:
            raise SystemExit("canonical module bytes differ")
    template = plan_row.get("canonical_augmentation_template", {})
    template_path = regular(template.get("path", ""), "canonical row template")
    template_locator = template.get('locator', '-')
    if template_locator == '-':
        template_digest = digest(template_path)
    else:
        match = re.fullmatch(r'([0-9]+):([1-9][0-9]*):([0-9a-f]{64})', str(template_locator))
        if not match:
            raise SystemExit('canonical template locator differs')
        start, length = map(int, match.groups()[:2])
        if length > 65536 or start + length > template_path.stat().st_size:
            raise SystemExit('canonical template slice exceeds bounds')
        with template_path.open('rb') as stream:
            stream.seek(start)
            template_bytes = stream.read(length)
        template_digest = hashlib.sha256(template_bytes).hexdigest()
        if (template_digest != match.group(3) or not template_bytes.endswith(b'\n')
                or template_bytes.count(b'\n') != 1):
            raise SystemExit('canonical template slice hash or framing differs')
    if not SHA.fullmatch(str(template.get("sha256", ""))) or template_digest != template["sha256"]:
        raise SystemExit("canonical row template bytes differ")
    product = {"isPP":"pp_data", "isAuAu":"auau_data"}.get(sample["dataset"])
    if product is None:
        raise SystemExit("canonical sample has no batch adapter")
    limits = [canonical.get(k) for k in ("max_seconds", "max_output_bytes", "max_rows")]
    if any(type(v) is not int or v <= 0 for v in limits) or limits[0] > 3600:
        raise SystemExit("canonical finite limits differ")
    codec = canonical.get('lossless_compression', 'inherited')
    if codec not in ('inherited', 'zstd5', 'lzma4'):
        raise SystemExit('canonical lossless compression policy differs')
    values = [str(entry), str(template_path), product, *map(str, limits), codec, template_locator]
    if any(any(c.isspace() for c in value) for value in values):
        raise SystemExit("canonical projection contains whitespace")
    canonical_projection = values

campaign = payload.get("campaign_tag")
if not isinstance(campaign, str) or not re.fullmatch(r"[A-Za-z0-9_.:-]+", campaign):
    raise SystemExit("campaign identity differs")
task = payload.get("task_id")
workstream = payload.get("workstream_id")
if any(not isinstance(v,str) or not re.fullmatch(r"[A-Za-z0-9_.:-]+",v) for v in (task,workstream)):
    raise SystemExit("campaign owner identity differs")
row = (
    system, sample["lane"], sample["dataset"], sample["sample"],
    runtime[f"{prefix}_wrapper_path"], runtime[f"{prefix}_macro_path"],
    runtime[f"{prefix}_config_path"], runtime["pp_install_prefix"],
    runtime["auau_install_prefix"], runtime["dependency_runtime_root"],
    runtime["release"], runtime["offline_main"], campaign,
    sample["period"], sample["si_di_role"], str(run), str(segment), task, workstream,
)
print("ROW\t" + "\t".join(row))
for key, value in sorted(environment.items()):
    print(f"ENV\t{key}\t{value}")
if canonical_projection:
    print("CANONICAL\t" + "\t".join(canonical_projection))
PY
fi

if [[ "${1:-}" == "__validate_terminal" ]]; then
  shift
  exec python3 - "$@" <<'PY'
import json, re, sys
from pathlib import Path
if len(sys.argv) != 11:
    raise SystemExit("terminal validator argument count differs")
metadata_path, progress_path, planned_raw, wrapper_raw, task, workstream, campaign, row_id, executor, skip_raw = sys.argv[1:]
if not re.fullmatch(r"[0-9]+", skip_raw): raise SystemExit("invalid skip count")
skip = int(skip_raw)
if not re.fullmatch(r"[1-9][0-9]*", planned_raw):
    raise SystemExit("planned count is invalid")
planned = int(planned_raw)
if not re.fullmatch(r"[0-9]+", wrapper_raw) or int(wrapper_raw) != 0:
    raise SystemExit("worker/wrapper success is absent")
metadata_file = Path(metadata_path)
progress_file = Path(progress_path)
if not metadata_file.is_file() or metadata_file.is_symlink():
    raise SystemExit("ROOT metadata projection is absent")
if not progress_file.is_file() or progress_file.is_symlink():
    raise SystemExit("live progress record is absent")
metadata = json.loads(metadata_file.read_text(encoding="utf-8"))
actual = metadata.get("rjevent_entries")
processed = metadata.get("processed_events")
if metadata.get("schema_version") != 14:
    raise SystemExit("replay schema differs")
if metadata.get("replay_complete") != 1:
    raise SystemExit("replay completion marker differs")
if any(isinstance(value, bool) or not isinstance(value, int) or value < 0
       for value in (actual, processed)):
    raise SystemExit("terminal event count is empty or invalid")
rejected=metadata.get("upstream_rejected_events")
if (type(rejected) is not int or rejected<0 or metadata.get("input_events")!=planned or
    actual!=processed or actual+rejected!=planned or
    metadata.get("source_entry_contract")!="ORIGINAL_PAIRED_DST_CURSOR_V1" or
    metadata.get("upstream_rejection_module")!="CaloStatusSkimmer" or
    metadata.get("upstream_rejection_return")!="ABORTEVENT" or
    metadata.get("source_identity_validated") is not True):
    raise SystemExit("input/retained/documented rejection or original-source identity differs")
progress = json.loads(progress_file.read_text(encoding="utf-8"))
expected = {
    "schema": "LongRunningCanaryProgressV1", "task_id": task,
    "workstream_id": workstream, "campaign_tag": campaign, "row_id": row_id,
    "executor": executor, "phase": "event_loop_complete", "processed": skip + planned,
    "total": skip + planned, "remaining": 0, "row_processed": planned,
    "row_total": planned, "row_remaining": 0, "skip_processed": skip,
    "skip_total": skip, "skip_remaining": 0,
    "deterministic_skip_processed": skip, "deterministic_skip_total": skip,
    "deterministic_skip_remaining": 0, "latest_error": "",
}
for key, value in expected.items():
    if progress.get(key) != value:
        raise SystemExit(f"live progress terminal field differs: {key}")
for key in ("rate_per_second", "eta_seconds", "updated_at"):
    if key not in progress:
        raise SystemExit(f"live progress terminal field is absent: {key}")
print(json.dumps({"actual_event_count": actual, "processed_events": processed,
                  "input_events":planned,"upstream_rejected_events":rejected,
                  "schema_version": 14, "replay_complete": 1}, sort_keys=True,
                 separators=(",", ":")))
PY
fi

if [[ "${1:-}" == "__build_receipt" ]]; then
  shift
  exec python3 - "$@" <<'PY'
import hashlib, json, math, os, re, sys
from pathlib import Path
if len(sys.argv) != 26:
    raise SystemExit("receipt builder argument count differs")
(time_raw, receipt_raw, candidate_raw, output_raw, request_raw, row_id, sample_id,
 cluster_raw, process_raw, manifest_raw, manifest_sha, records_sha, chunk_sha,
 plan_raw, plan_sha, planned_raw, actual_raw, processed_raw, original_raw,
 source_ordinal_raw, source_entry_raw, global_entry_raw, progress_raw, plan_locator, metadata_raw) = sys.argv[1:]
raw = Path(time_raw); receipt = Path(receipt_raw); candidate = Path(candidate_raw)
progress = Path(progress_raw); manifest = Path(manifest_raw)
for path, label in ((raw, "resource record"), (candidate, "candidate"),
                    (progress, "progress"), (manifest, "manifest"),
                    (Path(plan_raw), "tuple plan")):
    if not path.is_file() or path.is_symlink(): raise SystemExit(f"{label} is absent")
peaks = [int(m.group(1)) for line in raw.read_text(encoding="utf-8").splitlines()
         if (m := re.fullmatch(r"\s*Maximum resident set size \(kbytes\):\s*([0-9]+)\s*", line))]
if len(peaks) != 1 or peaks[0] <= 0: raise SystemExit("peak memory is absent")
time_text = raw.read_text(encoding="utf-8")
elapsed = re.findall(r"(?m)^\s*Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([0-9:.]+)\s*$", time_text)
if len(elapsed) != 1: raise SystemExit("terminal elapsed time is absent")
parts = elapsed[0].split(":")
if len(parts) not in (2, 3): raise SystemExit("terminal elapsed time is invalid")
runtime_seconds = 0.0
for part in parts: runtime_seconds = runtime_seconds * 60.0 + float(part)
if not math.isfinite(runtime_seconds) or runtime_seconds <= 0:
    raise SystemExit("terminal elapsed time is invalid")
cpu = {}
for name in ("User", "System"):
    values = re.findall(r"(?m)^\s*" + name + r" time \(seconds\):\s*([0-9.]+)\s*$", time_text)
    if len(values) != 1: raise SystemExit("terminal CPU time is absent")
    cpu[name.lower() + "_cpu_seconds"] = float(values[0])
    if not math.isfinite(cpu[name.lower() + "_cpu_seconds"]):
        raise SystemExit("terminal CPU time is invalid")
request = int(request_raw); peak_mb = math.ceil(peaks[0] / 1024)
if request <= 0 or request > 4096 or peak_mb > request: raise SystemExit("memory contract differs")
planned = int(planned_raw); actual = int(actual_raw); processed = int(processed_raw)
metadata=json.loads(Path(metadata_raw).read_text())
rejected=metadata.get('upstream_rejected_events')
if (planned<=0 or actual<0 or processed!=actual or type(rejected) is not int or rejected<0 or
    metadata.get('input_events')!=planned or actual+rejected!=planned or
    metadata.get('source_identity_validated') is not True or
    metadata.get('source_entry_contract')!='ORIGINAL_PAIRED_DST_CURSOR_V1'):
    raise SystemExit("receipt event counts differ")
payload_manifest = json.loads(manifest.read_text(encoding="utf-8"))
sample = payload_manifest["samples"][sample_id]
output_digest = hashlib.sha256()
with candidate.open("rb") as stream:
    for block in iter(lambda: stream.read(1024 * 1024), b""):
        output_digest.update(block)
payload = {
    "schema": "THE327Schema14DATAProductionRowReceiptV1",
    "status": "MEASURED_TERMINAL", "row_id": row_id, "sample_id": sample_id,
    "system": sample["system"], "period": sample["period"],
    "normalization": sample["normalization"], "source_roles": sample["source_roles"],
    "cluster_id": int(cluster_raw), "process_id": int(process_raw),
    "request_memory_mb": request, "raw_peak_memory_kb": peaks[0],
    "runtime_seconds": runtime_seconds, **cpu,
    "resource_time_sha256": hashlib.sha256(raw.read_bytes()).hexdigest(),
    "resource_time_raw": time_text,
    "peak_memory_mb": peak_mb, "production_manifest_path": str(manifest),
    "production_manifest_sha256": manifest_sha, "tuple_records_sha256": records_sha,
    "tuple_chunk_sha256": chunk_sha, "tuple_event_plan_path": plan_raw,
    "tuple_event_plan_sha256": plan_sha, "planned_event_count": planned,
    "tuple_event_plan_row_locator": plan_locator,
    "expected_event_count": planned, "actual_event_count": actual,
    "processed_events": processed,
    "event_count_contract": "EXACT_INPUT_EQUALS_RETAINED_PLUS_DOCUMENTED_UPSTREAM_REJECTED_V2",
    "input_events": planned, "upstream_rejected_events":rejected,
    "source_entry_contract":"ORIGINAL_PAIRED_DST_CURSOR_V1",
    "source_identity_validated":True,
    "upstream_rejection_module":"CaloStatusSkimmer", "upstream_rejection_return":"ABORTEVENT",
    "rejection_identities":"ReplayFoundationV1/RJUpstreamRejectedEventV1",
    "original_ordinal": int(original_raw), "source_file_ordinal": int(source_ordinal_raw),
    "source_entry_begin": int(source_entry_raw),
    "source_global_entry_begin": int(global_entry_raw),
    "scientific_output_path": output_raw,
    "scientific_output_sha256": output_digest.hexdigest(),
    "scientific_output_size_bytes": candidate.stat().st_size,
    "source_schema_version": 14, "source_schema_sha256": payload_manifest["source_schema_sha256"],
    "wrapper_exit_code": 0, "replay_complete": 1,
    "live_progress": {"path": str(progress),
        "sha256": hashlib.sha256(progress.read_bytes()).hexdigest(),
        "schema": "LongRunningCanaryProgressV1", "phase": "event_loop_complete",
        "processed": int(source_entry_raw)+planned, "total": int(source_entry_raw)+planned, "remaining": 0},
    "publication": "ATOMIC_NO_CLOBBER_LINK_RECEIPT_COMMIT_V1", "receipt_written_before_success": True,
}
tmp = receipt.with_name(receipt.name + f".tmp.{os.getpid()}")
tmp.write_text(json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n", encoding="utf-8")
os.chmod(tmp, 0o644); os.replace(tmp, receipt)
PY
fi

manifest="${1:-}"
manifest_sha="${2:-}"
row_id="${3:-}"
sample_id="${4:-}"
tuple_records_path="${5:-}"
tuple_records_sha="${6:-}"
byte_offset="${7:-}"
byte_count="${8:-}"
chunk_sha="${9:-}"
planned_event_count="${10:-}"
output_path="${11:-}"
measurement_path="${12:-}"
request_memory_mb="${13:-}"
original_ordinal="${14:-}"
tuple_plan_path="${15:-}"
tuple_plan_sha="${16:-}"
source_file_ordinal="${17:-}"
source_entry_begin="${18:-}"
source_global_entry_begin="${19:-}"
tuple_plan_locator="${20:-}"

progress_path=""
progress_enabled=0
write_failure_progress(){
  [[ "$progress_enabled" -eq 1 && -n "$progress_path" ]] || return 0
  python3 - "$progress_path" "$1" "$row_id" "$planned_event_count" <<'PY' || true
import datetime, json, os, sys
from pathlib import Path
path=Path(sys.argv[1]); message=sys.argv[2]; row_id=sys.argv[3]; row_total=int(sys.argv[4]); skip=int(os.environ["RJ_EVENT_OFFSET"]); total=row_total+skip
payload={"schema":"LongRunningCanaryProgressV1","processed":0,"total":total,
 "remaining":total,"row_processed":0,"row_total":row_total,"row_remaining":row_total,
 "skip_processed":0,"skip_total":skip,"skip_remaining":skip,
 "deterministic_skip_processed":0,"deterministic_skip_total":skip,
 "deterministic_skip_remaining":skip,"rate_per_second":0.0,"eta_seconds":None}
if path.is_file() and not path.is_symlink():
    try: payload=json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, json.JSONDecodeError): pass
processed=payload.get("processed",0); row_processed=payload.get("row_processed",0)
if isinstance(processed,bool) or not isinstance(processed,int) or processed<0: processed=0
if isinstance(row_processed,bool) or not isinstance(row_processed,int) or row_processed<0: row_processed=0
payload.update({"task_id":os.environ["RJ_PROGRESS_TASK_ID"],
 "workstream_id":os.environ["RJ_PROGRESS_WORKSTREAM_ID"],
 "campaign_tag":os.environ["RJ_PROGRESS_CAMPAIGN_TAG"],"row_id":row_id,
 "executor":os.environ["RJ_PROGRESS_EXECUTOR"],"phase":"worker_failed",
 "processed":min(processed,total),"total":total,"remaining":max(0,total-processed),
 "row_processed":min(row_processed,row_total),"row_total":row_total,
 "row_remaining":max(0,row_total-row_processed),"latest_error":message[:1000],
 "updated_at":datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds").replace("+00:00","Z")})
tmp=path.with_name(path.name+f".tmp.worker.{os.getpid()}")
tmp.write_text(json.dumps(payload,sort_keys=True,separators=(",",":"))+"\n",encoding="utf-8")
os.replace(tmp,path)
PY
}
die(){ message="$*"; write_failure_progress "$message"; printf '[THE327-SCHEMA14-DATA][ERROR] %s\n' "$message" >&2; exit 2; }

[[ "$#" -eq 20 ]] || die "worker argument count differs"
[[ -n "${_CONDOR_SCRATCH_DIR:-}" && -d "$_CONDOR_SCRATCH_DIR" ]] || die "worker sandbox is absent"
[[ "$output_path" == /sphenix/tg/tg01/bulk/*.root && ! -e "$output_path" && ! -L "$output_path" ]] || die "output path is unsafe or occupied"
[[ "$measurement_path" == /sphenix/tg/tg01/bulk/*.json && ! -e "$measurement_path" && ! -L "$measurement_path" ]] || die "measurement path is unsafe or occupied"
[[ "$request_memory_mb" =~ ^[1-9][0-9]*$ && "$request_memory_mb" -le 4096 ]] || die "memory request differs"
[[ "${ClusterId:-}" =~ ^[0-9]+$ && "${ProcId:-}" =~ ^[0-9]+$ ]] || die "scheduler identity is absent"

worker_self="$(cd "$(dirname "$0")" && pwd -P)/$(basename "$0")"
row_contract="$({ "$worker_self" __validate_row "$manifest" "$manifest_sha" "$row_id" "$sample_id" \
  "$tuple_records_path" "$tuple_records_sha" "$byte_offset" "$byte_count" "$chunk_sha" \
  "$planned_event_count" "$tuple_plan_path" "$tuple_plan_sha" "$original_ordinal" \
  "$source_file_ordinal" "$source_entry_begin" "$source_global_entry_begin" "$worker_self" "$tuple_plan_locator"; } 2>&1)" || \
  die "row/runtime contract validation failed: $row_contract"
first_line="${row_contract%%$'\n'*}"
IFS=$'\t' read -r marker system lane dataset sample wrapper macro config pp_install auau_install dependency_root release offline_main campaign period si_di_role run_number segment task_id workstream_id <<<"$first_line"
[[ "$marker" == ROW ]] || die "row/runtime projection differs"
environment_lines=""; [[ "$row_contract" == *$'\n'* ]] && environment_lines="${row_contract#*$'\n'}"
canonical_entry=""; canonical_template=""; canonical_product=""
canonical_seconds=""; canonical_bytes=""; canonical_rows=""; canonical_codec="inherited"; canonical_locator="-"
while IFS=$'\t' read -r env_marker key value rest; do
  [[ -z "$env_marker" ]] && continue
  if [[ "$env_marker" == CANONICAL ]]; then
    [[ -z "$canonical_entry" ]] || die "repeated canonical assembly projection"
    canonical_entry="$key"; canonical_template="$value"
    IFS=$'\t' read -r canonical_product canonical_seconds canonical_bytes canonical_rows canonical_codec canonical_locator <<<"$rest"
    continue
  fi
  [[ "$env_marker" == ENV && "$key" == RJ_* ]] || die "replay environment projection differs"
  [[ -z "$rest" ]] || die "replay environment contains unexpected fields"
  export "$key=$value"
done <<<"$environment_lines"

# Keep the native-tested external library directory exported by the frozen
# environment. The copied producer runtime does not contain shared/lib.
if [[ -z "$canonical_entry" ]]; then
  export RJ_RUNTIME_BASE="$dependency_root/shared" RJ_SNAPSHOT_LIB_PRECEDENCE=append
fi
export RJ_MACRO_PATH="$macro" RJ_CONFIG_YAML="$config"
export RJ_PINNED_RELEASE_NAME="$release" RJ_PINNED_OFFLINE_MAIN="$offline_main"
if [[ "$system" == pp ]]; then
  export RJ_PP_INSTALL_PREFIX="$pp_install"; unset RJ_AUAU_INSTALL_PREFIX RJ_SHARED_HEADER_INCLUDE
else
  export RJ_AUAU_INSTALL_PREFIX="$auau_install" RJ_SHARED_HEADER_INCLUDE="$pp_install/include/caloana"
  unset RJ_PP_INSTALL_PREFIX
fi
export RJ_DATASET="$dataset" RJ_IS_SIM=0 RJ_REPLAY_FOUNDATION_V1=1
unset RJ_SIM_SAMPLE
export RJ_EVENT_OFFSET="$source_entry_begin" RJ_REPLAY_RUN="$run_number" RJ_CDB_TIMESTAMP="$run_number"
# Both DATA producers expose the same audited paired-input cursor. AuAu's
# complete per-run frozen triplet was verified before exporting this row.
# Actual production admission remains bound to the sealed manifest/runtime.
export RJ_REQUIRE_ORIGINAL_SOURCE_CURSOR=1
export RJ_REPLAY_DATASET="$dataset" RJ_REPLAY_LANE="$lane" RJ_REPLAY_SAMPLE="$sample"
export RJ_REPLAY_SI_DI_ROLE="$si_di_role" RJ_REPLAY_PERIOD="$period" RJ_REPLAY_OWNERSHIP_STATE="$campaign"
export RJ_REPLAY_SOURCE_FILE_ORDINAL="$source_file_ordinal"
export RJ_REPLAY_SOURCE_ENTRY_BEGIN="$source_entry_begin"
export RJ_REPLAY_SOURCE_GLOBAL_ENTRY_BEGIN="$source_global_entry_begin"
export RJ_REPLAY_INPUT_URI_SHA256="$chunk_sha" RJ_REPLAY_INPUT_FILE_SHA256="$tuple_records_sha"
export RJ_REPLAY_SOURCE_MANIFEST_SHA256="$manifest_sha" RJ_REPLAY_SOURCE_SHA256="$manifest_sha"
export RJ_REPLAY_PROVENANCE_MANIFEST_SHA256="$manifest_sha" RJ_REPLAY_SEGMENT="$segment"
# Bind semantics to this frozen runtime/selection manifest, not the old
# two-event canary's source manifest. Runtime initialization requires this pin.
export RJ_REPLAY_SEMANTIC_SHA256="$manifest_sha"
# Force the existing macro's ordinary synchronized SIM manager route for every
# non-NONE source. These toggles do not pre-open or pre-read a DST.
export RJ_REQUIRE_G4=0 RJ_AUAU_CANDIDATE_SKIM_ONLY=0
export RJ_PPG12_PPSIM_G4_ONLY=0 RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4=0 RJ_PPG12_FIG11_SB_DIAGNOSTIC=0
export RJ_REQUIRE_NON_TINY_OUTPUT=1 RJ_MIN_OUTPUT_BYTES=50000 RJ_JOB_HEARTBEAT_SECONDS=0 RJ_PROFILE_JOB=0
export RJ_VERBOSITY=1 RJ_REQUEST_MEMORY_MB="$request_memory_mb"

progress_path="${measurement_path%.json}.progress.json"
[[ ! -e "$progress_path" && ! -L "$progress_path" ]] || die "live progress path is occupied"
export RJ_PROGRESS_TASK_ID="$task_id" RJ_PROGRESS_WORKSTREAM_ID="$workstream_id"
export RJ_PROGRESS_CAMPAIGN_TAG="$campaign" RJ_PROGRESS_ROW_ID="$row_id"
export RJ_PROGRESS_PATH="$progress_path" RJ_PROGRESS_SKIP_TOTAL="$source_entry_begin" RJ_PROGRESS_FD=-1
export RJ_PROGRESS_EXECUTOR="condor:${ClusterId}.${ProcId}@${HOSTNAME:-unknown}"
[[ "$RJ_PROGRESS_EXECUTOR" =~ ^[A-Za-z0-9_.:@-]+$ ]] || die "live progress executor identity differs"
progress_enabled=1
trap 'rc=$?; if [[ "$rc" -ne 0 ]]; then write_failure_progress "worker_exit_${rc}"; fi' EXIT

candidate_path="${output_path}.part"
measurement_candidate="${measurement_path}.part"
if [[ -n "$canonical_entry" ]]; then
  candidate_path="${_CONDOR_SCRATCH_DIR}/canonical-native.part.root"
  measurement_candidate="${_CONDOR_SCRATCH_DIR}/canonical-native.receipt.json"
fi
[[ ! -e "$candidate_path" && ! -L "$candidate_path" ]] || die "candidate path is occupied"
[[ ! -e "$measurement_candidate" && ! -L "$measurement_candidate" ]] || die "measurement candidate path is occupied"
chunk_path="${_CONDOR_SCRATCH_DIR}/input_tuple.tsv"
scratch_output="${_CONDOR_SCRATCH_DIR}/output"
/bin/dd if="$tuple_records_path" of="$chunk_path" iflag=skip_bytes,count_bytes skip="$byte_offset" count="$byte_count" status=none
[[ -s "$chunk_path" && "$(sha_file "$chunk_path")" == "$chunk_sha" ]] || die "tuple slice differs after extraction"
mkdir -p "$scratch_output" "$(dirname "$output_path")" "$(dirname "$measurement_path")"
time_raw="${_CONDOR_SCRATCH_DIR}/resource.time"
set +e
/usr/bin/time -v -o "$time_raw" bash "$wrapper" "$run_number" "$chunk_path" "$dataset" "$ClusterId" "$planned_event_count" 1 NONE "$scratch_output"
wrapper_status=$?
set -e
[[ "$wrapper_status" -eq 0 ]] || die "analysis wrapper exited ${wrapper_status}"
# Reuse the exact before-event validator to reject external/runtime drift
# during execution, before any scientific output is published.
"$worker_self" __validate_row "$manifest" "$manifest_sha" "$row_id" "$sample_id" \
  "$tuple_records_path" "$tuple_records_sha" "$byte_offset" "$byte_count" "$chunk_sha" \
  "$planned_event_count" "$tuple_plan_path" "$tuple_plan_sha" "$original_ordinal" \
  "$source_file_ordinal" "$source_entry_begin" "$source_global_entry_begin" "$worker_self" "$tuple_plan_locator" \
  >/dev/null || die "runtime or external dependency changed during execution"
mapfile -t roots < <(find "$scratch_output" -type f -name '*.root' -print)
[[ "${#roots[@]}" == 1 && -s "${roots[0]}" ]] || die "scientific ROOT output is absent or ambiguous"
cp -- "${roots[0]}" "$candidate_path"
chmod 0644 "$candidate_path"
[[ -s "$candidate_path" ]] || die "durable scientific ROOT candidate is absent"

metadata_path="${_CONDOR_SCRATCH_DIR}/root_metadata.json"
if ! root_error="$({
  set +u
  source /opt/sphenix/core/bin/sphenix_setup.sh -n "$release" >/dev/null
  setup_status=$?
  set -u
  [[ "$setup_status" -eq 0 ]] || exit "$setup_status"
  python3 - "$candidate_path" "$metadata_path" <<'PY'
import json, os, re, sys
from pathlib import Path
import ROOT
source=Path(sys.argv[1]); output=Path(sys.argv[2])
root=ROOT.TFile.Open(str(source),"READ")
if not root or root.IsZombie() or root.TestBit(ROOT.TFile.kRecovered): raise SystemExit("ROOT candidate open failed")
replay=root.GetDirectory("ReplayFoundationV1")
if not replay: raise SystemExit("ReplayFoundationV1 is absent")
def title(name):
    obj=replay.Get(name)
    if not obj or not obj.InheritsFrom("TNamed"): raise SystemExit(f"metadata absent: {name}")
    return str(obj.GetTitle())
events=replay.Get("RJEventV1")
if not events or not events.InheritsFrom("TTree"): raise SystemExit("RJEventV1 is absent")
expected_run=int(os.environ["RJ_REPLAY_RUN"])
expected_offset=int(os.environ["RJ_REPLAY_SOURCE_ENTRY_BEGIN"])
expected_source=int(os.environ["RJ_REPLAY_SOURCE_FILE_ORDINAL"])
expected_global=int(os.environ["RJ_REPLAY_SOURCE_GLOBAL_ENTRY_BEGIN"])
if title('source_entry_contract')!='ORIGINAL_PAIRED_DST_CURSOR_V1':
    raise SystemExit('original source-cursor receipt is absent')
input_count=int(title('input_events')); rejected_count=int(title('upstream_rejected_events'))
if title('upstream_rejection_module')!='CaloStatusSkimmer' or title('upstream_rejection_return')!='ABORTEVENT':
    raise SystemExit('unclassified upstream rejection')
rejected=replay.Get('RJUpstreamRejectedEventV1')
if not rejected or not rejected.InheritsFrom('TTree') or rejected.GetEntries()!=rejected_count:
    raise SystemExit('upstream rejection identities are absent')
rejected_entries=set(); physical_ids=set()
for row in rejected:
    entry=int(row.source_entry_ordinal); physical=int(row.physical_event_sequence)
    if (int(row.run)!=expected_run or not expected_offset<=entry<expected_offset+input_count or
        entry in rejected_entries or physical<0 or physical in physical_ids):
        raise SystemExit('invalid/duplicate rejected input identity')
    rejected_entries.add(entry);physical_ids.add(physical)
retained_entries=[entry for entry in range(expected_offset,expected_offset+input_count)
                  if entry not in rejected_entries]
if len(retained_entries)!=events.GetEntries():raise SystemExit('retained/rejected coverage differs')
events.SetBranchStatus("*",0)
for name in ("run", "source_entry_ordinal", "source_file_ordinal", "source_global_entry_ordinal",
             "physical_event_sequence", "physical_event_sequence_valid"):
    if not events.GetLeaf(name): raise SystemExit(f"DATA identity leaf absent: {name}")
    events.SetBranchStatus(name,1)
for index in range(int(events.GetEntriesFast())):
    if events.GetEntry(index) <= 0: raise SystemExit("DATA event identity read failed")
    if (int(events.GetLeaf("run").GetValueLong64()) != expected_run or
        int(events.GetLeaf("source_file_ordinal").GetValueLong64()) != expected_source or
        int(events.GetLeaf("source_entry_ordinal").GetValueLong64()) != retained_entries[index] or
        int(events.GetLeaf("source_global_entry_ordinal").GetValueLong64()) !=
            expected_global+retained_entries[index]-expected_offset or
        int(events.GetLeaf("physical_event_sequence_valid").GetValueLong64()) != 1):
        raise SystemExit("DATA run/source/entry identity differs")
    physical=int(events.GetLeaf("physical_event_sequence").GetValueLong64())
    if physical<0 or physical in physical_ids:raise SystemExit('invalid/duplicate retained physical identity')
    physical_ids.add(physical)
schema=title("rj_replay_schema_version"); complete=title("rj_replay_complete"); processed=title("processed_events")
if not re.fullmatch(r"[0-9]+",schema+complete+processed): raise SystemExit("ROOT integer metadata differs")
payload={"schema_version":int(schema),"replay_complete":int(complete),
 "processed_events":int(processed),"rjevent_entries":int(events.GetEntriesFast()),
 "input_events":input_count,"upstream_rejected_events":rejected_count,
 "source_entry_contract":title('source_entry_contract'),"source_identity_validated":True,
 "upstream_rejection_module":title('upstream_rejection_module'),
 "upstream_rejection_return":title('upstream_rejection_return')}
tmp=output.with_name(output.name+f".tmp.{os.getpid()}")
tmp.write_text(json.dumps(payload,sort_keys=True,separators=(",",":"))+"\n",encoding="utf-8")
os.replace(tmp,output); root.Close()
PY
} 2>&1)"; then
  printf '%s\n' "$root_error" >&2
  die "ROOT metadata validation failed; candidate preserved at ${candidate_path}"
fi
terminal="$({ "$worker_self" __validate_terminal "$metadata_path" "$progress_path" "$planned_event_count" "$wrapper_status" \
  "$RJ_PROGRESS_TASK_ID" "$RJ_PROGRESS_WORKSTREAM_ID" "$RJ_PROGRESS_CAMPAIGN_TAG" "$row_id" "$RJ_PROGRESS_EXECUTOR" "$source_entry_begin"; } 2>&1)" || \
  die "terminal contract failed; candidate preserved at ${candidate_path}: $terminal"
actual_events="$(python3 -c 'import json,sys; print(json.loads(sys.argv[1])["actual_event_count"])' "$terminal")"
processed_events="$(python3 -c 'import json,sys; print(json.loads(sys.argv[1])["processed_events"])' "$terminal")"

"$worker_self" __build_receipt "$time_raw" "$measurement_candidate" "$candidate_path" "$output_path" \
  "$request_memory_mb" "$row_id" "$sample_id" "$ClusterId" "$ProcId" "$manifest" "$manifest_sha" \
  "$tuple_records_sha" "$chunk_sha" "$tuple_plan_path" "$tuple_plan_sha" "$planned_event_count" \
  "$actual_events" "$processed_events" "$original_ordinal" "$source_file_ordinal" "$source_entry_begin" \
  "$source_global_entry_begin" "$progress_path" "$tuple_plan_locator" "$metadata_path" || die "terminal receipt failed; candidate preserved at ${candidate_path}"

if [[ -n "$canonical_entry" ]]; then
  (
    canonical_locator_args=()
    if [[ "$canonical_locator" != "-" ]]; then
      canonical_locator_args=(--augmentation-template-locator "$canonical_locator")
    fi
    set +u
    source /opt/sphenix/core/bin/sphenix_setup.sh -n "$release" >/dev/null || exit "$?"
    set -u
    /usr/bin/timeout --signal=TERM --kill-after=30s "${canonical_seconds}s" \
      python3 "$canonical_entry" "$candidate_path" "$output_path" \
      --request "$manifest" --producer-receipt "$measurement_candidate" \
      --receipt-output "$measurement_path" --product "$canonical_product" \
      --augmentation-template "$canonical_template" "${canonical_locator_args[@]}" \
      --lossless-compression "$canonical_codec" \
      --progress-output "${measurement_path%.json}.assembly-progress.json" \
      --max-seconds "$canonical_seconds" --max-output-bytes "$canonical_bytes" --max-rows "$canonical_rows"
  ) || die "canonical finalization failed; no committed output without matching receipt"
else
  # Same-filesystem atomic, no-clobber publication. Receipt is the commit marker.
  /bin/ln -- "$candidate_path" "$output_path" || die "final ROOT path occupied"
  /bin/ln -- "$measurement_candidate" "$measurement_path" || die "final receipt path occupied"
  /bin/rm -- "$candidate_path" "$measurement_candidate"
fi
trap - EXIT
printf '[THE327-SCHEMA14-DATA] PASS row=%s planned=%s actual=%s output=%s\n' "$row_id" "$planned_event_count" "$actual_events" "$output_path"
