#!/usr/bin/env bash
# Foreground-only implementation for run_ppg12_recoiljets_paired_oracle.sh.
# This script is intentionally not a submitter and contains no merge/promotion
# path.  It writes only below a new, explicitly authorized output directory.

set -Eeuo pipefail

die() {
  printf 'PPG12_PAIRED_ORACLE_WORKER_FAIL: %s\n' "$*" >&2
  exit 2
}

[[ $# -eq 21 ]] || die "internal argument-count mismatch"

expected_token="$1"
repo_root="$2"
output_dir="$3"
setup_script="$4"
lane_id="$5"
sample="$6"
period="$7"
interaction="$8"
ppg_macro="$9"
g4_full_list="${10}"
truthjet_full_list="${11}"
apply_bdt="${12}"
apply_config="${13}"
base_e_model="${14}"
base_v3e_model="${15}"
npb_model="${16}"
tower_mask="${17}"
recoil_runtime_manifest="${18}"
recoil_config="${19}"
reuse_ppg_raw_root="${20}"
reuse_ppg_raw_contract="${21}"
expected_roounfold_library_sha256="d135771391ae250bcb64c0889571825abe9924649485890e7a9c64648ee99062"
expected_roounfold_pcm_sha256="2d91962a7b42acf246c7a80339eee71ca2f7e6df18ef76051d24a83bc61d4244"
expected_roounfold_header_tree_sha256="ea9b923a8f6bc57b28027b7183b10e87246810c326208b1e36bf2b4b7491a458"
expected_recoeff_source_sha256="e9b25fdb6dd8a6bfbbad029cb90aaddc9489fdf2846c630ea63c8c41ac771eee"
expected_recoeff_roounfold_compat_sha256="f5a12905952a0f49a7521868935e7868eca7cf8de1facae12c26e0dd9b712891"
expected_estimator_revision="29f8223bd9b36dffab07961b597afa94185bbdf1"
recoeff_roounfold_compat_needle=', Form("response_matrix_full_%d", ieta), "", false));'
recoeff_roounfold_compat_replacement=', Form("response_matrix_full_%d", ieta), ""));'

case "$sample" in Photon5|Photon10|Photon20) ;; *) die "unsupported sample: $sample" ;; esac
case "$period" in
  0mrad) recoeff_var_type_suffix="0rad" ;;
  1p5mrad) recoeff_var_type_suffix="1p5mrad" ;;
  *) die "unsupported period: $period" ;;
esac
case "$interaction" in SI|DI) ;; *) die "unsupported interaction: $interaction" ;; esac
sample_lower="$(printf '%s' "$sample" | tr '[:upper:]' '[:lower:]')"
interaction_lower="$(printf '%s' "$interaction" | tr '[:upper:]' '[:lower:]')"
expected_lane_id="photon:${sample_lower}:${period}:${interaction_lower}"
[[ "$lane_id" == "$expected_lane_id" ]] || \
  die "embedded lane identity differs from physical lane: ${expected_lane_id}"
if [[ "$interaction" == DI ]]; then
  oracle_sample="${sample_lower}_double"
  recoil_sample="run28_photonjet${sample#Photon}_double"
else
  oracle_sample="$sample_lower"
  recoil_sample="run28_photonjet${sample#Photon}"
fi

worker_self="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)/$(basename "${BASH_SOURCE[0]}")"
driver_script="$(dirname "$worker_self")/run_ppg12_recoiljets_paired_oracle.sh"
ppg_wrapper="${repo_root}/macros/diagnostics/pp_currentian/Fun4All_ppg12_fixed_seed_oracle.C"
recoil_wrapper="${repo_root}/macros/diagnostics/pp_currentian/Fun4All_recoiljets_fixed_seed_oracle.C"
auditor="${repo_root}/scripts/diagnostics/pp_currentian/audit_ppg12_recoiljets_paired_oracle.py"
comparator="${repo_root}/scripts/diagnostics/pp_currentian/compare_ppg12_recoiljets_same_cluster_features.py"
aggregate_extractor="${repo_root}/scripts/diagnostics/pp_currentian/extract_ppg12_recoeff_executable_aggregate.py"

if [[ "$reuse_ppg_raw_root" == NONE && "$reuse_ppg_raw_contract" == NONE ]]; then
  reuse_mode="disabled"
elif [[ "$reuse_ppg_raw_root" != NONE && "$reuse_ppg_raw_contract" != NONE ]]; then
  reuse_mode="exact_contract_bound"
else
  die "raw PPG12 reuse requires both source and contract"
fi

manifest_role_path() {
  python3 - "$recoil_runtime_manifest" "$1" <<'PY'
from pathlib import Path
import hashlib
import json
import re
import sys

manifest = Path(sys.argv[1])
role = sys.argv[2]
try:
    data = json.loads(manifest.read_text())
except (OSError, json.JSONDecodeError) as exc:
    raise SystemExit(f"cannot read runtime manifest {manifest}: {exc}")
matches = [
    item
    for item in data.get("files", [])
    if isinstance(item, dict) and item.get("role") == role
]
if len(matches) != 1:
    raise SystemExit(f"runtime manifest role {role} has {len(matches)} matches")
row = matches[0]
path = Path(str(row.get("path", "")))
if not path.is_absolute() or not path.is_file() or path.stat().st_size <= 0:
    raise SystemExit(f"runtime manifest role {role} is not an absolute nonempty file: {path}")
expected = str(row.get("sha256", ""))
if not re.fullmatch(r"[0-9a-f]{64}", expected):
    raise SystemExit(f"runtime manifest role {role} has an invalid sha256")
actual = hashlib.sha256(path.read_bytes()).hexdigest()
if actual != expected:
    raise SystemExit(f"runtime manifest role {role} hash drifted")
print(path)
PY
}

recoeff_period_role="ppg_recoeff_period_config_${period}"
recoeff_period_config="$(manifest_role_path "$recoeff_period_role")" || \
  die "failed to resolve sealed ${period} RecoEff config"
recoeff_truth_vertex_role="ppg_recoeff_truth_vertex_reweight_${period}"
recoeff_truth_vertex_reweight="$(manifest_role_path "$recoeff_truth_vertex_role")" || \
  die "failed to resolve sealed ${period} truth-vertex reweight ROOT"
recoeff_yaml_cpp_header_receipt="$(manifest_role_path ppg_recoeff_yaml_cpp_header_tree_receipt)" || \
  die "failed to resolve sealed yaml-cpp header-tree receipt"

required_keys=(
  output_dir setup_script ppg_macro g4_full_list truthjet_full_list
  apply_bdt apply_config base_e_model base_v3e_model npb_model tower_mask
  recoil_runtime_manifest recoil_config ppg_recoeff_period_config
  ppg_recoeff_truth_vertex_reweight
  ppg_recoeff_yaml_cpp_header_tree_receipt
)
required_values=(
  "$output_dir" "$setup_script" "$ppg_macro" "$g4_full_list"
  "$truthjet_full_list" "$apply_bdt" "$apply_config" "$base_e_model"
  "$base_v3e_model" "$npb_model" "$tower_mask" "$recoil_runtime_manifest"
  "$recoil_config" "$recoeff_period_config" "$recoeff_truth_vertex_reweight"
  "$recoeff_yaml_cpp_header_receipt"
)
token_file_keys=(
  setup_script ppg_macro g4_full_list truthjet_full_list apply_bdt
  apply_config base_e_model base_v3e_model npb_model tower_mask
  recoil_runtime_manifest recoil_config driver_script worker_script
  ppg_wrapper recoil_wrapper comparator auditor aggregate_extractor
  ppg_recoeff_period_config ppg_recoeff_truth_vertex_reweight
  ppg_recoeff_yaml_cpp_header_tree_receipt
)
token_file_values=(
  "$setup_script" "$ppg_macro" "$g4_full_list" "$truthjet_full_list"
  "$apply_bdt" "$apply_config" "$base_e_model" "$base_v3e_model"
  "$npb_model" "$tower_mask" "$recoil_runtime_manifest" "$recoil_config"
  "$driver_script" "$worker_self" "$ppg_wrapper" "$recoil_wrapper"
  "$comparator" "$auditor" "$aggregate_extractor" "$recoeff_period_config"
  "$recoeff_truth_vertex_reweight" "$recoeff_yaml_cpp_header_receipt"
)
if [[ "$reuse_mode" == exact_contract_bound ]]; then
  token_file_keys+=(reuse_ppg_raw_root reuse_ppg_raw_contract)
  token_file_values+=("$reuse_ppg_raw_root" "$reuse_ppg_raw_contract")
fi
for ((index = 0; index < ${#token_file_keys[@]}; ++index)); do
  key="${token_file_keys[$index]}"
  path="${token_file_values[$index]}"
  [[ "$path" == /* && -f "$path" && -s "$path" ]] || \
    die "token-bound file is missing or empty for $key: $path"
done

sha256_file() {
  python3 - "$1" <<'PY'
from pathlib import Path
import hashlib
import sys
print(hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest())
PY
}

validate_pinned_roounfold_library() {
  local library="$1"
  [[ "$(sha256_file "$library")" == "$expected_roounfold_library_sha256" ]] || \
    die "sealed RooUnfold library differs from pinned historical digest"
}

preserve_reused_ppg_evidence() {
  local source_contract="$1"
  local source_raw_root="$2"
  local copied_raw_root="$3"
  local copied_ppg_log="$4"
  local receipt="$5"
  python3 - "$source_contract" "$source_raw_root" "$copied_raw_root" \
    "$copied_ppg_log" "$receipt" <<'PY'
from pathlib import Path
import hashlib
import json
import os
import shutil
import sys

source_contract, source_raw, copied_raw, copied_log, receipt = map(
    Path, sys.argv[1:]
)

def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()

try:
    contract = json.loads(source_contract.read_text())
except (OSError, json.JSONDecodeError) as exc:
    raise SystemExit(f"cannot read reused PPG12 contract: {exc}")
source_log = Path(str(contract.get("paths", {}).get("ppg_log", "")))
for label, path in (
    ("source contract", source_contract),
    ("source raw ROOT", source_raw),
    ("copied raw ROOT", copied_raw),
    ("source PPG12 log", source_log),
):
    if not path.is_absolute() or not path.is_file() or path.stat().st_size <= 0:
        raise SystemExit(f"reused PPG12 {label} is not an absolute nonempty file: {path}")

source_raw_hash = digest(source_raw)
copied_raw_hash = digest(copied_raw)
if copied_raw_hash != source_raw_hash:
    raise SystemExit("copied PPG12 raw ROOT differs from its reuse source")

copied_log.parent.mkdir(parents=True, exist_ok=True)
log_tmp = copied_log.with_name(f".{copied_log.name}.tmp.{os.getpid()}")
receipt_tmp = receipt.with_name(f".{receipt.name}.tmp.{os.getpid()}")
try:
    shutil.copy2(source_log, log_tmp)
    os.replace(log_tmp, copied_log)
    source_log_hash = digest(source_log)
    copied_log_hash = digest(copied_log)
    if copied_log_hash != source_log_hash:
        raise SystemExit("copied PPG12 reconstruction log differs from its reuse source")
    evidence = {
        "schema_version": 1,
        "mode": "exact_contract_bound",
        "source_contract": {
            "path": str(source_contract),
            "sha256": digest(source_contract),
        },
        "source_raw_root": {
            "path": str(source_raw),
            "sha256": source_raw_hash,
        },
        "copied_raw_root": {
            "path": str(copied_raw),
            "sha256": copied_raw_hash,
        },
        "source_ppg_log": {
            "path": str(source_log),
            "sha256": source_log_hash,
        },
        "copied_ppg_log": {
            "path": str(copied_log),
            "sha256": copied_log_hash,
        },
    }
    receipt_tmp.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n")
    os.replace(receipt_tmp, receipt)
finally:
    log_tmp.unlink(missing_ok=True)
    receipt_tmp.unlink(missing_ok=True)
PY
}

token_first_ph_seed=2991264730
token_pedestal_seed=4256268992
token_pedestal_sequence=534
token_ph_seed_sequence=2991264730,4256268992,2394322166,874466025,2240380304
derived_token="$({
  printf '%s\n' \
    'schema_version=3' \
    "lane=${sample}:${period}:${interaction}" \
    "lane_id=${lane_id}" \
    'rows=5' \
    'rng_mode=historical_fifo_replay_v2' \
    "first_ph_seed=${token_first_ph_seed}" \
    "pedestal_seed=${token_pedestal_seed}" \
    "pedestal=${token_pedestal_sequence}" \
    "ph_seed_sequence=${token_ph_seed_sequence}" \
    'source_graph=NONE,g4,truthjet,NONE,NONE' \
    'runtime=new.17' \
    "reuse_mode=${reuse_mode}"
  for ((index = 0; index < ${#required_keys[@]}; ++index)); do
    printf '%s=%s\n' "${required_keys[$index]}" "${required_values[$index]}"
  done
  for ((index = 0; index < ${#token_file_keys[@]}; ++index)); do
    printf '%s=%s sha256=%s\n' \
      "${token_file_keys[$index]}" "${token_file_values[$index]}" \
      "$(sha256_file "${token_file_values[$index]}")"
  done
} | LC_ALL=C sort | python3 -c \
  'import hashlib,sys; print("ppg12-oracle:" + hashlib.sha256(sys.stdin.buffer.read()).hexdigest())')"

[[ "${RJ_PPG12_PAIRED_ORACLE_RUN_TOKEN:-}" == "$expected_token" ]] || \
  die "missing exact authorization token from plan/run driver"
[[ "$expected_token" == ppg12-oracle:* ]] || die "malformed authorization token"
[[ "$derived_token" == "$expected_token" ]] || \
  die "token-bound path or file content changed after plan authorization"
[[ -z "${RANDOMSEED+x}" ]] || die "inherited shell RANDOMSEED is forbidden"
[[ "$output_dir" == /* && ! -e "$output_dir" ]] || \
  die "output directory must be absolute and absent: $output_dir"

# Start the scientific runtime from a genuinely empty environment.  Merely
# changing OFFLINE_MAIN is insufficient on SDCC: an inherited login shell can
# retain release_ana and user-install routes in PATH/LD_LIBRARY_PATH.
if [[ "${RJ_PPG12_PAIRED_ORACLE_CLEAN_ENV:-0}" != 1 ]]; then
  exec /usr/bin/env -i \
    HOME="${HOME:-/tmp}" \
    USER="${USER:-unknown}" \
    LOGNAME="${LOGNAME:-${USER:-unknown}}" \
    SHELL=/bin/bash \
    PATH=/usr/bin:/bin:/usr/sbin:/sbin \
    RJ_PPG12_PAIRED_ORACLE_CLEAN_ENV=1 \
    RJ_PPG12_PAIRED_ORACLE_RUN_TOKEN="$expected_token" \
    /bin/bash --noprofile --norc "$worker_self" "$@"
fi

for path in \
  "$setup_script" "$ppg_macro" "$g4_full_list" \
  "$truthjet_full_list" "$apply_bdt" "$apply_config" "$base_e_model" \
  "$base_v3e_model" "$npb_model" "$tower_mask" \
  "$recoil_runtime_manifest" "$recoil_config" "$recoeff_period_config" \
  "$recoeff_truth_vertex_reweight" "$recoeff_yaml_cpp_header_receipt" "$ppg_wrapper" \
  "$recoil_wrapper" "$auditor" "$comparator" "$aggregate_extractor"; do
  [[ "$path" == /* && -f "$path" && -s "$path" ]] || die "missing input: $path"
done

umask 0022
mkdir -p "$output_dir"
state_file="${output_dir}/RUN_STATE"
failure_report="${output_dir}/worker_failure.log"
failure_stderr_spool="${output_dir}/.worker_stderr.spool"
failure_stderr_fifo="${output_dir}/.worker_stderr.pipe"
failure_stage="worker_initialization"
failure_stage_log=""
failure_exit_status=""
failure_line=0
failure_command=""
printf 'RUNNING\n' > "$state_file"
completed=0

capture_worker_error() {
  local status="$1"
  local line="$2"
  local command="$3"
  failure_exit_status="$status"
  failure_line="$line"
  # BASH_COMMAND is source text rather than expanded argument values.  Bound it
  # anyway so a malformed compound command cannot make the receipt unbounded.
  failure_command="${command:0:2048}"
}

bounded_log_tail() {
  local path="$1"
  [[ -f "$path" ]] || return 0
  # Bound by both lines and bytes.  The byte cap is applied last so the final
  # receipt stays small even when one diagnostic line is unusually long.
  tail -n 160 -- "$path" 2>/dev/null | tail -c 32768 || true
}

finish_state() {
  local status=$?
  local report_tmp="${failure_report}.tmp.$$"
  trap - ERR EXIT
  if [[ $completed -eq 0 ]]; then
    printf 'FAILED\n' > "$state_file"
    # Restore the original stderr before waiting for the capture tee.  Closing
    # fd 2 is what delivers EOF to the FIFO, guaranteeing the spool is flushed
    # before its bounded tail is copied into the durable failure receipt.
    exec 2>&9
    exec 9>&-
    wait "$failure_tee_pid" || true
    {
      printf 'PPG12_PAIRED_ORACLE_WORKER_FAILURE_V1\n'
      printf 'stage=%s\n' "$failure_stage"
      printf 'exit_status=%s\n' "${failure_exit_status:-$status}"
      printf 'line=%s\n' "$failure_line"
      printf 'command='
      printf '%q\n' "$failure_command"
      printf 'stage_log=%s\n' "${failure_stage_log:-NONE}"
      printf 'stderr_tail_begin\n'
      bounded_log_tail "$failure_stderr_spool"
      printf '\nstderr_tail_end\n'
      if [[ -n "$failure_stage_log" && -f "$failure_stage_log" ]]; then
        printf 'stage_log_tail_begin\n'
        bounded_log_tail "$failure_stage_log"
        printf '\nstage_log_tail_end\n'
      fi
    } > "$report_tmp"
    mv -f "$report_tmp" "$failure_report"
  else
    exec 2>&9
    exec 9>&-
    wait "$failure_tee_pid" || true
    rm -f "$failure_report"
  fi
  rm -f "$failure_stderr_spool" "$failure_stderr_fifo" "$report_tmp"
  exit "$status"
}

# Preserve the worker's stderr in the foreground while keeping a temporary
# local spool for the bounded failure receipt.  The spool is removed on every
# ordinary success/failure exit and is never part of scientific evidence.
rm -f "$failure_stderr_spool" "$failure_stderr_fifo" "$failure_report"
mkfifo "$failure_stderr_fifo"
exec 9>&2
tee "$failure_stderr_spool" < "$failure_stderr_fifo" >&9 &
failure_tee_pid=$!
exec 2> "$failure_stderr_fifo"
rm -f "$failure_stderr_fifo"
trap 'capture_worker_error "$?" "$LINENO" "$BASH_COMMAND"' ERR
trap finish_state EXIT

input_dir="${output_dir}/inputs"
runtime_dir="${output_dir}/runtime"
ppg_dir="${output_dir}/ppg12"
recoil_dir="${output_dir}/recoiljets"
report_dir="${output_dir}/comparison"
recoeff_dir="${output_dir}/ppg12_recoeff"
baseline_layout="${recoeff_dir}/baseline"
trace_layout="${recoeff_dir}/trace"
mkdir -p "$input_dir" "$runtime_dir" "$ppg_dir" "$recoil_dir" "$report_dir" \
  "$baseline_layout/input/${oracle_sample}" "$baseline_layout/output" \
  "$trace_layout/input/${oracle_sample}" "$trace_layout/output"

# Collapse both source-locked period role families to the exact config and
# truth-vertex weight ROOT executed by this lane.  Downstream closure evidence
# hashes this selected manifest, so a 0mrad/1p5mrad swap cannot hide behind a
# static source manifest that contains both periods.
paired_runtime_manifest="${runtime_dir}/paired_runtime_manifest.json"
failure_stage="runtime_manifest_selection"
python3 - "$recoil_runtime_manifest" "$paired_runtime_manifest" "$period" \
  "$recoeff_period_config" "$recoeff_truth_vertex_reweight" <<'PY'
from pathlib import Path
import hashlib
import json
import re
import sys

source_path, output_path, period, selected_config_path, selected_vertex_path = sys.argv[1:]
source = json.loads(Path(source_path).read_text())
files = source.get("files")
if not isinstance(files, list):
    raise SystemExit("source runtime manifest files must be a list")
period_role_families = {
    "ppg_recoeff_period_config": {
        "0mrad": "ppg_recoeff_period_config_0mrad",
        "1p5mrad": "ppg_recoeff_period_config_1p5mrad",
    },
    "ppg_recoeff_truth_vertex_reweight": {
        "0mrad": "ppg_recoeff_truth_vertex_reweight_0mrad",
        "1p5mrad": "ppg_recoeff_truth_vertex_reweight_1p5mrad",
    },
}
if period not in {"0mrad", "1p5mrad"}:
    raise SystemExit(f"unsupported selected period: {period}")
by_role = {}
seen_paths = set()
for row in files:
    if not isinstance(row, dict) or not isinstance(row.get("role"), str):
        raise SystemExit("source runtime manifest contains malformed file metadata")
    role = row["role"]
    if role in by_role:
        raise SystemExit(f"source runtime manifest duplicates role {role}")
    by_role[role] = row
for family in period_role_families.values():
    for role in family.values():
        if role not in by_role:
            raise SystemExit(f"source runtime manifest lacks period role {role}")
selected_paths = {
    "ppg_recoeff_period_config": selected_config_path,
    "ppg_recoeff_truth_vertex_reweight": selected_vertex_path,
}
selected_rows = []
for selected_role, selected_path in selected_paths.items():
    source_role = period_role_families[selected_role][period]
    selected = dict(by_role[source_role])
    if Path(str(selected.get("path", ""))).resolve() != Path(selected_path).resolve():
        raise SystemExit(f"selected {selected_role} path differs from resolved runtime role")
    actual = hashlib.sha256(Path(selected_path).read_bytes()).hexdigest()
    if selected.get("sha256") != actual:
        raise SystemExit(f"selected {selected_role} digest differs from source runtime manifest")
    selected["role"] = selected_role
    selected_rows.append(selected)
period_roles = {
    role
    for family in period_role_families.values()
    for role in family.values()
}
selected_files = [
    dict(row) for row in files if row.get("role") not in period_roles
]
selected_files.extend(selected_rows)
for row in selected_files:
    path = str(Path(str(row.get("path", ""))).resolve())
    if path in seen_paths:
        raise SystemExit(f"selected runtime manifest aliases physical path {path}")
    seen_paths.add(path)
output = dict(source)
output["selected_period"] = period
output["source_runtime_manifest"] = {
    "path": str(Path(source_path).resolve()),
    "sha256": hashlib.sha256(Path(source_path).read_bytes()).hexdigest(),
}
output["files"] = sorted(selected_files, key=lambda row: row["role"])
Path(output_path).write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
PY

# Remove inherited custom-release routing before sourcing the one common
# runtime.  Keep only the minimal operating-system path needed to run the
# setup script; ordinary SDCC login shells can otherwise retain ana.560 and a
# user install even when OFFLINE_MAIN is later changed to new.17.
failure_stage="new17_runtime_setup"
export PATH=/usr/bin:/bin:/usr/sbin:/sbin
unset OFFLINE_MAIN MYINSTALL ROOT_INCLUDE_PATH LD_LIBRARY_PATH PYTHONPATH \
  CMAKE_PREFIX_PATH CPATH CPLUS_INCLUDE_PATH LIBRARY_PATH PKG_CONFIG_PATH
set +e
set +u
# shellcheck disable=SC1090
source "$setup_script" -n new.17
setup_status=$?
set -u
set -e
[[ $setup_status -eq 0 ]] || die "sPHENIX new.17 setup failed"
expected_offline="/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_new/new.17"
[[ "${OFFLINE_MAIN:-}" == "$expected_offline" ]] || \
  die "setup did not resolve exact new.17 OFFLINE_MAIN: ${OFFLINE_MAIN:-<unset>}"
runtime_routes="${PATH:-}:${LD_LIBRARY_PATH:-}:${ROOT_INCLUDE_PATH:-}:${PYTHONPATH:-}:${CMAKE_PREFIX_PATH:-}"
if [[ "$runtime_routes" == *'/release/release_ana/'* || \
      "$runtime_routes" == *'/sphenix/user/patsfan753/install'* || \
      "$runtime_routes" == *'/sphenix/u/patsfan753/install'* ]]; then
  die "new.17 setup retained a forbidden analysis-release or user-install route"
fi
command -v root >/dev/null 2>&1 || die "ROOT is unavailable after new.17 setup"
command -v python3 >/dev/null 2>&1 || die "python3 is unavailable after new.17 setup"
base_ld_library_path="${LD_LIBRARY_PATH:-}"
base_root_include_path="${ROOT_INCLUDE_PATH:-}"

read -r first_ph_seed pedestal_seed pedestal_sequence ph_seed_sequence < <(
  python3 "$auditor" rng-contract
)
[[ -n "$first_ph_seed" && -n "$pedestal_seed" && -n "$pedestal_sequence" && \
   -n "$ph_seed_sequence" ]] || die "failed to derive reconstruction RNG contract"

g4_slice="${input_dir}/g4hits_first5.list"
truthjet_slice="${input_dir}/dst_truth_jet_first5.list"
combined_list="${input_dir}/recoil_first5.list"
source_pair_receipt="${input_dir}/source_pair_receipt.json"
failure_stage="source_pair_materialization"
python3 - "$g4_full_list" "$truthjet_full_list" "$g4_slice" \
  "$truthjet_slice" "$combined_list" "$source_pair_receipt" "$sample" \
  "$interaction" <<'PY'
from pathlib import Path
import hashlib
import json
import re
import sys

g4_full, truth_full, g4_out, truth_out, combined_out, receipt = map(
    Path, sys.argv[1:7]
)
sample = sys.argv[7]
interaction = sys.argv[8]

def rows(path: Path) -> list[str]:
    return [
        line.strip()
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]

g4 = rows(g4_full)[:5]
truth = rows(truth_full)[:5]
if len(g4) != 5 or len(truth) != 5:
    raise SystemExit("source lists do not contain five non-comment rows")
if len(set(g4)) != 5 or len(set(truth)) != 5:
    raise SystemExit("source lists contain duplicate first-five rows")

sample_number = sample.removeprefix("Photon")
if sample_number not in {"5", "10", "20"}:
    raise SystemExit(f"unsupported photon source sample: {sample}")
if interaction not in {"SI", "DI"}:
    raise SystemExit(f"unsupported interaction mode: {interaction}")
identity_stem = (
    f"PhotonJet{sample_number}"
    if interaction == "SI"
    else f"PhotonJet{sample_number}_pythia8_Detroit"
)
identity_pattern = re.compile(
    rf"{re.escape(identity_stem)}-\d{{10}}-\d{{6}}\.root"
)

def normalized_identity(value: str, prefix: str, label: str) -> str:
    basename = Path(value).name
    if not basename.startswith(prefix):
        raise SystemExit(
            f"{label} basename does not start with exact recognized prefix {prefix}: {basename}"
        )
    identity = basename[len(prefix):]
    if identity_pattern.fullmatch(identity) is None:
        raise SystemExit(f"{label} basename has unexpected sample/run/segment identity: {basename}")
    return identity

identities = []
canonical_rows = []
for index, (g4_row, truth_row) in enumerate(zip(g4, truth), start=1):
    g4_identity = normalized_identity(g4_row, "G4Hits_pythia8_", f"G4 row {index}")
    truth_identity = normalized_identity(
        truth_row, "DST_TRUTH_JET_pythia8_", f"truth row {index}"
    )
    if g4_identity != truth_identity:
        raise SystemExit(f"row {index} G4/truth identity mismatch")
    identities.append(g4_identity)
    canonical_rows.append(f"{g4_row}\t{truth_row}\t{g4_identity}")
if len(set(identities)) != 5:
    raise SystemExit("source pair identities are not unique")
pair_hash = hashlib.sha256(("\n".join(canonical_rows) + "\n").encode()).hexdigest()
g4_out.write_text("\n".join(g4) + "\n")
truth_out.write_text("\n".join(truth) + "\n")
combined_out.write_text(
    "\n".join(f"NONE {left} {right} NONE NONE" for left, right in zip(g4, truth))
    + "\n"
)
receipt.write_text(
    json.dumps(
        {
            "schema_version": 1,
            "normalization": {
                "g4_prefix": "G4Hits_pythia8_",
                "truthjet_prefix": "DST_TRUTH_JET_pythia8_",
                "basename_only": True,
            },
            "sample": sample,
            "interaction": interaction,
            "row_count": 5,
            "identities": identities,
            "pair_identity_sha256": pair_hash,
        },
        indent=2,
        sort_keys=True,
    )
    + "\n"
)
PY

failure_stage="runtime_asset_resolution"
recoil_macro="$(manifest_role_path recoil_macro)"
recoil_lib="$(manifest_role_path libRecoilJets.so)"
ppg_lib="$(manifest_role_path libCaloAna24.so)"
calo_reco_lib="$(manifest_role_path libcalo_reco.so)"
clusteriso_lib="$(manifest_role_path libclusteriso.so)"
jetbase_lib="$(manifest_role_path libjetbase.so)"
photon_builder_header="$(manifest_role_path PhotonClusterBuilder.h)"
photon_builder_include_root="$(dirname "$(dirname "$photon_builder_header")")"
recoeff_macro="$(manifest_role_path ppg_recoeff_macro)"
recoeff_trace_macro="$(manifest_role_path ppg_recoeff_trace_macro)"
recoeff_source_macro="$(manifest_role_path ppg_recoeff_source_macro)"
recoeff_compat_macro="$(manifest_role_path ppg_recoeff_roounfold_compat_macro)"
recoeff_compat_receipt="$(manifest_role_path ppg_recoeff_roounfold_compat_transform_receipt)"
recoeff_trace_receipt="$(manifest_role_path ppg_recoeff_trace_transform_receipt)"
recoeff_cross_section_header="$(manifest_role_path ppg_recoeff_cross_section_header)"
recoeff_truth_vertex_header="$(manifest_role_path ppg_recoeff_truth_vertex_header)"
recoeff_canonical_config="$(manifest_role_path ppg_recoeff_canonical_config)"
recoeff_yaml_cpp="$(manifest_role_path ppg_recoeff_yaml_cpp)"
recoeff_yaml_cpp_header_receipt="$(manifest_role_path ppg_recoeff_yaml_cpp_header_tree_receipt)"
recoeff_yaml_cpp_include_root="$(python3 - "$recoeff_yaml_cpp_header_receipt" <<'PY'
from pathlib import Path
import hashlib
import json
import sys

receipt_path = Path(sys.argv[1])
data = json.loads(receipt_path.read_text())
if data.get("schema_version") != 1 or data.get("role") != "ppg_recoeff_yaml_cpp_header_tree":
    raise SystemExit("invalid sealed yaml-cpp header-tree receipt")
include_root = Path(str(data.get("include_root", "")))
tree = Path(str(data.get("staged_tree", "")))
if not include_root.is_absolute() or tree.resolve() != (include_root / "yaml-cpp").resolve():
    raise SystemExit("yaml-cpp receipt has inconsistent include root")
rows = data.get("files")
if not isinstance(rows, list) or not rows:
    raise SystemExit("yaml-cpp receipt has no header inventory")
expected = {}
for row in rows:
    relative = str(row.get("relative_path", "")) if isinstance(row, dict) else ""
    digest = str(row.get("sha256", "")) if isinstance(row, dict) else ""
    if not relative or len(digest) != 64 or relative in expected:
        raise SystemExit("yaml-cpp receipt contains malformed header metadata")
    expected[relative] = digest
actual_paths = sorted(path for path in tree.rglob("*") if path.is_file())
actual = {str(path.relative_to(tree)): path for path in actual_paths}
if set(actual) != set(expected) or "yaml.h" not in actual:
    raise SystemExit("sealed yaml-cpp header-tree membership drifted")
tree_digest = hashlib.sha256()
for relative in sorted(actual):
    path = actual[relative]
    if tree.resolve() not in path.resolve().parents:
        raise SystemExit("sealed yaml-cpp header escapes its tree")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != expected[relative]:
        raise SystemExit(f"sealed yaml-cpp header drifted: {relative}")
    tree_digest.update(relative.encode())
    tree_digest.update(b"\0")
    tree_digest.update(bytes.fromhex(digest))
if tree_digest.hexdigest() != data.get("tree_sha256"):
    raise SystemExit("sealed yaml-cpp header-tree digest drifted")
print(include_root.resolve())
PY
)" || die "failed to validate sealed yaml-cpp header tree"
recoeff_roounfold="$(manifest_role_path ppg_recoeff_roounfold)"
recoeff_roounfold_pcm="$(manifest_role_path ppg_recoeff_roounfold_pcm)"
recoeff_roounfold_header_receipt="$(manifest_role_path ppg_recoeff_roounfold_header_tree_receipt)"
recoeff_roounfold_response_header="$(manifest_role_path ppg_recoeff_roounfold_response_header)"
recoeff_roounfold_bayes_header="$(manifest_role_path ppg_recoeff_roounfold_bayes_header)"
recoeff_vertex_scan_data="$(manifest_role_path ppg_recoeff_vertex_scan_data)"
recoeff_mbd_correction="$(manifest_role_path ppg_recoeff_mbd_correction)"

python3 - "$recoeff_source_macro" "$recoeff_compat_macro" \
  "$recoeff_compat_receipt" "$recoeff_macro" "$recoeff_yaml_cpp" \
  "$recoeff_mbd_correction" "$expected_estimator_revision" \
  "$expected_recoeff_source_sha256" \
  "$expected_recoeff_roounfold_compat_sha256" \
  "$recoeff_roounfold_compat_needle" \
  "$recoeff_roounfold_compat_replacement" \
  "$expected_roounfold_library_sha256" "$expected_roounfold_pcm_sha256" \
  "$expected_roounfold_header_tree_sha256" <<'PY'
from pathlib import Path
import hashlib
import json
import sys

(
    source_raw, compat_raw, receipt_raw, baseline_raw, yaml_cpp_raw,
    mbd_correction_raw, expected_revision,
    expected_source_hash, expected_compat_hash, needle, replacement,
    expected_library_hash, expected_pcm_hash, expected_header_tree_hash,
) = sys.argv[1:]
source = Path(source_raw).resolve()
compat = Path(compat_raw).resolve()
receipt = Path(receipt_raw).resolve()
baseline = Path(baseline_raw).resolve()
yaml_cpp = Path(yaml_cpp_raw).resolve()
mbd_correction = Path(mbd_correction_raw).resolve()

def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

if digest(source) != expected_source_hash:
    raise SystemExit("canonical RecoEff source differs from pinned digest")
if digest(compat) != expected_compat_hash:
    raise SystemExit("RooUnfold compatibility macro differs from pinned digest")
data = json.loads(receipt.read_text())
if data.get("schema_version") != 1:
    raise SystemExit("invalid RooUnfold compatibility receipt schema")
if data.get("transform") != "ppg12_recoeff_roounfold_constructor_compat_v1":
    raise SystemExit("unexpected RooUnfold compatibility transform")
if data.get("source_revision") != expected_revision:
    raise SystemExit("RooUnfold compatibility transform uses the wrong source revision")
if Path(str(data.get("input_path", ""))).resolve() != source:
    raise SystemExit("RooUnfold compatibility input path differs from source role")
if data.get("input_sha256") != expected_source_hash:
    raise SystemExit("RooUnfold compatibility input digest differs")
if Path(str(data.get("output_path", ""))).resolve() != compat:
    raise SystemExit("RooUnfold compatibility output path differs from macro role")
if data.get("output_sha256") != expected_compat_hash:
    raise SystemExit("RooUnfold compatibility output digest differs")
operation = data.get("operation", {})
if not isinstance(operation, dict):
    raise SystemExit("RooUnfold compatibility receipt lacks its operation")
if operation.get("label") != "remove_unsupported_explicit_false_constructor_argument":
    raise SystemExit("RooUnfold compatibility operation label differs")
if operation.get("expected_count") != 1 or operation.get("observed_count") != 1:
    raise SystemExit("RooUnfold compatibility operation is not exact-once")
if operation.get("needle") != needle or operation.get("replacement") != replacement:
    raise SystemExit("RooUnfold compatibility operation text differs")
derived = source.read_bytes().replace(needle.encode(), replacement.encode())
if source.read_bytes().count(needle.encode()) != 1 or derived != compat.read_bytes():
    raise SystemExit("RooUnfold compatibility macro cannot be re-derived exactly")
runtime = data.get("historical_roounfold_contract", {})
if not isinstance(runtime, dict):
    raise SystemExit("RooUnfold compatibility receipt lacks runtime contract")
for key, expected in (
    ("library_sha256", expected_library_hash),
    ("pcm_sha256", expected_pcm_hash),
    ("header_tree_sha256", expected_header_tree_hash),
):
    if runtime.get(key) != expected:
        raise SystemExit(f"RooUnfold compatibility {key} differs")
if runtime.get("runtime_smoke_requires_default_overflow_false") is not True:
    raise SystemExit("RooUnfold compatibility does not require overflow-default smoke")
for key in (
    "canonical_source_unchanged",
    "constructor_default_overflow_equals_explicit_false",
    "response_object_setup_only",
):
    if data.get(key) is not True:
        raise SystemExit(f"RooUnfold compatibility invariant is false: {key}")
for key in (
    "selection_or_fill_expression_replaced",
    "purity_estimator_expression_replaced",
):
    if data.get(key) is not False:
        raise SystemExit(f"RooUnfold compatibility altered forbidden semantics: {key}")

if compat.parent.name != "macros" or compat.parent.parent.name != "estimator":
    raise SystemExit("RooUnfold compatibility macro has an unexpected runtime layout")
runtime_root = compat.parent.parent.parent
expected_layout = {
    "baseline": compat.parent / "RecoEffCalculator_TTreeReader.C",
    "yaml_cpp": runtime_root / "lib" / "libyaml-cpp.so",
    "mbd_correction": runtime_root / "estimator" / "data" / "MbdOut.corr",
}
actual_layout = {
    "baseline": baseline,
    "yaml_cpp": yaml_cpp,
    "mbd_correction": mbd_correction,
}
for role, expected in expected_layout.items():
    if actual_layout[role].resolve() != expected.resolve():
        raise SystemExit(f"executable RecoEff runtime layout differs for {role}")
text = compat.read_text()
path_rewrites = (
    (
        "/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so",
        str(yaml_cpp),
        "yaml_cpp",
    ),
    (
        "/sphenix/user/shuhangli/ppg12/efficiencytool/MbdOut.corr",
        str(mbd_correction),
        "mbd_correction",
    ),
)
for old, new, label in path_rewrites:
    observed = text.count(old)
    if observed != 1:
        raise SystemExit(
            "executable RecoEff path rewrite is not exact-once: "
            f"role={label} observed={observed}"
        )
    text = text.replace(old, new)
if text.encode() != baseline.read_bytes():
    raise SystemExit(
        "executable RecoEff baseline cannot be re-derived from the pinned "
        "compatibility macro and two allowed path rewrites"
    )
PY
sealed_apply_bdt="$(manifest_role_path ppg_apply_bdt_macro)"
sealed_apply_config="$(manifest_role_path ppg_apply_bdt_config)"
sealed_base_e_model="$(manifest_role_path ppg_apply_model_base_E)"
sealed_base_v3e_model="$(manifest_role_path ppg_apply_model_base_v3E)"
sealed_npb_model="$(manifest_role_path ppg_apply_npb_model)"
recoeff_include_root="$(dirname "$recoeff_cross_section_header")"
validate_pinned_roounfold_library "$recoeff_roounfold"
[[ "$(basename "$recoeff_roounfold_pcm")" == RooUnfoldDict_rdict.pcm ]] || \
  die "sealed RooUnfold PCM has the wrong basename"
[[ "$(sha256_file "$recoeff_roounfold_pcm")" == "$expected_roounfold_pcm_sha256" ]] || \
  die "sealed RooUnfold PCM differs from pinned historical digest"
[[ "$(cd "$(dirname "$recoeff_roounfold")" && pwd -P)" == \
   "$(cd "$(dirname "$recoeff_roounfold_pcm")" && pwd -P)" ]] || \
  die "sealed RooUnfold library and PCM are not co-located"
recoeff_roounfold_include_root="$(python3 - "$recoeff_roounfold_header_receipt" \
  "$recoeff_roounfold_response_header" "$recoeff_roounfold_bayes_header" \
  "$expected_roounfold_header_tree_sha256" <<'PY'
from pathlib import Path
import hashlib
import json
import re
import sys

receipt, response_role, bayes_role = map(Path, sys.argv[1:4])
expected_tree_digest = sys.argv[4]
try:
    data = json.loads(receipt.read_text())
except (OSError, json.JSONDecodeError) as exc:
    raise SystemExit(f"cannot read RooUnfold header receipt: {exc}")
expected_names = (
    "RooUnfold.h",
    "RooUnfoldResponse.h",
    "RooUnfoldBayes.h",
    "RooUnfoldBinByBin.h",
    "RooUnfoldErrors.h",
    "RooUnfoldInvert.h",
    "RooUnfoldParms.h",
    "RooUnfoldSvd.h",
    "RooUnfoldTUnfold.h",
)
if data.get("schema_version") != 1 or data.get("role") != "ppg_recoeff_roounfold_header_tree":
    raise SystemExit("invalid RooUnfold header-tree receipt")
include_root = Path(str(data.get("include_root", "")))
rows = data.get("files")
if not include_root.is_absolute() or not include_root.is_dir() or not isinstance(rows, list):
    raise SystemExit("invalid RooUnfold header-tree root or inventory")
names = tuple(
    str(row.get("relative_path", "")) if isinstance(row, dict) else ""
    for row in rows
)
if names != expected_names or len(set(names)) != len(names):
    raise SystemExit("RooUnfold header inventory differs from exact nine-file contract")
tree_digest = hashlib.sha256()
paths = {}
for row, name in zip(rows, names):
    header = include_root / name
    digest = str(row.get("sha256", ""))
    if (
        not header.is_file()
        or header.stat().st_size <= 0
        or header.resolve().parent != include_root.resolve()
        or not re.fullmatch(r"[0-9a-f]{64}", digest)
    ):
        raise SystemExit(f"invalid sealed RooUnfold header metadata: {name}")
    observed = hashlib.sha256(header.read_bytes()).hexdigest()
    if observed != digest:
        raise SystemExit(f"sealed RooUnfold header drifted: {name}")
    paths[name] = header.resolve()
    tree_digest.update(name.encode())
    tree_digest.update(b"\0")
    tree_digest.update(bytes.fromhex(observed))
if data.get("tree_sha256") != expected_tree_digest:
    raise SystemExit("sealed RooUnfold header tree differs from pinned digest")
if tree_digest.hexdigest() != data.get("tree_sha256"):
    raise SystemExit("sealed RooUnfold header-tree digest drifted")
if response_role.resolve() != paths["RooUnfoldResponse.h"]:
    raise SystemExit("sealed RooUnfoldResponse role differs from header tree")
if bayes_role.resolve() != paths["RooUnfoldBayes.h"]:
    raise SystemExit("sealed RooUnfoldBayes role differs from header tree")
print(include_root.resolve())
PY
)" || die "failed to validate sealed RooUnfold header tree"

# The public driver keeps the three historical asset arguments for backward
# compatibility.  Execution is allowed only when they are byte-identical to
# the sealed runtime roles; the mutable paths are never executed directly.
cmp -s "$apply_bdt" "$sealed_apply_bdt" || \
  die "provided apply_BDT macro differs from sealed source-locked runtime"
cmp -s "$apply_config" "$sealed_apply_config" || \
  die "provided apply_BDT config differs from sealed source-locked runtime"
cmp -s "$base_e_model" "$sealed_base_e_model" || \
  die "provided base_E model differs from sealed split model"
cmp -s "$base_v3e_model" "$sealed_base_v3e_model" || \
  die "provided base_v3E model differs from sealed split model"
cmp -s "$npb_model" "$sealed_npb_model" || \
  die "provided NPB model differs from sealed split model"
apply_bdt="$sealed_apply_bdt"
apply_config="$sealed_apply_config"
base_e_model="$sealed_base_e_model"
base_v3e_model="$sealed_base_v3e_model"
npb_model="$sealed_npb_model"

# The new.17 setup does not place OFFLINE_MAIN/rootmacros on ROOT's default
# macro path.  Resolve the release-owned calibration macro explicitly and add
# that exact directory to both Cling's include path and ROOT's macro path for
# each executable below.  Dynamic discovery would otherwise fail before the
# oracle runs despite the builder having already sealed this same file.
calo_macro_dir="${OFFLINE_MAIN}/rootmacros"
calo_calib="${calo_macro_dir}/Calo_Calib.C"
[[ -f "$calo_calib" && -s "$calo_calib" ]] || \
  die "exact new.17 Calo_Calib.C is missing: $calo_calib"

# CaloAna24's executable output contract is the cwd-local `caloana.root`;
# the Fun4All output argument controls the DST output, not the slimtree file.
# Bind the worker to that real executable filename instead of assuming the
# wrapper argument renames it.
ppg_raw_root="${ppg_dir}/caloana.root"
ppg_scored_root="${ppg_dir}/caloana_with_bdt_split.root"
recoil_root="${recoil_dir}/recoil.root"
ppg_log="${ppg_dir}/ppg12.log"
recoil_log="${recoil_dir}/recoiljets.log"
apply_log="${ppg_dir}/apply_bdt.log"
apply_evidence="${ppg_dir}/apply_bdt_stage_evidence.json"
apply_runtime_config="${ppg_dir}/config_nom_split_oracle.yaml"
baseline_recoeff_config="${baseline_layout}/config_bdt_nom_oracle.yaml"
trace_recoeff_config="${trace_layout}/config_bdt_nom_oracle.yaml"
baseline_recoeff_log="${baseline_layout}/recoeff.log"
trace_recoeff_log="${trace_layout}/recoeff.log"
baseline_recoeff_scan_log="${baseline_layout}/vertex_scan.log"
trace_recoeff_scan_log="${trace_layout}/vertex_scan.log"
baseline_eff_root="${baseline_layout}/output/MC_efficiency_${oracle_sample}_paired_oracle.root"
trace_eff_root="${trace_layout}/output/MC_efficiency_${oracle_sample}_paired_oracle.root"
baseline_response_root="${baseline_layout}/output/MC_response_${oracle_sample}_paired_oracle.root"
trace_response_root="${trace_layout}/output/MC_response_${oracle_sample}_paired_oracle.root"
baseline_vtxscan_root="${baseline_layout}/output/MC_efficiency_${oracle_sample}_paired_oracle_vtxscan.root"
trace_vtxscan_root="${trace_layout}/output/MC_efficiency_${oracle_sample}_paired_oracle_vtxscan.root"
ppg_candidate_trace="${trace_layout}/candidate_trace.csv"
ppg_response_trace="${trace_layout}/response_trace.csv"
aggregate_report="${report_dir}/executable_aggregate.json"
root_equivalence_report="${report_dir}/recoeff_root_equivalence.json"
report_md="${report_dir}/paired_oracle_report.md"
summary_csv="${report_dir}/paired_oracle_summary.csv"
candidate_csv="${report_dir}/paired_oracle_candidates.csv"
contract="${output_dir}/paired_oracle_contract.json"

failure_stage="estimator_config_materialization"
python3 - "$apply_config" "$apply_runtime_config" <<'PY'
from pathlib import Path
import sys

source, output = map(Path, sys.argv[1:])
text = source.read_text()
block = '    - node: "CLUSTERINFO_CEMC_NO_SPLIT"\n      model_suffix: "_nosplit"\n'
if text.count(block) != 1:
    raise SystemExit("canonical apply_BDT config no-split block changed")
text = text.replace(block, "")
if text.count('model_suffix: "_split"') != 1 or "_nosplit" in text:
    raise SystemExit("derived apply_BDT config is not the single canonical split lane")
output.write_text(text)
PY

for layout in "$baseline_layout" "$trace_layout"; do
  ln -s "$ppg_scored_root" "${layout}/input/${oracle_sample}/bdt_split.root"
  ln -s "$recoeff_vertex_scan_data" "${layout}/input/data_histo_bdt_nom_vtxscan.root"
  ln -s "$tower_mask" "${layout}/input/tower_masks_bdt_nom.root"
  ln -s "$recoeff_truth_vertex_reweight" \
    "${layout}/input/truth_vertex_reweight.root"
done

# Produce two byte-identical, relocatable estimator configs.  Their only
# changes from the source-locked period file are sealed I/O paths and the
# unique output suffix; all scientific analysis keys remain byte-for-byte.
python3 - "$recoeff_period_config" "$baseline_recoeff_config" \
  "$trace_recoeff_config" "$recoeff_var_type_suffix" "$period" <<'PY'
from pathlib import Path
import sys

source, baseline, trace = map(Path, sys.argv[1:4])
config_suffix = sys.argv[4]
period = sys.argv[5]
text = source.read_text()
truth_vertex_source = (
    "/sphenix/user/shuhangli/ppg12/efficiencytool/"
    f"truth_vertex_reweight/output/{period}/reweight.root"
)
rewrites = (
    ('photon_jet_file_root_dir: "/sphenix/user/shuhangli/ppg12/FunWithxgboost/"',
     'photon_jet_file_root_dir: "input/"'),
    ('eff_outfile: "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency"',
     'eff_outfile: "output/MC_efficiency"'),
    ('response_outfile: "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_response"',
     'response_outfile: "output/MC_response"'),
    ('data_outfile: "/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo"',
     'data_outfile: "output/data_histo"'),
    (f'var_type: "bdt_nom_{config_suffix}"', 'var_type: "paired_oracle"'),
    ('vertex_scan_data_file: ""',
     'vertex_scan_data_file: "input/data_histo_bdt_nom_vtxscan.root"'),
    ('tower_mask_file: "/sphenix/user/shuhangli/ppg12/efficiencytool/tower_masks_bdt_nom.root"',
     'tower_mask_file: "input/tower_masks_bdt_nom.root"'),
    (f'truth_vertex_reweight_file: "{truth_vertex_source}"',
     'truth_vertex_reweight_file: "input/truth_vertex_reweight.root"'),
)
for old, new in rewrites:
    observed = text.count(old)
    if observed != 1:
        raise SystemExit(f"estimator config rewrite count for {old!r}: {observed}")
    text = text.replace(old, new)
baseline.write_text(text)
trace.write_text(text)
if baseline.read_bytes() != trace.read_bytes():
    raise SystemExit("baseline and trace estimator configs differ")
PY

failure_stage="contract_materialization"
python3 - \
  "$contract" "$lane_id" "$sample" "$period" "$interaction" \
  "$setup_script" "$calo_calib" "$ppg_macro" "$ppg_lib" \
  "$g4_full_list" "$truthjet_full_list" "$g4_slice" "$truthjet_slice" \
  "$combined_list" "$source_pair_receipt" "$apply_bdt" "$apply_config" "$base_e_model" \
  "$base_v3e_model" "$npb_model" "$tower_mask" \
  "$recoil_runtime_manifest" "$paired_runtime_manifest" "$recoil_macro" "$recoil_config" \
  "$ppg_wrapper" "$recoil_wrapper" "$comparator" "$auditor" "$ppg_log" \
  "$recoil_log" "$ppg_raw_root" "$ppg_scored_root" "$recoil_root" \
  "$candidate_csv" "$expected_token" "$first_ph_seed" \
  "$pedestal_seed" "$pedestal_sequence" "$ph_seed_sequence" \
  "$recoeff_source_macro" "$recoeff_compat_macro" "$recoeff_compat_receipt" \
  "$recoeff_macro" "$recoeff_trace_macro" "$recoeff_trace_receipt" \
  "$recoeff_canonical_config" "$recoeff_period_config" \
  "$recoeff_truth_vertex_reweight" "$recoeff_yaml_cpp_header_receipt" \
  "$baseline_recoeff_config" "$trace_recoeff_config" \
  "$recoeff_cross_section_header" "$recoeff_truth_vertex_header" \
  "$recoeff_yaml_cpp" "$recoeff_roounfold" "$recoeff_roounfold_pcm" \
  "$recoeff_roounfold_header_receipt" \
  "$recoeff_roounfold_response_header" "$recoeff_roounfold_bayes_header" \
  "$recoeff_vertex_scan_data" "$recoeff_mbd_correction" \
  "$baseline_recoeff_log" "$trace_recoeff_log" \
  "$baseline_recoeff_scan_log" "$trace_recoeff_scan_log" \
  "$baseline_eff_root" "$trace_eff_root" \
  "$baseline_response_root" "$trace_response_root" \
  "$baseline_vtxscan_root" "$trace_vtxscan_root" \
  "$ppg_candidate_trace" "$ppg_response_trace" "$aggregate_report" \
  "$aggregate_extractor" "$apply_evidence" "$apply_runtime_config" \
  "$driver_script" "$worker_self" "$reuse_ppg_raw_root" \
  "$reuse_ppg_raw_contract" "$reuse_mode" <<'PY'
from pathlib import Path
import hashlib
import json
import sys

(
    contract, lane_id, sample, period, interaction,
    setup, calo, ppg_macro, ppg_lib, g4_full, truth_full, g4_slice,
    truth_slice, combined, source_pair_receipt, apply_bdt, apply_config, base_e,
    base_v3e, npb, mask, recoil_manifest, paired_runtime_manifest, recoil_macro,
    recoil_config, ppg_wrapper,
    recoil_wrapper, comparator, auditor, ppg_log, recoil_log, ppg_raw, ppg_scored,
    recoil_root, candidate_csv, token, first_ph_seed,
    pedestal_seed, pedestal_sequence, ph_seed_sequence,
    recoeff_source, recoeff_compat_macro, recoeff_compat_receipt,
    recoeff_macro, recoeff_trace_macro, recoeff_trace_receipt,
    recoeff_canonical_config, recoeff_period_config, recoeff_truth_vertex_reweight,
    recoeff_yaml_cpp_header_receipt, baseline_recoeff_config, trace_recoeff_config,
    cross_section_header, truth_vertex_header, yaml_cpp, roounfold, roounfold_pcm,
    roounfold_header_receipt,
    roounfold_response_header, roounfold_bayes_header, vertex_scan_data,
    mbd_correction, baseline_recoeff_log, trace_recoeff_log,
    baseline_scan_log, trace_scan_log, baseline_eff_root, trace_eff_root,
    baseline_response_root, trace_response_root, baseline_vtxscan_root,
    trace_vtxscan_root, ppg_candidate_trace, ppg_response_trace,
    aggregate_report, aggregate_extractor, apply_evidence, apply_runtime_config,
    driver_script, worker_script, reuse_raw_root, reuse_raw_contract, reuse_mode,
) = sys.argv[1:]

def digest(path: str) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

paths = {
    "setup_script": setup,
    "calo_calib": calo,
    "ppg_macro": ppg_macro,
    "ppg_caloana24": ppg_lib,
    "g4_full_list": g4_full,
    "truthjet_full_list": truth_full,
    "g4_slice": g4_slice,
    "truthjet_slice": truth_slice,
    "combined_list": combined,
    "source_pair_receipt": source_pair_receipt,
    "apply_bdt": apply_bdt,
    "apply_config": apply_config,
    "base_e_model": base_e,
    "base_v3e_model": base_v3e,
    "npb_model": npb,
    "tower_mask": mask,
    "recoil_runtime_manifest": recoil_manifest,
    "paired_runtime_manifest": paired_runtime_manifest,
    "recoil_macro": recoil_macro,
    "recoil_config": recoil_config,
    "ppg_wrapper": ppg_wrapper,
    "recoil_wrapper": recoil_wrapper,
    "comparator": comparator,
    "auditor": auditor,
    "ppg_log": ppg_log,
    "recoil_log": recoil_log,
    "ppg_raw_root": ppg_raw,
    "ppg_scored_root": ppg_scored,
    "recoil_root": recoil_root,
    "candidate_csv": candidate_csv,
    "ppg_recoeff_source_macro": recoeff_source,
    "ppg_recoeff_roounfold_compat_macro": recoeff_compat_macro,
    "ppg_recoeff_roounfold_compat_transform_receipt": recoeff_compat_receipt,
    "ppg_recoeff_macro": recoeff_macro,
    "ppg_recoeff_trace_macro": recoeff_trace_macro,
    "ppg_recoeff_trace_transform_receipt": recoeff_trace_receipt,
    "ppg_recoeff_canonical_config": recoeff_canonical_config,
    "ppg_recoeff_period_config": recoeff_period_config,
    "ppg_recoeff_truth_vertex_reweight": recoeff_truth_vertex_reweight,
    "ppg_recoeff_yaml_cpp_header_tree_receipt": recoeff_yaml_cpp_header_receipt,
    "ppg_recoeff_baseline_config": baseline_recoeff_config,
    "ppg_recoeff_trace_config": trace_recoeff_config,
    "ppg_recoeff_cross_section_header": cross_section_header,
    "ppg_recoeff_truth_vertex_header": truth_vertex_header,
    "ppg_recoeff_yaml_cpp": yaml_cpp,
    "ppg_recoeff_roounfold": roounfold,
    "ppg_recoeff_roounfold_pcm": roounfold_pcm,
    "ppg_recoeff_roounfold_header_tree_receipt": roounfold_header_receipt,
    "ppg_recoeff_roounfold_response_header": roounfold_response_header,
    "ppg_recoeff_roounfold_bayes_header": roounfold_bayes_header,
    "ppg_recoeff_vertex_scan_data": vertex_scan_data,
    "ppg_recoeff_mbd_correction": mbd_correction,
    "ppg_recoeff_baseline_log": baseline_recoeff_log,
    "ppg_recoeff_trace_log": trace_recoeff_log,
    "ppg_recoeff_baseline_scan_log": baseline_scan_log,
    "ppg_recoeff_trace_scan_log": trace_scan_log,
    "ppg_recoeff_baseline_eff_root": baseline_eff_root,
    "ppg_recoeff_trace_eff_root": trace_eff_root,
    "ppg_recoeff_baseline_response_root": baseline_response_root,
    "ppg_recoeff_trace_response_root": trace_response_root,
    "ppg_recoeff_baseline_vtxscan_root": baseline_vtxscan_root,
    "ppg_recoeff_trace_vtxscan_root": trace_vtxscan_root,
    "ppg12_executable_trace": ppg_candidate_trace,
    "ppg12_executable_response_trace": ppg_response_trace,
    "ppg12_executable_aggregate": aggregate_report,
    "ppg12_recoeff_root_equivalence": str(Path(aggregate_report).with_name("recoeff_root_equivalence.json")),
    "aggregate_extractor": aggregate_extractor,
    "apply_bdt_stage_evidence": apply_evidence,
    "apply_bdt_runtime_config": apply_runtime_config,
    "first_divergence_json": str(Path(candidate_csv).with_name("first_divergence.json")),
}
token_bound_paths = {
    "setup_script": setup,
    "ppg_macro": ppg_macro,
    "g4_full_list": g4_full,
    "truthjet_full_list": truth_full,
    "apply_bdt": apply_bdt,
    "apply_config": apply_config,
    "base_e_model": base_e,
    "base_v3e_model": base_v3e,
    "npb_model": npb,
    "tower_mask": mask,
    "recoil_runtime_manifest": recoil_manifest,
    "recoil_config": recoil_config,
    "driver_script": driver_script,
    "worker_script": worker_script,
    "ppg_wrapper": ppg_wrapper,
    "recoil_wrapper": recoil_wrapper,
    "comparator": comparator,
    "auditor": auditor,
    "aggregate_extractor": aggregate_extractor,
    "ppg_recoeff_period_config": recoeff_period_config,
    "ppg_recoeff_truth_vertex_reweight": recoeff_truth_vertex_reweight,
    "ppg_recoeff_yaml_cpp_header_tree_receipt": recoeff_yaml_cpp_header_receipt,
}
if reuse_mode == "exact_contract_bound":
    token_bound_paths.update({
        "reuse_ppg_raw_root": reuse_raw_root,
        "reuse_ppg_raw_contract": reuse_raw_contract,
    })
lane_data = {
    "lane_id": lane_id, "sample": sample, "period": period,
    "interaction": interaction, "rows": 5,
}
rng_data = {
    "mode": "historical_fifo_replay_v2",
    "reco_consts_randomseed": "absent",
    "ph_seed_sequence": [int(value) for value in ph_seed_sequence.split(",")],
    "first_ph_seed": int(first_ph_seed),
    "pedestal_seed": int(pedestal_seed),
    "pedestal_sequence": int(pedestal_sequence),
    "pedestal_file": f"pedestal-54256-0{int(pedestal_sequence):04d}.root",
    "source_log": "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/photon5/condorout/OutDir0/test.out",
    "source_log_call_count": 5,
}
runtime_data = {
    "profile": "new.17",
    "offline_main": "/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_new/new.17",
    "actual_offline_main": "/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_new/new.17",
}
source_pairs_data = json.loads(Path(source_pair_receipt).read_text())
raw_reconstruction_roles = (
    "setup_script", "calo_calib", "ppg_macro", "ppg_caloana24",
    "ppg_wrapper", "g4_full_list", "truthjet_full_list", "g4_slice",
    "truthjet_slice", "source_pair_receipt",
)
raw_reconstruction_dependencies = {
    "schema_version": 1,
    "lane": lane_data,
    "rng": rng_data,
    "runtime": runtime_data,
    "source_pairs": source_pairs_data,
    "files": {
        role: {"path": paths[role], "sha256": digest(paths[role])}
        for role in raw_reconstruction_roles
    },
}
data = {
    "schema_version": 3,
    "authorization_token": token,
    "authorization": {
        "token_schema_version": 3,
        "token_bound_files": {
            role: {"path": path, "sha256": digest(path)}
            for role, path in token_bound_paths.items()
        },
        "raw_ppg12_reuse": {
            "mode": reuse_mode,
            "source": None if reuse_mode == "disabled" else {
                "path": reuse_raw_root,
                "sha256": digest(reuse_raw_root),
            },
            "contract": None if reuse_mode == "disabled" else {
                "path": reuse_raw_contract,
                "sha256": digest(reuse_raw_contract),
            },
        },
    },
    "lane": lane_data,
    "rng": rng_data,
    "runtime": runtime_data,
    "runtime_hashes": {
        role: digest(paths[role])
        for role in (
            "setup_script", "calo_calib", "recoil_config", "ppg_wrapper",
            "recoil_wrapper", "comparator",
            "auditor", "aggregate_extractor",
            "source_pair_receipt", "paired_runtime_manifest",
            "ppg_recoeff_period_config", "ppg_recoeff_truth_vertex_reweight",
            "ppg_recoeff_yaml_cpp_header_tree_receipt",
            "ppg_recoeff_baseline_config", "ppg_recoeff_trace_config",
            "apply_bdt_runtime_config",
        )
    },
    "source_pairs": source_pairs_data,
    "raw_reconstruction_dependencies": raw_reconstruction_dependencies,
    "paths": paths,
}
Path(contract).write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
PY

failure_stage="contract_preflight_audit"
python3 "$auditor" preflight --contract "$contract"

# Use symlink-only library views so neither side inherits a broad user install.
ppg_lib_view="${runtime_dir}/ppg_lib"
recoil_lib_view="${runtime_dir}/recoil_lib"
mkdir -p "$ppg_lib_view" "$recoil_lib_view"
ln -s "$ppg_lib" "${ppg_lib_view}/libCaloAna24.so"
ln -s "$calo_reco_lib" "${recoil_lib_view}/libcalo_reco.so"
ln -s "$recoil_lib" "${recoil_lib_view}/libRecoilJets.so"
ln -s "$clusteriso_lib" "${recoil_lib_view}/libclusteriso.so"
ln -s "$jetbase_lib" "${recoil_lib_view}/libjetbase.so"

if [[ "$reuse_mode" == exact_contract_bound ]]; then
  failure_stage="ppg12_raw_reuse_validation"
  reuse_validation_mode="$(python3 - "$reuse_ppg_raw_contract" "$reuse_ppg_raw_root" \
    "$g4_slice" "$truthjet_slice" "$ppg_macro" "$ppg_wrapper" \
    "$setup_script" "$g4_full_list" "$truthjet_full_list" "$lane_id" \
    "$sample" "$period" "$interaction" "$ppg_lib" "$calo_calib" \
    "$source_pair_receipt" <<'PY'
from pathlib import Path
import hashlib
import json
import re
import sys

contract_path, raw_root, g4_slice, truth_slice = map(Path, sys.argv[1:5])
ppg_macro, ppg_wrapper, setup_script, g4_full, truth_full = map(
    Path, sys.argv[5:10]
)
lane_id, sample, period, interaction = sys.argv[10:14]
ppg_lib, calo_calib, source_pair_receipt = map(Path, sys.argv[14:17])
data = json.loads(contract_path.read_text())
authorization = data.get("authorization", {})
if authorization.get("token_schema_version") != 3:
    raise SystemExit(
        "reuse contract predates content-bound schema 3 and is scientifically inadmissible"
    )
token_files = authorization.get("token_bound_files", {})
if not isinstance(token_files, dict):
    raise SystemExit("reuse contract lacks token-bound file metadata")

# Recompute the *prior* plan/run authorization token before trusting any of
# the receipt fields below.  A schema-version marker plus self-consistent file
# hashes is not an authorization receipt: without this replay a hand-authored
# JSON wrapper could make an arbitrary raw ROOT look like an accepted run.
required_path_roles = [
    "output_dir", "setup_script", "ppg_macro", "g4_full_list",
    "truthjet_full_list", "apply_bdt", "apply_config", "base_e_model",
    "base_v3e_model", "npb_model", "tower_mask",
    "recoil_runtime_manifest", "recoil_config",
]
token_file_roles = [
    "setup_script", "ppg_macro", "g4_full_list", "truthjet_full_list",
    "apply_bdt", "apply_config", "base_e_model", "base_v3e_model",
    "npb_model", "tower_mask", "recoil_runtime_manifest", "recoil_config",
    "driver_script", "worker_script", "ppg_wrapper", "recoil_wrapper",
    "comparator", "auditor", "aggregate_extractor",
]
prior_paths = data.get("paths", {})
has_period_path = bool(prior_paths.get("ppg_recoeff_period_config"))
has_period_token = "ppg_recoeff_period_config" in token_files
if has_period_path != has_period_token:
    raise SystemExit("reuse contract has an incomplete period-config token binding")
if has_period_path:
    required_path_roles.append("ppg_recoeff_period_config")
    token_file_roles.append("ppg_recoeff_period_config")
has_vertex_path = bool(prior_paths.get("ppg_recoeff_truth_vertex_reweight"))
has_vertex_token = "ppg_recoeff_truth_vertex_reweight" in token_files
if has_vertex_path != has_vertex_token:
    raise SystemExit("reuse contract has an incomplete truth-vertex-weight token binding")
if has_vertex_path:
    required_path_roles.append("ppg_recoeff_truth_vertex_reweight")
    token_file_roles.append("ppg_recoeff_truth_vertex_reweight")
has_yaml_headers_path = bool(
    prior_paths.get("ppg_recoeff_yaml_cpp_header_tree_receipt")
)
has_yaml_headers_token = (
    "ppg_recoeff_yaml_cpp_header_tree_receipt" in token_files
)
if has_yaml_headers_path != has_yaml_headers_token:
    raise SystemExit("reuse contract has an incomplete yaml-cpp-header token binding")
if has_yaml_headers_path:
    required_path_roles.append("ppg_recoeff_yaml_cpp_header_tree_receipt")
    token_file_roles.append("ppg_recoeff_yaml_cpp_header_tree_receipt")
prior_reuse = authorization.get("raw_ppg12_reuse", {})
prior_reuse_mode = prior_reuse.get("mode")
if prior_reuse_mode not in {"disabled", "exact_contract_bound"}:
    raise SystemExit("reuse contract records an invalid prior raw-reuse mode")
if prior_reuse_mode == "exact_contract_bound":
    token_file_roles.extend(["reuse_ppg_raw_root", "reuse_ppg_raw_contract"])
if set(token_files) != set(token_file_roles):
    raise SystemExit("reuse contract token-bound role set is not a recognized schema-3 variant")

prior_output_dir = str(contract_path.parent.resolve())
required_values = {"output_dir": prior_output_dir}
for role in required_path_roles[1:]:
    value = str(prior_paths.get(role, ""))
    if not value:
        raise SystemExit(f"reuse contract lacks token-bound path role {role}")
    required_values[role] = value

token_lines = [
    "schema_version=3",
    f"lane={sample}:{period}:{interaction}",
    f"lane_id={lane_id}",
    "rows=5",
    "rng_mode=historical_fifo_replay_v2",
    "first_ph_seed=2991264730",
    "pedestal_seed=4256268992",
    "pedestal=534",
    "ph_seed_sequence=2991264730,4256268992,2394322166,874466025,2240380304",
    "source_graph=NONE,g4,truthjet,NONE,NONE",
    "runtime=new.17",
    f"reuse_mode={prior_reuse_mode}",
]
token_lines.extend(f"{role}={required_values[role]}" for role in required_path_roles)
for role in token_file_roles:
    item = token_files.get(role)
    if not isinstance(item, dict):
        raise SystemExit(f"reuse contract lacks token-bound {role}")
    path = str(item.get("path", ""))
    recorded = str(item.get("sha256", ""))
    if not path or len(recorded) != 64:
        raise SystemExit(f"reuse contract has malformed token-bound {role}")
    token_lines.append(f"{role}={path} sha256={recorded}")
token_payload = "".join(f"{line}\n" for line in sorted(token_lines)).encode()
recomputed_token = "ppg12-oracle:" + hashlib.sha256(token_payload).hexdigest()
if data.get("authorization_token") != recomputed_token:
    raise SystemExit("reuse contract authorization token does not replay exactly")

lane = data.get("lane", {})
if lane != {
    "lane_id": lane_id, "sample": sample, "period": period,
    "interaction": interaction, "rows": 5,
}:
    raise SystemExit("reuse contract lane differs from the exact paired oracle lane")
rng = data.get("rng", {})
if rng.get("mode") != "historical_fifo_replay_v2" or rng.get("ph_seed_sequence") != [
    2991264730, 4256268992, 2394322166, 874466025, 2240380304,
]:
    raise SystemExit("reuse contract RNG differs from the historical FIFO replay")
runtime = data.get("runtime", {})
if runtime.get("profile") != "new.17" or runtime.get("actual_offline_main") != (
    "/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_new/new.17"
):
    raise SystemExit("reuse contract runtime differs from exact new.17")
paths = data.get("paths", {})

run_state = contract_path.parent / "RUN_STATE"
if not run_state.is_file():
    raise SystemExit("reuse contract has no RUN_STATE")
run_state_value = run_state.read_text().strip()
legacy_failed_apply = run_state_value == "FAILED"
if run_state_value not in {"PASS", "FAILED"}:
    raise SystemExit("reuse contract has an inadmissible RUN_STATE")
if legacy_failed_apply:
    expected_legacy_sha256 = (
        "32adf2501ba593c24f380c55e89b072a07dbfea45391d5ed4c796c0bafe358ec"
    )
    if raw_root.stat().st_size != 7701676 or hashlib.sha256(
        raw_root.read_bytes()
    ).hexdigest() != expected_legacy_sha256:
        raise SystemExit("FAILED-run reuse is limited to the known attempt6 raw ROOT")
    ppg_log = Path(str(paths.get("ppg_log", "")))
    if not ppg_log.is_file() or re.search(
        r"(?im)^.*(?:Fun4AllServer::run|processed|event).*\b5000\b.*$",
        ppg_log.read_text(errors="replace"),
    ) is None:
        raise SystemExit("known FAILED-run reuse lacks the 5000-event completion log")
    apply_log = contract_path.parent / "ppg12" / "apply_bdt.log"
    apply_text = apply_log.read_text(errors="replace") if apply_log.is_file() else ""
    if "yaml-cpp/yaml.h" not in apply_text or not any(
        marker in apply_text for marker in ("file not found", "No such file or directory")
    ):
        raise SystemExit("FAILED-run reuse is not the known missing-yaml-header failure")
    if Path(str(paths.get("ppg_scored_root", ""))).exists() or Path(
        str(paths.get("apply_bdt_stage_evidence", ""))
    ).exists():
        raise SystemExit("FAILED-run reuse unexpectedly contains downstream apply_BDT output")

def rows(path: Path) -> list[str]:
    return [
        line.strip() for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]

for key, current in (("g4_slice", g4_slice), ("truthjet_slice", truth_slice)):
    previous = Path(str(paths.get(key, "")))
    if not previous.is_file() or rows(previous) != rows(current):
        raise SystemExit(f"reuse contract {key} differs from current five-file source slice")

def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def recorded_token_digest(key: str) -> str:
    item = token_files.get(key, {})
    if not isinstance(item, dict):
        raise SystemExit(f"reuse contract lacks token-bound {key}")
    previous = Path(str(item.get("path", "")))
    recorded = str(item.get("sha256", ""))
    if not previous.is_file() or digest(previous) != recorded:
        raise SystemExit(f"reuse contract recorded {key} no longer matches its source")
    return recorded

# Only inputs that can affect the raw CaloAna24 reconstruction are compared
# here.  Estimator configs, apply_BDT assets, RecoilJets code, and this worker
# may change without forcing the expensive raw PPG12 stage to run again.
for key, current in (
    ("setup_script", setup_script),
    ("g4_full_list", g4_full),
    ("truthjet_full_list", truth_full),
    ("ppg_macro", ppg_macro),
    ("ppg_wrapper", ppg_wrapper),
):
    recorded = recorded_token_digest(key)
    if not current.is_file() or digest(current) != recorded:
        raise SystemExit(f"reuse contract {key} differs from current reconstruction asset")

runtime_hashes = data.get("runtime_hashes", {})
prior_calo = Path(str(paths.get("calo_calib", "")))
prior_calo_hash = str(runtime_hashes.get("calo_calib", ""))
if (
    not prior_calo.is_file()
    or digest(prior_calo) != prior_calo_hash
    or not calo_calib.is_file()
    or digest(calo_calib) != prior_calo_hash
):
    raise SystemExit("reuse contract calo_calib differs from current reconstruction asset")

prior_runtime_manifest = Path(str(paths.get("recoil_runtime_manifest", "")))
prior_runtime_manifest_hash = recorded_token_digest("recoil_runtime_manifest")
if digest(prior_runtime_manifest) != prior_runtime_manifest_hash:
    raise SystemExit("reuse contract runtime manifest no longer matches its token")
manifest_data = json.loads(prior_runtime_manifest.read_text())
ppg_rows = [
    row for row in manifest_data.get("files", [])
    if isinstance(row, dict) and row.get("role") == "libCaloAna24.so"
]
if len(ppg_rows) != 1:
    raise SystemExit("reuse runtime manifest lacks one source-locked libCaloAna24.so")
prior_ppg_lib = Path(str(ppg_rows[0].get("path", "")))
prior_ppg_lib_hash = str(ppg_rows[0].get("sha256", ""))
if (
    not prior_ppg_lib.is_file()
    or digest(prior_ppg_lib) != prior_ppg_lib_hash
    or not ppg_lib.is_file()
    or digest(ppg_lib) != prior_ppg_lib_hash
):
    raise SystemExit("reuse contract libCaloAna24.so differs from current reconstruction asset")

def source_pair_evidence(g4_rows: list[str], truth_rows: list[str]) -> dict:
    import re

    if len(g4_rows) != 5 or len(truth_rows) != 5:
        raise SystemExit("raw reconstruction reuse requires exactly five source pairs")
    sample_number = sample.removeprefix("Photon")
    stem = (
        f"PhotonJet{sample_number}"
        if interaction == "SI"
        else f"PhotonJet{sample_number}_pythia8_Detroit"
    )
    pattern = re.compile(rf"{re.escape(stem)}-\d{{10}}-\d{{6}}\.root")
    identities = []
    canonical_rows = []
    for index, (g4_value, truth_value) in enumerate(
        zip(g4_rows, truth_rows), start=1
    ):
        g4_name = Path(g4_value).name
        truth_name = Path(truth_value).name
        g4_prefix = "G4Hits_pythia8_"
        truth_prefix = "DST_TRUTH_JET_pythia8_"
        if not g4_name.startswith(g4_prefix) or not truth_name.startswith(truth_prefix):
            raise SystemExit(f"reuse source pair {index} has an unrecognized prefix")
        g4_identity = g4_name[len(g4_prefix):]
        truth_identity = truth_name[len(truth_prefix):]
        if (
            g4_identity != truth_identity
            or pattern.fullmatch(g4_identity) is None
        ):
            raise SystemExit(f"reuse source pair {index} has a cross-mode identity mismatch")
        identities.append(g4_identity)
        canonical_rows.append(f"{g4_value}\t{truth_value}\t{g4_identity}")
    if len(set(identities)) != 5:
        raise SystemExit("reuse source pairs are not unique")
    return {
        "schema_version": 1,
        "normalization": {
            "g4_prefix": "G4Hits_pythia8_",
            "truthjet_prefix": "DST_TRUTH_JET_pythia8_",
            "basename_only": True,
        },
        "sample": sample,
        "interaction": interaction,
        "row_count": 5,
        "identities": identities,
        "pair_identity_sha256": hashlib.sha256(
            ("\n".join(canonical_rows) + "\n").encode()
        ).hexdigest(),
    }

prior_source_pairs = source_pair_evidence(rows(Path(str(paths["g4_slice"]))), rows(Path(str(paths["truthjet_slice"]))))
current_source_pairs = source_pair_evidence(rows(g4_slice), rows(truth_slice))
if prior_source_pairs != current_source_pairs:
    raise SystemExit("reuse source-pair identities differ from current reconstruction inputs")
if json.loads(source_pair_receipt.read_text()) != current_source_pairs:
    raise SystemExit("current source-pair receipt differs from exact reconstruction inputs")

prior_raw_dependencies = data.get("raw_reconstruction_dependencies")
if prior_raw_dependencies is not None:
    if prior_raw_dependencies.get("schema_version") != 1:
        raise SystemExit("reuse contract raw-reconstruction dependency schema differs")
    if prior_raw_dependencies.get("lane") != lane:
        raise SystemExit("reuse contract raw-reconstruction lane differs")
    if prior_raw_dependencies.get("rng") != rng:
        raise SystemExit("reuse contract raw-reconstruction RNG differs")
    if prior_raw_dependencies.get("runtime") != runtime:
        raise SystemExit("reuse contract raw-reconstruction runtime differs")
    if prior_raw_dependencies.get("source_pairs") != prior_source_pairs:
        raise SystemExit("reuse contract raw-reconstruction source receipt differs")
    dependency_files = prior_raw_dependencies.get("files", {})
    expected_dependency_roles = {
        "setup_script", "calo_calib", "ppg_macro", "ppg_caloana24",
        "ppg_wrapper", "g4_full_list", "truthjet_full_list", "g4_slice",
        "truthjet_slice", "source_pair_receipt",
    }
    if set(dependency_files) != expected_dependency_roles:
        raise SystemExit("reuse contract raw-reconstruction role set differs")
    for role, metadata in dependency_files.items():
        prior_path = Path(str(metadata.get("path", "")))
        recorded_hash = str(metadata.get("sha256", ""))
        if not prior_path.is_file() or digest(prior_path) != recorded_hash:
            raise SystemExit(
                f"reuse contract raw-reconstruction dependency changed: {role}"
            )
if raw_root.name != "caloana.root":
    raise SystemExit("reuse source is not the executable CaloAna24 caloana.root artifact")
recorded_raw = Path(str(paths.get("ppg_raw_root", "")))
if recorded_raw.resolve() != raw_root.resolve():
    raise SystemExit("reuse source path differs from the prior schema-3 contract")
apply_evidence = Path(str(paths.get("apply_bdt_stage_evidence", "")))
if not legacy_failed_apply:
    if not apply_evidence.is_file():
        raise SystemExit("reuse contract lacks completed apply_BDT evidence for the raw input")
    apply_data = json.loads(apply_evidence.read_text())
    raw_input = apply_data.get("input", {})
    if Path(str(raw_input.get("path", ""))).resolve() != raw_root.resolve():
        raise SystemExit("reuse source differs from the prior apply_BDT input")
    if str(raw_input.get("sha256", "")) != digest(raw_root):
        raise SystemExit("reuse source hash differs from the prior apply_BDT input hash")
print("legacy_failed_apply" if legacy_failed_apply else "completed_pass")
PY
  )" || die "raw PPG12 reuse contract validation failed"
  reuse_root_audit="${runtime_dir}/audit_reused_ppg12_raw.C"
  python3 - "$reuse_root_audit" "$reuse_ppg_raw_root" \
    "$reuse_validation_mode" <<'PY'
from pathlib import Path
import json
import sys

output, raw_root = map(Path, sys.argv[1:3])
mode = sys.argv[3]
expected_entries = 4993 if mode == "legacy_failed_apply" else -1
output.write_text(f'''#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <iostream>
{{
  TFile input({json.dumps(str(raw_root))}, "READ");
  if (input.IsZombie() || input.TestBit(TFile::kRecovered)) gSystem->Exit(91);
  const char *required[] = {{"slimtree", "sim_cross_counting", "tracking_radiograph"}};
  for (const char *name : required)
    if (!input.GetListOfKeys()->FindObject(name)) gSystem->Exit(92);
  TTree *tree = dynamic_cast<TTree *>(input.Get("slimtree"));
  if (!tree || tree->GetEntries() <= 0) gSystem->Exit(93);
  if ({expected_entries} >= 0 && tree->GetEntries() != {expected_entries}) gSystem->Exit(94);
  std::cout << "PPG12_RAW_REUSE_AUDIT mode={mode} entries="
            << tree->GetEntries() << std::endl;
  gSystem->Exit(0);
}}
''')
PY
  failure_stage="ppg12_raw_reuse_root_audit"
  failure_stage_log="${runtime_dir}/audit_reused_ppg12_raw.log"
  root -l -b -q "$reuse_root_audit" \
    >"$failure_stage_log" 2>&1 || {
    tail -n 80 "$failure_stage_log" >&2 || true
    die "raw PPG12 reuse ROOT structure audit failed"
  }
  cp -p "$reuse_ppg_raw_root" "$ppg_raw_root"
  [[ "$(sha256_file "$ppg_raw_root")" == "$(sha256_file "$reuse_ppg_raw_root")" ]] || \
    die "copied PPG12 raw reuse source differs from its token-bound input"
  reuse_evidence_receipt="${runtime_dir}/ppg12_raw_reuse_receipt.json"
  preserve_reused_ppg_evidence \
    "$reuse_ppg_raw_contract" "$reuse_ppg_raw_root" "$ppg_raw_root" \
    "$ppg_log" "$reuse_evidence_receipt" || \
    die "failed to preserve exact reused PPG12 reconstruction evidence"
else
  failure_stage="ppg12_raw_reconstruction"
  failure_stage_log="$ppg_log"
  ppg_runner="${runtime_dir}/run_ppg12.C"
  python3 - "$ppg_runner" "$ppg_macro" "$ppg_wrapper" "$g4_slice" \
    "$truthjet_slice" "$ppg_raw_root" "$calo_macro_dir" <<'PY'
from pathlib import Path
import json
import sys
runner, macro, wrapper, g4, truth, output, calo_macro_dir = sys.argv[1:]
call = f"Fun4All_ppg12_fixed_seed_oracle({json.dumps(g4)},{json.dumps(truth)},{json.dumps(output)})"
Path(runner).write_text(f'''{{
  Int_t error = 0;
  gROOT->SetMacroPath((std::string({json.dumps(calo_macro_dir + ':')}) + gROOT->GetMacroPath()).c_str());
  if (gROOT->LoadMacro({json.dumps(macro)}) < 0) throw std::runtime_error("PPG12 macro load failed");
  if (gROOT->LoadMacro({json.dumps(wrapper)}) < 0) throw std::runtime_error("PPG12 wrapper load failed");
  gROOT->ProcessLine({json.dumps(call)}, &error);
  if (error != TInterpreter::kNoError) throw std::runtime_error("PPG12 oracle call failed");
}}
''')
PY
  (
    cd "$ppg_dir"
    export RJ_PPG12_ORACLE_RUNTIME_PROFILE=new.17
    export LD_LIBRARY_PATH="${ppg_lib_view}:${base_ld_library_path}"
    export ROOT_INCLUDE_PATH="${calo_macro_dir}:${base_root_include_path}"
    root -l -b -q "$ppg_runner"
  ) 2>&1 | tee "$ppg_log"
fi
[[ -s "$ppg_raw_root" ]] || \
  die "PPG12 executable did not write its exact cwd-local caloana.root"

failure_stage="ppg12_apply_bdt"
failure_stage_log="$apply_log"
apply_runner="${runtime_dir}/apply_ppg12_bdt.C"
python3 - "$apply_runner" "$apply_bdt" "$apply_runtime_config" "$ppg_raw_root" \
  "$recoeff_yaml_cpp" <<'PY'
from pathlib import Path
import json
import sys
runner, macro, config, root_path, yaml_cpp = sys.argv[1:]
call = f"apply_BDT({json.dumps(config)},\"data\",{json.dumps(root_path)})"
Path(runner).write_text(f'''{{
  Int_t error = 0;
  if (gSystem->Load({json.dumps(yaml_cpp)}) < 0) throw std::runtime_error("sealed yaml-cpp load failed");
  if (gROOT->LoadMacro({json.dumps(macro)}) < 0) throw std::runtime_error("apply_BDT load failed");
  gROOT->ProcessLine({json.dumps(call)}, &error);
  if (error != TInterpreter::kNoError) throw std::runtime_error("apply_BDT call failed");
}}
''')
PY
(
  cd "$(dirname "$apply_bdt")"
  export LD_LIBRARY_PATH="$(dirname "$recoeff_yaml_cpp"):${ppg_lib_view}:${base_ld_library_path}"
  export ROOT_INCLUDE_PATH="${recoeff_yaml_cpp_include_root}:${base_root_include_path}"
  root -l -b -q "$apply_runner"
) 2>&1 | tee "$apply_log"
[[ -s "$ppg_scored_root" ]] || die "apply_BDT did not write expected scored ROOT"
python3 - "$apply_evidence" "$apply_bdt" "$apply_config" \
  "$apply_runtime_config" "$ppg_raw_root" "$ppg_scored_root" \
  "$recoil_runtime_manifest" <<'PY'
from pathlib import Path
import hashlib
import json
import sys

output, macro, canonical_config, runtime_config, raw_root, scored_root, manifest = map(Path, sys.argv[1:])

def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

runtime = json.loads(manifest.read_text())
assets = {
    item["role"]: {"path": item["path"], "sha256": item["sha256"]}
    for item in runtime["files"]
    if str(item.get("role", "")).startswith("ppg_apply_model_")
    or item.get("role") == "ppg_apply_npb_model"
}
data = {
    "schema_version": 1,
    "stage": "source_locked_canonical_apply_bdt_split",
    "revision": "29f8223bd9b36dffab07961b597afa94185bbdf1",
    "macro": {"path": str(macro), "sha256": digest(macro)},
    "canonical_config": {"path": str(canonical_config), "sha256": digest(canonical_config)},
    "runtime_config": {"path": str(runtime_config), "sha256": digest(runtime_config)},
    "runtime_config_change": "omit unused CLUSTERINFO_CEMC_NO_SPLIT lane only",
    "model_assets": assets,
    "input": {"path": str(raw_root), "sha256": digest(raw_root)},
    "output": {"path": str(scored_root), "sha256": digest(scored_root)},
    "branches_required": [
        "cluster_bdt_CLUSTERINFO_CEMC_base_E",
        "cluster_bdt_CLUSTERINFO_CEMC_base_v3E",
        "cluster_npb_score_CLUSTERINFO_CEMC",
    ],
    "shared_by_recoeff_baseline_and_trace": True,
}
output.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
PY

run_recoeff() {
  local label="$1"
  local layout="$2"
  local macro="$3"
  local config="$4"
  local scan_log="$5"
  local analysis_log="$6"
  local vtxscan_root="$7"
  local trace_mode="$8"
  local scan_runner="${runtime_dir}/run_recoeff_${label}_scan.C"
  local analysis_runner="${runtime_dir}/run_recoeff_${label}.C"

  python3 - "$scan_runner" "$analysis_runner" "$macro" "$config" \
    "$vtxscan_root" "$oracle_sample" "$recoeff_roounfold" <<'PY'
from pathlib import Path
import json
import sys

scan_runner, analysis_runner, macro, config, vtxscan, sample, roounfold = sys.argv[1:]

def runner(call: str) -> str:
    return f'''#include <TUnfold.h>
{{
  Int_t error = 0;
  if (gSystem->Load({json.dumps(roounfold)}) < 0) throw std::runtime_error("sealed RooUnfold load failed");
  if (gROOT->LoadMacro({json.dumps(macro)}) < 0) throw std::runtime_error("RecoEff macro load failed");
  gROOT->ProcessLine({json.dumps(call)}, &error);
  if (error != TInterpreter::kNoError) throw std::runtime_error("RecoEff call failed");
}}
'''

scan_call = f'RecoEffCalculator_TTreeReader({json.dumps(config)},{json.dumps(sample)},true,1.0,"")'
analysis_call = f'RecoEffCalculator_TTreeReader({json.dumps(config)},{json.dumps(sample)},false,1.0,{json.dumps(vtxscan)})'
Path(scan_runner).write_text(runner(scan_call))
Path(analysis_runner).write_text(runner(analysis_call))
PY

  failure_stage="ppg12_recoeff_${label}_vertex_scan"
  failure_stage_log="$scan_log"
  (
    cd "$layout"
    export LD_LIBRARY_PATH="$(dirname "$recoeff_yaml_cpp"):$(dirname "$recoeff_roounfold"):${base_ld_library_path}"
    export ROOT_INCLUDE_PATH="${recoeff_roounfold_include_root}:${recoeff_include_root}:${recoeff_yaml_cpp_include_root}:${base_root_include_path}"
    unset RJ_PPG12_EXEC_TRACE_CSV RJ_PPG12_EXEC_RESPONSE_TRACE_CSV
    root -l -b -q "$scan_runner"
  ) 2>&1 | tee "$scan_log"
  [[ -s "$vtxscan_root" ]] || die "$label RecoEff vertex scan is missing"

  failure_stage="ppg12_recoeff_${label}_analysis"
  failure_stage_log="$analysis_log"
  (
    cd "$layout"
    export LD_LIBRARY_PATH="$(dirname "$recoeff_yaml_cpp"):$(dirname "$recoeff_roounfold"):${base_ld_library_path}"
    export ROOT_INCLUDE_PATH="${recoeff_roounfold_include_root}:${recoeff_include_root}:${recoeff_yaml_cpp_include_root}:${base_root_include_path}"
    if [[ "$trace_mode" == 1 ]]; then
      export RJ_PPG12_EXEC_TRACE_CSV="$ppg_candidate_trace"
      export RJ_PPG12_EXEC_RESPONSE_TRACE_CSV="$ppg_response_trace"
    else
      unset RJ_PPG12_EXEC_TRACE_CSV RJ_PPG12_EXEC_RESPONSE_TRACE_CSV
    fi
    root -l -b -q "$analysis_runner"
  ) 2>&1 | tee "$analysis_log"
}

# Both estimator executions consume the exact same scored slimtree and
# byte-identical config.  The first is scientifically unmodified; the second
# adds only the candidate/response side channels.  Exact ROOT payload equality
# is checked by the aggregate extractor before the candidate trace is trusted.
failure_stage="ppg12_recoeff_baseline"
failure_stage_log="$baseline_recoeff_log"
run_recoeff baseline "$baseline_layout" "$recoeff_macro" \
  "$baseline_recoeff_config" "$baseline_recoeff_scan_log" \
  "$baseline_recoeff_log" "$baseline_vtxscan_root" 0
failure_stage="ppg12_recoeff_trace"
failure_stage_log="$trace_recoeff_log"
run_recoeff trace "$trace_layout" "$recoeff_trace_macro" \
  "$trace_recoeff_config" "$trace_recoeff_scan_log" \
  "$trace_recoeff_log" "$trace_vtxscan_root" 1

for required in "$baseline_eff_root" "$trace_eff_root" \
  "$baseline_response_root" "$trace_response_root" \
  "$ppg_candidate_trace" "$ppg_response_trace"; do
  [[ -s "$required" ]] || die "RecoEff executable evidence is missing: $required"
done

failure_stage="ppg12_recoeff_root_equivalence"
failure_stage_log=""
python3 "$aggregate_extractor" \
  --baseline-root "$baseline_eff_root" \
  --instrumented-root "$trace_eff_root" \
  --baseline-response-root "$baseline_response_root" \
  --instrumented-response-root "$trace_response_root" \
  --runtime-manifest "$paired_runtime_manifest" \
  --root-equivalence-only \
  --out-json "$root_equivalence_report"

failure_stage="recoiljets_reconstruction"
failure_stage_log="$recoil_log"
recoil_runner="${runtime_dir}/run_recoiljets.C"
python3 - "$recoil_runner" "$recoil_macro" "$recoil_wrapper" \
  "$combined_list" "$recoil_root" "$calo_macro_dir" <<'PY'
from pathlib import Path
import json
import sys
runner, macro, wrapper, combined, output, calo_macro_dir = sys.argv[1:]
call = f"Fun4All_recoiljets_fixed_seed_oracle({json.dumps(combined)},{json.dumps(output)})"
Path(runner).write_text(f'''{{
  Int_t error = 0;
  gROOT->SetMacroPath((std::string({json.dumps(calo_macro_dir + ':')}) + gROOT->GetMacroPath()).c_str());
  if (gROOT->LoadMacro({json.dumps(macro)}) < 0) throw std::runtime_error("RecoilJets macro load failed");
  if (gROOT->LoadMacro({json.dumps(wrapper)}) < 0) throw std::runtime_error("RecoilJets wrapper load failed");
  gROOT->ProcessLine({json.dumps(call)}, &error);
  if (error != TInterpreter::kNoError) throw std::runtime_error("RecoilJets oracle call failed");
}}
''')
PY
(
  cd "$recoil_dir"
  export RJ_PPG12_ORACLE_RUNTIME_PROFILE=new.17
  export LD_LIBRARY_PATH="${recoil_lib_view}:${base_ld_library_path}"
  # The implementation includes <caloreco/PhotonClusterBuilder.h>; therefore
  # ROOT_INCLUDE_PATH must point at the isolated prefix's include directory,
  # not the nested include/caloreco directory itself.
  export ROOT_INCLUDE_PATH="${calo_macro_dir}:${photon_builder_include_root}:${base_root_include_path}"
  export RJ_CONFIG_YAML="$recoil_config"
  export RJ_DATASET=isSim
  export RJ_IS_SIM=1
  export RJ_SIM_SAMPLE="$recoil_sample"
  export RJ_PPG12_CLOSURE_CANARY=1
  export RJ_PPG12_CLOSURE_CANARY_ID="${lane_id}:pairedoracle"
  export RJ_PPG12_PPSIM_REPLAY_SEEDS="$ph_seed_sequence"
  export RJ_PPG12_PPSIM_EXPECT_PEDESTAL_SEQUENCE="$pedestal_sequence"
  unset RJ_PPG12_PPSIM_FIXED_RANDOMSEED
  unset RJ_PPG12_PPSIM_FIXED_PEDESTAL_SEQUENCE
  export RJ_PPG12_PHOTON_YIELD=1
  if [[ "$interaction" == DI ]]; then
    export RJ_SIM_SIGNAL_SAMPLE_SET=ppg12_double
    export RJ_PPG12_PHOTON_YIELD_DOUBLE=1
    export RJ_PPG12_PERIOD_STRICT_DI=1
  else
    unset RJ_SIM_SIGNAL_SAMPLE_SET RJ_PPG12_PERIOD_STRICT_DI
    export RJ_PPG12_PHOTON_YIELD_DOUBLE=0
  fi
  # No blanket 4% classification smear. RecoilJets records its explicit
  # response-only additive-resolution ET in a separate trace branch.
  unset RJ_PPG12_PHOTON_YIELD_CLUSTER_ERES
  export RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4=1
  export RJ_PPG12_PPSIM_G4_ONLY=1
  export RJ_SIM_ALLOW_NONE_LISTS=1
  export RJ_PPG12_PERIOD="$period"
  export RJ_PPG12_CROSSING_PERIOD="$period"
  export RJ_PPG12_PERIOD_USE_LUMI_WEIGHT=1
  export RJ_PPG12_PERIOD_ALLOW_ALL_SIM=0
  export RJ_PPG12_PERIOD_ALLOW_MIX_OVERRIDE=0
  export RJ_PPG12_PERIOD_ALLOW_VERTEX_FILE_OVERRIDE=0
  export RJ_PPG12_TABLE_QA=1
  export RJ_PPG12_TABLE_QA_NPB_DATA_TAGGING=1
  export RJ_PPG12_FIG13_PARITY_QA=1
  export RJ_PPG12_FIG11_SB_DIAGNOSTIC=1
  export RJ_PP_PHOTONID_TRAINING_TREE=1
  export RJ_PP_PHOTONID_TRAINING_TREE_MAX_ENTRIES=0
  export RJ_PP_PHOTONID_SOURCE_ROLE=signal
  # Trace the entire accepted candidate population.  Signal status is an
  # audited field, never a pre-comparison filter.
  export RJ_PP_PHOTONID_PPG12_FILTER=0
  export RJ_PP_PHOTONID_REQUIRE_PRESELECTION=0
  export RJ_PHOTON_ID_ROW_MATCH='newPPG12|newPPG12|newPPG12'
  export RJ_DISABLE_ID_FANOUT=1
  export RJ_ID_FANOUT_MAX_ROWS=1
  export RJ_DISABLE_JET_PT_INTERNALIZATION=1
  export RJ_DISABLE_DPHI_INTERNALIZATION=1
  export RJ_DIRECT_DST_DOALL=1
  export RJ_DIRECT_NEVENTS=0
  export RJ_AUTO_MERGE=0
  export RJ_INTERNAL_JET_PT_MINS='5.0,7.0,10.0,12.0'
  export RJ_INTERNAL_DPHI_PI_FRACTIONS='0.5,0.875'
  export RJ_DISABLE_JES_CDB_AUDIT=1
  export RJ_FAIL_ON_MISSING_CALO_INPUT=0
  root -l -b -q "$recoil_runner"
) 2>&1 | tee "$recoil_log"
[[ -s "$recoil_root" ]] || die "RecoilJets executable did not write output ROOT"

failure_stage="candidate_comparison"
failure_stage_log=""
python3 "$comparator" \
  --rj-root "$recoil_root" \
  --ppg12-root "$ppg_scored_root" \
  --mask-root "$tower_mask" \
  --ppg12-max-events 5000 \
  --events-per-segment 1000 \
  --base-e-model "$base_e_model" \
  --base-v3e-model "$base_v3e_model" \
  --ppg12-executable-trace "$ppg_candidate_trace" \
  --ppg12-executable-response-trace "$ppg_response_trace" \
  --lane-id "$lane_id" \
  --runtime-contract "$contract" \
  --out-md "$report_md" \
  --out-csv "$summary_csv" \
  --out-candidates-csv "$candidate_csv"

failure_stage="executable_aggregate_extraction"
python3 "$aggregate_extractor" \
  --baseline-root "$baseline_eff_root" \
  --instrumented-root "$trace_eff_root" \
  --baseline-response-root "$baseline_response_root" \
  --instrumented-response-root "$trace_response_root" \
  --trace-csv "$ppg_candidate_trace" \
  --response-trace-csv "$ppg_response_trace" \
  --candidate-csv "$candidate_csv" \
  --lane-id "$lane_id" \
  --runtime-contract "$contract" \
  --runtime-manifest "$paired_runtime_manifest" \
  --asset apply_bdt_stage_evidence="$apply_evidence" \
  --asset estimator_source="$recoeff_source_macro" \
  --asset estimator_roounfold_compatibility="$recoeff_compat_receipt" \
  --asset estimator_trace_transform="$recoeff_trace_receipt" \
  --out-json "$aggregate_report"

failure_stage="postrun_audit"
python3 "$auditor" postrun --contract "$contract"
printf 'PASS\n' > "$state_file"
completed=1
failure_stage="complete"
printf 'PPG12_PAIRED_ORACLE_PASS output=%s contract=%s report=%s\n' \
  "$output_dir" "$contract" "$report_md"
