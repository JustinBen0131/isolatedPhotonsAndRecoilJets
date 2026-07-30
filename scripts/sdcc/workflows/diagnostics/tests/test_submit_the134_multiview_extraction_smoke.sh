#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd -P)"
controller="${repo_root}/scripts/sdcc/workflows/diagnostics/submit_the134_multiview_extraction_smoke.sh"
macro="${repo_root}/macros/Fun4All_recoilJets_unified_impl.C"
pp_config="${repo_root}/macros/analysis_config_the119_pp_replay_foundation.yaml"

bash -n "$controller"
grep -Fq 'RJ_THE134_PP_BASE_E_MODEL is required for the PPG12 p+p-SIM route' "$controller"
grep -Fq 'p+p config/base_E model path mismatch' "$controller"
grep -Fq 'RJ_THE134_CAPACITY_SELECTED_ROW' "$controller"
grep -Fq 'RJ_THE134_CAPACITY_COMBINE_TAG' "$controller"
grep -Fq 'capacity repair matrix must contain exactly its selected frozen witness' "$controller"
grep -Fq 'full capacity mode must not name a preserved partner campaign' "$controller"
grep -Fq 'RJ_THE134_MULTIVIEW_TRAINING_FILE=${sidecar};RJ_THE134_MULTIVIEW_SIDECAR_ONLY_V1=1;RJ_THE134_EPHEMERAL_ANALYSIS_OUTPUT=1;RJ_THE134_EXPECTED_SOURCE_ROLE=${role}' "$controller"
grep -Fq 'sidecar="${row_output}/training_views/${sample}/RJPhotonTrainingViewV1.root"' "$controller"
grep -Fq 'fast_extraction_mode="${RJ_THE134_FAST_EXTRACTION_V1:-0}"' "$controller"
grep -Fq 'RJ_THE134_FAST_EXTRACTION_V1 must be exactly 0 or 1' "$controller"
grep -Fq 'extra="${extra};RJ_THE134_FAST_EXTRACTION_V1=1"' "$controller"
grep -Fq 'def parse_ephemeral_analysis_receipt(' "$controller"
grep -Fq 'RECOILJETS_THE134_EPHEMERAL_ANALYSIS_V1 ' "$controller"
grep -Fq 'THE134_MULTIVIEW_SIDECAR_ONLY_ARTIFACT_PROFILE_V2' "$controller"
grep -Fq '"analysis_artifact_state": "EPHEMERAL_VALIDATED_NOT_RETAINED"' "$controller"
grep -Fq '"analysis_artifact_state": "DURABLE_VALIDATED"' "$controller"
grep -Fq '"full_extraction_plan": {' "$controller"
grep -Fq 'analysis_path.exists()' "$controller"
grep -Fq 'root_health_identity_join_certificate_capacity_v6.json' "$controller"

invalid_fast_log="${TMPDIR:-/tmp}/the134_invalid_fast_extraction.$$"
if RJ_THE134_FAST_EXTRACTION_V1=2 "$controller" inventory >"$invalid_fast_log" 2>&1; then
  printf 'invalid fast-extraction mode was accepted\n' >&2
  exit 1
fi
grep -Fq 'RJ_THE134_FAST_EXTRACTION_V1 must be exactly 0 or 1' "$invalid_fast_log"
rm -f "$invalid_fast_log"
python3 - "$controller" <<'PY'
from pathlib import Path
import ast
import hashlib
import json
import re
import sys
import tempfile

source = Path(sys.argv[1]).read_text(encoding="utf-8")
root_validator = source.split(
    '"$capacity_full_plan" "$capacity_full_plan_sha" "$fast_extraction_mode" <<\'PY\'\n',
    1,
)[1].split("\nPY\n}", 1)[0]
compile(root_validator, "<validate_root_health_and_joins>", "exec")
for required in (
    "ephemeral analysis health receipt cardinality differs",
    "ephemeral worker receipt does not bind the retained sidecar",
    "fast extraction requires the exact frozen sidecar-only ephemeral artifact profile",
    "analysis_writer_root_ephemeral_receipt_v1",
    "analysis_writer_root_v1",
):
    if required not in root_validator:
        raise ValueError(f"missing ephemeral/durable validation contract: {required}")
tree = ast.parse(root_validator)
selected = [
    node
    for node in tree.body
    if isinstance(node, ast.FunctionDef)
    and node.name in {"file_sha256", "parse_ephemeral_analysis_receipt"}
]
namespace = {
    "Path": Path,
    "hashlib": hashlib,
    "json": json,
    "HEX64": re.compile(r"^[0-9a-f]{64}$"),
    "ANALYSIS_MINIMUM_BYTES": 50_000,
    "expected_fast_extraction": False,
}
exec(compile(ast.Module(body=selected, type_ignores=[]), "<receipt-parser>", "exec"), namespace)
payload = {
    "schema": "THE134_EPHEMERAL_ANALYSIS_HEALTH_V1",
    "status": "PASS",
    "mode": "EPHEMERAL_CONDOR_SCRATCH",
    "analysis_size_bytes": 50_000,
    "analysis_minimum_bytes": 50_000,
    "analysis_key_inventory_sha256": "a" * 64,
    "analysis_config_present": True,
    "analysis_directory_present": True,
    "analysis_histogram_present": True,
    "analysis_root_non_zombie": True,
    "analysis_root_non_recovered": True,
    "analysis_root_retained": False,
    "sidecar_size_bytes": 12_345,
    "sidecar_tree_entries": 7,
}
prefix = "RECOILJETS_THE134_EPHEMERAL_ANALYSIS_V1 "
with tempfile.TemporaryDirectory() as temporary:
    stdout = Path(temporary) / "worker.out"
    stdout.write_text(prefix + json.dumps(payload, sort_keys=True) + "\n")
    parsed_path, parsed_sha, parsed_payload = namespace[
        "parse_ephemeral_analysis_receipt"
    ]({"condor_stdout": str(stdout)})
    assert parsed_path == stdout.resolve()
    assert parsed_sha == hashlib.sha256(stdout.read_bytes()).hexdigest()
    assert parsed_payload == payload
    for mutation in ("duplicate", "schema", "missing_key"):
        changed = dict(payload)
        if mutation == "schema":
            changed["schema"] = "DRIFT"
        elif mutation == "missing_key":
            changed.pop("analysis_config_present")
        text = prefix + json.dumps(changed, sort_keys=True) + "\n"
        if mutation == "duplicate":
            text += text
        stdout.write_text(text)
        try:
            namespace["parse_ephemeral_analysis_receipt"](
                {"condor_stdout": str(stdout)}
            )
        except ValueError:
            pass
        else:
            raise AssertionError(f"ephemeral receipt mutation accepted: {mutation}")
required_tail = [
    "scripts/sdcc/workflows/diagnostics/submit_the134_multiview_extraction_smoke.sh",
    "scripts/sdcc/workflows/diagnostics/materialize_the134_full_multiview_extraction.py",
    "scripts/sdcc/workflows/diagnostics/project_the134_preextraction_storage_quota.py",
    "scripts/sdcc/workflows/diagnostics/resolve_the134_full_multiview_extraction.py",
    "scripts/sdcc/workflows/diagnostics/the134_full_extraction_controller.py",
]


def validate(text: str) -> None:
    body = text.split("compute_code_sha() {", 1)[1].split("\nEOF\n}", 1)[0]
    logical = [
        line.split("|", 1)[0]
        for line in body.splitlines()
        if "|" in line and not line.lstrip().startswith(("while ", "done "))
    ]
    if logical[-len(required_tail) :] != required_tail:
        raise ValueError("aggregate-code inputs do not match the frozen V10/V11 manifest order")


validate(source)
mutated = source.replace(required_tail[1] + "|", "removed-authority-input|", 1)
try:
    validate(mutated)
except ValueError:
    pass
else:
    raise SystemExit("aggregate-code membership mutation was not rejected")
PY
grep -Fq 'std::string ppg12_base_e_model_file = "";' "$macro"
grep -Fq 'cfg.ppg12_base_e_model_file = detail::trim(AfterColon(line));' "$macro"
grep -Fq 'std::string baseEModelFile = cfg.ppg12_base_e_model_file;' "$macro"
grep -Fq 'ppg12_base_e_model_file: /sphenix/user/shuhangli/ppg12/FunWithxgboost/binned_models/model_base_E_split_single_tmva.root' "$pp_config"
if RJ_THE134_CAPACITY_SELECTED_ROW=not_a_frozen_row \
  "$controller" capacity-preflight >"${TMPDIR:-/tmp}/the134_invalid_capacity_selector.$$" 2>&1; then
  printf 'invalid capacity selector was accepted\n' >&2
  exit 1
fi
grep -Fq 'capacity selector is not a frozen witness row' \
  "${TMPDIR:-/tmp}/the134_invalid_capacity_selector.$$"
rm -f "${TMPDIR:-/tmp}/the134_invalid_capacity_selector.$$"

if RJ_THE134_CAPACITY_SELECTED_ROW=pp_background_jet8 \
  "$controller" capacity-preflight >"${TMPDIR:-/tmp}/the134_missing_combine_tag.$$" 2>&1; then
  printf 'one-row capacity repair without a combine tag was accepted\n' >&2
  exit 1
fi
grep -Fq 'one-row capacity repair requires a safe RJ_THE134_CAPACITY_COMBINE_TAG' \
  "${TMPDIR:-/tmp}/the134_missing_combine_tag.$$"
rm -f "${TMPDIR:-/tmp}/the134_missing_combine_tag.$$"

if RJ_THE134_CAPACITY_COMBINE_TAG=unexpected_partner \
  "$controller" capacity-preflight >"${TMPDIR:-/tmp}/the134_unexpected_combine_tag.$$" 2>&1; then
  printf 'full capacity mode with a preserved-partner tag was accepted\n' >&2
  exit 1
fi
grep -Fq 'full capacity mode must not name a preserved partner campaign' \
  "${TMPDIR:-/tmp}/the134_unexpected_combine_tag.$$"
rm -f "${TMPDIR:-/tmp}/the134_unexpected_combine_tag.$$"
if grep -Eq '(^|[^[:alnum:]_])(/usr/bin/)?(say|afplay)([^[:alnum:]_]|$)|osascript|NSSound' "$controller"; then
  printf 'forbidden local audio command in controller: %s\n' "$controller" >&2
  exit 1
fi

tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/the134-tuple-contract.XXXXXX")"
trap 'chmod -R u+w "$tmpdir" 2>/dev/null || true; rm -rf "$tmpdir"' EXIT
die() { return 2; }
# Extracted controller helpers call the controller's text logger. Keep that
# logger local to the test and use a name that cannot resolve to a macOS
# text-to-speech command.
log_line() { :; }
sha_file() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}
require_sha() {
  local _name="$1" value="$2"
  [[ "$value" =~ ^[0-9a-f]{64}$ ]] || return 2
}
require_file_hash() {
  local _label="$1" path="$2" expected="$3"
  require_sha "$_label" "$expected"
  [[ -s "$path" && "$(sha_file "$path")" == "$expected" ]]
}
eval "$(
  sed -n '/^validate_five_file_tuple_count()/,/^validate_five_field_fanout_contract()/p' \
    "$controller" | sed '$d'
)"
eval "$(sed -n '/^validate_pp_sim_weight_contract()/,/^}/p' "$controller")"
eval "$(
  sed -n '/^validate_capacity_preflight_contract()/,/^yaml_value()/p' \
    "$controller" | sed '$d'
)"
eval "$(sed -n '/^yaml_value()/,/^}/p' "$controller")"
eval "$(sed -n '/^validate_five_field_fanout_contract()/,/^}/p' "$controller")"
eval "$(sed -n '/^validate_snapshot_macro_provider()/,/^}/p' "$controller")"
eval "$(
  sed -n '/^condor_field()/,/^analysis_tag_for_dataset()/p' \
    "$controller" | sed '$d'
)"
eval "$(
  sed -n '/^verify_materialization_attempt_seal()/,/^submit_row()/p' \
    "$controller" | sed '$d'
)"
eval "$(
  sed -n '/^write_runtime_authority_manifest()/,/^require_inputs_and_hashes()/p' \
    "$controller" | sed '$d'
)"
eval "$(
  sed -n '/^validate_system_matrix_or_capacity_audit()/,/^write_capacity_resource_certificate()/p' \
    "$controller" | sed '$d'
)"
(
  eval "$(
    sed -n '/^common_extra_env()/,/^system_extra_env()/p' \
      "$controller" | sed '$d'
  )"
  manifest_field() {
    case "$2" in
      10) printf '%064d\n' 1 ;;
      12) printf '/tmp/config.yaml\n' ;;
      13) printf '%064d\n' 2 ;;
      17) printf '%064d\n' 3 ;;
      28) printf '/sphenix/tg/example/training_views/Jet8/RJPhotonTrainingViewV1.root\n' ;;
      *) return 2 ;;
    esac
  }
  pp_period=0mrad
  capture_et_min=5.0
  capacity_mode=0
  RJ_THE134_CODE_SHA256="$(printf '%064d' 4)"
  RJ_THE134_REPLAY_SCHEMA_SHA256="$(printf '%064d' 5)"
  RJ_THE134_SEMANTIC_SHA256="$(printf '%064d' 6)"
  fast_extraction_mode=0
  ordinary_extra="$(common_extra_env pp_background_jet8 pp pp_inclusive_sim Jet8 Jet8 BACKGROUND)"
  if [[ "$ordinary_extra" == *RJ_THE134_FAST_EXTRACTION_V1=* ]]; then
    printf 'fast-extraction flag leaked into the default capacity environment\n' >&2
    exit 1
  fi
  fast_extraction_mode=1
  fast_extra="$(common_extra_env pp_background_jet8 pp pp_inclusive_sim Jet8 Jet8 BACKGROUND)"
  [[ "$fast_extra" == *';RJ_THE134_FAST_EXTRACTION_V1=1'* ]] || {
    printf 'fast-extraction flag was not propagated into the candidate environment\n' >&2
    exit 1
  }
)

validate_pp_sim_weight_contract 0mrad
validate_pp_sim_weight_contract 1p5mrad
for invalid_pp_period in "" run28 1.5mrad; do
  if ( validate_pp_sim_weight_contract "$invalid_pp_period" ) >/dev/null 2>&1; then
    printf 'invalid p+p period authority was accepted: %s\n' "${invalid_pp_period:-<empty>}" >&2
    exit 1
  fi
done

capacity_mode=1
capacity_canary_id=the134_capacity_fixture
capacity_selected_row=
capacity_combine_tag=
capacity_pp_row=pp_background_jet8
capacity_auau_row=auau_background_jet12
output_root=/sphenix/tg/example/replay_foundation/the134_capacity_fixture
evidence_root=/sphenix/u/example/evidence/qa/the134_capacity_fixture
submit_root=/sphenix/u/example/condor/the134_capacity_fixture
capacity_full_plan="${tmpdir}/capacity_plan.json"
capacity_preflight_receipt="${tmpdir}/capacity_receipt.json"
python3 - "$capacity_full_plan" <<'PY'
from pathlib import Path
import json
import sys

path = Path(sys.argv[1])
payload = {
    "campaign": {
        "output_root": "/sphenix/tg/example/training/the134_full_fixture",
        "submit_root": "/sphenix/u/example/condor/the134_full_fixture",
    },
    "execution_partition": {
        "capacity_authority_earned": False,
        "capacity_canary_required_before_submission": True,
        "group_size": 7,
        "schema": "THE134_FULL_EXTRACTION_PARTITION_CONTRACT_V2",
    },
    "execution_state": "PREFLIGHT_ONLY_NO_CONDOR_MUTATION",
    "full_training_authority": 0,
    "input_manifests": {
        "bundle": {"sha256": "b" * 64},
        "materialization": {
            "readback": "PASS_SYMLINK_FREE_READONLY_CONTENT_EXACT",
            "sha256": "c" * 64,
        },
    },
    "rows": [
        {
            "full_training_authority": 0,
            "input_contract": {"group_size": 7},
            "row_id": "pp_background_jet8",
            "sample": "run28_jet8",
            "system": "pp",
        },
        {
            "full_training_authority": 0,
            "input_contract": {"group_size": 7},
            "row_id": "auau_background_jet12",
            "sample": "run28_embeddedJet12",
            "system": "auau",
        },
    ],
    "schema": "THE134_FULL_MULTIVIEW_EXTRACTION_PLAN_V2",
    "status": "PREFLIGHT_PASS",
    "submission_performed": False,
}
path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
capacity_full_plan_sha="$(sha_file "$capacity_full_plan")"
python3 - "$capacity_preflight_receipt" "$capacity_full_plan_sha" <<'PY'
from pathlib import Path
import json
import sys

path = Path(sys.argv[1])
payload = {
    "artifacts": {"plan": {"sha256": sys.argv[2]}},
    "authority_state": "PREFLIGHT_RESOLVED_NOT_EARNED",
    "bundle_manifest_sha256": "b" * 64,
    "execution_partition_sha256": "a" * 64,
    "full_training_authority": 0,
    "materialization_receipt_sha256": "c" * 64,
    "schema": "THE134_FULL_MULTIVIEW_EXTRACTION_PREFLIGHT_RECEIPT_V2",
    "status": "PASS",
    "submission_performed": False,
}
path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
capacity_preflight_receipt_sha="$(sha_file "$capacity_preflight_receipt")"
validate_capacity_preflight_contract
python3 - "$capacity_full_plan" <<'PY'
from pathlib import Path
import json
import sys

path = Path(sys.argv[1])
payload = json.loads(path.read_text())
payload["artifact_profile"] = {
    "analysis_root_role": "EPHEMERAL_CONDOR_SCRATCH_VALIDATED_NOT_RETAINED",
    "artifact_profile": "THE134_MULTIVIEW_SIDECAR_ONLY_V2",
    "retained_analysis_root_count_per_job": 0,
    "schema": "THE134_MULTIVIEW_SIDECAR_ONLY_ARTIFACT_PROFILE_V2",
    "training_sidecar_role": "RJPhotonTrainingViewV1",
}
payload["execution_partition"].update(
    {
        "expected_durable_root_artifact_count": 18577,
        "expected_job_count": 18577,
        "expected_retained_analysis_output_count": 0,
        "expected_sidecar_output_count": 18577,
        "schema": "THE134_FULL_EXTRACTION_PARTITION_CONTRACT_V3",
    }
)
path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
capacity_full_plan_sha="$(sha_file "$capacity_full_plan")"
python3 - "$capacity_preflight_receipt" "$capacity_full_plan_sha" <<'PY'
from pathlib import Path
import json
import sys

path = Path(sys.argv[1])
payload = json.loads(path.read_text())
payload["artifacts"]["plan"]["sha256"] = sys.argv[2]
path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
capacity_preflight_receipt_sha="$(sha_file "$capacity_preflight_receipt")"
validate_capacity_preflight_contract
python3 - "$capacity_full_plan" <<'PY'
from pathlib import Path
import json
import sys

path = Path(sys.argv[1])
payload = json.loads(path.read_text())
payload["execution_partition"]["expected_retained_analysis_output_count"] = 1
path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
capacity_full_plan_sha="$(sha_file "$capacity_full_plan")"
python3 - "$capacity_preflight_receipt" "$capacity_full_plan_sha" <<'PY'
from pathlib import Path
import json
import sys

path = Path(sys.argv[1])
payload = json.loads(path.read_text())
payload["artifacts"]["plan"]["sha256"] = sys.argv[2]
path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
capacity_preflight_receipt_sha="$(sha_file "$capacity_preflight_receipt")"
if ( validate_capacity_preflight_contract ) >/dev/null 2>&1; then
  printf 'V3 partition with retained analysis outputs was accepted\n' >&2
  exit 1
fi
python3 - "$capacity_full_plan" <<'PY'
from pathlib import Path
import json
import sys

path = Path(sys.argv[1])
payload = json.loads(path.read_text())
payload["execution_partition"]["expected_retained_analysis_output_count"] = 0
payload["artifact_profile"]["artifact_profile"] = "THE134_MULTIVIEW_SIDECAR_ONLY_V1"
path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
capacity_full_plan_sha="$(sha_file "$capacity_full_plan")"
python3 - "$capacity_preflight_receipt" "$capacity_full_plan_sha" <<'PY'
from pathlib import Path
import json
import sys

path = Path(sys.argv[1])
payload = json.loads(path.read_text())
payload["artifacts"]["plan"]["sha256"] = sys.argv[2]
path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
capacity_preflight_receipt_sha="$(sha_file "$capacity_preflight_receipt")"
if ( validate_capacity_preflight_contract ) >/dev/null 2>&1; then
  printf 'V3 partition with a non-V2 sidecar artifact profile was accepted\n' >&2
  exit 1
fi
python3 - "$capacity_full_plan" <<'PY'
from pathlib import Path
import json
import sys

path = Path(sys.argv[1])
payload = json.loads(path.read_text())
payload.pop("artifact_profile")
payload["execution_partition"] = {
    "capacity_authority_earned": False,
    "capacity_canary_required_before_submission": True,
    "group_size": 7,
    "schema": "THE134_FULL_EXTRACTION_PARTITION_CONTRACT_V2",
}
path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
capacity_full_plan_sha="$(sha_file "$capacity_full_plan")"
python3 - "$capacity_preflight_receipt" "$capacity_full_plan_sha" <<'PY'
from pathlib import Path
import json
import sys

path = Path(sys.argv[1])
payload = json.loads(path.read_text())
payload["artifacts"]["plan"]["sha256"] = sys.argv[2]
path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
capacity_preflight_receipt_sha="$(sha_file "$capacity_preflight_receipt")"
validate_capacity_preflight_contract
python3 - "$capacity_preflight_receipt" <<'PY'
from pathlib import Path
import json
import sys

path = Path(sys.argv[1])
payload = json.loads(path.read_text())
payload["schema"] = "THE134_FULL_MULTIVIEW_EXTRACTION_PREFLIGHT_RECEIPT_V1"
path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
capacity_preflight_receipt_sha="$(sha_file "$capacity_preflight_receipt")"
if ( validate_capacity_preflight_contract ) >/dev/null 2>&1; then
  printf 'hybrid V1 receipt plus V2 plan/partition was accepted\n' >&2
  exit 1
fi
python3 - "$capacity_full_plan" <<'PY'
from pathlib import Path
import json
import sys

path = Path(sys.argv[1])
payload = json.loads(path.read_text())
payload["schema"] = "THE134_FULL_MULTIVIEW_EXTRACTION_PLAN_V1"
payload["execution_partition"]["schema"] = (
    "THE134_FULL_EXTRACTION_PARTITION_CONTRACT_V1"
)
path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
capacity_full_plan_sha="$(sha_file "$capacity_full_plan")"
python3 - "$capacity_preflight_receipt" "$capacity_full_plan_sha" <<'PY'
from pathlib import Path
import json
import sys

path = Path(sys.argv[1])
payload = json.loads(path.read_text())
payload["artifacts"]["plan"]["sha256"] = sys.argv[2]
path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
capacity_preflight_receipt_sha="$(sha_file "$capacity_preflight_receipt")"
validate_capacity_preflight_contract
python3 - "$capacity_full_plan" "$capacity_preflight_receipt" <<'PY'
from pathlib import Path
import json
import sys

plan_path = Path(sys.argv[1])
receipt_path = Path(sys.argv[2])
plan = json.loads(plan_path.read_text())
receipt = json.loads(receipt_path.read_text())
plan["schema"] = "THE134_FULL_MULTIVIEW_EXTRACTION_PLAN_V2"
plan["execution_partition"]["schema"] = (
    "THE134_FULL_EXTRACTION_PARTITION_CONTRACT_V1"
)
plan_path.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
receipt["schema"] = "THE134_FULL_MULTIVIEW_EXTRACTION_PREFLIGHT_RECEIPT_V2"
receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
PY
capacity_full_plan_sha="$(sha_file "$capacity_full_plan")"
python3 - "$capacity_preflight_receipt" "$capacity_full_plan_sha" <<'PY'
from pathlib import Path
import json
import sys

path = Path(sys.argv[1])
payload = json.loads(path.read_text())
payload["artifacts"]["plan"]["sha256"] = sys.argv[2]
path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
capacity_preflight_receipt_sha="$(sha_file "$capacity_preflight_receipt")"
if ( validate_capacity_preflight_contract ) >/dev/null 2>&1; then
  printf 'hybrid V2 plan plus V1 partition was accepted\n' >&2
  exit 1
fi
python3 - "$capacity_full_plan" <<'PY'
from pathlib import Path
import json
import sys
path = Path(sys.argv[1])
payload = json.loads(path.read_text())
payload["execution_partition"]["schema"] = (
    "THE134_FULL_EXTRACTION_PARTITION_CONTRACT_V2"
)
payload["execution_partition"]["group_size"] = 6
path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
capacity_full_plan_sha="$(sha_file "$capacity_full_plan")"
python3 - "$capacity_preflight_receipt" "$capacity_full_plan_sha" <<'PY'
from pathlib import Path
import json
import sys
path = Path(sys.argv[1])
payload = json.loads(path.read_text())
payload["artifacts"]["plan"]["sha256"] = sys.argv[2]
path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
capacity_preflight_receipt_sha="$(sha_file "$capacity_preflight_receipt")"
if ( validate_capacity_preflight_contract ) >/dev/null 2>&1; then
  printf 'semantic group-size drift was accepted by the capacity preflight contract\n' >&2
  exit 1
fi

for required_fail_closed_text in \
  'validate_capacity_preflight_contract ||' \
  'validate_calo_reco_build_authority ||' \
  'if ! hashes="$(require_inputs_and_hashes)"; then'
do
  if ! grep -F "$required_fail_closed_text" "$controller" >/dev/null; then
    printf 'missing fail-closed input-authority propagation: %s\n' \
      "$required_fail_closed_text" >&2
    exit 1
  fi
done

fake_validator="${tmpdir}/prepare_the134_h70_matrix.py"
fake_contract="${tmpdir}/the134_h70_contract.py"
cat > "$fake_validator" <<'PY'
from pathlib import Path

TREE_NAME = "RJPhotonTrainingViewV1"
CORE_BRANCHES = ("definition_name",)


def read_tree(path: Path, _tree_name: str):
    state, source_id = path.read_text().strip().split("|", 1)
    populated = state == "POPULATED"
    arrays = {
        "definition_name": ["H70"] * (7 if populated else 0),
        "_source_id": source_id,
    }
    metadata = {
        "tree_num_entries": len(arrays["definition_name"]),
        "schema_sha256": "1" * 64,
        "pp_feature_contract_sha256": "2" * 64,
        "auau_feature_contract_sha256": "3" * 64,
    }
    return arrays, ["definition_name"], metadata


def validate_artifact_metadata(_metadata, _arrays, *, external):
    return []


def validate_file_arrays(arrays, *, system, source, view_name):
    return [], [], {
        "status": "PASS",
        "rows": len(arrays["definition_name"]),
        "candidates": 1,
        "selected_view": view_name,
        "selected_view_rows": 1,
        "selected_view_training_rows": 1,
        "definition_counts": {"H70": 1},
        "source": source,
        "source_occurrence_ids": [arrays["_source_id"]],
        "failures": [],
    }
PY
printf 'fixture contract\n' > "$fake_contract"

make_capacity_audit_fixture() {
  local base="$1" state="$2" system="$3" row_id="$4" sample="$5"
  mkdir -p "$base"
  python3 - "$base" "$state" "$system" "$row_id" "$sample" <<'PY'
from pathlib import Path
import csv
import hashlib
import json
import sys

base = Path(sys.argv[1])
state, system, row_id, sample = sys.argv[2:]
lane = "background"
dataset = "run28" if system == "pp" else "run28auau"
period = "0mrad" if system == "pp" else "AUAU_RUN24"
role = "SI" if system == "pp" else "EMBEDDED"
chunk_path = base / "chunk.list"
chunk_path.write_text("fixture input tuple\n")
chunk_sha = hashlib.sha256(chunk_path.read_bytes()).hexdigest()
args_path = base / "job.args"
args_line = f"{sample} {chunk_path} {dataset} 123 0 1 NONE {base / 'output'}"
args_path.write_text(args_line + "\n")
submitted_args_sha = hashlib.sha256(f"0 {args_line}\n".encode()).hexdigest()
manifest_sha = "c" * 64
config_sha = "d" * 64
code_sha = "e" * 64
source_contract = {
    "lane": lane,
    "dataset": dataset,
    "sample": sample,
    "period": period,
    "run": 28,
    "segment": 1,
    "si_di_role": role,
    "ownership_state": "source_role_frozen",
    "input_uri_hash": chunk_sha,
    "input_file_sha256": chunk_sha,
    "source_manifest_sha256": manifest_sha,
}
canonical = "|".join(
    str(source_contract[field])
    for field in (
        "lane",
        "dataset",
        "sample",
        "period",
        "run",
        "segment",
        "input_uri_hash",
        "input_file_sha256",
        "source_manifest_sha256",
    )
)
identity_sha = hashlib.sha256(canonical.encode()).hexdigest()
source_id = identity_sha[:32] if identity_sha[:32] != "0" * 32 else "0" * 31 + "1"
sidecar = base / "sidecar.root"
sidecar.write_text(f"{state}|{source_id}")
sidecar_sha = hashlib.sha256(sidecar.read_bytes()).hexdigest()
entries = 7 if state == "POPULATED" else 0
candidates = 1 if state == "POPULATED" else 0
report = {
    "row_id": row_id,
    "status": "PASS",
    "population_state": state,
    "replay_sources": 1,
    "sidecar": str(sidecar),
    "sidecar_bytes": sidecar.stat().st_size,
    "sidecar_sha256": sidecar_sha,
    "sidecar_entries": entries,
    "sidecar_candidates": candidates,
    "joined_sidecar_rows": entries,
    "source_contract": source_contract,
    "source_execution_contract": {
        "args_file": str(args_path.resolve()),
        "chunk_index": 1,
        "run": 28,
        "staged_chunk_list": str(chunk_path.resolve()),
        "submitted_args_sha256": submitted_args_sha,
    },
    "source_occurrence_id_hex": source_id,
    "source_identity_canonical_sha256": identity_sha,
}
(base / "root.json").write_text(json.dumps({
    "schema": "THE134_SMOKE_ROOT_HEALTH_IDENTITY_JOIN_V1",
    "status": "PASS",
    "scope": "capacity",
    "capacity_authority_earned": True,
    "full_training_authority": 0,
    "execution_group_size": 7,
    "source_occurrences_per_output": 1,
    "rows": [report],
}, indent=2, sort_keys=True) + "\n")
with (base / "manifest.tsv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=[
        "row_id", "system", "full_training_authority", "lane", "dataset",
        "sample", "source_manifest_sha256", "resolved_config_sha256",
        "code_sha256",
    ], delimiter="\t")
    writer.writeheader()
    writer.writerow({
        "row_id": row_id,
        "system": system,
        "full_training_authority": "0",
        "lane": lane,
        "dataset": dataset,
        "sample": sample,
        "source_manifest_sha256": manifest_sha,
        "resolved_config_sha256": config_sha,
        "code_sha256": code_sha,
    })
with (base / "receipt.tsv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=[
        "row_id", "multiview_sidecar", "staged_chunk_sha256",
        "args_file", "staged_chunk_list", "submitted_args_sha256",
    ], delimiter="\t")
    writer.writeheader()
    writer.writerow({
        "row_id": row_id,
        "multiview_sidecar": str(sidecar),
        "staged_chunk_sha256": chunk_sha,
        "args_file": str(args_path),
        "staged_chunk_list": str(chunk_path),
        "submitted_args_sha256": submitted_args_sha,
    })
(base / "provenance.json").write_text(json.dumps({
    "schema": "THE134_SOURCE_PROVENANCE_V1",
    "inputs": [{
        "path": str(sidecar),
        "row_id": row_id,
        "system": system,
        "source_sample": sample,
        "input_uri_sha256": chunk_sha,
        "input_file_sha256": chunk_sha,
        "source_manifest_sha256": manifest_sha,
        "config_sha256": config_sha,
        "code_sha256": code_sha,
    }],
}, indent=2, sort_keys=True) + "\n")
(base / "sidecars.list").write_text(str(sidecar) + "\n")
PY
}

capacity_mode=1
canonical_validator="$fake_validator"
validator="$fake_validator"
canonical_validator_contract="$fake_contract"
capacity_file_validator_sha="$(sha_file "$fake_validator")"
capacity_contract_validator_sha="$(sha_file "$fake_contract")"
pp_period=0mrad

empty_fixture="${tmpdir}/capacity-empty"
make_capacity_audit_fixture \
  "$empty_fixture" VALID_EMPTY pp pp_background_jet8 run28_jet8
root_health_join_certificate="${empty_fixture}/root.json"
submission_receipt="${empty_fixture}/receipt.tsv"
submission_manifest="${empty_fixture}/manifest.tsv"
validate_system_matrix_or_capacity_audit \
  pp pp_background_jet8 "${empty_fixture}/sidecars.list" \
  "${empty_fixture}/provenance.json" "${empty_fixture}/matrix.npz" \
  "${empty_fixture}/audit.json"
python3 - "${empty_fixture}/audit.json" <<'PY'
import json
import sys
audit = json.load(open(sys.argv[1]))
assert audit["status"] == "PASS"
assert audit["scope"] == "capacity"
assert audit["population_state"] == "VALID_EMPTY"
assert audit["matrix_materialized"] is False
assert audit["file_report"]["source_occurrence_ids"] == []
assert audit["file_report"]["source_occurrence_binding"] == "PAIRED_REPLAY_SOURCE_ROW"
PY
chmod u+w "${empty_fixture}/audit.json"
python3 - "${empty_fixture}/audit.json" <<'PY'
from pathlib import Path
import json
import sys
path = Path(sys.argv[1])
payload = json.loads(path.read_text())
payload["status"] = "FAIL"
path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
if ( validate_system_matrix_or_capacity_audit \
  pp pp_background_jet8 "${empty_fixture}/sidecars.list" \
  "${empty_fixture}/provenance.json" "${empty_fixture}/matrix.npz" \
  "${empty_fixture}/audit.json" ) >/dev/null 2>&1; then
  printf 'differing existing capacity audit was overwritten\n' >&2
  exit 1
fi

populated_fixture="${tmpdir}/capacity-populated"
make_capacity_audit_fixture \
  "$populated_fixture" POPULATED auau auau_background_jet12 run28_embeddedJet12
root_health_join_certificate="${populated_fixture}/root.json"
submission_receipt="${populated_fixture}/receipt.tsv"
submission_manifest="${populated_fixture}/manifest.tsv"
validate_system_matrix_or_capacity_audit \
  auau auau_background_jet12 "${populated_fixture}/sidecars.list" \
  "${populated_fixture}/provenance.json" "${populated_fixture}/matrix.npz" \
  "${populated_fixture}/audit.json"
python3 - "${populated_fixture}/audit.json" <<'PY'
import json
import sys
audit = json.load(open(sys.argv[1]))
assert audit["status"] == "PASS"
assert audit["population_state"] == "POPULATED"
assert audit["file_report"]["rows"] == 7
assert audit["file_report"]["source_occurrence_ids"] == [
    audit["source_occurrence_id_hex"]
]
PY
cp "${populated_fixture}/sidecar.root" "${populated_fixture}/sidecar.original"
python3 - "${populated_fixture}/sidecar.root" <<'PY'
from pathlib import Path
import sys
path = Path(sys.argv[1])
text = path.read_text()
path.write_text(("X" if text[0] != "X" else "Y") + text[1:])
assert len(path.read_text()) == len(text)
PY
if ( validate_system_matrix_or_capacity_audit \
  auau auau_background_jet12 "${populated_fixture}/sidecars.list" \
  "${populated_fixture}/provenance.json" "${populated_fixture}/matrix.npz" \
  "${populated_fixture}/audit.json" ) >/dev/null 2>&1; then
  printf 'same-size post-audit sidecar mutation was accepted\n' >&2
  exit 1
fi
cp "${populated_fixture}/sidecar.original" "${populated_fixture}/sidecar.root"
validate_system_matrix_or_capacity_audit \
  auau auau_background_jet12 "${populated_fixture}/sidecars.list" \
  "${populated_fixture}/provenance.json" "${populated_fixture}/matrix.npz" \
  "${populated_fixture}/audit.json"

capacity_mode=0

printf '%s\n%s\n' \
  'getenv = False' \
  'environment = RJ_PPG12_PHOTON_YIELD=1;RJ_PPG12_PHOTON_YIELD_DOUBLE=0;RJ_PPG12_PERIOD=0mrad;RJ_PPG12_PERIOD_USE_LUMI_WEIGHT=1;RJ_PPG12_PERIOD_STRICT_DI=0;RJ_PPG12_PERIOD_ALLOW_ALL_SIM=0;RJ_PPG12_PERIOD_ALLOW_MIX_OVERRIDE=0;RJ_PPG12_PERIOD_ALLOW_VERTEX_FILE_OVERRIDE=0' \
  > "${tmpdir}/valid-pp-weight.sub"
[[ "$(condor_field "${tmpdir}/valid-pp-weight.sub" getenv)" == False ]] || {
  printf 'sealed descriptor getenv=False was not recovered exactly\n' >&2
  exit 1
}
require_descriptor_env_exact unit "${tmpdir}/valid-pp-weight.sub" RJ_PPG12_PHOTON_YIELD 1
require_descriptor_env_exact unit "${tmpdir}/valid-pp-weight.sub" RJ_PPG12_PHOTON_YIELD_DOUBLE 0
require_descriptor_env_exact unit "${tmpdir}/valid-pp-weight.sub" RJ_PPG12_PERIOD 0mrad
require_descriptor_env_absent unit "${tmpdir}/valid-pp-weight.sub" RJ_PPG12_PHOTON_YIELD_MIX_WEIGHT
require_descriptor_env_absent unit "${tmpdir}/valid-pp-weight.sub" RJ_PP_VERTEX_REWEIGHT_FILE

sed 's/RJ_PPG12_PERIOD=0mrad/RJ_PPG12_PERIOD=0mrad;RJ_PPG12_PERIOD=1p5mrad/' \
  "${tmpdir}/valid-pp-weight.sub" > "${tmpdir}/duplicate-pp-period.sub"
if ( require_descriptor_env_exact unit "${tmpdir}/duplicate-pp-period.sub" RJ_PPG12_PERIOD 0mrad ) \
  >/dev/null 2>&1; then
  printf 'duplicate p+p period authority was accepted in the descriptor\n' >&2
  exit 1
fi

sed 's/$/;RJ_PPG12_PHOTON_YIELD_MIX_WEIGHT=0.776/' \
  "${tmpdir}/valid-pp-weight.sub" > "${tmpdir}/manual-pp-mix.sub"
if ( require_descriptor_env_absent unit "${tmpdir}/manual-pp-mix.sub" RJ_PPG12_PHOTON_YIELD_MIX_WEIGHT ) \
  >/dev/null 2>&1; then
  printf 'manual p+p mix-weight override was accepted in the descriptor\n' >&2
  exit 1
fi

printf '/calo\t/g4\t/jets\t/global\t/mbd\n' > "${tmpdir}/valid-pp.list"
validate_one_five_file_tuple unit pp "${tmpdir}/valid-pp.list"

printf '/calo\t/g4\t/jets\t/global\tNONE\n' > "${tmpdir}/valid-auau.list"
validate_one_five_file_tuple unit auau "${tmpdir}/valid-auau.list"

if validate_one_five_file_tuple unit pp "${tmpdir}/valid-auau.list"; then
  printf 'Au+Au NONE tuple leaked into p+p authorization\n' >&2
  exit 1
fi

if validate_one_five_file_tuple unit auau "${tmpdir}/valid-pp.list"; then
  printf 'real Au+Au column-five path was not rejected\n' >&2
  exit 1
fi

printf '/calo\t/g4\tNONE\t/global\tNONE\n' > "${tmpdir}/bad-auau-none-column.list"
if validate_one_five_file_tuple unit auau "${tmpdir}/bad-auau-none-column.list"; then
  printf 'Au+Au NONE outside column five was not rejected\n' >&2
  exit 1
fi

printf 'calo\tg4\tjets\tglobal\n' > "${tmpdir}/four-columns.list"
if validate_one_five_file_tuple unit pp "${tmpdir}/four-columns.list"; then
  printf 'four-column tuple was not rejected\n' >&2
  exit 1
fi

printf 'calo\tg4\tjets\tglobal\tmbd\ncalo2\tg42\tjets2\tglobal2\tmbd2\n' > "${tmpdir}/two-rows.list"
if validate_one_five_file_tuple unit pp "${tmpdir}/two-rows.list"; then
  printf 'two-row tuple ownership was not rejected\n' >&2
  exit 1
fi

: > "${tmpdir}/valid-seven-pp.list"
for index in 1 2 3 4 5 6 7; do
  printf '/calo%s\t/g4%s\t/jets%s\t/global%s\t/mbd%s\n' \
    "$index" "$index" "$index" "$index" "$index" \
    >> "${tmpdir}/valid-seven-pp.list"
done
validate_five_file_tuple_count unit pp "${tmpdir}/valid-seven-pp.list" 7
if validate_five_file_tuple_count unit pp "${tmpdir}/valid-seven-pp.list" 1; then
  printf 'seven-row capacity chunk leaked into one-row smoke authority\n' >&2
  exit 1
fi
head -n 6 "${tmpdir}/valid-seven-pp.list" > "${tmpdir}/six-pp.list"
if validate_five_file_tuple_count unit pp "${tmpdir}/six-pp.list" 7; then
  printf 'six-row chunk was accepted as a group-of-seven capacity witness\n' >&2
  exit 1
fi

printf 'calo\tg4\t\tglobal\tmbd\n' > "${tmpdir}/empty-column.list"
if validate_one_five_file_tuple unit pp "${tmpdir}/empty-column.list"; then
  printf 'empty tuple column was not rejected\n' >&2
  exit 1
fi

printf '%s\n' 'coneR: 0.40' > "${tmpdir}/valid.yaml"
printf '%s\n' 'dest|cfg|pre|tight|nonTight' > "${tmpdir}/valid.fanout"
validate_five_field_fanout_contract unit "${tmpdir}/valid.fanout" "${tmpdir}/valid.yaml" 0.40

mkdir -p "${tmpdir}/physical-snapshot/lib"
printf 'provider\n' > "${tmpdir}/physical-snapshot/lib/libcalo_reco.so"
ln -s "${tmpdir}/physical-snapshot" "${tmpdir}/logical-snapshot"
printf 'R__LOAD_LIBRARY(%s)\n' \
  "${tmpdir}/logical-snapshot/lib/libcalo_reco.so" \
  > "${tmpdir}/valid-calo-macro.C"
validate_snapshot_macro_provider \
  unit "${tmpdir}/valid-calo-macro.C" \
  "${tmpdir}/physical-snapshot/lib/libcalo_reco.so"

printf 'provider\n' > "${tmpdir}/foreign-libcalo_reco.so"
if validate_snapshot_macro_provider \
  unit "${tmpdir}/valid-calo-macro.C" \
  "${tmpdir}/foreign-libcalo_reco.so" >/dev/null 2>&1; then
  printf 'same-byte foreign CaloReco provider was accepted as the snapshot provider\n' >&2
  exit 1
fi

printf 'R__LOAD_LIBRARY(%s)\nR__LOAD_LIBRARY(%s)\n' \
  "${tmpdir}/logical-snapshot/lib/libcalo_reco.so" \
  "${tmpdir}/physical-snapshot/lib/libcalo_reco.so" \
  > "${tmpdir}/duplicate-calo-macro.C"
if validate_snapshot_macro_provider \
  unit "${tmpdir}/duplicate-calo-macro.C" \
  "${tmpdir}/physical-snapshot/lib/libcalo_reco.so" >/dev/null 2>&1; then
  printf 'duplicate CaloReco macro providers were accepted\n' >&2
  exit 1
fi

printf 'R__LOAD_LIBRARY(libcalo_reco.so)\n' > "${tmpdir}/relative-calo-macro.C"
if validate_snapshot_macro_provider \
  unit "${tmpdir}/relative-calo-macro.C" \
  "${tmpdir}/physical-snapshot/lib/libcalo_reco.so" >/dev/null 2>&1; then
  printf 'relative CaloReco macro provider was accepted\n' >&2
  exit 1
fi

fanout_mutations=(
  'dest|cfg|pre|tight'
  'dest|cfg|pre|tight|nonTight|0.40'
  'dest|cfg|pre|tight|nonTight|0.40|false|2.0'
  'dest|cfg|pre||nonTight'
)
for index in "${!fanout_mutations[@]}"; do
  printf '%s\n' "${fanout_mutations[$index]}" > "${tmpdir}/bad-fanout-${index}.txt"
  if validate_five_field_fanout_contract unit "${tmpdir}/bad-fanout-${index}.txt" "${tmpdir}/valid.yaml" 0.40; then
    printf 'fanout schema mutation %s was not rejected\n' "$index" >&2
    exit 1
  fi
done

printf '%s\n%s\n' 'dest|cfg|pre|tight|nonTight' 'dest2|cfg2|pre2|tight2|nonTight2' > "${tmpdir}/two-row.fanout"
if validate_five_field_fanout_contract unit "${tmpdir}/two-row.fanout" "${tmpdir}/valid.yaml" 0.40; then
  printf 'two-row fanout ownership was not rejected\n' >&2
  exit 1
fi

if validate_five_field_fanout_contract unit "${tmpdir}/valid.fanout" "${tmpdir}/missing.yaml" 0.40; then
  printf 'missing materialized config was not rejected\n' >&2
  exit 1
fi

printf '%s\n' 'coneR: 0.30' > "${tmpdir}/wrong-cone.yaml"
if validate_five_field_fanout_contract unit "${tmpdir}/valid.fanout" "${tmpdir}/wrong-cone.yaml" 0.40; then
  printf 'wrong materialized cone was not rejected\n' >&2
  exit 1
fi

printf '%s\n%s\n' 'coneR: 0.40' 'coneR: 0.40' > "${tmpdir}/duplicate-cone.yaml"
if validate_five_field_fanout_contract unit "${tmpdir}/valid.fanout" "${tmpdir}/duplicate-cone.yaml" 0.40; then
  printf 'duplicate materialized cone authority was not rejected\n' >&2
  exit 1
fi

# The production controller's die helper exits.  Run negative attempt-state
# transitions in subshells so the test can assert those exact hard failures.
die() { printf 'EXPECTED_ATTEMPT_REJECTION: %s\n' "$*" >&2; exit 2; }

runtime_root="${tmpdir}/runtime-authority"
mkdir -p "$runtime_root"
runtime_authority_manifest="${runtime_root}/authority.json"
runtime_authority_fingerprint="${runtime_authority_manifest}.sha256"
frozen_sha="$(printf 'a%.0s' {1..64})"
calo_reco_build_receipt="/frozen/build_receipt.json"
RJ_THE134_CALO_RECO_BUILD_RECEIPT_SHA256="$frozen_sha"
calo_reco_source_manifest="/frozen/source_manifest.json"
RJ_THE134_CALO_RECO_SOURCE_MANIFEST_SHA256="$frozen_sha"
calo_reco_library="/frozen/libcalo_reco.so"
RJ_THE134_CALO_RECO_LIBRARY_SHA256="$frozen_sha"
release_core_lib_dir="/release/lib"
release_core_lib64_dir="/release/lib64"
release_calo_io="/release/lib64/libcalo_io.so"
RJ_THE134_RELEASE_CALO_IO_SHA256="$frozen_sha"
release_clusteriso="/release/lib64/libclusteriso.so"
RJ_THE134_RELEASE_CLUSTERISO_SHA256="$frozen_sha"
release_jetbase="/release/lib64/libjetbase.so"
RJ_THE134_RELEASE_JETBASE_SHA256="$frozen_sha"
pp_base_e_model="/frozen/model_base_E.root"
RJ_THE134_PP_BASE_E_MODEL_SHA256="$frozen_sha"
pinned_release_name="ana.560"
pinned_offline_main="/release/ana.560"
pinned_calo_reco_soname="libcalo_reco.so.0"
pp_period="0mrad"

write_runtime_authority_manifest ensure
write_runtime_authority_manifest verify
pp_period="1p5mrad"
if ( write_runtime_authority_manifest verify ) \
  >"${tmpdir}/runtime-period-drift.stdout" 2>"${tmpdir}/runtime-period-drift.stderr"; then
  printf 'p+p period drift was accepted by the frozen runtime authority\n' >&2
  exit 1
fi
pp_period="0mrad"
chmod a-w "$runtime_authority_manifest" "$runtime_authority_fingerprint" "$runtime_root"
write_runtime_authority_manifest ensure
chmod u+w "$runtime_root" "$runtime_authority_manifest" "$runtime_authority_fingerprint"
python3 - "$runtime_authority_manifest" <<'PY'
from pathlib import Path
import json
import sys

path = Path(sys.argv[1])
payload = json.loads(path.read_text())
payload["status"] = "DRIFT"
path.write_text(json.dumps(payload, sort_keys=True) + "\n")
PY
if ( write_runtime_authority_manifest ensure ) \
  >"${tmpdir}/runtime-drift.stdout" 2>"${tmpdir}/runtime-drift.stderr"; then
  printf 'runtime-authority content drift was accepted\n' >&2
  exit 1
fi
grep -Fq 'existing runtime authority manifest differs' \
  "${tmpdir}/runtime-drift.stderr"

split_root="${tmpdir}/runtime-authority-split"
mkdir -p "$split_root"
runtime_authority_manifest="${split_root}/authority.json"
runtime_authority_fingerprint="${runtime_authority_manifest}.sha256"
printf '{}\n' > "$runtime_authority_manifest"
if ( write_runtime_authority_manifest ensure ) \
  >"${tmpdir}/runtime-split.stdout" 2>"${tmpdir}/runtime-split.stderr"; then
  printf 'incomplete runtime-authority manifest/fingerprint pair was accepted\n' >&2
  exit 1
fi
grep -Fq 'runtime-authority manifest/fingerprint pair is incomplete' \
  "${tmpdir}/runtime-split.stderr"

attempt_base="${tmpdir}/attempt-reuse"
selection="$(select_materialization_attempt "$attempt_base")"
IFS=$'\t' read -r disposition attempt_dir <<< "$selection"
[[ "$disposition" == new && "$(basename "$attempt_dir")" == attempt_01 ]] || {
  printf 'first materialization attempt was not attempt_01\n' >&2
  exit 1
}
printf 'descriptor\n' > "${attempt_dir}/unit.sub"
seal_materialization_attempt "$attempt_dir" $'contract\tvalue'
selection="$(select_materialization_attempt "$attempt_base")"
IFS=$'\t' read -r disposition reused_dir <<< "$selection"
[[ "$disposition" == reuse && "$reused_dir" == "$attempt_dir" ]] || {
  printf 'sealed materialization attempt was not reused exactly\n' >&2
  exit 1
}
[[ "$(find "$attempt_base" -type f -name '*.sub' | wc -l | tr -d ' ')" == 1 ]] || {
  printf 'exact reuse duplicated a materialized submit descriptor\n' >&2
  exit 1
}

chmod u+w "${attempt_dir}/unit.sub"
if ( verify_materialization_attempt_seal "$attempt_dir" ) \
  2>"${tmpdir}/writable-attempt.stderr"; then
  printf 'writable materialization mutation was accepted\n' >&2
  exit 1
fi
grep -Fq 'materialization attempt contains writable state' \
  "${tmpdir}/writable-attempt.stderr"
chmod a-w "${attempt_dir}/unit.sub"
verify_materialization_attempt_seal "$attempt_dir"

incomplete_base="${tmpdir}/attempt-incomplete"
mkdir -p "${incomplete_base}/attempt_01"
selection="$(select_materialization_attempt "$incomplete_base")"
IFS=$'\t' read -r disposition attempt_dir <<< "$selection"
[[ "$disposition" == new && "$(basename "$attempt_dir")" == attempt_02 ]] || {
  printf 'an incomplete attempt was overwritten instead of preserved\n' >&2
  exit 1
}
printf 'descriptor\n' > "${attempt_dir}/unit.sub"
seal_materialization_attempt "$attempt_dir" $'contract\trepaired'
selection="$(select_materialization_attempt "$incomplete_base")"
IFS=$'\t' read -r disposition reused_dir <<< "$selection"
[[ "$disposition" == reuse && "$reused_dir" == "$attempt_dir" ]] || {
  printf 'sealed replacement was not reused with its older incomplete evidence preserved\n' >&2
  exit 1
}

exhausted_base="${tmpdir}/attempt-exhausted"
mkdir -p "${exhausted_base}/attempt_01" "${exhausted_base}/attempt_02" "${exhausted_base}/attempt_03"
if ( select_materialization_attempt "$exhausted_base" ) >/dev/null 2>&1; then
  printf 'three incomplete materialization attempts did not exhaust the retry budget\n' >&2
  exit 1
fi

mixed_base="${tmpdir}/attempt-mixed"
selection="$(select_materialization_attempt "$mixed_base")"
IFS=$'\t' read -r _disposition mixed_attempt <<< "$selection"
printf 'descriptor\n' > "${mixed_attempt}/unit.sub"
seal_materialization_attempt "$mixed_attempt" $'contract\tmixed'
chmod u+w "$mixed_base"
mkdir "${mixed_base}/attempt_02"
if ( select_materialization_attempt "$mixed_base" ) >/dev/null 2>&1; then
  printf 'sealed and incomplete attempts were accepted together\n' >&2
  exit 1
fi

python3 - \
  "$controller" \
  "${repo_root}/src/RecoilJets.cc" \
  "${repo_root}/scripts/sdcc/runtime/condor/RecoilJets_Condor_submit.sh" <<'PY'
from pathlib import Path
import re
import sys

expansion = '"${materialize_contract_env[@]}"'
source = Path(sys.argv[1]).read_text()
analysis_source = Path(sys.argv[2]).read_text()
submitter_source = Path(sys.argv[3]).read_text()

def validate(text: str) -> None:
    required = (
        'local -a materialize_contract_env=(',
        'RJ_REPLAY_FOUNDATION_CANARY=1',
        'RJ_REPLAY_FOUNDATION_CAPACITY_CANARY=1',
        'RJ_REPLAY_FOUNDATION_CAPACITY_CANARY_ID="$capacity_canary_id"',
        'RJ_REPLAY_FOUNDATION_EXECUTION_PARTITION_SHA256',
        'readonly source_occurrences_per_output="1"',
        'capacity-preflight|capacity-submit|capacity-resume-submit|capacity-status|capacity-validate',
        '"$submitter" "$dataset" condorDoAllSmoke groupSize "$execution_group_size" maxJobs 1',
        '"$execution_row_count" "$execution_group_size"',
        '"$source_occurrences_per_output" "$capacity_mode"',
        '"source_occurrences_per_output": expected_source_count',
        'root_certificate.get("source_occurrences_per_output", -1)',
        'ANALYSIS_MINIMUM_BYTES = 50_000',
        'SIDECAR_GROSS_TRUNCATION_FLOOR_BYTES = 4_096',
        '"name": "analysis_writer_root_v1"',
        '"name": "photon_training_multiview_v1"',
        '"minimum_bytes_mode": "diagnostic_only"',
        'minimum_bytes=ANALYSIS_MINIMUM_BYTES',
        'minimum_bytes=0',
        'root_certificate.get("analysis_health_profile")',
        'root_certificate.get("sidecar_health_profile")',
        '"POPULATED" if candidate_definitions else "VALID_EMPTY"',
        '"populated_rows": populated_rows',
        '"valid_empty_rows": valid_empty_rows',
        '"capacity_population_witness:no populated selected row"',
        'root_certificate.get("populated_row_count", -1)',
        'root_certificate.get("valid_empty_row_count", -1)',
        '"source_identity_canonical_sha256": source_identity_sha256',
        '"source_occurrence_id_hex": (',
        '"analysis_sha256": file_sha256(analysis_path)',
        '"sidecar_sha256": file_sha256(sidecar_path)',
        'validate_system_matrix_or_capacity_audit',
        'RJ_THE134_CAPACITY_FILE_VALIDATOR_SHA256',
        'RJ_THE134_CAPACITY_CONTRACT_VALIDATOR_SHA256',
        'RJ_THE134_CAPACITY_CONTROLLER_SHA256',
        'RJ_THE134_CAPACITY_VALIDATION_COMMIT',
        'RJ_THE134_CAPACITY_SUBMISSION_MANIFEST_SHA256',
        'RJ_THE134_CAPACITY_SUBMISSION_RECEIPT_SHA256',
        'RJ_THE134_CAPACITY_SUBMISSION_JOURNAL_SHA256',
        'RJ_THE134_CAPACITY_RUNTIME_AUTHORITY_SHA256',
        'validate_capacity_postrun_authority',
        'write_runtime_authority_manifest verify',
        'root_health_identity_join_certificate_capacity_v6.json',
        'pp_capacity_multiview_audit_v2.json',
        'auau_capacity_multiview_audit_v2.json',
        '"validation_authority": validation_authority',
        'submitted argument identity differs from the sealed receipt',
        'source sample lacks the frozen run prefix',
        '"source_execution_contract": source_execution',
        'capacity source execution contract differs from the sealed root/receipt authority',
        'capacity validation forbids an alternate matrix-preparer path',
        '"file_validator": {',
        '"schema": "THE134_CAPACITY_MULTIVIEW_MATRIX_AUDIT_V1"',
        '"scope": "capacity"',
        '"class_balance_state": "NOT_APPLICABLE_CAPACITY_NON_TRAINING"',
        '"training_matrix_authority_earned": False',
        '"matrix_materialized": False',
        'validator.validate_artifact_metadata(',
        'validator.validate_file_arrays(',
        'existing capacity audit differs; preserve it before retry',
        'existing capacity resource certificate differs; preserve it before retry',
        'def resolve_classad_memory_usage_mb',
        '/Expr(((ResidentSetSize + 1023) / 1024))/',
        '"EVALUATED_FROM_RESIDENT_SET_SIZE_KB"',
        '"memory_usage_resolution": memory_usage_resolution',
        '"multiview_audits": audit_reports',
        '"$capacity_pp_multiview_audit" "$capacity_auau_multiview_audit"',
        'RJ_REPLAY_LANE="$lane"',
        'RJ_REPLAY_SCHEMA_SHA256="$RJ_THE134_REPLAY_SCHEMA_SHA256"',
        'validate_pp_sim_weight_contract "$pp_period"',
        'period="$pp_period"',
        'RJ_PPG12_PHOTON_YIELD=1',
        'RJ_PPG12_PERIOD="$pp_period"',
        'RJ_PINNED_CALO_RECO_RELEASE_COMPANIONS=1',
        'RJ_PINNED_CALO_RECO_SONAME="$pinned_calo_reco_soname"',
        'RJ_PINNED_RELEASE_NAME="$pinned_release_name"',
        'RJ_PINNED_OFFLINE_MAIN="$pinned_offline_main"',
        'RJ_PINNED_RELEASE_CALO_IO_PATH="$release_calo_io"',
        'RJ_PINNED_RELEASE_CLUSTERISO_PATH="$release_clusteriso"',
        'RJ_PINNED_RELEASE_JETBASE_PATH="$release_jetbase"',
        'RJ_RELEASE_CORE_LIB_DIR="$release_core_lib_dir"',
        'RJ_RELEASE_CORE_LIB64_DIR="$release_core_lib64_dir"',
        'RJ_CONDOR_SEALED_ENVIRONMENT=1',
        'env -u RJ_FORCE_RELEASE_CORE_LIBS -u RJ_FORCE_RELEASE_CALO_IO -u RJ_RELEASE_CALO_IO_PATH',
    )
    for token in required:
        if token not in text:
            raise ValueError(f"missing submit-shell replay-canary identity: {token}")
    helper_start = text.index("def resolve_classad_memory_usage_mb")
    helper_end = text.index("\n\nreports = []", helper_start)
    namespace = {}
    exec(text[helper_start:helper_end], namespace)
    resolve_memory = namespace["resolve_classad_memory_usage_mb"]
    expression = "/Expr(((ResidentSetSize + 1023) / 1024))/"
    if resolve_memory(
        {"MemoryUsage": expression, "ResidentSetSize": 1_250_000}
    ) != (1221, "EVALUATED_FROM_RESIDENT_SET_SIZE_KB"):
        raise ValueError("Condor expression-valued MemoryUsage resolution drifted")
    if resolve_memory(
        {"MemoryUsage": 1221, "ResidentSetSize": 1_250_000}
    ) != (1221, "NUMERIC_CLASSAD"):
        raise ValueError("numeric MemoryUsage resolution drifted")
    for invalid in (
        {"MemoryUsage": "/Expr(0)/", "ResidentSetSize": 1_250_000},
        {"MemoryUsage": expression, "ResidentSetSize": 0},
    ):
        try:
            resolve_memory(invalid)
        except ValueError:
            pass
        else:
            raise ValueError(f"invalid MemoryUsage fixture was accepted: {invalid}")
    if text.count(expansion) != 2:
        raise ValueError(
            "the shared submit-shell replay-canary identity must guard exactly "
            "the p+p and Au+Au dry-materialization calls"
        )
    worker_contract = (
        'RJ_REPLAY_FOUNDATION_V1=1;RJ_REPLAY_FOUNDATION_CANARY=1;'
        'RJ_REPLAY_TRACE=0;RJ_REPLAY_LANE=${lane}'
    )
    if worker_contract not in text:
        raise ValueError("worker descriptor replay-canary identity was lost")
    receipt_contract = (
        'materialized_config_values=()',
        'descriptor_env_values "$submit_file" RJ_CONFIG_YAML',
        '"${#materialized_config_values[@]}" == 1',
        'materialized_config_sha256',
        'file_sha256(materialized_config) != materialized_config_sha',
        'descriptor_env_values "$submit_file" RJ_SIM_ALLOW_NONE_LISTS',
        'p+p descriptor must reject NONE lists',
        'Au+Au descriptor must authorize only its typed optional MBD list',
        'RJ_SIM_ALLOW_NONE_LISTS=0',
        'RJ_SIM_ALLOW_NONE_LISTS=1',
        'sealed_getenv="$(condor_field "$submit_file" getenv)"',
        '"$sealed_getenv" == False',
        'descriptor must disable submit-host environment inheritance',
        'require_descriptor_env_exact "$row_id" "$submit_file" RJ_PPG12_PHOTON_YIELD 1',
        'require_descriptor_env_exact "$row_id" "$submit_file" RJ_PPG12_PHOTON_YIELD_DOUBLE 0',
        'require_descriptor_env_exact "$row_id" "$submit_file" RJ_PPG12_PERIOD "$pp_period"',
        'require_descriptor_env_exact "$row_id" "$submit_file" RJ_PPG12_PERIOD_USE_LUMI_WEIGHT 1',
        'require_descriptor_env_exact "$row_id" "$submit_file" RJ_PPG12_PERIOD_STRICT_DI 0',
        'require_descriptor_env_exact "$row_id" "$submit_file" RJ_PPG12_PERIOD_ALLOW_ALL_SIM 0',
        'require_descriptor_env_exact "$row_id" "$submit_file" RJ_PPG12_PERIOD_ALLOW_MIX_OVERRIDE 0',
        'require_descriptor_env_exact "$row_id" "$submit_file" RJ_PPG12_PERIOD_ALLOW_VERTEX_FILE_OVERRIDE 0',
        'require_descriptor_env_absent "$row_id" "$submit_file" RJ_PPG12_PHOTON_YIELD_MIX_WEIGHT',
        'require_descriptor_env_absent "$row_id" "$submit_file" RJ_PP_VERTEX_REWEIGHT_FILE',
    )
    for token in receipt_contract:
        if token not in text:
            raise ValueError(f"missing materialized-config receipt protection: {token}")
    single_provider_contract = (
        'readonly pinned_offline_main="/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_ana/ana.560"',
        '"$release_core_lib_dir" == "${pinned_offline_main}/lib"',
        '"$release_core_lib64_dir" == "${pinned_offline_main}/lib64"',
        'resolve_release_companion RELEASE_CALO_IO libcalo_io.so',
        'resolve_release_companion RELEASE_CLUSTERISO libclusteriso.so',
        'resolve_release_companion RELEASE_JETBASE libjetbase.so',
        'validate_calo_reco_build_authority',
        'THE134_ANA560_CALORECO_BUILD_RECEIPT_V3',
        'abi.get("rpath_runpath_exact_match") is not True',
        'THE134_ANA560_CALORECO_SOURCE_MANIFEST_V2',
        'readonly pinned_calo_reco_soname="libcalo_reco.so.0"',
        'snapshotted CaloReco SONAME alias does not resolve to its one provider',
        'RJ_THE134_CALO_RECO_BUILD_RECEIPT_SHA256',
        'RJ_THE134_CALO_RECO_SOURCE_MANIFEST_SHA256',
        'controller environment must not inherit RJ_FORCE_RELEASE_CORE_LIBS',
        'controller environment must not inherit RJ_FORCE_RELEASE_CALO_IO',
        'controller environment must not inherit RJ_RELEASE_CALO_IO_PATH',
        'snapshot illegally duplicates release-owned',
        'frozen executor lacks the exact ana.560 loader suffix',
        'frozen executor lacks the exact ${pinned_release_name} runtime witness',
        'validate_snapshot_macro_provider',
        'Calo_Calib macro must load exactly one absolute snapshotted CaloReco provider',
        'Calo_Calib macro does not load the same snapshotted CaloReco provider',
        'verify_sealed_snapshot_receipts',
        'snapshot_loader_receipt_sha256',
        'snapshot_manifest_sha256',
        'snapshot symlink must be relative',
        'snapshot symlink parent is writable',
        'snapshot symlink target is writable',
        'writable frozen snapshot root survived seal',
        'write_runtime_authority_manifest ensure',
        'write_runtime_authority_manifest verify',
        'runtime-authority manifest/fingerprint pair is incomplete',
        'existing runtime authority manifest differs from current frozen authority',
        '"schema": "THE134_SINGLE_PROVIDER_RUNTIME_AUTHORITY_V3"',
        '"interaction": "SI"',
        '"mix_weight": "period_auto"',
        '"vertex_reweight": "period_auto"',
        'verify_materialization_attempt_seal',
        'select_materialization_attempt',
        'seal_materialization_attempt',
        'materialization retry budget exhausted with three preserved incomplete attempts',
        'sealed materialization contract differs from exact readback',
        'resume ana.560 release-companion authority differs from frozen manifest',
        "NF != 46 {exit 1}",
        "NF != 28 {exit 1}",
    )
    for token in single_provider_contract:
        if token not in text:
            raise ValueError(f"missing single-provider protection: {token}")
    if text.count('RJ_PINNED_CALO_RECO_RELEASE_COMPANIONS=1') != 2:
        raise ValueError("single-provider mode must guard exactly the p+p and Au+Au materializers")
    if text.count('RJ_PINNED_CALO_RECO_SONAME="$pinned_calo_reco_soname"') != 2:
        raise ValueError("CaloReco SONAME authority must guard both materializers")
    if text.count('RJ_CONDOR_SEALED_ENVIRONMENT=1') != 2:
        raise ValueError("sealed Condor environment mode must guard both materializers")
    if text.count(
        'env -u RJ_FORCE_RELEASE_CORE_LIBS -u RJ_FORCE_RELEASE_CALO_IO '
        '-u RJ_RELEASE_CALO_IO_PATH'
    ) != 2:
        raise ValueError("stale runtime overrides must be stripped at both materializers")
    if text.count("\n  validate_calo_reco_build_authority ||\n") != 1:
        raise ValueError("CaloReco build authority must be invoked exactly once in preflight")
    submit_all_block = text.split("submit_all() {", 1)[1].split(
        "resume_submit() {", 1
    )[0]
    if not re.search(
        r"submit_all\(\) \{\n  preflight\n  assert_fresh_submission\n",
        "submit_all() {" + submit_all_block,
    ):
        raise ValueError("submit must perform idempotent preflight before freshness checks")

    manifest_block = text.split("write_submission_manifest() {", 1)[1].split(
        "write_runtime_authority_manifest() {", 1
    )[0]
    row_formats = re.findall(r"printf '([^']*%s[^']*)'", manifest_block)
    row_format = next((value for value in row_formats if value.count("%s") == 45), None)
    if row_format is None:
        raise ValueError("submission-manifest row format was not found")
    if row_format.count(r"\t") != 45:
        raise ValueError("46-field manifest format must contain 45 substitutions and 1 literal")

    receipt_block = text.split("submit_row() {", 1)[1].split(
        "assert_fresh_submission() {", 1
    )[0]
    receipt_formats = re.findall(r"printf '([^']*%s[^']*)'", receipt_block)
    receipt_format = next(
        (value for value in receipt_formats if value.count("%s") == 28), None
    )
    if receipt_format is None:
        raise ValueError("submission-receipt row format was not found")
    if receipt_format.count(r"\t") != 27:
        raise ValueError("submission receipt format must contain exactly 28 fields")

def validate_si_auto_weight_authority(analysis_text: str, submitter_text: str) -> None:
    analysis_contract = (
        'const bool mixWeightExplicit = (std::getenv("RJ_PPG12_PHOTON_YIELD_MIX_WEIGHT") != nullptr);',
        'const bool vertexFileExplicit = (std::getenv("RJ_PP_VERTEX_REWEIGHT_FILE") != nullptr);',
        'm_ppg12PhotonYieldDoubleInteraction ? m_ppg12PeriodFDouble : m_ppg12PeriodFSingle;',
        'if (!mixWeightExplicit)',
        'm_ppg12PhotonYieldMixWeight = expectedMix;',
        'm_ppg12PeriodMixWeightAuto = true;',
        'if (!vertexFileExplicit)',
        'm_vertexReweightFile = m_ppg12PeriodExpectedVertexFile;',
        'm_vertexReweightOn = true;',
        'm_ppg12PeriodVertexFileAuto = true;',
        '!m_ppg12PeriodContractEnabled ||',
        '!m_ppg12PeriodUseLumiWeight ||',
        '!m_vertexReweightOn ||',
        '!m_vertexReweightH)',
    )
    for token in analysis_contract:
        if token not in analysis_text:
            raise ValueError(f"missing ordinary SI automatic-weight authority: {token}")
    submitter_contract = (
        'extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_PHOTON_YIELD)"',
        'extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_PERIOD)"',
        'extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_PERIOD_USE_LUMI_WEIGHT)"',
        'extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_PERIOD_ALLOW_MIX_OVERRIDE)"',
        'extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_PERIOD_ALLOW_VERTEX_FILE_OVERRIDE)"',
        'extra="$(append_submit_extra_env_var "$extra" RJ_PP_VERTEX_REWEIGHT_FILE)"',
    )
    for token in submitter_contract:
        if token not in submitter_text:
            raise ValueError(f"missing submitter propagation authority: {token}")

validate(source)
validate_si_auto_weight_authority(analysis_source, submitter_source)

# Deliberately mutate one call site and prove this validator rejects the drift.
mutated = source.replace(expansion, "env", 1)
try:
    validate(mutated)
except ValueError:
    pass
else:
    raise SystemExit("one-call-site mutation was not rejected")

mutated = source.replace('RJ_PINNED_CALO_RECO_RELEASE_COMPANIONS=1', '', 1)
try:
    validate(mutated)
except ValueError:
    pass
else:
    raise SystemExit("single-provider mutation was not rejected")

mutated = source.replace("NF != 46 {exit 1}", "NF != 45 {exit 1}", 1)
try:
    validate(mutated)
except ValueError:
    pass
else:
    raise SystemExit("46-field manifest mutation was not rejected")

mutated = source.replace(
    "\n  validate_calo_reco_build_authority ||\n",
    "\n  : ||\n",
    1,
)
try:
    validate(mutated)
except ValueError:
    pass
else:
    raise SystemExit("CaloReco build-authority mutation was not rejected")

mutated = source.replace('RJ_CONDOR_SEALED_ENVIRONMENT=1', '', 1)
try:
    validate(mutated)
except ValueError:
    pass
else:
    raise SystemExit("sealed Condor environment mutation was not rejected")

mutated = source.replace('write_runtime_authority_manifest ensure', '', 1)
try:
    validate(mutated)
except ValueError:
    pass
else:
    raise SystemExit("runtime-authority ensure mutation was not rejected")

mutated = source.replace('"$sealed_getenv" == False', '"$sealed_getenv" == True', 1)
try:
    validate(mutated)
except ValueError:
    pass
else:
    raise SystemExit("descriptor getenv mutation was not rejected")

mutated = source.replace(
    '"schema": "THE134_CAPACITY_MULTIVIEW_MATRIX_AUDIT_V1"',
    '"schema": "THE134_FACTORIAL_VIEW_TRAINING_MATRIX_AUDIT_V1"',
    1,
)
try:
    validate(mutated)
except ValueError:
    pass
else:
    raise SystemExit("capacity-only audit schema mutation was not rejected")

mutated = source.replace(
    '"matrix_materialized": False',
    '"matrix_materialized": True',
    1,
)
try:
    validate(mutated)
except ValueError:
    pass
else:
    raise SystemExit("capacity matrix-materialization mutation was not rejected")

mutated = source.replace(
    'RJ_PPG12_PHOTON_YIELD=1 \\\n',
    '',
    1,
)
try:
    validate(mutated)
except ValueError:
    pass
else:
    raise SystemExit("p+p photon-yield weight-contract mutation was not rejected")

mutated = source.replace(
    'RJ_PPG12_PERIOD="$pp_period" \\\n',
    '',
    1,
)
try:
    validate(mutated)
except ValueError:
    pass
else:
    raise SystemExit("p+p period weight-contract mutation was not rejected")

mutated_analysis = analysis_source.replace(
    'm_ppg12PeriodMixWeightAuto = true;',
    'm_ppg12PeriodMixWeightAuto = false;',
    1,
)
try:
    validate_si_auto_weight_authority(mutated_analysis, submitter_source)
except ValueError:
    pass
else:
    raise SystemExit("automatic SI mix-weight authority mutation was not rejected")

mutated_analysis = analysis_source.replace(
    'm_ppg12PeriodVertexFileAuto = true;',
    'm_ppg12PeriodVertexFileAuto = false;',
    1,
)
try:
    validate_si_auto_weight_authority(mutated_analysis, submitter_source)
except ValueError:
    pass
else:
    raise SystemExit("automatic SI vertex-weight authority mutation was not rejected")

print("THE134_SUBMITTER_CANARY_IDENTITY_WIRING_PASS guarded_calls=2 mutations_rejected=11 tuple_mutations=3 fanout_mutations=8 pp_period_mutations=3 descriptor_weight_mutations=2 science_authority_mutations=2 runtime_authority_transitions=6 materialization_state_transitions=6")
PY
