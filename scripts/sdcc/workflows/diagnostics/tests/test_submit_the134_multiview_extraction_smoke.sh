#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd -P)"
controller="${repo_root}/scripts/sdcc/workflows/diagnostics/submit_the134_multiview_extraction_smoke.sh"

bash -n "$controller"

tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/the134-tuple-contract.XXXXXX")"
trap 'chmod -R u+w "$tmpdir" 2>/dev/null || true; rm -rf "$tmpdir"' EXIT
die() { return 2; }
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
        "schema": "THE134_FULL_EXTRACTION_PARTITION_CONTRACT_V1",
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
    "schema": "THE134_FULL_MULTIVIEW_EXTRACTION_PLAN_V1",
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
    "schema": "THE134_FULL_MULTIVIEW_EXTRACTION_PREFLIGHT_RECEIPT_V1",
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
        'capacity-preflight|capacity-submit|capacity-resume-submit|capacity-status|capacity-validate',
        '"$submitter" "$dataset" condorDoAllSmoke groupSize "$execution_group_size" maxJobs 1',
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
        '"schema": "THE134_SINGLE_PROVIDER_RUNTIME_AUTHORITY_V2"',
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
    if text.count("\n  validate_calo_reco_build_authority\n") != 1:
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
    "\n  validate_calo_reco_build_authority\n",
    "\n  :\n",
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

print("THE134_SUBMITTER_CANARY_IDENTITY_WIRING_PASS guarded_calls=2 mutations_rejected=9 tuple_mutations=3 fanout_mutations=8 pp_period_mutations=3 descriptor_weight_mutations=2 science_authority_mutations=2 runtime_authority_transitions=6 materialization_state_transitions=6")
PY
