#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
driver="$(cd "${script_dir}/.." && pwd -P)/submit_ppg12_stitched_purity_photon_canaries.sh"
repo_root="$(cd "${script_dir}/../../../../.." && pwd -P)"

tmp="$(mktemp -d "${TMPDIR:-/tmp}/ppg12-photon-canary-test.XXXXXX")"
trap 'rm -rf "$tmp"' EXIT

fail() { printf 'FAIL: %s\n' "$*" >&2; exit 1; }
pass() { printf 'PASS: %s\n' "$*"; }

mock_submitter="${tmp}/mock_submitter.sh"
cat > "$mock_submitter" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
{
  printf 'CALL\t'
  printf '%q ' "$@"
  printf '\n'
  for key in \
    RJ_CONFIG_YAML RJ_DEST_BASE_OVERRIDE RJ_SUBMISSION_NAMESPACE \
    RJ_PPG12_CLOSURE_CANARY RJ_PPG12_CLOSURE_CANARY_ID \
    RJ_PPG12_PPSIM_REPLAY_SEEDS RJ_PPG12_PPSIM_EXPECT_PEDESTAL_SEQUENCE \
    RJ_DISABLE_JES_CDB_AUDIT RJ_REQUIRE_SIM_GLOBAL \
    RJ_PPG12_PHOTON_YIELD RJ_PPG12_PHOTON_YIELD_CLUSTER_ERES \
    RJ_PPG12_PERIOD RJ_PPG12_CROSSING_PERIOD \
    RJ_PPG12_PERIOD_USE_LUMI_WEIGHT RJ_PPG12_PERIOD_ALLOW_ALL_SIM \
    RJ_PPG12_PERIOD_ALLOW_MIX_OVERRIDE RJ_PPG12_PERIOD_ALLOW_VERTEX_FILE_OVERRIDE \
    RJ_PPG12_TABLE_QA RJ_PPG12_TABLE_QA_NPB_DATA_TAGGING \
    RJ_PPG12_FIG13_PARITY_QA RJ_PPG12_FIG11_SB_DIAGNOSTIC \
    RJ_PP_PHOTONID_TRAINING_TREE RJ_PP_PHOTONID_TRAINING_TREE_MAX_ENTRIES \
    RJ_PP_PHOTONID_SOURCE_ROLE RJ_PP_PHOTONID_PPG12_FILTER \
    RJ_PP_PHOTONID_REQUIRE_PRESELECTION RJ_PHOTON_ID_ROW_MATCH \
    RJ_DISABLE_ID_FANOUT RJ_ID_FANOUT_MAX_ROWS RJ_DISABLE_JET_PT_INTERNALIZATION \
    RJ_DISABLE_DPHI_INTERNALIZATION RJ_DIRECT_DST_DOALL RJ_DIRECT_NEVENTS \
    RJ_VALIDATE_SIM_INPUT_MAX_LINES \
    RJ_INTERNAL_JET_PT_MINS RJ_INTERNAL_DPHI_PI_FRACTIONS RJ_REQUEST_MEMORY_MB \
    RJ_AUTO_MERGE RJ_REQUIRE_NON_TINY_OUTPUT RJ_MIN_OUTPUT_BYTES \
    RJ_FAIL_ON_MISSING_CALO_INPUT \
    RJ_SIM_SIGNAL_SAMPLE_SET RJ_PPG12_PHOTON_YIELD_DOUBLE \
    RJ_PPG12_PERIOD_STRICT_DI RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4 \
    RJ_PPG12_PPSIM_G4_ONLY RJ_SIM_ALLOW_NONE_LISTS \
    RJ_CODEX_CHAT_NAME RJ_CODEX_THREAD_ID
  do
    if [[ -n "${!key+x}" ]]; then
      printf 'ENV\t%s\t%s\n' "$key" "${!key}"
    else
      printf 'ABSENT\t%s\n' "$key"
    fi
  done
  printf 'END\n'
} >> "$MOCK_CALL_LOG"
EOF
chmod +x "$mock_submitter"

campaign="ppg12_photon_canary_unit"
evidence="${tmp}/evidence"
output_root="${tmp}/output"
common_env=(
  "RJ_PPG12_PHOTON_CANARY_CAMPAIGN_TAG=${campaign}"
  "RJ_PPG12_PHOTON_CANARY_EVIDENCE_DIR=${evidence}"
  "RJ_PPG12_PHOTON_CANARY_OUTPUT_ROOT=${output_root}"
  "RJ_PPG12_PHOTON_CANARY_SUBMITTER=${mock_submitter}"
  "RJ_PPG12_PHOTON_CANARY_CONFIG_YAML=${repo_root}/macros/analysis_config.yaml"
)

# Default invocation is a non-mutating plan and must not call the submitter.
env "${common_env[@]}" MOCK_CALL_LOG="${tmp}/calls.log" "$driver" > "${tmp}/plan.out"
[[ ! -e "${tmp}/calls.log" ]] || fail "plan mode invoked the canonical submitter"
[[ -s "${evidence}/photon_canary_plan.json" ]] || fail "plan JSON was not emitted"

python3 - "${evidence}/photon_canary_plan.json" <<'PY'
import json
import sys

with open(sys.argv[1]) as stream:
    doc = json.load(stream)

assert doc["schema"] == "ppg12-stitched-purity-photon-canary-plan/v4"
assert doc["status"] == "PLANNED"
assert doc["bounded_contract"]["lane_count"] == 12
assert doc["bounded_contract"]["group_size"] == 5
assert doc["bounded_contract"]["max_jobs_per_lane"] == 1
assert doc["bounded_contract"]["input_path_validation_rows"] == 5
assert doc["bounded_contract"]["automatic_merge"] is False
assert doc["bounded_contract"]["candidate_tree_max_entries"] == 0
assert doc["bounded_contract"]["photon_id_row"] == ["newPPG12"] * 3
assert doc["bounded_contract"]["classification_energy"] == "unsmeared_calibrated_et"
assert doc["bounded_contract"]["legacy_multiplicative_cluster_smearing"] == "disabled"
assert doc["bounded_contract"]["response_smearing"] == "ppg12_additive_response_only"
assert doc["bounded_contract"]["rebuild_calo_from_g4"] == "all lanes"
assert doc["bounded_contract"]["g4_only"] == "all lanes"
assert doc["bounded_contract"]["rng_contract"] == "historical_fifo_replay_v2"
assert doc["bounded_contract"]["replay_seed_sequence"] == [
    2991264730, 4256268992, 2394322166, 874466025, 2240380304
]
assert doc["bounded_contract"]["expected_pedestal_sequence"] == 534
assert doc["bounded_contract"]["estimator_random_seed"] == 42
assert doc["bounded_contract"]["estimator_toy_count"] == 20000
assert doc["bounded_contract"]["jes_cdb_audit_disabled"] is True
assert doc["bounded_contract"]["random_seed_policy"] == (
    "identical_historical_fifo_replay_for_both_executables_v2"
)
assert doc["bounded_contract"]["si_five_column_input_graph"] == [
    "NONE", "G4Hits", "DST_TRUTH_JET", "NONE", "NONE"
]
assert doc["bounded_contract"]["di_five_column_input_graph"] == [
    "NONE", "G4Hits", "DST_TRUTH_JET", "NONE", "NONE"
]
assert doc["bounded_contract"]["si_registered_input_graph"] == [
    "G4Hits", "DST_TRUTH_JET"
]
assert doc["bounded_contract"]["di_registered_input_graph"] == [
    "G4Hits", "DST_TRUTH_JET"
]
assert doc["bounded_contract"]["dst_global_registered"] is False

lanes = doc["lanes"]
assert len(lanes) == 12
assert len({lane["lane_id"] for lane in lanes}) == 12
assert len({lane["output_base"] for lane in lanes}) == 12
assert len({lane["submission_namespace"] for lane in lanes}) == 12
assert {lane["sample_key"] for lane in lanes} == {"photon5", "photon10", "photon20"}
assert {lane["period"] for lane in lanes} == {"0mrad", "1p5mrad"}
assert {lane["interaction"] for lane in lanes} == {"si", "di"}
for lane in lanes:
    suffix = "_double" if lane["interaction"] == "di" else ""
    number = lane["sample_key"].removeprefix("photon")
    assert lane["recoil_sample"] == f"run28_photonjet{number}{suffix}"
    assert lane["ppg12_oracle_root"].endswith(f"/photon{number}{suffix}/bdt_split.root")
    assert lane["ppg12_tree"] == "slimtree"
    assert lane["recoil_tree"] == "AuAuPhotonIDTrainingTree"
    assert lane["replay_seed_sequence"] == [
        2991264730, 4256268992, 2394322166, 874466025, 2240380304
    ]
    assert lane["expected_pedestal_sequence"] == 534
    assert lane["jes_cdb_audit_disabled"] is True
    assert lane["five_column_input_graph"] == [
        "NONE", "G4Hits", "DST_TRUTH_JET", "NONE", "NONE"
    ]
    if lane["interaction"] == "di":
        assert lane["source_contract"] == (
            "run28_double_none_g4_truthjet_none_none_g4_only_rebuild"
        )
        assert lane["registered_input_graph"] == ["G4Hits", "DST_TRUTH_JET"]
    else:
        assert lane["source_contract"] == (
            "run28_single_interaction_none_g4_truthjet_none_none_g4_only_rebuild"
        )
        assert lane["registered_input_graph"] == ["G4Hits", "DST_TRUTH_JET"]
    assert lane["dst_global_registered"] is False
    assert lane["submit_argv"] == [
        "isSim", "condorDoAllDirect", "groupSize", "5", "maxJobs", "1",
        f"SAMPLE={lane['recoil_sample']}",
    ]
PY
pass "plan emits exact 12-lane oracle mapping without submission"

# Submission needs both independent approvals: --submit and exact token.
if env "${common_env[@]}" MOCK_CALL_LOG="${tmp}/calls.log" \
  RJ_CODEX_CHAT_NAME=test RJ_CODEX_THREAD_ID=test \
  "$driver" --token "SUBMIT_${campaign}" > "${tmp}/missing-submit.out" 2>&1; then
  fail "submit succeeded without --submit"
fi
if env "${common_env[@]}" MOCK_CALL_LOG="${tmp}/calls.log" \
  RJ_CODEX_CHAT_NAME=test RJ_CODEX_THREAD_ID=test \
  "$driver" --submit --token WRONG > "${tmp}/bad-token.out" 2>&1; then
  fail "submit succeeded with the wrong token"
fi
[[ ! -e "${tmp}/calls.log" ]] || fail "failed authorization reached the submitter"
pass "submission authorization fails closed"

# Every deterministic/legacy control is driver-owned and must fail closed if
# inherited from the caller, even when the inherited value happens to match.
for conflict in \
  RANDOMSEED=7 \
  RJ_PPG12_CLOSURE_CANARY=1 \
  RJ_PPG12_CLOSURE_CANARY_ID=stale \
  RJ_PPG12_PPSIM_REPLAY_SEEDS=1,2,3,4,5 \
  RJ_PPG12_PPSIM_EXPECT_PEDESTAL_SEQUENCE=1 \
  RJ_PPG12_PPSIM_FIXED_RANDOMSEED=42 \
  RJ_PPG12_PPSIM_FIXED_PEDESTAL_SEQUENCE=1482 \
  RJ_DISABLE_JES_CDB_AUDIT=1 \
  RJ_REQUIRE_SIM_GLOBAL=1 \
  RJ_PPG12_PEDESTAL_OVERRIDE=pedestal.root \
  RJ_PPG12_DI_ARCHIVED_EXPECT_PEDESTAL=pedestal.root
do
  key="${conflict%%=*}"
  if env "${common_env[@]}" MOCK_CALL_LOG="${tmp}/calls.log" "$conflict" \
    "$driver" plan > "${tmp}/conflict.out" 2>&1; then
    fail "plan accepted inherited control ${key}"
  fi
  grep -Fq "$key" "${tmp}/conflict.out" || fail "conflict rejection did not name ${key}"
done
if env "${common_env[@]}" MOCK_CALL_LOG="${tmp}/calls.log" \
  RJ_SUBMIT_EXTRA_ENV='RJ_PPG12_PPSIM_REPLAY_SEEDS=1,2,3,4,5' \
  "$driver" plan > "${tmp}/extra-conflict.out" 2>&1; then
  fail "plan accepted a driver-owned control through RJ_SUBMIT_EXTRA_ENV"
fi
grep -Fq 'RJ_PPG12_PPSIM_REPLAY_SEEDS' "${tmp}/extra-conflict.out" || \
  fail "RJ_SUBMIT_EXTRA_ENV rejection did not name the conflicting control"
pass "inherited deterministic and legacy controls are rejected"

call_log="${tmp}/calls.log"
env "${common_env[@]}" MOCK_CALL_LOG="$call_log" \
  RJ_CODEX_CHAT_NAME='THE-97 | unit test' RJ_CODEX_THREAD_ID='unit-thread' \
  "$driver" --submit --token "SUBMIT_${campaign}" > "${tmp}/submit.out"

[[ "$(grep -c '^CALL' "$call_log")" -eq 12 ]] || fail "expected exactly 12 submitter calls"
[[ "$(grep -c '^END' "$call_log")" -eq 12 ]] || fail "submitter call records are incomplete"
[[ "$(wc -l < "${evidence}/submission_receipts.tsv" | tr -d ' ')" -eq 13 ]] || fail "receipt table is not 12 lanes plus header"

python3 - "$call_log" <<'PY'
import sys

blocks = []
current = None
for raw in open(sys.argv[1]):
    line = raw.rstrip("\n")
    if line.startswith("CALL\t"):
        current = {"call": line.split("\t", 1)[1], "env": {}, "absent": set()}
    elif line == "END":
        blocks.append(current)
        current = None
    elif line.startswith("ENV\t"):
        _, key, value = line.split("\t", 2)
        current["env"][key] = value
    elif line.startswith("ABSENT\t"):
        current["absent"].add(line.split("\t", 1)[1])

assert len(blocks) == 12
seen = set()
for block in blocks:
    call = block["call"]
    env = block["env"]
    assert call.startswith("isSim condorDoAllDirect groupSize 5 maxJobs 1 SAMPLE=run28_photonjet")
    sample = call.split("SAMPLE=", 1)[1].strip()
    period = env["RJ_PPG12_PERIOD"]
    interaction = "di" if sample.endswith("_double") else "si"
    seen.add((sample, period, interaction))
    assert env["RJ_PPG12_CROSSING_PERIOD"] == period
    assert env["RJ_PPG12_CLOSURE_CANARY"] == "1"
    assert env["RJ_PPG12_CLOSURE_CANARY_ID"] == (
        f"photon:photon{sample.split('photonjet', 1)[1].removesuffix('_double')}:{period}:{interaction}"
    )
    assert env["RJ_PPG12_PPSIM_REPLAY_SEEDS"] == (
        "2991264730,4256268992,2394322166,874466025,2240380304"
    )
    assert env["RJ_PPG12_PPSIM_EXPECT_PEDESTAL_SEQUENCE"] == "534"
    assert env["RJ_DISABLE_JES_CDB_AUDIT"] == "1"
    assert "RJ_REQUIRE_SIM_GLOBAL" in block["absent"]
    assert env["RJ_PPG12_PHOTON_YIELD"] == "1"
    assert "RJ_PPG12_PHOTON_YIELD_CLUSTER_ERES" in block["absent"]
    assert env["RJ_PPG12_TABLE_QA"] == "1"
    assert env["RJ_PPG12_TABLE_QA_NPB_DATA_TAGGING"] == "1"
    assert env["RJ_PPG12_FIG13_PARITY_QA"] == "1"
    assert env["RJ_PPG12_FIG11_SB_DIAGNOSTIC"] == "1"
    assert env["RJ_PP_PHOTONID_TRAINING_TREE"] == "1"
    assert env["RJ_PP_PHOTONID_TRAINING_TREE_MAX_ENTRIES"] == "0"
    assert env["RJ_PHOTON_ID_ROW_MATCH"] == "newPPG12|newPPG12|newPPG12"
    assert env["RJ_AUTO_MERGE"] == "0"
    assert env["RJ_INTERNAL_JET_PT_MINS"] == "5.0,7.0,10.0,12.0"
    assert env["RJ_INTERNAL_DPHI_PI_FRACTIONS"] == "0.5,0.875"
    assert env["RJ_REQUEST_MEMORY_MB"] == "6000"
    assert env["RJ_DIRECT_NEVENTS"] == "0"
    assert env["RJ_VALIDATE_SIM_INPUT_MAX_LINES"] == "5"
    assert env["RJ_DISABLE_ID_FANOUT"] == "1"
    assert env["RJ_CODEX_CHAT_NAME"] == "THE-97 | unit test"
    assert env["RJ_CODEX_THREAD_ID"] == "unit-thread"
    assert env["RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4"] == "1"
    assert env["RJ_PPG12_PPSIM_G4_ONLY"] == "1"
    assert env["RJ_SIM_ALLOW_NONE_LISTS"] == "1"
    assert env["RJ_FAIL_ON_MISSING_CALO_INPUT"] == "0"
    if interaction == "di":
        assert env["RJ_PPG12_PHOTON_YIELD_DOUBLE"] == "1"
        assert env["RJ_PPG12_PERIOD_STRICT_DI"] == "1"
        assert env["RJ_SIM_SIGNAL_SAMPLE_SET"] == "ppg12_double"
    else:
        assert env["RJ_PPG12_PHOTON_YIELD_DOUBLE"] == "0"
        for key in (
            "RJ_PPG12_PERIOD_STRICT_DI",
            "RJ_SIM_SIGNAL_SAMPLE_SET",
        ):
            assert key in block["absent"]

assert len(seen) == 12
PY
pass "authorized submit calls canonical submitter once per exact lane"

# A reused output namespace must stop before the first submitter call.
rm -f "$call_log"
mkdir -p "$output_root"
if env "${common_env[@]}" MOCK_CALL_LOG="$call_log" \
  RJ_CODEX_CHAT_NAME=test RJ_CODEX_THREAD_ID=test \
  "$driver" --submit --token "SUBMIT_${campaign}" > "${tmp}/reuse.out" 2>&1; then
  fail "submit accepted an existing output namespace"
fi
[[ ! -e "$call_log" ]] || fail "existing-output rejection reached the submitter"
pass "existing output namespace is rejected before mutation"

printf 'ALL TESTS PASSED\n'
