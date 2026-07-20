#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd -P)"
driver="${repo_root}/scripts/sdcc/workflows/diagnostics/run_ppg12_recoiljets_paired_oracle.sh"
worker="${repo_root}/scripts/sdcc/workflows/diagnostics/run_ppg12_recoiljets_paired_oracle_worker.sh"
auditor="${repo_root}/scripts/diagnostics/pp_currentian/audit_ppg12_recoiljets_paired_oracle.py"

bash -n "$driver"
bash -n "$worker"
python3 -m py_compile "$auditor"

for invariant in \
  'calo_calib="${calo_macro_dir}/Calo_Calib.C"' \
  'gROOT->SetMacroPath((std::string(' \
  'export ROOT_INCLUDE_PATH="${calo_macro_dir}:${base_root_include_path}"' \
  'export ROOT_INCLUDE_PATH="${calo_macro_dir}:${photon_builder_include_root}:${base_root_include_path}"'; do
  grep -Fq "$invariant" "$worker" || {
    echo "missing exact Calo_Calib routing invariant: $invariant" >&2
    exit 1
  }
done
for invariant in \
  'ppg_raw_root="${ppg_dir}/caloana.root"' \
  'ppg_scored_root="${ppg_dir}/caloana_with_bdt_split.root"' \
  'run_recoeff baseline' \
  'run_recoeff trace' \
  '--root-equivalence-only' \
  '--ppg12-executable-trace "$ppg_candidate_trace"' \
  '--ppg12-executable-response-trace "$ppg_response_trace"' \
  'export RJ_PP_PHOTONID_PPG12_FILTER=0' \
  'token-bound path or file content changed after plan authorization' \
  'reuse_mode="exact_contract_bound"' \
  'reuse contract predates content-bound schema 3 and is scientifically inadmissible' \
  'reuse contract authorization token does not replay exactly' \
  'reuse source hash differs from the prior apply_BDT input hash' \
  'FAILED-run reuse is limited to the known attempt6 raw ROOT' \
  'known FAILED-run reuse lacks the 5000-event completion log' \
  'FAILED-run reuse is not the known missing-yaml-header failure' \
  'input.TestBit(TFile::kRecovered)' \
  'ppg_recoeff_truth_vertex_reweight' \
  'ppg_recoeff_yaml_cpp_header_tree_receipt' \
  'export ROOT_INCLUDE_PATH="${recoeff_yaml_cpp_include_root}:${base_root_include_path}"' \
  '--out-json "$aggregate_report"'; do
  grep -Fq -- "$invariant" "$worker" || {
    echo "missing executable-estimator invariant: $invariant" >&2
    exit 1
  }
done
if grep -Fq 'ppg_raw_root="${ppg_dir}/output_sim.root"' "$worker"; then
  echo "worker still assumes the non-executable output_sim.root name" >&2
  exit 1
fi
if grep -Fq 'gSystem->Which(gROOT->GetMacroPath(), "Calo_Calib.C")' "$worker"; then
  echo "worker still relies on dynamic Calo_Calib discovery" >&2
  exit 1
fi

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
out="${tmp}/oracle-output"

for asset in \
  new17-setup.sh Fun4All_run_sim.C g4.list truth.list apply_BDT.C \
  config_nom.yaml base_E.root base_v3E.root npb.root mask.root \
  recoil_runtime.json recoil_config.yaml; do
  printf 'token-bound test asset: %s\n' "$asset" > "${tmp}/${asset}"
done

# The caller-facing estimator paths need not be the executable paths.  The
# source-locked runtime manifest is authoritative, so exercise the real SDCC
# topology: five distinct sealed files with byte-identical caller aliases.
sealed_dir="${tmp}/sealed-runtime/estimator"
mkdir -p "$sealed_dir"
cp "${tmp}/apply_BDT.C" "${sealed_dir}/apply_BDT.C"
cp "${tmp}/config_nom.yaml" "${sealed_dir}/config_nom.yaml"
cp "${tmp}/base_E.root" "${sealed_dir}/model_base_E.root"
cp "${tmp}/base_v3E.root" "${sealed_dir}/model_base_v3E.root"
cp "${tmp}/npb.root" "${sealed_dir}/npb_score.root"
printf 'period estimator config: 0mrad\n' > "${sealed_dir}/config_bdt_nom_0rad.yaml"
printf 'period estimator config: 1p5mrad\n' > "${sealed_dir}/config_bdt_nom_1p5mrad.yaml"
printf 'sealed truth-vertex weights: 0mrad\n' > "${sealed_dir}/truth_vertex_reweight_0mrad.root"
printf 'sealed truth-vertex weights: 1p5mrad\n' > "${sealed_dir}/truth_vertex_reweight_1p5mrad.root"
printf '{"schema_version":1,"role":"ppg_recoeff_yaml_cpp_header_tree"}\n' \
  > "${sealed_dir}/yaml_cpp_header_tree_receipt.json"
python3 - "${tmp}/recoil_runtime.json" "$sealed_dir" <<'PY'
from pathlib import Path
import json
import sys

manifest = Path(sys.argv[1])
sealed = Path(sys.argv[2])
roles = {
    "ppg_apply_bdt_macro": sealed / "apply_BDT.C",
    "ppg_apply_bdt_config": sealed / "config_nom.yaml",
    "ppg_apply_model_base_E": sealed / "model_base_E.root",
    "ppg_apply_model_base_v3E": sealed / "model_base_v3E.root",
    "ppg_apply_npb_model": sealed / "npb_score.root",
    "ppg_recoeff_period_config_0mrad": sealed / "config_bdt_nom_0rad.yaml",
    "ppg_recoeff_period_config_1p5mrad": sealed / "config_bdt_nom_1p5mrad.yaml",
    "ppg_recoeff_truth_vertex_reweight_0mrad": sealed / "truth_vertex_reweight_0mrad.root",
    "ppg_recoeff_truth_vertex_reweight_1p5mrad": sealed / "truth_vertex_reweight_1p5mrad.root",
    "ppg_recoeff_yaml_cpp_header_tree_receipt": sealed / "yaml_cpp_header_tree_receipt.json",
}
manifest.write_text(
    json.dumps(
        {"files": [{"role": role, "path": str(path)} for role, path in roles.items()]},
        indent=2,
        sort_keys=True,
    )
    + "\n"
)
PY

common_args=(
  --lane-id photon:photon5:1p5mrad:si
  --sample Photon5
  --period 1p5mrad
  --interaction SI
  --output-dir "$out"
  --setup-script "${tmp}/new17-setup.sh"
  --ppg-macro "${tmp}/Fun4All_run_sim.C"
  --g4-full-list "${tmp}/g4.list"
  --truthjet-full-list "${tmp}/truth.list"
  --apply-bdt "${tmp}/apply_BDT.C"
  --apply-config "${tmp}/config_nom.yaml"
  --base-e-model "${tmp}/base_E.root"
  --base-v3e-model "${tmp}/base_v3E.root"
  --npb-model "${tmp}/npb.root"
  --tower-mask "${tmp}/mask.root"
  --recoil-runtime-manifest "${tmp}/recoil_runtime.json"
  --recoil-config "${tmp}/recoil_config.yaml"
)

plan="$($driver "${common_args[@]}")"
grep -Fq 'lane: Photon5 | 1p5mrad | SI' <<<"$plan"
grep -Fq 'lane_id: photon:photon5:1p5mrad:si' <<<"$plan"
grep -Fq 'source_graph: NONE,g4,truthjet,NONE,NONE' <<<"$plan"
grep -Fq 'RNG: historical FIFO replay; five PH seeds=2991264730,4256268992,2394322166,874466025,2240380304; pedestal=534' <<<"$plan"
grep -Fq 'PPG12 + RecoilJets: one source-locked isolated new.17 runtime manifest required' <<<"$plan"
grep -Fq "PPG12 estimator config: ppg_recoeff_period_config_1p5mrad (${sealed_dir}/config_bdt_nom_1p5mrad.yaml)" <<<"$plan"
grep -Fq "PPG12 truth-vertex weights: ppg_recoeff_truth_vertex_reweight_1p5mrad (${sealed_dir}/truth_vertex_reweight_1p5mrad.root)" <<<"$plan"
grep -Fq 'execution: foreground only; no Condor; no merge' <<<"$plan"
token="$(sed -n 's/^  run_token: //p' <<<"$plan")"
[[ "$token" =~ ^ppg12-oracle:[0-9a-f]{64}$ ]]

# Passing the sealed paths directly must reproduce the caller-alias plan token:
# authorization binds the executable runtime assets, not a byte-identical
# lexical alias that the worker later replaces.
sealed_args=("${common_args[@]}")
for ((i = 0; i < ${#sealed_args[@]}; ++i)); do
  case "${sealed_args[$i]}" in
    "${tmp}/apply_BDT.C") sealed_args[$i]="${sealed_dir}/apply_BDT.C" ;;
    "${tmp}/config_nom.yaml") sealed_args[$i]="${sealed_dir}/config_nom.yaml" ;;
    "${tmp}/base_E.root") sealed_args[$i]="${sealed_dir}/model_base_E.root" ;;
    "${tmp}/base_v3E.root") sealed_args[$i]="${sealed_dir}/model_base_v3E.root" ;;
    "${tmp}/npb.root") sealed_args[$i]="${sealed_dir}/npb_score.root" ;;
  esac
done
sealed_plan="$($driver "${sealed_args[@]}")"
sealed_token="$(sed -n 's/^  run_token: //p' <<<"$sealed_plan")"
[[ "$sealed_token" == "$token" ]]

# Replay the emitted token from a worker-style contract containing only the
# sealed executable paths.  This is the exact regression for the first real
# lane, where the old driver authorized caller aliases and the worker recorded
# their manifest-selected replacements.
python3 - "$auditor" "$tmp" "$repo_root" "$token" "$out" "$driver" "$worker" <<'PY'
from pathlib import Path
import hashlib
import importlib.util
import sys

auditor_path = Path(sys.argv[1])
tmp = Path(sys.argv[2])
repo = Path(sys.argv[3])
token = sys.argv[4]
contract_path = Path(sys.argv[5]) / "paired_oracle_contract.json"
driver = Path(sys.argv[6]).resolve()
worker = Path(sys.argv[7]).resolve()

spec = importlib.util.spec_from_file_location("paired_oracle_token_replay", auditor_path)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)

sealed = tmp / "sealed-runtime" / "estimator"
paths = {
    "setup_script": str(tmp / "new17-setup.sh"),
    "ppg_macro": str(tmp / "Fun4All_run_sim.C"),
    "g4_full_list": str(tmp / "g4.list"),
    "truthjet_full_list": str(tmp / "truth.list"),
    "apply_bdt": str(sealed / "apply_BDT.C"),
    "apply_config": str(sealed / "config_nom.yaml"),
    "base_e_model": str(sealed / "model_base_E.root"),
    "base_v3e_model": str(sealed / "model_base_v3E.root"),
    "npb_model": str(sealed / "npb_score.root"),
    "tower_mask": str(tmp / "mask.root"),
    "recoil_runtime_manifest": str(tmp / "recoil_runtime.json"),
    "recoil_config": str(tmp / "recoil_config.yaml"),
    "ppg_recoeff_period_config": str(sealed / "config_bdt_nom_1p5mrad.yaml"),
    "ppg_recoeff_truth_vertex_reweight": str(sealed / "truth_vertex_reweight_1p5mrad.root"),
    "ppg_recoeff_yaml_cpp_header_tree_receipt": str(sealed / "yaml_cpp_header_tree_receipt.json"),
}
token_paths = {
    **paths,
    "driver_script": str(driver),
    "worker_script": str(worker),
    "ppg_wrapper": str(repo / "macros/diagnostics/pp_currentian/Fun4All_ppg12_fixed_seed_oracle.C"),
    "recoil_wrapper": str(repo / "macros/diagnostics/pp_currentian/Fun4All_recoiljets_fixed_seed_oracle.C"),
    "comparator": str(repo / "scripts/diagnostics/pp_currentian/compare_ppg12_recoiljets_same_cluster_features.py"),
    "auditor": str(auditor_path.resolve()),
    "aggregate_extractor": str(repo / "scripts/diagnostics/pp_currentian/extract_ppg12_recoeff_executable_aggregate.py"),
}

def digest(path: str) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

contract = {
    "authorization_token": token,
    "authorization": {
        "token_bound_files": {
            role: {"path": path, "sha256": digest(path)}
            for role, path in token_paths.items()
        },
        "raw_ppg12_reuse": {"mode": "disabled", "source": None, "contract": None},
    },
    "lane": {
        "lane_id": "photon:photon5:1p5mrad:si",
        "sample": "Photon5",
        "period": "1p5mrad",
        "interaction": "SI",
        "rows": 5,
    },
    "rng": module.historical_rng_contract(),
    "paths": paths,
}
module.validate_authorization_token(
    contract_path,
    contract,
    lane_id="photon:photon5:1p5mrad:si",
)
PY

printf 'caller alias is no longer byte-identical\n' >> "${tmp}/apply_BDT.C"
if "$driver" "${common_args[@]}" \
    >"${tmp}/caller-alias-mismatch.out" 2>"${tmp}/caller-alias-mismatch.err"; then
  echo "non-identical caller estimator alias unexpectedly passed" >&2
  exit 1
fi
grep -Fq 'differs from sealed' "${tmp}/caller-alias-mismatch.err"
cp "${sealed_dir}/apply_BDT.C" "${tmp}/apply_BDT.C"

di_args=("${common_args[@]}")
for ((i = 0; i < ${#di_args[@]}; ++i)); do
  case "${di_args[$i]}" in
    photon:photon5:1p5mrad:si) di_args[$i]=photon:photon20:0mrad:di ;;
    Photon5) di_args[$i]=Photon20 ;;
    1p5mrad) di_args[$i]=0mrad ;;
    SI) di_args[$i]=DI ;;
  esac
done
di_plan="$($driver "${di_args[@]}")"
grep -Fq 'lane: Photon20 | 0mrad | DI' <<<"$di_plan"
grep -Fq 'lane_id: photon:photon20:0mrad:di' <<<"$di_plan"
grep -Fq "PPG12 estimator config: ppg_recoeff_period_config_0mrad (${sealed_dir}/config_bdt_nom_0rad.yaml)" <<<"$di_plan"
grep -Fq "PPG12 truth-vertex weights: ppg_recoeff_truth_vertex_reweight_0mrad (${sealed_dir}/truth_vertex_reweight_0mrad.root)" <<<"$di_plan"
di_token="$(sed -n 's/^  run_token: //p' <<<"$di_plan")"
[[ "$di_token" =~ ^ppg12-oracle:[0-9a-f]{64}$ ]]
[[ "$di_token" != "$token" ]]

bad_lane_args=("${common_args[@]}")
bad_lane_args[1]=photon:photon10:1p5mrad:si
if "$driver" "${bad_lane_args[@]}" >"${tmp}/bad-lane.out" 2>"${tmp}/bad-lane.err"; then
  echo "mismatched physical lane identity unexpectedly passed" >&2
  exit 1
fi
grep -Fq 'lane-id differs from physical lane' "${tmp}/bad-lane.err"

printf 'mutated after plan\n' >> "${tmp}/recoil_config.yaml"
mutated_plan="$($driver "${common_args[@]}")"
mutated_token="$(sed -n 's/^  run_token: //p' <<<"$mutated_plan")"
[[ "$mutated_token" =~ ^ppg12-oracle:[0-9a-f]{64}$ ]]
[[ "$mutated_token" != "$token" ]]
if "$driver" --run --token "$token" "${common_args[@]}" \
    >"${tmp}/stale-content.out" 2>"${tmp}/stale-content.err"; then
  echo "content-stale token unexpectedly passed" >&2
  exit 1
fi
grep -Fq 'run token does not match' "${tmp}/stale-content.err"
printf 'token-bound test asset: recoil_config.yaml\n' > "${tmp}/recoil_config.yaml"

printf 'mutated after plan\n' >> "${sealed_dir}/truth_vertex_reweight_1p5mrad.root"
mutated_vertex_plan="$($driver "${common_args[@]}")"
mutated_vertex_token="$(sed -n 's/^  run_token: //p' <<<"$mutated_vertex_plan")"
[[ "$mutated_vertex_token" =~ ^ppg12-oracle:[0-9a-f]{64}$ ]]
[[ "$mutated_vertex_token" != "$token" ]]
if "$driver" --run --token "$token" "${common_args[@]}" \
    >"${tmp}/stale-vertex.out" 2>"${tmp}/stale-vertex.err"; then
  echo "truth-vertex-weight-stale token unexpectedly passed" >&2
  exit 1
fi
grep -Fq 'run token does not match' "${tmp}/stale-vertex.err"
printf 'sealed truth-vertex weights: 1p5mrad\n' \
  > "${sealed_dir}/truth_vertex_reweight_1p5mrad.root"

reuse_raw="${tmp}/caloana.root"
reuse_contract="${tmp}/prior_contract.json"
printf 'raw reuse bytes\n' > "$reuse_raw"
printf '{"schema_version":2}\n' > "$reuse_contract"
reuse_plan="$($driver "${common_args[@]}" \
  --reuse-ppg-raw-root "$reuse_raw" \
  --reuse-ppg-raw-contract "$reuse_contract")"
reuse_token="$(sed -n 's/^  run_token: //p' <<<"$reuse_plan")"
[[ "$reuse_token" =~ ^ppg12-oracle:[0-9a-f]{64}$ ]]
[[ "$reuse_token" != "$token" ]]
grep -Fq 'raw PPG12 reuse: exact_contract_bound' <<<"$reuse_plan"
if "$driver" "${common_args[@]}" --reuse-ppg-raw-root "$reuse_raw" \
    >"${tmp}/partial-reuse.out" 2>"${tmp}/partial-reuse.err"; then
  echo "partial reuse contract unexpectedly passed" >&2
  exit 1
fi
grep -Fq 'requires both' "${tmp}/partial-reuse.err"

if "$driver" --reco-seed 5709 "${common_args[@]}" \
    >"${tmp}/synthetic-seed.out" 2>"${tmp}/synthetic-seed.err"; then
  echo "synthetic reconstruction seed unexpectedly passed" >&2
  exit 1
fi
grep -Fq 'unknown argument: --reco-seed' "${tmp}/synthetic-seed.err"

if "$driver" --run --token ppg12-oracle:wrong "${common_args[@]}" \
    >"${tmp}/wrong.out" 2>"${tmp}/wrong.err"; then
  echo "wrong token unexpectedly passed" >&2
  exit 1
fi
grep -Fq 'run token does not match' "${tmp}/wrong.err"
[[ ! -e "$out" ]]

if RANDOMSEED=17 "$driver" --run --token "$token" "${common_args[@]}" \
    >"${tmp}/seed.out" 2>"${tmp}/seed.err"; then
  echo "inherited RANDOMSEED unexpectedly passed" >&2
  exit 1
fi
grep -Fq 'refusing inherited shell RANDOMSEED' "${tmp}/seed.err"
[[ ! -e "$out" ]]

if "$driver" --token "$token" "${common_args[@]}" \
    >"${tmp}/plan-token.out" 2>"${tmp}/plan-token.err"; then
  echo "plan mode unexpectedly accepted --token" >&2
  exit 1
fi
grep -Fq -- '--token is valid only with --run' "${tmp}/plan-token.err"

python3 - "$auditor" "$tmp" <<'PY'
from pathlib import Path
from unittest import mock
import importlib.util
import sys

module_path = Path(sys.argv[1])
tmp = Path(sys.argv[2])
spec = importlib.util.spec_from_file_location("paired_oracle_audit", module_path)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)

assert module.historical_rng_contract() == {
    "mode": "historical_fifo_replay_v2",
    "reco_consts_randomseed": "absent",
    "ph_seed_sequence": [
        2991264730, 4256268992, 2394322166, 874466025, 2240380304,
    ],
    "first_ph_seed": 2991264730,
    "pedestal_seed": 4256268992,
    "pedestal_sequence": 534,
    "pedestal_file": "pedestal-54256-00534.root",
    "source_log": module.HISTORICAL_RNG_SOURCE,
    "source_log_call_count": 5,
}

g4 = tmp / "g4.list"
truth = tmp / "truth.list"
combined = tmp / "combined.list"
g4_rows = [
    f"/source/G4Hits_pythia8_PhotonJet5-0000000028-{i:06d}.root"
    for i in range(5)
]
truth_rows = [
    f"/source/DST_TRUTH_JET_pythia8_PhotonJet5-0000000028-{i:06d}.root"
    for i in range(5)
]
g4.write_text("\n".join(g4_rows) + "\n")
truth.write_text("\n".join(truth_rows) + "\n")
combined.write_text(
    "\n".join(
        f"NONE {left} {right} NONE NONE"
        for left, right in zip(g4_rows, truth_rows)
    )
    + "\n"
)
si_evidence = module.validate_source_graph(
    g4, truth, combined, sample="Photon5", interaction="SI"
)
assert si_evidence["sample"] == "Photon5"
assert si_evidence["interaction"] == "SI"

bad = tmp / "four-stream.list"
bad.write_text(
    "\n".join(
        f"/calo {left} {right} /global /mbd"
        for left, right in zip(g4_rows, truth_rows)
    )
    + "\n"
)
try:
    module.validate_source_graph(
        g4, truth, bad, sample="Photon5", interaction="SI"
    )
except module.AuditFailure as exc:
    assert "exact oracle graph" in str(exc)
else:
    raise AssertionError("four-stream source graph unexpectedly passed")

di_g4 = tmp / "di-g4.list"
di_truth = tmp / "di-truth.list"
di_combined = tmp / "di-combined.list"
di_g4_rows = [
    f"/source/G4Hits_pythia8_PhotonJet20_pythia8_Detroit-0000000028-{i:06d}.root"
    for i in range(5)
]
di_truth_rows = [
    f"/source/DST_TRUTH_JET_pythia8_PhotonJet20_pythia8_Detroit-0000000028-{i:06d}.root"
    for i in range(5)
]
di_g4.write_text("\n".join(di_g4_rows) + "\n")
di_truth.write_text("\n".join(di_truth_rows) + "\n")
di_combined.write_text(
    "\n".join(
        f"NONE {left} {right} NONE NONE"
        for left, right in zip(di_g4_rows, di_truth_rows)
    )
    + "\n"
)
di_evidence = module.validate_source_graph(
    di_g4, di_truth, di_combined, sample="Photon20", interaction="DI"
)
assert di_evidence["interaction"] == "DI"
for wrong_mode, left, right, combined_path, sample in (
    ("DI", g4, truth, combined, "Photon5"),
    ("SI", di_g4, di_truth, di_combined, "Photon20"),
):
    try:
        module.validate_source_graph(
            left, right, combined_path, sample=sample, interaction=wrong_mode
        )
    except module.AuditFailure as exc:
        assert "unexpected sample/run/segment identity" in str(exc)
    else:
        raise AssertionError(f"{wrong_mode} accepted the opposite source identity mode")

macro = tmp / "Fun4All_run_sim.C"
macro.write_text(
    "INPUTREADHITS::listfile[0] = inputFile0;\n"
    "// INPUTREADHITS::listfile[1] = inputFile1;\n"
    "// INPUTREADHITS::listfile[2] = inputFile2;\n"
    "// INPUTREADHITS::listfile[3] = inputFile3;\n"
    "INPUTREADHITS::listfile[4] = inputFile4;\n"
)
module.validate_frozen_macro_graph(macro)
macro.write_text(macro.read_text() + "INPUTREADHITS::listfile[1] = inputFile1;\n")
try:
    module.validate_frozen_macro_graph(macro)
except module.AuditFailure as exc:
    assert "index 1" in str(exc)
else:
    raise AssertionError("active four-stream macro unexpectedly passed")

arbitrary = tmp / "arbitrary"
arbitrary.write_text("x")
stale = next(iter(module.STALE_LOCAL_PPG_MACRO_HASHES))
with mock.patch.object(module, "sha256", return_value=stale):
    try:
        module.require_hash(arbitrary, "0" * 64, "preserved macro")
    except module.AuditFailure as exc:
        assert "known stale local PPG macro" in str(exc)
    else:
        raise AssertionError("known stale macro digest unexpectedly passed")

log = tmp / "wrong-seed.log"
log.write_text(
    "ORACLE_RUNTIME side=ppg12 profile=new.17 "
    f"offline_main={module.EXPECTED_OFFLINE_MAIN}\n"
    "ORACLE_SEED_CONTRACT side=ppg12 mode=historical_fifo_replay_v2 "
    "rc_randomseed=absent "
    "ph_seed_sequence=2991264730,4256268992,2394322166,874466025,2240380304 "
    "pedestal_sequence=534\n"
    "PHRandomSeed::GetSeed() seed: 7\n"
    "PHRandomSeed::GetSeed() seed: 3421126067\n"
    "PHRandomSeed::GetSeed() seed: 4083286876\n"
    "PHRandomSeed::GetSeed() seed: 787846414\n"
    "PHRandomSeed::GetSeed() seed: 3143890026\n"
)
try:
    module.validate_log(
        log, "ppg12", arbitrary, module.historical_rng_contract()
    )
except module.AuditFailure as exc:
    assert "PHRandomSeed call sequence" in str(exc)
else:
    raise AssertionError("wrong first PH seed unexpectedly passed")

report = tmp / "candidate_report.csv"
feature_names = (
    "cluster_Et_score_input", "cluster_weta_cogx", "cluster_wphi_cogx",
    "vertexz", "cluster_Eta", "e11_over_e33", "cluster_et1",
    "cluster_et2", "cluster_et3", "cluster_et4", "e32_over_e35",
)
row = {
    "match_status": "matched",
    "candidate_identity": "seg0:evt1:trk7:rj1:ppg1",
    "identity_match": "1",
    "model_route_agree": "1",
    "common_agree": "1",
    "tag_agree": "1",
    "signal_status_agree": "1",
    "abcd_agree": "1",
    "rj_inferred_stored_model": "base_v3E",
    "ppg12_selected_model": "base_v3E",
    "ppg12_tag_evidence_source": module.PRESERVED_EXECUTABLE_EVIDENCE,
    "ppg12_isolation_abcd_evidence_source": module.PRESERVED_EXECUTABLE_EVIDENCE,
    "ppg12_truth_response_fill_evidence_source": module.PRESERVED_EXECUTABLE_EVIDENCE,
    "ppg12_weight_evidence_source": module.PRESERVED_EXECUTABLE_EVIDENCE,
    "base_E_score_delta": "0",
    "base_v3E_score_delta": "0",
    "selected_score_delta": "0",
    "stored_minus_routed_score": "0",
    "rj_stored_bdt_score": "0.8",
    "ppg12_selected_bdt_score": "0.8",
    "stored_common_agree": "1",
    "stored_tag_agree": "1",
    "stored_abcd_agree": "1",
    "rj_is_iso": "1",
    "ppg12_is_iso": "1",
    "rj_is_noniso": "0",
    "ppg12_is_noniso": "0",
    "rj_raw_eiso": "0.2",
    "ppg12_raw_eiso": "0.2",
    "rj_corrected_eiso": "0.3",
    "ppg12_corrected_eiso": "0.3",
    "rj_iso_threshold": "1.0",
    "ppg12_iso_threshold": "1.0",
    "rj_noniso_threshold": "3.0",
    "ppg12_noniso_threshold": "3.0",
    "truth_class_agree": "1",
    "rj_truth_class": "1",
    "truth_class": "1",
    "ppg12_is_signal": "1",
    "rj_logical_abcd_region": "1",
    "ppg12_logical_abcd_region": "1",
    "rj_analysis_window_pass": "1",
    "ppg12_analysis_window_pass": "1",
    "rj_response_Et": "20.0",
    "ppg12_response_Et": "20.0",
    "truth_pt": "20.0",
    "rj_response_window_pass": "1",
    "ppg12_response_window_pass": "1",
    "rj_signal_fill_A": "1",
    "rj_signal_fill_B": "0",
    "rj_signal_fill_C": "0",
    "rj_signal_fill_D": "0",
    "ppg12_signal_fill_A": "1",
    "ppg12_signal_fill_B": "0",
    "ppg12_signal_fill_C": "0",
    "ppg12_signal_fill_D": "0",
    "rj_signal_fill_multiplicity": "1",
    "ppg12_fill_multiplicity": "1",
    "rj_weight_lane_code": "1",
    "rj_weight_component_code": "1",
    "rj_weight_slice": "2.0",
    "ppg12_weight_sample": "2.0",
    "rj_weight_mix": "1.0",
    "ppg12_weight_mix": "1.0",
    "rj_weight_period": "1.0",
    "ppg12_weight_lumi": "1.0",
    "ppg12_weight_cross": "2.0",
    "rj_weight_vertex": "1.0",
    "ppg12_weight_vertex": "1.0",
    "ppg12_weight_truth_vertex": "1.0",
    "ppg12_weight_trigger": "1.0",
    "ppg12_weight_event": "2.0",
    "rj_weight_final": "2.0",
    "ppg12_weight_final": "2.0",
    "rj_weight_product_delta": "0",
    "rj_event_weight_delta": "0",
}
for feature in feature_names:
    row[f"{feature}_rj"] = "1.0"
    row[f"{feature}_ppg12"] = "1.0"
    row[f"{feature}_delta"] = "0"

def write_report(value):
    with report.open("w", newline="") as stream:
        writer = __import__("csv").DictWriter(stream, fieldnames=list(value))
        writer.writeheader()
        writer.writerow(value)

write_report(row)
assert module.validate_candidate_report(report) == 1
first_divergence = __import__("json").loads((tmp / "first_divergence.json").read_text())
assert first_divergence["status"] == "PASS"
assert first_divergence["first_divergence"] is None

for field, bad_value, expected in (
    ("rj_inferred_stored_model", "base_E", "stored model-route divergence"),
    ("stored_minus_routed_score", "0.01", "stored score divergence"),
    ("rj_stored_bdt_score", "0.80001", "stored-to-oracle score divergence"),
    ("stored_common_agree", "0", "stored_common_agree divergence"),
    ("stored_tag_agree", "0", "stored_tag_agree divergence"),
    ("stored_abcd_agree", "0", "stored_abcd_agree divergence"),
):
    mutated = dict(row)
    mutated[field] = bad_value
    write_report(mutated)
    try:
        module.validate_candidate_report(report)
    except module.AuditFailure as exc:
        assert expected in str(exc), (field, str(exc))
    else:
        raise AssertionError(f"stored Recoil evidence mutation passed: {field}")

mutated = dict(row)
mutated["cluster_Et_score_input_delta"] = "0.1"
write_report(mutated)
try:
    module.validate_candidate_report(report)
except module.AuditFailure:
    first_divergence = __import__("json").loads((tmp / "first_divergence.json").read_text())
    assert first_divergence["first_divergence"]["stage"] == "features"
    assert first_divergence["first_divergence"]["field"] == "cluster_Et_score_input"
else:
    raise AssertionError("feature mutation unexpectedly passed")

mutated = dict(row)
mutated["selected_score_delta"] = "nan"
write_report(mutated)
try:
    module.validate_candidate_report(report)
except module.AuditFailure as exc:
    assert "non-finite selected_score_delta" in str(exc)
else:
    raise AssertionError("non-finite candidate score unexpectedly passed")
PY

if rg -n 'condor_submit|condor_q|condor_rm|addChunks|current\.json' "$driver" "$worker"; then
  echo "foreground harness contains a forbidden submit/merge/promotion primitive" >&2
  exit 1
fi

if grep -Fq 'export RJ_PPG12_PHOTON_YIELD_CLUSTER_ERES=0.04' "$worker"; then
  echo "paired oracle still applies forbidden blanket 4 percent classification smear" >&2
  exit 1
fi

printf 'PPG12_PAIRED_ORACLE_TEST_PASS\n'
