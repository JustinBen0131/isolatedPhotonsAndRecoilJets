#!/usr/bin/env bash
# Plan/run front end for one physical same-event PPG12/RecoilJets oracle lane.
# Default mode is read-only plan.  Execution requires the printed contract
# token and delegates to a foreground-only worker; it never submits Condor.

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
repo_root="$(cd "${script_dir}/../../../.." && pwd -P)"
worker="${script_dir}/run_ppg12_recoiljets_paired_oracle_worker.sh"

mode=plan
provided_token=""
output_dir=""
lane_id=""
sample=""
period=""
interaction=""
setup_script="/opt/sphenix/core/bin/sphenix_setup.sh"
ppg_macro=""
g4_full_list=""
truthjet_full_list=""
apply_bdt="/sphenix/user/shuhangli/ppg12/FunWithxgboost/apply_BDT.C"
apply_config=""
base_e_model="/sphenix/user/shuhangli/ppg12/FunWithxgboost/binned_models/model_base_E_split_single_tmva.root"
base_v3e_model="/sphenix/user/shuhangli/ppg12/FunWithxgboost/binned_models/model_base_v3E_split_single_tmva.root"
npb_model="/sphenix/user/shuhangli/ppg12/FunWithxgboost/npb_models/npb_score_split_tmva.root"
tower_mask="/sphenix/user/shuhangli/ppg12/efficiencytool/tower_masks_bdt_nom.root"
recoil_runtime_manifest=""
recoil_config=""
reuse_ppg_raw_root=""
reuse_ppg_raw_contract=""

usage() {
  cat <<'EOF'
Usage:
  run_ppg12_recoiljets_paired_oracle.sh [--run --token TOKEN] \
    --lane-id photon:photonN:PERIOD:MODE --sample PhotonN \
    --period {0mrad|1p5mrad} --interaction {SI|DI} \
    --output-dir ABS --ppg-macro ABS --g4-full-list ABS --truthjet-full-list ABS \
    --apply-config ABS --recoil-runtime-manifest ABS --recoil-config ABS \
    [--reuse-ppg-raw-root ABS --reuse-ppg-raw-contract ABS] \
    [path overrides]

Default mode is a non-mutating plan.  The supported physical lanes are
Photon5/10/20 x 0mrad/1p5mrad x SI/DI, first five source rows.  Both
executables replay the exact five
PHRandomSeed values captured in historical PPG12 OutDir0; this reconstruction
contract is independent of downstream estimator-toy seed 42.  --run executes
in the foreground only and requires the exact token emitted by plan mode.
EOF
}

die() {
  printf 'PPG12_PAIRED_ORACLE_DRIVER_FAIL: %s\n' "$*" >&2
  exit 2
}

while (($#)); do
  case "$1" in
    --run) mode=run; shift ;;
    --token) [[ $# -ge 2 ]] || die "--token requires a value"; provided_token="$2"; shift 2 ;;
    --lane-id) lane_id="$2"; shift 2 ;;
    --sample) sample="$2"; shift 2 ;;
    --period) period="$2"; shift 2 ;;
    --interaction) interaction="$2"; shift 2 ;;
    --output-dir) output_dir="$2"; shift 2 ;;
    --setup-script) setup_script="$2"; shift 2 ;;
    --ppg-macro) ppg_macro="$2"; shift 2 ;;
    --g4-full-list) g4_full_list="$2"; shift 2 ;;
    --truthjet-full-list) truthjet_full_list="$2"; shift 2 ;;
    --apply-bdt) apply_bdt="$2"; shift 2 ;;
    --apply-config) apply_config="$2"; shift 2 ;;
    --base-e-model) base_e_model="$2"; shift 2 ;;
    --base-v3e-model) base_v3e_model="$2"; shift 2 ;;
    --npb-model) npb_model="$2"; shift 2 ;;
    --tower-mask) tower_mask="$2"; shift 2 ;;
    --recoil-runtime-manifest) recoil_runtime_manifest="$2"; shift 2 ;;
    --recoil-config) recoil_config="$2"; shift 2 ;;
    --reuse-ppg-raw-root) reuse_ppg_raw_root="$2"; shift 2 ;;
    --reuse-ppg-raw-contract) reuse_ppg_raw_contract="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

case "$sample" in Photon5|Photon10|Photon20) ;; *) die "unsupported --sample: $sample" ;; esac
case "$period" in 0mrad|1p5mrad) ;; *) die "unsupported --period: $period" ;; esac
case "$interaction" in SI|DI) ;; *) die "unsupported --interaction: $interaction" ;; esac
sample_lower="$(printf '%s' "$sample" | tr '[:upper:]' '[:lower:]')"
interaction_lower="$(printf '%s' "$interaction" | tr '[:upper:]' '[:lower:]')"
expected_lane_id="photon:${sample_lower}:${period}:${interaction_lower}"
[[ "$lane_id" == "$expected_lane_id" ]] || \
  die "--lane-id differs from physical lane; expected ${expected_lane_id}"

auditor="${repo_root}/scripts/diagnostics/pp_currentian/audit_ppg12_recoiljets_paired_oracle.py"
[[ -f "$auditor" && -s "$auditor" ]] || die "RNG-contract auditor is missing: $auditor"
read -r first_ph_seed pedestal_seed pedestal_sequence ph_seed_sequence < <(
  python3 "$auditor" rng-contract
)
[[ -n "$first_ph_seed" && -n "$pedestal_seed" && -n "$pedestal_sequence" && \
   -n "$ph_seed_sequence" ]] || die "failed to derive reconstruction RNG contract"

[[ -x "$worker" ]] || die "worker is missing or not executable: $worker"

ppg_wrapper="${repo_root}/macros/diagnostics/pp_currentian/Fun4All_ppg12_fixed_seed_oracle.C"
recoil_wrapper="${repo_root}/macros/diagnostics/pp_currentian/Fun4All_recoiljets_fixed_seed_oracle.C"
comparator="${repo_root}/scripts/diagnostics/pp_currentian/compare_ppg12_recoiljets_same_cluster_features.py"
aggregate_extractor="${repo_root}/scripts/diagnostics/pp_currentian/extract_ppg12_recoeff_executable_aggregate.py"
driver_self="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)/$(basename "${BASH_SOURCE[0]}")"

if [[ -n "$reuse_ppg_raw_root" || -n "$reuse_ppg_raw_contract" ]]; then
  [[ -n "$reuse_ppg_raw_root" && -n "$reuse_ppg_raw_contract" ]] || \
    die "raw PPG12 reuse requires both --reuse-ppg-raw-root and --reuse-ppg-raw-contract"
  reuse_mode="exact_contract_bound"
else
  reuse_mode="disabled"
fi

required_keys=(
  output_dir setup_script ppg_macro g4_full_list truthjet_full_list
  apply_bdt apply_config base_e_model base_v3e_model npb_model tower_mask
  recoil_runtime_manifest recoil_config
)
required_values=(
  "$output_dir" "$setup_script" "$ppg_macro" "$g4_full_list"
  "$truthjet_full_list" "$apply_bdt" "$apply_config" "$base_e_model"
  "$base_v3e_model" "$npb_model" "$tower_mask" "$recoil_runtime_manifest"
  "$recoil_config"
)
for ((index = 0; index < ${#required_keys[@]}; ++index)); do
  key="${required_keys[$index]}"
  value="${required_values[$index]}"
  [[ -n "$value" ]] || die "missing required --${key//_/-}"
  [[ "$value" == /* ]] || die "--${key//_/-} must be absolute: $value"
  [[ "$value" != *$'\n'* && "$value" != *$'\r'* ]] || die "unsafe newline in $key"
done

# The caller-facing PPG12 paths are assertions about the intended canonical
# assets; execution itself is source-locked to the copies sealed by the
# runtime builder.  Resolve those copies before constructing the plan token so
# the authorized paths, worker arguments, and contract receipt are identical.
# The byte comparisons retain the caller assertion and fail closed if either
# provenance chain drifts.
manifest_role_path() {
  python3 - "$recoil_runtime_manifest" "$1" <<'PY'
from pathlib import Path
import json
import sys

manifest = Path(sys.argv[1])
role = sys.argv[2]
try:
    data = json.loads(manifest.read_text())
except (OSError, json.JSONDecodeError) as exc:
    raise SystemExit(f"cannot read runtime manifest {manifest}: {exc}")
matches = [
    str(item.get("path", ""))
    for item in data.get("files", [])
    if item.get("role") == role
]
if len(matches) != 1:
    raise SystemExit(f"runtime manifest role {role} has {len(matches)} matches")
path = Path(matches[0])
if not path.is_absolute() or not path.is_file() or path.stat().st_size <= 0:
    raise SystemExit(f"runtime manifest role {role} is not an absolute nonempty file: {path}")
print(path)
PY
}

asserted_apply_bdt="$apply_bdt"
asserted_apply_config="$apply_config"
asserted_base_e_model="$base_e_model"
asserted_base_v3e_model="$base_v3e_model"
asserted_npb_model="$npb_model"
sealed_apply_bdt="$(manifest_role_path ppg_apply_bdt_macro)" || \
  die "failed to resolve sealed apply_BDT macro"
sealed_apply_config="$(manifest_role_path ppg_apply_bdt_config)" || \
  die "failed to resolve sealed apply_BDT config"
sealed_base_e_model="$(manifest_role_path ppg_apply_model_base_E)" || \
  die "failed to resolve sealed base_E model"
sealed_base_v3e_model="$(manifest_role_path ppg_apply_model_base_v3E)" || \
  die "failed to resolve sealed base_v3E model"
sealed_npb_model="$(manifest_role_path ppg_apply_npb_model)" || \
  die "failed to resolve sealed NPB model"
recoeff_period_role="ppg_recoeff_period_config_${period}"
recoeff_period_config="$(manifest_role_path "$recoeff_period_role")" || \
  die "failed to resolve sealed ${period} RecoEff config"
recoeff_truth_vertex_role="ppg_recoeff_truth_vertex_reweight_${period}"
recoeff_truth_vertex_reweight="$(manifest_role_path "$recoeff_truth_vertex_role")" || \
  die "failed to resolve sealed ${period} truth-vertex reweight ROOT"
recoeff_yaml_cpp_header_receipt="$(manifest_role_path ppg_recoeff_yaml_cpp_header_tree_receipt)" || \
  die "failed to resolve sealed yaml-cpp header-tree receipt"

asserted_asset_roles=(apply_bdt apply_config base_e_model base_v3e_model npb_model)
asserted_asset_paths=(
  "$asserted_apply_bdt" "$asserted_apply_config" "$asserted_base_e_model"
  "$asserted_base_v3e_model" "$asserted_npb_model"
)
sealed_asset_paths=(
  "$sealed_apply_bdt" "$sealed_apply_config" "$sealed_base_e_model"
  "$sealed_base_v3e_model" "$sealed_npb_model"
)
for ((index = 0; index < ${#asserted_asset_roles[@]}; ++index)); do
  role="${asserted_asset_roles[$index]}"
  asserted="${asserted_asset_paths[$index]}"
  sealed="${sealed_asset_paths[$index]}"
  [[ -f "$asserted" && -s "$asserted" ]] || \
    die "caller assertion for $role is missing or empty: $asserted"
  cmp -s "$asserted" "$sealed" || \
    die "caller assertion for $role differs from sealed source-locked runtime"
done

apply_bdt="$sealed_apply_bdt"
apply_config="$sealed_apply_config"
base_e_model="$sealed_base_e_model"
base_v3e_model="$sealed_base_v3e_model"
npb_model="$sealed_npb_model"
required_keys+=(
  ppg_recoeff_period_config ppg_recoeff_truth_vertex_reweight
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
  "$driver_self" "$worker" "$ppg_wrapper" "$recoil_wrapper" "$comparator"
  "$auditor" "$aggregate_extractor" "$recoeff_period_config"
  "$recoeff_truth_vertex_reweight" "$recoeff_yaml_cpp_header_receipt"
)
if [[ "$reuse_mode" == exact_contract_bound ]]; then
  token_file_keys+=(reuse_ppg_raw_root reuse_ppg_raw_contract)
  token_file_values+=("$reuse_ppg_raw_root" "$reuse_ppg_raw_contract")
fi
for ((index = 0; index < ${#token_file_keys[@]}; ++index)); do
  key="${token_file_keys[$index]}"
  value="${token_file_values[$index]}"
  [[ "$value" == /* && -f "$value" && -s "$value" ]] || \
    die "token-bound file is missing or empty for $key: $value"
done

sha256_file() {
  python3 - "$1" <<'PY'
from pathlib import Path
import hashlib
import sys
print(hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest())
PY
}

contract_token="$({
  printf '%s\n' \
    'schema_version=3' \
    "lane=${sample}:${period}:${interaction}" \
    "lane_id=${lane_id}" \
    'rows=5' \
    'rng_mode=historical_fifo_replay_v2' \
    "first_ph_seed=${first_ph_seed}" \
    "pedestal_seed=${pedestal_seed}" \
    "pedestal=${pedestal_sequence}" \
    "ph_seed_sequence=${ph_seed_sequence}" \
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

cat <<EOF
PPG12_PAIRED_ORACLE_PLAN
  lane: ${sample} | ${period} | ${interaction}
  lane_id: ${lane_id}
  rows: first five run-28 G4Hits + DST_TRUTH_JET identities
  source_graph: NONE,g4,truthjet,NONE,NONE
  runtime: new.17 for both executables
  RNG: historical FIFO replay; five PH seeds=${ph_seed_sequence}; pedestal=${pedestal_sequence}
  PPG12 + RecoilJets: one source-locked isolated new.17 runtime manifest required
  PPG12 estimator config: ${recoeff_period_role} (${recoeff_period_config})
  PPG12 truth-vertex weights: ${recoeff_truth_vertex_role} (${recoeff_truth_vertex_reweight})
  yaml-cpp headers: sealed tree receipt (${recoeff_yaml_cpp_header_receipt})
  raw PPG12 reuse: ${reuse_mode}${reuse_ppg_raw_root:+ (${reuse_ppg_raw_root})}
  execution: foreground only; no Condor; no merge; no current-pointer mutation
  output_dir: ${output_dir}
  run_token: ${contract_token}
EOF

if [[ "$mode" == plan ]]; then
  [[ -z "$provided_token" ]] || die "--token is valid only with --run"
  exit 0
fi

[[ "$provided_token" == "$contract_token" ]] || die "run token does not match this exact contract"
[[ ! -e "$output_dir" ]] || die "output path already exists: $output_dir"
[[ -z "${RANDOMSEED+x}" ]] || die "refusing inherited shell RANDOMSEED; historical replay requires it absent"

export RJ_PPG12_PAIRED_ORACLE_RUN_TOKEN="$contract_token"
exec "$worker" \
  "$contract_token" "$repo_root" "$output_dir" "$setup_script" \
  "$lane_id" "$sample" "$period" "$interaction" \
  "$ppg_macro" "$g4_full_list" "$truthjet_full_list" \
  "$apply_bdt" "$apply_config" "$base_e_model" "$base_v3e_model" \
  "$npb_model" "$tower_mask" "$recoil_runtime_manifest" "$recoil_config" \
  "${reuse_ppg_raw_root:-NONE}" "${reuse_ppg_raw_contract:-NONE}"
