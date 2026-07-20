#!/usr/bin/env bash
# Plan/run front end for the one-lane, same-event PPG12/RecoilJets oracle.
# Default mode is read-only plan.  Execution requires the printed contract
# token and delegates to a foreground-only worker; it never submits Condor.

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
repo_root="$(cd "${script_dir}/../../../.." && pwd -P)"
worker="${script_dir}/run_ppg12_recoiljets_paired_oracle_worker.sh"

mode=plan
provided_token=""
output_dir=""
setup_script="/opt/sphenix/core/bin/sphenix_setup.sh"
ppg_macro="/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/photon5/Fun4All_run_sim.C"
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

usage() {
  cat <<'EOF'
Usage:
  run_ppg12_recoiljets_paired_oracle.sh [--run --token TOKEN] \
    --output-dir ABS --g4-full-list ABS --truthjet-full-list ABS \
    --apply-config ABS --recoil-runtime-manifest ABS --recoil-config ABS \
    [path overrides]

Default mode is a non-mutating plan.  The only supported lane is Photon5,
1p5mrad, SI, first five source rows.  Both executables replay the exact five
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
    -h|--help) usage; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

auditor="${repo_root}/scripts/diagnostics/pp_currentian/audit_ppg12_recoiljets_paired_oracle.py"
[[ -f "$auditor" && -s "$auditor" ]] || die "RNG-contract auditor is missing: $auditor"
read -r first_ph_seed pedestal_seed pedestal_sequence ph_seed_sequence < <(
  python3 "$auditor" rng-contract
)
[[ -n "$first_ph_seed" && -n "$pedestal_seed" && -n "$pedestal_sequence" && \
   -n "$ph_seed_sequence" ]] || die "failed to derive reconstruction RNG contract"

[[ -x "$worker" ]] || die "worker is missing or not executable: $worker"

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

contract_token="$({
  printf '%s\n' \
    'schema_version=2' \
    'lane=Photon5:1p5mrad:SI' \
    'rows=5' \
    'rng_mode=historical_fifo_replay_v2' \
    "first_ph_seed=${first_ph_seed}" \
    "pedestal_seed=${pedestal_seed}" \
    "pedestal=${pedestal_sequence}" \
    "ph_seed_sequence=${ph_seed_sequence}" \
    'source_graph=NONE,g4,truthjet,NONE,NONE' \
    'runtime=new.17'
  for ((index = 0; index < ${#required_keys[@]}; ++index)); do
    printf '%s=%s\n' "${required_keys[$index]}" "${required_values[$index]}"
  done
} | LC_ALL=C sort | python3 -c \
  'import hashlib,sys; print("ppg12-oracle:" + hashlib.sha256(sys.stdin.buffer.read()).hexdigest())')"

cat <<EOF
PPG12_PAIRED_ORACLE_PLAN
  lane: Photon5 | 1p5mrad | SI
  rows: first five run-28 G4Hits + DST_TRUTH_JET identities
  source_graph: NONE,g4,truthjet,NONE,NONE
  runtime: new.17 for both executables
  RNG: historical FIFO replay; five PH seeds=${ph_seed_sequence}; pedestal=${pedestal_sequence}
  PPG12 + RecoilJets: one source-locked isolated new.17 runtime manifest required
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
  "$ppg_macro" "$g4_full_list" "$truthjet_full_list" \
  "$apply_bdt" "$apply_config" "$base_e_model" "$base_v3e_model" \
  "$npb_model" "$tower_mask" "$recoil_runtime_manifest" "$recoil_config"
