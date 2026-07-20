#!/usr/bin/env bash
# Foreground-only implementation for run_ppg12_recoiljets_paired_oracle.sh.
# This script is intentionally not a submitter and contains no merge/promotion
# path.  It writes only below a new, explicitly authorized output directory.

set -euo pipefail

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

case "$sample" in Photon5|Photon10|Photon20) ;; *) die "unsupported sample: $sample" ;; esac
case "$period" in 0mrad|1p5mrad) ;; *) die "unsupported period: $period" ;; esac
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
token_file_keys=(
  setup_script ppg_macro g4_full_list truthjet_full_list apply_bdt
  apply_config base_e_model base_v3e_model npb_model tower_mask
  recoil_runtime_manifest recoil_config driver_script worker_script
  ppg_wrapper recoil_wrapper comparator auditor aggregate_extractor
)
token_file_values=(
  "$setup_script" "$ppg_macro" "$g4_full_list" "$truthjet_full_list"
  "$apply_bdt" "$apply_config" "$base_e_model" "$base_v3e_model"
  "$npb_model" "$tower_mask" "$recoil_runtime_manifest" "$recoil_config"
  "$driver_script" "$worker_self" "$ppg_wrapper" "$recoil_wrapper"
  "$comparator" "$auditor" "$aggregate_extractor"
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
  "$recoil_runtime_manifest" "$recoil_config" "$ppg_wrapper" \
  "$recoil_wrapper" "$auditor" "$comparator" "$aggregate_extractor"; do
  [[ "$path" == /* && -f "$path" && -s "$path" ]] || die "missing input: $path"
done

umask 077
mkdir -p "$output_dir"
state_file="${output_dir}/RUN_STATE"
printf 'RUNNING\n' > "$state_file"
completed=0
finish_state() {
  if [[ $completed -eq 0 ]]; then
    printf 'FAILED\n' > "$state_file"
  fi
}
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

# Remove inherited custom-release routing before sourcing the one common
# runtime.  Keep only the minimal operating-system path needed to run the
# setup script; ordinary SDCC login shells can otherwise retain ana.560 and a
# user install even when OFFLINE_MAIN is later changed to new.17.
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
python3 - "$g4_full_list" "$truthjet_full_list" "$g4_slice" \
  "$truthjet_slice" "$combined_list" <<'PY'
from pathlib import Path
import sys

g4_full, truth_full, g4_out, truth_out, combined_out = map(Path, sys.argv[1:])

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
g4_out.write_text("\n".join(g4) + "\n")
truth_out.write_text("\n".join(truth) + "\n")
combined_out.write_text(
    "\n".join(f"NONE {left} {right} NONE NONE" for left, right in zip(g4, truth))
    + "\n"
)
PY

manifest_role_path() {
  python3 - "$recoil_runtime_manifest" "$1" <<'PY'
import json
import sys
data = json.load(open(sys.argv[1]))
matches = [str(item.get("path", "")) for item in data.get("files", []) if item.get("role") == sys.argv[2]]
if len(matches) != 1:
    raise SystemExit(f"runtime manifest role {sys.argv[2]} has {len(matches)} matches")
print(matches[0])
PY
}

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
recoeff_trace_receipt="$(manifest_role_path ppg_recoeff_trace_transform_receipt)"
recoeff_cross_section_header="$(manifest_role_path ppg_recoeff_cross_section_header)"
recoeff_truth_vertex_header="$(manifest_role_path ppg_recoeff_truth_vertex_header)"
recoeff_canonical_config="$(manifest_role_path ppg_recoeff_canonical_config)"
recoeff_yaml_cpp="$(manifest_role_path ppg_recoeff_yaml_cpp)"
recoeff_roounfold="$(manifest_role_path ppg_recoeff_roounfold)"
recoeff_roounfold_response_header="$(manifest_role_path ppg_recoeff_roounfold_response_header)"
recoeff_roounfold_bayes_header="$(manifest_role_path ppg_recoeff_roounfold_bayes_header)"
recoeff_vertex_scan_data="$(manifest_role_path ppg_recoeff_vertex_scan_data)"
recoeff_mbd_correction="$(manifest_role_path ppg_recoeff_mbd_correction)"
sealed_apply_bdt="$(manifest_role_path ppg_apply_bdt_macro)"
sealed_apply_config="$(manifest_role_path ppg_apply_bdt_config)"
sealed_base_e_model="$(manifest_role_path ppg_apply_model_base_E)"
sealed_base_v3e_model="$(manifest_role_path ppg_apply_model_base_v3E)"
sealed_npb_model="$(manifest_role_path ppg_apply_npb_model)"
recoeff_include_root="$(dirname "$recoeff_cross_section_header")"

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
done

# Produce two byte-identical, relocatable estimator configs.  Their only
# changes from the source-locked canonical file are sealed I/O paths and the
# unique output suffix; all scientific analysis keys remain byte-for-byte.
python3 - "$recoeff_canonical_config" "$baseline_recoeff_config" \
  "$trace_recoeff_config" <<'PY'
from pathlib import Path
import sys

source, baseline, trace = map(Path, sys.argv[1:])
text = source.read_text()
rewrites = (
    ('photon_jet_file_root_dir: "/sphenix/user/shuhangli/ppg12/FunWithxgboost/"',
     'photon_jet_file_root_dir: "input/"'),
    ('eff_outfile: "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency"',
     'eff_outfile: "output/MC_efficiency"'),
    ('response_outfile: "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_response"',
     'response_outfile: "output/MC_response"'),
    ('data_outfile: "/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo"',
     'data_outfile: "output/data_histo"'),
    ('var_type: "bdt_nom"', 'var_type: "paired_oracle"'),
    ('vertex_scan_data_file: ""',
     'vertex_scan_data_file: "input/data_histo_bdt_nom_vtxscan.root"'),
    ('tower_mask_file: "/sphenix/user/shuhangli/ppg12/efficiencytool/tower_masks_bdt_nom.root"',
     'tower_mask_file: "input/tower_masks_bdt_nom.root"'),
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

python3 - \
  "$contract" "$lane_id" "$sample" "$period" "$interaction" \
  "$setup_script" "$calo_calib" "$ppg_macro" "$ppg_lib" \
  "$g4_full_list" "$truthjet_full_list" "$g4_slice" "$truthjet_slice" \
  "$combined_list" "$apply_bdt" "$apply_config" "$base_e_model" \
  "$base_v3e_model" "$npb_model" "$tower_mask" \
  "$recoil_runtime_manifest" "$recoil_macro" "$recoil_config" \
  "$ppg_wrapper" "$recoil_wrapper" "$comparator" "$auditor" "$ppg_log" \
  "$recoil_log" "$ppg_raw_root" "$ppg_scored_root" "$recoil_root" \
  "$candidate_csv" "$expected_token" "$first_ph_seed" \
  "$pedestal_seed" "$pedestal_sequence" "$ph_seed_sequence" \
  "$recoeff_source_macro" "$recoeff_macro" "$recoeff_trace_macro" \
  "$recoeff_trace_receipt" "$recoeff_canonical_config" \
  "$baseline_recoeff_config" "$trace_recoeff_config" \
  "$recoeff_cross_section_header" "$recoeff_truth_vertex_header" \
  "$recoeff_yaml_cpp" "$recoeff_roounfold" \
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
    truth_slice, combined, apply_bdt, apply_config, base_e, base_v3e, npb,
    mask, recoil_manifest, recoil_macro, recoil_config, ppg_wrapper,
    recoil_wrapper, comparator, auditor, ppg_log, recoil_log, ppg_raw, ppg_scored,
    recoil_root, candidate_csv, token, first_ph_seed,
    pedestal_seed, pedestal_sequence, ph_seed_sequence,
    recoeff_source, recoeff_macro, recoeff_trace_macro, recoeff_trace_receipt,
    recoeff_canonical_config, baseline_recoeff_config, trace_recoeff_config,
    cross_section_header, truth_vertex_header, yaml_cpp, roounfold,
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
    "apply_bdt": apply_bdt,
    "apply_config": apply_config,
    "base_e_model": base_e,
    "base_v3e_model": base_v3e,
    "npb_model": npb,
    "tower_mask": mask,
    "recoil_runtime_manifest": recoil_manifest,
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
    "ppg_recoeff_macro": recoeff_macro,
    "ppg_recoeff_trace_macro": recoeff_trace_macro,
    "ppg_recoeff_trace_transform_receipt": recoeff_trace_receipt,
    "ppg_recoeff_canonical_config": recoeff_canonical_config,
    "ppg_recoeff_baseline_config": baseline_recoeff_config,
    "ppg_recoeff_trace_config": trace_recoeff_config,
    "ppg_recoeff_cross_section_header": cross_section_header,
    "ppg_recoeff_truth_vertex_header": truth_vertex_header,
    "ppg_recoeff_yaml_cpp": yaml_cpp,
    "ppg_recoeff_roounfold": roounfold,
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
}
if reuse_mode == "exact_contract_bound":
    token_bound_paths.update({
        "reuse_ppg_raw_root": reuse_raw_root,
        "reuse_ppg_raw_contract": reuse_raw_contract,
    })
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
    "lane": {
        "lane_id": lane_id, "sample": sample, "period": period,
        "interaction": interaction, "rows": 5,
    },
    "rng": {
        "mode": "historical_fifo_replay_v2",
        "reco_consts_randomseed": "absent",
        "ph_seed_sequence": [int(value) for value in ph_seed_sequence.split(",")],
        "first_ph_seed": int(first_ph_seed),
        "pedestal_seed": int(pedestal_seed),
        "pedestal_sequence": int(pedestal_sequence),
        "pedestal_file": f"pedestal-54256-0{int(pedestal_sequence):04d}.root",
        "source_log": "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/photon5/condorout/OutDir0/test.out",
        "source_log_call_count": 5,
    },
    "runtime": {
        "profile": "new.17",
        "offline_main": "/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_new/new.17",
        "actual_offline_main": "/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_new/new.17",
    },
    "runtime_hashes": {
        role: digest(paths[role])
        for role in (
            "setup_script", "calo_calib", "recoil_config", "ppg_wrapper",
            "recoil_wrapper", "comparator",
            "auditor", "aggregate_extractor",
            "ppg_recoeff_baseline_config", "ppg_recoeff_trace_config",
            "apply_bdt_runtime_config",
        )
    },
    "paths": paths,
}
Path(contract).write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
PY

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
  python3 - "$reuse_ppg_raw_contract" "$reuse_ppg_raw_root" \
    "$g4_slice" "$truthjet_slice" "$recoil_runtime_manifest" \
    "$ppg_macro" "$ppg_wrapper" "$setup_script" "$g4_full_list" \
    "$truthjet_full_list" "$worker_self" "$lane_id" "$sample" "$period" \
    "$interaction" <<'PY'
from pathlib import Path
import hashlib
import json
import sys

(
    contract_path, raw_root, g4_slice, truth_slice, runtime_manifest, ppg_macro,
    ppg_wrapper, setup_script, g4_full, truth_full, worker_script,
    lane_id, sample, period, interaction,
) = map(Path, sys.argv[1:])
lane_id, sample, period, interaction = map(
    str, (lane_id, sample, period, interaction)
)
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
prior_reuse = authorization.get("raw_ppg12_reuse", {})
prior_reuse_mode = prior_reuse.get("mode")
if prior_reuse_mode not in {"disabled", "exact_contract_bound"}:
    raise SystemExit("reuse contract records an invalid prior raw-reuse mode")
if prior_reuse_mode == "exact_contract_bound":
    token_file_roles.extend(["reuse_ppg_raw_root", "reuse_ppg_raw_contract"])

prior_paths = data.get("paths", {})
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
if not run_state.is_file() or run_state.read_text().strip() != "PASS":
    raise SystemExit("reuse contract does not belong to a completed PASS run")

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

for key, current in (
    ("setup_script", setup_script),
    ("g4_full_list", g4_full),
    ("truthjet_full_list", truth_full),
    ("recoil_runtime_manifest", runtime_manifest),
    ("ppg_macro", ppg_macro),
    ("ppg_wrapper", ppg_wrapper),
    ("worker_script", worker_script),
):
    item = token_files.get(key, {})
    if not isinstance(item, dict):
        raise SystemExit(f"reuse contract lacks token-bound {key}")
    previous = Path(str(item.get("path", "")))
    recorded = str(item.get("sha256", ""))
    if not previous.is_file() or digest(previous) != recorded:
        raise SystemExit(f"reuse contract recorded {key} no longer matches its source")
    if not current.is_file() or digest(current) != recorded:
        raise SystemExit(f"reuse contract {key} differs from current reconstruction asset")
if raw_root.name != "caloana.root":
    raise SystemExit("reuse source is not the executable CaloAna24 caloana.root artifact")
recorded_raw = Path(str(paths.get("ppg_raw_root", "")))
if recorded_raw.resolve() != raw_root.resolve():
    raise SystemExit("reuse source path differs from the prior schema-3 contract")
apply_evidence = Path(str(paths.get("apply_bdt_stage_evidence", "")))
if not apply_evidence.is_file():
    raise SystemExit("reuse contract lacks completed apply_BDT evidence for the raw input")
apply_data = json.loads(apply_evidence.read_text())
raw_input = apply_data.get("input", {})
if Path(str(raw_input.get("path", ""))).resolve() != raw_root.resolve():
    raise SystemExit("reuse source differs from the prior apply_BDT input")
if str(raw_input.get("sha256", "")) != digest(raw_root):
    raise SystemExit("reuse source hash differs from the prior apply_BDT input hash")
PY
  cp -p "$reuse_ppg_raw_root" "$ppg_raw_root"
  [[ "$(sha256_file "$ppg_raw_root")" == "$(sha256_file "$reuse_ppg_raw_root")" ]] || \
    die "copied PPG12 raw reuse source differs from its token-bound input"
  printf 'REUSED token-bound raw CaloAna24 source %s\n' "$reuse_ppg_raw_root" > "$ppg_log"
else
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

apply_runner="${runtime_dir}/apply_ppg12_bdt.C"
python3 - "$apply_runner" "$apply_bdt" "$apply_runtime_config" "$ppg_raw_root" <<'PY'
from pathlib import Path
import json
import sys
runner, macro, config, root_path = sys.argv[1:]
call = f"apply_BDT({json.dumps(config)},\"data\",{json.dumps(root_path)})"
Path(runner).write_text(f'''{{
  Int_t error = 0;
  if (gROOT->LoadMacro({json.dumps(macro)}) < 0) throw std::runtime_error("apply_BDT load failed");
  gROOT->ProcessLine({json.dumps(call)}, &error);
  if (error != TInterpreter::kNoError) throw std::runtime_error("apply_BDT call failed");
}}
''')
PY
(
  cd "$(dirname "$apply_bdt")"
  export LD_LIBRARY_PATH="${ppg_lib_view}:${base_ld_library_path}"
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
    "$vtxscan_root" "$oracle_sample" <<'PY'
from pathlib import Path
import json
import sys

scan_runner, analysis_runner, macro, config, vtxscan, sample = sys.argv[1:]

def runner(call: str) -> str:
    return f'''{{
  Int_t error = 0;
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

  (
    cd "$layout"
    export LD_LIBRARY_PATH="$(dirname "$recoeff_yaml_cpp"):$(dirname "$recoeff_roounfold"):${base_ld_library_path}"
    export ROOT_INCLUDE_PATH="${recoeff_include_root}:${base_root_include_path}"
    unset RJ_PPG12_EXEC_TRACE_CSV RJ_PPG12_EXEC_RESPONSE_TRACE_CSV
    root -l -b -q "$scan_runner"
  ) 2>&1 | tee "$scan_log"
  [[ -s "$vtxscan_root" ]] || die "$label RecoEff vertex scan is missing"

  (
    cd "$layout"
    export LD_LIBRARY_PATH="$(dirname "$recoeff_yaml_cpp"):$(dirname "$recoeff_roounfold"):${base_ld_library_path}"
    export ROOT_INCLUDE_PATH="${recoeff_include_root}:${base_root_include_path}"
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
run_recoeff baseline "$baseline_layout" "$recoeff_macro" \
  "$baseline_recoeff_config" "$baseline_recoeff_scan_log" \
  "$baseline_recoeff_log" "$baseline_vtxscan_root" 0
run_recoeff trace "$trace_layout" "$recoeff_trace_macro" \
  "$trace_recoeff_config" "$trace_recoeff_scan_log" \
  "$trace_recoeff_log" "$trace_vtxscan_root" 1

for required in "$baseline_eff_root" "$trace_eff_root" \
  "$baseline_response_root" "$trace_response_root" \
  "$ppg_candidate_trace" "$ppg_response_trace"; do
  [[ -s "$required" ]] || die "RecoEff executable evidence is missing: $required"
done

python3 "$aggregate_extractor" \
  --baseline-root "$baseline_eff_root" \
  --instrumented-root "$trace_eff_root" \
  --baseline-response-root "$baseline_response_root" \
  --instrumented-response-root "$trace_response_root" \
  --runtime-manifest "$recoil_runtime_manifest" \
  --root-equivalence-only \
  --out-json "$root_equivalence_report"

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
  --runtime-manifest "$recoil_runtime_manifest" \
  --asset apply_bdt_stage_evidence="$apply_evidence" \
  --asset estimator_source="$recoeff_source_macro" \
  --asset estimator_trace_transform="$recoeff_trace_receipt" \
  --out-json "$aggregate_report"

python3 "$auditor" postrun --contract "$contract"
printf 'PASS\n' > "$state_file"
completed=1
trap - EXIT
printf 'PPG12_PAIRED_ORACLE_PASS output=%s contract=%s report=%s\n' \
  "$output_dir" "$contract" "$report_md"
