#!/usr/bin/env bash
# Foreground-only implementation for run_ppg12_recoiljets_paired_oracle.sh.
# This script is intentionally not a submitter and contains no merge/promotion
# path.  It writes only below a new, explicitly authorized output directory.

set -euo pipefail

die() {
  printf 'PPG12_PAIRED_ORACLE_WORKER_FAIL: %s\n' "$*" >&2
  exit 2
}

[[ $# -eq 16 ]] || die "internal argument-count mismatch"

expected_token="$1"
repo_root="$2"
output_dir="$3"
setup_script="$4"
ppg_macro="$5"
ppg_lib="$6"
g4_full_list="$7"
truthjet_full_list="$8"
apply_bdt="$9"
apply_config="${10}"
base_e_model="${11}"
base_v3e_model="${12}"
npb_model="${13}"
tower_mask="${14}"
recoil_runtime_manifest="${15}"
recoil_config="${16}"

[[ "${RJ_PPG12_PAIRED_ORACLE_RUN_TOKEN:-}" == "$expected_token" ]] || \
  die "missing exact authorization token from plan/run driver"
[[ "$expected_token" == ppg12-oracle:* ]] || die "malformed authorization token"
[[ -z "${RANDOMSEED+x}" ]] || die "inherited shell RANDOMSEED is forbidden"
[[ "$output_dir" == /* && ! -e "$output_dir" ]] || \
  die "output directory must be absolute and absent: $output_dir"

# Start the scientific runtime from a genuinely empty environment.  Merely
# changing OFFLINE_MAIN is insufficient on SDCC: an inherited login shell can
# retain release_ana and user-install routes in PATH/LD_LIBRARY_PATH.
if [[ "${RJ_PPG12_PAIRED_ORACLE_CLEAN_ENV:-0}" != 1 ]]; then
  worker_self="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)/$(basename "${BASH_SOURCE[0]}")"
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

ppg_wrapper="${repo_root}/macros/diagnostics/pp_currentian/Fun4All_ppg12_fixed_seed_oracle.C"
recoil_wrapper="${repo_root}/macros/diagnostics/pp_currentian/Fun4All_recoiljets_fixed_seed_oracle.C"
auditor="${repo_root}/scripts/diagnostics/pp_currentian/audit_ppg12_recoiljets_paired_oracle.py"
comparator="${repo_root}/scripts/diagnostics/pp_currentian/compare_ppg12_recoiljets_same_cluster_features.py"

for path in \
  "$setup_script" "$ppg_macro" "$ppg_lib" "$g4_full_list" \
  "$truthjet_full_list" "$apply_bdt" "$apply_config" "$base_e_model" \
  "$base_v3e_model" "$npb_model" "$tower_mask" \
  "$recoil_runtime_manifest" "$recoil_config" "$ppg_wrapper" \
  "$recoil_wrapper" "$auditor" "$comparator"; do
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
mkdir -p "$input_dir" "$runtime_dir" "$ppg_dir" "$recoil_dir" "$report_dir"

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
calo_reco_lib="$(manifest_role_path libcalo_reco.so)"
clusteriso_lib="$(manifest_role_path libclusteriso.so)"
jetbase_lib="$(manifest_role_path libjetbase.so)"
photon_builder_header="$(manifest_role_path PhotonClusterBuilder.h)"
photon_builder_include_root="$(dirname "$(dirname "$photon_builder_header")")"

# The new.17 setup does not place OFFLINE_MAIN/rootmacros on ROOT's default
# macro path.  Resolve the release-owned calibration macro explicitly and add
# that exact directory to both Cling's include path and ROOT's macro path for
# each executable below.  Dynamic discovery would otherwise fail before the
# oracle runs despite the builder having already sealed this same file.
calo_macro_dir="${OFFLINE_MAIN}/rootmacros"
calo_calib="${calo_macro_dir}/Calo_Calib.C"
[[ -f "$calo_calib" && -s "$calo_calib" ]] || \
  die "exact new.17 Calo_Calib.C is missing: $calo_calib"

ppg_raw_root="${ppg_dir}/output_sim.root"
ppg_scored_root="${ppg_dir}/output_sim_with_bdt_split.root"
recoil_root="${recoil_dir}/recoil.root"
ppg_log="${ppg_dir}/ppg12.log"
recoil_log="${recoil_dir}/recoiljets.log"
apply_log="${ppg_dir}/apply_bdt.log"
report_md="${report_dir}/paired_oracle_report.md"
summary_csv="${report_dir}/paired_oracle_summary.csv"
candidate_csv="${report_dir}/paired_oracle_candidates.csv"
contract="${output_dir}/paired_oracle_contract.json"

python3 - \
  "$contract" "$setup_script" "$calo_calib" "$ppg_macro" "$ppg_lib" \
  "$g4_full_list" "$truthjet_full_list" "$g4_slice" "$truthjet_slice" \
  "$combined_list" "$apply_bdt" "$apply_config" "$base_e_model" \
  "$base_v3e_model" "$npb_model" "$tower_mask" \
  "$recoil_runtime_manifest" "$recoil_macro" "$recoil_config" \
  "$ppg_wrapper" "$recoil_wrapper" "$comparator" "$ppg_log" \
  "$recoil_log" "$ppg_raw_root" "$ppg_scored_root" "$recoil_root" \
  "$candidate_csv" "$expected_token" "$first_ph_seed" \
  "$pedestal_seed" "$pedestal_sequence" "$ph_seed_sequence" <<'PY'
from pathlib import Path
import hashlib
import json
import sys

(
    contract, setup, calo, ppg_macro, ppg_lib, g4_full, truth_full, g4_slice,
    truth_slice, combined, apply_bdt, apply_config, base_e, base_v3e, npb,
    mask, recoil_manifest, recoil_macro, recoil_config, ppg_wrapper,
    recoil_wrapper, comparator, ppg_log, recoil_log, ppg_raw, ppg_scored,
    recoil_root, candidate_csv, token, first_ph_seed,
    pedestal_seed, pedestal_sequence, ph_seed_sequence,
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
    "ppg_log": ppg_log,
    "recoil_log": recoil_log,
    "ppg_raw_root": ppg_raw,
    "ppg_scored_root": ppg_scored,
    "recoil_root": recoil_root,
    "candidate_csv": candidate_csv,
}
data = {
    "schema_version": 2,
    "authorization_token": token,
    "lane": {
        "sample": "Photon5", "period": "1p5mrad", "interaction": "SI",
        "rows": 5,
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
[[ -s "$ppg_raw_root" ]] || die "PPG12 executable did not write output_sim.root"

apply_runner="${runtime_dir}/apply_ppg12_bdt.C"
python3 - "$apply_runner" "$apply_bdt" "$apply_config" "$ppg_raw_root" <<'PY'
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
  export RJ_SIM_SAMPLE=run28_photonjet5
  export RJ_PPG12_CLOSURE_CANARY=1
  export RJ_PPG12_CLOSURE_CANARY_ID=photon:photon5:1p5mrad:si:pairedoracle
  export RJ_PPG12_PPSIM_REPLAY_SEEDS="$ph_seed_sequence"
  export RJ_PPG12_PPSIM_EXPECT_PEDESTAL_SEQUENCE="$pedestal_sequence"
  unset RJ_PPG12_PPSIM_FIXED_RANDOMSEED
  unset RJ_PPG12_PPSIM_FIXED_PEDESTAL_SEQUENCE
  export RJ_PPG12_PHOTON_YIELD=1
  export RJ_PPG12_PHOTON_YIELD_DOUBLE=0
  export RJ_PPG12_PHOTON_YIELD_CLUSTER_ERES=0.04
  export RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4=1
  export RJ_PPG12_PPSIM_G4_ONLY=1
  export RJ_SIM_ALLOW_NONE_LISTS=1
  export RJ_PPG12_PERIOD=1p5mrad
  export RJ_PPG12_CROSSING_PERIOD=1p5mrad
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
  export RJ_PP_PHOTONID_PPG12_FILTER=1
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
  --out-md "$report_md" \
  --out-csv "$summary_csv" \
  --out-candidates-csv "$candidate_csv"

python3 "$auditor" postrun --contract "$contract"
printf 'PASS\n' > "$state_file"
completed=1
trap - EXIT
printf 'PPG12_PAIRED_ORACLE_PASS output=%s contract=%s report=%s\n' \
  "$output_dir" "$contract" "$report_md"
