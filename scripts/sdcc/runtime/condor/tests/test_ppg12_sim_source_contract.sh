#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd)"
submitter="${repo_root}/scripts/sdcc/runtime/condor/RecoilJets_Condor_submit.sh"
worker="${repo_root}/scripts/sdcc/runtime/condor/RecoilJets_Condor.sh"
tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/ppg12-source-contract.XXXXXX")"
trap 'rm -rf "$tmpdir"' EXIT

err() { printf 'ERROR: %s\n' "$*" >&2; }
say() { printf '%s\n' "$*"; }
env_truthy() {
  case "${1:-}" in
    1|true|TRUE|yes|YES|on|ON) return 0 ;;
    *) return 1 ;;
  esac
}

# Exercise the production function itself without executing the submitter's
# top-level command dispatcher.
eval "$(sed -n '/^validate_ppg12_closure_canary_controls() {/,/^}/p' "$submitter")"
eval "$(sed -n '/^validate_ppg12_sim_source_contract() {/,/^}/p' "$submitter")"
eval "$(sed -n '/^validate_archived_ppg12_di_chunk_list() {/,/^}/p' "$worker")"

single_list="${tmpdir}/single.list"
closure_single_list="${tmpdir}/closure_single.list"
dual_list="${tmpdir}/dual.list"
closure_dual_list="${tmpdir}/closure_dual.list"
mixed_list="${tmpdir}/mixed.list"

printf '%s\t%s\t%s\t%s\t%s\n' \
  /single/calo.root \
  /sphenix/sim/js_pp200_signal/g4hits/a.root \
  /sphenix/sim/js_pp200_signal/nopileup/jets/a.root \
  /sphenix/sim/js_pp200_signal/nopileup/global/a.root \
  /single/mbd.root > "$single_list"

printf '%s\t%s\t%s\t%s\t%s\n' \
  NONE \
  /sphenix/sim/js_pp200_signal/g4hits/a.root \
  /sphenix/sim/js_pp200_signal/nopileup/jets/a.root \
  NONE \
  NONE > "$closure_single_list"

printf '%s\t%s\t%s\t%s\t%s\n' \
  NONE \
  /sphenix/sim/js_pp200_signal_dual/g4hits/a.root \
  /sphenix/sim/js_pp200_signal_dual/nopileup/jets/a.root \
  /sphenix/sim/js_pp200_signal_dual/nopileup/global/a.root \
  NONE > "$dual_list"

printf '%s\t%s\t%s\t%s\t%s\n' \
  NONE \
  /sphenix/sim/js_pp200_signal_dual/g4hits/a.root \
  /sphenix/sim/js_pp200_signal_dual/nopileup/jets/a.root \
  NONE \
  NONE > "$closure_dual_list"

worker_dual_list="${tmpdir}/worker_dual.list"
worker_dual_with_global_list="${tmpdir}/worker_dual_with_global.list"
printf '%s\t%s\t%s\t%s\t%s\n' \
  NONE \
  /sphenix/sim/js_pp200_signal_dual/g4hits/run0028/jet12/G4Hits_a.root \
  /sphenix/sim/js_pp200_signal_dual/nopileup/jets/run0028/jet12/DST_TRUTH_JET_a.root \
  NONE \
  NONE > "$worker_dual_list"
printf '%s\t%s\t%s\t%s\t%s\n' \
  NONE \
  /sphenix/sim/js_pp200_signal_dual/g4hits/run0028/jet12/G4Hits_a.root \
  /sphenix/sim/js_pp200_signal_dual/nopileup/jets/run0028/jet12/DST_TRUTH_JET_a.root \
  /sphenix/sim/js_pp200_signal_dual/nopileup/global/run0028/jet12/DST_GLOBAL_a.root \
  NONE > "$worker_dual_with_global_list"

validate_archived_ppg12_di_chunk_list "$worker_dual_list" jet12 >/dev/null 2>&1 || {
  printf 'FAIL archived DI worker rejected the frozen G4-only source graph\n' >&2
  exit 1
}
if validate_archived_ppg12_di_chunk_list "$worker_dual_with_global_list" jet12 >/dev/null 2>&1; then
  printf 'FAIL archived DI worker accepted a forbidden prebuilt global lane\n' >&2
  exit 1
fi

cp "$dual_list" "$mixed_list"
printf '%s\t%s\t%s\t%s\t%s\n' \
  NONE \
  /sphenix/sim/js_pp200_signal/g4hits/b.root \
  /sphenix/sim/js_pp200_signal/nopileup/jets/b.root \
  /sphenix/sim/js_pp200_signal/nopileup/global/b.root \
  NONE >> "$mixed_list"

clear_flags() {
  unset RJ_PPG12_PHOTON_YIELD RJ_PPG12_PHOTON_YIELD_DOUBLE
  unset RJ_PPG12_PERIOD_STRICT_DI RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4
  unset RJ_PPG12_PPSIM_G4_ONLY
  unset RJ_PPG12_CLOSURE_CANARY RJ_PPG12_CLOSURE_CANARY_ID
  unset RJ_PPG12_DIRECT_SMEAR_REPAIR
  unset RJ_PPG12_PPSIM_REPLAY_SEEDS
  unset RJ_PPG12_PPSIM_EXPECT_PEDESTAL_SEQUENCE
  unset RJ_PPG12_PPSIM_FIXED_RANDOMSEED
  unset RJ_PPG12_PPSIM_FIXED_PEDESTAL_SEQUENCE
  unset RANDOMSEED RJ_PPG12_PEDESTAL_OVERRIDE
  unset RJ_PPG12_DI_ARCHIVED_EXPECT_PEDESTAL
  unset RJ_DISABLE_JES_CDB_AUDIT RJ_SUBMIT_EXTRA_ENV
  unset RJ_REQUIRE_SIM_GLOBAL
}

expect_pass() {
  local label="$1"
  if ! validate_ppg12_sim_source_contract >/dev/null 2>&1; then
    printf 'FAIL expected pass: %s\n' "$label" >&2
    exit 1
  fi
}

expect_fail() {
  local label="$1"
  if validate_ppg12_sim_source_contract >/dev/null 2>&1; then
    printf 'FAIL expected rejection: %s\n' "$label" >&2
    exit 1
  fi
}

clear_flags
SIM_SAMPLE=run28_jet12
SIM_CLEAN_LIST="$single_list"
RJ_PPG12_PHOTON_YIELD=1
expect_pass 'valid SI sample and single-interaction streams'

RJ_SUBMIT_EXTRA_ENV='RANDOMSEED=7'
expect_pass 'ordinary SI behavior is unchanged outside closure-canary mode'

clear_flags
SIM_SAMPLE=run28_jet12_double
SIM_CLEAN_LIST="$dual_list"
RJ_PPG12_PHOTON_YIELD=1
RJ_PPG12_PHOTON_YIELD_DOUBLE=1
RJ_PPG12_PERIOD_STRICT_DI=1
RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4=1
RJ_PPG12_PPSIM_G4_ONLY=1
expect_pass 'valid DI sample, flags, and dual streams'

clear_flags
SIM_SAMPLE=run28_photonjet10
SIM_CLEAN_LIST="$closure_single_list"
RJ_PPG12_PHOTON_YIELD=1
RJ_PPG12_CLOSURE_CANARY=1
RJ_PPG12_CLOSURE_CANARY_ID='photon:photon10:0mrad:si'
RJ_PPG12_PPSIM_REPLAY_SEEDS='2991264730,4256268992,2394322166,874466025,2240380304'
RJ_PPG12_PPSIM_EXPECT_PEDESTAL_SEQUENCE=534
RJ_DISABLE_JES_CDB_AUDIT=1
RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4=1
RJ_PPG12_PPSIM_G4_ONLY=1
expect_pass 'closure SI uses only single-interaction G4Hits and truth-jet inputs'

unset RJ_PPG12_PHOTON_YIELD
expect_fail 'closure canary requires the PPG12 photon-yield path'

clear_flags
SIM_SAMPLE=run28_photonjet10_double
SIM_CLEAN_LIST="$closure_dual_list"
RJ_PPG12_PHOTON_YIELD=1
RJ_PPG12_PHOTON_YIELD_DOUBLE=1
RJ_PPG12_PERIOD_STRICT_DI=1
RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4=1
RJ_PPG12_PPSIM_G4_ONLY=1
RJ_PPG12_CLOSURE_CANARY=1
RJ_PPG12_CLOSURE_CANARY_ID='photon:photon10:0mrad:di'
RJ_PPG12_PPSIM_REPLAY_SEEDS='2991264730,4256268992,2394322166,874466025,2240380304'
RJ_PPG12_PPSIM_EXPECT_PEDESTAL_SEQUENCE=534
RJ_DISABLE_JES_CDB_AUDIT=1
expect_pass 'closure DI uses only dual-interaction G4Hits and truth-jet inputs'

clear_flags
SIM_SAMPLE=run28_photonjet10
SIM_CLEAN_LIST="$closure_single_list"
RJ_PPG12_PHOTON_YIELD=1
RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4=1
RJ_PPG12_PPSIM_G4_ONLY=1
RJ_PPG12_DIRECT_SMEAR_REPAIR=1
expect_pass 'direct repair SI uses the exact deployed G4-only source graph'

clear_flags
SIM_SAMPLE=run28_photonjet10_double
SIM_CLEAN_LIST="$closure_dual_list"
RJ_PPG12_PHOTON_YIELD=1
RJ_PPG12_PHOTON_YIELD_DOUBLE=1
RJ_PPG12_PERIOD_STRICT_DI=1
RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4=1
RJ_PPG12_PPSIM_G4_ONLY=1
RJ_PPG12_DIRECT_SMEAR_REPAIR=1
expect_pass 'direct repair DI uses the exact deployed dual G4-only source graph'

SIM_CLEAN_LIST="$dual_list"
expect_fail 'direct repair rejects a prebuilt global lane'

clear_flags
SIM_SAMPLE=run28_photonjet10
SIM_CLEAN_LIST="$single_list"
RJ_PPG12_PHOTON_YIELD=1
RJ_PPG12_CLOSURE_CANARY=1
RJ_PPG12_CLOSURE_CANARY_ID='photon:photon10:0mrad:si'
RJ_PPG12_PPSIM_REPLAY_SEEDS='2991264730,4256268992,2394322166,874466025,2240380304'
RJ_PPG12_PPSIM_EXPECT_PEDESTAL_SEQUENCE=534
RJ_DISABLE_JES_CDB_AUDIT=1
RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4=1
RJ_PPG12_PPSIM_G4_ONLY=1
expect_fail 'closure SI rejects prebuilt calo, global, and MBD lanes'

clear_flags
SIM_SAMPLE=run28_photonjet10
SIM_CLEAN_LIST="$closure_single_list"
RJ_PPG12_PHOTON_YIELD=1
RJ_PPG12_PPSIM_REPLAY_SEEDS='2991264730,4256268992,2394322166,874466025,2240380304'
RJ_PPG12_PPSIM_EXPECT_PEDESTAL_SEQUENCE=534
expect_fail 'historical replay controls without closure canary'

RJ_PPG12_CLOSURE_CANARY=1
RJ_PPG12_CLOSURE_CANARY_ID='photon:photon10:0mrad:si'
RJ_DISABLE_JES_CDB_AUDIT=1
RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4=1
RJ_PPG12_PPSIM_G4_ONLY=1
RJ_PPG12_PPSIM_REPLAY_SEEDS='1,2,3,4,5'
expect_fail 'closure canary rejects a noncanonical replay sequence'

RJ_PPG12_PPSIM_REPLAY_SEEDS='2991264730,4256268992,2394322166,874466025,2240380304'
unset RJ_DISABLE_JES_CDB_AUDIT
expect_fail 'closure canary requires the pre-graph JES CDB audit to be disabled'

RJ_DISABLE_JES_CDB_AUDIT=1
RANDOMSEED=42
expect_fail 'closure canary rejects inherited RANDOMSEED'

unset RANDOMSEED
RJ_PPG12_CLOSURE_CANARY_ID='unsafe id'
expect_fail 'closure canary rejects an unsafe canary id'

RJ_PPG12_CLOSURE_CANARY_ID='photon:photon10:0mrad:si'
RJ_SUBMIT_EXTRA_ENV='RJ_PPG12_PPSIM_REPLAY_SEEDS=1,2,3,4,5'
expect_fail 'closure canary rejects controls inherited through RJ_SUBMIT_EXTRA_ENV'

clear_flags
SIM_SAMPLE=run28_jet12
SIM_CLEAN_LIST="$single_list"
RJ_PPG12_PHOTON_YIELD=1
RJ_PPG12_PHOTON_YIELD_DOUBLE=1
RJ_PPG12_PERIOD_STRICT_DI=1
RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4=1
RJ_PPG12_PPSIM_G4_ONLY=1
expect_fail 'exact incident: unsuffixed SI sample with DI flags'

SIM_SAMPLE=run28_jet12_double
SIM_CLEAN_LIST="$single_list"
expect_fail 'double sample backed by single-interaction streams'

SIM_CLEAN_LIST="$mixed_list"
expect_fail 'one contaminated row in an otherwise dual list'

clear_flags
SIM_SAMPLE=run28_jet12_double
SIM_CLEAN_LIST="$dual_list"
RJ_PPG12_PHOTON_YIELD=1
RJ_PPG12_PHOTON_YIELD_DOUBLE=1
expect_fail 'double sample missing strict DI rebuild flags'

clear_flags
SIM_SAMPLE=run28_jet12
SIM_CLEAN_LIST="$dual_list"
RJ_PPG12_PHOTON_YIELD=1
expect_fail 'SI sample backed by dual-interaction streams'

# The deterministic closure controls must be copied explicitly into the worker
# environment; inherited submit-shell state alone is not an acceptable proof.
eval "$(sed -n '/^dataset_is_sim_like() {/,/^}/p' "$submitter")"
eval "$(sed -n '/^append_submit_extra_env_var() {/,/^}/p' "$submitter")"
eval "$(sed -n '/^remove_submit_extra_env_var() {/,/^}/p' "$submitter")"
eval "$(sed -n '/^submit_extra_env_var_is_truthy() {/,/^}/p' "$submitter")"
eval "$(sed -n '/^ppg12_period_sim_uses_auto_mix_weight() {/,/^}/p' "$submitter")"
eval "$(sed -n '/^build_submit_extra_env_fragment() {/,/^}/p' "$submitter")"
eval "$(sed -n '/^sim_requires_global_lane() {/,/^}/p' "$submitter")"

clear_flags
DATASET=isSim
RJ_SUBMIT_EXTRA_ENV=''
RJ_PPG12_PHOTON_YIELD=1
RJ_PPG12_CLOSURE_CANARY=1
RJ_PPG12_CLOSURE_CANARY_ID='photon:photon10:0mrad:si'
RJ_PPG12_PPSIM_REPLAY_SEEDS='2991264730,4256268992,2394322166,874466025,2240380304'
RJ_PPG12_PPSIM_EXPECT_PEDESTAL_SEQUENCE=534
RJ_DISABLE_JES_CDB_AUDIT=1
fragment="$(build_submit_extra_env_fragment)"
for expected in \
  RJ_PPG12_CLOSURE_CANARY=1 \
  RJ_PPG12_CLOSURE_CANARY_ID=photon:photon10:0mrad:si \
  'RJ_PPG12_PPSIM_REPLAY_SEEDS=2991264730,4256268992,2394322166,874466025,2240380304' \
  RJ_PPG12_PPSIM_EXPECT_PEDESTAL_SEQUENCE=534 \
  RJ_DISABLE_JES_CDB_AUDIT=1
do
  printf '%s' "$fragment" | tr ';' '\n' | grep -Fxq "$expected" || {
    printf 'FAIL missing explicit worker environment entry: %s\n' "$expected" >&2
    exit 1
  }
done
if sim_requires_global_lane; then
  printf 'FAIL closure canary incorrectly requires DST_GLOBAL\n' >&2
  exit 1
fi
unset RJ_PPG12_CLOSURE_CANARY RJ_PPG12_CLOSURE_CANARY_ID
unset RJ_PPG12_PPSIM_REPLAY_SEEDS RJ_PPG12_PPSIM_EXPECT_PEDESTAL_SEQUENCE
if ! sim_requires_global_lane; then
  printf 'FAIL ordinary PPG12 photon-yield production no longer requires DST_GLOBAL\n' >&2
  exit 1
fi
RJ_PPG12_DIRECT_SMEAR_REPAIR=1
if sim_requires_global_lane; then
  printf 'FAIL direct deployed-smear repair incorrectly requires DST_GLOBAL\n' >&2
  exit 1
fi

printf 'PASS ppg12_sim_source_contract\n'
