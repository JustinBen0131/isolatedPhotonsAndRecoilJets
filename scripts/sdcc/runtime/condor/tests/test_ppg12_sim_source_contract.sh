#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd)"
submitter="${repo_root}/scripts/sdcc/runtime/condor/RecoilJets_Condor_submit.sh"
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
eval "$(sed -n '/^validate_ppg12_sim_source_contract() {/,/^}/p' "$submitter")"

single_list="${tmpdir}/single.list"
dual_list="${tmpdir}/dual.list"
mixed_list="${tmpdir}/mixed.list"

printf '%s\t%s\t%s\t%s\t%s\n' \
  /single/calo.root \
  /sphenix/sim/js_pp200_signal/g4hits/a.root \
  /sphenix/sim/js_pp200_signal/nopileup/jets/a.root \
  /sphenix/sim/js_pp200_signal/nopileup/global/a.root \
  /single/mbd.root > "$single_list"

printf '%s\t%s\t%s\t%s\t%s\n' \
  NONE \
  /sphenix/sim/js_pp200_signal_dual/g4hits/a.root \
  /sphenix/sim/js_pp200_signal_dual/nopileup/jets/a.root \
  /sphenix/sim/js_pp200_signal_dual/nopileup/global/a.root \
  NONE > "$dual_list"

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

clear_flags
SIM_SAMPLE=run28_jet12_double
SIM_CLEAN_LIST="$dual_list"
RJ_PPG12_PHOTON_YIELD=1
RJ_PPG12_PHOTON_YIELD_DOUBLE=1
RJ_PPG12_PERIOD_STRICT_DI=1
RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4=1
RJ_PPG12_PPSIM_G4_ONLY=1
expect_pass 'valid DI sample, flags, and dual streams'

SIM_SAMPLE=run28_jet12
SIM_CLEAN_LIST="$single_list"
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

printf 'PASS ppg12_sim_source_contract\n'
