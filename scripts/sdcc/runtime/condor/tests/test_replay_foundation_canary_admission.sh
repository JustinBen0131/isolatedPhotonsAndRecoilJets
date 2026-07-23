#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd -P)"
submitter="${repo_root}/scripts/sdcc/runtime/condor/RecoilJets_Condor_submit.sh"

err() { printf 'ERROR: %s\n' "$*" >&2; }
say() { printf '%s\n' "$*"; }
env_truthy() {
  case "${1:-}" in 1|true|TRUE|yes|YES|on|ON) return 0 ;; *) return 1 ;; esac
}
auto_merge_enabled() {
  case "${RJ_AUTO_MERGE:-1}" in 0|false|FALSE|no|NO|off|OFF) return 1 ;; esac
  return 0
}

# Load only the admission function. The replay-foundation branch returns
# before the production-admission Python validator and its private contracts.
eval "$(sed -n '/^validate_ppg12_stitched_purity_admission()/,/^# Initializes paths for isSim mode/p' "$submitter" | sed '$d')"

reset_contract() {
  ACTION=condorDoAll
  SIM_SAMPLE=run28_jet8
  GROUP_SIZE=1
  GROUP_SIZE_EXPLICIT=1
  MAX_JOBS=1
  MAX_JOBS_EXPLICIT=1
  RJ_AUTO_MERGE=0
  RJ_REPLAY_FOUNDATION_CANARY=1
  RJ_REPLAY_LANE=pp_inclusive_jet8
  RJ_REPLAY_SCHEMA_SHA256=3d32edd4b1a093994c4d70ba7b4480532cc727dc45b05f4ee327aa2797e98a1c
  RJ_DEST_BASE_OVERRIDE=/sphenix/tg/tg01/bulk/example/replay_foundation/the134-test
}

expect_pass() {
  local label="$1"
  validate_ppg12_stitched_purity_admission || {
    printf 'FAIL expected pass: %s\n' "$label" >&2
    exit 1
  }
}

expect_fail() {
  local label="$1"
  if validate_ppg12_stitched_purity_admission; then
    printf 'FAIL expected rejection: %s\n' "$label" >&2
    exit 1
  fi
}

reset_contract
expect_pass 'exact bounded replay-foundation canary'

GROUP_SIZE=5
expect_fail 'group-size widening is rejected'
reset_contract
MAX_JOBS=2
expect_fail 'multiple jobs are rejected'
reset_contract
RJ_AUTO_MERGE=1
expect_fail 'automatic merge is rejected'
reset_contract
unset RJ_REPLAY_LANE
expect_fail 'missing lane identity is rejected'
reset_contract
RJ_REPLAY_SCHEMA_SHA256=not-a-sha256
expect_fail 'malformed schema identity is rejected'
reset_contract
RJ_DEST_BASE_OVERRIDE=/sphenix/tg/tg01/bulk/example/not_replay/the134-test
expect_fail 'non-replay output namespace is rejected'

printf 'PASS replay_foundation_canary_admission mutations=6\n'
