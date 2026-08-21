#!/usr/bin/env bash

set -euo pipefail

test_dir="$(mktemp -d)"
trap 'chmod -R u+w "$test_dir" 2>/dev/null || true; rm -rf "$test_dir"' EXIT
script="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)/launch_supervised_remote_phase.sh"

wait_terminal() {
  local path="$1"
  for _ in $(seq 1 200); do
    [[ -s "$path" ]] && return 0
    sleep 0.05
  done
  return 1
}

success_control="${test_dir}/success.control"
success_output="${test_dir}/success.output"
success_receipt="${success_output}/receipt.json"
"$script" \
  --control-root "$success_control" \
  --output-root "$success_output" \
  --expected-receipt "$success_receipt" \
  --label success_fixture \
  -- /bin/bash -c 'mkdir "$1"; printf "{\"status\":\"PASS\"}\n" >"$2"' \
  fixture "$success_output" "$success_receipt" \
  >"${test_dir}/success.launch"
wait_terminal "$success_control/terminal.json"
python3 - "$success_control/start.json" "$success_control/terminal.json" <<'PY'
import json
import sys

start = json.load(open(sys.argv[1], encoding="utf-8"))
terminal = json.load(open(sys.argv[2], encoding="utf-8"))
assert start["schema"] == "SUPERVISED_REMOTE_PHASE_START_V1"
assert terminal["schema"] == "SUPERVISED_REMOTE_PHASE_TERMINAL_V1"
assert terminal["exit_code"] == 0
assert terminal["signal_number"] == 0
assert terminal["expected_receipt_present"] is True
assert len(terminal["expected_receipt_sha256"]) == 64
assert start["argv_sha256"] == terminal["argv_sha256"]
assert start["wrapper_sha256"] == terminal["wrapper_sha256"]
PY

failure_control="${test_dir}/failure.control"
failure_output="${test_dir}/failure.output"
failure_receipt="${failure_output}/receipt.json"
"$script" \
  --control-root "$failure_control" \
  --output-root "$failure_output" \
  --expected-receipt "$failure_receipt" \
  --label failure_fixture \
  -- /bin/bash -c 'mkdir "$1"; printf "first bad\n" >&2; exit 7' \
  fixture "$failure_output" \
  >"${test_dir}/failure.launch"
wait_terminal "$failure_control/terminal.json"
python3 - "$failure_control/terminal.json" <<'PY'
import json
import sys

terminal = json.load(open(sys.argv[1], encoding="utf-8"))
assert terminal["exit_code"] == 7
assert terminal["expected_receipt_present"] is False
assert terminal["command_log_sha256"] is not None
PY
grep -Fx 'first bad' "$failure_control/command.log" >/dev/null

signal_control="${test_dir}/signal.control"
signal_output="${test_dir}/signal.output"
"$script" \
  --control-root "$signal_control" \
  --output-root "$signal_output" \
  --expected-receipt "$signal_output/receipt.json" \
  --label signal_fixture \
  -- /bin/bash -c 'mkdir "$1"; kill -TERM $$' fixture "$signal_output" \
  >"${test_dir}/signal.launch"
wait_terminal "$signal_control/terminal.json"
python3 - "$signal_control/terminal.json" <<'PY'
import json
import sys

terminal = json.load(open(sys.argv[1], encoding="utf-8"))
assert terminal["exit_code"] == 143
assert terminal["signal_number"] == 15
assert terminal["expected_receipt_present"] is False
PY

forward_control="${test_dir}/forward.control"
forward_output="${test_dir}/forward.output"
"$script" \
  --control-root "$forward_control" \
  --output-root "$forward_output" \
  --expected-receipt "$forward_output/receipt.json" \
  --label forward_fixture \
  -- /bin/bash -c 'mkdir "$1"; sleep 30' fixture "$forward_output" \
  >"${test_dir}/forward.launch"
python3 - "$forward_control/start.json" <<'PY'
import json
import os
import signal
import sys

start = json.load(open(sys.argv[1], encoding="utf-8"))
os.kill(start["supervisor"]["pid"], signal.SIGTERM)
PY
wait_terminal "$forward_control/terminal.json"
python3 - "$forward_control/terminal.json" <<'PY'
import json
import sys

terminal = json.load(open(sys.argv[1], encoding="utf-8"))
assert terminal["exit_code"] == 143
assert terminal["signal_number"] == 15
assert terminal["expected_receipt_present"] is False
PY

if "$script" \
  --control-root "$success_control" \
  --output-root "${test_dir}/duplicate.output" \
  --expected-receipt "${test_dir}/duplicate.output/receipt.json" \
  --label duplicate_fixture -- /usr/bin/true >/dev/null 2>&1; then
  echo "duplicate control root was accepted" >&2
  exit 1
fi

existing_output="${test_dir}/existing.output"
mkdir "$existing_output"
if "$script" \
  --control-root "${test_dir}/existing.control" \
  --output-root "$existing_output" \
  --expected-receipt "$existing_output/receipt.json" \
  --label existing_fixture -- /usr/bin/true >/dev/null 2>&1; then
  echo "existing output root was accepted" >&2
  exit 1
fi

printf 'SUPERVISED_REMOTE_PHASE_LAUNCHER_TEST_PASS\n'
