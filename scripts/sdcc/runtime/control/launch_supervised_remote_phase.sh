#!/usr/bin/env bash
# Launch one long-running SDCC phase with durable, atomic process receipts.

set -euo pipefail
umask 022

die() {
  printf 'SUPERVISED_REMOTE_PHASE_FAIL: %s\n' "$*" >&2
  exit 2
}

usage() {
  cat <<'EOF'
Usage:
  launch_supervised_remote_phase.sh \
    --control-root ABS --output-root ABS --expected-receipt ABS \
    --label SAFE_ID -- command [args ...]

The control root and output root must not exist. The expected receipt must be
strictly below the output root. The command is supervised in a detached
session; start.json and terminal.json are written atomically below the control
root. No existing namespace is overwritten or removed.
EOF
}

sha256_file() {
  python3 - "$1" <<'PY'
import hashlib
from pathlib import Path
import sys

digest = hashlib.sha256()
with Path(sys.argv[1]).open("rb") as stream:
    for chunk in iter(lambda: stream.read(1024 * 1024), b""):
        digest.update(chunk)
print(digest.hexdigest())
PY
}

canonical_self() {
  python3 - "$1" <<'PY'
from pathlib import Path
import sys
print(Path(sys.argv[1]).resolve(strict=True))
PY
}

proc_start_ticks() {
  local pid="$1"
  if [[ -r "/proc/${pid}/stat" ]]; then
    python3 - "/proc/${pid}/stat" <<'PY' 2>/dev/null || printf 'unavailable\n'
from pathlib import Path
import sys

payload = Path(sys.argv[1]).read_text(encoding="utf-8")
# /proc/<pid>/stat field 2 is parenthesized and may contain spaces. Fields 3+
# begin after the final ") "; starttime is field 22, i.e. tail index 19.
tail = payload.rsplit(") ", 1)[1].split()
print(tail[19])
PY
  else
    printf 'unavailable\n'
  fi
}

boot_id() {
  if [[ -r /proc/sys/kernel/random/boot_id ]]; then
    tr -d '\n' </proc/sys/kernel/random/boot_id
  else
    printf 'unavailable'
  fi
}

mode="launch"
control_root=""
output_root=""
expected_receipt=""
label=""

if [[ "${1:-}" == "--internal-supervise" ]]; then
  mode="supervise"
  shift
fi

while (($#)); do
  case "$1" in
    --control-root)
      [[ $# -ge 2 ]] || die "--control-root requires a value"
      control_root="$2"
      shift 2
      ;;
    --output-root)
      [[ $# -ge 2 ]] || die "--output-root requires a value"
      output_root="$2"
      shift 2
      ;;
    --expected-receipt)
      [[ $# -ge 2 ]] || die "--expected-receipt requires a value"
      expected_receipt="$2"
      shift 2
      ;;
    --label)
      [[ $# -ge 2 ]] || die "--label requires a value"
      label="$2"
      shift 2
      ;;
    --)
      shift
      break
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "unknown argument before --: $1"
      ;;
  esac
done

[[ "$control_root" == /* && "$output_root" == /* && "$expected_receipt" == /* ]] || \
  die "control, output, and receipt paths must be absolute"
[[ "$control_root" != *$'\n'* && "$control_root" != *$'\r'* && "$control_root" != *' '* ]] || \
  die "control root contains forbidden whitespace"
[[ "$output_root" != *$'\n'* && "$output_root" != *$'\r'* && "$output_root" != *' '* ]] || \
  die "output root contains forbidden whitespace"
[[ "$expected_receipt" != *$'\n'* && "$expected_receipt" != *$'\r'* && "$expected_receipt" != *' '* ]] || \
  die "expected receipt contains forbidden whitespace"
[[ "$label" =~ ^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$ ]] || die "label is not a safe identifier"
[[ $# -ge 1 ]] || die "a command is required after --"
[[ "$control_root" != "$output_root" ]] || die "control and output roots must differ"

python3 - "$output_root" "$expected_receipt" <<'PY' || \
  die "expected receipt must be strictly below the output root"
from pathlib import Path
import sys

root = Path(sys.argv[1])
receipt = Path(sys.argv[2])
if any(part in {".", ".."} for part in root.parts + receipt.parts):
    raise SystemExit(1)
try:
    relative = receipt.relative_to(root)
except ValueError:
    raise SystemExit(1)
if not relative.parts:
    raise SystemExit(1)
PY

self_path="$(canonical_self "$0")"
wrapper_sha256="$(sha256_file "$self_path")"

if [[ "$mode" == "launch" ]]; then
  [[ ! -e "$control_root" ]] || die "control root already exists: ${control_root}"
  [[ ! -e "$output_root" ]] || die "output root already exists: ${output_root}"
  [[ -d "$(dirname "$control_root")" ]] || die "control-root parent is missing"
  [[ -d "$(dirname "$output_root")" ]] || die "output-root parent is missing"
  mkdir "$control_root"

  python3 - "$control_root/argv.json" "$@" <<'PY'
import hashlib
import json
import os
from pathlib import Path
import sys

target = Path(sys.argv[1])
argv = sys.argv[2:]
digest = hashlib.sha256(b"\0".join(item.encode("utf-8") for item in argv)).hexdigest()
payload = {"argv": argv, "argv_sha256": digest}
temporary = target.with_suffix(".json.tmp")
temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
os.replace(temporary, target)
PY
  chmod a-w "$control_root/argv.json"

  supervisor_args=(
    --internal-supervise
    --control-root "$control_root"
    --output-root "$output_root"
    --expected-receipt "$expected_receipt"
    --label "$label"
    -- "$@"
  )
  if command -v setsid >/dev/null 2>&1; then
    nohup setsid "$self_path" "${supervisor_args[@]}" </dev/null \
      >"${control_root}/supervisor.log" 2>&1 &
  else
    nohup "$self_path" "${supervisor_args[@]}" </dev/null \
      >"${control_root}/supervisor.log" 2>&1 &
  fi
  supervisor_pid=$!
  printf '%s\n' "$supervisor_pid" >"${control_root}/supervisor.pid.tmp"
  mv "${control_root}/supervisor.pid.tmp" "${control_root}/supervisor.pid"
  chmod a-w "${control_root}/supervisor.pid"

  for _ in $(seq 1 300); do
    if [[ -s "${control_root}/start.json" || -s "${control_root}/terminal.json" ]]; then
      break
    fi
    if ! kill -0 "$supervisor_pid" 2>/dev/null; then
      break
    fi
    sleep 0.1
  done
  [[ -s "${control_root}/start.json" ]] || {
    if kill -0 "$supervisor_pid" 2>/dev/null; then
      python3 - "${control_root}/launch_indeterminate.json" "$label" \
        "$supervisor_pid" "$control_root" "$output_root" <<'PY'
import json
import os
from pathlib import Path
import sys

target, label, supervisor_pid, control_root, output_root = sys.argv[1:]
payload = {
    "schema": "SUPERVISED_REMOTE_PHASE_LAUNCH_INDETERMINATE_V1",
    "label": label,
    "reason": "supervisor_alive_without_start_receipt_after_30_seconds",
    "supervisor_pid": int(supervisor_pid),
    "control_root": control_root,
    "output_root": output_root,
    "retry_permitted": False,
}
target = Path(target)
temporary = target.with_suffix(".json.tmp")
temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
os.replace(temporary, target)
PY
      chmod a-w "${control_root}/launch_indeterminate.json"
    fi
    [[ -f "${control_root}/supervisor.log" ]] && tail -n 80 "${control_root}/supervisor.log" >&2
    die "supervisor did not materialize start.json; inspect the existing control root and never retry blindly"
  }

  cat <<EOF
SUPERVISED_REMOTE_PHASE_START_PASS
  label: ${label}
  supervisor_pid: ${supervisor_pid}
  supervisor_start_ticks: $(proc_start_ticks "$supervisor_pid")
  control_root: ${control_root}
  output_root: ${output_root}
  expected_receipt: ${expected_receipt}
  wrapper_sha256: ${wrapper_sha256}
EOF
  exit 0
fi

[[ -d "$control_root" ]] || die "control root is missing in supervisor mode"
[[ -f "$control_root/argv.json" ]] || die "argv.json is missing in supervisor mode"
mkdir "$control_root/.internal_supervisor_claim" 2>/dev/null || \
  die "internal supervisor was already claimed"
[[ ! -e "$control_root/start.json" && ! -e "$control_root/terminal.json" ]] || \
  die "supervisor receipt already exists"
[[ ! -e "$output_root" ]] || die "output root existed before supervised command start"

argv_sha256="$({ python3 - "$@" <<'PY'
import hashlib
import sys
print(hashlib.sha256(b"\0".join(item.encode("utf-8") for item in sys.argv[1:])).hexdigest())
PY
} )"
recorded_argv_sha256="$(python3 - "$control_root/argv.json" <<'PY'
import json
import sys
print(json.load(open(sys.argv[1], encoding="utf-8"))["argv_sha256"])
PY
)"
[[ "$argv_sha256" == "$recorded_argv_sha256" ]] || die "supervisor argv differs from launch receipt"

started_at="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
host_name="$(hostname -f 2>/dev/null || hostname)"
host_boot_id="$(boot_id)"
supervisor_pid="$$"
supervisor_start_ticks="$(proc_start_ticks "$supervisor_pid")"
command_log="${control_root}/command.log"
forwarded_signal=0
child_pid=""

forward_signal() {
  local signal_name="$1"
  local signal_number_value="$2"
  forwarded_signal="$signal_number_value"
  if [[ -n "$child_pid" ]]; then
    kill -"$signal_name" -- "-${child_pid}" 2>/dev/null || \
      kill -"$signal_name" "$child_pid" 2>/dev/null || true
  fi
}

trap 'forward_signal TERM 15' TERM
trap 'forward_signal INT 2' INT

set +e
if command -v setsid >/dev/null 2>&1; then
  setsid "$@" >"$command_log" 2>&1 &
else
  "$@" >"$command_log" 2>&1 &
fi
child_pid=$!
child_start_ticks="$(proc_start_ticks "$child_pid")"
set -e

python3 - "$control_root/start.json" "$label" "$started_at" "$host_name" \
  "$host_boot_id" "$supervisor_pid" "$supervisor_start_ticks" "$child_pid" \
  "$child_start_ticks" "$wrapper_sha256" "$argv_sha256" "$control_root" \
  "$output_root" "$expected_receipt" "$command_log" <<'PY'
import json
import os
from pathlib import Path
import sys

(target, label, started_at, host, boot_id, supervisor_pid,
 supervisor_ticks, child_pid, child_ticks, wrapper_sha256, argv_sha256,
 control_root, output_root, expected_receipt, command_log) = sys.argv[1:]
payload = {
    "schema": "SUPERVISED_REMOTE_PHASE_START_V1",
    "label": label,
    "started_at": started_at,
    "host": host,
    "boot_id": boot_id,
    "supervisor": {"pid": int(supervisor_pid), "proc_start_ticks": supervisor_ticks},
    "child": {"pid": int(child_pid), "proc_start_ticks": child_ticks},
    "wrapper_sha256": wrapper_sha256,
    "argv_sha256": argv_sha256,
    "control_root": control_root,
    "output_root": output_root,
    "expected_receipt": expected_receipt,
    "command_log": command_log,
}
target = Path(target)
temporary = target.with_suffix(".json.tmp")
temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
os.replace(temporary, target)
PY
chmod a-w "$control_root/start.json"

set +e
while true; do
  wait "$child_pid"
  exit_code=$?
  if ! kill -0 "$child_pid" 2>/dev/null; then
    break
  fi
done
set -e
trap - TERM INT
if ((forwarded_signal > 0 && exit_code < 128)); then
  exit_code=$((128 + forwarded_signal))
fi
ended_at="$(date -u +%Y-%m-%dT%H:%M:%SZ)"

inventory="${control_root}/output_inventory.tsv"
python3 - "$output_root" "${inventory}.tmp" <<'PY'
from pathlib import Path
import sys

root = Path(sys.argv[1])
target = Path(sys.argv[2])
rows = []
if root.is_dir():
    for path in root.rglob("*"):
        if path.is_file():
            stat = path.stat()
            rows.append((str(path), stat.st_size, stat.st_mtime_ns))
with target.open("w", encoding="utf-8") as stream:
    for path, size, mtime_ns in sorted(rows):
        stream.write(f"{size}\t{mtime_ns}\t{path}\n")
PY
mv "${inventory}.tmp" "$inventory"
chmod a-w "$inventory"
inventory_sha256="$(sha256_file "$inventory")"

command_log_sha256="null"
if [[ -f "$command_log" ]]; then
  command_log_sha256="$(sha256_file "$command_log")"
fi
receipt_sha256="null"
receipt_present="false"
if [[ -f "$expected_receipt" ]]; then
  receipt_present="true"
  receipt_sha256="$(sha256_file "$expected_receipt")"
fi
signal_number=0
if ((exit_code >= 128)); then
  signal_number=$((exit_code - 128))
fi

python3 - "$control_root/terminal.json" "$label" "$started_at" "$ended_at" \
  "$host_name" "$host_boot_id" "$supervisor_pid" "$supervisor_start_ticks" \
  "$child_pid" "$child_start_ticks" "$exit_code" "$signal_number" \
  "$wrapper_sha256" "$argv_sha256" "$control_root" "$output_root" \
  "$expected_receipt" "$receipt_present" "$receipt_sha256" "$command_log" \
  "$command_log_sha256" "$inventory" "$inventory_sha256" <<'PY'
import json
import os
from pathlib import Path
import sys

(target, label, started_at, ended_at, host, boot_id, supervisor_pid,
 supervisor_ticks, child_pid, child_ticks, exit_code, signal_number,
 wrapper_sha256, argv_sha256, control_root, output_root, expected_receipt,
 receipt_present, receipt_sha256, command_log, command_log_sha256,
 inventory, inventory_sha256) = sys.argv[1:]
payload = {
    "schema": "SUPERVISED_REMOTE_PHASE_TERMINAL_V1",
    "label": label,
    "started_at": started_at,
    "ended_at": ended_at,
    "host": host,
    "boot_id": boot_id,
    "supervisor": {"pid": int(supervisor_pid), "proc_start_ticks": supervisor_ticks},
    "child": {"pid": int(child_pid), "proc_start_ticks": child_ticks},
    "exit_code": int(exit_code),
    "signal_number": int(signal_number),
    "wrapper_sha256": wrapper_sha256,
    "argv_sha256": argv_sha256,
    "control_root": control_root,
    "output_root": output_root,
    "output_root_present": Path(output_root).is_dir(),
    "expected_receipt": expected_receipt,
    "expected_receipt_present": receipt_present == "true",
    "expected_receipt_sha256": None if receipt_sha256 == "null" else receipt_sha256,
    "command_log": command_log,
    "command_log_sha256": None if command_log_sha256 == "null" else command_log_sha256,
    "output_inventory": inventory,
    "output_inventory_sha256": inventory_sha256,
}
target = Path(target)
temporary = target.with_suffix(".json.tmp")
temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
os.replace(temporary, target)
PY
chmod a-w "$control_root/terminal.json"
chmod -R a-w "$control_root"
exit "$exit_code"
