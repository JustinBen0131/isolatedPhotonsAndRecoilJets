#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd -P)"
controller="${repo_root}/scripts/sdcc/workflows/diagnostics/submit_the119_replay_foundation_canaries.sh"
wrapper="${repo_root}/scripts/sdcc/workflows/diagnostics/submit_the134_shower_factorial_canaries.sh"

bash -n "$controller"
bash -n "$wrapper"
if grep -Eq '(^|[^[:alnum:]_])(/usr/bin/)?(say|afplay)([^[:alnum:]_]|$)|osascript|NSSound' \
  "$controller" "$wrapper"; then
  printf 'forbidden local audio command in THE-119/THE-134 canary surface\n' >&2
  exit 1
fi

die() {
  printf 'fixture failure: %s\n' "$*" >&2
  return 2
}
eval "$(sed -n '/^validate_writer_extra_template(){/,/^}/p' "$controller")"
eval "$(sed -n '/^render_writer_extra(){/,/^}/p' "$controller")"
eval "$(sed -n '/^render_arm_extra(){/,/^}/p' "$controller")"

extra_common_template='SHARED=1'
extra_pp_template='SYSTEM=pp'
extra_auau_template='SYSTEM=auau'
writer_extra_common_template='A=1;B=__OUTPUT__/sidecar.root'
writer_extra_pp_template='WRITER_SYSTEM=pp'
writer_extra_auau_template='WRITER_SYSTEM=auau'

[[ "$(render_arm_extra pp writer /tmp/pp)" == \
  'SHARED=1;SYSTEM=pp;A=1;B=/tmp/pp/sidecar.root;WRITER_SYSTEM=pp' ]]
[[ "$(render_arm_extra auau writer /tmp/auau)" == \
  'SHARED=1;SYSTEM=auau;A=1;B=/tmp/auau/sidecar.root;WRITER_SYSTEM=auau' ]]
[[ "$(render_arm_extra pp direct /tmp/direct)" == \
  'SHARED=1;SYSTEM=pp' ]]

if ( validate_writer_extra_template invalid 'A=1;A=2' ) >/dev/null 2>&1; then
  printf 'duplicate writer environment keys were accepted\n' >&2
  exit 1
fi
if ( validate_writer_extra_template invalid 'A=__UNKNOWN__' ) >/dev/null 2>&1; then
  printf 'unknown writer environment placeholder was accepted\n' >&2
  exit 1
fi
if ( render_writer_extra pp writer '/tmp/unsafe path' ) >/dev/null 2>&1; then
  printf 'unsafe writer output path was accepted\n' >&2
  exit 1
fi
writer_extra_common_template='SHARED=writer-duplicate'
if ( render_arm_extra pp writer /tmp/pp ) >/dev/null 2>&1; then
  printf 'cross-arm duplicate environment keys were accepted\n' >&2
  exit 1
fi
writer_extra_common_template='A=1;B=__OUTPUT__/sidecar.root'

for required in \
  'RJ_THE134_MULTIVIEW_SIDECAR_ONLY_V1=1' \
  'RJ_THE134_MULTIVIEW_TRAINING_FILE=__OUTPUT__/RJPhotonTrainingViewV1.root' \
  "RJ_THE119_PP_EXTRA_ENV_TEMPLATE='RJ_REPLAY_PERIOD=0mrad" \
  "RJ_THE119_AUAU_EXTRA_ENV_TEMPLATE='RJ_REPLAY_PERIOD=AUAU_RUN24" \
  'for arm in direct writer' \
  '${arm}:pp_inclusive_sim:run28_jet8' \
  '${arm}:auau_inclusive_embedded:run28_embeddedJet12' \
  '--terminal-gate-receipt "$RJ_THE119_EVIDENCE_ROOT/terminal_gate_receipt.json"'; do
  grep -F -- "$required" "$wrapper" >/dev/null
done

tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/the134-sidecar-selector.XXXXXX")"
trap 'rm -rf "$tmpdir"' EXIT
eval "$(sed -n '/^require_replacement_namespaces()/,/^}/p' "$wrapper")"
eval "$(sed -n '/^build_sidecar_only_keys()/,/^}/p' "$wrapper")"
eval "$(sed -n '/^require_exact_sidecar_keys()/,/^}/p' "$wrapper")"
eval "$(sed -n '/^validate_terminal_rows()/,/^}/p' "$wrapper")"

replacement_keys="$(build_sidecar_only_keys 0 | paste -sd, -)"
[[ "$replacement_keys" == \
  'direct:pp_inclusive_sim:run28_jet8,writer:pp_inclusive_sim:run28_jet8' ]]
full_keys="$(build_sidecar_only_keys 1 | paste -sd, -)"
[[ "$full_keys" == \
  'direct:pp_inclusive_sim:run28_jet8,direct:auau_inclusive_embedded:run28_embeddedJet12,writer:pp_inclusive_sim:run28_jet8,writer:auau_inclusive_embedded:run28_embeddedJet12' ]]
require_exact_sidecar_keys "$replacement_keys" "$replacement_keys" replacement
if ( require_exact_sidecar_keys \
  'direct:pp_inclusive_sim:run28_jet8' \
  "$replacement_keys" replacement ) >/dev/null 2>&1; then
  printf 'partial p+p replacement selector was accepted\n' >&2
  exit 1
fi
if ( require_exact_sidecar_keys \
  "${replacement_keys},writer:auau_inclusive_embedded:run28_embeddedJet12" \
  "$replacement_keys" replacement ) >/dev/null 2>&1; then
  printf 'widened p+p replacement selector was accepted\n' >&2
  exit 1
fi

replacement_tag=the134_pp_replacement_fixture
RJ_THE134_TAG="$replacement_tag"
RJ_THE134_OUTPUT_ROOT="${tmpdir}/output/${replacement_tag}"
RJ_THE134_EVIDENCE_ROOT="${tmpdir}/evidence/${replacement_tag}"
require_replacement_namespaces preflight
for special_tag in . ..; do
  RJ_THE134_TAG="$special_tag"
  RJ_THE134_OUTPUT_ROOT="${tmpdir}/special/output/${special_tag}"
  RJ_THE134_EVIDENCE_ROOT="${tmpdir}/special/evidence/${special_tag}"
  if ( require_replacement_namespaces preflight ) >/dev/null 2>&1; then
    printf 'special replacement tag %s was accepted\n' "$special_tag" >&2
    exit 1
  fi
done
RJ_THE134_TAG="$replacement_tag"
RJ_THE134_OUTPUT_ROOT="${tmpdir}/canonical/${replacement_tag}"
RJ_THE134_EVIDENCE_ROOT="${tmpdir}/canonical/sub/../${replacement_tag}"
if ( require_replacement_namespaces preflight ) >/dev/null 2>&1; then
  printf 'canonically identical replacement namespaces were accepted\n' >&2
  exit 1
fi
RJ_THE134_OUTPUT_ROOT="${tmpdir}/output/${replacement_tag}"
RJ_THE134_EVIDENCE_ROOT="${tmpdir}/evidence/${replacement_tag}"
mkdir -p "$RJ_THE134_EVIDENCE_ROOT"
if ( require_replacement_namespaces preflight ) >/dev/null 2>&1; then
  printf 'pre-existing replacement evidence namespace was accepted\n' >&2
  exit 1
fi
require_replacement_namespaces status
require_replacement_namespaces validate
require_replacement_namespaces aggregate
rm -rf "$RJ_THE134_EVIDENCE_ROOT"
mkdir -p "$RJ_THE134_OUTPUT_ROOT"
if ( require_replacement_namespaces submit ) >/dev/null 2>&1; then
  printf 'pre-existing replacement output namespace was accepted\n' >&2
  exit 1
fi
rm -rf "$RJ_THE134_OUTPUT_ROOT"
RJ_THE134_OUTPUT_ROOT="${tmpdir}/output/not-the-tag"
if ( require_replacement_namespaces preflight ) >/dev/null 2>&1; then
  printf 'mismatched replacement tag/output namespace was accepted\n' >&2
  exit 1
fi
if (
  unset RJ_THE134_EVIDENCE_ROOT
  require_replacement_namespaces preflight
) >/dev/null 2>&1; then
  printf 'replacement preflight without explicit evidence namespace was accepted\n' >&2
  exit 1
fi

RJ_THE119_EVIDENCE_ROOT="${tmpdir}/terminal"
RJ_THE119_OUTPUT_ROOT="${tmpdir}/terminal-output"
RJ_THE119_TAG=the134_pp_replacement_terminal_fixture
mkdir -p "$RJ_THE119_EVIDENCE_ROOT"
mkdir -p "$RJ_THE119_OUTPUT_ROOT"
printf '%s\n' 'fixture-preflight' > \
  "${RJ_THE119_EVIDENCE_ROOT}/preflight_receipt.txt"
sidecar_pp_replacement_mode=1
condor_q() { return 0; }
condor_history() { printf '%s\n' '4 0'; }
printf '%s\n' \
  "1001 0 2 --output ${RJ_THE119_OUTPUT_ROOT}/direct/pp_inclusive_sim/run28_jet8" \
  "1002 0 2 --output ${RJ_THE119_OUTPUT_ROOT}/writer/pp_inclusive_sim/run28_jet8" \
  > "${RJ_THE119_EVIDENCE_ROOT}/initial_queue.tsv"
validate_terminal_rows
terminal_receipt="${RJ_THE119_EVIDENCE_ROOT}/terminal_gate_receipt.json"
[[ -s "$terminal_receipt" ]]
python3 - "$terminal_receipt" <<'PY'
import json
import pathlib
import sys

payload = json.loads(pathlib.Path(sys.argv[1]).read_text(encoding="utf-8"))
assert payload["schema"] == "THE134_PP_REPLACEMENT_TERMINAL_GATE_V1"
assert payload["row_count"] == 2
assert {row["role"] for row in payload["rows"]} == {"direct", "writer"}
assert len({row["cluster_proc"] for row in payload["rows"]}) == 2
assert all(row["job_status"] == 4 for row in payload["rows"])
assert all(row["exit_code"] == 0 for row in payload["rows"])
PY
printf '%s\n' '1001 0 2 direct-row' > "${RJ_THE119_EVIDENCE_ROOT}/initial_queue.tsv"
if ( validate_terminal_rows ) >/dev/null 2>&1; then
  printf 'one-row replacement terminal receipt was accepted\n' >&2
  exit 1
fi
printf '%s\n' \
  "1001 0 2 --output ${RJ_THE119_OUTPUT_ROOT}/direct/pp_inclusive_sim/run28_jet8" \
  "1002 0 2 --output ${RJ_THE119_OUTPUT_ROOT}/writer/pp_inclusive_sim/run28_jet8" \
  '1003 0 2 duplicate-row' > "${RJ_THE119_EVIDENCE_ROOT}/initial_queue.tsv"
if ( validate_terminal_rows ) >/dev/null 2>&1; then
  printf 'three-row replacement terminal receipt was accepted\n' >&2
  exit 1
fi
printf '%s\n' \
  "1001 0 2 --output ${RJ_THE119_OUTPUT_ROOT}/direct/pp_inclusive_sim/run28_jet8" \
  "1001 0 2 --output ${RJ_THE119_OUTPUT_ROOT}/writer/pp_inclusive_sim/run28_jet8" \
  > "${RJ_THE119_EVIDENCE_ROOT}/initial_queue.tsv"
if ( validate_terminal_rows ) >/dev/null 2>&1; then
  printf 'duplicate Condor identity was accepted\n' >&2
  exit 1
fi
printf '%s\n' \
  "1001 0 2 --output ${RJ_THE119_OUTPUT_ROOT}/direct/pp_inclusive_sim/run28_jet8" \
  "1002 0 2 --output ${RJ_THE119_OUTPUT_ROOT}/direct/pp_inclusive_sim/run28_jet8" \
  > "${RJ_THE119_EVIDENCE_ROOT}/initial_queue.tsv"
if ( validate_terminal_rows ) >/dev/null 2>&1; then
  printf 'duplicate replacement role was accepted\n' >&2
  exit 1
fi
printf '%s\n' \
  "1001 0 2 --output ${RJ_THE119_OUTPUT_ROOT}/direct/pp_inclusive_sim/run28_jet8_evil" \
  "1002 0 2 --output ${RJ_THE119_OUTPUT_ROOT}/writer/pp_inclusive_sim/run28_jet8_evil" \
  > "${RJ_THE119_EVIDENCE_ROOT}/initial_queue.tsv"
if ( validate_terminal_rows ) >/dev/null 2>&1; then
  printf 'suffixed replacement output paths were accepted\n' >&2
  exit 1
fi
printf '%s\n' \
  "1001 0 2 --note=${RJ_THE119_OUTPUT_ROOT}/direct/pp_inclusive_sim/run28_jet8 unrelated-destination" \
  "1002 0 2 --note=${RJ_THE119_OUTPUT_ROOT}/writer/pp_inclusive_sim/run28_jet8 unrelated-destination" \
  > "${RJ_THE119_EVIDENCE_ROOT}/initial_queue.tsv"
if ( validate_terminal_rows ) >/dev/null 2>&1; then
  printf 'embedded replacement output paths were accepted\n' >&2
  exit 1
fi
printf '%s\n' \
  "1001 0 2 ${RJ_THE119_OUTPUT_ROOT}/writer/pp_inclusive_sim/run28_jet8 ${RJ_THE119_OUTPUT_ROOT}/direct/pp_inclusive_sim/run28_jet8" \
  "1002 0 2 --output ${RJ_THE119_OUTPUT_ROOT}/writer/pp_inclusive_sim/run28_jet8" \
  > "${RJ_THE119_EVIDENCE_ROOT}/initial_queue.tsv"
if ( validate_terminal_rows ) >/dev/null 2>&1; then
  printf 'one replacement row containing both role paths was accepted\n' >&2
  exit 1
fi

printf 'THE119_WRITER_EXTRA_ENVIRONMENT_TEST_PASS\n'
