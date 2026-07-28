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
  '${arm}:auau_inclusive_embedded:run28_embeddedJet12'; do
  grep -F "$required" "$wrapper" >/dev/null
done

printf 'THE119_WRITER_EXTRA_ENVIRONMENT_TEST_PASS\n'
