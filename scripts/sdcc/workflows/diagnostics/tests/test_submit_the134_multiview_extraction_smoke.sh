#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd -P)"
controller="${repo_root}/scripts/sdcc/workflows/diagnostics/submit_the134_multiview_extraction_smoke.sh"

bash -n "$controller"

tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/the134-tuple-contract.XXXXXX")"
trap 'rm -rf "$tmpdir"' EXIT
die() { return 2; }
eval "$(sed -n '/^validate_one_five_file_tuple()/,/^}/p' "$controller")"
eval "$(sed -n '/^yaml_value()/,/^}/p' "$controller")"
eval "$(sed -n '/^validate_five_field_fanout_contract()/,/^}/p' "$controller")"

printf 'calo\tg4\tjets\tglobal\tmbd\n' > "${tmpdir}/valid.list"
validate_one_five_file_tuple unit "${tmpdir}/valid.list"

printf 'calo\tg4\tjets\tglobal\n' > "${tmpdir}/four-columns.list"
if validate_one_five_file_tuple unit "${tmpdir}/four-columns.list"; then
  printf 'four-column tuple was not rejected\n' >&2
  exit 1
fi

printf 'calo\tg4\tjets\tglobal\tmbd\ncalo2\tg42\tjets2\tglobal2\tmbd2\n' > "${tmpdir}/two-rows.list"
if validate_one_five_file_tuple unit "${tmpdir}/two-rows.list"; then
  printf 'two-row tuple ownership was not rejected\n' >&2
  exit 1
fi

printf 'calo\tg4\t\tglobal\tmbd\n' > "${tmpdir}/empty-column.list"
if validate_one_five_file_tuple unit "${tmpdir}/empty-column.list"; then
  printf 'empty tuple column was not rejected\n' >&2
  exit 1
fi

printf '%s\n' 'coneR: 0.40' > "${tmpdir}/valid.yaml"
printf '%s\n' 'dest|cfg|pre|tight|nonTight' > "${tmpdir}/valid.fanout"
validate_five_field_fanout_contract unit "${tmpdir}/valid.fanout" "${tmpdir}/valid.yaml" 0.40

fanout_mutations=(
  'dest|cfg|pre|tight'
  'dest|cfg|pre|tight|nonTight|0.40'
  'dest|cfg|pre|tight|nonTight|0.40|false|2.0'
  'dest|cfg|pre||nonTight'
)
for index in "${!fanout_mutations[@]}"; do
  printf '%s\n' "${fanout_mutations[$index]}" > "${tmpdir}/bad-fanout-${index}.txt"
  if validate_five_field_fanout_contract unit "${tmpdir}/bad-fanout-${index}.txt" "${tmpdir}/valid.yaml" 0.40; then
    printf 'fanout schema mutation %s was not rejected\n' "$index" >&2
    exit 1
  fi
done

printf '%s\n%s\n' 'dest|cfg|pre|tight|nonTight' 'dest2|cfg2|pre2|tight2|nonTight2' > "${tmpdir}/two-row.fanout"
if validate_five_field_fanout_contract unit "${tmpdir}/two-row.fanout" "${tmpdir}/valid.yaml" 0.40; then
  printf 'two-row fanout ownership was not rejected\n' >&2
  exit 1
fi

if validate_five_field_fanout_contract unit "${tmpdir}/valid.fanout" "${tmpdir}/missing.yaml" 0.40; then
  printf 'missing materialized config was not rejected\n' >&2
  exit 1
fi

printf '%s\n' 'coneR: 0.30' > "${tmpdir}/wrong-cone.yaml"
if validate_five_field_fanout_contract unit "${tmpdir}/valid.fanout" "${tmpdir}/wrong-cone.yaml" 0.40; then
  printf 'wrong materialized cone was not rejected\n' >&2
  exit 1
fi

printf '%s\n%s\n' 'coneR: 0.40' 'coneR: 0.40' > "${tmpdir}/duplicate-cone.yaml"
if validate_five_field_fanout_contract unit "${tmpdir}/valid.fanout" "${tmpdir}/duplicate-cone.yaml" 0.40; then
  printf 'duplicate materialized cone authority was not rejected\n' >&2
  exit 1
fi

python3 - "$controller" <<'PY'
from pathlib import Path
import sys

expansion = 'env "${materialize_contract_env[@]}"'
source = Path(sys.argv[1]).read_text()

def validate(text: str) -> None:
    required = (
        'local -a materialize_contract_env=(',
        'RJ_REPLAY_FOUNDATION_CANARY=1',
        'RJ_REPLAY_LANE="$lane"',
        'RJ_REPLAY_SCHEMA_SHA256="$RJ_THE134_REPLAY_SCHEMA_SHA256"',
    )
    for token in required:
        if token not in text:
            raise ValueError(f"missing submit-shell replay-canary identity: {token}")
    if text.count(expansion) != 2:
        raise ValueError(
            "the shared submit-shell replay-canary identity must guard exactly "
            "the p+p and Au+Au dry-materialization calls"
        )
    worker_contract = (
        'RJ_REPLAY_FOUNDATION_V1=1;RJ_REPLAY_FOUNDATION_CANARY=1;'
        'RJ_REPLAY_TRACE=0;RJ_REPLAY_LANE=${lane}'
    )
    if worker_contract not in text:
        raise ValueError("worker descriptor replay-canary identity was lost")
    receipt_contract = (
        'materialized_config_values=()',
        'descriptor_env_values "$submit_file" RJ_CONFIG_YAML',
        '"${#materialized_config_values[@]}" == 1',
        'materialized_config_sha256',
        'file_sha256(materialized_config) != materialized_config_sha',
    )
    for token in receipt_contract:
        if token not in text:
            raise ValueError(f"missing materialized-config receipt protection: {token}")

validate(source)

# Deliberately mutate one call site and prove this validator rejects the drift.
mutated = source.replace(expansion, "env", 1)
try:
    validate(mutated)
except ValueError:
    pass
else:
    raise SystemExit("one-call-site mutation was not rejected")

print("THE134_SUBMITTER_CANARY_IDENTITY_WIRING_PASS guarded_calls=2 mutation_rejected=1 tuple_mutations=3 fanout_mutations=8")
PY
