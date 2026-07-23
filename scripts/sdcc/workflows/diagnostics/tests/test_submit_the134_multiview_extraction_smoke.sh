#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd -P)"
controller="${repo_root}/scripts/sdcc/workflows/diagnostics/submit_the134_multiview_extraction_smoke.sh"

bash -n "$controller"

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

validate(source)

# Deliberately mutate one call site and prove this validator rejects the drift.
mutated = source.replace(expansion, "env", 1)
try:
    validate(mutated)
except ValueError:
    pass
else:
    raise SystemExit("one-call-site mutation was not rejected")

print("THE134_SUBMITTER_CANARY_IDENTITY_WIRING_PASS guarded_calls=2 mutation_rejected=1")
PY
