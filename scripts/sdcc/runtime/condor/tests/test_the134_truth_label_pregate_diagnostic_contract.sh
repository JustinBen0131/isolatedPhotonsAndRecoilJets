#!/usr/bin/env bash
set -Eeuo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
repo_root="$(cd "${script_dir}/../../../../.." && pwd -P)"
source_file="${repo_root}/src/RecoilJets.cc"

python3 - "$source_file" <<'PY'
from pathlib import Path
import sys

source_path = Path(sys.argv[1])
source = source_path.read_text()

start_token = "// THE-134 V13 diagnostic-only pre-weight truth-label observer."
end_token = "// End THE-134 V13 diagnostic-only pre-weight truth-label observer."
assert source.count(start_token) == 1
assert source.count(end_token) == 1

start = source.index(start_token)
end = source.index(end_token, start)
block = source[start:end]
code_only = "\n".join(line.split("//", 1)[0] for line in block.splitlines())

stitch_accept = source.index("[pp photon stitch] process event")
weight_gate = source.index(
    "if (!configurePPG12SimEventWeight(sliceFactor, laneCode))",
    start,
)
assert stitch_accept < start < end < weight_gate

required = (
    '"RJ_THE134_TRUTH_LABEL_PREWEIGHT_DIAGNOSTIC_V1"',
    "THE134_PPG12_TRUTH_LABEL_PREWEIGHT_DIAG_V1",
    "ppPhotonSlice == PPG12PhotonSlice::kPhoton20",
    'stage=PRE_WEIGHT_GATE',
    'record=EVENT_CONTEXT',
    'record=CANDIDATE',
    'weight_gate_bypass=0',
    'new CaloRawClusterEval(topNode, "CEMC")',
    "buildPPG12TruthSignalPhotonMap(selectedEvent)",
    "if (haveCaloEval && rc)",
    "classifyRecoPhotonWithPPG12TruthTrack(rc,",
    "std::ostringstream diagnosticLine;",
    'reason = "EMPTY_SIGNAL_MAP"',
    'reason = "NO_MAX_PRIMARY"',
    'reason = "MATCH"',
    'reason = "PRIMARY_NOT_IN_MAP"',
)
missing = [token for token in required if token not in block]
assert not missing, f"missing V13 pre-weight diagnostic contracts: {missing}"

for forbidden in (
    "->Fill(",
    ".Fill(",
    "bumpHistFill(",
    "processCandidates(",
    "configurePPG12SimEventWeight(",
    "Fun4AllReturnCodes::",
):
    assert forbidden not in code_only, (
        f"diagnostic-only pre-weight observer contains forbidden mutation/control "
        f"token: {forbidden}"
    )
assert "return" not in code_only

weight_gate_contract = """if (!configurePPG12SimEventWeight(sliceFactor, laneCode))
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }"""
assert source.count(weight_gate_contract) == 1
assert source[weight_gate:weight_gate + len(weight_gate_contract)] == weight_gate_contract

weight_block = source[weight_gate:source.index("\n  }\n", weight_gate)]
assert "return Fun4AllReturnCodes::ABORTEVENT;" in weight_block

print("PASS: THE-134 truth-label observer is opt-in, pre-weight, and diagnostic-only")
PY
