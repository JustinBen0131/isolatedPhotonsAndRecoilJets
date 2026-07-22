#!/usr/bin/env bash
set -Eeuo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
repo_root="$(cd "${script_dir}/../../../../.." && pwd -P)"
macro="${repo_root}/macros/Fun4All_recoilJets_unified_impl.C"
the119_launcher="${repo_root}/scripts/sdcc/workflows/diagnostics/submit_the119_replay_foundation_canaries.sh"

python3 - "$macro" "$the119_launcher" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])
text = path.read_text()
launcher = Path(sys.argv[2]).read_text()

required = (
    'RJ_PPG12_CLOSURE_CANARY',
    'RJ_PPG12_CLOSURE_CANARY_ID',
    'RJ_PPG12_PPSIM_REPLAY_SEEDS',
    'RJ_PPG12_PPSIM_EXPECT_PEDESTAL_SEQUENCE',
    'historical PPG12 seed replay requires recoConsts RANDOMSEED to be absent',
    '2991264730,4256268992,2394322166,874466025,2240380304',
    'ppg12ExpectedPedestalSequence = 534',
    'RJ_DISABLE_JES_CDB_AUDIT=1',
    'ppg12ClosureCanaryId.size() <= 128',
    'PPG12 closure canary forbids inherited override',
    'const int naturalSequence = randGen.Integer(3260);',
    'const int sequence = naturalSequence;',
    'naturalSequence != ppg12ExpectedPedestalSequence',
    'ppg12_closure_rng_contract',
    'ppg12_closure_pedestal_natural_sequence',
    'ppg12_closure_phrandomseed_call_consumed',
    'RJ_REPLAY_FOUNDATION_DI_NEUTRALITY_CANARY',
    'RJ_REPLAY_FOUNDATION_DI_NEUTRALITY_CANARY_ID',
    'replay-foundation DI-neutrality canary requires the exact ',
    'PHRandomSeed::LoadSeed(seed);',
    'replay_foundation_di_neutrality_rng_contract',
)
missing = [token for token in required if token not in text]
assert not missing, f"missing closure RNG/provenance tokens: {missing}"

launcher_required = (
    'RJ_REPLAY_FOUNDATION_DI_NEUTRALITY_CANARY=1',
    'RJ_REPLAY_FOUNDATION_DI_NEUTRALITY_CANARY_ID="the119:${row_sample}:${period}"',
    'RJ_DISABLE_JES_CDB_AUDIT=1',
    'RJ_PPG12_PPSIM_REPLAY_SEEDS=2991264730,4256268992,2394322166,874466025,2240380304',
    'RJ_PPG12_PPSIM_EXPECT_PEDESTAL_SEQUENCE=534',
)
missing_launcher = [token for token in launcher_required if token not in launcher]
assert not missing_launcher, f"THE-119 DI launcher omits neutrality controls: {missing_launcher}"
di_start = launcher.index('elif [[ "$interaction" == di ]]')
di_end = launcher.index('\n  fi', di_start)
di_block = launcher[di_start:di_end]
assert 'RJ_PPG12_CLOSURE_CANARY' not in di_block

start = text.index('const bool usePPG12PPSimAuxInputs')
end = text.index('\n#else', start)
block = text[start:end]

ordered = (
    'InputInit();',
    'InputRegister();',
    'new FlagHandler()',
    'Enable::MBDRECO = true;',
    'Mbd_Reco();',
    'std::make_unique<GlobalVertexReco>()',
    'RunSettings(28);',
    'CEMC_Towers();',
    'HCALInner_Towers();',
    'HCALOuter_Towers();',
    'new TimerStats()',
    'InputManagers();',
    'randGen.SetSeed(PHRandomSeed());',
    'const int naturalSequence = randGen.Integer(3260);',
)
closure_start = block.index('if (ppg12ClosureCanary)', block.index('InputRegister();'))
closure_end = block.index('\n        else\n', closure_start)
closure = block[closure_start:closure_end]
prefix_positions = [block.index(token) for token in ordered[:2]]
positions = [closure.index(token) for token in ordered[2:-3]]
tail_positions = [block.index(token) for token in ordered[-3:]]
assert prefix_positions == sorted(prefix_positions)
assert positions == sorted(positions), f"PPG12 graph order changed: {positions}"
assert prefix_positions[-1] < block.index('if (ppg12ClosureCanary)', block.index('InputRegister();'))
assert tail_positions == sorted(tail_positions)
assert 'setenv("RJ_SKIP_CALO_TOWER_STATUS"' not in closure
assert 'if (!ppg12ClosureCanary)' in text

calib_start = text.index('if (usePPG12PPSimRebuildCaloFromG4)\n    {', end)
calib_end = text.index('\n    else if (isSim && caloInputMode == "simdst")', calib_start)
calib = text[calib_start:calib_end]
process_pos = calib.index('Process_Calo_Calib();')
nosplit_pos = calib.index('new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate_PPG12OracleNoSplit")')
assert process_pos < nosplit_pos
assert 'setSubclusterSplitting(false)' in calib
assert 'setOutputClusterNodeName("CLUSTERINFO_CEMC_NO_SPLIT")' in calib
assert 'set_UseTowerInfo(1)' in calib
assert 'if (ppg12ClosureCanary || !usePPG12Fig11G4OnlyRebuild)' in calib

assert 'CALO, GLOBAL, and MBD lanes must be NONE' in text
assert '!env_truthy_local("RJ_PPG12_PPSIM_G4_ONLY")' in text

print('PASS: deterministic PPG12 closure RNG and exact G4 graph are source-guarded')
PY
