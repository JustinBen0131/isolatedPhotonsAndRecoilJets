# Replay foundation schema validation

`RJReplayFoundationV1.h` defines the normalized, identity-linked ROOT tables
used by the THE-114 campaign.  The direct RecoilJets histograms remain the
authoritative detector execution; these tables are selection-neutral companion
data for exact offline replay.

Run the synthetic contract gate with the analysis ROOT runtime:

```bash
clang++ $(/Users/patsfan753/Desktop/analysis/env/bin/root-config --cflags) \
  -std=c++20 -Wall -Wextra -Werror \
  scripts/diagnostics/replay_foundation/the118_schema_canary.cc \
  $(/Users/patsfan753/Desktop/analysis/env/bin/root-config --libs) \
  -o /tmp/the118_schema_canary

/tmp/the118_schema_canary /tmp/the118_schema.root
/Users/patsfan753/Desktop/analysis/env/bin/python3 \
  scripts/diagnostics/replay_foundation/validate_the118_schema.py \
  /tmp/the118_schema.root --output /tmp/the118_validation.json
```

The gate requires all 16 table types, including the seven-definition
`RJShowerFeatureViewV1` payload, ZSTD compression, exact SHA-256-derived
128-bit identity closure, required branches, and fail-closed duplicate/orphan
fixtures.  Runtime population and direct/writer/replay scientific closure are
separate THE-119 canary gates; a passing synthetic schema does not certify a
physics lane.

## THE-134 shower-definition factorial gate

The H70 production gate stores a union of the H0/H70-centered absolute CEMC
tower neighborhoods and reconstructs seven typed definitions (`H70`, `H0`,
`G70`, `G0`, `O70`, `O0`, and `R70`) for every retained loose photon.  The
nominal analysis definition is H70; the other definitions are controls and do
not alter the 15--35 GeV reporting or validated-model domain.

Run the controlled center/phi-seam contract test locally:

```bash
clang++ $(/Users/patsfan753/Desktop/analysis/env/bin/root-config --cflags) \
  -std=c++20 -Wall -Wextra -Werror \
  scripts/diagnostics/replay_foundation/the134_shower_factorial_canary.cc \
  $(/Users/patsfan753/Desktop/analysis/env/bin/root-config --libs) \
  -o /tmp/the134_shower_factorial_canary
/tmp/the134_shower_factorial_canary
```

After the frozen pp and Au+Au libraries are built on SDCC, use
`submit_the134_shower_factorial_canaries.sh preflight|submit|status`.  Its
12-row direct/writer matrix covers one data, photon-simulation, and inclusive-
simulation source per system.  `validate_the134_shower_factorial.py` then
independently rebuilds every view from absolute tower rows, verifies active
model/view identity and TMVA score parity, and requires direct/writer histogram
and `Sumw2` neutrality.  Passing this bounded gate is required before any
THE-121/THE-122 broad submission.
