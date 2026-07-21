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

The gate requires all 15 table types, ZSTD compression, exact SHA-256-derived
128-bit identity closure, required branches, and fail-closed duplicate/orphan
fixtures.  Runtime population and direct/writer/replay scientific closure are
separate THE-119 canary gates; a passing synthetic schema does not certify a
physics lane.
