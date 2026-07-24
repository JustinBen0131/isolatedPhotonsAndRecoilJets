# Replay foundation schema validation

`RJReplayFoundationV1.h` defines the normalized, identity-linked ROOT tables
used by the THE-114 campaign.  The direct RecoilJets histograms remain the
authoritative detector execution; these tables are selection-neutral companion
data for exact offline replay.

Run the synthetic contract gate with the analysis ROOT runtime:

```bash
clang++ $(root-config --cflags) \
  -std=c++20 -Wall -Wextra -Werror \
  scripts/diagnostics/replay_foundation/the118_schema_canary.cc \
  $(root-config --libs) \
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
tower neighborhoods plus every RawCluster-owned tower needed to reproduce the
native shower denominator, centroid, and `et1`--`et4`.  Cells outside the 7x7
union carry a zero grid-membership bit and are replay provenance only; they do
not enter the frozen rectangular sums or moments.  The payload reconstructs
seven typed definitions (`H70`, `H0`, `G70`, `G0`, `O70`, `O0`, and `R70`) for
every retained loose photon.  The nominal analysis definition is H70; the
other definitions are controls and do not alter the 15--35 GeV reporting or
validated-model domain.

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

### Opt-in labeled multiview extraction

`RJPhotonTrainingViewV1.h` joins each accepted legacy photon-training label to
all seven retained shower definitions without changing
`AuAuPhotonIDTrainingTree` or the 16-tree replay inventory.  It is fail-closed
and writes a separate ROOT artifact only when all three settings are explicit:

```text
RJ_REPLAY_FOUNDATION_V1=1
RJ_THE134_MULTIVIEW_TRAINING_V1=1
RJ_THE134_MULTIVIEW_TRAINING_FILE=/separate/path/training_views.root
```

The corresponding legacy p+p or Au+Au training-tree mode must also be active.
The artifact records source/event/candidate/definition identities, exact
ordered pp11 or AuAu14 features, source-role labels, truth provenance, and a
weight-component ledger.  It never owns a working-point or tag decision;
below-15 rows are diagnostic with null WP/tag state.  Run its local schema and
default-inventory gate with:

```bash
clang++ $(/Users/patsfan753/Desktop/analysis/env/bin/root-config --cflags) \
  -std=c++20 -Wall -Wextra -Werror \
  scripts/diagnostics/replay_foundation/the134_multiview_training_canary.cc \
  $(/Users/patsfan753/Desktop/analysis/env/bin/root-config --libs) \
  -o /tmp/the134_multiview_training_canary
/tmp/the134_multiview_training_canary \
  /tmp/the134_default_replay.root /tmp/the134_pp_views.root /tmp/the134_auau_views.root
```

For a populated p+p extraction sidecar, run
`validate_the134_training_sidecar.py` with the exact expected source and hash
identity.  It implements
`RJ_ARTIFACT_HEALTH_PROFILE_V1 / photon_training_multiview_v1`: sidecar byte
size is diagnostic-only, while ROOT health, exact key and branch inventories,
completion, stable identities, seven-view population, feature finiteness,
weight-once closure, and model-domain safety are mandatory.  This does not
change the separate 50 kB gate for analysis/writer replay ROOT outputs.
