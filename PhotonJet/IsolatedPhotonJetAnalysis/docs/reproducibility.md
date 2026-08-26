# Reproducibility contract

A physics result is reproducible only when its inputs, software, resolved
scientific configuration, models, intermediate data contracts, and output
identities are all known. The current V1 is a downstream implementation
candidate, not yet a complete physics-result reproduction package.

## Required identities

- dataset manifest and ordered source identities;
- software commit and build environment;
- resolved analysis configuration;
- accepted TTree and output-object contract versions;
- training dataset and split identity;
- model bytes, variable order, preprocessing, score direction, and working
  point;
- selection, binning, weights, response, unfolding, and uncertainty settings;
- output hashes and validator status.

Paths are locations, not identities. A moved file remains the same input only
when its declared content identity is unchanged.

## Training lifecycle

Training must expose independently runnable stages for dataset construction,
feature validation, split construction, fitting, held-out evaluation, model
export, model-card generation, and working-point derivation. The exported model
manifest binds all inputs and the exact inference contract used by production
and offline validation.

## Run manifests

The histogram, skim, plot, and generated-annotation paths write the
stage-specific receipts documented with those commands. The generic
`PhotonJetRunManifestV1` helper is not wired into the CLI, and current commands
do not claim complete resolved-configuration, training, model, working-point,
or environment provenance. These remain explicit release gates.

`requirements-constraints.txt` records the direct dependency versions exercised
for this V1. It is a constraints file, not a transitive, hash-locked environment
or a sPHENIX/ROOT runtime specification. A clean sPHENIX environment rehearsal
remains not run.
