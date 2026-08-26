# Validation policy

The package can advance only after the complete production backend reproduces
the accepted reference and this public implementation reproduces every
capability in its declared surface. An analysis-repository projection must then
be byte- and mode-identical to the clean public-source commit selected for it.

## Public tests and migration checks

Two kinds of validation are intentionally separated:

- This repository contains durable unit tests, accepted-contract checks, small
  redistribution-safe fixtures, and compact integration tests useful to every
  contributor.
- Exhaustive golden/canonical and canonical/collaborator launch adapters,
  private production manifests, large data comparisons, and migration
  forensics remain outside the release. Their aggregate certificates may be
  published when useful.

The final package therefore stays clean without losing continuity proof.

## Equality levels

Deterministic source, configuration, manifest, model, and numeric payload files
are compared byte-for-byte.

ROOT container bytes can include UUIDs, timestamps, compression, basket layout,
and library-version details that do not alter scientific content. ROOT outputs
are therefore converted to a deterministic canonical record containing every
declared object path, class, title, axis, bin edge, label, entry count, bin
content, flow bin, error, variance/`Sumw2`, tree branch type/shape/order,
identity, and declared metadata. That canonical record is compared and hashed
byte-for-byte. Raw ROOT byte equality is additionally required where the frozen
environment demonstrates it is stable.

Exact floating comparison is the default. Any tolerance must be declared for
one field or object before the candidate run and include a numerical and
physics justification.

## Required layers

1. Build, load, dictionary, and dependency validation.
2. Resolved configuration equivalence.
3. Event and source identity equivalence.
4. Feature, model-score, and working-point equivalence.
5. Direct-reference ROOT equivalence.
6. Comprehensive TTree contract and value equivalence.
7. TTree-only regeneration of every declared direct scientific object.
8. Optional skim or sufficient-statistics closure.
9. Purity, response, unfolding, covariance, and plot-payload equivalence.
10. Final-render content review.
11. Clean-environment build and end-to-end rehearsal using only public docs.

The matrix covers p+p and Au+Au, data and simulation, centrality, weighting,
truth matching, fakes, misses, boundary migrations, and valid empty-object
cases.

## Cleanup rule

After the differential passes, migration-only adapters and compatibility
scaffolding are removed from the release candidate. Comments and configuration
are rewritten around the final public concepts. The complete applicable matrix
must pass again on that cleaned tree before the release state advances.
