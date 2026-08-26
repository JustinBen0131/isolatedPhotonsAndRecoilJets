# Contributing

This package values scientific traceability and direct readability over clever
abstraction. A change should make the physics behavior, configuration, or
reproduction path easier to audit.

## Before changing behavior

- State which stage, dataset class, accepted contract, and output objects are
  affected.
- Add or update the smallest fixture that exercises the change.
- Preserve the prior output as a reference unless the scientific change is
  intentional and documented.
- Record changes to selections, binning, weights, model variables, working
  points, response semantics, or uncertainty handling in `CHANGELOG.md`.

## Comment style

Comments should explain:

- the physics or detector reason for an operation;
- units, coordinate conventions, ranges, and sentinel values;
- identity, ordering, and one-entry-per-object invariants;
- why a non-obvious selection, match, weight, or numerical guard is required;
- ownership of a ROOT object or branch.

Comments should not record development chronology, personal machine paths, or
temporary debugging narratives.

## Configuration

Every configuration key must have one name, type, unit, default, validation
rule, and contract entry. Environment variables may provide deployment-specific
locations, but they must not silently change scientific defaults.

## Validation

Run the relevant unit and integration tests for every change. A refactor that
claims unchanged physics must also pass the differential checks appropriate to
the affected stage. See [`docs/validation.md`](docs/validation.md).

## Plot semantics

Do not type physics cut values, regions, working points, collision labels, or
sample status into a renderer or ROOT macro. Filtering and visible annotation
must share the executable selection program, and dataset text must come from a
typed manifest. Every promoted image requires its passing render receipt. See
[`docs/plot-semantics.md`](docs/plot-semantics.md).
