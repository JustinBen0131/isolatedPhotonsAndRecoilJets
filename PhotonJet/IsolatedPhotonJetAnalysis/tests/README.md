# Tests

The public V1 test suite contains:

- pure unit tests for identities, geometry, selections, weights, features, and
  numerical helpers;
- accepted tree, selection, dataset, plot, and receipt mutation tests;
- tiny ROOT fixtures with expected canonical inventories;
- small tree-starting offline integration tests;
- model-feature and hash-only reference checks;
- repository/projection boundary and machine-path leakage gates.

Producer integration, resolved full-analysis configuration, model inference,
and a clean sPHENIX end-to-end rehearsal are not part of the passing V1 suite;
they remain listed as `NOT_RUN` in the release matrix.

Large production comparisons and migration-only adapters remain outside the
release. See [`../docs/validation.md`](../docs/validation.md).
