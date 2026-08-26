# sPHENIX isolated-photon and recoil-jet analysis

> **Status: downstream differential candidate.** The public tree validator,
> selected offline stages, and provenance-compiled plot path are executable.
> The mini-DST producer, direct-reference emitter, and complete physics
> equivalence matrix are not yet implemented here. This repository is not
> approved for physics production or result reproduction.

This package is the canonical public source for the currently implemented
collaborator-facing stages of the sPHENIX isolated-photon and recoil-jet
workflow for p+p and Au+Au collisions. Exact, receipt-bound projections of a
clean commit can be generated for an analysis-repository subtree or a larger
integration tree. Those projections do not create a second implementation.

The complete production backend remains a separate integration reference until
the producer-to-tree and full downstream equivalence gates have run. This V1
therefore exposes only the supported, evidence-backed boundary described below.

The target workflow is:

```text
DST inputs -> comprehensive TTrees + direct references
           -> optional documented skims
           -> photon purity and background corrections
           -> detector response and unfolding
           -> numeric plot payloads and final figures
```

Photon-identification training, model validation, working-point derivation,
data-product contracts, and provenance are first-class parts of the package.

The current executable boundary is deliberately narrower:

```text
validated PhotonJetTrees_v1
  -> 1D xJgamma payload or documented recoil skim
  -> aggregate ABCD occupancy
  -> provenance-compiled, on-canvas figure annotations
```

The development-only adapter used to compare the accepted tree contract with
older production records is not shipped. It remains private validation
machinery, not a collaborator dependency.

## Design goals

- Preserve the validated physics behavior of the production implementation.
- Give each scientific stage one typed configuration and one documented entry
  point.
- Keep Fun4All and ROOT entry points directly usable.
- Provide a thin `photonjet` command for a guided end-to-end workflow.
- Add resolved configuration, model, and complete run-manifest binding only as
  each stage is implemented and validated; current receipts are stage-specific.
- Test deterministic artifacts byte-for-byte and ROOT content through a
  canonical object representation.
- Keep raw data, production logs, site-specific paths, and migration-only
  machinery outside the release.

## Planned repository map

| Path | Responsibility |
| --- | --- |
| `src/` | Compiled Fun4All producer and ROOT output contracts |
| `macros/` | Direct Fun4All entry points |
| `python/photonjet/` | Configuration, I/O, training, selections, analysis, plotting, and provenance |
| `config/` | Dataset manifests consumed by the current histogram/plot path |
| `contracts/` | Machine-readable definitions of current accepted data products |
| `workflows/` | Inspectable stage graph and stage-level commands |
| `models/` | Model manifests and documentation; model files only when redistribution is appropriate |
| `offline/` | Simple ROOT/Python-facing analysis entry points |
| `examples/` | Minimal p+p and Au+Au walkthroughs |
| `tests/` | Redistribution-safe unit and integration fixtures |
| `docs/` | Physics definitions, tree-contract reference, reproducibility, validation, and troubleshooting |

Directories are populated only when their implementation and tests are ready;
the map is not a promise that every target stage is already executable.

`contracts/` does not preserve superseded internal TTree variants. The first
release will define one accepted collaborator tree contract and only the
corrected nominal trigger and event-selection behavior. Historical adapters
used during differential validation remain outside this repository.

## Release states

The current state is `DOWNSTREAM_DIFFERENTIAL_CANDIDATE`.

1. `ARCHITECTURE_CANDIDATE` — package boundaries and contracts are being
   established; no behavior claim.
2. `DOWNSTREAM_DIFFERENTIAL_CANDIDATE` — selected vertical slices have a
   declared executable comparison matrix, while upstream integration remains
   explicitly not run.
3. `DIFFERENTIAL_IN_PROGRESS` — selected vertical slices reproduce the
   reference implementation.
4. `EQUIVALENCE_CANDIDATE` — the complete declared comparison matrix has run.
5. `RELEASE_CANDIDATE` — the cleaned package passes the matrix.
6. `TRANSFER_READY` — a clean-environment rehearsal and collaborator review
   also pass.

The exact validation policy is documented in
[`docs/validation.md`](docs/validation.md).

## Intended collaboration integration

The standalone repository is organized so its release tree can later be
placed under `PhotonJet/IsolatedPhotonJetAnalysis/` in
[`sPHENIX-Collaboration/analysis`](https://github.com/sPHENIX-Collaboration/analysis)
without carrying development-only material.
