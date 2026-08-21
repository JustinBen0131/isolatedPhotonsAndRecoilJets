# ML Helpers

Local and SDCC-transferable ML code is organized into `contracts/`,
`training/`, `validation/`, `working_points/`, `stacking/`, and `audits/`.

THE-134 factorial-view tooling uses one shared contract in
`contracts/the134_h70_contract.py`. The same commands are invoked with
`--view H70|H0|G70|G0|O70|O0|R70`; outputs and certificates are view-qualified.
For the source-complete extraction, invoke `prepare_the134_h70_matrix.py` once
per system with `--view H70 --factorial-audit-directory <fresh-dir>`.  That
single process opens each sidecar once, validates all seven views in memory,
writes one shared matrix, and emits seven view-qualified audits bound to the
same matrix SHA-256.  `run_the134_h70_model.py` then projects the requested
typed view into a private working matrix before any p+p row cap or training.
This replaces seven full GPFS scans without changing rows, ordering, labels,
weights, feature values, splits, model parameters, or working-point authority.
The supported sequence is:

```text
the134_full_extraction_controller.py
  -> prepare_the134_h70_matrix.py
  -> run_the134_h70_model.py
  -> materialize_the134_h70_holdout.py
  -> validate_the134_h70_model.py
  -> derive_the134_h70_wp.py
  -> build_the134_h70_model_registry.py
  -> build_the134_source_complete_replay_certificates.py
  -> assemble_the134_science_freeze_manifest.py
  -> build_the134_science_freeze_certificate.py
  -> project_the134_storage_quota.py
```

Before full-scope extraction, `build_the134_source_authority_manifest.py`
measures the exact 13-row, 65-list source inventory and emits a pinned
`THE134_FULL_EXTRACTION_SOURCE_AUTHORITY_V1` candidate receipt.
`build_the134_immutable_bundle_receipt.py` independently rehashes the complete
runtime artifact/dependency inventory. The non-submitting
`resolve_the134_full_multiview_extraction.py` accepts only those two pinned
receipts and produces the exact row descriptors; none of these preflight
surfaces grants scientific completion or production authority.

`the134_full_extraction_controller.py` is the guarded remote execution surface.
It reconstructs the frozen resolver and materializer inputs, seals each
submission attempt, binds every terminal proc to its exact staged chunk, and
emits the source-complete logical aggregate only after the duplicate, quota,
immutable-bundle, job, artifact, provenance, and exact-input-once gates pass.

The extraction executor then requires one row-specific
`THE134_SOURCE_PROVENANCE_V1` JSON derived from the resolved authority. It pins
every input's URI/file/manifest/config/code hashes and the SHA-256 of the
corresponding training-view ROOT file. Each required source must also declare
its full manifest hash, exact expected input and occurrence counts, and the
canonical hash of the complete input-record set. The extractor checks those
bindings independently, including the encoded source role/sample code and the
system-specific weight-component ledger (`event_weight == weight_final`, one
application, and exact declared component product). A required source with
zero 15--35 GeV rows is accepted only through a source-complete
`SOURCE_COMPLETE_ZERO_IN_DOMAIN` authority entry; a smoke tuple cannot satisfy
that gate.

p+p H70 is reuse-only from THE-116, and Au+Au H0 is reuse-only from THE-111.
`audit_the134_model_reuse.py` hashes the referenced immutable artifacts and
emits the reuse certificate consumed by `run_the134_h70_model.py`; every other
system/view pair is trained under THE-134.  Reuse is accepted only when every
artifact matches the system-specific frozen predecessor hash set; the p+p
reuse lane still invokes the established trainer in cache-only mode to
materialize and audit exact PPG12 weights without retraining.  Frozen
comparison controls are p+p H70 and Au+Au H0. Au+Au working points use the full combined weighted source
population via `materialize_the134_h70_wp_sample.py`; p+p uses the exact
event-group holdout.

Every downstream holdout, model-validation, working-point, and registry step
requires the model-completion receipt.  That receipt binds the XGBoost model,
TMVA export, metadata, weighted matrix, and full extraction authority; the
registry rejects any cross-model working-point or validation substitution.

Reuse preserves each accepted predecessor split boundary: THE-116 p+p H70 is
event-grouped 50/50 with seed 42, while THE-111 Au+Au H0 is an immutable
legacy candidate-row 90/10 control with seed 13. The latter is never
reconstructed or presented as new event-grouped H70 training; new Au+Au H70
and all other newly trained Au+Au views remain event-grouped 90/10.

The runner is a local/worker command surface, not a Condor submitter.  Its
default is a non-executing plan; `--execute` must be supplied by an already
authorized, duplicate-guarded campaign controller.

After both system matrices pass full source-population closure,
`build_the134_factorial_execution_plan.py` binds the two single-read matrix
manifests, all fourteen view-qualified extraction audits, and the exact p+p
H70/Au+Au H0 reuse certificates into one deterministic fourteen-lane plan.
It emits argv arrays only and never trains or submits anything.  Missing
views, early authority, hash drift, implicit retraining of a reuse lane, or a
parallel-view request outside 1--7 fails before a plan is written.
The default execution envelope is two concurrent views; any later increase
must remain a separately admitted executor decision.

`execute_the134_factorial_model_plan.py` is that bounded executor.  It accepts
exactly one system at a time, rehashes the plan and all bound artifacts,
refuses known SDCC login hosts and existing output namespaces, cannot widen
the plan's concurrency, never retries, and stops starting new views after the
first failed lane.  The externally recorded `--plan-sha256` is mandatory and
the plan is rebuilt deterministically before admission.  Without `--execute`
it performs admission only.

`build_the134_science_freeze_certificate.py` consumes one immutable manifest
that binds all seven paired p+p/Au+Au model registries and all fourteen
system/view replay certificates. Every replay certificate must carry the
complete system-specific source inventory with direct, writer, and cache
witnesses. Missing views, systems, sources, hashes, exact working points, or
closure gates fail closed.

`build_the134_source_complete_replay_certificates.py` independently rehashes
all seven source registries and all source/view validation receipts before
emitting the fourteen replay certificates.
`assemble_the134_science_freeze_manifest.py` then binds those certificates to
the seven exact paired model registries. These two deterministic builders are
the only supported bridge from source-complete extraction/model evidence to
the aggregate science-freeze certificate.

`project_the134_storage_quota.py` binds the aggregate science certificate to
source-complete byte projections. It budgets the direct-reference histogram
and normalized replay-TTree families produced by the same DST pass, final
merged products, evidence, and two simultaneous copies of the largest merge
workspace. The receipt passes only with at least 20 percent authoritative
quota headroom and the protected local-hot reserve intact. Neither receipt
promotes physics artifacts to `CANONICAL` or opens production by itself.
