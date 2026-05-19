# Storage And Cleanup

## Core Rule

Quota and cleanup hygiene are part of active analysis work, not an afterthought.
During Condor watches, training, validation, repeated heartbeats, failed runs,
or major plot pulls, periodically check storage pressure where relevant.
After major submissions, failed attempts, resubmits, plot batches, pulls, and
model validations, remove temporary manifests/dry-run clutter Codex created,
record protected timestamped artifacts, and remind Justin which stale artifacts
are ready for review/removal.

## Important Areas

- Local Mac: `/Users/patsfan753/Desktop/ThesisAnalysis`, especially
  `dataOutput/`, `InputFiles/`, large generated plot/model folders.
- SDCC checkout: `/sphenix/u/patsfan753/scratch/thesisAnalysis` and related
  GPFS model roots.
- Bulk area: `/sphenix/tg/tg01/bulk/...`, especially training-tree staging,
  RecoilJets fanout roots, merge trees, and superseded production attempts.

## Failure Recovery

Default behavior: keep small evidence, remove bulky failed payloads after a
fixed smoke/rerun succeeds or a clean replacement has been submitted.

Keep: command, log tail, traceback/held reason, manifest/metadata JSON, small
PNG/CSV/TXT diagnostics.

Remove after safe: failed smoke roots, partial validation caches, failed merge
outputs, `.failed_*`, `.quarantine_*`, `smoke_*`, `tmp_*`, partial model roots.

Quarantine is only for uncertainty windows. Once the corrected path exists and
evidence is preserved, delete quarantined bulky payloads when safe and record
the deletion/evidence path in `STATUS_DASHBOARD.md`.

## Cleanup Classification

Before deleting/offlining anything, classify:

- `KEEP_PROTECTED`: current analysis paths or live jobs.
- `KEEP_SMALL_MANIFESTS`: CSV/JSON/TXT/PNG evidence worth preserving.
- `OFFLINE_OR_COMPRESS`: bulky but potentially useful model/output packs.
- `SAFE_TO_REMOVE_AFTER_APPROVAL`: stale logs, scratch, failed smokes, and
  superseded partial outputs.

Ask Justin before deleting scientific outputs or trained models.

## Common Cleanup Candidates

- `stdout/`, `error/`, `log/` after evidence is consumed;
- `tmp_recoil_merge_*`, dry-run manifests, temporary merge notes;
- failed/superseded `output_*` production attempts;
- old `condor_sub/`, `condor_snapshots/`, `condor_recovery/`,
  `condor_segments/`, `condor_generated_configs/`, `condor_yaml_overrides/`
  once no live DAG/job references them;
- superseded `bdt_models/`, `mlp_models/`, `logreg_models/`,
  `local_bdt_training_outputs/`, `local_ml_pipeline_tests/`, tarballs;
- duplicate local/remote `dataOutput/`, `OutDir/`, `output/`, `outputSmoke/`,
  `transfer/`, `transfer_lite/` products.

## Guarded Helpers

For RecoilJets cleanup, start with dry runs:

```bash
./scripts/recoiljets_cleanup.sh dryrun dataset <dataset>
./scripts/recoiljets_cleanup.sh apply dataset <dataset>
./scripts/recoiljets_cleanup.sh dryrun smoke
./scripts/recoiljets_cleanup.sh apply smoke
```

Use broad `dryrun all` / `apply all` only for intentional broad fresh passes.
Never imply cleanup is safe while live Condor jobs may reference outputs.

For SDCC checkout logs, after consuming useful evidence and confirming live jobs
do not need them, the guarded log cleanup from the intended checkout is:

```bash
find stdout error log -type f -delete
```

Run it only from the intended SDCC `ThesisAnalysis` checkout.

## Bulk Vs Final Areas

Large intermediate pipeline products belong on `/sphenix/tg/tg01/bulk/...`.
Final merged outputs, local analysis-ready ROOT files, pulled products, plots,
and slide inputs belong in `InputFiles/...`, `dataOutput/...`, or explicit
Drive destinations.

When changing pipeline paths, preserve that split: bulk for scalable production
intermediates, user/local analysis areas for final products and
presentation-ready outputs. Do not redirect large per-segment products into the
user's home/basis directories unless Justin explicitly asks.

## AuAu Tight-BDT Sidecar Lifecycle

- Extraction smoke/full roots under
  `/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_<timestamp>/`
  are disposable once row counts, model training, and needed diagnostics have
  been consumed.
- Matching SDCC DAG/submit roots under `condor_sub/auauTightBDT_<timestamp>/`
  are disposable after the DAG has finished and its logs are no longer needed.
- Local SDCC checkout test products under `local_bdt_training_outputs/` are
  disposable after a better smoke/full extraction exists.
- Trained model directories under `bdt_models/tight_<timestamp>/` are not
  disposable unless a newer model is validated and the old one is explicitly
  superseded.

Keep the working checkout tidy during diagnostics: avoid leaving dry-run
`condor_sub/auauTightBDT_*` folders, temporary manifests, or fake-list test
files in the repo after local checks. Prefer timestamped roots for real SDCC
artifacts and remove only generated test clutter Codex created.
