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
- SDCC checkout and bulk areas: use the local-only paths recorded in
  `agent_context/local/SDCC_LOCAL_RULES.md` or current campaign status files.
  Do not commit site-specific remote paths or queue snapshots.

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

## Anomalous Path Artifacts

Treat blank-looking, quote-looking, or leading-space paths in the SDCC checkout
as suspicious artifacts, not normal project state. Examples include a path whose
printed name is only whitespace, or a zero-byte file such as ` 2`.

Before cleanup, inspect the exact path with escaped names, type, size, mtime,
and contents if it is a directory. If it is empty or clearly accidental, ask
Justin to remove it or remove it when the current request explicitly authorizes
cleanup. Never use wildcards for these paths; delete only the exact verified
path and run a post-cleanup check.

Treat front-facing agent/control-plane directories on SDCC the same way.
`agent_context/`, `.codex/`, `codex*` archive names, and local-only migration
evidence do not belong in the remote checkout. Keep these records locally under
`agent_context/local/...`; if they are found on SDCC, inspect them with escaped
names and sizes, classify them as cleanup candidates, and remove only after
explicit cleanup approval. Do not move them into another visible agent-named
remote directory.

Before running `mkdir`, `mv`, `cp`, or `rm` on SDCC from Codex-generated
commands, validate every path variable with a guard equivalent to:

```bash
case "$path" in
  ""|[[:space:]]*|*[[:space:]]) exit 2 ;;
esac
```

This prevents accidental directories such as a path named exactly two spaces.

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

Large intermediate pipeline products belong in the SDCC bulk area named by the
local-only SDCC rules or current campaign status.
Final merged outputs, local analysis-ready ROOT files, pulled products, plots,
and slide inputs belong in `InputFiles/...`, `dataOutput/...`, or explicit
Drive destinations.

When changing pipeline paths, preserve that split: bulk for scalable production
intermediates, user/local analysis areas for final products and
presentation-ready outputs. Do not redirect large per-segment products into the
user's home/basis directories unless Justin explicitly asks.

## SDCC Base Architecture

After THE-23 stage 6, SDCC base cleanup should preserve the root as a small
runtime surface. New durable outputs and evidence should use the typed domains
below instead of creating new loose root entries:

- RecoilJets output roots: `runs/recoiljets/{current,archive,smoke}/`.
- Condor/runtime state: `state/condor/{sub,snapshots,recovery,segments,generated_configs,yaml_overrides,lists,logs}/`.
- Inputs and calibration lists: `inputs/{dst_lists,sim_lists,grl,z_vertex,xsec}/`.
- Models: `models/{bdt,mlp,logreg,stack,promotion_pulls}/`.
- Small durable evidence: `evidence/{tables,xsec,qa,cleanup}/`.
- Transient scratch/debug material: `scratch/{tmp,debug,slide_ready,local_sim}/`.

Root-level `output_*`, `tmp_*`, loose `.csv`/`.txt`, `config.log`,
`pythia_xsec_firstPass*`, and cleanup stamps are cleanup candidates unless an
audited local-only manifest records them as an active compatibility hold. Move
historical scientific outputs into canonical homes; do not delete them without
separate explicit approval.

For new campaign-scoped SIM production after THE-23 stage 7, prefer the typed
I/O contract over flat roots:

- SDCC merge staging: `runs/recoiljets/current/<campaign>/`.
- TG bulk signal/background: `thesisAna/recoiljets/<family>/<campaign>/signal`
  and `.../background`.
- Local pulls: ROOT inputs under `InputFiles/`; compact QA, plots, and slide
  assets under `dataOutput/<domain>/<campaign>/`.

Treat flat TG roots such as `thesisAna/simembedded_<campaign>` and
`thesisAna/simembeddedinclusive_<campaign>` as legacy read-compatible paths,
not future defaults. Do not move historical TG bulk roots without a separate
live-job/reference/rollback audit.

## AuAu Tight-BDT Sidecar Lifecycle

- Extraction smoke/full roots in the timestamped SDCC bulk training area are
  disposable once row counts, model training, and needed diagnostics have been
  consumed.
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
