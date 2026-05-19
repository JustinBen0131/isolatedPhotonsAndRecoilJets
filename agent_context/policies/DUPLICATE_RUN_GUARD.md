# Duplicate Run Guard

This policy exists because repeated equivalent training/production runs waste
time, Condor slots, disk quota, and user attention.

## When This Applies

Apply before generating, submitting, or preparing:

- AuAu/pp BDT, MLP, stacker, logreg, NPB, sidecar, or ablation training;
- validation on sim, embedded signal/background, data, or score caches;
- RecoilJets production or production-style fanout;
- target-WP/WP80 processing, merge, final ROOT creation, or plot pipelines;
- any "new analysis" that could reuse an existing artifact.

## Run Fingerprint

Build this fingerprint before acting:

- campaign purpose and requested deliverable;
- dataset/mode, signal/background samples, trigger slices, and label semantics;
- model/product name, algorithm, feature list, WP/threshold, isolation policy;
- kinematic range, centrality bins, route/binning scheme, group size;
- config YAML, driver script, source code revision if known;
- input roots/manifests and output roots/tags;
- timestamped campaign tag, cluster/DAG IDs, Gmail subjects;
- expected small products: CSV/JSON/PNG/TXT/ROOT final names.

## Search Order

1. `agent_context/STATUS_DASHBOARD.md`
2. `agent_context/TASK_BOARD.md`
3. `codex_notes/*.md` for older dataset-validity state
4. Local `dataOutput/`, `InputFiles/`, `condor_generated_configs/`
5. Existing model registries under `agent_context/auau_model_registry/`
6. Gmail `RecoilJets Pipeline`, if the run was active/recent
7. Read-only SDCC path checks, if local evidence is insufficient

Use `rg` first. Prefer exact tags, product names, feature-family names, config
names, and output roots over broad prose search.

## Match Classes

- `IDENTICAL`: same dataset, model/features, cuts, labels, WP, route, and output
  intent. Hard stop.
- `MATERIAL_EQUIVALENT`: same science answer with only timestamp/log/root
  location changed. Hard stop unless Justin asks for a deliberate rerun.
- `SUPERSEDED`: older run replaced by a newer validated result. Do not rerun it
  unless the user wants a comparison.
- `RELATED_CONTROL`: useful comparison but not a duplicate. Continue only after
  naming the distinction.
- `UNCLEAR`: gather one more compact evidence check; do not submit while
  unclear.

## Required Response On Duplicate

Say plainly:

> This looks like a duplicate of `<tag/path>`. I am stopping before generating
> or submitting anything.

Then provide the evidence and the useful alternatives:

- reuse existing output;
- compare existing output to a new variant;
- rerun deliberately with a new reason, tag, and cleanup plan.

## Recording

When a new non-duplicate campaign is actually launched or a duplicate is
identified, record the fingerprint/evidence in `STATUS_DASHBOARD.md` or
`TASK_BOARD.md` so future sessions can short-circuit faster.
