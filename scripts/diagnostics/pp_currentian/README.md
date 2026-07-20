# PPG12 current-IAN diagnostic interfaces

The stitched-purity closure workflow is fail-closed and has four producer
layers:

1. `materialize_ppg12_paired_oracle_lane_inputs.py` converts one completed
   five-file photon paired oracle into separate reference and candidate
   lane-input sidecars.  Exact fill counts come from executable candidate-row
   flags, are cross-checked against the weighted executable aggregate, and are
   never inferred from TH1 moments.  Both paired sides share the RecoilJets
   orchestration configuration, while the exact period-specific preserved
   estimator configuration named by `ppg_recoeff_period_config` is hash-bound
   through the selected paired runtime manifest and bundle provenance.  The
   tool is intentionally photon-only;
   inclusive A/B/C/D lanes require their own exact fill ledger.
2. `extract_ppg12_stitched_purity_lane.py` creates one hash-bound lane extract.
3. `produce_ppg12_stitched_purity_evidence.py purity-repeat` executes the
   complete truth/raw/corrected estimator twice with one fresh `TRandom3(42)`
   stream per execution and exactly 20,000 toys per bin.
4. `assemble_ppg12_stitched_purity_manifest.py` assembles the lane index,
   global-minimality evidence, manifests, merge audit, and production wrapper consumed by
   `ppg12_stitched_purity_closure_gate.py`.

Materialize a completed paired photon oracle with the analysis Python runtime:

```bash
/Users/patsfan753/Desktop/analysis/env/bin/python3 \
  scripts/diagnostics/pp_currentian/materialize_ppg12_paired_oracle_lane_inputs.py \
  --paired-contract /absolute/path/to/paired_oracle_contract.json \
  --output-dir /absolute/path/to/materialized_lane
```

The output contains `reference/` and `candidate/` copies of
`lane_input.json`, `fill_evidence.json`, and the canonical extracted
`lane.json`, plus one hash-bound `materialization_manifest.json`.  The command
refuses non-PASS paired runs, stale provenance, a non-photon lane, missing raw
fill flags, weighted-cell drift, or an existing output directory.

Run the purity producer with the analysis Python environment:

```bash
/Users/patsfan753/Desktop/analysis/env/bin/python3 \
  scripts/diagnostics/pp_currentian/produce_ppg12_stitched_purity_evidence.py \
  purity-repeat \
  --lane-index /absolute/path/to/lane_index.json \
  --output /absolute/path/to/purity.json
```

Before assembling any manifest, generate the mandatory global-prefix evidence:

```bash
/Users/patsfan753/Desktop/analysis/env/bin/python3 \
  scripts/diagnostics/pp_currentian/assemble_ppg12_stitched_purity_manifest.py \
  coverage-evidence \
  --lane-index /absolute/path/to/lane_index.json \
  --output /absolute/path/to/global_minimality.json
```

All 32 lanes contribute group zero first. Additional five-file groups are
ordered by `(group_index, lane_id, group_id)`; every addition must reduce the
globally missing required reported-cell set, and no group may remain after the
first complete prefix. The assembler re-runs the canonical extractor for each
linked group sidecar and requires the group moment sums to equal the lane
aggregate. Manifests require this file via `--global-minimality-json`; edited,
stale, hand-authored, or non-minimal evidence fails closed.

The historical-production metrics must be generated from the actual ROOT
objects, including the upstream PPG12 MC-efficiency ROOT that contains the
historical unsuffixed B/C/D histograms:

```bash
/Users/patsfan753/Desktop/analysis/env/bin/python3 \
  scripts/diagnostics/pp_currentian/produce_ppg12_stitched_purity_evidence.py \
  historical-comparison \
  --candidate-purity /absolute/path/to/purity.json \
  --historical-final-root /absolute/path/to/Photon_final_bdt_nom_mc.root \
  --historical-abcd-root /absolute/path/to/MC_efficiency_jet_bdt_nom.root \
  --candidate-inclusive-root /absolute/path/to/inclusive.root \
  --candidate-photon-root /absolute/path/to/photon.root \
  --output /absolute/path/to/historical_comparison.json
```

The final PPG12 ROOT contains `gpurity_leak` and `h_leak_B/C/D`, but it does
not contain the upstream unsuffixed B/C/D inputs. Those counts must never be
inferred from the final graph.

## Canonical merge-evidence requirement

Only the assembler's `merge-audit` command is canonical admission evidence.
It re-resolves every candidate lane's `merge_input`, verifies the file SHA-256,
requires 32 unique inputs in deterministic 20-inclusive/12-photon lane order,
and emits the complete `input_links` list plus its set hash. A hand-written
family summary, an arithmetic report without `input_links`, or a report whose
input hashes no longer resolve is not sufficient for admission or production
promotion.

No command in this directory submits Condor jobs, changes `current.json`, or
publishes an artifact.
