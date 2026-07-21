# THE-106 C0-R diagnostic writer and cache-only replay

This directory is an additive, canary-only runtime surface for
`AU-AU-PHOTON-C0-RAW-QA`. It has no production registration and no ability to
submit jobs, run Condor, move pointers, merge files, or promote authority.

The writer consumes the existing disabled-by-default `the106::c0h2` and
`the106::c0rh` callbacks. It persists exact source/pair/delivery observations,
event finalizations, raw-QA candidate values and bits, ordered TMVA inputs and
score bits, and direct fill witnesses. Callback failures are caught inside the
writer and recorded; they never escape into the authoritative computation.

The replay binary reads only the sealed diagnostic cache and a separately
frozen inventory. It does not include or call Fun4All and cannot read DSTs.
It reconstructs the nine-variable PPG12 `pre` and applicable
`tight`/`nonTight` templates from candidate raw values, active triggers, slices,
and tight/tag state. Fill witnesses are compared one-for-one but are not used
as a transaction log or source of bin contents.

## Files

- `the106_c0r_diagnostic_writer.{h,cc}` — RAII registrations and immutable TSV
  cache writer.
- `Fun4All_the106_c0r_AuAu.C` — Arm B wrapper around the unchanged unified
  steering implementation.
- `the106_c0r_replay.cc` — cache-only replay and exact ROOT comparator.
- `the106_c0r_score_replay.cc` — exact same-runtime `RBDT("myBDT",
  model_reference)` evaluation from cached float32 input bits.
- `the106_c0r_freeze_inventory.cc` — freezes the complete bounded raw-QA
  object inventory from Arm A1 after validating every discovered object
  against the authoritative booking contract and before Arm C exists.
- `inventory_header.tsv` — schema header only; it is not an authoritative
  inventory.
- `test_the106_c0r_writer.cc` — in-memory paired-source/event/candidate test.
- `test_the106_c0r_replay.py` — synthetic end-to-end replay, comparison, and
  negative witness-mutation test.

## Isolated build

Build only against the isolated C0-H2/C0-RH framework prefix. Never point this
at or install it into the production prefix.

```sh
make \
  THE106_FRAMEWORK_PREFIX=/path/to/isolated/framework-prefix \
  ROOT_CONFIG=/path/to/root-config
make check \
  THE106_FRAMEWORK_PREFIX=/path/to/isolated/framework-prefix \
  ROOT_CONFIG=/path/to/root-config
```

The resulting `libTHE106C0R.so` is loaded only by the Arm B macro. Arm A1 uses
the normal Au+Au macro with all observers unregistered. Arm A0 uses the sealed
current-science snapshot with the THE-106 surface removed.

## Arm B wrapper

The wrapper takes the normal event limit, paired list, and ROOT destination,
followed by an already-created empty cache directory and sealed provenance
strings. The exact invocation is generated in the C0-R evidence namespace so
that credentials and SDCC paths are not stored here. The writer refuses to
overwrite any cache member and fails before the event loop when activation is
incomplete.

The wrapper includes the same `Fun4All_recoilJets_unified_impl.C` used by the
direct control. Its only semantic difference is the lexical lifetime of the
two neutral registrations.

## Frozen inventory contract

Before Arm C, create a new TSV beginning with `inventory_header.tsv`. Every
row binds one required ROOT object to direct booking/fill source evidence and
contains:

- percent-encoded directory, name, class, title, variable, trigger, tag,
  canonical view, and axis titles;
- exact photon-pT and centrality slice indices;
- exact bin count and binary64 `xmin`/`xmax` bits;
- actual direct `Sumw2` state;
- `required=1`; and
- `content_comparator=BITWISE_SINGLE_PROCESS_UNWEIGHTED`.

The inventory must be sealed from authoritative direct source plus the Arm A1
object contract before Arm C is inspected. The replay program intentionally
does not generate or broaden this inventory.

## Cache-only replay and closure

```sh
the106_c0r_replay replay \
  /isolated/C0-R/cache frozen_inventory.tsv armC.root armC_replay.json

the106_c0r_replay compare \
  frozen_inventory.tsv armA1.root armC.root armA1_vs_armC.json

the106_c0r_replay compare \
  frozen_inventory.tsv armB.root armC.root armB_vs_armC.json

the106_c0r_score_replay \
  /isolated/C0-R/cache same_runtime_score_closure.json
```

Score closure independently hashes the exact model bytes, requires one model
reference, and requires every score row to match the sealed model,
preprocessing, feature-order, and runtime authorities.  The runtime identity
is exactly
`ROOT/<gROOT-version>|TMVA::Experimental::RBDT|model_key=myBDT`.  The feature
order SHA-256 is the hash of the percent-decoded ordered feature names, each
followed by one newline (including after the last name).  Ordered float32
values and the returned score are compared by their exact IEEE-754 bits.
Every score joins to a finalized event, each admitted canonical raw-QA
candidate joins to a score observation, and per-event scored-candidate counts
must close.  A run with no valid `Compute` observation is rejected as vacuous.
The helper writes deterministic JSON for both pass and fail outcomes.

Replay fails closed on incomplete finalization, non-closed-world source state,
observer/invariant/serialization failures, duplicate identities, missing value
or trigger rows, noncanonical candidates, missing inventory contracts,
non-unit witness weights, missing or excess fill witnesses, or any exact ROOT
content/`Sumw2` mismatch.

## Deliberate limits

This helper supports only the nine raw shower-shape variables and the
data-only canonical PPG12 `pre`/`tight`/`nonTight` fanout. It does not support
isolation, ranking, ABCD, purity, efficiency, weights, truth, jets, pairs,
response, unfolding, production authority, or general replay. Score closure is
restricted to the exact `RBDT("myBDT", model_reference)` boundary used by the
frozen configuration; runtime and library identity remain sealed by the arm
manifests.
