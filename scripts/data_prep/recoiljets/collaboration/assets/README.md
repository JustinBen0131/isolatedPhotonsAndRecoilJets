# PhotonJetTrees

PhotonJetTrees is a ROOT TTree release for p+p and Au+Au photon+jet studies.
It contains data, photon+jet simulation, and inclusive simulation in a common
layout intended for ordinary ROOT event loops.

## Start here

For recoil studies, make a `TChain` named `photonJets` from the product's
`files.txt`. Each entry contains the event, photon, reconstructed jet, pair
kinematics, isolation, photon-identification score, and event weight for one
stored photon-jet pair.

```cpp
TChain pairs("photonJets");
// Add every path listed in data/auau/files.txt.
```

For an event-style loop with photon, jet, and pair arrays, open `eventTree`.
See `examples/quick_start.C` for a complete recoil loop that reads `files.txt`
directly.

## Weight contract and Au+Au inclusive limitation

`event_weight` is a verbatim value copied from the producer TTree. It is not a
package-level promise of a complete analysis weight. In particular, this
release does not apply either of the following downstream factors:

1. the Au+Au embedded-inclusive source-slice stitch factor
   `ownership_effective_cross_section / Npass` in the matching half-open
   maximum-truth-jet ownership window; or
2. the Au+Au centrality factor.

The complete nominal Au+Au embedded-inclusive weight has the exact order
`producer event_weight -> source-slice stitch -> centrality`. Its statistical
second moment uses the square of that complete final weight. Never interpret
the stored `event_weight` alone as that complete weight, and never substitute
`cross section / generated events` for the accepted source-slice factor.

`simulation/auau_embedded_inclusive` combines the Jet12, Jet20, Jet30, and
Jet40 source productions. The v1 package does not carry a globally stable
per-row source-slice identifier or a hash-bound source-stitch/centrality
sidecar. `source_file_index` is local to a part and is not such an identifier.
Therefore normalized Au+Au embedded-inclusive physics output is
`NOT_READY_FROM_PACKAGE_ALONE`. Use these Trees for raw topology, selection,
and relationship studies, or consume a separately supplied, hash-validated
composite analysis-weight payload. If that payload is absent, stop rather than
reconstructing a factor from filenames, part order, or `event_weight`.

## Files

```text
data/pp/files.txt
data/pp/parts/part_000000.root ...
data/auau/files.txt
data/auau/parts/part_000000.root ...
simulation/pp_photonjet/files.txt
simulation/pp_photonjet/parts/part_000000.root ...
simulation/pp_inclusive/files.txt
simulation/pp_inclusive/parts/part_000000.root ...
simulation/auau_embedded_photonjet/files.txt
simulation/auau_embedded_photonjet/parts/part_000000.root ...
simulation/auau_embedded_inclusive/files.txt
simulation/auau_embedded_inclusive/parts/part_000000.root ...
bdt/pp_model.root
bdt/auau_model.root
bdt/variables.txt
BRANCH_REFERENCE.md
examples/quick_start.C
examples/purity_regions.C
examples/truth_matching.C
```

Each ROOT file contains:

- `events`: one entry per stored event, with event-level quantities;
- `eventTree`: one entry per event, with photon, jet, and pair arrays;
- `photons`: one entry per stored photon candidate;
- `jets`: one entry per stored reconstructed jet;
- `photonJets`: one entry per stored photon-jet pair, with the associated event,
  photon, and jet values repeated for direct analysis;
- `truthPhotons`, `truthJets`, and `recoTruthLinks`: simulation truth objects
  and their stored reconstruction links.

The truth trees are empty in data files. See `BRANCH_REFERENCE.md` for branch
names and types.

The part names and `files.txt` order are stable. All parts in one product have
the same TTree schemas and can be added directly to a `TChain`. The examples
include the current chain file number in multi-file object keys, so events from
different parts cannot be mixed accidentally.

## Stored analysis range

The Au+Au `centrality` branch is continuous, so an analysis may define 0-20%,
20-50%, 50-80%, or another centrality interval within the coverage of the
chosen file. `vertex_z` is also stored so the analysis may choose its vertex
cut.

This release contains anti-kT R=0.4 reconstructed jets. R=0.3 and R=0.4 refer
to photon-isolation cones, not additional jet radii. Additional jet radii
require TTree regeneration. Every consumer must select `jet_radius` explicitly;
this remains mandatory if a later additive release carries more than one jet
radius.

## Photon-identification score

Use `bdt_score` as the photon-identification score. It has one meaning in every
file: the included model evaluated from the stored shower measurements built
with the analysis 70 MeV minimum cell-energy threshold. That threshold applies
to the shower measurements used by the BDT; it is not an isolation-tower cut.
Only photon candidates with complete, finite inputs for the included model are
exposed in the photon and photon-jet trees. Dependent photon-jet pairs and
reconstruction links are omitted when their photon is not exposed; no missing
input is replaced or assigned a placeholder score.

- p+p uses the included PPG12 11-input model in `bdt/pp_model.root`;
- Au+Au uses the included centrality-aware 14-input model in
  `bdt/auau_model.root`.

The exact ordered inputs are stored as `bdt_input_00` through
`bdt_input_13` and listed in `bdt/variables.txt`. The per-candidate branches
`bdt_tight_threshold`, `bdt_nontight_low_threshold`, and
`bdt_nontight_high_threshold` contain the established analysis score cuts.
The corresponding `bdt_is_tight`, `bdt_is_nontight`, and
`bdt_is_not_tight` branches are ready-to-use selection flags. For p+p these
are the PPG12 pT-dependent cuts. For Au+Au they are the centrality-dependent
70% tight and 90%-to-80% bounded non-tight cuts.

## Isolation

`iso_r03` and `iso_r04` store the R=0.3 and R=0.4 photon-isolation values.
Each has an isolated threshold and a non-isolated threshold. The interval
between those thresholds is the isolation gap.

## Analysis choices

The release stores objects and relationships without selecting a final physics
working point. An analysis must state its photon pT, vertex, centrality, BDT,
isolation, jet, and recoil cuts.

For an ABCD purity measurement, first apply the stated photon and event
acceptance, then select the event-leading photon independently in A, B, C, and
D. Also state whether non-tight means the full complement of the tight
selection or a bounded score band. `examples/purity_regions.C` demonstrates
both definitions explicitly with visible pT, eta, vertex, and centrality cuts.

## Examples

- `quick_start.C`: selects an event-leading tight, R=0.4-isolated photon and
  fills an explicitly raw/unweighted anti-kT R=0.4 reconstructed xJgamma
  distribution using the complete
  event identity, including the source file index. This is a package example,
  not the frozen PPG12 photon/jet-arbitration contract;
- `purity_regions.C`: computes raw/unweighted event-leading A/B/C/D counts with an
  explicit complement or bounded non-tight definition;
- `truth_matching.C`: joins stored photon-to-truth-photon links and plots the
  reconstructed-to-truth photon-pT ratio. It is link-level QA and does not
  silently choose a unique best match.

## Release scope

The source TTrees remain unchanged. PhotonJetTrees reorganizes the stored
physics quantities into analysis-friendly TTrees and does not apply a final
photon, isolation, centrality, jet, or recoil working point.

Au+Au DATA is an event-gated input sample: every retained event is terminally
valid, has direct GL1 ScaledVector bit 22, and passes
`MinimumBiasInfo::isAuAuMinimumBias()`. Its `run` and true
`scaled_trigger_bits` values come from the exact source-bound event-gate
witness, and the branches `scaled_bit22`, `minimum_bias_pass`, and
`nominal_event_selection_pass` make that nominal Schema10 event gate explicit.
The release receipt and exact provenance are recorded in
`validation/AUAU_RUN_GATE_PROVENANCE.json`.
