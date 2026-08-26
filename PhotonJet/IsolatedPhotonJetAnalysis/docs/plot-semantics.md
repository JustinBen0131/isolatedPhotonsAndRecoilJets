# Provenance-compiled plot semantics

Plot labels are analysis products in this repository. They are never copied
from a note, filename, old macro, or a second configuration file.

The recoil reducer compiles each `RecoilSelection` into a small executable
program. Its predicates perform the numerical comparisons in the event loop.
The histogram command writes the unchanged numeric payload plus a sidecar that
binds those payload bytes to that exact program and to the input-tree hashes.
The plot command verifies both hashes and compiles the human-readable labels
from the program. There is no command-line or Python argument for a cut label.

```text
typed selection
      |
      +---- executable predicates ----> histogram bins
      |
      +---- canonical serialization --> histogram receipt
                                           |
typed dataset manifest --------------------+
                                           v
                                   verified plot contract
                                      |             |
                                      v             v
                              matplotlib PNG   ROOT TLatex include
                                      |
                                      v
                                  render receipt
```

Only participating nodes are displayed. A region-A product records a tight,
isolated event-leading photon and omits the unused non-tight definition. A
region-C or region-D product includes the exact bounded-band or complement
definition because that choice changes the selected rows. Changing any cut
changes the program hash and the compiled label together.

Dataset text is also structured. Collision system, sample class, energy,
centrality, event count, and experiment status come from a validated
`PhotonJetDatasetManifestV1`. Unknown or free-form fields are rejected. Future
trigger, exposure, normalization, correction, and unfolding labels will be
admitted only with their upstream receipt bindings; they will not be inferred
from filenames.

## Canonical style and geometry

ROOT output uses the exact TLatex experiment mark:

```text
#it{#bf{sPHENIX}} Internal
```

Matplotlib renders bold italic `sPHENIX` followed by upright `Internal`.
Individual plots keep all semantic text inside the plotting frame; an external
whitespace header is forbidden. The renderer expands the displayed y-range to
create an internal annotation band and verifies that the complete data/error
envelope remains below it. Comparison tables and deliberately composed
multi-panel summaries may use a separate shared information region. The
renderer fails if a label leaves the frame, overlaps another annotation,
overlaps the data, or falls below the audience-readable font floor. It does
not hide cuts or silently shrink them to make a crowded plot pass.

Every accepted PNG has a `PhotonJetPlotRenderReceiptV1` sidecar containing the
exact visible strings, selection/dataset/style identities, image hash, and
layout results. A PNG without its passing receipt is not a canonical product.

ROOT macros use a generated annotation include plus
`macros/PhotonJetPlotStyle.h`. Run `photonjet plot verify-root-header` against
the same histogram receipt and dataset manifest before using that include. The
verification compares its bytes with a fresh contract compilation and fails if
the include was edited or is stale. Supplying separate handwritten physics
strings remains outside the accepted plotting path.
