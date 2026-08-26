# Architecture

The package follows one rule: advanced capability is welcome when its boundary
is scientifically meaningful and easy to inspect.

This repository is the canonical public source for its declared V1 surface.
`tools/project_release.py` copies that clean, committed source into exact
analysis and integration subtrees while binding file hashes and executable
modes. The complete production backend remains a separate reference until its
upstream and full downstream equivalence gates have run.

## Stages

```text
dataset manifest
      |
      v
Fun4All producer --------------------> direct-reference ROOT
      |
      v
comprehensive versioned TTrees
      |
      +------------------+
      |                  |
      v                  v
direct-equivalent     documented skim
reduction                or cache
      |                  |
      +---------+--------+
                |
                v
training / selections / purity / response / unfolding
                |
                v
canonical numeric plot payload
                |
                v
selection + dataset contract
                |
                v
receipted final figure / generated ROOT annotation
```

Each currently implemented stage can be invoked directly. The `photonjet`
command transparently parses its documented arguments and dispatches to the
same Python functions. V1 does not yet resolve the unimplemented full-analysis
configuration or emit a generic run manifest.

## C++ and Python boundary

The C++ layer owns detector-node access, event construction, occurrence-safe
identities, direct-reference ROOT objects, and comprehensive TTree writing.

The Python layer owns portable tree I/O, training and model evaluation,
selections, purity, response construction, unfolding orchestration, plot
payloads, figure rendering, and provenance inspection. ROOT macros remain
available where they are the clearest collaboration interface.

Pure definitions—bin edges, region semantics, model variable order, matching
rules, and contract versions—must not be independently redefined on both sides.
They are generated from or validated against one versioned contract.

The public tree contract starts from the accepted collaborator release. Known
wrong or incomplete historical tree layouts and pre-fix trigger/event-gate
behavior are migration inputs only and are not public modes.

## System specialization

p+p and Au+Au share pure physics calculations, accepted output contracts,
histogram declarations, identifiers, configuration types, and offline modules.
Their different detector nodes, event gates, centrality handling, background
subtraction, and model choices remain explicit system adapters.

## Provenance

The complete target provenance contract will bind:

- software commit and release state;
- resolved configuration hash;
- input manifest and content identities;
- accepted data-product contract versions;
- model and working-point identities;
- output paths and hashes;
- validator results and known caveats.

Current V1 commands instead write stage-specific histogram, skim, plot, and
annotation receipts where documented. The generic `PhotonJetRunManifestV1`
helper is an unwired contract primitive; configuration/model/working-point
binding through it is explicitly not run. Scheduler operations and site
forensics are outside this package.

Plot annotations follow the same rule. The selection program that executes the
cuts is serialized beside the numeric payload; a renderer can compile labels
from that program but cannot accept parallel physics prose. This keeps the
implementation small while making data/label drift a validation failure.
