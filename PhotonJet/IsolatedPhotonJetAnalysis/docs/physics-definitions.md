# Physics definitions

The default recoil selection uses an event-leading region-A photon with
(15 leq p_T^gamma < 35) GeV, (|eta^gamma|<0.7), anti-(k_T)
(R=0.4) jets with (p_T^mathrm{jet}>5) GeV and
(|eta^mathrm{jet}|<0.7), and
(Deltaphi>7pi/8). Every boundary is encoded once in
`RecoilSelection` and recorded in the output payload.

ABCD occupancy is event-leading independently in each region:

- A: tight and isolated;
- B: tight and non-isolated;
- C: non-tight and isolated;
- D: non-tight and non-isolated.

Isolated is strictly below the isolated threshold. Non-isolated is strictly
above the non-isolated threshold. A candidate on either boundary is excluded
from both adjacent regions. The non-tight score state must be chosen as either
the bounded score band or the full tight-score complement.

The response module distinguishes measured support from classification
support and represents matched, fake, miss, photon-feed, and xJ boundary
categories explicitly. The unfolding module uses the frozen low-iteration
candidate range and keeps high iterations diagnostic-only.

The public collision-system encoding is singular: p+p stores
`centrality = -1` as the not-applicable state, while Au+Au stores a percentile
in `[0, 100)`. Mixed or alternative encodings fail validation.
