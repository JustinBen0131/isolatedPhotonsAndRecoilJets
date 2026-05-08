# Reference Map

This is a compact decision map for local agents. It should point agents to the
right source of truth without copying large PDF content into the workspace
memory.

## Mission Hierarchy

1. PPG19 Au+Au gamma-jet / `x_{J#gamma}` is the ultimate target.
2. PPG18 pp gamma-jet is the validated baseline and comparison anchor.
3. PPG12 provides the photon-ID, BDT/NPB, ABCD, purity, and pp infrastructure
   reference.
4. ATLAS gamma-jet provides the final analysis-shape target.
5. PPG08 dijet xJ provides Au+Au xJ correction/style reference.

## Local References

| Reference | Path | Use For |
|---|---|---|
| PPG12 codebase | `ppg12codeGit/` | Implementation precedent for pp photon ID, BDT/NPB, ABCD, purity, stitching, systematics, and plotting. |
| PPG12 box-cut note | `usefulDocs/PPG12_analysis_note.pdf` | Reference photon-ID box cuts and isolated-photon pp baseline. |
| PPG12 BDT note | `usefulDocs/PPG12_analysis_note_withBDT.pdf` | BDT/NPB features, thresholds, training/application logic, and BDT-based photon-ID variants. |
| ATLAS gamma-jet note | `usefulDocs/Gamma_Jet_Analysis_Note (1).pdf` | Target gamma-jet analysis structure, xJgamma presentation, and final-analysis logic. |
| PPG08 dijet xJ note | `usefulDocs/PPG_08_dijet_xJ_in_Au_Au___draft_Conference_note (7).pdf` | Au+Au xJ correction strategy, presentation style, and dijet-analysis analogy. |
| GL1/GTM manual | `usefulDocs/Gl1-gtm_user_manual_v53.pdf` | Trigger semantics, live/scaled/raw/scaledown interpretation. |

## PPG12 Code Areas To Inspect First

- `ppg12codeGit/CLAUDE.md`: compact overview of the PPG12 pipeline.
- `ppg12codeGit/README.md`: broader PPG12 motivation and workflow.
- `ppg12codeGit/FunWithxgboost/`: BDT training/application infrastructure.
- `ppg12codeGit/efficiencytool/`: ABCD, efficiency, unfolding, and yield
  workflow.
- `ppg12codeGit/plotting/`: final plotting conventions and systematic
  aggregation.
- `ppg12codeGit/simcrosssection/` and related config/scripts: simulation
  cross-section and stitching references.

## Decision Rules

- For photon-ID cuts, start with PPG12 box-cut behavior before Au+Au changes.
- For BDT/NPB cuts, start with PPG12 BDT note and `FunWithxgboost`.
- For pp baseline comparisons, preserve PPG18/PPG12 behavior unless there is a
  documented reason to diverge.
- For Au+Au photon ID, treat new BDT/NPB/JetML variants as extensions that
  require embedded validation and pp comparison where practical.
- For xJgamma, start from the existing RecoilJets pp pipeline, then extend to
  Au+Au with embedded response/closure checks.
- For final presentation, make outputs ATLAS-like in analysis logic while
  retaining sPHENIX style and PPG constraints.
