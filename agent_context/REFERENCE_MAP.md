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
| PPG12 current IAN | `usefulDocs/PPG12_analysis_note_2026-05-21_v4_current_IAN.pdf` | Preferred current PPG12 internal analysis note. Use first for photon-ID, isolation, BDT/NPB, purity, pp baseline, and analysis-note logic. |
| PPG12 current paper draft | `usefulDocs/sPHENIX_PPG12_Paper_2026-05-21_current_draft.pdf` | Preferred current paper-draft reference for final PPG12 wording, figures, and physics framing. |
| PPG12 legacy box-cut note | `usefulDocs/PPG12_analysis_note_2026-01-07_legacy_boxcut.pdf` | Historical reference for photon-ID box cuts and isolated-photon pp baseline when older context is needed. |
| PPG12 legacy BDT note | `usefulDocs/PPG12_analysis_note_2026-05-03_legacy_withBDT.pdf` | Historical BDT/NPB reference; use only when comparing against older slide/code discussions or when v4 does not contain the needed detail. |
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

## Canonical Embedded Stitching Rules

Use these rules for future corrected embedded RecoilJets products unless
Justin explicitly starts a new Jet40 or alternate-boundary study.

- `isSimEmbedded` PhotonJet12+20: PhotonJet12 owns `12-21` and PhotonJet20
  owns `>=21`. Use matching effective cross sections for exactly those
  ownership windows: PhotonJet12 `2663.51030 pb`, PhotonJet20 `92.0516019 pb`.
- `isSimEmbeddedInclusive` Jet12+20+30: Jet12 owns `12-21`, Jet20 owns
  `21-31`, and Jet30 owns `>=31`. Use matching effective cross sections:
  Jet12 `1.22941038e6 pb`, Jet20 `3.87846122e4 pb`, Jet30 `1.79304310e3 pb`.
- Lower-only variants are diagnostics only because the samples overlap above
  the threshold and cannot be directly summed as a final stitched spectrum.
- Pre-fix Jet12+20+30 embedded-inclusive RecoilJets backgrounds made with the
  old worker gates are stale/tainted for stitched-background plots. Historical
  training caches are not invalidated by this specific RecoilJets worker-gate
  bug.
- Jet40 is not part of the current validated embedded-inclusive production
  because the embedded Jet40 sample is not available. Revisit only if a Jet40
  sample is added.
- Fine reco-cluster `E_T` leakage artifacts from the 12-40 GeV diagnostic
  should produce two final PNGs: a full-slide candidate to replace slide 7 and
  a plot-only PNG for Blair. Put the `#it{#bf{sPHENIX}} Internal` and PYTHIA8
  labels inside the plot canvas/axes so screenshots of the plot retain the
  internal-analysis label; do not leave those labels detached in the slide
  header.
- Completed fine reco-cluster `E_T` diagnostic tag
  `clusterEtFine40_stitch_20260521_212953` produced the canonical 1 GeV
  `12-40 GeV` artifacts:
  `dataOutput/stitchDiagnostics/focus21_clusterEt_leakage_20260521/inclusive_jet123_reco_cluster_et_weighted_components_blair.png`
  and
  `dataOutput/stitchDiagnostics/focus21_clusterEt_leakage_20260521/inclusive_jet123_reco_cluster_et_leakage_slide12_followup.png`.
  The corrected weighted diagnostic has no Jet12 contribution above 22 GeV;
  the 35-40 GeV bins are present but have zero total entries in this run.

## Scaled Trigger Study Pipeline

Use this entry when Justin refers to the scaled-trigger analysis that produced
slide 3 of `WP_GammaJets_5_20_26`
(`1-2l9IWNJMpFfunoShVrYk3ywtrIocMWMyePYnR85_zc`, slide
`g3e29c967769_0_16`).

- Run-list/scaler source: `./scripts/make_dstListsData.sh auau QA
  scaledTriggerAna` builds the common AuAu GRL subset
  `MBD_NS_geq_2_vtx_lt_150__Pho10_12`, requiring scaled-positive MBD
  vtx<150, Photon 10, and Photon 12 in the same run. The selected sample was
  620 runs.
- Production source: `scripts/RecoilJets_Condor_submit.sh` action
  `isAuAu scaledTriggerStudy` prepares a single config tagged
  `jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionReference_tightReference_nonTightReference_scaledTriggerStudy`
  and runs with `RJ_SCALED_TRIGGER_STUDY_ONLY=1`.
- Event loop source: `macros/Fun4All_recoilJets_unified_impl.C`
  `ScaledTriggerStudyReco` keeps geometry, calibration/clustering, MBD reco,
  and vertex reco, then fills max-cluster-energy histograms directly from
  `ScaledVector` bits 14, 22, and 23 for MBD vtx<150, Photon 10, and
  Photon 12.
- Merge source: `scripts/mergeRecoilJets.sh` applies per-run live/scaled
  factors from the scaled-trigger config before combining runs, marking ROOT
  chunks with `scaledTrigQA_perRunCorrected_applied`.
- Local artifacts: pulled ROOT
  `InputFiles/auau25/RecoilJets_auau_ALL_jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionReference_tightReference_nonTightReference_scaledTriggerStudy.root`;
  slide PNGs in `dataOutput/auau/scaledTriggerStudy/`; plotting macro
  `macros/PlotScaledTriggerStudy.C`.
- Slide-3 result: left plot is live/scaled-corrected max-cluster spectra for
  MBD, Photon 10, and Photon 12; right plot is Photon/MBD turn-on. Current
  plotted values are Photon 10 plateau about 0.911 with 90% turn-on near
  7.43 GeV, and Photon 12 plateau about 0.892 with 90% turn-on near 8.59 GeV.
  The unresolved caveat is the high-energy plateau below unity, to be checked
  versus centrality and run-by-run before treating it as fully understood.

## Decision Rules

- For photon-ID cuts, start with the current PPG12 IAN and paper draft before
  Au+Au changes.
- For BDT/NPB cuts, start with the current PPG12 IAN, paper draft, and
  `FunWithxgboost`.
- For pp baseline comparisons, preserve PPG18/PPG12 behavior unless there is a
  documented reason to diverge.
- For Au+Au photon ID, treat new BDT/NPB/JetML variants as extensions that
  require embedded validation and pp comparison where practical.
- For xJgamma, start from the existing RecoilJets pp pipeline, then extend to
  Au+Au with embedded response/closure checks.
- For final presentation, make outputs ATLAS-like in analysis logic while
  retaining sPHENIX style and PPG constraints.
