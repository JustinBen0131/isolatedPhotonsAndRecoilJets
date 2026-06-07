# Science Reference Hierarchy

## North Star

- Ultimate target: PPG19 Au+Au gamma-jet / `x_{J#gamma}`.
- PPG18 pp gamma-jet is the validated reference baseline for Au+Au extensions.
- PPG12 is the photon-ID, isolation, BDT/NPB, purity, and pp infrastructure
  reference.
- ATLAS gamma-jet is the final-analysis note style and correction target,
  adapted to sPHENIX/PPG constraints.

## Local References

- PPG12 implementation reference: `ppg12codeGit/`
- PPG12 current IAN:
  `usefulDocs/PPG12_analysis_note_2026-05-21_v4_current_IAN.pdf`
- PPG12 current paper draft:
  `usefulDocs/sPHENIX_PPG12_Paper_2026-05-21_current_draft.pdf`
- PPG12 legacy box-cut note:
  `usefulDocs/PPG12_analysis_note_2026-01-07_legacy_boxcut.pdf`
- PPG12 legacy BDT/NPB note:
  `usefulDocs/PPG12_analysis_note_2026-05-03_legacy_withBDT.pdf`
- ATLAS gamma-jet target: `usefulDocs/Gamma_Jet_Analysis_Note (1).pdf`
- Au+Au xJ/dijet correction/style: `usefulDocs/PPG_08_dijet_xJ_in_Au_Au___draft_Conference_note (7).pdf`
- Trigger semantics: `usefulDocs/Gl1-gtm_user_manual_v53.pdf`
- Compact map: `agent_context/REFERENCE_MAP.md`

## Decision Rules

- Do not invent photon-ID, purity, unfolding, xJ, BDT/NPB, stitching, or
  systematic conventions when a PPG12/PPG18/ATLAS/PPG08 precedent exists.
- For pp inclusive-jet stitching and stitched-source diagnostics, the baseline
  contract is canonical candidate rows, `ppg12_truth_window_pass_r04 > 0.5`,
  PPG12 sample/event cross-section weights, and recorded sample ownership caps.
  Keep trees broad enough for diagnostics, but require downstream consumers to
  apply the truth-window and ownership-cap gates before claiming a coherent
  stitched pp baseline. For ABCD-facing reco-cluster `E_T` diagnostics, also
  match the reference ABCD row scope: fixed-iso signal/sideband rows, reference
  tight or non-tight candidates, and exclusion of isolation-gap/neither rows.
- For PPG12 photon-ID, BDT/NPB, purity, isolation, pp baseline, and paper
  wording, prefer
  `PPG12_analysis_note_2026-05-21_v4_current_IAN.pdf` and
  `sPHENIX_PPG12_Paper_2026-05-21_current_draft.pdf` over older local PPG12
  PDFs unless the task is explicitly about historical comparisons.
- Reasoning order for physics changes:
  1. validate pp baseline behavior;
  2. extend to AuAu/embedded context;
  3. validate with embedded signal and embedded inclusive/background;
  4. promote to final products.
- For corrections, purity, ABCD, unfolding, normalization, stitching, response
  matrices, and surprising outputs, enter skeptic mode: check math, stale
  provenance, histogram-family mismatches, data/MC definitions, and pipeline
  stage.
- For trigger logic, GL1 bits, `ScaledVector`, live/scaled counts, or scaledown
  interpretation, check `usefulDocs/Gl1-gtm_user_manual_v53.pdf` before making
  trigger-semantics claims or code changes.
- For routine operational/code fixes, stay surgical and follow existing
  architecture.

## Current AuAu Photon-ID Baseline Bias

For the current AuAu/JSTG BDT story, default wording to the PPG12 photon-ID
25-feature family plus one AuAu context input, `centrality`, when that is the
chosen model. Treat older 32-feature baselines as controls/comparisons unless
Justin asks otherwise.

Recent known evidence: the Jet12+20+30 sixpack ablation had
`globalEtCent1535_bdt_ppg12PlusCent` READY with AUC `0.841561`, slightly above
the 32-feature baseline while being simpler and closer to PPG12.
