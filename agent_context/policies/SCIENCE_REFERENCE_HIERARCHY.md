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
- PPG12 box-cut note: `usefulDocs/PPG12_analysis_note.pdf`
- PPG12 BDT/NPB note: `usefulDocs/PPG12_analysis_note_withBDT.pdf`
- ATLAS gamma-jet target: `usefulDocs/Gamma_Jet_Analysis_Note (1).pdf`
- Au+Au xJ/dijet correction/style: `usefulDocs/PPG_08_dijet_xJ_in_Au_Au___draft_Conference_note (7).pdf`
- Trigger semantics: `usefulDocs/Gl1-gtm_user_manual_v53.pdf`
- Compact map: `agent_context/REFERENCE_MAP.md`

## Decision Rules

- Do not invent photon-ID, purity, unfolding, xJ, BDT/NPB, stitching, or
  systematic conventions when a PPG12/PPG18/ATLAS/PPG08 precedent exists.
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
