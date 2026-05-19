# Current Project State

This file captures durable current-state caveats that future Codex should check
before dataset-validity claims.

## ABCD / Purity Tight-Axis Fix

Local git commit `507eabd3` (`2026-05-05 21:17 EDT`) updated:

- `src/RecoilJets.cc`
- `src_AuAu/RecoilJets_AuAu.cc`

The fix changes `fillIsoSSTagCounters(...)` to use the caller-computed
`TightTag` from the active configured photon-ID variants.

Do not consider affected output families production-valid until the source is
transferred to SDCC, affected jobs rerun, fresh ROOT outputs pulled, and
`codex_notes/DATASET_STATUS.md` updated with evidence.

## Affected Families

Treat pre-fix ROOT products using non-reference tight or non-tight photon-ID
axes as stale only for affected ABCD/purity-derived histogram families:

- candidate ABCD counts:
  `h_isIsolated_isTight*`, `h_notIsolated_isTight*`,
  `h_isIsolated_notTight*`, `h_notIsolated_notTight*`
- candidate ABCD distributions:
  `h_Eiso_ABCD_[A-D]*`, `h_pTgamma_ABCD_[A-D]*`
- ABCD-tagged shower-shape templates:
  `h_ss_<var>_isIsolated_isTight*`,
  `h_ss_<var>_notIsolated_isTight*`,
  `h_ss_<var>_isIsolated_notTight*`,
  `h_ss_<var>_notIsolated_notTight*`
- truth-signal leakage:
  `h_sigABCD_MC*`
- event-leading xJ purity counts:
  `h_xJpurityLead_isIsolated_isTight*`,
  `h_xJpurityLead_notIsolated_isTight*`,
  `h_xJpurityLead_isIsolated_notTight*`,
  `h_xJpurityLead_notIsolated_notTight*`
- event-leading xJ leakage:
  `h_xJpurityLead_sigABCD_MC*`

Non-ABCD histograms, trees, event QA, jet QA, trigger QA, and unrelated outputs
can still be used if the requested plot/result does not consume affected
families. Check the macro/input before warning broadly.

`tightReference/nonTightReference` outputs are lower risk for this specific bug
because the old hard-coded ABCD axis matched reference behavior, but still need
normal provenance checks.

When fresh rerun outputs are pulled into `InputFiles/`, update
`codex_notes/DATASET_STATUS.md` and retire obsolete stale warnings.
