# Known Issues

Last updated: 2026-05-05

## ABCD/Purity Tight-Axis Mismatch

Status: Source fix present in local git commit `507eabd3`; SDCC transfer/rerun evidence not recorded.

Evidence:

- Local commit `507eabd3` changes `src/RecoilJets.cc` and `src_AuAu/RecoilJets_AuAu.cc` so `fillIsoSSTagCounters(...)` receives the caller-computed `TightTag`.
- The caller-computed `TightTag` uses the active configured tight/non-tight photon-ID variant, including `tightVariantA` / `newPPG12`-style behavior.
- `git blame` now attributes the key call-site fix lines to `507eabd3` at `2026-05-05 21:17:40 -0400`.

Impact:

- A/B/C/D region counts may be wrong for affected non-reference tight configurations produced before the fix.
- Purity histograms, leakage/correction factors, and downstream corrected xJ results derived from the affected ABCD families may be wrong.
- The highest-risk outputs are those using non-reference tight or non-tight photon-ID axes, such as `tightVariantA`, `nonTightVariantA`, renamed/new `newPPG12`-style settings, or later AuAu BDT sideband settings.
- `tightReference/nonTightReference` outputs are lower risk for this specific bug because the old hard-coded ABCD tight axis matched reference behavior.
- Outputs that only changed preselection while keeping `tightReference/nonTightReference` are not automatically stale for this tight-axis mismatch.
- The ROOT files are not globally invalid. Histograms, trees, QA, and plots that do not read or derive from the affected families can still be used for this specific issue.

Affected histogram families:

- Candidate ABCD counts: `h_isIsolated_isTight*`, `h_notIsolated_isTight*`, `h_isIsolated_notTight*`, `h_notIsolated_notTight*`.
- Candidate ABCD distributions: `h_Eiso_ABCD_[A-D]*`, `h_pTgamma_ABCD_[A-D]*`.
- ABCD-tagged shower-shape templates: `h_ss_<var>_isIsolated_isTight*`, `h_ss_<var>_notIsolated_isTight*`, `h_ss_<var>_isIsolated_notTight*`, `h_ss_<var>_notIsolated_notTight*`.
- Truth-signal ABCD leakage: `h_sigABCD_MC*`.
- Event-leading xJ purity counts: `h_xJpurityLead_isIsolated_isTight*`, `h_xJpurityLead_notIsolated_isTight*`, `h_xJpurityLead_isIsolated_notTight*`, `h_xJpurityLead_notIsolated_notTight*`.
- Event-leading xJ leakage: `h_xJpurityLead_sigABCD_MC*`.

Required closure criteria:

- Source fix transferred to SDCC from local commit `507eabd3` or a later commit containing the same fix.
- `RecoilJets.cc` and `RecoilJets_AuAu.cc` transferred to SDCC.
- Affected jobs rerun after transfer.
- Fresh ROOT outputs pulled into `InputFiles/`.
- `codex_notes/DATASET_STATUS.md` updated with timestamps, job IDs or SDCC pasted evidence, and final status.
- Obsolete stale warnings removed once post-fix output evidence exists.

Current recommended transfer command when ready:

```bash
./scripts/sftp_push_recoiljets.sh RecoilJets.cc RecoilJets_AuAu.cc
```
