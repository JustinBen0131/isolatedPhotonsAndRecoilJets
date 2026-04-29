# Analysis Environment

For any ROOT macro, `root`, `root-config`, ACLiC build, or ROOT-dependent script in this repo, use the analysis environment that matches the user's `use_root_analysis` fish function.

Canonical wrapper:

```bash
./scripts/root_in_analysis_env.sh /Users/patsfan753/Desktop/analysis/env/bin/root -l -q 'macros/MyMacro.C()'
```

Notes:

- This wrapper reproduces the `use_root_analysis` setup from `~/.config/fish/config.fish`.
- It forces the known-good RooUnfold install at `external/RooUnfold` to win over partial env copies.
- Do not rely on interactive shell startup files for ROOT in this repo; use the wrapper explicitly.

# Change Control

- Before making any file edit, state exactly which file(s) will change and what the change is intended to do.
- Do not make the edit until the user explicitly approves it, unless the user has already directly requested that exact edit.

# sPHENIX Pipeline Paths

Changes to any of the following files are changes the user will copy into the
SSH analysis checkout and run from:

```bash
/sphenix/u/patsfan753/scratch/thesisAnalysis
```

Treat these as sPHENIX-side pipeline files, not local Mac-only files:

- `/Users/patsfan753/Desktop/ThesisAnalysis/scripts/mergeRecoilJets.sh`
- `/Users/patsfan753/Desktop/ThesisAnalysis/scripts/RecoilJets_Condor_AuAu.sh`
- `/Users/patsfan753/Desktop/ThesisAnalysis/scripts/RecoilJets_Condor_submit.sh`
- `/Users/patsfan753/Desktop/ThesisAnalysis/scripts/RecoilJets_Condor.sh`
- `/Users/patsfan753/Desktop/ThesisAnalysis/src_AuAu/RecoilJets_AuAu.h`
- `/Users/patsfan753/Desktop/ThesisAnalysis/src_AuAu/RecoilJets_AuAu.cc`
- `/Users/patsfan753/Desktop/ThesisAnalysis/src/PhotonClusterBuilder.cc`
- `/Users/patsfan753/Desktop/ThesisAnalysis/src/PhotonClusterBuilder.h`
- `/Users/patsfan753/Desktop/ThesisAnalysis/src/RecoilJets.cc`
- `/Users/patsfan753/Desktop/ThesisAnalysis/src/RecoilJets.h`
- `/Users/patsfan753/Desktop/ThesisAnalysis/macros/analysis_config.yaml`
- `/Users/patsfan753/Desktop/ThesisAnalysis/macros/Fun4All_recoilJets_AuAu.C`
- `/Users/patsfan753/Desktop/ThesisAnalysis/macros/Fun4All_recoilJets_unified_impl.C`
- `/Users/patsfan753/Desktop/ThesisAnalysis/macros/Fun4All_recoilJets.C`

Only offline macros that the user runs after the online pipeline completes
should be treated as local-only code.
