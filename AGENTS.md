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
