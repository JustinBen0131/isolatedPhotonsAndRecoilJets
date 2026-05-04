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

# Input Files

- `InputFiles/` contains local analysis input ROOT files and should be treated as a first-class project directory like `macros/`, `scripts/`, `src/`, and `src_AuAu/`.
- When the user refers to input files or local ROOT inputs, first check `InputFiles/` from the repository base, especially `InputFiles/simEmbedded/`, `InputFiles/simPhotonJet/`, and data-specific subdirectories.
- The `@InputFiles` UI mention may not autocomplete in all clients; still use the literal path `InputFiles/...` when the user refers to those files.

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

# SDCC Transfer Workflow

Use `scripts/sftp_push_recoiljets.sh` as the only supported transfer helper for
copying known SSH pipeline files from the local Mac checkout to the SDCC
analysis checkout.

Use `scripts/sftp_get_recoiljets_outputs.sh` as the local-only helper for
pulling merged ROOT outputs from SDCC into `InputFiles/...`. It opens
interactive `sftp`, lists exact remote final ROOT files for the requested
dataset, prints the local overwrite preview, and then downloads only those
listed files after user confirmation. Do not SSH directly, do not run raw
`sftp` manually from Codex, and do not ask for or handle the SDCC password.

Local repo base:

```bash
/Users/patsfan753/Desktop/ThesisAnalysis
```

Remote SDCC repo base:

```bash
/sphenix/u/patsfan753/scratch/thesisAnalysis
```

Remote host:

```bash
patsfan753@sftp.sdcc.bnl.gov
```

Rules for future Codex runs:

- Treat the mapped files below as SSH pipeline files.
- When Codex edits any mapped file, its final response must include the exact
  `./scripts/sftp_push_recoiljets.sh ...` command the user should run.
- When the user wants the local GitHub repo updated after a successful SDCC
  transfer, include the opt-in commit form:

```bash
./scripts/sftp_push_recoiljets.sh <changed-files-or-groups> --commit-push -m "Concise change summary"
```

- The commit message should describe the actual code/config change, not just say
  "changes".
- The uploader commits only the selected transferable files after successful
  SFTP. It must never use `git add .` for this workflow.
- Do not SSH directly.
- Do not run `sftp` directly.
- Do not ask for or handle the SDCC password.
- The user runs the uploader manually and enters their password.
- Prefer the smallest upload command that covers the changed files.
- Prefer explicit changed file names for source files. Use `RecoilJets.cc`
  and `RecoilJets.h` directly for `src/`; use `RecoilJets_AuAu.cc` and
  `RecoilJets_AuAu.h` directly for `src_AuAu/`. Do not use `pp`, `auau`,
  or `both` source groups.
- If one known file changed, give a one-file command, for example:

```bash
./scripts/sftp_push_recoiljets.sh mergeRecoilJets.sh
```

- If several known files changed, give one command with all file names, for
  example:

```bash
./scripts/sftp_push_recoiljets.sh mergeRecoilJets.sh analysis_config.yaml RecoilJets_AuAu.cc
```

- If a whole non-source logical group changed, use the group, for example:

```bash
./scripts/sftp_push_recoiljets.sh condor
```

- If Codex creates or edits a new file that is needed on SDCC but is not in the
  uploader allowlist, explicitly say that the uploader allowlist needs to be
  extended. Ask before editing the uploader/AGENTS.md unless the user already
  requested that exact update.
- Do not use wildcard uploads for `macros/`, `scripts/`, `src/`, or `src_AuAu/`.
- New SSH-needed macro files must be added explicitly to the uploader allowlist
  with their intended remote path.
- `scripts/sftp_push_recoiljets.sh` is local-only. Do not include it in SDCC
  transfer commands.
- `scripts/sftp_get_recoiljets_outputs.sh` is local-only. Do not include it in
  SDCC transfer commands.
- Because the uploader/getter are local-only, do not rely on `--commit-push` to
  commit edits to those helper scripts or to `AGENTS.md`. Use normal local git
  commands for those files if needed.

Known transfer map:

```text
scripts/audit_auau_grl_projection.sh      -> scripts/audit_auau_grl_projection.sh
scripts/estimateEmbeddedPhotonXsec.sh     -> scripts/estimateEmbeddedPhotonXsec.sh
scripts/make_dstListsData.sh              -> scripts/make_dstListsData.sh
scripts/makeThesisSimLists.sh             -> scripts/makeThesisSimLists.sh
scripts/mergeRecoilJets.sh                -> scripts/mergeRecoilJets.sh
scripts/root_in_analysis_env.sh           -> scripts/root_in_analysis_env.sh
scripts/RecoilJets_Condor_AuAu.sh         -> RecoilJets_Condor_AuAu.sh
scripts/RecoilJets_Condor_submit.sh       -> RecoilJets_Condor_submit.sh
scripts/RecoilJets_Condor.sh              -> RecoilJets_Condor.sh
macros/analysis_config.yaml               -> macros/analysis_config.yaml
macros/Fun4All_recoilJets.C               -> macros/Fun4All_recoilJets.C
macros/Fun4All_recoilJets_AuAu.C          -> macros/Fun4All_recoilJets_AuAu.C
macros/Fun4All_recoilJets_unified_impl.C  -> macros/Fun4All_recoilJets_unified_impl.C
macros/PrintPPStitchDiagnostics.C         -> macros/PrintPPStitchDiagnostics.C
src/RecoilJets.cc                         -> src/RecoilJets.cc
src/RecoilJets.h                          -> src/RecoilJets.h
src_AuAu/RecoilJets_AuAu.cc               -> src_AuAu/RecoilJets_AuAu.cc
src_AuAu/RecoilJets_AuAu.h                -> src_AuAu/RecoilJets_AuAu.h
```

Uploader groups:

```text
condor    scripts/RecoilJets_Condor.sh + scripts/RecoilJets_Condor_AuAu.sh + scripts/RecoilJets_Condor_submit.sh
macros    known pipeline macros/configs
scripts   known pipeline helper scripts
pipeline  all known transferable pipeline files except local-only sftp helpers
```

Important mapping detail: most files preserve their relative path on the remote
side, but these three local files are under `scripts/` locally and live at the
remote repo base on SDCC:

```text
scripts/RecoilJets_Condor.sh       -> RecoilJets_Condor.sh
scripts/RecoilJets_Condor_AuAu.sh  -> RecoilJets_Condor_AuAu.sh
scripts/RecoilJets_Condor_submit.sh -> RecoilJets_Condor_submit.sh
```
