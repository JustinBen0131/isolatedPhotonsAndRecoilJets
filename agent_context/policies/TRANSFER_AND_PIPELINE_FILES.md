# Transfer And Pipeline Files

## Helpers

Use `scripts/sftp_push_recoiljets.sh` for mapped local-to-SDCC pipeline uploads
and read-only status/diff checks.

Use `scripts/sftp_get_recoiljets_outputs.sh` for pulling merged ROOT outputs
from SDCC into `InputFiles/...` and for local-only SIM merge helpers.

Do not SSH directly for transfer, do not run raw `sftp` manually from Codex,
and do not ask for or handle the SDCC password.

Before each helper run, state:

- exact command;
- action class: read-only status/diff, upload, download, local-only merge, or
  commit/push after upload;
- expected local and SDCC paths/files.

Routine low-risk helper runs can proceed after that preflight notice when the
task clearly calls for them. Ask for explicit confirmation before broad/risky
actions: groups like `pipeline`, `all`, or `changed`; ambiguous or overwrite
prone downloads; substantive `--commit-push`; allowlist extensions; or anything
that could delete/submit/hold/release/mutate beyond selected copied files.

Use the smallest explicit file list that does the job. Avoid wildcard uploads
and avoid broad groups unless the change genuinely spans that group. If a new
SDCC-needed macro/script/config/source file is not in the helper allowlist,
pause and ask before extending the allowlist.

If an SFTP helper opens an interactive password prompt, Justin enters it
manually. Codex must not ask for, store, or type the password.

## Read-Only Checks

```bash
./scripts/sftp_push_recoiljets.sh status
./scripts/sftp_push_recoiljets.sh status changed
./scripts/sftp_push_recoiljets.sh status condor
./scripts/sftp_push_recoiljets.sh diff <file-or-basename>
./scripts/sftp_push_recoiljets.sh diff changed
```

Interpretation:

- `MATCH`: local and SDCC byte-identical.
- `DIFFER`: both exist but differ.
- `MISSING`: absent on SDCC at mapped path.

## Upload Final-Response Rule

When Codex edits a mapped SDCC-side file and has not already run the upload
helper, final response must include the exact smallest upload command. For
substantive SDCC pipeline edits, include selected-file `--commit-push` form by
default unless Justin asks upload-only:

```bash
./scripts/sftp_push_recoiljets.sh <files-or-groups> --commit-push -m "Concise change summary"
```

The uploader commits only selected transferable files. It must never use
`git add .`.

For clean requested edits to mapped SDCC pipeline files, Codex should offer to
run the exact selected-file uploader command after tests pass, not merely leave
it for Justin. Ask first and explain which local paths will be uploaded, which
SDCC paths will be overwritten, and whether `makeProject` or another rebuild is
required afterward.

If Justin explicitly asks for upload-only, omit `--commit-push`. Commit
messages should describe the actual code/config change, not just say "changes".
Do not rely on the uploader's `--commit-push` to commit edits to the local-only
SFTP helpers or to `AGENTS.md`; use normal local git workflow for those files
if needed.

## Local And Remote Bases

Local repo:

```bash
/Users/patsfan753/Desktop/ThesisAnalysis
```

Remote SDCC repo:

```bash
/sphenix/u/patsfan753/scratch/thesisAnalysis
```

SFTP host:

```bash
patsfan753@sftp.sdcc.bnl.gov
```

## SDCC-Side Pipeline Files

Treat these as sPHENIX-side pipeline files Justin will copy/run from SDCC, not
local Mac-only files:

- `scripts/mergeRecoilJets.sh`
- `scripts/recoiljets_cleanup.sh`
- `scripts/RecoilJets_Condor_AuAu.sh`
- `scripts/RecoilJets_Condor_submit.sh`
- `scripts/RecoilJets_Condor.sh`
- `src_AuAu/RecoilJets_AuAu.h`
- `src_AuAu/RecoilJets_AuAu.cc`
- `src/PhotonClusterBuilder.cc`
- `src/PhotonClusterBuilder.h`
- `src/RecoilJets.cc`
- `src/RecoilJets.h`
- `macros/analysis_config.yaml`
- `macros/Fun4All_recoilJets_AuAu.C`
- `macros/Fun4All_recoilJets_unified_impl.C`
- `macros/Fun4All_recoilJets.C`

`macros/AnalyzeRecoilJets.h` is offline only. Do not make remote
production/merge/final-stitch scripts depend on it.

Only offline macros that Justin runs after the online pipeline completes should
be treated as local-only code.

If compiled source/headers change, tell Justin the correct rebuild directory:

```bash
cd /sphenix/u/patsfan753/scratch/thesisAnalysis/src_AuAu
make clean
makeProject

cd /sphenix/u/patsfan753/scratch/thesisAnalysis/src
make clean
makeProject
```

For shell/YAML/macro-only changes that do not require relinking, say no
`makeProject` rebuild is needed.

## Known Transfer Map

```text
scripts/audit_auau_grl_projection.sh      -> scripts/audit_auau_grl_projection.sh
scripts/audit_auau_ml_training_smoke.py   -> scripts/audit_auau_ml_training_smoke.py
scripts/estimateEmbeddedPhotonXsec.sh     -> scripts/estimateEmbeddedPhotonXsec.sh
scripts/make_dstListsData.sh              -> scripts/make_dstListsData.sh
scripts/makeThesisSimLists.sh             -> scripts/makeThesisSimLists.sh
scripts/mergeRecoilJets.sh                -> scripts/mergeRecoilJets.sh
scripts/recoiljets_cleanup.sh             -> scripts/recoiljets_cleanup.sh
scripts/root_in_analysis_env.sh           -> scripts/root_in_analysis_env.sh
scripts/submit_auau_bdt_widthstudy_pt1530_wp080.sh -> scripts/submit_auau_bdt_widthstudy_pt1530_wp080.sh
scripts/RecoilJets_Condor_AuAu.sh         -> RecoilJets_Condor_AuAu.sh
scripts/RecoilJets_Condor_submit.sh       -> RecoilJets_Condor_submit.sh
scripts/RecoilJets_Condor.sh              -> RecoilJets_Condor.sh
scripts/train_auau_jet_residual_bdt.py    -> scripts/train_auau_jet_residual_bdt.py
scripts/train_auau_photon_bdt.py          -> scripts/train_auau_photon_bdt.py
scripts/validate_auau_tight_bdt_on_sim.py -> scripts/validate_auau_tight_bdt_on_sim.py
macros/analysis_config.yaml               -> macros/analysis_config.yaml
macros/analysis_config_auau_bdt_widthstudy_pt1530_wp080.yaml -> macros/analysis_config_auau_bdt_widthstudy_pt1530_wp080.yaml
macros/Fun4All_recoilJets.C               -> macros/Fun4All_recoilJets.C
macros/Fun4All_recoilJets_AuAu.C          -> macros/Fun4All_recoilJets_AuAu.C
macros/Fun4All_recoilJets_unified_impl.C  -> macros/Fun4All_recoilJets_unified_impl.C
macros/PrintPPStitchDiagnostics.C         -> macros/PrintPPStitchDiagnostics.C
src/RecoilJets.cc                         -> src/RecoilJets.cc
src/RecoilJets.h                          -> src/RecoilJets.h
src_AuAu/RecoilJets_AuAu.cc               -> src_AuAu/RecoilJets_AuAu.cc
src_AuAu/RecoilJets_AuAu.h                -> src_AuAu/RecoilJets_AuAu.h
```

Groups:

```text
condor    scripts/RecoilJets_Condor.sh + scripts/RecoilJets_Condor_AuAu.sh + scripts/RecoilJets_Condor_submit.sh
macros    known pipeline macros/configs
scripts   known pipeline helper scripts
pipeline  all known transferable pipeline files except local-only sftp helpers
```

Important mapping detail: `scripts/RecoilJets_Condor*.sh` live under
`scripts/` locally but at the remote repo base on SDCC.

The transfer helpers themselves are local-only. Do not include them in SDCC
transfer commands, and do not rely on helper `--commit-push` for edits to the
helpers or `AGENTS.md`.
