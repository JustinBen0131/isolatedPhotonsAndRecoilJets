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

Never transfer local Codex control-plane state to SDCC. `agent_context/`,
`.codex/`, `codex*` migration/archive files, and local-only OS evidence must
remain local. SDCC transfer or remote-migration helpers must also reject empty,
whitespace-only, or leading/trailing-space path strings before creating
directories or moving files.

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

## Remote Submit Provenance

When a transferred or paste-ready SDCC command will launch
`RecoilJets_Condor_submit.sh`, `mergeRecoilJets.sh`, or another RecoilJets
workflow that emits pipeline emails, include Codex provenance in the remote
environment:

```bash
RJ_CODEX_CHAT_NAME="<short-workstream-label>"
RJ_CODEX_THREAD_ID="<codex-thread-id>"
export RJ_CODEX_CHAT_NAME RJ_CODEX_THREAD_ID
```

These assignments must be evaluated on SDCC before the DAG/meta files are
created. Passing only local `CODEX_*` variables on macOS is not enough, because
nested SSH and Condor scheduler nodes do not inherit them automatically.
If the thread id is unknown in the current shell, inspect the local environment
for `CODEX_THREAD_ID` before composing the remote command; otherwise stop and
record why the metadata could not be attached.

## Local And Remote Bases

Local repo:

```bash
/Users/patsfan753/Desktop/ThesisAnalysis
```

Remote SDCC repo and SFTP host are local-only access details. Load
`agent_context/local/SDCC_LOCAL_RULES.md` in Justin's local checkout or ask
Justin for the current local rule. Do not commit hostnames, usernames, or
remote checkout paths to GitHub-tracked policy files.

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

These SDCC-side pipeline files are protected path anchors for repo organization
purposes. Do not move or rename them as part of scripts/macros cleanup unless a
specific refactor is approved with a compatibility plan, reference audit, SDCC
upload plan, and smoke tests. New local side helpers should be organized around
them, not by moving these anchors opportunistically.

If compiled source/headers change, tell Justin the correct rebuild directory:

```bash
cd <REMOTE_THESIS_ANALYSIS>/src_AuAu
make clean
makeProject

cd <REMOTE_THESIS_ANALYSIS>/src
make clean
makeProject
```

For shell/YAML/macro-only changes that do not require relinking, say no
`makeProject` rebuild is needed.

## Known Transfer Map

```text
Exhaustive script mappings live in `scripts/sdcc/TRANSFER_MAP.tsv`,
`scripts/sdcc/SDCC_RUNTIME_INDEX.yaml`, and the generated
`LOCAL_FILES`/`REMOTE_FILES` arrays in
`scripts/sdcc/transfer/sftp_push_recoiljets.sh`.

Key script examples after THE-23 stage 5:

scripts/compat/local/audit_auau_grl_projection.sh -> scripts/sdcc/runtime/audit/audit_auau_grl_projection.sh
scripts/compat/local/audit_auau_ml_training_smoke.py -> scripts/ml/audits/audit_auau_ml_training_smoke.py
scripts/estimateEmbeddedPhotonXsec.sh     -> scripts/sdcc/runtime/xsec/estimateEmbeddedPhotonXsec.sh
scripts/make_dstListsData.sh              -> scripts/sdcc/runtime/lists/make_dstListsData.sh
scripts/makeThesisSimLists.sh             -> scripts/sdcc/runtime/lists/makeThesisSimLists.sh
scripts/mergeRecoilJets.sh                -> scripts/sdcc/runtime/merge/mergeRecoilJets.sh
scripts/recoiljets_cleanup.sh             -> scripts/sdcc/runtime/cleanup/recoiljets_cleanup.sh
scripts/root_in_analysis_env.sh           -> scripts/env/root_in_analysis_env.sh
scripts/compat/local/submit_auau_bdt_widthstudy_pt1530_wp080.sh -> scripts/sdcc/workflows/width_study/submit_auau_bdt_widthstudy_pt1530_wp080.sh
scripts/RecoilJets_Condor_AuAu.sh         -> RecoilJets_Condor_AuAu.sh
scripts/RecoilJets_Condor_submit.sh       -> RecoilJets_Condor_submit.sh
scripts/RecoilJets_Condor.sh              -> RecoilJets_Condor.sh
scripts/compat/local/train_auau_jet_residual_bdt.py -> scripts/ml/stacking/train_auau_jet_residual_bdt.py
scripts/compat/local/train_auau_photon_bdt.py -> scripts/ml/training/train_auau_photon_bdt.py
scripts/compat/local/validate_auau_tight_bdt_on_sim.py -> scripts/ml/validation/validate_auau_tight_bdt_on_sim.py
macros/analysis_config.yaml               -> macros/analysis_config.yaml
macros/analysis_config_auau_bdt_widthstudy_pt1530_wp080.yaml -> macros/analysis_config_auau_bdt_widthstudy_pt1530_wp080.yaml
macros/Fun4All_recoilJets.C               -> macros/Fun4All_recoilJets.C
macros/Fun4All_recoilJets_AuAu.C          -> macros/Fun4All_recoilJets_AuAu.C
macros/Fun4All_recoilJets_unified_impl.C  -> macros/Fun4All_recoilJets_unified_impl.C
macros/PrintPPStitchDiagnostics.C         -> macros/diagnostics/stitching/PrintPPStitchDiagnostics.C
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

Generated campaign configs should not be added to the static transfer allowlist
unless they exist locally and need an intentional upload. For example,
`analysis_config_auau_bdt_mlp_stack_wp080.yaml` is produced by the stack driver
from `analysis_config_auau_bdt_mlp_stack_template.yaml`; upload the template or
the driver, not a missing generated output path.

## Local Organization Detail

After THE-23 stage 5, local and SDCC top-level `scripts/` paths are minimal
command surfaces. Local canonical source lives under purpose folders such as
`scripts/sdcc/runtime/{condor,merge,lists,cleanup,xsec,audit}/`,
`scripts/sdcc/workflows/{submit,target_wp,width_study,stacking,scan,diagnostics}/`,
`scripts/sdcc/pipelines/{auau,pp}/`, `scripts/ml/{training,validation,working_points,stacking,audits}/`,
`scripts/plotting/{auau_bdt,pp_currentian,stitching,efficiency,trigger,truth_purity}/`,
`scripts/diagnostics/`, `scripts/data_prep/`, `scripts/slides/`, and
`scripts/os/`. On SDCC, canonical script source lives in matching subfolders;
historical remote `scripts/foo` paths are symlinks for operator/command
compatibility.

For SDCC uploads, continue using the historical local command path or group
name with `scripts/sftp_push_recoiljets.sh`. The helper resolves local symlink
targets before upload and compares/stages the canonical source. For scripts,
the remote write path is now usually the canonical SDCC subfolder path, while
the old remote `scripts/foo` symlink remains executable. Do not raw-copy
canonical folders to SDCC.

After THE-23 stage 6, future SDCC output/evidence paths should also be
canonical-path aware. New RecoilJets outputs belong under `runs/recoiljets/`,
Condor state under `state/condor/`, inputs under `inputs/`, model state under
`models/`, small evidence under `evidence/`, and transient material under
`scratch/`. Do not add new transfer mappings or helper defaults that create
root-level `output_*`, `tmp_*`, loose `.csv`/`.txt`, `config.log`,
`pythia_xsec_firstPass*`, or agent-named paths unless a local-only migration
manifest records the compatibility reason.

After THE-23 stage 7, campaign-aware SIM production should write future
checkout merge outputs under `runs/recoiljets/current/<campaign>/` and TG bulk
products under `thesisAna/recoiljets/<family>/<campaign>/{signal,background}/`.
Historical `output_<campaign>` and `thesisAna/simembedded*_<campaign>` paths
remain read-compatible legacy paths for old runs. Use
`scripts/sdcc/runtime/io/resolve_io_contract.py <campaign> --family <family>`
to inspect canonical and legacy paths, and use `SFTP_GET_CAMPAIGN_TAG=<tag>`
with `scripts/sftp_get_recoiljets_outputs.sh` when pulling a new canonical
campaign.
