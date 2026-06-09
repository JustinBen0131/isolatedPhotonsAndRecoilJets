# Scripts Directory

`scripts/` is a maintained command subsystem. The parent directory is a
minimal command surface: README/index files plus hard operator/runtime aliases.
Canonical source lives in purpose folders, old local names live in
`scripts/compat/local/`, and `scripts/bin/thesis-script` resolves IDs,
basenames, old paths, and canonical paths.

## Top-Level Contract

Top-level script names are hard contract symlinks unless the file is this
README or an index. Only active SDCC/runtime/transfer command contracts stay
front-facing. Historical aliases moved to
`scripts/compat/local/` and are indexed in `SCRIPT_INDEX.yaml`,
`HELPER_INDEX.yaml`, `slides/INDEX.yaml`, and `compat/ALIAS_INDEX.tsv`.

Keep these command paths available unless a coordinated compatibility migration
is explicitly approved:

- `RecoilJets_Condor*.sh`
- `mergeRecoilJets.sh`
- `recoiljets_cleanup.sh`
- `make_dstListsData.sh`
- `makeThesisSimLists.sh`
- `estimateEmbeddedPhotonXsec.sh`
- `root_in_analysis_env.sh`
- `sftp_push_recoiljets.sh`
- `sftp_get_recoiljets_outputs.sh`
- `auau_tight_*_pipeline.sh`
- `pp_photon_ml_pipeline.sh`

These are tied to SDCC transfer, Condor, ROOT setup, or existing operator
muscle memory. Edit their canonical source under the indexed purpose folder,
not through an assumed flat source layout.

## Canonical Homes

When adding new local-side code, use the narrowest purpose-specific folder:

- `scripts/bin/`: stable command gateways.
- `scripts/compat/local/`: old local names that should remain resolvable but
  should not clutter the parent.
- `scripts/sdcc/runtime/{condor,merge,lists,cleanup,xsec,audit}/`
- `scripts/sdcc/workflows/{submit,target_wp,width_study,stacking,scan,diagnostics}/`
- `scripts/sdcc/pipelines/{auau,pp}/`
- `scripts/plotting/{auau_bdt,pp_currentian,stitching,efficiency,trigger,truth_purity}/`
- `scripts/slides/<campaign>/<family>/`
- `scripts/ml/{training,validation,working_points,stacking,audits}/`
- `scripts/diagnostics/{pp_shuhang,auau_split,row_contracts,sdcc_audit,ml_validation}/`
- `scripts/data_prep/{manifests,ppg12,recoiljets,stitching}/`

Do not add new source to the top-level compatibility surface.

## SDCC Runtime Scripts

`scripts/sdcc/` is the canonical source home for scripts that may run on SDCC:

- `scripts/sdcc/runtime/{condor,merge,cleanup,lists,xsec,audit}/`: fixed-path
  Condor, merge, cleanup, DST-list, xsec, and audit entrypoints.
- `scripts/sdcc/pipelines/{auau,pp}/`: AuAu/pp ML pipeline drivers.
- `scripts/sdcc/workflows/{submit,target_wp,width_study,stacking,scan,diagnostics}/`:
  submit, merge, queue-gated, and campaign drivers.
- `scripts/sdcc/transfer/`: local SFTP push/get helpers.
- `scripts/sdcc/local_merge/root/`: local ROOT merge helper used after pulls.

Only SDCC-needed files should be selected by `sftp_push_recoiljets.sh`.
`scripts/sdcc/TRANSFER_MAP.tsv` is the readable transfer source of truth; the
shell helper arrays are regenerated from it. The remote checkout keeps hard
historical command names as symlinks or runtime aliases where needed.

## Slide Generators

`scripts/slides/` is the canonical slide-generator subsystem. Its `INDEX.yaml`
maps canonical source files to historical compatibility paths, campaign
families, slide keys, and default output folders. Use
`scripts/bin/thesis-script list slides` or `scripts/slides/INDEX.yaml` to find
slide code. Do not add slide builders to the parent directory.

## Indexes

- `scripts/COMMAND_INDEX.tsv`: command gateway lookup table.
- `scripts/SCRIPT_INDEX.yaml`: generated full script surface map.
- `scripts/HELPER_INDEX.yaml`: generated non-slide helper map.
- `scripts/slides/INDEX.yaml`: generated slide-generator map.
- `scripts/sdcc/SDCC_RUNTIME_INDEX.yaml`: generated SDCC script map.
- `scripts/sdcc/TRANSFER_MAP.tsv`: structured local-to-SDCC transfer map.
- `scripts/compat/ALIAS_INDEX.tsv`: hard and compat alias map.

Edit canonical source under the purpose folders; use old names only through
hard contracts, `scripts/compat/local/`, or `scripts/bin/thesis-script`.

## Codex Rule

Before creating another top-level script, search for existing helpers with
`rg --files scripts macros`, choose the narrowest reusable location, and avoid
ambiguous names like `make_plot.py` or `test.py`.
