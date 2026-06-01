# Macros Directory

`macros/` is now a small ROOT/Fun4All contract surface plus indexed local/offline macro domains.

## Protected Top-Level Anchors

Keep these paths fixed unless a separate audited compatibility migration is approved:

- `Fun4All_recoilJets*.C` and `Fun4All_auauTightBDTTraining.C`
- `analysis_config*.yaml`
- `Calo_Calib.C` and `HIJetReco.C`
- `AnalyzeRecoilJets*`, `AnalyzeTriggerGroupings.h`
- `sPhenixStyle.C` and `sPhenixStyle.h`

These are path-sensitive for SDCC, ROOT, existing plotting includes, or transfer helpers.

## Canonical Local/Offline Homes

- `plotting/auau_bdt/`
- `plotting/target_wp/`
- `plotting/width_study/`
- `plotting/stitching/`
- `plotting/pp_currentian/`
- `plotting/ssqa/`
- `diagnostics/auau_bdt/`
- `diagnostics/stitching/`
- `diagnostics/preselection/`
- `diagnostics/ssqa/`
- `diagnostics/model_checks/`
- `diagnostics/feasibility/`

Top-level wrappers for moved `.C` macros are hard compatibility aliases only
where old paths are actively referenced by policies, transfer maps, status
ledgers, or likely operator muscle memory. Lower-risk compatibility wrappers
live under `compat/local/`. Edit canonical files in the domain folders.

## Find A Macro

```bash
macros/bin/thesis-macro path PlotTarget80FoodChain
macros/bin/thesis-macro list plotting/width_study
```

`MACRO_INDEX.yaml` is the durable agent-facing index. `MACRO_INDEX.tsv` is the shell-friendly resolver index.

## Generated Artifacts

Do not keep ACLiC products, `.DS_Store`, AppleDouble files, `__pycache__`, or other generated build/cache artifacts in the front-facing macro tree. Use external cache/build locations or remove generated products after checks.
