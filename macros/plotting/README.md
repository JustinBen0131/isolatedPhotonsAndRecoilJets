# ROOT Plotting Macros

Canonical offline plotting macro implementations live in domain subfolders:

- `auau_bdt/`: Au+Au BDT, efficiency, centrality, and target-family plots.
- `target_wp/`: target working-point comparison plots.
- `width_study/`: BDT shower-width study plots.
- `stitching/`: embedded photon/jet stitching, xsec, and slide-style plot macros.
- `pp_currentian/`: pp CurrentIAN and pp-vs-AuAu comparison macros.
- `ssqa/`: SSQA summary plotting macros.

Top-level `macros/*.C` files for these entries are compatibility wrappers, not implementation sources. Use `macros/MACRO_INDEX.yaml` or `macros/bin/thesis-macro path <id-or-basename>` to find canonical source.
