# ROOT Diagnostic Macros

Canonical diagnostic macro implementations live in domain subfolders:

- `auau_bdt/`: Au+Au BDT comparison and box/NPB checks.
- `stitching/`: stitching and do-not-scale histogram inspections.
- `preselection/`: fail-fast preselection checks.
- `ssqa/`: SSQA quick-run diagnostics.
- `model_checks/`: model sanity summaries.
- `feasibility/`: exploratory feasibility diagnostics.

Top-level compatibility wrappers exist only to keep old ROOT commands runnable.
