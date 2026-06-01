#!/usr/bin/env python3
from pathlib import Path
import csv, sys
repo = Path(__file__).resolve().parents[2]
scripts_root = repo / "scripts"
map_path = scripts_root / "sdcc" / "TRANSFER_MAP.tsv"
errors = []
with map_path.open(newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))
if not rows:
    errors.append("transfer map is empty")
for i, row in enumerate(rows, 2):
    local = Path(row["local_canonical"])
    if row["local_canonical"].startswith("scripts/") or row["local_canonical"].startswith("macros/") or row["local_canonical"].startswith("src"):
        if not local.exists() and not local.is_symlink():
            errors.append(f"line {i}: missing local canonical {local}")
    alias = row.get("local_alias", "")
    if alias.startswith("scripts/") and alias != row["local_canonical"]:
        compat = Path("scripts/compat/local") / Path(alias).name
        hard = Path(alias)
        if not hard.exists() and not hard.is_symlink() and not compat.exists() and not compat.is_symlink():
            errors.append(f"line {i}: missing alias/compat for {alias}")
    if not row["upload_target"]:
        errors.append(f"line {i}: empty upload target")
if errors:
    print("[FAIL] transfer map validation")
    for e in errors:
        print(" -", e)
    sys.exit(1)
print(f"[OK] transfer map rows={len(rows)}")
