#!/usr/bin/env python3
"""Read-only hygiene probe for the SDCC thesisAnalysis checkout."""

from __future__ import annotations

import argparse
import os
from pathlib import Path


SKIP_DIR_NAMES = {
    "output",
    "outputSmoke",
    "stdout",
    "error",
    "log",
    "condor_sub",
    "condor_snapshots",
    "condor_recovery",
    "condor_segments",
    "bdt_models",
    "mlp_models",
    "logreg_models",
    "dataOutput",
    "InputFiles",
    "OutDir",
    "transfer",
    "transfer_lite",
    "local_bdt_training_outputs",
    "local_ml_pipeline_tests",
}

BASE_FACING_ROOTS = (
    ".",
    "scripts",
    "macros",
    "coresoftware_local",
)

STAGE6_ALLOWED_ROOT_OUTPUTS = {
    # Held at SDCC root during THE-23 stage 6 because current docs/macros still
    # reference these exact paths. New output_* roots should be born under runs/.
    "output_inclusive3_fixed_stitch_20260519_203527",
    "output_focus21_stitch_20260521_123833_exclusive31",
}


def is_bad_visible_name(name: str) -> bool:
    return (
        name.strip() != name
        or name.strip() == ""
        or name == "agent_context"
        or name == ".codex"
        or name.startswith("codex")
        or name.startswith("THE-")
    )


def is_root_contract_clutter(path: Path) -> bool:
    name = path.name
    if name in STAGE6_ALLOWED_ROOT_OUTPUTS:
        return False
    if name.startswith("output_"):
        return True
    if name.startswith("tmp_"):
        return True
    if name.endswith((".csv", ".txt")):
        return True
    if name in {"config.log"}:
        return True
    if name.startswith((".RJ_OUTPUTS_WIPED_", "pythia_xsec_firstPass")):
        return True
    return False


def is_generated_trash(path: Path) -> bool:
    name = path.name
    return (
        name == "__pycache__"
        or name == ".DS_Store"
        or name.startswith("._")
        or name.endswith((".pyc", ".pyo"))
        or ".codex_" in name
    )


def walk_base_facing(root: Path) -> list[Path]:
    findings: list[Path] = []
    for root_name in BASE_FACING_ROOTS:
        start = root / root_name
        if not start.exists():
            continue
        for dirpath, dirnames, filenames in os.walk(start):
            current = Path(dirpath)
            rel = current.relative_to(root)
            if root_name == ".":
                dirnames[:] = [
                    d
                    for d in dirnames
                    if d not in SKIP_DIR_NAMES and not d.startswith("output_")
                ]
                if str(rel) != "." and len(rel.parts) >= 2:
                    dirnames[:] = []
            for dirname in list(dirnames):
                candidate = current / dirname
                if is_generated_trash(candidate):
                    findings.append(candidate.relative_to(root))
            for filename in filenames:
                candidate = current / filename
                if is_generated_trash(candidate):
                    findings.append(candidate.relative_to(root))
    return sorted(set(findings), key=lambda p: str(p))


def broken_script_symlinks(root: Path) -> list[Path]:
    scripts = root / "scripts"
    if not scripts.exists():
        return []
    broken: list[Path] = []
    for path in scripts.rglob("*"):
        if path.is_symlink() and not path.exists():
            broken.append(path.relative_to(root))
    return sorted(broken, key=lambda p: str(p))


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Report SDCC checkout hygiene issues without mutating files."
    )
    parser.add_argument(
        "checkout",
        nargs="?",
        default=".",
        help="Checkout directory to inspect; defaults to current directory.",
    )
    parser.add_argument(
        "--report-only",
        action="store_true",
        help="Always exit 0 after printing findings.",
    )
    args = parser.parse_args()

    root = Path(args.checkout).resolve()
    if not root.is_dir():
        raise SystemExit(f"[ERROR] Not a directory: {root}")

    top_entries = sorted(root.iterdir(), key=lambda p: p.name)
    suspicious_top = [p for p in top_entries if is_bad_visible_name(p.name)]
    root_contract_clutter = [
        p for p in top_entries if is_root_contract_clutter(p)
    ]
    generated = walk_base_facing(root)
    broken = broken_script_symlinks(root)

    print(f"CHECKOUT {root}")
    print(f"TOPLEVEL_COUNT {len(top_entries)}")
    print(f"SUSPICIOUS_TOPLEVEL {len(suspicious_top)}")
    for path in suspicious_top:
        kind = "dir" if path.is_dir() else "file"
        print(f"  {kind}\t{path.name!r}")
    print(f"ROOT_CONTRACT_CLUTTER {len(root_contract_clutter)}")
    for path in root_contract_clutter:
        kind = "dir" if path.is_dir() else "file"
        print(f"  {kind}\t{path.name!r}")
    print(f"GENERATED_TRASH {len(generated)}")
    for path in generated[:200]:
        print(f"  {path}")
    if len(generated) > 200:
        print(f"  ... {len(generated) - 200} more")
    print(f"BROKEN_SCRIPTS_SYMLINKS {len(broken)}")
    for path in broken[:200]:
        print(f"  {path}")
    if len(broken) > 200:
        print(f"  ... {len(broken) - 200} more")

    issues = (
        len(suspicious_top)
        + len(root_contract_clutter)
        + len(generated)
        + len(broken)
    )
    if issues:
        print(f"RESULT FAIL issues={issues}")
        return 0 if args.report_only else 1
    print("RESULT PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
