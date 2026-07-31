#!/usr/bin/env python3
"""Bounded read-only guard for shared SDCC scientific paths.

The guard never walks a tree.  It checks at most 256 explicitly named paths
and their ancestor directories, preventing a permission preflight from
becoming another high-cardinality login-node workload.
"""

from __future__ import annotations

import argparse
import json
import os
import stat
import sys
from pathlib import Path
from typing import Iterable


SCHEMA = "SDCC_SCIENTIFIC_PERMISSION_GUARD_V1"
MAX_PATHS = 256


class PermissionContractError(RuntimeError):
    """A requested path is unsafe or outside the declared scientific root."""


def _mode(path: Path) -> int:
    return stat.S_IMODE(path.stat().st_mode)


def _beneath(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
    except ValueError:
        return False
    return True


def _existing_target(path: Path, *, allow_missing_leaf: bool) -> tuple[Path, bool]:
    if path.exists():
        return path.resolve(strict=True), True
    if not allow_missing_leaf:
        raise PermissionContractError(f"path does not exist: {path}")
    parent = path.parent.resolve(strict=True)
    return parent, False


def inspect_paths(
    root: Path,
    paths: Iterable[Path],
    *,
    allow_missing_leaf: bool = False,
) -> dict[str, object]:
    requested = list(paths)
    if len(requested) > MAX_PATHS:
        raise PermissionContractError(
            f"explicit path count exceeds bounded limit: {len(requested)} > {MAX_PATHS}"
        )
    resolved_root = root.resolve(strict=True)
    if not resolved_root.is_dir():
        raise PermissionContractError(f"scientific root is not a directory: {root}")
    targets = requested or [root]
    failures: list[dict[str, object]] = []
    checked: set[Path] = set()
    for raw in targets:
        if not raw.is_absolute() or str(raw) != str(raw).strip():
            raise PermissionContractError(f"path must be a clean absolute path: {raw!s}")
        target, target_exists = _existing_target(
            raw,
            allow_missing_leaf=allow_missing_leaf,
        )
        if not _beneath(target, resolved_root):
            raise PermissionContractError(
                f"path escapes scientific root: path={raw}, root={root}"
            )
        relative = target.relative_to(resolved_root)
        chain = [resolved_root]
        cursor = resolved_root
        for part in relative.parts:
            cursor = cursor / part
            chain.append(cursor)
        for candidate in chain:
            if candidate in checked:
                continue
            checked.add(candidate)
            mode = _mode(candidate)
            if candidate.is_dir():
                required = 0o050
                ok = mode & required == required
                kind = "directory"
            elif candidate.is_file():
                required = 0o040
                ok = mode & required == required
                kind = "file"
            else:
                ok = False
                required = 0
                kind = "unsupported"
            if not ok:
                failures.append(
                    {
                        "path": str(candidate),
                        "kind": kind,
                        "mode": f"{mode:04o}",
                        "required_group_bits": f"{required:04o}",
                    }
                )
        if target_exists and target.is_file() and target not in checked:
            mode = _mode(target)
            if not mode & 0o040:
                failures.append(
                    {
                        "path": str(target),
                        "kind": "file",
                        "mode": f"{mode:04o}",
                        "required_group_bits": "0040",
                    }
                )
    return {
        "schema": SCHEMA,
        "status": "PASS" if not failures else "FAIL",
        "root": str(resolved_root),
        "requested_path_count": len(targets),
        "checked_inode_count": len(checked),
        "bounded_path_limit": MAX_PATHS,
        "tree_walk_performed": False,
        "failures": failures,
    }


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Check explicit SDCC scientific paths for group visibility."
    )
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--path", type=Path, action="append", default=[])
    parser.add_argument("--allow-missing-leaf", action="store_true")
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    try:
        args = parse_args(argv)
        report = inspect_paths(
            args.root,
            args.path,
            allow_missing_leaf=args.allow_missing_leaf,
        )
    except (PermissionContractError, FileNotFoundError, OSError) as exc:
        print(f"[SDCC-PERMISSION-GUARD][ERROR] {exc}", file=sys.stderr)
        return 2
    print(json.dumps(report, sort_keys=True))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
