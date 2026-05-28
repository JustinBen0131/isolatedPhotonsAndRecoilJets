#!/usr/bin/env python3
"""Create and restore private local snapshots of ThesisAnalysis OS state."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import tarfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SNAPSHOT_DIR = Path("agent_context/local/os_snapshots")
EVENT_LOG = Path("agent_context/local/os_events.jsonl")

CORE_FILES = [
    Path("AGENTS.md"),
    Path(".gitignore"),
    Path("agent_context/CODEX_WORK_REGISTER.yaml"),
    Path("agent_context/ARTIFACT_REGISTRY.yaml"),
    Path("agent_context/REFERENCE_MAP.md"),
    Path("agent_context/THESIS_NARRATIVE_MAP.md"),
    Path("agent_context/templates/OS_POSTMORTEM_TEMPLATE.md"),
]

CORE_GLOBS = [
    "agent_context/policies/*",
    "scripts/codex_os_*.py",
    "scripts/codex_work_register_*.py",
]


def now_id() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def slugify(text: str) -> str:
    slug = re.sub(r"[^a-zA-Z0-9_.-]+", "-", text.strip()).strip("-").lower()
    return slug[:60] or "snapshot"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def iter_snapshot_paths(extra_paths: list[str] | None = None) -> list[Path]:
    paths: set[Path] = set()
    for path in CORE_FILES:
        if path.exists() and path.is_file():
            paths.add(path)
    for pattern in CORE_GLOBS:
        for path in Path(".").glob(pattern):
            if path.is_file() and "agent_context/local" not in path.as_posix():
                paths.add(path)
    for raw in extra_paths or []:
        path = Path(raw)
        if path.is_file():
            paths.add(path)
        elif path.is_dir():
            for child in path.rglob("*"):
                if child.is_file() and "agent_context/local" not in child.as_posix():
                    paths.add(child)
    return sorted(paths, key=lambda p: p.as_posix())


def manifest_for(snapshot_id: str, label: str, paths: list[Path]) -> dict[str, Any]:
    return {
        "snapshot_id": snapshot_id,
        "label": label,
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "file_count": len(paths),
        "files": [
            {
                "path": path.as_posix(),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }
            for path in paths
        ],
    }


def create_snapshot(label: str, extra_paths: list[str] | None = None) -> dict[str, Any]:
    SNAPSHOT_DIR.mkdir(parents=True, exist_ok=True)
    snapshot_id = f"{now_id()}-{slugify(label)}"
    tar_path = SNAPSHOT_DIR / f"{snapshot_id}.tar.gz"
    manifest_path = SNAPSHOT_DIR / f"{snapshot_id}.manifest.json"
    paths = iter_snapshot_paths(extra_paths)
    manifest = manifest_for(snapshot_id, label, paths)
    with tarfile.open(tar_path, "w:gz") as tar:
        for path in paths:
            tar.add(path, arcname=path.as_posix())
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    payload = {**manifest, "tar_path": tar_path.as_posix(), "manifest_path": manifest_path.as_posix()}
    append_event(
        {
            "event_type": "os_snapshot",
            "snapshot_id": snapshot_id,
            "label": label,
            "file_count": len(paths),
            "tar_path": tar_path.as_posix(),
            "manifest_path": manifest_path.as_posix(),
        }
    )
    return payload


def append_event(event: dict[str, Any]) -> None:
    EVENT_LOG.parent.mkdir(parents=True, exist_ok=True)
    event.setdefault("timestamp", datetime.now(timezone.utc).isoformat(timespec="seconds"))
    with EVENT_LOG.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(event, sort_keys=True, separators=(",", ":")) + "\n")


def command_create(args: argparse.Namespace) -> int:
    manifest = create_snapshot(args.label, args.include_path)
    if args.json:
        print(json.dumps(manifest, indent=2, sort_keys=True))
    else:
        print(f"OK: snapshot {manifest['snapshot_id']}")
        print(f"tar={manifest['tar_path']}")
        print(f"manifest={manifest['manifest_path']}")
        print(f"files={manifest['file_count']}")
    return 0


def load_manifests() -> list[dict[str, Any]]:
    manifests: list[dict[str, Any]] = []
    if not SNAPSHOT_DIR.exists():
        return manifests
    for path in sorted(SNAPSHOT_DIR.glob("*.manifest.json")):
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            data = {"snapshot_id": path.stem.replace(".manifest", ""), "manifest_path": path.as_posix(), "invalid": True}
        if isinstance(data, dict):
            data.setdefault("manifest_path", path.as_posix())
            manifests.append(data)
    return manifests


def command_list(args: argparse.Namespace) -> int:
    manifests = load_manifests()
    if args.json:
        print(json.dumps({"snapshots": manifests}, indent=2, sort_keys=True))
    elif not manifests:
        print("No OS snapshots found.")
    else:
        for item in manifests:
            print(
                f"{item.get('snapshot_id')} files={item.get('file_count')} "
                f"created={item.get('created_at')} label={item.get('label')}"
            )
    return 0


def resolve_snapshot(snapshot_id: str) -> tuple[Path, Path]:
    base = SNAPSHOT_DIR / snapshot_id
    tar_path = Path(snapshot_id) if Path(snapshot_id).exists() else base.with_suffix(".tar.gz")
    manifest_path = base.with_suffix(".manifest.json")
    if not tar_path.exists():
        tar_path = SNAPSHOT_DIR / f"{snapshot_id}.tar.gz"
    if not manifest_path.exists():
        manifest_path = SNAPSHOT_DIR / f"{snapshot_id}.manifest.json"
    return tar_path, manifest_path


def safe_members(tar: tarfile.TarFile) -> list[tarfile.TarInfo]:
    members: list[tarfile.TarInfo] = []
    for member in tar.getmembers():
        path = Path(member.name)
        if path.is_absolute() or ".." in path.parts:
            raise RuntimeError(f"unsafe path in snapshot: {member.name}")
        members.append(member)
    return members


def command_inspect(args: argparse.Namespace) -> int:
    tar_path, manifest_path = resolve_snapshot(args.snapshot_id)
    if not tar_path.exists() or not manifest_path.exists():
        print(f"ERROR: snapshot not found: {args.snapshot_id}")
        return 1
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if args.json:
        print(json.dumps(manifest, indent=2, sort_keys=True))
    else:
        print(f"snapshot_id={manifest.get('snapshot_id')}")
        print(f"label={manifest.get('label')}")
        print(f"created_at={manifest.get('created_at')}")
        print(f"file_count={manifest.get('file_count')}")
        for item in manifest.get("files", []):
            print(f"- {item.get('path')} bytes={item.get('bytes')} sha256={item.get('sha256')}")
    return 0


def command_restore(args: argparse.Namespace) -> int:
    tar_path, manifest_path = resolve_snapshot(args.snapshot_id)
    if not tar_path.exists() or not manifest_path.exists():
        print(f"ERROR: snapshot not found: {args.snapshot_id}")
        return 1
    with tarfile.open(tar_path, "r:gz") as tar:
        members = safe_members(tar)
        if not args.apply:
            print(f"DRY-RUN: would restore {len(members)} files from {tar_path}")
            for member in members:
                print(f"- {member.name}")
            return 0
        pre = create_snapshot(f"pre-restore-{args.snapshot_id}")
        tar.extractall(path=Path("."), members=members)
    print(f"OK: restored {len(members)} files from {tar_path}")
    print(f"pre_restore_snapshot={pre['snapshot_id']}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    create = subparsers.add_parser("create", help="create a private local OS snapshot")
    create.add_argument("--label", required=True)
    create.add_argument("--include-path", action="append")
    create.add_argument("--json", action="store_true")
    create.set_defaults(func=command_create)

    list_cmd = subparsers.add_parser("list", help="list private local OS snapshots")
    list_cmd.add_argument("--json", action="store_true")
    list_cmd.set_defaults(func=command_list)

    inspect = subparsers.add_parser("inspect", help="show a snapshot manifest")
    inspect.add_argument("snapshot_id")
    inspect.add_argument("--json", action="store_true")
    inspect.set_defaults(func=command_inspect)

    restore = subparsers.add_parser("restore", help="restore a snapshot; dry-run unless --apply is passed")
    restore.add_argument("snapshot_id")
    restore.add_argument("--apply", action="store_true")
    restore.set_defaults(func=command_restore)

    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
