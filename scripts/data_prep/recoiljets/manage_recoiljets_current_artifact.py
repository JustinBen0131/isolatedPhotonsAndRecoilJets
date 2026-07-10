#!/usr/bin/env python3
"""Register and resolve current RecoilJets production artifacts.

The registry is intentionally local and non-destructive: campaign output
directories stay immutable, while `current/<sample_key>/current.json` records
which final ROOT(s) plotting code should use by default.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from typing import Any


REPO = Path(__file__).resolve().parents[3]
REGISTRY_DIR = REPO / "dataOutput/current_recoiljets_artifacts"
REGISTRY_PATH = REGISTRY_DIR / "registry.json"


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def load_registry(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"schema_version": 1, "artifacts": [], "current": {}}
    with path.open() as handle:
        return json.load(handle)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def write_entry_snapshots(registry_path: Path, entry: dict[str, Any]) -> None:
    sample_key = entry["sample_key"]
    entry_id = entry["id"]
    write_json(registry_path.parent / "entries" / sample_key / f"{entry_id}.json", entry)
    if entry.get("status") in {"superseded", "bad", "paused"}:
        write_json(registry_path.parent / "backlog" / sample_key / f"{entry_id}.json", entry)


def write_current_pointer(registry_path: Path, entry: dict[str, Any]) -> None:
    sample_key = entry["sample_key"]
    root_paths = entry.get("root_paths", [])
    current_dir = registry_path.parent / "current" / sample_key
    pointer = {
        "schema_version": 1,
        "sample_key": sample_key,
        "current_entry_id": entry["id"],
        "campaign_tag": entry.get("campaign_tag", ""),
        "root_paths": root_paths,
        "artifact_kind": entry.get("artifact_kind", ""),
        "role": entry.get("role", ""),
        "produced_at": entry.get("produced_at", ""),
        "registered_at": entry.get("registered_at", ""),
        "registry": str(registry_path),
        "plot_policy": entry.get("plot_policy", ""),
        "sample_lane": entry.get("sample_lane", ""),
        "canonical_status": entry.get("canonical_status", "not_canonical"),
        "contract_report": entry.get("contract_report", ""),
        "promotion_basis": entry.get("promotion_basis", ""),
        "waivers": entry.get("waivers", []),
        "notes": entry.get("notes", ""),
    }
    write_json(current_dir / "current.json", pointer)
    link = current_dir / "current.root"
    if link.exists() or link.is_symlink():
        link.unlink()
    if len(root_paths) == 1:
        link.symlink_to(Path(root_paths[0]))


def normalize_path(path: str) -> str:
    p = Path(path).expanduser()
    if not p.is_absolute():
        p = (REPO / p).resolve()
    return str(p)


def verify_roots(paths: list[str]) -> list[dict[str, Any]]:
    evidence = []
    for item in paths:
        p = Path(item)
        evidence.append(
            {
                "path": str(p),
                "exists": p.exists(),
                "bytes": p.stat().st_size if p.exists() else None,
                "mtime": datetime.fromtimestamp(p.stat().st_mtime, timezone.utc).replace(microsecond=0).isoformat()
                if p.exists()
                else None,
            }
        )
    return evidence


def command_register(args: argparse.Namespace) -> int:
    registry = load_registry(args.registry)
    timestamp = args.produced_at or now_iso()
    root_paths = [normalize_path(p) for p in args.root_path]
    evidence = verify_roots(root_paths)
    missing = [row["path"] for row in evidence if not row["exists"]]
    if missing and not args.allow_missing:
        for path in missing:
            print(f"[ERROR] root path does not exist: {path}", file=sys.stderr)
        return 2

    entry_id = f"{timestamp}_{args.sample_key}_{args.campaign_tag}".replace(":", "").replace("+", "p")
    entry = {
        "id": entry_id,
        "sample_key": args.sample_key,
        "sample_family": args.sample_family,
        "artifact_kind": args.artifact_kind,
        "campaign_tag": args.campaign_tag,
        "status": args.status,
        "role": args.role,
        "root_paths": root_paths,
        "remote_paths": args.remote_path,
        "produced_at": timestamp,
        "registered_at": now_iso(),
        "plot_policy": args.plot_policy,
        "sample_lane": args.sample_lane,
        "canonical_status": args.canonical_status,
        "contract_report": args.contract_report,
        "promotion_basis": args.promotion_basis,
        "waivers": args.waiver,
        "notes": args.notes,
        "evidence": evidence,
    }

    previous_current = registry.get("current", {}).get(args.sample_key)
    if args.status == "current" and previous_current:
        for old in registry.get("artifacts", []):
            if old.get("id") == previous_current:
                old["status"] = "superseded"
                old["superseded_at"] = now_iso()
                old["superseded_by"] = entry_id
                write_entry_snapshots(args.registry, old)
                break
        entry["previous_current"] = previous_current

    registry.setdefault("artifacts", []).append(entry)
    if args.status == "current":
        registry.setdefault("current", {})[args.sample_key] = entry_id

    write_json(args.registry, registry)
    write_entry_snapshots(args.registry, entry)

    if args.status == "current":
        write_current_pointer(args.registry, entry)

    print(entry_id)
    return 0


def find_current(registry: dict[str, Any], sample_key: str) -> dict[str, Any]:
    entry_id = registry.get("current", {}).get(sample_key)
    if not entry_id:
        raise KeyError(f"no current artifact registered for sample_key={sample_key}")
    for entry in registry.get("artifacts", []):
        if entry.get("id") == entry_id:
            return entry
    raise KeyError(f"current entry id {entry_id} is missing from registry")


def command_resolve(args: argparse.Namespace) -> int:
    registry = load_registry(args.registry)
    try:
        entry = find_current(registry, args.sample_key)
    except KeyError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 2
    if args.json:
        print(json.dumps(entry, indent=2, sort_keys=True))
    else:
        for path in entry.get("root_paths", []):
            print(path)
    return 0


def command_list(args: argparse.Namespace) -> int:
    registry = load_registry(args.registry)
    rows = registry.get("artifacts", [])
    if args.sample_key:
        rows = [row for row in rows if row.get("sample_key") == args.sample_key]
    if args.json:
        print(json.dumps(rows, indent=2, sort_keys=True))
        return 0
    for row in rows:
        current_marker = " current" if registry.get("current", {}).get(row.get("sample_key")) == row.get("id") else ""
        print(
            f"{row.get('sample_key')} {row.get('status')}{current_marker} "
            f"{row.get('campaign_tag')} {row.get('produced_at')} {row.get('id')}"
        )
    return 0


def command_sync_snapshots(args: argparse.Namespace) -> int:
    registry = load_registry(args.registry)
    entry_by_id = {entry.get("id"): entry for entry in registry.get("artifacts", [])}
    for entry in registry.get("artifacts", []):
        write_entry_snapshots(args.registry, entry)
    for sample_key, entry_id in registry.get("current", {}).items():
        entry = entry_by_id.get(entry_id)
        if not entry:
            print(f"[ERROR] current entry id {entry_id} for {sample_key} is missing", file=sys.stderr)
            return 2
        if entry.get("status") != "current":
            print(f"[ERROR] current entry id {entry_id} for {sample_key} has status={entry.get('status')}", file=sys.stderr)
            return 2
        write_current_pointer(args.registry, entry)
    print(f"synced {len(registry.get('artifacts', []))} entries")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, default=REGISTRY_PATH)
    sub = parser.add_subparsers(dest="command", required=True)

    reg = sub.add_parser("register")
    reg.add_argument("--sample-key", required=True, help="Stable key, e.g. pp_sim_photonjet_merged")
    reg.add_argument("--sample-family", default="pp")
    reg.add_argument("--artifact-kind", default="final_root")
    reg.add_argument("--campaign-tag", required=True)
    reg.add_argument("--status", choices=["current", "candidate", "paused", "superseded", "bad"], default="current")
    reg.add_argument("--role", required=True)
    reg.add_argument("--root-path", action="append", required=True)
    reg.add_argument("--remote-path", action="append", default=[])
    reg.add_argument("--produced-at")
    reg.add_argument("--plot-policy", default="")
    reg.add_argument("--sample-lane", default="")
    reg.add_argument(
        "--canonical-status",
        choices=["not_canonical", "candidate_evidence", "canonical"],
        default="not_canonical",
    )
    reg.add_argument("--contract-report", default="")
    reg.add_argument("--promotion-basis", default="")
    reg.add_argument("--waiver", action="append", default=[])
    reg.add_argument("--notes", default="")
    reg.add_argument("--allow-missing", action="store_true")
    reg.set_defaults(func=command_register)

    res = sub.add_parser("resolve")
    res.add_argument("--sample-key", required=True)
    res.add_argument("--json", action="store_true")
    res.set_defaults(func=command_resolve)

    ls = sub.add_parser("list")
    ls.add_argument("--sample-key")
    ls.add_argument("--json", action="store_true")
    ls.set_defaults(func=command_list)

    sync = sub.add_parser("sync-snapshots")
    sync.set_defaults(func=command_sync_snapshots)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
