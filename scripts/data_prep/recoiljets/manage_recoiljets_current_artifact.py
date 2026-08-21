#!/usr/bin/env python3
"""Register and resolve current RecoilJets production artifacts.

The registry is intentionally local and non-destructive: campaign output
directories stay immutable, while `current/<sample_key>/current.json` records
which final ROOT(s) plotting code should use by default.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from typing import Any


REPO = Path(__file__).resolve().parents[3]
REGISTRY_DIR = REPO / "dataOutput/current_recoiljets_artifacts"
REGISTRY_PATH = REGISTRY_DIR / "registry.json"

PPG12_PRODUCTION_GATED_SAMPLE_FAMILIES = {
    "pp_sim_inclusivejet_merged": "inclusive",
    "pp_sim_photonjet_merged": "photon",
}
PPG12_PRODUCTION_GATE_SCHEMA = "ppg12-stitched-purity-production-gate/v1"
PPG12_GATE_REPORT_SCHEMA = "ppg12-stitched-purity-gate-report/v1"
PPG12_ADMISSION_SCHEMA = "ppg12-stitched-purity-admission/v1"
PPG12_PRODUCTION_WRAPPER_SCHEMA = "ppg12-stitched-purity-production/v1"
PPG12_MERGE_AUDIT_SCHEMA = "ppg12-stitched-purity-merge-audit/v1"
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
PPG12_CLOSURE_GATE = (
    REPO
    / "scripts/diagnostics/pp_currentian/ppg12_stitched_purity_closure_gate.py"
)


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
        "production_gate_report": entry.get("production_gate_report", ""),
        "production_gate_report_sha256": entry.get(
            "production_gate_report_sha256", ""
        ),
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


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_payload_sha256(payload: Any) -> str:
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def replay_ppg12_production_gate(
    contract_path: Path,
    admission_path: Path,
    production_wrapper_path: Path,
    merge_path: Path,
) -> dict[str, Any]:
    """Regenerate the production receipt from canonical executable evidence."""

    if not PPG12_CLOSURE_GATE.is_file():
        raise ValueError(f"canonical PPG12 closure gate is missing: {PPG12_CLOSURE_GATE}")
    with tempfile.TemporaryDirectory(prefix="ppg12-current-promotion-replay-") as temp_dir:
        outdir = Path(temp_dir)
        command = [
            sys.executable,
            str(PPG12_CLOSURE_GATE),
            "--contract",
            str(contract_path),
            "verify-production",
            "--admission-manifest",
            str(admission_path),
            "--production-manifest",
            str(production_wrapper_path),
            "--merge-audit",
            str(merge_path),
            "--outdir",
            str(outdir),
        ]
        completed = subprocess.run(
            command,
            text=True,
            capture_output=True,
            check=False,
        )
        replayed_path = outdir / "production_gate_report.json"
        if completed.returncode != 0 or not replayed_path.is_file():
            detail = completed.stderr.strip() or completed.stdout.strip()
            raise ValueError(
                "canonical PPG12 production-gate replay failed"
                + (f": {detail}" if detail else "")
            )
        return read_json_object(replayed_path, "replayed production gate report")


def read_json_object(path: Path, label: str) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text())
    except FileNotFoundError as exc:
        raise ValueError(f"{label} does not exist: {path}") from exc
    except json.JSONDecodeError as exc:
        raise ValueError(f"{label} is not valid JSON: {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"{label} must contain a JSON object: {path}")
    return payload


def require_sha256(value: Any, label: str) -> str:
    if not isinstance(value, str) or not SHA256_RE.fullmatch(value):
        raise ValueError(f"{label} is not a lowercase SHA-256")
    return value


def resolve_evidence_path(owner: Path, raw_path: Any, label: str) -> Path:
    if not isinstance(raw_path, str) or not raw_path:
        raise ValueError(f"{label} path is missing")
    path = Path(raw_path).expanduser()
    if not path.is_absolute():
        path = owner.parent / path
    try:
        return path.resolve(strict=True)
    except (FileNotFoundError, OSError) as exc:
        raise ValueError(f"{label} path cannot be resolved: {path}: {exc}") from exc


def verify_file_sha256(path: Path, expected: Any, label: str) -> str:
    expected_sha = require_sha256(expected, f"{label} SHA-256")
    actual_sha = file_sha256(path)
    if actual_sha != expected_sha:
        raise ValueError(
            f"{label} hash drift: path={path} expected={expected_sha} "
            f"actual={actual_sha}"
        )
    return actual_sha


def resolve_json_link(
    owner: Path,
    link: Any,
    label: str,
) -> tuple[Path, dict[str, Any]]:
    if not isinstance(link, dict):
        raise ValueError(f"{label} link is missing")
    path = resolve_evidence_path(owner, link.get("path"), label)
    verify_file_sha256(path, link.get("sha256"), label)
    return path, read_json_object(path, label)


def find_merge_audit_by_sha256(
    expected_sha: str,
    search_dirs: list[Path],
    candidate_manifest_path: Path,
) -> tuple[Path, dict[str, Any]]:
    """Resolve a legacy production merge hash to its local audit sidecar.

    The v1 production-gate emitter records the production merge SHA but not
    its path.  Preserve that workflow by searching only the direct evidence
    directories already linked by the gate, and accept a match only when it
    is a real merge-audit JSON bound to the exact production candidate.
    Newer reports may provide an explicit ``merge_audit`` link and bypass this
    compatibility resolver.
    """

    matches: list[tuple[Path, dict[str, Any]]] = []
    for directory in sorted({path.resolve() for path in search_dirs}):
        if not directory.is_dir():
            continue
        for candidate in sorted(directory.glob("*.json")):
            if not candidate.is_file() or file_sha256(candidate) != expected_sha:
                continue
            payload = read_json_object(candidate, "production merge audit")
            if payload.get("schema") != PPG12_MERGE_AUDIT_SCHEMA:
                continue
            link = payload.get("candidate_manifest")
            try:
                linked_path, _ = resolve_json_link(
                    candidate, link, "production merge-audit candidate_manifest"
                )
            except ValueError:
                continue
            if linked_path == candidate_manifest_path:
                matches.append((candidate.resolve(), payload))
    if not matches:
        raise ValueError(
            "production merge_audit_sha256 does not resolve to a linked "
            "merge-audit sidecar"
        )
    # Multiple byte-identical copies are equivalent because the report binds
    # the content hash and each accepted copy binds the same candidate file.
    return sorted(matches, key=lambda row: str(row[0]))[0]


def validate_ppg12_production_gate(
    report_path: Path,
    sample_key: str,
    root_paths: list[str],
) -> dict[str, Any]:
    """Bind a pp-SIM current-pointer mutation to the exact gated ROOT.

    The complete closure gate emits one receipt containing exactly one photon
    and one inclusive candidate artifact.  A current-pointer update is valid
    only when the receipt resolves to the zero-failure gate detail, admission,
    production wrapper, manifests, merge audit, historical evidence, and exact
    candidate ROOT bytes.  Merely constructing a PASS-shaped receipt or
    supplying syntactically valid 64-character hashes is insufficient.
    """

    family = PPG12_PRODUCTION_GATED_SAMPLE_FAMILIES[sample_key]
    report_path = report_path.expanduser().resolve()
    payload = read_json_object(report_path, "production gate report")
    if payload.get("schema") != PPG12_PRODUCTION_GATE_SCHEMA or payload.get("status") != "PASS":
        raise ValueError(
            "production gate report is not a passing "
            f"{PPG12_PRODUCTION_GATE_SCHEMA} report"
        )

    required_sha_fields = (
        "contract_sha256",
        "gate_report_payload_sha256",
        "admission_sha256",
        "production_manifest_sha256",
        "reference_manifest_sha256",
        "candidate_manifest_sha256",
        "merge_audit_sha256",
    )
    for key in required_sha_fields:
        require_sha256(payload.get(key), f"production gate report {key}")

    # The production_gate_report is an authorization receipt, not an
    # authority by itself.  Recompute the hash of the sibling gate report and
    # then walk every path+SHA link that the gate evaluated.
    gate_report_path = report_path.parent / "gate_report.json"
    gate_report = read_json_object(gate_report_path, "production gate detail report")
    if canonical_payload_sha256(gate_report) != payload["gate_report_payload_sha256"]:
        raise ValueError("production gate detail-report payload hash is inconsistent")
    if (
        gate_report.get("schema") != PPG12_GATE_REPORT_SCHEMA
        or gate_report.get("mode") != "verify-production"
        or gate_report.get("status") != "PASS"
        or gate_report.get("failure_count") != 0
        or gate_report.get("failures") != []
    ):
        raise ValueError(
            "production gate detail report is not a zero-failure "
            "verify-production PASS"
        )

    contract_link = gate_report.get("contract")
    if not isinstance(contract_link, dict):
        raise ValueError("production gate detail report contract link is missing")
    contract_path = resolve_evidence_path(
        gate_report_path, contract_link.get("path"), "closure contract"
    )
    contract = read_json_object(contract_path, "closure contract")
    contract_sha = canonical_payload_sha256(contract)
    if contract_link.get("sha256") != contract_sha or payload["contract_sha256"] != contract_sha:
        raise ValueError("production gate closure-contract hash is inconsistent")

    summary = gate_report.get("summary")
    if not isinstance(summary, dict):
        raise ValueError("production gate detail report summary is missing")
    admission_path = resolve_evidence_path(
        gate_report_path, summary.get("admission_manifest"), "admission manifest"
    )
    verify_file_sha256(
        admission_path, payload["admission_sha256"], "admission manifest"
    )
    admission = read_json_object(admission_path, "admission manifest")
    if (
        admission.get("schema") != PPG12_ADMISSION_SCHEMA
        or admission.get("status") != "PASS"
        or admission.get("contract_sha256") != contract_sha
    ):
        raise ValueError("linked admission manifest is not a passing contract match")
    for label in ("reference_manifest", "candidate_manifest", "merge_audit"):
        resolve_json_link(
            admission_path, admission.get(label), f"admission {label}"
        )
    admission_gate_path = admission_path.parent / "gate_report.json"
    admission_gate = read_json_object(admission_gate_path, "admission gate report")
    expected_admission_gate_sha = require_sha256(
        admission.get("gate_report_payload_sha256"),
        "admission gate_report_payload_sha256",
    )
    if canonical_payload_sha256(admission_gate) != expected_admission_gate_sha:
        raise ValueError("admission gate-report payload hash is inconsistent")
    if (
        admission_gate.get("schema") != PPG12_GATE_REPORT_SCHEMA
        or admission_gate.get("mode") != "admit"
        or admission_gate.get("status") != "PASS"
        or admission_gate.get("failure_count") != 0
        or admission_gate.get("failures") != []
    ):
        raise ValueError("linked admission gate report is not a zero-failure admit PASS")

    production_wrapper_path = resolve_evidence_path(
        gate_report_path, summary.get("production_wrapper"), "production wrapper"
    )
    verify_file_sha256(
        production_wrapper_path,
        payload["production_manifest_sha256"],
        "production wrapper",
    )
    wrapper = read_json_object(production_wrapper_path, "production wrapper")
    if wrapper.get("schema") != PPG12_PRODUCTION_WRAPPER_SCHEMA:
        raise ValueError("linked production wrapper has the wrong schema")
    if wrapper.get("admission_sha256") != payload["admission_sha256"]:
        raise ValueError("production wrapper is bound to a different admission")

    reference_path = resolve_evidence_path(
        gate_report_path,
        gate_report.get("reference_manifest"),
        "production reference manifest",
    )
    candidate_path = resolve_evidence_path(
        gate_report_path,
        gate_report.get("candidate_manifest"),
        "production candidate manifest",
    )
    verify_file_sha256(
        reference_path,
        payload["reference_manifest_sha256"],
        "production reference manifest",
    )
    verify_file_sha256(
        candidate_path,
        payload["candidate_manifest_sha256"],
        "production candidate manifest",
    )
    wrapper_reference_path, _ = resolve_json_link(
        production_wrapper_path,
        wrapper.get("reference_manifest"),
        "production-wrapper reference_manifest",
    )
    wrapper_candidate_path, _ = resolve_json_link(
        production_wrapper_path,
        wrapper.get("candidate_manifest"),
        "production-wrapper candidate_manifest",
    )
    if wrapper_reference_path != reference_path or wrapper_candidate_path != candidate_path:
        raise ValueError("production wrapper and gate report bind different manifests")

    merge_sha = payload["merge_audit_sha256"]
    merge_link = payload.get("merge_audit")
    if isinstance(merge_link, dict):
        merge_path, merge_audit = resolve_json_link(
            report_path, merge_link, "production merge audit"
        )
        if file_sha256(merge_path) != merge_sha:
            raise ValueError("production merge-audit link and SHA field disagree")
    else:
        merge_path, merge_audit = find_merge_audit_by_sha256(
            merge_sha,
            [report_path.parent, production_wrapper_path.parent, candidate_path.parent],
            candidate_path,
        )
    if merge_audit.get("schema") != PPG12_MERGE_AUDIT_SCHEMA:
        raise ValueError("linked production merge audit has the wrong schema")
    wrapper_merge_path, _ = resolve_json_link(
        production_wrapper_path,
        wrapper.get("merge_audit"),
        "production-wrapper merge_audit",
    )
    if wrapper_merge_path != merge_path or file_sha256(wrapper_merge_path) != merge_sha:
        raise ValueError("production wrapper and gate receipt bind different merge audits")
    merge_candidate_path, _ = resolve_json_link(
        merge_path,
        merge_audit.get("candidate_manifest"),
        "production merge-audit candidate_manifest",
    )
    if merge_candidate_path != candidate_path:
        raise ValueError("production merge audit is bound to a different candidate manifest")
    audits = merge_audit.get("audits")
    if not isinstance(audits, list) or not audits or any(
        not isinstance(row, dict)
        or str(row.get("status", "")).upper() != "PASS"
        or row.get("failures") not in ([], None)
        for row in audits
    ):
        raise ValueError("linked production merge audit is not passing")

    replayed_receipt = replay_ppg12_production_gate(
        contract_path,
        admission_path,
        production_wrapper_path,
        merge_path,
    )
    if replayed_receipt != payload:
        raise ValueError(
            "production gate receipt is not the canonical replay-derived result"
        )

    artifacts = payload.get("candidate_artifacts")
    if not isinstance(artifacts, list):
        raise ValueError("production gate report candidate_artifacts must be a list")
    artifact_payload: list[dict[str, str]] = []
    normalized: list[dict[str, str]] = []
    seen: set[str] = set()
    for index, row in enumerate(artifacts):
        if not isinstance(row, dict):
            raise ValueError(f"candidate_artifacts[{index}] must be an object")
        artifact_family = row.get("family")
        raw_path = row.get("path")
        sha256 = row.get("sha256")
        if artifact_family not in {"inclusive", "photon"}:
            raise ValueError(f"candidate_artifacts[{index}] has invalid family")
        if artifact_family in seen:
            raise ValueError(f"duplicate candidate artifact family: {artifact_family}")
        if not isinstance(raw_path, str):
            raise ValueError(f"candidate_artifacts[{index}] has invalid path/SHA-256")
        sha256 = require_sha256(
            sha256, f"candidate_artifacts[{index}] SHA-256"
        )
        seen.add(artifact_family)
        artifact_payload.append(
            {"family": artifact_family, "path": raw_path, "sha256": sha256}
        )
        normalized.append(
            {
                "family": artifact_family,
                "path": str(Path(raw_path).expanduser().resolve()),
                "sha256": sha256,
            }
        )
    artifact_payload.sort(key=lambda row: row["family"])
    normalized.sort(key=lambda row: row["family"])
    if seen != {"inclusive", "photon"}:
        raise ValueError("production gate report must bind both inclusive and photon artifacts")
    if payload.get("candidate_artifact_set_sha256") != canonical_payload_sha256(artifact_payload):
        raise ValueError("production gate report candidate artifact-set SHA-256 is inconsistent")
    wrapper_artifacts = wrapper.get("candidate_artifacts")
    if wrapper_artifacts != artifact_payload:
        raise ValueError("production wrapper and gate report bind different candidate artifacts")
    wrapper_assembly = wrapper.get("assembly")
    if not isinstance(wrapper_assembly, dict):
        raise ValueError("production wrapper assembly evidence is missing")
    if wrapper_assembly.get("contract_sha256") != contract_sha:
        raise ValueError("production wrapper closure-contract hash is stale")
    if wrapper_assembly.get("candidate_artifact_set_sha256") != canonical_payload_sha256(artifact_payload):
        raise ValueError("production wrapper candidate artifact-set hash is stale")
    resolve_json_link(
        production_wrapper_path,
        wrapper_assembly.get("historical_comparison"),
        "production historical comparison",
    )

    summary_artifacts = summary.get("candidate_artifacts")
    if summary_artifacts != artifact_payload:
        raise ValueError("gate detail report and receipt bind different candidate artifacts")
    if summary.get("candidate_artifact_set_sha256") != canonical_payload_sha256(artifact_payload):
        raise ValueError("gate detail report candidate artifact-set hash is stale")

    requested = [str(Path(path).expanduser().resolve()) for path in root_paths]
    if len(requested) != 1:
        raise ValueError(f"gated pp-SIM current promotion requires exactly one ROOT, got {len(requested)}")
    matches = [row for row in normalized if row["family"] == family]
    if len(matches) != 1 or matches[0]["path"] != requested[0]:
        raise ValueError(
            f"production gate report does not bind the requested {family} ROOT: {requested[0]}"
        )
    root = Path(requested[0])
    if not root.exists() or file_sha256(root) != matches[0]["sha256"]:
        raise ValueError(f"requested {family} ROOT does not match the production gate SHA-256: {root}")
    return payload


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

    production_gate_report = ""
    production_gate_report_sha256 = ""
    if args.status == "current" and args.sample_key in PPG12_PRODUCTION_GATED_SAMPLE_FAMILIES:
        if not args.production_gate_report:
            print(
                "[ERROR] pp-SIM current promotion is blocked without "
                "--production-gate-report",
                file=sys.stderr,
            )
            return 3
        report_path = Path(normalize_path(args.production_gate_report))
        try:
            validate_ppg12_production_gate(report_path, args.sample_key, root_paths)
        except ValueError as exc:
            print(f"[ERROR] pp-SIM production gate rejected promotion: {exc}", file=sys.stderr)
            return 3
        production_gate_report = str(report_path)
        production_gate_report_sha256 = file_sha256(report_path)

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
        "production_gate_report": production_gate_report,
        "production_gate_report_sha256": production_gate_report_sha256,
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
    resolved_current: list[tuple[str, dict[str, Any]]] = []
    for sample_key, entry_id in registry.get("current", {}).items():
        entry = entry_by_id.get(entry_id)
        if not entry:
            print(f"[ERROR] current entry id {entry_id} for {sample_key} is missing", file=sys.stderr)
            return 2
        if entry.get("status") != "current":
            print(f"[ERROR] current entry id {entry_id} for {sample_key} has status={entry.get('status')}", file=sys.stderr)
            return 2
        if sample_key in PPG12_PRODUCTION_GATED_SAMPLE_FAMILIES:
            raw_report = entry.get("production_gate_report")
            if not isinstance(raw_report, str) or not raw_report:
                print(
                    f"[ERROR] pp-SIM current snapshot sync is blocked without a production gate report: {sample_key}",
                    file=sys.stderr,
                )
                return 3
            report_path = Path(normalize_path(raw_report))
            recorded_sha = entry.get("production_gate_report_sha256")
            if (
                not report_path.exists()
                or not isinstance(recorded_sha, str)
                or file_sha256(report_path) != recorded_sha
            ):
                print(
                    f"[ERROR] pp-SIM current snapshot sync found a stale production gate report: {sample_key}",
                    file=sys.stderr,
                )
                return 3
            try:
                validate_ppg12_production_gate(
                    report_path, sample_key, list(entry.get("root_paths", []))
                )
            except ValueError as exc:
                print(
                    f"[ERROR] pp-SIM current snapshot sync rejected promotion evidence: {exc}",
                    file=sys.stderr,
                )
                return 3
        resolved_current.append((sample_key, entry))
    for entry in registry.get("artifacts", []):
        write_entry_snapshots(args.registry, entry)
    for _, entry in resolved_current:
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
    reg.add_argument(
        "--production-gate-report",
        default="",
        help=(
            "Required for current pp_sim_photonjet_merged or "
            "pp_sim_inclusivejet_merged promotion; must bind the exact ROOT SHA-256"
        ),
    )
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
