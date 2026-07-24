#!/usr/bin/env python3
"""Build and revalidate the exact THE-134 13-row source-authority manifest.

The output is consumed directly by
``resolve_the134_full_multiview_extraction.py``.  It freezes source, period,
SI/embedded role, ownership, five-list content, and aligned tuple identity
without submitting jobs or copying payloads.

The JSON is deliberately timestamp-free and canonical: identical inputs
produce identical bytes.  ``verify`` requires the previously recorded manifest
SHA-256, then rehashes every source list and rebuilds the complete payload from
the current paths.  A mutable path, duplicate row, missing list, or authority
drift therefore fails before submission.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
from pathlib import Path
from typing import Any, Iterable, Mapping


REPO_ROOT = Path(__file__).resolve().parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.ml.contracts.the134_h70_contract import (  # noqa: E402
    BACKGROUND_SOURCES,
    SIGNAL_SOURCES,
)


SOURCE_SCHEMA = "THE134_FULL_EXTRACTION_SOURCE_AUTHORITY_V1"
READBACK_SCHEMA = "THE134_FULL_EXTRACTION_SOURCE_AUTHORITY_READBACK_V1"
AUTHORITY_STATE = "MEASURED_CANDIDATE_NOT_SCIENTIFIC_AUTHORITY"
SOURCE_OWNERSHIP_STATE = "source_role_frozen"
AUAU_PERIOD = "AUAU_RUN24"
PP_SI_DI_ROLE = "SI"
AUAU_SI_DI_ROLE = "EMBEDDED"
PP_PERIODS = ("0mrad", "1p5mrad")
LIST_ROLES = ("calo_cluster", "g4hits", "jets", "global", "mbd_epd")
LIST_FILENAMES = {
    "calo_cluster": "DST_CALO_CLUSTER.matched.list",
    "g4hits": "G4Hits.matched.list",
    "jets": "DST_JETS.matched.list",
    "global": "DST_GLOBAL.matched.list",
    "mbd_epd": "DST_MBD_EPD.matched.list",
}
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")

# row_id, system, lane, dataset, sample, source_role, minimum_bias_gate,
# photon_id_row_match
SOURCE_ROWS = (
    (
        "pp_signal_photon5",
        "pp",
        "pp_photon_sim",
        "isSim",
        "run28_photonjet5",
        "signal",
        "not_applicable",
        "newPPG12",
    ),
    (
        "pp_signal_photon10",
        "pp",
        "pp_photon_sim",
        "isSim",
        "run28_photonjet10",
        "signal",
        "not_applicable",
        "newPPG12",
    ),
    (
        "pp_signal_photon20",
        "pp",
        "pp_photon_sim",
        "isSim",
        "run28_photonjet20",
        "signal",
        "not_applicable",
        "newPPG12",
    ),
    (
        "pp_background_jet8",
        "pp",
        "pp_inclusive_sim",
        "isSimInclusive",
        "run28_jet8",
        "background",
        "not_applicable",
        "newPPG12",
    ),
    (
        "pp_background_jet12",
        "pp",
        "pp_inclusive_sim",
        "isSimInclusive",
        "run28_jet12",
        "background",
        "not_applicable",
        "newPPG12",
    ),
    (
        "pp_background_jet20",
        "pp",
        "pp_inclusive_sim",
        "isSimInclusive",
        "run28_jet20",
        "background",
        "not_applicable",
        "newPPG12",
    ),
    (
        "pp_background_jet30",
        "pp",
        "pp_inclusive_sim",
        "isSimInclusive",
        "run28_jet30",
        "background",
        "not_applicable",
        "newPPG12",
    ),
    (
        "auau_signal_photon12",
        "auau",
        "auau_photon_embedded",
        "isSimEmbedded",
        "run28_embeddedPhoton12",
        "signal",
        "required_pass",
        "auauBDTSideband",
    ),
    (
        "auau_signal_photon20",
        "auau",
        "auau_photon_embedded",
        "isSimEmbedded",
        "run28_embeddedPhoton20",
        "signal",
        "required_pass",
        "auauBDTSideband",
    ),
    (
        "auau_background_jet12",
        "auau",
        "auau_inclusive_embedded",
        "isSimEmbeddedInclusive",
        "run28_embeddedJet12",
        "background",
        "required_pass",
        "auauBDTSideband",
    ),
    (
        "auau_background_jet20",
        "auau",
        "auau_inclusive_embedded",
        "isSimEmbeddedInclusive",
        "run28_embeddedJet20",
        "background",
        "required_pass",
        "auauBDTSideband",
    ),
    (
        "auau_background_jet30",
        "auau",
        "auau_inclusive_embedded",
        "isSimEmbeddedInclusive",
        "run28_embeddedJet30",
        "background",
        "required_pass",
        "auauBDTSideband",
    ),
    (
        "auau_background_jet40",
        "auau",
        "auau_inclusive_embedded",
        "isSimEmbeddedInclusive",
        "run28_embeddedJet40",
        "background",
        "required_pass",
        "auauBDTSideband",
    ),
)


class ManifestError(ValueError):
    """Fail-closed source-authority contract violation."""


def canonical_json_bytes(payload: Any) -> bytes:
    return (
        json.dumps(
            payload,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        )
        + "\n"
    ).encode("utf-8")


def canonical_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(
            payload,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        ).encode("utf-8")
    ).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_sha256(label: str, value: object) -> str:
    text = str(value)
    if not SHA256_RE.fullmatch(text):
        raise ManifestError(f"{label} must be a lowercase 64-character SHA-256")
    return text


def require_absolute_dir(label: str, value: object) -> Path:
    path = Path(str(value)).expanduser()
    if not path.is_absolute():
        raise ManifestError(f"{label} must be an absolute path: {path}")
    if not path.is_dir():
        raise ManifestError(f"{label} is not an existing directory: {path}")
    return path


def inventory_rows() -> list[dict[str, str]]:
    fields = (
        "row_id",
        "system",
        "lane",
        "dataset",
        "sample",
        "source_role",
        "minimum_bias_gate",
        "photon_id_row_match",
    )
    return [dict(zip(fields, values)) for values in SOURCE_ROWS]


def validate_frozen_inventory(rows: Iterable[Mapping[str, str]]) -> None:
    materialized = [dict(row) for row in rows]
    row_ids = [row["row_id"] for row in materialized]
    samples = [row["sample"] for row in materialized]
    if (
        len(materialized) != 13
        or len(set(row_ids)) != 13
        or len(set(samples)) != 13
    ):
        raise ManifestError("THE-134 source inventory must be exactly 13 unique rows")
    expected_sources = {
        "pp": set(SIGNAL_SOURCES["pp"]) | set(BACKGROUND_SOURCES["pp"]),
        "auau": set(SIGNAL_SOURCES["auau"]) | set(BACKGROUND_SOURCES["auau"]),
    }
    observed_sources = {
        system: {
            row["sample"] for row in materialized if row["system"] == system
        }
        for system in ("pp", "auau")
    }
    if observed_sources != expected_sources:
        raise ManifestError(
            "THE-134 source inventory differs from the frozen ML contract: "
            f"{observed_sources}"
        )
    if "run28_jet40" in observed_sources["pp"]:
        raise ManifestError("p+p Jet40 is diagnostic-only and cannot enter training")
    family_counts = {
        (system, role): sum(
            row["system"] == system and row["source_role"] == role
            for row in materialized
        )
        for system in ("pp", "auau")
        for role in ("signal", "background")
    }
    if family_counts != {
        ("pp", "signal"): 3,
        ("pp", "background"): 4,
        ("auau", "signal"): 2,
        ("auau", "background"): 4,
    }:
        raise ManifestError(f"source-family population drift: {family_counts}")


def executable_line(raw: str) -> str | None:
    stripped = raw.strip()
    if not stripped or stripped.startswith("#"):
        return None
    return stripped


def source_period(row: Mapping[str, str], pp_period: str) -> str:
    return pp_period if row["system"] == "pp" else AUAU_PERIOD


def source_si_di_role(row: Mapping[str, str]) -> str:
    return PP_SI_DI_ROLE if row["system"] == "pp" else AUAU_SI_DI_ROLE


def inspect_source_row(
    sim_list_root: Path,
    row: Mapping[str, str],
    *,
    pp_period: str,
) -> dict[str, Any]:
    sample_root = sim_list_root / row["sample"]
    if not sample_root.is_dir():
        raise ManifestError(
            f"{row['row_id']} sample directory is missing: {sample_root}"
        )

    list_records: list[dict[str, Any]] = []
    lines_by_role: dict[str, list[str]] = {}
    for role in LIST_ROLES:
        path = sample_root / LIST_FILENAMES[role]
        if not path.is_file():
            raise ManifestError(f"{row['row_id']} {role} list is missing: {path}")
        lines = path.read_text(encoding="utf-8").splitlines()
        executable_count = sum(executable_line(line) is not None for line in lines)
        if executable_count <= 0:
            raise ManifestError(
                f"{row['row_id']} {role} list has no executable entries"
            )
        resolved = path.resolve(strict=True)
        list_records.append(
            {
                "role": role,
                "path": str(path),
                "resolved_path": str(resolved),
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
                "line_count": len(lines),
                "executable_count": executable_count,
            }
        )
        lines_by_role[role] = lines

    physical_counts = {len(lines_by_role[role]) for role in LIST_ROLES}
    if len(physical_counts) != 1:
        raise ManifestError(
            f"{row['row_id']} five-list physical line counts differ"
        )

    tuples: list[dict[str, Any]] = []
    for line_index in range(next(iter(physical_counts))):
        inputs = {
            role: executable_line(lines_by_role[role][line_index])
            for role in LIST_ROLES
        }
        active = sum(value is not None for value in inputs.values())
        if active == 0:
            continue
        if active != len(LIST_ROLES):
            raise ManifestError(
                f"{row['row_id']} has a partial five-file tuple at physical "
                f"line {line_index + 1}"
            )
        tuples.append(
            {
                "tuple_index": len(tuples),
                "physical_line": line_index + 1,
                "inputs": inputs,
            }
        )
    if not tuples:
        raise ManifestError(f"{row['row_id']} has no complete input tuples")
    if {
        record["executable_count"] for record in list_records
    } != {len(tuples)}:
        raise ManifestError(
            f"{row['row_id']} executable list counts differ from tuple count"
        )

    tuple_records_sha256 = canonical_sha256(tuples)
    tuple_input_sha256s = [
        canonical_sha256(record["inputs"]) for record in tuples
    ]
    resolver_semantic_payload = {
        "row_id": row["row_id"],
        "system": row["system"],
        "sample": row["sample"],
        "lists": [
            {
                "role": record["role"],
                "sha256": record["sha256"],
                "line_count": record["line_count"],
                "executable_count": record["executable_count"],
            }
            for record in list_records
        ],
        "tuple_count": len(tuples),
        "tuple_records_sha256": tuple_records_sha256,
    }
    return {
        **dict(row),
        "source_period": source_period(row, pp_period),
        "source_si_di_role": source_si_di_role(row),
        "source_ownership_state": SOURCE_OWNERSHIP_STATE,
        "input_lists": list_records,
        "tuple_count": len(tuples),
        "tuple_records_sha256": tuple_records_sha256,
        "full_source_manifest_sha256": canonical_sha256(
            resolver_semantic_payload
        ),
        "expected_input_count": len(tuples),
        "expected_occurrence_count": len(tuples),
        "first_tuple_sha256": canonical_sha256(tuples[0]),
        "last_tuple_sha256": canonical_sha256(tuples[-1]),
        "_tuple_input_sha256s": tuple_input_sha256s,
    }


def row_semantic_record(row: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "row_id": row["row_id"],
        "system": row["system"],
        "lane": row["lane"],
        "dataset": row["dataset"],
        "sample": row["sample"],
        "source_role": row["source_role"],
        "source_period": row["source_period"],
        "source_si_di_role": row["source_si_di_role"],
        "source_ownership_state": row["source_ownership_state"],
        "minimum_bias_gate": row["minimum_bias_gate"],
        "photon_id_row_match": row["photon_id_row_match"],
        "input_lists": [
            {
                key: record[key]
                for key in (
                    "role",
                    "sha256",
                    "size_bytes",
                    "line_count",
                    "executable_count",
                )
            }
            for record in row["input_lists"]
        ],
        "tuple_count": row["tuple_count"],
        "tuple_records_sha256": row["tuple_records_sha256"],
        "full_source_manifest_sha256": row["full_source_manifest_sha256"],
    }


def build_manifest(sim_list_root: Path, pp_period: str) -> dict[str, Any]:
    root = require_absolute_dir("sim_list_root", sim_list_root)
    if pp_period not in PP_PERIODS:
        raise ManifestError(
            f"pp_period must be one of {PP_PERIODS}, observed {pp_period!r}"
        )
    frozen = inventory_rows()
    validate_frozen_inventory(frozen)
    rows = [
        inspect_source_row(root, row, pp_period=pp_period) for row in frozen
    ]
    tuple_identity_records = []
    all_tuple_identities = []
    for row in rows:
        identities = row.pop("_tuple_input_sha256s")
        tuple_identity_records.append(
            {
                "row_id": row["row_id"],
                "tuple_input_sha256s": identities,
            }
        )
        all_tuple_identities.extend(identities)
    unique_tuple_count = len(set(all_tuple_identities))
    if unique_tuple_count != len(all_tuple_identities):
        raise ManifestError(
            "duplicate five-file source tuples exist across THE-134 rows"
        )
    global_five_tuple_check = {
        "total": len(all_tuple_identities),
        "unique": unique_tuple_count,
        "duplicate_count": len(all_tuple_identities) - unique_tuple_count,
        "tuple_identity_records_sha256": canonical_sha256(
            tuple_identity_records
        ),
    }

    resolved_paths = [
        record["resolved_path"]
        for row in rows
        for record in row["input_lists"]
    ]
    if len(resolved_paths) != len(set(resolved_paths)):
        raise ManifestError("one physical source-list path is assigned twice")

    authority = {
        "scope": "THE134_FULL_13_ROW_TRAINING_EXTRACTION",
        "row_count": 13,
        "pp_period": pp_period,
        "pp_si_di_role": PP_SI_DI_ROLE,
        "auau_period": AUAU_PERIOD,
        "auau_si_di_role": AUAU_SI_DI_ROLE,
        "source_ownership_state": SOURCE_OWNERSHIP_STATE,
        "diagnostic_sources_excluded": ["run28_jet40"],
        "scientific_completion_granted": False,
    }
    semantic_rows = [row_semantic_record(row) for row in rows]
    payload = {
        "schema": SOURCE_SCHEMA,
        "status": "PASS",
        "authority_state": AUTHORITY_STATE,
        "authority": authority,
        "sim_list_root": str(root),
        "resolved_sim_list_root": str(root.resolve(strict=True)),
        "row_count": len(rows),
        "rows": rows,
        "global_five_tuple_check": global_five_tuple_check,
        "source_rows_semantic_sha256": canonical_sha256(semantic_rows),
        "manifest_semantic_sha256": canonical_sha256(
            {
                "authority_state": AUTHORITY_STATE,
                "authority": authority,
                "rows": semantic_rows,
                "global_five_tuple_check": global_five_tuple_check,
            }
        ),
    }
    validate_manifest_payload(payload, rehash=False)
    return payload


def validate_manifest_payload(
    payload: Mapping[str, Any],
    *,
    rehash: bool,
) -> dict[str, Any]:
    required_top = {
        "schema",
        "status",
        "authority_state",
        "authority",
        "sim_list_root",
        "resolved_sim_list_root",
        "row_count",
        "rows",
        "global_five_tuple_check",
        "source_rows_semantic_sha256",
        "manifest_semantic_sha256",
    }
    if set(payload) != required_top:
        raise ManifestError(
            "source manifest top-level inventory differs: "
            f"missing={sorted(required_top - set(payload))} "
            f"extra={sorted(set(payload) - required_top)}"
        )
    if payload["schema"] != SOURCE_SCHEMA or payload["status"] != "PASS":
        raise ManifestError("source manifest schema/status is invalid")
    if payload["authority_state"] != AUTHORITY_STATE:
        raise ManifestError("source manifest cannot grant scientific authority")
    rows = payload["rows"]
    if not isinstance(rows, list):
        raise ManifestError("source manifest rows must be a list")
    row_ids = [str(row.get("row_id", "")) for row in rows if isinstance(row, dict)]
    if len(rows) != 13 or len(row_ids) != 13 or len(set(row_ids)) != 13:
        raise ManifestError("source manifest rows are missing or duplicated")
    expected = inventory_rows()
    if row_ids != [row["row_id"] for row in expected]:
        raise ManifestError("source manifest row order/identity differs")
    if int(payload["row_count"]) != 13:
        raise ManifestError("source manifest row_count must be exactly 13")
    global_check = payload["global_five_tuple_check"]
    if not isinstance(global_check, dict) or set(global_check) != {
        "total",
        "unique",
        "duplicate_count",
        "tuple_identity_records_sha256",
    }:
        raise ManifestError("global five-tuple check inventory differs")
    tuple_total = sum(int(row.get("tuple_count", 0)) for row in rows)
    if (
        global_check["total"] != tuple_total
        or global_check["unique"] != tuple_total
        or global_check["duplicate_count"] != 0
    ):
        raise ManifestError("global five-tuple uniqueness closure failed")
    require_sha256(
        "global tuple identity records SHA-256",
        global_check["tuple_identity_records_sha256"],
    )

    authority = payload["authority"]
    if not isinstance(authority, dict):
        raise ManifestError("source manifest authority must be an object")
    expected_authority = {
        "scope": "THE134_FULL_13_ROW_TRAINING_EXTRACTION",
        "row_count": 13,
        "pp_period": authority.get("pp_period"),
        "pp_si_di_role": PP_SI_DI_ROLE,
        "auau_period": AUAU_PERIOD,
        "auau_si_di_role": AUAU_SI_DI_ROLE,
        "source_ownership_state": SOURCE_OWNERSHIP_STATE,
        "diagnostic_sources_excluded": ["run28_jet40"],
        "scientific_completion_granted": False,
    }
    if authority != expected_authority or authority["pp_period"] not in PP_PERIODS:
        raise ManifestError("source manifest authority fields differ")

    semantic_rows: list[dict[str, Any]] = []
    for observed, frozen in zip(rows, expected):
        if not isinstance(observed, dict):
            raise ManifestError("source manifest row must be an object")
        for field, value in frozen.items():
            if observed.get(field) != value:
                raise ManifestError(
                    f"{frozen['row_id']} frozen field {field} differs"
                )
        if observed.get("source_period") != source_period(
            frozen, authority["pp_period"]
        ):
            raise ManifestError(f"{frozen['row_id']} source period differs")
        if observed.get("source_si_di_role") != source_si_di_role(frozen):
            raise ManifestError(f"{frozen['row_id']} SI/embedded role differs")
        if observed.get("source_ownership_state") != SOURCE_OWNERSHIP_STATE:
            raise ManifestError(f"{frozen['row_id']} ownership state differs")
        records = observed.get("input_lists")
        if not isinstance(records, list) or [
            record.get("role") for record in records if isinstance(record, dict)
        ] != list(LIST_ROLES):
            raise ManifestError(f"{frozen['row_id']} list-role inventory differs")
        for record in records:
            if not isinstance(record, dict):
                raise ManifestError(f"{frozen['row_id']} list record is invalid")
            require_sha256(
                f"{frozen['row_id']} {record.get('role')} sha256",
                record.get("sha256", ""),
            )
            if rehash:
                path = Path(str(record.get("path", "")))
                if not path.is_absolute() or not path.is_file():
                    raise ManifestError(
                        f"{frozen['row_id']} source list is missing: {path}"
                    )
                resolved = str(path.resolve(strict=True))
                if resolved != record.get("resolved_path"):
                    raise ManifestError(
                        f"{frozen['row_id']} mutable path resolution drift: {path}"
                    )
                observed_sha = sha256_file(path)
                if observed_sha != record["sha256"]:
                    raise ManifestError(
                        f"{frozen['row_id']} source list hash drift: "
                        f"role={record['role']} expected={record['sha256']} "
                        f"observed={observed_sha}"
                    )
                if path.stat().st_size != int(record["size_bytes"]):
                    raise ManifestError(
                        f"{frozen['row_id']} source list size drift: "
                        f"role={record['role']}"
                    )
        semantic_rows.append(row_semantic_record(observed))

    expected_rows_sha = canonical_sha256(semantic_rows)
    if require_sha256(
        "source_rows_semantic_sha256",
        payload["source_rows_semantic_sha256"],
    ) != expected_rows_sha:
        raise ManifestError("source row semantic digest differs")
    expected_manifest_sha = canonical_sha256(
        {
            "authority_state": AUTHORITY_STATE,
            "authority": authority,
            "rows": semantic_rows,
            "global_five_tuple_check": global_check,
        }
    )
    if require_sha256(
        "manifest_semantic_sha256", payload["manifest_semantic_sha256"]
    ) != expected_manifest_sha:
        raise ManifestError("source manifest semantic digest differs")

    if rehash:
        rebuilt = build_manifest(
            Path(str(payload["sim_list_root"])), authority["pp_period"]
        )
        if canonical_json_bytes(rebuilt) != canonical_json_bytes(dict(payload)):
            raise ManifestError(
                "source manifest current-path rebuild differs from frozen payload"
            )
    return dict(payload)


def load_manifest(path: Path, expected_sha256: str) -> dict[str, Any]:
    expected = require_sha256("expected manifest SHA-256", expected_sha256)
    if not path.is_absolute() or not path.is_file():
        raise ManifestError(f"manifest must be an existing absolute file: {path}")
    observed = sha256_file(path)
    if observed != expected:
        raise ManifestError(
            f"manifest file hash drift: expected={expected} observed={observed}"
        )
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ManifestError(f"manifest is not valid JSON: {path}") from exc
    if not isinstance(payload, dict):
        raise ManifestError("manifest must contain one JSON object")
    return validate_manifest_payload(payload, rehash=True)


def atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    if not path.is_absolute():
        raise ManifestError(f"output path must be absolute: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    temporary.write_bytes(canonical_json_bytes(payload))
    os.replace(temporary, path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    build = subparsers.add_parser("build", help="build the deterministic manifest")
    build.add_argument("--sim-list-root", type=Path, required=True)
    build.add_argument("--pp-period", choices=PP_PERIODS, required=True)
    build.add_argument("--output", type=Path, required=True)

    verify = subparsers.add_parser(
        "verify", help="rehash and rebuild one pinned manifest"
    )
    verify.add_argument("--manifest", type=Path, required=True)
    verify.add_argument("--expected-sha256", required=True)
    verify.add_argument("--output-json", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        if args.command == "build":
            payload = build_manifest(args.sim_list_root, args.pp_period)
            atomic_write_json(args.output, payload)
            summary = {
                "schema": SOURCE_SCHEMA,
                "status": "PASS",
                "authority_state": AUTHORITY_STATE,
                "path": str(args.output),
                "sha256": sha256_file(args.output),
                "row_count": 13,
                "manifest_semantic_sha256": payload[
                    "manifest_semantic_sha256"
                ],
            }
        else:
            payload = load_manifest(args.manifest, args.expected_sha256)
            summary = {
                "schema": READBACK_SCHEMA,
                "status": "PASS",
                "authority_state": AUTHORITY_STATE,
                "manifest": str(args.manifest),
                "manifest_sha256": args.expected_sha256,
                "row_count": 13,
                "source_rows_rehashed": 13,
                "manifest_semantic_sha256": payload[
                    "manifest_semantic_sha256"
                ],
            }
            if args.output_json is not None:
                atomic_write_json(args.output_json, summary)
        print(json.dumps(summary, sort_keys=True))
        return 0
    except ManifestError as exc:
        print(f"THE134_SOURCE_AUTHORITY_ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
