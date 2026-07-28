#!/usr/bin/env python3
"""Certify the bounded THE-134 direct/sidecar-only differential canary.

The canary owns exactly p+p Jet8 and Au+Au embedded Jet12, with one direct and
one writer arm per system.  Analysis ROOT histograms must remain exactly
neutral.  Writer analysis files must contain the frozen validation-only
markers and no serialized ReplayFoundation directory.  The compact
RJPhotonTrainingViewV1 sidecars use their artifact-specific semantic health
contract rather than the analysis ROOT 50 kB size rule.  A bounded replacement
mode may certify only a corrected p+p pair and then bind it to the independently
validated, preserved Au+Au row through distinct immutable receipts.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import re
import shlex
import sys
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import uproot


HERE = Path(__file__).resolve()
REPOSITORY = HERE.parents[3]
PREPARER_PATH = (
    REPOSITORY / "scripts" / "ml" / "training" / "prepare_the134_h70_matrix.py"
)

CERTIFICATE_SCHEMA = "THE134_SIDECAR_DIFFERENTIAL_CERTIFICATE_V1"
COMPONENT_CERTIFICATE_SCHEMA = (
    "THE134_SIDECAR_DIFFERENTIAL_COMPONENT_CERTIFICATE_V1"
)
AGGREGATE_CERTIFICATE_SCHEMA = (
    "THE134_SIDECAR_DIFFERENTIAL_AGGREGATE_CERTIFICATE_V1"
)
TERMINAL_GATE_SCHEMA = "THE134_PP_REPLACEMENT_TERMINAL_GATE_V1"
PROFILE = "THE134_MULTIVIEW_SIDECAR_ONLY_V1"
SIDECAR_NAME = "RJPhotonTrainingViewV1.root"
ANALYSIS_MIN_BYTES = 50_000
SIDECAR_GROSS_TRUNCATION_BYTES = 4_096
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
HISTOGRAM_PREFIXES = ("TH1", "TH2", "TH3", "TProfile")
FORBIDDEN_REPLAY_PREFIX = "ReplayFoundationV1/"

ROWS = (
    {
        "system": "pp",
        "lane": "pp_inclusive_sim",
        "dataset": "isSimInclusive",
        "sample": "run28_jet8",
        "period": "0mrad",
        "si_di_role": "SI",
        "source": "run28_jet8",
    },
    {
        "system": "auau",
        "lane": "auau_inclusive_embedded",
        "dataset": "isSimEmbeddedInclusive",
        "sample": "run28_embeddedJet12",
        "period": "AUAU_RUN24",
        "si_di_role": "EMBEDDED",
        "source": "run28_embeddedJet12",
    },
)
PP_ROWS = tuple(row for row in ROWS if row["system"] == "pp")
AUAU_ROWS = tuple(row for row in ROWS if row["system"] == "auau")

MARKER_ENV = {
    "rj_replay_schema_sha256": "schema_sha",
    "rj_replay_semantic_sha256": "semantic_sha",
    "rj_replay_code_sha256": "code_sha256",
}
FIXED_MARKERS = {
    "rj_the134_multiview_sidecar_only_v1": "1",
    "rj_replay_transaction_state": "CONSTRUCTED_AND_VALIDATED",
    "rj_replay_serialization_state": "DISABLED",
    "rj_replay_cache_applicability": "NOT_APPLICABLE",
}


class ValidationError(RuntimeError):
    """Fail-closed differential-canary validation error."""


def load_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ValidationError(f"cannot load validator module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_sha256(tag: str, row: dict[str, str]) -> str:
    payload = (
        f"{tag}|{row['lane']}|{row['dataset']}|{row['sample']}|"
        "accepted-canary-source-v1"
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def parse_preflight(
    path: Path,
    ownership_rows: Iterable[dict[str, str]] = ROWS,
) -> dict[str, Any]:
    if not path.is_file():
        raise ValidationError(f"preflight receipt is missing: {path}")
    fields: dict[str, str] = {}
    pinned: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        if "=" in raw and not raw.startswith((" ", "\t")):
            key, value = raw.split("=", 1)
            if key in fields:
                raise ValidationError(
                    f"preflight receipt contains duplicate field: {key}"
                )
            fields[key] = value
            continue
        parts = raw.split()
        if len(parts) == 2 and SHA256_RE.fullmatch(parts[0]):
            basename = Path(parts[1]).name
            if basename in pinned:
                raise ValidationError(
                    f"preflight receipt contains duplicate pinned basename: {basename}"
                )
            pinned[basename] = parts[0]
    required = (
        "tag",
        "base",
        "code_sha256",
        "schema_sha",
        "semantic_sha",
        "only_keys",
    )
    missing = [name for name in required if not fields.get(name)]
    if missing:
        raise ValidationError(f"preflight receipt lacks fields: {missing}")
    for name in ("code_sha256", "schema_sha", "semantic_sha"):
        if SHA256_RE.fullmatch(fields[name]) is None:
            raise ValidationError(f"preflight {name} is not SHA-256")
    if not Path(fields["base"]).is_absolute():
        raise ValidationError("preflight base is not absolute")
    owned_rows = tuple(ownership_rows)
    expected_only_keys = {
        f"{arm}:{row['lane']}:{row['sample']}"
        for row in owned_rows
        for arm in ("direct", "writer")
    }
    observed_only_keys = fields["only_keys"].split(",")
    if (
        any(not key for key in observed_only_keys)
        or len(observed_only_keys) != len(expected_only_keys)
        or len(set(observed_only_keys)) != len(observed_only_keys)
        or set(observed_only_keys) != expected_only_keys
    ):
        raise ValidationError("preflight exact row ownership differs")
    if fields.get("writer_extra_common_template", "").find(
        "RJ_THE134_MULTIVIEW_SIDECAR_ONLY_V1=1"
    ) < 0:
        raise ValidationError("preflight lacks the sidecar-only writer profile")
    owned_systems = {row["system"] for row in owned_rows}
    if "pp" in owned_systems and fields.get("extra_pp_template", "").find(
        "RJ_PP_PHOTONID_EXTRACT_ONLY=1"
    ) < 0:
        raise ValidationError("preflight lacks the shared p+p extraction profile")
    if "auau" in owned_systems and fields.get("extra_auau_template", "").find(
        "RJ_AUAU_BDT_EXTRACT_ONLY=1"
    ) < 0:
        raise ValidationError("preflight lacks the shared Au+Au extraction profile")
    if fields.get("writer_extra_pp_template") or fields.get(
        "writer_extra_auau_template"
    ):
        raise ValidationError("preflight contains unexpected system-specific writer extensions")
    return {
        "fields": fields,
        "pinned": pinned,
        "path": str(path.resolve()),
        "sha256": sha256_file(path),
    }


def preflight_source_sha256(
    preflight: dict[str, Any],
    row: dict[str, str],
) -> str:
    override = preflight["fields"].get("source_sha_override", "")
    if override:
        if SHA256_RE.fullmatch(override) is None:
            raise ValidationError("preflight source_sha_override is not SHA-256")
        return override
    return source_sha256(preflight["fields"]["tag"], row)


def preflight_provenance(preflight: dict[str, Any]) -> dict[str, str]:
    return {
        "tag": preflight["fields"]["tag"],
        "code_sha256": preflight["fields"]["code_sha256"],
        "schema_sha256": preflight["fields"]["schema_sha"],
        "semantic_sha256": preflight["fields"]["semantic_sha"],
        "only_keys": preflight["fields"]["only_keys"],
    }


def parse_terminal_gate_receipt(
    path: Path,
    *,
    preflight: dict[str, Any],
    output_root: Path,
) -> dict[str, Any]:
    resolved = path.resolve()
    if not resolved.is_file():
        raise ValidationError(
            f"p+p terminal-gate receipt is missing: {resolved}"
        )
    expected_evidence_root = Path(preflight["path"]).resolve().parent
    if resolved.parent != expected_evidence_root:
        raise ValidationError(
            "p+p terminal-gate receipt is outside its preflight namespace"
        )
    try:
        payload = json.loads(resolved.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValidationError(
            f"p+p terminal-gate receipt is unreadable: {resolved}"
        ) from exc
    expected_fields = {
        "schema",
        "tag",
        "output_root",
        "preflight_receipt",
        "preflight_receipt_sha256",
        "initial_queue_tsv",
        "initial_queue_tsv_sha256",
        "row_count",
        "rows",
    }
    if not isinstance(payload, dict) or set(payload) != expected_fields:
        raise ValidationError("p+p terminal-gate field inventory differs")
    if payload.get("schema") != TERMINAL_GATE_SCHEMA:
        raise ValidationError("p+p terminal-gate receipt schema differs")
    if payload.get("tag") != preflight["fields"]["tag"]:
        raise ValidationError("p+p terminal-gate tag differs")
    if Path(str(payload.get("output_root", ""))).resolve() != output_root.resolve():
        raise ValidationError("p+p terminal-gate output root differs")
    preflight_path = Path(str(payload.get("preflight_receipt", ""))).resolve()
    preflight_sha = str(payload.get("preflight_receipt_sha256", ""))
    if (
        preflight_path != Path(preflight["path"]).resolve()
        or SHA256_RE.fullmatch(preflight_sha) is None
        or preflight_sha != preflight["sha256"]
        or not preflight_path.is_file()
        or sha256_file(preflight_path) != preflight_sha
    ):
        raise ValidationError("p+p terminal-gate preflight binding differs")
    initial_queue = Path(str(payload.get("initial_queue_tsv", ""))).resolve()
    initial_queue_sha = str(payload.get("initial_queue_tsv_sha256", ""))
    if (
        not initial_queue.is_file()
        or SHA256_RE.fullmatch(initial_queue_sha) is None
        or sha256_file(initial_queue) != initial_queue_sha
    ):
        raise ValidationError("p+p terminal-gate initial queue binding differs")
    if (
        initial_queue.parent != expected_evidence_root
        or initial_queue.name != "initial_queue.tsv"
    ):
        raise ValidationError(
            "p+p terminal-gate initial queue is outside its evidence namespace"
        )
    queue_bindings: set[tuple[str, str]] = set()
    queue_output_roots = {
        str(output_root.resolve()),
        str(Path(preflight["fields"]["base"])),
        str(Path(preflight["fields"]["base"]).resolve()),
    }
    queue_lines = initial_queue.read_text(encoding="utf-8").splitlines()
    if len(queue_lines) != 2:
        raise ValidationError("p+p bound initial queue cardinality differs")
    for raw in queue_lines:
        parts = raw.split(maxsplit=3)
        if (
            len(parts) != 4
            or not parts[0].isdigit()
            or not parts[1].isdigit()
            or not parts[2].isdigit()
        ):
            raise ValidationError("p+p bound initial queue row is malformed")
        cluster_proc = f"{int(parts[0])}.{int(parts[1])}"
        args = parts[3]
        try:
            arg_tokens = shlex.split(args)
        except ValueError as exc:
            raise ValidationError(
                "p+p bound initial queue arguments are malformed"
            ) from exc
        expected_role_tokens = {
            role: {
                (
                    f"{queue_root.rstrip('/')}/{role}/"
                    f"{PP_ROWS[0]['lane']}/{PP_ROWS[0]['sample']}"
                )
                for queue_root in queue_output_roots
            }
            for role in ("direct", "writer")
        }
        token_matches = [
            (index, role)
            for index, token in enumerate(arg_tokens)
            for role, expected_tokens in expected_role_tokens.items()
            if token in expected_tokens
        ]
        if (
            len(token_matches) != 1
            or token_matches[0][0] != len(arg_tokens) - 1
        ):
            raise ValidationError(
                "p+p bound initial queue row lacks one exact final "
                "direct/writer destination token"
            )
        role = token_matches[0][1]
        binding = (cluster_proc, role)
        if binding in queue_bindings or any(
            existing[0] == cluster_proc or existing[1] == role
            for existing in queue_bindings
        ):
            raise ValidationError(
                "p+p bound initial queue identities or roles are duplicated"
            )
        queue_bindings.add(binding)
    rows = payload.get("rows")
    if (
        payload.get("row_count") != 2
        or not isinstance(rows, list)
        or len(rows) != 2
    ):
        raise ValidationError("p+p terminal-gate cardinality differs")
    identities: set[str] = set()
    roles: set[str] = set()
    receipt_bindings: set[tuple[str, str]] = set()
    for row in rows:
        if not isinstance(row, dict):
            raise ValidationError("p+p terminal-gate row is not an object")
        if set(row) != {
            "cluster_id",
            "proc_id",
            "cluster_proc",
            "role",
            "lane",
            "sample",
            "job_status",
            "exit_code",
        }:
            raise ValidationError("p+p terminal-gate row field inventory differs")
        cluster_proc = str(row.get("cluster_proc", ""))
        cluster_id = row.get("cluster_id")
        proc_id = row.get("proc_id")
        role = str(row.get("role", ""))
        if (
            type(cluster_id) is not int
            or type(proc_id) is not int
            or cluster_proc != f"{cluster_id}.{proc_id}"
            or cluster_proc in identities
        ):
            raise ValidationError(
                "p+p terminal-gate Condor identities are malformed or duplicated"
            )
        if (
            role not in {"direct", "writer"}
            or role in roles
            or row.get("lane") != PP_ROWS[0]["lane"]
            or row.get("sample") != PP_ROWS[0]["sample"]
            or row.get("job_status") != 4
            or row.get("exit_code") != 0
        ):
            raise ValidationError("p+p terminal-gate role or terminal state differs")
        identities.add(cluster_proc)
        roles.add(role)
        receipt_bindings.add((cluster_proc, role))
    if roles != {"direct", "writer"}:
        raise ValidationError("p+p terminal-gate exact direct/writer roles differ")
    if receipt_bindings != queue_bindings:
        raise ValidationError(
            "p+p terminal-gate rows disagree with the bound initial queue"
        )
    return {
        "path": str(resolved),
        "sha256": sha256_file(resolved),
        "initial_queue_tsv": str(initial_queue),
        "initial_queue_tsv_sha256": initial_queue_sha,
        "row_count": 2,
        "cluster_procs": sorted(identities),
        "roles": sorted(roles),
        "queue_bindings": [
            {"cluster_proc": cluster_proc, "role": role}
            for cluster_proc, role in sorted(queue_bindings)
        ],
    }


def root_health(path: Path, minimum_bytes: int) -> dict[str, Any]:
    result = {
        "path": str(path),
        "bytes": path.stat().st_size if path.is_file() else 0,
        "sha256": sha256_file(path) if path.is_file() else "",
        "readable": False,
        "zombie": None,
        "recovered": None,
    }
    if not path.is_file() or result["bytes"] < minimum_bytes:
        return result
    try:
        import ROOT

        ROOT.gROOT.SetBatch(True)
        handle = ROOT.TFile.Open(str(path), "READ")
        if not handle:
            return result
        result["zombie"] = bool(handle.IsZombie())
        result["recovered"] = bool(handle.TestBit(ROOT.TFile.kRecovered))
        result["readable"] = not result["zombie"] and not result["recovered"]
        handle.Close()
    except Exception as exc:  # pragma: no cover - environment detail
        result["health_error"] = str(exc)
    return result


def discover_row_files(
    output_root: Path, row: dict[str, str]
) -> dict[str, Any]:
    roots: dict[str, Any] = {}
    for arm in ("direct", "writer"):
        base = output_root / arm / row["lane"] / row["sample"]
        files = sorted(path for path in base.rglob("*.root") if path.is_file())
        sidecars = [path for path in files if path.name == SIDECAR_NAME]
        analyses = [path for path in files if path.name != SIDECAR_NAME]
        roots[arm] = {
            "base": base,
            "analysis": analyses,
            "sidecars": sidecars,
        }
    if roots["direct"]["sidecars"]:
        raise ValidationError(
            f"{row['system']} direct arm unexpectedly contains a sidecar"
        )
    if len(roots["writer"]["sidecars"]) != 1:
        raise ValidationError(
            f"{row['system']} writer must contain exactly one sidecar"
        )
    if not roots["direct"]["analysis"] or not roots["writer"]["analysis"]:
        raise ValidationError(f"{row['system']} analysis ROOT coverage is empty")
    direct_relative = {
        path.relative_to(roots["direct"]["base"]): path
        for path in roots["direct"]["analysis"]
    }
    writer_relative = {
        path.relative_to(roots["writer"]["base"]): path
        for path in roots["writer"]["analysis"]
    }
    if set(direct_relative) != set(writer_relative):
        raise ValidationError(
            f"{row['system']} direct/writer analysis inventory differs"
        )
    roots["pairs"] = [
        (relative, direct_relative[relative], writer_relative[relative])
        for relative in sorted(direct_relative, key=str)
    ]
    return roots


def histogram_keys(root_file: uproot.ReadOnlyDirectory) -> dict[str, str]:
    return {
        str(key): str(class_name)
        for key, class_name in root_file.classnames(
            recursive=True, cycle=False
        ).items()
        if str(class_name).startswith(HISTOGRAM_PREFIXES)
    }


def equal_array(left: Any, right: Any) -> bool:
    left_array = np.asarray(left)
    right_array = np.asarray(right)
    if left_array.shape != right_array.shape:
        return False
    if left_array.dtype.kind in "fc" or right_array.dtype.kind in "fc":
        return bool(
            np.all(
                (left_array == right_array)
                | (np.isnan(left_array) & np.isnan(right_array))
            )
        )
    return bool(np.array_equal(left_array, right_array))


def histogram_signature(histogram: Any) -> dict[str, Any]:
    axes = []
    for axis in histogram.axes:
        axes.append(
            {
                "edges": np.asarray(axis.edges(flow=True)),
                "label": str(getattr(axis, "label", "")),
            }
        )
    signature: dict[str, Any] = {
        "class": str(histogram.classname),
        "title": str(histogram.title),
        "axes": axes,
        "values": np.asarray(histogram.values(flow=True)),
        "errors": np.asarray(histogram.errors(flow=True)),
        "variances": np.asarray(histogram.variances(flow=True)),
    }
    for name in (
        "fEntries",
        "fTsumw",
        "fTsumw2",
        "fTsumwx",
        "fTsumwx2",
        "fTsumwy",
        "fTsumwy2",
        "fTsumwxy",
    ):
        try:
            signature[name] = histogram.member(name)
        except (KeyError, AttributeError):
            pass
    try:
        signature["fSumw2"] = np.asarray(histogram.member("fSumw2"))
    except (KeyError, AttributeError, TypeError):
        signature["fSumw2"] = np.asarray([])
    return signature


def compare_signatures(left: dict[str, Any], right: dict[str, Any]) -> list[str]:
    failures: list[str] = []
    if set(left) != set(right):
        return ["signature_key_inventory"]
    for key in sorted(left):
        lvalue, rvalue = left[key], right[key]
        if key == "axes":
            if len(lvalue) != len(rvalue):
                failures.append("axis_count")
                continue
            for index, (laxis, raxis) in enumerate(zip(lvalue, rvalue)):
                if laxis["label"] != raxis["label"]:
                    failures.append(f"axis_{index}_label")
                if not equal_array(laxis["edges"], raxis["edges"]):
                    failures.append(f"axis_{index}_edges")
        elif isinstance(lvalue, np.ndarray) or isinstance(rvalue, np.ndarray):
            if not equal_array(lvalue, rvalue):
                failures.append(key)
        elif isinstance(lvalue, float) or isinstance(rvalue, float):
            if not (
                lvalue == rvalue
                or (
                    isinstance(lvalue, float)
                    and isinstance(rvalue, float)
                    and math.isnan(lvalue)
                    and math.isnan(rvalue)
                )
            ):
                failures.append(key)
        elif lvalue != rvalue:
            failures.append(key)
    return failures


def compare_histograms(direct: Path, writer: Path) -> dict[str, Any]:
    failures: list[str] = []
    with uproot.open(direct) as direct_file, uproot.open(writer) as writer_file:
        direct_keys = histogram_keys(direct_file)
        writer_keys = histogram_keys(writer_file)
        if direct_keys != writer_keys:
            failures.append("histogram_key_or_class_inventory")
        for key in sorted(set(direct_keys) & set(writer_keys)):
            mismatches = compare_signatures(
                histogram_signature(direct_file[key]),
                histogram_signature(writer_file[key]),
            )
            failures.extend(f"{key}:{item}" for item in mismatches)
    return {
        "status": "PASS" if not failures else "FAIL",
        "histograms": len(direct_keys),
        "failures": failures,
    }


def read_named(root_file: uproot.ReadOnlyDirectory, name: str) -> str:
    if name not in root_file:
        raise ValidationError(f"writer analysis ROOT lacks marker {name}")
    try:
        return str(root_file[name].member("fTitle"))
    except (KeyError, AttributeError) as exc:
        raise ValidationError(f"marker is not a TNamed: {name}") from exc


def validate_writer_markers(
    path: Path,
    row: dict[str, str],
    preflight: dict[str, Any],
) -> dict[str, Any]:
    expected = dict(FIXED_MARKERS)
    fields = preflight["fields"]
    for marker, field in MARKER_ENV.items():
        expected[marker] = fields[field]
    expected["rj_replay_source_sha256"] = preflight_source_sha256(
        preflight, row
    )
    pinned = preflight["pinned"]
    expected["rj_replay_config_sha256"] = pinned[
        (
            "analysis_config_the119_pp_replay_foundation.yaml"
            if row["system"] == "pp"
            else "analysis_config_the112_auau_combined_bdt_triplet.yaml"
        )
    ]
    model_candidates = [
        (name, value)
        for name, value in pinned.items()
        if (
            ("pp_tight_bdt" in name and row["system"] == "pp")
            or ("auau_tight_bdt" in name and row["system"] == "auau")
        )
    ]
    if len(model_candidates) != 1:
        raise ValidationError(
            f"{row['system']} preflight model identity is ambiguous"
        )
    expected["rj_replay_model_sha256"] = model_candidates[0][1]
    observed: dict[str, str] = {}
    with uproot.open(path) as root_file:
        classnames = root_file.classnames(recursive=True, cycle=False)
        forbidden = [
            str(key)
            for key in classnames
            if str(key).startswith(FORBIDDEN_REPLAY_PREFIX)
        ]
        if forbidden:
            raise ValidationError(
                f"sidecar-only writer serialized replay payload: {forbidden[:3]}"
            )
        for name, value in expected.items():
            observed[name] = read_named(root_file, name)
            if observed[name] != value:
                raise ValidationError(
                    f"{row['system']} marker {name} differs"
                )
    return {"status": "PASS", "expected": expected, "observed": observed}


def validate_sidecar(
    path: Path,
    row: dict[str, str],
    preflight: dict[str, Any],
    preparer: Any,
) -> dict[str, Any]:
    arrays, branches, metadata = preparer.read_tree(path, preparer.TREE_NAME)
    source_hash = preflight_source_sha256(preflight, row)
    config_name = (
        "analysis_config_the119_pp_replay_foundation.yaml"
        if row["system"] == "pp"
        else "analysis_config_the112_auau_combined_bdt_triplet.yaml"
    )
    entries = int(metadata["tree_num_entries"])
    input_uri_values = {
        str(value) for value in arrays["input_uri_sha256"].tolist()
    }
    input_file_values = {
        str(value) for value in arrays["input_file_sha256"].tolist()
    }
    if entries:
        external = {
            "input_uri_sha256": next(iter(input_uri_values), ""),
            "input_file_sha256": next(iter(input_file_values), ""),
            "source_manifest_sha256": source_hash,
            "config_sha256": preflight["pinned"][config_name],
            "code_sha256": preflight["fields"]["code_sha256"],
        }
        failures = preparer.validate_artifact_metadata(
            metadata, arrays, external=external
        )
        _selected, _vectors, report = preparer.validate_file_arrays(
            arrays,
            system=row["system"],
            source=row["source"],
            view_name="H70",
        )
        failures.extend(report["failures"])
    else:
        # Sparse canaries may legitimately produce a complete zero-entry
        # sidecar.  Validate its frozen metadata and exact empty-tree
        # completion state without fabricating source or candidate rows.
        failures = preparer.validate_artifact_metadata(
            metadata, arrays, external=None
        )
        report = {
            "candidates": 0,
            "definition_counts": {},
            "failures": [],
        }
        expected_empty_metadata = {
            "source_manifest_sha256": source_hash,
            "config_sha256": preflight["pinned"][config_name],
            "code_sha256": preflight["fields"]["code_sha256"],
            "rj_photon_training_entries": "0",
        }
        for name, expected in expected_empty_metadata.items():
            if str(metadata.get(name, "")) != expected:
                failures.append(f"valid_empty_metadata_{name}_mismatch")
    expected_constants = {
        "source_lane": row["lane"],
        "source_dataset": row["dataset"],
        "source_sample": row["sample"],
        "source_period": row["period"],
        "source_si_di_role": row["si_di_role"],
        "source_ownership_state": "source_role_frozen",
        "source_manifest_sha256": source_hash,
    }
    for branch, expected in expected_constants.items():
        observed = {str(value) for value in arrays[branch].tolist()}
        if observed and observed != {expected}:
            failures.append(f"{branch}_mismatch")
    if len(input_uri_values) > 1 or len(input_file_values) > 1:
        failures.append("input_hash_population_not_constant")
    return {
        "status": "PASS" if not failures else "FAIL",
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "branches": len(branches),
        "entries": entries,
        "candidates": int(report["candidates"]),
        "definition_counts": report["definition_counts"],
        "failures": failures,
    }


def require_healthy(
    paths: Iterable[Path], minimum_bytes: int
) -> list[dict[str, Any]]:
    reports = [root_health(path, minimum_bytes) for path in paths]
    failures = [report for report in reports if not report["readable"]]
    if failures:
        raise ValidationError(f"ROOT artifact health failed: {failures}")
    return reports


def validate_rows(
    output_root: Path,
    preflight: dict[str, Any],
    selected_rows: Iterable[dict[str, str]],
    preparer: Any,
) -> list[dict[str, Any]]:
    row_reports: list[dict[str, Any]] = []
    for row in selected_rows:
        discovered = discover_row_files(output_root, row)
        analysis_paths = [
            path
            for _relative, direct, writer in discovered["pairs"]
            for path in (direct, writer)
        ]
        analysis_health = require_healthy(analysis_paths, ANALYSIS_MIN_BYTES)
        sidecar = discovered["writer"]["sidecars"][0]
        sidecar_health = require_healthy(
            [sidecar], SIDECAR_GROSS_TRUNCATION_BYTES
        )[0]
        comparisons = []
        marker_reports = []
        for relative, direct, writer in discovered["pairs"]:
            comparison = compare_histograms(direct, writer)
            comparison["relative_path"] = str(relative)
            comparison["direct"] = str(direct)
            comparison["writer"] = str(writer)
            comparisons.append(comparison)
            marker_reports.append(
                validate_writer_markers(writer, row, preflight)
            )
        sidecar_report = validate_sidecar(
            sidecar, row, preflight, preparer
        )
        failures = [
            failure
            for comparison in comparisons
            for failure in comparison["failures"]
        ]
        failures.extend(sidecar_report["failures"])
        row_reports.append(
            {
                "system": row["system"],
                "lane": row["lane"],
                "sample": row["sample"],
                "status": "PASS" if not failures else "FAIL",
                "analysis_health": analysis_health,
                "sidecar_health": sidecar_health,
                "histogram_comparisons": comparisons,
                "writer_markers": marker_reports,
                "sidecar": sidecar_report,
                "failures": failures,
            }
        )
    return row_reports


def require_exact_component_root_inventory(
    output_root: Path,
    row_reports: Iterable[dict[str, Any]],
) -> None:
    expected: set[Path] = set()
    for row in row_reports:
        expected.update(
            Path(str(report["path"])).resolve()
            for report in row.get("analysis_health", [])
        )
        sidecar = row.get("sidecar_health")
        if isinstance(sidecar, dict):
            expected.add(Path(str(sidecar["path"])).resolve())
    observed = {
        path.resolve()
        for path in output_root.rglob("*.root")
        if path.is_file()
    }
    if observed != expected:
        unexpected = sorted(str(path) for path in observed - expected)
        missing = sorted(str(path) for path in expected - observed)
        raise ValidationError(
            "p+p component exact ROOT namespace differs: "
            f"unexpected={unexpected[:3]} missing={missing[:3]}"
        )


def certify_single_receipt(
    output_root: Path,
    preflight_receipt: Path,
    *,
    selected_rows: Iterable[dict[str, str]],
    ownership_rows: Iterable[dict[str, str]],
    schema: str,
    component_system: str | None = None,
    terminal_gate_receipt: Path | None = None,
) -> dict[str, Any]:
    resolved_output_root = output_root.resolve()
    resolved_receipt = preflight_receipt.resolve()
    selected = tuple(selected_rows)
    owned = tuple(ownership_rows)
    preflight = parse_preflight(resolved_receipt, owned)
    if Path(preflight["fields"]["base"]).resolve() != resolved_output_root:
        raise ValidationError(
            "single-receipt output root disagrees with preflight base"
        )
    terminal_gate = None
    if component_system == "pp":
        if terminal_gate_receipt is None:
            raise ValidationError(
                "p+p component certification requires a terminal-gate receipt"
            )
        terminal_gate = parse_terminal_gate_receipt(
            terminal_gate_receipt,
            preflight=preflight,
            output_root=resolved_output_root,
        )
    elif terminal_gate_receipt is not None:
        raise ValidationError(
            "terminal-gate receipt is only valid for p+p component certification"
        )
    preparer = load_module("the134_sidecar_differential_preparer", PREPARER_PATH)
    row_reports = validate_rows(
        resolved_output_root, preflight, selected, preparer
    )
    if component_system == "pp":
        require_exact_component_root_inventory(
            resolved_output_root, row_reports
        )
    failures = [
        f"{row['system']}:{failure}"
        for row in row_reports
        for failure in row["failures"]
    ]
    payload = {
        "schema": schema,
        "status": "PASS" if not failures else "FAIL",
        "artifact_profile": PROFILE,
        "output_root": str(resolved_output_root),
        "preflight_receipt": str(resolved_receipt),
        "preflight_receipt_sha256": preflight["sha256"],
        "receipt_ownership_systems": sorted(
            {row["system"] for row in owned}
        ),
        "provenance": preflight_provenance(preflight),
        "analysis_minimum_bytes": ANALYSIS_MIN_BYTES,
        "sidecar_minimum_bytes_mode": "ARTIFACT_SPECIFIC_GROSS_TRUNCATION_ONLY",
        "cache_replay_applicability": "NOT_APPLICABLE",
        "full_training_authority": 0,
        "full_extraction_authority": False,
        "broad_production_authority": False,
        "rows": row_reports,
        "failures": failures,
    }
    if component_system is not None:
        payload["component_system"] = component_system
    if terminal_gate is not None:
        payload["terminal_gate"] = terminal_gate
    return payload


def verify_component_artifact_bindings(payload: dict[str, Any]) -> None:
    output_root = Path(str(payload.get("output_root", ""))).resolve()
    if not output_root.is_dir():
        raise ValidationError(
            f"component output root is unavailable: {output_root}"
        )
    for row in payload["rows"]:
        analysis_health = row.get("analysis_health")
        sidecar_health = row.get("sidecar_health")
        if (
            not isinstance(analysis_health, list)
            or not analysis_health
            or not isinstance(sidecar_health, dict)
        ):
            raise ValidationError("component certificate lacks artifact bindings")
        report_groups = [
            (analysis_health, ANALYSIS_MIN_BYTES),
            ([sidecar_health], SIDECAR_GROSS_TRUNCATION_BYTES),
        ]
        for reports, minimum_bytes in report_groups:
            for report in reports:
                path = Path(str(report.get("path", "")))
                try:
                    path.resolve().relative_to(output_root)
                except ValueError as exc:
                    raise ValidationError(
                        f"component artifact is outside its output root: {path}"
                    ) from exc
                expected = str(report.get("sha256", ""))
                if (
                    not path.is_file()
                    or SHA256_RE.fullmatch(expected) is None
                    or report.get("readable") is not True
                    or report.get("zombie") is not False
                    or report.get("recovered") is not False
                    or int(report.get("bytes", 0)) < minimum_bytes
                ):
                    raise ValidationError(
                        f"component artifact binding is incomplete: {path}"
                    )
                if sha256_file(path) != expected:
                    raise ValidationError(
                        f"component artifact hash drift: {path}"
                    )


def verify_component_receipt_binding(
    payload: dict[str, Any],
    ownership_rows: Iterable[dict[str, str]],
) -> None:
    receipt = Path(str(payload.get("preflight_receipt", "")))
    expected_receipt_sha = str(payload.get("preflight_receipt_sha256", ""))
    if (
        not receipt.is_file()
        or SHA256_RE.fullmatch(expected_receipt_sha) is None
        or sha256_file(receipt) != expected_receipt_sha
    ):
        raise ValidationError("component preflight receipt binding differs")
    parsed = parse_preflight(receipt, ownership_rows)
    if (
        Path(str(payload.get("output_root", ""))).resolve()
        != Path(parsed["fields"]["base"]).resolve()
    ):
        raise ValidationError(
            "component output root disagrees with its preflight base"
        )
    if payload.get("provenance") != preflight_provenance(parsed):
        raise ValidationError(
            "component provenance disagrees with its preflight receipt"
        )


def verify_pp_terminal_binding(payload: dict[str, Any]) -> None:
    terminal_gate = payload.get("terminal_gate")
    if not isinstance(terminal_gate, dict):
        raise ValidationError("p+p component lacks terminal-gate binding")
    preflight = parse_preflight(
        Path(str(payload.get("preflight_receipt", ""))),
        PP_ROWS,
    )
    regenerated = parse_terminal_gate_receipt(
        Path(str(terminal_gate.get("path", ""))),
        preflight=preflight,
        output_root=Path(str(payload.get("output_root", ""))),
    )
    if terminal_gate != regenerated:
        raise ValidationError("p+p component terminal-gate binding differs")


def verify_component_shape(
    payload: dict[str, Any],
    *,
    expected_system: str,
    expected_receipt_systems: list[str],
    expected_row: dict[str, str],
) -> None:
    if payload.get("schema") != COMPONENT_CERTIFICATE_SCHEMA:
        raise ValidationError(f"{expected_system} component certificate schema differs")
    if payload.get("status") != "PASS":
        raise ValidationError(f"{expected_system} component certificate is not PASS")
    if payload.get("artifact_profile") != PROFILE:
        raise ValidationError(f"{expected_system} component artifact profile differs")
    if payload.get("component_system") != expected_system:
        raise ValidationError(f"component certificate is not {expected_system}")
    if payload.get("receipt_ownership_systems") != expected_receipt_systems:
        raise ValidationError(
            f"{expected_system} component receipt ownership differs"
        )
    rows = payload.get("rows")
    if (
        not isinstance(rows, list)
        or len(rows) != 1
        or rows[0].get("system") != expected_system
        or rows[0].get("lane") != expected_row["lane"]
        or rows[0].get("sample") != expected_row["sample"]
        or rows[0].get("status") != "PASS"
        or rows[0].get("failures")
        or payload.get("failures")
    ):
        raise ValidationError(
            f"{expected_system} component row population is invalid"
        )
    if (
        payload.get("full_training_authority") != 0
        or payload.get("full_extraction_authority") is not False
        or payload.get("broad_production_authority") is not False
    ):
        raise ValidationError(
            f"{expected_system} component authority boundary differs"
        )
    provenance = payload.get("provenance", {})
    for name in ("code_sha256", "schema_sha256", "semantic_sha256"):
        if SHA256_RE.fullmatch(str(provenance.get(name, ""))) is None:
            raise ValidationError(
                f"{expected_system} component provenance lacks {name}"
            )
    if expected_system == "pp" and not isinstance(
        payload.get("terminal_gate"), dict
    ):
        raise ValidationError("p+p component lacks terminal-gate receipt")


def load_pp_component_certificate(path: Path) -> dict[str, Any]:
    resolved = path.resolve()
    if not resolved.is_file():
        raise ValidationError(f"p+p component certificate is missing: {resolved}")
    try:
        payload = json.loads(resolved.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValidationError(
            f"p+p component certificate is unreadable: {resolved}"
        ) from exc
    verify_component_shape(
        payload,
        expected_system="pp",
        expected_receipt_systems=["pp"],
        expected_row=PP_ROWS[0],
    )
    verify_component_receipt_binding(payload, PP_ROWS)
    verify_pp_terminal_binding(payload)
    terminal_gate = payload["terminal_gate"]
    regenerated = certify_single_receipt(
        Path(str(payload["output_root"])),
        Path(str(payload["preflight_receipt"])),
        selected_rows=PP_ROWS,
        ownership_rows=PP_ROWS,
        schema=COMPONENT_CERTIFICATE_SCHEMA,
        component_system="pp",
        terminal_gate_receipt=Path(str(terminal_gate["path"])),
    )
    if payload != regenerated:
        raise ValidationError(
            "p+p component certificate differs from fresh full revalidation"
        )
    return {
        "certificate_path": str(resolved),
        "certificate_sha256": sha256_file(resolved),
        "payload": payload,
    }


def bind_aggregate_certificate(
    pp_binding: dict[str, Any],
    auau_component: dict[str, Any],
) -> dict[str, Any]:
    pp_component = pp_binding["payload"]
    verify_component_shape(
        pp_component,
        expected_system="pp",
        expected_receipt_systems=["pp"],
        expected_row=PP_ROWS[0],
    )
    verify_component_receipt_binding(pp_component, PP_ROWS)
    verify_pp_terminal_binding(pp_component)
    pp_regenerated = certify_single_receipt(
        Path(str(pp_component["output_root"])),
        Path(str(pp_component["preflight_receipt"])),
        selected_rows=PP_ROWS,
        ownership_rows=PP_ROWS,
        schema=COMPONENT_CERTIFICATE_SCHEMA,
        component_system="pp",
        terminal_gate_receipt=Path(
            str(pp_component["terminal_gate"]["path"])
        ),
    )
    if pp_component != pp_regenerated:
        raise ValidationError(
            "p+p component differs from aggregate-time full revalidation"
        )
    verify_component_shape(
        auau_component,
        expected_system="auau",
        expected_receipt_systems=["auau", "pp"],
        expected_row=AUAU_ROWS[0],
    )
    verify_component_receipt_binding(auau_component, ROWS)
    verify_component_artifact_bindings(auau_component)
    components = {"pp": pp_component, "auau": auau_component}
    all_rows = [
        row
        for system in ("pp", "auau")
        for row in components[system].get("rows", [])
    ]
    observed_systems = [row.get("system") for row in all_rows]
    if observed_systems != ["pp", "auau"]:
        raise ValidationError(
            "aggregate component rows must be exactly one p+p then one Au+Au"
        )
    failures = [
        f"{system}:{failure}"
        for system, component in components.items()
        for failure in component.get("failures", [])
    ]
    for row in all_rows:
        failures.extend(
            f"{row['system']}:{failure}" for failure in row.get("failures", [])
        )
    bindings = {
        "pp": {
            "component_certificate": pp_binding["certificate_path"],
            "component_certificate_sha256": pp_binding["certificate_sha256"],
            "output_root": pp_component["output_root"],
            "preflight_receipt": pp_component["preflight_receipt"],
            "preflight_receipt_sha256": pp_component[
                "preflight_receipt_sha256"
            ],
            "provenance": pp_component["provenance"],
        },
        "auau": {
            "component_certificate": None,
            "component_certificate_sha256": None,
            "output_root": auau_component["output_root"],
            "preflight_receipt": auau_component["preflight_receipt"],
            "preflight_receipt_sha256": auau_component[
                "preflight_receipt_sha256"
            ],
            "provenance": auau_component["provenance"],
        },
    }
    if (
        bindings["pp"]["provenance"]["code_sha256"]
        == bindings["auau"]["provenance"]["code_sha256"]
    ):
        raise ValidationError(
            "replacement aggregate requires corrected p+p and preserved Au+Au "
            "to retain distinct code provenance"
        )
    for identity in ("schema_sha256", "semantic_sha256"):
        if (
            bindings["pp"]["provenance"][identity]
            != bindings["auau"]["provenance"][identity]
        ):
            raise ValidationError(
                f"replacement aggregate requires identical {identity} provenance"
            )
    return {
        "schema": AGGREGATE_CERTIFICATE_SCHEMA,
        "status": "PASS" if not failures else "FAIL",
        "artifact_profile": PROFILE,
        "binding_mode": "INDEPENDENT_PP_COMPONENT_PLUS_PRESERVED_AUAU_RECEIPT",
        "component_bindings": bindings,
        "independent_code_provenance": True,
        "analysis_minimum_bytes": ANALYSIS_MIN_BYTES,
        "sidecar_minimum_bytes_mode": "ARTIFACT_SPECIFIC_GROSS_TRUNCATION_ONLY",
        "cache_replay_applicability": "NOT_APPLICABLE",
        "full_training_authority": 0,
        "full_extraction_authority": False,
        "broad_production_authority": False,
        "rows": all_rows,
        "failures": failures,
    }


def validate_aggregate(args: argparse.Namespace) -> dict[str, Any]:
    pp_binding = load_pp_component_certificate(
        args.pp_component_certificate
    )
    auau_component = certify_single_receipt(
        args.auau_output_root,
        args.auau_preflight_receipt,
        selected_rows=AUAU_ROWS,
        ownership_rows=ROWS,
        schema=COMPONENT_CERTIFICATE_SCHEMA,
        component_system="auau",
    )
    return bind_aggregate_certificate(pp_binding, auau_component)


def validate(args: argparse.Namespace) -> dict[str, Any]:
    aggregate_values = (
        getattr(args, "pp_component_certificate", None),
        getattr(args, "auau_output_root", None),
        getattr(args, "auau_preflight_receipt", None),
    )
    if any(aggregate_values):
        if not all(aggregate_values):
            raise ValidationError(
                "aggregate mode requires p+p certificate plus Au+Au root and receipt"
            )
        if (
            getattr(args, "output_root", None) is not None
            or getattr(args, "preflight_receipt", None) is not None
            or getattr(args, "component_system", None) is not None
            or getattr(args, "terminal_gate_receipt", None) is not None
        ):
            raise ValidationError(
                "aggregate mode cannot mix single-receipt arguments"
            )
        return validate_aggregate(args)
    if (
        getattr(args, "output_root", None) is None
        or getattr(args, "preflight_receipt", None) is None
    ):
        raise ValidationError(
            "single-receipt mode requires output root and preflight receipt"
        )
    component_system = getattr(args, "component_system", None)
    if component_system == "pp":
        return certify_single_receipt(
            args.output_root,
            args.preflight_receipt,
            selected_rows=PP_ROWS,
            ownership_rows=PP_ROWS,
            schema=COMPONENT_CERTIFICATE_SCHEMA,
            component_system="pp",
            terminal_gate_receipt=getattr(
                args, "terminal_gate_receipt", None
            ),
        )
    if component_system is not None:
        raise ValidationError(f"unsupported component system: {component_system}")
    if getattr(args, "terminal_gate_receipt", None) is not None:
        raise ValidationError(
            "terminal-gate receipt requires p+p component mode"
        )
    return certify_single_receipt(
        args.output_root,
        args.preflight_receipt,
        selected_rows=ROWS,
        ownership_rows=ROWS,
        schema=CERTIFICATE_SCHEMA,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--preflight-receipt", type=Path)
    parser.add_argument("--component-system", choices=("pp",))
    parser.add_argument("--terminal-gate-receipt", type=Path)
    parser.add_argument("--pp-component-certificate", type=Path)
    parser.add_argument("--auau-output-root", type=Path)
    parser.add_argument("--auau-preflight-receipt", type=Path)
    parser.add_argument("--output-json", type=Path, required=True)
    return parser.parse_args()


def write_certificate(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def main() -> int:
    args = parse_args()
    if any(
        (
            args.pp_component_certificate,
            args.auau_output_root,
            args.auau_preflight_receipt,
        )
    ):
        failure_schema = AGGREGATE_CERTIFICATE_SCHEMA
    elif args.component_system is not None:
        failure_schema = COMPONENT_CERTIFICATE_SCHEMA
    else:
        failure_schema = CERTIFICATE_SCHEMA
    try:
        payload = validate(args)
    except Exception as exc:
        payload = {
            "schema": failure_schema,
            "status": "FAIL",
            "artifact_profile": PROFILE,
            "failures": [f"{type(exc).__name__}:{exc}"],
            "full_training_authority": 0,
            "full_extraction_authority": False,
            "broad_production_authority": False,
        }
    write_certificate(args.output_json, payload)
    print(json.dumps(payload, sort_keys=True))
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
