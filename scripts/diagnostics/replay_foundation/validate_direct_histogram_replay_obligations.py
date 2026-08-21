#!/usr/bin/env python3
"""Inventory and certify direct-histogram TTree replay obligations.

The direct RecoilJets ROOT output is the detector-execution reference.  This
tool makes every TH1/TH2/TH3/TProfile object an explicit obligation: a
registered normalized-tree recipe must recreate the same key, class, title,
axes, entries, contributors, contents, errors, variances, Sumw2, and flow
bins.  Missing recipes, exemptions, ``DST_REQUIRED`` states, and partial
inventories fail closed.

The inventory command never invents replay recipes.  It records the exact
direct ROOT surface.  A separately authored registry must account for every
entry before ``validate`` can pass.  ``compare`` then checks a TTree-only
offline output against the direct reference without granting the replay
process access to direct histogram contents.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import uproot


INVENTORY_SCHEMA = "DirectHistogramInventoryV1"
REGISTRY_SCHEMA = "DirectHistogramReplayObligationV1"
CERTIFICATE_SCHEMA = "DirectHistogramReplayCertificateV1"
SCHEMA_VERSION = 1
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
FORBIDDEN_STATES = frozenset({"DST_REQUIRED", "EXEMPT", "EXEMPTED", "WAIVED"})
REQUIRED_COMPARE_FIELDS = frozenset(
    {
        "root_key",
        "object_class",
        "title",
        "axes",
        "entries",
        "contributors",
        "contents",
        "errors",
        "variances",
        "sumw2",
        "underflow_overflow",
    }
)


def canonical_json(value: Any) -> bytes:
    return (
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
        .encode("utf-8")
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_json(value: Any) -> str:
    return hashlib.sha256(canonical_json(value)).hexdigest()


def is_sha256(value: Any) -> bool:
    return isinstance(value, str) and SHA256_RE.fullmatch(value) is not None


def strip_cycle(key: str) -> str:
    return key.rsplit(";", 1)[0]


def is_declared_histogram_class(class_name: str) -> bool:
    return class_name.startswith(("TH1", "TH2", "TH3", "TProfile"))


def _member_text(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def _axis_record(axis: Any) -> dict[str, Any]:
    def member(name: str) -> str:
        try:
            return _member_text(axis.member(name))
        except (AttributeError, KeyError, uproot.exceptions.KeyInFileError):
            return ""

    return {
        "name": member("fName"),
        "title": member("fTitle"),
        "edges": [float(value) for value in axis.edges(flow=False)],
    }


def object_inventory_record(root_key: str, class_name: str, obj: Any) -> dict[str, Any]:
    path = Path(root_key)
    directory = "" if str(path.parent) == "." else str(path.parent)
    return {
        "root_key": root_key,
        "object_class": class_name,
        "directory": directory,
        "name": path.name,
        "title": _member_text(getattr(obj, "title", "")),
        "axes": [_axis_record(axis) for axis in obj.axes],
    }


def build_inventory(
    direct_path: Path,
    *,
    artifact_id: str,
    system: str,
    lane: str,
) -> dict[str, Any]:
    if system not in {"pp", "auau"}:
        raise ValueError(f"unsupported_system={system!r}")
    with uproot.open(direct_path) as root_file:
        classnames = root_file.classnames(recursive=True)
        entries = []
        for cycled_key, class_name in sorted(classnames.items()):
            root_key = strip_cycle(cycled_key)
            if not is_declared_histogram_class(class_name):
                continue
            entries.append(
                object_inventory_record(root_key, class_name, root_file[root_key])
            )
    inventory: dict[str, Any] = {
        "schema": INVENTORY_SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "artifact_id": artifact_id,
        "system": system,
        "lane": lane,
        "root_sha256": sha256_file(direct_path),
        "histogram_count": len(entries),
        "entries": entries,
    }
    inventory["inventory_sha256"] = sha256_json(entries)
    return inventory


def obligation_recipe_identity(obligation: Mapping[str, Any]) -> str:
    payload = {
        key: value
        for key, value in obligation.items()
        if key != "recipe_identity_sha256"
    }
    return sha256_json(payload)


def _find_forbidden(value: Any, path: str = "$") -> list[str]:
    failures: list[str] = []
    if isinstance(value, Mapping):
        for key, nested in value.items():
            failures.extend(_find_forbidden(nested, f"{path}.{key}"))
    elif isinstance(value, list):
        for index, nested in enumerate(value):
            failures.extend(_find_forbidden(nested, f"{path}[{index}]"))
    elif isinstance(value, str):
        normalized = value.upper()
        if any(state in normalized for state in FORBIDDEN_STATES):
            failures.append(f"forbidden_state:{path}:{value}")
    return failures


def _nonempty_text(
    record: Mapping[str, Any], field: str, prefix: str, failures: list[str]
) -> None:
    value = record.get(field)
    if not isinstance(value, str) or not value.strip():
        failures.append(f"{prefix}:invalid_{field}")


def _validate_axis(axis: Any, prefix: str, failures: list[str]) -> None:
    if not isinstance(axis, Mapping):
        failures.append(f"{prefix}:axis_not_object")
        return
    if set(axis) != {"name", "title", "edges"}:
        failures.append(f"{prefix}:axis_fields={sorted(axis)}")
    edges = axis.get("edges")
    if not isinstance(edges, list) or len(edges) < 2:
        failures.append(f"{prefix}:invalid_edges")
        return
    try:
        numeric = [float(value) for value in edges]
    except (TypeError, ValueError):
        failures.append(f"{prefix}:nonnumeric_edges")
        return
    if not all(math.isfinite(value) for value in numeric):
        failures.append(f"{prefix}:nonfinite_edges")
    if any(right <= left for left, right in zip(numeric, numeric[1:])):
        failures.append(f"{prefix}:nonmonotonic_edges")


def validate_inventory(inventory: Mapping[str, Any]) -> list[str]:
    failures: list[str] = []
    if inventory.get("schema") != INVENTORY_SCHEMA:
        failures.append(f"inventory_schema={inventory.get('schema')!r}")
    if inventory.get("schema_version") != SCHEMA_VERSION:
        failures.append(f"inventory_schema_version={inventory.get('schema_version')!r}")
    for field in ("artifact_id", "lane"):
        _nonempty_text(inventory, field, "inventory", failures)
    if inventory.get("system") not in {"pp", "auau"}:
        failures.append(f"inventory_system={inventory.get('system')!r}")
    if not is_sha256(inventory.get("root_sha256")):
        failures.append("inventory_invalid_root_sha256")
    entries = inventory.get("entries")
    if not isinstance(entries, list):
        failures.append("inventory_entries_not_array")
        return failures
    if inventory.get("histogram_count") != len(entries):
        failures.append(
            "inventory_count_mismatch:"
            f"{inventory.get('histogram_count')!r}!={len(entries)}"
        )
    if inventory.get("inventory_sha256") != sha256_json(entries):
        failures.append("inventory_sha256_mismatch")
    keys: list[str] = []
    for index, entry in enumerate(entries):
        prefix = f"inventory[{index}]"
        if not isinstance(entry, Mapping):
            failures.append(f"{prefix}:not_object")
            continue
        required = {
            "root_key",
            "object_class",
            "directory",
            "name",
            "title",
            "axes",
        }
        if set(entry) != required:
            failures.append(f"{prefix}:fields={sorted(entry)}")
        for field in ("root_key", "object_class", "name"):
            _nonempty_text(entry, field, prefix, failures)
        if not is_declared_histogram_class(str(entry.get("object_class", ""))):
            failures.append(f"{prefix}:unsupported_class={entry.get('object_class')!r}")
        root_key = str(entry.get("root_key", ""))
        keys.append(root_key)
        path = Path(root_key)
        expected_directory = "" if str(path.parent) == "." else str(path.parent)
        if entry.get("directory") != expected_directory:
            failures.append(f"{prefix}:directory_mismatch")
        if entry.get("name") != path.name:
            failures.append(f"{prefix}:name_mismatch")
        axes = entry.get("axes")
        if not isinstance(axes, list) or not axes:
            failures.append(f"{prefix}:invalid_axes")
        else:
            for axis_index, axis in enumerate(axes):
                _validate_axis(axis, f"{prefix}.axes[{axis_index}]", failures)
    duplicates = sorted({key for key in keys if keys.count(key) > 1})
    if duplicates:
        failures.append(f"inventory_duplicate_keys={duplicates}")
    return failures


def validate_registry(
    registry: Mapping[str, Any], inventory: Mapping[str, Any]
) -> dict[str, Any]:
    failures = validate_inventory(inventory)
    failures.extend(_find_forbidden(registry))
    if registry.get("schema") != REGISTRY_SCHEMA:
        failures.append(f"registry_schema={registry.get('schema')!r}")
    if registry.get("schema_version") != SCHEMA_VERSION:
        failures.append(f"registry_schema_version={registry.get('schema_version')!r}")

    reference = registry.get("direct_reference")
    if not isinstance(reference, Mapping):
        failures.append("direct_reference_not_object")
        reference = {}
    required_reference = {
        "artifact_id",
        "system",
        "lane",
        "root_sha256",
        "inventory_sha256",
        "bundle_sha256",
        "code_sha256",
        "schema_sha256",
        "config_sha256",
        "model_registry_sha256",
        "source_manifest_sha256",
    }
    if set(reference) != required_reference:
        failures.append(f"direct_reference_fields={sorted(reference)}")
    for field in (
        "root_sha256",
        "inventory_sha256",
        "bundle_sha256",
        "code_sha256",
        "schema_sha256",
        "config_sha256",
        "model_registry_sha256",
        "source_manifest_sha256",
    ):
        if not is_sha256(reference.get(field)):
            failures.append(f"direct_reference_invalid_{field}")
    for field in ("artifact_id", "system", "lane", "root_sha256", "inventory_sha256"):
        if reference.get(field) != inventory.get(field):
            failures.append(
                f"direct_reference_mismatch:{field}:"
                f"{reference.get(field)!r}!={inventory.get(field)!r}"
            )

    obligations = registry.get("obligations")
    if not isinstance(obligations, list) or not obligations:
        failures.append("obligations_missing_or_empty")
        obligations = []
    obligation_by_key: dict[str, Mapping[str, Any]] = {}
    for index, obligation in enumerate(obligations):
        prefix = f"obligation[{index}]"
        if not isinstance(obligation, Mapping):
            failures.append(f"{prefix}:not_object")
            continue
        root_key = obligation.get("root_key")
        if not isinstance(root_key, str) or not root_key:
            failures.append(f"{prefix}:invalid_root_key")
            continue
        if root_key in obligation_by_key:
            failures.append(f"duplicate_obligation:{root_key}")
        obligation_by_key[root_key] = obligation
        required_fields = {
            "root_key",
            "object_class",
            "directory",
            "name",
            "title",
            "axes",
            "source_inputs",
            "selection_rule",
            "fill_rule",
            "ordering_policy",
            "tie_policy",
            "weight_contract",
            "flow_policy",
            "response_classification",
            "equality_contract",
            "analysis_consumers",
            "ian_consumers",
            "recipe_identity_sha256",
        }
        if set(obligation) != required_fields:
            failures.append(f"{prefix}:fields={sorted(obligation)}")
        for field in (
            "selection_rule",
            "fill_rule",
            "ordering_policy",
            "tie_policy",
            "response_classification",
        ):
            _nonempty_text(obligation, field, prefix, failures)
        inputs = obligation.get("source_inputs")
        if not isinstance(inputs, list) or not inputs:
            failures.append(f"{prefix}:missing_source_inputs")
        else:
            for input_index, source_input in enumerate(inputs):
                source_prefix = f"{prefix}.source_inputs[{input_index}]"
                if not isinstance(source_input, Mapping):
                    failures.append(f"{source_prefix}:not_object")
                    continue
                if set(source_input) != {"tree", "branches"}:
                    failures.append(f"{source_prefix}:fields={sorted(source_input)}")
                _nonempty_text(source_input, "tree", source_prefix, failures)
                branches = source_input.get("branches")
                if (
                    not isinstance(branches, list)
                    or not branches
                    or any(not isinstance(value, str) or not value for value in branches)
                    or len(set(branches)) != len(branches)
                ):
                    failures.append(f"{source_prefix}:invalid_branches")
        weight = obligation.get("weight_contract")
        if not isinstance(weight, Mapping):
            failures.append(f"{prefix}:weight_contract_not_object")
        else:
            if set(weight) != {"components", "application_count", "rule"}:
                failures.append(f"{prefix}:weight_contract_fields={sorted(weight)}")
            if weight.get("application_count") != 1:
                failures.append(f"{prefix}:weight_application_count_not_one")
            components = weight.get("components")
            if (
                not isinstance(components, list)
                or not components
                or any(not isinstance(value, str) or not value for value in components)
                or len(set(components)) != len(components)
            ):
                failures.append(f"{prefix}:invalid_weight_components")
            _nonempty_text(weight, "rule", f"{prefix}.weight_contract", failures)
        flow = obligation.get("flow_policy")
        if not isinstance(flow, Mapping) or set(flow) != {"underflow", "overflow"}:
            failures.append(f"{prefix}:invalid_flow_policy")
        else:
            _nonempty_text(flow, "underflow", f"{prefix}.flow_policy", failures)
            _nonempty_text(flow, "overflow", f"{prefix}.flow_policy", failures)
        equality = obligation.get("equality_contract")
        if not isinstance(equality, Mapping):
            failures.append(f"{prefix}:equality_contract_not_object")
        else:
            if set(equality) != {"status", "numerical_mode", "compared_fields"}:
                failures.append(f"{prefix}:equality_contract_fields={sorted(equality)}")
            if equality.get("status") != "REPLAY_EXACT":
                failures.append(f"{prefix}:status_not_replay_exact")
            if equality.get("numerical_mode") != "EXACT":
                failures.append(f"{prefix}:numerical_mode_not_exact")
            compared = equality.get("compared_fields")
            if not isinstance(compared, list) or set(compared) != REQUIRED_COMPARE_FIELDS:
                failures.append(f"{prefix}:incomplete_compared_fields")
        for field in ("analysis_consumers", "ian_consumers"):
            consumers = obligation.get(field)
            if not isinstance(consumers, list) or any(
                not isinstance(value, str) or not value for value in consumers
            ):
                failures.append(f"{prefix}:invalid_{field}")
        expected_identity = obligation_recipe_identity(obligation)
        if obligation.get("recipe_identity_sha256") != expected_identity:
            failures.append(f"{prefix}:recipe_identity_sha256_mismatch")

    inventory_by_key = {
        entry["root_key"]: entry
        for entry in inventory.get("entries", [])
        if isinstance(entry, Mapping) and "root_key" in entry
    }
    missing = sorted(set(inventory_by_key) - set(obligation_by_key))
    extra = sorted(set(obligation_by_key) - set(inventory_by_key))
    if missing:
        failures.append(f"unregistered_direct_histograms={missing}")
    if extra:
        failures.append(f"registry_keys_not_in_direct_inventory={extra}")
    for root_key in sorted(set(inventory_by_key) & set(obligation_by_key)):
        inventory_surface = inventory_by_key[root_key]
        obligation_surface = {
            field: obligation_by_key[root_key].get(field)
            for field in (
                "root_key",
                "object_class",
                "directory",
                "name",
                "title",
                "axes",
            )
        }
        if obligation_surface != inventory_surface:
            failures.append(f"direct_surface_mismatch:{root_key}")

    return {
        "schema": CERTIFICATE_SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "status": "PASS" if not failures else "FAIL",
        "direct_artifact_id": inventory.get("artifact_id"),
        "direct_root_sha256": inventory.get("root_sha256"),
        "inventory_sha256": inventory.get("inventory_sha256"),
        "direct_histogram_count": len(inventory_by_key),
        "registered_obligation_count": len(obligation_by_key),
        "coverage_fraction": (
            1.0
            if inventory_by_key and set(inventory_by_key) == set(obligation_by_key)
            else (
                len(set(inventory_by_key) & set(obligation_by_key))
                / len(inventory_by_key)
                if inventory_by_key
                else 0.0
            )
        ),
        "unregistered_count": len(missing),
        "extra_registry_count": len(extra),
        "forbidden_state_count": sum(
            failure.startswith("forbidden_state:") for failure in failures
        ),
        "registry_sha256": sha256_json(registry),
        "failures": failures,
    }


def _safe_member(obj: Any, name: str) -> Any:
    try:
        return obj.member(name)
    except (AttributeError, KeyError, uproot.exceptions.KeyInFileError):
        return None


def _array_payload(value: Any) -> dict[str, Any] | None:
    if value is None:
        return None
    if hasattr(value, "member"):
        nested = _safe_member(value, "fArray")
        if nested is not None:
            value = nested
    try:
        array = np.asarray(value)
    except (TypeError, ValueError):
        return {"text": _member_text(value)}
    if array.dtype.kind == "O":
        return {"text": _member_text(value)}
    return {
        "shape": list(array.shape),
        "dtype": str(array.dtype),
        "bytes_sha256": hashlib.sha256(array.tobytes(order="C")).hexdigest(),
    }


def histogram_semantic_record(root_key: str, class_name: str, obj: Any) -> dict[str, Any]:
    surface = object_inventory_record(root_key, class_name, obj)
    values = np.asarray(obj.values(flow=True))
    try:
        variances = obj.variances(flow=True)
    except (AttributeError, TypeError, ValueError):
        variances = None
    variances_array = None if variances is None else np.asarray(variances)
    try:
        errors = obj.errors(flow=True)
    except (AttributeError, TypeError, ValueError):
        errors = None
    errors_array = None if errors is None else np.asarray(errors)
    return {
        **surface,
        "entries": _member_text(_safe_member(obj, "fEntries")),
        "contributors": {
            "bin_entries": _array_payload(_safe_member(obj, "fBinEntries")),
            "bin_sumw2": _array_payload(_safe_member(obj, "fBinSumw2")),
        },
        "contents": _array_payload(values),
        "errors": _array_payload(errors_array),
        "variances": _array_payload(variances_array),
        "sumw2": _array_payload(_safe_member(obj, "fSumw2")),
        "underflow_overflow": "included_in_flow_arrays",
    }


def compare_outputs(
    registry: Mapping[str, Any],
    direct_path: Path,
    replay_path: Path,
) -> dict[str, Any]:
    reference = registry.get("direct_reference", {})
    if not isinstance(reference, Mapping):
        reference = {}
    try:
        direct_inventory = build_inventory(
            direct_path,
            artifact_id=str(reference.get("artifact_id", "")),
            system=str(reference.get("system", "")),
            lane=str(reference.get("lane", "")),
        )
    except (OSError, ValueError, KeyError) as error:
        return {
            "schema": CERTIFICATE_SCHEMA,
            "schema_version": SCHEMA_VERSION,
            "status": "FAIL",
            "failures": [f"direct_inventory_error:{type(error).__name__}:{error}"],
        }
    registry_certificate = validate_registry(registry, direct_inventory)
    if registry_certificate["status"] != "PASS":
        return {
            "schema": CERTIFICATE_SCHEMA,
            "schema_version": SCHEMA_VERSION,
            "status": "FAIL",
            "registry_certificate": registry_certificate,
            "failures": ["registry_not_certified_for_direct_reference"],
        }
    obligations = registry.get("obligations")
    failures: list[str] = []
    if not isinstance(obligations, list) or not obligations:
        return {
            "schema": CERTIFICATE_SCHEMA,
            "schema_version": SCHEMA_VERSION,
            "status": "FAIL",
            "failures": ["obligations_missing_or_empty"],
        }
    expected_keys = {str(item.get("root_key")) for item in obligations}
    comparisons: list[dict[str, Any]] = []
    with uproot.open(direct_path) as direct_file, uproot.open(replay_path) as replay_file:
        direct_classes = {
            strip_cycle(key): value
            for key, value in direct_file.classnames(recursive=True).items()
            if is_declared_histogram_class(value)
        }
        replay_classes = {
            strip_cycle(key): value
            for key, value in replay_file.classnames(recursive=True).items()
            if is_declared_histogram_class(value)
        }
        if set(direct_classes) != expected_keys:
            failures.append(
                "direct_key_inventory_mismatch:"
                f"missing={sorted(expected_keys - set(direct_classes))}:"
                f"extra={sorted(set(direct_classes) - expected_keys)}"
            )
        if set(replay_classes) != expected_keys:
            failures.append(
                "replay_key_inventory_mismatch:"
                f"missing={sorted(expected_keys - set(replay_classes))}:"
                f"extra={sorted(set(replay_classes) - expected_keys)}"
            )
        for root_key in sorted(expected_keys & set(direct_classes) & set(replay_classes)):
            direct_record = histogram_semantic_record(
                root_key, direct_classes[root_key], direct_file[root_key]
            )
            replay_record = histogram_semantic_record(
                root_key, replay_classes[root_key], replay_file[root_key]
            )
            equal = direct_record == replay_record
            comparisons.append(
                {
                    "root_key": root_key,
                    "status": "PASS" if equal else "FAIL",
                    "direct_semantic_sha256": sha256_json(direct_record),
                    "replay_semantic_sha256": sha256_json(replay_record),
                }
            )
            if not equal:
                failures.append(f"semantic_mismatch:{root_key}")
    return {
        "schema": CERTIFICATE_SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "status": "PASS" if not failures else "FAIL",
        "direct_root_sha256": sha256_file(direct_path),
        "replay_root_sha256": sha256_file(replay_path),
        "registry_sha256": sha256_json(registry),
        "compared_histogram_count": len(comparisons),
        "comparisons": comparisons,
        "failures": failures,
    }


def load_json(path: Path) -> Any:
    with path.open(encoding="utf-8") as stream:
        return json.load(stream)


def write_report(report: Mapping[str, Any], output: Path | None) -> None:
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if output is None:
        print(text, end="")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(text, encoding="utf-8")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    inventory_parser = subparsers.add_parser(
        "inventory", help="Record the complete direct TH*/TProfile inventory"
    )
    inventory_parser.add_argument("--direct", required=True, type=Path)
    inventory_parser.add_argument("--artifact-id", required=True)
    inventory_parser.add_argument("--system", required=True, choices=("pp", "auau"))
    inventory_parser.add_argument("--lane", required=True)
    inventory_parser.add_argument("--output", required=True, type=Path)

    validate_parser = subparsers.add_parser(
        "validate", help="Require one exact registered replay recipe per direct object"
    )
    validate_parser.add_argument("--inventory", required=True, type=Path)
    validate_parser.add_argument("--registry", required=True, type=Path)
    validate_parser.add_argument("--output", type=Path)

    compare_parser = subparsers.add_parser(
        "compare", help="Compare direct and independent TTree-only replay ROOT outputs"
    )
    compare_parser.add_argument("--registry", required=True, type=Path)
    compare_parser.add_argument("--direct", required=True, type=Path)
    compare_parser.add_argument("--replay", required=True, type=Path)
    compare_parser.add_argument("--output", type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.command == "inventory":
        report = build_inventory(
            args.direct,
            artifact_id=args.artifact_id,
            system=args.system,
            lane=args.lane,
        )
        write_report(report, args.output)
        return 0
    if args.command == "validate":
        report = validate_registry(load_json(args.registry), load_json(args.inventory))
        write_report(report, args.output)
        return 0 if report["status"] == "PASS" else 1
    if args.command == "compare":
        report = compare_outputs(
            load_json(args.registry),
            args.direct,
            args.replay,
        )
        write_report(report, args.output)
        return 0 if report["status"] == "PASS" else 1
    raise AssertionError(f"unhandled_command={args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
