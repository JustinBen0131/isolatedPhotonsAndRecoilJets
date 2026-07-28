#!/usr/bin/env python3
"""Certify the bounded THE-134 direct/sidecar-only differential canary.

The canary owns exactly p+p Jet8 and Au+Au embedded Jet12, with one direct and
one writer arm per system.  Analysis ROOT histograms must remain exactly
neutral.  Writer analysis files must contain the frozen validation-only
markers and no serialized ReplayFoundation directory.  The compact
RJPhotonTrainingViewV1 sidecars use their artifact-specific semantic health
contract rather than the analysis ROOT 50 kB size rule.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import re
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


def parse_preflight(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise ValidationError(f"preflight receipt is missing: {path}")
    fields: dict[str, str] = {}
    pinned: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        if "=" in raw and not raw.startswith((" ", "\t")):
            key, value = raw.split("=", 1)
            fields[key] = value
            continue
        parts = raw.split()
        if len(parts) == 2 and SHA256_RE.fullmatch(parts[0]):
            pinned[Path(parts[1]).name] = parts[0]
    required = ("tag", "code_sha256", "schema_sha", "semantic_sha", "only_keys")
    missing = [name for name in required if not fields.get(name)]
    if missing:
        raise ValidationError(f"preflight receipt lacks fields: {missing}")
    for name in ("code_sha256", "schema_sha", "semantic_sha"):
        if SHA256_RE.fullmatch(fields[name]) is None:
            raise ValidationError(f"preflight {name} is not SHA-256")
    expected_only_keys = {
        f"{arm}:{row['lane']}:{row['sample']}"
        for row in ROWS
        for arm in ("direct", "writer")
    }
    if set(fields["only_keys"].split(",")) != expected_only_keys:
        raise ValidationError("preflight exact four-row ownership differs")
    if fields.get("writer_extra_common_template", "").find(
        "RJ_THE134_MULTIVIEW_SIDECAR_ONLY_V1=1"
    ) < 0:
        raise ValidationError("preflight lacks the sidecar-only writer profile")
    if fields.get("extra_pp_template", "").find(
        "RJ_PP_PHOTONID_EXTRACT_ONLY=1"
    ) < 0:
        raise ValidationError("preflight lacks the shared p+p extraction profile")
    if fields.get("extra_auau_template", "").find(
        "RJ_AUAU_BDT_EXTRACT_ONLY=1"
    ) < 0:
        raise ValidationError("preflight lacks the shared Au+Au extraction profile")
    if fields.get("writer_extra_pp_template") or fields.get(
        "writer_extra_auau_template"
    ):
        raise ValidationError("preflight contains unexpected system-specific writer extensions")
    return {"fields": fields, "pinned": pinned, "sha256": sha256_file(path)}


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
    expected["rj_replay_source_sha256"] = source_sha256(fields["tag"], row)
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
    source_hash = source_sha256(preflight["fields"]["tag"], row)
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


def validate(args: argparse.Namespace) -> dict[str, Any]:
    output_root = args.output_root.resolve()
    preflight = parse_preflight(args.preflight_receipt.resolve())
    preparer = load_module("the134_sidecar_differential_preparer", PREPARER_PATH)
    row_reports: list[dict[str, Any]] = []
    for row in ROWS:
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
    failures = [
        f"{row['system']}:{failure}"
        for row in row_reports
        for failure in row["failures"]
    ]
    return {
        "schema": CERTIFICATE_SCHEMA,
        "status": "PASS" if not failures else "FAIL",
        "artifact_profile": PROFILE,
        "output_root": str(output_root),
        "preflight_receipt": str(args.preflight_receipt.resolve()),
        "preflight_receipt_sha256": preflight["sha256"],
        "analysis_minimum_bytes": ANALYSIS_MIN_BYTES,
        "sidecar_minimum_bytes_mode": "ARTIFACT_SPECIFIC_GROSS_TRUNCATION_ONLY",
        "cache_replay_applicability": "NOT_APPLICABLE",
        "full_training_authority": 0,
        "full_extraction_authority": False,
        "broad_production_authority": False,
        "rows": row_reports,
        "failures": failures,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--preflight-receipt", type=Path, required=True)
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
    try:
        payload = validate(args)
    except Exception as exc:
        payload = {
            "schema": CERTIFICATE_SCHEMA,
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
