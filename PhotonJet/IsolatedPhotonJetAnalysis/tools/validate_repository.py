#!/usr/bin/env python3
"""Fail-closed public-boundary and reproducibility checks for this repository."""

from __future__ import annotations

import argparse
import ast
import gzip
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import re
import stat
import subprocess
import sys
import tarfile
import tempfile
from typing import Any
import unicodedata
import zipfile


REPO = Path(__file__).resolve().parents[1]
TEXT_SUFFIXES = {
    "",
    ".C",
    ".cc",
    ".cfg",
    ".cmake",
    ".cpp",
    ".h",
    ".hpp",
    ".html",
    ".ini",
    ".ipynb",
    ".json",
    ".log",
    ".md",
    ".py",
    ".rst",
    ".sh",
    ".tex",
    ".toml",
    ".tsv",
    ".txt",
    ".yaml",
    ".yml",
}
SENSITIVE_PATTERNS = {
    "task-style identifier": re.compile(
        r"(?<![A-Za-z0-9])(?:"
        r"(?!(?:SHA|INT|GPL|UTF|ISO|RFC|BDT|CHI)[-_/ ])[A-Z]{3}[-_/ ][0-9]{1,7}(?![A-Za-z0-9])"
        r"|(?!(?:SHA|INT|GPL|UTF|ISO|RFC|BDT|CHI)[0-9])[A-Z]{3}[0-9]{1,7}(?![A-Za-z0-9])"
        r"|(?i:THE)(?:[-_/ ]?)[0-9]{1,7}(?![A-Za-z0-9])"
        r")",
    ),
    "internal versioned adapter": re.compile(
        r"\b[A-Za-z0-9_]*(?:Foundation|Replay|Adapter)V[0-9]+\b"
    ),
    "workspace control path": re.compile(
        r"(?:^|[\s'\"`/])[A-Za-z0-9_.-]+_(?:context|notes)/"
        r"|(?:^|[\s'\"`/])scripts/(?:[a-z]{2}|control|governance|operations)(?:/|\b)",
        re.IGNORECASE,
    ),
    "development-only tree state": re.compile(
        r"\b[A-Z]+_TREE(?:__[A-Z0-9_]+)+\b"
    ),
    "absolute account path": re.compile(
        "(?:"
        + r"/(?:" + "|".join(("Users", "home")) + r")/[^/\s\x00]+/"
        + "|/" + "sphenix" + r"/[^/\s\x00]+/[^\s\x00]+"
        + ")",
        re.IGNORECASE,
    ),
    "ephemeral path": re.compile(
        "/" + r"(?:private/)?tmp/[^\s\x00]+|/" + "var" + r"/folders/[^\s\x00]+",
        re.IGNORECASE,
    ),
    "site filesystem path": re.compile(
        r"/(?:lustre|scratch|work|gpfs[0-9]*)/[^\s\x00]+", re.IGNORECASE
    ),
    "thread or UUID identifier": re.compile(
        r"\b[0-9a-f]{8}-[0-9a-f]{4}-[1-8][0-9a-f]{3}-[89ab0-9][0-9a-f]{3}-[0-9a-f]{12}\b",
        re.IGNORECASE,
    ),
    "automated authorship marker": re.compile(
        "|".join(("co" + "-authored-by", "assisted" + "-by", "generated" + "-by", "language" + " model")),
        re.IGNORECASE,
    ),
    "workflow governance vocabulary": re.compile(
        "(?:"
        + "agent" + r"[-_ ]?(?:control|context|task|session|workstream)"
        + "|" + "control" + r"[-_ ]?plane"
        + ")",
        re.IGNORECASE,
    ),
    "assistant or vendor provenance": re.compile(
        "|".join((
            "co" + "dex", "cl" + "aude", "chat" + "gpt", "open" + "ai", "anth" + "ropic",
            "g" + r"pt[-_ ]?[0-9]+", "a" + r"i[-_ ]+assistant",
            "large" + r"\s+language\s+model",
            "generated" + r"\s+(?:with|by)\s+(?:an?\s+)?(?:assistant|gpt|model)",
        )),
        re.IGNORECASE,
    ),
    "private credential material": re.compile(
        "(?:"
        + "-----BEGIN " + r"(?:RSA |EC |OPENSSH )?PRIVATE KEY-----"
        + "|" + "github" + r"_pat_[A-Za-z0-9_]{20,}"
        + "|" + "gh" + r"[pousr]_[A-Za-z0-9]{20,}"
        + "|" + "AK" + r"IA[0-9A-Z]{16}"
        + "|" + "AI" + r"za[0-9A-Za-z_-]{30,}"
        + "|" + "xox" + r"[abprs]-[0-9A-Za-z-]{10,}"
        + "|(?:" + "pass" + r"word|" + "to" + r"ken|" + "se" + r"cret)\s*[:=]\s*['\"]?[^\s'\";,]{8,}"
        + "|[A-Za-z][A-Za-z0-9+.-]*://[^\s/@:]+:[^\s/@]+@"
        + ")",
        re.IGNORECASE,
    ),
    "private contact identifier": re.compile(
        r"\b[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}\b",
        re.IGNORECASE,
    ),
    "scheduler or session identifier": re.compile(
        r"\b(?:cluster|job|session)[-_: ]?[0-9]{3,}\b",
        re.IGNORECASE,
    ),
    "private host or account identifier": re.compile(
        r"\b(?:user(?:name)?|account|host(?:name)?)\s*[:=]\s*['\"]?[A-Za-z0-9_.-]{2,}"
        r"|\b(?:login|submit|compute|batch)[-_]?[0-9]{1,5}\.[A-Za-z0-9.-]+\b",
        re.IGNORECASE,
    ),
    "nonportable URI or platform path": re.compile(
        r"(?:\b[A-Za-z]:\\|\\\\[A-Za-z0-9_.-]+\\|\b(?:file|root)://[^\s'\"`]+)",
        re.IGNORECASE,
    ),
}
RAW_PLOT_SEMANTIC = re.compile(
    r"(?:DrawLatex|\.text|\.annotate)\s*\([^;]*?[\"'][^\"']*"
    r"(?:GeV|#eta|p_\{?T|centrality|ABCD|sPHENIX)",
    re.IGNORECASE | re.DOTALL,
)
PLOT_CODE_EXCEPTIONS = {
    "python/photonjet/plotting/contract.py",
    "python/photonjet/plotting/root_bridge.py",
    "macros/PhotonJetPlotStyle.h",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _public_files() -> list[Path]:
    result: list[Path] = []
    for path in REPO.rglob("*"):
        if path.is_symlink():
            result.append(path)
            continue
        if not path.is_file():
            continue
        relative = path.relative_to(REPO)
        if relative.parts[0] in {".git", "build"} or "__pycache__" in relative.parts:
            continue
        result.append(path)
    return sorted(result, key=lambda item: item.relative_to(REPO).as_posix())


def _text(path: Path) -> str | None:
    if path.suffix not in TEXT_SUFFIXES:
        return None
    try:
        return unicodedata.normalize("NFKC", path.read_text(encoding="utf-8"))
    except UnicodeDecodeError:
        return None


def _python_raw_annotation(text: str) -> tuple[int, str] | None:
    try:
        tree = ast.parse(text)
    except SyntaxError:
        return None
    physics = re.compile(r"(?:GeV|#eta|p_\{?T|centrality|ABCD|sPHENIX)", re.IGNORECASE)
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        name = ""
        if isinstance(node.func, ast.Attribute):
            name = node.func.attr
        elif isinstance(node.func, ast.Name):
            name = node.func.id
        if name not in {"text", "annotate", "DrawLatex", "DrawLatexNDC"}:
            continue
        for argument in [*node.args, *(item.value for item in node.keywords)]:
            try:
                value = ast.literal_eval(argument)
            except (ValueError, TypeError):
                continue
            if isinstance(value, str) and physics.search(value):
                return int(getattr(argument, "lineno", node.lineno)), value
    return None


MAX_CONTAINER_BYTES = 64 * 1024 * 1024
MAX_CONTAINER_MEMBERS = 4096


def _decoded_binary_strings(data: bytes) -> list[str]:
    """Extract portable printable strings without treating random bytes as text."""

    result = [
        match.group().decode("ascii")
        for match in re.finditer(rb"[\x20-\x7e]{5,}", data)
    ]
    for pattern, encoding in (
        (rb"(?:[\x20-\x7e]\x00){4,}", "utf-16le"),
        (rb"(?:\x00[\x20-\x7e]){4,}", "utf-16be"),
    ):
        result.extend(
            match.group().decode(encoding)
            for match in re.finditer(pattern, data)
        )
    return result


def _sensitive_labels(strings: list[str]) -> set[str]:
    normalized = [unicodedata.normalize("NFKC", value) for value in strings]
    return {
        label
        for label, pattern in SENSITIVE_PATTERNS.items()
        if any(pattern.search(value) for value in normalized)
    }


def _root_metadata_strings(source: Any) -> tuple[list[str], list[str]]:
    """Decode compressed ROOT object metadata and any explicit string branches."""

    try:
        import uproot
    except ImportError:
        return [], ["ROOT metadata inspection requires uproot"]
    strings: list[str] = []
    failures: list[str] = []
    try:
        with uproot.open(source) as root:
            classnames = root.classnames(recursive=True)
            if len(classnames) > MAX_CONTAINER_MEMBERS:
                return [], ["ROOT object count exceeds the inspection bound"]
            for key, classname in classnames.items():
                strings.extend((str(key), str(classname)))
                try:
                    obj = root[key]
                except Exception as error:
                    failures.append(f"ROOT object metadata cannot be read: {key}: {error}")
                    continue
                for attribute in ("name", "title", "classname"):
                    value = getattr(obj, attribute, None)
                    if isinstance(value, str):
                        strings.append(value)
                if classname == "TObjString":
                    strings.append(str(obj))
                for member_name in getattr(obj, "member_names", ()):
                    try:
                        value = obj.member(member_name)
                    except Exception:
                        continue
                    if isinstance(value, str):
                        strings.append(value)
                    elif isinstance(value, bytes):
                        strings.extend(_decoded_binary_strings(value))
                for branch in getattr(obj, "branches", ()):
                    strings.extend(
                        str(value)
                        for value in (
                            getattr(branch, "name", ""),
                            getattr(branch, "title", ""),
                            getattr(branch, "typename", ""),
                        )
                        if value
                    )
                    branch_type = str(getattr(branch, "typename", "")).lower()
                    if "string" not in branch_type:
                        continue
                    values = branch.array(library="np")
                    if getattr(values, "size", len(values)) > MAX_CONTAINER_MEMBERS:
                        failures.append(f"ROOT string branch exceeds inspection bound: {key}/{branch.name}")
                        continue
                    for value in values:
                        if isinstance(value, str):
                            strings.append(value)
                        elif isinstance(value, bytes):
                            strings.extend(_decoded_binary_strings(value))
    except Exception as error:
        failures.append(f"ROOT metadata inspection failed: {error}")
    return strings, failures


def _bounded_read(stream: Any) -> bytes:
    data = stream.read(MAX_CONTAINER_BYTES + 1)
    if len(data) > MAX_CONTAINER_BYTES:
        raise ValueError("decompressed member exceeds inspection bound")
    return data


def _inspect_embedded_payload(name: str, data: bytes) -> tuple[set[str], list[str]]:
    labels = _sensitive_labels([name, *_decoded_binary_strings(data)])
    failures: list[str] = []
    if name.lower().endswith(".root"):
        root_strings, root_failures = _root_metadata_strings(io.BytesIO(data))
        labels.update(_sensitive_labels(root_strings))
        failures.extend(root_failures)
    return labels, failures


def _archive_findings(path: Path) -> tuple[set[str], list[str]] | None:
    labels: set[str] = set()
    failures: list[str] = []
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as archive:
            members = archive.infolist()
            if len(members) > MAX_CONTAINER_MEMBERS:
                return labels, ["archive member count exceeds inspection bound"]
            for member in members:
                if member.is_dir():
                    continue
                if member.flag_bits & 0x1:
                    failures.append(f"encrypted archive member is not inspectable: {member.filename}")
                    continue
                if member.file_size > MAX_CONTAINER_BYTES:
                    failures.append(f"archive member exceeds inspection bound: {member.filename}")
                    continue
                with archive.open(member) as stream:
                    found, errors = _inspect_embedded_payload(member.filename, _bounded_read(stream))
                labels.update(found)
                failures.extend(errors)
        return labels, failures
    if tarfile.is_tarfile(path):
        with tarfile.open(path, mode="r:*") as archive:
            members = archive.getmembers()
            if len(members) > MAX_CONTAINER_MEMBERS:
                return labels, ["archive member count exceeds inspection bound"]
            for member in members:
                if not member.isfile():
                    if member.issym() or member.islnk():
                        failures.append(f"archive link member is forbidden: {member.name}")
                    continue
                if member.size > MAX_CONTAINER_BYTES:
                    failures.append(f"archive member exceeds inspection bound: {member.name}")
                    continue
                stream = archive.extractfile(member)
                if stream is None:
                    failures.append(f"archive member cannot be read: {member.name}")
                    continue
                with stream:
                    found, errors = _inspect_embedded_payload(member.name, _bounded_read(stream))
                labels.update(found)
                failures.extend(errors)
        return labels, failures
    if path.suffix.lower() == ".gz":
        with gzip.open(path, "rb") as stream:
            found, errors = _inspect_embedded_payload(path.stem, _bounded_read(stream))
        labels.update(found)
        failures.extend(errors)
        return labels, failures
    return None


def _binary_sensitive_findings(path: Path) -> list[str]:
    data = path.read_bytes()
    labels = _sensitive_labels(_decoded_binary_strings(data))
    failures: list[str] = []
    if path.suffix.lower() == ".root":
        root_strings, root_failures = _root_metadata_strings(path)
        labels.update(_sensitive_labels(root_strings))
        failures.extend(root_failures)
    else:
        archive = _archive_findings(path)
        if archive is None:
            failures.append(f"unsupported binary format: {path.suffix or '<none>'}")
        else:
            archive_labels, archive_failures = archive
            labels.update(archive_labels)
            failures.extend(archive_failures)
    return [*sorted(labels), *failures]


def _json_object(path: Path, label: str, findings: list[str]) -> dict[str, Any] | None:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        findings.append(f"{label} cannot be decoded: {error}")
        return None
    if not isinstance(value, dict):
        findings.append(f"{label} must be a JSON object")
        return None
    return value


def _safe_tracked_file(
    raw_path: Any,
    tracked: set[str],
    label: str,
    findings: list[str],
) -> tuple[str, Path] | None:
    relative = str(raw_path)
    candidate = Path(relative)
    if not relative or candidate.is_absolute() or ".." in candidate.parts:
        findings.append(f"unsafe {label} path: {relative!r}")
        return None
    target = REPO / candidate
    try:
        target.resolve(strict=True).relative_to(REPO.resolve())
    except (OSError, ValueError):
        findings.append(f"{label} path leaves the repository: {relative}")
        return None
    if target.is_symlink() or not target.is_file():
        findings.append(f"{label} is missing or not a regular file: {relative}")
        return None
    if relative not in tracked:
        findings.append(f"{label} is not tracked by the release surface: {relative}")
        return None
    return relative, target


def _valid_sha256(value: Any) -> bool:
    return bool(re.fullmatch(r"[0-9a-f]{64}", str(value)))


def _validate_branch_contract(findings: list[str]) -> tuple[dict[str, Any] | None, int]:
    path = REPO / "contracts/photonjet_trees_v1_branches.json"
    contract = _json_object(path, "branch contract", findings)
    if contract is None:
        return None, 0
    if contract.get("schema") != "PhotonJetTreeBranchContractV1":
        findings.append("branch contract schema is invalid")
    if contract.get("contract") != "PhotonJetTrees_v1":
        findings.append("branch contract release name is invalid")
    order = contract.get("tree_order")
    branches = contract.get("branches")
    if not isinstance(order, list) or not isinstance(branches, dict):
        findings.append("branch contract tree order or branches are invalid")
        return contract, 0
    if order != list(branches):
        findings.append("branch contract trees are not unique and ordered")
    total = 0
    for tree_name in order:
        rows = branches.get(tree_name)
        if not isinstance(rows, list):
            findings.append(f"branch rows are invalid: {tree_name}")
            continue
        names = [row.get("name") for row in rows if isinstance(row, dict)]
        if len(names) != len(rows) or len(names) != len(set(names)):
            findings.append(f"branch names are invalid or duplicated: {tree_name}")
        for row in rows:
            if not isinstance(row, dict) or set(row) != {"name", "typename"}:
                findings.append(f"branch row shape is invalid: {tree_name}")
        total += len(rows)
    return contract, total


def _validate_root_fixture(
    target: Path,
    contract: dict[str, Any] | None,
    label: str,
    findings: list[str],
) -> None:
    if target.read_bytes()[:4] != b"root":
        findings.append(f"fixture does not have ROOT magic: {label}")
        return
    if contract is None:
        return
    try:
        import uproot

        with uproot.open(target) as root:
            observed_trees = list(root.keys(cycle=False))
            expected_trees = contract["tree_order"]
            if observed_trees != expected_trees:
                findings.append(f"fixture tree inventory differs: {label}")
                return
            for tree_name in expected_trees:
                expected = [
                    (row["name"], row["typename"])
                    for row in contract["branches"][tree_name]
                ]
                observed = list(root[tree_name].typenames().items())
                if observed != expected:
                    findings.append(f"fixture branch inventory differs: {label}:{tree_name}")
    except Exception as error:
        findings.append(f"fixture ROOT inventory cannot be verified: {label}: {error}")


def _validate_fixture_manifest(
    tracked: set[str],
    contract: dict[str, Any] | None,
    findings: list[str],
) -> None:
    manifest_path = REPO / "tests/fixtures/MANIFEST.json"
    manifest = _json_object(manifest_path, "fixture manifest", findings)
    if manifest is None:
        return
    if manifest.get("schema") != "PhotonJetEngineeringFixtureManifestV1":
        findings.append("fixture manifest schema is invalid")
    rows = manifest.get("fixtures")
    if not isinstance(rows, list):
        findings.append("fixture manifest rows must be a list")
        return
    row_paths = [str(row.get("path", "")) for row in rows if isinstance(row, dict)]
    if len(row_paths) != len(rows) or row_paths != sorted(row_paths) or len(row_paths) != len(set(row_paths)):
        findings.append("fixture manifest paths are not unique and sorted")
    required_shape = {"path", "collision_system", "kind", "sha256", "size_bytes"}
    for row in rows:
        if not isinstance(row, dict) or set(row) != required_shape:
            findings.append("fixture manifest row shape is invalid")
            continue
        resolved = _safe_tracked_file(row["path"], tracked, "fixture artifact", findings)
        if resolved is None:
            continue
        relative, target = resolved
        if row.get("collision_system") not in {"pp", "auau"}:
            findings.append(f"fixture collision system is invalid: {relative}")
        if row.get("kind") != "synthetic_public_tree_fixture":
            findings.append(f"fixture kind is invalid: {relative}")
        if not _valid_sha256(row.get("sha256")) or _sha256(target) != row.get("sha256"):
            findings.append(f"fixture artifact hash mismatch: {relative}")
        if target.stat().st_size != row.get("size_bytes"):
            findings.append(f"fixture artifact size mismatch: {relative}")
        _validate_root_fixture(target, contract, relative, findings)

    receipt_resolved = _safe_tracked_file(
        manifest.get("generation_receipt"), tracked, "fixture generation receipt", findings
    )
    if receipt_resolved is None:
        return
    _, receipt_path = receipt_resolved
    generation = _json_object(receipt_path, "fixture generation receipt", findings)
    if generation is None:
        return
    if generation.get("schema") != "PhotonJetEngineeringFixtureGenerationReceiptV1":
        findings.append("fixture generation receipt schema is invalid")
    if generation.get("fixtures") != rows:
        findings.append("fixture generation receipt differs from fixture manifest")
    if generation.get("public_regeneration") != "NOT_INCLUDED_IN_V1":
        findings.append("fixture regeneration status is not explicit")
    if generation.get("physics_use") != "FORBIDDEN_ENGINEERING_ONLY":
        findings.append("fixture physics-use boundary is not explicit")
    if generation.get("container_metadata_policy") != "NO_MACHINE_PATHS":
        findings.append("fixture container metadata policy is invalid")
    if generation.get("negative_metadata_scan") != "REQUIRED_BY_REPOSITORY_VALIDATOR":
        findings.append("fixture negative metadata scan is not required")


def _feature_cardinalities(findings: list[str]) -> dict[str, int]:
    path = REPO / "python/photonjet/training/features.py"
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"))
    except (OSError, SyntaxError) as error:
        findings.append(f"feature contract cannot be decoded: {error}")
        return {}
    result: dict[str, int] = {}
    for node in tree.body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        target = node.targets[0]
        if not isinstance(target, ast.Name) or target.id not in {"PP_H70_FEATURES", "AUAU_H70_FEATURES"}:
            continue
        try:
            values = ast.literal_eval(node.value)
        except (ValueError, TypeError):
            continue
        if not isinstance(values, tuple) or not values or len(values) != len(set(values)):
            findings.append(f"feature contract is invalid: {target.id}")
            continue
        result["pp" if target.id.startswith("PP_") else "auau"] = len(values)
    if set(result) != {"pp", "auau"}:
        findings.append("feature contract does not define both collision systems")
    return result


def _validate_model_manifests(tracked: set[str], findings: list[str]) -> None:
    feature_counts = _feature_cardinalities(findings)
    manifest_paths = sorted((REPO / "models/manifests").glob("*.json"))
    if [path.name for path in manifest_paths] != ["auau_h70.json", "pp_h70.json"]:
        findings.append("model manifest inventory is invalid")
    for manifest_path in manifest_paths:
        manifest = _json_object(manifest_path, f"model manifest {manifest_path.name}", findings)
        if manifest is None:
            continue
        system = manifest.get("system")
        expected_name = manifest_path.stem
        if manifest.get("schema") != "PhotonJetModelReferenceV1":
            findings.append(f"model manifest schema is invalid: {manifest_path.name}")
        if set(manifest) != {
            "schema", "name", "system", "artifact_filename", "artifact_sha256",
            "artifact_size_bytes", "feature_contract", "feature_count",
            "score_direction", "distribution", "usage", "inference_validation", "status",
        }:
            findings.append(f"model manifest row shape is invalid: {manifest_path.name}")
        if system not in {"pp", "auau"} or manifest.get("name") != expected_name:
            findings.append(f"model manifest identity is invalid: {manifest_path.name}")
        expected_contract = {
            "pp": "python/photonjet/training/features.py#PP_H70_FEATURES",
            "auau": "python/photonjet/training/features.py#AUAU_H70_FEATURES",
        }.get(str(system))
        if manifest.get("feature_contract") != expected_contract:
            findings.append(f"model feature-contract pointer is invalid: {manifest_path.name}")
        if manifest.get("feature_count") != feature_counts.get(str(system)):
            findings.append(f"model feature count differs: {manifest_path.name}")
        for key, expected in {
            "score_direction": "higher_is_signal",
            "distribution": "repository_hash_bound",
            "usage": "not_consumed_by_public_v1",
            "inference_validation": "NOT_RUN",
            "status": "HASH_BOUND_REFERENCE_ONLY",
        }.items():
            if manifest.get(key) != expected:
                findings.append(f"model manifest {key} is invalid: {manifest_path.name}")
        resolved = _safe_tracked_file(
            manifest.get("artifact_filename"), tracked, "model artifact", findings
        )
        if resolved is None:
            continue
        relative, artifact_path = resolved
        if not _valid_sha256(manifest.get("artifact_sha256")) or _sha256(artifact_path) != manifest.get("artifact_sha256"):
            findings.append(f"model artifact hash mismatch: {relative}")
        if artifact_path.stat().st_size != manifest.get("artifact_size_bytes"):
            findings.append(f"model artifact size mismatch: {relative}")
        if artifact_path.read_bytes()[:4] != b"root":
            findings.append(f"model artifact does not have ROOT magic: {relative}")
            continue
        try:
            import numpy as np
            import uproot

            with uproot.open(artifact_path) as root:
                if root.classnames(recursive=True) != {"myBDT;1": "TMVA::Experimental::RBDT"}:
                    findings.append(f"model ROOT inventory differs: {relative}")
                    continue
                model = root["myBDT"]
                cut_indices = np.asarray(model.member("fCutIndices"), dtype="uint32")
                expected_count = feature_counts.get(str(system))
                if expected_count is None or set(map(int, cut_indices)) != set(range(expected_count)):
                    findings.append(f"model RBDT feature indices differ: {relative}")
        except Exception as error:
            findings.append(f"model ROOT contract cannot be verified: {relative}: {error}")


REFERENCE_ARTIFACT_CLAIMS = {
    "pp_simulation_photonjet": ("c7c826f9e829b47c20aa4ca2e0455895787aa0bd8efb1096bd8a406f2d4c7821", 1000),
    "pp_simulation_inclusive": ("e4a7d07f3093b9019f47b042a8a8cc7366f729ed6fba03e66d13d927ba40fa2a", 1000),
    "auau_embedding_photonjet": ("fd7502eb6f3ee29173ffe7709da8a28589d1b8ed412ef85aa36dc884cee0e0ca", 1000),
    "auau_embedding_inclusive": ("d84c754898a900c4d245f4f50246346680e9bd2a7467079e9689142b273feccb", 1000),
}
BAYES_OUTPUT_CLAIMS = {
    "pp_1d": "71205444bdd9d9c25f933729a176afa2ec00e6940df0743daa3203d02392c4bb",
    "pp_2d": "558d1cbf894e8fa3cf4ced4f460f304e52b3c3ae8e0c74707530459bd16c6289",
    "auau_1d": "7648c7d06c0509951c01de6a6d5a365da3d5134e7b13e401ae63ec344eddf6fa",
    "auau_2d": "83a949a9ffef7cd31954fd0f54aa984d0ec0e019fdb51c77f6413b8cbc5ef2d4",
}
REQUIRED_NOT_RUN_GATES = {
    "Actual p+p and Au+Au mini-DST producer and direct-reference emission",
    "Golden producer-to-public branch-value differential from mini-DST inputs",
    "Data trigger, vertex, minimum-bias, and event-gate behavior",
    "Resolved full-analysis configuration consumption and generic per-stage run-manifest wiring",
    "Model inference from the packaged hash-only reference artifacts",
    "Independent model-score reevaluation and full training/working-point parity",
    "pTgamma-by-xJgamma Region A/C direct histograms with Sumw2",
    "pT-binned purity with leakage and closure corrections",
    "PhotonJetTrees-to-response construction, fakes, misses, and boundary categories",
    "Full unfolding toys, covariance, closure, refolding, and iteration selection",
    "Full-statistics or source-complete physics equivalence",
    "Clean sPHENIX environment rehearsal and collaboration transfer",
}


def _validate_release_matrix(branch_total: int, findings: list[str]) -> None:
    path = REPO / "RELEASE_MATRIX.json"
    matrix = _json_object(path, "release matrix", findings)
    if matrix is None:
        return
    if matrix.get("schema") != "PhotonJetReleaseMatrixV1":
        findings.append("release matrix schema is invalid")
    if matrix.get("release_state") != "DOWNSTREAM_DIFFERENTIAL_CANDIDATE":
        findings.append("release matrix state oversteps or differs from the V1 boundary")
    if matrix.get("claim_boundary") != (
        "Validated PhotonJetTrees_v1 through selected downstream artifacts; "
        "no mini-DST producer or complete physics-result claim"
    ):
        findings.append("release matrix claim boundary is invalid")

    references = matrix.get("reference_artifacts")
    if not isinstance(references, dict) or set(references) != set(REFERENCE_ARTIFACT_CLAIMS):
        findings.append("release matrix reference-artifact inventory is invalid")
    else:
        for name, (expected_sha256, expected_events) in REFERENCE_ARTIFACT_CLAIMS.items():
            row = references.get(name)
            if not isinstance(row, dict) or set(row) != {"sha256", "events"}:
                findings.append(f"release matrix reference row is invalid: {name}")
                continue
            if not _valid_sha256(row.get("sha256")) or row.get("sha256") != expected_sha256:
                findings.append(f"release matrix reference hash differs: {name}")
            if row.get("events") != expected_events:
                findings.append(f"release matrix reference event count differs: {name}")

    gates = matrix.get("passed_gates")
    if not isinstance(gates, list) or not all(isinstance(row, dict) for row in gates):
        findings.append("release matrix passed gates must be object rows")
        gates = []
    gate_names = [str(row.get("gate", "")) for row in gates]
    expected_gate_names = [
        "public_tree_contract",
        "xjgamma_histogram_differential",
        "recoil_skim_differential",
        "aggregate_abcd_differential",
        "iterative_bayes_kernel_differential",
        "plot_semantics_and_geometry",
        "redistribution_safe_tests",
        "wheel_install_contract",
    ]
    if gate_names != expected_gate_names or len(gate_names) != len(set(gate_names)):
        findings.append("release matrix passed-gate inventory is invalid")
    by_name = {str(row.get("gate", "")): row for row in gates}
    tree_gate = by_name.get("public_tree_contract", {})
    if tree_gate.get("status") != "PASS_FULL_PREVIEW" or f"{branch_total} ordered branch/type entries" not in str(tree_gate.get("scope", "")):
        findings.append("release matrix branch-contract claim differs from the executable contract")
    bayes_gate = by_name.get("iterative_bayes_kernel_differential", {})
    if bayes_gate.get("status") != "PASS_BYTE_EXACT" or bayes_gate.get("canonical_output_sha256") != BAYES_OUTPUT_CLAIMS:
        findings.append("release matrix iterative-Bayes claim is invalid")
    test_gate = by_name.get("redistribution_safe_tests", {})
    if test_gate.get("status") != "PASS" or re.search(r"\b[0-9]+\s+tests?\b", str(test_gate.get("scope", "")), re.IGNORECASE):
        findings.append("release matrix test claim is invalid or encodes an unbound raw count")
    for name in expected_gate_names:
        row = by_name.get(name, {})
        if not isinstance(row.get("scope"), str) or not row.get("scope", "").strip():
            findings.append(f"release matrix gate scope is missing: {name}")
        if not isinstance(row.get("status"), str) or not row.get("status", "").startswith("PASS"):
            findings.append(f"release matrix gate status is invalid: {name}")

    not_run = matrix.get("not_run_gates")
    if not isinstance(not_run, list) or not all(isinstance(item, str) for item in not_run):
        findings.append("release matrix NOT_RUN gates must be strings")
    elif set(not_run) != REQUIRED_NOT_RUN_GATES or len(not_run) != len(set(not_run)):
        findings.append("release matrix NOT_RUN gate inventory is incomplete or duplicated")
    evidence = matrix.get("evidence_policy")
    if evidence != {
        "collaborator_package": "Durable contracts, redistribution-safe fixtures, executable validators, and negative tests",
        "migration_checks": "Forensic comparators and site-specific receipts are retained outside the collaborator package",
    }:
        findings.append("release matrix evidence boundary is invalid")


def _validate_git_metadata(findings: list[str]) -> None:
    if not (REPO / ".git").exists():
        return
    try:
        branch = subprocess.run(
            ["git", "-C", str(REPO), "branch", "--show-current"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        history = subprocess.run(
            [
                "git", "-C", str(REPO), "log", "HEAD",
                "--format=%H%x00%an%x00%ae%x00%cn%x00%ce%x00%B%x1e",
            ],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
        tags = subprocess.run(
            ["git", "-C", str(REPO), "tag", "--merged", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.splitlines()
    except (OSError, subprocess.CalledProcessError) as error:
        findings.append(f"Git metadata cannot be inspected: {error}")
        return
    records = [item.lstrip("\n") for item in history.split("\x1e") if item.strip("\n")]
    if not records:
        findings.append("Git history has no reachable commit")
    for index, raw in enumerate(records):
        record = raw.split("\x00", 5)
        if len(record) != 6:
            findings.append(f"Git metadata record is malformed at reachable commit {index}")
            continue
        commit, _, _, _, _, message = record
        for label in sorted(_sensitive_labels([message]) - {"private contact identifier"}):
            findings.append(f"{label} in reachable Git commit {commit[:12]}")
    for surface, value in (("branch", branch), *(("tag", tag) for tag in tags)):
        labels = _sensitive_labels([value]) - {"private contact identifier"}
        for label in sorted(labels):
            findings.append(f"{label} in Git {surface}")
    _validate_git_history_trees([record.split("\x00", 1)[0] for record in records], findings)


def _historical_blob_findings(relative: str, data: bytes) -> list[str]:
    suffix = PurePosixPath(relative).suffix
    if suffix in TEXT_SUFFIXES:
        try:
            text = unicodedata.normalize("NFKC", data.decode("utf-8"))
        except UnicodeDecodeError:
            return ["historical text cannot be decoded as UTF-8"]
        return sorted(
            label for label, pattern in SENSITIVE_PATTERNS.items() if pattern.search(text)
        )
    with tempfile.NamedTemporaryFile(suffix=suffix) as handle:
        handle.write(data)
        handle.flush()
        return _binary_sensitive_findings(Path(handle.name))


def _validate_git_history_trees(commits: list[str], findings: list[str]) -> None:
    blob_findings: dict[str, list[str]] = {}
    for commit in commits:
        try:
            tree = subprocess.run(
                ["git", "-C", str(REPO), "ls-tree", "-rz", "--full-tree", commit],
                check=True,
                capture_output=True,
            ).stdout
        except (OSError, subprocess.CalledProcessError) as error:
            findings.append(f"reachable Git tree cannot be inspected: {commit[:12]}: {error}")
            continue
        for encoded in tree.split(b"\x00"):
            if not encoded:
                continue
            try:
                metadata, raw_path = encoded.split(b"\t", 1)
                mode, object_type, blob = metadata.decode("ascii").split(" ", 2)
                relative = unicodedata.normalize("NFKC", raw_path.decode("utf-8"))
            except (UnicodeDecodeError, ValueError):
                findings.append(f"reachable Git tree row is malformed: {commit[:12]}")
                continue
            for label in sorted(_sensitive_labels([relative])):
                findings.append(f"{label} in reachable Git path: {commit[:12]}:{relative}")
            if mode == "120000":
                findings.append(f"symbolic link in reachable Git tree: {commit[:12]}:{relative}")
            if object_type != "blob":
                continue
            if blob not in blob_findings:
                try:
                    data = subprocess.run(
                        ["git", "-C", str(REPO), "cat-file", "blob", blob],
                        check=True,
                        capture_output=True,
                    ).stdout
                    blob_findings[blob] = _historical_blob_findings(relative, data)
                except (OSError, subprocess.CalledProcessError) as error:
                    blob_findings[blob] = [f"historical Git blob cannot be read: {error}"]
            for label in blob_findings[blob]:
                findings.append(f"{label} in reachable Git content: {commit[:12]}:{relative}")


def _validate_projection_receipt(receipt_path: Path, findings: list[str]) -> set[str]:
    try:
        relative_receipt = receipt_path.resolve().relative_to(REPO.resolve()).as_posix()
    except ValueError:
        findings.append("projection receipt is outside the projected repository")
        return set()
    try:
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        findings.append(f"invalid projection receipt: {error}")
        return set()
    if receipt.get("schema") != "PhotonJetProjectionReceiptV2":
        findings.append("projection receipt schema is not PhotonJetProjectionReceiptV2")
    projection_target = receipt.get("target")
    expected_projection_paths = {
        "analysis": "PhotonJet/IsolatedPhotonJetAnalysis",
        "internal": "canonical/photonjet",
    }
    if projection_target not in expected_projection_paths:
        findings.append("projection receipt target is invalid")
    elif receipt.get("projection_path") != expected_projection_paths[projection_target]:
        findings.append("projection receipt path differs from its target")
    if receipt.get("source_tree_clean") is not True:
        findings.append("projection receipt is not bound to a clean source checkpoint")
    if not re.fullmatch(r"[0-9a-f]{40,64}", str(receipt.get("source_commit", ""))):
        findings.append("projection receipt source commit is invalid")
    rows = receipt.get("files")
    if not isinstance(rows, list):
        findings.append("projection receipt files must be a list")
        return set()
    paths: list[str] = []
    normalized_rows: list[dict[str, Any]] = []
    for row in rows:
        if not isinstance(row, dict):
            findings.append("projection receipt contains a non-object file row")
            continue
        relative = str(row.get("path", ""))
        candidate = Path(relative)
        if not relative or candidate.is_absolute() or ".." in candidate.parts:
            findings.append(f"unsafe path in projection receipt: {relative!r}")
            continue
        paths.append(relative)
        file_target = REPO / candidate
        if not file_target.is_file() or file_target.is_symlink():
            findings.append(f"projection receipt target is missing or not a regular file: {relative}")
            continue
        observed_size = file_target.stat().st_size
        observed_sha256 = _sha256(file_target)
        expected_mode = row.get("mode")
        if expected_mode not in {"100644", "100755"}:
            findings.append(f"projection receipt mode is invalid: {relative}")
        else:
            observed_mode = stat.S_IMODE(file_target.stat().st_mode)
            required_mode = 0o755 if expected_mode == "100755" else 0o644
            if observed_mode != required_mode:
                findings.append(
                    f"projection receipt mode mismatch: {relative} "
                    f"{observed_mode:04o} != {required_mode:04o}"
                )
        if observed_size != row.get("size_bytes"):
            findings.append(f"projection receipt size mismatch: {relative}")
        if observed_sha256 != row.get("sha256"):
            findings.append(f"projection receipt hash mismatch: {relative}")
        normalized_rows.append(
            {
                "path": relative,
                "mode": expected_mode,
                "sha256": row.get("sha256"),
                "size_bytes": row.get("size_bytes"),
            }
        )
    if paths != sorted(paths) or len(paths) != len(set(paths)):
        findings.append("projection receipt paths are not unique and sorted")
    if receipt.get("file_count") != len(rows):
        findings.append("projection receipt file count mismatch")
    content_sha256 = hashlib.sha256(
        json.dumps(normalized_rows, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    if content_sha256 != receipt.get("projection_content_sha256"):
        findings.append("projection receipt content digest mismatch")
    if not re.fullmatch(r"[0-9a-f]{64}", str(receipt.get("projection_manifest_sha256", ""))):
        findings.append("projection manifest digest is invalid")
    actual = {
        path.relative_to(REPO).as_posix()
        for path in _public_files()
        if path.resolve() != receipt_path.resolve()
    }
    expected = set(paths)
    always_forbidden = {".gitignore", "docs/draft_pull_request.md"}
    analysis_forbidden = always_forbidden | {
        "docs/analysis-integration.md",
        "docs/internal_integration.md",
        "manifests/projections.json",
        "tools/project_release.py",
    }
    internal_forbidden = always_forbidden | {"docs/analysis-integration.md"}
    internal_required = {
        "docs/internal_integration.md",
        "manifests/projections.json",
        "tools/project_release.py",
    }
    forbidden = analysis_forbidden if projection_target == "analysis" else internal_forbidden
    for extra in sorted(expected & forbidden):
        findings.append(f"forbidden file for projection target {projection_target}: {extra}")
    if projection_target == "internal":
        for missing in sorted(internal_required - expected):
            findings.append(f"required internal projection file is absent: {missing}")
        projection_manifest = REPO / "manifests/projections.json"
        if projection_manifest.is_file() and _sha256(projection_manifest) != receipt.get("projection_manifest_sha256"):
            findings.append("internal projection manifest hash differs from its receipt")
    for missing in sorted(expected - actual):
        findings.append(f"projection receipt file is absent: {missing}")
    for extra in sorted(actual - expected):
        findings.append(f"unreceipted file in projection: {extra}")
    if relative_receipt != "PROJECTION_RECEIPT.json":
        findings.append("projection receipt must be at repository root")
    return expected


EXPECTED_PROJECTION_MANIFEST = {
    "schema": "PhotonJetProjectionManifestV1",
    "include": [
        "README.md", "RELEASE_MATRIX.json", "CHANGELOG.md", "CONTRIBUTING.md",
        "LICENSE", "pyproject.toml", "requirements-constraints.txt", "bin", "config",
        "contracts", "docs", "examples", "macros", "manifests", "models", "offline",
        "python", "src", "tests", "tools", "workflows",
    ],
    "exclude": ["**/__pycache__/**", "**/.DS_Store", "docs/draft_pull_request.md"],
    "targets": {
        "analysis": {
            "path": "PhotonJet/IsolatedPhotonJetAnalysis",
            "exclude": [
                "docs/analysis-integration.md",
                "docs/internal_integration.md",
                "manifests/projections.json",
                "tools/project_release.py",
            ],
        },
        "internal": {
            "path": "canonical/photonjet",
            "exclude": ["docs/analysis-integration.md"],
        },
    },
}


def _projection_selected_paths(
    manifest: dict[str, Any], tracked: set[str], target: str
) -> set[str]:
    includes = manifest["include"]
    excludes = [*manifest["exclude"], *manifest["targets"][target]["exclude"]]
    result: set[str] = set()
    for relative in tracked:
        included = any(relative == item or relative.startswith(item + "/") for item in includes)
        if not included:
            continue
        path = PurePosixPath(relative)
        if any(path.match(pattern) for pattern in excludes):
            continue
        result.add(relative)
    return result


def _validate_projection_manifest(tracked: set[str], findings: list[str]) -> None:
    path = REPO / "manifests/projections.json"
    manifest = _json_object(path, "projection manifest", findings)
    if manifest is None:
        return
    if manifest != EXPECTED_PROJECTION_MANIFEST:
        findings.append("projection manifest differs from the canonical target contract")
        return
    for collection_name in ("include", "exclude"):
        rows = manifest[collection_name]
        if len(rows) != len(set(rows)):
            findings.append(f"projection manifest {collection_name} rows are duplicated")
    for target, spec in manifest["targets"].items():
        candidate = Path(spec["path"])
        if candidate.is_absolute() or ".." in candidate.parts:
            findings.append(f"unsafe projection target path: {target}")
        if len(spec["exclude"]) != len(set(spec["exclude"])):
            findings.append(f"projection target exclusions are duplicated: {target}")
    selected = {
        target: _projection_selected_paths(manifest, tracked, target)
        for target in ("analysis", "internal")
    }
    intended_internal_extras = {
        "docs/internal_integration.md",
        "manifests/projections.json",
        "tools/project_release.py",
    }
    if selected["internal"] - selected["analysis"] != intended_internal_extras:
        findings.append("internal projection extras differ from the canonical target contract")
    if selected["analysis"] - selected["internal"]:
        findings.append("analysis projection contains files absent from the internal projection")
    for target, paths in selected.items():
        for relative in sorted(paths):
            _safe_tracked_file(relative, tracked, f"{target} projection source", findings)


def validate(projection_receipt: Path | None = None) -> dict[str, Any]:
    findings: list[str] = []
    receipt_tracked = (
        _validate_projection_receipt(projection_receipt, findings)
        if projection_receipt is not None
        else None
    )
    files = _public_files()
    if receipt_tracked is not None:
        tracked = receipt_tracked
    else:
        try:
            tracked = set(
                subprocess.run(
                    ["git", "-C", str(REPO), "ls-files"],
                    check=True,
                    capture_output=True,
                    text=True,
                ).stdout.splitlines()
            )
        except (OSError, subprocess.CalledProcessError) as error:
            findings.append(f"release tracking metadata cannot be inspected: {error}")
            tracked = set()
    for path in files:
        relative = unicodedata.normalize("NFKC", path.relative_to(REPO).as_posix())
        if path.is_symlink():
            findings.append(f"symbolic link is forbidden in release projection: {relative}")
            continue
        for label, pattern in SENSITIVE_PATTERNS.items():
            if pattern.search(relative):
                findings.append(f"{label} in path: {relative}")
        text = _text(path)
        if text is None:
            for label in _binary_sensitive_findings(path):
                findings.append(f"{label} in binary: {relative}")
            continue
        for label, pattern in SENSITIVE_PATTERNS.items():
            match = pattern.search(text)
            if match:
                line = text.count("\n", 0, match.start()) + 1
                findings.append(f"{label}: {relative}:{line}")
        is_plot_code = (
            "plot" in relative.lower()
            or relative.startswith("macros/")
            or relative.startswith("examples/")
        ) and path.suffix in {".py", ".C", ".cc", ".cpp", ".h", ".hpp"}
        if is_plot_code and relative not in PLOT_CODE_EXCEPTIONS:
            python_match = _python_raw_annotation(text) if path.suffix == ".py" else None
            match = RAW_PLOT_SEMANTIC.search(text) if python_match is None else None
            if python_match is not None:
                line = python_match[0]
                findings.append(f"raw plot physics annotation: {relative}:{line}")
            elif match:
                line = text.count("\n", 0, match.start()) + 1
                findings.append(f"raw plot physics annotation: {relative}:{line}")
            if "#it{#bf{sPHENIX}}" in text:
                line = text[: text.index("#it{#bf{sPHENIX}}")].count("\n") + 1
                findings.append(f"experiment label outside canonical style source: {relative}:{line}")

    branch_contract, branch_total = _validate_branch_contract(findings)
    _validate_fixture_manifest(tracked, branch_contract, findings)
    _validate_model_manifests(tracked, findings)
    _validate_release_matrix(branch_total, findings)
    _validate_git_metadata(findings)

    if projection_receipt is None:
        _validate_projection_manifest(tracked, findings)

    cli = (REPO / "python/photonjet/cli.py").read_text(encoding="utf-8")
    for forbidden in ("--cut-label", "--physics-label", "--sample-label", "--sphenix-label"):
        if forbidden in cli:
            findings.append(f"free-form plot-label CLI option is forbidden: {forbidden}")

    for required in (
        "tests/fixtures/photonjet_trees_pp.root",
        "tests/fixtures/photonjet_trees_auau.root",
        "models/artifacts/pp_h70.root",
        "models/artifacts/auau_h70.root",
    ):
        if required not in tracked:
            findings.append(f"required release artifact is not tracked: {required}")

    return {
        "schema": "PhotonJetRepositoryValidationV1",
        "status": "PASS" if not findings else "FAIL",
        "file_count": len(files),
        "checks": {
            "public_boundary": "PASS" if not any(
                any(label in item for label in SENSITIVE_PATTERNS) for item in findings
            ) else "FAIL",
            "hash_bound_fixtures_and_models": "PASS" if not any(
                any(token in item for token in (
                    "fixture", "model", "branch contract", "branch row", "branch names",
                    "required release artifact", "feature contract",
                ))
                for item in findings
            ) else "FAIL",
            "plot_annotation_boundary": "PASS" if not any(
                item.startswith(("raw plot physics annotation", "experiment label", "free-form plot-label"))
                for item in findings
            ) else "FAIL",
            "projection_surface": "PASS" if not any("projection" in item for item in findings) else "FAIL",
            "release_claims": "PASS" if not any("release matrix" in item for item in findings) else "FAIL",
            "git_metadata": "PASS" if not any("Git " in item for item in findings) else "FAIL",
        },
        "findings": findings,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--projection-receipt",
        type=Path,
        help="validate a generated projection against its root receipt instead of source Git metadata",
    )
    args = parser.parse_args()
    result = validate(args.projection_receipt)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    sys.stdout.write(rendered)
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
