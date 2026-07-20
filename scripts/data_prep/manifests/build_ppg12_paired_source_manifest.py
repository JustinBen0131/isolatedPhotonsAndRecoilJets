#!/usr/bin/env python3
"""Build the exact 12-lane PPG12 paired-oracle source manifest.

The paired-oracle campaign driver consumes a compact path-only JSON document.
This builder is the fail-closed bridge from one sealed new.17 runtime plus an
explicit 12-row lane table to that document.  It validates every runtime file
against the SHA-256 recorded by the runtime builder, derives the five model
assets from their sealed roles, and refuses incomplete or duplicated lane
coverage.

The lane TSV must have this exact header::

    lane_id  sample  period  interaction  ppg_macro  g4_full_list  truthjet_full_list

Paths in both the TSV and CLI must already be absolute.  The output is written
with canonical lane ordering and deterministic JSON serialization.  Re-running
against an identical existing output is a no-op; overwriting different content
is refused.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
from pathlib import Path
from typing import Any


SCHEMA = "ppg12-paired-source-manifest/v1"
RUNTIME_PROFILE = "new.17"
RUNTIME_SCHEMA_VERSION = 1
ESTIMATOR_REVISION = "29f8223bd9b36dffab07961b597afa94185bbdf1"
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
GIT_COMMIT_RE = re.compile(r"^[0-9a-f]{40}$")

LANE_COLUMNS = (
    "lane_id",
    "sample",
    "period",
    "interaction",
    "ppg_macro",
    "g4_full_list",
    "truthjet_full_list",
)
LANE_PATH_FIELDS = ("ppg_macro", "g4_full_list", "truthjet_full_list")
COMMON_ROLE_FIELDS = {
    "apply_bdt": "ppg_apply_bdt_macro",
    "apply_config": "ppg_apply_bdt_config",
    "base_e_model": "ppg_apply_model_base_E",
    "base_v3e_model": "ppg_apply_model_base_v3E",
    "npb_model": "ppg_apply_npb_model",
}

# These are the roles the paired worker and the downstream closure contract
# actually consume.  Extra sealed roles are allowed, but omission of any of
# these makes the runtime unsuitable for the campaign.
REQUIRED_RUNTIME_ROLES = frozenset(
    {
        "recoil_macro",
        "recoil_impl",
        "libRecoilJets.so",
        "libCaloAna24.so",
        "libcalo_reco.so",
        "libclusteriso.so",
        "libjetbase.so",
        "PhotonClusterBuilder.h",
        "ppg_recoeff_source_macro",
        "ppg_recoeff_macro",
        "ppg_recoeff_trace_macro",
        "ppg_recoeff_trace_transform_receipt",
        "ppg_recoeff_cross_section_header",
        "ppg_recoeff_truth_vertex_header",
        "ppg_recoeff_period_config_0mrad",
        "ppg_recoeff_period_config_1p5mrad",
        "ppg_recoeff_truth_vertex_reweight_0mrad",
        "ppg_recoeff_truth_vertex_reweight_1p5mrad",
        "ppg_calculate_photon_yield",
        "ppg_apply_bdt_macro",
        "ppg_apply_bdt_config",
        "ppg_recoeff_yaml_cpp",
        "ppg_recoeff_yaml_cpp_header_tree_receipt",
        "ppg_recoeff_roounfold",
        "ppg_recoeff_roounfold_response_header",
        "ppg_recoeff_roounfold_bayes_header",
        "ppg_recoeff_vertex_scan_data",
        "ppg_recoeff_mbd_correction",
        "ppg_apply_model_base_E",
        "ppg_apply_model_base_v3E",
        "ppg_apply_npb_model",
    }
)


class ManifestBuildError(RuntimeError):
    """The requested manifest cannot be built without weakening provenance."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _require_sha256(value: str, label: str) -> str:
    normalized = value.strip().lower()
    if not SHA256_RE.fullmatch(normalized):
        raise ManifestBuildError(f"{label} must be one lowercase SHA-256 digest")
    return normalized


def _require_git_commit(value: str, label: str) -> str:
    normalized = value.strip().lower()
    if not GIT_COMMIT_RE.fullmatch(normalized):
        raise ManifestBuildError(f"{label} must be one lowercase 40-hex Git commit")
    return normalized


def _absolute_file(value: object, label: str) -> Path:
    if not isinstance(value, (str, Path)) or not str(value).strip():
        raise ManifestBuildError(f"{label} must be a non-empty absolute path")
    if str(value) != str(value).strip() or "\n" in str(value) or "\r" in str(value):
        raise ManifestBuildError(f"{label} contains unsafe whitespace: {value!r}")
    raw = Path(value)
    if not raw.is_absolute():
        raise ManifestBuildError(f"{label} must be absolute: {raw}")
    try:
        resolved = raw.resolve(strict=True)
    except FileNotFoundError as exc:
        raise ManifestBuildError(f"{label} is missing: {raw}") from exc
    if not resolved.is_file() or resolved.stat().st_size <= 0:
        raise ManifestBuildError(f"{label} must be a non-empty file: {resolved}")
    return resolved


def _read_json(path: Path, label: str) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ManifestBuildError(f"{label} is invalid JSON: {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ManifestBuildError(f"{label} must contain one JSON object: {path}")
    return payload


def _expected_lanes() -> list[tuple[str, str, str, str]]:
    return [
        (
            f"photon:photon{photon}:{period}:{interaction}",
            f"Photon{photon}",
            period,
            interaction.upper(),
        )
        for photon in (5, 10, 20)
        for period in ("0mrad", "1p5mrad")
        for interaction in ("si", "di")
    ]


def _runtime_roles(
    runtime_manifest: Path,
    expected_runtime_sha256: str,
    expected_estimator_revision: str,
) -> dict[str, Path]:
    observed_manifest_sha = _sha256(runtime_manifest)
    if observed_manifest_sha != expected_runtime_sha256:
        raise ManifestBuildError(
            "runtime manifest SHA-256 differs from the explicitly authorized attempt: "
            f"expected={expected_runtime_sha256}, observed={observed_manifest_sha}"
        )
    payload = _read_json(runtime_manifest, "runtime manifest")
    if payload.get("schema_version") != RUNTIME_SCHEMA_VERSION:
        raise ManifestBuildError("runtime manifest schema_version must be 1")
    if payload.get("runtime_profile") != RUNTIME_PROFILE:
        raise ManifestBuildError("runtime manifest must seal the new.17 profile")
    if payload.get("isolated_build") is not True:
        raise ManifestBuildError("runtime manifest must record isolated_build=true")
    if payload.get("estimator_revision") != expected_estimator_revision:
        raise ManifestBuildError(
            "runtime estimator revision differs from the closure contract: "
            f"expected={expected_estimator_revision}, "
            f"observed={payload.get('estimator_revision')!r}"
        )

    receipt = _absolute_file(payload.get("build_receipt", ""), "runtime build receipt")
    receipt_sha = _require_sha256(
        str(payload.get("build_receipt_sha256", "")),
        "runtime build_receipt_sha256",
    )
    if _sha256(receipt) != receipt_sha:
        raise ManifestBuildError("runtime build receipt SHA-256 is stale")

    raw_files = payload.get("files")
    if not isinstance(raw_files, list) or not raw_files:
        raise ManifestBuildError("runtime manifest files must be a non-empty list")
    roles: dict[str, Path] = {}
    # The preserved attempt-6 photon sources are period-independent: the
    # 0mrad and 1p5mrad lanes for one sample/interaction deliberately reuse
    # the same macro and paired event lists.  Reuse is valid only for that
    # exact cross-period pair.  Any sharing across samples, interactions, or
    # source roles remains a hard failure.
    seen_paths: dict[Path, tuple[tuple[str, str], str, str]] = {}
    for index, row in enumerate(raw_files):
        if not isinstance(row, dict):
            raise ManifestBuildError(f"runtime files[{index}] must be an object")
        raw_role = row.get("role")
        if not isinstance(raw_role, str) or not raw_role.strip():
            raise ManifestBuildError(f"runtime files[{index}].role is missing")
        role = raw_role.strip()
        if role != raw_role:
            raise ManifestBuildError(
                f"runtime files[{index}].role contains surrounding whitespace"
            )
        if role in roles:
            raise ManifestBuildError(f"runtime manifest duplicates role {role}")
        path = _absolute_file(row.get("path", ""), f"runtime role {role}")
        if path in seen_paths:
            raise ManifestBuildError(
                f"runtime roles {seen_paths[path]} and {role} duplicate path {path}"
            )
        declared_sha = _require_sha256(
            str(row.get("sha256", "")), f"runtime role {role} sha256"
        )
        observed_sha = _sha256(path)
        if declared_sha != observed_sha:
            raise ManifestBuildError(
                f"runtime role {role} SHA-256 is stale: "
                f"expected={declared_sha}, observed={observed_sha}"
            )
        roles[role] = path
        seen_paths[path] = role

    missing_roles = sorted(REQUIRED_RUNTIME_ROLES - set(roles))
    if missing_roles:
        raise ManifestBuildError(
            "runtime manifest lacks paired-oracle roles: " + ", ".join(missing_roles)
        )
    return roles


def _source_rows(path: Path, label: str) -> list[str]:
    rows = [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if len(rows) < 5:
        raise ManifestBuildError(f"{label} must contain at least five source rows")
    if len(rows) != len(set(rows)):
        raise ManifestBuildError(f"{label} contains duplicate source rows")
    return rows


def _normalized_event_identity(row: str, kind: str, label: str) -> str:
    prefixes = {
        "g4": "G4Hits_",
        "truthjet": "DST_TRUTH_JET_",
    }
    prefix = prefixes[kind]
    if any(char.isspace() for char in row):
        raise ManifestBuildError(f"{label} contains whitespace: {row!r}")
    path = Path(row)
    basename = path.name
    # The preserved lists use relative basenames.  Absolute rows are also
    # valid, but relative directory traversal is not part of this contract.
    if not path.is_absolute() and str(path) != basename:
        raise ManifestBuildError(
            f"{label} must be an absolute path or preserved basename: {row}"
        )
    if not basename.startswith(prefix) or not basename.endswith(".root"):
        raise ManifestBuildError(
            f"{label} must use the recognized {prefix}*.root naming contract: {row}"
        )
    identity = basename.removeprefix(prefix)
    if not identity or identity == ".root":
        raise ManifestBuildError(f"{label} has an empty event identity: {row}")
    return identity


def _paired_source_contract(
    g4_rows: list[str], truth_rows: list[str], label: str
) -> str:
    if len(g4_rows) != len(truth_rows):
        raise ManifestBuildError(f"{label}: G4 and TRUTH_JET lists have different row counts")
    g4_ids = [
        _normalized_event_identity(row, "g4", f"{label} G4 row {index}")
        for index, row in enumerate(g4_rows, start=1)
    ]
    truth_ids = [
        _normalized_event_identity(row, "truthjet", f"{label} TRUTH_JET row {index}")
        for index, row in enumerate(truth_rows, start=1)
    ]
    if len(g4_ids) != len(set(g4_ids)) or len(truth_ids) != len(set(truth_ids)):
        raise ManifestBuildError(f"{label}: normalized event identities must be unique")
    for index, (g4_id, truth_id) in enumerate(zip(g4_ids, truth_ids), start=1):
        if g4_id != truth_id:
            raise ManifestBuildError(
                f"{label}: event identity mismatch at row {index}: "
                f"G4={g4_id!r}, TRUTH_JET={truth_id!r}"
            )
    first_five = "\n".join(g4_ids[:5]).encode()
    paired_first_five = "\n".join(truth_ids[:5]).encode()
    g4_hash = hashlib.sha256(first_five).hexdigest()
    truth_hash = hashlib.sha256(paired_first_five).hexdigest()
    if g4_hash != truth_hash:
        raise ManifestBuildError(f"{label}: first-five event identity hash differs")
    return g4_hash


def _active_cpp(text: str) -> str:
    without_blocks = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
    return re.sub(r"//[^\n]*", "", without_blocks)


def _validate_ppg_macro(path: Path, interaction: str, label: str) -> None:
    active = _active_cpp(path.read_text(encoding="utf-8"))
    indices = {
        int(value)
        for value in re.findall(
            r"INPUTREADHITS\s*::\s*listfile\s*\[\s*([0-9]+)\s*\]\s*=",
            active,
        )
    }
    if indices != {0, 4}:
        raise ManifestBuildError(
            f"{label}: active INPUTREADHITS list indices must be exactly {{0, 4}}, "
            f"observed={sorted(indices)}"
        )
    embedding_two = bool(
        re.search(r"\badd_embedding_flag\s*\(\s*2\s*\)", active)
    )
    truth_jet_reconstruction = bool(
        re.search(r"\bnew\s+TruthJetInput\b", active)
    )
    if interaction == "DI":
        if not embedding_two:
            raise ManifestBuildError(
                f"{label}: DI macro must actively call add_embedding_flag(2)"
            )
        if not truth_jet_reconstruction:
            raise ManifestBuildError(
                f"{label}: DI macro must actively reconstruct truth jets with TruthJetInput"
            )
    elif embedding_two:
        raise ManifestBuildError(
            f"{label}: SI macro must not actively call add_embedding_flag(2)"
        )


def _lane_rows(lane_tsv: Path) -> list[dict[str, str]]:
    try:
        with lane_tsv.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            if tuple(reader.fieldnames or ()) != LANE_COLUMNS:
                raise ManifestBuildError(
                    "lane TSV header must be exact: " + "\t".join(LANE_COLUMNS)
                )
            raw_rows = list(reader)
    except csv.Error as exc:
        raise ManifestBuildError(f"lane TSV is malformed: {lane_tsv}: {exc}") from exc

    expected = _expected_lanes()
    expected_by_id = {
        lane_id: (sample, period, interaction)
        for lane_id, sample, period, interaction in expected
    }
    by_id: dict[str, dict[str, str]] = {}
    seen_paths: dict[Path, tuple[tuple[str, str], str, str]] = {}
    normalized_rows: dict[str, dict[str, str]] = {}
    for index, row in enumerate(raw_rows, start=2):
        normalized = {key: str(row.get(key, "")).strip() for key in LANE_COLUMNS}
        lane_id = normalized["lane_id"]
        if lane_id in by_id:
            raise ManifestBuildError(f"lane TSV duplicates lane_id {lane_id}")
        by_id[lane_id] = normalized
        if lane_id not in expected_by_id:
            raise ManifestBuildError(f"lane TSV contains noncanonical lane_id {lane_id!r}")
        expected_identity = expected_by_id[lane_id]
        observed_identity = (
            normalized["sample"],
            normalized["period"],
            normalized["interaction"],
        )
        if observed_identity != expected_identity:
            raise ManifestBuildError(
                f"{lane_id}: physical identity mismatch: "
                f"expected={expected_identity}, observed={observed_identity}"
            )
        normalized_paths: dict[str, str] = {}
        list_rows: dict[str, list[str]] = {}
        for field in LANE_PATH_FIELDS:
            path = _absolute_file(normalized[field], f"{lane_id} {field}")
            source_identity = (normalized["sample"], normalized["interaction"])
            prior = seen_paths.get(path)
            if prior is not None:
                prior_identity, prior_field, prior_lane = prior
                if prior_identity != source_identity or prior_field != field:
                    raise ManifestBuildError(
                        f"{lane_id} {field} illegally shares {prior_lane} "
                        f"{prior_field} path {path}"
                    )
            else:
                seen_paths[path] = (source_identity, field, lane_id)
            normalized_paths[field] = str(path)
            if field.endswith("_list"):
                list_rows[field] = _source_rows(path, f"{lane_id} {field}")
        _paired_source_contract(
            list_rows["g4_full_list"], list_rows["truthjet_full_list"], lane_id
        )
        _validate_ppg_macro(
            Path(normalized_paths["ppg_macro"]), normalized["interaction"], lane_id
        )
        normalized_rows[lane_id] = {
            "lane_id": lane_id,
            "sample": normalized["sample"],
            "period": normalized["period"],
            "interaction": normalized["interaction"],
            **normalized_paths,
        }

    missing = sorted(set(expected_by_id) - set(by_id))
    if missing or len(raw_rows) != len(expected):
        raise ManifestBuildError(
            "lane TSV must contain the exact 12 lanes; missing="
            + repr(missing)
            + f", rows={len(raw_rows)}"
        )

    for photon in (5, 10, 20):
        sample = f"Photon{photon}"
        for interaction in ("SI", "DI"):
            zero = normalized_rows[f"photon:photon{photon}:0mrad:{interaction.lower()}"]
            shifted = normalized_rows[
                f"photon:photon{photon}:1p5mrad:{interaction.lower()}"
            ]
            for field in LANE_PATH_FIELDS:
                if zero[field] != shifted[field]:
                    raise ManifestBuildError(
                        f"{sample} {interaction} {field} must be exactly shared "
                        "between 0mrad and 1p5mrad"
                    )

        si = normalized_rows[f"photon:photon{photon}:0mrad:si"]
        di = normalized_rows[f"photon:photon{photon}:0mrad:di"]
        for field in ("g4_full_list", "truthjet_full_list"):
            si_rows = _source_rows(Path(si[field]), f"{sample} SI {field}")
            di_rows = _source_rows(Path(di[field]), f"{sample} DI {field}")
            if si_rows == di_rows:
                raise ManifestBuildError(
                    f"{sample} SI and DI {field} source rows must be distinct"
                )
    return [normalized_rows[lane_id] for lane_id, *_ in expected]


def build_manifest(
    *,
    runtime_manifest: Path,
    expected_runtime_sha256: str,
    setup_script: Path,
    tower_mask: Path,
    recoil_config: Path,
    lane_tsv: Path,
    expected_estimator_revision: str = ESTIMATOR_REVISION,
) -> dict[str, Any]:
    runtime_manifest = _absolute_file(runtime_manifest, "runtime manifest")
    lane_tsv = _absolute_file(lane_tsv, "lane TSV")
    expected_runtime_sha256 = _require_sha256(
        expected_runtime_sha256, "expected runtime manifest SHA-256"
    )
    expected_estimator_revision = _require_git_commit(
        expected_estimator_revision, "expected estimator revision"
    )
    roles = _runtime_roles(
        runtime_manifest,
        expected_runtime_sha256,
        expected_estimator_revision,
    )

    explicit = {
        "setup_script": _absolute_file(setup_script, "setup script"),
        "tower_mask": _absolute_file(tower_mask, "tower mask"),
        "recoil_config": _absolute_file(recoil_config, "RecoilJets config"),
    }
    common_paths: dict[str, Path] = {
        **explicit,
        **{field: roles[role] for field, role in COMMON_ROLE_FIELDS.items()},
        "recoil_runtime_manifest": runtime_manifest,
    }
    reverse: dict[Path, str] = {}
    for field, path in common_paths.items():
        if path in reverse:
            raise ManifestBuildError(
                f"common assets {reverse[path]} and {field} duplicate path {path}"
            )
        reverse[path] = field

    lanes = _lane_rows(lane_tsv)
    lane_paths = {
        Path(lane[field])
        for lane in lanes
        for field in LANE_PATH_FIELDS
    }
    overlap = sorted(str(path) for path in lane_paths & set(common_paths.values()))
    if overlap:
        raise ManifestBuildError(
            "lane assets overlap common assets: " + ", ".join(overlap)
        )

    common_order = (
        "setup_script",
        "apply_bdt",
        "apply_config",
        "base_e_model",
        "base_v3e_model",
        "npb_model",
        "tower_mask",
        "recoil_runtime_manifest",
        "recoil_config",
    )
    return {
        "schema": SCHEMA,
        "common": {field: str(common_paths[field]) for field in common_order},
        "lanes": lanes,
    }


def _serialized(payload: dict[str, Any]) -> str:
    return json.dumps(
        payload,
        indent=2,
        sort_keys=True,
        ensure_ascii=False,
        allow_nan=False,
    ) + "\n"


def write_manifest(path: Path, payload: dict[str, Any]) -> str:
    if not path.is_absolute():
        raise ManifestBuildError(f"output must be absolute: {path}")
    rendered = _serialized(payload)
    if path.exists():
        if not path.is_file():
            raise ManifestBuildError(f"output exists and is not a file: {path}")
        if path.read_text(encoding="utf-8") == rendered:
            return "UNCHANGED"
        raise ManifestBuildError(f"refusing to overwrite different output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        with temporary.open("x", encoding="utf-8") as stream:
            stream.write(rendered)
            stream.flush()
            os.fsync(stream.fileno())
        temporary.replace(path)
    finally:
        if temporary.exists():
            temporary.unlink()
    return "WROTE"


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runtime-manifest", required=True, type=Path)
    parser.add_argument("--expected-runtime-sha256", required=True)
    parser.add_argument("--setup-script", required=True, type=Path)
    parser.add_argument("--tower-mask", required=True, type=Path)
    parser.add_argument("--recoil-config", required=True, type=Path)
    parser.add_argument("--lane-tsv", required=True, type=Path)
    parser.add_argument(
        "--expected-estimator-revision",
        default=ESTIMATOR_REVISION,
        help="frozen PPG12 estimator Git revision (default: closure contract revision)",
    )
    parser.add_argument("--output", required=True, type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        payload = build_manifest(
            runtime_manifest=args.runtime_manifest,
            expected_runtime_sha256=args.expected_runtime_sha256,
            setup_script=args.setup_script,
            tower_mask=args.tower_mask,
            recoil_config=args.recoil_config,
            lane_tsv=args.lane_tsv,
            expected_estimator_revision=args.expected_estimator_revision,
        )
        status = write_manifest(args.output, payload)
    except (ManifestBuildError, OSError) as exc:
        print(f"PPG12_PAIRED_SOURCE_MANIFEST_ERROR: {exc}", file=sys.stderr)
        return 2
    print(
        f"PPG12_PAIRED_SOURCE_MANIFEST_{status} "
        f"output={args.output} sha256={_sha256(args.output)} lanes=12"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
