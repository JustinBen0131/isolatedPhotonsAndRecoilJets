#!/usr/bin/env python3
"""Resolve a source-complete THE-134 multi-view extraction without submitting.

This controller is deliberately distinct from
``submit_the134_multiview_extraction_smoke.sh``.  The smoke owns one tuple per
source and can never acquire training authority.  This resolver instead proves
the immutable inputs for every tuple in the thirteen-source training matrix and
emits typed row descriptors with ``requested_scope=full`` and
``full_training_authority=0``.

The command has only two actions:

``inventory``
    Print the frozen thirteen-row source matrix.

``preflight``
    Verify immutable bundle and source-authority manifests, then atomically
    write a deterministic plan, row-descriptor JSONL, duplicate fingerprint,
    and receipt.  It does not invoke SSH, Condor, or an existing submitter.

Full authority is still earned only after the planned outputs pass
``prepare_the134_h70_matrix.py --scope full`` with exact source-population
closure.  A preflight receipt is therefore execution authority input, never a
scientific PASS certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import re
import shutil
import sys
import tempfile
from pathlib import Path, PurePosixPath
from typing import Any, Iterable


HERE = Path(__file__).resolve().parent
MATERIALIZER_PATH = HERE / "materialize_the134_immutable_bundle.py"
MATERIALIZER_SPEC = importlib.util.spec_from_file_location(
    "the134_immutable_bundle_materializer", MATERIALIZER_PATH
)
if MATERIALIZER_SPEC is None or MATERIALIZER_SPEC.loader is None:
    raise RuntimeError(
        f"cannot load immutable-bundle materializer: {MATERIALIZER_PATH}"
    )
bundle_materializer = importlib.util.module_from_spec(MATERIALIZER_SPEC)
MATERIALIZER_SPEC.loader.exec_module(bundle_materializer)


BUNDLE_SCHEMA = "THE134_FULL_EXTRACTION_IMMUTABLE_BUNDLE_V1"
MATERIALIZATION_SCHEMA = (
    "THE134_IMMUTABLE_BUILD_BUNDLE_MATERIALIZATION_V1"
)
SOURCE_SCHEMA = "THE134_FULL_EXTRACTION_SOURCE_AUTHORITY_V2"
SOURCE_AUTHORITY_SCOPE = "THE134_FULL_13_ROW_TRAINING_EXTRACTION"
SOURCE_AUTHORITY_ROW_COUNT = 13
SOURCE_AUTHORITY_AUAU_PERIOD = "AUAU_RUN24"
SOURCE_AUTHORITY_AUAU_SI_DI_ROLE = "EMBEDDED"
SOURCE_AUTHORITY_OWNERSHIP_STATE = "source_role_frozen"
SOURCE_AUTHORITY_DIAGNOSTIC_EXCLUSIONS = ("run28_jet40",)
PLAN_SCHEMA = "THE134_FULL_MULTIVIEW_EXTRACTION_PLAN_V2"
ROW_SCHEMA = "THE134_FULL_MULTIVIEW_EXTRACTION_ROW_V2"
RECEIPT_SCHEMA = "THE134_FULL_MULTIVIEW_EXTRACTION_PREFLIGHT_RECEIPT_V2"
DUPLICATE_SCHEMA = "THE134_FULL_MULTIVIEW_EXTRACTION_DUPLICATE_CONTRACT_V1"
EXECUTION_SCHEMA = "THE134_FULL_MULTIVIEW_EXTRACTION_EXECUTION_CONTRACT_V2"
PARTITION_SCHEMA = "THE134_FULL_EXTRACTION_PARTITION_CONTRACT_V3"
ROW_PARTITION_SCHEMA = "THE134_FULL_EXTRACTION_ROW_PARTITION_V1"
CHUNK_SCHEMA = "THE134_FULL_EXTRACTION_CHUNK_V1"
CHUNK_KEYS = frozenset(
    {
        "schema",
        "row_id",
        "full_source_manifest_sha256",
        "group_size",
        "chunk_index",
        "segment",
        "tuple_index_start",
        "tuple_index_end_exclusive",
        "tuple_count",
        "physical_lines",
        "tuple_input_sha256s",
        "tuple_records_sha256",
        "chunk_fingerprint_sha256",
    }
)
ROW_PARTITION_KEYS = frozenset(
    {
        "schema",
        "row_id",
        "full_source_manifest_sha256",
        "group_size",
        "source_tuple_count",
        "expected_chunk_count",
        "expected_job_count",
        "expected_output_pair_count",
        "expected_analysis_output_count",
        "expected_sidecar_output_count",
        "expected_source_occurrence_count",
        "source_occurrences_per_output_pair",
        "tail_tuple_count",
        "first_chunk_sha256",
        "last_chunk_sha256",
        "chunk_fingerprints",
        "chunk_records_sha256",
        "partition_sha256",
    }
)
EXECUTION_CHUNK_KEYS = CHUNK_KEYS | frozenset(
    {"global_chunk_index", "execution_chunk_sha256"}
)
TRAINING_SOURCE_SCHEMA = "THE134_PP_TRAINING_PERIOD_SI_CONTRACT_V1"
CLOSURE_BOUNDARY_SCHEMA = "THE134_FULL_TRAINING_CLOSURE_WITNESS_BOUNDARY_V1"

REQUESTED_SCOPE = "full"
PREFLIGHT_FULL_TRAINING_AUTHORITY = 0
PREFLIGHT_AUTHORITY_STATE = "PREFLIGHT_RESOLVED_NOT_EARNED"
AUTHORITY_EARNING_TOOL = "prepare_the134_h70_matrix.py"
AUTHORITY_EARNING_SCOPE = "full"
MODEL_DOMAIN_GEV = (15.0, 35.0)
CAPTURE_ET_MIN_GEV = 5.0
EXTRACTION_CONE_R = 0.4
LEGACY_TRAINING_TREE_MAX_ENTRIES = 0
EVENT_LIMIT_PER_JOB = 0
GROUP_SIZE = 7
REQUEST_MEMORY_MB = 3000
SHOWER_VIEWS = ("H70", "H0", "G70", "G0", "O70", "O0", "R70")
FROZEN_AUAU_FULL_EXTRACTION_CONTROLS = {
    "RJ_AUAU_BDT_EXTRACT_ONLY": "1",
    "RJ_AUAU_BDT_TRAINING_TREE": "1",
    "RJ_AUAU_BDT_TRAINING_TREE_MAX_ENTRIES": "0",
    "RJ_AUAU_BDT_NPB_DATA_TAGGING": "0",
    "RJ_REQUIRE_EMBEDDED_MINBIAS_CLASSIFIER": "1",
    "RJ_AUAU_BUILD_TOPOCLUSTER_ISOLATION": "0",
    "RJ_AUAU_USE_TOPOCLUSTER_ISOLATION": "0",
    "RJ_SIM_ALLOW_NONE_LISTS": "1",
}
COMMON_SUBMITTER_ADMISSION_KEYS = (
    "RJ_THE134_MULTIVIEW_TRAINING_V1",
    "RJ_THE134_MULTIVIEW_SIDECAR_ONLY_V1",
    "RJ_THE134_EPHEMERAL_ANALYSIS_OUTPUT",
    "RJ_REPLAY_SOURCE_MANIFEST_SHA256",
    "RJ_REPLAY_CONFIG_SHA256",
    "RJ_REPLAY_CODE_SHA256",
    "RJ_PROFILE_LABEL",
)
PP_SUBMITTER_ADMISSION_KEYS = (
    "RJ_PP_PHOTONID_EXTRACT_ONLY",
    "RJ_PP_PHOTONID_TRAINING_TREE",
    "RJ_PP_PHOTONID_TRAINING_TREE_MAX_ENTRIES",
    "RJ_PPG12_PHOTON_YIELD",
    "RJ_PPG12_PHOTON_YIELD_DOUBLE",
)
SIDECAR_ONLY_ARTIFACT_PROFILE = {
    "schema": "THE134_MULTIVIEW_SIDECAR_ONLY_ARTIFACT_PROFILE_V2",
    "artifact_profile": "THE134_MULTIVIEW_SIDECAR_ONLY_V2",
    "analysis_root_role": "EPHEMERAL_CONDOR_SCRATCH_VALIDATED_NOT_RETAINED",
    "analysis_health_record": "WORKER_STDOUT_STRUCTURED_V1",
    "retained_analysis_root_count_per_job": 0,
    "training_sidecar_role": "RJPhotonTrainingViewV1",
    "replay_transaction": "CONSTRUCTED_AND_VALIDATED",
    "replay_serialization": "DISABLED",
    "cache_replay_applicability": "NOT_APPLICABLE",
    "top_level_identity_markers": {
        "rj_replay_schema_sha256": "RJ_REPLAY_SCHEMA_SHA256",
        "rj_replay_semantic_sha256": "RJ_REPLAY_SEMANTIC_SHA256",
        "rj_replay_source_sha256": "RJ_REPLAY_SOURCE_SHA256",
        "rj_replay_model_sha256": "RJ_REPLAY_MODEL_SHA256",
        "rj_replay_config_sha256": "RJ_REPLAY_CONFIG_SHA256",
        "rj_replay_code_sha256": "RJ_REPLAY_CODE_SHA256",
    },
    "full_training_authority": 0,
}
LIST_ROLES = ("calo_cluster", "g4hits", "jets", "global", "mbd_epd")
LIST_FILENAMES = {
    "calo_cluster": "DST_CALO_CLUSTER.matched.list",
    "g4hits": "G4Hits.matched.list",
    "jets": "DST_JETS.matched.list",
    "global": "DST_GLOBAL.matched.list",
    "mbd_epd": "DST_MBD_EPD.matched.list",
}

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

REQUIRED_ARTIFACT_ROLES = (
    "submitter",
    "pp_executor",
    "auau_executor",
    "pp_config",
    "auau_config",
    "pp_library",
    "auau_library",
    "pp_model",
    "auau_model",
    "code_manifest",
    "replay_schema_header",
    "training_schema_header",
    "shower_factorial_header",
    "photon_cluster_header",
    "photon_cluster_builder_header",
    "calo_reco_library",
    "calo_reco_build_receipt",
    "calo_reco_source_manifest",
    "runtime_authority_manifest",
    "release_calo_io",
    "release_clusteriso",
    "release_jetbase",
)

RELEASE_PROVIDER_ROLES = {
    "libcalo_io.so": "release_calo_io",
    "libclusteriso.so": "release_clusteriso",
    "libjetbase.so": "release_jetbase",
}

SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
GIT_COMMIT_RE = re.compile(r"^[0-9a-f]{40}$")
TAG_RE = re.compile(r"^[a-z0-9][a-z0-9._-]{7,127}$")
ROLE_RE = re.compile(r"^[a-z][a-z0-9_]{1,63}$")


class ControllerError(ValueError):
    """Fail-closed contract violation."""


def canonical_json_bytes(payload: Any) -> bytes:
    return json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")


def canonical_sha256(payload: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(payload)).hexdigest()


def build_source_partition(
    row_id: str,
    tuples: list[dict[str, Any]],
    *,
    full_source_manifest_sha256: str,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Build the exact ordered seven-tuple execution partition for one row."""

    expected_indices = list(range(len(tuples)))
    observed_indices = [record.get("tuple_index") for record in tuples]
    if observed_indices != expected_indices:
        raise ControllerError(
            f"{row_id} tuple indices are not the exact ordered range"
        )
    chunks: list[dict[str, Any]] = []
    for chunk_index, start in enumerate(range(0, len(tuples), GROUP_SIZE)):
        stop = min(start + GROUP_SIZE, len(tuples))
        members = tuples[start:stop]
        record = {
            "schema": CHUNK_SCHEMA,
            "row_id": row_id,
            "full_source_manifest_sha256": full_source_manifest_sha256,
            "group_size": GROUP_SIZE,
            "chunk_index": chunk_index,
            "segment": chunk_index + 1,
            "tuple_index_start": start,
            "tuple_index_end_exclusive": stop,
            "tuple_count": len(members),
            "physical_lines": [member["physical_line"] for member in members],
            "tuple_input_sha256s": [
                canonical_sha256(member["inputs"]) for member in members
            ],
            "tuple_records_sha256": canonical_sha256(members),
        }
        record["chunk_fingerprint_sha256"] = canonical_sha256(record)
        chunks.append(record)

    chunk_count = len(chunks)
    partition_payload = {
        "schema": ROW_PARTITION_SCHEMA,
        "row_id": row_id,
        "full_source_manifest_sha256": full_source_manifest_sha256,
        "group_size": GROUP_SIZE,
        "source_tuple_count": len(tuples),
        "expected_chunk_count": chunk_count,
        "expected_job_count": chunk_count,
        "expected_output_pair_count": chunk_count,
        "expected_analysis_output_count": chunk_count,
        "expected_sidecar_output_count": chunk_count,
        "expected_source_occurrence_count": chunk_count,
        "source_occurrences_per_output_pair": 1,
        "tail_tuple_count": chunks[-1]["tuple_count"],
        "first_chunk_sha256": chunks[0]["chunk_fingerprint_sha256"],
        "last_chunk_sha256": chunks[-1]["chunk_fingerprint_sha256"],
        "chunk_fingerprints": [
            record["chunk_fingerprint_sha256"] for record in chunks
        ],
    }
    summary = {
        **partition_payload,
        "chunk_records_sha256": canonical_sha256(chunks),
        "partition_sha256": canonical_sha256(partition_payload),
    }
    validate_source_partition(summary, chunks, tuples)
    return summary, chunks


def validate_source_partition(
    summary: dict[str, Any],
    chunks: list[dict[str, Any]],
    tuples: list[dict[str, Any]],
) -> None:
    """Reject gaps, overlaps, reordering, count drift, and digest drift."""

    row_id = str(summary.get("row_id", ""))
    if set(summary) != ROW_PARTITION_KEYS:
        raise ControllerError("row partition field inventory differs")
    if summary.get("schema") != ROW_PARTITION_SCHEMA or not row_id:
        raise ControllerError("row partition schema or row identity differs")
    source_manifest_sha256 = require_sha256(
        f"{row_id} partition full_source_manifest_sha256",
        summary.get("full_source_manifest_sha256", ""),
    )
    expected_chunk_count = (len(tuples) + GROUP_SIZE - 1) // GROUP_SIZE
    count_fields = {
        "expected_chunk_count": expected_chunk_count,
        "expected_job_count": expected_chunk_count,
        "expected_output_pair_count": expected_chunk_count,
        "expected_analysis_output_count": expected_chunk_count,
        "expected_sidecar_output_count": expected_chunk_count,
        "expected_source_occurrence_count": expected_chunk_count,
    }
    if summary.get("group_size") != GROUP_SIZE:
        raise ControllerError(f"{row_id} partition group_size differs")
    if summary.get("source_tuple_count") != len(tuples):
        raise ControllerError(f"{row_id} partition source-tuple count differs")
    for field, expected in count_fields.items():
        if summary.get(field) != expected:
            raise ControllerError(f"{row_id} partition {field} differs")
    if summary.get("source_occurrences_per_output_pair") != 1:
        raise ControllerError(
            f"{row_id} must have one source occurrence per output pair"
        )
    if len(chunks) != expected_chunk_count:
        raise ControllerError(f"{row_id} partition chunk count differs")

    flattened_tuple_sha256s: list[str] = []
    expected_start = 0
    for chunk_index, chunk in enumerate(chunks):
        if set(chunk) != CHUNK_KEYS:
            raise ControllerError(f"{row_id} chunk field inventory differs")
        if chunk.get("schema") != CHUNK_SCHEMA:
            raise ControllerError(f"{row_id} chunk schema differs")
        if chunk.get("row_id") != row_id:
            raise ControllerError(f"{row_id} chunk row identity differs")
        if chunk.get("full_source_manifest_sha256") != source_manifest_sha256:
            raise ControllerError(f"{row_id} chunk source-manifest binding differs")
        start = chunk.get("tuple_index_start")
        stop = chunk.get("tuple_index_end_exclusive")
        if (
            chunk.get("group_size") != GROUP_SIZE
            or chunk.get("chunk_index") != chunk_index
            or chunk.get("segment") != chunk_index + 1
            or start != expected_start
            or not isinstance(stop, int)
            or stop <= start
            or stop > len(tuples)
        ):
            raise ControllerError(
                f"{row_id} chunk ordering or tuple bounds differ"
            )
        expected_members = tuples[start:stop]
        expected_count = len(expected_members)
        if expected_count > GROUP_SIZE or (
            chunk_index + 1 < expected_chunk_count
            and expected_count != GROUP_SIZE
        ):
            raise ControllerError(f"{row_id} chunk size contract differs")
        if chunk.get("tuple_count") != expected_count:
            raise ControllerError(f"{row_id} chunk tuple count differs")
        expected_physical_lines = [
            member["physical_line"] for member in expected_members
        ]
        expected_tuple_sha256s = [
            canonical_sha256(member["inputs"]) for member in expected_members
        ]
        if chunk.get("physical_lines") != expected_physical_lines:
            raise ControllerError(f"{row_id} chunk physical-line mapping differs")
        if chunk.get("tuple_input_sha256s") != expected_tuple_sha256s:
            raise ControllerError(f"{row_id} chunk tuple identity mapping differs")
        if chunk.get("tuple_records_sha256") != canonical_sha256(
            expected_members
        ):
            raise ControllerError(f"{row_id} chunk tuple-record digest differs")
        fingerprint_payload = {
            key: value
            for key, value in chunk.items()
            if key != "chunk_fingerprint_sha256"
        }
        if chunk.get("chunk_fingerprint_sha256") != canonical_sha256(
            fingerprint_payload
        ):
            raise ControllerError(f"{row_id} chunk fingerprint differs")
        flattened_tuple_sha256s.extend(expected_tuple_sha256s)
        expected_start = stop

    expected_tuple_sha256s = [
        canonical_sha256(record["inputs"]) for record in tuples
    ]
    if (
        expected_start != len(tuples)
        or flattened_tuple_sha256s != expected_tuple_sha256s
    ):
        raise ControllerError(
            f"{row_id} chunk partition is not exact ordered tuple coverage"
        )
    chunk_fingerprints = [
        chunk["chunk_fingerprint_sha256"] for chunk in chunks
    ]
    partition_payload = {
        key: value
        for key, value in summary.items()
        if key not in {"chunk_records_sha256", "partition_sha256"}
    }
    if summary.get("chunk_fingerprints") != chunk_fingerprints:
        raise ControllerError(f"{row_id} chunk fingerprint inventory differs")
    if summary.get("tail_tuple_count") != chunks[-1]["tuple_count"]:
        raise ControllerError(f"{row_id} partition tail tuple count differs")
    if summary.get("first_chunk_sha256") != chunk_fingerprints[0]:
        raise ControllerError(f"{row_id} first chunk digest differs")
    if summary.get("last_chunk_sha256") != chunk_fingerprints[-1]:
        raise ControllerError(f"{row_id} last chunk digest differs")
    if summary.get("chunk_records_sha256") != canonical_sha256(chunks):
        raise ControllerError(f"{row_id} chunk-record inventory digest differs")
    if summary.get("partition_sha256") != canonical_sha256(partition_payload):
        raise ControllerError(f"{row_id} row-partition digest differs")


def training_period_si_contract(pp_period: str) -> dict[str, str]:
    if pp_period not in {"0mrad", "1p5mrad"}:
        raise ControllerError(f"unsupported p+p training period: {pp_period!r}")
    return {
        "schema": TRAINING_SOURCE_SCHEMA,
        "system": "pp",
        "period": pp_period,
        "si_di_role": "SI",
        "row_scope": "ALL_PP_TRAINING_ROWS",
    }


def closure_witness_boundary() -> dict[str, Any]:
    return {
        "schema": CLOSURE_BOUNDARY_SCHEMA,
        "preflight_contains_closure_witness": False,
        "sole_authority_earning_tool": AUTHORITY_EARNING_TOOL,
        "required_argv_suffix": ["--scope", AUTHORITY_EARNING_SCOPE],
        "required_audit_schema": "THE134_FACTORIAL_VIEW_TRAINING_MATRIX_AUDIT_V1",
        "required_audit_status": "PASS",
        "required_source_population_closure_status": "PASS",
    }


def validate_preflight_authority_payload(payload: Any, *, label: str) -> None:
    """Reject any preflight surface that claims authority it cannot earn."""

    def walk(value: Any, path: str) -> Iterable[tuple[str, Any, str]]:
        if isinstance(value, dict):
            for key, nested in value.items():
                nested_path = f"{path}.{key}" if path else str(key)
                yield str(key), nested, nested_path
                yield from walk(nested, nested_path)
        elif isinstance(value, list):
            for index, nested in enumerate(value):
                yield from walk(nested, f"{path}[{index}]")

    for key, value, path in walk(payload, label):
        if key == "requested_scope" and value != REQUESTED_SCOPE:
            raise ControllerError(
                f"{path} must remain requested_scope={REQUESTED_SCOPE}"
            )
        if (
            key == "full_training_authority"
            and value != PREFLIGHT_FULL_TRAINING_AUTHORITY
        ):
            raise ControllerError(
                f"{path} must remain full_training_authority="
                f"{PREFLIGHT_FULL_TRAINING_AUTHORITY} during preflight"
            )
        if key == "authority_state" and value != PREFLIGHT_AUTHORITY_STATE:
            raise ControllerError(
                f"{path} must remain authority_state={PREFLIGHT_AUTHORITY_STATE}"
            )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_sha256(label: str, value: object) -> str:
    text = str(value)
    if not SHA256_RE.fullmatch(text):
        raise ControllerError(
            f"{label} must be a lowercase 64-character SHA-256"
        )
    return text


def require_positive_int(label: str, value: object) -> int:
    if isinstance(value, bool):
        raise ControllerError(f"{label} must be a positive integer")
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise ControllerError(f"{label} must be a positive integer") from exc
    if result <= 0 or str(result) != str(value):
        raise ControllerError(f"{label} must be a canonical positive integer")
    return result


def require_safe_text(label: str, value: object) -> str:
    text = str(value)
    if not text or any(character in text for character in ("\x00", "\n", "\r", "\t")):
        raise ControllerError(f"{label} is empty or contains a control character")
    return text


def require_local_file(label: str, value: object) -> Path:
    text = require_safe_text(label, value)
    path = Path(text).expanduser()
    if not path.is_absolute():
        raise ControllerError(f"{label} must be an absolute path: {path}")
    if not path.is_file():
        raise ControllerError(f"{label} does not exist as a regular file: {path}")
    return path


def require_local_dir(label: str, value: object) -> Path:
    text = require_safe_text(label, value)
    path = Path(text).expanduser()
    if not path.is_absolute():
        raise ControllerError(f"{label} must be an absolute path: {path}")
    if not path.is_dir():
        raise ControllerError(f"{label} does not exist as a directory: {path}")
    return path


def require_remote_root(label: str, value: object) -> str:
    text = require_safe_text(label, value)
    if ";" in text:
        raise ControllerError(f"{label} cannot contain a semicolon")
    path = PurePosixPath(text)
    if not path.is_absolute() or ".." in path.parts or str(path) != text.rstrip("/"):
        raise ControllerError(
            f"{label} must be a normalized absolute POSIX path without a trailing slash"
        )
    if str(path) == "/":
        raise ControllerError(f"{label} cannot be the filesystem root")
    return str(path)


def load_pinned_json(path: Path, expected_sha256: str, label: str) -> dict[str, Any]:
    expected = require_sha256(f"{label} expected hash", expected_sha256)
    if not path.is_file():
        raise ControllerError(f"{label} is missing: {path}")
    observed = sha256_file(path)
    if observed != expected:
        raise ControllerError(
            f"{label} hash drift: expected={expected} observed={observed} path={path}"
        )
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ControllerError(f"{label} is not valid JSON: {path}") from exc
    if not isinstance(payload, dict):
        raise ControllerError(f"{label} must contain one JSON object")
    return payload


def inventory_rows() -> list[dict[str, str]]:
    keys = (
        "row_id",
        "system",
        "lane",
        "dataset",
        "sample",
        "source_role",
        "minimum_bias_gate",
        "photon_id_row_match",
    )
    return [dict(zip(keys, row)) for row in SOURCE_ROWS]


def validate_inventory() -> None:
    rows = inventory_rows()
    row_ids = {row["row_id"] for row in rows}
    samples = {row["sample"] for row in rows}
    if len(rows) != 13 or len(row_ids) != 13 or len(samples) != 13:
        raise ControllerError("THE-134 source matrix is not exactly 13 unique rows")
    if "run28_jet40" in samples:
        raise ControllerError("p+p Jet40 is diagnostic-only and cannot enter training")
    expected_counts = {
        ("pp", "signal"): 3,
        ("pp", "background"): 4,
        ("auau", "signal"): 2,
        ("auau", "background"): 4,
    }
    observed_counts = {
        key: sum(
            row["system"] == key[0] and row["source_role"] == key[1]
            for row in rows
        )
        for key in expected_counts
    }
    if observed_counts != expected_counts:
        raise ControllerError(
            f"THE-134 source-family closure differs: {observed_counts}"
        )
    if any(
        row["minimum_bias_gate"] != "required_pass"
        for row in rows
        if row["system"] == "auau"
    ):
        raise ControllerError("every Au+Au training row must require MinimumBias")


def resolve_release_providers(
    *,
    release_core_lib_dir: Path,
    release_core_lib64_dir: Path,
    artifacts: dict[str, dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    """Resolve the one pinned runtime provider for each bundled companion.

    The immutable bundle retains byte-identical copies as provenance artifacts,
    while execution intentionally uses the versioned CVMFS release provider.
    Search order must match the frozen wrapper exactly: lib64, then lib.
    """

    providers: dict[str, dict[str, Any]] = {}
    search_roots = (release_core_lib64_dir, release_core_lib_dir)
    for family, role in RELEASE_PROVIDER_ROLES.items():
        authority = artifacts.get(role)
        if not isinstance(authority, dict):
            raise ControllerError(
                f"bundle is missing release-provider authority: {role}"
            )
        expected_sha = require_sha256(
            f"bundle artifact {role}", authority.get("sha256", "")
        )
        candidates = [root / family for root in search_roots]
        existing = [candidate for candidate in candidates if candidate.is_file()]
        if not existing:
            raise ControllerError(
                f"pinned release provider is missing: {family}"
            )
        selected = existing[0]
        observed_sha = sha256_file(selected)
        if observed_sha != expected_sha:
            raise ControllerError(
                f"pinned release provider hash differs from immutable bundle: "
                f"family={family} expected={expected_sha} "
                f"observed={observed_sha} path={selected}"
            )
        providers[family] = {
            "path": str(selected),
            "resolved_path": str(selected.resolve(strict=True)),
            "sha256": observed_sha,
            "bundle_artifact_path": str(authority["path"]),
        }
    return providers


def validate_bundle(payload: dict[str, Any]) -> dict[str, Any]:
    if payload.get("schema") != BUNDLE_SCHEMA:
        raise ControllerError(f"bundle schema must be {BUNDLE_SCHEMA}")
    if payload.get("status") != "PASS":
        raise ControllerError("immutable bundle receipt is not PASS")
    public_commit = str(payload.get("public_commit", ""))
    if not GIT_COMMIT_RE.fullmatch(public_commit):
        raise ControllerError("bundle public_commit must be a 40-character Git SHA")
    hashes = {
        name: require_sha256(f"bundle {name}", payload.get(name, ""))
        for name in (
            "code_sha256",
            "replay_schema_sha256",
            "training_schema_sha256",
            "semantic_sha256",
        )
    }
    bundle_identity = require_sha256(
        "bundle identity", payload.get("bundle_identity_sha256", "")
    )
    if payload.get("semantic_fingerprint_sha256") != bundle_identity:
        raise ControllerError(
            "bundle semantic fingerprint differs from bundle identity"
        )
    runtime = payload.get("runtime")
    if not isinstance(runtime, dict):
        raise ControllerError("bundle runtime contract must be an object")
    release = require_safe_text("bundle runtime release", runtime.get("release", ""))
    if release != "ana.560":
        raise ControllerError("bundle runtime release must remain ana.560")
    offline_main = require_remote_root(
        "bundle runtime offline_main", runtime.get("offline_main", "")
    )
    soname = require_safe_text(
        "bundle runtime calo_reco_soname", runtime.get("calo_reco_soname", "")
    )
    if soname != "libcalo_reco.so.0":
        raise ControllerError("bundle runtime must pin libcalo_reco.so.0")
    if require_positive_int(
        "bundle runtime request_memory_mb", runtime.get("request_memory_mb", 0)
    ) != REQUEST_MEMORY_MB:
        raise ControllerError(
            f"bundle runtime request_memory_mb must remain {REQUEST_MEMORY_MB}"
        )
    release_core_lib_dir = require_local_dir(
        "bundle runtime release_core_lib_dir",
        runtime.get("release_core_lib_dir", ""),
    )
    release_core_lib64_dir = require_local_dir(
        "bundle runtime release_core_lib64_dir",
        runtime.get("release_core_lib64_dir", ""),
    )

    raw_artifacts = payload.get("artifacts")
    if not isinstance(raw_artifacts, list):
        raise ControllerError("bundle artifacts must be a list")
    artifacts: dict[str, dict[str, Any]] = {}
    for raw in raw_artifacts:
        if not isinstance(raw, dict):
            raise ControllerError("bundle artifact records must be objects")
        role = str(raw.get("role", ""))
        if not ROLE_RE.fullmatch(role):
            raise ControllerError(f"invalid bundle artifact role: {role!r}")
        if role in artifacts:
            raise ControllerError(f"duplicate bundle artifact role: {role}")
        path = require_local_file(f"bundle artifact {role}", raw.get("path", ""))
        expected_sha = require_sha256(
            f"bundle artifact {role}", raw.get("sha256", "")
        )
        observed_sha = sha256_file(path)
        if observed_sha != expected_sha:
            raise ControllerError(
                f"bundle artifact hash drift: role={role} expected={expected_sha} "
                f"observed={observed_sha} path={path}"
            )
        expected_size = require_positive_int(
            f"bundle artifact {role} size_bytes", raw.get("size_bytes", 0)
        )
        observed_size = path.stat().st_size
        if observed_size != expected_size:
            raise ControllerError(
                f"bundle artifact size drift: role={role} expected={expected_size} "
                f"observed={observed_size} path={path}"
            )
        if role in {"submitter", "pp_executor", "auau_executor"} and not os.access(
            path, os.X_OK
        ):
            raise ControllerError(
                f"bundle artifact must be executable: role={role} path={path}"
            )
        artifacts[role] = {
            "role": role,
            "path": str(path),
            "resolved_path": str(path.resolve(strict=True)),
            "sha256": observed_sha,
            "size_bytes": observed_size,
        }
    missing = sorted(set(REQUIRED_ARTIFACT_ROLES) - set(artifacts))
    if missing:
        raise ControllerError(f"bundle is missing required artifact roles: {missing}")
    release_providers = resolve_release_providers(
        release_core_lib_dir=release_core_lib_dir,
        release_core_lib64_dir=release_core_lib64_dir,
        artifacts=artifacts,
    )
    normalized_artifacts = [artifacts[role] for role in sorted(artifacts)]
    return {
        "schema": BUNDLE_SCHEMA,
        "status": "PASS",
        "public_commit": public_commit,
        "bundle_identity_sha256": bundle_identity,
        **hashes,
        "runtime": {
            "release": release,
            "offline_main": offline_main,
            "calo_reco_soname": soname,
            "request_memory_mb": REQUEST_MEMORY_MB,
            "release_core_lib_dir": str(release_core_lib_dir),
            "release_core_lib64_dir": str(release_core_lib64_dir),
            "release_providers": release_providers,
        },
        "artifacts": normalized_artifacts,
        "artifact_by_role": artifacts,
        "semantic_fingerprint_sha256": canonical_sha256(
            {
                "public_commit": public_commit,
                **hashes,
                "runtime": {
                    "release": release,
                    "offline_main": offline_main,
                    "calo_reco_soname": soname,
                    "request_memory_mb": REQUEST_MEMORY_MB,
                },
                "artifacts": [
                    {
                        "role": artifact["role"],
                        "sha256": artifact["sha256"],
                        "size_bytes": artifact["size_bytes"],
                    }
                    for artifact in normalized_artifacts
                ],
            }
        ),
    }


def validate_materialization_binding(
    payload: dict[str, Any],
    *,
    materialization_path: Path,
    bundle_path: Path,
    bundle_file_sha256: str,
    bundle: dict[str, Any],
) -> dict[str, Any]:
    """Prove the resolver receipt belongs to one rehashed materialized bundle."""

    if payload.get("schema") != MATERIALIZATION_SCHEMA:
        raise ControllerError(
            f"materialization schema must be {MATERIALIZATION_SCHEMA}"
        )
    try:
        validated = bundle_materializer.validate_materialization_payload(
            payload,
            materialization_receipt_path=materialization_path,
            rehash=True,
        )
    except bundle_materializer.MaterializationError as exc:
        raise ControllerError(
            f"immutable bundle materialization is invalid: {exc}"
        ) from exc

    bound_bundle_path = Path(
        str(validated["resolver_bundle_receipt"])
    ).resolve(strict=True)
    if bound_bundle_path != bundle_path:
        raise ControllerError(
            "materialization receipt is not bound to the supplied resolver "
            f"bundle receipt: materialized={bound_bundle_path} "
            f"supplied={bundle_path}"
        )
    bound_bundle_sha256 = require_sha256(
        "materialization resolver bundle receipt",
        validated["resolver_bundle_receipt_sha256"],
    )
    if bound_bundle_sha256 != bundle_file_sha256:
        raise ControllerError(
            "materialization resolver bundle receipt hash differs: "
            f"materialized={bound_bundle_sha256} supplied={bundle_file_sha256}"
        )
    if str(validated["public_commit"]) != bundle["public_commit"]:
        raise ControllerError(
            "materialization public commit differs from resolver bundle"
        )
    bundle_identity = require_sha256(
        "materialization bundle identity",
        validated["bundle_identity_sha256"],
    )
    if bundle_identity != bundle["bundle_identity_sha256"]:
        raise ControllerError(
            "materialization bundle identity differs from resolver bundle"
        )
    return {
        "schema": MATERIALIZATION_SCHEMA,
        "status": "PASS",
        "materialization_receipt": str(materialization_path),
        "resolver_bundle_receipt": str(bound_bundle_path),
        "resolver_bundle_receipt_sha256": bound_bundle_sha256,
        "bundle_identity_sha256": bundle_identity,
        "digest_named_bundle_path": str(
            validated["digest_named_bundle_path"]
        ),
        "artifact_count": int(validated["artifact_count"]),
        "total_artifact_bytes": int(validated["total_artifact_bytes"]),
        "readback": "PASS_SYMLINK_FREE_READONLY_CONTENT_EXACT",
    }


def executable_line(raw: str) -> str | None:
    stripped = raw.strip()
    if not stripped or stripped.startswith("#"):
        return None
    return stripped


def inspect_source_entry(
    raw: dict[str, Any],
    expected_row: dict[str, str],
) -> dict[str, Any]:
    for field in ("row_id", "system", "sample"):
        if str(raw.get(field, "")) != expected_row[field]:
            raise ControllerError(
                f"source authority {expected_row['row_id']} {field} differs: "
                f"{raw.get(field)!r}"
            )
    raw_lists = raw.get("input_lists")
    if not isinstance(raw_lists, list):
        raise ControllerError(
            f"{expected_row['row_id']} input_lists must be a list"
        )
    by_role: dict[str, dict[str, Any]] = {}
    lines_by_role: dict[str, list[str]] = {}
    for record in raw_lists:
        if not isinstance(record, dict):
            raise ControllerError(
                f"{expected_row['row_id']} input-list records must be objects"
            )
        role = str(record.get("role", ""))
        if role not in LIST_ROLES or role in by_role:
            raise ControllerError(
                f"{expected_row['row_id']} has invalid or duplicate list role: {role!r}"
            )
        path = require_local_file(
            f"{expected_row['row_id']} {role} list", record.get("path", "")
        )
        if path.name != LIST_FILENAMES[role]:
            raise ControllerError(
                f"{expected_row['row_id']} {role} list must be named "
                f"{LIST_FILENAMES[role]}"
            )
        expected_sha = require_sha256(
            f"{expected_row['row_id']} {role} list", record.get("sha256", "")
        )
        observed_sha = sha256_file(path)
        if observed_sha != expected_sha:
            raise ControllerError(
                f"{expected_row['row_id']} {role} list hash drift: "
                f"expected={expected_sha} observed={observed_sha}"
            )
        lines = path.read_text(encoding="utf-8").splitlines()
        declared_line_count = require_positive_int(
            f"{expected_row['row_id']} {role} line_count",
            record.get("line_count", 0),
        )
        declared_executable_count = require_positive_int(
            f"{expected_row['row_id']} {role} executable_count",
            record.get("executable_count", 0),
        )
        executable_count = sum(executable_line(line) is not None for line in lines)
        if len(lines) != declared_line_count:
            raise ControllerError(
                f"{expected_row['row_id']} {role} physical line count drift"
            )
        if executable_count != declared_executable_count:
            raise ControllerError(
                f"{expected_row['row_id']} {role} executable line count drift"
            )
        by_role[role] = {
            "role": role,
            "path": str(path),
            "resolved_path": str(path.resolve(strict=True)),
            "sha256": observed_sha,
            "line_count": len(lines),
            "executable_count": executable_count,
        }
        lines_by_role[role] = lines
    if set(by_role) != set(LIST_ROLES):
        raise ControllerError(
            f"{expected_row['row_id']} must bind exactly {list(LIST_ROLES)}"
        )
    list_parents = {Path(record["path"]).parent for record in by_role.values()}
    if len(list_parents) != 1:
        raise ControllerError(
            f"{expected_row['row_id']} five source lists do not share one sample root"
        )
    sample_root = next(iter(list_parents))
    if sample_root.name != expected_row["sample"]:
        raise ControllerError(
            f"{expected_row['row_id']} source-list parent must be named "
            f"{expected_row['sample']}"
        )
    physical_counts = {len(lines_by_role[role]) for role in LIST_ROLES}
    if len(physical_counts) != 1:
        raise ControllerError(
            f"{expected_row['row_id']} five-list physical row counts differ"
        )

    tuples: list[dict[str, Any]] = []
    physical_count = next(iter(physical_counts))
    for line_index in range(physical_count):
        values = {
            role: executable_line(lines_by_role[role][line_index])
            for role in LIST_ROLES
        }
        active = sum(value is not None for value in values.values())
        if active == 0:
            continue
        if active != len(LIST_ROLES):
            raise ControllerError(
                f"{expected_row['row_id']} has a partial five-file tuple at "
                f"physical line {line_index + 1}"
            )
        tuples.append(
            {
                "tuple_index": len(tuples),
                "physical_line": line_index + 1,
                "inputs": values,
            }
        )
    if not tuples:
        raise ControllerError(f"{expected_row['row_id']} has no executable tuples")
    executable_counts = {
        by_role[role]["executable_count"] for role in LIST_ROLES
    }
    if executable_counts != {len(tuples)}:
        raise ControllerError(
            f"{expected_row['row_id']} executable tuple/list counts differ"
        )

    tuple_records_sha = canonical_sha256(tuples)
    tuple_input_sha256s = [
        canonical_sha256(record["inputs"]) for record in tuples
    ]
    normalized_lists = [by_role[role] for role in LIST_ROLES]
    source_semantic_payload = {
        "row_id": expected_row["row_id"],
        "system": expected_row["system"],
        "sample": expected_row["sample"],
        "lists": [
            {
                "role": record["role"],
                "sha256": record["sha256"],
                "line_count": record["line_count"],
                "executable_count": record["executable_count"],
            }
            for record in normalized_lists
        ],
        "tuple_count": len(tuples),
        "tuple_records_sha256": tuple_records_sha,
    }
    full_source_manifest_sha = canonical_sha256(source_semantic_payload)
    if require_positive_int(
        f"{expected_row['row_id']} tuple_count",
        raw.get("tuple_count", 0),
    ) != len(tuples):
        raise ControllerError(
            f"{expected_row['row_id']} source tuple_count differs from tuples"
        )
    if any(
        field in raw
        for field in ("expected_input_count", "expected_occurrence_count")
    ):
        raise ControllerError(
            f"{expected_row['row_id']} V2 source authority contains "
            "ambiguous execution-count aliases"
        )
    if require_sha256(
        f"{expected_row['row_id']} tuple_records_sha256",
        raw.get("tuple_records_sha256", ""),
    ) != tuple_records_sha:
        raise ControllerError(
            f"{expected_row['row_id']} tuple-record authority differs"
        )
    if require_sha256(
        f"{expected_row['row_id']} full_source_manifest_sha256",
        raw.get("full_source_manifest_sha256", ""),
    ) != full_source_manifest_sha:
        raise ControllerError(
            f"{expected_row['row_id']} full-source manifest authority differs"
        )
    first_tuple_sha = canonical_sha256(tuples[0])
    last_tuple_sha = canonical_sha256(tuples[-1])
    if require_sha256(
        f"{expected_row['row_id']} first_tuple_sha256",
        raw.get("first_tuple_sha256", ""),
    ) != first_tuple_sha:
        raise ControllerError(
            f"{expected_row['row_id']} first-tuple authority differs"
        )
    if require_sha256(
        f"{expected_row['row_id']} last_tuple_sha256",
        raw.get("last_tuple_sha256", ""),
    ) != last_tuple_sha:
        raise ControllerError(
            f"{expected_row['row_id']} last-tuple authority differs"
        )
    partition, chunks = build_source_partition(
        expected_row["row_id"],
        tuples,
        full_source_manifest_sha256=full_source_manifest_sha,
    )
    return {
        "row_id": expected_row["row_id"],
        "system": expected_row["system"],
        "sample": expected_row["sample"],
        "input_lists": normalized_lists,
        "tuple_count": len(tuples),
        "tuple_records_sha256": tuple_records_sha,
        "full_source_manifest_sha256": full_source_manifest_sha,
        "first_tuple_sha256": first_tuple_sha,
        "last_tuple_sha256": last_tuple_sha,
        "partition_contract": partition,
        "sim_list_root": str(sample_root.parent),
        "_tuple_input_sha256s": tuple_input_sha256s,
        "_chunks": chunks,
    }


def validate_source_manifest(payload: dict[str, Any]) -> list[dict[str, Any]]:
    if payload.get("schema") != SOURCE_SCHEMA:
        raise ControllerError(f"source-authority schema must be {SOURCE_SCHEMA}")
    if payload.get("status") != "PASS":
        raise ControllerError("source-authority manifest is not PASS")
    raw_rows = payload.get("rows")
    if not isinstance(raw_rows, list):
        raise ControllerError("source-authority rows must be a list")
    indexed: dict[str, dict[str, Any]] = {}
    for raw in raw_rows:
        if not isinstance(raw, dict):
            raise ControllerError("source-authority row records must be objects")
        row_id = str(raw.get("row_id", ""))
        if not row_id or row_id in indexed:
            raise ControllerError(
                f"source-authority row id is missing or duplicated: {row_id!r}"
            )
        indexed[row_id] = raw
    expected_rows = inventory_rows()
    expected_ids = {row["row_id"] for row in expected_rows}
    if set(indexed) != expected_ids:
        raise ControllerError(
            "source-authority row closure differs: "
            f"missing={sorted(expected_ids - set(indexed))} "
            f"extra={sorted(set(indexed) - expected_ids)}"
        )
    observed = [
        inspect_source_entry(indexed[row["row_id"]], row) for row in expected_rows
    ]
    sim_roots = {record["sim_list_root"] for record in observed}
    if len(sim_roots) != 1:
        raise ControllerError(
            f"source-authority rows do not share one RJ_SIM_ROOT_OVERRIDE: {sorted(sim_roots)}"
        )
    all_tuple_identities = [
        identity
        for record in observed
        for identity in record.pop("_tuple_input_sha256s")
    ]
    if len(all_tuple_identities) != len(set(all_tuple_identities)):
        raise ControllerError(
            "duplicate five-file source tuples exist across THE-134 rows"
        )
    return observed


def validate_source_period_authority(
    payload: dict[str, Any], *, pp_period: str
) -> dict[str, str]:
    authority = payload.get("authority")
    if not isinstance(authority, dict):
        raise ControllerError(
            "source-authority manifest must contain an authority object"
        )
    frozen_period = require_safe_text(
        "source-authority pp_period", authority.get("pp_period", "")
    )
    if frozen_period != pp_period:
        raise ControllerError(
            "requested p+p period differs from frozen source authority: "
            f"requested={pp_period} frozen={frozen_period}"
        )
    frozen_si_di_role = require_safe_text(
        "source-authority pp_si_di_role",
        authority.get("pp_si_di_role", ""),
    )
    if frozen_si_di_role != "SI":
        raise ControllerError(
            "source-authority pp_si_di_role must remain SI"
        )
    expected_authority = {
        "scope": SOURCE_AUTHORITY_SCOPE,
        "row_count": SOURCE_AUTHORITY_ROW_COUNT,
        "pp_period": frozen_period,
        "pp_si_di_role": frozen_si_di_role,
        "auau_period": SOURCE_AUTHORITY_AUAU_PERIOD,
        "auau_si_di_role": SOURCE_AUTHORITY_AUAU_SI_DI_ROLE,
        "source_ownership_state": SOURCE_AUTHORITY_OWNERSHIP_STATE,
        "diagnostic_sources_excluded": list(
            SOURCE_AUTHORITY_DIAGNOSTIC_EXCLUSIONS
        ),
        "scientific_completion_granted": False,
    }
    if authority != expected_authority:
        raise ControllerError(
            "source-authority V2 field inventory differs from the canonical "
            "13-row builder contract"
        )
    return {
        "pp_period": frozen_period,
        "pp_si_di_role": frozen_si_di_role,
    }


def artifact_for_system(
    artifacts: dict[str, dict[str, Any]], system: str, kind: str
) -> dict[str, Any]:
    role = f"{system}_{kind}"
    return artifacts[role]


def validate_frozen_auau_runtime_environment(
    environment: dict[str, str],
) -> None:
    for key, expected_value in FROZEN_AUAU_FULL_EXTRACTION_CONTROLS.items():
        actual_value = environment.get(key)
        if actual_value != expected_value:
            raise ControllerError(
                "Au+Au full-extraction runtime control "
                f"{key} must remain {expected_value!r}; got {actual_value!r}"
            )


def validate_sidecar_only_contract(
    environment: dict[str, str],
    artifact_profile: dict[str, Any],
) -> None:
    if environment.get("RJ_THE134_MULTIVIEW_SIDECAR_ONLY_V1") != "1":
        raise ControllerError(
            "RJ_THE134_MULTIVIEW_SIDECAR_ONLY_V1 must remain '1'"
        )
    if environment.get("RJ_THE134_EPHEMERAL_ANALYSIS_OUTPUT") != "1":
        raise ControllerError(
            "RJ_THE134_EPHEMERAL_ANALYSIS_OUTPUT must remain '1'"
        )
    if artifact_profile != SIDECAR_ONLY_ARTIFACT_PROFILE:
        raise ControllerError(
            "THE-134 sidecar-only artifact profile differs from the frozen "
            "validation-only contract"
        )


def runtime_environment(
    row: dict[str, str],
    *,
    tag: str,
    pp_period: str,
    bundle: dict[str, Any],
    source: dict[str, Any],
    sidecar_template: str,
) -> dict[str, str]:
    system = row["system"]
    artifacts = bundle["artifact_by_role"]
    model = artifact_for_system(artifacts, system, "model")
    config = artifact_for_system(artifacts, system, "config")
    environment = {
        "RJ_REPLAY_FOUNDATION_V1": "1",
        "RJ_REPLAY_FOUNDATION_CANARY": "0",
        "RJ_REPLAY_TRACE": "0",
        "RJ_REPLAY_LANE": row["lane"],
        "RJ_REPLAY_DATASET": row["dataset"],
        "RJ_REPLAY_SAMPLE": row["sample"],
        "RJ_REPLAY_PERIOD": pp_period if system == "pp" else "AUAU_RUN24",
        "RJ_REPLAY_SI_DI_ROLE": "SI" if system == "pp" else "EMBEDDED",
        "RJ_REPLAY_OWNERSHIP_STATE": "source_role_frozen",
        "RJ_REPLAY_SOURCE_MANIFEST_SHA256": source[
            "full_source_manifest_sha256"
        ],
        "RJ_REPLAY_SOURCE_SHA256": source["full_source_manifest_sha256"],
        "RJ_REPLAY_MODEL_SHA256": model["sha256"],
        "RJ_REPLAY_MODEL_SCORE_NAME": (
            "tight_bdt_score" if system == "pp" else "auau_tight_bdt_score"
        ),
        "RJ_REPLAY_MODEL_SHOWER_DEFINITION": "H70" if system == "pp" else "H0",
        "RJ_REPLAY_CONFIG_SHA256": config["sha256"],
        "RJ_REPLAY_CODE_SHA256": bundle["code_sha256"],
        "RJ_REPLAY_SCHEMA_SHA256": bundle["replay_schema_sha256"],
        "RJ_REPLAY_SEMANTIC_SHA256": bundle["semantic_sha256"],
        "RJ_REPLAY_PHOTON_CAPTURE_ET_MIN": str(CAPTURE_ET_MIN_GEV),
        "RJ_REPLAY_JET_CONSTITUENT_PT_MIN": str(CAPTURE_ET_MIN_GEV),
        "RJ_THE134_MULTIVIEW_TRAINING_V1": "1",
        "RJ_THE134_MULTIVIEW_TRAINING_FILE": sidecar_template,
        "RJ_THE134_MULTIVIEW_SIDECAR_ONLY_V1": "1",
        "RJ_THE134_EPHEMERAL_ANALYSIS_OUTPUT": "1",
        "RJ_THE134_EXPECTED_SOURCE_ROLE": row["source_role"],
        "RJ_ID_FANOUT_MAX_ROWS": "1",
        "RJ_AUTO_MERGE": "0",
        "RJ_CLEAN_OUTPUT_BASE": "0",
        "RJ_REQUEST_MEMORY": f"{REQUEST_MEMORY_MB}MB",
        "RJ_REQUIRE_NON_TINY_OUTPUT": "1",
        "RJ_MIN_OUTPUT_BYTES": "50000",
        "RJ_FAIL_ON_MISSING_CALO_INPUT": "1",
        "RJ_VALIDATE_SIM_INPUT_PATHS": "1",
        "RJ_PROFILE_JOB": "1",
        "RJ_JOB_HEARTBEAT_SECONDS": "120",
        "RJ_PROFILE_LABEL": f"{tag}_{row['row_id']}",
        "RJ_PHOTON_ID_ROW_MATCH": row["photon_id_row_match"],
    }
    if system == "pp":
        environment.update(
            {
                "RJ_SIM_ALLOW_NONE_LISTS": "0",
                "RJ_PPG12_PHOTON_YIELD": "1",
                "RJ_PPG12_PHOTON_YIELD_DOUBLE": "0",
                "RJ_PPG12_PERIOD": pp_period,
                "RJ_PPG12_PERIOD_USE_LUMI_WEIGHT": "1",
                "RJ_PPG12_PERIOD_STRICT_DI": "0",
                "RJ_PPG12_PERIOD_ALLOW_ALL_SIM": "0",
                "RJ_PPG12_PERIOD_ALLOW_MIX_OVERRIDE": "0",
                "RJ_PPG12_PERIOD_ALLOW_VERTEX_FILE_OVERRIDE": "0",
                "RJ_PP_PHOTONID_EXTRACT_ONLY": "1",
                "RJ_PP_PHOTONID_TRAINING_TREE": "1",
                "RJ_PP_PHOTONID_TRAINING_TREE_MAX_ENTRIES": str(
                    LEGACY_TRAINING_TREE_MAX_ENTRIES
                ),
                "RJ_PP_PHOTONID_SOURCE_ROLE": row["source_role"],
                "RJ_PP_PHOTONID_PPG12_FILTER": "1",
                "RJ_PP_PHOTONID_REQUIRE_PRESELECTION": "0",
            }
        )
    else:
        environment.update(FROZEN_AUAU_FULL_EXTRACTION_CONTROLS)
        validate_frozen_auau_runtime_environment(environment)
    validate_sidecar_only_contract(
        environment, dict(SIDECAR_ONLY_ARTIFACT_PROFILE)
    )
    return dict(sorted(environment.items()))


def serialized_environment(environment: dict[str, str]) -> str:
    for key, value in environment.items():
        if not re.fullmatch(r"[A-Z][A-Z0-9_]*", key):
            raise ControllerError(f"environment key is not canonical: {key!r}")
        if any(character in value for character in (";", "\n", "\r", "\x00")):
            raise ControllerError(
                f"environment value cannot be serialized safely: {key}"
            )
    return ";".join(f"{key}={value}" for key, value in sorted(environment.items()))


def submitter_admission_environment(
    system: str, worker_environment: dict[str, str]
) -> dict[str, str]:
    """Return only the worker bindings the submit shell must validate itself.

    These values remain in ``RJ_SUBMIT_EXTRA_ENV`` for the worker.  The sealed
    submit process also needs an explicit copy so its pre-Condor admission can
    verify the exact execution manifest.  p+p photon-yield controls must not
    leak into Au+Au materialization.
    """

    if system not in {"pp", "auau"}:
        raise ControllerError(
            f"unsupported system for submitter admission: {system}"
        )
    keys = list(COMMON_SUBMITTER_ADMISSION_KEYS)
    if system == "pp":
        keys.extend(
            key
            for key in worker_environment
            if key.startswith("RJ_PPG12_")
            or key.startswith("RJ_PP_PHOTONID_")
        )
        missing_required = [
            key
            for key in PP_SUBMITTER_ADMISSION_KEYS
            if key not in worker_environment
        ]
        if missing_required:
            raise ControllerError(
                "p+p worker environment is missing required submitter bindings: "
                + ",".join(missing_required)
            )
    keys = list(dict.fromkeys(keys))
    missing = [key for key in keys if key not in worker_environment]
    if missing:
        raise ControllerError(
            "worker environment is missing submitter admission bindings: "
            + ",".join(missing)
        )
    return {key: worker_environment[key] for key in keys}


def build_descriptors(
    *,
    tag: str,
    pp_period: str,
    training_source_contract_sha256: str,
    output_root: str,
    evidence_root: str,
    submit_root: str,
    partition_artifact_name: str,
    partition_artifact_sha256: str,
    bundle: dict[str, Any],
    sources: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    source_by_row = {source["row_id"]: source for source in sources}
    artifacts = bundle["artifact_by_role"]
    descriptors: list[dict[str, Any]] = []
    for row in inventory_rows():
        source = source_by_row[row["row_id"]]
        system = row["system"]
        sidecar_template = (
            f"{output_root}/{row['row_id']}/training_views/"
            "$(Cluster).$(Process).root"
        )
        analysis_namespace = f"{output_root}/{row['row_id']}"
        worker_environment = runtime_environment(
            row,
            tag=tag,
            pp_period=pp_period,
            bundle=bundle,
            source=source,
            sidecar_template=sidecar_template,
        )
        artifact_profile = dict(SIDECAR_ONLY_ARTIFACT_PROFILE)
        validate_sidecar_only_contract(worker_environment, artifact_profile)
        materialization_environment = {
            "RJ_DAG_DRYRUN": "1",
            "RJ_CONDOR_SEALED_ENVIRONMENT": "1",
            "RJ_SIM_ROOT_OVERRIDE": source["sim_list_root"],
            "RJ_CONFIG_YAML": artifact_for_system(
                artifacts, system, "config"
            )["path"],
            (
                "RJ_PP_LIBRARY_OVERRIDE"
                if system == "pp"
                else "RJ_AUAU_LIBRARY_OVERRIDE"
            ): artifact_for_system(artifacts, system, "library")["path"],
            "RJ_PHOTON_CLUSTER_BUILDER_HEADER_OVERRIDE": artifacts[
                "photon_cluster_builder_header"
            ]["path"],
            "RJ_CALO_RECO_LIBRARY_OVERRIDE": artifacts["calo_reco_library"][
                "path"
            ],
            "RJ_PHOTON_CLUSTER_BUILDER_LIBRARY_OVERRIDE": "",
            "RJ_PINNED_CALO_RECO_RELEASE_COMPANIONS": "1",
            "RJ_PINNED_CALO_RECO_SONAME": bundle["runtime"][
                "calo_reco_soname"
            ],
            "RJ_PINNED_RELEASE_NAME": bundle["runtime"]["release"],
            "RJ_PINNED_OFFLINE_MAIN": bundle["runtime"]["offline_main"],
            "RJ_PINNED_RELEASE_CALO_IO_PATH": bundle["runtime"][
                "release_providers"
            ]["libcalo_io.so"]["path"],
            "RJ_PINNED_RELEASE_CALO_IO_SHA256": bundle["runtime"][
                "release_providers"
            ]["libcalo_io.so"]["sha256"],
            "RJ_PINNED_RELEASE_CLUSTERISO_PATH": bundle["runtime"][
                "release_providers"
            ]["libclusteriso.so"]["path"],
            "RJ_PINNED_RELEASE_CLUSTERISO_SHA256": bundle["runtime"][
                "release_providers"
            ]["libclusteriso.so"]["sha256"],
            "RJ_PINNED_RELEASE_JETBASE_PATH": bundle["runtime"][
                "release_providers"
            ]["libjetbase.so"]["path"],
            "RJ_PINNED_RELEASE_JETBASE_SHA256": bundle["runtime"][
                "release_providers"
            ]["libjetbase.so"]["sha256"],
            "RJ_RELEASE_CORE_LIB_DIR": bundle["runtime"][
                "release_core_lib_dir"
            ],
            "RJ_RELEASE_CORE_LIB64_DIR": bundle["runtime"][
                "release_core_lib64_dir"
            ],
            "RJ_AUTO_MERGE": "0",
            "RJ_STAGE_EMAIL_MODE": "none",
            "RJ_CLEAN_OUTPUT_BASE": "0",
            "RJ_REQUEST_MEMORY": f"{REQUEST_MEMORY_MB}MB",
            "RJ_REQUIRE_NON_TINY_OUTPUT": "1",
            "RJ_MIN_OUTPUT_BYTES": "50000",
            "RJ_FAIL_ON_MISSING_CALO_INPUT": "1",
            "RJ_VALIDATE_SIM_INPUT_PATHS": "1",
            "RJ_PROFILE_JOB": "1",
            "RJ_JOB_HEARTBEAT_SECONDS": "120",
            "RJ_DEST_BASE_OVERRIDE": analysis_namespace,
            "RJ_SUBMISSION_NAMESPACE": row["row_id"],
            "RJ_CONDOR_SUB_DIR": f"{submit_root}/{row['row_id']}",
            "RJ_PHOTON_ID_ROW_MATCH": row["photon_id_row_match"],
            "RJ_ID_FANOUT_MAX_ROWS": "1",
            "RJ_SIM_ALLOW_NONE_LISTS": "0" if system == "pp" else "1",
            "RJ_SUBMIT_EXTRA_ENV": serialized_environment(worker_environment),
        }
        materialization_environment.update(
            submitter_admission_environment(system, worker_environment)
        )
        descriptor = {
            "schema": ROW_SCHEMA,
            **row,
            "source_period": pp_period if system == "pp" else "AUAU_RUN24",
            "source_si_di_role": "SI" if system == "pp" else "EMBEDDED",
            "training_period_si_contract_sha256": (
                training_source_contract_sha256 if system == "pp" else None
            ),
            "requested_scope": REQUESTED_SCOPE,
            "full_training_authority": PREFLIGHT_FULL_TRAINING_AUTHORITY,
            "authority_state": PREFLIGHT_AUTHORITY_STATE,
            "input_contract": {
                "group_size": GROUP_SIZE,
                "event_limit_per_job": EVENT_LIMIT_PER_JOB,
                "source_tuple_count": source["tuple_count"],
                "expected_chunk_count": source["partition_contract"][
                    "expected_chunk_count"
                ],
                "expected_job_count": source["partition_contract"][
                    "expected_job_count"
                ],
                "expected_output_pair_count": source["partition_contract"][
                    "expected_output_pair_count"
                ],
                "expected_analysis_output_count": source[
                    "partition_contract"
                ]["expected_analysis_output_count"],
                "expected_sidecar_output_count": source[
                    "partition_contract"
                ]["expected_sidecar_output_count"],
                "expected_source_occurrence_count": source[
                    "partition_contract"
                ]["expected_source_occurrence_count"],
                "source_occurrences_per_output_pair": source[
                    "partition_contract"
                ]["source_occurrences_per_output_pair"],
                "row_partition_sha256": source["partition_contract"][
                    "partition_sha256"
                ],
                "chunk_records_sha256": source["partition_contract"][
                    "chunk_records_sha256"
                ],
                "partition_artifact_name": partition_artifact_name,
                "partition_artifact_sha256": partition_artifact_sha256,
                "full_source_manifest_sha256": source[
                    "full_source_manifest_sha256"
                ],
                "tuple_records_sha256": source["tuple_records_sha256"],
                "first_tuple_sha256": source["first_tuple_sha256"],
                "last_tuple_sha256": source["last_tuple_sha256"],
                "input_lists": source["input_lists"],
            },
            "science_contract": {
                "nominal_model_domain_gev": list(MODEL_DOMAIN_GEV),
                "loose_capture_et_min_gev": CAPTURE_ET_MIN_GEV,
                "extraction_cone_r": EXTRACTION_CONE_R,
                "legacy_training_tree_max_entries": (
                    LEGACY_TRAINING_TREE_MAX_ENTRIES
                ),
                "shower_views": list(SHOWER_VIEWS),
            },
            "bundle_contract": {
                "public_commit": bundle["public_commit"],
                "code_sha256": bundle["code_sha256"],
                "replay_schema_sha256": bundle["replay_schema_sha256"],
                "training_schema_sha256": bundle["training_schema_sha256"],
                "semantic_sha256": bundle["semantic_sha256"],
                "library": artifact_for_system(
                    artifacts, system, "library"
                ),
                "model": artifact_for_system(artifacts, system, "model"),
                "config": artifact_for_system(artifacts, system, "config"),
                "submitter": artifacts["submitter"],
                "executor": artifacts[
                    "pp_executor" if system == "pp" else "auau_executor"
                ],
            },
            "execution_contract": {
                "schema": EXECUTION_SCHEMA,
                "state": "PREFLIGHT_ONLY_NO_CONDOR_MUTATION",
                "existing_submitter_argv": [
                    artifacts["submitter"]["path"],
                    row["dataset"],
                    "condorDoAll",
                    "groupSize",
                    str(GROUP_SIZE),
                    f"SAMPLE={row['sample']}",
                ],
                "materialization_requires": "RJ_DAG_DRYRUN=1",
                "submission_requires": (
                    "separate explicit duplicate-guarded campaign authorization"
                ),
                "required_execution_inputs": [
                    "RJ_CODEX_CHAT_NAME",
                    "RJ_CODEX_THREAD_ID",
                ],
                "partition_artifact_name": partition_artifact_name,
                "partition_artifact_sha256": partition_artifact_sha256,
                "row_partition_sha256": source["partition_contract"][
                    "partition_sha256"
                ],
                "analysis_output_namespace": analysis_namespace,
                "multiview_sidecar_template": sidecar_template,
                "submit_namespace": f"{submit_root}/{row['row_id']}",
                "evidence_namespace": f"{evidence_root}/{row['row_id']}",
                "artifact_profile": artifact_profile,
                "materialization_environment": dict(
                    sorted(materialization_environment.items())
                ),
                "worker_environment": worker_environment,
                "forbidden_ambient_environment": [
                    "RJ_PPG12_CROSSING_PERIOD",
                    "RJ_PPG12_PHOTON_YIELD_MIX_WEIGHT",
                    "RJ_PP_VERTEX_REWEIGHT_FILE",
                    "RJ_PPG12_PHOTON_YIELD_TRUTH_VERTEX",
                    "RJ_PPG12_PHOTON_YIELD_BUILDER_TRUTH_VERTEX",
                    "RJ_PPG12_PHOTON_YIELD_RECO_TRUTH_VERTEX",
                    "RJ_THE134_TRUTH_LABEL_DIAGNOSTIC_V1",
                    "RJ_THE134_TRUTH_LABEL_COUNTER_DIAGNOSTIC_V1",
                    "RJ_THE134_TRUTH_LABEL_PREWEIGHT_DIAGNOSTIC_V1",
                ],
            },
        }
        descriptor["row_fingerprint_sha256"] = canonical_sha256(descriptor)
        descriptors.append(descriptor)
    return descriptors


def duplicate_contract(
    *,
    training_source_contract_sha256: str,
    bundle: dict[str, Any],
    sources: list[dict[str, Any]],
) -> dict[str, Any]:
    return {
        "schema": DUPLICATE_SCHEMA,
        "purpose": "THE134_SOURCE_COMPLETE_FACTORIAL_VIEW_TRAINING_EXTRACTION",
        "requested_scope": REQUESTED_SCOPE,
        "full_training_authority": PREFLIGHT_FULL_TRAINING_AUTHORITY,
        "authority_state": PREFLIGHT_AUTHORITY_STATE,
        "training_period_si_contract_sha256": training_source_contract_sha256,
        "science_contract": {
            "model_domain_gev": list(MODEL_DOMAIN_GEV),
            "capture_et_min_gev": CAPTURE_ET_MIN_GEV,
            "extraction_cone_r": EXTRACTION_CONE_R,
            "shower_views": list(SHOWER_VIEWS),
        },
        "bundle_semantic_fingerprint_sha256": bundle[
            "semantic_fingerprint_sha256"
        ],
        "sources": [
            {
                "row_id": source["row_id"],
                "system": source["system"],
                "sample": source["sample"],
                "tuple_count": source["tuple_count"],
                "tuple_records_sha256": source["tuple_records_sha256"],
                "full_source_manifest_sha256": source[
                    "full_source_manifest_sha256"
                ],
            }
            for source in sources
        ],
    }


def ordered_partition_chunks(
    sources: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    expected_row_ids = [row["row_id"] for row in inventory_rows()]
    observed_row_ids = [source["row_id"] for source in sources]
    if observed_row_ids != expected_row_ids:
        raise ControllerError("partition source row order differs")
    chunks: list[dict[str, Any]] = []
    global_chunk_index = 0
    for source in sources:
        row_chunks = source.get("_chunks")
        if not isinstance(row_chunks, list):
            raise ControllerError(
                f"{source['row_id']} partition chunks are unavailable"
            )
        for chunk in row_chunks:
            record = {
                **chunk,
                "global_chunk_index": global_chunk_index,
            }
            record["execution_chunk_sha256"] = canonical_sha256(record)
            chunks.append(record)
            global_chunk_index += 1
    validate_ordered_partition_chunks(chunks, sources)
    return chunks


def validate_ordered_partition_chunks(
    chunks: list[dict[str, Any]],
    sources: list[dict[str, Any]],
) -> None:
    expected: list[tuple[str, dict[str, Any]]] = [
        (source["row_id"], chunk)
        for source in sources
        for chunk in source["_chunks"]
    ]
    if len(chunks) != len(expected):
        raise ControllerError("global partition chunk count differs")
    for global_index, (observed, (row_id, source_chunk)) in enumerate(
        zip(chunks, expected)
    ):
        if set(observed) != EXECUTION_CHUNK_KEYS:
            raise ControllerError(
                f"{row_id} global partition field inventory differs"
            )
        if observed.get("global_chunk_index") != global_index:
            raise ControllerError("global partition chunk index differs")
        for field, value in source_chunk.items():
            if observed.get(field) != value:
                raise ControllerError(
                    f"{row_id} global partition chunk payload differs"
                )
        fingerprint_payload = {
            key: value
            for key, value in observed.items()
            if key != "execution_chunk_sha256"
        }
        if observed.get("execution_chunk_sha256") != canonical_sha256(
            fingerprint_payload
        ):
            raise ControllerError(
                f"{row_id} global execution-chunk digest differs"
            )


def aggregate_partition_contract(
    sources: list[dict[str, Any]],
    *,
    partition_artifact_name: str,
    partition_artifact_sha256: str,
) -> dict[str, Any]:
    rows = [source["partition_contract"] for source in sources]
    source_tuple_count = sum(row["source_tuple_count"] for row in rows)
    expected_chunk_count = sum(row["expected_chunk_count"] for row in rows)
    expected_job_count = sum(row["expected_job_count"] for row in rows)
    expected_output_pair_count = sum(
        row["expected_output_pair_count"] for row in rows
    )
    expected_analysis_output_count = sum(
        row["expected_analysis_output_count"] for row in rows
    )
    expected_sidecar_output_count = sum(
        row["expected_sidecar_output_count"] for row in rows
    )
    expected_source_occurrence_count = sum(
        row["expected_source_occurrence_count"] for row in rows
    )
    equal_execution_counts = {
        expected_chunk_count,
        expected_job_count,
        expected_output_pair_count,
        expected_analysis_output_count,
        expected_sidecar_output_count,
        expected_source_occurrence_count,
    }
    if len(equal_execution_counts) != 1:
        raise ControllerError(
            "aggregate chunk/job/output/source-occurrence counts differ"
        )
    partition_sha256 = require_sha256(
        "partition artifact SHA-256", partition_artifact_sha256
    )
    return {
        "schema": PARTITION_SCHEMA,
        "group_size": GROUP_SIZE,
        "source_tuple_count": source_tuple_count,
        "expected_chunk_count": expected_chunk_count,
        "expected_job_count": expected_job_count,
        "expected_output_pair_count": expected_output_pair_count,
        "expected_analysis_output_count": expected_analysis_output_count,
        "expected_sidecar_output_count": expected_sidecar_output_count,
        "expected_physical_root_artifact_count": (
            expected_analysis_output_count + expected_sidecar_output_count
        ),
        "expected_retained_analysis_output_count": 0,
        "expected_durable_root_artifact_count": (
            expected_sidecar_output_count
        ),
        "expected_source_occurrence_count": (
            expected_source_occurrence_count
        ),
        "source_occurrences_per_output_pair": 1,
        "row_partition_records_sha256": canonical_sha256(rows),
        "row_partition_sha256s": [
            row["partition_sha256"] for row in rows
        ],
        "partition_artifact": {
            "name": partition_artifact_name,
            "sha256": partition_sha256,
            "record_count": expected_chunk_count,
        },
        "basis": (
            "ordered_disjoint_seven_tuple_partition_one_execution_per_chunk"
        ),
        "capacity_canary_required_before_submission": True,
        "capacity_authority_earned": False,
    }


def write_json(path: Path, payload: Any) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )


def write_preflight_outputs(
    *,
    out_dir: Path,
    materialization_path: Path,
    materialization_file_sha256: str,
    materialization: dict[str, Any],
    bundle_path: Path,
    bundle_file_sha256: str,
    source_path: Path,
    source_file_sha256: str,
    tag: str,
    pp_period: str,
    output_root: str,
    evidence_root: str,
    submit_root: str,
    bundle: dict[str, Any],
    sources: list[dict[str, Any]],
) -> dict[str, Any]:
    if out_dir.exists():
        raise ControllerError(f"preflight output directory already exists: {out_dir}")
    out_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(prefix=f".{out_dir.name}.tmp.", dir=str(out_dir.parent))
    )
    try:
        training_source_contract = training_period_si_contract(pp_period)
        training_source_contract_sha256 = canonical_sha256(
            training_source_contract
        )
        partition_path = temporary / "the134_full_extraction_partition.jsonl"
        partition_chunks = ordered_partition_chunks(sources)
        with partition_path.open("w", encoding="utf-8") as stream:
            for chunk in partition_chunks:
                stream.write(
                    json.dumps(
                        chunk,
                        sort_keys=True,
                        separators=(",", ":"),
                        ensure_ascii=True,
                    )
                    + "\n"
                )
        partition_artifact_sha256 = sha256_file(partition_path)
        execution_partition = aggregate_partition_contract(
            sources,
            partition_artifact_name=partition_path.name,
            partition_artifact_sha256=partition_artifact_sha256,
        )
        descriptors = build_descriptors(
            tag=tag,
            pp_period=pp_period,
            training_source_contract_sha256=training_source_contract_sha256,
            output_root=output_root,
            evidence_root=evidence_root,
            submit_root=submit_root,
            partition_artifact_name=partition_path.name,
            partition_artifact_sha256=partition_artifact_sha256,
            bundle=bundle,
            sources=sources,
        )
        duplicate = duplicate_contract(
            training_source_contract_sha256=training_source_contract_sha256,
            bundle=bundle,
            sources=sources,
        )
        duplicate_fingerprint = canonical_sha256(duplicate)
        execution_fingerprint = canonical_sha256(
            {
                "schema": EXECUTION_SCHEMA,
                "tag": tag,
                "output_root": output_root,
                "evidence_root": evidence_root,
                "submit_root": submit_root,
                "materialization_receipt_sha256": (
                    materialization_file_sha256
                ),
                "bundle_manifest_sha256": bundle_file_sha256,
                "source_manifest_sha256": source_file_sha256,
                "duplicate_fingerprint_sha256": duplicate_fingerprint,
                "partition_artifact_sha256": partition_artifact_sha256,
                "execution_partition_sha256": canonical_sha256(
                    execution_partition
                ),
                "row_fingerprints": [
                    descriptor["row_fingerprint_sha256"]
                    for descriptor in descriptors
                ],
            }
        )
        plan = {
            "schema": PLAN_SCHEMA,
            "status": "PREFLIGHT_PASS",
            "execution_state": "PREFLIGHT_ONLY_NO_CONDOR_MUTATION",
            "submission_performed": False,
            "campaign": {
                "tag": tag,
                "output_root": output_root,
                "evidence_root": evidence_root,
                "submit_root": submit_root,
            },
            "training_period_si_contract": training_source_contract,
            "training_period_si_contract_sha256": (
                training_source_contract_sha256
            ),
            "authority": {
                "requested_scope": REQUESTED_SCOPE,
                "full_training_authority": PREFLIGHT_FULL_TRAINING_AUTHORITY,
                "authority_state": PREFLIGHT_AUTHORITY_STATE,
            },
            "artifact_profile": dict(SIDECAR_ONLY_ARTIFACT_PROFILE),
            "closure_witness_boundary": closure_witness_boundary(),
            "source_family_closure": {
                "row_count": len(descriptors),
                "pp_signal": 3,
                "pp_background": 4,
                "auau_signal": 2,
                "auau_background": 4,
                "pp_jet40_excluded_as_diagnostic_only": True,
            },
            "execution_partition": execution_partition,
            "input_manifests": {
                "materialization": {
                    "path": str(materialization_path),
                    "sha256": materialization_file_sha256,
                    "bundle_identity_sha256": materialization[
                        "bundle_identity_sha256"
                    ],
                    "digest_named_bundle_path": materialization[
                        "digest_named_bundle_path"
                    ],
                    "readback": materialization["readback"],
                },
                "bundle": {
                    "path": str(bundle_path),
                    "sha256": bundle_file_sha256,
                    "semantic_fingerprint_sha256": bundle[
                        "semantic_fingerprint_sha256"
                    ],
                },
                "sources": {
                    "path": str(source_path),
                    "sha256": source_file_sha256,
                },
            },
            "duplicate_contract": duplicate,
            "duplicate_fingerprint_sha256": duplicate_fingerprint,
            "execution_fingerprint_sha256": execution_fingerprint,
            "rows": descriptors,
        }
        plan_path = temporary / "the134_full_extraction_plan.json"
        rows_path = temporary / "the134_full_extraction_rows.jsonl"
        duplicate_path = (
            temporary / "the134_full_extraction_duplicate_fingerprint.sha256"
        )
        write_json(plan_path, plan)
        with rows_path.open("w", encoding="utf-8") as stream:
            for descriptor in descriptors:
                stream.write(
                    json.dumps(
                        descriptor,
                        sort_keys=True,
                        separators=(",", ":"),
                        ensure_ascii=True,
                    )
                    + "\n"
                )
        duplicate_path.write_text(duplicate_fingerprint + "\n", encoding="ascii")
        receipt = {
            "schema": RECEIPT_SCHEMA,
            "status": "PASS",
            "submission_performed": False,
            "row_count": len(descriptors),
            "requested_scope": REQUESTED_SCOPE,
            "full_training_authority": PREFLIGHT_FULL_TRAINING_AUTHORITY,
            "authority_state": PREFLIGHT_AUTHORITY_STATE,
            "artifact_profile_sha256": canonical_sha256(
                plan["artifact_profile"]
            ),
            "training_period_si_contract_sha256": (
                training_source_contract_sha256
            ),
            "closure_witness_boundary_sha256": canonical_sha256(
                plan["closure_witness_boundary"]
            ),
            "execution_partition_sha256": canonical_sha256(
                plan["execution_partition"]
            ),
            "bundle_manifest_sha256": bundle_file_sha256,
            "materialization_receipt_sha256": (
                materialization_file_sha256
            ),
            "source_manifest_sha256": source_file_sha256,
            "duplicate_fingerprint_sha256": duplicate_fingerprint,
            "execution_fingerprint_sha256": execution_fingerprint,
            "artifacts": {
                "plan": {
                    "name": plan_path.name,
                    "sha256": sha256_file(plan_path),
                },
                "rows": {
                    "name": rows_path.name,
                    "sha256": sha256_file(rows_path),
                },
                "duplicate_fingerprint": {
                    "name": duplicate_path.name,
                    "sha256": sha256_file(duplicate_path),
                },
                "partition": {
                    "name": partition_path.name,
                    "sha256": partition_artifact_sha256,
                    "record_count": len(partition_chunks),
                },
            },
        }
        validate_preflight_authority_payload(plan, label="plan")
        validate_preflight_authority_payload(receipt, label="receipt")
        receipt_path = temporary / "preflight_receipt.json"
        write_json(receipt_path, receipt)
        os.replace(temporary, out_dir)
        result = {
            "status": "PASS",
            "row_count": len(descriptors),
            "requested_scope": REQUESTED_SCOPE,
            "full_training_authority": PREFLIGHT_FULL_TRAINING_AUTHORITY,
            "authority_state": PREFLIGHT_AUTHORITY_STATE,
            "submission_performed": False,
            "out_dir": str(out_dir),
            "materialization_receipt_sha256": (
                materialization_file_sha256
            ),
            "execution_partition_sha256": canonical_sha256(
                plan["execution_partition"]
            ),
            "partition_artifact_sha256": partition_artifact_sha256,
            "duplicate_fingerprint_sha256": duplicate_fingerprint,
            "execution_fingerprint_sha256": execution_fingerprint,
        }
        validate_preflight_authority_payload(result, label="result")
        return result
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="action", required=True)
    inventory = subparsers.add_parser(
        "inventory", help="print the frozen thirteen-source matrix"
    )
    inventory.add_argument("--format", choices=("json", "tsv"), default="json")
    preflight = subparsers.add_parser(
        "preflight", help="verify inputs and emit a non-submitting full plan"
    )
    preflight.add_argument("--bundle-manifest", type=Path, required=True)
    preflight.add_argument("--bundle-sha256", required=True)
    preflight.add_argument(
        "--materialization-receipt", type=Path, required=True
    )
    preflight.add_argument("--materialization-receipt-sha256", required=True)
    preflight.add_argument("--source-manifest", type=Path, required=True)
    preflight.add_argument("--source-manifest-sha256", required=True)
    preflight.add_argument("--tag", required=True)
    preflight.add_argument("--pp-period", choices=("0mrad", "1p5mrad"), required=True)
    preflight.add_argument("--output-root", required=True)
    preflight.add_argument("--evidence-root", required=True)
    preflight.add_argument("--submit-root", required=True)
    preflight.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args(argv)


def run_inventory(output_format: str) -> int:
    validate_inventory()
    rows = inventory_rows()
    if output_format == "json":
        print(
            json.dumps(
                {
                    "schema": "THE134_FULL_EXTRACTION_SOURCE_INVENTORY_V1",
                    "row_count": len(rows),
                    "rows": rows,
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        fields = list(rows[0])
        print("\t".join(fields))
        for row in rows:
            print("\t".join(row[field] for field in fields))
    return 0


def run_preflight(args: argparse.Namespace) -> int:
    validate_inventory()
    if not TAG_RE.fullmatch(args.tag):
        raise ControllerError(
            "tag must be 8-128 lowercase alphanumeric/period/underscore/hyphen characters"
        )
    output_root = require_remote_root("output root", args.output_root)
    evidence_root = require_remote_root("evidence root", args.evidence_root)
    submit_root = require_remote_root("submit root", args.submit_root)
    if len({output_root, evidence_root, submit_root}) != 3:
        raise ControllerError("output, evidence, and submit roots must be distinct")
    for label, root in (
        ("output root", output_root),
        ("evidence root", evidence_root),
        ("submit root", submit_root),
    ):
        if PurePosixPath(root).name != args.tag:
            raise ControllerError(
                f"{label} basename must equal the exact campaign tag: {args.tag}"
            )
    bundle_path = args.bundle_manifest.resolve(strict=True)
    materialization_path = args.materialization_receipt.absolute()
    if not materialization_path.is_file():
        raise ControllerError(
            "immutable bundle materialization receipt is missing: "
            f"{materialization_path}"
        )
    source_path = args.source_manifest.resolve(strict=True)
    bundle_payload = load_pinned_json(
        bundle_path, args.bundle_sha256, "immutable bundle manifest"
    )
    source_payload = load_pinned_json(
        source_path, args.source_manifest_sha256, "source-authority manifest"
    )
    bundle = validate_bundle(bundle_payload)
    materialization_payload = load_pinned_json(
        materialization_path,
        args.materialization_receipt_sha256,
        "immutable bundle materialization receipt",
    )
    materialization = validate_materialization_binding(
        materialization_payload,
        materialization_path=materialization_path,
        bundle_path=bundle_path,
        bundle_file_sha256=args.bundle_sha256,
        bundle=bundle,
    )
    validate_source_period_authority(
        source_payload, pp_period=args.pp_period
    )
    sources = validate_source_manifest(source_payload)
    result = write_preflight_outputs(
        out_dir=args.out_dir,
        materialization_path=materialization_path,
        materialization_file_sha256=args.materialization_receipt_sha256,
        materialization=materialization,
        bundle_path=bundle_path,
        bundle_file_sha256=args.bundle_sha256,
        source_path=source_path,
        source_file_sha256=args.source_manifest_sha256,
        tag=args.tag,
        pp_period=args.pp_period,
        output_root=output_root,
        evidence_root=evidence_root,
        submit_root=submit_root,
        bundle=bundle,
        sources=sources,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


def main(argv: Iterable[str] | None = None) -> int:
    try:
        args = parse_args(argv)
        if args.action == "inventory":
            return run_inventory(args.format)
        if args.action == "preflight":
            return run_preflight(args)
        raise ControllerError(f"unsupported action: {args.action}")
    except (ControllerError, FileNotFoundError) as exc:
        print(f"[THE134-FULL-EXTRACT][ERROR] {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
