#!/usr/bin/env python3
"""Fail-closed checks for the one-lane PPG12/RecoilJets paired oracle.

The supported contract is intentionally narrow: Photon5, 1.5 mrad, SI,
five source rows, the exact five-value historical PHRandomSeed sequence, and
the G4Hits + truth-jet source graph used by the live preserved executable.
This utility never submits work or changes a current-artifact pointer.
Estimator toy seed 42 is a separate downstream contract.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable


SCHEMA_VERSION = 3
VALID_SAMPLES = {"Photon5", "Photon10", "Photon20"}
VALID_PERIODS = {"0mrad", "1p5mrad"}
VALID_INTERACTIONS = {"SI", "DI"}
LEGACY_PINNED_LANE = {
    "sample": "Photon5",
    "period": "1p5mrad",
    "interaction": "SI",
    "rows": 5,
}
HISTORICAL_PH_SEED_SEQUENCE = [
    2991264730,
    4256268992,
    2394322166,
    874466025,
    2240380304,
]
HISTORICAL_PEDESTAL_SEQUENCE = 534
HISTORICAL_PEDESTAL_FILE = "pedestal-54256-00534.root"
HISTORICAL_RNG_SOURCE = (
    "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/"
    "run28/photon5/condorout/OutDir0/test.out"
)
EXPECTED_RUNTIME_PROFILE = "new.17"
EXPECTED_OFFLINE_MAIN = (
    "/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/"
    "release/release_new/new.17"
)
EXPECTED_PPG_SOURCE_REVISION = "1c0ff86bf0ebabfba63a1abc4512cbe59fe48e31"
EXPECTED_ESTIMATOR_REVISION = "29f8223bd9b36dffab07961b597afa94185bbdf1"
EXPECTED_RECOEFF_SOURCE_HASH = "e9b25fdb6dd8a6bfbbad029cb90aaddc9489fdf2846c630ea63c8c41ac771eee"
EXPECTED_RECOEFF_CONFIG_HASH = "42b7be1628843d5b7607ab988ffb58c6d019d8ade01b3c4d528611498db95732"
EXPECTED_HASHES = {
    "ppg_macro": "95f12ffa8e283dec5add217200b9c36a8c2084459daa4efb59d3fe2a2c21ab3c",
    "g4_full_list": "47978eac14253516d5dbcbfb9381372c49c80879eddae420b4e65de087d7bbd1",
    "truthjet_full_list": "dacbd05167f31109786a7eb38cea89bd33c28432fb97ea2af3bddedbf627ea3e",
    "apply_bdt": "bd6e7c5bc9858ddad9bc835552d818c00290bb7de3f5036f44d8bf4804734366",
    "apply_config": "b8d1bc359a647cc913f213777fc42958b532b30c37a63bb318680130eb6e321b",
    "base_e_model": "7e2d5ed9d1216ab30b92a137482a5bfa111d78e86d4f5a459179c38bc7993997",
    "base_v3e_model": "7679e634260402fb3815b2733767182690eec7587f9e09bffc307a05d00d59df",
    "npb_model": "d6086dadac534013cda15cdfb69c1683776d3456d9e439903589653e8ac19eab",
    "tower_mask": "86a48919e6e0bc00fe9aa7f055a98fb202f51cdd6f8f7cc5e6d2677c922915cd",
}
EXPECTED_APPLY_MODEL_HASHES = {
    "ppg_apply_model_base": "0ec432081df6bd5cbdc68c948c4c325331bc0220834b4720429c2049c86880e1",
    "ppg_apply_model_base_vr": "8b70b9bda2430fa7147694ebb8f51de72122ee525a8b5b59dca0d50c4c86000f",
    "ppg_apply_model_base_v0": "d71acf56911f4648baecf17c9bd567fbfa6f646e31ce20300ea1801fe77f15d4",
    "ppg_apply_model_base_v1": "56f59bb0d80f47726cf425423169bd6556f15353c0b2204e3124d84f9d3c579f",
    "ppg_apply_model_base_v2": "d863f545a2d9d8557243ddea38a3a3165c83c11d756eb65e1e998e90d68c75bf",
    "ppg_apply_model_base_v3": "3a722fb2c16f0120d62e74710963187f8d2efad9b0868ea411cff33055052393",
    "ppg_apply_model_base_E": "7e2d5ed9d1216ab30b92a137482a5bfa111d78e86d4f5a459179c38bc7993997",
    "ppg_apply_model_base_v0E": "6a302d7ebece4f5a592a38edb8de934ebb3a2935fdedff4d7bf012d54c7870dd",
    "ppg_apply_model_base_v1E": "b75c9e3c3c4a6e6333c79567b8de0faacf8ece6591e085813be974af6f0571e9",
    "ppg_apply_model_base_v2E": "c5a14d44b3655516692b012f15f2f84419f1acba70479d23c82421b6eeb49f27",
    "ppg_apply_model_base_v3E": "7679e634260402fb3815b2733767182690eec7587f9e09bffc307a05d00d59df",
    "ppg_apply_npb_model": "d6086dadac534013cda15cdfb69c1683776d3456d9e439903589653e8ac19eab",
}
STALE_LOCAL_PPG_MACRO_HASHES = {
    "28f9f27f36dde1b3ce64b2e6f1262efbcc9144b25d261a114693596d7a135442",
    "28c65527599e3bc03d0938cbe01ec3997955208037a59726349904512aadcdf3",
}
FORBIDDEN_RECOIL_PREFIXES = (
    "/sphenix/u/patsfan753/thesisAnalysis/install",
    "/sphenix/user/patsfan753/install",
)
REQUIRED_RECOIL_ROLES = {
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
    "ppg_recoeff_canonical_config",
    "ppg_calculate_photon_yield",
    "ppg_apply_bdt_macro",
    "ppg_apply_bdt_config",
    "ppg_recoeff_yaml_cpp",
    "ppg_recoeff_roounfold",
    "ppg_recoeff_roounfold_response_header",
    "ppg_recoeff_roounfold_bayes_header",
    "ppg_recoeff_vertex_scan_data",
    "ppg_recoeff_mbd_correction",
} | set(EXPECTED_APPLY_MODEL_HASHES)
EXPECTED_TRACE_OPERATIONS = (
    "event_identity_reader",
    "trace_stream_setup",
    "capture_raw_isolation",
    "truth_vertex_factor_scope",
    "capture_truth_vertex_factor",
    "trace_model_scope",
    "reuse_scoped_model",
    "reuse_scoped_score",
    "candidate_trace_emit",
    "response_trace_emit",
)


class AuditFailure(RuntimeError):
    pass


class RootTRandom3:
    """Minimal, dependency-free ROOT TRandom3 integer stream.

    The paired oracle needs the exact PHRandomSeed call sequence before ROOT
    is imported.  ROOT's fixed-seed PHRandomSeed path is a TRandom3 stream of
    ``Integer(UINT_MAX) + 1`` values.  Keeping this implementation here makes
    the seed/pedestal contract derivable and testable instead of hard-coding a
    guessed pedestal index.
    """

    _N = 624
    _M = 397
    _MATRIX_A = 0x9908B0DF
    _UPPER_MASK = 0x80000000
    _LOWER_MASK = 0x7FFFFFFF

    def __init__(self, seed: int):
        if not 1 <= seed <= 0xFFFFFFFF:
            raise AuditFailure("TRandom3 seed must be in [1, 4294967295]")
        self._mt = [0] * self._N
        self._mt[0] = seed & 0xFFFFFFFF
        for index in range(1, self._N):
            previous = self._mt[index - 1]
            self._mt[index] = (
                1812433253 * (previous ^ (previous >> 30)) + index
            ) & 0xFFFFFFFF
        self._index = self._N

    def _twist(self) -> None:
        for index in range(self._N):
            value = (
                (self._mt[index] & self._UPPER_MASK)
                | (self._mt[(index + 1) % self._N] & self._LOWER_MASK)
            )
            twisted = value >> 1
            if value & 1:
                twisted ^= self._MATRIX_A
            self._mt[index] = self._mt[(index + self._M) % self._N] ^ twisted
        self._index = 0

    def _uint32(self) -> int:
        if self._index >= self._N:
            self._twist()
        value = self._mt[self._index]
        self._index += 1
        value ^= value >> 11
        value ^= (value << 7) & 0x9D2C5680
        value ^= (value << 15) & 0xEFC60000
        value ^= value >> 18
        return value & 0xFFFFFFFF

    def integer(self, maximum: int) -> int:
        if maximum <= 0:
            return 0
        while True:
            value = self._uint32()
            if value != 0:
                break
        return int(maximum * (value * (1.0 / 4294967296.0)))


def historical_rng_contract() -> dict[str, Any]:
    return {
        "mode": "historical_fifo_replay_v2",
        "reco_consts_randomseed": "absent",
        "ph_seed_sequence": HISTORICAL_PH_SEED_SEQUENCE,
        "first_ph_seed": HISTORICAL_PH_SEED_SEQUENCE[0],
        "pedestal_seed": HISTORICAL_PH_SEED_SEQUENCE[1],
        "pedestal_sequence": HISTORICAL_PEDESTAL_SEQUENCE,
        "pedestal_file": HISTORICAL_PEDESTAL_FILE,
        "source_log": HISTORICAL_RNG_SOURCE,
        "source_log_call_count": 5,
    }


def lane_id_from_contract(lane: dict[str, Any]) -> str:
    sample = str(lane.get("sample", ""))
    period = str(lane.get("period", ""))
    interaction = str(lane.get("interaction", ""))
    if (
        sample not in VALID_SAMPLES
        or period not in VALID_PERIODS
        or interaction not in VALID_INTERACTIONS
        or lane.get("rows") != 5
    ):
        raise AuditFailure(f"unsupported physical lane contract: {lane}")
    return f"photon:{sample.lower()}:{period}:{interaction.lower()}"


def validate_authorization_token(
    contract_path: Path,
    contract: dict[str, Any],
    *,
    lane_id: str,
) -> None:
    """Replay the plan token from the exact physical lane and file digests."""
    lane = contract["lane"]
    paths = contract.get("paths", {})
    authorization = contract.get("authorization", {})
    token_files = authorization.get("token_bound_files", {})
    reuse = authorization.get("raw_ppg12_reuse", {})
    reuse_mode = reuse.get("mode")
    required_path_roles = (
        "output_dir", "setup_script", "ppg_macro", "g4_full_list",
        "truthjet_full_list", "apply_bdt", "apply_config", "base_e_model",
        "base_v3e_model", "npb_model", "tower_mask",
        "recoil_runtime_manifest", "recoil_config",
    )
    required_values = {"output_dir": str(contract_path.parent)}
    for role in required_path_roles[1:]:
        value = str(paths.get(role, ""))
        if not value:
            raise AuditFailure(f"contract lacks token-bound path role {role}")
        required_values[role] = value
    rng = contract.get("rng", {})
    ph_sequence = ",".join(str(value) for value in rng.get("ph_seed_sequence", []))
    token_lines = [
        "schema_version=3",
        f"lane={lane['sample']}:{lane['period']}:{lane['interaction']}",
        f"lane_id={lane_id}",
        "rows=5",
        "rng_mode=historical_fifo_replay_v2",
        f"first_ph_seed={rng.get('first_ph_seed')}",
        f"pedestal_seed={rng.get('pedestal_seed')}",
        f"pedestal={rng.get('pedestal_sequence')}",
        f"ph_seed_sequence={ph_sequence}",
        "source_graph=NONE,g4,truthjet,NONE,NONE",
        "runtime=new.17",
        f"reuse_mode={reuse_mode}",
    ]
    token_lines.extend(
        f"{role}={required_values[role]}" for role in required_path_roles
    )
    for role in sorted(token_files):
        item = token_files[role]
        token_lines.append(
            f"{role}={item.get('path', '')} sha256={item.get('sha256', '')}"
        )
    payload = "".join(f"{line}\n" for line in sorted(token_lines)).encode()
    expected = "ppg12-oracle:" + hashlib.sha256(payload).hexdigest()
    if contract.get("authorization_token") != expected:
        raise AuditFailure(
            "paired-oracle authorization token does not replay from its lane and assets"
        )


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    try:
        data = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise AuditFailure(f"cannot read JSON {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise AuditFailure(f"JSON root must be an object: {path}")
    return data


def require_file(path_value: str, label: str) -> Path:
    if not path_value:
        raise AuditFailure(f"missing path for {label}")
    path = Path(path_value)
    if not path.is_absolute():
        raise AuditFailure(f"{label} must be absolute: {path}")
    if not path.is_file() or path.stat().st_size <= 0:
        raise AuditFailure(f"{label} is missing or empty: {path}")
    return path


def require_hash(path: Path, expected: str, label: str) -> str:
    actual = sha256(path)
    if actual in STALE_LOCAL_PPG_MACRO_HASHES:
        raise AuditFailure(f"{label} resolves to a known stale local PPG macro: {actual}")
    if actual != expected:
        raise AuditFailure(
            f"{label} hash mismatch: expected {expected}, observed {actual}, path={path}"
        )
    return actual


def noncomment_rows(path: Path) -> list[str]:
    return [
        line.strip()
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]


def source_identity(value: str) -> tuple[str, str]:
    match = re.search(r"-(\d{10})-(\d{6})\.root$", Path(value).name)
    if not match:
        raise AuditFailure(f"cannot recover run/segment identity from source: {value}")
    return match.group(1), match.group(2)


def validate_source_graph(g4_slice: Path, truth_slice: Path, combined: Path) -> None:
    g4_rows = noncomment_rows(g4_slice)
    truth_rows = noncomment_rows(truth_slice)
    combined_rows = noncomment_rows(combined)
    if len(g4_rows) != 5 or len(truth_rows) != 5 or len(combined_rows) != 5:
        raise AuditFailure(
            "paired oracle requires exactly five G4, truth-jet, and combined rows"
        )
    if len(set(g4_rows)) != 5 or len(set(truth_rows)) != 5:
        raise AuditFailure("source slices contain duplicate rows")
    for index, (g4, truth, line) in enumerate(
        zip(g4_rows, truth_rows, combined_rows), start=1
    ):
        if source_identity(g4) != source_identity(truth):
            raise AuditFailure(f"row {index} G4/truth-jet identity mismatch")
        fields = line.split()
        if len(fields) != 5:
            raise AuditFailure(f"combined row {index} does not have exactly five columns")
        expected = ["NONE", g4, truth, "NONE", "NONE"]
        if fields != expected:
            raise AuditFailure(
                f"combined row {index} violates exact oracle graph; expected {expected}, got {fields}"
            )


def validate_frozen_macro_graph(path: Path) -> None:
    text = path.read_text(errors="replace")
    required_active = (
        "INPUTREADHITS::listfile[0] = inputFile0;",
        "INPUTREADHITS::listfile[4] = inputFile4;",
    )
    for line in required_active:
        if line not in text:
            raise AuditFailure(f"preserved macro is missing active source line: {line}")
    for index in (1, 2, 3):
        active = re.search(
            rf"^(?!\s*//)\s*INPUTREADHITS::listfile\[{index}\]\s*=",
            text,
            flags=re.MULTILINE,
        )
        if active:
            raise AuditFailure(
                f"preserved macro unexpectedly activates source index {index}; "
                "four-stream SI is not the executable oracle"
            )


def is_below(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
        return True
    except ValueError:
        return False


def validate_recoil_manifest(path: Path, expected_offline: str) -> dict[str, Path]:
    data = load_json(path)
    if data.get("schema_version") != 1:
        raise AuditFailure("Recoil runtime manifest schema_version must be 1")
    if data.get("runtime_profile") != EXPECTED_RUNTIME_PROFILE:
        raise AuditFailure("Recoil runtime manifest is not pinned to new.17")
    if data.get("offline_main") != expected_offline:
        raise AuditFailure("Recoil runtime manifest OFFLINE_MAIN differs from oracle runtime")
    if data.get("isolated_build") is not True:
        raise AuditFailure("Recoil runtime manifest must declare isolated_build=true")
    if data.get("estimator_revision") != EXPECTED_ESTIMATOR_REVISION:
        raise AuditFailure("Recoil runtime manifest uses the wrong estimator revision")
    build_receipt = require_file(
        data.get("build_receipt", ""), "Recoil build receipt"
    )
    receipt_digest = str(data.get("build_receipt_sha256", ""))
    if not re.fullmatch(r"[0-9a-f]{64}", receipt_digest):
        raise AuditFailure("Recoil runtime manifest lacks build_receipt_sha256")
    require_hash(build_receipt, receipt_digest, "Recoil build receipt")
    receipt = load_json(build_receipt)
    ppg_source = receipt.get("ppg12_source", {})
    if not isinstance(ppg_source, dict):
        raise AuditFailure("runtime receipt lacks source-locked PPG12 provenance")
    if ppg_source.get("revision") != EXPECTED_PPG_SOURCE_REVISION:
        raise AuditFailure("runtime receipt uses the wrong PPG12 source revision")
    if ppg_source.get("working_tree_ignored") is not True:
        raise AuditFailure("runtime receipt does not exclude mutable PPG12 working-tree edits")
    if ppg_source.get("rebuilt_against_common_runtime") is not True:
        raise AuditFailure("runtime receipt does not rebuild PPG12 against the common runtime")
    rewrites = receipt.get("staged_rewrites", {})
    if not isinstance(rewrites, dict):
        raise AuditFailure("runtime receipt lacks staged-rewrite provenance")
    if rewrites.get("archived_ppg12_binary_reused") is not False:
        raise AuditFailure("runtime receipt reuses the ABI-incompatible archived PPG12 binary")
    if rewrites.get("ppg12_source_locked_rebuild") is not True:
        raise AuditFailure("runtime receipt lacks the source-locked PPG12 rebuild gate")
    estimator = receipt.get("ppg12_estimator", {})
    if not isinstance(estimator, dict):
        raise AuditFailure("runtime receipt lacks preserved estimator provenance")
    if estimator.get("revision") != EXPECTED_ESTIMATOR_REVISION:
        raise AuditFailure("runtime receipt uses the wrong estimator revision")
    if estimator.get("working_tree_ignored") is not True:
        raise AuditFailure("runtime receipt does not exclude mutable estimator edits")
    if estimator.get("reconstruction_revision_is_distinct") is not True:
        raise AuditFailure("runtime receipt conflates reconstruction and estimator revisions")
    if rewrites.get("estimator_revision_separate_from_reconstruction") is not True:
        raise AuditFailure("runtime receipt lacks the separate-estimator-revision gate")
    if rewrites.get("estimator_trace_requires_exact_root_equivalence") is not True:
        raise AuditFailure("runtime receipt does not require dual ROOT equivalence")
    isolated_root = build_receipt.parent
    files = data.get("files")
    if not isinstance(files, list):
        raise AuditFailure("Recoil runtime manifest files must be a list")
    by_role: dict[str, Path] = {}
    for item in files:
        if not isinstance(item, dict):
            raise AuditFailure("invalid Recoil runtime file entry")
        role = str(item.get("role", ""))
        file_path = require_file(str(item.get("path", "")), f"Recoil {role}")
        expected = str(item.get("sha256", ""))
        if not re.fullmatch(r"[0-9a-f]{64}", expected):
            raise AuditFailure(f"invalid expected digest for Recoil {role}")
        require_hash(file_path, expected, f"Recoil {role}")
        if role in by_role:
            raise AuditFailure(f"duplicate Recoil runtime role: {role}")
        by_role[role] = file_path
        resolved = str(file_path.resolve())
        if any(resolved.startswith(prefix + "/") for prefix in FORBIDDEN_RECOIL_PREFIXES):
            raise AuditFailure(f"Recoil oracle file escapes isolated build: {file_path}")
        if not is_below(file_path, isolated_root):
            raise AuditFailure(
                f"Recoil oracle file is outside build-receipt root {isolated_root}: {file_path}"
            )
    missing = sorted(REQUIRED_RECOIL_ROLES - set(by_role))
    if missing:
        raise AuditFailure(f"Recoil runtime manifest lacks required roles: {missing}")

    source_meta = estimator.get("canonical_source", {})
    baseline_meta = estimator.get("staged_uninstrumented_macro", {})
    trace_meta = estimator.get("instrumented_macro", {})
    config_meta = estimator.get("canonical_config", {})
    yield_meta = estimator.get("calculate_photon_yield", {})
    for label, metadata, role in (
        ("canonical estimator source", source_meta, "ppg_recoeff_source_macro"),
        ("uninstrumented estimator macro", baseline_meta, "ppg_recoeff_macro"),
        ("instrumented estimator macro", trace_meta, "ppg_recoeff_trace_macro"),
        ("canonical estimator config", config_meta, "ppg_recoeff_canonical_config"),
        ("CalculatePhotonYield source", yield_meta, "ppg_calculate_photon_yield"),
    ):
        if not isinstance(metadata, dict):
            raise AuditFailure(f"runtime receipt lacks {label} metadata")
        if Path(str(metadata.get("path", ""))).resolve() != by_role[role].resolve():
            raise AuditFailure(f"runtime receipt {label} path differs from manifest role")
        if str(metadata.get("sha256", "")) != sha256(by_role[role]):
            raise AuditFailure(f"runtime receipt {label} digest differs from manifest role")
    apply_meta = estimator.get("apply_bdt_stage", {})
    if not isinstance(apply_meta, dict):
        raise AuditFailure("runtime receipt lacks source-locked apply_BDT metadata")
    if apply_meta.get("revision") != EXPECTED_ESTIMATOR_REVISION:
        raise AuditFailure("apply_BDT stage uses the wrong source revision")
    if apply_meta.get("working_tree_ignored") is not True:
        raise AuditFailure("apply_BDT stage does not exclude mutable checkout edits")
    if apply_meta.get("scored_root_shared_by_both_recoeff_runs") is not True:
        raise AuditFailure("apply_BDT stage does not require one shared scored ROOT")
    for key, role in (
        ("macro", "ppg_apply_bdt_macro"),
        ("config", "ppg_apply_bdt_config"),
    ):
        metadata = apply_meta.get(key, {})
        if not isinstance(metadata, dict):
            raise AuditFailure(f"apply_BDT receipt lacks {key} metadata")
        if Path(str(metadata.get("path", ""))).resolve() != by_role[role].resolve():
            raise AuditFailure(f"apply_BDT receipt {key} path differs from manifest role")
        if str(metadata.get("sha256", "")) != sha256(by_role[role]):
            raise AuditFailure(f"apply_BDT receipt {key} digest differs from manifest role")
    model_meta = apply_meta.get("model_assets", [])
    if not isinstance(model_meta, list) or len(model_meta) != len(EXPECTED_APPLY_MODEL_HASHES):
        raise AuditFailure("apply_BDT receipt does not contain exactly 11 split models plus NPB")
    receipt_model_hashes = {
        str(item.get("sha256", ""))
        for item in model_meta
        if isinstance(item, dict)
    }
    if receipt_model_hashes != set(EXPECTED_APPLY_MODEL_HASHES.values()):
        raise AuditFailure("apply_BDT receipt model hashes differ from the canonical split assets")
    receipt_model_paths = {
        Path(str(item.get("path", ""))).resolve()
        for item in model_meta
        if isinstance(item, dict)
    }
    manifest_model_paths = {
        by_role[role].resolve() for role in EXPECTED_APPLY_MODEL_HASHES
    }
    if receipt_model_paths != manifest_model_paths:
        raise AuditFailure("apply_BDT receipt model paths differ from manifest roles")
    for role, expected_hash in EXPECTED_APPLY_MODEL_HASHES.items():
        require_hash(by_role[role], expected_hash, role)
    if baseline_meta.get("scientific_expression_changes") is not False:
        raise AuditFailure("uninstrumented estimator macro declares scientific changes")
    if baseline_meta.get("sealed_path_rewrites_only") is not True:
        raise AuditFailure("uninstrumented estimator macro is not path-rewrite-only")
    if trace_meta.get("candidate_trace_side_channel_only") is not True:
        raise AuditFailure("instrumented estimator trace is not declared side-channel-only")
    if trace_meta.get("requires_exact_root_equivalence") is not True:
        raise AuditFailure("instrumented estimator does not require exact ROOT equivalence")
    require_hash(
        by_role["ppg_recoeff_source_macro"],
        EXPECTED_RECOEFF_SOURCE_HASH,
        "canonical RecoEff source",
    )
    require_hash(
        by_role["ppg_recoeff_canonical_config"],
        EXPECTED_RECOEFF_CONFIG_HASH,
        "canonical RecoEff config",
    )

    transform_path = by_role["ppg_recoeff_trace_transform_receipt"]
    if Path(str(trace_meta.get("transform_receipt", ""))).resolve() != transform_path.resolve():
        raise AuditFailure("instrumented estimator transform receipt path differs")
    if str(trace_meta.get("transform_receipt_sha256", "")) != sha256(transform_path):
        raise AuditFailure("instrumented estimator transform receipt digest differs")
    transform = load_json(transform_path)
    if transform.get("schema_version") != 1:
        raise AuditFailure("estimator transform receipt schema_version must be 1")
    if transform.get("transform") != "ppg12_recoeff_candidate_trace_v1":
        raise AuditFailure("unexpected estimator trace transform")
    if transform.get("source_revision") != EXPECTED_ESTIMATOR_REVISION:
        raise AuditFailure("estimator trace transform uses the wrong source revision")
    if transform.get("trace_side_channel_only") is not True:
        raise AuditFailure("estimator transform is not side-channel-only")
    if transform.get("selection_or_fill_expression_replaced") is not False:
        raise AuditFailure("estimator transform replaced a selection or fill expression")
    if Path(str(transform.get("input_path", ""))).resolve() != by_role["ppg_recoeff_macro"].resolve():
        raise AuditFailure("estimator transform input path differs from baseline macro")
    if transform.get("input_sha256") != sha256(by_role["ppg_recoeff_macro"]):
        raise AuditFailure("estimator transform input digest differs from baseline macro")
    if Path(str(transform.get("output_path", ""))).resolve() != by_role["ppg_recoeff_trace_macro"].resolve():
        raise AuditFailure("estimator transform output path differs from trace macro")
    if transform.get("output_sha256") != sha256(by_role["ppg_recoeff_trace_macro"]):
        raise AuditFailure("estimator transform output digest differs from trace macro")
    operations = transform.get("operations", [])
    if not isinstance(operations, list) or tuple(
        str(item.get("label", "")) for item in operations if isinstance(item, dict)
    ) != EXPECTED_TRACE_OPERATIONS:
        raise AuditFailure("estimator trace transform operation sequence differs")
    header = by_role["PhotonClusterBuilder.h"]
    if header.parent.name != "caloreco" or header.parent.parent.name != "include":
        raise AuditFailure(
            "isolated PhotonClusterBuilder.h must be staged as "
            "<isolated-prefix>/include/caloreco/PhotonClusterBuilder.h"
        )
    return by_role


def validate_ldd(library: Path, expected_offline: str) -> None:
    proc = subprocess.run(
        ["ldd", str(library)], text=True, capture_output=True, check=False
    )
    if proc.returncode != 0:
        raise AuditFailure(f"ldd failed for {library}: {proc.stderr.strip()}")
    for line in proc.stdout.splitlines():
        if "not found" in line:
            raise AuditFailure(f"unresolved dependency for {library}: {line.strip()}")
        match = re.search(r"=>\s+(/\S+)", line)
        if not match:
            continue
        dependency = match.group(1)
        if "/release/release_" in dependency and not dependency.startswith(
            expected_offline + "/"
        ):
            raise AuditFailure(
                f"{library} resolves a dependency outside new.17: {dependency}"
            )


def validate_contract(contract_path: Path, run_ldd: bool = True) -> dict[str, Any]:
    contract = load_json(contract_path)
    if contract.get("schema_version") != SCHEMA_VERSION:
        raise AuditFailure("paired-oracle contract schema_version must be 3")
    lane = contract.get("lane")
    if not isinstance(lane, dict):
        raise AuditFailure("paired-oracle lane contract must be an object")
    lane_id = lane_id_from_contract(lane)
    if lane.get("lane_id") != lane_id:
        raise AuditFailure("paired-oracle embedded lane_id differs from physical lane")
    if "seed" in lane:
        raise AuditFailure(
            "historical replay must not contain a synthetic reconstruction seed"
        )
    expected_rng = historical_rng_contract()
    if contract.get("rng") != expected_rng:
        raise AuditFailure("paired-oracle RNG metadata is not derived from its seed")
    runtime = contract.get("runtime", {})
    if runtime.get("profile") != EXPECTED_RUNTIME_PROFILE:
        raise AuditFailure("paired oracle runtime profile must be new.17")
    if runtime.get("offline_main") != EXPECTED_OFFLINE_MAIN:
        raise AuditFailure("paired oracle OFFLINE_MAIN is not the frozen new.17 path")
    if runtime.get("actual_offline_main") != EXPECTED_OFFLINE_MAIN:
        raise AuditFailure("active OFFLINE_MAIN does not resolve to frozen new.17")

    paths = contract.get("paths", {})
    authorization = contract.get("authorization", {})
    if not isinstance(authorization, dict) or authorization.get(
        "token_schema_version"
    ) != 3:
        raise AuditFailure("paired-oracle authorization is not content-bound schema 3")
    token_files = authorization.get("token_bound_files", {})
    required_token_roles = {
        "setup_script", "ppg_macro", "g4_full_list", "truthjet_full_list",
        "apply_bdt", "apply_config", "base_e_model", "base_v3e_model",
        "npb_model", "tower_mask", "recoil_runtime_manifest", "recoil_config",
        "driver_script", "worker_script", "ppg_wrapper", "recoil_wrapper",
        "comparator", "auditor", "aggregate_extractor",
    }
    reuse = authorization.get("raw_ppg12_reuse", {})
    if not isinstance(reuse, dict) or reuse.get("mode") not in {
        "disabled", "exact_contract_bound",
    }:
        raise AuditFailure("raw PPG12 reuse authorization is invalid")
    if reuse.get("mode") == "exact_contract_bound":
        required_token_roles.update({"reuse_ppg_raw_root", "reuse_ppg_raw_contract"})
        for key in ("source", "contract"):
            item = reuse.get(key)
            if not isinstance(item, dict):
                raise AuditFailure(f"raw PPG12 reuse lacks token-bound {key}")
            file_path = require_file(str(item.get("path", "")), f"reuse {key}")
            require_hash(file_path, str(item.get("sha256", "")), f"reuse {key}")
    elif reuse.get("source") is not None or reuse.get("contract") is not None:
        raise AuditFailure("disabled raw PPG12 reuse unexpectedly names an artifact")
    if not isinstance(token_files, dict) or set(token_files) != required_token_roles:
        raise AuditFailure("paired-oracle token-bound file role set differs")
    for role, item in token_files.items():
        if not isinstance(item, dict):
            raise AuditFailure(f"invalid token-bound metadata for {role}")
        file_path = require_file(str(item.get("path", "")), f"token-bound {role}")
        expected_hash = str(item.get("sha256", ""))
        if not re.fullmatch(r"[0-9a-f]{64}", expected_hash):
            raise AuditFailure(f"invalid token-bound digest for {role}")
        require_hash(file_path, expected_hash, f"token-bound {role}")
        if role in paths and file_path.resolve() != Path(str(paths[role])).resolve():
            raise AuditFailure(f"token-bound {role} path differs from contract paths")
    validate_authorization_token(contract_path, contract, lane_id=lane_id)

    ppg_macro = require_file(paths.get("ppg_macro", ""), "preserved PPG12 macro")
    if {key: lane.get(key) for key in LEGACY_PINNED_LANE} == LEGACY_PINNED_LANE:
        require_hash(ppg_macro, EXPECTED_HASHES["ppg_macro"], "preserved PPG12 macro")
    validate_frozen_macro_graph(ppg_macro)
    g4_full = require_file(paths.get("g4_full_list", ""), "G4 full list")
    truth_full = require_file(paths.get("truthjet_full_list", ""), "truth-jet full list")
    if {key: lane.get(key) for key in LEGACY_PINNED_LANE} == LEGACY_PINNED_LANE:
        require_hash(g4_full, EXPECTED_HASHES["g4_full_list"], "G4 full list")
        require_hash(
            truth_full, EXPECTED_HASHES["truthjet_full_list"], "truth-jet full list"
        )
    g4_slice = require_file(paths.get("g4_slice", ""), "G4 five-row slice")
    truth_slice = require_file(
        paths.get("truthjet_slice", ""), "truth-jet five-row slice"
    )
    combined = require_file(paths.get("combined_list", ""), "Recoil combined list")
    if noncomment_rows(g4_slice) != noncomment_rows(g4_full)[:5]:
        raise AuditFailure("G4 slice is not exactly the first five full-list rows")
    if noncomment_rows(truth_slice) != noncomment_rows(truth_full)[:5]:
        raise AuditFailure(
            "truth-jet slice is not exactly the first five full-list rows"
        )
    expected_identities = [("0000000028", f"{segment:06d}") for segment in range(5)]
    observed_identities = [source_identity(row) for row in noncomment_rows(g4_slice)]
    if observed_identities != expected_identities:
        raise AuditFailure(
            "source slice is not the frozen run-28 segments 000000--000004"
        )
    validate_source_graph(g4_slice, truth_slice, combined)

    for role in (
        "apply_bdt",
        "apply_config",
        "base_e_model",
        "base_v3e_model",
        "npb_model",
        "tower_mask",
    ):
        asset = require_file(paths.get(role, ""), role)
        require_hash(asset, EXPECTED_HASHES[role], role)

    recoil_manifest = require_file(
        paths.get("recoil_runtime_manifest", ""), "Recoil runtime manifest"
    )
    recoil_files = validate_recoil_manifest(recoil_manifest, EXPECTED_OFFLINE_MAIN)
    if Path(paths["apply_bdt"]).resolve() != recoil_files["ppg_apply_bdt_macro"].resolve():
        raise AuditFailure("contract apply_BDT macro differs from sealed runtime role")
    if Path(paths["apply_config"]).resolve() != recoil_files["ppg_apply_bdt_config"].resolve():
        raise AuditFailure("contract apply_BDT config differs from sealed runtime role")
    for contract_role, manifest_role in (
        ("base_e_model", "ppg_apply_model_base_E"),
        ("base_v3e_model", "ppg_apply_model_base_v3E"),
        ("npb_model", "ppg_apply_npb_model"),
    ):
        if Path(paths[contract_role]).resolve() != recoil_files[manifest_role].resolve():
            raise AuditFailure(f"contract {contract_role} differs from sealed runtime role")

    apply_runtime_config = require_file(
        paths.get("apply_bdt_runtime_config", ""), "derived apply_BDT runtime config"
    )
    canonical_apply_text = Path(paths["apply_config"]).read_text()
    no_split_block = (
        '    - node: "CLUSTERINFO_CEMC_NO_SPLIT"\n'
        '      model_suffix: "_nosplit"\n'
    )
    if canonical_apply_text.count(no_split_block) != 1:
        raise AuditFailure("canonical apply_BDT no-split block changed")
    if apply_runtime_config.read_text() != canonical_apply_text.replace(no_split_block, ""):
        raise AuditFailure("apply_BDT runtime config contains changes beyond omitting no-split")

    canonical_recoeff_config = recoil_files["ppg_recoeff_canonical_config"]
    baseline_config = require_file(
        paths.get("ppg_recoeff_baseline_config", ""), "baseline RecoEff config"
    )
    trace_config = require_file(
        paths.get("ppg_recoeff_trace_config", ""), "trace RecoEff config"
    )
    if baseline_config.read_bytes() != trace_config.read_bytes():
        raise AuditFailure("baseline and instrumented RecoEff configs differ")
    expected_config = canonical_recoeff_config.read_text()
    for old, new in (
        ('photon_jet_file_root_dir: "/sphenix/user/shuhangli/ppg12/FunWithxgboost/"',
         'photon_jet_file_root_dir: "input/"'),
        ('eff_outfile: "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency"',
         'eff_outfile: "output/MC_efficiency"'),
        ('response_outfile: "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_response"',
         'response_outfile: "output/MC_response"'),
        ('data_outfile: "/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo"',
         'data_outfile: "output/data_histo"'),
        ('var_type: "bdt_nom"', 'var_type: "paired_oracle"'),
        ('vertex_scan_data_file: ""',
         'vertex_scan_data_file: "input/data_histo_bdt_nom_vtxscan.root"'),
        ('tower_mask_file: "/sphenix/user/shuhangli/ppg12/efficiencytool/tower_masks_bdt_nom.root"',
         'tower_mask_file: "input/tower_masks_bdt_nom.root"'),
    ):
        if expected_config.count(old) != 1:
            raise AuditFailure(f"canonical RecoEff config rewrite marker changed: {old}")
        expected_config = expected_config.replace(old, new)
    if baseline_config.read_text() != expected_config:
        raise AuditFailure("RecoEff runtime config contains non-path scientific changes")
    ppg_lib = require_file(paths.get("ppg_caloana24", ""), "source-locked libCaloAna24")
    if recoil_files["libCaloAna24.so"].resolve() != ppg_lib.resolve():
        raise AuditFailure("contract PPG12 library differs from isolated runtime manifest")
    recoil_macro = require_file(paths.get("recoil_macro", ""), "Recoil oracle macro")
    if recoil_files["recoil_macro"].resolve() != recoil_macro.resolve():
        raise AuditFailure("contract Recoil macro differs from isolated runtime manifest")
    if run_ldd:
        validate_ldd(ppg_lib, EXPECTED_OFFLINE_MAIN)
        for role in (
            "libcalo_reco.so",
            "libRecoilJets.so",
            "libclusteriso.so",
            "libjetbase.so",
        ):
            validate_ldd(recoil_files[role], EXPECTED_OFFLINE_MAIN)

    calo_calib = require_file(paths.get("calo_calib", ""), "resolved Calo_Calib.C")
    runtime_hashes = contract.get("runtime_hashes", {})
    for role, file_path in (
        ("setup_script", require_file(paths.get("setup_script", ""), "setup script")),
        ("calo_calib", calo_calib),
        ("recoil_config", require_file(paths.get("recoil_config", ""), "Recoil config")),
        ("ppg_wrapper", require_file(paths.get("ppg_wrapper", ""), "PPG12 wrapper")),
        ("recoil_wrapper", require_file(paths.get("recoil_wrapper", ""), "Recoil wrapper")),
        ("comparator", require_file(paths.get("comparator", ""), "oracle comparator")),
        ("auditor", require_file(paths.get("auditor", ""), "oracle auditor")),
        ("aggregate_extractor", require_file(paths.get("aggregate_extractor", ""), "aggregate extractor")),
        ("ppg_recoeff_baseline_config", baseline_config),
        ("ppg_recoeff_trace_config", trace_config),
        ("apply_bdt_runtime_config", apply_runtime_config),
    ):
        expected = str(runtime_hashes.get(role, ""))
        if not re.fullmatch(r"[0-9a-f]{64}", expected):
            raise AuditFailure(f"missing runtime fingerprint for {role}")
        require_hash(file_path, expected, role)
    return contract


def ph_seed_sequence(log_text: str) -> list[int]:
    values = [
        int(value)
        for value in re.findall(r"PHRandomSeed::GetSeed\(\) seed:\s*(\d+)", log_text)
    ]
    if len(values) != 5:
        raise AuditFailure(
            f"historical replay requires exactly five PHRandomSeed records; got {len(values)}"
        )
    return values


def actual_pedestal(log_text: str) -> str:
    for line in log_text.splitlines():
        if line.startswith("ORACLE_SEED_CONTRACT"):
            continue
        match = re.search(r"(pedestal-54256-0\d{4}\.root)", line)
        if match:
            return match.group(1)
    raise AuditFailure("log does not expose the actual pedestal file")


def parse_keyed_line(log_text: str, prefix: str, side: str) -> dict[str, str]:
    for line in log_text.splitlines():
        if not line.startswith(prefix + " "):
            continue
        fields = dict(re.findall(r"([A-Za-z0-9_]+)=([^\s]+)", line))
        if fields.get("side") == side:
            return fields
    raise AuditFailure(f"missing {prefix} line for {side}")


def library_block(log_text: str, side: str) -> str:
    pattern = re.compile(
        rf"ORACLE_LIBRARIES_BEGIN side={re.escape(side)}\n(.*?)\n"
        rf"ORACLE_LIBRARIES_END side={re.escape(side)}",
        flags=re.DOTALL,
    )
    match = pattern.search(log_text)
    if not match:
        raise AuditFailure(f"missing loaded-library block for {side}")
    return match.group(1)


def validate_log(
    path: Path,
    side: str,
    calo_calib: Path,
    rng: dict[str, Any],
    manifest_files: dict[str, Path] | None = None,
) -> None:
    text = path.read_text(errors="replace")
    runtime = parse_keyed_line(text, "ORACLE_RUNTIME", side)
    if runtime.get("profile") != EXPECTED_RUNTIME_PROFILE:
        raise AuditFailure(f"{side} log does not report new.17")
    if runtime.get("offline_main") != EXPECTED_OFFLINE_MAIN:
        raise AuditFailure(f"{side} log reports the wrong OFFLINE_MAIN")
    seed_line = parse_keyed_line(text, "ORACLE_SEED_CONTRACT", side)
    if seed_line.get("mode") != rng["mode"]:
        raise AuditFailure(f"{side} log reports the wrong RNG replay mode")
    if seed_line.get("rc_randomseed") != "absent":
        raise AuditFailure(f"{side} replay unexpectedly sets recoConsts RANDOMSEED")
    if seed_line.get("ph_seed_sequence") != ",".join(
        str(value) for value in rng["ph_seed_sequence"]
    ):
        raise AuditFailure(f"{side} wrapper reports the wrong replay sequence")
    if seed_line.get("pedestal_sequence") != str(rng["pedestal_sequence"]):
        raise AuditFailure(f"{side} wrapper reports the wrong pedestal sequence")
    if ph_seed_sequence(text) != rng["ph_seed_sequence"]:
        raise AuditFailure(f"{side} PHRandomSeed call sequence differs from contract")
    if actual_pedestal(text) != rng["pedestal_file"]:
        raise AuditFailure(f"{side} actual pedestal differs from contract")
    graph = parse_keyed_line(text, "ORACLE_SOURCE_GRAPH", side)
    if graph.get("columns") != "NONE,g4,truthjet,NONE,NONE":
        raise AuditFailure(f"{side} log reports a non-oracle source graph")
    calo = parse_keyed_line(text, "ORACLE_CALO_CALIB", side)
    if Path(calo.get("path", "")).resolve() != calo_calib.resolve():
        raise AuditFailure(f"{side} loaded a different Calo_Calib.C")
    libraries = library_block(text, side)
    if "/release/release_ana/" in libraries:
        raise AuditFailure(f"{side} loaded an analysis release in the new.17 oracle")
    for other in re.findall(r"/release/release_new/(new\.[0-9]+)", libraries):
        if other != EXPECTED_RUNTIME_PROFILE:
            raise AuditFailure(f"{side} loaded mixed release {other}")

    if side == "ppg12":
        assert manifest_files is not None
        dynamic = parse_keyed_line(text, "ORACLE_DYNAMIC_LIBRARY", side)
        lib_path = require_file(dynamic.get("path", ""), "loaded libCaloAna24")
        if sha256(lib_path) != sha256(manifest_files["libCaloAna24.so"]):
            raise AuditFailure("loaded libCaloAna24 differs from isolated runtime manifest")
    else:
        assert manifest_files is not None
        for role in (
            "libcalo_reco.so",
            "libRecoilJets.so",
            "libclusteriso.so",
            "libjetbase.so",
        ):
            line = None
            for candidate in text.splitlines():
                if (
                    candidate.startswith("ORACLE_DYNAMIC_LIBRARY side=recoiljets ")
                    and f"name={role}" in candidate
                ):
                    line = candidate
                    break
            if line is None:
                raise AuditFailure(f"RecoilJets log lacks actual library path for {role}")
            fields = dict(re.findall(r"([A-Za-z0-9_.]+)=([^\s]+)", line))
            actual = require_file(fields.get("path", ""), f"loaded {role}")
            if sha256(actual) != sha256(manifest_files[role]):
                raise AuditFailure(f"loaded {role} differs from isolated runtime manifest")


def require_root_file(path: Path, tree_name: str, branches: Iterable[str]) -> None:
    try:
        import ROOT  # type: ignore
    except ImportError as exc:
        raise AuditFailure("PyROOT is required for postrun ROOT validation") from exc
    root_file = ROOT.TFile.Open(str(path))
    if not root_file or root_file.IsZombie():
        raise AuditFailure(f"unreadable or zombie ROOT file: {path}")
    if root_file.TestBit(ROOT.TFile.kRecovered):
        raise AuditFailure(f"recovered ROOT file is not admissible: {path}")
    tree = root_file.Get(tree_name)
    if not tree:
        raise AuditFailure(f"missing tree {tree_name} in {path}")
    present = {branch.GetName() for branch in tree.GetListOfBranches()}
    missing = sorted(set(branches) - present)
    if missing:
        raise AuditFailure(f"missing {tree_name} branches in {path}: {missing}")
    root_file.Close()


def as_float(row: dict[str, str], field: str) -> float:
    try:
        value = float(row[field])
    except (KeyError, ValueError) as exc:
        raise AuditFailure(f"candidate report has invalid {field}") from exc
    if not math.isfinite(value):
        raise AuditFailure(f"candidate report has non-finite {field}")
    return value


def as_int(row: dict[str, str], field: str) -> int:
    try:
        return int(row[field])
    except (KeyError, ValueError) as exc:
        raise AuditFailure(f"candidate report has invalid {field}") from exc


FIRST_DIVERGENCE_SCHEMA_VERSION = 1
PRESERVED_EXECUTABLE_EVIDENCE = "preserved_ppg12_executable"
FIRST_DIVERGENCE_STAGE_ORDER = (
    "population",
    "features",
    "scores",
    "route",
    "tags",
    "isolation_abcd",
    "truth_response_fills",
    "weights",
)
ORACLE_FEATURE_NAMES = (
    "cluster_Et_score_input",
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "vertexz",
    "cluster_Eta",
    "e11_over_e33",
    "cluster_et1",
    "cluster_et2",
    "cluster_et3",
    "cluster_et4",
    "e32_over_e35",
)


def _first_divergence(
    stage: str,
    csv_row: int,
    row: dict[str, str],
    field: str,
    reason: str,
    *,
    observed: Any = None,
    reference: Any = None,
    delta: float | None = None,
    tolerance: float | None = None,
) -> dict[str, Any]:
    finding: dict[str, Any] = {
        "stage": stage,
        "csv_row": csv_row,
        "candidate_identity": row.get("candidate_identity", ""),
        "field": field,
        "reason": reason,
        "observed": observed,
        "reference": reference,
    }
    if delta is not None:
        finding["delta"] = delta
    if tolerance is not None:
        finding["tolerance"] = tolerance
    return finding


def analyze_candidate_report(
    path: Path,
    *,
    expected_lane_id: str | None = None,
    expected_runtime_contract_sha256: str | None = None,
) -> dict[str, Any]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    stage_counts = {
        stage: {"evaluated": 0, "failed": 0}
        for stage in FIRST_DIVERGENCE_STAGE_ORDER
    }
    if not rows:
        stage_counts["population"]["failed"] = 1
        return {
            "schema_version": FIRST_DIVERGENCE_SCHEMA_VERSION,
            "candidate_report": str(path.resolve()),
            "status": "FAIL",
            "row_count": 0,
            "stage_order": list(FIRST_DIVERGENCE_STAGE_ORDER),
            "stage_counts": stage_counts,
            "first_divergence": {
                "stage": "population",
                "csv_row": 1,
                "candidate_identity": "",
                "field": "candidate_report",
                "reason": "candidate report contains no rows",
                "observed": 0,
                "reference": ">0",
            },
        }
    first: dict[str, Any] | None = None

    def fail(stage: str, index: int, row: dict[str, str], field: str, reason: str,
             **evidence: Any) -> bool:
        nonlocal first
        stage_counts[stage]["failed"] += 1
        first = _first_divergence(stage, index, row, field, reason, **evidence)
        return False

    identities: set[str] = set()
    stage = "population"
    for index, row in enumerate(rows, start=2):
        stage_counts[stage]["evaluated"] += 1
        if expected_lane_id is not None and row.get("lane_id") != expected_lane_id:
            fail(
                stage,
                index,
                row,
                "lane_id",
                "candidate row belongs to a different physical lane",
                observed=row.get("lane_id"),
                reference=expected_lane_id,
            )
            break
        if (
            expected_runtime_contract_sha256 is not None
            and row.get("runtime_contract_sha256")
            != expected_runtime_contract_sha256
        ):
            fail(
                stage,
                index,
                row,
                "runtime_contract_sha256",
                "candidate row is not bound to the exact runtime contract",
                observed=row.get("runtime_contract_sha256"),
                reference=expected_runtime_contract_sha256,
            )
            break
        identity = row.get("candidate_identity", "")
        if not identity:
            fail(stage, index, row, "candidate_identity", "missing stable candidate identity")
            break
        if identity in identities:
            fail(stage, index, row, "candidate_identity", "duplicate stable candidate identity",
                 observed=identity)
            break
        identities.add(identity)
        try:
            identity_match = as_int(row, "identity_match")
        except AuditFailure as exc:
            fail(stage, index, row, "identity_match", str(exc))
            break
        if row.get("match_status") != "matched" or identity_match != 1:
            fail(
                stage, index, row, "match_status", "candidate-population mismatch",
                observed=row.get("match_status"), reference="matched",
            )
            break

    if first is None:
        stage = "features"
        for index, row in enumerate(rows, start=2):
            stage_counts[stage]["evaluated"] += 1
            try:
                for feature in ORACLE_FEATURE_NAMES:
                    ppg = as_float(row, f"{feature}_ppg12")
                    delta = abs(as_float(row, f"{feature}_delta"))
                    tolerance = max(1.0e-6, 1.0e-5 * abs(ppg))
                    if delta > tolerance:
                        fail(
                            stage, index, row, feature, "feature divergence",
                            observed=as_float(row, f"{feature}_rj"), reference=ppg,
                            delta=delta, tolerance=tolerance,
                        )
                        break
            except AuditFailure as exc:
                fail(stage, index, row, "feature_schema", str(exc))
            if first is not None:
                break

    if first is None:
        stage = "scores"
        for index, row in enumerate(rows, start=2):
            stage_counts[stage]["evaluated"] += 1
            try:
                for field in (
                    "base_E_score_delta", "base_v3E_score_delta",
                    "selected_score_delta", "stored_minus_routed_score",
                ):
                    delta = abs(as_float(row, field))
                    if delta > 1.0e-6:
                        reason = (
                            "stored score divergence"
                            if field == "stored_minus_routed_score"
                            else "score divergence"
                        )
                        fail(stage, index, row, field, reason,
                             delta=delta, tolerance=1.0e-6)
                        break
                if first is None:
                    delta = abs(
                        as_float(row, "rj_stored_bdt_score")
                        - as_float(row, "ppg12_selected_bdt_score")
                    )
                    if delta > 1.0e-6:
                        fail(stage, index, row, "rj_stored_bdt_score",
                             "stored-to-oracle score divergence",
                             delta=delta, tolerance=1.0e-6)
            except AuditFailure as exc:
                fail(stage, index, row, "score_schema", str(exc))
            if first is not None:
                break

    if first is None:
        stage = "route"
        for index, row in enumerate(rows, start=2):
            stage_counts[stage]["evaluated"] += 1
            try:
                if as_int(row, "model_route_agree") != 1:
                    fail(stage, index, row, "model_route_agree", "model-route divergence",
                         observed=row.get("rj_selected_model"),
                         reference=row.get("ppg12_selected_model"))
                elif row.get("rj_inferred_stored_model") != row.get("ppg12_selected_model"):
                    fail(stage, index, row, "rj_inferred_stored_model",
                         "stored model-route divergence",
                         observed=row.get("rj_inferred_stored_model"),
                         reference=row.get("ppg12_selected_model"))
            except AuditFailure as exc:
                fail(stage, index, row, "route_schema", str(exc))
            if first is not None:
                break

    if first is None:
        stage = "tags"
        for index, row in enumerate(rows, start=2):
            stage_counts[stage]["evaluated"] += 1
            try:
                if row.get("ppg12_tag_evidence_source") != PRESERVED_EXECUTABLE_EVIDENCE:
                    fail(
                        stage,
                        index,
                        row,
                        "ppg12_tag_evidence_source",
                        "preserved-executable tag evidence is missing",
                        observed=row.get("ppg12_tag_evidence_source"),
                        reference=PRESERVED_EXECUTABLE_EVIDENCE,
                    )
                else:
                    for field in (
                        "signal_status_agree", "common_agree", "tag_agree",
                        "stored_common_agree", "stored_tag_agree",
                    ):
                        if as_int(row, field) != 1:
                            fail(stage, index, row, field, f"{field} divergence")
                            break
            except AuditFailure as exc:
                fail(stage, index, row, "tag_schema", str(exc))
            if first is not None:
                break

    if first is None:
        stage = "isolation_abcd"
        for index, row in enumerate(rows, start=2):
            stage_counts[stage]["evaluated"] += 1
            try:
                if (
                    row.get("ppg12_isolation_abcd_evidence_source")
                    != PRESERVED_EXECUTABLE_EVIDENCE
                ):
                    fail(
                        stage,
                        index,
                        row,
                        "ppg12_isolation_abcd_evidence_source",
                        "preserved-executable isolation/ABCD evidence is missing",
                        observed=row.get("ppg12_isolation_abcd_evidence_source"),
                        reference=PRESERVED_EXECUTABLE_EVIDENCE,
                    )
                for left, right in (
                    ("rj_is_iso", "ppg12_is_iso"),
                    ("rj_is_noniso", "ppg12_is_noniso"),
                ):
                    if first is not None:
                        break
                    if as_int(row, left) != as_int(row, right):
                        fail(stage, index, row, left, "isolation assignment divergence",
                             observed=row.get(left), reference=row.get(right))
                        break
                if first is None:
                    for left, right in (
                        ("rj_raw_eiso", "ppg12_raw_eiso"),
                        ("rj_corrected_eiso", "ppg12_corrected_eiso"),
                        ("rj_iso_threshold", "ppg12_iso_threshold"),
                        ("rj_noniso_threshold", "ppg12_noniso_threshold"),
                    ):
                        delta = abs(as_float(row, left) - as_float(row, right))
                        if delta > 1.0e-6:
                            fail(stage, index, row, left, "isolation-value divergence",
                                 delta=delta, tolerance=1.0e-6)
                            break
                if first is None:
                    for field in ("abcd_agree", "stored_abcd_agree"):
                        if as_int(row, field) != 1:
                            fail(stage, index, row, field, f"{field} divergence")
                            break
            except AuditFailure as exc:
                fail(stage, index, row, "isolation_schema", str(exc))
            if first is not None:
                break

    if first is None:
        stage = "truth_response_fills"
        for index, row in enumerate(rows, start=2):
            stage_counts[stage]["evaluated"] += 1
            try:
                if (
                    row.get("ppg12_truth_response_fill_evidence_source")
                    != PRESERVED_EXECUTABLE_EVIDENCE
                ):
                    fail(
                        stage,
                        index,
                        row,
                        "ppg12_truth_response_fill_evidence_source",
                        "preserved-executable truth/response/fill evidence is missing",
                        observed=row.get("ppg12_truth_response_fill_evidence_source"),
                        reference=PRESERVED_EXECUTABLE_EVIDENCE,
                    )
                is_signal = as_int(row, "ppg12_is_signal")
                if is_signal not in (0, 1):
                    fail(stage, index, row, "ppg12_is_signal", "non-binary signal status")
                for field in ("truth_class_agree",):
                    if first is not None:
                        break
                    if as_int(row, field) != 1:
                        fail(stage, index, row, field, "truth-class divergence",
                             observed=row.get("rj_truth_class"), reference=row.get("truth_class"))
                        break
                for left, right in (
                    ("rj_logical_abcd_region", "ppg12_logical_abcd_region"),
                    ("rj_analysis_window_pass", "ppg12_analysis_window_pass"),
                    ("rj_signal_fill_multiplicity", "ppg12_fill_multiplicity"),
                    ("rj_signal_fill_A", "ppg12_signal_fill_A"),
                    ("rj_signal_fill_B", "ppg12_signal_fill_B"),
                    ("rj_signal_fill_C", "ppg12_signal_fill_C"),
                    ("rj_signal_fill_D", "ppg12_signal_fill_D"),
                ):
                    if first is not None:
                        break
                    if as_int(row, left) != as_int(row, right):
                        fail(stage, index, row, left, "truth/response/fill divergence",
                             observed=row.get(left), reference=row.get(right))
                if first is None and is_signal == 1:
                    response_delta = abs(
                        as_float(row, "rj_response_Et")
                        - as_float(row, "ppg12_response_Et")
                    )
                    if response_delta > 1.0e-6:
                        fail(
                            stage,
                            index,
                            row,
                            "rj_response_Et",
                            "response-ET divergence",
                            observed=row.get("rj_response_Et"),
                            reference=row.get("ppg12_response_Et"),
                            delta=response_delta,
                            tolerance=1.0e-6,
                        )
                if first is None and is_signal == 1 and (
                    as_int(row, "rj_response_window_pass")
                    != as_int(row, "ppg12_response_window_pass")
                ):
                    fail(
                        stage,
                        index,
                        row,
                        "rj_response_window_pass",
                        "response-window divergence",
                        observed=row.get("rj_response_window_pass"),
                        reference=row.get("ppg12_response_window_pass"),
                    )
            except AuditFailure as exc:
                fail(stage, index, row, "truth_response_fill_schema", str(exc))
            if first is not None:
                break

    if first is None:
        stage = "weights"
        for index, row in enumerate(rows, start=2):
            stage_counts[stage]["evaluated"] += 1
            try:
                if row.get("ppg12_weight_evidence_source") != PRESERVED_EXECUTABLE_EVIDENCE:
                    fail(
                        stage,
                        index,
                        row,
                        "ppg12_weight_evidence_source",
                        "preserved-executable event-weight evidence is missing",
                        observed=row.get("ppg12_weight_evidence_source"),
                        reference=PRESERVED_EXECUTABLE_EVIDENCE,
                    )
                elif as_int(row, "rj_weight_lane_code") != 1:
                    fail(stage, index, row, "rj_weight_lane_code",
                         "paired photon oracle must use photon+jet lane code 1",
                         observed=row.get("rj_weight_lane_code"), reference=1)
                elif as_int(row, "rj_weight_component_code") != 1:
                    fail(stage, index, row, "rj_weight_component_code",
                         "paired SI oracle must use SI photon component code 1",
                         observed=row.get("rj_weight_component_code"), reference=1)
                else:
                    final = as_float(row, "rj_weight_final")
                    tolerance = max(1.0e-6, 1.0e-6 * abs(final))
                    factor_pairs = (
                        ("rj_weight_slice", "ppg12_weight_sample"),
                        ("rj_weight_mix", "ppg12_weight_mix"),
                        ("rj_weight_period", "ppg12_weight_lumi"),
                    )
                    for left, right in factor_pairs:
                        left_value = as_float(row, left)
                        right_value = as_float(row, right)
                        delta = abs(left_value - right_value)
                        pair_tolerance = max(1.0e-9, 1.0e-6 * abs(right_value))
                        if delta > pair_tolerance:
                            fail(
                                stage, index, row, left, "weight-factor divergence",
                                observed=left_value, reference=right_value,
                                delta=delta, tolerance=pair_tolerance,
                            )
                            break
                    ppg_cross = (
                        as_float(row, "ppg12_weight_sample")
                        * as_float(row, "ppg12_weight_mix")
                        * as_float(row, "ppg12_weight_lumi")
                    )
                    if first is None:
                        delta = abs(ppg_cross - as_float(row, "ppg12_weight_cross"))
                        if delta > max(1.0e-9, 1.0e-6 * abs(ppg_cross)):
                            fail(stage, index, row, "ppg12_weight_cross",
                                 "executable cross-section weight does not factorize")
                    ppg_vertex = (
                        as_float(row, "ppg12_weight_vertex")
                        * as_float(row, "ppg12_weight_truth_vertex")
                    )
                    if first is None:
                        delta = abs(as_float(row, "rj_weight_vertex") - ppg_vertex)
                        if delta > max(1.0e-9, 1.0e-6 * abs(ppg_vertex)):
                            fail(stage, index, row, "rj_weight_vertex",
                                 "vertex-weight divergence", delta=delta,
                                 tolerance=max(1.0e-9, 1.0e-6 * abs(ppg_vertex)))
                    trigger = as_float(row, "ppg12_weight_trigger")
                    if first is None and abs(trigger - 1.0) > 1.0e-9:
                        fail(stage, index, row, "ppg12_weight_trigger",
                             "simulation trigger factor must be unity",
                             observed=trigger, reference=1.0)
                    expected_event = ppg_cross * ppg_vertex * trigger
                    if first is None:
                        delta = abs(expected_event - as_float(row, "ppg12_weight_event"))
                        if delta > max(1.0e-9, 1.0e-6 * abs(expected_event)):
                            fail(stage, index, row, "ppg12_weight_event",
                                 "executable event weight does not factorize")
                    cross_delta = abs(final - as_float(row, "ppg12_weight_final"))
                    if first is None and cross_delta > tolerance:
                        fail(
                            stage,
                            index,
                            row,
                            "rj_weight_final",
                            "preserved-executable event-weight divergence",
                            observed=row.get("rj_weight_final"),
                            reference=row.get("ppg12_weight_final"),
                            delta=cross_delta,
                            tolerance=tolerance,
                        )
                    for field in ("rj_weight_product_delta", "rj_event_weight_delta"):
                        if first is not None:
                            break
                        delta = abs(as_float(row, field))
                        if delta > tolerance:
                            fail(stage, index, row, field, "weight-factor divergence",
                                 delta=delta, tolerance=tolerance)
                            break
            except AuditFailure as exc:
                fail(stage, index, row, "weight_schema", str(exc))
            if first is not None:
                break

    return {
        "schema_version": FIRST_DIVERGENCE_SCHEMA_VERSION,
        "candidate_report": str(path.resolve()),
        "status": "FAIL" if first is not None else "PASS",
        "row_count": len(rows),
        "stage_order": list(FIRST_DIVERGENCE_STAGE_ORDER),
        "stage_counts": stage_counts,
        "first_divergence": first,
    }


def validate_candidate_report(
    path: Path,
    first_divergence_path: Path | None = None,
    *,
    expected_lane_id: str | None = None,
    expected_runtime_contract_sha256: str | None = None,
) -> int:
    report = analyze_candidate_report(
        path,
        expected_lane_id=expected_lane_id,
        expected_runtime_contract_sha256=expected_runtime_contract_sha256,
    )
    output = first_divergence_path or path.with_name("first_divergence.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if report["status"] != "PASS":
        finding = report["first_divergence"] or {}
        raise AuditFailure(
            f"{finding.get('stage', 'unknown')} first divergence at CSV row "
            f"{finding.get('csv_row', '?')}: {finding.get('reason', 'unknown')} "
            f"field={finding.get('field', '?')}"
        )
    return int(report["row_count"])


def validate_apply_evidence(path: Path, paths: dict[str, Any]) -> None:
    data = load_json(path)
    if data.get("schema_version") != 1:
        raise AuditFailure("apply_BDT stage evidence schema_version must be 1")
    if data.get("stage") != "source_locked_canonical_apply_bdt_split":
        raise AuditFailure("unexpected apply_BDT stage identity")
    if data.get("revision") != EXPECTED_ESTIMATOR_REVISION:
        raise AuditFailure("apply_BDT evidence uses the wrong source revision")
    if data.get("shared_by_recoeff_baseline_and_trace") is not True:
        raise AuditFailure("apply_BDT output is not shared by both RecoEff runs")
    for key, contract_role in (
        ("macro", "apply_bdt"),
        ("canonical_config", "apply_config"),
        ("runtime_config", "apply_bdt_runtime_config"),
        ("input", "ppg_raw_root"),
        ("output", "ppg_scored_root"),
    ):
        item = data.get(key, {})
        if not isinstance(item, dict):
            raise AuditFailure(f"apply_BDT evidence lacks {key}")
        target = require_file(paths.get(contract_role, ""), contract_role)
        if Path(str(item.get("path", ""))).resolve() != target.resolve():
            raise AuditFailure(f"apply_BDT evidence {key} path differs from contract")
        if str(item.get("sha256", "")) != sha256(target):
            raise AuditFailure(f"apply_BDT evidence {key} hash differs from contract")
    required_branches = {
        "cluster_bdt_CLUSTERINFO_CEMC_base_E",
        "cluster_bdt_CLUSTERINFO_CEMC_base_v3E",
        "cluster_npb_score_CLUSTERINFO_CEMC",
    }
    if set(data.get("branches_required", [])) != required_branches:
        raise AuditFailure("apply_BDT evidence branch contract differs")
    assets = data.get("model_assets", {})
    if not isinstance(assets, dict) or set(assets) != set(EXPECTED_APPLY_MODEL_HASHES):
        raise AuditFailure("apply_BDT evidence lacks the exact 12 model roles")
    for role, expected_hash in EXPECTED_APPLY_MODEL_HASHES.items():
        item = assets.get(role, {})
        if not isinstance(item, dict) or item.get("sha256") != expected_hash:
            raise AuditFailure(f"apply_BDT evidence has wrong model hash for {role}")


def validate_aggregate_evidence(
    path: Path,
    *,
    mode: str,
    contract_path: Path,
    lane_id: str,
    contract_paths: dict[str, Any],
) -> None:
    data = load_json(path)
    if data.get("schema_version") != 1:
        raise AuditFailure("RecoEff aggregate evidence schema_version must be 1")
    if data.get("evidence_source") != "preserved_ppg12_executable_aggregate":
        raise AuditFailure("RecoEff aggregate is not preserved-executable evidence")
    if data.get("mode") != mode or data.get("status") != "PASS":
        raise AuditFailure(f"RecoEff aggregate evidence did not pass in {mode} mode")
    equivalence = data.get("root_equivalence", {})
    if not isinstance(equivalence, dict) or equivalence.get("pass") is not True:
        raise AuditFailure("instrumented and uninstrumented RecoEff ROOT files differ")
    if equivalence.get("efficiency_root_exact") is not True:
        raise AuditFailure("RecoEff efficiency ROOT payloads differ")
    if equivalence.get("response_root_exact") is not True:
        raise AuditFailure("RecoEff response ROOT payloads differ")
    if mode == "full":
        comparison = data.get("aggregate_comparison", {})
        if not isinstance(comparison, dict) or comparison.get("pass") is not True:
            raise AuditFailure("executable/RecoilJets weighted leakage cells differ")
        required_scopes = {"tags", "isolation_abcd", "truth_abcd_fills", "weights"}
        if set(data.get("certified_scopes", [])) != required_scopes:
            raise AuditFailure("full RecoEff aggregate does not certify required scopes")
        lane_identity = data.get("lane_identity", {})
        if not isinstance(lane_identity, dict):
            raise AuditFailure("full RecoEff aggregate lacks lane identity")
        if lane_identity.get("lane_id") != lane_id:
            raise AuditFailure("full RecoEff aggregate belongs to a different lane")
        if lane_identity.get("runtime_contract_sha256") != sha256(contract_path):
            raise AuditFailure("full RecoEff aggregate is bound to a stale runtime contract")
        provenance = data.get("provenance", {})
        if not isinstance(provenance, dict):
            raise AuditFailure("full RecoEff aggregate lacks provenance")
        for aggregate_role, contract_role in (
            ("runtime_contract", None),
            ("candidate_csv", "candidate_csv"),
            ("trace_csv", "ppg12_executable_trace"),
            ("response_trace_csv", "ppg12_executable_response_trace"),
        ):
            item = provenance.get(aggregate_role, {})
            if not isinstance(item, dict):
                raise AuditFailure(
                    f"full RecoEff aggregate lacks {aggregate_role} provenance"
                )
            expected_path = (
                contract_path
                if contract_role is None
                else require_file(
                    str(contract_paths.get(contract_role, "")), contract_role
                )
            )
            if Path(str(item.get("path", ""))).resolve() != expected_path.resolve():
                raise AuditFailure(
                    f"full RecoEff aggregate {aggregate_role} path differs from contract"
                )
            if item.get("sha256") != sha256(expected_path):
                raise AuditFailure(
                    f"full RecoEff aggregate {aggregate_role} hash differs from contract"
                )


def validate_postrun(contract_path: Path) -> int:
    contract = validate_contract(contract_path)
    paths = contract["paths"]
    lane_id = lane_id_from_contract(contract["lane"])
    contract_digest = sha256(contract_path)
    recoil_files = validate_recoil_manifest(
        Path(paths["recoil_runtime_manifest"]), EXPECTED_OFFLINE_MAIN
    )
    calo_calib = Path(paths["calo_calib"])
    rng = contract["rng"]
    validate_log(
        Path(paths["ppg_log"]), "ppg12", calo_calib, rng, recoil_files
    )
    validate_log(
        Path(paths["recoil_log"]), "recoiljets", calo_calib, rng, recoil_files
    )
    require_root_file(
        Path(paths["ppg_scored_root"]),
        "slimtree",
        (
            "cluster_bdt_CLUSTERINFO_CEMC_base_E",
            "cluster_bdt_CLUSTERINFO_CEMC_base_v3E",
            "cluster_npb_score_CLUSTERINFO_CEMC",
        ),
    )
    require_root_file(
        Path(paths["recoil_root"]),
        "AuAuPhotonIDTrainingTree",
        (
            "cluster_Et_score_input",
            "tight_bdt_score",
            "ppg12_common_pass",
            "ppg12_tight_tag",
            "ppg12_raw_eiso",
            "ppg12_reco_eiso",
            "ppg12_is_iso",
            "ppg12_is_noniso",
            "ppg12_response_Et",
            "ppg12_truth_class",
            "ppg12_logical_abcd_region",
            "ppg12_analysis_window_pass",
            "ppg12_response_window_pass",
            "ppg12_signal_fill_A",
            "ppg12_signal_fill_B",
            "ppg12_signal_fill_C",
            "ppg12_signal_fill_D",
            "ppg12_signal_fill_multiplicity",
            "ppg12_weight_lane_code",
            "ppg12_weight_component_code",
            "ppg12_weight_slice",
            "ppg12_weight_vertex",
            "ppg12_weight_mix",
            "ppg12_weight_period",
            "ppg12_weight_final",
        ),
    )
    validate_apply_evidence(Path(paths["apply_bdt_stage_evidence"]), paths)
    validate_aggregate_evidence(
        Path(paths["ppg12_recoeff_root_equivalence"]),
        mode="root_equivalence_only",
        contract_path=contract_path,
        lane_id=lane_id,
        contract_paths=paths,
    )
    validate_aggregate_evidence(
        Path(paths["ppg12_executable_aggregate"]),
        mode="full",
        contract_path=contract_path,
        lane_id=lane_id,
        contract_paths=paths,
    )
    candidate_path = Path(paths["candidate_csv"])
    divergence_path = Path(
        paths.get("first_divergence_json", candidate_path.with_name("first_divergence.json"))
    )
    return validate_candidate_report(
        candidate_path,
        divergence_path,
        expected_lane_id=lane_id,
        expected_runtime_contract_sha256=contract_digest,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=("preflight", "postrun", "rng-contract"))
    parser.add_argument("--contract", type=Path)
    parser.add_argument("--seed", type=int)
    parser.add_argument(
        "--skip-ldd",
        action="store_true",
        help="Unit-test only: skip live library dependency resolution.",
    )
    args = parser.parse_args()
    try:
        if args.phase == "rng-contract":
            if args.seed is not None:
                raise AuditFailure(
                    "historical rng-contract forbids a synthetic --seed"
                )
            rng = historical_rng_contract()
            print(
                rng["first_ph_seed"],
                rng["pedestal_seed"],
                rng["pedestal_sequence"],
                ",".join(str(value) for value in rng["ph_seed_sequence"]),
            )
        elif args.contract is None:
            raise AuditFailure(f"{args.phase} requires --contract")
        elif args.phase == "preflight":
            validate_contract(args.contract, run_ldd=not args.skip_ldd)
            print(f"PPG12_PAIRED_ORACLE_PREFLIGHT_PASS contract={args.contract}")
        else:
            count = validate_postrun(args.contract)
            print(
                "PPG12_PAIRED_ORACLE_POSTRUN_PASS "
                f"contract={args.contract} candidates={count}"
            )
    except AuditFailure as exc:
        print(f"PPG12_PAIRED_ORACLE_FAIL: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
