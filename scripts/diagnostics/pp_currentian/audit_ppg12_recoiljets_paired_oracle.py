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


SCHEMA_VERSION = 2
EXPECTED_LANE = {
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
}


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
        raise AuditFailure("paired-oracle contract schema_version must be 2")
    lane = contract.get("lane")
    if not isinstance(lane, dict):
        raise AuditFailure("paired-oracle lane contract must be an object")
    for key, expected in EXPECTED_LANE.items():
        if lane.get(key) != expected:
            raise AuditFailure(f"unsupported lane contract: {lane}")
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
    ppg_macro = require_file(paths.get("ppg_macro", ""), "preserved PPG12 macro")
    require_hash(ppg_macro, EXPECTED_HASHES["ppg_macro"], "preserved PPG12 macro")
    validate_frozen_macro_graph(ppg_macro)
    g4_full = require_file(paths.get("g4_full_list", ""), "G4 full list")
    truth_full = require_file(paths.get("truthjet_full_list", ""), "truth-jet full list")
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


def validate_candidate_report(path: Path) -> int:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise AuditFailure("candidate report contains no rows")
    feature_names = (
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
    for index, row in enumerate(rows, start=2):
        if row.get("match_status") != "matched" or as_int(row, "identity_match") != 1:
            raise AuditFailure(f"candidate-population mismatch at CSV row {index}")
        for field in ("model_route_agree", "common_agree", "tag_agree", "abcd_agree"):
            if as_int(row, field) != 1:
                raise AuditFailure(f"{field} divergence at CSV row {index}")
        if row.get("rj_inferred_stored_model") != row.get("ppg12_selected_model"):
            raise AuditFailure(f"stored model-route divergence at CSV row {index}")
        for feature in feature_names:
            ppg = as_float(row, f"{feature}_ppg12")
            delta = abs(as_float(row, f"{feature}_delta"))
            tolerance = max(1.0e-6, 1.0e-5 * abs(ppg))
            if delta > tolerance:
                raise AuditFailure(
                    f"feature divergence {feature} at row {index}: {delta} > {tolerance}"
                )
        for field in ("base_E_score_delta", "base_v3E_score_delta", "selected_score_delta"):
            if abs(as_float(row, field)) > 1.0e-6:
                raise AuditFailure(f"score divergence {field} at CSV row {index}")
        if abs(as_float(row, "stored_minus_routed_score")) > 1.0e-6:
            raise AuditFailure(f"stored score divergence at CSV row {index}")
        if abs(
            as_float(row, "rj_stored_bdt_score")
            - as_float(row, "ppg12_selected_bdt_score")
        ) > 1.0e-6:
            raise AuditFailure(f"stored-to-oracle score divergence at CSV row {index}")
        for field in ("stored_common_agree", "stored_tag_agree", "stored_abcd_agree"):
            if as_int(row, field) != 1:
                raise AuditFailure(f"{field} divergence at CSV row {index}")
        for left, right in (
            ("rj_is_iso", "ppg12_is_iso"),
            ("rj_is_noniso", "ppg12_is_noniso"),
        ):
            if as_int(row, left) != as_int(row, right):
                raise AuditFailure(f"isolation assignment divergence at CSV row {index}")
        for left, right in (
            ("rj_raw_eiso", "ppg12_raw_eiso"),
            ("rj_corrected_eiso", "ppg12_corrected_eiso"),
            ("rj_iso_threshold", "ppg12_iso_threshold"),
            ("rj_noniso_threshold", "ppg12_noniso_threshold"),
        ):
            if abs(as_float(row, left) - as_float(row, right)) > 1.0e-6:
                raise AuditFailure(f"isolation-value divergence at CSV row {index}")
    return len(rows)


def validate_postrun(contract_path: Path) -> int:
    contract = validate_contract(contract_path)
    paths = contract["paths"]
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
        ),
    )
    return validate_candidate_report(Path(paths["candidate_csv"]))


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
