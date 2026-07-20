#!/usr/bin/env python3
"""Freeze and plan the 20 inclusive-jet same-event PPG12 oracle lanes.

This helper is intentionally orchestration-only.  It never starts ROOT,
Condor, a merge, or a current-pointer update.  ``freeze`` converts an explicit
source specification into a content-addressed manifest.  ``plan`` revalidates
every byte before emitting deterministic argv vectors for the foreground
paired-oracle driver.

The paired driver must explicitly advertise the ``inclusive`` capability with
this source comment before a plan is execution-ready::

    # PPG12_PAIRED_ORACLE_CAPABILITIES: photon,inclusive

Until the inclusive executable path is admitted, the planner still emits the
complete reviewable 20-lane contract but fails closed when
``--require-execution-ready`` is requested.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import sys
from pathlib import Path
from typing import Any, Iterable


SOURCE_SPEC_SCHEMA = "ppg12-inclusive-paired-source-spec/v1"
SOURCE_MANIFEST_SCHEMA = "ppg12-inclusive-paired-source-manifest/v1"
PLAN_SCHEMA = "ppg12-inclusive-paired-oracle-plan/v1"

SAMPLES = ("Jet8", "Jet12", "Jet20", "Jet30", "Jet40")
PERIODS = ("0mrad", "1p5mrad")
INTERACTIONS = ("SI", "DI")
ROWS_PER_LANE = 5

OWNERSHIP_WINDOWS: dict[str, dict[str, Any]] = {
    "Jet8": {"lower_gev": 9.0, "upper_gev": 14.0},
    "Jet12": {"lower_gev": 14.0, "upper_gev": 21.0},
    "Jet20": {"lower_gev": 21.0, "upper_gev": 32.0},
    "Jet30": {"lower_gev": 32.0, "upper_gev": 42.0},
    "Jet40": {"lower_gev": 42.0, "upper_gev": 100.0},
}
CANDIDATE_ET_CAPS_GEV: dict[str, float] = {
    "Jet8": 15.0,
    "Jet12": 23.0,
    "Jet20": 35.0,
    "Jet30": 45.0,
    "Jet40": 100.0,
}

DRIVER_ARG_KEYS = (
    "setup_script",
    "ppg_macro",
    "g4_full_list",
    "truthjet_full_list",
    "apply_bdt",
    "apply_config",
    "base_e_model",
    "base_v3e_model",
    "npb_model",
    "tower_mask",
    "recoil_runtime_manifest",
    "recoil_config",
)
LANE_KEYS = {
    "lane_id",
    "sample",
    "period",
    "interaction",
    "sample_key",
    "source_family",
    "driver_args",
}
CAPABILITY_RE = re.compile(
    r"^#\s*PPG12_PAIRED_ORACLE_CAPABILITIES:\s*([^\n]+)$", re.MULTILINE
)
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


class ContractError(RuntimeError):
    """Raised when a source or plan contract is not exact."""


def canonical_bytes(payload: Any) -> bytes:
    return (
        json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
        + "\n"
    ).encode("utf-8")


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ContractError(f"cannot read JSON {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ContractError(f"top-level JSON object required: {path}")
    return payload


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(json.dumps(payload, indent=2, sort_keys=True).encode("utf-8") + b"\n")


def expected_lane_id(sample: str, period: str, interaction: str) -> str:
    return f"inclusive:{sample.lower()}:{period}:{interaction.lower()}"


def expected_lane_ids() -> list[str]:
    return [
        expected_lane_id(sample, period, interaction)
        for period in PERIODS
        for sample in SAMPLES
        for interaction in INTERACTIONS
    ]


def expected_sample_key(sample: str, interaction: str) -> str:
    base = f"run28_jet{sample.removeprefix('Jet')}"
    return f"{base}_double" if interaction == "DI" else base


def expected_source_family(interaction: str) -> str:
    return "js_pp200_signal_dual" if interaction == "DI" else "js_pp200_signal"


def canonical_ownership_contract() -> dict[str, Any]:
    return {
        "variable": "max_truth_jet_pt_r04",
        "placement": "before_downstream_inclusive_physics_fills",
        "lower_bound": "inclusive",
        "upper_bound": "inclusive",
        "zero_truth_jet": "reject",
        "windows": OWNERSHIP_WINDOWS,
    }


def canonical_candidate_et_cap_contract() -> dict[str, Any]:
    return {
        "variable": "cluster_et",
        "placement": "before_downstream_inclusive_candidate_physics_fills",
        "upper_bound": "inclusive",
        "caps_gev": CANDIDATE_ET_CAPS_GEV,
    }


def canonical_global_contract() -> dict[str, Any]:
    return {
        "rows_per_lane": ROWS_PER_LANE,
        "jet5": "excluded",
        "global_scale": 1.0,
        "jet8_scale": 1.0,
        "external_scale": 1.0,
        "ownership_gate": canonical_ownership_contract(),
        "candidate_et_cap": canonical_candidate_et_cap_contract(),
        "execution": "foreground_only_no_scheduler_no_merge_no_promotion",
    }


def require_unit_scale(name: str, value: Any) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ContractError(f"{name} must be numeric 1.0")
    if not math.isfinite(float(value)) or float(value) != 1.0:
        raise ContractError(f"{name} must be exactly 1.0; observed {value!r}")


def require_exact_contract(contract: Any) -> None:
    if not isinstance(contract, dict):
        raise ContractError("contract must be an object")
    expected = canonical_global_contract()
    if set(contract) != set(expected):
        raise ContractError(
            f"contract keys differ: expected {sorted(expected)}, observed {sorted(contract)}"
        )
    for key in ("global_scale", "jet8_scale", "external_scale"):
        require_unit_scale(key, contract.get(key))
    for key in (
        "rows_per_lane",
        "jet5",
        "ownership_gate",
        "candidate_et_cap",
        "execution",
    ):
        if contract.get(key) != expected[key]:
            raise ContractError(
                f"contract {key} differs from the frozen inclusive contract"
            )


def require_absolute_file(path_text: Any, role: str) -> Path:
    if not isinstance(path_text, str) or not path_text:
        raise ContractError(f"{role} must be a nonempty path string")
    if "\n" in path_text or "\r" in path_text:
        raise ContractError(f"{role} contains a newline")
    path = Path(path_text)
    if not path.is_absolute():
        raise ContractError(f"{role} must be absolute: {path}")
    if not path.is_file() or path.stat().st_size <= 0:
        raise ContractError(f"{role} is missing or empty: {path}")
    return path


def noncomment_rows(path: Path) -> list[str]:
    rows = [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if len(rows) < ROWS_PER_LANE:
        raise ContractError(
            f"source list {path} has {len(rows)} rows; need at least {ROWS_PER_LANE}"
        )
    return rows


def normalized_event_identity(row: str, kind: str, label: str) -> str:
    prefix = {"g4": "G4Hits_", "truthjet": "DST_TRUTH_JET_"}[kind]
    if any(character.isspace() for character in row):
        raise ContractError(f"{label} contains whitespace: {row!r}")
    path = Path(row)
    basename = path.name
    if not path.is_absolute() and str(path) != basename:
        raise ContractError(
            f"{label} must be an absolute path or preserved basename: {row}"
        )
    if not basename.startswith(prefix) or not basename.endswith(".root"):
        raise ContractError(
            f"{label} must use the recognized {prefix}*.root naming contract: {row}"
        )
    identity = basename.removeprefix(prefix)
    if not identity or identity == ".root":
        raise ContractError(f"{label} has an empty event identity: {row}")
    return identity


def active_cpp(text: str) -> str:
    without_blocks = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
    return re.sub(r"//[^\n]*", "", without_blocks)


def validate_ppg_macro(path: Path, interaction: str, label: str) -> dict[str, Any]:
    active = active_cpp(path.read_text(encoding="utf-8"))
    active_indices = sorted(
        {
            int(value)
            for value in re.findall(
                r"INPUTREADHITS\s*::\s*listfile\s*\[\s*([0-9]+)\s*\]\s*=",
                active,
            )
        }
    )
    if active_indices != [0, 4]:
        raise ContractError(
            f"{label}: active INPUTREADHITS list indices must be exactly {{0, 4}}, "
            f"observed={active_indices}"
        )
    embedding_flag_two = bool(
        re.search(r"\badd_embedding_flag\s*\(\s*2\s*\)", active)
    )
    truth_jet_input = bool(re.search(r"\bnew\s+TruthJetInput\b", active))
    if interaction == "DI":
        if not truth_jet_input:
            raise ContractError(
                f"{label}: DI macro must actively reconstruct truth jets with TruthJetInput"
            )
        if not embedding_flag_two:
            raise ContractError(
                f"{label}: DI macro must actively call add_embedding_flag(2)"
            )
    elif embedding_flag_two:
        raise ContractError(
            f"{label}: SI macro must not actively call add_embedding_flag(2)"
        )
    return {
        "active_input_indices": active_indices,
        "truth_jet_input": truth_jet_input,
        "embedding_flag_2": embedding_flag_two,
    }


def file_record(path: Path) -> dict[str, Any]:
    return {
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def validate_source_rows(
    *, sample: str, interaction: str, g4_rows: list[str], truth_rows: list[str]
) -> list[str]:
    if len(g4_rows) != len(truth_rows):
        raise ContractError(
            f"{sample} {interaction} G4/truth list lengths differ: "
            f"{len(g4_rows)} vs {len(truth_rows)}"
        )
    family = expected_source_family(interaction)
    family_fragment = f"/{family}/"
    bad_g4 = next((row for row in g4_rows if family_fragment not in row), None)
    bad_truth = next((row for row in truth_rows if family_fragment not in row), None)
    if bad_g4 is not None or bad_truth is not None:
        raise ContractError(
            f"{sample} {interaction} source row does not belong to {family}"
        )
    jet = sample.removeprefix("Jet")
    truth_slice_re = re.compile(rf"/(?:run0028|run28)/jet{re.escape(jet)}/")
    if any(not truth_slice_re.search(row) for row in truth_rows):
        raise ContractError(
            f"{sample} {interaction} truth-jet list is not uniformly the jet{jet} slice"
        )
    g4_identities = [
        normalized_event_identity(
            row, "g4", f"{sample} {interaction} G4 row {index}"
        )
        for index, row in enumerate(g4_rows, start=1)
    ]
    truth_identities = [
        normalized_event_identity(
            row, "truthjet", f"{sample} {interaction} TRUTH_JET row {index}"
        )
        for index, row in enumerate(truth_rows, start=1)
    ]
    if len(g4_identities) != len(set(g4_identities)):
        raise ContractError(
            f"{sample} {interaction} normalized G4 event identities are not unique"
        )
    if len(truth_identities) != len(set(truth_identities)):
        raise ContractError(
            f"{sample} {interaction} normalized TRUTH_JET event identities are not unique"
        )
    for index, (g4_identity, truth_identity) in enumerate(
        zip(g4_identities, truth_identities), start=1
    ):
        if g4_identity != truth_identity:
            raise ContractError(
                f"{sample} {interaction} event identity mismatch at row {index}: "
                f"G4={g4_identity!r}, TRUTH_JET={truth_identity!r}"
            )
    return g4_identities


def freeze_lane(lane: Any) -> dict[str, Any]:
    if not isinstance(lane, dict):
        raise ContractError("each lane must be an object")
    if set(lane) != LANE_KEYS:
        raise ContractError(
            f"lane keys differ: expected {sorted(LANE_KEYS)}, observed {sorted(lane)}"
        )
    sample = lane["sample"]
    period = lane["period"]
    interaction = lane["interaction"]
    if sample not in SAMPLES:
        raise ContractError(f"unsupported inclusive sample: {sample!r}")
    if period not in PERIODS:
        raise ContractError(f"unsupported period: {period!r}")
    if interaction not in INTERACTIONS:
        raise ContractError(f"unsupported interaction: {interaction!r}")
    lane_id = expected_lane_id(sample, period, interaction)
    if lane["lane_id"] != lane_id:
        raise ContractError(f"lane_id mismatch: expected {lane_id}, got {lane['lane_id']}")
    sample_key = expected_sample_key(sample, interaction)
    if lane["sample_key"] != sample_key:
        raise ContractError(
            f"sample_key mismatch for {lane_id}: expected {sample_key}, got {lane['sample_key']}"
        )
    family = expected_source_family(interaction)
    if lane["source_family"] != family:
        raise ContractError(
            f"source family mismatch for {lane_id}: expected {family}, got {lane['source_family']}"
        )

    driver_args = lane["driver_args"]
    if not isinstance(driver_args, dict) or set(driver_args) != set(DRIVER_ARG_KEYS):
        observed = sorted(driver_args) if isinstance(driver_args, dict) else type(driver_args).__name__
        raise ContractError(
            f"driver_args for {lane_id} must contain exactly {list(DRIVER_ARG_KEYS)}; "
            f"observed {observed}"
        )
    paths = {
        key: require_absolute_file(driver_args[key], f"{lane_id}:{key}")
        for key in DRIVER_ARG_KEYS
    }
    g4_rows = noncomment_rows(paths["g4_full_list"])
    truth_rows = noncomment_rows(paths["truthjet_full_list"])
    event_identities = validate_source_rows(
        sample=sample,
        interaction=interaction,
        g4_rows=g4_rows,
        truth_rows=truth_rows,
    )
    macro_contract = validate_ppg_macro(
        paths["ppg_macro"], interaction, f"{lane_id}:ppg_macro"
    )
    first_five_pairs = [
        {"row": index + 1, "g4": g4, "truthjet": truth}
        for index, (g4, truth) in enumerate(
            zip(g4_rows[:ROWS_PER_LANE], truth_rows[:ROWS_PER_LANE])
        )
    ]
    first_five_identities = event_identities[:ROWS_PER_LANE]
    return {
        "lane_id": lane_id,
        "sample": sample,
        "period": period,
        "interaction": interaction,
        "sample_key": sample_key,
        "source_family": family,
        "ownership_window": OWNERSHIP_WINDOWS[sample],
        "candidate_et_cap_gev": CANDIDATE_ET_CAPS_GEV[sample],
        "driver_args": {key: file_record(paths[key]) for key in DRIVER_ARG_KEYS},
        "ppg_macro_contract": macro_contract,
        "source_binding": {
            "row_count": len(g4_rows),
            "first_five_pairs": first_five_pairs,
            "first_five_pairs_sha256": sha256_bytes(canonical_bytes(first_five_pairs)),
            "first_five_event_identities": first_five_identities,
            "first_five_event_identity_sha256": sha256_bytes(
                "\n".join(first_five_identities).encode("utf-8")
            ),
            "all_event_identity_sha256": sha256_bytes(
                "\n".join(event_identities).encode("utf-8")
            ),
        },
    }


def freeze_spec(spec_path: Path) -> dict[str, Any]:
    spec = load_json(spec_path)
    if spec.get("schema") != SOURCE_SPEC_SCHEMA:
        raise ContractError(
            f"source spec schema must be {SOURCE_SPEC_SCHEMA}; observed {spec.get('schema')!r}"
        )
    if set(spec) != {"schema", "contract", "lanes"}:
        raise ContractError("source spec must contain exactly schema, contract, and lanes")
    require_exact_contract(spec["contract"])
    if not isinstance(spec["lanes"], list):
        raise ContractError("lanes must be a list")
    frozen = [freeze_lane(lane) for lane in spec["lanes"]]
    by_id: dict[str, dict[str, Any]] = {}
    for lane in frozen:
        lane_id = lane["lane_id"]
        if lane_id in by_id:
            raise ContractError(f"duplicate lane: {lane_id}")
        by_id[lane_id] = lane
    expected = expected_lane_ids()
    if set(by_id) != set(expected):
        missing = sorted(set(expected) - set(by_id))
        extra = sorted(set(by_id) - set(expected))
        raise ContractError(f"20-lane coverage mismatch; missing={missing}, extra={extra}")
    payload = {
        "schema": SOURCE_MANIFEST_SCHEMA,
        "status": "source_locked",
        "contract": canonical_global_contract(),
        "lanes": [by_id[lane_id] for lane_id in expected],
    }
    payload["manifest_payload_sha256"] = sha256_bytes(canonical_bytes(payload))
    return payload


def verify_file_record(record: Any, role: str) -> Path:
    if not isinstance(record, dict) or set(record) != {"path", "bytes", "sha256"}:
        raise ContractError(f"invalid file record for {role}")
    if not isinstance(record["sha256"], str) or not SHA256_RE.fullmatch(record["sha256"]):
        raise ContractError(f"invalid sha256 for {role}")
    path = require_absolute_file(record["path"], role)
    if path.stat().st_size != record["bytes"]:
        raise ContractError(f"byte size changed for {role}: {path}")
    observed = sha256_file(path)
    if observed != record["sha256"]:
        raise ContractError(f"sha256 changed for {role}: {path}")
    return path


def verify_manifest(manifest_path: Path) -> dict[str, Any]:
    manifest = load_json(manifest_path)
    if manifest.get("schema") != SOURCE_MANIFEST_SCHEMA:
        raise ContractError(f"manifest schema must be {SOURCE_MANIFEST_SCHEMA}")
    if set(manifest) != {
        "schema",
        "status",
        "contract",
        "lanes",
        "manifest_payload_sha256",
    }:
        raise ContractError("source manifest contains unexpected or missing top-level keys")
    if manifest["status"] != "source_locked":
        raise ContractError("source manifest is not source_locked")
    require_exact_contract(manifest["contract"])
    expected_digest = manifest.pop("manifest_payload_sha256")
    if not isinstance(expected_digest, str) or not SHA256_RE.fullmatch(expected_digest):
        raise ContractError("source manifest payload hash is malformed")
    observed_digest = sha256_bytes(canonical_bytes(manifest))
    manifest["manifest_payload_sha256"] = expected_digest
    if observed_digest != expected_digest:
        raise ContractError("source manifest payload hash changed")

    lanes = manifest["lanes"]
    if not isinstance(lanes, list) or len(lanes) != 20:
        raise ContractError("source manifest must contain exactly 20 lanes")
    observed_ids = [lane.get("lane_id") for lane in lanes if isinstance(lane, dict)]
    if observed_ids != expected_lane_ids():
        raise ContractError("source manifest lane order or coverage differs from canonical order")
    for lane in lanes:
        lane_id = lane["lane_id"]
        expected_fields = {
            "lane_id",
            "sample",
            "period",
            "interaction",
            "sample_key",
            "source_family",
            "ownership_window",
            "candidate_et_cap_gev",
            "driver_args",
            "ppg_macro_contract",
            "source_binding",
        }
        if set(lane) != expected_fields:
            raise ContractError(f"manifest lane fields changed for {lane_id}")
        sample = lane["sample"]
        period = lane["period"]
        interaction = lane["interaction"]
        if lane_id != expected_lane_id(sample, period, interaction):
            raise ContractError(f"manifest lane identity changed for {lane_id}")
        if lane["sample_key"] != expected_sample_key(sample, interaction):
            raise ContractError(f"manifest sample key changed for {lane_id}")
        if lane["source_family"] != expected_source_family(interaction):
            raise ContractError(f"manifest source family changed for {lane_id}")
        if lane["ownership_window"] != OWNERSHIP_WINDOWS[sample]:
            raise ContractError(f"manifest ownership window changed for {lane_id}")
        if lane["candidate_et_cap_gev"] != CANDIDATE_ET_CAPS_GEV[sample]:
            raise ContractError(f"manifest candidate ET cap changed for {lane_id}")
        args = lane["driver_args"]
        if not isinstance(args, dict) or set(args) != set(DRIVER_ARG_KEYS):
            raise ContractError(f"manifest driver args changed for {lane_id}")
        paths = {
            key: verify_file_record(args[key], f"{lane_id}:{key}")
            for key in DRIVER_ARG_KEYS
        }
        g4_rows = noncomment_rows(paths["g4_full_list"])
        truth_rows = noncomment_rows(paths["truthjet_full_list"])
        event_identities = validate_source_rows(
            sample=sample,
            interaction=interaction,
            g4_rows=g4_rows,
            truth_rows=truth_rows,
        )
        macro_contract = validate_ppg_macro(
            paths["ppg_macro"], interaction, f"{lane_id}:ppg_macro"
        )
        if lane["ppg_macro_contract"] != macro_contract:
            raise ContractError(f"PPG macro contract changed for {lane_id}")
        binding = lane["source_binding"]
        if not isinstance(binding, dict) or set(binding) != {
            "row_count",
            "first_five_pairs",
            "first_five_pairs_sha256",
            "first_five_event_identities",
            "first_five_event_identity_sha256",
            "all_event_identity_sha256",
        }:
            raise ContractError(f"manifest source binding changed for {lane_id}")
        pairs = [
            {"row": index + 1, "g4": g4, "truthjet": truth}
            for index, (g4, truth) in enumerate(
                zip(g4_rows[:ROWS_PER_LANE], truth_rows[:ROWS_PER_LANE])
            )
        ]
        if binding["row_count"] != len(g4_rows) or binding["first_five_pairs"] != pairs:
            raise ContractError(f"first-five source binding changed for {lane_id}")
        pair_digest = sha256_bytes(canonical_bytes(pairs))
        if binding["first_five_pairs_sha256"] != pair_digest:
            raise ContractError(f"first-five pair hash changed for {lane_id}")
        first_five_identities = event_identities[:ROWS_PER_LANE]
        if binding["first_five_event_identities"] != first_five_identities:
            raise ContractError(f"first-five event identities changed for {lane_id}")
        identity_digest = sha256_bytes(
            "\n".join(first_five_identities).encode("utf-8")
        )
        if binding["first_five_event_identity_sha256"] != identity_digest:
            raise ContractError(f"first-five event identity hash changed for {lane_id}")
        all_identity_digest = sha256_bytes(
            "\n".join(event_identities).encode("utf-8")
        )
        if binding["all_event_identity_sha256"] != all_identity_digest:
            raise ContractError(f"all-row event identity hash changed for {lane_id}")
    return manifest


def driver_capabilities(driver: Path) -> set[str]:
    text = driver.read_text(encoding="utf-8")
    match = CAPABILITY_RE.search(text)
    if not match:
        return set()
    return {item.strip().lower() for item in match.group(1).split(",") if item.strip()}


def option_for_key(key: str) -> str:
    return "--" + key.replace("_", "-")


def build_plan(
    *, manifest_path: Path, driver_path: Path, output_root: Path
) -> dict[str, Any]:
    manifest = verify_manifest(manifest_path)
    driver = require_absolute_file(str(driver_path), "paired oracle driver")
    if not os.access(driver, os.X_OK):
        raise ContractError(f"paired oracle driver is not executable: {driver}")
    if not output_root.is_absolute() or str(output_root) == "/":
        raise ContractError(f"output root must be an absolute non-root path: {output_root}")
    if output_root.exists():
        raise ContractError(f"output root must be absent at plan time: {output_root}")

    capabilities = sorted(driver_capabilities(driver))
    blockers: list[str] = []
    if "inclusive" not in capabilities:
        blockers.append(
            "paired driver has not advertised the source-reviewed inclusive executable capability"
        )
    lanes: list[dict[str, Any]] = []
    for lane in manifest["lanes"]:
        lane_output = output_root / lane["lane_id"].replace(":", "__")
        argv = [
            str(driver),
            "--lane-id",
            lane["lane_id"],
            "--sample",
            lane["sample"],
            "--period",
            lane["period"],
            "--interaction",
            lane["interaction"],
            "--output-dir",
            str(lane_output),
        ]
        for key in DRIVER_ARG_KEYS:
            argv.extend([option_for_key(key), lane["driver_args"][key]["path"]])
        if "--run" in argv or "--token" in argv:
            raise ContractError("planner generated a mutation-bearing driver argument")
        lane_payload = {
            "lane_id": lane["lane_id"],
            "sample_key": lane["sample_key"],
            "source_family": lane["source_family"],
            "ownership_window": lane["ownership_window"],
            "candidate_et_cap_gev": lane["candidate_et_cap_gev"],
            "first_five_pairs_sha256": lane["source_binding"][
                "first_five_pairs_sha256"
            ],
            "first_five_event_identity_sha256": lane["source_binding"][
                "first_five_event_identity_sha256"
            ],
            "ppg_macro_contract": lane["ppg_macro_contract"],
            "output_dir": str(lane_output),
            "argv": argv,
        }
        lane_payload["lane_contract_sha256"] = sha256_bytes(
            canonical_bytes(lane_payload)
        )
        lanes.append(lane_payload)

    source_manifest_record = file_record(manifest_path.resolve())
    payload = {
        "schema": PLAN_SCHEMA,
        "status": "execution_ready" if not blockers else "plan_only_blocked",
        "execution_ready": not blockers,
        "blockers": blockers,
        "execution_contract": "foreground_only_no_scheduler_no_merge_no_promotion",
        "source_manifest": {
            **source_manifest_record,
            "manifest_payload_sha256": manifest["manifest_payload_sha256"],
        },
        "driver": {**file_record(driver), "capabilities": capabilities},
        "output_root": str(output_root),
        "lane_count": len(lanes),
        "rows_per_lane": ROWS_PER_LANE,
        "jet5": "excluded",
        "scales": {"global": 1.0, "jet8": 1.0, "external": 1.0},
        "lanes": lanes,
    }
    payload["plan_payload_sha256"] = sha256_bytes(canonical_bytes(payload))
    return payload


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    freeze = subparsers.add_parser("freeze", help="freeze a 20-lane source spec")
    freeze.add_argument("--spec", type=Path, required=True)
    freeze.add_argument("--out", type=Path, required=True)

    verify = subparsers.add_parser("verify", help="revalidate a frozen source manifest")
    verify.add_argument("--manifest", type=Path, required=True)

    plan = subparsers.add_parser("plan", help="emit a deterministic foreground plan")
    plan.add_argument("--manifest", type=Path, required=True)
    plan.add_argument("--driver", type=Path, required=True)
    plan.add_argument("--output-root", type=Path, required=True)
    plan.add_argument("--out", type=Path)
    plan.add_argument("--require-execution-ready", action="store_true")
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        if args.command == "freeze":
            payload = freeze_spec(args.spec.resolve())
            write_json(args.out.resolve(), payload)
            print(
                "PPG12_INCLUSIVE_SOURCE_FREEZE_PASS "
                f"lanes={len(payload['lanes'])} manifest={args.out.resolve()} "
                f"sha256={payload['manifest_payload_sha256']}"
            )
            return 0
        if args.command == "verify":
            payload = verify_manifest(args.manifest.resolve())
            print(
                "PPG12_INCLUSIVE_SOURCE_VERIFY_PASS "
                f"lanes={len(payload['lanes'])} "
                f"sha256={payload['manifest_payload_sha256']}"
            )
            return 0
        if args.command == "plan":
            payload = build_plan(
                manifest_path=args.manifest.resolve(),
                driver_path=args.driver.resolve(),
                output_root=args.output_root,
            )
            if args.require_execution_ready and not payload["execution_ready"]:
                raise ContractError("; ".join(payload["blockers"]))
            if args.out:
                write_json(args.out.resolve(), payload)
            else:
                sys.stdout.buffer.write(json.dumps(payload, indent=2, sort_keys=True).encode("utf-8") + b"\n")
            print(
                "PPG12_INCLUSIVE_PLAN_PASS "
                f"lanes={payload['lane_count']} status={payload['status']} "
                f"sha256={payload['plan_payload_sha256']}",
                file=sys.stderr,
            )
            return 0
    except ContractError as exc:
        print(f"PPG12_INCLUSIVE_ORCHESTRATION_FAIL: {exc}", file=sys.stderr)
        return 2
    raise AssertionError(f"unhandled command: {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
