#!/usr/bin/env python3
"""Fail-closed 12-lane Condor fanout for the PPG12 paired executable oracle.

This helper is invoked by ``submit_ppg12_stitched_purity_photon_canaries.sh``.
It generates one deterministic plan, queue, wrapper, and submit description;
submits exactly one job per canonical photon lane only after a matching token;
runs one token-bound lane inside Condor; and provides a read-only JSON audit.

There is deliberately no merge, promotion, retry, release, removal, or other
job-control interface here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any


PLAN_SCHEMA = "ppg12-paired-oracle-condor-plan/v1"
GENERATION_SCHEMA = "ppg12-paired-oracle-condor-generation/v1"
RECEIPT_SCHEMA = "ppg12-paired-oracle-condor-lane-receipt/v1"
AUDIT_SCHEMA = "ppg12-paired-oracle-condor-audit/v1"
SUBMISSION_SCHEMA = "ppg12-paired-oracle-condor-submission/v1"
SOURCE_SCHEMA = "ppg12-paired-source-manifest/v1"
TOKEN_PREFIX = "RUN_PPG12_PAIRED_CONDOR_"
TOKEN_RE = re.compile(r"^RUN_PPG12_PAIRED_CONDOR_[0-9a-f]{64}$")
LANE_TOKEN_RE = re.compile(r"^ppg12-oracle:[0-9a-f]{64}$")
SAFE_TAG_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,159}$")
SAFE_CHAT_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9 ._:+/@|()-]{0,255}$")
SAFE_THREAD_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:+/@-]{0,255}$")
COMMON_NAMES = (
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
LANE_SOURCE_NAMES = ("ppg_macro", "g4_full_list", "truthjet_full_list")


class FanoutError(RuntimeError):
    """The fanout contract is incomplete, stale, or unsafe."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def payload_sha256(payload: object) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def read_json(path: Path, label: str) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise FanoutError(f"missing {label}: {path}") from exc
    except json.JSONDecodeError as exc:
        raise FanoutError(f"invalid {label}: {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise FanoutError(f"{label} must contain one JSON object: {path}")
    return payload


def absolute_file(value: object, label: str, *, executable: bool = False) -> Path:
    if not isinstance(value, (str, Path)) or not str(value).strip():
        raise FanoutError(f"{label} must be a non-empty absolute path")
    raw_text = str(value)
    if raw_text != raw_text.strip() or any(char in raw_text for char in "\r\n\t"):
        raise FanoutError(f"{label} contains unsafe whitespace: {raw_text!r}")
    raw = Path(raw_text)
    if not raw.is_absolute():
        raise FanoutError(f"{label} must be absolute: {raw}")
    try:
        path = raw.resolve(strict=True)
    except FileNotFoundError as exc:
        raise FanoutError(f"missing {label}: {raw}") from exc
    if not path.is_file() or path.stat().st_size <= 0:
        raise FanoutError(f"{label} must be a non-empty file: {path}")
    if executable and not os.access(path, os.X_OK):
        raise FanoutError(f"{label} must be executable: {path}")
    return path


def absolute_directory_string(value: object, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise FanoutError(f"{label} must be a non-empty absolute path")
    if value != value.strip() or any(char.isspace() for char in value):
        raise FanoutError(f"{label} cannot contain whitespace in Condor mode: {value!r}")
    if not Path(value).is_absolute():
        raise FanoutError(f"{label} must be absolute: {value}")
    return str(Path(value))


def file_link(path: Path) -> dict[str, str]:
    return {"path": str(path), "sha256": sha256(path)}


def validate_link(value: object, label: str, *, executable: bool = False) -> dict[str, str]:
    if not isinstance(value, dict) or set(value) != {"path", "sha256"}:
        raise FanoutError(f"{label} must be an exact path/SHA-256 link")
    path = absolute_file(value["path"], label, executable=executable)
    observed = sha256(path)
    if value["sha256"] != observed:
        raise FanoutError(
            f"{label} SHA-256 drift: expected={value['sha256']}, observed={observed}"
        )
    return {"path": str(path), "sha256": observed}


def expected_lanes() -> list[tuple[str, str, str, str]]:
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


def _list_rows(path: Path, label: str) -> list[str]:
    rows = [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if len(rows) < 5:
        raise FanoutError(f"{label} must contain at least five source rows")
    if len(rows) != len(set(rows)):
        raise FanoutError(f"{label} contains duplicate source rows")
    return rows


def _event_identity(row: str, kind: str, label: str) -> str:
    prefix = {"g4": "G4Hits_", "truthjet": "DST_TRUTH_JET_"}[kind]
    if any(char.isspace() for char in row):
        raise FanoutError(f"{label} contains whitespace: {row!r}")
    path = Path(row)
    basename = path.name
    if not path.is_absolute() and str(path) != basename:
        raise FanoutError(
            f"{label} must be an absolute path or preserved basename: {row}"
        )
    if not basename.startswith(prefix) or not basename.endswith(".root"):
        raise FanoutError(f"{label} must match recognized {prefix}*.root")
    identity = basename.removeprefix(prefix)
    if not identity or identity == ".root":
        raise FanoutError(f"{label} has an empty event identity")
    return identity


def _paired_rows(g4_rows: list[str], truth_rows: list[str], label: str) -> str:
    if len(g4_rows) != len(truth_rows):
        raise FanoutError(f"{label}: G4 and TRUTH_JET list lengths differ")
    g4_ids = [
        _event_identity(row, "g4", f"{label} G4 row {index}")
        for index, row in enumerate(g4_rows, start=1)
    ]
    truth_ids = [
        _event_identity(row, "truthjet", f"{label} TRUTH_JET row {index}")
        for index, row in enumerate(truth_rows, start=1)
    ]
    if len(g4_ids) != len(set(g4_ids)) or len(truth_ids) != len(set(truth_ids)):
        raise FanoutError(f"{label}: normalized event identities must be unique")
    for index, (g4_id, truth_id) in enumerate(zip(g4_ids, truth_ids), start=1):
        if g4_id != truth_id:
            raise FanoutError(
                f"{label}: event identity mismatch at row {index}: "
                f"G4={g4_id!r}, TRUTH_JET={truth_id!r}"
            )
    g4_hash = hashlib.sha256("\n".join(g4_ids[:5]).encode()).hexdigest()
    truth_hash = hashlib.sha256("\n".join(truth_ids[:5]).encode()).hexdigest()
    if g4_hash != truth_hash:
        raise FanoutError(f"{label}: first-five event identity hash differs")
    return g4_hash


def _validate_macro(path: Path, interaction: str, label: str) -> None:
    text = path.read_text(encoding="utf-8")
    active = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
    active = re.sub(r"//[^\n]*", "", active)
    indices = {
        int(value)
        for value in re.findall(
            r"INPUTREADHITS\s*::\s*listfile\s*\[\s*([0-9]+)\s*\]\s*=",
            active,
        )
    }
    if indices != {0, 4}:
        raise FanoutError(
            f"{label}: active INPUTREADHITS list indices must be exactly {{0, 4}}, "
            f"observed={sorted(indices)}"
        )
    embedding_two = bool(
        re.search(r"\badd_embedding_flag\s*\(\s*2\s*\)", active)
    )
    if interaction == "DI":
        if not embedding_two:
            raise FanoutError(f"{label}: DI macro must call add_embedding_flag(2)")
        if not re.search(r"\bnew\s+TruthJetInput\b", active):
            raise FanoutError(f"{label}: DI macro must reconstruct TruthJetInput")
    elif embedding_two:
        raise FanoutError(f"{label}: SI macro must not call add_embedding_flag(2)")


def source_contract(source_manifest: Path, output_root: str) -> dict[str, Any]:
    source_manifest = absolute_file(source_manifest, "source manifest")
    source = read_json(source_manifest, "source manifest")
    if source.get("schema") != SOURCE_SCHEMA or set(source) != {"schema", "common", "lanes"}:
        raise FanoutError(f"source manifest must use the exact {SOURCE_SCHEMA} schema")
    raw_common = source.get("common")
    if not isinstance(raw_common, dict) or set(raw_common) != set(COMMON_NAMES):
        raise FanoutError("source manifest common asset role set differs")
    common = {
        name: file_link(absolute_file(raw_common[name], f"common {name}"))
        for name in COMMON_NAMES
    }

    raw_lanes = source.get("lanes")
    if not isinstance(raw_lanes, list):
        raise FanoutError("source manifest lanes must be a list")
    # Validate raw cardinality and identifiers before constructing a mapping.
    # A dict comprehension would otherwise collapse a 13th duplicate row into
    # an apparently valid 12-lane manifest.
    if len(raw_lanes) != 12:
        raise FanoutError(
            "source manifest must contain exactly 12 raw canonical lane rows; "
            f"observed={len(raw_lanes)}"
        )
    raw_lane_ids: list[str] = []
    expected_fields = {
        "lane_id", "sample", "period", "interaction", *LANE_SOURCE_NAMES
    }
    for index, row in enumerate(raw_lanes):
        if not isinstance(row, dict) or set(row) != expected_fields:
            raise FanoutError(
                f"source manifest lane row {index} has unexpected fields"
            )
        lane_id = row.get("lane_id")
        if not isinstance(lane_id, str) or not lane_id:
            raise FanoutError(
                f"source manifest lane row {index} omits a valid lane_id"
            )
        raw_lane_ids.append(lane_id)
    duplicate_ids = sorted(
        lane_id for lane_id in set(raw_lane_ids) if raw_lane_ids.count(lane_id) > 1
    )
    if duplicate_ids:
        raise FanoutError(
            "source manifest duplicates raw lane_id values: "
            + ", ".join(duplicate_ids)
        )
    by_id: dict[str, dict[str, Any]] = {}
    seen_paths: dict[Path, tuple[tuple[str, str], str, str]] = {}
    for row in raw_lanes:
        lane_id = row["lane_id"]
        by_id[lane_id] = row

    expected = expected_lanes()
    if set(by_id) != {row[0] for row in expected} or len(by_id) != 12:
        raise FanoutError("source manifest must contain the exact 12 canonical lanes")
    lanes: list[dict[str, Any]] = []
    for lane_id, sample, period, interaction in expected:
        raw = by_id[lane_id]
        if (raw.get("sample"), raw.get("period"), raw.get("interaction")) != (
            sample, period, interaction
        ):
            raise FanoutError(f"{lane_id}: embedded physical identity differs")
        sources: dict[str, dict[str, str]] = {}
        list_rows: dict[str, list[str]] = {}
        for name in LANE_SOURCE_NAMES:
            path = absolute_file(raw[name], f"{lane_id} {name}")
            if any(char.isspace() for char in str(path)):
                raise FanoutError(f"{lane_id} {name} path contains whitespace")
            source_identity = (sample, interaction)
            prior = seen_paths.get(path)
            if prior is not None:
                prior_identity, prior_name, prior_lane = prior
                if prior_identity != source_identity or prior_name != name:
                    raise FanoutError(
                        f"{lane_id} {name} illegally shares {prior_lane} "
                        f"{prior_name} path {path}"
                    )
            else:
                seen_paths[path] = (source_identity, name, lane_id)
            sources[name] = file_link(path)
            if name.endswith("_list"):
                list_rows[name] = _list_rows(path, f"{lane_id} {name}")
        first_five_identity_sha256 = _paired_rows(
            list_rows["g4_full_list"], list_rows["truthjet_full_list"], lane_id
        )
        _validate_macro(Path(sources["ppg_macro"]["path"]), interaction, lane_id)
        sample_key = sample.lower()
        output_base = f"{output_root}/{sample_key}_{period}_{interaction.lower()}"
        lanes.append(
            {
                "lane_id": lane_id,
                "lane_key": f"{sample_key}_{period}_{interaction.lower()}",
                "sample": sample,
                "period": period,
                "interaction": interaction,
                "rows": 5,
                "first_five_event_identity_sha256": first_five_identity_sha256,
                "output_base": output_base,
                "sources": sources,
                "expected_evidence": {
                    "run_state": output_base + "/RUN_STATE",
                    "runtime_contract": output_base + "/paired_oracle_contract.json",
                    "candidate_csv": output_base + "/comparison/paired_oracle_candidates.csv",
                    "executable_aggregate": output_base + "/comparison/executable_aggregate.json",
                    "lane_receipt": output_base + "/condor_lane_receipt.json",
                },
            }
        )

    by_normalized_id = {lane["lane_id"]: lane for lane in lanes}
    for photon in (5, 10, 20):
        sample = f"Photon{photon}"
        for interaction in ("SI", "DI"):
            zero = by_normalized_id[
                f"photon:photon{photon}:0mrad:{interaction.lower()}"
            ]
            shifted = by_normalized_id[
                f"photon:photon{photon}:1p5mrad:{interaction.lower()}"
            ]
            for name in LANE_SOURCE_NAMES:
                if zero["sources"][name] != shifted["sources"][name]:
                    raise FanoutError(
                        f"{sample} {interaction} {name} must be exactly shared "
                        "between 0mrad and 1p5mrad"
                    )

        si = by_normalized_id[f"photon:photon{photon}:0mrad:si"]
        di = by_normalized_id[f"photon:photon{photon}:0mrad:di"]
        for name in ("g4_full_list", "truthjet_full_list"):
            si_rows = _list_rows(Path(si["sources"][name]["path"]), f"{sample} SI {name}")
            di_rows = _list_rows(Path(di["sources"][name]["path"]), f"{sample} DI {name}")
            if si_rows == di_rows:
                raise FanoutError(
                    f"{sample} SI and DI {name} source rows must be distinct"
                )
    return {
        "source_manifest": file_link(source_manifest),
        "common": common,
        "lanes": lanes,
    }


def provenance(chat_name: str, thread_id: str) -> dict[str, str]:
    values = {
        "RJ_CODEX_CHAT_NAME": chat_name,
        "RJ_CODEX_THREAD_ID": thread_id,
    }
    if (
        chat_name != chat_name.strip()
        or not SAFE_CHAT_NAME_RE.fullmatch(chat_name)
        or any(token in chat_name for token in ('$', '"', "'", "\\", ";"))
    ):
        raise FanoutError(
            "RJ_CODEX_CHAT_NAME must be an explicit Condor-safe title; "
            "spaces and | are allowed, submit syntax is not"
        )
    if thread_id != thread_id.strip() or not SAFE_THREAD_ID_RE.fullmatch(thread_id):
        raise FanoutError(
            "RJ_CODEX_THREAD_ID must be explicit and Condor-safe; "
            "use letters, digits, ._:+/@-"
        )
    return values


def condor_string(value: str, label: str) -> str:
    """Render one validated value as an HTCondor double-quoted string.

    Submit-file macro expansion happens before ClassAd parsing, so provenance
    validation rejects ``$`` and submit delimiters.  This renderer then owns
    the surrounding quotes, preserving spaces and ``|`` in human task names.
    """

    if any(char in value for char in ('\r', '\n', '$', '"', "'", "\\", ";")):
        raise FanoutError(f"{label} cannot be represented safely in a Condor string")
    return f'"{value}"'


def condor_environment(assignments: dict[str, str]) -> str:
    """Render validated key/value assignments in HTCondor's new syntax."""

    rendered: list[str] = []
    for name, value in assignments.items():
        if not re.fullmatch(r"[A-Z][A-Z0-9_]*", name):
            raise FanoutError(f"unsafe Condor environment variable name: {name!r}")
        if any(char in value for char in ('\r', '\n', '$', '"', "'", "\\", ";", "=")):
            raise FanoutError(
                f"{name} cannot be represented safely in the Condor environment"
            )
        rendered.append(f"{name}={value}")
    return f'"{";".join(rendered)}"'


def build_authorization(
    *,
    source_manifest: Path,
    paired_driver: Path,
    campaign_driver: Path,
    helper: Path,
    campaign_tag: str,
    output_root: str,
    evidence_dir: str,
    chat_name: str,
    thread_id: str,
) -> dict[str, Any]:
    if not SAFE_TAG_RE.fullmatch(campaign_tag):
        raise FanoutError(f"unsafe campaign tag: {campaign_tag!r}")
    output_root = absolute_directory_string(output_root, "output root")
    evidence_dir = absolute_directory_string(evidence_dir, "evidence directory")
    if output_root == evidence_dir:
        raise FanoutError("output root and evidence directory must differ")
    paired_driver = absolute_file(paired_driver, "paired driver", executable=True)
    campaign_driver = absolute_file(campaign_driver, "campaign driver", executable=True)
    helper = absolute_file(helper, "Condor fanout helper")
    source = source_contract(source_manifest, output_root)
    return {
        "campaign_tag": campaign_tag,
        "output_root": output_root,
        "evidence_dir": evidence_dir,
        **source,
        "paired_driver": file_link(paired_driver),
        "campaign_driver": file_link(campaign_driver),
        "fanout_helper": file_link(helper),
        "provenance": provenance(chat_name, thread_id),
        "bounded_contract": {
            "lane_count": 12,
            "jobs_per_lane": 1,
            "rows_per_lane": 5,
            "execution": "condor_source_locked_paired_executable",
            "automatic_merge": False,
            "promotion": False,
            "external_normalization": "forbidden",
            "automatic_retries": False,
            "job_control_interface": False,
        },
    }


def canonical_plan(authorization: dict[str, Any]) -> dict[str, Any]:
    authorization_sha = payload_sha256(authorization)
    return {
        "schema": PLAN_SCHEMA,
        "mode": "condor_fanout",
        "authorization": authorization,
        "authorization_sha256": authorization_sha,
        "submission_token": TOKEN_PREFIX + authorization_sha,
    }


def _rendered_json(payload: object) -> str:
    return json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"


def write_same_or_fail(path: Path, content: str, *, executable: bool = False) -> str:
    if not path.is_absolute():
        raise FanoutError(f"generated path must be absolute: {path}")
    if path.exists():
        if not path.is_file() or path.read_text(encoding="utf-8") != content:
            raise FanoutError(f"refusing to overwrite different generated file: {path}")
        if executable and not os.access(path, os.X_OK):
            path.chmod(0o500)
        return "UNCHANGED"
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        with temporary.open("x", encoding="utf-8") as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        temporary.chmod(0o500 if executable else 0o400)
        temporary.replace(path)
    finally:
        if temporary.exists():
            temporary.unlink()
    return "WROTE"


def generation_paths(evidence_dir: Path) -> dict[str, Path]:
    return {
        "plan": evidence_dir / "photon_condor_plan.json",
        "queue": evidence_dir / "photon_condor_queue.tsv",
        "wrapper": evidence_dir / "run_photon_condor_lane.sh",
        "submit": evidence_dir / "photon_condor.submit",
        "generation": evidence_dir / "photon_condor_generation.json",
        "submission_receipt": evidence_dir / "photon_condor_submission.json",
    }


def render_queue(plan: dict[str, Any]) -> str:
    # HTCondor's ``queue ... from`` consumes every non-comment row; a TSV
    # header would silently become a thirteenth job.
    rows = [
        f"{lane['lane_id']}\t{lane['lane_key']}"
        for lane in plan["authorization"]["lanes"]
    ]
    return "\n".join(rows) + "\n"


def render_wrapper(plan_path: Path, token: str, campaign_driver: Path) -> str:
    return "\n".join(
        (
            "#!/usr/bin/env bash",
            "set -Eeuo pipefail",
            "IFS=$'\\n\\t'",
            "umask 077",
            '[[ "$#" -eq 1 ]] || { echo "expected one lane_id" >&2; exit 64; }',
            "exec "
            + shlex.quote(str(campaign_driver))
            + " --condor-lane --plan "
            + shlex.quote(str(plan_path))
            + " --token "
            + shlex.quote(token)
            + ' --lane-id "$1"',
            "",
        )
    )


def render_submit(
    *,
    plan: dict[str, Any],
    wrapper: Path,
    queue: Path,
    log_dir: Path,
) -> str:
    auth = plan["authorization"]
    prov = auth["provenance"]
    return "\n".join(
        (
            "universe = vanilla",
            f"executable = {wrapper}",
            'arguments = "$(lane_id)"',
            f"initialdir = {Path(auth['campaign_driver']['path']).parents[4]}",
            "getenv = False",
            "environment = " + condor_environment(prov),
            "should_transfer_files = NO",
            "notification = Never",
            "request_cpus = 1",
            "request_memory = 4GB",
            "request_disk = 8GB",
            f"output = {log_dir}/$(lane_key).out",
            f"error = {log_dir}/$(lane_key).err",
            f"log = {log_dir}/fanout.log",
            "on_exit_hold = (ExitBySignal =?= True) || (ExitCode =!= 0)",
            "+JobBatchName = "
            + condor_string(auth["campaign_tag"], "JobBatchName"),
            '+RJ_PPG12_FANOUT = "exact_12_lane_paired_oracle"',
            "+RJ_CODEX_CHAT_NAME = "
            + condor_string(prov["RJ_CODEX_CHAT_NAME"], "RJ_CODEX_CHAT_NAME"),
            "+RJ_CODEX_THREAD_ID = "
            + condor_string(prov["RJ_CODEX_THREAD_ID"], "RJ_CODEX_THREAD_ID"),
            '+RJ_AUTOMATIC_RETRIES = False',
            f"queue lane_id,lane_key from {queue}",
            "",
        )
    )


def generate(args: argparse.Namespace) -> dict[str, Any]:
    helper = Path(__file__).resolve()
    authorization = build_authorization(
        source_manifest=args.source_manifest,
        paired_driver=args.paired_driver,
        campaign_driver=args.campaign_driver,
        helper=helper,
        campaign_tag=args.campaign_tag,
        output_root=args.output_root,
        evidence_dir=args.evidence_dir,
        chat_name=args.codex_chat_name,
        thread_id=args.codex_thread_id,
    )
    plan = canonical_plan(authorization)
    evidence_dir = Path(authorization["evidence_dir"])
    paths = generation_paths(evidence_dir)
    log_dir = evidence_dir / "condor_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    write_same_or_fail(paths["plan"], _rendered_json(plan))
    write_same_or_fail(paths["queue"], render_queue(plan))
    write_same_or_fail(
        paths["wrapper"],
        render_wrapper(
            paths["plan"],
            plan["submission_token"],
            Path(authorization["campaign_driver"]["path"]),
        ),
        executable=True,
    )
    write_same_or_fail(
        paths["submit"],
        render_submit(
            plan=plan,
            wrapper=paths["wrapper"],
            queue=paths["queue"],
            log_dir=log_dir,
        ),
    )
    generated_links = {
        name: file_link(paths[name]) for name in ("plan", "queue", "wrapper", "submit")
    }
    generation = {
        "schema": GENERATION_SCHEMA,
        "authorization_sha256": plan["authorization_sha256"],
        "submission_token_sha256": hashlib.sha256(
            plan["submission_token"].encode()
        ).hexdigest(),
        "authorized_job_count": 12,
        "artifacts": generated_links,
        "contract": {
            "one_job_per_lane": True,
            "automatic_merge": False,
            "promotion": False,
            "automatic_retries": False,
            "job_control_interface": False,
        },
    }
    generation["verification_payload_sha256"] = payload_sha256(generation)
    write_same_or_fail(paths["generation"], _rendered_json(generation))
    return {
        "plan": plan,
        "paths": {name: str(path) for name, path in paths.items()},
    }


def validate_plan(path: Path) -> dict[str, Any]:
    path = absolute_file(path, "Condor plan")
    plan = read_json(path, "Condor plan")
    if set(plan) != {
        "schema", "mode", "authorization", "authorization_sha256", "submission_token"
    } or plan.get("schema") != PLAN_SCHEMA or plan.get("mode") != "condor_fanout":
        raise FanoutError("Condor plan schema or field coverage differs")
    auth = plan.get("authorization")
    if not isinstance(auth, dict):
        raise FanoutError("Condor plan authorization is missing")
    rebuilt = build_authorization(
        source_manifest=Path(auth["source_manifest"]["path"]),
        paired_driver=Path(auth["paired_driver"]["path"]),
        campaign_driver=Path(auth["campaign_driver"]["path"]),
        helper=Path(auth["fanout_helper"]["path"]),
        campaign_tag=auth["campaign_tag"],
        output_root=auth["output_root"],
        evidence_dir=auth["evidence_dir"],
        chat_name=auth["provenance"]["RJ_CODEX_CHAT_NAME"],
        thread_id=auth["provenance"]["RJ_CODEX_THREAD_ID"],
    )
    if auth != rebuilt:
        raise FanoutError("Condor plan authorization differs from current source assets")
    auth_sha = payload_sha256(auth)
    if plan["authorization_sha256"] != auth_sha:
        raise FanoutError("Condor plan authorization hash is stale")
    expected_token = TOKEN_PREFIX + auth_sha
    if plan["submission_token"] != expected_token:
        raise FanoutError("Condor plan driver token is stale")
    return plan


def validate_generation(plan_path: Path, plan: dict[str, Any]) -> dict[str, Any]:
    paths = generation_paths(Path(plan["authorization"]["evidence_dir"]))
    if plan_path.resolve() != paths["plan"].resolve():
        raise FanoutError("Condor plan is outside its authorized evidence directory")
    generation = read_json(absolute_file(paths["generation"], "generation manifest"), "generation manifest")
    if generation.get("schema") != GENERATION_SCHEMA:
        raise FanoutError("generation manifest schema differs")
    unhashed = dict(generation)
    declared = unhashed.pop("verification_payload_sha256", None)
    if declared != payload_sha256(unhashed):
        raise FanoutError("generation manifest payload hash is stale")
    if generation.get("authorization_sha256") != plan["authorization_sha256"]:
        raise FanoutError("generation manifest binds a different authorization")
    artifacts = generation.get("artifacts")
    if not isinstance(artifacts, dict) or set(artifacts) != {"plan", "queue", "wrapper", "submit"}:
        raise FanoutError("generation manifest artifact coverage differs")
    for name, link in artifacts.items():
        observed = validate_link(link, f"generated {name}", executable=name == "wrapper")
        if observed["path"] != str(paths[name].resolve()):
            raise FanoutError(f"generated {name} path differs from the authorized path")
    if paths["queue"].read_text(encoding="utf-8") != render_queue(plan):
        raise FanoutError("generated Condor queue differs from the canonical 12 lanes")
    if paths["wrapper"].read_text(encoding="utf-8") != render_wrapper(
        paths["plan"],
        plan["submission_token"],
        Path(plan["authorization"]["campaign_driver"]["path"]),
    ):
        raise FanoutError("generated Condor wrapper differs")
    expected_submit = render_submit(
        plan=plan,
        wrapper=paths["wrapper"],
        queue=paths["queue"],
        log_dir=Path(plan["authorization"]["evidence_dir"]) / "condor_logs",
    )
    if paths["submit"].read_text(encoding="utf-8") != expected_submit:
        raise FanoutError("generated Condor submit description differs")
    return generation


def require_token(plan: dict[str, Any], token: str) -> None:
    if not TOKEN_RE.fullmatch(token) or token != plan["submission_token"]:
        raise FanoutError("provided token does not authorize this exact Condor plan")


def require_live_provenance(plan: dict[str, Any]) -> None:
    expected = plan["authorization"]["provenance"]
    observed = {
        "RJ_CODEX_CHAT_NAME": os.environ.get("RJ_CODEX_CHAT_NAME", ""),
        "RJ_CODEX_THREAD_ID": os.environ.get("RJ_CODEX_THREAD_ID", ""),
    }
    if observed != expected:
        raise FanoutError(
            "live RJ_CODEX provenance differs from the token-bound Condor plan"
        )


def require_live_submission_contract(plan: dict[str, Any]) -> None:
    auth = plan["authorization"]
    source = os.environ.get("RJ_PPG12_PAIRED_SOURCE_MANIFEST", "")
    if not source:
        raise FanoutError("live submission lacks RJ_PPG12_PAIRED_SOURCE_MANIFEST")
    observed = {
        "source_manifest": str(
            absolute_file(source, "live paired source manifest")
        ),
        "campaign_tag": os.environ.get("RJ_PPG12_PHOTON_CANARY_CAMPAIGN_TAG", ""),
        "output_root": os.environ.get("RJ_PPG12_PHOTON_CANARY_OUTPUT_ROOT", ""),
        "evidence_dir": os.environ.get("RJ_PPG12_PHOTON_CANARY_EVIDENCE_DIR", ""),
    }
    expected = {
        "source_manifest": auth["source_manifest"]["path"],
        "campaign_tag": auth["campaign_tag"],
        "output_root": auth["output_root"],
        "evidence_dir": auth["evidence_dir"],
    }
    if observed != expected:
        raise FanoutError("live submission contract differs from the token-bound plan")


def lane_args(plan: dict[str, Any], lane_id: str) -> tuple[dict[str, Any], list[str]]:
    matches = [lane for lane in plan["authorization"]["lanes"] if lane["lane_id"] == lane_id]
    if len(matches) != 1:
        raise FanoutError(f"lane_id is not one exact authorized lane: {lane_id!r}")
    lane = matches[0]
    common = plan["authorization"]["common"]
    args = [
        "--lane-id", lane["lane_id"],
        "--sample", lane["sample"],
        "--period", lane["period"],
        "--interaction", lane["interaction"],
        "--output-dir", lane["output_base"],
        "--ppg-macro", lane["sources"]["ppg_macro"]["path"],
        "--g4-full-list", lane["sources"]["g4_full_list"]["path"],
        "--truthjet-full-list", lane["sources"]["truthjet_full_list"]["path"],
        "--setup-script", common["setup_script"]["path"],
        "--apply-bdt", common["apply_bdt"]["path"],
        "--apply-config", common["apply_config"]["path"],
        "--base-e-model", common["base_e_model"]["path"],
        "--base-v3e-model", common["base_v3e_model"]["path"],
        "--npb-model", common["npb_model"]["path"],
        "--tower-mask", common["tower_mask"]["path"],
        "--recoil-runtime-manifest", common["recoil_runtime_manifest"]["path"],
        "--recoil-config", common["recoil_config"]["path"],
    ]
    return lane, args


def execute_lane(args: argparse.Namespace) -> int:
    plan_path = absolute_file(args.plan, "Condor plan")
    plan = validate_plan(plan_path)
    require_token(plan, args.token)
    require_live_provenance(plan)
    lane, driver_args = lane_args(plan, args.lane_id)
    output = Path(lane["output_base"])
    if output.exists():
        raise FanoutError(f"lane output already exists; automatic retry forbidden: {output}")
    paired_driver = Path(plan["authorization"]["paired_driver"]["path"])
    lane_plan = subprocess.run(
        [str(paired_driver), *driver_args],
        text=True,
        capture_output=True,
        check=False,
    )
    if lane_plan.returncode != 0:
        raise FanoutError(
            f"paired driver planning failed for {args.lane_id}: "
            f"{lane_plan.stderr.strip() or lane_plan.stdout.strip()}"
        )
    match = re.search(r"^  run_token: (ppg12-oracle:[0-9a-f]{64})$", lane_plan.stdout, re.M)
    if not match or not LANE_TOKEN_RE.fullmatch(match.group(1)):
        raise FanoutError(f"paired driver emitted no valid lane token for {args.lane_id}")
    completed = subprocess.run(
        [str(paired_driver), "--run", "--token", match.group(1), *driver_args],
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        raise FanoutError(f"paired executable failed for {args.lane_id}")
    evidence = lane["expected_evidence"]
    run_state = absolute_file(evidence["run_state"], f"{args.lane_id} RUN_STATE")
    if run_state.read_text(encoding="utf-8").strip() != "PASS":
        raise FanoutError(f"{args.lane_id} did not write RUN_STATE=PASS")
    links = {
        name: file_link(absolute_file(evidence[name], f"{args.lane_id} {name}"))
        for name in ("runtime_contract", "candidate_csv", "executable_aggregate", "run_state")
    }
    receipt = {
        "schema": RECEIPT_SCHEMA,
        "status": "PASS",
        "lane_id": args.lane_id,
        "authorization_sha256": plan["authorization_sha256"],
        "plan": file_link(plan_path),
        "paired_driver": plan["authorization"]["paired_driver"],
        "evidence": links,
    }
    receipt["verification_payload_sha256"] = payload_sha256(receipt)
    receipt_path = Path(evidence["lane_receipt"])
    write_same_or_fail(receipt_path, _rendered_json(receipt))
    print(f"PPG12_PAIRED_CONDOR_LANE_PASS lane_id={args.lane_id} receipt={receipt_path}")
    return 0


def validate_lane_receipt(
    plan_path: Path, plan: dict[str, Any], lane: dict[str, Any]
) -> dict[str, Any]:
    receipt_path = Path(lane["expected_evidence"]["lane_receipt"])
    if not receipt_path.is_file():
        if Path(lane["output_base"]).exists():
            return {
                "lane_id": lane["lane_id"],
                "status": "FAIL",
                "receipt": str(receipt_path),
                "reason": "lane output exists without a valid completion receipt",
            }
        return {"lane_id": lane["lane_id"], "status": "MISSING", "receipt": str(receipt_path)}
    try:
        receipt = read_json(receipt_path, f"{lane['lane_id']} receipt")
        if receipt.get("schema") != RECEIPT_SCHEMA or receipt.get("lane_id") != lane["lane_id"]:
            raise FanoutError("receipt schema or lane identity differs")
        unhashed = dict(receipt)
        declared = unhashed.pop("verification_payload_sha256", None)
        if declared != payload_sha256(unhashed):
            raise FanoutError("receipt verification hash is stale")
        if receipt.get("authorization_sha256") != plan["authorization_sha256"]:
            raise FanoutError("receipt binds a different authorization")
        if validate_link(receipt.get("plan"), "receipt plan") != file_link(plan_path):
            raise FanoutError("receipt binds a different plan")
        if receipt.get("paired_driver") != plan["authorization"]["paired_driver"]:
            raise FanoutError("receipt binds a different paired driver")
        evidence = receipt.get("evidence")
        if not isinstance(evidence, dict) or set(evidence) != {
            "runtime_contract", "candidate_csv", "executable_aggregate", "run_state"
        }:
            raise FanoutError("receipt evidence coverage differs")
        for name, link in evidence.items():
            observed = validate_link(link, f"receipt {name}")
            if observed["path"] != str(Path(lane["expected_evidence"][name]).resolve()):
                raise FanoutError(f"receipt {name} path differs")
        if Path(evidence["run_state"]["path"]).read_text(encoding="utf-8").strip() != "PASS":
            raise FanoutError("RUN_STATE is not PASS")
        return {
            "lane_id": lane["lane_id"],
            "status": "PASS",
            "receipt": str(receipt_path),
            "receipt_sha256": sha256(receipt_path),
        }
    except (FanoutError, OSError) as exc:
        return {
            "lane_id": lane["lane_id"],
            "status": "FAIL",
            "receipt": str(receipt_path),
            "reason": str(exc),
        }


def audit(args: argparse.Namespace) -> int:
    plan_path = absolute_file(args.plan, "Condor plan")
    plan = validate_plan(plan_path)
    validate_generation(plan_path, plan)
    rows = [
        validate_lane_receipt(plan_path, plan, lane)
        for lane in plan["authorization"]["lanes"]
    ]
    counts = {status: sum(row["status"] == status for row in rows) for status in ("PASS", "MISSING", "FAIL")}
    status = "PASS" if counts["PASS"] == 12 else "FAIL" if counts["FAIL"] else "INCOMPLETE"
    report = {
        "schema": AUDIT_SCHEMA,
        "status": status,
        "authorization_sha256": plan["authorization_sha256"],
        "plan": file_link(plan_path),
        "counts": counts,
        "lanes": rows,
        "read_only": True,
    }
    report["verification_payload_sha256"] = payload_sha256(report)
    print(_rendered_json(report), end="")
    return 0 if status == "PASS" else 3 if status == "INCOMPLETE" else 2


def submit(args: argparse.Namespace) -> int:
    plan_path = absolute_file(args.plan, "Condor plan")
    plan = validate_plan(plan_path)
    require_token(plan, args.token)
    require_live_provenance(plan)
    require_live_submission_contract(plan)
    validate_generation(plan_path, plan)
    output_root = Path(plan["authorization"]["output_root"])
    if output_root.exists():
        raise FanoutError(f"output root already exists; refusing duplicate fanout: {output_root}")
    paths = generation_paths(Path(plan["authorization"]["evidence_dir"]))
    if paths["submission_receipt"].exists():
        raise FanoutError(
            "submission receipt already exists; refusing a second condor_submit invocation: "
            f"{paths['submission_receipt']}"
        )
    command = absolute_file(args.condor_submit, "condor_submit", executable=True)
    completed = subprocess.run(
        [str(command), "-terse", str(paths["submit"])],
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        raise FanoutError(
            "condor_submit failed: " + (completed.stderr.strip() or completed.stdout.strip())
        )
    receipt = {
        "schema": SUBMISSION_SCHEMA,
        "status": "SUBMITTED",
        "authorized_job_count": 12,
        "authorization_sha256": plan["authorization_sha256"],
        "plan": file_link(plan_path),
        "submit_file": file_link(paths["submit"]),
        "condor_submit": file_link(command),
        "stdout": completed.stdout.strip(),
        "stderr": completed.stderr.strip(),
        "automatic_retries": False,
        "automatic_merge": False,
    }
    write_same_or_fail(paths["submission_receipt"], _rendered_json(receipt))
    print(
        "PPG12_PAIRED_CONDOR_SUBMITTED jobs=12 "
        f"receipt={paths['submission_receipt']} readback={completed.stdout.strip()}"
    )
    return 0


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(description=__doc__)
    commands = root.add_subparsers(dest="command", required=True)
    plan = commands.add_parser("plan")
    plan.add_argument("--source-manifest", required=True, type=Path)
    plan.add_argument("--paired-driver", required=True, type=Path)
    plan.add_argument("--campaign-driver", required=True, type=Path)
    plan.add_argument("--campaign-tag", required=True)
    plan.add_argument("--output-root", required=True)
    plan.add_argument("--evidence-dir", required=True)
    plan.add_argument("--codex-chat-name", required=True)
    plan.add_argument("--codex-thread-id", required=True)
    submit_parser = commands.add_parser("submit")
    submit_parser.add_argument("--plan", required=True, type=Path)
    submit_parser.add_argument("--token", required=True)
    submit_parser.add_argument("--condor-submit", required=True, type=Path)
    lane = commands.add_parser("lane")
    lane.add_argument("--plan", required=True, type=Path)
    lane.add_argument("--token", required=True)
    lane.add_argument("--lane-id", required=True)
    audit_parser = commands.add_parser("audit")
    audit_parser.add_argument("--plan", required=True, type=Path)
    return root


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        if args.command == "plan":
            result = generate(args)
            print(
                "PPG12_PAIRED_CONDOR_PLAN jobs=12 "
                f"plan={result['paths']['plan']} submit={result['paths']['submit']} "
                f"run_token={result['plan']['submission_token']}"
            )
            return 0
        if args.command == "submit":
            return submit(args)
        if args.command == "lane":
            return execute_lane(args)
        if args.command == "audit":
            return audit(args)
        raise FanoutError(f"unsupported command {args.command}")
    except (FanoutError, OSError) as exc:
        print(f"PPG12_PAIRED_CONDOR_ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
