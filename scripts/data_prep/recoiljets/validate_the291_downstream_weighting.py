#!/usr/bin/env python3
"""Fail closed on drift in the THE-291 Au+Au inclusive downstream surface."""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
from typing import Any, Mapping, Sequence


SCRIPTS_ROOT = Path(__file__).resolve().parents[2]
if str(SCRIPTS_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_ROOT))


SCHEMA = "THE291AuAuEmbeddedInclusiveDownstreamHardeningValidationV1"
POLICY_SCHEMA = "THE291AuAuEmbeddedInclusiveDownstreamSurfacePolicyV1"
POLICY_STATUS = "PASS_FAIL_CLOSED_SURFACE_SEALED"
POLICY_SCOPE = (
    "all repository source files that mention both Au+Au embedded-inclusive inputs "
    "and downstream weighting or output behavior, collaborator weight-contract "
    "notes/examples, and the tracked CI/allow-list gate"
)
POLICY_CHANGE_RULE = (
    "any added, removed, renamed, or byte-changed surface file requires explicit "
    "policy review and resealing"
)
POLICY_PROHIBITIONS = (
    "no Au+Au embedded-inclusive nominal consumer may reconstruct source or centrality factors",
    "no Au+Au embedded-inclusive nominal plot may accept a centrality-only receipt",
    "no nominal PNG or table may omit the transitive downstream artifact receipt",
    "no collaborator example may treat Tree event_weight as a complete Au+Au inclusive weight",
)
DEFAULT_POLICY = Path(__file__).with_name("the291_downstream_surface_v1.json")
DEFAULT_REPO = Path(__file__).resolve().parents[3]
DEFAULT_SOURCE_ASSEMBLY = DEFAULT_REPO / (
    "dataOutput/the243_golden_ppg_analysis_closure_deck_20260821/"
    "stitching_updates_current_20260820/"
    "THE248_SCHEMA10_EMBEDDED_INCLUSIVE_STITCHING_ASSEMBLED_V1.json"
)
DEFAULT_SOURCE_RECEIPT = DEFAULT_SOURCE_ASSEMBLY.with_name(
    "THE248_SCHEMA10_EMBEDDED_INCLUSIVE_STITCHING_ASSEMBLY_RECEIPT_V1.json"
)
EXPECTED_SOURCE_ASSEMBLY_SHA256 = (
    "a79bae1b17a1c921e9a80b470c887340e6c76e497f11ad3c84222a57761b5206"
)
EXPECTED_SOURCE_RECEIPT_SHA256 = (
    "9b694d37b6cf85ccc39d4d0c2696168490328a08e6a378e7e2dd7b494ad3b9c5"
)

SOURCE_SUFFIXES = {".py", ".C", ".cc", ".cpp", ".h", ".sh"}
COLLABORATOR_CONTRACT_PATHS = {
    "scripts/data_prep/recoiljets/collaboration/assets/README.md",
    "scripts/data_prep/recoiljets/collaboration/assets/BRANCH_REFERENCE.md",
    "scripts/data_prep/recoiljets/collaboration/assets/examples/quick_start.C",
    "scripts/data_prep/recoiljets/collaboration/assets/examples/purity_regions.C",
}
GOVERNANCE_CONTRACT_PATHS = {
    ".gitignore",
    ".github/workflows/auau-inclusive-weight-contract.yml",
}
ROLE_BY_PATH = {
    "scripts/data_prep/recoiljets/auau_embedded_inclusive_schema10_weighting.py":
        "canonical_weight_provider",
    "scripts/data_prep/recoiljets/build_auau_embedded_inclusive_reco_cluster_source_fraction_payload.py":
        "canonical_payload_builder",
    "scripts/data_prep/recoiljets/assemble_the248_schema10_stitching_products.py":
        "canonical_source_stitch_assembler",
    "scripts/data_prep/recoiljets/assemble_the243_schema10_response_products.py":
        "canonical_response_router",
    "scripts/data_prep/recoiljets/pp_schema10_weighting.py":
        "pp_only_rejects_auau_inclusive",
    "scripts/data_prep/recoiljets/reduce_the248_schema10_stitching.py":
        "raw_source_reducer",
    "scripts/data_prep/recoiljets/the243_h70_corrections.py":
        "photon_signal_only_router",
    "scripts/data_prep/recoiljets/build_the255_auau_offline_package.py":
        "photon_signal_only_router",
    "scripts/plotting/canonical_schema10_baseline.py": "canonical_stitch_guard",
    "scripts/plotting/plot_label_contract.py": "canonical_plot_contract",
    "scripts/plotting/truth_purity/make_auau_current_reference_inclusive_candidate_yields_by_centrality.py":
        "canonical_nominal_plot_consumer",
    "scripts/slides/the243/stitching/make_the243_current_stitch_contract_slides.py":
        "canonical_stitch_diagnostic",
    "scripts/slides/the243/stitching/make_the248_current_schema10_stitching_closure_slides.py":
        "canonical_stitch_diagnostic",
    "scripts/slides/wp_gammajets/stitching/make_embedded_stitching_closure_slide.py":
        "legacy_diagnostic_only",
    "scripts/data_prep/recoiljets/validate_the291_downstream_weighting.py":
        "surface_enforcement",
    "scripts/data_prep/recoiljets/tests/test_validate_the291_downstream_weighting.py":
        "contract_test",
}

REQUIRED_TOKENS = {
    ".gitignore": (
        "!/.github/workflows/auau-inclusive-weight-contract.yml",
    ),
    ".github/workflows/auau-inclusive-weight-contract.yml": (
        "numpy==1.26.4",
        "validate_the291_downstream_weighting.py",
        "test_auau_embedded_inclusive_schema10_weighting",
        "test_validate_the291_downstream_weighting",
        "test_plot_label_contract",
    ),
    "scripts/data_prep/recoiljets/auau_embedded_inclusive_schema10_weighting.py": (
        "CanonicalSchema10AuAuEmbeddedInclusiveAnalysisWeightReceiptV1",
        "CanonicalSchema10AuAuEmbeddedInclusiveDownstreamArtifactReceiptV1",
        "validate_analysis_weight_payload",
        "write_downstream_artifact_receipt",
        "load_downstream_artifact_receipt",
        "sigma_eff/Npass",
        "sigma_eff/Ngen",
    ),
    "scripts/plotting/plot_label_contract.py": (
        "simulation_family",
        "inclusive_background",
        "centrality-only receipt is forbidden",
        "load_analysis_weight_receipt",
    ),
    "scripts/plotting/truth_purity/make_auau_current_reference_inclusive_candidate_yields_by_centrality.py": (
        'simulation_family="inclusive_background"',
        "def read_candidate_block",
        "allow_pickle=False",
        "write_downstream_artifact_receipt",
        "load_downstream_artifact_receipt",
        'path.open("x"',
        "os.link",
    ),
    "scripts/data_prep/recoiljets/assemble_the248_schema10_stitching_products.py": (
        "load_source_stitch_artifact",
        "CANONICAL_SOURCE_STITCH_ONLY__CENTRALITY_NOT_APPLIED",
        "nominal_downstream_analysis_ready",
    ),
    "scripts/data_prep/recoiljets/assemble_the243_schema10_response_products.py": (
        "load_source_stitch_artifact",
        "sigma_eff/Npass",
    ),
    "scripts/data_prep/recoiljets/pp_schema10_weighting.py": (
        "pp_only_rejects_auau_inclusive",
        "do not normalize Au+Au embedded inclusive jets",
    ),
    "scripts/slides/the243/stitching/make_the248_current_schema10_stitching_closure_slides.py": (
        "THE248Schema10StitchingAssemblyTerminalReceiptV4",
        "load_source_stitch_artifact",
        "legacy cross_section/generated_events weight is forbidden",
    ),
    "scripts/slides/wp_gammajets/stitching/make_embedded_stitching_closure_slide.py": (
        'THE291_CONTRACT_ROLE = "legacy_diagnostic_only"',
        "--historical-diagnostic-only",
        "HISTORICAL NON-NOMINAL",
    ),
}

FORBIDDEN_TOKENS_BY_PATH = {
    "scripts/plotting/truth_purity/make_auau_current_reference_inclusive_candidate_yields_by_centrality.py": (
        "agent_context",
        "unified_response",
    ),
}


class HardeningValidationError(ValueError):
    """Raised when a sealed surface or its semantics differs."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as error:
        raise HardeningValidationError(f"cannot read protected surface file: {path}") from error


def _is_surface_source(path: Path, text: str) -> bool:
    if path.suffix not in SOURCE_SUFFIXES:
        return False
    lowered = text.lower()
    has_auau = "auau" in lowered or "au+au" in lowered
    has_inclusive = any(token in lowered for token in (
        "embedded_inclusive",
        "embedded inclusive",
        "embedded-inclusive",
        "auau_jet12",
        "auau_jet20",
        "auau_jet30",
        "auau_jet40",
    ))
    has_downstream_behavior = any(token in lowered for token in (
        "weight", "stitch", "sumw", "hist", "plot", "render", "response", "correction", "unfold"
    ))
    return has_auau and has_inclusive and has_downstream_behavior


def discover_surface(repo: Path) -> dict[str, str]:
    surface: dict[str, str] = {}
    scripts = repo / "scripts"
    for path in sorted(scripts.rglob("*")):
        if (
            not path.is_file()
            or path.is_symlink()
            or path.suffix not in SOURCE_SUFFIXES
            or "__pycache__" in path.parts
        ):
            continue
        relative = path.relative_to(repo).as_posix()
        text = _read(path)
        if _is_surface_source(path, text):
            surface[relative] = text
    for relative in sorted(COLLABORATOR_CONTRACT_PATHS | GOVERNANCE_CONTRACT_PATHS):
        path = repo / relative
        if not path.is_file() or path.is_symlink():
            raise HardeningValidationError(
                f"required explicit contract file is missing: {relative}"
            )
        surface[relative] = _read(path)
    return dict(sorted(surface.items()))


def role_for_path(relative: str) -> str:
    if relative in ROLE_BY_PATH:
        return ROLE_BY_PATH[relative]
    if relative in COLLABORATOR_CONTRACT_PATHS:
        return "collaborator_contract"
    if relative in GOVERNANCE_CONTRACT_PATHS:
        return "governance_contract"
    if "/tests/" in relative or Path(relative).name.startswith("test_"):
        return "contract_or_legacy_test_frozen"
    if relative.endswith(".sh"):
        return "legacy_infrastructure_frozen"
    return "legacy_analysis_surface_frozen"


def build_policy(repo: Path) -> dict[str, Any]:
    surface = discover_surface(repo)
    return {
        "schema": POLICY_SCHEMA,
        "status": POLICY_STATUS,
        "scope": POLICY_SCOPE,
        "change_rule": POLICY_CHANGE_RULE,
        "entries": [
            {
                "path": relative,
                "role": role_for_path(relative),
                "sha256": sha256_file(repo / relative),
            }
            for relative in surface
        ],
        "prohibitions": list(POLICY_PROHIBITIONS),
    }


def _load_json(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise HardeningValidationError(f"invalid JSON object: {path}") from error
    if not isinstance(payload, dict):
        raise HardeningValidationError(f"JSON object required: {path}")
    return payload


def _direct_bad_divisions(path: Path, text: str) -> list[str]:
    if path.suffix != ".py":
        pattern = re.compile(
            r"\b(?:cross_section\w*|xsec\w*|sigma_eff\w*)\s*/\s*"
            r"(?:generated_events\w*|ngen\w*)\b",
            re.IGNORECASE,
        )
        return [
            f"forbidden cross-section/Ngen division at {path}:{text.count(chr(10), 0, match.start()) + 1}"
            for match in pattern.finditer(text)
        ]
    try:
        tree = ast.parse(text, filename=str(path))
    except SyntaxError as error:
        return [f"protected Python source is not parseable: {path}: {error}"]
    findings: list[str] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.BinOp) or not isinstance(node.op, ast.Div):
            continue
        left = ast.unparse(node.left).lower()
        right = ast.unparse(node.right).lower()
        if (
            any(token in left for token in ("cross_section", "xsec", "sigma_eff"))
            and any(token in right for token in ("generated_events", "ngen"))
        ):
            findings.append(
                f"forbidden cross-section/Ngen division at {path}:{getattr(node, 'lineno', 0)}"
            )
    return findings


def validate_repository(
    repo: Path,
    policy_path: Path,
    *,
    require_canonical_surface: bool = True,
) -> tuple[dict[str, Any], list[str]]:
    errors: list[str] = []
    try:
        policy = _load_json(policy_path)
        surface = discover_surface(repo)
    except HardeningValidationError as error:
        return {}, [str(error)]
    if (
        set(policy) != {"schema", "status", "scope", "change_rule", "entries", "prohibitions"}
        or policy.get("schema") != POLICY_SCHEMA
        or policy.get("status") != POLICY_STATUS
        or not isinstance(policy.get("entries"), list)
    ):
        return {}, ["THE-291 surface policy schema/status/fields differ"]
    if policy.get("scope") != POLICY_SCOPE:
        errors.append("THE-291 surface policy scope differs")
    if policy.get("change_rule") != POLICY_CHANGE_RULE:
        errors.append("THE-291 surface policy change rule differs")
    if policy.get("prohibitions") != list(POLICY_PROHIBITIONS):
        errors.append("THE-291 surface policy prohibitions differ")
    entries: dict[str, Mapping[str, Any]] = {}
    for entry in policy["entries"]:
        if not isinstance(entry, Mapping) or set(entry) != {"path", "role", "sha256"}:
            errors.append("THE-291 surface policy contains an invalid entry")
            continue
        relative = entry.get("path")
        if not isinstance(relative, str) or relative in entries:
            errors.append("THE-291 surface policy contains a duplicate or invalid path")
            continue
        entries[relative] = entry
    if list(entries) != sorted(entries):
        errors.append("THE-291 surface policy entries are not canonically sorted")
    discovered = set(surface)
    sealed = set(entries)
    if discovered != sealed:
        errors.append(
            "THE-291 protected surface inventory drifted: added={} removed={}".format(
                sorted(discovered - sealed), sorted(sealed - discovered)
            )
        )
    for relative in sorted(discovered & sealed):
        entry = entries[relative]
        expected_role = role_for_path(relative)
        if entry.get("role") != expected_role:
            errors.append(f"THE-291 surface role differs for {relative}")
        actual_hash = sha256_file(repo / relative)
        if entry.get("sha256") != actual_hash:
            errors.append(f"THE-291 protected surface hash drifted for {relative}")

    if require_canonical_surface:
        missing_required = sorted(set(ROLE_BY_PATH) - discovered)
        if missing_required:
            errors.append(f"THE-291 required canonical surface is missing: {missing_required}")
    for relative, tokens in REQUIRED_TOKENS.items():
        if relative not in surface:
            continue
        for token in tokens:
            if token not in surface[relative]:
                errors.append(f"THE-291 required token missing from {relative}: {token}")
    for relative, tokens in FORBIDDEN_TOKENS_BY_PATH.items():
        if relative not in surface:
            continue
        for token in tokens:
            if token in surface[relative]:
                errors.append(f"THE-291 forbidden private dependency in {relative}: {token}")

    division_exempt_roles = {
        "canonical_weight_provider",
        "canonical_response_router",
        "pp_only_rejects_auau_inclusive",
        "photon_signal_only_router",
        "contract_test",
        "contract_or_legacy_test_frozen",
    }
    for relative, text in surface.items():
        if role_for_path(relative) not in division_exempt_roles:
            errors.extend(_direct_bad_divisions(repo / relative, text))

    readme = surface.get(
        "scripts/data_prep/recoiljets/collaboration/assets/README.md", ""
    )
    branch_reference = surface.get(
        "scripts/data_prep/recoiljets/collaboration/assets/BRANCH_REFERENCE.md", ""
    )
    quick_start = surface.get(
        "scripts/data_prep/recoiljets/collaboration/assets/examples/quick_start.C", ""
    )
    purity = surface.get(
        "scripts/data_prep/recoiljets/collaboration/assets/examples/purity_regions.C", ""
    )
    for label, text in (("README", readme), ("BRANCH_REFERENCE", branch_reference)):
        for token in (
            "NOT_READY_FROM_PACKAGE_ALONE",
            "ownership_effective_cross_section / Npass",
            "centrality",
            "event_weight",
        ):
            if token not in text:
                errors.append(f"collaborator {label} lacks required fail-closed weight note: {token}")
    forbidden_example_patterns = {
        "quick_start.C": (
            r"SetBranchAddress\(\s*\"event_weight\"",
            r"Fill\(\s*xjgamma\s*,",
            r"weighted pairs",
        ),
        "purity_regions.C": (
            r"SetBranchAddress\(\s*\"event_weight\"",
            r"\.event_weight",
            r"(?<!un)weighted event-leading",
        ),
    }
    for label, text in (("quick_start.C", quick_start), ("purity_regions.C", purity)):
        for pattern in forbidden_example_patterns[label]:
            if re.search(pattern, text, re.IGNORECASE):
                errors.append(f"collaborator {label} still applies or advertises Tree event_weight")
        if "raw/unweighted" not in text.lower():
            errors.append(f"collaborator {label} is not visibly raw/unweighted")

    return {
        "policy_path": str(policy_path.resolve()),
        "policy_sha256": sha256_file(policy_path),
        "surface_file_count": len(surface),
        "role_counts": {
            role: sum(1 for path in surface if role_for_path(path) == role)
            for role in sorted({role_for_path(path) for path in surface})
        },
        "collaborator_contract_hashes": {
            relative: sha256_file(repo / relative)
            for relative in sorted(COLLABORATOR_CONTRACT_PATHS)
        },
    }, errors


def validate_live_source(
    assembly_path: Path,
    receipt_path: Path,
) -> tuple[dict[str, Any], list[str]]:
    errors: list[str] = []
    if not assembly_path.is_file() or not receipt_path.is_file():
        return {}, ["accepted embedded-inclusive source-stitch artifact or receipt is missing"]
    assembly_hash = sha256_file(assembly_path)
    receipt_hash = sha256_file(receipt_path)
    if assembly_hash != EXPECTED_SOURCE_ASSEMBLY_SHA256:
        errors.append("accepted embedded-inclusive source assembly SHA-256 differs")
    if receipt_hash != EXPECTED_SOURCE_RECEIPT_SHA256:
        errors.append("accepted embedded-inclusive source receipt SHA-256 differs")
    samples: dict[str, Any] = {}
    if not errors:
        try:
            from data_prep.recoiljets.auau_embedded_inclusive_schema10_weighting import (
                load_source_stitch_artifact,
            )

            source = load_source_stitch_artifact(
                assembly_path, receipt_path, verify_authority_file=True
            )
            samples = {
                sample.sample_id: {
                    "normalization_denominator_events": sample.normalization_denominator_events,
                    "stitching_weight_pb_per_owned_event": sample.stitching_weight_pb_per_owned_event,
                    "ownership_window_gev": [
                        sample.ownership_low_gev,
                        sample.ownership_high_gev,
                    ],
                }
                for sample in source.samples
            }
        except (ImportError, OSError, ValueError) as error:
            errors.append(f"accepted embedded-inclusive source chain failed deep validation: {error}")
    return {
        "assembly_path": str(assembly_path.resolve()),
        "assembly_sha256": assembly_hash,
        "receipt_path": str(receipt_path.resolve()),
        "receipt_sha256": receipt_hash,
        "samples": samples,
    }, errors


def git_state(repo: Path) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, command in (
        ("head", ["git", "rev-parse", "HEAD"]),
        ("status_porcelain", ["git", "status", "--short"]),
    ):
        completed = subprocess.run(
            command, cwd=repo, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        result[key] = completed.stdout.strip() if completed.returncode == 0 else None
    return result


def write_receipt(path: Path, payload: Mapping[str, Any]) -> None:
    path = path.resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    parser.add_argument("--emit-policy", action="store_true")
    parser.add_argument("--receipt", type=Path)
    parser.add_argument("--source-assembly", type=Path, default=DEFAULT_SOURCE_ASSEMBLY)
    parser.add_argument("--source-receipt", type=Path, default=DEFAULT_SOURCE_RECEIPT)
    parser.add_argument("--skip-live-source", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    repo = args.repo.resolve()
    if args.emit_policy:
        print(json.dumps(build_policy(repo), indent=2, sort_keys=False))
        return 0
    repository, errors = validate_repository(repo, args.policy.resolve())
    live_source: dict[str, Any] = {"status": "NOT_RUN"}
    if not args.skip_live_source:
        live_source, findings = validate_live_source(
            args.source_assembly.resolve(), args.source_receipt.resolve()
        )
        errors.extend(findings)
    payload = {
        "schema": SCHEMA,
        "status": "pass" if not errors else "fail",
        "scope": "THE-291 downstream code, plot enforcement, collaborator notes/examples, and accepted source-stitch live readback",
        "repository": str(repo),
        "git": git_state(repo),
        "surface_validation": repository,
        "private_control_plane_validation": {
            "status": "NOT_APPLICABLE_PUBLIC_REPOSITORY_BOUNDARY"
        },
        "accepted_source_stitch_validation": live_source,
        "centrality_weight_payload_materialization": "NOT_RUN_NO_READY_CANONICAL_CENTRALITY_RECEIPT",
        "nominal_downstream_plot_materialization": "NOT_RUN_NO_COMPLETE_COMPOSITE_WEIGHT_PAYLOAD",
        "published_collaborator_package_update": "NOT_RUN_REQUIRES_EXPLICIT_REMOTE_MUTATION_AUTHORITY",
        "errors": errors,
    }
    if args.receipt is not None:
        write_receipt(args.receipt, payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
