#!/usr/bin/env python3
"""Preflight THE-79 PHENIX-like R_AA campaign inputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path
from typing import Any


EXPECTED_MODEL_ROOT = Path(
    "/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/"
    "the79_fullsrc_binned14_ppg12pt5to40_20260627_2052"
)
EXPECTED_WP = "auauEtFineCent7BDT|centlinear|0.5113985411|0.0011136541|8|40|1"
EXPECTED_AUAU_ROW = [["newPPG12", "auauEtFineCent7BDT", "auauBDTSideband"]]
EXPECTED_SIDE_MODE = "relativeToTight"
EXPECTED_SIDE_MIN_OFFSET = -0.20
EXPECTED_SIDE_MAX_OFFSET = -0.03
EXPECTED_PT_EDGES = [5, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36, 40]
FINAL_PT_BINS = [8, 8.5, 9, 9.5, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36, 40]
EXPECTED_CENT7 = [0, 10, 20, 30, 40, 50, 60, 80]
EXPECTED_AUAU_FEATURES = [
    "cluster_Et",
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "cluster_weta33_cogx",
    "cluster_wphi33_cogx",
    "vertexz",
    "cluster_Eta",
    "e11_over_e33",
    "cluster_et1",
    "cluster_et2",
    "cluster_et3",
    "cluster_et4",
    "e32_over_e35",
    "centrality",
]
EXPECTED_PP_FEATURES = [
    "cluster_Et",
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
]
REQUIRED_TAA_COLUMNS = [
    "cent_low",
    "cent_high",
    "N_coll",
    "TAA_mb_inv",
    "TAA_unc",
    "sigmaNN_mb",
    "source",
    "version",
]


def strip_comment(line: str) -> str:
    in_quote = False
    quote_char = ""
    out = []
    for ch in line:
        if ch in {"'", '"'}:
            if in_quote and ch == quote_char:
                in_quote = False
            elif not in_quote:
                in_quote = True
                quote_char = ch
        if ch == "#" and not in_quote:
            break
        out.append(ch)
    return "".join(out).rstrip()


def parse_scalar(value: str) -> Any:
    value = value.strip()
    if value == "":
        return ""
    if value in {"true", "True"}:
        return True
    if value in {"false", "False"}:
        return False
    if (value.startswith('"') and value.endswith('"')) or (
        value.startswith("'") and value.endswith("'")
    ):
        return value[1:-1]
    try:
        if any(ch in value for ch in [".", "e", "E"]):
            return float(value)
        return int(value)
    except ValueError:
        return value


def parse_inline_list(value: str) -> list[Any]:
    value = value.strip()
    if not (value.startswith("[") and value.endswith("]")):
        return []
    body = value[1:-1].strip()
    if not body:
        return []
    parts = []
    token = []
    in_quote = False
    quote_char = ""
    for ch in body:
        if ch in {"'", '"'}:
            if in_quote and ch == quote_char:
                in_quote = False
            elif not in_quote:
                in_quote = True
                quote_char = ch
            token.append(ch)
            continue
        if ch == "," and not in_quote:
            parts.append(parse_scalar("".join(token).strip()))
            token = []
            continue
        token.append(ch)
    if token:
        parts.append(parse_scalar("".join(token).strip()))
    return parts


def parse_simple_yaml(path: Path) -> dict[str, Any]:
    lines = path.read_text().splitlines()
    data: dict[str, Any] = {}
    i = 0
    while i < len(lines):
        line = strip_comment(lines[i])
        if not line.strip() or line.startswith(" ") or ":" not in line:
            i += 1
            continue
        key, rhs = line.split(":", 1)
        key = key.strip()
        rhs = rhs.strip()
        if rhs.startswith("["):
            data[key] = parse_inline_list(rhs)
        elif rhs:
            data[key] = parse_scalar(rhs)
        else:
            rows = []
            j = i + 1
            while j < len(lines):
                child = strip_comment(lines[j])
                if not child.strip():
                    j += 1
                    continue
                if not child.startswith(" "):
                    break
                stripped = child.strip()
                if stripped.startswith("-"):
                    item = stripped[1:].strip()
                    if item.startswith("["):
                        rows.append(parse_inline_list(item))
                    else:
                        rows.append(parse_scalar(item))
                j += 1
            data[key] = rows
            i = j - 1
        i += 1
    return data


def as_numbers(values: Any) -> list[float]:
    if not isinstance(values, list):
        return []
    out = []
    for value in values:
        try:
            out.append(float(value))
        except (TypeError, ValueError):
            return []
    return out


def assert_equal(errors: list[str], label: str, got: Any, expected: Any) -> None:
    if got != expected:
        errors.append(f"{label}: got {got!r}, expected {expected!r}")


def validate_yaml(pp_yaml: Path, auau_yaml: Path) -> tuple[list[str], list[str], dict[str, Any]]:
    errors: list[str] = []
    warnings: list[str] = []
    pp = parse_simple_yaml(pp_yaml)
    auau = parse_simple_yaml(auau_yaml)

    assert_equal(errors, "pp photon_id_sets", pp.get("photon_id_sets"), [["newPPG12", "newPPG12", "newPPG12"]])
    assert_equal(errors, "pp tight_bdt_features", pp.get("tight_bdt_features"), EXPECTED_PP_FEATURES)
    if str(pp.get("tight_bdt_model_file", "")).find("model_base_v3E_split_single_tmva.root") < 0:
        errors.append(f"pp tight_bdt_model_file is not the PPG12 baseV3E split model: {pp.get('tight_bdt_model_file')}")
    assert_equal(errors, "pp final pT bins", as_numbers(pp.get("unfold_reco_photon_pt_bins")), FINAL_PT_BINS)

    assert_equal(errors, "AuAu photon_id_sets", auau.get("photon_id_sets"), EXPECTED_AUAU_ROW)
    if auau.get("photon_id_sets") == [["newPPG12", "auauEtFineCent7BDT", "auauBDTComplement"]]:
        errors.append("AuAu R_AA row still uses broad auauBDTComplement; use bounded auauBDTSideband for THE-79")
    assert_equal(errors, "AuAu model dir", str(auau.get("auau_tight_bdt_expanded_model_dir")), str(EXPECTED_MODEL_ROOT))
    assert_equal(errors, "AuAu BDT product", auau.get("auau_tight_bdt_etFineCent7_product"), "base14_perEtCent7")
    assert_equal(errors, "AuAu score label", auau.get("auau_tight_bdt_score_label"), "score_base14_perEtCent7")
    assert_equal(errors, "AuAu WP entries", auau.get("auau_tight_bdt_working_point_entries"), [EXPECTED_WP])
    assert_equal(errors, "AuAu non-tight sideband mode", auau.get("auau_nontight_bdt_sideband_mode"), EXPECTED_SIDE_MODE)
    try:
        if not math.isclose(float(auau.get("auau_nontight_bdt_relative_min_offset")), EXPECTED_SIDE_MIN_OFFSET):
            errors.append(
                "AuAu non-tight relative min offset: "
                f"got {auau.get('auau_nontight_bdt_relative_min_offset')!r}, expected {EXPECTED_SIDE_MIN_OFFSET}"
            )
        if not math.isclose(float(auau.get("auau_nontight_bdt_relative_max_offset")), EXPECTED_SIDE_MAX_OFFSET):
            errors.append(
                "AuAu non-tight relative max offset: "
                f"got {auau.get('auau_nontight_bdt_relative_max_offset')!r}, expected {EXPECTED_SIDE_MAX_OFFSET}"
            )
    except (TypeError, ValueError):
        errors.append("AuAu non-tight relative sideband offsets are not numeric")
    assert_equal(errors, "AuAu etfine pT edges", as_numbers(auau.get("auau_tight_bdt_etfine_pt_bin_edges")), [float(x) for x in EXPECTED_PT_EDGES])
    assert_equal(errors, "AuAu cent7 edges", as_numbers(auau.get("auau_tight_bdt_cent7_edges")), [float(x) for x in EXPECTED_CENT7])
    assert_equal(errors, "AuAu final pT bins", as_numbers(auau.get("unfold_reco_photon_pt_bins")), FINAL_PT_BINS)
    for feature_key in ["auau_tight_bdt_features", "auau_tight_bdt_centDep_features", "auau_tight_bdt_centAsFeatBase3x3_features"]:
        assert_equal(errors, f"AuAu {feature_key}", auau.get(feature_key), EXPECTED_AUAU_FEATURES)
    try:
        if not (math.isclose(float(auau.get("auau_tight_bdt_apply_pt_min")), 8.0) and math.isclose(float(auau.get("auau_tight_bdt_apply_pt_max")), 40.0)):
            errors.append("AuAu apply pT range is not 8-40 GeV")
    except (TypeError, ValueError):
        errors.append("AuAu apply pT range is not numeric")

    manifest = {
        "pp_yaml": str(pp_yaml),
        "auau_yaml": str(auau_yaml),
        "pp_path": "PPG12 newPPG12 tight/sideband baseV3E",
        "auau_model_product": "base14_perEtCent7",
        "auau_model_count_expected": 98,
        "auau_score_label": "score_base14_perEtCent7",
        "auau_wp": EXPECTED_WP,
        "auau_sideband": {
            "mode": EXPECTED_SIDE_MODE,
            "relative_min_offset": EXPECTED_SIDE_MIN_OFFSET,
            "relative_max_offset": EXPECTED_SIDE_MAX_OFFSET,
        },
        "final_pt_bins": FINAL_PT_BINS,
        "routing_pt_edges": EXPECTED_PT_EDGES,
        "centrality_edges": EXPECTED_CENT7,
    }
    return errors, warnings, manifest


def validate_model_root(model_root: Path) -> tuple[list[str], dict[str, Any]]:
    errors: list[str] = []
    registry_path = model_root / "model_registry.json"
    report: dict[str, Any] = {"model_root": str(model_root), "registry": str(registry_path)}
    if not registry_path.exists():
        return [f"model registry missing: {registry_path}"], report
    registry = json.loads(registry_path.read_text())
    report["registry_status"] = registry.get("status")
    report["trained_model_count"] = registry.get("trained_model_count")
    models = registry.get("models") or registry.get("trained_models") or []
    if registry.get("status") != "READY":
        errors.append(f"model registry status is not READY: {registry.get('status')}")
    base14 = [m for m in models if m.get("product") == "base14_perEtCent7"]
    report["base14_perEtCent7_records"] = len(base14)
    if len(base14) != 98:
        errors.append(f"base14_perEtCent7 registry records = {len(base14)}, expected 98")
    bad_features = [
        m.get("model_id")
        for m in base14
        if list(m.get("features") or []) != EXPECTED_AUAU_FEATURES
    ]
    if bad_features:
        errors.append(f"base14_perEtCent7 models with wrong features: {bad_features[:5]}")
    missing = []
    for plo, phi in zip(EXPECTED_PT_EDGES[:-1], EXPECTED_PT_EDGES[1:]):
        for clo, chi in zip(EXPECTED_CENT7[:-1], EXPECTED_CENT7[1:]):
            path = model_root / (
                f"auau_tight_bdt_base14_perEtCent7_pt_{plo:03d}_{phi:03d}_"
                f"cent_{clo:03d}_{chi:03d}_tmva.root"
            )
            if not path.exists():
                missing.append(str(path))
    report["base14_perEtCent7_tmva_files_expected"] = 98
    report["base14_perEtCent7_tmva_files_missing"] = len(missing)
    if missing:
        errors.append(f"missing base14_perEtCent7 TMVA files: {missing[:5]}")
    return errors, report


def validate_taa_table(path: Path | None, allow_pending: bool) -> tuple[list[str], dict[str, Any]]:
    if path is None:
        status = "TAA_PENDING"
        msg = "No official TAA table supplied; normalized R_AA plotting must remain blocked."
        return ([] if allow_pending else [msg]), {"taa_status": status, "message": msg}
    errors: list[str] = []
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        missing_cols = [c for c in REQUIRED_TAA_COLUMNS if c not in (reader.fieldnames or [])]
        if missing_cols:
            errors.append(f"TAA table missing required columns: {missing_cols}")
            return errors, {"taa_status": "INVALID", "path": str(path)}
        rows = list(reader)
    if not rows:
        errors.append("TAA table has no rows")
    for idx, row in enumerate(rows, start=1):
        for col in ["cent_low", "cent_high", "N_coll", "TAA_mb_inv", "TAA_unc", "sigmaNN_mb"]:
            try:
                value = float(row[col])
            except ValueError:
                errors.append(f"TAA row {idx} column {col} is not numeric: {row[col]!r}")
                continue
            if col in {"N_coll", "TAA_mb_inv", "sigmaNN_mb"} and value <= 0:
                errors.append(f"TAA row {idx} column {col} must be positive: {value}")
    return errors, {"taa_status": "READY" if not errors else "INVALID", "path": str(path), "rows": len(rows)}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pp-yaml", required=True, type=Path)
    parser.add_argument("--auau-yaml", required=True, type=Path)
    parser.add_argument("--model-root", type=Path, default=None)
    parser.add_argument("--taa-table", type=Path, default=None)
    parser.add_argument("--allow-taa-pending", action="store_true")
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()

    errors, warnings, manifest = validate_yaml(args.pp_yaml, args.auau_yaml)
    if args.model_root is not None:
        model_errors, model_report = validate_model_root(args.model_root)
        errors.extend(model_errors)
        manifest["model_report"] = model_report
    taa_errors, taa_report = validate_taa_table(args.taa_table, args.allow_taa_pending)
    errors.extend(taa_errors)
    manifest["taa_report"] = taa_report
    manifest["status"] = "PASS" if not errors else "FAIL"
    manifest["errors"] = errors
    manifest["warnings"] = warnings

    text = json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text)
    print(text, end="")
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
