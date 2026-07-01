#!/usr/bin/env python3
"""Preflight pp/AuAu final xJgamma production-axis contract."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


FINAL_PHOTON_PT_BINS = [15, 17, 19, 21, 23, 26, 35]
FINAL_UNFOLD_RECO_PHOTON_PT_BINS = [10, 15, 17, 19, 21, 23, 26, 35, 40]
FINAL_UNFOLD_TRUTH_PHOTON_PT_BINS = [5, 10, 15, 17, 19, 21, 23, 26, 35, 40]
PPG12_YIELD_RECO_BINS = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36]
OLD_AUAU_RESPONSE_BINS = [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 35, 40]


def strip_comment(line: str) -> str:
    in_quote = False
    quote_char = ""
    out: list[str] = []
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
    parts: list[Any] = []
    token: list[str] = []
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
    data: dict[str, Any] = {}
    lines = path.read_text().splitlines()
    for raw in lines:
        line = strip_comment(raw)
        if not line.strip() or line.startswith(" ") or ":" not in line:
            continue
        key, rhs = line.split(":", 1)
        rhs = rhs.strip()
        if rhs.startswith("["):
            data[key.strip()] = parse_inline_list(rhs)
        elif rhs:
            data[key.strip()] = parse_scalar(rhs)
    return data


def as_numbers(values: Any) -> list[float]:
    if not isinstance(values, list):
        return []
    out: list[float] = []
    for value in values:
        try:
            out.append(float(value))
        except (TypeError, ValueError):
            return []
    return out


def normalize(values: list[float]) -> list[float | int]:
    out: list[float | int] = []
    for value in values:
        if abs(value - round(value)) < 1e-9:
            out.append(int(round(value)))
        else:
            out.append(value)
    return out


def env_truthy(value: str | None) -> bool:
    return value is not None and value.strip().lower() in {"1", "true", "yes", "on"}


def validate_config(label: str, path: Path, data: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    gamma_bins = normalize(as_numbers(data.get("jes3_photon_pt_bins")))
    reco_bins = normalize(as_numbers(data.get("unfold_reco_photon_pt_bins")))
    truth_bins = normalize(as_numbers(data.get("unfold_truth_photon_pt_bins")))

    if gamma_bins != FINAL_PHOTON_PT_BINS:
        errors.append(
            f"{label} {path}: jes3_photon_pt_bins={gamma_bins}, expected {FINAL_PHOTON_PT_BINS}"
        )
    if reco_bins != FINAL_UNFOLD_RECO_PHOTON_PT_BINS:
        errors.append(
            f"{label} {path}: unfold_reco_photon_pt_bins={reco_bins}, expected {FINAL_UNFOLD_RECO_PHOTON_PT_BINS}"
        )
    if truth_bins != FINAL_UNFOLD_TRUTH_PHOTON_PT_BINS:
        errors.append(
            f"{label} {path}: unfold_truth_photon_pt_bins={truth_bins}, expected {FINAL_UNFOLD_TRUTH_PHOTON_PT_BINS}"
        )
    if reco_bins == PPG12_YIELD_RECO_BINS:
        errors.append(
            f"{label} {path}: final xJgamma reco axis is the PPG12 photon-yield axis; decouple RJ_PPG12_PHOTON_YIELD from unfolding bins"
        )
    if reco_bins == OLD_AUAU_RESPONSE_BINS:
        errors.append(
            f"{label} {path}: final xJgamma reco axis is the old AuAu response axis, not the matched pp/AuAu final-analysis axis"
        )
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pp-yaml", required=True, type=Path)
    parser.add_argument("--auau-yaml", required=True, type=Path)
    parser.add_argument(
        "--env",
        action="append",
        default=[],
        help="Environment assignment to validate, e.g. RJ_PPG12_PHOTON_YIELD_APPLY_BINNING=0",
    )
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()

    pp = parse_simple_yaml(args.pp_yaml)
    auau = parse_simple_yaml(args.auau_yaml)
    errors = []
    errors.extend(validate_config("pp", args.pp_yaml, pp))
    errors.extend(validate_config("AuAu", args.auau_yaml, auau))

    env = {}
    for assignment in args.env:
        if "=" not in assignment:
            errors.append(f"Malformed --env assignment: {assignment!r}")
            continue
        key, value = assignment.split("=", 1)
        env[key.strip()] = value.strip()

    if env_truthy(env.get("RJ_PPG12_PHOTON_YIELD_APPLY_BINNING")):
        errors.append(
            "RJ_PPG12_PHOTON_YIELD_APPLY_BINNING=1 is forbidden for final xJgamma production; "
            "PPG12 photon-yield mode may be on, but it must not override xJgamma axes."
        )

    manifest = {
        "pp_yaml": str(args.pp_yaml),
        "auau_yaml": str(args.auau_yaml),
        "final_photon_pt_bins": FINAL_PHOTON_PT_BINS,
        "unfold_reco_photon_pt_bins": FINAL_UNFOLD_RECO_PHOTON_PT_BINS,
        "unfold_truth_photon_pt_bins": FINAL_UNFOLD_TRUTH_PHOTON_PT_BINS,
        "ppg12_yield_apply_binning_allowed": False,
        "status": "PASS" if not errors else "FAIL",
        "errors": errors,
    }
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(manifest, indent=2) + "\n")

    print(json.dumps(manifest, indent=2))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
