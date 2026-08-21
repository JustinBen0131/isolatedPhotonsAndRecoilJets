#!/usr/bin/env python3
"""Validate paired THE-112 bounded/complement canary ROOT outputs.

The validator is intentionally stricter than the historical THE-100 check:
it requires versioned, cone-separated score/isolation surfaces and proves that
each candidate is represented once in canonical R=0.4 and once in the R=0.3
robustness view.  It never selects the nominal sideband; it only validates the
transport needed by the separate offline ranking step.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
from typing import Any

import numpy as np
import uproot
import yaml


CENTRALITY_TOKENS = ("cent_0_20", "cent_20_50", "cent_50_80")
PT_EDGES = np.asarray([15.0, 17.0, 19.0, 21.0, 23.0, 26.0, 35.0])
SCORE_STEP = 0.01
V2_TOKEN = "h3_auauSidebandScanV2_scoreMinusT80_vs_EisoMinusCut_vs_pT_"
# Reconstructed Au+Au isolation is sliding-only.  R=0.4 is canonical and
# R=0.3 is a robustness view; fixed reconstructed-isolation views are absent
# unless Justin separately authorizes such a diagnostic.
ISO_VIEWS = (
    "isoR40_isSliding",
    "isoR30_isSliding",
)
EXPECTED_INTERNAL_ISO_VIEWS = (
    "isoR40_isSliding:0.40:true:0.0,"
    "isoR30_isSliding:0.30:true:0.0"
)
EXPECTED_CONE_RADII = (0.40, 0.30)


def clean_keys(root_file: uproot.ReadOnlyDirectory) -> list[str]:
    return sorted(str(key).split(";", 1)[0] for key in root_file.keys(recursive=True))


def values(root_file: uproot.ReadOnlyDirectory, key: str) -> np.ndarray:
    return np.asarray(root_file[key].values(flow=False), dtype=np.float64)


def matching(keys: list[str], token: str) -> list[str]:
    return [key for key in keys if token in key]


def require(condition: bool, message: str, failures: list[str]) -> None:
    if not condition:
        failures.append(message)


def summarize(path: Path) -> tuple[uproot.ReadOnlyDirectory, list[str]]:
    if not path.is_file() or path.stat().st_size < 50_000:
        raise RuntimeError(f"missing or tiny ROOT file: {path}")
    root_file = uproot.open(path)
    keys = clean_keys(root_file)
    if not keys:
        raise RuntimeError(f"ROOT file has no keys: {path}")
    return root_file, keys


def expected_surface_basename(category: str, cent: str, view: str) -> str:
    return f"{V2_TOKEN}{category}_{view}_{cent}"


def key_scope_and_basename(key: str) -> tuple[str, str]:
    scope, separator, basename = key.rpartition("/")
    return (scope if separator else "", basename if separator else key)


def v2_scopes(keys: list[str]) -> tuple[str, ...]:
    return tuple(sorted({key_scope_and_basename(key)[0] for key in matching(keys, V2_TOKEN)}))


def surface_map(keys: list[str], category: str) -> dict[tuple[str, str, str], str]:
    found: dict[tuple[str, str, str], str] = {}
    expected = {
        expected_surface_basename(category, cent, view): (cent, view)
        for cent in CENTRALITY_TOKENS
        for view in ISO_VIEWS
    }
    for key in matching(keys, V2_TOKEN + category):
        scope, basename = key_scope_and_basename(key)
        if basename in expected:
            cent, view = expected[basename]
            found[(scope, cent, view)] = key
    return found


def validate_v2_inventory(
    keys: list[str],
    sample_kind: str,
    failures: list[str],
) -> dict[str, Any]:
    """Require one and only one histogram for every allowed V2 namespace."""

    categories = (
        ("all",)
        if sample_kind == "data"
        else ("all", "truthPrompt", "truthSignal", "truthBackground")
    )
    expected = {
        expected_surface_basename(category, cent, view)
        for category in categories
        for cent in CENTRALITY_TOKENS
        for view in ISO_VIEWS
    }
    actual_keys = matching(keys, V2_TOKEN)
    scoped_counts: dict[str, Counter[str]] = {}
    for key in actual_keys:
        scope, basename = key_scope_and_basename(key)
        scoped_counts.setdefault(scope, Counter())[basename] += 1

    require(bool(scoped_counts), "V2 inventory: no scan namespace found", failures)
    missing: dict[str, list[str]] = {}
    extra: dict[str, list[str]] = {}
    duplicates: dict[str, dict[str, int]] = {}
    for scope, counts in sorted(scoped_counts.items()):
        scope_label = scope or "<root>"
        scope_missing = sorted(expected - set(counts))
        scope_extra = sorted(set(counts) - expected)
        scope_duplicates = {name: count for name, count in sorted(counts.items()) if count != 1}
        if scope_missing:
            missing[scope_label] = scope_missing
        if scope_extra:
            extra[scope_label] = scope_extra
        if scope_duplicates:
            duplicates[scope_label] = scope_duplicates

    require(not missing, f"V2 inventory: missing exact namespaces by trigger scope: {missing}", failures)
    require(not extra, f"V2 inventory: forbidden extra namespaces by trigger scope: {extra}", failures)
    require(not duplicates, f"V2 inventory: duplicate namespaces by trigger scope: {duplicates}", failures)

    return {
        "scope_count": len(scoped_counts),
        "scopes": [scope or "<root>" for scope in sorted(scoped_counts)],
        "expected_count_per_scope": len(expected),
        "expected_count_total": len(expected) * len(scoped_counts),
        "observed_count": len(actual_keys),
        "categories": list(categories),
        "views": list(ISO_VIEWS),
        "missing": missing,
        "extra": extra,
        "duplicates": duplicates,
    }


def validate_no_fixed_namespaces(
    keys: list[str],
    label: str,
    failures: list[str],
) -> dict[str, Any]:
    fixed_keys = sorted(key for key in keys if "_fixedIso" in key)
    require(
        not fixed_keys,
        f"{label}: fixed reconstructed-isolation histogram namespaces are forbidden: {fixed_keys}",
        failures,
    )
    return {"fixed_namespace_count": len(fixed_keys), "fixed_namespaces": fixed_keys}


def validate_embedded_config(
    root_file: uproot.ReadOnlyDirectory,
    label: str,
    failures: list[str],
) -> dict[str, Any]:
    """Audit the stamped config stored in a canary ROOT file."""

    if "analysis_config_yaml" not in root_file:
        failures.append(f"{label}: analysis_config_yaml missing")
        return {"present": False, "decoded": False}

    try:
        config_text = str(root_file["analysis_config_yaml"])
        config = yaml.safe_load(config_text)
    except Exception as exc:  # noqa: BLE001 - report corrupt provenance as validation failure
        failures.append(f"{label}: analysis_config_yaml could not be decoded: {exc}")
        return {"present": True, "decoded": False}
    if not isinstance(config, dict):
        failures.append(f"{label}: analysis_config_yaml did not decode to a mapping")
        return {"present": True, "decoded": False}

    sliding = config.get("isSlidingIso")
    sliding_and_fixed = config.get("isSlidingAndFixed")
    fixed_gev = config.get("fixedGeV")
    internal_views = config.get("internal_iso_cone_views")
    cone_r = config.get("coneR")
    isolation_wp = config.get("isolation_wp")
    truth_iso = isolation_wp.get("truthIsoGeV") if isinstance(isolation_wp, dict) else None

    require(sliding is True, f"{label}: isSlidingIso must be true, observed {sliding!r}", failures)
    require(
        sliding_and_fixed is False,
        f"{label}: isSlidingAndFixed must be false, observed {sliding_and_fixed!r}",
        failures,
    )
    fixed_is_zero = (
        isinstance(fixed_gev, (int, float))
        and not isinstance(fixed_gev, bool)
        and np.isclose(float(fixed_gev), 0.0, rtol=0.0, atol=1e-12)
    )
    require(fixed_is_zero, f"{label}: fixedGeV must be the 0 GeV sentinel, observed {fixed_gev!r}", failures)
    require(
        internal_views == EXPECTED_INTERNAL_ISO_VIEWS,
        f"{label}: internal_iso_cone_views mismatch: expected "
        f"{EXPECTED_INTERNAL_ISO_VIEWS!r}, observed {internal_views!r}",
        failures,
    )
    try:
        if isinstance(cone_r, list):
            cone_values = tuple(float(value) for value in cone_r)
        elif isinstance(cone_r, (int, float)) and not isinstance(cone_r, bool):
            # The source campaign YAML carries the ordered [0.4, 0.3] list.
            # The Condor matrix expander stamps the nominal outer lane as the
            # scalar 0.4 while preserving both exact diagnostic views in
            # internal_iso_cone_views.  Accept that stamped representation,
            # but never a scalar robustness-only or fixed-isolation lane.
            cone_values = (float(cone_r),) if internal_views == EXPECTED_INTERNAL_ISO_VIEWS else ()
        else:
            cone_values = ()
    except (TypeError, ValueError):
        cone_values = ()
    cone_contract_ok = (
        len(cone_values) == len(EXPECTED_CONE_RADII)
        and np.allclose(cone_values, EXPECTED_CONE_RADII, rtol=0.0, atol=1e-12)
    ) or (
        len(cone_values) == 1
        and np.isclose(cone_values[0], EXPECTED_CONE_RADII[0], rtol=0.0, atol=1e-12)
        and internal_views == EXPECTED_INTERNAL_ISO_VIEWS
    )
    require(
        cone_contract_ok,
        f"{label}: coneR must be ordered [0.4, 0.3], or stamped nominal 0.4 with the exact "
        f"ordered internal-view contract; observed {cone_r!r}",
        failures,
    )
    truth_is_four = (
        isinstance(truth_iso, (int, float))
        and not isinstance(truth_iso, bool)
        and np.isclose(float(truth_iso), 4.0, rtol=0.0, atol=1e-12)
    )
    require(
        truth_is_four,
        f"{label}: isolation_wp.truthIsoGeV must be 4 GeV, observed {truth_iso!r}",
        failures,
    )

    return {
        "present": True,
        "decoded": True,
        "isSlidingIso": sliding,
        "isSlidingAndFixed": sliding_and_fixed,
        "fixedGeV": fixed_gev,
        "internal_iso_cone_views": internal_views,
        "coneR": cone_r,
        "truthIsoGeV": truth_iso,
    }


def validate_surface_family(
    root_file: uproot.ReadOnlyDirectory,
    keys: list[str],
    category: str,
    failures: list[str],
    require_nonzero: bool,
) -> dict[str, Any]:
    scopes = v2_scopes(keys)
    surfaces = surface_map(keys, category)
    expected = {
        (scope, cent, view)
        for scope in scopes
        for cent in CENTRALITY_TOKENS
        for view in ISO_VIEWS
    }
    missing = sorted(expected - set(surfaces))
    require(not missing, f"{category}: missing cone/centrality surfaces: {missing}", failures)

    integrals: dict[str, float] = {}
    max_pair_delta = 0.0
    for (scope, cent, view), key in sorted(surfaces.items()):
        histogram = root_file[key]
        # Include underflow/overflow in the candidate-weight identity check.
        # The sideband scan uses the sign of the isolation coordinate, so a
        # high-isolation overflow is still a valid non-isolated candidate and
        # must not disappear from transport validation.
        array = np.asarray(histogram.values(flow=True), dtype=np.float64)
        require(np.all(np.isfinite(array)), f"{category}: non-finite bins in {key}", failures)
        require(np.all(array >= -1e-12), f"{category}: negative bins in {key}", failures)
        x_edges = np.asarray(histogram.axis(0).edges(), dtype=np.float64)
        z_edges = np.asarray(histogram.axis(2).edges(), dtype=np.float64)
        require(
            x_edges.size >= 2 and np.allclose(np.diff(x_edges), SCORE_STEP, rtol=0.0, atol=1e-12),
            f"{category}: score axis is not uniformly {SCORE_STEP:g} in {key}",
            failures,
        )
        require(
            np.array_equal(z_edges, PT_EDGES),
            f"{category}: wrong photon-pT edges in {key}: {z_edges.tolist()}",
            failures,
        )
        scope_label = scope or "<root>"
        integrals[f"{scope_label}:{cent}_{view}"] = float(np.sum(array))

    for scope in scopes:
        scope_label = scope or "<root>"
        for cent in CENTRALITY_TOKENS:
            view_integrals = [
                integrals.get(f"{scope_label}:{cent}_{view}", 0.0)
                for view in ISO_VIEWS
            ]
            reference = view_integrals[0]
            for view, integral in zip(ISO_VIEWS[1:], view_integrals[1:]):
                max_pair_delta = max(max_pair_delta, abs(reference - integral))
                require(
                    np.isclose(reference, integral, rtol=0.0, atol=1e-9),
                    f"{category}: internal-view candidate-weight mismatch in "
                    f"{scope_label}/{cent}: {ISO_VIEWS[0]}={reference} vs {view}={integral}",
                    failures,
                )
    total = float(sum(integrals.values()))
    if require_nonzero:
        require(total > 0.0, f"{category}: all V2 surfaces are empty", failures)
    return {
        "surface_count": len(surfaces),
        "trigger_scopes": [scope or "<root>" for scope in scopes],
        "integrals": integrals,
        "max_internal_view_integral_delta": max_pair_delta,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bounded", type=Path, required=True)
    parser.add_argument("--complement", type=Path, required=True)
    parser.add_argument("--sample-kind", choices=("data", "signal", "inclusive"), required=True)
    parser.add_argument(
        "--allow-empty-tight",
        action="store_true",
        help="Allow sparse canaries with no common tight histogram; never use for merged production QA.",
    )
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()

    failures: list[str] = []
    bounded, bounded_keys = summarize(args.bounded)
    complement, complement_keys = summarize(args.complement)

    config_reports = {
        "bounded": validate_embedded_config(bounded, "bounded", failures),
        "complement": validate_embedded_config(complement, "complement", failures),
    }
    namespace_reports = {
        "bounded": validate_no_fixed_namespaces(bounded_keys, "bounded", failures),
        "complement": validate_no_fixed_namespaces(complement_keys, "complement", failures),
    }

    for label, keys in (("bounded", bounded_keys), ("complement", complement_keys)):
        require(any(key.endswith("h_centrality") for key in keys), f"{label}: h_centrality missing", failures)
        for cent in CENTRALITY_TOKENS:
            require(any(cent in key for key in keys), f"{label}: no content for {cent}", failures)
        require(bool(matching(keys, "h_xJpurityLead_isIsolated_notTight")),
                f"{label}: Region-C leading-photon counter missing", failures)
        require(bool(matching(keys, "sidebandC")),
                f"{label}: sideband-C recoil output missing", failures)

    require(
        not matching(bounded_keys, "h3_auauDualView_scoreMinusT80_vs_Eiso_vs_pT"),
        "bounded: obsolete mixed-cone THE-100 surface is present",
        failures,
    )
    v2_inventory = validate_v2_inventory(bounded_keys, args.sample_kind, failures)
    require(
        not matching(complement_keys, V2_TOKEN),
        "complement: V2 scan surfaces must be owned only by the bounded scan lane",
        failures,
    )
    surface_reports: dict[str, Any] = {
        "all": validate_surface_family(bounded, bounded_keys, "all", failures, True)
    }
    if args.sample_kind != "data":
        surface_reports["truthPrompt"] = validate_surface_family(
            bounded, bounded_keys, "truthPrompt", failures, args.sample_kind == "signal"
        )
        surface_reports["truthSignal"] = validate_surface_family(
            bounded, bounded_keys, "truthSignal", failures, args.sample_kind == "signal"
        )
        surface_reports["truthBackground"] = validate_surface_family(
            bounded,
            bounded_keys,
            "truthBackground",
            failures,
            args.sample_kind == "inclusive",
        )

    common_tight = sorted(
        set(matching(bounded_keys, "h_Eiso_tight"))
        & set(matching(complement_keys, "h_Eiso_tight"))
    )
    require(bool(common_tight) or args.allow_empty_tight, "no common tight-isolation histograms", failures)
    max_tight_delta = 0.0
    for key in common_tight:
        left = values(bounded, key)
        right = values(complement, key)
        require(left.shape == right.shape, f"tight shape mismatch: {key}", failures)
        if left.shape == right.shape:
            delta = float(np.max(np.abs(left - right))) if left.size else 0.0
            max_tight_delta = max(max_tight_delta, delta)
            require(np.allclose(left, right, rtol=0.0, atol=1e-9),
                    f"tight output differs between views: {key} max_delta={delta}", failures)

    common_nontight = sorted(
        set(matching(bounded_keys, "h_Eiso_nonTight"))
        & set(matching(complement_keys, "h_Eiso_nonTight"))
    )
    require(bool(common_nontight), "no common non-tight isolation histograms", failures)
    subset_violations = 0
    for key in common_nontight:
        left = values(bounded, key)
        right = values(complement, key)
        if left.shape != right.shape:
            failures.append(f"non-tight shape mismatch: {key}")
            continue
        subset_violations += int(np.count_nonzero(left > right + 1e-9))
    require(subset_violations == 0,
            f"bounded non-tight is not a subset of complement: {subset_violations} bins", failures)

    report: dict[str, Any] = {
        "schema": "THE112_AUAU_SIDEBAND_CANARY_VALIDATION_V1",
        "status": "PASS" if not failures else "FAIL",
        "sample_kind": args.sample_kind,
        "bounded": str(args.bounded),
        "complement": str(args.complement),
        "bounded_keys": len(bounded_keys),
        "complement_keys": len(complement_keys),
        "embedded_config": config_reports,
        "reconstructed_isolation_namespaces": namespace_reports,
        "v2_inventory": v2_inventory,
        "surface_reports": surface_reports,
        "common_tight_histograms": len(common_tight),
        "allow_empty_tight": args.allow_empty_tight,
        "max_tight_bin_delta": max_tight_delta,
        "common_nontight_histograms": len(common_nontight),
        "bounded_subset_violations": subset_violations,
        "failures": failures,
    }
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(payload)
    print(payload, end="")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
