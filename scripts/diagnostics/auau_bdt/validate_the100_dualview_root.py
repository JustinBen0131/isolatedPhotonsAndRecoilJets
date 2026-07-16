#!/usr/bin/env python3
"""Validate paired THE-100 bounded/complement RecoilJets outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import uproot


CENTRALITY_TOKENS = ("cent_0_20", "cent_20_50", "cent_50_80")
PT_EDGES = np.asarray([15.0, 17.0, 19.0, 21.0, 23.0, 26.0, 35.0])


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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bounded", type=Path, required=True)
    parser.add_argument("--complement", type=Path, required=True)
    parser.add_argument("--sample-kind", choices=("data", "signal", "inclusive"), required=True)
    parser.add_argument(
        "--allow-empty-tight",
        action="store_true",
        help="Allow a statistically sparse canary with no booked tight histogram; never use for merged production QA.",
    )
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()

    failures: list[str] = []
    bounded, bounded_keys = summarize(args.bounded)
    complement, complement_keys = summarize(args.complement)

    for label, root_file, keys in (
        ("bounded", bounded, bounded_keys),
        ("complement", complement, complement_keys),
    ):
        require(any(key.endswith("h_centrality") for key in keys),
                f"{label}: h_centrality missing", failures)
        for cent in CENTRALITY_TOKENS:
            require(any(cent in key for key in keys),
                    f"{label}: no content for {cent}", failures)
        require(bool(matching(keys, "h_xJpurityLead_isIsolated_notTight")),
                f"{label}: Region-C leading-photon counter missing", failures)
        require(bool(matching(keys, "sidebandC")),
                f"{label}: sideband-C recoil output missing", failures)

    all_surfaces = matching(
        bounded_keys, "h3_auauDualView_scoreMinusT80_vs_Eiso_vs_pT_all"
    )
    require(bool(all_surfaces), "bounded: dual-view all-candidate surface missing", failures)
    surface_integral = 0.0
    for key in all_surfaces:
        array = values(bounded, key)
        require(np.all(np.isfinite(array)), f"bounded: non-finite surface bins in {key}", failures)
        require(np.all(array >= -1e-12), f"bounded: negative surface bins in {key}", failures)
        surface_integral += float(np.sum(array))
        z_edges = np.asarray(bounded[key].axis(2).edges(), dtype=np.float64)
        require(np.array_equal(z_edges, PT_EDGES),
                f"bounded: wrong photon-pT edges in {key}: {z_edges.tolist()}", failures)
    require(surface_integral > 0.0, "bounded: dual-view surface is empty", failures)

    signal_surface_integral = 0.0
    signal_surfaces = matching(
        bounded_keys, "h3_auauDualView_scoreMinusT80_vs_Eiso_vs_pT_truthSignal"
    )
    for key in signal_surfaces:
        signal_surface_integral += float(np.sum(values(bounded, key)))
    if args.sample_kind == "signal":
        require(signal_surface_integral > 0.0,
                "bounded signal: truth-signal surface is empty or missing", failures)
        require(signal_surface_integral <= surface_integral + 1e-9,
                "bounded signal: truth-signal exceeds all-candidate surface", failures)

    common_tight = sorted(
        set(matching(bounded_keys, "h_Eiso_tight"))
        & set(matching(complement_keys, "h_Eiso_tight"))
    )
    require(
        bool(common_tight) or args.allow_empty_tight,
        "no common tight-isolation histograms",
        failures,
    )
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
        "schema": "THE100_AUAU_DUALVIEW_ROOT_VALIDATION_V1",
        "status": "PASS" if not failures else "FAIL",
        "sample_kind": args.sample_kind,
        "bounded": str(args.bounded),
        "complement": str(args.complement),
        "bounded_keys": len(bounded_keys),
        "complement_keys": len(complement_keys),
        "all_surface_count": len(all_surfaces),
        "all_surface_integral": surface_integral,
        "truth_signal_surface_count": len(signal_surfaces),
        "truth_signal_surface_integral": signal_surface_integral,
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
