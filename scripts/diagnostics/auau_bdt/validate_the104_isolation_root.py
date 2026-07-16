#!/usr/bin/env python3
"""Validate THE-104 truth-matched AuAu isolation histogram contracts."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import uproot


PT_EDGES = (15, 17, 19, 21, 23, 26, 30, 35)
CENT_EDGES = tuple(range(0, 85, 5))
CONES = ("R30", "R40")
PATTERN = re.compile(
    r"^SIM/h_EisoReco_truthSigMatched_iso(?P<cone>R30|R40)_"
    r"pT_(?P<ptlo>\d+)_(?P<pthi>\d+)_cent_(?P<clo>\d+)_(?P<chi>\d+)$"
)


def clean_keys(root_file: uproot.ReadOnlyDirectory) -> list[str]:
    return sorted(str(key).split(";", 1)[0] for key in root_file.keys(recursive=True))


def expected_keys() -> set[str]:
    return {
        f"SIM/h_EisoReco_truthSigMatched_iso{cone}_pT_{ptlo}_{pthi}_cent_{clo}_{chi}"
        for cone in CONES
        for ptlo, pthi in zip(PT_EDGES[:-1], PT_EDGES[1:])
        for clo, chi in zip(CENT_EDGES[:-1], CENT_EDGES[1:])
    }


def require(condition: bool, message: str, failures: list[str]) -> None:
    if not condition:
        failures.append(message)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument(
        "--allow-sparse",
        action="store_true",
        help="Canary mode: require valid populated objects for both cones, not the full matrix.",
    )
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()

    failures: list[str] = []
    require(args.root.is_file(), f"missing ROOT file: {args.root}", failures)
    require(
        args.root.stat().st_size >= 50_000 if args.root.is_file() else False,
        f"tiny ROOT file: {args.root}",
        failures,
    )
    if failures:
        print(json.dumps({"status": "FAIL", "failures": failures}, indent=2))
        return 1

    root_file = uproot.open(args.root)
    keys = clean_keys(root_file)
    require("SIM" in {key.split("/", 1)[0] for key in keys}, "missing SIM directory", failures)

    observed = {key for key in keys if PATTERN.match(key)}
    expected = expected_keys()
    unexpected = sorted(observed - expected)
    missing = sorted(expected - observed)
    require(not unexpected, f"unexpected target histogram names: {unexpected[:5]}", failures)
    if args.allow_sparse:
        for cone in CONES:
            require(
                any(f"_iso{cone}_" in key for key in observed),
                f"canary has no {cone} truth-matched isolation object",
                failures,
            )
    else:
        require(not missing, f"missing {len(missing)} of {len(expected)} target histograms", failures)

    cone_integrals = {cone: 0.0 for cone in CONES}
    populated = 0
    for key in sorted(observed):
        hist = root_file[key]
        values = np.asarray(hist.values(flow=False), dtype=np.float64)
        variances = hist.variances(flow=False)
        edges = np.asarray(hist.axis().edges(), dtype=np.float64)
        require(values.shape == (700,), f"{key}: wrong shape {values.shape}", failures)
        require(np.all(np.isfinite(values)), f"{key}: non-finite contents", failures)
        require(np.all(values >= -1e-12), f"{key}: negative weighted contents", failures)
        require(variances is not None, f"{key}: Sumw2 variances missing", failures)
        if variances is not None:
            variances = np.asarray(variances, dtype=np.float64)
            require(np.all(np.isfinite(variances)), f"{key}: non-finite variances", failures)
            require(np.all(variances >= -1e-12), f"{key}: negative variances", failures)
        require(
            np.allclose(edges[[0, -1]], [-20.0, 50.0], atol=1e-12),
            f"{key}: wrong isolation range {edges[[0, -1]].tolist()}",
            failures,
        )
        integral = float(np.sum(values))
        if integral > 0.0:
            populated += 1
        match = PATTERN.match(key)
        assert match is not None
        cone_integrals[match.group("cone")] += integral

    for cone in CONES:
        require(cone_integrals[cone] > 0.0, f"{cone}: zero aggregate population", failures)

    report = {
        "schema": "THE104_AUAU_CANONICAL_MINBIAS_ISOLATION_ROOT_VALIDATION_V1",
        "status": "PASS" if not failures else "FAIL",
        "root": str(args.root),
        "allow_sparse": args.allow_sparse,
        "expected_histograms": len(expected),
        "observed_histograms": len(observed),
        "populated_histograms": populated,
        "missing_histograms": len(missing),
        "cone_integrals": cone_integrals,
        "physics_contract": {
            "sample": "AuAu embedded Photon12+Photon20",
            "minimum_bias_classifier": "required before physics histogram filling",
            "photon_pt_edges_gev": list(PT_EDGES),
            "centrality_edges_percent": list(CENT_EDGES),
            "cones": list(CONES),
        },
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
