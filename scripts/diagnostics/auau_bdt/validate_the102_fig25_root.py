#!/usr/bin/env python3
"""Validate THE-102 joint background-correlation surfaces."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import uproot


CENTRALITY_TOKENS = ("cent_0_20", "cent_20_50", "cent_50_80")
AXES = ("e11e33", "bdtScore")
BASE = "h2_auauFig25_{axis}_vs_Eiso_background_pT15to35"


def clean_keys(root_file: uproot.ReadOnlyDirectory) -> list[str]:
    return sorted(str(key).split(";", 1)[0] for key in root_file.keys(recursive=True))


def require(condition: bool, message: str, failures: list[str]) -> None:
    if not condition:
        failures.append(message)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--allow-empty-centrality", action="store_true")
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()

    failures: list[str] = []
    require(args.root.is_file(), f"missing ROOT file: {args.root}", failures)
    require(args.root.stat().st_size >= 50_000 if args.root.is_file() else False,
            f"tiny ROOT file: {args.root}", failures)
    if failures:
        print(json.dumps({"status": "FAIL", "failures": failures}, indent=2))
        return 1

    root_file = uproot.open(args.root)
    keys = clean_keys(root_file)
    report_surfaces: dict[str, dict[str, object]] = {}

    for axis in AXES:
        token = BASE.format(axis=axis)
        matches = [key for key in keys if token in key]
        require(len(matches) == 3,
                f"{axis}: expected 3 centrality surfaces, found {len(matches)}", failures)
        for cent in CENTRALITY_TOKENS:
            cent_matches = [key for key in matches if cent in key]
            require(len(cent_matches) == 1,
                    f"{axis} {cent}: expected exactly one surface, found {len(cent_matches)}", failures)
            if len(cent_matches) != 1:
                continue
            key = cent_matches[0]
            hist = root_file[key]
            values = np.asarray(hist.values(flow=False), dtype=np.float64)
            variances = hist.variances(flow=False)
            x_edges = np.asarray(hist.axis(0).edges(), dtype=np.float64)
            y_edges = np.asarray(hist.axis(1).edges(), dtype=np.float64)
            integral = float(np.sum(values))
            require(values.shape == (50, 160),
                    f"{key}: wrong shape {values.shape}, expected (50, 160)", failures)
            require(np.all(np.isfinite(values)), f"{key}: non-finite contents", failures)
            require(np.all(values >= -1e-12), f"{key}: negative weighted contents", failures)
            require(np.allclose(x_edges[[0, -1]], [0.0, 1.000001], atol=1e-9),
                    f"{key}: wrong x range {x_edges[[0, -1]].tolist()}", failures)
            require(np.allclose(y_edges[[0, -1]], [-20.0, 60.0], atol=1e-9),
                    f"{key}: wrong isolation range {y_edges[[0, -1]].tolist()}", failures)
            require(variances is not None, f"{key}: Sumw2 variances missing", failures)
            if not args.allow_empty_centrality:
                require(integral > 0.0, f"{key}: empty surface", failures)
            report_surfaces[f"{axis}:{cent}"] = {
                "key": key,
                "shape": list(values.shape),
                "integral": integral,
            }

    classifier_keys = [key for key in keys if "pmtDiag" in key or "embeddedMinBias" in key]
    require(not classifier_keys,
            "embedded MinimumBiasClassifier diagnostic/filter objects unexpectedly present", failures)

    report = {
        "schema": "THE102_AUAU_FIG25_ROOT_VALIDATION_V1",
        "status": "PASS" if not failures else "FAIL",
        "root": str(args.root),
        "allow_empty_centrality": args.allow_empty_centrality,
        "surfaces": report_surfaces,
        "classifier_object_count": len(classifier_keys),
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
