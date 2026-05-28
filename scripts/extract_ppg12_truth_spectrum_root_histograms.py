#!/usr/bin/env python3
"""Export compact CSV/JSON from PPG12 truth-spectrum ROOT histograms."""

from __future__ import annotations

import csv
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import uproot


BASE = Path("/sphenix/user/shuhangli/ppg12/efficiencytool/results")


@dataclass(frozen=True)
class Sample:
    group: str
    sample: str
    xsec_pb: float
    win_lo: float
    win_hi: float
    color: str
    hist_name: str


PHOTON_SAMPLES = [
    Sample("photon", "photon5", 146359.3, 0.0, 14.0, "#E7298A", "h_max_truth_photon_pt_filtered"),
    Sample("photon", "photon10", 6944.675, 14.0, 22.0, "#33A02C", "h_max_truth_photon_pt_filtered"),
    Sample("photon", "photon20", 130.4461, 22.0, 200.0, "#1F78B4", "h_max_truth_photon_pt_filtered"),
]

JET_SAMPLES = [
    Sample("jet", "jet8", 1.3013e7, 9.0, 14.0, "#E7298A", "h_max_truth_jet_pt_filtered"),
    Sample("jet", "jet12", 1.4903e6, 14.0, 21.0, "#33A02C", "h_max_truth_jet_pt_filtered"),
    Sample("jet", "jet20", 6.2623e4, 21.0, 32.0, "#1F78B4", "h_max_truth_jet_pt_filtered"),
    Sample("jet", "jet30", 2.5298e3, 32.0, 42.0, "#FF7F00", "h_max_truth_jet_pt_filtered"),
    Sample("jet", "jet40", 1.3553e2, 42.0, 200.0, "#E7298A", "h_max_truth_jet_pt_filtered"),
]


def read_scalar_string(f: uproot.ReadOnlyDirectory, key: str) -> str | None:
    try:
        obj = f[key]
    except Exception:
        return None
    try:
        return str(obj.member("fTitle") or obj.member("fName"))
    except Exception:
        return str(obj)


def emit_sample(spec: Sample, writer: csv.DictWriter) -> dict[str, object]:
    path = BASE / f"truth_spectrum_{spec.sample}.root"
    with uproot.open(path) as f:
        h = f[spec.hist_name]
        values, edges = h.to_numpy(flow=False)
        variances = h.variances(flow=False)
        if variances is None:
            errors = np.sqrt(np.clip(values, 0, None))
        else:
            errors = np.sqrt(np.clip(variances, 0, None))
        nevt_read = read_scalar_string(f, "nevt_read")
        per_entry_weight = read_scalar_string(f, "per_entry_weight_pb_per_GeV")
        bin_width = read_scalar_string(f, "bin_width_GeV")

    centers = 0.5 * (edges[:-1] + edges[1:])
    for lo, hi, center, val, err in zip(edges[:-1], edges[1:], centers, values, errors):
        writer.writerow(
            {
                "group": spec.group,
                "sample": spec.sample,
                "bin_low": f"{lo:.8g}",
                "bin_high": f"{hi:.8g}",
                "bin_center": f"{center:.8g}",
                "density_pb_per_gev": f"{float(val):.12g}",
                "density_err_pb_per_gev": f"{float(err):.12g}",
                "xsec_pb": f"{spec.xsec_pb:.9g}",
                "stitch_window_low": f"{spec.win_lo:.6g}",
                "stitch_window_high": f"{spec.win_hi:.6g}",
                "used_in_stitch": int((center >= spec.win_lo) and (center < spec.win_hi)),
                "color": spec.color,
                "source_root": str(path),
                "hist_name": spec.hist_name,
            }
        )
    return {
        "group": spec.group,
        "sample": spec.sample,
        "source_root": str(path),
        "hist_name": spec.hist_name,
        "xsec_pb": spec.xsec_pb,
        "stitch_window": [spec.win_lo, spec.win_hi],
        "nevt_read": nevt_read,
        "per_entry_weight_pb_per_GeV": per_entry_weight,
        "bin_width_GeV": bin_width,
    }


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: extract_ppg12_truth_spectrum_root_histograms.py OUT_CSV OUT_JSON", file=sys.stderr)
        return 2
    out_csv = Path(sys.argv[1])
    out_json = Path(sys.argv[2])
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "group",
        "sample",
        "bin_low",
        "bin_high",
        "bin_center",
        "density_pb_per_gev",
        "density_err_pb_per_gev",
        "xsec_pb",
        "stitch_window_low",
        "stitch_window_high",
        "used_in_stitch",
        "color",
        "source_root",
        "hist_name",
    ]
    summaries: list[dict[str, object]] = []
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for spec in PHOTON_SAMPLES + JET_SAMPLES:
            print(f"[extract] {spec.sample} {spec.hist_name}", file=sys.stderr)
            summaries.append(emit_sample(spec, writer))
    out_json.write_text(
        json.dumps(
            {
                "schema": "PPG12_CURRENT_IAN_TRUTH_SPECTRUM_ROOT_HISTOGRAMS_V1",
                "weight_definition": "PPG12 truth_spectrum ROOT histograms; per-entry weight recorded in each ROOT file metadata",
                "photon_windows": {"photon5": [0, 14], "photon10": [14, 22], "photon20": [22, 200]},
                "jet_windows": {"jet8": [9, 14], "jet12": [14, 21], "jet20": [21, 32], "jet30": [32, 42], "jet40": [42, 200]},
                "samples": summaries,
            },
            indent=2,
        )
    )
    print(out_csv)
    print(out_json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
