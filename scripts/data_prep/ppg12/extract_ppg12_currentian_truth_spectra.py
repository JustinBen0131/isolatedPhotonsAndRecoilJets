#!/usr/bin/env python3
"""Extract compact PPG12 current-IAN truth stitching spectra from pp slimtrees."""

from __future__ import annotations
# Keep purpose-folder helpers runnable when invoked directly.
import sys as _codex_sys
from pathlib import Path as _CodexPath
_CODEX_THIS_FILE = _CodexPath(__file__).resolve()
_CODEX_SCRIPTS_DIR = next((p for p in _CODEX_THIS_FILE.parents if p.name == "scripts"), _CODEX_THIS_FILE.parent)
_CODEX_SCRIPTS_DIR_STR = str(_CODEX_SCRIPTS_DIR)
if _CODEX_SCRIPTS_DIR_STR not in _codex_sys.path:
    _codex_sys.path.append(_CODEX_SCRIPTS_DIR_STR)
del _CODEX_THIS_FILE, _CODEX_SCRIPTS_DIR, _CODEX_SCRIPTS_DIR_STR

import csv
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import awkward as ak
import numpy as np
import uproot


SLIMTREE_TEMPLATE = "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/{sample}/condorout/combined.root"


@dataclass(frozen=True)
class Sample:
    group: str
    sample: str
    xsec_pb: float
    win_lo: float
    win_hi: float
    color: str


PHOTON_SAMPLES = [
    Sample("photon", "photon5", 146359.3, 0.0, 14.0, "#E7298A"),
    Sample("photon", "photon10", 6944.675, 14.0, 22.0, "#33A02C"),
    Sample("photon", "photon20", 130.4461, 22.0, 200.0, "#1F78B4"),
]

JET_SAMPLES = [
    Sample("jet", "jet8", 1.3013e7, 9.0, 14.0, "#E7298A"),
    Sample("jet", "jet12", 1.4903e6, 14.0, 21.0, "#33A02C"),
    Sample("jet", "jet20", 6.2623e4, 21.0, 32.0, "#1F78B4"),
    Sample("jet", "jet30", 2.5298e3, 32.0, 42.0, "#FF7F00"),
    Sample("jet", "jet40", 1.3553e2, 42.0, 200.0, "#E7298A"),
]


def leading_photon_pt(tree: uproot.TTree, edges: np.ndarray) -> tuple[np.ndarray, np.ndarray, int, int]:
    hist = np.zeros(len(edges) - 1, dtype=np.float64)
    hist2 = np.zeros_like(hist)
    n_seen = 0
    n_filled = 0
    branches = ["particle_pid", "particle_Pt", "particle_Eta"]
    for arr in tree.iterate(branches, step_size="500 MB"):
        pids = arr["particle_pid"]
        pts = arr["particle_Pt"]
        etas = arr["particle_Eta"]
        mask = (pids == 22) & (abs(etas) < 0.7)
        sel = pts[mask]
        has = ak.num(sel) > 0
        max_pt = ak.to_numpy(ak.fill_none(ak.max(sel, axis=1), -1.0))
        values = max_pt[ak.to_numpy(has)]
        h, _ = np.histogram(values, bins=edges)
        hist += h
        hist2 += h
        n_seen += len(pids)
        n_filled += len(values)
    return hist, hist2, n_seen, n_filled


def leading_jet_pt(tree: uproot.TTree, edges: np.ndarray) -> tuple[np.ndarray, np.ndarray, int, int]:
    hist = np.zeros(len(edges) - 1, dtype=np.float64)
    hist2 = np.zeros_like(hist)
    n_seen = 0
    n_filled = 0
    branch = "jet_truth_Pt_AntiKt_Truth_r04"
    for arr in tree.iterate([branch], step_size="500 MB"):
        pts = arr[branch]
        has = ak.num(pts) > 0
        max_pt = ak.to_numpy(ak.fill_none(ak.max(pts, axis=1), -1.0))
        values = max_pt[ak.to_numpy(has)]
        h, _ = np.histogram(values, bins=edges)
        hist += h
        hist2 += h
        n_seen += len(pts)
        n_filled += len(values)
    return hist, hist2, n_seen, n_filled


def extract_group(samples: list[Sample], edges: np.ndarray, writer: csv.DictWriter) -> list[dict[str, object]]:
    summary: list[dict[str, object]] = []
    for spec in samples:
        path = Path(SLIMTREE_TEMPLATE.format(sample=spec.sample))
        print(f"[extract] {spec.group}:{spec.sample} {path}", file=sys.stderr, flush=True)
        with uproot.open(path) as f:
            tree = f["slimtree"]
            if spec.group == "photon":
                raw, raw_sumw2, n_events, n_filled = leading_photon_pt(tree, edges)
            else:
                raw, raw_sumw2, n_events, n_filled = leading_jet_pt(tree, edges)
        width = np.diff(edges)
        per_event_weight = spec.xsec_pb / float(n_events)
        density = raw * per_event_weight / width
        density_err = np.sqrt(raw_sumw2) * per_event_weight / width
        centers = 0.5 * (edges[:-1] + edges[1:])
        in_window = (centers >= spec.win_lo) & (centers < spec.win_hi)
        for lo, hi, c, n, y, ey, keep in zip(edges[:-1], edges[1:], centers, raw, density, density_err, in_window):
            writer.writerow(
                {
                    "group": spec.group,
                    "sample": spec.sample,
                    "bin_low": f"{lo:.6g}",
                    "bin_high": f"{hi:.6g}",
                    "bin_center": f"{c:.6g}",
                    "raw_events": int(n),
                    "xsec_pb": f"{spec.xsec_pb:.9g}",
                    "n_events": int(n_events),
                    "per_event_weight_pb": f"{per_event_weight:.12g}",
                    "density_pb_per_gev": f"{y:.12g}",
                    "density_err_pb_per_gev": f"{ey:.12g}",
                    "stitch_window_low": f"{spec.win_lo:.6g}",
                    "stitch_window_high": f"{spec.win_hi:.6g}",
                    "used_in_stitch": int(bool(keep)),
                    "color": spec.color,
                }
            )
        summary.append(
            {
                "group": spec.group,
                "sample": spec.sample,
                "path": str(path),
                "xsec_pb": spec.xsec_pb,
                "n_events": n_events,
                "n_with_leading_object": n_filled,
                "per_event_weight_pb": per_event_weight,
                "stitch_window": [spec.win_lo, spec.win_hi],
            }
        )
        print(
            f"[extract] done {spec.sample}: n_events={n_events} n_object={n_filled} weight={per_event_weight:.4e}",
            file=sys.stderr,
            flush=True,
        )
    return summary


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: extract_ppg12_currentian_truth_spectra.py OUT_CSV OUT_JSON", file=sys.stderr)
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
        "raw_events",
        "xsec_pb",
        "n_events",
        "per_event_weight_pb",
        "density_pb_per_gev",
        "density_err_pb_per_gev",
        "stitch_window_low",
        "stitch_window_high",
        "used_in_stitch",
        "color",
    ]
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        photon_summary = extract_group(PHOTON_SAMPLES, np.arange(0.0, 50.0001, 0.5), writer)
        jet_summary = extract_group(JET_SAMPLES, np.arange(8.0, 55.0001, 0.5), writer)

    out_json.write_text(
        json.dumps(
            {
                "schema": "PPG12_CURRENT_IAN_TRUTH_STITCH_SPECTRA_V1",
                "source": SLIMTREE_TEMPLATE,
                "eta_cut_for_leading_photon": "|particle_Eta| < 0.7",
                "photon_windows": {"photon5": [0, 14], "photon10": [14, 22], "photon20": [22, 200]},
                "jet_windows": {"jet8": [9, 14], "jet12": [14, 21], "jet20": [21, 32], "jet30": [32, 42], "jet40": [42, 200]},
                "weight_definition": "density_pb_per_gev = raw_events * sigma_pb / n_events / bin_width",
                "samples": photon_summary + jet_summary,
            },
            indent=2,
        )
    )
    print(out_csv)
    print(out_json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
