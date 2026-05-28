#!/usr/bin/env python3
"""Extract summed scaled-trigger histogram bins for selected clean runs.

Run this on SDCC where PyROOT can read the per-run ROOT files. The run list is
provided through SCALED_TRIGGER_RUNS as a comma-separated list. Output is a
small CSV to stdout so the plotting step can stay local.
"""

from __future__ import annotations

import csv
import os
import sys

import ROOT  # type: ignore


DEFAULT_INPUT_DIR = (
    "/sphenix/u/patsfan753/scratch/thesisAnalysis/output/auau/perRun/"
    "jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionReference_"
    "tightReference_nonTightReference_scaledTriggerStudy"
)

HISTS = {
    "mbd": (
        "MBD_NS_geq_2_vtx_lt_150/"
        "h_maxEnergyClus_NewTriggerFilling_perRunCorrected_MBD_NS_geq_2_vtx_lt_150"
    ),
    "p10": "Photon_10/h_maxEnergyClus_NewTriggerFilling_perRunCorrected_Photon_10",
    "p12": "Photon_12/h_maxEnergyClus_NewTriggerFilling_perRunCorrected_Photon_12",
}

DEFAULT_CLEAN_RUNS = [
    68502,
    71224,
    71252,
    71256,
    71278,
    71329,
    71346,
    71378,
    71406,
    71419,
    71455,
    71460,
    71467,
    71801,
    71812,
    71844,
    71943,
    71954,
    71957,
    71985,
    72019,
    72086,
    72138,
    72149,
    72151,
    72170,
    72171,
    72200,
    72220,
    72272,
    72286,
    72369,
    72370,
    72372,
    72387,
    72434,
    72478,
    72496,
    72505,
    72508,
    72526,
    72569,
    72573,
    72593,
    72621,
    72627,
    72641,
    72664,
    72665,
    72862,
    72864,
    72867,
    72968,
    72969,
    73004,
    73038,
    73056,
    73057,
    74399,
    74413,
    74518,
    74523,
    74550,
]


def parse_runs() -> list[int]:
    raw = os.environ.get("SCALED_TRIGGER_RUNS", "").strip()
    if not raw:
        return DEFAULT_CLEAN_RUNS
    return [int(token) for token in raw.replace(" ", "").split(",") if token]


def main() -> int:
    ROOT.gROOT.SetBatch(True)
    runs = parse_runs()
    input_dir = os.environ.get("SCALED_TRIGGER_INPUT_DIR", DEFAULT_INPUT_DIR)

    summed: dict[str, object] = {}
    missing: list[str] = []
    used: list[int] = []

    for run in runs:
        path = os.path.join(input_dir, f"chunkMerge_run_{run:08d}.root")
        root_file = ROOT.TFile.Open(path)
        if not root_file or root_file.IsZombie():
            missing.append(f"{run}: could not open {path}")
            continue
        hists = {name: root_file.Get(hist_path) for name, hist_path in HISTS.items()}
        if any(not hist for hist in hists.values()):
            missing.append(f"{run}: missing one or more expected histograms")
            root_file.Close()
            continue

        if not summed:
            for name, hist in hists.items():
                clone = hist.Clone(f"{name}_sum_clean_runs")
                clone.SetDirectory(0)
                clone.Reset()
                summed[name] = clone

        for name, hist in hists.items():
            summed[name].Add(hist)  # type: ignore[union-attr]
        used.append(run)
        root_file.Close()

    if missing:
        for line in missing:
            print(f"WARNING,{line}", file=sys.stderr)
    if not summed:
        raise RuntimeError("No histograms were summed.")

    writer = csv.writer(sys.stdout)
    writer.writerow(["BIN", "bin", "center", "width", "mbd", "p10", "p12"])
    mbd = summed["mbd"]
    p10 = summed["p10"]
    p12 = summed["p12"]
    xaxis = mbd.GetXaxis()  # type: ignore[union-attr]
    for ibin in range(1, mbd.GetNbinsX() + 1):  # type: ignore[union-attr]
        writer.writerow(
            [
                "BIN",
                ibin,
                xaxis.GetBinCenter(ibin),
                xaxis.GetBinWidth(ibin),
                mbd.GetBinContent(ibin),  # type: ignore[union-attr]
                p10.GetBinContent(ibin),  # type: ignore[union-attr]
                p12.GetBinContent(ibin),  # type: ignore[union-attr]
            ]
        )
    writer.writerow(["SUMMARY", "requested_runs", len(runs)])
    writer.writerow(["SUMMARY", "used_runs", len(used)])
    writer.writerow(["SUMMARY", "first_run", min(used)])
    writer.writerow(["SUMMARY", "last_run", max(used)])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
