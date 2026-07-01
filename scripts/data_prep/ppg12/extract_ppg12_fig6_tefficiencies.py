#!/usr/bin/env python3
"""Print PPG12 Figure 6 TEfficiency values as CSV.

This script is safe to stream into a read-only SDCC Python/ROOT session:

    python3 - < scripts/data_prep/ppg12/extract_ppg12_fig6_tefficiencies.py
"""

import ROOT


SOURCE_ROOT = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_bdt_nom.root"
OBJECTS = [
    ("reco", "eff_reco_eta_0"),
    ("id", "eff_id_eta_0"),
    ("iso", "eff_iso_eta_0"),
    ("all", "eff_all_eta_0"),
]


def main() -> None:
    f = ROOT.TFile.Open(SOURCE_ROOT)
    print("stage,pt_low,pt_high,pt_mid,efficiency,err_low,err_high,source_root,source_object")
    if not f or f.IsZombie():
        raise SystemExit(f"OPEN_FAILED {SOURCE_ROOT}")

    for stage, obj_name in OBJECTS:
        eff = f.Get(obj_name)
        if not eff:
            raise SystemExit(f"MISSING {obj_name}")
        hist = eff.GetTotalHistogram()
        axis = hist.GetXaxis()
        for i in range(1, hist.GetNbinsX() + 1):
            lo = axis.GetBinLowEdge(i)
            hi = axis.GetBinUpEdge(i)
            mid = axis.GetBinCenter(i)
            val = eff.GetEfficiency(i)
            err_low = eff.GetEfficiencyErrorLow(i)
            err_high = eff.GetEfficiencyErrorUp(i)
            print(
                "{},{:.12g},{:.12g},{:.12g},{:.17g},{:.17g},{:.17g},{},{}".format(
                    stage, lo, hi, mid, val, err_low, err_high, SOURCE_ROOT, obj_name
                )
            )


if __name__ == "__main__":
    main()
