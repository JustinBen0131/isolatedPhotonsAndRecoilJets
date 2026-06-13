#!/usr/bin/env python3
"""Extract pp PPG12_TABLE_QA_V1 inclusive sample histograms to JSON.

Run this on a machine where the sample-level RecoilJets ROOT files are
directly visible. The output is a small histogram cache for local plotting; it
does not rewrite or merge ROOT files.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import ROOT


ROOT.gROOT.SetBatch(True)


VARS = (
    "weta_cogx",
    "wphi_cogx",
    "et1",
    "et2",
    "et3",
    "et4",
    "e11_to_e33",
    "e17_to_e77",
    "e32_to_e35",
    "bdt",
    "npb_score",
)
PT_TOKENS = ("0", "2", "3", "1535")
CUTS = ("cut0", "cut1", "cut2", "cut3", "cut4")
DEFAULT_SAMPLES = ("jet5", "jet8", "jet12", "jet20", "jet30", "jet40")
CURRENT_IAN_SAMPLES = ("jet8", "jet12", "jet20", "jet30", "jet40")
CFG_TAG = "preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12"


def hist_payload(hist: ROOT.TH1) -> dict:
    nbins = hist.GetNbinsX()
    x = [float(hist.GetBinCenter(i)) for i in range(1, nbins + 1)]
    y = [float(hist.GetBinContent(i)) for i in range(1, nbins + 1)]
    e = [float(hist.GetBinError(i)) for i in range(1, nbins + 1)]
    sumw = float(sum(y))
    sumw2 = float(sum(err * err for err in e))
    max_frac = float(max(y) / sumw) if sumw > 0 else 0.0
    neff = float((sumw * sumw) / sumw2) if sumw2 > 0 else 0.0
    return {
        "x": x,
        "y": y,
        "e": e,
        "entries": float(hist.GetEntries()),
        "sumw": sumw,
        "sumw2": sumw2,
        "neff": neff,
        "root_effective_entries": float(hist.GetEffectiveEntries()),
        "max_bin_fraction": max_frac,
        "nbins": int(nbins),
        "xlow": float(hist.GetXaxis().GetXmin()),
        "xhigh": float(hist.GetXaxis().GetXmax()),
    }


def combine_payloads(payloads: list[dict]) -> dict | None:
    nonempty = [p for p in payloads if p and p.get("sumw", 0.0) > 0.0]
    if not nonempty:
        return None
    x = nonempty[0]["x"]
    y = [0.0 for _ in x]
    e2 = [0.0 for _ in x]
    entries = 0.0
    sources = []
    for payload in nonempty:
        if len(payload["x"]) != len(x) or any(abs(a - b) > 1e-9 for a, b in zip(payload["x"], x)):
            raise RuntimeError("Cannot combine histograms with different x binning")
        y = [a + b for a, b in zip(y, payload["y"])]
        e2 = [a + b * b for a, b in zip(e2, payload["e"])]
        entries += float(payload.get("entries", 0.0))
        sources.extend(payload.get("source_samples", []))
    e = [v ** 0.5 for v in e2]
    sumw = float(sum(y))
    sumw2 = float(sum(e2))
    max_frac = float(max(y) / sumw) if sumw > 0 else 0.0
    neff = float((sumw * sumw) / sumw2) if sumw2 > 0 else 0.0
    return {
        "x": x,
        "y": y,
        "e": e,
        "entries": entries,
        "sumw": sumw,
        "sumw2": sumw2,
        "neff": neff,
        "max_bin_fraction": max_frac,
        "nbins": len(x),
        "source_samples": sorted(set(sources)),
    }


def open_root(path: Path) -> ROOT.TFile:
    f = ROOT.TFile.Open(str(path))
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open ROOT file: {path}")
    return f


def sample_file(base_dir: Path, sample: str, cfg_tag: str) -> Path:
    return base_dir / f"RecoilJets_{sample}_ALL_{cfg_tag}.root"


def project_x_payload(hist: ROOT.TH2, *, source_histogram: str) -> dict:
    proj = hist.ProjectionX(f"{hist.GetName()}_px_for_json")
    if not proj:
        raise RuntimeError(f"ProjectionX failed for {source_histogram}")
    proj.SetDirectory(0)
    payload = hist_payload(proj)
    payload["source_histogram"] = source_histogram
    payload["source_projection"] = "TH2::ProjectionX, no rebin/no range; renderer applies PPG12 display transform"
    return payload


def read_sample(base_dir: Path, sample: str, cfg_tag: str, source: str) -> dict:
    path = sample_file(base_dir, sample, cfg_tag)
    f = open_root(path)
    hists: dict[str, dict] = {}
    try:
        for var in VARS:
            for pt_token in PT_TOKENS:
                for cut in CUTS:
                    prefix = "h2d" if source == "h2d-projectx" else "h1d"
                    key = f"{prefix}_{var}_eta0_pt{pt_token}_{cut}"
                    hist = f.Get(f"SIM/{key}")
                    if not hist:
                        continue
                    if source == "h2d-projectx":
                        if not hist.InheritsFrom("TH2"):
                            continue
                        payload = project_x_payload(hist, source_histogram=f"SIM/{key}")
                    else:
                        if not hist.InheritsFrom("TH1"):
                            continue
                        payload = hist_payload(hist)
                        payload["source_histogram"] = f"SIM/{key}"
                        payload["source_projection"] = "direct TH1 payload"
                    payload["source_samples"] = [sample]
                    hists[key] = payload
    finally:
        f.Close()
    return {"path": str(path), "hists": hists}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-dir", required=True, type=Path)
    parser.add_argument("--cfg-tag", default=CFG_TAG)
    parser.add_argument("--samples", nargs="+", default=list(DEFAULT_SAMPLES))
    parser.add_argument(
        "--source",
        choices=("h2d-projectx", "h1d"),
        default="h2d-projectx",
        help="Use PPG12-style h2d ProjectionX payloads by default; h1d is retained for old diagnostics.",
    )
    args = parser.parse_args()

    samples: dict[str, dict] = {}
    for sample in args.samples:
        samples[sample] = read_sample(args.base_dir, sample, args.cfg_tag, args.source)

    sample_sets = {
        "current_ian_jet8to40": [s for s in CURRENT_IAN_SAMPLES if s in samples],
        "available_jet5to40": list(samples),
    }
    combined: dict[str, dict] = {}
    all_keys = sorted({key for sample in samples.values() for key in sample["hists"]})
    for set_name, set_samples in sample_sets.items():
        combined[set_name] = {}
        for key in all_keys:
            payloads = [samples[sample]["hists"].get(key) for sample in set_samples]
            combo = combine_payloads([p for p in payloads if p is not None])
            if combo is not None:
                combined[set_name][key] = combo

    output = {
        "schema": (
            "THE42_PPG12_TABLE_QA_INCLUSIVE_SAMPLE_HIST_CACHE_V2"
            if args.source == "h2d-projectx"
            else "THE42_PPG12_TABLE_QA_INCLUSIVE_SAMPLE_HIST_CACHE_V1"
        ),
        "note": (
            "Sample-level PPG12_TABLE_QA histograms are already weighted in "
            "RecoilJets by ppg12InclusiveJetSliceXsecPb(slice)/7.3113. "
            "These combined caches intentionally apply no finalStitch weight."
        ),
        "source": args.source,
        "topdir": "SIM",
        "cfg_tag": args.cfg_tag,
        "base_dir": str(args.base_dir),
        "samples": samples,
        "sample_sets": sample_sets,
        "combined": combined,
    }
    print(json.dumps(output, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
