#!/usr/bin/env python3
"""Extract fine 1 GeV reco-cluster ET component spectra from ROOT outputs.

This is intentionally a narrow diagnostic helper for the Blair leakage check.
It sums ABCD regions and all centrality slices for the opt-in histogram family
created with RJ_RECO_CLUSTER_ET_FINE_DIAG=1.
"""

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

import argparse
import csv
import json
import math
from pathlib import Path

import ROOT


REPO = Path(__file__).resolve().parents[1]
OUTDIR = REPO / "dataOutput/stitchDiagnostics/focus21_clusterEt_leakage_jet1234_20260526"
DEFAULT_CSV = OUTDIR / "inclusive_jet1234_reco_cluster_et_weighted_components_fine1gev12to50.csv"
DEFAULT_SUMMARY = OUTDIR / "inclusive_jet1234_reco_cluster_et_weighted_components_fine1gev12to50_summary.json"

SAMPLES = ("Jet12", "Jet20", "Jet30", "Jet40")
REGIONS = ("A", "B", "C", "D")
DEFAULT_GATES = {
    "Jet12": "12 <= pT_truth_jet < 21",
    "Jet20": "21 <= pT_truth_jet < 31",
    "Jet30": "31 <= pT_truth_jet < 41",
    "Jet40": "pT_truth_jet >= 41",
}
DEFAULT_SIGMAS = {
    "Jet12": 1.22772477e6,
    "Jet20": 3.88117850e4,
    "Jet30": 1.73665908e3,
    "Jet40": 1.00642312e2,
}
DEFAULT_SCALES = {
    "Jet12": 12199.8889732,
    "Jet20": 385.6418704,
    "Jet30": 17.2567536,
    "Jet40": 1.0,
}
DEFAULT_ROOT_NAMES = {
    "Jet12": "RecoilJets_embeddedJet12_ALL_preselectionReference_tightReference_nonTightReference_baseVariant.root",
    "Jet20": "RecoilJets_embeddedJet20_ALL_preselectionReference_tightReference_nonTightReference_baseVariant.root",
    "Jet30": "RecoilJets_embeddedJet30_ALL_preselectionReference_tightReference_nonTightReference_baseVariant.root",
    "Jet40": "RecoilJets_embeddedJet40_ALL_preselectionReference_tightReference_nonTightReference_baseVariant.root",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-base", type=Path, help="Remote/local output_<tag> base containing simembeddedinclusive sample ALL ROOTs.")
    parser.add_argument("--chunk-base", type=Path, help="Directory containing chunkMerge_embeddedJet*_grp*.root files.")
    parser.add_argument("--jet12-root", type=Path)
    parser.add_argument("--jet20-root", type=Path)
    parser.add_argument("--jet30-root", type=Path)
    parser.add_argument("--jet40-root", type=Path)
    parser.add_argument("--jet12-scale", type=float, default=DEFAULT_SCALES["Jet12"])
    parser.add_argument("--jet20-scale", type=float, default=DEFAULT_SCALES["Jet20"])
    parser.add_argument("--jet30-scale", type=float, default=DEFAULT_SCALES["Jet30"])
    parser.add_argument("--jet40-scale", type=float, default=DEFAULT_SCALES["Jet40"])
    parser.add_argument("--jet12-sigma", type=float, default=DEFAULT_SIGMAS["Jet12"])
    parser.add_argument("--jet20-sigma", type=float, default=DEFAULT_SIGMAS["Jet20"])
    parser.add_argument("--jet30-sigma", type=float, default=DEFAULT_SIGMAS["Jet30"])
    parser.add_argument("--jet40-sigma", type=float, default=DEFAULT_SIGMAS["Jet40"])
    parser.add_argument("--iso-tag", default="isoR30_fixedIso4GeV")
    parser.add_argument("--hist-token", default="fine1GeV_12to50")
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    return parser.parse_args()


def sample_roots(args: argparse.Namespace) -> dict[str, list[Path]]:
    roots = {
        "Jet12": args.jet12_root,
        "Jet20": args.jet20_root,
        "Jet30": args.jet30_root,
        "Jet40": args.jet40_root,
    }
    if args.output_base is not None:
        for sample, name in DEFAULT_ROOT_NAMES.items():
            roots[sample] = args.output_base / "simembeddedinclusive" / name
    if args.chunk_base is not None:
        chunk_names = {
            "Jet12": "chunkMerge_embeddedJet12_grp*.root",
            "Jet20": "chunkMerge_embeddedJet20_grp*.root",
            "Jet30": "chunkMerge_embeddedJet30_grp*.root",
            "Jet40": "chunkMerge_embeddedJet40_grp*.root",
        }
        out: dict[str, list[Path]] = {}
        for sample, pattern in chunk_names.items():
            paths = sorted(args.chunk_base.glob(pattern))
            if not paths:
                raise SystemExit(f"no chunk ROOTs found for {sample}: {args.chunk_base / pattern}")
            out[sample] = paths
        return out
    missing = [sample for sample, path in roots.items() if path is None]
    if missing:
        raise SystemExit(f"missing ROOT path(s) for: {', '.join(missing)}")
    return {sample: [Path(path)] for sample, path in roots.items()}


def matching_hists(root_file: ROOT.TFile, region: str, hist_token: str, iso_tag: str) -> list[ROOT.TH1]:
    sim = root_file.Get("SIM")
    if not sim:
        raise KeyError("missing SIM directory")
    prefix = f"h_recoClusterEt_ABCD_{region}_{hist_token}"
    if iso_tag:
        exact = f"{prefix}_{iso_tag}"
        obj = sim.Get(exact)
        if obj and obj.InheritsFrom("TH1"):
            hist = obj.Clone(f"{exact}_clone")
            hist.SetDirectory(0)
            return [hist]
    matches: list[tuple[str, ROOT.TH1]] = []
    for key in sim.GetListOfKeys():
        name = key.GetName()
        if not name.startswith(prefix):
            continue
        if iso_tag and iso_tag not in name:
            continue
        obj = sim.Get(name)
        if obj and obj.InheritsFrom("TH1"):
            hist = obj.Clone(f"{name}_clone")
            hist.SetDirectory(0)
            matches.append((name, hist))
    if not matches:
        raise KeyError(f"no SIM/{prefix}* histograms found for iso tag {iso_tag}")
    inclusive = [hist for name, hist in matches if "_cent_" not in name]
    if inclusive:
        return inclusive
    return [hist for _, hist in matches]


def abcd_cent_sum(paths: list[Path], hist_token: str, iso_tag: str) -> ROOT.TH1:
    total = None
    for path in paths:
        root_file = ROOT.TFile.Open(str(path), "READ")
        if not root_file or root_file.IsZombie():
            raise OSError(f"could not open ROOT file: {path}")
        for region in REGIONS:
            for hist in matching_hists(root_file, region, hist_token, iso_tag):
                if total is None:
                    total = hist.Clone("abcd_cent_sum")
                    total.SetDirectory(0)
                else:
                    total.Add(hist)
        root_file.Close()
    if total is None:
        raise RuntimeError(f"no histograms accumulated from {paths}")
    return total


def row_records(sample: str, hist: ROOT.TH1, scale: float, sigma: float, paths: list[Path]) -> list[dict[str, object]]:
    axis = hist.GetXaxis()
    rows: list[dict[str, object]] = []
    for ibin in range(1, hist.GetNbinsX() + 1):
        lo = float(axis.GetBinLowEdge(ibin))
        hi = float(axis.GetBinUpEdge(ibin))
        raw = float(hist.GetBinContent(ibin))
        raw_err = float(hist.GetBinError(ibin))
        rows.append(
            {
                "observable": "ABCD-summed reco photon-cluster ET, fine 1 GeV bins",
                "sample": sample,
                "truth_jet_gate": DEFAULT_GATES[sample],
                "sigma_eff_pb": sigma,
                "merge_scale": scale,
                "hist_family": "sum_cent h_recoClusterEt_ABCD_{A,B,C,D}_{hist token}_{iso view}_cent_*",
                "region_sum": "A+B+C+D",
                "centrality_sum": "0-80%",
                "bin_index": ibin,
                "x_low": lo,
                "x_high": hi,
                "x_center": 0.5 * (lo + hi),
                "x_err_low": 0.5 * (hi - lo),
                "x_err_high": 0.5 * (hi - lo),
                "is_overflow": 0,
                "raw_entries": raw,
                "raw_error": raw_err,
                "weighted_entries": raw * scale,
                "weighted_error": raw_err * scale,
                "root_path": ";".join(str(path) for path in paths),
            }
        )
    return rows


def summarize_range(rows: list[dict[str, object]], lo: float, hi: float) -> dict[str, object]:
    by_sample = {sample: 0.0 for sample in SAMPLES}
    for row in rows:
        xlo = float(row["x_low"])
        xhi = float(row["x_high"])
        if xlo < lo or xhi > hi:
            continue
        by_sample[str(row["sample"])] += float(row["weighted_entries"])
    total = sum(by_sample.values())
    fractions = {
        sample: (by_sample[sample] / total if total > 0 else 0.0)
        for sample in SAMPLES
    }
    return {
        "x_range": f"{lo:g}-{hi:g}",
        "weighted_sum": total,
        "fractions": fractions,
        "weighted_entries": by_sample,
    }


def write_outputs(args: argparse.Namespace) -> None:
    roots = sample_roots(args)
    scales = {"Jet12": args.jet12_scale, "Jet20": args.jet20_scale, "Jet30": args.jet30_scale, "Jet40": args.jet40_scale}
    sigmas = {"Jet12": args.jet12_sigma, "Jet20": args.jet20_sigma, "Jet30": args.jet30_sigma, "Jet40": args.jet40_sigma}

    all_rows: list[dict[str, object]] = []
    for sample in SAMPLES:
        hist = abcd_cent_sum(roots[sample], args.hist_token, args.iso_tag)
        all_rows.extend(row_records(sample, hist, scales[sample], sigmas[sample], roots[sample]))

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(all_rows[0].keys()))
        writer.writeheader()
        writer.writerows(all_rows)

    summary = {
        "input_roots": {sample: [str(path) for path in paths] for sample, paths in roots.items()},
        "output_csv": str(args.csv),
        "observable": "ABCD-summed reco photon-cluster ET, fine 1 GeV bins",
        "centrality_sum": "0-80%",
        "iso_tag": args.iso_tag,
        "truth_jet_gates": DEFAULT_GATES,
        "sample_metadata": [
            {
                "sample": sample,
                "truth_jet_gate": DEFAULT_GATES[sample],
                "sigma_eff_pb": sigmas[sample],
                "merge_scale": scales[sample],
                "root_path": ";".join(str(path) for path in roots[sample]),
            }
            for sample in SAMPLES
        ],
        "bins_of_interest": [
            summarize_range(all_rows, 22.0, 24.0),
            summarize_range(all_rows, 24.0, 26.0),
            summarize_range(all_rows, 26.0, 35.0),
            summarize_range(all_rows, 35.0, 40.0),
            summarize_range(all_rows, 40.0, 50.0),
        ],
        "plot_caveat": "Dedicated diagnostic histogram with true 1 GeV bins from 12 to 50 GeV.",
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Wrote {args.csv}")
    print(f"Wrote {args.summary}")


def main() -> None:
    ROOT.gROOT.SetBatch(True)
    write_outputs(parse_args())


if __name__ == "__main__":
    main()
