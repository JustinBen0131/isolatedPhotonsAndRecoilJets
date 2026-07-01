#!/usr/bin/env python3
"""Replay compact AuAu photon-candidate skim rows into QA histograms.

This is the fast offline side of THE-83/THE-69C.  The skim is produced once by
Fun4All, then this script regenerates the photon-ID rows needed for QA without
another DST pass.
"""

from __future__ import annotations

import argparse
import json
import math
from array import array
from pathlib import Path
from typing import Any

import ROOT


DEFAULT_PT_BINS = [15, 16, 18, 20, 22, 24, 26, 28, 30, 32, 35]
DEFAULT_CENT_BINS = [(0, 20), (20, 50), (50, 80)]
DEFAULT_ROWS = (
    "preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement",
    "preselectionNewPPG12_tightReference_nonTightReference",
)
TREE_NAME = "AuAuPhotonCandidateSkim"
WP80_INTERCEPT = 0.53471108
WP80_SLOPE = 0.0012284143


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", help="Input ROOT files containing AuAuPhotonCandidateSkim")
    parser.add_argument("--output", required=True, help="Output replay ROOT file")
    parser.add_argument("--manifest", help="Optional replay manifest JSON")
    parser.add_argument("--campaign-tag", default="", help="Optional campaign tag recorded in the manifest")
    parser.add_argument("--tree", default=TREE_NAME)
    parser.add_argument("--pt-bins", default=",".join(str(x) for x in DEFAULT_PT_BINS))
    parser.add_argument(
        "--rows",
        default=",".join(DEFAULT_ROWS),
        help="Comma-separated replay rows. Supported: default BDT and reference box rows.",
    )
    parser.add_argument("--max-entries", type=int, default=0, help="Debug cap; 0 means all entries")
    return parser.parse_args()


def parse_float_list(text: str) -> list[float]:
    return [float(x) for x in text.split(",") if x.strip()]


def cent_label(lo: int, hi: int) -> str:
    return f"cent{lo}_{hi}"


def row_key(row: str) -> str:
    if "tightAuAuCentInputBase3x3BDT" in row:
        return "baseline_bdt"
    if "tightReference" in row:
        return "box_reference"
    raise ValueError(f"Unsupported replay row: {row}")


def branch_for_row(row: str) -> tuple[str, str, str]:
    key = row_key(row)
    if key == "baseline_bdt":
        return "baseline_bdt_abcd_region", "baseline_bdt_tight", "baseline_bdt_nontight"
    return "box_abcd_region", "box_tight", "box_nontight"


def make_h1(name: str, title: str, bins: list[float]) -> ROOT.TH1F:
    arr = array("d", [float(x) for x in bins])
    h = ROOT.TH1F(name, title, len(bins) - 1, arr)
    h.Sumw2()
    return h


def safe_get(entry: Any, name: str, default: float = math.nan) -> float:
    try:
        return getattr(entry, name)
    except AttributeError:
        return default


def fill_ratio(num: ROOT.TH1F, den: ROOT.TH1F, out: ROOT.TH1F) -> None:
    for ib in range(1, out.GetNbinsX() + 1):
        n = num.GetBinContent(ib)
        d = den.GetBinContent(ib)
        if d <= 0:
            out.SetBinContent(ib, 0.0)
            out.SetBinError(ib, 0.0)
            continue
        r = n / d
        # Binomial-like uncertainty is sufficient for a replay QA product.
        err = math.sqrt(max(r * (1.0 - r), 0.0) / d) if 0.0 <= r <= 1.0 else 0.0
        out.SetBinContent(ib, r)
        out.SetBinError(ib, err)


def main() -> int:
    args = parse_args()
    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.SetDefaultSumw2(True)

    pt_bins = parse_float_list(args.pt_bins)
    rows = [r.strip() for r in args.rows.split(",") if r.strip()]

    chain = ROOT.TChain(args.tree)
    for path in args.inputs:
        chain.Add(path)
    entries = int(chain.GetEntries())
    if entries <= 0:
        raise RuntimeError(f"No entries found in tree {args.tree} from inputs: {args.inputs}")

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fout = ROOT.TFile(str(output), "RECREATE")
    if not fout or fout.IsZombie():
        raise RuntimeError(f"Could not open output ROOT: {output}")

    histograms: dict[str, ROOT.TH1] = {}
    purity_inputs: dict[tuple[str, str, str], ROOT.TH1F] = {}

    for row in rows:
        key = row_key(row)
        row_dir = fout.mkdir(key)
        row_dir.cd()
        for lo, hi in DEFAULT_CENT_BINS:
            clabel = cent_label(lo, hi)
            histograms[f"{key}/h_bdt_score_{clabel}"] = ROOT.TH1F(
                f"h_bdt_score_{clabel}", ";AuAu BDT score;Candidates", 100, 0.0, 1.0
            )
            histograms[f"{key}/h_e11_over_e33_all_{clabel}"] = ROOT.TH1F(
                f"h_e11_over_e33_all_{clabel}", ";E_{1x1}/E_{3x3};Candidates", 100, 0.0, 1.2
            )
            histograms[f"{key}/h_e11_over_e33_tight_iso_{clabel}"] = ROOT.TH1F(
                f"h_e11_over_e33_tight_iso_{clabel}", ";E_{1x1}/E_{3x3};Candidates", 100, 0.0, 1.2
            )
            histograms[f"{key}/h_eiso_all_{clabel}"] = ROOT.TH1F(
                f"h_eiso_all_{clabel}", ";E_{T}^{iso} [GeV];Candidates", 120, -20.0, 40.0
            )
            histograms[f"{key}/h_eiso_tight_{clabel}"] = ROOT.TH1F(
                f"h_eiso_tight_{clabel}", ";E_{T}^{iso} [GeV];Candidates", 120, -20.0, 40.0
            )
            histograms[f"{key}/h_tight_iso_yield_{clabel}"] = make_h1(
                f"h_tight_iso_yield_{clabel}", ";cluster E_{T} [GeV];A-region candidates", pt_bins
            )
            histograms[f"{key}/h_lead_xj_tight_iso_{clabel}"] = ROOT.TH1F(
                f"h_lead_xj_tight_iso_{clabel}", ";x_{J#gamma}^{lead};Candidates", 60, 0.0, 3.0
            )
            for region in ("A", "B", "C", "D"):
                h = make_h1(
                    f"h_abcd_{region}_{clabel}",
                    f";cluster E_{{T}} [GeV];region {region} candidates",
                    pt_bins,
                )
                histograms[f"{key}/{h.GetName()}"] = h
                purity_inputs[(key, clabel, region)] = h
            for h in list(histograms.values()):
                h.Sumw2()
        fout.cd()

    max_entries = args.max_entries if args.max_entries > 0 else entries
    for i, entry in enumerate(chain):
        if i >= max_entries:
            break
        pt = float(safe_get(entry, "cluster_Et"))
        cent = float(safe_get(entry, "centrality"))
        if not (math.isfinite(pt) and math.isfinite(cent)):
            continue
        if pt < pt_bins[0] or pt >= pt_bins[-1]:
            continue

        clabel = None
        for lo, hi in DEFAULT_CENT_BINS:
            if cent >= lo and cent < hi:
                clabel = cent_label(lo, hi)
                break
        if clabel is None:
            continue

        weight = float(safe_get(entry, "event_weight", 1.0))
        if not math.isfinite(weight):
            weight = 1.0
        bdt_score = float(safe_get(entry, "auau_tight_bdt_score", -2.0))
        e11e33 = float(safe_get(entry, "e11_over_e33", math.nan))
        eiso = float(safe_get(entry, "reco_eiso", math.nan))
        xj = float(safe_get(entry, "lead_xj", math.nan))

        for row in rows:
            key = row_key(row)
            region_branch, tight_branch, _ = branch_for_row(row)
            region_code = int(safe_get(entry, region_branch, 0))
            tight_pass = int(safe_get(entry, tight_branch, 0)) == 1

            if 0.0 <= bdt_score <= 1.0:
                histograms[f"{key}/h_bdt_score_{clabel}"].Fill(bdt_score, weight)
            if math.isfinite(e11e33):
                histograms[f"{key}/h_e11_over_e33_all_{clabel}"].Fill(e11e33, weight)
            if math.isfinite(eiso):
                histograms[f"{key}/h_eiso_all_{clabel}"].Fill(eiso, weight)
                if tight_pass:
                    histograms[f"{key}/h_eiso_tight_{clabel}"].Fill(eiso, weight)

            if region_code in (1, 2, 3, 4):
                region = "ABCD"[region_code - 1]
                purity_inputs[(key, clabel, region)].Fill(pt, weight)
                if region == "A":
                    histograms[f"{key}/h_tight_iso_yield_{clabel}"].Fill(pt, weight)
                    if math.isfinite(e11e33):
                        histograms[f"{key}/h_e11_over_e33_tight_iso_{clabel}"].Fill(e11e33, weight)
                    if math.isfinite(xj) and xj >= 0.0:
                        histograms[f"{key}/h_lead_xj_tight_iso_{clabel}"].Fill(xj, weight)

    for row in rows:
        key = row_key(row)
        row_dir = fout.GetDirectory(key)
        row_dir.cd()
        for lo, hi in DEFAULT_CENT_BINS:
            clabel = cent_label(lo, hi)
            h_a = purity_inputs[(key, clabel, "A")]
            h_b = purity_inputs[(key, clabel, "B")]
            h_c = purity_inputs[(key, clabel, "C")]
            h_d = purity_inputs[(key, clabel, "D")]
            h_fake = make_h1(f"h_abcd_fake_estimate_BC_over_D_{clabel}", ";cluster E_{T} [GeV];BC/D", pt_bins)
            h_raw_purity = make_h1(f"h_raw_abcd_purity_{clabel}", ";cluster E_{T} [GeV];1 - BC/(AD)", pt_bins)
            for ib in range(1, h_a.GetNbinsX() + 1):
                a = h_a.GetBinContent(ib)
                b = h_b.GetBinContent(ib)
                c = h_c.GetBinContent(ib)
                d = h_d.GetBinContent(ib)
                fake = (b * c / d) if d > 0.0 else 0.0
                purity = (a - fake) / a if a > 0.0 else 0.0
                h_fake.SetBinContent(ib, fake)
                h_raw_purity.SetBinContent(ib, purity)
            histograms[f"{key}/{h_fake.GetName()}"] = h_fake
            histograms[f"{key}/{h_raw_purity.GetName()}"] = h_raw_purity
        fout.cd()

    for path, hist in histograms.items():
        dirname = str(Path(path).parent)
        fout.cd(dirname if dirname != "." else "")
        hist.Write("", ROOT.TObject.kOverwrite)

    fout.Close()

    manifest = {
        "inputs": args.inputs,
        "output": str(output),
        "campaign_tag": args.campaign_tag,
        "tree": args.tree,
        "entries_seen": entries,
        "entries_processed": min(max_entries, entries),
        "rows": rows,
        "row_keys": {row: row_key(row) for row in rows},
        "pt_bins": pt_bins,
        "centrality_bins": DEFAULT_CENT_BINS,
        "skim_contract": {
            "tree": TREE_NAME,
            "baseline_bdt_fields": [
                "auau_tight_bdt_score",
                "baseline_wp80_threshold",
                "baseline_bdt_tight",
                "baseline_bdt_nontight",
                "baseline_bdt_abcd_region",
            ],
            "note": "baseline fields use the fixed THE-57/THE-69 default WP80 contract; replay does not retune the cut",
        },
        "default_wp80_contract": {
            "model_id": "centAsFeatBase3x3_pt15to35",
            "mode": "centlinear",
            "formula": "T80(c)=0.53471108+0.0012284143*c",
            "intercept": WP80_INTERCEPT,
            "slope": WP80_SLOPE,
            "pt_window": [15.0, 35.0],
        },
    }
    if args.manifest:
        manifest_path = Path(args.manifest)
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
