#!/usr/bin/env python3
"""Audit THE-76 Fig.24/Fig.3 high-pT Eiso tail discrepancies.

This is a read-only diagnostic: it inspects existing ROOT/CSV/JSON artifacts
and writes compact audit tables.  It does not submit, merge, transfer, or alter
analysis code.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
import textwrap
from pathlib import Path
from typing import Any

import numpy as np
import uproot


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
CAMPAIGN = "the76_ppg12_fig24_photonjet_fix_20260702_014217"
DEFAULT_CURRENT_ROOT = (
    REPO
    / "dataOutput/ppg12Parity"
    / CAMPAIGN
    / "final_roots/photonjet/RecoilJets_photonjet5plus10plus20_MERGED.root"
)
DEFAULT_PPG12_FIG24_REMOTE_ROOT = (
    "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_bdt_nom.root"
)
DEFAULT_PPG12_FIG24_OBJECT = "h_singal_reco_isoET_0"
DEFAULT_FIG24_CUTOFF_CSV = (
    REPO
    / "dataOutput/ppg12Parity"
    / CAMPAIGN
    / "fig24_iso_efficiency/fig24_ppg12_cutoff_points_from_sdcc_nom_root.csv"
)
DEFAULT_FIG24_VARIANT_CSV = (
    REPO
    / "dataOutput/ppg12Parity"
    / CAMPAIGN
    / "fig24_iso_efficiency/fig24_weight_variant_localization.csv"
)
DEFAULT_FIG3_RATIO_CSV = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_parity_full_20260701_003024"
    / "data_iso_template_fig3/ppg12_fig3_iso_template_pt16_22_sdcc_vs_current_PPG12_scaledtrigger30_ratio.csv"
)
DEFAULT_FIG3_JSON = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_parity_full_20260701_003024"
    / "data_iso_template_fig3/ppg12_sdcc_fig3_iso_template_pt16_22_tight_data.json"
)
DEFAULT_OUTDIR = (
    REPO
    / "dataOutput/ppg12Parity/control_plane/audits/the76_fig24_fig3_eiso_tail_20260706"
)

PT_SLICES = [
    (16.0, 22.0),
    (22.0, 26.0),
    (26.0, 30.0),
    (30.0, 34.0),
    (34.0, 36.0),
]
FIXED_EISO_VALUES = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
EFFS = [0.70, 0.80, 0.90]


def root_find_bin(edges: np.ndarray, value: float) -> int:
    nbins = len(edges) - 1
    if value < edges[0]:
        return 1
    if value >= edges[-1]:
        return nbins
    return int(np.searchsorted(edges, value, side="right"))


def weighted_quantile_from_bins(edges: np.ndarray, values: np.ndarray, frac: float, *, low: float | None = None) -> float:
    if values.size == 0:
        return math.nan
    total = float(np.nansum(values))
    if total <= 0:
        return math.nan
    start = 1
    if low is not None:
        start = root_find_bin(edges, low)
    current = 0.0
    for ibin in range(start, len(edges)):
        current += float(values[ibin - 1])
        if current >= total * frac:
            return float(edges[ibin])
    return float(edges[-1])


def distribution_stats(edges: np.ndarray, values: np.ndarray, variances: np.ndarray | None = None) -> dict[str, float]:
    values = np.asarray(values, dtype=float)
    centers = 0.5 * (edges[:-1] + edges[1:])
    total = float(np.nansum(values))
    if variances is None:
        variances = np.clip(values, 0.0, None)
    var_sum = float(np.nansum(variances))
    if total <= 0:
        out = {
            "entries": total,
            "effective_entries": 0.0,
            "mean": math.nan,
            "median": math.nan,
            "rms": math.nan,
            "q70": math.nan,
            "q80": math.nan,
            "q90": math.nan,
            "cut70_isolow_minus1": math.nan,
            "cut80_isolow_minus1": math.nan,
            "cut90_isolow_minus1": math.nan,
        }
    else:
        mean = float(np.nansum(values * centers) / total)
        rms = float(math.sqrt(max(0.0, np.nansum(values * (centers - mean) ** 2) / total)))
        out = {
            "entries": total,
            "effective_entries": float(total * total / var_sum) if var_sum > 0 else math.nan,
            "mean": mean,
            "median": weighted_quantile_from_bins(edges, values, 0.50),
            "rms": rms,
            "q70": weighted_quantile_from_bins(edges, values, 0.70),
            "q80": weighted_quantile_from_bins(edges, values, 0.80),
            "q90": weighted_quantile_from_bins(edges, values, 0.90),
            "cut70_isolow_minus1": weighted_quantile_from_bins(edges, values, 0.70, low=-1.0),
            "cut80_isolow_minus1": weighted_quantile_from_bins(edges, values, 0.80, low=-1.0),
            "cut90_isolow_minus1": weighted_quantile_from_bins(edges, values, 0.90, low=-1.0),
        }
    for eiso in FIXED_EISO_VALUES:
        mask = edges[1:] <= eiso
        out[f"cdf_le_{eiso:g}"] = float(np.nansum(values[mask]) / total) if total > 0 else math.nan
    return out


def th2_slice(values: np.ndarray, variances: np.ndarray | None, xedges: np.ndarray, yedges: np.ndarray, lo: float, hi: float) -> tuple[np.ndarray, np.ndarray | None]:
    xlo = root_find_bin(xedges, lo)
    xhi = root_find_bin(xedges, hi)
    xlo = max(1, min(values.shape[0], xlo))
    xhi = max(1, min(values.shape[0], xhi))
    if xhi < xlo:
        xlo, xhi = xhi, xlo
    vals = np.nansum(values[xlo - 1 : xhi, :], axis=0)
    vars_out = None
    if variances is not None:
        vars_out = np.nansum(variances[xlo - 1 : xhi, :], axis=0)
    return vals, vars_out


def load_th2_from_root(path: Path, key: str) -> dict[str, Any]:
    with uproot.open(path) as handle:
        if key not in handle:
            raise KeyError(f"{key} missing in {path}")
        hist = handle[key]
        values, xedges, yedges = hist.to_numpy(flow=False)
        variances = hist.variances(flow=False)
    return {
        "source": str(path),
        "object": key,
        "values": values.astype(float),
        "variances": None if variances is None else variances.astype(float),
        "xedges": xedges.astype(float),
        "yedges": yedges.astype(float),
    }


def remote_extract_ppg12_fig24(remote_root: str, object_name: str, login_host: str, worker_host: str) -> dict[str, Any]:
    remote_py = textwrap.dedent(
        f"""
        import json
        import ROOT
        ROOT.gROOT.SetBatch(True)
        path = {remote_root!r}
        obj = {object_name!r}
        f = ROOT.TFile.Open(path)
        if not f or f.IsZombie():
            raise SystemExit("could not open " + path)
        h = f.Get(obj)
        if not h:
            raise SystemExit("missing " + obj)
        payload = {{
            "source": path,
            "object": obj,
            "xedges": [],
            "yedges": [],
            "values": [],
            "variances": [],
        }}
        for ix in range(1, h.GetNbinsX() + 1):
            if ix == 1:
                payload["xedges"].append(float(h.GetXaxis().GetBinLowEdge(ix)))
            payload["xedges"].append(float(h.GetXaxis().GetBinUpEdge(ix)))
        for iy in range(1, h.GetNbinsY() + 1):
            if iy == 1:
                payload["yedges"].append(float(h.GetYaxis().GetBinLowEdge(iy)))
            payload["yedges"].append(float(h.GetYaxis().GetBinUpEdge(iy)))
        for ix in range(1, h.GetNbinsX() + 1):
            row = []
            erow = []
            for iy in range(1, h.GetNbinsY() + 1):
                row.append(float(h.GetBinContent(ix, iy)))
                erow.append(float(h.GetBinError(ix, iy)) ** 2)
            payload["values"].append(row)
            payload["variances"].append(erow)
        print("JSON_BEGIN")
        print(json.dumps(payload))
        print("JSON_END")
        """
    )
    auth_sock = subprocess.run(
        ["launchctl", "getenv", "SSH_AUTH_SOCK"],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout.strip()
    env = os.environ.copy()
    if auth_sock:
        env["SSH_AUTH_SOCK"] = auth_sock
    command = [
        "ssh",
        login_host,
        (
            "ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null "
            f"{worker_host} 'python3 -'"
        ),
    ]
    proc = subprocess.run(
        command,
        input=remote_py,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=env,
        check=False,
        timeout=180,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"remote Fig24 extraction failed rc={proc.returncode}\n{proc.stdout}")
    if "JSON_BEGIN" not in proc.stdout or "JSON_END" not in proc.stdout:
        raise RuntimeError(f"remote Fig24 extraction did not emit JSON markers\n{proc.stdout}")
    payload = json.loads(proc.stdout.split("JSON_BEGIN", 1)[1].split("JSON_END", 1)[0].strip())
    payload["values"] = np.asarray(payload["values"], dtype=float)
    payload["variances"] = np.asarray(payload["variances"], dtype=float)
    payload["xedges"] = np.asarray(payload["xedges"], dtype=float)
    payload["yedges"] = np.asarray(payload["yedges"], dtype=float)
    return payload


def write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_ppg12_cutoff_points(path: Path) -> dict[tuple[float, int], dict[str, Any]]:
    lines = path.read_text().splitlines()
    if "POINTS_BEGIN" not in lines:
        raise RuntimeError(f"{path} missing POINTS_BEGIN")
    start = lines.index("POINTS_BEGIN") + 1
    out: dict[tuple[float, int], dict[str, Any]] = {}
    for row in csv.DictReader(lines[start:]):
        if row.get("kind") != "POINT":
            continue
        eff = float(row["eff"])
        idx = int(row["bin"])
        out[(eff, idx)] = row
    return out


def fig24_shape_rows(label: str, th2: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    values = th2["values"]
    variances = th2["variances"]
    xedges = th2["xedges"]
    yedges = th2["yedges"]
    for lo, hi in PT_SLICES:
        vals, vars_out = th2_slice(values, variances, xedges, yedges, lo, hi)
        stats = distribution_stats(yedges, vals, vars_out)
        row = {
            "source": label,
            "pt_low": lo,
            "pt_high": hi,
            "object": th2["object"],
            "root": th2["source"],
        }
        row.update(stats)
        rows.append(row)
    return rows


def fig24_cutoff_comparison_rows(cutoff_csv: Path, current_th2: dict[str, Any]) -> list[dict[str, Any]]:
    ppg12 = parse_ppg12_cutoff_points(cutoff_csv)
    rows: list[dict[str, Any]] = []
    values = current_th2["values"]
    variances = current_th2["variances"]
    xedges = current_th2["xedges"]
    yedges = current_th2["yedges"]
    pt_edges = np.linspace(10.0, 36.0, 14)
    for eff in EFFS:
        for idx in range(13):
            lo = float(pt_edges[idx])
            hi = float(pt_edges[idx + 1])
            vals, vars_out = th2_slice(values, variances, xedges, yedges, lo, hi)
            stats = distribution_stats(yedges, vals, vars_out)
            p = ppg12.get((eff, idx + 1))
            current_cut = stats[f"cut{int(round(100 * eff))}_isolow_minus1"]
            ppg12_cut = float(p["cut"]) if p else math.nan
            rows.append(
                {
                    "eff": eff,
                    "bin": idx + 1,
                    "pt_low": lo,
                    "pt_high": hi,
                    "pt_center": 0.5 * (lo + hi),
                    "ppg12_cut": ppg12_cut,
                    "current_cut": current_cut,
                    "sdcc_over_current": ppg12_cut / current_cut if current_cut and math.isfinite(ppg12_cut) else math.nan,
                    "current_entries": stats["entries"],
                    "current_mean": stats["mean"],
                    "current_rms": stats["rms"],
                }
            )
    return rows


def hist1d_stats_from_root(path: Path, key: str) -> dict[str, Any] | None:
    with uproot.open(path) as handle:
        if key not in handle:
            return None
        hist = handle[key]
        values, edges = hist.to_numpy(flow=False)
        variances = hist.variances(flow=False)
    stats = distribution_stats(edges.astype(float), values.astype(float), None if variances is None else variances.astype(float))
    return {"key": key, **stats}


def current_component_rows(root: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    patterns = [
        ("exact_fig24_corrected", "SIM/h_ppg12_fig24_EisoReco_truthSigMatched_pT_26_28"),
        ("exact_fig24_corrected", "SIM/h_ppg12_fig24_EisoReco_truthSigMatched_pT_28_30"),
        ("exact_fig24_corrected", "SIM/h_ppg12_fig24_EisoReco_truthSigMatched_pT_30_32"),
        ("exact_fig24_corrected", "SIM/h_ppg12_fig24_EisoReco_truthSigMatched_pT_32_34"),
        ("exact_fig24_corrected", "SIM/h_ppg12_fig24_EisoReco_truthSigMatched_pT_34_36"),
        ("legacy_corrected", "SIM/h_EisoReco_truthSigMatched_pT_26_35"),
        ("raw_all", "SIM/h_Eiso_ppg12_topo_raw_pT_26_35"),
        ("raw_tight", "SIM/h_Eiso_ppg12_topo_raw_tight_pT_26_35"),
        ("raw_nonTight", "SIM/h_Eiso_ppg12_topo_raw_nonTight_pT_26_35"),
    ]
    families = ["g4Stored", "g4NoEmbed", "hepmcAnyPhoton", "hepmcParentStable", "hepmcParent"]
    for fam in families:
        for vtx in ["vtxW", "noVtxW"]:
            patterns.append((f"stitch_{fam}_{vtx}", f"SIM/h_EisoReco_truthSigMatched_stitch_{fam}_{vtx}_pT_26_35"))
    for label, key in patterns:
        stats = hist1d_stats_from_root(root, key)
        if stats is None:
            rows.append({"source": label, "key": key, "status": "missing"})
        else:
            rows.append({"source": label, "status": "present", **stats})
    return rows


def fig3_tail_rows(path: Path) -> list[dict[str, Any]]:
    rows_raw = list(csv.DictReader(path.open()))
    windows = [(-0.5, 6.0), (0.0, 2.0), (2.0, 6.0), (6.0, 12.0), (12.0, 20.0)]
    out: list[dict[str, Any]] = []
    for lo, hi in windows:
        vals = []
        centers = []
        for row in rows_raw:
            try:
                center = float(row["bin_center"])
                ratio = float(row["current_over_ppg12"])
                ppg12 = float(row["ppg12_sdcc_counts_per_width"])
            except Exception:
                continue
            if lo <= center <= hi and ppg12 > 10.0 and math.isfinite(ratio):
                vals.append(ratio)
                centers.append(center)
        if vals:
            idx = int(np.argmax(np.abs(np.asarray(vals) - 1.0)))
            out.append(
                {
                    "window_low": lo,
                    "window_high": hi,
                    "n_bins": len(vals),
                    "mean_current_over_ppg12": float(np.mean(vals)),
                    "min_current_over_ppg12": float(np.min(vals)),
                    "max_current_over_ppg12": float(np.max(vals)),
                    "max_abs_deviation": float(np.max(np.abs(np.asarray(vals) - 1.0))),
                    "max_deviation_center": centers[idx],
                }
            )
        else:
            out.append({"window_low": lo, "window_high": hi, "n_bins": 0})
    return out


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True))


def write_source_chain_tables(outdir: Path, *, ppg12_fig24_remote_root: str, current_root: Path, fig3_csv: Path) -> None:
    rows = [
        {
            "target": "Fig24",
            "PPG12 ROOT": ppg12_fig24_remote_root,
            "object": "h_singal_reco_isoET_0",
            "macro": "ppg12codeGit/efficiencytool/FindETCut.C:7-105",
            "observable": "cluster E_T vs corrected reco topo Eiso",
            "selection": "truth-matched signal cluster, eta bin 0; see RecoEffCalculator_TTreeReader.C",
            "correction": "cluster_iso_topo_04 then recoisoET = recoisoET * mc_iso_scale + mc_iso_shift; config_bdt_nom.yaml uses 1.2, 0.1",
            "binning": "400 x bins 0-50, 4400 y bins -5-50; cutoff scan 13 bins 10-36 and isolow=-1",
            "evidence": "FindETCut.C:78-105; config_bdt_nom.yaml:61,81-83; RecoEffCalculator_TTreeReader.C:690,1213-1215,2151-2164,2261-2266,2916",
        },
        {
            "target": "Fig3",
            "PPG12 ROOT": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo_bdt_nom.root",
            "object": "h_tight_isoET_0_3 + h_tight_isoET_0_4 + h_tight_isoET_0_5",
            "macro": "ppg12codeGit/plotting/plot_isoET.C / plot_sideband.C family",
            "observable": "tight-data isolation template in 16-22 GeV pT window",
            "selection": "tight data clusters in eta bin 0, pT indexes 3-5",
            "correction": "data path uses the stored recoisoET chosen by the efficiency tool; MC scale/shift is only under issim",
            "binning": "source 0.1 GeV bins; local reproduction applies PPG12 variable rebinning",
            "evidence": f"{fig3_csv}; make_ppg12_fig3_iso_data_overlay.py:24-41,106-152",
        },
    ]
    write_rows(outdir / "ppg12_source_chain_table.csv", rows)
    md = ["# PPG12 Source Chain Table", ""]
    for row in rows:
        md.append(f"## {row['target']}")
        for key, value in row.items():
            if key != "target":
                md.append(f"- **{key}**: {value}")
        md.append("")
    (outdir / "ppg12_source_chain_table.md").write_text("\n".join(md))

    current_rows = [
        {
            "target": "Fig24",
            "current ROOT": str(current_root),
            "object": "SIM/h_singal_reco_isoET_0 and SIM/h_ppg12_fig24_signal_reco_isoET_eta0",
            "fill line": "src/RecoilJets.cc:9157,9180-9244",
            "raw isolation input": "ppg12PhotonYieldRawEiso(recoMatch, topNode)",
            "correction": "ppg12PhotonYieldEiso(raw) = 1.2*raw + 0.1 by default",
            "cluster selection": "findRecoPhotonMatchedToTruthSignal plus PhotonClusterv1 and tower mask, |eta|<0.7, Fig24 reco pT bin, photon-sample ownership",
            "truth/data selection": "truth-isolated signal loop, reco truth match",
            "weight": "Fig24 compatibility hist fill is unweighted in the dedicated fill call",
            "evidence": "src/RecoilJets.cc:9115-9157,9180-9244,13487-13627,19128-19153",
        },
        {
            "target": "Fig3",
            "current ROOT": "dataOutput/current_recoiljets_artifacts/current/pp_data_merged/current.root",
            "object": "PPG12_scaledtrigger30/h_tight_isoET_0_3,4,5",
            "fill line": "src/RecoilJets.cc tight-data iso template family, exact line not changed in this pass",
            "raw isolation input": "same PPG12 photon-yield iso path when enabled for pp PPG12 output",
            "correction": "data is not MC-shifted; the Fig3 path checks the stored/current tight-data Eiso template",
            "cluster selection": "tight data clusters in eta bin 0, pT indexes 3-5",
            "truth/data selection": "data tight template; no truth match",
            "weight": "data trigger/directory normalization handled in plotting script",
            "evidence": "make_ppg12_fig3_iso_data_overlay.py; existing ratio CSV/manifest",
        },
    ]
    write_rows(outdir / "current_source_chain_table.csv", current_rows)
    md = ["# Current Source Chain Table", ""]
    for row in current_rows:
        md.append(f"## {row['target']}")
        for key, value in row.items():
            if key != "target":
                md.append(f"- **{key}**: {value}")
        md.append("")
    (outdir / "current_source_chain_table.md").write_text("\n".join(md))


def write_cluster_ladder(outdir: Path) -> None:
    rows = [
        {
            "stage": "raw isolation branch/object",
            "PPG12 behavior": "Reads cluster_iso_topo_04_<clusternodename> when use_topo_iso == 2.",
            "RecoilJets behavior": "Uses PhotonClusterv1 stored ppg12_topo_raw_eiso_04 if vertex-compatible, else recomputes topo-cone sum from CLUSTERINFO_CEMC_NO_SPLIT-like topo cluster node and subtracts candidate ET.",
            "same/different/unknown": "different_or_fallback_dependent",
            "could cause high-pT Eiso tail?": "yes",
            "evidence": "RecoEffCalculator_TTreeReader.C:690,2151-2164; RecoilJets.cc:13487-13612",
        },
        {
            "stage": "MC scale/offset correction",
            "PPG12 behavior": "Applies recoisoET = recoisoET*1.2 + 0.1 for simulation in nominal config.",
            "RecoilJets behavior": "Applies ppg12PhotonYieldEiso(raw) = m_scale*raw + m_shift, defaults 1.2 and 0.1.",
            "same/different/unknown": "same_by_code_path_defaults",
            "could cause high-pT Eiso tail?": "unlikely unless exact artifact env overrode defaults",
            "evidence": "config_bdt_nom.yaml:81-82; RecoilJets.h:1621-1622; RecoilJets.cc:13614-13627",
        },
        {
            "stage": "truth-match/signal cluster",
            "PPG12 behavior": "Efficiency tool fills h_singal_reco_isoET for reconstructed clusters matched to signal truth particle after its cluster/photon loops.",
            "RecoilJets behavior": "findRecoPhotonMatchedToTruthSignal, PhotonClusterv1 dynamic cast, tower mask, eta/reco-pT/sample ownership gate before Fig24 fill.",
            "same/different/unknown": "unknown_exact_identity",
            "could cause high-pT Eiso tail?": "yes",
            "evidence": "RecoEffCalculator_TTreeReader.C:2860-2917; RecoilJets.cc:9115-9244",
        },
        {
            "stage": "component/weight mixture",
            "PPG12 behavior": "MC_efficiency_bdt_nom.root is the nominal merged efficiency tool product.",
            "RecoilJets behavior": "Current ROOT is downstream photon+jet merged artifact; dedicated Fig24 fill itself is unweighted, but the event population is the merged production population.",
            "same/different/unknown": "unknown_component_identity",
            "could cause high-pT Eiso tail?": "possible but not leading from current local variants",
            "evidence": "Fig24 variant localization and current ROOT key families",
        },
        {
            "stage": "cutoff extraction",
            "PPG12 behavior": "FindETCut.C sums TH2 pT interval and scans y from isolow=-1.",
            "RecoilJets behavior": "make_the76_fig24_iso_cut_overlay.py replicates that scan.",
            "same/different/unknown": "same",
            "could cause high-pT Eiso tail?": "no",
            "evidence": "FindETCut.C:7-68; make_the76_fig24_iso_cut_overlay.py:71-103",
        },
    ]
    write_rows(outdir / "cluster_selection_truth_match_ladder.csv", rows)
    md = ["# Cluster Selection / Truth-Match Ladder", ""]
    for row in rows:
        md.append(f"- **{row['stage']}**: {row['same/different/unknown']}. {row['PPG12 behavior']} / {row['RecoilJets behavior']} Evidence: {row['evidence']}")
    (outdir / "cluster_selection_truth_match_ladder.md").write_text("\n".join(md))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--current-root", type=Path, default=DEFAULT_CURRENT_ROOT)
    ap.add_argument("--fig24-cutoff-csv", type=Path, default=DEFAULT_FIG24_CUTOFF_CSV)
    ap.add_argument("--fig3-ratio-csv", type=Path, default=DEFAULT_FIG3_RATIO_CSV)
    ap.add_argument("--fig3-json", type=Path, default=DEFAULT_FIG3_JSON)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    ap.add_argument("--refresh-sdcc-fig24", action="store_true")
    ap.add_argument("--ppg12-fig24-remote-root", default=DEFAULT_PPG12_FIG24_REMOTE_ROOT)
    ap.add_argument("--ppg12-fig24-object", default=DEFAULT_PPG12_FIG24_OBJECT)
    ap.add_argument("--login-host", default="patsfan753@ssh.sdcc.bnl.gov")
    ap.add_argument("--worker-host", default="sphnxuser05.sdcc.bnl.gov")
    args = ap.parse_args()

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    cache_path = outdir / "ppg12_fig24_sdcc_h_singal_reco_isoET_0_payload.json"

    current_th2 = load_th2_from_root(args.current_root, "SIM/h_singal_reco_isoET_0")
    current_shape_rows = fig24_shape_rows("current_SIM_h_singal_reco_isoET_0", current_th2)
    write_rows(outdir / "fig24_current_eiso_shape_by_pt.csv", current_shape_rows)

    ppg12_th2: dict[str, Any] | None = None
    ppg12_extract_error: str | None = None
    try:
        if args.refresh_sdcc_fig24 or not cache_path.exists():
            payload = remote_extract_ppg12_fig24(
                args.ppg12_fig24_remote_root,
                args.ppg12_fig24_object,
                args.login_host,
                args.worker_host,
            )
            json_payload = {
                "source": payload["source"],
                "object": payload["object"],
                "xedges": payload["xedges"].tolist(),
                "yedges": payload["yedges"].tolist(),
                "values": payload["values"].tolist(),
                "variances": payload["variances"].tolist(),
            }
            write_json(cache_path, json_payload)
        else:
            json_payload = json.loads(cache_path.read_text())
            payload = {
                "source": json_payload["source"],
                "object": json_payload["object"],
                "xedges": np.asarray(json_payload["xedges"], dtype=float),
                "yedges": np.asarray(json_payload["yedges"], dtype=float),
                "values": np.asarray(json_payload["values"], dtype=float),
                "variances": np.asarray(json_payload["variances"], dtype=float),
            }
        ppg12_th2 = payload
    except Exception as exc:
        ppg12_extract_error = str(exc)

    shape_comparison_rows: list[dict[str, Any]] = []
    if ppg12_th2 is not None:
        ppg12_shape_rows = fig24_shape_rows("ppg12_sdcc_h_singal_reco_isoET_0", ppg12_th2)
        write_rows(outdir / "fig24_ppg12_sdcc_eiso_shape_by_pt.csv", ppg12_shape_rows)
        for ppg, cur in zip(ppg12_shape_rows, current_shape_rows):
            row = {
                "pt_low": cur["pt_low"],
                "pt_high": cur["pt_high"],
                "ppg12_entries": ppg["entries"],
                "current_entries": cur["entries"],
                "ppg12_mean": ppg["mean"],
                "current_mean": cur["mean"],
                "current_minus_ppg12_mean": cur["mean"] - ppg["mean"],
                "ppg12_rms": ppg["rms"],
                "current_rms": cur["rms"],
                "current_minus_ppg12_rms": cur["rms"] - ppg["rms"],
                "ppg12_cut70": ppg["cut70_isolow_minus1"],
                "current_cut70": cur["cut70_isolow_minus1"],
                "sdcc_over_current_cut70": ppg["cut70_isolow_minus1"] / cur["cut70_isolow_minus1"],
                "ppg12_cut80": ppg["cut80_isolow_minus1"],
                "current_cut80": cur["cut80_isolow_minus1"],
                "sdcc_over_current_cut80": ppg["cut80_isolow_minus1"] / cur["cut80_isolow_minus1"],
                "ppg12_cut90": ppg["cut90_isolow_minus1"],
                "current_cut90": cur["cut90_isolow_minus1"],
                "sdcc_over_current_cut90": ppg["cut90_isolow_minus1"] / cur["cut90_isolow_minus1"],
            }
            for eiso in FIXED_EISO_VALUES:
                key = f"cdf_le_{eiso:g}"
                row[f"ppg12_{key}"] = ppg[key]
                row[f"current_{key}"] = cur[key]
                row[f"current_minus_ppg12_{key}"] = cur[key] - ppg[key]
            shape_comparison_rows.append(row)
        write_rows(outdir / "fig24_eiso_shape_comparison_by_pt.csv", shape_comparison_rows)

    cutoff_rows = fig24_cutoff_comparison_rows(args.fig24_cutoff_csv, current_th2)
    write_rows(outdir / "fig24_cutoff_comparison_from_current_th2.csv", cutoff_rows)

    component_rows = current_component_rows(args.current_root)
    write_rows(outdir / "fig24_current_component_and_raw_localization.csv", component_rows)

    fig3_rows = fig3_tail_rows(args.fig3_ratio_csv)
    write_rows(outdir / "fig3_data_tail_summary.csv", fig3_rows)

    write_source_chain_tables(
        outdir,
        ppg12_fig24_remote_root=args.ppg12_fig24_remote_root,
        current_root=args.current_root,
        fig3_csv=args.fig3_ratio_csv,
    )
    write_cluster_ladder(outdir)

    raw_contract_rows = [
        {
            "source": "PPG12 Fig24",
            "raw observable": "cluster_iso_topo_04_<clusternodename>",
            "correction location": "RecoEffCalculator_TTreeReader.C after raw recoisoET selection, before h_singal_reco_isoET fill",
            "correction formula": "recoisoET = recoisoET * mc_iso_scale + mc_iso_shift",
            "applied count": "one for simulation nominal config",
            "evidence": "config_bdt_nom.yaml:61,81-82; RecoEffCalculator_TTreeReader.C:2151-2164,2261-2266,2916",
            "consequence": "No evidence that Fig24 cutoff residual is caused by a missing or double MC scale/shift.",
        },
        {
            "source": "RecoilJets Fig24",
            "raw observable": "PhotonClusterv1 ppg12_topo_raw_eiso_04 if valid, otherwise recomputed topo cone sum minus candidate ET",
            "correction location": "ppg12PhotonYieldEiso called before Fig24 fill",
            "correction formula": "return m_ppg12PhotonYieldMcIsoScale * raw + m_ppg12PhotonYieldMcIsoShift",
            "applied count": "one by code path defaults; exact artifact env override not found in current ROOT metadata",
            "evidence": "RecoilJets.cc:9157,13487-13627; RecoilJets.h:1621-1622",
            "consequence": "Remaining leading uncertainty is raw topo isolation or cluster selection, not the documented linear correction.",
        },
    ]
    write_rows(outdir / "raw_vs_corrected_eiso_contract.csv", raw_contract_rows)

    summary = {
        "audit": "THE-76 Fig24/Fig3 Eiso tail",
        "current_root": str(args.current_root),
        "ppg12_fig24_remote_root": args.ppg12_fig24_remote_root,
        "ppg12_fig24_object": args.ppg12_fig24_object,
        "ppg12_fig24_payload_cache": str(cache_path) if cache_path.exists() else None,
        "ppg12_fig24_extract_error": ppg12_extract_error,
        "fig24_high_pt_cutoff_summary": {
            "n_rows": len(cutoff_rows),
            "high_pt_min_sdcc_over_current": min(
                r["sdcc_over_current"] for r in cutoff_rows if r["pt_low"] >= 26.0 and math.isfinite(r["sdcc_over_current"])
            ),
            "high_pt_mean_sdcc_over_current": float(
                np.mean([r["sdcc_over_current"] for r in cutoff_rows if r["pt_low"] >= 26.0 and math.isfinite(r["sdcc_over_current"])])
            ),
        },
        "fig3_tail_summary": fig3_rows,
        "classification": "candidate: data and MC both indicate PhotonClusterBuilder isolation tail",
        "classification_basis": [
            "Fig24 current cutoffs are systematically larger than PPG12 at high pT, so current Eiso is harder/broader.",
            "The PPG12 and RecoilJets documented linear MC correction is nominally the same.",
            "Existing no-vtx/coarse variants do not remove the high-pT behavior.",
            "Existing Fig3 data template comparison shows current/PPG12 ratios above one in the high-isolation tail windows, consistent with a shared isolation-template tail issue.",
            "Current final ROOT lacks enough per-sample/period component histograms to prove a weighted-combination cause.",
        ],
        "canonicalization_allowed": False,
        "minimal_next_diagnostic": "Run a tiny same-cluster/topo-isolation canary on photon20 high-pT events: compare PPG12 cluster_iso_topo_04, RecoilJets stored ppg12_topo_raw_eiso_04, and RecoilJets recomputed topo sum before applying 1.2*raw+0.1.",
        "outputs": {
            "fig24_current_eiso_shape_by_pt": str(outdir / "fig24_current_eiso_shape_by_pt.csv"),
            "fig24_ppg12_sdcc_eiso_shape_by_pt": str(outdir / "fig24_ppg12_sdcc_eiso_shape_by_pt.csv") if ppg12_th2 is not None else None,
            "fig24_shape_comparison_by_pt": str(outdir / "fig24_eiso_shape_comparison_by_pt.csv") if ppg12_th2 is not None else None,
            "fig24_cutoff_comparison": str(outdir / "fig24_cutoff_comparison_from_current_th2.csv"),
            "component_localization": str(outdir / "fig24_current_component_and_raw_localization.csv"),
            "fig3_data_tail_summary": str(outdir / "fig3_data_tail_summary.csv"),
            "raw_vs_corrected_contract": str(outdir / "raw_vs_corrected_eiso_contract.csv"),
            "cluster_ladder": str(outdir / "cluster_selection_truth_match_ladder.md"),
        },
    }
    write_json(outdir / "fig24_fig3_eiso_tail_summary.json", summary)

    md = [
        "# THE-76 Fig24/Fig3 Eiso Tail Audit",
        "",
        f"- Current ROOT: `{args.current_root}`",
        f"- PPG12 Fig24 source: `{args.ppg12_fig24_remote_root}:{args.ppg12_fig24_object}`",
        f"- PPG12 Fig24 extraction error: `{ppg12_extract_error}`" if ppg12_extract_error else "- PPG12 Fig24 TH2 extracted read-only from SDCC/cache.",
        f"- Classification: `{summary['classification']}`",
        f"- Canonicalization allowed: `{summary['canonicalization_allowed']}`",
        "",
        "## Key Evidence",
        "",
    ]
    for item in summary["classification_basis"]:
        md.append(f"- {item}")
    md += [
        "",
        "## Minimal Next Diagnostic",
        "",
        summary["minimal_next_diagnostic"],
        "",
    ]
    (outdir / "fig24_fig3_eiso_tail_summary.md").write_text("\n".join(md))

    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
