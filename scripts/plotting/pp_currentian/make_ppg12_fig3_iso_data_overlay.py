#!/usr/bin/env python3
"""Overlay PPG12 Fig. 3 isolation-template data with current pp output.

This script compares the tight-data isolation distribution used in the PPG12
paper-draft Fig. 3 contract against the current RecoilJets pp data output for
the same pT bin, using the same variable rebinning and bin-width scaling from
``ppg12codeGit/plotting/CONF_plots.C``.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import uproot


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
CAMPAIGN = "the76_ppg12_parity_full_20260701_003024"
DEFAULT_CURRENT_ROOT = REPO / "dataOutput/current_recoiljets_artifacts/current/pp_data_merged/current.root"
DEFAULT_OUTDIR = REPO / "dataOutput/ppg12Parity" / CAMPAIGN / "data_iso_template_fig3"
DEFAULT_REMOTE_PPG12_ROOT = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo_bdt_nom.root"
DEFAULT_LOGIN_HOST = "patsfan753@ssh.sdcc.bnl.gov"
DEFAULT_WORKER_HOST = "sphnxuser05.sdcc.bnl.gov"
DEFAULT_CURRENT_LABEL = "Current pp data"

PT_INDEXES = (3, 4, 5)
PT_LABEL = r"$16 < E_T^\gamma < 22\ \mathrm{GeV}$"
HIST_BASE = "h_tight_isoET_0_"


@dataclass(frozen=True)
class HistSeries:
    label: str
    source: str
    keys: list[str]
    raw_values: np.ndarray
    raw_variances: np.ndarray
    raw_edges: np.ndarray
    values: np.ndarray
    variances: np.ndarray
    edges: np.ndarray
    density: np.ndarray
    density_variance: np.ndarray

    @property
    def centers(self) -> np.ndarray:
        return self.edges[:-1] + 0.5 * np.diff(self.edges)

    @property
    def density_errors(self) -> np.ndarray:
        return np.sqrt(np.clip(self.density_variance, 0.0, None))

    @property
    def raw_integral(self) -> float:
        return float(np.sum(self.raw_values))


def read_current_series(path: Path, directory: str, label: str) -> HistSeries:
    values_sum: np.ndarray | None = None
    variances_sum: np.ndarray | None = None
    edges_ref: np.ndarray | None = None
    keys: list[str] = []
    with uproot.open(path) as fh:
        for idx in PT_INDEXES:
            key = f"{directory}/{HIST_BASE}{idx}"
            if key not in fh:
                raise KeyError(f"missing current histogram {key} in {path}")
            hist = fh[key]
            values, edges = hist.to_numpy(flow=False)
            variances = hist.variances(flow=False)
            if variances is None:
                variances = np.clip(values, 0.0, None)
            values = values.astype(float)
            variances = variances.astype(float)
            edges = edges.astype(float)
            if values_sum is None:
                values_sum = np.zeros_like(values, dtype=float)
                variances_sum = np.zeros_like(variances, dtype=float)
                edges_ref = edges
            elif len(edges) != len(edges_ref) or float(np.max(np.abs(edges - edges_ref))) > 1.0e-9:
                raise ValueError(f"incompatible current binning for {key}")
            values_sum += values
            variances_sum += variances
            keys.append(key)
    assert values_sum is not None and variances_sum is not None and edges_ref is not None
    return make_processed_series(
        label=label,
        source=str(path),
        keys=keys,
        raw_values=values_sum,
        raw_variances=variances_sum,
        raw_edges=edges_ref,
    )


def extract_remote_ppg12_payload(remote_root: str, login_host: str, worker_host: str) -> dict[str, Any]:
    remote_py = textwrap.dedent(
        f"""
        import json
        import ROOT
        ROOT.gROOT.SetBatch(True)
        path = {remote_root!r}
        f = ROOT.TFile.Open(path)
        if not f or f.IsZombie():
            raise SystemExit("could not open " + path)
        payload = {{"source_root": path, "hist_base": {HIST_BASE!r}, "pt_indexes": {list(PT_INDEXES)!r}, "hists": []}}
        for idx in payload["pt_indexes"]:
            key = payload["hist_base"] + str(idx)
            h = f.Get(key)
            if not h:
                raise SystemExit("missing " + key)
            values = []
            errors = []
            edges = []
            for ib in range(1, h.GetNbinsX() + 1):
                if ib == 1:
                    edges.append(float(h.GetXaxis().GetBinLowEdge(ib)))
                values.append(float(h.GetBinContent(ib)))
                errors.append(float(h.GetBinError(ib)))
                edges.append(float(h.GetXaxis().GetBinLowEdge(ib) + h.GetXaxis().GetBinWidth(ib)))
            payload["hists"].append({{"key": key, "values": values, "errors": errors, "edges": edges}})
        print("JSON_BEGIN")
        print(json.dumps(payload))
        print("JSON_END")
        """
    )
    command = [
        "ssh",
        login_host,
        (
            "ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null "
            f"{worker_host} 'python3 -'"
        ),
    ]
    env = os.environ.copy()
    auth_sock = subprocess.run(
        ["launchctl", "getenv", "SSH_AUTH_SOCK"],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout.strip()
    if auth_sock:
        env["SSH_AUTH_SOCK"] = auth_sock
    proc = subprocess.run(
        command,
        input=remote_py,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=env,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"remote SDCC extraction failed with rc={proc.returncode}\n{proc.stdout}")
    if "JSON_BEGIN" not in proc.stdout or "JSON_END" not in proc.stdout:
        raise RuntimeError(f"remote SDCC extraction did not emit JSON markers\n{proc.stdout}")
    return json.loads(proc.stdout.split("JSON_BEGIN", 1)[1].split("JSON_END", 1)[0].strip())


def read_sdcc_series(json_path: Path, *, refresh: bool, remote_root: str, login_host: str, worker_host: str) -> HistSeries:
    if refresh or not json_path.exists():
        payload = extract_remote_ppg12_payload(remote_root, login_host, worker_host)
        json_path.parent.mkdir(parents=True, exist_ok=True)
        json_path.write_text(json.dumps(payload, indent=2, sort_keys=True))
    else:
        payload = json.loads(json_path.read_text())
    values_sum: np.ndarray | None = None
    variances_sum: np.ndarray | None = None
    edges_ref: np.ndarray | None = None
    keys: list[str] = []
    for hist in payload["hists"]:
        values = np.asarray(hist["values"], dtype=float)
        variances = np.square(np.asarray(hist["errors"], dtype=float))
        edges = np.asarray(hist["edges"], dtype=float)
        if values_sum is None:
            values_sum = np.zeros_like(values, dtype=float)
            variances_sum = np.zeros_like(variances, dtype=float)
            edges_ref = edges
        elif len(edges) != len(edges_ref) or float(np.max(np.abs(edges - edges_ref))) > 1.0e-9:
            raise ValueError(f"incompatible SDCC binning for {hist['key']}")
        values_sum += values
        variances_sum += variances
        keys.append(hist["key"])
    assert values_sum is not None and variances_sum is not None and edges_ref is not None
    return make_processed_series(
        label="PPG12 SDCC ROOT data",
        source=str(payload["source_root"]),
        keys=keys,
        raw_values=values_sum,
        raw_variances=variances_sum,
        raw_edges=edges_ref,
    )


def variable_rebin(values: np.ndarray, variances: np.ndarray, edges: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    new_edges = [float(edges[0])]
    out_values: list[float] = []
    out_variances: list[float] = []
    i = 0
    while i < len(values):
        group = 1 if edges[i] < 2.5 else 5
        j = min(i + group, len(values))
        out_values.append(float(np.sum(values[i:j])))
        out_variances.append(float(np.sum(variances[i:j])))
        new_edges.append(float(edges[j]))
        i = j
    return np.asarray(out_values), np.asarray(out_variances), np.asarray(new_edges)


def make_processed_series(
    *,
    label: str,
    source: str,
    keys: list[str],
    raw_values: np.ndarray,
    raw_variances: np.ndarray,
    raw_edges: np.ndarray,
) -> HistSeries:
    values, variances, edges = variable_rebin(raw_values, raw_variances, raw_edges)
    widths = np.diff(edges)
    density = values / widths
    density_variance = variances / np.square(widths)
    return HistSeries(
        label=label,
        source=source,
        keys=keys,
        raw_values=raw_values,
        raw_variances=raw_variances,
        raw_edges=raw_edges,
        values=values,
        variances=variances,
        edges=edges,
        density=density,
        density_variance=density_variance,
    )


def ratio_payload(current: HistSeries, sdcc: HistSeries) -> dict[str, np.ndarray | float]:
    if len(current.edges) != len(sdcc.edges) or float(np.max(np.abs(current.edges - sdcc.edges))) > 1.0e-9:
        raise ValueError("current and SDCC rebinned edges do not match")
    mask = sdcc.density > 0
    ratio = np.full_like(sdcc.density, np.nan, dtype=float)
    ratio_err = np.full_like(sdcc.density, np.nan, dtype=float)
    ratio[mask] = current.density[mask] / sdcc.density[mask]
    cur_rel = np.divide(
        current.density_errors,
        current.density,
        out=np.full_like(current.density, np.nan, dtype=float),
        where=current.density > 0,
    )
    sdcc_rel = np.divide(
        sdcc.density_errors,
        sdcc.density,
        out=np.full_like(sdcc.density, np.nan, dtype=float),
        where=sdcc.density > 0,
    )
    ratio_err[mask] = ratio[mask] * np.sqrt(np.square(cur_rel[mask]) + np.square(sdcc_rel[mask]))
    finite = np.isfinite(ratio)
    max_dev = float(np.nanmax(np.abs(ratio[finite] - 1.0))) if np.any(finite) else math.nan
    max_idx = int(np.nanargmax(np.abs(ratio - 1.0))) if np.any(finite) else -1
    main_mask = finite & (current.centers >= 0.0) & (current.centers <= 6.0) & (sdcc.density > 10.0)
    main_dev = float(np.nanmax(np.abs(ratio[main_mask] - 1.0))) if np.any(main_mask) else math.nan
    main_idx = int(np.nanargmax(np.where(main_mask, np.abs(ratio - 1.0), np.nan))) if np.any(main_mask) else -1
    main_mean = float(np.nanmean(ratio[main_mask])) if np.any(main_mask) else math.nan
    return {
        "ratio": ratio,
        "ratio_err": ratio_err,
        "max_abs_deviation": max_dev,
        "max_abs_deviation_percent": 100.0 * max_dev if math.isfinite(max_dev) else math.nan,
        "max_deviation_center": float(current.centers[max_idx]) if max_idx >= 0 else math.nan,
        "main_window": "0 <= Eiso <= 6 GeV and PPG12 density > 10",
        "main_mean_ratio": main_mean,
        "main_max_abs_deviation": main_dev,
        "main_max_abs_deviation_percent": 100.0 * main_dev if math.isfinite(main_dev) else math.nan,
        "main_max_deviation_center": float(current.centers[main_idx]) if main_idx >= 0 else math.nan,
    }


def write_csv(path: Path, current: HistSeries, sdcc: HistSeries, ratios: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "bin_low",
                "bin_high",
                "bin_center",
                "ppg12_sdcc_counts_per_width",
                "ppg12_sdcc_error",
                "current_counts_per_width",
                "current_error",
                "current_over_ppg12",
                "current_over_ppg12_error",
            ]
        )
        for lo, hi, x, s, se, c, ce, r, re in zip(
            current.edges[:-1],
            current.edges[1:],
            current.centers,
            sdcc.density,
            sdcc.density_errors,
            current.density,
            current.density_errors,
            ratios["ratio"],
            ratios["ratio_err"],
        ):
            writer.writerow([lo, hi, x, s, se, c, ce, r, re])


def draw_overlay(path: Path, current: HistSeries, sdcc: HistSeries, ratios: dict[str, Any], *, current_dir: str) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "mathtext.fontset": "dejavusans",
            "axes.linewidth": 1.45,
            "xtick.major.width": 1.35,
            "ytick.major.width": 1.35,
            "xtick.minor.width": 1.0,
            "ytick.minor.width": 1.0,
        }
    )
    fig = plt.figure(figsize=(7.72, 9.98), dpi=160)
    gs = fig.add_gridspec(2, 1, height_ratios=[3.1, 1.0], hspace=0.035)
    ax = fig.add_subplot(gs[0])
    rax = fig.add_subplot(gs[1], sharex=ax)

    ax.errorbar(
        sdcc.centers,
        sdcc.density,
        yerr=sdcc.density_errors,
        fmt="o",
        ms=5.2,
        color="black",
        ecolor="black",
        elinewidth=1.1,
        capsize=0,
        label="PPG12 SDCC ROOT data",
        zorder=3,
    )
    ax.errorbar(
        current.centers,
        current.density,
        yerr=current.density_errors,
        fmt="s",
        ms=4.7,
        mfc="white",
        mec="#1f77b4",
        mew=1.2,
        color="#1f77b4",
        ecolor="#58a5df",
        elinewidth=1.0,
        capsize=0,
        label=current.label,
        zorder=4,
    )

    ymax = max(float(np.nanmax(sdcc.density + sdcc.density_errors)), float(np.nanmax(current.density + current.density_errors)))
    ax.set_xlim(-1.0, 15.0)
    ax.set_ylim(0.0, ymax * 1.22)
    ax.set_ylabel("Counts / Bin Width", fontsize=22)
    ax.tick_params(axis="both", which="both", direction="in", top=True, right=True, labelsize=17, length=7)
    ax.tick_params(axis="both", which="minor", length=4)
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    plt.setp(ax.get_xticklabels(), visible=False)

    ax.text(0.93, 0.93, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax.transAxes, ha="right", va="top", fontsize=22)
    ax.text(0.93, 0.855, r"$p{+}p\ \sqrt{s}=200\ \mathrm{GeV}$", transform=ax.transAxes, ha="right", va="top", fontsize=20)
    ax.text(0.93, 0.785, PT_LABEL, transform=ax.transAxes, ha="right", va="top", fontsize=20)
    ax.text(0.93, 0.715, r"$|\eta^\gamma| < 0.7$", transform=ax.transAxes, ha="right", va="top", fontsize=20)

    legend = ax.legend(loc="upper right", bbox_to_anchor=(0.93, 0.63), frameon=False, fontsize=17, handlelength=1.5)
    for lh in legend.legend_handles:
        try:
            lh.set_linewidth(1.3)
        except Exception:
            pass

    stats_text = (
        f"PPG12 N = {sdcc.raw_integral:.1f}\n"
        f"Current N = {current.raw_integral:.0f}\n"
        f"main max |R-1| = {ratios['main_max_abs_deviation_percent']:.1f}%\n"
        f"near $E_T^{{iso}}$ = {ratios['main_max_deviation_center']:.2g} GeV"
    )
    ax.text(0.56, 0.47, stats_text, transform=ax.transAxes, ha="left", va="top", fontsize=16)
    rax.axhline(1.0, color="gray", linestyle=(0, (4, 4)), linewidth=1.1)
    rax.errorbar(
        current.centers,
        ratios["ratio"],
        yerr=ratios["ratio_err"],
        fmt="o",
        ms=4.3,
        color="#1f77b4",
        ecolor="#58a5df",
        elinewidth=1.0,
        capsize=0,
    )
    finite = np.asarray(ratios["ratio"])[np.isfinite(ratios["ratio"])]
    if finite.size:
        lo = float(np.nanpercentile(finite, 2)) - 0.18
        hi = float(np.nanpercentile(finite, 98)) + 0.18
        lo = min(lo, 0.65)
        hi = max(hi, 1.35)
        if hi - lo > 2.5:
            lo, hi = 0.0, min(max(float(np.nanmax(finite)) * 1.12, 1.6), 4.0)
        rax.set_ylim(max(0.0, lo), hi)
    else:
        rax.set_ylim(0.5, 1.5)
    rax.set_xlabel(r"$E_T^{\mathrm{iso,reco}}\ \mathrm{[GeV]}$", fontsize=22, loc="right")
    rax.set_ylabel("Current / PPG12", fontsize=16)
    rax.tick_params(axis="both", which="both", direction="in", top=True, right=True, labelsize=16, length=7)
    rax.tick_params(axis="both", which="minor", length=4)
    rax.xaxis.set_minor_locator(AutoMinorLocator(5))
    rax.yaxis.set_minor_locator(AutoMinorLocator(5))
    fig.subplots_adjust(left=0.13, right=0.965, bottom=0.095, top=0.965)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)


def write_manifest(path: Path, current: HistSeries, sdcc: HistSeries, ratios: dict[str, Any], args: argparse.Namespace, csv_path: Path, png_path: Path, sdcc_json: Path) -> None:
    manifest = {
        "campaign": args.campaign_tag,
        "comparison": "PPG12 Fig. 3 tight-data isolation template, SDCC ROOT vs current RecoilJets pp data",
        "pt_bin": "16 < E_T^gamma < 22 GeV",
        "eta_selection": "|eta^gamma| < 0.7",
        "ppg12_reference_code": "ppg12codeGit/plotting/CONF_plots.C",
        "processing_contract": {
            "summed_pt_indexes": list(PT_INDEXES),
            "histogram_base": HIST_BASE,
            "variable_rebin": "group size 1 for original bins with low edge < 2.5 GeV, else group size 5",
            "scaling": "counts divided by rebinned bin width",
            "ratio": "current density / PPG12 SDCC density",
        },
        "ppg12_sdcc": {
            "source_root": sdcc.source,
            "keys": sdcc.keys,
            "sdcc_json": str(sdcc_json),
            "raw_integral": sdcc.raw_integral,
        },
        "current": {
            "source_root": current.source,
            "label": current.label,
            "directory": args.current_dir,
            "keys": current.keys,
            "raw_integral": current.raw_integral,
            "artifact_pointer": str(DEFAULT_CURRENT_ROOT),
        },
        "metrics": {
            "max_abs_deviation_percent": ratios["max_abs_deviation_percent"],
            "max_deviation_center": ratios["max_deviation_center"],
            "main_window": ratios["main_window"],
            "main_mean_ratio": ratios["main_mean_ratio"],
            "main_max_abs_deviation_percent": ratios["main_max_abs_deviation_percent"],
            "main_max_deviation_center": ratios["main_max_deviation_center"],
        },
        "outputs": {
            "png": str(png_path),
            "csv": str(csv_path),
            "manifest": str(path),
        },
    }
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--current-root", type=Path, default=DEFAULT_CURRENT_ROOT)
    parser.add_argument("--current-dir", default="PPG12_scaledtrigger30")
    parser.add_argument("--current-label", default=DEFAULT_CURRENT_LABEL)
    parser.add_argument("--campaign-tag", default=CAMPAIGN)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--remote-ppg12-root", default=DEFAULT_REMOTE_PPG12_ROOT)
    parser.add_argument("--login-host", default=DEFAULT_LOGIN_HOST)
    parser.add_argument("--worker-host", default=DEFAULT_WORKER_HOST)
    parser.add_argument("--refresh-sdcc", action="store_true")
    parser.add_argument("--tag", default=None, help="output tag suffix; defaults from current-dir")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    tag = args.tag or args.current_dir.replace("/", "_")
    args.outdir.mkdir(parents=True, exist_ok=True)
    sdcc_json = args.outdir / "ppg12_sdcc_fig3_iso_template_pt16_22_tight_data.json"
    sdcc = read_sdcc_series(
        sdcc_json,
        refresh=args.refresh_sdcc,
        remote_root=args.remote_ppg12_root,
        login_host=args.login_host,
        worker_host=args.worker_host,
    )
    current = read_current_series(args.current_root.resolve(), args.current_dir, args.current_label)
    ratios = ratio_payload(current, sdcc)

    stem = f"ppg12_fig3_iso_template_pt16_22_sdcc_vs_current_{tag}_ratio"
    png_path = args.outdir / f"{stem}.png"
    csv_path = args.outdir / f"{stem}.csv"
    manifest_path = args.outdir / f"{stem}_manifest.json"
    write_csv(csv_path, current, sdcc, ratios)
    draw_overlay(png_path, current, sdcc, ratios, current_dir=args.current_dir)
    write_manifest(manifest_path, current, sdcc, ratios, args, csv_path, png_path, sdcc_json)
    print(json.dumps({"png": str(png_path), "csv": str(csv_path), "manifest": str(manifest_path)}, indent=2))


if __name__ == "__main__":
    main()
