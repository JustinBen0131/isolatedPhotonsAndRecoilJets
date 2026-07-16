#!/usr/bin/env python3
"""Render a PPG12 IAN Fig. 18-style 0-mrad inclusive-MC NPB SI/DI overlay.

The PPG12 source macro, ``plot_npb_score_di.C``, sums the NPB-score spectra
of jet8/12/20/30/40 independently for SI and DI, then unit-normalizes each
interaction component.  This helper extracts the equivalent unweighted
TableQA score spectra from the already merged THE97 source components on SDCC
and applies the PPG12 per-sample cross-section/event-count factors locally.

The TableQA family intentionally stores raw candidate counts.  Its source
histogram is therefore a source-level parity diagnostic: PPG12's candidate
level truth-vertex reweight cannot be recovered from an already aggregated raw
score histogram.  That limitation is recorded in the manifest and never
hidden in an absolute parity claim.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import LogLocator, MultipleLocator


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
DEFAULT_REMOTE_BASE = (
    "/sphenix/u/patsfan753/scratch/thesisAnalysis/runs/recoiljets/current/"
    "the97_ppg12_final_parity_full_20260709_2230/merge_components"
)
DEFAULT_PPG12_RESULTS = "/sphenix/user/shuhangli/ppg12/efficiencytool/results"
DEFAULT_OUTDIR = REPO / "dataOutput/ppg12Parity/the97_ppg12_final_parity_full_20260709_2230/fig18_npb_si_di"
CFG = "jetMinPtScan_dphiScan_vz60_isoR40_isSliding_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12"
SAMPLES = ("jet8", "jet12", "jet20", "jet30", "jet40")

# CrossSectionWeights.h / PPG12 ShowerShapeCheck.C convention, normalized to
# jet50.  The common jet50 denominator cancels after each component is unit
# normalized, but retaining it makes the source convention explicit.
XSEC_PB = {
    "jet8": 1.15e7,
    "jet12": 1.4903e6,
    "jet20": 6.2623e4,
    "jet30": 2.5298e3,
    "jet40": 1.3553e2,
}
JET50_PB = 7.3113
MARKER = "__PPG12_FIG18_PAYLOAD__"


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--remote-base", default=DEFAULT_REMOTE_BASE)
    ap.add_argument("--ppg12-results", default=DEFAULT_PPG12_RESULTS)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    ap.add_argument("--output-name", default="ppg12_fig18_current_0mrad_npb_si_di_overlay.png")
    ap.add_argument(
        "--with-ppg12-reference",
        action="store_true",
        help="Read the exact SDCC PPG12 source histograms and render an SDCC/current ratio panel.",
    )
    return ap.parse_args()


def remote_payload(remote_base: str, ppg12_results: str, include_reference: bool) -> dict:
    """Extract only the compact 100-bin source histograms over SSH."""
    remote_code = f'''import json
import numpy as np
import uproot

base = {remote_base!r}
ppg12_results = {ppg12_results!r}
include_reference = {include_reference!r}
cfg = {CFG!r}
samples = {SAMPLES!r}
payload = {{"components": {{}}, "source_roots": {{}}}}
for component, suffix in (("si", ""), ("di", "_double")):
    out = {{}}
    roots = {{}}
    for sample in samples:
        root = (f"{{base}}/inclusivejet_0mrad_{{component}}/siminclusive/"
                f"RecoilJets_{{sample}}{{suffix}}_ALL_{{cfg}}.root")
        with uproot.open(root) as fin:
            hsum = None
            edges = None
            for ptbin in range(5):
                hist = fin[f"SIM/h1d_npb_score_eta0_pt{{ptbin}}_cut0"]
                values, this_edges = hist.to_numpy(flow=False)
                hsum = np.asarray(values, dtype=float) if hsum is None else hsum + values
                edges = np.asarray(this_edges, dtype=float)
            meta = fin["SIM/h_ppInclusiveJetStitch_ppg12Fig6EventWeighted_metadata"]
            mvals, medges = meta.to_numpy(flow=False)
        out[sample] = {{"values": hsum.tolist(), "events_processed": float(mvals[0])}}
        roots[sample] = root
    payload["components"][component] = {{"edges": edges.tolist(), "samples": out}}
    payload["source_roots"][component] = roots
if include_reference:
    reference = {{}}
    reference_roots = {{}}
    for component, mix in (("si", "nom"), ("di", "double")):
        out = {{}}
        roots = {{}}
        for sample in samples:
            root = f"{{ppg12_results}}/MC_efficiencyshower_shape_{{sample}}_{{mix}}_inclusive_showershape_0rad.root"
            with uproot.open(root) as fin:
                hsum = None
                edges = None
                for ptbin in range(5):
                    for bdtbin in range(3):
                        hist = fin[f"h_npb_score_eta0_pt{{ptbin}}_bdt{{bdtbin}}"]
                        values, this_edges = hist.to_numpy(flow=False)
                        hsum = np.asarray(values, dtype=float) if hsum is None else hsum + values
                        edges = np.asarray(this_edges, dtype=float)
            out[sample] = hsum.tolist()
            roots[sample] = root
        reference[component] = {{"edges": edges.tolist(), "samples": out}}
        reference_roots[component] = roots
    payload["ppg12_reference"] = reference
    payload["ppg12_reference_roots"] = reference_roots
print({MARKER!r})
print(json.dumps(payload, separators=(",", ":")))
'''
    sock = os.environ.get("SSH_AUTH_SOCK", "")
    if not sock:
        launchctl = subprocess.run(
            ["launchctl", "getenv", "SSH_AUTH_SOCK"], text=True, capture_output=True, check=False
        )
        sock = launchctl.stdout.strip()
    if not sock:
        raise RuntimeError("SSH_AUTH_SOCK is unavailable; cannot perform the read-only SDCC extraction")

    inner = (
        "ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null "
        "-o BatchMode=yes sphnxuser05.sdcc.bnl.gov " + shlex.quote("python3 - <<'PY'\n" + remote_code + "\nPY")
    )
    env = os.environ.copy()
    env["SSH_AUTH_SOCK"] = sock
    run = subprocess.run(
        ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=20", "patsfan753@ssh.sdcc.bnl.gov", inner],
        text=True,
        capture_output=True,
        check=False,
        env=env,
    )
    if run.returncode != 0:
        raise RuntimeError(f"Read-only SDCC extraction failed:\n{run.stderr[-2000:]}")
    if MARKER not in run.stdout:
        raise RuntimeError("SDCC extraction returned no marked JSON payload")
    return json.loads(run.stdout.rsplit(MARKER, 1)[1].strip())


def weighted_component(payload: dict, component: str) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    item = payload["components"][component]
    edges = np.asarray(item["edges"], dtype=float)
    total = np.zeros(len(edges) - 1, dtype=float)
    scales: dict[str, float] = {}
    for sample in SAMPLES:
        entry = item["samples"][sample]
        events = float(entry["events_processed"])
        if not np.isfinite(events) or events <= 0:
            raise RuntimeError(f"Invalid processed-event metadata for {component}/{sample}: {events}")
        scale = XSEC_PB[sample] / JET50_PB / events
        total += scale * np.asarray(entry["values"], dtype=float)
        scales[sample] = scale
    area = float(np.sum(total))
    if not np.isfinite(area) or area <= 0:
        raise RuntimeError(f"Empty weighted NPB spectrum for {component}")
    return edges, total / area, scales


def ppg12_component(payload: dict, component: str) -> tuple[np.ndarray, np.ndarray]:
    item = payload["ppg12_reference"][component]
    edges = np.asarray(item["edges"], dtype=float)
    total = np.zeros(len(edges) - 1, dtype=float)
    for sample in SAMPLES:
        total += np.asarray(item["samples"][sample], dtype=float)
    area = float(np.sum(total))
    if not np.isfinite(area) or area <= 0:
        raise RuntimeError(f"Empty PPG12 reference NPB spectrum for {component}")
    return edges, total / area


def render(edges: np.ndarray, si: np.ndarray, di: np.ndarray, output: Path, delta: float) -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "axes.linewidth": 1.1,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
    })
    fig, ax = plt.subplots(figsize=(6.2, 5.1), dpi=180)
    ax.stairs(si, edges, color="#1f4dff", linewidth=1.4, label="Single Interaction (current)")
    ax.stairs(di, edges, color="#ff4a4a", linewidth=1.4, linestyle=(0, (1.2, 1.8)), label="Double Interaction (current)")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(1e-4, 1.0)
    ax.set_yscale("log")
    ax.set_xlabel("NPB Score", fontsize=11.5)
    ax.set_ylabel("Normalized Counts", fontsize=11.5)
    ax.xaxis.set_major_locator(MultipleLocator(0.1))
    ax.xaxis.set_minor_locator(MultipleLocator(0.05))
    ax.yaxis.set_major_locator(LogLocator(base=10, numticks=5))
    ax.tick_params(labelsize=10)
    ax.text(0.03, 0.97, r"$\bf{\it{sPHENIX}}$ Internal", transform=ax.transAxes,
            ha="left", va="top", fontsize=10)
    ax.text(0.03, 0.90, "Inclusive MC (PYTHIA8)", transform=ax.transAxes,
            ha="left", va="top", fontsize=9.5)
    ax.text(0.03, 0.83, r"$|\eta| < 0.7$", transform=ax.transAxes,
            ha="left", va="top", fontsize=9.5)
    ax.text(0.03, 0.76, "0 mrad crossing", transform=ax.transAxes,
            ha="left", va="top", fontsize=9.5)
    ax.text(0.03, 0.69, rf"$\Delta$(DI$-$SI) below NPB<0.5: {delta:.3f}", transform=ax.transAxes,
            ha="left", va="top", fontsize=8.8)
    leg = ax.legend(loc="upper right", frameon=False, fontsize=9.3, handlelength=2.5)
    for line in leg.get_lines():
        line.set_linewidth(1.5)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output, facecolor="white")
    plt.close(fig)


def render_comparison(
    edges: np.ndarray,
    current_si: np.ndarray,
    current_di: np.ndarray,
    ppg12_si: np.ndarray,
    ppg12_di: np.ndarray,
    output: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "axes.linewidth": 1.1,
        "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
    })
    fig, (ax, ratio_ax) = plt.subplots(
        2, 1, figsize=(6.4, 6.3), dpi=180, sharex=True,
        gridspec_kw={"height_ratios": (3.0, 1.0), "hspace": 0.05},
    )
    # Color identifies SI/DI; line style identifies the source. This makes the
    # four curves readable without visually changing PPG12's SI/DI story.
    ax.stairs(ppg12_si, edges, color="#174ea6", linewidth=1.2, linestyle="--", label="PPG12 SDCC SI")
    ax.stairs(current_si, edges, color="#1f4dff", linewidth=1.55, label="Current SI")
    ax.stairs(ppg12_di, edges, color="#b31412", linewidth=1.2, linestyle="--", label="PPG12 SDCC DI")
    ax.stairs(current_di, edges, color="#ff4a4a", linewidth=1.55, linestyle=(0, (1.2, 1.8)), label="Current DI")
    ax.set_yscale("log")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(1e-4, 1.0)
    ax.set_ylabel("Normalized Counts", fontsize=11.5)
    ax.yaxis.set_major_locator(LogLocator(base=10, numticks=5))
    ax.tick_params(labelsize=9.5)
    ax.text(0.03, 0.96, r"$\bf{\it{sPHENIX}}$ Internal", transform=ax.transAxes,
            ha="left", va="top", fontsize=10)
    ax.text(0.03, 0.88, "Inclusive MC (PYTHIA8), 0 mrad", transform=ax.transAxes,
            ha="left", va="top", fontsize=9.3)
    ax.text(0.03, 0.80, r"$|\eta| < 0.7$, $10 < E_T^{\gamma,\mathrm{reco}} < 30$ GeV", transform=ax.transAxes,
            ha="left", va="top", fontsize=8.8)
    ax.legend(loc="upper right", frameon=False, fontsize=8.6, ncol=2, handlelength=2.5,
              columnspacing=1.2)

    centers = 0.5 * (edges[:-1] + edges[1:])
    valid_si = (current_si >= 1e-4) & (ppg12_si >= 1e-4)
    valid_di = (current_di >= 1e-4) & (ppg12_di >= 1e-4)
    ratio_si = np.divide(ppg12_si, current_si, out=np.full_like(ppg12_si, np.nan), where=current_si > 0)
    ratio_di = np.divide(ppg12_di, current_di, out=np.full_like(ppg12_di, np.nan), where=current_di > 0)
    ratio_ax.axhline(1.0, color="0.4", linewidth=0.9, linestyle=(0, (3, 3)))
    ratio_ax.plot(centers[valid_si], ratio_si[valid_si], "o", color="#1f4dff", markersize=2.8, label="SI")
    ratio_ax.plot(centers[valid_di], ratio_di[valid_di], "o", color="#ff4a4a", markersize=2.8, label="DI")
    ratio_values = np.concatenate((ratio_si[valid_si], ratio_di[valid_di]))
    ymax = max(1.5, float(np.max(ratio_values)) * 1.12) if ratio_values.size else 1.5
    ratio_ax.set_ylim(0.0, ymax)
    ratio_ax.set_ylabel("PPG12 /\nCurrent", fontsize=9.7)
    ratio_ax.set_xlabel("NPB Score", fontsize=11.5)
    ratio_ax.xaxis.set_major_locator(MultipleLocator(0.1))
    ratio_ax.xaxis.set_minor_locator(MultipleLocator(0.05))
    ratio_ax.tick_params(labelsize=9.5)
    ratio_ax.legend(loc="upper left", frameon=False, fontsize=8.3, ncol=2, handletextpad=0.3)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return ratio_si, ratio_di, np.asarray([valid_si, valid_di])


def main() -> int:
    args = parse_args()
    payload = remote_payload(args.remote_base, args.ppg12_results, args.with_ppg12_reference)
    edges_si, si, si_scales = weighted_component(payload, "si")
    edges_di, di, di_scales = weighted_component(payload, "di")
    if not np.allclose(edges_si, edges_di):
        raise RuntimeError("SI and DI NPB bin edges differ")
    below = edges_si[1:] <= 0.5 + 1e-12
    delta = float(np.sum(di[below]) - np.sum(si[below]))

    args.outdir.mkdir(parents=True, exist_ok=True)
    output = args.outdir / args.output_name
    if args.with_ppg12_reference and args.output_name == "ppg12_fig18_current_0mrad_npb_si_di_overlay.png":
        output = args.outdir / "ppg12_fig18_sdcc_vs_current_0mrad_npb_si_di_overlay_ratio.png"
    bins = args.outdir / (output.stem + "_bins.csv")
    manifest = args.outdir / (output.stem + ".manifest.json")
    if args.with_ppg12_reference:
        ppg_edges_si, ppg_si = ppg12_component(payload, "si")
        ppg_edges_di, ppg_di = ppg12_component(payload, "di")
        if not (np.allclose(edges_si, ppg_edges_si) and np.allclose(edges_si, ppg_edges_di)):
            raise RuntimeError("PPG12 and current NPB bin edges differ")
        ratio_si, ratio_di, valid_ratios = render_comparison(edges_si, si, di, ppg_si, ppg_di, output)
    else:
        ppg_si = ppg_di = ratio_si = ratio_di = valid_ratios = None
        render(edges_si, si, di, output, delta)
    with bins.open("w", newline="") as stream:
        fields = ["npb_lo", "npb_hi", "single_interaction", "double_interaction"]
        if args.with_ppg12_reference:
            fields += ["ppg12_single_interaction", "ppg12_double_interaction", "ppg12_over_current_si", "ppg12_over_current_di"]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for idx, (lo, hi, one, two) in enumerate(zip(edges_si[:-1], edges_si[1:], si, di)):
            row = {"npb_lo": lo, "npb_hi": hi, "single_interaction": one, "double_interaction": two}
            if args.with_ppg12_reference:
                row.update({"ppg12_single_interaction": ppg_si[idx], "ppg12_double_interaction": ppg_di[idx],
                            "ppg12_over_current_si": ratio_si[idx], "ppg12_over_current_di": ratio_di[idx]})
            writer.writerow(row)
    record = {
        "output_png": str(output),
        "bin_table": str(bins),
        "source_roots": payload["source_roots"],
        "source_histograms": "SIM/h1d_npb_score_eta0_pt{0..4}_cut0",
        "selection": "inclusive-jet MC, 0 mrad, eta0 (-0.7<eta<0.7), pT 10-30 GeV, TableQA cut0",
        "ppg12_reference_macro": "ppg12codeGit/plotting/plot_npb_score_di.C",
        "aggregation": "sum jet8/12/20/30/40 separately for SI and DI; scale each source by PPG12 xsec/jet50 divided by events_processed; unit-normalize each component",
        "ppg12_mix_note": "The PPG12 SI/DI mix factor is constant within each curve and cancels after separate unit normalization.",
        "vertex_weight_caveat": "Current TableQA spectra are raw candidate counts. Candidate-level truth-vertex weights in the historical ShowerShapeCheck source cannot be reconstructed after aggregation; this is a source-level visual diagnostic, not an absolute vertex-weighted parity claim.",
        "si_source_scales": si_scales,
        "di_source_scales": di_scales,
        "probability_below_npb_0p5": {"single": float(np.sum(si[below])), "double": float(np.sum(di[below])), "double_minus_single": delta},
    }
    if args.with_ppg12_reference:
        record.update({
            "ppg12_reference_roots": payload["ppg12_reference_roots"],
            "ppg12_reference_histograms": "h_npb_score_eta0_pt{0..4}_bdt{0..2}",
            "ratio_definition": "PPG12 SDCC unit-normalized component / current unit-normalized component; points shown only where both densities are >=1e-4",
            "ratio_valid_bins": {"single": int(np.count_nonzero(valid_ratios[0])), "double": int(np.count_nonzero(valid_ratios[1]))},
        })
    manifest.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(output)
    print(bins)
    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
