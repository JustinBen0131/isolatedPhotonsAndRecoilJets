#!/usr/bin/env python3
"""Regenerate THE-76 photon Fig.5 overlay from the fixed pre-vz campaign.

This script is intentionally separate from the older single-source diagnostic.
It reads the fresh Condor output for the source-stage pre-firstEventCuts/vz
histogram family and applies the PPG12 Fig.5 convention externally:

    current bin value = raw pre-vz count * XSEC[sample] / slimtree.num_entries

No RecoilJets metadata denominator, period/SI-DI duplication, global fitted
scale, or merge-stage artifact is used for the photon panel.
"""

from __future__ import annotations

import csv
import json
import math
import os
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
CAMPAIGN_FULL = "the76_ppg12_parity_full_20260701_003024"
CAMPAIGN_PREVZ = "the76_ppg12_fig5_prevz_strict_20260705_224243"
BASE_FULL = REPO / "dataOutput/ppg12Parity" / CAMPAIGN_FULL
BASE_PREVZ = REPO / "dataOutput/ppg12Parity" / CAMPAIGN_PREVZ
OUT_DIR = BASE_PREVZ / "strict_stitched_photonjet_prevz"
SOURCE_CSV = BASE_FULL / "reference_audit/sdcc_pull/live_ppg12_ian_source_points.csv"
JET_POINTS = BASE_FULL / "strict_stitched_inclusivejet/jet_source_scope_common_exposure_corrected_points.csv"
ROOT_MACRO = REPO / "macros/plotting/pp_currentian/MakeTHE76PhotonFig5OverlayRootStyle.C"
ROOT_EXE = Path("/Users/patsfan753/Desktop/analysis/env/bin/root")
POINTS_CSV = OUT_DIR / "photon_data_over_fit_sdcc_vs_current_overlay_prevz_strict_points.csv"
POINTS_MANIFEST = OUT_DIR / "photon_data_over_fit_sdcc_vs_current_overlay_prevz_strict_manifest.json"

REMOTE_BASE = f"/sphenix/tg/tg01/bulk/jbennett/thesisAna/sim/{CAMPAIGN_PREVZ}"
SETTING = "preselectionReference_tightReference_nonTightReference"
HIST = "SIM/h_ppPhotonStitch_ppg12Fig5PreVz_maxPhotonPt_kept"

SAMPLES = {
    "photon5": {
        "remote_dir": "run28_photonjet5",
        "xsec_pb": 146359.3,
        "lo": 0.0,
        "hi": 14.0,
        "ppg12_slimtree_entries": 9986593,
    },
    "photon10": {
        "remote_dir": "run28_photonjet10",
        "xsec_pb": 6944.675,
        "lo": 14.0,
        "hi": 22.0,
        "ppg12_slimtree_entries": 9998561,
    },
    "photon20": {
        "remote_dir": "run28_photonjet20",
        "xsec_pb": 130.4461,
        "lo": 22.0,
        "hi": 200.0,
        "ppg12_slimtree_entries": 9999988,
    },
}

COLORS = {
    "photon5": "#e83e8c",
    "photon10": "#2ca02c",
    "photon20": "#1da1f2",
    "jet8": "#e83e8c",
    "jet12": "#2ca02c",
    "jet20": "#1da1f2",
    "jet30": "#ff6f00",
    "jet40": "#d65ad1",
}


def _ratio_error(num: float, num_err: float, den: float, den_err: float) -> float:
    if not all(math.isfinite(v) for v in (num, num_err, den, den_err)) or num <= 0 or den <= 0:
        return math.nan
    return (num / den) * math.hypot(num_err / num, den_err / den)


def run_remote_aggregate() -> dict:
    remote_code = f"""
import json
from pathlib import Path
import ROOT

base = Path({REMOTE_BASE!r})
setting = {SETTING!r}
hist_name = {HIST!r}
samples = {json.dumps(SAMPLES)}
out = {{}}
for sample, spec in samples.items():
    d = base / setting / spec["remote_dir"]
    files = sorted(d.glob("*.root"))
    values = None
    errors2 = None
    edges = None
    missing = 0
    zombie = 0
    for p in files:
        f = ROOT.TFile.Open(str(p))
        if not f or f.IsZombie():
            zombie += 1
            continue
        h = f.Get(hist_name)
        if not h:
            missing += 1
            f.Close()
            continue
        if values is None:
            nb = h.GetNbinsX()
            values = [0.0] * nb
            errors2 = [0.0] * nb
            edges = [h.GetXaxis().GetBinLowEdge(i) for i in range(1, nb + 1)]
            edges.append(h.GetXaxis().GetBinUpEdge(nb))
        for i in range(1, h.GetNbinsX() + 1):
            values[i - 1] += h.GetBinContent(i)
            err = h.GetBinError(i)
            errors2[i - 1] += err * err
        f.Close()
    out[sample] = {{
        "source_dir": str(d),
        "files": len(files),
        "missing_histogram": missing,
        "zombie": zombie,
        "edges": edges,
        "values": values,
        "errors": [x ** 0.5 for x in errors2] if errors2 is not None else None,
        "histogram": hist_name,
    }}
print(json.dumps(out))
"""
    ssh_sock = os.environ.get("SSH_AUTH_SOCK") or subprocess.check_output(
        ["launchctl", "getenv", "SSH_AUTH_SOCK"], text=True
    ).strip()
    cmd = [
        "ssh",
        "patsfan753@ssh.sdcc.bnl.gov",
        "ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null sphnxuser05.sdcc.bnl.gov 'python3 -u -'",
    ]
    proc = subprocess.run(
        cmd,
        input=remote_code,
        text=True,
        capture_output=True,
        check=True,
        env={**os.environ, "SSH_AUTH_SOCK": ssh_sock},
    )
    lines = [line for line in proc.stdout.splitlines() if line.strip().startswith("{")]
    if not lines:
        raise RuntimeError(f"remote aggregate did not emit JSON; stdout={proc.stdout!r}; stderr={proc.stderr!r}")
    return json.loads(lines[-1])


def read_source_rows() -> dict[tuple[str, float], dict[str, float | str]]:
    out: dict[tuple[str, float], dict[str, float | str]] = {}
    with SOURCE_CSV.open() as f:
        for row in csv.DictReader(f):
            if row["group"] != "photon" or row["sample"] not in SAMPLES:
                continue
            x = float(row["bin_center"])
            if 10.0 <= x <= 40.0:
                out[(row["sample"], round(x, 6))] = {
                    "sample": row["sample"],
                    "bin_low": float(row["bin_low"]),
                    "bin_high": float(row["bin_high"]),
                    "bin_center": x,
                    "value": float(row["value"]),
                    "error": float(row["error"]),
                    "fit_value": float(row["fit_value"]),
                }
    return out


def write_points(aggregate: dict) -> tuple[Path, Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    source = read_source_rows()
    rows: list[dict[str, object]] = []
    for sample, spec in SAMPLES.items():
        rec = aggregate[sample]
        if rec["values"] is None:
            raise RuntimeError(f"No usable {HIST} objects found for {sample}: {rec}")
        edges = [float(x) for x in rec["edges"]]
        values = [float(x) for x in rec["values"]]
        errors = [float(x) for x in rec["errors"]]
        xsec = float(spec["xsec_pb"])
        denom = float(spec["ppg12_slimtree_entries"])
        for i, raw in enumerate(values):
            lo, hi = edges[i], edges[i + 1]
            x = 0.5 * (lo + hi)
            if not (10.0 <= x <= 40.0 and spec["lo"] <= x < spec["hi"]):
                continue
            srow = source.get((sample, round(x, 6)))
            if not srow:
                continue
            cur = raw * xsec / denom
            cur_err = errors[i] * xsec / denom
            fit = float(srow["fit_value"])
            sdcc = float(srow["value"])
            sdcc_err = float(srow["error"])
            rows.append(
                {
                    "sample": sample,
                    "bin_low": lo,
                    "bin_high": hi,
                    "bin_center": x,
                    "ppg12_sdcc_value": sdcc,
                    "ppg12_sdcc_error": sdcc_err,
                    "ppg12_fit_value": fit,
                    "ppg12_sdcc_over_fit": sdcc / fit,
                    "ppg12_sdcc_over_fit_error": sdcc_err / fit,
                    "current_value": cur,
                    "current_error": cur_err,
                    "current_over_fit": cur / fit,
                    "current_over_fit_error": cur_err / fit,
                    "current_over_sdcc": cur / sdcc,
                    "sdcc_over_current": sdcc / cur if cur > 0 else math.nan,
                    "current_value_mode": "prevz_raw_count_times_xsec_over_ppg12_slimtree_entries",
                    "current_source_files": int(rec["files"]),
                    "current_missing_histogram_files": int(rec["missing_histogram"]),
                    "current_zombie_files": int(rec["zombie"]),
                    "ppg12_slimtree_entries": int(denom),
                    "comparison_scope": "strict pre-vz photon5+10+20 stitch; no period/SI-DI aggregate; no global scale",
                }
            )
    fields = [
        "sample", "bin_low", "bin_high", "bin_center",
        "ppg12_sdcc_value", "ppg12_sdcc_error", "ppg12_fit_value",
        "ppg12_sdcc_over_fit", "ppg12_sdcc_over_fit_error",
        "current_value", "current_error", "current_over_fit",
        "current_over_fit_error", "current_over_sdcc", "sdcc_over_current",
        "current_value_mode", "current_source_files", "current_missing_histogram_files",
        "current_zombie_files", "ppg12_slimtree_entries", "comparison_scope",
    ]
    with POINTS_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    POINTS_MANIFEST.write_text(json.dumps({
        "status": "ok",
        "campaign_tag": CAMPAIGN_PREVZ,
        "points_csv": str(POINTS_CSV),
        "ppg12_source_csv": str(SOURCE_CSV),
        "remote_base": REMOTE_BASE,
        "setting": SETTING,
        "histogram": HIST,
        "samples": aggregate,
        "normalization": "current_value = raw pre-vz count * xsec / PPG12 slimtree.num_entries; no bin-width division and no global scale",
        "ppg12_slimtree_entries": {k: v["ppg12_slimtree_entries"] for k, v in SAMPLES.items()},
        "stitch_windows": {k: [v["lo"], v["hi"]] for k, v in SAMPLES.items()},
        "important_caveat": "This is the strict Fig.5 source-stage diagnostic. It intentionally bypasses RecoilJets post-vz metadata denominators and period/SI-DI aggregate physics weights.",
    }, indent=2) + "\n")
    return POINTS_CSV, POINTS_MANIFEST


def run_root_plot(points_csv: Path) -> tuple[Path, Path]:
    png = OUT_DIR / "photon_data_over_fit_sdcc_vs_current_overlay_root_ppg12_style_prevz_strict.png"
    manifest = OUT_DIR / "photon_data_over_fit_sdcc_vs_current_overlay_root_ppg12_style_prevz_strict_manifest.json"
    scope = (
        "Strict Fig.5 pre-vz diagnostic: PPG12 SDCC open markers vs RecoilJets "
        "pre-firstEventCuts source-stage photon5+10+20 stitch; XSEC/PPG12 slimTree entries; no global scale"
    )
    call = f'{ROOT_MACRO}("{points_csv}","{png}","{manifest}","{scope}")'
    subprocess.run(
        [str(REPO / "scripts/root_in_analysis_env.sh"), str(ROOT_EXE), "-b", "-l", "-q", call],
        check=True,
        cwd=REPO,
    )
    return png, manifest


def write_ratio_summary(points_csv: Path) -> tuple[Path, Path, Path]:
    out_png = OUT_DIR / "ppg12_over_current_stitched_spectra_two_panel_summary_prevz_strict_photon.png"
    out_csv = OUT_DIR / "ppg12_over_current_stitched_spectra_two_panel_summary_prevz_strict_photon_points.csv"
    out_manifest = OUT_DIR / "ppg12_over_current_stitched_spectra_two_panel_summary_prevz_strict_photon_manifest.json"
    rows: list[dict[str, object]] = []
    with points_csv.open() as f:
        for row in csv.DictReader(f):
            sdcc = float(row["ppg12_sdcc_value"])
            cur = float(row["current_value"])
            sdcc_err = float(row["ppg12_sdcc_error"])
            cur_err = float(row["current_error"])
            rows.append(
                {
                    "panel": "photon+jet",
                    "sample": row["sample"],
                    "bin_low": float(row["bin_low"]),
                    "bin_high": float(row["bin_high"]),
                    "bin_center": float(row["bin_center"]),
                    "ppg12_value": sdcc,
                    "ppg12_error": sdcc_err,
                    "current_value": cur,
                    "current_error": cur_err,
                    "ratio_ppg12_over_current": sdcc / cur,
                    "ratio_error": _ratio_error(sdcc, sdcc_err, cur, cur_err),
                    "current_definition": "strict pre-vz RecoilJets source-stage Fig.5, XSEC/PPG12 slimTree entries",
                }
            )
    if JET_POINTS.exists():
        with JET_POINTS.open() as f:
            for row in csv.DictReader(f):
                sdcc = float(row["ppg12_source_value"])
                cur = float(row["current_common_exposure_value"])
                sdcc_err = float(row["ppg12_source_error"])
                cur_err = float(row["current_common_exposure_error"])
                rows.append(
                    {
                        "panel": "inclusive jet",
                        "sample": row["sample"],
                        "bin_low": float(row["bin_low"]),
                        "bin_high": float(row["bin_high"]),
                        "bin_center": float(row["bin_center"]),
                        "ppg12_value": sdcc,
                        "ppg12_error": sdcc_err,
                        "current_value": cur,
                        "current_error": cur_err,
                        "ratio_ppg12_over_current": sdcc / cur,
                        "ratio_error": _ratio_error(sdcc, sdcc_err, cur, cur_err),
                        "current_definition": "existing inclusive-jet source-exposure corrected points; unchanged in this regeneration",
                    }
                )
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "mathtext.fontset": "dejavuserif",
        "axes.unicode_minus": False,
    })
    fig, axes = plt.subplots(2, 1, figsize=(16, 9), sharex=False)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.93, bottom=0.105, hspace=0.28)
    panels = [
        ("photon+jet", ["photon5", "photon10", "photon20"], (10, 40), (0.94, 1.06), "photon+jet stitched spectra, strict pre-vz diagnostic"),
        ("inclusive jet", ["jet8", "jet12", "jet20", "jet30", "jet40"], (9, 50), (0.86, 1.08), "inclusive jet stitched spectra, source-count corrected"),
    ]
    for ax, (panel, samples, xlim, ylim, title) in zip(axes, panels, strict=True):
        ax.axhline(1.0, color="0.45", lw=1.5, ls=(0, (5, 5)), zorder=0)
        for sample in samples:
            pts = [r for r in rows if r["panel"] == panel and r["sample"] == sample]
            if not pts:
                continue
            ax.errorbar(
                [float(r["bin_center"]) for r in pts],
                [float(r["ratio_ppg12_over_current"]) for r in pts],
                yerr=[float(r["ratio_error"]) for r in pts],
                fmt="o",
                ms=6.5,
                lw=1.2,
                elinewidth=1.1,
                capsize=0,
                color=COLORS[sample],
                markeredgecolor=COLORS[sample],
                markerfacecolor=COLORS[sample],
                label=sample,
            )
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_ylabel("PPG12 / Current", fontsize=24)
        ax.set_title(title, loc="left", fontsize=25, fontweight="bold", pad=8)
        ax.tick_params(axis="both", which="major", labelsize=20, direction="in", top=True, right=True, length=8, width=1.4)
        ax.tick_params(axis="both", which="minor", direction="in", top=True, right=True, length=4, width=1.0)
        ax.minorticks_on()
        ax.legend(loc="upper right", ncol=len(samples), frameon=False, fontsize=17, handletextpad=0.4, columnspacing=1.0, borderpad=0.1)
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)
    axes[0].tick_params(labelbottom=False)
    axes[1].set_xlabel("Leading object $p_T$ or $E_T$ [GeV]", fontsize=26)
    fig.savefig(out_png, dpi=160)
    plt.close(fig)
    ratios = [float(r["ratio_ppg12_over_current"]) for r in rows]
    out_manifest.write_text(json.dumps({
        "status": "ok",
        "plot_png": str(out_png),
        "points_csv": str(out_csv),
        "photon_points_csv": str(points_csv),
        "inclusive_jet_points_csv": str(JET_POINTS) if JET_POINTS.exists() else None,
        "ratio_definition": "PPG12 SDCC stitched spectrum divided by current analysis stitched spectrum",
        "photon_current_definition": "strict pre-vz RecoilJets source-stage Fig.5, XSEC/PPG12 slimTree entries",
        "inclusive_jet_current_definition": "unchanged existing inclusive-jet source-exposure corrected points",
        "ratio_min": min(ratios),
        "ratio_max": max(ratios),
    }, indent=2) + "\n")
    return out_png, out_csv, out_manifest


def main() -> None:
    aggregate = run_remote_aggregate()
    points_csv, points_manifest = write_points(aggregate)
    root_png, root_manifest = run_root_plot(points_csv)
    ratio_png, ratio_csv, ratio_manifest = write_ratio_summary(points_csv)
    print(points_csv)
    print(points_manifest)
    print(root_png)
    print(root_manifest)
    print(ratio_png)
    print(ratio_csv)
    print(ratio_manifest)


if __name__ == "__main__":
    main()
