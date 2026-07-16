#!/usr/bin/env python3
"""Overlay PPG12 Fig. 25 BDT profiles with the current combined output.

The two PPG12 profiles are extracted read-only from their SDCC ROOT files with
ROOT's own ``RebinX(2)`` and ``ProfileX`` implementation.  The current profile
is extracted from the locally registered current inclusive-SIM artifact using
the identical operations.  The lower panel shows each SDCC source divided by
the current profile with independent statistical errors propagated in
quadrature.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

try:
    import ROOT
except ModuleNotFoundError as exc:
    raise SystemExit(
        "PyROOT is required. Run with /Users/patsfan753/Desktop/analysis/env/bin/python3."
    ) from exc


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
CURRENT_POINTER = (
    REPO
    / "dataOutput/current_recoiljets_artifacts/current/pp_sim_inclusivejet_merged/current.json"
)
DEFAULT_OUTDIR = (
    REPO
    / "dataOutput/ppg12Parity/the97_ppg12_final_parity_full_20260709_2230/"
      "fig25_bdt_iso_profile"
)
HIST_CURRENT = "SIM/h2d_bdt_eta0_pt1_cut1"
HIST_SDCC = "h2d_bdt_eta0_pt1_cut1"
SDCC_BASE = "/sphenix/user/shuhangli/ppg12/efficiencytool/results"
SDCC_V1 = f"{SDCC_BASE}/MC_efficiencyshower_shape_jet_inclusive_combined_showershape.root"
SDCC_V3 = f"{SDCC_BASE}/MC_efficiencyshower_shape_jet_inclusive_showershape_etbin_v3E_v3E.root"
MARKER = "__PPG12_FIG25_PROFILE_PAYLOAD__"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--current-pointer", type=Path, default=CURRENT_POINTER)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--output-name",
        default="ppg12_fig25_bdt_sdcc_basev1_basev3e_vs_current_overlay_ratio.png",
    )
    return parser.parse_args()


def profile_from_th2(hist: Any, profile_name: str) -> dict[str, Any]:
    hist.SetDirectory(0)
    hist.RebinX(2)
    hist.GetXaxis().SetRangeUser(0.0, 1.0)
    profile = hist.ProfileX(profile_name, 1, -1, "")
    profile.SetDirectory(0)

    rows: list[dict[str, float]] = []
    axis = profile.GetXaxis()
    for index in range(1, profile.GetNbinsX() + 1):
        low = float(axis.GetBinLowEdge(index))
        high = float(axis.GetBinUpEdge(index))
        if high <= 0.0 or low >= 1.0:
            continue
        value = float(profile.GetBinContent(index))
        error = float(profile.GetBinError(index))
        entries = float(profile.GetBinEntries(index))
        if entries <= 0.0 or not np.isfinite(value) or not np.isfinite(error):
            continue
        rows.append(
            {
                "x_low": max(0.0, low),
                "x_high": min(1.0, high),
                "x": 0.5 * (max(0.0, low) + min(1.0, high)),
                "value": value,
                "error": error,
                "entries": entries,
            }
        )
    return {
        "correlation": float(hist.GetCorrelationFactor()),
        "rows": rows,
    }


def extract_local(root_path: Path) -> dict[str, Any]:
    ROOT.gROOT.SetBatch(True)
    source = ROOT.TFile.Open(str(root_path), "READ")
    if not source or source.IsZombie():
        raise RuntimeError(f"Cannot open current ROOT: {root_path}")
    obj = source.Get(HIST_CURRENT)
    if not obj:
        raise RuntimeError(f"Missing {HIST_CURRENT} in {root_path}")
    hist = obj.Clone("h2_current_fig25_extract")
    payload = profile_from_th2(hist, "pfx_current_fig25_extract")
    source.Close()
    payload.update(
        {
            "path": str(root_path),
            "histogram": HIST_CURRENT,
            "root_version": str(ROOT.gROOT.GetVersion()),
        }
    )
    return payload


def extract_sdcc() -> dict[str, Any]:
    remote_code = f'''import json, os
from datetime import datetime, timezone
import ROOT

ROOT.gROOT.SetBatch(True)
sources = {{"base_v1E": {SDCC_V1!r}, "base_v3E": {SDCC_V3!r}}}
payload = {{"root_version": str(ROOT.gROOT.GetVersion()), "profiles": {{}}}}
for label, path in sources.items():
    source = ROOT.TFile.Open(path, "READ")
    if not source or source.IsZombie():
        raise RuntimeError(f"Cannot open {{path}}")
    obj = source.Get({HIST_SDCC!r})
    if not obj:
        raise RuntimeError(f"Missing {HIST_SDCC} in {{path}}")
    hist = obj.Clone(f"h2_{{label}}_fig25_extract")
    hist.SetDirectory(0)
    hist.RebinX(2)
    hist.GetXaxis().SetRangeUser(0.0, 1.0)
    profile = hist.ProfileX(f"pfx_{{label}}_fig25_extract", 1, -1, "")
    profile.SetDirectory(0)
    rows = []
    axis = profile.GetXaxis()
    for index in range(1, profile.GetNbinsX() + 1):
        low = float(axis.GetBinLowEdge(index))
        high = float(axis.GetBinUpEdge(index))
        if high <= 0.0 or low >= 1.0:
            continue
        value = float(profile.GetBinContent(index))
        error = float(profile.GetBinError(index))
        entries = float(profile.GetBinEntries(index))
        if entries <= 0.0:
            continue
        rows.append({{
            "x_low": max(0.0, low), "x_high": min(1.0, high),
            "x": 0.5 * (max(0.0, low) + min(1.0, high)),
            "value": value, "error": error, "entries": entries,
        }})
    stat = os.stat(path)
    payload["profiles"][label] = {{
        "path": path,
        "histogram": {HIST_SDCC!r},
        "size_bytes": int(stat.st_size),
        "mtime_iso": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
        "correlation": float(hist.GetCorrelationFactor()),
        "rows": rows,
    }}
    source.Close()
print({MARKER!r})
print(json.dumps(payload, separators=(",", ":")))
'''
    sock = os.environ.get("SSH_AUTH_SOCK", "")
    if not sock:
        found = subprocess.run(
            ["launchctl", "getenv", "SSH_AUTH_SOCK"],
            text=True,
            capture_output=True,
            check=False,
        )
        sock = found.stdout.strip()
    if not sock:
        raise RuntimeError("SSH_AUTH_SOCK is unavailable")

    inner = (
        "ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null "
        "-o BatchMode=yes sphnxuser05.sdcc.bnl.gov "
        + shlex.quote("python3 - <<'PY'\n" + remote_code + "\nPY")
    )
    env = os.environ.copy()
    env["SSH_AUTH_SOCK"] = sock
    result = subprocess.run(
        [
            "ssh",
            "-o",
            "BatchMode=yes",
            "-o",
            "ConnectTimeout=20",
            "patsfan753@ssh.sdcc.bnl.gov",
            inner,
        ],
        text=True,
        capture_output=True,
        check=False,
        env=env,
    )
    if result.returncode != 0:
        raise RuntimeError(f"Read-only SDCC extraction failed:\n{result.stderr[-3000:]}")
    if MARKER not in result.stdout:
        raise RuntimeError("SDCC extraction returned no marked JSON payload")
    return json.loads(result.stdout.rsplit(MARKER, 1)[1].strip())


def arrays(profile: dict[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rows = profile["rows"]
    return (
        np.asarray([row["x"] for row in rows], dtype=float),
        np.asarray([row["value"] for row in rows], dtype=float),
        np.asarray([row["error"] for row in rows], dtype=float),
        np.asarray([row["x_high"] - row["x_low"] for row in rows], dtype=float),
    )


def aligned_ratio(
    numerator: dict[str, Any], denominator: dict[str, Any]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xn, yn, en, _ = arrays(numerator)
    xd, yd, ed, _ = arrays(denominator)
    if len(xn) != len(xd) or not np.allclose(xn, xd, rtol=0.0, atol=1e-10):
        raise RuntimeError("Profile bin centers are not aligned")
    valid = np.isfinite(yn) & np.isfinite(en) & np.isfinite(yd) & np.isfinite(ed) & (yd != 0.0)
    ratio = np.full_like(yn, np.nan)
    error = np.full_like(yn, np.nan)
    ratio[valid] = yn[valid] / yd[valid]
    nonzero_num = valid & (yn != 0.0)
    error[nonzero_num] = np.abs(ratio[nonzero_num]) * np.sqrt(
        (en[nonzero_num] / yn[nonzero_num]) ** 2
        + (ed[nonzero_num] / yd[nonzero_num]) ** 2
    )
    zero_num = valid & (yn == 0.0)
    error[zero_num] = en[zero_num] / np.abs(yd[zero_num])
    return xn, ratio, error


def render(
    historical: dict[str, Any],
    v3e: dict[str, Any],
    current: dict[str, Any],
    output: Path,
) -> dict[str, float]:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.linewidth": 1.1,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig = plt.figure(figsize=(6.5, 7.0), dpi=180)
    grid = fig.add_gridspec(2, 1, height_ratios=(3.15, 1.15), hspace=0.035)
    top = fig.add_subplot(grid[0])
    ratio_ax = fig.add_subplot(grid[1], sharex=top)

    series = [
        (historical, "PPG12 SDCC base_v1E (Fig. 25)", "black", "o", "none"),
        (v3e, "PPG12 SDCC base_v3E", "#d62728", "s", "none"),
        (current, "Current combined base_v3E", "#1559d6", "o", "#1559d6"),
    ]
    upper_extent = 0.0
    for profile, label, color, marker, fill in series:
        x, y, error, _ = arrays(profile)
        upper_extent = max(upper_extent, float(np.nanmax(y + error)))
        top.errorbar(
            x,
            y,
            yerr=error,
            color=color,
            marker=marker,
            markerfacecolor=fill,
            markeredgecolor=color,
            markersize=3.1,
            linewidth=1.05,
            elinewidth=0.75,
            capsize=0.0,
            drawstyle="steps-mid",
            label=label,
            zorder=3,
        )

    top.set_xlim(0.0, 1.0)
    top.set_ylim(0.0, max(8.0, 1.08 * upper_extent))
    top.set_ylabel(r"$\langle E_T^{\mathrm{iso}}\rangle$ [GeV]", fontsize=12)
    top.tick_params(labelbottom=False, labelsize=10)
    top.yaxis.set_minor_locator(AutoMinorLocator(5))
    top.text(
        0.025,
        0.965,
        r"$\bf{\it{sPHENIX}}$ Internal",
        transform=top.transAxes,
        ha="left",
        va="top",
        fontsize=11.5,
    )
    top.text(
        0.025,
        0.905,
        r"$p+p\ \sqrt{s}=200$ GeV" + "\n" + r"$|\eta|<0.7$",
        transform=top.transAxes,
        ha="left",
        va="top",
        fontsize=10.3,
    )
    top.text(
        0.975,
        0.715,
        r"$14<p_T<18$ GeV, w/ nbkg cut" + "\nBackground MC",
        transform=top.transAxes,
        ha="right",
        va="top",
        fontsize=9.5,
    )
    top.legend(loc="upper right", frameon=False, fontsize=9.0, handlelength=2.1)

    ratios = []
    ratio_specs = [
        (historical, "base_v1E / Current", "black", "o"),
        (v3e, "base_v3E / Current", "#d62728", "s"),
    ]
    for profile, label, color, marker in ratio_specs:
        x, ratio, error = aligned_ratio(profile, current)
        ratios.append((x, ratio, error))
        ratio_ax.errorbar(
            x,
            ratio,
            yerr=error,
            color=color,
            marker=marker,
            markerfacecolor="none",
            markeredgecolor=color,
            markersize=3.0,
            linewidth=0.9,
            elinewidth=0.7,
            capsize=0.0,
            drawstyle="steps-mid",
            label=label,
        )

    finite_lows: list[float] = []
    finite_highs: list[float] = []
    for _, ratio, error in ratios:
        valid = np.isfinite(ratio) & np.isfinite(error)
        finite_lows.extend((ratio[valid] - error[valid]).tolist())
        finite_highs.extend((ratio[valid] + error[valid]).tolist())
    low = min(finite_lows + [1.0])
    high = max(finite_highs + [1.0])
    span = max(high - low, 0.2)
    ratio_ax.set_ylim(max(0.0, low - 0.08 * span), high + 0.08 * span)
    ratio_ax.axhline(1.0, color="0.35", linewidth=0.9, linestyle=(0, (4, 3)))
    ratio_ax.set_xlabel("BDT score", fontsize=11.5)
    ratio_ax.set_ylabel("SDCC / Current", fontsize=10.5)
    ratio_ax.xaxis.set_major_locator(MultipleLocator(0.1))
    ratio_ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ratio_ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ratio_ax.tick_params(labelsize=9.5)
    ratio_ax.legend(loc="lower left", frameon=False, fontsize=8.5, handlelength=2.0)

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.subplots_adjust(left=0.16, right=0.975, bottom=0.105, top=0.985)
    fig.savefig(output)
    plt.close(fig)
    return {
        "ratio_y_min": float(ratio_ax.get_ylim()[0]),
        "ratio_y_max": float(ratio_ax.get_ylim()[1]),
        "upper_y_max": max(8.0, 1.08 * upper_extent),
    }


def main() -> int:
    args = parse_args()
    pointer = json.loads(args.current_pointer.read_text())
    root_path = Path(pointer["root_paths"][0])
    current = extract_local(root_path)
    sdcc = extract_sdcc()
    historical = sdcc["profiles"]["base_v1E"]
    v3e = sdcc["profiles"]["base_v3E"]

    args.outdir.mkdir(parents=True, exist_ok=True)
    output = args.outdir / args.output_name
    csv_path = args.outdir / (output.stem + ".csv")
    payload_path = args.outdir / (output.stem + ".source_profiles.json")
    manifest_path = args.outdir / (output.stem + ".manifest.json")

    bounds = render(historical, v3e, current, output)
    x_h, r_h, e_h = aligned_ratio(historical, current)
    x_v, r_v, e_v = aligned_ratio(v3e, current)
    x_c, y_c, e_c, _ = arrays(current)
    _, y_h, ey_h, _ = arrays(historical)
    _, y_v, ey_v, _ = arrays(v3e)
    with csv_path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "bdt_center",
                "ppg12_base_v1E_mean_iso",
                "ppg12_base_v1E_error",
                "ppg12_base_v3E_mean_iso",
                "ppg12_base_v3E_error",
                "current_base_v3E_mean_iso",
                "current_base_v3E_error",
                "ppg12_base_v1E_over_current",
                "ppg12_base_v1E_over_current_error",
                "ppg12_base_v3E_over_current",
                "ppg12_base_v3E_over_current_error",
            ]
        )
        for index in range(len(x_c)):
            writer.writerow(
                [
                    f"{x_c[index]:.6g}",
                    f"{y_h[index]:.9g}",
                    f"{ey_h[index]:.9g}",
                    f"{y_v[index]:.9g}",
                    f"{ey_v[index]:.9g}",
                    f"{y_c[index]:.9g}",
                    f"{e_c[index]:.9g}",
                    f"{r_h[index]:.9g}",
                    f"{e_h[index]:.9g}",
                    f"{r_v[index]:.9g}",
                    f"{e_v[index]:.9g}",
                ]
            )

    payload = {
        "sdcc": sdcc,
        "current": current,
    }
    payload_path.write_text(json.dumps(payload, indent=2) + "\n")
    sha = ""
    basis = pointer.get("promotion_basis", "")
    if "sha256 " in basis:
        sha = basis.split("sha256 ", 1)[1].split()[0]
    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "figure": "PPG12 Fig. 25 BDT-isolation profile: SDCC base_v1E and base_v3E versus current",
        "current_pointer": str(args.current_pointer),
        "current_root": str(root_path),
        "current_sha256": sha,
        "sdcc_extraction": "read-only nested SSH; compact ROOT ProfileX arrays only",
        "histograms": {"sdcc": HIST_SDCC, "current": HIST_CURRENT},
        "operation": "RebinX(2), x range 0-1, ProfileX over all isolation bins",
        "ratio_definition": "SDCC profile / current profile; independent statistical errors propagated in quadrature",
        "correlations": {
            "ppg12_base_v1E": historical["correlation"],
            "ppg12_base_v3E": v3e["correlation"],
            "current_base_v3E": current["correlation"],
        },
        "axis_bounds": bounds,
        "output_png": str(output),
        "source_profiles_json": str(payload_path),
        "bin_csv": str(csv_path),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
