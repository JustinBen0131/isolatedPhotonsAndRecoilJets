#!/usr/bin/env python3
"""Compare pp-data non-tight Eiso sideband shapes against PPG12.

This diagnostic is deliberately read-only.  It compares the RecoilJets data
ABCD isolation distributions with Shuhang's PPG12 `data_histo_bdt_nom.root`
histograms to isolate whether the remaining purity mismatch is driven by the
data-side non-tight isolation split.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import textwrap
from pathlib import Path
from typing import Any

import ROOT


ROOT.gROOT.SetBatch(True)

TRIGGER_DIR = "Photon_4_GeV_plus_MBD_NS_geq_1"
PT_BINS = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36]
DEFAULT_CURRENT_ROOT = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/InputFiles/pp24/"
    "ppg12_photon_yield_v1_data_20260620/pp/"
    "RecoilJets_pp_ALL_jetMinPtScan_dphiScan_vz60_isoR40_isSliding_"
    "preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12.root"
)
DEFAULT_REMOTE_PPG12_ROOT = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo_bdt_nom.root"
DEFAULT_OUTDIR = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/ppg12PhotonYield/"
    "ppg12_photon_yield_v1_data_20260620/purity_fig29_comparison/"
    "globalmbd_mbddigi_componentmix_20260629_current_pp"
)


def root_open(path: str | Path) -> ROOT.TFile:
    f = ROOT.TFile.Open(str(path))
    if not f or f.IsZombie():
        raise RuntimeError(f"could not open ROOT file: {path}")
    return f


def hist_stats(h: ROOT.TH1, iso_thr: float, noniso_thr: float) -> dict[str, float]:
    nb = h.GetNbinsX()
    total_all = 0.0
    total_regular = 0.0
    underflow = float(h.GetBinContent(0))
    overflow = float(h.GetBinContent(nb + 1))
    total_all += underflow + overflow
    mean_num = 0.0
    iso = 0.0
    gap = 0.0
    noniso_window = 0.0
    over20 = 0.0
    below_min = 0.0
    for ib in range(1, nb + 1):
        y = float(h.GetBinContent(ib))
        x = float(h.GetXaxis().GetBinCenter(ib))
        total_regular += y
        total_all += y
        mean_num += x * y
        if x <= -20.0:
            below_min += y
        elif x < iso_thr:
            iso += y
        elif x <= noniso_thr:
            gap += y
        elif x < 20.0:
            noniso_window += y
        else:
            over20 += y
    mean = mean_num / total_regular if total_regular else math.nan
    return {
        "total_all": total_all,
        "total_regular": total_regular,
        "underflow": underflow,
        "overflow": overflow,
        "mean_regular": mean,
        "below_min_center": below_min,
        "iso_center": iso,
        "gap_center": gap,
        "noniso_window_center": noniso_window,
        "over20_center": over20,
        "frac_iso_center": iso / total_regular if total_regular else math.nan,
        "frac_gap_center": gap / total_regular if total_regular else math.nan,
        "frac_noniso_window_center": noniso_window / total_regular if total_regular else math.nan,
        "frac_overflow_all": overflow / total_all if total_all else math.nan,
    }


def compact_values(obj: Any, names: list[str]) -> dict[str, list[float]]:
    out: dict[str, list[float]] = {}
    for name in names:
        h = obj.Get(name)
        if not h:
            raise RuntimeError(f"missing histogram {name}")
        out[name] = [float(h.GetBinContent(i)) for i in range(1, h.GetNbinsX() + 1)]
    return out


def current_payload(path: Path, directory: str) -> dict[str, Any]:
    f = root_open(path)
    d = f.Get(directory)
    if not d:
        raise RuntimeError(f"missing directory {directory} in {path}")
    compact = compact_values(
        d,
        ["h_nontight_iso_cluster_0", "h_nontight_noniso_cluster_0"],
    )
    rows: list[dict[str, Any]] = []
    for ipt, (lo, hi) in enumerate(zip(PT_BINS[:-1], PT_BINS[1:])):
        mid = 0.5 * (lo + hi)
        iso_thr = 0.490 + 0.037 * mid
        noniso_thr = iso_thr + 0.8
        hname = f"h_Eiso_nonTight_pT_{lo}_{hi}"
        h = d.Get(hname)
        if not h:
            raise RuntimeError(f"missing current histogram {directory}/{hname}")
        c_name = f"h_Eiso_ABCD_C_pT_{lo}_{hi}"
        d_name = f"h_Eiso_ABCD_D_pT_{lo}_{hi}"
        h_c = d.Get(c_name)
        h_d = d.Get(d_name)
        if not h_c or not h_d:
            raise RuntimeError(f"missing current ABCD Eiso histogram for {lo}-{hi}")
        c_int = sum(float(h_c.GetBinContent(i)) for i in range(0, h_c.GetNbinsX() + 2))
        d_int = sum(float(h_d.GetBinContent(i)) for i in range(0, h_d.GetNbinsX() + 2))
        row = {
            "bin": f"{lo:g}-{hi:g}",
            "pt_lo": lo,
            "pt_hi": hi,
            "iso_thr_mid": iso_thr,
            "noniso_thr_mid": noniso_thr,
            **hist_stats(h, iso_thr, noniso_thr),
            "compact_C": compact["h_nontight_iso_cluster_0"][ipt],
            "compact_D": compact["h_nontight_noniso_cluster_0"][ipt],
            "eiso_C_integral_all": c_int,
            "eiso_D_integral_all": d_int,
        }
        rows.append(row)
    return {
        "path": str(path),
        "directory": directory,
        "rows": rows,
    }


def remote_ppg12_payload(remote_root: str, login_host: str, worker_host: str) -> dict[str, Any]:
    script = textwrap.dedent(
        f"""
        import json, math
        import ROOT
        ROOT.gROOT.SetBatch(True)
        pt_bins = {PT_BINS!r}
        path = {remote_root!r}
        f = ROOT.TFile.Open(path)
        if not f or f.IsZombie():
            raise SystemExit("could not open " + path)

        def hist_stats(h, iso_thr, noniso_thr):
            nb = h.GetNbinsX()
            total_all = float(h.GetBinContent(0) + h.GetBinContent(nb + 1))
            total_regular = 0.0
            mean_num = 0.0
            underflow = float(h.GetBinContent(0))
            overflow = float(h.GetBinContent(nb + 1))
            iso = gap = noniso_window = over20 = below_min = 0.0
            for ib in range(1, nb + 1):
                y = float(h.GetBinContent(ib))
                x = float(h.GetXaxis().GetBinCenter(ib))
                total_regular += y
                total_all += y
                mean_num += x * y
                if x <= -20.0:
                    below_min += y
                elif x < iso_thr:
                    iso += y
                elif x <= noniso_thr:
                    gap += y
                elif x < 20.0:
                    noniso_window += y
                else:
                    over20 += y
            return {{
                "total_all": total_all,
                "total_regular": total_regular,
                "underflow": underflow,
                "overflow": overflow,
                "mean_regular": mean_num / total_regular if total_regular else float("nan"),
                "below_min_center": below_min,
                "iso_center": iso,
                "gap_center": gap,
                "noniso_window_center": noniso_window,
                "over20_center": over20,
                "frac_iso_center": iso / total_regular if total_regular else float("nan"),
                "frac_gap_center": gap / total_regular if total_regular else float("nan"),
                "frac_noniso_window_center": noniso_window / total_regular if total_regular else float("nan"),
                "frac_overflow_all": overflow / total_all if total_all else float("nan"),
            }}

        compact = {{}}
        for name in ["h_nontight_iso_cluster_0", "h_nontight_noniso_cluster_0"]:
            h = f.Get(name)
            if not h:
                raise SystemExit("missing " + name)
            compact[name] = [float(h.GetBinContent(i)) for i in range(1, h.GetNbinsX() + 1)]

        rows = []
        for ipt, (lo, hi) in enumerate(zip(pt_bins[:-1], pt_bins[1:])):
            mid = 0.5 * (lo + hi)
            iso_thr = 0.490 + 0.037 * mid
            noniso_thr = iso_thr + 0.8
            hname = f"h_nontight_isoET_0_{{ipt}}"
            h = f.Get(hname)
            if not h:
                raise SystemExit("missing " + hname)
            row = {{
                "bin": f"{{lo:g}}-{{hi:g}}",
                "pt_lo": lo,
                "pt_hi": hi,
                "iso_thr_mid": iso_thr,
                "noniso_thr_mid": noniso_thr,
                **hist_stats(h, iso_thr, noniso_thr),
                "compact_C": compact["h_nontight_iso_cluster_0"][ipt],
                "compact_D": compact["h_nontight_noniso_cluster_0"][ipt],
            }}
            rows.append(row)
        print("JSON_PAYLOAD_BEGIN")
        print(json.dumps({{"path": path, "rows": rows}}))
        print("JSON_PAYLOAD_END")
        """
    ).strip()
    env = os.environ.copy()
    sock = subprocess.run(["launchctl", "getenv", "SSH_AUTH_SOCK"], text=True, capture_output=True, check=False)
    if sock.stdout.strip():
        env["SSH_AUTH_SOCK"] = sock.stdout.strip()
    cmd = [
        "ssh",
        login_host,
        f"ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null {worker_host} 'python3 -'",
    ]
    result = subprocess.run(cmd, input=script, text=True, capture_output=True, env=env, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"remote PPG12 read failed rc={result.returncode}\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}")
    text = result.stdout
    begin = text.index("JSON_PAYLOAD_BEGIN") + len("JSON_PAYLOAD_BEGIN")
    end = text.index("JSON_PAYLOAD_END")
    return json.loads(text[begin:end].strip())


def div(num: float, den: float) -> float:
    return num / den if den else math.nan


def fmt(v: float) -> str:
    return "nan" if not math.isfinite(v) else f"{v:.6g}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--current-root", type=Path, default=DEFAULT_CURRENT_ROOT)
    ap.add_argument("--current-dir", default=TRIGGER_DIR)
    ap.add_argument("--remote-ppg12-root", default=DEFAULT_REMOTE_PPG12_ROOT)
    ap.add_argument("--login-host", default="patsfan753@ssh.sdcc.bnl.gov")
    ap.add_argument("--worker-host", default="sphnxuser05.sdcc.bnl.gov")
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    ap.add_argument("--tag", default="20260629")
    args = ap.parse_args()

    current = current_payload(args.current_root, args.current_dir)
    ppg12 = remote_ppg12_payload(args.remote_ppg12_root, args.login_host, args.worker_host)

    args.outdir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    for c, p in zip(current["rows"], ppg12["rows"]):
        if c["bin"] != p["bin"]:
            raise RuntimeError(f"bin mismatch current={c['bin']} ppg12={p['bin']}")
        row = {
            "bin": c["bin"],
            "pt_lo": c["pt_lo"],
            "pt_hi": c["pt_hi"],
            "current_total_regular": c["total_regular"],
            "ppg12_total_regular": p["total_regular"],
            "current_mean_eiso": c["mean_regular"],
            "ppg12_mean_eiso": p["mean_regular"],
            "current_frac_iso_center": c["frac_iso_center"],
            "ppg12_frac_iso_center": p["frac_iso_center"],
            "current_frac_gap_center": c["frac_gap_center"],
            "ppg12_frac_gap_center": p["frac_gap_center"],
            "current_frac_noniso_window_center": c["frac_noniso_window_center"],
            "ppg12_frac_noniso_window_center": p["frac_noniso_window_center"],
            "current_C": c["compact_C"],
            "ppg12_C": p["compact_C"],
            "current_D": c["compact_D"],
            "ppg12_D": p["compact_D"],
            "current_D_over_CplusD": div(c["compact_D"], c["compact_C"] + c["compact_D"]),
            "ppg12_D_over_CplusD": div(p["compact_D"], p["compact_C"] + p["compact_D"]),
            "current_Eiso_D_integral_all": c["eiso_D_integral_all"],
            "current_Eiso_C_integral_all": c["eiso_C_integral_all"],
        }
        row["ratio_D_over_CplusD"] = div(row["current_D_over_CplusD"], row["ppg12_D_over_CplusD"])
        row["ratio_noniso_window_frac"] = div(
            row["current_frac_noniso_window_center"],
            row["ppg12_frac_noniso_window_center"],
        )
        rows.append(row)

    payload = {"current": current, "ppg12": ppg12, "rows": rows}
    json_path = args.outdir / f"the76_data_side_eiso_sideband_current_vs_ppg12_{args.tag}.json"
    csv_path = args.outdir / f"the76_data_side_eiso_sideband_current_vs_ppg12_{args.tag}.csv"
    md_path = args.outdir / f"the76_data_side_eiso_sideband_current_vs_ppg12_{args.tag}.md"
    json_path.write_text(json.dumps(payload, indent=2))

    import csv

    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "# THE-76 pp-data non-tight Eiso sideband parity",
        "",
        f"Current ROOT: `{args.current_root}` / `{args.current_dir}`",
        f"PPG12 ROOT: `{args.remote_ppg12_root}`",
        "",
        "## What this checks",
        "",
        "This compares non-tight Eiso spectra and exact compact C/D counts. The window fractions use the PPG12 mid-bin thresholds only as a shape diagnostic; the compact C/D columns are the authoritative region counts from each file.",
        "",
        "## Non-tight sideband readout",
        "",
        "| bin | current mean Eiso | PPG12 mean Eiso | current noniso frac | PPG12 noniso frac | noniso-frac ratio | current D/(C+D) | PPG12 D/(C+D) | D-split ratio |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            f"| {row['bin']} | {fmt(row['current_mean_eiso'])} | {fmt(row['ppg12_mean_eiso'])} | "
            f"{fmt(row['current_frac_noniso_window_center'])} | {fmt(row['ppg12_frac_noniso_window_center'])} | "
            f"{fmt(row['ratio_noniso_window_frac'])} | {fmt(row['current_D_over_CplusD'])} | "
            f"{fmt(row['ppg12_D_over_CplusD'])} | {fmt(row['ratio_D_over_CplusD'])} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- If the current non-tight Eiso spectrum is shifted lower than PPG12 and D/(C+D) is suppressed, the remaining purity mismatch is data-side isolation construction/window parity, not the final purity solver.",
            "- PPG12 `config_bdt_nom.yaml` uses `use_topo_iso: 2`, i.e. `cluster_iso_topo_04`, with `-20 < Eiso < 0.490+0.037*ET` for iso and `0.490+0.037*ET+0.8 < Eiso < 20` for noniso.",
            "- Current `fillIsoSSTagCounters` is still fed by the generic PhotonClusterBuilder `eiso()` path, while the PPG12 signal-yield path has a separate topo-isolation helper.",
            "",
        ]
    )
    md_path.write_text("\n".join(lines))

    print(f"wrote {json_path}")
    print(f"wrote {csv_path}")
    print(f"wrote {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
