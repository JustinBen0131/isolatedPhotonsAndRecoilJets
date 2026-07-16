#!/usr/bin/env python3
"""Compare PPG12 and RecoilJets pp-data stage histograms for Fig. 29 parity.

The goal is to identify where the data-side ABCD populations first diverge
after the signal-leakage comparison has been made compatible with PPG12.
The PPG12 ROOT can be read directly if local, or read-only over the standard
SDCC nested SSH path.
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

TRIGGER_DIR = "PPG12_scaledtrigger30"
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

HISTS = [
    "h_all_cluster_0",
    "h_common_cluster_0",
    "h_tight_cluster_0",
    "h_tight_iso_cluster_0",
    "h_tight_noniso_cluster_0",
    "h_nontight_iso_cluster_0",
    "h_nontight_noniso_cluster_0",
]


def root_open(path: str | Path) -> ROOT.TFile:
    f = ROOT.TFile.Open(str(path))
    if not f or f.IsZombie():
        raise RuntimeError(f"could not open ROOT file: {path}")
    return f


def hist_payload(path: str | Path, directory: str | None) -> dict[str, Any]:
    f = root_open(path)
    d = f.Get(directory) if directory else f
    if not d:
        raise RuntimeError(f"missing directory {directory} in {path}")
    payload: dict[str, Any] = {"path": str(path), "directory": directory or "", "hists": {}}
    for name in HISTS:
        h = d.Get(name)
        if not h:
            payload["hists"][name] = {"present": False, "values": [], "errors": [], "edges": []}
            continue
        values = [float(h.GetBinContent(i)) for i in range(1, h.GetNbinsX() + 1)]
        errors = [float(h.GetBinError(i)) for i in range(1, h.GetNbinsX() + 1)]
        edges = [float(h.GetXaxis().GetBinLowEdge(i)) for i in range(1, h.GetNbinsX() + 1)]
        edges.append(float(h.GetXaxis().GetBinUpEdge(h.GetNbinsX())))
        payload["hists"][name] = {"present": True, "values": values, "errors": errors, "edges": edges}
    return payload


def remote_ppg12_payload(remote_root: str, login_host: str, worker_host: str) -> dict[str, Any]:
    script = textwrap.dedent(
        f"""
        import json
        import ROOT
        ROOT.gROOT.SetBatch(True)
        hists = {HISTS!r}
        path = {remote_root!r}
        f = ROOT.TFile.Open(path)
        if not f or f.IsZombie():
            raise SystemExit("could not open " + path)
        payload = {{"path": path, "directory": "", "hists": {{}}}}
        for name in hists:
            h = f.Get(name)
            if not h:
                payload["hists"][name] = {{"present": False, "values": [], "errors": [], "edges": []}}
                continue
            values = [float(h.GetBinContent(i)) for i in range(1, h.GetNbinsX() + 1)]
            errors = [float(h.GetBinError(i)) for i in range(1, h.GetNbinsX() + 1)]
            edges = [float(h.GetXaxis().GetBinLowEdge(i)) for i in range(1, h.GetNbinsX() + 1)]
            edges.append(float(h.GetXaxis().GetBinUpEdge(h.GetNbinsX())))
            payload["hists"][name] = {{"present": True, "values": values, "errors": errors, "edges": edges}}
        print("JSON_PAYLOAD_BEGIN")
        print(json.dumps(payload))
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
    if den == 0.0:
        return math.nan
    return num / den


def fmt(value: float) -> str:
    if not math.isfinite(value):
        return "nan"
    return f"{value:.6g}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--current-root", type=Path, default=DEFAULT_CURRENT_ROOT)
    ap.add_argument("--current-dir", default=TRIGGER_DIR)
    ap.add_argument("--ppg12-root", type=Path)
    ap.add_argument("--remote-ppg12-root", default=DEFAULT_REMOTE_PPG12_ROOT)
    ap.add_argument("--login-host", default="patsfan753@ssh.sdcc.bnl.gov")
    ap.add_argument("--worker-host", default="sphnxuser05.sdcc.bnl.gov")
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    ap.add_argument("--tag", default="20260629")
    args = ap.parse_args()

    current = hist_payload(args.current_root, args.current_dir)
    if args.ppg12_root:
        ppg12 = hist_payload(args.ppg12_root, None)
        ppg12_label = str(args.ppg12_root)
    else:
        ppg12 = remote_ppg12_payload(args.remote_ppg12_root, args.login_host, args.worker_host)
        ppg12_label = args.remote_ppg12_root

    args.outdir.mkdir(parents=True, exist_ok=True)
    json_path = args.outdir / f"the76_data_side_stage_hists_current_vs_ppg12_{args.tag}.json"
    json_path.write_text(json.dumps({"current": current, "ppg12": ppg12}, indent=2))

    edges = current["hists"]["h_all_cluster_0"]["edges"]
    ppg_edges = ppg12["hists"]["h_all_cluster_0"]["edges"]
    if len(edges) != len(ppg_edges) or any(abs(a - b) > 1.0e-6 for a, b in zip(edges, ppg_edges)):
        raise RuntimeError(f"bin edge mismatch: current={edges}, ppg12={ppg_edges}")

    rows: list[dict[str, float | str]] = []
    for i in range(len(edges) - 1):
        c = {name: current["hists"][name]["values"][i] for name in HISTS}
        p = {name: ppg12["hists"][name]["values"][i] for name in HISTS}
        row: dict[str, float | str] = {
            "pt_lo": edges[i],
            "pt_hi": edges[i + 1],
            "bin": f"{edges[i]:g}-{edges[i + 1]:g}",
        }
        for name in HISTS:
            short = name.removeprefix("h_").removesuffix("_0")
            row[f"current_{short}"] = c[name]
            row[f"ppg12_{short}"] = p[name]
            row[f"ratio_{short}"] = div(c[name], p[name])
        row.update(
            {
                "current_common_over_all": div(c["h_common_cluster_0"], c["h_all_cluster_0"]),
                "ppg12_common_over_all": div(p["h_common_cluster_0"], p["h_all_cluster_0"]),
                "current_tight_over_common": div(c["h_tight_cluster_0"], c["h_common_cluster_0"]),
                "ppg12_tight_over_common": div(p["h_tight_cluster_0"], p["h_common_cluster_0"]),
                "current_A_over_tight": div(c["h_tight_iso_cluster_0"], c["h_tight_cluster_0"]),
                "ppg12_A_over_tight": div(p["h_tight_iso_cluster_0"], p["h_tight_cluster_0"]),
                "current_B_over_tight": div(c["h_tight_noniso_cluster_0"], c["h_tight_cluster_0"]),
                "ppg12_B_over_tight": div(p["h_tight_noniso_cluster_0"], p["h_tight_cluster_0"]),
                "current_C_over_common": div(c["h_nontight_iso_cluster_0"], c["h_common_cluster_0"]),
                "ppg12_C_over_common": div(p["h_nontight_iso_cluster_0"], p["h_common_cluster_0"]),
                "current_D_over_common": div(c["h_nontight_noniso_cluster_0"], c["h_common_cluster_0"]),
                "ppg12_D_over_common": div(p["h_nontight_noniso_cluster_0"], p["h_common_cluster_0"]),
                "current_D_over_CplusD": div(
                    c["h_nontight_noniso_cluster_0"],
                    c["h_nontight_iso_cluster_0"] + c["h_nontight_noniso_cluster_0"],
                ),
                "ppg12_D_over_CplusD": div(
                    p["h_nontight_noniso_cluster_0"],
                    p["h_nontight_iso_cluster_0"] + p["h_nontight_noniso_cluster_0"],
                ),
            }
        )
        row["ratio_D_over_common_fraction"] = div(float(row["current_D_over_common"]), float(row["ppg12_D_over_common"]))
        row["ratio_D_over_CplusD_fraction"] = div(float(row["current_D_over_CplusD"]), float(row["ppg12_D_over_CplusD"]))
        rows.append(row)

    csv_path = args.outdir / f"the76_data_side_stage_hists_current_vs_ppg12_{args.tag}.csv"
    with csv_path.open("w", newline="") as handle:
        import csv

        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    md_path = args.outdir / f"the76_data_side_stage_hists_current_vs_ppg12_{args.tag}.md"
    lines = [
        "# THE-76 data-side stage histogram parity",
        "",
        f"Current ROOT: `{args.current_root}` / `{args.current_dir}`",
        f"PPG12 ROOT: `{ppg12_label}`",
        "",
        "## What this checks",
        "",
        "This compares aggregate pp-data histograms before the final `CalculatePhotonYield` solve: all clusters, common clusters, tight clusters, and A/B/C/D cluster populations.",
        "",
        "## Key ratios",
        "",
        "| bin | current common/all | PPG12 common/all | current tight/common | PPG12 tight/common | current D/common | PPG12 D/common | D/common ratio | current D/(C+D) | PPG12 D/(C+D) |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            f"| {row['bin']} | {fmt(float(row['current_common_over_all']))} | {fmt(float(row['ppg12_common_over_all']))} | "
            f"{fmt(float(row['current_tight_over_common']))} | {fmt(float(row['ppg12_tight_over_common']))} | "
            f"{fmt(float(row['current_D_over_common']))} | {fmt(float(row['ppg12_D_over_common']))} | "
            f"{fmt(float(row['ratio_D_over_common_fraction']))} | {fmt(float(row['current_D_over_CplusD']))} | "
            f"{fmt(float(row['ppg12_D_over_CplusD']))} |"
        )
    lines.extend(
        [
            "",
            "## Immediate readout",
            "",
            "- If common/all diverges first, the problem is event/candidate acceptance before photon ID.",
            "- If tight/common diverges first, the problem is BDT/tight-vs-nontight classification.",
            "- If D/common and D/(C+D) diverge after C is comparable, the problem is the non-tight isolation sideband definition or per-candidate isolation quantity.",
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
