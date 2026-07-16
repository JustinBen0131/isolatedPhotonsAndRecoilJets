#!/usr/bin/env python3
"""Read-only THE-94 jet8 PPG12 reference provenance audit.

This script does not submit, merge, transfer, or mutate SDCC state. It probes
small metadata from PPG12 reference ROOT files and source lists, then writes
local CSV/JSON evidence tables for the jet8-only Fig.6 discrepancy.
"""

from __future__ import annotations

import base64
import csv
import json
import os
import re
import shlex
import subprocess
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[3]
OUTDIR = (
    REPO
    / "dataOutput/ppg12Parity/THE-94_ppg12_pp_inclusivejet_sim_parity/"
    "jet8_reference_provenance_20260704"
)

HIST_NAME = "h_max_truth_jet_pT"
LO = 9.0
HI = 14.0
JET8_XSEC = 1.15e7
JET8_OLD_XSEC = 1.3013e7
JET50_XSEC = 7.3113
OBSERVED = 2.339817922469992
FORCED_XSEC = JET8_XSEC / OBSERVED

PPG12_REFS = {
    "true_period_combined": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_jet_bdt_nom.root",
    "period_combined_0mrad": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_jet_bdt_nom_0rad.root",
    "period_combined_1p5mrad": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_jet_bdt_nom_1p5mrad.root",
    "ian_no_suffix_jet8": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_jet8_bdt_nom.root",
    "jet8_component_0mrad_SI": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_jet8_nom_bdt_nom_0rad.root",
    "jet8_component_0mrad_DI": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_jet8_double_bdt_nom_0rad.root",
    "jet8_component_1p5mrad_SI": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_jet8_nom_bdt_nom_1p5mrad.root",
    "jet8_component_1p5mrad_DI": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_jet8_double_bdt_nom_1p5mrad.root",
    "jet8_nom_all_periods": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_jet8_nom_bdt_nom.root",
    "jet8_double_all_periods": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_jet8_double_bdt_nom.root",
}

PPG12_SLIMTREES = {
    "jet8_SI": "/sphenix/user/shuhangli/ppg12/FunWithxgboost/jet8/bdt_split.root",
    "jet8_DI": "/sphenix/user/shuhangli/ppg12/FunWithxgboost/jet8_double/bdt_split.root",
    "jet12_SI_control": "/sphenix/user/shuhangli/ppg12/FunWithxgboost/jet12/bdt_split.root",
    "jet12_DI_control": "/sphenix/user/shuhangli/ppg12/FunWithxgboost/jet12_double/bdt_split.root",
}

PPG12_SOURCE_DIRS = {
    "jet8_SI_run28": "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/jet8",
    "jet8_DI_run28": "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/jet8_double",
    "jet12_SI_run28_control": "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/jet12",
    "jet12_DI_run28_control": "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/jet12_double",
}


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def remote_python(payload: str) -> dict[str, Any]:
    encoded = base64.b64encode(payload.encode()).decode()
    remote_cmd = f"python3 -c 'import base64; exec(base64.b64decode(\"{encoded}\"))'"
    cmd = [
        "ssh",
        "patsfan753@ssh.sdcc.bnl.gov",
        (
            "ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null "
            f"sphnxuser05.sdcc.bnl.gov {shlex.quote(remote_cmd)}"
        ),
    ]
    env = os.environ.copy()
    sock = subprocess.check_output(["launchctl", "getenv", "SSH_AUTH_SOCK"], text=True).strip()
    env["SSH_AUTH_SOCK"] = sock
    proc = subprocess.run(cmd, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    (OUTDIR / "remote_stdout.txt").write_text(proc.stdout)
    (OUTDIR / "remote_stderr.txt").write_text(proc.stderr)
    if proc.returncode != 0:
        raise RuntimeError(f"remote probe failed rc={proc.returncode}; see {OUTDIR/'remote_stderr.txt'}")
    marker = "JSON_RESULT_START\n"
    idx = proc.stdout.rfind(marker)
    if idx < 0:
        raise RuntimeError(f"remote probe did not emit JSON marker; see {OUTDIR/'remote_stdout.txt'}")
    return json.loads(proc.stdout[idx + len(marker):])


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    payload = f"""
import json, os, re
from ctypes import c_double
import ROOT
ROOT.gROOT.SetBatch(True)

HIST_NAME = {HIST_NAME!r}
LO = {LO!r}
HI = {HI!r}
ppg12_refs = {json.dumps(PPG12_REFS)}
slimtrees = {json.dumps(PPG12_SLIMTREES)}
source_dirs = {json.dumps(PPG12_SOURCE_DIRS)}

def hist_metrics(label, path):
    row = dict(label=label, path=path, object=HIST_NAME, exists=False)
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        row["error"] = "zombie_or_missing"
        return row
    row["exists"] = True
    row["file_size_bytes"] = os.path.getsize(path) if os.path.exists(path) else None
    row["file_mtime"] = os.path.getmtime(path) if os.path.exists(path) else None
    h = f.Get(HIST_NAME)
    if not h:
        row["error"] = "hist_missing"
        f.Close()
        return row
    b1 = h.GetXaxis().FindBin(LO + 1e-9)
    b2 = h.GetXaxis().FindBin(HI - 1e-9)
    err = c_double(0.0)
    integral = h.IntegralAndError(b1, b2, err)
    row.update(
        nbins=int(h.GetNbinsX()),
        x_min=float(h.GetXaxis().GetXmin()),
        x_max=float(h.GetXaxis().GetXmax()),
        bin1=int(b1),
        bin2=int(b2),
        integral_9_14=float(integral),
        error_9_14=float(err.value),
        entries=float(h.GetEntries()),
        sumw_all=float(h.Integral(0, h.GetNbinsX()+1)),
        mean=float(h.GetMean()),
    )
    cfg = f.Get("config")
    if cfg:
        try:
            s = str(cfg.GetString().Data())
            row["config_lumi"] = _find_yaml_scalar(s, "lumi")
            row["config_lumi_target"] = _find_yaml_scalar(s, "lumi_target")
            row["config_truth_vertex_reweight_on"] = _find_yaml_scalar(s, "truth_vertex_reweight_on")
            row["config_truth_vertex_reweight_file"] = _find_yaml_scalar(s, "truth_vertex_reweight_file")
            row["config_var_type"] = _find_yaml_scalar(s, "var_type")
        except Exception as exc:
            row["config_error"] = repr(exc)
    for key in ["merge_lumi_period0_bdt_nom_0rad", "merge_lumi_period1_bdt_nom_1p5mrad", "merge_lumi_target"]:
        obj = f.Get(key)
        if obj:
            try:
                row[key] = str(obj.GetTitle())
            except Exception:
                pass
    f.Close()
    return row

def _find_yaml_scalar(text, key):
    pat = re.compile(r"^\\s*" + re.escape(key) + r":\\s*(.+?)\\s*$", re.M)
    m = pat.search(text)
    if not m:
        return ""
    return m.group(1).strip().strip('"')

def tree_metrics(label, path):
    row = dict(label=label, path=path, exists=False)
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        row["error"] = "zombie_or_missing"
        return row
    row["exists"] = True
    row["file_size_bytes"] = os.path.getsize(path) if os.path.exists(path) else None
    row["file_mtime"] = os.path.getmtime(path) if os.path.exists(path) else None
    t = f.Get("slimtree")
    if not t:
        row["error"] = "slimtree_missing"
        f.Close()
        return row
    row["entries"] = int(t.GetEntries())
    row["branches"] = ",".join([br.GetName() for br in t.GetListOfBranches()][:80])
    f.Close()
    return row

def list_metrics(label, d):
    rows = []
    if not os.path.isdir(d):
        return [dict(label=label, dir=d, exists=False)]
    names = sorted([n for n in os.listdir(d) if n.endswith(".list") or n.endswith(".C") or n.endswith(".sh")])
    for name in names:
        path = os.path.join(d, name)
        row = dict(label=label, dir=d, file=name, path=path, exists=True)
        if os.path.isfile(path):
            try:
                with open(path, "r", errors="replace") as fh:
                    lines = [ln.rstrip("\\n") for ln in fh]
                nonblank = [ln for ln in lines if ln.strip() and not ln.strip().startswith("#")]
                row["line_count"] = len(lines)
                row["nonblank_noncomment_count"] = len(nonblank)
                row["first_noncomment"] = nonblank[0] if nonblank else ""
                row["mentions_dataset28"] = any("0000000028" in ln for ln in lines)
                row["mentions_dataset29"] = any("0000000029" in ln for ln in lines)
                row["mentions_bseidlitz"] = any("bseidlitz" in ln for ln in lines)
                row["mentions_g4hits"] = any("G4Hits" in ln or "g4hits" in ln for ln in lines)
                row["mentions_dst_allreco"] = any("DST_ALLRECO" in ln for ln in lines)
                if name == "Fun4All_run_sim.C":
                    row["active_input_lines"] = "; ".join([ln.strip() for ln in lines if "inputFile" in ln or "listfile" in ln][:16])
            except Exception as exc:
                row["error"] = repr(exc)
        rows.append(row)
    return rows

hist_rows = [hist_metrics(k, v) for k, v in ppg12_refs.items()]
tree_rows = [tree_metrics(k, v) for k, v in slimtrees.items()]
list_rows = []
for k, v in source_dirs.items():
    list_rows.extend(list_metrics(k, v))

print("JSON_RESULT_START")
print(json.dumps(dict(hist_rows=hist_rows, tree_rows=tree_rows, list_rows=list_rows), sort_keys=True))
"""
    data = remote_python(payload)

    hist_rows = data["hist_rows"]
    tree_rows = data["tree_rows"]
    list_rows = data["list_rows"]

    write_csv(
        OUTDIR / "ppg12_jet8_reference_root_metrics.csv",
        hist_rows,
        [
            "label",
            "path",
            "object",
            "exists",
            "integral_9_14",
            "error_9_14",
            "entries",
            "sumw_all",
            "nbins",
            "x_min",
            "x_max",
            "config_lumi",
            "config_lumi_target",
            "config_truth_vertex_reweight_on",
            "config_truth_vertex_reweight_file",
            "config_var_type",
            "merge_lumi_period0_bdt_nom_0rad",
            "merge_lumi_period1_bdt_nom_1p5mrad",
            "merge_lumi_target",
            "file_size_bytes",
            "file_mtime",
            "error",
        ],
    )
    write_csv(
        OUTDIR / "ppg12_jet8_slimtree_metrics.csv",
        tree_rows,
        ["label", "path", "exists", "entries", "file_size_bytes", "file_mtime", "branches", "error"],
    )
    write_csv(
        OUTDIR / "ppg12_jet8_source_list_metrics.csv",
        list_rows,
        [
            "label",
            "dir",
            "file",
            "path",
            "exists",
            "line_count",
            "nonblank_noncomment_count",
            "first_noncomment",
            "mentions_dataset28",
            "mentions_dataset29",
            "mentions_bseidlitz",
            "mentions_g4hits",
            "mentions_dst_allreco",
            "active_input_lines",
            "error",
        ],
    )

    def by_label(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
        return {str(r["label"]): r for r in rows}

    h = by_label(hist_rows)
    t = by_label(tree_rows)
    ratios: list[dict[str, Any]] = []
    comp_pairs = [
        ("0mrad_SI", "jet8_component_0mrad_SI", "jet8_SI"),
        ("0mrad_DI", "jet8_component_0mrad_DI", "jet8_DI"),
        ("1p5mrad_SI", "jet8_component_1p5mrad_SI", "jet8_SI"),
        ("1p5mrad_DI", "jet8_component_1p5mrad_DI", "jet8_DI"),
    ]
    for comp, hist_label, tree_label in comp_pairs:
        hr = h.get(hist_label, {})
        tr = t.get(tree_label, {})
        integ = float(hr.get("integral_9_14") or 0.0)
        entries = float(hr.get("entries") or 0.0)
        source_entries = float(tr.get("entries") or 0.0)
        ratios.append(
            {
                "quantity": f"ppg12_{comp}_hist_entries_over_slimtree_entries",
                "value": entries / source_entries if source_entries else None,
                "explains_2p338": False,
                "note": "Raw in-window fill population fraction in PPG12 source.",
            }
        )
        ratios.append(
            {
                "quantity": f"ppg12_{comp}_integral_per_hist_entry",
                "value": integ / entries if entries else None,
                "explains_2p338": False,
                "note": "Average weighted fill per raw PPG12 in-window entry.",
            }
        )
    ratios.extend(
        [
            {
                "quantity": "observed_current_over_ppg12_true_period",
                "value": OBSERVED,
                "explains_2p338": True,
                "note": "Measured final jet8 failure.",
            },
            {
                "quantity": "known_old_xsec_over_ppg12_xsec",
                "value": JET8_OLD_XSEC / JET8_XSEC,
                "explains_2p338": False,
                "note": "Known old jet8 xsec issue is too small.",
            },
            {
                "quantity": "forced_effective_jet8_xsec_pb",
                "value": FORCED_XSEC,
                "explains_2p338": True,
                "note": "Numerical-only implied xsec; not found as a PPG12 source constant.",
            },
        ]
    )
    write_csv(OUTDIR / "ppg12_jet8_real_ratio_search.csv", ratios, ["quantity", "value", "explains_2p338", "note"])

    manifest = {
        "task": "THE-94 jet8 PPG12 reference provenance audit",
        "output_dir": str(OUTDIR),
        "root_metrics": str(OUTDIR / "ppg12_jet8_reference_root_metrics.csv"),
        "slimtree_metrics": str(OUTDIR / "ppg12_jet8_slimtree_metrics.csv"),
        "source_list_metrics": str(OUTDIR / "ppg12_jet8_source_list_metrics.csv"),
        "real_ratio_search": str(OUTDIR / "ppg12_jet8_real_ratio_search.csv"),
        "observed_factor": OBSERVED,
        "known_old_xsec_factor": JET8_OLD_XSEC / JET8_XSEC,
        "forced_effective_xsec_pb": FORCED_XSEC,
    }
    (OUTDIR / "ppg12_jet8_reference_provenance_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
