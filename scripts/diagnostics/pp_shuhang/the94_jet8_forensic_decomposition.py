#!/usr/bin/env python3
"""THE-94 jet8-only inclusive Fig.6 forensic decomposition.

This is a read-only diagnostic.  It compares the accepted THE-94 jet8 output
against PPG12 component/reference products and writes local CSV/JSON evidence.
It does not submit Condor, merge, transfer payloads, edit SDCC files, or apply
any fitted scale.
"""

from __future__ import annotations

import base64
import csv
import json
import math
import os
import shlex
import subprocess
from pathlib import Path
from statistics import mean, median
from typing import Any

import numpy as np
import uproot


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
PYTHON = Path("/Users/patsfan753/Desktop/analysis/env/bin/python")
OUT_DIR = (
    REPO
    / "dataOutput/ppg12Parity/THE-94_ppg12_pp_inclusivejet_sim_parity/jet8_forensics_20260704"
)

FINAL_ROOT = (
    REPO
    / "dataOutput/ppg12Parity/the94_ppg12_inclusivejet_fig6_fixed_20260702_204032/"
    "final_roots/inclusivejet/"
    "RecoilJets_jet8_ALL_jetMinPtScan_dphiScan_vz60_isoR40_isSliding_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12.root"
)
TRUE_PERIOD_CSV = (
    REPO
    / "dataOutput/ppg12Parity/the94_ppg12_inclusivejet_fig6_fixed_20260702_204032/"
    "contract_canary/fig6_inclusive/fig6_inclusive_current_kept_vs_ppg12_true_period_combined_points.csv"
)
NO_SUFFIX_CSV = (
    REPO
    / "dataOutput/ppg12Parity/the94_ppg12_inclusivejet_fig6_fixed_20260702_204032/"
    "contract_canary/fig6_inclusive/fig6_inclusive_current_kept_vs_ppg12_source_points.csv"
)

HIST_KEPT = "SIM/h_ppInclusiveJetStitch_ppg12Fig6EventWeighted_r04_maxTruthJetPt_kept"
HIST_META = "SIM/h_ppInclusiveJetStitch_ppg12Fig6EventWeighted_metadata"

JET8_XSEC = 1.15e7
JET8_OLD_XSEC = 1.3013e7
JET50_XSEC = 7.3113
OBSERVED_FACTOR = 2.339817922469992
EFFECTIVE_XSEC_IF_FORCED = JET8_XSEC / OBSERVED_FACTOR

# Current all-z nominal period weights from the PPG12 merge-feeder configs
# (config-schema.md: 0 mrad 47.2076, 1.5 mrad 17.1642, target 64.3718 pb^-1).
LUMI_0MRAD = 47.2076
LUMI_1P5MRAD = 17.1642
LUMI_TARGET = 64.3718

# Historical 60 cm factors from PPG12 sample_combining_check.tex. Kept separate
# because that report explains the old Blair jet8 SI anomaly, not the current
# all-z true-period reference used by THE-94.
LUMI_0MRAD_60CM = 32.66
LUMI_1P5MRAD_60CM = 16.27
LUMI_TARGET_60CM = LUMI_0MRAD_60CM + LUMI_1P5MRAD_60CM
MIX_DOUBLE_0MRAD = 0.224
MIX_DOUBLE_1P5MRAD = 0.079
JET8_SI_VTX_MEAN_0MRAD_PPG12_REPORT = 1.192
JET8_SI_VTX_MEAN_1P5MRAD_PPG12_REPORT = 2.948


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def summarize(values: list[float]) -> dict[str, Any]:
    clean = [float(v) for v in values if math.isfinite(float(v))]
    if not clean:
        return {"n": 0, "median": None, "mean": None, "min": None, "max": None}
    return {
        "n": len(clean),
        "median": median(clean),
        "mean": mean(clean),
        "min": min(clean),
        "max": max(clean),
    }


def local_hist_integral(path: Path, hist_name: str, lo: float = 9.0, hi: float = 14.0) -> dict[str, Any]:
    with uproot.open(path) as root_file:
        hist = root_file[hist_name]
        values, edges = hist.to_numpy(flow=False)
        variances = hist.variances(flow=False)
        centers = 0.5 * (edges[:-1] + edges[1:])
        mask = (centers >= lo) & (centers < hi)
        err2 = np.zeros_like(values, dtype=float) if variances is None else np.clip(variances, 0.0, None)
        return {
            "integral": float(values[mask].sum()),
            "error": float(np.sqrt(err2[mask].sum())),
            "entries": float(getattr(hist, "member", lambda _: np.nan)("fEntries")),
            "nonzero_bins": int(np.count_nonzero(values[mask])),
        }


def local_metadata(path: Path) -> dict[str, float]:
    with uproot.open(path) as root_file:
        values, _ = root_file[HIST_META].to_numpy(flow=False)
    labels = [
        "events_seen_raw",
        "window_low_GeV",
        "window_high_GeV",
        "upper_edge_inclusive",
        "bin_width_GeV",
        "xsec_pb",
        "xsec_over_jet50",
        "sample_bin",
        "truth_def_code",
        "includes_current_weight",
        "includes_slice_weight",
        "value_mode_code",
    ]
    return {label: float(values[i]) for i, label in enumerate(labels) if i < len(values)}


def jet8_ratios_from_csv(path: Path, ratio_col: str) -> dict[str, Any]:
    ratios: list[float] = []
    current_integral = 0.0
    ref_integral = 0.0
    for row in read_csv(path):
        if row.get("sample") != "jet8":
            continue
        ratios.append(float(row[ratio_col]))
        current_integral += float(row.get("current_kept_value", row.get("current_value", 0.0)))
        if "ppg12_true_period_value" in row:
            ref_integral += float(row["ppg12_true_period_value"])
        elif "ppg12_source_value" in row:
            ref_integral += float(row["ppg12_source_value"])
    stats = summarize(ratios)
    stats.update({"current_integral": current_integral, "reference_integral": ref_integral})
    return stats


def run_remote_probe() -> dict[str, Any]:
    cached = OUT_DIR / "remote_probe.json"
    if cached.exists() and os.environ.get("THE94_JET8_FORCE_REMOTE_PROBE", "0") != "1":
        return json.loads(cached.read_text())

    remote_code = r'''
import glob
import json
import math
import os
from pathlib import Path

import ROOT

ROOT.gROOT.SetBatch(True)

LO = 9.0
HI = 14.0
RESULTS = Path("/sphenix/user/shuhangli/ppg12/efficiencytool/results")
ORIGINAL_BASE = Path("/sphenix/tg/tg01/bulk/jbennett/thesisAna/siminclusive/the94_ppg12_inclusivejet_fig6_fixed_20260702_204032")
STYLE = "jetMinPtScan_dphiScan_vz60_isoR40_isSliding_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12"
CURRENT_HIST = "SIM/h_ppInclusiveJetStitch_ppg12Fig6EventWeighted_r04_maxTruthJetPt_kept"
CURRENT_META = "SIM/h_ppInclusiveJetStitch_ppg12Fig6EventWeighted_metadata"
REF_HIST = "h_max_truth_jet_pT"


def open_hist(path, hist):
    f = ROOT.TFile.Open(str(path), "READ")
    if not f or f.IsZombie():
        return None, None, "zombie_or_unopenable"
    h = f.Get(hist)
    if not h:
        f.Close()
        return None, None, "missing_hist"
    return f, h, ""


def integrate_file(path, hist):
    f, h, err = open_hist(path, hist)
    if err:
        return {"path": str(path), "ok": False, "error": err}
    total = 0.0
    err2 = 0.0
    nonzero = 0
    for i in range(1, h.GetNbinsX() + 1):
        center = h.GetXaxis().GetBinCenter(i)
        if LO <= center < HI:
            v = float(h.GetBinContent(i))
            e = float(h.GetBinError(i))
            total += v
            err2 += e * e
            if v != 0.0:
                nonzero += 1
    entries = float(h.GetEntries())
    f.Close()
    return {
        "path": str(path),
        "ok": True,
        "integral": total,
        "error": math.sqrt(err2),
        "entries": entries,
        "nonzero_bins": nonzero,
    }


def metadata_first(path):
    f, h, err = open_hist(path, CURRENT_META)
    if err:
        return {"error": err}
    labels = [
        "events_seen_raw",
        "window_low_GeV",
        "window_high_GeV",
        "upper_edge_inclusive",
        "bin_width_GeV",
        "xsec_pb",
        "xsec_over_jet50",
        "sample_bin",
        "truth_def_code",
        "includes_current_weight",
        "includes_slice_weight",
        "value_mode_code",
    ]
    out = {}
    for idx, label in enumerate(labels, start=1):
        if idx <= h.GetNbinsX():
            out[label] = float(h.GetBinContent(idx))
    f.Close()
    return out


def sum_current_component(period, component, sample_dir):
    d = ORIGINAL_BASE / f"{period}_{component}" / STYLE / sample_dir
    files = sorted(glob.glob(str(d / "*.root")))
    total = 0.0
    err2 = 0.0
    entries = 0.0
    nonzero_files = 0
    missing = 0
    zombie = 0
    min_size = None
    max_size = None
    first_meta = {}
    for n, path in enumerate(files):
        try:
            size = os.path.getsize(path)
            min_size = size if min_size is None else min(min_size, size)
            max_size = size if max_size is None else max(max_size, size)
        except OSError:
            pass
        rec = integrate_file(path, CURRENT_HIST)
        if not rec["ok"]:
            if rec["error"] == "missing_hist":
                missing += 1
            else:
                zombie += 1
            continue
        total += rec["integral"]
        err2 += rec["error"] * rec["error"]
        entries += rec["entries"]
        if rec["integral"] != 0.0:
            nonzero_files += 1
        if n == 0:
            first_meta = metadata_first(path)
    return {
        "path": str(d),
        "files": len(files),
        "nonzero_files": nonzero_files,
        "missing_hist_files": missing,
        "zombie_or_unopenable_files": zombie,
        "min_size": min_size,
        "max_size": max_size,
        "integral": total,
        "error": math.sqrt(err2),
        "entries": entries,
        "metadata_first_file": first_meta,
    }


def read_text(path, max_lines=16):
    try:
        with open(path, "r", errors="replace") as handle:
            return "".join(handle.readlines()[:max_lines])
    except Exception as exc:
        return f"UNREADABLE: {exc}"


components = {
    "0mrad_SI": ("0mrad", "single", "run28_jet8", RESULTS / "MC_efficiency_jet8_nom_bdt_nom_0rad.root"),
    "0mrad_DI": ("0mrad", "double", "run28_jet8_double", RESULTS / "MC_efficiency_jet8_double_bdt_nom_0rad.root"),
    "1p5mrad_SI": ("1p5mrad", "single", "run28_jet8", RESULTS / "MC_efficiency_jet8_nom_bdt_nom_1p5mrad.root"),
    "1p5mrad_DI": ("1p5mrad", "double", "run28_jet8_double", RESULTS / "MC_efficiency_jet8_double_bdt_nom_1p5mrad.root"),
}

component_rows = []
for label, (period, comp, sample_dir, ref_path) in components.items():
    cur = sum_current_component(period, comp, sample_dir)
    ref = integrate_file(ref_path, REF_HIST)
    row = {
        "component_label": label,
        "period": period,
        "component": "SI" if comp == "single" else "DI",
        "sample_dir": sample_dir,
        "current": cur,
        "ppg12_reference": ref,
        "current_over_ppg12": (cur["integral"] / ref["integral"]) if ref.get("ok") and ref.get("integral") else None,
        "ppg12_reference_path": str(ref_path),
    }
    component_rows.append(row)

reference_names = {
    "true_period_combined_product": RESULTS / "MC_efficiency_jet_bdt_nom.root",
    "period_combined_0mrad": RESULTS / "MC_efficiency_jet_bdt_nom_0rad.root",
    "period_combined_1p5mrad": RESULTS / "MC_efficiency_jet_bdt_nom_1p5mrad.root",
    "ian_no_suffix_per_sample": RESULTS / "MC_efficiency_jet8_bdt_nom.root",
    "jet8_nom_all_periods": RESULTS / "MC_efficiency_jet8_nom_bdt_nom.root",
    "jet8_double_all_periods": RESULTS / "MC_efficiency_jet8_double_bdt_nom.root",
    "jet8_nom_0mrad": RESULTS / "MC_efficiency_jet8_nom_bdt_nom_0rad.root",
    "jet8_nom_1p5mrad": RESULTS / "MC_efficiency_jet8_nom_bdt_nom_1p5mrad.root",
    "jet8_double_0mrad": RESULTS / "MC_efficiency_jet8_double_bdt_nom_0rad.root",
    "jet8_double_1p5mrad": RESULTS / "MC_efficiency_jet8_double_bdt_nom_1p5mrad.root",
    "jet8_nom_1p5mrad_norewt": RESULTS / "MC_efficiency_jet8_nom_bdt_nom_1p5mrad_norewt.root",
    "jet8_nom_1p5mrad_newrewt": RESULTS / "MC_efficiency_jet8_nom_bdt_nom_1p5mrad_newrewt.root",
    "jet8_nomold_0mrad": RESULTS / "MC_efficiency_jet8_bdt_nomold_0rad.root",
    "jet8_nomold_1p5mrad": RESULTS / "MC_efficiency_jet8_bdt_nomold_1p5mrad.root",
}
reference_rows = []
for kind, path in reference_names.items():
    rec = integrate_file(path, REF_HIST)
    rec["reference_kind"] = kind
    rec["object_name"] = REF_HIST
    rec["status"] = (
        "main_target" if kind == "true_period_combined_product"
        else "historical_ian_reference" if kind == "ian_no_suffix_per_sample"
        else "component_or_variant_reference"
    )
    reference_rows.append(rec)

remote_ppg12_jet8_dir = Path("/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/jet8")
source_probe = {
    "remote_ppg12_jet8_dir": str(remote_ppg12_jet8_dir),
    "exists": remote_ppg12_jet8_dir.exists(),
    "fun4all_excerpt": read_text(remote_ppg12_jet8_dir / "Fun4All_run_sim.C", 120),
    "test_list_excerpt": read_text(remote_ppg12_jet8_dir / "test.list", 8),
    "g4hits_list_excerpt": read_text(remote_ppg12_jet8_dir / "g4hits.list", 3),
}

out = {
    "components": component_rows,
    "references": reference_rows,
    "source_probe": source_probe,
}
print("__THE94_JSON_START__")
print(json.dumps(out, sort_keys=True))
print("__THE94_JSON_END__")
'''
    encoded = base64.b64encode(remote_code.encode()).decode()
    inner = (
        "cd /sphenix/u/patsfan753/scratch/thesisAnalysis && "
        "python3 - <<'PY'\n"
        "import base64\n"
        f"code = base64.b64decode('{encoded}').decode()\n"
        "exec(code)\n"
        "PY"
    )
    remote_gateway_cmd = (
        "ssh -o BatchMode=yes -o StrictHostKeyChecking=no "
        "-o UserKnownHostsFile=/dev/null sphnxuser05.sdcc.bnl.gov "
        + shlex.quote(inner)
    )
    env = os.environ.copy()
    if not env.get("SSH_AUTH_SOCK"):
        sock = subprocess.run(
            ["launchctl", "getenv", "SSH_AUTH_SOCK"],
            check=False,
            text=True,
            capture_output=True,
        ).stdout.strip()
        if sock:
            env["SSH_AUTH_SOCK"] = sock
    result = subprocess.run(
        ["ssh", "-o", "BatchMode=yes", "patsfan753@ssh.sdcc.bnl.gov", remote_gateway_cmd],
        check=False,
        text=True,
        capture_output=True,
        env=env,
    )
    (OUT_DIR / "remote_probe_stdout.txt").write_text(result.stdout)
    (OUT_DIR / "remote_probe_stderr.txt").write_text(result.stderr)
    if result.returncode != 0:
        raise RuntimeError(
            f"remote probe failed with code {result.returncode}; "
            f"see {OUT_DIR / 'remote_probe_stdout.txt'} and {OUT_DIR / 'remote_probe_stderr.txt'}"
        )
    start = result.stdout.find("__THE94_JSON_START__")
    end = result.stdout.find("__THE94_JSON_END__")
    if start < 0 or end < 0 or end <= start:
        raise RuntimeError(f"remote JSON markers missing; see {OUT_DIR / 'remote_probe_stdout.txt'}")
    payload = result.stdout[start + len("__THE94_JSON_START__") : end].strip()
    return json.loads(payload)


def component_factor(period: str, component: str) -> tuple[float, float]:
    if period == "0mrad":
        lumi = LUMI_0MRAD / LUMI_TARGET
        mix = 1.0 - MIX_DOUBLE_0MRAD if component == "SI" else MIX_DOUBLE_0MRAD
    elif period == "1p5mrad":
        lumi = LUMI_1P5MRAD / LUMI_TARGET
        mix = 1.0 - MIX_DOUBLE_1P5MRAD if component == "SI" else MIX_DOUBLE_1P5MRAD
    else:
        raise ValueError(period)
    return lumi, mix


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    remote = run_remote_probe()
    (OUT_DIR / "remote_probe.json").write_text(json.dumps(remote, indent=2, sort_keys=True))

    final = local_hist_integral(FINAL_ROOT, HIST_KEPT)
    final_meta = local_metadata(FINAL_ROOT)
    true_stats = jet8_ratios_from_csv(TRUE_PERIOD_CSV, "current_over_ppg12_true_period")
    source_stats = jet8_ratios_from_csv(NO_SUFFIX_CSV, "current_kept_over_ppg12")

    component_rows: list[dict[str, Any]] = []
    for rec in remote["components"]:
        label = rec["component_label"]
        period = rec["period"]
        component = rec["component"]
        current = rec["current"]
        ref = rec["ppg12_reference"]
        lumi_factor, mix_factor = component_factor(period, component)
        vtx_summary = "not_recorded_in_current_root"
        if label == "0mrad_SI":
            vtx_summary = f"PPG12_report_jet8_SI_mean_truth_vtx_weight={JET8_SI_VTX_MEAN_0MRAD_PPG12_REPORT}"
        elif label == "1p5mrad_SI":
            vtx_summary = f"PPG12_report_jet8_SI_mean_truth_vtx_weight={JET8_SI_VTX_MEAN_1P5MRAD_PPG12_REPORT}"
        elif component == "DI":
            vtx_summary = "PPG12_report_DI_broad_truth_z_mean_reweight_approximately_1"
        meta = current.get("metadata_first_file") or {}
        component_rows.append(
            {
                "component": label,
                "period": period,
                "si_di": component,
                "source_input_path_family": current["path"],
                "source_family_classification": (
                    "current_RecoilJets_standard_run28_dataset28_nominal"
                    if component == "SI"
                    else "current_RecoilJets_run28_double_DI_nominal"
                ),
                "current_integral": current["integral"],
                "current_error": current["error"],
                "ppg12_reference_integral": ref.get("integral"),
                "ppg12_reference_error": ref.get("error"),
                "current_over_ppg12_component": rec.get("current_over_ppg12"),
                "raw_entries_hist_getentries": current["entries"],
                "weighted_entries_integral": current["integral"],
                "xsec_used_pb": meta.get("xsec_pb", JET8_XSEC),
                "xsec_over_jet50": meta.get("xsec_over_jet50", JET8_XSEC / JET50_XSEC),
                "lumi_factor": lumi_factor,
                "mix_factor": mix_factor,
                "vertex_weight_summary": vtx_summary,
                "root_files_merged": current["files"],
                "nonzero_files": current["nonzero_files"],
                "missing_hist_files": current["missing_hist_files"],
                "zombie_or_unopenable_files": current["zombie_or_unopenable_files"],
                "ppg12_reference_path": rec["ppg12_reference_path"],
                "ppg12_reference_object": "h_max_truth_jet_pT",
            }
        )

    reference_rows: list[dict[str, Any]] = []
    for rec in remote["references"]:
        kind = rec["reference_kind"]
        path = rec["path"]
        source_family = "unknown"
        if "nom_bdt_nom" in path:
            source_family = "jet8_nom_SI_component"
        if "double" in path:
            source_family = "jet8_double_DI_component"
        if kind == "ian_no_suffix_per_sample":
            source_family = "historical_IAN_no_suffix_per_sample"
        if kind == "true_period_combined_product":
            source_family = "true_period_combined_jet8plus12plus20plus30plus40"
        reference_rows.append(
            {
                "reference_kind": kind,
                "path": path,
                "object_name": rec.get("object_name", "h_max_truth_jet_pT"),
                "source_sample_family": source_family,
                "xsec_convention": "PPG12 CrossSectionWeights.h jet8cross=1.15e7, jet50cross=7.3113",
                "event_or_source_count": "not encoded in histogram-only reference; component integral recorded",
                "owned_window_integral_9_14": rec.get("integral"),
                "owned_window_error": rec.get("error"),
                "entries_hist_getentries": rec.get("entries"),
                "nonzero_bins": rec.get("nonzero_bins"),
                "status": rec.get("status"),
                "nominal_or_historical": (
                    "nominal_main_target"
                    if kind == "true_period_combined_product"
                    else "historical_ian_reference"
                    if kind == "ian_no_suffix_per_sample"
                    else "component_or_variant_diagnostic"
                ),
            }
        )

    refs = {row["reference_kind"]: row for row in reference_rows}
    true_integral = refs["true_period_combined_product"]["owned_window_integral_9_14"]
    no_suffix_integral = refs["ian_no_suffix_per_sample"]["owned_window_integral_9_14"]
    nom0 = refs["jet8_nom_0mrad"]["owned_window_integral_9_14"]
    nom1 = refs["jet8_nom_1p5mrad"]["owned_window_integral_9_14"]
    dbl0 = refs["jet8_double_0mrad"]["owned_window_integral_9_14"]
    dbl1 = refs["jet8_double_1p5mrad"]["owned_window_integral_9_14"]
    comp_current_sum = sum(float(row["current_integral"]) for row in component_rows)

    predicted_jet8_si_period_ratio = (
        (LUMI_1P5MRAD * (1.0 - MIX_DOUBLE_1P5MRAD))
        / (LUMI_0MRAD * (1.0 - MIX_DOUBLE_0MRAD))
        * (JET8_SI_VTX_MEAN_1P5MRAD_PPG12_REPORT / JET8_SI_VTX_MEAN_0MRAD_PPG12_REPORT)
    )
    expected_wide_si_period_ratio = (
        (LUMI_1P5MRAD * (1.0 - MIX_DOUBLE_1P5MRAD))
        / (LUMI_0MRAD * (1.0 - MIX_DOUBLE_0MRAD))
    )
    expected_di_period_ratio = (
        (LUMI_1P5MRAD * MIX_DOUBLE_1P5MRAD)
        / (LUMI_0MRAD * MIX_DOUBLE_0MRAD)
    )
    historical_wide_si_period_ratio_60cm = (
        (LUMI_1P5MRAD_60CM * (1.0 - MIX_DOUBLE_1P5MRAD))
        / (LUMI_0MRAD_60CM * (1.0 - MIX_DOUBLE_0MRAD))
    )
    historical_predicted_jet8_si_period_ratio_60cm = (
        historical_wide_si_period_ratio_60cm
        * (JET8_SI_VTX_MEAN_1P5MRAD_PPG12_REPORT / JET8_SI_VTX_MEAN_0MRAD_PPG12_REPORT)
    )
    scalar_rows = [
        {
            "candidate_factor": "observed_current_over_true_period",
            "value": OBSERVED_FACTOR,
            "explains_2p3398": True,
            "defensible": "observed_failure_not_correction",
            "note": "Measured flat jet8 current/PPG12 factor from true-period-combined no-scale check.",
        },
        {
            "candidate_factor": "known_old_xsec_over_ppg12_xsec",
            "value": JET8_OLD_XSEC / JET8_XSEC,
            "explains_2p3398": False,
            "defensible": "real_known_old_xsec_issue_but_too_small",
            "note": "1.3013e7 / 1.15e7 = 1.1316, not 2.3398.",
        },
        {
            "candidate_factor": "forced_effective_xsec_pb_needed",
            "value": EFFECTIVE_XSEC_IF_FORCED,
            "explains_2p3398": True,
            "defensible": "not_defensible_without_ppg12_source",
            "note": "Would require jet8 xsec about 4.91e6 pb; not found in PPG12 CrossSectionWeights.h.",
        },
        {
            "candidate_factor": "ppg12_true_over_no_suffix_reference_integral",
            "value": true_integral / no_suffix_integral if no_suffix_integral else None,
            "explains_2p3398": False,
            "defensible": "real_reference_identity_factor",
            "note": "Shows no-suffix and true-period-combined are not the same reference target.",
        },
        {
            "candidate_factor": "ppg12_nom_1p5_over_0mrad_measured",
            "value": nom1 / nom0 if nom0 else None,
            "explains_2p3398": False,
            "defensible": "real_ppg12_jet8_SI_anomaly",
            "note": "PPG12 jet8 SI period ratio is anomalous due to documented narrow truth-z source.",
        },
        {
            "candidate_factor": "ppg12_expected_wide_SI_1p5_over_0mrad",
            "value": expected_wide_si_period_ratio,
            "explains_2p3398": False,
            "defensible": "real_allz_policy_factor_for_normal_samples",
            "note": "Expected all-z period ratio for wide-z SI samples.",
        },
        {
            "candidate_factor": "historical_60cm_ppg12_report_predicted_jet8_SI_anomalous_period_ratio",
            "value": historical_predicted_jet8_si_period_ratio_60cm,
            "explains_2p3398": False,
            "defensible": "real_historical_exception_factor_not_current_allz_target",
            "note": "Old 60 cm report factor: 0.591 * 2.948/1.192 = 1.462. Historical jet8 SI anomaly, not the current all-z scalar.",
        },
        {
            "candidate_factor": "current_allz_if_old_jet8_narrow_z_means_applied",
            "value": predicted_jet8_si_period_ratio,
            "explains_2p3398": False,
            "defensible": "diagnostic_only",
            "note": "All-z analogue of the old narrow-z boost. Current PPG12 component refs instead track the wide-z all-z period ratio.",
        },
        {
            "candidate_factor": "ppg12_double_1p5_over_0mrad_measured",
            "value": dbl1 / dbl0 if dbl0 else None,
            "explains_2p3398": False,
            "defensible": "real_DI_factor",
            "note": "DI is not the source of the observed flat 2.34 excess if this is near 0.176.",
        },
        {
            "candidate_factor": "ppg12_expected_DI_1p5_over_0mrad",
            "value": expected_di_period_ratio,
            "explains_2p3398": False,
            "defensible": "real_policy_factor_for_DI_samples",
            "note": "Expected period ratio for DI samples.",
        },
        {
            "candidate_factor": "current_component_sum_over_final_root",
            "value": comp_current_sum / final["integral"] if final["integral"] else None,
            "explains_2p3398": False,
            "defensible": "merge_composition_check",
            "note": "Checks whether local final root is just the sum of the four current components.",
        },
    ]

    source_probe = remote["source_probe"]
    source_probe_summary = {
        "remote_ppg12_jet8_dir": source_probe.get("remote_ppg12_jet8_dir"),
        "exists": source_probe.get("exists"),
        "fun4all_contains_test_list": "test.list" in source_probe.get("fun4all_excerpt", ""),
        "test_list_mentions_bseidlitz": "bseidlitz" in source_probe.get("test_list_excerpt", ""),
        "g4hits_list_available": bool(source_probe.get("g4hits_list_excerpt", "").strip())
        and "UNREADABLE" not in source_probe.get("g4hits_list_excerpt", ""),
        "note": "If remote source files are not readable, use ppg12codeGit/reports/sample_combining_check.tex as the source identity evidence.",
    }

    write_csv(
        OUT_DIR / "jet8_component_table.csv",
        component_rows,
        [
            "component",
            "period",
            "si_di",
            "source_input_path_family",
            "source_family_classification",
            "current_integral",
            "current_error",
            "ppg12_reference_integral",
            "ppg12_reference_error",
            "current_over_ppg12_component",
            "raw_entries_hist_getentries",
            "weighted_entries_integral",
            "xsec_used_pb",
            "xsec_over_jet50",
            "lumi_factor",
            "mix_factor",
            "vertex_weight_summary",
            "root_files_merged",
            "nonzero_files",
            "missing_hist_files",
            "zombie_or_unopenable_files",
            "ppg12_reference_path",
            "ppg12_reference_object",
        ],
    )
    write_csv(
        OUT_DIR / "jet8_reference_identity_table.csv",
        reference_rows,
        [
            "reference_kind",
            "path",
            "object_name",
            "source_sample_family",
            "xsec_convention",
            "event_or_source_count",
            "owned_window_integral_9_14",
            "owned_window_error",
            "entries_hist_getentries",
            "nonzero_bins",
            "status",
            "nominal_or_historical",
        ],
    )
    write_csv(
        OUT_DIR / "jet8_scalar_factor_diagnosis.csv",
        scalar_rows,
        ["candidate_factor", "value", "explains_2p3398", "defensible", "note"],
    )

    manifest = {
        "task": "THE-94 jet8-only inclusive Fig.6 forensic decomposition",
        "output_dir": str(OUT_DIR),
        "final_root": str(FINAL_ROOT),
        "final_integral_9_14": final,
        "final_metadata_raw": final_meta,
        "ratio_stats_true_period": true_stats,
        "ratio_stats_no_suffix": source_stats,
        "source_probe_summary": source_probe_summary,
        "component_table": str(OUT_DIR / "jet8_component_table.csv"),
        "reference_identity_table": str(OUT_DIR / "jet8_reference_identity_table.csv"),
        "scalar_factor_diagnosis": str(OUT_DIR / "jet8_scalar_factor_diagnosis.csv"),
        "remote_probe_json": str(OUT_DIR / "remote_probe.json"),
        "remote_probe_stdout": str(OUT_DIR / "remote_probe_stdout.txt"),
        "remote_probe_stderr": str(OUT_DIR / "remote_probe_stderr.txt"),
        "conclusion_hint": (
            "jet8-only factor is not explained by the old cross-section correction alone; "
            "component/reference rows determine whether PPG12's historical jet8 SI anomaly is the target."
        ),
    }
    (OUT_DIR / "jet8_forensic_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
