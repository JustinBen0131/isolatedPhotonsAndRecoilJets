#!/usr/bin/env python3
"""Audit whether a ROOT file satisfies the PPG12 photon-yield contract.

This is intentionally an auditor, not a converter.  It checks for the object
families consumed by PPG12 efficiencytool/CalculatePhotonYield.C and reports
nearby RecoilJets diagnostic objects separately so we do not relabel a proxy as
the PPG12 photon-yield product.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Iterable


PPG12_RECO_BINS = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36]
PPG12_TRUTH_BINS = [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36, 45]

ROOT_WRAPPER = Path("scripts/root_in_analysis_env.sh")
ROOT_BIN = Path("/Users/patsfan753/Desktop/analysis/env/bin/root")


REQUIRED_OBJECTS = {
    "data_abcd": [
        "h_tight_iso_cluster_0",
        "h_tight_noniso_cluster_0",
        "h_nontight_iso_cluster_0",
        "h_nontight_noniso_cluster_0",
        "h_common_cluster_0",
    ],
    "signal_leakage": [
        "h_tight_iso_cluster_signal_0",
        "h_tight_noniso_cluster_signal_0",
        "h_nontight_iso_cluster_signal_0",
        "h_nontight_noniso_cluster_signal_0",
    ],
    "background_notmatch": [
        "h_tight_iso_cluster_background_0",
        "h_tight_noniso_cluster_background_0",
        "h_nontight_iso_cluster_background_0",
        "h_nontight_noniso_cluster_background_0",
        "h_tight_iso_cluster_notmatch_0",
        "h_tight_noniso_cluster_notmatch_0",
        "h_nontight_iso_cluster_notmatch_0",
        "h_nontight_noniso_cluster_notmatch_0",
    ],
    "efficiency_truth": [
        "eff_reco_eta_0",
        "eff_iso_eta_0",
        "eff_id_eta_0",
        "h_truth_pT_0",
        "h_truth_pT_novtx_0",
        "h_truth_pT_vertexcut_0",
        "h_truth_pT_vertexcut_mbd_cut_0",
        "h_truth_pT_vertexcut_mbd_north_cut_0",
        "h_truth_pT_vertexcut_mbd_south_cut_0",
    ],
    "response": [
        "h_pT_truth_response_0",
        "h_pT_reco_response_0",
        "h_pT_reco_fake_0",
        "h_response_full_0",
        "response_matrix_full_0",
        "h_pT_truth_half_response_0",
        "h_pT_reco_half_response_0",
        "h_response_half_0",
        "response_matrix_half_0",
    ],
}

RAW_RECOILJETS_REQUIRED_OBJECTS = {
    "data_abcd": [
        "h_tight_iso_cluster_0",
        "h_tight_noniso_cluster_0",
        "h_nontight_iso_cluster_0",
        "h_nontight_noniso_cluster_0",
        "h_common_cluster_0",
        "h_all_cluster_0",
        "h_tight_cluster_0",
    ],
    "signal_leakage": [
        "h_tight_iso_cluster_signal_0",
        "h_tight_noniso_cluster_signal_0",
        "h_nontight_iso_cluster_signal_0",
        "h_nontight_noniso_cluster_signal_0",
        "h_all_cluster_signal_0",
        "h_tight_cluster_signal_0",
    ],
    "background_notmatch": [
        "h_tight_iso_cluster_notmatch_0",
        "h_tight_noniso_cluster_notmatch_0",
        "h_nontight_iso_cluster_notmatch_0",
        "h_nontight_noniso_cluster_notmatch_0",
    ],
    "truth_response_histograms": [
        "h_truth_pT_0",
        "h_truth_pT_novtx_0",
        "h_truth_pT_vertexcut_0",
        "h_truth_pT_vertexcut_mbd_cut_0",
        "h_pT_truth_response_0",
        "h_pT_reco_response_0",
        "h_pT_reco_fake_0",
        "h_response_full_0",
        "h_pT_truth_half_response_0",
        "h_pT_reco_half_response_0",
        "h_response_half_0",
        "h_pT_truth_secondhalf_response_0",
        "h_pT_reco_secondhalf_response_0",
    ],
}

RECOILJETS_PROXY_PATTERNS = [
    r"h_pTgamma_ABCD_[ABCD]",
    r"h_isIsolated_isTight",
    r"h_notIsolated_isTight",
    r"h_isIsolated_notTight",
    r"h_notIsolated_notTight",
    r"h_sigABCD_MC",
    r"h_xJpurityLead_",
    r"h_unfoldRecoPho_pTgamma",
    r"h_unfoldTruthPho_pTgamma",
    r"h2_unfoldResponsePho_pTgamma",
]


def run_root_key_dump(root_file: Path) -> list[dict[str, str]]:
    macro = r'''
#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TClass.h>
#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TArrayD.h>
#include <TSystem.h>
#include <iostream>
#include <sstream>

std::string axis_edges(TAxis* axis)
{
  if (!axis) return "";
  std::ostringstream os;
  const int n = axis->GetNbins();
  const TArrayD* bins = axis->GetXbins();
  if (bins && bins->GetSize() == n + 1)
  {
    for (int i = 0; i <= n; ++i)
    {
      if (i) os << ",";
      os << bins->At(i);
    }
  }
  else
  {
    for (int i = 1; i <= n; ++i)
    {
      if (i > 1) os << ",";
      os << axis->GetBinLowEdge(i);
    }
    if (n > 0) os << "," << axis->GetBinUpEdge(n);
  }
  return os.str();
}

void walk(TDirectory* dir, const std::string& prefix)
{
  TIter next(dir->GetListOfKeys());
  TKey* key = nullptr;
  while ((key = static_cast<TKey*>(next())))
  {
    const std::string name = key->GetName();
    const std::string cls = key->GetClassName();
    const std::string path = prefix.empty() ? name : prefix + "/" + name;
    std::cout << "ROOTKEY\t" << cls << "\t" << path << "\n";
    TObject* obj = dir->Get(name.c_str());
    if (auto* h1 = dynamic_cast<TH1*>(obj))
    {
      std::cout << "ROOTHIST\t" << cls << "\t" << path << "\t"
                << h1->GetNbinsX() << "\t" << axis_edges(h1->GetXaxis()) << "\t";
      if (auto* h2 = dynamic_cast<TH2*>(obj))
      {
        std::cout << h2->GetNbinsY() << "\t" << axis_edges(h2->GetYaxis());
      }
      else
      {
        std::cout << "0\t";
      }
      std::cout << "\n";
    }
    TClass* cl = TClass::GetClass(cls.c_str());
    if (cl && cl->InheritsFrom(TDirectory::Class()))
    {
      TDirectory* sub = dynamic_cast<TDirectory*>(dir->Get(name.c_str()));
      if (sub) walk(sub, path);
    }
  }
}

void dump_keys(const char* filename)
{
  TFile f(filename, "READ");
  if (f.IsZombie())
  {
    std::cerr << "ERROR zombie file: " << filename << "\n";
    gSystem->Exit(2);
  }
  walk(&f, "");
}
'''
    with tempfile.TemporaryDirectory() as td:
        macro_path = Path(td) / "dump_keys.C"
        macro_path.write_text(macro)
        cmd = [
            str(ROOT_WRAPPER),
            str(ROOT_BIN),
            "-l",
            "-b",
            "-q",
            f"{macro_path}+(\"{root_file}\")",
        ]
        proc = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if proc.returncode != 0:
            raise RuntimeError(proc.stdout)

    rows: list[dict[str, str]] = []
    by_path: dict[str, dict[str, str]] = {}
    for line in proc.stdout.splitlines():
        if line.startswith("ROOTKEY\t"):
            _, cls, path = line.split("\t", 2)
            row = {"class": cls, "path": path, "basename": path.rsplit("/", 1)[-1]}
            rows.append(row)
            by_path[path] = row
        elif line.startswith("ROOTHIST\t"):
            parts = line.split("\t", 6)
            if len(parts) != 7:
                continue
            _, cls, path, nx, x_edges, ny, y_edges = parts
            row = by_path.get(path)
            if row is None:
                row = {"class": cls, "path": path, "basename": path.rsplit("/", 1)[-1]}
                rows.append(row)
                by_path[path] = row
            row.update({"nx": nx, "x_edges": x_edges, "ny": ny, "y_edges": y_edges})
    return rows


def _matching_paths(rows: Iterable[dict[str, str]], basename: str) -> list[str]:
    return sorted(row["path"] for row in rows if row["basename"] == basename)


def _matching_rows(rows: Iterable[dict[str, str]], basename: str) -> list[dict[str, str]]:
    return sorted((row for row in rows if row["basename"] == basename), key=lambda row: row["path"])


def _pattern_paths(rows: Iterable[dict[str, str]], pattern: str) -> list[str]:
    rx = re.compile(pattern)
    return sorted(row["path"] for row in rows if rx.search(row["basename"]))


def _parse_edges(text: str) -> list[float]:
    if not text:
        return []
    return [float(tok) for tok in text.split(",") if tok]


def _edges_match(found: list[float], expected: list[float], tol: float = 1.0e-6) -> bool:
    if len(found) != len(expected):
        return False
    return all(abs(a - b) <= tol for a, b in zip(found, expected))


def _expected_axes_for_required_object(name: str) -> tuple[list[float] | None, list[float] | None]:
    if name.startswith("h_response_"):
        return PPG12_RECO_BINS, PPG12_TRUTH_BINS
    if name.startswith("h_truth_pT") or name.startswith("h_pT_truth"):
        return PPG12_TRUTH_BINS, None
    if name.startswith("h_") and (
        "cluster" in name
        or name.startswith("h_pT_reco")
        or name.startswith("h_all_")
        or name.startswith("h_common_")
    ):
        return PPG12_RECO_BINS, None
    return None, None


def _binning_mismatch(name: str, row: dict[str, str]) -> str | None:
    expected_x, expected_y = _expected_axes_for_required_object(name)
    if expected_x is None:
        return None
    found_x = _parse_edges(row.get("x_edges", ""))
    if not _edges_match(found_x, expected_x):
        return f"{row['path']}: x_edges={found_x} expected={expected_x}"
    if expected_y is not None:
        found_y = _parse_edges(row.get("y_edges", ""))
        if not _edges_match(found_y, expected_y):
            return f"{row['path']}: y_edges={found_y} expected={expected_y}"
    return None


def _contract_for(mode: str, sample_role: str) -> dict[str, list[str]]:
    base = RAW_RECOILJETS_REQUIRED_OBJECTS if mode == "raw-recoiljets" else REQUIRED_OBJECTS
    if sample_role == "all":
        return base
    if sample_role == "data":
        return {"data_abcd": base["data_abcd"]}
    if sample_role == "signal-mc":
        keep = ["data_abcd", "signal_leakage", "truth_response_histograms"]
        return {key: base[key] for key in keep if key in base}
    if sample_role == "inclusive-mc":
        keep = ["data_abcd", "signal_leakage", "background_notmatch", "truth_response_histograms"]
        return {key: base[key] for key in keep if key in base}
    raise ValueError(f"unknown sample_role: {sample_role}")


def build_report(root_file: Path, mode: str, sample_role: str) -> dict:
    rows = run_root_key_dump(root_file)
    contract = _contract_for(mode, sample_role)
    required = {}
    missing = []
    binning_mismatches = []
    for family, names in contract.items():
        family_rows = {}
        for name in names:
            matching_rows = _matching_rows(rows, name)
            matches = [row["path"] for row in matching_rows]
            family_rows[name] = matches
            if not matches:
                missing.append(name)
            else:
                for row in matching_rows:
                    mismatch = _binning_mismatch(name, row)
                    if mismatch:
                        binning_mismatches.append(mismatch)
        required[family] = family_rows

    proxies = {pat: _pattern_paths(rows, pat)[:40] for pat in RECOILJETS_PROXY_PATTERNS}
    verdict = "pass" if not missing and not binning_mismatches else "fail"

    return {
        "root_file": str(root_file),
        "mode": mode,
        "sample_role": sample_role,
        "object_count": len(rows),
        "ppg12_reco_bins": PPG12_RECO_BINS,
        "ppg12_truth_bins": PPG12_TRUTH_BINS,
        "required": required,
        "missing_required": sorted(missing),
        "missing_required_count": len(missing),
        "binning_mismatches": sorted(binning_mismatches),
        "binning_mismatch_count": len(binning_mismatches),
        "recoiljets_proxy_patterns": proxies,
        "verdict": verdict,
        "interpretation": (
            "This file satisfies the PPG12 Photon_final/CalculatePhotonYield object contract."
            if verdict == "pass"
            else (
                "This file does not satisfy the raw RecoilJets PPG12_PHOTON_YIELD_V1 contract."
                if mode == "raw-recoiljets"
                else "This file does not satisfy the PPG12 Photon_final/CalculatePhotonYield object contract; "
                     "RecoilJets proxy ABCD/response objects must not be relabeled as Fig.29 photon-yield inputs."
            )
        ),
    }


def write_markdown(report: dict, path: Path) -> None:
    lines = [
        "# THE-76 PPG12 Photon-Yield Contract Audit",
        "",
        f"- ROOT file: `{report['root_file']}`",
        f"- Mode: `{report['mode']}`",
        f"- Sample role: `{report['sample_role']}`",
        f"- Verdict: **{report['verdict'].upper()}**",
        f"- Required objects missing: `{report['missing_required_count']}`",
        f"- Required-object binning mismatches: `{report['binning_mismatch_count']}`",
        f"- PPG12 reco bins: `{report['ppg12_reco_bins']}`",
        f"- PPG12 truth bins: `{report['ppg12_truth_bins']}`",
        "",
        "## Required PPG12 Objects",
        "",
        "| family | object | status | matching paths |",
        "| --- | --- | --- | --- |",
    ]
    for family, objects in report["required"].items():
        for name, matches in objects.items():
            status = "present" if matches else "missing"
            match_text = "<br>".join(f"`{m}`" for m in matches[:5]) if matches else ""
            lines.append(f"| {family} | `{name}` | {status} | {match_text} |")

    if report["binning_mismatches"]:
        lines.extend([
            "",
            "## Binning Mismatches",
            "",
        ])
        for mismatch in report["binning_mismatches"]:
            lines.append(f"- `{mismatch}`")

    lines.extend([
        "",
        "## RecoilJets Proxy Objects Found",
        "",
        "These are useful diagnostics or raw ingredients, but they are not the PPG12",
        "`Photon_final` contract consumed by `CalculatePhotonYield.C` unless an",
        "explicit converter/provenance layer maps them with the correct bins and",
        "semantics.",
        "",
        "| pattern | example matching paths |",
        "| --- | --- |",
    ])
    for pattern, matches in report["recoiljets_proxy_patterns"].items():
        match_text = "<br>".join(f"`{m}`" for m in matches[:8]) if matches else ""
        lines.append(f"| `{pattern}` | {match_text} |")

    lines.extend([
        "",
        "## Interpretation",
        "",
        report["interpretation"],
        "",
    ])
    path.write_text("\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path, help="ROOT file to audit")
    parser.add_argument("--out-dir", required=True, type=Path, help="Directory for JSON/Markdown reports")
    parser.add_argument(
        "--mode",
        choices=("photon-final", "raw-recoiljets"),
        default="photon-final",
        help="photon-final is the strict PPG12 CalculatePhotonYield contract; raw-recoiljets checks the pre-export RecoilJets contract.",
    )
    parser.add_argument(
        "--sample-role",
        choices=("all", "data", "signal-mc", "inclusive-mc"),
        default="all",
        help="Limit required object families to what the sample type can physically contain.",
    )
    args = parser.parse_args()

    root_file = args.root.resolve()
    if not root_file.exists():
        raise FileNotFoundError(root_file)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    report = build_report(root_file, args.mode, args.sample_role)
    json_path = args.out_dir / "the76_ppg12_photon_yield_contract_audit.json"
    md_path = args.out_dir / "the76_ppg12_photon_yield_contract_audit.md"
    json_path.write_text(json.dumps(report, indent=2, sort_keys=True))
    write_markdown(report, md_path)

    print(f"wrote {json_path}")
    print(f"wrote {md_path}")
    print(f"verdict={report['verdict']} missing={report['missing_required_count']}")
    return 0 if report["verdict"] == "pass" else 3


if __name__ == "__main__":
    raise SystemExit(main())
