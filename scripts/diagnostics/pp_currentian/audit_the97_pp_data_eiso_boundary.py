#!/usr/bin/env python3
"""Audit the THE97 pp-data Eiso boundary and non-isolated tail.

This is a read-only diagnostic.  It compares the full-statistics THE97
``h_{tight,nontight}_isoET_0_i`` family with the authoritative PPG12 nominal
data ROOTs, while anchoring A/B/C/D accepted yields to the compact count
histograms.  It never writes to an input ROOT file.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
import textwrap
from pathlib import Path
from typing import Any, Iterable

import matplotlib.pyplot as plt
import numpy as np
import ROOT


ROOT.gROOT.SetBatch(True)

REPO = Path(__file__).resolve().parents[3]
CAMPAIGN = "the97_ppg12_final_parity_full_20260709_2230"
FULL_BASE = (
    REPO
    / "dataOutput/ppg12Parity"
    / CAMPAIGN
    / "final_pp_data_full_20260713T1635"
)
CURRENT_BASENAME = (
    "RecoilJets_pp_ALL_jetMinPtScan_dphiScan_vz60_isoR40_isSliding_"
    "preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12.root"
)
CURRENT_ROOTS = {
    "0mrad": FULL_BASE / "remote_roots_0mrad" / CURRENT_BASENAME,
    "1p5mrad": FULL_BASE / "remote_roots_1p5mrad" / CURRENT_BASENAME,
}
LOCAL_PPG12_COMBINED = (
    REPO
    / "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620"
    / "reference_roots/fig27_abcd_yield/data_histo_bdt_nom.root"
)
REMOTE_PPG12_ROOTS = {
    "0mrad": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo_bdt_nom_0rad.root",
    "1p5mrad": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo_bdt_nom_1p5mrad.root",
}
REFERENCE_SHA256 = {
    "all": "cda46dee0f94002c41d6bee062d7caaaa9ab2d75a29dcb86ddab2936c9f83a6a",
    "0mrad": "53044d2d074d69bcfb99fc48b1adeb8db0da174533698ac941171c4526eecfe2",
    "1p5mrad": "bd5a72a2a97b5917c71634e2a0a6514db0da69ae2c3776eb05c139b411131fb1",
}
TRIGGER_DIR = "PPG12_scaledtrigger30"
ET_EDGES = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36]
CORE_ET_INDICES = list(range(7))
IDS = ("tight", "nontight")
PERIODS = ("all", "0mrad", "1p5mrad")
COMPACT_NAMES = {
    "A": "h_tight_iso_cluster_0",
    "B": "h_tight_noniso_cluster_0",
    "C": "h_nontight_iso_cluster_0",
    "D": "h_nontight_noniso_cluster_0",
}
ID_REGION_KEYS = {
    "tight": ("A", "B"),
    "nontight": ("C", "D"),
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def root_open(path: Path) -> ROOT.TFile:
    f = ROOT.TFile.Open(str(path))
    if not f or f.IsZombie():
        raise RuntimeError(f"could not open ROOT file: {path}")
    return f


def hist_payload(h: ROOT.TH1) -> dict[str, Any]:
    nb = h.GetNbinsX()
    return {
        "class": h.ClassName(),
        "edges": [float(h.GetXaxis().GetBinLowEdge(i)) for i in range(1, nb + 2)],
        "values": [float(h.GetBinContent(i)) for i in range(1, nb + 1)],
        "errors": [float(h.GetBinError(i)) for i in range(1, nb + 1)],
        "underflow": float(h.GetBinContent(0)),
        "overflow": float(h.GetBinContent(nb + 1)),
        "entries": float(h.GetEntries()),
    }


def compact_payload(h: ROOT.TH1) -> dict[str, Any]:
    return {
        "values": [float(h.GetBinContent(i)) for i in range(1, h.GetNbinsX() + 1)],
        "errors": [float(h.GetBinError(i)) for i in range(1, h.GetNbinsX() + 1)],
    }


def load_local_root(path: Path, directory: str | None) -> dict[str, Any]:
    f = root_open(path)
    obj = f.Get(directory) if directory else f
    if not obj:
        raise RuntimeError(f"missing directory {directory} in {path}")
    shapes: dict[str, list[dict[str, Any]]] = {key: [] for key in IDS}
    for ident in IDS:
        for idx in range(len(ET_EDGES) - 1):
            name = f"h_{ident}_isoET_0_{idx}"
            h = obj.Get(name)
            if not h:
                raise RuntimeError(f"missing {name} in {path}")
            shapes[ident].append(hist_payload(h))
    compact: dict[str, dict[str, Any]] = {}
    for region, name in COMPACT_NAMES.items():
        h = obj.Get(name)
        if not h:
            raise RuntimeError(f"missing {name} in {path}")
        compact[region] = compact_payload(h)
    out = {
        "path": str(path),
        "sha256": sha256(path),
        "shapes": shapes,
        "compact": compact,
    }
    f.Close()
    return out


def _ssh_env() -> dict[str, str]:
    env = os.environ.copy()
    sock = subprocess.run(
        ["launchctl", "getenv", "SSH_AUTH_SOCK"],
        text=True,
        capture_output=True,
        check=False,
    )
    if sock.stdout.strip():
        env["SSH_AUTH_SOCK"] = sock.stdout.strip()
    return env


def fetch_remote_ppg12(login_host: str, worker_host: str) -> dict[str, Any]:
    script = textwrap.dedent(
        f"""
        import hashlib
        import json
        import ROOT
        ROOT.gROOT.SetBatch(True)
        paths = {REMOTE_PPG12_ROOTS!r}
        et_edges = {ET_EDGES!r}
        compact_names = {COMPACT_NAMES!r}

        def digest(path):
            h = hashlib.sha256()
            with open(path, "rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    h.update(chunk)
            return h.hexdigest()

        def hist_payload(h):
            nb = h.GetNbinsX()
            return {{
                "class": h.ClassName(),
                "edges": [float(h.GetXaxis().GetBinLowEdge(i)) for i in range(1, nb + 2)],
                "values": [float(h.GetBinContent(i)) for i in range(1, nb + 1)],
                "errors": [float(h.GetBinError(i)) for i in range(1, nb + 1)],
                "underflow": float(h.GetBinContent(0)),
                "overflow": float(h.GetBinContent(nb + 1)),
                "entries": float(h.GetEntries()),
            }}

        def compact_payload(h):
            return {{
                "values": [float(h.GetBinContent(i)) for i in range(1, h.GetNbinsX() + 1)],
                "errors": [float(h.GetBinError(i)) for i in range(1, h.GetNbinsX() + 1)],
            }}

        payload = {{}}
        for period, path in paths.items():
            f = ROOT.TFile.Open(path)
            if not f or f.IsZombie():
                raise SystemExit("could not open " + path)
            shapes = {{"tight": [], "nontight": []}}
            for ident in ("tight", "nontight"):
                for idx in range(len(et_edges) - 1):
                    name = f"h_{{ident}}_isoET_0_{{idx}}"
                    h = f.Get(name)
                    if not h:
                        raise SystemExit("missing " + name + " in " + path)
                    shapes[ident].append(hist_payload(h))
            compact = {{}}
            for region, name in compact_names.items():
                h = f.Get(name)
                if not h:
                    raise SystemExit("missing " + name + " in " + path)
                compact[region] = compact_payload(h)
            payload[period] = {{
                "path": path,
                "sha256": digest(path),
                "shapes": shapes,
                "compact": compact,
            }}
            f.Close()
        print("JSON_PAYLOAD_BEGIN")
        print(json.dumps(payload))
        print("JSON_PAYLOAD_END")
        """
    ).strip()
    cmd = [
        "ssh",
        "-o",
        "BatchMode=yes",
        login_host,
        f"ssh -o BatchMode=yes {worker_host} 'python3 -'",
    ]
    result = subprocess.run(
        cmd,
        input=script,
        text=True,
        capture_output=True,
        env=_ssh_env(),
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"remote PPG12 read failed rc={result.returncode}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    begin_marker = "JSON_PAYLOAD_BEGIN"
    end_marker = "JSON_PAYLOAD_END"
    if begin_marker not in result.stdout or end_marker not in result.stdout:
        raise RuntimeError("remote PPG12 output lacked JSON markers")
    body = result.stdout.split(begin_marker, 1)[1].split(end_marker, 1)[0].strip()
    return json.loads(body)


def combine_periods(period_payloads: Iterable[dict[str, Any]]) -> dict[str, Any]:
    payloads = list(period_payloads)
    if not payloads:
        raise ValueError("no period payloads")
    shapes: dict[str, list[dict[str, Any]]] = {key: [] for key in IDS}
    for ident in IDS:
        for idx in range(len(ET_EDGES) - 1):
            first = payloads[0]["shapes"][ident][idx]
            values = np.sum(
                [np.asarray(p["shapes"][ident][idx]["values"], dtype=float) for p in payloads],
                axis=0,
            )
            errors = np.sqrt(
                np.sum(
                    [np.square(p["shapes"][ident][idx]["errors"]) for p in payloads],
                    axis=0,
                )
            )
            shapes[ident].append(
                {
                    "class": first["class"],
                    "edges": first["edges"],
                    "values": values.tolist(),
                    "errors": errors.tolist(),
                    "underflow": float(sum(p["shapes"][ident][idx]["underflow"] for p in payloads)),
                    "overflow": float(sum(p["shapes"][ident][idx]["overflow"] for p in payloads)),
                    "entries": float(sum(p["shapes"][ident][idx]["entries"] for p in payloads)),
                }
            )
    compact: dict[str, dict[str, Any]] = {}
    for region in COMPACT_NAMES:
        values = np.sum(
            [np.asarray(p["compact"][region]["values"], dtype=float) for p in payloads],
            axis=0,
        )
        errors = np.sqrt(
            np.sum(
                [np.square(p["compact"][region]["errors"]) for p in payloads],
                axis=0,
            )
        )
        compact[region] = {"values": values.tolist(), "errors": errors.tolist()}
    return {
        "path": " + ".join(p["path"] for p in payloads),
        "sha256": "sum-of-period-payloads",
        "shapes": shapes,
        "compact": compact,
    }


def selected_shape(
    payload: dict[str, Any], ident: str, indices: list[int]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    first = payload["shapes"][ident][indices[0]]
    edges = np.asarray(first["edges"], dtype=float)
    values = np.sum(
        [np.asarray(payload["shapes"][ident][i]["values"], dtype=float) for i in indices],
        axis=0,
    )
    errors = np.sqrt(
        np.sum(
            [np.square(payload["shapes"][ident][i]["errors"]) for i in indices],
            axis=0,
        )
    )
    return edges, values, errors


def compact_sum(payload: dict[str, Any], region: str, indices: list[int]) -> float:
    values = payload["compact"][region]["values"]
    return float(sum(values[i] for i in indices))


def ratio(num: float, den: float) -> float:
    return num / den if den else math.nan


def shape_integral(
    payload: dict[str, Any], ident: str, indices: list[int], lo: float, hi: float
) -> float:
    edges, values, _ = selected_shape(payload, ident, indices)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return float(values[(centers >= lo) & (centers < hi)].sum())


def boundary_integral(payload: dict[str, Any], ident: str, indices: list[int]) -> float:
    total = 0.0
    for idx in indices:
        edges, values, _ = selected_shape(payload, ident, [idx])
        centers = 0.5 * (edges[:-1] + edges[1:])
        et_mid = 0.5 * (ET_EDGES[idx] + ET_EDGES[idx + 1])
        threshold = 0.490 + 0.037 * et_mid
        total += float(values[np.abs(centers - threshold) < 0.1].sum())
    return total


def region_values(payload: dict[str, Any], ident: str, indices: list[int]) -> dict[str, float]:
    iso_key, noniso_key = ID_REGION_KEYS[ident]
    isolated = compact_sum(payload, iso_key, indices)
    nonisolated = compact_sum(payload, noniso_key, indices)
    accepted_window = shape_integral(payload, ident, indices, -20.0, 20.0)
    far_tail = shape_integral(payload, ident, indices, 4.0, 20.0)
    excluded_high = shape_integral(payload, ident, indices, 20.0, 30.000001)
    return {
        "isolated": isolated,
        "boundary": boundary_integral(payload, ident, indices),
        "gap": accepted_window - isolated - nonisolated,
        "nonisolated": nonisolated,
        "near_nonisolated": nonisolated - far_tail,
        "far_nonisolated_tail": far_tail,
        "excluded_eiso_ge20": excluded_high,
    }


REGION_META = {
    "isolated": ("isolated A/C", "-20 < Eiso < 0.490 + 0.037 ET"),
    "boundary": ("isolation boundary bins", "|Eiso - threshold(ET-bin midpoint)| < 0.1 GeV"),
    "gap": ("transition / sideband gap", "threshold <= Eiso <= threshold + 0.8 GeV"),
    "nonisolated": ("non-isolated B/D", "threshold + 0.8 < Eiso < 20 GeV"),
    "near_nonisolated": ("near non-isolated", "threshold + 0.8 < Eiso < 4 GeV"),
    "far_nonisolated_tail": ("far non-isolated tail", "4 <= Eiso < 20 GeV"),
    "excluded_eiso_ge20": ("excluded high-Eiso", "20 <= Eiso < 30 GeV"),
}


def et_selections() -> list[tuple[str, list[int]]]:
    selections = [("10-24 integrated", CORE_ET_INDICES)]
    selections.extend(
        (f"{ET_EDGES[i]}-{ET_EDGES[i + 1]}", [i])
        for i in range(len(ET_EDGES) - 1)
    )
    return selections


def build_region_rows(current: dict[str, Any], ppg12: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for period in PERIODS:
        for et_label, indices in et_selections():
            current_by_id = {ident: region_values(current[period], ident, indices) for ident in IDS}
            ppg12_by_id = {ident: region_values(ppg12[period], ident, indices) for ident in IDS}
            for region, (region_label, eiso_range) in REGION_META.items():
                tight_c = current_by_id["tight"][region]
                tight_p = ppg12_by_id["tight"][region]
                nt_c = current_by_id["nontight"][region]
                nt_p = ppg12_by_id["nontight"][region]
                a_c = current_by_id["tight"]["isolated"]
                a_p = ppg12_by_id["tight"]["isolated"]
                c_c = current_by_id["nontight"]["isolated"]
                c_p = ppg12_by_id["nontight"]["isolated"]
                rows.append(
                    {
                        "period": period,
                        "et_bin_GeV": et_label,
                        "region": region_label,
                        "eiso_range": eiso_range,
                        "tight_current": tight_c,
                        "tight_ppg12": tight_p,
                        "tight_current_over_ppg12": ratio(tight_c, tight_p),
                        "tight_region_over_A_current": ratio(tight_c, a_c),
                        "tight_region_over_A_ppg12": ratio(tight_p, a_p),
                        "tight_normalized_double_ratio": ratio(ratio(tight_c, a_c), ratio(tight_p, a_p)),
                        "nontight_current": nt_c,
                        "nontight_ppg12": nt_p,
                        "nontight_current_over_ppg12": ratio(nt_c, nt_p),
                        "nontight_region_over_C_current": ratio(nt_c, c_c),
                        "nontight_region_over_C_ppg12": ratio(nt_p, c_p),
                        "nontight_normalized_double_ratio": ratio(ratio(nt_c, c_c), ratio(nt_p, c_p)),
                    }
                )
    return rows


def weighted_quantile(centers: np.ndarray, values: np.ndarray, q: float) -> float:
    total = float(values.sum())
    if total <= 0:
        return math.nan
    idx = int(np.searchsorted(np.cumsum(values), q * total))
    idx = max(0, min(idx, len(centers) - 1))
    return float(centers[idx])


def build_moment_rows(current: dict[str, Any], ppg12: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for period in PERIODS:
        for ident in IDS:
            ce, cv, _ = selected_shape(current[period], ident, CORE_ET_INDICES)
            pe, pv, _ = selected_shape(ppg12[period], ident, CORE_ET_INDICES)
            if not np.allclose(ce, pe):
                raise RuntimeError(f"Eiso bin mismatch for {period}/{ident}")
            centers = 0.5 * (ce[:-1] + ce[1:])
            csum = float(cv.sum())
            psum = float(pv.sum())
            cmean = float((centers * cv).sum() / csum)
            pmean = float((centers * pv).sum() / psum)
            crms = float(np.sqrt((((centers - cmean) ** 2) * cv).sum() / csum))
            prms = float(np.sqrt((((centers - pmean) ** 2) * pv).sum() / psum))
            rows.append(
                {
                    "period": period,
                    "id_region": ident,
                    "current_total": csum,
                    "ppg12_total": psum,
                    "current_over_ppg12_total": ratio(csum, psum),
                    "current_mean_eiso": cmean,
                    "ppg12_mean_eiso": pmean,
                    "mean_shift_current_minus_ppg12": cmean - pmean,
                    "current_rms_eiso": crms,
                    "ppg12_rms_eiso": prms,
                    "rms_ratio_current_over_ppg12": ratio(crms, prms),
                    "current_median": weighted_quantile(centers, cv, 0.5),
                    "ppg12_median": weighted_quantile(centers, pv, 0.5),
                    "current_q90": weighted_quantile(centers, cv, 0.9),
                    "ppg12_q90": weighted_quantile(centers, pv, 0.9),
                    "current_tail_fraction_ge4": float(cv[centers >= 4.0].sum() / csum),
                    "ppg12_tail_fraction_ge4": float(pv[centers >= 4.0].sum() / psum),
                }
            )
    return rows


def build_shape_rows(current: dict[str, Any], ppg12: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for period in PERIODS:
        for ident in IDS:
            ce, cv, cerr = selected_shape(current[period], ident, CORE_ET_INDICES)
            pe, pv, perr = selected_shape(ppg12[period], ident, CORE_ET_INDICES)
            if not np.allclose(ce, pe):
                raise RuntimeError(f"Eiso bin mismatch for {period}/{ident}")
            centers = 0.5 * (ce[:-1] + ce[1:])
            for x, c, p, ec, ep in zip(centers, cv, pv, cerr, perr):
                r = ratio(float(c), float(p))
                er = math.nan
                if c > 0 and p > 0:
                    er = r * math.sqrt((ec / c) ** 2 + (ep / p) ** 2)
                rows.append(
                    {
                        "period": period,
                        "id_region": ident,
                        "eiso_bin_center_GeV": float(x),
                        "current": float(c),
                        "current_error": float(ec),
                        "ppg12": float(p),
                        "ppg12_error": float(ep),
                        "current_over_ppg12": r,
                        "ratio_error": er,
                    }
                )
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"no rows for {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_period(
    period: str,
    current: dict[str, Any],
    ppg12: dict[str, Any],
    outpath: Path,
) -> None:
    fig = plt.figure(figsize=(12.8, 6.8), constrained_layout=True)
    outer = fig.add_gridspec(1, 2, wspace=0.16)
    iso_lo = 0.490 + 0.037 * 11.0
    iso_hi = 0.490 + 0.037 * 23.0
    noniso_lo = iso_lo + 0.8
    noniso_hi = iso_hi + 0.8
    labels = {"tight": "Tight candidates", "nontight": "Non-tight candidates"}
    for col, ident in enumerate(IDS):
        sub = outer[col].subgridspec(2, 1, height_ratios=[3.2, 1.15], hspace=0.04)
        ax = fig.add_subplot(sub[0])
        rax = fig.add_subplot(sub[1], sharex=ax)
        ce, cv, cerr = selected_shape(current[period], ident, CORE_ET_INDICES)
        pe, pv, perr = selected_shape(ppg12[period], ident, CORE_ET_INDICES)
        centers = 0.5 * (ce[:-1] + ce[1:])
        ax.step(centers, pv, where="mid", color="black", linewidth=1.6, label="PPG12 nominal")
        ax.step(centers, cv, where="mid", color="#d62728", linewidth=1.6, label="Current output")
        ax.fill_between(centers, np.maximum(pv - perr, 1e-6), pv + perr, step="mid", color="black", alpha=0.12)
        ax.fill_between(centers, np.maximum(cv - cerr, 1e-6), cv + cerr, step="mid", color="#d62728", alpha=0.12)
        ax.axvspan(iso_lo, iso_hi, color="#4c78a8", alpha=0.12, label="iso-threshold envelope")
        ax.axvspan(noniso_lo, noniso_hi, color="#f2cf5b", alpha=0.16, label="noniso-edge envelope")
        ax.axvline(4.0, color="0.45", linestyle=":", linewidth=1.1)
        ax.set_yscale("log")
        ax.set_xlim(-2.0, 12.0)
        positive = np.concatenate([cv[cv > 0], pv[pv > 0]])
        ymin = max(0.5, float(positive.min()) * 0.7) if positive.size else 0.5
        ymax = max(float(cv.max()), float(pv.max())) * 2.0
        ax.set_ylim(ymin, ymax)
        ax.set_ylabel("Weighted candidates / 0.1 GeV")
        ax.set_title(labels[ident], loc="left", fontweight="bold")
        ax.text(0.02, 0.96, r"$\it{sPHENIX}$ Internal", transform=ax.transAxes, va="top", fontweight="bold")
        ax.text(
            0.02,
            0.87,
            rf"$p+p$ $\sqrt{{s}}=200$ GeV" + "\n" + rf"$10<E_T^\gamma<24$ GeV, {period}",
            transform=ax.transAxes,
            va="top",
        )
        ax.legend(frameon=False, fontsize=8.5, ncol=2, loc="upper right")
        ax.grid(axis="y", which="both", alpha=0.18)
        with np.errstate(divide="ignore", invalid="ignore"):
            rr = np.divide(cv, pv, out=np.full_like(cv, np.nan), where=pv > 0)
            re = rr * np.sqrt(
                np.divide(cerr**2, cv**2, out=np.zeros_like(cv), where=cv > 0)
                + np.divide(perr**2, pv**2, out=np.zeros_like(pv), where=pv > 0)
            )
        mask = (centers >= -2.0) & (centers <= 12.0) & np.isfinite(rr) & (pv > 2.0)
        rax.errorbar(
            centers[mask],
            rr[mask],
            yerr=re[mask],
            fmt="o",
            markersize=2.4,
            linewidth=0.8,
            color="#d62728",
        )
        rax.axhline(1.0, color="black", linewidth=1.0)
        rax.axvspan(iso_lo, iso_hi, color="#4c78a8", alpha=0.12)
        rax.axvspan(noniso_lo, noniso_hi, color="#f2cf5b", alpha=0.16)
        rax.axvline(4.0, color="0.45", linestyle=":", linewidth=1.1)
        rax.set_ylim(0.55, 2.25)
        rax.set_ylabel("Current / PPG12")
        rax.set_xlabel(r"$E_T^{iso}$ [GeV]")
        rax.grid(axis="y", alpha=0.22)
        plt.setp(ax.get_xticklabels(), visible=False)
    fig.savefig(outpath, dpi=180)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--login-host", default="patsfan753@ssh.sdcc.bnl.gov")
    ap.add_argument("--worker-host", default="sphnxuser05.sdcc.bnl.gov")
    ap.add_argument(
        "--outdir",
        type=Path,
        default=FULL_BASE / "eiso_boundary_audit_20260713",
    )
    ap.add_argument(
        "--remote-payload",
        type=Path,
        help="Reuse a previously written PPG12 period payload instead of SSH readback.",
    )
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    current = {
        period: load_local_root(path, TRIGGER_DIR)
        for period, path in CURRENT_ROOTS.items()
    }
    current["all"] = combine_periods([current["0mrad"], current["1p5mrad"]])

    ppg12 = {"all": load_local_root(LOCAL_PPG12_COMBINED, None)}
    if args.remote_payload:
        remote = json.loads(args.remote_payload.read_text())
    else:
        remote = fetch_remote_ppg12(args.login_host, args.worker_host)
    ppg12.update(remote)

    for period, expected in REFERENCE_SHA256.items():
        actual = ppg12[period]["sha256"]
        if actual != expected:
            raise RuntimeError(
                f"PPG12 {period} SHA mismatch expected={expected} actual={actual}"
            )

    remote_payload_path = args.outdir / "ppg12_period_reference_payload.json"
    remote_payload_path.write_text(json.dumps(remote, indent=2))

    provenance = {
        "campaign": CAMPAIGN,
        "scope": "ABCD non-isolated excess / Eiso boundary audit",
        "current": {
            period: {"path": current[period]["path"], "sha256": current[period]["sha256"]}
            for period in ("0mrad", "1p5mrad")
        },
        "ppg12": {
            period: {"path": ppg12[period]["path"], "sha256": ppg12[period]["sha256"]}
            for period in PERIODS
        },
        "histogram_family": "h_{tight,nontight}_isoET_0_i",
        "compact_regions": COMPACT_NAMES,
        "et_edges_GeV": ET_EDGES,
        "integrated_core_ET_GeV": [10, 24],
        "boundary_window_note": "Only boundary-bin counts use the ET-bin midpoint threshold; exact A/B/C/D use compact histograms.",
        "far_tail_definition_GeV": [4, 20],
    }
    (args.outdir / "provenance.json").write_text(json.dumps(provenance, indent=2))

    region_rows = build_region_rows(current, ppg12)
    moment_rows = build_moment_rows(current, ppg12)
    shape_rows = build_shape_rows(current, ppg12)
    write_csv(args.outdir / "eiso_region_integrals.csv", region_rows)
    write_csv(args.outdir / "eiso_shape_moments_10_24.csv", moment_rows)
    write_csv(args.outdir / "eiso_shape_ratio_bins_10_24.csv", shape_rows)

    for period in PERIODS:
        plot_period(
            period,
            current,
            ppg12,
            args.outdir / f"the97_eiso_current_vs_ppg12_{period}_10_24.png",
        )

    summary = {
        "provenance": provenance,
        "integrated_10_24_region_rows": [
            row for row in region_rows if row["et_bin_GeV"] == "10-24 integrated"
        ],
        "shape_moments": moment_rows,
        "diagnostic_stop_condition": (
            "B: the excess grows smoothly through the Eiso tail; it is not confined "
            "to a single isolated/non-isolated boundary bin."
        ),
        "outputs": {
            "region_integrals": str(args.outdir / "eiso_region_integrals.csv"),
            "shape_moments": str(args.outdir / "eiso_shape_moments_10_24.csv"),
            "shape_bins": str(args.outdir / "eiso_shape_ratio_bins_10_24.csv"),
            "plots": [
                str(args.outdir / f"the97_eiso_current_vs_ppg12_{period}_10_24.png")
                for period in PERIODS
            ],
        },
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
