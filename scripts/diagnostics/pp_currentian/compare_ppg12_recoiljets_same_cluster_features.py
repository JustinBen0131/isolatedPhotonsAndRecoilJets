#!/usr/bin/env python3
"""Compare PPG12 slimtree clusters against RecoilJets pp photon-ID rows.

This is a read-only THE-76 forensic helper.  It matches signal clusters by
eventnumber + truth track id + nearest eta/phi, then reports where the BDT
classification inputs first diverge.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import ROOT  # type: ignore


ROOT.gROOT.SetBatch(True)

RECO_BINS = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36]
FEATURES = [
    ("cluster_Et_score_input", "cluster_Et"),
    ("cluster_weta_cogx", "cluster_weta_cogx"),
    ("cluster_wphi_cogx", "cluster_wphi_cogx"),
    ("vertexz", "vertexz"),
    ("cluster_Eta", "cluster_Eta"),
    ("e11_over_e33", "e11_over_e33"),
    ("cluster_et1", "cluster_et1"),
    ("cluster_et2", "cluster_et2"),
    ("cluster_et3", "cluster_et3"),
    ("cluster_et4", "cluster_et4"),
    ("e32_over_e35", "e32_over_e35"),
]


def find_bin(x: float) -> int:
    for i in range(len(RECO_BINS) - 1):
        if RECO_BINS[i] < x < RECO_BINS[i + 1]:
            return i
    return -1


def bin_label(i: int) -> str:
    return f"{RECO_BINS[i]}-{RECO_BINS[i + 1]}"


def finite(x: Any) -> bool:
    try:
        return math.isfinite(float(x))
    except Exception:
        return False


def open_interval(x: float, lo: float, hi: float) -> bool:
    return finite(x) and lo < x < hi


def safe_div(n: float, d: float) -> float:
    return n / d if d else float("nan")


def dphi(a: float, b: float) -> float:
    x = a - b
    while x > math.pi:
        x -= 2.0 * math.pi
    while x <= -math.pi:
        x += 2.0 * math.pi
    return x


def selected_ppg12_model(smeared_et: float) -> str:
    # config_bdt_nom.yaml: fallback base_E, ET bins [8, 15, 35] both base_v3E.
    if 8.0 <= smeared_et < 35.0:
        return "base_v3E"
    return "base_E"


def ppg12_thresholds(et: float) -> tuple[float, float, float]:
    tight_min = 0.8333333333333334 - 0.003333333333333336 * et
    nt_min = 0.7333333333333333 - 0.01333333333333333 * et
    nt_max = 0.6666666666666666 + 0.003333333333333336 * et
    return tight_min, nt_min, nt_max


def classify(row: dict[str, float], score: float, et_for_cuts: float) -> tuple[int, dict[str, bool]]:
    """Return PPG12 tag: 0 preselection fail, 1 tight, 2 nontight, 3 neither."""
    e11e33 = row["e11_over_e33"]
    e32e35 = row["e32_over_e35"]
    wr = row["cluster_wphi_cogx"] / row["cluster_weta_cogx"] if row["cluster_weta_cogx"] != 0 else float("nan")
    common = (
        open_interval(row["cluster_prob"], 0.0, 1.0)
        and open_interval(e11e33, 0.0, 0.98)
        and finite(wr)
        and wr > 0.0
        and finite(row["cluster_weta_cogx"])
        and row["cluster_weta_cogx"] < 2.0
        and row["npb_score"] > 0.5
    )
    flags: dict[str, bool] = {"common": common}
    if not common:
        return 0, flags

    tight_min, nt_min, nt_max = ppg12_thresholds(et_for_cuts)
    tight_prob = open_interval(row["cluster_prob"], 0.0, 1.0)
    tight_weta = open_interval(row["cluster_weta_cogx"], 0.0, 1.0)
    tight_wphi = open_interval(row["cluster_wphi_cogx"], 0.0, 1.0)
    tight_et1 = open_interval(row["cluster_et1"], 0.5, 1.0)
    tight_et2 = open_interval(row["cluster_et2"], 0.0, 1.0)
    tight_et3 = open_interval(row["cluster_et3"], 0.0, 1.0)
    tight_et4 = open_interval(row["cluster_et4"], 0.0, 1.0)
    tight_e11e33 = open_interval(e11e33, 0.0, 1.0)
    tight_e32e35 = open_interval(e32e35, 0.8, 1.0)
    tight_bdt = finite(score) and score > tight_min and score < 1.0

    flags.update(
        tight_prob=tight_prob,
        tight_weta=tight_weta,
        tight_wphi=tight_wphi,
        tight_et1=tight_et1,
        tight_et2=tight_et2,
        tight_et3=tight_et3,
        tight_et4=tight_et4,
        tight_e11e33=tight_e11e33,
        tight_e32e35=tight_e32e35,
        tight_bdt=tight_bdt,
        nt_bdt=finite(score) and score > nt_min and score < nt_max,
    )
    tight = all(
        flags[k]
        for k in (
            "tight_prob",
            "tight_weta",
            "tight_wphi",
            "tight_et1",
            "tight_et2",
            "tight_et3",
            "tight_et4",
            "tight_e11e33",
            "tight_e32e35",
            "tight_bdt",
        )
    )
    if tight:
        return 1, flags

    nt_shape = (
        open_interval(row["cluster_prob"], 0.0, 1.0)
        and open_interval(row["cluster_weta_cogx"], 0.0, 1.0)
        and open_interval(row["cluster_wphi_cogx"], 0.0, 1.0)
        and open_interval(row["cluster_et1"], 0.6, 1.0)
        and open_interval(row["cluster_et4"], 0.0, 1.0)
        and open_interval(e11e33, 0.0, 1.0)
        and open_interval(e32e35, 0.8, 1.0)
    )
    nfail = 0
    if not tight_weta:
        nfail += 1
    if not tight_prob:
        nfail += 1
    if not tight_bdt:
        nfail += 1
    flags["nt_shape"] = nt_shape
    flags["nfail_any"] = nfail > 0
    if nt_shape and flags["nt_bdt"] and nfail > 0:
        return 2, flags
    return 3, flags


def tower_masked(mask: Any, ieta_value: float, iphi_value: float) -> bool:
    if not mask:
        return False
    ieta = int(ieta_value)
    iphi = int(iphi_value)
    if ieta < 0 or ieta >= mask.GetNbinsX():
        return False
    if iphi < 0 or iphi >= mask.GetNbinsY():
        return False
    return mask.GetBinContent(ieta + 1, iphi + 1) > 0


@dataclass
class Agg:
    n: int = 0
    rj_common: int = 0
    ppg12_common: int = 0
    rj_tight: int = 0
    rj_nontight: int = 0
    ppg12_tight: int = 0
    ppg12_nontight: int = 0
    model_base_e: int = 0
    model_base_v3e: int = 0
    score_diff_sum: float = 0.0
    score_abs_sum: float = 0.0
    rj_score_sum: float = 0.0
    ppg12_score_sum: float = 0.0
    feature_abs: dict[str, float] = field(default_factory=lambda: defaultdict(float))
    feature_signed: dict[str, float] = field(default_factory=lambda: defaultdict(float))
    feature_n: dict[str, int] = field(default_factory=lambda: defaultdict(int))

    def add_feature(self, name: str, rj: float, ppg12: float) -> None:
        if not (finite(rj) and finite(ppg12)):
            return
        diff = float(rj) - float(ppg12)
        self.feature_abs[name] += abs(diff)
        self.feature_signed[name] += diff
        self.feature_n[name] += 1


def read_rj_rows(path: str, events_per_segment: int) -> tuple[list[dict[str, float]], set[tuple[int, int]]]:
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open RecoilJets ROOT: {path}")
    t = f.Get("AuAuPhotonIDTrainingTree")
    if not t:
        raise RuntimeError(f"Missing AuAuPhotonIDTrainingTree in {path}")
    branches = {branch.GetName() for branch in t.GetListOfBranches()}
    has_score_input_et = "cluster_Et_score_input" in branches
    has_cluster_prob = "cluster_prob" in branches

    rows: list[dict[str, float]] = []
    events: set[tuple[int, int]] = set()
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        if int(getattr(t, "is_signal")) != 1:
            continue
        evt = int(getattr(t, "evt"))
        global_eventnumber = int(getattr(t, "eventnumber"))
        if events_per_segment > 0 and global_eventnumber > 0:
            segment = int((global_eventnumber - 1) // events_per_segment)
            eventnumber = int((global_eventnumber - 1) % events_per_segment) + 1
        else:
            segment = int((evt - 1) // events_per_segment) if events_per_segment > 0 and evt > 0 else 0
            eventnumber = global_eventnumber
        row = {
            "segment": segment,
            "evt": evt,
            "eventnumber": eventnumber,
            "global_eventnumber": global_eventnumber,
            "truth_track_id": int(getattr(t, "truth_track_id")),
            "cluster_Et": float(getattr(t, "cluster_Et")),
            "cluster_Et_score_input": (
                float(getattr(t, "cluster_Et_score_input")) if has_score_input_et else float("nan")
            ),
            "cluster_Eta": float(getattr(t, "cluster_Eta")),
            "cluster_Phi": float(getattr(t, "cluster_Phi")),
            "vertexz": float(getattr(t, "ppg12_kin_vertexz")),
            "cluster_weta_cogx": float(getattr(t, "cluster_weta_cogx")),
            "cluster_wphi_cogx": float(getattr(t, "cluster_wphi_cogx")),
            "e11_over_e33": float(getattr(t, "e11_over_e33")),
            "e32_over_e35": float(getattr(t, "e32_over_e35")),
            "cluster_et1": float(getattr(t, "cluster_et1")),
            "cluster_et2": float(getattr(t, "cluster_et2")),
            "cluster_et3": float(getattr(t, "cluster_et3")),
            "cluster_et4": float(getattr(t, "cluster_et4")),
            "cluster_prob": float(getattr(t, "cluster_prob")) if has_cluster_prob else float("nan"),
            "npb_score": float(getattr(t, "npb_score")),
            "tight_bdt_score": float(getattr(t, "tight_bdt_score")),
            "ppg12_common_pass": int(getattr(t, "ppg12_common_pass")),
            "ppg12_tight_tag": int(getattr(t, "ppg12_tight_tag")),
        }
        rows.append(row)
        events.add((segment, eventnumber))
    return rows, events


def read_ppg12_rows(
    path: str,
    wanted_events: set[tuple[int, int]],
    max_entries: int,
    mask_path: str,
    events_per_segment: int,
) -> dict[tuple[int, int], list[dict[str, float]]]:
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open PPG12 ROOT: {path}")
    t = f.Get("slimtree")
    if not t:
        raise RuntimeError(f"Missing slimtree in {path}")

    mf = ROOT.TFile.Open(mask_path) if mask_path else None
    mask = mf.Get("mask_phisymm_tight") if mf else None
    rand = ROOT.TRandom3(0)
    out: dict[tuple[int, int], list[dict[str, float]]] = defaultdict(list)
    stop = t.GetEntries()
    if max_entries > 0:
        stop = min(stop, max_entries)

    for ie in range(stop):
        t.GetEntry(ie)
        eventnumber = int(getattr(t, "eventnumber"))
        segment = int(ie // events_per_segment) if events_per_segment > 0 else 0
        event_key = (segment, eventnumber)
        if event_key not in wanted_events:
            continue
        vertexz = float(getattr(t, "vertexz"))
        if abs(vertexz) > 60.0:
            continue

        ncluster = min(int(getattr(t, "ncluster_CLUSTERINFO_CEMC")), 20000)
        cluster_Et = getattr(t, "cluster_Et_CLUSTERINFO_CEMC")
        cluster_Eta = getattr(t, "cluster_Eta_CLUSTERINFO_CEMC")
        cluster_Phi = getattr(t, "cluster_Phi_CLUSTERINFO_CEMC")
        cluster_prob = getattr(t, "cluster_prob_CLUSTERINFO_CEMC")
        cluster_truth = getattr(t, "cluster_truthtrkID_CLUSTERINFO_CEMC")
        cluster_ietacent = getattr(t, "cluster_ietacent_CLUSTERINFO_CEMC")
        cluster_iphicent = getattr(t, "cluster_iphicent_CLUSTERINFO_CEMC")
        cluster_weta = getattr(t, "cluster_weta_cogx_CLUSTERINFO_CEMC")
        cluster_wphi = getattr(t, "cluster_wphi_cogx_CLUSTERINFO_CEMC")
        cluster_et1 = getattr(t, "cluster_et1_CLUSTERINFO_CEMC")
        cluster_et2 = getattr(t, "cluster_et2_CLUSTERINFO_CEMC")
        cluster_et3 = getattr(t, "cluster_et3_CLUSTERINFO_CEMC")
        cluster_et4 = getattr(t, "cluster_et4_CLUSTERINFO_CEMC")
        cluster_e11 = getattr(t, "cluster_e11_CLUSTERINFO_CEMC")
        cluster_e33 = getattr(t, "cluster_e33_CLUSTERINFO_CEMC")
        cluster_e32 = getattr(t, "cluster_e32_CLUSTERINFO_CEMC")
        cluster_e35 = getattr(t, "cluster_e35_CLUSTERINFO_CEMC")
        cluster_npb = getattr(t, "cluster_npb_score_CLUSTERINFO_CEMC")
        bdt_base_e = getattr(t, "cluster_bdt_CLUSTERINFO_CEMC_base_E")
        bdt_base_v3e = getattr(t, "cluster_bdt_CLUSTERINFO_CEMC_base_v3E")

        # Build signal truth-track set for this event.
        nparticles = min(int(getattr(t, "nparticles")), 20000)
        particle_pid = getattr(t, "particle_pid")
        particle_trkid = getattr(t, "particle_trkid")
        particle_class = getattr(t, "particle_photonclass")
        particle_iso03 = getattr(t, "particle_truth_iso_03")
        signal_tracks = set()
        for ip in range(nparticles):
            if int(particle_pid[ip]) == 22 and int(particle_class[ip]) < 3 and float(particle_iso03[ip]) < 4.0:
                signal_tracks.add(int(particle_trkid[ip]))

        for ic in range(ncluster):
            if tower_masked(mask, float(cluster_ietacent[ic]), float(cluster_iphicent[ic])):
                continue
            raw_et = float(cluster_Et[ic])
            smeared_et = raw_et * rand.Gaus(1.0, 0.04)
            truth_track = int(cluster_truth[ic])
            if truth_track not in signal_tracks:
                continue
            if smeared_et < 5.0:
                continue
            e11 = float(cluster_e11[ic])
            e33 = float(cluster_e33[ic])
            e32 = float(cluster_e32[ic])
            e35 = float(cluster_e35[ic])
            model = selected_ppg12_model(smeared_et)
            score = float(bdt_base_v3e[ic]) if model == "base_v3E" else float(bdt_base_e[ic])
            row = {
                "segment": segment,
                "eventnumber": eventnumber,
                "truth_track_id": truth_track,
                "cluster_Et": raw_et,
                "smeared_cluster_Et": smeared_et,
                "cluster_Eta": float(cluster_Eta[ic]),
                "cluster_Phi": float(cluster_Phi[ic]),
                "vertexz": vertexz,
                "cluster_weta_cogx": float(cluster_weta[ic]),
                "cluster_wphi_cogx": float(cluster_wphi[ic]),
                "e11_over_e33": e11 / e33 if e33 > 0 else 0.0,
                "e32_over_e35": e32 / e35 if e35 > 0 else 0.0,
                "cluster_et1": float(cluster_et1[ic]),
                "cluster_et2": float(cluster_et2[ic]),
                "cluster_et3": float(cluster_et3[ic]),
                "cluster_et4": float(cluster_et4[ic]),
                "cluster_prob": float(cluster_prob[ic]),
                "npb_score": float(cluster_npb[ic]),
                "bdt_base_E": float(bdt_base_e[ic]),
                "bdt_base_v3E": float(bdt_base_v3e[ic]),
                "selected_model": model,
                "selected_bdt_score": score,
            }
            tag, flags = classify(row, score, smeared_et)
            row["ppg12_recomputed_tag"] = tag
            row["ppg12_common_pass"] = 1 if flags.get("common") else 0
            out[event_key].append(row)
    return out


def match_rows(
    rj_rows: list[dict[str, float]],
    ppg12_by_event: dict[tuple[int, int], list[dict[str, float]]],
) -> tuple[list[tuple[dict[str, float], dict[str, float], float]], int, list[dict[str, float]]]:
    matches = []
    misses = 0
    used: set[tuple[int, int, int]] = set()
    for rj in rj_rows:
        event_key = (int(rj["segment"]), int(rj["eventnumber"]))
        candidates = [
            (idx, p)
            for idx, p in enumerate(ppg12_by_event.get(event_key, []))
            if int(p["truth_track_id"]) == int(rj["truth_track_id"])
        ]
        best = None
        best_dist = float("inf")
        for idx, p in candidates:
            key = (event_key[0], event_key[1], idx)
            if key in used:
                continue
            deta = float(rj["cluster_Eta"]) - float(p["cluster_Eta"])
            dph = dphi(float(rj["cluster_Phi"]), float(p["cluster_Phi"]))
            dist = math.hypot(deta, dph)
            if dist < best_dist:
                best_dist = dist
                best = (idx, p)
        if best is None or best_dist > 0.02:
            misses += 1
            continue
        used.add((event_key[0], event_key[1], best[0]))
        matches.append((rj, best[1], best_dist))
    ppg12_only = []
    for (segment, eventnumber), rows in ppg12_by_event.items():
        for idx, row in enumerate(rows):
            if (int(segment), int(eventnumber), idx) not in used:
                ppg12_only.append(row)
    return matches, misses, ppg12_only


def write_reports(
    matches: list[tuple[dict[str, float], dict[str, float], float]],
    misses: int,
    ppg12_only: list[dict[str, float]],
    out_md: Path,
    out_csv: Path,
) -> None:
    aggs = [Agg() for _ in range(len(RECO_BINS) - 1)]
    ppg12_only_by_bin = [
        {"total": 0, "common": 0, "tight": 0, "nontight": 0, "neither": 0, "preselection_fail": 0}
        for _ in range(len(RECO_BINS) - 1)
    ]
    examples = []
    for rj, ppg, dist in matches:
        ib = find_bin(float(ppg["smeared_cluster_Et"]))
        if ib < 0:
            continue
        agg = aggs[ib]
        agg.n += 1
        agg.rj_common += int(rj["ppg12_common_pass"] == 1)
        agg.ppg12_common += int(ppg["ppg12_common_pass"] == 1)
        agg.rj_tight += int(rj["ppg12_tight_tag"] == 1)
        agg.rj_nontight += int(rj["ppg12_tight_tag"] == 2)
        agg.ppg12_tight += int(ppg["ppg12_recomputed_tag"] == 1)
        agg.ppg12_nontight += int(ppg["ppg12_recomputed_tag"] == 2)
        agg.model_base_e += int(ppg["selected_model"] == "base_E")
        agg.model_base_v3e += int(ppg["selected_model"] == "base_v3E")
        score_diff = float(rj["tight_bdt_score"]) - float(ppg["selected_bdt_score"])
        agg.score_diff_sum += score_diff
        agg.score_abs_sum += abs(score_diff)
        agg.rj_score_sum += float(rj["tight_bdt_score"])
        agg.ppg12_score_sum += float(ppg["selected_bdt_score"])
        for rj_name, ppg_name in FEATURES:
            rj_val = float(rj.get(rj_name, float("nan")))
            if rj_name == "cluster_Et_score_input" and not finite(rj_val):
                # Legacy foreground trees did not store this branch; keep the
                # old same-cluster proxy so historical reports remain readable.
                rj_val = float(ppg["cluster_Et"])
            agg.add_feature(rj_name, rj_val, float(ppg[ppg_name]))
        if len(examples) < 30 or abs(score_diff) > min(abs(e[0]) for e in examples):
            examples.append((score_diff, dist, rj, ppg))
            examples = sorted(examples, key=lambda x: abs(x[0]), reverse=True)[:30]

    for ppg in ppg12_only:
        ib = find_bin(float(ppg["smeared_cluster_Et"]))
        if ib < 0:
            continue
        bucket = ppg12_only_by_bin[ib]
        bucket["total"] += 1
        bucket["common"] += int(ppg["ppg12_common_pass"] == 1)
        tag = int(ppg["ppg12_recomputed_tag"])
        if tag == 0:
            bucket["preselection_fail"] += 1
        elif tag == 1:
            bucket["tight"] += 1
        elif tag == 2:
            bucket["nontight"] += 1
        else:
            bucket["neither"] += 1

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "bin", "matches", "rj_common_frac", "ppg12_common_frac",
            "rj_nontight_common", "ppg12_nontight_common",
            "rj_tight_common", "ppg12_tight_common",
            "mean_rj_score", "mean_ppg12_score", "mean_score_diff", "mean_abs_score_diff",
            "base_v3E_rows", "base_E_rows",
        ])
        for i, agg in enumerate(aggs):
            writer.writerow([
                bin_label(i), agg.n,
                safe_div(agg.rj_common, agg.n), safe_div(agg.ppg12_common, agg.n),
                safe_div(agg.rj_nontight, agg.rj_common), safe_div(agg.ppg12_nontight, agg.ppg12_common),
                safe_div(agg.rj_tight, agg.rj_common), safe_div(agg.ppg12_tight, agg.ppg12_common),
                safe_div(agg.rj_score_sum, agg.n), safe_div(agg.ppg12_score_sum, agg.n),
                safe_div(agg.score_diff_sum, agg.n), safe_div(agg.score_abs_sum, agg.n),
                agg.model_base_v3e, agg.model_base_e,
            ])

    lines = [
        "# THE-76 same-cluster PPG12/RecoilJets feature parity",
        "",
        f"- Matched rows: {len(matches)}",
        f"- Unmatched RecoilJets signal rows: {misses}",
        f"- PPG12-only signal rows in same events: {len(ppg12_only)}",
        f"- CSV: `{out_csv}`",
        "",
        "## Stage Summary",
        "",
        "| reco ET bin | matches | RJ common | PPG12 common | RJ NT/common | PPG12 NT/common | RJ tight/common | PPG12 tight/common | mean RJ BDT | mean PPG12 BDT | mean score diff | mean abs score diff | base_v3E | base_E |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for i, agg in enumerate(aggs):
        lines.append(
            f"| {bin_label(i)} | {agg.n} | {safe_div(agg.rj_common, agg.n):.4f} | "
            f"{safe_div(agg.ppg12_common, agg.n):.4f} | "
            f"{safe_div(agg.rj_nontight, agg.rj_common):.4f} | "
            f"{safe_div(agg.ppg12_nontight, agg.ppg12_common):.4f} | "
            f"{safe_div(agg.rj_tight, agg.rj_common):.4f} | "
            f"{safe_div(agg.ppg12_tight, agg.ppg12_common):.4f} | "
            f"{safe_div(agg.rj_score_sum, agg.n):.4f} | "
            f"{safe_div(agg.ppg12_score_sum, agg.n):.4f} | "
            f"{safe_div(agg.score_diff_sum, agg.n):+.4f} | "
            f"{safe_div(agg.score_abs_sum, agg.n):.4f} | "
            f"{agg.model_base_v3e} | {agg.model_base_e} |"
        )

    lines += [
        "",
        "## PPG12-Only Signal Candidates In Same Events",
        "",
        "| reco ET bin | total | common | tight | nontight | neither | preselection fail |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for i, bucket in enumerate(ppg12_only_by_bin):
        lines.append(
            f"| {bin_label(i)} | {bucket['total']} | {bucket['common']} | "
            f"{bucket['tight']} | {bucket['nontight']} | {bucket['neither']} | "
            f"{bucket['preselection_fail']} |"
        )
    lines += ["", "## Mean Absolute Feature Differences", ""]
    for i, agg in enumerate(aggs):
        if agg.n == 0:
            continue
        lines.append(f"### {bin_label(i)}")
        lines.append("")
        lines.append("| feature | mean signed RJ-PPG12 | mean abs |")
        lines.append("| --- | ---: | ---: |")
        for name, _ in FEATURES:
            n = agg.feature_n.get(name, 0)
            lines.append(
                f"| {name} | {safe_div(agg.feature_signed.get(name, 0.0), n):+.6g} | "
                f"{safe_div(agg.feature_abs.get(name, 0.0), n):.6g} |"
            )
        lines.append("")

    lines += ["## Largest Score-Shift Examples", "", "| score diff | dR | segment | event | truth_track | RJ tag | PPG12 tag | RJ score | PPG12 score | raw ET | RJ smeared ET | PPG12 smeared ET | weta RJ/PPG12 | e11e33 RJ/PPG12 | et2 RJ/PPG12 | et3 RJ/PPG12 | et4 RJ/PPG12 | model |", "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- | --- | --- | --- |"]
    for score_diff, dist, rj, ppg in examples:
        lines.append(
            f"| {score_diff:+.5f} | {dist:.5g} | {int(rj['segment'])} | {int(rj['eventnumber'])} | {int(rj['truth_track_id'])} | "
            f"{int(rj['ppg12_tight_tag'])} | {int(ppg['ppg12_recomputed_tag'])} | "
            f"{float(rj['tight_bdt_score']):.5f} | {float(ppg['selected_bdt_score']):.5f} | "
            f"{float(ppg['cluster_Et']):.4f} | {float(rj['cluster_Et']):.4f} | {float(ppg['smeared_cluster_Et']):.4f} | "
            f"{float(rj['cluster_weta_cogx']):.4f}/{float(ppg['cluster_weta_cogx']):.4f} | "
            f"{float(rj['e11_over_e33']):.4f}/{float(ppg['e11_over_e33']):.4f} | "
            f"{float(rj['cluster_et2']):.4f}/{float(ppg['cluster_et2']):.4f} | "
            f"{float(rj['cluster_et3']):.4f}/{float(ppg['cluster_et3']):.4f} | "
            f"{float(rj['cluster_et4']):.4f}/{float(ppg['cluster_et4']):.4f} | {ppg['selected_model']} |"
        )
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rj-root", required=True)
    parser.add_argument("--ppg12-root", default="/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon20/bdt_split.root")
    parser.add_argument("--mask-root", default="/sphenix/user/shuhangli/ppg12/efficiencytool/tower_masks_bdt_nom.root")
    parser.add_argument("--ppg12-max-events", type=int, default=5000)
    parser.add_argument("--events-per-segment", type=int, default=1000)
    parser.add_argument("--out-md", required=True)
    parser.add_argument("--out-csv", required=True)
    args = parser.parse_args()

    rj_rows, wanted_events = read_rj_rows(args.rj_root, args.events_per_segment)
    ppg12_by_event = read_ppg12_rows(
        args.ppg12_root,
        wanted_events,
        args.ppg12_max_events,
        args.mask_root,
        args.events_per_segment,
    )
    matches, misses, ppg12_only = match_rows(rj_rows, ppg12_by_event)
    write_reports(matches, misses, ppg12_only, Path(args.out_md), Path(args.out_csv))
    print(
        f"matched={len(matches)} unmatched={misses} ppg12_only={len(ppg12_only)} "
        f"out_md={args.out_md} out_csv={args.out_csv}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
