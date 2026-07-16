#!/usr/bin/env python3
"""Decompose RecoilJets signal candidates missing from the PPG12 slimtree view.

This local/read-only THE-76 helper compares a small RecoilJets foreground ROOT
against Shuhang's PPG12 bdt_split.root.  It is intentionally stage-oriented:
match by corrected source-file event key, truth track id, and nearest reco
eta/phi, then classify unmatched RecoilJets rows by the first PPG12 gate that
cannot reproduce them.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import ROOT  # type: ignore


ROOT.gROOT.SetBatch(True)

RECO_BINS = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36]
FEATURE_FIELDS = [
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "cluster_et1",
    "cluster_et2",
    "cluster_et3",
    "cluster_et4",
    "e11_over_e33",
    "e32_over_e35",
    "cluster_prob",
    "npb_score",
]


def finite(x: Any) -> bool:
    try:
        return math.isfinite(float(x))
    except Exception:
        return False


def dphi(a: float, b: float) -> float:
    x = a - b
    while x > math.pi:
        x -= 2.0 * math.pi
    while x <= -math.pi:
        x += 2.0 * math.pi
    return x


def find_bin(x: float) -> int:
    for i in range(len(RECO_BINS) - 1):
        if RECO_BINS[i] < x < RECO_BINS[i + 1]:
            return i
    return -1


def bin_label(i: int) -> str:
    return f"{RECO_BINS[i]}-{RECO_BINS[i + 1]}"


def safe_div(num: float, den: float) -> float:
    return num / den if den else float("nan")


def open_interval(x: float, lo: float, hi: float) -> bool:
    return finite(x) and lo < x < hi


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


def selected_ppg12_model(et: float) -> str:
    return "base_v3E" if 8.0 <= et < 35.0 else "base_E"


def ppg12_thresholds(et: float) -> tuple[float, float, float]:
    tight_min = 0.815625 - 0.0015625 * et
    nt_min = 0.7333333333333333 - 0.01333333333333333 * et
    nt_max = 0.684375 + 0.0015625 * et
    return tight_min, nt_min, nt_max


def classify_ppg12(row: dict[str, float], score: float, et_for_cuts: float) -> tuple[str, dict[str, bool]]:
    e11e33 = row["e11_over_e33"]
    e32e35 = row["e32_over_e35"]
    weta = row["cluster_weta_cogx"]
    wphi = row["cluster_wphi_cogx"]
    wr = wphi / weta if weta != 0.0 else float("nan")
    common = (
        open_interval(row["cluster_prob"], 0.0, 1.0)
        and open_interval(e11e33, 0.0, 0.98)
        and finite(wr)
        and wr > 0.0
        and finite(weta)
        and weta < 2.0
        and row["npb_score"] > 0.5
    )
    flags = {"common": common}
    if not common:
        return "preselection_fail", flags

    tight_min, nt_min, nt_max = ppg12_thresholds(et_for_cuts)
    tight_weta = open_interval(weta, 0.0, 1.0)
    tight_wphi = open_interval(wphi, 0.0, 1.0)
    tight_et1 = open_interval(row["cluster_et1"], 0.5, 1.0)
    tight_et2 = open_interval(row["cluster_et2"], 0.0, 1.0)
    tight_et3 = open_interval(row["cluster_et3"], 0.0, 1.0)
    tight_et4 = open_interval(row["cluster_et4"], 0.0, 1.0)
    tight_prob = open_interval(row["cluster_prob"], 0.0, 1.0)
    tight_e11e33 = open_interval(e11e33, 0.0, 1.0)
    tight_e32e35 = open_interval(e32e35, 0.8, 1.0)
    tight_bdt = finite(score) and score > tight_min and score < 1.0
    tight = (
        tight_weta
        and tight_wphi
        and tight_et1
        and tight_et2
        and tight_et3
        and tight_et4
        and tight_prob
        and tight_e11e33
        and tight_e32e35
        and tight_bdt
    )
    flags.update(
        tight=tight,
        tight_weta=tight_weta,
        tight_prob=tight_prob,
        tight_bdt=tight_bdt,
        nt_bdt=finite(score) and score > nt_min and score < nt_max,
    )
    if tight:
        return "tight", flags

    nt_shape = (
        open_interval(weta, 0.0, 1.0)
        and open_interval(wphi, 0.0, 1.0)
        and open_interval(row["cluster_prob"], 0.0, 1.0)
        and open_interval(e11e33, 0.0, 1.0)
        and open_interval(e32e35, 0.8, 1.0)
        and open_interval(row["cluster_et1"], 0.6, 1.0)
        and open_interval(row["cluster_et4"], 0.0, 1.0)
    )
    nfail = 0
    if not tight_weta:
        nfail += 1
    if not tight_prob:
        nfail += 1
    if not tight_bdt:
        nfail += 1
    flags.update(nt_shape=nt_shape, nfail_gt0=nfail > 0)
    if nt_shape and flags["nt_bdt"] and nfail > 0:
        return "nontight", flags
    return "neither", flags


@dataclass
class EventInfo:
    vertexz: float
    pass_vertex: bool
    ppg12_entry: int


def read_rj_rows(path: str, events_per_segment: int) -> tuple[list[dict[str, float]], set[tuple[int, int]]]:
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open RecoilJets ROOT: {path}")
    t = f.Get("AuAuPhotonIDTrainingTree")
    if not t:
        raise RuntimeError(f"Missing AuAuPhotonIDTrainingTree in {path}")
    branches = {b.GetName() for b in t.GetListOfBranches()}
    rows: list[dict[str, float]] = []
    events: set[tuple[int, int]] = set()
    for ir in range(t.GetEntries()):
        t.GetEntry(ir)
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
            "row": ir,
            "evt": evt,
            "global_eventnumber": global_eventnumber,
            "segment": segment,
            "eventnumber": eventnumber,
            "vertexz": float(getattr(t, "vertexz")) if "vertexz" in branches else float("nan"),
            "truth_track_id": int(getattr(t, "truth_track_id")),
            "cluster_Et": float(getattr(t, "cluster_Et")),
            "cluster_Eta": float(getattr(t, "cluster_Eta")),
            "cluster_Phi": float(getattr(t, "cluster_Phi")),
            "ppg12_kin_vertexz": float(getattr(t, "ppg12_kin_vertexz")) if "ppg12_kin_vertexz" in branches else float("nan"),
            "common": int(getattr(t, "ppg12_common_pass")) if "ppg12_common_pass" in branches else -1,
            "tag": int(getattr(t, "ppg12_tight_tag")) if "ppg12_tight_tag" in branches else -1,
            "iso": int(getattr(t, "ppg12_is_iso")) if "ppg12_is_iso" in branches else -1,
            "noniso": int(getattr(t, "ppg12_is_noniso")) if "ppg12_is_noniso" in branches else -1,
            "truth_window": int(getattr(t, "ppg12_truth_window_pass_r04")) if "ppg12_truth_window_pass_r04" in branches else -1,
            "bdt_score": float(getattr(t, "tight_bdt_score")) if "tight_bdt_score" in branches else float("nan"),
            "score_input_et": float(getattr(t, "score_input_et")) if "score_input_et" in branches else float("nan"),
        }
        for feature in FEATURE_FIELDS:
            row[feature] = float(getattr(t, feature)) if feature in branches else float("nan")
        rows.append(row)
        events.add((segment, eventnumber))
    return rows, events


def read_ppg12_rows(
    path: str,
    mask_path: str,
    entries: int,
    events_per_segment: int,
) -> tuple[dict[tuple[int, int], EventInfo], dict[tuple[int, int], list[dict[str, float]]], dict[tuple[int, int], list[dict[str, float]]]]:
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open PPG12 ROOT: {path}")
    t = f.Get("slimtree")
    if not t:
        raise RuntimeError(f"Missing slimtree in {path}")
    mf = ROOT.TFile.Open(mask_path)
    if not mf or mf.IsZombie():
        raise RuntimeError(f"Could not open mask ROOT: {mask_path}")
    mask = mf.Get("mask_phisymm_tight")
    if not mask:
        raise RuntimeError("Missing mask_phisymm_tight")

    event_info: dict[tuple[int, int], EventInfo] = {}
    pre_by_event: dict[tuple[int, int], list[dict[str, float]]] = defaultdict(list)
    pass_by_event: dict[tuple[int, int], list[dict[str, float]]] = defaultdict(list)
    rand = ROOT.TRandom3(0)
    stop = t.GetEntries() if entries <= 0 else min(entries, t.GetEntries())
    for ie in range(stop):
        t.GetEntry(ie)
        segment = int(ie // events_per_segment) if events_per_segment > 0 else 0
        eventnumber = int(getattr(t, "eventnumber"))
        event_key = (segment, eventnumber)
        vertexz = float(getattr(t, "vertexz"))
        pass_vertex = abs(vertexz) <= 60.0
        event_info[event_key] = EventInfo(vertexz=vertexz, pass_vertex=pass_vertex, ppg12_entry=ie)

        nparticles = min(int(getattr(t, "nparticles")), 20000)
        particle_pid = getattr(t, "particle_pid")
        particle_trkid = getattr(t, "particle_trkid")
        particle_class = getattr(t, "particle_photonclass")
        particle_iso03 = getattr(t, "particle_truth_iso_03")
        signal_tracks = set()
        signal_truth_pt: dict[int, float] = {}
        for ip in range(nparticles):
            if int(particle_pid[ip]) == 22 and int(particle_class[ip]) < 3 and float(particle_iso03[ip]) < 4.0:
                track_id = int(particle_trkid[ip])
                signal_tracks.add(track_id)
                eta = float(getattr(t, "particle_Eta")[ip])
                energy = float(getattr(t, "particle_E")[ip])
                signal_truth_pt[track_id] = energy / math.cosh(eta) if finite(eta) else float("nan")

        ncluster = min(int(getattr(t, "ncluster_CLUSTERINFO_CEMC")), 20000)
        cluster_Et = getattr(t, "cluster_Et_CLUSTERINFO_CEMC")
        cluster_Eta = getattr(t, "cluster_Eta_CLUSTERINFO_CEMC")
        cluster_Phi = getattr(t, "cluster_Phi_CLUSTERINFO_CEMC")
        cluster_prob = getattr(t, "cluster_prob_CLUSTERINFO_CEMC")
        cluster_truth = getattr(t, "cluster_truthtrkID_CLUSTERINFO_CEMC")
        cluster_ietacent = getattr(t, "cluster_ietacent_CLUSTERINFO_CEMC")
        cluster_iphicent = getattr(t, "cluster_iphicent_CLUSTERINFO_CEMC")
        cluster_iso_topo_04 = getattr(t, "cluster_iso_topo_04_CLUSTERINFO_CEMC")
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

        smeared_ets: list[float | None] = [None] * ncluster
        if pass_vertex:
            for ic in range(ncluster):
                if tower_masked(mask, float(cluster_ietacent[ic]), float(cluster_iphicent[ic])):
                    continue
                smeared_ets[ic] = float(cluster_Et[ic]) * rand.Gaus(1.0, 0.04)

        for ic in range(ncluster):
            if tower_masked(mask, float(cluster_ietacent[ic]), float(cluster_iphicent[ic])):
                continue
            truth_track = int(cluster_truth[ic])
            if truth_track not in signal_tracks:
                continue
            e11 = float(cluster_e11[ic])
            e33 = float(cluster_e33[ic])
            e32 = float(cluster_e32[ic])
            e35 = float(cluster_e35[ic])
            raw_et = float(cluster_Et[ic])
            base_row = {
                "segment": segment,
                "eventnumber": eventnumber,
                "ppg12_entry": ie,
                "icluster": ic,
                "truth_track_id": truth_track,
                "raw_et": raw_et,
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
                "eiso": 1.2 * float(cluster_iso_topo_04[ic]) + 0.1,
                "truth_pt": signal_truth_pt.get(truth_track, float("nan")),
            }
            pre_by_event[event_key].append(dict(base_row))
            smeared_et = smeared_ets[ic]
            if not pass_vertex or smeared_et is None or smeared_et < 5.0:
                continue
            model = selected_ppg12_model(smeared_et)
            score = float(bdt_base_v3e[ic]) if model == "base_v3E" else float(bdt_base_e[ic])
            row = dict(base_row)
            row.update(
                smeared_et=smeared_et,
                selected_model=model,
                selected_bdt_score=score,
            )
            tag, flags = classify_ppg12(row, score, smeared_et)
            iso_max = 0.490 + 0.037 * smeared_et
            noniso_min = iso_max + 0.8
            is_iso = row["eiso"] > -20.0 and row["eiso"] < iso_max
            is_noniso = row["eiso"] > noniso_min and row["eiso"] < 20.0
            in_response_window = (
                finite(row["truth_pt"])
                and 8.0 < row["truth_pt"] < 45.0
                and 10.0 < smeared_et < 36.0
            )
            row["tag_name"] = tag
            row["common"] = 1 if flags.get("common") else 0
            row["tag"] = {"preselection_fail": 0, "tight": 1, "nontight": 2, "neither": 3}[tag]
            row["is_iso"] = int(is_iso)
            row["is_noniso"] = int(is_noniso)
            row["in_response_window"] = int(in_response_window)
            pass_by_event[event_key].append(row)

    return event_info, pre_by_event, pass_by_event


def nearest_same_track(rj: dict[str, float], rows: list[dict[str, float]]) -> tuple[dict[str, float] | None, float]:
    best = None
    best_dr = float("inf")
    for p in rows:
        if int(p["truth_track_id"]) != int(rj["truth_track_id"]):
            continue
        dr = math.hypot(
            float(rj["cluster_Eta"]) - float(p["cluster_Eta"]),
            dphi(float(rj["cluster_Phi"]), float(p["cluster_Phi"])),
        )
        if dr < best_dr:
            best_dr = dr
            best = p
    return best, best_dr


def row_region(row: dict[str, float]) -> str:
    tag = int(row.get("tag", -1))
    iso = int(row.get("iso", -1))
    noniso = int(row.get("noniso", -1))
    if tag == 1 and iso == 1:
        return "A_like_tight_iso"
    if tag == 1 and noniso == 1:
        return "B_like_tight_noniso"
    if tag == 2 and iso == 1:
        return "C_like_nontight_iso"
    if tag == 2 and noniso == 1:
        return "D_like_nontight_noniso"
    if tag == 1:
        return "tight_gap"
    if tag == 2:
        return "nontight_gap"
    if tag == 3:
        return "neither"
    if tag == 0:
        return "preselection_fail"
    return "unknown"


def decompose(
    rj_rows: list[dict[str, float]],
    event_info: dict[tuple[int, int], EventInfo],
    pre_by_event: dict[tuple[int, int], list[dict[str, float]]],
    pass_by_event: dict[tuple[int, int], list[dict[str, float]]],
) -> tuple[list[dict[str, Any]], Counter, dict[int, Counter], set[tuple[int, int, int]]]:
    used_pass: set[tuple[int, int, int]] = set()
    out: list[dict[str, Any]] = []
    category_counter: Counter = Counter()
    by_bin: dict[int, Counter] = defaultdict(Counter)

    for rj in rj_rows:
        event_key = (int(rj["segment"]), int(rj["eventnumber"]))
        info = event_info.get(event_key)
        pass_rows = pass_by_event.get(event_key, [])
        pre_rows = pre_by_event.get(event_key, [])

        best_pass, dr_pass = nearest_same_track(rj, pass_rows)
        matched = False
        if best_pass is not None and dr_pass <= 0.02:
            key = (event_key[0], event_key[1], int(best_pass["icluster"]))
            if key not in used_pass:
                used_pass.add(key)
                matched = True

        best_pre, dr_pre = nearest_same_track(rj, pre_rows)
        if matched:
            category = "matched_ppg12_pass_dr002"
        elif info is None:
            category = "ppg12_event_not_in_span"
        elif not info.pass_vertex:
            category = "ppg12_reco_vertex_fail"
        elif best_pre is None:
            category = "no_ppg12_signal_track_in_event_pre_vtx"
        elif best_pass is None:
            category = "same_track_pre_vtx_but_removed_before_pass_rows"
        elif dr_pass <= 0.10:
            category = "same_track_pass_cluster_dr002_to_dr010"
        else:
            category = "same_track_pass_cluster_dr_gt010"

        ib = find_bin(float(rj["cluster_Et"]))
        region = row_region(rj)
        category_counter[category] += 1
        by_bin[ib][category] += 1
        by_bin[ib][f"region:{region}"] += 1
        if int(rj.get("common", -1)) == 1:
            by_bin[ib]["common"] += 1
        if int(rj.get("tag", -1)) == 1:
            by_bin[ib]["tight"] += 1
        if int(rj.get("tag", -1)) == 2:
            by_bin[ib]["nontight"] += 1

        out.append(
            {
                **rj,
                "category": category,
                "ppg12_vertexz": info.vertexz if info else float("nan"),
                "ppg12_vertex_pass": int(info.pass_vertex) if info else -1,
                "ppg12_icluster": int(best_pass["icluster"]) if best_pass is not None else -1,
                "ppg12_raw_et": best_pass["raw_et"] if best_pass is not None else float("nan"),
                "ppg12_smeared_et": best_pass["smeared_et"] if best_pass is not None else float("nan"),
                "ppg12_selected_model": best_pass["selected_model"] if best_pass is not None else "",
                "ppg12_selected_bdt_score": best_pass["selected_bdt_score"] if best_pass is not None else float("nan"),
                "ppg12_tag": best_pass["tag"] if best_pass is not None else -1,
                "ppg12_tag_name": best_pass["tag_name"] if best_pass is not None else "",
                "ppg12_common": best_pass["common"] if best_pass is not None else -1,
                "nearest_pre_dr": dr_pre if best_pre is not None else float("nan"),
                "nearest_pass_dr": dr_pass if best_pass is not None else float("nan"),
                "matched_ppg12_pass": int(matched),
                "region": region,
                **{
                    f"ppg12_{feature}": best_pass[feature]
                    if best_pass is not None and feature in best_pass else float("nan")
                    for feature in FEATURE_FIELDS
                },
            }
        )
    return out, category_counter, by_bin, used_pass


def summarize_ppg12_pass_rows(
    pass_by_event: dict[tuple[int, int], list[dict[str, float]]],
    used_pass: set[tuple[int, int, int]],
) -> tuple[Counter, dict[int, Counter]]:
    overall: Counter = Counter()
    by_bin: dict[int, Counter] = defaultdict(Counter)
    for event_key, rows in pass_by_event.items():
        for row in rows:
            key = (event_key[0], event_key[1], int(row["icluster"]))
            tag_name = str(row.get("tag_name", "unknown"))
            ib = find_bin(float(row.get("smeared_et", float("nan"))))
            matched = key in used_pass
            state = "matched" if matched else "unmatched"
            overall["total"] += 1
            overall[f"tag:{tag_name}"] += 1
            overall[state] += 1
            overall[f"{state}:{tag_name}"] += 1
            by_bin[ib]["total"] += 1
            by_bin[ib][f"tag:{tag_name}"] += 1
            by_bin[ib][state] += 1
            by_bin[ib][f"{state}:{tag_name}"] += 1
            if row.get("tag") == 1 and row.get("is_iso") == 1 and row.get("in_response_window") == 1:
                by_bin[ib]["region:A"] += 1
            if row.get("tag") == 1 and row.get("is_noniso") == 1:
                by_bin[ib]["region:B"] += 1
            if row.get("tag") == 2 and row.get("is_iso") == 1:
                by_bin[ib]["region:C"] += 1
            if row.get("tag") == 2 and row.get("is_noniso") == 1:
                by_bin[ib]["region:D"] += 1
    return overall, by_bin


def write_outputs(
    rows: list[dict[str, Any]],
    category_counter: Counter,
    by_bin: dict[int, Counter],
    ppg12_overall: Counter,
    ppg12_by_bin: dict[int, Counter],
    out_md: Path,
    out_csv: Path,
) -> None:
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "row", "segment", "eventnumber", "global_eventnumber", "truth_track_id",
        "cluster_Et", "cluster_Eta", "cluster_Phi", "vertexz", "ppg12_kin_vertexz",
        "ppg12_vertexz", "ppg12_vertex_pass", "common", "tag", "iso", "noniso",
        "truth_window", "bdt_score", "score_input_et", "region", "category",
        "ppg12_icluster", "ppg12_raw_et", "ppg12_smeared_et",
        "ppg12_selected_model", "ppg12_selected_bdt_score", "ppg12_tag",
        "ppg12_tag_name", "ppg12_common", "nearest_pre_dr", "nearest_pass_dr",
        "matched_ppg12_pass",
    ]
    for feature in FEATURE_FIELDS:
        fieldnames.extend([feature, f"ppg12_{feature}"])
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    lines = [
        "# THE-76 unmatched RecoilJets signal-candidate decomposition",
        "",
        f"- RecoilJets signal rows: {len(rows)}",
        f"- CSV: `{out_csv}`",
        "",
        "## Overall Categories",
        "",
        "| category | rows | fraction |",
        "| --- | ---: | ---: |",
    ]
    for cat, n in category_counter.most_common():
        lines.append(f"| {cat} | {n} | {safe_div(n, len(rows)):.4f} |")

    category_names = [cat for cat, _ in category_counter.most_common()]
    lines += [
        "",
        "## By RecoilJets ET Bin",
        "",
        "| bin | total | common | tight | nontight | " + " | ".join(category_names) + " |",
        "| --- | ---: | ---: | ---: | ---: | " + " | ".join(["---:"] * len(category_names)) + " |",
    ]
    for ib in range(len(RECO_BINS) - 1):
        c = by_bin.get(ib, Counter())
        total = sum(c.get(cat, 0) for cat in category_names)
        vals = [str(c.get(cat, 0)) for cat in category_names]
        lines.append(
            f"| {bin_label(ib)} | {total} | {c.get('common', 0)} | {c.get('tight', 0)} | {c.get('nontight', 0)} | "
            + " | ".join(vals)
            + " |"
        )

    region_names = sorted({key.split(":", 1)[1] for c in by_bin.values() for key in c if key.startswith("region:")})
    lines += [
        "",
        "## Region Makeup By RecoilJets ET Bin",
        "",
        "| bin | " + " | ".join(region_names) + " |",
        "| --- | " + " | ".join(["---:"] * len(region_names)) + " |",
    ]
    for ib in range(len(RECO_BINS) - 1):
        c = by_bin.get(ib, Counter())
        vals = [str(c.get(f"region:{name}", 0)) for name in region_names]
        lines.append(f"| {bin_label(ib)} | " + " | ".join(vals) + " |")

    lines += [
        "",
        "## PPG12 Pass Rows Reverse Match",
        "",
        f"- PPG12 pass signal rows in inspected span: {ppg12_overall.get('total', 0)}",
        f"- Exact RecoilJets matches used: {ppg12_overall.get('matched', 0)}",
        f"- PPG12 pass rows unmatched by exact RecoilJets cluster match: {ppg12_overall.get('unmatched', 0)}",
        "",
        "| tag | total | matched | unmatched |",
        "| --- | ---: | ---: | ---: |",
    ]
    for tag_name in ["tight", "nontight", "neither", "preselection_fail"]:
        total = ppg12_overall.get(f"tag:{tag_name}", 0)
        matched = ppg12_overall.get(f"matched:{tag_name}", 0)
        unmatched = ppg12_overall.get(f"unmatched:{tag_name}", 0)
        if total or matched or unmatched:
            lines.append(f"| {tag_name} | {total} | {matched} | {unmatched} |")

    lines += [
        "",
        "## PPG12 Pass Rows Reverse Match By ET Bin",
        "",
        "| bin | total | matched | unmatched | tight | nontight | A | B | C | D | unmatched_tight | unmatched_nontight |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for ib in range(len(RECO_BINS) - 1):
        c = ppg12_by_bin.get(ib, Counter())
        lines.append(
            f"| {bin_label(ib)} | {c.get('total', 0)} | {c.get('matched', 0)} | {c.get('unmatched', 0)} | "
            f"{c.get('tag:tight', 0)} | {c.get('tag:nontight', 0)} | "
            f"{c.get('region:A', 0)} | {c.get('region:B', 0)} | "
            f"{c.get('region:C', 0)} | {c.get('region:D', 0)} | "
            f"{c.get('unmatched:tight', 0)} | {c.get('unmatched:nontight', 0)} |"
        )

    out_md.write_text("\n".join(lines) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rj-root", required=True)
    ap.add_argument("--ppg12-root", default="/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon20/bdt_split.root")
    ap.add_argument("--mask-root", default="/sphenix/user/shuhangli/ppg12/efficiencytool/tower_masks_bdt_nom.root")
    ap.add_argument("--entries", type=int, default=2000)
    ap.add_argument("--events-per-segment", type=int, default=1000)
    ap.add_argument("--out-md", required=True)
    ap.add_argument("--out-csv", required=True)
    args = ap.parse_args()

    rj_rows, _ = read_rj_rows(args.rj_root, args.events_per_segment)
    event_info, pre_by_event, pass_by_event = read_ppg12_rows(
        args.ppg12_root, args.mask_root, args.entries, args.events_per_segment
    )
    rows, category_counter, by_bin, used_pass = decompose(
        rj_rows, event_info, pre_by_event, pass_by_event
    )
    ppg12_overall, ppg12_by_bin = summarize_ppg12_pass_rows(pass_by_event, used_pass)
    write_outputs(
        rows,
        category_counter,
        by_bin,
        ppg12_overall,
        ppg12_by_bin,
        Path(args.out_md),
        Path(args.out_csv),
    )
    print(
        f"rj_signal_rows={len(rows)} categories={dict(category_counter)} "
        f"ppg12_pass_rows={ppg12_overall.get('total', 0)} "
        f"ppg12_unmatched={ppg12_overall.get('unmatched', 0)} "
        f"out_md={args.out_md} out_csv={args.out_csv}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
