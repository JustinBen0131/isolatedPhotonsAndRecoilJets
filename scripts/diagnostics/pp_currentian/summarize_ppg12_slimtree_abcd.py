#!/usr/bin/env python3
"""Recompute PPG12 Fig.29 signal ABCD leakage from a PPG12 slimtree.

This is a read-only THE-76 forensic helper.  It mirrors the relevant
RecoEffCalculator_TTreeReader.C logic closely enough to compare a small
foreground RecoilJets sample against the same PPG12 slimtree entry span.
"""

from __future__ import annotations

import argparse
import math
from collections import defaultdict
from dataclasses import dataclass
from typing import Any

import ROOT  # type: ignore


ROOT.gROOT.SetBatch(True)

RECO_BINS = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36]
TRUTH_BINS = [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36, 45]


def finite(x: Any) -> bool:
    try:
        return math.isfinite(float(x))
    except Exception:
        return False


def open_interval(x: float, lo: float, hi: float) -> bool:
    return finite(x) and lo < x < hi


def find_bin(x: float, bins: list[float] = RECO_BINS) -> int:
    for i in range(len(bins) - 1):
        if bins[i] < x < bins[i + 1]:
            return i
    return -1


def bin_label(i: int) -> str:
    return f"{RECO_BINS[i]}-{RECO_BINS[i + 1]}"


def safe_div(num: float, den: float) -> float:
    return num / den if den else float("nan")


def selected_ppg12_model(et: float) -> str:
    # config_bdt_nom.yaml: fallback base_E, ET-binned base_v3E for [8,15), [15,35).
    return "base_v3E" if 8.0 <= et < 35.0 else "base_E"


def ppg12_thresholds(et: float) -> tuple[float, float, float]:
    tight_min = 0.8333333333333334 - 0.003333333333333336 * et
    nt_min = 0.7333333333333333 - 0.01333333333333333 * et
    nt_max = 0.6666666666666666 + 0.003333333333333336 * et
    return tight_min, nt_min, nt_max


def classify(row: dict[str, float], score: float, et: float) -> tuple[str, dict[str, bool]]:
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

    tight_min, nt_min, nt_max = ppg12_thresholds(et)
    tight_weta = open_interval(weta, 0.0, 1.0)
    tight_wphi = open_interval(wphi, 0.0, 1.0)
    tight_et1 = open_interval(row["cluster_et1"], 0.5, 1.0)
    tight_et2 = open_interval(row["cluster_et2"], 0.0, 1.0)
    tight_et3 = open_interval(row["cluster_et3"], 0.0, 1.0)
    tight_e11e33 = open_interval(e11e33, 0.0, 1.0)
    tight_e32e35 = open_interval(e32e35, 0.8, 1.0)
    tight_et4 = open_interval(row["cluster_et4"], 0.0, 1.0)
    tight_prob = open_interval(row["cluster_prob"], 0.0, 1.0)
    tight_bdt = finite(score) and score > tight_min and score < 1.0
    tight = (
        tight_weta
        and tight_wphi
        and tight_et1
        and tight_et2
        and tight_et3
        and tight_e11e33
        and tight_e32e35
        and tight_et4
        and tight_prob
        and tight_bdt
    )
    flags.update(
        tight=tight,
        tight_weta=tight_weta,
        tight_wphi=tight_wphi,
        tight_et1=tight_et1,
        tight_et2=tight_et2,
        tight_et3=tight_et3,
        tight_e11e33=tight_e11e33,
        tight_e32e35=tight_e32e35,
        tight_et4=tight_et4,
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
    # config_bdt_nom.yaml increments nfail for weta and BDT toggles; prob is
    # hard-coded in PPG12 but cannot fail after the broad prob prefilter.
    nfail = 0
    if not tight_weta:
        nfail += 1
    if not tight_prob:
        nfail += 1
    if not tight_bdt:
        nfail += 1
    nontight = nt_shape and flags["nt_bdt"] and nfail > 0
    flags.update(nt_shape=nt_shape, nontight=nontight, nfail_gt0=nfail > 0)
    if nontight:
        return "nontight", flags
    return "neither", flags


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
class Counts:
    signal: float = 0.0
    common: float = 0.0
    tight: float = 0.0
    nontight: float = 0.0
    neither: float = 0.0
    iso: float = 0.0
    noniso: float = 0.0
    A: float = 0.0
    B: float = 0.0
    C: float = 0.0
    D: float = 0.0
    response_signal: float = 0.0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ppg12-root", required=True)
    ap.add_argument("--mask-root", required=True)
    ap.add_argument("--rj-root", default="", help="optional RecoilJets ROOT used to restrict counted events")
    ap.add_argument("--events-per-segment", type=int, default=1000)
    ap.add_argument("--entries", type=int, default=2000, help="slimtree entries to process; <=0 means all")
    ap.add_argument("--weighted", action="store_true", help="use event_weight branch if present")
    args = ap.parse_args()

    f = ROOT.TFile.Open(args.ppg12_root)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open {args.ppg12_root}")
    t = f.Get("slimtree")
    if not t:
        raise RuntimeError(f"Missing slimtree in {args.ppg12_root}")
    branches = {b.GetName() for b in t.GetListOfBranches()}

    mf = ROOT.TFile.Open(args.mask_root)
    if not mf or mf.IsZombie():
        raise RuntimeError(f"Could not open mask {args.mask_root}")
    mask = mf.Get("mask_phisymm_tight")
    if not mask:
        raise RuntimeError(f"Missing mask_phisymm_tight in {args.mask_root}")

    rand = ROOT.TRandom3(0)
    counts: list[Counts] = [Counts() for _ in range(len(RECO_BINS) - 1)]
    totals = defaultdict(float)
    wanted_events: set[tuple[int, int]] | None = None
    if args.rj_root:
        rf = ROOT.TFile.Open(args.rj_root)
        if not rf or rf.IsZombie():
            raise RuntimeError(f"Could not open RecoilJets event-mask ROOT: {args.rj_root}")
        rt = rf.Get("AuAuPhotonIDTrainingTree")
        if not rt:
            raise RuntimeError(f"Missing AuAuPhotonIDTrainingTree in {args.rj_root}")
        wanted_events = set()
        for ir in range(rt.GetEntries()):
            rt.GetEntry(ir)
            evt = int(getattr(rt, "evt"))
            global_eventnumber = int(getattr(rt, "eventnumber"))
            if args.events_per_segment > 0 and global_eventnumber > 0:
                segment = int((global_eventnumber - 1) // args.events_per_segment)
                eventnumber = int((global_eventnumber - 1) % args.events_per_segment) + 1
            else:
                segment = int((evt - 1) // args.events_per_segment) if args.events_per_segment > 0 and evt > 0 else 0
                eventnumber = global_eventnumber
            wanted_events.add((segment, eventnumber))

    nentries = t.GetEntries()
    stop = nentries if args.entries <= 0 else min(args.entries, nentries)

    for ie in range(stop):
        t.GetEntry(ie)
        eventnumber = int(getattr(t, "eventnumber"))
        segment = int(ie // args.events_per_segment) if args.events_per_segment > 0 else 0
        count_this_event = wanted_events is None or (segment, eventnumber) in wanted_events
        vertexz = float(getattr(t, "vertexz"))
        event_weight = 1.0
        if args.weighted and "weight" in branches:
            event_weight = float(getattr(t, "weight"))
        elif args.weighted and "event_weight" in branches:
            event_weight = float(getattr(t, "event_weight"))
        ncluster = int(getattr(t, "ncluster_CLUSTERINFO_CEMC"))
        cluster_Et = getattr(t, "cluster_Et_CLUSTERINFO_CEMC")
        cluster_Eta = getattr(t, "cluster_Eta_CLUSTERINFO_CEMC")
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

        smeared_et = [float("nan")] * ncluster
        unmasked = [False] * ncluster
        for ic in range(ncluster):
            if tower_masked(mask, float(cluster_ietacent[ic]), float(cluster_iphicent[ic])):
                continue
            unmasked[ic] = True
            smeared_et[ic] = float(cluster_Et[ic]) * rand.Gaus(1.0, 0.04)

        if not count_this_event or abs(vertexz) > 60.0:
            continue

        nparticles = int(getattr(t, "nparticles"))
        particle_pid = getattr(t, "particle_pid")
        particle_trkid = getattr(t, "particle_trkid")
        particle_pt = getattr(t, "particle_Pt")
        particle_class = getattr(t, "particle_photonclass")
        particle_iso03 = getattr(t, "particle_truth_iso_03")
        signal_by_track: dict[int, float] = {}
        for ip in range(nparticles):
            if (
                int(particle_pid[ip]) == 22
                and int(particle_class[ip]) < 3
                and float(particle_iso03[ip]) < 4.0
            ):
                signal_by_track[int(particle_trkid[ip])] = float(particle_pt[ip])

        for ic in range(ncluster):
            if not unmasked[ic]:
                continue
            et = smeared_et[ic]
            if et < 5.0:
                continue
            eta = float(cluster_Eta[ic])
            if not (-0.7 < eta < 0.7):
                continue
            reco_bin = find_bin(et)
            if reco_bin < 0:
                continue

            truth_track = int(cluster_truth[ic])
            if truth_track not in signal_by_track:
                continue
            truth_pt = signal_by_track[truth_track]
            row = {
                "cluster_prob": float(cluster_prob[ic]),
                "cluster_weta_cogx": float(cluster_weta[ic]),
                "cluster_wphi_cogx": float(cluster_wphi[ic]),
                "cluster_et1": float(cluster_et1[ic]),
                "cluster_et2": float(cluster_et2[ic]),
                "cluster_et3": float(cluster_et3[ic]),
                "cluster_et4": float(cluster_et4[ic]),
                "e11_over_e33": float(cluster_e11[ic]) / float(cluster_e33[ic]) if float(cluster_e33[ic]) > 0 else 0.0,
                "e32_over_e35": float(cluster_e32[ic]) / float(cluster_e35[ic]) if float(cluster_e35[ic]) > 0 else 0.0,
                "npb_score": float(cluster_npb[ic]),
            }
            model = selected_ppg12_model(et)
            score = float(bdt_base_v3e[ic]) if model == "base_v3E" else float(bdt_base_e[ic])
            tag, flags = classify(row, score, et)
            raw_eiso = float(cluster_iso_topo_04[ic])
            eiso = 1.2 * raw_eiso + 0.1
            iso_thr = 0.490 + 0.037 * et
            noniso_thr = iso_thr + 0.8
            iso = eiso > -20.0 and eiso < iso_thr
            noniso = eiso > noniso_thr and eiso < 20.0
            in_response = (
                TRUTH_BINS[0] < truth_pt < TRUTH_BINS[-1]
                and RECO_BINS[0] < et < RECO_BINS[-1]
            )

            c = counts[reco_bin]
            c.signal += event_weight
            totals["signal"] += event_weight
            if in_response:
                c.response_signal += event_weight
            if flags.get("common"):
                c.common += event_weight
            if tag == "tight":
                c.tight += event_weight
            elif tag == "nontight":
                c.nontight += event_weight
            elif tag == "neither":
                c.neither += event_weight
            if iso:
                c.iso += event_weight
            if noniso:
                c.noniso += event_weight

            if tag == "tight" and iso and in_response:
                c.A += event_weight
            elif tag == "tight" and noniso:
                c.B += event_weight
            elif tag == "nontight" and iso:
                c.C += event_weight
            elif tag == "nontight" and noniso:
                c.D += event_weight

    print(f"source={args.ppg12_root}")
    print(f"entries_processed={stop}")
    if args.rj_root:
        print(f"event_mask={args.rj_root}")
        print(f"event_mask_size={len(wanted_events or [])}")
    print(f"weighted={args.weighted}")
    print("| bin | signal | response | common | tight | nontight | neither | A | B | C | D | B/A | C/A | D/A | (C+D)/A | NT/common |")
    print("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for i, c in enumerate(counts):
        print(
            f"| {bin_label(i)} | {c.signal:.6g} | {c.response_signal:.6g} | {c.common:.6g} | "
            f"{c.tight:.6g} | {c.nontight:.6g} | {c.neither:.6g} | "
            f"{c.A:.6g} | {c.B:.6g} | {c.C:.6g} | {c.D:.6g} | "
            f"{safe_div(c.B, c.A):.6g} | {safe_div(c.C, c.A):.6g} | "
            f"{safe_div(c.D, c.A):.6g} | {safe_div(c.C + c.D, c.A):.6g} | "
            f"{safe_div(c.nontight, c.common):.6g} |"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
