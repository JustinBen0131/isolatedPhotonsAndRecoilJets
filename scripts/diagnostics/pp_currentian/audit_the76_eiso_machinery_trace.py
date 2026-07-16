#!/usr/bin/env python3
"""Write THE-76 Fig24/Fig3 Eiso machinery trace artifacts.

This is a source/procedure ledger, not a production runner.  It records the
first isolation-machinery operation that is not yet proven equivalent between
PPG12 and RecoilJets, and defines the smallest diagnostic canary needed next.
"""

from __future__ import annotations

import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


REPO = Path(__file__).resolve().parents[3]
AUDIT_DIR = (
    REPO
    / "dataOutput/ppg12Parity/control_plane/audits/the76_eiso_machinery_trace_20260706"
)


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def md_table(rows: list[dict[str, str]]) -> str:
    if not rows:
        return ""
    cols = list(rows[0].keys())
    lines = [
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join(["---"] * len(cols)) + " |",
    ]
    for row in rows:
        vals = [str(row.get(c, "")).replace("\n", "<br>").replace("|", "\\|") for c in cols]
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines) + "\n"


def write_md(path: Path, title: str, rows: list[dict[str, str]], intro: Iterable[str] = ()) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    body = [f"# {title}", ""]
    body.extend(intro)
    if intro:
        body.append("")
    body.append(md_table(rows))
    path.write_text("\n".join(body), encoding="utf-8")


def rel(path: str) -> str:
    return path


def main() -> int:
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    generated_at = datetime.now(timezone.utc).isoformat()

    ppg12_rows = [
        {
            "stage": "topo node lookup",
            "PPG12 operation": "anatreemaker reads TOPOCLUSTER_ALLCALO and TOPOCLUSTER_ALLCALO_SOFT from the node tree",
            "exact source path/line": rel("ppg12codeGit/anatreemaker/source/CaloAna24.cc:1115-1120"),
            "input": "topNode",
            "output": "topoClusterContainer, topoClusterContainer_soft",
            "cut/condition": "if missing, topo iso is left as zero/warned",
            "can affect high-Eiso tail yes/no": "yes",
        },
        {
            "stage": "candidate axis and ET",
            "PPG12 operation": "uses RawClusterUtility eta/phi at m_vertex and candidate ET = cluster energy / cosh(eta)",
            "exact source path/line": rel("ppg12codeGit/anatreemaker/source/CaloAna24.cc:1190-1203"),
            "input": "recoCluster, m_vertex",
            "output": "eta, phi, ET",
            "cut/condition": "requires ET >= clusterpTmin",
            "can affect high-Eiso tail yes/no": "yes",
        },
        {
            "stage": "topo cone sum",
            "PPG12 operation": "calculateET_topo_6cones sums topo-cluster ET for dR thresholds including R < 0.4",
            "exact source path/line": rel("ppg12codeGit/anatreemaker/source/CaloAna24.cc:1320-1325,2553-2589"),
            "input": "iso_axis_eta, iso_axis_phi, TOPOCLUSTER_ALLCALO",
            "output": "topoET_04",
            "cut/condition": "RawClusterUtility eta/phi at vertexz; topo ET = energy / cosh(eta); dR with phi wrapping",
            "can affect high-Eiso tail yes/no": "yes",
        },
        {
            "stage": "raw branch construction",
            "PPG12 operation": "stores raw topo isolation as cluster_iso_topo_04 = topoET_04 - ET",
            "exact source path/line": rel("ppg12codeGit/anatreemaker/source/CaloAna24.cc:1934-1939"),
            "input": "topoET_04, candidate ET",
            "output": "cluster_iso_topo_04",
            "cut/condition": "candidate cluster is not separately excluded by ID; one candidate ET scalar is subtracted",
            "can affect high-Eiso tail yes/no": "yes",
        },
        {
            "stage": "branch write",
            "PPG12 operation": "writes cluster_iso_topo_04_<clusterNode> to slimtree",
            "exact source path/line": rel("ppg12codeGit/anatreemaker/source/CaloAna24.cc:256-261"),
            "input": "cluster_iso_topo_04 array",
            "output": "slimtree branch cluster_iso_topo_04_*",
            "cut/condition": "per-cluster branch in saved slimtree entry",
            "can affect high-Eiso tail yes/no": "yes",
        },
        {
            "stage": "Fig24/Fig3 branch read",
            "PPG12 operation": "efficiency tool reads cluster_iso_topo_04_<clusterNode> with TTreeReaderArray",
            "exact source path/line": rel("ppg12codeGit/efficiencytool/RecoEffCalculator_TTreeReader.C:689-690"),
            "input": "slimtree branch",
            "output": "cluster_iso_topo_04[icluster]",
            "cut/condition": "selected by use_topo_iso = 2",
            "can affect high-Eiso tail yes/no": "yes",
        },
        {
            "stage": "raw-to-corrected convention",
            "PPG12 operation": "for nominal MC, recoisoET = cluster_iso_topo_04 then recoisoET = recoisoET * 1.2 + 0.1",
            "exact source path/line": rel("ppg12codeGit/efficiencytool/config_bdt_nom.yaml:61,81-82; ppg12codeGit/efficiencytool/RecoEffCalculator_TTreeReader.C:2152-2154,2261-2265"),
            "input": "cluster_iso_topo_04",
            "output": "corrected recoisoET",
            "cut/condition": "MC only for scale/shift; data keeps the selected branch convention without MC scale/shift",
            "can affect high-Eiso tail yes/no": "yes",
        },
        {
            "stage": "signal fill",
            "PPG12 operation": "fills h_singal_reco_isoET_0 with cluster_Et and corrected recoisoET for truth-associated signal clusters",
            "exact source path/line": rel("ppg12codeGit/efficiencytool/RecoEffCalculator_TTreeReader.C:1213-1215,2732-2825,2916"),
            "input": "cluster_Et, recoisoET, truth association by cluster_truthtrkID/photon_reco",
            "output": "h_singal_reco_isoET_0",
            "cut/condition": "eta bin must be valid; truth track must map to photon_reco for signal fill",
            "can affect high-Eiso tail yes/no": "yes",
        },
        {
            "stage": "cutoff extraction",
            "PPG12 operation": "FindETCut.C integrates h_singal_reco_isoET_0 to obtain 70/80/90 percent cutoffs",
            "exact source path/line": rel("ppg12codeGit/efficiencytool/FindETCut.C:86-105"),
            "input": "h_singal_reco_isoET_0",
            "output": "cutoff points",
            "cut/condition": "previous audits ruled out cutoff extraction as the source of the filled-shape discrepancy",
            "can affect high-Eiso tail yes/no": "no",
        },
    ]

    rj_rows = [
        {
            "stage": "builder config",
            "RecoilJets operation": "Fun4All config enables PPG12 topo isolation on PhotonClusterBuilder with TOPOCLUSTER_ALLCALO and R = 0.4",
            "exact source path/line": rel("macros/Fun4All_recoilJets_unified_impl.C:5470-5484"),
            "input": "macro config flags",
            "output": "PhotonClusterBuilder PPG12 topo settings",
            "cut/condition": "candidate-exclusion flag controlled by ppg12ExcludeCandidateTopo",
            "can affect high-Eiso tail yes/no": "yes",
        },
        {
            "stage": "topo node lookup",
            "RecoilJets operation": "PhotonClusterBuilder resolves m_ppg12_topocluster_container from the configured topo node",
            "exact source path/line": rel("src/PhotonClusterBuilder.cc:309-315,725-727; src/PhotonClusterBuilder.h:174-177"),
            "input": "topNode, m_ppg12_topocluster_node",
            "output": "m_ppg12_topocluster_container",
            "cut/condition": "missing container disables valid stored branch",
            "can affect high-Eiso tail yes/no": "yes",
        },
        {
            "stage": "builder vertex",
            "RecoilJets operation": "PhotonClusterBuilder chooses m_vertex from truth/global-MBD/MBD/global vertex configuration and stores it with the PPG12 topo branch",
            "exact source path/line": rel("src/PhotonClusterBuilder.cc:746-832,2339-2345"),
            "input": "vertex nodes/config",
            "output": "m_vertex and ppg12_topo_vertex_z",
            "cut/condition": "stored branch can later be rejected if RecoilJets expected vertex differs by >= 1e-3",
            "can affect high-Eiso tail yes/no": "yes",
        },
        {
            "stage": "stored raw construction",
            "RecoilJets operation": "PhotonClusterBuilder computes ppg12_topo_raw_eiso_04 = topo sum ET - candidate ET",
            "exact source path/line": rel("src/PhotonClusterBuilder.cc:2069-2080,2554-2605"),
            "input": "iso_seed_eta/phi, candidate ET, m_ppg12_topocluster_container",
            "output": "ppg12_topo_raw_eiso_04, ppg12_topo_valid_04",
            "cut/condition": "finite inputs, candidate ET > 0, topo ET > 0, dR < radius; optional candidate extra subtraction if exclude flag is on",
            "can affect high-Eiso tail yes/no": "yes",
        },
        {
            "stage": "stored branch write",
            "RecoilJets operation": "PhotonClusterBuilder writes ppg12_topo_raw_eiso_04, sumET, valid flag, radius, vertex, and exclude-candidate flag as shower-shape parameters",
            "exact source path/line": rel("src/PhotonClusterBuilder.cc:2337-2345"),
            "input": "computed PPG12 topo values",
            "output": "PhotonClusterv1 shower-shape parameters",
            "cut/condition": "only when m_use_ppg12_topocluster_isolation is enabled",
            "can affect high-Eiso tail yes/no": "yes",
        },
        {
            "stage": "stored branch use",
            "RecoilJets operation": "RecoilJets first tries the stored ppg12_topo_raw_eiso_04 if valid and vertex-compatible",
            "exact source path/line": rel("src/RecoilJets.cc:13487-13525"),
            "input": "PhotonClusterv1 stored raw Eiso, stored vertex, expected vertex",
            "output": "stored raw Eiso",
            "cut/condition": "stored valid > 0.5, raw finite, raw < 1e8, vertex mismatch < 1e-3 or expected vertex nonfinite",
            "can affect high-Eiso tail yes/no": "yes",
        },
        {
            "stage": "fallback recompute",
            "RecoilJets operation": "if stored branch is invalid or vertex-incompatible, RecoilJets recomputes topo-cone isolation from m_ppg12TopoClusters",
            "exact source path/line": rel("src/RecoilJets.cc:13526-13611"),
            "input": "axis eta/phi, candidate ET, m_ppg12TopoClusters, ppg12PhotonYieldKinematicVertexZ()",
            "output": "raw iso = topoEt04 - candidateEt, optional extra candidate subtraction",
            "cut/condition": "silent fallback path is not persisted in the merged ROOT",
            "can affect high-Eiso tail yes/no": "yes",
        },
        {
            "stage": "raw-to-corrected convention",
            "RecoilJets operation": "ppg12PhotonYieldEiso applies m_ppg12PhotonYieldMcIsoScale * raw + m_ppg12PhotonYieldMcIsoShift for pp SIM",
            "exact source path/line": rel("src/RecoilJets.cc:13614-13626; src/RecoilJets.h:1621-1622"),
            "input": "stored or recomputed raw Eiso",
            "output": "corrected Eiso",
            "cut/condition": "scale=1.2 shift=0.1 by default; prior audit did not find correction as leading mismatch",
            "can affect high-Eiso tail yes/no": "yes",
        },
        {
            "stage": "Fig24 fill",
            "RecoilJets operation": "fills SIM/h_singal_reco_isoET_0 and named Fig24 variants with rPt and corrected Eiso for matched signal clusters",
            "exact source path/line": rel("src/RecoilJets.cc:9157,9180-9244"),
            "input": "recoMatch, corrected Eiso, sample ownership, active triggers",
            "output": "SIM/h_singal_reco_isoET_0",
            "cut/condition": "requires eta/pT bins, finite Eiso, photon sample context/window/owner truth",
            "can affect high-Eiso tail yes/no": "yes",
        },
        {
            "stage": "Fig3 data path",
            "RecoilJets operation": "data-side isolation calls the same ppg12PhotonYieldRawEiso/ppg12PhotonYieldEiso helpers when PPG12 yield iso is enabled",
            "exact source path/line": rel("src/RecoilJets.cc:8143-8148,10812-10818"),
            "input": "data PhotonClusterv1 / topo containers",
            "output": "current data Eiso values",
            "cut/condition": "Fig3 data tail sharing the discrepancy points to the common isolation path, not MC-only weighting",
            "can affect high-Eiso tail yes/no": "yes",
        },
    ]

    comparison_rows = [
        {
            "operation": "cluster collection source",
            "PPG12 behavior": "TOPOCLUSTER_ALLCALO resolved in CaloAna24",
            "RecoilJets behavior": "TOPOCLUSTER_ALLCALO configured for PhotonClusterBuilder and used by RecoilJets fallback",
            "same/different/unknown": "same by name, not proven same event content",
            "exact evidence": "CaloAna24.cc:1115; Fun4All_recoilJets_unified_impl.C:5482; PhotonClusterBuilder.h:174",
            "expected symptom if different": "tail or global shift in raw Eiso",
            "next proof": "same-event/topo-container sum comparison",
        },
        {
            "operation": "topo-cluster definition",
            "PPG12 behavior": "consumes existing TOPOCLUSTER_ALLCALO; producer settings not audited here",
            "RecoilJets behavior": "consumes existing TOPOCLUSTER_ALLCALO; same producer assumed for stored/fallback",
            "same/different/unknown": "unknown",
            "exact evidence": "both source paths consume the node, not producer provenance",
            "expected symptom if different": "broad Eiso tail in both data and MC",
            "next proof": "topo node provenance and same-event topo list checksum",
        },
        {
            "operation": "candidate cluster exclusion",
            "PPG12 behavior": "cluster_iso_topo_04 = topoET_04 - ET; no separate topo-cluster ID exclusion visible",
            "RecoilJets behavior": "stored branch does topo_sum - candidate_et; fallback can optionally subtract candidate ET a second time if exclude flag is on",
            "same/different/unknown": "same if exclude flag off; different if on",
            "exact evidence": "CaloAna24.cc:1939; PhotonClusterBuilder.cc:2600-2604; RecoilJets.cc:13594-13602",
            "expected symptom if different": "large Eiso scale offset, likely not just high-pT tail",
            "next proof": "record ppg12_topo_exclude_candidate and RecoilJets flag per run",
        },
        {
            "operation": "cone radius",
            "PPG12 behavior": "R < 0.4 for cluster_iso_topo_04",
            "RecoilJets behavior": "R=0.4 configured for stored branch; fallback uses kPPG12YieldRecoIsoConeR=0.4",
            "same/different/unknown": "same",
            "exact evidence": "CaloAna24.cc:2568-2570; Fun4All_recoilJets_unified_impl.C:5483; RecoilJets.cc:378,13591",
            "expected symptom if different": "tail and cutoff mismatch",
            "next proof": "none unless runtime config differs",
        },
        {
            "operation": "eta/phi distance metric",
            "PPG12 behavior": "RawClusterUtility eta/phi at vertexz, deltaR with phi wrap",
            "RecoilJets behavior": "RawClusterUtility eta/phi at selected vertex; TVector2 phi wrap in fallback; builder deltaR wrap",
            "same/different/unknown": "formula same; vertex choice/path unknown per cluster",
            "exact evidence": "CaloAna24.cc:2558-2567; CaloAna24.h:400-408; PhotonClusterBuilder.cc:2569-2595; RecoilJets.cc:13573-13591",
            "expected symptom if different": "fallback-only high-pT tail if vertex compatibility pushes clusters into different vertex path",
            "next proof": "record stored vertex, expected vertex, and path choice",
        },
        {
            "operation": "cluster/tower energy variable",
            "PPG12 behavior": "topo ET = cluster energy / cosh(topo eta); candidate ET = reco energy / cosh(candidate eta)",
            "RecoilJets behavior": "same ET formula for topo and candidate ET from PhotonClusterv1 shower-shape parameter",
            "same/different/unknown": "same formula; candidate source not proven same",
            "exact evidence": "CaloAna24.cc:1193-1198,2564-2567; PhotonClusterBuilder.cc:2587; RecoilJets.cc:13560-13565,13585",
            "expected symptom if different": "raw Eiso shift proportional to cluster ET",
            "next proof": "same-cluster candidate ET comparison",
        },
        {
            "operation": "neighbor object threshold",
            "PPG12 behavior": "topo-cluster loop excludes nonpositive/nonexistent only implicitly through topo clusters; no explicit topo ET threshold in calculateET_topo_6cones",
            "RecoilJets behavior": "builder/fallback skip nonfinite or topo ET <= 0",
            "same/different/unknown": "probably same for positive topo clusters, not proven",
            "exact evidence": "CaloAna24.cc:2559-2589; PhotonClusterBuilder.cc:2587-2591; RecoilJets.cc:13585-13587",
            "expected symptom if different": "small tail/normalization differences",
            "next proof": "topo list and summed ET canary",
        },
        {
            "operation": "raw isolation branch",
            "PPG12 behavior": "single branch cluster_iso_topo_04 from slimtree",
            "RecoilJets behavior": "stored ppg12_topo_raw_eiso_04 if valid, otherwise recomputed raw Eiso",
            "same/different/unknown": "different path policy",
            "exact evidence": "RecoEffCalculator_TTreeReader.C:689-690,2152-2154; RecoilJets.cc:13503-13611",
            "expected symptom if different": "mixed distribution can be broader/tail-heavy if fallback path differs",
            "next proof": "path-split canary",
        },
        {
            "operation": "validity checks",
            "PPG12 behavior": "branch is consumed directly once event/cluster exists",
            "RecoilJets behavior": "stored branch requires valid flag, finite raw, finite vertex compatibility",
            "same/different/unknown": "different",
            "exact evidence": "RecoEffCalculator_TTreeReader.C:2152-2154; RecoilJets.cc:13503-13525",
            "expected symptom if different": "clusters leave PPG12 branch path and enter recompute path",
            "next proof": "stored_valid and vertex-compatible counters vs pT",
        },
        {
            "operation": "vertex compatibility",
            "PPG12 behavior": "topo branch constructed at CaloAna vertexz and then consumed without a later vertex compatibility gate",
            "RecoilJets behavior": "stored branch rejected if stored vertex and expected vertex differ by >= 1e-3",
            "same/different/unknown": "different",
            "exact evidence": "CaloAna24.cc:2558; RecoilJets.cc:13510-13515",
            "expected symptom if different": "high-pT clusters may silently use recomputed isolation under a different vertex",
            "next proof": "per-cluster vertex compatibility/path-choice canary",
        },
        {
            "operation": "fallback/recompute path",
            "PPG12 behavior": "no equivalent fallback visible in Fig24/Fig3 consumer",
            "RecoilJets behavior": "recomputes from m_ppg12TopoClusters if stored value fails checks",
            "same/different/unknown": "different",
            "exact evidence": "RecoilJets.cc:13526-13611",
            "expected symptom if different": "tail-heavy final mixed Eiso if recompute differs from slimtree branch",
            "next proof": "stored-only vs recomputed-only q70/q80/q90",
        },
        {
            "operation": "correction formula",
            "PPG12 behavior": "MC nominal applies 1.2*raw + 0.1",
            "RecoilJets behavior": "pp SIM applies 1.2*raw + 0.1",
            "same/different/unknown": "same",
            "exact evidence": "config_bdt_nom.yaml:81-82; RecoilJets.h:1621-1622; RecoilJets.cc:13614-13626",
            "expected symptom if different": "global scale/offset, already not leading",
            "next proof": "none unless runtime env overrides scale/shift",
        },
        {
            "operation": "cluster pT definition",
            "PPG12 behavior": "uses cluster_Et from anatreemaker",
            "RecoilJets behavior": "uses rPt/recoMatch cluster pT/ET in Fig24 fill",
            "same/different/unknown": "unknown",
            "exact evidence": "RecoEffCalculator_TTreeReader.C:2283-2288,2916; RecoilJets.cc:9180-9188,9230",
            "expected symptom if different": "pT-bin migration, especially high pT",
            "next proof": "same-cluster pT/eta canary",
        },
        {
            "operation": "truth-match/signal selection",
            "PPG12 behavior": "requires cluster truth track maps to a particle in photon_reco; dR cut removed/commented",
            "RecoilJets behavior": "uses recoMatch and sample owner truth context before Fig24 fill",
            "same/different/unknown": "similar intent, not same-cluster proven",
            "exact evidence": "RecoEffCalculator_TTreeReader.C:2732-2825,2916; RecoilJets.cc:9157,9180-9221",
            "expected symptom if different": "selected cluster set differs, can change tails",
            "next proof": "same-cluster truth-match canary",
        },
        {
            "operation": "Fig24 fill object",
            "PPG12 behavior": "h_singal_reco_isoET_0",
            "RecoilJets behavior": "SIM/h_singal_reco_isoET_0 compatibility object and named variants",
            "same/different/unknown": "same target shape, different upstream machinery",
            "exact evidence": "RecoEffCalculator_TTreeReader.C:1213-1215,2916; RecoilJets.cc:9227-9244",
            "expected symptom if different": "not causal by itself",
            "next proof": "source path canary",
        },
        {
            "operation": "Fig3 fill object",
            "PPG12 behavior": "efficiency-tool data isolation templates use the same topo branch convention without MC scale/shift",
            "RecoilJets behavior": "data path invokes the same PPG12 yield iso helpers",
            "same/different/unknown": "same helper family but not same raw path proven",
            "exact evidence": "RecoEffCalculator_TTreeReader.C:2152-2154; RecoilJets.cc:8143-8148,10812-10818",
            "expected symptom if different": "same high-isolation tail in data and MC",
            "next proof": "data stored-vs-recompute counters",
        },
    ]

    canary_rows = [
        {
            "diagnostic": "stored-vs-recomputed path split",
            "scope": "photon20, 30-36 GeV Fig24-matched signal clusters",
            "status": "not runnable from final merged ROOT",
            "exact blocker": "current final ROOT does not store ppg12_topo_raw_eiso_04, stored-valid flag, stored vertex, expected vertex, recomputed raw Eiso, or path choice per cluster",
            "smallest instrumentation": "diagnostic-only foreground canary around RecoilJets::ppg12PhotonYieldRawEiso that records stored_valid, vertex_compatible, stored_raw, recomputed_raw, final_raw, path_choice, corrected_eiso, rPt, rEta, truth-match, pT/Eiso bin, and weight",
            "decision enabled": "stored-only vs recomputed-only q70/q80/q90 and fallback fraction vs pT",
        },
        {
            "diagnostic": "same-cluster PPG12 comparison",
            "scope": "matched photon20 high-pT clusters in PPG12 slimtree and RecoilJets event",
            "status": "blocked until event/cluster keys and RecoilJets per-cluster raw-path outputs exist",
            "exact blocker": "PPG12 has cluster_iso_topo_04 in slimtree, but current RecoilJets products do not persist comparable per-cluster stored/recomputed raw Eiso or path labels",
            "smallest instrumentation": "same foreground canary plus event/run/source id, cluster id if available, truth track id, cluster pT/eta/phi; then compare to PPG12 slimtree row",
            "decision enabled": "distinguish stored branch mismatch from recompute mismatch from truth-match/selection mismatch",
        },
        {
            "diagnostic": "Fig3 data linkage",
            "scope": "data clusters entering current Fig3 isolation overlay",
            "status": "blocked by same missing path labels",
            "exact blocker": "data final histograms show the high tail but do not persist isolation path choice",
            "smallest instrumentation": "same path-choice counters for data candidates, without MC truth columns",
            "decision enabled": "if data and MC tails are driven by the same recompute path, isolation machinery is favored over MC weighting",
        },
    ]

    summary = {
        "generated_at": generated_at,
        "current_state": (
            "Fig24/Fig3 high-Eiso discrepancy is real in filled Eiso distributions; "
            "final merged ROOTs do not persist per-cluster isolation path labels."
        ),
        "first_different_or_unproven_operation": (
            "raw isolation branch/path policy: PPG12 consumes one slimtree branch "
            "cluster_iso_topo_04, while RecoilJets uses PhotonClusterBuilder stored "
            "ppg12_topo_raw_eiso_04 only if valid and vertex-compatible, otherwise "
            "silently recomputes topo-cone isolation."
        ),
        "most_likely_controllable_mechanism": (
            "recomputed topo-cone fallback or vertex-compatibility path selection "
            "makes the final RecoilJets Eiso distribution broader than the PPG12 "
            "cluster_iso_topo_04 distribution."
        ),
        "classification": "blocked: exact source-level operation found but diagnostic labels missing",
        "canonicalization_allowed": False,
        "exact_fix_or_rerun_needed": (
            "Add diagnostic-only foreground instrumentation around "
            "RecoilJets::ppg12PhotonYieldRawEiso and, if needed, PhotonClusterBuilder "
            "to record stored raw, recomputed raw, vertex compatibility, and path choice "
            "for bounded photon20 30-36 GeV signal clusters. Do not broad-rerun until "
            "the path split identifies whether the fix is stored-branch construction, "
            "vertex compatibility, fallback recomputation, or cluster selection."
        ),
        "validation_target_after_fix": (
            "Path-split q70/q80/q90 must identify the tail source; after the controlled "
            "fix, Fig24 h_singal_reco_isoET_0 cutoff ratios should no longer show the "
            "34-36 GeV SDCC/current drop, and Fig3 data tail should improve if the same "
            "path is causal."
        ),
    }

    artifacts = {
        "ppg12_ledger_csv": AUDIT_DIR / "ppg12_isolation_machinery_ledger.csv",
        "ppg12_ledger_json": AUDIT_DIR / "ppg12_isolation_machinery_ledger.json",
        "ppg12_ledger_md": AUDIT_DIR / "ppg12_isolation_machinery_ledger.md",
        "recoiljets_ledger_csv": AUDIT_DIR / "recoiljets_isolation_machinery_ledger.csv",
        "recoiljets_ledger_json": AUDIT_DIR / "recoiljets_isolation_machinery_ledger.json",
        "recoiljets_ledger_md": AUDIT_DIR / "recoiljets_isolation_machinery_ledger.md",
        "equivalence_csv": AUDIT_DIR / "operation_equivalence_table.csv",
        "equivalence_json": AUDIT_DIR / "operation_equivalence_table.json",
        "equivalence_md": AUDIT_DIR / "operation_equivalence_table.md",
        "canary_csv": AUDIT_DIR / "diagnostic_canary_plan.csv",
        "canary_json": AUDIT_DIR / "diagnostic_canary_plan.json",
        "canary_md": AUDIT_DIR / "diagnostic_canary_plan.md",
        "summary_json": AUDIT_DIR / "eiso_machinery_trace_summary.json",
        "summary_md": AUDIT_DIR / "eiso_machinery_trace_summary.md",
        "audit_index": AUDIT_DIR / "audit_index.json",
    }

    write_csv(artifacts["ppg12_ledger_csv"], ppg12_rows)
    write_json(artifacts["ppg12_ledger_json"], {"rows": ppg12_rows, "generated_at": generated_at})
    write_md(
        artifacts["ppg12_ledger_md"],
        "THE-76 PPG12 Isolation Machinery Ledger",
        ppg12_rows,
        [
            "This ledger traces the PPG12 source path for Fig24/Fig3 Eiso.",
            "The key raw observable is `cluster_iso_topo_04`, then MC Fig24 applies `1.2 * raw + 0.1`.",
        ],
    )

    write_csv(artifacts["recoiljets_ledger_csv"], rj_rows)
    write_json(artifacts["recoiljets_ledger_json"], {"rows": rj_rows, "generated_at": generated_at})
    write_md(
        artifacts["recoiljets_ledger_md"],
        "THE-76 RecoilJets Isolation Machinery Ledger",
        rj_rows,
        [
            "This ledger traces the current RecoilJets/PhotonClusterBuilder source path for Fig24/Fig3 Eiso.",
            "The critical policy difference is stored-branch use with vertex compatibility followed by silent recompute fallback.",
        ],
    )

    write_csv(artifacts["equivalence_csv"], comparison_rows)
    write_json(artifacts["equivalence_json"], {"rows": comparison_rows, "generated_at": generated_at})
    write_md(
        artifacts["equivalence_md"],
        "THE-76 Eiso Operation Equivalence Table",
        comparison_rows,
        [
            "Rows marked different or unknown are not automatically causal.",
            "The first operation both different/unproven and capable of the observed tail is the raw isolation branch/path policy.",
        ],
    )

    write_csv(artifacts["canary_csv"], canary_rows)
    write_json(artifacts["canary_json"], {"rows": canary_rows, "generated_at": generated_at})
    write_md(
        artifacts["canary_md"],
        "THE-76 Eiso Diagnostic Canary Plan",
        canary_rows,
        [
            "No production or broad rerun is part of this canary.",
            "The current final ROOT cannot run this canary because it lacks per-cluster path labels.",
        ],
    )

    write_json(artifacts["summary_json"], summary)
    summary_md = f"""# THE-76 Eiso Machinery Trace Summary

- Generated: {generated_at}
- Current state: {summary["current_state"]}
- First different/unproven operation: {summary["first_different_or_unproven_operation"]}
- Most likely controllable mechanism: {summary["most_likely_controllable_mechanism"]}
- Classification: `{summary["classification"]}`
- Canonicalization allowed now: `{summary["canonicalization_allowed"]}`

## Exact Fix / Rerun Needed

{summary["exact_fix_or_rerun_needed"]}

## Validation Target After Fix

{summary["validation_target_after_fix"]}
"""
    artifacts["summary_md"].write_text(summary_md, encoding="utf-8")

    index_payload = {
        "generated_at": generated_at,
        "audit_dir": str(AUDIT_DIR),
        "artifacts": {k: str(v) for k, v in artifacts.items() if k != "audit_index"},
        "classification": summary["classification"],
        "first_different_or_unproven_operation": summary["first_different_or_unproven_operation"],
    }
    write_json(artifacts["audit_index"], index_payload)

    print(str(AUDIT_DIR))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
