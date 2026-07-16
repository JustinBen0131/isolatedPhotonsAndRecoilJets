#!/usr/bin/env python3
"""Build the THE-76 factorized Eiso trace ledger.

This is a local, read-only diagnostic artifact writer. It does not inspect or
mutate SDCC state. The goal is to separate event/cluster-population questions
from raw topo-isolation machinery questions for the PPG12 Fig.24/Fig.3
high-Eiso tail discrepancy.
"""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable


REPO = Path(__file__).resolve().parents[3]
OUTDIR = (
    REPO
    / "dataOutput/ppg12Parity/control_plane/audits"
    / "the76_eiso_factorized_trace_20260706"
)


@dataclass(frozen=True)
class OperationRow:
    operation: str
    ppg12_behavior: str
    ppg12_evidence: str
    recoiljets_behavior: str
    recoiljets_evidence: str
    same_different_unknown: str
    can_produce_high_eiso_tail: str
    expected_symptom_if_different: str
    next_proof: str


@dataclass(frozen=True)
class CanaryRow:
    field: str
    status: str
    evidence_or_blocker: str
    smallest_next_step: str


def write_table(base: Path, stem: str, rows: Iterable[object]) -> None:
    rows = list(rows)
    if not rows:
        raise ValueError(f"no rows for {stem}")
    dicts = [asdict(row) for row in rows]

    with (base / f"{stem}.json").open("w") as f:
        json.dump(dicts, f, indent=2, sort_keys=True)
        f.write("\n")

    with (base / f"{stem}.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(dicts[0].keys()))
        writer.writeheader()
        writer.writerows(dicts)

    with (base / f"{stem}.md").open("w") as f:
        headers = list(dicts[0].keys())
        f.write("| " + " | ".join(headers) + " |\n")
        f.write("| " + " | ".join(["---"] * len(headers)) + " |\n")
        for row in dicts:
            vals = [str(row[h]).replace("\n", "<br>").replace("|", "\\|") for h in headers]
            f.write("| " + " | ".join(vals) + " |\n")


def operation_rows() -> list[OperationRow]:
    return [
        OperationRow(
            operation="event population before isolation",
            ppg12_behavior=(
                "Anatreemaker writes a slimtree row only if saveevent is true: any "
                "stored cluster container has ncluster>0 or nparticles>0. Fig.24 then "
                "runs from that slimtree population."
            ),
            ppg12_evidence="ppg12codeGit/anatreemaker/source/CaloAna24.cc:2140-2153",
            recoiljets_behavior=(
                "Full RecoilJets event loop calls firstEventCuts before downstream "
                "analysis; for SIM it skips trigger logic but still applies |vz| if enabled."
            ),
            recoiljets_evidence="src/RecoilJets.cc:3672-3724 and src/RecoilJets.cc:4759-4764",
            same_different_unknown="different/unknown",
            can_produce_high_eiso_tail="yes",
            expected_symptom_if_different=(
                "Different vertex/event-activity population changes topo-cluster "
                "multiplicity and broadens or softens isolation before any raw-sum bug."
            ),
            next_proof=(
                "Foreground event-population canary: trigger/run/vz/nTopoClusters/sumET "
                "for the same photon20 source events before comparing raw Eiso."
            ),
        ),
        OperationRow(
            operation="candidate cluster selection",
            ppg12_behavior=(
                "Candidate clusters are read from the configured cluster container; "
                "E, eta, phi, ET are computed at m_vertex and clusters with ET < "
                "clusterpTmin are skipped."
            ),
            ppg12_evidence="ppg12codeGit/anatreemaker/source/CaloAna24.cc:1183-1203",
            recoiljets_behavior=(
                "PhotonClusterBuilder/RecoilJets use PhotonClusterv1 candidates and "
                "Fig.24 signal matching before filling h_singal_reco_isoET_0."
            ),
            recoiljets_evidence=(
                "src/RecoilJets.cc:9188-9208; src/RecoilJets.cc:13551-13558"
            ),
            same_different_unknown="unknown",
            can_produce_high_eiso_tail="yes",
            expected_symptom_if_different=(
                "A different set of signal clusters, especially at high pT, changes "
                "the sampled event environment and Eiso quantiles."
            ),
            next_proof=(
                "Canary must print cluster ET/eta/phi/truth-match before raw isolation."
            ),
        ),
        OperationRow(
            operation="topo-cluster collection source",
            ppg12_behavior="Uses RawClusterContainer node TOPOCLUSTER_ALLCALO.",
            ppg12_evidence=(
                "ppg12codeGit/anatreemaker/source/CaloAna24.cc:1320-1325 and "
                "ppg12codeGit/anatreemaker/source/CaloAna24.cc:2553-2560"
            ),
            recoiljets_behavior=(
                "Registers/uses TOPOCLUSTER_ALLCALO for PPG12 photon-yield topo isolation."
            ),
            recoiljets_evidence=(
                "macros/Fun4All_recoilJets_unified_impl.C:5469-5484 and "
                "macros/Fun4All_recoilJets_unified_impl.C:5513-5516"
            ),
            same_different_unknown="same name, content provenance unknown",
            can_produce_high_eiso_tail="yes",
            expected_symptom_if_different=(
                "Different topo-cluster builder settings or tower inputs directly alter "
                "neighbor ET sums and high-isolation tails."
            ),
            next_proof=(
                "Print topocluster count, ET spectrum, and neighbor list for same high-pT event."
            ),
        ),
        OperationRow(
            operation="candidate eta/phi/ET definition",
            ppg12_behavior=(
                "Uses RawClusterUtility eta/phi at vertex and ET = energy / cosh(eta)."
            ),
            ppg12_evidence="ppg12codeGit/anatreemaker/source/CaloAna24.cc:1190-1198",
            recoiljets_behavior=(
                "Stored path uses PhotonClusterBuilder iso seed and candidate ET; fallback "
                "uses stored iso axis or cluster eta/phi and cluster_et/cluster_pt."
            ),
            recoiljets_evidence=(
                "src/PhotonClusterBuilder.cc:2074-2075 and "
                "src/RecoilJets.cc:13602-13623"
            ),
            same_different_unknown="unknown for exact event/cluster",
            can_produce_high_eiso_tail="yes",
            expected_symptom_if_different=(
                "A shifted cone axis or candidate ET changes both cone membership and "
                "the final subtraction."
            ),
            next_proof=(
                "Canary print PPG12 candidate ET/eta/phi and RecoilJets axis/ET for same cluster."
            ),
        ),
        OperationRow(
            operation="neighbor topo-cluster collection",
            ppg12_behavior=(
                "Loops every topo cluster in TOPOCLUSTER_ALLCALO and computes eta/phi/ET at vertexz."
            ),
            ppg12_evidence="ppg12codeGit/anatreemaker/source/CaloAna24.cc:2558-2567",
            recoiljets_behavior=(
                "Stored path loops builder topocluster container at m_vertex; fallback loops "
                "m_ppg12TopoClusters at ppg12PhotonYieldKinematicVertexZ()."
            ),
            recoiljets_evidence=(
                "src/PhotonClusterBuilder.cc:2569-2588 and src/RecoilJets.cc:13630-13643"
            ),
            same_different_unknown="unknown",
            can_produce_high_eiso_tail="yes",
            expected_symptom_if_different=(
                "The same candidate can receive different neighbor ET when vertex or "
                "topocluster container differs."
            ),
            next_proof=(
                "High-pT cluster canary must print all neighbor topo clusters in R<0.4 "
                "with ET and dR for PPG12-style and RecoilJets-style sums."
            ),
        ),
        OperationRow(
            operation="DeltaR < 0.4 cone rule",
            ppg12_behavior="Uses sqrt(deta^2+dphi^2) with phi wrapping and dR < 0.4.",
            ppg12_evidence=(
                "ppg12codeGit/anatreemaker/source/CaloAna24.h:400-408 and "
                "ppg12codeGit/anatreemaker/source/CaloAna24.cc:2567-2570"
            ),
            recoiljets_behavior=(
                "Uses deltaR/dEta/dPhi wrapping and radius 0.4 in builder and fallback."
            ),
            recoiljets_evidence=(
                "src/PhotonClusterBuilder.h:174-177; src/PhotonClusterBuilder.cc:2593-2596; "
                "src/RecoilJets.cc:13645-13649"
            ),
            same_different_unknown="same",
            can_produce_high_eiso_tail="no unless axis/input differs",
            expected_symptom_if_different="Not a leading mechanism from source alone.",
            next_proof="No separate proof needed until axis/neighbors are aligned.",
        ),
        OperationRow(
            operation="candidate exclusion / candidate ET subtraction",
            ppg12_behavior=(
                "Adds candidate cluster if inside cone, then stores cluster_iso_topo_04 = topoET_04 - ET."
            ),
            ppg12_evidence=(
                "ppg12codeGit/anatreemaker/source/CaloAna24.cc:2568-2570 and "
                "ppg12codeGit/anatreemaker/source/CaloAna24.cc:1934-1939"
            ),
            recoiljets_behavior=(
                "Builder/fallback compute raw = topo_sum - candidate_et; optional extra "
                "exclude-candidate subtraction is controlled by a flag."
            ),
            recoiljets_evidence=(
                "src/PhotonClusterBuilder.cc:2599-2605; src/RecoilJets.cc:13651-13659; "
                "src/PhotonClusterBuilder.h:174-177"
            ),
            same_different_unknown="same if exclude flag false; verify runtime",
            can_produce_high_eiso_tail="yes if flag differs or candidate ET differs",
            expected_symptom_if_different=(
                "An extra subtraction would soften, not harden, Eiso; a missing subtraction "
                "or candidate ET mismatch could harden."
            ),
            next_proof="Canary print exclude flag and candidate ET used in raw sum.",
        ),
        OperationRow(
            operation="summed ET definition",
            ppg12_behavior="Uses topo ET = cluster energy / cosh(topo eta); no explicit ET<=0 skip in shown loop.",
            ppg12_evidence="ppg12codeGit/anatreemaker/source/CaloAna24.cc:2564-2567",
            recoiljets_behavior=(
                "Uses topo ET = energy / cosh(topo eta); RecoilJets skips nonfinite or ET<=0 neighbors."
            ),
            recoiljets_evidence="src/PhotonClusterBuilder.cc:2580-2591 and src/RecoilJets.cc:13638-13644",
            same_different_unknown="slightly different guard; likely small unless negative/zero topo clusters exist",
            can_produce_high_eiso_tail="possible but not leading",
            expected_symptom_if_different=(
                "Skipping nonpositive neighbors cannot create a harder positive tail; container "
                "content/axis differences remain more plausible."
            ),
            next_proof="Neighbor-list canary should include nonpositive/nonfinite neighbors if present.",
        ),
        OperationRow(
            operation="stored branch versus recompute policy",
            ppg12_behavior="Consumes one stored slimtree value cluster_iso_topo_04 in efficiency code.",
            ppg12_evidence="ppg12codeGit/efficiencytool/RecoEffCalculator_TTreeReader.C:2151-2155",
            recoiljets_behavior=(
                "Accepts stored ppg12_topo_raw_eiso_04 only if valid and vertex-compatible; otherwise "
                "silently recomputes from topoclusters."
            ),
            recoiljets_evidence="src/RecoilJets.cc:13560-13591 and src/RecoilJets.cc:13593-13668",
            same_different_unknown="different",
            can_produce_high_eiso_tail="yes",
            expected_symptom_if_different=(
                "Fallback subset can be tail-heavy if recompute uses a different vertex/axis/container "
                "than the builder or PPG12 slimtree value."
            ),
            next_proof=(
                "Diagnostic-only foreground canary: stored valid, stored vertex, expected vertex, "
                "stored raw, recomputed raw, final path choice, q70/q80/q90 by path."
            ),
        ),
        OperationRow(
            operation="linear Eiso correction",
            ppg12_behavior="For MC, corrected Eiso = 1.2 * raw + 0.1.",
            ppg12_evidence=(
                "ppg12codeGit/efficiencytool/config_bdt_nom.yaml:61,81-82 and "
                "ppg12codeGit/efficiencytool/RecoEffCalculator_TTreeReader.C:2261-2265"
            ),
            recoiljets_behavior="For pp SIM, corrected Eiso = m_scale * raw + m_shift.",
            recoiljets_evidence="src/RecoilJets.cc:13671-13683 and src/RecoilJets.h:1621-1622",
            same_different_unknown="same by source/default",
            can_produce_high_eiso_tail="no, mostly ruled out",
            expected_symptom_if_different="A scale/offset mismatch would move all pT bins coherently.",
            next_proof="Already secondary; keep in validation table but do not lead with it.",
        ),
        OperationRow(
            operation="Fig24 fill population",
            ppg12_behavior=(
                "Signal fill requires cluster truth-track id to map to photon_reco; fills "
                "h_singal_reco_isoET with cluster_Et and corrected recoisoET."
            ),
            ppg12_evidence="ppg12codeGit/efficiencytool/RecoEffCalculator_TTreeReader.C:2732-2825 and 2916-2917",
            recoiljets_behavior=(
                "Fills h_singal_reco_isoET_0 with matched reco cluster and corrected "
                "ppg12PhotonYieldEiso(ppg12PhotonYieldRawEiso(...))."
            ),
            recoiljets_evidence="src/RecoilJets.cc:9188-9208",
            same_different_unknown="unknown",
            can_produce_high_eiso_tail="yes",
            expected_symptom_if_different=(
                "Even identical raw-isolation machinery gives different distributions if the "
                "signal-cluster set differs."
            ),
            next_proof=(
                "Before raw-Eiso identity, canary must print selected cluster/truth-match set "
                "and event activity."
            ),
        ),
    ]


def canary_rows() -> list[CanaryRow]:
    return [
        CanaryRow(
            field="local final ROOT per-cluster carrier",
            status="blocked in current final ROOT",
            evidence_or_blocker=(
                "Current merged ROOT contains histogram directories, but no TTree/event carrier "
                "with cluster IDs, neighbor lists, stored/recomputed path labels, or raw per-cluster values."
            ),
            smallest_next_step=(
                "Run a diagnostic-only foreground event canary, not a final-histogram audit."
            ),
        ),
        CanaryRow(
            field="event population canary",
            status="required before raw isolation comparison",
            evidence_or_blocker=(
                "Fig.3 data can be biased by trigger/run/vz/event-quality/activity; Fig.24 MC can "
                "be biased by vz/event-saving/source-composition/activity. These are upstream of Eiso."
            ),
            smallest_next_step=(
                "For bounded photon20 30-36 GeV events, print trigger/run/vz/nTopoClusters/sumET "
                "and selected cluster count before raw Eiso."
            ),
        ),
        CanaryRow(
            field="cluster selection canary",
            status="required before raw isolation comparison",
            evidence_or_blocker=(
                "Different high-pT signal cluster set can reproduce a high-Eiso tail without "
                "any isolation formula difference."
            ),
            smallest_next_step=(
                "Print cluster ET/eta/phi, truth-match ID/status, and whether each side selects it."
            ),
        ),
        CanaryRow(
            field="neighbor-list raw isolation canary",
            status="required after population alignment",
            evidence_or_blocker=(
                "The first raw-number operation that can differ for an equivalent cluster is the "
                "neighbor-topocluster collection and ET sum in R<0.4 under the chosen vertex/path."
            ),
            smallest_next_step=(
                "Print all neighbors in R<0.4, neighbor ET, dR, candidate exclusion decision, "
                "raw sum, corrected Eiso, stored/recomputed path choice."
            ),
        ),
    ]


def write_summary(base: Path) -> None:
    summary = {
        "campaign": "THE-76 | PPG12 PhotonJet SIM Parity",
        "target": "Fig24/Fig3 high-Eiso tail factorized trace",
        "classification": "blocked: source trace identifies first unproven operation, but per-cluster diagnostic labels are missing",
        "upstream_fork": (
            "Need decide whether the tail is from event/cluster population or raw Eiso "
            "for equivalent clusters."
        ),
        "first_raw_isolation_operation_that_can_change_value": (
            "neighbor topo-cluster collection and summed ET in R<0.4, evaluated with the "
            "selected candidate axis, vertex, topocluster container, and stored-vs-recompute path"
        ),
        "current_first_factorized_unknown": (
            "event and cluster population equivalence before isolation"
        ),
        "canonicalization_allowed": False,
        "minimal_next_diagnostic": (
            "diagnostic-only foreground photon20 30-36 GeV canary printing event activity, "
            "selected cluster identity, all R<0.4 neighbor topo clusters, stored raw Eiso, "
            "recomputed raw Eiso, final path choice, and corrected Eiso"
        ),
    }
    with (base / "factorized_trace_summary.json").open("w") as f:
        json.dump(summary, f, indent=2, sort_keys=True)
        f.write("\n")
    with (base / "factorized_trace_summary.md").open("w") as f:
        f.write("# THE-76 Fig24/Fig3 Eiso Factorized Trace Summary\n\n")
        for key, value in summary.items():
            f.write(f"- **{key}**: {value}\n")


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    write_table(OUTDIR, "event_cluster_isolation_factorization_table", operation_rows())
    write_table(OUTDIR, "high_pt_cluster_canary_requirements", canary_rows())
    write_summary(OUTDIR)
    index = {
        "audit_dir": str(OUTDIR),
        "tables": [
            "event_cluster_isolation_factorization_table.md",
            "event_cluster_isolation_factorization_table.csv",
            "event_cluster_isolation_factorization_table.json",
            "high_pt_cluster_canary_requirements.md",
            "high_pt_cluster_canary_requirements.csv",
            "high_pt_cluster_canary_requirements.json",
            "factorized_trace_summary.md",
            "factorized_trace_summary.json",
            "local_root_event_carrier_check.txt",
        ],
    }
    with (OUTDIR / "audit_index.json").open("w") as f:
        json.dump(index, f, indent=2, sort_keys=True)
        f.write("\n")
    print(OUTDIR)


if __name__ == "__main__":
    main()
