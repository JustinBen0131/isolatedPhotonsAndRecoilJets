#!/usr/bin/env python3
"""THE-76 photon+jet Fig.5 slimTree-equivalence audit.

This is a read-only local diagnostic. It does not run production, submit
Condor, merge, transfer, mutate Slides, or canonicalize any contract.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

try:
    import uproot
except ModuleNotFoundError as exc:  # pragma: no cover
    raise SystemExit(
        "uproot is required; use /Users/patsfan753/Desktop/analysis/env/bin/python3"
    ) from exc


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
OUT_DIR = (
    REPO
    / "dataOutput/ppg12Parity/control_plane/audits/"
    / "the76_photonjet_slimtree_equivalence_20260702"
)

CURRENT_ROOT = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_fig24_photonjet_fix_20260702_014217/"
    / "final_roots/photonjet/RecoilJets_photonjet5plus10plus20_MERGED.root"
)
CURRENT_POINTS_CSV = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_fig24_photonjet_fix_20260702_014217/"
    / "strict_stitched_photonjet_current/"
    / "photon_data_over_fit_sdcc_vs_current_overlay_current_root_points.csv"
)
CURRENT_POINTS_MANIFEST = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_fig24_photonjet_fix_20260702_014217/"
    / "strict_stitched_photonjet_current/"
    / "photon_data_over_fit_sdcc_vs_current_overlay_current_root_points_manifest.json"
)
PREVIOUS_AUDIT = (
    REPO
    / "dataOutput/ppg12Parity/control_plane/audits/"
    / "the76_photonjet_fig5_identity_20260702"
)

PPG12_BUILD = REPO / "ppg12codeGit/plotting/build_h_max_photon_pT_uncut.py"
PPG12_PLOT = REPO / "ppg12codeGit/plotting/plot_combine_uncut.C"
PPG12_CALOANA = REPO / "ppg12codeGit/anatreemaker/source/CaloAna24.cc"
PPG12_CALOANA_H = REPO / "ppg12codeGit/anatreemaker/source/CaloAna24.h"
PPG12_RUN28 = REPO / "ppg12codeGit/anatreemaker/macro_maketree/sim/run28"
RECOILJETS = REPO / "src/RecoilJets.cc"

CURRENT_HIST = "SIM/h_ppPhotonStitch_ppg12Fig5_maxPhotonPt_kept"


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_json(path: Path) -> Any:
    return json.loads(path.read_text())


def snippet(path: Path, start: int, end: int) -> str:
    lines = path.read_text(errors="replace").splitlines()
    out = []
    for line_no in range(start, min(end, len(lines)) + 1):
        out.append(f"{path.relative_to(REPO)}:{line_no}: {lines[line_no - 1]}")
    return "\n".join(out)


def hist_stats(path: Path, names: list[str]) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    with uproot.open(path) as f:
        for name in names:
            if name not in f:
                out[name] = {"exists": False}
                continue
            h = f[name]
            vals, _edges = h.to_numpy(flow=False)
            stat: dict[str, Any] = {
                "exists": True,
                "classname": h.classname,
                "integral_visible_bins": float(np.sum(vals)),
                "nonzero_bins": int(np.count_nonzero(vals)),
                "max_bin_content": float(np.max(vals)) if len(vals) else 0.0,
            }
            for member in ("fEntries", "fTsumw", "fTsumw2"):
                try:
                    stat[member] = float(h.member(member))
                except Exception:
                    stat[member] = None
            if stat.get("fEntries") not in (None, 0.0):
                stat["fTsumw_over_fEntries"] = stat["fTsumw"] / stat["fEntries"]
            try:
                stat["weighted"] = bool(h.weighted)
            except Exception:
                stat["weighted"] = None
            out[name] = stat
    return out


def row(
    axis: str,
    ppg12: str,
    recoil: str,
    status: str,
    could_explain: str,
    evidence: str,
    next_proof: str = "",
) -> dict[str, str]:
    return {
        "transformation_axis": axis,
        "ppg12_anatreemaker_slimtree_behavior": ppg12,
        "recoiljets_direct_dst_behavior": recoil,
        "same_different_unknown": status,
        "could_explain_per_mille_to_percent_residual": could_explain,
        "exact_evidence_path_or_command": evidence,
        "next_proof_needed_if_unknown": next_proof,
    }


def markdown_table(rows: list[dict[str, str]]) -> str:
    fields = [
        "transformation_axis",
        "ppg12_anatreemaker_slimtree_behavior",
        "recoiljets_direct_dst_behavior",
        "same_different_unknown",
        "could_explain_per_mille_to_percent_residual",
        "exact_evidence_path_or_command",
        "next_proof_needed_if_unknown",
    ]
    text = ["| " + " | ".join(fields) + " |", "| " + " | ".join(["---"] * len(fields)) + " |"]
    for r in rows:
        text.append("| " + " | ".join(str(r.get(f, "")).replace("\n", "<br>") for f in fields) + " |")
    return "\n".join(text) + "\n"


def load_previous() -> dict[str, Any]:
    out: dict[str, Any] = {}
    for name in (
        "final_classification.json",
        "raw_root_bin_identity.json",
        "component_prestitch_identity.json",
    ):
        path = PREVIOUS_AUDIT / name
        out[name] = read_json(path) if path.exists() else {"missing": str(path)}
    if CURRENT_POINTS_MANIFEST.exists():
        out["current_points_manifest.json"] = read_json(CURRENT_POINTS_MANIFEST)
    return out


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    previous = load_previous()
    root_stats = hist_stats(
        CURRENT_ROOT,
        [
            CURRENT_HIST,
            "SIM/h_ppPhotonStitch_ppg12Fig5_maxPhotonPt_all",
            "SIM/h_ppPhotonStitch_ppg12Fig5_maxPhotonPt_rejected",
            "SIM/h_ppg12_vtxqa_sim_truth_z_unweighted_pre_vzcut",
            "SIM/h_ppg12_vtxqa_sim_truth_z_unweighted_post_vzcut",
            "SIM/h_ppg12_vtxqa_sim_period_event_weight_pre_vzcut",
            "SIM/h_ppg12_vtxqa_sim_period_event_weight_post_vzcut",
            "SIM/h_ppg12_vtxqa_sim_vertex_weight_pre_vzcut",
            "SIM/h_ppg12_vtxqa_sim_vertex_weight_post_vzcut",
        ],
    )

    ppg12_summary = {
        "fig5_builder": str(PPG12_BUILD),
        "source_slimtree_template": "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/{sample}/condorout/combined.root",
        "output_root": "/sphenix/user/shuhangli/ppg12/plotting/photon_max_pT_uncut.root",
        "objects": [
            "h_max_photon_pT_photon5",
            "h_max_photon_pT_photon10",
            "h_max_photon_pT_photon20",
            "h_max_photon_pT_photon5_sumw2",
            "h_max_photon_pT_photon10_sumw2",
            "h_max_photon_pT_photon20_sumw2",
        ],
        "builder_evidence": [
            snippet(PPG12_BUILD, 20, 24),
            snippet(PPG12_BUILD, 27, 37),
            snippet(PPG12_BUILD, 39, 53),
            snippet(PPG12_BUILD, 65, 75),
        ],
        "slimtree_event_save_evidence": [
            snippet(PPG12_CALOANA, 760, 770),
            snippet(PPG12_CALOANA, 882, 888),
            snippet(PPG12_CALOANA, 975, 1002),
            snippet(PPG12_CALOANA, 2140, 2153),
        ],
        "run28_production_evidence": [
            snippet(PPG12_RUN28 / "photon5/run_condor.sh", 5, 13),
            snippet(PPG12_RUN28 / "hadd_combined.sh", 13, 18),
            snippet(PPG12_RUN28 / "hadd_combined.sh", 46, 56),
        ],
    }

    recoil_summary = {
        "current_root": str(CURRENT_ROOT),
        "current_histogram": CURRENT_HIST,
        "current_points_csv": str(CURRENT_POINTS_CSV),
        "artifact_pointer": str(REPO / "dataOutput/current_recoiljets_artifacts/current/pp_sim_photonjet_merged/current.json"),
        "source_evidence": [
            snippet(RECOILJETS, 111, 130),
            snippet(RECOILJETS, 714, 762),
            snippet(RECOILJETS, 1016, 1055),
            snippet(RECOILJETS, 3669, 3738),
            snippet(RECOILJETS, 4688, 4715),
            snippet(RECOILJETS, 4999, 5024),
            snippet(RECOILJETS, 5078, 5112),
        ],
        "root_histogram_stats": root_stats,
    }

    pre = root_stats["SIM/h_ppg12_vtxqa_sim_truth_z_unweighted_pre_vzcut"].get("fEntries")
    post = root_stats["SIM/h_ppg12_vtxqa_sim_truth_z_unweighted_post_vzcut"].get("fEntries")
    vz_fraction = (post / pre) if pre else math.nan

    rows = [
        row(
            "DST file list",
            "run_condor.sh slices sample lists from /sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/{sample}.",
            "Current final ROOT does not contain the source DST lists; registered artifact records only the merged ROOT.",
            "unknown",
            "yes",
            "ppg12codeGit/anatreemaker/macro_maketree/sim/run28/photon*/run_condor.sh; dataOutput/current_recoiljets_artifacts/current/pp_sim_photonjet_merged/current.json",
            "Read-only SDCC file-list digest for photon5/10/20 plus the July2 RecoilJets input manifest.",
        ),
        row(
            "event/entry preservation",
            "slimtree->Fill() only when saveevent is true; saveevent is set by ncluster>0 or nparticles>0. Fig.5 builder then uses slimtree.num_entries as n_total.",
            "RecoilJets Fig.5 fills after firstEventCuts(); exact ROOT QA shows truth-z entries change from pre_vzcut to post_vzcut.",
            "different",
            "yes",
            f"{PPG12_CALOANA}:2140-2153; {PPG12_BUILD}:30-37; {RECOILJETS}:3669-3738; current ROOT pre/post fEntries={pre}/{post} (post/pre={vz_fraction:.6f})",
            "",
        ),
        row(
            "duplicate handling",
            "PPG12 hadd merges OutDir*/caloana.root into combined.root with hadd -f -j 1 -n 50.",
            "RecoilJets current final ROOT is a merged artifact; per-job duplicate accounting is not stored in the final ROOT.",
            "unknown",
            "yes",
            "ppg12codeGit/anatreemaker/macro_maketree/sim/run28/hadd_combined.sh:46-56; current artifact manifest only records final merged ROOT.",
            "Compare per-job lists and merge manifests on SDCC read-only.",
        ),
        row(
            "sample labels photon5/photon10/photon20",
            "Fig.5 builder loops photon5, photon10, photon20.",
            "RecoilJets detects photon sample from RJ_PPG12_PHOTON_SAMPLE/RJ_PHOTONJET_SAMPLE/RJ_SIM_SAMPLE or outfile text.",
            "same",
            "no",
            f"{PPG12_BUILD}:68-70; {RECOILJETS}:714-727",
            "",
        ),
        row(
            "generated-event/source-count denominator",
            "w_evt = XSEC[sample] / slimtree.num_entries, not PHGenIntegral and not an external generated-event count.",
            "Plot manifest uses a constant validated source-scope factor on the current ROOT; exact denominator lineage is not encoded in the final ROOT.",
            "different",
            "yes",
            f"{PPG12_BUILD}:30-37; {CURRENT_POINTS_MANIFEST}",
            "Read-only SDCC check of slimtree entries, h_sim_cross_counting, PHGenIntegral, and current RecoilJets source denominators.",
        ),
        row(
            "truth photon collection",
            "CaloAna24 seeds primary_particles from GetPrimaryParticleRange() and requires trutheval->get_embed(truth) >= 1.",
            "RecoilJets ppg12MaxStoredTruthPhotonPt uses G4TruthInfo primary range and requireEmbed=true.",
            "same",
            "not by itself",
            f"{PPG12_CALOANA}:760-770; {RECOILJETS}:1016-1036",
            "Same-DST canary still needed to prove row-level equality.",
        ),
        row(
            "photon PID/status",
            "Fig.5 builder uses branch particle_pid == 22. CaloAna24 stores pid for embedded primary particles after pT/E/eta cuts.",
            "RecoilJets requires p->get_pid() == 22 on G4Truth primary particles.",
            "same",
            "not by itself",
            f"{PPG12_BUILD}:39-48; {PPG12_CALOANA}:988-993; {RECOILJETS}:1035-1049",
            "Same-DST canary needed for exact value comparison.",
        ),
        row(
            "photon promptness/source category",
            "Fig.5 builder does not use particle_photonclass or mother PID.",
            "RecoilJets Fig.5 does not require promptness/source category.",
            "same",
            "no",
            f"{PPG12_BUILD}:39-48; {RECOILJETS}:1035-1055",
            "",
        ),
        row(
            "max photon pT definition",
            "max over slimTree particle_Pt for particle_pid==22 and |particle_Eta|<0.7.",
            "max over direct G4 primary embedded photons with pid 22, pT>1, E>1, finite eta, |eta|<0.7.",
            "unknown",
            "yes",
            f"{PPG12_BUILD}:41-49; {RECOILJETS}:1035-1055",
            "Same-DST event canary comparing PPG12 branch value to direct RecoilJets computed max pT.",
        ),
        row(
            "eta/acceptance cuts",
            "slimTree stores particles with |eta|<=1.5 and Fig.5 builder applies |eta|<0.7.",
            "RecoilJets Fig.5 call applies etaAbsMax=0.7 after G4 primary scan.",
            "same",
            "not by itself",
            f"{PPG12_CALOANA}:886-888; {PPG12_BUILD}:42-48; {RECOILJETS}:1049-1050",
            "",
        ),
        row(
            "z-vertex cuts",
            "No Fig.5 z-vertex cut is applied in build_h_max_photon_pT_uncut.py; no explicit |z| cut appears before slimtree->Fill() in the inspected CaloAna24 path.",
            "RecoilJets Fig.5 is after firstEventCuts(); for sim this applies global |vz| if enabled. Exact current ROOT has pre_vzcut fEntries 119,994,000 and post_vzcut fEntries 77,280,278.",
            "different",
            "yes",
            f"{PPG12_BUILD}:39-53; {PPG12_CALOANA}:2140-2153; {RECOILJETS}:3700-3721; current ROOT h_ppg12_vtxqa_sim_truth_z_unweighted_pre/post_vzcut",
            "A no-vz-cut RecoilJets canary or same-DST event table would quantify the exact bin shifts.",
        ),
        row(
            "event-quality cuts",
            "Anatreemaker requires required nodes and saves only events with clusters or particles; no Fig.5-specific trigger/vz gate in the builder.",
            "RecoilJets requires fetchNodes() and firstEventCuts() before Fig.5 fill.",
            "different",
            "yes",
            f"{PPG12_CALOANA}:2140-2153; {RECOILJETS}:4450-4455 and 4605-4636",
            "Same-DST canary should print row inclusion/exclusion reason.",
        ),
        row(
            "period split",
            "PPG12 Fig.5 source path in build script points to run28/{sample} SI-like sample names only.",
            "Current RecoilJets ROOT is the July2 photonjet merged artifact; final ROOT does not expose per-period component identity in the Fig.5 object.",
            "unknown",
            "yes",
            f"{PPG12_BUILD}:20-22; {CURRENT_POINTS_MANIFEST}",
            "Read-only current merge manifest/per-run inventory.",
        ),
        row(
            "SI/DI handling",
            "PPG12 Fig.5 builder loops photon5/photon10/photon20, not *_double sample names.",
            "RecoilJets current artifact may include period/SI-DI production context; exact component composition is not encoded in the single final Fig.5 histogram.",
            "unknown",
            "yes",
            f"{PPG12_BUILD}:68-70; {CURRENT_ROOT}",
            "Read current production manifest and per-component final ROOTs, if present.",
        ),
        row(
            "cross-section weight",
            "XSEC constants are photon5=146359.3, photon10=6944.675, photon20=130.4461 and enter as XSEC/n_total.",
            "RecoilJets has matching xsec constants, but current plot uses post-hoc source-scope factor on the final merged ROOT.",
            "different",
            "yes",
            f"{PPG12_BUILD}:20 and 36; {RECOILJETS}:770-777; {CURRENT_POINTS_MANIFEST}",
            "Recover exact RecoilJets per-sample normalization denominator.",
        ),
        row(
            "vertex reweight",
            "PPG12 Fig.5 builder does not read or apply vertex-reweight branches/files.",
            "Current RecoilJets ROOT contains weighted vertex QA and Fig.5 histogram is weighted (fTsumw != fEntries); local source line order alone is not exact-artifact proof of when the factor enters Fig.5.",
            "different",
            "yes",
            f"{PPG12_BUILD}:39-53; current ROOT {CURRENT_HIST} fEntries={root_stats[CURRENT_HIST]['fEntries']} fTsumw={root_stats[CURRENT_HIST]['fTsumw']}; {RECOILJETS}:5078-5112",
            "Foreground canary should print final fill weight at the Fig.5 fill call.",
        ),
        row(
            "luminosity/period weight",
            "Not used by build_h_max_photon_pT_uncut.py.",
            "RecoilJets has period/mix/lumi machinery for PPG12 photon-yield context; exact Fig.5 fill-stage proof is missing, but current ROOT weighted QA is present.",
            "different",
            "yes",
            f"{PPG12_BUILD}:39-53; {RECOILJETS}:5078-5112; current ROOT h_ppg12_vtxqa_sim_period_event_weight_pre/post_vzcut",
            "Same-DST canary/fill-stage printout.",
        ),
        row(
            "SI/DI mix weight",
            "Not used by build_h_max_photon_pT_uncut.py.",
            "RecoilJets multiplies m_ppg12PhotonYieldMixWeight in PPG12 MC context; exact Fig.5 fill-stage linkage remains unproven.",
            "different",
            "yes",
            f"{RECOILJETS}:5104-5112",
            "Same-DST canary/fill-stage printout.",
        ),
        row(
            "final event fill weight",
            "One constant per sample, XSEC/slimtree_entries.",
            "Current ROOT Fig.5 histogram is weighted; current plotting multiplies by a constant source-scope factor. This is not the same as the PPG12 builder's direct XSEC/slimtree_entries construction.",
            "different",
            "yes",
            f"{PPG12_BUILD}:36 and 51-53; current ROOT {CURRENT_HIST} weighted={root_stats[CURRENT_HIST]['weighted']} fTsumw_over_fEntries={root_stats[CURRENT_HIST]['fTsumw_over_fEntries']}",
            "Instrument RecoilJets Fig.5 fill to print factorized weight, or inspect exact run logs if already emitted.",
        ),
        row(
            "binning",
            "0.5 GeV edges from 0 to 50.",
            "Current plotted extraction uses matched 0.5 GeV centers from 10 to 40; current ROOT histogram is 0 to 60 with 0.5 GeV bins.",
            "same for plotted range",
            "no",
            f"{PPG12_BUILD}:23; {RECOILJETS}:5011-5013; {CURRENT_POINTS_CSV}",
            "",
        ),
        row(
            "underflow/overflow",
            "Builder uses numpy.histogram over explicit 0-50 edges; values outside are not included in H.",
            "RecoilJets ROOT histogram has 0-60 range; plotted range uses 10-40 bins only.",
            "different outside plotted range",
            "unlikely for 10-40 central residual",
            f"{PPG12_BUILD}:23 and 51-53; {RECOILJETS}:5011-5013",
            "Only relevant if normalization includes overflow; current plot uses visible bins, so not primary.",
        ),
        row(
            "stitch/sample ownership",
            "PPG12 plot_combine_uncut uses photon5 c<14, photon10 14<=c<22, photon20 c>=22.",
            "RecoilJets ppg12PhotonSliceWindow uses [0,14), [14,22), [22,inf).",
            "same",
            "no",
            f"{PPG12_PLOT}:95-118; {RECOILJETS}:751-762",
            "",
        ),
    ]

    write_csv(
        OUT_DIR / "slimtree_transformation_table.csv",
        rows,
        [
            "transformation_axis",
            "ppg12_anatreemaker_slimtree_behavior",
            "recoiljets_direct_dst_behavior",
            "same_different_unknown",
            "could_explain_per_mille_to_percent_residual",
            "exact_evidence_path_or_command",
            "next_proof_needed_if_unknown",
        ],
    )
    write_json(OUT_DIR / "slimtree_transformation_table.json", rows)
    write_text(
        OUT_DIR / "slimtree_transformation_table.md",
        "# THE-76 Photon+Jet Fig.5 SlimTree Transformation Table\n\n"
        + markdown_table(rows),
    )

    same_dst_canary = {
        "status": "blocked_missing_local_ppg12_slimtree_and_same_event_manifest",
        "central_question": "For the same DST event, do PPG12 and RecoilJets fill the same Fig.5 bin with the same weight?",
        "missing_local_artifacts": [
            "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/{photon5,photon10,photon20}/condorout/combined.root",
            "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/{sample}/{dst,g4hits,dst_calo_cluster,dst_mbd_epd,dst_truth_jet}.list",
            "Exact July2 RecoilJets per-component input manifest or per-event diagnostic table",
        ],
        "smallest_read_only_sdcc_check": [
            "For photon5/photon10/photon20, open PPG12 combined.root with uproot and print slimtree.num_entries, h_sim_cross_counting if present, branch names, and the first 100 rows of event identifiers plus leading particle_pid==22 |eta|<0.7 pT.",
            "Open the corresponding RecoilJets per-component ROOTs or rerun a <=100-event foreground canary in scratch, with no Condor, printing event id, vz pass, sample window, max photon pT, bin, and factorized fill weight.",
            "Compare same source-file/event identifiers before any broad rerun.",
        ],
        "required_columns": [
            "sample",
            "source_file",
            "entry/run/event identifier",
            "appears_in_ppg12_slimtree",
            "appears_in_recoiljets_direct_read",
            "ppg12_photon_pt_branch_value",
            "recoiljets_computed_max_photon_pt",
            "ppg12_selected",
            "recoiljets_selected",
            "ppg12_bin",
            "recoiljets_bin",
            "ppg12_weight",
            "recoiljets_cross_section_or_source_weight",
            "period_lumi_weight",
            "si_di_mix_weight",
            "vertex_reweight",
            "final_fill_weight",
            "mismatch_reason",
        ],
    }

    component_localization = {
        "status": "residual visible in all three owned sample regions; true pre-stitch component identity is not locally testable because current final ROOT lacks separate per-sample component histograms",
        "previous_component_summary": previous.get("component_prestitch_identity.json", {}),
        "localization_candidates": [
            "event population and |vz| gate mismatch",
            "source-count/slimtree-entry denominator mismatch",
            "weighted RecoilJets Fig.5 object versus PPG12 XSEC/slimtree_entries builder",
            "same-DST max-photon-pT value mismatch not yet tested",
        ],
    }

    source_count_status = {
        "ppg12_builder_denominator": "slimtree.num_entries from /sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/{sample}/condorout/combined.root",
        "recoiljets_denominator": "not encoded as a matching per-sample source denominator in the current final ROOT; plot uses source_scope_factor from previous strict extraction",
        "current_manifest": previous.get("current_points_manifest.json", {}),
        "status": "not proven identical",
        "consequence": "Exact central-value identity cannot be claimed until the PPG12 slimtree entry counts and the RecoilJets source normalization denominators are locked and shown identical or intentionally transformed.",
    }

    vertex_status = {
        "classification": "ambiguous",
        "reason": "PPG12 Fig.5 builder has no vertex-reweight input, while the current RecoilJets final ROOT is weighted and contains vertex QA. However the exact Fig.5 fill-stage factor cannot be proven from local source alone because the ROOT object is weighted and the local source line order is not runtime proof for the already-produced artifact.",
        "positive_evidence": {
            "current_root_weighted_fig5": root_stats[CURRENT_HIST],
            "vtx_pre_post": {
                "truth_z_unweighted_pre_vzcut": root_stats["SIM/h_ppg12_vtxqa_sim_truth_z_unweighted_pre_vzcut"],
                "truth_z_unweighted_post_vzcut": root_stats["SIM/h_ppg12_vtxqa_sim_truth_z_unweighted_post_vzcut"],
                "period_event_weight_pre_vzcut": root_stats["SIM/h_ppg12_vtxqa_sim_period_event_weight_pre_vzcut"],
                "period_event_weight_post_vzcut": root_stats["SIM/h_ppg12_vtxqa_sim_period_event_weight_post_vzcut"],
            },
        },
        "next_proof": "A <=100-event foreground canary must print CurrentWeight/final fill weight at the Fig.5 Fill call for photon5/10/20.",
    }

    final = {
        "classification": "candidate cause: PPG12 slimTree event filtering differs",
        "canonicalization_allowed_now": False,
        "first_causal_identity_break_candidate": "PPG12 Fig.5 builder reads slimTree rows without a Fig.5 |vz| gate and normalizes by slimtree.num_entries; RecoilJets current Fig.5 object is filled after firstEventCuts() and the exact ROOT shows a large pre/post vz-cut event population change.",
        "quantitative_anchor": {
            "current_root_truth_z_unweighted_pre_vzcut_fEntries": pre,
            "current_root_truth_z_unweighted_post_vzcut_fEntries": post,
            "post_over_pre": vz_fraction,
            "previous_max_relative_fig5_residual": previous.get("raw_root_bin_identity.json", {}).get("max_extracted_relative_difference"),
        },
        "candidate_causal_mechanisms": [
            "Event population mismatch from RecoilJets SIM |vz|/firstEventCuts gate versus PPG12 Fig.5 slimTree builder lacking a corresponding Fig.5 cut.",
            "Denominator mismatch because PPG12 uses slimtree.num_entries while current plotting uses a constant source-scope factor on a weighted RecoilJets final ROOT.",
            "Weight-handling mismatch because PPG12 Fig.5 applies a constant XSEC/slimtree_entries sample weight, while the current RecoilJets Fig.5 ROOT object is weighted (fTsumw != fEntries).",
            "Same-DST max-photon-pT equality remains unproven without the PPG12 slimTree rows and event identifiers.",
        ],
        "ruled_out": [
            "plotting/CSV extraction as the source of non-unity central values",
            "floating-point, ROOT serialization, or CSV precision as an explanation for per-mille/percent residuals",
            "sample ownership window mismatch for the plotted [0,14), [14,22), [22,inf) windows",
            "a pure final-plot style problem",
        ],
        "ambiguous": [
            "exact PPG12 slimTree source file list and event identifiers",
            "exact RecoilJets per-component input population for the current merged artifact",
            "row-level equality of PPG12 stored particle_Pt and RecoilJets direct max photon pT",
            "which exact event-weight factor entered the current Fig.5 histogram fill in the already-produced ROOT",
        ],
        "minimal_next_diagnostic": "Run a tiny read-only SDCC foreground same-DST canary for photon5/10/20 that prints slimTree row presence, direct RecoilJets event presence, max photon pT, vz gate, bin, PPG12 XSEC/slimtree_entries weight, RecoilJets factorized fill weight, and mismatch reason for <=100 shared events.",
    }

    write_json(OUT_DIR / "ppg12_slimtree_chain_summary.json", ppg12_summary)
    write_text(
        OUT_DIR / "ppg12_slimtree_chain_summary.md",
        "# PPG12 anatreemaker/slimTree Fig.5 source chain\n\n"
        "The Fig.5 source path is `build_h_max_photon_pT_uncut.py`, which reads "
        "`run28/{sample}/condorout/combined.root:slimtree`, uses "
        "`particle_pid`, `particle_Pt`, and `particle_Eta`, and weights each "
        "sample with `XSEC / slimtree.num_entries`.\n\n"
        "## Evidence\n\n"
        + "\n\n".join(f"```text\n{x}\n```" for x in ppg12_summary["builder_evidence"] + ppg12_summary["slimtree_event_save_evidence"])
    )
    write_json(OUT_DIR / "recoiljets_direct_chain_summary.json", recoil_summary)
    write_text(
        OUT_DIR / "recoiljets_direct_chain_summary.md",
        "# RecoilJets direct-DST Fig.5 source chain\n\n"
        "The current artifact is the registered `pp_sim_photonjet_merged` ROOT. "
        "The Fig.5 object is `SIM/h_ppPhotonStitch_ppg12Fig5_maxPhotonPt_kept`. "
        "The exact ROOT object is weighted and the vertex QA shows a pre/post "
        "event-population change before the final object.\n\n"
        "## ROOT stats\n\n"
        f"```json\n{json.dumps(root_stats, indent=2, sort_keys=True)}\n```\n\n"
        "## Code evidence\n\n"
        + "\n\n".join(f"```text\n{x}\n```" for x in recoil_summary["source_evidence"])
    )
    write_json(OUT_DIR / "same_dst_event_canary_plan.json", same_dst_canary)
    write_text(
        OUT_DIR / "same_dst_event_canary_plan.md",
        "# Same-DST event canary plan\n\n"
        f"Status: `{same_dst_canary['status']}`\n\n"
        "Central question: "
        + same_dst_canary["central_question"]
        + "\n\n## Missing local artifacts\n\n"
        + "\n".join(f"- `{x}`" for x in same_dst_canary["missing_local_artifacts"])
        + "\n\n## Smallest read-only SDCC check\n\n"
        + "\n".join(f"- {x}" for x in same_dst_canary["smallest_read_only_sdcc_check"])
        + "\n\n## Required columns\n\n"
        + "\n".join(f"- `{x}`" for x in same_dst_canary["required_columns"])
        + "\n",
    )
    write_json(OUT_DIR / "component_level_localization.json", component_localization)
    write_text(
        OUT_DIR / "component_level_localization.md",
        "# Component-level localization\n\n"
        f"{component_localization['status']}\n\n"
        "Candidate locations:\n"
        + "\n".join(f"- {x}" for x in component_localization["localization_candidates"])
        + "\n",
    )
    write_json(OUT_DIR / "source_count_denominator_status.json", source_count_status)
    write_text(
        OUT_DIR / "source_count_denominator_status.md",
        "# Source-count / generated-denominator status\n\n"
        f"- PPG12: {source_count_status['ppg12_builder_denominator']}\n"
        f"- RecoilJets/current: {source_count_status['recoiljets_denominator']}\n"
        f"- Status: `{source_count_status['status']}`\n\n"
        f"Consequence: {source_count_status['consequence']}\n",
    )
    write_json(OUT_DIR / "vertex_reweight_slimtree_context.json", vertex_status)
    write_text(
        OUT_DIR / "vertex_reweight_slimtree_context.md",
        "# Vertex-reweight status in the slimTree context\n\n"
        f"Classification: `{vertex_status['classification']}`\n\n"
        f"{vertex_status['reason']}\n\n"
        f"Next proof: {vertex_status['next_proof']}\n",
    )
    write_json(OUT_DIR / "final_classification.json", final)
    write_text(
        OUT_DIR / "final_classification.md",
        "# THE-76 Fig.5 slimTree-layer classification\n\n"
        f"Classification: `{final['classification']}`\n\n"
        f"First causal identity-break candidate: {final['first_causal_identity_break_candidate']}\n\n"
        "## Quantitative anchor\n\n"
        f"```json\n{json.dumps(final['quantitative_anchor'], indent=2, sort_keys=True)}\n```\n\n"
        "## Candidate causal mechanisms\n\n"
        + "\n".join(f"- {x}" for x in final["candidate_causal_mechanisms"])
        + "\n\n## Ruled out\n\n"
        + "\n".join(f"- {x}" for x in final["ruled_out"])
        + "\n\n## Ambiguous\n\n"
        + "\n".join(f"- {x}" for x in final["ambiguous"])
        + "\n\n"
        f"Canonicalization allowed now: `{str(final['canonicalization_allowed_now']).lower()}`\n\n"
        f"Minimal next diagnostic: {final['minimal_next_diagnostic']}\n",
    )
    write_json(
        OUT_DIR / "audit_index.json",
        {
            "out_dir": str(OUT_DIR),
            "current_root": str(CURRENT_ROOT),
            "current_histogram": CURRENT_HIST,
            "previous_identity_audit": str(PREVIOUS_AUDIT),
            "artifacts": sorted(str(p) for p in OUT_DIR.iterdir()),
            "classification": final["classification"],
            "canonicalization_allowed_now": False,
        },
    )

    print(f"Wrote slimTree equivalence audit to {OUT_DIR}")
    print(f"Classification: {final['classification']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
