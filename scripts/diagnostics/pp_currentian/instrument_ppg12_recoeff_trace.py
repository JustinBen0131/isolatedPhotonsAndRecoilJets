#!/usr/bin/env python3
"""Deterministically instrument the preserved PPG12 RecoEff macro.

The transformation is deliberately narrow and fail closed.  It adds a CSV
side channel for the candidate decisions already computed by the executable;
it does not replace any selection, weight, histogram fill, or RNG expression.
The paired harness must run both this copy and the uninstrumented staged copy
and prove exact ROOT-object/content/Sumw2 equivalence before using the trace as
preserved-executable evidence.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


SCHEMA_VERSION = 1
TRANSFORM_NAME = "ppg12_recoeff_candidate_trace_v1"


class TransformFailure(RuntimeError):
    pass


def digest_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def replace_exact(text: str, old: str, new: str, label: str) -> str:
    count = text.count(old)
    if count != 1:
        raise TransformFailure(
            f"{label}: expected exactly one insertion marker, observed {count}"
        )
    return text.replace(old, new)


def instrument(text: str) -> tuple[str, list[dict[str, object]]]:
    operations: list[dict[str, object]] = []

    def apply(old: str, new: str, label: str) -> None:
        nonlocal text
        text = replace_exact(text, old, new, label)
        operations.append(
            {
                "label": label,
                "expected_count": 1,
                "marker_sha256": digest_bytes(old.encode()),
                "replacement_sha256": digest_bytes(new.encode()),
            }
        )

    apply(
        "    TTreeReaderValue<int> runnumber(reader, \"runnumber\");\n",
        "    TTreeReaderValue<int> runnumber(reader, \"runnumber\");\n"
        "    TTreeReaderValue<int> oracle_eventnumber(reader, \"eventnumber\");\n",
        "event_identity_reader",
    )
    apply(
        "    int ientry = 0;\n\n    while (reader.Next())\n",
        r'''    int ientry = 0;

    const char *oracle_trace_path = gSystem->Getenv("RJ_PPG12_EXEC_TRACE_CSV");
    const char *oracle_response_trace_path =
        gSystem->Getenv("RJ_PPG12_EXEC_RESPONSE_TRACE_CSV");
    std::ofstream oracle_trace;
    std::ofstream oracle_response_trace;
    if (oracle_trace_path && oracle_trace_path[0] != '\0')
    {
        oracle_trace.open(oracle_trace_path, std::ios::out | std::ios::trunc);
        if (!oracle_trace.is_open())
        {
            std::cerr << "[PPG12OracleTrace] ERROR: cannot open "
                      << oracle_trace_path << std::endl;
            return;
        }
        oracle_trace
            << "tree_entry,chain_file_index,local_tree_entry,runnumber,eventnumber,"
            << "cluster_index,truth_track_id,"
            << "cluster_Et,raw_eiso,"
            << "corrected_eiso,iso_threshold,noniso_threshold,sample_weight,mix_weight,"
            << "lumi_weight,cross_weight,vertex_weight,truth_vertex_weight,"
            << "trigger_weight,event_weight,weight,selected_model,"
            << "base_E_score,base_v3E_score,selected_score,cluster_weta_cogx,"
            << "cluster_wphi_cogx,vertexz,cluster_Eta,e11_over_e33,cluster_et1,"
            << "cluster_et2,cluster_et3,cluster_et4,e32_over_e35,common_pass,"
            << "tight,nontight,is_iso,is_noniso,logical_abcd_region,is_signal,"
            << "truth_particle_index,truth_class,truth_pt,"
            << "analysis_window_pass,signal_fill_A,signal_fill_B,signal_fill_C,"
            << "signal_fill_D,fill_multiplicity\n";
        oracle_trace.precision(17);
    }
    if (oracle_response_trace_path && oracle_response_trace_path[0] != '\0')
    {
        oracle_response_trace.open(
            oracle_response_trace_path, std::ios::out | std::ios::trunc);
        if (!oracle_response_trace.is_open())
        {
            std::cerr << "[PPG12OracleTrace] ERROR: cannot open response trace "
                      << oracle_response_trace_path << std::endl;
            return;
        }
        oracle_response_trace
            << "tree_entry,chain_file_index,local_tree_entry,runnumber,eventnumber,"
            << "cluster_index,response_Et,"
            << "response_window_pass\n";
        oracle_response_trace.precision(17);
    }

    while (reader.Next())
''',
        "trace_stream_setup",
    )
    apply(
        "            // fudge the MC isoET\n            if (issim)\n",
        "            const float oracle_raw_recoisoET = recoisoET;\n"
        "            // fudge the MC isoET\n            if (issim)\n",
        "capture_raw_isolation",
    )
    apply(
        "        weight = cross_weight;\n        vertex_weight = 1.0;\n",
        "        weight = cross_weight;\n"
        "        vertex_weight = 1.0;\n"
        "        float oracle_truth_vertex_weight = 1.0;\n",
        "truth_vertex_factor_scope",
    )
    apply(
        "                weight *= vertexz_truth_mb_ptr\n"
        "                    ? TruthVertexWeight(h_truth_vtx_reweight, *vertexz_truth, **vertexz_truth_mb_ptr)\n"
        "                    : TruthVertexWeight(h_truth_vtx_reweight, *vertexz_truth);\n",
        "                weight *= vertexz_truth_mb_ptr\n"
        "                    ? TruthVertexWeight(h_truth_vtx_reweight, *vertexz_truth, **vertexz_truth_mb_ptr)\n"
        "                    : TruthVertexWeight(h_truth_vtx_reweight, *vertexz_truth);\n"
        "                oracle_truth_vertex_weight = vertexz_truth_mb_ptr\n"
        "                    ? TruthVertexWeight(h_truth_vtx_reweight, *vertexz_truth, **vertexz_truth_mb_ptr)\n"
        "                    : TruthVertexWeight(h_truth_vtx_reweight, *vertexz_truth);\n",
        "capture_truth_vertex_factor",
    )
    apply(
        "            bool common_pass = false;\n"
        "            bool tight = false;\n"
        "            bool nontight = false;\n"
        "            bool iso = false;\n"
        "            bool noniso = false;\n",
        "            bool common_pass = false;\n"
        "            bool tight = false;\n"
        "            bool nontight = false;\n"
        "            bool iso = false;\n"
        "            bool noniso = false;\n"
        "            std::string selected_bdt_model = bdt_model_name;\n"
        "            float bdt_score = std::numeric_limits<float>::quiet_NaN();\n",
        "trace_model_scope",
    )
    apply(
        "                std::string selected_bdt_model = bdt_model_name;\n",
        "                selected_bdt_model = bdt_model_name;\n",
        "reuse_scoped_model",
    )
    apply(
        "                float bdt_score = (*bdt_arrays[selected_bdt_model])[icluster];\n",
        "                bdt_score = (*bdt_arrays[selected_bdt_model])[icluster];\n",
        "reuse_scoped_score",
    )
    trace_block = r'''            if (oracle_trace.is_open())
            {
                int oracle_truth_index = -1;
                int oracle_truth_class = -1;
                float oracle_truth_pt = std::numeric_limits<float>::quiet_NaN();
                bool oracle_is_signal = false;
                auto oracle_track_it = particle_trkidmap.find(cluster_truthtrkID[icluster]);
                if (oracle_track_it != particle_trkidmap.end())
                {
                    oracle_truth_index = oracle_track_it->second;
                    oracle_truth_class = particle_photonclass[oracle_truth_index];
                    oracle_truth_pt = particle_Pt[oracle_truth_index];
                    oracle_is_signal = photon_reco.find(oracle_truth_index) != photon_reco.end();
                }
                const int oracle_region = tight && iso ? 1 :
                                          tight && noniso ? 2 :
                                          nontight && iso ? 3 :
                                          nontight && noniso ? 4 : 0;
                const int oracle_analysis_window = oracle_is_signal &&
                    oracle_truth_pt > pTmin_truth && oracle_truth_pt < pTmax_truth &&
                    cluster_Et[icluster] > pTmin && cluster_Et[icluster] < pTmax;
                const int oracle_fill_a =
                    oracle_is_signal && oracle_region == 1 && oracle_analysis_window;
                const int oracle_fill_b = oracle_is_signal && oracle_region == 2;
                const int oracle_fill_c = oracle_is_signal && oracle_region == 3;
                const int oracle_fill_d = oracle_is_signal && oracle_region == 4;
                const int oracle_fill_multiplicity =
                    oracle_fill_a + oracle_fill_b + oracle_fill_c + oracle_fill_d;
                std::string oracle_selected_model = bdt_model_name;
                if (use_et_binned_bdt)
                {
                    const float oracle_et = cluster_Et[icluster];
                    for (int ib = 0; ib < (int)bdt_et_bin_models.size(); ++ib)
                    {
                        if (oracle_et >= bdt_et_bin_edges[ib] &&
                            oracle_et < bdt_et_bin_edges[ib + 1])
                        {
                            oracle_selected_model = bdt_et_bin_models[ib];
                            break;
                        }
                    }
                }
                const float oracle_base_e_score = (*bdt_arrays["base_E"])[icluster];
                const float oracle_base_v3e_score = (*bdt_arrays["base_v3E"])[icluster];
                const float oracle_selected_score =
                    (*bdt_arrays[oracle_selected_model])[icluster];
                const float oracle_lumi_weight = issim ? lumi / lumi_target : 1.0;
                const float oracle_trigger_weight =
                    event_weight != 0.0 ? weight / event_weight : 1.0;
                oracle_trace
                    << ientry << ',' << chain.GetTreeNumber() << ','
                    << chain.GetTree()->GetReadEntry() << ',' << *runnumber << ','
                    << *oracle_eventnumber << ',' << icluster << ','
                    << cluster_truthtrkID[icluster] << ',' << cluster_Et[icluster] << ','
                    << oracle_raw_recoisoET << ',' << recoisoET << ','
                    << recoiso_max << ',' << recononiso_min << ',' << sc.weight << ','
                    << mix_weight << ',' << oracle_lumi_weight << ',' << cross_weight << ','
                    << vertex_weight << ',' << oracle_truth_vertex_weight << ','
                    << oracle_trigger_weight << ',' << event_weight << ',' << weight << ','
                    << oracle_selected_model << ',' << oracle_base_e_score << ','
                    << oracle_base_v3e_score << ',' << oracle_selected_score << ','
                    << cluster_weta_cogx[icluster] << ','
                    << cluster_wphi_cogx[icluster] << ',' << *vertexz << ','
                    << cluster_Eta[icluster] << ',' << e11_over_e33 << ','
                    << cluster_et1[icluster] << ',' << cluster_et2[icluster] << ','
                    << cluster_et3[icluster] << ',' << cluster_et4[icluster] << ','
                    << e32_over_e35 << ',' << int(common_pass) << ',' << int(tight) << ','
                    << int(nontight) << ',' << int(iso) << ',' << int(noniso) << ','
                    << oracle_region << ',' << int(oracle_is_signal) << ','
                    << oracle_truth_index << ',' << oracle_truth_class << ','
                    << oracle_truth_pt << ',' << oracle_analysis_window << ','
                    << oracle_fill_a << ',' << oracle_fill_b << ',' << oracle_fill_c << ','
                    << oracle_fill_d << ',' << oracle_fill_multiplicity << '\n';
            }

            if (tight && iso)
'''
    apply(
        "            if (tight && iso)\n",
        trace_block,
        "candidate_trace_emit",
    )
    apply(
        "                    double cluster_Et_smear = cluster_Et[icluster]\n"
        "                        + ((sigma_extra_GeV > 0) ? rand->Gaus(0, sigma_extra_GeV) : 0.0);\n",
        "                    double cluster_Et_smear = cluster_Et[icluster]\n"
        "                        + ((sigma_extra_GeV > 0) ? rand->Gaus(0, sigma_extra_GeV) : 0.0);\n"
        "                    if (oracle_response_trace.is_open())\n"
        "                    {\n"
        "                        const int oracle_response_window =\n"
        "                            particle_Pt[iparticle] > pTmin_truth &&\n"
        "                            particle_Pt[iparticle] < pTmax_truth &&\n"
        "                            cluster_Et_smear > pTmin && cluster_Et_smear < pTmax;\n"
        "                        oracle_response_trace << ientry << ',' << chain.GetTreeNumber() << ','\n"
        "                            << chain.GetTree()->GetReadEntry() << ',' << *runnumber << ','\n"
        "                            << *oracle_eventnumber << ',' << icluster << ','\n"
        "                            << cluster_Et_smear << ',' << oracle_response_window << '\\n';\n"
        "                    }\n",
        "response_trace_emit",
    )
    return text, operations


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--source-revision", required=True)
    args = parser.parse_args()

    try:
        source = args.input.read_bytes()
        transformed, operations = instrument(source.decode())
        output = transformed.encode()
        if output == source:
            raise TransformFailure("instrumentation made no change")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.receipt.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_bytes(output)
        receipt = {
            "schema_version": SCHEMA_VERSION,
            "transform": TRANSFORM_NAME,
            "source_revision": args.source_revision,
            "input_path": str(args.input.resolve()),
            "input_sha256": digest_bytes(source),
            "output_path": str(args.output.resolve()),
            "output_sha256": digest_bytes(output),
            "operations": operations,
            "selection_or_fill_expression_replaced": False,
            "trace_side_channel_only": True,
        }
        args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    except (OSError, UnicodeDecodeError, TransformFailure) as exc:
        print(f"PPG12_RECOEFF_TRACE_INSTRUMENT_FAIL: {exc}", file=__import__("sys").stderr)
        return 2
    print(
        "PPG12_RECOEFF_TRACE_INSTRUMENT_PASS "
        f"output={args.output} receipt={args.receipt}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
