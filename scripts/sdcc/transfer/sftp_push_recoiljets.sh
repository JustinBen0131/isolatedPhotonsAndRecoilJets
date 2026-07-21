#!/usr/bin/env bash
set -euo pipefail

LOCAL_BASE="${RJ_SFTP_LOCAL_BASE:-/Users/patsfan753/Desktop/ThesisAnalysis}"
REMOTE_BASE="/sphenix/u/patsfan753/scratch/thesisAnalysis"
REMOTE_HOST="${RJ_SFTP_REMOTE_HOST:-patsfan753@sftp.sdcc.bnl.gov}"
REMOTE_TRANSPORT="${RJ_SFTP_TRANSPORT:-sftp}"
REMOTE_SSH_GATEWAY="${RJ_SDCC_SSH_GATEWAY:-patsfan753@ssh.sdcc.bnl.gov}"
REMOTE_SSH_TARGET="${RJ_SDCC_SSH_TARGET:-sphnxuser05.sdcc.bnl.gov}"

LOCAL_FILES=(
  "scripts/sdcc/runtime/audit/audit_auau_grl_projection.sh"
  "scripts/ml/audits/audit_auau_ml_training_smoke.py"
  "scripts/diagnostics/auau_split/audit_auau_truth_tags.py"
  "scripts/sdcc/pipelines/auau/auau_tight_bdt_pipeline.sh"
  "scripts/sdcc/pipelines/auau/auau_tight_logreg_pipeline.sh"
  "scripts/sdcc/pipelines/auau/auau_tight_mlp_pipeline.sh"
  "scripts/sdcc/pipelines/pp/pp_photon_ml_pipeline.sh"
  "scripts/sdcc/workflows/training/submit_the116_pp_bdt_15to35.sh"
  "scripts/data_prep/manifests/build_currentian_fast_manifests.py"
  "scripts/diagnostics/pp_shuhang/compare_pp_currentian_insitu_stitch_contract.py"
  "scripts/slides/pp_currentian/stitching/make_pp_currentian_insitu_contract_stitch_slides.py"
  "scripts/plotting/pp_currentian/make_ppg12_fig19_bdt_overlay.py"
  "scripts/diagnostics/pp_shuhang/audit_pp_basev3e_shuhang_equivalence.py"
  "scripts/plotting/pp_currentian/render_pp_currentian_shuhang_overlay_from_summary.py"
  "scripts/ml/validation/validate_pp_photon_ml_tables.py"
  "scripts/sdcc/workflows/submit/submit_auau_logreg_full_chain.sh"
  "scripts/sdcc/runtime/xsec/estimateEmbeddedPhotonXsec.sh"
  "scripts/data_prep/stitching/extract_focus21_fine_cluster_et_components.py"
  "scripts/slides/working_point/make_focus21_cluster_et_leakage_slide.py"
  "scripts/plotting/stitching/make_focus21_inclusive_reco_cluster_et_blair_plot.py"
  "scripts/sdcc/runtime/lists/make_dstListsData.sh"
  "scripts/sdcc/runtime/lists/makeThesisSimLists.sh"
  "scripts/sdcc/runtime/lists/makePPG12DoubleSimLists.sh"
  "scripts/sdcc/runtime/merge/mergeRecoilJets.sh"
  "scripts/sdcc/runtime/io/recoiljets_io_paths.sh"
  "scripts/sdcc/runtime/cleanup/recoiljets_cleanup.sh"
  "scripts/env/root_in_analysis_env.sh"
  "scripts/sdcc/workflows/width_study/submit_auau_bdt_widthstudy_pt1530_wp080.sh"
  "scripts/sdcc/workflows/width_study/submit_auau_bdt_widthstudy_windows_wp050.sh"
  "scripts/sdcc/workflows/width_study/merge_auau_bdt_widthstudy_windows_wp050_staged.sh"
  "scripts/sdcc/workflows/target_wp/merge_auau_bdt_target80_ready.sh"
  "scripts/sdcc/workflows/target_wp/submit_auau_bdt_etfine_centstudy_wp050.sh"
  "scripts/sdcc/workflows/target_wp/submit_auau_bdt_etfine_centstudy_target80.sh"
  "scripts/sdcc/workflows/target_wp/submit_auau_bdt_targetwp_pair.sh"
  "scripts/sdcc/workflows/target_wp/submit_auau_logreg_targetwp_pair.sh"
  "scripts/sdcc/workflows/target_wp/submit_auau_bdt_target80_config_dir.sh"
  "scripts/sdcc/workflows/target_wp/submit_auau_mlp_targetwp_pair.sh"
  "scripts/sdcc/workflows/submit/submit_the79_phenix_like_raa_campaign.sh"
  "scripts/diagnostics/ml_validation/check_the79_phenix_like_raa_inputs.py"
  "scripts/sdcc/workflows/submit/submit_auau_mlp_highpt_sweep.sh"
  "scripts/sdcc/workflows/submit/submit_auau_mlp_kitchensink.sh"
  "scripts/sdcc/workflows/submit/submit_auau_mlp_finept_distilled_sweep.sh"
  "scripts/sdcc/workflows/stacking/submit_auau_stacked_bdt_mlp_calibrator.sh"
  "scripts/sdcc/workflows/stacking/submit_auau_stacked_bdt_mlp_full_feature_chain.sh"
  "scripts/sdcc/workflows/stacking/submit_auau_stacked_bdt_mlp_sweep.sh"
  "scripts/sdcc/workflows/stacking/submit_auau_stack_matrix_wave.sh"
  "scripts/sdcc/workflows/diagnostics/submit_auau_iso_visible_diagnostic_chain.sh"
  "scripts/sdcc/workflows/diagnostics/submit_the38_tree_depth_capacity_campaign.sh"
  "scripts/sdcc/workflows/diagnostics/submit_the8_corrected_baseline_diagnostic_expansion.sh"
  "scripts/sdcc/workflows/stacking/submit_auau_oof_residual_superstacker.sh"
  "scripts/sdcc/workflows/stacking/submit_auau_global_sixpack_oof_stack.sh"
  "scripts/sdcc/workflows/stacking/submit_fresh_pp_auau_oof_stack_campaign.sh"
  "scripts/sdcc/workflows/stacking/submit_fresh_pp_auau_oof_stack_dag_campaign.sh"
  "scripts/sdcc/workflows/stacking/auau_bdt_mlp_stack_production_driver.sh"
  "scripts/ml/stacking/promote_auau_stacked_bdt_mlp.py"
  "scripts/plotting/auau_bdt/make_auau_bdt_training_closure.py"
  "scripts/plotting/auau_bdt/stacking/make_auau_bdt_mlp_stack_roc_overlay.py"
  "scripts/plotting/auau_bdt/stacking/make_auau_bdt_mlp_stack_score_separation.py"
  "scripts/plotting/auau_bdt/stacking/make_auau_mlp_expert_validation_plots.py"
  "scripts/plotting/auau_bdt/stacking/make_auau_mlp_training_curves.py"
  "scripts/plotting/auau_bdt/stacking/make_auau_stacked_training_curves.py"
  "scripts/diagnostics/ml_validation/make_auau_iso_visible_diagnostic_summary.py"
  "scripts/plotting/auau_bdt/make_auau_isolation_feature_correlations.py"
  "scripts/ml/stacking/train_auau_stacked_bdt_mlp_calibrator.py"
  "scripts/ml/stacking/train_auau_stacked_bdt_mlp_sweep.py"
  "scripts/ml/stacking/train_auau_oof_residual_superstacker.py"
  "scripts/ml/stacking/train_photon_bdt_mlp_oof_stack.py"
  "scripts/ml/stacking/train_photon_bdt_mlp_oof_stack_staged.py"
  "scripts/sdcc/workflows/stacking/auau_mlp_bdt_beating_driver.sh"
  "scripts/sdcc/workflows/target_wp/prepare_auau_bdt_target80_available_campaigns.sh"
  "scripts/ml/working_points/make_auau_bdt_target_wp_config.py"
  "scripts/ml/working_points/make_auau_logreg_target_wp_config.py"
  "scripts/ml/working_points/make_auau_mlp_target_wp_config.py"
  "scripts/ml/working_points/derive_the57_full_weighted_wp80.py"
  "scripts/ml/working_points/derive_the107_paired_wp.py"
  "scripts/ml/working_points/derive_corrected_auau_shower_wp.py"
  "scripts/ml/validation/validate_xgb_tmva_runtime_parity.py"
  "scripts/ml/validation/audit_the107_auau_bdt_extraction.py"
  "scripts/ml/validation/compare_the107_auau_bdt_label_contracts.py"
  "scripts/ml/validation/audit_corrected_auau_shower_contract_roots.py"
  "scripts/plotting/auau_bdt/make_corrected_shower_contract_comparison.py"
  "scripts/sdcc/workflows/validation/run_corrected_auau_shower_contract_audits.sh"
  "scripts/sdcc/runtime/condor/RecoilJets_Condor_AuAu.sh"
  "scripts/sdcc/runtime/condor/RecoilJets_Condor_submit.sh"
  "scripts/sdcc/runtime/condor/RecoilJets_Condor.sh"
  "scripts/ml/stacking/train_auau_jet_residual_bdt.py"
  "scripts/ml/training/train_auau_photon_bdt.py"
  "scripts/ml/training/train_auau_photon_logreg.py"
  "scripts/ml/training/train_auau_photon_mlp.py"
  "scripts/ml/validation/validate_auau_tight_bdt_on_sim.py"
  "scripts/ml/validation/validate_auau_tight_logreg_on_sim.py"
  "scripts/ml/validation/validate_auau_tight_mlp_on_sim.py"
  "macros/analysis_config.yaml"
  "macros/analysis_config_auau_bdt_validation.yaml"
  "macros/analysis_config_auau_bdt_validation_wp080.yaml"
  "macros/analysis_config_auau_bdt_validation_wp080_no3x3.yaml"
  "macros/analysis_config_the88_bounded_sideband_default_bdt.yaml"
  "macros/analysis_config_the88a_embedded_truthiso_diagnostics.yaml"
  "macros/analysis_config_the95_embedded_pmt_low_calo_diagnostics.yaml"
  "macros/analysis_config_the100_auau_dualview_comparison.yaml"
  "macros/analysis_config_the102_auau_fig25_correlations.yaml"
  "macros/analysis_config_the104_canonical_minbias_isolation_base.yaml"
  "macros/analysis_config_the104_canonical_minbias_isolation_phosub.yaml"
  "macros/analysis_config_auau_mlp_validation.yaml"
  "macros/analysis_config_auau_mlp_v2_validation.yaml"
  "macros/analysis_config_auau_bdt_widthstudy_pt1530_wp080.yaml"
  "macros/analysis_config_auau_bdt_widthstudy_pt1530_wp050.yaml"
  "macros/analysis_config_auau_bdt_etfine_centstudy_wp050.yaml"
  "macros/analysis_config_auau_bdt_mlp_stack_template.yaml"
  "macros/analysis_config_auau_bdt_base3x3_pt5to40_targetwp_template.yaml"
  "macros/analysis_config_the42_wp80_centlinear_ss_overlay.yaml"
  "macros/analysis_config_the79_phenix_like_raa_pp_ppg12.yaml"
  "macros/analysis_config_the79_phenix_like_raa_auau_bdt98_wp80.yaml"
  "scripts/sdcc/workflows/submit/submit_the96_auau_iso_baseline5pct_campaign.sh"
  "scripts/sdcc/workflows/diagnostics/submit_the95_embedded_pmt_low_calo_diagnostics.sh"
  "scripts/sdcc/workflows/submit/submit_the98_auau_emcal_phosub_iso_campaign.sh"
  "scripts/sdcc/workflows/diagnostics/submit_the100_auau_dualview_campaign.sh"
  "scripts/diagnostics/auau_bdt/validate_the100_dualview_root.py"
  "scripts/sdcc/workflows/diagnostics/submit_the102_auau_fig25_correlations.sh"
  "scripts/diagnostics/auau_bdt/validate_the102_fig25_root.py"
  "scripts/sdcc/workflows/diagnostics/submit_the104_canonical_minbias_isolation_slides.sh"
  "scripts/diagnostics/auau_bdt/validate_the104_isolation_root.py"
  "scripts/sdcc/workflows/submit/heartbeat_the96_auau_iso_baseline5pct_campaign.sh"
  "macros/Calo_Calib.C"
  "macros/Fun4All_recoilJets.C"
  "macros/Fun4All_recoilJets_AuAu.C"
  "macros/Fun4All_auauTightBDTTraining.C"
  "macros/Fun4All_recoilJets_unified_impl.C"
  "macros/diagnostics/stitching/PrintPPStitchDiagnostics.C"
  "src/PhotonClusterBuilder.cc"
  "src/PhotonClusterBuilder.h"
  "coresoftware_local/offline/packages/CaloReco/CaloTowerStatus.cc"
  "coresoftware_local/offline/packages/CaloReco/CaloTowerStatus.h"
  "coresoftware_local/offline/packages/CaloReco/BEmcRecCEMC.cc"
  "coresoftware_local/offline/packages/CaloReco/BEmcRecCEMC.h"
  "coresoftware_local/offline/packages/CaloReco/RawClusterBuilderTemplate.cc"
  "coresoftware_local/offline/packages/CaloReco/RawClusterBuilderTopo.cc"
  "coresoftware_local/offline/packages/CaloBase/RawTowerDefs.h"
  "src/RecoilJets.cc"
  "src/RecoilJets.h"
  "src/PPG12SimWeight.h"
  "src_AuAu/configure.ac"
  "src_AuAu/Makefile.am"
  "src_AuAu/RecoilJets_AuAu.cc"
  "src_AuAu/RecoilJets_AuAu.h"
  "scripts/sdcc/workflows/diagnostics/submit_the105_auau_shower_contract_factorial.sh"
  "scripts/sdcc/workflows/diagnostics/the105_preserve_invalid_shower_shapes.patch"
  "scripts/sdcc/workflows/diagnostics/the105_skip_invalid_rawcluster_tower_coordinates.patch"
  "macros/analysis_config_the112_auau_combined_bdt_triplet.yaml"
  "scripts/sdcc/workflows/diagnostics/submit_the112_auau_sideband_campaign.sh"
  "scripts/diagnostics/auau_bdt/rank_the112_sideband_scan.py"
  "scripts/diagnostics/auau_bdt/validate_the112_sideband_canary.py"
  "scripts/diagnostics/pp_currentian/compare_ppg12_recoiljets_same_cluster_features.py"
  "scripts/diagnostics/pp_currentian/run_ppg12_photon_oracle_canary_audit.py"
  "scripts/diagnostics/pp_currentian/instrument_ppg12_recoeff_trace.py"
  "scripts/diagnostics/pp_currentian/extract_ppg12_recoeff_executable_aggregate.py"
  "scripts/diagnostics/pp_currentian/produce_ppg12_stitched_purity_evidence.py"
  "scripts/diagnostics/pp_currentian/ppg12_stitched_purity_closure_gate.py"
  "scripts/diagnostics/pp_currentian/assemble_ppg12_stitched_purity_manifest.py"
  "scripts/diagnostics/pp_currentian/extract_ppg12_stitched_purity_lane.py"
  "scripts/data_prep/manifests/build_ppg12_paired_source_manifest.py"
  "scripts/sdcc/workflows/diagnostics/submit_ppg12_stitched_purity_photon_canaries.sh"
  "scripts/sdcc/workflows/diagnostics/ppg12_paired_oracle_condor_fanout.py"
  "scripts/sdcc/workflows/diagnostics/build_ppg12_oracle_new17_runtime.sh"
  "scripts/sdcc/workflows/diagnostics/run_ppg12_recoiljets_paired_oracle.sh"
  "scripts/sdcc/workflows/diagnostics/run_ppg12_recoiljets_paired_oracle_worker.sh"
  "scripts/diagnostics/pp_currentian/audit_ppg12_recoiljets_paired_oracle.py"
  "macros/diagnostics/pp_currentian/Fun4All_ppg12_fixed_seed_oracle.C"
  "macros/diagnostics/pp_currentian/Fun4All_recoiljets_fixed_seed_oracle.C"
  "src/configure.ac"
  "src/Makefile.am"
  "src/autogen.sh"
)

REMOTE_FILES=(
  "scripts/sdcc/runtime/audit/audit_auau_grl_projection.sh"
  "scripts/ml/audits/audit_auau_ml_training_smoke.py"
  "scripts/diagnostics/auau_split/audit_auau_truth_tags.py"
  "scripts/sdcc/pipelines/auau/auau_tight_bdt_pipeline.sh"
  "scripts/sdcc/pipelines/auau/auau_tight_logreg_pipeline.sh"
  "scripts/sdcc/pipelines/auau/auau_tight_mlp_pipeline.sh"
  "scripts/sdcc/pipelines/pp/pp_photon_ml_pipeline.sh"
  "scripts/sdcc/workflows/training/submit_the116_pp_bdt_15to35.sh"
  "scripts/data_prep/manifests/build_currentian_fast_manifests.py"
  "scripts/diagnostics/pp_shuhang/compare_pp_currentian_insitu_stitch_contract.py"
  "scripts/slides/pp_currentian/stitching/make_pp_currentian_insitu_contract_stitch_slides.py"
  "scripts/plotting/pp_currentian/make_ppg12_fig19_bdt_overlay.py"
  "scripts/diagnostics/pp_shuhang/audit_pp_basev3e_shuhang_equivalence.py"
  "scripts/plotting/pp_currentian/render_pp_currentian_shuhang_overlay_from_summary.py"
  "scripts/ml/validation/validate_pp_photon_ml_tables.py"
  "scripts/sdcc/workflows/submit/submit_auau_logreg_full_chain.sh"
  "scripts/sdcc/runtime/xsec/estimateEmbeddedPhotonXsec.sh"
  "scripts/data_prep/stitching/extract_focus21_fine_cluster_et_components.py"
  "scripts/slides/working_point/make_focus21_cluster_et_leakage_slide.py"
  "scripts/plotting/stitching/make_focus21_inclusive_reco_cluster_et_blair_plot.py"
  "scripts/sdcc/runtime/lists/make_dstListsData.sh"
  "scripts/sdcc/runtime/lists/makeThesisSimLists.sh"
  "scripts/sdcc/runtime/lists/makePPG12DoubleSimLists.sh"
  "scripts/sdcc/runtime/merge/mergeRecoilJets.sh"
  "scripts/sdcc/runtime/io/recoiljets_io_paths.sh"
  "scripts/sdcc/runtime/cleanup/recoiljets_cleanup.sh"
  "scripts/env/root_in_analysis_env.sh"
  "scripts/sdcc/workflows/width_study/submit_auau_bdt_widthstudy_pt1530_wp080.sh"
  "scripts/sdcc/workflows/width_study/submit_auau_bdt_widthstudy_windows_wp050.sh"
  "scripts/sdcc/workflows/width_study/merge_auau_bdt_widthstudy_windows_wp050_staged.sh"
  "scripts/sdcc/workflows/target_wp/merge_auau_bdt_target80_ready.sh"
  "scripts/sdcc/workflows/target_wp/submit_auau_bdt_etfine_centstudy_wp050.sh"
  "scripts/sdcc/workflows/target_wp/submit_auau_bdt_etfine_centstudy_target80.sh"
  "scripts/sdcc/workflows/target_wp/submit_auau_bdt_targetwp_pair.sh"
  "scripts/sdcc/workflows/target_wp/submit_auau_logreg_targetwp_pair.sh"
  "scripts/sdcc/workflows/target_wp/submit_auau_bdt_target80_config_dir.sh"
  "scripts/sdcc/workflows/target_wp/submit_auau_mlp_targetwp_pair.sh"
  "scripts/sdcc/workflows/submit/submit_the79_phenix_like_raa_campaign.sh"
  "scripts/diagnostics/ml_validation/check_the79_phenix_like_raa_inputs.py"
  "scripts/sdcc/workflows/submit/submit_auau_mlp_highpt_sweep.sh"
  "scripts/sdcc/workflows/submit/submit_auau_mlp_kitchensink.sh"
  "scripts/sdcc/workflows/submit/submit_auau_mlp_finept_distilled_sweep.sh"
  "scripts/sdcc/workflows/stacking/submit_auau_stacked_bdt_mlp_calibrator.sh"
  "scripts/sdcc/workflows/stacking/submit_auau_stacked_bdt_mlp_full_feature_chain.sh"
  "scripts/sdcc/workflows/stacking/submit_auau_stacked_bdt_mlp_sweep.sh"
  "scripts/sdcc/workflows/stacking/submit_auau_stack_matrix_wave.sh"
  "scripts/sdcc/workflows/diagnostics/submit_auau_iso_visible_diagnostic_chain.sh"
  "scripts/sdcc/workflows/diagnostics/submit_the38_tree_depth_capacity_campaign.sh"
  "scripts/sdcc/workflows/diagnostics/submit_the8_corrected_baseline_diagnostic_expansion.sh"
  "scripts/sdcc/workflows/stacking/submit_auau_oof_residual_superstacker.sh"
  "scripts/sdcc/workflows/stacking/submit_auau_global_sixpack_oof_stack.sh"
  "scripts/sdcc/workflows/stacking/submit_fresh_pp_auau_oof_stack_campaign.sh"
  "scripts/sdcc/workflows/stacking/submit_fresh_pp_auau_oof_stack_dag_campaign.sh"
  "scripts/sdcc/workflows/stacking/auau_bdt_mlp_stack_production_driver.sh"
  "scripts/ml/stacking/promote_auau_stacked_bdt_mlp.py"
  "scripts/plotting/auau_bdt/make_auau_bdt_training_closure.py"
  "scripts/plotting/auau_bdt/stacking/make_auau_bdt_mlp_stack_roc_overlay.py"
  "scripts/plotting/auau_bdt/stacking/make_auau_bdt_mlp_stack_score_separation.py"
  "scripts/plotting/auau_bdt/stacking/make_auau_mlp_expert_validation_plots.py"
  "scripts/plotting/auau_bdt/stacking/make_auau_mlp_training_curves.py"
  "scripts/plotting/auau_bdt/stacking/make_auau_stacked_training_curves.py"
  "scripts/diagnostics/ml_validation/make_auau_iso_visible_diagnostic_summary.py"
  "scripts/plotting/auau_bdt/make_auau_isolation_feature_correlations.py"
  "scripts/ml/stacking/train_auau_stacked_bdt_mlp_calibrator.py"
  "scripts/ml/stacking/train_auau_stacked_bdt_mlp_sweep.py"
  "scripts/ml/stacking/train_auau_oof_residual_superstacker.py"
  "scripts/ml/stacking/train_photon_bdt_mlp_oof_stack.py"
  "scripts/ml/stacking/train_photon_bdt_mlp_oof_stack_staged.py"
  "scripts/sdcc/workflows/stacking/auau_mlp_bdt_beating_driver.sh"
  "scripts/sdcc/workflows/target_wp/prepare_auau_bdt_target80_available_campaigns.sh"
  "scripts/ml/working_points/make_auau_bdt_target_wp_config.py"
  "scripts/ml/working_points/make_auau_logreg_target_wp_config.py"
  "scripts/ml/working_points/make_auau_mlp_target_wp_config.py"
  "scripts/ml/working_points/derive_the57_full_weighted_wp80.py"
  "scripts/ml/working_points/derive_the107_paired_wp.py"
  "scripts/ml/working_points/derive_corrected_auau_shower_wp.py"
  "scripts/ml/validation/validate_xgb_tmva_runtime_parity.py"
  "scripts/ml/validation/audit_the107_auau_bdt_extraction.py"
  "scripts/ml/validation/compare_the107_auau_bdt_label_contracts.py"
  "scripts/ml/validation/audit_corrected_auau_shower_contract_roots.py"
  "scripts/plotting/auau_bdt/make_corrected_shower_contract_comparison.py"
  "scripts/sdcc/workflows/validation/run_corrected_auau_shower_contract_audits.sh"
  "RecoilJets_Condor_AuAu.sh"
  "RecoilJets_Condor_submit.sh"
  "RecoilJets_Condor.sh"
  "scripts/ml/stacking/train_auau_jet_residual_bdt.py"
  "scripts/ml/training/train_auau_photon_bdt.py"
  "scripts/ml/training/train_auau_photon_logreg.py"
  "scripts/ml/training/train_auau_photon_mlp.py"
  "scripts/ml/validation/validate_auau_tight_bdt_on_sim.py"
  "scripts/ml/validation/validate_auau_tight_logreg_on_sim.py"
  "scripts/ml/validation/validate_auau_tight_mlp_on_sim.py"
  "macros/analysis_config.yaml"
  "macros/analysis_config_auau_bdt_validation.yaml"
  "macros/analysis_config_auau_bdt_validation_wp080.yaml"
  "macros/analysis_config_auau_bdt_validation_wp080_no3x3.yaml"
  "macros/analysis_config_the88_bounded_sideband_default_bdt.yaml"
  "macros/analysis_config_the88a_embedded_truthiso_diagnostics.yaml"
  "macros/analysis_config_the95_embedded_pmt_low_calo_diagnostics.yaml"
  "macros/analysis_config_the100_auau_dualview_comparison.yaml"
  "macros/analysis_config_the102_auau_fig25_correlations.yaml"
  "macros/analysis_config_the104_canonical_minbias_isolation_base.yaml"
  "macros/analysis_config_the104_canonical_minbias_isolation_phosub.yaml"
  "macros/analysis_config_auau_mlp_validation.yaml"
  "macros/analysis_config_auau_mlp_v2_validation.yaml"
  "macros/analysis_config_auau_bdt_widthstudy_pt1530_wp080.yaml"
  "macros/analysis_config_auau_bdt_widthstudy_pt1530_wp050.yaml"
  "macros/analysis_config_auau_bdt_etfine_centstudy_wp050.yaml"
  "macros/analysis_config_auau_bdt_mlp_stack_template.yaml"
  "macros/analysis_config_auau_bdt_base3x3_pt5to40_targetwp_template.yaml"
  "macros/analysis_config_the42_wp80_centlinear_ss_overlay.yaml"
  "macros/analysis_config_the79_phenix_like_raa_pp_ppg12.yaml"
  "macros/analysis_config_the79_phenix_like_raa_auau_bdt98_wp80.yaml"
  "scripts/sdcc/workflows/submit/submit_the96_auau_iso_baseline5pct_campaign.sh"
  "scripts/sdcc/workflows/diagnostics/submit_the95_embedded_pmt_low_calo_diagnostics.sh"
  "scripts/sdcc/workflows/submit/submit_the98_auau_emcal_phosub_iso_campaign.sh"
  "scripts/sdcc/workflows/diagnostics/submit_the100_auau_dualview_campaign.sh"
  "scripts/diagnostics/auau_bdt/validate_the100_dualview_root.py"
  "scripts/sdcc/workflows/diagnostics/submit_the102_auau_fig25_correlations.sh"
  "scripts/diagnostics/auau_bdt/validate_the102_fig25_root.py"
  "scripts/sdcc/workflows/diagnostics/submit_the104_canonical_minbias_isolation_slides.sh"
  "scripts/diagnostics/auau_bdt/validate_the104_isolation_root.py"
  "scripts/sdcc/workflows/submit/heartbeat_the96_auau_iso_baseline5pct_campaign.sh"
  "macros/Calo_Calib.C"
  "macros/Fun4All_recoilJets.C"
  "macros/Fun4All_recoilJets_AuAu.C"
  "macros/Fun4All_auauTightBDTTraining.C"
  "macros/Fun4All_recoilJets_unified_impl.C"
  "macros/diagnostics/stitching/PrintPPStitchDiagnostics.C"
  "coresoftware_local/offline/packages/CaloReco/PhotonClusterBuilder.cc"
  "coresoftware_local/offline/packages/CaloReco/PhotonClusterBuilder.h"
  "coresoftware_local/offline/packages/CaloReco/CaloTowerStatus.cc"
  "coresoftware_local/offline/packages/CaloReco/CaloTowerStatus.h"
  "coresoftware_local/offline/packages/CaloReco/BEmcRecCEMC.cc"
  "coresoftware_local/offline/packages/CaloReco/BEmcRecCEMC.h"
  "coresoftware_local/offline/packages/CaloReco/RawClusterBuilderTemplate.cc"
  "coresoftware_local/offline/packages/CaloReco/RawClusterBuilderTopo.cc"
  "coresoftware_local/offline/packages/CaloBase/RawTowerDefs.h"
  "src/RecoilJets.cc"
  "src/RecoilJets.h"
  "src/PPG12SimWeight.h"
  "src_AuAu/configure.ac"
  "src_AuAu/Makefile.am"
  "src_AuAu/RecoilJets_AuAu.cc"
  "src_AuAu/RecoilJets_AuAu.h"
  "scripts/sdcc/workflows/diagnostics/submit_the105_auau_shower_contract_factorial.sh"
  "scripts/sdcc/workflows/diagnostics/the105_preserve_invalid_shower_shapes.patch"
  "scripts/sdcc/workflows/diagnostics/the105_skip_invalid_rawcluster_tower_coordinates.patch"
  "macros/analysis_config_the112_auau_combined_bdt_triplet.yaml"
  "scripts/sdcc/workflows/diagnostics/submit_the112_auau_sideband_campaign.sh"
  "scripts/diagnostics/auau_bdt/rank_the112_sideband_scan.py"
  "scripts/diagnostics/auau_bdt/validate_the112_sideband_canary.py"
  "scripts/diagnostics/pp_currentian/compare_ppg12_recoiljets_same_cluster_features.py"
  "scripts/diagnostics/pp_currentian/run_ppg12_photon_oracle_canary_audit.py"
  "scripts/diagnostics/pp_currentian/instrument_ppg12_recoeff_trace.py"
  "scripts/diagnostics/pp_currentian/extract_ppg12_recoeff_executable_aggregate.py"
  "scripts/diagnostics/pp_currentian/produce_ppg12_stitched_purity_evidence.py"
  "scripts/diagnostics/pp_currentian/ppg12_stitched_purity_closure_gate.py"
  "scripts/diagnostics/pp_currentian/assemble_ppg12_stitched_purity_manifest.py"
  "scripts/diagnostics/pp_currentian/extract_ppg12_stitched_purity_lane.py"
  "scripts/data_prep/manifests/build_ppg12_paired_source_manifest.py"
  "scripts/sdcc/workflows/diagnostics/submit_ppg12_stitched_purity_photon_canaries.sh"
  "scripts/sdcc/workflows/diagnostics/ppg12_paired_oracle_condor_fanout.py"
  "scripts/sdcc/workflows/diagnostics/build_ppg12_oracle_new17_runtime.sh"
  "scripts/sdcc/workflows/diagnostics/run_ppg12_recoiljets_paired_oracle.sh"
  "scripts/sdcc/workflows/diagnostics/run_ppg12_recoiljets_paired_oracle_worker.sh"
  "scripts/diagnostics/pp_currentian/audit_ppg12_recoiljets_paired_oracle.py"
  "macros/diagnostics/pp_currentian/Fun4All_ppg12_fixed_seed_oracle.C"
  "macros/diagnostics/pp_currentian/Fun4All_recoiljets_fixed_seed_oracle.C"
  "src/configure.ac"
  "src/Makefile.am"
  "src/autogen.sh"
)

GROUP_CONDOR=(
  "scripts/sdcc/runtime/condor/RecoilJets_Condor_AuAu.sh"
  "scripts/sdcc/runtime/condor/RecoilJets_Condor_submit.sh"
  "scripts/sdcc/runtime/condor/RecoilJets_Condor.sh"
)

GROUP_MACROS=(
  "macros/analysis_config.yaml"
  "macros/analysis_config_auau_bdt_validation.yaml"
  "macros/analysis_config_auau_bdt_validation_wp080.yaml"
  "macros/analysis_config_auau_bdt_validation_wp080_no3x3.yaml"
  "macros/analysis_config_the88_bounded_sideband_default_bdt.yaml"
  "macros/analysis_config_the88a_embedded_truthiso_diagnostics.yaml"
  "macros/analysis_config_the95_embedded_pmt_low_calo_diagnostics.yaml"
  "macros/analysis_config_the100_auau_dualview_comparison.yaml"
  "macros/analysis_config_the102_auau_fig25_correlations.yaml"
  "macros/analysis_config_the104_canonical_minbias_isolation_base.yaml"
  "macros/analysis_config_the104_canonical_minbias_isolation_phosub.yaml"
  "macros/analysis_config_auau_mlp_validation.yaml"
  "macros/analysis_config_auau_mlp_v2_validation.yaml"
  "macros/analysis_config_auau_bdt_widthstudy_pt1530_wp080.yaml"
  "macros/analysis_config_auau_bdt_widthstudy_pt1530_wp050.yaml"
  "macros/analysis_config_auau_bdt_etfine_centstudy_wp050.yaml"
  "macros/analysis_config_auau_bdt_mlp_stack_template.yaml"
  "macros/analysis_config_auau_bdt_base3x3_pt5to40_targetwp_template.yaml"
  "macros/analysis_config_the42_wp80_centlinear_ss_overlay.yaml"
  "macros/analysis_config_the79_phenix_like_raa_pp_ppg12.yaml"
  "macros/analysis_config_the79_phenix_like_raa_auau_bdt98_wp80.yaml"
  "macros/analysis_config_the112_auau_combined_bdt_triplet.yaml"
  "macros/Calo_Calib.C"
  "macros/Fun4All_recoilJets.C"
  "macros/Fun4All_recoilJets_AuAu.C"
  "macros/Fun4All_auauTightBDTTraining.C"
  "macros/Fun4All_recoilJets_unified_impl.C"
  "macros/diagnostics/stitching/PrintPPStitchDiagnostics.C"
)

GROUP_SCRIPTS=(
  "scripts/sdcc/runtime/audit/audit_auau_grl_projection.sh"
  "scripts/ml/audits/audit_auau_ml_training_smoke.py"
  "scripts/diagnostics/auau_split/audit_auau_truth_tags.py"
  "scripts/sdcc/pipelines/auau/auau_tight_bdt_pipeline.sh"
  "scripts/sdcc/pipelines/auau/auau_tight_logreg_pipeline.sh"
  "scripts/sdcc/pipelines/auau/auau_tight_mlp_pipeline.sh"
  "scripts/sdcc/pipelines/pp/pp_photon_ml_pipeline.sh"
  "scripts/data_prep/manifests/build_currentian_fast_manifests.py"
  "scripts/diagnostics/pp_shuhang/compare_pp_currentian_insitu_stitch_contract.py"
  "scripts/slides/pp_currentian/stitching/make_pp_currentian_insitu_contract_stitch_slides.py"
  "scripts/plotting/pp_currentian/make_ppg12_fig19_bdt_overlay.py"
  "scripts/diagnostics/pp_shuhang/audit_pp_basev3e_shuhang_equivalence.py"
  "scripts/plotting/pp_currentian/render_pp_currentian_shuhang_overlay_from_summary.py"
  "scripts/ml/validation/validate_pp_photon_ml_tables.py"
  "scripts/sdcc/workflows/submit/submit_auau_logreg_full_chain.sh"
  "scripts/sdcc/runtime/xsec/estimateEmbeddedPhotonXsec.sh"
  "scripts/data_prep/stitching/extract_focus21_fine_cluster_et_components.py"
  "scripts/slides/working_point/make_focus21_cluster_et_leakage_slide.py"
  "scripts/plotting/stitching/make_focus21_inclusive_reco_cluster_et_blair_plot.py"
  "scripts/sdcc/runtime/lists/make_dstListsData.sh"
  "scripts/sdcc/runtime/lists/makeThesisSimLists.sh"
  "scripts/sdcc/runtime/lists/makePPG12DoubleSimLists.sh"
  "scripts/sdcc/runtime/merge/mergeRecoilJets.sh"
  "scripts/sdcc/runtime/io/recoiljets_io_paths.sh"
  "scripts/sdcc/runtime/cleanup/recoiljets_cleanup.sh"
  "scripts/env/root_in_analysis_env.sh"
  "scripts/sdcc/workflows/width_study/submit_auau_bdt_widthstudy_pt1530_wp080.sh"
  "scripts/sdcc/workflows/width_study/submit_auau_bdt_widthstudy_windows_wp050.sh"
  "scripts/sdcc/workflows/width_study/merge_auau_bdt_widthstudy_windows_wp050_staged.sh"
  "scripts/sdcc/workflows/target_wp/merge_auau_bdt_target80_ready.sh"
  "scripts/sdcc/workflows/target_wp/submit_auau_bdt_etfine_centstudy_wp050.sh"
  "scripts/sdcc/workflows/target_wp/submit_auau_bdt_etfine_centstudy_target80.sh"
  "scripts/sdcc/workflows/target_wp/submit_auau_bdt_targetwp_pair.sh"
  "scripts/sdcc/workflows/target_wp/submit_auau_logreg_targetwp_pair.sh"
  "scripts/sdcc/workflows/target_wp/submit_auau_bdt_target80_config_dir.sh"
  "scripts/sdcc/workflows/target_wp/submit_auau_mlp_targetwp_pair.sh"
  "scripts/sdcc/workflows/submit/submit_the79_phenix_like_raa_campaign.sh"
  "scripts/diagnostics/ml_validation/check_the79_phenix_like_raa_inputs.py"
  "scripts/sdcc/workflows/submit/submit_auau_mlp_highpt_sweep.sh"
  "scripts/sdcc/workflows/submit/submit_auau_mlp_kitchensink.sh"
  "scripts/sdcc/workflows/submit/submit_auau_mlp_finept_distilled_sweep.sh"
  "scripts/sdcc/workflows/diagnostics/submit_the95_embedded_pmt_low_calo_diagnostics.sh"
  "scripts/sdcc/workflows/diagnostics/submit_the100_auau_dualview_campaign.sh"
  "scripts/diagnostics/auau_bdt/validate_the100_dualview_root.py"
  "scripts/sdcc/workflows/diagnostics/submit_the102_auau_fig25_correlations.sh"
  "scripts/diagnostics/auau_bdt/validate_the102_fig25_root.py"
  "scripts/sdcc/workflows/diagnostics/submit_the104_canonical_minbias_isolation_slides.sh"
  "scripts/diagnostics/auau_bdt/validate_the104_isolation_root.py"
  "scripts/sdcc/workflows/diagnostics/submit_the105_auau_shower_contract_factorial.sh"
  "scripts/sdcc/workflows/diagnostics/the105_preserve_invalid_shower_shapes.patch"
  "scripts/sdcc/workflows/diagnostics/the105_skip_invalid_rawcluster_tower_coordinates.patch"
  "scripts/sdcc/workflows/diagnostics/submit_the112_auau_sideband_campaign.sh"
  "scripts/diagnostics/auau_bdt/rank_the112_sideband_scan.py"
  "scripts/diagnostics/auau_bdt/validate_the112_sideband_canary.py"
  "scripts/sdcc/workflows/stacking/submit_auau_stacked_bdt_mlp_calibrator.sh"
  "scripts/sdcc/workflows/stacking/submit_auau_stacked_bdt_mlp_full_feature_chain.sh"
  "scripts/sdcc/workflows/stacking/submit_auau_stacked_bdt_mlp_sweep.sh"
  "scripts/sdcc/workflows/stacking/submit_auau_stack_matrix_wave.sh"
  "scripts/sdcc/workflows/diagnostics/submit_auau_iso_visible_diagnostic_chain.sh"
  "scripts/sdcc/workflows/diagnostics/submit_the38_tree_depth_capacity_campaign.sh"
  "scripts/sdcc/workflows/diagnostics/submit_the8_corrected_baseline_diagnostic_expansion.sh"
  "scripts/sdcc/workflows/stacking/submit_auau_oof_residual_superstacker.sh"
  "scripts/sdcc/workflows/stacking/submit_auau_global_sixpack_oof_stack.sh"
  "scripts/sdcc/workflows/stacking/submit_fresh_pp_auau_oof_stack_campaign.sh"
  "scripts/sdcc/workflows/stacking/auau_bdt_mlp_stack_production_driver.sh"
  "scripts/ml/stacking/promote_auau_stacked_bdt_mlp.py"
  "scripts/plotting/auau_bdt/make_auau_bdt_training_closure.py"
  "scripts/plotting/auau_bdt/stacking/make_auau_bdt_mlp_stack_roc_overlay.py"
  "scripts/plotting/auau_bdt/stacking/make_auau_bdt_mlp_stack_score_separation.py"
  "scripts/plotting/auau_bdt/stacking/make_auau_mlp_expert_validation_plots.py"
  "scripts/plotting/auau_bdt/stacking/make_auau_mlp_training_curves.py"
  "scripts/plotting/auau_bdt/stacking/make_auau_stacked_training_curves.py"
  "scripts/diagnostics/ml_validation/make_auau_iso_visible_diagnostic_summary.py"
  "scripts/plotting/auau_bdt/make_auau_isolation_feature_correlations.py"
  "scripts/ml/stacking/train_auau_stacked_bdt_mlp_calibrator.py"
  "scripts/ml/stacking/train_auau_stacked_bdt_mlp_sweep.py"
  "scripts/ml/stacking/train_auau_oof_residual_superstacker.py"
  "scripts/ml/stacking/train_photon_bdt_mlp_oof_stack.py"
  "scripts/sdcc/workflows/stacking/auau_mlp_bdt_beating_driver.sh"
  "scripts/sdcc/workflows/target_wp/prepare_auau_bdt_target80_available_campaigns.sh"
  "scripts/ml/working_points/make_auau_bdt_target_wp_config.py"
  "scripts/ml/working_points/make_auau_logreg_target_wp_config.py"
  "scripts/ml/working_points/make_auau_mlp_target_wp_config.py"
  "scripts/ml/working_points/derive_the57_full_weighted_wp80.py"
  "scripts/ml/working_points/derive_the107_paired_wp.py"
  "scripts/ml/working_points/derive_corrected_auau_shower_wp.py"
  "scripts/ml/validation/validate_xgb_tmva_runtime_parity.py"
  "scripts/ml/validation/audit_the107_auau_bdt_extraction.py"
  "scripts/ml/validation/compare_the107_auau_bdt_label_contracts.py"
  "scripts/ml/validation/audit_corrected_auau_shower_contract_roots.py"
  "scripts/plotting/auau_bdt/make_corrected_shower_contract_comparison.py"
  "scripts/sdcc/workflows/validation/run_corrected_auau_shower_contract_audits.sh"
  "scripts/sdcc/runtime/condor/RecoilJets_Condor_AuAu.sh"
  "scripts/sdcc/runtime/condor/RecoilJets_Condor_submit.sh"
  "scripts/sdcc/runtime/condor/RecoilJets_Condor.sh"
  "scripts/ml/stacking/train_auau_jet_residual_bdt.py"
  "scripts/ml/training/train_auau_photon_bdt.py"
  "scripts/ml/training/train_auau_photon_logreg.py"
  "scripts/ml/training/train_auau_photon_mlp.py"
  "scripts/ml/validation/validate_auau_tight_bdt_on_sim.py"
  "scripts/ml/validation/validate_auau_tight_logreg_on_sim.py"
  "scripts/ml/validation/validate_auau_tight_mlp_on_sim.py"
)

usage() {
  cat <<'EOF'
Usage:
  ./scripts/sftp_push_recoiljets.sh <group-or-file> [group-or-file ...] [--commit-push -m "message"]
  ./scripts/sftp_push_recoiljets.sh status [group-or-file ...]
  ./scripts/sftp_push_recoiljets.sh diff [group-or-file ...]

Uploads selected known files from the local Mac checkout to the SDCC analysis
checkout via interactive sftp. No password is stored; sftp prompts normally.
Upload mode prints the overwrite preview and then starts sftp directly; it does
not ask for an extra y/N confirmation.
Set RJ_SFTP_LOCAL_BASE to an absolute clean worktree path when a focused branch
must be deployed without reading files from the default live checkout.
If the SFTP endpoint rejects the available key but the SDCC SSH gateway works,
set RJ_SFTP_TRANSPORT=ssh-tar. That transport still uses this mapped-file
allowlist and path validation, then streams an exact tar payload through the
approved SSH gateway.

Local organization note:
  Top-level scripts are hard command aliases after THE-23 stage 5. This helper
  uses canonical local paths generated from scripts/sdcc/TRANSFER_MAP.tsv and
  also resolves retired scripts/foo names through scripts/compat/local. Script
  uploads target canonical SDCC subfolder paths where applicable; hard Condor
  runtime aliases still upload to their protected SDCC runtime paths.

Read-only modes:
  status
      Fetch selected mapped remote files into a temp directory and print whether
      each SDCC file matches the local file. With no selection, checks pipeline.
  diff
      Fetch selected mapped remote files and print unified diffs for files that
      differ. With no selection, checks pipeline.

Groups:
  condor    RecoilJets_Condor*.sh submit/wrapper files
  macros    known pipeline macros/configs
  scripts   known pipeline helper scripts
  pipeline  all known transferable pipeline files except local-only sftp helpers
  all       alias for pipeline
  changed   all changed known transferable files in this checkout

Examples:
  ./scripts/sftp_push_recoiljets.sh status
  ./scripts/sftp_push_recoiljets.sh status changed
  ./scripts/sftp_push_recoiljets.sh status pipeline
  ./scripts/sftp_push_recoiljets.sh diff RecoilJets_Condor_submit.sh
  ./scripts/sftp_push_recoiljets.sh RecoilJets.cc RecoilJets.h
  ./scripts/sftp_push_recoiljets.sh condor
  ./scripts/sftp_push_recoiljets.sh changed
  ./scripts/sftp_push_recoiljets.sh mergeRecoilJets.sh analysis_config.yaml RecoilJets_AuAu.cc
  ./scripts/sftp_push_recoiljets.sh scripts/mergeRecoilJets.sh
  ./scripts/sftp_push_recoiljets.sh RecoilJets.cc RecoilJets.h RecoilJets_Condor_submit.sh --commit-push -m "Update PP recoil jet stitching diagnostics"

Options:
  --commit-push
      After a successful SFTP upload, run git add/commit/push for exactly the
      selected transferable local files. Never uses git add .
  -m, --message
      Commit message required with --commit-push.

Known files (canonical local path -> SDCC upload target; see scripts/sdcc/TRANSFER_MAP.tsv for aliases and canonical SDCC paths):
EOF
  local i
  for (( i=0; i<${#LOCAL_FILES[@]}; i++ )); do
    printf '  %-45s -> %s\n' "${LOCAL_FILES[$i]}" "${REMOTE_FILES[$i]}"
  done
}

die_with_usage() {
  echo "[ERROR] $*" >&2
  echo >&2
  usage >&2
  exit 2
}

normalize_arg() {
  local x="$1"
  x="${x#./}"
  echo "$x"
}

validate_remote_base() {
  local path="$1"
  case "$path" in
    ""|*[[:space:]]*|*agent_context*|*.codex*|*codex*|*THE-*)
      echo "[ERROR] Unsafe SDCC remote base path: ${path}" >&2
      exit 2
      ;;
  esac
}

validate_remote_relative_path() {
  local path="$1"
  local label="${2:-remote path}"
  case "$path" in
    ""|*[[:space:]]*|/*|.|..|../*|*/../*|*/..|*//*|*agent_context*|*.codex*|*codex*|*THE-*)
      echo "[ERROR] Unsafe ${label}: ${path}" >&2
      echo "[ERROR] Refusing empty, absolute, whitespace, traversal, agent-context, codex, or THE-* SDCC targets." >&2
      exit 2
      ;;
  esac
}

validate_selected_remote_paths() {
  local remote
  validate_remote_base "$REMOTE_BASE"
  for remote in "${selected_remote[@]}"; do
    validate_remote_relative_path "$remote" "selected SDCC upload target"
  done
}

real_abs_path() {
  python3 -c 'import pathlib, sys; print(pathlib.Path(sys.argv[1]).resolve())' "$1"
}

real_rel_to_local_base() {
  local abs="$1"
  python3 - "$LOCAL_BASE" "$abs" <<'PY'
import pathlib
import sys

base = pathlib.Path(sys.argv[1]).resolve()
path = pathlib.Path(sys.argv[2]).resolve()
try:
    print(path.relative_to(base))
except ValueError:
    print(path)
PY
}

upload_source_path() {
  local rel="$1"
  local abs="${LOCAL_BASE}/${rel}"
  if [[ -L "$abs" ]]; then
    real_abs_path "$abs"
  else
    printf '%s\n' "$rel"
  fi
}

stage_paths_for_selected() {
  local rel abs target target_rel existing
  local -a out=()
  for rel in "${selected_local[@]}"; do
    out+=( "$rel" )
    abs="${LOCAL_BASE}/${rel}"
    if [[ -L "$abs" ]]; then
      target="$(real_abs_path "$abs")"
      target_rel="$(real_rel_to_local_base "$target")"
      out+=( "$target_rel" )
    fi
  done

  local -a unique=()
  for rel in "${out[@]}"; do
    for existing in "${unique[@]:-}"; do
      [[ "$existing" == "$rel" ]] && continue 2
    done
    unique+=( "$rel" )
  done
  printf '%s\n' "${unique[@]}"
}

run_commit_push() {
  echo
  echo "Git commit/push requested."
  echo "Repository : ${LOCAL_BASE}"
  echo "Message    : ${commit_message}"
  echo
  echo "Current git status:"
  (cd "$LOCAL_BASE" && git status --short)
  echo
  echo "Staging exactly the uploaded local files and their canonical symlink targets:"
  local f
  local -a stage_files=()
  while IFS= read -r f; do
    [[ -n "$f" ]] || continue
    stage_files+=( "$f" )
    echo "  $f"
  done < <(stage_paths_for_selected)

  (cd "$LOCAL_BASE" && git add -- "${stage_files[@]}")

  if (cd "$LOCAL_BASE" && git diff --cached --quiet -- "${stage_files[@]}"); then
    echo
    echo "[OK] No staged changes in selected uploaded files; skipping git commit/push."
    return 0
  fi

  echo
  echo "Staged diff summary:"
  (cd "$LOCAL_BASE" && git diff --cached --stat -- "${stage_files[@]}")
  echo
  (cd "$LOCAL_BASE" && git commit -m "$commit_message" -- "${stage_files[@]}")
  (cd "$LOCAL_BASE" && git push)
  echo
  echo "[OK] Git commit and push complete."
}

find_index_for_local() {
  local rel="$1"
  local i
  for (( i=0; i<${#LOCAL_FILES[@]}; i++ )); do
    if [[ "${LOCAL_FILES[$i]}" == "$rel" ]]; then
      echo "$i"
      return 0
    fi
  done

  if [[ "$rel" == macros/* ]]; then
    local macro_matches=()
    for (( i=0; i<${#LOCAL_FILES[@]}; i++ )); do
      if [[ "${LOCAL_FILES[$i]}" == macros/* && "$(basename "${LOCAL_FILES[$i]}")" == "$(basename "$rel")" ]]; then
        macro_matches+=( "$i" )
      fi
    done
    if (( ${#macro_matches[@]} == 1 )); then
      echo "${macro_matches[0]}"
      return 0
    fi
  fi

  local arg_abs="${LOCAL_BASE}/${rel}"
  local compat_abs=""
  if [[ "$rel" == scripts/* ]]; then
    compat_abs="${LOCAL_BASE}/scripts/compat/local/$(basename "$rel")"
  fi
  if [[ -e "$arg_abs" || -L "$arg_abs" || -e "$compat_abs" || -L "$compat_abs" ]]; then
    local arg_real local_abs local_real
    if [[ -e "$arg_abs" || -L "$arg_abs" ]]; then
      arg_real="$(real_abs_path "$arg_abs")"
    else
      arg_real="$(real_abs_path "$compat_abs")"
    fi
    for (( i=0; i<${#LOCAL_FILES[@]}; i++ )); do
      local_abs="${LOCAL_BASE}/${LOCAL_FILES[$i]}"
      if [[ -e "$local_abs" || -L "$local_abs" ]]; then
        local_real="$(real_abs_path "$local_abs")"
        if [[ "$local_real" == "$arg_real" ]]; then
          echo "$i"
          return 0
        fi
      fi
    done
  fi
  return 1
}

selected_local=()
selected_remote=()
selection_args=()
commit_push=0
commit_message=""

add_index() {
  local idx="$1"
  local rel="${LOCAL_FILES[$idx]}"
  local remote="${REMOTE_FILES[$idx]}"
  local existing
  for existing in "${selected_local[@]:-}"; do
    if [[ "$existing" == "$rel" ]]; then
      return 0
    fi
  done
  selected_local+=( "$rel" )
  selected_remote+=( "$remote" )
}

add_local_rel() {
  local rel="$1"
  local idx
  idx="$(find_index_for_local "$rel")" || die_with_usage "Unknown file in group definition: $rel"
  add_index "$idx"
}

add_changed_known_files() {
  local line path idx any=0
  while IFS= read -r line; do
    [[ -n "$line" ]] || continue
    path="${line:3}"
    if [[ "$path" == *" -> "* ]]; then
      path="${path##* -> }"
    fi
    path="${path#\"}"
    path="${path%\"}"
    path="$(normalize_arg "$path")"
    if idx="$(find_index_for_local "$path")"; then
      add_index "$idx"
      any=1
    fi
  done < <(cd "$LOCAL_BASE" && git status --porcelain --untracked-files=all)

  (( any )) || die_with_usage "No changed known transferable files found in ${LOCAL_BASE}."
}

add_group() {
  local group="$1"
  local f
  case "$group" in
    condor)
      for f in "${GROUP_CONDOR[@]}"; do add_local_rel "$f"; done
      ;;
    macros)
      for f in "${GROUP_MACROS[@]}"; do add_local_rel "$f"; done
      ;;
    scripts)
      for f in "${GROUP_SCRIPTS[@]}"; do add_local_rel "$f"; done
      ;;
    pipeline|all)
      add_local_rel "src/PhotonClusterBuilder.cc"
      add_local_rel "src/PhotonClusterBuilder.h"
      add_local_rel "coresoftware_local/offline/packages/CaloReco/CaloTowerStatus.cc"
      add_local_rel "coresoftware_local/offline/packages/CaloReco/CaloTowerStatus.h"
      add_local_rel "coresoftware_local/offline/packages/CaloReco/BEmcRecCEMC.cc"
      add_local_rel "coresoftware_local/offline/packages/CaloReco/BEmcRecCEMC.h"
      add_local_rel "coresoftware_local/offline/packages/CaloReco/RawClusterBuilderTemplate.cc"
      add_local_rel "coresoftware_local/offline/packages/CaloReco/RawClusterBuilderTopo.cc"
      add_local_rel "coresoftware_local/offline/packages/CaloBase/RawTowerDefs.h"
      add_local_rel "src/RecoilJets.cc"
      add_local_rel "src/RecoilJets.h"
      add_local_rel "src/PPG12SimWeight.h"
      add_local_rel "src_AuAu/RecoilJets_AuAu.cc"
      add_local_rel "src_AuAu/RecoilJets_AuAu.h"
      add_group condor
      add_group macros
      add_group scripts
      ;;
    changed|changed-pipeline|changedPipeline|changed-known)
      add_changed_known_files
      ;;
    *)
      return 1
      ;;
  esac
}

resolve_file_arg() {
  local arg
  arg="$(normalize_arg "$1")"

  local idx
  if idx="$(find_index_for_local "$arg")"; then
    add_index "$idx"
    return 0
  fi

  if [[ "$arg" == */* ]]; then
    die_with_usage "Unknown relative path: $arg"
  fi

  local matches=()
  local i
  for (( i=0; i<${#LOCAL_FILES[@]}; i++ )); do
    if [[ "$(basename "${LOCAL_FILES[$i]}")" == "$arg" ]]; then
      matches+=( "$i" )
    fi
  done

  if (( ${#matches[@]} == 0 )); then
    die_with_usage "Unknown group, relative path, or basename: $arg"
  fi
  if (( ${#matches[@]} > 1 )); then
    echo "[ERROR] Ambiguous basename: $arg" >&2
    echo "Use one of these relative paths:" >&2
    for idx in "${matches[@]}"; do
      echo "  ${LOCAL_FILES[$idx]}" >&2
    done
    exit 2
  fi

  add_index "${matches[0]}"
}

hash_file() {
  local file="$1"
  if command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$file" | awk '{print $1}'
  elif command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$file" | awk '{print $1}'
  else
    echo "[ERROR] Need shasum or sha256sum for status/diff mode." >&2
    exit 5
  fi
}

run_remote_compare() {
  local compare_mode="$1"
  local tmp_dir
  tmp_dir="$(mktemp -d "${TMPDIR:-/tmp}/sftp_push_recoiljets_check.XXXXXX")"
  local batch
  batch="${tmp_dir}/fetch.sftp"

  {
    printf 'cd %s\n' "$REMOTE_BASE"
    local i
    for (( i=0; i<${#selected_remote[@]}; i++ )); do
      printf -- '-get %s %s\n' "${selected_remote[$i]}" "${tmp_dir}/remote_${i}"
    done
  } > "$batch"

  echo
  echo "Remote host : ${REMOTE_HOST}"
  echo "Local base  : ${LOCAL_BASE}"
  echo "Remote base : ${REMOTE_BASE}"
  echo
  echo "Read-only ${compare_mode}: fetching ${#selected_local[@]} mapped file(s) into a temp directory."
  echo "No remote files will be modified."
  echo "Opening interactive sftp. Enter your SDCC password when prompted."

  if sftp \
      -oBatchMode=no \
      -oPreferredAuthentications=password,keyboard-interactive,publickey \
      -b "$batch" \
      "$REMOTE_HOST"; then
    :
  else
    local status=$?
    rm -rf "$tmp_dir"
    echo
    echo "[ERROR] sftp fetch failed with exit code ${status}." >&2
    exit "$status"
  fi

  echo
  printf '%-8s  %-45s  %s\n' "STATUS" "LOCAL" "REMOTE"
  printf '%-8s  %-45s  %s\n' "--------" "---------------------------------------------" "---------------------------------------------"

  local matches=0 differs=0 missing=0 checked=0
  local i local_file remote_copy local_hash remote_hash result
  for (( i=0; i<${#selected_local[@]}; i++ )); do
    local_file="${LOCAL_BASE}/${selected_local[$i]}"
    remote_copy="${tmp_dir}/remote_${i}"
    checked=$(( checked + 1 ))
    if [[ ! -f "$remote_copy" ]]; then
      result="MISSING"
      missing=$(( missing + 1 ))
    else
      local_hash="$(hash_file "$local_file")"
      remote_hash="$(hash_file "$remote_copy")"
      if [[ "$local_hash" == "$remote_hash" ]]; then
        result="MATCH"
        matches=$(( matches + 1 ))
      else
        result="DIFFER"
        differs=$(( differs + 1 ))
      fi
    fi
    printf '%-8s  %-45s  %s\n' "$result" "${selected_local[$i]}" "${selected_remote[$i]}"

    if [[ "$compare_mode" == "diff" && "$result" == "DIFFER" ]]; then
      echo
      echo "Diff for ${selected_local[$i]}:"
      diff -u "$local_file" "$remote_copy" | sed \
        -e "1s|.*|--- local:${selected_local[$i]}|" \
        -e "2s|.*|+++ sdcc:${selected_remote[$i]}|" || true
      echo
    elif [[ "$compare_mode" == "diff" && "$result" == "MISSING" ]]; then
      echo
      echo "Remote file missing for ${selected_local[$i]} -> ${selected_remote[$i]}"
      echo
    fi
  done

  rm -rf "$tmp_dir"
  echo
  echo "Summary: checked=${checked} match=${matches} differ=${differs} missing=${missing}"
  if (( differs > 0 || missing > 0 )); then
    exit 1
  fi
}

run_ssh_tar_upload() {
  local stage_dir tar_file
  stage_dir="$(mktemp -d "${TMPDIR:-/tmp}/sftp_push_recoiljets_stage.XXXXXX")"
  tar_file="$(mktemp "${TMPDIR:-/tmp}/sftp_push_recoiljets_payload.XXXXXX.tar")"
  local i src dst dst_dir
  for (( i=0; i<${#selected_local[@]}; i++ )); do
    src="${LOCAL_BASE}/${selected_local[$i]}"
    dst="${stage_dir}/${selected_remote[$i]}"
    dst_dir="$(dirname "$dst")"
    mkdir -p "$dst_dir"
    cp -p "$(real_abs_path "$src")" "$dst"
  done

  (cd "$stage_dir" && tar -cf "$tar_file" .)

  if [[ -z "${SSH_AUTH_SOCK:-}" ]] && command -v launchctl >/dev/null 2>&1; then
    export SSH_AUTH_SOCK="$(launchctl getenv SSH_AUTH_SOCK 2>/dev/null || true)"
  fi

  echo
  echo "Using RJ_SFTP_TRANSPORT=ssh-tar"
  echo "Gateway     : ${REMOTE_SSH_GATEWAY}"
  echo "Target host : ${REMOTE_SSH_TARGET}"
  echo "Remote base : ${REMOTE_BASE}"
  echo "Payload tar : ${tar_file}"

  if ssh "${REMOTE_SSH_GATEWAY}" \
      "ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null ${REMOTE_SSH_TARGET} 'mkdir -p ${REMOTE_BASE} && cd ${REMOTE_BASE} && tar -xf -'" \
      < "$tar_file"; then
    rm -rf "$stage_dir" "$tar_file"
    return 0
  fi

  local status=$?
  rm -rf "$stage_dir" "$tar_file"
  return "$status"
}

if (( $# == 0 )); then
  usage
  exit 2
fi

case "${1:-}" in
  -h|--help|help)
    usage
    exit 0
    ;;
esac

mode="upload"
case "${1:-}" in
  status|check)
    mode="status"
    shift
    ;;
  diff)
    mode="diff"
    shift
    ;;
  upload|push)
    mode="upload"
    shift
    ;;
esac

while (( $# > 0 )); do
  case "$1" in
    --commit-push)
      commit_push=1
      shift
      ;;
    -m|--message)
      opt="$1"
      shift
      if (( $# == 0 )); then
        die_with_usage "$opt requires a commit message"
      fi
      commit_message="$1"
      shift
      ;;
    --)
      shift
      while (( $# > 0 )); do
        selection_args+=( "$1" )
        shift
      done
      ;;
    -*)
      die_with_usage "Unknown option: $1"
      ;;
    *)
      selection_args+=( "$1" )
      shift
      ;;
  esac
done

if [[ "$mode" != "upload" && "$commit_push" -eq 1 ]]; then
  die_with_usage "--commit-push is only valid for upload mode"
fi

if (( commit_push )) && [[ -z "$commit_message" ]]; then
  die_with_usage "--commit-push requires -m \"commit message\""
fi

if (( ${#selection_args[@]} == 0 )) && [[ "$mode" != "upload" ]]; then
  selection_args+=( "pipeline" )
fi

if (( ${#selection_args[@]} == 0 )); then
  die_with_usage "No groups or files selected."
fi

for arg in "${selection_args[@]}"; do
  arg="$(normalize_arg "$arg")"
  if add_group "$arg"; then
    continue
  fi
  resolve_file_arg "$arg"
done

if (( ${#selected_local[@]} == 0 )); then
  die_with_usage "No files selected."
fi

validate_selected_remote_paths

if [[ ! -d "$LOCAL_BASE" ]]; then
  echo "[ERROR] Local base does not exist: $LOCAL_BASE" >&2
  exit 3
fi

for f in "${selected_local[@]}"; do
  if [[ ! -f "${LOCAL_BASE}/${f}" ]]; then
    echo "[ERROR] Missing local file: ${LOCAL_BASE}/${f}" >&2
    exit 4
  fi
done

if [[ "$mode" == "status" || "$mode" == "diff" ]]; then
  run_remote_compare "$mode"
  exit 0
fi

echo
echo "Remote host : ${REMOTE_HOST}"
echo "Local base  : ${LOCAL_BASE}"
echo "Remote base : ${REMOTE_BASE}"
echo
echo "Files to upload:"
for (( i=0; i<${#selected_local[@]}; i++ )); do
  echo "  ${LOCAL_BASE}/${selected_local[$i]}"
  echo "    -> ${REMOTE_BASE}/${selected_remote[$i]}"
done
echo
echo "This will overwrite the remote files listed above."
echo "Proceeding without an extra y/N prompt."

batch_file="$(mktemp "${TMPDIR:-/tmp}/sftp_push_recoiljets.XXXXXX")"
cleanup() {
  rm -f "$batch_file"
}
trap cleanup EXIT

{
  printf 'lcd %s\n' "$LOCAL_BASE"
  printf 'cd %s\n' "$REMOTE_BASE"
  for (( i=0; i<${#selected_local[@]}; i++ )); do
    printf 'put %s %s\n' "$(upload_source_path "${selected_local[$i]}")" "${selected_remote[$i]}"
  done
} > "$batch_file"

echo
echo "sftp batch commands:"
sed 's/^/  /' "$batch_file"
echo
echo "Opening interactive sftp. Enter your SDCC password when prompted."
if [[ "$REMOTE_TRANSPORT" == "ssh-tar" ]]; then
  if run_ssh_tar_upload; then
    :
  else
    status=$?
    echo "[ERROR] ssh-tar upload failed with exit code ${status}." >&2
    exit "$status"
  fi
elif sftp \
    -oBatchMode=no \
    -oPreferredAuthentications=password,keyboard-interactive,publickey \
    -b "$batch_file" \
    "$REMOTE_HOST"; then
  echo
  echo "[OK] Upload complete."
else
  status=$?
  echo
  echo "[ERROR] sftp upload failed with exit code ${status}." >&2
  echo "[ERROR] No success confirmation was received from sftp." >&2
  exit "$status"
fi

echo
echo "Uploaded ${#selected_local[@]} file(s)."

if (( commit_push )); then
  run_commit_push
fi
