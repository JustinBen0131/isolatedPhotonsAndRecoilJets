#include "../../../src/RJReplayFoundationV1.h"

#include <TFile.h>

#include <iostream>
#include <string>

using namespace RJReplayFoundationV1;

int main(int argc, char** argv)
{
  if (argc != 2)
  {
    std::cerr << "usage: the118_schema_canary OUTPUT.root\n";
    return 64;
  }

  TFile output(argv[1], "RECREATE", "THE-118 normalized replay schema canary");
  const std::string h(64, 'a');
  Writer writer;
  std::string error;
  if (!writer.initialize(&output, Metadata{h,h,h,h,h,h}, &error))
  {
    std::cerr << error << '\n';
    return 1;
  }

  const auto sourceId = makeIdentity("source|pp_data|run1|seg2");
  const auto eventId = makeIdentity(sourceId.hex() + "|event|3");
  const auto candidateId = makeIdentity(eventId.hex() + "|candidate|0");
  const auto modelId = makeIdentity("model|the116|sha");
  const auto jetId = makeIdentity(eventId.hex() + "|jet|r04|0");
  const auto pairId = makeIdentity(candidateId.hex() + "|" + jetId.hex());
  const auto truthPhotonId = makeIdentity(eventId.hex() + "|truth_photon|0");
  const auto truthJetId = makeIdentity(eventId.hex() + "|truth_jet|r04|0");

  SourceOccurrenceRow source;
  source.id=sourceId; source.lane="pp_data"; source.dataset="synthetic"; source.sample="data";
  source.period="0mrad"; source.run=1; source.segment=2; source.si_di_role="DATA";
  source.ownership_state="OWNED"; source.input_uri_hash=h; source.input_file_sha256=h; source.source_manifest_sha256=h;
  if(!writer.fill(source,&error)) return 2;

  EventRow event;
  event.id=eventId; event.source_id=sourceId; event.run=1; event.event_sequence=3; event.trigger_bits=4;
  event.vertex_z=0.25; event.centrality=-1; event.event_weight=1; event.terminal_status=0;
  event.candidate_count=1; event.tag_count=1; event.recoil_count=1;
  if(!writer.fill(event,&error)) return 3;

  PhotonCandidateRow candidate;
  candidate.id=candidateId; candidate.event_id=eventId; candidate.encounter_ordinal=0;
  candidate.rank_keys={18.0,0.1,0.2}; candidate.cluster_et=18; candidate.eta=0.1; candidate.phi=0.2;
  candidate.ordered_features=std::vector<float>(11,0.5F); candidate.preselection_bitmask=7;
  candidate.shower_definition_views={"H70"};
  candidate.finite_feature_state=1; candidate.below15_retention_state=0;
  if(!writer.fill(candidate,&error)) return 4;

  ModelEvaluationRow model;
  model.candidate_id=candidateId; model.model_id=modelId; model.model_sha256=h;
  model.shower_definition_id="H70"; model.shower_semantic_sha256=h;
  model.ordered_input_witnesses=candidate.ordered_features; model.raw_score=0.7; model.finite_score=1;
  model.applicability_state=static_cast<int>(ModelApplicability::VALIDATED_DOMAIN);
  model.wp70=0.8; model.wp80=0.6; model.wp90=0.4; model.delta_wp70=-0.1; model.delta_wp80=0.1; model.delta_wp90=0.3;
  if(!writer.fill(model,&error)) return 5;

  ShowerCellRow cell;
  cell.candidate_id=candidateId; cell.local_eta_index=0; cell.local_phi_index=0;
  cell.tower_eta_index=48; cell.tower_phi_index=128; cell.tower_key=11;
  cell.calibrated_energy=2.0; cell.rawcluster_map_value=1.9; cell.is_good=1; cell.seed_state=1;
  cell.denominator_membership=1; cell.rawcluster_owned=1; cell.rawcluster_value_present=1;
  cell.floor0_membership=1; cell.floor70_membership=1; cell.grid_membership_bitmask=3;
  if(!writer.fill(cell,&error)) return 6;

  ShowerFeatureViewRow view;
  view.candidate_id=candidateId; view.definition_id=makeIdentity("shower-definition|H70");
  view.definition_name="H70"; view.semantic_sha256=h; view.ordered_features=candidate.ordered_features;
  view.floor_gev=0.070; view.cog_eta=3.25; view.cog_phi=2.75;
  view.raw_center_eta=48.75; view.raw_center_phi=128.25;
  view.center_eta_index=48; view.center_phi_index=128;
  view.e11=2.0; view.e33=2.0; view.e32=2.0; view.e35=2.0;
  view.e11_over_e33=1.0; view.e32_over_e35=1.0;
  view.weta_cogx=0.1; view.wphi_cogx=0.2; view.weta33_cogx=0.1; view.wphi33_cogx=0.2;
  view.native_et1=2.0; view.native_et2=2.0; view.native_et3=2.0; view.native_et4=2.0;
  view.moment_eta_numerator=0.2; view.moment_phi_numerator=0.4; view.moment_denominator=2.0;
  view.moment33_eta_numerator=0.2; view.moment33_phi_numerator=0.4; view.moment33_denominator=2.0;
  view.energy_source=0; view.rectangular_membership=0; view.moment_membership=1;
  view.finite_feature_state=1; view.good_cell_count=1; view.owned_cell_count=1;
  view.active_sum_cell_count=1; view.active_moment_cell_count=1;
  if(!writer.fill(view,&error)) return 25;

  IsolationConstituentRow iso;
  iso.candidate_id=candidateId; iso.constituent_id=makeIdentity(candidateId.hex()+"|iso|0");
  iso.delta_eta=0.1; iso.delta_phi=0.1; iso.delta_r=0.1414; iso.subsystem=1;
  iso.raw_energy=0.5; iso.calibrated_energy=0.51; iso.sub1_energy=0.3; iso.phosub_residual=0.02;
  iso.quality_state=1; iso.mask_state=0; iso.candidate_removal_state=0;
  if(!writer.fill(iso,&error)) return 7;

  IsolationWitnessRow witness;
  witness.candidate_id=candidateId; witness.isolation_id=makeIdentity("standard_sub1|r04|reco");
  witness.radius=0.4; witness.subtraction_method=1; witness.reconstructed_or_truth=0;
  witness.cone_sum=1.2; witness.threshold=2.0; witness.pass_state=1; witness.constituent_count=1;
  if(!writer.fill(witness,&error)) return 8;

  JetRow jet;
  jet.id=jetId; jet.event_id=eventId; jet.algorithm="antikt"; jet.radius=0.4;
  jet.input_identity="tower"; jet.subtraction_identity="sub1"; jet.raw_pt=11; jet.corrected_pt=10;
  jet.eta=-0.2; jet.phi=3.0; jet.deterministic_order=0;
  if(!writer.fill(jet,&error)) return 9;

  JetConstituentRow jc;
  jc.jet_id=jetId; jc.constituent_id=makeIdentity(jetId.hex()+"|constituent|0"); jc.constituent_ordinal=0;
  jc.subsystem=1; jc.energy=3; jc.eta=-0.2; jc.phi=3.0; jc.quality_state=1;
  if(!writer.fill(jc,&error)) return 10;

  PhotonJetPairRow pair;
  pair.id=pairId; pair.event_id=eventId; pair.candidate_id=candidateId; pair.jet_id=jetId;
  pair.delta_phi=2.8; pair.xjgamma=10.0/18.0; pair.recoil_state=1;
  if(!writer.fill(pair,&error)) return 11;

  TruthPhotonRow truthPhoton;
  truthPhoton.id=truthPhotonId; truthPhoton.event_id=eventId; truthPhoton.pt=18.2;
  truthPhoton.eta=0.1; truthPhoton.phi=0.2; truthPhoton.prompt_class=1; truthPhoton.source_role=1;
  truthPhoton.truth_isolation_witness=0.2; truthPhoton.reporting_guard_state=0;
  if(!writer.fill(truthPhoton,&error)) return 12;

  TruthJetRow truthJet;
  truthJet.id=truthJetId; truthJet.event_id=eventId; truthJet.algorithm="antikt"; truthJet.radius=0.4;
  truthJet.pt=10.5; truthJet.eta=-0.2; truthJet.phi=3.0; truthJet.ownership_state="OWNED";
  if(!writer.fill(truthJet,&error)) return 13;

  RecoTruthLinkRow photonLink;
  photonLink.id=makeIdentity(candidateId.hex()+"|"+truthPhotonId.hex());
  photonLink.reco_type=static_cast<int>(RecoTruthType::PHOTON); photonLink.reco_id=candidateId;
  photonLink.truth_type=static_cast<int>(RecoTruthType::PHOTON); photonLink.truth_id=truthPhotonId;
  photonLink.match_metric=0.01; photonLink.link_class=static_cast<int>(LinkClass::MATCH);
  if(!writer.fill(photonLink,&error)) return 14;

  const auto jetLinks=buildDeterministicJetLinks(
      {{jetId,jet.corrected_pt,jet.eta,jet.phi,jet.radius}},
      {{truthJetId,truthJet.pt,truthJet.eta,truthJet.phi,truthJet.radius}},0.3);
  if(jetLinks.size()!=2 ||
     jetLinks[0].link_class!=static_cast<int>(LinkClass::MATCH_CANDIDATE) ||
     jetLinks[1].link_class!=static_cast<int>(LinkClass::MATCH)) return 15;
  for(const auto& link:jetLinks) if(!writer.fill(link,&error)) return 15;

  // Deterministic one-to-many/many-to-one stress: the lower-dR edge wins;
  // the remaining reco and truth objects form explicit fake/miss rows while
  // every viable edge remains retained as MATCH_CANDIDATE.
  const auto r0=makeIdentity("jet-match-r0"),r1=makeIdentity("jet-match-r1");
  const auto t0=makeIdentity("jet-match-t0"),t1=makeIdentity("jet-match-t1");
  const auto stress=buildDeterministicJetLinks(
      {{r0,12.0,0.00,0.00,0.4},{r1,20.0,0.02,0.00,0.4}},
      {{t0,11.0,0.01,0.00,0.4},{t1,10.0,1.00,1.00,0.4}},0.3);
  int candidates=0,matches=0,fakes=0,misses=0;
  for(const auto& link:stress)
  {
    candidates+=link.link_class==static_cast<int>(LinkClass::MATCH_CANDIDATE);
    matches+=link.link_class==static_cast<int>(LinkClass::MATCH);
    fakes+=link.link_class==static_cast<int>(LinkClass::RECO_FAKE);
    misses+=link.link_class==static_cast<int>(LinkClass::TRUTH_MISS);
  }
  if(candidates!=2||matches!=1||fakes!=1||misses!=1) return 18;

  WeightComponentRow weight;
  weight.target_id=eventId; weight.component_type="event_final"; weight.final_weight=1; weight.application_count=1;
  if(!writer.fill(weight,&error)) return 19;

  EventDisplaySnapshotRow snapshot;
  snapshot.id=makeIdentity(eventId.hex()+"|snapshot|0"); snapshot.event_id=eventId;
  snapshot.selection_reason="synthetic"; snapshot.quota_class="bounded"; snapshot.serialized_payload_hash=h;
  if(!writer.fill(snapshot,&error)) return 20;

  // Negative fixtures must fail closed without changing the valid population.
  EventRow orphan=event; orphan.id=makeIdentity("orphan-event"); orphan.source_id=makeIdentity("missing-source");
  if(writer.fill(orphan,&error)) return 21;
  if(writer.fill(candidate,&error)) return 22; // duplicate candidate
  if(writer.fill(cell,&error)) return 26; // duplicate candidate/local-cell key
  if(writer.fill(view,&error)) return 27; // duplicate candidate/definition key
  RecoTruthLinkRow invalid=photonLink; invalid.id=makeIdentity("invalid-link"); invalid.truth_id=makeIdentity("missing-truth");
  if(writer.fill(invalid,&error)) return 23;

  if(!writer.finish(&error))
  {
    std::cerr << error << '\n';
    return 24;
  }
  output.Close();
  std::cout << "THE118_SCHEMA_CANARY_PASS\n";
  return 0;
}
