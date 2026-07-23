#include "../../../src/RJPhotonTrainingViewV1.h"

#include <TDirectory.h>
#include <TFile.h>
#include <TKey.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <cstdint>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
using namespace RJReplayFoundationV1;

void require(bool state,const std::string& message)
{
  if(!state)throw std::runtime_error(message);
}

SourceOccurrenceRow source(const std::string& lane,
                           const std::string& sample="training")
{
  const std::string hash(64,'a');
  SourceOccurrenceRow row;
  row.lane=lane;row.dataset="synthetic";row.sample=sample;row.period="test";
  row.si_di_role="SI";row.ownership_state="OWNED";row.run=7;row.segment=11;
  row.input_uri_hash=hash;row.input_file_sha256=hash;row.source_manifest_sha256=hash;
  row.id=RJPhotonTrainingViewV1::sourceOccurrenceIdentity(row);
  return row;
}

void rejectSourceSampleCodeMutation(const std::string& path)
{
  const std::string hash(64,'f');
  const auto sourceRow=source("pp_sim","run28_photonjet20");
  RJPhotonTrainingViewV1::Runtime runtime;std::string error;
  require(runtime.initialize(path,RJPhotonTrainingViewV1::kSystemPP,
                             sourceRow,hash,hash,&error),error);
  require(runtime.sourceSampleCode()==20,
          "numeric source sample code was not derived from the frozen source identity");
  RJPhotonTrainingViewV1::Label label;
  label.cluster_index=0;label.training_label=1;label.is_signal=1;
  label.label_authority="PPG12_SOURCE_ROLE";label.source_role=1;
  label.source_sample_code=3;label.ppg12_source_role_label=1;
  label.weight_application_count=1;
  require(!runtime.recordLabel(1,0,label,&error),
          "legacy sample-bin ordinal was accepted as the source threshold");
  require(error=="training label sample code does not match the frozen source identity",
          "source-sample mutation failed for an unexpected reason: "+error);
}

std::vector<ShowerFeatureViewRow> views(const Identity128& candidateId)
{
  std::vector<ShowerFeatureViewRow> result;
  int index=0;
  for(const auto& definition:RJShowerFactorialV1::definitions())
  {
    ShowerFeatureViewRow row;
    row.candidate_id=candidateId;
    row.definition_id=makeIdentity(std::string("shower-definition|")+definition.name+"|"+
                                   RJShowerFactorialV1::semanticText(definition));
    row.definition_name=definition.name;
    row.semantic_sha256=RJShowerFactorialV1::semanticSha256(definition);
    row.weta_cogx=0.01+0.001*index;
    row.wphi_cogx=0.02+0.001*index;
    row.weta33_cogx=0.03+0.001*index;
    row.wphi33_cogx=0.04+0.001*index;
    row.e11_over_e33=0.8;
    row.native_et1=1.1;row.native_et2=1.2;row.native_et3=1.3;row.native_et4=1.4;
    row.e32_over_e35=0.9;
    row.finite_feature_state=1;
    result.push_back(row);
    ++index;
  }
  return result;
}

void validateDefaultReplayInventory(const std::string& path)
{
  const std::string hash(64,'b');
  {
    TFile output(path.c_str(),"RECREATE");
    Writer writer;std::string error;
    require(writer.initialize(&output,Metadata{hash,hash,hash,hash,hash,hash},&error),error);
    require(writer.finish(&error),error);
    output.Write("",TObject::kOverwrite);
  }
  TFile input(path.c_str(),"READ");
  require(input.IsOpen()&&!input.IsZombie(),"default replay artifact is unreadable");
  auto* directory=input.GetDirectory("ReplayFoundationV1");
  require(directory!=nullptr,"ReplayFoundationV1 directory is missing");
  int treeCount=0;
  TIter next(directory->GetListOfKeys());
  while(auto* key=dynamic_cast<TKey*>(next()))
    if(std::string(key->GetClassName())=="TTree")++treeCount;
  require(treeCount==16,"default ReplayFoundationV1 inventory is not exactly 16 trees");
  require(input.Get(RJPhotonTrainingViewV1::kTreeName)==nullptr,
          "opt-in training tree leaked into the default replay artifact");
}

void writeTrainingArtifact(const std::string& path,int systemCode,double et,double centrality,int labelValue)
{
  const std::string hash(64,'c');
  const auto sourceRow=source(systemCode==RJPhotonTrainingViewV1::kSystemPP?"pp_sim":"auau_embed");
  const auto eventId=RJPhotonTrainingViewV1::eventIdentity(sourceRow,sourceRow.run,3);
  constexpr std::uint64_t clusterMapKey=17;
  const auto candidateId=makeIdentity(eventId.hex()+"|candidate|0|"+std::to_string(clusterMapKey));
  RJPhotonTrainingViewV1::Runtime runtime;std::string error;
  require(runtime.initialize(path,systemCode,sourceRow,hash,hash,&error),error);
  RJPhotonTrainingViewV1::Label label;
  label.cluster_index=0;label.training_label=labelValue;label.is_signal=labelValue==1?1:0;
  label.label_authority="PPG12_SOURCE_ROLE";label.source_role=labelValue==1?1:2;
  label.ppg12_source_role_label=labelValue;label.weight_application_count=1;
  if(systemCode==RJPhotonTrainingViewV1::kSystemPP)
  {
    label.weight_slice=2.0;label.weight_cross_section=2.0;
    label.weight_vertex=0.5;label.weight_si_di=0.8;
    label.weight_period=1.25;label.weight_exposure=1.25;
    label.weight_final=1.0;
  }
  else
  {
    label.weight_vertex=1.25;label.weight_exposure=0.8;
    label.weight_final=1.0;
  }
  require(runtime.recordLabel(3,0,label,&error),error);
  RJPhotonTrainingViewV1::CandidateContext candidate;
  candidate.event_id=eventId;candidate.candidate_id=candidateId;candidate.run=7;
  candidate.event_sequence=3;candidate.encounter_ordinal=0;candidate.cluster_et=et;
  candidate.cluster_map_key=clusterMapKey;
  candidate.eta=0.1;candidate.phi=0.2;candidate.vertex_z=2.5;
  candidate.centrality=centrality;candidate.event_weight=1.0;
  auto candidateViews=views(candidateId);
  if(systemCode==RJPhotonTrainingViewV1::kSystemPP)
  {
    for(auto& view:candidateViews)
    {
      // Exercise the system-specific finiteness boundary: these fields are
      // Au+Au-only diagnostics and are not members of the canonical p+p 11.
      view.weta33_cogx=std::numeric_limits<double>::quiet_NaN();
      view.wphi33_cogx=std::numeric_limits<double>::quiet_NaN();
      view.finite_feature_state=0;
    }
  }
  require(runtime.appendCandidate(candidate,candidateViews,&error),error);
  require(runtime.finishEvent(3,&error),error);
  require(runtime.entries()==7,"training runtime did not buffer exactly seven rows");
  require(runtime.finish(&error),error);
}

void rejectCandidateIdentityMutation(const std::string& path)
{
  const std::string hash(64,'d');
  const auto sourceRow=source("pp_sim");
  const auto eventId=RJPhotonTrainingViewV1::eventIdentity(sourceRow,sourceRow.run,9);
  const auto wrongCandidateId=makeIdentity(eventId.hex()+"|candidate|mutated");
  RJPhotonTrainingViewV1::Runtime runtime;std::string error;
  require(runtime.initialize(path,RJPhotonTrainingViewV1::kSystemPP,
                             sourceRow,hash,hash,&error),error);
  RJPhotonTrainingViewV1::Label label;
  label.cluster_index=0;label.training_label=1;label.is_signal=1;
  label.label_authority="PPG12_SOURCE_ROLE";label.source_role=1;
  label.ppg12_source_role_label=1;label.weight_application_count=1;
  require(runtime.recordLabel(9,0,label,&error),error);
  RJPhotonTrainingViewV1::CandidateContext candidate;
  candidate.event_id=eventId;candidate.candidate_id=wrongCandidateId;
  candidate.run=sourceRow.run;candidate.event_sequence=9;candidate.encounter_ordinal=0;
  candidate.cluster_map_key=17;candidate.cluster_et=20.0;
  candidate.eta=0.1;candidate.phi=0.2;candidate.vertex_z=2.5;
  require(!runtime.appendCandidate(candidate,views(wrongCandidateId),&error),
          "candidate-ID mutation was not rejected");
  require(error=="training candidate stable identity input mismatch",
          "candidate-ID mutation failed for an unexpected reason: "+error);
}

void rejectSourceAndEventIdentityMutations(const std::string& sourcePath,
                                           const std::string& eventPath)
{
  const std::string hash(64,'e');
  auto mutatedSource=source("pp_sim");
  mutatedSource.id=makeIdentity("mutated-source-identity");
  RJPhotonTrainingViewV1::Runtime sourceRuntime;std::string error;
  require(!sourceRuntime.initialize(sourcePath,RJPhotonTrainingViewV1::kSystemPP,
                                    mutatedSource,hash,hash,&error),
          "source-ID mutation was not rejected");
  require(error=="training-view source stable identity or input/manifest hash is invalid",
          "source-ID mutation failed for an unexpected reason: "+error);

  const auto sourceRow=source("pp_sim");
  RJPhotonTrainingViewV1::Runtime eventRuntime;error.clear();
  require(eventRuntime.initialize(eventPath,RJPhotonTrainingViewV1::kSystemPP,
                                  sourceRow,hash,hash,&error),error);
  RJPhotonTrainingViewV1::Label label;
  label.cluster_index=0;label.training_label=1;label.is_signal=1;
  label.label_authority="PPG12_SOURCE_ROLE";label.source_role=1;
  label.ppg12_source_role_label=1;label.weight_application_count=1;
  require(eventRuntime.recordLabel(12,0,label,&error),error);
  RJPhotonTrainingViewV1::CandidateContext candidate;
  candidate.event_id=makeIdentity("mutated-event-identity");
  candidate.run=sourceRow.run;candidate.event_sequence=12;
  candidate.encounter_ordinal=0;candidate.cluster_map_key=17;
  candidate.candidate_id=makeIdentity(
      candidate.event_id.hex()+"|candidate|0|"+std::to_string(candidate.cluster_map_key));
  candidate.cluster_et=20.0;candidate.eta=0.1;candidate.phi=0.2;
  candidate.vertex_z=2.5;
  require(!eventRuntime.appendCandidate(candidate,views(candidate.candidate_id),&error),
          "event-ID mutation was not rejected");
  require(error=="training candidate event stable identity input mismatch",
          "event-ID mutation failed for an unexpected reason: "+error);
}

void rejectWeightLedgerMutations(const std::string& aliasPath,
                                 const std::string& productPath,
                                 const std::string& eventWeightPath,
                                 const std::string& ordinalPath)
{
  const std::string hash(64,'f');
  const auto sourceRow=source("pp_sim");
  {
    RJPhotonTrainingViewV1::Runtime runtime;std::string error;
    require(runtime.initialize(aliasPath,RJPhotonTrainingViewV1::kSystemPP,
                               sourceRow,hash,hash,&error),error);
    RJPhotonTrainingViewV1::Label label;
    label.cluster_index=0;label.training_label=1;label.is_signal=1;
    label.label_authority="PPG12_SOURCE_ROLE";label.source_role=1;
    label.ppg12_source_role_label=1;label.weight_application_count=1;
    label.weight_slice=2.0;label.weight_cross_section=1.0;
    label.weight_period=1.25;label.weight_exposure=1.25;
    label.weight_final=2.5;
    require(!runtime.recordLabel(1,0,label,&error),
            "p+p weight-alias mutation was not rejected");
    require(error=="p+p training weight provenance aliases do not close",
            "p+p weight-alias mutation failed for an unexpected reason: "+error);
  }
  {
    RJPhotonTrainingViewV1::Runtime runtime;std::string error;
    require(runtime.initialize(productPath,RJPhotonTrainingViewV1::kSystemPP,
                               sourceRow,hash,hash,&error),error);
    RJPhotonTrainingViewV1::Label label;
    label.cluster_index=0;label.training_label=1;label.is_signal=1;
    label.label_authority="PPG12_SOURCE_ROLE";label.source_role=1;
    label.ppg12_source_role_label=1;label.weight_application_count=1;
    label.weight_slice=2.0;label.weight_cross_section=2.0;
    label.weight_vertex=0.5;label.weight_si_di=0.8;
    label.weight_period=1.25;label.weight_exposure=1.25;
    label.weight_final=2.0;
    require(!runtime.recordLabel(1,0,label,&error),
            "p+p weight-product mutation was not rejected");
    require(error=="training label weight components do not close to final weight",
            "p+p weight-product mutation failed for an unexpected reason: "+error);
  }
  {
    RJPhotonTrainingViewV1::Runtime runtime;std::string error;
    require(runtime.initialize(eventWeightPath,RJPhotonTrainingViewV1::kSystemPP,
                               sourceRow,hash,hash,&error),error);
    RJPhotonTrainingViewV1::Label label;
    label.cluster_index=0;label.training_label=1;label.is_signal=1;
    label.label_authority="PPG12_SOURCE_ROLE";label.source_role=1;
    label.ppg12_source_role_label=1;label.weight_application_count=1;
    label.weight_slice=2.0;label.weight_cross_section=2.0;
    label.weight_vertex=0.5;label.weight_si_di=0.8;
    label.weight_period=1.25;label.weight_exposure=1.25;
    label.weight_final=1.0;
    require(runtime.recordLabel(1,0,label,&error),error);
    const auto eventId=RJPhotonTrainingViewV1::eventIdentity(sourceRow,sourceRow.run,1);
    constexpr std::uint64_t clusterMapKey=17;
    RJPhotonTrainingViewV1::CandidateContext candidate;
    candidate.event_id=eventId;
    candidate.candidate_id=makeIdentity(
        eventId.hex()+"|candidate|0|"+std::to_string(clusterMapKey));
    candidate.run=sourceRow.run;candidate.event_sequence=1;
    candidate.encounter_ordinal=0;candidate.cluster_map_key=clusterMapKey;
    candidate.cluster_et=20.0;candidate.eta=0.1;candidate.phi=0.2;
    candidate.vertex_z=2.5;candidate.centrality=-1.0;
    candidate.event_weight=1.25;
    require(!runtime.appendCandidate(candidate,views(candidate.candidate_id),&error),
            "candidate/final event-weight mutation was not rejected");
    require(error=="training candidate event weight does not match final weight ledger",
            "candidate/final event-weight mutation failed for an unexpected reason: "+error);
  }
  {
    RJPhotonTrainingViewV1::Runtime runtime;std::string error;
    require(runtime.initialize(ordinalPath,RJPhotonTrainingViewV1::kSystemPP,
                               sourceRow,hash,hash,&error),error);
    RJPhotonTrainingViewV1::Label label;
    label.cluster_index=5;label.training_label=1;label.is_signal=1;
    label.label_authority="PPG12_SOURCE_ROLE";label.source_role=1;
    label.ppg12_source_role_label=1;label.weight_application_count=1;
    require(!runtime.recordLabel(1,0,label,&error),
            "label cluster-index/ordinal mutation was not rejected");
    require(error=="training-label cluster index does not match encounter ordinal",
            "label cluster-index/ordinal mutation failed for an unexpected reason: "+error);
  }
}

void validateTrainingArtifact(const std::string& path,int expectedFeatureCount,
                              int expectedDomain,int expectedBelow15,
                              int expectedNominal,double expectedCentrality,
                              int expectedSystem,int expectedLabel)
{
  TFile input(path.c_str(),"READ");
  require(input.IsOpen()&&!input.IsZombie(),"training artifact is unreadable");
  auto* tree=dynamic_cast<TTree*>(input.Get(RJPhotonTrainingViewV1::kTreeName));
  require(tree!=nullptr,"training tree is missing");
  require(tree->GetEntries()==7,"training tree does not have seven definition rows");
  require(input.GetDirectory("ReplayFoundationV1")==nullptr,
          "separate training artifact unexpectedly owns a replay directory");

  std::string* definition=nullptr;
  std::string* showerSemantic=nullptr;
  std::string* featureContract=nullptr;
  std::vector<float>* features=nullptr;
  int featureCount=0,finite=0,system=0,domain=0,below15=0,nominal=0,wp=0,tag=0;
  int trainingLabel=0,weightApplicationCount=0;
  std::uint64_t clusterMapKey=0;
  double centrality=0.0,eventWeight=0.0,weightSlice=0.0,weightCrossSection=0.0;
  double weightVertex=0.0,weightSiDi=0.0,weightPeriod=0.0;
  double weightExposure=0.0,weightFinal=0.0;
  std::uint64_t trainingHi=0,trainingLo=0,sourceHi=0,sourceLo=0;
  std::uint64_t eventHi=0,eventLo=0,candidateHi=0,candidateLo=0,definitionHi=0,definitionLo=0;
  tree->SetBranchAddress("definition_name",&definition);
  tree->SetBranchAddress("shower_semantic_sha256",&showerSemantic);
  tree->SetBranchAddress("feature_contract_sha256",&featureContract);
  tree->SetBranchAddress("ordered_features",&features);
  tree->SetBranchAddress("feature_count",&featureCount);
  tree->SetBranchAddress("finite_feature_state",&finite);
  tree->SetBranchAddress("system_code",&system);
  tree->SetBranchAddress("model_domain_state",&domain);
  tree->SetBranchAddress("below15_retention_state",&below15);
  tree->SetBranchAddress("nominal_training_eligible",&nominal);
  tree->SetBranchAddress("working_point_state",&wp);
  tree->SetBranchAddress("tag_state",&tag);
  tree->SetBranchAddress("training_label",&trainingLabel);
  tree->SetBranchAddress("weight_application_count",&weightApplicationCount);
  tree->SetBranchAddress("weight_slice",&weightSlice);
  tree->SetBranchAddress("weight_cross_section",&weightCrossSection);
  tree->SetBranchAddress("weight_vertex",&weightVertex);
  tree->SetBranchAddress("weight_si_di",&weightSiDi);
  tree->SetBranchAddress("weight_period",&weightPeriod);
  tree->SetBranchAddress("weight_exposure",&weightExposure);
  tree->SetBranchAddress("weight_final",&weightFinal);
  tree->SetBranchAddress("cluster_map_key",&clusterMapKey);
  tree->SetBranchAddress("centrality",&centrality);
  tree->SetBranchAddress("event_weight",&eventWeight);
  tree->SetBranchAddress("training_view_id_hi",&trainingHi);
  tree->SetBranchAddress("training_view_id_lo",&trainingLo);
  tree->SetBranchAddress("source_occurrence_id_hi",&sourceHi);
  tree->SetBranchAddress("source_occurrence_id_lo",&sourceLo);
  tree->SetBranchAddress("event_id_hi",&eventHi);
  tree->SetBranchAddress("event_id_lo",&eventLo);
  tree->SetBranchAddress("candidate_id_hi",&candidateHi);
  tree->SetBranchAddress("candidate_id_lo",&candidateLo);
  tree->SetBranchAddress("definition_id_hi",&definitionHi);
  tree->SetBranchAddress("definition_id_lo",&definitionLo);
  std::set<std::string> definitions;
  std::set<std::pair<std::uint64_t,std::uint64_t>> trainingIds;
  std::pair<std::uint64_t,std::uint64_t> expectedSource{0,0};
  std::pair<std::uint64_t,std::uint64_t> expectedEvent{0,0};
  std::pair<std::uint64_t,std::uint64_t> expectedCandidate{0,0};
  for(Long64_t entry=0;entry<tree->GetEntries();++entry)
  {
    tree->GetEntry(entry);
    require(definition!=nullptr&&definitions.insert(*definition).second,
            "definition identity is missing or duplicated");
    require(features!=nullptr&&static_cast<int>(features->size())==expectedFeatureCount,
            "ordered feature vector has the wrong length");
    require(featureCount==expectedFeatureCount&&finite==1,
            "feature count/finite state is invalid");
    require(system==expectedSystem&&trainingLabel==expectedLabel&&weightApplicationCount==1,
            "system, source-role label, or weight application count changed");
    require(std::fabs(eventWeight-weightFinal)<1.0e-12,
            "candidate event weight and final weight ledger disagree");
    if(expectedSystem==RJPhotonTrainingViewV1::kSystemPP)
    {
      require(weightSlice==2.0&&weightCrossSection==2.0&&weightVertex==0.5&&
              weightSiDi==0.8&&weightPeriod==1.25&&weightExposure==1.25&&
              weightFinal==1.0,
              "non-unit p+p weight aliases or unique-factor product changed");
    }
    else
    {
      require(weightSlice==1.0&&weightCrossSection==1.0&&weightVertex==1.25&&
              weightSiDi==1.0&&weightPeriod==1.0&&weightExposure==0.8&&
              weightFinal==1.0,
              "non-unit AuAu vertex/exposure product changed");
    }
    require(featureContract!=nullptr&&
            *featureContract==RJPhotonTrainingViewV1::featureContractSha256(expectedSystem),
            "feature contract hash changed");
    require(showerSemantic!=nullptr&&definition!=nullptr&&
            *showerSemantic==RJShowerFactorialV1::semanticSha256(
                RJShowerFactorialV1::definition(*definition)),
            "shower semantic hash changed");
    require(domain==expectedDomain&&below15==expectedBelow15&&nominal==expectedNominal,
            "model-domain or nominal-training state is invalid");
    require(wp==-1&&tag==-1,"training artifact owns a forbidden WP/tag state");
    require((trainingHi|trainingLo)!=0&&(sourceHi|sourceLo)!=0&&(eventHi|eventLo)!=0&&
            (candidateHi|candidateLo)!=0&&(definitionHi|definitionLo)!=0,
            "source/event/candidate/definition/training identity is null");
    const auto rebuiltCandidate=makeIdentity(
        Identity128{eventHi,eventLo}.hex()+"|candidate|0|"+std::to_string(clusterMapKey));
    require(rebuiltCandidate==Identity128{candidateHi,candidateLo},
            "candidate identity cannot be rebuilt from its stable inputs");
    require(trainingIds.emplace(trainingHi,trainingLo).second,
            "training-view identity is duplicated");
    const std::pair<std::uint64_t,std::uint64_t> sourceId{sourceHi,sourceLo};
    const std::pair<std::uint64_t,std::uint64_t> eventId{eventHi,eventLo};
    const std::pair<std::uint64_t,std::uint64_t> candidateId{candidateHi,candidateLo};
    if(entry==0){expectedSource=sourceId;expectedEvent=eventId;expectedCandidate=candidateId;}
    require(sourceId==expectedSource&&eventId==expectedEvent&&candidateId==expectedCandidate,
            "factorial rows do not share exact source/event/candidate identities");
    require(std::fabs(centrality-expectedCentrality)<1.0e-12,
            "centrality provenance changed");
    if(expectedFeatureCount==14)
    {
      require(std::fabs(features->at(3)-static_cast<float>(0.03+0.001*entry))<1.0e-6F,
              "AuAu weta33 is not feature index 3");
      require(std::fabs(features->at(13)-static_cast<float>(expectedCentrality))<1.0e-6F,
              "AuAu centrality is not the final feature");
    }
  }
  require(definitions==std::set<std::string>({"H70","H0","G70","G0","O70","O0","R70"}),
          "training artifact definition set is incomplete");
}
} // namespace

int main(int argc,char** argv)
{
  if(argc!=4)
  {
    std::cerr<<"usage: the134_multiview_training_canary DEFAULT.root PP.root AUAU.root\n";
    return 64;
  }
  try
  {
    validateDefaultReplayInventory(argv[1]);
    writeTrainingArtifact(argv[2],RJPhotonTrainingViewV1::kSystemPP,20.0,-1.0,1);
    validateTrainingArtifact(argv[2],11,static_cast<int>(ModelApplicability::VALIDATED_DOMAIN),0,1,-1.0,RJPhotonTrainingViewV1::kSystemPP,1);
    writeTrainingArtifact(argv[3],RJPhotonTrainingViewV1::kSystemAuAu,12.0,30.25,-1);
    validateTrainingArtifact(argv[3],14,static_cast<int>(ModelApplicability::DIAGNOSTIC_EXTRAPOLATION),1,0,30.25,RJPhotonTrainingViewV1::kSystemAuAu,-1);
    const std::string mutationPath=std::string(argv[2])+".candidate_identity_mutation.root";
    rejectCandidateIdentityMutation(mutationPath);
    std::remove(mutationPath.c_str());
    const std::string sourceMutationPath=std::string(argv[2])+".source_identity_mutation.root";
    const std::string eventMutationPath=std::string(argv[2])+".event_identity_mutation.root";
    rejectSourceAndEventIdentityMutations(sourceMutationPath,eventMutationPath);
    std::remove(sourceMutationPath.c_str());
    std::remove(eventMutationPath.c_str());
    const std::string sampleCodeMutationPath=
        std::string(argv[2])+".source_sample_code_mutation.root";
    rejectSourceSampleCodeMutation(sampleCodeMutationPath);
    std::remove(sampleCodeMutationPath.c_str());
    const std::string aliasMutationPath=std::string(argv[2])+".weight_alias_mutation.root";
    const std::string productMutationPath=std::string(argv[2])+".weight_product_mutation.root";
    const std::string eventWeightMutationPath=std::string(argv[2])+".event_weight_mutation.root";
    const std::string ordinalMutationPath=std::string(argv[2])+".label_ordinal_mutation.root";
    rejectWeightLedgerMutations(aliasMutationPath,productMutationPath,
                                eventWeightMutationPath,ordinalMutationPath);
    std::remove(aliasMutationPath.c_str());
    std::remove(productMutationPath.c_str());
    std::remove(eventWeightMutationPath.c_str());
    std::remove(ordinalMutationPath.c_str());
    std::cout<<"THE134_MULTIVIEW_TRAINING_CANARY_PASS\n";
    return 0;
  }
  catch(const std::exception& error)
  {
    std::cerr<<"THE134_MULTIVIEW_TRAINING_CANARY_FAIL: "<<error.what()<<'\n';
    return 1;
  }
}
