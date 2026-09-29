#ifndef RJ_PRODUCER_CAPTURE_V1_H
#define RJ_PRODUCER_CAPTURE_V1_H

// Additive raw witnesses only. No selection, reconstruction, scoring, unit
// conversion, or matching is performed here. Templates keep focused fixtures
// independent of Fun4All/ROOT and use the existing native object accessors.
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace RJProducerCaptureV1
{
inline float missingFloat() { return std::numeric_limits<float>::quiet_NaN(); }
inline double missingDouble() { return std::numeric_limits<double>::quiet_NaN(); }

inline const std::vector<std::string>& npbFeatureNames()
{
  static const std::vector<std::string> names={
    "cluster_Et","cluster_Eta","vertexz","e11_over_e33","e32_over_e35",
    "e11_over_e22","e11_over_e13","e11_over_e15","e11_over_e17",
    "e11_over_e31","e11_over_e51","e11_over_e71","e22_over_e33",
    "e22_over_e35","e22_over_e37","e22_over_e53","cluster_weta_cogx",
    "cluster_wphi_cogx","cluster_et1","cluster_et2","cluster_et3",
    "cluster_et4","cluster_w32","cluster_w52","cluster_w72"};
  return names;
}

// State records what the existing invocation did, not physics applicability.
enum EvaluationState : std::int32_t
{
  NOT_CAPTURED=0, MODEL_UNAVAILABLE=1, OUTSIDE_PRODUCER_DOMAIN=2,
  INPUT_NONFINITE=3, EVALUATED_FINITE=4, EVALUATED_NONFINITE=5
};

struct NPBWitness
{
  std::int32_t capture_version=0,evaluation_state=NOT_CAPTURED;
  // 1=named score with configured min/max/eta; 2=primary BDT invocation (no
  // additional NPB domain gate). Neither value certifies a validated domain.
  std::int32_t domain_mode=0,in_producer_domain=-1,model_available=-1;
  std::int32_t input_order_is_npb25=0,inputs_captured=0,inputs_finite=0;
  float min_et=missingFloat(),max_et=missingFloat(),max_abs_eta=missingFloat();
  float raw_score=missingFloat();
  std::uint32_t source_photon_id=0;
  float source_cluster_et=missingFloat(),source_eta=missingFloat(),source_phi=missingFloat();
  // Object id and container key are distinct; only an observed container entry
  // establishes the latter. Policy 1=same object, 2=native first-match result.
  std::int32_t binding_policy=0,source_identity_available=0,kinematic_input_bits_available=0;
  std::uint32_t source_native_cluster_key=0,input_et_bits=0,input_eta_bits=0;
  std::int64_t source_encounter_ordinal=-1;
  std::string source_cluster_node;
  std::vector<std::string> feature_names;
  std::vector<float> ordered_inputs;
  std::vector<std::int32_t> ordered_input_finite;
  // Exact configured path only. Resolve/hash through the bound runtime packet;
  // this string is deliberately NOT labelled a verified model hash.
  std::string model_file;
};

struct NativeTiming
{
  std::int32_t capture_version=0,contributing_towers=-1,raw_time_finite=0;
  double mean_time_samples=missingDouble(),energy_denominator=missingDouble();
  double time_energy_numerator=missingDouble(),producer_vertex_z=missingDouble();
  std::string input_cluster_node,calibrated_tower_node;
};

struct MbdTimes
{
  std::int32_t capture_version=0,available=0,native_is_valid=-1;
  std::int32_t t0_finite=0,south_time_finite=0,north_time_finite=0;
  std::int32_t south_npmt=-1,north_npmt=-1;
  double t0_ns=missingDouble(),south_time_ns=missingDouble(),north_time_ns=missingDouble();
};

struct TowerTimeStatus
{
  std::int32_t capture_version=0,available=0,time_finite=0;
  double time_samples=missingDouble();
  std::uint32_t status=0; // full native uint8_t word, not reconstructed from isGood
};

// Same option interpretation as the existing AuAu event-QA path. This does
// not configure or change tower status/calibration or photon selection.
inline bool eventCaloRequireGood(const char* raw)
{
  if(!raw)return true;
  std::string flag(raw);
  std::transform(flag.begin(),flag.end(),flag.begin(),
                 [](unsigned char c){return std::tolower(c);});
  return !(flag=="0"||flag=="false"||flag=="no"||flag=="off");
}

// Observe current-event TOWERINFO_CALIB_{CEMC,HCALIN,HCALOUT}, never SUB1 or
// retowered inputs. Match native AuAu QA: signed finite energies, double
// accumulation, then binary32 storage promoted to the existing double branches.
// PPG12's float-per-addition totals are NOT claimed bitwise equivalent.
// Masks use bits 0/1/2 for CEMC/HCALIN/HCALOUT. Validity is numeric/container
// evidence, not certification of the detector geometry or calibration payload.
// Absent containers are NaN, not fabricated zero. Incomplete present containers
// retain the finite partial sum but have their valid bit cleared. Empty present
// containers also remain invalid. Existing histogram/selection state is untouched.
template<class Row,class Container>
void captureCalorimeterSums(Row& row,const std::array<Container*,3>& containers,bool requireGood)
{
  row.event_calo_capture_version=1;
  row.event_calo_available_mask=0;row.event_calo_valid_mask=0;
  row.event_calo_require_isgood=requireGood?1:0;
  std::array<double,3> sums{{missingDouble(),missingDouble(),missingDouble()}};
  for(std::size_t layer=0;layer<containers.size();++layer)
  {
    auto* towers=containers[layer];
    if(!towers)continue;
    row.event_calo_available_mask|=1U<<layer;
    double sum=0.;bool valid=towers->size()>0;
    for(unsigned int channel=0;channel<towers->size();++channel)
    {
      auto* tower=towers->get_tower_at_channel(channel);
      if(!tower){valid=false;continue;}
      if(requireGood&&!tower->get_isGood())continue;
      const double energy=tower->get_energy();
      if(!std::isfinite(energy)){valid=false;continue;}
      sum+=energy;
    }
    sums[layer]=sum;
    if(valid&&std::isfinite(static_cast<float>(sum)))row.event_calo_valid_mask|=1U<<layer;
  }
  row.emcal_total_energy=static_cast<float>(sums[0]);
  row.ihcal_total_energy=static_cast<float>(sums[1]);
  row.ohcal_total_energy=static_cast<float>(sums[2]);
  row.total_calo_energy=static_cast<float>(sums[0]+sums[1]+sums[2]);
}

template<class Photon> float parameter(const Photon* photon,const std::string& name)
{
  if(!photon)return missingFloat();
  const auto& values=photon->get_all_shower_shapes();
  const auto it=values.find(name);
  return it==values.end()?missingFloat():it->second;
}

template<class Photon> std::string textParameter(const Photon* photon,const std::string& prefix)
{
  if(!photon)return {};
  for(const auto& item:photon->get_all_shower_shapes())
    if(item.first.compare(0,prefix.size(),prefix)==0)return item.first.substr(prefix.size());
  return {};
}

template<class Photon> std::int32_t integerParameter(const Photon* photon,const std::string& name)
{
  const float value=parameter(photon,name);
  return std::isfinite(value)&&value>=-1&&value<=100000&&std::floor(value)==value
      ?static_cast<std::int32_t>(value):-1;
}

inline std::string npbPrefix(const std::string& scoreName)
{ return "rj_capture_"+scoreName+"_v1_"; }

inline std::uint32_t floatBits(float value)
{
  static_assert(sizeof(float)==sizeof(std::uint32_t),"capture requires binary32 float");
  static_assert(std::numeric_limits<float>::is_iec559,"capture requires IEEE float");
  std::uint32_t bits=0;std::memcpy(&bits,&value,sizeof(bits));return bits;
}

inline void bindNPBObject(NPBWitness& witness,std::uint32_t key,std::int64_t ordinal,
                          const std::string& node,std::int32_t policy)
{
  witness.source_identity_available=1;witness.source_native_cluster_key=key;
  witness.source_encounter_ordinal=ordinal;witness.source_cluster_node=node;
  witness.binding_policy=policy;
}

// Observe the pointer returned by the existing native matching function. This
// does not choose a match or change its order; get_id() is never used as a key.
template<class Container,class Photon>
void bindNPBContainerObject(NPBWitness& witness,const Container* container,
                           const Photon* selected,const std::string& node,std::int32_t policy)
{
  if(!container||!selected)return;
  const auto range=container->getClusters();std::int64_t ordinal=0;
  for(auto it=range.first;it!=range.second;++it,++ordinal)
    if(it->second==selected)
    {
      bindNPBObject(witness,static_cast<std::uint32_t>(it->first),ordinal,node,policy);
      return;
    }
}

template<class Photon>
void recordNPB(Photon* photon,const std::string& scoreName,const std::string& modelFile,
               const std::vector<std::string>& names,const std::vector<float>& inputs,
               bool inputsCaptured,bool modelAvailable,bool inDomain,bool evaluated,
               int domainMode,float minEt,float maxEt,float maxAbsEta,float score)
{
  if(!photon)return;
  const std::string prefix=npbPrefix(scoreName);
  auto put=[&](const std::string& name,float value)
  { photon->set_shower_shape_parameter(prefix+name,value); };
  const bool finite=inputsCaptured&&inputs.size()==names.size()&&
      std::all_of(inputs.begin(),inputs.end(),[](float value){return std::isfinite(value);});
  const int state=evaluated?(std::isfinite(score)?EVALUATED_FINITE:EVALUATED_NONFINITE):
      (!modelAvailable?MODEL_UNAVAILABLE:(!inDomain?OUTSIDE_PRODUCER_DOMAIN:INPUT_NONFINITE));
  put("capture_version",1);put("evaluation_state",state);put("domain_mode",domainMode);
  put("in_producer_domain",inDomain?1:0);put("model_available",modelAvailable?1:0);
  put("input_order_is_npb25",names==npbFeatureNames()?1:0);
  put("inputs_captured",inputsCaptured?1:0);put("inputs_finite",finite?1:0);
  put("min_et",minEt);put("max_et",maxEt);put("max_abs_eta",maxAbsEta);put("raw_score",score);
  std::string featureOrder;
  for(const auto& name:names){if(!featureOrder.empty())featureOrder+='|';featureOrder+=name;}
  put("feature_names:"+featureOrder,0);put("model_file:"+modelFile,0);
  // Copy the exact vector passed to the existing Compute(), or its exact
  // pre-compute invalid-input witness. Never recalculate ratios/vertex here.
  if(inputsCaptured)for(std::size_t i=0;i<inputs.size();++i)put("input_"+std::to_string(i),inputs[i]);
}

template<class Photon> NPBWitness readNPB(const Photon* photon,const std::string& scoreName)
{
  NPBWitness result;
  const std::string prefix=npbPrefix(scoreName);
  if(parameter(photon,prefix+"capture_version")!=1)return result;
  result.capture_version=1;
  result.source_photon_id=static_cast<std::uint32_t>(photon->get_id());
  result.source_cluster_et=parameter(photon,"cluster_pt");
  result.source_eta=parameter(photon,"cluster_eta");result.source_phi=parameter(photon,"cluster_phi");
  result.evaluation_state=integerParameter(photon,prefix+"evaluation_state");
  result.domain_mode=integerParameter(photon,prefix+"domain_mode");
  result.in_producer_domain=integerParameter(photon,prefix+"in_producer_domain");
  result.model_available=integerParameter(photon,prefix+"model_available");
  result.inputs_captured=integerParameter(photon,prefix+"inputs_captured");
  result.min_et=parameter(photon,prefix+"min_et");result.max_et=parameter(photon,prefix+"max_et");
  result.max_abs_eta=parameter(photon,prefix+"max_abs_eta");result.raw_score=parameter(photon,prefix+"raw_score");
  result.model_file=textParameter(photon,prefix+"model_file:");
  std::istringstream order(textParameter(photon,prefix+"feature_names:"));
  std::string name;
  while(std::getline(order,name,'|'))result.feature_names.push_back(name);
  result.input_order_is_npb25=result.feature_names==npbFeatureNames()?1:0;
  if(result.inputs_captured==1)for(std::size_t i=0;i<result.feature_names.size();++i)
  {
    const float value=parameter(photon,prefix+"input_"+std::to_string(i));
    result.ordered_inputs.push_back(value);
    result.ordered_input_finite.push_back(std::isfinite(value)?1:0);
  }
  result.inputs_finite=result.inputs_captured==1&&!result.feature_names.empty()&&
      std::all_of(result.ordered_input_finite.begin(),result.ordered_input_finite.end(),
                  [](std::int32_t value){return value==1;})?1:0;
  if(result.input_order_is_npb25==1&&result.ordered_inputs.size()==25)
  {
    result.kinematic_input_bits_available=1;
    result.input_et_bits=floatBits(result.ordered_inputs[0]);
    result.input_eta_bits=floatBits(result.ordered_inputs[1]);
  }
  return result;
}

template<class Photon>
void recordNativeTiming(Photon* photon,float mean,float denominator,float numerator,int towers,
                        const std::string& inputNode,const std::string& towerNode)
{
  photon->set_shower_shape_parameter("rj_native_timing_v1",1);
  photon->set_shower_shape_parameter("rj_native_timing_mean",mean);
  photon->set_shower_shape_parameter("rj_native_timing_denominator",denominator);
  photon->set_shower_shape_parameter("rj_native_timing_numerator",numerator);
  photon->set_shower_shape_parameter("rj_native_timing_towers",static_cast<float>(towers));
  photon->set_shower_shape_parameter("rj_native_timing_input_node:"+inputNode,0);
  photon->set_shower_shape_parameter("rj_native_timing_tower_node:"+towerNode,0);
}

template<class Photon> NativeTiming readNativeTiming(const Photon* photon)
{
  NativeTiming result;
  if(parameter(photon,"rj_native_timing_v1")!=1)return result;
  result.capture_version=1;
  result.mean_time_samples=parameter(photon,"rj_native_timing_mean");
  result.energy_denominator=parameter(photon,"rj_native_timing_denominator");
  result.time_energy_numerator=parameter(photon,"rj_native_timing_numerator");
  result.contributing_towers=integerParameter(photon,"rj_native_timing_towers");
  result.raw_time_finite=std::isfinite(result.mean_time_samples)?1:0;
  result.producer_vertex_z=parameter(photon,"vertex_z");
  result.input_cluster_node=textParameter(photon,"rj_native_timing_input_node:");
  result.calibrated_tower_node=textParameter(photon,"rj_native_timing_tower_node:");
  return result;
}

template<class Mbd> MbdTimes captureMbd(Mbd* mbd)
{
  MbdTimes result;result.capture_version=1;
  if(!mbd)return result;
  result.available=1;result.native_is_valid=mbd->isValid();
  result.t0_ns=mbd->get_t0();result.south_time_ns=mbd->get_time(0);result.north_time_ns=mbd->get_time(1);
  result.south_npmt=mbd->get_npmt(0);result.north_npmt=mbd->get_npmt(1);
  result.t0_finite=std::isfinite(result.t0_ns)?1:0;
  result.south_time_finite=std::isfinite(result.south_time_ns)?1:0;
  result.north_time_finite=std::isfinite(result.north_time_ns)?1:0;
  return result;
}

template<class Tower> TowerTimeStatus captureTower(Tower* tower)
{
  TowerTimeStatus result;result.capture_version=1;
  if(!tower)return result;
  result.available=1;result.time_samples=tower->get_time();
  result.status=static_cast<std::uint32_t>(tower->get_status());
  result.time_finite=std::isfinite(result.time_samples)?1:0;
  return result;
}
} // namespace RJProducerCaptureV1
#endif
