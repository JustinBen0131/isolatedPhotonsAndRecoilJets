#ifndef RJ_REPLAY_FOUNDATION_V1_H
#define RJ_REPLAY_FOUNDATION_V1_H

// Versioned, selection-neutral normalized replay tables for the synchronized
// pp/AuAu photon+jet foundation.  This header is intentionally self-contained
// so the pp and AuAu libraries can still be staged and built independently.

#include <TDirectory.h>
#include <TFile.h>
#include <TNamed.h>
#include <TTree.h>
#include <Compression.h>
#include <RtypesCore.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace RJReplayFoundationV1
{
constexpr const char* kSchemaName = "RJ_REPLAY_FOUNDATION_V1";
constexpr int kSchemaVersion = 2;

// ROOT leaf-list codes '/l' and '/L' are defined in terms of ULong64_t and
// Long64_t, respectively.  std::uint64_t/std::int64_t are intentionally kept
// for the detector-neutral API and identity arithmetic; on LP64 Linux those
// standard types are not necessarily the same C++ types ROOT expects when a
// branch address is bound.
using SerializedUInt64 = ULong64_t;
using SerializedInt64 = Long64_t;
static_assert(sizeof(SerializedUInt64) == sizeof(std::uint64_t),
              "ROOT unsigned 64-bit persistence type changed width");
static_assert(sizeof(SerializedInt64) == sizeof(std::int64_t),
              "ROOT signed 64-bit persistence type changed width");

struct Identity128
{
  std::uint64_t hi = 0;
  std::uint64_t lo = 0;

  bool isNull() const { return hi == 0 && lo == 0; }
  bool operator==(const Identity128& other) const { return hi == other.hi && lo == other.lo; }
  bool operator!=(const Identity128& other) const { return !(*this == other); }
  std::string hex() const
  {
    std::ostringstream os;
    os << std::hex << std::setfill('0') << std::setw(16) << hi << std::setw(16) << lo;
    return os.str();
  }
};

struct SerializedIdentity128
{
  SerializedUInt64 hi = 0;
  SerializedUInt64 lo = 0;
};

struct IdentityHash
{
  std::size_t operator()(const Identity128& id) const noexcept
  {
    return static_cast<std::size_t>(id.hi ^ (id.lo + 0x9e3779b97f4a7c15ULL + (id.hi << 6U) + (id.hi >> 2U)));
  }
};

namespace detail
{
inline std::uint32_t rotr(std::uint32_t x, std::uint32_t n) { return (x >> n) | (x << (32U - n)); }

inline std::array<std::uint8_t, 32> sha256(const std::string& input)
{
  static constexpr std::uint32_t k[64] = {
    0x428a2f98U,0x71374491U,0xb5c0fbcfU,0xe9b5dba5U,0x3956c25bU,0x59f111f1U,0x923f82a4U,0xab1c5ed5U,
    0xd807aa98U,0x12835b01U,0x243185beU,0x550c7dc3U,0x72be5d74U,0x80deb1feU,0x9bdc06a7U,0xc19bf174U,
    0xe49b69c1U,0xefbe4786U,0x0fc19dc6U,0x240ca1ccU,0x2de92c6fU,0x4a7484aaU,0x5cb0a9dcU,0x76f988daU,
    0x983e5152U,0xa831c66dU,0xb00327c8U,0xbf597fc7U,0xc6e00bf3U,0xd5a79147U,0x06ca6351U,0x14292967U,
    0x27b70a85U,0x2e1b2138U,0x4d2c6dfcU,0x53380d13U,0x650a7354U,0x766a0abbU,0x81c2c92eU,0x92722c85U,
    0xa2bfe8a1U,0xa81a664bU,0xc24b8b70U,0xc76c51a3U,0xd192e819U,0xd6990624U,0xf40e3585U,0x106aa070U,
    0x19a4c116U,0x1e376c08U,0x2748774cU,0x34b0bcb5U,0x391c0cb3U,0x4ed8aa4aU,0x5b9cca4fU,0x682e6ff3U,
    0x748f82eeU,0x78a5636fU,0x84c87814U,0x8cc70208U,0x90befffaU,0xa4506cebU,0xbef9a3f7U,0xc67178f2U};
  std::vector<std::uint8_t> msg(input.begin(), input.end());
  const std::uint64_t bitLen = static_cast<std::uint64_t>(msg.size()) * 8ULL;
  msg.push_back(0x80U);
  while ((msg.size() % 64U) != 56U) msg.push_back(0U);
  for (int i = 7; i >= 0; --i) msg.push_back(static_cast<std::uint8_t>((bitLen >> (i * 8)) & 0xffU));

  std::uint32_t h[8] = {0x6a09e667U,0xbb67ae85U,0x3c6ef372U,0xa54ff53aU,
                        0x510e527fU,0x9b05688cU,0x1f83d9abU,0x5be0cd19U};
  for (std::size_t off = 0; off < msg.size(); off += 64U)
  {
    std::uint32_t w[64] = {};
    for (int i = 0; i < 16; ++i)
      w[i] = (static_cast<std::uint32_t>(msg[off + 4*i]) << 24U) |
             (static_cast<std::uint32_t>(msg[off + 4*i + 1]) << 16U) |
             (static_cast<std::uint32_t>(msg[off + 4*i + 2]) << 8U) |
             static_cast<std::uint32_t>(msg[off + 4*i + 3]);
    for (int i = 16; i < 64; ++i)
    {
      const std::uint32_t s0 = rotr(w[i-15],7U) ^ rotr(w[i-15],18U) ^ (w[i-15] >> 3U);
      const std::uint32_t s1 = rotr(w[i-2],17U) ^ rotr(w[i-2],19U) ^ (w[i-2] >> 10U);
      w[i] = w[i-16] + s0 + w[i-7] + s1;
    }
    std::uint32_t a=h[0],b=h[1],c=h[2],d=h[3],e=h[4],f=h[5],g=h[6],hh=h[7];
    for (int i = 0; i < 64; ++i)
    {
      const std::uint32_t s1 = rotr(e,6U) ^ rotr(e,11U) ^ rotr(e,25U);
      const std::uint32_t ch = (e & f) ^ ((~e) & g);
      const std::uint32_t t1 = hh + s1 + ch + k[i] + w[i];
      const std::uint32_t s0 = rotr(a,2U) ^ rotr(a,13U) ^ rotr(a,22U);
      const std::uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
      const std::uint32_t t2 = s0 + maj;
      hh=g; g=f; f=e; e=d+t1; d=c; c=b; b=a; a=t1+t2;
    }
    h[0]+=a; h[1]+=b; h[2]+=c; h[3]+=d; h[4]+=e; h[5]+=f; h[6]+=g; h[7]+=hh;
  }
  std::array<std::uint8_t,32> out{};
  for (int i=0;i<8;++i) for (int j=0;j<4;++j) out[4*i+j]=static_cast<std::uint8_t>((h[i]>>(24-8*j))&0xffU);
  return out;
}
} // namespace detail

inline Identity128 makeIdentity(const std::string& canonicalInput)
{
  const auto digest = detail::sha256(canonicalInput);
  Identity128 id;
  for (int i=0;i<8;++i) id.hi=(id.hi<<8U)|digest[i];
  for (int i=8;i<16;++i) id.lo=(id.lo<<8U)|digest[i];
  if (id.isNull()) id.lo=1; // reserved null identity; SHA-256 collision guard.
  return id;
}

inline std::string sha256Hex(const std::string& canonicalInput)
{
  const auto digest = detail::sha256(canonicalInput);
  std::ostringstream os;
  os << std::hex << std::setfill('0');
  for (const auto byte : digest) os << std::setw(2) << static_cast<unsigned int>(byte);
  return os.str();
}

enum class ModelApplicability : std::int32_t { VALIDATED_DOMAIN=0, DIAGNOSTIC_EXTRAPOLATION=1, MODEL_NOT_APPLICABLE=2, INPUT_INVALID=3 };
enum class LinkClass : std::int32_t { MATCH=0, RECO_FAKE=1, TRUTH_MISS=2, WRONG_PHOTON=3, WRONG_RECOIL=4, MATCH_CANDIDATE=5 };
enum class RecoTruthType : std::int32_t { NONE=0, PHOTON=1, JET=2 };

struct SourceOccurrenceRow { Identity128 id; std::string lane,dataset,sample,period,si_di_role,ownership_state,input_uri_hash,input_file_sha256,source_manifest_sha256; std::int32_t run=0,segment=0; };
struct EventRow { Identity128 id,source_id; std::int32_t run=0; std::int64_t event_sequence=0; std::uint64_t trigger_bits=0; double vertex_z=0,centrality=-1,event_weight=1; std::int32_t terminal_status=0,candidate_count=0,tag_count=0,recoil_count=0; };
struct PhotonCandidateRow { Identity128 id,event_id; std::int32_t encounter_ordinal=0,finite_feature_state=0,below15_retention_state=0; std::vector<double> rank_keys; double cluster_et=0,eta=0,phi=0; std::vector<float> ordered_features; std::vector<std::string> shower_definition_views; std::uint64_t preselection_bitmask=0; };
struct ModelEvaluationRow { Identity128 candidate_id,model_id; std::string model_sha256,shower_definition_id,shower_semantic_sha256; std::vector<float> ordered_input_witnesses; double raw_score=std::numeric_limits<double>::quiet_NaN(); std::int32_t finite_score=0,applicability_state=2; double wp70=std::numeric_limits<double>::quiet_NaN(),wp80=std::numeric_limits<double>::quiet_NaN(),wp90=std::numeric_limits<double>::quiet_NaN(),delta_wp70=std::numeric_limits<double>::quiet_NaN(),delta_wp80=std::numeric_limits<double>::quiet_NaN(),delta_wp90=std::numeric_limits<double>::quiet_NaN(); };
struct ShowerCellRow { Identity128 candidate_id; std::int32_t local_eta_index=0,local_phi_index=0,tower_eta_index=-1,tower_phi_index=-1; std::uint64_t tower_key=0; double calibrated_energy=0,rawcluster_map_value=std::numeric_limits<double>::quiet_NaN(); std::int32_t is_good=0,is_zero=0,is_negative=0,is_nonfinite=0,seed_state=0,denominator_membership=0,rawcluster_owned=0,rawcluster_value_present=0,floor0_membership=0,floor70_membership=0,grid_membership_bitmask=0; };
struct ShowerFeatureViewRow
{
  Identity128 candidate_id,definition_id;
  std::string definition_name,semantic_sha256;
  std::vector<float> ordered_features;
  double floor_gev=0,cog_eta=std::numeric_limits<double>::quiet_NaN(),cog_phi=std::numeric_limits<double>::quiet_NaN();
  double raw_center_eta=std::numeric_limits<double>::quiet_NaN(),raw_center_phi=std::numeric_limits<double>::quiet_NaN();
  std::int32_t center_eta_index=-1,center_phi_index=-1;
  double e11=0,e33=0,e32=0,e35=0,e11_over_e33=std::numeric_limits<double>::quiet_NaN(),e32_over_e35=std::numeric_limits<double>::quiet_NaN();
  double weta_cogx=std::numeric_limits<double>::quiet_NaN(),wphi_cogx=std::numeric_limits<double>::quiet_NaN(),weta33_cogx=std::numeric_limits<double>::quiet_NaN(),wphi33_cogx=std::numeric_limits<double>::quiet_NaN();
  double native_et1=std::numeric_limits<double>::quiet_NaN(),native_et2=std::numeric_limits<double>::quiet_NaN(),native_et3=std::numeric_limits<double>::quiet_NaN(),native_et4=std::numeric_limits<double>::quiet_NaN();
  double moment_eta_numerator=0,moment_phi_numerator=0,moment_denominator=0,moment33_eta_numerator=0,moment33_phi_numerator=0,moment33_denominator=0;
  std::int32_t energy_source=0,rectangular_membership=0,moment_membership=0,finite_feature_state=0;
  std::int32_t good_cell_count=0,owned_cell_count=0,active_sum_cell_count=0,active_moment_cell_count=0,exact_zero_count=0,negative_count=0,nonfinite_count=0;
};
struct IsolationConstituentRow { Identity128 candidate_id,constituent_id; double delta_eta=0,delta_phi=0,delta_r=0,raw_energy=0,calibrated_energy=0,sub1_energy=0,phosub_residual=0; std::int32_t subsystem=0,quality_state=0,mask_state=0,candidate_removal_state=0; };
struct IsolationWitnessRow { Identity128 candidate_id,isolation_id; double radius=0,cone_sum=0,threshold=0; std::int32_t subtraction_method=0,reconstructed_or_truth=0,pass_state=0,constituent_count=0; };
struct JetRow { Identity128 id,event_id; std::string algorithm,input_identity,subtraction_identity; double radius=0,raw_pt=0,corrected_pt=0,eta=0,phi=0; std::uint64_t quality_bitmask=0; std::int32_t deterministic_order=0; };
struct JetConstituentRow { Identity128 jet_id,constituent_id; std::int32_t constituent_ordinal=0,subsystem=0,quality_state=0; double energy=0,eta=0,phi=0; };
struct PhotonJetPairRow { Identity128 id,event_id,candidate_id,jet_id; double delta_phi=0,xjgamma=0; std::int32_t recoil_state=0,photon_rank=0,jet_rank=0,wrong_photon_class=0,wrong_recoil_class=0; };
struct TruthPhotonRow { Identity128 id,event_id; double pt=0,eta=0,phi=0,truth_isolation_witness=0; std::int32_t prompt_class=0,source_role=0,reporting_guard_state=0; };
struct TruthJetRow { Identity128 id,event_id; std::string algorithm,ownership_state; double radius=0,pt=0,eta=0,phi=0; std::int32_t reporting_guard_state=0; };
struct RecoTruthLinkRow { Identity128 id,reco_id,truth_id; std::int32_t reco_type=0,truth_type=0,link_class=0; double match_metric=std::numeric_limits<double>::quiet_NaN(); };
struct WeightComponentRow { Identity128 target_id; std::string component_type; double slice_weight=1,cross_section_weight=1,vertex_weight=1,si_di_weight=1,period_weight=1,exposure_weight=1,final_weight=1; std::int32_t application_count=0; };
struct EventDisplaySnapshotRow { Identity128 id,event_id; std::string selection_reason,quota_class,serialized_payload_hash; };

// Selection-neutral jet matching shared by the p+p and Au+Au replay writers.
// Every same-radius edge inside the frozen candidate gate is retained as a
// MATCH_CANDIDATE row.  MATCH/RECO_FAKE/TRUTH_MISS rows then form an exact,
// disjoint final partition.  Candidate rows are diagnostics and are excluded
// from that partition.
struct JetMatchObject
{
  Identity128 id;
  double pt=0,eta=0,phi=0,radius=0;
};

inline std::vector<RecoTruthLinkRow> buildDeterministicJetLinks(
    const std::vector<JetMatchObject>& reco,
    const std::vector<JetMatchObject>& truth,
    double deltaRMax)
{
  if (!std::isfinite(deltaRMax) || deltaRMax <= 0.0)
    throw std::invalid_argument("jet match deltaR must be finite and positive");

  struct Candidate
  {
    std::size_t recoIndex=0,truthIndex=0;
    double deltaR=0;
  };
  std::vector<Candidate> candidates;
  std::vector<double> nearestReco(reco.size(),std::numeric_limits<double>::quiet_NaN());
  std::vector<double> nearestTruth(truth.size(),std::numeric_limits<double>::quiet_NaN());
  auto identityLess=[](const Identity128& a,const Identity128& b)
  { return a.hi<b.hi || (a.hi==b.hi && a.lo<b.lo); };
  auto updateNearest=[](double& current,double value)
  { if(!std::isfinite(current)||value<current) current=value; };
  auto wrappedDeltaPhi=[](double a,double b)
  { return std::atan2(std::sin(a-b),std::cos(a-b)); };

  for(std::size_t ir=0;ir<reco.size();++ir)
  {
    if(reco[ir].id.isNull()||!std::isfinite(reco[ir].eta)||!std::isfinite(reco[ir].phi))continue;
    for(std::size_t it=0;it<truth.size();++it)
    {
      if(truth[it].id.isNull()||!std::isfinite(truth[it].eta)||!std::isfinite(truth[it].phi))continue;
      if(std::fabs(reco[ir].radius-truth[it].radius)>1.0e-6)continue;
      const double dr=std::hypot(reco[ir].eta-truth[it].eta,
                                 wrappedDeltaPhi(reco[ir].phi,truth[it].phi));
      if(!std::isfinite(dr))continue;
      updateNearest(nearestReco[ir],dr);
      updateNearest(nearestTruth[it],dr);
      if(dr<deltaRMax)candidates.push_back({ir,it,dr});
    }
  }
  std::sort(candidates.begin(),candidates.end(),[&](const Candidate& a,const Candidate& b)
  {
    if(a.deltaR!=b.deltaR)return a.deltaR<b.deltaR;
    if(reco[a.recoIndex].pt!=reco[b.recoIndex].pt)
      return reco[a.recoIndex].pt>reco[b.recoIndex].pt;
    if(reco[a.recoIndex].id!=reco[b.recoIndex].id)
      return identityLess(reco[a.recoIndex].id,reco[b.recoIndex].id);
    return identityLess(truth[a.truthIndex].id,truth[b.truthIndex].id);
  });

  std::vector<RecoTruthLinkRow> links;
  links.reserve(candidates.size()+reco.size()+truth.size());
  for(const Candidate& edge:candidates)
  {
    RecoTruthLinkRow row;
    row.id=makeIdentity(reco[edge.recoIndex].id.hex()+"|"+
                        truth[edge.truthIndex].id.hex()+"|jet_match_candidate");
    row.reco_type=static_cast<int>(RecoTruthType::JET);
    row.reco_id=reco[edge.recoIndex].id;
    row.truth_type=static_cast<int>(RecoTruthType::JET);
    row.truth_id=truth[edge.truthIndex].id;
    row.match_metric=edge.deltaR;
    row.link_class=static_cast<int>(LinkClass::MATCH_CANDIDATE);
    links.push_back(row);
  }

  std::unordered_set<Identity128,IdentityHash> matchedReco,matchedTruth;
  for(const Candidate& edge:candidates)
  {
    const auto& recoObject=reco[edge.recoIndex];
    const auto& truthObject=truth[edge.truthIndex];
    if(matchedReco.count(recoObject.id)||matchedTruth.count(truthObject.id))continue;
    matchedReco.insert(recoObject.id);matchedTruth.insert(truthObject.id);
    RecoTruthLinkRow row;
    row.id=makeIdentity(recoObject.id.hex()+"|"+truthObject.id.hex()+"|jet_match");
    row.reco_type=static_cast<int>(RecoTruthType::JET);row.reco_id=recoObject.id;
    row.truth_type=static_cast<int>(RecoTruthType::JET);row.truth_id=truthObject.id;
    row.match_metric=edge.deltaR;row.link_class=static_cast<int>(LinkClass::MATCH);
    links.push_back(row);
  }
  for(std::size_t ir=0;ir<reco.size();++ir)
  {
    if(matchedReco.count(reco[ir].id))continue;
    RecoTruthLinkRow row;
    row.id=makeIdentity(reco[ir].id.hex()+"|jet_fake");
    row.reco_type=static_cast<int>(RecoTruthType::JET);row.reco_id=reco[ir].id;
    row.truth_type=static_cast<int>(RecoTruthType::NONE);
    row.match_metric=nearestReco[ir];row.link_class=static_cast<int>(LinkClass::RECO_FAKE);
    links.push_back(row);
  }
  for(std::size_t it=0;it<truth.size();++it)
  {
    if(matchedTruth.count(truth[it].id))continue;
    RecoTruthLinkRow row;
    row.id=makeIdentity(truth[it].id.hex()+"|jet_miss");
    row.reco_type=static_cast<int>(RecoTruthType::NONE);
    row.truth_type=static_cast<int>(RecoTruthType::JET);row.truth_id=truth[it].id;
    row.match_metric=nearestTruth[it];row.link_class=static_cast<int>(LinkClass::TRUTH_MISS);
    links.push_back(row);
  }
  return links;
}

struct Metadata
{
  std::string schema_sha256,semantic_sha256,source_sha256,model_sha256,config_sha256,code_sha256;
};

enum class WriterMode
{
  SERIALIZE,
  VALIDATE_ONLY
};

class Writer
{
 public:
  Writer() = default;
  Writer(const Writer&) = delete;
  Writer& operator=(const Writer&) = delete;

  bool initialize(TFile* file, const Metadata& metadata, std::string* error=nullptr)
  {
    return initialize(file,metadata,WriterMode::SERIALIZE,error);
  }

  bool initialize(TFile* file, const Metadata& metadata, WriterMode mode, std::string* error=nullptr)
  {
    if (!file || !file->IsOpen()) return fail(error,"output file is not open");
    if (!validMetadata(metadata)) return fail(error,"metadata hashes must be 64 lowercase hexadecimal characters");
    if(mode!=WriterMode::SERIALIZE&&mode!=WriterMode::VALIDATE_ONLY)return fail(error,"unsupported writer mode");
    m_file=file; m_metadata=metadata; m_mode=mode;
    m_file->SetCompressionAlgorithm(static_cast<int>(ROOT::RCompressionSetting::EAlgorithm::kZSTD));
    m_file->SetCompressionLevel(5);
    if(m_mode==WriterMode::VALIDATE_ONLY){m_initialized=true;return true;}
    TDirectory* saved=gDirectory;
    m_dir=m_file->GetDirectory("ReplayFoundationV1");
    if (!m_dir) m_dir=m_file->mkdir("ReplayFoundationV1");
    if (!m_dir) return fail(error,"failed to create ReplayFoundationV1 directory");
    m_dir->cd();
    bookTrees();
    TNamed schema("rj_replay_schema",kSchemaName); schema.Write("rj_replay_schema",TObject::kOverwrite);
    TNamed schemaVersion("rj_replay_schema_version","2"); schemaVersion.Write("rj_replay_schema_version",TObject::kOverwrite);
    writeMeta("schema_sha256",metadata.schema_sha256); writeMeta("semantic_sha256",metadata.semantic_sha256);
    writeMeta("source_sha256",metadata.source_sha256); writeMeta("model_sha256",metadata.model_sha256);
    writeMeta("config_sha256",metadata.config_sha256); writeMeta("code_sha256",metadata.code_sha256);
    if (saved) saved->cd();
    m_initialized=true;
    return true;
  }

  bool fill(const SourceOccurrenceRow& r,std::string* e=nullptr){ if(!ready(e)||r.id.isNull()||!insert(m_sourceIds,r.id))return fail(e,"invalid or duplicate source identity"); if(validationOnly())return true; m_source=r; serialize(m_sourceId,r.id); m_tSource->Fill(); return true; }
  bool fill(const EventRow& r,std::string* e=nullptr){ if(!ready(e)||r.id.isNull()||!contains(m_sourceIds,r.source_id)||!insert(m_eventIds,r.id))return fail(e,"event foreign key or identity failure"); if(validationOnly())return true; m_event=r; serialize(m_eventId,r.id); serialize(m_eventSourceId,r.source_id); m_eventSequence=serialize(r.event_sequence); m_eventTriggerBits=serialize(r.trigger_bits); m_tEvent->Fill(); return true; }
  bool fill(const PhotonCandidateRow& r,std::string* e=nullptr){ if(!ready(e)||r.id.isNull()||!contains(m_eventIds,r.event_id)||!insert(m_candidateIds,r.id))return fail(e,"candidate foreign key or identity failure"); if(validationOnly())return true; m_candidate=r; serialize(m_candidateId,r.id); serialize(m_candidateEventId,r.event_id); m_candidatePreselectionBitmask=serialize(r.preselection_bitmask); m_tCandidate->Fill(); return true; }
  bool fill(const ModelEvaluationRow& r,std::string* e=nullptr){ if(!ready(e)||r.model_id.isNull()||r.shower_definition_id.empty()||!validHex64(r.shower_semantic_sha256)||!contains(m_candidateIds,r.candidate_id)||!insert(m_modelEvalIds,makeIdentity(r.candidate_id.hex()+"|"+r.model_id.hex())))return fail(e,"model evaluation shower semantics, foreign key, or duplicate failure"); if(validationOnly())return true; m_model=r; serialize(m_modelCandidateId,r.candidate_id); serialize(m_modelId,r.model_id); m_tModel->Fill(); return true; }
  bool fill(const ShowerCellRow& r,std::string* e=nullptr){ const auto key=makeIdentity(r.candidate_id.hex()+"|cell|"+std::to_string(r.tower_key)); const bool gridOrOwnedProvenance=(r.grid_membership_bitmask>0&&r.grid_membership_bitmask<=3)||(r.grid_membership_bitmask==0&&r.rawcluster_owned!=0&&r.rawcluster_value_present!=0); if(!ready(e)||r.tower_eta_index<0||r.tower_eta_index>=96||r.tower_phi_index<0||r.tower_phi_index>=256||!gridOrOwnedProvenance||!contains(m_candidateIds,r.candidate_id)||!insert(m_showerCellIds,key))return fail(e,"shower-cell tower identity, grid-or-owned provenance, candidate foreign key, or duplicate failure"); if(validationOnly())return true; m_shower=r; serialize(m_showerCandidateId,r.candidate_id); m_showerTowerKey=serialize(r.tower_key); m_tShower->Fill(); return true; }
  bool fill(const ShowerFeatureViewRow& r,std::string* e=nullptr){ const auto key=makeIdentity(r.candidate_id.hex()+"|view|"+r.definition_id.hex()); if(!ready(e)||r.definition_id.isNull()||r.definition_name.empty()||!validHex64(r.semantic_sha256)||r.center_eta_index<0||r.center_eta_index>=96||r.center_phi_index<0||r.center_phi_index>=256||!std::isfinite(r.raw_center_eta)||!std::isfinite(r.raw_center_phi)||!contains(m_candidateIds,r.candidate_id)||!insert(m_showerViewIds,key))return fail(e,"shower-feature-view identity, center, semantic hash, foreign key, or duplicate failure"); if(validationOnly())return true; m_showerView=r; serialize(m_showerViewCandidateId,r.candidate_id); serialize(m_showerViewDefinitionId,r.definition_id); m_tShowerView->Fill(); return true; }
  bool fill(const IsolationConstituentRow& r,std::string* e=nullptr){ if(!ready(e)||r.constituent_id.isNull()||!contains(m_candidateIds,r.candidate_id))return fail(e,"isolation constituent foreign key failure"); if(validationOnly())return true; m_isoConstituent=r; serialize(m_isoConstituentCandidateId,r.candidate_id); serialize(m_isoConstituentId,r.constituent_id); m_tIsoConstituent->Fill(); return true; }
  bool fill(const IsolationWitnessRow& r,std::string* e=nullptr){ if(!ready(e)||r.isolation_id.isNull()||!contains(m_candidateIds,r.candidate_id))return fail(e,"isolation witness foreign key failure"); if(validationOnly())return true; m_isoWitness=r; serialize(m_isoWitnessCandidateId,r.candidate_id); serialize(m_isoWitnessId,r.isolation_id); m_tIsoWitness->Fill(); return true; }
  bool fill(const JetRow& r,std::string* e=nullptr){ if(!ready(e)||r.id.isNull()||!contains(m_eventIds,r.event_id)||!insert(m_jetIds,r.id))return fail(e,"jet foreign key or identity failure"); if(validationOnly())return true; m_jet=r; serialize(m_jetId,r.id); serialize(m_jetEventId,r.event_id); m_jetQualityBitmask=serialize(r.quality_bitmask); m_tJet->Fill(); return true; }
  bool fill(const JetConstituentRow& r,std::string* e=nullptr){ if(!ready(e)||r.constituent_id.isNull()||!contains(m_jetIds,r.jet_id))return fail(e,"jet constituent foreign key failure"); if(validationOnly())return true; m_jetConstituent=r; serialize(m_jetConstituentJetId,r.jet_id); serialize(m_jetConstituentId,r.constituent_id); m_tJetConstituent->Fill(); return true; }
  bool fill(const PhotonJetPairRow& r,std::string* e=nullptr){ if(!ready(e)||r.id.isNull()||!contains(m_eventIds,r.event_id)||!contains(m_candidateIds,r.candidate_id)||!contains(m_jetIds,r.jet_id)||!insert(m_pairIds,r.id))return fail(e,"pair foreign key or identity failure"); if(validationOnly())return true; m_pair=r; serialize(m_pairId,r.id); serialize(m_pairEventId,r.event_id); serialize(m_pairCandidateId,r.candidate_id); serialize(m_pairJetId,r.jet_id); m_tPair->Fill(); return true; }
  bool fill(const TruthPhotonRow& r,std::string* e=nullptr){ if(!ready(e)||r.id.isNull()||!contains(m_eventIds,r.event_id)||!insert(m_truthPhotonIds,r.id))return fail(e,"truth photon foreign key or identity failure"); if(validationOnly())return true; m_truthPhoton=r; serialize(m_truthPhotonId,r.id); serialize(m_truthPhotonEventId,r.event_id); m_tTruthPhoton->Fill(); return true; }
  bool fill(const TruthJetRow& r,std::string* e=nullptr){ if(!ready(e)||r.id.isNull()||!contains(m_eventIds,r.event_id)||!insert(m_truthJetIds,r.id))return fail(e,"truth jet foreign key or identity failure"); if(validationOnly())return true; m_truthJet=r; serialize(m_truthJetId,r.id); serialize(m_truthJetEventId,r.event_id); m_tTruthJet->Fill(); return true; }
  bool fill(const RecoTruthLinkRow& r,std::string* e=nullptr){ if(!ready(e)||r.id.isNull()||!validRecoTruthLink(r)||!insert(m_linkIds,r.id))return fail(e,"reco-truth link identity/type/foreign-key failure"); if(validationOnly())return true; m_link=r; serialize(m_linkId,r.id); serialize(m_linkRecoId,r.reco_id); serialize(m_linkTruthId,r.truth_id); m_tLink->Fill(); return true; }
  bool fill(const WeightComponentRow& r,std::string* e=nullptr){ if(!ready(e)||r.target_id.isNull()||r.component_type.empty()||r.application_count<0)return fail(e,"weight component identity/type/count failure"); if(validationOnly())return true; m_weight=r; serialize(m_weightTargetId,r.target_id); m_tWeight->Fill(); return true; }
  bool fill(const EventDisplaySnapshotRow& r,std::string* e=nullptr){ if(!ready(e)||r.id.isNull()||!contains(m_eventIds,r.event_id)||!insert(m_snapshotIds,r.id))return fail(e,"snapshot foreign key or identity failure"); if(validationOnly())return true; m_snapshot=r; serialize(m_snapshotId,r.id); serialize(m_snapshotEventId,r.event_id); m_tSnapshot->Fill(); return true; }

  bool finish(std::string* error=nullptr)
  {
    if(!ready(error)) return false;
    if(validationOnly()){m_finished=true;return true;}
    TDirectory* saved=gDirectory; m_dir->cd();
    for(TTree* tree:m_trees) if(!tree||tree->Write("",TObject::kOverwrite)<=0){ if(saved)saved->cd(); return fail(error,"tree write failure"); }
    TNamed complete("rj_replay_complete","1"); complete.Write("rj_replay_complete",TObject::kOverwrite);
    if(saved) saved->cd();
    m_finished=true;
    return true;
  }

  const Metadata& metadata() const { return m_metadata; }

 private:
  static bool fail(std::string* error,const std::string& message){ if(error)*error=message; return false; }
  bool ready(std::string* e) const { return m_initialized&&!m_finished?true:fail(e,"writer is not active"); }
  bool validationOnly() const { return m_mode==WriterMode::VALIDATE_ONLY; }
  static bool validHex64(const std::string& s){ if(s.size()!=64)return false; for(char c:s)if(!((c>='0'&&c<='9')||(c>='a'&&c<='f')))return false; return true; }
  static bool validMetadata(const Metadata& m){ return validHex64(m.schema_sha256)&&validHex64(m.semantic_sha256)&&validHex64(m.source_sha256)&&validHex64(m.model_sha256)&&validHex64(m.config_sha256)&&validHex64(m.code_sha256); }
  static bool insert(std::unordered_set<Identity128,IdentityHash>& s,const Identity128& id){ return s.insert(id).second; }
  static bool contains(const std::unordered_set<Identity128,IdentityHash>& s,const Identity128& id){ return s.find(id)!=s.end(); }
  static SerializedUInt64 serialize(std::uint64_t value){ return static_cast<SerializedUInt64>(value); }
  static SerializedInt64 serialize(std::int64_t value){ return static_cast<SerializedInt64>(value); }
  static void serialize(SerializedIdentity128& output,const Identity128& input)
  { output.hi=serialize(input.hi); output.lo=serialize(input.lo); }
  void writeMeta(const char* key,const std::string& value){ TNamed obj(key,value.c_str()); obj.Write(key,TObject::kOverwrite); }
  TTree* tree(const char* name){ auto* t=new TTree(name,name); t->SetAutoFlush(-5000000); m_trees.push_back(t); return t; }
  static void unsigned64(TTree* t,const char* n,SerializedUInt64* p){ const std::string leaf=std::string(n)+"/l"; t->Branch(n,p,leaf.c_str()); }
  static void signed64(TTree* t,const char* n,SerializedInt64* p){ const std::string leaf=std::string(n)+"/L"; t->Branch(n,p,leaf.c_str()); }
  static void idBranches(TTree* t,const char* prefix,SerializedIdentity128* id){ const std::string hi=std::string(prefix)+"_hi"; const std::string lo=std::string(prefix)+"_lo"; unsigned64(t,hi.c_str(),&id->hi); unsigned64(t,lo.c_str(),&id->lo); }
  template<class T> static void scalar(TTree* t,const char* n,T* p,const char* type){ const std::string leaf=std::string(n)+"/"+type; t->Branch(n,p,leaf.c_str()); }
  static void str(TTree* t,const char* n,std::string* p){ t->Branch(n,p); }
  template<class T> static void vec(TTree* t,const char* n,std::vector<T>* p){ t->Branch(n,p); }

  void bookTrees()
  {
    m_tSource=tree("RJSourceOccurrenceV1"); idBranches(m_tSource,"source_occurrence_id",&m_sourceId); str(m_tSource,"lane",&m_source.lane); str(m_tSource,"dataset",&m_source.dataset); str(m_tSource,"sample",&m_source.sample); str(m_tSource,"period",&m_source.period); scalar(m_tSource,"run",&m_source.run,"I"); scalar(m_tSource,"segment",&m_source.segment,"I"); str(m_tSource,"si_di_role",&m_source.si_di_role); str(m_tSource,"ownership_state",&m_source.ownership_state); str(m_tSource,"input_uri_hash",&m_source.input_uri_hash); str(m_tSource,"input_file_sha256",&m_source.input_file_sha256); str(m_tSource,"source_manifest_sha256",&m_source.source_manifest_sha256);
    m_tEvent=tree("RJEventV1"); idBranches(m_tEvent,"event_id",&m_eventId); idBranches(m_tEvent,"source_occurrence_id",&m_eventSourceId); scalar(m_tEvent,"run",&m_event.run,"I"); signed64(m_tEvent,"event_sequence",&m_eventSequence); unsigned64(m_tEvent,"trigger_bits",&m_eventTriggerBits); scalar(m_tEvent,"vertex_z",&m_event.vertex_z,"D"); scalar(m_tEvent,"centrality",&m_event.centrality,"D"); scalar(m_tEvent,"event_weight",&m_event.event_weight,"D"); scalar(m_tEvent,"terminal_status",&m_event.terminal_status,"I"); scalar(m_tEvent,"candidate_count",&m_event.candidate_count,"I"); scalar(m_tEvent,"tag_count",&m_event.tag_count,"I"); scalar(m_tEvent,"recoil_count",&m_event.recoil_count,"I");
    m_tCandidate=tree("RJPhotonCandidateV1"); idBranches(m_tCandidate,"candidate_id",&m_candidateId); idBranches(m_tCandidate,"event_id",&m_candidateEventId); scalar(m_tCandidate,"encounter_ordinal",&m_candidate.encounter_ordinal,"I"); vec(m_tCandidate,"rank_keys",&m_candidate.rank_keys); scalar(m_tCandidate,"cluster_et",&m_candidate.cluster_et,"D"); scalar(m_tCandidate,"eta",&m_candidate.eta,"D"); scalar(m_tCandidate,"phi",&m_candidate.phi,"D"); vec(m_tCandidate,"ordered_features",&m_candidate.ordered_features); vec(m_tCandidate,"shower_definition_views",&m_candidate.shower_definition_views); unsigned64(m_tCandidate,"preselection_bitmask",&m_candidatePreselectionBitmask); scalar(m_tCandidate,"finite_feature_state",&m_candidate.finite_feature_state,"I"); scalar(m_tCandidate,"below15_retention_state",&m_candidate.below15_retention_state,"I");
    m_tModel=tree("RJModelEvaluationV1"); idBranches(m_tModel,"candidate_id",&m_modelCandidateId); idBranches(m_tModel,"model_id",&m_modelId); str(m_tModel,"model_sha256",&m_model.model_sha256); str(m_tModel,"shower_definition_id",&m_model.shower_definition_id); str(m_tModel,"shower_semantic_sha256",&m_model.shower_semantic_sha256); vec(m_tModel,"ordered_input_witnesses",&m_model.ordered_input_witnesses); scalar(m_tModel,"raw_score",&m_model.raw_score,"D"); scalar(m_tModel,"finite_score",&m_model.finite_score,"I"); scalar(m_tModel,"applicability_state",&m_model.applicability_state,"I"); scalar(m_tModel,"wp70",&m_model.wp70,"D"); scalar(m_tModel,"wp80",&m_model.wp80,"D"); scalar(m_tModel,"wp90",&m_model.wp90,"D"); scalar(m_tModel,"delta_wp70",&m_model.delta_wp70,"D"); scalar(m_tModel,"delta_wp80",&m_model.delta_wp80,"D"); scalar(m_tModel,"delta_wp90",&m_model.delta_wp90,"D");
    m_tShower=tree("RJShowerCellV1"); idBranches(m_tShower,"candidate_id",&m_showerCandidateId); scalar(m_tShower,"local_eta_index",&m_shower.local_eta_index,"I"); scalar(m_tShower,"local_phi_index",&m_shower.local_phi_index,"I"); scalar(m_tShower,"tower_eta_index",&m_shower.tower_eta_index,"I"); scalar(m_tShower,"tower_phi_index",&m_shower.tower_phi_index,"I"); unsigned64(m_tShower,"tower_key",&m_showerTowerKey); scalar(m_tShower,"calibrated_energy",&m_shower.calibrated_energy,"D"); scalar(m_tShower,"rawcluster_map_value",&m_shower.rawcluster_map_value,"D"); scalar(m_tShower,"is_good",&m_shower.is_good,"I"); scalar(m_tShower,"is_zero",&m_shower.is_zero,"I"); scalar(m_tShower,"is_negative",&m_shower.is_negative,"I"); scalar(m_tShower,"is_nonfinite",&m_shower.is_nonfinite,"I"); scalar(m_tShower,"seed_state",&m_shower.seed_state,"I"); scalar(m_tShower,"denominator_membership",&m_shower.denominator_membership,"I"); scalar(m_tShower,"rawcluster_owned",&m_shower.rawcluster_owned,"I"); scalar(m_tShower,"rawcluster_value_present",&m_shower.rawcluster_value_present,"I"); scalar(m_tShower,"floor0_membership",&m_shower.floor0_membership,"I"); scalar(m_tShower,"floor70_membership",&m_shower.floor70_membership,"I"); scalar(m_tShower,"grid_membership_bitmask",&m_shower.grid_membership_bitmask,"I");
    m_tShowerView=tree("RJShowerFeatureViewV1"); idBranches(m_tShowerView,"candidate_id",&m_showerViewCandidateId); idBranches(m_tShowerView,"definition_id",&m_showerViewDefinitionId); str(m_tShowerView,"definition_name",&m_showerView.definition_name); str(m_tShowerView,"semantic_sha256",&m_showerView.semantic_sha256); vec(m_tShowerView,"ordered_features",&m_showerView.ordered_features); scalar(m_tShowerView,"floor_gev",&m_showerView.floor_gev,"D"); scalar(m_tShowerView,"cog_eta",&m_showerView.cog_eta,"D"); scalar(m_tShowerView,"cog_phi",&m_showerView.cog_phi,"D"); scalar(m_tShowerView,"raw_center_eta",&m_showerView.raw_center_eta,"D"); scalar(m_tShowerView,"raw_center_phi",&m_showerView.raw_center_phi,"D"); scalar(m_tShowerView,"center_eta_index",&m_showerView.center_eta_index,"I"); scalar(m_tShowerView,"center_phi_index",&m_showerView.center_phi_index,"I"); scalar(m_tShowerView,"e11",&m_showerView.e11,"D"); scalar(m_tShowerView,"e33",&m_showerView.e33,"D"); scalar(m_tShowerView,"e32",&m_showerView.e32,"D"); scalar(m_tShowerView,"e35",&m_showerView.e35,"D"); scalar(m_tShowerView,"e11_over_e33",&m_showerView.e11_over_e33,"D"); scalar(m_tShowerView,"e32_over_e35",&m_showerView.e32_over_e35,"D"); scalar(m_tShowerView,"weta_cogx",&m_showerView.weta_cogx,"D"); scalar(m_tShowerView,"wphi_cogx",&m_showerView.wphi_cogx,"D"); scalar(m_tShowerView,"weta33_cogx",&m_showerView.weta33_cogx,"D"); scalar(m_tShowerView,"wphi33_cogx",&m_showerView.wphi33_cogx,"D"); scalar(m_tShowerView,"native_et1",&m_showerView.native_et1,"D"); scalar(m_tShowerView,"native_et2",&m_showerView.native_et2,"D"); scalar(m_tShowerView,"native_et3",&m_showerView.native_et3,"D"); scalar(m_tShowerView,"native_et4",&m_showerView.native_et4,"D"); scalar(m_tShowerView,"moment_eta_numerator",&m_showerView.moment_eta_numerator,"D"); scalar(m_tShowerView,"moment_phi_numerator",&m_showerView.moment_phi_numerator,"D"); scalar(m_tShowerView,"moment_denominator",&m_showerView.moment_denominator,"D"); scalar(m_tShowerView,"moment33_eta_numerator",&m_showerView.moment33_eta_numerator,"D"); scalar(m_tShowerView,"moment33_phi_numerator",&m_showerView.moment33_phi_numerator,"D"); scalar(m_tShowerView,"moment33_denominator",&m_showerView.moment33_denominator,"D"); scalar(m_tShowerView,"energy_source",&m_showerView.energy_source,"I"); scalar(m_tShowerView,"rectangular_membership",&m_showerView.rectangular_membership,"I"); scalar(m_tShowerView,"moment_membership",&m_showerView.moment_membership,"I"); scalar(m_tShowerView,"finite_feature_state",&m_showerView.finite_feature_state,"I"); scalar(m_tShowerView,"good_cell_count",&m_showerView.good_cell_count,"I"); scalar(m_tShowerView,"owned_cell_count",&m_showerView.owned_cell_count,"I"); scalar(m_tShowerView,"active_sum_cell_count",&m_showerView.active_sum_cell_count,"I"); scalar(m_tShowerView,"active_moment_cell_count",&m_showerView.active_moment_cell_count,"I"); scalar(m_tShowerView,"exact_zero_count",&m_showerView.exact_zero_count,"I"); scalar(m_tShowerView,"negative_count",&m_showerView.negative_count,"I"); scalar(m_tShowerView,"nonfinite_count",&m_showerView.nonfinite_count,"I");
    m_tIsoConstituent=tree("RJIsolationConstituentV1"); idBranches(m_tIsoConstituent,"candidate_id",&m_isoConstituentCandidateId); idBranches(m_tIsoConstituent,"constituent_id",&m_isoConstituentId); scalar(m_tIsoConstituent,"delta_eta",&m_isoConstituent.delta_eta,"D"); scalar(m_tIsoConstituent,"delta_phi",&m_isoConstituent.delta_phi,"D"); scalar(m_tIsoConstituent,"delta_r",&m_isoConstituent.delta_r,"D"); scalar(m_tIsoConstituent,"subsystem",&m_isoConstituent.subsystem,"I"); scalar(m_tIsoConstituent,"raw_energy",&m_isoConstituent.raw_energy,"D"); scalar(m_tIsoConstituent,"calibrated_energy",&m_isoConstituent.calibrated_energy,"D"); scalar(m_tIsoConstituent,"sub1_energy",&m_isoConstituent.sub1_energy,"D"); scalar(m_tIsoConstituent,"phosub_residual",&m_isoConstituent.phosub_residual,"D"); scalar(m_tIsoConstituent,"quality_state",&m_isoConstituent.quality_state,"I"); scalar(m_tIsoConstituent,"mask_state",&m_isoConstituent.mask_state,"I"); scalar(m_tIsoConstituent,"candidate_removal_state",&m_isoConstituent.candidate_removal_state,"I");
    m_tIsoWitness=tree("RJIsolationWitnessV1"); idBranches(m_tIsoWitness,"candidate_id",&m_isoWitnessCandidateId); idBranches(m_tIsoWitness,"isolation_identity",&m_isoWitnessId); scalar(m_tIsoWitness,"radius",&m_isoWitness.radius,"D"); scalar(m_tIsoWitness,"subtraction_method",&m_isoWitness.subtraction_method,"I"); scalar(m_tIsoWitness,"reconstructed_or_truth",&m_isoWitness.reconstructed_or_truth,"I"); scalar(m_tIsoWitness,"cone_sum",&m_isoWitness.cone_sum,"D"); scalar(m_tIsoWitness,"threshold",&m_isoWitness.threshold,"D"); scalar(m_tIsoWitness,"pass_state",&m_isoWitness.pass_state,"I"); scalar(m_tIsoWitness,"constituent_count",&m_isoWitness.constituent_count,"I");
    m_tJet=tree("RJJetV1"); idBranches(m_tJet,"jet_id",&m_jetId); idBranches(m_tJet,"event_id",&m_jetEventId); str(m_tJet,"algorithm",&m_jet.algorithm); scalar(m_tJet,"radius",&m_jet.radius,"D"); str(m_tJet,"input_identity",&m_jet.input_identity); str(m_tJet,"subtraction_identity",&m_jet.subtraction_identity); scalar(m_tJet,"raw_pt",&m_jet.raw_pt,"D"); scalar(m_tJet,"corrected_pt",&m_jet.corrected_pt,"D"); scalar(m_tJet,"eta",&m_jet.eta,"D"); scalar(m_tJet,"phi",&m_jet.phi,"D"); unsigned64(m_tJet,"quality_bitmask",&m_jetQualityBitmask); scalar(m_tJet,"deterministic_order",&m_jet.deterministic_order,"I");
    m_tJetConstituent=tree("RJJetConstituentV1"); idBranches(m_tJetConstituent,"jet_id",&m_jetConstituentJetId); idBranches(m_tJetConstituent,"constituent_identity",&m_jetConstituentId); scalar(m_tJetConstituent,"constituent_ordinal",&m_jetConstituent.constituent_ordinal,"I"); scalar(m_tJetConstituent,"subsystem",&m_jetConstituent.subsystem,"I"); scalar(m_tJetConstituent,"energy",&m_jetConstituent.energy,"D"); scalar(m_tJetConstituent,"eta",&m_jetConstituent.eta,"D"); scalar(m_tJetConstituent,"phi",&m_jetConstituent.phi,"D"); scalar(m_tJetConstituent,"quality_state",&m_jetConstituent.quality_state,"I");
    m_tPair=tree("RJPhotonJetPairV1"); idBranches(m_tPair,"pair_id",&m_pairId); idBranches(m_tPair,"event_id",&m_pairEventId); idBranches(m_tPair,"candidate_id",&m_pairCandidateId); idBranches(m_tPair,"jet_id",&m_pairJetId); scalar(m_tPair,"delta_phi",&m_pair.delta_phi,"D"); scalar(m_tPair,"xjgamma",&m_pair.xjgamma,"D"); scalar(m_tPair,"recoil_state",&m_pair.recoil_state,"I"); scalar(m_tPair,"photon_rank",&m_pair.photon_rank,"I"); scalar(m_tPair,"jet_rank",&m_pair.jet_rank,"I"); scalar(m_tPair,"wrong_photon_class",&m_pair.wrong_photon_class,"I"); scalar(m_tPair,"wrong_recoil_class",&m_pair.wrong_recoil_class,"I");
    m_tTruthPhoton=tree("RJTruthPhotonV1"); idBranches(m_tTruthPhoton,"truth_photon_id",&m_truthPhotonId); idBranches(m_tTruthPhoton,"event_id",&m_truthPhotonEventId); scalar(m_tTruthPhoton,"pt",&m_truthPhoton.pt,"D"); scalar(m_tTruthPhoton,"eta",&m_truthPhoton.eta,"D"); scalar(m_tTruthPhoton,"phi",&m_truthPhoton.phi,"D"); scalar(m_tTruthPhoton,"prompt_class",&m_truthPhoton.prompt_class,"I"); scalar(m_tTruthPhoton,"source_role",&m_truthPhoton.source_role,"I"); scalar(m_tTruthPhoton,"truth_isolation_witness",&m_truthPhoton.truth_isolation_witness,"D"); scalar(m_tTruthPhoton,"reporting_guard_state",&m_truthPhoton.reporting_guard_state,"I");
    m_tTruthJet=tree("RJTruthJetV1"); idBranches(m_tTruthJet,"truth_jet_id",&m_truthJetId); idBranches(m_tTruthJet,"event_id",&m_truthJetEventId); str(m_tTruthJet,"algorithm",&m_truthJet.algorithm); scalar(m_tTruthJet,"radius",&m_truthJet.radius,"D"); scalar(m_tTruthJet,"pt",&m_truthJet.pt,"D"); scalar(m_tTruthJet,"eta",&m_truthJet.eta,"D"); scalar(m_tTruthJet,"phi",&m_truthJet.phi,"D"); str(m_tTruthJet,"ownership_state",&m_truthJet.ownership_state); scalar(m_tTruthJet,"reporting_guard_state",&m_truthJet.reporting_guard_state,"I");
    m_tLink=tree("RJRecoTruthLinkV1"); idBranches(m_tLink,"link_id",&m_linkId); scalar(m_tLink,"reco_type",&m_link.reco_type,"I"); idBranches(m_tLink,"reco_id",&m_linkRecoId); scalar(m_tLink,"truth_type",&m_link.truth_type,"I"); idBranches(m_tLink,"truth_id",&m_linkTruthId); scalar(m_tLink,"match_metric",&m_link.match_metric,"D"); scalar(m_tLink,"link_class",&m_link.link_class,"I");
    m_tWeight=tree("RJWeightComponentV1"); idBranches(m_tWeight,"target_id",&m_weightTargetId); str(m_tWeight,"component_type",&m_weight.component_type); scalar(m_tWeight,"slice_weight",&m_weight.slice_weight,"D"); scalar(m_tWeight,"cross_section_weight",&m_weight.cross_section_weight,"D"); scalar(m_tWeight,"vertex_weight",&m_weight.vertex_weight,"D"); scalar(m_tWeight,"si_di_weight",&m_weight.si_di_weight,"D"); scalar(m_tWeight,"period_weight",&m_weight.period_weight,"D"); scalar(m_tWeight,"exposure_weight",&m_weight.exposure_weight,"D"); scalar(m_tWeight,"final_weight",&m_weight.final_weight,"D"); scalar(m_tWeight,"application_count",&m_weight.application_count,"I");
    m_tSnapshot=tree("RJEventDisplaySnapshotV1"); idBranches(m_tSnapshot,"snapshot_id",&m_snapshotId); idBranches(m_tSnapshot,"event_id",&m_snapshotEventId); str(m_tSnapshot,"selection_reason",&m_snapshot.selection_reason); str(m_tSnapshot,"quota_class",&m_snapshot.quota_class); str(m_tSnapshot,"serialized_payload_hash",&m_snapshot.serialized_payload_hash);
  }

  bool validRecoTruthLink(const RecoTruthLinkRow& r) const
  {
    const auto rt=static_cast<RecoTruthType>(r.reco_type), tt=static_cast<RecoTruthType>(r.truth_type);
    const auto lc=static_cast<LinkClass>(r.link_class);
    const bool recoOk=(rt==RecoTruthType::PHOTON&&contains(m_candidateIds,r.reco_id))||(rt==RecoTruthType::JET&&contains(m_jetIds,r.reco_id))||(rt==RecoTruthType::NONE&&r.reco_id.isNull());
    const bool truthOk=(tt==RecoTruthType::PHOTON&&contains(m_truthPhotonIds,r.truth_id))||(tt==RecoTruthType::JET&&contains(m_truthJetIds,r.truth_id))||(tt==RecoTruthType::NONE&&r.truth_id.isNull());
    if(!recoOk||!truthOk) return false;
    if(lc==LinkClass::RECO_FAKE) return rt!=RecoTruthType::NONE&&tt==RecoTruthType::NONE;
    if(lc==LinkClass::TRUTH_MISS) return rt==RecoTruthType::NONE&&tt!=RecoTruthType::NONE;
    return rt!=RecoTruthType::NONE&&tt!=RecoTruthType::NONE;
  }

  TFile* m_file=nullptr; TDirectory* m_dir=nullptr; bool m_initialized=false,m_finished=false; WriterMode m_mode=WriterMode::SERIALIZE; Metadata m_metadata;
  std::vector<TTree*> m_trees;
  TTree *m_tSource=nullptr,*m_tEvent=nullptr,*m_tCandidate=nullptr,*m_tModel=nullptr,*m_tShower=nullptr,*m_tShowerView=nullptr,*m_tIsoConstituent=nullptr,*m_tIsoWitness=nullptr,*m_tJet=nullptr,*m_tJetConstituent=nullptr,*m_tPair=nullptr,*m_tTruthPhoton=nullptr,*m_tTruthJet=nullptr,*m_tLink=nullptr,*m_tWeight=nullptr,*m_tSnapshot=nullptr;
  SourceOccurrenceRow m_source; EventRow m_event; PhotonCandidateRow m_candidate; ModelEvaluationRow m_model; ShowerCellRow m_shower; ShowerFeatureViewRow m_showerView; IsolationConstituentRow m_isoConstituent; IsolationWitnessRow m_isoWitness; JetRow m_jet; JetConstituentRow m_jetConstituent; PhotonJetPairRow m_pair; TruthPhotonRow m_truthPhoton; TruthJetRow m_truthJet; RecoTruthLinkRow m_link; WeightComponentRow m_weight; EventDisplaySnapshotRow m_snapshot;
  SerializedIdentity128 m_sourceId;
  SerializedIdentity128 m_eventId,m_eventSourceId;
  SerializedIdentity128 m_candidateId,m_candidateEventId;
  SerializedIdentity128 m_modelCandidateId,m_modelId;
  SerializedIdentity128 m_showerCandidateId;
  SerializedIdentity128 m_showerViewCandidateId,m_showerViewDefinitionId;
  SerializedIdentity128 m_isoConstituentCandidateId,m_isoConstituentId;
  SerializedIdentity128 m_isoWitnessCandidateId,m_isoWitnessId;
  SerializedIdentity128 m_jetId,m_jetEventId;
  SerializedIdentity128 m_jetConstituentJetId,m_jetConstituentId;
  SerializedIdentity128 m_pairId,m_pairEventId,m_pairCandidateId,m_pairJetId;
  SerializedIdentity128 m_truthPhotonId,m_truthPhotonEventId;
  SerializedIdentity128 m_truthJetId,m_truthJetEventId;
  SerializedIdentity128 m_linkId,m_linkRecoId,m_linkTruthId;
  SerializedIdentity128 m_weightTargetId;
  SerializedIdentity128 m_snapshotId,m_snapshotEventId;
  SerializedInt64 m_eventSequence=0;
  SerializedUInt64 m_eventTriggerBits=0,m_candidatePreselectionBitmask=0;
  SerializedUInt64 m_showerTowerKey=0,m_jetQualityBitmask=0;
  std::unordered_set<Identity128,IdentityHash> m_sourceIds,m_eventIds,m_candidateIds,m_modelEvalIds,m_showerCellIds,m_showerViewIds,m_jetIds,m_pairIds,m_truthPhotonIds,m_truthJetIds,m_linkIds,m_snapshotIds;
};
} // namespace RJReplayFoundationV1

#endif
