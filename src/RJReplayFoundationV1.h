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

#include <array>
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
constexpr int kSchemaVersion = 1;

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

enum class ModelApplicability : std::int32_t { VALIDATED_DOMAIN=0, DIAGNOSTIC_EXTRAPOLATION=1, MODEL_NOT_APPLICABLE=2, INPUT_INVALID=3 };
enum class LinkClass : std::int32_t { MATCH=0, RECO_FAKE=1, TRUTH_MISS=2, WRONG_PHOTON=3, WRONG_RECOIL=4 };
enum class RecoTruthType : std::int32_t { NONE=0, PHOTON=1, JET=2 };

struct SourceOccurrenceRow { Identity128 id; std::string lane,dataset,sample,period,si_di_role,ownership_state,input_uri_hash,input_file_sha256,source_manifest_sha256; std::int32_t run=0,segment=0; };
struct EventRow { Identity128 id,source_id; std::int32_t run=0; std::int64_t event_sequence=0; std::uint64_t trigger_bits=0; double vertex_z=0,centrality=-1,event_weight=1; std::int32_t terminal_status=0,candidate_count=0,tag_count=0,recoil_count=0; };
struct PhotonCandidateRow { Identity128 id,event_id; std::int32_t encounter_ordinal=0,finite_feature_state=0,below15_retention_state=0; std::vector<double> rank_keys; double cluster_et=0,eta=0,phi=0; std::vector<float> ordered_features; std::uint64_t preselection_bitmask=0; };
struct ModelEvaluationRow { Identity128 candidate_id,model_id; std::string model_sha256; std::vector<float> ordered_input_witnesses; double raw_score=std::numeric_limits<double>::quiet_NaN(); std::int32_t finite_score=0,applicability_state=2; double wp70=std::numeric_limits<double>::quiet_NaN(),wp80=std::numeric_limits<double>::quiet_NaN(),wp90=std::numeric_limits<double>::quiet_NaN(),delta_wp70=std::numeric_limits<double>::quiet_NaN(),delta_wp80=std::numeric_limits<double>::quiet_NaN(),delta_wp90=std::numeric_limits<double>::quiet_NaN(); };
struct ShowerCellRow { Identity128 candidate_id; std::int32_t local_eta_index=0,local_phi_index=0; std::uint64_t tower_key=0; double calibrated_energy=0; std::int32_t is_good=0,is_zero=0,is_negative=0,is_nonfinite=0,seed_state=0,denominator_membership=0; };
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

struct Metadata
{
  std::string schema_sha256,semantic_sha256,source_sha256,model_sha256,config_sha256,code_sha256;
};

class Writer
{
 public:
  Writer() = default;
  Writer(const Writer&) = delete;
  Writer& operator=(const Writer&) = delete;

  bool initialize(TFile* file, const Metadata& metadata, std::string* error=nullptr)
  {
    if (!file || !file->IsOpen()) return fail(error,"output file is not open");
    if (!validMetadata(metadata)) return fail(error,"metadata hashes must be 64 lowercase hexadecimal characters");
    m_file=file; m_metadata=metadata;
    m_file->SetCompressionAlgorithm(static_cast<int>(ROOT::RCompressionSetting::EAlgorithm::kZSTD));
    m_file->SetCompressionLevel(5);
    TDirectory* saved=gDirectory;
    m_dir=m_file->GetDirectory("ReplayFoundationV1");
    if (!m_dir) m_dir=m_file->mkdir("ReplayFoundationV1");
    if (!m_dir) return fail(error,"failed to create ReplayFoundationV1 directory");
    m_dir->cd();
    bookTrees();
    TNamed schema("rj_replay_schema",kSchemaName); schema.Write("rj_replay_schema",TObject::kOverwrite);
    TNamed schemaVersion("rj_replay_schema_version","1"); schemaVersion.Write("rj_replay_schema_version",TObject::kOverwrite);
    writeMeta("schema_sha256",metadata.schema_sha256); writeMeta("semantic_sha256",metadata.semantic_sha256);
    writeMeta("source_sha256",metadata.source_sha256); writeMeta("model_sha256",metadata.model_sha256);
    writeMeta("config_sha256",metadata.config_sha256); writeMeta("code_sha256",metadata.code_sha256);
    if (saved) saved->cd();
    m_initialized=true;
    return true;
  }

  bool fill(const SourceOccurrenceRow& r,std::string* e=nullptr){ if(!ready(e)||r.id.isNull()||!insert(m_sourceIds,r.id))return fail(e,"invalid or duplicate source identity"); m_source=r; m_tSource->Fill(); return true; }
  bool fill(const EventRow& r,std::string* e=nullptr){ if(!ready(e)||r.id.isNull()||!contains(m_sourceIds,r.source_id)||!insert(m_eventIds,r.id))return fail(e,"event foreign key or identity failure"); m_event=r; m_tEvent->Fill(); return true; }
  bool fill(const PhotonCandidateRow& r,std::string* e=nullptr){ if(!ready(e)||r.id.isNull()||!contains(m_eventIds,r.event_id)||!insert(m_candidateIds,r.id))return fail(e,"candidate foreign key or identity failure"); m_candidate=r; m_tCandidate->Fill(); return true; }
  bool fill(const ModelEvaluationRow& r,std::string* e=nullptr){ if(!ready(e)||r.model_id.isNull()||!contains(m_candidateIds,r.candidate_id)||!insert(m_modelEvalIds,makeIdentity(r.candidate_id.hex()+"|"+r.model_id.hex())))return fail(e,"model evaluation foreign key or duplicate failure"); m_model=r; m_tModel->Fill(); return true; }
  bool fill(const ShowerCellRow& r,std::string* e=nullptr){ if(!ready(e)||!contains(m_candidateIds,r.candidate_id))return fail(e,"shower-cell candidate foreign key failure"); m_shower=r; m_tShower->Fill(); return true; }
  bool fill(const IsolationConstituentRow& r,std::string* e=nullptr){ if(!ready(e)||r.constituent_id.isNull()||!contains(m_candidateIds,r.candidate_id))return fail(e,"isolation constituent foreign key failure"); m_isoConstituent=r; m_tIsoConstituent->Fill(); return true; }
  bool fill(const IsolationWitnessRow& r,std::string* e=nullptr){ if(!ready(e)||r.isolation_id.isNull()||!contains(m_candidateIds,r.candidate_id))return fail(e,"isolation witness foreign key failure"); m_isoWitness=r; m_tIsoWitness->Fill(); return true; }
  bool fill(const JetRow& r,std::string* e=nullptr){ if(!ready(e)||r.id.isNull()||!contains(m_eventIds,r.event_id)||!insert(m_jetIds,r.id))return fail(e,"jet foreign key or identity failure"); m_jet=r; m_tJet->Fill(); return true; }
  bool fill(const JetConstituentRow& r,std::string* e=nullptr){ if(!ready(e)||r.constituent_id.isNull()||!contains(m_jetIds,r.jet_id))return fail(e,"jet constituent foreign key failure"); m_jetConstituent=r; m_tJetConstituent->Fill(); return true; }
  bool fill(const PhotonJetPairRow& r,std::string* e=nullptr){ if(!ready(e)||r.id.isNull()||!contains(m_eventIds,r.event_id)||!contains(m_candidateIds,r.candidate_id)||!contains(m_jetIds,r.jet_id)||!insert(m_pairIds,r.id))return fail(e,"pair foreign key or identity failure"); m_pair=r; m_tPair->Fill(); return true; }
  bool fill(const TruthPhotonRow& r,std::string* e=nullptr){ if(!ready(e)||r.id.isNull()||!contains(m_eventIds,r.event_id)||!insert(m_truthPhotonIds,r.id))return fail(e,"truth photon foreign key or identity failure"); m_truthPhoton=r; m_tTruthPhoton->Fill(); return true; }
  bool fill(const TruthJetRow& r,std::string* e=nullptr){ if(!ready(e)||r.id.isNull()||!contains(m_eventIds,r.event_id)||!insert(m_truthJetIds,r.id))return fail(e,"truth jet foreign key or identity failure"); m_truthJet=r; m_tTruthJet->Fill(); return true; }
  bool fill(const RecoTruthLinkRow& r,std::string* e=nullptr){ if(!ready(e)||r.id.isNull()||!validRecoTruthLink(r)||!insert(m_linkIds,r.id))return fail(e,"reco-truth link identity/type/foreign-key failure"); m_link=r; m_tLink->Fill(); return true; }
  bool fill(const WeightComponentRow& r,std::string* e=nullptr){ if(!ready(e)||r.target_id.isNull()||r.component_type.empty()||r.application_count<0)return fail(e,"weight component identity/type/count failure"); m_weight=r; m_tWeight->Fill(); return true; }
  bool fill(const EventDisplaySnapshotRow& r,std::string* e=nullptr){ if(!ready(e)||r.id.isNull()||!contains(m_eventIds,r.event_id)||!insert(m_snapshotIds,r.id))return fail(e,"snapshot foreign key or identity failure"); m_snapshot=r; m_tSnapshot->Fill(); return true; }

  bool finish(std::string* error=nullptr)
  {
    if(!ready(error)) return false;
    TDirectory* saved=gDirectory; m_dir->cd();
    for(TTree* tree:m_trees) if(!tree||tree->Write("",TObject::kOverwrite)<=0){ if(saved)saved->cd(); return fail(error,"tree write failure"); }
    TNamed complete("rj_replay_complete","1"); complete.Write("rj_replay_complete",TObject::kOverwrite);
    if(saved) saved->cd();
    m_finished=true;
    return true;
  }

 private:
  static bool fail(std::string* error,const std::string& message){ if(error)*error=message; return false; }
  bool ready(std::string* e) const { return m_initialized&&!m_finished?true:fail(e,"writer is not active"); }
  static bool validHex64(const std::string& s){ if(s.size()!=64)return false; for(char c:s)if(!((c>='0'&&c<='9')||(c>='a'&&c<='f')))return false; return true; }
  static bool validMetadata(const Metadata& m){ return validHex64(m.schema_sha256)&&validHex64(m.semantic_sha256)&&validHex64(m.source_sha256)&&validHex64(m.model_sha256)&&validHex64(m.config_sha256)&&validHex64(m.code_sha256); }
  static bool insert(std::unordered_set<Identity128,IdentityHash>& s,const Identity128& id){ return s.insert(id).second; }
  static bool contains(const std::unordered_set<Identity128,IdentityHash>& s,const Identity128& id){ return s.find(id)!=s.end(); }
  void writeMeta(const char* key,const std::string& value){ TNamed obj(key,value.c_str()); obj.Write(key,TObject::kOverwrite); }
  TTree* tree(const char* name){ auto* t=new TTree(name,name); t->SetAutoFlush(-5000000); m_trees.push_back(t); return t; }
  static void idBranches(TTree* t,const char* prefix,Identity128* id){ const std::string hi=std::string(prefix)+"_hi"; const std::string lo=std::string(prefix)+"_lo"; t->Branch(hi.c_str(),&id->hi,(hi+"/l").c_str()); t->Branch(lo.c_str(),&id->lo,(lo+"/l").c_str()); }
  template<class T> static void scalar(TTree* t,const char* n,T* p,const char* type){ const std::string leaf=std::string(n)+"/"+type; t->Branch(n,p,leaf.c_str()); }
  static void str(TTree* t,const char* n,std::string* p){ t->Branch(n,p); }
  template<class T> static void vec(TTree* t,const char* n,std::vector<T>* p){ t->Branch(n,p); }

  void bookTrees()
  {
    m_tSource=tree("RJSourceOccurrenceV1"); idBranches(m_tSource,"source_occurrence_id",&m_source.id); str(m_tSource,"lane",&m_source.lane); str(m_tSource,"dataset",&m_source.dataset); str(m_tSource,"sample",&m_source.sample); str(m_tSource,"period",&m_source.period); scalar(m_tSource,"run",&m_source.run,"I"); scalar(m_tSource,"segment",&m_source.segment,"I"); str(m_tSource,"si_di_role",&m_source.si_di_role); str(m_tSource,"ownership_state",&m_source.ownership_state); str(m_tSource,"input_uri_hash",&m_source.input_uri_hash); str(m_tSource,"input_file_sha256",&m_source.input_file_sha256); str(m_tSource,"source_manifest_sha256",&m_source.source_manifest_sha256);
    m_tEvent=tree("RJEventV1"); idBranches(m_tEvent,"event_id",&m_event.id); idBranches(m_tEvent,"source_occurrence_id",&m_event.source_id); scalar(m_tEvent,"run",&m_event.run,"I"); scalar(m_tEvent,"event_sequence",&m_event.event_sequence,"L"); scalar(m_tEvent,"trigger_bits",&m_event.trigger_bits,"l"); scalar(m_tEvent,"vertex_z",&m_event.vertex_z,"D"); scalar(m_tEvent,"centrality",&m_event.centrality,"D"); scalar(m_tEvent,"event_weight",&m_event.event_weight,"D"); scalar(m_tEvent,"terminal_status",&m_event.terminal_status,"I"); scalar(m_tEvent,"candidate_count",&m_event.candidate_count,"I"); scalar(m_tEvent,"tag_count",&m_event.tag_count,"I"); scalar(m_tEvent,"recoil_count",&m_event.recoil_count,"I");
    m_tCandidate=tree("RJPhotonCandidateV1"); idBranches(m_tCandidate,"candidate_id",&m_candidate.id); idBranches(m_tCandidate,"event_id",&m_candidate.event_id); scalar(m_tCandidate,"encounter_ordinal",&m_candidate.encounter_ordinal,"I"); vec(m_tCandidate,"rank_keys",&m_candidate.rank_keys); scalar(m_tCandidate,"cluster_et",&m_candidate.cluster_et,"D"); scalar(m_tCandidate,"eta",&m_candidate.eta,"D"); scalar(m_tCandidate,"phi",&m_candidate.phi,"D"); vec(m_tCandidate,"ordered_features",&m_candidate.ordered_features); scalar(m_tCandidate,"preselection_bitmask",&m_candidate.preselection_bitmask,"l"); scalar(m_tCandidate,"finite_feature_state",&m_candidate.finite_feature_state,"I"); scalar(m_tCandidate,"below15_retention_state",&m_candidate.below15_retention_state,"I");
    m_tModel=tree("RJModelEvaluationV1"); idBranches(m_tModel,"candidate_id",&m_model.candidate_id); idBranches(m_tModel,"model_id",&m_model.model_id); str(m_tModel,"model_sha256",&m_model.model_sha256); vec(m_tModel,"ordered_input_witnesses",&m_model.ordered_input_witnesses); scalar(m_tModel,"raw_score",&m_model.raw_score,"D"); scalar(m_tModel,"finite_score",&m_model.finite_score,"I"); scalar(m_tModel,"applicability_state",&m_model.applicability_state,"I"); scalar(m_tModel,"wp70",&m_model.wp70,"D"); scalar(m_tModel,"wp80",&m_model.wp80,"D"); scalar(m_tModel,"wp90",&m_model.wp90,"D"); scalar(m_tModel,"delta_wp70",&m_model.delta_wp70,"D"); scalar(m_tModel,"delta_wp80",&m_model.delta_wp80,"D"); scalar(m_tModel,"delta_wp90",&m_model.delta_wp90,"D");
    m_tShower=tree("RJShowerCellV1"); idBranches(m_tShower,"candidate_id",&m_shower.candidate_id); scalar(m_tShower,"local_eta_index",&m_shower.local_eta_index,"I"); scalar(m_tShower,"local_phi_index",&m_shower.local_phi_index,"I"); scalar(m_tShower,"tower_key",&m_shower.tower_key,"l"); scalar(m_tShower,"calibrated_energy",&m_shower.calibrated_energy,"D"); scalar(m_tShower,"is_good",&m_shower.is_good,"I"); scalar(m_tShower,"is_zero",&m_shower.is_zero,"I"); scalar(m_tShower,"is_negative",&m_shower.is_negative,"I"); scalar(m_tShower,"is_nonfinite",&m_shower.is_nonfinite,"I"); scalar(m_tShower,"seed_state",&m_shower.seed_state,"I"); scalar(m_tShower,"denominator_membership",&m_shower.denominator_membership,"I");
    m_tIsoConstituent=tree("RJIsolationConstituentV1"); idBranches(m_tIsoConstituent,"candidate_id",&m_isoConstituent.candidate_id); idBranches(m_tIsoConstituent,"constituent_id",&m_isoConstituent.constituent_id); scalar(m_tIsoConstituent,"delta_eta",&m_isoConstituent.delta_eta,"D"); scalar(m_tIsoConstituent,"delta_phi",&m_isoConstituent.delta_phi,"D"); scalar(m_tIsoConstituent,"delta_r",&m_isoConstituent.delta_r,"D"); scalar(m_tIsoConstituent,"subsystem",&m_isoConstituent.subsystem,"I"); scalar(m_tIsoConstituent,"raw_energy",&m_isoConstituent.raw_energy,"D"); scalar(m_tIsoConstituent,"calibrated_energy",&m_isoConstituent.calibrated_energy,"D"); scalar(m_tIsoConstituent,"sub1_energy",&m_isoConstituent.sub1_energy,"D"); scalar(m_tIsoConstituent,"phosub_residual",&m_isoConstituent.phosub_residual,"D"); scalar(m_tIsoConstituent,"quality_state",&m_isoConstituent.quality_state,"I"); scalar(m_tIsoConstituent,"mask_state",&m_isoConstituent.mask_state,"I"); scalar(m_tIsoConstituent,"candidate_removal_state",&m_isoConstituent.candidate_removal_state,"I");
    m_tIsoWitness=tree("RJIsolationWitnessV1"); idBranches(m_tIsoWitness,"candidate_id",&m_isoWitness.candidate_id); idBranches(m_tIsoWitness,"isolation_identity",&m_isoWitness.isolation_id); scalar(m_tIsoWitness,"radius",&m_isoWitness.radius,"D"); scalar(m_tIsoWitness,"subtraction_method",&m_isoWitness.subtraction_method,"I"); scalar(m_tIsoWitness,"reconstructed_or_truth",&m_isoWitness.reconstructed_or_truth,"I"); scalar(m_tIsoWitness,"cone_sum",&m_isoWitness.cone_sum,"D"); scalar(m_tIsoWitness,"threshold",&m_isoWitness.threshold,"D"); scalar(m_tIsoWitness,"pass_state",&m_isoWitness.pass_state,"I"); scalar(m_tIsoWitness,"constituent_count",&m_isoWitness.constituent_count,"I");
    m_tJet=tree("RJJetV1"); idBranches(m_tJet,"jet_id",&m_jet.id); idBranches(m_tJet,"event_id",&m_jet.event_id); str(m_tJet,"algorithm",&m_jet.algorithm); scalar(m_tJet,"radius",&m_jet.radius,"D"); str(m_tJet,"input_identity",&m_jet.input_identity); str(m_tJet,"subtraction_identity",&m_jet.subtraction_identity); scalar(m_tJet,"raw_pt",&m_jet.raw_pt,"D"); scalar(m_tJet,"corrected_pt",&m_jet.corrected_pt,"D"); scalar(m_tJet,"eta",&m_jet.eta,"D"); scalar(m_tJet,"phi",&m_jet.phi,"D"); scalar(m_tJet,"quality_bitmask",&m_jet.quality_bitmask,"l"); scalar(m_tJet,"deterministic_order",&m_jet.deterministic_order,"I");
    m_tJetConstituent=tree("RJJetConstituentV1"); idBranches(m_tJetConstituent,"jet_id",&m_jetConstituent.jet_id); idBranches(m_tJetConstituent,"constituent_identity",&m_jetConstituent.constituent_id); scalar(m_tJetConstituent,"constituent_ordinal",&m_jetConstituent.constituent_ordinal,"I"); scalar(m_tJetConstituent,"subsystem",&m_jetConstituent.subsystem,"I"); scalar(m_tJetConstituent,"energy",&m_jetConstituent.energy,"D"); scalar(m_tJetConstituent,"eta",&m_jetConstituent.eta,"D"); scalar(m_tJetConstituent,"phi",&m_jetConstituent.phi,"D"); scalar(m_tJetConstituent,"quality_state",&m_jetConstituent.quality_state,"I");
    m_tPair=tree("RJPhotonJetPairV1"); idBranches(m_tPair,"pair_id",&m_pair.id); idBranches(m_tPair,"event_id",&m_pair.event_id); idBranches(m_tPair,"candidate_id",&m_pair.candidate_id); idBranches(m_tPair,"jet_id",&m_pair.jet_id); scalar(m_tPair,"delta_phi",&m_pair.delta_phi,"D"); scalar(m_tPair,"xjgamma",&m_pair.xjgamma,"D"); scalar(m_tPair,"recoil_state",&m_pair.recoil_state,"I"); scalar(m_tPair,"photon_rank",&m_pair.photon_rank,"I"); scalar(m_tPair,"jet_rank",&m_pair.jet_rank,"I"); scalar(m_tPair,"wrong_photon_class",&m_pair.wrong_photon_class,"I"); scalar(m_tPair,"wrong_recoil_class",&m_pair.wrong_recoil_class,"I");
    m_tTruthPhoton=tree("RJTruthPhotonV1"); idBranches(m_tTruthPhoton,"truth_photon_id",&m_truthPhoton.id); idBranches(m_tTruthPhoton,"event_id",&m_truthPhoton.event_id); scalar(m_tTruthPhoton,"pt",&m_truthPhoton.pt,"D"); scalar(m_tTruthPhoton,"eta",&m_truthPhoton.eta,"D"); scalar(m_tTruthPhoton,"phi",&m_truthPhoton.phi,"D"); scalar(m_tTruthPhoton,"prompt_class",&m_truthPhoton.prompt_class,"I"); scalar(m_tTruthPhoton,"source_role",&m_truthPhoton.source_role,"I"); scalar(m_tTruthPhoton,"truth_isolation_witness",&m_truthPhoton.truth_isolation_witness,"D"); scalar(m_tTruthPhoton,"reporting_guard_state",&m_truthPhoton.reporting_guard_state,"I");
    m_tTruthJet=tree("RJTruthJetV1"); idBranches(m_tTruthJet,"truth_jet_id",&m_truthJet.id); idBranches(m_tTruthJet,"event_id",&m_truthJet.event_id); str(m_tTruthJet,"algorithm",&m_truthJet.algorithm); scalar(m_tTruthJet,"radius",&m_truthJet.radius,"D"); scalar(m_tTruthJet,"pt",&m_truthJet.pt,"D"); scalar(m_tTruthJet,"eta",&m_truthJet.eta,"D"); scalar(m_tTruthJet,"phi",&m_truthJet.phi,"D"); str(m_tTruthJet,"ownership_state",&m_truthJet.ownership_state); scalar(m_tTruthJet,"reporting_guard_state",&m_truthJet.reporting_guard_state,"I");
    m_tLink=tree("RJRecoTruthLinkV1"); idBranches(m_tLink,"link_id",&m_link.id); scalar(m_tLink,"reco_type",&m_link.reco_type,"I"); idBranches(m_tLink,"reco_id",&m_link.reco_id); scalar(m_tLink,"truth_type",&m_link.truth_type,"I"); idBranches(m_tLink,"truth_id",&m_link.truth_id); scalar(m_tLink,"match_metric",&m_link.match_metric,"D"); scalar(m_tLink,"link_class",&m_link.link_class,"I");
    m_tWeight=tree("RJWeightComponentV1"); idBranches(m_tWeight,"target_id",&m_weight.target_id); str(m_tWeight,"component_type",&m_weight.component_type); scalar(m_tWeight,"slice_weight",&m_weight.slice_weight,"D"); scalar(m_tWeight,"cross_section_weight",&m_weight.cross_section_weight,"D"); scalar(m_tWeight,"vertex_weight",&m_weight.vertex_weight,"D"); scalar(m_tWeight,"si_di_weight",&m_weight.si_di_weight,"D"); scalar(m_tWeight,"period_weight",&m_weight.period_weight,"D"); scalar(m_tWeight,"exposure_weight",&m_weight.exposure_weight,"D"); scalar(m_tWeight,"final_weight",&m_weight.final_weight,"D"); scalar(m_tWeight,"application_count",&m_weight.application_count,"I");
    m_tSnapshot=tree("RJEventDisplaySnapshotV1"); idBranches(m_tSnapshot,"snapshot_id",&m_snapshot.id); idBranches(m_tSnapshot,"event_id",&m_snapshot.event_id); str(m_tSnapshot,"selection_reason",&m_snapshot.selection_reason); str(m_tSnapshot,"quota_class",&m_snapshot.quota_class); str(m_tSnapshot,"serialized_payload_hash",&m_snapshot.serialized_payload_hash);
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

  TFile* m_file=nullptr; TDirectory* m_dir=nullptr; bool m_initialized=false,m_finished=false; Metadata m_metadata;
  std::vector<TTree*> m_trees;
  TTree *m_tSource=nullptr,*m_tEvent=nullptr,*m_tCandidate=nullptr,*m_tModel=nullptr,*m_tShower=nullptr,*m_tIsoConstituent=nullptr,*m_tIsoWitness=nullptr,*m_tJet=nullptr,*m_tJetConstituent=nullptr,*m_tPair=nullptr,*m_tTruthPhoton=nullptr,*m_tTruthJet=nullptr,*m_tLink=nullptr,*m_tWeight=nullptr,*m_tSnapshot=nullptr;
  SourceOccurrenceRow m_source; EventRow m_event; PhotonCandidateRow m_candidate; ModelEvaluationRow m_model; ShowerCellRow m_shower; IsolationConstituentRow m_isoConstituent; IsolationWitnessRow m_isoWitness; JetRow m_jet; JetConstituentRow m_jetConstituent; PhotonJetPairRow m_pair; TruthPhotonRow m_truthPhoton; TruthJetRow m_truthJet; RecoTruthLinkRow m_link; WeightComponentRow m_weight; EventDisplaySnapshotRow m_snapshot;
  std::unordered_set<Identity128,IdentityHash> m_sourceIds,m_eventIds,m_candidateIds,m_modelEvalIds,m_jetIds,m_pairIds,m_truthPhotonIds,m_truthJetIds,m_linkIds,m_snapshotIds;
};
} // namespace RJReplayFoundationV1

#endif
