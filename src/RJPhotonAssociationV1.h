#ifndef RJ_PHOTON_ASSOCIATION_V1_H
#define RJ_PHOTON_ASSOCIATION_V1_H

#include "RJDominantTruthWitnessV1.h"

#include <cmath>
#include <cstdint>
#include <iterator>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

// Pure, single-source/event photon association. No kinematic, class, isolation,
// vertex, or angular cut belongs here. The caller supplies the selected truth
// IDs separately from the FULL retained raw-photon inventory. Source-local row
// identities must never be pooled across calls/sources. Native dominant-primary
// witnesses remain the association provenance, including selected-domain fakes.
namespace RJPhotonAssociationV1
{
// Per-event evidence for an explicitly declared reconstructed-object domain.
// The caller owns the cuts: this census never changes a selection or equates
// nodes_ready/terminal success with completeness. Inventory every container
// key independently, then account for each serialized eligible key. A missing
// container must not masquerade as an inspected, legitimately empty one.
class RecoCaptureCensus
{
 public:
  enum class Disposition { ELIGIBLE, OUTSIDE_DOMAIN, UNCLASSIFIABLE };
  void inspect(std::uint32_t key,Disposition disposition)
  {
    if(disposition!=Disposition::ELIGIBLE && disposition!=Disposition::OUTSIDE_DOMAIN &&
       disposition!=Disposition::UNCLASSIFIABLE)fail("unknown reco capture disposition");
    if(m_finished || !m_seen.emplace(key,disposition).second)
      fail("duplicate or late reco capture inventory key");
    if(disposition==Disposition::ELIGIBLE)m_eligible.insert(key);
    else if(disposition==Disposition::UNCLASSIFIABLE)++m_unknown;
  }
  void serialized(std::uint32_t key)
  {
    if(m_finished || !m_eligible.count(key) || !m_written.insert(key).second)
      fail("foreign, excluded, duplicate or late serialized reco key");
  }
  void finish(bool containerAvailable,std::size_t expectedEntries)
  {
    if(m_finished)fail("reco capture census already sealed");
    m_finished=true;
    m_complete=!m_invalid && containerAvailable && m_seen.size()==expectedEntries &&
        m_unknown==0 && m_written==m_eligible;
  }
  bool complete() const {return m_finished && m_complete && !m_invalid;}
  std::size_t inspected() const {return m_seen.size();}
  std::size_t eligible() const {return m_eligible.size();}
  std::size_t written() const {return m_written.size();}
  std::size_t unclassifiable() const {return m_unknown;}
  template<class Event> void record(Event& event) const
  {
    if(!m_finished || m_seen.size()>static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
      throw std::invalid_argument("unsealed or oversized reco capture census");
    event.reco_photon_capture_version=1;
    event.reco_photon_capture_domain=1; // PhotonClusterv1, finite pt/eta/phi, 5<=pt<40, |eta|<0.7
    event.reco_photon_capture_state=complete()?1:2;
    event.reco_photon_inspected=static_cast<std::int32_t>(inspected());
    event.reco_photon_eligible=static_cast<std::int32_t>(eligible());
    event.reco_photon_written=static_cast<std::int32_t>(written());
    event.reco_photon_unclassifiable=static_cast<std::int32_t>(unclassifiable());
  }
 private:
  [[noreturn]] void fail(const char* message)
  {m_invalid=true;throw std::invalid_argument(message);}
  std::map<std::uint32_t,Disposition> m_seen;
  std::unordered_set<std::uint32_t> m_eligible,m_written;
  std::size_t m_unknown=0;
  bool m_finished=false,m_complete=false,m_invalid=false;
};

enum class Outcome : std::int32_t
{
  UNKNOWN=-1, MATCH=0, RECO_FAKE=1, TRUTH_MISS=2
};

enum class Reason
{
  DOMINANT_PRIMARY_IDENTITY,
  EVALUATOR_UNAVAILABLE_OR_INVALID_PRIMARY,
  EVALUATOR_PROVED_NO_PRIMARY,
  DOMINANT_PRIMARY_NOT_PHOTON,
  DOMINANT_PHOTON_ABSENT_FROM_TRUTH_INVENTORY,
  DOMINANT_PHOTON_OUTSIDE_TRUTH_SELECTION,
  NO_SELECTED_RECO_ASSOCIATION,
  RECO_ASSOCIATION_OR_CAPTURE_INCOMPLETE
};

template<class Identity>
struct Link
{
  Identity reco_id{},truth_id{},association_truth_id{};
  Outcome outcome=Outcome::UNKNOWN;
  Reason reason=Reason::EVALUATOR_UNAVAILABLE_OR_INVALID_PRIMARY;
};

template<class Identity>
struct Result
{
  std::vector<Link<Identity>> links;
  bool associations_complete=true;
  bool response_complete=false;
};

// Candidate rows expose id,event_id,dominant_truth. Truth rows expose
// id,event_id,native_track_id,embedding_id,g4_photon_valid. Identity supports
// == and a caller-supplied hash; its default value is the null reference.
//
// recoCaptureComplete is an EXTERNAL assertion that every selected reco object
// in this event was captured, including an explicitly known empty population.
// nodes_ready, terminal_status==0 and truth_denominator_complete do not prove
// it. Likewise truthCaptureComplete must cover the declared truth population.
// Known matches/fakes remain useful without those assertions, but a response
// cannot be certified. Any unresolved selected reco association prevents all
// unmatched selected truth objects in this event from becoming ordinary misses.
template<class Identity,class IdentityHash,class CandidateRows,class TruthRows>
Result<Identity> associate(
    const std::string& sourceScope,const Identity& eventId,
    const CandidateRows& candidates,const TruthRows& rawTruth,
    const std::vector<Identity>& selectedTruthIds,
    bool recoCaptureComplete=false,bool truthCaptureComplete=false)
{
  if(sourceScope.empty() || eventId==Identity{})
    throw std::invalid_argument("photon association needs source/event scope");
  using NativeKey=std::pair<std::int32_t,std::int32_t>;
  std::map<NativeKey,Identity> truthByNativeKey;
  std::unordered_set<Identity,IdentityHash> truthIds,recoIds,selected,matched;
  for(const auto& row:rawTruth)
  {
    if(!(row.event_id==eventId) || row.id==Identity{} ||
       row.native_track_id<=0 || row.g4_photon_valid!=1 ||
       !truthIds.insert(row.id).second ||
       !truthByNativeKey.emplace(NativeKey{row.embedding_id,row.native_track_id},row.id).second)
      throw std::invalid_argument("invalid/ambiguous scoped raw truth inventory");
  }
  for(const auto& id:selectedTruthIds)
    if(!truthIds.count(id) || !selected.insert(id).second)
      throw std::invalid_argument("selected truth must be a unique raw-inventory subset");

  Result<Identity> result;
  result.links.reserve(candidates.size()+selectedTruthIds.size());
  for(const auto& row:candidates)
  {
    if(!(row.event_id==eventId) || row.id==Identity{} || !recoIds.insert(row.id).second)
      throw std::invalid_argument("invalid/ambiguous scoped reco inventory");
    Link<Identity> link;
    link.reco_id=row.id;
    const auto& witness=row.dominant_truth;
    const bool evaluatorKnown=
        witness.evaluator_mode==RJDominantTruthWitnessV1::TOWERINFO ||
        witness.evaluator_mode==RJDominantTruthWitnessV1::LEGACY_RAWTOWER;
    if(evaluatorKnown && witness.state==RJDominantTruthWitnessV1::NO_PRIMARY)
    {
      link.outcome=Outcome::RECO_FAKE;
      link.reason=Reason::EVALUATOR_PROVED_NO_PRIMARY;
    }
    else if(evaluatorKnown && witness.state==RJDominantTruthWitnessV1::VALID_PRIMARY &&
            witness.track_id>0 && witness.pid!=0 &&
            std::isfinite(witness.energy_contribution) && witness.energy_contribution>=0.)
    {
      if(witness.pid!=22)
      {
        link.outcome=Outcome::RECO_FAKE;
        link.reason=Reason::DOMINANT_PRIMARY_NOT_PHOTON;
      }
      else
      {
        const auto found=truthByNativeKey.find({witness.embedding_id,witness.track_id});
        if(found==truthByNativeKey.end())
          link.reason=Reason::DOMINANT_PHOTON_ABSENT_FROM_TRUTH_INVENTORY;
        else
        {
          link.association_truth_id=found->second;
          if(selected.count(found->second))
          {
            link.truth_id=found->second;
            link.outcome=Outcome::MATCH;
            link.reason=Reason::DOMINANT_PRIMARY_IDENTITY;
            matched.insert(found->second);
          }
          else
          {
            link.outcome=Outcome::RECO_FAKE;
            link.reason=Reason::DOMINANT_PHOTON_OUTSIDE_TRUTH_SELECTION;
          }
        }
      }
    }
    if(link.outcome==Outcome::UNKNOWN)result.associations_complete=false;
    result.links.push_back(link);
  }
  // Follow caller order, not unordered-container iteration order.
  for(const auto& id:selectedTruthIds)
  {
    if(matched.count(id))continue;
    Link<Identity> link;
    link.truth_id=id;
    const bool known=recoCaptureComplete && result.associations_complete;
    link.outcome=known?Outcome::TRUTH_MISS:Outcome::UNKNOWN;
    link.reason=known?Reason::NO_SELECTED_RECO_ASSOCIATION:
        Reason::RECO_ASSOCIATION_OR_CAPTURE_INCOMPLETE;
    result.links.push_back(link);
  }
  result.response_complete=recoCaptureComplete && truthCaptureComplete &&
      result.associations_complete;
  return result;
}
}  // namespace RJPhotonAssociationV1

#endif
