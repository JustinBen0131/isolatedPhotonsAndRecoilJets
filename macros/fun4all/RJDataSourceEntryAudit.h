#ifndef RJ_DATA_SOURCE_ENTRY_AUDIT_H
#define RJ_DATA_SOURCE_ENTRY_AUDIT_H
#include "RJDataSourceEntryContract.h"
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/InputFileHandlerReturnCodes.h>
#include <fun4all/SubsysReco.h>
#include <phool/PHNodeIOManager.h>
#include <phool/getClass.h>
#include <ffaobjects/EventHeader.h>
#include <calostatusskimmer/CaloStatusSkimmer.h>
#include <TFile.h>
#include <TNamed.h>
#include <TTree.h>
#include <memory>
#include <string>

namespace rj_source_entry
{
// The PHNodeIOManager cursor identifies the original ROOT entry, not the
// number of events that survived reconstruction or detector-quality cuts.
class InputManager : public Fun4AllDstInputManager
{
 public:
  explicit InputManager(const std::string& name) : Fun4AllDstInputManager(name) {}
  void openForCursor()
  {
    // HasSyncObject() is zero before fileopen; Fun4AllServer::skip can then
    // return success through the base no-op SkipForThisManager implementation.
    // InputFileHandler uses SUCCESS=1, unlike Fun4All EVENT_OK=0.
    if (!IsOpen() && OpenNextFile()!=InputFileHandlerReturnCodes::SUCCESS)
      throw std::runtime_error("could not open exact paired DATA input before skip");
    if (!IManager()) throw std::runtime_error("opened DATA input has no cursor");
  }
  std::int64_t currentEntry()
  {
    if (!IManager()) throw std::runtime_error("original DST input manager is absent");
    return static_cast<std::int64_t>(IManager()->getEventNumber())-1;
  }
};
inline std::shared_ptr<Accounting>& active()
{ static std::shared_ptr<Accounting> value; return value; }
class Tracker : public SubsysReco
{
 public:
  Tracker() : SubsysReco("RJOriginalSourceEntryTracker") {}
  int process_event(PHCompositeNode* top) override
  {
    auto* server=Fun4AllServer::instance();
    auto* jet=dynamic_cast<InputManager*>(server->getInputManager("DST_JET_IN"));
    auto* calo=dynamic_cast<InputManager*>(server->getInputManager("DST_JETCALO_IN"));
    auto* header=findNode::getClass<EventHeader>(top,"EventHeader");
    if (!active() || !jet || !calo || !header)
      throw std::runtime_error("original paired DATA identity unavailable");
    active()->input(jet->currentEntry(),calo->currentEntry(),header->get_RunNumber(),header->get_EvtSequence());
    return Fun4AllReturnCodes::EVENT_OK;
  }
};
class Skimmer : public CaloStatusSkimmer
{
 public:
  explicit Skimmer(const std::string& name) : CaloStatusSkimmer(name) {}
  int process_event(PHCompositeNode* top) override
  {
    const int result=CaloStatusSkimmer::process_event(top);
    if (active())
    {
      if (result!=Fun4AllReturnCodes::EVENT_OK && result!=Fun4AllReturnCodes::ABORTEVENT)
        throw std::runtime_error("unexpected CaloStatusSkimmer return code");
      active()->skim(result==Fun4AllReturnCodes::ABORTEVENT);
    }
    return result; // Preserve the release skimmer's selection and defaults exactly.
  }
};
// RecoilJets' scope-exit writer retains even analysis-rejected events. Fun4All
// skips downstream process_event after ABORTEVENT, but calls ResetEvent for
// every subsystem. Commit once there; never treat an analysis cut as data loss.
class Retained : public SubsysReco
{
 public:
  Retained() : SubsysReco("RJOriginalSourceEntryRetained") {}
  int ResetEvent(PHCompositeNode*) override
  {
    if (!active())
      throw std::runtime_error("retained original DATA identity unavailable");
    active()->finishEvent();
    return Fun4AllReturnCodes::EVENT_OK;
  }
};
inline void writeMetadata(const std::string& path)
{
  if (!active()) return;
  active()->complete();
  TFile file(path.c_str(),"UPDATE");
  auto* directory=file.GetDirectory("ReplayFoundationV1");
  if (file.IsZombie() || !directory) throw std::runtime_error("source audit ROOT output unavailable");
  directory->cd();
  if (directory->Get("RJUpstreamRejectedEventV1"))
    throw std::runtime_error("source audit already present");
  TNamed("source_entry_contract","ORIGINAL_PAIRED_DST_CURSOR_V1").Write();
  TNamed("input_events",std::to_string(active()->attempted()).c_str()).Write();
  TNamed("upstream_rejected_events",std::to_string(active()->rejected().size()).c_str()).Write();
  TNamed("upstream_rejection_module","CaloStatusSkimmer").Write();
  TNamed("upstream_rejection_return","ABORTEVENT").Write();
  TTree rejected("RJUpstreamRejectedEventV1","Original input identities rejected before RecoilJets");
  Long64_t entry=-1,physical=-1; Int_t run=-1;
  rejected.Branch("source_entry_ordinal",&entry);
  rejected.Branch("physical_event_sequence",&physical);
  rejected.Branch("run",&run);
  for (const auto& row : active()->rejected())
  { entry=row.entry; physical=row.physical; run=row.run; rejected.Fill(); }
  rejected.Write();
  rejected.SetDirectory(nullptr);
  if (file.TestBit(TFile::kWriteError)) throw std::runtime_error("source audit write failed");
  file.Close();
}
}
#endif
