// Exact, self-contained reproduction of the event selection used by
// sPHENIX-Collaboration/analysis sEPD-Study TriggerQA at:
//   Fun4All_TriggerQA.C blob d132e704544300897568017a57e0756b8a411ad8
//   TriggerQA.cc       blob 4f5a77f28f61c6c7e42faf6c280a06e686d57789
//
// The pre-skimmer counter is diagnostic only.  The post-skimmer counters
// intentionally preserve Apurva's order and logic.

#include <calostatusskimmer/CaloStatusSkimmer.h>
#include <calotrigger/TriggerAnalyzer.h>
#include <calotrigger/TriggerRunInfoReco.h>
#include <ffaobjects/EventHeader.h>
#include <ffamodules/CDBInterface.h>
#include <ffamodules/FlagHandler.h>
#include <fun4all/Fun4AllBase.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>
#include <fun4all/SubsysReco.h>
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/GlobalVertexReco.h>
#include <mbd/MbdReco.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/recoConsts.h>

#include <TSystem.h>
#include <TH1F.h>

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libfun4allutils.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libglobalvertex.so)
R__LOAD_LIBRARY(libcalotrigger.so)
R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(libCaloStatusSkimmer.so)

namespace the126
{
  class InputCounter : public SubsysReco
  {
   public:
    InputCounter()
      : SubsysReco("THE126InputCounter")
    {
    }

    int process_event(PHCompositeNode*) override
    {
      ++m_events;
      return Fun4AllReturnCodes::EVENT_OK;
    }

    std::uint64_t events() const { return m_events; }

   private:
    std::uint64_t m_events{0};
  };

  class SelectionCounter : public SubsysReco
  {
   public:
    SelectionCounter(InputCounter* input_counter,
                     int run,
                     int chunk,
                     int input_files,
                     const std::string& input_list,
                     const std::string& output_tsv)
      : SubsysReco("THE126SelectionCounter")
      , m_input_counter(input_counter)
      , m_run(run)
      , m_chunk(chunk)
      , m_input_files(input_files)
      , m_input_list(input_list)
      , m_output_tsv(output_tsv)
    {
    }

    int Init(PHCompositeNode*) override
    {
      m_trigger_analyzer = std::make_unique<TriggerAnalyzer>();
      m_event_hist = new TH1F("hEvent", ";Type;Events", 8, 0, 8);
      const char* labels[8] = {
          "All",
          "Has Z",
          "|z| < 10 cm",
          "Trig 12",
          "Trig 14",
          "|z| < 10 cm & (Trig 12 | Trig 14)",
          "|z| < 10 cm & Trig 12",
          "|z| < 10 cm & Trig 14"};
      for (int i = 0; i < 8; ++i)
      {
        m_event_hist->GetXaxis()->SetBinLabel(i + 1, labels[i]);
      }
      Fun4AllServer::instance()->registerHisto(m_event_hist);
      return Fun4AllReturnCodes::EVENT_OK;
    }

    int process_event(PHCompositeNode* topNode) override
    {
      EventHeader* event_header =
          findNode::getClass<EventHeader>(topNode, "EventHeader");
      if (!event_header)
      {
        std::cerr << "[THE126][ERROR] Missing EventHeader" << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }

      ++m_after_calo_status;
      m_event_hist->Fill(0);
      m_trigger_analyzer->decodeTriggers(topNode);
      const bool trig12 = m_trigger_analyzer->didTriggerFire(12);
      const bool trig14 = m_trigger_analyzer->didTriggerFire(14);
      if (trig12)
      {
        ++m_trig12;
        m_event_hist->Fill(3);
      }
      if (trig14)
      {
        ++m_trig14;
        m_event_hist->Fill(4);
      }

      GlobalVertexMap* vertex_map =
          findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
      if (!vertex_map)
      {
        std::cerr << "[THE126][ERROR] Missing GlobalVertexMap" << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }
      if (vertex_map->empty())
      {
        return Fun4AllReturnCodes::ABORTEVENT;
      }

      ++m_has_z;
      m_event_hist->Fill(1);
      GlobalVertex* vertex = vertex_map->begin()->second;
      const double z = vertex->get_z();
      if (std::abs(z) < 10.0)
      {
        ++m_z10;
        m_event_hist->Fill(2);
        if (trig12 || trig14)
        {
          ++m_z10_trig12_or_trig14;
          m_event_hist->Fill(5);
        }
        if (trig12)
        {
          ++m_z10_trig12;
          m_event_hist->Fill(6);
        }
        if (trig14)
        {
          ++m_z10_trig14;
          m_event_hist->Fill(7);
        }
      }

      return Fun4AllReturnCodes::EVENT_OK;
    }

    int End(PHCompositeNode*) override
    {
      std::ofstream out(m_output_tsv);
      if (!out)
      {
        std::cerr << "[THE126][ERROR] Cannot write " << m_output_tsv
                  << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }

      constexpr double sigma_mbd_b = 6.324;
      const double lumi_ub_inv =
          static_cast<double>(m_z10_trig12_or_trig14) /
          (sigma_mbd_b * 1.0e6);
      const double lumi_trig12_ub_inv =
          static_cast<double>(m_z10_trig12) / (sigma_mbd_b * 1.0e6);
      const double lumi_trig14_ub_inv =
          static_cast<double>(m_z10_trig14) / (sigma_mbd_b * 1.0e6);

      out << "run\tchunk\tinput_files\tinput_events\tpass_calo_status"
             "\thas_z\tz10\ttrig12\ttrig14\tz10_trig12_or_trig14"
             "\tz10_trig12\tz10_trig14\tsigma_mbd_b\tlumi_ub_inv"
             "\tlumi_trig12_ub_inv\tlumi_trig14_ub_inv"
             "\tinput_list\n";
      out << m_run << '\t'
          << m_chunk << '\t'
          << m_input_files << '\t'
          << (m_input_counter ? m_input_counter->events() : 0) << '\t'
          << m_after_calo_status << '\t'
          << m_has_z << '\t'
          << m_z10 << '\t'
          << m_trig12 << '\t'
          << m_trig14 << '\t'
          << m_z10_trig12_or_trig14 << '\t'
          << m_z10_trig12 << '\t'
          << m_z10_trig14 << '\t'
          << std::fixed << std::setprecision(6) << sigma_mbd_b << '\t'
          << std::setprecision(10) << lumi_ub_inv << '\t'
          << lumi_trig12_ub_inv << '\t'
          << lumi_trig14_ub_inv << '\t'
          << m_input_list << '\n';
      out.close();

      std::cout << "[THE126][SUMMARY] run=" << m_run
                << " chunk=" << m_chunk
                << " input_events="
                << (m_input_counter ? m_input_counter->events() : 0)
                << " pass_calo_status=" << m_after_calo_status
                << " z10_trig12_or_trig14=" << m_z10_trig12_or_trig14
                << " z10_trig12=" << m_z10_trig12
                << " z10_trig14=" << m_z10_trig14
                << " lumi_ub_inv=" << std::setprecision(10) << lumi_ub_inv
                << std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }

   private:
    InputCounter* m_input_counter{nullptr};
    int m_run{0};
    int m_chunk{0};
    int m_input_files{0};
    std::string m_input_list;
    std::string m_output_tsv;
    std::unique_ptr<TriggerAnalyzer> m_trigger_analyzer;
    TH1F* m_event_hist{nullptr};
    std::uint64_t m_after_calo_status{0};
    std::uint64_t m_has_z{0};
    std::uint64_t m_z10{0};
    std::uint64_t m_trig12{0};
    std::uint64_t m_trig14{0};
    std::uint64_t m_z10_trig12_or_trig14{0};
    std::uint64_t m_z10_trig12{0};
    std::uint64_t m_z10_trig14{0};
  };
}

void Fun4All_AuAuLumi10(const std::string& input_list,
                        const std::string& output_root,
                        const std::string& output_tsv,
                        int run,
                        int chunk,
                        int n_events = 0,
                        const std::string& dbtag = "newcdbtag")
{
  if (run <= 0)
  {
    throw std::runtime_error("run must be positive");
  }

  std::ifstream list_stream(input_list);
  if (!list_stream)
  {
    throw std::runtime_error("cannot open input list " + input_list);
  }
  int input_files = 0;
  std::string line;
  while (std::getline(list_stream, line))
  {
    if (!line.empty()) ++input_files;
  }
  if (input_files == 0)
  {
    throw std::runtime_error("input list is empty: " + input_list);
  }

  std::cout << "[THE126] input_list=" << input_list
            << " output_root=" << output_root
            << " output_tsv=" << output_tsv
            << " run=" << run
            << " chunk=" << chunk
            << " input_files=" << input_files
            << " n_events=" << n_events
            << " dbtag=" << dbtag << std::endl;

  Fun4AllServer* server = Fun4AllServer::instance();
  recoConsts* rc = recoConsts::instance();
  rc->set_StringFlag("CDB_GLOBALTAG", dbtag);
  rc->set_uint64Flag("TIMESTAMP", static_cast<std::uint64_t>(run));
  CDBInterface::instance()->Verbosity(1);

  auto* flag_handler = new FlagHandler();
  server->registerSubsystem(flag_handler);

  auto* input_counter = new the126::InputCounter();
  server->registerSubsystem(input_counter);

  auto* calo_status = new CaloStatusSkimmer("CaloStatusSkimmer");
  server->registerSubsystem(calo_status);

  auto* mbd_reco = new MbdReco();
  server->registerSubsystem(mbd_reco);

  auto* vertex_reco = new GlobalVertexReco();
  server->registerSubsystem(vertex_reco);

  auto* trigger_info = new TriggerRunInfoReco();
  trigger_info->Verbosity(1);
  server->registerSubsystem(trigger_info);

  auto* selection_counter =
      new the126::SelectionCounter(input_counter, run, chunk, input_files,
                                   input_list, output_tsv);
  server->registerSubsystem(selection_counter);

  auto* input = new Fun4AllDstInputManager("THE126Input");
  input->AddListFile(input_list);
  server->registerInputManager(input);

  server->Verbosity(Fun4AllBase::VERBOSITY_QUIET);
  server->run(n_events);
  server->End();
  server->dumpHistos(output_root);
  CDBInterface::instance()->Print();
  server->PrintTimer();
  delete server;

  std::cout << "[THE126] All done" << std::endl;
  gSystem->Exit(0);
}
