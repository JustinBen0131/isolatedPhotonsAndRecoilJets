//======================================================================
//  Fun4All_recoilJets.C
//  --------------------------------------------------------------------
#pragma once
#if defined(__CINT__) || defined(__CLING__)
  R__ADD_INCLUDE_PATH(/sphenix/u/patsfan753/thesisAnalysis/install/include)
#endif
#if defined(__CLING__)
  #pragma cling add_include_path("/sphenix/u/patsfan753/thesisAnalysis/install/include")
#endif
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

//–––– Standard Fun4All / sPHENIX ––––––––––––––––––––––––––––––––––––––
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllUtils.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <calobase/TowerInfoContainer.h>
#include <mbd/MbdPmtContainer.h>
#include <caloreco/CaloTowerStatus.h>
#include <caloreco/CaloWaveformProcessing.h>
#include <caloreco/CaloTowerBuilder.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <calotrigger/MinimumBiasClassifier.h>
#include <ffamodules/FlagHandler.h>
#include <ffamodules/CDBInterface.h>
#include <clusteriso/ClusterIso.h>
#include <calotrigger/TriggerRunInfoReco.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <caloreco/CaloGeomMapping.h>
#include <caloreco/RawClusterPositionCorrection.h>
#include <calobase/RawTowerGeom.h>
#include <caloreco/RawTowerCalibration.h>
#include <caloreco/PhotonClusterBuilder.h>
#include <jetbase/Jet.h>
#include <g4jets/TruthJetInput.h>

#include <jetbase/FastJetOptions.h>
#include <jetbackground/FastJetAlgoSub.h>
#include <globalvertex/GlobalVertexReco.h>
#include <caloreco/CaloTowerCalib.h>

#include <centrality/CentralityReco.h>
#include <mbd/MbdEvent.h>
#include <mbd/MbdReco.h>
#include <zdcinfo/ZdcReco.h>
#include <phool/recoConsts.h>
#include <phool/PHRandomSeed.h>
#include <jetbase/FastJetOptions.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <jetbase/JetCalib.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/DetermineTowerBackground.h>
#include <jetbackground/SubtractTowers.h>
#include <jetbackground/CopyAndSubtractJets.h>
#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/src/RecoilJets.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <cstdlib>     // getenv
#include <algorithm>   // std::transform
#include <cctype>      // std::tolower
#include <iomanip>     // std::setw, std::setprecision
#include <dlfcn.h>     // dlopen RTLD_NOLOAD
#include <TSystem.h>   // gSystem, GetBuildArch/Compiler info
#include <Calo_Calib.C>

// Load your local overrides FIRST so their symbols/dictionaries win
R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libcalo_reco.so)
R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libcalo_io.so)
R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libRecoilJets.so)
R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libclusteriso.so)
R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libjetbase.so)


// Then load the rest of the environment stack
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffarawobjects.so)
R__LOAD_LIBRARY(libcaloTreeGen.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libjetbase.so)
R__LOAD_LIBRARY(libglobalvertex.so)
R__LOAD_LIBRARY(libcentrality.so)      // always
R__LOAD_LIBRARY(libcentrality_io.so)   // if you instantiate CentralityReco
R__LOAD_LIBRARY(libcalotrigger.so)
R__LOAD_LIBRARY(libzdcinfo.so)
R__LOAD_LIBRARY(libmbd.so)


//======================================================================
//  Convenience helpers
//======================================================================
namespace detail
{
  /// Throw a nicely formatted exception on unrecoverable error
  [[noreturn]] void bail(const std::string& msg)
  {
    std::ostringstream oss;
    oss << "\n[FATAL] Fun4All_recoilJets :: " << msg << '\n';
    throw std::runtime_error(oss.str());
  }

  /// Trim whitespace from both ends (for robust list-file parsing)
  inline std::string trim(std::string s)
  {
    const char* ws = " \t\r\n";
    s.erase(0, s.find_first_not_of(ws));
    s.erase(s.find_last_not_of(ws) + 1);
    return s;
  }
}

namespace detail
{
  inline FastJetAlgoSub* fjAlgo(const float R)
  {
    FastJetOptions o{};               // IMPORTANT: value-initialize ALL fields to safe defaults
    o.algo            = Jet::ANTIKT;  // algorithm
    o.jet_R           = R;            // jet radius
    o.use_jet_min_pt  = true;         // enable a ptmin
    o.jet_min_pt      = 0.0f;         // ptmin value
    o.verbosity       = 0;            // quiet
    return new FastJetAlgoSub(o);
  }
}

// ============================================================================
// Ensure the composite nodes JetCalib hard-requires exist.
// JetCalib::CreateNodeTree requires BOTH:
//   - PHCompositeNode "ANTIKT"
//   - PHCompositeNode "TOWER"
// Your pp DSTs often don't have "TOWER", so JetCalib aborts the run.
// ============================================================================
class EnsureJetCalibNodes final : public SubsysReco
{
 public:
  explicit EnsureJetCalibNodes(const std::string& name="EnsureJetCalibNodes")
    : SubsysReco(name) {}

  int InitRun(PHCompositeNode* topNode) override
  {
    PHNodeIterator iter(topNode);

    auto* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
    if (!dstNode)
    {
      std::cerr << "[EnsureJetCalibNodes] FATAL: DST node not found\n";
      return Fun4AllReturnCodes::ABORTRUN;
    }

    // JetCalib requires ANTIKT to exist
    auto* antiktNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "ANTIKT"));
    if (!antiktNode)
    {
      antiktNode = new PHCompositeNode("ANTIKT");
      dstNode->addNode(antiktNode);
      if (Verbosity() > 0)
        std::cout << "[EnsureJetCalibNodes] created DST/ANTIKT composite node\n";
    }

    // JetCalib requires TOWER to exist (this is what you're missing)
    auto* towerNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "TOWER"));
    if (!towerNode)
    {
      towerNode = new PHCompositeNode("TOWER");
      dstNode->addNode(towerNode);
      if (Verbosity() > 0)
        std::cout << "[EnsureJetCalibNodes] created DST/TOWER composite node\n";
    }

    return Fun4AllReturnCodes::EVENT_OK;
  }
};


class TowerAudit final : public SubsysReco
{
 public:
  explicit TowerAudit(const std::string& name,
                      const std::string& nodeCEMC,
                      const std::string& nodeHCIN,
                      const std::string& nodeHCOUT,
                      int maxPrintEvents = 5)
    : SubsysReco(name)
    , m_nodeCEMC(nodeCEMC)
    , m_nodeHCIN(nodeHCIN)
    , m_nodeHCOUT(nodeHCOUT)
    , m_maxPrint(maxPrintEvents)
  {}

  int InitRun(PHCompositeNode* topNode) override
  {
    // Just verify nodes exist at InitRun (they may still be filled per-event)
    auto* cemc = findNode::getClass<TowerInfoContainer>(topNode, m_nodeCEMC);
    auto* hcin = findNode::getClass<TowerInfoContainer>(topNode, m_nodeHCIN);
    auto* hcot = findNode::getClass<TowerInfoContainer>(topNode, m_nodeHCOUT);

    if (Verbosity() > 0)
    {
      std::cout << "[" << Name() << "::InitRun]"
                << " nodes:"
                << " CEMC=" << m_nodeCEMC << (cemc ? "(OK)" : "(MISSING)")
                << " HCALIN=" << m_nodeHCIN << (hcin ? "(OK)" : "(MISSING)")
                << " HCALOUT=" << m_nodeHCOUT << (hcot ? "(OK)" : "(MISSING)")
                << std::endl;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }

  int process_event(PHCompositeNode* topNode) override
  {
    ++m_evt;
    if (m_evt > m_maxPrint) return Fun4AllReturnCodes::EVENT_OK;

    auto auditOne = [&](const char* label, const std::string& node)
    {
      auto* cont = findNode::getClass<TowerInfoContainer>(topNode, node);
      if (!cont)
      {
        std::cout << "[" << Name() << "] evt=" << m_evt
                  << " " << label << " node=" << node << " MISSING\n";
        return;
      }

      const unsigned int nt = cont->size();
      unsigned int nNull = 0, nGood = 0, nBad = 0;
      unsigned int nEpos = 0, nEposGood = 0, nEposBad = 0;

      double sumE = 0.0, sumEgood = 0.0, sumEbad = 0.0;
      double sumT = 0.0, sumTgood = 0.0, sumTbad = 0.0;
      unsigned int nT = 0, nTgood = 0, nTbad = 0;

      // sample the whole container (24576 for CEMC is fine for a few events)
      for (unsigned int ch = 0; ch < nt; ++ch)
      {
        TowerInfo* t = cont->get_tower_at_channel(ch);
        if (!t) { ++nNull; continue; }

        const bool good = t->get_isGood();
        const float e   = t->get_energy();
        const float tim = t->get_time();

        if (good) ++nGood; else ++nBad;

        if (std::isfinite(e))
        {
          sumE += e;
          if (e > 0) ++nEpos;
          if (good) { sumEgood += e; if (e > 0) ++nEposGood; }
          else      { sumEbad  += e; if (e > 0) ++nEposBad;  }
        }

        if (std::isfinite(tim))
        {
          ++nT;
          sumT += tim;
          if (good) { ++nTgood; sumTgood += tim; }
          else      { ++nTbad;  sumTbad  += tim; }
        }
      }

      const double fracBad = (nt > 0) ? (double)nBad / (double)nt : 0.0;
      const double meanT   = (nT > 0) ? sumT / (double)nT : 0.0;
      const double meanTg  = (nTgood > 0) ? sumTgood / (double)nTgood : 0.0;
      const double meanTb  = (nTbad > 0) ? sumTbad / (double)nTbad : 0.0;

      std::cout << "[" << Name() << "] evt=" << m_evt << " " << label
                << " node=" << node
                << " ntowers=" << nt
                << " null=" << nNull
                << " good=" << nGood
                << " bad=" << nBad
                << " fracBad=" << std::fixed << std::setprecision(4) << fracBad
                << " | sumE=" << std::setprecision(3) << sumE
                << " (good=" << sumEgood << ", bad=" << sumEbad << ")"
                << " | E>0: all=" << nEpos << " good=" << nEposGood << " bad=" << nEposBad
                << " | meanTime=" << std::setprecision(3) << meanT
                << " (good=" << meanTg << ", bad=" << meanTb << ")"
                << "\n";
    };

    auditOne("CEMC",   m_nodeCEMC);
    auditOne("HCALIN", m_nodeHCIN);
    auditOne("HCALOUT",m_nodeHCOUT);

    return Fun4AllReturnCodes::EVENT_OK;
  }

 private:
  std::string m_nodeCEMC;
  std::string m_nodeHCIN;
  std::string m_nodeHCOUT;
  int m_maxPrint = 5;
  int m_evt = 0;
};



class JetCalibOneEventProbe final : public SubsysReco
{
 public:
  JetCalibOneEventProbe(const std::string& name,
                        const std::string& rawNode,
                        const std::string& calibNode,
                        int maxJetsToPrint = 12)
    : SubsysReco(name)
    , m_rawNode(rawNode)
    , m_calibNode(calibNode)
    , m_maxJets(maxJetsToPrint)
  {}

  int process_event(PHCompositeNode* topNode) override
  {
    // Hard gate: ONLY show anything at ultra-verbose
    if (Verbosity() < 20) return Fun4AllReturnCodes::EVENT_OK;

    // Print only once for the whole job
    if (m_fired) return Fun4AllReturnCodes::EVENT_OK;

    auto* raw   = findNode::getClass<JetContainer>(topNode, m_rawNode);
    auto* calib = findNode::getClass<JetContainer>(topNode, m_calibNode);

    // Only trigger when both exist and actually have jets
    if (!raw || !calib) return Fun4AllReturnCodes::EVENT_OK;
    if (raw->size() == 0 || calib->size() == 0) return Fun4AllReturnCodes::EVENT_OK;

    // Optional: grab z-vertex JetCalib uses (GlobalVertexMap/MBD z)
    float zvtx = -999.0f;
    auto* vtxmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
    if (vtxmap && !vtxmap->empty())
    {
      GlobalVertex* vtx = vtxmap->begin()->second;
      if (vtx)
      {
        auto mbdStart = vtx->find_vertexes(GlobalVertex::MBD);
        auto mbdEnd   = vtx->end_vertexes();
        for (auto it = mbdStart; it != mbdEnd; ++it)
        {
          const auto& [type, vec] = *it;
          if (type != GlobalVertex::MBD) continue;
          for (const auto* vv : vec)
          {
            if (!vv) continue;
            zvtx = vv->get_z();
          }
        }
      }
    }

    const int nRaw   = (int) raw->size();
    const int nCalib = (int) calib->size();
    const int nPrint = std::min({m_maxJets, nRaw, nCalib});

    std::cout << "\n====================================================================\n";
    std::cout << "[JES PROBE] One-event before/after JetCalib\n";
    std::cout << "  rawNode   = " << m_rawNode   << " (n=" << nRaw   << ")\n";
    std::cout << "  calibNode = " << m_calibNode << " (n=" << nCalib << ")\n";
    std::cout << "  zvtx(MBD) = " << std::fixed << std::setprecision(2) << zvtx << " cm\n";
    std::cout << "--------------------------------------------------------------------\n";
    std::cout << "  i |   pt_raw   eta_raw   phi_raw  ||   pt_cal   eta_cal   phi_cal  ||  scale\n";
    std::cout << "----+---------------------------------------------------------------+--------\n";

    for (int i = 0; i < nPrint; ++i)
    {
      const Jet* jr = raw->get_jet(i);
      const Jet* jc = calib->get_jet(i);
      if (!jr || !jc) continue;

      const float pr = jr->get_pt();
      const float pc = jc->get_pt();
      const float sc = (pr > 0.0f) ? (pc / pr) : -1.0f;

      std::cout
        << std::setw(3) << i << " | "
        << std::setw(8) << std::fixed << std::setprecision(3) << pr << " "
        << std::setw(8) << std::setprecision(3) << jr->get_eta() << " "
        << std::setw(8) << std::setprecision(3) << jr->get_phi() << " || "
        << std::setw(8) << std::setprecision(3) << pc << " "
        << std::setw(8) << std::setprecision(3) << jc->get_eta() << " "
        << std::setw(8) << std::setprecision(3) << jc->get_phi() << " || "
        << std::setw(6) << std::setprecision(3) << sc
        << "\n";
    }

    std::cout << "====================================================================\n\n";

    m_fired = true;
    return Fun4AllReturnCodes::EVENT_OK;
  }

 private:
  std::string m_rawNode;
  std::string m_calibNode;
  int m_maxJets = 12;
  bool m_fired = false;
};



//======================================================================
//  The actual steering macro
//======================================================================
void Fun4All_recoilJets(const int   nEvents   =  0,
                        const char* listFile  = "input_files.list",
                        const char* outRoot   = "TrigPlot.root",
                        const bool  verbose   = false)
{
  //--------------------------------------------------------------------
  // 0.  Banner & basic environment sanity
  //--------------------------------------------------------------------
  if (verbose) {
      std::cout << "\n>>> Fun4All_recoilJets – ana.495 driver <<<\n"
                << "    Input list : " << listFile  << '\n'
                << "    Output file: " << outRoot   << '\n'
                << "    nEvents    : " << nEvents   << (nEvents==0? " (all)\n":"\n");
  }

  Fun4AllServer* se = Fun4AllServer::instance();
  if (!se) detail::bail("unable to obtain Fun4AllServer instance!");

  //--------------------------------------------------------------------
  // 1.  Parse the file list & determine run / segment
  //--------------------------------------------------------------------
  std::ifstream list(listFile);
  if (!list.is_open())
      detail::bail("cannot open input list \"" + std::string(listFile) + "\"");

  // ---------------------------------------------------------------
    // Parse list file
    //   DATA:  1 column  -> calo DST
    //   SIM :  5 columns -> calo + G4Hits + (truth jets) + global + mbd_epd
    // ---------------------------------------------------------------
    std::vector<std::string> filesCalo;
    std::vector<std::string> filesG4;
    std::vector<std::string> filesJets;
    std::vector<std::string> filesGlobal;
    std::vector<std::string> filesMbd;

    for (std::string line; std::getline(list, line); )
    {
          line = detail::trim(line);
          if (line.empty()) continue;
          if (!line.empty() && line[0] == '#') continue;

          // Columns:
          //   DATA : <DST_CALO_CLUSTER>
          //   SIM  : <DST_CALO_CLUSTER> <G4Hits> <DST_JETS> <DST_GLOBAL> <DST_MBD_EPD>
          std::istringstream iss(line);
          std::string fCalo, fG4, fJets, fGlobal, fMbd;
          iss >> fCalo >> fG4 >> fJets >> fGlobal >> fMbd;

          if (fCalo.empty()) continue;

          filesCalo.emplace_back(fCalo);

          // Keep vectors key-aligned (same length as filesCalo):
          // empty string == "not provided on this line"
          filesG4.emplace_back((!fG4.empty() && fG4 != "NONE") ? fG4 : std::string{});
          filesJets.emplace_back((!fJets.empty() && fJets != "NONE") ? fJets : std::string{});
          filesGlobal.emplace_back((!fGlobal.empty() && fGlobal != "NONE") ? fGlobal : std::string{});
          filesMbd.emplace_back((!fMbd.empty() && fMbd != "NONE") ? fMbd : std::string{});
    }

  if (filesCalo.empty())
        detail::bail("input list \"" + std::string(listFile) + "\" is empty");

  const std::string& firstFile = filesCalo.front(); // used for GetRunSegment

  // Simulation detection:
  //  - primary: RJ_DATASET=isSim
  //  - fallback: RJ_IS_SIM=1 (wrapper sets this)
  bool isSim = false;
  if (const char* ds = std::getenv("RJ_DATASET"))
  {
      std::string s = detail::trim(std::string(ds));
      std::transform(s.begin(), s.end(), s.begin(),
                     [](unsigned char c){ return std::tolower(c); });
      if (s == "issim" || s == "sim") isSim = true;
    }
    if (!isSim)
    {
      if (const char* f = std::getenv("RJ_IS_SIM"))
        isSim = (std::atoi(f) != 0);
  }

  auto runSegPair = Fun4AllUtils::GetRunSegment(firstFile);
  int run = runSegPair.first;
  int seg = runSegPair.second;

  if (run <= 0)
  {
      if (isSim)
      {
        // If the filename doesn't encode a run number, use a safe CDB timestamp.
        run = 1;
        seg = 0;

        // Optional override if you ever want it:
        //   export RJ_CDB_TIMESTAMP=XXXX
        if (const char* ts = std::getenv("RJ_CDB_TIMESTAMP"))
        {
          try
          {
            const int t = std::stoi(ts);
            if (t > 0) run = t;
          }
          catch (...) { /* keep default */ }
        }

        if (verbose)
          std::cout << "[INFO] Simulation mode: run/seg not in filename → using TIMESTAMP=" << run << "\n";
      }
      else
      {
        detail::bail("failed to extract run number from first file: " + firstFile);
      }
  }

  if (verbose)
          std::cout << "[INFO] Run=" << run << "  Seg=" << seg
                    << "  (" << filesCalo.size() << " files)\n";

  // --------------------------------------------------------------------
  // Global verbosity control (RJ_VERBOSITY from env; defaults to 10;
  // Condor detection → 0). Also silences std::cout/cerr globally when 0.
  // --------------------------------------------------------------------
  int vlevel = 10;
  if (const char* venv = std::getenv("RJ_VERBOSITY"))
  {
      try { vlevel = std::stoi(venv); } catch (...) { vlevel = 10; }
  }
  else
  {
      // Detect Condor jobs defensively; default to quiet in batch.
      if (std::getenv("_CONDOR_SCRATCH_DIR") || std::getenv("_CONDOR_JOB_AD"))
        vlevel = 0;
  }

  // RAII silence for global std::cout/cerr if vlevel==0
  struct ScopedSilence {
    std::ofstream   sink;
    std::streambuf* cout_save = nullptr;
    std::streambuf* cerr_save = nullptr;
    bool active = false;
    void enable() {
      if (active) return;
      sink.open("/dev/null");
      cout_save = std::cout.rdbuf(sink.rdbuf());
      cerr_save = std::cerr.rdbuf(sink.rdbuf());
      active = true;
    }
    ~ScopedSilence() {
      if (active) {
        std::cout.rdbuf(cout_save);
        std::cerr.rdbuf(cerr_save);
      }
    }
  } _silence;

  if (vlevel == 0) _silence.enable();

  // --------------------------------------------------------------------
  // 2.  CDB + IO managers
  // --------------------------------------------------------------------
  recoConsts* rc = recoConsts::instance();
  // CDB_GLOBALTAG is REQUIRED for any CDBInterface::getUrl() call.
  // Allow override via env, otherwise default to your known-good tag.
  std::string gtag = "newcdbtag";
  if (const char* envGT = std::getenv("RJ_CDB_GLOBALTAG"))
  {
      std::string tmp = detail::trim(std::string(envGT));
      if (!tmp.empty()) gtag = tmp;
  }
  rc->set_StringFlag("CDB_GLOBALTAG", gtag);

  if (vlevel > 0)
      std::cout << "[INFO] CDB_GLOBALTAG=" << gtag << "\n";

  // TIMESTAMP:
  //  - DATA: use run number (must be > 1000 so Calo_Calib treats it as DATA)
  //  - SIM : use a known-good fixed timestamp (SIM does NOT run Process_Calo_Calib)
  unsigned long long cdbts = static_cast<unsigned long long>(run);

  if (!isSim)
  {
      // DATA must satisfy Calo_Calib's heuristic (TIMESTAMP > 1000 => data)
      if (cdbts <= 1000ULL)
      {
        std::cerr << "[FATAL] DATA run number " << cdbts
                  << " would make Calo_Calib treat this as SIM (TIMESTAMP<=1000).\n";
        throw std::runtime_error("Invalid DATA TIMESTAMP for Calo_Calib (must be > 1000).");
      }
    }
    else
    {
      cdbts = 47289ULL;  // keep your old working SIM timestamp
    }

    if (const char* ts = std::getenv("RJ_CDB_TIMESTAMP"))
    {
      char* end = nullptr;
      unsigned long long tmp = std::strtoull(ts, &end, 10);
      if (end != ts && tmp > 0ULL) cdbts = tmp;
  }

  rc->set_uint64Flag("TIMESTAMP", cdbts);

  if (vlevel > 0)
      std::cout << "[INFO] CDB TIMESTAMP=" << rc->get_uint64Flag("TIMESTAMP")
                << " (isSim=" << (isSim ? "true" : "false") << ")\n";

  CDBInterface::instance()->Verbosity(0);


  auto* flag = new FlagHandler();
  se->registerSubsystem(flag);

  // ------------------------------------------------------------------
  // Decide how to source truth jets in SIM:
  //
  //   RJ_TRUTH_JETS_MODE=AUTO  (default)
  //       - if list has 3rd column (DST_JETS): read jets directly from DST
  //       - otherwise: build truth jets from TRUTH particles (TruthJetInput)
  //
  //   RJ_TRUTH_JETS_MODE=DST
  //       - require 3rd column (DST_JETS) and read truth jets from DST
  //
  //   RJ_TRUTH_JETS_MODE=BUILD
  //       - ignore 3rd column and build truth jets from TRUTH particles
  //
  //   RJ_TRUTH_JETS_MODE=BOTH
  //       - read DST jets AND also build a second "from particles" truth-jet
  //         collection under a different node name for QA comparisons
  // ------------------------------------------------------------------
  bool useDSTTruthJets = false;
  bool buildTruthJetsFromParticles = false;
  bool buildTruthJetsAsAltNode = false;   // only true in BOTH mode

    auto all_nonempty = [](const std::vector<std::string>& v) -> bool
    {
        if (v.empty()) return false;
        for (const auto& s : v) if (s.empty()) return false;
        return true;
      };

      const bool listHasG4     = all_nonempty(filesG4);
      const bool listHasJets   = all_nonempty(filesJets);
      const bool listHasGlobal = all_nonempty(filesGlobal);
      const bool listHasMbd    = all_nonempty(filesMbd);

      // Normalize env token to lowercase
      auto env_lower = [](const char* key, const std::string& def) -> std::string
      {
        const char* raw = std::getenv(key);
        std::string s = raw ? detail::trim(std::string(raw)) : def;
        std::transform(s.begin(), s.end(), s.begin(),
                       [](unsigned char c){ return std::tolower(c); });
        return s;
      };

      const std::string truthMode = env_lower("RJ_TRUTH_JETS_MODE", "auto");
      if (truthMode == "dst")
      {
        useDSTTruthJets = true;
        buildTruthJetsFromParticles = false;
        buildTruthJetsAsAltNode = false;
      }
      else if (truthMode == "build")
      {
        useDSTTruthJets = false;
        buildTruthJetsFromParticles = true;
        buildTruthJetsAsAltNode = false;
      }
      else if (truthMode == "both")
      {
        useDSTTruthJets = true;
        buildTruthJetsFromParticles = true;
        buildTruthJetsAsAltNode = true;
      }
      else
      {
        // AUTO (or any unrecognized token): prefer DST truth jets if provided
        useDSTTruthJets = listHasJets;
        buildTruthJetsFromParticles = !listHasJets;
        buildTruthJetsAsAltNode = false;
    }

    // ------------------ Calo cluster DST (always) -------------------
    auto* inCalo = new Fun4AllDstInputManager("DST_CALO_CLUSTER_IN");
    for (const auto& f : filesCalo) inCalo->AddFile(f);
    se->registerInputManager(inCalo);

    if (isSim)
    {
        // For SIM we REQUIRE reco-vertex streams (GLOBAL + MBD_EPD).
        // G4Hits is OPTIONAL: if missing, we skip truth-photon matching QA.
        if (!listHasGlobal || !listHasMbd)
        {
          std::ostringstream os;
          os << "isSim requires reco-vertex DSTs (DST_GLOBAL + DST_MBD_EPD) paired 1:1 with calo files.\n"
             << "Input list must include these columns per line:\n"
             << "  <DST_CALO_CLUSTER> <G4_or_NONE> <DST_JETS_or_NONE> <DST_GLOBAL> <DST_MBD_EPD>\n"
             << "But your list is missing DST_GLOBAL and/or DST_MBD_EPD on at least one line.";
          detail::bail(os.str());
        }

        // G4 is OPTIONAL (photonjet productions may not provide it).
        // If you want to force it: export RJ_REQUIRE_G4=1
        bool requireG4 = false;
        if (const char* env = std::getenv("RJ_REQUIRE_G4")) requireG4 = (std::atoi(env) != 0);

        if (listHasG4)
        {
          auto* inG4 = new Fun4AllDstInputManager("DST_G4HITS_IN");
          for (const auto& f : filesG4) inG4->AddFile(f);
          se->registerInputManager(inG4);
        }
        else
        {
          if (requireG4)
          {
            detail::bail("RJ_REQUIRE_G4=1 but no G4Hits stream was provided in the input list.");
          }
          else
          {
            std::cout << "\033[31m[WARN] isSim: no G4Hits stream provided (or it's 'NONE'). "
                         "Continuing; truth-photon matching QA will be skipped.\033[0m\n";
          }
        }

        auto* inGlobal = new Fun4AllDstInputManager("DST_GLOBAL_IN");
        for (const auto& f : filesGlobal) inGlobal->AddFile(f);
        se->registerInputManager(inGlobal);

        auto* inMbd = new Fun4AllDstInputManager("DST_MBD_EPD_IN");
        for (const auto& f : filesMbd) inMbd->AddFile(f);
        se->registerInputManager(inMbd);

        if (verbose)
          std::cout << "[INFO] isSim: registered input managers (Calo + Global + MBD_EPD"
                    << (listHasG4 ? " + G4" : " (no G4)") << ")\n";
    }


      // ------------------ Jets DST (SIM optional; for truth jets) -------
      if (isSim && useDSTTruthJets)
      {
        if (!listHasJets)
        {
          std::ostringstream os;
          os << "RJ_TRUTH_JETS_MODE=" << truthMode << " requires a list with at least 3 columns:\n"
             << "  <DST_CALO_CLUSTER> <G4Hits> <DST_JETS> [<DST_GLOBAL> <DST_MBD_EPD>]\n"
             << "Use your staged 5-column master list built from the matched lists.";
          detail::bail(os.str());
        }

        auto* inJets = new Fun4AllDstInputManager("DST_JETS_IN");
        for (const auto& f : filesJets) inJets->AddFile(f);
        se->registerInputManager(inJets);

        if (verbose)
          std::cout << "[INFO] isSim: registered DST_JETS input manager (truth jets from DST)\n";
      }

      if (verbose && isSim)
      {
        std::cout << "[INFO] Truth-jet mode (RJ_TRUTH_JETS_MODE=" << truthMode << "): "
                  << (useDSTTruthJets ? "DST" : "")
                  << ((useDSTTruthJets && buildTruthJetsFromParticles) ? "+" : "")
                  << (buildTruthJetsFromParticles ? "BUILD" : "")
                  << "\n";
    }

  // --------------------------------------------------------------------
  // 3.  Geometry + status + calibration + clustering
  // --------------------------------------------------------------------
  //
  // IMPORTANT:
  // RawClusterBuilderTemplate::process_event() ALWAYS requires "TOWERGEOM_CEMC".
  // But CaloGeomMapping with UseDetailedGeometry(true) publishes ONLY
  // "TOWERGEOM_CEMC_DETAILED" for CEMC.
  // So we ALWAYS create the legacy node (TOWERGEOM_CEMC), and optionally
  // also create the detailed node (TOWERGEOM_CEMC_DETAILED).
  //
  bool useDetailedCemcGeom = true;  // keep your current behavior by default
  if (const char* env = std::getenv("RJ_DETAILED_CEMC_GEOM"))
  {
      useDetailedCemcGeom = (std::atoi(env) != 0);
    }

    // Always publish the legacy/simple CEMC geometry node: TOWERGEOM_CEMC
    {
      auto* geomCemcLegacy = new CaloGeomMapping("Geom_CEMC");
      geomCemcLegacy->set_detector_name("CEMC");
      geomCemcLegacy->set_UseDetailedGeometry(false);
      se->registerSubsystem(geomCemcLegacy);
    }

    // Optionally publish the detailed CEMC geometry node: TOWERGEOM_CEMC_DETAILED
    if (useDetailedCemcGeom)
    {
      auto* geomCemcDetailed = new CaloGeomMapping("Geom_CEMC_DETAILED");
      geomCemcDetailed->set_detector_name("CEMC");
      geomCemcDetailed->set_UseDetailedGeometry(true);
      se->registerSubsystem(geomCemcDetailed);
    }

    // HCAL nodes (detailed not supported; CaloGeomMapping will fall back internally)
    for (const std::string& det : {"HCALIN","HCALOUT"})
    {
      auto* geom = new CaloGeomMapping(("Geom_" + det).c_str());
      geom->set_detector_name(det);
      geom->set_UseDetailedGeometry(true);
      se->registerSubsystem(geom);
  }


    // ------------------------------------------------------------------
    // Calo calibration + clustering
    //
    // DATA: run the standard Calo_Calib reconstruction chain.
    // SIM : skip (SIM analysis DST already has calibrated towers/clusters; re-running
    //             the reco chain is what triggers CaloTowerStatus hotmap/default-map failure).
    // ------------------------------------------------------------------
    if (!isSim)
    {
      if (vlevel > 0) std::cout << "[DATA] running Process_Calo_Calib()\n";
      Process_Calo_Calib();
    }
    else
    {
      if (vlevel > 0)
      {
        std::cout << "[isSim] skipping Process_Calo_Calib() "
                     "(SIM DST already has TOWERINFO_CALIB and CLUSTERINFO_CEMC)\n";
      }
    }


  if (vlevel > 0) std::cout << "Calibrating MBD" << std::endl;
  std::unique_ptr<MbdReco> mbdreco = std::make_unique<MbdReco>();
  se->registerSubsystem(mbdreco.release());

  if (!isSim)
    {
      if (vlevel > 0) std::cout << "Calibrating ZDC" << std::endl;
      auto* zdcreco = new ZdcReco();
      zdcreco->set_zdc1_cut(0.0);
      zdcreco->set_zdc2_cut(0.0);
      se->registerSubsystem(zdcreco);
    }
    else
    {
      if (vlevel > 0) std::cout << "[isSim] skipping ZdcReco (sim DST has no TOWERS_ZDC)" << std::endl;
  }


  if (vlevel > 0) std::cout << "Retrieving Vtx Info" << std::endl;
  std::unique_ptr<GlobalVertexReco> gvertex = std::make_unique<GlobalVertexReco>();
  se->registerSubsystem(gvertex.release());

  // Decide dataset early so we can gate Au+Au-only modules (centrality)
  // IMPORTANT: isSim should behave like pp-style reconstruction here.
  bool isAuAuData = true;

  if (const char* env = std::getenv("RJ_DATASET"))
  {
        std::string s = detail::trim(std::string(env));
        std::string sLower = s;
        std::transform(sLower.begin(), sLower.end(), sLower.begin(),
                       [](unsigned char c){ return std::tolower(c); });

        if (sLower == "issim" || sLower == "sim")
        {
          isSim = true;        // ensure consistency even if RJ_IS_SIM wasn't set
          isAuAuData = false;  // sim uses pp-style chain
        }
        else if (sLower == "ispp" || sLower == "pp")
        {
          isAuAuData = false;
        }
        else if (sLower == "isauau" || sLower == "auau" || sLower == "aa")
        {
          isAuAuData = true;
        }
        else
        {
          // Unknown token: fallback heuristic based on known pp runs
          isAuAuData = (run > 53864);
        }
    }
    else
    {
        // Fallback heuristic based on known pp runs: ≤ 53864 → p+p, > 53864 → Au+Au
        isAuAuData = (run > 53864);
    }

  if (isAuAuData)
  {
      std::cout << "building minbias classifier" << std::endl;
      auto* mb = new MinimumBiasClassifier();
      mb->Verbosity(0);
      mb->setOverwriteScale(
            "/cvmfs/sphenix.sdcc.bnl.gov/calibrations/sphnxpro/cdb/CentralityScale/42/6b/426bc1b56ba544201b0213766bee9478_cdb_centrality_scale_54912.root");
      se->registerSubsystem(mb);
      
      if (vlevel > 0) std::cout << "building centrality classifier (Au+Au)" << std::endl;
      auto* cent = new CentralityReco();
      cent->Verbosity(0);
      cent->setOverwriteScale(
              "/cvmfs/sphenix.sdcc.bnl.gov/calibrations/sphnxpro/cdb/CentralityScale/42/6b/426bc1b56ba544201b0213766bee9478_cdb_centrality_scale_54912.root");
      se->registerSubsystem(cent);
  }
  else
  {
        if (vlevel > 0) std::cout << "[pp dataset] skipping CentralityReco" << std::endl;
  }
    
  // ---------------------- Reco jets + JES calibration (pp-style only) ----------------------
  //
  // IMPORTANT:
  // JetCalib::CreateNodeTree() requires a PHCompositeNode named "TOWER".
  // Many pp DSTs do NOT have it, so JetCalib aborts unless we create it.
  // We register EnsureJetCalibNodes once (only when we intend to run JetCalib).
  //
    // Apply pp JES calibration for ALL pp-style running (pp data AND isSim).
    // Keep Au+Au excluded.
    const bool doJetCalibAny = (!isAuAuData);
    if (doJetCalibAny)
    {
        auto* ensure = new EnsureJetCalibNodes("EnsureJetCalibNodes_forJES");
        // keep this modest; set RJ_JETCALIB_NODE_VERBOSE=1 for prints
        int nodeV = 0;
        if (const char* env = std::getenv("RJ_JETCALIB_NODE_VERBOSE")) nodeV = std::atoi(env);
        ensure->Verbosity(nodeV);
        se->registerSubsystem(ensure);

        if (vlevel > 0)
          std::cout << "[INFO] JES: enabling JetCalib (pp-style: pp data + isSim) => ensuring DST/TOWER exists\n";
    }

  // Optional: control JetCalib verbosity independently
  int jetcalV = 1; // default: show InitRun/process_event messages
  if (const char* env = std::getenv("RJ_JETCALIB_VERBOSITY"))
  {
      jetcalV = std::atoi(env);
  }

  for (const auto& jnm : RecoilJets::kJetRadii)
  {
      const std::string radKey = jnm.key;            // "r02" or "r04"
      const int   D = std::stoi(radKey.substr(1));   // 2 or 4
      const float R = 0.1f * D;                      // 0.2 or 0.4

      // Canonical node name that RecoilJets reads (from RecoilJets::kJetRadii)
      const std::string calibNode = jnm.pp_node;

      // Apply JES calibration for pp-like chains (pp data + pp-style SIM)
      // Run pp JES calibration in BOTH pp data and isSim (pp-style chains)
      const bool doJetCalib = (!isAuAuData);

      // If calibrating: build RAW jets to a separate node to avoid name collision
      const std::string rawNode = doJetCalib ? (calibNode + "_RAW") : calibNode;

      // ---------------------- JetReco (build jets) ----------------------
      const std::string recoName = std::string("JetsReco_") + (isSim ? "Sim_" : "Data_") + radKey;

      auto* jreco = new JetReco(recoName);
      jreco->add_input(new TowerJetInput(Jet::CEMC_TOWERINFO,    "TOWERINFO_CALIB"));
      jreco->add_input(new TowerJetInput(Jet::HCALIN_TOWERINFO,  "TOWERINFO_CALIB"));
      jreco->add_input(new TowerJetInput(Jet::HCALOUT_TOWERINFO, "TOWERINFO_CALIB"));

      jreco->add_algo(detail::fjAlgo(R), rawNode);
      jreco->set_algo_node("ANTIKT");
      jreco->set_input_node("TOWERINFO_CALIB");

      // If you want JetReco debug: export RJ_JETRECO_VERBOSITY=1
      int jetrecoV = 0;
      if (const char* env = std::getenv("RJ_JETRECO_VERBOSITY")) jetrecoV = std::atoi(env);
      jreco->Verbosity(jetrecoV);

      se->registerSubsystem(jreco);

      if (vlevel > 0)
      {
        std::cout << "[INFO] (" << (isSim ? "isSim" : "data") << ") reco jets: built "
                  << rawNode << " (R=" << R << ") from TOWERINFO_CALIB\n";
      }

      // ---------------------- JetCalib (apply JES) ----------------------
      if (doJetCalib)
      {
        auto* jcal = new JetCalib(std::string("JetCalib_") + radKey);

        // JetCalib reads RAW jets and writes CALIB jets into the canonical node
        jcal->set_InputNode(rawNode);
        jcal->set_OutputNode(calibNode);

        jcal->set_JetRadius(R);
        jcal->set_ZvrtxNode("GlobalVertexMap");  // what JetCalib expects

        // Full pp JES: Zvrtx + eta dependent
        jcal->set_ApplyZvrtxDependentCalib(true);
        jcal->set_ApplyEtaDependentCalib(true);

        jcal->Verbosity(jetcalV);
        se->registerSubsystem(jcal);

        // ------------------------------------------------------------
        // Targeted JES probe: prints ONCE, ONLY if Verbosity() >= 20
        // ------------------------------------------------------------
        auto* probe = new JetCalibOneEventProbe(std::string("JetCalibOneEventProbe_") + radKey,
                                               rawNode,
                                               calibNode,
                                               /*maxJetsToPrint=*/12);
        probe->Verbosity(vlevel);  // uses your RJ_VERBOSITY; probe prints only if >=20
        se->registerSubsystem(probe);

        if (vlevel > 0)
        {
          std::cout << "[INFO] JES calib enabled: " << rawNode << " -> " << calibNode
                    << " (R=" << R << ", Zvrtx+eta dependent, JetCalibVerbosity=" << jetcalV << ")\n";
        }
     }
  }
    
  // ---------------------- Truth jets -----------------------------------------
  // If useDSTTruthJets==true: they already exist from DST_JETS_IN.
  // If buildTruthJetsFromParticles==true: build them from TRUTH particles here.
  if (isSim && buildTruthJetsFromParticles)
    {
      if (vlevel > 0)
      {
        if (useDSTTruthJets)
          std::cout << "[INFO] (isSim) truth jets: ALSO building from TRUTH particles (QA 'both' mode)\n";
        else
          std::cout << "[INFO] (isSim) truth jets: building from TRUTH particles (no DST_JETS)\n";
      }

      for (const auto& jnm : RecoilJets::kJetRadii)
      {
        const std::string radKey = jnm.key;          // "r02" or "r04"
        const int   D = std::stoi(radKey.substr(1));
        const float R = 0.1f * D;

        auto* truthReco = new JetReco(std::string("TruthJetReco_FromParticles_") + radKey);
        auto* tji = new TruthJetInput(Jet::PARTICLE);
        tji->add_embedding_flag(1);  // pythia/herwig particles only
        truthReco->add_input(tji);

        // Node naming:
        //  - If we're NOT reading DST jets, keep canonical name "AntiKt_Truth_<rKey>"
        //  - If we ARE reading DST jets too (BOTH mode), avoid collisions
        const std::string truthNode = (useDSTTruthJets && buildTruthJetsAsAltNode)
                                    ? (std::string("AntiKt_TruthFromParticles_") + radKey)
                                    : (std::string("AntiKt_Truth_") + radKey);

          truthReco->add_algo(detail::fjAlgo(R), truthNode);
          truthReco->set_algo_node("ANTIKT");
          truthReco->set_input_node("TRUTH");

          // Same deal for truth-jet reco: keep JetReco quiet.
          truthReco->Verbosity(0);

          se->registerSubsystem(truthReco);

        if (vlevel > 0)
          std::cout << "[INFO] (isSim) truth jets: produced node " << truthNode
                    << " (R=" << R << ")\n";
      }
    }
    else if (isSim && useDSTTruthJets && vlevel > 0)
    {
      std::cout << "[INFO] (isSim) truth jets: using nodes from DST_JETS (no TruthJetInput reco)\n";
  }



  // --------------------------------------------------------------------
  // 5.  Run-information helper (optional but handy)
  // --------------------------------------------------------------------
  if (!isSim)
  {
      auto* trigInfo = new TriggerRunInfoReco();
      trigInfo->Verbosity(vlevel);
      se->registerSubsystem(trigInfo);
  }
  else
  {
      if (vlevel > 0) std::cout << "[isSim] skipping TriggerRunInfoReco" << std::endl;
  }
    
  // Build photon clusters
  auto* photonBuilder = new PhotonClusterBuilder("PhotonClusterBuilder");
  photonBuilder->set_input_cluster_node("CLUSTERINFO_CEMC");
  photonBuilder->set_output_photon_node("PHOTONCLUSTER_CEMC");
  photonBuilder->Verbosity(2);
  se->registerSubsystem(photonBuilder);

  auto* recoilJets = new RecoilJets(outRoot);

  recoilJets->setUseVzCut(true, 30.0);
  recoilJets->setIsolationWP(1.08128, 0.0299107, 1.0, 0.30, 0.0);

    
  // RecoilJets inherits SubsysReco::Verbosity(int)
  recoilJets->Verbosity(vlevel);
  if (verbose) std::cout << "[INFO] RJ_VERBOSITY → " << vlevel << '\n';
  // Pick analysis type for the module (isPP / isAuAu / isSim), case-insensitive.
  // Fallback: ≤ 53864 → isPP, > 53864 → isAuAu.
  std::string dtype = "isAuAu";
  if (const char* env = std::getenv("RJ_DATASET"))
  {
            std::string s = detail::trim(std::string(env));
            std::string sLower = s;
            std::transform(sLower.begin(), sLower.end(), sLower.begin(),
                           [](unsigned char c){ return std::tolower(c); });

            if (sLower == "issim" || sLower == "sim")      dtype = "isSim";
            else if (sLower == "ispp" || sLower == "pp")   dtype = "isPP";
            else                                           dtype = "isAuAu";
    }
    else
    {
            dtype = (run <= 53864) ? "isPP" : "isAuAu";
  }

  if (verbose) std::cout << "[INFO] RJ_DATASET → " << dtype << '\n';
  recoilJets->setDataType(dtype);

  se->registerSubsystem(recoilJets);

  //--------------------------------------------------------------------
  // 6.  Run
  //--------------------------------------------------------------------
  try
  {
    if (verbose) std::cout << "[INFO] Starting event loop …\n";
    se->run(nEvents);
    se->End();
    if (verbose) std::cout << "[INFO] Finished successfully.\n";
  }
  catch (const std::exception& e)
  {
    detail::bail(std::string("exception in Fun4All: ") + e.what());
  }

//  //--------------------------------------------------------------------
//  // 7.  Clean exit
//  //--------------------------------------------------------------------
//  gSystem->Exit(0);
}

#endif   // ROOT_VERSION guard
