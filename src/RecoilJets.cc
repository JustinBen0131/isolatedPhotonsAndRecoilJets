//==========================================================================
//  sPHENIX EMCal × sEPD × MBD correlator
//  Implementation file  – no duplicated definitions
//==========================================================================
#include "RecoilJets.h"
//––– Fun4All / PHOOL -------------------------------------------------------
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/getClass.h>
#include <phool/recoConsts.h>
#include <jetbase/JetContainer.h>
#include <array>
//––– ROOT & CLHEP ----------------------------------------------------------
#include <TProfile.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <TMath.h>
#include <TH2Poly.h>
#include <TKey.h>
#include <CLHEP/Vector/ThreeVector.h>
#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>
//––– sPHENIX objects -------------------------------------------------------
#include <globalvertex/GlobalVertex.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoDefs.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/RawCluster.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/RawClusterUtility.h>
#include <mbd/MbdPmtHit.h>
#include <mbd/MbdGeom.h>
#include <mbd/MbdOut.h>
#include <ffarawobjects/Gl1Packet.h>
#include <mbd/MbdPmtContainer.h>
#include <epd/EpdReco.h>
#include <epd/EpdGeom.h>
#include <centrality/CentralityInfo.h>
#include <calotrigger/MinimumBiasInfo.h>
#include <calotrigger/MinimumBiasClassifier.h>   // optional but handy

#include <eventplaneinfo/Eventplaneinfo.h>
#include <eventplaneinfo/Eventplaneinfov1.h>
#include <eventplaneinfo/EventplaneinfoMap.h>

// Standard C++ -------------------------------------------------------------
#include <atomic>
#include <algorithm>   // std::clamp
#include <cmath>       // std::cosh, std::hypot, std::fmod
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <regex>
#include <tuple>

#ifdef _OPENMP
  #include <omp.h>
#endif

//–––––––– helpers ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
#define CLR_BLUE   "\033[1;34m"
#define CLR_CYAN   "\033[1;36m"
#define CLR_GREEN  "\033[1;32m"
#define CLR_YELLOW "\033[1;33m"
#define CLR_MAGENTA "\033[1;35m"
#define CLR_RESET  "\033[0m"

#undef  LOG
#define LOG(lvl, colour, msg)                                           \
  do {                                                                  \
    if (static_cast<int>(Verbosity()) >= static_cast<int>(lvl))         \
      std::cout << colour << msg << CLR_RESET << std::endl;             \
  } while (false)

/** Always print, independent of Verbosity() */
#define PROGRESS(MSG)                                                     \
    do {                                                                  \
        if (static_cast<int>(Verbosity()) >= 1)                           \
            std::cout << CLR_CYAN << MSG << CLR_RESET << std::endl;       \
    } while (false)
//==========================================================================
//  ctor
//==========================================================================
RecoilJets::RecoilJets(const std::string& outFile)
: SubsysReco("RecoilJets"),
  Outfile(outFile)
{
  if (Outfile.empty())
  {
    std::cerr << "[FATAL] output filename is empty.\n";
    std::exit(EXIT_FAILURE);
  }
}

RecoilJets::~RecoilJets()
{
  
}

// ----------------------------------------------------------------------
//  getCentralitySlice  – returns true if the current event falls into
//                        one of the user–defined centrality intervals
//                        (lo ≤ cent < hi).  On success it also builds
//                        the “_lo_hi” tag that all QA routines use
//                        when addressing centrality‑specific histograms.
//
//  • lo,hi   – returned lower / upper edge (undefined if the function
//              returns false)
//  • tag     – empty if the function returns false
// ----------------------------------------------------------------------
bool RecoilJets::getCentralitySlice(int& lo,
                                              int& hi,
                                              std::string& tag) const
{
  lo = hi = -1;
  for (std::size_t i = 0; i + 1 < m_centEdges.size(); ++i)
    if (m_centBin >= m_centEdges[i] && m_centBin < m_centEdges[i + 1])
    { lo = m_centEdges[i]; hi = m_centEdges[i + 1]; break; }

  const bool hasSlice = (lo >= 0);
  tag.clear();
  if (hasSlice)
  {
    std::ostringstream ss;
    ss << '_' << lo << '_' << hi;
    tag = ss.str();
  }
  return hasSlice;
}


int RecoilJets::Init(PHCompositeNode* topNode)
{
  LOG(1, CLR_BLUE, "[Init] RecoilJets – starting");

  /* 0.  book-keeping & QA histograms --------------------------------- */
  out = new TFile(Outfile.c_str(), "RECREATE");
  LOG(1, CLR_GREEN, "[Init] opened output file: " << Outfile);

  trigAna = new TriggerAnalyzer();
  LOG(1, CLR_GREEN, "[Init] booking scalar QA histograms …");
  createHistos_Data();
    
  /* 1.  optional DST node-tree dump ---------------------------------- */
  if (Verbosity() >= 2)           // ← adjust threshold as desired
  {
    std::cout << CLR_CYAN
              << "\n[Init] ── DST node-tree dump ────────────────────────────"
              << CLR_RESET << std::endl;

    /* depth-first walk implemented with a std::function so that the
       lambda can recurse without shadowing problems                       */
    std::function<void(PHCompositeNode*, int)> dumpTree =
      [&](PHCompositeNode* node, int depth)
    {
      if (!node) return;

      /* print this node ------------------------------------------------ */
      const std::string indent(depth * 3, ' ');
      std::cout << indent << node->getName()
                << " (" << node->getType() << ")";

      /* for IO-data nodes also print the contained class name ---------- */
      if (node->getType() == "PHIODataNode")
        std::cout << " <" << node->getClass() << '>';

      std::cout << '\n';

      /* iterate over children with the *real* PHOOL API ---------------- */
      PHNodeIterator it(node);
      auto& kids = it.ls();                       // PHPointerList<PHNode>
      for (size_t i = 0; i < kids.length(); ++i)
      {
        PHNode* child = kids[i];
        if (!child) continue;

        if (auto* comp = dynamic_cast<PHCompositeNode*>(child))
        {
          dumpTree(comp, depth + 1);              // recurse
        }
        else
        {
          const std::string ind2((depth + 1) * 3, ' ');
          std::cout << ind2 << child->getName()
                    << " (" << child->getType() << ")";
          if (child->getType() == "PHIODataNode")
            std::cout << " <" << child->getClass() << '>';
          std::cout << '\n';
        }
      }
    };

    /* pick the correct root: use argument if non-null, else global ---- */
    PHCompositeNode* root = topNode
                              ? topNode
                              : Fun4AllServer::instance()->topNode();

    dumpTree(root, 0);

    std::cout << CLR_CYAN
              << "[Init] ──────────────────────────────────────────────────\n"
              << CLR_RESET << std::endl;
  }

  LOG(1, CLR_BLUE, "[Init] RecoilJets – done");
  return Fun4AllReturnCodes::EVENT_OK;
}


int RecoilJets::InitRun(PHCompositeNode* topNode)
{
  /* 0. banner -------------------------------------------------------- */
  const uint64_t run = recoConsts::instance()->get_uint64Flag("TIMESTAMP", 0);
  LOG(1, CLR_BLUE, "[InitRun] ------------------------------------------------------------");
  LOG(1, CLR_BLUE, "[InitRun] Starting InitRun  –  TIMESTAMP = " << run);

  /* 3. sanity‑check user centrality edges --------------------------- */
  if (m_centEdges.empty())
    LOG(0, CLR_YELLOW, "[InitRun] WARNING: centrality edges vector is EMPTY");
  else
  {
    bool mono = std::is_sorted(m_centEdges.begin(), m_centEdges.end());
    if (!mono)
      LOG(0, CLR_YELLOW, "[InitRun] WARNING: centrality edges not monotonic");
  }

  LOG(1, CLR_BLUE, "[InitRun] InitRun completed successfully");
  return Fun4AllReturnCodes::EVENT_OK;
}



void RecoilJets::createHistos_Data()
{
  for (const auto& kv : triggerNameMap)
  {
    const std::string trig = kv.second;
    if (Verbosity() > 1)
      std::cout << CLR_BLUE << "  ├─ trigger \"" << trig
                << "\" – booking scalar QA" << CLR_RESET << std::endl;

    out->mkdir(trig.c_str())->cd();
    HistMap& H = qaHistogramsByTrigger[trig];

    bookEventPlaneCentralityQA(trig, H);
    bookPhotonJetQA(trig, H);
      
    H["h_vertexZ"] = new TH1F(("h_vertexZ_" + trig).c_str(),
                                "Primary vertex z;z_{vtx} [cm]",
                                240, -60., 60.);
      
    H["h_centrality"] =
          new TH1F(("h_centrality_" + trig).c_str(),
                   "Centrality percentile (MBD);centrality [%];Events",
                   100, 0., 100.);
      
    H["cnt_"+trig+"_raw"] =
          new TH1I(("cnt_"+trig+"_raw").c_str(),
                   (trig+" – raw bit fired;flag;Events").c_str(),
                   1, 0.5, 1.5);

    H["cnt_"+trig+"_live"] =
          new TH1I(("cnt_"+trig+"_live").c_str(),
                   (trig+" – live bit fired;flag;Events").c_str(),
                   1, 0.5, 1.5);

    H["cnt_"+trig+"_scaled"] =
          new TH1I(("cnt_"+trig+"_scaled").c_str(),
                   (trig+" – scaled bit fired;flag;Events").c_str(),
                   1, 0.5, 1.5);
      
      /* ── vertex‑cut compliance histogram ─────────────────────────────── */
    if (int vCut = extractVtxCut(trig); vCut > 0)
    {
        auto* h = new TH1I(Form("h_vtxRelToCut_%s", trig.c_str()),
                           Form("|z_{vtx}| vs trigger‑cut = %d cm;relation;Events", vCut),
                           3, 0.5, 3.5);
        h->GetXaxis()->SetBinLabel(1, "< cut");
        h->GetXaxis()->SetBinLabel(2, "= cut");
        h->GetXaxis()->SetBinLabel(3, "> cut");
        H[h->GetName()] = h;
    }
      
    out->cd();
  }
}



// ======================================================================
//  process_event – one-event driver, pp/Au+Au aware trigger gating
// ======================================================================
int RecoilJets::process_event(PHCompositeNode* topNode)
{
  /* ------------------------------------------------------------------ */
  /* 0) Banner & counter                                                */
  /* ------------------------------------------------------------------ */
  ++event_count;
  PROGRESS("=================================   event "
           << std::setw(4) << event_count
           << "   =====================================  (Verb="
           << Verbosity() << ')');

  /* ------------------------------------------------------------------ */
  /* 1) Mandatory nodes                                                 */
  /* ------------------------------------------------------------------ */
  LOG(4, CLR_BLUE, "  [process_event] – node sanity");
  if (!fetchNodes(topNode))
  {
    LOG(4, CLR_YELLOW,
        "    mandatory node(s) missing → ABORTEVENT");
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  /* ------------------------------------------------------------------ */
  /* 2) Trigger gating (pp vs Au+Au)                                    */
  /*     - pp:  TriggerAnalyzer, require "MBD N&S >= 1"                 */
  /*     - Au+Au: use existing firstEventCuts() (MB+trigger gate)       */
  /* ------------------------------------------------------------------ */
  std::vector<std::string> activeTrig;  // short keys used for booking/filling

  if (!m_isAuAu)
  {
    // -------- p+p path -------------------------------------------------
    if (!trigAna)
    {
      std::cerr << "[ERROR] No TriggerAnalyzer pointer in pp mode!\n";
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    if (Verbosity() > 0)
      std::cout << "[DEBUG] pp mode: decoding GL1 triggers …" << std::endl;

    trigAna->decodeTriggers(topNode);

    // Build "short-name" list from DB-name map
    activeTrig.clear();
    activeTrig.reserve(triggerNameMap_pp.size());
    for (const auto& kv : triggerNameMap_pp)
    {
      const std::string& dbName   = kv.first;   // GL1/DB name
      const std::string& shortKey = kv.second;  // histogram-friendly short name
      if (trigAna->didTriggerFire(dbName))
      {
        activeTrig.push_back(shortKey);
        if (Verbosity() > 0)
          std::cout << "Trigger fired: \"" << dbName
                    << "\" → short \"" << shortKey << "\"\n";
      }
    }

    // Require "MBD N&S >= 1" in pp
    const char* requiredShort = "MBD_NandS_geq_1";
    bool haveRequired = false;
    for (const auto& s : activeTrig) { if (s == requiredShort) { haveRequired = true; break; } }

    if (!haveRequired)
    {
      if (Verbosity() > 0)
        std::cout << "[trigger-gate] pp event skipped: required \"MBD N&S >= 1\" not active.\n";
      // Soft-skip: do not fill anything for this event
      return Fun4AllReturnCodes::EVENT_OK;
    }

    // Ensure we always have at least one tag
    if (activeTrig.empty()) activeTrig.emplace_back("ALL");
  }
  else
  {
    // -------- Au+Au path ----------------------------------------------
    if (!firstEventCuts(topNode, activeTrig))
    {
      ++m_evtNoTrig;
      LOG(4, CLR_YELLOW, "    event rejected by MB/Trigger gate – skip");
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    // Some configurations of firstEventCuts may not supply names; keep consistent
    if (activeTrig.empty()) activeTrig.emplace_back("ALL");
  }

  /* ------------------------------------------------------------------ */
  /* 3) Vertex-z QA                                                     */
  /* ------------------------------------------------------------------ */
  for (const auto& t : activeTrig)
  {
    auto itTrig = qaHistogramsByTrigger.find(t);
    if (itTrig != qaHistogramsByTrigger.end())
    {
      auto& H = itTrig->second;
      if (auto hit = H.find("h_vertexZ"); hit != H.end())
        static_cast<TH1F*>(hit->second)->Fill(m_vz);
    }
  }

  /* ------------------------------------------------------------------ */
  /* 4) Centrality lookup (Au+Au only)                                  */
  /*     Keep m_centBin = -1 in pp (acts as minimum-bias)               */
  /* ------------------------------------------------------------------ */
  if (m_isAuAu)
  {
    CentralityInfo* central =
      findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");

    if (!central)
    {
      LOG(4, CLR_YELLOW,
          "    CentralityInfo node missing (Au+Au) – ABORTEVENT");
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    const float centile =
      central->get_centrality_bin(CentralityInfo::PROP::mbd_NS);

    if (!std::isfinite(centile) || centile < 0.f)
    {
      LOG(4, CLR_YELLOW,
          "    invalid mbd_NS centile – treating as minimum-bias (0–100%)");
      m_centBin = -1;
    }
    else
    {
      m_centBin = static_cast<int>(centile);
      LOG(5, CLR_GREEN, "    centrality bin = " << m_centBin << '%');
    }

    // Fill centrality histogram if booked under each active trigger
    if (centile >= 0.f && centile <= 100.f)
    {
      for (const auto& t : activeTrig)
      {
        auto itTrig = qaHistogramsByTrigger.find(t);
        if (itTrig != qaHistogramsByTrigger.end())
        {
          auto& H = itTrig->second;
          if (auto hc = H.find("h_centrality"); hc != H.end())
            static_cast<TH1F*>(hc->second)->Fill(centile);
        }
      }
    }
  }
  else
  {
    // pp: no centrality
    m_centBin = -1;
  }

  /* ------------------------------------------------------------------ */
  /* 5) Vertex-z guard (applies to both, primarily relevant for Au+Au)  */
  /* ------------------------------------------------------------------ */
  if (!std::isfinite(m_vz) || std::abs(m_vz) > 60.0)
  {
    LOG(4, CLR_YELLOW,
        "    Vertex-z (" << m_vz << " cm) outside bounds – skip event");
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  /* ------------------------------------------------------------------ */
  /* 6) Done – downstream (QA/correlations/etc.) will use activeTrig    */
  /* ------------------------------------------------------------------ */
  LOG(4, CLR_GREEN, "  [process_event] – completed OK");
  return Fun4AllReturnCodes::EVENT_OK;
}


//==========================================================================
//  fetchNodes – check presence of all required nodes & cache pointers
//==========================================================================
bool RecoilJets::fetchNodes(PHCompositeNode* top)
{
    
  /* ------------------------------------------------------------------ */
  /* 0.  Minimum‑bias information         */
  /* ------------------------------------------------------------------ */
  m_isMinBias = false;                                // reset per event

  if (auto* mbInfo = findNode::getClass<MinimumBiasInfo>(top,
                                                           "MinimumBiasInfo"))
        m_isMinBias = mbInfo->isAuAuMinimumBias();
  else
        LOG(1, CLR_YELLOW,
            "  – MinimumBiasInfo node missing (treating as !MB)");

  /* NOTE: Actual rejection is done in firstEventCuts(), so we only
   *       cache m_isMinBias here and keep processing. */

    
  /* ––– primary vertex –––––––––––––––––––––––––––––––––––––––––––––––– */
  GlobalVertexMap* vmap = findNode::getClass<GlobalVertexMap>(top,"GlobalVertexMap");
  m_vtx = nullptr; m_vx = m_vy = m_vz = 0.;

  if (!vmap)          { LOG(1, CLR_YELLOW, "  – GlobalVertexMap node **missing** → skip event"); return false; }
  if (vmap->empty())  { LOG(2, CLR_YELLOW, "  – GlobalVertexMap is **empty** → skip event");     return false; }

  m_vtx = vmap->begin()->second;
  if (!m_vtx)         { LOG(1, CLR_YELLOW, "  – vertex pointer null → skip event");              return false; }

  m_vx = m_vtx->get_x();
  m_vy = m_vtx->get_y();
  m_vz = m_vtx->get_z();

  /* ––– calorimeter towers & geometry –––––––––––––––––––––––––––––––– */
  m_calo.clear();
  for (const auto& ci : m_caloInfo)
  {
    const std::string node = std::get<0>(ci),
                      geo  = std::get<1>(ci),
                      lbl  = std::get<2>(ci);

    auto* tw = findNode::getClass<TowerInfoContainer>(top, node);
    auto* ge = findNode::getClass<RawTowerGeomContainer>(top, geo);
    if (!tw || !ge)
    { LOG(2, CLR_YELLOW, "  – missing " << lbl << " nodes → skip"); return false; }

    m_calo[lbl] = { tw, ge, 0. };
  }

  /* ––– remaining detectors –––––––––––––––––––––––––––––––––––––––––– */
  m_sepd     = findNode::getClass<TowerInfoContainer>(top, "TOWERINFO_CALIB_SEPD");
  m_mbdpmts  = findNode::getClass<MbdPmtContainer  >(top, "MbdPmtContainer");
  m_mbdgeom  = findNode::getClass<MbdGeom         >(top, "MbdGeom");
  m_epdgeom  = findNode::getClass<EpdGeom         >(top, "TOWERGEOM_EPD");
  m_epmap    = findNode::getClass<EventplaneinfoMap>(top, "EventplaneinfoMap");
  m_clus     = findNode::getClass<RawClusterContainer>(top, "CLUSTERINFO_CEMC");

  const bool ok_sepd = (m_sepd && m_epdgeom);
  const bool ok_mbd  = (m_mbdpmts && m_mbdgeom);

  if (!ok_sepd || !ok_mbd)
  {
    LOG(2, CLR_YELLOW, "  – missing mandatory SEPD and/or MBD nodes → skip event");
    return false;
  }
  return true;
}



int RecoilJets::ResetEvent(PHCompositeNode*)
{
  
  m_centBin = -1;

  return Fun4AllReturnCodes::EVENT_OK;
}

int RecoilJets::Reset(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

// --------------------------------------------------------------------------
//  End – enhanced diagnostics, robust against dangling pointers
// --------------------------------------------------------------------------
int RecoilJets::End(PHCompositeNode*)
{
  auto warn = [&](const std::string& m)
  { std::cerr << CLR_YELLOW << "[End] " << m << CLR_RESET << '\n'; };

  auto info = [&](int lvl, const std::string& m)
  { LOG(lvl, CLR_GREEN, "[End] " << m); };

  //--------------------------------------------------------------------
  // 1. Basic checks on the output file pointer
  //--------------------------------------------------------------------
  if (!out)           { warn("TFile* 'out' is nullptr – nothing to write");  return Fun4AllReturnCodes::ABORTEVENT; }
  if (!out->IsOpen()) { warn("Output file is *not* open ("+std::string(out->GetName())+")"); return Fun4AllReturnCodes::ABORTEVENT; }

  //--------------------------------------------------------------------
  // 2. Write histograms trigger‑by‑trigger
  //--------------------------------------------------------------------
  std::size_t nHistExpected = 0, nHistWritten = 0;

  for (auto& [trig, hMap] : qaHistogramsByTrigger)
  {
    TDirectory* dir = out->GetDirectory(trig.c_str());
    if (!dir) dir = out->mkdir(trig.c_str());
    dir->cd();

    std::size_t nExpThis = 0, nWrtThis = 0;

    for (auto& [key, obj] : hMap)
    {
      ++nHistExpected; ++nExpThis;

      TH1* h = dynamic_cast<TH1*>(obj);
      if (!h)              { warn("Object '"+key+"' (trigger "+trig+") is not TH1 – skipped"); continue; }
      if (h->GetEntries()==0)
      { if (Verbosity()>1) warn("Histogram '"+key+"' (trigger "+trig+") has 0 entries – skipped"); continue; }

      try
      {
        if (h->Write("", TObject::kOverwrite) > 0) { ++nHistWritten; ++nWrtThis; }
        else warn("Write() returned 0 for '"+key+"' (trigger "+trig+")");
      }
      catch (const std::exception& e)
      { warn("Exception while writing '"+key+"' – "+std::string(e.what())); }
    }
    info(1, "trigger '"+trig+"': "+std::to_string(nWrtThis)+" / "
               +std::to_string(nExpThis)+" histograms written");
    out->cd();
  }

  //--------------------------------------------------------------------
  // 3.  Human‑readable summary  (must run *before* the file is deleted)
  //--------------------------------------------------------------------
  if (Verbosity() > 0)
  {
    std::cout << "\n\033[1mHistogram summary\033[0m\n"
              << "\033[1mTrigger                        │ Histogram                           │  Entries\033[0m\n"
              << "-------------------------------------------------------------------------------\n";

    for (const auto& [trig, hMap] : qaHistogramsByTrigger)
      for (const auto& [key, obj]  : hMap)
        if (const TH1* h = dynamic_cast<const TH1*>(obj))
        {
            if (h->GetEntries() == 0) continue;          // <‑‑ only list filled ones
            std::cout << std::left << std::setw(30) << trig << " │ "
            << std::setw(32) << key  << " │ "
            << std::right<< std::setw(10)
            << static_cast<Long64_t>(h->GetEntries()) << '\n';
        }
    if (!m_totStat.empty())
    {
      std::cout << "-------------------------------------------------------------------------------\n"
                << "\n\033[1mπ0‑QA cut summary (all events)\033[0m\n"
                << "\033[1mcut‑key                                    │ tested        passed    eff[%]\033[0m\n"
                << "--------------------------------------------------------------------------\n";
      for (const auto& [key, s] : m_totStat)
      {
        const double eff = s.tested ? 100.*s.passed / s.tested : 0.;
        std::cout << std::left  << std::setw(44) << key << " │ "
                  << std::right << std::setw(12) << s.tested
                  << std::setw(12) << s.passed
                  << std::setw(9)  << std::fixed << std::setprecision(2) << eff << '\n';
      }
      std::cout << "--------------------------------------------------------------------------\n";
    }
  }
  out->cd("triggerQA");
  if (h_MBTrigCorr && h_MBTrigCorr->GetEntries() > 0) h_MBTrigCorr->Write();

  //--------------------------------------------------------------------
  // 4.  Write footer & close the file
  //--------------------------------------------------------------------
  if (Verbosity() >= 1)
      std::cout << "\nOutput ROOT file →  " << out->GetName() << "\n\n";
    
  info(1, "writing TFile footer and closing ("+std::to_string(nHistWritten)
           +" / "+std::to_string(nHistExpected)+" objects written)");

  try      { out->Close(); }
  catch (const std::exception& e)
  { warn("Exception during TFile::Write/Close – "+std::string(e.what())); }

  delete out; out = nullptr;          // safe: we no longer dereference histos
  info(0, "Done.");
  return Fun4AllReturnCodes::EVENT_OK;
}


void RecoilJets::Print(const std::string&) const
{ /* nothing to print beyond ROOT histograms */ }
