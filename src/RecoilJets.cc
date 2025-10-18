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
#include <jetbase/Jet.h>
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
#include <calobase/PhotonClusterContainer.h>
#include <calobase/PhotonClusterv1.h>
#include <clusteriso/ClusterIso.h>

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


using namespace PhoIDCuts;
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

  /* NOTE: Actual rejection is done in firstEventCuts(), so only
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
  m_photons  = findNode::getClass<PhotonClusterContainer>(top, "PHOTONCLUSTER_CEMC"); // from PhotonClusterBuilder

  // --- UE-subtracted jet containers (SUB1) --------------------------------
  m_jets.clear();
  for (const auto& rk : kJetRadii)
  {
      const std::string rKey = rk.first;   // "r02"
      const std::string node = rk.second;  // e.g. "AntiKt_TowerInfo_HIRecoSeedsSub_r02"
      auto* jc = findNode::getClass<JetContainer>(top, node);
      if (!jc)
        LOG(2, CLR_YELLOW, "  – jet node missing: " << node << " (xJ disabled for " << rKey << ")");
      m_jets[rKey] = jc; // may be nullptr; we check later
  }

  const bool ok_sepd = (m_sepd && m_epdgeom);
  const bool ok_mbd  = (m_mbdpmts && m_mbdgeom);

  if (!ok_sepd || !ok_mbd)
  {
        LOG(2, CLR_YELLOW, "  – missing mandatory SEPD and/or MBD nodes → skip event");
        return false;
  }
  return true;
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


int RecoilJets::InitRun(PHCompositeNode* /*topNode*/)
{
  const uint64_t run = recoConsts::instance()->get_uint64Flag("TIMESTAMP", 0);
  LOG(1, CLR_BLUE, "[InitRun] ------------------------------------------------------------");
  LOG(1, CLR_BLUE, "[InitRun] Starting InitRun  –  TIMESTAMP = " << run);

  // Dataset flag + configured binning
  LOG(1, CLR_GREEN, "[InitRun] Dataset      : " << (m_isAuAu ? "Au+Au" : "p+p"));
  {
    std::ostringstream os;
    os << "[InitRun] gamma-ET bins: {";
    for (std::size_t i = 0; i+1 < m_gammaEtBins.size(); ++i)
      os << (i? ", ":" ") << m_gammaEtBins[i] << "–" << m_gammaEtBins[i+1];
    os << " }";
    LOG(1, CLR_GREEN, os.str());
  }
  if (m_isAuAu)
  {
    if (m_centEdges.empty())
      LOG(0, CLR_YELLOW, "[InitRun] WARNING: centrality edges vector is EMPTY");
    else
    {
      bool mono = std::is_sorted(m_centEdges.begin(), m_centEdges.end());
      if (!mono)
        LOG(0, CLR_YELLOW, "[InitRun] WARNING: centrality edges not monotonic");
      std::ostringstream oc;
      oc << "[InitRun] centrality bins: {";
      for (std::size_t i = 0; i+1 < m_centEdges.size(); ++i)
        oc << (i? ", ":" ") << m_centEdges[i] << "–" << m_centEdges[i+1];
      oc << " }";
      LOG(1, CLR_GREEN, oc.str());
    }
  }

  LOG(1, CLR_BLUE, "[InitRun] InitRun completed successfully");
  return Fun4AllReturnCodes::EVENT_OK;
}


void RecoilJets::createHistos_Data()
{
  // Book per-trigger counters for whichever data type we’re running
  if (m_isAuAu)
  {
    for (const auto& kv : triggerNameMapAuAu)    // kv: std::pair<int,std::string>
    {
      const std::string trig = kv.second;
      if (Verbosity() > 1)
        std::cout << CLR_BLUE << "  ├─ trigger \"" << trig
                  << "\" – booking trigger counter" << CLR_RESET << std::endl;

      TDirectory* dir = out->GetDirectory(trig.c_str());
      if (!dir) dir = out->mkdir(trig.c_str());
      dir->cd();

      HistMap& H = qaHistogramsByTrigger[trig];
      const std::string hname = "cnt_" + trig;
      if (H.find(hname) == H.end())
      {
        auto* h = new TH1I(hname.c_str(), (hname+";count;entries").c_str(), 1, 0.5, 1.5);
        h->GetXaxis()->SetBinLabel(1, "count");
        H[hname] = h;
      }
      out->cd();
    }
  }
  else
  {
    for (const auto& kv : triggerNameMap_pp)     // kv: std::pair<std::string,std::string>
    {
      const std::string trig = kv.second;
      if (Verbosity() > 1)
        std::cout << CLR_BLUE << "  ├─ trigger \"" << trig
                  << "\" – booking trigger counter" << CLR_RESET << std::endl;

      TDirectory* dir = out->GetDirectory(trig.c_str());
      if (!dir) dir = out->mkdir(trig.c_str());
      dir->cd();

      HistMap& H = qaHistogramsByTrigger[trig];
      const std::string hname = "cnt_" + trig;
      if (H.find(hname) == H.end())
      {
        auto* h = new TH1I(hname.c_str(), (hname+";count;entries").c_str(), 1, 0.5, 1.5);
        h->GetXaxis()->SetBinLabel(1, "count");
        H[hname] = h;
      }
      out->cd();
    }
  }
}



bool RecoilJets::firstEventCuts(PHCompositeNode* topNode,
                                std::vector<std::string>& activeTrig)
{
  activeTrig.clear();

  if (!m_isAuAu)
  {
    // ---------- p+p path: require "MBD N&S >= 1" via TriggerAnalyzer ----------
    if (!trigAna) return false;

    trigAna->decodeTriggers(topNode);

    for (const auto& kv : triggerNameMap_pp)
    {
      const std::string& dbName   = kv.first;   // GL1/DB name
      const std::string& shortKey = kv.second;  // histogram-friendly short key
      if (trigAna->didTriggerFire(dbName))
        activeTrig.push_back(shortKey);
    }

    // Strict pp gate: must have the MB short key present
    const char* requiredShort = "MBD_NandS_geq_1";
    const bool haveRequired =
      std::find(activeTrig.begin(), activeTrig.end(), requiredShort) != activeTrig.end();
    if (!haveRequired) return false;
  }
  else
  {
    // ---------- Au+Au path: MB + any scaled bit from the map ----------
    if (!m_isMinBias) return false;

    uint64_t wScaled = 0;
    if (auto* gl1 = findNode::getClass<Gl1Packet>(topNode, "GL1Packet"))
      wScaled = gl1->lValue(0, "ScaledVector");
    else if (auto* gl1b = findNode::getClass<Gl1Packet>(topNode, "14001"))
      wScaled = gl1b->lValue(0, "ScaledVector");

    if (wScaled)
    {
      const auto bits = extractTriggerBits(wScaled, static_cast<int>(event_count));
      for (const auto& [bit, key] : triggerNameMapAuAu)
        if (checkTriggerCondition(bits, bit))
          activeTrig.push_back(key);
    }

    if (activeTrig.empty()) return false;
  }

  // ---------- Global |vz| veto (applies to both datasets) ----------
  if (m_useVzCut && std::fabs(m_vz) >= m_vzCut) return false;

  // Uniform downstream behavior
  if (activeTrig.empty()) activeTrig.emplace_back("ALL");

  return true;
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
  std::cout << "==================== processing event "
              << std::setw(6) << event_count
              << " ====================" << std::endl;

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
  /* 2) Trigger gating (pp & Au+Au) — unified in firstEventCuts()       */
  /* ------------------------------------------------------------------ */
  std::vector<std::string> activeTrig;
  if (!firstEventCuts(topNode, activeTrig))
  {
      ++m_evtNoTrig;
      LOG(4, CLR_YELLOW, "    event rejected by MB/Trigger/Vz gate – skip");
      return Fun4AllReturnCodes::ABORTEVENT;
  }

  /* ------------------------------------------------------------------ */
  /* 3) Trigger counters (one per trigger) + Vertex-z QA                */
  /* ------------------------------------------------------------------ */

  // Bump the per-trigger counter once per accepted event
  for (const auto& t : activeTrig)
  {
      auto itTrig = qaHistogramsByTrigger.find(t);
      if (itTrig != qaHistogramsByTrigger.end())
      {
        auto& H = itTrig->second;
        const std::string hname = "cnt_" + t;
        if (auto hc = H.find(hname); hc != H.end())
          static_cast<TH1I*>(hc->second)->Fill(1);

        // Optional existing vertex QA (kept intact)
        if (auto hvz = H.find("h_vertexZ"); hvz != H.end())
          static_cast<TH1F*>(hvz->second)->Fill(m_vz);
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

  processCandidates(topNode, activeTrig);
    
  LOG(4, CLR_GREEN, "  [process_event] – completed OK");
  return Fun4AllReturnCodes::EVENT_OK;
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
  }

  // ------------------ Analysis summary (verbosity-controlled) ------------------
  if (Verbosity() >= 1)
  {
      std::cout << "\n\033[1mSelection summary (dataset: " << (m_isAuAu ? "Au+Au" : "p+p")
                << ", events=" << event_count << ")\033[0m\n";
      for (const auto& kvT : m_catByTrig)
      {
        const std::string& trig = kvT.first;
        std::cout << "\n\033[1mTrigger: " << trig << "\033[0m\n";
        std::cout << "slice                                 |   seen |  tight | nonTgt | isoPass | isoFail |   SBtot |  SBpass\n";
        std::cout << "--------------------------------------+--------+--------+--------+---------+---------+--------+--------\n";

        // Ordered print by slice label
        std::vector<std::string> keys;
        keys.reserve(kvT.second.size());
        for (auto& s : kvT.second) keys.push_back(s.first);
        std::sort(keys.begin(), keys.end());

        for (const auto& sfx : keys)
        {
          const CatStat& S = kvT.second.at(sfx);
          std::ostringstream lab; lab << (sfx.empty() ? "<all>" : sfx);
          std::cout << std::left  << std::setw(38) << lab.str() << " | "
                    << std::right << std::setw(6)  << S.seen     << " | "
                    << std::setw(6)  << S.tight    << " | "
                    << std::setw(6)  << S.nonTight << " | "
                    << std::setw(7)  << S.isoPass  << " | "
                    << std::setw(7)  << S.isoFail  << " | "
                    << std::setw(6)  << S.idSB_total << " | "
                    << std::setw(6)  << S.idSB_pass  << "\n";
        }
      }

      // Histogram fill counts (analysis set)
      if (!m_histFill.empty())
      {
        std::cout << "\n\033[1mHistogram fills (analysis set)\033[0m\n";
        std::cout << "trigger::histogram                               | fills\n";
        std::cout << "-----------------------------------------------+-------\n";
        // Stable, alphabetical
        std::vector<std::pair<std::string,std::size_t>> vv(m_histFill.begin(), m_histFill.end());
        std::sort(vv.begin(), vv.end(),
                  [](auto& a, auto& b){ return a.first < b.first; });
        for (const auto& p : vv)
        {
          std::cout << std::left  << std::setw(47) << p.first
                    << " | " << std::right << std::setw(5) << p.second << "\n";
        }
      }
  }
  // ---------------------------------------------------------------------------

  // 4.  Write footer & close the file
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


// Build SSVars from PhotonClusterv1 (names come from PhotonClusterBuilder)
RecoilJets::SSVars RecoilJets::makeSSFromPhoton(const PhotonClusterv1* pho, double et) const
{
  SSVars v{};
  const double weta_cogx = pho->get_shower_shape_parameter("weta_cogx");
  const double wphi_cogx = pho->get_shower_shape_parameter("wphi_cogx");
  const double et1       = pho->get_shower_shape_parameter("et1");

  const double e11 = pho->get_shower_shape_parameter("e11");
  const double e33 = pho->get_shower_shape_parameter("e33");
  const double e32 = pho->get_shower_shape_parameter("e32");
  const double e35 = pho->get_shower_shape_parameter("e35");

  v.weta_cogx     = weta_cogx;
  v.wphi_cogx     = wphi_cogx;
  v.et1           = et1;
  v.e11_over_e33  = (e33 > 0.0) ? (e11 / e33) : 0.0;
  v.e32_over_e35  = (e35 > 0.0) ? (e32 / e35) : 0.0;
  v.et_gamma      = et;
  return v;
}

// Centralized candidate processing for the event
void RecoilJets::processCandidates(PHCompositeNode* topNode,
                                   const std::vector<std::string>& activeTrig)
{
  // centrality bin index (Au+Au only, else -1)
  const int centIdx = (m_isAuAu ? findCentBin(m_centBin) : -1);

    // Prefer photons (full shower shapes) if available
    if (m_photons)
    {
      PhotonClusterContainer::ConstRange prange = m_photons->getClusters();
      for (auto pit = prange.first; pit != prange.second; ++pit)
      {
        const auto* pho = dynamic_cast<const PhotonClusterv1*>(pit->second);
        if (!pho) continue;

        // If the underlying object is also a RawCluster, grab it
        const RawCluster* rc = dynamic_cast<const RawCluster*>(pit->second);

        // Eta: prefer builder-provided; else derive (if RawCluster available)
        double eta = pho->get_shower_shape_parameter("cluster_eta");
        if (!std::isfinite(eta) || eta == 0.0)
        {
            if (rc) eta = RawClusterUtility::GetPseudorapidity(*rc, CLHEP::Hep3Vector(0,0,m_vz));
            else     continue; // cannot determine eta → skip
        }

        // Global fiducial (PPG-12): |eta^gamma| < m_etaAbsMax (default 0.7)
        if (!std::isfinite(eta) || std::fabs(eta) >= m_etaAbsMax) continue;

        const double et = pho->get_energy() / std::cosh(eta);
        const int etIdx = findEtBin(et);
        if (etIdx < 0) continue;


        // Build SSVars
        const SSVars v = makeSSFromPhoton(pho, et);

        // Fill ID/ISO counters only if we can evaluate isolation (need RawCluster*)
        if (rc)
          for (const auto& trigShort : activeTrig)
            fillIsoSSTagCounters(trigShort, rc, v, et, centIdx, topNode);

        // xJ requires a valid φ; use builder-provided value only
        double phi_gamma = pho->get_shower_shape_parameter("cluster_phi");
        if (!std::isfinite(phi_gamma)) continue;

        const bool iso = (rc ? isIsolated(rc, et, topNode) : false);
        const TightTag tightTag = classifyPhotonTightness(v);

        if (iso && tightTag == TightTag::kTight)
        {
          JetContainer* jets = nullptr;
          if (auto it = m_jets.find("r02"); it != m_jets.end()) jets = it->second;

          if (jets)
          {
            double bestPt = -1.0;
            for (const Jet* j : *jets)
            {
              if (!j) continue;
              const double dphi = std::fabs(TVector2::Phi_mpi_pi(j->get_phi() - phi_gamma));
              if (dphi < m_minBackToBack) continue;   // back-to-back (default 7π/8)
              const double pt = j->get_pt();
              if (pt < m_minJetPt) continue;
              if (pt > bestPt) bestPt = pt;
            }

            if (bestPt > 0.0)
            {
              const double xJ = bestPt / et;
              for (const auto& trigShort : activeTrig)
              {
                auto* hx = getOrBookXJHist(trigShort, etIdx, (m_isAuAu ? centIdx : -1));
                hx->Fill(xJ);
                bumpHistFill(trigShort, std::string("h_xJ") + suffixForBins(etIdx, (m_isAuAu ? centIdx : -1)));
              }
            }
          }
        }
      }
      return;
    }


  // Fallback: RawCluster – isolation only (no SS windows available)
  if (m_clus)
  {
    RawClusterContainer::ConstRange range = m_clus->getClusters();
    for (auto it = range.first; it != range.second; ++it)
    {
      const RawCluster* clus = it->second;
      if (!clus) continue;

      const double eta = RawClusterUtility::GetPseudorapidity(*clus, CLHEP::Hep3Vector(0,0,m_vz));
      if (!std::isfinite(eta) || std::fabs(eta) >= m_etaAbsMax) continue;

      const double et  = clus->get_energy() / std::cosh(eta);

      const int etIdx = findEtBin(et);
      if (etIdx < 0) continue;


      const bool iso = isIsolated(clus, et, topNode);

      for (const auto& trigShort : activeTrig)
      {
        if (iso)
          getOrBookCountHist(trigShort, "h_isolatedCount_ET",    etIdx, (m_isAuAu ? centIdx : -1))->Fill(1);
        else
          getOrBookCountHist(trigShort, "h_nonIsolatedCount_ET", etIdx, (m_isAuAu ? centIdx : -1))->Fill(1);
      }
    }
  }
}



bool RecoilJets::passesPhotonPreselection(const SSVars& v)
{
  const bool pass_e11e33 = (v.e11_over_e33 < PRE_E11E33_MAX);
  const bool pass_et1    = in_open_interval(v.et1, PRE_ET1_MIN, PRE_ET1_MAX);
  const bool pass_e32e35 = in_open_interval(v.e32_over_e35, PRE_E32E35_MIN, PRE_E32E35_MAX);
  const bool pass_weta   = (v.weta_cogx < PRE_WETA_MAX);
  return pass_e11e33 && pass_et1 && pass_e32e35 && pass_weta;
}

RecoilJets::TightTag RecoilJets::classifyPhotonTightness(const SSVars& v)
{
  if (!passesPhotonPreselection(v)) return TightTag::kPreselectionFail;

  const double w_hi = tight_w_hi(v.et_gamma);
  const bool pass_weta   = in_open_interval(v.weta_cogx,  TIGHT_W_LO, w_hi);
  const bool pass_wphi   = in_open_interval(v.wphi_cogx,  TIGHT_W_LO, w_hi);
  const bool pass_e11e33 = in_open_interval(v.e11_over_e33, TIGHT_E11E33_MIN, TIGHT_E11E33_MAX);
  const bool pass_et1    = in_open_interval(v.et1,          TIGHT_ET1_MIN,    TIGHT_ET1_MAX);
  const bool pass_e32e35 = in_open_interval(v.e32_over_e35, TIGHT_E32E35_MIN, TIGHT_E32E35_MAX);

  int n_fail = (!pass_weta) + (!pass_wphi) + (!pass_e11e33) + (!pass_et1) + (!pass_e32e35);
  if (n_fail == 0) return TightTag::kTight;
  if (n_fail >= 2) return TightTag::kNonTight;
  return TightTag::kNeither;
}


void RecoilJets::setIsolationWP(double aGeV, double bPerGeV,
                                double sideGapGeV, double coneR, double towerMin)
{
  m_isoA      = aGeV;
  m_isoB      = bPerGeV;
  m_isoGap    = sideGapGeV;
  m_isoConeR  = std::max(0.05, coneR);
  m_isoTowMin = std::max(0.0, towerMin);
}

double RecoilJets::eiso(const RawCluster* clus, PHCompositeNode* /*topNode*/) const
{
  if (!clus) return 0.0;

  // ClusterIso stores isolation keyed by cone size in tenths (R=0.3 -> 3)
  int cone10 = static_cast<int>(std::lround(10.0 * m_isoConeR));
  if (cone10 < 1) cone10 = 1;

  // Prefer UE-subtracted isolation if present; fall back to unsubtracted.
  // The 3rd flag corresponds to the TowerInfo-based calculation.
  double iso = clus->get_et_iso(cone10, /*subtracted=*/true,  /*towerinfo=*/true);
  if (!std::isfinite(iso) || iso <= 0.0)
    iso = clus->get_et_iso(cone10, /*subtracted=*/false, /*towerinfo=*/true);

  return std::isfinite(iso) ? iso : 0.0;
}


bool RecoilJets::isIsolated(const RawCluster* clus, double et_gamma, PHCompositeNode* topNode) const
{
  const double thr  = m_isoA + m_isoB * et_gamma;
  const double eiso = this->eiso(clus, topNode);
  return (eiso < thr);
}

bool RecoilJets::isNonIsolated(const RawCluster* clus, double et_gamma, PHCompositeNode* topNode) const
{
  const double thr  = m_isoA + m_isoB * et_gamma + m_isoGap;
  const double eiso = this->eiso(clus, topNode);
  return (eiso >= thr);
}

// ---------- E_T / centrality bin helpers ----------
int RecoilJets::findEtBin(double et) const
{
  if (m_gammaEtBins.size() < 2) return -1;
  for (size_t i=0; i+1<m_gammaEtBins.size(); ++i)
    if (et >= m_gammaEtBins[i] && et < m_gammaEtBins[i+1]) return static_cast<int>(i);
  return -1;
}

int RecoilJets::findCentBin(int cent) const
{
  if (m_centEdges.size() < 2) return -1;
  for (size_t i=0; i+1<m_centEdges.size(); ++i)
    if (cent >= m_centEdges[i] && cent < m_centEdges[i+1]) return static_cast<int>(i);
  return -1;
}

// suffix: _ET_lo_hi  [ _cent_clo_chi only if isAuAu & centIdx>=0 ]
std::string RecoilJets::suffixForBins(int etIdx, int centIdx) const
{
  std::ostringstream s;
  if (etIdx >= 0) {
    const double lo = m_gammaEtBins[etIdx];
    const double hi = m_gammaEtBins[etIdx+1];
    s << "_ET_" << std::fixed << std::setprecision(0) << lo << '_' << hi;
  }
  if (m_isAuAu && centIdx >= 0) {
    const int clo = m_centEdges[centIdx];
    const int chi = m_centEdges[centIdx+1];
    s << "_cent_" << clo << '_' << chi;
  }
  return s.str();
}

TH1I* RecoilJets::getOrBookCountHist(const std::string& trig,
                                     const std::string& base,
                                     int etIdx, int centIdx)
{
  const std::string name = base + suffixForBins(etIdx, centIdx);
  auto& H = qaHistogramsByTrigger[trig];
  if (auto it = H.find(name); it != H.end())
    return dynamic_cast<TH1I*>(it->second);

  TDirectory* cur = gDirectory;
  TDirectory* dir = out->GetDirectory(trig.c_str());
  if (!dir) dir = out->mkdir(trig.c_str());
  dir->cd();

  auto* h = new TH1I(name.c_str(), (name+";count;entries").c_str(), 1, 0.5, 1.5);
  h->GetXaxis()->SetBinLabel(1,"count");
  H[name] = h;

  if (cur) cur->cd();
  return h;
}

TH1F* RecoilJets::getOrBookXJHist(const std::string& trig, int etIdx, int centIdx)
{
  const std::string base = "h_xJ";
  const std::string name = base + suffixForBins(etIdx, centIdx); // h_xJ_ET_lo_hi[_cent_clo_chi]

  auto& H = qaHistogramsByTrigger[trig];
  if (auto it = H.find(name); it != H.end())
    return dynamic_cast<TH1F*>(it->second);

  TDirectory* cur = gDirectory;
  TDirectory* dir = out->GetDirectory(trig.c_str());
  if (!dir) dir = out->mkdir(trig.c_str());
  dir->cd();

  auto* h = new TH1F(name.c_str(),
                     (name+";x_{J}=p_{T}^{jet}/E_{T}^{#gamma};Entries").c_str(),
                     60, 0.0, 3.0);
  H[name] = h;

  if (cur) cur->cd();
  return h;
}

TH1F* RecoilJets::getOrBookIsoHist(const std::string& trig, int etIdx, int centIdx)
{
  const std::string base = "h_Eiso";
  const std::string name = base + suffixForBins(etIdx, centIdx);

  auto& H = qaHistogramsByTrigger[trig];
  if (auto it = H.find(name); it != H.end())
    return dynamic_cast<TH1F*>(it->second);

  TDirectory* cur = gDirectory;
  TDirectory* dir = out->GetDirectory(trig.c_str());
  if (!dir) dir = out->mkdir(trig.c_str());
  dir->cd();

  auto* h = new TH1F(name.c_str(),
                     (name+";E_{T}^{iso} [GeV];Entries").c_str(),
                     120, 0.0, 12.0);
  H[name] = h;

  if (cur) cur->cd();
  return h;
}

// Record a histogram fill (for end-of-job diagnostics)
void RecoilJets::bumpHistFill(const std::string& trig, const std::string& hnameWithSuffix)
{
  m_histFill[trig + "::" + hnameWithSuffix] += 1;
}

// Book a shower-shape histogram for a given variable and tag
TH1F* RecoilJets::getOrBookSSHist(const std::string& trig,
                                  const std::string& varKey,
                                  const std::string& tagKey,
                                  int etIdx, int centIdx)
{
  const std::string base = "h_ss_" + varKey + "_" + tagKey;
  const std::string name = base + suffixForBins(etIdx, centIdx);

  auto& H = qaHistogramsByTrigger[trig];
  if (auto it = H.find(name); it != H.end())
    return dynamic_cast<TH1F*>(it->second);

  // Choose safe, generic ranges
  int    nb = 120;
  double lo = 0.0, hi = 1.2;
  if (varKey == "weta" || varKey == "wphi") { nb = 120; lo = 0.0; hi = 1.2; }
  else if (varKey == "et1")                 { nb = 120; lo = 0.0; hi = 1.2; }
  else if (varKey == "e11e33" || varKey == "e32e35") { nb = 120; lo = 0.0; hi = 1.2; }

  TDirectory* cur = gDirectory;
  TDirectory* dir = out->GetDirectory(trig.c_str());
  if (!dir) dir = out->mkdir(trig.c_str());
  dir->cd();

  const std::string title = base + ";"+varKey+";Entries";
  auto* h = new TH1F(name.c_str(), title.c_str(), nb, lo, hi);
  H[name] = h;

  if (cur) cur->cd();
  return h;
}




void RecoilJets::fillIsoSSTagCounters(const std::string& trig,
                                      const RawCluster* clus,
                                      const SSVars& v,
                                      double et_gamma,
                                      int centIdx,
                                      PHCompositeNode* topNode)
{
  const int etIdx = findEtBin(et_gamma);
  if (etIdx < 0 || !clus) return;

  // Slice key (used for stats & histogram names)
  const std::string slice = suffixForBins(etIdx, (m_isAuAu ? centIdx : -1));

  // Compute iso and record the Eiso distribution
  const double eiso_et = eiso(clus, topNode);
  {
    auto* h = getOrBookIsoHist(trig, etIdx, (m_isAuAu ? centIdx : -1));
    h->Fill(eiso_et);
    bumpHistFill(trig, "h_Eiso" + slice);
  }

  // Categorize per selection
  const bool iso = (eiso_et < (m_isoA + m_isoB * et_gamma));
  const TightTag tag = classifyPhotonTightness(v);

  // Update per-slice counters
  auto& S = m_catByTrig[trig][slice];
  S.seen += 1;
  (tag == TightTag::kTight ? S.tight : S.nonTight) += 1;
  (iso ? S.isoPass : S.isoFail) += 1;

  // Fill per-variable shower-shape histos for tight/nonTight only
  if (tag == TightTag::kTight || tag == TightTag::kNonTight)
  {
    const std::string tagKey = (tag == TightTag::kTight) ? "tight" : "nontight";
    auto fillSS = [&](const std::string& key, double val)
    {
      auto* h = getOrBookSSHist(trig, key, tagKey, etIdx, (m_isAuAu ? centIdx : -1));
      h->Fill(val);
      bumpHistFill(trig, "h_ss_" + key + "_" + tagKey + slice);
    };
    fillSS("weta",   v.weta_cogx);
    fillSS("wphi",   v.wphi_cogx);
    fillSS("et1",    v.et1);
    fillSS("e11e33", v.e11_over_e33);
    fillSS("e32e35", v.e32_over_e35);
  }

  // Fill (and track) the marginal counters
  if (tag == TightTag::kTight)
  {
    auto* h = getOrBookCountHist(trig, "h_tightCount_ET", etIdx, (m_isAuAu ? centIdx : -1));
    h->Fill(1);
    bumpHistFill(trig, "h_tightCount_ET" + slice);
  }
  else
  {
    auto* h = getOrBookCountHist(trig, "h_nonTightCount_ET", etIdx, (m_isAuAu ? centIdx : -1));
    h->Fill(1);
    bumpHistFill(trig, "h_nonTightCount_ET" + slice);
  }

  if (iso)
  {
    auto* h = getOrBookCountHist(trig, "h_isolatedCount_ET", etIdx, (m_isAuAu ? centIdx : -1));
    h->Fill(1);
    bumpHistFill(trig, "h_isolatedCount_ET" + slice);
  }
  else
  {
    auto* h = getOrBookCountHist(trig, "h_nonIsolatedCount_ET", etIdx, (m_isAuAu ? centIdx : -1));
    h->Fill(1);
    bumpHistFill(trig, "h_nonIsolatedCount_ET" + slice);
  }

  // ID-sideband (fail ≥2 tight cuts)
  if (tag == TightTag::kNonTight)
  {
    S.idSB_total += 1;
    {
      auto* h = getOrBookCountHist(trig, "h_idSB_total", etIdx, (m_isAuAu ? centIdx : -1));
      h->Fill(1);
      bumpHistFill(trig, "h_idSB_total" + slice);
    }
    if (iso)
    {
      S.idSB_pass += 1;
      auto* h = getOrBookCountHist(trig, "h_idSB_pass", etIdx, (m_isAuAu ? centIdx : -1));
      h->Fill(1);
      bumpHistFill(trig, "h_idSB_pass" + slice);
    }
  }
}
