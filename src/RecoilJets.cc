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
#include <calobase/RawTowerGeomContainer.h>           // <-- required base
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/RawClusterUtility.h>
#include <mbd/MbdPmtHit.h>
#include <mbd/MbdGeom.h>
#include <mbd/MbdOut.h>
#include <ffarawobjects/Gl1Packet.h>
#include <mbd/MbdPmtContainer.h>
#include <centrality/CentralityInfo.h>
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

bool RecoilJets::fetchNodes(PHCompositeNode* top)
{
  /* ––– primary vertex –––––––––––––––––––––––––––––––––––––––––––––––– */
  GlobalVertexMap* vmap = findNode::getClass<GlobalVertexMap>(top, "GlobalVertexMap");
  m_vtx = nullptr; m_vx = m_vy = m_vz = 0.;

  if (!vmap)         { LOG(1, CLR_YELLOW, "  – GlobalVertexMap node **missing** → skip event"); return false; }
  if (vmap->empty()) { LOG(2, CLR_YELLOW, "  – GlobalVertexMap is **empty** → skip event");     return false; }

  m_vtx = vmap->begin()->second;
  if (!m_vtx)        { LOG(1, CLR_YELLOW, "  – vertex pointer null → skip event");              return false; }

  m_vx = m_vtx->get_x();
  m_vy = m_vtx->get_y();
  m_vz = m_vtx->get_z();

  LOG(4, CLR_BLUE, "    [fetchNodes] dataset=" << (m_isAuAu ? "Au+Au" : "p+p")
                                              << "  vz=" << std::fixed << std::setprecision(2) << m_vz);

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
    {
      LOG(2, CLR_YELLOW, "  – missing " << lbl << " nodes → skip");
      return false;
    }
    m_calo[lbl] = { tw, ge, 0. };
    LOG(5, CLR_BLUE,  "    [fetchNodes] towers OK: " << lbl << "  node=" << node << "  geom=" << geo);
  }

  /* ––– clusters & photons ––––––––––––––––––––––––––––––––––––––––––– */
  m_clus    = findNode::getClass<RawClusterContainer>(top, "CLUSTERINFO_CEMC");
  m_photons = findNode::getClass<PhotonClusterContainer>(top, "PHOTONCLUSTER_CEMC");

  if (m_clus)
  {
    // RawClusterContainer does not expose size() directly in older interfaces,
    // so just print presence here.
    LOG(5, CLR_BLUE, "    [fetchNodes] CLUSTERINFO_CEMC present");
  }
  else
  {
    LOG(3, CLR_YELLOW, "    [fetchNodes] CLUSTERINFO_CEMC **missing**");
  }

  LOG(5, CLR_BLUE, std::string("    [fetchNodes] PHOTONCLUSTER_CEMC ") + (m_photons ? "present" : "missing"));

  /* ––– jet containers (dataset-aware) ––––––––––––––––––––––––––––––– */
  m_jets.clear();
  for (const auto& jnm : kJetRadii)
  {
    const std::string rKey = jnm.key;                              // "r02"
    const std::string node = m_isAuAu ? jnm.aa_node : jnm.pp_node; // dataset-aware name

    auto* jc = findNode::getClass<JetContainer>(top, node);
    if (!jc)
    {
      LOG(2, CLR_YELLOW, "  – jet node missing: " << node << " (xJ disabled for " << rKey << ")");
    }
    else
    {
      // Print current jet count in the container for this event
      LOG(4, CLR_GREEN, "    [fetchNodes] jet node found: " << node
                        << "  jets=" << jc->size());
    }
    m_jets[rKey] = jc; // may be nullptr; checked later
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
    // ---------- Au+Au path: use scaled GL1 bits only (no MB pre-gate) ----------
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
  // -------- guards & context ------------------------------------------------
  if (!topNode)
  {
    LOG(0, CLR_YELLOW, "  [processCandidates] topNode == nullptr – aborting this event");
    return;
  }

  if (activeTrig.empty())
  {
    LOG(0, CLR_YELLOW, "  [processCandidates] activeTrig is EMPTY – nothing to fill this event");
    return;
  }

  // Centrality bin index (Au+Au only, else -1)
  const int centIdx = (m_isAuAu ? findCentBin(m_centBin) : -1);
  LOG(4, CLR_BLUE, "  [processCandidates] dataset=" << (m_isAuAu ? "Au+Au" : "p+p")
                                                   << "  centIdx=" << centIdx
                                                   << "  m_centBin=" << m_centBin
                                                   << "  nTriggers=" << activeTrig.size());

  // If neither photons nor clusters are present, we cannot do anything
  if (!m_photons && !m_clus)
  {
    LOG(2, CLR_YELLOW, "  [processCandidates] PHOTONCLUSTER_CEMC and CLUSTERINFO_CEMC both MISSING – nothing to process");
    return;
  }

  // -------- PHOTON path (preferred; full shower shapes) ---------------------
  if (m_photons)
  {
    PhotonClusterContainer::ConstRange prange = m_photons->getClusters();

    // Quick emptiness check
    std::size_t nPho = 0;
    for (auto tmp = prange.first; tmp != prange.second; ++tmp) ++nPho;
    LOG(4, CLR_BLUE, "    [processCandidates] photons present: " << nPho);
    if (nPho == 0)
      LOG(3, CLR_YELLOW, "    [processCandidates] photon container is EMPTY for this event");

    std::size_t nUsed  = 0;
    std::size_t nSkipEta = 0, nSkipEtBin = 0, nSkipPhi = 0, nNoRC = 0, nNotIso = 0, nNotTight = 0;

    try
    {
      int iPho = 0;
      for (auto pit = prange.first; pit != prange.second; ++pit, ++iPho)
      {
        // Concrete type check
        const auto* pho = dynamic_cast<const PhotonClusterv1*>(pit->second);
        if (!pho)
        {
          LOG(5, CLR_YELLOW, "      [pho#" << iPho << "] pointer is not PhotonClusterv1 – skipping");
          continue;
        }

        // If the underlying object is also a RawCluster, grab it (needed for ISO)
        const RawCluster* rc = dynamic_cast<const RawCluster*>(pit->second);

        // Eta: prefer builder-provided; else derive (if RawCluster available)
        double eta = pho->get_shower_shape_parameter("cluster_eta");
        if (!std::isfinite(eta) || eta == 0.0)
        {
          if (rc) eta = RawClusterUtility::GetPseudorapidity(*rc, CLHEP::Hep3Vector(0,0,m_vz));
          else
          {
            ++nSkipEta;
            LOG(5, CLR_YELLOW, "      [pho#" << iPho << "] cannot determine eta (no rc and builder eta invalid) – skip");
            continue;
          }
        }

        // Global fiducial (PPG-12): |eta^gamma| < m_etaAbsMax (default 0.7)
        if (!std::isfinite(eta) || std::fabs(eta) >= m_etaAbsMax)
        {
          ++nSkipEta;
          if (Verbosity() >= 5)
            LOG(5, CLR_YELLOW, "      [pho#" << iPho << "] |eta|=" << std::fabs(eta)
                                             << " outside fiducial (max=" << m_etaAbsMax << ") – skip");
          continue;
        }

        const double energy = pho->get_energy();
        if (!std::isfinite(energy) || energy <= 0.0)
        {
          ++nSkipEtBin;
          LOG(5, CLR_YELLOW, "      [pho#" << iPho << "] non-positive or NaN energy=" << energy << " – skip");
          continue;
        }

        const double et = energy / std::cosh(eta);
        if (!std::isfinite(et) || et <= 0.0)
        {
          ++nSkipEtBin;
          LOG(5, CLR_YELLOW, "      [pho#" << iPho << "] invalid ET=" << et << " (E=" << energy << ", eta=" << eta << ") – skip");
          continue;
        }

        const int etIdx = findEtBin(et);
        if (etIdx < 0)
        {
          ++nSkipEtBin;
          if (Verbosity() >= 5)
            LOG(5, CLR_YELLOW, "      [pho#" << iPho << "] ET=" << et << " not in configured ET bins – skip");
          continue;
        }

        // Build SSVars (basic validation only; keep logic identical)
        const SSVars v = makeSSFromPhoton(pho, et);

        // Fill ID/ISO counters only if we can evaluate isolation (need RawCluster*)
        if (rc)
        {
          for (const auto& trigShort : activeTrig)
            fillIsoSSTagCounters(trigShort, rc, v, et, centIdx, topNode);
        }
        else
        {
          ++nNoRC;
          if (Verbosity() >= 5)
            LOG(5, CLR_YELLOW, "      [pho#" << iPho << "] RawCluster pointer is NULL – cannot compute isolation counters");
        }

        // xJ requires a valid φ (builder-provided only)
        double phi_gamma = pho->get_shower_shape_parameter("cluster_phi");
        if (!std::isfinite(phi_gamma))
        {
          ++nSkipPhi;
          if (Verbosity() >= 5)
            LOG(5, CLR_YELLOW, "      [pho#" << iPho << "] cluster_phi invalid – skip xJ step");
          continue;
        }

        const bool iso = (rc ? isIsolated(rc, et, topNode) : false);
        const TightTag tightTag = classifyPhotonTightness(v);

        if (!(iso && tightTag == TightTag::kTight))
        {
          if (!iso)      ++nNotIso;
          if (tightTag != TightTag::kTight) ++nNotTight;
          if (Verbosity() >= 6)
            LOG(6, CLR_BLUE, "      [pho#" << iPho << "] NOT used for xJ (iso=" << iso
                                           << ", tightTag=" << static_cast<int>(tightTag) << ")");
          continue;
        }

        // Dataset-aware r02 jet container
        JetContainer* jets = nullptr;
        if (auto it = m_jets.find("r02"); it != m_jets.end()) jets = it->second;

        if (!jets)
        {
          LOG(3, CLR_YELLOW, "      [pho#" << iPho << "] jet container r02 is nullptr – cannot form xJ");
          continue;
        }

        // Debug: show how many jets and print the first few
        const std::size_t nj = jets->size();
        LOG(4, CLR_BLUE, "      [pho#" << iPho << "] r02 jet container size = " << nj);
        if (Verbosity() >= 5 && nj > 0)
        {
          int shown = 0;
          for (const Jet* j : *jets)
          {
            if (!j) continue;
            if (shown < 5)
            {
              LOG(5, CLR_CYAN, "        jet " << shown
                    << "  pt="  << std::fixed << std::setprecision(2) << j->get_pt()
                    << "  eta=" << std::setprecision(3) << j->get_eta()
                    << "  phi=" << std::setprecision(3) << j->get_phi());
              ++shown;
            }
            else break;
          }
        }

        // Scan jets for away-side match
        int nPassDphi = 0, nPassPt = 0;
        double bestPt = -1.0;
        for (const Jet* j : *jets)
        {
          if (!j) continue;
          const double dphi = std::fabs(TVector2::Phi_mpi_pi(j->get_phi() - phi_gamma));
          if (dphi >= m_minBackToBack) ++nPassDphi; else continue;

          const double pt = j->get_pt();
          if (pt >= m_minJetPt)        ++nPassPt;   else continue;

          if (pt > bestPt) bestPt = pt;
        }

        if (bestPt > 0.0)
        {
          const double xJ = bestPt / et;
          if (Verbosity() >= 5)
            LOG(5, CLR_GREEN, "      [pho#" << iPho << "] xJ = " << std::fixed << std::setprecision(3)
                                            << xJ << "  (bestPt=" << bestPt
                                            << ", ET^γ=" << et
                                            << ", passed dphi=" << nPassDphi
                                            << ", passed pt=" << nPassPt << ")");

          for (const auto& trigShort : activeTrig)
          {
            auto* hx = getOrBookXJHist(trigShort, etIdx, (m_isAuAu ? centIdx : -1));
            hx->Fill(xJ);
            bumpHistFill(trigShort, std::string("h_xJ") + suffixForBins(etIdx, (m_isAuAu ? centIdx : -1)));
          }
          ++nUsed;
        }
        else
        {
          LOG(4, CLR_YELLOW, "      [pho#" << iPho << "] no away-side jet passes Δφ/pt cuts (bestPt ≤ 0)"
                                           << "  [passed dphi=" << nPassDphi
                                           << ", passed pt=" << nPassPt
                                           << ", minBackToBack=" << m_minBackToBack
                                           << ", minJetPt=" << m_minJetPt << "]");
        }
      } // photon loop
    }
    catch (const std::exception& e)
    {
      LOG(0, CLR_YELLOW, "    [processCandidates] EXCEPTION in photon path: " << e.what());
    }
    catch (...)
    {
      LOG(0, CLR_YELLOW, "    [processCandidates] UNKNOWN exception in photon path");
    }

    // Summary for photon path
    if (Verbosity() >= 4)
      LOG(4, CLR_BLUE, "    [processCandidates] photon summary: used=" << nUsed
                                  << "  skipEta=" << nSkipEta
                                  << "  skipEtBin=" << nSkipEtBin
                                  << "  skipPhi=" << nSkipPhi
                                  << "  noRC=" << nNoRC
                                  << "  notIso=" << nNotIso
                                  << "  notTight=" << nNotTight);

    return; // keep original behavior: prefer photon path and return
  }

  // -------- FALLBACK: RawCluster – isolation only (no SS windows) ----------
  if (m_clus)
  {
    RawClusterContainer::ConstRange range = m_clus->getClusters();

    // Emptiness check
    std::size_t nCl = 0;
    for (auto tmp = range.first; tmp != range.second; ++tmp) ++nCl;
    LOG(4, CLR_BLUE, "    [processCandidates] clusters present: " << nCl);
    if (nCl == 0)
      LOG(3, CLR_YELLOW, "    [processCandidates] cluster container is EMPTY for this event");

    std::size_t nIso = 0, nNonIso = 0, nSkipEta = 0, nSkipEtBin = 0;

    try
    {
      int iCl = 0;
      for (auto it = range.first; it != range.second; ++it, ++iCl)
      {
        const RawCluster* clus = it->second;
        if (!clus)
        {
          LOG(5, CLR_YELLOW, "      [clus#" << iCl << "] null pointer – skip");
          continue;
        }

        const double eta = RawClusterUtility::GetPseudorapidity(*clus, CLHEP::Hep3Vector(0,0,m_vz));
        if (!std::isfinite(eta) || std::fabs(eta) >= m_etaAbsMax)
        {
          ++nSkipEta;
          if (Verbosity() >= 5)
            LOG(5, CLR_YELLOW, "      [clus#" << iCl << "] |eta|=" << std::fabs(eta)
                                              << " outside fiducial (max=" << m_etaAbsMax << ") – skip");
          continue;
        }

        const double energy = clus->get_energy();
        if (!std::isfinite(energy) || energy <= 0.0)
        {
          ++nSkipEtBin;
          LOG(5, CLR_YELLOW, "      [clus#" << iCl << "] non-positive or NaN energy=" << energy << " – skip");
          continue;
        }

        const double et  = energy / std::cosh(eta);
        if (!std::isfinite(et) || et <= 0.0)
        {
          ++nSkipEtBin;
          LOG(5, CLR_YELLOW, "      [clus#" << iCl << "] invalid ET=" << et << " – skip");
          continue;
        }

        const int etIdx = findEtBin(et);
        if (etIdx < 0)
        {
          ++nSkipEtBin;
          if (Verbosity() >= 5)
            LOG(5, CLR_YELLOW, "      [clus#" << iCl << "] ET=" << et << " not in configured ET bins – skip");
          continue;
        }

        const bool iso = isIsolated(clus, et, topNode);
        for (const auto& trigShort : activeTrig)
        {
          if (iso)
          {
            getOrBookCountHist(trigShort, "h_isolatedCount_ET",    etIdx, (m_isAuAu ? centIdx : -1))->Fill(1);
          }
          else
          {
            getOrBookCountHist(trigShort, "h_nonIsolatedCount_ET", etIdx, (m_isAuAu ? centIdx : -1))->Fill(1);
          }
        }
        (iso ? ++nIso : ++nNonIso);
      } // cluster loop
    }
    catch (const std::exception& e)
    {
      LOG(0, CLR_YELLOW, "    [processCandidates] EXCEPTION in cluster fallback: " << e.what());
    }
    catch (...)
    {
      LOG(0, CLR_YELLOW, "    [processCandidates] UNKNOWN exception in cluster fallback");
    }

    if (Verbosity() >= 4)
      LOG(4, CLR_BLUE, "    [processCandidates] cluster summary: iso=" << nIso
                                  << "  nonIso=" << nNonIso
                                  << "  skipEta=" << nSkipEta
                                  << "  skipEtBin=" << nSkipEtBin);
  }
}





bool RecoilJets::passesPhotonPreselection(const SSVars& v)
{
  // Basic sanity on inputs
  const bool ok_vals =
      std::isfinite(v.weta_cogx) &&
      std::isfinite(v.wphi_cogx) &&
      std::isfinite(v.et1) &&
      std::isfinite(v.e11_over_e33) &&
      std::isfinite(v.e32_over_e35) &&
      std::isfinite(v.et_gamma);

  if (!ok_vals)
  {
    LOG(2, CLR_YELLOW,
        "  [passesPhotonPreselection] non-finite SSVars detected: "
        << "weta=" << v.weta_cogx << " wphi=" << v.wphi_cogx
        << " et1=" << v.et1 << " e11/e33=" << v.e11_over_e33
        << " e32/e35=" << v.e32_over_e35 << " ET^γ=" << v.et_gamma);
    return false;
  }

  const bool pass_e11e33 = (v.e11_over_e33 < PRE_E11E33_MAX);
  const bool pass_et1    = in_open_interval(v.et1, PRE_ET1_MIN, PRE_ET1_MAX);
  const bool pass_e32e35 = in_open_interval(v.e32_over_e35, PRE_E32E35_MIN, PRE_E32E35_MAX);
  const bool pass_weta   = (v.weta_cogx < PRE_WETA_MAX);

  if (Verbosity() >= 5)
  {
    LOG(5, CLR_BLUE,
        "  [passesPhotonPreselection] "
        << "weta=" << v.weta_cogx << " (<" << PRE_WETA_MAX << ") → " << pass_weta
        << " | et1=" << v.et1 << " ∈ (" << PRE_ET1_MIN << "," << PRE_ET1_MAX << ") → " << pass_et1
        << " | e11/e33=" << v.e11_over_e33 << " (<" << PRE_E11E33_MAX << ") → " << pass_e11e33
        << " | e32/e35=" << v.e32_over_e35 << " ∈ (" << PRE_E32E35_MIN << "," << PRE_E32E35_MAX << ") → " << pass_e32e35);
  }

  const bool pass_all = pass_e11e33 && pass_et1 && pass_e32e35 && pass_weta;

  if (!pass_all && Verbosity() >= 4)
  {
    LOG(4, CLR_YELLOW,
        "  [passesPhotonPreselection] FAILED: weta=" << v.weta_cogx
        << ", et1=" << v.et1
        << ", e11/e33=" << v.e11_over_e33
        << ", e32/e35=" << v.e32_over_e35);
  }

  return pass_all;
}


RecoilJets::TightTag RecoilJets::classifyPhotonTightness(const SSVars& v)
{
  // If preselection fails, short-circuit
  if (!passesPhotonPreselection(v))
  {
    if (Verbosity() >= 5)
      LOG(5, CLR_YELLOW, "  [classifyPhotonTightness] Preselection FAILED – returning kPreselectionFail");
    return TightTag::kPreselectionFail;
  }

  // Upper width cut is ET-dependent
  const double w_hi = tight_w_hi(v.et_gamma);
  if (!std::isfinite(w_hi))
  {
    LOG(2, CLR_YELLOW, "  [classifyPhotonTightness] non-finite tight_w_hi for ET^γ=" << v.et_gamma
                        << " – treating as failure");
    return TightTag::kNeither;
  }

  const bool pass_weta   = in_open_interval(v.weta_cogx,  TIGHT_W_LO, w_hi);
  const bool pass_wphi   = in_open_interval(v.wphi_cogx,  TIGHT_W_LO, w_hi);
  const bool pass_e11e33 = in_open_interval(v.e11_over_e33, TIGHT_E11E33_MIN, TIGHT_E11E33_MAX);
  const bool pass_et1    = in_open_interval(v.et1,          TIGHT_ET1_MIN,    TIGHT_ET1_MAX);
  const bool pass_e32e35 = in_open_interval(v.e32_over_e35, TIGHT_E32E35_MIN, TIGHT_E32E35_MAX);

  if (Verbosity() >= 5)
  {
    LOG(5, CLR_BLUE,
        "  [classifyPhotonTightness] ET^γ=" << v.et_gamma
        << " → w_hi=" << w_hi
        << " | weta=" << v.weta_cogx << " ok=" << pass_weta
        << " | wphi=" << v.wphi_cogx << " ok=" << pass_wphi
        << " | e11/e33=" << v.e11_over_e33 << " ok=" << pass_e11e33
        << " | et1=" << v.et1 << " ok=" << pass_et1
        << " | e32/e35=" << v.e32_over_e35 << " ok=" << pass_e32e35);
  }

  const int n_fail = (!pass_weta) + (!pass_wphi) + (!pass_e11e33) + (!pass_et1) + (!pass_e32e35);

  TightTag tag;
  if (n_fail == 0)      tag = TightTag::kTight;
  else if (n_fail >= 2) tag = TightTag::kNonTight;
  else                  tag = TightTag::kNeither;

  if (Verbosity() >= 4)
    LOG(4, CLR_GREEN, "  [classifyPhotonTightness] n_fail=" << n_fail
                        << " → tag=" << static_cast<int>(tag));

  return tag;
}


void RecoilJets::setIsolationWP(double aGeV, double bPerGeV,
                                double sideGapGeV, double coneR, double towerMin)
{
  // Keep original behavior (min floors), but report anything suspicious
  if (!std::isfinite(aGeV) || !std::isfinite(bPerGeV) ||
      !std::isfinite(sideGapGeV) || !std::isfinite(coneR) || !std::isfinite(towerMin))
  {
    LOG(2, CLR_YELLOW,
        "  [setIsolationWP] Non-finite input(s): "
        << "A=" << aGeV << " B=" << bPerGeV << " gap=" << sideGapGeV
        << " coneR=" << coneR << " towerMin=" << towerMin);
  }

  m_isoA      = aGeV;
  m_isoB      = bPerGeV;
  m_isoGap    = sideGapGeV;

  const double cone_before  = coneR;
  const double tower_before = towerMin;

  m_isoConeR  = std::max(0.05, coneR);
  m_isoTowMin = std::max(0.0,  towerMin);

  if (Verbosity() >= 4)
  {
    LOG(4, CLR_BLUE,
        "  [setIsolationWP] Isolation WP set:"
        << "  A=" << m_isoA
        << "  B=" << m_isoB
        << "  gap=" << m_isoGap
        << "  coneR=" << m_isoConeR << (cone_before != m_isoConeR ? " (clamped)" : "")
        << "  towerMin=" << m_isoTowMin << (tower_before != m_isoTowMin ? " (clamped)" : ""));
  }
}


double RecoilJets::eiso(const RawCluster* clus, PHCompositeNode* /*topNode*/) const
{
  if (!clus)
  {
    if (Verbosity() >= 3)
      LOG(3, CLR_YELLOW, "  [eiso] clus==nullptr → return 0");
  #ifdef __cpp_lib_isfinite
    return 0.0;
  #else
    return 0.0;
  #endif
  }

  // ClusterIso stores isolation keyed by cone size in tenths (R=0.3 → 3)
  int cone10 = static_cast<int>(std::lround(10.0 * m_isoConeR));
  if (cone10 < 1)
  {
    if (Verbosity() >= 4)
      LOG(4, CLR_YELLOW, "  [eiso] cone size < 0.1? cone10=" << cone10 << " – clamping to 1");
    cone10 = 1;
  }

  // Prefer UE-subtracted; fallback to unsubtracted if invalid
  double iso_sub = clus->get_et_iso(cone10, /*subtracted=*/true,  /*towerinfo=*/true);
  double iso_uns = clus->get_et_iso(cone10, /*subtracted=*/false, /*towerinfo=*/true);

  bool sub_ok = std::isfinite(iso_sub) && (iso_sub > 0.0 || iso_uns <= 0.0);
  double iso  = sub_ok ? iso_sub : iso_uns;

  if (!std::isfinite(iso))
  {
    if (Verbosity() >= 3)
      LOG(3, CLR_YELLOW, "  [eiso] non-finite iso (sub=" << iso_sub << ", uns=" << iso_uns << ") → return 0");
    return 0.0;
  }

  if (Verbosity() >= 5)
  {
    LOG(5, CLR_BLUE, "  [eiso] cone10=" << cone10
            << "  iso_sub=" << iso_sub << "  iso_uns=" << iso_uns
            << "  used=" << iso << (sub_ok ? " (sub)" : " (uns)"));
  }
  return (iso < 0.0 ? 0.0 : iso); // guard tiny negatives
}


bool RecoilJets::isIsolated(const RawCluster* clus, double et_gamma, PHCompositeNode* topNode) const
{
  if (!std::isfinite(et_gamma) || et_gamma <= 0.0)
  {
    if (Verbosity() >= 4)
      LOG(4, CLR_YELLOW, "  [isIsolated] invalid ET^γ=" << et_gamma << " → false");
    return false;
  }

  const double thr  = m_isoA + m_isoB * et_gamma;
  const double eiso_val = this->eiso(clus, topNode);

  if (Verbosity() >= 5)
    LOG(5, CLR_BLUE, "  [isIsolated] ET^γ=" << et_gamma << "  thr=" << thr << "  eiso=" << eiso_val
                        << "  pass=" << (eiso_val < thr));

  return (eiso_val < thr);
}


bool RecoilJets::isNonIsolated(const RawCluster* clus, double et_gamma, PHCompositeNode* topNode) const
{
  if (!std::isfinite(et_gamma) || et_gamma <= 0.0)
  {
    if (Verbosity() >= 4)
      LOG(4, CLR_YELLOW, "  [isNonIsolated] invalid ET^γ=" << et_gamma << " → false");
    return false;
  }

  const double thr  = m_isoA + m_isoB * et_gamma + m_isoGap;
  const double eiso_val = this->eiso(clus, topNode);

  if (Verbosity() >= 5)
    LOG(5, CLR_BLUE, "  [isNonIsolated] ET^γ=" << et_gamma << "  thr=" << thr << "  eiso=" << eiso_val
                        << "  pass=" << (eiso_val >= thr));

  return (eiso_val >= thr);
}


// ---------- E_T / centrality bin helpers ----------
int RecoilJets::findEtBin(double et) const
{
  if (m_gammaEtBins.size() < 2)
  {
    if (Verbosity() >= 3)
      LOG(3, CLR_YELLOW, "  [findEtBin] ET bin edges vector has size < 2 – returning -1");
    return -1;
  }

  // Optional monotonicity check (diagnostic only)
  if (Verbosity() >= 6)
  {
    bool mono = std::is_sorted(m_gammaEtBins.begin(), m_gammaEtBins.end());
    if (!mono)
      LOG(6, CLR_YELLOW, "  [findEtBin] ET bin edges are NOT monotonic");
  }

  if (!std::isfinite(et))
  {
    if (Verbosity() >= 4)
      LOG(4, CLR_YELLOW, "  [findEtBin] non-finite ET=" << et << " – returning -1");
    return -1;
  }

  for (size_t i = 0; i + 1 < m_gammaEtBins.size(); ++i)
  {
    if (et >= m_gammaEtBins[i] && et < m_gammaEtBins[i + 1])
      return static_cast<int>(i);
  }

  if (Verbosity() >= 5)
    LOG(5, CLR_BLUE, "  [findEtBin] ET=" << et << " out of configured range – returning -1");
  return -1;
}


int RecoilJets::findCentBin(int cent) const
{
  if (m_centEdges.size() < 2)
  {
    if (Verbosity() >= 3)
      LOG(3, CLR_YELLOW, "  [findCentBin] centrality edges vector has size < 2 – returning -1");
    return -1;
  }

  // Optional monotonicity check (diagnostic only)
  if (Verbosity() >= 6)
  {
    bool mono = std::is_sorted(m_centEdges.begin(), m_centEdges.end());
    if (!mono)
      LOG(6, CLR_YELLOW, "  [findCentBin] centrality edges are NOT monotonic");
  }

  for (size_t i = 0; i + 1 < m_centEdges.size(); ++i)
  {
    if (cent >= m_centEdges[i] && cent < m_centEdges[i + 1])
      return static_cast<int>(i);
  }

  if (Verbosity() >= 5)
    LOG(5, CLR_BLUE, "  [findCentBin] cent=" << cent << " out of configured range – returning -1");
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
  const std::string suffix = suffixForBins(etIdx, centIdx);
  const std::string name   = base + suffix;

  if (Verbosity() >= 5)
    LOG(5, CLR_BLUE, "  [getOrBookCountHist] trig=\"" << trig << "\" base=\"" << base
           << "\" etIdx=" << etIdx << " centIdx=" << centIdx
           << " → name=\"" << name << "\"");

  if (trig.empty() || base.empty())
  {
    LOG(2, CLR_YELLOW, "  [getOrBookCountHist] empty trig/base (trig=\"" << trig
           << "\", base=\"" << base << "\") – returning nullptr");
    return nullptr;
  }

  // Map slot
  auto& H = qaHistogramsByTrigger[trig];

  // If already booked, verify type and return
  if (auto it = H.find(name); it != H.end())
  {
    if (auto* h = dynamic_cast<TH1I*>(it->second))
    {
      if (Verbosity() >= 6)
        LOG(6, CLR_GREEN, "    [getOrBookCountHist] reusing existing TH1I \"" << name << "\"");
      return h;
    }
    LOG(2, CLR_YELLOW, "    [getOrBookCountHist] name clash: object \"" << name
                       << "\" exists but is not TH1I – replacing it");
    H.erase(it);
  }

  // Safety: output file & directory
  if (!out || !out->IsOpen())
  {
    LOG(1, CLR_YELLOW, "  [getOrBookCountHist] output TFile invalid/null – returning nullptr");
    return nullptr;
  }

  TDirectory* const prevDir = gDirectory;
  TDirectory* dir = out->GetDirectory(trig.c_str());
  if (!dir)
  {
    if (Verbosity() >= 4) LOG(4, CLR_BLUE, "    [getOrBookCountHist] creating directory \"" << trig << '"');
    dir = out->mkdir(trig.c_str());
  }
  if (!dir)
  {
    LOG(1, CLR_YELLOW, "  [getOrBookCountHist] failed to create/access directory \"" << trig << "\"");
    if (prevDir) prevDir->cd();
    return nullptr;
  }

  dir->cd();
  const std::string title = name + ";count;entries";

  if (Verbosity() >= 5)
    LOG(5, CLR_BLUE, "    [getOrBookCountHist] booking TH1I name=\"" << name << "\" title=\"" << title << '"');

  auto* h = new TH1I(name.c_str(), title.c_str(), 1, 0.5, 1.5);
  if (!h)
  {
    LOG(1, CLR_YELLOW, "  [getOrBookCountHist] new TH1I failed for \"" << name << '"');
    if (prevDir) prevDir->cd();
    return nullptr;
  }
  h->GetXaxis()->SetBinLabel(1,"count");
  H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}


TH1F* RecoilJets::getOrBookXJHist(const std::string& trig, int etIdx, int centIdx)
{
  const std::string base   = "h_xJ";
  const std::string suffix = suffixForBins(etIdx, centIdx);
  const std::string name   = base + suffix;

  if (Verbosity() >= 5)
    LOG(5, CLR_BLUE, "  [getOrBookXJHist] trig=\"" << trig
           << "\" etIdx=" << etIdx << " centIdx=" << centIdx
           << " → name=\"" << name << "\"");

  if (trig.empty())
  {
    LOG(2, CLR_YELLOW, "  [getOrBookXJHist] empty trig – returning nullptr");
    return nullptr;
  }

  auto& H = qaHistogramsByTrigger[trig];

  if (auto it = H.find(name); it != H.end())
  {
    if (auto* h = dynamic_cast<TH1F*>(it->second))
    {
      if (Verbosity() >= 6)
        LOG(6, CLR_GREEN, "    [getOrBookXJHist] reusing existing TH1F \"" << name << "\"");
      return h;
    }
    LOG(2, CLR_YELLOW, "    [getOrBookXJHist] name clash: object \"" << name
                       << "\" exists but is not TH1F – replacing it");
    H.erase(it);
  }

  if (!out || !out->IsOpen())
  {
    LOG(1, CLR_YELLOW, "  [getOrBookXJHist] output TFile invalid/null – returning nullptr");
    return nullptr;
  }

  TDirectory* const prevDir = gDirectory;
  TDirectory* dir = out->GetDirectory(trig.c_str());
  if (!dir) dir = out->mkdir(trig.c_str());
  if (!dir)
  {
    LOG(1, CLR_YELLOW, "  [getOrBookXJHist] failed to create/access directory \"" << trig << "\"");
    if (prevDir) prevDir->cd();
    return nullptr;
  }

  dir->cd();

  // Binning chosen for visibility: 60 bins in [0,3]
  const int    nbins = 60;
  const double xmin  = 0.0;
  const double xmax  = 3.0;

  const std::string title = name + ";x_{J}=p_{T}^{jet}/E_{T}^{#gamma};Entries";

  if (Verbosity() >= 5)
    LOG(5, CLR_BLUE, "    [getOrBookXJHist] booking TH1F name=\"" << name
         << "\" bins=" << nbins << " range=[" << xmin << "," << xmax
         << "] title=\"" << title << '"');

  auto* h = new TH1F(name.c_str(), title.c_str(), nbins, xmin, xmax);
  if (!h)
  {
    LOG(1, CLR_YELLOW, "  [getOrBookXJHist] new TH1F failed for \"" << name << '"');
    if (prevDir) prevDir->cd();
    return nullptr;
  }
  H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}


TH1F* RecoilJets::getOrBookIsoHist(const std::string& trig, int etIdx, int centIdx)
{
  const std::string base   = "h_Eiso";
  const std::string suffix = suffixForBins(etIdx, centIdx);
  const std::string name   = base + suffix;

  if (Verbosity() >= 5)
    LOG(5, CLR_BLUE, "  [getOrBookIsoHist] trig=\"" << trig
           << "\" etIdx=" << etIdx << " centIdx=" << centIdx
           << " → name=\"" << name << "\"");

  if (trig.empty())
  {
    LOG(2, CLR_YELLOW, "  [getOrBookIsoHist] empty trig – returning nullptr");
    return nullptr;
  }

  auto& H = qaHistogramsByTrigger[trig];

  if (auto it = H.find(name); it != H.end())
  {
    if (auto* h = dynamic_cast<TH1F*>(it->second))
    {
      if (Verbosity() >= 6)
        LOG(6, CLR_GREEN, "    [getOrBookIsoHist] reusing existing TH1F \"" << name << "\"");
      return h;
    }
    LOG(2, CLR_YELLOW, "    [getOrBookIsoHist] name clash: object \"" << name
                       << "\" exists but is not TH1F – replacing it");
    H.erase(it);
  }

  if (!out || !out->IsOpen())
  {
    LOG(1, CLR_YELLOW, "  [getOrBookIsoHist] output TFile invalid/null – returning nullptr");
    return nullptr;
  }

  TDirectory* const prevDir = gDirectory;
  TDirectory* dir = out->GetDirectory(trig.c_str());
  if (!dir) dir = out->mkdir(trig.c_str());
  if (!dir)
  {
    LOG(1, CLR_YELLOW, "  [getOrBookIsoHist] failed to create/access directory \"" << trig << "\"");
    if (prevDir) prevDir->cd();
    return nullptr;
  }

  dir->cd();

  // Binning chosen to match your previous layout
  const int    nbins = 120;
  const double xmin  = 0.0;
  const double xmax  = 12.0;

  const std::string title = name + ";E_{T}^{iso} [GeV];Entries";

  if (Verbosity() >= 5)
    LOG(5, CLR_BLUE, "    [getOrBookIsoHist] booking TH1F name=\"" << name
         << "\" bins=" << nbins << " range=[" << xmin << "," << xmax
         << "] title=\"" << title << '"');

  auto* h = new TH1F(name.c_str(), title.c_str(), nbins, xmin, xmax);
  if (!h)
  {
    LOG(1, CLR_YELLOW, "  [getOrBookIsoHist] new TH1F failed for \"" << name << '"');
    if (prevDir) prevDir->cd();
    return nullptr;
  }
  H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}


// Record a histogram fill (for end-of-job diagnostics)
void RecoilJets::bumpHistFill(const std::string& trig, const std::string& hnameWithSuffix)
{
  if (trig.empty() || hnameWithSuffix.empty())
  {
    LOG(2, CLR_YELLOW, "  [bumpHistFill] empty trig or hist name (\"" << trig
           << "\", \"" << hnameWithSuffix << "\") – ignoring");
    return;
  }
  const std::string key = trig + "::" + hnameWithSuffix;
  auto it = m_histFill.find(key);
  if (it == m_histFill.end())
  {
    m_histFill.emplace(key, 1);
    if (Verbosity() >= 6)
      LOG(6, CLR_BLUE, "    [bumpHistFill] first fill → \"" << key << "\" = 1");
  }
  else
  {
    ++(it->second);
    if (Verbosity() >= 7)
      LOG(7, CLR_BLUE, "    [bumpHistFill] ++ \"" << key << "\" → " << it->second);
  }
}


TH1F* RecoilJets::getOrBookSSHist(const std::string& trig,
                                  const std::string& varKey,
                                  const std::string& tagKey,
                                  int etIdx, int centIdx)
{
  const std::string base   = "h_ss_" + varKey + "_" + tagKey;
  const std::string suffix = suffixForBins(etIdx, centIdx);
  const std::string name   = base + suffix;

  if (Verbosity() >= 5)
    LOG(5, CLR_BLUE, "  [getOrBookSSHist] trig=\"" << trig << "\" varKey=\"" << varKey
           << "\" tagKey=\"" << tagKey << "\" etIdx=" << etIdx << " centIdx=" << centIdx
           << " → name=\"" << name << "\"");

  if (trig.empty() || varKey.empty() || tagKey.empty())
  {
    LOG(2, CLR_YELLOW, "  [getOrBookSSHist] empty trig/varKey/tagKey – returning nullptr");
    return nullptr;
  }

  auto& H = qaHistogramsByTrigger[trig];

  if (auto it = H.find(name); it != H.end())
  {
    if (auto* h = dynamic_cast<TH1F*>(it->second))
    {
      if (Verbosity() >= 6)
        LOG(6, CLR_GREEN, "    [getOrBookSSHist] reusing existing TH1F \"" << name << "\"");
      return h;
    }
    LOG(2, CLR_YELLOW, "    [getOrBookSSHist] name clash: object \"" << name
                       << "\" exists but is not TH1F – replacing it");
    H.erase(it);
  }

  if (!out || !out->IsOpen())
  {
    LOG(1, CLR_YELLOW, "  [getOrBookSSHist] output TFile invalid/null – returning nullptr");
    return nullptr;
  }

  // Choose sane, variable-specific ranges (matches your previous defaults)
  int    nb = 120;
  double lo = 0.0, hi = 1.2;
  if (varKey == "weta" || varKey == "wphi") { nb = 120; lo = 0.0; hi = 1.2; }
  else if (varKey == "et1")                 { nb = 120; lo = 0.0; hi = 1.2; }
  else if (varKey == "e11e33" || varKey == "e32e35") { nb = 120; lo = 0.0; hi = 1.2; }
  else
  {
    LOG(3, CLR_YELLOW, "  [getOrBookSSHist] unknown varKey \"" << varKey
         << "\" – using default nb=" << nb << " range=[" << lo << "," << hi << "]");
  }

  TDirectory* const prevDir = gDirectory;
  TDirectory* dir = out->GetDirectory(trig.c_str());
  if (!dir) dir = out->mkdir(trig.c_str());
  if (!dir)
  {
    LOG(1, CLR_YELLOW, "  [getOrBookSSHist] failed to create/access directory \"" << trig << "\"");
    if (prevDir) prevDir->cd();
    return nullptr;
  }

  dir->cd();

  const std::string title = base + ";" + varKey + ";Entries";

  if (Verbosity() >= 5)
    LOG(5, CLR_BLUE, "    [getOrBookSSHist] booking TH1F name=\"" << name
         << "\" bins=" << nb << " range=[" << lo << "," << hi
         << "] title=\"" << title << '"');

  auto* h = new TH1F(name.c_str(), title.c_str(), nb, lo, hi);
  if (!h)
  {
    LOG(1, CLR_YELLOW, "  [getOrBookSSHist] new TH1F failed for \"" << name << '"');
    if (prevDir) prevDir->cd();
    return nullptr;
  }
  H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}




void RecoilJets::fillIsoSSTagCounters(const std::string& trig,
                                      const RawCluster* clus,
                                      const SSVars& v,
                                      double et_gamma,
                                      int centIdx,
                                      PHCompositeNode* topNode)
{
  if (!clus)
  {
    LOG(3, CLR_YELLOW, "  [fillIsoSSTagCounters] clus==nullptr – abort slice fill");
    return;
  }

  const int etIdx = findEtBin(et_gamma);
  if (etIdx < 0)
  {
    if (Verbosity() >= 4)
      LOG(4, CLR_YELLOW, "  [fillIsoSSTagCounters] ET bin not found for ET^γ=" << et_gamma
                          << " – skipping fills");
    return;
  }

  const int effCentIdx = (m_isAuAu ? centIdx : -1);
  const std::string slice = suffixForBins(etIdx, effCentIdx);

  // Isolation value and its histogram
  const double eiso_et = eiso(clus, topNode);
  if (auto* hIso = getOrBookIsoHist(trig, etIdx, effCentIdx))
  {
    hIso->Fill(eiso_et);
    bumpHistFill(trig, "h_Eiso" + slice);
  }
  else
  {
    LOG(2, CLR_YELLOW, "  [fillIsoSSTagCounters] getOrBookIsoHist returned nullptr for slice \"" << slice << '"');
  }

  // Selections
  const double thrIso = (m_isoA + m_isoB * et_gamma);
  const bool   iso    = (eiso_et < thrIso);
  const TightTag tag  = classifyPhotonTightness(v);

  if (Verbosity() >= 5)
  {
    LOG(5, CLR_BLUE, "    [fillIsoSSTagCounters] slice=\"" << slice
        << "\" eiso=" << eiso_et << " thrIso=" << thrIso
        << " isoPass=" << iso
        << " tightTag=" << static_cast<int>(tag));
  }

  // Update per-slice counters (in-memory bookkeeping)
  auto& S = m_catByTrig[trig][slice];
  S.seen += 1;
  (tag == TightTag::kTight ? S.tight : S.nonTight) += 1;
  (iso ? S.isoPass : S.isoFail) += 1;

  // Shower-shape histograms for tight/nontight only
  if (tag == TightTag::kTight || tag == TightTag::kNonTight)
  {
    const std::string tagKey = (tag == TightTag::kTight) ? "tight" : "nontight";
    auto fillSS = [&](const std::string& key, double val)
    {
      if (auto* h = getOrBookSSHist(trig, key, tagKey, etIdx, effCentIdx))
      {
        h->Fill(val);
        bumpHistFill(trig, "h_ss_" + key + "_" + tagKey + slice);
      }
      else
      {
        LOG(2, CLR_YELLOW, "    [fillSS] getOrBookSSHist returned nullptr for key=\"" << key
            << "\" tagKey=\"" << tagKey << "\" slice=\"" << slice << '"');
      }
    };

    fillSS("weta",   v.weta_cogx);
    fillSS("wphi",   v.wphi_cogx);
    fillSS("et1",    v.et1);
    fillSS("e11e33", v.e11_over_e33);
    fillSS("e32e35", v.e32_over_e35);
  }

  // Marginal counters (h_tight/nonTight & h_iso/nonIso)
  {
    const char* name = (tag == TightTag::kTight) ? "h_tightCount_ET" : "h_nonTightCount_ET";
    if (auto* h = getOrBookCountHist(trig, name, etIdx, effCentIdx))
    {
      h->Fill(1);
      bumpHistFill(trig, std::string(name) + slice);
    }
    else
    {
      LOG(2, CLR_YELLOW, "  [fillIsoSSTagCounters] getOrBookCountHist returned nullptr for \"" << name
                         << "\" slice=\"" << slice << '"');
    }
  }
  {
    const char* name = (iso ? "h_isolatedCount_ET" : "h_nonIsolatedCount_ET");
    if (auto* h = getOrBookCountHist(trig, name, etIdx, effCentIdx))
    {
      h->Fill(1);
      bumpHistFill(trig, std::string(name) + slice);
    }
    else
    {
      LOG(2, CLR_YELLOW, "  [fillIsoSSTagCounters] getOrBookCountHist returned nullptr for \"" << name
                         << "\" slice=\"" << slice << '"');
    }
  }

  // ID sideband bookkeeping (fail ≥ 2 tight cuts)
  if (tag == TightTag::kNonTight)
  {
    S.idSB_total += 1;

    if (auto* h = getOrBookCountHist(trig, "h_idSB_total", etIdx, effCentIdx))
    {
      h->Fill(1);
      bumpHistFill(trig, "h_idSB_total" + slice);
    }
    else
    {
      LOG(2, CLR_YELLOW, "  [fillIsoSSTagCounters] getOrBookCountHist nullptr for h_idSB_total slice \"" << slice << '"');
    }

    if (iso)
    {
      S.idSB_pass += 1;
      if (auto* h = getOrBookCountHist(trig, "h_idSB_pass", etIdx, effCentIdx))
      {
        h->Fill(1);
        bumpHistFill(trig, "h_idSB_pass" + slice);
      }
      else
      {
        LOG(2, CLR_YELLOW, "  [fillIsoSSTagCounters] getOrBookCountHist nullptr for h_idSB_pass slice \"" << slice << '"');
      }
    }
  }
}
