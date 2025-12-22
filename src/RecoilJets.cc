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
#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMap.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoDefs.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/RawCluster.h>
#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/coresoftware_local/offline/packages/CaloBase/PhotonClusterContainer.h"
#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/coresoftware_local/offline/packages/CaloBase/PhotonClusterv1.h"
#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/coresoftware_local/offline/packages/CaloReco/PhotonClusterBuilder.h"
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

#define CLR_BLUE   "\033[1;34m"
#define CLR_CYAN   "\033[1;36m"
#define CLR_GREEN  "\033[1;32m"
#define CLR_YELLOW "\033[1;33m"
#define CLR_MAGENTA "\033[1;35m"
#define CLR_RED    "\033[1;31m"
#define CLR_RESET  "\033[0m"

#undef  LOG
#define LOG(lvl, colour, msg)                                           \
  do {                                                                  \
    if (static_cast<int>(Verbosity()) >= static_cast<int>(lvl))         \
      std::cout << (colour) << msg << CLR_RESET << std::endl;           \
  } while (false)


/** Always print, independent of Verbosity() */
#define PROGRESS(MSG)                                                     \
    do {                                                                  \
        if (static_cast<int>(Verbosity()) >= 1)                           \
            std::cout << CLR_CYAN << MSG << CLR_RESET << std::endl;       \
    } while (false)


using namespace PhoIDCuts;

// Friendly label for printing tight category
static const char* tightTagName(RecoilJets::TightTag t)
{
  switch (t) {
    case RecoilJets::TightTag::kPreselectionFail: return "PreselectionFail";
    case RecoilJets::TightTag::kTight:            return "Tight";
    case RecoilJets::TightTag::kNonTight:         return "NonTight(>=2)";
    case RecoilJets::TightTag::kNeither:          return "Neither(1 fail)";
    default:                                      return "UNKNOWN";
  }
}

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
  /* ––– primary vertex ––––––––––––––––––––––––––––––––––––––––––––––––
     *
     * PhotonClusterBuilder uses MbdVertexMap(z) to compute cluster eta/phi
     * and the ET threshold. To keep RecoilJets fully consistent, we prefer
     * the SAME MBD z here.
     *
     * We still optionally read GlobalVertexMap so m_vtx/m_vx/m_vy remain
     * available for any legacy/diagnostic usage, but m_vz (used everywhere
     * in RecoilJets kinematics) will match the builder when MBD exists.
  * ------------------------------------------------------------------ */
  GlobalVertexMap* gvmap  = findNode::getClass<GlobalVertexMap>(top, "GlobalVertexMap");
  MbdVertexMap*    mbdmap = findNode::getClass<MbdVertexMap>(top, "MbdVertexMap");

  m_vtx = nullptr; m_vx = m_vy = 0.0; m_vz = 0.0;

  // (Optional) keep GlobalVertex pointer for x/y diagnostics
  if (gvmap && !gvmap->empty())
    {
      m_vtx = gvmap->begin()->second;
      if (m_vtx)
      {
        m_vx = m_vtx->get_x();
        m_vy = m_vtx->get_y();
      }
    }

    // Prefer MBD z (matches PhotonClusterBuilder)
    bool haveVz = false;
    if (mbdmap && !mbdmap->empty())
    {
      MbdVertex* mv = mbdmap->begin()->second;
      if (mv)
      {
        m_vz = mv->get_z();
        haveVz = std::isfinite(m_vz);
      }
  }

  // Fallback: use GlobalVertex z ONLY if MBD is unavailable
  if (!haveVz)
  {
      if (!gvmap)         { LOG(1, CLR_YELLOW, "  – MbdVertexMap missing AND GlobalVertexMap missing → skip event"); return false; }
      if (gvmap->empty()) { LOG(2, CLR_YELLOW, "  – MbdVertexMap missing AND GlobalVertexMap empty → skip event");   return false; }

      m_vtx = gvmap->begin()->second;
      if (!m_vtx)         { LOG(1, CLR_YELLOW, "  – MbdVertexMap missing AND Global vertex pointer null → skip event"); return false; }

      m_vx = m_vtx->get_x();
      m_vy = m_vtx->get_y();
      m_vz = m_vtx->get_z();
      haveVz = std::isfinite(m_vz);
  }

  if (!haveVz)
  {
      LOG(1, CLR_YELLOW, "  – vertex z is non-finite → skip event");
      return false;
  }

  LOG(4, CLR_BLUE, "    [fetchNodes] dataset=" << (m_isAuAu ? "Au+Au" : "p+p")
                                                << "  vz(used)=" << std::fixed << std::setprecision(2) << m_vz);

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

  if (!m_photons)
  {
      LOG(0, CLR_YELLOW,
          "    [fetchNodes] PHOTONCLUSTER_CEMC is MISSING → ABORTEVENT. "
          "PhotonClusterBuilder likely did not run or PHOTONCLUSTER_CEMC node name mismatch.");
      return false;
  }

  LOG(5, CLR_BLUE, "    [fetchNodes] PHOTONCLUSTER_CEMC present");

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
      os << "[InitRun] gamma-pT bins: {";
      for (std::size_t i = 0; i+1 < m_gammaPtBins.size(); ++i)
        os << (i? ", ":" ") << m_gammaPtBins[i] << "–" << m_gammaPtBins[i+1];
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
  // Common vertex histogram binning (match your |vz| cut of 60 cm)
  const int    nbVz  = 120;
  const double vzMin = -60.0;
  const double vzMax =  60.0;

  // ------------------------------------------------------------------
  // SIMULATION MODE: do NOT depend on trigger names.
  // Everything goes under a single directory: /SIM/
  // ------------------------------------------------------------------
  if (m_isSim)
  {
    const std::string trig = "SIM";

    TDirectory* dir = out->GetDirectory(trig.c_str());
    if (!dir) dir = out->mkdir(trig.c_str());
    dir->cd();

    HistMap& H = qaHistogramsByTrigger[trig];

    // 1) per-job event counter (SIM)
    const std::string hcnt = "cnt_" + trig;
    if (H.find(hcnt) == H.end())
    {
      auto* h = new TH1I(hcnt.c_str(), (hcnt + ";count;entries").c_str(), 1, 0.5, 1.5);
      h->GetXaxis()->SetBinLabel(1, "count");
      H[hcnt] = h;
    }

    // 2) vertex-z QA
    if (H.find("h_vertexZ") == H.end())
    {
      auto* hvz = new TH1F("h_vertexZ", "h_vertexZ;v_{z} [cm];Entries", nbVz, vzMin, vzMax);
      H["h_vertexZ"] = hvz;
    }

    // Optional: if you ever run a HI-like sim with centrality available
    if (m_isAuAu && H.find("h_centrality") == H.end())
    {
      auto* hc = new TH1F("h_centrality", "h_centrality;Centrality [%];Entries", 100, 0.0, 100.0);
      H["h_centrality"] = hc;
    }

    out->cd();
    return;
  }

  // ------------------------------------------------------------------
  // DATA MODE (existing behavior): per-trigger directories
  // ------------------------------------------------------------------
  if (m_isAuAu)
  {
    for (const auto& kv : triggerNameMapAuAu) // kv: std::pair<int,std::string>
    {
      const std::string trig = kv.second;

      // Make sure the trigger directory exists
      TDirectory* dir = out->GetDirectory(trig.c_str());
      if (!dir) dir = out->mkdir(trig.c_str());
      dir->cd();

      HistMap& H = qaHistogramsByTrigger[trig];

      // 1) per-trigger event counter
      const std::string hcnt = "cnt_" + trig;
      if (H.find(hcnt) == H.end())
      {
        auto* h = new TH1I(hcnt.c_str(), (hcnt + ";count;entries").c_str(), 1, 0.5, 1.5);
        h->GetXaxis()->SetBinLabel(1, "count");
        H[hcnt] = h;
      }

      // 2) vertex-z QA
      if (H.find("h_vertexZ") == H.end())
      {
        auto* hvz = new TH1F("h_vertexZ", "h_vertexZ;v_{z} [cm];Entries", nbVz, vzMin, vzMax);
        H["h_vertexZ"] = hvz;
      }

      // 3) centrality QA (Au+Au only)
      if (H.find("h_centrality") == H.end())
      {
        auto* hc = new TH1F("h_centrality", "h_centrality;Centrality [%];Entries", 100, 0.0, 100.0);
        H["h_centrality"] = hc;
      }

      out->cd();
    }
  }
  else
  {
    for (const auto& kv : triggerNameMap_pp) // kv: std::pair<std::string,std::string>
    {
      const std::string trig = kv.second;

      // Make sure the trigger directory exists
      TDirectory* dir = out->GetDirectory(trig.c_str());
      if (!dir) dir = out->mkdir(trig.c_str());
      dir->cd();

      HistMap& H = qaHistogramsByTrigger[trig];

      // 1) per-trigger event counter
      const std::string hcnt = "cnt_" + trig;
      if (H.find(hcnt) == H.end())
      {
        auto* h = new TH1I(hcnt.c_str(), (hcnt + ";count;entries").c_str(), 1, 0.5, 1.5);
        h->GetXaxis()->SetBinLabel(1, "count");
        H[hcnt] = h;
      }

      // 2) vertex-z QA
      if (H.find("h_vertexZ") == H.end())
      {
        auto* hvz = new TH1F("h_vertexZ", "h_vertexZ;v_{z} [cm];Entries", nbVz, vzMin, vzMax);
        H["h_vertexZ"] = hvz;
      }

      out->cd();
    }
  }
}



bool RecoilJets::firstEventCuts(PHCompositeNode* topNode,
                                std::vector<std::string>& activeTrig)
{
  m_lastReject = EventReject::None;
  activeTrig.clear();

  // ------------------------------------------------------------------
  // SIMULATION MODE:
  //   - skip ALL trigger logic
  //   - fill everything under a single key: "SIM"
  //   - still apply the global |vz| cut if enabled
  // ------------------------------------------------------------------
  if (m_isSim)
  {
    activeTrig.emplace_back("SIM");

    if (m_useVzCut && std::fabs(m_vz) >= m_vzCut)
    {
      m_lastReject = EventReject::Vz;
      return false;
    }

    m_lastReject = EventReject::None;
    return true;
  }

  if (!m_isAuAu)
  {
    // ---------- p+p path: accept if ANY configured pp trigger fired ----------
    if (!trigAna) { m_lastReject = EventReject::Trigger; return false; }

    trigAna->decodeTriggers(topNode);

    for (const auto& kv : triggerNameMap_pp)
    {
      const std::string& dbName   = kv.first;   // GL1/DB name
      const std::string& shortKey = kv.second;  // friendly key
      if (trigAna->didTriggerFire(dbName))
        activeTrig.push_back(shortKey);
    }

    // Gate: at least one of the configured pp triggers must have fired
    if (activeTrig.empty()) { m_lastReject = EventReject::Trigger; return false; }
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

    if (activeTrig.empty()) { m_lastReject = EventReject::Trigger; return false; }
  }

  // ---------- Global |vz| veto (applies to both datasets) ----------
  if (m_useVzCut && std::fabs(m_vz) >= m_vzCut)
  {
    m_lastReject = EventReject::Vz;
    return false;
  }

  // Uniform downstream behavior
  if (activeTrig.empty()) activeTrig.emplace_back("ALL");

  m_lastReject = EventReject::None;
  return true;
}




int RecoilJets::process_event(PHCompositeNode* topNode)
{
  /* ------------------------------------------------------------------ */
  /* 0) Banner & counter                                                */
  /* ------------------------------------------------------------------ */
  ++event_count;
  ++m_bk.evt_seen;
  std::cout << "==================== processing event "
            << std::setw(6) << event_count
            << " ====================" << std::endl;

  /* ------------------------------------------------------------------ */
  /* 1) Mandatory nodes                                                 */
  /* ------------------------------------------------------------------ */
  LOG(4, CLR_BLUE, "  [process_event] – node sanity");
  if (!fetchNodes(topNode))
  {
    LOG(4, CLR_YELLOW, "    mandatory node(s) missing → ABORTEVENT");
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  /* ------------------------------------------------------------------ */
  /* 2) Trigger gating (pp & Au+Au) — unified in firstEventCuts()       */
  /* ------------------------------------------------------------------ */
  std::vector<std::string> activeTrig;
  if (!firstEventCuts(topNode, activeTrig))
  {
    ++m_evtNoTrig;  // keep your legacy counter
    if (m_lastReject == EventReject::Trigger) ++m_bk.evt_fail_trigger;
    else if (m_lastReject == EventReject::Vz) ++m_bk.evt_fail_vz;

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

  /* ------------------------------------------------------------------ */
  /* 6) Pure jet QA (independent of photon pipeline)                     */
  /*     Filled once per accepted event, after centrality/vz.            */
  /* ------------------------------------------------------------------ */
  const int centIdxForJets = (m_isAuAu ? findCentBin(m_centBin) : -1);
  for (const auto& jnm : kJetRadii)
  {
      fillInclusiveJetQA(activeTrig, centIdxForJets, jnm.key);
  }

  processCandidates(topNode, activeTrig);

  ++m_bk.evt_accepted;
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
            std::cout << "slice                                 | iso∧tight | ¬iso∧tight | iso∧¬tight | ¬iso∧¬tight |   seen | SBtot | SBpass\n";
            std::cout << "--------------------------------------+-----------+------------+------------+-------------+--------+-------+-------\n";
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
                        << std::right << std::setw(9)  << S.n_iso_tight        << " | "
                        << std::setw(10) << S.n_nonIso_tight     << " | "
                        << std::setw(10) << S.n_iso_nonTight     << " | "
                        << std::setw(11) << S.n_nonIso_nonTight  << " | "
                        << std::setw(6)  << S.seen               << " | "
                        << std::setw(5)  << S.idSB_total         << " | "
                        << std::setw(5)  << S.idSB_pass          << "\n";
            }

        }

        // NEW: job-wide meticulous cutflow
        printCutSummary();
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
RecoilJets::SSVars RecoilJets::makeSSFromPhoton(const PhotonClusterv1* pho, double pt_gamma) const
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
  v.pt_gamma      = pt_gamma;
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

  // ========================= PHOTON path (preferred) ========================
  if (m_photons)
  {
    PhotonClusterContainer::ConstRange prange = m_photons->getClusters();

    // Quick emptiness check
    std::size_t nPho = 0;
    for (auto tmp = prange.first; tmp != prange.second; ++tmp) ++nPho;
    LOG(4, CLR_BLUE, "    [processCandidates] photons present: " << nPho);
    if (nPho == 0)
      LOG(3, CLR_YELLOW, "    [processCandidates] photon container is EMPTY for this event");

    // Per-event counters
    std::size_t nUsed = 0;
    std::size_t nSkipEta = 0, nSkipEtBin = 0, nSkipPhi = 0, nNoRC = 0, nNotIso = 0, nNotTight = 0;

    try
    {
        // ---- Event-level de-duplication for jet matching / JES fills ----
        // If multiple photons pass (iso ∧ tight) in the same event, keep ONLY
        // the leading one in pT^gamma and do jet matching + JES3 fills once.
        bool   haveLeadIsoTight = false;
        int    leadPhoIndex     = -1;
        int    leadPtIdx        = -1;
        double leadPtGamma      = -1.0;
        double leadPhiGamma     = 0.0;

        int iPho = 0;
        for (auto pit = prange.first; pit != prange.second; ++pit, ++iPho)
        {
        // Concrete types
        const auto* pho = dynamic_cast<const PhotonClusterv1*>(pit->second);
        if (!pho)
        {
          LOG(5, CLR_YELLOW, "      [pho#" << iPho << "] pointer is not PhotonClusterv1 – skipping");
          continue;
        }

        // Underlying RawCluster is needed to build kinematics at the event vertex
        const RawCluster* rc = dynamic_cast<const RawCluster*>(pit->second);
        if (!rc)
        {
          ++nNoRC; ++m_bk.pho_noRC;
          if (Verbosity() >= 5)
            LOG(5, CLR_YELLOW, "      [pho#" << iPho << "] RawCluster pointer is NULL – cannot compute kinematics");
          continue;
        }

        ++m_bk.pho_total;

        // --------------------------------------------------------------
        // Use the EXACT kinematics produced by PhotonClusterBuilder:
        //   - vertex_z    (MBD z used in builder)
        //   - cluster_eta / cluster_phi
        //   - cluster_pt
        //
        // This prevents ET/pT mismatches from recomputing with a different
        // vertex choice or a different construction path.
        // --------------------------------------------------------------
        double eta      = pho->get_shower_shape_parameter("cluster_eta");
        double phi      = pho->get_shower_shape_parameter("cluster_phi");
        double pt_gamma = pho->get_shower_shape_parameter("cluster_pt");
        const double energy = rc->get_energy();

        // Strict mode: ONLY use PhotonClusterBuilder-provided kinematics.
        // If they are missing/non-finite, skip this candidate and print why.
        if (!std::isfinite(eta) || !std::isfinite(phi) || !std::isfinite(pt_gamma))
        {
            LOG(0, CLR_YELLOW,
                "      [pho#" << iPho << "] PhotonClusterBuilder kinematics are non-finite/missing: "
                << "eta=" << eta << " phi=" << phi << " pt=" << pt_gamma
                << " → skipping candidate (no RawClusterUtility fallback allowed).");
            continue;
        }

        // -------- early transverse-scale gate (fast reject) ---------------
        constexpr double kMinPtGamma = 2.0; // GeV (kept as a fast sanity floor)

        // quick reject: E < 2 GeV implies pT < 2 GeV for any η
        if (energy < kMinPtGamma)
        {
          ++nSkipEtBin; ++m_bk.pho_early_E;
          if (Verbosity() >= 6)
            LOG(6, CLR_YELLOW, "      [pho#" << iPho << "] early reject (E<2 GeV): E=" << energy);
          continue;
        }

        // |η| fiducial cut
        if (!std::isfinite(eta) || std::fabs(eta) >= m_etaAbsMax)
        {
          ++nSkipEta; ++m_bk.pho_eta_fail;
          if (Verbosity() >= 2)
          {
            std::cout << CLR_MAGENTA
                      << "      [skip:eta pho#" << iPho << "]"
                      << "  pT=" << std::fixed << std::setprecision(2) << pt_gamma
                      << "  E="  << std::fixed << std::setprecision(2) << energy
                      << "  eta="<< std::fixed << std::setprecision(3) << eta
                      << "  phi="<< std::fixed << std::setprecision(3) << phi
                      << "  (|eta|>=" << m_etaAbsMax << ")"
                      << CLR_RESET << std::endl;
          }
          continue;
        }

        // Enforce a minimal transverse scale sanity floor
        if (!std::isfinite(pt_gamma) || pt_gamma < kMinPtGamma)
        {
          ++nSkipEtBin; ++m_bk.pho_early_E;
          if (Verbosity() >= 6)
            LOG(6, CLR_YELLOW, "      [pho#" << iPho << "] early reject (pT<2 GeV): pT=" << pt_gamma);
          continue;
        }

        // Now do your configured pT binning (drives ALL pT-sliced histograms + JES objects)
        const int ptIdx = findPtBin(pt_gamma);
        if (ptIdx < 0)
        {
          ++nSkipEtBin; ++m_bk.pho_etbin_out;
          if (Verbosity() >= 2)
          {
            std::cout << CLR_MAGENTA
                      << "      [skip:pTbin pho#" << iPho << "]"
                      << "  pT=" << std::fixed << std::setprecision(2) << pt_gamma
                      << "  E="  << std::fixed << std::setprecision(2) << energy
                      << "  eta="<< std::fixed << std::setprecision(3) << eta
                      << "  phi="<< std::fixed << std::setprecision(3) << phi
                      << CLR_RESET << std::endl;
          }
          continue;
        }

        // φ to use for jet matching (used later for jet scan)
        const double phi_gamma = phi;

        // ---- PURE isolation QA (per pT/[cent] slice), BEFORE any pre-selection ----
        // This block is intentionally independent of:
        //   - preselection
        //   - tightness
        //   - any later xJ usability decision
        {
            const int effCentIdx = (m_isAuAu ? centIdx : -1);
            const std::string slice = suffixForBins(ptIdx, effCentIdx);

            // Total isolation (existing behavior)
            const double eiso_tot = eiso(rc, topNode);

            // Component isolation (PhotonClusterBuilder iso_* pieces)
            // Default to fail-safe (goes to overflow with your [-5,12] binning).
            double eiso_emcal  = 1e9;
            double eiso_hcalin = 1e9;
            double eiso_hcalout= 1e9;

            const int cone10 = static_cast<int>(std::lround(10.0 * m_isoConeR));
            const char* k_em = nullptr;
            const char* k_hi = nullptr;
            const char* k_ho = nullptr;

            if (cone10 == 3)
            {
              k_em = "iso_03_emcal";
              k_hi = "iso_03_hcalin";
              k_ho = "iso_03_hcalout";
            }
            else if (cone10 == 4)
            {
              k_em = "iso_04_emcal";
              k_hi = "iso_04_hcalin";
              k_ho = "iso_04_hcalout";
            }

            if (pho && k_em && k_hi && k_ho)
            {
              eiso_emcal   = pho->get_shower_shape_parameter(k_em);
              eiso_hcalin  = pho->get_shower_shape_parameter(k_hi);
              eiso_hcalout = pho->get_shower_shape_parameter(k_ho);

              // If any component is non-finite, treat as fail-safe (overflow)
              if (!std::isfinite(eiso_emcal) || !std::isfinite(eiso_hcalin) || !std::isfinite(eiso_hcalout))
              {
                eiso_emcal = eiso_hcalin = eiso_hcalout = 1e9;
              }
            }

            // Pure isolation threshold decision (signal line only)
            const double thrIso  = m_isoA + m_isoB * pt_gamma;
            const bool   isoPass = (eiso_tot < thrIso);

            for (const auto& trigShort : activeTrig)
            {
              // (A) Existing total isolation histogram (unchanged name/slicing)
              if (auto* hIso = getOrBookIsoHist(trigShort, ptIdx, effCentIdx))
              {
                hIso->Fill(eiso_tot);
                bumpHistFill(trigShort, std::string("h_Eiso") + slice);
              }

              // (B) NEW component isolation histograms (EMCal / IHCal / OHCal)
              if (auto* hEm = getOrBookIsoPartHist(trigShort, "h_Eiso_emcal",
                                                   "E_{T}^{iso,EMCal} [GeV]",
                                                   ptIdx, effCentIdx))
              {
                hEm->Fill(eiso_emcal);
                bumpHistFill(trigShort, std::string("h_Eiso_emcal") + slice);
              }

              if (auto* hHi = getOrBookIsoPartHist(trigShort, "h_Eiso_hcalin",
                                                   "E_{T}^{iso,IHCAL} [GeV]",
                                                   ptIdx, effCentIdx))
              {
                hHi->Fill(eiso_hcalin);
                bumpHistFill(trigShort, std::string("h_Eiso_hcalin") + slice);
              }

              if (auto* hHo = getOrBookIsoPartHist(trigShort, "h_Eiso_hcalout",
                                                   "E_{T}^{iso,OHCAL} [GeV]",
                                                   ptIdx, effCentIdx))
              {
                hHo->Fill(eiso_hcalout);
                bumpHistFill(trigShort, std::string("h_Eiso_hcalout") + slice);
              }

              // (C) NEW isolation decision counts (PASS vs FAIL), independent of SS
              if (auto* hDec = getOrBookIsoDecisionHist(trigShort, ptIdx, effCentIdx))
              {
                hDec->Fill(isoPass ? 1 : 2);   // bin 1 = PASS, bin 2 = FAIL
                bumpHistFill(trigShort, std::string("h_isoDecision") + slice);
              }
            }
        }

        // 1) Build shower-shape inputs (for preselection and tightness)
        const SSVars v = makeSSFromPhoton(pho, pt_gamma);
        ++m_bk.pho_reached_pre_iso;

        // ---------- Preselection breakdown (count by criterion) ----------
        const bool pass_e11e33 = (v.e11_over_e33 < PRE_E11E33_MAX);
        const bool pass_et1    = in_open_interval(v.et1, PRE_ET1_MIN, PRE_ET1_MAX);
        const bool pass_e32e35 = in_open_interval(v.e32_over_e35, PRE_E32E35_MIN, PRE_E32E35_MAX);
        const bool pass_weta   = (v.weta_cogx < PRE_WETA_MAX);

        bool pre_ok = pass_e11e33 && pass_et1 && pass_e32e35 && pass_weta;
        if (!pre_ok)
        {
          // 3) Book/fill per-cut, per-bin FAILURE histograms
          const int effCentIdx = (m_isAuAu ? centIdx : -1);
          const std::string slice = suffixForBins(ptIdx, effCentIdx);

          for (const auto& trigShort : activeTrig)
          {
            if (!pass_weta)
            {
              if (auto* h = getOrBookCountHist(trigShort, "h_preFail_weta", ptIdx, effCentIdx))
              { h->Fill(1); bumpHistFill(trigShort, std::string("h_preFail_weta") + slice); }
            }
            if (!pass_et1)
            {
              const char* base = (v.et1 <= PRE_ET1_MIN) ? "h_preFail_et1_low" : "h_preFail_et1_high";
              if (auto* h = getOrBookCountHist(trigShort, base, ptIdx, effCentIdx))
              { h->Fill(1); bumpHistFill(trigShort, std::string(base) + slice); }
            }
            if (!pass_e11e33)
            {
              if (auto* h = getOrBookCountHist(trigShort, "h_preFail_e11e33_high", ptIdx, effCentIdx))
              { h->Fill(1); bumpHistFill(trigShort, std::string("h_preFail_e11e33_high") + slice); }
            }
            if (!pass_e32e35)
            {
              const char* base = (v.e32_over_e35 <= PRE_E32E35_MIN) ? "h_preFail_e32e35_low" : "h_preFail_e32e35_high";
              if (auto* h = getOrBookCountHist(trigShort, base, ptIdx, effCentIdx))
              { h->Fill(1); bumpHistFill(trigShort, std::string(base) + slice); }
            }
          }

          // Maintain existing bookkeeping
          if (!pass_weta)   ++m_bk.pre_fail_weta;
          if (!pass_et1)    { if (v.et1 <= PRE_ET1_MIN) ++m_bk.pre_fail_et1_low; else ++m_bk.pre_fail_et1_high; }
          if (!pass_e11e33) ++m_bk.pre_fail_e11e33_high;
          if (!pass_e32e35) { if (v.e32_over_e35 <= PRE_E32E35_MIN) ++m_bk.pre_fail_e32e35_low; else ++m_bk.pre_fail_e32e35_high; }

          // 4) Human-readable, single-line diagnostic with ALL values and reasons
          if (Verbosity() >= 4)
          {
            const char* et1_state    = pass_et1    ? "PASS" : (v.et1 <= PRE_ET1_MIN ? "FAIL(low)" : "FAIL(high)");
            const char* e32e35_state = pass_e32e35 ? "PASS" : (v.e32_over_e35 <= PRE_E32E35_MIN ? "FAIL(low)" : "FAIL(high)");

            std::ostringstream msg;
            msg << "      [pho#" << iPho << "] preselection FAIL"
                << " | weta="    << std::fixed << std::setprecision(3) << v.weta_cogx
                << "  cut:<"     << PRE_WETA_MAX       << " → " << (pass_weta ? "PASS" : "FAIL")
                << " | et1="     << std::fixed << std::setprecision(3) << v.et1
                << "  cut:("     << PRE_ET1_MIN        << "," << PRE_ET1_MAX       << ") → " << et1_state
                << " | e11/e33=" << std::fixed << std::setprecision(3) << v.e11_over_e33
                << "  cut:<"     << PRE_E11E33_MAX     << " → " << (pass_e11e33 ? "PASS" : "FAIL")
                << " | e32/e35=" << std::fixed << std::setprecision(3) << v.e32_over_e35
                << "  cut:("     << PRE_E32E35_MIN     << "," << PRE_E32E35_MAX    << ") → " << e32e35_state
                << " | pT^{#gamma}=" << std::fixed << std::setprecision(2) << v.pt_gamma;
            LOG(4, CLR_MAGENTA, msg.str());
          }

          ++nNotTight; // keep legacy counter of “not usable” for xJ
          continue;
        }

        if (Verbosity() >= 4)
        {
          std::ostringstream msg;
          msg << "      [pho#" << iPho << "] preselection PASS"
              << " | weta="    << std::fixed << std::setprecision(3) << v.weta_cogx  << "  cut:<"  << PRE_WETA_MAX    << " → PASS"
              << " | et1="     << std::fixed << std::setprecision(3) << v.et1        << "  cut:("  << PRE_ET1_MIN     << "," << PRE_ET1_MAX     << ") → PASS"
              << " | e11/e33=" << std::fixed << std::setprecision(3) << v.e11_over_e33 << "  cut:<" << PRE_E11E33_MAX  << " → PASS"
              << " | e32/e35=" << std::fixed << std::setprecision(3) << v.e32_over_e35 << "  cut:(" << PRE_E32E35_MIN  << "," << PRE_E32E35_MAX  << ") → PASS"
              << " | pT^{#gamma}=" << std::fixed << std::setprecision(2) << v.pt_gamma;
          LOG(4, CLR_RED, msg.str());
        }
        ++m_bk.pre_pass;

        // ---------- Isolation (count pass/fail) ----------
        const bool iso = isIsolated(rc, pt_gamma, topNode);
        if (iso) ++m_bk.iso_pass; else ++m_bk.iso_fail;

        // ---------- Tight classification breakdown ----------
        const double w_hi = tight_w_hi(v.pt_gamma);
        const bool t_weta   = in_open_interval(v.weta_cogx,  TIGHT_W_LO, w_hi);
        const bool t_wphi   = in_open_interval(v.wphi_cogx,  TIGHT_W_LO, w_hi);
        const bool t_e11e33 = in_open_interval(v.e11_over_e33, TIGHT_E11E33_MIN, TIGHT_E11E33_MAX);
        const bool t_et1    = in_open_interval(v.et1,          TIGHT_ET1_MIN,    TIGHT_ET1_MAX);
        const bool t_e32e35 = in_open_interval(v.e32_over_e35, TIGHT_E32E35_MIN, TIGHT_E32E35_MAX);

        int tight_fails = 0;
        if (!t_weta)   { ++m_bk.tight_fail_weta;   ++tight_fails; }
        if (!t_wphi)   { ++m_bk.tight_fail_wphi;   ++tight_fails; }
        if (!t_e11e33) { ++m_bk.tight_fail_e11e33; ++tight_fails; }
        if (!t_et1)    { ++m_bk.tight_fail_et1;    ++tight_fails; }
        if (!t_e32e35) { ++m_bk.tight_fail_e32e35; ++tight_fails; }

        RecoilJets::TightTag tightTag;
        if (tight_fails == 0)      { tightTag = TightTag::kTight;      ++m_bk.tight_tight; }
        else if (tight_fails >= 2) { tightTag = TightTag::kNonTight;   ++m_bk.tight_nonTight; }
        else                       { tightTag = TightTag::kNeither;    ++m_bk.tight_neither; }

        // NEW: record the 2×2 (iso, tight) category + SS variable hists
        // This also fills h_Eiso once and prints a detailed decision line.
        for (const auto& trigShort : activeTrig)
        {
          fillIsoSSTagCounters(trigShort, rc, v, pt_gamma, centIdx, topNode);
        }

            // Your original xJ usable gate (unchanged)
            if (!(iso && tightTag == TightTag::kTight))
            {
              if (!iso) ++nNotIso;
              if (tightTag != TightTag::kTight) ++nNotTight;
              if (Verbosity() >= 6)
                LOG(6, CLR_BLUE, "      [pho#" << iPho << "] NOT used for xJ (iso=" << iso
                                               << ", tightTag=" << tightTagName(tightTag) << ")");
              continue;
            }

            // ---- IMPORTANT CHANGE ----
            // Do NOT jet-match here. If >1 photon passes iso∧tight in the same event,
            // jet-matching here would double-fill xJ/JES histograms.
            //
            // Instead: keep ONLY the event-leading iso∧tight photon (highest pT^gamma),
            // and do jet matching + JES3 fills ONCE after the photon loop.
            if (!haveLeadIsoTight || pt_gamma > leadPtGamma)
            {
              haveLeadIsoTight = true;
              leadPhoIndex     = iPho;
              leadPtIdx        = ptIdx;
              leadPtGamma      = pt_gamma;
              leadPhiGamma     = phi_gamma;

              if (Verbosity() >= 6)
                LOG(6, CLR_GREEN, "      [pho#" << iPho << "] marked as event-leading iso∧tight photon for jet matching (pT="
                                                << std::fixed << std::setprecision(2) << leadPtGamma << ")");
            }

            // Defer jet matching (prevents double-filling JES TH3s if multiple photons pass).
            continue;
        } // photon loop

        // ------------------------------------------------------------------
        // Jet matching + JES fills ONCE per event using the event-leading
        // iso∧tight photon (by pT^gamma).
        // ------------------------------------------------------------------
        if (haveLeadIsoTight)
        {
          // r02 jet container
            JetContainer* jets = nullptr;
            if (auto itJ = m_jets.find("r02"); itJ != m_jets.end()) jets = itJ->second;

            // Centrality index used for matching-QA (centrality-only suffix in Au+Au)
            const int effCentIdx_M = (m_isAuAu ? centIdx : -1);

            if (!jets)
            {
              LOG(3, CLR_YELLOW, "      [lead pho#" << leadPhoIndex << "] jet container r02 is nullptr – cannot form xJ");

              // Still fill matching-QA once per event-leading iso∧tight photon
              // (treat as "NoJetPt" and maxDphi underflow sentinel).
              const int    status = 1;      // NoJetPt
              const double maxDphiFill = -0.01; // underflow sentinel (no jets pT+eta)

              for (const auto& trigShort : activeTrig)
              {
                if (auto* hs = getOrBookMatchStatusVsPtGamma(trigShort, effCentIdx_M))
                { hs->Fill(leadPtGamma, status); bumpHistFill(trigShort, hs->GetName()); }

                if (auto* hm = getOrBookMatchMaxDphiVsPtGamma(trigShort, effCentIdx_M))
                { hm->Fill(leadPtGamma, maxDphiFill); bumpHistFill(trigShort, hm->GetName()); }

                if (auto* pn = getOrBookNRecoilJetsVsPtGamma(trigShort, effCentIdx_M))
                { pn->Fill(leadPtGamma, 0.0); bumpHistFill(trigShort, pn->GetName()); }
              }
            }
            else
            {
                // Scan jets (R=0.2):
                //   - baseline pT gate:        pT > m_minJetPt
                //   - fiducial containment:    |eta_jet| < m_jetEtaAbsMax
                //
                // Apply the SAME fiducial cut to:
                //   * recoil jet1 (away-side)
                //   * the jet2 candidate used for alpha
                // so alpha/xJ are not polluted by edge / partial-containment jets.
                int nPassPt   = 0;
                int nPassEta  = 0;
                int nPassDphi = 0;

                double bestPt = -1.0;
                const Jet* bestJet = nullptr;

                double all1Pt = -1.0;
                const Jet* all1Jet = nullptr;
                double all2Pt = -1.0;
                const Jet* all2Jet = nullptr;

                // NEW: max |Δφ| over jets that pass pT+eta (even if they fail the Δφ cut)
                double maxDphi = -1.0;

                for (const Jet* j : *jets)
                {
                  if (!j) continue;

                  const double jpt  = j->get_pt();
                  const double jeta = j->get_eta();
                  const double jphi = j->get_phi();

                  if (!std::isfinite(jpt) || !std::isfinite(jeta) || !std::isfinite(jphi)) continue;
                  if (jpt < m_minJetPt) continue;
                  ++nPassPt;

                  // Fiducial acceptance (jet-axis containment)
                  if (std::fabs(jeta) >= m_jetEtaAbsMax) continue;
                  ++nPassEta;

                  // Track top-2 jets in the event by pT (for "jet2" definition)
                  if (jpt > all1Pt)
                  {
                    all2Pt = all1Pt; all2Jet = all1Jet;
                    all1Pt = jpt;    all1Jet = j;
                  }
                  else if (jpt > all2Pt)
                  {
                    all2Pt = jpt; all2Jet = j;
                  }

                  // Δφ for diagnostics and recoil selection (away-side of LEADING photon)
                  const double dphiAbs = std::fabs(TVector2::Phi_mpi_pi(jphi - leadPhiGamma));
                  if (std::isfinite(dphiAbs) && dphiAbs > maxDphi) maxDphi = dphiAbs;

                  if (dphiAbs >= m_minBackToBack)
                  {
                    ++nPassDphi;

                    if (jpt > bestPt)
                    {
                      bestPt  = jpt;
                      bestJet = j;
                    }
                  }
                }

                // NEW: matching status category (filled once per event-leading iso∧tight photon)
                int status = 0;
                if (nPassPt == 0)                 status = 1; // NoJetPt
                else if (nPassEta == 0)           status = 2; // NoJetEta
                else if (nPassDphi == 0 || !bestJet) status = 3; // NoBackToBack
                else                               status = 4; // Matched

                // NEW: maxDphi fill value (underflow sentinel if no jets pass pT+eta)
                const double maxDphiFill = (maxDphi >= 0.0 ? maxDphi : -0.01);

                // NEW: fill matching-QA histograms (centrality-only suffix; pTgamma is axis)
                for (const auto& trigShort : activeTrig)
                {
                  if (auto* hs = getOrBookMatchStatusVsPtGamma(trigShort, effCentIdx_M))
                  { hs->Fill(leadPtGamma, status); bumpHistFill(trigShort, hs->GetName()); }

                  if (auto* hm = getOrBookMatchMaxDphiVsPtGamma(trigShort, effCentIdx_M))
                  { hm->Fill(leadPtGamma, maxDphiFill); bumpHistFill(trigShort, hm->GetName()); }

                  if (auto* pn = getOrBookNRecoilJetsVsPtGamma(trigShort, effCentIdx_M))
                  { pn->Fill(leadPtGamma, static_cast<double>(nPassDphi)); bumpHistFill(trigShort, pn->GetName()); }

                  // Fill selected recoil Δφ only when matched
                  if (bestJet)
                  {
                    const double dphiSel = std::fabs(TVector2::Phi_mpi_pi(bestJet->get_phi() - leadPhiGamma));
                    if (auto* hd = getOrBookMatchDphiVsPtGamma(trigShort, effCentIdx_M))
                    { hd->Fill(leadPtGamma, dphiSel); bumpHistFill(trigShort, hd->GetName()); }

                    // Optional: is selected recoil the leading jet among pT+eta jets?
                    const double isLeading = (all1Jet && bestJet == all1Jet) ? 1.0 : 0.0;
                    if (auto* hl = getOrBookRecoilIsLeadingVsPtGamma(trigShort, effCentIdx_M))
                    { hl->Fill(leadPtGamma, isLeading); bumpHistFill(trigShort, hl->GetName()); }
                  }
                }

                // Keep your physics-output fills unchanged (xJ, jet1/2 pT, alpha, JES3)
                if (bestPt > 0.0 && bestJet)
                {
                  // xJ uses photon pT
                  const double xJ = bestPt / leadPtGamma;

                  // Define jet2 as the highest-pT jet excluding the recoil jet1.
                  // If none exists, jet2Pt=0 and jet2Jet=nullptr, so alpha=0.
                  const Jet* jet2Jet = nullptr;
                  double jet2Pt = 0.0;
                  if (all1Jet && all1Jet != bestJet)      { jet2Jet = all1Jet; jet2Pt = all1Pt; }
                  else if (all2Jet && all2Jet != bestJet) { jet2Jet = all2Jet; jet2Pt = all2Pt; }
                  jet2Pt = std::max(0.0, jet2Pt);

                  const double alpha = (leadPtGamma > 0.0 ? (jet2Pt / leadPtGamma) : 0.0);

                  if (Verbosity() >= 5)
                    LOG(5, CLR_GREEN, "      [lead pho#" << leadPhoIndex << "] xJ = " << std::fixed << std::setprecision(3)
                                                    << xJ
                                                    << "  jet1Pt=" << bestPt
                                                    << "  jet2Pt=" << jet2Pt
                                                    << "  alpha=" << alpha
                                                    << "  (pT^{#gamma}=" << leadPtGamma
                                                    << ", passed dphi=" << nPassDphi
                                                    << ", passed pt=" << nPassPt << ")");

                  const int effCentIdx = (m_isAuAu ? centIdx : -1);
                  const std::string slice    = suffixForBins(leadPtIdx, effCentIdx);
                  const std::string centOnly = suffixForBins(-1,        effCentIdx);

                  // Jet-only QA for the recoil selection (jet1) and alpha jet (jet2)
                  // Uses the SAME pT^gamma/centrality slicing as the xJ/JES histograms.
                  fillSelectedJetQA(activeTrig, leadPtIdx, effCentIdx, "r02", bestJet, jet2Jet);

                  for (const auto& trigShort : activeTrig)
                  {
                    // 1) xJ spectrum (per pT bin)
                    if (auto* hx = getOrBookXJHist(trigShort, leadPtIdx, effCentIdx))
                    {
                      hx->Fill(xJ);
                      bumpHistFill(trigShort, std::string("h_xJ") + slice);
                    }

                    // 2) jet1 pT spectrum (per pT bin)
                    if (auto* hj1 = getOrBookJet1PtHist(trigShort, leadPtIdx, effCentIdx))
                    {
                      hj1->Fill(bestPt);
                      bumpHistFill(trigShort, std::string("h_jet1Pt") + slice);
                    }

                    // 3) jet2 pT spectrum (per pT bin)
                    if (auto* hj2 = getOrBookJet2PtHist(trigShort, leadPtIdx, effCentIdx))
                    {
                      hj2->Fill(jet2Pt);
                      bumpHistFill(trigShort, std::string("h_jet2Pt") + slice);
                    }

                    // 4) alpha spectrum (per pT bin)
                    if (auto* ha = getOrBookAlphaHist(trigShort, leadPtIdx, effCentIdx))
                    {
                      ha->Fill(alpha);
                      bumpHistFill(trigShort, std::string("h_alpha") + slice);
                    }

                    // 5) TH3 correlations (cent-only suffix; pT^gamma is an axis)
                    if (auto* h3x = getOrBookJES3_xJ_alphaHist(trigShort, effCentIdx))
                    {
                      h3x->Fill(leadPtGamma, xJ, alpha);
                      bumpHistFill(trigShort, std::string("h_JES3_pT_xJ_alpha") + centOnly);
                    }

                    if (auto* h3j = getOrBookJES3_jet1Pt_alphaHist(trigShort, effCentIdx))
                    {
                      h3j->Fill(leadPtGamma, bestPt, alpha);
                      bumpHistFill(trigShort, std::string("h_JES3_pT_jet1Pt_alpha") + centOnly);
                    }
                  }

                  ++nUsed;
                }
                else
                {
                  LOG(4, CLR_YELLOW, "      [lead pho#" << leadPhoIndex << "] no away-side jet passes Δφ/pt cuts (bestPt ≤ 0)"
                                                   << "  [passed dphi=" << nPassDphi
                                                   << ", passed pt=" << nPassPt
                                                   << ", minBackToBack=" << m_minBackToBack
                                                   << ", minJetPt=" << m_minJetPt << "]");
                }
            }
        }
        else
        {
          if (Verbosity() >= 5)
            LOG(5, CLR_BLUE, "      [processCandidates] no iso∧tight photon in this event → no jet matching / JES fills");
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
    }
    catch (const std::exception& e)
    {
      LOG(0, CLR_YELLOW, "    [processCandidates] EXCEPTION in photon path: " << e.what());
    }
    catch (...)
    {
      LOG(0, CLR_YELLOW, "    [processCandidates] UNKNOWN exception in photon path");
    }

    return; // keep original behavior: prefer photon path and return
  } // end photon path

  // =================== NO RawCluster fallback ===============================
  // We require PhotonClusterBuilder products for ALL photon kinematics.
  // If PHOTONCLUSTER_CEMC is missing, do NOT silently switch definitions.
  LOG(0, CLR_YELLOW,
        "  [processCandidates] PHOTONCLUSTER_CEMC is MISSING for this event. "
        "Refusing RawCluster fallback to prevent inconsistent photon kinematics. "
        "This usually means PhotonClusterBuilder was not run (or node name mismatch) in the macro/DST.");

  return;
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
      std::isfinite(v.pt_gamma);

  if (!ok_vals)
  {
    LOG(2, CLR_YELLOW,
        "  [passesPhotonPreselection] non-finite SSVars detected: "
        << "weta=" << v.weta_cogx << " wphi=" << v.wphi_cogx
        << " et1=" << v.et1 << " e11/e33=" << v.e11_over_e33
        << " e32/e35=" << v.e32_over_e35 << " pT^γ=" << v.pt_gamma);
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

  if (Verbosity() >= 4)
  {
      if (!pass_all)
      {
        LOG(4, CLR_YELLOW,
            "  [passesPhotonPreselection] FAILED: weta=" << v.weta_cogx
            << ", et1=" << v.et1
            << ", e11/e33=" << v.e11_over_e33
            << ", e32/e35=" << v.e32_over_e35);
      }
      else
      {
        LOG(4, CLR_RED,
            "  [passesPhotonPreselection] PASS: weta=" << v.weta_cogx
            << ", et1=" << v.et1
            << ", e11/e33=" << v.e11_over_e33
            << ", e32/e35=" << v.e32_over_e35);
      }
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
    const double w_hi = tight_w_hi(v.pt_gamma);
    if (!std::isfinite(w_hi))
    {
      LOG(2, CLR_YELLOW, "  [classifyPhotonTightness] non-finite tight_w_hi for pT^γ=" << v.pt_gamma
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
    {
      const char* tagName = tightTagName(tag);
      const char* colour  = (tag == TightTag::kTight) ? CLR_RED : CLR_GREEN;
      LOG(4, colour, "  [classifyPhotonTightness] n_fail=" << n_fail
                      << " → tag=" << tagName);
    }

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
      LOG(3, CLR_YELLOW, "  [eiso] clus==nullptr → return +inf (fail-safe)");
    return 1e9; // fail-safe: not isolated
  }

  // We require PhotonClusterBuilder-provided iso_* parameters stored on PhotonClusterv1.
  const auto* pho = dynamic_cast<const PhotonClusterv1*>(clus);
  if (!pho)
  {
    if (Verbosity() >= 2)
      LOG(2, CLR_YELLOW, "  [eiso] clus is not PhotonClusterv1 → return +inf (fail-safe)");
    return 1e9;
  }

  // PhotonClusterBuilder provides full (EMCal + HCALIN + HCALOUT) isolation for R=0.3 and R=0.4:
  //   iso_03_emcal (already subtracts ET^gamma), iso_03_hcalin, iso_03_hcalout
  //   iso_04_emcal (already subtracts ET^gamma), iso_04_hcalin, iso_04_hcalout
  const int cone10 = static_cast<int>(std::lround(10.0 * m_isoConeR));

  const char* k_em = nullptr;
  const char* k_hi = nullptr;
  const char* k_ho = nullptr;

  if (cone10 == 3)
  {
    k_em = "iso_03_emcal";
    k_hi = "iso_03_hcalin";
    k_ho = "iso_03_hcalout";
  }
  else if (cone10 == 4)
  {
    k_em = "iso_04_emcal";
    k_hi = "iso_04_hcalin";
    k_ho = "iso_04_hcalout";
  }
  else
  {
    if (Verbosity() >= 2)
      LOG(2, CLR_YELLOW,
          "  [eiso] m_isoConeR=" << m_isoConeR << " (cone10=" << cone10
          << ") not supported by PhotonClusterBuilder full calo iso → return +inf (fail-safe)");
    return 1e9;
  }

  const double em = pho->get_shower_shape_parameter(k_em);
  const double hi = pho->get_shower_shape_parameter(k_hi);
  const double ho = pho->get_shower_shape_parameter(k_ho);

  // Total isolation ET (PhotonClusterBuilder convention):
  //   iso_*_emcal already has (-ET^gamma) applied, HCAL terms are raw sums.
  const double iso = em + hi + ho;

  if (!std::isfinite(iso))
  {
    if (Verbosity() >= 2)
      LOG(2, CLR_YELLOW,
          "  [eiso] non-finite PhotonClusterBuilder iso: "
          << k_em << "=" << em << " " << k_hi << "=" << hi << " " << k_ho << "=" << ho
          << " → return +inf (fail-safe)");
    return 1e9;
  }

  if (Verbosity() >= 5)
  {
    LOG(5, CLR_BLUE,
        "  [eiso(builder)] cone10=" << cone10
        << "  " << k_em << "=" << em
        << "  " << k_hi << "=" << hi
        << "  " << k_ho << "=" << ho
        << "  total=" << iso);
  }

  return iso;
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
    const bool passIso = (eiso_val < thr);

    if (Verbosity() >= 5)
        LOG(5, (passIso ? CLR_RED : CLR_YELLOW),
            "  [isIsolated] ET^γ=" << et_gamma << "  thr=" << thr << "  eiso=" << eiso_val
            << "  pass=" << (passIso ? 1 : 0));

    return passIso;
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
int RecoilJets::findPtBin(double pt) const
{
  if (m_gammaPtBins.size() < 2)
  {
    if (Verbosity() >= 3)
      LOG(3, CLR_YELLOW, "  [findPtBin] pT bin edges vector has size < 2 – returning -1");
    return -1;
  }

  // Optional monotonicity check (diagnostic only)
  if (Verbosity() >= 6)
  {
    bool mono = std::is_sorted(m_gammaPtBins.begin(), m_gammaPtBins.end());
    if (!mono)
      LOG(6, CLR_YELLOW, "  [findPtBin] pT bin edges are NOT monotonic");
  }

  if (!std::isfinite(pt))
  {
    if (Verbosity() >= 4)
      LOG(4, CLR_YELLOW, "  [findPtBin] non-finite pT=" << pt << " – returning -1");
    return -1;
  }

  for (size_t i = 0; i + 1 < m_gammaPtBins.size(); ++i)
  {
    if (pt >= m_gammaPtBins[i] && pt < m_gammaPtBins[i + 1])
      return static_cast<int>(i);
  }

  if (Verbosity() >= 5)
    LOG(5, CLR_BLUE, "  [findPtBin] pT=" << pt << " out of configured range – returning -1");
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
std::string RecoilJets::suffixForBins(int ptIdx, int centIdx) const
{
  std::ostringstream s;

  // Photon pT slicing (uses m_gammaPtBins as the configured pT bin edges)
  if (ptIdx >= 0) {
    const double lo = m_gammaPtBins[ptIdx];
    const double hi = m_gammaPtBins[ptIdx+1];
    s << "_pT_" << std::fixed << std::setprecision(0) << lo << '_' << hi;
  }

  // Centrality slicing (Au+Au only, and only if centIdx>=0)
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

  // 60 bins in [0,3]
  const int    nbins = 60;
  const double xmin  = 0.0;
  const double xmax  = 3.0;

  // NOTE: now explicitly labeled as pTγ in the title
  const std::string title = name + ";x_{J}=p_{T}^{jet}/p_{T}^{#gamma};Entries";

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


// ------------------------------------------------------------------
// NEW: Leading recoil jet pT distribution (same slicing as h_xJ)
// ------------------------------------------------------------------
TH1F* RecoilJets::getOrBookJet1PtHist(const std::string& trig, int etIdx, int centIdx)
{
  const std::string base   = "h_jet1Pt";
  const std::string suffix = suffixForBins(etIdx, centIdx);
  const std::string name   = base + suffix;

  if (trig.empty()) return nullptr;

  auto& H = qaHistogramsByTrigger[trig];
  if (auto it = H.find(name); it != H.end())
    if (auto* h = dynamic_cast<TH1F*>(it->second)) return h;

  if (!out || !out->IsOpen()) return nullptr;

  TDirectory* const prevDir = gDirectory;
  TDirectory* dir = out->GetDirectory(trig.c_str());
  if (!dir) dir = out->mkdir(trig.c_str());
  if (!dir) { if (prevDir) prevDir->cd(); return nullptr; }

  dir->cd();

  const int    nb = 120;
  const double lo = 0.0;
  const double hi = 60.0;

  const std::string title = name + ";p_{T}^{jet1} [GeV];Entries";
  auto* h = new TH1F(name.c_str(), title.c_str(), nb, lo, hi);
  H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}


// ------------------------------------------------------------------
// NEW: Subleading jet pT distribution (same slicing as h_xJ)
// ------------------------------------------------------------------
TH1F* RecoilJets::getOrBookJet2PtHist(const std::string& trig, int etIdx, int centIdx)
{
  const std::string base   = "h_jet2Pt";
  const std::string suffix = suffixForBins(etIdx, centIdx);
  const std::string name   = base + suffix;

  if (trig.empty()) return nullptr;

  auto& H = qaHistogramsByTrigger[trig];
  if (auto it = H.find(name); it != H.end())
    if (auto* h = dynamic_cast<TH1F*>(it->second)) return h;

  if (!out || !out->IsOpen()) return nullptr;

  TDirectory* const prevDir = gDirectory;
  TDirectory* dir = out->GetDirectory(trig.c_str());
  if (!dir) dir = out->mkdir(trig.c_str());
  if (!dir) { if (prevDir) prevDir->cd(); return nullptr; }

  dir->cd();

  const int    nb = 120;
  const double lo = 0.0;
  const double hi = 60.0;

  const std::string title = name + ";p_{T}^{jet2} [GeV];Entries";
  auto* h = new TH1F(name.c_str(), title.c_str(), nb, lo, hi);
  H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}


// ------------------------------------------------------------------
// NEW: alpha = pTjet2 / pTgamma (same slicing as h_xJ)
// ------------------------------------------------------------------
TH1F* RecoilJets::getOrBookAlphaHist(const std::string& trig, int etIdx, int centIdx)
{
  const std::string base   = "h_alpha";
  const std::string suffix = suffixForBins(etIdx, centIdx);
  const std::string name   = base + suffix;

  if (trig.empty()) return nullptr;

  auto& H = qaHistogramsByTrigger[trig];
  if (auto it = H.find(name); it != H.end())
    if (auto* h = dynamic_cast<TH1F*>(it->second)) return h;

  if (!out || !out->IsOpen()) return nullptr;

  TDirectory* const prevDir = gDirectory;
  TDirectory* dir = out->GetDirectory(trig.c_str());
  if (!dir) dir = out->mkdir(trig.c_str());
  if (!dir) { if (prevDir) prevDir->cd(); return nullptr; }

  dir->cd();

  const int    nb = 60;
  const double lo = 0.0;
  const double hi = 2.0;

  const std::string title = name + ";#alpha=p_{T}^{jet2}/p_{T}^{#gamma};Entries";
  auto* h = new TH1F(name.c_str(), title.c_str(), nb, lo, hi);
  H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}


// ------------------------------------------------------------------
// NEW: TH3 for xJ vs alpha vs pTgamma (centrality suffix only)
// ------------------------------------------------------------------
TH3F* RecoilJets::getOrBookJES3_xJ_alphaHist(const std::string& trig, int centIdx)
{
  const std::string base   = "h_JES3_pT_xJ_alpha";
  const std::string suffix = suffixForBins(-1, centIdx); // cent-only (no pT suffix)
  const std::string name   = base + suffix;

  if (trig.empty()) return nullptr;

  auto& H = qaHistogramsByTrigger[trig];
  if (auto it = H.find(name); it != H.end())
  if (auto* h = dynamic_cast<TH3F*>(it->second)) return h;

  if (!out || !out->IsOpen()) return nullptr;
  if (m_gammaPtBins.size() < 2) return nullptr;

  TDirectory* const prevDir = gDirectory;
  TDirectory* dir = out->GetDirectory(trig.c_str());
  if (!dir) dir = out->mkdir(trig.c_str());
  if (!dir) { if (prevDir) prevDir->cd(); return nullptr; }

  dir->cd();

  const int nx = static_cast<int>(m_gammaPtBins.size()) - 1;
  const double* xbins = m_gammaPtBins.data();

  const int ny = 60;
  const double ylo = 0.0, yhi = 3.0;   // xJ

  const int nz = 40;
  const double zlo = 0.0, zhi = 2.0;   // alpha

  const std::string title =
    name + ";p_{T}^{#gamma} [GeV];x_{J}=p_{T}^{jet1}/p_{T}^{#gamma};#alpha=p_{T}^{jet2}/p_{T}^{#gamma}";

  auto* h = new TH3F(name.c_str(), title.c_str(),
                     nx, xbins,
                     ny, ylo, yhi,
                     nz, zlo, zhi);

  H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}


// ------------------------------------------------------------------
// NEW: TH3 for jet1Pt vs alpha vs pTgamma (centrality suffix only)
// ------------------------------------------------------------------
TH3F* RecoilJets::getOrBookJES3_jet1Pt_alphaHist(const std::string& trig, int centIdx)
{
  const std::string base   = "h_JES3_pT_jet1Pt_alpha";
  const std::string suffix = suffixForBins(-1, centIdx); // cent-only
  const std::string name   = base + suffix;

  if (trig.empty()) return nullptr;

  auto& H = qaHistogramsByTrigger[trig];
  if (auto it = H.find(name); it != H.end())
    if (auto* h = dynamic_cast<TH3F*>(it->second)) return h;

  if (!out || !out->IsOpen()) return nullptr;
  if (m_gammaPtBins.size() < 2) return nullptr;

  TDirectory* const prevDir = gDirectory;
  TDirectory* dir = out->GetDirectory(trig.c_str());
  if (!dir) dir = out->mkdir(trig.c_str());
  if (!dir) { if (prevDir) prevDir->cd(); return nullptr; }

  dir->cd();

  const int nx = static_cast<int>(m_gammaPtBins.size()) - 1;
  const double* xbins = m_gammaPtBins.data();

  const int ny = 120;
  const double ylo = 0.0, yhi = 60.0;  // jet1Pt

  const int nz = 40;
  const double zlo = 0.0, zhi = 2.0;   // alpha

  const std::string title =
    name + ";p_{T}^{#gamma} [GeV];p_{T}^{jet1} [GeV];#alpha=p_{T}^{jet2}/p_{T}^{#gamma}";

  auto* h = new TH3F(name.c_str(), title.c_str(),
                     nx, xbins,
                     ny, ylo, yhi,
                     nz, zlo, zhi);

  H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}


// ============================================================================
//  Matching-QA bookers (centrality suffix only; pT^gamma is an axis)
// ============================================================================

TH2F* RecoilJets::getOrBookMatchDphiVsPtGamma(const std::string& trig, int centIdx)
{
  const std::string base   = "h_match_dphi_vs_pTgamma";
  const std::string suffix = suffixForBins(-1, centIdx);   // centrality-only
  const std::string name   = base + suffix;

  if (trig.empty()) return nullptr;

  auto& H = qaHistogramsByTrigger[trig];
  if (auto it = H.find(name); it != H.end())
  {
    if (auto* h = dynamic_cast<TH2F*>(it->second)) return h;
    H.erase(it);
  }

  if (!out || !out->IsOpen()) return nullptr;
  if (m_gammaPtBins.size() < 2) return nullptr;

  TDirectory* const prevDir = gDirectory;
  TDirectory* dir = out->GetDirectory(trig.c_str());
  if (!dir) dir = out->mkdir(trig.c_str());
  if (!dir) { if (prevDir) prevDir->cd(); return nullptr; }

  dir->cd();

  const int nx = static_cast<int>(m_gammaPtBins.size()) - 1;
  const double* xbins = m_gammaPtBins.data();

  const int ny = 64;
  const double ylo = 0.0;
  const double yhi = M_PI;

  const std::string title =
    name + ";p_{T}^{#gamma} [GeV];|#Delta#phi(#gamma, recoil jet)| [rad]";

  auto* h = new TH2F(name.c_str(), title.c_str(), nx, xbins, ny, ylo, yhi);
  H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookMatchStatusVsPtGamma(const std::string& trig, int centIdx)
{
  const std::string base   = "h_match_status_vs_pTgamma";
  const std::string suffix = suffixForBins(-1, centIdx);   // centrality-only
  const std::string name   = base + suffix;

  if (trig.empty()) return nullptr;

  auto& H = qaHistogramsByTrigger[trig];
  if (auto it = H.find(name); it != H.end())
  {
    if (auto* h = dynamic_cast<TH2F*>(it->second)) return h;
    H.erase(it);
  }

  if (!out || !out->IsOpen()) return nullptr;
  if (m_gammaPtBins.size() < 2) return nullptr;

  TDirectory* const prevDir = gDirectory;
  TDirectory* dir = out->GetDirectory(trig.c_str());
  if (!dir) dir = out->mkdir(trig.c_str());
  if (!dir) { if (prevDir) prevDir->cd(); return nullptr; }

  dir->cd();

  const int nx = static_cast<int>(m_gammaPtBins.size()) - 1;
  const double* xbins = m_gammaPtBins.data();

  const int ny = 4;
  const double ylo = 0.5;
  const double yhi = 4.5;

  const std::string title =
    name + ";p_{T}^{#gamma} [GeV];match status";

  auto* h = new TH2F(name.c_str(), title.c_str(), nx, xbins, ny, ylo, yhi);

  h->GetYaxis()->SetBinLabel(1, "NoJetPt");
  h->GetYaxis()->SetBinLabel(2, "NoJetEta");
  h->GetYaxis()->SetBinLabel(3, "NoBackToBack");
  h->GetYaxis()->SetBinLabel(4, "Matched");

  H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookMatchMaxDphiVsPtGamma(const std::string& trig, int centIdx)
{
  const std::string base   = "h_match_maxdphi_vs_pTgamma";
  const std::string suffix = suffixForBins(-1, centIdx);   // centrality-only
  const std::string name   = base + suffix;

  if (trig.empty()) return nullptr;

  auto& H = qaHistogramsByTrigger[trig];
  if (auto it = H.find(name); it != H.end())
  {
    if (auto* h = dynamic_cast<TH2F*>(it->second)) return h;
    H.erase(it);
  }

  if (!out || !out->IsOpen()) return nullptr;
  if (m_gammaPtBins.size() < 2) return nullptr;

  TDirectory* const prevDir = gDirectory;
  TDirectory* dir = out->GetDirectory(trig.c_str());
  if (!dir) dir = out->mkdir(trig.c_str());
  if (!dir) { if (prevDir) prevDir->cd(); return nullptr; }

  dir->cd();

  const int nx = static_cast<int>(m_gammaPtBins.size()) - 1;
  const double* xbins = m_gammaPtBins.data();

  const int ny = 64;
  const double ylo = 0.0;
  const double yhi = M_PI;

  const std::string title =
    name + ";p_{T}^{#gamma} [GeV];max |#Delta#phi(#gamma, jet)| over jets passing p_{T}+|#eta| [rad]";

  auto* h = new TH2F(name.c_str(), title.c_str(), nx, xbins, ny, ylo, yhi);
  H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}

TProfile* RecoilJets::getOrBookNRecoilJetsVsPtGamma(const std::string& trig, int centIdx)
{
  const std::string base   = "p_nRecoilJets_vs_pTgamma";
  const std::string suffix = suffixForBins(-1, centIdx);   // centrality-only
  const std::string name   = base + suffix;

  if (trig.empty()) return nullptr;

  auto& H = qaHistogramsByTrigger[trig];
  if (auto it = H.find(name); it != H.end())
  {
    if (auto* p = dynamic_cast<TProfile*>(it->second)) return p;
    H.erase(it);
  }

  if (!out || !out->IsOpen()) return nullptr;
  if (m_gammaPtBins.size() < 2) return nullptr;

  TDirectory* const prevDir = gDirectory;
  TDirectory* dir = out->GetDirectory(trig.c_str());
  if (!dir) dir = out->mkdir(trig.c_str());
  if (!dir) { if (prevDir) prevDir->cd(); return nullptr; }

  dir->cd();

  const int nx = static_cast<int>(m_gammaPtBins.size()) - 1;
  const double* xbins = m_gammaPtBins.data();

  const std::string title =
    name + ";p_{T}^{#gamma} [GeV];#LT N_{recoil jets} #GT (jets passing p_{T}+|#eta|+#Delta#phi cut)";

  auto* p = new TProfile(name.c_str(), title.c_str(), nx, xbins, 0.0, 50.0);
  H[name] = p;

  if (prevDir) prevDir->cd();
  return p;
}

TH2F* RecoilJets::getOrBookRecoilIsLeadingVsPtGamma(const std::string& trig, int centIdx)
{
  const std::string base   = "h_recoilIsLeading_vs_pTgamma";
  const std::string suffix = suffixForBins(-1, centIdx);   // centrality-only
  const std::string name   = base + suffix;

  if (trig.empty()) return nullptr;

  auto& H = qaHistogramsByTrigger[trig];
  if (auto it = H.find(name); it != H.end())
  {
    if (auto* h = dynamic_cast<TH2F*>(it->second)) return h;
    H.erase(it);
  }

  if (!out || !out->IsOpen()) return nullptr;
  if (m_gammaPtBins.size() < 2) return nullptr;

  TDirectory* const prevDir = gDirectory;
  TDirectory* dir = out->GetDirectory(trig.c_str());
  if (!dir) dir = out->mkdir(trig.c_str());
  if (!dir) { if (prevDir) prevDir->cd(); return nullptr; }

  dir->cd();

  const int nx = static_cast<int>(m_gammaPtBins.size()) - 1;
  const double* xbins = m_gammaPtBins.data();

  const int ny = 2;
  const double ylo = -0.5;
  const double yhi =  1.5;

  const std::string title =
    name + ";p_{T}^{#gamma} [GeV];recoil jet is leading jet (p_{T}+|#eta| set)";

  auto* h = new TH2F(name.c_str(), title.c_str(), nx, xbins, ny, ylo, yhi);
  h->GetYaxis()->SetBinLabel(1, "NotLeading");
  h->GetYaxis()->SetBinLabel(2, "Leading");

  H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}


// ============================================================================
//  Your original getOrBookIsoHist (UNCHANGED) follows
// ============================================================================

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
  const int    nbins = 170;
  const double xmin  = -5.0;
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


// ------------------------------------------------------------------
// NEW: generic "iso-like-h_Eiso" booker for component isolation spectra
//      (EMCal / IHCal / OHCal), same binning & same slicing rules.
// ------------------------------------------------------------------
TH1F* RecoilJets::getOrBookIsoPartHist(const std::string& trig,
                                       const std::string& base,
                                       const std::string& xAxisTitle,
                                       int ptIdx, int centIdx)
{
  const std::string suffix = suffixForBins(ptIdx, centIdx);
  const std::string name   = base + suffix;

  if (Verbosity() >= 5)
    LOG(5, CLR_BLUE, "  [getOrBookIsoPartHist] trig=\"" << trig
           << "\" base=\"" << base << "\" ptIdx=" << ptIdx << " centIdx=" << centIdx
           << " → name=\"" << name << "\"");

  if (trig.empty() || base.empty())
  {
    LOG(2, CLR_YELLOW, "  [getOrBookIsoPartHist] empty trig/base – returning nullptr");
    return nullptr;
  }

  auto& H = qaHistogramsByTrigger[trig];

  if (auto it = H.find(name); it != H.end())
  {
    if (auto* h = dynamic_cast<TH1F*>(it->second)) return h;
    LOG(2, CLR_YELLOW, "    [getOrBookIsoPartHist] name clash: object \"" << name
                       << "\" exists but is not TH1F – replacing it");
    H.erase(it);
  }

  if (!out || !out->IsOpen())
  {
    LOG(1, CLR_YELLOW, "  [getOrBookIsoPartHist] output TFile invalid/null – returning nullptr");
    return nullptr;
  }

  TDirectory* const prevDir = gDirectory;
  TDirectory* dir = out->GetDirectory(trig.c_str());
  if (!dir) dir = out->mkdir(trig.c_str());
  if (!dir)
  {
    LOG(1, CLR_YELLOW, "  [getOrBookIsoPartHist] failed to create/access directory \"" << trig << "\"");
    if (prevDir) prevDir->cd();
    return nullptr;
  }

  dir->cd();

  // Match h_Eiso binning
  const int    nbins = 170;
  const double xmin  = -5.0;
  const double xmax  = 12.0;

  const std::string xlab  = (xAxisTitle.empty() ? "E_{T}^{iso} [GeV]" : xAxisTitle);
  const std::string title = name + ";" + xlab + ";Entries";

  auto* h = new TH1F(name.c_str(), title.c_str(), nbins, xmin, xmax);
  if (!h)
  {
    LOG(1, CLR_YELLOW, "  [getOrBookIsoPartHist] new TH1F failed for \"" << name << '"');
    if (prevDir) prevDir->cd();
    return nullptr;
  }

  H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}


// ------------------------------------------------------------------
// NEW: isolation PASS/FAIL counter histogram (2 bins), same slicing rules.
//   bin 1 = PASS  (Eiso < thr)
//   bin 2 = FAIL  (Eiso >= thr)
// ------------------------------------------------------------------
TH1I* RecoilJets::getOrBookIsoDecisionHist(const std::string& trig, int ptIdx, int centIdx)
{
  const std::string base   = "h_isoDecision";
  const std::string suffix = suffixForBins(ptIdx, centIdx);
  const std::string name   = base + suffix;

  if (Verbosity() >= 5)
    LOG(5, CLR_BLUE, "  [getOrBookIsoDecisionHist] trig=\"" << trig
           << "\" ptIdx=" << ptIdx << " centIdx=" << centIdx
           << " → name=\"" << name << "\"");

  if (trig.empty())
  {
    LOG(2, CLR_YELLOW, "  [getOrBookIsoDecisionHist] empty trig – returning nullptr");
    return nullptr;
  }

  auto& H = qaHistogramsByTrigger[trig];

  if (auto it = H.find(name); it != H.end())
  {
    if (auto* h = dynamic_cast<TH1I*>(it->second)) return h;
    LOG(2, CLR_YELLOW, "    [getOrBookIsoDecisionHist] name clash: object \"" << name
                       << "\" exists but is not TH1I – replacing it");
    H.erase(it);
  }

  if (!out || !out->IsOpen())
  {
    LOG(1, CLR_YELLOW, "  [getOrBookIsoDecisionHist] output TFile invalid/null – returning nullptr");
    return nullptr;
  }

  TDirectory* const prevDir = gDirectory;
  TDirectory* dir = out->GetDirectory(trig.c_str());
  if (!dir) dir = out->mkdir(trig.c_str());
  if (!dir)
  {
    LOG(1, CLR_YELLOW, "  [getOrBookIsoDecisionHist] failed to create/access directory \"" << trig << "\"");
    if (prevDir) prevDir->cd();
    return nullptr;
  }

  dir->cd();

  auto* h = new TH1I(name.c_str(), (name + ";Isolation decision;Entries").c_str(), 2, 0.5, 2.5);
  if (!h)
  {
    LOG(1, CLR_YELLOW, "  [getOrBookIsoDecisionHist] new TH1I failed for \"" << name << '"');
    if (prevDir) prevDir->cd();
    return nullptr;
  }

  h->GetXaxis()->SetBinLabel(1, "PASS");
  h->GetXaxis()->SetBinLabel(2, "FAIL");
  H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}



// ============================================================================
//  Jet QA helpers (generic bookers + fillers)
//  - Naming:  <base>_<rKey> + suffixForBins(ptIdx,centIdx)
//  - Pure jet QA uses ptIdx = -1 (centrality-only suffix in AuAu)
// ============================================================================

TH1I* RecoilJets::getOrBookJetQA1I(const std::string& trig,
                                  const std::string& base,
                                  const std::string& xAxisTitle,
                                  const std::string& rKey,
                                  int ptIdx, int centIdx,
                                  int nbins, double xmin, double xmax)
{
  const std::string suffix = suffixForBins(ptIdx, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || base.empty() || rKey.empty()) return nullptr;

  auto& H = qaHistogramsByTrigger[trig];

  if (auto it = H.find(name); it != H.end())
  {
    if (auto* h = dynamic_cast<TH1I*>(it->second)) return h;
    H.erase(it);
  }

  if (!out || !out->IsOpen()) return nullptr;

  TDirectory* const prevDir = gDirectory;
  TDirectory* dir = out->GetDirectory(trig.c_str());
  if (!dir) dir = out->mkdir(trig.c_str());
  if (!dir) { if (prevDir) prevDir->cd(); return nullptr; }

  dir->cd();

  const std::string xlab  = (xAxisTitle.empty() ? base : xAxisTitle);
  const std::string title = name + ";" + xlab + ";Entries";

  auto* h = new TH1I(name.c_str(), title.c_str(), nbins, xmin, xmax);
  H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}


TH1F* RecoilJets::getOrBookJetQA1F(const std::string& trig,
                                  const std::string& base,
                                  const std::string& xAxisTitle,
                                  const std::string& rKey,
                                  int ptIdx, int centIdx,
                                  int nbins, double xmin, double xmax)
{
  const std::string suffix = suffixForBins(ptIdx, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || base.empty() || rKey.empty()) return nullptr;

  auto& H = qaHistogramsByTrigger[trig];

  if (auto it = H.find(name); it != H.end())
  {
    if (auto* h = dynamic_cast<TH1F*>(it->second)) return h;
    H.erase(it);
  }

  if (!out || !out->IsOpen()) return nullptr;

  TDirectory* const prevDir = gDirectory;
  TDirectory* dir = out->GetDirectory(trig.c_str());
  if (!dir) dir = out->mkdir(trig.c_str());
  if (!dir) { if (prevDir) prevDir->cd(); return nullptr; }

  dir->cd();

  const std::string xlab  = (xAxisTitle.empty() ? base : xAxisTitle);
  const std::string title = name + ";" + xlab + ";Entries";

  auto* h = new TH1F(name.c_str(), title.c_str(), nbins, xmin, xmax);
  H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}


TH2F* RecoilJets::getOrBookJetQA2F(const std::string& trig,
                                  const std::string& base,
                                  const std::string& xAxisTitle,
                                  const std::string& yAxisTitle,
                                  const std::string& rKey,
                                  int ptIdx, int centIdx,
                                  int nxbins, double xmin, double xmax,
                                  int nybins, double ymin, double ymax)
{
  const std::string suffix = suffixForBins(ptIdx, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || base.empty() || rKey.empty()) return nullptr;

  auto& H = qaHistogramsByTrigger[trig];

  if (auto it = H.find(name); it != H.end())
  {
    if (auto* h = dynamic_cast<TH2F*>(it->second)) return h;
    H.erase(it);
  }

  if (!out || !out->IsOpen()) return nullptr;

  TDirectory* const prevDir = gDirectory;
  TDirectory* dir = out->GetDirectory(trig.c_str());
  if (!dir) dir = out->mkdir(trig.c_str());
  if (!dir) { if (prevDir) prevDir->cd(); return nullptr; }

  dir->cd();

  const std::string xlab  = (xAxisTitle.empty() ? "x" : xAxisTitle);
  const std::string ylab  = (yAxisTitle.empty() ? "y" : yAxisTitle);
  const std::string title = name + ";" + xlab + ";" + ylab;

  auto* h = new TH2F(name.c_str(), title.c_str(),
                     nxbins, xmin, xmax,
                     nybins, ymin, ymax);
  H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}


// ---------------------------- Pure jet QA filler ----------------------------
void RecoilJets::fillInclusiveJetQA(const std::vector<std::string>& activeTrig,
                                   int centIdx,
                                   const std::string& rKey)
{
  if (activeTrig.empty()) return;

  JetContainer* jets = nullptr;
  if (auto it = m_jets.find(rKey); it != m_jets.end()) jets = it->second;
  if (!jets) return;

  const int effCentIdx = (m_isAuAu ? centIdx : -1);

  // Event-level (fiducial jets only)
  int    nJetsFid = 0;
  double HT       = 0.0;
  double leadPt   = -1.0;
  double leadEta  = 0.0;
  double leadPhi  = 0.0;
  double subPt    = -1.0;

  constexpr int    nbPt  = 200;  constexpr double ptLo  = 0.0;   constexpr double ptHi  = 100.0;
  constexpr int    nbEta = 120;  constexpr double etaLo = -1.2;  constexpr double etaHi =  1.2;
  constexpr int    nbPhi = 128;  constexpr double phiLo = -M_PI; constexpr double phiHi =  M_PI;

  constexpr int    nbMass = 120; constexpr double mLo = 0.0;     constexpr double mHi = 30.0;

  // 2D occupancy and mass-vs-pt
  constexpr int    nbEta2 = 60;
  constexpr int    nbPhi2 = 64;
  constexpr int    nbPt2  = 100;
  constexpr int    nbM2   = 80;
  constexpr double m2Hi   = 20.0;

  // Loop over jets once; fill "all" first (pT cut only), then "incl" (fiducial)
  for (const Jet* j : *jets)
  {
    if (!j) continue;

    const double pt  = j->get_pt();
    const double eta = j->get_eta();

    if (!std::isfinite(pt) || !std::isfinite(eta)) continue;
    if (pt < m_minJetPt) continue;

    const double phi = TVector2::Phi_mpi_pi(j->get_phi());
    const double mass = j->get_mass();

    // ----------------- "all": pT cut only (no eta acceptance) -----------------
    for (const auto& trigShort : activeTrig)
    {
      if (auto* h = getOrBookJetQA1F(trigShort, "h_jetPt_all", "p_{T}^{jet} [GeV]", rKey, -1, effCentIdx, nbPt, ptLo, ptHi))
      { h->Fill(pt); bumpHistFill(trigShort, h->GetName()); }

      if (auto* h = getOrBookJetQA1F(trigShort, "h_jetEta_all", "#eta^{jet}", rKey, -1, effCentIdx, nbEta, etaLo, etaHi))
      { h->Fill(eta); bumpHistFill(trigShort, h->GetName()); }

      if (auto* h = getOrBookJetQA1F(trigShort, "h_jetPhi_all", "#phi^{jet} [rad]", rKey, -1, effCentIdx, nbPhi, phiLo, phiHi))
      { h->Fill(phi); bumpHistFill(trigShort, h->GetName()); }

      if (auto* h2 = getOrBookJetQA2F(trigShort, "h_jetEtaPhi_all", "#eta^{jet}", "#phi^{jet} [rad]",
                                      rKey, -1, effCentIdx, nbEta2, etaLo, etaHi, nbPhi2, phiLo, phiHi))
      { h2->Fill(eta, phi); bumpHistFill(trigShort, h2->GetName()); }

      if (std::isfinite(mass) && mass >= 0.0)
      {
        if (auto* hm = getOrBookJetQA1F(trigShort, "h_jetMass_all", "m^{jet} [GeV]", rKey, -1, effCentIdx, nbMass, mLo, mHi))
        { hm->Fill(mass); bumpHistFill(trigShort, hm->GetName()); }

        if (auto* h2m = getOrBookJetQA2F(trigShort, "h_jetMassVsPt_all", "p_{T}^{jet} [GeV]", "m^{jet} [GeV]",
                                         rKey, -1, effCentIdx, nbPt2, ptLo, ptHi, nbM2, mLo, m2Hi))
        { h2m->Fill(pt, mass); bumpHistFill(trigShort, h2m->GetName()); }
      }
    }

    // ----------------- "incl": fiducial jets used for analysis ----------------
    if (std::fabs(eta) >= m_jetEtaAbsMax) continue;

    ++nJetsFid;
    HT += pt;

    if (pt > leadPt)
    {
      subPt   = leadPt;
      leadPt  = pt;
      leadEta = eta;
      leadPhi = phi;
    }
    else if (pt > subPt)
    {
      subPt = pt;
    }

    for (const auto& trigShort : activeTrig)
    {
      if (auto* h = getOrBookJetQA1F(trigShort, "h_jetPt_incl", "p_{T}^{jet} [GeV]", rKey, -1, effCentIdx, nbPt, ptLo, ptHi))
      { h->Fill(pt); bumpHistFill(trigShort, h->GetName()); }

      if (auto* h = getOrBookJetQA1F(trigShort, "h_jetEta_incl", "#eta^{jet}", rKey, -1, effCentIdx, nbEta, etaLo, etaHi))
      { h->Fill(eta); bumpHistFill(trigShort, h->GetName()); }

      if (auto* h = getOrBookJetQA1F(trigShort, "h_jetPhi_incl", "#phi^{jet} [rad]", rKey, -1, effCentIdx, nbPhi, phiLo, phiHi))
      { h->Fill(phi); bumpHistFill(trigShort, h->GetName()); }

      if (auto* h2 = getOrBookJetQA2F(trigShort, "h_jetEtaPhi_incl", "#eta^{jet}", "#phi^{jet} [rad]",
                                      rKey, -1, effCentIdx, nbEta2, etaLo, etaHi, nbPhi2, phiLo, phiHi))
      { h2->Fill(eta, phi); bumpHistFill(trigShort, h2->GetName()); }

      if (std::isfinite(mass) && mass >= 0.0)
      {
        if (auto* hm = getOrBookJetQA1F(trigShort, "h_jetMass_incl", "m^{jet} [GeV]", rKey, -1, effCentIdx, nbMass, mLo, mHi))
        { hm->Fill(mass); bumpHistFill(trigShort, hm->GetName()); }

        if (auto* h2m = getOrBookJetQA2F(trigShort, "h_jetMassVsPt_incl", "p_{T}^{jet} [GeV]", "m^{jet} [GeV]",
                                         rKey, -1, effCentIdx, nbPt2, ptLo, ptHi, nbM2, mLo, m2Hi))
        { h2m->Fill(pt, mass); bumpHistFill(trigShort, h2m->GetName()); }
      }
    }
  }

  // Event-level histograms (fiducial jets only)
  for (const auto& trigShort : activeTrig)
  {
    if (auto* h = getOrBookJetQA1I(trigShort, "h_nJets", "N_{jets}", rKey, -1, effCentIdx, 100, 0.0, 100.0))
    { h->Fill(nJetsFid); bumpHistFill(trigShort, h->GetName()); }

    if (auto* h = getOrBookJetQA1F(trigShort, "h_HT", "H_{T} [GeV]", rKey, -1, effCentIdx, 300, 0.0, 150.0))
    { h->Fill(HT); bumpHistFill(trigShort, h->GetName()); }

    if (leadPt >= 0.0)
    {
      if (auto* h = getOrBookJetQA1F(trigShort, "h_leadJetPt", "p_{T}^{lead jet} [GeV]", rKey, -1, effCentIdx, 200, 0.0, 100.0))
      { h->Fill(leadPt); bumpHistFill(trigShort, h->GetName()); }

      if (auto* h = getOrBookJetQA1F(trigShort, "h_leadJetEta", "#eta^{lead jet}", rKey, -1, effCentIdx, 120, -1.2, 1.2))
      { h->Fill(leadEta); bumpHistFill(trigShort, h->GetName()); }

      if (auto* h = getOrBookJetQA1F(trigShort, "h_leadJetPhi", "#phi^{lead jet} [rad]", rKey, -1, effCentIdx, 128, -M_PI, M_PI))
      { h->Fill(leadPhi); bumpHistFill(trigShort, h->GetName()); }
    }

    if (subPt >= 0.0)
    {
      if (auto* h = getOrBookJetQA1F(trigShort, "h_subleadJetPt", "p_{T}^{sublead jet} [GeV]", rKey, -1, effCentIdx, 200, 0.0, 100.0))
      { h->Fill(subPt); bumpHistFill(trigShort, h->GetName()); }
    }
  }
}


// ----------------------- Selected jet QA (gamma+jet jets) -----------------------
void RecoilJets::fillSelectedJetQA(const std::vector<std::string>& activeTrig,
                                  int ptIdx, int centIdx,
                                  const std::string& rKey,
                                  const Jet* jet1,
                                  const Jet* jet2)
{
  if (activeTrig.empty() || !jet1) return;

  const int effCentIdx = (m_isAuAu ? centIdx : -1);

  const double eta1 = jet1->get_eta();
  const double phi1 = TVector2::Phi_mpi_pi(jet1->get_phi());
  const double m1   = jet1->get_mass();

  for (const auto& trigShort : activeTrig)
  {
    if (auto* h = getOrBookJetQA1F(trigShort, "h_jet1Eta_sel", "#eta^{jet1}", rKey, ptIdx, effCentIdx, 120, -1.2, 1.2))
    { h->Fill(eta1); bumpHistFill(trigShort, h->GetName()); }

    if (auto* h = getOrBookJetQA1F(trigShort, "h_jet1Phi_sel", "#phi^{jet1} [rad]", rKey, ptIdx, effCentIdx, 128, -M_PI, M_PI))
    { h->Fill(phi1); bumpHistFill(trigShort, h->GetName()); }

    if (std::isfinite(m1) && m1 >= 0.0)
    {
      if (auto* h = getOrBookJetQA1F(trigShort, "h_jet1Mass_sel", "m^{jet1} [GeV]", rKey, ptIdx, effCentIdx, 120, 0.0, 30.0))
      { h->Fill(m1); bumpHistFill(trigShort, h->GetName()); }
    }

    if (jet2)
    {
      const double eta2 = jet2->get_eta();
      const double phi2 = TVector2::Phi_mpi_pi(jet2->get_phi());

      if (auto* h = getOrBookJetQA1F(trigShort, "h_jet2Eta_sel", "#eta^{jet2}", rKey, ptIdx, effCentIdx, 120, -1.2, 1.2))
      { h->Fill(eta2); bumpHistFill(trigShort, h->GetName()); }

      if (auto* h = getOrBookJetQA1F(trigShort, "h_jet2Phi_sel", "#phi^{jet2} [rad]", rKey, ptIdx, effCentIdx, 128, -M_PI, M_PI))
      { h->Fill(phi2); bumpHistFill(trigShort, h->GetName()); }
    }
  }
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

  const int ptIdx = findPtBin(pt_gamma);
  if (ptIdx < 0)
  {
    if (Verbosity() >= 4)
      LOG(4, CLR_YELLOW, "  [fillIsoSSTagCounters] ET bin not found for ET^γ=" << pt_gamma
                          << " – skipping fills");
    return;
  }

  const int effCentIdx = (m_isAuAu ? centIdx : -1);
  const std::string slice = suffixForBins(ptIdx, effCentIdx);

  // --- Isolation decision (Eiso histogram is filled unconditionally in the photon loop) ---
  const double eiso_et = eiso(clus, topNode);
  const double thrIso  = (m_isoA + m_isoB * et_gamma);
  const bool   iso     = (eiso_et < thrIso);


  // --- Tight sub-cuts (for printing & category) ---
  const double w_hi = tight_w_hi(v.pt_gamma);
  const bool pass_weta   = in_open_interval(v.weta_cogx,  TIGHT_W_LO, w_hi);
  const bool pass_wphi   = in_open_interval(v.wphi_cogx,  TIGHT_W_LO, w_hi);
  const bool pass_e11e33 = in_open_interval(v.e11_over_e33, TIGHT_E11E33_MIN, TIGHT_E11E33_MAX);
  const bool pass_et1    = in_open_interval(v.et1,          TIGHT_ET1_MIN,    TIGHT_ET1_MAX);
  const bool pass_e32e35 = in_open_interval(v.e32_over_e35, TIGHT_E32E35_MIN, TIGHT_E32E35_MAX);

  const int tight_fails  = (!pass_weta) + (!pass_wphi) + (!pass_e11e33) + (!pass_et1) + (!pass_e32e35);
  TightTag tag;
  if (tight_fails == 0)      tag = TightTag::kTight;
  else if (tight_fails >= 2) tag = TightTag::kNonTight;
  else                       tag = TightTag::kNeither;

  // --- Update per-slice counters (bookkeeping only) ---
  auto& S = m_catByTrig[trig][slice];
  S.seen += 1;

  // --- Choose 2×2 category and fill its COUNT histogram ---
  const char* comboBase = nullptr;     // base name for count histogram
  const char* comboKeySS = nullptr;    // tagKey for SS histograms
  if (iso && tag == TightTag::kTight)           { comboBase = "h_isIsolated_isTight";    comboKeySS = "isIsolated_isTight";    ++S.n_iso_tight; }
  else if (!iso && tag == TightTag::kTight)     { comboBase = "h_notIsolated_isTight";   comboKeySS = "notIsolated_isTight";   ++S.n_nonIso_tight; }
  else if (iso && tag != TightTag::kTight)      { comboBase = "h_isIsolated_notTight";   comboKeySS = "isIsolated_notTight";   ++S.n_iso_nonTight; }
  else                                           { comboBase = "h_notIsolated_notTight"; comboKeySS = "notIsolated_notTight"; ++S.n_nonIso_nonTight; }

  if (auto* hc = getOrBookCountHist(trig, comboBase, ptIdx, effCentIdx))
  {
    hc->Fill(1);
    bumpHistFill(trig, std::string(comboBase) + slice);
  }
  else
  {
    LOG(2, CLR_YELLOW, "  [fillIsoSSTagCounters] getOrBookCountHist returned nullptr for \"" << comboBase
                       << "\" slice=\"" << slice << '"');
  }

  if (Verbosity() >= 5)
  {
      std::ostringstream os;
      os << "    [fill] trig=\"" << trig << "\" slice=\"" << slice << "\" → "
         << comboBase
         << " | Eiso=" << std::fixed << std::setprecision(3) << eiso_et
         << " thr=" << std::fixed << std::setprecision(3) << thrIso
         << " iso=" << (iso ? "PASS" : "FAIL")
         << " | tight=" << tightTagName(tag)
         << " [weta "    << (pass_weta   ? "PASS" : "FAIL")
         << ", wphi "    << (pass_wphi   ? "PASS" : "FAIL")
         << ", e11/e33 " << (pass_e11e33 ? "PASS" : "FAIL")
         << ", et1 "     << (pass_et1    ? "PASS" : "FAIL")
         << ", e32/e35 " << (pass_e32e35 ? "PASS" : "FAIL")
         << "]";
      const char* colour = (iso || tag == TightTag::kTight) ? CLR_RED : CLR_CYAN;
      LOG(5, colour, os.str());
  }


  // --- SS variable histograms for THIS category (always) ---
  auto fillSS = [&](const std::string& key, double val)
  {
    if (auto* h = getOrBookSSHist(trig, key, comboKeySS, ptIdx, effCentIdx))
    {
      h->Fill(val);
      bumpHistFill(trig, "h_ss_" + key + "_" + comboKeySS + slice);
    }
    else
    {
      LOG(2, CLR_YELLOW, "    [fillSS] getOrBookSSHist returned nullptr for key=\"" << key
          << "\" tagKey=\"" << comboKeySS << "\" slice=\"" << slice << '"');
    }
  };
  fillSS("weta",   v.weta_cogx);
  fillSS("wphi",   v.wphi_cogx);
  fillSS("et1",    v.et1);
  fillSS("e11e33", v.e11_over_e33);
  fillSS("e32e35", v.e32_over_e35);

  // --- Optional: ID sideband bookkeeping (unchanged behavior) ---
  if (tag == TightTag::kNonTight)
  {
    S.idSB_total += 1;
    if (auto* h = getOrBookCountHist(trig, "h_idSB_total", ptIdx, effCentIdx))
    {
      h->Fill(1);
      bumpHistFill(trig, "h_idSB_total" + slice);
    }
    if (iso)
    {
      S.idSB_pass += 1;
      if (auto* h = getOrBookCountHist(trig, "h_idSB_pass", ptIdx, effCentIdx))
      {
        h->Fill(1);
        bumpHistFill(trig, "h_idSB_pass" + slice);
      }
    }
  }
}


void RecoilJets::printCutSummary() const
{
  using std::cout; using std::endl; using std::setw; using std::left; using std::right;

  cout << "\n\033[1mJob-wide cutflow summary\033[0m\n";

  // Events
  cout << "\n\033[1mEvents\033[0m\n";
  cout << "metric              | value\n"
          "--------------------+-------\n";
  cout << left << setw(20) << "seen"            << " | " << right << setw(7) << m_bk.evt_seen         << '\n';
  cout << left << setw(20) << "fail(trigger)"   << " | " << right << setw(7) << m_bk.evt_fail_trigger << '\n';
  cout << left << setw(20) << "fail(|vz|)"      << " | " << right << setw(7) << m_bk.evt_fail_vz      << '\n';
  cout << left << setw(20) << "accepted"        << " | " << right << setw(7) << m_bk.evt_accepted     << '\n';

  // Photon candidates: high-level
  cout << "\n\033[1mPhoton candidates (PHOTONCLUSTER_CEMC)\033[0m\n";
  cout << "metric                         | value\n"
          "-------------------------------+-------\n";
  cout << left << setw(31) << "total (container)"          << " | " << right << setw(7) << m_bk.pho_total           << '\n';
  cout << left << setw(31) << "early reject: E<2 GeV"      << " | " << right << setw(7) << m_bk.pho_early_E         << '\n';
  cout << left << setw(31) << "fail |eta| fiducial"        << " | " << right << setw(7) << m_bk.pho_eta_fail        << '\n';
  cout << left << setw(31) << "fail ET bin (out of range)" << " | " << right << setw(7) << m_bk.pho_etbin_out       << '\n';
  cout << left << setw(31) << "no RawCluster ptr"          << " | " << right << setw(7) << m_bk.pho_noRC            << '\n';
  cout << left << setw(31) << "reached pre+iso step"       << " | " << right << setw(7) << m_bk.pho_reached_pre_iso << '\n';

  // Preselection breakdown
  cout << "\n\033[1mPreselection (among reached pre+iso)\033[0m\n";
  cout << "outcome/criterion              | count\n"
          "-------------------------------+-------\n";
  cout << left << setw(31) << "pass preselection"          << " | " << right << setw(7) << m_bk.pre_pass             << '\n';
  cout << left << setw(31) << "fail weta >= 0.60"          << " | " << right << setw(7) << m_bk.pre_fail_weta        << '\n';
  cout << left << setw(31) << "fail et1 < 0.60"            << " | " << right << setw(7) << m_bk.pre_fail_et1_low     << '\n';
  cout << left << setw(31) << "fail et1 > 1.00"            << " | " << right << setw(7) << m_bk.pre_fail_et1_high    << '\n';
  cout << left << setw(31) << "fail e11/e33 >= 0.98"       << " | " << right << setw(7) << m_bk.pre_fail_e11e33_high << '\n';
  cout << left << setw(31) << "fail e32/e35 < 0.80"        << " | " << right << setw(7) << m_bk.pre_fail_e32e35_low  << '\n';
  cout << left << setw(31) << "fail e32/e35 > 1.00"        << " | " << right << setw(7) << m_bk.pre_fail_e32e35_high << '\n';

  // Tight classification (for preselection-passed)
  cout << "\n\033[1mTight classification (preselection-passed)\033[0m\n";
  cout << "category                       | count\n"
          "-------------------------------+-------\n";
  cout << left << setw(31) << "Tight (all 5 pass)"         << " | " << right << setw(7) << m_bk.tight_tight     << '\n';
  cout << left << setw(31) << "Neither (exactly 1 fail)"   << " | " << right << setw(7) << m_bk.tight_neither   << '\n';
  cout << left << setw(31) << "Non-Tight (>=2 fail)"       << " | " << right << setw(7) << m_bk.tight_nonTight  << '\n';

  cout << "— per-cut fails —\n";
  cout << left << setw(31) << "fail weta tight"            << " | " << right << setw(7) << m_bk.tight_fail_weta   << '\n';
  cout << left << setw(31) << "fail wphi tight"            << " | " << right << setw(7) << m_bk.tight_fail_wphi   << '\n';
  cout << left << setw(31) << "fail et1 tight"             << " | " << right << setw(7) << m_bk.tight_fail_et1    << '\n';
  cout << left << setw(31) << "fail e11/e33 tight"         << " | " << right << setw(7) << m_bk.tight_fail_e11e33 << '\n';
  cout << left << setw(31) << "fail e32/e35 tight"         << " | " << right << setw(7) << m_bk.tight_fail_e32e35 << '\n';

  // Isolation
  cout << "\n\033[1mIsolation (preselection-passed)\033[0m\n";
  cout << "outcome                        | count\n"
          "-------------------------------+-------\n";
  cout << left << setw(31) << "isolated (pass)"            << " | " << right << setw(7) << m_bk.iso_pass         << '\n';
  cout << left << setw(31) << "non-isolated (fail)"        << " | " << right << setw(7) << m_bk.iso_fail         << '\n';
}
