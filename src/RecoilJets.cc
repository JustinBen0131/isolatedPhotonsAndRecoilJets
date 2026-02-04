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
#include <TProfile3D.h>
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
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4eval/CaloRawClusterEval.h>

#include <calobase/RawTowerDefs.h>
#include <ffaobjects/EventHeader.h>

#if defined(__GNUC__) && !defined(__clang__)
  #pragma GCC diagnostic push
  // HepMC2 headers trigger deprecated std::iterator warnings under modern libstdc++
  #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>

#if defined(__GNUC__) && !defined(__clang__)
  #pragma GCC diagnostic pop
#endif
#include <calobase/TowerInfoDefs.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/RawCluster.h>
#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/coresoftware_local/offline/packages/CaloBase/PhotonClusterv1.h"
#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/coresoftware_local/offline/packages/CaloReco/PhotonClusterBuilder.h"
#include <calobase/RawTowerGeomContainer.h>
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
#include <limits>
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

// ============================================================================
// Jet-radius helpers (single source of truth for acceptance cuts)
//   rKey format: "r02" -> R=0.2, "r04" -> R=0.4
//   Containment/fiducial policy: |eta_jet| < 1.1 - R
// ============================================================================
namespace
{
  inline double jetRFromKey(const std::string& rKey)
  {
    // Expect "r02", "r04", etc.
    if (rKey.size() >= 3 && (rKey[0] == 'r' || rKey[0] == 'R'))
    {
      try
      {
        const double d = std::stod(rKey.substr(1));   // "02"->2, "04"->4
        return 0.1 * d;                               // 0.2, 0.4
      }
      catch (...) {}
    }
    return std::numeric_limits<double>::quiet_NaN();
  }

  inline double jetEtaAbsMaxForRKey(const std::string& rKey)
  {
    const double R = jetRFromKey(rKey);
    if (!std::isfinite(R)) return 0.0;

    // sPHENIX calorimeter acceptance: |eta| < 1.1
    const double etaMax = 1.1 - R;
    return (etaMax > 0.0 ? etaMax : 0.0);
  }
}


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
  if (!top)
  {
    LOG(0, CLR_YELLOW, "  [fetchNodes] top == nullptr → skip");
    return false;
  }

  // Small local helpers (keep this file self-contained)
  auto trim = [](std::string s) -> std::string
  {
    const char* ws = " \t\r\n";
    s.erase(0, s.find_first_not_of(ws));
    s.erase(s.find_last_not_of(ws) + 1);
    return s;
  };

  auto toLower = [](std::string s) -> std::string
  {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return s;
  };

  // ------------------------------------------------------------------
  // Vertex maps
  // ------------------------------------------------------------------
  GlobalVertexMap* gvmap  = findNode::getClass<GlobalVertexMap>(top, "GlobalVertexMap");
  MbdVertexMap*    mbdmap = findNode::getClass<MbdVertexMap>(top, "MbdVertexMap");

  // ------------------------------------------------------------------
  // Dataset / sim mode determination
  //   - start from module flags (setDataType("isSim") sets m_isSim)
  //   - allow wrapper env override (RJ_DATASET / RJ_IS_SIM)
  // ------------------------------------------------------------------
  bool isSim  = m_isSim;
  bool isAuAu = m_isAuAu;

  if (const char* ds = std::getenv("RJ_DATASET"))
  {
    const std::string s = toLower(trim(std::string(ds)));

    if (s == "issim" || s == "sim")
    {
      isSim  = true;
      isAuAu = false;
    }
    else if (s == "isauau" || s == "auau" || s == "aa")
    {
      isSim  = false;
      isAuAu = true;
    }
    else if (s == "ispp" || s == "pp")
    {
      isSim  = false;
      isAuAu = false;
    }
  }

  if (const char* f = std::getenv("RJ_IS_SIM"))
  {
    // RJ_IS_SIM=1 forces sim; RJ_IS_SIM=0 forces "not sim"
    const bool envSim = (std::atoi(f) != 0);
    if (envSim)
    {
      isSim  = true;
      isAuAu = false;
    }
    else
    {
      isSim = false;
      // Do NOT force isAuAu here — RJ_DATASET should set that if needed.
    }
  }

  // Keep internal flags consistent with what this event will do
  m_isSim  = isSim;
  m_isAuAu = isAuAu;

  // ------------------------------------------------------------------
  // Defaults
  // ------------------------------------------------------------------
  m_vtx = nullptr;
  m_vx  = 0.0;
  m_vy  = 0.0;
  m_vz  = 0.0;

  double gv_z    = std::numeric_limits<double>::quiet_NaN();
  double mbd_z   = std::numeric_limits<double>::quiet_NaN();
  double truth_z = std::numeric_limits<double>::quiet_NaN();

  bool haveGVZ    = false;
  bool haveMBDZ   = false;
  bool haveTruthZ = false;

  // ------------------------------------------------------------------
  // SIM ONLY: truth vertex from G4TruthInfo (requires paired G4Hits DST)
  // ------------------------------------------------------------------
  if (isSim)
  {
    auto* truth = findNode::getClass<PHG4TruthInfoContainer>(top, "G4TruthInfo");
    if (!truth)
    {
      LOG(3, CLR_YELLOW,
          "    [fetchNodes] isSim: G4TruthInfo missing (did you load paired G4Hits DST?)");
    }
    else
    {
      auto vr = truth->GetPrimaryVtxRange();
      if (vr.first != vr.second && vr.first->second)
      {
        const PHG4VtxPoint* vtx = vr.first->second;
        truth_z = vtx->get_z();
        haveTruthZ = std::isfinite(truth_z);
      }
      else
      {
        LOG(3, CLR_YELLOW, "    [fetchNodes] isSim: G4TruthInfo has no primary vertex");
      }
    }
  }

  // QA-only: we may still compare truth_z vs reco_z elsewhere; never use truth_z as reco vertex
  (void) truth_z;
  (void) haveTruthZ;

  // GlobalVertex (keep pointer + x/y)
  if (gvmap && !gvmap->empty())
  {
    m_vtx = gvmap->begin()->second;
    if (m_vtx)
    {
      m_vx  = m_vtx->get_x();
      m_vy  = m_vtx->get_y();
      gv_z  = m_vtx->get_z();
      haveGVZ = std::isfinite(gv_z);
    }
  }

  // MBD z (preferred for consistency with PhotonClusterBuilder)
  if (mbdmap && !mbdmap->empty())
  {
    MbdVertex* mv = mbdmap->begin()->second;
    if (mv)
    {
      mbd_z    = mv->get_z();
      haveMBDZ = std::isfinite(mbd_z);
    }
  }

  // Choose the z we will USE (RECO vertex ONLY):
  //   • SIM + data: MBD → else Global → else SKIP
  //   • Never use TRUTH as a fallback for reco objects
  const char* vz_source = "none";

  if (haveMBDZ)
  {
    m_vz = static_cast<float>(mbd_z);
    vz_source = "MBD";
  }
  else if (haveGVZ)
  {
    m_vz = static_cast<float>(gv_z);
    vz_source = "Global";
  }
  else
  {
    LOG(1, CLR_YELLOW, "  [fetchNodes] no usable RECO vertex z (MBD/Global) → skip event");
    return false;
  }

  // Print comparison ONLY if it matters (data only)
  if (!isSim && haveMBDZ && haveGVZ)
  {
    const double dz = mbd_z - gv_z;
    constexpr double kDzPrint = 1e-3;

    if (std::fabs(dz) > kDzPrint || Verbosity() >= 7)
    {
      LOG(3, CLR_BLUE,
          "  – vertex compare: vz(MBD)=" << std::fixed << std::setprecision(3) << mbd_z
          << "  vz(Global)=" << std::fixed << std::setprecision(3) << gv_z
          << "  Δz(MBD-Global)=" << std::fixed << std::setprecision(3) << dz);
    }
  }

  LOG(4, CLR_BLUE,
      "    [fetchNodes] dataset=" << (isSim ? "SIM" : (isAuAu ? "Au+Au" : "p+p"))
      << "  vz(used)=" << std::fixed << std::setprecision(2) << m_vz
      << "  (source=" << vz_source << ")");

  // ------------------------------------------------------------------
  // SIM QA (PRE-CUT): truth vs reco-used vertex z
  //   - Fill here so it is recorded even if the event later fails |vz| cut
  //   - X = truth vz (from G4TruthInfo)
  //   - Y = reco-used vz (your chosen m_vz: MBD → Global → TRUTH fallback)
  // ------------------------------------------------------------------
  if (isSim && haveTruthZ)
  {
    auto itTrig = qaHistogramsByTrigger.find("SIM");
    if (itTrig != qaHistogramsByTrigger.end())
    {
      auto& H = itTrig->second;
      auto itH2 = H.find("h_vzTruthVsReco");
      if (itH2 != H.end())
      {
        static_cast<TH2F*>(itH2->second)->Fill(static_cast<float>(truth_z), m_vz);
        bumpHistFill("SIM", "h_vzTruthVsReco");
      }
    }
  }

  // ------------------------------------------------------------------
  // Calo towers & geometry
  // ------------------------------------------------------------------
  m_calo.clear();
  for (const auto& ci : m_caloInfo)
  {
    const std::string node = std::get<0>(ci);
    const std::string geo  = std::get<1>(ci);
    const std::string lbl  = std::get<2>(ci);

    auto* tw = findNode::getClass<TowerInfoContainer>(top, node);
    auto* ge = findNode::getClass<RawTowerGeomContainer>(top, geo);
    if (!tw || !ge)
    {
      LOG(2, CLR_YELLOW, "  – missing " << lbl << " nodes → skip");
      return false;
    }
    m_calo[lbl] = { tw, ge, 0.0 };
    LOG(5, CLR_BLUE, "    [fetchNodes] towers OK: " << lbl << "  node=" << node << "  geom=" << geo);
  }

  // ------------------------------------------------------------------
  // Clusters & photons
  // ------------------------------------------------------------------
  m_clus    = findNode::getClass<RawClusterContainer>(top, "CLUSTERINFO_CEMC");
  m_photons = findNode::getClass<RawClusterContainer>(top, "PHOTONCLUSTER_CEMC");

  if (Verbosity() >= 3)
  {
    auto countClusters = [](RawClusterContainer* c) -> size_t
    {
      if (!c) return 0;
      size_t n = 0;
      auto r = c->getClusters();
      for (auto it = r.first; it != r.second; ++it) ++n;
      return n;
    };

    LOG(3, CLR_CYAN,
        "    [fetchNodes] sizes: CLUSTERINFO_CEMC=" << countClusters(m_clus)
        << " | PHOTONCLUSTER_CEMC=" << countClusters(m_photons));
  }

  if (!m_photons)
  {
    LOG(0, CLR_YELLOW,
        "    [fetchNodes] PHOTONCLUSTER_CEMC is MISSING → ABORTEVENT. "
        "PhotonClusterBuilder likely did not run or node name mismatch.");
    return false;
  }

  // ------------------------------------------------------------------
  // Reco jets: cache all configured radii in parallel.
  //
  // Radii selection is driven by:
  //   - setActiveJetRKeys([...]) (typically from YAML in Fun4All_recoilJets.C)
  //   - fallback: RecoilJets::kJetRadii (baseline r02+r04)
  //
  // Legacy knob (still honored as a "primary" label in logs only):
  //   export RJ_RECO_JET_KEY=r02   (or r04, etc.)
  // ------------------------------------------------------------------
  std::vector<std::string> rKeysToRun;
  if (!m_activeJetRKeys.empty())
      rKeysToRun = m_activeJetRKeys;
  else
  {
      rKeysToRun.reserve(kJetRadii.size());
      for (const auto& jnm : kJetRadii) rKeysToRun.push_back(jnm.key);
  }

  // de-dup defensively (preserve order)
  {
      std::vector<std::string> uniq;
      uniq.reserve(rKeysToRun.size());
      for (const auto& k : rKeysToRun)
      {
        if (k.empty()) continue;
        if (std::find(uniq.begin(), uniq.end(), k) == uniq.end())
          uniq.push_back(k);
      }
      rKeysToRun.swap(uniq);
    }

    if (rKeysToRun.empty())
      rKeysToRun = {kJetRadii.front().key};

    auto recoJetNodeForRKey = [&](const std::string& rKey) -> std::string
    {
      for (const auto& jnm : kJetRadii)
      {
        if (jnm.key == rKey)
          return (isAuAu ? jnm.aa_node : jnm.pp_node);
      }
      return std::string("AntiKt_Tower_") + rKey;
    };

    std::string primaryRecoKey = trim(m_xjRecoJetKey);
    if (const char* rk = std::getenv("RJ_RECO_JET_KEY"))
      primaryRecoKey = trim(std::string(rk));

    if (primaryRecoKey.empty())
      primaryRecoKey = rKeysToRun.front();

    bool keyOK = false;
    for (const auto& k : rKeysToRun)
      if (primaryRecoKey == k) { keyOK = true; break; }

    if (!keyOK)
    {
      primaryRecoKey = rKeysToRun.front();
      LOG(1, CLR_YELLOW,
          "    [fetchNodes] requested RJ_RECO_JET_KEY not in active radii → using \"" << primaryRecoKey << "\"");
    }

    // Persist (legacy) — may be printed elsewhere
    m_xjRecoJetKey = primaryRecoKey;

    m_jets.clear();
    m_jetsRaw.clear();
    for (const auto& rKey : rKeysToRun)
    {
      const std::string node = recoJetNodeForRKey(rKey);
      auto* jc = findNode::getClass<JetContainer>(top, node);

      if (!jc)
      {
        LOG(2, CLR_YELLOW, "  – reco jet node missing: " << node << " (rKey=" << rKey << ")");
      }
      else
      {
        LOG(4, CLR_GREEN,
            "    [fetchNodes] reco jet node found: " << node
            << "  rKey=" << rKey << "  jets=" << jc->size());
      }

      // Store even if nullptr; downstream QA/matching will skip gracefully per-radius
      m_jets[rKey] = jc;

      // If JetCalib is enabled, jets are built to "<node>_RAW" then calibrated into "<node>".
      // RAW jets reliably retain comp_vec (tower constituents) needed for EventDisplay diagnostics.
      const std::string nodeRaw = node + "_RAW";
      auto* jcRaw = findNode::getClass<JetContainer>(top, nodeRaw);
      m_jetsRaw[rKey] = jcRaw;
    }

    if (Verbosity() >= 4)
    {
      std::ostringstream os;
      os << "    [fetchNodes] reco jet radii cached: {";
      for (std::size_t i = 0; i < rKeysToRun.size(); ++i)
      {
        if (i) os << ", ";
        const std::string rk = rKeysToRun[i];
        auto it = m_jets.find(rk);
        const bool ok = (it != m_jets.end() && it->second);
        os << rk << ":" << (ok ? "OK" : "MISSING");
      }
      os << "}  (primary=" << primaryRecoKey << ")";
      LOG(4, CLR_BLUE, os.str());
    }

    // ------------------------------------------------------------------
    // Truth jets: SIM only.
    // Cache truth jet containers for the same set of radii.
    //
    // Recommended override for parallel radii:
    //   export RJ_TRUTH_JETS_NODE=AntiKt_Truth_{rKey}
    //   export RJ_TRUTH_JETS_NODE=AntiKt_TruthFromParticles_{rKey}
    //
    // If RJ_TRUTH_JETS_NODE is set without "{rKey}", we also try "<env>_<rKey>"
    // (plus the standard canonical nodes) for each radius.
    // ------------------------------------------------------------------
    m_truthJetsByRKey.clear();
    m_truthJetsNodeByRKey.clear();

    if (isSim)
    {
      const char* tn_env_c = std::getenv("RJ_TRUTH_JETS_NODE");
      const std::string tn_env = tn_env_c ? trim(std::string(tn_env_c)) : std::string{};

      for (const auto& rKey : rKeysToRun)
      {
        const std::string canonicalTruth = std::string("AntiKt_Truth_") + rKey;
        const std::string altTruth       = std::string("AntiKt_TruthFromParticles_") + rKey;

        std::vector<std::string> candidates;
        candidates.reserve(6);

        if (!tn_env.empty())
        {
          const std::string placeholder = "{rKey}";
          if (tn_env.find(placeholder) != std::string::npos)
          {
            std::string expanded = tn_env;
            expanded.replace(expanded.find(placeholder), placeholder.size(), rKey);
            candidates.push_back(expanded);
          }
          else
          {
            candidates.push_back(tn_env);
            if (tn_env.find("_" + rKey) == std::string::npos)
              candidates.push_back(tn_env + "_" + rKey);
          }
        }

        // Always try standard truth nodes for this radius
        candidates.push_back(canonicalTruth);
        candidates.push_back(altTruth);

        JetContainer* tj = nullptr;
        std::string   usedNode;

        for (const auto& node : candidates)
        {
          if (node.empty()) continue;
          if (auto* tmp = findNode::getClass<JetContainer>(top, node))
          {
            tj = tmp;
            usedNode = node;
            break;
          }
        }

        m_truthJetsByRKey[rKey] = tj;
        m_truthJetsNodeByRKey[rKey] = usedNode;

        if (!tj)
        {
          LOG(2, CLR_YELLOW,
              "    [fetchNodes] isSim: truth jet container NOT found for rKey=" << rKey
              << " (tried canonical=" << canonicalTruth
              << ", alt=" << altTruth
              << (tn_env.empty() ? "" : ", plus RJ_TRUTH_JETS_NODE candidates")
              << "). Truth-jet QA disabled for this radius.");
        }
        else
        {
          LOG(4, CLR_GREEN,
              "    [fetchNodes] isSim: truth jet node found for rKey=" << rKey
              << ": " << usedNode << "  jets=" << tj->size());
        }
      }
    }

    // -------------------------------------------------------------------------
    // EventDisplay diagnostics payload nodes (SIM only). Optional: never affects physics.
    //
    // IMPORTANT SAFETY:
    //   - This block must NEVER cause the event to be skipped.
    //   - If anything is missing, we simply mark diagnostics as unavailable for this event.
    // -------------------------------------------------------------------------
    m_evtDiagNodesReady = false;

    // Always clear per-event pointers
    m_evtHeader          = nullptr;
    m_evtDispTowersCEMC  = nullptr;
    m_evtDispTowersIHCal = nullptr;
    m_evtDispTowersOHCal = nullptr;
    m_evtDispGeomCEMC    = nullptr;
    m_evtDispGeomIHCal   = nullptr;
    m_evtDispGeomOHCal   = nullptr;

    if (!(isSim && m_evtDiagEnabled))
    {
      if (Verbosity() >= 5)
      {
        LOG(5, CLR_CYAN,
            "[EventDisplayTree][nodes] skipped: "
            << "isSim=" << (isSim ? "true" : "false")
            << "  m_evtDiagEnabled=" << (m_evtDiagEnabled ? "true" : "false"));
      }
    }
    else
    {
      // NOTE: do NOT assume any node exists; only set m_evtDiagNodesReady=true if ALL required are present.
      m_evtHeader = findNode::getClass<EventHeader>(top, "EventHeader");

      // JetReco in Fun4All_recoilJets.C builds jets from TowerJetInput(..., "TOWERINFO_CALIB"),
      // so comp.second indices are defined against TOWERINFO_CALIB. Prefer that container.
      TowerInfoContainer* towersAll = findNode::getClass<TowerInfoContainer>(top, "TOWERINFO_CALIB");

      const bool usingAll = (towersAll != nullptr);
      const char* nodeTowersAll = "TOWERINFO_CALIB";
      const char* nodeTowersCEMC = usingAll ? "TOWERINFO_CALIB" : "TOWERINFO_CALIB_CEMC";
      const char* nodeTowersIHCal = usingAll ? "TOWERINFO_CALIB" : "TOWERINFO_CALIB_HCALIN";
      const char* nodeTowersOHCal = usingAll ? "TOWERINFO_CALIB" : "TOWERINFO_CALIB_HCALOUT";

      if (towersAll)
      {
        m_evtDispTowersCEMC  = towersAll;
        m_evtDispTowersIHCal = towersAll;
        m_evtDispTowersOHCal = towersAll;
      }
      else
      {
        m_evtDispTowersCEMC  = findNode::getClass<TowerInfoContainer>(top, "TOWERINFO_CALIB_CEMC");
        m_evtDispTowersIHCal = findNode::getClass<TowerInfoContainer>(top, "TOWERINFO_CALIB_HCALIN");
        m_evtDispTowersOHCal = findNode::getClass<TowerInfoContainer>(top, "TOWERINFO_CALIB_HCALOUT");
      }

      m_evtDispGeomCEMC    = findNode::getClass<RawTowerGeomContainer>(top, "TOWERGEOM_CEMC");
      m_evtDispGeomIHCal   = findNode::getClass<RawTowerGeomContainer>(top, "TOWERGEOM_HCALIN");
      m_evtDispGeomOHCal   = findNode::getClass<RawTowerGeomContainer>(top, "TOWERGEOM_HCALOUT");

      const bool missing =
          (!m_evtHeader ||
           !m_evtDispTowersCEMC || !m_evtDispTowersIHCal || !m_evtDispTowersOHCal ||
           !m_evtDispGeomCEMC   || !m_evtDispGeomIHCal   || !m_evtDispGeomOHCal);

      if (!missing)
      {
        m_evtDiagNodesReady = true;

        if (Verbosity() >= 4)
        {
          LOG(4, CLR_GREEN,
              "[EventDisplayTree][nodes][OK] "
              << "EventHeader=OK"
              << "  towers(" << (usingAll ? nodeTowersAll : "split") << "): "
              << "CEMC=" << (m_evtDispTowersCEMC ? "OK" : "MISSING")
              << " IHCal=" << (m_evtDispTowersIHCal ? "OK" : "MISSING")
              << " OHCal=" << (m_evtDispTowersOHCal ? "OK" : "MISSING")
              << "  geom: "
              << "CEMC=" << (m_evtDispGeomCEMC ? "OK" : "MISSING")
              << " IHCal=" << (m_evtDispGeomIHCal ? "OK" : "MISSING")
              << " OHCal=" << (m_evtDispGeomOHCal ? "OK" : "MISSING"));
        }

        if (Verbosity() >= 6)
        {
          const auto nC = (m_evtDispTowersCEMC  ? m_evtDispTowersCEMC->size()  : 0);
          const auto nI = (m_evtDispTowersIHCal ? m_evtDispTowersIHCal->size() : 0);
          const auto nO = (m_evtDispTowersOHCal ? m_evtDispTowersOHCal->size() : 0);
          LOG(6, CLR_CYAN,
              "[EventDisplayTree][nodes][sizes] "
              << "CEMC=" << nC << "  HCALIN=" << nI << "  HCALOUT=" << nO
              << "  (usingAll=" << (usingAll ? "YES" : "NO") << ")");
        }
      }
      else
      {
        static bool s_warned_once = false;
        static int  s_warned_count = 0;

        if (!s_warned_once)
        {
          LOG(1, CLR_YELLOW, "[EventDisplayTree] disabled for this event: missing required node(s) for diagnostics payload");
          LOG(1, CLR_YELLOW, "  EventHeader            : " << (m_evtHeader ? "OK" : "MISSING"));
          LOG(1, CLR_YELLOW, "  " << nodeTowersCEMC     << " : " << (m_evtDispTowersCEMC ? "OK" : "MISSING"));
          LOG(1, CLR_YELLOW, "  " << nodeTowersIHCal    << " : " << (m_evtDispTowersIHCal ? "OK" : "MISSING"));
          LOG(1, CLR_YELLOW, "  " << nodeTowersOHCal    << " : " << (m_evtDispTowersOHCal ? "OK" : "MISSING"));
          LOG(1, CLR_YELLOW, "  TOWERGEOM_CEMC         : " << (m_evtDispGeomCEMC ? "OK" : "MISSING"));
          LOG(1, CLR_YELLOW, "  TOWERGEOM_HCALIN       : " << (m_evtDispGeomIHCal ? "OK" : "MISSING"));
          LOG(1, CLR_YELLOW, "  TOWERGEOM_HCALOUT      : " << (m_evtDispGeomOHCal ? "OK" : "MISSING"));
          LOG(1, CLR_YELLOW, "  towersAll(" << nodeTowersAll << "): " << (towersAll ? "FOUND" : "MISSING"));
          s_warned_once = true;
        }

        if (Verbosity() >= 4 && s_warned_count < 5)
        {
          LOG(4, CLR_YELLOW,
              "[EventDisplayTree][nodes][MISSING] "
              << "EventHeader=" << (m_evtHeader ? "OK" : "MISSING")
              << "  Towers: "
              << nodeTowersCEMC << "=" << (m_evtDispTowersCEMC ? "OK" : "MISSING")
              << " " << nodeTowersIHCal << "=" << (m_evtDispTowersIHCal ? "OK" : "MISSING")
              << " " << nodeTowersOHCal << "=" << (m_evtDispTowersOHCal ? "OK" : "MISSING")
              << "  Geom: "
              << "TOWERGEOM_CEMC=" << (m_evtDispGeomCEMC ? "OK" : "MISSING")
              << " TOWERGEOM_HCALIN=" << (m_evtDispGeomIHCal ? "OK" : "MISSING")
              << " TOWERGEOM_HCALOUT=" << (m_evtDispGeomOHCal ? "OK" : "MISSING"));
          ++s_warned_count;
        }

        m_evtDiagNodesReady = false;
      }
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

  // ------------------------------------------------------------------
  // reset per-pT negative-isolation accounting (SS-independent)
  // ------------------------------------------------------------------
  const std::size_t nPtBins = (m_gammaPtBins.size() > 1 ? (m_gammaPtBins.size() - 1) : 0);
  m_nIsoBuilderByPt.assign(nPtBins, 0ULL);
  m_nIsoBuilderNegByPt.assign(nPtBins, 0ULL);

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

    // -------------------------------------------------------------------------
    // EventDisplay diagnostics payload (offline rendering; independent of Verbosity()).
    //
    // Enable via either:
    //   - RecoilJets::enableEventDisplayDiagnostics(true), or
    //   - environment: RJ_ENABLE_EVENTDISPLAY_DIAGNOSTICS=1
    //
    // Optional cap (per rKey × pTγ bin × {NUM,MissA,MissB}):
    //   - environment: RJ_EVTDIAG_MAX_PER_BIN=N   (0 => unlimited)
    // -------------------------------------------------------------------------
    if (const char* e = gSystem->Getenv("RJ_ENABLE_EVENTDISPLAY_DIAGNOSTICS"))
    {
      m_evtDiagEnabled = (std::atoi(e) != 0);
    }

    if (const char* n = gSystem->Getenv("RJ_EVTDIAG_MAX_PER_BIN"))
    {
      m_evtDiagMaxPerBin = std::max(0, std::atoi(n));
    }

    m_evtDiagSavedPerBin.clear();

    if (m_evtDiagEnabled)
    {
      initEventDisplayDiagnosticsTree();
    }


    LOG(1, CLR_BLUE, "[InitRun] InitRun completed successfully");
    return Fun4AllReturnCodes::EVENT_OK;
  }



void RecoilJets::createHistos_Data()
{
  // Vertex histogram binning tied to your configured vz cut (m_vzCut).
  // Keep the same resolution as before: ~0.5 cm per bin.
  const double vzAbs    = std::max(1.0, std::fabs(static_cast<double>(m_vzCut)));
  const double vzBinW   = 0.5; // cm/bin
  const int    nbVz     = std::max(1, static_cast<int>(std::lround((2.0 * vzAbs) / vzBinW)));
  const double vzMin    = -vzAbs;
  const double vzMax    =  vzAbs;

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

      // 3) SIM ONLY: truth-vs-reco vertex comparison (filled BEFORE |vz| cut)
      if (H.find("h_vzTruthVsReco") == H.end())
      {
        auto* h2 = new TH2F("h_vzTruthVsReco",
                            "h_vzTruthVsReco;v_{z}^{truth} [cm];v_{z}^{reco-used} [cm]",
                            nbVz, vzMin, vzMax,
                            nbVz, vzMin, vzMax);
        H["h_vzTruthVsReco"] = h2;
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
  // DATA MODE: per-trigger directories
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

  // Small helper for printing lists
  auto joinList = [](const std::vector<std::string>& v, const char* sep) -> std::string
  {
    std::ostringstream os;
    for (std::size_t i = 0; i < v.size(); ++i)
    {
      if (i) os << sep;
      os << v[i];
    }
    return os.str();
  };

  // Grab GL1 scaled vector for diagnostics (works for pp & Au+Au if node exists)
  uint64_t gl1Scaled = 0;
  {
    Gl1Packet* gl1 = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    if (!gl1) gl1 = findNode::getClass<Gl1Packet>(topNode, "14001");
    if (gl1) gl1Scaled = gl1->lValue(0, "ScaledVector");
  }

  // ------------------------------------------------------------------
  // SIMULATION MODE:
  //   - skip trigger logic
  //   - still apply global |vz| cut if enabled
  // ------------------------------------------------------------------
  if (m_isSim)
  {
    activeTrig.emplace_back("SIM");

    if (m_useVzCut && std::fabs(m_vz) >= m_vzCut)
    {
      m_lastReject = EventReject::Vz;
      if (Verbosity() >= 4)
      {
        std::ostringstream os;
        os << "    [firstEventCuts] REJECT (SIM |vz| cut)"
           << " | vz=" << std::fixed << std::setprecision(3) << m_vz
           << " | |vz|cut=" << m_vzCut;
        LOG(4, CLR_YELLOW, os.str());
      }
      return false;
    }

    if (Verbosity() >= 4)
    {
      std::ostringstream os;
      os << "    [firstEventCuts] ACCEPT (SIM)"
         << " | vz=" << std::fixed << std::setprecision(3) << m_vz;
      if (m_useVzCut) os << " (|vz|cut=" << m_vzCut << ")";
      LOG(4, CLR_GREEN, os.str());
    }

    m_lastReject = EventReject::None;
    return true;
  }

  // ------------------------------------------------------------------
  // DATA MODE
  // ------------------------------------------------------------------
  if (!m_isAuAu)
  {
    // ---------- p+p: accept if ANY configured pp trigger fired ----------
    if (!trigAna)
    {
      m_lastReject = EventReject::Trigger;
      LOG(3, CLR_YELLOW, "    [firstEventCuts] REJECT (p+p): trigAna == nullptr (TriggerAnalyzer missing)");
      return false;
    }

    trigAna->decodeTriggers(topNode);

    std::vector<std::string> firedDb;
    std::vector<std::string> checkedLines;
    firedDb.reserve(triggerNameMap_pp.size());
    checkedLines.reserve(triggerNameMap_pp.size());

    for (const auto& kv : triggerNameMap_pp)
    {
      const std::string& dbName   = kv.first;   // GL1/DB name
      const std::string& shortKey = kv.second;  // friendly key

      const bool fired = trigAna->didTriggerFire(dbName);
      checkedLines.push_back(dbName + "->" + (fired ? "1" : "0") + " (" + shortKey + ")");
      if (fired) activeTrig.push_back(shortKey);
      if (fired) firedDb.push_back(dbName);

      if (Verbosity() >= 6)
      {
        LOG(6, CLR_CYAN,
            "      [pp trigger] " << dbName
            << " fired=" << (fired ? 1 : 0)
            << "  shortKey=\"" << shortKey << "\"");
      }
    }

    if (activeTrig.empty())
    {
      m_lastReject = EventReject::Trigger;

      if (Verbosity() >= 4)
      {
        std::ostringstream os;
        os << "    [firstEventCuts] REJECT (p+p Trigger): no configured pp triggers fired"
           << " | checked={" << joinList(checkedLines, ", ") << "}"
           << " | vz=" << std::fixed << std::setprecision(3) << m_vz;
        if (m_useVzCut) os << " (|vz|cut=" << m_vzCut << ")";
        if (gl1Scaled) os << " | GL1Scaled=0x" << std::hex << gl1Scaled << std::dec;
        LOG(4, CLR_YELLOW, os.str());
      }
      return false;
    }

    if (Verbosity() >= 4)
    {
      std::ostringstream os;
      os << "    [firstEventCuts] Trigger PASS (p+p)"
         << " | activeTrig={" << joinList(activeTrig, ", ") << "}"
         << " | vz=" << std::fixed << std::setprecision(3) << m_vz;
      if (gl1Scaled) os << " | GL1Scaled=0x" << std::hex << gl1Scaled << std::dec;
      LOG(4, CLR_GREEN, os.str());
    }
  }
  else
  {
    // ---------- Au+Au: scaled GL1 bits only ----------
    uint64_t wScaled = 0;
    if (auto* gl1 = findNode::getClass<Gl1Packet>(topNode, "GL1Packet"))
      wScaled = gl1->lValue(0, "ScaledVector");
    else if (auto* gl1b = findNode::getClass<Gl1Packet>(topNode, "14001"))
      wScaled = gl1b->lValue(0, "ScaledVector");

    if (Verbosity() >= 6)
    {
      std::ostringstream os;
      os << "    [AuAu trigger] ScaledVector="
         << "0x" << std::hex << wScaled << std::dec
         << " (event_count=" << event_count << ")";
      LOG(6, CLR_CYAN, os.str());
    }

    if (wScaled)
    {
      const auto bits = extractTriggerBits(wScaled, static_cast<int>(event_count));
      for (const auto& [bit, key] : triggerNameMapAuAu)
        if (checkTriggerCondition(bits, bit))
          activeTrig.push_back(key);
    }

    if (activeTrig.empty())
    {
      m_lastReject = EventReject::Trigger;

      if (Verbosity() >= 4)
      {
        std::ostringstream os;
        os << "    [firstEventCuts] REJECT (Au+Au Trigger): no configured scaled triggers fired"
           << " | ScaledVector=0x" << std::hex << wScaled << std::dec
           << " | vz=" << std::fixed << std::setprecision(3) << m_vz;
        if (m_useVzCut) os << " (|vz|cut=" << m_vzCut << ")";
        LOG(4, CLR_YELLOW, os.str());
      }
      return false;
    }

    if (Verbosity() >= 4)
    {
      std::ostringstream os;
      os << "    [firstEventCuts] Trigger PASS (Au+Au)"
         << " | activeTrig={" << joinList(activeTrig, ", ") << "}"
         << " | vz=" << std::fixed << std::setprecision(3) << m_vz;
      LOG(4, CLR_GREEN, os.str());
    }
  }

  // ---------- Global |vz| veto ----------
  if (m_useVzCut && std::fabs(m_vz) >= m_vzCut)
  {
    m_lastReject = EventReject::Vz;

    if (Verbosity() >= 4)
    {
      std::ostringstream os;
      os << "    [firstEventCuts] REJECT (|vz| cut)"
         << " | vz=" << std::fixed << std::setprecision(3) << m_vz
         << " | |vz|cut=" << m_vzCut
         << " | triggers={" << joinList(activeTrig, ", ") << "}";
      LOG(4, CLR_YELLOW, os.str());
    }
    return false;
  }

  // Uniform downstream behavior
  if (activeTrig.empty()) activeTrig.emplace_back("ALL");

  m_lastReject = EventReject::None;

  if (Verbosity() >= 4)
  {
    std::ostringstream os;
    os << "    [firstEventCuts] ACCEPT"
       << " | triggers={" << joinList(activeTrig, ", ") << "}"
       << " | vz=" << std::fixed << std::setprecision(3) << m_vz;
    if (m_useVzCut) os << " (|vz|cut=" << m_vzCut << ")";
    LOG(4, CLR_GREEN, os.str());
  }

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

      const char* why = "UNKNOWN";
      if (m_lastReject == EventReject::Trigger) why = "Trigger";
      else if (m_lastReject == EventReject::Vz) why = "|vz|";

      std::ostringstream os;
      os << "    event rejected by " << why << " gate – skip"
         << " | vz=" << std::fixed << std::setprecision(3) << m_vz;
      if (m_useVzCut) os << " (|vz|cut=" << m_vzCut << ")";
      os << " | nActiveTrig=" << activeTrig.size();

      if (!activeTrig.empty())
      {
        os << " | triggers={";
        for (std::size_t i = 0; i < activeTrig.size(); ++i)
        {
          if (i) os << ", ";
          os << activeTrig[i];
        }
        os << "}";
      }

      LOG(4, CLR_YELLOW, os.str());
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  /* ------------------------------------------------------------------ */
  /* 3) Trigger counters (one per trigger) + Vertex-z QA                */
  /* ------------------------------------------------------------------ */

    // Bump the per-trigger counter once per accepted event (and LOG fills via bumpHistFill)
    for (const auto& t : activeTrig)
    {
      auto itTrig = qaHistogramsByTrigger.find(t);
      if (itTrig == qaHistogramsByTrigger.end()) continue;

      auto& H = itTrig->second;

      // (1) per-trigger event counter
      const std::string hcnt = "cnt_" + t;
      if (auto hc = H.find(hcnt); hc != H.end())
      {
        static_cast<TH1I*>(hc->second)->Fill(1);
        bumpHistFill(t, hcnt);
      }

      // (2) vertex-z QA
      if (auto hvz = H.find("h_vertexZ"); hvz != H.end())
      {
        static_cast<TH1F*>(hvz->second)->Fill(m_vz);
        bumpHistFill(t, "h_vertexZ");
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
          if (itTrig == qaHistogramsByTrigger.end()) continue;

          auto& H = itTrig->second;
          if (auto hc = H.find("h_centrality"); hc != H.end())
          {
            static_cast<TH1F*>(hc->second)->Fill(centile);
            bumpHistFill(t, "h_centrality");
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
  if (!std::isfinite(m_vz) || (m_useVzCut && std::fabs(m_vz) >= m_vzCut))
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
  for (const auto& kv : m_jets)
  {
        fillInclusiveJetQA(activeTrig, centIdxForJets, kv.first);
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
//  End – Diagnostics, robust against dangling pointers
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
  if (!out)
  {
    warn("TFile* 'out' is nullptr – nothing to write");
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if (!out->IsOpen())
  {
    warn("Output file is *not* open (" + std::string(out->GetName()) + ")");
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  //--------------------------------------------------------------------
  // 2. Write histograms trigger-by-trigger
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
      ++nHistExpected;
      ++nExpThis;

      TH1* h = dynamic_cast<TH1*>(obj);
      if (!h)
      {
        warn("Object '" + key + "' (trigger " + trig + ") is not TH1 – skipped");
        continue;
      }

      if (h->GetEntries() == 0)
      {
        if (Verbosity() > 1)
          warn("Histogram '" + key + "' (trigger " + trig + ") has 0 entries – skipped");
        continue;
      }

      try
      {
        if (h->Write("", TObject::kOverwrite) > 0)
        {
          ++nHistWritten;
          ++nWrtThis;
        }
        else
        {
          warn("Write() returned 0 for '" + key + "' (trigger " + trig + ")");
        }
      }
      catch (const std::exception& e)
      {
        warn("Exception while writing '" + key + "' – " + std::string(e.what()));
      }
    }

    info(1, "trigger '" + trig + "': " + std::to_string(nWrtThis) + " / "
               + std::to_string(nExpThis) + " histograms written");
    out->cd();
  }

  //--------------------------------------------------------------------
  // 3. Human-readable histogram summary (must run *before* file deletion)
  //--------------------------------------------------------------------
  if (Verbosity() > 0)
  {
    std::cout << "\n\033[1mHistogram summary\033[0m\n"
              << "\033[1mTrigger                        │ Histogram                           │  Entries\033[0m\n"
              << "-------------------------------------------------------------------------------\n";

    for (const auto& [trig, hMap] : qaHistogramsByTrigger)
    {
      for (const auto& [key, obj] : hMap)
      {
        const TH1* h = dynamic_cast<const TH1*>(obj);
        if (!h) continue;
        if (h->GetEntries() == 0) continue;  // only list filled ones

        std::cout << std::left  << std::setw(30) << trig << " │ "
                  << std::setw(32) << key  << " │ "
                  << std::right << std::setw(10)
                  << static_cast<Long64_t>(h->GetEntries()) << '\n';
      }
    }
  }

    // --------------------------------------------------------------------
    // 4. Analysis summary (verbosity-controlled)
    // --------------------------------------------------------------------
    if (Verbosity() >= 1)
    {
      const char* dsLabel = (m_isSim ? "SIM" : (m_isAuAu ? "Au+Au" : "p+p"));

      std::cout << "\n\033[1mSelection summary (dataset: " << dsLabel
                << ", events=" << event_count << ")\033[0m\n";

      // -------------------- ABCD summary (existing) --------------------
      for (const auto& kvT : m_catByTrig)
      {
        const std::string& trig = kvT.first;

        std::cout << "\n\033[1mTrigger: " << trig << "\033[0m\n";

        std::cout << "ABCD (PPG12): iso=Eiso<thrIso; nSB=Eiso>thrIso+gap; GAP excluded; "
                     "NT2=non-tight(>=2 fails); Neither(1 fail) excluded.\n";

        std::cout << "slice                                 |  A(iso,T) | B(nSB,T)   | C(iso,NT2) | D(nSB,NT2)  |  seen\n";
        std::cout << "--------------------------------------+-----------+-----------+------------+-------------+-------\n";

        std::vector<std::string> keys;
        keys.reserve(kvT.second.size());
        for (const auto& s : kvT.second) keys.push_back(s.first);
        std::sort(keys.begin(), keys.end());

        for (const auto& sfx : keys)
        {
          const CatStat& S = kvT.second.at(sfx);
          std::ostringstream lab;
          lab << (sfx.empty() ? "<all>" : sfx);

          std::cout << std::left  << std::setw(38) << lab.str() << " | "
                    << std::right << std::setw(9)  << S.n_iso_tight       << " | "
                    << std::setw(9)  << S.n_nonIso_tight    << " | "
                    << std::setw(10) << S.n_iso_nonTight    << " | "
                    << std::setw(11) << S.n_nonIso_nonTight << " | "
                    << std::setw(5)  << S.seen              << "\n";
        }
      }

      // -------------------- Photon cutflow (existing) --------------------
      printCutSummary();


      // =========================================================================
      //  Jet matching + unfolding summary (computed from already-filled hists)
      //
      // Notes:
      //  - "match status" is filled ONCE per event that has a leading iso∧tight γ
      //    (one fill per rKey), using:
      //      1=NoJetPt, 2=NoJetEta, 3=NoBackToBack, 4=Matched
      //  - Unfolding TH2 "reco" is filled once per (leading iso∧tight γ) × (each recoil jet).
      // =========================================================================
      auto starts_with = [](const std::string& s, const std::string& pfx) -> bool
      {
        return (s.size() >= pfx.size() && std::equal(pfx.begin(), pfx.end(), s.begin()));
      };

      auto fmt_frac = [](double num, double den) -> std::string
      {
        std::ostringstream os;
        if (den <= 0.0) { os << "n/a"; return os.str(); }
        os << std::fixed << std::setprecision(3) << (num / den);
        return os.str();
      };

      struct TH2Stats { double entries = 0.0; double inrange = 0.0; double meanY = std::numeric_limits<double>::quiet_NaN(); };

      auto th2_y_stats = [&](const TH2* h) -> TH2Stats
      {
        TH2Stats S;
        if (!h) return S;

        S.entries = h->GetEntries();

        const int nx = h->GetNbinsX();
        const int ny = h->GetNbinsY();

        double sumW = 0.0;
        double sumY = 0.0;

        for (int ix = 1; ix <= nx; ++ix)
        {
          for (int iy = 1; iy <= ny; ++iy)
          {
            const double w = h->GetBinContent(ix, iy);
            if (w == 0.0) continue;
            const double y = h->GetYaxis()->GetBinCenter(iy);
            sumW += w;
            sumY += w * y;
          }
        }

        S.inrange = sumW;
        if (sumW > 0.0) S.meanY = sumY / sumW;
        return S;
      };

      auto sum_th1_by_prefix = [&](const HistMap& H, const std::string& pfx,
                                   double& outEntries, double& outMean) -> void
      {
        double n = 0.0;
        double sumMean = 0.0;

        for (const auto& kv : H)
        {
          const std::string& name = kv.first;
          if (!starts_with(name, pfx)) continue;
          const TH1* h = dynamic_cast<const TH1*>(kv.second);
          if (!h) continue;

          const double e = static_cast<double>(h->GetEntries());
          if (e <= 0.0) continue;

          n += e;
          sumMean += h->GetMean() * e;
        }

        outEntries = n;
        outMean    = (n > 0.0 ? (sumMean / n) : std::numeric_limits<double>::quiet_NaN());
      };

      std::cout << "\n\033[1mJet matching summary (leading iso\u2227tight photons)\033[0m\n";
      std::cout << "Jet cuts: pT > " << std::fixed << std::setprecision(2) << m_minJetPt
                << " GeV, |Δφ(γ,jet)| >= " << std::setprecision(3) << m_minBackToBack
                << ", fiducial |η_jet| < 1.1 - R\n";

      for (const auto& trigPair : qaHistogramsByTrigger)
      {
        const std::string& trig = trigPair.first;
        const HistMap&      H   = trigPair.second;

        std::cout << "\n\033[1mTrigger: " << trig << "\033[0m\n";

        std::cout << "rKey  | |eta|< (1.1-R) | Nlead(iso\u2227tight) | NoJetPt | NoJetEta | NoBackToBack | Matched | fMatched | fBackToBack(given jets)\n";
        std::cout << "------+---------------+------------------+---------+----------+-------------+---------+---------+--------------------------\n";

        const bool useActiveRKeys = (!m_activeJetRKeys.empty());
        const std::size_t nRKeys = useActiveRKeys ? m_activeJetRKeys.size() : kJetRadii.size();

        for (std::size_t iR = 0; iR < nRKeys; ++iR)
        {
            const std::string rKey = useActiveRKeys ? m_activeJetRKeys[iR] : kJetRadii[iR].key;
            const double etaMax = jetEtaAbsMaxForRKey(rKey);

          // Collect ALL match-status histograms for this rKey (includes centrality-sliced ones in Au+Au)
          const std::string pfxStatus = std::string("h_match_status_vs_pTgamma_") + rKey;

          std::vector<const TH2*> statusHists;
          for (const auto& kv : H)
          {
            if (!starts_with(kv.first, pfxStatus)) continue;
            if (auto* h = dynamic_cast<const TH2*>(kv.second)) statusHists.push_back(h);
          }

          double c1 = 0.0, c2 = 0.0, c3 = 0.0, c4 = 0.0;
          std::vector<std::array<double, 5>> perX; // [ix][1..4]

          if (!statusHists.empty())
          {
            const int nx0 = statusHists.front()->GetNbinsX();
            perX.assign(static_cast<std::size_t>(nx0 + 1), std::array<double, 5>{0,0,0,0,0});

            for (const TH2* h : statusHists)
            {
              const int nx = h->GetNbinsX();
              const int ny = h->GetNbinsY();
              if (nx != nx0 || ny < 4) continue;

              for (int ix = 1; ix <= nx; ++ix)
              {
                const double v1 = h->GetBinContent(ix, 1);
                const double v2 = h->GetBinContent(ix, 2);
                const double v3 = h->GetBinContent(ix, 3);
                const double v4 = h->GetBinContent(ix, 4);

                c1 += v1; c2 += v2; c3 += v3; c4 += v4;

                perX[ix][1] += v1;
                perX[ix][2] += v2;
                perX[ix][3] += v3;
                perX[ix][4] += v4;
              }
            }
          }

          const double nLead = c1 + c2 + c3 + c4;
          const double nJetsExist = c3 + c4;

          std::cout << std::left << std::setw(4) << rKey << " | "
                    << std::right << std::setw(13) << std::fixed << std::setprecision(2) << etaMax << " | "
                    << std::setw(16) << static_cast<long long>(std::llround(nLead)) << " | "
                    << std::setw(7)  << static_cast<long long>(std::llround(c1)) << " | "
                    << std::setw(8)  << static_cast<long long>(std::llround(c2)) << " | "
                    << std::setw(11) << static_cast<long long>(std::llround(c3)) << " | "
                    << std::setw(7)  << static_cast<long long>(std::llround(c4)) << " | "
                    << std::setw(7)  << fmt_frac(c4, nLead) << " | "
                    << std::setw(24) << fmt_frac(c4, nJetsExist)
                    << "\n";

          // Optional per-pTγ-bin breakdown (very useful for debugging)
          if (Verbosity() >= 2 && !statusHists.empty())
          {
            const TH2* href = statusHists.front();
            const int nx = href->GetNbinsX();

            std::cout << "    per-pTγ match status (" << rKey << "):\n";
            std::cout << "    pTγ bin     | NoJetPt | NoJetEta | NoBackToBack | Matched | N\n";
            std::cout << "    ------------+---------+----------+-------------+---------+------\n";

            for (int ix = 1; ix <= nx; ++ix)
            {
              const double lo = href->GetXaxis()->GetBinLowEdge(ix);
              const double hi = href->GetXaxis()->GetBinUpEdge(ix);

              const double b1 = perX[ix][1];
              const double b2 = perX[ix][2];
              const double b3 = perX[ix][3];
              const double b4 = perX[ix][4];
              const double bt = b1 + b2 + b3 + b4;

              if (bt <= 0.0) continue;

              std::ostringstream lab;
              lab << std::fixed << std::setprecision(0) << lo << "–" << hi;

              std::cout << "    " << std::left << std::setw(12) << lab.str() << " | "
                        << std::right << std::setw(7)  << static_cast<long long>(std::llround(b1)) << " | "
                        << std::setw(8)  << static_cast<long long>(std::llround(b2)) << " | "
                        << std::setw(11) << static_cast<long long>(std::llround(b3)) << " | "
                        << std::setw(7)  << static_cast<long long>(std::llround(b4)) << " | "
                        << std::setw(6)  << static_cast<long long>(std::llround(bt)) << "\n";
            }
          }

          // -------------------- matched-only / jet-physics summaries --------------------
          // Mean max|dphi| over events that had >=1 jet passing pT+eta (NoBackToBack + Matched)
          const std::string pfxMaxDphi = std::string("h_match_maxdphi_vs_pTgamma_") + rKey;
          TH2Stats Smax;
          {
            double sumEntries = 0.0;
            double sumInRange = 0.0;
            double sumMeanYw  = 0.0;

            int nUsed = 0;
            for (const auto& kv : H)
            {
              if (!starts_with(kv.first, pfxMaxDphi)) continue;
              const TH2* h = dynamic_cast<const TH2*>(kv.second);
              if (!h) continue;

              TH2Stats s = th2_y_stats(h);
              sumEntries += s.entries;
              sumInRange += s.inrange;
              if (std::isfinite(s.meanY) && s.inrange > 0.0)
              {
                sumMeanYw += s.meanY * s.inrange;
                ++nUsed;
              }
            }

            Smax.entries = sumEntries;
            Smax.inrange = sumInRange;
            if (sumInRange > 0.0) Smax.meanY = sumMeanYw / sumInRange;
          }

          // Mean |dphi| for matched recoil jet1 (only filled when Matched)
          const std::string pfxDphi = std::string("h_match_dphi_vs_pTgamma_") + rKey;
          TH2Stats Sdphi;
          {
            double sumEntries = 0.0;
            double sumInRange = 0.0;
            double sumMeanYw  = 0.0;

            int nUsed = 0;
            for (const auto& kv : H)
            {
              if (!starts_with(kv.first, pfxDphi)) continue;
              const TH2* h = dynamic_cast<const TH2*>(kv.second);
              if (!h) continue;

              TH2Stats s = th2_y_stats(h);
              sumEntries += s.entries;
              sumInRange += s.inrange;
              if (std::isfinite(s.meanY) && s.inrange > 0.0)
              {
                sumMeanYw += s.meanY * s.inrange;
                ++nUsed;
              }
            }

            Sdphi.entries = sumEntries;
            Sdphi.inrange = sumInRange;
            if (sumInRange > 0.0) Sdphi.meanY = sumMeanYw / sumInRange;
          }

          // Mean <N recoil jets> from the profile (if present)
          const std::string pfxNrecoil = std::string("p_nRecoilJets_vs_pTgamma_") + rKey;
          double profEntries = 0.0;
          double profSum     = 0.0;
          for (const auto& kv : H)
          {
            if (!starts_with(kv.first, pfxNrecoil)) continue;
            const TProfile* p = dynamic_cast<const TProfile*>(kv.second);
            if (!p) continue;

            const int nx = p->GetNbinsX();
            for (int ix = 1; ix <= nx; ++ix)
            {
              const double e = p->GetBinEntries(ix);
              if (e <= 0.0) continue;
              profEntries += e;
              profSum     += p->GetBinContent(ix) * e;
            }
          }
          const double meanNrecoil = (profEntries > 0.0 ? profSum / profEntries : std::numeric_limits<double>::quiet_NaN());

          // Matched recoil-jet outputs (leading recoil jet only)
          double xjN=0.0, xjMean=std::numeric_limits<double>::quiet_NaN();
          double aN =0.0, aMean =std::numeric_limits<double>::quiet_NaN();
          double j1N=0.0, j1Mean=std::numeric_limits<double>::quiet_NaN();

          sum_th1_by_prefix(H, std::string("h_xJ_")     + rKey, xjN, xjMean);
          sum_th1_by_prefix(H, std::string("h_alpha_")  + rKey, aN,  aMean);
          sum_th1_by_prefix(H, std::string("h_jet1Pt_") + rKey, j1N, j1Mean);

            auto sum_th2_prefix = [&](const std::string& pfx) -> TH2Stats
            {
              TH2Stats stats;
              double sumEntries = 0.0, sumInRange = 0.0, sumMeanYw = 0.0;

              for (const auto& kv : H)
              {
                if (!starts_with(kv.first, pfx)) continue;
                const TH2* h = dynamic_cast<const TH2*>(kv.second);
                if (!h) continue;

                TH2Stats s = th2_y_stats(h);
                sumEntries += s.entries;
                sumInRange += s.inrange;
                if (std::isfinite(s.meanY) && s.inrange > 0.0)
                  sumMeanYw += s.meanY * s.inrange;
              }

              stats.entries = sumEntries;
              stats.inrange = sumInRange;
              if (sumInRange > 0.0) stats.meanY = sumMeanYw / sumInRange;
              return stats;
            };

          const TH2Stats Ureco  = sum_th2_prefix(std::string("h2_unfoldReco_pTgamma_xJ_incl_")       + rKey);
          const TH2Stats Utruth = sum_th2_prefix(std::string("h2_unfoldTruth_pTgamma_xJ_incl_")      + rKey);
          const TH2Stats Ufakes = sum_th2_prefix(std::string("h2_unfoldRecoFakes_pTgamma_xJ_incl_")  + rKey);
          const TH2Stats Umiss  = sum_th2_prefix(std::string("h2_unfoldTruthMisses_pTgamma_xJ_incl_")+ rKey);

          // Response is "global bin vs global bin" => only count fills
          double rspEntries = 0.0;
          {
            const std::string pfxRsp = std::string("h2_unfoldResponse_pTgamma_xJ_incl_") + rKey;
            for (const auto& kv : H)
            {
              if (!starts_with(kv.first, pfxRsp)) continue;
              const TH2* h = dynamic_cast<const TH2*>(kv.second);
              if (!h) continue;
              rspEntries += h->GetEntries();
            }
          }

          // Print the “extra” stats only if something exists
          if (nLead > 0.0 || xjN > 0.0 || Ureco.entries > 0.0 || Utruth.entries > 0.0 || rspEntries > 0.0)
          {
            std::cout << "    " << rKey << " details:\n";
            if (std::isfinite(meanNrecoil))
              std::cout << "      <N_recoil jets> (profile): " << std::fixed << std::setprecision(3) << meanNrecoil
                        << " (entries=" << static_cast<long long>(std::llround(profEntries)) << ")\n";

            if (std::isfinite(Smax.meanY))
              std::cout << "      <max|Δφ|> over jets passing pT+eta: " << std::fixed << std::setprecision(3) << Smax.meanY
                        << " rad (in-range fills=" << static_cast<long long>(std::llround(Smax.inrange)) << ")\n";

            if (std::isfinite(Sdphi.meanY))
              std::cout << "      <|Δφ(γ,jet1)|> for matched events: " << std::fixed << std::setprecision(3) << Sdphi.meanY
                        << " rad (fills=" << static_cast<long long>(std::llround(Sdphi.inrange)) << ")\n";

            if (xjN > 0.0)
              std::cout << "      xJ fills=" << static_cast<long long>(std::llround(xjN))
                        << "  <xJ>=" << std::fixed << std::setprecision(3) << xjMean << "\n";

            if (j1N > 0.0)
              std::cout << "      jet1 pT fills=" << static_cast<long long>(std::llround(j1N))
                        << "  <pT_jet1>=" << std::fixed << std::setprecision(3) << j1Mean << " GeV\n";

            if (aN > 0.0)
              std::cout << "      alpha fills=" << static_cast<long long>(std::llround(aN))
                        << "  <alpha>=" << std::fixed << std::setprecision(3) << aMean << "\n";

            std::cout << "      Unfolding (inclusive pairing):\n";
            std::cout << "        reco pairs:   entries=" << static_cast<long long>(std::llround(Ureco.entries))
                      << "  in-range=" << static_cast<long long>(std::llround(Ureco.inrange))
                      << "  <xJ>=" << (std::isfinite(Ureco.meanY) ? (static_cast<std::ostringstream&&>(std::ostringstream() << std::fixed << std::setprecision(3) << Ureco.meanY).str()) : std::string("n/a"))
                      << "\n";
            std::cout << "        truth pairs:  entries=" << static_cast<long long>(std::llround(Utruth.entries))
                      << "  in-range=" << static_cast<long long>(std::llround(Utruth.inrange))
                      << "  <xJ>=" << (std::isfinite(Utruth.meanY) ? (static_cast<std::ostringstream&&>(std::ostringstream() << std::fixed << std::setprecision(3) << Utruth.meanY).str()) : std::string("n/a"))
                      << "\n";
            std::cout << "        response fills (global-bin TH2): entries=" << static_cast<long long>(std::llround(rspEntries)) << "\n";
            std::cout << "        reco fakes:   entries=" << static_cast<long long>(std::llround(Ufakes.entries))
                      << "  in-range=" << static_cast<long long>(std::llround(Ufakes.inrange)) << "\n";
            std::cout << "        truth misses: entries=" << static_cast<long long>(std::llround(Umiss.entries))
                      << "  in-range=" << static_cast<long long>(std::llround(Umiss.inrange)) << "\n";
          }
        }
      }
    }


  //--------------------------------------------------------------------
  // 5. Write footer & close the file
  //--------------------------------------------------------------------
  if (Verbosity() >= 1)
  {
      std::cout << "\nOutput ROOT file →  " << out->GetName() << "\n\n";

      info(1, "writing TFile footer and closing (" + std::to_string(nHistWritten)
                 + " / " + std::to_string(nHistExpected) + " objects written)");
  }

    // EventDisplay diagnostics payload (offline rendering; independent of Verbosity()).
    if (m_evtDiagTree)
      {
        if (Verbosity() >= 1)
        {
          const Long64_t nED = m_evtDiagTree->GetEntries();

          std::cout << CLR_CYAN
                    << "\n[End][EventDisplayTree] entries=" << nED
                    << "  fillsRecorded=" << m_evtDiagNFill
                    << "  anyTowers=" << m_evtDiagNFillWithAnyTowers
                    << "  selTowers=" << m_evtDiagNFillWithSelTowers
                    << "  bestTowers=" << m_evtDiagNFillWithBestTowers
                    << CLR_RESET << "\n";

          std::cout << "  byCat: "
                    << "NUM fills=" << m_evtDiagNFillByCat[0] << " anyTowers=" << m_evtDiagNFillWithAnyTowersByCat[0]
                    << "  MissA fills=" << m_evtDiagNFillByCat[1] << " anyTowers=" << m_evtDiagNFillWithAnyTowersByCat[1]
                    << "  MissB fills=" << m_evtDiagNFillByCat[2] << " anyTowers=" << m_evtDiagNFillWithAnyTowersByCat[2]
                    << "\n";

          if (nED > 0 && m_evtDiagNFillWithAnyTowers == 0)
          {
            std::cout << CLR_YELLOW
                      << "  [WARN] EventDisplayTree has entries but ALL tower payloads are empty. "
                      << "Most likely: constituent decoding mismatch (TowerInfo node vs Jet SRC) or comp_vec stripped by calibration."
                      << CLR_RESET << "\n";
          }
        }

        out->cd();
        m_evtDiagTree->Write("", TObject::kOverwrite);
      }

    try
    {
      out->Close();
  }
  catch (const std::exception& e)
  {
    warn("Exception during TFile::Write/Close – " + std::string(e.what()));
  }

  delete out;
  out = nullptr;

  info(0, "Done.");
  return Fun4AllReturnCodes::EVENT_OK;
}



// Build SSVars from PhotonClusterv1 (names come from PhotonClusterBuilder)
RecoilJets::SSVars RecoilJets::makeSSFromPhoton(const PhotonClusterv1* pho, double pt_gamma) const
{
  SSVars v{};
  v.pt_gamma = pt_gamma;

  if (!pho)
  {
    if (Verbosity() >= 3)
      LOG(3, CLR_YELLOW, "  [makeSSFromPhoton] pho==nullptr → returning default SSVars");
    return v;
  }

  // Raw stored shower-shape parameters (PhotonClusterBuilder names)
  const double weta_cogx = pho->get_shower_shape_parameter("weta_cogx");
  const double wphi_cogx = pho->get_shower_shape_parameter("wphi_cogx");
  const double et1       = pho->get_shower_shape_parameter("et1");

  const double e11 = pho->get_shower_shape_parameter("e11");
  const double e33 = pho->get_shower_shape_parameter("e33");
  const double e32 = pho->get_shower_shape_parameter("e32");
  const double e35 = pho->get_shower_shape_parameter("e35");

  // Derived ratios (protect denominators)
  const double e11_over_e33_raw = (e33 > 0.0) ? (e11 / e33) : std::numeric_limits<double>::quiet_NaN();
  const double e32_over_e35_raw = (e35 > 0.0) ? (e32 / e35) : std::numeric_limits<double>::quiet_NaN();

  v.weta_cogx    = weta_cogx;
  v.wphi_cogx    = wphi_cogx;
  v.et1          = et1;

  // Preserve NaN for invalid denominators/inputs.
  // This ensures preselection/tight cuts fail cleanly for invalid SS inputs.
  v.e11_over_e33 = e11_over_e33_raw;
  v.e32_over_e35 = e32_over_e35_raw;

  // Print only when:
  //   - Verbosity >= 8 (always), OR
  //   - Verbosity >= 6 and something looks wrong (non-finite / denom<=0)
  const bool bad =
      !std::isfinite(weta_cogx) || !std::isfinite(wphi_cogx) || !std::isfinite(et1) ||
      !std::isfinite(e11) || !std::isfinite(e33) || !std::isfinite(e32) || !std::isfinite(e35) ||
      (e33 <= 0.0) || (e35 <= 0.0);

  if (Verbosity() >= 8 || (Verbosity() >= 6 && bad))
  {
    std::ostringstream os;
    os << "      [SS extract] pT^γ=" << std::fixed << std::setprecision(2) << pt_gamma
       << " | weta_cogx=" << std::setprecision(4) << weta_cogx
       << " wphi_cogx=" << std::setprecision(4) << wphi_cogx
       << " et1=" << std::setprecision(4) << et1
       << " | e11=" << std::setprecision(4) << e11
       << " e33=" << std::setprecision(4) << e33
       << " (e11/e33=" << std::setprecision(4) << (std::isfinite(e11_over_e33_raw) ? e11_over_e33_raw : -999.0) << ")"
       << " | e32=" << std::setprecision(4) << e32
       << " e35=" << std::setprecision(4) << e35
       << " (e32/e35=" << std::setprecision(4) << (std::isfinite(e32_over_e35_raw) ? e32_over_e35_raw : -999.0) << ")";
    if (bad) os << "  [WARN]";
    LOG(6, bad ? CLR_YELLOW : CLR_BLUE, os.str());
  }

  return v;
}



void RecoilJets::fillUnfoldResponseMatrixAndTruthDistributions(
    const std::vector<std::string>& activeTrig,
    const std::string& rKey,
    const int effCentIdx_M,
    const double leadPtGamma,
    const double tPt,
    const double tPhi,
    const std::vector<const Jet*>& recoJetsFid,
    const std::vector<char>& recoJetsFidIsRecoil,
    const Jet* recoil1Jet)
{
  // Keep local dR identical in behavior to the caller
  auto dR = [](double eta1, double phi1, double eta2, double phi2) -> double
  {
    const double dphi = TVector2::Phi_mpi_pi(phi1 - phi2);
    const double deta = eta1 - eta2;
    return std::sqrt(deta*deta + dphi*dphi);
  };

  JetContainer* truthJets = nullptr;
  if (auto itT = m_truthJetsByRKey.find(rKey); itT != m_truthJetsByRKey.end()) truthJets = itT->second;
  if (!truthJets) return;

  const double etaMaxTruth = jetEtaAbsMaxForRKey(rKey);

  std::vector<const Jet*> truthJetsFid;
  std::vector<char>       truthJetsFidIsRecoil;
  truthJetsFid.reserve(truthJets->size());
  truthJetsFidIsRecoil.reserve(truthJets->size());

  // 1) Build truth jet lists + fill truth inclusive distribution
  for (const Jet* tj : *truthJets)
  {
    if (!tj) continue;

    const double ptj  = tj->get_pt();
    const double etaj = tj->get_eta();
    const double phij = tj->get_phi();

    if (!std::isfinite(ptj) || !std::isfinite(etaj) || !std::isfinite(phij)) continue;
    if (ptj < m_minJetPt) continue;
    if (std::fabs(etaj) >= etaMaxTruth) continue;

    const double dphiAbs = std::fabs(TVector2::Phi_mpi_pi(phij - tPhi));
    const bool isRecoil = (dphiAbs >= m_minBackToBack);

    truthJetsFid.push_back(tj);
    truthJetsFidIsRecoil.push_back(isRecoil ? 1 : 0);

    if (isRecoil)
    {
      const double xJt = ptj / tPt;

      for (const auto& trigShort : activeTrig)
      {
        if (auto* h2t = getOrBookUnfoldTruthPtXJIncl(trigShort, rKey, effCentIdx_M))
        {
          h2t->Fill(tPt, xJt);
          bumpHistFill(trigShort, h2t->GetName());
        }

        // inclusive |Δphi(γ,jet)| for each TRUTH recoil jet passing cuts
        if (auto* h2dt = getOrBookUnfoldTruthPtDphiIncl(trigShort, rKey, effCentIdx_M))
        {
          h2dt->Fill(tPt, dphiAbs);
          bumpHistFill(trigShort, h2dt->GetName());
        }
      }
    }
  }

    // ------------------------------------------------------------------
    //  leading-truth recoil jet1 match bookkeeping vs truth pT^gamma
    // ------------------------------------------------------------------
    {
      const double kLeadMatchDR = m_jetMatchDRMax;

      int iTruthLead = -1;
      double ptTruthLead = -1.0;

      for (int it = 0; it < static_cast<int>(truthJetsFid.size()); ++it)
      {
        if (!truthJetsFidIsRecoil[it]) continue;

        const Jet* tj = truthJetsFid[it];
        if (!tj) continue;

        const double ptj = tj->get_pt();
        if (!std::isfinite(ptj) || ptj <= 0.0) continue;

        if (ptj > ptTruthLead)
        {
          ptTruthLead = ptj;
          iTruthLead  = it;
        }
      }

      if (iTruthLead >= 0)
      {
        const Jet* tjLead = truthJetsFid[iTruthLead];

        // DEN: truth-leading away-side recoil jet exists
        for (const auto& trigShort : activeTrig)
        {
          if (auto* hDen = getOrBookLeadTruthRecoilMatchDenPtGammaTruth(trigShort, rKey, effCentIdx_M))
          { hDen->Fill(tPt); bumpHistFill(trigShort, hDen->GetName()); }
        }

        // Best reco fid match to truth-leading recoil jet (independent of the selected recoil1Jet)
        const Jet* rjTruthBest = nullptr;
        double drTruthBest = 1e9;
        if (tjLead)
        {
          for (const Jet* rjAny : recoJetsFid)
          {
            if (!rjAny) continue;

            const double drAny = dR(rjAny->get_eta(), rjAny->get_phi(),
                                    tjLead->get_eta(), tjLead->get_phi());
            if (drAny < drTruthBest)
            {
              drTruthBest = drAny;
              rjTruthBest = rjAny;
            }
          }
        }
        const bool hasRecoMatchToTruthLead = (rjTruthBest && drTruthBest < kLeadMatchDR);

        // For MissA subtyping: does the reco jet matched to the truth-leading recoil jet
        // itself pass the recoil definition? (pT/eta are already satisfied by recoJetsFid,
        // but we evaluate sequentially to support a robust MissA2 cutflow histogram.)
        bool truthMatchPassesRecoil = false;
        int  truthMatchFailCut = 0; // 1=pTmin, 2=|eta|, 3=dphi
        if (hasRecoMatchToTruthLead && rjTruthBest)
        {
          const double ptm  = rjTruthBest->get_pt();
          const double etam = rjTruthBest->get_eta();
          const double phim = rjTruthBest->get_phi();

          const bool passPt  = (std::isfinite(ptm)  && ptm  >= m_minJetPt);
          const bool passEta = (std::isfinite(etam) && std::fabs(etam) < etaMaxTruth);

          const double dphiAbsM = std::fabs(TVector2::Phi_mpi_pi(phim - tPhi));
          const bool passDphi = (std::isfinite(dphiAbsM) && dphiAbsM >= m_minBackToBack);

          truthMatchPassesRecoil = (passPt && passEta && passDphi);

          if (!passPt)        truthMatchFailCut = 1;
          else if (!passEta)  truthMatchFailCut = 2;
          else if (!passDphi) truthMatchFailCut = 3;
          else                truthMatchFailCut = 0;
        }

        // Does the analysis-selected reco recoilJet1 match the truth-leading recoil jet?
        bool leadRecoMatches = false;
        double drLead = 1e9;
        if (recoil1Jet && tjLead)
        {
          drLead = dR(recoil1Jet->get_eta(), recoil1Jet->get_phi(),
                      tjLead->get_eta(), tjLead->get_phi());
          leadRecoMatches = (drLead < kLeadMatchDR);
        }

        // Fill the standard NUM / MissA / MissB bookkeeping (+ MissA subtypes)
        if (leadRecoMatches)
        {
            for (const auto& trigShort : activeTrig)
            {
              if (auto* hNum = getOrBookLeadTruthRecoilMatchNumPtGammaTruth(trigShort, rKey, effCentIdx_M))
              { hNum->Fill(tPt); bumpHistFill(trigShort, hNum->GetName()); }
            }
          }
          else if (hasRecoMatchToTruthLead)
          {
            for (const auto& trigShort : activeTrig)
            {
              if (auto* hA = getOrBookLeadTruthRecoilMatchMissA_PtGammaTruth(trigShort, rKey, effCentIdx_M))
              { hA->Fill(tPt); bumpHistFill(trigShort, hA->GetName()); }

              if (truthMatchPassesRecoil)
              {
                if (auto* hA1 = getOrBookLeadTruthRecoilMatchMissA1_PtGammaTruth(trigShort, rKey, effCentIdx_M))
                { hA1->Fill(tPt); bumpHistFill(trigShort, hA1->GetName()); }
              }
              else
              {
                if (auto* hA2 = getOrBookLeadTruthRecoilMatchMissA2_PtGammaTruth(trigShort, rKey, effCentIdx_M))
                { hA2->Fill(tPt); bumpHistFill(trigShort, hA2->GetName()); }

                if (truthMatchFailCut > 0)
                {
                  if (auto* hCF = getOrBookLeadTruthRecoilMatchMissA2_Cutflow(trigShort, rKey, effCentIdx_M))
                  { hCF->Fill((double)truthMatchFailCut); bumpHistFill(trigShort, hCF->GetName()); }
                }
              }
            }
          }
          else
          {
            for (const auto& trigShort : activeTrig)
            {
              if (auto* hB = getOrBookLeadTruthRecoilMatchMissB_PtGammaTruth(trigShort, rKey, effCentIdx_M))
              { hB->Fill(tPt); bumpHistFill(trigShort, hB->GetName()); }
            }
          }

           // -------------------------------------------------------------------------
            // EventDisplay diagnostics payload (offline rendering; independent of Verbosity()).
            // One entry per (event, rKey) when enabled.
            // -------------------------------------------------------------------------
            if (m_evtDiagEnabled)
            {
              if (!m_evtDiagNodesReady)
              {
                if (Verbosity() >= 3)
                {
                  LOG(3, CLR_YELLOW,
                      "[EventDisplayTree][skip] nodes not ready "
                      << "(evt=" << ecvt << " rKey=" << rKey << ")");
                }
              }
              else
              {
                const int ptBin = findPtBin(tPt);

                EventDisplayCat cat = EventDisplayCat::MissB;
                if (leadRecoMatches)
                {
                  cat = EventDisplayCat::NUM;
                }
                else if (hasRecoMatchToTruthLead)
                {
                  cat = EventDisplayCat::MissA;
                }

                if (ptBin < 0)
                {
                  if (Verbosity() >= 4)
                  {
                    LOG(4, CLR_YELLOW,
                        "[EventDisplayTree][skip] ptBin<0 (out of range) "
                        << "(evt=" << ecvt << " rKey=" << rKey
                        << " tPt=" << tPt << ")");
                  }
                }
                else if (!(recoil1Jet || rjTruthBest || tjLead))
                {
                  if (Verbosity() >= 4)
                  {
                    LOG(4, CLR_YELLOW,
                        "[EventDisplayTree][skip] no jets available for diagnostics "
                        << "(evt=" << ecvt << " rKey=" << rKey
                        << " ptBin=" << ptBin
                        << " cat=" << static_cast<int>(cat) << ")");
                  }
                }
                else
                {
                  const bool need = eventDisplayDiagnosticsNeed(rKey, ptBin, cat);

                  if (Verbosity() >= 5)
                  {
                    LOG(5, CLR_CYAN,
                        "[EventDisplayTree][check] "
                        << "(evt=" << ecvt << " rKey=" << rKey
                        << " ptBin=" << ptBin
                        << " cat=" << static_cast<int>(cat)
                        << " need=" << (need ? "YES" : "NO")
                        << " sel=" << (recoil1Jet ? "Y" : "N")
                        << " best=" << (rjTruthBest ? "Y" : "N")
                        << " truthLead=" << (tjLead ? "Y" : "N")
                        << " tPt=" << tPt
                        << " leadPtGamma=" << leadPtGamma << ")");
                  }

                  if (need)
                  {
                    if (Verbosity() >= 4)
                    {
                      LOG(4, CLR_CYAN,
                          "[EventDisplayTree][fill] dispatching fillEventDisplayDiagnostics "
                          << "(evt=" << ecvt << " rKey=" << rKey
                          << " ptBin=" << ptBin
                          << " cat=" << static_cast<int>(cat) << ")");
                    }

                    fillEventDisplayDiagnostics(rKey, ptBin, cat, tPt, tPhi, leadPtGamma, recoil1Jet, rjTruthBest, tjLead);

                    if (Verbosity() >= 4)
                    {
                      LOG(4, CLR_GREEN,
                          "[EventDisplayTree][fill] returned from fillEventDisplayDiagnostics "
                          << "(evt=" << ecvt << " rKey=" << rKey
                          << " ptBin=" << ptBin
                          << " cat=" << static_cast<int>(cat) << ")");
                    }
                  }
                }
              }
            }

        // ------------------------------------------------------------------
        // diagnostics (SIM-only): characterize the selected recoilJet1 in each class.
        //
        // Filled only when recoil1Jet exists (i.e. the event actually contributes to reco xJ).
        // ------------------------------------------------------------------
        if (recoil1Jet && tjLead)
        {
          const double ptReco1  = recoil1Jet->get_pt();
          const double ptTruth  = tjLead->get_pt();

          const double phiReco1 = recoil1Jet->get_phi();

          if (std::isfinite(ptReco1) && std::isfinite(ptTruth) && ptReco1 > 0.0 && ptTruth > 0.0)
          {
            // Use truth photon axis (tPhi) to keep this within the truth-conditioned block.
            // (haveTruthPho ensures reco and truth photon directions are aligned event-by-event.)
            const double dphiReco1 = std::fabs(TVector2::Phi_mpi_pi(phiReco1 - tPhi));
            const double dRReco1TruthLead = drLead;

            const double xJReco1 = (leadPtGamma > 0.0 ? (ptReco1 / leadPtGamma) : 0.0);

            for (const auto& trigShort : activeTrig)
            {
              if (leadRecoMatches)
              {
                // (A1) pT(recoilJet1^reco) vs pT(truth-leading recoil jet)
                if (auto* h = getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtTruthLead_num(trigShort, rKey, effCentIdx_M))
                { h->Fill(ptTruth, ptReco1); bumpHistFill(trigShort, h->GetName()); }

                // (A2) pT(recoilJet1^reco) vs pT(reco jet matched to truth-leading recoil jet)
                if (hasRecoMatchToTruthLead && rjTruthBest)
                {
                  const double ptRecoMatch = rjTruthBest->get_pt();
                  if (std::isfinite(ptRecoMatch) && ptRecoMatch > 0.0)
                  {
                    if (auto* h = getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtRecoTruthMatch_num(trigShort, rKey, effCentIdx_M))
                    { h->Fill(ptRecoMatch, ptReco1); bumpHistFill(trigShort, h->GetName()); }
                  }
                }

                // (B3) |Δphi(γ^truth, recoilJet1^reco)| vs pTγ,truth
                if (auto* h = getOrBookLeadTruthRecoilMatchDphiRecoJet1VsPtGammaTruth_num(trigShort, rKey, effCentIdx_M))
                { h->Fill(tPt, dphiReco1); bumpHistFill(trigShort, h->GetName()); }

                // (B4) ΔR(recoilJet1^reco, truth-leading recoil jet) vs pTγ,truth
                if (auto* h = getOrBookLeadTruthRecoilMatchDRRecoJet1VsPtGammaTruth_num(trigShort, rKey, effCentIdx_M))
                { h->Fill(tPt, dRReco1TruthLead); bumpHistFill(trigShort, h->GetName()); }

                // (C5) xJ(recoilJet1^reco) vs |Δphi(γ^truth, recoilJet1^reco)|
                if (auto* h = getOrBookLeadTruthRecoilMatchXJRecoJet1VsDphiRecoJet1_num(trigShort, rKey, effCentIdx_M))
                { h->Fill(dphiReco1, xJReco1); bumpHistFill(trigShort, h->GetName()); }
              }
              else if (hasRecoMatchToTruthLead)
              {
                // MissA: truth-leading recoil jet is reconstructed, but we selected the wrong reco recoilJet1

                if (auto* h = getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtTruthLead_missA(trigShort, rKey, effCentIdx_M))
                { h->Fill(ptTruth, ptReco1); bumpHistFill(trigShort, h->GetName()); }

                if (truthMatchPassesRecoil)
                {
                  if (auto* h = getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtTruthLead_missA1(trigShort, rKey, effCentIdx_M))
                  { h->Fill(ptTruth, ptReco1); bumpHistFill(trigShort, h->GetName()); }
                }
                else
                {
                  if (auto* h = getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtTruthLead_missA2(trigShort, rKey, effCentIdx_M))
                  { h->Fill(ptTruth, ptReco1); bumpHistFill(trigShort, h->GetName()); }
                }

                if (rjTruthBest)
                {
                  const double ptRecoMatch = rjTruthBest->get_pt();
                  if (std::isfinite(ptRecoMatch) && ptRecoMatch > 0.0)
                  {
                    if (auto* h = getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtRecoTruthMatch_missA(trigShort, rKey, effCentIdx_M))
                    { h->Fill(ptRecoMatch, ptReco1); bumpHistFill(trigShort, h->GetName()); }

                    if (truthMatchPassesRecoil)
                    {
                      if (auto* h = getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtRecoTruthMatch_missA1(trigShort, rKey, effCentIdx_M))
                      { h->Fill(ptRecoMatch, ptReco1); bumpHistFill(trigShort, h->GetName()); }
                    }
                    else
                    {
                      if (auto* h = getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtRecoTruthMatch_missA2(trigShort, rKey, effCentIdx_M))
                      { h->Fill(ptRecoMatch, ptReco1); bumpHistFill(trigShort, h->GetName()); }
                    }
                  }
                }

                if (auto* h = getOrBookLeadTruthRecoilMatchDphiRecoJet1VsPtGammaTruth_missA(trigShort, rKey, effCentIdx_M))
                { h->Fill(tPt, dphiReco1); bumpHistFill(trigShort, h->GetName()); }

                if (auto* h = getOrBookLeadTruthRecoilMatchDRRecoJet1VsPtGammaTruth_missA(trigShort, rKey, effCentIdx_M))
                { h->Fill(tPt, dRReco1TruthLead); bumpHistFill(trigShort, h->GetName()); }

                if (auto* h = getOrBookLeadTruthRecoilMatchXJRecoJet1VsDphiRecoJet1_missA(trigShort, rKey, effCentIdx_M))
                { h->Fill(dphiReco1, xJReco1); bumpHistFill(trigShort, h->GetName()); }
              }
              else
              {
                // MissB: truth-leading recoil jet has no reco match among recoJetsFid

                if (auto* h = getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtTruthLead_missB(trigShort, rKey, effCentIdx_M))
                { h->Fill(ptTruth, ptReco1); bumpHistFill(trigShort, h->GetName()); }

                if (auto* h = getOrBookLeadTruthRecoilMatchDphiRecoJet1VsPtGammaTruth_missB(trigShort, rKey, effCentIdx_M))
                { h->Fill(tPt, dphiReco1); bumpHistFill(trigShort, h->GetName()); }

                if (auto* h = getOrBookLeadTruthRecoilMatchDRRecoJet1VsPtGammaTruth_missB(trigShort, rKey, effCentIdx_M))
                { h->Fill(tPt, dRReco1TruthLead); bumpHistFill(trigShort, h->GetName()); }

                if (auto* h = getOrBookLeadTruthRecoilMatchXJRecoJet1VsDphiRecoJet1_missB(trigShort, rKey, effCentIdx_M))
                { h->Fill(dphiReco1, xJReco1); bumpHistFill(trigShort, h->GetName()); }
              }
            }
          }
        }
      }
    }

  // 2) Jet ΔR matching (fiducial jets only) to build response + fakes/misses
  struct CandPair { double dr; int iReco; int iTruth; };
  std::vector<CandPair> cands;
  cands.reserve(recoJetsFid.size() * truthJetsFid.size());

  for (int ir = 0; ir < static_cast<int>(recoJetsFid.size()); ++ir)
  {
    const Jet* rj = recoJetsFid[ir];
    if (!rj) continue;

    for (int it = 0; it < static_cast<int>(truthJetsFid.size()); ++it)
    {
      const Jet* tj = truthJetsFid[it];
      if (!tj) continue;

      const double dr = dR(rj->get_eta(), rj->get_phi(), tj->get_eta(), tj->get_phi());
      if (dr < 0.3) cands.push_back({dr, ir, it});
    }
  }

  std::sort(cands.begin(), cands.end(),
            [](const CandPair& a, const CandPair& b){ return a.dr < b.dr; });

  std::vector<char> recoMatched(recoJetsFid.size(), 0);
  std::vector<char> truthMatched(truthJetsFid.size(), 0);

  for (const auto& trigShort : activeTrig)
  {
    // IMPORTANT:
    // recoMatched/truthMatched must be per-trigger bookkeeping because we fill a
    // complete set of histograms per trigger. If activeTrig has multiple entries,
    // we must reset these flags each iteration to avoid cross-trigger state leakage.
    std::fill(recoMatched.begin(), recoMatched.end(), 0);
    std::fill(truthMatched.begin(), truthMatched.end(), 0);

    auto* h2Reco  = getOrBookUnfoldRecoPtXJIncl(trigShort, rKey, effCentIdx_M);
    auto* h2Truth = getOrBookUnfoldTruthPtXJIncl(trigShort, rKey, effCentIdx_M);
    auto* hRsp    = getOrBookUnfoldResponsePtXJIncl(trigShort, rKey, effCentIdx_M);
    auto* hFake   = getOrBookUnfoldRecoFakesPtXJIncl(trigShort, rKey, effCentIdx_M);
    auto* hMiss   = getOrBookUnfoldTruthMissesPtXJIncl(trigShort, rKey, effCentIdx_M);

    auto* h2RecoM  = getOrBookUnfoldRecoMatchedPtXJIncl(trigShort, rKey, effCentIdx_M);
    auto* h2TruthM = getOrBookUnfoldTruthMatchedPtXJIncl(trigShort, rKey, effCentIdx_M);

    auto* hFakeA   = getOrBookUnfoldRecoFakesPtXJIncl_typeA(trigShort, rKey, effCentIdx_M);
    auto* hFakeB   = getOrBookUnfoldRecoFakesPtXJIncl_typeB(trigShort, rKey, effCentIdx_M);
    auto* hMissA   = getOrBookUnfoldTruthMissesPtXJIncl_typeA(trigShort, rKey, effCentIdx_M);
    auto* hMissB   = getOrBookUnfoldTruthMissesPtXJIncl_typeB(trigShort, rKey, effCentIdx_M);

    auto* hDR           = getOrBookUnfoldJetMatchDR(trigShort, rKey, effCentIdx_M);
    auto* hPtResp       = getOrBookUnfoldJetPtResponsePtTruth(trigShort, rKey, effCentIdx_M);
    auto* hPtRespAll    = getOrBookUnfoldJetPtResponseAllPtTruth(trigShort, rKey, effCentIdx_M);
    auto* hPtRespLead   = getOrBookLeadRecoilJetPtResponsePtTruth(trigShort, rKey, effCentIdx_M);
    auto* hPtScatterLead= getOrBookLeadRecoilJetPtTruthPtReco(trigShort, rKey, effCentIdx_M);
    auto* hDRLead       = getOrBookLeadRecoilJetMatchDR(trigShort, rKey, effCentIdx_M);

    if (!h2Reco || !h2Truth || !hRsp || !hFake || !hMiss) continue;

    // --- greedy one-to-one matching ---
    for (const auto& cp : cands)
    {
        if (recoMatched[cp.iReco]) continue;
        if (truthMatched[cp.iTruth]) continue;

        recoMatched[cp.iReco]   = 1;
        truthMatched[cp.iTruth] = 1;

        const Jet* rj = recoJetsFid[cp.iReco];
        const Jet* tj = truthJetsFid[cp.iTruth];
        if (!rj || !tj) continue;

        const bool recoSel  = (recoJetsFidIsRecoil[cp.iReco] != 0);
        const bool truthSel = (truthJetsFidIsRecoil[cp.iTruth] != 0);

        const double xJr = rj->get_pt() / leadPtGamma;
        const double xJt = tj->get_pt() / tPt;

        // --- inclusive matched-fiducial pT response (no truthSel/recoSel gating) ---
        if (hPtRespAll && tj->get_pt() > 0.0)
        {
          hPtRespAll->Fill(tj->get_pt(), rj->get_pt() / tj->get_pt());
          bumpHistFill(trigShort, hPtRespAll->GetName());
        }

        // --- lead recoil jet1 response + ΔR sanity + pT scatter ---
        if (recoil1Jet && (rj == recoil1Jet))
        {
          if (hDRLead)
          { hDRLead->Fill(cp.dr); bumpHistFill(trigShort, hDRLead->GetName()); }

          if (tj->get_pt() > 0.0)
          {
            if (hPtRespLead)
            { hPtRespLead->Fill(tj->get_pt(), rj->get_pt() / tj->get_pt()); bumpHistFill(trigShort, hPtRespLead->GetName()); }

            if (hPtScatterLead)
            { hPtScatterLead->Fill(tj->get_pt(), rj->get_pt()); bumpHistFill(trigShort, hPtScatterLead->GetName()); }
          }
        }

        if (truthSel && recoSel)
        {
          const int gTruth = h2Truth->FindBin(tPt, xJt);
          const int gReco  = h2Reco->FindBin(leadPtGamma, xJr);

          hRsp->Fill(static_cast<double>(gTruth), static_cast<double>(gReco));
          bumpHistFill(trigShort, hRsp->GetName());

          if (h2TruthM)
          { h2TruthM->Fill(tPt, xJt); bumpHistFill(trigShort, h2TruthM->GetName()); }

          if (h2RecoM)
          { h2RecoM->Fill(leadPtGamma, xJr); bumpHistFill(trigShort, h2RecoM->GetName()); }

          if (hDR)
          { hDR->Fill(cp.dr); bumpHistFill(trigShort, hDR->GetName()); }

          if (hPtResp && tj->get_pt() > 0.0)
          { hPtResp->Fill(tj->get_pt(), rj->get_pt() / tj->get_pt()); bumpHistFill(trigShort, hPtResp->GetName()); }
      }
      else if (truthSel && !recoSel)
      {
        hMiss->Fill(tPt, xJt);
        bumpHistFill(trigShort, hMiss->GetName());

        if (hMissA)
        { hMissA->Fill(tPt, xJt); bumpHistFill(trigShort, hMissA->GetName()); }
      }
      else if (!truthSel && recoSel)
      {
        hFake->Fill(leadPtGamma, xJr);
        bumpHistFill(trigShort, hFake->GetName());

        if (hFakeA)
        { hFakeA->Fill(leadPtGamma, xJr); bumpHistFill(trigShort, hFakeA->GetName()); }
      }
    }

    // --- unmatched truth recoil jets => MISS type B ---
    for (std::size_t it = 0; it < truthJetsFid.size(); ++it)
    {
      if (truthMatched[it]) continue;
      if (!truthJetsFidIsRecoil[it]) continue;

      const Jet* tj = truthJetsFid[it];
      if (!tj) continue;

      const double xJt = tj->get_pt() / tPt;

      hMiss->Fill(tPt, xJt);
      bumpHistFill(trigShort, hMiss->GetName());

      if (hMissB)
      { hMissB->Fill(tPt, xJt); bumpHistFill(trigShort, hMissB->GetName()); }
    }

    // --- unmatched reco recoil jets => FAKE type B ---
    for (std::size_t ir = 0; ir < recoJetsFid.size(); ++ir)
    {
      if (recoMatched[ir]) continue;
      if (!recoJetsFidIsRecoil[ir]) continue;

      const Jet* rj = recoJetsFid[ir];
      if (!rj) continue;

      const double xJr = rj->get_pt() / leadPtGamma;

      hFake->Fill(leadPtGamma, xJr);
      bumpHistFill(trigShort, hFake->GetName());

      if (hFakeB)
      { hFakeB->Fill(leadPtGamma, xJr); bumpHistFill(trigShort, hFakeB->GetName()); }
    }
  }
}



void RecoilJets::fillRecoTruthJES3MatchingQA(const std::vector<std::string>& activeTrig,
                                            const std::string& rKey,
                                            const int effCentIdx_M,
                                            const double leadPtGamma,
                                            const double xJ,
                                            const double alpha,
                                            const double tPt,
                                            const double tPhi,
                                            const Jet* recoil1Jet)
{
  // Local dR helper (identical behavior to the lambda you use in the caller)
  auto dR = [](double eta1, double phi1, double eta2, double phi2) -> double
  {
    const double dphi = TVector2::Phi_mpi_pi(phi1 - phi2);
    const double deta = eta1 - eta2;
    return std::sqrt(deta*deta + dphi*dphi);
  };

  //  RECO truth-PHOTON-tagged subset (photon match only; no truth-jet match required)
  for (const auto& trigShort : activeTrig)
  {
    if (auto* hTagPho = getOrBookJES3RecoTruthPhoTagged_xJ_alphaHist(trigShort, rKey, effCentIdx_M))
    { hTagPho->Fill(leadPtGamma, xJ, alpha); bumpHistFill(trigShort, hTagPho->GetName()); }
  }

  JetContainer* truthJets = nullptr;
  if (auto itT = m_truthJetsByRKey.find(rKey); itT != m_truthJetsByRKey.end()) truthJets = itT->second;

  if (!truthJets)
  {
    if (Verbosity() >= 4)
      LOG(4, CLR_YELLOW, "      [truthQA] rKey=" << rKey << " truth jets missing → skip truth JES3 fills");
    return;
  }
  else if (tPt <= 0.0)
  {
    if (Verbosity() >= 4)
      LOG(4, CLR_YELLOW, "      [truthQA] rKey=" << rKey << " matched truth gamma has invalid pT → skip truth JES3 fills");
    return;
  }

  // truth recoil jets (same R), |eta| < 1.1 - R, away-side to truth gamma
  const double etaMaxTruth = jetEtaAbsMaxForRKey(rKey);

  double tj1Pt = -1.0;
  const Jet* tj1 = nullptr;

  // 1) truth recoil jet1: highest-pT fiducial jet back-to-back to truth gamma
  for (const Jet* tj : *truthJets)
  {
    if (!tj) continue;

    const double ptj  = tj->get_pt();
    const double etaj = tj->get_eta();
    const double phij = tj->get_phi();

    if (!std::isfinite(ptj) || !std::isfinite(etaj) || !std::isfinite(phij)) continue;
    if (ptj < m_minJetPt) continue;
    if (std::fabs(etaj) >= etaMaxTruth) continue;

    const double dphiAbs = std::fabs(TVector2::Phi_mpi_pi(phij - tPhi));
    if (dphiAbs < m_minBackToBack) continue;

    if (ptj > tj1Pt)
    {
      tj1Pt = ptj;
      tj1   = tj;
    }
  }

  if (!tj1 || tj1Pt <= 0.0)
  {
    if (Verbosity() >= 5)
      LOG(5, CLR_YELLOW, "      [truthQA] rKey=" << rKey << " no truth recoil jet1 found → skip truth JES3 fills");
    return;
  }

  // 2) truth jet2 for α: highest-pT fiducial truth jet excluding tj1 (NO Δφ requirement)
  double tj2Pt = -1.0;
  for (const Jet* tj : *truthJets)
  {
    if (!tj) continue;

    const double ptj  = tj->get_pt();
    const double etaj = tj->get_eta();
    const double phij = tj->get_phi();

    if (!std::isfinite(ptj) || !std::isfinite(etaj) || !std::isfinite(phij)) continue;
    if (ptj < m_minJetPt) continue;
    if (std::fabs(etaj) >= etaMaxTruth) continue;
    if (tj == tj1) continue;

    if (ptj > tj2Pt) tj2Pt = ptj;
  }

  const double xJt = tj1Pt / tPt;
  const double aT  = (tj2Pt > 0.0 ? (tj2Pt / tPt) : 0.0);

  // TRUTH reco-conditioned but NO reco↔truth jet match requirement
  for (const auto& trigShort : activeTrig)
  {
    if (auto* ht3_nm = getOrBookJES3TruthRecoCondNoJetMatch_xJ_alphaHist(trigShort, rKey, effCentIdx_M))
    { ht3_nm->Fill(tPt, xJt, aT); bumpHistFill(trigShort, ht3_nm->GetName()); }
  }

  // require reco jet1 ↔ truth jet1 match in ΔR for the "jet-matched truth" JES3
  const double recoJeta = recoil1Jet->get_eta();
  const double recoJphi = recoil1Jet->get_phi();
  const double drJet = dR(recoJeta, recoJphi, tj1->get_eta(), tj1->get_phi());

  if (drJet > 0.3)
  {
    if (Verbosity() >= 5)
      LOG(5, CLR_YELLOW, "      [truthQA] rKey=" << rKey
                      << " reco jet1 ↔ truth jet1 ΔR=" << drJet << " > 0.3 → skip jet-matched truth JES3 fills");
    return;
  }

  for (const auto& trigShort : activeTrig)
  {
    // existing: TRUTH reco-conditioned + jet-matched
    if (auto* ht3x = getOrBookJES3Truth_xJ_alphaHist(trigShort, rKey, effCentIdx_M))
    { ht3x->Fill(tPt, xJt, aT); bumpHistFill(trigShort, ht3x->GetName()); }

    if (auto* ht3j = getOrBookJES3Truth_jet1Pt_alphaHist(trigShort, rKey, effCentIdx_M))
    { ht3j->Fill(tPt, tj1Pt, aT); bumpHistFill(trigShort, ht3j->GetName()); }

    // RECO truth-tagged (photon match + jet match)
    if (auto* hRecoTag = getOrBookJES3RecoTruthTagged_xJ_alphaHist(trigShort, rKey, effCentIdx_M))
    { hRecoTag->Fill(leadPtGamma, xJ, alpha); bumpHistFill(trigShort, hRecoTag->GetName()); }
  }
}


bool RecoilJets::runLeadIsoTightPhotonJetLoopAllRadii(
    const std::vector<std::string>& activeTrig,
    const int effCentIdx_M,
    const int centIdx,
    const int leadPhoIndex,
    const int leadPtIdx,
    const double leadPtGamma,
    const double leadEtaGamma,
    const double leadPhiGamma,
    const bool haveTruthPho,
    const double tPt,
    const double tPhi,
    PHG4TruthInfoContainer* truth)
{
  bool filledAnyRadius = false;

    // Run the EXACT SAME jet logic for every configured radius in parallel
  for (const auto& kv : m_jets)
  {
    const std::string rKey = kv.first;

    JetContainer* jets = kv.second;

    // Radius-dependent containment cut: |eta_jet| < 1.1 - R
    const double Rjet            = jetRFromKey(rKey);
    const double jetEtaAbsMaxUse = jetEtaAbsMaxForRKey(rKey);

    if (!jets)
    {
      LOG(3, CLR_YELLOW,
          "      [lead pho#" << leadPhoIndex << "] rKey=" << rKey
          << " jet container is nullptr – cannot form xJ");

      // Still fill matching-QA once per event-leading iso∧tight photon (per rKey)
      const int    status = 1;          // NoJetPt
      const double maxDphiFill = -0.01; // underflow sentinel

      for (const auto& trigShort : activeTrig)
      {
        if (auto* hs = getOrBookMatchStatusVsPtGamma(trigShort, rKey, effCentIdx_M))
        { hs->Fill(leadPtGamma, status); bumpHistFill(trigShort, hs->GetName()); }

        if (auto* hm = getOrBookMatchMaxDphiVsPtGamma(trigShort, rKey, effCentIdx_M))
        { hm->Fill(leadPtGamma, maxDphiFill); bumpHistFill(trigShort, hm->GetName()); }

        if (auto* pn = getOrBookNRecoilJetsVsPtGamma(trigShort, rKey, effCentIdx_M))
        { pn->Fill(leadPtGamma, 0.0); bumpHistFill(trigShort, pn->GetName()); }
      }
      continue;
    }

    // Scan jets:
    //   - baseline pT gate:        pT > m_minJetPt
    //   - fiducial containment:    |eta_jet| < jetEtaAbsMaxUse
    //   - recoil selection:        |Δφ(γ,jet)| >= m_minBackToBack
    //
    // IMPORTANT:
    //   - recoil1/recoil2 are still the leading/subleading recoil jets (your existing alpha logic).
    //   - we also cache *all* fiducial jets and a recoil-flag per jet so we can do
    //     inclusive photon–jet pairing (ATLAS-like) + build the unfolding response.
    int nPassPt   = 0;
    int nPassEta  = 0;
    int nPassDphi = 0;

    // lists for inclusive pairing + truth matching response
    std::vector<const Jet*> recoJetsFid;          // jets passing pT+eta
    std::vector<char>       recoJetsFidIsRecoil;  // whether each fid jet passes the recoil Δφ cut
    recoJetsFid.reserve(jets->size());
    recoJetsFidIsRecoil.reserve(jets->size());

    // Leading recoil jet (pass pT+eta+dphi): defines xJgamma
    double recoil1Pt = -1.0;
    const Jet* recoil1Jet = nullptr;

    // Leading + subleading jets in the full fiducial set (pass pT+eta; NO Δφ requirement)
    // Used for:
    //   - recoil-is-leading QA (compare recoil1Jet vs all1Jet)
    //   - alpha jet2 definition: highest-pT fiducial jet excluding recoil1Jet
    double all1Pt = -1.0;
    const Jet* all1Jet = nullptr;
    double all2Pt = -1.0;
    [[maybe_unused]]const Jet* all2Jet = nullptr;

    // max |Δφ| over jets that pass pT+eta (even if they fail the Δφ cut)
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

      if (std::fabs(jeta) >= jetEtaAbsMaxUse) continue;
      ++nPassEta;

//      // Track leading+subleading jets in the pT+eta set (NO Δφ requirement)
//      if (jpt > all1Pt)
//      {
//        all2Pt  = all1Pt;
//        all2Jet = all1Jet;
//
//        all1Pt  = jpt;
//        all1Jet = j;
//      }
//      else if (jpt > all2Pt)
//      {
//        all2Pt  = jpt;
//        all2Jet = j;
//      }
        // Track leading+subleading jets in the pT+eta set (NO Δφ requirement)
        //   - ensure "leading jet" is not the photon: veto jets with ΔR(γ,jet) < 0.4
        // Photon–jet overlap veto: exclude jets near the photon direction (ΔR < 0.4)
                const double dphiPho = TVector2::Phi_mpi_pi(jphi - leadPhiGamma);
                const double detaPho = (jeta - leadEtaGamma);
                const double dRPho2  = (detaPho*detaPho + dphiPho*dphiPho);

                if (!std::isfinite(dRPho2) || (dRPho2 < (0.4 * 0.4)))
                {
                  continue;
                }

                // Track leading+subleading jets in the pT+eta set (NO Δφ requirement)
                if (jpt > all1Pt)
                {
                  all2Pt  = all1Pt;
                  all2Jet = all1Jet;

                  all1Pt  = jpt;
                  all1Jet = j;
                }
                else if (jpt > all2Pt)
                {
                  all2Pt  = jpt;
                  all2Jet = j;
                }


              const double dphiAbs = std::fabs(TVector2::Phi_mpi_pi(jphi - leadPhiGamma));
              if (std::isfinite(dphiAbs) && dphiAbs > maxDphi) maxDphi = dphiAbs;

              const bool isRecoil = (dphiAbs >= m_minBackToBack);

              // cache fiducial jets for inclusive pairing / response matching
              recoJetsFid.push_back(j);
              recoJetsFidIsRecoil.push_back(isRecoil ? 1 : 0);

              if (isRecoil)
              {
                ++nPassDphi;
            }
        }


      // Sam-style leading-jet definition:
      //   - pick the leading fiducial jet excluding photon overlap (all1Jet)
      //   - THEN require it be back-to-back to define recoil1Jet
      if (all1Jet)
      {
        const double dphiLeadAbs = std::fabs(TVector2::Phi_mpi_pi(all1Jet->get_phi() - leadPhiGamma));
        if (std::isfinite(dphiLeadAbs) && dphiLeadAbs >= m_minBackToBack)
        {
          recoil1Pt  = all1Pt;
          recoil1Jet = all1Jet;
        }
      }

    // matching status category
    int status = 0;
    if (nPassPt == 0)                        status = 1; // NoJetPt
    else if (nPassEta == 0)                  status = 2; // NoJetEta
    else if (nPassDphi == 0 || !recoil1Jet)  status = 3; // NoBackToBack
    else                                     status = 4; // Matched

    const double maxDphiFill = (maxDphi >= 0.0 ? maxDphi : -0.01);

    // fill matching-QA histograms (radius-tagged; centrality-only suffix; pTgamma is axis)
    for (const auto& trigShort : activeTrig)
    {
      if (auto* hs = getOrBookMatchStatusVsPtGamma(trigShort, rKey, effCentIdx_M))
      { hs->Fill(leadPtGamma, status); bumpHistFill(trigShort, hs->GetName()); }

      if (auto* hm = getOrBookMatchMaxDphiVsPtGamma(trigShort, rKey, effCentIdx_M))
      { hm->Fill(leadPtGamma, maxDphiFill); bumpHistFill(trigShort, hm->GetName()); }

      if (auto* pn = getOrBookNRecoilJetsVsPtGamma(trigShort, rKey, effCentIdx_M))
      { pn->Fill(leadPtGamma, static_cast<double>(nPassDphi)); bumpHistFill(trigShort, pn->GetName()); }

      if (recoil1Jet)
      {
        const double dphiSel = std::fabs(TVector2::Phi_mpi_pi(recoil1Jet->get_phi() - leadPhiGamma));
        if (auto* hd = getOrBookMatchDphiVsPtGamma(trigShort, rKey, effCentIdx_M))
        { hd->Fill(leadPtGamma, dphiSel); bumpHistFill(trigShort, hd->GetName()); }

        const double isLeading = (all1Jet && recoil1Jet == all1Jet) ? 1.0 : 0.0;
        if (auto* hl = getOrBookRecoilIsLeadingVsPtGamma(trigShort, rKey, effCentIdx_M))
        { hl->Fill(leadPtGamma, isLeading); bumpHistFill(trigShort, hl->GetName()); }
      }
    }

    // ------------------------------------------------------------------
    // Inclusive γ–jet pairing histograms for unfolding (ATLAS-style)
    // ------------------------------------------------------------------
    for (std::size_t irj = 0; irj < recoJetsFid.size(); ++irj)
    {
      if (!recoJetsFidIsRecoil[irj]) continue;

      const Jet* rj = recoJetsFid[irj];
      if (!rj) continue;

      const double jpt  = rj->get_pt();
      const double jphi = rj->get_phi();
      if (!std::isfinite(jpt) || jpt <= 0.0) continue;
      if (!std::isfinite(jphi)) continue;

      const double xJr     = jpt / leadPtGamma;
      const double dphiAbs = std::fabs(TVector2::Phi_mpi_pi(jphi - leadPhiGamma));

      for (const auto& trigShort : activeTrig)
      {
        if (auto* h2 = getOrBookUnfoldRecoPtXJIncl(trigShort, rKey, effCentIdx_M))
        {
          h2->Fill(leadPtGamma, xJr);
          bumpHistFill(trigShort, h2->GetName());
        }

        // inclusive |Δphi(γ,jet)| for each recoil jet passing cuts
        if (auto* h2d = getOrBookUnfoldRecoPtDphiIncl(trigShort, rKey, effCentIdx_M))
        {
          h2d->Fill(leadPtGamma, dphiAbs);
          bumpHistFill(trigShort, h2d->GetName());
        }
      }
    }

    // ------------------------------------------------------------------
    // response matrix + truth distribution for unfolding
    // ------------------------------------------------------------------
    if (m_isSim && haveTruthPho && (tPt > 0.0))
    {
      fillUnfoldResponseMatrixAndTruthDistributions(activeTrig,
                                                    rKey,
                                                    effCentIdx_M,
                                                    leadPtGamma,
                                                    tPt,
                                                    tPhi,
                                                    recoJetsFid,
                                                    recoJetsFidIsRecoil,
                                                    recoil1Jet);
    }

    // ------------------------------------------------------------------
    // Physics-output fills (xJ, alpha, JES3) + Jet13/Profile3D — radius-tagged
    // ------------------------------------------------------------------
    if (recoil1Pt > 0.0 && recoil1Jet)
    {
      filledAnyRadius = true;

      const double xJ = recoil1Pt / leadPtGamma;

        // α definition (correct + robust):
        //   jet2 = highest-pT *fiducial* jet (pT+eta) excluding the selected recoil jet,
        //          WITH photon–jet overlap removal so we do NOT count the photon clustered as a jet.
        //
        //   Overlap removal:
        //     require ΔR(γ, jet2) > max(0.30, R)  (0.30 ~ your isolation cone; R is jet radius)
        //
        //   No Δφ requirement on jet2 (standard α radiation proxy).
        const Jet* jet2Jet = nullptr;
        double jet2Pt = 0.0;

        auto dR = [](double eta1, double phi1, double eta2, double phi2) -> double
        {
          const double dphi = TVector2::Phi_mpi_pi(phi1 - phi2);
          const double deta = eta1 - eta2;
          return std::sqrt(deta*deta + dphi*dphi);
        };

        const double minPhoJetDR = std::max(0.30, Rjet);

        for (const Jet* j2 : recoJetsFid) // recoJetsFid already pass pT + eta
        {
          if (!j2) continue;
          if (j2 == recoil1Jet) continue;

          const double j2pt  = j2->get_pt();
          const double j2eta = j2->get_eta();
          const double j2phi = j2->get_phi();

          if (!std::isfinite(j2pt) || !std::isfinite(j2eta) || !std::isfinite(j2phi)) continue;

          // Photon–jet overlap removal (prevents alpha from being dominated by a jet at the photon direction)
          const double drPho = dR(j2eta, j2phi, leadEtaGamma, leadPhiGamma);
          if (drPho <= minPhoJetDR) continue;

          if (j2pt > jet2Pt)
          {
            jet2Pt  = j2pt;
            jet2Jet = j2;
          }
      }

      const double alpha = (leadPtGamma > 0.0 ? (jet2Pt / leadPtGamma) : 0.0);

      if (Verbosity() >= 5)
      {
        LOG(5, CLR_GREEN,
            "      [lead pho#" << leadPhoIndex << "] rKey=" << rKey
            << " (R=" << std::fixed << std::setprecision(2) << Rjet
            << ", |eta|<" << std::fixed << std::setprecision(2) << jetEtaAbsMaxUse << ")"
            << " | xJ=" << std::fixed << std::setprecision(3) << xJ
            << " jet1Pt=" << recoil1Pt
            << " jet2Pt(fid,no#Delta#phi)=" << jet2Pt
            << " alpha=" << alpha
            << " (passed dphi=" << nPassDphi
            << ", passed pt=" << nPassPt << ")");
      }

      const int effCentIdx = (m_isAuAu ? centIdx : -1);

      // Jet-only QA for selected jets (jet2 uses alpha-definition above)
      fillSelectedJetQA(activeTrig, leadPtIdx, effCentIdx, rKey, recoil1Jet, jet2Jet);

      // Jet13 + Balance3 (radius-tagged)
      {
        const double jeta1 = recoil1Jet->get_eta();
        const double jphi1 = TVector2::Phi_mpi_pi(recoil1Jet->get_phi());

        for (const auto& trigShort : activeTrig)
        {
          if (auto* h3j = getOrBookJet13RecoilJet1(trigShort, rKey))
          { h3j->Fill(leadPtGamma, jeta1, jphi1); bumpHistFill(trigShort, h3j->GetName()); }

          if (auto* p3 = getOrBookBalance3(trigShort, rKey))
          { p3->Fill(leadPtGamma, jeta1, jphi1, xJ); bumpHistFill(trigShort, p3->GetName()); }
        }
      }

      for (const auto& trigShort : activeTrig)
      {
        if (auto* hx = getOrBookXJHist(trigShort, rKey, leadPtIdx, effCentIdx))
        { hx->Fill(xJ); bumpHistFill(trigShort, hx->GetName()); }

        if (auto* hj1 = getOrBookJet1PtHist(trigShort, rKey, leadPtIdx, effCentIdx))
        { hj1->Fill(recoil1Pt); bumpHistFill(trigShort, hj1->GetName()); }

        if (auto* hj2 = getOrBookJet2PtHist(trigShort, rKey, leadPtIdx, effCentIdx))
        { hj2->Fill(jet2Pt); bumpHistFill(trigShort, hj2->GetName()); }

        if (auto* ha = getOrBookAlphaHist(trigShort, rKey, leadPtIdx, effCentIdx))
        { ha->Fill(alpha); bumpHistFill(trigShort, ha->GetName()); }

        if (auto* h3x = getOrBookJES3_xJ_alphaHist(trigShort, rKey, effCentIdx))
        { h3x->Fill(leadPtGamma, xJ, alpha); bumpHistFill(trigShort, h3x->GetName()); }

        if (auto* h3j = getOrBookJES3_jet1Pt_alphaHist(trigShort, rKey, effCentIdx))
        { h3j->Fill(leadPtGamma, recoil1Pt, alpha); bumpHistFill(trigShort, h3j->GetName()); }
      }

      // -------------------- SIM ONLY: truth matching QA (radius-tagged) --------------------
      if (m_isSim && truth && haveTruthPho)
      {
        fillRecoTruthJES3MatchingQA(activeTrig,
                                    rKey,
                                    effCentIdx_M,
                                    leadPtGamma,
                                    xJ,
                                    alpha,
                                    tPt,
                                    tPhi,
                                    recoil1Jet);
      }
    }
    else
    {
      if (Verbosity() >= 4)
      {
        LOG(4, CLR_YELLOW,
            "      [lead pho#" << leadPhoIndex << "] rKey=" << rKey
            << " no recoil jet1 passes Δφ/pt/eta cuts"
            << " [passed dphi=" << nPassDphi
            << ", passed pt=" << nPassPt
            << ", minBackToBack=" << m_minBackToBack
            << ", minJetPt=" << m_minJetPt
            << ", |eta|<" << std::fixed << std::setprecision(2) << jetEtaAbsMaxUse << "]");
      }
    }
  } // radii loop

  return filledAnyRadius;
}




bool RecoilJets::runLeadIsoTightPhotonJetMatchingAndUnfolding(
    const std::vector<std::string>& activeTrig,
    const int effCentIdx_M,
    const int centIdx,
    const int leadPhoIndex,
    const int leadPtIdx,
    const double leadPtGamma,
    const double leadEtaGamma,
    const double leadPhiGamma,
    const bool haveTruthSigPho,
    const double tPtSig,
    const bool haveTruthPho,
    const double tPt,
    const double tPhi,
    PHG4TruthInfoContainer* truth)
{
  // This function only does unfolding photon bookkeeping + delegates jet logic per radius.

  // -------------------- Photon-only unfolding fills (N_gamma) --------------------
  for (const auto& trigShort : activeTrig)
  {
    // Measured reco photon spectrum (DATA + SIM)
    if (auto* hR = getOrBookUnfoldRecoPhoPtGamma(trigShort, effCentIdx_M))
    { hR->Fill(leadPtGamma); bumpHistFill(trigShort, hR->GetName()); }
  
    if (m_isSim)
    {
      // Reco photon fakes (selected reco photon not matching truth signal photon)
      if (!haveTruthPho)
      {
        if (auto* hRF = getOrBookUnfoldRecoPhoFakesPtGamma(trigShort, effCentIdx_M))
        { hRF->Fill(leadPtGamma); bumpHistFill(trigShort, hRF->GetName()); }
      }
  
      // Photon response (truth -> reco) only for matched truth signal photon
      if (haveTruthSigPho && haveTruthPho)
      {
        if (auto* hResp = getOrBookUnfoldResponsePhoPtGamma(trigShort, effCentIdx_M))
        { hResp->Fill(tPtSig, leadPtGamma); bumpHistFill(trigShort, hResp->GetName()); }
      }
  
      // Truth misses (truth signal exists but does not appear as the selected reco photon)
      if (haveTruthSigPho && !haveTruthPho)
      {
        if (auto* hTM = getOrBookUnfoldTruthPhoMissesPtGamma(trigShort, effCentIdx_M))
        { hTM->Fill(tPtSig); bumpHistFill(trigShort, hTM->GetName()); }
      }
    }
  }
  
  
  if (Verbosity() >= 5)
  {
    LOG(5, CLR_CYAN,
        "      [lead pho#" << leadPhoIndex << "] jet matching for all radii"
        << " | pT^γ=" << std::fixed << std::setprecision(2) << leadPtGamma
        << " eta^γ=" << std::fixed << std::setprecision(3) << leadEtaGamma
        << " phi^γ=" << std::fixed << std::setprecision(3) << leadPhiGamma);
  }
  
  bool filledAnyRadius = false;
  
    filledAnyRadius =
        runLeadIsoTightPhotonJetLoopAllRadii(activeTrig,
                                             effCentIdx_M,
                                             centIdx,
                                             leadPhoIndex,
                                             leadPtIdx,
                                             leadPtGamma,
                                             leadEtaGamma,
                                             leadPhiGamma,
                                             haveTruthPho,
                                             tPt,
                                             tPhi,
                                             truth);

  return filledAnyRadius;
}

void RecoilJets::fillPureIsolationQA(PHCompositeNode* topNode,
                                     const std::vector<std::string>& activeTrig,
                                     const PhotonClusterv1* pho,
                                     const RawCluster* rc,
                                     const int ptIdx,
                                     const int centIdx,
                                     const double pt_gamma)
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

  // ----------------------------------------------------------------
  //  accumulate negative-iso fractions once per reco photon candidate
  // (NOT per trigger), and BEFORE any SS-based cuts.
  // ----------------------------------------------------------------
  if (ptIdx >= 0 && ptIdx < static_cast<int>(m_nIsoBuilderByPt.size()))
  {
    if (std::isfinite(eiso_tot) && eiso_tot < 1e8)
    {
      ++m_nIsoBuilderByPt[ptIdx];
      if (eiso_tot < 0.0) ++m_nIsoBuilderNegByPt[ptIdx];
    }
  }

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

    // (B) component isolation histograms (EMCal / IHCal / OHCal)
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


void RecoilJets::fillTruthSigABCDLeakageCounters(PHCompositeNode* topNode,
                                                 const std::vector<std::string>& activeTrig,
                                                 const int centIdx)
{
  const int effCentIdx_sig = (m_isAuAu ? centIdx : -1);

  // Need HepMC event for truth-signal definition
  PHHepMCGenEventMap* hepmcmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  PHHepMCGenEvent*    hepmc    = nullptr;
  HepMC::GenEvent*    evt      = nullptr;

  if (hepmcmap)
  {
    hepmc = hepmcmap->get(0);
    if (!hepmc) hepmc = hepmcmap->get(1);
    if (!hepmc && !hepmcmap->empty()) hepmc = hepmcmap->begin()->second;
    if (hepmc) evt = hepmc->getEvent();
  }

  if (!evt)
  {
    if (Verbosity() >= 4)
    {
      LOG(4, CLR_YELLOW,
          "    [sigABCD] SIM: HepMC event missing → cannot fill truth-signal leakage counters");
    }
    return;
  }

  // CaloRawClusterEval enforces the "BEST-MATCHED truth primary particle" condition
  CaloRawClusterEval clustereval(topNode, "CEMC");
  clustereval.set_usetowerinfo(true);
  clustereval.next_event(topNode);

  // Fallback if CLUSTERINFO_* is not available in the DST
  if (!clustereval.has_reduced_node_pointers())
  {
    if (Verbosity() >= 6)
    {
      LOG(6, CLR_YELLOW,
          "    [sigABCD] SIM: CaloRawClusterEval reduced nodes missing with towerinfo=true → retry with towerinfo=false");
    }

    clustereval.set_usetowerinfo(false);
    clustereval.next_event(topNode);
  }

  if (!clustereval.has_reduced_node_pointers())
  {
    if (Verbosity() >= 4)
    {
      LOG(4, CLR_YELLOW,
          "    [sigABCD] SIM: CaloRawClusterEval missing required nodes → cannot truth↔reco match photons → skip leakage counters");
    }
    return;
  }

  int nTruthSig        = 0;
  int nTruthSigMatched = 0;

  for (auto it = evt->particles_begin(); it != evt->particles_end(); ++it)
  {
    const HepMC::GenParticle* p = *it;
    if (!p) continue;

    // Truth isolation is computed inside isTruthPromptIsolatedSignalPhoton(...)
    // for prompt (direct/frag) final-state photons in acceptance.
    double isoEtTruth   = std::numeric_limits<double>::quiet_NaN();
    const bool passTruthIso = isTruthPromptIsolatedSignalPhoton(evt, p, isoEtTruth);

    // ---------------------------
    //  Truth isolation QA (SIM)
    // ---------------------------
    // Fill the isoEtTruth distribution for ALL prompt truth photons where isoEtTruth was computed
    // (this includes both PASS and FAIL of the iso cut).
    if (std::isfinite(isoEtTruth))
    {
      for (const auto& trigShort : activeTrig)
      {
        // Distribution
        if (auto* hIsoT = getOrBookTruthIsoHist(trigShort, "h_EisoTruth", 200, 0.0, 50.0))
        {
          hIsoT->Fill(isoEtTruth);
          bumpHistFill(trigShort, hIsoT->GetName());
        }

        // PASS/FAIL decision (bin1=PASS, bin2=FAIL)
        if (auto* hDecT = getOrBookTruthIsoDecisionHist(trigShort, "h_EisoTruthDecision"))
        {
          hDecT->Fill(passTruthIso ? 1 : 2);
          bumpHistFill(trigShort, hDecT->GetName());
        }
      }
    }

    // Keep your original leakage logic: only proceed for isolated truth-signal photons
    if (!passTruthIso) continue;
    ++nTruthSig;

    const RawCluster* recoMatch = nullptr;
    double rPt = 0.0, rEta = 0.0, rPhi = 0.0, drBest = 1e9;
    float  eBest = -1.0f;

    if (!findRecoPhotonMatchedToTruthSignal(evt, p, clustereval,
                                           recoMatch, rPt, rEta, rPhi, drBest, eBest))
    {
      if (Verbosity() >= 8)
      {
        LOG(8, CLR_YELLOW,
            "    [sigABCD] truth barcode=" << p->barcode()
            << " passTruthIso=1 but no reco match found");
      }
      continue;
    }
    ++nTruthSigMatched;

    const auto* recoPho = dynamic_cast<const PhotonClusterv1*>(recoMatch);
    if (!recoPho)
    {
      if (Verbosity() >= 8)
      {
        LOG(8, CLR_YELLOW,
            "    [sigABCD] truth barcode=" << p->barcode()
            << " matched reco cluster is not PhotonClusterv1 → skip");
      }
      continue;
    }

    const int ptIdx_sig = findPtBin(rPt);
    if (ptIdx_sig < 0)
    {
      if (Verbosity() >= 8)
      {
        LOG(8, CLR_YELLOW,
            "    [sigABCD] truth barcode=" << p->barcode()
            << " matched reco pT=" << std::fixed << std::setprecision(2) << rPt
            << " outside configured m_gammaPtBins → skip");
      }
      continue;
    }

    const double eiso_et = eiso(recoMatch, topNode);
    if (!std::isfinite(eiso_et) || eiso_et > 1e8)
    {
      if (Verbosity() >= 8)
      {
        LOG(8, CLR_YELLOW,
            "    [sigABCD] truth barcode=" << p->barcode()
            << " reco pT=" << std::fixed << std::setprecision(2) << rPt
            << " has invalid Eiso=" << eiso_et << " → skip");
      }
      continue;
    }

    // reco isolation distribution for TRUTH-ISOLATED signal photons (direct+frag),
    // independent of reco shower-shape classification.
    {
      const std::string slice_sig = suffixForBins(ptIdx_sig, effCentIdx_sig);

      for (const auto& trigShort : activeTrig)
      {
        if (auto* hIso = getOrBookIsoPartHist(trigShort,
                                             "h_EisoReco_truthSigMatched",
                                             "E_{T}^{iso,reco} [GeV] (truth iso signal match)",
                                             ptIdx_sig, effCentIdx_sig))
        {
          hIso->Fill(eiso_et);
          bumpHistFill(trigShort, std::string("h_EisoReco_truthSigMatched") + slice_sig);
        }
      }
    }

    // Reco-side ABCD classification (PPG12-equivalent)
    const SSVars   v   = makeSSFromPhoton(recoPho, rPt);
    const TightTag tag = classifyPhotonTightness(v);

    // Exclude: preselection fail and Neither(1 fail), consistent with your ABCD logic
    if (!(tag == TightTag::kTight || tag == TightTag::kNonTight))
    {
      if (Verbosity() >= 8)
      {
        LOG(8, CLR_YELLOW,
            "    [sigABCD] truth barcode=" << p->barcode()
            << " reco pT=" << std::fixed << std::setprecision(2) << rPt
            << " fails preselection/Neither → tag=" << tightTagName(tag) << " → skip");
      }
      continue;
    }

    const double thrIso    = m_isoA + m_isoB * rPt;
    const double thrNonIso = thrIso + m_isoGap;

    const bool iso    = (eiso_et < thrIso);
    const bool nonIso = (eiso_et > thrNonIso);

    // GAP excluded
    if (!iso && !nonIso)
    {
      if (Verbosity() >= 8)
      {
        LOG(8, CLR_YELLOW,
            "    [sigABCD] truth barcode=" << p->barcode()
            << " reco pT=" << std::fixed << std::setprecision(2) << rPt
            << " lies in ISO GAP (Eiso=" << std::fixed << std::setprecision(3) << eiso_et
            << ", thrIso=" << thrIso
            << ", thrNonIso=" << thrNonIso << ") → skip");
      }
      continue;
    }

    int regBin = 0;
    if      (iso    && tag == TightTag::kTight)    regBin = 1; // A
    else if (nonIso && tag == TightTag::kTight)    regBin = 2; // B
    else if (iso    && tag == TightTag::kNonTight) regBin = 3; // C
    else if (nonIso && tag == TightTag::kNonTight) regBin = 4; // D

    if (regBin <= 0) continue;

    for (const auto& trigShort : activeTrig)
    {
      if (auto* h = getOrBookSigABCDLeakageHist(trigShort, ptIdx_sig, effCentIdx_sig))
      {
        h->Fill(regBin);
        bumpHistFill(trigShort, h->GetName());
      }
    }

    if (Verbosity() >= 7)
    {
      LOG(7, CLR_CYAN,
          "    [sigABCD] truth barcode=" << p->barcode()
          << " reco pT=" << std::fixed << std::setprecision(2) << rPt
          << " Eiso=" << std::fixed << std::setprecision(3) << eiso_et
          << " thrIso=" << std::fixed << std::setprecision(3) << thrIso
          << " thrNonIso=" << std::fixed << std::setprecision(3) << thrNonIso
          << " tag=" << tightTagName(tag)
          << " → region=" << regBin
          << " (drBest=" << std::fixed << std::setprecision(4) << drBest
          << ", eContrib=" << eBest
          << ", isoEtTruth=" << std::fixed << std::setprecision(3) << isoEtTruth << ")");
    }
  } // end truth particle loop

  if (Verbosity() >= 5)
  {
    LOG(5, CLR_BLUE,
        "    [sigABCD] summary: truthSig=" << nTruthSig
        << " truthSigMatchedReco=" << nTruthSigMatched
        << " (ABCD filled only if reco passes preselection AND is not GAP AND not Neither)");
  }
}





// Centralized candidate processing for the event
void RecoilJets::processCandidates(PHCompositeNode* topNode,
                                   const std::vector<std::string>& activeTrig)
{
  // -------- guards & context ------------------------------------------------
  if (!topNode)
  {
    LOG(0, CLR_YELLOW, "  [processCandidates] topNode == nullptr – aborting this event :(");
    return;
  }
  if (activeTrig.empty())
  {
    LOG(0, CLR_YELLOW, "  [processCandidates] activeTrig is EMPTY – nothing to fill this event :(");
    return;
  }

  // Centrality bin index (Au+Au only, else -1)
  const int centIdx = (m_isAuAu ? findCentBin(m_centBin) : -1);

  if (Verbosity() >= 4)
  {
      std::ostringstream os;
      os << "  [processCandidates] dataset=" << (m_isAuAu ? "Au+Au" : "p+p")
         << "  centIdx=" << centIdx
         << "  m_centBin=" << m_centBin
         << "  nTriggers=" << activeTrig.size()
         << "  triggers={";

      for (std::size_t i = 0; i < activeTrig.size(); ++i)
      {
        if (i) os << ", ";
        os << activeTrig[i];
      }
      os << "}";

      LOG(4, CLR_BLUE, os.str());
    }

    // If neither photons nor clusters are present, we cannot do anything
    if (!m_photons && !m_clus)
    {
      LOG(2, CLR_YELLOW, "  [processCandidates] PHOTONCLUSTER_CEMC and CLUSTERINFO_CEMC both MISSING – nothing to process - sus");
      return;
    }

    // ==========================================================================
    //  Truth-signal → reco ABCD leakage counters
    //
    // For each truth "signal photon" (prompt + truth-iso) that has a matched reco
    // PhotonClusterBuilder cluster, classify the reco cluster into PPG12 ABCD
    // regions using the SAME reco definitions as your purity logic:
    //   - preselection must pass
    //   - GAP excluded
    //   - Neither(1 fail) excluded
    //
    // Output histograms (under /SIM/):
    //   h_sigABCD_MC_pT_lo_hi[_cent_lo_hi]  (TH1I, bins: 1=A, 2=B, 3=C, 4=D)
    // ==========================================================================
    if (m_isSim && m_photons)
    {
      fillTruthSigABCDLeakageCounters(topNode, activeTrig, centIdx);
    }

    // ==========================================================================
    // (SIM ONLY): PURE truth xJgamma distribution (NO reco gating)
    //
    // Books/Fills (per trigger, per radius, per centrality-suffix):
    //   h_JES3TruthPure_pT_xJ_alpha_<rKey><centSuffix>
    //
    // Truth photon:
    //   uses your encoded truth-signal definition:
    //     isTruthPromptIsolatedSignalPhoton(evt, p, isoEt)
    //
    // Truth recoil jets (per rKey truth jet container):
    //   - pT > m_minJetPt
    //   - |eta| < (1.1 - R)
    //   - |Δphi(truth gamma, truth jet)| >= m_minBackToBack
    //
    // Then fill:
    //   (tPt, xJt=tj1Pt/tPt, aT=tj2Pt/tPt)
    // ==========================================================================
    if (m_isSim)
    {
      const int effCentIdx_truth = (m_isAuAu ? centIdx : -1);

      // Need HepMC event for truth-signal definition
      PHHepMCGenEventMap* hepmcmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
      PHHepMCGenEvent*    hepmc    = nullptr;
      HepMC::GenEvent*    evt      = nullptr;

      if (hepmcmap)
      {
        hepmc = hepmcmap->get(0);
        if (!hepmc) hepmc = hepmcmap->get(1);
        if (!hepmc && !hepmcmap->empty()) hepmc = hepmcmap->begin()->second;
        if (hepmc) evt = hepmc->getEvent();
      }

      if (!evt)
      {
        if (Verbosity() >= 4)
        {
          LOG(4, CLR_YELLOW,
              "    [truthXJgamma] SIM: HepMC event missing → cannot fill h_JES3TruthPure_pT_xJ_alpha");
        }
      }
      else
      {
        // Pick the event-leading truth signal photon (highest pT)
        bool   haveTruthSigPho = false;
        double tPt  = -1.0;
        double tPhi = 0.0;

        for (auto it = evt->particles_begin(); it != evt->particles_end(); ++it)
        {
          const HepMC::GenParticle* p = *it;
          if (!p) continue;

          double isoEt = 0.0;
          if (!isTruthPromptIsolatedSignalPhoton(evt, p, isoEt)) continue;

          const double pt  = std::hypot(p->momentum().px(), p->momentum().py());
          const double eta = p->momentum().pseudoRapidity();
          const double phi = TVector2::Phi_mpi_pi(p->momentum().phi());

          if (!std::isfinite(pt) || !std::isfinite(eta) || !std::isfinite(phi) || pt <= 0.0) continue;

          if (!haveTruthSigPho || pt > tPt)
          {
            haveTruthSigPho = true;
            tPt  = pt;
            tPhi = phi;
          }
        }

        if (!haveTruthSigPho || tPt <= 0.0)
        {
          if (Verbosity() >= 5)
          {
            LOG(5, CLR_YELLOW,
                "    [truthXJgamma] no truth signal photon found → skip pure truth xJgamma fills");
          }
        }
        else
        {
          // Fill per jet radius (truth jets are radius-tagged)
          for (const auto& kvT : m_truthJetsByRKey)
          {
            const std::string rKey = kvT.first;

            JetContainer* truthJets = kvT.second;

            if (!truthJets)
            {
              if (Verbosity() >= 5)
              {
                LOG(5, CLR_YELLOW,
                    "    [truthXJgamma] rKey=" << rKey << " truth jet container missing → skip");
              }
              continue;
            }

            const double etaMaxTruth = jetEtaAbsMaxForRKey(rKey);

            double tj1Pt = -1.0;
            const Jet* tj1 = nullptr;

            // 1) truth recoil jet1: highest-pT fiducial jet back-to-back to truth gamma
            for (const Jet* tj : *truthJets)
            {
              if (!tj) continue;

              const double ptj  = tj->get_pt();
              const double etaj = tj->get_eta();
              const double phij = tj->get_phi();

              if (!std::isfinite(ptj) || !std::isfinite(etaj) || !std::isfinite(phij)) continue;
              if (ptj < m_minJetPt) continue;
              if (std::fabs(etaj) >= etaMaxTruth) continue;

              const double dphiAbs = std::fabs(TVector2::Phi_mpi_pi(phij - tPhi));
              if (dphiAbs < m_minBackToBack) continue;

              if (ptj > tj1Pt)
              {
                tj1Pt = ptj;
                tj1   = tj;
              }
            }

            if (!tj1 || tj1Pt <= 0.0)
            {
              if (Verbosity() >= 6)
              {
                LOG(6, CLR_YELLOW,
                    "    [truthXJgamma] rKey=" << rKey << " no truth recoil jet1 found → skip");
              }
              continue;
            }

            // 2) truth jet2 for α: highest-pT fiducial truth jet excluding tj1 (NO Δφ requirement)
            double tj2Pt = -1.0;
            for (const Jet* tj : *truthJets)
            {
              if (!tj) continue;

              const double ptj  = tj->get_pt();
              const double etaj = tj->get_eta();
              const double phij = tj->get_phi();

              if (!std::isfinite(ptj) || !std::isfinite(etaj) || !std::isfinite(phij)) continue;
              if (ptj < m_minJetPt) continue;
              if (std::fabs(etaj) >= etaMaxTruth) continue;
              if (tj == tj1) continue;

              if (ptj > tj2Pt) tj2Pt = ptj;
            }

            const double xJt = tj1Pt / tPt;
            const double aT  = (tj2Pt > 0.0 ? (tj2Pt / tPt) : 0.0);

            for (const auto& trigShort : activeTrig)
            {
              if (auto* ht3_pure = getOrBookJES3TruthPure_xJ_alphaHist(trigShort, rKey, effCentIdx_truth))
              {
                ht3_pure->Fill(tPt, xJt, aT);
                bumpHistFill(trigShort, ht3_pure->GetName());
              }
            }

            if (Verbosity() >= 7)
            {
              LOG(7, CLR_CYAN,
                  "    [truthXJgamma] rKey=" << rKey
                  << " tPt=" << std::fixed << std::setprecision(2) << tPt
                  << " tj1Pt=" << std::fixed << std::setprecision(2) << tj1Pt
                  << " tj2Pt(fid,no#Delta#phi)=" << std::fixed << std::setprecision(2) << (tj2Pt > 0.0 ? tj2Pt : 0.0)
                  << " → xJt=" << std::fixed << std::setprecision(3) << xJt
                  << " aT="  << std::fixed << std::setprecision(3) << aT);
            }
          } // end rKey loop
        }
      }
    }

    // ========================= PHOTON path ========================
    if (m_photons)
    {
    RawClusterContainer::ConstRange prange = m_photons->getClusters();

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
      double leadEtaGamma     = 0.0;
      double leadPhiGamma     = 0.0;

      //  leading reco photon candidate in pT^gamma within analysis pT bins,
      // WITHOUT requiring iso or tightness (used for jet-level cutflow QA vs pT^gamma)
      bool   haveLeadAnyPho  = false;
      int    leadAnyPtIdx    = -1;
      double leadAnyPtGamma  = -1.0;
      double leadAnyPhiGamma = 0.0;

      // Keep a pointer to the actual leading reco photon cluster (needed for SIM truth↔reco association)
      const RawCluster* leadRc = nullptr;
        
      // Event-level diagnostic: number of reco photon candidates that pass the SAME
      // iso∧tight gate used for the event-leading photon selection in this event.
      int nIsoTightPhoCand = 0;

      int iPho = 0;
      for (auto pit = prange.first; pit != prange.second; ++pit, ++iPho)
      {
        // Concrete type: PhotonClusterv1
        const auto* pho = dynamic_cast<const PhotonClusterv1*>(pit->second);
        if (!pho)
        {
              LOG(5, CLR_YELLOW, "      [pho#" << iPho << "] pointer is not PhotonClusterv1 – skipping");
              continue;
        }

        // Upcast through inheritance (no RTTI, no second cast)
        const RawCluster* rc = pho;



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
                << " → skipping candidate.");
            continue;
        }

           // -------- reco pT^gamma floor (PPG12 reco-matching requirement) ----
           // PPG12 requirement you stated: pT_reco > 5 GeV  (and for photons pT == ET).
           constexpr double kMinPtGamma = 5.0; // GeV

           if (pt_gamma < kMinPtGamma)
           {
             ++nSkipEtBin; ++m_bk.pho_early_E;

             if (Verbosity() >= 6)
             {
               const int dbgPtIdx = findPtBin(pt_gamma);

               std::ostringstream os;
               os << "      [pho#" << iPho << "] early reject (pT^#gamma floor)"
                  << " | pT="  << std::fixed << std::setprecision(2) << pt_gamma
                  << " < "     << kMinPtGamma
                  << " | E="   << std::fixed << std::setprecision(2) << energy
                  << " | eta=" << std::fixed << std::setprecision(3) << eta
                  << " | phi=" << std::fixed << std::setprecision(3) << phi
                  << " | pTbinIdx=" << dbgPtIdx;
               LOG(6, CLR_YELLOW, os.str());
             }
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

        //  choose event-leading reco photon candidate in pT^gamma (within analysis bins),
        // WITHOUT iso/tight requirement.
        if (!haveLeadAnyPho || pt_gamma > leadAnyPtGamma)
        {
                haveLeadAnyPho  = true;
                leadAnyPtIdx    = ptIdx;
                leadAnyPtGamma  = pt_gamma;
                leadAnyPhiGamma = phi_gamma;
        }

        fillPureIsolationQA(topNode, activeTrig, pho, rc, ptIdx, centIdx, pt_gamma);

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

        // record the 2×2 (iso, tight) category + SS variable hists
        // This also fills h_Eiso once and prints a detailed decision line.
        for (const auto& trigShort : activeTrig)
        {
          fillIsoSSTagCounters(trigShort, rc, v, pt_gamma, centIdx, topNode);
        }

           //  xJ usable gate (preselection already passed above):
           //    - for ALL iso∧tight photons: fill the reco photon 3D baseline TH3
           //    - for jet matching: defer to event-leading iso∧tight photon
          if (!(iso && tightTag == TightTag::kTight))
          {
            if (!iso) ++nNotIso;
            if (tightTag != TightTag::kTight) ++nNotTight;
            if (Verbosity() >= 6)
              LOG(6, CLR_BLUE, "      [pho#" << iPho << "] NOT used for xJ (iso=" << iso
                                             << ", tightTag=" << tightTagName(tightTag) << ")");
            continue;
          }

          // Count ALL iso∧tight photon candidates in this event (photon-side ambiguity diagnostic)
          ++nIsoTightPhoCand;

          for (const auto& trigShort : activeTrig)
           {
             if (auto* h3 = getOrBookPho3TightIso(trigShort))
             {
               h3->Fill(pt_gamma, eta, TVector2::Phi_mpi_pi(phi_gamma));
               bumpHistFill(trigShort, h3->GetName());
             }
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
             leadEtaGamma     = eta;
             leadPhiGamma     = phi_gamma;

             // Store the actual cluster pointer for later truth↔reco association
             leadRc           = rc;

             if (Verbosity() >= 6)
               LOG(6, CLR_GREEN, "      [pho#" << iPho << "] marked as event-leading iso∧tight photon for jet matching (pT="
                                               << std::fixed << std::setprecision(2) << leadPtGamma << ")");
           }

            // Defer jet matching (prevents double-filling JES TH3s if multiple photons pass).
            continue;
        } // photon loop
        
        // ------------------------------------------------------------------
        // Event-level photon-side diagnostic:
        //   N = number of iso∧tight photon candidates in this event.
        // Fill once per event (only when a leading iso∧tight photon exists),
        // binned by the selected leading-photon pT bin (leadPtIdx).
        // ------------------------------------------------------------------
        if (haveLeadIsoTight && leadPtIdx >= 0)
        {
          for (const auto& trigShort : activeTrig)
          {
            if (auto* hN = getOrBookNIsoTightPhoCandHist(trigShort, leadPtIdx, centIdx))
            {
              hN->Fill(nIsoTightPhoCand);
              bumpHistFill(trigShort, hN->GetName());
            }
          }
        }

        // ------------------------------------------------------------------
        // Jet matching + JES fills ONCE per event using the event-leading
        // iso∧tight photon (by pT^gamma), but run identically for ALL radii
        // in kJetRadii (r02 + r04) in the SAME loop.
        // ------------------------------------------------------------------

        // Centrality index used for unfolding/matching-QA suffixing (centrality-only suffix in Au+Au)
        const int effCentIdx_M = (m_isAuAu ? centIdx : -1);

        // -------------------- SIM: define the event-leading truth signal photon (for N_gamma unfolding) --------------------
        HepMC::GenEvent* evtHepMC = nullptr;

        bool   haveTruthSigPho = false;
        const HepMC::GenParticle* truthSigPho = nullptr;
        double tPtSig  = 0.0;
        double tEtaSig = 0.0;
        double tPhiSig = 0.0;
        double tIsoEtSig = 0.0;

        if (m_isSim)
        {
          // Grab HepMC event (required for prompt/direct/fragmentation classification)
          PHHepMCGenEventMap* hepmcmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
          PHHepMCGenEvent*    hepmc    = nullptr;

          if (hepmcmap)
          {
            // Try common keys then fallback to first entry
            hepmc = hepmcmap->get(0);
            if (!hepmc) hepmc = hepmcmap->get(1);
            if (!hepmc && !hepmcmap->empty()) hepmc = hepmcmap->begin()->second;
            if (hepmc) evtHepMC = hepmc->getEvent();
          }

          if (evtHepMC)
          {
            for (auto it = evtHepMC->particles_begin(); it != evtHepMC->particles_end(); ++it)
            {
              const HepMC::GenParticle* p = *it;
              if (!p) continue;

              double isoEt = 0.0;
              if (!isTruthPromptIsolatedSignalPhoton(evtHepMC, p, isoEt)) continue;

              const double pt  = std::hypot(p->momentum().px(), p->momentum().py());
              const double eta = p->momentum().pseudoRapidity();
              const double phi = p->momentum().phi();
              if (!std::isfinite(pt) || !std::isfinite(eta) || !std::isfinite(phi) || pt <= 0.0) continue;

              if (!haveTruthSigPho || pt > tPtSig)
              {
                haveTruthSigPho = true;
                truthSigPho = p;
                tPtSig  = pt;
                tEtaSig = eta;
                tPhiSig = phi;
                tIsoEtSig = isoEt;
              }
            }
          }
          else
          {
            if (Verbosity() >= 5)
              LOG(5, CLR_YELLOW,
                  "      [truthQA] PHHepMCGenEventMap/HepMC event missing → skip truth photon unfolding helpers");
          }

          // Fill truth photon spectrum once per event (SIM)
          if (haveTruthSigPho)
          {
            for (const auto& trigShort : activeTrig)
            {
              if (auto* hT = getOrBookUnfoldTruthPhoPtGamma(trigShort, effCentIdx_M))
              { hT->Fill(tPtSig); bumpHistFill(trigShort, hT->GetName()); }

              // If no reco leading iso∧tight photon exists, this truth photon is a MISS for N_gamma
              if (!haveLeadIsoTight)
              {
                if (auto* hTM = getOrBookUnfoldTruthPhoMissesPtGamma(trigShort, effCentIdx_M))
                { hTM->Fill(tPtSig); bumpHistFill(trigShort, hTM->GetName()); }
              }
            }
          }
        }

        if (haveLeadIsoTight)
        {
          // ΔR helper (used for SIM truth matching)
          auto dR = [](double eta1, double phi1, double eta2, double phi2) -> double
          {
            const double dphi = TVector2::Phi_mpi_pi(phi1 - phi2);
            const double deta = eta1 - eta2;
            return std::sqrt(deta*deta + dphi*dphi);
          };

          // -------------------- SIM: does the leading reco iso∧tight photon correspond to the event-leading truth signal photon? --------------------
          bool   haveTruthPho = false;
          double tPt  = 0.0, tEta = 0.0, tPhi = 0.0;

          PHG4TruthInfoContainer* truth = nullptr;

          if (m_isSim)
          {
            truth = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
            if (!truth)
            {
              LOG(4, CLR_YELLOW, "      [truthQA] G4TruthInfo missing → skip truth matching QA (all radii)");
            }
            else if (!evtHepMC || !haveTruthSigPho || !truthSigPho)
            {
              LOG(4, CLR_YELLOW,
                  "      [truthQA] no truth signal photon / HepMC event missing → skip truth photon association to lead reco photon");
            }
            else if (!leadRc)
            {
              LOG(4, CLR_YELLOW,
                  "      [truthQA] leadRc is nullptr (unexpected) → skip truth photon matching");
            }
            else
            {
              // CaloRawClusterEval enforces the "BEST-MATCHED truth primary particle" condition
              CaloRawClusterEval clustereval(topNode, "CEMC");
              clustereval.set_usetowerinfo(true);
              clustereval.next_event(topNode);

              // Fallback if CLUSTERINFO_* is not available in the DST
              if (!clustereval.has_reduced_node_pointers())
              {
                clustereval.set_usetowerinfo(false);
                clustereval.next_event(topNode);
              }

              if (!clustereval.has_reduced_node_pointers())
              {
                LOG(4, CLR_YELLOW,
                    "      [truthQA] CaloRawClusterEval missing required nodes (clusters/towers/truth) → skip truth↔reco photon association");
              }
              else
              {
                const double kPhoMatchDR = m_phoMatchDRMax;

                const RawCluster* recoMatch = nullptr;
                double rPt = 0.0, rEta = 0.0, rPhi = 0.0, drBest = 1e9;
                float  eBest = -1.0f;

                if (findRecoPhotonMatchedToTruthSignal(evtHepMC, truthSigPho, clustereval,
                                                      recoMatch, rPt, rEta, rPhi, drBest, eBest))
                {
                  // Require: the truth signal photon's matched reco is the event-leading iso∧tight photon
                  if (recoMatch == leadRc)
                  {
                    const double drPho = dR(leadEtaGamma, leadPhiGamma, tEtaSig, tPhiSig);
                    if (drPho < kPhoMatchDR)
                    {
                      haveTruthPho = true;
                      tPt  = tPtSig;
                      tEta = tEtaSig;
                      tPhi = tPhiSig;

                      if (Verbosity() >= 5)
                      {
                        LOG(5, CLR_GREEN,
                            "      [truthQA] lead reco photon is CaloRawClusterEval-best-matched to the event-leading truth signal photon:"
                            << " barcode=" << truthSigPho->barcode()
                            << " pT^truth=" << std::fixed << std::setprecision(2) << tPt
                            << " eta^truth=" << std::fixed << std::setprecision(3) << tEta
                            << " phi^truth=" << std::fixed << std::setprecision(3) << tPhi
                            << " ETiso^truth=" << std::fixed << std::setprecision(3) << tIsoEtSig
                            << " ΔR=" << std::fixed << std::setprecision(4) << drPho
                            << " Edep(best)=" << eBest);
                      }
                    }
                  }
                }

                if (!haveTruthPho && Verbosity() >= 5)
                {
                  LOG(5, CLR_YELLOW,
                      "      [truthQA] event-leading truth signal photon does NOT match the event-leading reco iso∧tight photon");
                }
              }
            }
          }
          const bool filledAnyRadius =
              runLeadIsoTightPhotonJetMatchingAndUnfolding(activeTrig,
                                                           effCentIdx_M,
                                                           centIdx,
                                                           leadPhoIndex,
                                                           leadPtIdx,
                                                           leadPtGamma,
                                                           leadEtaGamma,
                                                           leadPhiGamma,
                                                           haveTruthSigPho,
                                                           tPtSig,
                                                           haveTruthPho,
                                                           tPt,
                                                           tPhi,
                                                           truth);

            // Preserve your old "used" semantics: count once per event if any radius filled
          if (filledAnyRadius) ++nUsed;

        }
        else
        {
          if (Verbosity() >= 5)
            LOG(5, CLR_BLUE, "      [processCandidates] no iso∧tight photon in this event → no jet matching / JES fills");
        }

        // ------------------------------------------------------------------
        //  jet-by-jet cutflow status vs pT^gamma (NOT event-level).
        // Uses the event-leading reco photon candidate in pT^gamma WITHOUT requiring iso/tight.
        //
        // Fills: h_jetcutflow_status_vs_pTgamma_<rKey><centSuffix>
        // Y bins:
        //   1 = FailJetPt
        //   2 = FailJetEta
        //   3 = FailBackToBack
        //   4 = PassAll
        // ------------------------------------------------------------------
        if (haveLeadAnyPho && leadAnyPtIdx >= 0)
        {

          for (const auto& kvJ : m_jets)
          {
            const std::string& rKey = kvJ.first;

            JetContainer* jets = kvJ.second;
            if (!jets) continue;

            const double jetEtaAbsMaxUse = jetEtaAbsMaxForRKey(rKey);

            for (const Jet* j : *jets)
            {
              if (!j) continue;

              const double jpt  = j->get_pt();
              const double jeta = j->get_eta();
              const double jphi = j->get_phi();

              if (!std::isfinite(jpt) || !std::isfinite(jeta) || !std::isfinite(jphi)) continue;

              int jStatus = 0;
              if (jpt < m_minJetPt) jStatus = 1;
              else if (std::fabs(jeta) >= jetEtaAbsMaxUse) jStatus = 2;
              else
              {
                const double dphiAbs = std::fabs(TVector2::Phi_mpi_pi(jphi - leadAnyPhiGamma));
                jStatus = (dphiAbs >= m_minBackToBack) ? 4 : 3;
              }

              for (const auto& trigShort : activeTrig)
              {
                if (auto* hj = getOrBookJetCutflowStatusVsPtGamma(trigShort, rKey, effCentIdx_M))
                { hj->Fill(leadAnyPtGamma, jStatus); bumpHistFill(trigShort, hj->GetName()); }
              }
            }
          }
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
  // Require PhotonClusterBuilder products for ALL photon kinematics.
  LOG(0, CLR_YELLOW,
        "  [processCandidates] PHOTONCLUSTER_CEMC is MISSING for this event. ");

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

    const bool pass_e11e33 = (v.e11_over_e33 < m_phoid_pre_e11e33_max);
    const bool pass_et1    = in_open_interval(v.et1, m_phoid_pre_et1_min, m_phoid_pre_et1_max);
    const bool pass_e32e35 = in_open_interval(v.e32_over_e35, m_phoid_pre_e32e35_min, m_phoid_pre_e32e35_max);
    const bool pass_weta   = (v.weta_cogx < m_phoid_pre_weta_max);

    if (Verbosity() >= 5)
    {
      LOG(5, CLR_BLUE,
          "  [passesPhotonPreselection] "
          << "weta=" << v.weta_cogx << " (<" << m_phoid_pre_weta_max << ") → " << pass_weta
          << " | et1=" << v.et1 << " ∈ (" << m_phoid_pre_et1_min << "," << m_phoid_pre_et1_max << ") → " << pass_et1
          << " | e11/e33=" << v.e11_over_e33 << " (<" << m_phoid_pre_e11e33_max << ") → " << pass_e11e33
          << " | e32/e35=" << v.e32_over_e35 << " ∈ (" << m_phoid_pre_e32e35_min << "," << m_phoid_pre_e32e35_max << ") → " << pass_e32e35);
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
  const double w_hi = m_phoid_tight_w_hi_intercept + m_phoid_tight_w_hi_slope * v.pt_gamma;
  if (!std::isfinite(w_hi) || !std::isfinite(v.pt_gamma) || v.pt_gamma <= 0.0)
  {
        LOG(2, CLR_YELLOW, "  [classifyPhotonTightness] non-finite tight_w_hi for pT^γ=" << v.pt_gamma
                            << " – treating as failure");
        return TightTag::kNeither;
  }

  const bool pass_weta   = in_open_interval(v.weta_cogx,  m_phoid_tight_w_lo, w_hi);
  const bool pass_wphi   = in_open_interval(v.wphi_cogx,  m_phoid_tight_w_lo, w_hi);
  const bool pass_e11e33 = in_open_interval(v.e11_over_e33, m_phoid_tight_e11e33_min, m_phoid_tight_e11e33_max);
  const bool pass_et1    = in_open_interval(v.et1,          m_phoid_tight_et1_min,    m_phoid_tight_et1_max);
  const bool pass_e32e35 = in_open_interval(v.e32_over_e35, m_phoid_tight_e32e35_min, m_phoid_tight_e32e35_max);

  if (Verbosity() >= 5)
  {
      LOG(5, CLR_BLUE,
          "  [classifyPhotonTightness] pT^γ=" << v.pt_gamma
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

  // PPG12 isolation definition is fixed to ΔR = 0.3.
  // Ignore user-provided coneR for "exact PPG12" operation.
  constexpr double kPPG12ConeR = 0.3;

    if (Verbosity() >= 2 && std::fabs(coneR - kPPG12ConeR) > 1e-9)
    {
      LOG(2, CLR_YELLOW,
          "  [setIsolationWP] Requested coneR=" << coneR
          << " but PPG12 requires coneR=0.3. Forcing coneR=0.3.");
    }

  m_isoConeR  = kPPG12ConeR;
  m_isoTowMin = std::max(0.0, towerMin);

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

  return (eiso_val > thr);
}


// -----------------------------------------------------------------------------
// Unified truth-MC signal definition: isolated prompt photon (CaloAna-matching)
// -----------------------------------------------------------------------------
bool RecoilJets::isTruthPromptIsolatedSignalPhoton(const HepMC::GenEvent* evt,
                                                   const HepMC::GenParticle* pho,
                                                   double& isoEt) const
{
  isoEt = std::numeric_limits<double>::quiet_NaN();
  if (!evt || !pho) return false;

  constexpr double kTruthEtaAbsMax = 0.7;

  // CaloAna-style truth isolation parameters
  constexpr double kIsoConeR  = 0.3;
  constexpr double kMergerDR  = 0.001;
  constexpr double kIsoMaxGeV = 4.0;

  // Final-state photon requirement (your project requirement)
  if (pho->pdg_id() != 22) return false;
  if (pho->status() != 1)  return false;

  const double etaPho = pho->momentum().pseudoRapidity();
  const double phiPho = pho->momentum().phi();
  if (!std::isfinite(etaPho) || !std::isfinite(phiPho)) return false;
  if (std::fabs(etaPho) >= kTruthEtaAbsMax) return false;

  // -------------------------------------------------------------------------
  // 1) Prompt classification: walk back photon-in/photon-out vertices, then topology classify
  // -------------------------------------------------------------------------
  const HepMC::GenVertex* vertex = pho->production_vertex();
  if (!vertex) return false;  // CaloAna returns "can't identify" on null vertex

  // Incoming particles at production vertex
  std::vector<const HepMC::GenParticle*> incoming_particles;
  incoming_particles.reserve(4);
  for (auto inItr = vertex->particles_in_const_begin(); inItr != vertex->particles_in_const_end(); ++inItr)
  {
    if (*inItr) incoming_particles.push_back(*inItr);
  }

  // Walk back while there is exactly one incoming particle and it's a photon (pid==22)
  while (incoming_particles.size() == 1 && incoming_particles[0] && incoming_particles[0]->pdg_id() == 22)
  {
    vertex = incoming_particles[0]->production_vertex();
    if (!vertex) return false;

    incoming_particles.clear();
    for (auto inItr = vertex->particles_in_const_begin(); inItr != vertex->particles_in_const_end(); ++inItr)
    {
      if (*inItr) incoming_particles.push_back(*inItr);
    }
  }

  // Outgoing particles at the (possibly walked-back) vertex
  std::vector<const HepMC::GenParticle*> outgoing_particles;
  outgoing_particles.reserve(4);
  for (auto outItr = vertex->particles_out_const_begin(); outItr != vertex->particles_out_const_end(); ++outItr)
  {
    if (*outItr) outgoing_particles.push_back(*outItr);
  }

  // Make sure there is a photon in outgoing (CaloAna requirement)
  bool has_out_photon = false;
  for (const auto* op : outgoing_particles)
  {
    if (op && op->pdg_id() == 22) { has_out_photon = true; break; }
  }
  if (!has_out_photon) return false;

  // direct photon:1, fragmentation photon:2, decayed photon:3, can't identify:0
  int photonclass = 0;

  // direct photon 2->2: both incoming/outgoing are partons+photon (|pdg|<=22)
  if (incoming_particles.size() == 2 && outgoing_particles.size() == 2)
  {
    const int in0 = incoming_particles[0] ? incoming_particles[0]->pdg_id() : 0;
    const int in1 = incoming_particles[1] ? incoming_particles[1]->pdg_id() : 0;
    const int out0 = outgoing_particles[0] ? outgoing_particles[0]->pdg_id() : 0;
    const int out1 = outgoing_particles[1] ? outgoing_particles[1]->pdg_id() : 0;

    if (std::abs(in0) <= 22 && std::abs(in1) <= 22 &&
        std::abs(out0) <= 22 && std::abs(out1) <= 22)
    {
      photonclass = 1;
    }
  }
  // fragmentation / decay: one incoming
  else if (incoming_particles.size() == 1 && incoming_particles[0])
  {
    const int inpid = incoming_particles[0]->pdg_id();

    // fragmentation photon 1->2: incoming |pid|<=11, outgoing size==2, outgoing contains incoming pid
    if (std::abs(inpid) <= 11 && outgoing_particles.size() == 2)
    {
      bool has_inpid_out = false;
      for (const auto* op : outgoing_particles)
      {
        if (op && op->pdg_id() == inpid) { has_inpid_out = true; break; }
      }
      if (has_inpid_out)
      {
        photonclass = 2;
      }
    }

    // decayed photon: incoming |pid| > 37 (CaloAna sets this after frag check)
    if (std::abs(inpid) > 37)
    {
      photonclass = 3;
    }
  }

  // Accept only direct OR fragmentation (prompt definition used by CaloAna for signal)
  if (!(photonclass == 1 || photonclass == 2)) return false;

  // -------------------------------------------------------------------------
  // 2) Truth isolation: IDENTICAL to CaloAna truth iso method
  //    isoEt = sum_{ΔR<R} Et(p)  -  sum_{ΔR<merger} Et(p)
  //    where merger=0.001 excludes the photon (and any ultra-close merged stuff)
  // -------------------------------------------------------------------------
  auto et_of = [](const HepMC::FourVector& p) -> double
  {
    const double px = p.px();
    const double py = p.py();
    const double pz = p.pz();
    const double e  = p.e();

    const double pt = std::hypot(px, py);
    const double pmag = std::hypot(pt, pz);
    if (!std::isfinite(e) || !std::isfinite(pt) || !std::isfinite(pmag) || pmag <= 0.0) return 0.0;

    // TLorentzVector::Et equivalent: Et = E * pt / |p|
    const double et = e * pt / pmag;
    return (std::isfinite(et) ? et : 0.0);
  };

  double isoSumEt    = 0.0;
  double clusterEt   = 0.0;

  for (auto it = evt->particles_begin(); it != evt->particles_end(); ++it)
  {
    const HepMC::GenParticle* q = *it;
    if (!q) continue;

    // CaloAna uses a "primary_particles" list; operationally these are stable primaries.
    // In HepMC we approximate that with status==1 final-state.
    if (q->status() != 1) continue;

    const double qeta = q->momentum().pseudoRapidity();
    const double qphi = q->momentum().phi();
    if (!std::isfinite(qeta) || !std::isfinite(qphi)) continue;

    const double dphi = TVector2::Phi_mpi_pi(phiPho - qphi);
    const double deta = etaPho - qeta;
    const double dr   = std::sqrt(deta*deta + dphi*dphi);

    const double qEt = et_of(q->momentum());
    if (qEt <= 0.0) continue;

    if (dr < kIsoConeR) isoSumEt += qEt;
    if (dr < kMergerDR) clusterEt += qEt;
  }

  isoEt = isoSumEt - clusterEt;
  if (!std::isfinite(isoEt)) return false;

  return (isoEt < kIsoMaxGeV);
}

// -----------------------------------------------------------------------------
// Unified truth→reco photon matching using CaloRawClusterEval (SIM only)
// -----------------------------------------------------------------------------
bool RecoilJets::findRecoPhotonMatchedToTruthSignal(const HepMC::GenEvent* evt,
                                                    const HepMC::GenParticle* truthPho,
                                                    CaloRawClusterEval& clustereval,
                                                    const RawCluster*& recoPho,
                                                    double& recoPt,
                                                    double& recoEta,
                                                    double& recoPhi,
                                                    double& drBest,
                                                    float& eContribBest) const
{
  recoPho       = nullptr;
  recoPt        = std::numeric_limits<double>::quiet_NaN();
  recoEta       = std::numeric_limits<double>::quiet_NaN();
  recoPhi       = std::numeric_limits<double>::quiet_NaN();
  drBest        = std::numeric_limits<double>::infinity();
  eContribBest  = std::numeric_limits<float>::lowest();

  if (!evt || !truthPho) return false;
  if (!m_photons) return false;

  // Reco-side requirements (your spec)
  constexpr double kRecoEtaAbsMax = 0.7;
  constexpr double kRecoPtMinGeV  = 5.0;
  constexpr double kMatchDR       = 0.05;

  // Truth kinematics
  const double tEta = truthPho->momentum().pseudoRapidity();
  const double tPhi = TVector2::Phi_mpi_pi(truthPho->momentum().phi());
  if (!std::isfinite(tEta) || !std::isfinite(tPhi)) return false;

  const int truthBarcode = truthPho->barcode();

  const auto prange = m_photons->getClusters();
  int iPho = 0;
  for (auto pit = prange.first; pit != prange.second; ++pit, ++iPho)
  {
    const auto* pho = dynamic_cast<const PhotonClusterv1*>(pit->second);
    if (!pho) continue;

    const RawCluster* rc = pho;

    double eta = pho->get_shower_shape_parameter("cluster_eta");
    double phi = TVector2::Phi_mpi_pi(pho->get_shower_shape_parameter("cluster_phi"));
    double pt  = pho->get_shower_shape_parameter("cluster_pt");

    if (!std::isfinite(eta) || !std::isfinite(phi) || !std::isfinite(pt)) continue;
    if (std::fabs(eta) >= kRecoEtaAbsMax) continue;
    if (pt <= kRecoPtMinGeV) continue;

    const double dphi = TVector2::Phi_mpi_pi(tPhi - phi);
    const double deta = tEta - eta;
    const double dr   = std::sqrt(deta*deta + dphi*dphi);
    if (dr >= kMatchDR) continue;

    // BEST-MATCHED truth primary particle (by deposited energy) for this reco cluster
    RawCluster* rc_nc = const_cast<RawCluster*>(rc);
    PHG4Particle* primary = clustereval.max_truth_primary_particle_by_energy(rc_nc);
    if (!primary) continue;

    // Must match the truth photon barcode (truth↔reco association)
    if (primary->get_barcode() != truthBarcode) continue;

    // Defensive: require it's actually a photon primary
    if (primary->get_pid() != 22) continue;

    const float econtrib = clustereval.get_energy_contribution(rc_nc, primary);

    // Choose by largest energy contribution, tie-break by smallest ΔR
    bool take = false;
    if (!recoPho) take = true;
    else if (econtrib > eContribBest) take = true;
    else if (econtrib == eContribBest && dr < drBest) take = true;

    if (take)
    {
      recoPho       = rc;
      recoPt        = pt;
      recoEta       = eta;
      recoPhi       = phi;
      drBest        = dr;
      eContribBest  = econtrib;
    }
  }

  return (recoPho != nullptr);
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

void RecoilJets::initEventDisplayDiagnosticsTree()
{
  if (!m_evtDiagEnabled)
  {
    if (Verbosity() >= 4)
    {
      LOG(4, CLR_CYAN, "[EventDisplayTree][init] skipped: m_evtDiagEnabled=false");
    }
    return;
  }

  if (m_evtDiagTree)
  {
    if (Verbosity() >= 5)
    {
      LOG(5, CLR_CYAN, "[EventDisplayTree][init] already exists: entries=" << m_evtDiagTree->GetEntries());
    }
    return;
  }

  if (!out)
  {
    LOG(1, CLR_YELLOW, "[EventDisplayTree][init] cannot create tree: output TFile pointer is null");
    m_evtDiagEnabled = false;
    return;
  }

  if (!out->IsOpen())
  {
    LOG(1, CLR_YELLOW, "[EventDisplayTree][init] cannot create tree: output file is not open (" << out->GetName() << ")");
    m_evtDiagEnabled = false;
    return;
  }

  if (!out->cd())
  {
    LOG(1, CLR_YELLOW, "[EventDisplayTree][init] cannot create tree: out->cd() failed (" << out->GetName() << ")");
    m_evtDiagEnabled = false;
    return;
  }

  if (Verbosity() >= 2)
  {
    LOG(2, CLR_CYAN, "[EventDisplayTree][init] creating EventDisplayTree in file: " << out->GetName());
  }

  m_evtDiagTree = new TTree("EventDisplayTree",
                            "EventDisplay diagnostics per event per jet radius (offline rendering)");

  if (!m_evtDiagTree)
  {
    LOG(1, CLR_YELLOW, "[EventDisplayTree][init] failed to allocate TTree (nullptr)");
    m_evtDiagEnabled = false;
    return;
  }

  m_evtDiagTree->SetDirectory(out);

  bool ok = true;

  auto AddLeaf = [&](const char* name, void* addr, const char* leaflist)
  {
    TBranch* br = m_evtDiagTree->Branch(name, addr, leaflist);
    if (!br)
    {
      ok = false;
      LOG(1, CLR_YELLOW, "[EventDisplayTree][init] Branch FAILED: \"" << name << "\" leaflist=\"" << leaflist << "\"");
    }
    else if (Verbosity() >= 6)
    {
      LOG(6, CLR_GREEN, "[EventDisplayTree][init] Branch OK: \"" << name << "\" title=\"" << br->GetTitle() << "\"");
    }
  };

  auto AddObj = [&](const char* name, auto* addr)
  {
    TBranch* br = m_evtDiagTree->Branch(name, addr);
    if (!br)
    {
      ok = false;
      LOG(1, CLR_YELLOW, "[EventDisplayTree][init] Branch FAILED: \"" << name << "\" (object)");
    }
    else if (Verbosity() >= 6)
    {
      LOG(6, CLR_GREEN, "[EventDisplayTree][init] Branch OK: \"" << name << "\" class=\"" << br->GetClassName() << "\"");
    }
  };

  AddLeaf("run", &m_evtDiag_run, "run/I");
  AddLeaf("evt", &m_evtDiag_evt, "evt/I");
  AddLeaf("event_count", &m_evtDiag_eventCount, "event_count/L");
  AddLeaf("vz", &m_evtDiag_vz, "vz/F");

  AddObj("rKey", &m_evtDiag_rKey);
  AddLeaf("ptBin", &m_evtDiag_ptBin, "ptBin/I");
  AddLeaf("cat", &m_evtDiag_cat, "cat/I");

  AddLeaf("ptGammaTruth", &m_evtDiag_ptGammaTruth, "ptGammaTruth/F");
  AddLeaf("phiGammaTruth", &m_evtDiag_phiGammaTruth, "phiGammaTruth/F");
  AddLeaf("ptGammaReco", &m_evtDiag_ptGammaReco, "ptGammaReco/F");

  AddLeaf("sel_pt", &m_evtDiag_sel_pt, "sel_pt/F");
  AddLeaf("sel_eta", &m_evtDiag_sel_eta, "sel_eta/F");
  AddLeaf("sel_phi", &m_evtDiag_sel_phi, "sel_phi/F");

  AddLeaf("best_pt", &m_evtDiag_best_pt, "best_pt/F");
  AddLeaf("best_eta", &m_evtDiag_best_eta, "best_eta/F");
  AddLeaf("best_phi", &m_evtDiag_best_phi, "best_phi/F");

  AddLeaf("truthLead_pt", &m_evtDiag_truthLead_pt, "truthLead_pt/F");
  AddLeaf("truthLead_eta", &m_evtDiag_truthLead_eta, "truthLead_eta/F");
  AddLeaf("truthLead_phi", &m_evtDiag_truthLead_phi, "truthLead_phi/F");

  AddLeaf("drSelToTruthLead", &m_evtDiag_drSelToTruthLead, "drSelToTruthLead/F");
  AddLeaf("drBestToTruthLead", &m_evtDiag_drBestToTruthLead, "drBestToTruthLead/F");

  AddObj("sel_calo", &m_evtDiag_sel_calo);
  AddObj("sel_ieta", &m_evtDiag_sel_ieta);
  AddObj("sel_iphi", &m_evtDiag_sel_iphi);
  AddObj("sel_etaTower", &m_evtDiag_sel_etaTower);
  AddObj("sel_phiTower", &m_evtDiag_sel_phiTower);
  AddObj("sel_etTower", &m_evtDiag_sel_etTower);
  AddObj("sel_eTower", &m_evtDiag_sel_eTower);

  AddObj("best_calo", &m_evtDiag_best_calo);
  AddObj("best_ieta", &m_evtDiag_best_ieta);
  AddObj("best_iphi", &m_evtDiag_best_iphi);
  AddObj("best_etaTower", &m_evtDiag_best_etaTower);
  AddObj("best_phiTower", &m_evtDiag_best_phiTower);
  AddObj("best_etTower", &m_evtDiag_best_etTower);
  AddObj("best_eTower", &m_evtDiag_best_eTower);

  if (!ok)
  {
    LOG(1, CLR_YELLOW, "[EventDisplayTree][init] one or more Branch() calls failed; disabling EventDisplayTree (physics unaffected)");
    delete m_evtDiagTree;
    m_evtDiagTree = nullptr;
    m_evtDiagEnabled = false;
    return;
  }

  if (Verbosity() >= 2)
  {
    const int nBr = (m_evtDiagTree->GetListOfBranches() ? m_evtDiagTree->GetListOfBranches()->GetEntries() : -1);
    LOG(2, CLR_GREEN, "[EventDisplayTree][init] OK: booked EventDisplayTree with nBranches=" << nBr);
  }
}




void RecoilJets::resetEventDisplayDiagnosticsBuffers()
{
  m_evtDiag_run        = 0;
  m_evtDiag_evt        = 0;
  m_evtDiag_eventCount = 0;
  m_evtDiag_vz         = 0.0f;

  m_evtDiag_rKey.clear();
  m_evtDiag_ptBin = -1;
  m_evtDiag_cat   = -1;

  m_evtDiag_ptGammaTruth  = 0.0f;
  m_evtDiag_phiGammaTruth = 0.0f;
  m_evtDiag_ptGammaReco   = 0.0f;

  m_evtDiag_sel_pt  = 0.0f;
  m_evtDiag_sel_eta = 0.0f;
  m_evtDiag_sel_phi = 0.0f;

  m_evtDiag_best_pt  = 0.0f;
  m_evtDiag_best_eta = 0.0f;
  m_evtDiag_best_phi = 0.0f;

  m_evtDiag_truthLead_pt  = 0.0f;
  m_evtDiag_truthLead_eta = 0.0f;
  m_evtDiag_truthLead_phi = 0.0f;

  m_evtDiag_drSelToTruthLead  = -1.0f;
  m_evtDiag_drBestToTruthLead = -1.0f;

  m_evtDiag_sel_calo.clear();
  m_evtDiag_sel_ieta.clear();
  m_evtDiag_sel_iphi.clear();
  m_evtDiag_sel_etaTower.clear();
  m_evtDiag_sel_phiTower.clear();
  m_evtDiag_sel_etTower.clear();
  m_evtDiag_sel_eTower.clear();

  m_evtDiag_best_calo.clear();
  m_evtDiag_best_ieta.clear();
  m_evtDiag_best_iphi.clear();
  m_evtDiag_best_etaTower.clear();
  m_evtDiag_best_phiTower.clear();
  m_evtDiag_best_etTower.clear();
  m_evtDiag_best_eTower.clear();
}

bool RecoilJets::eventDisplayDiagnosticsNeed(const std::string& rKey, int ptBin, EventDisplayCat cat)
{
  if (m_evtDiagMaxPerBin <= 0)
  {
    return true;
  }

  std::ostringstream oss;
  oss << rKey << "_ptBin" << ptBin << "_cat" << static_cast<int>(cat);
  const std::string key = oss.str();

  int& n = m_evtDiagSavedPerBin[key];
  if (n >= m_evtDiagMaxPerBin)
  {
    return false;
  }
  ++n;
  return true;
}

void RecoilJets::appendEventDisplayDiagnosticsFromJet(const Jet* jet,
                                                      std::vector<int>& calo,
                                                      std::vector<int>& ieta,
                                                      std::vector<int>& iphi,
                                                      std::vector<float>& eta,
                                                      std::vector<float>& phi,
                                                      std::vector<float>& et,
                                                      std::vector<float>& e) const
{
  if (!jet)
  {
    return;
  }

    for (auto it = jet->begin_comp(); it != jet->end_comp(); ++it)
    {
      const Jet::SRC src = it->first;
      const unsigned int idx = it->second;

      TowerInfoContainer* towers = nullptr;
      RawTowerGeomContainer* geom = nullptr;
      RawTowerDefs::CalorimeterId caloId = RawTowerDefs::CalorimeterId::CEMC;
      int caloCode = -1;

      if (src == Jet::CEMC_TOWERINFO || src == Jet::CEMC_TOWERINFO_RETOWER)
      {
        towers = m_evtDispTowersCEMC;
        geom   = m_evtDispGeomCEMC;
        caloId = RawTowerDefs::CalorimeterId::CEMC;
        caloCode = 0;
      }
      else if (src == Jet::HCALIN_TOWERINFO)
      {
        towers = m_evtDispTowersIHCal;
        geom   = m_evtDispGeomIHCal;
        caloId = RawTowerDefs::CalorimeterId::HCALIN;
        caloCode = 1;
      }
      else if (src == Jet::HCALOUT_TOWERINFO)
      {
        towers = m_evtDispTowersOHCal;
        geom   = m_evtDispGeomOHCal;
        caloId = RawTowerDefs::CalorimeterId::HCALOUT;
        caloCode = 2;
      }
      else
      {
        continue;
      }

      if (!towers || !geom)
      {
        continue;
      }

      const auto nTowers = towers->size();
      if (nTowers == 0)
      {
        continue;
      }

      // TowerJetInput stores TowerInfo constituents as CHANNEL indices.
      // Treat compid as channel only; never decode_key here (can segfault if misused).
      if (idx >= nTowers)
      {
        continue;
      }

      TowerInfo* tower = towers->get_tower_at_channel(idx);
      if (!tower)
      {
        continue;
      }

      const unsigned int calokey = towers->encode_key(idx);
      const int twr_ieta = towers->getTowerEtaBin(calokey);
      const int twr_iphi = towers->getTowerPhiBin(calokey);

    RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(caloId, twr_ieta, twr_iphi);
    RawTowerGeom* tower_geom = geom->get_tower_geometry(key);
    if (!tower_geom)
    {
      continue;
    }

    const double x = tower_geom->get_center_x();
    const double y = tower_geom->get_center_y();
    const double r = std::sqrt(x * x + y * y);
    if (!(r > 0))
    {
      continue;
    }

    const double tower_phi = std::atan2(y, x);

    const double eta0 = tower_geom->get_eta();
    const double z0 = std::sinh(eta0) * r;
    const double z = z0 - m_vz;
    const double tower_eta = std::asinh(z / r);

    const double E = tower->get_energy();
    const double Et = E / std::cosh(tower_eta);

    if (!std::isfinite(Et))
    {
      continue;
    }

    calo.push_back(caloCode);
    ieta.push_back(twr_ieta);
    iphi.push_back(twr_iphi);
    eta.push_back(static_cast<float>(tower_eta));
    phi.push_back(static_cast<float>(tower_phi));
    et.push_back(static_cast<float>(Et));
    e.push_back(static_cast<float>(E));
  }
}

void RecoilJets::fillEventDisplayDiagnostics(const std::string& rKey,
                                             int ptBin,
                                             EventDisplayCat cat,
                                             double truthGammaPt,
                                             double truthGammaPhi,
                                             double recoGammaPt,
                                             const Jet* selectedRecoilJet,
                                             const Jet* recoTruthBest,
                                             const Jet* truthLeadRecoilJet)
{
  if (!m_evtDiagEnabled || !m_evtDiagNodesReady)
  {
    return;
  }

  initEventDisplayDiagnosticsTree();
  if (!m_evtDiagTree)
  {
    return;
  }

  resetEventDisplayDiagnosticsBuffers();

  if (m_evtHeader)
  {
    m_evtDiag_run = m_evtHeader->get_RunNumber();
    m_evtDiag_evt = m_evtHeader->get_EvtSequence();
  }
  m_evtDiag_eventCount = event_count;
  m_evtDiag_vz = static_cast<float>(m_vz);

  m_evtDiag_rKey   = rKey;
  m_evtDiag_ptBin  = ptBin;
  m_evtDiag_cat    = static_cast<int>(cat);

  m_evtDiag_ptGammaTruth  = static_cast<float>(truthGammaPt);
  m_evtDiag_phiGammaTruth = static_cast<float>(truthGammaPhi);
  m_evtDiag_ptGammaReco   = static_cast<float>(recoGammaPt);

  if (selectedRecoilJet)
  {
    m_evtDiag_sel_pt  = static_cast<float>(selectedRecoilJet->get_pt());
    m_evtDiag_sel_eta = static_cast<float>(selectedRecoilJet->get_eta());
    m_evtDiag_sel_phi = static_cast<float>(selectedRecoilJet->get_phi());
  }

  if (recoTruthBest)
  {
    m_evtDiag_best_pt  = static_cast<float>(recoTruthBest->get_pt());
    m_evtDiag_best_eta = static_cast<float>(recoTruthBest->get_eta());
    m_evtDiag_best_phi = static_cast<float>(recoTruthBest->get_phi());
  }

  if (truthLeadRecoilJet)
  {
    m_evtDiag_truthLead_pt  = static_cast<float>(truthLeadRecoilJet->get_pt());
    m_evtDiag_truthLead_eta = static_cast<float>(truthLeadRecoilJet->get_eta());
    m_evtDiag_truthLead_phi = static_cast<float>(truthLeadRecoilJet->get_phi());
  }

  auto dR = [](double eta1, double phi1, double eta2, double phi2)
  {
    const double dphi = TVector2::Phi_mpi_pi(phi1 - phi2);
    const double deta = eta1 - eta2;
    return std::sqrt(deta * deta + dphi * dphi);
  };

  if (selectedRecoilJet && truthLeadRecoilJet)
  {
    m_evtDiag_drSelToTruthLead = static_cast<float>(dR(selectedRecoilJet->get_eta(),
                                                       selectedRecoilJet->get_phi(),
                                                       truthLeadRecoilJet->get_eta(),
                                                       truthLeadRecoilJet->get_phi()));
  }

  if (recoTruthBest && truthLeadRecoilJet)
  {
    m_evtDiag_drBestToTruthLead = static_cast<float>(dR(recoTruthBest->get_eta(),
                                                        recoTruthBest->get_phi(),
                                                        truthLeadRecoilJet->get_eta(),
                                                        truthLeadRecoilJet->get_phi()));
  }

    appendEventDisplayDiagnosticsFromJet(selectedRecoilJet,
                                         m_evtDiag_sel_calo,
                                         m_evtDiag_sel_ieta,
                                         m_evtDiag_sel_iphi,
                                         m_evtDiag_sel_etaTower,
                                         m_evtDiag_sel_phiTower,
                                         m_evtDiag_sel_etTower,
                                         m_evtDiag_sel_eTower);

    appendEventDisplayDiagnosticsFromJet(recoTruthBest,
                                         m_evtDiag_best_calo,
                                         m_evtDiag_best_ieta,
                                         m_evtDiag_best_iphi,
                                         m_evtDiag_best_etaTower,
                                         m_evtDiag_best_phiTower,
                                         m_evtDiag_best_etTower,
                                         m_evtDiag_best_eTower);

    const bool hasSelTowers  = !m_evtDiag_sel_etTower.empty();
    const bool hasBestTowers = !m_evtDiag_best_etTower.empty();
    const bool hasAnyTowers  = (hasSelTowers || hasBestTowers);

    ++m_evtDiagNFill;

    const int icat = static_cast<int>(cat);
    if (icat >= 0 && icat < 3)
    {
      ++m_evtDiagNFillByCat[icat];
      if (hasAnyTowers) ++m_evtDiagNFillWithAnyTowersByCat[icat];
    }

    if (hasSelTowers)  ++m_evtDiagNFillWithSelTowers;
    if (hasBestTowers) ++m_evtDiagNFillWithBestTowers;
    if (hasAnyTowers)  ++m_evtDiagNFillWithAnyTowers;

    m_evtDiagTree->Fill();
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


TH1F* RecoilJets::getOrBookXJHist(const std::string& trig,
                                 const std::string& rKey,
                                 int etIdx, int centIdx)
{
  const std::string base   = "h_xJ";
  const std::string suffix = suffixForBins(etIdx, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (Verbosity() >= 5)
    LOG(5, CLR_BLUE, "  [getOrBookXJHist] trig=\"" << trig
           << "\" rKey=\"" << rKey
           << "\" etIdx=" << etIdx << " centIdx=" << centIdx
           << " → name=\"" << name << "\"");

  if (trig.empty() || rKey.empty())
  {
    LOG(2, CLR_YELLOW, "  [getOrBookXJHist] empty trig/rKey – returning nullptr");
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

  const int    nbins = 60;
  const double xmin  = 0.0;
  const double xmax  = 3.0;

  const std::string title = name + ";x_{J}=p_{T}^{jet1}/p_{T}^{#gamma};Entries";

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


TH1F* RecoilJets::getOrBookJet1PtHist(const std::string& trig,
                                     const std::string& rKey,
                                     int etIdx, int centIdx)
{
  const std::string base   = "h_jet1Pt";
  const std::string suffix = suffixForBins(etIdx, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

TH1F* RecoilJets::getOrBookJet2PtHist(const std::string& trig,
                                     const std::string& rKey,
                                     int etIdx, int centIdx)
{
  const std::string base   = "h_jet2Pt";
  const std::string suffix = suffixForBins(etIdx, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

TH1F* RecoilJets::getOrBookAlphaHist(const std::string& trig,
                                    const std::string& rKey,
                                    int etIdx, int centIdx)
{
  const std::string base   = "h_alpha";
  const std::string suffix = suffixForBins(etIdx, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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


// ============================================================================
//  Unfolding bookers for inclusive γ–jet pairing (ATLAS-style 2D unfolding)
//   - Uses PPG12 photon bin edges (reco + truth)
//   - Uses log-ish xJ edges with explicit under/overflow bins
//   - All are TH2F, so your End() writer will write them automatically.
// ============================================================================

TH2F* RecoilJets::getOrBookUnfoldRecoPtXJIncl(const std::string& trig,
                                              const std::string& rKey,
                                              int centIdx)
{
  const std::string base   = "h2_unfoldReco_pTgamma_xJ_incl";
  const std::string suffix = suffixForBins(-1, centIdx);   // centrality-only (Au+Au); empty in pp
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

    // Reco photon edges for unfolding (includes [10,15] underflow and [35,40] overflow)
    const std::vector<double>& kPtReco = m_unfoldRecoPhotonPtBins;

    // ATLAS-like log-ish xJ edges, with explicit under/overflow bins for unfolding stability
    const std::vector<double>& kXJ = m_unfoldXJBins;

    const int nx = static_cast<int>(kPtReco.size()) - 1;
    const int ny = static_cast<int>(kXJ.size()) - 1;

  const std::string title =
    name + ";p_{T}^{#gamma,reco} [GeV];x_{J#gamma}^{reco}=p_{T}^{jet,reco}/p_{T}^{#gamma,reco}";

  auto* h = new TH2F(name.c_str(), title.c_str(),
                     nx, kPtReco.data(),
                     ny, kXJ.data());
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookUnfoldRecoPtDphiIncl(const std::string& trig,
                                              const std::string& rKey,
                                              int centIdx)
{
    const std::string base   = "h2_unfoldReco_pTgamma_dphi_incl";
    const std::string suffix = suffixForBins(-1, centIdx);   // centrality-only (Au+Au); empty in pp
    const std::string name   = base + "_" + rKey + suffix;

    if (trig.empty() || rKey.empty()) return nullptr;

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

    // Must match reco pTγ unfolding binning (includes [10,15] underflow and [35,40] overflow)
    const std::vector<double>& kPtReco = m_unfoldRecoPhotonPtBins;

    const int nx = static_cast<int>(kPtReco.size()) - 1;
    const int ny = 64;

    const std::string title =
      name + ";p_{T}^{#gamma,reco} [GeV];|#Delta#phi(#gamma,jet)| [rad]";

    auto* h = new TH2F(name.c_str(), title.c_str(),
                       nx, kPtReco.data(),
                       ny, 0.0, M_PI);
    h->Sumw2();

    H[name] = h;
    if (prevDir) prevDir->cd();
    return h;
}

TH2F* RecoilJets::getOrBookUnfoldTruthPtDphiIncl(const std::string& trig,
                                               const std::string& rKey,
                                               int centIdx)
{
    const std::string base   = "h2_unfoldTruth_pTgamma_dphi_incl";
    const std::string suffix = suffixForBins(-1, centIdx);
    const std::string name   = base + "_" + rKey + suffix;

    if (trig.empty() || rKey.empty()) return nullptr;

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

    // Must match truth pTγ unfolding binning you already use
      const std::vector<double>& kPtTruth = m_unfoldTruthPhotonPtBins;

    const int nx = static_cast<int>(kPtTruth.size()) - 1;
    const int ny = 64;

    const std::string title =
      name + ";p_{T}^{#gamma,truth} [GeV];|#Delta#phi(#gamma^{truth},jet^{truth})| [rad]";

    auto* h = new TH2F(name.c_str(), title.c_str(),
                       nx, kPtTruth.data(),
                       ny, 0.0, M_PI);
    h->Sumw2();

    H[name] = h;
    if (prevDir) prevDir->cd();
    return h;
}

TH2F* RecoilJets::getOrBookUnfoldTruthPtXJIncl(const std::string& trig,
                                               const std::string& rKey,
                                               int centIdx)
{
  const std::string base   = "h2_unfoldTruth_pTgamma_xJ_incl";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

    // Truth photon edges for unfolding (includes [5,10] and [10,15] underflow, plus [35,40] overflow)
    const std::vector<double>& kPtTruth = m_unfoldTruthPhotonPtBins;

    // Must match reco xJ edges for a clean response definition
    const std::vector<double>& kXJ = m_unfoldXJBins;

    const int nx = static_cast<int>(kPtTruth.size()) - 1;
    const int ny = static_cast<int>(kXJ.size()) - 1;

  const std::string title =
    name + ";p_{T}^{#gamma,truth} [GeV];x_{J#gamma}^{truth}=p_{T}^{jet,truth}/p_{T}^{#gamma,truth}";

  auto* h = new TH2F(name.c_str(), title.c_str(),
                     nx, kPtTruth.data(),
                     ny, kXJ.data());
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookUnfoldResponsePtXJIncl(const std::string& trig,
                                                  const std::string& rKey,
                                                  int centIdx)
{
  const std::string base   = "h2_unfoldResponse_pTgamma_xJ_incl";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

    // These define the TH2 template bin counts (FindBin() global bin indexing includes +2 each axis)
    const std::vector<double>& kPtReco  = m_unfoldRecoPhotonPtBins;   // nbinsX(reco)=8
    const std::vector<double>& kPtTruth = m_unfoldTruthPhotonPtBins;  // nbinsX(truth)=9
    const std::vector<double>& kXJ      = m_unfoldXJBins;             // nbinsY=15

    const int nPtReco  = static_cast<int>(kPtReco.size())  - 1;
    const int nPtTruth = static_cast<int>(kPtTruth.size()) - 1;
    const int nXJ      = static_cast<int>(kXJ.size())      - 1;

  // Global bin count used by TH2::FindBin for a TH2 with (nPt, nXJ):
  //   globalBin in [0, (nPt+2)*(nXJ+2)-1]
  const int nGlobTruth = (nPtTruth + 2) * (nXJ + 2);
  const int nGlobReco  = (nPtReco  + 2) * (nXJ + 2);

  const std::string title =
    name + ";global bin (truth: p_{T}^{#gamma}, x_{J#gamma});global bin (reco: p_{T}^{#gamma}, x_{J#gamma})";

  auto* h = new TH2F(name.c_str(), title.c_str(),
                     nGlobTruth, -0.5, static_cast<double>(nGlobTruth) - 0.5,
                     nGlobReco,  -0.5, static_cast<double>(nGlobReco)  - 0.5);
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookUnfoldRecoFakesPtXJIncl(const std::string& trig,
                                                   const std::string& rKey,
                                                   int centIdx)
{
  const std::string base   = "h2_unfoldRecoFakes_pTgamma_xJ_incl";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

    const std::vector<double>& kPtReco = m_unfoldRecoPhotonPtBins;
    const std::vector<double>& kXJ = m_unfoldXJBins;

    const int nx = static_cast<int>(kPtReco.size()) - 1;
    const int ny = static_cast<int>(kXJ.size()) - 1;

  const std::string title =
    name + ";p_{T}^{#gamma,reco} [GeV];x_{J#gamma}^{reco} (FAKES)";

  auto* h = new TH2F(name.c_str(), title.c_str(),
                     nx, kPtReco.data(),
                     ny, kXJ.data());
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookUnfoldTruthMissesPtXJIncl(const std::string& trig,
                                                     const std::string& rKey,
                                                     int centIdx)
{
  const std::string base   = "h2_unfoldTruthMisses_pTgamma_xJ_incl";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

    const std::vector<double>& kPtTruth = m_unfoldTruthPhotonPtBins;
  const std::vector<double>& kXJ = m_unfoldXJBins;

  const int nx = static_cast<int>(kPtTruth.size()) - 1;
  const int ny = static_cast<int>(kXJ.size()) - 1;

  const std::string title =
    name + ";p_{T}^{#gamma,truth} [GeV];x_{J#gamma}^{truth} (MISSES)";

  auto* h = new TH2F(name.c_str(), title.c_str(),
                     nx, kPtTruth.data(),
                     ny, kXJ.data());
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

// -------------------------------------------------------------------------
// Photon-only unfolding (for N_gamma normalization)
// -------------------------------------------------------------------------
TH1F* RecoilJets::getOrBookUnfoldRecoPhoPtGamma(const std::string& trig,
                                                int centIdx)
{
  const std::string base   = "h_unfoldRecoPho_pTgamma";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + suffix;

  if (trig.empty()) return nullptr;

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

    const std::vector<double>& kPtReco = m_unfoldRecoPhotonPtBins;
    const int nb = static_cast<int>(kPtReco.size()) - 1;

  const std::string title =
    name + ";p_{T}^{#gamma,reco} [GeV];N_{#gamma}^{reco}";

  auto* h = new TH1F(name.c_str(), title.c_str(), nb, kPtReco.data());
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH1F* RecoilJets::getOrBookUnfoldTruthPhoPtGamma(const std::string& trig,
                                                 int centIdx)
{
  const std::string base   = "h_unfoldTruthPho_pTgamma";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + suffix;

  if (trig.empty()) return nullptr;

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

    const std::vector<double>& kPtTruth = m_unfoldTruthPhotonPtBins;
  const int nb = static_cast<int>(kPtTruth.size()) - 1;

  const std::string title =
    name + ";p_{T}^{#gamma,truth} [GeV];N_{#gamma}^{truth}";

  auto* h = new TH1F(name.c_str(), title.c_str(), nb, kPtTruth.data());
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookUnfoldResponsePhoPtGamma(const std::string& trig,
                                                    int centIdx)
{
  const std::string base   = "h2_unfoldResponsePho_pTgamma";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + suffix;

  if (trig.empty()) return nullptr;

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

    const std::vector<double>& kPtTruth = m_unfoldTruthPhotonPtBins;
    const std::vector<double>& kPtReco  = m_unfoldRecoPhotonPtBins;

    const int nx = static_cast<int>(kPtTruth.size()) - 1;
    const int ny = static_cast<int>(kPtReco.size())  - 1;

  const std::string title =
    name + ";p_{T}^{#gamma,truth} [GeV];p_{T}^{#gamma,reco} [GeV]";

  auto* h = new TH2F(name.c_str(), title.c_str(),
                     nx, kPtTruth.data(),
                     ny, kPtReco.data());
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH1F* RecoilJets::getOrBookUnfoldRecoPhoFakesPtGamma(const std::string& trig,
                                                     int centIdx)
{
  const std::string base   = "h_unfoldRecoPhoFakes_pTgamma";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + suffix;

  if (trig.empty()) return nullptr;

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

    const std::vector<double>& kPtReco = m_unfoldRecoPhotonPtBins;
    const int nb = static_cast<int>(kPtReco.size()) - 1;

  const std::string title =
    name + ";p_{T}^{#gamma,reco} [GeV];Reco fakes";

  auto* h = new TH1F(name.c_str(), title.c_str(), nb, kPtReco.data());
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH1F* RecoilJets::getOrBookUnfoldTruthPhoMissesPtGamma(const std::string& trig,
                                                       int centIdx)
{
  const std::string base   = "h_unfoldTruthPhoMisses_pTgamma";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + suffix;

  if (trig.empty()) return nullptr;

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

    const std::vector<double>& kPtTruth = m_unfoldTruthPhotonPtBins;
  const int nb = static_cast<int>(kPtTruth.size()) - 1;

  const std::string title =
    name + ";p_{T}^{#gamma,truth} [GeV];Truth misses";

  auto* h = new TH1F(name.c_str(), title.c_str(), nb, kPtTruth.data());
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

// -------------------------------------------------------------------------
// Unfolding QA helpers (SIM ONLY): explicit matched + type-split fakes/misses + jet-match quality
// -------------------------------------------------------------------------
TH2F* RecoilJets::getOrBookUnfoldTruthMatchedPtXJIncl(const std::string& trig,
                                                      const std::string& rKey,
                                                      int centIdx)
{
  const std::string base   = "h2_unfoldTruthMatched_pTgamma_xJ_incl";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

  const std::vector<double>& kPtTruth = m_unfoldTruthPhotonPtBins;
  const std::vector<double>& kXJ = m_unfoldXJBins;

  const int nx = static_cast<int>(kPtTruth.size()) - 1;
  const int ny = static_cast<int>(kXJ.size()) - 1;

  const std::string title =
    name + ";p_{T}^{#gamma,truth} [GeV];x_{J#gamma}^{truth} (MATCHED)";

  auto* h = new TH2F(name.c_str(), title.c_str(),
                     nx, kPtTruth.data(),
                     ny, kXJ.data());
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookUnfoldRecoMatchedPtXJIncl(const std::string& trig,
                                                     const std::string& rKey,
                                                     int centIdx)
{
  const std::string base   = "h2_unfoldRecoMatched_pTgamma_xJ_incl";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

  const std::vector<double>& kPtReco = m_unfoldRecoPhotonPtBins;
  const std::vector<double>& kXJ = m_unfoldXJBins;

  const int nx = static_cast<int>(kPtReco.size()) - 1;
  const int ny = static_cast<int>(kXJ.size()) - 1;

  const std::string title =
    name + ";p_{T}^{#gamma,reco} [GeV];x_{J#gamma}^{reco} (MATCHED)";

  auto* h = new TH2F(name.c_str(), title.c_str(),
                     nx, kPtReco.data(),
                     ny, kXJ.data());
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookUnfoldRecoFakesPtXJIncl_typeA(const std::string& trig,
                                                         const std::string& rKey,
                                                         int centIdx)
{
  const std::string base   = "h2_unfoldRecoFakes_typeA_pTgamma_xJ_incl";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

    const std::vector<double>& kPtReco = m_unfoldRecoPhotonPtBins;
    const std::vector<double>& kXJ = m_unfoldXJBins;

    const int nx = static_cast<int>(kPtReco.size()) - 1;
    const int ny = static_cast<int>(kXJ.size()) - 1;

  const std::string title =
    name + ";p_{T}^{#gamma,reco} [GeV];x_{J#gamma}^{reco} (FAKES typeA)";

  auto* h = new TH2F(name.c_str(), title.c_str(),
                     nx, kPtReco.data(),
                     ny, kXJ.data());
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookUnfoldRecoFakesPtXJIncl_typeB(const std::string& trig,
                                                         const std::string& rKey,
                                                         int centIdx)
{
  const std::string base   = "h2_unfoldRecoFakes_typeB_pTgamma_xJ_incl";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

    const std::vector<double>& kPtReco = m_unfoldRecoPhotonPtBins;
    const std::vector<double>& kXJ = m_unfoldXJBins;

    const int nx = static_cast<int>(kPtReco.size()) - 1;
    const int ny = static_cast<int>(kXJ.size()) - 1;

  const std::string title =
    name + ";p_{T}^{#gamma,reco} [GeV];x_{J#gamma}^{reco} (FAKES typeB)";

  auto* h = new TH2F(name.c_str(), title.c_str(),
                     nx, kPtReco.data(),
                     ny, kXJ.data());
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookUnfoldTruthMissesPtXJIncl_typeA(const std::string& trig,
                                                           const std::string& rKey,
                                                           int centIdx)
{
  const std::string base   = "h2_unfoldTruthMisses_typeA_pTgamma_xJ_incl";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

  const std::vector<double>& kPtTruth = m_unfoldTruthPhotonPtBins;
  const std::vector<double>& kXJ = m_unfoldXJBins;

  const int nx = static_cast<int>(kPtTruth.size()) - 1;
  const int ny = static_cast<int>(kXJ.size()) - 1;

  const std::string title =
    name + ";p_{T}^{#gamma,truth} [GeV];x_{J#gamma}^{truth} (MISSES typeA)";

  auto* h = new TH2F(name.c_str(), title.c_str(),
                     nx, kPtTruth.data(),
                     ny, kXJ.data());
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookUnfoldTruthMissesPtXJIncl_typeB(const std::string& trig,
                                                           const std::string& rKey,
                                                           int centIdx)
{
  const std::string base   = "h2_unfoldTruthMisses_typeB_pTgamma_xJ_incl";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

    const std::vector<double>& kPtTruth = m_unfoldTruthPhotonPtBins;
  const std::vector<double>& kXJ = m_unfoldXJBins;

  const int nx = static_cast<int>(kPtTruth.size()) - 1;
  const int ny = static_cast<int>(kXJ.size()) - 1;

  const std::string title =
    name + ";p_{T}^{#gamma,truth} [GeV];x_{J#gamma}^{truth} (MISSES typeB)";

  auto* h = new TH2F(name.c_str(), title.c_str(),
                     nx, kPtTruth.data(),
                     ny, kXJ.data());
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH1F* RecoilJets::getOrBookUnfoldJetMatchDR(const std::string& trig,
                                            const std::string& rKey,
                                            int centIdx)
{
  const std::string base   = "h_unfoldJetMatch_dR";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

  auto* h = new TH1F(name.c_str(),
                     (name + ";#DeltaR(reco jet, truth jet);Counts").c_str(),
                     60, 0.0, 0.30);
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookUnfoldJetPtResponsePtTruth(const std::string& trig,
                                                      const std::string& rKey,
                                                      int centIdx)
{
  const std::string base   = "h2_unfoldJetPtResponse_pTtruth_ratio";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

  TH2F* h = nullptr;

  std::vector<double> tmpJetPt;
  const std::vector<double>* pJetPt = &m_unfoldJetPtBins;
  if (pJetPt->size() < 2)
  {
      tmpJetPt.reserve(121);
      for (int i = 0; i <= 120; ++i) tmpJetPt.push_back(0.5 * (double)i);
      pJetPt = &tmpJetPt;
  }

  const int nb = static_cast<int>(pJetPt->size()) - 1;
  h = new TH2F(name.c_str(),
                (name + ";p_{T}^{jet,truth} [GeV];p_{T}^{jet,reco}/p_{T}^{jet,truth}").c_str(),
                nb, pJetPt->data(),
                120, 0.0, 2.0);
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookUnfoldJetPtResponseAllPtTruth(const std::string& trig,
                                                         const std::string& rKey,
                                                         int centIdx)
{
  const std::string base   = "h2_unfoldJetPtResponseAll_pTtruth_ratio";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

  TH2F* h = nullptr;

  std::vector<double> tmpJetPt;
  const std::vector<double>* pJetPt = &m_unfoldJetPtBins;
  if (pJetPt->size() < 2)
  {
      tmpJetPt.reserve(121);
      for (int i = 0; i <= 120; ++i) tmpJetPt.push_back(0.5 * (double)i);
      pJetPt = &tmpJetPt;
  }

  const int nb = static_cast<int>(pJetPt->size()) - 1;
  h = new TH2F(name.c_str(),
                (name + ";p_{T}^{jet,truth} [GeV];p_{T}^{jet,reco}/p_{T}^{jet,truth} (all matched fid jets)").c_str(),
                nb, pJetPt->data(),
                120, 0.0, 2.0);
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookLeadRecoilJetPtResponsePtTruth(const std::string& trig,
                                                          const std::string& rKey,
                                                          int centIdx)
{
  const std::string base   = "h2_leadRecoilJetPtResponse_pTtruth_ratio";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

  TH2F* h = nullptr;

  std::vector<double> tmpJetPt;
  tmpJetPt.reserve(121);
  for (int i = 0; i <= 120; ++i) tmpJetPt.push_back(0.5 * (double)i);
  const std::vector<double>* pJetPt = &tmpJetPt;

  const int nb = static_cast<int>(pJetPt->size()) - 1;
  h = new TH2F(name.c_str(),
                (name + ";p_{T}^{jet,truth} [GeV];p_{T}^{jet,reco}/p_{T}^{jet,truth} (lead recoil jet1)").c_str(),
                nb, pJetPt->data(),
                120, 0.0, 2.0);
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookLeadRecoilJetPtTruthPtReco(const std::string& trig,
                                                      const std::string& rKey,
                                                      int centIdx)
{
  const std::string base   = "h2_leadRecoilJet_pTtruth_pTreco";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

  TH2F* h = nullptr;

  std::vector<double> tmpJetPt;
  tmpJetPt.reserve(121);
  for (int i = 0; i <= 120; ++i) tmpJetPt.push_back(0.5 * (double)i);
  const std::vector<double>* pJetPt = &tmpJetPt;

  const int nb = static_cast<int>(pJetPt->size()) - 1;
  h = new TH2F(name.c_str(),
                (name + ";p_{T}^{jet,truth} [GeV];p_{T}^{jet,reco} [GeV] (lead recoil jet1)").c_str(),
                nb, pJetPt->data(),
                nb, pJetPt->data());
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH1F* RecoilJets::getOrBookLeadRecoilJetMatchDR(const std::string& trig,
                                                const std::string& rKey,
                                                int centIdx)
{
  const std::string base   = "h_leadRecoilJetMatch_dR";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

  auto* h = new TH1F(name.c_str(),
                     (name + ";#DeltaR(lead recoil jet1, matched truth jet);Counts").c_str(),
                     60, 0.0, 0.30);
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

// -------------------------------------------------------------------------
//  JES3-style *leading truth recoil jet1* match bookkeeping vs truth pT^gamma
//   names:
//     h_leadTruthRecoilMatch_den_pTgammaTruth_<rKey><centSuffix>
//     h_leadTruthRecoilMatch_num_pTgammaTruth_<rKey><centSuffix>
//     h_leadTruthRecoilMatch_missA_pTgammaTruth_<rKey><centSuffix>
//     h_leadTruthRecoilMatch_missB_pTgammaTruth_<rKey><centSuffix>
// -------------------------------------------------------------------------
TH1F* RecoilJets::getOrBookLeadTruthRecoilMatchDenPtGammaTruth(const std::string& trig,
                                                               const std::string& rKey,
                                                               int centIdx)
{
    const std::string base   = "h_leadTruthRecoilMatch_den_pTgammaTruth";
    const std::string suffix = suffixForBins(-1, centIdx);
    const std::string name   = base + "_" + rKey + suffix;

    if (trig.empty() || rKey.empty()) return nullptr;

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

    const std::vector<double>& kPtTruth = m_gammaPtBins;
    const int nb = static_cast<int>(kPtTruth.size()) - 1;

    const std::string title =
      name + ";p_{T}^{#gamma,truth} [GeV];Den (truth lead recoil jet1 exists)";

    auto* h = new TH1F(name.c_str(), title.c_str(), nb, kPtTruth.data());
    h->Sumw2();

    H[name] = h;
    if (prevDir) prevDir->cd();
    return h;
}

TH1F* RecoilJets::getOrBookLeadTruthRecoilMatchNumPtGammaTruth(const std::string& trig,
                                                               const std::string& rKey,
                                                               int centIdx)
{
    const std::string base   = "h_leadTruthRecoilMatch_num_pTgammaTruth";
    const std::string suffix = suffixForBins(-1, centIdx);
    const std::string name   = base + "_" + rKey + suffix;

    if (trig.empty() || rKey.empty()) return nullptr;

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

    const std::vector<double>& kPtTruth = m_gammaPtBins;
    const int nb = static_cast<int>(kPtTruth.size()) - 1;

    const std::string title =
      name + ";p_{T}^{#gamma,truth} [GeV];Num (reco recoil jet1 matches truth jet1)";

    auto* h = new TH1F(name.c_str(), title.c_str(), nb, kPtTruth.data());
    h->Sumw2();

    H[name] = h;
    if (prevDir) prevDir->cd();
    return h;
}

TH1F* RecoilJets::getOrBookLeadTruthRecoilMatchMissA_PtGammaTruth(const std::string& trig,
                                                                    const std::string& rKey,
                                                                    int centIdx)
  {
    const std::string base   = "h_leadTruthRecoilMatch_missA_pTgammaTruth";
    const std::string suffix = suffixForBins(-1, centIdx);
    const std::string name   = base + "_" + rKey + suffix;

    if (trig.empty() || rKey.empty()) return nullptr;

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

    const std::vector<double>& kPtTruth = m_gammaPtBins;
    const int nb = static_cast<int>(kPtTruth.size()) - 1;

    const std::string title =
      name + ";p_{T}^{#gamma,truth} [GeV];MissA (reco fid jet near truth jet1, but Num failed)";

    auto* h = new TH1F(name.c_str(), title.c_str(), nb, kPtTruth.data());
    h->Sumw2();

    H[name] = h;
    if (prevDir) prevDir->cd();
    return h;
}

TH1F* RecoilJets::getOrBookLeadTruthRecoilMatchMissA1_PtGammaTruth(const std::string& trig,
                                                                   const std::string& rKey,
                                                                   int centIdx)
{
    const std::string base   = "h_leadTruthRecoilMatch_missA1_pTgammaTruth";
    const std::string suffix = suffixForBins(-1, centIdx);
    const std::string name   = base + "_" + rKey + suffix;

    if (trig.empty() || rKey.empty()) return nullptr;

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

    const std::vector<double>& kPtTruth = m_gammaPtBins;
    const int nb = static_cast<int>(kPtTruth.size()) - 1;

    const std::string title =
      name + ";p_{T}^{#gamma,truth} [GeV];MissA1 (truth-match passes recoil cuts; competitor/ordering)";

    auto* h = new TH1F(name.c_str(), title.c_str(), nb, kPtTruth.data());
    h->Sumw2();

    H[name] = h;
    if (prevDir) prevDir->cd();
    return h;
}

TH1F* RecoilJets::getOrBookLeadTruthRecoilMatchMissA2_PtGammaTruth(const std::string& trig,
                                                                   const std::string& rKey,
                                                                   int centIdx)
{
    const std::string base   = "h_leadTruthRecoilMatch_missA2_pTgammaTruth";
    const std::string suffix = suffixForBins(-1, centIdx);
    const std::string name   = base + "_" + rKey + suffix;

    if (trig.empty() || rKey.empty()) return nullptr;

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

    const std::vector<double>& kPtTruth = m_gammaPtBins;
    const int nb = static_cast<int>(kPtTruth.size()) - 1;

    const std::string title =
      name + ";p_{T}^{#gamma,truth} [GeV];MissA2 (truth-match fails recoil cuts; gate-exclusion)";

    auto* h = new TH1F(name.c_str(), title.c_str(), nb, kPtTruth.data());
    h->Sumw2();

    H[name] = h;
    if (prevDir) prevDir->cd();
    return h;
}

TH1I* RecoilJets::getOrBookLeadTruthRecoilMatchMissA2_Cutflow(const std::string& trig,
                                                              const std::string& rKey,
                                                              int centIdx)
{
    const std::string base   = "h_leadTruthRecoilMatch_missA2_cutflow";
    const std::string suffix = suffixForBins(-1, centIdx);
    const std::string name   = base + "_" + rKey + suffix;

    if (trig.empty() || rKey.empty()) return nullptr;

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

    auto* h = new TH1I(name.c_str(),
                       (name + ";cut stage;Counts (MissA2 only)").c_str(),
                       3, 0.5, 3.5);
    h->Sumw2();
    h->GetXaxis()->SetBinLabel(1, "fail p_{T}^{min}");
    h->GetXaxis()->SetBinLabel(2, "fail |#eta|");
    h->GetXaxis()->SetBinLabel(3, "fail #Delta#phi");

    H[name] = h;
    if (prevDir) prevDir->cd();
    return h;
}

TH1F* RecoilJets::getOrBookLeadTruthRecoilMatchMissB_PtGammaTruth(const std::string& trig,
                                                                  const std::string& rKey,
                                                                  int centIdx)
{
    const std::string base   = "h_leadTruthRecoilMatch_missB_pTgammaTruth";
    const std::string suffix = suffixForBins(-1, centIdx);
    const std::string name   = base + "_" + rKey + suffix;

    if (trig.empty() || rKey.empty()) return nullptr;

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

    const std::vector<double>& kPtTruth = m_gammaPtBins;
    const int nb = static_cast<int>(kPtTruth.size()) - 1;

    const std::string title =
      name + ";p_{T}^{#gamma,truth} [GeV];MissB (no reco fid jet within #DeltaR<0.3 of truth jet1)";

    auto* h = new TH1F(name.c_str(), title.c_str(), nb, kPtTruth.data());
    h->Sumw2();

    H[name] = h;
    if (prevDir) prevDir->cd();
    return h;
}

// -------------------------------------------------------------------------
// (SIM ONLY): diagnostics for why reco xJ differs from truth-conditioned xJ
//   Filled inside the same DEN / NUM / MissA / MissB classification block.
// -------------------------------------------------------------------------

TH2F* RecoilJets::getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtTruthLead_num(const std::string& trig,
                                                                           const std::string& rKey,
                                                                           int centIdx)
{
  const std::string base   = "h2_leadTruthRecoilMatch_pTrecoJet1_vs_pTtruthLead_num";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

  std::vector<double> tmpJetPt;
  tmpJetPt.reserve(121);
  for (int i = 0; i <= 120; ++i) tmpJetPt.push_back(0.5 * (double)i);
  const std::vector<double>* pJetPt = &tmpJetPt;

  const int nb = static_cast<int>(pJetPt->size()) - 1;
  auto* h = new TH2F(name.c_str(),
                       (name + ";p_{T}^{lead recoil jet,truth} [GeV];p_{T}^{recoilJet1,reco} [GeV] (NUM)").c_str(),
                       nb, pJetPt->data(),
                       nb, pJetPt->data());
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtTruthLead_missA(const std::string& trig,
                                                                             const std::string& rKey,
                                                                             int centIdx)
{
  const std::string base   = "h2_leadTruthRecoilMatch_pTrecoJet1_vs_pTtruthLead_missA";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

  std::vector<double> tmpJetPt;
  tmpJetPt.reserve(121);
  for (int i = 0; i <= 120; ++i) tmpJetPt.push_back(0.5 * (double)i);
  const std::vector<double>* pJetPt = &tmpJetPt;

  const int nb = static_cast<int>(pJetPt->size()) - 1;
  auto* h = new TH2F(name.c_str(),
                       (name + ";p_{T}^{lead recoil jet,truth} [GeV];p_{T}^{recoilJet1,reco} [GeV] (MissA)").c_str(),
                       nb, pJetPt->data(),
                       nb, pJetPt->data());
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtTruthLead_missA1(const std::string& trig,
                                                                              const std::string& rKey,
                                                                              int centIdx)
{
  const std::string base   = "h2_leadTruthRecoilMatch_pTrecoJet1_vs_pTtruthLead_missA1";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

  std::vector<double> tmpJetPt;
  tmpJetPt.reserve(121);
  for (int i = 0; i <= 120; ++i) tmpJetPt.push_back(0.5 * (double)i);
  const std::vector<double>* pJetPt = &tmpJetPt;

  const int nb = static_cast<int>(pJetPt->size()) - 1;
  auto* h = new TH2F(name.c_str(),
                       (name + ";p_{T}^{lead recoil jet,truth} [GeV];p_{T}^{recoilJet1,reco} [GeV] (MissA1)").c_str(),
                       nb, pJetPt->data(),
                       nb, pJetPt->data());
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtTruthLead_missA2(const std::string& trig,
                                                                              const std::string& rKey,
                                                                              int centIdx)
{
  const std::string base   = "h2_leadTruthRecoilMatch_pTrecoJet1_vs_pTtruthLead_missA2";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

  std::vector<double> tmpJetPt;
  tmpJetPt.reserve(121);
  for (int i = 0; i <= 120; ++i) tmpJetPt.push_back(0.5 * (double)i);
  const std::vector<double>* pJetPt = &tmpJetPt;

  const int nb = static_cast<int>(pJetPt->size()) - 1;
  auto* h = new TH2F(name.c_str(),
                       (name + ";p_{T}^{lead recoil jet,truth} [GeV];p_{T}^{recoilJet1,reco} [GeV] (MissA2)").c_str(),
                       nb, pJetPt->data(),
                       nb, pJetPt->data());
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtTruthLead_missB(const std::string& trig,
                                                                             const std::string& rKey,
                                                                             int centIdx)
{
  const std::string base   = "h2_leadTruthRecoilMatch_pTrecoJet1_vs_pTtruthLead_missB";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

  std::vector<double> tmpJetPt;
  tmpJetPt.reserve(121);
  for (int i = 0; i <= 120; ++i) tmpJetPt.push_back(0.5 * (double)i);
  const std::vector<double>* pJetPt = &tmpJetPt;

  const int nb = static_cast<int>(pJetPt->size()) - 1;
  auto* h = new TH2F(name.c_str(),
                       (name + ";p_{T}^{lead recoil jet,truth} [GeV];p_{T}^{recoilJet1,reco} [GeV] (MissB)").c_str(),
                       nb, pJetPt->data(),
                       nb, pJetPt->data());
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

// (A2) pT(recoilJet1^reco) vs pT(reco jet matched to truth-leading recoil jet), for NUM / MissA

TH2F* RecoilJets::getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtRecoTruthMatch_num(const std::string& trig,
                                                                                const std::string& rKey,
                                                                                int centIdx)
{
  const std::string base   = "h2_leadTruthRecoilMatch_pTrecoJet1_vs_pTrecoTruthMatch_num";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

  std::vector<double> tmpJetPt;
  tmpJetPt.reserve(121);
  for (int i = 0; i <= 120; ++i) tmpJetPt.push_back(0.5 * (double)i);
  const std::vector<double>* pJetPt = &tmpJetPt;

  const int nb = static_cast<int>(pJetPt->size()) - 1;
  auto* h = new TH2F(name.c_str(),
                       (name + ";p_{T}^{reco match to truth-lead} [GeV];p_{T}^{recoilJet1,reco} [GeV] (NUM)").c_str(),
                       nb, pJetPt->data(),
                       nb, pJetPt->data());
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtRecoTruthMatch_missA(const std::string& trig,
                                                                                  const std::string& rKey,
                                                                                  int centIdx)
{
  const std::string base   = "h2_leadTruthRecoilMatch_pTrecoJet1_vs_pTrecoTruthMatch_missA";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

  std::vector<double> tmpJetPt;
  tmpJetPt.reserve(121);
  for (int i = 0; i <= 120; ++i) tmpJetPt.push_back(0.5 * (double)i);
  const std::vector<double>* pJetPt = &tmpJetPt;

  const int nb = static_cast<int>(pJetPt->size()) - 1;
  auto* h = new TH2F(name.c_str(),
                       (name + ";p_{T}^{reco match to truth-lead} [GeV];p_{T}^{recoilJet1,reco} [GeV] (MissA)").c_str(),
                       nb, pJetPt->data(),
                       nb, pJetPt->data());
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtRecoTruthMatch_missA1(const std::string& trig,
                                                                                   const std::string& rKey,
                                                                                   int centIdx)
{
  const std::string base   = "h2_leadTruthRecoilMatch_pTrecoJet1_vs_pTrecoTruthMatch_missA1";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

  std::vector<double> tmpJetPt;
  tmpJetPt.reserve(121);
  for (int i = 0; i <= 120; ++i) tmpJetPt.push_back(0.5 * (double)i);
  const std::vector<double>* pJetPt = &tmpJetPt;

  const int nb = static_cast<int>(pJetPt->size()) - 1;
  auto* h = new TH2F(name.c_str(),
                       (name + ";p_{T}^{reco match to truth-lead} [GeV];p_{T}^{recoilJet1,reco} [GeV] (MissA1)").c_str(),
                       nb, pJetPt->data(),
                       nb, pJetPt->data());
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtRecoTruthMatch_missA2(const std::string& trig,
                                                                                   const std::string& rKey,
                                                                                   int centIdx)
{
  const std::string base   = "h2_leadTruthRecoilMatch_pTrecoJet1_vs_pTrecoTruthMatch_missA2";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

  std::vector<double> tmpJetPt;
  tmpJetPt.reserve(121);
  for (int i = 0; i <= 120; ++i) tmpJetPt.push_back(0.5 * (double)i);
  const std::vector<double>* pJetPt = &tmpJetPt;

  const int nb = static_cast<int>(pJetPt->size()) - 1;
  auto* h = new TH2F(name.c_str(),
                       (name + ";p_{T}^{reco match to truth-lead} [GeV];p_{T}^{recoilJet1,reco} [GeV] (MissA2)").c_str(),
                       nb, pJetPt->data(),
                       nb, pJetPt->data());
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

// (B3) |Δphi(γ^truth, recoilJet1^reco)| vs pTγ,truth, split by class

TH2F* RecoilJets::getOrBookLeadTruthRecoilMatchDphiRecoJet1VsPtGammaTruth_num(const std::string& trig,
                                                                              const std::string& rKey,
                                                                              int centIdx)
{
  const std::string base   = "h2_leadTruthRecoilMatch_dphiRecoJet1_num_pTgammaTruth";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

  const std::vector<double>& kPtTruth = m_gammaPtBins;
  const int nb = static_cast<int>(kPtTruth.size()) - 1;

  auto* h = new TH2F(name.c_str(),
                     (name + ";p_{T}^{#gamma,truth} [GeV];|#Delta#phi(#gamma^{truth}, recoilJet1^{reco})| [rad] (NUM)").c_str(),
                     nb, kPtTruth.data(),
                     64, 0.0, M_PI);
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookLeadTruthRecoilMatchDphiRecoJet1VsPtGammaTruth_missA(const std::string& trig,
                                                                                const std::string& rKey,
                                                                                int centIdx)
{
  const std::string base   = "h2_leadTruthRecoilMatch_dphiRecoJet1_missA_pTgammaTruth";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

    const std::vector<double>& kPtTruth = m_gammaPtBins;
  const int nb = static_cast<int>(kPtTruth.size()) - 1;

  auto* h = new TH2F(name.c_str(),
                     (name + ";p_{T}^{#gamma,truth} [GeV];|#Delta#phi(#gamma^{truth}, recoilJet1^{reco})| [rad] (MissA)").c_str(),
                     nb, kPtTruth.data(),
                     64, 0.0, M_PI);
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookLeadTruthRecoilMatchDphiRecoJet1VsPtGammaTruth_missB(const std::string& trig,
                                                                                const std::string& rKey,
                                                                                int centIdx)
{
  const std::string base   = "h2_leadTruthRecoilMatch_dphiRecoJet1_missB_pTgammaTruth";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

  const std::vector<double>& kPtTruth = m_gammaPtBins;
  const int nb = static_cast<int>(kPtTruth.size()) - 1;

  auto* h = new TH2F(name.c_str(),
                     (name + ";p_{T}^{#gamma,truth} [GeV];|#Delta#phi(#gamma^{truth}, recoilJet1^{reco})| [rad] (MissB)").c_str(),
                     nb, kPtTruth.data(),
                     64, 0.0, M_PI);
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

// (B4) ΔR(recoilJet1^reco, truth-leading recoil jet) vs pTγ,truth, split by class

TH2F* RecoilJets::getOrBookLeadTruthRecoilMatchDRRecoJet1VsPtGammaTruth_num(const std::string& trig,
                                                                            const std::string& rKey,
                                                                            int centIdx)
{
  const std::string base   = "h2_leadTruthRecoilMatch_dRRecoJet1_vs_truthLead_num_pTgammaTruth";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

  const std::vector<double>& kPtTruth = m_gammaPtBins;
  const int nb = static_cast<int>(kPtTruth.size()) - 1;

  auto* h = new TH2F(name.c_str(),
                     (name + ";p_{T}^{#gamma,truth} [GeV];#DeltaR(recoilJet1^{reco}, truth lead recoil) (NUM)").c_str(),
                     nb, kPtTruth.data(),
                     70, 0.0, 3.5);
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookLeadTruthRecoilMatchDRRecoJet1VsPtGammaTruth_missA(const std::string& trig,
                                                                              const std::string& rKey,
                                                                              int centIdx)
{
  const std::string base   = "h2_leadTruthRecoilMatch_dRRecoJet1_vs_truthLead_missA_pTgammaTruth";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

  const std::vector<double>& kPtTruth = m_gammaPtBins;
  const int nb = static_cast<int>(kPtTruth.size()) - 1;

  auto* h = new TH2F(name.c_str(),
                     (name + ";p_{T}^{#gamma,truth} [GeV];#DeltaR(recoilJet1^{reco}, truth lead recoil) (MissA)").c_str(),
                     nb, kPtTruth.data(),
                     70, 0.0, 3.5);
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookLeadTruthRecoilMatchDRRecoJet1VsPtGammaTruth_missB(const std::string& trig,
                                                                              const std::string& rKey,
                                                                              int centIdx)
{
  const std::string base   = "h2_leadTruthRecoilMatch_dRRecoJet1_vs_truthLead_missB_pTgammaTruth";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

  const std::vector<double>& kPtTruth = m_gammaPtBins;
  const int nb = static_cast<int>(kPtTruth.size()) - 1;

  auto* h = new TH2F(name.c_str(),
                     (name + ";p_{T}^{#gamma,truth} [GeV];#DeltaR(recoilJet1^{reco}, truth lead recoil) (MissB)").c_str(),
                     nb, kPtTruth.data(),
                     70, 0.0, 3.5);
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

// (C5) xJ(recoilJet1^reco) vs |Δphi(γ^truth, recoilJet1^reco)|, split by class

TH2F* RecoilJets::getOrBookLeadTruthRecoilMatchXJRecoJet1VsDphiRecoJet1_num(const std::string& trig,
                                                                            const std::string& rKey,
                                                                            int centIdx)
{
  const std::string base   = "h2_leadTruthRecoilMatch_xJRecoJet1_vs_dphiRecoJet1_num";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

  auto* h = new TH2F(name.c_str(),
                     (name + ";|#Delta#phi(#gamma^{truth}, recoilJet1^{reco})| [rad];x_{J}^{reco}=p_{T}^{recoilJet1}/p_{T}^{#gamma} (NUM)").c_str(),
                     64, 0.0, M_PI,
                     60, 0.0, 3.0);
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookLeadTruthRecoilMatchXJRecoJet1VsDphiRecoJet1_missA(const std::string& trig,
                                                                              const std::string& rKey,
                                                                              int centIdx)
{
  const std::string base   = "h2_leadTruthRecoilMatch_xJRecoJet1_vs_dphiRecoJet1_missA";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

  auto* h = new TH2F(name.c_str(),
                     (name + ";|#Delta#phi(#gamma^{truth}, recoilJet1^{reco})| [rad];x_{J}^{reco}=p_{T}^{recoilJet1}/p_{T}^{#gamma} (MissA)").c_str(),
                     64, 0.0, M_PI,
                     60, 0.0, 3.0);
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookLeadTruthRecoilMatchXJRecoJet1VsDphiRecoJet1_missB(const std::string& trig,
                                                                              const std::string& rKey,
                                                                              int centIdx)
{
  const std::string base   = "h2_leadTruthRecoilMatch_xJRecoJet1_vs_dphiRecoJet1_missB";
  const std::string suffix = suffixForBins(-1, centIdx);
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

  auto* h = new TH2F(name.c_str(),
                     (name + ";|#Delta#phi(#gamma^{truth}, recoilJet1^{reco})| [rad];x_{J}^{reco}=p_{T}^{recoilJet1}/p_{T}^{#gamma} (MissB)").c_str(),
                     64, 0.0, M_PI,
                     60, 0.0, 3.0);
  h->Sumw2();

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}

// ------------------------------------------------------------------
// NEW / REQUIRED: TH3 for xJ vs alpha vs pTgamma (radius-tagged; centrality suffix only)
//   name pattern: h_JES3_pT_xJ_alpha_<rKey><centSuffix>
// ------------------------------------------------------------------
TH3F* RecoilJets::getOrBookJES3_xJ_alphaHist(const std::string& trig,
                                             const std::string& rKey,
                                             int centIdx)
{
  const std::string base   = "h_JES3_pT_xJ_alpha";
  const std::string suffix = suffixForBins(-1, centIdx); // cent-only (no pT suffix)
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

  auto& H = qaHistogramsByTrigger[trig];
  if (auto it = H.find(name); it != H.end())
  {
    if (auto* h = dynamic_cast<TH3F*>(it->second)) return h;
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

  const int ny = 60;
  const double ylo = 0.0, yhi = 3.0;   // xJ

  const int nz = 40;
  const double zlo = 0.0, zhi = 2.0;   // alpha

  const std::string title =
    name + ";p_{T}^{#gamma} [GeV];x_{J}=p_{T}^{jet1}/p_{T}^{#gamma};#alpha=p_{T}^{jet2}/p_{T}^{#gamma}";

  std::vector<double> ybins(ny + 1);
  for (int i = 0; i <= ny; ++i)
    ybins[i] = ylo + (yhi - ylo) * (static_cast<double>(i) / ny);

  std::vector<double> zbins(nz + 1);
  for (int i = 0; i <= nz; ++i)
    zbins[i] = zlo + (zhi - zlo) * (static_cast<double>(i) / nz);

  auto* h = new TH3F(name.c_str(), title.c_str(),
                     nx, xbins,
                     ny, ybins.data(),
                     nz, zbins.data());

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}


// ------------------------------------------------------------------
// NEW / REQUIRED: TH3 for jet1Pt vs alpha vs pTgamma (radius-tagged; centrality suffix only)
//   name pattern: h_JES3_pT_jet1Pt_alpha_<rKey><centSuffix>
// ------------------------------------------------------------------
TH3F* RecoilJets::getOrBookJES3_jet1Pt_alphaHist(const std::string& trig,
                                                const std::string& rKey,
                                                int centIdx)
{
  const std::string base   = "h_JES3_pT_jet1Pt_alpha";
  const std::string suffix = suffixForBins(-1, centIdx); // cent-only
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

  auto& H = qaHistogramsByTrigger[trig];
  if (auto it = H.find(name); it != H.end())
  {
    if (auto* h = dynamic_cast<TH3F*>(it->second)) return h;
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

  const int ny = 120;
  const double ylo = 0.0, yhi = 60.0;  // jet1Pt

  const int nz = 40;
  const double zlo = 0.0, zhi = 2.0;   // alpha

  const std::string title =
    name + ";p_{T}^{#gamma} [GeV];p_{T}^{jet1} [GeV];#alpha=p_{T}^{jet2}/p_{T}^{#gamma}";

  std::vector<double> ybins(ny + 1);
  for (int i = 0; i <= ny; ++i)
    ybins[i] = ylo + (yhi - ylo) * (static_cast<double>(i) / ny);

  std::vector<double> zbins(nz + 1);
  for (int i = 0; i <= nz; ++i)
    zbins[i] = zlo + (zhi - zlo) * (static_cast<double>(i) / nz);

  auto* h = new TH3F(name.c_str(), title.c_str(),
                     nx, xbins,
                     ny, ybins.data(),
                     nz, zbins.data());

  H[name] = h;
  if (prevDir) prevDir->cd();
  return h;
}


// ============================================================================
//  Matching-QA bookers (centrality suffix only; pT^gamma is an axis)
// ============================================================================

TH2F* RecoilJets::getOrBookMatchDphiVsPtGamma(const std::string& trig,
                                              const std::string& rKey,
                                              int centIdx)
{
  const std::string base   = "h_match_dphi_vs_pTgamma";
  const std::string suffix = suffixForBins(-1, centIdx);   // centrality-only
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

TH2F* RecoilJets::getOrBookMatchStatusVsPtGamma(const std::string& trig,
                                                const std::string& rKey,
                                                int centIdx)
{
  const std::string base   = "h_match_status_vs_pTgamma";
  const std::string suffix = suffixForBins(-1, centIdx);   // centrality-only
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

TH2F* RecoilJets::getOrBookJetCutflowStatusVsPtGamma(const std::string& trig,
                                                     const std::string& rKey,
                                                     int centIdx)
{
  const std::string base   = "h_jetcutflow_status_vs_pTgamma";
  const std::string suffix = suffixForBins(-1, centIdx);   // centrality-only
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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
    name + ";p_{T}^{#gamma} [GeV];jet cutflow status";

  auto* h = new TH2F(name.c_str(), title.c_str(), nx, xbins, ny, ylo, yhi);

  h->GetYaxis()->SetBinLabel(1, "FailJetPt");
  h->GetYaxis()->SetBinLabel(2, "FailJetEta");
  h->GetYaxis()->SetBinLabel(3, "FailBackToBack");
  h->GetYaxis()->SetBinLabel(4, "PassAll");

  H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}

TH2F* RecoilJets::getOrBookMatchMaxDphiVsPtGamma(const std::string& trig,
                                                 const std::string& rKey,
                                                 int centIdx)
{
  const std::string base   = "h_match_maxdphi_vs_pTgamma";
  const std::string suffix = suffixForBins(-1, centIdx);   // centrality-only
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

TProfile* RecoilJets::getOrBookNRecoilJetsVsPtGamma(const std::string& trig,
                                                    const std::string& rKey,
                                                    int centIdx)
{
  const std::string base   = "p_nRecoilJets_vs_pTgamma";
  const std::string suffix = suffixForBins(-1, centIdx);   // centrality-only
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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

TH2F* RecoilJets::getOrBookRecoilIsLeadingVsPtGamma(const std::string& trig,
                                                    const std::string& rKey,
                                                    int centIdx)
{
  const std::string base   = "h_recoilIsLeading_vs_pTgamma";
  const std::string suffix = suffixForBins(-1, centIdx);   // centrality-only
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

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
//  NEW 3D baseline bookers (NO centrality suffix)
// ============================================================================

TH3F* RecoilJets::getOrBookPho3TightIso(const std::string& trig)
{
  const std::string name = "h_Pho3_pT_eta_phi_tightIso";
  if (trig.empty()) return nullptr;

  auto& H = qaHistogramsByTrigger[trig];
  if (auto it = H.find(name); it != H.end())
  {
    if (auto* h = dynamic_cast<TH3F*>(it->second)) return h;
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

  const int ny = 56;
  const double ylo = -0.7, yhi = 0.7;

  const int nz = 64;
  const double zlo = -M_PI, zhi = M_PI;

    const std::string title = name + ";p_{T}^{#gamma} [GeV];#eta^{#gamma};#phi^{#gamma} [rad]";

    // TH3F does not support mixed (variable x bins, uniform y/z) constructors in ROOT 6.32.
    // Build uniform y/z bin edges explicitly and use the fully-variable-bins constructor.
    std::vector<double> ybins(ny + 1);
    for (int i = 0; i <= ny; ++i)
    {
      ybins[i] = ylo + (yhi - ylo) * (static_cast<double>(i) / static_cast<double>(ny));
    }

    std::vector<double> zbins(nz + 1);
    for (int i = 0; i <= nz; ++i)
    {
      zbins[i] = zlo + (zhi - zlo) * (static_cast<double>(i) / static_cast<double>(nz));
    }

    auto* h = new TH3F(name.c_str(), title.c_str(),
                       nx, xbins,
                       ny, ybins.data(),
                       nz, zbins.data());
    H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}

TH3F* RecoilJets::getOrBookJet13RecoilJet1(const std::string& trig,
                                          const std::string& rKey)
{
  const std::string base = "h_Jet13_pTgamma_eta_phi_recoilJet1";
  const std::string name = base + "_" + rKey;

  if (trig.empty() || rKey.empty()) return nullptr;

  auto& H = qaHistogramsByTrigger[trig];
  if (auto it = H.find(name); it != H.end())
  {
    if (auto* h = dynamic_cast<TH3F*>(it->second)) return h;
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

  const int ny = 60;
  const double ylo = -1.1, yhi = 1.1;

  const int nz = 64;
  const double zlo = -M_PI, zhi = M_PI;

    const std::string title = name + ";p_{T}^{#gamma} [GeV];#eta^{jet1};#phi^{jet1} [rad]";

    // TH3F does not support mixed (variable x bins, uniform y/z) constructors in ROOT 6.32.
    // Build uniform y/z bin edges explicitly and use the fully-variable-bins constructor.
    std::vector<double> ybins(ny + 1);
    for (int i = 0; i <= ny; ++i)
    {
      ybins[i] = ylo + (yhi - ylo) * (static_cast<double>(i) / static_cast<double>(ny));
    }

    std::vector<double> zbins(nz + 1);
    for (int i = 0; i <= nz; ++i)
    {
      zbins[i] = zlo + (zhi - zlo) * (static_cast<double>(i) / static_cast<double>(nz));
    }

    auto* h = new TH3F(name.c_str(), title.c_str(),
                       nx, xbins,
                       ny, ybins.data(),
                       nz, zbins.data());
    H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}

TProfile3D* RecoilJets::getOrBookBalance3(const std::string& trig,
                                          const std::string& rKey)
{
  const std::string base = "p_Balance3_pTgamma_eta_phi";
  const std::string name = base + "_" + rKey;

  if (trig.empty() || rKey.empty()) return nullptr;

  auto& H = qaHistogramsByTrigger[trig];
  if (auto it = H.find(name); it != H.end())
  {
    if (auto* p = dynamic_cast<TProfile3D*>(it->second)) return p;
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

  const int ny = 60;
  const double ylo = -1.1, yhi = 1.1;

  const int nz = 64;
  const double zlo = -M_PI, zhi = M_PI;

  const std::string title =
    name + ";p_{T}^{#gamma} [GeV];#eta^{jet1};#phi^{jet1} [rad];#LT x_{J} #GT";

    // TProfile3D does not support mixed (variable x bins, uniform y/z) constructors in ROOT 6.32.
    // Build uniform y/z bin edges explicitly and use the fully-variable-bins constructor.
    std::vector<double> ybins(ny + 1);
    for (int i = 0; i <= ny; ++i)
    {
      ybins[i] = ylo + (yhi - ylo) * (static_cast<double>(i) / static_cast<double>(ny));
    }

    std::vector<double> zbins(nz + 1);
    for (int i = 0; i <= nz; ++i)
    {
      zbins[i] = zlo + (zhi - zlo) * (static_cast<double>(i) / static_cast<double>(nz));
    }

    auto* p = new TProfile3D(name.c_str(), title.c_str(),
                             nx, xbins,
                             ny, ybins.data(),
                             nz, zbins.data(),
                             "");
    H[name] = p;

  if (prevDir) prevDir->cd();
  return p;
}


// ============================================================================
//  SIM-only truth JES3 bookers (centrality suffix only; pTgamma is axis)
// ============================================================================

TH3F* RecoilJets::getOrBookJES3Truth_xJ_alphaHist(const std::string& trig,
                                                 const std::string& rKey,
                                                 int centIdx)
{
  const std::string base   = "h_JES3Truth_pT_xJ_alpha";
  const std::string suffix = suffixForBins(-1, centIdx);   // cent-only
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

  auto& H = qaHistogramsByTrigger[trig];
  if (auto it = H.find(name); it != H.end())
  {
    if (auto* h = dynamic_cast<TH3F*>(it->second)) return h;
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

  const int ny = 60;
  const double ylo = 0.0, yhi = 3.0; // xJ

  const int nz = 40;
  const double zlo = 0.0, zhi = 2.0; // alpha

  std::vector<double> ybins(ny + 1);
  for (int i = 0; i <= ny; ++i)
    ybins[i] = ylo + (yhi - ylo) * (static_cast<double>(i) / ny);

  std::vector<double> zbins(nz + 1);
  for (int i = 0; i <= nz; ++i)
    zbins[i] = zlo + (zhi - zlo) * (static_cast<double>(i) / nz);

  const std::string title =
    name + ";p_{T}^{#gamma,truth} [GeV];x_{J}^{truth}=p_{T}^{jet1,truth}/p_{T}^{#gamma,truth};#alpha^{truth}=p_{T}^{jet2,truth}/p_{T}^{#gamma,truth}";

  auto* h = new TH3F(name.c_str(), title.c_str(),
                     nx, xbins,
                     ny, ybins.data(),
                     nz, zbins.data());
  H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}

// -----------------------------------------------------------------------------
// NEW (SIM-only): TRUTH reco-conditioned but NO reco↔truth jet match
// Name: h_JES3TruthRecoCondNoJetMatch_pT_xJ_alpha_<rKey><centSuffix>
// -----------------------------------------------------------------------------
TH3F* RecoilJets::getOrBookJES3TruthRecoCondNoJetMatch_xJ_alphaHist(const std::string& trig,
                                                                    const std::string& rKey,
                                                                    int centIdx)
{
  const std::string base   = "h_JES3TruthRecoCondNoJetMatch_pT_xJ_alpha";
  const std::string suffix = suffixForBins(-1, centIdx);   // cent-only
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

  auto& H = qaHistogramsByTrigger[trig];
  if (auto it = H.find(name); it != H.end())
  {
    if (auto* h = dynamic_cast<TH3F*>(it->second)) return h;
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

  const int ny = 60;
  const double ylo = 0.0, yhi = 3.0; // xJ
  const int nz = 40;
  const double zlo = 0.0, zhi = 2.0; // alpha

  std::vector<double> ybins(ny + 1);
  for (int i = 0; i <= ny; ++i) ybins[i] = ylo + (yhi - ylo) * (static_cast<double>(i) / ny);

  std::vector<double> zbins(nz + 1);
  for (int i = 0; i <= nz; ++i) zbins[i] = zlo + (zhi - zlo) * (static_cast<double>(i) / nz);

  const std::string title =
    name + ";p_{T}^{#gamma,truth} [GeV];x_{J}^{truth}=p_{T}^{jet1,truth}/p_{T}^{#gamma,truth};#alpha^{truth}=p_{T}^{jet2,truth}/p_{T}^{#gamma,truth}";

  auto* h = new TH3F(name.c_str(), title.c_str(),
                     nx, xbins,
                     ny, ybins.data(),
                     nz, zbins.data());
  H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}

// -----------------------------------------------------------------------------
// NEW (SIM-only): RECO xJgamma for truth-tagged PHOTON events (no truth-jet match required)
// Name: h_JES3RecoTruthPhoTagged_pT_xJ_alpha_<rKey><centSuffix>
// -----------------------------------------------------------------------------
TH3F* RecoilJets::getOrBookJES3RecoTruthPhoTagged_xJ_alphaHist(const std::string& trig,
                                                               const std::string& rKey,
                                                               int centIdx)
{
  const std::string base   = "h_JES3RecoTruthPhoTagged_pT_xJ_alpha";
  const std::string suffix = suffixForBins(-1, centIdx);   // cent-only
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

  auto& H = qaHistogramsByTrigger[trig];
  if (auto it = H.find(name); it != H.end())
  {
    if (auto* h = dynamic_cast<TH3F*>(it->second)) return h;
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

  const int ny = 60;
  const double ylo = 0.0, yhi = 3.0; // xJ
  const int nz = 40;
  const double zlo = 0.0, zhi = 2.0; // alpha

  std::vector<double> ybins(ny + 1);
  for (int i = 0; i <= ny; ++i) ybins[i] = ylo + (yhi - ylo) * (static_cast<double>(i) / ny);

  std::vector<double> zbins(nz + 1);
  for (int i = 0; i <= nz; ++i) zbins[i] = zlo + (zhi - zlo) * (static_cast<double>(i) / nz);

  const std::string title =
    name + ";p_{T}^{#gamma,reco} [GeV];x_{J}^{reco}=p_{T}^{jet1}/p_{T}^{#gamma};#alpha^{reco}=p_{T}^{jet2}/p_{T}^{#gamma}";

  auto* h = new TH3F(name.c_str(), title.c_str(),
                     nx, xbins,
                     ny, ybins.data(),
                     nz, zbins.data());
  H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}

// -----------------------------------------------------------------------------
// NEW (SIM-only): RECO xJgamma for truth-tagged pairs (photon truth-signal + reco jet1 matched to truth jet1)
// Name: h_JES3RecoTruthTagged_pT_xJ_alpha_<rKey><centSuffix>
// -----------------------------------------------------------------------------
TH3F* RecoilJets::getOrBookJES3RecoTruthTagged_xJ_alphaHist(const std::string& trig,
                                                            const std::string& rKey,
                                                            int centIdx)
{
  const std::string base   = "h_JES3RecoTruthTagged_pT_xJ_alpha";
  const std::string suffix = suffixForBins(-1, centIdx);   // cent-only
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

  auto& H = qaHistogramsByTrigger[trig];
  if (auto it = H.find(name); it != H.end())
  {
    if (auto* h = dynamic_cast<TH3F*>(it->second)) return h;
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

  const int ny = 60;
  const double ylo = 0.0, yhi = 3.0; // xJ
  const int nz = 40;
  const double zlo = 0.0, zhi = 2.0; // alpha

  std::vector<double> ybins(ny + 1);
  for (int i = 0; i <= ny; ++i) ybins[i] = ylo + (yhi - ylo) * (static_cast<double>(i) / ny);

  std::vector<double> zbins(nz + 1);
  for (int i = 0; i <= nz; ++i) zbins[i] = zlo + (zhi - zlo) * (static_cast<double>(i) / nz);

  const std::string title =
    name + ";p_{T}^{#gamma,reco} [GeV];x_{J}^{reco}=p_{T}^{jet1}/p_{T}^{#gamma};#alpha^{reco}=p_{T}^{jet2}/p_{T}^{#gamma}";

  auto* h = new TH3F(name.c_str(), title.c_str(),
                     nx, xbins,
                     ny, ybins.data(),
                     nz, zbins.data());
  H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}

// -----------------------------------------------------------------------------
// PURE truth xJgamma distribution booker (NO reco gating, NO reco↔truth jet match)
// Name: h_JES3TruthPure_pT_xJ_alpha_<rKey><centSuffix>
// -----------------------------------------------------------------------------
TH3F* RecoilJets::getOrBookJES3TruthPure_xJ_alphaHist(const std::string& trig,
                                                     const std::string& rKey,
                                                     int centIdx)
{
  const std::string base   = "h_JES3TruthPure_pT_xJ_alpha";
  const std::string suffix = suffixForBins(-1, centIdx);   // cent-only
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

  auto& H = qaHistogramsByTrigger[trig];
  if (auto it = H.find(name); it != H.end())
  {
    if (auto* h = dynamic_cast<TH3F*>(it->second)) return h;
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

  const int ny = 60;
  const double ylo = 0.0, yhi = 3.0; // xJ

  const int nz = 40;
  const double zlo = 0.0, zhi = 2.0; // alpha

  std::vector<double> ybins(ny + 1);
  for (int i = 0; i <= ny; ++i)
    ybins[i] = ylo + (yhi - ylo) * (static_cast<double>(i) / ny);

  std::vector<double> zbins(nz + 1);
  for (int i = 0; i <= nz; ++i)
    zbins[i] = zlo + (zhi - zlo) * (static_cast<double>(i) / nz);

  const std::string title =
    name + ";p_{T}^{#gamma,truth} [GeV];x_{J}^{truth}=p_{T}^{jet1,truth}/p_{T}^{#gamma,truth};#alpha^{truth}=p_{T}^{jet2,truth}/p_{T}^{#gamma,truth}";

  auto* h = new TH3F(name.c_str(), title.c_str(),
                     nx, xbins,
                     ny, ybins.data(),
                     nz, zbins.data());
  H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}

TH3F* RecoilJets::getOrBookJES3Truth_jet1Pt_alphaHist(const std::string& trig,
                                                     const std::string& rKey,
                                                     int centIdx)
{
  const std::string base   = "h_JES3Truth_pT_jet1Pt_alpha";
  const std::string suffix = suffixForBins(-1, centIdx);   // cent-only
  const std::string name   = base + "_" + rKey + suffix;

  if (trig.empty() || rKey.empty()) return nullptr;

  auto& H = qaHistogramsByTrigger[trig];
  if (auto it = H.find(name); it != H.end())
  {
    if (auto* h = dynamic_cast<TH3F*>(it->second)) return h;
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

  const int ny = 120;
  const double ylo = 0.0, yhi = 60.0; // jet1Pt

  const int nz = 40;
  const double zlo = 0.0, zhi = 2.0; // alpha

  std::vector<double> ybins(ny + 1);
  for (int i = 0; i <= ny; ++i)
    ybins[i] = ylo + (yhi - ylo) * (static_cast<double>(i) / ny);

  std::vector<double> zbins(nz + 1);
  for (int i = 0; i <= nz; ++i)
    zbins[i] = zlo + (zhi - zlo) * (static_cast<double>(i) / nz);

  const std::string title =
    name + ";p_{T}^{#gamma,truth} [GeV];p_{T}^{jet1,truth} [GeV];#alpha^{truth}=p_{T}^{jet2,truth}/p_{T}^{#gamma,truth}";

  auto* h = new TH3F(name.c_str(), title.c_str(),
                     nx, xbins,
                     ny, ybins.data(),
                     nz, zbins.data());
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
//  isolation PASS/FAIL counter histogram (2 bins), same slicing rules.
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

// Event-level reco-photon multiplicity diagnostic:
//   N = number of reco photon candidates passing the SAME iso∧tight gate used
//       for the event-leading photon selection (filled once per event; binned
//       by the selected leading-photon pT bin via suffixForBins()).
TH1I* RecoilJets::getOrBookNIsoTightPhoCandHist(const std::string& trig, int ptIdx, int centIdx)
{
  const std::string base   = "h_nIsoTightPhoCand";
  const std::string suffix = suffixForBins(ptIdx, centIdx);
  const std::string name   = base + suffix;

  if (Verbosity() >= 5)
    LOG(5, CLR_BLUE, "  [getOrBookNIsoTightPhoCandHist] trig=\"" << trig
           << "\" ptIdx=" << ptIdx << " centIdx=" << centIdx
           << " → name=\"" << name << "\"");

  if (trig.empty())
  {
    LOG(2, CLR_YELLOW, "  [getOrBookNIsoTightPhoCandHist] empty trig – returning nullptr");
    return nullptr;
  }

  auto& H = qaHistogramsByTrigger[trig];

  if (auto it = H.find(name); it != H.end())
  {
    if (auto* h = dynamic_cast<TH1I*>(it->second)) return h;
    LOG(2, CLR_YELLOW, "    [getOrBookNIsoTightPhoCandHist] name clash: object \"" << name
                       << "\" exists but is not TH1I – replacing it");
    H.erase(it);
  }

  if (!out || !out->IsOpen())
  {
    LOG(1, CLR_YELLOW, "  [getOrBookNIsoTightPhoCandHist] output TFile invalid/null – returning nullptr");
    return nullptr;
  }

  TDirectory* const prevDir = gDirectory;
  TDirectory* dir = out->GetDirectory(trig.c_str());
  if (!dir) dir = out->mkdir(trig.c_str());
  if (!dir)
  {
    LOG(1, CLR_YELLOW, "  [getOrBookNIsoTightPhoCandHist] failed to create/access directory \"" << trig << "\"");
    if (prevDir) prevDir->cd();
    return nullptr;
  }

  dir->cd();

  // Integer multiplicity: 0..10
  auto* h = new TH1I(name.c_str(),
                     (name + ";N_{#gamma}^{iso+tight} candidates;Entries").c_str(),
                     11, -0.5, 10.5);
  if (!h)
  {
    LOG(1, CLR_YELLOW, "  [getOrBookNIsoTightPhoCandHist] new TH1I failed for \"" << name << '"');
    if (prevDir) prevDir->cd();
    return nullptr;
  }

  H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}


// ============================================================================
//  SIM-only truth isolation QA bookers
//    - stored directly in the trigger directory (SIM => /SIM/)
//    - NOT sliced by pT/centrality (these are global truth QA histograms)
// ============================================================================
TH1F* RecoilJets::getOrBookTruthIsoHist(const std::string& trig,
                                       const std::string& name,
                                       int nbins, double xmin, double xmax)
{
  if (trig.empty() || name.empty()) return nullptr;

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

  const std::string title = name + ";E_{T}^{iso,truth} [GeV];Entries";
  auto* h = new TH1F(name.c_str(), title.c_str(), nbins, xmin, xmax);
  h->Sumw2();
  H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}

TH1I* RecoilJets::getOrBookTruthIsoDecisionHist(const std::string& trig,
                                               const std::string& name)
{
  if (trig.empty() || name.empty()) return nullptr;

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

  const std::string title = name + ";Truth isolation cut decision;Entries";
  auto* h = new TH1I(name.c_str(), title.c_str(), 2, 0.5, 2.5);
  h->GetXaxis()->SetBinLabel(1, "PASS");
  h->GetXaxis()->SetBinLabel(2, "FAIL");
  H[name] = h;

  if (prevDir) prevDir->cd();
  return h;
}

// ------------------------------------------------------------------
//  matched truth-signal → reco ABCD leakage histogram
//   - one histogram per (pT[/cent]) slice
//   - bins: 1=A, 2=B, 3=C, 4=D
// ------------------------------------------------------------------
TH1I* RecoilJets::getOrBookSigABCDLeakageHist(const std::string& trig, int ptIdx, int centIdx)
{
    const std::string base   = "h_sigABCD_MC";
    const std::string suffix = suffixForBins(ptIdx, centIdx);
    const std::string name   = base + suffix;

    if (trig.empty()) return nullptr;

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

    const std::string title =
      name + ";Reco ABCD region for matched truth-signal photons;Entries";

    auto* h = new TH1I(name.c_str(), title.c_str(), 4, 0.5, 4.5);
    h->GetXaxis()->SetBinLabel(1, "A");
    h->GetXaxis()->SetBinLabel(2, "B");
    h->GetXaxis()->SetBinLabel(3, "C");
    h->GetXaxis()->SetBinLabel(4, "D");

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

  // Fiducial/containment cut depends on jet radius key:
  //   |eta_jet| < 1.1 - R
  const double etaAbsMaxUse = jetEtaAbsMaxForRKey(rKey);

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
    // Fiducial/containment cut depends on jet radius key:
    //   |eta_jet| < 1.1 - R
    if (std::fabs(eta) >= etaAbsMaxUse) continue;

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

    // Per-event jet summary (printed once per rKey, not once per trigger)
    if (Verbosity() >= 6)
    {
      std::ostringstream os;
      os << "    [jetQA:event] rKey=" << rKey
         << " | jets(container)=" << jets->size()
         << " | fiducial jets (pT>" << m_minJetPt
         << ", |eta|<" << std::fixed << std::setprecision(2) << etaAbsMaxUse
         << ")=" << nJetsFid
         << " | HT=" << std::fixed << std::setprecision(2) << HT;

      if (leadPt >= 0.0)
      {
        os << " | lead=(pT=" << std::fixed << std::setprecision(2) << leadPt
           << ", eta=" << std::fixed << std::setprecision(3) << leadEta
           << ", phi=" << std::fixed << std::setprecision(3) << leadPhi << ")";
      }
      if (subPt >= 0.0)
      {
        os << " | sublead pT=" << std::fixed << std::setprecision(2) << subPt;
      }

      LOG(6, CLR_CYAN, os.str());
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

  // ------------------------------------------------------------------
  // (A) Job-wide fill counter (used by your End() summary)
  // ------------------------------------------------------------------
  const std::string key = trig + "::" + hnameWithSuffix;

  bool firstThisHist = false;
  long long jobCount = 0;

  auto it = m_histFill.find(key);
  if (it == m_histFill.end())
  {
    m_histFill.emplace(key, 1);
    firstThisHist = true;
    jobCount = 1;

    // Keep your old "first fill" behavior (useful at medium verbosity)
    if (Verbosity() >= 6)
      LOG(6, CLR_BLUE, "    [bumpHistFill] first fill → \"" << key << "\" = 1");
  }
  else
  {
    ++(it->second);
    jobCount = it->second;
  }

  // ------------------------------------------------------------------
  // (B) Per-fill logging (ONLY at high verbosity), rate-limited per event
  //
  //  Verbosity() >= 7: prints fill lines
  //  Verbosity() >= 8: larger budget (more lines per event)
  //  Verbosity() >= 9: essentially unlimited
  //
  //  Optional filter:
  //    export RJ_HIST_VERBOSE_REGEX='h_ss_|h_Eiso|h_xJ'
  //  ------------------------------------------------------------------
  const int v = static_cast<int>(Verbosity());
  if (v < 7) return;

  // ---- optional regex filter (match on "trig::hist") -----------------
  static bool cfgInit = false;
  static bool useFilter = false;
  static std::regex filter;

  if (!cfgInit)
  {
    cfgInit = true;
    if (gSystem)
    {
      const char* re = gSystem->Getenv("RJ_HIST_VERBOSE_REGEX");
      if (re && std::string(re).size() > 0)
      {
        try
        {
          filter = std::regex(re);
          useFilter = true;
          std::cout << CLR_CYAN
                    << "[RecoilJets] RJ_HIST_VERBOSE_REGEX enabled: \"" << re << "\""
                    << CLR_RESET << std::endl;
        }
        catch (...)
        {
          useFilter = false;
          std::cout << CLR_YELLOW
                    << "[RecoilJets] RJ_HIST_VERBOSE_REGEX is invalid — ignoring filter"
                    << CLR_RESET << std::endl;
        }
      }
    }
  }

  if (useFilter)
  {
    if (!std::regex_search(key, filter) && !std::regex_search(hnameWithSuffix, filter))
      return;
  }

  // ---- per-event print budget ---------------------------------------
  const long long evt = static_cast<long long>(event_count);

  static long long lastEvt = -1;
  static int linesThisEvt = 0;
  static bool warnedThisEvt = false;

  if (evt != lastEvt)
  {
    lastEvt = evt;
    linesThisEvt = 0;
    warnedThisEvt = false;
  }

  int cap = 200;
  if (v >= 9) cap = 100000000;   // effectively unlimited
  else if (v >= 8) cap = 2000;   // heavy but still bounded
  else             cap = 200;    // safe default

  if (linesThisEvt >= cap)
  {
    if (!warnedThisEvt)
    {
      warnedThisEvt = true;
      LOG(6, CLR_YELLOW,
          "    [HFill] (evt=" << evt << ") print budget reached (cap=" << cap
          << "). Suppressing further histogram-fill prints this event. "
          << "Raise Verbosity() to >=9 for unlimited, or set RJ_HIST_VERBOSE_REGEX to narrow output.");
    }
    return;
  }
  ++linesThisEvt;

  // ---- resolve TObject for type + entries (best-effort) --------------
  TObject* obj = nullptr;
  auto itTrig = qaHistogramsByTrigger.find(trig);
  if (itTrig != qaHistogramsByTrigger.end())
  {
    auto itObj = itTrig->second.find(hnameWithSuffix);
    if (itObj != itTrig->second.end()) obj = itObj->second;
  }

  const char* cls = obj ? obj->ClassName() : "unresolved";
  long long entries = -1;
  if (auto* h = dynamic_cast<TH1*>(obj)) entries = static_cast<long long>(h->GetEntries());

  std::ostringstream os;
  os << "    [HFill] evt=" << evt
     << " trig=\"" << trig << "\""
     << " hist=\"" << hnameWithSuffix << "\""
     << " type=" << cls;

  if (entries >= 0) os << " entries=" << entries;
  os << " jobCount=" << jobCount;

  if (firstThisHist) os << " (first-time)";

  const char* col = firstThisHist ? CLR_GREEN : CLR_CYAN;
  LOG(7, col, os.str());
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
                                      double pt_gamma,
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
      LOG(4, CLR_YELLOW,
          "  [fillIsoSSTagCounters] pT bin not found for pT^γ=" << pt_gamma
          << " – skipping fills");
    return;
  }

  const int effCentIdx = (m_isAuAu ? centIdx : -1);
  const std::string slice = suffixForBins(ptIdx, effCentIdx);

  // -------------------------------------------------------------------------
  // Isolation region definition for purity (PPG12)
  //
  //   ISO region:      Eiso < (A + B * pT)
  //   NONISO sideband: Eiso > (A + B * pT + gap)   with gap = +1 GeV
  //   GAP region:      [thrIso, thrIso+gap] is excluded from A–B–C–D
  // -------------------------------------------------------------------------
  const double eiso_et = eiso(clus, topNode);

  // If PhotonClusterBuilder iso_* is missing/non-finite, your eiso() returns a huge fail-safe.
  // Do NOT treat that as "non-isolated sideband" — skip purity classification entirely.
  if (!std::isfinite(eiso_et) || eiso_et > 1e8)
  {
    if (Verbosity() >= 4)
      LOG(4, CLR_YELLOW,
          "  [fillIsoSSTagCounters] invalid/missing Eiso=" << eiso_et
          << " → skip A–B–C–D fill (purity regions require valid iso)");
    return;
  }

  const double thrIso    = (m_isoA + m_isoB * pt_gamma);
  const double thrNonIso = thrIso + m_isoGap;

  const bool iso    = (eiso_et < thrIso);
  const bool nonIso = (eiso_et > thrNonIso);

  // Exclude the gap region from A–B–C–D (PPG12 definition)
  if (!iso && !nonIso)
  {
    if (Verbosity() >= 5)
    {
      std::ostringstream os;
      os << "    [SS+Iso] trig=\"" << trig << "\" slice=\"" << slice << "\""
         << " pT^γ=" << std::fixed << std::setprecision(2) << pt_gamma
         << " | Eiso=" << std::setprecision(3) << eiso_et
         << " in GAP: [" << std::setprecision(3) << thrIso
         << ", " << std::setprecision(3) << thrNonIso << "]"
         << " → excluded from A–B–C–D";
      LOG(5, CLR_YELLOW, os.str());
    }
    return;
  }

  // -------------------------------------------------------------------------
  // Tight sub-cuts (PPG12 Table 4) — applied AFTER preselection
  //
  // IMPORTANT PPG12 detail:
  //   Non-tight is defined as "fails at least TWO of the five tight requirements".
  //   Clusters failing exactly ONE of the five tight requirements are NOT "non-tight"
  //   and must be EXCLUDED from the ABCD purity regions.
  // -------------------------------------------------------------------------
  const double w_hi = tight_w_hi(v.pt_gamma);

  const bool pass_weta   = in_open_interval(v.weta_cogx,     TIGHT_W_LO,       w_hi);
  const bool pass_wphi   = in_open_interval(v.wphi_cogx,     TIGHT_W_LO,       w_hi);
  const bool pass_e11e33 = in_open_interval(v.e11_over_e33,  TIGHT_E11E33_MIN, TIGHT_E11E33_MAX);
  const bool pass_et1    = in_open_interval(v.et1,           TIGHT_ET1_MIN,    TIGHT_ET1_MAX);
  const bool pass_e32e35 = in_open_interval(v.e32_over_e35,  TIGHT_E32E35_MIN, TIGHT_E32E35_MAX);

  const int tight_fails =
      (!pass_weta) + (!pass_wphi) + (!pass_e11e33) + (!pass_et1) + (!pass_e32e35);

  TightTag tag;
  if (tight_fails == 0)      tag = TightTag::kTight;
  else if (tight_fails >= 2) tag = TightTag::kNonTight;   // PPG12 "non-tight"
  else                       tag = TightTag::kNeither;    // exactly 1 fail (EXCLUDE from ABCD)

  // -------------------------------------------------------------------------
  // PPG12-equivalent ABCD gating:
  //   - Only Tight and NonTight(>=2 fails) enter the A–B–C–D regions
  //   - kNeither (exactly 1 fail) is explicitly excluded
  // -------------------------------------------------------------------------
  if (tag == TightTag::kNeither)
  {
    if (Verbosity() >= 5)
    {
      std::ostringstream os;
      os << "    [SS+Iso] trig=\"" << trig << "\" slice=\"" << slice << "\""
         << " pT^γ=" << std::fixed << std::setprecision(2) << pt_gamma
         << " | Eiso=" << std::setprecision(3) << eiso_et
         << " thrIso=" << std::setprecision(3) << thrIso
         << " thrNonIso=" << std::setprecision(3) << thrNonIso
         << " → region=" << (iso ? "ISO" : "NONISO")
         << " | tight=" << tightTagName(tag) << " (fails=" << tight_fails << ")"
         << " | NOTE: exactly-1-fail is NOT PPG12 non-tight → excluded from ABCD";
      LOG(5, CLR_YELLOW, os.str());
    }
    return;
  }

  // If this ever happens, it indicates the caller violated the expected control flow
  // (this function assumes preselection already passed).
  if (tag == TightTag::kPreselectionFail)
  {
    if (Verbosity() >= 4)
      LOG(4, CLR_YELLOW,
          "  [fillIsoSSTagCounters] tag==kPreselectionFail (unexpected here) → skip ABCD fill");
    return;
  }

  // -------------------------------------------------------------------------
  // Update per-slice counters (only for candidates that ACTUALLY enter ABCD)
  // -------------------------------------------------------------------------
  auto& S = m_catByTrig[trig][slice];
  S.seen += 1;

  // -------------------------------------------------------------------------
  // Choose A–B–C–D category using ISO vs NONISO sideband AND Tight vs NonTight(>=2 fails)
  //
  // Keep histogram names for backward compatibility:
  //   - "notTight" in names below now means "PPG12 non-tight (>=2 fails)" ONLY.
  // -------------------------------------------------------------------------
  const char* comboBase  = nullptr;
  const char* comboKeySS = nullptr;
  char        region     = '?';

  if (iso && tag == TightTag::kTight)
  {
    // Region A: tight, isolated
    region    = 'A';
    comboBase = "h_isIsolated_isTight";
    comboKeySS= "isIsolated_isTight";
    ++S.n_iso_tight;
  }
  else if (nonIso && tag == TightTag::kTight)
  {
    // Region B: tight, non-isolated (strict NONISO sideband)
    region    = 'B';
    comboBase = "h_notIsolated_isTight";
    comboKeySS= "notIsolated_isTight";
    ++S.n_nonIso_tight;
  }
  else if (iso && tag == TightTag::kNonTight)
  {
    // Region C: NON-TIGHT (>=2 fails), isolated
    region    = 'C';
    comboBase = "h_isIsolated_notTight";     // legacy name retained
    comboKeySS= "isIsolated_notTight";       // legacy key retained
    ++S.n_iso_nonTight;
  }
  else if (nonIso && tag == TightTag::kNonTight)
  {
    // Region D: NON-TIGHT (>=2 fails), non-isolated (strict NONISO sideband)
    region    = 'D';
    comboBase = "h_notIsolated_notTight";    // legacy name retained
    comboKeySS= "notIsolated_notTight";      // legacy key retained
    ++S.n_nonIso_nonTight;
  }
  else
  {
    // Should be unreachable because:
    //   - gap removed earlier
    //   - kNeither removed earlier
    //   - tag must now be either kTight or kNonTight
    if (Verbosity() >= 4)
    {
      std::ostringstream os;
      os << "  [fillIsoSSTagCounters] UNREACHABLE mapping hit:"
         << " iso=" << iso << " nonIso=" << nonIso
         << " tag=" << tightTagName(tag) << " (fails=" << tight_fails << ")"
         << " | Eiso=" << eiso_et
         << " thrIso=" << thrIso
         << " thrNonIso=" << thrNonIso
         << " → skip fills to be safe";
      LOG(4, CLR_YELLOW, os.str());
    }
    return;
  }

  // -------------------------------------------------------------------------
  // Count histogram for this ABCD region (per pT/cent slice)
  // -------------------------------------------------------------------------
  if (auto* hc = getOrBookCountHist(trig, comboBase, ptIdx, effCentIdx))
  {
    hc->Fill(1);
    bumpHistFill(trig, std::string(comboBase) + slice);
  }
  else
  {
    LOG(2, CLR_YELLOW,
        "  [fillIsoSSTagCounters] getOrBookCountHist returned nullptr for \""
        << comboBase << "\" slice=\"" << slice << '"');
  }

  // -------------------------------------------------------------------------
  // Verbose diagnostics (optional)
  // -------------------------------------------------------------------------
  if (Verbosity() >= 5)
  {
    // Optional: break isolation into PhotonClusterBuilder components (best-effort)
    double iso_em = std::numeric_limits<double>::quiet_NaN();
    double iso_hi = std::numeric_limits<double>::quiet_NaN();
    double iso_ho = std::numeric_limits<double>::quiet_NaN();

    const int cone10 = static_cast<int>(std::lround(10.0 * m_isoConeR));
    const char* k_em = nullptr;
    const char* k_hi = nullptr;
    const char* k_ho = nullptr;

    if (cone10 == 3)      { k_em = "iso_03_emcal"; k_hi = "iso_03_hcalin"; k_ho = "iso_03_hcalout"; }
    else if (cone10 == 4) { k_em = "iso_04_emcal"; k_hi = "iso_04_hcalin"; k_ho = "iso_04_hcalout"; }

    if (k_em && k_hi && k_ho)
    {
      if (auto* pho = dynamic_cast<const PhotonClusterv1*>(clus))
      {
        iso_em = pho->get_shower_shape_parameter(k_em);
        iso_hi = pho->get_shower_shape_parameter(k_hi);
        iso_ho = pho->get_shower_shape_parameter(k_ho);
      }
    }

    // (1) Compact summary line
    {
      std::ostringstream os;
      os << "    [SS+Iso] trig=\"" << trig << "\" slice=\"" << slice << "\""
         << " pT^γ=" << std::fixed << std::setprecision(2) << pt_gamma
         << " | Eiso=" << std::setprecision(3) << eiso_et
         << " thrIso=" << std::setprecision(3) << thrIso
         << " thrNonIso=" << std::setprecision(3) << thrNonIso
         << " → region=" << region << "(" << (iso ? "ISO" : "NONISO") << ")"
         << " | tight=" << tightTagName(tag) << " (fails=" << tight_fails << ")"
         << " | w_hi=" << std::setprecision(3) << w_hi;

      // Highlight the true signal region A in red; others cyan
      const char* colour = (region == 'A') ? CLR_RED : CLR_CYAN;
      LOG(5, colour, os.str());
    }

    // (2) Structured details
    if (Verbosity() >= 6)
    {
      // Isolation breakdown
      {
        std::ostringstream os;
        os << "      [iso parts] coneR=" << std::fixed << std::setprecision(2) << m_isoConeR
           << " : ";
        if (std::isfinite(iso_em) && std::isfinite(iso_hi) && std::isfinite(iso_ho))
        {
          os << k_em << "=" << std::setprecision(3) << iso_em
             << " " << k_hi << "=" << std::setprecision(3) << iso_hi
             << " " << k_ho << "=" << std::setprecision(3) << iso_ho
             << " → total=" << std::setprecision(3) << (iso_em + iso_hi + iso_ho);
        }
        else
        {
          os << "(components unavailable) total(Eiso)=" << std::setprecision(3) << eiso_et;
        }
        os << " | thrIso=" << std::setprecision(3) << thrIso
           << " | thrNonIso=" << std::setprecision(3) << thrNonIso;
        LOG(6, CLR_BLUE, os.str());
      }

      // Shower-shape values + cut windows (tight)
      {
        std::ostringstream os;
        os << "      [SS vars] "
           << "weta=" << std::fixed << std::setprecision(3) << v.weta_cogx
           << " ∈ (" << TIGHT_W_LO << "," << std::setprecision(3) << w_hi << ") → " << (pass_weta ? "PASS" : "FAIL")
           << " | wphi=" << std::setprecision(3) << v.wphi_cogx
           << " ∈ (" << TIGHT_W_LO << "," << std::setprecision(3) << w_hi << ") → " << (pass_wphi ? "PASS" : "FAIL")
           << " | et1=" << std::setprecision(3) << v.et1
           << " ∈ (" << TIGHT_ET1_MIN << "," << TIGHT_ET1_MAX << ") → " << (pass_et1 ? "PASS" : "FAIL")
           << " | e11/e33=" << std::setprecision(3) << v.e11_over_e33
           << " ∈ (" << TIGHT_E11E33_MIN << "," << TIGHT_E11E33_MAX << ") → " << (pass_e11e33 ? "PASS" : "FAIL")
           << " | e32/e35=" << std::setprecision(3) << v.e32_over_e35
           << " ∈ (" << TIGHT_E32E35_MIN << "," << TIGHT_E32E35_MAX << ") → " << (pass_e32e35 ? "PASS" : "FAIL");
        LOG(6, CLR_BLUE, os.str());
      }
    }
  }

  // -------------------------------------------------------------------------
  // SS variable histograms for THIS ABCD category (only for ISO or NONISO sideband)
  // -------------------------------------------------------------------------
  auto fillSS = [&](const std::string& key, double val)
  {
    if (auto* h = getOrBookSSHist(trig, key, comboKeySS, ptIdx, effCentIdx))
    {
      h->Fill(val);
      bumpHistFill(trig, "h_ss_" + key + "_" + comboKeySS + slice);
    }
    else
    {
      LOG(2, CLR_YELLOW,
          "    [fillSS] getOrBookSSHist returned nullptr for key=\"" << key
          << "\" tagKey=\"" << comboKeySS << "\" slice=\"" << slice << '"');
    }
  };

  fillSS("weta",   v.weta_cogx);
  fillSS("wphi",   v.wphi_cogx);
  fillSS("et1",    v.et1);
  fillSS("e11e33", v.e11_over_e33);
  fillSS("e32e35", v.e32_over_e35);

  // -------------------------------------------------------------------------
  // No additional counting histograms beyond the four ABCD regions:
  //   A: h_isIsolated_isTight
  //   B: h_notIsolated_isTight
  //   C: h_isIsolated_notTight
  //   D: h_notIsolated_notTight
  // -------------------------------------------------------------------------
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
  cout << left << setw(31) << "early reject: pT^#gamma<5"  << " | " << right << setw(7) << m_bk.pho_early_E         << '\n';
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
