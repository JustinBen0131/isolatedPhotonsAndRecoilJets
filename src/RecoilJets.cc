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
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>
#include <calobase/TowerInfoDefs.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/RawCluster.h>
#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/coresoftware_local/offline/packages/CaloBase/PhotonClusterv1.h"
#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/coresoftware_local/offline/packages/CaloReco/PhotonClusterBuilder.h"
#include <clusteriso/ClusterIso.h>
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
  double truth_x = std::numeric_limits<double>::quiet_NaN();
  double truth_y = std::numeric_limits<double>::quiet_NaN();
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
        truth_x = vtx->get_x();
        truth_y = vtx->get_y();
        truth_z = vtx->get_z();
        haveTruthZ = std::isfinite(truth_z);
      }
      else
      {
        LOG(3, CLR_YELLOW, "    [fetchNodes] isSim: G4TruthInfo has no primary vertex");
      }
    }
  }

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

  // Choose the z we will USE:
  //   • isSim : TRUTH (if present) → else MBD → else Global
  //   • data  : MBD → else Global
  const char* vz_source = "none";
  if (isSim && haveTruthZ)
  {
    m_vx = static_cast<float>(truth_x);
    m_vy = static_cast<float>(truth_y);
    m_vz = static_cast<float>(truth_z);
    vz_source = "TRUTH";
  }
  else if (haveMBDZ)
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
    LOG(1, CLR_YELLOW, "  [fetchNodes] no usable vertex z → skip event");
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

  if (!m_photons)
  {
    LOG(0, CLR_YELLOW,
        "    [fetchNodes] PHOTONCLUSTER_CEMC is MISSING → ABORTEVENT. "
        "PhotonClusterBuilder likely did not run or node name mismatch.");
    return false;
  }

    // ------------------------------------------------------------------
    // Reco jets: cache ALL radii listed in kJetRadii (r02 + r04) in parallel.
    //
    // Legacy knob (still honored as a "primary" label in logs only):
    //   export RJ_RECO_JET_KEY=r02   (or r04)
    //
    // IMPORTANT: we DO NOT disable the other radius — all jet QA / matching
    //            runs for every entry in kJetRadii.
    // ------------------------------------------------------------------
    std::string primaryRecoKey = trim(m_xjRecoJetKey);
    if (const char* rk = std::getenv("RJ_RECO_JET_KEY"))
      primaryRecoKey = trim(std::string(rk));

    if (primaryRecoKey.empty())
      primaryRecoKey = kJetRadii.front().key;

    bool keyOK = false;
    for (const auto& jnm : kJetRadii)
      if (primaryRecoKey == jnm.key) { keyOK = true; break; }

    if (!keyOK)
    {
      primaryRecoKey = kJetRadii.front().key;
      LOG(1, CLR_YELLOW,
          "    [fetchNodes] requested RJ_RECO_JET_KEY not in kJetRadii → using \"" << primaryRecoKey << "\"");
    }

    // Persist (legacy) — may be printed elsewhere
    m_xjRecoJetKey = primaryRecoKey;

    m_jets.clear();
    for (const auto& jnm : kJetRadii)
    {
      const std::string rKey = jnm.key;

      const std::string node = (isAuAu ? jnm.aa_node : jnm.pp_node);
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
    }

    if (Verbosity() >= 4)
    {
      std::ostringstream os;
      os << "    [fetchNodes] reco jet radii cached: {";
      for (std::size_t i = 0; i < kJetRadii.size(); ++i)
      {
        if (i) os << ", ";
        const std::string rk = kJetRadii[i].key;
        auto it = m_jets.find(rk);
        const bool ok = (it != m_jets.end() && it->second);
        os << rk << ":" << (ok ? "OK" : "MISSING");
      }
      os << "}  (primary=" << primaryRecoKey << ")";
      LOG(4, CLR_BLUE, os.str());
    }

    // ------------------------------------------------------------------
    // Truth jets: SIM only.
    // Cache truth jet containers for ALL radii in kJetRadii (r02 + r04).
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

      for (const auto& jnm : kJetRadii)
      {
        const std::string rKey = jnm.key;

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

        // job-wide cutflow
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

  delete out; out = nullptr;          
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
    LOG(2, CLR_YELLOW, "  [processCandidates] PHOTONCLUSTER_CEMC and CLUSTERINFO_CEMC both MISSING – nothing to process");
    return;
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

              if (Verbosity() >= 6)
                LOG(6, CLR_GREEN, "      [pho#" << iPho << "] marked as event-leading iso∧tight photon for jet matching (pT="
                                                << std::fixed << std::setprecision(2) << leadPtGamma << ")");
            }

            // Defer jet matching (prevents double-filling JES TH3s if multiple photons pass).
            continue;
        } // photon loop

        // ------------------------------------------------------------------
        // Jet matching + JES fills ONCE per event using the event-leading
        // iso∧tight photon (by pT^gamma), but run identically for ALL radii
        // in kJetRadii (r02 + r04) in the SAME loop.
        // ------------------------------------------------------------------
        if (haveLeadIsoTight)
        {
          // Centrality index used for matching-QA (centrality-only suffix in Au+Au)
          const int effCentIdx_M = (m_isAuAu ? centIdx : -1);

          // ΔR helper (used for SIM truth matching)
          auto dR = [](double eta1, double phi1, double eta2, double phi2) -> double
          {
            const double dphi = TVector2::Phi_mpi_pi(phi1 - phi2);
            const double deta = eta1 - eta2;
            return std::sqrt(deta*deta + dphi*dphi);
          };

          // -------------------- SIM: match truth photon ONCE per event --------------------
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
            else
            {
                // 1) TRUTH PHOTON definition for matching (MC):
                //    |eta| < 0.7, PID=22, prompt (direct OR fragmentation) from HepMC history,
                //    and truth isolation ETiso < 4 GeV within ΔR = 0.4 (final-state particles,
                //    excluding neutrinos and the photon itself). Then match to reco with ΔR < 0.05.

                constexpr double kTruthEtaAbsMax = 0.7;
                constexpr double kTruthIsoConeR  = 0.4;
                constexpr double kTruthIsoMaxGeV = 4.0;
                constexpr double kPhoMatchDR     = 0.05;

                double bestDR = 1e9;

                // Grab HepMC event (required for prompt/direct/fragmentation classification)
                PHHepMCGenEventMap* hepmcmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
                PHHepMCGenEvent*    hepmc    = nullptr;
                HepMC::GenEvent*    evt      = nullptr;

                if (hepmcmap)
                {
                  // Try common keys then fallback to first entry
                  hepmc = hepmcmap->get(0);
                  if (!hepmc) hepmc = hepmcmap->get(1);
                  if (!hepmc && !hepmcmap->empty()) hepmc = hepmcmap->begin()->second;
                  if (hepmc) evt = hepmc->getEvent();
                }

                if (!evt)
                {
                  LOG(4, CLR_YELLOW,
                      "      [truthQA] PHHepMCGenEventMap/HepMC event missing → cannot apply prompt+truth-iso definition → skip truth photon matching");
                }
                else
                {
                  // Identify beam particles (so we don't reject prompt photons just because beams are hadrons)
                  const HepMC::GenParticle* beam1 = nullptr;
                  const HepMC::GenParticle* beam2 = nullptr;
                  if (evt->valid_beam_particles())
                  {
                    beam1 = evt->beam_particles().first;
                    beam2 = evt->beam_particles().second;
                  }

                  auto isBeamParticle = [&](const HepMC::GenParticle* gp) -> bool
                  {
                    if (!gp) return false;
                    if (gp == beam1 || gp == beam2) return true;

                    // If beams aren't flagged, treat incoming hadrons with no production vertex as "beam-like"
                    const int apdg = std::abs(gp->pdg_id());
                    if (apdg > 100 && gp->production_vertex() == nullptr) return true;

                    return false;
                  };

                  auto isNeutrino = [](int pdg) -> bool
                  {
                    const int apdg = std::abs(pdg);
                    return (apdg == 12 || apdg == 14 || apdg == 16);
                  };

                  auto deltaR = [&](double eta1, double phi1, double eta2, double phi2) -> double
                  {
                    const double dphi = TVector2::Phi_mpi_pi(phi1 - phi2);
                    const double deta = eta1 - eta2;
                    return std::sqrt(deta*deta + dphi*dphi);
                  };

                  // Prompt (direct or fragmentation) photon classifier:
                  // Reject only if there is a NON-beam hadron ancestor (pi0/eta/...).
                  auto isPromptPhoton = [&](const HepMC::GenParticle* pho) -> bool
                  {
                    if (!pho) return false;

                    std::vector<const HepMC::GenParticle*> stack;
                    stack.reserve(64);

                    const HepMC::GenVertex* vtx = pho->production_vertex();
                    if (!vtx) return true;  // best-effort: no ancestry info

                    for (auto pit = vtx->particles_in_const_begin(); pit != vtx->particles_in_const_end(); ++pit)
                      if (*pit) stack.push_back(*pit);

                    int steps = 0;
                    while (!stack.empty() && steps < 5000)
                    {
                      const HepMC::GenParticle* a = stack.back();
                      stack.pop_back();
                      ++steps;

                      if (!a) continue;

                      const int apdg = std::abs(a->pdg_id());

                      // Hadron ancestor (not beam) → NOT prompt (decay photon)
                      if (apdg > 100 && !isBeamParticle(a))
                        return false;

                      const HepMC::GenVertex* pv = a->production_vertex();
                      if (!pv) continue;

                      for (auto ip = pv->particles_in_const_begin(); ip != pv->particles_in_const_end(); ++ip)
                        if (*ip) stack.push_back(*ip);
                    }

                    return true;
                  };

                  // Loop over HepMC particles to find truth photons satisfying the NEW definition
                  for (auto it = evt->particles_begin(); it != evt->particles_end(); ++it)
                  {
                    const HepMC::GenParticle* p = *it;
                    if (!p) continue;
                    if (p->pdg_id() != 22) continue;

                    // Use final-state photons
                    if (p->status() != 1) continue;

                    const double pt = std::hypot(p->momentum().px(), p->momentum().py());
                    if (!std::isfinite(pt) || pt <= 0.0) continue;

                    const double eta = p->momentum().pseudoRapidity();
                    const double phi = p->momentum().phi();
                    if (!std::isfinite(eta) || !std::isfinite(phi)) continue;

                    if (std::fabs(eta) >= kTruthEtaAbsMax) continue;

                    // Prompt/direct/fragmentation requirement
                    if (!isPromptPhoton(p)) continue;

                    // Truth isolation: sum pT of FINAL-STATE particles in ΔR<0.4 (exclude neutrinos + this photon)
                    double iso = 0.0;
                    for (auto jt = evt->particles_begin(); jt != evt->particles_end(); ++jt)
                    {
                      const HepMC::GenParticle* q = *jt;
                      if (!q) continue;
                      if (q == p) continue;
                      if (q->status() != 1) continue;
                      if (isNeutrino(q->pdg_id())) continue;

                      const double qpt = std::hypot(q->momentum().px(), q->momentum().py());
                      if (!std::isfinite(qpt) || qpt <= 0.0) continue;

                      const double qeta = q->momentum().pseudoRapidity();
                      const double qphi = q->momentum().phi();
                      if (!std::isfinite(qeta) || !std::isfinite(qphi)) continue;

                      if (deltaR(eta, phi, qeta, qphi) < kTruthIsoConeR) iso += qpt;
                    }

                    if (!std::isfinite(iso) || iso >= kTruthIsoMaxGeV) continue;

                    // Match reco ↔ truth photon (ΔR < 0.05)
                    const double drPho = deltaR(leadEtaGamma, leadPhiGamma, eta, phi);
                    if (drPho < kPhoMatchDR && drPho < bestDR)
                    {
                      bestDR = drPho;
                      haveTruthPho = true;
                      tPt  = pt;
                      tEta = eta;
                      tPhi = phi;
                    }
                  }

                  if (!haveTruthPho && Verbosity() >= 5)
                  {
                    LOG(5, CLR_YELLOW,
                        "      [truthQA] no truth prompt photon matched (ΔR<0.05) with truth iso ETiso<4 GeV in ΔR=0.4");
                  }
                }
            }
          }
          // ---------------------------------------------------------------------------

          if (Verbosity() >= 5)
          {
            LOG(5, CLR_CYAN,
                "      [lead pho#" << leadPhoIndex << "] jet matching for all radii"
                << " | pT^γ=" << std::fixed << std::setprecision(2) << leadPtGamma
                << " eta^γ=" << std::fixed << std::setprecision(3) << leadEtaGamma
                << " phi^γ=" << std::fixed << std::setprecision(3) << leadPhiGamma);
          }

          bool filledAnyRadius = false;

          // Run the EXACT SAME jet logic for every radius in kJetRadii (r02 + r04)
          for (const auto& jnm : kJetRadii)
          {
            const std::string rKey = jnm.key;

            JetContainer* jets = nullptr;
            if (auto itJ = m_jets.find(rKey); itJ != m_jets.end()) jets = itJ->second;

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
            // IMPORTANT: jet2 for alpha is the SUBLEADING recoil jet (same recoil selection),
            // not just the 2nd jet in the event.
            int nPassPt   = 0;
            int nPassEta  = 0;
            int nPassDphi = 0;

            // Leading & subleading recoil jets (pass pT+eta+dphi)
            double recoil1Pt = -1.0;
            const Jet* recoil1Jet = nullptr;
            double recoil2Pt = -1.0;
            const Jet* recoil2Jet = nullptr;

            // Leading & subleading in full (pT+eta) set (for isLeading QA)
            double all1Pt = -1.0;
            const Jet* all1Jet = nullptr;
            double all2Pt = -1.0;
            const Jet* all2Jet = nullptr;

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

              // Track top-2 jets in the pT+eta set (for recoil-is-leading QA only)
              if (jpt > all1Pt)
              {
                all2Pt = all1Pt; all2Jet = all1Jet;
                all1Pt = jpt;    all1Jet = j;
              }
              else if (jpt > all2Pt)
              {
                all2Pt = jpt; all2Jet = j;
              }

              const double dphiAbs = std::fabs(TVector2::Phi_mpi_pi(jphi - leadPhiGamma));
              if (std::isfinite(dphiAbs) && dphiAbs > maxDphi) maxDphi = dphiAbs;

              if (dphiAbs >= m_minBackToBack)
              {
                ++nPassDphi;

                // Track leading/subleading recoil jets
                if (jpt > recoil1Pt)
                {
                  recoil2Pt  = recoil1Pt;
                  recoil2Jet = recoil1Jet;

                  recoil1Pt  = jpt;
                  recoil1Jet = j;
                }
                else if (jpt > recoil2Pt)
                {
                  recoil2Pt  = jpt;
                  recoil2Jet = j;
                }
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

            // Physics-output fills (xJ, alpha, JES3) + Jet13/Profile3D — radius-tagged
            if (recoil1Pt > 0.0 && recoil1Jet)
            {
              filledAnyRadius = true;

              const double xJ = recoil1Pt / leadPtGamma;

              // α uses SUBLEADING recoil jet; if none, α=0 (always fill TH3)
              const double jet2Pt = (recoil2Pt > 0.0 ? recoil2Pt : 0.0);
              const double alpha  = (leadPtGamma > 0.0 ? (jet2Pt / leadPtGamma) : 0.0);

              if (Verbosity() >= 5)
              {
                LOG(5, CLR_GREEN,
                    "      [lead pho#" << leadPhoIndex << "] rKey=" << rKey
                    << " (R=" << std::fixed << std::setprecision(2) << Rjet
                    << ", |eta|<" << std::fixed << std::setprecision(2) << jetEtaAbsMaxUse << ")"
                    << " | xJ=" << std::fixed << std::setprecision(3) << xJ
                    << " jet1Pt=" << recoil1Pt
                    << " jet2Pt(recoil)=" << jet2Pt
                    << " alpha=" << alpha
                    << " (passed dphi=" << nPassDphi
                    << ", passed pt=" << nPassPt << ")");
              }

              const int effCentIdx = (m_isAuAu ? centIdx : -1);

              // Jet-only QA for selected jets (jet2 is subleading recoil)
              fillSelectedJetQA(activeTrig, leadPtIdx, effCentIdx, rKey, recoil1Jet, recoil2Jet);

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
                JetContainer* truthJets = nullptr;
                if (auto itT = m_truthJetsByRKey.find(rKey); itT != m_truthJetsByRKey.end()) truthJets = itT->second;

                if (!truthJets)
                {
                  if (Verbosity() >= 4)
                    LOG(4, CLR_YELLOW, "      [truthQA] rKey=" << rKey << " truth jets missing → skip truth JES3 fills");
                }
                else if (tPt <= 0.0)
                {
                  if (Verbosity() >= 4)
                    LOG(4, CLR_YELLOW, "      [truthQA] rKey=" << rKey << " matched truth gamma has invalid pT → skip truth JES3 fills");
                }
                else
                {
                  // 2) truth recoil jets (same R), |eta| < 1.1 - R, away-side to truth gamma
                  const double etaMaxTruth = jetEtaAbsMaxForRKey(rKey);

                  double tj1Pt = -1.0; const Jet* tj1 = nullptr;
                  double tj2Pt = -1.0; const Jet* tj2 = nullptr;

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
                      tj2Pt = tj1Pt; tj2 = tj1;
                      tj1Pt = ptj;   tj1 = tj;
                    }
                    else if (ptj > tj2Pt)
                    {
                      tj2Pt = ptj; tj2 = tj;
                    }
                  }

                  if (!tj1 || tj1Pt <= 0.0)
                  {
                    if (Verbosity() >= 5)
                      LOG(5, CLR_YELLOW, "      [truthQA] rKey=" << rKey << " no truth recoil jet1 found → skip truth JES3 fills");
                  }
                  else
                  {
                    // 3) require reco jet1 ↔ truth jet1 match in ΔR
                    const double recoJeta = recoil1Jet->get_eta();
                    const double recoJphi = recoil1Jet->get_phi();
                    const double drJet = dR(recoJeta, recoJphi, tj1->get_eta(), tj1->get_phi());

                    if (drJet > 0.2)
                    {
                      if (Verbosity() >= 5)
                        LOG(5, CLR_YELLOW, "      [truthQA] rKey=" << rKey
                                        << " reco jet1 ↔ truth jet1 ΔR=" << drJet << " > 0.2 → skip truth JES3 fills");
                    }
                    else
                    {
                      const double xJt = tj1Pt / tPt;
                      const double aT  = (tj2Pt > 0.0 ? (tj2Pt / tPt) : 0.0);

                      for (const auto& trigShort : activeTrig)
                      {
                        if (auto* ht3x = getOrBookJES3Truth_xJ_alphaHist(trigShort, rKey, effCentIdx_M))
                        { ht3x->Fill(tPt, xJt, aT); bumpHistFill(trigShort, ht3x->GetName()); }

                        if (auto* ht3j = getOrBookJES3Truth_jet1Pt_alphaHist(trigShort, rKey, effCentIdx_M))
                        { ht3j->Fill(tPt, tj1Pt, aT); bumpHistFill(trigShort, ht3j->GetName()); }
                      }
                    }
                  }
                }
              }
              // -----------------------------------------------------------------------
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

          // Preserve your old "used" semantics: count once per event if any radius filled
          if (filledAnyRadius) ++nUsed;
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

  auto* h = new TH3F(name.c_str(), title.c_str(), nx, xbins, ny, ylo, yhi, nz, zlo, zhi);
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

  auto* h = new TH3F(name.c_str(), title.c_str(), nx, xbins, ny, ylo, yhi, nz, zlo, zhi);
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

  auto* p = new TProfile3D(name.c_str(), title.c_str(), nx, xbins, ny, ylo, yhi, nz, zlo, zhi);
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
  // Tight sub-cuts (PPG12 Table 4) — applied AFTER preselection (as you already do)
  // -------------------------------------------------------------------------
  const double w_hi = tight_w_hi(v.pt_gamma);

  const bool pass_weta   = in_open_interval(v.weta_cogx,     TIGHT_W_LO,       w_hi);
  const bool pass_wphi   = in_open_interval(v.wphi_cogx,     TIGHT_W_LO,       w_hi);
  const bool pass_e11e33 = in_open_interval(v.e11_over_e33,  TIGHT_E11E33_MIN, TIGHT_E11E33_MAX);
  const bool pass_et1    = in_open_interval(v.et1,          TIGHT_ET1_MIN,    TIGHT_ET1_MAX);
  const bool pass_e32e35 = in_open_interval(v.e32_over_e35, TIGHT_E32E35_MIN, TIGHT_E32E35_MAX);

  const int tight_fails  =
      (!pass_weta) + (!pass_wphi) + (!pass_e11e33) + (!pass_et1) + (!pass_e32e35);

  TightTag tag;
  if (tight_fails == 0)      tag = TightTag::kTight;
  else if (tight_fails >= 2) tag = TightTag::kNonTight;
  else                       tag = TightTag::kNeither;

  // -------------------------------------------------------------------------
  // Update per-slice counters (purity regions only)
  // -------------------------------------------------------------------------
  auto& S = m_catByTrig[trig][slice];
  S.seen += 1;

  // -------------------------------------------------------------------------
  // Choose A–B–C–D category using ISO vs NONISO sideband (NOT "fail iso")
  // Keep histogram names for backward compatibility.
  // -------------------------------------------------------------------------
  const char* comboBase  = nullptr;
  const char* comboKeySS = nullptr;

  if (iso && tag == TightTag::kTight)
  {
    // Region A
    comboBase  = "h_isIsolated_isTight";
    comboKeySS = "isIsolated_isTight";
    ++S.n_iso_tight;
  }
  else if (nonIso && tag == TightTag::kTight)
  {
    // Region B (strict NONISO sideband)
    comboBase  = "h_notIsolated_isTight";
    comboKeySS = "notIsolated_isTight";
    ++S.n_nonIso_tight;
  }
  else if (iso && tag != TightTag::kTight)
  {
    // Region C (includes kNonTight and kNeither as "not tight" in your 2×2)
    comboBase  = "h_isIsolated_notTight";
    comboKeySS = "isIsolated_notTight";
    ++S.n_iso_nonTight;
  }
  else
  {
    // Region D (strict NONISO sideband)
    comboBase  = "h_notIsolated_notTight";
    comboKeySS = "notIsolated_notTight";
    ++S.n_nonIso_nonTight;
  }

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
         << " → region=" << (iso ? "ISO" : "NONISO")
         << " | tight=" << tightTagName(tag) << " (fails=" << tight_fails << ")"
         << " | w_hi=" << std::setprecision(3) << w_hi;
      const char* colour = (iso && tag == TightTag::kTight) ? CLR_RED : CLR_CYAN;
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
  // SS variable histograms for THIS category (only for ISO or NONISO sideband)
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
  // Optional: ID sideband bookkeeping
  // Keep your existing behavior: count NonTight total, and NonTight+ISO pass
  // (this corresponds to C region for the ID sideband).
  // -------------------------------------------------------------------------
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
