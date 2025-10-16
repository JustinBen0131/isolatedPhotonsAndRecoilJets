// Tell Emacs this is C++   -*- C++ -*-
#ifndef RECOILJETS_H
#define RECOILJETS_H
//==========================================================================
//  EMCal × sEPD × MBD correlator – headers
//==========================================================================
//––– Framework ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
#include <fun4all/SubsysReco.h>
#include <calotrigger/TriggerAnalyzer.h>
#include <phool/PHCompositeNode.h>

//––– ROOT base ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
#include <TFile.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TH2F.h>
#include <TH2Poly.h>
#include <TLorentzVector.h>

//––– sPHENIX objects ––––––––––––––––––––––––––––––––––––––––––––––––––––––
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawClusterContainer.h>
#include <globalvertex/GlobalVertexMap.h>
#include <mbd/MbdGeom.h>
#include <mbd/MbdPmtContainer.h>
#include <epd/EpdGeom.h>
#include <epd/EpdReco.h>
#include <centrality/CentralityInfo.h>
#include <eventplaneinfo/EventplaneinfoMap.h>

//––– STL ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
#include <string>
#include <map>
#include <vector>
#include <tuple>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <array>
#include <unordered_map>
#include <bitset>          // ← for extractTriggerBits()
#include <TVector2.h>      // ← TVector2::Phi_mpi_pi in doCaloQA()
#include <jetbase/JetContainer.h>

// --------------------------------------------------------------------------
//  Helper alias: maps histogram name → ROOT object*
// --------------------------------------------------------------------------
using HistMap = std::map<std::string, TObject*>;
class CentralityInfo;
class MinimumBiasInfo;
class Fun4AllHistoManager;
class PHCompositeNode;
class RawCluster;

// ==========================================================================
//  CLASS DCaloTowerCalibECLARATION
// ==========================================================================
class RecoilJets : public SubsysReco
{
 public:
  static constexpr std::array<std::pair<const char*, const char*>, 2> kJetRadii {{
            {"r02", "AntiKt_TowerInfo_HIRecoSeedsSub_r02"},
            {"r04", "AntiKt_TowerInfo_HIRecoSeedsSub_r04"}
    }};
  // ---------- construction / destruction ----------------------------------
  explicit RecoilJets(const std::string& out = "caloTreeData.root");
  ~RecoilJets() override;

  // ---------- Fun4All hooks ------------------------------------------------
  int Init          (PHCompositeNode*) override;
  int InitRun       (PHCompositeNode*) override;
  int process_event (PHCompositeNode*) override;
  int ResetEvent    (PHCompositeNode*) override;
  int End           (PHCompositeNode*) override;
  int Reset         (PHCompositeNode*) override;
  void Print(const std::string& what = "ALL") const override;

  // ---------- user configuration ------------------------------------------
  void setVerbose  (int level)   { Verbosity(level);        }
  void setRunNumber(int r)       { m_runNumber = r;         }
  void setVzCut    (double c)    { m_vzCut = std::fabs(c);  }
  void enableVzCut (bool f = true) { m_useVzCut = f;        }
  void setCentralityEdges(const std::vector<int>& e) { m_centEdges = e; }

    // Photon-isolation configuration (PPG12-style baseline)
    //   coneR        : isolation cone ΔR
    //   isoFrac      : (ΣE_T in cone excluding tag) / E_T^γ  < isoFrac
    //   hadFrac      : (ΣE_T^HCAL near axis) / E_T^γ < hadFrac (leakage veto)
    //   dphiMinRad   : back-to-back requirement on |Δφ(γ,jet)|, default 7π/8
  void configurePhotonIsolation(double coneR = 0.4,
                                  double isoFrac = 0.10,
                                  double hadFrac = 0.10,
                                  double dphiMinRad = 7.0*M_PI/8.0)
  {
      m_isoConeR = std::max(0.05, coneR);
      m_isoFrac  = std::max(0.0, isoFrac);
      m_hadFrac  = std::max(0.0, hadFrac);
      m_minBackToBack = std::clamp(dphiMinRad, 0.0, M_PI);
  }

    // Configure E_T^γ bin edges (must be size>=2, monotonically increasing)
  void setGammaEtBins(const std::vector<double>& edges)
  {
      if (edges.size() >= 2) m_gammaEtBins = edges;
  }

 private:
  // ======================================================================
  // 1) One‑time booking helpers
  // ======================================================================
  void      createHistos_Data();                          // main booker
  TH2Poly*  makeMbdHitmap(const std::string&, const MbdGeom*, int arm);
  TH2F*     makeEpdHitmap(const std::string& name, EpdGeom* geom, int arm);
  bool      m_sepdMapReady {false};
  void      buildSepdChannelMap();
  bool m_isMinBias {false};
  // ======================================================================
  // 2) Per‑event helpers (called in process_event)
  // ======================================================================
  bool getCentralitySlice(int& lo, int& hi, std::string& tag) const;
  bool fetchNodes (PHCompositeNode*);                     // guards + cache
    
  void doCaloQA   (PHCompositeNode* topNode,
                   const std::vector<std::string>&);
  void doSepdQA   (const std::vector<std::string>&);
  void doMbdQA    (const std::vector<std::string>&);
  void fillCorrelations (const std::vector<std::string>&);

  // ======================================================================
  // 3) Configuration & run‑time state
  // ======================================================================
  // --- run‑wide -----------------------------------------------------------
  int         m_runNumber   = -1;
  bool        verbose       = true;
  double      m_vzCut       = 10.;        // [cm]
  double      m_towMinE   {0.050};
  bool        m_useVzCut    = true;
  const GlobalVertex* m_vtx {nullptr};
  double m_vx {0.}, m_vy {0.}, m_vz {0.};
  std::vector<unsigned> m_epdKey;
  std::vector<int> m_centEdges{0,10,20,30,40,60,80};
  int                         m_centBin   = -1;   // 0…99
  std::map<std::string,int>   m_centIdxCache;         // "0_10" → 0, etc.

  std::string       Outfile;              // ROOT output file name
  TFile*            out      = nullptr;
  TriggerAnalyzer*  trigAna  = nullptr;
  std::size_t       event_count = 0;

  // --- trigger bookkeeping if using TriggerAnalyzer package -----------------------------------------------
//  std::map<std::string, std::string> triggerNameMap {
//        {"MBD N&S >= 2", "MBD_NandS_geq_2"}};
    
  // --- trigger bookkeeping :  bit‑index  →  human‑readable key ----------
  std::map<int, std::string> triggerNameMap = {
      {10, "MBD_NS_geq_2"},
      {11, "MBD_NS_geq_1"},
      {12, "MBD_NS_geq_2_vtx_lt_10"},
      {13, "MBD_NS_geq_2_vtx_lt_30"},
      {14, "MBD_NS_geq_2_vtx_lt_150"},
      {15, "MBD_NS_geq_1_vtx_lt_10"},
      {16, "photon_6_plus_MBD_NS_geq_2_vtx_lt_10"},
      {17, "photon_8_plus_MBD_NS_geq_2_vtx_lt_10"},
      {18, "photon_10_plus_MBD_NS_geq_2_vtx_lt_10"},
      {19, "photon_12_plus_MBD_NS_geq_2_vtx_lt_10"},
      {20, "photon_6_plus_MBD_NS_geq_2_vtx_lt_150"},
      {21, "photon_8_plus_MBD_NS_geq_2_vtx_lt_150"},
      {22, "photon_10_plus_MBD_NS_geq_2_vtx_lt_150"},
      {23, "photon_12_plus_MBD_NS_geq_2_vtx_lt_150"}
  };
    
  std::map<std::string, HistMap>     qaHistogramsByTrigger;
  // ── trigger QA helpers ────────────────────────────────────────────────
  TH2I* h_MBTrigCorr   = nullptr;                 // 2‑D map: MinBias × Trigger
  std::unordered_map<std::string,int> m_trigBin;  // trigger‑key → x‑bin index

  // first‑event gate: MB + trigger selection (declared here, defined in .cc)
  bool firstEventCuts(PHCompositeNode*   topNode,
                        std::vector<std::string>& activeTrig);
  // --- analysis cuts ------------------------------------------------------
  const std::vector<float>               m_asymCuts {0.5f, 0.7f};
  const std::vector<float>               m_chi2Cuts {1.f, 4.f};
  const std::vector<float>               m_minClusE {2.f};
  const std::vector<std::pair<float,float>> m_ptBins {
        {2,4},{4,6},{6,8},{8,10},{10,12},{12,15},{15,20},{20,30} };

  // --- detector lists -----------------------------------------------------
  const std::vector<std::tuple<std::string,std::string,std::string>> m_caloInfo {
        {"TOWERINFO_CALIB_CEMC",   "TOWERGEOM_CEMC_DETAILED",   "CEMC"},
        {"TOWERINFO_CALIB_HCALIN", "TOWERGEOM_HCALIN", "IHCAL"},
        {"TOWERINFO_CALIB_HCALOUT","TOWERGEOM_HCALOUT","OHCAL"} };

  // --- run‑time caches ----------------------------------------------------
  struct CaloCache {
    TowerInfoContainer*    tw = nullptr;
    RawTowerGeomContainer* g  = nullptr;
    double                 sumE = 0.;
  };
  std::map<std::string, CaloCache> m_calo;         // keyed by "CEMC"/…

  TowerInfoContainer*  m_sepd     = nullptr;
  MbdPmtContainer*     m_mbdpmts  = nullptr;
  MbdGeom*             m_mbdgeom  = nullptr;
  EpdGeom*             m_epdgeom  = nullptr;
  EventplaneinfoMap*   m_epmap    = nullptr;   // pointer to EventplaneinfoMap
  RawClusterContainer* m_clus     = nullptr;


  // --------------------------------------------------------------------------

};
#endif  // RECOILJETS_H
