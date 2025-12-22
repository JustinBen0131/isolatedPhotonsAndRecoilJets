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
#include <TH1I.h>
#include <TH3F.h>
#include <TH2F.h>
#include <TH2Poly.h>
#include <TLorentzVector.h>

//––– sPHENIX objects ––––––––––––––––––––––––––––––––––––––––––––––––––––––
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>
#include <calobase/RawTowerGeomContainer.h>           // <-- required base
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawClusterContainer.h>
#include <globalvertex/GlobalVertexMap.h>
#include <mbd/MbdGeom.h>
#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/coresoftware_local/offline/packages/CaloReco/PhotonClusterBuilder.h"
#include <mbd/MbdPmtContainer.h>
#include <epd/EpdGeom.h>
#include <epd/EpdReco.h>
#include <centrality/CentralityInfo.h>
#include <eventplaneinfo/EventplaneinfoMap.h>
#include <jetbase/JetContainer.h>

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
#include <bitset>
#include <TVector2.h>     

// --------------------------------------------------------------------------
//  Helper alias: maps histogram name → ROOT object*
// --------------------------------------------------------------------------
using HistMap = std::map<std::string, TObject*>;
class CentralityInfo;
class MinimumBiasInfo;
class Fun4AllHistoManager;
class PHCompositeNode;
class RawCluster;
class PhotonClusterv1;
class PhotonClusterContainer;
class Jet;       // forward declaration for jet QA helpers
class TProfile;  // forward declaration (matching-QA profile)

// ---------------- PPG-12 working-point constants (Table 4) -----------------
namespace PhoIDCuts {
  // Pre-selection windows
  inline constexpr double PRE_E11E33_MAX = 0.98;
  inline constexpr double PRE_ET1_MIN    = 0.60;
  inline constexpr double PRE_ET1_MAX    = 1.00;
  inline constexpr double PRE_E32E35_MIN = 0.80;
  inline constexpr double PRE_E32E35_MAX = 1.00;
  inline constexpr double PRE_WETA_MAX   = 0.60;  // only η width is capped

  // Tight windows (two carry pT dependence; for photons pT == ET)
  inline constexpr double TIGHT_W_LO     = 0.0;   // 0 < w < 0.15 + 0.006 pT
  inline double tight_w_hi(double pt_gamma) { return 0.15 + 0.006 * pt_gamma; }

  inline constexpr double TIGHT_E11E33_MIN = 0.40;
  inline constexpr double TIGHT_E11E33_MAX = 0.98;
  inline constexpr double TIGHT_ET1_MIN    = 0.90;
  inline constexpr double TIGHT_ET1_MAX    = 1.00;
  inline constexpr double TIGHT_E32E35_MIN = 0.92;
  inline constexpr double TIGHT_E32E35_MAX = 1.00;

  // Inequality helpers
  inline bool in_open_interval (double x, double lo, double hi) { return (x>lo) && (x<hi); }
  inline bool in_closed_left   (double x, double lo, double hi) { return (x>=lo) && (x<hi); }
  inline bool in_closed_right  (double x, double lo, double hi) { return (x>lo) && (x<=hi); }
  inline bool in_semi_open     (double x, double lo, double hi) { return (x>lo) && (x<=hi); }
}

// ==========================================================================
//  CLASS DECLARATION
// ==========================================================================
class RecoilJets : public SubsysReco
{
 public:
  // Jet node names per radius for pp (raw) and Au+Au (UE-subtracted)
  struct JetNodeMap { const char* key; const char* pp_node; const char* aa_node; };
  static constexpr std::array<JetNodeMap, 1> kJetRadii {{
      {"r02", "AntiKt_TowerInfo_PPJetsRaw_r02", "AntiKt_TowerInfo_HIRecoSeedsSub_r02"}
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

  // ---------- user configuration ------------------------------------------
  void setVerbose  (int level)   { Verbosity(level);        }
  void setRunNumber(int r)       { m_runNumber = r;         }
  void setVzCut    (double c)    { m_vzCut = std::fabs(c);  }
  void enableVzCut (bool f = true) { m_useVzCut = f;        }
  void setCentralityEdges(const std::vector<int>& e) { m_centEdges = e; }

  // Global photon fiducial (PPG-12): |eta^gamma| < 0.7 by default
  void setEtaAbsMax(double a) { m_etaAbsMax = std::max(0.0, a); }

    // ---- data type (pp vs Au+Au) and sim flag ------------------------------
    enum class DataType { kPP, kAuAu };

    // string-friendly setter for Fun4All macros
    // Accepted tokens:
    //   isPP / pp
    //   isAuAu / auau
    //   isSim / sim   (simulation; pp-style logic but trigger-dependent code can be skipped)
    void setDataType(const std::string& s)
    {
          const std::string t = s;

          // Simulation token
          if (t=="isSim" || t=="issim" || t=="sim" || t=="SIM" || t=="Sim")
          {
            m_isSim = true;
            m_isAuAu = false;   // sim behaves like pp-type for Au+Au gating
            return;
          }

          // Data tokens
          m_isSim = false;
          if (t=="isAuAu" || t=="AuAu" || t=="AUAU" || t=="auau" || t=="AA")
            m_isAuAu = true;
          else
            m_isAuAu = false; // default: pp
    }

    // type-safe setter (optional)
    void setDataType(DataType t)
    {
          m_isSim = false;
          m_isAuAu = (t==DataType::kAuAu);
    }

    // queries
    bool isAuAu() const { return m_isAuAu; }
    bool isSim()  const { return m_isSim;  }


  // ---------- Shower‑shape selection (PPG‑12 style) ------------------------
  struct SSVars {
    double weta_cogx;       // w_η^cogX (seed-excluded 2nd moment in η)
    double wphi_cogx;       // w_φ^cogX (seed-excluded 2nd moment in φ)
    double et1;             // 2×2 core fraction proxy used in note
    double e11_over_e33;    // E11/E33 compactness
    double e32_over_e35;    // strip ratio (3×2 over 3×5, CoG-aligned)
    double pt_gamma;        // photon candidate pT (GeV) used everywhere for binning & width caps
  };

  enum class TightTag : uint8_t {
    kPreselectionFail = 0,  // fails pre-selection (do not consider as candidate)
    kTight            = 1,  // passes all 5 tight SS cuts
    kNonTight         = 2,  // passes pre-selection but fails ≥ 2 tight SS cuts
    kNeither          = 3   // passes pre-selection, fails exactly 1 tight cut
  };

  // API: simple, stable signatures you can call from your cluster loop
  bool     passesPhotonPreselection(const SSVars& v);
  TightTag classifyPhotonTightness (const SSVars& v);

  // ---------- Photon‑isolation configuration (PPG‑12) ----------------------
  //   coneR        : isolation cone ΔR
  //   isoFrac      : (ΣE_T in cone excluding tag) / E_T^γ  < isoFrac
  //   hadFrac      : (ΣE_T^HCAL near axis) / E_T^γ < hadFrac (leakage veto)
  //   dphiMinRad   : back-to-back requirement on |Δφ(γ,jet)|, default 7π/8
  void configurePhotonIsolation(double coneR = 0.4,
                                double isoFrac = 0.10,
                                double hadFrac = 0.10,
                                double dphiMinRad = 7.0*M_PI/8.0)
  {
      m_isoConeR      = std::max(0.05, coneR);
      m_isoFracLegacy = std::max(0.0, isoFrac);
      m_hadFrac       = std::max(0.0, hadFrac);
      m_minBackToBack = std::clamp(dphiMinRad, 0.0, M_PI);
  }

  // ---------- Isolation working point (signal line + sideband) -------------
  void setIsolationWP(double aGeV      = 1.08128,
                        double bPerGeV   = 0.0299107,
                        double sideGapGeV= 1.0,
                        double coneR     = 0.3,
                        double towerMin  = 0.060);

  // Jet selection baseline (applies to jet QA + away-side jet scans for xJ)
  // This is intentionally a reco-level technical floor.
  void setMinJetPt(double pt) { m_minJetPt = std::max(0.0, pt); }

  // Jet fiducial acceptance for the *jet axis*.
  // Defensible first-pass for R=0.2 barrel TowerInfo jets:
  //   |eta_jet| < (eta_max_calo - R) ≈ 1.1 - 0.2 = 0.9
  void setJetEtaAbsMax(double a) { m_jetEtaAbsMax = std::max(0.0, a); }

  // Compute Eiso for a cluster (EMCal+iHCal+oHCal in ΔR=coneR, excl. the cluster)
  double eiso(const RawCluster* clus, PHCompositeNode* topNode) const;

  // PPG-12: "isolated" (signal line)
  bool isIsolated(const RawCluster* clus, double et_gamma, PHCompositeNode* topNode) const;

  // PPG-12: "non-isolated sideband" (line + gap)
  bool isNonIsolated(const RawCluster* clus, double et_gamma, PHCompositeNode* topNode) const;

  // ---------- Binning controls (user-editable arrays) ----------------------
  // Configure photon pT bin edges (must be size>=2, monotonically increasing).
  void setGammaPtBins(const std::vector<double>& edges) { if (edges.size()>=2) m_gammaPtBins = edges; }

  // Backward-compatible alias (if any macros still call the old name)
  void setGammaEtBins(const std::vector<double>& edges) { setGammaPtBins(edges); }

    // Configure centrality edges (already available via setCentralityEdges()).

 private:
  // ======================================================================
  // 1) One‑time booking helpers
  // ======================================================================
  void      createHistos_Data();                          // main booker
  TH2Poly*  makeMbdHitmap(const std::string&, const MbdGeom*, int arm);
  TH2F*     makeEpdHitmap(const std::string& name, EpdGeom* geom, int arm);
  bool      m_sepdMapReady {false};
  void      buildSepdChannelMap();
  bool      m_isMinBias {false};

  // ======================================================================
  // 2) Per‑event helpers (called in process_event)
  // ======================================================================
  bool getCentralitySlice(int& lo, int& hi, std::string& tag) const;
  bool fetchNodes (PHCompositeNode*);                     // guards + cache
  void doCaloQA   (PHCompositeNode* topNode, const std::vector<std::string>&);
  void doSepdQA   (const std::vector<std::string>&);
  void doMbdQA    (const std::vector<std::string>&);
  void fillCorrelations (const std::vector<std::string>&);

  // ======================================================================
  // 3) Configuration & run‑time state
  // ======================================================================
  // --- run-wide -----------------------------------------------------------
  int         m_runNumber   = -1;
  bool        verbose       = true;
  double      m_vzCut       = 30.;        // [cm]
  double      m_towMinE     {0.050};
  bool        m_useVzCut    = true;

  // PPG-12 photon fiducial: |eta^gamma| < 0.7
  double      m_etaAbsMax   {0.7};

  // data-type flags:
  //   m_isAuAu: controls Au+Au-only logic (centrality path, UE-subtracted jets, etc.)
  //   m_isSim : marks simulation so trigger/GL1-dependent logic can be skipped
  bool        m_isAuAu      = false;      // default: p+p
  bool        m_isSim       = false;      // default: data

  // vertex & caches
  const GlobalVertex* m_vtx {nullptr};
  double m_vx {0.}, m_vy {0.}, m_vz {0.};
  std::vector<unsigned> m_epdKey;

  // User centrality edges and current bin (0…99)
  std::vector<int> m_centEdges{0,10,20,30,40,60,80};
  int              m_centBin   = -1;   // 0…99
  std::map<std::string,int> m_centIdxCache; // "0_10" → 0, etc.

  // Photon pT bin edges used everywhere for slicing & JES calibration
  std::vector<double> m_gammaPtBins {10, 15, 20, 50};

  // Output file & QA store
  std::string  Outfile;
  TFile*       out      = nullptr;
  TriggerAnalyzer*  trigAna  = nullptr;
  std::size_t   event_count = 0;
  std::size_t   m_evtNoTrig = 0;  // count events rejected by MB/trigger gate

  std::map<std::string, HistMap>     qaHistogramsByTrigger;
  TH2I* h_MBTrigCorr   = nullptr;                 // 2-D map: MinBias × Trigger
  std::unordered_map<std::string,int> m_trigBin;  // trigger-key → x-bin index

    struct CatStat {
        // photons that reached preselection (this slice)
        std::size_t seen{0};

        // LEGACY (kept for compatibility; not used in new table)
        std::size_t tight{0};
        std::size_t nonTight{0};
        std::size_t isoPass{0};
        std::size_t isoFail{0};

        // New 2×2 joint categories (exactly one increments per preselection-passed photon)
        std::size_t n_iso_tight{0};        // isolated ∧ tight
        std::size_t n_nonIso_tight{0};     // ¬isolated ∧ tight
        std::size_t n_iso_nonTight{0};     // isolated ∧ ¬tight (NonTight or Neither)
        std::size_t n_nonIso_nonTight{0};  // ¬isolated ∧ ¬tight

        // ID sideband bookkeeping (unchanged)
        std::size_t idSB_total{0};
        std::size_t idSB_pass{0};
    };

    
  // Per trigger -> per slice suffix (e.g., "_ET_10_12[_cent_0_10]") → counters
  std::map<std::string, std::map<std::string, CatStat>> m_catByTrig;
    // Histogram fill counts: "<trig>::<histNameWithSuffix>" → fills
    std::map<std::string, std::size_t> m_histFill;

    // ---------- Run-wide cutflow bookkeeping ----------
    enum class EventReject : uint8_t { None = 0, Trigger = 1, Vz = 2 };

    struct CutBookkeeping
    {
      // Events
      std::size_t evt_seen{0};
      std::size_t evt_fail_trigger{0};
      std::size_t evt_fail_vz{0};
      std::size_t evt_accepted{0};

      // Photon candidates from PHOTONCLUSTER_CEMC path
      std::size_t pho_total{0};
      std::size_t pho_early_E{0};
      std::size_t pho_eta_fail{0};
      std::size_t pho_etbin_out{0};
      std::size_t pho_noRC{0};
      std::size_t pho_reached_pre_iso{0};  // reached the preselection+isolation step

      // Preselection breakdown (fail reasons)
      std::size_t pre_fail_weta{0};
      std::size_t pre_fail_et1_low{0};
      std::size_t pre_fail_et1_high{0};
      std::size_t pre_fail_e11e33_high{0};
      std::size_t pre_fail_e32e35_low{0};
      std::size_t pre_fail_e32e35_high{0};
      std::size_t pre_pass{0};

      // Tight classification breakdown (among preselection-passed)
      std::size_t tight_tight{0};
      std::size_t tight_neither{0};    // exactly 1 tight cut failed
      std::size_t tight_nonTight{0};   // ≥2 tight cuts failed

      // Per-tight-cut failures (counted over all preselection-passed)
      std::size_t tight_fail_weta{0};
      std::size_t tight_fail_wphi{0};
      std::size_t tight_fail_et1{0};
      std::size_t tight_fail_e11e33{0};
      std::size_t tight_fail_e32e35{0};

      // Isolation
      std::size_t iso_pass{0};
      std::size_t iso_fail{0};
    };

    CutBookkeeping m_bk;
    EventReject    m_lastReject{EventReject::None};


    // Helper to record histogram fills
    void bumpHistFill(const std::string& trig, const std::string& hnameWithSuffix);

    // Human-readable, tabulated end-of-job summary
    void printCutSummary() const;

    // first-event gate: MB + trigger selection (declared here, defined in .cc)
    bool firstEventCuts(PHCompositeNode*   topNode,
                        std::vector<std::string>& activeTrig);

  // --- detector lists -----------------------------------------------------
  const std::vector<std::tuple<std::string,std::string,std::string>> m_caloInfo {
        {"TOWERINFO_CALIB_CEMC",   "TOWERGEOM_CEMC",   "CEMC"},
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
  EventplaneinfoMap*   m_epmap    = nullptr;
  RawClusterContainer* m_clus     = nullptr;
  PhotonClusterContainer* m_photons = nullptr;

  // Jet containers keyed by radius tag ("r02", "r04", …); we’ll use r02
  std::map<std::string, JetContainer*> m_jets;
    // Jet selection baseline pT floor (GeV) used for:
    //   - pure jet QA (to suppress pathological ultra-low-pT jets)
    //   - recoil jet scan for xJ
    double m_minJetPt {5.0};

    // Jet fiducial acceptance (jet axis).
    // Default is a containment-motivated first pass for R=0.2:
    //   |eta_jet| < 0.9
    double m_jetEtaAbsMax {0.9};

  // --------------------------------------------------------------------------
  // Triggers (as in your code) – pp & Au+Au name maps
  // --------------------------------------------------------------------------
  std::map<std::string,std::string> triggerNameMap_pp {
 //     {"MBD N&S >= 1",          "MBD_NandS_geq_1"},
//      {"Photon 3 GeV + MBD NS >= 1","Photon_3_GeV_plus_MBD_NS_geq_1"},
      {"Photon 4 GeV + MBD NS >= 1","Photon_4_GeV_plus_MBD_NS_geq_1"}
//      {"Photon 5 GeV + MBD NS >= 1","Photon_5_GeV_plus_MBD_NS_geq_1"}
  };
  std::map<int,std::string>* activeTriggerNameMap_pp{nullptr};

  std::map<int, std::string> triggerNameMapAuAu = {
//        {10, "MBD_NS_geq_2"},
//        {11, "MBD_NS_geq_1"},
//        {12, "MBD_NS_geq_2_vtx_lt_10"},
//        {13, "MBD_NS_geq_2_vtx_lt_30"},
        {14, "MBD_NS_geq_2_vtx_lt_150"}
//        {15, "MBD_NS_geq_1_vtx_lt_10"},
//        {16, "photon_6_plus_MBD_NS_geq_2_vtx_lt_10"},
//        {17, "photon_8_plus_MBD_NS_geq_2_vtx_lt_10"},
//        {18, "photon_10_plus_MBD_NS_geq_2_vtx_lt_10"},
//        {19, "photon_12_plus_MBD_NS_geq_2_vtx_lt_10"},
//        {20, "photon_6_plus_MBD_NS_geq_2_vtx_lt_150"},
//        {21, "photon_8_plus_MBD_NS_geq_2_vtx_lt_150"},
//        {22, "photon_10_plus_MBD_NS_geq_2_vtx_lt_150"},
//        {23, "photon_12_plus_MBD_NS_geq_2_vtx_lt_150"}
  };

  // Bit helpers
  std::string bitsetToList(uint64_t word) const;
  inline std::vector<int> extractTriggerBits(uint64_t b_gl1_scaledvec, [[maybe_unused]]int entry) {
    std::vector<int> trig_bits;
    std::bitset<64> bits(b_gl1_scaledvec);
    for (unsigned int bit = 0; bit < 64; bit++) if (((b_gl1_scaledvec >> bit) & 0x1U) == 0x1U) trig_bits.push_back(bit);
    return trig_bits;
  }
  inline bool checkTriggerCondition(const std::vector<int> &trig_bits, int inputBit) {
    for (const int &bit : trig_bits) if (bit == inputBit) return true;
    return false;
  }

  // ---------- Isolation WP config (PPG-12) ----------
  double m_isoA       {1.08128};   // GeV
  double m_isoB       {0.0299107}; // per GeV
  double m_isoGap     {1.0};       // GeV (sideband offset)
  double m_isoConeR   {0.3};       // ΔR
  double m_isoTowMin  {0.060};     // GeV
  double m_isoFracLegacy {0.10};   // legacy fractional definition (if needed)
  double m_hadFrac    {0.10};
  double m_minBackToBack {7.0*M_PI/8.0};


  // ---------- Internal helpers for photon pT(/centrality) binning & counters -------
  // NOTE: findPtBin + m_gammaPtBins define the photon pT slicing used everywhere (including JES calibration).
  int  findPtBin(double pt) const;                 // returns -1 if out-of-range
  int  findCentBin(int cent) const;                // returns -1 if out-of-range
  std::string suffixForBins(int ptIdx, int centIdx) const;

  TH1I* getOrBookCountHist(const std::string& trig,
                             const std::string& base,
                             int ptIdx, int centIdx);

    // 1D spectra (pT^gamma-sliced via suffixForBins)
    TH1F* getOrBookXJHist (const std::string& trig, int ptIdx, int centIdx);
    TH1F* getOrBookIsoHist(const std::string& trig, int ptIdx, int centIdx);

    // NEW: pure isolation QA (independent of preselection/tightness)
    // Component spectra use the same slicing rules as h_Eiso.
    TH1F* getOrBookIsoPartHist(const std::string& trig,
                               const std::string& base,
                               const std::string& xAxisTitle,
                               int ptIdx, int centIdx);

    // NEW: 2-bin (PASS/FAIL) isolation decision counter hist (same slicing rules)
    TH1I* getOrBookIsoDecisionHist(const std::string& trig, int ptIdx, int centIdx);


  // JES calibration primitives (filled only for iso∧tight γ + matched recoil jet)
  TH1F* getOrBookJet1PtHist(const std::string& trig, int ptIdx, int centIdx);
  TH1F* getOrBookJet2PtHist(const std::string& trig, int ptIdx, int centIdx);
  TH1F* getOrBookAlphaHist (const std::string& trig, int ptIdx, int centIdx);

  // Correlation-preserving TH3s (pT^gamma is an axis; centrality suffix only)
  TH3F* getOrBookJES3_xJ_alphaHist      (const std::string& trig, int centIdx);
  TH3F* getOrBookJES3_jet1Pt_alphaHist  (const std::string& trig, int centIdx);

   // ---------------------- Matching-QA (centrality suffix only) ----------------------
   // NOTE: pT^gamma is an axis (m_gammaPtBins). Do NOT pT-slice these in the name.
   TH2F*     getOrBookMatchDphiVsPtGamma       (const std::string& trig, int centIdx);
   TH2F*     getOrBookMatchStatusVsPtGamma     (const std::string& trig, int centIdx);
   TH2F*     getOrBookMatchMaxDphiVsPtGamma    (const std::string& trig, int centIdx);

   // Optional, still low-bloat:
   TProfile* getOrBookNRecoilJetsVsPtGamma     (const std::string& trig, int centIdx);
   TH2F*     getOrBookRecoilIsLeadingVsPtGamma (const std::string& trig, int centIdx);

   // Shower-shape histogram booker:
   //    h_ss_<varKey>_<tagKey> + suffixForBins(pT[,cent])
   TH1F* getOrBookSSHist(const std::string& trig,
                              const std::string& varKey,
                              const std::string& tagKey,
                              int ptIdx, int centIdx);

    void  fillIsoSSTagCounters(const std::string& trig,
                                const RawCluster* clus,
                                const SSVars& v,
                                double pt_gamma,
                                int centIdx,
                                PHCompositeNode* topNode);

    // ---------------------- Pure jet QA (photon-independent) ----------------------
    void fillInclusiveJetQA(const std::vector<std::string>& activeTrig,
                            int centIdx,
                            const std::string& rKey);

    // ---------------- Jet-only QA for jets used in the gamma+jet step -------------
    void fillSelectedJetQA(const std::vector<std::string>& activeTrig,
                           int ptIdx, int centIdx,
                           const std::string& rKey,
                           const Jet* jet1,
                           const Jet* jet2);

    // --------------- Generic jet-QA histogram bookers (uses suffixForBins) --------
    TH1I* getOrBookJetQA1I(const std::string& trig,
                           const std::string& base,
                           const std::string& xAxisTitle,
                           const std::string& rKey,
                           int ptIdx, int centIdx,
                           int nbins, double xmin, double xmax);

    TH1F* getOrBookJetQA1F(const std::string& trig,
                           const std::string& base,
                           const std::string& xAxisTitle,
                           const std::string& rKey,
                           int ptIdx, int centIdx,
                           int nbins, double xmin, double xmax);

    TH2F* getOrBookJetQA2F(const std::string& trig,
                           const std::string& base,
                           const std::string& xAxisTitle,
                           const std::string& yAxisTitle,
                           const std::string& rKey,
                           int ptIdx, int centIdx,
                           int nxbins, double xmin, double xmax,
                           int nybins, double ymin, double ymax);

    SSVars makeSSFromPhoton(const PhotonClusterv1* pho, double pt_gamma) const;
    void   processCandidates(PHCompositeNode* topNode,
                                 const std::vector<std::string>& activeTrig);
  };

#endif  // RECOILJETS_H
