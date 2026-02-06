// Tell Emacs this is C++   -*- C++ -*-
#ifndef RECOILJETS_H
#define RECOILJETS_H

// ============================================================================
// RecoilJets
//   - Photon (tight/iso) selection + recoil jet matching + xJ/alpha/JES outputs
//   - Supports pp, Au+Au, and simulation
//   - IMPORTANT: jet-dependent outputs are radius-tagged (r02, r04, ...)
// ============================================================================

// ---------------------------------------------------------------------------
// Fun4All / PHOOL
// ---------------------------------------------------------------------------
#include <fun4all/SubsysReco.h>
#include <calotrigger/TriggerAnalyzer.h>
#include <phool/PHCompositeNode.h>

// ---------------------------------------------------------------------------
// ROOT
// ---------------------------------------------------------------------------
#include <TFile.h>
#include <TTree.h>
#include <TObject.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2F.h>
#include <TH2Poly.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TProfile3D.h>
#include <TLorentzVector.h>
#include <TVector2.h>

// ---------------------------------------------------------------------------
// sPHENIX object headers (used as pointer types throughout)
// NOTE: These are kept here for "drop-in" builds since the .cc relies on
//       this header to provide these definitions in multiple places.
// ---------------------------------------------------------------------------
#include <calobase/RawClusterContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/RawTowerGeom.h>
#include <globalvertex/GlobalVertexMap.h>
#include <jetbase/JetContainer.h>

// You currently include PhotonClusterBuilder in the header; keep for drop-in.
// (This is not ideal for portability, but preserves your existing build.)
#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/coresoftware_local/offline/packages/CaloReco/PhotonClusterBuilder.h"

// ---------------------------------------------------------------------------
// STL
// ---------------------------------------------------------------------------
#include <algorithm>
#include <array>
#include <bitset>
#include <cstdint>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

// ---------------------------------------------------------------------------
// Helper alias: maps histogram name -> ROOT object*
// ---------------------------------------------------------------------------
using HistMap = std::map<std::string, TObject*>;

// Forward declarations (pointers only in the interface)
class CentralityInfo;
class Gl1Packet;
class EventHeader;
class PHG4TruthInfoContainer;
class PHG4Particle;
class PHG4VtxPoint;

class GlobalVertex;
class RawCluster;
class PhotonClusterv1;
class Jet;

// g4eval: used for truth↔reco association of EMCal clusters
class CaloRawClusterEval;

// HepMC forward declarations (used by unified truth signal classifier)
namespace HepMC
{
  class GenEvent;
  class GenParticle;
}

// ============================================================================
// Photon ID working points / helper functions  (PPG12 Table 4)
// ============================================================================
// NOTE:
//   .cc uses:  using namespace PhoIDCuts;
//   PPG12 defines:
//     - Preselection cuts applied first (fail => rejected; does NOT enter A–B–C–D)
//     - Tight cuts applied only after preselection
//     - Non-tight: passes preselection AND fails >=2 of the 5 tight cuts
//   Also: for photons in this context, pT^gamma ≡ ET^gamma (use pt_gamma everywhere).
// ============================================================================
namespace PhoIDCuts
{
  // ------------------- Pre-selection (Table 4) -------------------
  // Applied to EMCal photon-cluster candidates BEFORE tight/non-tight classification.
  inline constexpr double PRE_E11E33_MAX = 0.98;  // E11/E33 < 0.98
  inline constexpr double PRE_ET1_MIN    = 0.60;  // 0.6 < et1 < 1.0
  inline constexpr double PRE_ET1_MAX    = 1.00;
  inline constexpr double PRE_E32E35_MIN = 0.80;  // 0.8 < E32/E35 < 1.0
  inline constexpr double PRE_E32E35_MAX = 1.00;
  inline constexpr double PRE_WETA_MAX   = 0.60;  // weta_cogx < 0.6  (only η width capped in preselection)

  // ------------------- Tight ID (Table 4) -------------------
  // Applied ONLY to candidates that PASS preselection.
  // Widths use a pT(=ET) dependent upper bound:
  inline constexpr double TIGHT_W_LO = 0.0;       // 0 < w < (0.15 + 0.006 * ET^gamma)

  inline constexpr double TIGHT_E11E33_MIN = 0.40; // 0.4 < E11/E33 < 0.98
  inline constexpr double TIGHT_E11E33_MAX = 0.98;

  inline constexpr double TIGHT_ET1_MIN    = 0.90; // 0.9 < et1 < 1.0
  inline constexpr double TIGHT_ET1_MAX    = 1.00;

  inline constexpr double TIGHT_E32E35_MIN = 0.92; // 0.92 < E32/E35 < 1.0
  inline constexpr double TIGHT_E32E35_MAX = 1.00;

  // Open interval helper: (lo < x < hi)
  inline bool in_open_interval(double x, double lo, double hi)
  {
    return (std::isfinite(x) && std::isfinite(lo) && std::isfinite(hi) && (x > lo) && (x < hi));
  }

  // Tight upper width bound: w_hi(ET^gamma) = 0.15 + 0.006 * ET^gamma
  // In your code, pass pt_gamma (since pT^gamma == ET^gamma here).
  inline double tight_w_hi(double pt_gamma)
  {
    if (!std::isfinite(pt_gamma) || pt_gamma <= 0.0)
      return std::numeric_limits<double>::quiet_NaN();
    return 0.15 + 0.006 * pt_gamma;
  }

} // namespace PhoIDCuts

// ============================================================================
// RecoilJets analysis module
// ============================================================================
class RecoilJets : public SubsysReco
{
public:
  // -------------------------------------------------------------------------
  // Lightweight categorization used by photon ID QA
  // -------------------------------------------------------------------------
  enum class TightTag : std::uint8_t
  {
    kPreselectionFail = 0,
    kTight            = 1,
    kNonTight         = 2,  // >=2 fails
    kNeither          = 3   // exactly 1 fail
  };

  enum class EventReject : std::uint8_t
  {
    None    = 0,
    Trigger = 1,
    Vz      = 2
  };

  // Shower-shape variables extracted from PhotonClusterv1
  struct SSVars
  {
    double pt_gamma      = 0.0;
    double weta_cogx     = 0.0;
    double wphi_cogx     = 0.0;
    double et1           = 0.0;
    double e11_over_e33  = 0.0;
    double e32_over_e35  = 0.0;
  };

  // Per-(trigger,slice) category counters printed in End()
  struct CatStat
  {
    long long n_iso_tight        = 0;
    long long n_nonIso_tight     = 0;
    long long n_iso_nonTight     = 0;
    long long n_nonIso_nonTight  = 0;
    long long seen               = 0;
  };

  // Job-wide bookkeeping counters printed in printCutSummary()
  struct Bookkeeping
  {
    // Event-level
    long long evt_seen         = 0;
    long long evt_fail_trigger = 0;
    long long evt_fail_vz      = 0;
    long long evt_accepted     = 0;

    // Photon-level
    long long pho_total            = 0;
    long long pho_early_E          = 0;
    long long pho_eta_fail         = 0;
    long long pho_etbin_out        = 0;
    long long pho_noRC             = 0;
    long long pho_reached_pre_iso  = 0;

    // Preselection breakdown
    long long pre_pass              = 0;
    long long pre_fail_weta         = 0;
    long long pre_fail_et1_low      = 0;
    long long pre_fail_et1_high     = 0;
    long long pre_fail_e11e33_high  = 0;
    long long pre_fail_e32e35_low   = 0;
    long long pre_fail_e32e35_high  = 0;

    // Tight classification breakdown
    long long tight_tight       = 0;
    long long tight_neither     = 0;
    long long tight_nonTight    = 0;

    // Tight per-cut fails
    long long tight_fail_weta   = 0;
    long long tight_fail_wphi   = 0;
    long long tight_fail_et1    = 0;
    long long tight_fail_e11e33 = 0;
    long long tight_fail_e32e35 = 0;

    // Isolation
    long long iso_pass          = 0;
    long long iso_fail          = 0;
  };

  // -------------------------------------------------------------------------
  // Jet radius mapping (key + node names per dataset)
  // -------------------------------------------------------------------------
  struct JetRadiusDef
  {
    std::string key;      // e.g. "r02"
    std::string pp_node;  // jet node name used in pp
    std::string aa_node;  // jet node name used in Au+Au
  };

  // Radii to run in parallel. These are used in fetchNodes(), jet QA, matching.
  //
  // NOTE: Replace node names here if your DST uses different names.
    inline static const std::array<JetRadiusDef, 2> kJetRadii = {{
      // key   pp_node              aa_node
      { "r02", "AntiKt_Tower_r02",  "AntiKt_Tower_r02_Sub1"  },
      { "r04", "AntiKt_Tower_r04",  "AntiKt_Tower_r04_Sub1"  }
    }};

  // Trigger maps:
  //  - pp:   DB/GL1 name -> short key (directory name)
  //  - AuAu: bit index   -> short key (directory name)
  //
  // NOTE: These are placeholders if you don't already have your project-specific
  //       maps. Keep/restore your exact trigger mappings here.
  inline static const std::vector<std::pair<std::string, std::string>> triggerNameMap_pp = {
    //     {"MBD N&S >= 1",          "MBD_NandS_geq_1"},
    //      {"Photon 3 GeV + MBD NS >= 1","Photon_3_GeV_plus_MBD_NS_geq_1"},
    {"Photon 4 GeV + MBD NS >= 1","Photon_4_GeV_plus_MBD_NS_geq_1"}
    //      {"Photon 5 GeV + MBD NS >= 1","Photon_5_GeV_plus_MBD_NS_geq_1"}
  };
    
  inline static const std::vector<std::pair<int, std::string>> triggerNameMapAuAu = {
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

  // -------------------------------------------------------------------------
  // Construction / Fun4All hooks
  // -------------------------------------------------------------------------
  explicit RecoilJets(const std::string& outFile);
  ~RecoilJets() override;

  int Init(PHCompositeNode* topNode) override;
  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int ResetEvent(PHCompositeNode* topNode) override;
  int Reset(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  // -------------------------------------------------------------------------
  // Configuration helpers (keep names stable; safe defaults)
  // -------------------------------------------------------------------------
  void setDataType(const std::string& s)
  {
      // Supported user-facing keys (also overridden by env RJ_DATASET / RJ_IS_SIM)
      //   "isSim"/"sim"  -> sim
      //   "isAuAu"/"auau"-> Au+Au
      //   "ispp"/"pp"    -> pp
      const std::string t = s;
      if (t == "isSim" || t == "sim") { m_isSim = true;  m_isAuAu = false; }
      if (t == "isAuAu" || t == "auau"|| t == "aa") { m_isSim = false; m_isAuAu = true; }
      if (t == "ispp"  || t == "pp")  { m_isSim = false; m_isAuAu = false; }
  }

  // Optional: limit active jet radii (keys like "r02","r04"). Empty => all kJetRadii.
  void setActiveJetRKeys(const std::vector<std::string>& keys) { m_activeJetRKeys = keys; }

  void setPhotonEtaAbsMax(double v) { m_etaAbsMax = v; }
  void setMinJetPt(double v)        { m_minJetPt = v; }
  void setMinBackToBack(double v)   { m_minBackToBack = v; }

  void setUseVzCut(bool on, double vz = 60.0)
  {
      m_useVzCut = on;
      m_vzCut    = static_cast<float>(vz);
  }
  void setGammaPtBins(const std::vector<double>& bins)
          {
            // Canonical pT^gamma binning for ALL photon-binned histograms (reco + truth):
            //   [15,17,19,21,23,26,35] GeV  (6 bins, start at 15 GeV)
            //
            // NOTE: Unfolding histograms book their own extended pT^gamma binning:
            //   reco : [10,15] + (canonical bins) + [35,40]
            //   truth: [5,10] + [10,15] + (canonical bins) + [35,40]
            if (!bins.empty())
            {
              m_gammaPtBins = bins;
            }
  }

  // Phase-1 YAML knobs (matching thresholds)
  void setPhoMatchDRMax(double v) { m_phoMatchDRMax = v; }
  void setJetMatchDRMax(double v) { m_jetMatchDRMax = v; }

  // Phase-1 YAML knobs (unfolding explicit bin edges)
  void setUnfoldRecoPhotonPtBins(const std::vector<double>& bins)  { m_unfoldRecoPhotonPtBins = bins; }
  void setUnfoldTruthPhotonPtBins(const std::vector<double>& bins) { m_unfoldTruthPhotonPtBins = bins; }
  void setUnfoldJetPtBins(const std::vector<double>& bins)         { m_unfoldJetPtBins = bins; }
  void setUnfoldXJBins(const std::vector<double>& bins)            { m_unfoldXJBins = bins; }

  // Analysis provenance stamping (written once into output ROOT by RecoilJets.cc)
  void setAnalysisConfigYAML(const std::string& yamlText, const std::string& tag)
  {
        m_analysisConfigYAMLText = yamlText;
        m_analysisConfigTag      = tag;
  }


  void setCentEdges(const std::vector<int>& edges)     { m_centEdges = edges; }

  // Isolation WP (implemented in .cc)
  void setIsolationWP(double aGeV, double bPerGeV,
                              double sideGapGeV, double coneR, double towerMin);
  void setIsSlidingIso(bool on) { m_isSlidingIso = on; }

    // Photon ID cuts (PPG12 Table 4): allow YAML override while preserving baseline defaults
  void setPhotonIDCuts(double pre_e11e33_max,
                           double pre_et1_min,
                           double pre_et1_max,
                           double pre_e32e35_min,
                           double pre_e32e35_max,
                           double pre_weta_max,
                           double tight_w_lo,
                           double tight_w_hi_intercept,
                           double tight_w_hi_slope,
                           double tight_e11e33_min,
                           double tight_e11e33_max,
                           double tight_et1_min,
                           double tight_et1_max,
                           double tight_e32e35_min,
                           double tight_e32e35_max)
  {
        m_phoid_pre_e11e33_max = pre_e11e33_max;
        m_phoid_pre_et1_min    = pre_et1_min;
        m_phoid_pre_et1_max    = pre_et1_max;
        m_phoid_pre_e32e35_min = pre_e32e35_min;
        m_phoid_pre_e32e35_max = pre_e32e35_max;
        m_phoid_pre_weta_max   = pre_weta_max;

        m_phoid_tight_w_lo           = tight_w_lo;
        m_phoid_tight_w_hi_intercept = tight_w_hi_intercept;
        m_phoid_tight_w_hi_slope     = tight_w_hi_slope;

        m_phoid_tight_e11e33_min = tight_e11e33_min;
        m_phoid_tight_e11e33_max = tight_e11e33_max;

        m_phoid_tight_et1_min    = tight_et1_min;
        m_phoid_tight_et1_max    = tight_et1_max;

        m_phoid_tight_e32e35_min = tight_e32e35_min;
        m_phoid_tight_e32e35_max = tight_e32e35_max;
  }

    // EventDisplay diagnostics payload (offline rendering; independent of Verbosity()).
    // When enabled, a compact per-event-per-radius TTree ("EventDisplayTree") is written to the output ROOT file.
    void enableEventDisplayDiagnostics(bool on = true) { m_evtDiagEnabled = on; }
    void setEventDisplayDiagnosticsMaxPerBin(int n)    { m_evtDiagMaxPerBin = n; }

  // -------------------------------------------------------------------------
  // Read-only state access
  // -------------------------------------------------------------------------
  const Bookkeeping& bookkeeping() const { return m_bk; }
  EventReject lastReject() const         { return m_lastReject; }

private:
  // -------------------------------------------------------------------------
  // Internal helpers (nodes, event selection)
  // -------------------------------------------------------------------------
  bool jetRKeyActive(const std::string& rKey) const
  {
    if (m_activeJetRKeys.empty()) return true;
    return (std::find(m_activeJetRKeys.begin(), m_activeJetRKeys.end(), rKey) != m_activeJetRKeys.end());
  }

  bool fetchNodes(PHCompositeNode* topNode);
  bool firstEventCuts(PHCompositeNode* topNode, std::vector<std::string>& activeTrig);
  void createHistos_Data();

    void fillUnfoldResponseMatrixAndTruthDistributions(
          const std::vector<std::string>& activeTrig,
          const std::string& rKey,
          const int effCentIdx_M,
          const double leadPtGamma,
          const double tPt,
          const double tEta,
          const double tPhi,
          const std::vector<const Jet*>& recoJetsFid,
          const std::vector<char>& recoJetsFidIsRecoil,
          const Jet* recoil1Jet);

      
    void fillRecoTruthJES3MatchingQA(const std::vector<std::string>& activeTrig,
                                      const std::string& rKey,
                                      const int effCentIdx_M,
                                      const double leadPtGamma,
                                      const double xJ,
                                      const double alpha,
                                      const double tPt,
                                      const double tEta,
                                      const double tPhi,
                                      const Jet* recoil1Jet);

    bool runLeadIsoTightPhotonJetLoopAllRadii(
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
            const double tEta,
            const double tPhi,
            PHG4TruthInfoContainer* truth);

      
    bool runLeadIsoTightPhotonJetMatchingAndUnfolding(
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
            const double tEta,
            const double tPhi,
            PHG4TruthInfoContainer* truth);

  void fillPureIsolationQA(PHCompositeNode* topNode,
                               const std::vector<std::string>& activeTrig,
                               const PhotonClusterv1* pho,
                               const RawCluster* rc,
                               const int ptIdx,
                               const int centIdx,
                               const double pt_gamma);

  void fillTruthSigABCDLeakageCounters(PHCompositeNode* topNode,
                                         const std::vector<std::string>& activeTrig,
                                         const int centIdx);

  void processCandidates(PHCompositeNode* topNode, const std::vector<std::string>& activeTrig);

  bool getCentralitySlice(int& lo, int& hi, std::string& tag) const;

  // Au+Au scaled trigger helper (safe default behavior)
  static std::bitset<64> extractTriggerBits(std::uint64_t scaledVec, int /*eventNumber*/)
  {
    return std::bitset<64>(scaledVec);
  }
  static bool checkTriggerCondition(const std::bitset<64>& bits, int bit)
  {
    if (bit < 0 || bit >= 64) return false;
    return bits.test(static_cast<std::size_t>(bit));
  }

  // -------------------------------------------------------------------------
  // Photon ID helpers
  // -------------------------------------------------------------------------
  SSVars makeSSFromPhoton(const PhotonClusterv1* pho, double pt_gamma) const;
  bool   passesPhotonPreselection(const SSVars& v);
  TightTag classifyPhotonTightness(const SSVars& v);


  // Isolation helpers
  double eiso(const RawCluster* clus, PHCompositeNode* topNode) const;
  bool   isIsolated(const RawCluster* clus, double et_gamma, PHCompositeNode* topNode) const;
  bool   isNonIsolated(const RawCluster* clus, double et_gamma, PHCompositeNode* topNode) const;

    // Unified truth-MC signal definition for "isolated prompt photon" (SIM only)
    // Definition
    //   |eta| < 0.7, PID=22, status=1 (final state),
    //   prompt classification via CaloAna photon_type logic:
    //     - walk back photon-in/photon-out vertices
    //     - direct=1 if 2->2 with |pdg|<=22 on all legs
    //     - frag  =2 if 1->2 with |incoming pdg|<=11 and outgoing contains incoming pid (and photon)
    //   and truth isolation ETiso_truth < 4 GeV where (CaloAna truth-iso):
    //     ETiso = sum_{ΔR<0.3} Et(final-state)  -  sum_{ΔR<0.001} Et(final-state)
    //     (the ΔR<0.001 subtraction removes the photon itself, and any ultra-merged pieces).
  bool isTruthPromptIsolatedSignalPhoton(const HepMC::GenEvent* evt,
                                          const HepMC::GenParticle* pho,
                                          double& isoEt) const;
    
  // Unified truth→reco photon matching using CaloRawClusterEval.
  // Match requirements:
  //   - reco |eta| < 0.7, reco pT > 5 GeV
  //   - ΔR(truth,reco) < 0.05
  //   - reco cluster’s BEST-MATCHED truth primary (by deposited energy) has same barcode as truth photon
  // Chooses the candidate with the largest energy contribution (tie-breaker: smallest ΔR).
  bool findRecoPhotonMatchedToTruthSignal(const HepMC::GenEvent* evt,
                                           const HepMC::GenParticle* truthPho,
                                           CaloRawClusterEval& clustereval,
                                           const RawCluster*& recoPho,
                                           double& recoPt,
                                           double& recoEta,
                                           double& recoPhi,
                                           double& drBest,
                                           float& eContribBest) const;


  // -------------------------------------------------------------------------
  // Binning helpers
  // -------------------------------------------------------------------------
  int         findPtBin(double pt) const;
  int         findCentBin(int cent) const;
  std::string suffixForBins(int ptIdx, int centIdx) const;

  // EventDisplay categories (used by the diagnostics payload stored in EventDisplayTree)
  enum class EventDisplayCat : int { NUM = 0, MissA = 1, MissB = 2 };

  // EventDisplay diagnostics payload (offline rendering; independent of Verbosity()).
  void        initEventDisplayDiagnosticsTree();
  void        resetEventDisplayDiagnosticsBuffers();
  bool        eventDisplayDiagnosticsNeed(const std::string& rKey, int ptBin, EventDisplayCat cat);
  void        appendEventDisplayDiagnosticsFromJet(const Jet* jet,
                                                      std::vector<int>& calo,
                                                      std::vector<int>& ieta,
                                                      std::vector<int>& iphi,
                                                      std::vector<float>& eta,
                                                      std::vector<float>& phi,
                                                      std::vector<float>& et,
                                                      std::vector<float>& e) const;
    
    void        fillEventDisplayDiagnostics(const std::string& rKey,
                                               int ptBin,
                                               EventDisplayCat cat,
                                               double truthGammaPt,
                                               double truthGammaPhi,
                                               double recoGammaPt,
                                               double recoGammaEta,
                                               double recoGammaPhi,
                                               const Jet* selectedRecoilJet,
                                               const Jet* recoTruthBest,
                                               const Jet* truthLeadRecoilJet);
  // -------------------------------------------------------------------------
  // Histogram utilities (bookers)
  // -------------------------------------------------------------------------
  TH1I* getOrBookCountHist(const std::string& trig,
                             const std::string& base,
                             int ptIdx, int centIdx);

  // Isolation spectra (reco)
  TH1F* getOrBookIsoHist(const std::string& trig, int ptIdx, int centIdx);
  TH1F* getOrBookIsoPartHist(const std::string& trig,
                                     const std::string& base,
                                     const std::string& xAxisTitle,
                                     int ptIdx, int centIdx);
  TH1F* getOrBookPtGammaHist(const std::string& trig,
                                 const std::string& base,
                                 int centIdx);
  TH1I* getOrBookIsoDecisionHist(const std::string& trig, int ptIdx, int centIdx);

  // Event-level photon multiplicity diagnostic:
  // N = number of reco photon candidates passing the SAME tight+iso+fiducial cuts
  // used for the leading photon selection (filled once per event; binned by leading-photon pT bin).
  TH1I* getOrBookNIsoTightPhoCandHist(const std::string& trig, int ptIdx, int centIdx);

  TH2F* getOrBookIsoCompareHist(const std::string& trig, int ptIdx, int centIdx);

  // -------------------------------------------------------------------------
  // SIM ONLY: matched truth-signal → reco ABCD leakage counters (per pT[/cent] slice)
  //   bins: 1=A, 2=B, 3=C, 4=D
  // -------------------------------------------------------------------------
  TH1I* getOrBookSigABCDLeakageHist(const std::string& trig, int ptIdx, int centIdx);

  // Isolation QA (truth, SIM only)
  // These are NOT sliced; they live in the trigger directory (SIM => /SIM/).
  TH1F* getOrBookTruthIsoHist(const std::string& trig,
                                const std::string& name,
                                int nbins, double xmin, double xmax);

  TH1I* getOrBookTruthIsoDecisionHist(const std::string& trig,
                                        const std::string& name);

  // Shower-shape distributions
  TH1F* getOrBookSSHist(const std::string& trig,
                        const std::string& varKey,
                        const std::string& tagKey,
                        int ptIdx, int centIdx);

  // Physics outputs (already radius-tagged in your .cc)
  TH1F* getOrBookXJHist(const std::string& trig,
                          const std::string& rKey,
                          int ptIdx, int centIdx);
  TH1F* getOrBookJet1PtHist(const std::string& trig,
                              const std::string& rKey,
                              int ptIdx, int centIdx);
  TH1F* getOrBookJet2PtHist(const std::string& trig,
                              const std::string& rKey,
                              int ptIdx, int centIdx);
  TH1F* getOrBookAlphaHist(const std::string& trig,
                             const std::string& rKey,
                             int ptIdx, int centIdx);

  // -------------------------------------------------------------------------
  // inclusive γ–jet unfolding histograms (ATLAS-style 2D in (pT^γ, xJγ))
  //
  //  - Reco (DATA + MC):
  //      h2_unfoldReco_pTgamma_xJ_incl_<rKey><centSuffix>
  //
  //  - Truth (MC):
  //      h2_unfoldTruth_pTgamma_xJ_incl_<rKey><centSuffix>
  //
  //  - Response (MC, matched truth↔reco pairs):
  //      h2_unfoldResponse_pTgamma_xJ_incl_<rKey><centSuffix>
  //    stored as global-bin(truth) vs global-bin(reco), using TH2::FindBin().
  //
  //  - Fakes/Misses (MC, for closure / response completeness):
  //      h2_unfoldRecoFakes_pTgamma_xJ_incl_<rKey><centSuffix>
  //      h2_unfoldTruthMisses_pTgamma_xJ_incl_<rKey><centSuffix>
  // -------------------------------------------------------------------------
  TH2F* getOrBookUnfoldRecoPtXJIncl      (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookUnfoldTruthPtXJIncl     (const std::string& trig, const std::string& rKey, int centIdx);

  // inclusive |Δphi(gamma,jet)| per recoil jet that passes pT+eta+recoil Δphi cuts
  TH2F* getOrBookUnfoldRecoPtDphiIncl    (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookUnfoldTruthPtDphiIncl   (const std::string& trig, const std::string& rKey, int centIdx);

  TH2F* getOrBookUnfoldResponsePtXJIncl  (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookUnfoldRecoFakesPtXJIncl (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookUnfoldTruthMissesPtXJIncl(const std::string& trig, const std::string& rKey, int centIdx);

  // photon-only unfolding (for N_gamma normalization): DATA+SIM reco, SIM truth
  TH1F* getOrBookUnfoldRecoPhoPtGamma        (const std::string& trig, int centIdx);
  TH1F* getOrBookUnfoldTruthPhoPtGamma       (const std::string& trig, int centIdx);
  TH2F* getOrBookUnfoldResponsePhoPtGamma    (const std::string& trig, int centIdx);
  TH1F* getOrBookUnfoldRecoPhoFakesPtGamma   (const std::string& trig, int centIdx);
  TH1F* getOrBookUnfoldTruthPhoMissesPtGamma (const std::string& trig, int centIdx);

  // unfolding QA helpers (SIM only): explicit matched distributions + type-split fakes/misses
  TH2F* getOrBookUnfoldTruthMatchedPtXJIncl      (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookUnfoldRecoMatchedPtXJIncl       (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookUnfoldRecoFakesPtXJIncl_typeA   (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookUnfoldRecoFakesPtXJIncl_typeB   (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookUnfoldTruthMissesPtXJIncl_typeA (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookUnfoldTruthMissesPtXJIncl_typeB (const std::string& trig, const std::string& rKey, int centIdx);

  // jet-match QA (SIM only): match ΔR and pT response for matched recoil jets
  TH1F* getOrBookUnfoldJetMatchDR           (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookUnfoldJetPtResponsePtTruth (const std::string& trig, const std::string& rKey, int centIdx);

  // additional JES diagnostics (SIM only):
  //  - ALL matched fiducial jet pairs (no recoil/truth recoil-selection gating)
  //  - lead recoil jet1-only response + ΔR sanity
  //  - lead recoil jet1 scatter: pT(reco) vs pT(truth)
  TH2F* getOrBookUnfoldJetPtResponseAllPtTruth   (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadRecoilJetPtResponsePtTruth  (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadRecoilJetPtTruthPtReco      (const std::string& trig, const std::string& rKey, int centIdx);
  TH1F* getOrBookLeadRecoilJetMatchDR            (const std::string& trig, const std::string& rKey, int centIdx);

  // -------------------------------------------------------------------------
  // (SIM ONLY): JES3-style *leading truth recoil jet1* match bookkeeping vs truth pT^gamma
  //   Den  : truth leading recoil jet1 exists (truth recoil definition)
  //   Num  : Den + reco recoil jet1 matches truth jet1 (ΔR < 0.3)
  //   MissA: Den + some reco fid jet within ΔR < 0.3 of truth jet1, but Num failed
  //   MissB: Den + no reco fid jet within ΔR < 0.3 of truth jet1
  // -------------------------------------------------------------------------
  TH1F* getOrBookLeadTruthRecoilMatchDenPtGammaTruth    (const std::string& trig, const std::string& rKey, int centIdx);
  TH1F* getOrBookLeadTruthRecoilMatchNumPtGammaTruth    (const std::string& trig, const std::string& rKey, int centIdx);
  TH1F* getOrBookLeadTruthRecoilMatchMissA_PtGammaTruth (const std::string& trig, const std::string& rKey, int centIdx);
  TH1F* getOrBookLeadTruthRecoilMatchMissB_PtGammaTruth (const std::string& trig, const std::string& rKey, int centIdx);

  // NEW: MissA subtypes (SIM-only)
  //   MissA1: truth-matched reco jet passes recoil definition → competitor/ordering
  //   MissA2: truth-matched reco jet fails recoil definition  → gate-exclusion
  TH1F* getOrBookLeadTruthRecoilMatchMissA1_PtGammaTruth (const std::string& trig, const std::string& rKey, int centIdx);
  TH1F* getOrBookLeadTruthRecoilMatchMissA2_PtGammaTruth (const std::string& trig, const std::string& rKey, int centIdx);
  TH1I* getOrBookLeadTruthRecoilMatchMissA2_Cutflow      (const std::string& trig, const std::string& rKey, int centIdx);

  // -------------------------------------------------------------------------
  // (SIM ONLY): diagnostics for why reco xJ differs from truth-conditioned xJ
  //   Filled inside the same DEN / NUM / MissA / MissB classification block.
  //
  //   NOTE: all are radius-tagged and use centrality-only suffix (like the
  //   existing LeadTruthRecoilMatch bookkeeping).
  // -------------------------------------------------------------------------

  // (A1) pT(recoilJet1^reco) vs pT(truth-leading recoil jet), split by class
  TH2F* getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtTruthLead_num    (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtTruthLead_missA  (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtTruthLead_missA1 (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtTruthLead_missA2 (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtTruthLead_missB  (const std::string& trig, const std::string& rKey, int centIdx);

  // (A2) pT(recoilJet1^reco) vs pT(reco jet matched to truth-leading recoil jet), for NUM / MissA
  TH2F* getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtRecoTruthMatch_num    (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtRecoTruthMatch_missA  (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtRecoTruthMatch_missA1 (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtRecoTruthMatch_missA2 (const std::string& trig, const std::string& rKey, int centIdx);

  // (B3) |Δphi(γ^truth, recoilJet1^reco)| vs pTγ,truth, split by class
  TH2F* getOrBookLeadTruthRecoilMatchDphiRecoJet1VsPtGammaTruth_num   (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadTruthRecoilMatchDphiRecoJet1VsPtGammaTruth_missA (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadTruthRecoilMatchDphiRecoJet1VsPtGammaTruth_missB (const std::string& trig, const std::string& rKey, int centIdx);

  // (B4) ΔR(recoilJet1^reco, truth-leading recoil jet) vs pTγ,truth, split by class
  TH2F* getOrBookLeadTruthRecoilMatchDRRecoJet1VsPtGammaTruth_num   (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadTruthRecoilMatchDRRecoJet1VsPtGammaTruth_missA (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadTruthRecoilMatchDRRecoJet1VsPtGammaTruth_missB (const std::string& trig, const std::string& rKey, int centIdx);

  // (C5) xJ(recoilJet1^reco) vs |Δphi(γ^truth, recoilJet1^reco)|, split by class
  TH2F* getOrBookLeadTruthRecoilMatchXJRecoJet1VsDphiRecoJet1_num   (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadTruthRecoilMatchXJRecoJet1VsDphiRecoJet1_missA (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadTruthRecoilMatchXJRecoJet1VsDphiRecoJet1_missB (const std::string& trig, const std::string& rKey, int centIdx);

  // -------------------------------------------------------------------------
  //  radius-tagged matching-QA bookers
  //   - centrality suffix only (Au+Au)
  //   - pT^gamma is an axis
  //   - name pattern: <base>_<rKey><centSuffix>
  // -------------------------------------------------------------------------
  TH2F*     getOrBookMatchStatusVsPtGamma     (const std::string& trig, const std::string& rKey, int centIdx);

  // per-jet cutflow status vs pT^gamma (NOT event-level; not iso/tight-conditioned)
  // name: h_jetcutflow_status_vs_pTgamma_<rKey><centSuffix>
  // y bins:
  //   1 = FailJetPt, 2 = FailJetEta, 3 = FailBackToBack, 4 = PassAll
  TH2F*     getOrBookJetCutflowStatusVsPtGamma(const std::string& trig, const std::string& rKey, int centIdx);

  TH2F*     getOrBookMatchMaxDphiVsPtGamma    (const std::string& trig, const std::string& rKey, int centIdx);
  TProfile* getOrBookNRecoilJetsVsPtGamma     (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F*     getOrBookMatchDphiVsPtGamma       (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F*     getOrBookRecoilIsLeadingVsPtGamma (const std::string& trig, const std::string& rKey, int centIdx);

  // -------------------------------------------------------------------------
  // radius-tagged JES3 bookers (data)
  //   - centrality suffix only (Au+Au)
  //   - name pattern: <base>_<rKey><centSuffix>
  // -------------------------------------------------------------------------
  TH3F* getOrBookJES3_xJ_alphaHist     (const std::string& trig, const std::string& rKey, int centIdx);
  TH3F* getOrBookJES3_jet1Pt_alphaHist (const std::string& trig, const std::string& rKey, int centIdx);

  // -------------------------------------------------------------------------
  // radius-tagged Jet13 + Balance3 bookers
  //   - no centrality suffix
  //   - name pattern: <base>_<rKey>
  // -------------------------------------------------------------------------
  TH3F*       getOrBookPho3TightIso       (const std::string& trig);               // photon baseline
  TH3F*       getOrBookJet13RecoilJet1    (const std::string& trig, const std::string& rKey);
  TProfile3D* getOrBookBalance3           (const std::string& trig, const std::string& rKey);

  // -------------------------------------------------------------------------
  // radius-tagged JES3 bookers (truth)
  //   - centrality suffix only (Au+Au)
  //   - name pattern: <base>_<rKey><centSuffix>
  // -------------------------------------------------------------------------
  // TRUTH reco-conditioned + reco↔truth jet1 match (your existing "works too well" truth)
  TH3F* getOrBookJES3Truth_xJ_alphaHist         (const std::string& trig, const std::string& rKey, int centIdx);

  // TRUTH reco-conditioned, but WITHOUT requiring reco jet1 ↔ truth jet1 ΔR match
  TH3F* getOrBookJES3TruthRecoCondNoJetMatch_xJ_alphaHist(const std::string& trig, const std::string& rKey, int centIdx);

  TH3F* getOrBookJES3Truth_jet1Pt_alphaHist     (const std::string& trig, const std::string& rKey, int centIdx);

  // PURE truth xJgamma distribution (no reco gating, no reco↔truth jet matching)
  TH3F* getOrBookJES3TruthPure_xJ_alphaHist     (const std::string& trig, const std::string& rKey, int centIdx);

  // -------------------------------------------------------------------------
  // reco-side JES3 subsets to diagnose fakes / wrong-jet assignment
  // -------------------------------------------------------------------------
  // RECO xJ,alpha for events where the LEADING reco iso∧tight photon is truth-signal tagged (photon match only)
  TH3F* getOrBookJES3RecoTruthPhoTagged_xJ_alphaHist(const std::string& trig, const std::string& rKey, int centIdx);

  // RECO xJ,alpha for truth-tagged pairs (photon match + reco jet1 matched to truth jet1)
  TH3F* getOrBookJES3RecoTruthTagged_xJ_alphaHist(const std::string& trig, const std::string& rKey, int centIdx);

  // -------------------------------------------------------------------------
  // Jet QA: generic bookers + fillers (already radius-tagged in your .cc)
  // -------------------------------------------------------------------------
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

  void fillInclusiveJetQA(const std::vector<std::string>& activeTrig,
                          int centIdx,
                          const std::string& rKey);

  void fillSelectedJetQA(const std::vector<std::string>& activeTrig,
                         int ptIdx, int centIdx,
                         const std::string& rKey,
                         const Jet* jet1,
                         const Jet* jet2);

  // -------------------------------------------------------------------------
  // SS + Iso category accounting
  // -------------------------------------------------------------------------
  void fillIsoSSTagCounters(const std::string& trig,
                            const RawCluster* clus,
                            const SSVars& v,
                            double pt_gamma,
                            int centIdx,
                            PHCompositeNode* topNode);

  void bumpHistFill(const std::string& trig, const std::string& hnameWithSuffix);
  void printCutSummary() const;

  // -------------------------------------------------------------------------
  // Member data
  // -------------------------------------------------------------------------
  std::string Outfile;
  TFile*      out    = nullptr;

  TriggerAnalyzer* trigAna = nullptr;

  // Dataset flags (these may be overridden at runtime by env in fetchNodes())
  bool m_isSim  = false;
  bool m_isAuAu = false;

  // Vertex / event selection
  GlobalVertex* m_vtx = nullptr;
  float m_vx = 0.0f;
  float m_vy = 0.0f;
  float m_vz = 0.0f;

  bool  m_useVzCut = true;
  float m_vzCut    = 30.0f;

  // Centrality
  int m_centBin = -1;                 // 0..99 (Au+Au), or -1 in pp
  std::vector<int> m_centEdges;       // centrality bin edges, e.g. {0,10,20,...,100}

  // Photon fiducial + binning
  double m_etaAbsMax = 0.7;           // photon |eta| cut

  std::vector<double> m_gammaPtBins = {15,17,19,21,23,26,35};  // canonical photon pT bin edges (6 bins, start at 15 GeV)

  // Photon ID cuts (PPG12 Table 4) defaults (match PhoIDCuts namespace baseline)
  double m_phoid_pre_e11e33_max = 0.98;
  double m_phoid_pre_et1_min    = 0.60;
  double m_phoid_pre_et1_max    = 1.00;
  double m_phoid_pre_e32e35_min = 0.80;
  double m_phoid_pre_e32e35_max = 1.00;
  double m_phoid_pre_weta_max   = 0.60;

  double m_phoid_tight_w_lo           = 0.0;
  double m_phoid_tight_w_hi_intercept = 0.15;
  double m_phoid_tight_w_hi_slope     = 0.006;

  double m_phoid_tight_e11e33_min = 0.40;
  double m_phoid_tight_e11e33_max = 0.98;

  double m_phoid_tight_et1_min    = 0.90;
  double m_phoid_tight_et1_max    = 1.00;

  double m_phoid_tight_e32e35_min = 0.92;
  double m_phoid_tight_e32e35_max = 1.00;


  // Phase-1 YAML knobs (matching thresholds)
  double m_phoMatchDRMax = 0.05;
  double m_jetMatchDRMax = 0.3;

  // Phase-1 YAML knobs (explicit unfolding bin edges)
  std::vector<double> m_unfoldRecoPhotonPtBins  = {10,15,17,19,21,23,26,35,40};
  std::vector<double> m_unfoldTruthPhotonPtBins = {5,10,15,17,19,21,23,26,35,40};
  std::vector<double> m_unfoldJetPtBins;
  std::vector<double> m_unfoldXJBins = {0.0,0.20,0.24,0.29,0.35,0.41,0.50,0.60,0.72,0.86,1.03,1.24,1.49,1.78,2.14,3.0};

  // Analysis provenance stamping (written once into output ROOT by RecoilJets.cc)
  std::string m_analysisConfigYAMLText = "";
  std::string m_analysisConfigTag      = "";
  bool        m_analysisConfigStamped  = false;


  // Isolation WP (PPG12 nominal)
  double m_isoA      = 1.08128;
  double m_isoB      = 0.0299107;
  double m_isoGap    = 1.0;
  double m_isoConeR  = 0.3;
  double m_isoTowMin = 0.0;
  bool   m_isSlidingIso = true;

  // Jet selection WP
  //  double m_minJetPt      = 5.0;
  double m_minJetPt      = 10.0;
  double m_minBackToBack = 7.0 * M_PI / 8.0;    // radians

  // Legacy "primary" reco jet key (still used for printing/overrides only)
  std::string m_xjRecoJetKey = "r04";

  // Nodes: photons / clusters
  RawClusterContainer* m_clus    = nullptr;
  RawClusterContainer* m_photons = nullptr;

  // Calo tower bundles (node cache)
  struct CaloBundle
  {
    TowerInfoContainer*   towers = nullptr;
    RawTowerGeomContainer* geom  = nullptr;
    double                sumEt  = 0.0;
  };

  std::vector<std::tuple<std::string, std::string, std::string>> m_caloInfo;
  std::map<std::string, CaloBundle> m_calo;

    // -------------------------------------------------------------------------
    // parallel jet containers by radius key
    // -------------------------------------------------------------------------
    std::map<std::string, JetContainer*> m_jets;
    std::map<std::string, JetContainer*> m_jetsRaw;

  // Optional: limit active jet radii (keys like "r02","r04"). Empty => all kJetRadii.
  std::vector<std::string> m_activeJetRKeys;

  // -------------------------------------------------------------------------
  // NEW: truth jet containers by radius key (SIM only)
  // -------------------------------------------------------------------------
  std::map<std::string, JetContainer*> m_truthJetsByRKey;
  std::map<std::string, std::string>   m_truthJetsNodeByRKey;

  // -------------------------------------------------------------------------
  // EventDisplay diagnostics payload support nodes (optional; never affects physics)
  // -------------------------------------------------------------------------
  EventHeader*           m_evtHeader          = nullptr;
  TowerInfoContainer*    m_evtDispTowersCEMC  = nullptr;
  TowerInfoContainer*    m_evtDispTowersIHCal = nullptr;
  TowerInfoContainer*    m_evtDispTowersOHCal = nullptr;
  RawTowerGeomContainer* m_evtDispGeomCEMC    = nullptr;
  RawTowerGeomContainer* m_evtDispGeomIHCal   = nullptr;
  RawTowerGeomContainer* m_evtDispGeomOHCal   = nullptr;

  // -------------------------------------------------------------------------
  // EventDisplay diagnostics payload (offline rendering; independent of Verbosity()).
  //
  //  - One TTree entry per (event, rKey) when enabled.
  //  - Stores compact "ingredients" (jet kinematics + sparse tower constituent lists)
  //    so that single-event displays can be rendered OFFLINE without any online scanning
  //    or Verbosity() hacks.
  // -------------------------------------------------------------------------
  bool m_evtDiagEnabled    = false;  // user-controlled switch (macro setter or env)
  bool m_evtDiagNodesReady = false;  // per-event: true if required tower/geom nodes are available

  int  m_evtDiagMaxPerBin = 0;       // 0 => unlimited; else limit entries per (rKey,ptBin,cat)
  std::unordered_map<std::string, int> m_evtDiagSavedPerBin;

  long long m_evtDiagNFill = 0;
  long long m_evtDiagNFillWithSelTowers = 0;
  long long m_evtDiagNFillWithBestTowers = 0;
  long long m_evtDiagNFillWithAnyTowers = 0;
  long long m_evtDiagNFillByCat[3] = {0, 0, 0};
  long long m_evtDiagNFillWithAnyTowersByCat[3] = {0, 0, 0};

  TTree* m_evtDiagTree = nullptr;

  int       m_evtDiag_run        = 0;
  int       m_evtDiag_evt        = 0;
  long long m_evtDiag_eventCount = 0;
  float     m_evtDiag_vz         = 0.0f;

    std::string m_evtDiag_rKey;
    int         m_evtDiag_ptBin = -1;
    int         m_evtDiag_cat   = -1;
    int         m_evtDiag_isSim = 0;

    float m_evtDiag_ptGammaTruth  = 0.0f;
    float m_evtDiag_phiGammaTruth = 0.0f;
    float m_evtDiag_ptGammaReco   = 0.0f;
    float m_evtDiag_etaGammaReco  = 0.0f;
    float m_evtDiag_phiGammaReco  = 0.0f;

  float m_evtDiag_sel_pt  = 0.0f;
  float m_evtDiag_sel_eta = 0.0f;
  float m_evtDiag_sel_phi = 0.0f;

  float m_evtDiag_best_pt  = 0.0f;
  float m_evtDiag_best_eta = 0.0f;
  float m_evtDiag_best_phi = 0.0f;

  float m_evtDiag_truthLead_pt  = 0.0f;
  float m_evtDiag_truthLead_eta = 0.0f;
  float m_evtDiag_truthLead_phi = 0.0f;

  float m_evtDiag_drSelToTruthLead  = -1.0f;
  float m_evtDiag_drBestToTruthLead = -1.0f;

  std::vector<int>   m_evtDiag_sel_calo;
  std::vector<int>   m_evtDiag_sel_ieta;
  std::vector<int>   m_evtDiag_sel_iphi;
  std::vector<float> m_evtDiag_sel_etaTower;
  std::vector<float> m_evtDiag_sel_phiTower;
  std::vector<float> m_evtDiag_sel_etTower;
  std::vector<float> m_evtDiag_sel_eTower;

  std::vector<int>   m_evtDiag_best_calo;
  std::vector<int>   m_evtDiag_best_ieta;
  std::vector<int>   m_evtDiag_best_iphi;
  std::vector<float> m_evtDiag_best_etaTower;
  std::vector<float> m_evtDiag_best_phiTower;
  std::vector<float> m_evtDiag_best_etTower;
  std::vector<float> m_evtDiag_best_eTower;

  // -------------------------------------------------------------------------
  // Diagnostics / accounting
  // -------------------------------------------------------------------------
  EventReject m_lastReject = EventReject::None;

  long long event_count  = 0;
  long long m_evtNoTrig  = 0;

  Bookkeeping m_bk{};

  // For End() histogram summary (counts how many times each histogram was filled)
  std::unordered_map<std::string, long long> m_histFill;

  // -------------------------------------------------------------------------
  // NEW: per-pT-bin negative-isolation bookkeeping (independent of SS classification)
  //
  //  - "Builder"    = RecoilJets::eiso() using PhotonClusterBuilder iso_* layer sums
  //  - "ClusterIso" = RawCluster::get_et_iso(radiusx10,false,true) (UNSUBTRACTED)
  // -------------------------------------------------------------------------
  std::vector<unsigned long long> m_nIsoBuilderByPt;
  std::vector<unsigned long long> m_nIsoBuilderNegByPt;
  std::vector<unsigned long long> m_nIsoClusterIsoByPt;
  std::vector<unsigned long long> m_nIsoClusterIsoNegByPt;

  // Histogram storage: trigger -> (name -> object*)
  std::map<std::string, HistMap> qaHistogramsByTrigger;

  // Per-trigger slice counters printed in End()
  std::map<std::string, std::map<std::string, CatStat>> m_catByTrig;
};

#endif // RECOILJETS_H
