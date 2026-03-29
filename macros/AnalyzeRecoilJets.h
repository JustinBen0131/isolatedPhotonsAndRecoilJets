// AnalyzeRecoilJets.h

#ifndef ANALYZE_RECOIL_JETS_H
#define ANALYZE_RECOIL_JETS_H

// =============================================================================
// ROOT
// =============================================================================
#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TIterator.h>
#include <TObject.h>
#include <TNamed.h>
#include <TAxis.h>

#include <TH1.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2.h>
#include <TH2F.h>
#include <TH3.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>

#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMarker.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLine.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TMatrixDSym.h>
#include <TMultiGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TString.h>

// =============================================================================
// C++ STL
// =============================================================================
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <utility>
#include <limits>
#include <random>

namespace ARJ
{
  // ---------------------------------------------------------------------------
  // Convenience aliases (kept inside ARJ to avoid polluting global namespace)
  // ---------------------------------------------------------------------------
  using std::cout;
  using std::cerr;
  using std::endl;
  using std::string;
  using std::vector;
  using std::map;

  // #############################################################################
  // #                                                                           #
  // #                       CONFIGURATION PANEL                                 #
  // #                                                                           #
  // #  Edit ONLY this section. Everything below is derived automatically.       #
  // #                                                                           #
  // #############################################################################

  // ===========================================================================
  // 1. RUN MODE  (set exactly one group to true)
  // ===========================================================================
  //   PP DATA ONLY:   isPPdataOnly=true,  isSimAndDataPP=false
  //   SIM ONLY:       isPPdataOnly=false, isSimAndDataPP=false  + one SIM toggle
  //   SIM + PP DATA:  isPPdataOnly=false, isSimAndDataPP=true   + one SIM toggle
  //   AuAu ONLY:      isAuAuOnly=true    (all others false)
  // ---------------------------------------------------------------------------
  inline bool isPPdataOnly   = false;
  inline bool isSimAndDataPP = false;
  inline bool isAuAuOnly     = true;
  inline bool isPPdataAndAUAU = true;
  inline bool isRun25pp      = false;

  // ===========================================================================
  // 2. SIM SAMPLE  (set exactly ONE to true when SIM is included)
  // ===========================================================================
  //   Single-slice (raw counts):
  inline bool isPhotonJet5               = false;
  inline bool isPhotonJet10              = false;
  inline bool isPhotonJet20              = false;
  //   Weighted merged (y-axis ~ "Counts / pb^{-1}"):
  inline bool bothPhoton5and10sim        = false;
  inline bool bothPhoton5and20sim        = false;
  inline bool bothPhoton10and20sim       = false;
  inline bool allPhoton5and10and20sim    = false;
  //   Special SIM samples:
  inline bool isSimMB                    = false;   // MinBias DETROIT tune
  inline bool isSimJet5                  = false;   // inclusive jet5
  inline bool isSimEmbedded              = false;   // embedded photon20 in AuAu

  inline bool doPhotonJetMerge = false;

  //   RooUnfold: true = run both non-purity and purity-corrected passes + overlay.
  inline bool do_xJ_PPunfold = false;
  //   Internal: selects raw vs ABCD purity-corrected reco inputs per pass.
  inline bool gApplyPurityCorrectionForUnfolding = false;

  //   One-off Sam vs Justin unsmear comparison:
  inline bool doSamVsJustinUnsmearOverlays = false;

  // ===========================================================================
  // 3a. PP / SIM CUT DEFAULTS  (drives all PP, SIM, and SIM+PP paths)
  //     All input/output paths are derived automatically from these values.
  // ===========================================================================
  inline const int    kJetPtMin        = 5;            // GeV: 3, 5, or 10
  inline const string kB2BCut          = "7pi_8";      // "7pi_8" or "pi_2"
  inline const int    kVzCut           = 30;            // cm: 30 or 60
  inline const string kIsoConeR        = "isoR30";     // "isoR30" or "isoR40"
  inline const string kIsoMode         = "isSliding";  // "isSliding" or "fixedIso5GeV"
  inline const double kPhotonEtaAbsMax = 0.7;

  // ===========================================================================
  // 3b. Au+Au CUT DEFAULTS  (independent from PP/SIM — drives all AuAu and
  //     embedded-SIM paths; does NOT need to match the PP cuts above)
  // ===========================================================================
  inline const int    kAA_JetPtMin     = 5;            // GeV: 3, 5, or 10
  inline const string kAA_B2BCut       = "7pi_8";      // "7pi_8" or "pi_2"
  inline const int    kAA_VzCut        = 30;            // cm: 30 or 60
  inline const string kAA_IsoConeR     = "isoR40";     // "isoR30" or "isoR40"
  inline const string kAA_IsoMode      = "fixedIso5GeV";// "isSliding" or "fixedIso5GeV"
  inline const string kAA_UEVariant    = "variantB";      // "noSub","baseVariant","variantA","variantB"
  // Au+Au trigger directory name(s) inside the ROOT file.
  // Set one, two, or all three.  Analysis runs independently for each.
  inline const vector<string> kTriggersAuAu = {
  //      "MBD_NS_geq_2_vtx_lt_150"
       "photon_10_plus_MBD_NS_geq_2_vtx_lt_150"
      // ,"photon_12_plus_MBD_NS_geq_2_vtx_lt_150"
  };

  // ===========================================================================
  // 4. BINNING  (edit these arrays to change pT / xJ / centrality slicing)
  // ===========================================================================
  //   JES3 photon pT bin edges (drives pT-sliced QA tables + JES3 booking):
  inline const vector<double> kJES3PhotonPtBins         = {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 35};
  //   Unfolding photon pT bin edges (reco and truth, independent from JES3):
  inline const vector<double> kUnfoldRecoPhotonPtBins   = {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 35, 40};
  inline const vector<double> kUnfoldTruthPhotonPtBins  = {5, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 35, 40};
  //   xJ = pT(jet) / pT(gamma) bin edges:
  inline const vector<double> kUnfoldXjBins             = {0.0, 0.20, 0.24, 0.29, 0.35, 0.41, 0.50, 0.60, 0.72, 0.86, 1.03, 1.24, 1.49, 1.78, 2.14, 3.0};
  //   Uniform jet pT binning (expanded into edges at startup):
  inline constexpr double kUnfoldJetPtStart = 0.0;
  inline constexpr double kUnfoldJetPtStop  = 60.0;
  inline constexpr double kUnfoldJetPtStep  = 0.5;
  //   Au+Au centrality bin edges (percent):
  inline const vector<double> kCentralityEdges          = {0, 10, 20, 40, 60, 80};


  // ===========================================================================
  // 6. CROSS SECTIONS  (pb, for SIM slice weighting)
  // ===========================================================================
  inline constexpr double kSigmaPhoton5_pb  = 89266.571;
  inline constexpr double kSigmaPhoton10_pb = 6692.7611;
  inline constexpr double kSigmaPhoton20_pb = 105.79868;

  // #############################################################################
  // #                  END OF CONFIGURATION PANEL                               #
  // #  Everything below is derived — you should not need to edit it.            #
  // #############################################################################

  // True if the selected SIM sample is a weighted multi-slice merge (hist units become ~pb/bin)
  inline bool IsWeightedSIMSelected()
  {
        return (bothPhoton5and10sim || bothPhoton5and20sim || bothPhoton10and20sim || allPhoton5and10and20sim);
  }

  // --- Derived tag builders (DO NOT EDIT) ---
  // PP / SIM tag (Section 3a defaults)
  inline string CfgTag()
  {
        return "jetMinPt" + std::to_string(kJetPtMin) + "_" +
               kB2BCut + "_vz" + std::to_string(kVzCut) + "_" +
               kIsoConeR + "_" + kIsoMode;
  }

  // Au+Au tag (Section 3b defaults)
  inline string CfgTagAA()
  {
        return "jetMinPt" + std::to_string(kAA_JetPtMin) + "_" +
               kAA_B2BCut + "_vz" + std::to_string(kAA_VzCut) + "_" +
               kAA_IsoConeR + "_" + kAA_IsoMode;
  }

  inline string CfgTagWithUE_AA()
  {
        return CfgTagAA() + "_" + kAA_UEVariant;
  }

  // Convenience alias — all existing CfgTagWithUE() call sites are AuAu-context
  inline string CfgTagWithUE()
  {
        return CfgTagWithUE_AA();
  }

  // For overlay comparisons across different cut combos
  inline string CfgTagFor(int jetPtMin,
                          const string& b2bCut,
                          int vzCut,
                          const string& isoConeR,
                          const string& isoMode)
  {
      return "jetMinPt" + std::to_string(jetPtMin) + "_" +
             b2bCut + "_vz" + std::to_string(vzCut) + "_" +
             isoConeR + "_" + isoMode;
  }

  inline string CfgTagWithUEFor(int jetPtMin,
                                const string& b2bCut,
                                int vzCut,
                                const string& isoConeR,
                                const string& isoMode,
                                const string& ueVariant)
  {
      return CfgTagFor(jetPtMin, b2bCut, vzCut, isoConeR, isoMode) + "_" + ueVariant;
  }

  inline string B2BLabelFor(const string& b2bCut)
  {
      if (b2bCut == "7pi_8") return "7#pi/8";
      if (b2bCut == "pi_2")  return "#pi/2";
      return b2bCut;
  }

  // Back-to-back label for plot annotations (ROOT TLatex format)
  inline string B2BLabel()
  {
      if (kB2BCut == "7pi_8") return "7#pi/8";
      if (kB2BCut == "pi_2")  return "#pi/2";
      return kB2BCut;
  }

  // vzCutCm as a double for existing plot labels that use it
  inline const double vzCutCm = static_cast<double>(kVzCut);

  // =============================================================================
  // FIXED INPUTS / TRIGGERS
  // =============================================================================
  inline const string kTriggerPP       = "Photon_4_GeV_plus_MBD_NS_geq_1";
  // Backward-compatible alias: defaults to first entry in kTriggersAuAu
  inline const string kTriggerAuAuGold = kTriggersAuAu.front();
  inline const string kDirSIM          = "SIM";

  // =========================================================================
  // INPUT PATH BUILDERS — PP/SIM from CfgTag(), AuAu from CfgTagWithUE_AA()
  // =========================================================================
  inline const string kThesisBase = "/Users/patsfan753/Desktop/ThesisAnalysis";
  inline const string kInputBase  = kThesisBase + "/InputFiles";
  inline const string kOutputBase = kThesisBase + "/dataOutput";

  // --- Input paths (one function per dataset) ---
  inline string InputPP()
  {
      return kInputBase + "/pp24/RecoilJets_pp_ALL_" + CfgTag() + ".root";
  }
  inline string InputPP25()
  {
      return kInputBase + "/pp25/RecoilJets_pp25_ALL_" + CfgTag() + ".root";
  }
  inline string InputPP(bool run25)
  {
      return run25 ? InputPP25() : InputPP();
  }
  inline string InputAuAu()
  {
      return kInputBase + "/auau25/RecoilJets_auau_ALL_" + CfgTagWithUE() + ".root";
  }
  inline string InputSim(const string& sampleTag)
  {
      // sampleTag = "photonjet5", "photonjet10", "photonjet20"
      return kInputBase + "/simPhotonJet/RecoilJets_" + sampleTag + "_ALL_" + CfgTag() + ".root";
  }
  inline string InputSim(const string& sampleTag, const string& cfgTag)
  {
      return kInputBase + "/simPhotonJet/RecoilJets_" + sampleTag + "_ALL_" + cfgTag + ".root";
  }
  inline string InputSimJet5()
  {
      return kInputBase + "/InclusiveJetSIM/RecoilJets_jet5_ALL_" + CfgTag() + ".root";
  }
  inline string InputSimMB()
  {
      return kInputBase + "/MinBiasSIM_DETROITtune/RecoilJets_detroit_ALL_" + CfgTag() + ".root";
  }
  inline string InputSimEmbedded()
  {
      return kInputBase + "/simEmbedded/RecoilJets_embeddedPhoton20_ALL_" + CfgTagWithUE() + ".root";
  }
  // Variant override for embedded (e.g. loop over UE variants)
  inline string InputSimEmbedded(const string& ueVariant)
  {
    return kInputBase + "/simEmbedded/RecoilJets_embeddedPhoton20_ALL_" +
           CfgTagAA() + "_" + ueVariant + ".root";
  }
  // =========================================================================
  // OUTPUT PATH BUILDERS — tag-aware, never overwrite across cut combos
  // =========================================================================
  inline string OutputPP()
  {
      return kOutputBase + "/pp/" + CfgTag();
  }
  inline string OutputAuAu()
  {
      return kOutputBase + "/auau/" + CfgTagWithUE();
  }
  inline string OutputPPAuAu()
  {
      return kOutputBase + "/pp_auau/pp_" + CfgTag() + "__aa_" + CfgTagWithUE_AA();
  }
  inline string OutputIndividualSim(const string& simLabel)
  {
      // simLabel = "photonJet5_SIM", "photonJet10_SIM", "photonJet20_SIM"
      return kOutputBase + "/individualSim/" + CfgTag() + "/" + simLabel;
  }
  inline string OutputCombinedSimOnly(const string& comboLabel)
  {
      // comboLabel = "photonJet10and20merged_SIM", etc.
      return kOutputBase + "/combinedSimOnly/" + CfgTag() + "/" + comboLabel;
  }
  inline string OutputCombinedSimOnly(const string& cfgTag, const string& comboLabel)
  {
      return kOutputBase + "/combinedSimOnly/" + cfgTag + "/" + comboLabel;
  }
  inline string OutputSimAndData(const string& comboLabel)
  {
      // comboLabel = "photonJet10and20merged_SIM", etc.
      return kOutputBase + "/simAndData/" + CfgTag() + "/" + comboLabel;
  }
  inline string OutputSimMB()
  {
      return kOutputBase + "/simMBpp/" + CfgTag();
  }
  inline string OutputSimJet5()
  {
      return kOutputBase + "/simJet5ppINCLUSIVE/" + CfgTag();
  }
  inline string OutputSimEmbedded()
  {
      return kOutputBase + "/auau/embeddedPhoton20/" + CfgTagWithUE();
  }

  // Merged SIM ROOT file path (lives inside the combinedSimOnly output dir)
  inline string MergedSimPath(const string& comboLabel, const string& mergedFilename)
  {
      // e.g. MergedSimPath("photonJet10and20merged_SIM", "RecoilJets_photonjet10plus20_MERGED.root")
      return OutputCombinedSimOnly(comboLabel) + "/" + mergedFilename;
  }
  inline string MergedSimPath(const string& cfgTag, const string& comboLabel, const string& mergedFilename)
  {
      return OutputCombinedSimOnly(cfgTag, comboLabel) + "/" + mergedFilename;
  }

  // =============================================================================
  // Binning (initialized from Configuration Panel constants above)
  // =============================================================================
  struct BinningCfg
  {
        vector<double> jes3_photon_pt_bins         = {8, 10,12,14,16,18,20,22,24,26,35};
        vector<double> unfold_reco_photon_pt_bins  = {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 35, 40};
        vector<double> unfold_truth_photon_pt_bins = {5, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 35, 40};
        vector<double> unfold_xj_bins              = {0.0,0.20,0.24,0.29,0.35,0.41,0.50,0.60,0.72,0.86,1.03,1.24,1.49,1.78,2.14,3.0};

        // Au+Au centrality bin edges (percent) used for PP vs AuAu overlays
        vector<double> centrality_edges             = {0, 10, 20, 40, 60, 80};

        double unfold_jet_pt_start = 0.0;
        double unfold_jet_pt_stop  = 60.0;
        double unfold_jet_pt_step  = 0.5;
        vector<double> unfold_jet_pt_edges;  // expanded from start/stop/step
    };

    inline string Trim(std::string s)
    {
      const char* ws = " \t\r\n";
      s.erase(0, s.find_first_not_of(ws));
      s.erase(s.find_last_not_of(ws) + 1);
      return s;
    }

    inline void ExpandUniformEdges(vector<double>& edges, double start, double stop, double step)
    {
      edges.clear();
      if (!(std::isfinite(start) && std::isfinite(stop) && std::isfinite(step))) return;
      if (step <= 0.0) return;
      if (stop <= start) return;

      const int n = (int) std::llround((stop - start) / step);
      edges.reserve((std::size_t)n + 2);
      for (int i = 0; i <= n; ++i)
      {
        edges.push_back(start + step * (double)i);
      }
      if (edges.empty() || std::fabs(edges.back() - stop) > 1e-9) edges.push_back(stop);
    }

    inline const BinningCfg& Binning()
    {
        static BinningCfg cfg;
        static bool loaded = false;
        if (loaded) return cfg;
        loaded = true;

        ExpandUniformEdges(cfg.unfold_jet_pt_edges, cfg.unfold_jet_pt_start, cfg.unfold_jet_pt_stop, cfg.unfold_jet_pt_step);
        if (cfg.unfold_jet_pt_edges.size() < 2)
          ExpandUniformEdges(cfg.unfold_jet_pt_edges, 0.0, 60.0, 0.5);

        return cfg;
    }

  // JES3 photon pT bin edges (authoritative for pT-bin loops/labels)
  inline const vector<double>& kPtEdges = Binning().jes3_photon_pt_bins;
  inline const int kNPtBins = (int)kPtEdges.size() - 1;

  // Unfolding photon pT bin edges (explicit, independent from JES3 bins)
  inline const vector<double>& kUnfoldRecoPtEdges  = Binning().unfold_reco_photon_pt_bins;
  inline const vector<double>& kUnfoldTruthPtEdges = Binning().unfold_truth_photon_pt_bins;
  inline const int kNUnfoldRecoPtBins  = (int)kUnfoldRecoPtEdges.size()  - 1;
  inline const int kNUnfoldTruthPtBins = (int)kUnfoldTruthPtEdges.size() - 1;

    // Jet radii keys
    inline const vector<string> kRKeys = {"r02","r04","r06"};

  // =============================================================================
  // ANSI helpers
  // =============================================================================
  inline const string ANSI_RESET    = "\033[0m";
  inline const string ANSI_BOLD_RED = "\033[1;31m";
  inline const string ANSI_BOLD_GRN = "\033[1;32m";
  inline const string ANSI_BOLD_YEL = "\033[1;33m";
  inline const string ANSI_BOLD_CYN = "\033[1;36m";
  inline const string ANSI_BOLD_WHT = "\033[1;37m";
  inline const string ANSI_DIM      = "\033[2m";


  // =============================================================================
  // Small data structures
  // =============================================================================
  struct PtBin
  {
      int lo = 0;
      int hi = 0;
      string label;   // "10-12"
      string folder;  // "pT_10_12"
      string suffix;  // "_pT_10_12"
  };

    // Cached bin list
    inline const vector<PtBin>& PtBins()
    {
        static vector<PtBin> v;
        if (!v.empty()) return v;
        v.reserve((std::size_t)kNPtBins);

        for (int i = 0; i < kNPtBins; ++i)
        {
          PtBin b;
          b.lo = (int) std::llround(kPtEdges[(std::size_t)i]);
          b.hi = (int) std::llround(kPtEdges[(std::size_t)i+1]);
          {
            std::ostringstream s; s << b.lo << "-" << b.hi; b.label = s.str();
          }
          {
            std::ostringstream s; s << "pT_" << b.lo << "_" << b.hi; b.folder = s.str();
          }
          {
            std::ostringstream s; s << "_pT_" << b.lo << "_" << b.hi; b.suffix = s.str();
          }
          v.push_back(b);
        }
        return v;
    }

    // Cached unfolding RECO pT bin list (derived from BinningCfg::unfold_reco_photon_pt_bins)
    inline const vector<PtBin>& UnfoldRecoPtBins()
    {
        static vector<PtBin> v;
        if (!v.empty()) return v;

        const auto& edges = kUnfoldRecoPtEdges;
        const int n = (int)edges.size() - 1;
        if (n <= 0) return v;

        v.reserve((std::size_t)n);

        for (int i = 0; i < n; ++i)
        {
          PtBin b;
          b.lo = (int) std::llround(edges[(std::size_t)i]);
          b.hi = (int) std::llround(edges[(std::size_t)i+1]);
          {
            std::ostringstream s; s << b.lo << "-" << b.hi; b.label = s.str();
          }
          {
            std::ostringstream s; s << "pT_" << b.lo << "_" << b.hi; b.folder = s.str();
          }
          {
            std::ostringstream s; s << "_pT_" << b.lo << "_" << b.hi; b.suffix = s.str();
          }
          v.push_back(b);
        }

        return v;
    }

    // Analysis-visible unfolding RECO pT bins:
    //   - exclude the first reco bin  (8-10)  -> treat as reco underflow support bin
    //   - exclude the last reco bin   (35-40) -> treat as reco overflow support bin
    // This yields the 9 analysis bins:
    //   10-12, 12-14, 14-16, 16-18, 18-20, 20-22, 22-24, 24-26, 26-35
    inline const vector<PtBin>& UnfoldAnalysisRecoPtBins()
    {
        static vector<PtBin> v;
        if (!v.empty()) return v;

        const auto& all = UnfoldRecoPtBins();
        if (all.size() <= 2) return v;

        v.reserve(all.size() - 2);
        for (std::size_t i = 1; i + 1 < all.size(); ++i)
        {
          v.push_back(all[i]);
        }

        return v;
    }

    struct CentBin
    {
        int lo = 0;
        int hi = 0;
        string label;   // "0-10"
        string folder;  // "0_10"
        string suffix;  // "_cent_0_10"
    };

    // Cached centrality bin list (derived from YAML key: centrality_edges)
    inline const vector<CentBin>& CentBins()
    {
        static vector<CentBin> v;
        if (!v.empty()) return v;

        const auto& edges = Binning().centrality_edges;
        if (edges.size() < 2) return v;

        v.reserve(edges.size() - 1);

        for (std::size_t i = 0; i + 1 < edges.size(); ++i)
        {
          CentBin b;
          b.lo = (int) std::llround(edges[i]);
          b.hi = (int) std::llround(edges[i+1]);
          {
            std::ostringstream s; s << b.lo << "-" << b.hi; b.label = s.str();
          }
          {
            std::ostringstream s; s << b.lo << "_" << b.hi; b.folder = s.str();
          }
          {
            std::ostringstream s; s << "_cent_" << b.lo << "_" << b.hi; b.suffix = s.str();
          }
          v.push_back(b);
        }

        return v;
    }

    inline double RFromKey(const string& rKey)
    {
        if (rKey == "r02") return 0.2;
        if (rKey == "r04") return 0.4;
        if (rKey == "r06") return 0.6;
        return 0.0;
    }

  inline double FidEtaAbsFromKey(const string& rKey)
  {
    const double R = RFromKey(rKey);
    return 1.1 - R;
  }

  // Non-copyable but movable (safe for std::vector)
  struct Dataset
  {
        string label;       // "SIM" or "DATA"
        bool   isSim = false;

        string trigger;     // for data
        string topDirName;  // "SIM" or trigger

        string inFilePath;
        string outBase;

        string centFolder;  // e.g. "0_10" for Au+Au centrality-scoped plotting
        string centSuffix;  // e.g. "_cent_0_10"
        string centLabel;   // e.g. "Centrality: 0-10%"

        TFile*      file   = nullptr;
        TDirectory* topDir = nullptr;

        // Existing raw missing-object log (one line per failure occurrence)
        std::ofstream missingOut;
        int missingCount = 0;

        // NEW: Tracking for end-of-job histogram coverage summary
        // Key is ALWAYS the full path as printed in inventory:  "<topDirName>/<relName>"
        map<string, int>    requestCounts;   // fullpath -> #times GetObj was called
        map<string, int>    missingCounts;   // fullpath -> #times LogMissing fired
        map<string, string> missingReason;   // fullpath -> last recorded reason string

        Dataset() = default;
        Dataset(const Dataset&) = delete;
        Dataset& operator=(const Dataset&) = delete;
        Dataset(Dataset&&) noexcept = default;
        Dataset& operator=(Dataset&&) noexcept = default;
  };

  struct MatchCache
  {
    map<string, vector<double> > NphoLead;    // per rKey, size kNPtBins
    map<string, vector<double> > NphoMatched; // per rKey, size kNPtBins
  };

  struct LeakageFactors
  {
    vector<double> fB, fC, fD; // per pT bin
    bool available = false;
  };

  // Inventory item
  struct InvItem
  {
    string path;
    string cls;
    double entries;
  };

  // =============================================================================
  // Path helpers
  // =============================================================================
  inline void EnsureDir(const string& path)
  {
    if (path.empty()) return;
    gSystem->mkdir(path.c_str(), true);
  }

  inline string DirnameFromPath(const string& filepath)
  {
    const size_t pos = filepath.find_last_of('/');
    if (pos == string::npos) return "";
    return filepath.substr(0, pos);
  }

  inline void EnsureParentDirForFile(const string& filepath)
  {
    EnsureDir(DirnameFromPath(filepath));
  }

  inline string JoinPath(const string& a, const string& b)
  {
    if (a.empty()) return b;
    if (b.empty()) return a;
    if (a.back() == '/') return a + b;
    return a + "/" + b;
  }

  inline string FullPath(const Dataset& ds, const string& rel)
  {
    return ds.topDirName + "/" + rel;
  }

  inline void LogMissing(Dataset& ds, const string& fullpath, const string& reason)
  {
      if (ds.missingOut.is_open())
      {
        ds.missingOut << fullpath << " :: " << reason << "\n";
      }
      ds.missingCount++;

      // NEW: tabulate unique missing objects (and how often they were missing)
      ds.missingCounts[fullpath]++;
      ds.missingReason[fullpath] = reason;
  }

  // =============================================================================
  // Required-file opener (abort with clear message if missing)
  // =============================================================================
  inline TFile* OpenRequiredFile(const string& path, const string& context = "")
  {
      TFile* f = TFile::Open(path.c_str(), "READ");
      if (!f || f->IsZombie())
      {
          cerr << ANSI_BOLD_RED
               << "[FATAL] Required input file not found"
               << (context.empty() ? "" : " (" + context + ")")
               << ":\n  " << path << "\n"
               << ANSI_RESET;
          std::exit(1);
      }
      return f;
  }

  // =============================================================================
  // Robust ROOT getters (template must live in header)
  // =============================================================================
  template <class T>
  inline T* GetObj(Dataset& ds, const string& relName,
                     bool logMissing = true,
                     bool logZero = true,
                     bool treatZeroAsMissing = false)
    {
      if (!ds.topDir) return nullptr;

      vector<string> relNamesToTry;
      relNamesToTry.reserve(2);

      if (!ds.centSuffix.empty() && relName.find("_cent_") == string::npos)
      {
        relNamesToTry.push_back(relName + ds.centSuffix);
      }
      relNamesToTry.push_back(relName);

      TObject* obj = nullptr;
      string usedRelName;

      for (const auto& cand : relNamesToTry)
      {
        obj = ds.topDir->Get(cand.c_str());
        if (obj)
        {
          usedRelName = cand;
          break;
        }
      }

      if (!obj)
      {
        const string fp = FullPath(ds, relNamesToTry.front());
        ds.requestCounts[fp]++;
        if (logMissing) LogMissing(ds, fp, "MISSING");
        return nullptr;
      }

      const string fp = FullPath(ds, usedRelName);
      ds.requestCounts[fp]++;

      T* out = dynamic_cast<T*>(obj);
      if (!out)
      {
        if (logMissing) LogMissing(ds, fp, string("WRONG_TYPE actual=") + obj->ClassName());
        return nullptr;
      }

      // entries check (works for TH* and profiles; for others, skip)
      double ent = -1.0;
      if (auto* h = dynamic_cast<TH1*>(out)) ent = h->GetEntries();
      else if (auto* p = dynamic_cast<TProfile*>(out)) ent = p->GetEntries();
      else if (auto* p2 = dynamic_cast<TProfile2D*>(out)) ent = p2->GetEntries();
      else if (auto* p3 = dynamic_cast<TProfile3D*>(out)) ent = p3->GetEntries();

      if (ent >= 0.0 && ent <= 0.0)
      {
        if (logZero) LogMissing(ds, fp, "ZERO_ENTRIES");
        if (treatZeroAsMissing) return nullptr;
      }

      return out;
  }

  // =============================================================================
  // Small numeric + cloning helpers
  // =============================================================================
  inline double SafeDivide(double a, double b, double def = 0.0)
  {
    if (b == 0.0) return def;
    return a / b;
  }

  inline TH1* CloneTH1(const TH1* h, const string& newName)
  {
    if (!h) return nullptr;
    TH1* c = dynamic_cast<TH1*>(h->Clone(newName.c_str()));
    if (!c) return nullptr;
    c->SetDirectory(nullptr);
    return c;
  }

  inline TH2* CloneTH2(const TH2* h, const string& newName)
  {
    if (!h) return nullptr;
    TH2* c = dynamic_cast<TH2*>(h->Clone(newName.c_str()));
    if (!c) return nullptr;
    c->SetDirectory(nullptr);
    return c;
  }

  inline TH3* CloneTH3(const TH3* h, const string& newName)
  {
    if (!h) return nullptr;
    TH3* c = dynamic_cast<TH3*>(h->Clone(newName.c_str()));
    if (!c) return nullptr;
    c->SetDirectory(nullptr);
    return c;
  }

  inline void EnsureSumw2(TH1* h)
  {
    if (!h) return;
    if (h->GetSumw2N() == 0) h->Sumw2();
  }

  inline double SmallestPositiveBinContent(const TH1* h)
  {
    if (!h) return 0.0;
    double minPos = std::numeric_limits<double>::max();
    const int nb = h->GetNbinsX();
    for (int i = 1; i <= nb; ++i)
    {
      const double v = h->GetBinContent(i);
      if (v > 0.0 && v < minPos) minPos = v;
    }
    if (!std::isfinite(minPos) || minPos == std::numeric_limits<double>::max()) return 0.0;
    return minPos;
  }

  inline void NormalizeToUnitArea(TH1* h)
  {
    if (!h) return;
    const double integral = h->Integral(0, h->GetNbinsX() + 1);
    if (integral > 0.0) h->Scale(1.0 / integral);
  }

  inline void NormalizeToUnitAreaInRange(TH1* h, double xmin, double xmax)
  {
    if (!h) return;
    const int b1 = h->GetXaxis()->FindBin(xmin + 1e-9);
    const int b2 = h->GetXaxis()->FindBin(xmax - 1e-9);
    const double integral = h->Integral(b1, b2);
    if (integral > 0.0) h->Scale(1.0 / integral);
  }

  // Build an "inclusive over pTgamma bins" reco histogram by summing the pT-sliced family.
  inline TH1* SumRecoPtSlicedTH1(Dataset& ds, const string& histBase)
  {
    TH1* sum = nullptr;

    for (int i = 0; i < kNPtBins; ++i)
    {
      const PtBin& b = PtBins()[i];
      const string hname = histBase + b.suffix;

      TH1* h = GetObj<TH1>(ds, hname, false, false, true);
      if (!h) continue;

      if (!sum)
      {
        sum = CloneTH1(h, histBase + "_INCLUSIVE_overPt");
        if (sum)
        {
          sum->Reset("ICES");
          sum->SetDirectory(nullptr);
        }
      }
      if (sum) sum->Add(h);
    }

    return sum;
  }

  // =============================================================================
  // Style / plotting primitives
  // =============================================================================
  inline void SetupGlobalStyle()
  {
      gStyle->SetOptStat(0);
      gStyle->SetTitleBorderSize(0);
      gStyle->SetTitleFillColor(0);
      gStyle->SetPadLeftMargin(0.16);
      gStyle->SetPadRightMargin(0.05);
      gStyle->SetPadBottomMargin(0.14);
      gStyle->SetPadTopMargin(0.08);
      gStyle->SetTitleXOffset(0.95);
      gStyle->SetTitleYOffset(1.20);
      gStyle->SetTitleOffset(0.95,"x");
      gStyle->SetTitleOffset(1.20,"y");
      gStyle->SetTitleOffset(1.15,"z");
      gStyle->SetLegendBorderSize(0);
      gStyle->SetLegendFillColor(0);
  }

  inline void ApplyCanvasMargins1D(TCanvas& c)
  {
    c.SetLeftMargin(0.16);
    c.SetRightMargin(0.05);
    c.SetBottomMargin(0.14);
    c.SetTopMargin(0.08);
  }

  inline void ApplyCanvasMargins2D(TCanvas& c)
  {
    c.SetLeftMargin(0.18);
    c.SetRightMargin(0.15);
    c.SetBottomMargin(0.16);
    c.SetTopMargin(0.08);
  }

  inline void DrawLatexLines(double x, double y,
                             const vector<string>& lines,
                             double size = 0.035,
                             double dy   = 0.045)
  {
    TLatex latex;
    latex.SetNDC(true);
    latex.SetTextFont(42);
    latex.SetTextSize(size);
    double yy = y;
    for (const auto& s : lines)
    {
      latex.DrawLatex(x, yy, s.c_str());
      yy -= dy;
    }
  }

  inline vector<string> DefaultHeaderLines(const Dataset& ds)
  {
          vector<string> lines;

          if (ds.isSim)
          {
            std::string simLabel = "UNKNOWN";

            if (isPhotonJet5)                 simLabel = "photonJet5";
            else if (isPhotonJet10)           simLabel = "photonJet10";
            else if (isPhotonJet20)           simLabel = "photonJet20";
            else if (bothPhoton5and10sim)     simLabel = "photonJet5+10 merged";
            else if (bothPhoton5and20sim)     simLabel = "photonJet5+20 merged";
            else if (bothPhoton10and20sim)    simLabel = "photonJet10+20 merged";
            else if (allPhoton5and10and20sim) simLabel = "photonJet5+10+20 merged";

            lines.push_back(std::string("Dataset: SIM (") + simLabel + ")");
          }
          else
          {
            lines.push_back(std::string("Dataset: DATA (") + ds.trigger + ")");
          }

          if (!ds.centLabel.empty())
          {
            lines.push_back(ds.centLabel);
          }

          return lines;
  }

  inline void SaveCanvas(TCanvas& c, const string& filepath)
  {
    EnsureParentDirForFile(filepath);
    c.SaveAs(filepath.c_str());
  }

  inline void DrawFidEtaLines1D(double etaAbs)
  {
    if (!gPad) return;
    const double ymin = gPad->GetUymin();
    const double ymax = gPad->GetUymax();
    TLine l1(+etaAbs, ymin, +etaAbs, ymax);
    TLine l2(-etaAbs, ymin, -etaAbs, ymax);
    l1.SetLineStyle(2); l2.SetLineStyle(2);
    l1.SetLineWidth(2); l2.SetLineWidth(2);
    l1.Draw("same");
    l2.Draw("same");
  }

  // =============================================================================
  // ISO (Eiso) 1D label helpers (for per-pT single PNGs)
  // =============================================================================
  inline string IsoTitleForBase(const string& base)
  {
      if (base == "h_Eiso_emcal")   return "E_{T}^{iso, EMCal}";
      if (base == "h_Eiso_hcalin")  return "E_{T}^{iso, IHCal}";
      if (base == "h_Eiso_hcalout") return "E_{T}^{iso, OHCal}";
      if (base == "h_Eiso")         return "E_{T}^{iso, Total}";
      return "E_{T}^{iso}";
  }

  inline void DrawIsoCornerLabels1D(const Dataset& ds, const string& isoTitle, const PtBin& pb)
  {
        TLatex t;
        t.SetNDC(true);
        t.SetTextFont(42);

        // Title (top-left) — match your table styling
        t.SetTextAlign(13);     // left, top
        t.SetTextSize(0.08);
        t.DrawLatex(0.45, 0.985, isoTitle.c_str());

        // pT bin (top-right, larger)
        t.SetTextAlign(33);     // right, top
        t.SetTextSize(0.070);
        t.DrawLatex(0.95, 0.88,
          TString::Format("p_{T}^{#gamma}: %d-%d GeV", pb.lo, pb.hi).Data());

        // eta cut (under pT, smaller)
        t.SetTextSize(0.055);
        t.DrawLatex(0.85, 0.75, "|#eta^{#gamma}| < 0.7");

        if (!ds.centLabel.empty())
        {
          t.SetTextSize(0.050);
          t.DrawLatex(0.95, 0.66, ds.centLabel.c_str());
        }
  }

  inline void DrawAndSaveTH1_Iso(const Dataset& ds,
                                   TH1* h,
                                   const string& filepath,
                                   const string& xTitle,
                                   const string& yTitle,
                                   const string& isoTitle,
                                   const PtBin& pb,
                                   bool logy,
                                   const char* drawOpt = "hist")
  {
      if (!h) return;

      EnsureSumw2(h);

      TCanvas c("c_iso","c_iso",900,700);
      ApplyCanvasMargins1D(c);
      c.SetLogy(logy);

      const bool mergedSimWeightedToPb = (ds.isSim && IsWeightedSIMSelected());

      const string yTitleEff =
        (yTitle == "A.U." || yTitle == "A.U")
          ? "Fraction of entries"
          : ((mergedSimWeightedToPb && (yTitle == "Counts" || yTitle == "counts"))
              ? "Counts / pb^{-1}"
              : yTitle);

      h->SetLineWidth(2);
      h->SetTitle("");
      h->GetXaxis()->SetTitle(xTitle.c_str());
      h->GetYaxis()->SetTitle(yTitleEff.c_str());

      const string opt = (drawOpt ? string(drawOpt) : string("hist"));

      if (logy)
      {
        const double minPos = SmallestPositiveBinContent(h);
        if (minPos > 0.0) h->SetMinimum(0.5 * minPos);
        else              h->SetMinimum(1e-6);
      }

      h->Draw(opt.c_str());

      // ISO corner label styling (NO DefaultHeaderLines, NO multi-line block)
      DrawIsoCornerLabels1D(ds, isoTitle, pb);

      SaveCanvas(c, filepath);
  }

  inline void DrawAndSaveTH1_Common(const Dataset& ds, TH1* h,
                                   const string& filepath,
                                   const string& xTitle,
                                   const string& yTitle,
                                   const vector<string>& extraLines,
                                   bool logy,
                                   bool drawFidLines = false,
                                   double fidEtaAbs  = 0.0,
                                   const char* drawOpt = "hist")
  {
    if (!h) return;

    EnsureSumw2(h);

    TCanvas c("c1","c1",900,700);
    ApplyCanvasMargins1D(c);
    c.SetLogy(logy);

      const bool mergedSimWeightedToPb = (ds.isSim && bothPhoton10and20sim);

      const string yTitleEff =
        (yTitle == "A.U." || yTitle == "A.U")
          ? "Fraction of entries"
          : ((mergedSimWeightedToPb && (yTitle == "Counts" || yTitle == "counts"))
              ? "Counts / pb^{-1}"
              : yTitle);

    h->SetLineWidth(2);
    h->SetTitle("");
    h->GetXaxis()->SetTitle(xTitle.c_str());
    h->GetYaxis()->SetTitle(yTitleEff.c_str());

    const string opt = (drawOpt ? string(drawOpt) : string("hist"));
    const bool wantsErrors =
      (opt.find("E") != string::npos) || (opt.find("e") != string::npos) ||
      (opt.find("P") != string::npos) || (opt.find("p") != string::npos);

    if (wantsErrors)
    {
      h->SetMarkerStyle(20);
      h->SetMarkerSize(1.0);
    }

    if (logy)
    {
      const double minPos = SmallestPositiveBinContent(h);
      if (minPos > 0.0) h->SetMinimum(0.5 * minPos);
      else              h->SetMinimum(1e-6);
    }

    h->Draw(opt.c_str());

    if (drawFidLines && fidEtaAbs > 0.0)
    {
      DrawFidEtaLines1D(fidEtaAbs);
    }

      vector<string> lines = DefaultHeaderLines(ds);
      for (const auto& s : extraLines) lines.push_back(s);
      DrawLatexLines(0.14, 0.92, lines, 0.034, 0.045);

      SaveCanvas(c, filepath);
    }

    // ---------------------------------------------------------------------------
    // RECO xJ (integrated over alpha): per-pTbin PNG with xJ floors + legend + title
    //
    // Draws:
    //   - blue dashed line: x_{J,min}^{abs}  = jetPtMin / pTmax^gamma
    //   - red  dashed line: x_{J,min}^{full} = jetPtMin / pTmin^gamma
    //   - legend with both relations + numeric value
    //   - big title: "RECO (unconditional) x_{J#gamma}, p_{T}^{#gamma} = X - Y GeV, R = 0.Z"
    // ---------------------------------------------------------------------------
    inline void DrawAndSave_xJRecoIntegratedAlpha_WithFloors(const Dataset& ds,
                                                             TH1* h,
                                                             const string& filepath,
                                                             double ptMinGamma,
                                                             double ptMaxGamma,
                                                             double R,
                                                             double jetPtMin_GeV = 10.0,
                                                             bool logy = false)
    {
      if (!h) return;

      EnsureSumw2(h);

      TCanvas c("c_xJReco","c_xJReco",900,700);
      ApplyCanvasMargins1D(c);
      c.SetLogy(logy);

      const bool mergedSimWeightedToPb = (ds.isSim && bothPhoton10and20sim);

      const string yTitleEff =
        (mergedSimWeightedToPb ? "Counts / pb^{-1}" : "Counts");

      h->SetTitle("");
      h->SetLineWidth(2);
      h->SetMarkerStyle(20);
      h->SetMarkerSize(1.0);

        h->GetXaxis()->SetTitle("x_{J#gamma}");
        h->GetXaxis()->SetRangeUser(0.0, 2.0);
        h->GetYaxis()->SetTitle(yTitleEff.c_str());

      if (logy)
      {
        const double minPos = SmallestPositiveBinContent(h);
        if (minPos > 0.0) h->SetMinimum(0.5 * minPos);
        else              h->SetMinimum(1e-6);
      }

        h->Draw("E1");
        gPad->Update();

        if (jetPtMin_GeV <= 0.0) jetPtMin_GeV = static_cast<double>(kJetPtMin);
        const string bbLabel = B2BLabel();

        const double xAbs  = (ptMaxGamma > 0.0) ? (jetPtMin_GeV / ptMaxGamma) : -1.0;
        const double xFull = (ptMinGamma > 0.0) ? (jetPtMin_GeV / ptMinGamma) : -1.0;

        const double yMin = gPad->GetUymin();
        const double yMax = gPad->GetUymax();

        TLine* lnAbs = new TLine(xAbs,  yMin, xAbs,  yMax);
        lnAbs->SetLineColor(kBlue + 1);
        lnAbs->SetLineStyle(2);
        lnAbs->SetLineWidth(2);

        TLine* lnFull = new TLine(xFull, yMin, xFull, yMax);
        lnFull->SetLineColor(kRed + 1);
        lnFull->SetLineStyle(2);
        lnFull->SetLineWidth(2);

        if (xAbs > 0.0)  lnAbs->Draw("same");
        if (xFull > 0.0) lnFull->Draw("same");

        TLegend* leg = new TLegend(0.52, 0.70, 0.95, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.032);

        if (xAbs > 0.0)
          leg->AddEntry(lnAbs,
            TString::Format("x_{J, min}^{abs} = #frac{%.0f}{p_{T, max}^{#gamma}} = %.3f", jetPtMin_GeV, xAbs),
            "l");
        if (xFull > 0.0)
          leg->AddEntry(lnFull,
            TString::Format("x_{J, min}^{full} = #frac{%.0f}{p_{T, min}^{#gamma}} = %.3f", jetPtMin_GeV, xFull),
            "l");

        leg->Draw();

        {
          TLatex tCuts;
          tCuts.SetNDC(true);
          tCuts.SetTextFont(42);
          tCuts.SetTextAlign(33);
          tCuts.SetTextSize(0.038);
          tCuts.DrawLatex(0.92, 0.62, TString::Format("|#Delta#phi(#gamma,jet)| > %s", bbLabel.c_str()).Data());
          tCuts.DrawLatex(0.92, 0.54, TString::Format("p_{T}^{jet} > %.0f GeV", jetPtMin_GeV).Data());
        }

        TLatex ttl;
        ttl.SetNDC(true);
        ttl.SetTextFont(42);
        ttl.SetTextSize(0.052);
        ttl.DrawLatex(0.12, 0.94,
          TString::Format("RECO (unconditional) x_{J#gamma}, p_{T}^{#gamma} = %.0f - %.0f GeV, R = %.1f",
            ptMinGamma, ptMaxGamma, R).Data());

      SaveCanvas(c, filepath);
      delete leg;
      delete lnFull;
      delete lnAbs;
  }

    


  inline void DrawAndSaveTH2_Common(const Dataset& ds, TH2* h,
                                   const string& filepath,
                                   const string& xTitle,
                                   const string& yTitle,
                                   const string& zTitle,
                                   const vector<string>& extraLines,
                                   bool logz,
                                   bool drawFidLines = false,
                                   double fidEtaAbs  = 0.0)
  {
    if (!h) return;

    TCanvas c("c2","c2",950,780);
    ApplyCanvasMargins2D(c);
    c.SetLogz(logz);

    h->SetTitle("");
    h->GetXaxis()->SetTitle(xTitle.c_str());
    h->GetYaxis()->SetTitle(yTitle.c_str());
    h->GetZaxis()->SetTitle(zTitle.c_str());

    h->Draw("COLZ");

    if (drawFidLines && fidEtaAbs > 0.0)
    {
      const double ymin = h->GetYaxis()->GetXmin();
      const double ymax = h->GetYaxis()->GetXmax();
      TLine l1(+fidEtaAbs, ymin, +fidEtaAbs, ymax);
      TLine l2(-fidEtaAbs, ymin, -fidEtaAbs, ymax);
      l1.SetLineStyle(2); l2.SetLineStyle(2);
      l1.SetLineWidth(2); l2.SetLineWidth(2);
      l1.Draw("same"); l2.Draw("same");
    }

    vector<string> lines = DefaultHeaderLines(ds);
    for (const auto& s : extraLines) lines.push_back(s);
    DrawLatexLines(0.14, 0.92, lines, 0.034, 0.045);

    SaveCanvas(c, filepath);
  }

  inline void DrawOverlayTwoTH1(const Dataset& ds,
                               TH1* h1, TH1* h2,
                               const string& label1, const string& label2,
                               const string& filepath,
                               const string& xTitle, const string& yTitle,
                               const vector<string>& extraLines,
                               bool logy)
  {
    if (!h1 || !h2) return;

    TCanvas c("c_ov","c_ov",900,700);
    ApplyCanvasMargins1D(c);
    c.SetLogy(logy);

    const string yTitleEff =
      (yTitle == "A.U." || yTitle == "A.U")
        ? "Fraction of entries"
        : yTitle;

    h1->SetTitle("");
    h1->GetXaxis()->SetTitle(xTitle.c_str());
    h1->GetYaxis()->SetTitle(yTitleEff.c_str());

    h1->SetLineWidth(2);
    h2->SetLineWidth(2);
    h1->SetLineColor(1);
    h2->SetLineColor(2);

    const double maxv = std::max(h1->GetMaximum(), h2->GetMaximum());
    h1->SetMaximum(maxv * 1.25);

    h1->Draw("hist");
    h2->Draw("hist same");

    TLegend leg(0.62, 0.73, 0.92, 0.88);
    leg.SetTextFont(42);
    leg.SetTextSize(0.033);
    leg.AddEntry(h1, label1.c_str(), "l");
    leg.AddEntry(h2, label2.c_str(), "l");
    leg.Draw();

    vector<string> lines = DefaultHeaderLines(ds);
    for (const auto& s : extraLines) lines.push_back(s);
    DrawLatexLines(0.14, 0.92, lines, 0.034, 0.045);

    SaveCanvas(c, filepath);
  }

  // =============================================================================
  // Vertex Z rebinner (0.5 cm bins over [-vzCut,+vzCut])
  // =============================================================================
  inline TH1F* RebinToFixedBinWidthVertexZ(const TH1* hOrig, double vzCut)
  {
    if (!hOrig) return nullptr;

    const double L = std::fabs(vzCut);
    const int nbVz = (int)std::lround(4.0 * L); // (2L)/0.5 = 4L
    if (nbVz <= 0) return nullptr;

    TH1F* hNew = new TH1F("h_vertexZ_fixed", "h_vertexZ_fixed", nbVz, -L, +L);
    hNew->Sumw2();

    const TAxis* ax = hOrig->GetXaxis();
    const int nbo = ax->GetNbins();

    for (int i = 1; i <= nbo; ++i)
    {
      const double xlo = ax->GetBinLowEdge(i);
      const double xhi = ax->GetBinUpEdge(i);
      const double w   = xhi - xlo;
      if (w <= 0) continue;

      const double content = hOrig->GetBinContent(i);
      const double err     = hOrig->GetBinError(i);

      const double olo = std::max(xlo, -L);
      const double ohi = std::min(xhi, +L);
      if (ohi <= olo) continue;

      const int jlo = hNew->GetXaxis()->FindBin(olo + 1e-9);
      const int jhi = hNew->GetXaxis()->FindBin(ohi - 1e-9);

      for (int j = jlo; j <= jhi; ++j)
      {
        const double nlo = hNew->GetXaxis()->GetBinLowEdge(j);
        const double nhi = hNew->GetXaxis()->GetBinUpEdge(j);

        const double ilo = std::max(olo, nlo);
        const double ihi = std::min(ohi, nhi);
        if (ihi <= ilo) continue;

        const double frac = (ihi - ilo) / w;
        const double add  = content * frac;
        const double addE = err * frac; // approximate

        const double prev  = hNew->GetBinContent(j);
        const double prevE = hNew->GetBinError(j);

        hNew->SetBinContent(j, prev + add);
        hNew->SetBinError(j, std::sqrt(prevE*prevE + addE*addE));
      }
    }

    hNew->GetXaxis()->SetTitle("v_{z} [cm]");
    hNew->GetYaxis()->SetTitle("Counts");
    return hNew;
  }

  // =============================================================================
  // Inventory (recursive)
  // =============================================================================
  inline void CollectInventoryRecursive(TDirectory* dir, const string& prefix, vector<InvItem>& out)
  {
    if (!dir) return;
    TIter next(dir->GetListOfKeys());
    while (TKey* key = (TKey*)next())
    {
      const string name = key->GetName();
      const string cls  = key->GetClassName();

      const string full = prefix + name;

      if (cls == "TDirectoryFile" || cls == "TDirectory")
      {
        TDirectory* sub = dynamic_cast<TDirectory*>(dir->Get(name.c_str()));
        if (sub) CollectInventoryRecursive(sub, full + "/", out);
        continue;
      }

      TObject* obj = dir->Get(name.c_str());
      double ent = -1.0;

      if (auto* h = dynamic_cast<TH1*>(obj)) ent = h->GetEntries();
      else if (auto* p = dynamic_cast<TProfile*>(obj)) ent = p->GetEntries();
      else if (auto* p2 = dynamic_cast<TProfile2D*>(obj)) ent = p2->GetEntries();
      else if (auto* p3 = dynamic_cast<TProfile3D*>(obj)) ent = p3->GetEntries();

      out.push_back({full, cls, ent});
    }
  }

  inline void PrintInventoryToTerminal(const Dataset& ds, const vector<InvItem>& items)
  {
    cout << ANSI_BOLD_CYN << "\n[INVENTORY] " << ds.label << " topDir=" << ds.topDirName
         << " (#objects=" << items.size() << ")" << ANSI_RESET << "\n";

    const int wPath = 92;
    const int wCls  = 18;
    const int wEnt  = 14;

    cout << std::left
         << std::setw(wPath) << "path"
         << std::setw(wCls)  << "class"
         << std::right << std::setw(wEnt)  << "entries"
         << "\n";
    cout << string(wPath + wCls + wEnt, '-') << "\n";

    for (const auto& it : items)
    {
      std::ostringstream ent;
      if (it.entries < 0) ent << "n/a";
      else ent << std::fixed << std::setprecision(0) << it.entries;

      string path = it.path;
      if ((int)path.size() > wPath-1) path = "..." + path.substr(path.size() - (wPath-4));

      cout << std::left
           << std::setw(wPath) << path
           << std::setw(wCls)  << it.cls
           << std::right << std::setw(wEnt) << ent.str()
           << "\n";
    }
  }

  // =============================================================================
  // Shared 3x3 table for pT-sliced TH1 families
  // =============================================================================
  inline void Make3x3Table_TH1(Dataset& ds,
                              const string& histBase,
                              const string& outDir,
                              const string& outName,
                              const string& xTitle,
                              const string& yTitle,
                              bool logy,
                              bool normalizeShape,
                              const vector<string>& commonLines)
  {
    EnsureDir(outDir);

    TCanvas c("c_tbl","c_tbl",1500,800);
    c.Divide(3,2, 0.001, 0.001);

    const auto& bins = PtBins();

    std::vector<TH1*> keepAlive;
    keepAlive.reserve(kNPtBins);

    for (int i = 0; i < kNPtBins; ++i)
    {
      c.cd(i+1);
      gPad->SetLeftMargin(0.14);
      gPad->SetRightMargin(0.05);
      gPad->SetBottomMargin(0.14);
      gPad->SetTopMargin(0.10);
      gPad->SetLogy(logy);

      const PtBin& b = bins[i];
      const string hname = histBase + b.suffix;

        // ------------------------------------------------------------
        // ISO tables get a compact label style:
        //   - Top-left: title (E_T^{iso, ...})
        //   - Top-right: big pT bin
        //   - Under it: smaller |eta| cut
        // Non-ISO tables keep the old multi-line header behavior.
        // ------------------------------------------------------------
        const bool isIsoTable =
          (histBase.rfind("h_Eiso", 0) == 0);  // starts with "h_Eiso"

        auto IsoTitleForBase = [&](const std::string& base)->std::string
        {
          if (base == "h_Eiso_emcal")  return "E_{T}^{iso, EMCal}";
          if (base == "h_Eiso_hcalin") return "E_{T}^{iso, IHCal}";
          if (base == "h_Eiso_hcalout")return "E_{T}^{iso, OHCal}";
          if (base == "h_Eiso")        return "E_{T}^{iso, Total}";
          return "E_{T}^{iso}";
        };

        auto DrawIsoPadLabels = [&](const PtBin& pb)->void
        {
          TLatex t;
          t.SetNDC(true);
          t.SetTextFont(42);

          // Title (top-left)
          t.SetTextAlign(13);          // left, top
          t.SetTextSize(0.08);
          t.DrawLatex(0.45, 0.985, IsoTitleForBase(histBase).c_str());

          // pT bin (top-right, larger)
          t.SetTextAlign(33);          // right, top
          t.SetTextSize(0.070);
          t.DrawLatex(0.95, 0.88,
            TString::Format("p_{T}^{#gamma}: %d-%d GeV", pb.lo, pb.hi).Data());

          // eta cut (under pT, smaller)
          t.SetTextSize(0.055);
          t.DrawLatex(0.85, 0.75, "|#eta^{#gamma}| < 0.7");

          if (!ds.centLabel.empty())
          {
            t.SetTextSize(0.050);
            t.DrawLatex(0.95, 0.66, ds.centLabel.c_str());
          }
        };

        TH1* h = GetObj<TH1>(ds, hname, true, true, true);
        if (!h)
        {
          // Keep your "MISSING" message, but use the new compact label style for iso tables
          TLatex t;
          t.SetNDC(true);
          t.SetTextFont(42);
          t.SetTextAlign(22);          // center, center
          t.SetTextSize(0.075);
          t.DrawLatex(0.50, 0.55, "MISSING");

          if (isIsoTable)
          {
            DrawIsoPadLabels(b);
          }
          else
          {
              // old behavior for non-iso tables
              vector<string> lines = commonLines;
              if (!ds.centLabel.empty()) lines.push_back(ds.centLabel);
              {
                std::ostringstream s;
                s << histBase << "   pT^{#gamma}: " << b.lo << "-" << b.hi << " GeV";
                lines.push_back(s.str());
              }
              DrawLatexLines(0.16, 0.90, lines, 0.040, 0.050);
          }

          continue;
        }

        TH1* hc = CloneTH1(h, TString::Format("%s_tbl_%d", histBase.c_str(), i).Data());
        if (!hc) continue;

        // Keep the drawn hist alive until after SaveCanvas()
        hc->SetDirectory(nullptr);

        if (normalizeShape) NormalizeToUnitArea(hc);

        hc->SetLineWidth(2);
        hc->SetTitle("");

        // Axis titles
        hc->GetXaxis()->SetTitle(xTitle.c_str());
        hc->GetYaxis()->SetTitle(yTitle.c_str());

        // ---- Make axis text larger (especially X) and "zoom in" by reducing offsets ----
        // X axis: bigger tick labels + bigger title, closer to axis
        hc->GetXaxis()->SetLabelSize(0.060);
        hc->GetXaxis()->SetTitleSize(0.070);
        hc->GetXaxis()->SetTitleOffset(0.85);

        // Y axis: also slightly larger for consistency on small pads
        hc->GetYaxis()->SetLabelSize(0.055);
        hc->GetYaxis()->SetTitleSize(0.065);
        hc->GetYaxis()->SetTitleOffset(1.05);

        // (Optional) clearer tick marks on small pads
        hc->GetXaxis()->SetTickLength(0.03);
        hc->GetYaxis()->SetTickLength(0.03);

        // Force x-range ONLY for the total isolation 3x3 table:
        //   table3x3_Eiso_total.png  (histBase == "h_Eiso")
        const bool isTotalIsoTable = (isIsoTable && histBase == "h_Eiso");
        if (isTotalIsoTable)
        {
          hc->GetXaxis()->SetRangeUser(-2.0, 6.0);
        }

        hc->Draw("hist");

        // ------------------------------------------------------------------
        // For the total isolation table: draw ONE iso-cut line at the pT-bin CENTER
        //   cut(ptCenter) = 1.08128 + 0.0299107 * ptCenter
        // and integrate counts left/right of that ONE threshold.
        // ------------------------------------------------------------------
        if (isTotalIsoTable)
        {
          constexpr double kIsoA = 1.08128;
          constexpr double kIsoB = 0.0299107;

          const double ptCenter = 0.5 * (static_cast<double>(b.lo) + static_cast<double>(b.hi));
          const double cut = kIsoA + kIsoB * ptCenter;

          // y-range for the vertical line (avoid y=0 on log plots)
          double yBot = 0.0;
          if (logy)
          {
            const double minPos = SmallestPositiveBinContent(hc);
            yBot = (minPos > 0.0) ? (0.5 * minPos) : 1e-6;
          }

          const double yTop = (hc->GetMaximum() > 0.0) ? (hc->GetMaximum() * 1.25) : 1.0;

          // Draw ONE line only (clone so it persists until SaveCanvas)
          {
            TLine l(cut, yBot, cut, yTop);
            l.SetLineWidth(2);
            l.SetLineStyle(2); // dashed
            l.SetLineColor(kRed + 1);
            l.DrawClone("same");
          }

          // Integrals (under/overflow included) relative to ONE cut
          const int nbx = hc->GetNbinsX();
          const double nTot = hc->Integral(0, nbx + 1);

          const int bCut = hc->GetXaxis()->FindBin(cut);

          const double nLeft  = hc->Integral(0, bCut);            // <= cut bin
          const double nRight = hc->Integral(bCut + 1, nbx + 1);  // > cut bin

          // Compact annotations (left side so it won't collide with pT label)
          TLatex t;
          t.SetNDC(true);
          t.SetTextFont(42);
          t.SetTextAlign(13);

          t.SetTextSize(0.052);
          t.DrawLatex(
            0.16, 0.70,
            TString::Format("iso cut @ pT center: %.2f GeV", cut).Data()
          );

          t.SetTextSize(0.048);
          t.DrawLatex(
            0.16, 0.62,
            TString::Format("N(<=cut)=%.0f  N(>cut)=%.0f", nLeft, nRight).Data()
          );

          if (nTot > 0.0)
          {
            t.DrawLatex(
              0.16, 0.54,
              TString::Format("f(<=cut)=%.3f  f(>cut)=%.3f", nLeft / nTot, nRight / nTot).Data()
            );
          }
        }

        if (isIsoTable)
        {
          // new compact iso labels
          DrawIsoPadLabels(b);
        }
        else
        {
            // old behavior for everything else
            vector<string> lines = commonLines;
            if (!ds.centLabel.empty()) lines.push_back(ds.centLabel);
            {
              std::ostringstream s;
              s << histBase << "   pT^{#gamma}: " << b.lo << "-" << b.hi << " GeV";
              lines.push_back(s.str());
            }
            DrawLatexLines(0.16, 0.90, lines, 0.040, 0.050);
        }

        keepAlive.push_back(hc);
    }

      SaveCanvas(c, JoinPath(outDir, outName));

      // Now it is safe to delete the histograms drawn into pads
      for (TH1* htmp : keepAlive) delete htmp;
      keepAlive.clear();
  }

  // =============================================================================
  // Leakage corrected solver
  // =============================================================================
  inline bool SolveLeakageCorrectedSA(double A, double B, double C, double D,
                                     double fB, double fC, double fD,
                                     double& outSA)
  {
    outSA = 0.0;
    if (A <= 0.0) { outSA = 0.0; return true; }

    double S = A;
    if (D != 0.0)
    {
      const double Asig_raw = A - B*(C/D);
      S = std::min(std::max(Asig_raw, 0.0), A);
    }

    auto F = [&](double s)->double
    {
      const double denom = (D - fD*s);
      if (denom == 0.0) return std::numeric_limits<double>::quiet_NaN();
      return A - (B - fB*s)*(C - fC*s)/denom;
    };

    const int maxIter = 200;
    const double tol = 1e-6;
    double lambda = 0.25;

    for (int it = 0; it < maxIter; ++it)
    {
      if (fD > 0.0)
      {
        const double sMax = (D / fD) * 0.999;
        if (std::isfinite(sMax)) S = std::min(S, std::max(0.0, sMax));
      }

      const double fval = F(S);
      if (!std::isfinite(fval)) return false;

      double Snew = (1.0 - lambda)*S + lambda*fval;
      if (!std::isfinite(Snew)) return false;

      Snew = std::min(std::max(Snew, 0.0), A);

      const double delta = std::fabs(Snew - S);
      S = Snew;

      if (delta < tol)
      {
        outSA = S;
        return true;
      }

      if (it > 10 && delta > 0.5*A && lambda > 0.05) lambda *= 0.5;
    }

    outSA = S;
    return false;
  }

  // =============================================================================
  // Tiny shared I/O helpers (used in several sections)
  // =============================================================================
  inline void WriteTextFile(const string& filepath, const vector<string>& lines)
  {
    EnsureParentDirForFile(filepath);
    std::ofstream out(filepath.c_str());
    for (const auto& s : lines) out << s << "\n";
  }

  inline void WriteJetSummaryTxt(const string& filepath,
                                const string& rKey,
                                double Nevt,
                                const map<string,double>& scalars,
                                const vector<string>& notes)
  {
    EnsureParentDirForFile(filepath);
    std::ofstream out(filepath.c_str());
    out << "rKey: " << rKey << "\n";
    out << "Nevt (from h_HT entries): " << std::fixed << std::setprecision(0) << Nevt << "\n\n";
    out << "[SCALARS]\n";
    for (const auto& kv : scalars)
    {
      out << " " << kv.first << " = " << std::fixed << std::setprecision(6) << kv.second << "\n";
    }
    if (!notes.empty())
    {
      out << "\n[NOTES]\n";
      for (const auto& s : notes) out << " - " << s << "\n";
    }
  }

  // =============================================================================
  // Tiny shared counters
  // =============================================================================
  inline double ReadEventCount(Dataset& ds)
  {
    const string cntName = "cnt_" + ds.topDirName;
    TH1* cnt = GetObj<TH1>(ds, cntName, true, true, false);
    if (!cnt) return 0.0;
    return cnt->GetBinContent(1);
  }

  inline double Read1BinCount(Dataset& ds, const string& hname)
  {
    TH1* h = GetObj<TH1>(ds, hname, true, true, false);
    if (!h) return 0.0;
    return h->GetBinContent(1);
  }

  inline double DetermineNevtForRKey(Dataset& ds, const string& rKey, double fallback)
  {
    TH1* h = GetObj<TH1>(ds, "h_HT_" + rKey, true, true, false);
    if (!h) return fallback;
    const double Nevt = h->GetEntries();
    return (Nevt > 0.0 ? Nevt : fallback);
  }

  // =============================================================================
  // MatchCache helpers
  // =============================================================================
  inline void InitMatchCache(MatchCache& mc)
  {
    mc.NphoLead.clear();
    mc.NphoMatched.clear();
    for (const auto& rKey : kRKeys)
    {
      mc.NphoLead[rKey]    = vector<double>(kNPtBins, 0.0);
      mc.NphoMatched[rKey] = vector<double>(kNPtBins, 0.0);
    }
  }

  inline int FindPtBinIndexByEdges(int lo, int hi)
  {
    const auto& bins = PtBins();
    for (int i = 0; i < kNPtBins; ++i)
    {
      if (bins[i].lo == lo && bins[i].hi == hi) return i;
    }
    return -1;
  }

  // =============================================================================
  // TH3 / Profile3D projection helpers
  // =============================================================================
  inline string AxisBinLabel(const TAxis* ax, int ibin, const string& unit = "", int prec = 0)
  {
    if (!ax) return "";
    const double lo = ax->GetBinLowEdge(ibin);
    const double hi = ax->GetBinUpEdge(ibin);

    std::ostringstream s;
    s << std::fixed << std::setprecision(prec) << lo << "-" << hi;
    if (!unit.empty()) s << " " << unit;
    return s.str();
  }

  inline TH2* ProjectYZ_AtXbin_TH3(const TH3* h3, int xbin, const string& newName)
  {
    if (!h3) return nullptr;
    TH3* hc = CloneTH3(h3, newName + "_clone3");
    if (!hc) return nullptr;
    hc->GetXaxis()->SetRange(xbin, xbin);

    TH2* h2 = dynamic_cast<TH2*>(hc->Project3D("zy"));
    if (h2) h2->SetDirectory(nullptr);

    delete hc;
    return h2;
  }

  inline TH1* ProjectY_AtXbin_TH3(const TH3* h3, int xbin, const string& newName)
  {
    if (!h3) return nullptr;
    TH3* hc = CloneTH3(h3, newName + "_clone3y");
    if (!hc) return nullptr;
    hc->GetXaxis()->SetRange(xbin, xbin);

    TH1* h1 = dynamic_cast<TH1*>(hc->Project3D("y"));
    if (h1) h1->SetDirectory(nullptr);

    delete hc;
    return h1;
  }

  // ---------------------------------------------------------------------------
  // Helper: build overlap weights mapping old X-axis binning onto a target [xLo,xHi]
  // ---------------------------------------------------------------------------
  inline bool XaxisOverlapWeights(const TAxis* ax,
                                    double xLo,
                                    double xHi,
                                    int& outBinLo,
                                    int& outBinHi,
                                    std::vector<double>& outW)
  {
      outW.clear();
      outBinLo = -1;
      outBinHi = -1;
      if (!ax) return false;
      if (xHi <= xLo) return false;

      const int nb = ax->GetNbins();
      int bLo = -1, bHi = -1;

      // find all old bins that overlap [xLo,xHi]
      for (int ib = 1; ib <= nb; ++ib)
      {
        const double lo = ax->GetBinLowEdge(ib);
        const double hi = ax->GetBinUpEdge(ib);
        const double oLo = std::max(lo, xLo);
        const double oHi = std::min(hi, xHi);
        if (oHi > oLo)
        {
          if (bLo < 0) bLo = ib;
          bHi = ib;
        }
      }
      if (bLo < 0 || bHi < 0) return false;

      outBinLo = bLo;
      outBinHi = bHi;

      const double widthTarget = (xHi - xLo);
      outW.reserve(bHi - bLo + 1);

      for (int ib = bLo; ib <= bHi; ++ib)
      {
        const double lo = ax->GetBinLowEdge(ib);
        const double hi = ax->GetBinUpEdge(ib);
        const double oLo = std::max(lo, xLo);
        const double oHi = std::min(hi, xHi);
        const double overlap = std::max(0.0, oHi - oLo);
        outW.push_back(overlap / widthTarget);
      }

      return true;
  }

  // ---------------------------------------------------------------------------
  // Helper: sum Y projections over an X-bin range with explicit weights (weights are per-bin)
  // ---------------------------------------------------------------------------
  inline TH1* ProjectY_AtXrange_TH3_Weighted(const TH3* h3,
                                              int xbinLo,
                                              int xbinHi,
                                              const std::vector<double>& w,
                                              const std::string& newName)
  {
      if (!h3) return nullptr;
      if (xbinHi < xbinLo) return nullptr;

      TH1* sum = nullptr;

      for (int xb = xbinLo; xb <= xbinHi; ++xb)
      {
        TH1* h = ProjectY_AtXbin_TH3(h3, xb, newName + TString::Format("_xb%d", xb).Data());
        if (!h) continue;

        const int iw = xb - xbinLo;
        const double ww = (iw >= 0 && iw < (int)w.size()) ? w[iw] : 1.0;
        h->Scale(ww);

        if (!sum)
        {
          sum = CloneTH1(h, newName);
          if (sum)
          {
            sum->Reset("ICES");
            sum->SetDirectory(nullptr);
          }
        }
        if (sum) sum->Add(h);
        delete h;
      }

      return sum;
  }


  inline TH1* ProjectY_AtXbin_AndZRange_TH3(const TH3* h3,
                                           int xbin,
                                           int zbinLo,
                                           int zbinHi,
                                           const string& newName)
  {
    if (!h3) return nullptr;
    if (zbinHi < zbinLo) return nullptr;

    TH3* hc = CloneTH3(h3, newName + "_clone3y_zrng");
    if (!hc) return nullptr;

    hc->GetXaxis()->SetRange(xbin, xbin);
    hc->GetZaxis()->SetRange(zbinLo, zbinHi);

    TH1* h1 = dynamic_cast<TH1*>(hc->Project3D("y"));
    if (h1) h1->SetDirectory(nullptr);

    delete hc;
    return h1;
  }

  inline TH1* ProjectY_AtXbin_AndAlphaMax_TH3(const TH3* h3,
                                             int xbin,
                                             double alphaMax,
                                             const string& newName)
  {
    if (!h3) return nullptr;

    const TAxis* az = h3->GetZaxis();
    if (!az) return nullptr;

    const int nb = az->GetNbins();
    if (nb <= 0) return nullptr;

    const double zmin = az->GetXmin();
    const double zmax = az->GetXmax();
    const double eps  = 1e-9;

    if (alphaMax <= zmin) return nullptr;

    double a = alphaMax;
    if (a > zmax) a = zmax;

    int zLo = 1;
    int zHi = az->FindBin(a - eps);

    if (zHi < 1)  zHi = 1;
    if (zHi > nb) zHi = nb;
    if (zHi < zLo) return nullptr;

    return ProjectY_AtXbin_AndZRange_TH3(h3, xbin, zLo, zHi, newName);
  }

  inline TH1* ProjectZ_AtXbin_TH3(const TH3* h3, int xbin, const string& newName)
  {
    if (!h3) return nullptr;
    TH3* hc = CloneTH3(h3, newName + "_clone3z");
    if (!hc) return nullptr;
    hc->GetXaxis()->SetRange(xbin, xbin);

    TH1* h1 = dynamic_cast<TH1*>(hc->Project3D("z"));
    if (h1) h1->SetDirectory(nullptr);

    delete hc;
    return h1;
  }

  inline TProfile2D* ProjectYZ_AtXbin_Profile3D(const TProfile3D* p3, int xbin, const string& newName)
  {
    if (!p3) return nullptr;
    TProfile3D* pc = (TProfile3D*)p3->Clone((newName + "_cloneP3").c_str());
    if (!pc) return nullptr;
    pc->SetDirectory(nullptr);
    pc->GetXaxis()->SetRange(xbin, xbin);

    TProfile2D* p2 = dynamic_cast<TProfile2D*>(pc->Project3DProfile("zy"));
    if (p2) p2->SetDirectory(nullptr);

    delete pc;
    return p2;
  }

  // =============================================================================
  // SIM merge utilities (photonJet slices: 5 / 10 / 20)
  // =============================================================================
  inline double ReadEventCountFromFile(TFile* f, const string& topDirName)
  {
      if (!f) return 0.0;
      TDirectory* d = f->GetDirectory(topDirName.c_str());
      if (!d) return 0.0;

      TH1* cnt = dynamic_cast<TH1*>(d->Get(("cnt_" + topDirName).c_str()));
      if (!cnt) return 0.0;
      return cnt->GetBinContent(1);
  }

    // Existing 2-slice recursive merge helper (kept for backwards compatibility)
  inline void CopyAndScaleAddRecursive(TDirectory* outDir,
                                         TDirectory* d10, double w10,
                                         TDirectory* d20, double w20)
  {
      if (!outDir) return;
      if (!d10 && !d20) return;

      auto IsDirClass = [](const std::string& cls)->bool
      {
        return (cls == "TDirectoryFile" || cls == "TDirectory");
      };

      std::set<std::string> seen;

      // Pass 1: keys in d10
      if (d10)
      {
        TIter next(d10->GetListOfKeys());
        while (TKey* key = (TKey*)next())
        {
          const std::string name = key->GetName();
          const std::string cls  = key->GetClassName();
          seen.insert(name);

          if (IsDirClass(cls))
          {
            TDirectory* sub10 = dynamic_cast<TDirectory*>(d10->Get(name.c_str()));
            TDirectory* sub20 = (d20 ? dynamic_cast<TDirectory*>(d20->Get(name.c_str())) : nullptr);

            outDir->cd();
            TDirectory* subOut = outDir->mkdir(name.c_str());
            if (!subOut) continue;

            CopyAndScaleAddRecursive(subOut, sub10, w10, sub20, w20);
            continue;
          }

          TObject* o10 = d10->Get(name.c_str());
          TObject* o20 = (d20 ? d20->Get(name.c_str()) : nullptr);

          if (auto* h10 = dynamic_cast<TH1*>(o10))
          {
            outDir->cd();
            TH1* h = dynamic_cast<TH1*>(h10->Clone(name.c_str()));
            if (!h) continue;
            h->SetDirectory(outDir);
            if (h->GetSumw2N() == 0) h->Sumw2();

            if (w10 != 0.0) h->Scale(w10);

            if (auto* h20 = dynamic_cast<TH1*>(o20))
            {
              TH1* tmp = dynamic_cast<TH1*>(h20->Clone((name + "_tmp20").c_str()));
              if (tmp)
              {
                tmp->SetDirectory(nullptr);
                if (tmp->GetSumw2N() == 0) tmp->Sumw2();
                if (w20 != 0.0) tmp->Scale(w20);
                h->Add(tmp);
                delete tmp;
              }
            }

            h->Write(name.c_str(), TObject::kOverwrite);
            continue;
          }

          outDir->cd();
          if (o10)
          {
            TObject* c = o10->Clone(name.c_str());
            if (c) c->Write(name.c_str(), TObject::kOverwrite);
          }
        }
      }

      // Pass 2: keys in d20 not in d10
      if (d20)
      {
        TIter next(d20->GetListOfKeys());
        while (TKey* key = (TKey*)next())
        {
          const std::string name = key->GetName();
          const std::string cls  = key->GetClassName();
          if (seen.count(name)) continue;

          if (IsDirClass(cls))
          {
            TDirectory* sub20 = dynamic_cast<TDirectory*>(d20->Get(name.c_str()));
            if (!sub20) continue;

            outDir->cd();
            TDirectory* subOut = outDir->mkdir(name.c_str());
            if (!subOut) continue;

            CopyAndScaleAddRecursive(subOut, nullptr, 0.0, sub20, w20);
            continue;
          }

          TObject* o20 = d20->Get(name.c_str());

          if (auto* h20 = dynamic_cast<TH1*>(o20))
          {
            outDir->cd();
            TH1* h = dynamic_cast<TH1*>(h20->Clone(name.c_str()));
            if (!h) continue;
            h->SetDirectory(outDir);
            if (h->GetSumw2N() == 0) h->Sumw2();
            if (w20 != 0.0) h->Scale(w20);
            h->Write(name.c_str(), TObject::kOverwrite);
            continue;
          }

          outDir->cd();
          if (o20)
          {
            TObject* c = o20->Clone(name.c_str());
            if (c) c->Write(name.c_str(), TObject::kOverwrite);
          }
        }
      }
  }

  // Add one weighted slice into an *existing* output directory tree (supports N-slice merges)
  inline void AddScaledRecursive(TDirectory* outDir, TDirectory* inDir, double w)
  {
      if (!outDir || !inDir) return;

      auto IsDirClass = [](const std::string& cls)->bool
      {
        return (cls == "TDirectoryFile" || cls == "TDirectory");
      };

      TIter next(inDir->GetListOfKeys());
      while (TKey* key = (TKey*)next())
      {
        const std::string name = key->GetName();
        const std::string cls  = key->GetClassName();

        if (IsDirClass(cls))
        {
          TDirectory* subIn = dynamic_cast<TDirectory*>(inDir->Get(name.c_str()));
          if (!subIn) continue;

          outDir->cd();
          TDirectory* subOut = outDir->GetDirectory(name.c_str());
          if (!subOut) subOut = outDir->mkdir(name.c_str());
          if (!subOut) continue;

          AddScaledRecursive(subOut, subIn, w);
          continue;
        }

        TObject* objIn = inDir->Get(name.c_str());
        if (!objIn) continue;

        if (auto* hIn = dynamic_cast<TH1*>(objIn))
        {
          outDir->cd();
          TH1* hOut = dynamic_cast<TH1*>(outDir->Get(name.c_str()));

          if (!hOut)
          {
            TH1* hNew = dynamic_cast<TH1*>(hIn->Clone(name.c_str()));
            if (!hNew) continue;
            hNew->SetDirectory(outDir);
            if (hNew->GetSumw2N() == 0) hNew->Sumw2();
            if (w != 0.0) hNew->Scale(w);
            hNew->Write(name.c_str(), TObject::kOverwrite);
          }
          else
          {
            TH1* tmp = dynamic_cast<TH1*>(hIn->Clone((name + "_tmpAdd").c_str()));
            if (!tmp) continue;
            tmp->SetDirectory(nullptr);
            if (tmp->GetSumw2N() == 0) tmp->Sumw2();
            if (w != 0.0) tmp->Scale(w);
            hOut->Add(tmp);
            hOut->Write(name.c_str(), TObject::kOverwrite);
            delete tmp;
          }

          continue;
        }

        // Non-hist objects: only copy if missing (avoid clobbering metadata repeatedly)
        outDir->cd();
        if (!outDir->Get(name.c_str()))
        {
          TObject* c = objIn->Clone(name.c_str());
          if (c) c->Write(name.c_str(), TObject::kOverwrite);
        }
      }
  }

  // Generic N-slice merge (use for 2-slice or 3-slice weighted merges)
  inline bool BuildMergedSIMFile_PhotonSlices(const vector<string>& inFiles,
                                                  const vector<double>& sigmas_pb,
                                                  const string& outMerged,
                                                  const string& topDirName,
                                                  const vector<string>& sliceLabels = {})
      {
        cout << ANSI_BOLD_CYN
             << "\n[MERGE SIM] Building merged SIM file with cross-section weights\n"
             << "  out    = " << outMerged << "\n"
             << "  topDir = " << topDirName << "\n"
             << ANSI_RESET;

        const size_t n = inFiles.size();
        if (n < 2 || sigmas_pb.size() != n)
        {
          cout << ANSI_BOLD_YEL
               << "[MERGE SIM][WARN] BuildMergedSIMFile_PhotonSlices needs >=2 inputs and a matching sigma list. Skipping this merge target.\n"
               << ANSI_RESET;
          return false;
        }

        vector<TFile*> fin(n, nullptr);
        vector<TDirectory*> din(n, nullptr);
        vector<double> N(n, 0.0);
        vector<double> w(n, 0.0);

        auto CloseAll = [&](){
          for (auto* f : fin) { if (f) f->Close(); }
        };

        for (size_t i = 0; i < n; ++i)
        {
          cout << "  in[" << i << "] = " << inFiles[i]
               << "   sigma_pb=" << std::setprecision(12) << sigmas_pb[i] << "\n";

          fin[i] = TFile::Open(inFiles[i].c_str(), "READ");
          if (!fin[i] || fin[i]->IsZombie())
          {
            cout << ANSI_BOLD_YEL
                 << "[MERGE SIM][WARN] Cannot open input SIM file. Skipping this merge target:\n"
                 << "  " << inFiles[i] << "\n"
                 << ANSI_RESET;
            CloseAll();
            return false;
          }

          din[i] = fin[i]->GetDirectory(topDirName.c_str());
          if (!din[i])
          {
            cout << ANSI_BOLD_YEL
                 << "[MERGE SIM][WARN] Missing topDir '" << topDirName << "' in input SIM file. Skipping this merge target:\n"
                 << "  " << inFiles[i] << "\n"
                 << ANSI_RESET;
            CloseAll();
            return false;
          }

          N[i] = ReadEventCountFromFile(fin[i], topDirName);
          if (N[i] <= 0.0)
          {
            cout << ANSI_BOLD_YEL
                 << "[MERGE SIM][WARN] Naccepted <= 0 for input. Skipping this merge target:\n"
                 << "  " << inFiles[i] << "\n"
                 << ANSI_RESET;
            CloseAll();
            return false;
          }

            w[i] = sigmas_pb[i] / N[i];
          }

          // Normalize by a common factor so merged histograms remain "count-like",
          // while preserving correct relative (sigma/N) weighting across slices.
          double wRef = 0.0;
          for (size_t i = 0; i < n; ++i)
          {
            if (w[i] > 0.0) { wRef = w[i]; break; }
          }
          if (wRef <= 0.0)
          {
            cout << ANSI_BOLD_YEL
                 << "[MERGE SIM][WARN] wRef <= 0 (cannot normalize slice weights). Skipping this merge target.\n"
                 << ANSI_RESET;
            CloseAll();
            return false;
          }
          for (size_t i = 0; i < n; ++i)
          {
            w[i] /= wRef;
          }

          cout << ANSI_BOLD_YEL << "[MERGE SIM] Slice weights (relative): w = (sigma/N)/(sigma0/N0)\n" << ANSI_RESET;
          for (size_t i = 0; i < n; ++i)
          {
            const string lab = (!sliceLabels.empty() && sliceLabels.size() == n) ? sliceLabels[i] : std::to_string(i);
            cout << "  [" << lab << "]  N=" << std::fixed << std::setprecision(0) << N[i]
                 << "   sigma_pb=" << std::setprecision(12) << sigmas_pb[i]
                 << "   w=" << std::setprecision(12) << w[i] << "\n";
          }

        EnsureParentDirForFile(outMerged);
        TFile* fout = TFile::Open(outMerged.c_str(), "RECREATE");
        if (!fout || fout->IsZombie())
        {
          cout << ANSI_BOLD_YEL << "[MERGE SIM][WARN] Cannot create merged output file. Skipping this merge target." << ANSI_RESET << "\n";
          CloseAll();
          return false;
        }

        fout->cd();
        TDirectory* outTop = fout->mkdir(topDirName.c_str());
        if (!outTop)
        {
          cout << ANSI_BOLD_YEL << "[MERGE SIM][WARN] Cannot create topDir in merged output. Skipping this merge target." << ANSI_RESET << "\n";
          fout->Close();
          CloseAll();
          return false;
        }

        // Build output by adding each slice into the same output tree
        for (size_t i = 0; i < n; ++i)
        {
          AddScaledRecursive(outTop, din[i], w[i]);
        }

        // Carry over / merge EventDisplayTree (diagnostic TTree) so merged SIM keeps eventDisplay capability.
        // NOTE: this is UNWEIGHTED: we simply concatenate entries across slices (diagnostics only).
        {
          TTree* tOutED = nullptr;

          for (size_t i = 0; i < n; ++i)
          {
            TTree* tInED = dynamic_cast<TTree*>(fin[i] ? fin[i]->Get("EventDisplayTree") : nullptr);
            if (!tInED && din[i])
            {
              tInED = dynamic_cast<TTree*>(din[i]->Get("EventDisplayTree"));
            }
            if (!tInED) continue;

            if (!tOutED)
            {
              outTop->cd();
              tOutED = tInED->CloneTree(0);
              if (tOutED) tOutED->SetDirectory(outTop);
            }

            if (tOutED)
            {
              tOutED->CopyEntries(tInED);
            }
          }

          if (tOutED)
          {
            outTop->cd();
            tOutED->Write("EventDisplayTree", TObject::kOverwrite);
          }
        }

        // Metadata
        std::ostringstream oss;
        oss << "Merged photonJet slices. Nslices=" << n << " ";
        for (size_t i = 0; i < n; ++i)
        {
          const string lab = (!sliceLabels.empty() && sliceLabels.size() == n) ? sliceLabels[i] : std::to_string(i);
          oss << "[" << lab
              << " N=" << std::fixed << std::setprecision(0) << N[i]
              << " sigma_pb=" << std::setprecision(12) << sigmas_pb[i]
              << " w=" << std::setprecision(12) << w[i]
              << "] ";
        }

        outTop->cd();
        TNamed meta("MERGE_INFO", oss.str().c_str());
        meta.Write("MERGE_INFO", TObject::kOverwrite);

        fout->Write();
        fout->Close();
        CloseAll();

        cout << ANSI_BOLD_CYN
             << "[MERGE SIM] Done. Merged file written: " << outMerged << "\n"
             << "[MERGE SIM] NOTE: histograms are weighted by relative (sigma/N) slice weights (overall normalization arbitrary).\n"
             << ANSI_RESET;

        return true;
  }

  // Backwards-compatible wrapper (existing call sites remain valid)
  inline bool BuildMergedSIMFile_Photon10And20(const string& in10,
                                                 const string& in20,
                                                 const string& outMerged,
                                                 const string& topDirName,
                                                 double sigma10_pb,
                                                 double sigma20_pb)
  {
      return BuildMergedSIMFile_PhotonSlices({in10, in20},
                                             {sigma10_pb, sigma20_pb},
                                             outMerged,
                                             topDirName,
                                             {"photonJet10", "photonJet20"});
  }

  // =============================================================================
  // Run mode + SIM sample helpers
  // =============================================================================
  enum class SimSample
  {
        kNone,
        kPhotonJet5,
        kPhotonJet10,
        kPhotonJet20,
        kPhotonJet5And10Merged,
        kPhotonJet5And20Merged,
        kPhotonJet10And20Merged,
        kPhotonJet5And10And20Merged,
        kSimMB,
        kSimJet5,
        kSimEmbedded,
        kInvalid
  };

  inline bool IsMergedSimSample(SimSample s)
  {
        return (s == SimSample::kPhotonJet5And10Merged ||
                s == SimSample::kPhotonJet5And20Merged ||
                s == SimSample::kPhotonJet10And20Merged ||
                s == SimSample::kPhotonJet5And10And20Merged);
  }

  inline SimSample CurrentSimSample()
  {
        const int nTrue =
          (isPhotonJet5 ? 1 : 0) +
          (isPhotonJet10 ? 1 : 0) +
          (isPhotonJet20 ? 1 : 0) +
          (bothPhoton5and10sim ? 1 : 0) +
          (bothPhoton5and20sim ? 1 : 0) +
          (bothPhoton10and20sim ? 1 : 0) +
          (allPhoton5and10and20sim ? 1 : 0) +
          (isSimMB ? 1 : 0) +
          (isSimJet5 ? 1 : 0) +
          (isSimEmbedded ? 1 : 0);

        if (nTrue == 0) return SimSample::kNone;
        if (nTrue != 1) return SimSample::kInvalid;

        if (isPhotonJet5)            return SimSample::kPhotonJet5;
        if (isPhotonJet10)           return SimSample::kPhotonJet10;
        if (isPhotonJet20)           return SimSample::kPhotonJet20;
        if (bothPhoton5and10sim)     return SimSample::kPhotonJet5And10Merged;
        if (bothPhoton5and20sim)     return SimSample::kPhotonJet5And20Merged;
        if (bothPhoton10and20sim)    return SimSample::kPhotonJet10And20Merged;
        if (allPhoton5and10and20sim) return SimSample::kPhotonJet5And10And20Merged;
        if (isSimMB)                 return SimSample::kSimMB;
        if (isSimJet5)               return SimSample::kSimJet5;
        if (isSimEmbedded)           return SimSample::kSimEmbedded;

        return SimSample::kInvalid;
  }

  inline string SimSampleLabel(SimSample s)
  {
          switch (s)
          {
            case SimSample::kNone:                      return "NONE";
            case SimSample::kPhotonJet5:                return "photonJet5";
            case SimSample::kPhotonJet10:               return "photonJet10";
            case SimSample::kPhotonJet20:               return "photonJet20";
            case SimSample::kPhotonJet5And10Merged:     return "photonJet5and10merged";
            case SimSample::kPhotonJet5And20Merged:     return "photonJet5and20merged";
            case SimSample::kPhotonJet10And20Merged:    return "photonJet10and20merged";
            case SimSample::kPhotonJet5And10And20Merged:return "photonJet5and10and20merged";
            case SimSample::kSimMB:                     return "simMB";
            case SimSample::kSimJet5:                   return "simJet5";
            case SimSample::kSimEmbedded:               return "embeddedPhoton20_ALL_" + CfgTagWithUE();
            default:                                    return "INVALID";
          }
    }

    inline string SimInputPathForSample(SimSample s)
    {
        switch (s)
        {
          case SimSample::kPhotonJet5:                return InputSim("photonjet5");
          case SimSample::kPhotonJet10:               return InputSim("photonjet10");
          case SimSample::kPhotonJet20:               return InputSim("photonjet20");
          case SimSample::kPhotonJet5And10Merged:     return MergedSimPath("photonJet5and10merged_SIM", "RecoilJets_photonjet5plus10_MERGED.root");
          case SimSample::kPhotonJet5And20Merged:     return MergedSimPath("photonJet5and20merged_SIM", "RecoilJets_photonjet5plus20_MERGED.root");
          case SimSample::kPhotonJet10And20Merged:    return MergedSimPath("photonJet10and20merged_SIM", "RecoilJets_photonjet10plus20_MERGED.root");
          case SimSample::kPhotonJet5And10And20Merged:return MergedSimPath("photonJet5and10and20merged_SIM", "RecoilJets_photonjet5plus10plus20_MERGED.root");
          case SimSample::kSimMB:                     return InputSimMB();
          case SimSample::kSimJet5:                   return InputSimJet5();
          case SimSample::kSimEmbedded:               return InputSimEmbedded();
          default:                                    return "";
        }
    }

    inline string SimOutBaseForSample(SimSample s)
    {
          switch (s)
          {
            case SimSample::kPhotonJet5:                return OutputIndividualSim("photonJet5_SIM");
            case SimSample::kPhotonJet10:               return OutputIndividualSim("photonJet10_SIM");
            case SimSample::kPhotonJet20:               return OutputIndividualSim("photonJet20_SIM");
            case SimSample::kPhotonJet5And10Merged:     return OutputCombinedSimOnly("photonJet5and10merged_SIM");
            case SimSample::kPhotonJet5And20Merged:     return OutputCombinedSimOnly("photonJet5and20merged_SIM");
            case SimSample::kPhotonJet10And20Merged:    return OutputCombinedSimOnly("photonJet10and20merged_SIM");
            case SimSample::kPhotonJet5And10And20Merged:return OutputCombinedSimOnly("photonJet5and10and20merged_SIM");
            case SimSample::kSimMB:                     return OutputSimMB();
            case SimSample::kSimJet5:                   return OutputSimJet5();
            case SimSample::kSimEmbedded:               return OutputSimEmbedded();
            default:                                    return "";
          }
    }

    inline bool ValidateRunConfig(string* errMsg = nullptr)
    {
        // Disallow contradictory run-mode toggles.
        if (isPPdataOnly && isSimAndDataPP)
        {
          if (errMsg) *errMsg = "Both isPPdataOnly and isSimAndDataPP are true. Choose only one.";
          return false;
        }

        if (isAuAuOnly && (isPPdataOnly || isSimAndDataPP))
        {
          if (errMsg) *errMsg = "isAuAuOnly=true is mutually exclusive with isPPdataOnly and isSimAndDataPP. Choose only one run mode.";
          return false;
        }

        const SimSample ss = CurrentSimSample();

        // AuAu-only run: SIM sample toggles must be OFF.
        if (isAuAuOnly)
        {
          if (ss != SimSample::kNone)
          {
            if (errMsg)
            {
              *errMsg =
                "AuAu-only mode selected, but a SIM sample toggle is set. "
                "Set isPhotonJet5=false, isPhotonJet10=false, isPhotonJet20=false, "
                "bothPhoton5and10sim=false, bothPhoton5and20sim=false, bothPhoton10and20sim=false, "
                "allPhoton5and10and20sim=false, isSimMB=false, isSimJet5=false, isSimEmbedded=false.";
            }
            return false;
          }
          return true;
        }

        // PP-only run: SIM sample toggles must be OFF.
        if (isPPdataOnly)
        {
          if (ss != SimSample::kNone)
          {
            if (errMsg)
            {
              *errMsg =
                "PP-data-only mode selected, but a SIM sample toggle is set. "
                "Set isPhotonJet5=false, isPhotonJet10=false, isPhotonJet20=false, "
                "bothPhoton5and10sim=false, bothPhoton5and20sim=false, bothPhoton10and20sim=false, "
                "allPhoton5and10and20sim=false, isSimMB=false, isSimJet5=false, isSimEmbedded=false.";
            }
            return false;
          }
          return true;
        }

        // Any non-PP-only run includes SIM (either SIM-only or SIM+PP).
        if (ss == SimSample::kNone)
        {
          if (errMsg)
          {
            *errMsg =
              "SIM is required (SIM-only or SIM+PP), but no SIM sample was selected. "
              "Set exactly one of: isPhotonJet5, isPhotonJet10, isPhotonJet20, "
              "bothPhoton5and10sim, bothPhoton5and20sim, bothPhoton10and20sim, allPhoton5and10and20sim, "
              "isSimMB, isSimJet5, isSimEmbedded.";
          }
          return false;
        }

        if (ss == SimSample::kInvalid)
        {
          if (errMsg)
          {
            *errMsg =
              "Invalid SIM sample toggle combination. "
              "Set EXACTLY ONE of: isPhotonJet5, isPhotonJet10, isPhotonJet20, "
              "bothPhoton5and10sim, bothPhoton5and20sim, bothPhoton10and20sim, allPhoton5and10and20sim, "
              "isSimMB, isSimJet5, isSimEmbedded.";
          }
          return false;
        }

        // SIM+DATA mode contract: in-situ calibration + combined steps are defined for
        // the merged photonJet10+20 or photonJet5+10+20 SIM sample.
        if (isSimAndDataPP &&
            ss != SimSample::kPhotonJet10And20Merged &&
            ss != SimSample::kPhotonJet5And10And20Merged &&
            ss != SimSample::kSimMB &&
            ss != SimSample::kSimJet5)
        {
        if (errMsg)
        {
              *errMsg =
              "SIM+DATA (isSimAndDataPP=true) requires a merged photonJet10+20 or photonJet5+10+20 SIM sample, "
              "or isSimMB=true or isSimJet5=true. "
              "Set: bothPhoton10and20sim=true or allPhoton5and10and20sim=true or isSimMB=true or isSimJet5=true "
              "(all other SIM sample toggles false).";
          }
          return false;
        }

        return true;
    }

    // Backwards-compatible name used by older guard code paths.
    inline bool ExactlyOneModeSet()
    {
        return ValidateRunConfig(nullptr);
    }

    // =============================================================================
    // Run mode helpers
    // =============================================================================
    enum class RunMode
    {
          kPPDataOnly,
          kAuAuOnly,
          kSimOnly,
          kSimAndDataPP,
          kInvalid
    };

    inline RunMode CurrentRunMode()
    {
          if (!ValidateRunConfig(nullptr)) return RunMode::kInvalid;
          if (isAuAuOnly)    return RunMode::kAuAuOnly;
          if (isPPdataOnly)  return RunMode::kPPDataOnly;
          if (isSimAndDataPP) return RunMode::kSimAndDataPP;
          return RunMode::kSimOnly;
    }

    inline string RunModeLabel(RunMode m)
    {
        switch (m)
         {
            case RunMode::kPPDataOnly:   return "PP_DATA_ONLY";
            case RunMode::kAuAuOnly:     return "AUAU_ONLY";
            case RunMode::kSimOnly:      return "SIM_ONLY";
            case RunMode::kSimAndDataPP: return "SIM_AND_DATA_PP";
            default:                     return "INVALID";
        }
    }


} // namespace ARJ

#endif // ANALYZE_RECOIL_JETS_H
