// AnalyzeRecoilJets.h
// -----------------------------------------------------------------------------
// Toolbox + shared infrastructure for AnalyzeRecoilJets.cpp
// -----------------------------------------------------------------------------

#ifndef ANALYZE_RECOIL_JETS_H
#define ANALYZE_RECOIL_JETS_H

// =============================================================================
// ROOT
// =============================================================================
#include <TFile.h>
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
#include <TStyle.h>
#include <TSystem.h>
#include <TLine.h>
#include <TGraph.h>
#include <TGraphErrors.h>
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

// =============================================================================
// HOW TO USE THESE TOGGLES (pp data + photonJet SIM: 5 / 10 / 20)
// =============================================================================
//
// You control TWO things:
//
// (A) RUN MODE (PP-only vs SIM-only vs SIM+PP):
//
//   1) PP DATA ONLY (no SIM at all)
//        isPPdataOnly   = true;
//        isSimAndDataPP = false;
//        // IMPORTANT: all SIM sample toggles below MUST be false.
//
//   2) SIM ONLY (no PP data)
//        isPPdataOnly   = false;
//        isSimAndDataPP = false;
//        // choose EXACTLY ONE SIM sample toggle below.
//
//   3) SIM + PP DATA (run both; enables SIM/Data combined steps where implemented)
//        isPPdataOnly   = false;
//        isSimAndDataPP = true;
//        // choose EXACTLY ONE SIM sample toggle below.
//
// (B) WHICH SIM SAMPLE TO USE (choose EXACTLY ONE when SIM is included):
//
//   --- Single-slice SIM (raw event-count histograms) ---
//   - photonJet5 only:
//        isPhotonJet5 = true;
//        // all other SIM toggles false
//        -> SIM input:  kInSIM5
//        -> SIM output: kOutSIM5Base  (.../photonJet5_SIM)
//
//   - photonJet10 only:
//        isPhotonJet10 = true;
//        // all other SIM toggles false
//        -> SIM input:  kInSIM10
//        -> SIM output: kOutSIM10Base (.../photonJet10_SIM)
//
//   - photonJet20 only:
//        isPhotonJet20 = true;
//        // all other SIM toggles false
//        -> SIM input:  kInSIM20
//        -> SIM output: kOutSIM20Base (.../photonJet20_SIM)
//
//   --- Weighted merged SIM (histograms become weighted; y-axis ~ "Counts / pb^{-1}") ---
//   When a merged mode is selected, the code builds/uses a merged ROOT file whose histograms are:
//        H_merged = sum_i ( (sigma_i / N_i) * H_i )
//   where N_i is read from cnt_SIM bin 1 in each slice file (accepted event count).
//
//   - merged photonJet5+10:
//        bothPhoton5and10sim = true;
//        // all other SIM toggles false
//        -> merged output file: kMergedSIMOut_5and10
//        -> SIM output base:    kOutSIM5and10MergedBase (.../photonJet5and10merged_SIM)
//
//   - merged photonJet5+20:
//        bothPhoton5and20sim = true;
//        // all other SIM toggles false
//        -> merged output file: kMergedSIMOut_5and20
//        -> SIM output base:    kOutSIM5and20MergedBase (.../photonJet5and20merged_SIM)
//
//   - merged photonJet10+20:
//        bothPhoton10and20sim = true;
//        // all other SIM toggles false
//        -> merged output file: kMergedSIMOut
//        -> SIM output base:    kOutSIMMergedBase (.../photonJet10and20merged_SIM)
//
//   - merged photonJet5+10+20:
//        allPhoton5and10and20sim = true;
//        // all other SIM toggles false
//        -> merged output file: kMergedSIMOut_5and10and20
//        -> SIM output base:    kOutSIM5and10and20MergedBase (.../photonJet5and10and20merged_SIM)
//
// Cross sections used for weights (pb):
//   - photonJet5  : kSigmaPhoton5_pb  = 89266.571
//   - photonJet10 : kSigmaPhoton10_pb = 6692.7611
//   - photonJet20 : kSigmaPhoton20_pb = 105.79868
//
// INVALID COMBINATIONS (hard error):
//   - Setting more than one SIM sample toggle true.
//   - Setting any SIM sample toggle true while isPPdataOnly=true.
//   - Setting isPPdataOnly=true AND isSimAndDataPP=true.
//
// NOTES:
//   - For merged modes, IsWeightedSIMSelected() returns true, so plotting helpers label y-axes as "Counts / pb^{-1}".
//   - For single-slice SIM, y-axes remain raw "Counts" (unweighted).
// =============================================================================

  inline bool isPPdataOnly   = false;
  inline bool isSimAndDataPP = false;

  // SIM sample selection toggles (choose EXACTLY ONE for any SIM-including run)
  inline bool isPhotonJet5               = false;
  inline bool isPhotonJet10              = true;
  inline bool isPhotonJet20              = false;

  inline bool bothPhoton5and10sim        = false;
  inline bool bothPhoton5and20sim        = false;
  inline bool bothPhoton10and20sim       = false;

  inline bool allPhoton5and10and20sim    = false;

  // True if the selected SIM sample is a weighted multi-slice merge (hist units become ~pb/bin)
  inline bool IsWeightedSIMSelected()
  {
      return (bothPhoton5and10sim || bothPhoton5and20sim || bothPhoton10and20sim || allPhoton5and10and20sim);
  }

  // Displayed range [-vzCutCm,+vzCutCm] and 0.5 cm display bin width
  inline double vzCutCm = 30.0;

  // =============================================================================
  // FIXED INPUTS (pp + photonJet5/10/20 SIM)
  // =============================================================================
  inline const string kTriggerPP = "Photon_4_GeV_plus_MBD_NS_geq_1";
  inline const string kDirSIM    = "SIM";

  inline const string kInPP =
        "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/pp/RecoilJets_pp_ALL.root";

  inline const string kInSIM5 =
        "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/photonJet5_SIM/RecoilJets_photonjet5_ALL.root";

  inline const string kInSIM10 =
        "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/photonJet10_SIM/RecoilJets_photonjet10_ALL.root";

  inline const string kInSIM20 =
        "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/photonJet20_SIM/RecoilJets_photonjet20_ALL.root";

  inline const string kOutPPBase =
        "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/pp";

  // SIM outputs are routed by the SIM sample toggle(s) above:
  inline const string kOutSIM5Base =
        "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/photonJet5_SIM";

  inline const string kOutSIM10Base =
        "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/photonJet10_SIM";

  inline const string kOutSIM20Base =
        "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/photonJet20_SIM";

  // Keep the existing name for 10+20 merged (so older code stays readable):
  inline const string kOutSIMMergedBase =
        "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/photonJet10and20merged_SIM";

  inline const string kOutSIM5and10MergedBase =
        "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/photonJet5and10merged_SIM";

  inline const string kOutSIM5and20MergedBase =
        "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/photonJet5and20merged_SIM";

  inline const string kOutSIM5and10and20MergedBase =
        "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/photonJet5and10and20merged_SIM";

  // Merged SIM ROOT outputs (weighted merges)
  inline const string kMergedSIMOut =
        "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/photonJet10and20merged_SIM/RecoilJets_photonjet10plus20_MERGED.root";

  inline const string kMergedSIMOut_5and10 =
        "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/photonJet5and10merged_SIM/RecoilJets_photonjet5plus10_MERGED.root";

  inline const string kMergedSIMOut_5and20 =
        "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/photonJet5and20merged_SIM/RecoilJets_photonjet5plus20_MERGED.root";

  inline const string kMergedSIMOut_5and10and20 =
        "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/photonJet5and10and20merged_SIM/RecoilJets_photonjet5plus10plus20_MERGED.root";

  // Cross sections (pb) used for per-event weights w = sigma / Naccepted
  inline constexpr double kSigmaPhoton5_pb  = 89266.571;
  inline constexpr double kSigmaPhoton10_pb = 6692.7611;
  inline constexpr double kSigmaPhoton20_pb = 105.79868;

  // =============================================================================
  // Binning
  // =============================================================================
  inline constexpr int kNPtBins = 9;
  inline constexpr int kPtEdges[kNPtBins + 1] = {10,12,14,16,18,20,22,24,26,35};

  // Jet radii keys
  inline const vector<string> kRKeys = {"r02","r04"};

  // =============================================================================
  // ANSI helpers
  // =============================================================================
  inline const string ANSI_RESET    = "\033[0m";
  inline const string ANSI_BOLD_RED = "\033[1;31m";
  inline const string ANSI_BOLD_GRN = "\033[1;32m";
  inline const string ANSI_BOLD_YEL = "\033[1;33m";
  inline const string ANSI_BOLD_CYN = "\033[1;36m";
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
    v.reserve(kNPtBins);

    for (int i = 0; i < kNPtBins; ++i)
    {
      PtBin b;
      b.lo = kPtEdges[i];
      b.hi = kPtEdges[i+1];
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

  inline double RFromKey(const string& rKey)
  {
    if (rKey == "r02") return 0.2;
    if (rKey == "r04") return 0.4;
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
  // Robust ROOT getters (template must live in header)
  // =============================================================================
  template <class T>
  inline T* GetObj(Dataset& ds, const string& relName,
                   bool logMissing = true,
                   bool logZero = true,
                   bool treatZeroAsMissing = false)
  {
    if (!ds.topDir) return nullptr;

    // PP-only pass: ignore centrality-tagged objects if caller accidentally asks.
    if (relName.find("_cent_") != string::npos) return nullptr;

    const string fp = FullPath(ds, relName);
    ds.requestCounts[fp]++;

    TObject* obj = ds.topDir->Get(relName.c_str());

    if (!obj)
    {
      if (logMissing) LogMissing(ds, fp, "MISSING");
      return nullptr;
    }

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
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadTopMargin(0.08);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
  }

  inline void ApplyCanvasMargins1D(TCanvas& c)
  {
    c.SetLeftMargin(0.12);
    c.SetRightMargin(0.05);
    c.SetBottomMargin(0.12);
    c.SetTopMargin(0.08);
  }

  inline void ApplyCanvasMargins2D(TCanvas& c)
  {
    c.SetLeftMargin(0.12);
    c.SetRightMargin(0.15);
    c.SetBottomMargin(0.12);
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

  inline void DrawIsoCornerLabels1D(const string& isoTitle, const PtBin& pb)
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
      DrawIsoCornerLabels1D(isoTitle, pb);

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

      if (name.find("_cent_") != string::npos) continue;

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

    TCanvas c("c_tbl","c_tbl",1500,1200);
    c.Divide(3,3, 0.001, 0.001);

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
        cout << ANSI_BOLD_RED
             << "[MERGE SIM][FATAL] BuildMergedSIMFile_PhotonSlices needs >=2 inputs and matching sigma list.\n"
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
          cout << ANSI_BOLD_RED << "[MERGE SIM][FATAL] Cannot open input SIM file: " << inFiles[i] << ANSI_RESET << "\n";
          CloseAll();
          return false;
        }

        din[i] = fin[i]->GetDirectory(topDirName.c_str());
        if (!din[i])
        {
          cout << ANSI_BOLD_RED << "[MERGE SIM][FATAL] Missing topDir in SIM file: " << inFiles[i] << ANSI_RESET << "\n";
          CloseAll();
          return false;
        }

        N[i] = ReadEventCountFromFile(fin[i], topDirName);
        if (N[i] <= 0.0)
        {
          cout << ANSI_BOLD_RED << "[MERGE SIM][FATAL] Naccepted <= 0 for input: " << inFiles[i] << ANSI_RESET << "\n";
          CloseAll();
          return false;
        }

        w[i] = sigmas_pb[i] / N[i];
      }

      cout << ANSI_BOLD_YEL << "[MERGE SIM] Slice weights (w = sigma/N) [pb/event]:\n" << ANSI_RESET;
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
        cout << ANSI_BOLD_RED << "[MERGE SIM][FATAL] Cannot create merged output file." << ANSI_RESET << "\n";
        CloseAll();
        return false;
      }

      fout->cd();
      TDirectory* outTop = fout->mkdir(topDirName.c_str());
      if (!outTop)
      {
        cout << ANSI_BOLD_RED << "[MERGE SIM][FATAL] Cannot create topDir in merged output." << ANSI_RESET << "\n";
        fout->Close();
        CloseAll();
        return false;
      }

      // Build output by adding each slice into the same output tree
      for (size_t i = 0; i < n; ++i)
      {
        AddScaledRecursive(outTop, din[i], w[i]);
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
           << "[MERGE SIM] NOTE: histograms are weighted (units ~ pb per bin). Entries are no longer raw event counts.\n"
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
          (allPhoton5and10and20sim ? 1 : 0);

        if (nTrue == 0) return SimSample::kNone;
        if (nTrue != 1) return SimSample::kInvalid;

        if (isPhotonJet5)            return SimSample::kPhotonJet5;
        if (isPhotonJet10)           return SimSample::kPhotonJet10;
        if (isPhotonJet20)           return SimSample::kPhotonJet20;
        if (bothPhoton5and10sim)     return SimSample::kPhotonJet5And10Merged;
        if (bothPhoton5and20sim)     return SimSample::kPhotonJet5And20Merged;
        if (bothPhoton10and20sim)    return SimSample::kPhotonJet10And20Merged;
        if (allPhoton5and10and20sim) return SimSample::kPhotonJet5And10And20Merged;

        return SimSample::kInvalid;
  }

  inline string SimSampleLabel(SimSample s)
  {
        switch (s)
        {
          case SimSample::kNone:                     return "NONE";
          case SimSample::kPhotonJet5:               return "photonJet5";
          case SimSample::kPhotonJet10:              return "photonJet10";
          case SimSample::kPhotonJet20:              return "photonJet20";
          case SimSample::kPhotonJet5And10Merged:    return "photonJet5and10merged";
          case SimSample::kPhotonJet5And20Merged:    return "photonJet5and20merged";
          case SimSample::kPhotonJet10And20Merged:   return "photonJet10and20merged";
          case SimSample::kPhotonJet5And10And20Merged:return "photonJet5and10and20merged";
          default:                                   return "INVALID";
        }
  }

  inline string SimInputPathForSample(SimSample s)
  {
        switch (s)
        {
          case SimSample::kPhotonJet5:                return kInSIM5;
          case SimSample::kPhotonJet10:               return kInSIM10;
          case SimSample::kPhotonJet20:               return kInSIM20;
          case SimSample::kPhotonJet5And10Merged:     return kMergedSIMOut_5and10;
          case SimSample::kPhotonJet5And20Merged:     return kMergedSIMOut_5and20;
          case SimSample::kPhotonJet10And20Merged:    return kMergedSIMOut;
          case SimSample::kPhotonJet5And10And20Merged:return kMergedSIMOut_5and10and20;
          default:                                    return "";
        }
  }

  inline string SimOutBaseForSample(SimSample s)
  {
        switch (s)
        {
          case SimSample::kPhotonJet5:                return kOutSIM5Base;
          case SimSample::kPhotonJet10:               return kOutSIM10Base;
          case SimSample::kPhotonJet20:               return kOutSIM20Base;
          case SimSample::kPhotonJet5And10Merged:     return kOutSIM5and10MergedBase;
          case SimSample::kPhotonJet5And20Merged:     return kOutSIM5and20MergedBase;
          case SimSample::kPhotonJet10And20Merged:    return kOutSIMMergedBase;
          case SimSample::kPhotonJet5And10And20Merged:return kOutSIM5and10and20MergedBase;
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

        const SimSample ss = CurrentSimSample();

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
                "allPhoton5and10and20sim=false.";
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
              "bothPhoton5and10sim, bothPhoton5and20sim, bothPhoton10and20sim, allPhoton5and10and20sim.";
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
              "bothPhoton5and10sim, bothPhoton5and20sim, bothPhoton10and20sim, allPhoton5and10and20sim.";
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
      kSimOnly,
      kSimAndDataPP,
      kInvalid
  };

  inline RunMode CurrentRunMode()
  {
      if (!ValidateRunConfig(nullptr)) return RunMode::kInvalid;
      if (isPPdataOnly)   return RunMode::kPPDataOnly;
      if (isSimAndDataPP) return RunMode::kSimAndDataPP;
      return RunMode::kSimOnly;
  }

  inline string RunModeLabel(RunMode m)
  {
      switch (m)
      {
        case RunMode::kPPDataOnly:   return "PP_DATA_ONLY";
        case RunMode::kSimOnly:      return "SIM_ONLY";
        case RunMode::kSimAndDataPP: return "SIM_AND_DATA_PP";
        default:                     return "INVALID";
      }
  }


} // namespace ARJ

#endif // ANALYZE_RECOIL_JETS_H
