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
  // GLOBAL USER TOGGLES (EDIT HERE)
  // =============================================================================
  inline bool isPPdataOnly               = false;
  inline bool isSimOnly                  = true;
  inline bool isSimAndDataPP             = false;
  inline bool isSimOnlyWithPhoton10and20 = false;

  // Displayed range [-vzCutCm,+vzCutCm] and 0.5 cm display bin width
  inline double vzCutCm = 30.0;

  // =============================================================================
  // FIXED INPUTS (pp + photonJet10 SIM only)
  // =============================================================================
  inline const string kTriggerPP = "Photon_4_GeV_plus_MBD_NS_geq_1";
  inline const string kDirSIM    = "SIM";

  inline const string kInPP =
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/pp/RecoilJets_pp_ALL.root";

  inline const string kInSIM10 =
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/photonJet10_SIM/RecoilJets_photonjet10_ALL.root";

  inline const string kInSIM20 =
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/photonJet20_SIM/RecoilJets_photonjet20_ALL.root";

  inline const string kOutPPBase =
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/pp";

  inline const string kOutSIMBase =
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/photonJet10_SIM";

  inline const string kMergedSIMOut =
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/photonJet10_SIM/RecoilJets_photonjet10plus20_MERGED.root";

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
    if (ds.isSim) lines.push_back("Dataset: SIM (photonJet10)");
    else          lines.push_back(string("Dataset: DATA (") + ds.trigger + ")");
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

    const string yTitleEff =
      (yTitle == "A.U." || yTitle == "A.U")
        ? "Fraction of entries"
        : yTitle;

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

      TH1* h = GetObj<TH1>(ds, hname, true, true, true);
      if (!h)
      {
        TLatex t;
        t.SetNDC(true);
        t.SetTextFont(42);
        t.SetTextSize(0.06);
        t.DrawLatex(0.15, 0.55, "MISSING");
        t.SetTextSize(0.05);
        std::ostringstream s;
        s << histBase << "  pT^{#gamma}: " << b.lo << "-" << b.hi;
        t.DrawLatex(0.15, 0.45, s.str().c_str());
        continue;
      }

      TH1* hc = CloneTH1(h, TString::Format("%s_tbl_%d", histBase.c_str(), i).Data());
      if (!hc) continue;
      if (normalizeShape) NormalizeToUnitArea(hc);

      hc->SetLineWidth(2);
      hc->SetTitle("");
      hc->GetXaxis()->SetTitle(xTitle.c_str());
      hc->GetYaxis()->SetTitle(yTitle.c_str());
      hc->Draw("hist");

      vector<string> lines = commonLines;
      {
        std::ostringstream s;
        s << histBase << "   pT^{#gamma}: " << b.lo << "-" << b.hi << " GeV";
        lines.push_back(s.str());
      }
      DrawLatexLines(0.16, 0.90, lines, 0.040, 0.050);

      delete hc;
    }

    SaveCanvas(c, JoinPath(outDir, outName));
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
  // SIM merge utilities (photonJet10 + photonJet20)
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

  inline bool BuildMergedSIMFile_Photon10And20(const string& in10,
                                               const string& in20,
                                               const string& outMerged,
                                               const string& topDirName,
                                               double sigma10_pb,
                                               double sigma20_pb)
  {
    cout << ANSI_BOLD_CYN
         << "\n[MERGE SIM] Building merged SIM file with cross-section weights\n"
         << "  in10 = " << in10 << "\n"
         << "  in20 = " << in20 << "\n"
         << "  out  = " << outMerged << "\n"
         << "  topDir = " << topDirName << "\n"
         << ANSI_RESET;

    TFile* f10 = TFile::Open(in10.c_str(), "READ");
    TFile* f20 = TFile::Open(in20.c_str(), "READ");
    if (!f10 || f10->IsZombie() || !f20 || f20->IsZombie())
    {
      cout << ANSI_BOLD_RED << "[MERGE SIM][FATAL] Cannot open one of the input SIM files." << ANSI_RESET << "\n";
      if (f10) f10->Close();
      if (f20) f20->Close();
      return false;
    }

    TDirectory* d10 = f10->GetDirectory(topDirName.c_str());
    TDirectory* d20 = f20->GetDirectory(topDirName.c_str());
    if (!d10 || !d20)
    {
      cout << ANSI_BOLD_RED << "[MERGE SIM][FATAL] Missing topDir in one of the SIM files." << ANSI_RESET << "\n";
      f10->Close(); f20->Close();
      return false;
    }

    const double N10 = ReadEventCountFromFile(f10, topDirName);
    const double N20 = ReadEventCountFromFile(f20, topDirName);

    cout << ANSI_BOLD_YEL
         << "[MERGE SIM] Event counts from cnt_" << topDirName << " bin1:\n"
         << "  N10 = " << std::fixed << std::setprecision(0) << N10 << "\n"
         << "  N20 = " << std::fixed << std::setprecision(0) << N20 << "\n"
         << ANSI_RESET;

    if (N10 <= 0.0 || N20 <= 0.0)
    {
      cout << ANSI_BOLD_RED << "[MERGE SIM][FATAL] N10 or N20 <= 0 (cannot compute weights)." << ANSI_RESET << "\n";
      f10->Close(); f20->Close();
      return false;
    }

    const double w10 = sigma10_pb / N10;
    const double w20 = sigma20_pb / N20;

    const double sigmaTot = sigma10_pb + sigma20_pb;
    const double f10frac  = (sigmaTot > 0.0) ? (sigma10_pb / sigmaTot) : 0.0;
    const double f20frac  = (sigmaTot > 0.0) ? (sigma20_pb / sigmaTot) : 0.0;

    const double Leff = (sigma10_pb > 0.0 ? (N10 / sigma10_pb) : 0.0)
                      + (sigma20_pb > 0.0 ? (N20 / sigma20_pb) : 0.0);

    cout << ANSI_BOLD_YEL
         << "[MERGE SIM] Cross sections (pb): sigma10=" << sigma10_pb
         << "  sigma20=" << sigma20_pb
         << "  sigmaTot=" << sigmaTot << "\n"
         << "[MERGE SIM] Slice fractions:     f10=" << std::setprecision(6) << f10frac
         << "  f20=" << std::setprecision(6) << f20frac << "\n"
         << "[MERGE SIM] Per-event weights:   w10=" << std::setprecision(12) << w10
         << "  w20=" << std::setprecision(12) << w20 << "   [pb/event]\n"
         << "[MERGE SIM] Effective Lumi:      Leff=" << std::setprecision(6) << Leff << "   [pb^{-1}]\n"
         << ANSI_RESET;

    EnsureParentDirForFile(outMerged);
    TFile* fout = TFile::Open(outMerged.c_str(), "RECREATE");
    if (!fout || fout->IsZombie())
    {
      cout << ANSI_BOLD_RED << "[MERGE SIM][FATAL] Cannot create merged output file." << ANSI_RESET << "\n";
      f10->Close(); f20->Close();
      return false;
    }

    fout->cd();
    TDirectory* outTop = fout->mkdir(topDirName.c_str());
    if (!outTop)
    {
      cout << ANSI_BOLD_RED << "[MERGE SIM][FATAL] Cannot create topDir in merged output." << ANSI_RESET << "\n";
      fout->Close(); f10->Close(); f20->Close();
      return false;
    }

    CopyAndScaleAddRecursive(outTop, d10, w10, d20, w20);

    outTop->cd();
    TNamed meta("MERGE_INFO",
      TString::Format("Merged photonJet10+20. N10=%.0f N20=%.0f sigma10=%.6f pb sigma20=%.6f pb w10=%.12g w20=%.12g",
                      N10, N20, sigma10_pb, sigma20_pb, w10, w20).Data());
    meta.Write("MERGE_INFO", TObject::kOverwrite);

    fout->Write();
    fout->Close();
    f10->Close();
    f20->Close();

    cout << ANSI_BOLD_CYN
         << "[MERGE SIM] Done. Merged file written: " << outMerged << "\n"
         << "[MERGE SIM] NOTE: histograms are now weighted (units ~ pb per bin). Entries are no longer raw event counts.\n"
         << ANSI_RESET;

    return true;
  }

  // =============================================================================
  // Run-mode guard (moved out of cpp)
  // =============================================================================
  inline bool ExactlyOneModeSet()
  {
    const int nTrue =
      (isPPdataOnly ? 1 : 0) +
      (isSimOnly ? 1 : 0) +
      (isSimAndDataPP ? 1 : 0) +
      (isSimOnlyWithPhoton10and20 ? 1 : 0);
    return (nTrue == 1);
  }


  // =============================================================================
  // Run mode helpers
  // =============================================================================
  enum class RunMode
  {
    kPPDataOnly,
    kSimOnly,
    kSimAndDataPP,
    kSimOnlyMergedPhoton10And20,
    kInvalid
  };

  inline RunMode CurrentRunMode()
  {
    if (!ExactlyOneModeSet()) return RunMode::kInvalid;
    if (isPPdataOnly)               return RunMode::kPPDataOnly;
    if (isSimOnly)                  return RunMode::kSimOnly;
    if (isSimAndDataPP)             return RunMode::kSimAndDataPP;
    if (isSimOnlyWithPhoton10and20) return RunMode::kSimOnlyMergedPhoton10And20;
    return RunMode::kInvalid;
  }

  inline string RunModeLabel(RunMode m)
  {
    switch (m)
    {
      case RunMode::kPPDataOnly:               return "PP_DATA_ONLY";
      case RunMode::kSimOnly:                  return "SIM_ONLY";
      case RunMode::kSimAndDataPP:             return "SIM_AND_DATA_PP";
      case RunMode::kSimOnlyMergedPhoton10And20:return "SIM_ONLY_MERGED_PHOTON10AND20";
      default:                                 return "INVALID";
    }
  }


} // namespace ARJ

#endif // ANALYZE_RECOIL_JETS_H
