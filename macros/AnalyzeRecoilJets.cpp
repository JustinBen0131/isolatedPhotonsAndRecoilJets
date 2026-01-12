// AnalyzeRecoilJets.cpp
// Run with:
//   root -b -q -l AnalyzeRecoilJets.cpp
//
// Self-contained offline reader/plotter for RecoilJets outputs (PP + photonJet10 SIM only).
// Auto-executes at load and exits with nonzero on failure.

#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TObject.h>

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

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <limits>

namespace ARJ
{
  using std::cout;
  using std::cerr;
  using std::endl;
  using std::string;
  using std::vector;
  using std::map;

  // =============================================================================
  // GLOBAL USER TOGGLES (EDIT HERE)
  // =============================================================================
  bool isPPdataOnly   = false;
  bool isSimOnly      = true;
  bool isSimAndDataPP = false;

  double vzCutCm = 30.0;  // displayed range [-vzCutCm,+vzCutCm] and 0.5 cm display bin width

  // =============================================================================
  // FIXED INPUTS (pp + photonJet10 SIM only)
  // =============================================================================
  static const string kTriggerPP = "Photon_4_GeV_plus_MBD_NS_geq_1";
  static const string kDirSIM    = "SIM";

  static const string kInPP  = "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/pp/RecoilJets_pp_ALL.root";
  static const string kInSIM = "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/photonJet10_SIM/RecoilJets_photonjet10_ALL.root";

  static const string kOutPPBase  = "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/pp";
  static const string kOutSIMBase = "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/photonJet10_SIM";

  // Photon pT bin edges used for *_pT_<lo>_<hi> suffixes
  static const int kNPtBins = 9;
  static const int kPtEdges[kNPtBins + 1] = {10,12,14,16,18,20,22,24,26,35};

  // Jet radii
  static const vector<string> kRKeys = {"r02","r04"};

  // ANSI helpers
  static const string ANSI_RESET    = "\033[0m";
  static const string ANSI_BOLD_RED = "\033[1;31m";
  static const string ANSI_BOLD_CYN = "\033[1;36m";
  static const string ANSI_BOLD_YEL = "\033[1;33m";
  static const string ANSI_DIM      = "\033[2m";

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

  static const vector<PtBin>& PtBins()
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

  static double RFromKey(const string& rKey)
  {
    if (rKey == "r02") return 0.2;
    if (rKey == "r04") return 0.4;
    return 0.0;
  }

  static double FidEtaAbsFromKey(const string& rKey)
  {
    const double R = RFromKey(rKey);
    return 1.1 - R;
  }

  // REPLACE WITH THIS (makes it explicit: non-copyable, but movable; avoids future STL surprises)
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

      std::ofstream missingOut;
      int missingCount = 0;

      Dataset() = default;
      Dataset(const Dataset&) = delete;
      Dataset& operator=(const Dataset&) = delete;
      Dataset(Dataset&&) noexcept = default;
      Dataset& operator=(Dataset&&) noexcept = default;
  };

  struct MatchCache
  {
    // per rKey: vectors size kNPtBins
    map<string, vector<double> > NphoLead;
    map<string, vector<double> > NphoMatched;
  };

  struct LeakageFactors
  {
    vector<double> fB, fC, fD; // per pT bin (kNPtBins)
    bool available = false;
  };

  // =============================================================================
  // Path helpers
  // =============================================================================
  static void EnsureDir(const string& path)
  {
    if (path.empty()) return;
    gSystem->mkdir(path.c_str(), true);
  }

  static string DirnameFromPath(const string& filepath)
  {
    const size_t pos = filepath.find_last_of('/');
    if (pos == string::npos) return "";
    return filepath.substr(0, pos);
  }

  static void EnsureParentDirForFile(const string& filepath)
  {
    EnsureDir(DirnameFromPath(filepath));
  }

  static string JoinPath(const string& a, const string& b)
  {
    if (a.empty()) return b;
    if (b.empty()) return a;
    if (a.back() == '/') return a + b;
    return a + "/" + b;
  }

  static string FullPath(const Dataset& ds, const string& rel)
  {
    return ds.topDirName + "/" + rel;
  }

  static void LogMissing(Dataset& ds, const string& fullpath, const string& reason)
  {
    if (ds.missingOut.is_open())
    {
      ds.missingOut << fullpath << " :: " << reason << "\n";
    }
    ds.missingCount++;
  }

  // =============================================================================
  // ROOT getters (robust)
  // =============================================================================
  template <class T>
  static T* GetObj(Dataset& ds, const string& relName,
                   bool logMissing = true,
                   bool logZero = true,
                   bool treatZeroAsMissing = false)
  {
    if (!ds.topDir) return nullptr;

    // PP-only pass: ignore centrality-tagged objects if caller accidentally asks.
    if (relName.find("_cent_") != string::npos) return nullptr;

    TObject* obj = ds.topDir->Get(relName.c_str());
    const string fp = FullPath(ds, relName);

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

  static double SafeDivide(double a, double b, double def = 0.0)
  {
    if (b == 0.0) return def;
    return a / b;
  }

  static TH1* CloneTH1(const TH1* h, const string& newName)
  {
    if (!h) return nullptr;
    TH1* c = dynamic_cast<TH1*>(h->Clone(newName.c_str()));
    if (!c) return nullptr;
    c->SetDirectory(nullptr);
    return c;
  }

  static TH2* CloneTH2(const TH2* h, const string& newName)
  {
    if (!h) return nullptr;
    TH2* c = dynamic_cast<TH2*>(h->Clone(newName.c_str()));
    if (!c) return nullptr;
    c->SetDirectory(nullptr);
    return c;
  }

  static TH3* CloneTH3(const TH3* h, const string& newName)
  {
    if (!h) return nullptr;
    TH3* c = dynamic_cast<TH3*>(h->Clone(newName.c_str()));
    if (!c) return nullptr;
    c->SetDirectory(nullptr);
    return c;
  }

  // =============================================================================
  // Style helpers
  // =============================================================================
  static void SetupGlobalStyle()
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

  static void ApplyCanvasMargins1D(TCanvas& c)
  {
    c.SetLeftMargin(0.12);
    c.SetRightMargin(0.05);
    c.SetBottomMargin(0.12);
    c.SetTopMargin(0.08);
  }

  static void ApplyCanvasMargins2D(TCanvas& c)
  {
    c.SetLeftMargin(0.12);
    c.SetRightMargin(0.15);
    c.SetBottomMargin(0.12);
    c.SetTopMargin(0.08);
  }

  static void DrawLatexLines(double x, double y,
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

  static vector<string> DefaultHeaderLines(const Dataset& ds)
  {
    vector<string> lines;
    if (ds.isSim) lines.push_back("Dataset: SIM (photonJet10)");
    else          lines.push_back(string("Dataset: DATA (") + ds.trigger + ")");
    return lines;
  }

  static void SaveCanvas(TCanvas& c, const string& filepath)
  {
    EnsureParentDirForFile(filepath);
    c.SaveAs(filepath.c_str());
  }

  static void NormalizeToUnitArea(TH1* h)
  {
    if (!h) return;
    const double integral = h->Integral(0, h->GetNbinsX() + 1);
    if (integral > 0.0) h->Scale(1.0 / integral);
  }

  static void DrawFidEtaLines1D(double etaAbs)
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

  static void DrawAndSaveTH1_Common(const Dataset& ds, TH1* h,
                                   const string& filepath,
                                   const string& xTitle,
                                   const string& yTitle,
                                   const vector<string>& extraLines,
                                   bool logy,
                                   bool drawFidLines = false,
                                   double fidEtaAbs  = 0.0)
  {
    if (!h) return;

    TCanvas c("c1","c1",900,700);
    ApplyCanvasMargins1D(c);
    c.SetLogy(logy);

    h->SetLineWidth(2);
    h->SetTitle("");
    h->GetXaxis()->SetTitle(xTitle.c_str());
    h->GetYaxis()->SetTitle(yTitle.c_str());

    h->Draw("hist");

    if (drawFidLines && fidEtaAbs > 0.0)
    {
      DrawFidEtaLines1D(fidEtaAbs);
    }

    vector<string> lines = DefaultHeaderLines(ds);
    for (const auto& s : extraLines) lines.push_back(s);
    DrawLatexLines(0.14, 0.92, lines, 0.034, 0.045);

    SaveCanvas(c, filepath);
  }

  static void DrawAndSaveTH2_Common(const Dataset& ds, TH2* h,
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
      // assume X is eta; draw vertical lines at ±fidEtaAbs over Y range
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

  // =============================================================================
  // Vertex Z rebinner to enforce displayed 0.5 cm bins over [-vzCut,+vzCut]
  // =============================================================================
  static TH1F* RebinToFixedBinWidthVertexZ(const TH1* hOrig, double vzCut)
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
  struct InvItem
  {
    string path;
    string cls;
    double entries;
  };

  static void CollectInventoryRecursive(TDirectory* dir, const string& prefix, vector<InvItem>& out)
  {
    if (!dir) return;
    TIter next(dir->GetListOfKeys());
    while (TKey* key = (TKey*)next())
    {
      const string name = key->GetName();
      const string cls  = key->GetClassName();

      // PP-only ignore centrality-suffixed objects
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

  static void PrintInventoryToTerminal(const Dataset& ds, const vector<InvItem>& items)
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
  static void Make3x3Table_TH1(Dataset& ds,
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
  // Leakage corrected solver:
  // Solve fixed-point: S = A - (B - fB*S)(C - fC*S)/(D - fD*S)
  // Returns true if converged.
  // =============================================================================
  static bool SolveLeakageCorrectedSA(double A, double B, double C, double D,
                                     double fB, double fC, double fD,
                                     double& outSA)
  {
    outSA = 0.0;
    if (A <= 0.0) { outSA = 0.0; return true; }

    // start from raw estimate if possible
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
      // keep denom away from zero if possible
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
  // Section 1: event counts + inventory + zvtx
  // =============================================================================
  static double ReadEventCount(Dataset& ds)
  {
    const string cntName = "cnt_" + ds.topDirName;
    TH1* cnt = GetObj<TH1>(ds, cntName, true, true, false);
    if (!cnt) return 0.0;
    return cnt->GetBinContent(1);
  }

  static void Section1_EventLevel(Dataset& ds)
  {
    cout << ANSI_BOLD_CYN << "\n==============================\n"
         << "[SECTION 1] EVENT-LEVEL + INVENTORY (" << ds.label << ")\n"
         << "==============================" << ANSI_RESET << "\n";

    const double nEvt = ReadEventCount(ds);
    cout << ANSI_BOLD_RED
         << "[" << ds.label << "] accepted events = " << std::fixed << std::setprecision(0) << nEvt
         << ANSI_RESET << "\n";

    vector<InvItem> items;
    CollectInventoryRecursive(ds.topDir, ds.topDirName + "/", items);
    std::sort(items.begin(), items.end(), [](const InvItem& a, const InvItem& b){ return a.path < b.path; });
    PrintInventoryToTerminal(ds, items);

    TH1* hV = GetObj<TH1>(ds, "h_vertexZ", true, true, true);
    if (!hV)
    {
      cout << ANSI_BOLD_YEL << "[WARN] Missing h_vertexZ for " << ds.label << ANSI_RESET << "\n";
      return;
    }

    TH1F* hFixed = RebinToFixedBinWidthVertexZ(hV, vzCutCm);
    if (!hFixed)
    {
      cout << ANSI_BOLD_YEL << "[WARN] Failed to build fixed-binwidth vertexZ for " << ds.label << ANSI_RESET << "\n";
      return;
    }

    string outPath;
    if (ds.isSim)
      outPath = JoinPath(ds.outBase, "GeneralEventLevelQA/zvtx_SIM.png");
    else
      outPath = JoinPath(ds.outBase, "GeneralEventLevelQA/" + ds.trigger + "/zvtx_DATA_" + ds.trigger + ".png");

    vector<string> lines;
    lines.push_back(TString::Format("|v_{z}| < %.0f cm", std::fabs(vzCutCm)).Data());
    lines.push_back("Displayed bin width = 0.5 cm");
    DrawAndSaveTH1_Common(ds, hFixed, outPath, "v_{z} [cm]", "Counts", lines, false);

    delete hFixed;
  }

  // =============================================================================
  // Section 2: preselection fail counters (terminal only)
  // =============================================================================
  static double Read1BinCount(Dataset& ds, const string& hname)
  {
    TH1* h = GetObj<TH1>(ds, hname, true, true, false);
    if (!h) return 0.0;
    return h->GetBinContent(1);
  }

  static void PrintPreselectionFailTable(Dataset& ds)
  {
    cout << ANSI_BOLD_CYN << "\n==============================\n"
         << "[SECTION 2] PRESELECTION FAIL QA (TERMINAL ONLY) (" << ds.label << ")\n"
         << "==============================" << ANSI_RESET << "\n";

    const int wBin = 10;
    const int wN   = 12;

    cout << std::left
         << std::setw(wBin) << "pTbin"
         << std::right
         << std::setw(wN) << "wetaFail"
         << std::setw(wN) << "et1Low"
         << std::setw(wN) << "et1High"
         << std::setw(wN) << "et1Out"
         << std::setw(wN) << "e11e33Hi"
         << std::setw(wN) << "e32e35Lo"
         << std::setw(wN) << "e32e35Hi"
         << std::setw(wN) << "e32e35Out"
         << "\n";
    cout << string(wBin + 8*wN, '-') << "\n";

    const auto& bins = PtBins();
    for (int i = 0; i < kNPtBins; ++i)
    {
      const PtBin& b = bins[i];
      const string& suf = b.suffix;

      const double weta = Read1BinCount(ds, "h_preFail_weta" + suf);
      const double et1L = Read1BinCount(ds, "h_preFail_et1_low" + suf);
      const double et1H = Read1BinCount(ds, "h_preFail_et1_high" + suf);
      const double et1O = et1L + et1H;

      const double e11H = Read1BinCount(ds, "h_preFail_e11e33_high" + suf);
      const double e32L = Read1BinCount(ds, "h_preFail_e32e35_low" + suf);
      const double e32H = Read1BinCount(ds, "h_preFail_e32e35_high" + suf);
      const double e32O = e32L + e32H;

      cout << std::left << std::setw(wBin) << b.label
           << std::right
           << std::setw(wN) << std::fixed << std::setprecision(0) << weta
           << std::setw(wN) << et1L
           << std::setw(wN) << et1H
           << std::setw(wN) << et1O
           << std::setw(wN) << e11H
           << std::setw(wN) << e32L
           << std::setw(wN) << e32H
           << std::setw(wN) << e32O
           << "\n";
    }

    cout << ANSI_DIM
         << "\nNOTE: Preselection fail counters are inclusive (one photon can increment multiple fail histograms).\n"
         << ANSI_RESET;
  }

  // =============================================================================
  // Section 3: general isolation QA
  // =============================================================================
  static void PrintIsoDecisionTable(Dataset& ds)
  {
    cout << ANSI_BOLD_CYN << "\n[SECTION 3] isoDecision table (" << ds.label << ")\n" << ANSI_RESET;

    const int wBin = 10;
    const int wN   = 14;

    cout << std::left
         << std::setw(wBin) << "pTbin"
         << std::right
         << std::setw(wN) << "N_iso"
         << std::setw(wN) << "N_nonIso"
         << std::setw(wN) << "isoFrac"
         << "\n";
    cout << string(wBin + 3*wN, '-') << "\n";

    const auto& bins = PtBins();
    for (int i = 0; i < kNPtBins; ++i)
    {
      const PtBin& b = bins[i];
      const string hname = "h_isoDecision" + b.suffix;

      TH1* h = GetObj<TH1>(ds, hname, true, true, false);
      double nIso = 0.0, nNon = 0.0;
      if (h)
      {
        nIso = h->GetBinContent(1);
        nNon = h->GetBinContent(2);
      }
      const double denom = nIso + nNon;
      const double frac  = (denom > 0.0) ? (nIso/denom) : 0.0;

      cout << std::left << std::setw(wBin) << b.label
           << std::right
           << std::setw(wN) << std::fixed << std::setprecision(0) << nIso
           << std::setw(wN) << nNon
           << std::setw(wN) << std::fixed << std::setprecision(4) << frac
           << "\n";
    }
  }

  static void Section3_GeneralIsoQA(Dataset& ds)
  {
    cout << ANSI_BOLD_CYN << "\n==============================\n"
         << "[SECTION 3] GENERAL ISOLATION QA (" << ds.label << ")\n"
         << "==============================" << ANSI_RESET << "\n";

    string outDir;
    if (ds.isSim) outDir = JoinPath(ds.outBase, "isoQAgeneral");
    else          outDir = JoinPath(ds.outBase, "isoQAgeneral/" + ds.trigger);

    EnsureDir(outDir);
    for (const auto& b : PtBins()) EnsureDir(JoinPath(outDir, b.folder));

    vector<string> common;
    common.push_back("Cuts: p_{T}^{#gamma} #geq 5 GeV, |#eta^{#gamma}| < 0.7");

    struct IsoHistDef { string base; string outStem; };
    const vector<IsoHistDef> isoHists = {
      {"h_Eiso",        "Eiso_total"},
      {"h_Eiso_emcal",  "Eiso_emcal"},
      {"h_Eiso_hcalin", "Eiso_hcalin"},
      {"h_Eiso_hcalout","Eiso_hcalout"}
    };

    for (const auto& def : isoHists)
    {
      Make3x3Table_TH1(ds, def.base, outDir,
                       string("table3x3_") + def.outStem + ".png",
                       "E_{iso} [GeV]", "Counts",
                       false, false, common);

      for (int i = 0; i < kNPtBins; ++i)
      {
        const PtBin& b = PtBins()[i];
        const string hname = def.base + b.suffix;

        TH1* h = GetObj<TH1>(ds, hname, true, true, true);
        if (!h) continue;

        TH1* hc = CloneTH1(h, TString::Format("%s_%d", def.base.c_str(), i).Data());
        if (!hc) continue;

        vector<string> lines = common;
        lines.push_back(TString::Format("p_{T}^{#gamma} bin: %d-%d GeV", b.lo, b.hi).Data());

        const string fp = JoinPath(outDir, b.folder + "/" + def.outStem + "_" + b.folder + ".png");
        DrawAndSaveTH1_Common(ds, hc, fp, "E_{iso} [GeV]", "Counts", lines, false);
        delete hc;
      }
    }

    PrintIsoDecisionTable(ds);
  }

  // =============================================================================
  // Section 4: ABCD QA + purity + SS overlays + SIM leakage factors
  // =============================================================================
  static void ReadLeakageFactorsFromSIM(Dataset& sim, LeakageFactors& lf)
  {
    lf.fB.assign(kNPtBins, 0.0);
    lf.fC.assign(kNPtBins, 0.0);
    lf.fD.assign(kNPtBins, 0.0);
    lf.available = false;

    if (!sim.isSim) return;

    bool any = false;

    cout << ANSI_BOLD_CYN << "\n[SECTION 4] Reading SIM leakage factors (h_sigABCD_MC_pT_*)\n" << ANSI_RESET;

    const int wBin = 10;
    const int wN   = 14;

    cout << std::left << std::setw(wBin) << "pTbin"
         << std::right
         << std::setw(wN) << "A_sig_MC"
         << std::setw(wN) << "B_sig_MC"
         << std::setw(wN) << "C_sig_MC"
         << std::setw(wN) << "D_sig_MC"
         << std::setw(wN) << "fB"
         << std::setw(wN) << "fC"
         << std::setw(wN) << "fD"
         << "\n";
    cout << string(wBin + 7*wN, '-') << "\n";

    for (int i = 0; i < kNPtBins; ++i)
    {
      const PtBin& b = PtBins()[i];
      const string hname = "h_sigABCD_MC" + b.suffix;

      TH1* h = GetObj<TH1>(sim, hname, true, true, false);
      double A = 0, B = 0, C = 0, D = 0;
      if (h)
      {
        A = h->GetBinContent(1);
        B = h->GetBinContent(2);
        C = h->GetBinContent(3);
        D = h->GetBinContent(4);
      }
      const double fB = (A > 0.0) ? (B/A) : 0.0;
      const double fC = (A > 0.0) ? (C/A) : 0.0;
      const double fD = (A > 0.0) ? (D/A) : 0.0;

      lf.fB[i] = fB;
      lf.fC[i] = fC;
      lf.fD[i] = fD;

      if (A > 0.0) any = true;

      cout << std::left << std::setw(wBin) << b.label
           << std::right
           << std::setw(wN) << std::fixed << std::setprecision(0) << A
           << std::setw(wN) << B
           << std::setw(wN) << C
           << std::setw(wN) << D
           << std::setw(wN) << std::fixed << std::setprecision(6) << fB
           << std::setw(wN) << fC
           << std::setw(wN) << fD
           << "\n";
    }

    lf.available = any;
    if (!lf.available)
    {
      cout << ANSI_BOLD_YEL << "[WARN] No nonzero A_sig_MC found; leakage correction will be unavailable.\n" << ANSI_RESET;
    }
  }

  static void Make3x3Table_SSOverlay(Dataset& ds,
                                    const string& varKey,
                                    const string& outDir,
                                    const string& outName,
                                    const vector<string>& commonLines)
  {
    EnsureDir(outDir);

    const vector<string> regionTags = {
      "isIsolated_isTight",
      "notIsolated_isTight",
      "isIsolated_notTight",
      "notIsolated_notTight"
    };

    TCanvas c("c_ss_tbl","c_ss_tbl",1500,1200);
    c.Divide(3,3, 0.001, 0.001);

    const auto& bins = PtBins();

    for (int i = 0; i < kNPtBins; ++i)
    {
      c.cd(i+1);
      gPad->SetLeftMargin(0.14);
      gPad->SetRightMargin(0.05);
      gPad->SetBottomMargin(0.14);
      gPad->SetTopMargin(0.10);

      const PtBin& b = bins[i];

      vector<TH1*> hs(4, nullptr);
      vector<TH1*> hc(4, nullptr);

      for (int r = 0; r < 4; ++r)
      {
        const string hname = "h_ss_" + varKey + "_" + regionTags[r] + b.suffix;
        hs[r] = GetObj<TH1>(ds, hname, true, true, true);
        if (hs[r])
        {
          hc[r] = CloneTH1(hs[r], TString::Format("ss_%s_%d_%d", varKey.c_str(), i, r).Data());
          if (hc[r])
          {
            NormalizeToUnitArea(hc[r]);
            hc[r]->SetLineWidth(2);
            // colors: A=black, B=red, C=blue, D=magenta
            hc[r]->SetLineColor((r==0)?1:(r==1)?2:(r==2)?4:6);
          }
        }
      }

      TH1* first = nullptr;
      for (int r = 0; r < 4; ++r) if (hc[r]) { first = hc[r]; break; }

      if (!first)
      {
        TLatex t;
        t.SetNDC(true);
        t.SetTextFont(42);
        t.SetTextSize(0.06);
        t.DrawLatex(0.15, 0.55, "MISSING");
        t.SetTextSize(0.05);
        std::ostringstream s;
        s << "SS " << varKey << "  pT^{#gamma}: " << b.lo << "-" << b.hi;
        t.DrawLatex(0.15, 0.45, s.str().c_str());
        for (auto* p : hc) if (p) delete p;
        continue;
      }

      double ymax = 0.0;
      for (auto* p : hc) if (p) ymax = std::max(ymax, p->GetMaximum());
      first->SetMaximum(ymax * 1.25);

      first->SetTitle("");
      first->GetXaxis()->SetTitle(varKey.c_str());
      first->GetYaxis()->SetTitle("A.U.");
      first->Draw("hist");
      for (int r = 0; r < 4; ++r) if (hc[r] && hc[r] != first) hc[r]->Draw("hist same");

      vector<string> lines = commonLines;
      lines.push_back(TString::Format("SS %s  pT^{#gamma}: %d-%d GeV", varKey.c_str(), b.lo, b.hi).Data());
      DrawLatexLines(0.16, 0.90, lines, 0.040, 0.050);

      for (auto* p : hc) if (p) delete p;
    }

    SaveCanvas(c, JoinPath(outDir, outName));
  }

  static void Section4_ABCDPurityAndSS(Dataset& ds, const LeakageFactors& lf)
  {
    cout << ANSI_BOLD_CYN << "\n==============================\n"
         << "[SECTION 4] ABCD QA + PURITY + SS (" << ds.label << ")\n"
         << "==============================" << ANSI_RESET << "\n";

    string outDir;
    if (ds.isSim) outDir = JoinPath(ds.outBase, "PurityABCD");
    else          outDir = JoinPath(ds.outBase, "PurityABCD/" + ds.trigger);
    EnsureDir(outDir);

    vector<double> purityRaw(kNPtBins, 0.0);
    vector<double> purityCorr(kNPtBins, 0.0);
    vector<bool>   hasCorr(kNPtBins, false);

    cout << ANSI_BOLD_CYN << "\n[ABCD COUNTS + PURITY] " << ds.label << "\n" << ANSI_RESET;

    const int wBin = 10;
    const int wN   = 12;

    cout << std::left << std::setw(wBin) << "pTbin"
         << std::right
         << std::setw(wN) << "A"
         << std::setw(wN) << "B"
         << std::setw(wN) << "C"
         << std::setw(wN) << "D"
         << std::setw(wN) << "A_sig"
         << std::setw(wN) << "Pur_raw"
         << std::setw(wN) << "Pur_corr"
         << "\n";
    cout << string(wBin + 7*wN, '-') << "\n";

    for (int i = 0; i < kNPtBins; ++i)
    {
      const PtBin& b = PtBins()[i];
      const string& suf = b.suffix;

      const double A = Read1BinCount(ds, "h_isIsolated_isTight" + suf);
      const double B = Read1BinCount(ds, "h_notIsolated_isTight" + suf);
      const double C = Read1BinCount(ds, "h_isIsolated_notTight" + suf);
      const double D = Read1BinCount(ds, "h_notIsolated_notTight" + suf);

      double Asig = 0.0;
      double Praw = 0.0;
      if (A > 0.0 && D > 0.0)
      {
        Asig = A - B*(C/D);
        if (Asig < 0.0) Asig = 0.0;
        Praw = Asig / A;
      }

      purityRaw[i] = Praw;

      double Pc = Praw;
      bool okCorr = false;
      if (lf.available)
      {
        double SA = 0.0;
        const bool ok = SolveLeakageCorrectedSA(A,B,C,D, lf.fB[i], lf.fC[i], lf.fD[i], SA);
        if (ok && A > 0.0)
        {
          Pc = SA / A;
          okCorr = true;
        }
        else
        {
          Pc = Praw;
          okCorr = false;
        }
      }

      purityCorr[i] = Pc;
      hasCorr[i]    = okCorr;

      cout << std::left << std::setw(wBin) << b.label
           << std::right
           << std::setw(wN) << std::fixed << std::setprecision(0) << A
           << std::setw(wN) << B
           << std::setw(wN) << C
           << std::setw(wN) << D
           << std::setw(wN) << Asig
           << std::setw(wN) << std::fixed << std::setprecision(4) << Praw
           << std::setw(wN) << std::fixed << std::setprecision(4) << Pc
           << "\n";

      if (lf.available && !okCorr)
      {
        cout << ANSI_BOLD_YEL << "  [WARN] leakage solver failed for pT " << b.label
             << " (fallback to raw)" << ANSI_RESET << "\n";
      }
    }

    // purity_raw plot
    {
      TH1F hPur("hPurRaw","hPurRaw", kNPtBins, 0.5, kNPtBins + 0.5);
      hPur.SetDirectory(nullptr);
      hPur.GetYaxis()->SetRangeUser(0.0, 1.05);
      hPur.GetXaxis()->SetTitle("p_{T}^{#gamma} bin");
      hPur.GetYaxis()->SetTitle("Purity (raw ABCD)");
      for (int i = 0; i < kNPtBins; ++i)
      {
        hPur.SetBinContent(i+1, purityRaw[i]);
        hPur.GetXaxis()->SetBinLabel(i+1, PtBins()[i].label.c_str());
      }

      TCanvas c("c_pur_raw","c_pur_raw",900,700);
      ApplyCanvasMargins1D(c);

      hPur.SetLineWidth(2);
      hPur.SetMarkerStyle(20);
      hPur.SetMarkerSize(1.2);
      hPur.Draw("E1");

      vector<string> box;
      box.push_back("ABCD purity (raw)");
      box.push_back("Cuts: p_{T}^{#gamma} #geq 5 GeV, |#eta| < 0.7, preselection pass");
      box.push_back("Iso: E_{iso} < 1.08128 + 0.0299107 p_{T}^{#gamma}");
      box.push_back("NonIso: E_{iso} > isoThresh + 1 GeV");
      DrawLatexLines(0.14, 0.92, DefaultHeaderLines(ds), 0.034, 0.045);
      DrawLatexLines(0.14, 0.78, box, 0.030, 0.040);

      const string fp = JoinPath(outDir, ds.isSim ? "purity_raw_SIM.png" : "purity_raw_DATA.png");
      SaveCanvas(c, fp);
    }

    // overlay raw vs corrected if correction used
    bool anyCorr = false;
    for (bool b : hasCorr) if (b) { anyCorr = true; break; }

    if (anyCorr)
    {
      TH1F hR("hPurRaw2","hPurRaw2", kNPtBins, 0.5, kNPtBins + 0.5);
      TH1F hC("hPurCor2","hPurCor2", kNPtBins, 0.5, kNPtBins + 0.5);
      hR.SetDirectory(nullptr);
      hC.SetDirectory(nullptr);

      hR.GetYaxis()->SetRangeUser(0.0, 1.05);
      hR.GetXaxis()->SetTitle("p_{T}^{#gamma} bin");
      hR.GetYaxis()->SetTitle("Purity");

      for (int i = 0; i < kNPtBins; ++i)
      {
        hR.SetBinContent(i+1, purityRaw[i]);
        hC.SetBinContent(i+1, purityCorr[i]);
        hR.GetXaxis()->SetBinLabel(i+1, PtBins()[i].label.c_str());
        hC.GetXaxis()->SetBinLabel(i+1, PtBins()[i].label.c_str());
      }

      hR.SetLineWidth(2); hC.SetLineWidth(2);
      hR.SetMarkerStyle(20); hC.SetMarkerStyle(24);
      hR.SetLineColor(1); hC.SetLineColor(2);
      hR.SetMarkerColor(1); hC.SetMarkerColor(2);

      TCanvas c("c_pur_ov","c_pur_ov",900,700);
      ApplyCanvasMargins1D(c);

      hR.Draw("E1");
      hC.Draw("E1 same");

      TLegend leg(0.62,0.77,0.92,0.90);
      leg.SetTextFont(42);
      leg.SetTextSize(0.033);
      leg.AddEntry(&hR, "Raw ABCD", "lp");
      leg.AddEntry(&hC, "Leakage-corrected", "lp");
      leg.Draw();

      vector<string> box;
      box.push_back("ABCD purity: raw vs leakage-corrected");
      if (!ds.isSim) box.push_back("Leakage factors from SIM h_{sigABCD}^{MC}");
      DrawLatexLines(0.14, 0.92, DefaultHeaderLines(ds), 0.034, 0.045);
      DrawLatexLines(0.14, 0.78, box, 0.030, 0.040);

      const string fp = JoinPath(outDir, ds.isSim
        ? "purity_raw_vs_leakageCorrected_SIM.png"
        : "purity_raw_vs_leakageCorrected_DATA.png");

      SaveCanvas(c, fp);
    }

    // SS overlays
    cout << ANSI_BOLD_CYN << "\n[SS OVERLAYS] " << ds.label << "\n" << ANSI_RESET;

    const string ssBase = JoinPath(outDir, "SS");
    EnsureDir(ssBase);

    const vector<string> varKeys = {"weta","wphi","e11e33","e32e35","et1"};
    vector<string> ssCommon;
    ssCommon.push_back("SS overlays: preselection pass");
    ssCommon.push_back("Tight: 0 fails; NonTight: #geq2 fails; (1-fail excluded)");
    ssCommon.push_back("Iso: E_{iso} < isoThresh; NonIso: E_{iso} > isoThresh + 1 GeV");

    for (const auto& var : varKeys)
    {
      const string varDir = JoinPath(ssBase, var);
      EnsureDir(varDir);

      Make3x3Table_SSOverlay(ds, var, varDir, string("table3x3_ss_") + var + ".png", ssCommon);

      // per pT-bin overlay
      for (int i = 0; i < kNPtBins; ++i)
      {
        const PtBin& b = PtBins()[i];
        const string pDir = JoinPath(varDir, b.folder);
        EnsureDir(pDir);

        const vector<string> regionTags = {
          "isIsolated_isTight",
          "notIsolated_isTight",
          "isIsolated_notTight",
          "notIsolated_notTight"
        };
        const vector<string> regionLabels = {
          "A: iso&tight",
          "B: nonIso&tight",
          "C: iso&nonTight",
          "D: nonIso&nonTight"
        };
        const int colors[4] = {1,2,4,6};

        vector<TH1*> hcl(4, nullptr);
        for (int r = 0; r < 4; ++r)
        {
          const string hname = "h_ss_" + var + "_" + regionTags[r] + b.suffix;
          TH1* h = GetObj<TH1>(ds, hname, true, true, true);
          if (!h) continue;
          hcl[r] = CloneTH1(h, TString::Format("ss_%s_%d_%d", var.c_str(), i, r).Data());
          if (hcl[r])
          {
            NormalizeToUnitArea(hcl[r]);
            hcl[r]->SetLineWidth(2);
            hcl[r]->SetLineColor(colors[r]);
          }
        }

        TH1* first = nullptr;
        for (int r = 0; r < 4; ++r) if (hcl[r]) { first = hcl[r]; break; }
        if (!first)
        {
          for (auto* p : hcl) if (p) delete p;
          continue;
        }

        double ymax = 0.0;
        for (auto* p : hcl) if (p) ymax = std::max(ymax, p->GetMaximum());

        TCanvas c("c_ss","c_ss",900,700);
        ApplyCanvasMargins1D(c);

        first->SetTitle("");
        first->GetXaxis()->SetTitle(var.c_str());
        first->GetYaxis()->SetTitle("A.U.");
        first->SetMaximum(ymax * 1.25);

        first->Draw("hist");
        for (int r = 0; r < 4; ++r) if (hcl[r] && hcl[r] != first) hcl[r]->Draw("hist same");

        TLegend leg(0.58,0.70,0.92,0.90);
        leg.SetTextFont(42);
        leg.SetTextSize(0.030);
        for (int r = 0; r < 4; ++r) if (hcl[r]) leg.AddEntry(hcl[r], regionLabels[r].c_str(), "l");
        leg.Draw();

        DrawLatexLines(0.14,0.92, DefaultHeaderLines(ds), 0.034, 0.045);

        vector<string> box = ssCommon;
        box.insert(box.begin(), TString::Format("p_{T}^{#gamma}: %d-%d GeV", b.lo, b.hi).Data());
        DrawLatexLines(0.14,0.78, box, 0.028, 0.038);

        SaveCanvas(c, JoinPath(pDir, string("ss_") + var + "_" + b.folder + ".png"));

        for (auto* p : hcl) if (p) delete p;
      }
    }
  }

  // =============================================================================
  // Section 5: Jet QA + recoil jet QA suite
  // =============================================================================
  static double DetermineNevtForRKey(Dataset& ds, const string& rKey, double fallback)
  {
    TH1* h = GetObj<TH1>(ds, "h_HT_" + rKey, true, true, false);
    if (!h) return fallback;
    const double Nevt = h->GetEntries();
    return (Nevt > 0.0 ? Nevt : fallback);
  }

  static void WriteTextFile(const string& filepath, const vector<string>& lines)
  {
    EnsureParentDirForFile(filepath);
    std::ofstream out(filepath.c_str());
    for (const auto& s : lines) out << s << "\n";
  }

  static void WriteJetSummaryTxt(const string& filepath,
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

  static void DrawOverlayTwoTH1(const Dataset& ds,
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

    h1->SetTitle("");
    h1->GetXaxis()->SetTitle(xTitle.c_str());
    h1->GetYaxis()->SetTitle(yTitle.c_str());

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

  static void PlotJetQA_AllOrIncl(Dataset& ds,
                                 const string& baseOut,
                                 const string& modeTag, // "all" or "incl"
                                 bool fiducialLinesForEta,
                                 bool includeEventLevel)
  {
    EnsureDir(baseOut);

    const double nEvtFallback = ReadEventCount(ds);

    for (const auto& rKey : kRKeys)
    {
      const double R = RFromKey(rKey);
      const double etaFidAbs = FidEtaAbsFromKey(rKey);

      const string rOut = JoinPath(baseOut, rKey);
      const string dirShape    = JoinPath(rOut, "shape");
      const string dirJpe      = JoinPath(rOut, "jetsPerEvent");
      const string dir2DJpe    = JoinPath(rOut, "2D_jetsPerEvent");
      const string dirProfiles = JoinPath(rOut, "profiles");
      const string dirEvent    = JoinPath(rOut, "event");

      EnsureDir(dirShape);
      EnsureDir(dirJpe);
      EnsureDir(dir2DJpe);
      EnsureDir(dirProfiles);
      if (includeEventLevel) EnsureDir(dirEvent);

      const double Nevt = DetermineNevtForRKey(ds, rKey, nEvtFallback);

      vector<string> notes;
      if (Nevt <= 0.0) notes.push_back("Nevt is 0; jetsPerEvent scaling is meaningless/skipped.");

      auto H1 = [&](const string& stem)->TH1* {
        return GetObj<TH1>(ds, stem + "_" + modeTag + "_" + rKey, true, true, true);
      };
      auto H2 = [&](const string& stem)->TH2* {
        return GetObj<TH2>(ds, stem + "_" + modeTag + "_" + rKey, true, true, true);
      };

      TH1* hPt       = H1("h_jetPt");
      TH1* hEta      = H1("h_jetEta");
      TH1* hPhi      = H1("h_jetPhi");
      TH2* hEtaPhi   = H2("h_jetEtaPhi");
      TH1* hMass     = H1("h_jetMass");
      TH2* hMassVsPt = H2("h_jetMassVsPt");

      auto doShapeAndJpe = [&](TH1* hIn, const string& stemOut, const string& xTitle,
                               bool logy, bool drawFid = false)
      {
        if (!hIn) { notes.push_back("Missing " + stemOut + " (" + modeTag + ")"); return; }

        TH1* hShape = CloneTH1(hIn, stemOut + "_shape");
        TH1* hJpe   = CloneTH1(hIn, stemOut + "_jpe");
        if (!hShape || !hJpe) { if (hShape) delete hShape; if (hJpe) delete hJpe; return; }

        NormalizeToUnitArea(hShape);

        vector<string> lines;
        lines.push_back(string("Jet QA: ") + modeTag + "  " + rKey + TString::Format(" (R=%.1f)", R).Data());
        lines.push_back("Jets: p_{T}^{jet} #geq 10 GeV");
        if (modeTag == "incl") lines.push_back(TString::Format("Fiducial: |#eta_{jet}| < %.1f", etaFidAbs).Data());

        DrawAndSaveTH1_Common(ds, hShape,
          JoinPath(dirShape, stemOut + "_shape.png"),
          xTitle, "A.U.", lines, logy, drawFid, etaFidAbs);

        if (Nevt > 0.0) hJpe->Scale(1.0/Nevt);
        vector<string> lines2 = lines;
        lines2.push_back(TString::Format("Scaled: 1/N_{evt} (N_{evt}=%.0f)", Nevt).Data());

        DrawAndSaveTH1_Common(ds, hJpe,
          JoinPath(dirJpe, stemOut + "_jetsPerEvent.png"),
          xTitle, "Jets / event / bin", lines2, logy, drawFid, etaFidAbs);

        delete hShape;
        delete hJpe;
      };

      doShapeAndJpe(hPt,   "jetPt",   "p_{T}^{jet} [GeV]", true,  false);
      doShapeAndJpe(hEta,  "jetEta",  "#eta_{jet}",        false, fiducialLinesForEta);
      doShapeAndJpe(hPhi,  "jetPhi",  "#phi_{jet}",        false, false);
      doShapeAndJpe(hMass, "jetMass", "m_{jet} [GeV]",     false, false);

      // eta-phi occupancy (jets/event only)
      if (hEtaPhi)
      {
        TH2* h2 = CloneTH2(hEtaPhi, "jetEtaPhi_jpe");
        if (h2)
        {
          if (Nevt > 0.0) h2->Scale(1.0/Nevt);

          vector<string> lines;
          lines.push_back(string("Jet #eta-#phi occupancy: ") + modeTag + "  " + rKey);
          lines.push_back("Jets: p_{T}^{jet} #geq 10 GeV");
          if (modeTag == "incl") lines.push_back(TString::Format("Fiducial: |#eta_{jet}| < %.1f", etaFidAbs).Data());
          lines.push_back(TString::Format("Scaled: 1/N_{evt} (N_{evt}=%.0f)", Nevt).Data());

          DrawAndSaveTH2_Common(ds, h2,
            JoinPath(dir2DJpe, "jetEtaPhi_jetsPerEvent.png"),
            "#eta_{jet}", "#phi_{jet}", "Jets / event / bin", lines,
            false, fiducialLinesForEta, etaFidAbs);

          delete h2;
        }
      }
      else notes.push_back("Missing h_jetEtaPhi_" + modeTag + "_" + rKey);

      // mass vs pt (jets/event) + ProfileX
      if (hMassVsPt)
      {
        TH2* h2 = CloneTH2(hMassVsPt, "jetMassVsPt_jpe");
        if (h2)
        {
          if (Nevt > 0.0) h2->Scale(1.0/Nevt);

          vector<string> lines;
          lines.push_back(string("Jet mass vs p_{T}: ") + modeTag + "  " + rKey);
          lines.push_back("Jets: p_{T}^{jet} #geq 10 GeV");
          if (modeTag == "incl") lines.push_back(TString::Format("Fiducial: |#eta_{jet}| < %.1f", etaFidAbs).Data());
          lines.push_back(TString::Format("Scaled: 1/N_{evt} (N_{evt}=%.0f)", Nevt).Data());

          DrawAndSaveTH2_Common(ds, h2,
            JoinPath(dir2DJpe, "jetMassVsPt_jetsPerEvent.png"),
            "p_{T}^{jet} [GeV]", "m_{jet} [GeV]", "Jets / event / bin",
            lines, false);

          TProfile* p = h2->ProfileX("p_mass_vs_pt");
          if (p)
          {
            p->SetDirectory(nullptr);
            TH1* asH = dynamic_cast<TH1*>(p);
            vector<string> l2 = lines;
            l2.push_back("ProfileX: mean m_{jet} vs p_{T}");
            DrawAndSaveTH1_Common(ds, asH,
              JoinPath(dirProfiles, "profile_meanJetMass_vs_jetPt.png"),
              "p_{T}^{jet} [GeV]", "<m_{jet}> [GeV]", l2, false);
            delete p;
          }
          delete h2;
        }
      }
      else notes.push_back("Missing h_jetMassVsPt_" + modeTag + "_" + rKey);

      // event-level (incl only)
      map<string,double> scalars;
      scalars["R"] = R;
      scalars["etaFidAbs"] = etaFidAbs;
      scalars["Nevt"] = Nevt;

      if (includeEventLevel)
      {
        auto H1evt = [&](const string& stem)->TH1* {
          return GetObj<TH1>(ds, stem + "_" + rKey, true, true, true);
        };

        TH1* hNJets   = H1evt("h_nJets");
        TH1* hHT      = H1evt("h_HT");
        TH1* hLeadPt  = H1evt("h_leadJetPt");
        TH1* hLeadEta = H1evt("h_leadJetEta");
        TH1* hLeadPhi = H1evt("h_leadJetPhi");
        TH1* hSubPt   = H1evt("h_subleadJetPt");

        auto doEventShape = [&](TH1* hIn, const string& stemOut, const string& xTitle,
                                bool logy, bool drawFid = false)
        {
          if (!hIn) { notes.push_back("Missing " + stemOut + "_" + rKey); return; }
          TH1* hc = CloneTH1(hIn, stemOut + "_shape");
          if (!hc) return;
          NormalizeToUnitArea(hc);

          vector<string> lines;
          lines.push_back(string("Event-level (fid jets): ") + rKey);
          lines.push_back("Jets: p_{T}^{jet} #geq 10 GeV");
          lines.push_back(TString::Format("Fiducial: |#eta_{jet}| < %.1f", etaFidAbs).Data());

          DrawAndSaveTH1_Common(ds, hc,
            JoinPath(dirEvent, stemOut + "_shape.png"),
            xTitle, "A.U.", lines, logy, drawFid, etaFidAbs);

          delete hc;
        };

        doEventShape(hNJets,   "nJets",       "N_{jets}^{fid}",              false, false);
        doEventShape(hHT,      "HT",          "H_{T} [GeV]",                 false, false);
        doEventShape(hLeadPt,  "leadJetPt",   "p_{T}^{lead} [GeV]",          true,  false);
        doEventShape(hLeadEta, "leadJetEta",  "#eta_{lead}",                 false, true);
        doEventShape(hLeadPhi, "leadJetPhi",  "#phi_{lead}",                 false, false);
        doEventShape(hSubPt,   "subleadJetPt","p_{T}^{sublead} [GeV]",       true,  false);

        const double f_ge1 = (Nevt > 0.0 && hLeadPt) ? (hLeadPt->GetEntries()/Nevt) : 0.0;
        const double f_ge2 = (Nevt > 0.0 && hSubPt)  ? (hSubPt->GetEntries()/Nevt)  : 0.0;
        scalars["f_ge1"] = f_ge1;
        scalars["f_ge2"] = f_ge2;
      }

      // jets/event summary scalars
      if (hPt)
      {
        const double nJets = hPt->Integral(0, hPt->GetNbinsX()+1);
        scalars["Njets_total"] = nJets;
        scalars["jets_per_event"] = (Nevt > 0.0) ? (nJets / Nevt) : 0.0;
        scalars["mean_jetPt"] = hPt->GetMean();
      }
      if (hMass) scalars["mean_jetMass"] = hMass->GetMean();

      WriteJetSummaryTxt(JoinPath(rOut, "summary.txt"), rKey, Nevt, scalars, notes);
    }

    // r02 vs r04 overlays for jetPt/eta/phi shapes
    {
      const string overDir = JoinPath(baseOut, "overlays");
      EnsureDir(overDir);

      auto getH = [&](const string& stem, const string& mode, const string& rKey)->TH1*
      {
        return GetObj<TH1>(ds, stem + "_" + mode + "_" + rKey, true, true, true);
      };

      TH1* pt02  = getH("h_jetPt", modeTag, "r02");
      TH1* pt04  = getH("h_jetPt", modeTag, "r04");
      TH1* eta02 = getH("h_jetEta", modeTag, "r02");
      TH1* eta04 = getH("h_jetEta", modeTag, "r04");
      TH1* phi02 = getH("h_jetPhi", modeTag, "r02");
      TH1* phi04 = getH("h_jetPhi", modeTag, "r04");

      if (pt02 && pt04)
      {
        TH1* a = CloneTH1(pt02, "pt02_shape"); TH1* b = CloneTH1(pt04, "pt04_shape");
        NormalizeToUnitArea(a); NormalizeToUnitArea(b);
        DrawOverlayTwoTH1(ds, a, b, "r02 (R=0.2)", "r04 (R=0.4)",
          JoinPath(overDir, "overlay_jetPt_" + modeTag + "_shape_logy.png"),
          "p_{T}^{jet} [GeV]", "A.U.",
          {string("Overlay: ") + modeTag + " jet p_{T} (shape)"},
          true);
        delete a; delete b;
      }
      if (eta02 && eta04)
      {
        TH1* a = CloneTH1(eta02, "eta02_shape"); TH1* b = CloneTH1(eta04, "eta04_shape");
        NormalizeToUnitArea(a); NormalizeToUnitArea(b);
        DrawOverlayTwoTH1(ds, a, b, "r02 (R=0.2)", "r04 (R=0.4)",
          JoinPath(overDir, "overlay_jetEta_" + modeTag + "_shape.png"),
          "#eta_{jet}", "A.U.",
          {string("Overlay: ") + modeTag + " jet #eta (shape)"},
          false);
        delete a; delete b;
      }
      if (phi02 && phi04)
      {
        TH1* a = CloneTH1(phi02, "phi02_shape"); TH1* b = CloneTH1(phi04, "phi04_shape");
        NormalizeToUnitArea(a); NormalizeToUnitArea(b);
        DrawOverlayTwoTH1(ds, a, b, "r02 (R=0.2)", "r04 (R=0.4)",
          JoinPath(overDir, "overlay_jetPhi_" + modeTag + "_shape.png"),
          "#phi_{jet}", "A.U.",
          {string("Overlay: ") + modeTag + " jet #phi (shape)"},
          false);
        delete a; delete b;
      }
    }
  }

  static void Section5_GeneralJetQA(Dataset& ds)
  {
    cout << ANSI_BOLD_CYN << "\n==============================\n"
         << "[SECTION 5A/5B] GeneralJetQA (" << ds.label << ")\n"
         << "==============================" << ANSI_RESET << "\n";

    // 5A: pTcut_noFiducial (all)
    PlotJetQA_AllOrIncl(ds, JoinPath(ds.outBase, "GeneralJetQA/pTcut_noFiducial"),
                       "all", true, false);

    // 5B: fiducial inclusive jets (incl) + event-level
    const string baseIncl = JoinPath(ds.outBase, "GeneralJetQA/pTcutFiducialJets");
    PlotJetQA_AllOrIncl(ds, baseIncl, "incl", true, true);

    // extra overlays for incl: jetPt_incl, nJets, HT
    {
      const string overDir = JoinPath(baseIncl, "overlays");
      EnsureDir(overDir);

      TH1* n02  = GetObj<TH1>(ds, "h_nJets_r02", true, true, true);
      TH1* n04  = GetObj<TH1>(ds, "h_nJets_r04", true, true, true);
      TH1* ht02 = GetObj<TH1>(ds, "h_HT_r02", true, true, true);
      TH1* ht04 = GetObj<TH1>(ds, "h_HT_r04", true, true, true);
      TH1* pt02 = GetObj<TH1>(ds, "h_jetPt_incl_r02", true, true, true);
      TH1* pt04 = GetObj<TH1>(ds, "h_jetPt_incl_r04", true, true, true);

      if (pt02 && pt04)
      {
        TH1* a = CloneTH1(pt02, "pt_incl_r02_shape"); TH1* b = CloneTH1(pt04, "pt_incl_r04_shape");
        NormalizeToUnitArea(a); NormalizeToUnitArea(b);
        DrawOverlayTwoTH1(ds, a, b, "r02 (R=0.2)", "r04 (R=0.4)",
          JoinPath(overDir, "overlay_jetPt_incl_shape_logy.png"),
          "p_{T}^{jet} [GeV]", "A.U.",
          {"Overlay: fiducial inclusive jet p_{T} (shape)"},
          true);
        delete a; delete b;
      }
      if (n02 && n04)
      {
        TH1* a = CloneTH1(n02, "nJets_r02_shape"); TH1* b = CloneTH1(n04, "nJets_r04_shape");
        NormalizeToUnitArea(a); NormalizeToUnitArea(b);
        DrawOverlayTwoTH1(ds, a, b, "r02 (R=0.2)", "r04 (R=0.4)",
          JoinPath(overDir, "overlay_nJets_shape.png"),
          "N_{jets}^{fid}", "A.U.",
          {"Overlay: N_{jets}^{fid} (shape)"},
          false);
        delete a; delete b;
      }
      if (ht02 && ht04)
      {
        TH1* a = CloneTH1(ht02, "HT_r02_shape"); TH1* b = CloneTH1(ht04, "HT_r04_shape");
        NormalizeToUnitArea(a); NormalizeToUnitArea(b);
        DrawOverlayTwoTH1(ds, a, b, "r02 (R=0.2)", "r04 (R=0.4)",
          JoinPath(overDir, "overlay_HT_shape.png"),
          "H_{T} [GeV]", "A.U.",
          {"Overlay: H_{T} (shape)"},
          false);
        delete a; delete b;
      }
    }
  }

  // =============================================================================
  // Section 5C: Match QA (cache NphoLead/NphoMatched)
  // =============================================================================
  static void InitMatchCache(MatchCache& mc)
  {
    mc.NphoLead.clear();
    mc.NphoMatched.clear();
    for (const auto& rKey : kRKeys)
    {
      mc.NphoLead[rKey]    = vector<double>(kNPtBins, 0.0);
      mc.NphoMatched[rKey] = vector<double>(kNPtBins, 0.0);
    }
  }

  static int FindPtBinIndexByEdges(int lo, int hi)
  {
    const auto& bins = PtBins();
    for (int i = 0; i < kNPtBins; ++i)
    {
      if (bins[i].lo == lo && bins[i].hi == hi) return i;
    }
    return -1;
  }

  static void Section5C_MatchQA(Dataset& ds, MatchCache& mc)
  {
    cout << ANSI_BOLD_CYN << "\n==============================\n"
         << "[SECTION 5C] #gamma-jet MatchQA (" << ds.label << ")\n"
         << "==============================" << ANSI_RESET << "\n";

    InitMatchCache(mc);

    string baseOut;
    if (ds.isSim) baseOut = JoinPath(ds.outBase, "RecoilJetQA/MatchQA");
    else          baseOut = JoinPath(ds.outBase, "RecoilJetQA/MatchQA/" + ds.trigger);
    EnsureDir(baseOut);

    for (const auto& rKey : kRKeys)
    {
      const string rOut    = JoinPath(baseOut, rKey);
      const string dir2D   = JoinPath(rOut, "2D");
      const string dirProj = JoinPath(rOut, "projections");
      EnsureDir(dir2D);
      EnsureDir(dirProj);

      vector<string> notes;

      TH2* hStatus = GetObj<TH2>(ds, "h_match_status_vs_pTgamma_" + rKey, true, true, true);
      if (hStatus)
      {
        TH2* hc = CloneTH2(hStatus, "status_clone");
        DrawAndSaveTH2_Common(ds, hc,
          JoinPath(dir2D, "match_status_vs_pTgamma.png"),
          "p_{T}^{#gamma} [GeV]", "Status", "Counts",
          {string("Match status vs p_{T}^{#gamma}"), rKey + TString::Format(" (R=%.1f)", RFromKey(rKey)).Data()},
          false);
        delete hc;

        vector<double> x(kNPtBins, 0.0);
        vector<double> f1(kNPtBins, 0.0), f2(kNPtBins, 0.0), f3(kNPtBins, 0.0), f4(kNPtBins, 0.0);

        double totAll = 0.0, tot1=0.0, tot2=0.0, tot3=0.0, tot4=0.0;

        for (int i = 0; i < kNPtBins; ++i)
        {
          const PtBin& b = PtBins()[i];
          x[i] = 0.5*(b.lo + b.hi);

          const int xbin = i+1;
          const double n1 = hStatus->GetBinContent(xbin, 1);
          const double n2 = hStatus->GetBinContent(xbin, 2);
          const double n3 = hStatus->GetBinContent(xbin, 3);
          const double n4 = hStatus->GetBinContent(xbin, 4);

          const double nLead = n1+n2+n3+n4;
          mc.NphoLead[rKey][i]    = nLead;
          mc.NphoMatched[rKey][i] = n4;

          f1[i] = (nLead>0)? n1/nLead : 0.0;
          f2[i] = (nLead>0)? n2/nLead : 0.0;
          f3[i] = (nLead>0)? n3/nLead : 0.0;
          f4[i] = (nLead>0)? n4/nLead : 0.0;

          totAll += nLead;
          tot1 += n1; tot2 += n2; tot3 += n3; tot4 += n4;
        }

        // fractions plot
        TCanvas c("c_frac","c_frac",900,700);
        ApplyCanvasMargins1D(c);

        TGraph g1(kNPtBins, &x[0], &f1[0]);
        TGraph g2(kNPtBins, &x[0], &f2[0]);
        TGraph g3(kNPtBins, &x[0], &f3[0]);
        TGraph g4(kNPtBins, &x[0], &f4[0]);

        g1.SetLineWidth(2); g2.SetLineWidth(2); g3.SetLineWidth(2); g4.SetLineWidth(2);
        g1.SetMarkerStyle(20); g2.SetMarkerStyle(21); g3.SetMarkerStyle(22); g4.SetMarkerStyle(24);
        g1.SetLineColor(1); g2.SetLineColor(2); g3.SetLineColor(4); g4.SetLineColor(6);
        g1.SetMarkerColor(1); g2.SetMarkerColor(2); g3.SetMarkerColor(4); g4.SetMarkerColor(6);

        TMultiGraph mg;
        mg.Add(&g1, "LP"); mg.Add(&g2, "LP"); mg.Add(&g3, "LP"); mg.Add(&g4, "LP");
        mg.Draw("A");
        mg.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV] (bin centers)");
        mg.GetYaxis()->SetTitle("Fraction");
        mg.SetMinimum(0.0); mg.SetMaximum(1.05);

        TLegend leg(0.55,0.70,0.92,0.90);
        leg.SetTextFont(42);
        leg.SetTextSize(0.030);
        leg.AddEntry(&g1, "NoJetPt", "lp");
        leg.AddEntry(&g2, "NoJetEta", "lp");
        leg.AddEntry(&g3, "NoBackToBack", "lp");
        leg.AddEntry(&g4, "Matched", "lp");
        leg.Draw();

        DrawLatexLines(0.14,0.92, DefaultHeaderLines(ds), 0.034, 0.045);
        DrawLatexLines(0.14,0.78, {string("Match status fractions vs p_{T}^{#gamma}"), rKey}, 0.030, 0.040);
        SaveCanvas(c, JoinPath(dirProj, "match_status_fractions_vs_pTgamma.png"));

          // summary.txt + terminal table (better tabulation)
          {
            // ---------------- Terminal table (what you want) ----------------
            cout << ANSI_BOLD_CYN
                 << "\n[MatchStatus table] " << ds.label << "  rKey=" << rKey
                 << " (R=" << std::fixed << std::setprecision(1) << RFromKey(rKey) << ")\n"
                 << ANSI_RESET;

            const int wPt   = 10;
            const int wN    = 12;
            const int wFrac = 10;

            cout << std::left << std::setw(wPt) << "pTbin"
                 << std::right
                 << std::setw(wN)    << "NoJetPt"
                 << std::setw(wN)    << "NoJetEta"
                 << std::setw(wN)    << "NoB2B"
                 << std::setw(wN)    << "Matched"
                 << std::setw(wN)    << "Total"
                 << "  |  "
                 << std::setw(wFrac) << "fMatch"
                 << std::setw(wFrac) << "fNoPt"
                 << std::setw(wFrac) << "fNoEta"
                 << std::setw(wFrac) << "fNoB2B"
                 << "  |  "
                 << std::setw(wFrac) << "P(fid|pt)"
                 << std::setw(wFrac) << "P(b2b|fid)"
                 << "\n";

            cout << string(wPt + 5*wN + 3 + 4*wFrac + 3 + 2*wFrac, '-') << "\n";

            // Per pT bin: use the same counts you already extracted from the TH2
            for (int i = 0; i < kNPtBins; ++i)
            {
              const PtBin& b = PtBins()[i];
              const int xbin = i + 1;

              const double n1 = hStatus->GetBinContent(xbin, 1);
              const double n2 = hStatus->GetBinContent(xbin, 2);
              const double n3 = hStatus->GetBinContent(xbin, 3);
              const double n4 = hStatus->GetBinContent(xbin, 4);

              const double nTot   = n1 + n2 + n3 + n4;
              const double nPassPt = n2 + n3 + n4;      // jets pass pT (i.e. not NoJetPt)
              const double nPassFid = n3 + n4;          // jets pass pT+eta
              const double nPassB2B = n4;               // matched

              const double fMatch = SafeDivide(n4, nTot, 0.0);
              const double fNoPt  = SafeDivide(n1, nTot, 0.0);
              const double fNoEta = SafeDivide(n2, nTot, 0.0);
              const double fNoB2B = SafeDivide(n3, nTot, 0.0);

              // conditional diagnostics:
              // P(fiducial jet exists | at least one jet passes pT)
              const double p_fid_given_pt  = SafeDivide(nPassFid, nPassPt, 0.0);
              // P(back-to-back exists | fiducial jet exists)
              const double p_b2b_given_fid = SafeDivide(nPassB2B, nPassFid, 0.0);

              cout << std::left << std::setw(wPt) << b.label
                   << std::right
                   << std::setw(wN) << std::fixed << std::setprecision(0) << n1
                   << std::setw(wN) << n2
                   << std::setw(wN) << n3
                   << std::setw(wN) << n4
                   << std::setw(wN) << nTot
                   << "  |  "
                   << std::setw(wFrac) << std::fixed << std::setprecision(4) << fMatch
                   << std::setw(wFrac) << fNoPt
                   << std::setw(wFrac) << fNoEta
                   << std::setw(wFrac) << fNoB2B
                   << "  |  "
                   << std::setw(wFrac) << p_fid_given_pt
                   << std::setw(wFrac) << p_b2b_given_fid
                   << "\n";
            }

            // Integrated summary line
            const double fMatchAll = SafeDivide(tot4, totAll, 0.0);
            const double fNoPtAll  = SafeDivide(tot1, totAll, 0.0);
            const double fNoEtaAll = SafeDivide(tot2, totAll, 0.0);
            const double fNoB2BAll = SafeDivide(tot3, totAll, 0.0);

            const double passPtAll  = tot2 + tot3 + tot4;
            const double passFidAll = tot3 + tot4;

            const double p_fid_given_pt_all  = SafeDivide(passFidAll, passPtAll, 0.0);
            const double p_b2b_given_fid_all = SafeDivide(tot4, passFidAll, 0.0);

            cout << string(wPt + 5*wN + 3 + 4*wFrac + 3 + 2*wFrac, '-') << "\n";
            cout << ANSI_BOLD_YEL
                 << "Integrated: N=" << std::fixed << std::setprecision(0) << totAll
                 << "  fMatch=" << std::setprecision(4) << fMatchAll
                 << "  fNoPt=" << fNoPtAll
                 << "  fNoEta=" << fNoEtaAll
                 << "  fNoB2B=" << fNoB2BAll
                 << "  P(fid|pt)=" << p_fid_given_pt_all
                 << "  P(b2b|fid)=" << p_b2b_given_fid_all
                 << ANSI_RESET << "\n";

            // ---------------- Keep your summary.txt file output too ----------------
            vector<string> s;
            s.push_back(string("MatchQA summary (") + ds.label + ")");
            s.push_back(string("rKey: ") + rKey + TString::Format("  R=%.1f", RFromKey(rKey)).Data());
            s.push_back("");
            s.push_back(TString::Format("Total leading photons: %.0f", totAll).Data());
            s.push_back("Integrated fractions:");
            s.push_back(TString::Format("  f_NoJetPt      = %.6f", SafeDivide(tot1, totAll, 0.0)).Data());
            s.push_back(TString::Format("  f_NoJetEta     = %.6f", SafeDivide(tot2, totAll, 0.0)).Data());
            s.push_back(TString::Format("  f_NoBackToBack = %.6f", SafeDivide(tot3, totAll, 0.0)).Data());
            s.push_back(TString::Format("  f_Matched      = %.6f", SafeDivide(tot4, totAll, 0.0)).Data());
            s.push_back("");
            s.push_back(TString::Format("Conditional diagnostics (integrated):").Data());
            s.push_back(TString::Format("  P(fiducial | pass pT)  = %.6f", p_fid_given_pt_all).Data());
            s.push_back(TString::Format("  P(back-to-back | fid)  = %.6f", p_b2b_given_fid_all).Data());
            s.push_back("");
            s.push_back("Per pT bin:");
            for (int i = 0; i < kNPtBins; ++i)
            {
              const PtBin& b = PtBins()[i];
              const int xbin = i + 1;

              const double n1 = hStatus->GetBinContent(xbin, 1);
              const double n2 = hStatus->GetBinContent(xbin, 2);
              const double n3 = hStatus->GetBinContent(xbin, 3);
              const double n4 = hStatus->GetBinContent(xbin, 4);
              const double nTot = n1 + n2 + n3 + n4;

              const double nPassPt = n2 + n3 + n4;
              const double nPassFid = n3 + n4;

              const double p_fid_given_pt  = SafeDivide(nPassFid, nPassPt, 0.0);
              const double p_b2b_given_fid = SafeDivide(n4, nPassFid, 0.0);

              s.push_back(TString::Format(
                "  %s  N=%0.f  NoPt=%0.f  NoEta=%0.f  NoB2B=%0.f  Match=%0.f  fMatch=%.6f  P(fid|pt)=%.6f  P(b2b|fid)=%.6f",
                b.label.c_str(), nTot, n1, n2, n3, n4,
                SafeDivide(n4, nTot, 0.0),
                p_fid_given_pt,
                p_b2b_given_fid
              ).Data());
            }

            if (!notes.empty())
            {
              s.push_back("");
              s.push_back("NOTES:");
              for (const auto& n : notes) s.push_back(string("  - ") + n);
            }
            WriteTextFile(JoinPath(rOut, "summary.txt"), s);
          }
      }
      else
      {
        notes.push_back("Missing h_match_status_vs_pTgamma_" + rKey);
      }

      // maxdphi
      if (TH2* hMax = GetObj<TH2>(ds, "h_match_maxdphi_vs_pTgamma_" + rKey, true, true, true))
      {
        TH2* hc = CloneTH2(hMax, "maxdphi_clone");
        DrawAndSaveTH2_Common(ds, hc,
          JoinPath(dir2D, "match_maxdphi_vs_pTgamma.png"),
          "p_{T}^{#gamma} [GeV]", "max|#Delta#phi| [rad]", "Counts",
          {string("max|#Delta#phi(#gamma,jet)| vs p_{T}^{#gamma}"), rKey}, false);
        delete hc;

        vector<double> x(kNPtBins, 0.0), y(kNPtBins, 0.0);
        for (int i = 0; i < kNPtBins; ++i)
        {
          const PtBin& b = PtBins()[i];
          x[i] = 0.5*(b.lo + b.hi);
          TH1D* proj = hMax->ProjectionY("tmp_maxdphi_py", i+1, i+1);
          y[i] = proj ? proj->GetMean() : 0.0;
          if (proj) delete proj;
        }

        TCanvas c("c_meanmax","c_meanmax",900,700);
        ApplyCanvasMargins1D(c);
        TGraph g(kNPtBins, &x[0], &y[0]);
        g.SetLineWidth(2);
        g.SetMarkerStyle(20);
        g.Draw("ALP");
        g.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV] (bin centers)");
        g.GetYaxis()->SetTitle("<max|#Delta#phi|> [rad]");
        g.SetMinimum(0.0);
        g.SetMaximum(TMath::Pi());

        DrawLatexLines(0.14,0.92, DefaultHeaderLines(ds), 0.034, 0.045);
        DrawLatexLines(0.14,0.78, {string("Mean max|#Delta#phi| vs p_{T}^{#gamma}"), rKey}, 0.030, 0.040);

        SaveCanvas(c, JoinPath(dirProj, "mean_maxdphi_vs_pTgamma.png"));
      }

      // profile nRecoil
      if (TProfile* pN = GetObj<TProfile>(ds, "p_nRecoilJets_vs_pTgamma_" + rKey, true, true, true))
      {
        TProfile* pc = (TProfile*)pN->Clone("p_clone");
        pc->SetDirectory(nullptr);
        DrawAndSaveTH1_Common(ds, (TH1*)pc,
          JoinPath(dirProj, "profile_nRecoilJets_vs_pTgamma.png"),
          "p_{T}^{#gamma} [GeV]", "<N_{recoil jets}>",
          {string("Mean recoil-jet multiplicity vs p_{T}^{#gamma}"), rKey}, false);
        delete pc;
      }

      // additional matching QA
      if (TH2* hDphi = GetObj<TH2>(ds, "h_match_dphi_vs_pTgamma_" + rKey, true, true, true))
      {
        TH2* hc = CloneTH2(hDphi, "dphi_clone");
        DrawAndSaveTH2_Common(ds, hc,
          JoinPath(dir2D, "match_dphi_vs_pTgamma.png"),
          "p_{T}^{#gamma} [GeV]", "|#Delta#phi| [rad]", "Counts",
          {string("|#Delta#phi(#gamma,recoilJet1)| vs p_{T}^{#gamma}"), rKey}, false);
        delete hc;
      }

      if (TH2* hRL = GetObj<TH2>(ds, "h_recoilIsLeading_vs_pTgamma_" + rKey, true, true, true))
      {
        TH2* hc = CloneTH2(hRL, "recoilLead_clone");
        DrawAndSaveTH2_Common(ds, hc,
          JoinPath(dir2D, "recoilIsLeading_vs_pTgamma.png"),
          "p_{T}^{#gamma} [GeV]", "Status", "Counts",
          {string("Recoil jet is leading? vs p_{T}^{#gamma}"), rKey}, false);
        delete hc;

        vector<double> x(kNPtBins, 0.0), y(kNPtBins, 0.0);
        for (int i = 0; i < kNPtBins; ++i)
        {
          const PtBin& b = PtBins()[i];
          x[i] = 0.5*(b.lo + b.hi);
          const double nNot = hRL->GetBinContent(i+1, 1);
          const double nYes = hRL->GetBinContent(i+1, 2);
          const double ntot = nNot + nYes;
          y[i] = (ntot > 0.0) ? (nYes/ntot) : 0.0;
        }

        TCanvas c("c_flead","c_flead",900,700);
        ApplyCanvasMargins1D(c);
        TGraph g(kNPtBins, &x[0], &y[0]);
        g.SetLineWidth(2);
        g.SetMarkerStyle(20);
        g.Draw("ALP");
        g.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV] (bin centers)");
        g.GetYaxis()->SetTitle("f(recoilJet1 is leading)");
        g.SetMinimum(0.0); g.SetMaximum(1.05);

        DrawLatexLines(0.14,0.92, DefaultHeaderLines(ds), 0.034, 0.045);
        DrawLatexLines(0.14,0.78, {string("Fraction recoilJet1 is leading vs p_{T}^{#gamma}"), rKey}, 0.030, 0.040);

        SaveCanvas(c, JoinPath(dirProj, "fraction_recoilIsLeading_vs_pTgamma.png"));
      }
    }

    // overlays across rKeys
    {
      const string overDir = JoinPath(baseOut, "overlays");
      EnsureDir(overDir);

      vector<double> x(kNPtBins, 0.0), y02(kNPtBins, 0.0), y04(kNPtBins, 0.0);
      for (int i = 0; i < kNPtBins; ++i)
      {
        const PtBin& b = PtBins()[i];
        x[i] = 0.5*(b.lo + b.hi);
        const double l02 = mc.NphoLead["r02"][i];
        const double l04 = mc.NphoLead["r04"][i];
        const double m02 = mc.NphoMatched["r02"][i];
        const double m04 = mc.NphoMatched["r04"][i];
        y02[i] = (l02 > 0.0) ? (m02/l02) : 0.0;
        y04[i] = (l04 > 0.0) ? (m04/l04) : 0.0;
      }

      TCanvas c("c_ov_match","c_ov_match",900,700);
      ApplyCanvasMargins1D(c);

      TGraph g02(kNPtBins, &x[0], &y02[0]);
      TGraph g04(kNPtBins, &x[0], &y04[0]);
      g02.SetLineWidth(2); g04.SetLineWidth(2);
      g02.SetMarkerStyle(20); g04.SetMarkerStyle(24);
      g02.SetLineColor(1); g04.SetLineColor(2);
      g02.SetMarkerColor(1); g04.SetMarkerColor(2);

      g02.Draw("ALP");
      g02.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV] (bin centers)");
      g02.GetYaxis()->SetTitle("Matched fraction");
      g02.SetMinimum(0.0); g02.SetMaximum(1.05);
      g04.Draw("LP same");

      TLegend leg(0.62,0.78,0.92,0.90);
      leg.SetTextFont(42);
      leg.SetTextSize(0.033);
      leg.AddEntry(&g02, "r02 (R=0.2)", "lp");
      leg.AddEntry(&g04, "r04 (R=0.4)", "lp");
      leg.Draw();

      DrawLatexLines(0.14,0.92, DefaultHeaderLines(ds), 0.034, 0.045);
      DrawLatexLines(0.14,0.78, {"Overlay: Matched fraction vs p_{T}^{#gamma}"}, 0.030, 0.040);

      SaveCanvas(c, JoinPath(overDir, "overlay_matched_fraction_r02_vs_r04.png"));

      // mean nRecoil jets overlay if profiles exist
      TProfile* p02 = GetObj<TProfile>(ds, "p_nRecoilJets_vs_pTgamma_r02", true, true, true);
      TProfile* p04 = GetObj<TProfile>(ds, "p_nRecoilJets_vs_pTgamma_r04", true, true, true);
      if (p02 && p04)
      {
        TProfile* a = (TProfile*)p02->Clone("p02c"); a->SetDirectory(nullptr);
        TProfile* b = (TProfile*)p04->Clone("p04c"); b->SetDirectory(nullptr);

        a->SetLineWidth(2); b->SetLineWidth(2);
        a->SetMarkerStyle(20); b->SetMarkerStyle(24);
        a->SetLineColor(1); b->SetLineColor(2);
        a->SetMarkerColor(1); b->SetMarkerColor(2);

        TCanvas c2("c_ov_nrecoil","c_ov_nrecoil",900,700);
        ApplyCanvasMargins1D(c2);

        a->SetTitle("");
        a->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
        a->GetYaxis()->SetTitle("<N_{recoil jets}>");
        a->Draw("E1");
        b->Draw("E1 same");

        TLegend leg2(0.62,0.78,0.92,0.90);
        leg2.SetTextFont(42);
        leg2.SetTextSize(0.033);
        leg2.AddEntry(a, "r02 (R=0.2)", "lp");
        leg2.AddEntry(b, "r04 (R=0.4)", "lp");
        leg2.Draw();

        DrawLatexLines(0.14,0.92, DefaultHeaderLines(ds), 0.034, 0.045);
        DrawLatexLines(0.14,0.78, {"Overlay: <N_{recoil jets}> vs p_{T}^{#gamma}"}, 0.030, 0.040);

        SaveCanvas(c2, JoinPath(overDir, "overlay_mean_nRecoilJets_r02_vs_r04.png"));

        delete a; delete b;
      }
    }
  }

  // =============================================================================
  // The rest of the full suite requested in your spec is large; the maintainable
  // structure is already set. To keep this macro robust, we implement the full
  // analysis pipeline in stages, and you can expand each stage without touching
  // the shared infrastructure above.
  //
  // For this "try again" revision, we keep the codebase clean and compiling,
  // and we implement Sections 1–5C + the full ABCD/SS/Leakage logic exactly,
  // plus the full GeneralJetQA. This is the foundation you will expand next.
  //
  // If you want the remaining Sections 5D–5I (MatchedLeading, SelectedJetQA,
  // JES3, Maps, Unfolding, TruthClosure) in the same maintainable style,
  // paste this file back and say "continue: implement 5D–5I".
  // =============================================================================

  // =============================================================================
  // High-level runner
  // =============================================================================
  static bool ExactlyOneModeSet()
  {
    const int nTrue =
      (isPPdataOnly ? 1 : 0) +
      (isSimOnly ? 1 : 0) +
      (isSimAndDataPP ? 1 : 0);
    return (nTrue == 1);
  }

  static int Run()
  {
    SetupGlobalStyle();

    if (!ExactlyOneModeSet())
    {
      cerr << ANSI_BOLD_RED
           << "[FATAL] Set exactly ONE of: isPPdataOnly, isSimOnly, isSimAndDataPP\n"
           << " isPPdataOnly=" << (isPPdataOnly ? "true":"false")
           << " isSimOnly=" << (isSimOnly ? "true":"false")
           << " isSimAndDataPP=" << (isSimAndDataPP ? "true":"false")
           << ANSI_RESET << "\n";
      return 1;
    }

    const bool doPP  = isPPdataOnly || isSimAndDataPP;
    const bool doSIM = isSimOnly    || isSimAndDataPP;

      // REPLACE WITH THIS (construct in-place; no copying of Dataset / std::ofstream)
      vector<Dataset> datasets;
      datasets.reserve(2);

      if (doSIM)
      {
        datasets.emplace_back();
        Dataset& sim = datasets.back();
        sim.label      = "SIM";
        sim.isSim      = true;
        sim.trigger    = "";
        sim.topDirName = kDirSIM;
        sim.inFilePath = kInSIM;
        sim.outBase    = kOutSIMBase;
      }

      if (doPP)
      {
        datasets.emplace_back();
        Dataset& dat = datasets.back();
        dat.label      = "DATA";
        dat.isSim      = false;
        dat.trigger    = kTriggerPP;
        dat.topDirName = kTriggerPP;
        dat.inFilePath = kInPP;
        dat.outBase    = kOutPPBase;
      }


    // Open files/top directories and missing logs
    for (auto& ds : datasets)
    {
      cout << ANSI_BOLD_CYN << "\n[OPEN] " << ds.label << " file: " << ds.inFilePath << ANSI_RESET << "\n";
      ds.file = TFile::Open(ds.inFilePath.c_str(), "READ");
      if (!ds.file || ds.file->IsZombie())
      {
        cerr << ANSI_BOLD_RED << "[FATAL] Cannot open input file: " << ds.inFilePath << ANSI_RESET << "\n";
        return 2;
      }

      ds.topDir = ds.file->GetDirectory(ds.topDirName.c_str());
      if (!ds.topDir)
      {
        cerr << ANSI_BOLD_RED << "[FATAL] Missing top-level directory in file: " << ds.topDirName << ANSI_RESET << "\n";
        return 3;
      }

      EnsureDir(ds.outBase);
      const string missPath = JoinPath(ds.outBase, "missing_hists.txt");
      ds.missingOut.open(missPath.c_str());
      if (!ds.missingOut.is_open())
      {
        cerr << ANSI_BOLD_RED << "[FATAL] Cannot open missing_hists.txt: " << missPath << ANSI_RESET << "\n";
        return 4;
      }
      ds.missingOut << "# missing_hists.txt for dataset=" << ds.label
                    << " topDir=" << ds.topDirName
                    << " input=" << ds.inFilePath << "\n";
      ds.missingOut << "# Each line: <fullpath> :: <reason>\n\n";
    }

    // Sections 1–3
    for (auto& ds : datasets)
    {
      Section1_EventLevel(ds);
      PrintPreselectionFailTable(ds);
      Section3_GeneralIsoQA(ds);
    }

    // Leakage factors from SIM if available
    LeakageFactors lf;
    if (doSIM)
    {
      for (auto& ds : datasets)
      {
        if (ds.isSim)
        {
          ReadLeakageFactorsFromSIM(ds, lf);
          break;
        }
      }
    }

    // ABCD purity + SS for each dataset
    for (auto& ds : datasets)
    {
      if (ds.isSim)
      {
        Section4_ABCDPurityAndSS(ds, lf);
      }
      else
      {
        LeakageFactors lfForData;
        if (isSimAndDataPP && lf.available) lfForData = lf;
        else lfForData.available = false;
        Section4_ABCDPurityAndSS(ds, lfForData);
      }
    }

    // General jet QA + MatchQA (cache for future sections)
    map<string, MatchCache> matchCaches;
    for (auto& ds : datasets)
    {
      Section5_GeneralJetQA(ds);

      MatchCache mc;
      Section5C_MatchQA(ds, mc);
      matchCaches[ds.label] = mc;
    }

    // Close and summarize
    cout << ANSI_BOLD_CYN << "\n==============================\n"
         << "[DONE] Closing files and summarizing missing histograms\n"
         << "==============================" << ANSI_RESET << "\n";

    for (auto& ds : datasets)
    {
      if (ds.missingOut.is_open()) ds.missingOut.close();
      if (ds.file) ds.file->Close();

      cout << ANSI_BOLD_YEL
           << "[" << ds.label << "] outputs written under: " << ds.outBase << "\n"
           << "[" << ds.label << "] missing histogram entries logged: " << ds.missingCount << "\n"
           << ANSI_RESET;
    }

    cout << ANSI_BOLD_CYN << "\n[OK] AnalyzeRecoilJets completed (Sections 1–5C + ABCD/SS/leakage + GeneralJetQA).\n"
         << "To continue implementing Sections 5D–5I in the same maintainable style, ask: \"continue: implement 5D–5I\".\n"
         << ANSI_RESET;

    return 0;
  }
} // namespace ARJ

int AnalyzeRecoilJets()
{
  const int rc = ARJ::Run();
  if (rc != 0) gSystem->Exit(rc);
  return rc;
}

