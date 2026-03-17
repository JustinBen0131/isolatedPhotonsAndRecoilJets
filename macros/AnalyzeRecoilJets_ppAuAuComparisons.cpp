// =============================================================================
// OPTIONAL: PP vs Au+Au (gold-gold) photon-ID deliverables
//
// Deliverable Set A (SS vars): per pT × centrality
//   - inclusive (no isolation requirement)
//   - iso pass
//   - iso fail (nonIso)
//
// Deliverable Set B (Eiso total): per pT × centrality
//   - inclusive (no shower-shape / tight requirement)
//   - tight pass
//   - tight fail
//
// Output base:
//   kOutPPAuAuBase/
//     {noIsoRequired, isoPassSSplots, isoFailSSplots,
//      noSS_isoSpectra, tightIsoSpectra, nonTightIsoSpectra}/
//        <centLo>_<centHi>/
//          (SS folders only) <ssVar>/
//            table2x3_overlay_pp_vs_auau.png
//            overlay_pp_vs_auau_pT_<lo>_<hi>.png  (first 6 pT bins)
// =============================================================================
static TH1* GetTH1FromTopDir(TDirectory* topDir, const string& hname)
{
if (!topDir) return nullptr;
TObject* obj = topDir->Get(hname.c_str());
if (!obj) return nullptr;
return dynamic_cast<TH1*>(obj);
}

static void StyleOverlayHist(TH1* h, int color, int mstyle)
{
if (!h) return;
h->SetLineWidth(2);
h->SetLineColor(color);
h->SetMarkerStyle(mstyle);
h->SetMarkerSize(0.95);
h->SetMarkerColor(color);
}

static void DrawMissingPad(const string& titleLine)
{
TLatex t;
t.SetNDC(true);
t.SetTextFont(42);
t.SetTextAlign(22);
t.SetTextSize(0.080);
t.DrawLatex(0.50, 0.55, "MISSING");

t.SetTextSize(0.050);
t.DrawLatex(0.50, 0.42, titleLine.c_str());
}

static TH1* CloneNormalizeStyle(TH1* hIn,
                              const string& newName,
                              int color,
                              int mstyle)
{
if (!hIn) return nullptr;

TH1* h = CloneTH1(hIn, newName);
if (!h) return nullptr;

EnsureSumw2(h);

const double I = h->Integral(0, h->GetNbinsX() + 1);
if (I > 0.0) h->Scale(1.0 / I);

StyleOverlayHist(h, color, mstyle);
h->SetTitle("");

return h;
}

static void DrawOverlayPad_PPvsAuAu(TDirectory* ppTop,
                                  TDirectory* aaTop,
                                  const string& histBase,
                                  const string& centSuffix,
                                  const string& centLabel,
                                  const PtBin& pb,
                                  const string& xTitle,
                                  const string& topLeftTitle,
                                  bool forceIsoXRange,
                                  vector<TH1*>& keepAlive)
{
const string hPPName = histBase + pb.suffix;
const string hAAName = histBase + pb.suffix + centSuffix;

TH1* rawPP = GetTH1FromTopDir(ppTop, hPPName);
TH1* rawAA = GetTH1FromTopDir(aaTop, hAAName);

if (!rawPP || !rawAA)
{
  std::ostringstream s;
  s << "pT: " << pb.lo << "-" << pb.hi << "  " << centLabel;
  DrawMissingPad(s.str());
  return;
}

// ------------------------------------------------------------------
// Verbosity (PP vs AuAu): print raw + normalization diagnostics ONCE
// per (histBase, cent, pT bin), even though we draw it twice (table + per-bin).
// This is specifically to debug cases where the normalized AuAu curve looks ~0.
// ------------------------------------------------------------------
{
  static std::set<std::string> s_printed;
  const std::string key = TString::Format("PP_AuAu|%s|%s|%s", histBase.c_str(), centSuffix.c_str(), pb.folder.c_str()).Data();

  if (!s_printed.count(key))
  {
    s_printed.insert(key);

    auto integralAll = [&](TH1* h)->double
    {
      if (!h) return 0.0;
      return h->Integral(0, h->GetNbinsX() + 1);
    };

    auto integralInRange = [&](TH1* h, double xlo, double xhi)->double
    {
      if (!h) return 0.0;
      const int b1 = h->GetXaxis()->FindBin(xlo);
      const int b2 = h->GetXaxis()->FindBin(xhi);
      return h->Integral(b1, b2);
    };

    const double ppIall = integralAll(rawPP);
    const double aaIall = integralAll(rawAA);

    const double xloDbg = forceIsoXRange ? -2.0 : rawPP->GetXaxis()->GetXmin();
    const double xhiDbg = forceIsoXRange ?  6.0 : rawPP->GetXaxis()->GetXmax();

    const double ppIwin = integralInRange(rawPP, xloDbg, xhiDbg);
    const double aaIwin = integralInRange(rawAA, xloDbg, xhiDbg);

    cout << ANSI_BOLD_CYN
         << "\n[PP_AuAu DEBUG] " << histBase
         << " | " << centLabel
         << " | pT " << pb.lo << "-" << pb.hi << " (" << pb.folder << ")"
         << ANSI_RESET << "\n";

    cout << "  PP hist: " << hPPName << "\n"
         << "    Entries=" << std::fixed << std::setprecision(0) << rawPP->GetEntries()
         << "  Integral(all)=" << ppIall
         << "  Integral(win [" << xloDbg << "," << xhiDbg << "])=" << ppIwin
         << "  Mean=" << std::setprecision(6) << rawPP->GetMean()
         << "  RMS="  << rawPP->GetRMS()
         << "  MaxBinContent=" << std::setprecision(6) << rawPP->GetMaximum()
         << "\n";

    cout << "  AuAu hist: " << hAAName << "\n"
         << "    Entries=" << std::fixed << std::setprecision(0) << rawAA->GetEntries()
         << "  Integral(all)=" << aaIall
         << "  Integral(win [" << xloDbg << "," << xhiDbg << "])=" << aaIwin
         << "  Mean=" << std::setprecision(6) << rawAA->GetMean()
         << "  RMS="  << rawAA->GetRMS()
         << "  MaxBinContent=" << std::setprecision(6) << rawAA->GetMaximum()
         << "\n";

    const double ppScale = (ppIall > 0.0) ? (1.0/ppIall) : 0.0;
    const double aaScale = (aaIall > 0.0) ? (1.0/aaIall) : 0.0;

    cout << "  Normalization factors (unit area over ALL bins including under/overflow):\n"
         << "    PP  scale = " << std::setprecision(12) << ppScale << "\n"
         << "    AuAu scale = " << std::setprecision(12) << aaScale << "\n";

    // If the “looks like 0 everywhere” effect is happening, it’s usually because:
    // - AuAu integral is enormous relative to the peak region (tiny per-bin probabilities)
    // - or the histogram is extremely wide / has far tails
    // - or max after scaling is orders-of-magnitude smaller than PP’s max.
    const double ppMaxNorm = (ppIall > 0.0) ? (rawPP->GetMaximum() * ppScale) : 0.0;
    const double aaMaxNorm = (aaIall > 0.0) ? (rawAA->GetMaximum() * aaScale) : 0.0;

    cout << "  Peak heights after scaling:\n"
         << "    PP  max(norm) = " << std::setprecision(12) << ppMaxNorm << "\n"
         << "    AuAu max(norm) = " << std::setprecision(12) << aaMaxNorm << "\n";

    if (ppMaxNorm > 0.0)
    {
      cout << "  Ratio AuAu/PP of normalized peak = " << std::setprecision(12)
           << ((ppMaxNorm > 0.0) ? (aaMaxNorm/ppMaxNorm) : 0.0) << "\n";
    }
  }
}

TH1* hPP = CloneNormalizeStyle(rawPP,
  TString::Format("pp_%s_%s", histBase.c_str(), pb.folder.c_str()).Data(),
  kBlack, 20);

TH1* hAA = CloneNormalizeStyle(rawAA,
  TString::Format("aa_%s_%s%s", histBase.c_str(), pb.folder.c_str(), centSuffix.c_str()).Data(),
  kRed + 1, 24);

if (!hPP || !hAA)
{
  if (hPP) delete hPP;
  if (hAA) delete hAA;
  std::ostringstream s;
  s << "pT: " << pb.lo << "-" << pb.hi << "  " << centLabel;
  DrawMissingPad(s.str());
  return;
}

// Iso spectra: keep the same x-range as your iso tables (most informative region)
if (forceIsoXRange)
{
  hPP->GetXaxis()->SetRangeUser(-2.0, 6.0);
  hAA->GetXaxis()->SetRangeUser(-2.0, 6.0);
}

hPP->GetXaxis()->SetTitle(xTitle.c_str());
hPP->GetYaxis()->SetTitle("Normalized counts");

const double ymax = std::max(hPP->GetMaximum(), hAA->GetMaximum());
hPP->SetMaximum(ymax * 1.35);

hPP->Draw("E1");
hAA->Draw("E1 same");

  // Always label BOTH curves explicitly (PP vs Au+Au) for every overlay canvas/pad.
  // Move legend slightly left/down to avoid collisions with the pT/cent text.
  TLegend leg(0.52, 0.70, 0.90, 0.86);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.036);
  leg.AddEntry(hPP, "PP data", "ep");
  leg.AddEntry(hAA, "Au+Au data", "ep");
  leg.Draw();

  TLatex t;
  t.SetNDC(true);
  t.SetTextFont(42);

  // ------------------------------------------------------------------
  // For SS variables: use centered top title + #gamma-ID cut annotation + cut lines
  // For non-SS (e.g. Eiso): keep the legacy label scheme.
  // ------------------------------------------------------------------
  const bool isSS =
    (xTitle == "weta" || xTitle == "wphi" || xTitle == "et1" || xTitle == "e11e33" || xTitle == "e32e35");

  if (isSS)
  {
    // Centered title (var + cent + pT)
    t.SetTextAlign(22);
    t.SetTextSize(0.045);
    t.DrawLatex(0.50, 0.94,
      TString::Format("%s, %s, p_{T}^{#gamma} = %d-%d GeV",
                      xTitle.c_str(), centLabel.c_str(), pb.lo, pb.hi).Data());

    // #gamma-ID annotation (top-left)
    TLatex tcut;
    tcut.SetNDC(true);
    tcut.SetTextFont(42);
    tcut.SetTextAlign(13);
    tcut.SetTextSize(0.040);

    bool drawCuts = false;
    double cutLo = 0.0;
    double cutHi = 0.0;
    std::string cutText;

    if (xTitle == "e11e33")
    {
      cutText = "Tight #gamma-ID: 0.4 < #frac{E_{11}}{E_{33}} < 0.98";
      drawCuts = true;
      cutLo = 0.4;
      cutHi = 0.98;
    }
    else if (xTitle == "e32e35")
    {
      cutText = "#gamma-ID: 0.92 < #frac{E_{32}}{E_{35}} < 1.0";
      drawCuts = true;
      cutLo = 0.92;
      cutHi = 1.0;
    }
    else if (xTitle == "et1")
    {
      cutText = "#gamma-ID: 0.9 < et1 < 1.0";
      drawCuts = true;
      cutLo = 0.9;
      cutHi = 1.0;
    }
    else if (xTitle == "weta")
    {
      cutText = "#gamma-ID: 0 < w_{#eta}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
    }
    else if (xTitle == "wphi")
    {
      cutText = "#gamma-ID: 0 < w_{#phi}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
    }

    if (!cutText.empty())
    {
      tcut.DrawLatex(0.16, 0.86, cutText.c_str());
    }

    // Draw cut lines using pad y-range (robust, and works for table pads + single canvases)
    if (drawCuts)
    {
      gPad->Update();
      const double yMin = gPad->GetUymin();
      const double yMax = gPad->GetUymax();

      TLine* l1 = new TLine(cutLo, yMin, cutLo, yMax);
      l1->SetLineColor(kGreen + 2);
      l1->SetLineWidth(2);
      l1->SetLineStyle(2);
      l1->Draw("same");

      TLine* l2 = new TLine(cutHi, yMin, cutHi, yMax);
      l2->SetLineColor(kOrange + 7);
      l2->SetLineWidth(2);
      l2->SetLineStyle(2);
      l2->Draw("same");

      gPad->RedrawAxis();
    }
  }
  else
  {
    // Legacy label scheme (for non-SS, e.g. Eiso)
    t.SetTextAlign(13);
    t.SetTextSize(0.055);
    t.DrawLatex(0.16, 0.92, topLeftTitle.c_str());

    t.SetTextAlign(33);
    t.SetTextSize(0.055);
    t.DrawLatex(0.95, 0.92,
      TString::Format("p_{T}^{#gamma}: %d-%d GeV", pb.lo, pb.hi).Data());
    t.SetTextSize(0.048);
    t.DrawLatex(0.95, 0.84, centLabel.c_str());
  }

keepAlive.push_back(hPP);
keepAlive.push_back(hAA);
}

static void Make2x3Table_PPvsAuAu(TDirectory* ppTop,
                                TDirectory* aaTop,
                                const string& histBase,
                                const string& centSuffix,
                                const string& centLabel,
                                const string& outPng,
                                const string& xTitle,
                                const string& topLeftTitle,
                                bool forceIsoXRange)
{
const int nPads = std::min(6, kNPtBins);
if (nPads <= 0) return;

TCanvas c(
  TString::Format("c_ppauau_tbl_%s%s", histBase.c_str(), centSuffix.c_str()).Data(),
  "c_ppauau_tbl", 1500, 800
);
c.Divide(3,2, 0.001, 0.001);

const auto& bins = PtBins();

vector<TH1*> keepAlive;
keepAlive.reserve((std::size_t)nPads * 2);

for (int i = 0; i < nPads; ++i)
{
  c.cd(i+1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.10);
  gPad->SetLogy(false);

  const PtBin& pb = bins[i];

  DrawOverlayPad_PPvsAuAu(ppTop, aaTop, histBase, centSuffix, centLabel, pb,
                          xTitle, topLeftTitle, forceIsoXRange, keepAlive);
}

SaveCanvas(c, outPng);

for (TH1* h : keepAlive) delete h;
keepAlive.clear();
}

// =============================================================================
// FULL-RANGE UNNORMALIZED QA for a (pp, AuAu) histogram family:
//   - per-pT PNGs in:
//       <outDir>/AuAu_unNormalized/<pTbin>/AuAu_unNormalized_<histBase>_<pTbin>.png
//       <outDir>/pp_unNormalized/<pTbin>/pp_unNormalized_<histBase>_<pTbin>.png
//   - plus 2x3 tables (first 6 pT bins) in:
//       <outDir>/AuAu_unNormalized/table2x3_AuAu_unNormalized.png
//       <outDir>/pp_unNormalized/table2x3_pp_unNormalized.png
//
// NOTE:
//   - FULL RANGE ONLY: no SetRangeUser() anywhere; we UnZoom() explicitly.
// =============================================================================
// =============================================================================
// NEW: per-pT-bin PP+AuAu overlay across ALL centrality bins
//   Outputs (inside <outDir>/noSS_isoSpectra/<pb.folder>/):
//     ppAuAu_allCent_overlay_<histBase>_<pb.folder>.png        (unnormalized counts)
//     ppAuAu_allCent_overlay_norm_<histBase>_<pb.folder>.png   (unit-area normalized)
//
//   AuAu curves: one per centrality bin, distinct closed-circle colors
//   PP curve: closed blue circles (kBlue+1, marker 20), drawn last on top
// =============================================================================
static void MakePPAuAuAllCent_OverlayByCent(TDirectory* ppTop,
                                            TDirectory* aaTop,
                                            const string& outDir,
                                            const string& histBase,
                                            const PtBin& pb,
                                            const string& xTitle)
{
  const auto& centBins = CentBins();
  if (centBins.empty() || !aaTop) return;

  const string qaDirOverlay = JoinPath(JoinPath(outDir, "noSS_isoSpectra"), pb.folder);
  EnsureDir(qaDirOverlay);

  // Distinct colors for centrality bins (closed circles throughout)
  const int centColors[] = {
    kRed + 1, kOrange + 7, kGreen + 2, kMagenta + 1, kCyan + 2, kViolet + 1,
    kAzure + 2, kSpring + 5, kPink + 7, kTeal + 3, kGray + 2, kYellow + 2
  };
  const int nCentColors = (int)(sizeof(centColors) / sizeof(centColors[0]));

  // Fetch PP histogram once
  const string hPPName = histBase + pb.suffix;
  TH1* rawPP = GetTH1FromTopDir(ppTop, hPPName);

  // ---------------------------------------------------------------
  // Inner lambda: build and save one overlay canvas
  //   normalize=false -> raw counts
  //   normalize=true  -> unit-area shape
  // ---------------------------------------------------------------
  auto MakeOneCanvas = [&](bool normalize)
  {
    const string cName = TString::Format("c_allCent_ov%s_%s_%s",
      normalize ? "_norm" : "",
      histBase.c_str(), pb.folder.c_str()).Data();

    TCanvas cOv(cName.c_str(), cName.c_str(), 900, 700);
    ApplyCanvasMargins1D(cOv);
    cOv.SetLogy(false);

    double yMax = 0.0;
    TH1* hFirst = nullptr;

    TLegend leg(0.52, 0.52, 0.90, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextFont(42);
    leg.SetTextSize(0.030);

    vector<TH1*> keepAlive;
    keepAlive.reserve(centBins.size() + 1);

    // Draw AuAu per centrality
    for (int ic = 0; ic < (int)centBins.size(); ++ic)
    {
      const auto& cb = centBins[ic];
      const string hAAName = histBase + pb.suffix + cb.suffix;
      TH1* rawAA = GetTH1FromTopDir(aaTop, hAAName);
      if (!rawAA) continue;

      TH1* hAA = CloneTH1(rawAA,
        TString::Format("aa_allCent_ov%s_%s_%s%s",
          normalize ? "_norm" : "",
          histBase.c_str(), pb.folder.c_str(), cb.suffix.c_str()).Data());
      if (!hAA) continue;

      EnsureSumw2(hAA);
      hAA->GetXaxis()->UnZoom();
      hAA->SetTitle("");
      hAA->GetXaxis()->SetTitle(xTitle.c_str());
      hAA->GetYaxis()->SetTitle(normalize ? "Normalized counts" : "Counts");

      if (normalize)
      {
        const double integral = hAA->Integral(0, hAA->GetNbinsX() + 1);
        if (integral > 0.0) hAA->Scale(1.0 / integral);
      }

      const int col = centColors[ic % nCentColors];
      StyleOverlayHist(hAA, col, 20);

      yMax = std::max(yMax, hAA->GetMaximum());

      if (!hFirst)
      {
        hFirst = hAA;
        hFirst->Draw("E1");
      }
      else
      {
        hAA->Draw("E1 same");
      }

      leg.AddEntry(hAA,
        TString::Format("Au+Au %d-%d%%", cb.lo, cb.hi).Data(), "ep");

      keepAlive.push_back(hAA);
    }

    // Draw PP on top (closed blue circles)
    if (rawPP)
    {
      TH1* hPP = CloneTH1(rawPP,
        TString::Format("pp_allCent_ov%s_%s_%s",
          normalize ? "_norm" : "",
          histBase.c_str(), pb.folder.c_str()).Data());
      if (hPP)
      {
        EnsureSumw2(hPP);
        hPP->GetXaxis()->UnZoom();
        hPP->SetTitle("");
        hPP->GetXaxis()->SetTitle(xTitle.c_str());
        hPP->GetYaxis()->SetTitle(normalize ? "Normalized counts" : "Counts");

        if (normalize)
        {
          const double integral = hPP->Integral(0, hPP->GetNbinsX() + 1);
          if (integral > 0.0) hPP->Scale(1.0 / integral);
        }

        StyleOverlayHist(hPP, kBlue + 1, 20);

        yMax = std::max(yMax, hPP->GetMaximum());

        if (!hFirst)
        {
          hFirst = hPP;
          hFirst->Draw("E1");
        }
        else
        {
          hPP->Draw("E1 same");
        }

        leg.AddEntry(hPP, "PP data", "ep");
        keepAlive.push_back(hPP);
      }
    }

    if (hFirst)
    {
      hFirst->SetMaximum(yMax * 1.35);

      leg.Draw();

      TLatex tOv;
      tOv.SetNDC(true);
      tOv.SetTextFont(42);
      tOv.SetTextAlign(22);
      tOv.SetTextSize(0.042);
      tOv.DrawLatex(0.50, 0.94,
        TString::Format("E_{T}^{Iso} PP+AuAu%s, p_{T}^{#gamma} = %d-%d GeV, all cent",
                        normalize ? " (norm)" : "", pb.lo, pb.hi).Data());

      const string outPng = JoinPath(qaDirOverlay,
        TString::Format("ppAuAu_allCent_overlay%s_%s_%s.png",
          normalize ? "_norm" : "",
          histBase.c_str(), pb.folder.c_str()).Data());

      SaveCanvas(cOv, outPng);
    }

    for (TH1* h : keepAlive) delete h;
  };

  MakeOneCanvas(false); // unnormalized counts
  MakeOneCanvas(true);  // unit-area normalized
}

static void MakeUnnormalizedQA_PPvsAuAu(TDirectory* ppTop,
                                     TDirectory* aaTop,
                                     const string& outDir,
                                     const string& histBase,
                                     const string& centSuffix,
                                     const string& centLabel,
                                     const string& xTitle,
                                     const string& topLeftTitle)
{
const int nPads = std::min(6, kNPtBins);
if (nPads <= 0) return;

const auto& bins = PtBins();

const string qaBaseAA = JoinPath(outDir, "AuAu_unNormalized");
const string qaBasePP = JoinPath(outDir, "pp_unNormalized");
EnsureDir(qaBaseAA);
EnsureDir(qaBasePP);

// -------------------------
// 2x3 table: AuAu counts (FULL RANGE)
// -------------------------
{
    TCanvas cAAtbl(
      TString::Format("c_aa_unNorm_tbl_%s%s", histBase.c_str(), centSuffix.c_str()).Data(),
      "c_aa_unNorm_tbl", 1500, 800
    );
    cAAtbl.Divide(3,2, 0.001, 0.001);

    vector<TH1*> keepAliveAA;
    keepAliveAA.reserve((std::size_t)nPads);

    for (int i = 0; i < nPads; ++i)
    {
      const PtBin& pb = bins[i];
      const string hAAName = histBase + pb.suffix + centSuffix;
      TH1* rawAA = GetTH1FromTopDir(aaTop, hAAName);

      cAAtbl.cd(i+1);
      gPad->SetLeftMargin(0.14);
      gPad->SetRightMargin(0.05);
      gPad->SetBottomMargin(0.14);
      gPad->SetTopMargin(0.12);
      gPad->SetLogy(false);

      if (!rawAA)
      {
        std::ostringstream s;
        s << "pT: " << pb.lo << "-" << pb.hi << "  " << centLabel;
        DrawMissingPad(s.str());
        continue;
      }

      TH1* hAAc = CloneTH1(rawAA,
        TString::Format("aa_tbl_counts_%s_%s%s", histBase.c_str(), pb.folder.c_str(), centSuffix.c_str()).Data());
      if (!hAAc) continue;

      EnsureSumw2(hAAc);
      hAAc->GetXaxis()->UnZoom();
      hAAc->SetTitle("");
      hAAc->GetXaxis()->SetTitle(xTitle.c_str());
      hAAc->GetYaxis()->SetTitle("Counts");
      StyleOverlayHist(hAAc, kRed + 1, 24);

      const double ymax = hAAc->GetMaximum();
      hAAc->SetMaximum(ymax * 1.35);

      hAAc->Draw("E1");

      TLegend leg(0.52, 0.70, 0.90, 0.86);
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);
      leg.SetTextFont(42);
      leg.SetTextSize(0.036);
      leg.AddEntry(hAAc, "Au+Au", "ep");
      leg.Draw();

      // Centered pad title: what is plotted
      TLatex t;
      t.SetNDC(true);
      t.SetTextFont(42);
      t.SetTextAlign(22);
      t.SetTextSize(0.042);
      t.DrawLatex(0.50, 0.93,
        TString::Format("%s, %s, p_{T}^{#gamma} = %d-%d GeV",
                        xTitle.c_str(), centLabel.c_str(), pb.lo, pb.hi).Data());

      // SS cut text + cut lines (where applicable)
      TLatex tcut;
      tcut.SetNDC(true);
      tcut.SetTextFont(42);
      tcut.SetTextAlign(13);
      tcut.SetTextSize(0.038);

      bool drawCuts = false;
      double cutLo = 0.0;
      double cutHi = 0.0;
      std::string cutText;

      if (xTitle == "e11e33")
      {
        cutText = "Tight #gamma-ID: 0.4 < #frac{E_{11}}{E_{33}} < 0.98";
        drawCuts = true;
        cutLo = 0.4;
        cutHi = 0.98;
      }
      else if (xTitle == "e32e35")
      {
        cutText = "#gamma-ID: 0.92 < #frac{E_{32}}{E_{35}} < 1.0";
        drawCuts = true;
        cutLo = 0.92;
        cutHi = 1.0;
      }
      else if (xTitle == "et1")
      {
        cutText = "#gamma-ID: 0.9 < et1 < 1.0";
        drawCuts = true;
        cutLo = 0.9;
        cutHi = 1.0;
      }
      else if (xTitle == "weta")
      {
        cutText = "#gamma-ID: 0 < w_{#eta}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
      }
      else if (xTitle == "wphi")
      {
        cutText = "#gamma-ID: 0 < w_{#phi}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
      }

      if (!cutText.empty())
      {
        tcut.DrawLatex(0.16, 0.86, cutText.c_str());
      }

      if (drawCuts)
      {
        gPad->Update();
        const double yMin = gPad->GetUymin();
        const double yMax = gPad->GetUymax();

        TLine* l1 = new TLine(cutLo, yMin, cutLo, yMax);
        l1->SetLineColor(kGreen + 2);
        l1->SetLineWidth(2);
        l1->SetLineStyle(2);
        l1->Draw("same");

        TLine* l2 = new TLine(cutHi, yMin, cutHi, yMax);
        l2->SetLineColor(kOrange + 7);
        l2->SetLineWidth(2);
        l2->SetLineStyle(2);
        l2->Draw("same");

        gPad->RedrawAxis();
      }

      keepAliveAA.push_back(hAAc);
    }
    SaveCanvas(cAAtbl, JoinPath(qaBaseAA, "table2x3_AuAu_unNormalized.png"));

    for (TH1* h : keepAliveAA) delete h;
    keepAliveAA.clear();
  }
  // -------------------------
  // 2x3 table: PP counts (FULL RANGE)
  // -------------------------
  {
    TCanvas cPPtbl(
      TString::Format("c_pp_unNorm_tbl_%s", histBase.c_str()).Data(),
      "c_pp_unNorm_tbl", 1500, 800
    );
    cPPtbl.Divide(3,2, 0.001, 0.001);

    vector<TH1*> keepAlivePP;
    keepAlivePP.reserve((std::size_t)nPads);

    for (int i = 0; i < nPads; ++i)
    {
      const PtBin& pb = bins[i];
      const string hPPName = histBase + pb.suffix;
      TH1* rawPP = GetTH1FromTopDir(ppTop, hPPName);

      cPPtbl.cd(i+1);
      gPad->SetLeftMargin(0.14);
      gPad->SetRightMargin(0.05);
      gPad->SetBottomMargin(0.14);
      gPad->SetTopMargin(0.12);
      gPad->SetLogy(false);

      if (!rawPP)
      {
        std::ostringstream s;
        s << "pT: " << pb.lo << "-" << pb.hi;
        DrawMissingPad(s.str());
        continue;
      }

      TH1* hPPc = CloneTH1(rawPP,
        TString::Format("pp_tbl_counts_%s_%s", histBase.c_str(), pb.folder.c_str()).Data());
      if (!hPPc) continue;

      EnsureSumw2(hPPc);
      hPPc->GetXaxis()->UnZoom();
      hPPc->SetTitle("");
      hPPc->GetXaxis()->SetTitle(xTitle.c_str());
      hPPc->GetYaxis()->SetTitle("Counts");
      StyleOverlayHist(hPPc, kBlack, 20);

      const double ymax = hPPc->GetMaximum();
      hPPc->SetMaximum(ymax * 1.35);

      hPPc->Draw("E1");

      TLegend leg(0.52, 0.70, 0.90, 0.86);
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);
      leg.SetTextFont(42);
      leg.SetTextSize(0.036);
      leg.AddEntry(hPPc, "PP data (counts)", "ep");
      leg.Draw();

      // Centered pad title: what is plotted
      TLatex t;
      t.SetNDC(true);
      t.SetTextFont(42);
      t.SetTextAlign(22);
      t.SetTextSize(0.044);
      t.DrawLatex(0.50, 0.93,
          TString::Format("%s, PP, p_{T}^{#gamma} = %d-%d GeV, Run24pp",
                          xTitle.c_str(), pb.lo, pb.hi).Data());

      // For E_{T}^{iso} plots: draw the isolation cut line and print the cut value per pT bin
      if (xTitle.find("E_{T}^{iso}") != std::string::npos)
      {
          const double pTmid  = 0.5 * ((double)pb.lo + (double)pb.hi);
          const double cutIso = 1.08128 + 0.0299107 * pTmid;

          gPad->Update();
          const double yMinIso = gPad->GetUymin();
          const double yMaxIso = gPad->GetUymax();

          TLine* lIso = new TLine(cutIso, yMinIso, cutIso, yMaxIso);
          lIso->SetLineColor(kRed + 1);
          lIso->SetLineWidth(2);
          lIso->SetLineStyle(2);
          lIso->Draw("same");

          TLatex txIso;
          txIso.SetNDC(true);
          txIso.SetTextFont(42);
          txIso.SetTextAlign(31);
          txIso.SetTextSize(0.04);
          txIso.DrawLatex(0.89, 0.8,
            TString::Format("E_{T}^{iso} < 1.08128 + 0.0299107*p_{T}^{#gamma} #approx %.7f", cutIso).Data());

          gPad->RedrawAxis();
      }

      // SS cut text + cut lines (where applicable)
      TLatex tcut;
      tcut.SetNDC(true);
      tcut.SetTextFont(42);
      tcut.SetTextAlign(13);
      tcut.SetTextSize(0.038);

      bool drawCuts = false;
      double cutLo = 0.0;
      double cutHi = 0.0;
      std::string cutText;

      if (xTitle == "e11e33")
      {
        cutText = "Tight #gamma-ID: 0.4 < #frac{E_{11}}{E_{33}} < 0.98";
        drawCuts = true;
        cutLo = 0.4;
        cutHi = 0.98;
      }
      else if (xTitle == "e32e35")
      {
        cutText = "#gamma-ID: 0.92 < #frac{E_{32}}{E_{35}} < 1.0";
        drawCuts = true;
        cutLo = 0.92;
        cutHi = 1.0;
      }
      else if (xTitle == "et1")
      {
        cutText = "#gamma-ID: 0.9 < et1 < 1.0";
        drawCuts = true;
        cutLo = 0.9;
        cutHi = 1.0;
      }
      else if (xTitle == "weta")
      {
        cutText = "#gamma-ID: 0 < w_{#eta}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
      }
      else if (xTitle == "wphi")
      {
        cutText = "#gamma-ID: 0 < w_{#phi}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
      }

      if (!cutText.empty())
      {
        tcut.DrawLatex(0.16, 0.86, cutText.c_str());
      }

      if (drawCuts)
      {
        gPad->Update();
        const double yMin = gPad->GetUymin();
        const double yMax = gPad->GetUymax();

        TLine* l1 = new TLine(cutLo, yMin, cutLo, yMax);
        l1->SetLineColor(kGreen + 2);
        l1->SetLineWidth(2);
        l1->SetLineStyle(2);
        l1->Draw("same");

        TLine* l2 = new TLine(cutHi, yMin, cutHi, yMax);
        l2->SetLineColor(kOrange + 7);
        l2->SetLineWidth(2);
        l2->SetLineStyle(2);
        l2->Draw("same");

        gPad->RedrawAxis();
      }

      keepAlivePP.push_back(hPPc);
    }

    SaveCanvas(cPPtbl, JoinPath(qaBasePP, "table2x3_pp_unNormalized.png"));

    for (TH1* h : keepAlivePP) delete h;
    keepAlivePP.clear();
}

// -------------------------
// Per-pT PNGs (FULL RANGE ONLY; exactly one per pT bin)
// -------------------------
for (int i = 0; i < nPads; ++i)
{
  const PtBin& pb = bins[i];

  const string hPPName = histBase + pb.suffix;
  const string hAAName = histBase + pb.suffix + centSuffix;

  TH1* rawPP = GetTH1FromTopDir(ppTop, hPPName);
  TH1* rawAA = GetTH1FromTopDir(aaTop, hAAName);

  const string qaDirAA  = JoinPath(qaBaseAA, pb.folder);
  const string qaDirPP  = JoinPath(qaBasePP, pb.folder);
  EnsureDir(qaDirAA);
  EnsureDir(qaDirPP);

  if (rawAA)
  {
      // Baseline (no UE subtraction): from current aaTop
      TH1* hNoUE = CloneTH1(rawAA,
        TString::Format("aa_counts_noUE_%s_%s%s", histBase.c_str(), pb.folder.c_str(), centSuffix.c_str()).Data());

      // UE-subtracted: open kInAuAuGoldNew and fetch the same hist name
      static TFile* fUE = nullptr;
      static TDirectory* aaTopUE = nullptr;

      if (!fUE)
      {
        fUE = TFile::Open(kInAuAuGoldNew.c_str(), "READ");
        if (fUE && aaTop)
        {
          aaTopUE = fUE->GetDirectory(aaTop->GetName());
          if (!aaTopUE) aaTopUE = fUE;
        }
      }

      TH1* rawUE = nullptr;
      if (aaTopUE)
      {
        const string hUEName = histBase + pb.suffix + centSuffix;
        rawUE = GetTH1FromTopDir(aaTopUE, hUEName);
      }

      TH1* hUE = nullptr;
      if (rawUE)
      {
        hUE = CloneTH1(rawUE,
          TString::Format("aa_counts_UEsub_%s_%s%s", histBase.c_str(), pb.folder.c_str(), centSuffix.c_str()).Data());
      }

      if (hNoUE)
      {
        EnsureSumw2(hNoUE);
        if (hUE) EnsureSumw2(hUE);

        TCanvas cAA(
          TString::Format("c_aa_counts_%s_%s%s", histBase.c_str(), pb.folder.c_str(), centSuffix.c_str()).Data(),
          "c_aa_counts", 900, 700
        );
        ApplyCanvasMargins1D(cAA);
        cAA.SetLogy(false);

        hNoUE->GetXaxis()->UnZoom();
        hNoUE->SetTitle("");
        hNoUE->GetXaxis()->SetTitle(xTitle.c_str());
        hNoUE->GetYaxis()->SetTitle("Counts");

        // Match table styling: black = no UE subtraction, red = UE subtracted
        StyleOverlayHist(hNoUE, kBlack, 20);
        if (hUE) StyleOverlayHist(hUE, kRed + 1, 24);

        const double ymax = std::max(hNoUE->GetMaximum(), (hUE ? hUE->GetMaximum() : 0.0));
        hNoUE->SetMaximum(ymax * 1.35);
        if (hUE) hUE->SetMaximum(ymax * 1.35);

        hNoUE->Draw("E1");
        if (hUE) hUE->Draw("E1 same");

        TLegend leg(0.52, 0.70, 0.90, 0.86);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.SetTextFont(42);
        leg.SetTextSize(0.036);

        if (hUE)
        {
          leg.AddEntry(hNoUE, "no UE subtraction", "ep");
          leg.AddEntry(hUE,   "UE subtracted",    "ep");
        }
        else
        {
          leg.AddEntry(hNoUE, "Au+Au", "ep");
        }
        leg.Draw();

        // Single centered title only (no extra on-canvas printouts, no "(inclusive)" label).
        std::string centPretty = centLabel;
        if (centPretty.rfind("cent:", 0) == 0)
        {
            centPretty = "Cent = " + centPretty.substr(5);
            while (!centPretty.empty() && centPretty[0] == ' ') centPretty.erase(centPretty.begin());
        }

        TLatex t;
        t.SetNDC(true);
        t.SetTextFont(42);
        t.SetTextAlign(22);
        t.SetTextSize(0.045);

        t.DrawLatex(0.50, 0.94,
            TString::Format("E_{T}^{Iso}, run3auau, %s, p_{T}^{#gamma} = %d-%d GeV",
                            centPretty.c_str(), pb.lo, pb.hi).Data());

        SaveCanvas(cAA, JoinPath(qaDirAA,
            TString::Format("AuAu_unNormalized_%s_%s.png", histBase.c_str(), pb.folder.c_str()).Data()));

          // ---------------------------------------------------------------
          // NEW: PP + AuAu overlay per centrality per pT bin
          //   -> <outDir>/noSS_isoSpectra/<pb.folder>/ppAuAu_overlay_<histBase>_<pb.folder>.png
          // ---------------------------------------------------------------
          if (rawPP)
          {
            const string qaDirOverlay = JoinPath(JoinPath(outDir, "noSS_isoSpectra"), pb.folder);
            EnsureDir(qaDirOverlay);

            TH1* hPPov = CloneTH1(rawPP,
              TString::Format("pp_ov_%s_%s%s", histBase.c_str(), pb.folder.c_str(), centSuffix.c_str()).Data());

            if (hPPov)
            {
              EnsureSumw2(hPPov);

              TCanvas cOv(
                TString::Format("c_ov_%s_%s%s", histBase.c_str(), pb.folder.c_str(), centSuffix.c_str()).Data(),
                "c_ov", 900, 700
              );
              ApplyCanvasMargins1D(cOv);
              cOv.SetLogy(false);

              hNoUE->GetXaxis()->UnZoom();
              hPPov->GetXaxis()->UnZoom();

              hNoUE->SetTitle("");
              hPPov->SetTitle("");
              hNoUE->GetXaxis()->SetTitle(xTitle.c_str());
              hNoUE->GetYaxis()->SetTitle("Counts");

              StyleOverlayHist(hNoUE, kBlack,  20);
              StyleOverlayHist(hPPov, kBlue+1, 20);

              const double ymaxOv = std::max(hNoUE->GetMaximum(), hPPov->GetMaximum());
              hNoUE->SetMaximum(ymaxOv * 1.35);

              hNoUE->Draw("E1");
              hPPov->Draw("E1 same");

              TLegend legOv(0.52, 0.70, 0.90, 0.86);
              legOv.SetBorderSize(0);
              legOv.SetFillStyle(0);
              legOv.SetTextFont(42);
              legOv.SetTextSize(0.036);
              legOv.AddEntry(hNoUE, "Au+Au (no UE sub)", "ep");
              legOv.AddEntry(hPPov, "PP data",           "ep");
              legOv.Draw();

              std::string centPrettyOv = centLabel;
              if (centPrettyOv.rfind("cent:", 0) == 0)
              {
                centPrettyOv = "Cent = " + centPrettyOv.substr(5);
                while (!centPrettyOv.empty() && centPrettyOv[0] == ' ') centPrettyOv.erase(centPrettyOv.begin());
              }

              TLatex tOv;
              tOv.SetNDC(true);
              tOv.SetTextFont(42);
              tOv.SetTextAlign(22);
              tOv.SetTextSize(0.045);
              tOv.DrawLatex(0.50, 0.94,
              TString::Format("E_{T}^{Iso} PP+AuAu, %s, p_{T}^{#gamma} = %d-%d GeV",
                                centPrettyOv.c_str(), pb.lo, pb.hi).Data());

              SaveCanvas(cOv, JoinPath(qaDirOverlay,
                                TString::Format("ppAuAu_overlay_%s_%s.png", histBase.c_str(), pb.folder.c_str()).Data()));

              delete hPPov;
            }

            // ---------------------------------------------------------------
            // NEW: PP+AuAu overlay across ALL centrality bins (unnorm + norm)
            //   -> noSS_isoSpectra/<pb.folder>/ppAuAu_allCent_overlay_<histBase>_<pb.folder>.png
            //   -> noSS_isoSpectra/<pb.folder>/ppAuAu_allCent_overlay_norm_<histBase>_<pb.folder>.png
            // ---------------------------------------------------------------
            MakePPAuAuAllCent_OverlayByCent(ppTop, aaTop, outDir, histBase, pb, xTitle);
          }
          // ---------------------------------------------------------------

          delete hNoUE;
          if (hUE) delete hUE;
      }
      else
      {
        if (hNoUE) delete hNoUE;
        if (hUE) delete hUE;
      }
  }

  if (rawPP)
  {
    TH1* hPPc = CloneTH1(rawPP,
      TString::Format("pp_counts_%s_%s", histBase.c_str(), pb.folder.c_str()).Data());
    if (hPPc)
    {
      EnsureSumw2(hPPc);

      TCanvas cPP(
        TString::Format("c_pp_counts_%s_%s", histBase.c_str(), pb.folder.c_str()).Data(),
        "c_pp_counts", 900, 700
      );
      ApplyCanvasMargins1D(cPP);
      cPP.SetLogy(false);

      hPPc->GetXaxis()->UnZoom();
      hPPc->SetTitle("");
      hPPc->GetXaxis()->SetTitle(xTitle.c_str());
      hPPc->GetYaxis()->SetTitle("Counts");

      StyleOverlayHist(hPPc, kBlack, 20);
      hPPc->Draw("E1");

      TLegend leg(0.52, 0.70, 0.90, 0.86);
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);
      leg.SetTextFont(42);
      leg.SetTextSize(0.036);
      leg.AddEntry(hPPc, "PP data (counts)", "ep");
      leg.Draw();

      TLatex t;
      t.SetNDC(true);
      t.SetTextFont(42);

      t.SetTextAlign(13);
      t.SetTextSize(0.055);
      t.DrawLatex(0.16, 0.92, topLeftTitle.c_str());

      t.SetTextAlign(33);
      t.SetTextSize(0.055);
      t.DrawLatex(0.95, 0.92,
        TString::Format("p_{T}^{#gamma}: %d-%d GeV", pb.lo, pb.hi).Data());

      SaveCanvas(cPP, JoinPath(qaDirPP,
        TString::Format("pp_unNormalized_%s_%s.png", histBase.c_str(), pb.folder.c_str()).Data()));

      delete hPPc;
    }
  }
}
}

// =============================================================================
// NEW: 2x3 AuAu unnormalized (counts) table overlaying TWO AuAu inputs
//   - "no UE-sub"          : aaTop
//   - "with UE-sub nodes"  : aaTopNew
// Output:
//   <outDir>/AuAu_unNormalized/table2x3_AuAu_unNormalized_overlay_UEsub.png
// =============================================================================
static void Make2x3Table_AuAuUnNormalized_OverlayUEsub(TDirectory* aaTop,
                                                       TDirectory* aaTopNew,
                                                       const string& outDir,
                                                       const string& histBase,
                                                       const string& centSuffix,
                                                       const string& centLabel,
                                                       const string& xTitle)
{
  if (!aaTop || !aaTopNew) return;

  const int nPads = std::min(6, kNPtBins);
  if (nPads <= 0) return;

  const auto& bins = PtBins();

  const string qaBaseAA = JoinPath(outDir, "AuAu_unNormalized");
  EnsureDir(qaBaseAA);

  TCanvas cAAtbl(
    TString::Format("c_aa_unNorm_tbl_overlayUEsub_%s%s", histBase.c_str(), centSuffix.c_str()).Data(),
    "c_aa_unNorm_tbl_overlayUEsub", 1500, 800
  );
  cAAtbl.Divide(3,2, 0.001, 0.001);

    vector<TH1*> keepAliveAA;
    keepAliveAA.reserve((std::size_t)nPads * 2);

    vector<TLegend*> keepAliveLeg;
    keepAliveLeg.reserve((std::size_t)nPads);

    for (int i = 0; i < nPads; ++i)
    {
      const PtBin& pb = bins[i];
      const string hName = histBase + pb.suffix + centSuffix;

      TH1* rawOld = GetTH1FromTopDir(aaTop,    hName);
      TH1* rawNew = GetTH1FromTopDir(aaTopNew, hName);

      cAAtbl.cd(i+1);
      gPad->SetLeftMargin(0.14);
      gPad->SetRightMargin(0.05);
      gPad->SetBottomMargin(0.14);
      gPad->SetTopMargin(0.12);
      gPad->SetLogy(false);

      if (!rawOld || !rawNew)
      {
        std::ostringstream s;
        s << "pT: " << pb.lo << "-" << pb.hi << "  " << centLabel;
        DrawMissingPad(s.str());
        continue;
      }

      TH1* hOld = CloneTH1(rawOld,
        TString::Format("aa_tbl_counts_old_%s_%s%s", histBase.c_str(), pb.folder.c_str(), centSuffix.c_str()).Data());
      TH1* hNew = CloneTH1(rawNew,
        TString::Format("aa_tbl_counts_new_%s_%s%s", histBase.c_str(), pb.folder.c_str(), centSuffix.c_str()).Data());

      if (!hOld || !hNew)
      {
        if (hOld) delete hOld;
        if (hNew) delete hNew;
        std::ostringstream s;
        s << "pT: " << pb.lo << "-" << pb.hi << "  " << centLabel;
        DrawMissingPad(s.str());
        continue;
      }

      EnsureSumw2(hOld);
      EnsureSumw2(hNew);
      hOld->GetXaxis()->UnZoom();
      hNew->GetXaxis()->UnZoom();

      hOld->SetTitle("");
      hNew->SetTitle("");
      hOld->GetXaxis()->SetTitle(xTitle.c_str());
      hOld->GetYaxis()->SetTitle("Counts");

      StyleOverlayHist(hOld, kBlack,   20);
      StyleOverlayHist(hNew, kRed + 1, 24);

      const double ymax = std::max(hOld->GetMaximum(), hNew->GetMaximum());
      hOld->SetMaximum(ymax * 1.35);

      hOld->Draw("E1");
      hNew->Draw("E1 same");

      TLegend* leg = new TLegend(0.52, 0.68, 0.90, 0.86);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.032);
      leg->AddEntry(hOld, "no UE subtraction", "ep");
      leg->AddEntry(hNew, "UE subtracted", "ep");
      leg->Draw();

      keepAliveLeg.push_back(leg);

      TLatex t;
      t.SetNDC(true);
      t.SetTextFont(42);
      t.SetTextAlign(22);
      t.SetTextSize(0.042);
      t.DrawLatex(0.50, 0.93,
        TString::Format("%s, %s, p_{T}^{#gamma} = %d-%d GeV",
                        xTitle.c_str(), centLabel.c_str(), pb.lo, pb.hi).Data());

      gPad->RedrawAxis();

      keepAliveAA.push_back(hOld);
      keepAliveAA.push_back(hNew);
    }

    SaveCanvas(cAAtbl, JoinPath(qaBaseAA, "table2x3_AuAu_unNormalized_overlay_UEsub.png"));

    for (TLegend* l : keepAliveLeg) delete l;
    keepAliveLeg.clear();

    for (TH1* h : keepAliveAA) delete h;
    keepAliveAA.clear();
}

static void MakePerPtOverlays_PPvsAuAu(TDirectory* ppTop,
                                     TDirectory* aaTop,
                                     const string& histBase,
                                     const string& centSuffix,
                                     const string& centLabel,
                                     const string& outDir,
                                     const string& xTitle,
                                     const string& topLeftTitle,
                                     bool forceIsoXRange)
{
const int nPads = std::min(6, kNPtBins);
if (nPads <= 0) return;

const auto& bins = PtBins();

// For Eiso families, keep the existing behavior:
// - generate FULL-range unnormalized QA outputs (per-pT + 2x3 tables)
if (forceIsoXRange)
{
  MakeUnnormalizedQA_PPvsAuAu(ppTop, aaTop, outDir, histBase, centSuffix, centLabel, xTitle, topLeftTitle);
}

for (int i = 0; i < nPads; ++i)
{
  const PtBin& pb = bins[i];

  TCanvas c(
    TString::Format("c_ppauau_%s_%s%s", histBase.c_str(), pb.folder.c_str(), centSuffix.c_str()).Data(),
    "c_ppauau_one", 900, 700
  );
  ApplyCanvasMargins1D(c);
  c.SetLogy(false);

  vector<TH1*> keepAlive;
  keepAlive.reserve(2);

  DrawOverlayPad_PPvsAuAu(ppTop, aaTop, histBase, centSuffix, centLabel, pb,
                          xTitle, topLeftTitle, forceIsoXRange, keepAlive);

  SaveCanvas(c, JoinPath(outDir,
    TString::Format("overlay_pp_vs_auau_%s.png", pb.folder.c_str()).Data()));

  for (TH1* h : keepAlive) delete h;
  keepAlive.clear();
}
}

static void ProduceFamily_PPvsAuAu(TDirectory* ppTop,
                                 TDirectory* aaTop,
                                 const string& outDir,
                                 const string& histBase,
                                 const string& centSuffix,
                                 const string& centLabel,
                                 const string& xTitle,
                                 const string& topLeftTitle,
                                 bool forceIsoXRange)
{
EnsureDir(outDir);

const string tablePng = JoinPath(outDir, "table2x3_overlay_pp_vs_auau.png");

Make2x3Table_PPvsAuAu(ppTop, aaTop, histBase, centSuffix, centLabel, tablePng,
                      xTitle, topLeftTitle, forceIsoXRange);

MakePerPtOverlays_PPvsAuAu(ppTop, aaTop, histBase, centSuffix, centLabel, outDir,
                           xTitle, topLeftTitle, forceIsoXRange);

if (histBase == "h_Eiso")
{
  static std::set<string> s_donePeripheralVsPP;

  const auto& centBins = CentBins();
  if (centBins.empty()) return;

  const auto& ptBins = PtBins();
  if (ptBins.empty()) return;

  const CentBin& peripheralBin = centBins.back();

  const string peripheralDir = JoinPath(DirnameFromPath(outDir), peripheralBin.folder);
  EnsureDir(peripheralDir);

  const string peripheralTablePng = JoinPath(peripheralDir, "table3x3_overlay_peripheral_vs_pp.png");
  if (s_donePeripheralVsPP.count(peripheralTablePng)) return;
  s_donePeripheralVsPP.insert(peripheralTablePng);

  TCanvas cTbl("c_peripheral_vs_pp_iso_tbl", "c_peripheral_vs_pp_iso_tbl", 1500, 1200);
  cTbl.Divide(3, 3, 0.001, 0.001);

  vector<TH1*> keepAlive;
  keepAlive.reserve((std::size_t)std::min(9, kNPtBins) * 2);

  vector<TLegend*> keepAliveLeg;
  keepAliveLeg.reserve((std::size_t)std::min(9, kNPtBins));

  const int nPads = std::min(9, kNPtBins);

  for (int i = 0; i < nPads; ++i)
  {
    const PtBin& pb = ptBins[i];

    const string hPPName = histBase + pb.suffix;
    const string hAAName = histBase + pb.suffix + peripheralBin.suffix;

    TH1* rawPP = GetTH1FromTopDir(ppTop, hPPName);
    TH1* rawAA = GetTH1FromTopDir(aaTop, hAAName);

    cTbl.cd(i + 1);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.14);
    gPad->SetTopMargin(0.12);
    gPad->SetLogy(false);

    if (!rawPP || !rawAA)
    {
      DrawMissingPad(TString::Format("Peripheral vs pp, p_{T}^{#gamma} = %d-%d GeV", pb.lo, pb.hi).Data());
      continue;
    }

    TH1* hPP = CloneTH1(rawPP,
      TString::Format("h_peripheralVsPP_pp_%s", pb.folder.c_str()).Data());

    TH1* hAA = CloneTH1(rawAA,
      TString::Format("h_peripheralVsPP_auau_%s", pb.folder.c_str()).Data());

    if (!hPP || !hAA)
    {
      if (hPP) delete hPP;
      if (hAA) delete hAA;
      DrawMissingPad(TString::Format("Peripheral vs pp, p_{T}^{#gamma} = %d-%d GeV", pb.lo, pb.hi).Data());
      continue;
    }

    hPP->SetLineWidth(2);
    hPP->SetLineColor(kBlack);
    hPP->SetMarkerStyle(24);
    hPP->SetMarkerSize(0.85);
    hPP->SetMarkerColor(kBlack);
    hPP->SetFillStyle(0);

    hAA->SetLineWidth(2);
    hAA->SetLineColor(kRed + 1);
    hAA->SetMarkerStyle(20);
    hAA->SetMarkerSize(0.85);
    hAA->SetMarkerColor(kRed + 1);
    hAA->SetFillStyle(0);

    hPP->GetXaxis()->SetTitle(xTitle.c_str());
    hPP->GetYaxis()->SetTitle("Counts");
    hAA->GetXaxis()->SetTitle(xTitle.c_str());
    hAA->GetYaxis()->SetTitle("Counts");

    hPP->GetXaxis()->SetRangeUser(-2.0, 6.0);
    hAA->GetXaxis()->SetRangeUser(-2.0, 6.0);

    const double yMax = std::max(hPP->GetMaximum(), hAA->GetMaximum());
    hPP->SetMaximum((yMax > 0.0) ? (1.25 * yMax) : 1.0);
    hPP->SetMinimum(0.0);

    hPP->Draw("E1");
    hAA->Draw("E1 same");

    auto* leg = new TLegend(0.62, 0.72, 0.90, 0.86);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.040);
    leg->AddEntry(hPP, "Run24pp", "ep");
    leg->AddEntry(hAA, "Run3auau", "ep");
    leg->Draw();
    keepAliveLeg.push_back(leg);

    TLatex t;
    t.SetNDC(true);
    t.SetTextFont(42);
    t.SetTextAlign(23);
    t.SetTextSize(0.052);
    t.DrawLatex(0.50, 0.958,
      TString::Format("Peripheral vs pp for p_{T}^{#gamma} = %d - %d GeV", pb.lo, pb.hi).Data()
    );

    keepAlive.push_back(hPP);
    keepAlive.push_back(hAA);
  }

  SaveCanvas(cTbl, peripheralTablePng);

  for (TH1* h : keepAlive) delete h;
  for (TLegend* leg : keepAliveLeg) delete leg;
}
}

// =============================================================================
// Tight vs nonTight isolation overlays (PP and AuAu)
//
// Output base:
//   kOutPPAuAuBase/tight_nonTight_isoSpectraOverlay/
//     pp/
//       overlay_tight_vs_nonTight_pT_<lo>_<hi>.png
//     AuAu/<centLo>_<centHi>/
//       overlay_tight_vs_nonTight_pT_<lo>_<hi>.png
// =============================================================================
static void ProduceTightNonTightIsoOverlays(TDirectory* ppTop,
                                          TDirectory* aaTop,
                                          const string& outBase,
                                          const string& centFolder,
                                          const string& centSuffix,
                                          const string& centLabel)
{
const int nPads = std::min(6, kNPtBins);
if (nPads <= 0) return;

const auto& bins = PtBins();

// -------------------------
// PP (no centrality)
// -------------------------
  {
    const string outPP = JoinPath(outBase, "pp");
    EnsureDir(outPP);

    // -------------------------
    // NEW: 2x3 pT summary table (first 6 pT bins)
    // tight vs nonTight (normalized isolation spectra)
    // -------------------------
    {
      TCanvas cTbl("c_pp_tight_nonTight_tbl", "c_pp_tight_nonTight_tbl", 1500, 800);
      cTbl.Divide(3,2, 0.001, 0.001);

      vector<TH1*> keepAlive;
      keepAlive.reserve((std::size_t)nPads * 2);

      vector<TLegend*> keepAliveLeg;
      keepAliveLeg.reserve((std::size_t)nPads);

      for (int i = 0; i < nPads; ++i)
      {
        const PtBin& pb = bins[i];

        const string hT = string("h_Eiso_tight")    + pb.suffix;
        const string hN = string("h_Eiso_nonTight") + pb.suffix;

        TH1* rawT = GetTH1FromTopDir(ppTop, hT);
        TH1* rawN = GetTH1FromTopDir(ppTop, hN);

        cTbl.cd(i+1);
        gPad->SetLeftMargin(0.14);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.14);
        gPad->SetTopMargin(0.12);
        gPad->SetLogy(false);

        if (!rawT || !rawN)
        {
          std::ostringstream s;
          s << "pT: " << pb.lo << "-" << pb.hi;
          DrawMissingPad(s.str());
          continue;
        }

        TH1* ht = CloneNormalizeStyle(rawT,
          TString::Format("pp_tight_tbl_%s", pb.folder.c_str()).Data(),
          kBlack, 20);

        TH1* hn = CloneNormalizeStyle(rawN,
          TString::Format("pp_nontight_tbl_%s", pb.folder.c_str()).Data(),
          kRed + 1, 24);

        if (!ht || !hn)
        {
          if (ht) delete ht;
          if (hn) delete hn;
          std::ostringstream s;
          s << "pT: " << pb.lo << "-" << pb.hi;
          DrawMissingPad(s.str());
          continue;
        }

        ht->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
        ht->GetYaxis()->SetTitle("Normalized counts");

        ht->GetXaxis()->SetRangeUser(-2.0, 6.0);
        hn->GetXaxis()->SetRangeUser(-2.0, 6.0);

        const double ymax = std::max(ht->GetMaximum(), hn->GetMaximum());
        ht->SetMaximum(ymax * 1.35);

        ht->Draw("E1");
        hn->Draw("E1 same");

        TLegend* leg = new TLegend(0.52, 0.70, 0.90, 0.86);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.036);
        leg->AddEntry(ht, "run24pp data (tight)", "ep");
        leg->AddEntry(hn, "run24pp data (nonTight)", "ep");
        leg->Draw();
        keepAliveLeg.push_back(leg);

        TLatex t;
        t.SetNDC(true);
        t.SetTextFont(42);
        t.SetTextAlign(22);
        t.SetTextSize(0.050);
        t.DrawLatex(0.50, 0.93,
          TString::Format("tight vs nonTight, p_{T}^{#gamma} = %d-%d GeV", pb.lo, pb.hi).Data());

        gPad->RedrawAxis();

        keepAlive.push_back(ht);
        keepAlive.push_back(hn);
      }

      SaveCanvas(cTbl, JoinPath(outPP, "table2x3_overlay_tight_vs_nonTight.png"));

      for (TLegend* l : keepAliveLeg) delete l;
      keepAliveLeg.clear();

      for (TH1* h : keepAlive) delete h;
      keepAlive.clear();
    }

    for (int i = 0; i < nPads; ++i)
    {
      const PtBin& pb = bins[i];

    const string hT = string("h_Eiso_tight")    + pb.suffix;
    const string hN = string("h_Eiso_nonTight") + pb.suffix;

    TH1* rawT = GetTH1FromTopDir(ppTop, hT);
    TH1* rawN = GetTH1FromTopDir(ppTop, hN);

    if (!rawT || !rawN)
    {
      cout << ANSI_BOLD_YEL
           << "[PP_AuAu DEBUG][WARN] Missing PP tight/nonTight hists: " << hT << " or " << hN
           << " (" << pb.folder << ")"
           << ANSI_RESET << "\n";
      continue;
    }

    TH1* ht = CloneNormalizeStyle(rawT,
      TString::Format("pp_tight_%s", pb.folder.c_str()).Data(),
      kBlack, 20);

    TH1* hn = CloneNormalizeStyle(rawN,
      TString::Format("pp_nontight_%s", pb.folder.c_str()).Data(),
      kRed + 1, 24);

    if (!ht || !hn)
    {
      if (ht) delete ht;
      if (hn) delete hn;
      continue;
    }

    TCanvas c(
      TString::Format("c_pp_tight_nonTight_%s", pb.folder.c_str()).Data(),
      "c_pp_tight_nonTight", 900, 700
    );
    ApplyCanvasMargins1D(c);
    c.SetLogy(false);

    ht->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
    ht->GetYaxis()->SetTitle("Normalized counts");

    ht->GetXaxis()->SetRangeUser(-2.0, 6.0);
    hn->GetXaxis()->SetRangeUser(-2.0, 6.0);

    const double ymax = std::max(ht->GetMaximum(), hn->GetMaximum());
    ht->SetMaximum(ymax * 1.35);

    ht->Draw("E1");
    hn->Draw("E1 same");

    // Always label BOTH curves explicitly: dataset + tightness.
    TLegend leg(0.52, 0.70, 0.90, 0.86);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextFont(42);
    leg.SetTextSize(0.036);
    leg.AddEntry(ht, "run24pp data (tight)", "ep");
    leg.AddEntry(hn, "run24pp data (nonTight)", "ep");
    leg.Draw();

    TLatex t;
    t.SetNDC(true);
    t.SetTextFont(42);
    t.SetTextAlign(33);
    t.SetTextSize(0.055);
    t.DrawLatex(0.95, 0.92,
      TString::Format("p_{T}^{#gamma}: %d-%d GeV", pb.lo, pb.hi).Data());

    SaveCanvas(c, JoinPath(outPP,
      TString::Format("overlay_tight_vs_nonTight_%s.png", pb.folder.c_str()).Data()));

    delete ht;
    delete hn;
  }
}

// -------------------------
// AuAu (per centrality)
// -------------------------
{
  const string outAA = JoinPath(JoinPath(outBase, "AuAu"), centFolder);
  EnsureDir(JoinPath(outBase, "AuAu"));
  EnsureDir(outAA);

  for (int i = 0; i < nPads; ++i)
  {
    const PtBin& pb = bins[i];

    const string hT = string("h_Eiso_tight")    + pb.suffix + centSuffix;
    const string hN = string("h_Eiso_nonTight") + pb.suffix + centSuffix;

    TH1* rawT = GetTH1FromTopDir(aaTop, hT);
    TH1* rawN = GetTH1FromTopDir(aaTop, hN);

    if (!rawT || !rawN)
    {
      cout << ANSI_BOLD_YEL
           << "[PP_AuAu DEBUG][WARN] Missing AuAu tight/nonTight hists: " << hT << " or " << hN
           << " (" << centLabel << ", " << pb.folder << ")"
           << ANSI_RESET << "\n";
      continue;
    }

    TH1* ht = CloneNormalizeStyle(rawT,
      TString::Format("aa_tight_%s%s", pb.folder.c_str(), centSuffix.c_str()).Data(),
      kBlack, 20);

    TH1* hn = CloneNormalizeStyle(rawN,
      TString::Format("aa_nontight_%s%s", pb.folder.c_str(), centSuffix.c_str()).Data(),
      kRed + 1, 24);

    if (!ht || !hn)
    {
      if (ht) delete ht;
      if (hn) delete hn;
      continue;
    }

    TCanvas c(
      TString::Format("c_aa_tight_nonTight_%s%s", pb.folder.c_str(), centSuffix.c_str()).Data(),
      "c_aa_tight_nonTight", 900, 700
    );
    ApplyCanvasMargins1D(c);
    c.SetLogy(false);

    ht->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
    ht->GetYaxis()->SetTitle("Normalized counts");

    ht->GetXaxis()->SetRangeUser(-2.0, 6.0);
    hn->GetXaxis()->SetRangeUser(-2.0, 6.0);

    const double ymax = std::max(ht->GetMaximum(), hn->GetMaximum());
    ht->SetMaximum(ymax * 1.35);

    ht->Draw("E1");
    hn->Draw("E1 same");

      // Always label BOTH curves explicitly: dataset + tightness.
      TLegend leg(0.52, 0.70, 0.90, 0.86);
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);
      leg.SetTextFont(42);
      leg.SetTextSize(0.036);
      leg.AddEntry(ht, "Au+Au data (tight)", "ep");
      leg.AddEntry(hn, "Au+Au data (nonTight)", "ep");
      leg.Draw();

    TLatex t;
    t.SetNDC(true);
    t.SetTextFont(42);

    t.SetTextAlign(33);
    t.SetTextSize(0.055);
    t.DrawLatex(0.95, 0.92,
      TString::Format("p_{T}^{#gamma}: %d-%d GeV", pb.lo, pb.hi).Data());
    t.SetTextSize(0.048);
    t.DrawLatex(0.95, 0.84, centLabel.c_str());

    SaveCanvas(c, JoinPath(outAA,
      TString::Format("overlay_tight_vs_nonTight_%s.png", pb.folder.c_str()).Data()));

    delete ht;
    delete hn;
  }
}
}

// =============================================================================
// IsoFail vs IsoPass overlays for shower-shape (SS) spectra (PP and AuAu)
//
// Output base:
//   kOutPPAuAuBase/isoFail_isoPass_SSspectra/
//     pp/<ssVar>/
//       overlay_isoFail_vs_isoPass_pT_<lo>_<hi>.png
//     AuAu/<centLo>_<centHi>/<ssVar>/
//       overlay_isoFail_vs_isoPass_pT_<lo>_<hi>.png
//
// Uses:
//   isoPass:  h_ss_<var>_iso
//   isoFail:  h_ss_<var>_nonIso
// =============================================================================
static void ProduceIsoFailIsoPassSSOverlays(TDirectory* ppTop,
                                          TDirectory* aaTop,
                                          const string& outBase,
                                          const string& centFolder,
                                          const string& centSuffix,
                                          const string& centLabel,
                                          const vector<string>& ssVars)
{
const int nPads = std::min(6, kNPtBins);
if (nPads <= 0) return;

const auto& bins = PtBins();

  // -------------------------
  // PP (no centrality)
  // -------------------------
  {
    const string outPP = JoinPath(outBase, "pp");
    EnsureDir(outPP);

    for (const auto& var : ssVars)
    {
      const string outPPvar = JoinPath(outPP, var);
      EnsureDir(outPPvar);

      // -------------------------
      // NEW: 2x3 pT summary table (first 6 pT bins)
      // -------------------------
      {
        TCanvas cTbl(
          TString::Format("c_pp_isoFail_isoPass_tbl_%s", var.c_str()).Data(),
          "c_pp_isoFail_isoPass_tbl", 1500, 800
        );
        cTbl.Divide(3,2, 0.001, 0.001);

          vector<TH1*> keepAlive;
          keepAlive.reserve((std::size_t)nPads * 2);

          vector<TLegend*> keepAliveLeg;
          keepAliveLeg.reserve((std::size_t)nPads);

          for (int i = 0; i < nPads; ++i)
          {
            const PtBin& pb = bins[i];

            const string hPass = string("h_ss_") + var + string("_iso")    + pb.suffix;
            const string hFail = string("h_ss_") + var + string("_nonIso") + pb.suffix;

            TH1* rawPass = GetTH1FromTopDir(ppTop, hPass);
            TH1* rawFail = GetTH1FromTopDir(ppTop, hFail);

            cTbl.cd(i+1);
            gPad->SetLeftMargin(0.14);
            gPad->SetRightMargin(0.05);
            gPad->SetBottomMargin(0.14);
            gPad->SetTopMargin(0.12);
            gPad->SetLogy(false);

            if (!rawPass || !rawFail)
            {
              std::ostringstream s;
              s << "pT: " << pb.lo << "-" << pb.hi;
              DrawMissingPad(s.str());
              continue;
            }

            TH1* hP = CloneNormalizeStyle(rawPass,
              TString::Format("pp_isoPass_tbl_%s_%s", var.c_str(), pb.folder.c_str()).Data(),
              kBlack, 20);

            TH1* hF = CloneNormalizeStyle(rawFail,
              TString::Format("pp_isoFail_tbl_%s_%s", var.c_str(), pb.folder.c_str()).Data(),
              kRed + 1, 24);

            if (!hP || !hF)
            {
              if (hP) delete hP;
              if (hF) delete hF;
              std::ostringstream s;
              s << "pT: " << pb.lo << "-" << pb.hi;
              DrawMissingPad(s.str());
              continue;
            }

              hP->GetXaxis()->SetTitle(var.c_str());
              hP->GetYaxis()->SetTitle("Normalized counts");

              if (var == "e32e35")
              {
                hP->GetXaxis()->SetRangeUser(0.8, 1.2);
                hF->GetXaxis()->SetRangeUser(0.8, 1.2);
              }

              const double ymax = std::max(hP->GetMaximum(), hF->GetMaximum());
              hP->SetMaximum(ymax * 1.35);

            hP->Draw("E1");
            hF->Draw("E1 same");

            // Legend (top-right): isolated vs non-isolated (keep alive until SaveCanvas)
            TLegend* leg = new TLegend(0.58, 0.70, 0.92, 0.86);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextFont(42);
            leg->SetTextSize(0.036);
            leg->AddEntry(hP, "isolated", "ep");
            leg->AddEntry(hF, "non-isolated", "ep");
            leg->Draw();

            keepAliveLeg.push_back(leg);

            // Centered title above each pad
            TLatex t;
            t.SetNDC(true);
            t.SetTextFont(42);
            t.SetTextAlign(22);
            t.SetTextSize(0.042);
            t.DrawLatex(0.50, 0.93,
              TString::Format("%s, PP, p_{T}^{#gamma} = %d-%d GeV",
                              var.c_str(), pb.lo, pb.hi).Data());

            // SS cut label + cut lines (where applicable)
            TLatex tcut;
            tcut.SetNDC(true);
            tcut.SetTextFont(42);
            tcut.SetTextAlign(13);
            tcut.SetTextSize(0.038);

            bool drawCuts = false;
            double cutLo = 0.0;
            double cutHi = 0.0;
            std::string cutText;

            if (var == "e11e33")
            {
              cutText = "Tight #gamma-ID: 0.4 < #frac{E_{11}}{E_{33}} < 0.98";
              drawCuts = true;
              cutLo = 0.4;
              cutHi = 0.98;
            }
            else if (var == "e32e35")
            {
              cutText = "#gamma-ID: 0.92 < #frac{E_{32}}{E_{35}} < 1.0";
              drawCuts = true;
              cutLo = 0.92;
              cutHi = 1.0;
            }
            else if (var == "et1")
            {
              cutText = "#gamma-ID: 0.9 < et1 < 1.0";
              drawCuts = true;
              cutLo = 0.9;
              cutHi = 1.0;
            }
            else if (var == "weta")
            {
              cutText = "#gamma-ID: 0 < w_{#eta}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
            }
            else if (var == "wphi")
            {
              cutText = "#gamma-ID: 0 < w_{#phi}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
            }

            if (!cutText.empty())
            {
              tcut.DrawLatex(0.16, 0.86, cutText.c_str());
            }

            if (drawCuts)
            {
              gPad->Update();
              const double yMin = gPad->GetUymin();
              const double yMax = gPad->GetUymax();

              TLine* l1 = new TLine(cutLo, yMin, cutLo, yMax);
              l1->SetLineColor(kGreen + 2);
              l1->SetLineWidth(2);
              l1->SetLineStyle(2);
              l1->Draw("same");

              TLine* l2 = new TLine(cutHi, yMin, cutHi, yMax);
              l2->SetLineColor(kOrange + 7);
              l2->SetLineWidth(2);
              l2->SetLineStyle(2);
              l2->Draw("same");

              gPad->RedrawAxis();
            }

            keepAlive.push_back(hP);
            keepAlive.push_back(hF);
          }

          SaveCanvas(cTbl, JoinPath(outPPvar, "table2x3_isoFail_vs_isoPass.png"));

          for (TLegend* l : keepAliveLeg) delete l;
          keepAliveLeg.clear();

          for (TH1* h : keepAlive) delete h;
          keepAlive.clear();
      }

      // -------------------------
      // Individual pT-bin overlays (first 6 pT bins)
      // -------------------------
      for (int i = 0; i < nPads; ++i)
      {
        const PtBin& pb = bins[i];

        const string hPass = string("h_ss_") + var + string("_iso")    + pb.suffix;
        const string hFail = string("h_ss_") + var + string("_nonIso") + pb.suffix;

        TH1* rawPass = GetTH1FromTopDir(ppTop, hPass);
        TH1* rawFail = GetTH1FromTopDir(ppTop, hFail);

        if (!rawPass || !rawFail)
        {
          cout << ANSI_BOLD_YEL
               << "[PP_AuAu DEBUG][WARN] Missing PP isoPass/isoFail SS hists: " << hPass << " or " << hFail
               << " (" << var << ", " << pb.folder << ")"
               << ANSI_RESET << "\n";
          continue;
        }

        TH1* hP = CloneNormalizeStyle(rawPass,
          TString::Format("pp_isoPass_%s_%s", var.c_str(), pb.folder.c_str()).Data(),
          kBlack, 20);

        TH1* hF = CloneNormalizeStyle(rawFail,
          TString::Format("pp_isoFail_%s_%s", var.c_str(), pb.folder.c_str()).Data(),
          kRed + 1, 24);

        if (!hP || !hF)
        {
          if (hP) delete hP;
          if (hF) delete hF;
          continue;
        }

        TCanvas c(
          TString::Format("c_pp_isoFail_isoPass_%s_%s", var.c_str(), pb.folder.c_str()).Data(),
          "c_pp_isoFail_isoPass", 900, 700
        );
        ApplyCanvasMargins1D(c);
        c.SetLogy(false);

        hP->GetXaxis()->SetTitle(var.c_str());
        hP->GetYaxis()->SetTitle("Normalized counts");

        const double ymax = std::max(hP->GetMaximum(), hF->GetMaximum());
        hP->SetMaximum(ymax * 1.35);

        hP->Draw("E1");
        hF->Draw("E1 same");

        // Legend (top-right): isolated vs non-isolated
        TLegend leg(0.58, 0.70, 0.92, 0.86);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.SetTextFont(42);
        leg.SetTextSize(0.036);
        leg.AddEntry(hP, "isolated", "ep");
        leg.AddEntry(hF, "non-isolated", "ep");
        leg.Draw();

        // Centered title above canvas
        TLatex t;
        t.SetNDC(true);
        t.SetTextFont(42);
        t.SetTextAlign(22);
        t.SetTextSize(0.045);
        t.DrawLatex(0.50, 0.94,
          TString::Format("%s, PP, p_{T}^{#gamma} = %d-%d GeV",
                          var.c_str(), pb.lo, pb.hi).Data());

        // SS cut label + cut lines (where applicable)
        TLatex tcut;
        tcut.SetNDC(true);
        tcut.SetTextFont(42);
        tcut.SetTextAlign(13);
        tcut.SetTextSize(0.040);

        bool drawCuts = false;
        double cutLo = 0.0;
        double cutHi = 0.0;
        std::string cutText;

        if (var == "e11e33")
        {
          cutText = "Tight #gamma-ID: 0.4 < #frac{E_{11}}{E_{33}} < 0.98";
          drawCuts = true;
          cutLo = 0.4;
          cutHi = 0.98;
        }
        else if (var == "e32e35")
        {
          cutText = "#gamma-ID: 0.92 < #frac{E_{32}}{E_{35}} < 1.0";
          drawCuts = true;
          cutLo = 0.92;
          cutHi = 1.0;
        }
        else if (var == "et1")
        {
          cutText = "#gamma-ID: 0.9 < et1 < 1.0";
          drawCuts = true;
          cutLo = 0.9;
          cutHi = 1.0;
        }
        else if (var == "weta")
        {
          cutText = "#gamma-ID: 0 < w_{#eta}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
        }
        else if (var == "wphi")
        {
          cutText = "#gamma-ID: 0 < w_{#phi}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
        }

        if (!cutText.empty())
        {
          tcut.DrawLatex(0.16, 0.86, cutText.c_str());
        }

        if (drawCuts)
        {
          gPad->Update();
          const double yMin = gPad->GetUymin();
          const double yMax = gPad->GetUymax();

          TLine* l1 = new TLine(cutLo, yMin, cutLo, yMax);
          l1->SetLineColor(kGreen + 2);
          l1->SetLineWidth(2);
          l1->SetLineStyle(2);
          l1->Draw("same");

          TLine* l2 = new TLine(cutHi, yMin, cutHi, yMax);
          l2->SetLineColor(kOrange + 7);
          l2->SetLineWidth(2);
          l2->SetLineStyle(2);
          l2->Draw("same");

          gPad->RedrawAxis();
        }

        SaveCanvas(c, JoinPath(outPPvar,
          TString::Format("overlay_isoFail_vs_isoPass_%s.png", pb.folder.c_str()).Data()));

        delete hP;
        delete hF;
      }
    }
  }
  
  // -------------------------
  // AuAu (per centrality)
  // -------------------------
  {
    const string outAA = JoinPath(JoinPath(outBase, "AuAu"), centFolder);
    EnsureDir(JoinPath(outBase, "AuAu"));
    EnsureDir(outAA);

    for (const auto& var : ssVars)
    {
      const string outAAvar = JoinPath(outAA, var);
      EnsureDir(outAAvar);

      // -------------------------
      // NEW: 2x3 pT summary table (first 6 pT bins)
      // -------------------------
      {
        TCanvas cTbl(
          TString::Format("c_aa_isoFail_isoPass_tbl_%s%s", var.c_str(), centSuffix.c_str()).Data(),
          "c_aa_isoFail_isoPass_tbl", 1500, 800
        );
        cTbl.Divide(3,2, 0.001, 0.001);

        vector<TH1*> keepAlive;
        keepAlive.reserve((std::size_t)nPads * 2);

        vector<TLegend*> keepAliveLeg;
        keepAliveLeg.reserve((std::size_t)nPads);

        for (int i = 0; i < nPads; ++i)
        {
            const PtBin& pb = bins[i];

            const string hPass = string("h_ss_") + var + string("_iso")    + pb.suffix + centSuffix;
            const string hFail = string("h_ss_") + var + string("_nonIso") + pb.suffix + centSuffix;

            TH1* rawPass = GetTH1FromTopDir(aaTop, hPass);
            TH1* rawFail = GetTH1FromTopDir(aaTop, hFail);

            cTbl.cd(i+1);
            gPad->SetLeftMargin(0.14);
            gPad->SetRightMargin(0.05);
            gPad->SetBottomMargin(0.14);
            gPad->SetTopMargin(0.12);
            gPad->SetLogy(false);

            if (!rawPass || !rawFail)
            {
              std::ostringstream s;
              s << "pT: " << pb.lo << "-" << pb.hi << "  " << centLabel;
              DrawMissingPad(s.str());
              continue;
            }

            TH1* hP = CloneNormalizeStyle(rawPass,
              TString::Format("aa_isoPass_tbl_%s_%s%s", var.c_str(), pb.folder.c_str(), centSuffix.c_str()).Data(),
              kBlack, 20);

            TH1* hF = CloneNormalizeStyle(rawFail,
              TString::Format("aa_isoFail_tbl_%s_%s%s", var.c_str(), pb.folder.c_str(), centSuffix.c_str()).Data(),
              kRed + 1, 24);

            if (!hP || !hF)
            {
              if (hP) delete hP;
              if (hF) delete hF;
              std::ostringstream s;
              s << "pT: " << pb.lo << "-" << pb.hi << "  " << centLabel;
              DrawMissingPad(s.str());
              continue;
            }

            hP->GetXaxis()->SetTitle(var.c_str());
            hP->GetYaxis()->SetTitle("Normalized counts");

            const double ymax = std::max(hP->GetMaximum(), hF->GetMaximum());
            hP->SetMaximum(ymax * 1.35);

            hP->Draw("E1");
            hF->Draw("E1 same");

            // Legend (top-right): isolated vs non-isolated (keep alive until SaveCanvas)
            TLegend* leg = new TLegend(0.58, 0.70, 0.92, 0.86);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextFont(42);
            leg->SetTextSize(0.036);
            leg->AddEntry(hP, "isolated", "ep");
            leg->AddEntry(hF, "non-isolated", "ep");
            leg->Draw();

            keepAliveLeg.push_back(leg);

            // Centered title above each pad
            TLatex t;
            t.SetNDC(true);
            t.SetTextFont(42);
            t.SetTextAlign(22);
            t.SetTextSize(0.042);
            t.DrawLatex(0.50, 0.93,
              TString::Format("%s, %s, p_{T}^{#gamma} = %d-%d GeV",
                              var.c_str(), centLabel.c_str(), pb.lo, pb.hi).Data());

            // SS cut label + cut lines (where applicable)
            TLatex tcut;
            tcut.SetNDC(true);
            tcut.SetTextFont(42);
            tcut.SetTextAlign(13);
            tcut.SetTextSize(0.038);

            bool drawCuts = false;
            double cutLo = 0.0;
            double cutHi = 0.0;
            std::string cutText;

            if (var == "e11e33")
            {
              cutText = "Tight #gamma-ID: 0.4 < #frac{E_{11}}{E_{33}} < 0.98";
              drawCuts = true;
              cutLo = 0.4;
              cutHi = 0.98;
            }
            else if (var == "e32e35")
            {
              cutText = "#gamma-ID: 0.92 < #frac{E_{32}}{E_{35}} < 1.0";
              drawCuts = true;
              cutLo = 0.92;
              cutHi = 1.0;
            }
            else if (var == "et1")
            {
              cutText = "#gamma-ID: 0.9 < et1 < 1.0";
              drawCuts = true;
              cutLo = 0.9;
              cutHi = 1.0;
            }
            else if (var == "weta")
            {
              cutText = "#gamma-ID: 0 < w_{#eta}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
            }
            else if (var == "wphi")
            {
              cutText = "#gamma-ID: 0 < w_{#phi}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
            }

            if (!cutText.empty())
            {
              tcut.DrawLatex(0.16, 0.86, cutText.c_str());
            }

            if (drawCuts)
            {
              gPad->Update();
              const double yMin = gPad->GetUymin();
              const double yMax = gPad->GetUymax();

              TLine* l1 = new TLine(cutLo, yMin, cutLo, yMax);
              l1->SetLineColor(kGreen + 2);
              l1->SetLineWidth(2);
              l1->SetLineStyle(2);
              l1->Draw("same");

              TLine* l2 = new TLine(cutHi, yMin, cutHi, yMax);
              l2->SetLineColor(kOrange + 7);
              l2->SetLineWidth(2);
              l2->SetLineStyle(2);
              l2->Draw("same");

              gPad->RedrawAxis();
            }

            keepAlive.push_back(hP);
            keepAlive.push_back(hF);
          }

          SaveCanvas(cTbl, JoinPath(outAAvar, "table2x3_isoFail_vs_isoPass.png"));

          for (TLegend* l : keepAliveLeg) delete l;
          keepAliveLeg.clear();

          for (TH1* h : keepAlive) delete h;
          keepAlive.clear();
      }

      // -------------------------
      // Individual pT-bin overlays (first 6 pT bins)
      // -------------------------
      for (int i = 0; i < nPads; ++i)
      {
        const PtBin& pb = bins[i];

        const string hPass = string("h_ss_") + var + string("_iso")    + pb.suffix + centSuffix;
        const string hFail = string("h_ss_") + var + string("_nonIso") + pb.suffix + centSuffix;

        TH1* rawPass = GetTH1FromTopDir(aaTop, hPass);
        TH1* rawFail = GetTH1FromTopDir(aaTop, hFail);

        if (!rawPass || !rawFail)
        {
          cout << ANSI_BOLD_YEL
               << "[PP_AuAu DEBUG][WARN] Missing AuAu isoPass/isoFail SS hists: " << hPass << " or " << hFail
               << " (" << centLabel << ", " << var << ", " << pb.folder << ")"
               << ANSI_RESET << "\n";
          continue;
        }

        TH1* hP = CloneNormalizeStyle(rawPass,
          TString::Format("aa_isoPass_%s_%s%s", var.c_str(), pb.folder.c_str(), centSuffix.c_str()).Data(),
          kBlack, 20);

        TH1* hF = CloneNormalizeStyle(rawFail,
          TString::Format("aa_isoFail_%s_%s%s", var.c_str(), pb.folder.c_str(), centSuffix.c_str()).Data(),
          kRed + 1, 24);

        if (!hP || !hF)
        {
          if (hP) delete hP;
          if (hF) delete hF;
          continue;
        }

        TCanvas c(
          TString::Format("c_aa_isoFail_isoPass_%s_%s%s", var.c_str(), pb.folder.c_str(), centSuffix.c_str()).Data(),
          "c_aa_isoFail_isoPass", 900, 700
        );
        ApplyCanvasMargins1D(c);
        c.SetLogy(false);

        hP->GetXaxis()->SetTitle(var.c_str());
        hP->GetYaxis()->SetTitle("Normalized counts");

        const double ymax = std::max(hP->GetMaximum(), hF->GetMaximum());
        hP->SetMaximum(ymax * 1.35);

        hP->Draw("E1");
        hF->Draw("E1 same");

        // Legend (top-right): isolated vs non-isolated
        TLegend leg(0.58, 0.70, 0.92, 0.86);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.SetTextFont(42);
        leg.SetTextSize(0.036);
        leg.AddEntry(hP, "isolated", "ep");
        leg.AddEntry(hF, "non-isolated", "ep");
        leg.Draw();

        // Centered title above canvas
        TLatex t;
        t.SetNDC(true);
        t.SetTextFont(42);
        t.SetTextAlign(22);
        t.SetTextSize(0.045);
        t.DrawLatex(0.50, 0.94,
          TString::Format("%s, %s, p_{T}^{#gamma} = %d-%d GeV",
                          var.c_str(), centLabel.c_str(), pb.lo, pb.hi).Data());

        // SS cut label + cut lines (where applicable)
        TLatex tcut;
        tcut.SetNDC(true);
        tcut.SetTextFont(42);
        tcut.SetTextAlign(13);
        tcut.SetTextSize(0.040);

        bool drawCuts = false;
        double cutLo = 0.0;
        double cutHi = 0.0;
        std::string cutText;

        if (var == "e11e33")
        {
          cutText = "Tight #gamma-ID: 0.4 < #frac{E_{11}}{E_{33}} < 0.98";
          drawCuts = true;
          cutLo = 0.4;
          cutHi = 0.98;
        }
        else if (var == "e32e35")
        {
          cutText = "#gamma-ID: 0.92 < #frac{E_{32}}{E_{35}} < 1.0";
          drawCuts = true;
          cutLo = 0.92;
          cutHi = 1.0;
        }
        else if (var == "et1")
        {
          cutText = "#gamma-ID: 0.9 < et1 < 1.0";
          drawCuts = true;
          cutLo = 0.9;
          cutHi = 1.0;
        }
        else if (var == "weta")
        {
          cutText = "#gamma-ID: 0 < w_{#eta}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
        }
        else if (var == "wphi")
        {
          cutText = "#gamma-ID: 0 < w_{#phi}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
        }

        if (!cutText.empty())
        {
          tcut.DrawLatex(0.16, 0.86, cutText.c_str());
        }

        if (drawCuts)
        {
          gPad->Update();
          const double yMin = gPad->GetUymin();
          const double yMax = gPad->GetUymax();

          TLine* l1 = new TLine(cutLo, yMin, cutLo, yMax);
          l1->SetLineColor(kGreen + 2);
          l1->SetLineWidth(2);
          l1->SetLineStyle(2);
          l1->Draw("same");

          TLine* l2 = new TLine(cutHi, yMin, cutHi, yMax);
          l2->SetLineColor(kOrange + 7);
          l2->SetLineWidth(2);
          l2->SetLineStyle(2);
          l2->Draw("same");

          gPad->RedrawAxis();
        }

        SaveCanvas(c, JoinPath(outAAvar,
          TString::Format("overlay_isoFail_vs_isoPass_%s.png", pb.folder.c_str()).Data()));

        delete hP;
        delete hF;
      }
    }
  }
}

// =============================================================================
// NEW: PP-only unnormalized ZOOMED 2x3 table for e32e35 (per centrality folder)
//
// Output:
//   <outDir>/pp_unNormalized/table2x3_pp_unNormalized_zoom.png
//
// Notes:
//   - Uses PP hist only (same regardless of centrality, but written under each cent folder).
//   - Unnormalized counts (no scaling).
//   - Zoomed to peak region for e32e35.
// =============================================================================
static void MakePPUnNormalizedZoomTable_e32e35(TDirectory* ppTop,
                                            const string& outDir,
                                            const string& centLabel)
{
if (!ppTop) return;

const int nPads = std::min(6, kNPtBins);
if (nPads <= 0) return;

const auto& bins = PtBins();

const string ppUnDir = JoinPath(outDir, "pp_unNormalized");
EnsureDir(ppUnDir);

TCanvas c(
  TString::Format("c_pp_unNorm_zoom_e32e35_%s", centLabel.c_str()).Data(),
  "c_pp_unNorm_zoom_e32e35", 1500, 800
);
c.Divide(3,2, 0.001, 0.001);

vector<TH1*> keepAlive;
keepAlive.reserve((std::size_t)nPads);

for (int i = 0; i < nPads; ++i)
{
  const PtBin& pb = bins[i];

  // PP inclusive SS hist name for e32e35
  const string hPPName = string("h_ss_e32e35_inclusive") + pb.suffix;
  TH1* rawPP = GetTH1FromTopDir(ppTop, hPPName);

  c.cd(i+1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.12);
  gPad->SetLogy(false);

  if (!rawPP)
  {
    std::ostringstream s;
    s << "pT: " << pb.lo << "-" << pb.hi << "  " << centLabel;
    DrawMissingPad(s.str());
    continue;
  }

  TH1* h = CloneTH1(rawPP,
    TString::Format("pp_unNorm_zoom_e32e35_%s", pb.folder.c_str()).Data());
  if (!h) continue;

  EnsureSumw2(h);

  h->SetTitle("");
  h->GetXaxis()->SetTitle("e32e35");
  h->GetYaxis()->SetTitle("Counts");

  // Zoom to peak region for e32e35
  h->GetXaxis()->SetRangeUser(0.85, 1.05);

  StyleOverlayHist(h, kBlack, 20);

  const double ymax = h->GetMaximum();
  h->SetMaximum(ymax * 1.35);

  h->Draw("E1");

  // Legend (top-right)
  TLegend leg(0.60, 0.70, 0.92, 0.86);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.036);
  leg.AddEntry(h, "PP data (counts)", "ep");
  leg.Draw();

  // Centered title
  TLatex t;
  t.SetNDC(true);
  t.SetTextFont(42);
  t.SetTextAlign(22);
  t.SetTextSize(0.042);
  t.DrawLatex(0.50, 0.93,
    TString::Format("e32e35, PP (zoom), %s, p_{T}^{#gamma} = %d-%d GeV",
                    centLabel.c_str(), pb.lo, pb.hi).Data());

  // SS cut label (top-left)
  TLatex tcut;
  tcut.SetNDC(true);
  tcut.SetTextFont(42);
  tcut.SetTextAlign(13);
  tcut.SetTextSize(0.038);
  tcut.DrawLatex(0.16, 0.86, "#gamma-ID: 0.92 < #frac{E_{32}}{E_{35}} < 1.0");

  // Vertical cut lines at 0.92 and 1.0 (use pad y-range to guarantee visibility)
  gPad->Update();
  const double yMin = gPad->GetUymin();
  const double yMax = gPad->GetUymax();

  TLine* l1 = new TLine(0.92, yMin, 0.92, yMax);
  l1->SetLineColor(kGreen + 2);
  l1->SetLineWidth(2);
  l1->SetLineStyle(2);
  l1->Draw("same");

  TLine* l2 = new TLine(1.0, yMin, 1.0, yMax);
  l2->SetLineColor(kOrange + 7);
  l2->SetLineWidth(2);
  l2->SetLineStyle(2);
  l2->Draw("same");

  gPad->RedrawAxis();

  keepAlive.push_back(h);
}

SaveCanvas(c, JoinPath(ppUnDir, "table2x3_pp_unNormalized_zoom.png"));

for (TH1* h : keepAlive) delete h;
keepAlive.clear();
}

// =============================================================================
// NEW: AuAu unnormalized Eiso (counts) summary by centrality (first 6)
//
// Output:
//   <kOutPPAuAuBase>/noSS_isoSpectra/table2x3_AuAu_unNormalized.png
//
// For each pT bin (first 6 pads), overlay AuAu Eiso counts for each
// centrality bin (first 6) on the same pad.
// =============================================================================
static void Make2x3Table_AuAuUnNormalized_ByCent(TDirectory* aaTop,
                                               const string& outDir,
                                               const string& histBase,
                                               const string& xTitle)
{
if (!aaTop) return;

const auto& ptBins   = PtBins();
const auto& centBins = CentBins();

const int nPads  = std::min(6, kNPtBins);
const int nCents = std::min(6, (int)centBins.size());
if (nPads <= 0 || nCents <= 0) return;

EnsureDir(outDir);

const int colors[] = {
  kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 1,
  kViolet + 1, kAzure + 2, kSpring + 5, kPink + 7, kTeal + 3, kGray + 2
};
const int nColors = (int)(sizeof(colors)/sizeof(colors[0]));

TCanvas c(
  TString::Format("c_aa_unNorm_byCent_%s", histBase.c_str()).Data(),
  "c_aa_unNorm_byCent", 1500, 800
);
c.Divide(3,2, 0.001, 0.001);

vector<TH1*> keepAlive;
keepAlive.reserve((std::size_t)nPads * (std::size_t)nCents);

vector<TLegend*> keepAliveLeg;
keepAliveLeg.reserve((std::size_t)nPads);

for (int ipt = 0; ipt < nPads; ++ipt)
{
  const PtBin& pb = ptBins[ipt];

  c.cd(ipt+1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.12);
  gPad->SetLogy(false);

  vector<TH1*> histsPad;
  histsPad.reserve((std::size_t)nCents);

  vector<std::string> labelsPad;
  labelsPad.reserve((std::size_t)nCents);

  for (int ic = 0; ic < nCents; ++ic)
  {
    const auto& cb = centBins[ic];
    const string hAAName = histBase + pb.suffix + cb.suffix;

    TH1* rawAA = GetTH1FromTopDir(aaTop, hAAName);
    if (!rawAA) continue;

    TH1* hAAc = CloneTH1(rawAA,
      TString::Format("aa_unNorm_byCent_%s_%s%s",
                      histBase.c_str(), pb.folder.c_str(), cb.suffix.c_str()).Data());
    if (!hAAc) continue;

    EnsureSumw2(hAAc);
    hAAc->GetXaxis()->UnZoom();
    hAAc->SetTitle("");
    hAAc->GetXaxis()->SetTitle(xTitle.c_str());
    hAAc->GetYaxis()->SetTitle("Counts");

    const int col = colors[ic % nColors];
    StyleOverlayHist(hAAc, col, 20);
    hAAc->SetMarkerStyle(20);

    histsPad.push_back(hAAc);
    labelsPad.push_back(TString::Format("%d-%d%%", cb.lo, cb.hi).Data());
  }

  if (histsPad.empty())
  {
    std::ostringstream s;
    s << "pT: " << pb.lo << "-" << pb.hi << "  AuAu by cent";
    DrawMissingPad(s.str());
    continue;
  }

  double yMax = 0.0;
  for (TH1* h : histsPad) yMax = std::max(yMax, (double)h->GetMaximum());

  histsPad[0]->SetMaximum(yMax * 1.35);
  histsPad[0]->Draw("E1");

  for (std::size_t j = 1; j < histsPad.size(); ++j)
  {
    histsPad[j]->Draw("E1 same");
  }

  TLegend* leg = new TLegend(0.52, 0.60, 0.90, 0.86);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.030);

  for (std::size_t j = 0; j < histsPad.size(); ++j)
  {
    leg->AddEntry(histsPad[j], labelsPad[j].c_str(), "ep");
  }
  leg->Draw();
  keepAliveLeg.push_back(leg);

  TLatex t;
  t.SetNDC(true);
  t.SetTextFont(42);
  t.SetTextAlign(22);
  t.SetTextSize(0.042);
  t.DrawLatex(0.50, 0.93,
    TString::Format("Au+Au (counts), centrality overlays, p_{T}^{#gamma} = %d-%d GeV",
                    pb.lo, pb.hi).Data());

  gPad->RedrawAxis();

  for (TH1* h : histsPad) keepAlive.push_back(h);
}

SaveCanvas(c, JoinPath(outDir, "table2x3_AuAu_unNormalized.png"));

for (TLegend* l : keepAliveLeg) delete l;
keepAliveLeg.clear();

for (TH1* h : keepAlive) delete h;
keepAlive.clear();
}

// =============================================================================
// NEW: AuAu unnormalized Eiso (counts) summary by pT overlays (all pT bins)
//
// Output:
//   <kOutPPAuAuBase>/noSS_isoSpectra/table2x3_AuAu_unNormalized_byPtOverlays.png
//
// For each centrality bin (first 6 pads), overlay AuAu Eiso counts for each
// pT bin (ALL available pT bins) on the same pad.
// =============================================================================
static void Make2x3Table_AuAuUnNormalized_ByPtOverlaysPerCent(TDirectory* aaTop,
                                                            const string& outDir,
                                                            const string& histBase,
                                                            const string& xTitle)
{
if (!aaTop) return;

const auto& ptBins   = PtBins();
const auto& centBins = CentBins();

const int nPads  = std::min(6, (int)centBins.size());
const int nPt    = (int)ptBins.size();
if (nPads <= 0 || nPt <= 0) return;

EnsureDir(outDir);

const int colors[] = {
  kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 1,
  kViolet + 1, kAzure + 2, kSpring + 5, kPink + 7, kTeal + 3, kGray + 2
};
const int nColors = (int)(sizeof(colors)/sizeof(colors[0]));

TCanvas c(
  TString::Format("c_aa_unNorm_byPtOverlays_%s", histBase.c_str()).Data(),
  "c_aa_unNorm_byPtOverlays", 1500, 800
);
c.Divide(3,2, 0.001, 0.001);

vector<TH1*> keepAlive;
keepAlive.reserve((std::size_t)nPads * (std::size_t)nPt);

vector<TLegend*> keepAliveLeg;
keepAliveLeg.reserve((std::size_t)nPads);

for (int ic = 0; ic < nPads; ++ic)
{
  const auto& cb = centBins[ic];

  c.cd(ic+1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.12);
  gPad->SetLogy(false);

  vector<TH1*> histsPad;
  histsPad.reserve((std::size_t)nPt);

  vector<std::string> labelsPad;
  labelsPad.reserve((std::size_t)nPt);

  for (int ipt = 0; ipt < nPt; ++ipt)
  {
    const PtBin& pb = ptBins[ipt];
    const string hAAName = histBase + pb.suffix + cb.suffix;

    TH1* rawAA = GetTH1FromTopDir(aaTop, hAAName);
    if (!rawAA) continue;

    TH1* hAAc = CloneTH1(rawAA,
      TString::Format("aa_unNorm_byPt_%s_%s%s",
                      histBase.c_str(), pb.folder.c_str(), cb.suffix.c_str()).Data());
    if (!hAAc) continue;

    EnsureSumw2(hAAc);
    hAAc->GetXaxis()->UnZoom();
    hAAc->SetTitle("");
    hAAc->GetXaxis()->SetTitle(xTitle.c_str());
    hAAc->GetYaxis()->SetTitle("Counts");

    const int col = colors[ipt % nColors];
    StyleOverlayHist(hAAc, col, 20);
    hAAc->SetMarkerStyle(20);

    histsPad.push_back(hAAc);
    labelsPad.push_back(TString::Format("%d-%d GeV", pb.lo, pb.hi).Data());
  }

  if (histsPad.empty())
  {
    std::ostringstream s;
    s << "cent: " << cb.lo << "-" << cb.hi << "  AuAu by pT";
    DrawMissingPad(s.str());
    continue;
  }

  double yMax = 0.0;
  for (TH1* h : histsPad) yMax = std::max(yMax, (double)h->GetMaximum());

  histsPad[0]->SetMaximum(yMax * 1.35);
  histsPad[0]->Draw("E1");

  for (std::size_t j = 1; j < histsPad.size(); ++j)
  {
    histsPad[j]->Draw("E1 same");
  }

  TLegend* leg = new TLegend(0.52, 0.60, 0.90, 0.86);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.028);

  for (std::size_t j = 0; j < histsPad.size(); ++j)
  {
    leg->AddEntry(histsPad[j], labelsPad[j].c_str(), "ep");
  }
  leg->Draw();
  keepAliveLeg.push_back(leg);

    TLatex t;
    t.SetNDC(true);
    t.SetTextFont(42);
    t.SetTextAlign(22);
    t.SetTextSize(0.042);
    t.DrawLatex(0.50, 0.93,
      TString::Format("Au+Au (counts), p_{T}^{#gamma} overlays, cent = %d-%d%%",
                      cb.lo, cb.hi).Data());

    gPad->RedrawAxis();

    // Also write the exact same pad as a standalone PNG inside:
    //   <outDir>/<centLo>_<centHi>/AuAu_unNormalized/AuAu_unNormalized_byPtOverlays.png
    {
      const string outCent  = JoinPath(outDir, cb.folder);
      const string qaBaseAA = JoinPath(outCent, "AuAu_unNormalized");
      EnsureDir(qaBaseAA);

      TCanvas cOne(
        TString::Format("c_aa_unNorm_byPtOverlays_single_%s_%s", histBase.c_str(), cb.folder.c_str()).Data(),
        "c_aa_unNorm_byPtOverlays_single", 500, 400
      );
      cOne.cd();
      gPad->SetLeftMargin(0.14);
      gPad->SetRightMargin(0.05);
      gPad->SetBottomMargin(0.14);
      gPad->SetTopMargin(0.12);
      gPad->SetLogy(false);

      histsPad[0]->SetMaximum(yMax * 1.35);
      histsPad[0]->Draw("E1");

      for (std::size_t j = 1; j < histsPad.size(); ++j)
      {
        histsPad[j]->Draw("E1 same");
      }

      TLegend legOne(0.52, 0.60, 0.90, 0.86);
      legOne.SetBorderSize(0);
      legOne.SetFillStyle(0);
      legOne.SetTextFont(42);
      legOne.SetTextSize(0.028);

      for (std::size_t j = 0; j < histsPad.size(); ++j)
      {
        legOne.AddEntry(histsPad[j], labelsPad[j].c_str(), "ep");
      }
      legOne.Draw();

      TLatex tOne;
      tOne.SetNDC(true);
      tOne.SetTextFont(42);
      tOne.SetTextAlign(22);
      tOne.SetTextSize(0.042);
      tOne.DrawLatex(0.50, 0.93,
        TString::Format("Au+Au (counts), p_{T}^{#gamma} overlays, cent = %d-%d%%",
                        cb.lo, cb.hi).Data());

      gPad->RedrawAxis();

      SaveCanvas(cOne, JoinPath(qaBaseAA, "AuAu_unNormalized_byPtOverlays.png"));
    }

    for (TH1* h : histsPad) keepAlive.push_back(h);
}

SaveCanvas(c, JoinPath(outDir, "table2x3_AuAu_unNormalized_byPtOverlays.png"));

for (TLegend* l : keepAliveLeg) delete l;
keepAliveLeg.clear();

  for (TH1* h : keepAlive) delete h;
  keepAlive.clear();
}

// =============================================================================
// NEW: AuAu-only SS (counts) summary by centrality (first 6) — WITH SS cut label + vlines
//
// Output (inside noIsoRequired/<ssVar>/):
//   table2x3_AuAu_unNormalized.png
//
// For each pT bin (first 6 pads), overlay AuAu SS counts for each
// centrality bin (first 6) on the same pad.
// =============================================================================
static void Make2x3Table_AuAuUnNormalized_SS_ByCent(TDirectory* aaTop,
                                                   const string& outDir,
                                                   const string& histBase,
                                                   const string& ssVar)
{
  if (!aaTop) return;

  const auto& ptBins   = PtBins();
  const auto& centBins = CentBins();

  const int nPads  = std::min(6, kNPtBins);
  const int nCents = std::min(6, (int)centBins.size());
  if (nPads <= 0 || nCents <= 0) return;

  EnsureDir(outDir);

  const int colors[] = {
    kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 1,
    kViolet + 1, kAzure + 2, kSpring + 5, kPink + 7, kTeal + 3, kGray + 2
  };
  const int nColors = (int)(sizeof(colors)/sizeof(colors[0]));

  TCanvas c(
    TString::Format("c_aa_unNorm_byCent_SS_%s", histBase.c_str()).Data(),
    "c_aa_unNorm_byCent_SS", 1500, 800
  );
  c.Divide(3,2, 0.001, 0.001);

  vector<TH1*> keepAlive;
  keepAlive.reserve((std::size_t)nPads * (std::size_t)nCents);

  vector<TLegend*> keepAliveLeg;
  keepAliveLeg.reserve((std::size_t)nPads);

  for (int ipt = 0; ipt < nPads; ++ipt)
  {
    const PtBin& pb = ptBins[ipt];

    c.cd(ipt+1);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.14);
    gPad->SetTopMargin(0.12);
    gPad->SetLogy(false);

    vector<TH1*> histsPad;
    histsPad.reserve((std::size_t)nCents);

    vector<std::string> labelsPad;
    labelsPad.reserve((std::size_t)nCents);

    for (int ic = 0; ic < nCents; ++ic)
    {
      const auto& cb = centBins[ic];
      const string hAAName = histBase + pb.suffix + cb.suffix;

      TH1* rawAA = GetTH1FromTopDir(aaTop, hAAName);
      if (!rawAA) continue;

      TH1* hAAc = CloneTH1(rawAA,
        TString::Format("aa_unNorm_byCent_SS_%s_%s%s",
                        histBase.c_str(), pb.folder.c_str(), cb.suffix.c_str()).Data());
      if (!hAAc) continue;

      EnsureSumw2(hAAc);
      hAAc->GetXaxis()->UnZoom();
      hAAc->SetTitle("");
      hAAc->GetXaxis()->SetTitle(ssVar.c_str());
      hAAc->GetYaxis()->SetTitle("Counts");

      const int col = colors[ic % nColors];
      StyleOverlayHist(hAAc, col, 20);
      hAAc->SetMarkerStyle(20);

      histsPad.push_back(hAAc);
      labelsPad.push_back(TString::Format("%d-%d%%", cb.lo, cb.hi).Data());
    }

    if (histsPad.empty())
    {
      std::ostringstream s;
      s << "pT: " << pb.lo << "-" << pb.hi << "  AuAu by cent";
      DrawMissingPad(s.str());
      continue;
    }

    double yMax = 0.0;
    for (TH1* h : histsPad) yMax = std::max(yMax, (double)h->GetMaximum());

    histsPad[0]->SetMaximum(yMax * 1.35);
    histsPad[0]->Draw("E1");

    for (std::size_t j = 1; j < histsPad.size(); ++j)
    {
      histsPad[j]->Draw("E1 same");
    }

    TLegend* leg = new TLegend(0.52, 0.60, 0.90, 0.86);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.030);

    for (std::size_t j = 0; j < histsPad.size(); ++j)
    {
      leg->AddEntry(histsPad[j], labelsPad[j].c_str(), "ep");
    }
    leg->Draw();
    keepAliveLeg.push_back(leg);

    // Centered title (per pad)
    TLatex t;
    t.SetNDC(true);
    t.SetTextFont(42);
    t.SetTextAlign(22);
    t.SetTextSize(0.042);
    t.DrawLatex(0.50, 0.93,
      TString::Format("Au+Au (counts), %s, centrality overlays, p_{T}^{#gamma} = %d-%d GeV",
                      ssVar.c_str(), pb.lo, pb.hi).Data());

    // SS cut label (top-left) + cut lines (where applicable)
    TLatex tcut;
    tcut.SetNDC(true);
    tcut.SetTextFont(42);
    tcut.SetTextAlign(13);
    tcut.SetTextSize(0.038);

    bool drawCuts = false;
    double cutLo = 0.0;
    double cutHi = 0.0;
    std::string cutText;

    if (ssVar == "e11e33")
    {
      cutText = "Tight #gamma-ID: 0.4 < #frac{E_{11}}{E_{33}} < 0.98";
      drawCuts = true;
      cutLo = 0.4;
      cutHi = 0.98;
    }
    else if (ssVar == "e32e35")
    {
      cutText = "#gamma-ID: 0.92 < #frac{E_{32}}{E_{35}} < 1.0";
      drawCuts = true;
      cutLo = 0.92;
      cutHi = 1.0;
    }
    else if (ssVar == "et1")
    {
      cutText = "#gamma-ID: 0.9 < et1 < 1.0";
      drawCuts = true;
      cutLo = 0.9;
      cutHi = 1.0;
    }
    else if (ssVar == "weta")
    {
      cutText = "#gamma-ID: 0 < w_{#eta}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
    }
    else if (ssVar == "wphi")
    {
      cutText = "#gamma-ID: 0 < w_{#phi}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
    }

    if (!cutText.empty())
    {
      tcut.DrawLatex(0.16, 0.86, cutText.c_str());
    }

    if (drawCuts)
    {
      gPad->Update();
      const double yMin = gPad->GetUymin();
      const double yMaxPad = gPad->GetUymax();

      TLine* l1 = new TLine(cutLo, yMin, cutLo, yMaxPad);
      l1->SetLineColor(kGreen + 2);
      l1->SetLineWidth(2);
      l1->SetLineStyle(2);
      l1->Draw("same");

      TLine* l2 = new TLine(cutHi, yMin, cutHi, yMaxPad);
      l2->SetLineColor(kOrange + 7);
      l2->SetLineWidth(2);
      l2->SetLineStyle(2);
      l2->Draw("same");

      gPad->RedrawAxis();
    }

    gPad->RedrawAxis();

    for (TH1* h : histsPad) keepAlive.push_back(h);
  }

    SaveCanvas(c, JoinPath(outDir, "table2x3_AuAu_unNormalized.png"));

    // Zoomed tables for selected SS vars:
    //   e11e33, e32e35, et1 -> [0.9, 1.1]
    //   weta, wphi          -> [0.0, 0.1]
    bool doZoom = false;
    double zxLo = 0.0;
    double zxHi = 0.0;
    std::string zoomTag;

    if (ssVar == "e11e33" || ssVar == "e32e35" || ssVar == "et1")
    {
      doZoom = true;
      zxLo = 0.9;
      zxHi = 1.1;
      zoomTag = "0p9_1p1";
    }
    else if (ssVar == "weta" || ssVar == "wphi")
    {
      doZoom = true;
      zxLo = 0.0;
      zxHi = 0.1;
      zoomTag = "0p0_0p1";
    }

    if (doZoom)
    {
      TCanvas cz(
        TString::Format("c_aa_unNorm_byCent_SS_%s_zoom_%s", histBase.c_str(), zoomTag.c_str()).Data(),
        "c_aa_unNorm_byCent_SS_zoom", 1500, 800
      );
      cz.Divide(3,2, 0.001, 0.001);

      vector<TH1*> keepAliveZ;
      keepAliveZ.reserve((std::size_t)nPads * (std::size_t)nCents);

      vector<TLegend*> keepAliveLegZ;
      keepAliveLegZ.reserve((std::size_t)nPads);

      for (int ipt = 0; ipt < nPads; ++ipt)
      {
        const PtBin& pb = ptBins[ipt];

        cz.cd(ipt+1);
        gPad->SetLeftMargin(0.14);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.14);
        gPad->SetTopMargin(0.12);
        gPad->SetLogy(false);

        vector<TH1*> histsPad;
        histsPad.reserve((std::size_t)nCents);

        vector<std::string> labelsPad;
        labelsPad.reserve((std::size_t)nCents);

        for (int ic = 0; ic < nCents; ++ic)
        {
          const auto& cb = centBins[ic];
          const string hAAName = histBase + pb.suffix + cb.suffix;

          TH1* rawAA = GetTH1FromTopDir(aaTop, hAAName);
          if (!rawAA) continue;

          TH1* hAAc = CloneTH1(rawAA,
            TString::Format("aa_unNorm_byCent_SS_%s_%s%s_zoom_%s",
                            histBase.c_str(), pb.folder.c_str(), cb.suffix.c_str(), zoomTag.c_str()).Data());
          if (!hAAc) continue;

          EnsureSumw2(hAAc);
          hAAc->GetXaxis()->UnZoom();
          hAAc->SetTitle("");
          hAAc->GetXaxis()->SetTitle(ssVar.c_str());
          hAAc->GetYaxis()->SetTitle("Counts");
          hAAc->GetXaxis()->SetRangeUser(zxLo, zxHi);

          const int col = colors[ic % nColors];
          StyleOverlayHist(hAAc, col, 20);
          hAAc->SetMarkerStyle(20);

          histsPad.push_back(hAAc);
          labelsPad.push_back(TString::Format("%d-%d%%", cb.lo, cb.hi).Data());
        }

        if (histsPad.empty())
        {
          std::ostringstream s;
          s << "pT: " << pb.lo << "-" << pb.hi << "  AuAu by cent (zoom)";
          DrawMissingPad(s.str());
          continue;
        }

        double yMax = 0.0;
        for (TH1* h : histsPad) yMax = std::max(yMax, (double)h->GetMaximum());

        histsPad[0]->SetMaximum(yMax * 1.35);
        histsPad[0]->Draw("E1");

        for (std::size_t j = 1; j < histsPad.size(); ++j)
        {
          histsPad[j]->Draw("E1 same");
        }

        TLegend* leg = new TLegend(0.52, 0.60, 0.90, 0.86);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.030);

        for (std::size_t j = 0; j < histsPad.size(); ++j)
        {
          leg->AddEntry(histsPad[j], labelsPad[j].c_str(), "ep");
        }
        leg->Draw();
        keepAliveLegZ.push_back(leg);

        TLatex t;
        t.SetNDC(true);
        t.SetTextFont(42);
        t.SetTextAlign(22);
        t.SetTextSize(0.042);
        t.DrawLatex(0.50, 0.93,
          TString::Format("Au+Au (counts), %s, centrality overlays, p_{T}^{#gamma} = %d-%d GeV",
                          ssVar.c_str(), pb.lo, pb.hi).Data());

        TLatex tcut;
        tcut.SetNDC(true);
        tcut.SetTextFont(42);
        tcut.SetTextAlign(13);
        tcut.SetTextSize(0.038);

        bool drawCuts = false;
        double cutLo = 0.0;
        double cutHi = 0.0;
        std::string cutText;

        if (ssVar == "e11e33")
        {
          cutText = "Tight #gamma-ID: 0.4 < #frac{E_{11}}{E_{33}} < 0.98";
          drawCuts = true;
          cutLo = 0.4;
          cutHi = 0.98;
        }
        else if (ssVar == "e32e35")
        {
          cutText = "#gamma-ID: 0.92 < #frac{E_{32}}{E_{35}} < 1.0";
          drawCuts = true;
          cutLo = 0.92;
          cutHi = 1.0;
        }
        else if (ssVar == "et1")
        {
          cutText = "#gamma-ID: 0.9 < et1 < 1.0";
          drawCuts = true;
          cutLo = 0.9;
          cutHi = 1.0;
        }
        else if (ssVar == "weta")
        {
          cutText = "#gamma-ID: 0 < w_{#eta}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
        }
        else if (ssVar == "wphi")
        {
          cutText = "#gamma-ID: 0 < w_{#phi}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
        }

        if (!cutText.empty())
        {
          tcut.DrawLatex(0.16, 0.86, cutText.c_str());
        }

        if (drawCuts)
        {
          gPad->Update();
          const double yMin = gPad->GetUymin();
          const double yMaxPad = gPad->GetUymax();

          TLine* l1 = new TLine(cutLo, yMin, cutLo, yMaxPad);
          l1->SetLineColor(kGreen + 2);
          l1->SetLineWidth(2);
          l1->SetLineStyle(2);
          l1->Draw("same");

          TLine* l2 = new TLine(cutHi, yMin, cutHi, yMaxPad);
          l2->SetLineColor(kOrange + 7);
          l2->SetLineWidth(2);
          l2->SetLineStyle(2);
          l2->Draw("same");

          gPad->RedrawAxis();
        }

        gPad->RedrawAxis();

        for (TH1* h : histsPad) keepAliveZ.push_back(h);
      }

      SaveCanvas(cz, JoinPath(outDir,
        TString::Format("table2x3_AuAu_unNormalized_zoom_%s.png", zoomTag.c_str()).Data()));

      for (TLegend* l : keepAliveLegZ) delete l;
      keepAliveLegZ.clear();

      for (TH1* h : keepAliveZ) delete h;
      keepAliveZ.clear();
    }

    for (TLegend* l : keepAliveLeg) delete l;
    keepAliveLeg.clear();

    for (TH1* h : keepAlive) delete h;
    keepAlive.clear();
}

// =============================================================================
// NEW: AuAu-only SS (counts) summary by pT overlays (all pT bins) — WITH SS cut label + vlines
//
// Output (inside noIsoRequired/<ssVar>/):
//   table2x3_AuAu_unNormalized_byPtOverlays.png
//
// For each centrality bin (first 6 pads), overlay AuAu SS counts for each
// pT bin (ALL available pT bins) on the same pad.
// =============================================================================
static void Make2x3Table_AuAuUnNormalized_SS_ByPtOverlaysPerCent(TDirectory* aaTop,
                                                                 const string& outDir,
                                                                 const string& histBase,
                                                                 const string& ssVar)
{
  if (!aaTop) return;

  const auto& ptBins   = PtBins();
  const auto& centBins = CentBins();

  const int nPads  = std::min(6, (int)centBins.size());
  const int nPt    = (int)ptBins.size();
  if (nPads <= 0 || nPt <= 0) return;

  EnsureDir(outDir);

  const int colors[] = {
    kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 1,
    kViolet + 1, kAzure + 2, kSpring + 5, kPink + 7, kTeal + 3, kGray + 2
  };
  const int nColors = (int)(sizeof(colors)/sizeof(colors[0]));

  TCanvas c(
    TString::Format("c_aa_unNorm_byPtOverlays_SS_%s", histBase.c_str()).Data(),
    "c_aa_unNorm_byPtOverlays_SS", 1500, 800
  );
  c.Divide(3,2, 0.001, 0.001);

  vector<TH1*> keepAlive;
  keepAlive.reserve((std::size_t)nPads * (std::size_t)nPt);

  vector<TLegend*> keepAliveLeg;
  keepAliveLeg.reserve((std::size_t)nPads);

  for (int ic = 0; ic < nPads; ++ic)
  {
    const auto& cb = centBins[ic];

    c.cd(ic+1);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.14);
    gPad->SetTopMargin(0.12);
    gPad->SetLogy(false);

    vector<TH1*> histsPad;
    histsPad.reserve((std::size_t)nPt);

    vector<std::string> labelsPad;
    labelsPad.reserve((std::size_t)nPt);

    for (int ipt = 0; ipt < nPt; ++ipt)
    {
      const PtBin& pb = ptBins[ipt];
      const string hAAName = histBase + pb.suffix + cb.suffix;

      TH1* rawAA = GetTH1FromTopDir(aaTop, hAAName);
      if (!rawAA) continue;

      TH1* hAAc = CloneTH1(rawAA,
        TString::Format("aa_unNorm_byPt_SS_%s_%s%s",
                        histBase.c_str(), pb.folder.c_str(), cb.suffix.c_str()).Data());
      if (!hAAc) continue;

      EnsureSumw2(hAAc);
      hAAc->GetXaxis()->UnZoom();
      hAAc->SetTitle("");
      hAAc->GetXaxis()->SetTitle(ssVar.c_str());
      hAAc->GetYaxis()->SetTitle("Counts");

      const int col = colors[ipt % nColors];
      StyleOverlayHist(hAAc, col, 20);
      hAAc->SetMarkerStyle(20);

      histsPad.push_back(hAAc);
      labelsPad.push_back(TString::Format("%d-%d GeV", pb.lo, pb.hi).Data());
    }

    if (histsPad.empty())
    {
      std::ostringstream s;
      s << "cent: " << cb.lo << "-" << cb.hi << "  AuAu by pT";
      DrawMissingPad(s.str());
      continue;
    }

    double yMax = 0.0;
    for (TH1* h : histsPad) yMax = std::max(yMax, (double)h->GetMaximum());

    histsPad[0]->SetMaximum(yMax * 1.35);
    histsPad[0]->Draw("E1");

    for (std::size_t j = 1; j < histsPad.size(); ++j)
    {
      histsPad[j]->Draw("E1 same");
    }

    TLegend* leg = new TLegend(0.52, 0.60, 0.90, 0.86);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.028);

    for (std::size_t j = 0; j < histsPad.size(); ++j)
    {
      leg->AddEntry(histsPad[j], labelsPad[j].c_str(), "ep");
    }
    leg->Draw();
    keepAliveLeg.push_back(leg);

    // Centered title (per pad)
    TLatex t;
    t.SetNDC(true);
    t.SetTextFont(42);
    t.SetTextAlign(22);
    t.SetTextSize(0.042);
    t.DrawLatex(0.50, 0.93,
      TString::Format("Au+Au (counts), %s, p_{T}^{#gamma} overlays, cent = %d-%d%%",
                      ssVar.c_str(), cb.lo, cb.hi).Data());

    // SS cut label (top-left) + cut lines (where applicable)
    TLatex tcut;
    tcut.SetNDC(true);
    tcut.SetTextFont(42);
    tcut.SetTextAlign(13);
    tcut.SetTextSize(0.038);

    bool drawCuts = false;
    double cutLo = 0.0;
    double cutHi = 0.0;
    std::string cutText;

    if (ssVar == "e11e33")
    {
      cutText = "Tight #gamma-ID: 0.4 < #frac{E_{11}}{E_{33}} < 0.98";
      drawCuts = true;
      cutLo = 0.4;
      cutHi = 0.98;
    }
    else if (ssVar == "e32e35")
    {
      cutText = "#gamma-ID: 0.92 < #frac{E_{32}}{E_{35}} < 1.0";
      drawCuts = true;
      cutLo = 0.92;
      cutHi = 1.0;
    }
    else if (ssVar == "et1")
    {
      cutText = "#gamma-ID: 0.9 < et1 < 1.0";
      drawCuts = true;
      cutLo = 0.9;
      cutHi = 1.0;
    }
    else if (ssVar == "weta")
    {
      cutText = "#gamma-ID: 0 < w_{#eta}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
    }
    else if (ssVar == "wphi")
    {
      cutText = "#gamma-ID: 0 < w_{#phi}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
    }

    if (!cutText.empty())
    {
      tcut.DrawLatex(0.16, 0.86, cutText.c_str());
    }

    if (drawCuts)
    {
      gPad->Update();
      const double yMin = gPad->GetUymin();
      const double yMaxPad = gPad->GetUymax();

      TLine* l1 = new TLine(cutLo, yMin, cutLo, yMaxPad);
      l1->SetLineColor(kGreen + 2);
      l1->SetLineWidth(2);
      l1->SetLineStyle(2);
      l1->Draw("same");

      TLine* l2 = new TLine(cutHi, yMin, cutHi, yMaxPad);
      l2->SetLineColor(kOrange + 7);
      l2->SetLineWidth(2);
      l2->SetLineStyle(2);
      l2->Draw("same");

      gPad->RedrawAxis();
    }

    gPad->RedrawAxis();

    for (TH1* h : histsPad) keepAlive.push_back(h);
  }

  SaveCanvas(c, JoinPath(outDir, "table2x3_AuAu_unNormalized_byPtOverlays.png"));

  // Zoomed tables for selected SS vars:
  //   e11e33, e32e35, et1 -> [0.9, 1.1]
  //   weta, wphi          -> [0.0, 0.1]
  bool doZoom = false;
  double zxLo = 0.0;
  double zxHi = 0.0;
  std::string zoomTag;

  if (ssVar == "e11e33" || ssVar == "e32e35" || ssVar == "et1")
  {
      doZoom = true;
      zxLo = 0.9;
      zxHi = 1.1;
      zoomTag = "0p9_1p1";
  }
  else if (ssVar == "weta" || ssVar == "wphi")
    {
      doZoom = true;
      zxLo = 0.0;
      zxHi = 0.1;
      zoomTag = "0p0_0p1";
  }

  if (doZoom)
  {
      TCanvas cz(
        TString::Format("c_aa_unNorm_byPtOverlays_SS_%s_zoom_%s", histBase.c_str(), zoomTag.c_str()).Data(),
        "c_aa_unNorm_byPtOverlays_SS_zoom", 1500, 800
      );
      cz.Divide(3,2, 0.001, 0.001);

      vector<TH1*> keepAliveZ;
      keepAliveZ.reserve((std::size_t)nPads * (std::size_t)nPt);

      vector<TLegend*> keepAliveLegZ;
      keepAliveLegZ.reserve((std::size_t)nPads);

      for (int ic = 0; ic < nPads; ++ic)
      {
        const auto& cb = centBins[ic];

        cz.cd(ic+1);
        gPad->SetLeftMargin(0.14);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.14);
        gPad->SetTopMargin(0.12);
        gPad->SetLogy(false);

        vector<TH1*> histsPad;
        histsPad.reserve((std::size_t)nPt);

        vector<std::string> labelsPad;
        labelsPad.reserve((std::size_t)nPt);

        for (int ipt = 0; ipt < nPt; ++ipt)
        {
          const PtBin& pb = ptBins[ipt];
          const string hAAName = histBase + pb.suffix + cb.suffix;

          TH1* rawAA = GetTH1FromTopDir(aaTop, hAAName);
          if (!rawAA) continue;

          TH1* hAAc = CloneTH1(rawAA,
            TString::Format("aa_unNorm_byPt_SS_%s_%s%s_zoom_%s",
                            histBase.c_str(), pb.folder.c_str(), cb.suffix.c_str(), zoomTag.c_str()).Data());
          if (!hAAc) continue;

          EnsureSumw2(hAAc);
          hAAc->GetXaxis()->UnZoom();
          hAAc->SetTitle("");
          hAAc->GetXaxis()->SetTitle(ssVar.c_str());
          hAAc->GetYaxis()->SetTitle("Counts");
          hAAc->GetXaxis()->SetRangeUser(zxLo, zxHi);

          const int col = colors[ipt % nColors];
          StyleOverlayHist(hAAc, col, 20);
          hAAc->SetMarkerStyle(20);

          histsPad.push_back(hAAc);
          labelsPad.push_back(TString::Format("%d-%d GeV", pb.lo, pb.hi).Data());
        }

        if (histsPad.empty())
        {
          std::ostringstream s;
          s << "cent: " << cb.lo << "-" << cb.hi << "  AuAu by pT (zoom)";
          DrawMissingPad(s.str());
          continue;
        }

        double yMax = 0.0;
        for (TH1* h : histsPad) yMax = std::max(yMax, (double)h->GetMaximum());

        histsPad[0]->SetMaximum(yMax * 1.35);
        histsPad[0]->Draw("E1");

        for (std::size_t j = 1; j < histsPad.size(); ++j)
        {
          histsPad[j]->Draw("E1 same");
        }

        TLegend* leg = new TLegend(0.52, 0.60, 0.90, 0.86);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.028);

        for (std::size_t j = 0; j < histsPad.size(); ++j)
        {
          leg->AddEntry(histsPad[j], labelsPad[j].c_str(), "ep");
        }
        leg->Draw();
        keepAliveLegZ.push_back(leg);

        TLatex t;
        t.SetNDC(true);
        t.SetTextFont(42);
        t.SetTextAlign(22);
        t.SetTextSize(0.042);
        t.DrawLatex(0.50, 0.93,
          TString::Format("Au+Au (counts), %s, p_{T}^{#gamma} overlays, cent = %d-%d%%",
                          ssVar.c_str(), cb.lo, cb.hi).Data());

        TLatex tcut;
        tcut.SetNDC(true);
        tcut.SetTextFont(42);
        tcut.SetTextAlign(13);
        tcut.SetTextSize(0.038);

        bool drawCuts = false;
        double cutLo = 0.0;
        double cutHi = 0.0;
        std::string cutText;

        if (ssVar == "e11e33")
        {
          cutText = "Tight #gamma-ID: 0.4 < #frac{E_{11}}{E_{33}} < 0.98";
          drawCuts = true;
          cutLo = 0.4;
          cutHi = 0.98;
        }
        else if (ssVar == "e32e35")
        {
          cutText = "#gamma-ID: 0.92 < #frac{E_{32}}{E_{35}} < 1.0";
          drawCuts = true;
          cutLo = 0.92;
          cutHi = 1.0;
        }
        else if (ssVar == "et1")
        {
          cutText = "#gamma-ID: 0.9 < et1 < 1.0";
          drawCuts = true;
          cutLo = 0.9;
          cutHi = 1.0;
        }
        else if (ssVar == "weta")
        {
          cutText = "#gamma-ID: 0 < w_{#eta}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
        }
        else if (ssVar == "wphi")
        {
          cutText = "#gamma-ID: 0 < w_{#phi}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
        }

        if (!cutText.empty())
        {
          tcut.DrawLatex(0.16, 0.86, cutText.c_str());
        }

        if (drawCuts)
        {
          gPad->Update();
          const double yMin = gPad->GetUymin();
          const double yMaxPad = gPad->GetUymax();

          TLine* l1 = new TLine(cutLo, yMin, cutLo, yMaxPad);
          l1->SetLineColor(kGreen + 2);
          l1->SetLineWidth(2);
          l1->SetLineStyle(2);
          l1->Draw("same");

          TLine* l2 = new TLine(cutHi, yMin, cutHi, yMaxPad);
          l2->SetLineColor(kOrange + 7);
          l2->SetLineWidth(2);
          l2->SetLineStyle(2);
          l2->Draw("same");

          gPad->RedrawAxis();
        }

        gPad->RedrawAxis();

        for (TH1* h : histsPad) keepAliveZ.push_back(h);
      }

      SaveCanvas(cz, JoinPath(outDir,
        TString::Format("table2x3_AuAu_unNormalized_byPtOverlays_zoom_%s.png", zoomTag.c_str()).Data()));

      for (TLegend* l : keepAliveLegZ) delete l;
      keepAliveLegZ.clear();

      for (TH1* h : keepAliveZ) delete h;
      keepAliveZ.clear();
    }

    for (TLegend* l : keepAliveLeg) delete l;
    keepAliveLeg.clear();

    for (TH1* h : keepAlive) delete h;
    keepAlive.clear();
}

void RunPPvsAuAuDeliverables(Dataset& dsPP)
{
  cout << ANSI_BOLD_CYN << "\n[EXTRA] PP vs Au+Au (gold-gold) photon-ID deliverables\n" << ANSI_RESET;

// --- Open AuAu gold-gold file(s) ---
//   (1) kInAuAuGold    : "no UE-sub" (current)
//   (2) kInAuAuGoldNew : "with UE-sub nodes" (new)
TFile* fAA = TFile::Open(kInAuAuGold.c_str(), "READ");
if (!fAA || fAA->IsZombie())
  {
    cout << ANSI_BOLD_RED
         << "[ERROR] Cannot open AuAu gold file: " << kInAuAuGold
         << ANSI_RESET << "\n";
    if (fAA) { fAA->Close(); delete fAA; }
    return;
  }

  TDirectory* aaTop = fAA->GetDirectory(kTriggerAuAuGold.c_str());
  if (!aaTop)
  {
    cout << ANSI_BOLD_RED
         << "[ERROR] Missing AuAu trigger directory '" << kTriggerAuAuGold
         << "' in file: " << kInAuAuGold
         << ANSI_RESET << "\n";
    fAA->Close(); delete fAA;
    return;
  }

  TFile* fAANew = TFile::Open(kInAuAuGoldNew.c_str(), "READ");
  if (!fAANew || fAANew->IsZombie())
  {
    cout << ANSI_BOLD_YEL
         << "[WARN] Cannot open AuAu gold NEW file (UE-sub overlay disabled): " << kInAuAuGoldNew
         << ANSI_RESET << "\n";
    if (fAANew) { fAANew->Close(); delete fAANew; }
    fAANew = nullptr;
  }

  TDirectory* aaTopNew = nullptr;
  if (fAANew)
  {
    aaTopNew = fAANew->GetDirectory(kTriggerAuAuGold.c_str());
    if (!aaTopNew)
    {
      cout << ANSI_BOLD_YEL
           << "[WARN] Missing AuAu trigger directory '" << kTriggerAuAuGold
           << "' in NEW file (UE-sub overlay disabled): " << kInAuAuGoldNew
           << ANSI_RESET << "\n";
    }
}

if (!dsPP.topDir)
{
  cout << ANSI_BOLD_RED
       << "[ERROR] PP dataset topDir is null (cannot read PP histograms)."
       << ANSI_RESET << "\n";
  fAA->Close(); delete fAA;
  return;
}

// Output base (dedicated PP vs AuAu folder)
const string outBase = kOutPPAuAuBase;
EnsureDir(outBase);

const auto& centBins = CentBins();
if (centBins.empty())
{
  cout << ANSI_BOLD_YEL
       << "[WARN] centrality_edges missing/invalid (no centrality bins). Nothing to do."
       << ANSI_RESET << "\n";
  fAA->Close(); delete fAA;
  return;
}

// ---------------------------------------------------------------------
// Deliverable Set A: SS spectra (inclusive / iso / nonIso)
// ---------------------------------------------------------------------
const vector<string> ssVars = {"weta","wphi","et1","e11e33","e32e35"};

struct SSDef { string folder; string tag; string label; };
const vector<SSDef> ssDefs = {
  {"noIsoRequired",  "inclusive", "Inclusive"},
  {"isoPassSSplots", "iso",       "Iso pass"},
  {"isoFailSSplots", "nonIso",    "Iso fail"}
};

for (const auto& cb : centBins)
{
    const string centFolder = cb.folder;
    const string centSuffix = cb.suffix;
    const string centLabel  = TString::Format("cent: %d-%d%%", cb.lo, cb.hi).Data();

    cout << ANSI_BOLD_CYN
         << "\n[PP_AuAu] =========================================================\n"
         << "[PP_AuAu] Centrality bin: " << centLabel << "  (suffix=" << centSuffix << ", folder=" << centFolder << ")\n"
         << "[PP_AuAu] ========================================================="
         << ANSI_RESET << "\n";

    cout << ANSI_BOLD_CYN
         << "\n[PP_AuAu] --- Shower-shape tables (2x3) + per-pT overlays: Inclusive / Iso pass / Iso fail ---"
         << ANSI_RESET << "\n";

    for (const auto& def : ssDefs)
    {
      cout << ANSI_BOLD_GRN
           << "\n[PP_AuAu]   SS category: " << def.label
           << "  (folder=" << def.folder << ", tag=" << def.tag << ")"
           << ANSI_RESET << "\n";

      for (const auto& var : ssVars)
      {
        cout << "  [PP_AuAu]     SS var: " << var
             << "  -> building table2x3 + per-bin overlays" << "\n";

        const string outDir = JoinPath(outBase, JoinPath(def.folder, JoinPath(centFolder, var)));
        const string histBase = "h_ss_" + var + "_" + def.tag;

        const string topLeft = TString::Format("SS %s (%s)", var.c_str(), def.label.c_str()).Data();

        ProduceFamily_PPvsAuAu(dsPP.topDir, aaTop, outDir, histBase, centSuffix, centLabel,
                                 var, topLeft, false);

        // ------------------------------------------------------------
        // NEW: For noIsoRequired/e32e35 only, write an additional PP-only
        // unnormalized ZOOMED 2x3 table under:
        //   <outDir>/pp_unNormalized/table2x3_pp_unNormalized_zoom.png
        // Does NOT affect existing outputs.
        // ------------------------------------------------------------
        if (def.folder == "noIsoRequired" && var == "e32e35")
        {
            MakePPUnNormalizedZoomTable_e32e35(dsPP.topDir, outDir, centLabel);
        }

        // ---------------------------------------------------------------------
        // NEW: In noIsoRequired only, also write FULL-range UNNORMALIZED SS
        // distributions (PP + AuAu) per pT bin per centrality, plus 2x3 tables.
        // Output inside:
        //   .../noIsoRequired/<cent>/<ssVar>/{AuAu_unNormalized,pp_unNormalized}/...
        // ---------------------------------------------------------------------
        if (def.folder == "noIsoRequired")
        {
            MakeUnnormalizedQA_PPvsAuAu(dsPP.topDir, aaTop, outDir, histBase, centSuffix, centLabel, var, topLeft);
        }
      }
    }

    // ---------------------------------------------------------------------
    // Deliverable Set B: total Eiso spectra split by tightness
    // ---------------------------------------------------------------------
    cout << ANSI_BOLD_CYN
         << "\n[PP_AuAu] --- Isolation spectra (2x3) + per-pT overlays: inclusive / tight / nonTight ---"
         << ANSI_RESET << "\n";

    struct IsoDef { string folder; string base; string label; };
    const vector<IsoDef> isoDefs = {
      {"noSS_isoSpectra",    "h_Eiso",          "E_{T}^{iso, Total} (inclusive)"},
      {"tightIsoSpectra",    "h_Eiso_tight",    "E_{T}^{iso, Total} (tight pass)"},
      {"nonTightIsoSpectra", "h_Eiso_nonTight", "E_{T}^{iso, Total} (tight fail)"}
    };

    for (const auto& idef : isoDefs)
    {
      cout << ANSI_BOLD_GRN
           << "\n[PP_AuAu]   Iso category: " << idef.label
           << "  (folder=" << idef.folder << ", base=" << idef.base << ")"
           << ANSI_RESET << "\n";

      const string outDir = JoinPath(outBase, JoinPath(idef.folder, centFolder));
      ProduceFamily_PPvsAuAu(dsPP.topDir, aaTop, outDir, idef.base, centSuffix, centLabel,
                               "E_{T}^{iso} [GeV]", idef.label, true);

      // ---------------------------------------------------------------------
      // NEW: For inclusive noSS_isoSpectra only, also write an AuAu (counts)
      // 2x3 table overlaying the current AuAu file vs the UE-subtracted file.
      // Output:
      //   <outDir>/AuAu_unNormalized/table2x3_AuAu_unNormalized_overlay_UEsub.png
      // ---------------------------------------------------------------------
      if (idef.folder == "noSS_isoSpectra" && aaTopNew)
      {
          Make2x3Table_AuAuUnNormalized_OverlayUEsub(aaTop, aaTopNew, outDir,
                                                    idef.base, centSuffix, centLabel,
                                                    "E_{T}^{iso} [GeV]");
        }
    }

    // ---------------------------------------------------------------------
    // NEW: Tight vs nonTight overlays (PP and AuAu per centrality)
    // Output:
    //   kOutPPAuAuBase/tight_nonTight_isoSpectraOverlay/pp/...
    //   kOutPPAuAuBase/tight_nonTight_isoSpectraOverlay/AuAu/<cent>/...
    // ---------------------------------------------------------------------
    {
      const string tnBase = JoinPath(outBase, "tight_nonTight_isoSpectraOverlay");
      EnsureDir(tnBase);
      ProduceTightNonTightIsoOverlays(dsPP.topDir, aaTop, tnBase, centFolder, centSuffix, centLabel);
    }

    // ---------------------------------------------------------------------
    // NEW: isoFail vs isoPass SS overlays (PP and AuAu per centrality)
    // Output:
    //   kOutPPAuAuBase/isoFail_isoPass_SSspectra/pp/<ssVar>/...
    //   kOutPPAuAuBase/isoFail_isoPass_SSspectra/AuAu/<cent>/<ssVar>/...
    // ---------------------------------------------------------------------
    {
      const string outSS = JoinPath(outBase, "isoFail_isoPass_SSspectra");
      EnsureDir(outSS);
      ProduceIsoFailIsoPassSSOverlays(dsPP.topDir, aaTop, outSS, centFolder, centSuffix, centLabel, ssVars);
    }
  }

  // ---------------------------------------------------------------------
  // NEW (noIsoRequired summary outside per-centrality folders):
  //
  // Output to:
  //   <kOutPPAuAuBase>/noIsoRequired/<ssVar>/<pTbin>/table2x3_overlay_pp_vs_auau_byCent.png
  //
  // For each SS variable and each pT bin, produce a 2x3 table where each pad
  // corresponds to a centrality bin (first 6), overlaying:
  //   - PP (same curve in every pad for that pT bin)
  //   - AuAu (centrality-dependent curve)
  // with legend + centered title including SS var, centrality, and pT bin.
  // ---------------------------------------------------------------------
  {
    const int nPtPads = std::min(6, kNPtBins);
    const int nCentPads = std::min(6, (int)centBins.size());

    if (nPtPads > 0 && nCentPads > 0)
    {
      const auto& ptBins = PtBins();

      const string baseNoIso = JoinPath(outBase, "noIsoRequired");
      EnsureDir(baseNoIso);

      for (const auto& var : ssVars)
      {
            const string varBase = JoinPath(baseNoIso, var);
            EnsureDir(varBase);

            // ---------------------------------------------------------------------
            // NEW: In noIsoRequired/<ssVar>/, also write TWO AuAu-only summary tables:
            //   (A) 2x3 table of first 6 pT bins, each pad overlays centralities (AuAu only)
            //       -> noIsoRequired/<ssVar>/table2x3_AuAu_unNormalized.png
            //   (B) 2x3 table of first 6 centralities, each pad overlays ALL pT bins (AuAu only)
            //       -> noIsoRequired/<ssVar>/table2x3_AuAu_unNormalized_byPtOverlays.png
            // Both include centered title + SS cut label + vertical cut lines.
            // ---------------------------------------------------------------------
            {
              const string histBaseAuAu = string("h_ss_") + var + string("_inclusive");
              Make2x3Table_AuAuUnNormalized_SS_ByCent(aaTop, varBase, histBaseAuAu, var);
              Make2x3Table_AuAuUnNormalized_SS_ByPtOverlaysPerCent(aaTop, varBase, histBaseAuAu, var);
            }

            // ---------------------------------------------------------------------
            // NEW: In noIsoRequired/<ssVar>/ (OUTSIDE pT folders), write:
            //   (1) PP inclusive overlay across ALL pT bins (from YAML PtBins)
            //       -> noIsoRequired/<ssVar>/overlay_pp_byPt.png
            //   (2) AuAu inclusive overlay across ALL pT bins, PER centrality bin
            //       -> noIsoRequired/<ssVar>/overlay_auau_<centFolder>_byPt.png
            // All curves: closed circle markers, distinct colors, legend = pT bins.
            // ---------------------------------------------------------------------
          {
            const int nAllPt = kNPtBins;
            if (nAllPt > 0)
            {
              const int colors[] = {
                kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 1,
                kViolet + 1, kAzure + 2, kSpring + 5, kPink + 7, kTeal + 3, kGray + 2
              };
              const int nColors = (int)(sizeof(colors)/sizeof(colors[0]));

              // -------------------------
              // (1) PP: overlay ALL pT bins
              // -------------------------
              {
                TCanvas cPPall(
                  TString::Format("c_noIso_pp_byPt_%s", var.c_str()).Data(),
                  "c_noIso_pp_byPt", 900, 700
                );
                ApplyCanvasMargins1D(cPPall);
                cPPall.SetLogy(false);

                double yMax = 0.0;
                TH1* hFirst = nullptr;

                  TLegend leg(0.47, 0.65, 0.92, 0.91);
                  leg.SetBorderSize(0);
                  leg.SetFillStyle(0);
                  leg.SetTextFont(42);
                  leg.SetTextSize(0.030);
                  leg.SetNColumns(2);

                  // For e32e35 only: place the pT-bin legend in the middle-left of the canvas
                  // (keeps the rest of the overlay_pp_byPt plots unchanged).
                  if (var == "e32e35")
                  {
                    leg.SetX1NDC(0.14);
                    leg.SetX2NDC(0.44);
                    leg.SetY1NDC(0.35);
                    leg.SetY2NDC(0.72);
                    leg.SetTextSize(0.032);
                    leg.SetNColumns(1);
                  }
                  // For et1 only: place the pT-bin legend in the middle-left of the canvas
                  if (var == "et1")
                  {
                    leg.SetX1NDC(0.14);
                    leg.SetX2NDC(0.44);
                    leg.SetY1NDC(0.35);
                    leg.SetY2NDC(0.72);
                    leg.SetTextSize(0.032);
                    leg.SetNColumns(1);
                  }

                vector<TH1*> keepAlivePP;
                keepAlivePP.reserve((std::size_t)nAllPt);

                for (int iptAll = 0; iptAll < nAllPt; ++iptAll)
                {
                  const PtBin& pbAll = ptBins[iptAll];
                  const string hPPName = string("h_ss_") + var + string("_inclusive") + pbAll.suffix;
                  TH1* rawPP = GetTH1FromTopDir(dsPP.topDir, hPPName);
                  if (!rawPP) continue;

                  const int col = colors[iptAll % nColors];

                  TH1* hPPc = CloneNormalizeStyle(rawPP,
                    TString::Format("pp_noIso_byPt_%s_%s", var.c_str(), pbAll.folder.c_str()).Data(),
                    col, 20);

                  if (!hPPc) continue;

                  hPPc->SetMarkerStyle(20);
                  hPPc->SetMarkerSize(1.05);

                  yMax = std::max(yMax, hPPc->GetMaximum());

                  if (!hFirst)
                  {
                    hFirst = hPPc;
                    hFirst->GetXaxis()->SetTitle(var.c_str());
                    hFirst->GetYaxis()->SetTitle("Normalized counts");
                    hFirst->Draw("E1");
                  }
                  else
                  {
                    hPPc->Draw("E1 same");
                  }

                  leg.AddEntry(hPPc,
                    TString::Format("p_{T}^{#gamma}: %d-%d GeV", pbAll.lo, pbAll.hi).Data(),
                    "ep");

                  keepAlivePP.push_back(hPPc);
                }

                  if (hFirst)
                  {
                    hFirst->SetMaximum(yMax * 1.35);

                    TLatex t;
                    t.SetNDC(true);
                    t.SetTextFont(42);
                    t.SetTextAlign(22);
                    t.SetTextSize(0.045);
                    t.DrawLatex(0.50, 0.93,
                      TString::Format("pp (no Iso Req), %s, overlay across p_{T}^{#gamma} bins", var.c_str()).Data());

                    // ------------------------------------------------------------
                    // Add SS cut label (top-left) + cut lines (where applicable)
                    // ------------------------------------------------------------
                    TLatex tcut;
                    tcut.SetNDC(true);
                    tcut.SetTextFont(42);
                    tcut.SetTextAlign(13);
                    tcut.SetTextSize(0.038);

                    bool drawCuts = false;
                    double cutLo = 0.0;
                    double cutHi = 0.0;
                    std::string cutText;

                    if (var == "e11e33")
                    {
                      cutText = "Tight #gamma-ID: 0.4 < #frac{E_{11}}{E_{33}} < 0.98";
                      drawCuts = true;
                      cutLo = 0.4;
                      cutHi = 0.98;
                    }
                    else if (var == "e32e35")
                    {
                      cutText = "#gamma-ID: 0.92 < #frac{E_{32}}{E_{35}} < 1.0";
                      drawCuts = true;
                      cutLo = 0.92;
                      cutHi = 1.0;
                    }
                    else if (var == "et1")
                    {
                      cutText = "#gamma-ID: 0.9 < et1 < 1.0";
                      drawCuts = true;
                      cutLo = 0.9;
                      cutHi = 1.0;
                    }
                    else if (var == "weta")
                    {
                      cutText = "#gamma-ID: 0 < w_{#eta}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
                    }
                    else if (var == "wphi")
                    {
                      cutText = "#gamma-ID: 0 < w_{#phi}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
                    }

                    if (!cutText.empty())
                    {
                      tcut.DrawLatex(0.16, 0.89, cutText.c_str());
                    }

                    leg.Draw();

                    if (drawCuts)
                    {
                      gPad->Update();
                      const double yMin = gPad->GetUymin();
                      const double yMaxPad = gPad->GetUymax();

                      TLine* l1 = new TLine(cutLo, yMin, cutLo, yMaxPad);
                      l1->SetLineColor(kGreen + 2);
                      l1->SetLineWidth(2);
                      l1->SetLineStyle(2);
                      l1->Draw("same");

                      TLine* l2 = new TLine(cutHi, yMin, cutHi, yMaxPad);
                      l2->SetLineColor(kOrange + 7);
                      l2->SetLineWidth(2);
                      l2->SetLineStyle(2);
                      l2->Draw("same");
                    }

                    gPad->RedrawAxis();

                    SaveCanvas(cPPall, JoinPath(varBase, "overlay_pp_byPt.png"));
                  }

                for (TH1* h : keepAlivePP) delete h;
                keepAlivePP.clear();
              }

              // -------------------------
              // (2) AuAu: per-centrality overlay ALL pT bins
              // -------------------------
              for (size_t ic = 0; ic < centBins.size(); ++ic)
              {
                const auto& cb2 = centBins[ic];
                const string centSuffix2 = cb2.suffix;
                const string centFolder2 = cb2.folder;
                const string centLabel2  = TString::Format("Cent = %d-%d%%", cb2.lo, cb2.hi).Data();

                TCanvas cAAall(
                  TString::Format("c_noIso_aa_byPt_%s_%s", var.c_str(), centFolder2.c_str()).Data(),
                  "c_noIso_aa_byPt", 900, 700
                );
                ApplyCanvasMargins1D(cAAall);
                cAAall.SetLogy(false);

                double yMax = 0.0;
                TH1* hFirst = nullptr;

                TLegend leg(0.52, 0.58, 0.90, 0.86);
                leg.SetBorderSize(0);
                leg.SetFillStyle(0);
                leg.SetTextFont(42);
                leg.SetTextSize(0.034);

                vector<TH1*> keepAliveAA;
                keepAliveAA.reserve((std::size_t)nAllPt);

                for (int iptAll = 0; iptAll < nAllPt; ++iptAll)
                {
                  const PtBin& pbAll = ptBins[iptAll];
                  const string hAAName = string("h_ss_") + var + string("_inclusive") + pbAll.suffix + centSuffix2;
                  TH1* rawAA = GetTH1FromTopDir(aaTop, hAAName);
                  if (!rawAA) continue;

                  const int col = colors[iptAll % nColors];

                  TH1* hAAc = CloneNormalizeStyle(rawAA,
                    TString::Format("aa_noIso_byPt_%s_%s%s", var.c_str(), pbAll.folder.c_str(), centSuffix2.c_str()).Data(),
                    col, 20);

                  if (!hAAc) continue;

                  hAAc->SetMarkerStyle(20);
                  hAAc->SetMarkerSize(1.05);

                  yMax = std::max(yMax, hAAc->GetMaximum());

                  if (!hFirst)
                  {
                    hFirst = hAAc;
                    hFirst->GetXaxis()->SetTitle(var.c_str());
                    hFirst->GetYaxis()->SetTitle("Normalized counts");
                    hFirst->Draw("E1");
                  }
                  else
                  {
                    hAAc->Draw("E1 same");
                  }

                  leg.AddEntry(hAAc,
                    TString::Format("p_{T}^{#gamma}: %d-%d GeV", pbAll.lo, pbAll.hi).Data(),
                    "ep");

                  keepAliveAA.push_back(hAAc);
                }

                if (hFirst)
                {
                  hFirst->SetMaximum(yMax * 1.35);

                  TLatex t;
                  t.SetNDC(true);
                  t.SetTextFont(42);
                  t.SetTextAlign(22);
                  t.SetTextSize(0.045);
                  t.DrawLatex(0.50, 0.93,
                    TString::Format("Au+Au (inclusive), %s, %s, overlay across p_{T}^{#gamma} bins",
                                    var.c_str(), centLabel2.c_str()).Data());

                  leg.Draw();
                  gPad->RedrawAxis();

                  SaveCanvas(cAAall, JoinPath(varBase,
                    TString::Format("overlay_auau_%s_byPt.png", centFolder2.c_str()).Data()));
                }

                for (TH1* h : keepAliveAA) delete h;
                keepAliveAA.clear();
              }
            }
        }

        for (int ipt = 0; ipt < nPtPads; ++ipt)
        {
          const PtBin& pb = ptBins[ipt];

          const string ptDir = JoinPath(varBase, pb.folder);
          EnsureDir(ptDir);

          TCanvas c(
            TString::Format("c_noIso_centSummary_%s_%s", var.c_str(), pb.folder.c_str()).Data(),
            "c_noIso_centSummary", 1500, 800
          );
          c.Divide(3,2, 0.001, 0.001);

          vector<TH1*> keepAlive;
          keepAlive.reserve((std::size_t)nCentPads * 2);

          // PP inclusive SS hist name (same for all centrality pads)
          const string hPPName = string("h_ss_") + var + string("_inclusive") + pb.suffix;
          TH1* rawPP = GetTH1FromTopDir(dsPP.topDir, hPPName);

          for (int ic = 0; ic < nCentPads; ++ic)
          {
            const auto& cb2 = centBins[ic];
            const string centSuffix2 = cb2.suffix;
            const string centLabel2  = TString::Format("Cent = %d-%d%%", cb2.lo, cb2.hi).Data();

            c.cd(ic+1);
            gPad->SetLeftMargin(0.14);
            gPad->SetRightMargin(0.05);
            gPad->SetBottomMargin(0.14);
            gPad->SetTopMargin(0.12);
            gPad->SetLogy(false);

            // AuAu inclusive SS hist name for this centrality
            const string hAAName = string("h_ss_") + var + string("_inclusive") + pb.suffix + centSuffix2;
            TH1* rawAA = GetTH1FromTopDir(aaTop, hAAName);

            if (!rawPP || !rawAA)
            {
              std::ostringstream s;
              s << "pT: " << pb.lo << "-" << pb.hi << "  " << centLabel2;
              DrawMissingPad(s.str());
              continue;
            }

            TH1* hPP = CloneNormalizeStyle(rawPP,
              TString::Format("pp_noIso_centSummary_%s_%s_%d", var.c_str(), pb.folder.c_str(), ic).Data(),
              kBlack, 20);

            TH1* hAA = CloneNormalizeStyle(rawAA,
              TString::Format("aa_noIso_centSummary_%s_%s%s_%d", var.c_str(), pb.folder.c_str(), centSuffix2.c_str(), ic).Data(),
              kRed + 1, 24);

            if (!hPP || !hAA)
            {
              if (hPP) delete hPP;
              if (hAA) delete hAA;
              std::ostringstream s;
              s << "pT: " << pb.lo << "-" << pb.hi << "  " << centLabel2;
              DrawMissingPad(s.str());
              continue;
            }

            hPP->GetXaxis()->SetTitle(var.c_str());
            hPP->GetYaxis()->SetTitle("Normalized counts");

            const double ymax = std::max(hPP->GetMaximum(), hAA->GetMaximum());
            hPP->SetMaximum(ymax * 1.35);

            hPP->Draw("E1");
            hAA->Draw("E1 same");

            // Always label BOTH curves explicitly (PP vs Au+Au).
            TLegend leg(0.52, 0.68, 0.90, 0.84);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            leg.SetTextFont(42);
            leg.SetTextSize(0.034);
            leg.AddEntry(hPP, "PP data", "ep");
            leg.AddEntry(hAA, "Au+Au data", "ep");
            leg.Draw();

              // Centered title for each pad
              TLatex t;
              t.SetNDC(true);
              t.SetTextFont(42);
              t.SetTextAlign(22);
              t.SetTextSize(0.040);

              t.DrawLatex(0.50, 0.93,
                TString::Format("%s, %s, p_{T}^{#gamma} = %d-%d GeV",
                                var.c_str(), centLabel2.c_str(), pb.lo, pb.hi).Data());

              // ------------------------------------------------------------
              // Add #gamma-ID annotation (top-left) + draw cut lines where applicable
              // ------------------------------------------------------------
              TLatex tcut;
              tcut.SetNDC(true);
              tcut.SetTextFont(42);
              tcut.SetTextAlign(13);
              tcut.SetTextSize(0.036);

              bool drawCuts = false;
              double cutLo = 0.0;
              double cutHi = 0.0;
              std::string cutText;

              if (var == "e11e33")
              {
                cutText = "Tight #gamma-ID: 0.4 < #frac{E_{11}}{E_{33}} < 0.98";
                drawCuts = true;
                cutLo = 0.4;
                cutHi = 0.98;
              }
              else if (var == "e32e35")
              {
                cutText = "#gamma-ID: 0.92 < #frac{E_{32}}{E_{35}} < 1.0";
                drawCuts = true;
                cutLo = 0.92;
                cutHi = 1.0;
              }
              else if (var == "et1")
              {
                cutText = "#gamma-ID: 0.9 < et1 < 1.0";
                drawCuts = true;
                cutLo = 0.9;
                cutHi = 1.0;
              }
              else if (var == "weta")
              {
                cutText = "#gamma-ID: 0 < w_{#eta}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
              }
              else if (var == "wphi")
              {
                cutText = "#gamma-ID: 0 < w_{#phi}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
              }

              if (!cutText.empty())
              {
                tcut.DrawLatex(0.16, 0.84, cutText.c_str());
              }

              if (drawCuts)
              {
                // Use pad y-range (not hist max) and force "same" draw so lines persist/appear.
                gPad->Update();
                const double yMin = gPad->GetUymin();
                const double yMax = gPad->GetUymax();

                TLine* l1 = new TLine(cutLo, yMin, cutLo, yMax);
                l1->SetLineColor(kGreen + 2);
                l1->SetLineWidth(2);
                l1->SetLineStyle(2);
                l1->Draw("same");

                TLine* l2 = new TLine(cutHi, yMin, cutHi, yMax);
                l2->SetLineColor(kOrange + 7);
                l2->SetLineWidth(2);
                l2->SetLineStyle(2);
                l2->Draw("same");

                gPad->RedrawAxis();
              }

            keepAlive.push_back(hPP);
            keepAlive.push_back(hAA);
          }

            SaveCanvas(c, JoinPath(ptDir, "table2x3_overlay_pp_vs_auau_byCent.png"));

            // ------------------------------------------------------------------
            // NEW: AuAu-only overlay by centrality (all closed circles)
            // Output to:
            //   <kOutPPAuAuBase>/noIsoRequired/<ssVar>/<pTbin>/overlay_auau_byCent.png
            // This is produced for e11e33 and ALL shower-shape variables (ssVars).
            // ------------------------------------------------------------------
            {
              TCanvas cCent(
                TString::Format("c_auau_byCent_%s_%s", var.c_str(), pb.folder.c_str()).Data(),
                "c_auau_byCent", 900, 700
              );
              ApplyCanvasMargins1D(cCent);
              cCent.SetLogy(false);

              const int colors[] = {
                kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 1,
                kViolet + 1, kAzure + 2, kSpring + 5, kPink + 7, kTeal + 3, kGray + 2
              };
              const int nColors = (int)(sizeof(colors)/sizeof(colors[0]));

              double yMax = 0.0;
              TH1* hFirst = nullptr;

              TLegend leg(0.52, 0.58, 0.90, 0.86);
              leg.SetBorderSize(0);
              leg.SetFillStyle(0);
              leg.SetTextFont(42);
              leg.SetTextSize(0.034);

              for (size_t ic = 0; ic < centBins.size(); ++ic)
              {
                const auto& cb2 = centBins[ic];
                const string centSuffix2 = cb2.suffix;
                const string centLabel2  = TString::Format("%d-%d%%", cb2.lo, cb2.hi).Data();

                const string hAAName = string("h_ss_") + var + string("_inclusive") + pb.suffix + centSuffix2;
                TH1* rawAA = GetTH1FromTopDir(aaTop, hAAName);
                if (!rawAA) continue;

                const int col = colors[(int)(ic % (size_t)nColors)];

                TH1* hAAc = CloneNormalizeStyle(rawAA,
                  TString::Format("aa_byCent_%s_%s%s_%zu", var.c_str(), pb.folder.c_str(), centSuffix2.c_str(), ic).Data(),
                  col, 20);

                if (!hAAc) continue;

                // Force closed circles for every centrality bin
                hAAc->SetMarkerStyle(20);
                hAAc->SetMarkerSize(1.05);

                yMax = std::max(yMax, hAAc->GetMaximum());

                if (!hFirst)
                {
                  hFirst = hAAc;
                  hFirst->GetXaxis()->SetTitle(var.c_str());
                  hFirst->GetYaxis()->SetTitle("Normalized counts");
                  hFirst->Draw("E1");
                }
                else
                {
                  hAAc->Draw("E1 same");
                }

                leg.AddEntry(hAAc, centLabel2.c_str(), "ep");
                keepAlive.push_back(hAAc);
              }

              if (hFirst)
              {
                hFirst->SetMaximum(yMax * 1.35);

                TLatex t;
                t.SetNDC(true);
                t.SetTextFont(42);
                t.SetTextAlign(22);
                t.SetTextSize(0.045);
                t.DrawLatex(0.50, 0.93,
                  TString::Format("Au+Au by centrality, %s, p_{T}^{#gamma} = %d-%d GeV",
                                  var.c_str(), pb.lo, pb.hi).Data());

                leg.Draw();
                gPad->RedrawAxis();

                SaveCanvas(cCent, JoinPath(ptDir, "overlay_auau_byCent.png"));
              }
            }

            // ------------------------------------------------------------------
            // NEW: PP-only zoomed distribution for this SS variable + pT bin
            // Output to:
            //   <kOutPPAuAuBase>/noIsoRequired/<ssVar>/<pTbin>/pp_zoom.png
            // ------------------------------------------------------------------
            if (rawPP)
            {
              TH1* hPPz = CloneNormalizeStyle(rawPP,
                TString::Format("pp_zoom_%s_%s", var.c_str(), pb.folder.c_str()).Data(),
                kBlack, 20);

              if (hPPz)
              {
                TCanvas cPPz(
                  TString::Format("c_pp_zoom_%s_%s", var.c_str(), pb.folder.c_str()).Data(),
                  "c_pp_zoom", 900, 700
                );
                ApplyCanvasMargins1D(cPPz);
                cPPz.SetLogy(false);

                hPPz->GetXaxis()->SetTitle(var.c_str());
                hPPz->GetYaxis()->SetTitle("Normalized counts");

                // Default zoom windows (reasonable + stable)
                // - For ratio-like vars: focus around the tight-id region
                // - For widths: show low region where cut lives
                if (var == "e11e33")
                {
                  hPPz->GetXaxis()->SetRangeUser(0.25, 1.05);
                }
                else if (var == "e32e35")
                {
                  hPPz->GetXaxis()->SetRangeUser(0.85, 1.05);
                }
                else if (var == "et1")
                {
                  hPPz->GetXaxis()->SetRangeUser(0.80, 1.05);
                }
                else if (var == "weta" || var == "wphi")
                {
                  hPPz->GetXaxis()->SetRangeUser(0.0, 0.40);
                }

                hPPz->Draw("E1");

                // Add the same gamma-ID annotation + lines (where applicable)
                TLatex tcut;
                tcut.SetNDC(true);
                tcut.SetTextFont(42);
                tcut.SetTextAlign(13);
                tcut.SetTextSize(0.040);

                bool drawCuts = false;
                double cutLo = 0.0;
                double cutHi = 0.0;
                std::string cutText;

                if (var == "e11e33")
                {
                  cutText = "Tight #gamma-ID: 0.4 < #frac{E_{11}}{E_{33}} < 0.98";
                  drawCuts = true;
                  cutLo = 0.4;
                  cutHi = 0.98;
                }
                else if (var == "e32e35")
                {
                  cutText = "#gamma-ID: 0.92 < #frac{E_{32}}{E_{35}} < 1.0";
                  drawCuts = true;
                  cutLo = 0.92;
                  cutHi = 1.0;
                }
                else if (var == "et1")
                {
                  cutText = "#gamma-ID: 0.9 < et1 < 1.0";
                  drawCuts = true;
                  cutLo = 0.9;
                  cutHi = 1.0;
                }
                else if (var == "weta")
                {
                  cutText = "#gamma-ID: 0 < w_{#eta}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
                }
                else if (var == "wphi")
                {
                  cutText = "#gamma-ID: 0 < w_{#phi}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
                }

                if (!cutText.empty())
                {
                  tcut.DrawLatex(0.16, 0.86, cutText.c_str());
                }

                if (drawCuts)
                {
                  gPad->Update();
                  const double yMin = gPad->GetUymin();
                  const double yMax = gPad->GetUymax();

                  TLine* l1 = new TLine(cutLo, yMin, cutLo, yMax);
                  l1->SetLineColor(kGreen + 2);
                  l1->SetLineWidth(2);
                  l1->SetLineStyle(2);
                  l1->Draw("same");

                  TLine* l2 = new TLine(cutHi, yMin, cutHi, yMax);
                  l2->SetLineColor(kOrange + 7);
                  l2->SetLineWidth(2);
                  l2->SetLineStyle(2);
                  l2->Draw("same");

                  gPad->RedrawAxis();
                }

                // Centered title
                TLatex t;
                t.SetNDC(true);
                t.SetTextFont(42);
                t.SetTextAlign(22);
                t.SetTextSize(0.045);
                t.DrawLatex(0.50, 0.94,
                  TString::Format("%s, PP, p_{T}^{#gamma} = %d-%d GeV",
                                  var.c_str(), pb.lo, pb.hi).Data());

                SaveCanvas(cPPz, JoinPath(ptDir, "pp_zoom.png"));

                delete hPPz;
              }
            }

            for (TH1* h : keepAlive) delete h;
            keepAlive.clear();
        }
      }
    }
  }

  // ---------------------------------------------------------------------
  // NEW: Per-pT unnormalized PP/AuAu overlays across centrality (skip 80-100),
  // rendered as a 2x3 canvas containing ALL SS variables.
  //
  // Layout:
  //   [ weta ] [ wphi ] [ blank ]
  //   [ e11e33 ] [ et1 ] [ e32e35 ]
  //
  // Output (saved into each SS variable's pT folder):
  //   <outBase>/noIsoRequired/<ssVar>/<pTbin>/table2x3_ppAuAu_unNormalized_SS_overlaysByCent.png
  // ---------------------------------------------------------------------
  {
      if (aaTop)
      {
        const string baseNoIso = JoinPath(outBase, "noIsoRequired");
        EnsureDir(baseNoIso);

        const auto& ptBinsLocal   = PtBins();
        const auto& centBinsLocal = CentBins();

        if (!ptBinsLocal.empty() && !centBinsLocal.empty())
        {
          struct VarDef { std::string var; std::string label; };
          const std::vector<VarDef> vars =
          {
            {"weta",   "w_{#eta}"},
            {"wphi",   "w_{#phi}"},
            {"e11e33", "E_{11}/E_{33}"},
            {"et1",    "et1"},
            {"e32e35", "E_{32}/E_{35}"}
          };

          auto LabelForVar = [&](const std::string& v) -> std::string
          {
            for (const auto& vd : vars) { if (vd.var == v) return vd.label; }
            return v;
          };

          const int centColors[] = { kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 1, kViolet + 1, kAzure + 2 };
          const int nCentColors = (int)(sizeof(centColors)/sizeof(centColors[0]));

          for (int ipt = 0; ipt < (int)ptBinsLocal.size(); ++ipt)
          {
            const PtBin& pb = ptBinsLocal[ipt];

            TCanvas cAll(
              TString::Format("c_noIso_ppAuAu_unNorm_byCent_allSS_%s", pb.folder.c_str()).Data(),
              "c_noIso_ppAuAu_unNorm_byCent_allSS", 2100, 1400
            );
            cAll.Divide(3, 2, 0.001, 0.001);

            vector<TH1*> keepH;
            keepH.reserve((std::size_t)vars.size() * (std::size_t)centBinsLocal.size());

            vector<TLegend*> keepL;
            keepL.reserve((std::size_t)vars.size());

            vector<TLine*> keepLines;
            keepLines.reserve((std::size_t)vars.size() * 2);

            const std::string padVars[6] = {"", "weta", "wphi", "e11e33", "et1", "e32e35"};

            for (int ipad = 0; ipad < 6; ++ipad)
            {
              const std::string v = padVars[ipad];
              cAll.cd(ipad + 1);

              gPad->SetLeftMargin(0.14);
              gPad->SetRightMargin(0.05);
              gPad->SetBottomMargin(0.14);
              gPad->SetTopMargin(0.18);

              if (v.empty())
              {
                gPad->SetLogy(false);
                gPad->Clear();
                continue;
              }

              gPad->SetLogy(true);

              const std::string vlabel = LabelForVar(v);
              const string histBase = string("h_ss_") + v + string("_inclusive");

              TH1* rawPP = GetTH1FromTopDir(dsPP.topDir, histBase + pb.suffix);

              TH1* hPP = nullptr;
              if (rawPP)
              {
                hPP = CloneTH1(rawPP,
                  TString::Format("pp_unNorm_byCent_%s_%s", v.c_str(), pb.folder.c_str()).Data());
                if (hPP)
                {
                  EnsureSumw2(hPP);
                  hPP->GetXaxis()->UnZoom();
                  hPP->SetTitle("");
                  hPP->GetXaxis()->SetTitle(vlabel.c_str());
                  hPP->GetYaxis()->SetTitle("Counts");
                  StyleOverlayHist(hPP, kBlack, 20);
                  hPP->SetMarkerStyle(20);
                  hPP->SetMarkerSize(1.05);
                  keepH.push_back(hPP);
                }
              }

              std::vector<TH1*> hAAs;
              std::vector<std::string> aaLabels;

              for (int ic = 0, icDraw = 0; ic < (int)centBinsLocal.size(); ++ic)
              {
                const auto& cb = centBinsLocal[(std::size_t)ic];
                if (cb.lo == 80 && cb.hi == 100) continue;

                TH1* rawAA = GetTH1FromTopDir(aaTop, histBase + pb.suffix + cb.suffix);
                if (!rawAA) continue;

                const int col = centColors[icDraw % nCentColors];
                ++icDraw;

                TH1* hAA = CloneTH1(rawAA,
                  TString::Format("aa_unNorm_byCent_%s_%s%s", v.c_str(), pb.folder.c_str(), cb.suffix.c_str()).Data());
                if (!hAA) continue;

                EnsureSumw2(hAA);
                hAA->GetXaxis()->UnZoom();
                hAA->SetTitle("");
                hAA->GetXaxis()->SetTitle(vlabel.c_str());
                hAA->GetYaxis()->SetTitle("Counts");
                StyleOverlayHist(hAA, col, 24);
                hAA->SetMarkerStyle(24);
                hAA->SetMarkerSize(1.05);

                hAAs.push_back(hAA);
                aaLabels.push_back(TString::Format("AuAu %d-%d%%", cb.lo, cb.hi).Data());
                keepH.push_back(hAA);
              }

              if (!hPP && hAAs.empty())
              {
                DrawMissingPad(TString::Format("%s, p_{T}^{#gamma} = %d-%d GeV", vlabel.c_str(), pb.lo, pb.hi).Data());
                continue;
              }

              double yMax = 0.0;
              double minPos = 1e99;

              auto ScanHist = [&](TH1* h)
              {
                if (!h) return;
                yMax = std::max(yMax, (double)h->GetMaximum());
                const int nb = h->GetNbinsX();
                for (int ib = 1; ib <= nb; ++ib)
                {
                  const double y = h->GetBinContent(ib);
                  if (y > 0.0 && y < minPos) minPos = y;
                }
              };

              ScanHist(hPP);
              for (TH1* h : hAAs) ScanHist(h);

              if (!(minPos < 1e98)) minPos = 0.5;

              const double yMin = minPos * 0.5;
              const double yMaxDraw = (yMax > 0.0 ? (yMax * 3.0) : 1.0);

              TH1* hFirst = hPP ? hPP : hAAs[0];
              hFirst->SetMinimum(yMin);
              hFirst->SetMaximum(yMaxDraw);
              hFirst->Draw("E1");

              if (hPP && hFirst != hPP) hPP->Draw("E1 same");
              for (TH1* h : hAAs) { if (h != hFirst) h->Draw("E1 same"); }

              const bool isTopRow = (ipad < 3);
              const bool shiftTopTwoLeft = (ipad == 1 || ipad == 2);

              const double lx1 = isTopRow ? (shiftTopTwoLeft ? 0.48 : 0.52) : 0.16;
              const double lx2 = isTopRow ? (shiftTopTwoLeft ? 0.89 : 0.93) : 0.57;
              const double ly1 = 0.52;
              const double ly2 = 0.78;

              TLegend* leg = new TLegend(lx1, ly1, lx2, ly2);
              leg->SetBorderSize(0);
              leg->SetFillStyle(0);
              leg->SetTextFont(42);
              leg->SetTextSize(0.030);
              leg->SetNColumns(2);
              leg->SetEntrySeparation(0.10);
              leg->SetColumnSeparation(0.12);

              if (hPP) leg->AddEntry(hPP, "PP (Run24pp)", "ep");
              for (std::size_t j = 0; j < hAAs.size(); ++j)
                  leg->AddEntry(hAAs[j], aaLabels[j].c_str(), "ep");

              leg->Draw();
              keepL.push_back(leg);

              {
                TLatex t;
                t.SetNDC(true);
                t.SetTextFont(42);
                t.SetTextAlign(22);
                t.SetTextSize(0.045);
                t.DrawLatex(0.50, 0.955,
                  TString::Format("%s, p_{T}^{#gamma} = %d-%d GeV, pp/AuAu overlays",
                                  vlabel.c_str(), pb.lo, pb.hi).Data());
              }

              {
                TLatex tcut;
                tcut.SetNDC(true);
                tcut.SetTextFont(42);
                tcut.SetTextAlign(13);
                tcut.SetTextSize(0.034);

                bool drawCuts = false;
                double cutLo = 0.0;
                double cutHi = 0.0;
                std::string cutText;

                if (v == "e11e33")
                {
                  cutText = "#gamma-ID: 0.4 < #frac{E_{11}}{E_{33}} < 0.98";
                  drawCuts = true;
                  cutLo = 0.4;
                  cutHi = 0.98;
                }
                else if (v == "e32e35")
                {
                  cutText = "#gamma-ID: 0.92 < #frac{E_{32}}{E_{35}} < 1.0";
                  drawCuts = true;
                  cutLo = 0.92;
                  cutHi = 1.0;
                }
                else if (v == "et1")
                {
                  cutText = "#gamma-ID: 0.9 < et1 < 1.0";
                  drawCuts = true;
                  cutLo = 0.9;
                  cutHi = 1.0;
                }
                else if (v == "weta")
                {
                  cutText = "#gamma-ID: 0 < w_{#eta}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
                }
                else if (v == "wphi")
                {
                  cutText = "#gamma-ID: 0 < w_{#phi}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
                }

                if (!cutText.empty())
                {
                  tcut.DrawLatex(0.16, 0.88, cutText.c_str());
                }

                  if (drawCuts)
                  {
                    const double yMinPad = hFirst->GetMinimum();
                    const double yMaxPad = hFirst->GetMaximum();

                    TLine* l1 = new TLine(cutLo, yMinPad, cutLo, yMaxPad);
                    l1->SetLineColor(kBlack);
                    l1->SetLineWidth(3);
                    l1->SetLineStyle(2);
                    l1->Draw("same");
                    keepLines.push_back(l1);

                    TLine* l2 = new TLine(cutHi, yMinPad, cutHi, yMaxPad);
                    l2->SetLineColor(kBlack);
                    l2->SetLineWidth(3);
                    l2->SetLineStyle(2);
                    l2->Draw("same");
                    keepLines.push_back(l2);

                    gPad->RedrawAxis();
                  }

                  if (v == "weta" || v == "wphi")
                  {
                    const double yMinPad = hFirst->GetMinimum();
                    const double yMaxPad = hFirst->GetMaximum();

                    TLine* lw = new TLine(0.3, yMinPad, 0.3, yMaxPad);
                    lw->SetLineColor(kBlack);
                    lw->SetLineWidth(3);
                    lw->SetLineStyle(2);
                    lw->Draw("same");
                    keepLines.push_back(lw);

                    gPad->RedrawAxis();
                  }

                gPad->RedrawAxis();
              }
            }

            for (const auto& vd : vars)
            {
              const string varBase = JoinPath(baseNoIso, vd.var);
              EnsureDir(varBase);

              const string ptDir = JoinPath(varBase, pb.folder);
              EnsureDir(ptDir);

              SaveCanvas(cAll, JoinPath(ptDir, "table2x3_ppAuAu_unNormalized_SS_overlaysByCent.png"));
            }

            for (TLegend* l : keepL) delete l;
            keepL.clear();

            for (TLine* ln : keepLines) delete ln;
            keepLines.clear();

            for (TH1* h : keepH) delete h;
            keepH.clear();
          }
        }
      }
    }

    // ---------------------------------------------------------------------
    // NEW: noIsoRequired/0_10  (AuAu counts) — 1x5 SS pT-overlay summary table
    //
    // Purpose:
    //   For the 0-10% centrality bin, make a single 1x5 canvas where each pad
    //   is one SS variable with ALL pT bins overlaid (same legend in each pad),
    //   including the SAME SS cut label + vertical cut lines used elsewhere.
    //
    // Output:
    //   <kOutPPAuAuBase>/noIsoRequired/0_10/table1x5_AuAu_unNormalized_SS_byPtOverlays.png
    //
    // Ordering:
    //   weta, wphi, e11e33, et1, e32e35
    // Titles per pad:
    //   <SS label>, cent = 0-10%, Run25auau
    // ---------------------------------------------------------------------
    {
      if (aaTop)
      {
        const string baseNoIso = JoinPath(outBase, "noIsoRequired");
        EnsureDir(baseNoIso);

        const auto& centBinsLocal = CentBins();
        const auto& ptBinsLocal   = PtBins();

        const CentBin* cb0 = nullptr;
        for (const auto& cb : centBinsLocal)
        {
          if (cb.lo == 0 && cb.hi == 10) { cb0 = &cb; break; }
        }
        if (!cb0 && !centBinsLocal.empty()) cb0 = &centBinsLocal[0];

        if (cb0 && !ptBinsLocal.empty())
        {
          const string outCent = JoinPath(baseNoIso, cb0->folder);
          EnsureDir(outCent);

          struct VarDef { std::string var; std::string label; };
          const std::vector<VarDef> vars =
          {
            {"weta",   "w_{#eta}"},
            {"wphi",   "w_{#phi}"},
            {"e11e33", "E_{11}/E_{33}"},
            {"et1",    "et1"},
            {"e32e35", "E_{32}/E_{35}"}
          };

          const int colors[] = {
            kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 1,
            kViolet + 1, kAzure + 2, kSpring + 5, kPink + 7, kTeal + 3, kGray + 2
          };
          const int nColors = (int)(sizeof(colors)/sizeof(colors[0]));

          TCanvas c("c_noIso_0_10_SS_byPt_1x5", "c_noIso_0_10_SS_byPt_1x5", 2600, 750);
          c.Divide(5, 1, 0.001, 0.001);

          std::vector<TH1*> keepAlive;
          keepAlive.reserve(vars.size() * ptBinsLocal.size());

          std::vector<TLegend*> keepLeg;
          keepLeg.reserve(vars.size());

          for (int iv = 0; iv < (int)vars.size(); ++iv)
          {
            const std::string& var    = vars[iv].var;
            const std::string& vlabel = vars[iv].label;

            c.cd(iv + 1);
            gPad->SetLeftMargin(0.14);
            gPad->SetRightMargin(0.05);
            gPad->SetBottomMargin(0.14);
            gPad->SetTopMargin(0.18);
            gPad->SetLogy(true);

            const bool isW = (var == "weta" || var == "wphi");

            const string histBaseAuAu = string("h_ss_") + var + string("_inclusive");

            std::vector<TH1*> histsPad;
            histsPad.reserve(ptBinsLocal.size());

            std::vector<std::string> labelsPad;
            labelsPad.reserve(ptBinsLocal.size());

            for (int ipt = 0; ipt < (int)ptBinsLocal.size(); ++ipt)
            {
              const PtBin& pb = ptBinsLocal[ipt];
              const string hAAName = histBaseAuAu + pb.suffix + cb0->suffix;

              TH1* rawAA = GetTH1FromTopDir(aaTop, hAAName);
              if (!rawAA) continue;

              TH1* hAAc = CloneTH1(rawAA,
                TString::Format("aa_noIso_0_10_SS_byPt_%s_%s%s",
                                var.c_str(), pb.folder.c_str(), cb0->suffix.c_str()).Data());
              if (!hAAc) continue;

              EnsureSumw2(hAAc);
              hAAc->GetXaxis()->UnZoom();
              hAAc->SetTitle("");
              hAAc->GetXaxis()->SetTitle(vlabel.c_str());
              hAAc->GetYaxis()->SetTitle("Counts");

              const int col = colors[ipt % nColors];
              StyleOverlayHist(hAAc, col, 20);
              hAAc->SetMarkerStyle(20);

              histsPad.push_back(hAAc);
              labelsPad.push_back(TString::Format("%d-%d GeV", pb.lo, pb.hi).Data());
            }

            if (histsPad.empty())
            {
              DrawMissingPad(TString::Format("%s, cent %d-%d%%", vlabel.c_str(), cb0->lo, cb0->hi).Data());
              continue;
            }

            double yMax = 0.0;
            double minPos = 1e99;
            for (TH1* h : histsPad)
            {
              yMax = std::max(yMax, (double)h->GetMaximum());
              const int nb = h->GetNbinsX();
              for (int ib = 1; ib <= nb; ++ib)
              {
                const double y = h->GetBinContent(ib);
                if (y > 0.0 && y < minPos) minPos = y;
              }
            }
            if (!(minPos < 1e98)) minPos = 0.5;

            histsPad[0]->SetMinimum(minPos * 0.5);
            histsPad[0]->SetMaximum((yMax > 0.0) ? (yMax * 3.0) : 1.0);
            histsPad[0]->Draw("E1");
            for (std::size_t j = 1; j < histsPad.size(); ++j) histsPad[j]->Draw("E1 same");

            // Same pT legend on every pad, but position depends on variable group
            //   weta/wphi: center-right
            //   e11e33/et1/e32e35: middle-left
            TLegend* leg = nullptr;
            if (isW)
              leg = new TLegend(0.55, 0.45, 0.93, 0.80);
            else
              leg = new TLegend(0.16, 0.45, 0.54, 0.80);

            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextFont(42);
            leg->SetTextSize(0.032);
            leg->SetNColumns(2);

            for (std::size_t j = 0; j < histsPad.size(); ++j)
              leg->AddEntry(histsPad[j], labelsPad[j].c_str(), "ep");

            leg->Draw();
            keepLeg.push_back(leg);

            // Centered per-pad title (larger)
            {
              TLatex t;
              t.SetNDC(true);
              t.SetTextFont(42);
              t.SetTextAlign(22);
              t.SetTextSize(0.048);
              t.DrawLatex(0.50, 0.955,
                TString::Format("%s, cent = %d-%d%%, Run25auau", vlabel.c_str(), cb0->lo, cb0->hi).Data());
            }

            // SS cut label (top-left) + cut lines (where applicable) — match existing per-var styling
            {
              TLatex tcut;
              tcut.SetNDC(true);
              tcut.SetTextFont(42);
              tcut.SetTextAlign(13);
              tcut.SetTextSize(0.034);

              bool drawCuts = false;
              double cutLo = 0.0;
              double cutHi = 0.0;
              std::string cutText;

              if (var == "e11e33")
              {
                cutText = "Tight #gamma-ID: 0.4 < #frac{E_{11}}{E_{33}} < 0.98";
                drawCuts = true;
                cutLo = 0.4;
                cutHi = 0.98;
              }
              else if (var == "e32e35")
              {
                cutText = "#gamma-ID: 0.92 < #frac{E_{32}}{E_{35}} < 1.0";
                drawCuts = true;
                cutLo = 0.92;
                cutHi = 1.0;
              }
              else if (var == "et1")
              {
                cutText = "#gamma-ID: 0.9 < et1 < 1.0";
                drawCuts = true;
                cutLo = 0.9;
                cutHi = 1.0;
              }
              else if (var == "weta")
              {
                cutText = "#gamma-ID: 0 < w_{#eta}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
              }
              else if (var == "wphi")
              {
                cutText = "#gamma-ID: 0 < w_{#phi}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
              }

              if (!cutText.empty())
              {
                tcut.DrawLatex(0.16, 0.88, cutText.c_str());
              }

              if (drawCuts)
              {
                gPad->Update();
                const double yMin = gPad->GetUymin();
                const double yMaxPad = gPad->GetUymax();

                TLine* l1 = new TLine(cutLo, yMin, cutLo, yMaxPad);
                l1->SetLineColor(kGreen + 2);
                l1->SetLineWidth(2);
                l1->SetLineStyle(2);
                l1->Draw("same");

                TLine* l2 = new TLine(cutHi, yMin, cutHi, yMaxPad);
                l2->SetLineColor(kOrange + 7);
                l2->SetLineWidth(2);
                l2->SetLineStyle(2);
                l2->Draw("same");

                gPad->RedrawAxis();
              }

              gPad->RedrawAxis();
            }

            for (TH1* h : histsPad) keepAlive.push_back(h);
          }

          SaveCanvas(c, JoinPath(outCent, "table1x5_AuAu_unNormalized_SS_byPtOverlays.png"));

          // -------------------------------------------------------------------
          // NEW: 2x5 table: top row = AuAu 0-10% (same as 1x5), bottom row = PP (Run24pp)
          // Output (same folder as 1x5):
          //   <outCent>/table2x5_AuAu0_10_and_PP_unNormalized_SS_byPtOverlays.png
          //
          // Requirements:
          //   - bottom row: PP unnormalized counts, log-y, pT overlays, NO legend
          //   - titles: top row uses Run25auau; bottom row uses Run24pp
          // -------------------------------------------------------------------
          {
              TCanvas c2("c_noIso_0_10_SS_byPt_2x5", "c_noIso_0_10_SS_byPt_2x5", 2600, 1450);
              c2.Divide(5, 2, 0.001, 0.001);

              std::vector<TH1*> keepAlive2;
              keepAlive2.reserve((std::size_t)vars.size() * (std::size_t)ptBinsLocal.size() * 2);

              std::vector<TLegend*> keepLeg2;
              keepLeg2.reserve((std::size_t)vars.size());

              for (int iv = 0; iv < (int)vars.size(); ++iv)
              {
                const std::string& var    = vars[iv].var;
                const std::string& vlabel = vars[iv].label;

                const bool isW = (var == "weta" || var == "wphi");
                const string histBase = string("h_ss_") + var + string("_inclusive");

                // -----------------------------
                // TOP ROW: AuAu 0-10% by pT overlays (with legend)
                // -----------------------------
                c2.cd(iv + 1);
                gPad->SetLeftMargin(0.14);
                gPad->SetRightMargin(0.05);
                gPad->SetBottomMargin(0.14);
                gPad->SetTopMargin(0.18);
                gPad->SetLogy(true);

                std::vector<TH1*> histsAA;
                histsAA.reserve(ptBinsLocal.size());

                std::vector<std::string> labelsAA;
                labelsAA.reserve(ptBinsLocal.size());

                for (int ipt = 0; ipt < (int)ptBinsLocal.size(); ++ipt)
                {
                  const PtBin& pb = ptBinsLocal[ipt];
                  const string hAAName = histBase + pb.suffix + cb0->suffix;

                  TH1* rawAA = GetTH1FromTopDir(aaTop, hAAName);
                  if (!rawAA) continue;

                  TH1* hAAc = CloneTH1(rawAA,
                    TString::Format("aa_noIso_0_10_SS_byPt_2x5_top_%s_%s%s",
                                    var.c_str(), pb.folder.c_str(), cb0->suffix.c_str()).Data());
                  if (!hAAc) continue;

                  EnsureSumw2(hAAc);
                  hAAc->GetXaxis()->UnZoom();
                  hAAc->SetTitle("");
                  hAAc->GetXaxis()->SetTitle(vlabel.c_str());
                  hAAc->GetYaxis()->SetTitle("Counts");

                  const int col = colors[ipt % nColors];
                  StyleOverlayHist(hAAc, col, 20);
                  hAAc->SetMarkerStyle(20);

                  histsAA.push_back(hAAc);
                  labelsAA.push_back(TString::Format("%d-%d GeV", pb.lo, pb.hi).Data());
                }

                if (histsAA.empty())
                {
                  DrawMissingPad(TString::Format("%s, cent %d-%d%%", vlabel.c_str(), cb0->lo, cb0->hi).Data());
                }
                else
                {
                  double yMaxAA = 0.0;
                  double minPosAA = 1e99;
                  for (TH1* h : histsAA)
                  {
                    yMaxAA = std::max(yMaxAA, (double)h->GetMaximum());
                    const int nb = h->GetNbinsX();
                    for (int ib = 1; ib <= nb; ++ib)
                    {
                      const double y = h->GetBinContent(ib);
                      if (y > 0.0 && y < minPosAA) minPosAA = y;
                    }
                  }
                  if (!(minPosAA < 1e98)) minPosAA = 0.5;

                  histsAA[0]->SetMinimum(minPosAA * 0.5);
                  histsAA[0]->SetMaximum((yMaxAA > 0.0) ? (yMaxAA * 3.0) : 1.0);
                  histsAA[0]->Draw("E1");
                  for (std::size_t j = 1; j < histsAA.size(); ++j) histsAA[j]->Draw("E1 same");

                  TLegend* legTop = nullptr;
                  if (isW) legTop = new TLegend(0.55, 0.45, 0.93, 0.80);
                  else     legTop = new TLegend(0.16, 0.45, 0.54, 0.80);

                  legTop->SetBorderSize(0);
                  legTop->SetFillStyle(0);
                  legTop->SetTextFont(42);
                  legTop->SetTextSize(0.032);
                  legTop->SetNColumns(2);

                  for (std::size_t j = 0; j < histsAA.size(); ++j)
                    legTop->AddEntry(histsAA[j], labelsAA[j].c_str(), "ep");

                  legTop->Draw();
                  keepLeg2.push_back(legTop);

                  // Title (AuAu row)
                  {
                    TLatex t;
                    t.SetNDC(true);
                    t.SetTextFont(42);
                    t.SetTextAlign(22);
                    t.SetTextSize(0.052);
                    t.DrawLatex(0.50, 0.95,
                      TString::Format("%s, cent = %d-%d%%, Run25auau", vlabel.c_str(), cb0->lo, cb0->hi).Data());
                  }

                  // Cut label + cut lines (same behavior as 1x5)
                  {
                    TLatex tcut;
                    tcut.SetNDC(true);
                    tcut.SetTextFont(42);
                    tcut.SetTextAlign(13);
                    tcut.SetTextSize(0.037);

                    bool drawCuts = false;
                    double cutLo = 0.0;
                    double cutHi = 0.0;
                    std::string cutText;

                    if (var == "e11e33")
                    {
                      cutText = "#gamma-ID: 0.4 < #frac{E_{11}}{E_{33}} < 0.98";
                      drawCuts = true;
                      cutLo = 0.4;
                      cutHi = 0.98;
                    }
                    else if (var == "e32e35")
                    {
                      cutText = "#gamma-ID: 0.92 < #frac{E_{32}}{E_{35}} < 1.0";
                      drawCuts = true;
                      cutLo = 0.92;
                      cutHi = 1.0;
                    }
                    else if (var == "et1")
                    {
                      cutText = "#gamma-ID: 0.9 < et1 < 1.0";
                      drawCuts = true;
                      cutLo = 0.9;
                      cutHi = 1.0;
                    }
                    else if (var == "weta")
                    {
                      cutText = "#gamma-ID: 0 < w_{#eta}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
                    }
                    else if (var == "wphi")
                    {
                      cutText = "#gamma-ID: 0 < w_{#phi}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
                    }

                    if (!cutText.empty())
                    {
                      tcut.DrawLatex(0.16, 0.88, cutText.c_str());
                    }

                      if (drawCuts)
                      {
                        gPad->Update();
                        const double yMin = histsAA[0]->GetMinimum();
                        const double yMax = histsAA[0]->GetMaximum();

                        TLine* l1 = new TLine(cutLo, yMin, cutLo, yMax);
                        l1->SetLineColor(kBlack);
                        l1->SetLineWidth(3);
                        l1->SetLineStyle(2);
                        l1->Draw("same");

                        TLine* l2 = new TLine(cutHi, yMin, cutHi, yMax);
                        l2->SetLineColor(kBlack);
                        l2->SetLineWidth(3);
                        l2->SetLineStyle(2);
                        l2->Draw("same");

                        TLegend* legCut = new TLegend(0.62, 0.80, 0.93, 0.90);
                        legCut->SetBorderSize(0);
                        legCut->SetFillStyle(0);
                        legCut->SetTextFont(42);
                        legCut->SetTextSize(0.040);
                        legCut->AddEntry(l1, "cut bounds", "l");
                        legCut->Draw();
                        keepLeg2.push_back(legCut);

                        gPad->RedrawAxis();
                      }

                      if (isW)
                      {
                        gPad->Update();
                        const double yMin = histsAA[0]->GetMinimum();
                        const double yMax = histsAA[0]->GetMaximum();

                        TLine* lw = new TLine(0.3, yMin, 0.3, yMax);
                        lw->SetLineColor(kBlack);
                        lw->SetLineWidth(3);
                        lw->SetLineStyle(2);
                        lw->Draw("same");

                          TLegend* legW = new TLegend(0.16, 0.89, 0.54, 0.93);
                          legW->SetBorderSize(0);
                          legW->SetFillStyle(0);
                          legW->SetTextFont(42);
                          legW->SetTextSize(0.038);
                          legW->AddEntry(lw, "w^{cogX, high} for p_{T}^{#gamma} = 25 GeV", "l");
                          legW->Draw();
                          keepLeg2.push_back(legW);

                        gPad->RedrawAxis();
                      }

                      gPad->RedrawAxis();
                  }
                }

                for (TH1* h : histsAA) keepAlive2.push_back(h);

                // -----------------------------
                // BOTTOM ROW: PP by pT overlays (NO legend)
                // -----------------------------
                c2.cd((iv + 1) + 5);
                gPad->SetLeftMargin(0.14);
                gPad->SetRightMargin(0.05);
                gPad->SetBottomMargin(0.14);
                gPad->SetTopMargin(0.18);
                gPad->SetLogy(true);

                std::vector<TH1*> histsPP;
                histsPP.reserve(ptBinsLocal.size());

                for (int ipt = 0; ipt < (int)ptBinsLocal.size(); ++ipt)
                {
                  const PtBin& pb = ptBinsLocal[ipt];
                  const string hPPName = histBase + pb.suffix;

                  TH1* rawPP = nullptr;
                  rawPP = GetTH1FromTopDir(dsPP.topDir, hPPName);
                  if (!rawPP) continue;

                  TH1* hPPc = CloneTH1(rawPP,
                    TString::Format("pp_SS_byPt_2x5_bot_%s_%s",
                                    var.c_str(), pb.folder.c_str()).Data());
                  if (!hPPc) continue;

                  EnsureSumw2(hPPc);
                  hPPc->GetXaxis()->UnZoom();
                  hPPc->SetTitle("");
                  hPPc->GetXaxis()->SetTitle(vlabel.c_str());
                  hPPc->GetYaxis()->SetTitle("Counts");

                  const int col = colors[ipt % nColors];
                  StyleOverlayHist(hPPc, col, 20);
                  hPPc->SetMarkerStyle(20);

                  histsPP.push_back(hPPc);
                }

                if (histsPP.empty())
                {
                  DrawMissingPad(TString::Format("%s, Run24pp", vlabel.c_str()).Data());
                }
                else
                {
                  double yMaxPP = 0.0;
                  double minPosPP = 1e99;
                  for (TH1* h : histsPP)
                  {
                    yMaxPP = std::max(yMaxPP, (double)h->GetMaximum());
                    const int nb = h->GetNbinsX();
                    for (int ib = 1; ib <= nb; ++ib)
                    {
                      const double y = h->GetBinContent(ib);
                      if (y > 0.0 && y < minPosPP) minPosPP = y;
                    }
                  }
                  if (!(minPosPP < 1e98)) minPosPP = 0.5;

                  histsPP[0]->SetMinimum(minPosPP * 0.5);
                  histsPP[0]->SetMaximum((yMaxPP > 0.0) ? (yMaxPP * 3.0) : 1.0);
                  histsPP[0]->Draw("E1");
                  for (std::size_t j = 1; j < histsPP.size(); ++j) histsPP[j]->Draw("E1 same");

                  // Title (PP row)
                  {
                    TLatex t;
                    t.SetNDC(true);
                    t.SetTextFont(42);
                    t.SetTextAlign(22);
                    t.SetTextSize(0.052);
                    t.DrawLatex(0.50, 0.95,
                      TString::Format("%s, Run24pp", vlabel.c_str()).Data());
                  }

                  // Cut label + cut lines (same as AuAu row)
                  {
                    TLatex tcut;
                    tcut.SetNDC(true);
                    tcut.SetTextFont(42);
                    tcut.SetTextAlign(13);
                    tcut.SetTextSize(0.037);

                    bool drawCuts = false;
                    double cutLo = 0.0;
                    double cutHi = 0.0;
                    std::string cutText;

                    if (var == "e11e33")
                    {
                      cutText = "#gamma-ID: 0.4 < #frac{E_{11}}{E_{33}} < 0.98";
                      drawCuts = true;
                      cutLo = 0.4;
                      cutHi = 0.98;
                    }
                    else if (var == "e32e35")
                    {
                      cutText = "#gamma-ID: 0.92 < #frac{E_{32}}{E_{35}} < 1.0";
                      drawCuts = true;
                      cutLo = 0.92;
                      cutHi = 1.0;
                    }
                    else if (var == "et1")
                    {
                      cutText = "#gamma-ID: 0.9 < et1 < 1.0";
                      drawCuts = true;
                      cutLo = 0.9;
                      cutHi = 1.0;
                    }
                    else if (var == "weta")
                    {
                      cutText = "#gamma-ID: 0 < w_{#eta}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
                    }
                    else if (var == "wphi")
                    {
                      cutText = "#gamma-ID: 0 < w_{#phi}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
                    }

                    if (!cutText.empty())
                    {
                      tcut.DrawLatex(0.16, 0.88, cutText.c_str());
                    }

                    if (drawCuts)
                    {
                      gPad->Update();
                      const double yMin = histsPP[0]->GetMinimum();
                      const double yMax = histsPP[0]->GetMaximum();

                      TLine* l1 = new TLine(cutLo, yMin, cutLo, yMax);
                      l1->SetLineColor(kBlack);
                      l1->SetLineWidth(3);
                      l1->SetLineStyle(2);
                      l1->Draw("same");

                      TLine* l2 = new TLine(cutHi, yMin, cutHi, yMax);
                      l2->SetLineColor(kBlack);
                      l2->SetLineWidth(3);
                      l2->SetLineStyle(2);
                      l2->Draw("same");

                        TLegend* legCut = new TLegend(0.62, 0.80, 0.93, 0.90);
                        legCut->SetBorderSize(0);
                        legCut->SetFillStyle(0);
                        legCut->SetTextFont(42);
                        legCut->SetTextSize(0.040);
                        legCut->AddEntry(l1, "cut bounds", "l");
                        legCut->Draw();
                        keepLeg2.push_back(legCut);

                      gPad->RedrawAxis();
                    }

                    if (isW)
                    {
                      gPad->Update();
                      const double yMin = histsPP[0]->GetMinimum();
                      const double yMax = histsPP[0]->GetMaximum();

                      TLine* lw = new TLine(0.3, yMin, 0.3, yMax);
                      lw->SetLineColor(kBlack);
                      lw->SetLineWidth(3);
                      lw->SetLineStyle(2);
                      lw->Draw("same");

                        TLegend* legW = new TLegend(0.16, 0.89, 0.54, 0.93);
                        legW->SetBorderSize(0);
                        legW->SetFillStyle(0);
                        legW->SetTextFont(42);
                        legW->SetTextSize(0.038);
                        legW->AddEntry(lw, "w^{cogX, high} for p_{T}^{#gamma} = 25 GeV", "l");
                        legW->Draw();
                        keepLeg2.push_back(legW);
                        
                      gPad->RedrawAxis();
                    }

                    gPad->RedrawAxis();
                  }
                }

                for (TH1* h : histsPP) keepAlive2.push_back(h);
              }

              SaveCanvas(c2, JoinPath(outCent, "table2x5_AuAu0_10_and_PP_unNormalized_SS_byPtOverlays.png"));

              for (int ipad = 1; ipad <= 10; ++ipad)
              {
                TPad* p = dynamic_cast<TPad*>(c2.cd(ipad));
                if (!p) continue;
                p->SetLogy(false);
                p->Modified();
              }
              c2.cd();
              c2.Modified();
              c2.Update();

              SaveCanvas(c2, JoinPath(outCent, "table2x5_AuAu0_10_and_PP_unNormalized_SS_byPtOverlays_linearY.png"));

              // -------------------------------------------------------------------
              // NEW: PP-only (noIsoRequired) — per-pT-bin 1x5 SS tables
              //
              // Outputs (per pT folder):
              //   (A) table1x5_PP_SS.png                             (inclusive, raw counts)
              //   (B) table1x5_PP_SS_pre_DataSigBkg.png              (PPG12: preselection, Data vs MC)
              //   (C) table1x5_PP_SS_tight_DataSigBkg.png            (PPG12: tight,       Data vs MC)
              //   (D) table1x5_PP_SS_nonTight_DataSigBkg.png         (PPG12: non-tight,   Data vs MC)
              //
              // Pads (L->R):
              //   weta, wphi, e11e33, et1, e32e35
              //
              // Notes:
              //   - (B-D) expect RecoilJets to have filled:
              //       DATA: h_ss_<var>_<tag>_pT_*
              //       SIM : h_ss_<var>_<tag>_sig_pT_* and ..._bkg_pT_*
              //   - (B-D) use unit-area normalization (shape comparison), NOT log-y
              // -------------------------------------------------------------------
              {
                  const string outNoIso = JoinPath(kOutPPAuAuBase, "noIsoRequired");
                  EnsureDir(outNoIso);

                  // SIM input (for *_sig / *_bkg SS templates)
                  TFile* fSimSS = nullptr;
                  TDirectory* simTopSS = nullptr;
                  {
                    const string simPath = SimInputPathForSample(CurrentSimSample());
                    if (!simPath.empty())
                    {
                      fSimSS = TFile::Open(simPath.c_str(), "READ");
                      if (fSimSS && !fSimSS->IsZombie())
                      {
                        simTopSS = fSimSS->GetDirectory(kDirSIM.c_str());
                        if (!simTopSS)
                        {
                          cout << ANSI_BOLD_YEL
                               << "[WARN] SS templates: missing topDir '" << kDirSIM
                               << "' in SIM file: " << simPath
                               << ANSI_RESET << "\n";
                        }
                      }
                      else
                      {
                        cout << ANSI_BOLD_YEL
                             << "[WARN] SS templates: cannot open SIM file for *_sig/*_bkg: " << simPath
                             << ANSI_RESET << "\n";
                        if (fSimSS) { fSimSS->Close(); delete fSimSS; }
                        fSimSS = nullptr;
                      }
                    }
                  }

                  const std::vector<std::string> ppg12Tags = {"pre", "tight", "nonTight"};

                  for (int ipt = 0; ipt < (int)ptBinsLocal.size(); ++ipt)
                  {
                    const PtBin& pb = ptBinsLocal[ipt];
                    const string outPt = JoinPath(outNoIso, pb.folder);
                    EnsureDir(outPt);

                    // ------------------------------------------------------------
                    // (A) Inclusive PP-only raw-counts table (existing output)
                    // ------------------------------------------------------------
                    {
                      TCanvas cPP(
                        TString::Format("c_noIso_pp_SS_1x5_%s", pb.folder.c_str()).Data(),
                        "c_noIso_pp_SS_1x5", 2600, 750
                      );
                      cPP.Divide(5, 1, 0.001, 0.001);

                      std::vector<TH1*> keepAlivePP;
                      keepAlivePP.reserve(vars.size());

                      for (int iv = 0; iv < (int)vars.size(); ++iv)
                      {
                        const std::string& var    = vars[iv].var;
                        const std::string& vlabel = vars[iv].label;

                        cPP.cd(iv + 1);
                        gPad->SetLeftMargin(0.14);
                        gPad->SetRightMargin(0.05);
                        gPad->SetBottomMargin(0.14);
                        gPad->SetTopMargin(0.18);
                        gPad->SetLogy(false);

                        const string hPPName = string("h_ss_") + var + string("_inclusive") + pb.suffix;
                        TH1* rawPP = GetTH1FromTopDir(dsPP.topDir, hPPName);

                        if (!rawPP)
                        {
                          TLatex tm;
                          tm.SetNDC(true);
                          tm.SetTextFont(42);
                          tm.SetTextAlign(22);
                          tm.SetTextSize(0.055);
                          tm.DrawLatex(0.50, 0.55, "MISSING");
                          continue;
                        }

                        TH1* hPPc = (TH1*)rawPP->Clone(
                          TString::Format("pp_noIso_SS_%s_%s", var.c_str(), pb.folder.c_str()).Data()
                        );
                        if (!hPPc) continue;
                        hPPc->SetDirectory(nullptr);

                        hPPc->SetLineColor(kBlack);
                        hPPc->SetMarkerColor(kBlack);

                        hPPc->SetMarkerStyle(20);
                        hPPc->SetMarkerSize(1.05);
                        hPPc->GetXaxis()->SetTitle(vlabel.c_str());
                        hPPc->GetYaxis()->SetTitle("Counts");
                        hPPc->SetMinimum(0.0);

                        const double yBuf = 1.25;
                        const double ymax = yBuf * hPPc->GetMaximum();
                        if (ymax > 0.0) hPPc->SetMaximum(ymax);

                        hPPc->Draw("E1");

                        // Pad header (SS label + pT bin)
                        {
                          TLatex th;
                          th.SetNDC(true);
                          th.SetTextFont(42);
                          th.SetTextAlign(22);
                          th.SetTextSize(0.055);
                          th.DrawLatex(0.50, 0.94,
                            TString::Format("%s, p_{T}^{#gamma}: %d-%d GeV", vlabel.c_str(), pb.lo, pb.hi).Data());
                        }

                        // ------------------------------------------------------------
                        // Add SS cut label (top-left) + cut lines (where applicable)
                        // ------------------------------------------------------------
                        {
                          TLatex tcut;
                          tcut.SetNDC(true);
                          tcut.SetTextFont(42);
                          tcut.SetTextAlign(13);
                          tcut.SetTextSize(0.045);

                          bool drawCuts = false;
                          double cutLo = 0.0;
                          double cutHi = 0.0;
                          std::string cutText;

                          if (var == "e11e33")
                          {
                            cutText = "Tight #gamma-ID: 0.4 < #frac{E_{11}}{E_{33}} < 0.98";
                            drawCuts = true;
                            cutLo = 0.4;
                            cutHi = 0.98;
                          }
                          else if (var == "e32e35")
                          {
                            cutText = "#gamma-ID: 0.92 < #frac{E_{32}}{E_{35}} < 1.0";
                            drawCuts = true;
                            cutLo = 0.92;
                            cutHi = 1.0;
                          }
                          else if (var == "et1")
                          {
                            cutText = "#gamma-ID: 0.9 < et1 < 1.0";
                            drawCuts = true;
                            cutLo = 0.9;
                            cutHi = 1.0;
                          }
                          else if (var == "weta")
                          {
                            cutText = "#gamma-ID: 0 < w_{#eta}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
                          }
                          else if (var == "wphi")
                          {
                            cutText = "#gamma-ID: 0 < w_{#phi}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
                          }

                          if (!cutText.empty())
                          {
                            tcut.DrawLatex(0.16, 0.86, cutText.c_str());
                          }

                          if (drawCuts)
                          {
                            gPad->Update();
                            const double yMin = gPad->GetUymin();
                            const double yMaxPad = gPad->GetUymax();

                            TLine* l1 = new TLine(cutLo, yMin, cutLo, yMaxPad);
                            l1->SetLineColor(kGreen + 2);
                            l1->SetLineWidth(2);
                            l1->SetLineStyle(2);
                            l1->Draw("same");

                            TLine* l2 = new TLine(cutHi, yMin, cutHi, yMaxPad);
                            l2->SetLineColor(kOrange + 7);
                            l2->SetLineWidth(2);
                            l2->SetLineStyle(2);
                            l2->Draw("same");
                          }

                          gPad->RedrawAxis();
                        }

                        keepAlivePP.push_back(hPPc);
                      }

                      SaveCanvas(cPP, JoinPath(outPt, "table1x5_PP_SS.png"));

                      for (TH1* h : keepAlivePP) delete h;
                      keepAlivePP.clear();
                    }

                    // ------------------------------------------------------------
                    // (B-D) PPG12-style tables: Data vs Signal MC vs Background MC
                    // ------------------------------------------------------------
                    for (const auto& tag : ppg12Tags)
                    {
                      TCanvas cPP(
                        TString::Format("c_noIso_pp_SS_1x5_%s_%s", tag.c_str(), pb.folder.c_str()).Data(),
                        "c_noIso_pp_SS_1x5", 2600, 750
                      );
                      cPP.Divide(5, 1, 0.001, 0.001);

                      std::vector<TH1*> keepAlivePP;
                      keepAlivePP.reserve(vars.size() * 3);

                      std::vector<TLegend*> keepLegPP;
                      keepLegPP.reserve(vars.size());

                      bool anyPad = false;

                      for (int iv = 0; iv < (int)vars.size(); ++iv)
                      {
                        const std::string& var    = vars[iv].var;
                        const std::string& vlabel = vars[iv].label;

                        cPP.cd(iv + 1);
                        gPad->SetLeftMargin(0.14);
                        gPad->SetRightMargin(0.05);
                        gPad->SetBottomMargin(0.14);
                        gPad->SetTopMargin(0.18);
                        gPad->SetLogy(false);

                        const bool isW = (var == "weta" || var == "wphi");

                        const string hDataName = string("h_ss_") + var + string("_") + tag + pb.suffix;
                        TH1* rawData = GetTH1FromTopDir(dsPP.topDir, hDataName);

                        TH1* rawSig = nullptr;
                        TH1* rawBkg = nullptr;

                        if (simTopSS)
                        {
                          const string hSigName = string("h_ss_") + var + string("_") + tag + string("_sig") + pb.suffix;
                          const string hBkgName = string("h_ss_") + var + string("_") + tag + string("_bkg") + pb.suffix;
                          rawSig = GetTH1FromTopDir(simTopSS, hSigName);
                          rawBkg = GetTH1FromTopDir(simTopSS, hBkgName);
                        }

                        if (!rawData && !rawSig && !rawBkg)
                        {
                          TLatex tm;
                          tm.SetNDC(true);
                          tm.SetTextFont(42);
                          tm.SetTextAlign(22);
                          tm.SetTextSize(0.055);
                          tm.DrawLatex(0.50, 0.55, "MISSING");
                          continue;
                        }

                        anyPad = true;

                        TH1* hData = nullptr;
                        TH1* hSig  = nullptr;
                        TH1* hBkg  = nullptr;

                        if (rawData)
                          hData = CloneNormalizeStyle(rawData,
                            TString::Format("pp_%s_%s_%s_data", tag.c_str(), var.c_str(), pb.folder.c_str()).Data(),
                            kBlack, 20);

                        if (rawSig)
                          hSig = CloneNormalizeStyle(rawSig,
                            TString::Format("pp_%s_%s_%s_sig", tag.c_str(), var.c_str(), pb.folder.c_str()).Data(),
                            kRed + 1, 24);

                        if (rawBkg)
                          hBkg = CloneNormalizeStyle(rawBkg,
                            TString::Format("pp_%s_%s_%s_bkg", tag.c_str(), var.c_str(), pb.folder.c_str()).Data(),
                            kBlue + 1, 25);

                        TH1* hFirst = (hData ? hData : (hSig ? hSig : hBkg));
                        if (!hFirst)
                        {
                          TLatex tm;
                          tm.SetNDC(true);
                          tm.SetTextFont(42);
                          tm.SetTextAlign(22);
                          tm.SetTextSize(0.055);
                          tm.DrawLatex(0.50, 0.55, "MISSING");
                          continue;
                        }

                        hFirst->GetXaxis()->SetTitle(vlabel.c_str());
                        hFirst->GetYaxis()->SetTitle("Normalized counts");

                        double yMax = 0.0;
                        if (hData) yMax = std::max(yMax, (double)hData->GetMaximum());
                        if (hSig)  yMax = std::max(yMax, (double)hSig->GetMaximum());
                        if (hBkg)  yMax = std::max(yMax, (double)hBkg->GetMaximum());

                        hFirst->SetMinimum(0.0);
                        hFirst->SetMaximum((yMax > 0.0) ? (yMax * 1.35) : 1.0);

                        hFirst->Draw("E1");
                        if (hData && hData != hFirst) hData->Draw("E1 same");
                        if (hSig  && hSig  != hFirst) hSig->Draw("E1 same");
                        if (hBkg  && hBkg  != hFirst) hBkg->Draw("E1 same");

                        // Legend: Data / Signal MC / Background MC
                        {
                          TLegend* leg = nullptr;
                          if (isW) leg = new TLegend(0.55, 0.62, 0.93, 0.86);
                          else     leg = new TLegend(0.16, 0.62, 0.54, 0.86);

                          leg->SetBorderSize(0);
                          leg->SetFillStyle(0);
                          leg->SetTextFont(42);
                          leg->SetTextSize(0.036);

                          if (hData) leg->AddEntry(hData, "Data", "ep");
                          if (hSig)  leg->AddEntry(hSig,  "Signal MC", "ep");
                          if (hBkg)  leg->AddEntry(hBkg,  "Background MC", "ep");

                          leg->Draw();
                          keepLegPP.push_back(leg);
                        }

                        // Pad header (SS label + category + pT bin)
                        {
                          std::string tagLabel = "Preselection";
                          if (tag == "tight") tagLabel = "Tight";
                          else if (tag == "nonTight") tagLabel = "Non-tight";

                          TLatex th;
                          th.SetNDC(true);
                          th.SetTextFont(42);
                          th.SetTextAlign(22);
                          th.SetTextSize(0.050);
                          th.DrawLatex(0.50, 0.94,
                            TString::Format("%s, %s, p_{T}^{#gamma}: %d-%d GeV",
                              vlabel.c_str(), tagLabel.c_str(), pb.lo, pb.hi).Data());
                        }

                        // ------------------------------------------------------------
                        // Add SS cut label (top-left) + cut lines (where applicable)
                        // ------------------------------------------------------------
                        {
                          TLatex tcut;
                          tcut.SetNDC(true);
                          tcut.SetTextFont(42);
                          tcut.SetTextAlign(13);
                          tcut.SetTextSize(0.040);

                          bool drawCuts = false;
                          double cutLo = 0.0;
                          double cutHi = 0.0;
                          std::string cutText;

                          if (var == "e11e33")
                          {
                            cutText = "Tight #gamma-ID: 0.4 < #frac{E_{11}}{E_{33}} < 0.98";
                            drawCuts = true;
                            cutLo = 0.4;
                            cutHi = 0.98;
                          }
                          else if (var == "e32e35")
                          {
                            cutText = "#gamma-ID: 0.92 < #frac{E_{32}}{E_{35}} < 1.0";
                            drawCuts = true;
                            cutLo = 0.92;
                            cutHi = 1.0;
                          }
                          else if (var == "et1")
                          {
                            cutText = "#gamma-ID: 0.9 < et1 < 1.0";
                            drawCuts = true;
                            cutLo = 0.9;
                            cutHi = 1.0;
                          }
                          else if (var == "weta")
                          {
                            cutText = "#gamma-ID: 0 < w_{#eta}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
                          }
                          else if (var == "wphi")
                          {
                            cutText = "#gamma-ID: 0 < w_{#phi}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
                          }

                          if (!cutText.empty())
                          {
                            tcut.DrawLatex(0.16, 0.86, cutText.c_str());
                          }

                          gPad->Update();
                          const double yMin = gPad->GetUymin();
                          const double yMaxPad = gPad->GetUymax();

                          if (drawCuts)
                          {
                            TLine* l1 = new TLine(cutLo, yMin, cutLo, yMaxPad);
                            l1->SetLineColor(kGreen + 2);
                            l1->SetLineWidth(2);
                            l1->SetLineStyle(2);
                            l1->Draw("same");

                            TLine* l2 = new TLine(cutHi, yMin, cutHi, yMaxPad);
                            l2->SetLineColor(kOrange + 7);
                            l2->SetLineWidth(2);
                            l2->SetLineStyle(2);
                            l2->Draw("same");
                          }

                          // For weta/wphi: show the pT=25 GeV "high" line used elsewhere
                          if (isW)
                          {
                            TLine* lw = new TLine(0.3, yMin, 0.3, yMaxPad);
                            lw->SetLineColor(kBlack);
                            lw->SetLineWidth(2);
                            lw->SetLineStyle(2);
                            lw->Draw("same");
                          }

                          gPad->RedrawAxis();
                        }

                        if (hData) keepAlivePP.push_back(hData);
                        if (hSig)  keepAlivePP.push_back(hSig);
                        if (hBkg)  keepAlivePP.push_back(hBkg);
                      }

                      if (anyPad)
                      {
                        SaveCanvas(cPP, JoinPath(outPt,
                          TString::Format("table1x5_PP_SS_%s_DataSigBkg.png", tag.c_str()).Data()));
                      }

                      for (TLegend* l : keepLegPP) delete l;
                      keepLegPP.clear();

                      for (TH1* h : keepAlivePP) delete h;
                      keepAlivePP.clear();
                    }
                  }

                  if (fSimSS) { fSimSS->Close(); delete fSimSS; fSimSS = nullptr; }
                }

              for (TLegend* l : keepLeg2) delete l;
              keepLeg2.clear();

              for (TH1* h : keepAlive2) delete h;
              keepAlive2.clear();
            }

            for (TLegend* l : keepLeg) delete l;
            keepLeg.clear();

            for (TH1* h : keepAlive) delete h;
            keepAlive.clear();
        }
     }
  }

  // ---------------------------------------------------------------------
  // NEW: noSS_isoSpectra summary (AuAu counts by centrality per pT bin)
  // Output:
  //   <kOutPPAuAuBase>/noSS_isoSpectra/table2x3_AuAu_unNormalized.png
  // ---------------------------------------------------------------------
  {
      const string outNoSS = JoinPath(outBase, "noSS_isoSpectra");
      EnsureDir(outNoSS);

      const string histBaseIso = "h_Eiso";
      const string xTitleIso   = "E_{T}^{iso} [GeV]";

      Make2x3Table_AuAuUnNormalized_ByCent(aaTop, outNoSS, histBaseIso, xTitleIso);
      Make2x3Table_AuAuUnNormalized_ByPtOverlaysPerCent(aaTop, outNoSS, histBaseIso, xTitleIso);

      // -------------------------------------------------------------------
      // NEW: Terminal numeric summaries for the AuAu noSS isolation QA block
      //
      // (1) table2x3_AuAu_unNormalized_byPtOverlays.png
      //     -> h_Eiso shape summary (per centrality, per pT bin):
      //        xMin, xMax, yMax, and x-bin location of yMax
      //
      // (2) Iso-decision acceptance comparison (no UE vs UE-sub)
      //     -> h_isoDecision_*  (bin 1 = PASS, bin 2 = FAIL)
      //        Print PASS / FAIL / pass fraction / gain in PASS yield
      // -------------------------------------------------------------------
      auto PrintByPtOverlaySummary = [&](TDirectory* top,
                                         const std::string& tag,
                                         const std::string& filePath)
      {
        if (!top) return;

        const auto& ptBins   = PtBins();
        const auto& centBins = CentBins();

        const int nPads = std::min(6, (int)centBins.size());
        const int nPt   = (int)ptBins.size();

        std::ios::fmtflags f = cout.flags();
        std::streamsize    p = cout.precision();

        cout << ANSI_BOLD_CYN
             << "\n[AuAu SUMMARY] table2x3_AuAu_unNormalized_byPtOverlays.png  (" << tag << ")\n"
             << ANSI_RESET;
        cout << "  file: " << filePath << "\n";
        cout << "  histBase: " << histBaseIso << "   xTitle: " << xTitleIso << "\n";
        cout << "  pads: first " << nPads << " centrality bins (matches 2x3 table)\n\n";

        cout << ANSI_BOLD_WHT;
        cout << std::left
             << std::setw(10) << "cent"
             << std::setw(10) << "pTgamma"
             << std::setw(10) << "xMin"
             << std::setw(10) << "xMax"
             << std::setw(14) << "yMax"
             << std::setw(10) << "x@Max"
             << std::setw(12) << "binLo"
             << std::setw(12) << "binHi"
             << ANSI_RESET << "\n";

        cout << std::string(88, '-') << "\n";

        for (int ic = 0; ic < nPads; ++ic)
        {
          const auto& cb = centBins[ic];
          const std::string centStr = TString::Format("%d-%d%%", cb.lo, cb.hi).Data();

          for (int ipt = 0; ipt < nPt; ++ipt)
          {
            const PtBin& pb = ptBins[ipt];
            const std::string ptStr = TString::Format("%d-%d", pb.lo, pb.hi).Data();

            const string hName = histBaseIso + pb.suffix + cb.suffix;
            TH1* h = GetTH1FromTopDir(top, hName);

            if (!h)
            {
              cout << std::left
                   << std::setw(10) << centStr
                   << std::setw(10) << ptStr
                   << ANSI_BOLD_RED << "MISSING" << ANSI_RESET << "\n";
              continue;
            }

            const int nb = h->GetNbinsX();

            const double xMin = h->GetXaxis()->GetBinLowEdge(1);
            const double xMax = h->GetXaxis()->GetBinUpEdge(nb);

            const int ibMax = h->GetMaximumBin();
            const double yMax = h->GetBinContent(ibMax);

            const double xAtMax = h->GetXaxis()->GetBinCenter(ibMax);
            const double binLo  = h->GetXaxis()->GetBinLowEdge(ibMax);
            const double binHi  = h->GetXaxis()->GetBinUpEdge(ibMax);

            cout << std::left
                 << std::setw(10) << centStr
                 << std::setw(10) << ptStr
                 << std::setw(10) << std::fixed << std::setprecision(3) << xMin
                 << std::setw(10) << std::fixed << std::setprecision(3) << xMax
                 << std::setw(14) << std::fixed << std::setprecision(3) << yMax
                 << std::setw(10) << std::fixed << std::setprecision(3) << xAtMax
                 << std::setw(12) << std::fixed << std::setprecision(3) << binLo
                 << std::setw(12) << std::fixed << std::setprecision(3) << binHi
                 << "\n";
          }
        }

        cout << "\n";
        cout.flags(f);
        cout.precision(p);
      };

      auto PrintIsoDecisionSummary = [&](TDirectory* topNoUE,
                                         TDirectory* topUE,
                                         const std::string& fileNoUE,
                                         const std::string& fileUE)
      {
        if (!topNoUE || !topUE) return;

        const auto& ptBins   = PtBins();
        const auto& centBins = CentBins();

        const int nPads = std::min(6, (int)centBins.size());
        const int nPt   = (int)ptBins.size();
        const string histBaseIsoDecision = "h_isoDecision";

        struct IsoCounts
        {
          bool ok = false;
          double nPass = 0.0;
          double nFail = 0.0;
          double nTot  = 0.0;
          double fPass = 0.0;
        };

        struct GainEntry
        {
          double gain = 0.0;
          std::string cent;
          std::string pt;
          double passNoUE = 0.0;
          double passUE = 0.0;
          double fNoUE = 0.0;
          double fUE = 0.0;
        };

        auto ReadIsoCounts = [&](TDirectory* top, const std::string& hName) -> IsoCounts
        {
          IsoCounts out;
          if (!top) return out;

          TH1* h = GetTH1FromTopDir(top, hName);
          if (!h) return out;

          out.ok = true;
          out.nPass = h->GetBinContent(1);
          out.nFail = h->GetBinContent(2);
          out.nTot  = out.nPass + out.nFail;
          out.fPass = (out.nTot > 0.0 ? out.nPass / out.nTot : 0.0);
          return out;
        };

        std::vector<GainEntry> gainEntries;
        std::vector<GainEntry> revivedEntries;

        std::ios::fmtflags f = cout.flags();
        std::streamsize    p = cout.precision();

        cout << ANSI_BOLD_CYN
             << "\n[AuAu ISO ACCEPTANCE SUMMARY] with / without UE subtraction\n"
             << ANSI_RESET;
        cout << "  Histogram base : " << histBaseIsoDecision << "   (bin 1 = PASS, bin 2 = FAIL)\n";
        cout << "  Quantity       : photons passing isolation cut, independent of tight photon ID\n";
        cout << "  File (no UE)   : " << fileNoUE << "\n";
        cout << "  File (UE-sub)  : " << fileUE << "\n";
        cout << "  Pads compared  : first " << nPads << " centrality bins (matches 2x3 layout)\n";

        for (int ic = 0; ic < nPads; ++ic)
        {
          const auto& cb = centBins[ic];
          const std::string centStr = TString::Format("%d-%d%%", cb.lo, cb.hi).Data();

          cout << "\n" << ANSI_BOLD_WHT
               << "==============================================================================================\n"
               << "CENTRALITY = " << centStr << "\n"
               << "=============================================================================================="
               << ANSI_RESET << "\n";

          cout << ANSI_BOLD_WHT
               << std::left
               << std::setw(10) << "pTgamma"
               << std::setw(12) << "PASS(noUE)"
               << std::setw(12) << "FAIL(noUE)"
               << std::setw(14) << "fpass(noUE)"
               << std::setw(12) << "PASS(UE)"
               << std::setw(12) << "FAIL(UE)"
               << std::setw(12) << "fpass(UE)"
               << std::setw(12) << "Gain(PASS)"
               << std::setw(12) << "Delta fpass"
               << ANSI_RESET << "\n";

          cout << std::string(108, '-') << "\n";

          double sumPassNoUE = 0.0;
          double sumFailNoUE = 0.0;
          double sumPassUE   = 0.0;
          double sumFailUE   = 0.0;

          for (int ipt = 0; ipt < nPt; ++ipt)
          {
            const PtBin& pb = ptBins[ipt];
            const std::string ptStr = TString::Format("%d-%d", pb.lo, pb.hi).Data();

            const string hName = histBaseIsoDecision + pb.suffix + cb.suffix;

            const IsoCounts cNoUE = ReadIsoCounts(topNoUE, hName);
            const IsoCounts cUE   = ReadIsoCounts(topUE,   hName);

            const std::string sPassNoUE = cNoUE.ok ? std::string(TString::Format("%.0f", cNoUE.nPass).Data()) : "MISSING";
            const std::string sFailNoUE = cNoUE.ok ? std::string(TString::Format("%.0f", cNoUE.nFail).Data()) : "MISSING";
            const std::string sFNoUE    = cNoUE.ok ? std::string(TString::Format("%.3f", cNoUE.fPass).Data()) : "--";

            const std::string sPassUE = cUE.ok ? std::string(TString::Format("%.0f", cUE.nPass).Data()) : "MISSING";
            const std::string sFailUE = cUE.ok ? std::string(TString::Format("%.0f", cUE.nFail).Data()) : "MISSING";
            const std::string sFUE    = cUE.ok ? std::string(TString::Format("%.3f", cUE.fPass).Data()) : "--";

            std::string sGain = "--";
            std::string sDeltaF = "--";

            if (cNoUE.ok && cUE.ok)
            {
              sumPassNoUE += cNoUE.nPass;
              sumFailNoUE += cNoUE.nFail;
              sumPassUE   += cUE.nPass;
              sumFailUE   += cUE.nFail;

              sDeltaF = TString::Format("%+.3f", cUE.fPass - cNoUE.fPass).Data();

              if (cNoUE.nPass > 0.0)
              {
                const double gain = cUE.nPass / cNoUE.nPass;
                sGain = TString::Format("%.3f", gain).Data();

                GainEntry ge;
                ge.gain = gain;
                ge.cent = centStr;
                ge.pt = ptStr;
                ge.passNoUE = cNoUE.nPass;
                ge.passUE = cUE.nPass;
                ge.fNoUE = cNoUE.fPass;
                ge.fUE = cUE.fPass;
                gainEntries.push_back(ge);
              }
              else if (cUE.nPass > 0.0)
              {
                sGain = "INF";

                GainEntry ge;
                ge.gain = 1e99;
                ge.cent = centStr;
                ge.pt = ptStr;
                ge.passNoUE = cNoUE.nPass;
                ge.passUE = cUE.nPass;
                ge.fNoUE = cNoUE.fPass;
                ge.fUE = cUE.fPass;
                revivedEntries.push_back(ge);
              }
            }

            cout << std::left
                 << std::setw(10) << ptStr
                 << std::setw(12) << sPassNoUE
                 << std::setw(12) << sFailNoUE
                 << std::setw(14) << sFNoUE
                 << std::setw(12) << sPassUE
                 << std::setw(12) << sFailUE
                 << std::setw(12) << sFUE
                 << std::setw(12) << sGain
                 << std::setw(12) << sDeltaF
                 << "\n";
          }

          const double sumTotNoUE = sumPassNoUE + sumFailNoUE;
          const double sumTotUE   = sumPassUE   + sumFailUE;

          const double fTotNoUE = (sumTotNoUE > 0.0 ? sumPassNoUE / sumTotNoUE : 0.0);
          const double fTotUE   = (sumTotUE   > 0.0 ? sumPassUE   / sumTotUE   : 0.0);

          std::string sGainTot = "--";
          if (sumPassNoUE > 0.0) sGainTot = TString::Format("%.3f", sumPassUE / sumPassNoUE).Data();
          else if (sumPassUE > 0.0) sGainTot = "INF";

          cout << std::string(108, '-') << "\n";
          cout << ANSI_BOLD_YEL
               << std::left
               << std::setw(10) << "Subtotal"
               << std::setw(12) << std::string(TString::Format("%.0f", sumPassNoUE).Data())
               << std::setw(12) << std::string(TString::Format("%.0f", sumFailNoUE).Data())
               << std::setw(14) << std::string(TString::Format("%.3f", fTotNoUE).Data())
               << std::setw(12) << std::string(TString::Format("%.0f", sumPassUE).Data())
               << std::setw(12) << std::string(TString::Format("%.0f", sumFailUE).Data())
               << std::setw(12) << std::string(TString::Format("%.3f", fTotUE).Data())
               << std::setw(12) << sGainTot
               << std::setw(12) << std::string(TString::Format("%+.3f", fTotUE - fTotNoUE).Data())
               << ANSI_RESET << "\n";
        }

        if (!gainEntries.empty())
        {
          std::sort(gainEntries.begin(), gainEntries.end(),
                    [](const GainEntry& a, const GainEntry& b)
                    {
                      return a.gain > b.gain;
                    });

          const int nShow = std::min(3, (int)gainEntries.size());

          cout << "\n" << ANSI_BOLD_GRN
               << "[Top " << nShow << " finite-gain bins by PASS(UE) / PASS(noUE)]"
               << ANSI_RESET << "\n";
          cout << ANSI_BOLD_WHT
               << std::left
               << std::setw(10) << "cent"
               << std::setw(10) << "pTgamma"
               << std::setw(14) << "PASS(noUE)"
               << std::setw(12) << "PASS(UE)"
               << std::setw(14) << "fpass(noUE)"
               << std::setw(12) << "fpass(UE)"
               << std::setw(12) << "Gain"
               << ANSI_RESET << "\n";
          cout << std::string(84, '-') << "\n";

          for (int i = 0; i < nShow; ++i)
          {
            const auto& ge = gainEntries[i];
            cout << std::left
                 << std::setw(10) << ge.cent
                 << std::setw(10) << ge.pt
                 << std::setw(14) << std::string(TString::Format("%.0f", ge.passNoUE).Data())
                 << std::setw(12) << std::string(TString::Format("%.0f", ge.passUE).Data())
                 << std::setw(14) << std::string(TString::Format("%.3f", ge.fNoUE).Data())
                 << std::setw(12) << std::string(TString::Format("%.3f", ge.fUE).Data())
                 << std::setw(12) << std::string(TString::Format("%.3f", ge.gain).Data())
                 << "\n";
          }
        }

        if (!revivedEntries.empty())
        {
          cout << "\n" << ANSI_BOLD_GRN
               << "[Bins revived by UE subtraction: PASS(noUE) = 0 and PASS(UE) > 0]"
               << ANSI_RESET << "\n";
          cout << ANSI_BOLD_WHT
               << std::left
               << std::setw(10) << "cent"
               << std::setw(10) << "pTgamma"
               << std::setw(14) << "PASS(noUE)"
               << std::setw(12) << "PASS(UE)"
               << std::setw(14) << "fpass(noUE)"
               << std::setw(12) << "fpass(UE)"
               << ANSI_RESET << "\n";
          cout << std::string(72, '-') << "\n";

          for (const auto& ge : revivedEntries)
          {
            cout << std::left
                 << std::setw(10) << ge.cent
                 << std::setw(10) << ge.pt
                 << std::setw(14) << std::string(TString::Format("%.0f", ge.passNoUE).Data())
                 << std::setw(12) << std::string(TString::Format("%.0f", ge.passUE).Data())
                 << std::setw(14) << std::string(TString::Format("%.3f", ge.fNoUE).Data())
                 << std::setw(12) << std::string(TString::Format("%.3f", ge.fUE).Data())
                 << "\n";
          }
        }

        cout << "\n";
        cout.flags(f);
        cout.precision(p);
      };

      PrintByPtOverlaySummary(aaTop, "non UE-subtracted", kInAuAuGold);

      if (aaTopNew)
      {
        PrintByPtOverlaySummary(aaTopNew, "UE-subtracted", kInAuAuGoldNew);
      }
      else
      {
        cout << ANSI_BOLD_YEL
             << "[WARN] UE-subtracted Eiso-shape summary skipped (aaTopNew is null): " << kInAuAuGoldNew
             << ANSI_RESET << "\n";
      }

      if (aaTop && aaTopNew)
      {
        PrintIsoDecisionSummary(aaTop, aaTopNew, kInAuAuGold, kInAuAuGoldNew);

        // -------------------------------------------------------------------
        // NEW (presentation-ready): 0-10% isolation pass fraction vs pTgamma
        // Output:
        //   <outBase>/noSS_isoSpectra/0_10/pho_isoPassFraction_vs_pTgamma_noUE_vs_UEsub.png
        //
        // y = f_pass = PASS/(PASS+FAIL), with binomial uncertainty:
        //   sigma_f = sqrt( f(1-f) / Ntot )
        // -------------------------------------------------------------------
        {
          const auto& ptBins   = PtBins();
          const auto& centBins = CentBins();

          if (!centBins.empty())
          {
            const auto& cb = centBins[0]; // 0-10% expected first bin
            const string centDir = JoinPath(outNoSS, TString::Format("%d_%d", cb.lo, cb.hi).Data());
            EnsureDir(centDir);

            const string histBaseIsoDecision = "h_isoDecision";

            vector<double> x, ex, yNoUE, eyNoUE, yUE, eyUE;

            for (int ipt = 0; ipt < (int)ptBins.size(); ++ipt)
            {
              const PtBin& pb = ptBins[ipt];

              const double xcen = 0.5 * ((double)pb.lo + (double)pb.hi);
              const double xerr = 0.5 * ((double)pb.hi - (double)pb.lo);

              const string hName = histBaseIsoDecision + pb.suffix + cb.suffix;

              TH1* hNo = GetTH1FromTopDir(aaTop,    hName);
              TH1* hU  = GetTH1FromTopDir(aaTopNew, hName);

              if (!hNo || !hU) continue;

              const double passNo  = hNo->GetBinContent(1);
              const double failNo  = hNo->GetBinContent(2);
              const double totNo   = passNo + failNo;
              const double fNo     = (totNo > 0.0 ? passNo / totNo : 0.0);
              const double efNo    = (totNo > 0.0 ? std::sqrt(fNo * (1.0 - fNo) / totNo) : 0.0);

              const double passUE  = hU->GetBinContent(1);
              const double failUE  = hU->GetBinContent(2);
              const double totUE   = passUE + failUE;
              const double fUEv    = (totUE > 0.0 ? passUE / totUE : 0.0);
              const double efUE    = (totUE > 0.0 ? std::sqrt(fUEv * (1.0 - fUEv) / totUE) : 0.0);

              x.push_back(xcen);
              ex.push_back(xerr);

              yNoUE.push_back(fNo);
              eyNoUE.push_back(efNo);

              yUE.push_back(fUEv);
              eyUE.push_back(efUE);
            }

            if (!x.empty())
            {
              TCanvas c("c_pho_isoPassFrac_0_10", "c_pho_isoPassFrac_0_10", 950, 700);
              ApplyCanvasMargins1D(c);

              TH1F frame("frame_isoPassFrac", "", 1, 12.0, 36.0);
              frame.SetMinimum(0.0);
              frame.SetMaximum(1.05);
              frame.SetTitle("");
              frame.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV] (bin centers)");
              frame.GetYaxis()->SetTitle("Isolation pass fraction  f_{pass} = PASS/(PASS+FAIL)");
              frame.Draw("axis");

              TGraphErrors gNo((int)x.size(), &x[0], &yNoUE[0], &ex[0], &eyNoUE[0]);
              gNo.SetMarkerStyle(20);
              gNo.SetMarkerSize(1.1);
              gNo.SetMarkerColor(kBlack);
              gNo.SetLineColor(kBlack);
              gNo.SetLineWidth(2);
              gNo.Draw("P same");

              TGraphErrors gU((int)x.size(), &x[0], &yUE[0], &ex[0], &eyUE[0]);
              gU.SetMarkerStyle(24);
              gU.SetMarkerSize(1.1);
              gU.SetMarkerColor(kRed + 1);
              gU.SetLineColor(kRed + 1);
              gU.SetLineWidth(2);
              gU.Draw("P same");

              // Reference line at f_pass = 0.5
              {
                TLine l(12.0, 0.5, 36.0, 0.5);
                l.SetLineStyle(2);
                l.SetLineWidth(2);
                l.Draw("same");
              }

              TLegend leg(0.55, 0.52, 0.88, 0.67);
              leg.SetBorderSize(0);
              leg.SetFillStyle(0);
              leg.SetTextFont(42);
              leg.SetTextSize(0.035);
              leg.AddEntry(&gNo, "no UE subtraction", "pe");
              leg.AddEntry(&gU,  "UE-subtracted iso", "pe");
              leg.Draw();

              // Centered title (presentation-ready)
              {
                TLatex tx;
                tx.SetNDC();
                tx.SetTextFont(42);
                tx.SetTextAlign(22);
                tx.SetTextSize(0.042);
                tx.DrawLatex(0.50, 0.965,
                             TString::Format("AuAu isolation acceptance vs p_{T}^{#gamma} (cent = %d-%d%%)", cb.lo, cb.hi).Data());
              }

              SaveCanvas(c, JoinPath(centDir, "pho_isoPassFraction_vs_pTgamma_noUE_vs_UEsub.png"));
            }
          }
        }
      }
      else
      {
        cout << ANSI_BOLD_YEL
             << "[WARN] Iso-decision acceptance summary skipped because one or both AuAu top directories are null."
             << ANSI_RESET << "\n";
      }
  }

  cout << ANSI_BOLD_GRN
       << "  -> Wrote PP vs AuAu deliverables under: " << outBase
       << ANSI_RESET << "\n";

  if (fAANew)
  {
    fAANew->Close();
    delete fAANew;
  }

  fAA->Close();
  delete fAA;
}

