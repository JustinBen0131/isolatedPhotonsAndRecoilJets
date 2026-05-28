#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace
{
const std::string kCfg = "preselectionReference_tightReference_nonTightReference_baseVariant";
const std::string kOutDir = "dataOutput/diagnostics/stitch_boundary_scan_20260520/correctness_scan";

struct Slice
{
  std::string key;
  std::string label;
  std::string path;
  std::string histName;
  double nominalLo = 0.0;
  double nominalHi = 1.0e9;
  double sigmaNominalPb = 0.0;
  int color = kBlack;
  int marker = 20;
};

struct Window
{
  int slice = -1;
  std::string label;
  double lo = 0.0;
  double hi = 1.0e9;
};

struct Scenario
{
  std::string key;
  std::string label;
  std::string shortLabel;
  std::vector<Window> windows;
  std::vector<double> boundaries;
  bool productionExact = false;
};

struct BuiltWindow
{
  Window window;
  double rawWindow = 0.0;
  double rawNominal = 0.0;
  double sigmaProxyPb = 0.0;
  double perEntryPb = 0.0;
  std::unique_ptr<TH1> hist;
};

struct BuiltScenario
{
  Scenario scenario;
  std::vector<BuiltWindow> components;
  std::unique_ptr<TH1> total;
  std::unique_ptr<TF1> fit;
  std::unique_ptr<TH1> ratio;
  bool ok = false;
  double logRms = std::numeric_limits<double>::infinity();
  double maxAbs = std::numeric_limits<double>::infinity();
};

std::string SafeKey(std::string s)
{
  for (char& c : s)
  {
    if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-')) c = '_';
  }
  return s;
}

std::string Fixed(double v, int p = 3)
{
  std::ostringstream os;
  os << std::fixed << std::setprecision(p) << v;
  return os.str();
}

std::string Sci(double v, int p = 3)
{
  std::ostringstream os;
  os << std::scientific << std::setprecision(p) << v;
  return os.str();
}

double IntegralWindow(const TH1* h, double lo, double hi)
{
  if (!h) return 0.0;
  double sum = 0.0;
  for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
  {
    const double x = h->GetBinCenter(ib);
    if (x < lo || x >= hi) continue;
    sum += h->GetBinContent(ib);
  }
  return sum;
}

std::unique_ptr<TH1> CloneWindow(const TH1* src, const std::string& name, double lo, double hi)
{
  if (!src) return nullptr;
  std::unique_ptr<TH1> out(dynamic_cast<TH1*>(src->Clone(name.c_str())));
  if (!out) return nullptr;
  out->SetDirectory(nullptr);
  if (out->GetSumw2N() == 0) out->Sumw2();
  for (int ib = 1; ib <= out->GetNbinsX(); ++ib)
  {
    const double x = out->GetBinCenter(ib);
    if (x < lo || x >= hi)
    {
      out->SetBinContent(ib, 0.0);
      out->SetBinError(ib, 0.0);
    }
  }
  return out;
}

void StyleHist(TH1* h, int color, int marker, double markerSize = 0.75)
{
  if (!h) return;
  h->SetStats(false);
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  h->SetMarkerSize(markerSize);
  h->SetLineWidth(2);
}

std::unique_ptr<TH1> LoadWindowHist(const Slice& slice,
                                   const Window& window,
                                   const std::string& name,
                                   double& rawWindow,
                                   double& rawNominal,
                                   double& sigmaProxyPb,
                                   double& perEntryPb)
{
  rawWindow = 0.0;
  rawNominal = 0.0;
  sigmaProxyPb = 0.0;
  perEntryPb = 0.0;

  std::unique_ptr<TFile> f(TFile::Open(slice.path.c_str(), "READ"));
  if (!f || f->IsZombie())
  {
    std::cerr << "[ERROR] Cannot open " << slice.path << "\n";
    return nullptr;
  }
  TDirectory* d = f->GetDirectory("SIM");
  if (!d)
  {
    std::cerr << "[ERROR] Missing SIM directory in " << slice.path << "\n";
    return nullptr;
  }
  TH1* hAll = dynamic_cast<TH1*>(d->Get(slice.histName.c_str()));
  if (!hAll)
  {
    std::cerr << "[ERROR] Missing SIM/" << slice.histName << " in " << slice.path << "\n";
    return nullptr;
  }

  rawWindow = IntegralWindow(hAll, window.lo, window.hi);
  rawNominal = IntegralWindow(hAll, slice.nominalLo, slice.nominalHi);
  if (!(rawWindow > 0.0) || !(rawNominal > 0.0) || !(slice.sigmaNominalPb > 0.0))
  {
    std::cerr << "[WARN] Bad counts for " << slice.key << " window " << window.lo
              << "-" << window.hi << ": rawWindow=" << rawWindow
              << " rawNominal=" << rawNominal << "\n";
    return nullptr;
  }

  sigmaProxyPb = slice.sigmaNominalPb * rawWindow / rawNominal;
  perEntryPb = sigmaProxyPb / rawWindow;
  auto h = CloneWindow(hAll, name, window.lo, window.hi);
  if (!h) return nullptr;
  h->Scale(perEntryPb);
  h->Rebin(2);
  StyleHist(h.get(), slice.color, slice.marker, 0.78);
  h->SetTitle("");
  return h;
}

double PositiveNear(const TH1* h, double x0)
{
  if (!h) return 1.0;
  const int ib0 = h->GetXaxis()->FindBin(x0);
  for (int delta = 0; delta <= 6; ++delta)
  {
    for (int sign : {-1, 1})
    {
      const int ib = ib0 + sign * delta;
      if (ib < 1 || ib > h->GetNbinsX()) continue;
      const double y = h->GetBinContent(ib);
      if (y > 0.0) return y;
    }
  }
  return 1.0;
}

double EstimateSlope(const TH1* h, double xmin, double xmax)
{
  double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
  for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
  {
    const double x = h->GetBinCenter(ib);
    const double y = h->GetBinContent(ib);
    if (x < xmin || x > xmax || y <= 0.0) continue;
    if (x1 <= 0.0)
    {
      x1 = x;
      y1 = y;
    }
    x2 = x;
    y2 = y;
  }
  if (!(x1 > 0.0) || !(x2 > x1) || !(y1 > 0.0) || !(y2 > 0.0)) return 6.0;
  return std::clamp(-std::log(y2 / y1) / std::log(x2 / x1), 0.5, 40.0);
}

std::unique_ptr<TF1> FitModifiedPowerLaw(const TH1* h,
                                        const std::string& name,
                                        double xmin,
                                        double xmax,
                                        int& fitStatus)
{
  fitStatus = -1;
  if (!h) return nullptr;
  TGraphErrors g;
  int ip = 0;
  for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
  {
    const double x = h->GetBinCenter(ib);
    const double y = h->GetBinContent(ib);
    if (x < xmin || x > xmax || y <= 0.0) continue;
    g.SetPoint(ip, x, std::log(y));
    g.SetPointError(ip, 0.0, 1.0);
    ++ip;
  }
  if (g.GetN() < 8) return nullptr;

  TF1 logFit((name + "_logfit").c_str(),
             "[0] + ([1] + [2]*TMath::Log(x) + [3]*x)*TMath::Log(1.0/x)",
             xmin, xmax);
  const double y20 = std::max(PositiveNear(h, 20.0), 1.0e-30);
  const double slope = EstimateSlope(h, xmin, xmax);
  const double c1 = 2.0;
  const double c2 = 0.01;
  const double logA = std::log(y20) - (slope + c1 * std::log(20.0) + c2 * 20.0) * std::log(1.0 / 20.0);
  logFit.SetParameters(logA, slope, c1, c2);
  g.Fit(&logFit, "QNR");
  TFitResultPtr result = g.Fit(&logFit, "QNR S");
  fitStatus = static_cast<int>(result);

  auto fit = std::make_unique<TF1>((name + "_modified_power_law").c_str(),
                                   "TMath::Exp([0]) * TMath::Power(1.0/x, [1] + [2]*TMath::Log(x) + [3]*x)",
                                   xmin, xmax);
  for (int i = 0; i < 4; ++i) fit->SetParameter(i, logFit.GetParameter(i));
  fit->SetLineColor(kGray + 2);
  fit->SetLineStyle(2);
  fit->SetLineWidth(2);
  return fit;
}

std::unique_ptr<TH1> RatioToFit(const TH1* h, const TF1* fit, const std::string& name)
{
  if (!h || !fit) return nullptr;
  std::unique_ptr<TH1> r(dynamic_cast<TH1*>(h->Clone(name.c_str())));
  if (!r) return nullptr;
  r->SetDirectory(nullptr);
  r->Reset("ICES");
  if (r->GetSumw2N() == 0) r->Sumw2();
  for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
  {
    const double x = h->GetBinCenter(ib);
    const double y = h->GetBinContent(ib);
    const double den = fit->Eval(x);
    if (y <= 0.0 || den <= 0.0) continue;
    r->SetBinContent(ib, y / den);
    r->SetBinError(ib, h->GetBinError(ib) / den);
  }
  StyleHist(r.get(), kBlack, 20, 0.70);
  return r;
}

double MeanRatio(const TH1* ratio, double lo, double hi)
{
  if (!ratio) return std::numeric_limits<double>::quiet_NaN();
  double sum = 0.0;
  double wsum = 0.0;
  for (int ib = 1; ib <= ratio->GetNbinsX(); ++ib)
  {
    const double x = ratio->GetBinCenter(ib);
    const double y = ratio->GetBinContent(ib);
    if (x < lo || x >= hi || y <= 0.0) continue;
    const double err = ratio->GetBinError(ib);
    const double w = (err > 0.0) ? 1.0 / (err * err) : 1.0;
    sum += w * y;
    wsum += w;
  }
  return (wsum > 0.0) ? sum / wsum : std::numeric_limits<double>::quiet_NaN();
}

BuiltScenario BuildScenario(const std::vector<Slice>& slices,
                            const Scenario& scenario,
                            const std::string& prefix,
                            double fitMin,
                            double fitMax)
{
  BuiltScenario built;
  built.scenario = scenario;
  for (std::size_t iw = 0; iw < scenario.windows.size(); ++iw)
  {
    const Window& w = scenario.windows[iw];
    if (w.slice < 0 || w.slice >= static_cast<int>(slices.size())) continue;
    BuiltWindow bw;
    bw.window = w;
    bw.hist = LoadWindowHist(slices[w.slice], w,
                             prefix + "_" + scenario.key + "_" + std::to_string(iw),
                             bw.rawWindow, bw.rawNominal, bw.sigmaProxyPb, bw.perEntryPb);
    if (!bw.hist) continue;
    built.components.push_back(std::move(bw));
  }
  if (built.components.size() != scenario.windows.size())
  {
    std::cerr << "[INVALID] " << prefix << " " << scenario.key
              << ": built " << built.components.size() << "/"
              << scenario.windows.size()
              << " components; skipping incomplete scenario." << std::endl;
    built.components.clear();
    return built;
  }

  built.total.reset(dynamic_cast<TH1*>(built.components.front().hist->Clone((prefix + "_" + scenario.key + "_total").c_str())));
  built.total->SetDirectory(nullptr);
  built.total->Reset("ICES");
  if (built.total->GetSumw2N() == 0) built.total->Sumw2();
  for (const auto& comp : built.components) built.total->Add(comp.hist.get());
  StyleHist(built.total.get(), kBlack, 20, 0.72);

  int fitStatus = -1;
  built.fit = FitModifiedPowerLaw(built.total.get(), prefix + "_" + scenario.key, fitMin, fitMax, fitStatus);
  built.ratio = RatioToFit(built.total.get(), built.fit.get(), prefix + "_" + scenario.key + "_ratio");
  built.ok = static_cast<bool>(built.fit) && static_cast<bool>(built.ratio);
  if (!built.ok) return built;

  double sumLog2 = 0.0;
  int n = 0;
  double maxAbs = 0.0;
  for (int ib = 1; ib <= built.ratio->GetNbinsX(); ++ib)
  {
    const double x = built.ratio->GetBinCenter(ib);
    const double r = built.ratio->GetBinContent(ib);
    if (x < fitMin || x > fitMax || r <= 0.0) continue;
    sumLog2 += std::pow(std::log(r), 2);
    maxAbs = std::max(maxAbs, std::abs(r - 1.0));
    ++n;
  }
  built.logRms = (n > 0) ? std::sqrt(sumLog2 / n) : std::numeric_limits<double>::infinity();
  built.maxAbs = maxAbs;
  return built;
}

void DrawFullScenario(const std::vector<Slice>& slices,
                      const BuiltScenario& built,
                      const std::string& sampleTitle,
                      const std::string& xTitle,
                      const std::string& outputPrefix,
                      double xMin,
                      double xMax)
{
  if (!built.ok) return;

  const std::string outPng = kOutDir + "/" + outputPrefix + "_" + SafeKey(built.scenario.key) + ".png";
  TCanvas c(("c_" + outputPrefix + "_" + built.scenario.key).c_str(), "", 1120, 860);
  TPad top("top", "top", 0.0, 0.31, 1.0, 1.0);
  TPad bot("bot", "bot", 0.0, 0.0, 1.0, 0.32);
  top.SetLeftMargin(0.115);
  top.SetRightMargin(0.035);
  top.SetBottomMargin(0.025);
  top.SetLogy();
  bot.SetLeftMargin(0.115);
  bot.SetRightMargin(0.035);
  bot.SetTopMargin(0.03);
  bot.SetBottomMargin(0.27);
  top.Draw();
  bot.Draw();

  top.cd();
  built.total->SetTitle("");
  built.total->GetXaxis()->SetRangeUser(xMin, xMax);
  built.total->GetYaxis()->SetTitle("proxy #sigma_{eff}-scaled entries [pb / bin]");
  built.total->GetYaxis()->SetTitleOffset(1.05);
  built.total->SetMinimum(std::max(1.0e-12, built.total->GetMaximum() * 2.0e-5));
  built.total->SetMaximum(std::max(1.0e-9, built.total->GetMaximum() * 90.0));
  built.total->GetXaxis()->SetLabelSize(0.0);
  built.total->GetXaxis()->SetTitleSize(0.0);
  built.total->Draw("AXIS");
  for (const auto& comp : built.components) comp.hist->Draw("E1 SAME");
  built.fit->Draw("SAME");
  built.total->Draw("E1 SAME");
  gPad->RedrawAxis();

  TLegend leg(0.15, 0.13, 0.49, 0.13 + 0.055 * (built.components.size() + 2));
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.032);
  for (const auto& comp : built.components)
  {
    leg.AddEntry(comp.hist.get(), comp.window.label.c_str(), "lep");
  }
  leg.AddEntry(built.total.get(), "combined sum", "lep");
  leg.AddEntry(built.fit.get(), "Fit: modified power law", "l");
  leg.Draw();

  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(0.039);
  lat.DrawLatex(0.15, 0.88, sampleTitle.c_str());
  lat.SetTextSize(0.031);
  lat.DrawLatex(0.15, 0.825, built.scenario.label.c_str());
  lat.SetTextSize(0.027);
  lat.DrawLatex(0.15, 0.780, built.scenario.productionExact
                             ? "Nominal window: uses existing generator xsec constants"
                             : "Shifted-window diagnostic: xsecs are local count-fraction proxies; rerun generator xsec before production");
  TLatex sph;
  sph.SetNDC(true);
  sph.SetTextFont(42);
  sph.SetTextAlign(33);
  sph.SetTextSize(0.040);
  sph.DrawLatex(0.92, 0.88, "#it{#bf{sPHENIX}} Internal");
  sph.SetTextSize(0.032);
  sph.DrawLatex(0.92, 0.835, "Pythia Overlay #sqrt{s_{NN}} = 200 GeV");

  for (double b : built.scenario.boundaries)
  {
    TLine line(b, built.total->GetMinimum(), b, built.total->GetMaximum() / 4.0);
    line.SetLineColor(kGray + 1);
    line.SetLineStyle(3);
    line.Draw("SAME");
  }

  bot.cd();
  built.ratio->SetTitle("");
  built.ratio->GetXaxis()->SetRangeUser(xMin, xMax);
  built.ratio->GetXaxis()->SetTitle(xTitle.c_str());
  built.ratio->GetYaxis()->SetTitle("combined / fit");
  built.ratio->GetYaxis()->SetRangeUser(0.85, 1.15);
  built.ratio->GetYaxis()->SetNdivisions(505);
  built.ratio->GetYaxis()->SetTitleSize(0.085);
  built.ratio->GetYaxis()->SetLabelSize(0.075);
  built.ratio->GetYaxis()->SetTitleOffset(0.48);
  built.ratio->GetXaxis()->SetTitleSize(0.085);
  built.ratio->GetXaxis()->SetLabelSize(0.075);
  built.ratio->Draw("E1");
  TLine one(xMin, 1.0, xMax, 1.0);
  one.SetLineColor(kGray + 1);
  one.SetLineStyle(2);
  one.Draw("SAME");
  for (double b : built.scenario.boundaries)
  {
    TLine line(b, 0.85, b, 1.15);
    line.SetLineColor(kGray + 1);
    line.SetLineStyle(3);
    line.Draw("SAME");
  }
  c.SaveAs(outPng.c_str());
  std::cout << "[DONE] Wrote " << outPng << "\n";
}

void DrawRatioComparison(const std::vector<BuiltScenario>& built,
                         const std::string& outputPrefix,
                         const std::string& title,
                         const std::string& xTitle,
                         double xMin,
                         double xMax)
{
  const std::string outPng = kOutDir + "/" + outputPrefix + "_ratio_comparison.png";
  int nOk = 0;
  for (const auto& b : built)
  {
    if (b.ok) ++nOk;
  }
  if (nOk == 0)
  {
    std::cerr << "[INVALID] " << outputPrefix
              << ": no complete scenarios; writing invalid-placeholder PNG." << std::endl;
    TCanvas cInvalid(("c_" + outputPrefix + "_invalid").c_str(), "", 1050, 640);
    cInvalid.SetMargin(0.08, 0.05, 0.12, 0.10);
    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextFont(42);
    lat.SetTextSize(0.045);
    lat.DrawLatex(0.12, 0.65, "INVALID diagnostic input");
    lat.SetTextSize(0.032);
    lat.DrawLatex(0.12, 0.56, "At least one required slice histogram is missing.");
    lat.DrawLatex(0.12, 0.49, "Do not use this output as a stitched-spectrum comparison.");
    cInvalid.SaveAs(outPng.c_str());
    return;
  }

  TCanvas c(("c_" + outputPrefix + "_ratio_comparison").c_str(), "", 1050, 640);
  c.SetLeftMargin(0.105);
  c.SetRightMargin(0.04);
  c.SetBottomMargin(0.13);
  c.SetTopMargin(0.10);

  TH1F frame(("frame_" + outputPrefix).c_str(), "", 10, xMin, xMax);
  frame.SetStats(false);
  frame.GetYaxis()->SetRangeUser(0.85, 1.15);
  frame.GetXaxis()->SetTitle(xTitle.c_str());
  frame.GetYaxis()->SetTitle("combined / own modified-power-law fit");
  frame.GetYaxis()->SetTitleOffset(0.85);
  frame.Draw("AXIS");

  const int colors[] = {kBlack, kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 2, kOrange + 7};
  TLegend leg(0.15, 0.69, 0.50, 0.90);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.031);
  int idx = 0;
  for (const auto& b : built)
  {
    if (!b.ok) continue;
    b.ratio->SetLineColor(colors[idx % 6]);
    b.ratio->SetMarkerColor(colors[idx % 6]);
    b.ratio->SetMarkerStyle(20 + idx);
    b.ratio->SetMarkerSize(0.62);
    b.ratio->Draw("E1 SAME");
    leg.AddEntry(b.ratio.get(), b.scenario.shortLabel.c_str(), "lep");
    ++idx;
  }
  TLine one(xMin, 1.0, xMax, 1.0);
  one.SetLineColor(kGray + 1);
  one.SetLineStyle(2);
  one.Draw("SAME");
  leg.Draw();

  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(0.039);
  lat.DrawLatex(0.15, 0.935, title.c_str());
  lat.SetTextSize(0.028);
  lat.DrawLatex(0.15, 0.640, "Each curve uses the same fit family; shifted-window normalizations are diagnostic proxies.");

  c.SaveAs(outPng.c_str());
  std::cout << "[DONE] Wrote " << outPng << "\n";
}

void WriteMetrics(const std::vector<BuiltScenario>& built,
                  const std::string& outputPrefix)
{
  const std::string csv = kOutDir + "/" + outputPrefix + "_metrics.csv";
  std::ofstream out(csv);
  out << "scenario,production_exact,log_rms,max_abs,boundary,left_mean,right_mean,right_over_left,component,window_lo,window_hi,raw_window,raw_nominal,sigma_proxy_pb,per_entry_pb\n";
  for (const auto& b : built)
  {
    if (!b.ok) continue;
    for (double boundary : b.scenario.boundaries)
    {
      const double left = MeanRatio(b.ratio.get(), boundary - 2.0, boundary);
      const double right = MeanRatio(b.ratio.get(), boundary, boundary + 2.0);
      const double jump = (std::isfinite(left) && left != 0.0) ? right / left : std::numeric_limits<double>::quiet_NaN();
      for (const auto& comp : b.components)
      {
        out << b.scenario.key << ","
            << (b.scenario.productionExact ? 1 : 0) << ","
            << b.logRms << ","
            << b.maxAbs << ","
            << boundary << ","
            << left << ","
            << right << ","
            << jump << ","
            << comp.window.label << ","
            << comp.window.lo << ","
            << comp.window.hi << ","
            << comp.rawWindow << ","
            << comp.rawNominal << ","
            << comp.sigmaProxyPb << ","
            << comp.perEntryPb << "\n";
      }
    }
  }
  std::cout << "[DONE] Wrote " << csv << "\n";
}

void WriteReadme()
{
  const std::string path = kOutDir + "/README.md";
  std::ofstream out(path);
  out << "# Embedded Stitch Boundary Correctness Scan\n\n";
  out << "Generated by `macros/MakeEmbeddedStitchBoundaryCorrectnessScan.C`.\n\n";
  out << "Purpose: test whether visible kinks are driven by hard sample ownership boundaries, "
      << "without pretending the shifted-window normalizations are final production constants.\n\n";
  out << "The nominal scenarios use existing generator effective cross sections. Shifted scenarios "
      << "use local RecoilJets `h_*_all` count fractions to proxy the changed effective xsec; "
      << "those must be replaced by generator-only `estimateEmbeddedPhotonXsec.sh` reruns before "
      << "production use.\n\n";
  out << "Inputs are per-slice local ROOT files, not final physics claims from tainted merged "
      << "Jet12+20+30 outputs.\n";
  std::cout << "[DONE] Wrote " << path << "\n";
}
}

void MakeEmbeddedStitchBoundaryCorrectnessScan()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0.5);
  gSystem->mkdir(kOutDir.c_str(), true);

  const std::vector<Slice> photonSlices = {
      {"photon12", "PhotonJet12",
       "InputFiles/simEmbedded/RecoilJets_embeddedPhoton12_ALL_" + kCfg + ".root",
       "h_embedStitch_filterPhotonPt_all", 12.0, 20.0, 2598.12425, kBlue + 1, 20},
      {"photon20", "PhotonJet20",
       "InputFiles/simEmbedded/RecoilJets_embeddedPhoton20_ALL_" + kCfg + ".root",
       "h_embedStitch_filterPhotonPt_all", 20.0, 1.0e9, 133.317866, kRed + 1, 21},
  };
  const std::vector<Scenario> photonScenarios = {
      {"boundary20", "Nominal ownership: PhotonJet12 12-20, PhotonJet20 >=20", "boundary 20 GeV",
       {{0, "PhotonJet12: 12-20 GeV", 12.0, 20.0}, {1, "PhotonJet20: >=20 GeV", 20.0, 1.0e9}},
       {20.0}, true},
      {"boundary21", "Diagnostic ownership: PhotonJet12 12-21, PhotonJet20 >=21", "boundary 21 GeV",
       {{0, "PhotonJet12: 12-21 GeV", 12.0, 21.0}, {1, "PhotonJet20: >=21 GeV", 21.0, 1.0e9}},
       {21.0}, false},
      {"boundary22", "PPG12-inspired diagnostic: PhotonJet12 12-22, PhotonJet20 >=22", "boundary 22 GeV",
       {{0, "PhotonJet12: 12-22 GeV", 12.0, 22.0}, {1, "PhotonJet20: >=22 GeV", 22.0, 1.0e9}},
       {22.0}, false},
  };

  std::vector<BuiltScenario> builtPhoton;
  for (const auto& s : photonScenarios)
  {
    auto b = BuildScenario(photonSlices, s, "embedded_photon", 12.0, 45.0);
    DrawFullScenario(photonSlices, b, "Embedded PhotonJet12+20 generator stitching",
                     "p_{T,filter}^{#gamma} [GeV]", "embedded_photon12_20", 10.0, 45.0);
    builtPhoton.push_back(std::move(b));
  }
  DrawRatioComparison(builtPhoton, "embedded_photon12_20",
                      "PhotonJet12+20 boundary scan", "p_{T,filter}^{#gamma} [GeV]", 10.0, 45.0);
  WriteMetrics(builtPhoton, "embedded_photon12_20");

  const std::string incBase = "InputFiles/InclusiveJetSIM_EMBEDDED/inclusive3_20260515_172346/";
  const std::vector<Slice> jetSlices = {
      {"jet12", "Jet12", incBase + "RecoilJets_embeddedJet12_ALL_" + kCfg + ".root",
       "h_embedInclusiveStitch_filterJetPt_all", 12.0, 20.0, 1.21692467e6, kBlue + 1, 20},
      {"jet20", "Jet20", incBase + "RecoilJets_embeddedJet20_ALL_" + kCfg + ".root",
       "h_embedInclusiveStitch_filterJetPt_all", 20.0, 30.0, 5.44464934e4, kOrange + 7, 21},
      {"jet30", "Jet30", incBase + "RecoilJets_embeddedJet30_ALL_" + kCfg + ".root",
       "h_embedInclusiveStitch_filterJetPt_all", 30.0, 1.0e9, 2.40291630e3, kMagenta + 2, 22},
  };
  const std::vector<Scenario> jetScenarios = {
      {"nominal_20_30", "Nominal ownership: Jet12 12-20, Jet20 20-30, Jet30 >=30", "12-20 / 20-30 / >=30",
       {{0, "Jet12: 12-20 GeV", 12.0, 20.0}, {1, "Jet20: 20-30 GeV", 20.0, 30.0}, {2, "Jet30: >=30 GeV", 30.0, 1.0e9}},
       {20.0, 30.0}, true},
      {"boundary21", "Blair/minimal diagnostic: Jet12 12-21, Jet20 21-30, Jet30 >=30", "12-21 / 21-30 / >=30",
       {{0, "Jet12: 12-21 GeV", 12.0, 21.0}, {1, "Jet20: 21-30 GeV", 21.0, 30.0}, {2, "Jet30: >=30 GeV", 30.0, 1.0e9}},
       {21.0, 30.0}, false},
      {"boundary22", "Stronger diagnostic: Jet12 12-22, Jet20 22-30, Jet30 >=30", "12-22 / 22-30 / >=30",
       {{0, "Jet12: 12-22 GeV", 12.0, 22.0}, {1, "Jet20: 22-30 GeV", 22.0, 30.0}, {2, "Jet30: >=30 GeV", 30.0, 1.0e9}},
       {22.0, 30.0}, false},
      {"ppg12_style", "PPG12-style diagnostic: Jet12 14-21, Jet20 21-32, Jet30 >=32", "14-21 / 21-32 / >=32",
       {{0, "Jet12: 14-21 GeV", 14.0, 21.0}, {1, "Jet20: 21-32 GeV", 21.0, 32.0}, {2, "Jet30: >=32 GeV", 32.0, 1.0e9}},
       {21.0, 32.0}, false},
  };

  std::vector<BuiltScenario> builtJets;
  for (const auto& s : jetScenarios)
  {
    auto b = BuildScenario(jetSlices, s, "embedded_inclusive_jet", 12.0, 50.0);
    DrawFullScenario(jetSlices, b, "Embedded inclusive Jet12+20+30 generator stitching",
                     "max p_{T}^{jet,truth} [GeV]", "embedded_inclusive_jet12_20_30", 10.0, 50.0);
    builtJets.push_back(std::move(b));
  }
  DrawRatioComparison(builtJets, "embedded_inclusive_jet12_20_30",
                      "Embedded inclusive Jet12+20+30 boundary scan",
                      "max p_{T}^{jet,truth} [GeV]", 10.0, 50.0);
  WriteMetrics(builtJets, "embedded_inclusive_jet12_20_30");
  WriteReadme();
}
