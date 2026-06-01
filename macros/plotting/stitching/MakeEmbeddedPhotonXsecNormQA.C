#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMath.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TH1.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace
{
// IMPORTANT:
//   After RecoilJets_AuAu stitching, PhotonJet12 must use the EXCLUSIVE
//   generator cross section for 12 <= pT_filter^gamma < 20 GeV.
constexpr double kSigmaEmbeddedPhoton12To20_pb = 2598.12425;
constexpr double kSigmaEmbeddedPhoton20Plus_pb = 133.317866;
constexpr double kSigmaEmbeddedInclusiveJet12_pb = 1.21692467e6;
constexpr double kSigmaEmbeddedInclusiveJet20_pb = 5.56198698e4;

std::string gConfigTag =
    "jetMinPt5_7pi_8_vz60_isoR40_fixedIso2GeV_baseVariant_"
    "preselectionReference_tightReference_nonTightReference";
std::string gOutDirOverride;

std::string ConfigTag()
{
  return gConfigTag;
}

std::string File12()
{
  return "InputFiles/simEmbedded/RecoilJets_embeddedPhoton12_ALL_" + ConfigTag() + ".root";
}

std::string File20()
{
  return "InputFiles/simEmbedded/RecoilJets_embeddedPhoton20_ALL_" + ConfigTag() + ".root";
}

std::string OutDir()
{
  if (!gOutDirOverride.empty()) return gOutDirOverride;
  return "dataOutput/combinedSimOnlyEMBEDDED/" + ConfigTag() + "/photonJet12and20merged_SIM";
}

std::string InclusiveFile12()
{
  return "InputFiles/InclusiveJetSIM_EMBEDDED/RecoilJets_embeddedJet12_ALL_" + ConfigTag() + ".root";
}

std::string InclusiveFile20()
{
  return "InputFiles/InclusiveJetSIM_EMBEDDED/RecoilJets_embeddedJet20_ALL_" + ConfigTag() + ".root";
}

std::string InclusiveOutDir()
{
  return "dataOutput/combinedSimOnlyEMBEDDED/" + ConfigTag() + "/embeddedJet12and20merged_SIM";
}

int VzCutFromConfigTag()
{
  const std::string cfg = ConfigTag();
  if (cfg.find("vz60") != std::string::npos) return 60;
  if (cfg.find("vz30") != std::string::npos) return 30;
  std::cerr << "[WARN] Could not infer vz cut from config tag: " << cfg
            << ". Falling back to 30 cm." << std::endl;
  return 30;
}

double IsoConeRFromConfigTag()
{
  const std::string cfg = ConfigTag();
  if (cfg.find("isoR40") != std::string::npos) return 0.4;
  if (cfg.find("isoR30") != std::string::npos) return 0.3;
  std::cerr << "[WARN] Could not infer isolation cone radius from config tag: " << cfg
            << ". Falling back to 0.3." << std::endl;
  return 0.3;
}

std::string VzEtaLabel()
{
  return TString::Format("|v_{z}| < %d cm,  |#eta^{#gamma}| < 0.7", VzCutFromConfigTag()).Data();
}

std::string IsoConeLabel()
{
  return TString::Format("#DeltaR_{cone} < %.1f", IsoConeRFromConfigTag()).Data();
}

double ReadEventCount(TFile* f)
{
  if (!f) return 0.0;
  TDirectory* d = f->GetDirectory("SIM");
  if (!d) return 0.0;

  TH1* cnt = dynamic_cast<TH1*>(d->Get("cnt_SIM"));
  if (!cnt) return 0.0;

  const double bin1 = cnt->GetBinContent(1);
  if (bin1 > 0.0) return bin1;

  const double integral = cnt->Integral(0, cnt->GetNbinsX() + 1);
  if (integral > 0.0) return integral;

  return cnt->GetEntries();
}

TH1* CloneFromDir(TDirectory* d, const std::string& name, const std::string& cloneName)
{
  if (!d) return nullptr;
  TH1* h = dynamic_cast<TH1*>(d->Get(name.c_str()));
  if (!h) return nullptr;

  TH1* c = dynamic_cast<TH1*>(h->Clone(cloneName.c_str()));
  if (!c) return nullptr;
  c->SetDirectory(nullptr);
  if (c->GetSumw2N() == 0) c->Sumw2();
  return c;
}

std::unique_ptr<TH1> BuildRegionPtSpectrum(TFile* f,
                                           const std::vector<std::string>& regions,
                                           const std::string& cloneName)
{
  if (!f) return nullptr;
  TDirectory* d = f->GetDirectory("SIM");
  if (!d) return nullptr;

  auto addSet = [&](const std::vector<std::string>& suffixes) -> std::unique_ptr<TH1>
  {
    std::unique_ptr<TH1> sum;
    for (const std::string& suffix : suffixes)
    {
      for (const std::string& region : regions)
      {
        const std::string name = "h_pTgamma_ABCD_" + region + suffix;
        std::unique_ptr<TH1> h(CloneFromDir(d, name, cloneName + "_" + region + suffix));
        if (!h) continue;

        if (!sum)
        {
          sum.reset(dynamic_cast<TH1*>(h->Clone(cloneName.c_str())));
          if (!sum) return nullptr;
          sum->SetDirectory(nullptr);
          sum->Reset("ICES");
          if (sum->GetSumw2N() == 0) sum->Sumw2();
        }
        sum->Add(h.get());
      }
    }
    return sum;
  };

  std::unique_ptr<TH1> inclusive = addSet({""});
  if (inclusive && inclusive->Integral(0, inclusive->GetNbinsX() + 1) > 0.0)
  {
    return inclusive;
  }

  return addSet({"_cent_0_20", "_cent_20_50", "_cent_50_80"});
}

std::unique_ptr<TH1> BuildKeptFilterPtSpectrum(TFile* f, const std::string& cloneName)
{
  if (!f) return nullptr;
  TDirectory* d = f->GetDirectory("SIM");
  return std::unique_ptr<TH1>(CloneFromDir(d, "h_embedStitch_filterPhotonPt_kept", cloneName));
}

std::unique_ptr<TH1> BuildKeptInclusiveJetFilterPtSpectrum(TFile* f, const std::string& cloneName)
{
  if (!f) return nullptr;
  TDirectory* d = f->GetDirectory("SIM");
  return std::unique_ptr<TH1>(CloneFromDir(d, "h_embedInclusiveStitch_filterJetPt_kept", cloneName));
}

void StyleHist(TH1* h, int color, int marker)
{
  if (!h) return;
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  h->SetMarkerSize(0.85);
  h->SetLineWidth(2);
}

std::string Sci(double v, int p = 3)
{
  std::ostringstream os;
  os << std::scientific << std::setprecision(p) << v;
  return os.str();
}

std::string Fixed(double v, int p = 3)
{
  std::ostringstream os;
  os << std::fixed << std::setprecision(p) << v;
  return os.str();
}

enum class FitFamily
{
  kPowerLaw,
  kHagedorn,
  kModifiedPowerLaw,
  kCurvedPowerLaw,
  kLogPoly4,
  kLogPoly6,
  kLogPoly8,
  kLogPoly10,
  kLogPoly12,
  kLogCubicSpline,
  kPiecewisePowerLaw,
  kPiecewiseLogQuad,
  kExponential
};

struct FitReference
{
  FitFamily family = FitFamily::kPowerLaw;
  std::string tag;
  std::string label;
  std::string formulaLabel;
  std::unique_ptr<TF1> curve;
  std::unique_ptr<TH1> ratio;
  bool ok = false;
  int status = -1;
  int nPoints = 0;
  double chi2Ndf = std::numeric_limits<double>::infinity();
  double logRms = std::numeric_limits<double>::infinity();
  double maxAbs = std::numeric_limits<double>::infinity();
  double boundaryMaxAbs = std::numeric_limits<double>::infinity();
  double score = std::numeric_limits<double>::infinity();
};

std::string FitFamilyTag(FitFamily family)
{
  switch (family)
  {
    case FitFamily::kPowerLaw: return "powerlaw";
    case FitFamily::kHagedorn: return "hagedorn";
    case FitFamily::kModifiedPowerLaw: return "modified_powerlaw";
    case FitFamily::kCurvedPowerLaw: return "curved_powerlaw";
    case FitFamily::kLogPoly4: return "logpoly4";
    case FitFamily::kLogPoly6: return "logpoly6";
    case FitFamily::kLogPoly8: return "logpoly8";
    case FitFamily::kLogPoly10: return "logpoly10";
    case FitFamily::kLogPoly12: return "logpoly12";
    case FitFamily::kLogCubicSpline: return "log_cubic_spline";
    case FitFamily::kPiecewisePowerLaw: return "piecewise_powerlaw";
    case FitFamily::kPiecewiseLogQuad: return "piecewise_logquad";
    case FitFamily::kExponential: return "exponential";
  }
  return "unknown";
}

std::string FitFamilyLabel(FitFamily family)
{
  switch (family)
  {
    case FitFamily::kPowerLaw: return "Fit reference";
    case FitFamily::kHagedorn: return "Fit reference";
    case FitFamily::kModifiedPowerLaw: return "Modified power-law fit";
    case FitFamily::kCurvedPowerLaw: return "Fit reference";
    case FitFamily::kLogPoly4: return "Fit reference";
    case FitFamily::kLogPoly6: return "Fit reference";
    case FitFamily::kLogPoly8: return "Fit reference";
    case FitFamily::kLogPoly10: return "Fit reference";
    case FitFamily::kLogPoly12: return "Fit reference";
    case FitFamily::kLogCubicSpline: return "Fit reference";
    case FitFamily::kPiecewisePowerLaw: return "Fit reference";
    case FitFamily::kPiecewiseLogQuad: return "Fit reference";
    case FitFamily::kExponential: return "Fit reference";
  }
  return "Fit reference";
}

std::string FitFamilyFormulaLabel(FitFamily family)
{
  switch (family)
  {
    case FitFamily::kPowerLaw: return "Fit: A (p_{T}/20)^{-n}";
    case FitFamily::kHagedorn: return "Fit: A (1+p_{T}/p_{0})^{-n}";
    case FitFamily::kModifiedPowerLaw: return "Fit: modified power law";
    case FitFamily::kCurvedPowerLaw: return "Fit: curved power law";
    case FitFamily::kLogPoly4: return "Fit: log-polynomial continuum";
    case FitFamily::kLogPoly6: return "Fit: log-polynomial continuum";
    case FitFamily::kLogPoly8: return "Fit: log-polynomial continuum";
    case FitFamily::kLogPoly10: return "Fit: log-polynomial continuum";
    case FitFamily::kLogPoly12: return "Fit: log-polynomial continuum";
    case FitFamily::kLogCubicSpline: return "Fit: smooth log-cubic spline";
    case FitFamily::kPiecewisePowerLaw: return "Fit: local power law by stitch region";
    case FitFamily::kPiecewiseLogQuad: return "Fit: piecewise log-quadratic continuum";
    case FitFamily::kExponential: return "Fit: A e^{-bp_{T}}";
  }
  return "Fit reference";
}

int LogPolyOrder(FitFamily family)
{
  switch (family)
  {
    case FitFamily::kLogPoly4: return 4;
    case FitFamily::kLogPoly6: return 6;
    case FitFamily::kLogPoly8: return 8;
    case FitFamily::kLogPoly10: return 10;
    case FitFamily::kLogPoly12: return 12;
    default: return -1;
  }
}

std::string LogPolyExpression(int order, int offset = 0)
{
  std::ostringstream formula;
  formula << "[" << offset << "] + [" << offset + 1 << "]*TMath::Log(x/20.0)";
  for (int ip = 2; ip <= order; ++ip)
  {
    formula << " + [" << offset + ip << "]*TMath::Power(TMath::Log(x/20.0), " << ip << ")";
  }
  return formula.str();
}

std::string PiecewiseLogPolyExpression(int order, const std::vector<double>& boundaries)
{
  const int nSegments = static_cast<int>(boundaries.size()) + 1;
  auto segmentExpr = [&](int segment) {
    return LogPolyExpression(order, segment * (order + 1));
  };

  std::string expr = segmentExpr(nSegments - 1);
  for (int segment = nSegments - 2; segment >= 0; --segment)
  {
    std::ostringstream wrapped;
    wrapped << "(x < " << boundaries[segment] << " ? (" << segmentExpr(segment)
            << ") : (" << expr << "))";
    expr = wrapped.str();
  }
  return expr;
}

std::string LogCubicSplineExpression(const std::vector<double>& boundaries)
{
  std::ostringstream formula;
  formula << "[0] + [1]*TMath::Log(x/20.0)"
          << " + [2]*TMath::Power(TMath::Log(x/20.0), 2)"
          << " + [3]*TMath::Power(TMath::Log(x/20.0), 3)";
  for (std::size_t ib = 0; ib < boundaries.size(); ++ib)
  {
    const double boundary = boundaries[ib];
    const double knot = std::log(boundary / 20.0);
    formula << " + [" << (4 + static_cast<int>(ib)) << "]"
            << "*(x > " << std::setprecision(16) << boundary
            << " ? TMath::Power(TMath::Log(x/20.0) - " << knot << ", 3) : 0.0)";
  }
  return formula.str();
}

double PositiveContentNear(const TH1* h, double x)
{
  if (!h) return 1.0;
  const int ib = h->GetXaxis()->FindBin(x);
  const int lo = std::max(1, ib - 2);
  const int hi = std::min(h->GetNbinsX(), ib + 2);
  for (int j = 0; j <= 2; ++j)
  {
    const int left = ib - j;
    const int right = ib + j;
    if (left >= lo && left <= hi && h->GetBinContent(left) > 0.0) return h->GetBinContent(left);
    if (right >= lo && right <= hi && h->GetBinContent(right) > 0.0) return h->GetBinContent(right);
  }
  for (int i = 1; i <= h->GetNbinsX(); ++i)
  {
    if (h->GetBinContent(i) > 0.0) return h->GetBinContent(i);
  }
  return 1.0;
}

double EstimatePowerSlope(const TH1* h, double fitXMin, double fitXMax)
{
  double xFirst = 0.0;
  double yFirst = 0.0;
  double xLast = 0.0;
  double yLast = 0.0;
  for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
  {
    const double x = h->GetBinCenter(ib);
    const double y = h->GetBinContent(ib);
    if (x < fitXMin || x > fitXMax || y <= 0.0) continue;
    if (xFirst <= 0.0)
    {
      xFirst = x;
      yFirst = y;
    }
    xLast = x;
    yLast = y;
  }
  if (xFirst <= 0.0 || xLast <= xFirst || yFirst <= 0.0 || yLast <= 0.0) return 6.0;
  const double slope = -std::log(yLast / yFirst) / std::log(xLast / xFirst);
  return std::min(20.0, std::max(0.5, slope));
}

std::unique_ptr<TGraphErrors> BuildLogSpectrumGraph(const TH1* h,
                                                    const std::string& name,
                                                    double fitXMin,
                                                    double fitXMax)
{
  auto graph = std::make_unique<TGraphErrors>();
  graph->SetName(name.c_str());
  int ip = 0;
  for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
  {
    const double x = h->GetBinCenter(ib);
    const double y = h->GetBinContent(ib);
    if (x < fitXMin || x > fitXMax || y <= 0.0) continue;
    graph->SetPoint(ip, x, std::log(y));
    graph->SetPointError(ip, 0.0, 1.0);
    ++ip;
  }
  return graph;
}

std::unique_ptr<TH1> MakeRatioToFit(const TH1* h, const TF1* fit, const std::string& name)
{
  if (!h || !fit) return nullptr;
  std::unique_ptr<TH1> ratio(dynamic_cast<TH1*>(h->Clone(name.c_str())));
  if (!ratio) return nullptr;
  ratio->SetDirectory(nullptr);
  ratio->Reset("ICES");
  if (ratio->GetSumw2N() == 0) ratio->Sumw2();
  for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
  {
    const double x = h->GetBinCenter(ib);
    const double den = fit->Eval(x);
    const double num = h->GetBinContent(ib);
    if (den <= 0.0 || num <= 0.0) continue;
    ratio->SetBinContent(ib, num / den);
    ratio->SetBinError(ib, h->GetBinError(ib) / den);
  }
  return ratio;
}

FitReference BuildFitReference(const TH1* h,
                               FitFamily family,
                               const std::string& name,
                               double drawXMin,
                               double drawXMax,
                               double fitXMin,
                               double fitXMax,
                               const std::vector<double>& stitchBoundaries)
{
  FitReference ref;
  ref.family = family;
  ref.tag = FitFamilyTag(family);
  ref.label = FitFamilyLabel(family);
  ref.formulaLabel = FitFamilyFormulaLabel(family);
  if (!h) return ref;

  auto graph = BuildLogSpectrumGraph(h, name + "_logGraph_" + ref.tag, fitXMin, fitXMax);
  ref.nPoints = graph ? graph->GetN() : 0;
  const bool isPiecewisePowerLaw = (family == FitFamily::kPiecewisePowerLaw);
  const bool isPiecewiseLogQuad = (family == FitFamily::kPiecewiseLogQuad);
  const bool isModifiedPowerLaw = (family == FitFamily::kModifiedPowerLaw);
  const bool isPiecewise = (isPiecewisePowerLaw || isPiecewiseLogQuad);
  const bool isLogCubicSpline = (family == FitFamily::kLogCubicSpline);
  const int logPolyOrder = isPiecewisePowerLaw ? 1 :
                           isPiecewiseLogQuad ? 2 : LogPolyOrder(family);
  const int nPar = isLogCubicSpline ? (4 + static_cast<int>(stitchBoundaries.size())) :
                   isPiecewise ? ((static_cast<int>(stitchBoundaries.size()) + 1) * (logPolyOrder + 1)) :
                   (logPolyOrder > 0) ? (logPolyOrder + 1) :
                   isModifiedPowerLaw ? 4 :
                   (family == FitFamily::kHagedorn || family == FitFamily::kCurvedPowerLaw) ? 3 : 2;
  if (!graph || ref.nPoints < nPar + 2) return ref;

  std::string logFormula;
  if (family == FitFamily::kPowerLaw) logFormula = "[0] - [1]*TMath::Log(x/20.0)";
  if (family == FitFamily::kHagedorn) logFormula = "[0] - [2]*TMath::Log(1.0 + x/[1])";
  if (isModifiedPowerLaw) logFormula = "[0] + ([1] + [2]*TMath::Log(x) + [3]*x)*TMath::Log(1.0/x)";
  if (family == FitFamily::kCurvedPowerLaw) logFormula = "[0] + [1]*TMath::Log(x/20.0) + [2]*TMath::Power(TMath::Log(x/20.0), 2)";
  if (isLogCubicSpline) logFormula = LogCubicSplineExpression(stitchBoundaries);
  if (isPiecewise) logFormula = PiecewiseLogPolyExpression(logPolyOrder, stitchBoundaries);
  if (!isPiecewise && logPolyOrder > 0) logFormula = LogPolyExpression(logPolyOrder);
  if (family == FitFamily::kExponential) logFormula = "[0] - [1]*(x-20.0)";

  auto logFit = std::make_unique<TF1>((name + "_logFit_" + ref.tag).c_str(), logFormula.c_str(), fitXMin, fitXMax);
  const double y20 = std::max(PositiveContentNear(h, 20.0), 1.0e-30);
  const double slope = EstimatePowerSlope(h, fitXMin, fitXMax);
  if (family == FitFamily::kPowerLaw)
  {
    logFit->SetParameters(std::log(y20), slope);
    logFit->SetParLimits(1, 0.05, 40.0);
  }
  else if (family == FitFamily::kHagedorn)
  {
    const double p0 = 12.0;
    const double n = std::max(2.0, slope + 2.0);
    logFit->SetParameters(std::log(y20) + n * std::log(1.0 + 20.0 / p0), p0, n);
    logFit->SetParLimits(1, 1.0, 120.0);
    logFit->SetParLimits(2, 0.05, 80.0);
  }
  else if (isModifiedPowerLaw)
  {
    const double n = slope;
    const double c1 = 2.0;
    const double c2 = 0.01;
    const double logA = std::log(y20) -
                        (n + c1 * std::log(20.0) + c2 * 20.0) *
                        std::log(1.0 / 20.0);
    logFit->SetParameters(logA, n, c1, c2);
    logFit->SetParNames("logA", "n", "c1", "c2");
  }
  else if (family == FitFamily::kCurvedPowerLaw)
  {
    logFit->SetParameters(std::log(y20), -slope, 0.0);
    logFit->SetParLimits(1, -50.0, 5.0);
    logFit->SetParLimits(2, -50.0, 50.0);
  }
  else if (isLogCubicSpline)
  {
    logFit->SetParameter(0, std::log(y20));
    logFit->SetParameter(1, -slope);
    for (int ip = 2; ip < nPar; ++ip) logFit->SetParameter(ip, 0.0);
  }
  else if (isPiecewise)
  {
    std::vector<double> edges;
    edges.push_back(fitXMin);
    for (double boundary : stitchBoundaries) edges.push_back(boundary);
    edges.push_back(fitXMax);
    for (int segment = 0; segment + 1 < static_cast<int>(edges.size()); ++segment)
    {
      const double center = 0.5 * (edges[segment] + edges[segment + 1]);
      const double yCenter = std::max(PositiveContentNear(h, center), 1.0e-30);
      const double zCenter = std::log(center / 20.0);
      const int offset = segment * (logPolyOrder + 1);
      logFit->SetParameter(offset, std::log(yCenter) + slope * zCenter);
      logFit->SetParameter(offset + 1, -slope);
      for (int ip = 2; ip <= logPolyOrder; ++ip) logFit->SetParameter(offset + ip, 0.0);
    }
  }
  else if (logPolyOrder > 0)
  {
    logFit->SetParameter(0, std::log(y20));
    logFit->SetParameter(1, -slope);
    for (int ip = 2; ip < nPar; ++ip) logFit->SetParameter(ip, 0.0);
  }
  else
  {
    const double b = std::max(0.001, slope / 20.0);
    logFit->SetParameters(std::log(y20), b);
    logFit->SetParLimits(1, 0.0001, 5.0);
  }

  if (isModifiedPowerLaw) graph->Fit(logFit.get(), "QNR");
  TFitResultPtr result = graph->Fit(logFit.get(), "QNR S");
  ref.status = static_cast<int>(result);
  if (ref.status != 0)
  {
    std::cout << "[FIT] " << name << " " << ref.tag
              << " status=" << ref.status << " entering ratio scoring with penalty" << std::endl;
  }

  std::string curveFormula;
  if (family == FitFamily::kPowerLaw) curveFormula = "TMath::Exp([0]) * TMath::Power(x/20.0, -[1])";
  if (family == FitFamily::kHagedorn) curveFormula = "TMath::Exp([0]) * TMath::Power(1.0 + x/[1], -[2])";
  if (isModifiedPowerLaw) curveFormula = "TMath::Exp([0]) * TMath::Power(1.0/x, [1] + [2]*TMath::Log(x) + [3]*x)";
  if (family == FitFamily::kCurvedPowerLaw) curveFormula = "TMath::Exp([0] + [1]*TMath::Log(x/20.0) + [2]*TMath::Power(TMath::Log(x/20.0), 2))";
  if (isLogCubicSpline) curveFormula = "TMath::Exp(" + LogCubicSplineExpression(stitchBoundaries) + ")";
  if (isPiecewise) curveFormula = "TMath::Exp(" + PiecewiseLogPolyExpression(logPolyOrder, stitchBoundaries) + ")";
  if (!isPiecewise && logPolyOrder > 0) curveFormula = "TMath::Exp(" + LogPolyExpression(logPolyOrder) + ")";
  if (family == FitFamily::kExponential) curveFormula = "TMath::Exp([0] - [1]*(x-20.0))";
  ref.curve = std::make_unique<TF1>((name + "_curve_" + ref.tag).c_str(), curveFormula.c_str(), fitXMin, drawXMax);
  for (int ip = 0; ip < nPar; ++ip) ref.curve->SetParameter(ip, logFit->GetParameter(ip));
  ref.curve->SetLineColor(kGray + 2);
  ref.curve->SetLineStyle(2);
  ref.curve->SetLineWidth(2);

  ref.ratio = MakeRatioToFit(h, ref.curve.get(), name + "_ratio_" + ref.tag);
  if (!ref.ratio) return ref;
  ref.ratio->SetStats(false);
  StyleHist(ref.ratio.get(), kBlack, 20);

  double sumLog2 = 0.0;
  double maxAbs = 0.0;
  double boundaryMaxAbs = 0.0;
  int nUsed = 0;
  for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
  {
    const double x = h->GetBinCenter(ib);
    const double y = h->GetBinContent(ib);
    const double den = ref.curve->Eval(x);
    if (x < fitXMin || x > fitXMax || y <= 0.0 || den <= 0.0) continue;
    const double r = y / den;
    if (!std::isfinite(r) || r <= 0.0) continue;
    const double abs = std::abs(r - 1.0);
    sumLog2 += std::pow(std::log(r), 2);
    maxAbs = std::max(maxAbs, abs);
    for (double boundary : stitchBoundaries)
    {
      if (std::abs(x - boundary) <= 2.0) boundaryMaxAbs = std::max(boundaryMaxAbs, abs);
    }
    ++nUsed;
  }
  if (nUsed <= nPar) return ref;

  ref.chi2Ndf = logFit->GetNDF() > 0 ? logFit->GetChisquare() / logFit->GetNDF() : std::numeric_limits<double>::infinity();
  ref.logRms = std::sqrt(sumLog2 / nUsed);
  ref.maxAbs = maxAbs;
  ref.boundaryMaxAbs = boundaryMaxAbs;
  const double complexityPenalty = (family == FitFamily::kPowerLaw) ? 0.0 :
                                   (family == FitFamily::kHagedorn) ? 0.004 :
                                   (family == FitFamily::kModifiedPowerLaw) ? 0.006 :
                                   (family == FitFamily::kCurvedPowerLaw) ? 0.005 :
                                   (family == FitFamily::kLogPoly4) ? 0.006 :
                                   (family == FitFamily::kLogPoly6) ? 0.007 :
                                   (family == FitFamily::kLogPoly8) ? 0.008 :
                                   (family == FitFamily::kLogPoly10) ? 0.010 :
                                   (family == FitFamily::kLogPoly12) ? 0.012 :
                                   (family == FitFamily::kLogCubicSpline) ? 0.010 :
                                   (family == FitFamily::kPiecewisePowerLaw) ? 0.010 :
                                   (family == FitFamily::kPiecewiseLogQuad) ? 0.016 : 0.020;
  const double statusPenalty = (ref.status == 0) ? 0.0 : 0.010;
  ref.score = ref.logRms + 0.35 * ref.boundaryMaxAbs + 0.10 * ref.maxAbs + complexityPenalty + statusPenalty;
  ref.ok = std::isfinite(ref.score);
  std::cout << "[FIT] " << name << " " << ref.tag
            << " status=" << ref.status
            << " n=" << nUsed
            << " chi2ndf=" << ref.chi2Ndf
            << " logRms=" << ref.logRms
            << " maxAbs=" << ref.maxAbs
            << " boundaryMaxAbs=" << ref.boundaryMaxAbs
            << " score=" << ref.score << std::endl;
  return ref;
}

int SelectFitReference(const std::vector<FitReference>& refs)
{
  int best = -1;
  for (std::size_t i = 0; i < refs.size(); ++i)
  {
    if (!refs[i].ok) continue;
    if (best < 0 || refs[i].score < refs[best].score) best = static_cast<int>(i);
  }
  return best;
}

int SelectFitReference(const std::vector<FitReference>& refs, FitFamily requiredFamily)
{
  for (std::size_t i = 0; i < refs.size(); ++i)
  {
    if (refs[i].ok && refs[i].family == requiredFamily) return static_cast<int>(i);
  }
  return -1;
}

std::unique_ptr<TH1> WeightedClone(std::unique_ptr<TH1> h, double weight)
{
  if (!h) return nullptr;
  h->Scale(weight);
  return h;
}

std::unique_ptr<TH1> BuildWeightedPairHistogram(TFile* f12,
                                                TFile* f20,
                                                const std::string& histName,
                                                const std::string& cloneName,
                                                double w12,
                                                double w20)
{
  TDirectory* d12 = f12 ? f12->GetDirectory("SIM") : nullptr;
  TDirectory* d20 = f20 ? f20->GetDirectory("SIM") : nullptr;

  std::unique_ptr<TH1> h12(CloneFromDir(d12, histName, cloneName + "_12"));
  std::unique_ptr<TH1> h20(CloneFromDir(d20, histName, cloneName + "_20"));
  if (!h12 && !h20) return nullptr;

  std::unique_ptr<TH1> hSum;
  if (h12)
  {
    hSum.reset(dynamic_cast<TH1*>(h12->Clone(cloneName.c_str())));
    if (!hSum) return nullptr;
    hSum->SetDirectory(nullptr);
    hSum->Scale(w12);
  }
  else if (h20)
  {
    hSum.reset(dynamic_cast<TH1*>(h20->Clone(cloneName.c_str())));
    if (!hSum) return nullptr;
    hSum->SetDirectory(nullptr);
    hSum->Reset("ICES");
    if (hSum->GetSumw2N() == 0) hSum->Sumw2();
  }

  if (h20)
  {
    h20->Scale(w20);
    hSum->Add(h20.get());
  }

  return hSum;
}

std::unique_ptr<TH1> BuildWeightedPairHistogramSum(TFile* f12,
                                                   TFile* f20,
                                                   const std::vector<std::string>& histNames,
                                                   const std::string& cloneName,
                                                   double w12,
                                                   double w20)
{
  std::unique_ptr<TH1> hSum;
  int nAdded = 0;
  for (const std::string& histName : histNames)
  {
    std::unique_ptr<TH1> h = BuildWeightedPairHistogram(
        f12, f20, histName, cloneName + "_" + std::to_string(nAdded), w12, w20);
    if (!h) continue;

    if (!hSum)
    {
      hSum.reset(dynamic_cast<TH1*>(h->Clone(cloneName.c_str())));
      if (!hSum) return nullptr;
      hSum->SetDirectory(nullptr);
      hSum->Reset("ICES");
      if (hSum->GetSumw2N() == 0) hSum->Sumw2();
    }

    hSum->Add(h.get());
    ++nAdded;
  }

  return hSum;
}

bool FindEfficiencyCut(TH1* hIn, double eff, double& cut, double& cutErr)
{
  cut = 0.0;
  cutErr = 0.0;
  if (!hIn) return false;

  const int nb = hIn->GetNbinsX();
  const double total = hIn->Integral(1, nb);
  if (!(total > 0.0)) return false;

  double running = 0.0;
  int ibCut = nb;
  for (int ib = 1; ib <= nb; ++ib)
  {
    running += hIn->GetBinContent(ib);
    if ((running / total) >= eff)
    {
      ibCut = ib;
      break;
    }
  }

  const double binLo = hIn->GetXaxis()->GetBinLowEdge(ibCut);
  const double binHi = hIn->GetXaxis()->GetBinUpEdge(ibCut);
  const double prev = running - hIn->GetBinContent(ibCut);
  const double binC = hIn->GetBinContent(ibCut);
  const double target = eff * total;

  if (binC > 0.0)
  {
    const double frac = std::min(1.0, std::max(0.0, (target - prev) / binC));
    cut = binLo + frac * (binHi - binLo);
  }
  else
  {
    cut = 0.5 * (binLo + binHi);
  }

  cutErr = 0.5 * (binHi - binLo);
  return std::isfinite(cut);
}

void StyleEffGraph(TGraphErrors& g, int color, int marker)
{
  g.SetLineWidth(2);
  g.SetLineColor(color);
  g.SetMarkerColor(color);
  g.SetMarkerStyle(marker);
  g.SetMarkerSize(marker == 22 ? 1.45 : 1.35);
}

struct IsoFlatCentResult
{
  int lo = 0;
  int hi = 0;
  double center = 0.0;
  double halfWidth = 0.0;
  bool have70 = false;
  bool have80 = false;
  bool have90 = false;
  double flat70 = 0.0;
  double err70 = 0.0;
  double flat80 = 0.0;
  double err80 = 0.0;
  double flat90 = 0.0;
  double err90 = 0.0;
};

bool ComputeFlatGraphAverage(TGraphErrors& g, double xLo, double xHi, double& value, double& error)
{
  value = 0.0;
  error = 0.0;

  double sumW = 0.0;
  double sumWY = 0.0;
  int nUsed = 0;

  for (int ip = 0; ip < g.GetN(); ++ip)
  {
    double x = 0.0;
    double y = 0.0;
    g.GetPoint(ip, x, y);
    if (x < xLo || x > xHi) continue;

    const double ey = g.GetErrorY(ip);
    const double w = (ey > 0.0) ? (1.0 / (ey * ey)) : 1.0;
    sumW += w;
    sumWY += w * y;
    ++nUsed;
  }

  if (nUsed <= 0 || !(sumW > 0.0)) return false;
  value = sumWY / sumW;
  error = std::sqrt(1.0 / sumW);
  return std::isfinite(value);
}

std::vector<IsoFlatCentResult> DrawIsoEfficiencyCutoffQA(TFile* f12, TFile* f20, double w12, double w20)
{
  const std::vector<double> ptEdges = {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 35};
  const double fitXLo = 10.0;
  const double fitXHi = 35.0;
  struct CentBin
  {
    int lo;
    int hi;
    std::vector<std::string> suffixes;
    std::string tag;
  };
  const std::vector<CentBin> centBins = {
      {0, 20, {"_cent_0_20"}, "cent_0_20"},
      {20, 50, {"_cent_20_50"}, "cent_20_50"},
      {50, 80, {"_cent_50_80"}, "cent_50_80"},
  };

  std::vector<IsoFlatCentResult> flatResults;

  for (const auto& cb : centBins)
  {
    std::vector<double> x;
    std::vector<double> ex;
    std::vector<double> y70, ey70;
    std::vector<double> y80, ey80;
    std::vector<double> y90, ey90;
    double yMin = std::numeric_limits<double>::max();
    double yMax = -std::numeric_limits<double>::max();

    for (std::size_t ipt = 1; ipt + 1 < ptEdges.size(); ++ipt)
    {
      const int lo = static_cast<int>(std::llround(ptEdges[ipt]));
      const int hi = static_cast<int>(std::llround(ptEdges[ipt + 1]));
      std::vector<std::string> histNames;
      histNames.reserve(cb.suffixes.size());
      for (const std::string& suffix : cb.suffixes)
      {
        histNames.push_back("h_EisoReco_truthSigMatched_pT_" +
                            std::to_string(lo) + "_" + std::to_string(hi) + suffix);
      }

      std::unique_ptr<TH1> h = BuildWeightedPairHistogramSum(
          f12, f20, histNames, "h_weightedIso_" + cb.tag + "_" + std::to_string(lo) + "_" + std::to_string(hi), w12, w20);
      if (!h || h->Integral(1, h->GetNbinsX()) <= 0.0) continue;

      double c70 = 0.0, e70 = 0.0;
      double c80 = 0.0, e80 = 0.0;
      double c90 = 0.0, e90 = 0.0;
      if (!FindEfficiencyCut(h.get(), 0.70, c70, e70)) continue;
      if (!FindEfficiencyCut(h.get(), 0.80, c80, e80)) continue;
      if (!FindEfficiencyCut(h.get(), 0.90, c90, e90)) continue;

      x.push_back(0.5 * (ptEdges[ipt] + ptEdges[ipt + 1]));
      ex.push_back(0.0);
      y70.push_back(c70);
      ey70.push_back(e70);
      y80.push_back(c80);
      ey80.push_back(e80);
      y90.push_back(c90);
      ey90.push_back(e90);

      yMin = std::min(yMin, std::min(c70 - e70, std::min(c80 - e80, c90 - e90)));
      yMax = std::max(yMax, std::max(c70 + e70, std::max(c80 + e80, c90 + e90)));
    }

    if (x.empty())
    {
      std::cerr << "[WARN] No iso-efficiency points for " << cb.tag << std::endl;
      continue;
    }

    const double pad = (yMax > yMin) ? 0.70 * (yMax - yMin) : 0.5;
    const double yLo = std::max(0.0, yMin - pad);
    const double yHi = yMax + pad;

    TCanvas c(("c_embeddedPhoton_stitchedIsoCutEfficiency_" + cb.tag).c_str(),
              "stitched embedded photon iso-cut efficiency", 1100, 800);
    c.SetLeftMargin(0.14);
    c.SetRightMargin(0.04);
    c.SetBottomMargin(0.13);
    c.SetTopMargin(0.10);

    TH1F frame(("hFrame_embeddedPhoton_stitchedIsoCutEfficiency_" + cb.tag).c_str(),
               "", 100, 10.0, 35.0);
    frame.SetDirectory(nullptr);
    frame.SetStats(0);
    frame.SetMinimum(yLo);
    frame.SetMaximum(yHi);
    frame.GetXaxis()->SetTitle("Cluster p_{T} [GeV]");
    frame.GetYaxis()->SetTitle("E_{T}^{iso} Cutoff [GeV]");
    frame.GetXaxis()->SetTitleSize(0.060);
    frame.GetYaxis()->SetTitleSize(0.060);
    frame.GetXaxis()->SetLabelSize(0.050);
    frame.GetYaxis()->SetLabelSize(0.050);
    frame.GetYaxis()->SetTitleOffset(1.05);
    frame.Draw();

    TGraphErrors g90(static_cast<int>(x.size()), x.data(), y90.data(), ex.data(), ey90.data());
    TGraphErrors g80(static_cast<int>(x.size()), x.data(), y80.data(), ex.data(), ey80.data());
    TGraphErrors g70(static_cast<int>(x.size()), x.data(), y70.data(), ex.data(), ey70.data());
    StyleEffGraph(g90, kMagenta + 1, 20);
    StyleEffGraph(g80, kGreen + 2, 21);
    StyleEffGraph(g70, kBlue + 1, 22);

    IsoFlatCentResult flat;
    flat.lo = cb.lo;
    flat.hi = cb.hi;
    flat.center = 0.5 * (cb.lo + cb.hi);
    flat.halfWidth = 0.5 * (cb.hi - cb.lo);
    flat.have70 = ComputeFlatGraphAverage(g70, fitXLo, fitXHi, flat.flat70, flat.err70);
    flat.have80 = ComputeFlatGraphAverage(g80, fitXLo, fitXHi, flat.flat80, flat.err80);
    flat.have90 = ComputeFlatGraphAverage(g90, fitXLo, fitXHi, flat.flat90, flat.err90);
    flatResults.push_back(flat);

    TF1 f90(("f_embeddedIsoFlat90_" + cb.tag).c_str(), "[0]", fitXLo, fitXHi);
    TF1 f80(("f_embeddedIsoFlat80_" + cb.tag).c_str(), "[0]", fitXLo, fitXHi);
    TF1 f70(("f_embeddedIsoFlat70_" + cb.tag).c_str(), "[0]", fitXLo, fitXHi);
    if (flat.have90) f90.SetParameter(0, flat.flat90);
    if (flat.have80) f80.SetParameter(0, flat.flat80);
    if (flat.have70) f70.SetParameter(0, flat.flat70);
    f90.SetLineColor(kMagenta + 1);
    f80.SetLineColor(kGreen + 2);
    f70.SetLineColor(kBlue + 1);
    f90.SetLineWidth(3);
    f80.SetLineWidth(3);
    f70.SetLineWidth(3);
    f90.SetLineStyle(2);
    f80.SetLineStyle(2);
    f70.SetLineStyle(2);

    g90.Draw("PE1 SAME");
    g80.Draw("PE1 SAME");
    g70.Draw("PE1 SAME");
    if (flat.have90) f90.Draw("SAME");
    if (flat.have80) f80.Draw("SAME");
    if (flat.have70) f70.Draw("SAME");

    TLatex cutLabel;
    cutLabel.SetTextFont(42);
    cutLabel.SetTextAlign(31);
    cutLabel.SetTextSize(0.030);
    const double labelX = 34.45;
    const double labelYOffset = 0.015 * (yHi - yLo);
    auto drawCutValue = [&](bool have, double value, int color)
    {
      if (!have) return;
      cutLabel.SetTextColor(color);
      cutLabel.DrawLatex(labelX, value + labelYOffset, TString::Format("%.2f", value).Data());
    };
    drawCutValue(flat.have90, flat.flat90, kMagenta + 1);
    drawCutValue(flat.have80, flat.flat80, kGreen + 2);
    drawCutValue(flat.have70, flat.flat70, kBlue + 1);
    cutLabel.SetTextColor(kBlack);

    TLegend leg(0.20, 0.16, 0.78, 0.24);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextFont(42);
    leg.SetTextSize(0.030);
    leg.AddEntry(&g70, "70% Efficiency", "ep");
    leg.AddEntry(&g80, "80% Efficiency", "ep");
    leg.AddEntry(&g90, "90% Efficiency", "ep");
    leg.SetNColumns(3);
    leg.Draw();

    TLatex info;
    info.SetNDC(true);
    info.SetTextFont(42);
    info.SetTextAlign(13);
    info.SetTextSize(0.030);
    info.DrawLatex(0.18, 0.88, "Photon+Jet 12+20 Embedded SIM");
    info.DrawLatex(0.18, 0.83, TString::Format("%d-%d%% centrality", cb.lo, cb.hi).Data());
    info.DrawLatex(0.18, 0.78, VzEtaLabel().c_str());
    info.DrawLatex(0.18, 0.73, IsoConeLabel().c_str());

    TLatex sph;
    sph.SetNDC(true);
    sph.SetTextFont(42);
    sph.SetTextAlign(33);
    sph.SetTextSize(0.042);
    sph.DrawLatex(0.92, 0.88, "#bf{sPHENIX} #it{Internal}");
    sph.SetTextSize(0.034);
    sph.DrawLatex(0.92, 0.83, "Pythia Overlay  #sqrt{s_{NN}} = 200 GeV");

    const std::string outPng = OutDir() + "/embeddedPhoton_stitchedIsoCutEfficiency_" + cb.tag + ".png";
    c.SaveAs(outPng.c_str());
    std::cout << "[DONE] Wrote " << outPng << std::endl;
  }

  return flatResults;
}

void DrawIsoFlatCutoffVsCentralityQA(const std::vector<IsoFlatCentResult>& flatResults)
{
  if (flatResults.empty())
  {
    std::cerr << "[WARN] No flat iso cutoff values available for centrality summary" << std::endl;
    return;
  }

  std::vector<double> x70, y70, ex70, ey70;
  std::vector<double> x80, y80, ex80, ey80;
  std::vector<double> x90, y90, ex90, ey90;

  double yMin = std::numeric_limits<double>::max();
  double yMax = -std::numeric_limits<double>::max();
  auto addPoint = [&](std::vector<double>& x,
                      std::vector<double>& y,
                      std::vector<double>& ex,
                      std::vector<double>& ey,
                      const IsoFlatCentResult& r,
                      double value,
                      double err)
  {
    x.push_back(r.center);
    y.push_back(value);
    ex.push_back(0.0);
    ey.push_back(err);
    yMin = std::min(yMin, value - err);
    yMax = std::max(yMax, value + err);
  };

  for (const auto& r : flatResults)
  {
    if (r.have70) addPoint(x70, y70, ex70, ey70, r, r.flat70, r.err70);
    if (r.have80) addPoint(x80, y80, ex80, ey80, r, r.flat80, r.err80);
    if (r.have90) addPoint(x90, y90, ex90, ey90, r, r.flat90, r.err90);
  }

  if (x70.empty() && x80.empty() && x90.empty()) return;

  const double pad = (yMax > yMin) ? 0.65 * (yMax - yMin) : 0.75;
  TCanvas c("c_embeddedPhoton_stitchedIsoCutEfficiencyFits_vsCentrality",
            "stitched embedded photon flat iso-cut fits vs centrality", 900, 700);
  c.SetLeftMargin(0.14);
  c.SetRightMargin(0.04);
  c.SetBottomMargin(0.13);
  c.SetTopMargin(0.10);

  TH1F frame("hFrame_embeddedPhoton_stitchedIsoCutEfficiencyFits_vsCentrality",
             "", 100, 0.0, 80.0);
  frame.SetDirectory(nullptr);
  frame.SetStats(0);
  frame.SetMinimum(std::max(0.0, yMin - pad));
  frame.SetMaximum(yMax + pad);
  frame.GetXaxis()->SetTitle("Centrality [%]");
  frame.GetYaxis()->SetTitle("E_{T}^{iso} Cutoff [GeV]");
  frame.GetXaxis()->SetTitleSize(0.055);
  frame.GetYaxis()->SetTitleSize(0.055);
  frame.GetXaxis()->SetLabelSize(0.045);
  frame.GetYaxis()->SetLabelSize(0.045);
  frame.GetYaxis()->SetTitleOffset(1.15);
  frame.Draw();

  TLegend leg(0.48, 0.62, 0.92, 0.78);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.028);

  auto drawFitGraph = [&](std::vector<double>& x,
                          std::vector<double>& y,
                          std::vector<double>& ex,
                          std::vector<double>& ey,
                          const char* name,
                          const char* label,
                          int color,
                          int marker)
  {
    if (x.empty()) return;
    TGraphErrors* g = new TGraphErrors(static_cast<int>(x.size()), x.data(), y.data(), ex.data(), ey.data());
    g->SetLineWidth(2);
    g->SetLineColor(color);
    g->SetMarkerColor(color);
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(marker == 22 ? 1.5 : 1.2);
    g->Draw("P SAME");

    TF1* fit = new TF1((std::string("fit_") + name).c_str(), "pol1", 0.0, 80.0);
    fit->SetLineColor(color);
    fit->SetLineWidth(2);
    fit->SetLineStyle(2);
    g->Fit(fit, "QNR");
    fit->Draw("SAME");
    leg.AddEntry(g, TString::Format("%s: y = %.4fx %+.2f", label, fit->GetParameter(1), fit->GetParameter(0)).Data(), "lp");
  };

  drawFitGraph(x90, y90, ex90, ey90, "90", "90% Eff", kMagenta + 1, 20);
  drawFitGraph(x80, y80, ex80, ey80, "80", "80% Eff", kGreen + 2, 21);
  drawFitGraph(x70, y70, ex70, ey70, "70", "70% Eff", kBlue + 1, 22);
  leg.Draw();

  TLatex title;
  title.SetNDC(true);
  title.SetTextFont(42);
  title.SetTextAlign(23);
  title.SetTextSize(0.042);
  title.DrawLatex(0.50, 0.98, "Flat E_{T}^{iso} cutoff vs centrality");

  TLatex info;
  info.SetNDC(true);
  info.SetTextFont(42);
  info.SetTextAlign(13);
  info.SetTextSize(0.030);
  info.DrawLatex(0.18, 0.88, "Photon+Jet 12+20 Embedded SIM");
  info.DrawLatex(0.18, 0.83, "constant fit over full plotted p_{T}^{#gamma} range");
  info.DrawLatex(0.18, 0.78, VzEtaLabel().c_str());
  info.DrawLatex(0.18, 0.73, IsoConeLabel().c_str());

  TLatex sph;
  sph.SetNDC(true);
  sph.SetTextFont(42);
  sph.SetTextAlign(33);
  sph.SetTextSize(0.042);
  sph.DrawLatex(0.92, 0.88, "#bf{sPHENIX} #it{Internal}");
  sph.SetTextSize(0.034);
  sph.DrawLatex(0.92, 0.83, "Pythia Overlay  #sqrt{s_{NN}} = 200 GeV");

  const std::string outPng = OutDir() + "/embeddedPhoton_stitchedIsoCutEfficiencyFits_vsCentrality.png";
  c.SaveAs(outPng.c_str());
  std::cout << "[DONE] Wrote " << outPng << std::endl;
}

void DrawIsoCentralityOverlay10To12QA(TFile* f12, TFile* f20, double w12, double w20)
{
  struct CentBin
  {
    int lo;
    int hi;
    std::vector<std::string> suffixes;
    std::string tag;
    int color;
  };

  const std::vector<CentBin> centBins = {
      {0, 20, {"_cent_0_20"}, "cent_0_20", kBlack},
      {20, 50, {"_cent_20_50"}, "cent_20_50", kBlue + 1},
      {50, 80, {"_cent_50_80"}, "cent_50_80", kOrange + 1},
  };

  std::vector<std::unique_ptr<TH1>> hOwned;
  std::vector<TH1*> hCents;
  std::vector<std::string> labels;
  double yMax = 0.0;

  for (const auto& cb : centBins)
  {
    std::vector<std::string> histNames;
    histNames.reserve(cb.suffixes.size());
    for (const std::string& suffix : cb.suffixes)
    {
      histNames.push_back("h_Eiso_pT_10_12" + suffix);
    }

    std::unique_ptr<TH1> h = BuildWeightedPairHistogramSum(
        f12, f20, histNames, "h_combinedIsoCentOverlay_10_12_" + cb.tag, w12, w20);
    if (!h || h->Integral(1, h->GetNbinsX()) <= 0.0)
    {
      std::cerr << "[WARN] Missing or empty 10-12 GeV isolation histograms for "
                << cb.tag << " centrality overlay" << std::endl;
      continue;
    }

    h->Rebin(10);
    if (h->GetSumw2N() == 0) h->Sumw2();
    const double integral = h->Integral(1, h->GetNbinsX());
    if (integral <= 0.0) continue;
    h->Scale(1.0 / integral);

    h->SetTitle("");
    h->SetStats(0);
    h->SetLineColor(cb.color);
    h->SetMarkerColor(cb.color);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.9);
    h->SetLineWidth(2);
    h->SetFillStyle(0);

    yMax = std::max(yMax, h->GetMaximum());
    labels.push_back(std::to_string(cb.lo) + "-" + std::to_string(cb.hi) + "%");
    hCents.push_back(h.get());
    hOwned.push_back(std::move(h));
  }

  if (hCents.empty())
  {
    std::cerr << "[WARN] No centrality histograms available for embeddedPhoton_stitchedIsoCentralityOverlay_pT_10_12.png" << std::endl;
    return;
  }

  TCanvas c("c_embeddedPhoton_stitchedIsoCentralityOverlay_pT_10_12",
            "stitched embedded photon isolation centrality overlay", 900, 700);
  c.SetLeftMargin(0.14);
  c.SetRightMargin(0.04);
  c.SetBottomMargin(0.13);
  c.SetTopMargin(0.10);

  hCents[0]->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
  hCents[0]->GetYaxis()->SetTitle("Normalized to Unit Area");
  hCents[0]->GetXaxis()->SetTitleSize(0.048);
  hCents[0]->GetYaxis()->SetTitleSize(0.050);
  hCents[0]->GetXaxis()->SetLabelSize(0.040);
  hCents[0]->GetYaxis()->SetLabelSize(0.050);
  hCents[0]->GetYaxis()->SetTitleOffset(1.15);
  hCents[0]->GetXaxis()->SetRangeUser(-10.0, 50.0);
  hCents[0]->SetMinimum(0.0);
  hCents[0]->SetMaximum((yMax > 0.0) ? 1.25 * yMax : 1.0);
  hCents[0]->Draw("E1");
  for (std::size_t ih = 1; ih < hCents.size(); ++ih)
  {
    hCents[ih]->Draw("E1 SAME");
  }

  TLegend leg(0.56, 0.39, 0.89, 0.61);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.038);
  leg.SetNColumns(2);
  for (std::size_t ih = 0; ih < hCents.size(); ++ih)
  {
    leg.AddEntry(hCents[ih], labels[ih].c_str(), "ep");
  }
  leg.Draw();

  TLatex title;
  title.SetNDC(true);
  title.SetTextFont(42);
  title.SetTextAlign(23);
  title.SetTextSize(0.052);
  title.DrawLatex(0.50, 0.99, "Embedded photon centrality overlays, p_{T}^{#gamma} = 10-12 GeV");

  TLatex info;
  info.SetNDC(true);
  info.SetTextFont(42);
  info.SetTextAlign(33);
  info.SetTextSize(0.040);
  info.DrawLatex(0.90, 0.88, "PhotonJet12+20 stitched SIM");
  info.DrawLatex(0.90, 0.84, "with UE Sub");
  info.DrawLatex(0.90, 0.80, IsoConeLabel().c_str());

  TLatex sph;
  sph.SetNDC(true);
  sph.SetTextFont(42);
  sph.SetTextAlign(33);
  sph.SetTextSize(0.052);
  sph.DrawLatex(0.90, 0.32, "#bf{sPHENIX} #it{Internal}");
  sph.SetTextSize(0.042);
  sph.DrawLatex(0.90, 0.26, "Au+Au  #sqrt{s_{NN}} = 200 GeV");

  const std::string outPng = OutDir() + "/embeddedPhoton_stitchedIsoCentralityOverlay_pT_10_12.png";
  c.SaveAs(outPng.c_str());
  std::cout << "[DONE] Wrote " << outPng << std::endl;
}

void DrawSpectrumSmoothQA(std::unique_ptr<TH1> h12,
                          std::unique_ptr<TH1> h20,
                          const std::string& outputName,
                          const std::string& title,
                          const std::string& xTitle,
                          const std::string& note,
                          double w12,
                          double w20)
{
  if (!h12 || !h20)
  {
    std::cerr << "[WARN] Missing histogram(s) for " << outputName << std::endl;
    return;
  }

  const bool isPhotonFilterPtQA = (outputName == "embeddedPhoton_stitchedTruthFilterPtSpectrum");
  const bool isInclusiveJetFilterPtQA = (outputName == "embeddedInclusiveJet_stitchedTruthFilterPtSpectrum");
  const bool isFilterPtQA = (isPhotonFilterPtQA || isInclusiveJetFilterPtQA);
  const double xMin = isFilterPtQA ? 10.0 : h12->GetXaxis()->GetXmin();
  const double xMax = isFilterPtQA ? 45.0 : h12->GetXaxis()->GetXmax();
  const double ratioMin = isFilterPtQA ? 0.85 : 0.45;
  const double ratioMax = isFilterPtQA ? 1.15 : 1.55;

  h12->SetTitle("");
  h20->SetTitle("");
  if (isFilterPtQA)
  {
    h12->Rebin(2);
    h20->Rebin(2);
  }

  StyleHist(h12.get(), kBlue + 1, 20);
  StyleHist(h20.get(), kRed + 1, 21);

  std::unique_ptr<TH1> hSum(dynamic_cast<TH1*>(h12->Clone((outputName + "_sum").c_str())));
  hSum->SetDirectory(nullptr);
  hSum->Add(h20.get());
  hSum->SetTitle("");
  StyleHist(hSum.get(), kBlack, 24);

  const double fitXMin = isFilterPtQA ? 12.0 : xMin;
  const double fitXMax = xMax;
  const std::vector<double> stitchBoundaries = isFilterPtQA ? std::vector<double>{20.0} : std::vector<double>{};
  std::vector<FitReference> fitRefs;
  fitRefs.push_back(BuildFitReference(hSum.get(), FitFamily::kPowerLaw, outputName, xMin, xMax, fitXMin, fitXMax, stitchBoundaries));
  fitRefs.push_back(BuildFitReference(hSum.get(), FitFamily::kHagedorn, outputName, xMin, xMax, fitXMin, fitXMax, stitchBoundaries));
  fitRefs.push_back(BuildFitReference(hSum.get(), FitFamily::kModifiedPowerLaw, outputName, xMin, xMax, fitXMin, fitXMax, stitchBoundaries));
  fitRefs.push_back(BuildFitReference(hSum.get(), FitFamily::kCurvedPowerLaw, outputName, xMin, xMax, fitXMin, fitXMax, stitchBoundaries));
  fitRefs.push_back(BuildFitReference(hSum.get(), FitFamily::kLogPoly4, outputName, xMin, xMax, fitXMin, fitXMax, stitchBoundaries));
  fitRefs.push_back(BuildFitReference(hSum.get(), FitFamily::kLogPoly6, outputName, xMin, xMax, fitXMin, fitXMax, stitchBoundaries));
  fitRefs.push_back(BuildFitReference(hSum.get(), FitFamily::kLogPoly8, outputName, xMin, xMax, fitXMin, fitXMax, stitchBoundaries));
  fitRefs.push_back(BuildFitReference(hSum.get(), FitFamily::kLogPoly10, outputName, xMin, xMax, fitXMin, fitXMax, stitchBoundaries));
  fitRefs.push_back(BuildFitReference(hSum.get(), FitFamily::kLogPoly12, outputName, xMin, xMax, fitXMin, fitXMax, stitchBoundaries));
  fitRefs.push_back(BuildFitReference(hSum.get(), FitFamily::kLogCubicSpline, outputName, xMin, xMax, fitXMin, fitXMax, stitchBoundaries));
  fitRefs.push_back(BuildFitReference(hSum.get(), FitFamily::kPiecewisePowerLaw, outputName, xMin, xMax, fitXMin, fitXMax, stitchBoundaries));
  fitRefs.push_back(BuildFitReference(hSum.get(), FitFamily::kPiecewiseLogQuad, outputName, xMin, xMax, fitXMin, fitXMax, stitchBoundaries));
  fitRefs.push_back(BuildFitReference(hSum.get(), FitFamily::kExponential, outputName, xMin, xMax, fitXMin, fitXMax, stitchBoundaries));
  int selectedFit = -1;
  if (isFilterPtQA) selectedFit = SelectFitReference(fitRefs, FitFamily::kModifiedPowerLaw);
  if (selectedFit < 0) selectedFit = SelectFitReference(fitRefs);
  if (selectedFit < 0)
  {
    std::cerr << "[WARN] No acceptable fit reference for " << outputName << std::endl;
    return;
  }

  hSum->GetXaxis()->SetTitle(xTitle.c_str());
  hSum->GetXaxis()->SetRangeUser(xMin, xMax);
  hSum->GetYaxis()->SetTitle("#sigma_{eff}/N scaled entries [pb / bin]");
  hSum->GetYaxis()->SetTitleOffset(1.12);

  const double ymax = std::max({hSum->GetMaximum(), h12->GetMaximum(), h20->GetMaximum()});
  h12->GetXaxis()->SetRangeUser(xMin, xMax);
  h20->GetXaxis()->SetRangeUser(xMin, xMax);
  hSum->SetMinimum(std::max(1.0e-8, ymax * 2.0e-5));
  hSum->SetMaximum(std::max(1.0e-6, ymax * (isFilterPtQA ? 85.0 : 18.0)));

  auto drawFitPlot = [&](FitReference& fitRef, const std::string& outBaseName)
  {
    fitRef.ratio->SetTitle("");
    fitRef.ratio->GetXaxis()->SetTitle(xTitle.c_str());
    fitRef.ratio->GetXaxis()->SetRangeUser(xMin, xMax);
    fitRef.ratio->GetYaxis()->SetTitle("stitched / fit");
    fitRef.ratio->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
    fitRef.ratio->GetYaxis()->SetNdivisions(505);
    fitRef.ratio->GetYaxis()->SetTitleSize(0.085);
    fitRef.ratio->GetYaxis()->SetLabelSize(0.075);
    fitRef.ratio->GetYaxis()->SetTitleOffset(0.50);
    fitRef.ratio->GetXaxis()->SetTitleSize(0.090);
    fitRef.ratio->GetXaxis()->SetLabelSize(0.080);

    std::unique_ptr<TH1> r12 = MakeRatioToFit(h12.get(), fitRef.curve.get(), outBaseName + "_ratio12");
    std::unique_ptr<TH1> r20 = MakeRatioToFit(h20.get(), fitRef.curve.get(), outBaseName + "_ratio20");
    if (r12)
    {
      r12->SetStats(false);
      StyleHist(r12.get(), kBlue + 1, 20);
    }
    if (r20)
    {
      r20->SetStats(false);
      StyleHist(r20.get(), kRed + 1, 21);
    }

    TCanvas c(("c_" + outBaseName).c_str(), outBaseName.c_str(), 1050, 850);
    TPad top("top", "top", 0.0, 0.30, 1.0, 1.0);
    TPad bot("bot", "bot", 0.0, 0.0, 1.0, 0.31);
    top.SetBottomMargin(0.025);
    top.SetLeftMargin(0.12);
    top.SetRightMargin(0.04);
    top.SetLogy();
    bot.SetTopMargin(0.03);
    bot.SetBottomMargin(0.28);
    bot.SetLeftMargin(0.12);
    bot.SetRightMargin(0.04);
    top.Draw();
    bot.Draw();

    top.cd();
    hSum->GetXaxis()->SetLabelSize(0.0);
    hSum->GetXaxis()->SetTitleSize(0.0);
    hSum->Draw("AXIS");
    h12->Draw("E1 SAME");
    h20->Draw("E1 SAME");
    fitRef.curve->Draw("SAME");
    gPad->RedrawAxis();

    TLegend leg(isFilterPtQA ? 0.15 : 0.55,
                isFilterPtQA ? 0.16 : 0.60,
                isFilterPtQA ? 0.50 : 0.91,
                isFilterPtQA ? 0.41 : 0.87);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.038);
    leg.AddEntry(h12.get(), isInclusiveJetFilterPtQA ? "Jet12 stitched" : (isFilterPtQA ? "PhotonJet12 stitched" : "weighted PhotonJet12"), "lep");
    leg.AddEntry(h20.get(), isInclusiveJetFilterPtQA ? "Jet20 stitched" : (isFilterPtQA ? "PhotonJet20 stitched" : "weighted PhotonJet20"), "lep");
    leg.AddEntry(fitRef.curve.get(), fitRef.label.c_str(), "l");
    leg.Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.040);
    lat.DrawLatex(0.15, 0.86, title.c_str());
    lat.SetTextSize(0.031);
    lat.DrawLatex(0.15, 0.80, note.c_str());
    if (isPhotonFilterPtQA)
    {
      lat.DrawLatex(0.15, 0.75, "PhotonJet12: 12 #leq p_{T,filter}^{#gamma} < 20 GeV; PhotonJet20: p_{T,filter}^{#gamma} #geq 20 GeV");
      lat.DrawLatex(0.15, 0.70, ("w_{12#rightarrow20}=" + Sci(w12) + " pb/event, w_{20+}=" + Sci(w20) + " pb/event").c_str());
    }
    else if (isInclusiveJetFilterPtQA)
    {
      lat.DrawLatex(0.15, 0.75, "Jet12: 12 #leq max p_{T}^{jet,truth} < 20 GeV; Jet20: max p_{T}^{jet,truth} #geq 20 GeV");
      lat.DrawLatex(0.15, 0.70, ("w_{12}=" + Sci(w12) + " pb/event, w_{20}=" + Sci(w20) + " pb/event").c_str());
    }
    else
    {
      lat.DrawLatex(0.15, 0.75, ("w_{12#rightarrow20}=" + Sci(w12) + " pb/event, w_{20+}=" + Sci(w20) + " pb/event").c_str());
    }
    if (isFilterPtQA)
    {
      TLatex tSphM;
      tSphM.SetNDC(true);
      tSphM.SetTextFont(42);
      tSphM.SetTextAlign(33);
      tSphM.SetTextSize(0.042);
      tSphM.DrawLatex(0.92, 0.58, "#bf{sPHENIX} #it{Internal}");
      tSphM.SetTextSize(0.034);
      tSphM.DrawLatex(0.92, 0.53, "Pythia Overlay #sqrt{s_{NN}} = 200 GeV");
    }

    bot.cd();
    fitRef.ratio->Draw("AXIS");
    if (r12) r12->Draw("E1 SAME");
    if (r20) r20->Draw("E1 SAME");
    TLine one(xMin, 1.0, xMax, 1.0);
    one.SetLineColor(kGray + 1);
    one.SetLineStyle(2);
    one.Draw("SAME");

    const std::string outPng = OutDir() + "/" + outBaseName + ".png";
    c.SaveAs(outPng.c_str());
    std::cout << "[DONE] Wrote " << outPng << std::endl;
  };

  for (auto& fitRef : fitRefs)
  {
    if (!fitRef.ok) continue;
    drawFitPlot(fitRef, outputName + "_fit_" + fitRef.tag);
  }
  drawFitPlot(fitRefs[selectedFit], outputName);
  std::cout << "[FIT] selected " << outputName << " -> " << fitRefs[selectedFit].tag << std::endl;
}

void DrawCompositionQA(std::unique_ptr<TH1> h12,
                       std::unique_ptr<TH1> h20,
                       double w12,
                       double w20)
{
  if (!h12 || !h20)
  {
    std::cerr << "[WARN] Missing ABCD histograms for composition QA" << std::endl;
    return;
  }

  StyleHist(h12.get(), kBlue + 1, 20);
  StyleHist(h20.get(), kRed + 1, 21);

  std::unique_ptr<TH1> hSum(dynamic_cast<TH1*>(h12->Clone("h_composition_sum")));
  hSum->SetDirectory(nullptr);
  hSum->Add(h20.get());
  StyleHist(hSum.get(), kBlack, 24);

  std::unique_ptr<TH1> frac12(dynamic_cast<TH1*>(h12->Clone("h_fractionPhoton12")));
  std::unique_ptr<TH1> frac20(dynamic_cast<TH1*>(h20->Clone("h_fractionPhoton20")));
  frac12->SetDirectory(nullptr);
  frac20->SetDirectory(nullptr);
  frac12->Divide(hSum.get());
  frac20->Divide(hSum.get());
  StyleHist(frac12.get(), kBlue + 1, 20);
  StyleHist(frac20.get(), kRed + 1, 21);

  hSum->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
  hSum->GetYaxis()->SetTitle("#sigma_{eff}/N scaled entries [pb / bin]");
  hSum->GetYaxis()->SetTitleOffset(1.12);

  frac12->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
  frac12->GetYaxis()->SetTitle("sample fraction");
  frac12->GetYaxis()->SetRangeUser(0.0, 1.08);
  frac12->GetYaxis()->SetNdivisions(505);
  frac12->GetYaxis()->SetTitleSize(0.085);
  frac12->GetYaxis()->SetLabelSize(0.075);
  frac12->GetYaxis()->SetTitleOffset(0.50);
  frac12->GetXaxis()->SetTitleSize(0.090);
  frac12->GetXaxis()->SetLabelSize(0.080);

  TCanvas c("c_embeddedPhoton_stitchedSampleComposition_ABCDsum",
            "embedded photon stitched sample composition ABCD sum", 1050, 850);
  TPad top("top", "top", 0.0, 0.30, 1.0, 1.0);
  TPad bot("bot", "bot", 0.0, 0.0, 1.0, 0.31);
  top.SetBottomMargin(0.025);
  top.SetLeftMargin(0.12);
  top.SetRightMargin(0.04);
  top.SetLogy();
  bot.SetTopMargin(0.03);
  bot.SetBottomMargin(0.28);
  bot.SetLeftMargin(0.12);
  bot.SetRightMargin(0.04);
  top.Draw();
  bot.Draw();

  top.cd();
  const double ymax = std::max({hSum->GetMaximum(), h12->GetMaximum(), h20->GetMaximum()});
  hSum->SetMinimum(std::max(1.0e-8, ymax * 2.0e-5));
  hSum->SetMaximum(std::max(1.0e-6, ymax * 18.0));
  hSum->Draw("E1");
  h12->Draw("E1 SAME");
  h20->Draw("E1 SAME");

  TLegend leg(0.55, 0.63, 0.91, 0.87);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.040);
  leg.AddEntry(hSum.get(), "weighted ABCD sum", "lep");
  leg.AddEntry(h12.get(), "weighted PhotonJet12 ABCD", "lep");
  leg.AddEntry(h20.get(), "weighted PhotonJet20 ABCD", "lep");
  leg.Draw();

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.040);
  lat.DrawLatex(0.15, 0.86, "Embedded photon stitched sample composition");
  lat.SetTextSize(0.032);
  lat.DrawLatex(0.15, 0.80, "Reference ID; ABCD-summed reco photon spectrum");
  lat.DrawLatex(0.15, 0.75, ("w_{12#rightarrow20}=" + Sci(w12) + " pb/event, w_{20+}=" + Sci(w20) + " pb/event").c_str());

  bot.cd();
  frac12->Draw("E1");
  frac20->Draw("E1 SAME");
  TLine half(frac12->GetXaxis()->GetXmin(), 0.5, frac12->GetXaxis()->GetXmax(), 0.5);
  half.SetLineColor(kGray + 1);
  half.SetLineStyle(2);
  half.Draw("SAME");

  TLegend leg2(0.58, 0.70, 0.91, 0.93);
  leg2.SetBorderSize(0);
  leg2.SetFillStyle(0);
  leg2.SetTextSize(0.075);
  leg2.AddEntry(frac12.get(), "PhotonJet12 / sum", "lep");
  leg2.AddEntry(frac20.get(), "PhotonJet20 / sum", "lep");
  leg2.Draw();

  const std::string outPng = OutDir() + "/embeddedPhoton_stitchedSampleComposition_ABCDsum.png";
  c.SaveAs(outPng.c_str());
  std::cout << "[DONE] Wrote " << outPng << std::endl;

  const std::string legacyNormPng = OutDir() + "/embeddedPhoton12to20_plusPhoton20_xsecNormalizationQA.png";
  c.SaveAs(legacyNormPng.c_str());
  std::cout << "[DONE] Wrote " << legacyNormPng << std::endl;
}
}

void MakeEmbeddedPhotonXsecNormQA()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0.5);

  const std::string file12 = File12();
  const std::string file20 = File20();
  std::unique_ptr<TFile> f12(TFile::Open(file12.c_str(), "READ"));
  std::unique_ptr<TFile> f20(TFile::Open(file20.c_str(), "READ"));
  if (!f12 || f12->IsZombie())
  {
    std::cerr << "[ERROR] Cannot open " << file12 << std::endl;
    return;
  }
  if (!f20 || f20->IsZombie())
  {
    std::cerr << "[ERROR] Cannot open " << file20 << std::endl;
    return;
  }

  if (kSigmaEmbeddedPhoton12To20_pb <= 0.0 || kSigmaEmbeddedPhoton20Plus_pb <= 0.0)
  {
    std::cerr << "[ERROR] Embedded photon cross sections are not set." << std::endl;
    return;
  }

  const double n12 = ReadEventCount(f12.get());
  const double n20 = ReadEventCount(f20.get());
  if (n12 <= 0.0 || n20 <= 0.0)
  {
    std::cerr << "[ERROR] Bad event counts: N12=" << n12 << " N20=" << n20 << std::endl;
    return;
  }

  const double w12 = kSigmaEmbeddedPhoton12To20_pb / n12;
  const double w20 = kSigmaEmbeddedPhoton20Plus_pb / n20;

  gSystem->mkdir(OutDir().c_str(), true);

  DrawSpectrumSmoothQA(
      WeightedClone(BuildKeptFilterPtSpectrum(f12.get(), "h_filterPt12_weighted"), w12),
      WeightedClone(BuildKeptFilterPtSpectrum(f20.get(), "h_filterPt20_weighted"), w20),
      "embeddedPhoton_stitchedTruthFilterPtSpectrum",
      "Embedded photon generator stitching spectrum",
      "p_{T,filter}^{#gamma} [GeV]",
      "Uses SIM/h_embedStitch_filterPhotonPt_kept",
      w12,
      w20);

  const std::vector<IsoFlatCentResult> flatIsoResults = DrawIsoEfficiencyCutoffQA(f12.get(), f20.get(), w12, w20);
  DrawIsoFlatCutoffVsCentralityQA(flatIsoResults);
  DrawIsoCentralityOverlay10To12QA(f12.get(), f20.get(), w12, w20);

  std::cout << "[INFO] N12=" << n12 << " N20=" << n20
            << " w12=" << w12 << " pb/event"
            << " w20=" << w20 << " pb/event" << std::endl;
}

void MakeEmbeddedPhotonXsecNormQA(const char* configTag)
{
  if (!configTag || std::string(configTag).empty())
  {
    std::cerr << "[ERROR] Empty config tag passed to MakeEmbeddedPhotonXsecNormQA" << std::endl;
    return;
  }
  gConfigTag = configTag;
  MakeEmbeddedPhotonXsecNormQA();
}

void MakeEmbeddedInclusiveJetXsecNormQA()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0.5);

  const std::string file12 = InclusiveFile12();
  const std::string file20 = InclusiveFile20();
  std::unique_ptr<TFile> f12(TFile::Open(file12.c_str(), "READ"));
  std::unique_ptr<TFile> f20(TFile::Open(file20.c_str(), "READ"));
  if (!f12 || f12->IsZombie())
  {
    std::cerr << "[ERROR] Cannot open " << file12 << std::endl;
    return;
  }
  if (!f20 || f20->IsZombie())
  {
    std::cerr << "[ERROR] Cannot open " << file20 << std::endl;
    return;
  }

  const double n12 = ReadEventCount(f12.get());
  const double n20 = ReadEventCount(f20.get());
  if (n12 <= 0.0 || n20 <= 0.0)
  {
    std::cerr << "[ERROR] Bad inclusive-jet event counts: N12=" << n12 << " N20=" << n20 << std::endl;
    return;
  }

  const double w12 = kSigmaEmbeddedInclusiveJet12_pb / n12;
  const double w20 = kSigmaEmbeddedInclusiveJet20_pb / n20;

  const std::string oldOutDir = gOutDirOverride;
  gOutDirOverride = InclusiveOutDir();
  gSystem->mkdir(OutDir().c_str(), true);

  DrawSpectrumSmoothQA(
      WeightedClone(BuildKeptInclusiveJetFilterPtSpectrum(f12.get(), "h_inclusiveFilterJetPt12_weighted"), w12),
      WeightedClone(BuildKeptInclusiveJetFilterPtSpectrum(f20.get(), "h_inclusiveFilterJetPt20_weighted"), w20),
      "embeddedInclusiveJet_stitchedTruthFilterPtSpectrum",
      "Embedded inclusive-jet generator stitching spectrum",
      "max p_{T}^{jet,truth} [GeV]",
      "Uses SIM/h_embedInclusiveStitch_filterJetPt_kept",
      w12,
      w20);

  gOutDirOverride = oldOutDir;

  std::cout << "[INFO] Inclusive N12=" << n12 << " N20=" << n20
            << " w12=" << w12 << " pb/event"
            << " w20=" << w20 << " pb/event" << std::endl;
}

void MakeEmbeddedInclusiveJetXsecNormQA(const char* configTag)
{
  if (!configTag || std::string(configTag).empty())
  {
    std::cerr << "[ERROR] Empty config tag passed to MakeEmbeddedInclusiveJetXsecNormQA" << std::endl;
    return;
  }
  gConfigTag = configTag;
  MakeEmbeddedInclusiveJetXsecNormQA();
}

std::vector<std::string> EmbeddedPhotonConfigTagsFromLocalInputs()
{
  std::set<std::string> tags;
  const std::string dir = "InputFiles/simEmbedded";

  const std::string prefix = "RecoilJets_embeddedPhoton12_ALL_";
  const std::string suffix = ".root";
  const TString listing = gSystem->GetFromPipe(("ls " + dir + "/" + prefix + "*" + suffix + " 2>/dev/null").c_str());
  std::istringstream in(listing.Data());
  std::string path;
  while (std::getline(in, path))
  {
    if (path.empty()) continue;
    const std::size_t slash = path.find_last_of('/');
    const std::string name = (slash == std::string::npos) ? path : path.substr(slash + 1);
    if (name.size() <= prefix.size() + suffix.size()) continue;
    if (name.rfind(prefix, 0) != 0) continue;
    if (name.substr(name.size() - suffix.size()) != suffix) continue;
    const std::string tag = name.substr(prefix.size(), name.size() - prefix.size() - suffix.size());

    const std::string file20 = dir + "/RecoilJets_embeddedPhoton20_ALL_" + tag + ".root";
    if (!gSystem->AccessPathName(file20.c_str()))
    {
      tags.insert(tag);
    }
    else
    {
      std::cerr << "[WARN] Skipping tag without matching embeddedPhoton20 file: " << tag << std::endl;
    }
  }

  return std::vector<std::string>(tags.begin(), tags.end());
}

void MakeEmbeddedPhotonXsecNormQA_AllLocalSimEmbedded()
{
  const std::vector<std::string> tags = EmbeddedPhotonConfigTagsFromLocalInputs();
  if (tags.empty())
  {
    std::cerr << "[ERROR] No complete embeddedPhoton12+20 cfg tags found in InputFiles/simEmbedded" << std::endl;
    return;
  }

  std::cout << "[INFO] Regenerating embedded photon stitching QA for "
            << tags.size() << " cfg tag(s)." << std::endl;
  for (const std::string& tag : tags)
  {
    std::cout << "\n[CFG] " << tag << std::endl;
    MakeEmbeddedPhotonXsecNormQA(tag.c_str());
  }
}
