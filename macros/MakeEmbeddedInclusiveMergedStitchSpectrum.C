#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TMath.h>
#include <TNamed.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

namespace
{
const std::string kDefaultInput =
    "dataOutput/combinedSimOnlyEMBEDDED/"
    "preselectionReference_tightReference_nonTightReference_baseVariant/"
    "embeddedJet12and20merged_SIM/RecoilJets_embeddedJet12plus20_MERGED.root";

std::string ParentDir(const std::string& path)
{
  const std::size_t slash = path.find_last_of('/');
  if (slash == std::string::npos) return ".";
  return path.substr(0, slash);
}

std::string Sci(double value, int precision = 3)
{
  std::ostringstream os;
  os.setf(std::ios::scientific);
  os.precision(precision);
  os << value;
  return os.str();
}

bool ExtractRefScalePbPerEvent(const std::string& mergeInfo, const std::string& label, double& scale)
{
  scale = 0.0;
  const std::regex pattern("\\[" + label + R"( Nraw=([0-9.eE+\-]+) sigma_pb=([0-9.eE+\-]+) w=([0-9.eE+\-]+)\])");
  std::smatch match;
  if (!std::regex_search(mergeInfo, match, pattern) || match.size() < 4) return false;

  const double nRaw = std::stod(match[1].str());
  const double sigmaPb = std::stod(match[2].str());
  if (!(nRaw > 0.0) || !(sigmaPb > 0.0)) return false;

  scale = sigmaPb / nRaw;
  return std::isfinite(scale) && scale > 0.0;
}

std::string Fixed(double value, int precision = 3)
{
  std::ostringstream os;
  os << std::fixed << std::setprecision(precision) << value;
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

void StyleRatio(TH1* h, int color = kBlack, int markerStyle = 20)
{
  if (!h) return;
  h->SetStats(false);
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(markerStyle);
  h->SetMarkerSize(0.85);
  h->SetLineWidth(2);
}

double PositiveContentNear(const TH1* h, double x)
{
  if (!h) return 1.0;
  const int ib = h->GetXaxis()->FindBin(x);
  for (int j = 0; j <= 2; ++j)
  {
    const int left = ib - j;
    const int right = ib + j;
    if (left >= 1 && left <= h->GetNbinsX() && h->GetBinContent(left) > 0.0) return h->GetBinContent(left);
    if (right >= 1 && right <= h->GetNbinsX() && h->GetBinContent(right) > 0.0) return h->GetBinContent(right);
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

std::unique_ptr<TH1> RatioToFit(const TH1* h, const TF1* fit, const std::string& name)
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

  ref.ratio = RatioToFit(h, ref.curve.get(), name + "_ratio_" + ref.tag);
  if (!ref.ratio) return ref;
  StyleRatio(ref.ratio.get());

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

std::unique_ptr<TH1> CloneRangeHist(const TH1* src, const char* name, double xLo, double xHi,
                                    int color, int markerStyle = 20)
{
  if (!src) return nullptr;
  std::unique_ptr<TH1> h(dynamic_cast<TH1*>(src->Clone(name)));
  if (!h) return nullptr;
  h->SetDirectory(nullptr);
  h->SetStats(false);
  for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
  {
    const double x = h->GetBinCenter(ib);
    const bool keep = (x >= xLo) && (x < xHi);
    if (!keep)
    {
      h->SetBinContent(ib, 0.0);
      h->SetBinError(ib, 0.0);
    }
  }
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(markerStyle);
  h->SetMarkerSize(0.90);
  h->SetLineWidth(2);
  return h;
}
}

void MakeEmbeddedInclusiveMergedStitchSpectrum(const char* inputPath = kDefaultInput.c_str())
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0.5);

  const std::string inPath = inputPath ? inputPath : kDefaultInput;
  std::unique_ptr<TFile> f(TFile::Open(inPath.c_str(), "READ"));
  if (!f || f->IsZombie())
  {
    std::cerr << "[ERROR] Cannot open " << inPath << std::endl;
    return;
  }

  TDirectory* d = f->GetDirectory("SIM");
  if (!d)
  {
    std::cerr << "[ERROR] Missing SIM directory in " << inPath << std::endl;
    return;
  }

  TH1* hIn = dynamic_cast<TH1*>(d->Get("h_embedInclusiveStitch_filterJetPt_kept"));
  if (!hIn)
  {
    std::cerr << "[ERROR] Missing SIM/h_embedInclusiveStitch_filterJetPt_kept in " << inPath << std::endl;
    return;
  }

  TNamed* mergeInfoObj = dynamic_cast<TNamed*>(d->Get("MERGE_INFO"));
  const std::string mergeInfo = mergeInfoObj ? mergeInfoObj->GetTitle() : "";
  const bool isThreeSlice = mergeInfo.find("embeddedJet30") != std::string::npos ||
                            inPath.find("embeddedJet12plus20plus30") != std::string::npos;
  const std::string refLabel = isThreeSlice ? "embeddedJet30" : "embeddedJet20";
  double refScalePbPerEvent = 1.0;
  const bool haveAbsScale = ExtractRefScalePbPerEvent(mergeInfo, refLabel, refScalePbPerEvent);
  if (!haveAbsScale)
  {
    std::cerr << "[WARN] Could not parse " << refLabel << " scale from MERGE_INFO; plotting relative merged entries." << std::endl;
  }

  std::unique_ptr<TH1> h(dynamic_cast<TH1*>(hIn->Clone("h_embeddedInclusiveJetFinalStitch_pb")));
  h->SetDirectory(nullptr);
  if (h->GetSumw2N() == 0) h->Sumw2();
  if (haveAbsScale) h->Scale(refScalePbPerEvent);
  h->Rebin(2);
  h->SetStats(false);
  h->SetTitle("");
  h->SetLineColor(kBlack);
  h->SetMarkerColor(kBlack);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.9);
  h->SetLineWidth(2);

  const double xMin = 10.0;
  const double xMax = 45.0;
  const double fitXMin = 12.0;
  const double fitXMax = xMax;
  const std::vector<double> stitchBoundaries = isThreeSlice ? std::vector<double>{20.0, 30.0}
                                                            : std::vector<double>{20.0};
  const std::string plotBaseName = "embeddedInclusiveJet_finalMergedStitchedTruthJetPtSpectrum";
  std::vector<FitReference> fitRefs;
  fitRefs.push_back(BuildFitReference(h.get(), FitFamily::kPowerLaw, plotBaseName, xMin, xMax, fitXMin, fitXMax, stitchBoundaries));
  fitRefs.push_back(BuildFitReference(h.get(), FitFamily::kHagedorn, plotBaseName, xMin, xMax, fitXMin, fitXMax, stitchBoundaries));
  fitRefs.push_back(BuildFitReference(h.get(), FitFamily::kModifiedPowerLaw, plotBaseName, xMin, xMax, fitXMin, fitXMax, stitchBoundaries));
  fitRefs.push_back(BuildFitReference(h.get(), FitFamily::kCurvedPowerLaw, plotBaseName, xMin, xMax, fitXMin, fitXMax, stitchBoundaries));
  fitRefs.push_back(BuildFitReference(h.get(), FitFamily::kLogPoly4, plotBaseName, xMin, xMax, fitXMin, fitXMax, stitchBoundaries));
  fitRefs.push_back(BuildFitReference(h.get(), FitFamily::kLogPoly6, plotBaseName, xMin, xMax, fitXMin, fitXMax, stitchBoundaries));
  fitRefs.push_back(BuildFitReference(h.get(), FitFamily::kLogPoly8, plotBaseName, xMin, xMax, fitXMin, fitXMax, stitchBoundaries));
  fitRefs.push_back(BuildFitReference(h.get(), FitFamily::kLogPoly10, plotBaseName, xMin, xMax, fitXMin, fitXMax, stitchBoundaries));
  fitRefs.push_back(BuildFitReference(h.get(), FitFamily::kLogPoly12, plotBaseName, xMin, xMax, fitXMin, fitXMax, stitchBoundaries));
  fitRefs.push_back(BuildFitReference(h.get(), FitFamily::kLogCubicSpline, plotBaseName, xMin, xMax, fitXMin, fitXMax, stitchBoundaries));
  fitRefs.push_back(BuildFitReference(h.get(), FitFamily::kPiecewisePowerLaw, plotBaseName, xMin, xMax, fitXMin, fitXMax, stitchBoundaries));
  fitRefs.push_back(BuildFitReference(h.get(), FitFamily::kPiecewiseLogQuad, plotBaseName, xMin, xMax, fitXMin, fitXMax, stitchBoundaries));
  fitRefs.push_back(BuildFitReference(h.get(), FitFamily::kExponential, plotBaseName, xMin, xMax, fitXMin, fitXMax, stitchBoundaries));
  int selectedFit = SelectFitReference(fitRefs, FitFamily::kModifiedPowerLaw);
  if (selectedFit < 0) selectedFit = SelectFitReference(fitRefs);
  if (selectedFit < 0)
  {
    std::cerr << "[WARN] No acceptable fit reference for " << plotBaseName << std::endl;
    return;
  }

  constexpr int kJet12Color = kBlue + 1;
  constexpr int kJet20Color = kOrange + 7;
  constexpr int kJet30Color = kPink + 6;
  std::unique_ptr<TH1> hJet12;
  std::unique_ptr<TH1> hJet20;
  std::unique_ptr<TH1> hJet30;
  if (isThreeSlice)
  {
    hJet12 = CloneRangeHist(h.get(), "h_embeddedInclusiveJetFinalStitch_jet12", 12.0, 20.0, kJet12Color);
    hJet20 = CloneRangeHist(h.get(), "h_embeddedInclusiveJetFinalStitch_jet20to30", 20.0, 30.0, kJet20Color);
    hJet30 = CloneRangeHist(h.get(), "h_embeddedInclusiveJetFinalStitch_jet30", 30.0, 1.0e9, kJet30Color);
  }

  h->GetXaxis()->SetRangeUser(xMin, xMax);
  h->GetYaxis()->SetTitle(haveAbsScale ? "#sigma_{eff}/N scaled entries [pb / bin]"
                                       : "weighted merged entries / bin");
  h->GetYaxis()->SetTitleOffset(1.12);
  h->SetMinimum(std::max(1.0e-8, h->GetMaximum() * 1.0e-5));
  h->SetMaximum(std::max(1.0e-6, h->GetMaximum() * 70.0));

  auto drawFitPlot = [&](FitReference& fitRef, const std::string& outBaseName)
  {
    fitRef.ratio->GetXaxis()->SetRangeUser(xMin, xMax);
    fitRef.ratio->GetXaxis()->SetTitle("max p_{T}^{jet,truth} [GeV]");
    fitRef.ratio->GetYaxis()->SetTitle("stitched / fit");
    fitRef.ratio->GetYaxis()->SetRangeUser(0.85, 1.15);
    fitRef.ratio->GetYaxis()->SetNdivisions(505);
    fitRef.ratio->GetYaxis()->SetTitleSize(0.085);
    fitRef.ratio->GetYaxis()->SetLabelSize(0.075);
    fitRef.ratio->GetYaxis()->SetTitleOffset(0.50);
    fitRef.ratio->GetXaxis()->SetTitleSize(0.090);
    fitRef.ratio->GetXaxis()->SetLabelSize(0.080);

    std::unique_ptr<TH1> rJet12;
    std::unique_ptr<TH1> rJet20;
    std::unique_ptr<TH1> rJet30;
    if (isThreeSlice)
    {
      rJet12 = CloneRangeHist(fitRef.ratio.get(), ("h_" + outBaseName + "_ratio_jet12").c_str(), 12.0, 20.0, kJet12Color);
      rJet20 = CloneRangeHist(fitRef.ratio.get(), ("h_" + outBaseName + "_ratio_jet20to30").c_str(), 20.0, 30.0, kJet20Color);
      rJet30 = CloneRangeHist(fitRef.ratio.get(), ("h_" + outBaseName + "_ratio_jet30").c_str(), 30.0, 1.0e9, kJet30Color);
    }

    TCanvas c(("c_" + outBaseName).c_str(), "embedded inclusive jet final merged stitch spectrum", 1050, 850);
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
    h->GetXaxis()->SetLabelSize(0.0);
    h->GetXaxis()->SetTitleSize(0.0);
    h->Draw("E1");
    fitRef.curve->Draw("SAME");
    if (isThreeSlice)
    {
      hJet12->Draw("E1 SAME");
      hJet20->Draw("E1 SAME");
      hJet30->Draw("E1 SAME");
    }
    else
    {
      h->Draw("E1 SAME");
    }

    TLegend leg(0.15, 0.14, 0.52, isThreeSlice ? 0.38 : 0.34);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextFont(42);
    leg.SetTextSize(isThreeSlice ? 0.032 : 0.038);
    if (isThreeSlice)
    {
      leg.AddEntry(hJet12.get(), "Jet12: 12-20 GeV", "lep");
      leg.AddEntry(hJet20.get(), "Jet20-to-30: 20-30 GeV", "lep");
      leg.AddEntry(hJet30.get(), "Jet30: #geq 30 GeV", "lep");
    }
    else
    {
      leg.AddEntry(h.get(), "Jet12+20 stitched", "lep");
    }
    leg.AddEntry(fitRef.curve.get(), fitRef.label.c_str(), "l");
    leg.Draw();

    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextFont(42);
    lat.SetTextSize(0.040);
    lat.DrawLatex(0.15, 0.86, "Embedded inclusive-jet generator stitching spectrum");
    lat.SetTextSize(0.031);
    lat.DrawLatex(0.15, 0.80, isThreeSlice ? "Final RecoilJets_embeddedJet12plus20plus30_MERGED.root"
                                            : "Final RecoilJets_embeddedJet12plus20_MERGED.root");
    if (isThreeSlice)
    {
      lat.DrawLatex(0.15, 0.75, "Jet12: 12 #leq max p_{T}^{jet,truth} < 20 GeV; Jet20: 20-30 GeV; Jet30: #geq 30 GeV");
    }
    else
    {
      lat.DrawLatex(0.15, 0.75, "Jet12: 12 #leq max p_{T}^{jet,truth} < 20 GeV; Jet20: max p_{T}^{jet,truth} #geq 20 GeV");
    }
    if (haveAbsScale)
    {
      lat.DrawLatex(0.15, 0.70, ("Final stitch scale = " + Sci(refScalePbPerEvent) + " pb per " +
                                 (isThreeSlice ? "Jet30" : "Jet20") + "-weighted event").c_str());
    }
    else
    {
      lat.DrawLatex(0.15, 0.70, "MERGE_INFO scale unavailable; shown in relative final-stitch units");
    }
    TLatex sph;
    sph.SetNDC(true);
    sph.SetTextFont(42);
    sph.SetTextAlign(33);
    sph.SetTextSize(0.042);
    sph.DrawLatex(0.92, 0.58, "#it{#bf{sPHENIX}} Internal");
    sph.SetTextSize(0.034);
    sph.DrawLatex(0.92, 0.53, "Pythia Overlay #sqrt{s_{NN}} = 200 GeV");

    bot.cd();
    fitRef.ratio->Draw("E1");
    if (isThreeSlice)
    {
      rJet12->Draw("E1 SAME");
      rJet20->Draw("E1 SAME");
      rJet30->Draw("E1 SAME");
    }
    TLine one(xMin, 1.0, xMax, 1.0);
    one.SetLineColor(kGray + 1);
    one.SetLineStyle(2);
    one.Draw("SAME");

    const std::string outPng = ParentDir(inPath) + "/" + outBaseName + ".png";
    c.SaveAs(outPng.c_str());
    std::cout << "[DONE] Wrote " << outPng << std::endl;
  };

  for (auto& fitRef : fitRefs)
  {
    if (!fitRef.ok) continue;
    drawFitPlot(fitRef, plotBaseName + "_fit_" + fitRef.tag);
  }
  drawFitPlot(fitRefs[selectedFit], plotBaseName);

  std::cout << "[FIT] selected " << plotBaseName << " -> " << fitRefs[selectedFit].tag << std::endl;
  std::cout << "[INFO] input=" << inPath << std::endl;
  std::cout << "[INFO] mergeInfo=" << mergeInfo << std::endl;
  std::cout << "[INFO] refScalePbPerEvent=" << refScalePbPerEvent << std::endl;
}
