#include "sPhenixStyle.C"

#include <TBox.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace
{
struct Moments
{
  double mean = 0.0;
  double sigma = 0.0;
};

struct FeatureSep
{
  std::string feature;
  std::string label;
  double sigMean = 0.0;
  double bkgMean = 0.0;
  double separation = 0.0;
};

std::vector<std::string> SplitCsvLine(const std::string& line)
{
  std::vector<std::string> out;
  std::string token;
  std::stringstream ss(line);
  while (std::getline(ss, token, ',')) out.push_back(token);
  return out;
}

double ToDouble(const std::string& s)
{
  try
  {
    return std::stod(s);
  }
  catch (...)
  {
    return 0.0;
  }
}

bool StartsWith(const std::string& s, const std::string& prefix)
{
  return s.size() >= prefix.size() && s.compare(0, prefix.size(), prefix) == 0;
}

bool Contains(const std::string& s, const std::string& needle)
{
  return s.find(needle) != std::string::npos;
}

std::map<std::string, std::map<std::string, Moments>>
ReadInclusiveMoments(const std::string& csv)
{
  std::ifstream in(csv);
  std::map<std::string, std::map<std::string, Moments>> out;
  std::string line;
  std::getline(in, line);
  while (std::getline(in, line))
  {
    auto cols = SplitCsvLine(line);
    if (cols.size() < 7) continue;
    const std::string& feature = cols[0];
    const std::string& scope = cols[1];
    const std::string& cls = cols[2];
    if (scope != "inclusive") continue;
    out[feature][cls] = {ToDouble(cols[4]), ToDouble(cols[5])};
  }
  return out;
}

std::unique_ptr<TH1> SumInclusiveTemplates(TDirectory* dir,
                                           const std::string& var,
                                           const std::string& truthTag,
                                           const std::string& name)
{
  if (!dir) return nullptr;
  const std::string prefix = "h_ss_" + var + "_pre_" + truthTag + "_pT_";
  std::unique_ptr<TH1> sum;

  for (auto* keyObj : *dir->GetListOfKeys())
  {
    const std::string key = keyObj->GetName();
    if (!StartsWith(key, prefix) || Contains(key, "_cent_")) continue;
    TH1* h = dynamic_cast<TH1*>(dir->Get(key.c_str()));
    if (!h) continue;
    if (!sum)
    {
      sum.reset(dynamic_cast<TH1*>(h->Clone(name.c_str())));
      sum->SetDirectory(nullptr);
      sum->Reset("ICES");
    }
    sum->Add(h);
  }
  return sum;
}

Moments HistMoments(TDirectory* dir, const std::string& var, const std::string& truthTag)
{
  auto h = SumInclusiveTemplates(dir, var, truthTag,
                                 Form("h_sum_%s_%s", var.c_str(), truthTag.c_str()));
  if (!h || h->Integral() <= 0.0) return {};
  return {h->GetMean(), h->GetStdDev()};
}

double Separation(const Moments& sig, const Moments& bkg)
{
  const double denom = std::sqrt(std::max(1.0e-12, 0.5 * (sig.sigma * sig.sigma + bkg.sigma * bkg.sigma)));
  return std::abs(sig.mean - bkg.mean) / denom;
}

std::string PrettyLabel(const std::string& feature)
{
  if (feature == "cluster_Et") return "cluster E_{T}";
  if (feature == "cluster_et1") return "leading tower fraction";
  if (feature == "cluster_et2") return "second tower fraction";
  if (feature == "cluster_et3") return "third tower fraction";
  if (feature == "cluster_et4") return "fourth tower fraction";
  if (feature == "e11_over_e33") return "E_{1x1}/E_{3x3}";
  if (feature == "e32_over_e35") return "E_{3x2}/E_{3x5}";
  if (feature == "cluster_weta33_cogx") return "3x3 #eta shower width";
  if (feature == "cluster_wphi33_cogx") return "3x3 #phi shower width";
  return feature;
}
}

void MakeAuAu3x3FeatureSeparation(
  const char* reportDir = "dataOutput/auauTightBDTValidation/model_validation_condor_cent3x3_20260510_221419",
  const char* outDir = "dataOutput/auauTightBDTValidation/model_validation_condor_cent3x3_20260510_221419/money_plots_3x3")
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  const std::string report(reportDir);
  const std::string out(outDir);
  auto csvMoments = ReadInclusiveMoments(report + "/validation_feature_summary.csv");

  const std::string sigPath =
    "dataOutput/combinedSimOnlyEMBEDDED/preselectionNewPPG12_tightAuAuCentInputBDT_nonTightAuAuBDTComplement_baseVariant/photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root";
  const std::string bkgPath =
    "dataOutput/combinedSimOnlyEMBEDDED/preselectionNewPPG12_tightAuAuCentInputBDT_nonTightAuAuBDTComplement_baseVariant/embeddedJet12and20merged_SIM/RecoilJets_embeddedJet12plus20_MERGED.root";

  std::unique_ptr<TFile> sigFile(TFile::Open(sigPath.c_str(), "READ"));
  std::unique_ptr<TFile> bkgFile(TFile::Open(bkgPath.c_str(), "READ"));
  TDirectory* sigDir = sigFile ? dynamic_cast<TDirectory*>(sigFile->Get("SIM")) : nullptr;
  TDirectory* bkgDir = bkgFile ? dynamic_cast<TDirectory*>(bkgFile->Get("SIM")) : nullptr;

  std::vector<std::string> features = {
    "cluster_Et",
    "e32_over_e35",
    "cluster_et1",
    "cluster_wphi33_cogx",
    "cluster_weta33_cogx",
    "e11_over_e33",
    "cluster_et3",
    "cluster_et2",
    "cluster_et4",
    "centrality",
  };

  std::vector<FeatureSep> rows;
  for (const auto& feature : features)
  {
    Moments sig;
    Moments bkg;
    if (feature == "cluster_weta33_cogx")
    {
      sig = HistMoments(sigDir, "weta33", "sig");
      bkg = HistMoments(bkgDir, "weta33", "bkg");
    }
    else if (feature == "cluster_wphi33_cogx")
    {
      sig = HistMoments(sigDir, "wphi33", "sig");
      bkg = HistMoments(bkgDir, "wphi33", "bkg");
    }
    else
    {
      sig = csvMoments[feature]["signal"];
      bkg = csvMoments[feature]["background"];
    }
    rows.push_back({feature, PrettyLabel(feature), sig.mean, bkg.mean, Separation(sig, bkg)});
  }

  std::sort(rows.begin(), rows.end(), [](const FeatureSep& a, const FeatureSep& b) {
    return a.separation > b.separation;
  });

  gSystem->mkdir(out.c_str(), true);
  std::ofstream table(out + "/bdt3x3_feature_separation_inputs.csv");
  table << "feature,signal_mean,background_mean,separation\n";
  for (const auto& row : rows)
  {
    table << row.feature << "," << row.sigMean << "," << row.bkgMean << "," << row.separation << "\n";
  }

  const int n = std::min<int>(10, rows.size());
  double xmax = 0.0;
  for (int i = 0; i < n; ++i) xmax = std::max(xmax, rows[i].separation);
  xmax *= 1.14;

  TCanvas c("c_bdt3x3_feature_separation", "", 1600, 1000);
  c.SetMargin(0.27, 0.04, 0.14, 0.12);
  TH2F frame("frame", ";Signal/background separation of input distribution;", 1, 0.0, xmax, n, 0.0, n);
  frame.SetStats(false);
  frame.GetXaxis()->SetTitleSize(0.043);
  frame.GetXaxis()->SetLabelSize(0.038);
  frame.GetYaxis()->SetLabelSize(0.0);
  frame.GetYaxis()->SetTickLength(0.0);
  frame.GetYaxis()->SetNdivisions(0);
  frame.Draw("AXIS");

  TLatex label;
  label.SetNDC(false);
  label.SetTextFont(42);
  label.SetTextAlign(32);
  label.SetTextSize(0.035);

  const int barColor = kAzure + 2;
  std::vector<std::unique_ptr<TBox>> bars;
  for (int ir = 0; ir < n; ++ir)
  {
    const auto& row = rows[n - 1 - ir];
    const double y0 = ir + 0.18;
    const double y1 = ir + 0.82;
    auto box = std::make_unique<TBox>(0.0, y0, row.separation, y1);
    box->SetFillColor(barColor);
    box->SetLineColor(barColor);
    box->Draw();
    bars.push_back(std::move(box));
    label.DrawLatex(-0.03 * xmax, ir + 0.50, row.label.c_str());
  }

  TLine grid;
  grid.SetLineColorAlpha(kGray + 1, 0.25);
  grid.SetLineStyle(1);
  for (double x = 0.5; x < xmax; x += 0.5) grid.DrawLine(x, 0.0, x, n);

  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(23);
  tx.SetTextSize(0.040);
  tx.DrawLatex(0.54, 0.955, "Most separating inputs before the 3#times3-width BDT");
  tx.SetTextAlign(13);
  tx.SetTextSize(0.033);
  tx.DrawLatex(0.70, 0.27, "#it{#bf{sPHENIX}} Internal");
  tx.SetTextSize(0.028);
  tx.DrawLatex(0.70, 0.22, "Pythia overlay, #sqrt{s_{NN}} = 200 GeV");

  const std::string png = out + "/bdt3x3_feature_separation_inputs.png";
  c.SaveAs(png.c_str());
  std::cout << "[MakeAuAu3x3FeatureSeparation] wrote " << png << "\n";
}
