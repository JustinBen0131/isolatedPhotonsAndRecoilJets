#include "sPhenixStyle.C"

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

namespace
{
struct Product
{
  std::string key;
  std::string label;
  int color;
  int marker;
};

const std::vector<Product> kProducts = {
  {"centINDcontrol", "Cluster-shape only", kAzure + 1, 20},
  {"centAsFeat", "Cluster-shape + centrality", kOrange + 7, 21},
  {"centDepBDTs", "Centrality-specific models", kGreen + 2, 22},
};

const std::vector<std::pair<std::string, std::string>> kCentBins = {
  {"0_20", "0-20%"},
  {"20_50", "20-50%"},
  {"50_80", "50-80%"},
};

std::string Slurp(const std::string& path)
{
  std::ifstream in(path);
  std::ostringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

double ExtractMetric(const std::string& json, const std::string& product, const std::string& field)
{
  const std::string needle = "\"" + product + "\"";
  const auto pos = json.find(needle);
  if (pos == std::string::npos) return std::numeric_limits<double>::quiet_NaN();
  auto end = json.size();
  for (const auto& p : kProducts)
  {
    if (p.key == product) continue;
    const auto next = json.find("\"" + p.key + "\"", pos + needle.size());
    if (next != std::string::npos) end = std::min(end, next);
  }
  const std::string body = json.substr(pos, end - pos);
  const std::regex valueRe("\"" + field + "\"\\s*:\\s*([-+0-9.eE]+)");
  std::smatch value;
  if (!std::regex_search(body, value, valueRe)) return std::numeric_limits<double>::quiet_NaN();
  return std::stod(value.str(1));
}

double ExtractCentralityAuc(const std::string& json, const std::string& product, const std::string& centKey)
{
  const std::string needle = "\"" + product + "\"";
  const auto pos = json.find(needle);
  if (pos == std::string::npos) return std::numeric_limits<double>::quiet_NaN();
  auto end = json.size();
  for (const auto& p : kProducts)
  {
    if (p.key == product) continue;
    const auto next = json.find("\"" + p.key + "\"", pos + needle.size());
    if (next != std::string::npos) end = std::min(end, next);
  }
  const std::string productBody = json.substr(pos, end - pos);
  const std::regex centRe("\"auc_by_centrality\"\\s*:\\s*\\{([\\s\\S]*?)\\}");
  std::smatch centBlock;
  if (!std::regex_search(productBody, centBlock, centRe)) return std::numeric_limits<double>::quiet_NaN();
  const std::string centBody = centBlock.str(1);
  const std::regex valueRe("\"" + centKey + "\"\\s*:\\s*([-+0-9.eE]+)");
  std::smatch value;
  if (!std::regex_search(centBody, value, valueRe)) return std::numeric_limits<double>::quiet_NaN();
  return std::stod(value.str(1));
}

void DrawSphenixLabel(double x, double y, const char* line2 = "Pythia overlay, #sqrt{s_{NN}} = 200 GeV")
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextSize(0.040);
  tx.DrawLatex(x, y, "#it{#bf{sPHENIX}} Internal");
  tx.SetTextSize(0.034);
  tx.DrawLatex(x, y - 0.050, line2);
  tx.DrawLatex(x, y - 0.092, "Embedded #gamma+jet signal and inclusive-jet background");
}

void DrawSphenixHeader(double x, double y)
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextSize(0.034);
  tx.DrawLatex(x, y, "#it{#bf{sPHENIX}} Internal  |  Pythia overlay, #sqrt{s_{NN}} = 200 GeV");
  tx.SetTextSize(0.030);
  tx.DrawLatex(x, y - 0.040, "Embedded #gamma+jet signal and inclusive-jet background");
}

void Save(TCanvas& c, const std::string& path)
{
  c.SaveAs(path.c_str());
  std::cout << "[MakeAuAuTightBDTValidationSlidePlots] wrote " << path << "\n";
}
}

void MakeAuAuTightBDTValidationSlidePlots(
    const char* reportDir = "dataOutput/auauTightBDTValidation/model_validation_20260508_115524",
    const char* outDir = "dataOutput/auauTightBDTValidation/model_validation_20260508_115524/slideReady")
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);

  const std::string report(reportDir);
  const std::string out(outDir);
  gSystem->mkdir(out.c_str(), true);

  std::unique_ptr<TFile> f(TFile::Open((report + "/validation_curves.root").c_str(), "READ"));
  if (!f || f->IsZombie())
  {
    std::cerr << "Cannot open " << report << "/validation_curves.root\n";
    return;
  }
  const std::string metrics = Slurp(report + "/validation_metrics.json");

  {
    TCanvas c("c_auauTightBDT_roc", "", 950, 760);
    c.SetLeftMargin(0.12);
    c.SetRightMargin(0.04);
    c.SetBottomMargin(0.12);
    c.SetTopMargin(0.16);

    TH2F frame("frame_roc", ";Background fake rate;Signal efficiency", 100, 0, 1, 100, 0, 1.02);
    frame.GetXaxis()->SetTitleOffset(1.05);
    frame.GetYaxis()->SetTitleOffset(1.10);
    frame.Draw();

    TLine diag(0, 0, 1, 1);
    diag.SetLineStyle(7);
    diag.SetLineColor(kGray + 2);
    diag.Draw();

    TLegend leg(0.43, 0.20, 0.90, 0.40);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextFont(42);
    leg.SetTextSize(0.033);

    for (const auto& p : kProducts)
    {
      auto* g = dynamic_cast<TGraph*>(f->Get(("gROC_" + p.key + "_signalEfficiency_vs_backgroundFakeRate").c_str()));
      if (!g) continue;
      g->SetLineColor(p.color);
      g->SetLineWidth(3);
      g->SetMarkerColor(p.color);
      g->SetMarkerStyle(p.marker);
      g->SetMarkerSize(0.45);
      g->Draw("L same");
      const double auc = ExtractMetric(metrics, p.key, "auc_inclusive");
      leg.AddEntry(g, Form("%s  AUC = %.3f", p.label.c_str(), auc), "l");
    }
    leg.Draw();
    DrawSphenixHeader(0.14, 0.955);

    TLatex title;
    title.SetNDC();
    title.SetTextFont(42);
    title.SetTextSize(0.033);
    title.DrawLatex(0.14, 0.855, "ROC scan for Au+Au tight photon-ID BDT candidates");

    Save(c, out + "/auau_tight_bdt_roc_inclusive_sphenix.png");
  }

  {
    TCanvas c("c_auauTightBDT_score", "", 980, 760);
    c.SetLogy();
    c.SetLeftMargin(0.12);
    c.SetRightMargin(0.04);
    c.SetBottomMargin(0.12);
    c.SetTopMargin(0.07);

    auto* hSig = dynamic_cast<TH1*>(f->Get("hScore_centDepBDTs_signal_unitArea"));
    auto* hBkg = dynamic_cast<TH1*>(f->Get("hScore_centDepBDTs_background_unitArea"));
    if (hSig && hBkg)
    {
      hSig->SetDirectory(nullptr);
      hBkg->SetDirectory(nullptr);
      hSig->SetStats(false);
      hBkg->SetStats(false);
      hSig->SetLineColor(kOrange + 7);
      hSig->SetLineWidth(3);
      hBkg->SetLineColor(kAzure + 1);
      hBkg->SetFillColorAlpha(kAzure + 1, 0.30);
      hBkg->SetLineWidth(2);
      hBkg->SetTitle(";BDT score;Unit-normalized candidates");
      hBkg->GetYaxis()->SetRangeUser(3e-3, 20);
      hBkg->GetXaxis()->SetTitleOffset(1.05);
      hBkg->GetYaxis()->SetTitleOffset(1.10);
      hBkg->Draw("hist");
      hSig->Draw("hist same");

      TLegend leg(0.58, 0.18, 0.91, 0.31);
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);
      leg.SetTextFont(42);
      leg.SetTextSize(0.035);
      leg.AddEntry(hSig, "Truth-matched photons", "l");
      leg.AddEntry(hBkg, "Inclusive-jet candidates", "f");
      leg.Draw();

      DrawSphenixLabel(0.17, 0.88);

      TLatex note;
      note.SetNDC();
      note.SetTextFont(42);
      note.SetTextSize(0.035);
      note.DrawLatex(0.17, 0.72, "Centrality-specific model output");
      note.DrawLatex(0.17, 0.67, Form("#LTscore#GT_{sig}=%.3f,  #LTscore#GT_{bkg}=%.3f",
                                        ExtractMetric(metrics, "centDepBDTs", "signal_score_mean"),
                                        ExtractMetric(metrics, "centDepBDTs", "background_score_mean")));

      Save(c, out + "/auau_tight_bdt_score_separation_centDep_sphenix.png");
    }
  }

  {
    TCanvas c("c_auauTightBDT_aucByCent", "", 980, 760);
    c.SetLeftMargin(0.12);
    c.SetRightMargin(0.04);
    c.SetBottomMargin(0.12);
    c.SetTopMargin(0.07);

    TH2F frame("frame_auc_cent", ";Centrality;ROC AUC", 3, 0.5, 3.5, 100, 0.84, 0.91);
    for (int i = 0; i < 3; ++i) frame.GetXaxis()->SetBinLabel(i + 1, kCentBins[i].second.c_str());
    frame.GetXaxis()->SetTitleOffset(1.05);
    frame.GetYaxis()->SetTitleOffset(1.10);
    frame.Draw();

    TLegend leg(0.43, 0.19, 0.90, 0.37);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextFont(42);
    leg.SetTextSize(0.033);

    std::vector<std::unique_ptr<TGraph>> keep;
    for (const auto& p : kProducts)
    {
      auto g = std::make_unique<TGraph>();
      for (int i = 0; i < 3; ++i)
      {
        const double auc = ExtractCentralityAuc(metrics, p.key, kCentBins[i].first);
        g->SetPoint(i, i + 1, auc);
      }
      g->SetLineColor(p.color);
      g->SetMarkerColor(p.color);
      g->SetLineWidth(3);
      g->SetMarkerStyle(p.marker);
      g->SetMarkerSize(1.3);
      g->Draw("LP same");
      leg.AddEntry(g.get(), p.label.c_str(), "lp");
      keep.push_back(std::move(g));
    }
    leg.Draw();
    DrawSphenixLabel(0.17, 0.88);

    TLatex title;
    title.SetNDC();
    title.SetTextFont(42);
    title.SetTextSize(0.038);
    title.DrawLatex(0.17, 0.74, "Model separation versus centrality");

    Save(c, out + "/auau_tight_bdt_auc_by_centrality_sphenix.png");
  }
}
