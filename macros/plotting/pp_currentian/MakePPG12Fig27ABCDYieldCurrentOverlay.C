#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace
{
struct Region
{
  const char *letter;
  const char *hist;
  const char *label;
  int color;
  int marker;
  double offset;
};

const Region kRegions[] = {
    {"A", "h_tight_iso_cluster_0", "A: tight iso", kBlack, 20, -0.18},
    {"B", "h_tight_noniso_cluster_0", "B: tight noniso", kRed + 1, 21, -0.06},
    {"C", "h_nontight_iso_cluster_0", "C: nontight iso", kBlue + 1, 22, 0.06},
    {"D", "h_nontight_noniso_cluster_0", "D: nontight noniso", kMagenta + 1, 23, 0.18},
};

TH1 *requireHist(TDirectory *dir, const char *name, const char *label)
{
  if (!dir)
  {
    Error("MakePPG12Fig27ABCDYieldCurrentOverlay", "Null directory for %s", label);
    return nullptr;
  }
  auto *hist = dynamic_cast<TH1 *>(dir->Get(name));
  if (!hist)
  {
    Error("MakePPG12Fig27ABCDYieldCurrentOverlay", "Missing %s histogram %s", label, name);
    return nullptr;
  }
  auto *clone = dynamic_cast<TH1 *>(hist->Clone(Form("%s_%s_clone", label, name)));
  clone->SetDirectory(nullptr);
  clone->SetStats(false);
  return clone;
}

bool sameBinning(const TH1 *a, const TH1 *b)
{
  if (!a || !b || a->GetNbinsX() != b->GetNbinsX())
    return false;
  for (int i = 1; i <= a->GetNbinsX(); ++i)
  {
    if (std::fabs(a->GetXaxis()->GetBinLowEdge(i) - b->GetXaxis()->GetBinLowEdge(i)) > 1.0e-6)
      return false;
    if (std::fabs(a->GetXaxis()->GetBinUpEdge(i) - b->GetXaxis()->GetBinUpEdge(i)) > 1.0e-6)
      return false;
  }
  return true;
}

TGraphErrors *makeYieldGraph(const TH1 *hist, const Region &region, bool current)
{
  const int n = hist->GetNbinsX();
  auto *graph = new TGraphErrors(n);
  graph->SetName(Form("g_%s_%s", current ? "current" : "ppg12", region.letter));
  for (int i = 1; i <= n; ++i)
  {
    const double x = hist->GetXaxis()->GetBinCenter(i);
    const double y = hist->GetBinContent(i);
    graph->SetPoint(i - 1, x, y);
    graph->SetPointError(i - 1, 0.0, hist->GetBinError(i));
  }
  graph->SetLineColor(region.color);
  graph->SetMarkerColor(region.color);
  graph->SetMarkerStyle(current ? region.marker : 24);
  graph->SetMarkerSize(current ? 1.05 : 1.15);
  graph->SetLineWidth(2);
  return graph;
}

TGraphErrors *makeRatioGraph(const TH1 *current, const TH1 *ppg12, const Region &region)
{
  const int n = current->GetNbinsX();
  auto *graph = new TGraphErrors(n);
  graph->SetName(Form("g_ratio_%s", region.letter));
  for (int i = 1; i <= n; ++i)
  {
    const double width = current->GetXaxis()->GetBinWidth(i);
    const double x = current->GetXaxis()->GetBinCenter(i) + region.offset * width;
    const double c = current->GetBinContent(i);
    const double p = ppg12->GetBinContent(i);
    const double ec = current->GetBinError(i);
    const double ep = ppg12->GetBinError(i);
    double r = std::numeric_limits<double>::quiet_NaN();
    double er = 0.0;
    if (p > 0.0)
    {
      r = c / p;
      const double t1 = ec / p;
      const double t2 = c * ep / (p * p);
      er = std::sqrt(t1 * t1 + t2 * t2);
    }
    graph->SetPoint(i - 1, x, r);
    graph->SetPointError(i - 1, 0.0, er);
  }
  graph->SetLineColor(region.color);
  graph->SetMarkerColor(region.color);
  graph->SetMarkerStyle(region.marker);
  graph->SetMarkerSize(1.0);
  graph->SetLineWidth(2);
  return graph;
}

void setupGraphFrame(TH1 *frame)
{
  frame->SetStats(false);
  frame->GetXaxis()->SetRangeUser(10.0, 36.0);
  frame->GetXaxis()->SetTitle("E_{T}^{#gamma,reco} [GeV]");
  frame->GetYaxis()->SetTitle("raw ABCD yield");
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetXaxis()->SetTitleSize(0);
  frame->GetYaxis()->SetTitleSize(0.062);
  frame->GetYaxis()->SetTitleOffset(0.82);
  frame->GetYaxis()->SetLabelSize(0.052);
}

double ratioPanelYMax(const std::vector<TGraphErrors *> &graphs, double ymin)
{
  double ymax = 1.20;
  for (auto *graph : graphs)
  {
    if (!graph)
      continue;
    for (int i = 0; i < graph->GetN(); ++i)
    {
      double x = 0.0;
      double y = 0.0;
      graph->GetPoint(i, x, y);
      if (!std::isfinite(y))
        continue;
      const double upper = y + graph->GetErrorY(i);
      if (std::isfinite(upper))
        ymax = std::max(ymax, upper);
    }
  }
  const double buffer = std::max(0.06, 0.08 * (ymax - ymin));
  return std::ceil((ymax + buffer) * 20.0) / 20.0;
}

void setupRatioFrame(TH1 *frame, double ymax)
{
  frame->SetStats(false);
  frame->GetXaxis()->SetRangeUser(10.0, 36.0);
  frame->GetYaxis()->SetRangeUser(0.20, ymax);
  frame->GetXaxis()->SetTitle("E_{T}^{#gamma,reco} [GeV]");
  frame->GetYaxis()->SetTitle("Current / PPG12");
  frame->GetXaxis()->SetTitleSize(0.13);
  frame->GetXaxis()->SetTitleOffset(0.92);
  frame->GetXaxis()->SetLabelSize(0.105);
  frame->GetYaxis()->SetTitleSize(0.095);
  frame->GetYaxis()->SetTitleOffset(0.54);
  frame->GetYaxis()->SetLabelSize(0.085);
  frame->GetYaxis()->SetNdivisions(505);
}

std::string fmt(double value, int precision = 6)
{
  std::ostringstream os;
  os << std::setprecision(precision) << std::defaultfloat << value;
  return os.str();
}
}  // namespace

void MakePPG12Fig27ABCDYieldCurrentOverlay(
    const char *current_root =
        "dataOutput/ppg12Parity/the76_ppg12_parity_full_20260701_003024/final_roots/pp_data_hierarchical_v1/RecoilJets_pp_ALL_PERIOD_COMBINED.root",
    const char *ppg12_root =
        "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620/reference_roots/fig27_abcd_yield/data_histo_bdt_nom.root",
    const char *current_dir_name = "PPG12_scaledtrigger30",
    const char *output_dir =
        "dataOutput/ppg12Parity/the76_ppg12_parity_full_20260701_003024/data_abcd_fig27",
    const char *current_label = "July pp output")
{
  gSystem->mkdir(output_dir, kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLineWidth(2);
  gStyle->SetEndErrorSize(3);

  TFile fCurrent(current_root, "READ");
  TFile fPpg12(ppg12_root, "READ");
  if (fCurrent.IsZombie() || fCurrent.TestBit(TFile::kRecovered))
  {
    Error("MakePPG12Fig27ABCDYieldCurrentOverlay",
          "Current ROOT is zombie or recovered: %s", current_root);
    return;
  }
  if (fPpg12.IsZombie() || fPpg12.TestBit(TFile::kRecovered))
  {
    Error("MakePPG12Fig27ABCDYieldCurrentOverlay",
          "PPG12 ROOT is zombie or recovered: %s", ppg12_root);
    return;
  }
  auto *currentDir = dynamic_cast<TDirectory *>(fCurrent.Get(current_dir_name));
  if (!currentDir)
  {
    Error("MakePPG12Fig27ABCDYieldCurrentOverlay", "Missing current directory %s", current_dir_name);
    return;
  }

  std::vector<TH1 *> currentHists;
  std::vector<TH1 *> ppg12Hists;
  for (const auto &region : kRegions)
  {
    TH1 *hc = requireHist(currentDir, region.hist, Form("current_%s", region.letter));
    TH1 *hp = requireHist(&fPpg12, region.hist, Form("ppg12_%s", region.letter));
    if (!hc || !hp)
      return;
    if (!sameBinning(hc, hp))
    {
      Error("MakePPG12Fig27ABCDYieldCurrentOverlay", "Binning mismatch for region %s", region.letter);
      return;
    }
    currentHists.push_back(hc);
    ppg12Hists.push_back(hp);
  }

  std::vector<TGraphErrors *> currentGraphs;
  std::vector<TGraphErrors *> ppg12Graphs;
  std::vector<TGraphErrors *> ratioGraphs;
  for (int i = 0; i < 4; ++i)
  {
    currentGraphs.push_back(makeYieldGraph(currentHists[i], kRegions[i], true));
    ppg12Graphs.push_back(makeYieldGraph(ppg12Hists[i], kRegions[i], false));
    ratioGraphs.push_back(makeRatioGraph(currentHists[i], ppg12Hists[i], kRegions[i]));
  }

  auto *canvas = new TCanvas("c_fig27_abcd_current_vs_ppg12", "c_fig27_abcd_current_vs_ppg12", 900, 1050);
  canvas->SetMargin(0, 0, 0, 0);
  auto *padTop = new TPad("padTop", "padTop", 0, 0.31, 1, 1);
  auto *padBot = new TPad("padBot", "padBot", 0, 0, 1, 0.31);
  padTop->SetBottomMargin(0.015);
  padTop->SetTopMargin(0.045);
  padTop->SetLeftMargin(0.15);
  padTop->SetRightMargin(0.035);
  padBot->SetTopMargin(0.015);
  padBot->SetBottomMargin(0.30);
  padBot->SetLeftMargin(0.15);
  padBot->SetRightMargin(0.035);
  padTop->SetTicks(1, 1);
  padBot->SetTicks(1, 1);
  padTop->Draw();
  padBot->Draw();

  padTop->cd();
  padTop->SetLogy();
  auto *frameTop = dynamic_cast<TH1 *>(ppg12Hists[0]->Clone("frameTop"));
  frameTop->Reset("ICES");
  frameTop->GetYaxis()->SetRangeUser(1.0, 2.5e5);
  setupGraphFrame(frameTop);
  frameTop->Draw("AXIS");
  for (auto *g : currentGraphs)
    g->Draw("PZ SAME");
  for (auto *g : ppg12Graphs)
    g->Draw("PZ SAME");

  TLatex text;
  text.SetNDC(true);
  text.SetTextFont(42);
  text.SetTextSize(0.042);
  text.DrawLatex(0.31, 0.91, "#it{#bf{sPHENIX}} Internal");
  text.SetTextSize(0.035);
  text.DrawLatex(0.31, 0.855, "p+p #sqrt{s} = 200 GeV");
  text.DrawLatex(0.31, 0.810, "|#eta^{#gamma}| < 0.7");
  text.DrawLatex(0.31, 0.765, "Data, bdt_nom");
  text.DrawLatex(0.31, 0.720, "raw ABCD yields; no scale factor");

  auto *legSource = new TLegend(0.66, 0.79, 0.93, 0.91);
  legSource->SetBorderSize(0);
  legSource->SetFillStyle(0);
  legSource->SetTextFont(42);
  legSource->SetTextSize(0.033);
  legSource->AddEntry(ppg12Graphs[0], "PPG12 SDCC ROOT", "pe");
  legSource->AddEntry(currentGraphs[0], current_label, "pe");
  legSource->Draw();

  auto *legRegion = new TLegend(0.66, 0.55, 0.93, 0.77);
  legRegion->SetBorderSize(0);
  legRegion->SetFillStyle(0);
  legRegion->SetTextFont(42);
  legRegion->SetTextSize(0.030);
  for (int i = 0; i < 4; ++i)
    legRegion->AddEntry(currentGraphs[i], kRegions[i].label, "pe");
  legRegion->Draw();
  padTop->RedrawAxis();

  padBot->cd();
  auto *frameRatio = dynamic_cast<TH1 *>(ppg12Hists[0]->Clone("frameRatio"));
  frameRatio->Reset("ICES");
  const double ratioYMax = ratioPanelYMax(ratioGraphs, 0.20);
  setupRatioFrame(frameRatio, ratioYMax);
  frameRatio->Draw("AXIS");
  TLine lineOne(10.0, 1.0, 36.0, 1.0);
  lineOne.SetLineStyle(2);
  lineOne.SetLineColor(kGray + 2);
  lineOne.SetLineWidth(2);
  lineOne.Draw();
  for (auto *g : ratioGraphs)
    g->Draw("PZ SAME");
  padBot->RedrawAxis();

  const TString pngPath = Form("%s/fig27_abcd_yield_current_vs_ppg12_overlay_ratio.png", output_dir);
  canvas->SaveAs(pngPath);

  const TString csvPath = Form("%s/fig27_abcd_yield_current_vs_ppg12_overlay_ratio.csv", output_dir);
  std::ofstream csv(csvPath.Data());
  csv << "region,label,pt_lo,pt_hi,ppg12_yield,ppg12_error,current_yield,current_error,current_over_ppg12,current_over_ppg12_error\n";
  double maxAbsDev = -1.0;
  std::string maxDevWhere;
  for (int r = 0; r < 4; ++r)
  {
    TH1 *hc = currentHists[r];
    TH1 *hp = ppg12Hists[r];
    for (int i = 1; i <= hp->GetNbinsX(); ++i)
    {
      const double c = hc->GetBinContent(i);
      const double p = hp->GetBinContent(i);
      const double ec = hc->GetBinError(i);
      const double ep = hp->GetBinError(i);
      double ratio = std::numeric_limits<double>::quiet_NaN();
      double eratio = std::numeric_limits<double>::quiet_NaN();
      if (p > 0.0)
      {
        ratio = c / p;
        const double t1 = ec / p;
        const double t2 = c * ep / (p * p);
        eratio = std::sqrt(t1 * t1 + t2 * t2);
        const double dev = std::fabs(ratio - 1.0);
        if (dev > maxAbsDev)
        {
          maxAbsDev = dev;
          maxDevWhere = Form("%s %.0f-%.0f GeV", kRegions[r].letter,
                             hp->GetXaxis()->GetBinLowEdge(i), hp->GetXaxis()->GetBinUpEdge(i));
        }
      }
      csv << kRegions[r].letter << "," << kRegions[r].label << ","
          << hp->GetXaxis()->GetBinLowEdge(i) << "," << hp->GetXaxis()->GetBinUpEdge(i) << ","
          << p << "," << ep << "," << c << "," << ec << "," << ratio << "," << eratio << "\n";
    }
  }
  csv.close();

  const TString manifestPath = Form("%s/fig27_abcd_yield_current_vs_ppg12_overlay_ratio_manifest.json", output_dir);
  std::ofstream manifest(manifestPath.Data());
  manifest << "{\n";
  manifest << "  \"schema\": \"THE93_FIG27_ABCD_YIELD_CURRENT_VS_PPG12_V1\",\n";
  manifest << "  \"png\": \"" << pngPath.Data() << "\",\n";
  manifest << "  \"csv\": \"" << csvPath.Data() << "\",\n";
  manifest << "  \"current_root\": \"" << current_root << "\",\n";
  manifest << "  \"current_label\": \"" << current_label << "\",\n";
  manifest << "  \"current_dir\": \"" << current_dir_name << "\",\n";
  manifest << "  \"ppg12_root\": \"" << ppg12_root << "\",\n";
  manifest << "  \"ppg12_remote_source\": \"/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo_bdt_nom.root\",\n";
  manifest << "  \"histograms\": [\n";
  for (int i = 0; i < 4; ++i)
  {
    manifest << "    {\"region\": \"" << kRegions[i].letter << "\", \"hist\": \"" << kRegions[i].hist
             << "\", \"label\": \"" << kRegions[i].label << "\", \"current_integral\": "
             << currentHists[i]->Integral() << ", \"ppg12_integral\": " << ppg12Hists[i]->Integral()
             << ", \"integral_ratio\": " << currentHists[i]->Integral() / ppg12Hists[i]->Integral() << "}";
    manifest << (i == 3 ? "\n" : ",\n");
  }
  manifest << "  ],\n";
  manifest << "  \"normalization\": \"raw counts, no scale factor\",\n";
  manifest << "  \"binning\": \"validated identical 10,12,14,16,18,20,22,24,26,28,32,36 GeV edges\",\n";
  manifest << "  \"ratio_y_range\": [0.2, " << ratioYMax << "],\n";
  manifest << "  \"ratio_y_range_note\": \"upper bound auto-fit to max(Current/PPG12 + statistical error) with small visual buffer\",\n";
  manifest << "  \"max_abs_ratio_deviation\": " << maxAbsDev << ",\n";
  manifest << "  \"max_abs_ratio_deviation_where\": \"" << maxDevWhere << "\"\n";
  manifest << "}\n";
  manifest.close();

  printf("wrote_png=%s\n", pngPath.Data());
  printf("wrote_csv=%s\n", csvPath.Data());
  printf("wrote_manifest=%s\n", manifestPath.Data());
  printf("max_abs_ratio_deviation=%g at %s\n", maxAbsDev, maxDevWhere.c_str());
}
