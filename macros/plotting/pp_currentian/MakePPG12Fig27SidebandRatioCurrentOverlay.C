#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TPad.h>
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
struct RatioSpec
{
  const char *key;
  const char *numHist;
  const char *denHist;
  const char *label;
  int color;
  int marker;
  double offset;
};

const RatioSpec kRatios[] = {
    {"BoverA", "h_tight_noniso_cluster_0", "h_tight_iso_cluster_0",
     "B/A: tight noniso", kBlack, 20, -0.18},
    {"CoverA", "h_nontight_iso_cluster_0", "h_tight_iso_cluster_0",
     "C/A: nontight iso", kRed + 1, 21, 0.0},
    {"DoverA", "h_nontight_noniso_cluster_0", "h_tight_iso_cluster_0",
     "D/A: nontight noniso", kBlue + 1, 22, 0.18},
};

TH1 *requireHist(TDirectory *dir, const char *name, const char *label)
{
  if (!dir)
  {
    Error("MakePPG12Fig27SidebandRatioCurrentOverlay", "Null directory for %s", label);
    return nullptr;
  }
  auto *hist = dynamic_cast<TH1 *>(dir->Get(name));
  if (!hist)
  {
    Error("MakePPG12Fig27SidebandRatioCurrentOverlay", "Missing %s histogram %s", label, name);
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

TH1 *makeRootDivideRatio(const TH1 *num, const TH1 *den, const RatioSpec &spec, const char *prefix)
{
  auto *ratio = dynamic_cast<TH1 *>(num->Clone(Form("%s_%s", prefix, spec.key)));
  ratio->SetDirectory(nullptr);
  ratio->SetStats(false);
  ratio->Divide(den);
  ratio->SetLineColor(spec.color);
  ratio->SetMarkerColor(spec.color);
  ratio->SetLineWidth(2);
  return ratio;
}

TGraphErrors *histToGraph(const TH1 *hist, const RatioSpec &spec, bool current)
{
  const int n = hist->GetNbinsX();
  auto *graph = new TGraphErrors(n);
  graph->SetName(Form("g_%s_%s", current ? "current" : "ppg12", spec.key));
  for (int i = 1; i <= n; ++i)
  {
    const double width = hist->GetXaxis()->GetBinWidth(i);
    const double x = hist->GetXaxis()->GetBinCenter(i) + spec.offset * width;
    graph->SetPoint(i - 1, x, hist->GetBinContent(i));
    graph->SetPointError(i - 1, 0.0, hist->GetBinError(i));
  }
  graph->SetLineColor(spec.color);
  graph->SetMarkerColor(spec.color);
  graph->SetMarkerStyle(current ? spec.marker : 24);
  graph->SetMarkerSize(current ? 1.02 : 1.15);
  graph->SetLineWidth(2);
  return graph;
}

TGraphErrors *makeRatioOfRatiosGraph(const TH1 *current, const TH1 *ppg12, const RatioSpec &spec)
{
  const int n = current->GetNbinsX();
  auto *graph = new TGraphErrors(n);
  graph->SetName(Form("g_current_over_ppg12_%s", spec.key));
  for (int i = 1; i <= n; ++i)
  {
    const double width = current->GetXaxis()->GetBinWidth(i);
    const double x = current->GetXaxis()->GetBinCenter(i) + spec.offset * width;
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
  graph->SetLineColor(spec.color);
  graph->SetMarkerColor(spec.color);
  graph->SetMarkerStyle(spec.marker);
  graph->SetMarkerSize(0.95);
  graph->SetLineWidth(2);
  return graph;
}

void setupTopFrame(TH1 *frame)
{
  frame->SetStats(false);
  frame->GetXaxis()->SetRangeUser(10.0, 36.0);
  frame->GetYaxis()->SetRangeUser(0.0, 2.0);
  frame->GetYaxis()->SetTitle("ABCD yield ratio");
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetXaxis()->SetTitleSize(0);
  frame->GetYaxis()->SetTitleSize(0.062);
  frame->GetYaxis()->SetTitleOffset(0.88);
  frame->GetYaxis()->SetLabelSize(0.052);
  frame->GetYaxis()->SetNdivisions(505);
}

std::pair<double, double> ratioOfRatiosYRange(const std::vector<TGraphErrors *> &graphs)
{
  double ymin = 0.45;
  double ymax = 1.75;
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
      const double lower = y - graph->GetErrorY(i);
      const double upper = y + graph->GetErrorY(i);
      if (std::isfinite(lower))
        ymin = std::min(ymin, lower);
      if (std::isfinite(upper))
        ymax = std::max(ymax, upper);
    }
  }
  const double span = std::max(0.1, ymax - ymin);
  const double buffer = std::max(0.06, 0.08 * span);
  ymin = std::max(0.0, std::floor((ymin - buffer) * 20.0) / 20.0);
  ymax = std::ceil((ymax + buffer) * 20.0) / 20.0;
  return {ymin, ymax};
}

void setupBottomFrame(TH1 *frame, double ymin, double ymax)
{
  frame->SetStats(false);
  frame->GetXaxis()->SetRangeUser(10.0, 36.0);
  frame->GetYaxis()->SetRangeUser(ymin, ymax);
  frame->GetXaxis()->SetTitle("E_{T}^{#gamma,reco} [GeV]");
  frame->GetYaxis()->SetTitle("Current / PPG12");
  frame->GetXaxis()->SetTitleSize(0.13);
  frame->GetXaxis()->SetTitleOffset(0.92);
  frame->GetXaxis()->SetLabelSize(0.105);
  frame->GetYaxis()->SetTitleSize(0.09);
  frame->GetYaxis()->SetTitleOffset(0.57);
  frame->GetYaxis()->SetLabelSize(0.083);
  frame->GetYaxis()->SetNdivisions(505);
}

std::string fmt(double value, int precision = 6)
{
  std::ostringstream os;
  os << std::setprecision(precision) << std::defaultfloat << value;
  return os.str();
}
}  // namespace

void MakePPG12Fig27SidebandRatioCurrentOverlay(
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
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLineWidth(2);
  gStyle->SetEndErrorSize(3);

  TFile fCurrent(current_root, "READ");
  TFile fPpg12(ppg12_root, "READ");
  if (fCurrent.IsZombie() || fCurrent.TestBit(TFile::kRecovered))
  {
    Error("MakePPG12Fig27SidebandRatioCurrentOverlay",
          "Current ROOT is zombie or recovered: %s", current_root);
    return;
  }
  if (fPpg12.IsZombie() || fPpg12.TestBit(TFile::kRecovered))
  {
    Error("MakePPG12Fig27SidebandRatioCurrentOverlay",
          "PPG12 ROOT is zombie or recovered: %s", ppg12_root);
    return;
  }

  auto *currentDir = dynamic_cast<TDirectory *>(fCurrent.Get(current_dir_name));
  if (!currentDir)
  {
    Error("MakePPG12Fig27SidebandRatioCurrentOverlay", "Missing current directory %s", current_dir_name);
    return;
  }

  std::vector<TH1 *> currentRatios;
  std::vector<TH1 *> ppg12Ratios;
  std::vector<TGraphErrors *> currentGraphs;
  std::vector<TGraphErrors *> ppg12Graphs;
  std::vector<TGraphErrors *> ratioOfRatiosGraphs;

  for (const auto &spec : kRatios)
  {
    TH1 *cNum = requireHist(currentDir, spec.numHist, Form("current_num_%s", spec.key));
    TH1 *cDen = requireHist(currentDir, spec.denHist, Form("current_den_%s", spec.key));
    TH1 *pNum = requireHist(&fPpg12, spec.numHist, Form("ppg12_num_%s", spec.key));
    TH1 *pDen = requireHist(&fPpg12, spec.denHist, Form("ppg12_den_%s", spec.key));
    if (!cNum || !cDen || !pNum || !pDen)
      return;
    if (!sameBinning(cNum, cDen) || !sameBinning(cNum, pNum) || !sameBinning(cNum, pDen))
    {
      Error("MakePPG12Fig27SidebandRatioCurrentOverlay", "Binning mismatch for %s", spec.key);
      return;
    }
    TH1 *cRatio = makeRootDivideRatio(cNum, cDen, spec, "current");
    TH1 *pRatio = makeRootDivideRatio(pNum, pDen, spec, "ppg12");
    currentRatios.push_back(cRatio);
    ppg12Ratios.push_back(pRatio);
    currentGraphs.push_back(histToGraph(cRatio, spec, true));
    ppg12Graphs.push_back(histToGraph(pRatio, spec, false));
    ratioOfRatiosGraphs.push_back(makeRatioOfRatiosGraph(cRatio, pRatio, spec));
  }

  auto *canvas = new TCanvas("c_fig27_sideband_ratio_current_vs_ppg12",
                             "c_fig27_sideband_ratio_current_vs_ppg12", 900, 1050);
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
  auto *frameTop = dynamic_cast<TH1 *>(ppg12Ratios[0]->Clone("frameTop"));
  frameTop->Reset("ICES");
  setupTopFrame(frameTop);
  frameTop->Draw("AXIS");
  for (auto *g : currentGraphs)
    g->Draw("PZ SAME");
  for (auto *g : ppg12Graphs)
    g->Draw("PZ SAME");

  TLatex text;
  text.SetNDC(true);
  text.SetTextFont(42);
  text.SetTextSize(0.039);
  text.DrawLatex(0.19, 0.91, "#bf{sPHENIX} Internal");
  text.DrawLatex(0.19, 0.855, "#it{p}+#it{p} #sqrt{s} = 200 GeV");
  text.DrawLatex(0.19, 0.800, "|#eta^{#gamma}| < 0.7");
  text.DrawLatex(0.19, 0.745, "Data, bdt_nom");
  text.DrawLatex(0.19, 0.690, "PPG12 Fig. 27 sideband ratios");

  auto *legSource = new TLegend(0.55, 0.79, 0.84, 0.93);
  legSource->SetBorderSize(0);
  legSource->SetFillStyle(0);
  legSource->SetTextFont(42);
  legSource->SetTextSize(0.037);
  legSource->SetHeader("source", "C");
  legSource->AddEntry(ppg12Graphs[0], "PPG12 SDCC", "p");
  legSource->AddEntry(currentGraphs[0], current_label, "p");
  legSource->Draw();

  auto *legRatio = new TLegend(0.55, 0.56, 0.91, 0.76);
  legRatio->SetBorderSize(0);
  legRatio->SetFillStyle(0);
  legRatio->SetTextFont(42);
  legRatio->SetTextSize(0.037);
  legRatio->SetHeader("ratio definition", "C");
  for (int i = 0; i < 3; ++i)
    legRatio->AddEntry(currentGraphs[i], kRatios[i].label, "p");
  legRatio->Draw();

  padBot->cd();
  auto *frameBottom = dynamic_cast<TH1 *>(ppg12Ratios[0]->Clone("frameBottom"));
  frameBottom->Reset("ICES");
  const auto bottomRange = ratioOfRatiosYRange(ratioOfRatiosGraphs);
  setupBottomFrame(frameBottom, bottomRange.first, bottomRange.second);
  frameBottom->Draw("AXIS");
  TLine one(10.0, 1.0, 36.0, 1.0);
  one.SetLineColor(kGray + 2);
  one.SetLineStyle(2);
  one.SetLineWidth(2);
  one.Draw();
  for (auto *g : ratioOfRatiosGraphs)
    g->Draw("PZ SAME");

  const std::string png = std::string(output_dir) + "/fig27_sideband_ratios_current_vs_ppg12_overlay_ratio.png";
  canvas->SaveAs(png.c_str());

  const std::string csv = std::string(output_dir) + "/fig27_sideband_ratios_current_vs_ppg12_overlay_ratio.csv";
  std::ofstream out(csv);
  out << "ratio,bin,low,high,center,ppg12,ppg12_err,current,current_err,current_over_ppg12,current_over_ppg12_err\n";
  double maxAbsDev = -1.0;
  std::string maxDesc;
  for (int ir = 0; ir < 3; ++ir)
  {
    const TH1 *p = ppg12Ratios[ir];
    const TH1 *c = currentRatios[ir];
    for (int i = 1; i <= p->GetNbinsX(); ++i)
    {
      const double pv = p->GetBinContent(i);
      const double cv = c->GetBinContent(i);
      const double pe = p->GetBinError(i);
      const double ce = c->GetBinError(i);
      double rr = std::numeric_limits<double>::quiet_NaN();
      double re = 0.0;
      if (pv > 0)
      {
        rr = cv / pv;
        re = std::sqrt((ce / pv) * (ce / pv) + (cv * pe / (pv * pv)) * (cv * pe / (pv * pv)));
        const double dev = std::fabs(rr - 1.0);
        if (dev > maxAbsDev)
        {
          maxAbsDev = dev;
          maxDesc = std::string(kRatios[ir].key) + " " + fmt(p->GetXaxis()->GetBinLowEdge(i), 3) + "-" +
                    fmt(p->GetXaxis()->GetBinUpEdge(i), 3) + " GeV";
        }
      }
      out << kRatios[ir].key << "," << i << "," << p->GetXaxis()->GetBinLowEdge(i) << ","
          << p->GetXaxis()->GetBinUpEdge(i) << "," << p->GetXaxis()->GetBinCenter(i) << ","
          << pv << "," << pe << "," << cv << "," << ce << "," << rr << "," << re << "\n";
    }
  }

  const std::string manifest = std::string(output_dir) + "/fig27_sideband_ratios_current_vs_ppg12_overlay_ratio_manifest.json";
  std::ofstream meta(manifest);
  meta << "{\n";
  meta << "  \"plot\": \"" << png << "\",\n";
  meta << "  \"csv\": \"" << csv << "\",\n";
  meta << "  \"ppg12_root\": \"" << ppg12_root << "\",\n";
  meta << "  \"current_root\": \"" << current_root << "\",\n";
  meta << "  \"current_directory\": \"" << current_dir_name << "\",\n";
  meta << "  \"current_label\": \"" << current_label << "\",\n";
  meta << "  \"source_code_reference\": \"ppg12codeGit/plotting/plot_sideband_selection.C\",\n";
  meta << "  \"ratio_definitions\": [\"B/A: h_tight_noniso_cluster_0 / h_tight_iso_cluster_0\", "
       << "\"C/A: h_nontight_iso_cluster_0 / h_tight_iso_cluster_0\", "
       << "\"D/A: h_nontight_noniso_cluster_0 / h_tight_iso_cluster_0\"],\n";
  meta << "  \"normalization\": \"none; direct ROOT TH1::Divide sideband-yield ratios\",\n";
  meta << "  \"bottom_panel\": \"current ratio divided by PPG12 ratio, propagated from ROOT ratio errors\",\n";
  meta << "  \"bottom_y_range\": [" << bottomRange.first << ", " << bottomRange.second << "],\n";
  meta << "  \"max_abs_current_over_ppg12_minus_one\": " << maxAbsDev << ",\n";
  meta << "  \"max_abs_deviation_bin\": \"" << maxDesc << "\"\n";
  meta << "}\n";

  Info("MakePPG12Fig27SidebandRatioCurrentOverlay", "Wrote %s", png.c_str());
}
