#include <TCanvas.h>
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
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace {

struct PeriodConfig {
  const char *period;
  const char *title;
  const char *ppg12_root;
  const char *ppg12_hist;
  const char *current_hist;
  int color;
  double z_cut;
};

const char *kCurrentRoot =
    "dataOutput/ppg12Parity/the76_ppg12_parity_full_20260701_003024/"
    "final_roots/pp_data_hierarchical_v1/RecoilJets_pp_ALL_PERIOD_COMBINED.root";

const char *kOutDir =
    "dataOutput/ppg12Parity/the76_ppg12_parity_full_20260701_003024/"
    "data_vertex_fig81";

const char *kPpg12DataGlob =
    "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/data/ana521/"
    "condorout/part_*_with_bdt_split.root";

void SetStyle() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetEndErrorSize(0);
}

std::unique_ptr<TH1> LoadHist(TFile &file, const char *path,
                              const char *clone_name) {
  auto *obj = file.Get(path);
  auto *hist = dynamic_cast<TH1 *>(obj);
  if (!hist) {
    std::cerr << "ERROR: missing histogram " << path << std::endl;
    return nullptr;
  }
  auto out = std::unique_ptr<TH1>(dynamic_cast<TH1 *>(hist->Clone(clone_name)));
  out->SetDirectory(nullptr);
  out->Sumw2(false);
  out->Sumw2(true);
  return out;
}

std::unique_ptr<TH1> LoadPpg12(const PeriodConfig &cfg) {
  TFile file(cfg.ppg12_root, "READ");
  if (file.IsZombie() || file.TestBit(TFile::kRecovered)) {
    std::cerr << "ERROR: invalid PPG12 ROOT " << cfg.ppg12_root << std::endl;
    return nullptr;
  }
  return LoadHist(file, cfg.ppg12_hist, Form("ppg12_%s", cfg.period));
}

std::unique_ptr<TH1> LoadCurrent(TFile &file, const PeriodConfig &cfg) {
  auto hist = LoadHist(file, cfg.current_hist, Form("current_%s_1cm", cfg.period));
  if (!hist) return nullptr;
  if (hist->GetNbinsX() != 400) {
    std::cerr << "ERROR: expected 400 current bins, got " << hist->GetNbinsX()
              << " for " << cfg.current_hist << std::endl;
    return nullptr;
  }
  auto *rebinned_raw = hist->Rebin(5, Form("current_%s", cfg.period));
  auto out = std::unique_ptr<TH1>(dynamic_cast<TH1 *>(rebinned_raw));
  out->SetDirectory(nullptr);
  return out;
}

void NormalizeShape(TH1 *hist) {
  const double integral = hist->Integral(1, hist->GetNbinsX());
  if (integral > 0.0) hist->Scale(1.0 / integral);
}

std::unique_ptr<TGraphErrors> MakeRatioGraph(const TH1 *num, const TH1 *den,
                                             const char *name) {
  auto graph = std::make_unique<TGraphErrors>();
  graph->SetName(name);
  int ip = 0;
  for (int ib = 1; ib <= den->GetNbinsX(); ++ib) {
    const double d = den->GetBinContent(ib);
    const double n = num->GetBinContent(ib);
    if (d <= 0.0 || n <= 0.0) continue;
    const double ratio = n / d;
    const double en = num->GetBinError(ib);
    const double ed = den->GetBinError(ib);
    const double eratio =
        ratio * std::sqrt(std::pow(en / n, 2) + std::pow(ed / d, 2));
    graph->SetPoint(ip, den->GetXaxis()->GetBinCenter(ib), ratio);
    graph->SetPointError(ip, 0.0, eratio);
    ++ip;
  }
  return graph;
}

void ConfigureAxis(TH1 *hist, double label_scale = 1.0) {
  hist->GetXaxis()->SetTitleSize(0.075 * label_scale);
  hist->GetYaxis()->SetTitleSize(0.072 * label_scale);
  hist->GetXaxis()->SetLabelSize(0.060 * label_scale);
  hist->GetYaxis()->SetLabelSize(0.060 * label_scale);
  hist->GetYaxis()->SetTitleOffset(0.78 / label_scale);
  hist->GetXaxis()->SetTitleOffset(0.93);
}

void DrawPeriod(TPad *top, TPad *bottom, const PeriodConfig &cfg, TH1 *ppg12,
                TH1 *current, const char *currentRoot, const std::string &stage,
                std::ofstream &csv,
                std::ofstream &manifest) {
  auto ppg12_shape = std::unique_ptr<TH1>(
      dynamic_cast<TH1 *>(ppg12->Clone(Form("ppg12_shape_%s", cfg.period))));
  auto current_shape = std::unique_ptr<TH1>(dynamic_cast<TH1 *>(
      current->Clone(Form("current_shape_%s", cfg.period))));
  ppg12_shape->SetDirectory(nullptr);
  current_shape->SetDirectory(nullptr);
  NormalizeShape(ppg12_shape.get());
  NormalizeShape(current_shape.get());

  ppg12_shape->SetLineColor(kBlack);
  ppg12_shape->SetMarkerColor(kBlack);
  ppg12_shape->SetMarkerStyle(24);
  ppg12_shape->SetMarkerSize(0.45);
  ppg12_shape->SetLineWidth(2);
  current_shape->SetLineColor(cfg.color);
  current_shape->SetMarkerColor(cfg.color);
  current_shape->SetMarkerStyle(20);
  current_shape->SetMarkerSize(0.42);
  current_shape->SetLineWidth(2);

  top->cd();
  top->SetLogy();
  top->SetGrid(1, 1);
  top->SetBottomMargin(0.02);
  top->SetTopMargin(0.08);
  top->SetLeftMargin(0.13);
  top->SetRightMargin(0.04);

  auto frame = std::unique_ptr<TH1>(
      dynamic_cast<TH1 *>(ppg12_shape->Clone(Form("frame_%s", cfg.period))));
  frame->Reset("ICESM");
  frame->SetMinimum(1e-6);
  frame->SetMaximum(0.30);
  frame->GetXaxis()->SetRangeUser(-200.0, 200.0);
  frame->GetYaxis()->SetTitle("shape-normalized events");
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetXaxis()->SetTitleSize(0);
  frame->GetYaxis()->SetTitleSize(0.048);
  frame->GetYaxis()->SetLabelSize(0.044);
  frame->GetYaxis()->SetTitleOffset(1.15);
  frame->Draw("AXIS");
  ppg12_shape->Draw("E1 SAME");
  current_shape->Draw("E1 SAME");

  for (double x : {-cfg.z_cut, cfg.z_cut}) {
    auto *line = new TLine(x, 1e-6, x, 0.30);
    line->SetLineColor(cfg.color);
    line->SetLineStyle(3);
    line->SetLineWidth(1);
    line->Draw();
  }

  TLatex tex;
  tex.SetNDC();
  tex.SetTextFont(42);
  tex.SetTextSize(0.042);
  tex.DrawLatex(0.18, 0.86, cfg.title);
  tex.SetTextSize(0.032);
  tex.DrawLatex(0.18, 0.78, "#it{sPHENIX} Internal");
  tex.DrawLatex(0.18, 0.715, "p+p #sqrt{s} = 200 GeV");
  tex.DrawLatex(0.18, 0.655,
                stage == "pre_vzcut" ? "trigger bit 30, pre-vzcut"
                                      : "trigger bit 30, post-vzcut");

  auto *leg = new TLegend(0.47, 0.72, 0.91, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.031);
  leg->AddEntry(ppg12_shape.get(), "PPG12 SDCC h_{D}", "lep");
  leg->AddEntry(current_shape.get(), "This analysis", "lep");
  leg->Draw();

  bottom->cd();
  bottom->SetGrid(1, 1);
  bottom->SetTopMargin(0.04);
  bottom->SetBottomMargin(0.30);
  bottom->SetLeftMargin(0.13);
  bottom->SetRightMargin(0.04);

  auto ratio_frame = std::unique_ptr<TH1>(dynamic_cast<TH1 *>(
      ppg12_shape->Clone(Form("ratio_frame_%s", cfg.period))));
  ratio_frame->Reset("ICESM");
  ratio_frame->SetMinimum(0.40);
  ratio_frame->SetMaximum(6.20);
  ratio_frame->GetXaxis()->SetRangeUser(-200.0, 200.0);
  ratio_frame->GetXaxis()->SetTitle("z_{reco} (cm)");
  ratio_frame->GetYaxis()->SetTitle("Current / PPG12");
  ConfigureAxis(ratio_frame.get(), 0.84);
  ratio_frame->GetYaxis()->SetNdivisions(505);
  ratio_frame->Draw("AXIS");

  auto *unit = new TLine(-200.0, 1.0, 200.0, 1.0);
  unit->SetLineStyle(2);
  unit->SetLineColor(kGray + 2);
  unit->Draw();

  auto ratio = MakeRatioGraph(current_shape.get(), ppg12_shape.get(),
                              Form("ratio_%s", cfg.period));
  ratio->SetMarkerStyle(20);
  ratio->SetMarkerColor(cfg.color);
  ratio->SetLineColor(cfg.color);
  ratio->SetMarkerSize(0.45);
  ratio->Draw("P SAME");

  for (int ib = 1; ib <= ppg12_shape->GetNbinsX(); ++ib) {
    const double p = ppg12_shape->GetBinContent(ib);
    const double c = current_shape->GetBinContent(ib);
    const double r = (p > 0.0 && c > 0.0) ? c / p : std::numeric_limits<double>::quiet_NaN();
    csv << cfg.period << "," << ib << ","
        << ppg12_shape->GetXaxis()->GetBinLowEdge(ib) << ","
        << ppg12_shape->GetXaxis()->GetBinUpEdge(ib) << ","
        << ppg12_shape->GetXaxis()->GetBinCenter(ib) << ","
        << ppg12->GetBinContent(ib) << "," << current->GetBinContent(ib)
        << "," << p << "," << c << "," << r << "\n";
  }

  manifest << "    {\n";
  manifest << "      \"period\": \"" << cfg.period << "\",\n";
  manifest << "      \"ppg12_root\": \"" << cfg.ppg12_root << "\",\n";
  manifest << "      \"ppg12_hist\": \"" << cfg.ppg12_hist << "\",\n";
  manifest << "      \"current_root\": \"" << currentRoot << "\",\n";
  manifest << "      \"current_hist\": \"" << cfg.current_hist << "\",\n";
  manifest << "      \"current_rebin_factor\": 5,\n";
  manifest << "      \"ppg12_raw_integral\": " << ppg12->Integral(0, ppg12->GetNbinsX() + 1) << ",\n";
  manifest << "      \"current_raw_integral\": " << current->Integral(0, current->GetNbinsX() + 1) << ",\n";
  manifest << "      \"z_cut_cm\": " << cfg.z_cut << "\n";
  manifest << "    }";

  // Keep drawn ROOT objects alive until the canvas is printed.
  frame.release();
  ppg12_shape.release();
  current_shape.release();
  ratio_frame.release();
  ratio.release();
}

}  // namespace

void MakePPG12Fig81DataVertexCurrentOverlay() {
  SetStyle();
  const char *currentRoot = gSystem->Getenv("RJ_FIG81_CURRENT_ROOT");
  if (!currentRoot || std::string(currentRoot).empty()) currentRoot = kCurrentRoot;
  const char *outDir = gSystem->Getenv("RJ_FIG81_OUT_DIR");
  if (!outDir || std::string(outDir).empty()) outDir = kOutDir;
  const char *stageEnv = gSystem->Getenv("RJ_FIG81_STAGE");
  std::string stage = (stageEnv && std::string(stageEnv).size()) ? stageEnv : "pre_vzcut";
  if (stage != "pre_vzcut" && stage != "post_vzcut") {
    std::cerr << "ERROR: RJ_FIG81_STAGE must be pre_vzcut or post_vzcut, got "
              << stage << std::endl;
    return;
  }
  gSystem->mkdir(outDir, true);

  std::vector<PeriodConfig> configs = {
      {"1p5mrad", "1.5 mrad reco vertex",
       "ppg12codeGit/efficiencytool/truth_vertex_reweight/output/1p5mrad/"
       "reweight.root.photon10_backup",
       "h_D",
       "",
       kBlue + 1, 83.0},
      {"0mrad", "0 mrad reco vertex",
       "ppg12codeGit/efficiencytool/truth_vertex_reweight/output/0mrad/"
       "reweight.root.photon10_backup",
       "h_D",
       "",
       kRed + 1, 144.0},
  };
  std::vector<std::string> currentHistPaths;
  currentHistPaths.reserve(configs.size());
  for (auto &cfg : configs) {
    currentHistPaths.emplace_back(
        std::string("PPG12_scaledtrigger30/h_ppg12_vtxqa_data_reco_z_triggered_") +
        cfg.period + "_" + stage);
    cfg.current_hist = currentHistPaths.back().c_str();
  }

  TFile current_file(currentRoot, "READ");
  if (current_file.IsZombie() || current_file.TestBit(TFile::kRecovered)) {
    std::cerr << "ERROR: invalid current ROOT " << currentRoot << std::endl;
    return;
  }

  std::ofstream csv(std::string(outDir) +
                    "/fig81_data_reco_vertex_sdcc_vs_current_shape_ratio.csv");
  csv << std::setprecision(12);
  csv << "period,bin,low,high,center,ppg12_raw,current_raw_rebinned,"
         "ppg12_shape,current_shape,current_over_ppg12_shape\n";

  std::ofstream manifest(std::string(outDir) +
                         "/fig81_data_reco_vertex_sdcc_vs_current_shape_ratio_manifest.json");
  manifest << std::setprecision(12);
  manifest << "{\n";
  manifest << "  \"artifact\": \"PPG12 Fig.81 data reco vertex current overlay\",\n";
  manifest << "  \"ratio_definition\": \"shape-normalized current / shape-normalized PPG12 SDCC\",\n";
  manifest << "  \"source_data_glob_sdcc\": \"" << kPpg12DataGlob << "\",\n";
  manifest << "  \"source_config\": \"ppg12codeGit/efficiencytool/truth_vertex_reweight/config.yaml\",\n";
  manifest << "  \"source_plot_code\": \"ppg12codeGit/efficiencytool/truth_vertex_reweight/plot_comparisons.py\",\n";
  manifest << "  \"current_note\": \"Current RecoilJets histogram is trigger-bit-30 period QA at current_stage; it is rebinned from 1 cm to PPG12 5 cm bins and shape-normalized before ratio.\",\n";
  manifest << "  \"current_root_override\": \"" << currentRoot << "\",\n";
  manifest << "  \"current_stage\": \"" << stage << "\",\n";
  manifest << "  \"periods\": [\n";

  auto canvas = std::make_unique<TCanvas>("c_fig81_vertex_overlay",
                                          "Fig81 vertex current overlay", 1500,
                                          850);

  for (size_t i = 0; i < configs.size(); ++i) {
    const double x0 = i == 0 ? 0.00 : 0.50;
    const double x1 = i == 0 ? 0.50 : 1.00;
    canvas->cd();
    auto *top = new TPad(Form("top_%zu", i), "", x0, 0.32, x1, 1.00);
    auto *bot = new TPad(Form("bot_%zu", i), "", x0, 0.00, x1, 0.32);
    top->Draw();
    bot->Draw();

    auto ppg12 = LoadPpg12(configs[i]);
    auto current = LoadCurrent(current_file, configs[i]);
    if (!ppg12 || !current) return;
    DrawPeriod(top, bot, configs[i], ppg12.get(), current.get(), currentRoot, stage, csv, manifest);
    manifest << (i + 1 == configs.size() ? "\n" : ",\n");
  }
  manifest << "  ],\n";
  manifest << "  \"outputs\": {\n";
  manifest << "    \"png\": \"" << outDir << "/fig81_data_reco_vertex_sdcc_vs_current_shape_ratio.png\",\n";
  manifest << "    \"csv\": \"" << outDir << "/fig81_data_reco_vertex_sdcc_vs_current_shape_ratio.csv\"\n";
  manifest << "  }\n";
  manifest << "}\n";

  const std::string png = std::string(outDir) +
                          "/fig81_data_reco_vertex_sdcc_vs_current_shape_ratio.png";
  canvas->SaveAs(png.c_str());
  std::cout << "wrote " << png << std::endl;
}
