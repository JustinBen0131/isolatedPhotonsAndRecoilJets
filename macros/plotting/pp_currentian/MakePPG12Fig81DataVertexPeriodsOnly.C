#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace {

struct PeriodInput {
  const char *period;
  const char *label;
  const char *root_path;
  int color;
  double z_cut;
  int run_min;
  int run_max;
};

const char *kOutDir =
    "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620/"
    "shower_shape_reference_validation/fig81_vertex_period_data";

const char *kSourceGlob =
    "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/data/ana521/"
    "condorout/part_*_with_bdt_split.root";

void SetPlainStyle() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045, "XYZ");
  gStyle->SetTitleSize(0.055, "XYZ");
  gStyle->SetTitleOffset(1.05, "X");
  gStyle->SetTitleOffset(1.15, "Y");
}

std::unique_ptr<TH1> LoadHist(const PeriodInput &input) {
  auto file = std::make_unique<TFile>(input.root_path, "READ");
  if (!file || file->IsZombie() || file->TestBit(TFile::kRecovered)) {
    std::cerr << "ERROR: invalid ROOT file " << input.root_path << std::endl;
    return nullptr;
  }

  auto hist = dynamic_cast<TH1 *>(file->Get("h_D"));
  if (!hist) {
    std::cerr << "ERROR: missing h_D in " << input.root_path << std::endl;
    return nullptr;
  }

  auto clone = std::unique_ptr<TH1>(dynamic_cast<TH1 *>(hist->Clone(
      Form("h_D_%s", input.period))));
  clone->SetDirectory(nullptr);
  clone->SetLineColor(input.color);
  clone->SetMarkerColor(input.color);
  clone->SetLineWidth(3);
  clone->SetMarkerStyle(20);
  clone->SetMarkerSize(0.1);
  return clone;
}

void WriteCsv(const std::vector<PeriodInput> &inputs,
              const std::vector<std::unique_ptr<TH1>> &hists) {
  std::ofstream csv(std::string(kOutDir) +
                    "/ppg12_fig81_data_reco_vertex_periods_only.csv");
  csv << "period,run_min,run_max,bin,low,high,center,counts,error\n";
  csv << std::setprecision(12);
  for (size_t ih = 0; ih < hists.size(); ++ih) {
    const auto *h = hists[ih].get();
    const auto &p = inputs[ih];
    for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
      csv << p.period << "," << p.run_min << "," << p.run_max << "," << ib
          << "," << h->GetXaxis()->GetBinLowEdge(ib) << ","
          << h->GetXaxis()->GetBinUpEdge(ib) << ","
          << h->GetXaxis()->GetBinCenter(ib) << "," << h->GetBinContent(ib)
          << "," << h->GetBinError(ib) << "\n";
    }
  }
}

void WriteManifest(const std::vector<PeriodInput> &inputs,
                   const std::vector<std::unique_ptr<TH1>> &hists) {
  std::ofstream js(std::string(kOutDir) +
                   "/ppg12_fig81_data_reco_vertex_periods_only_manifest.json");
  js << std::setprecision(12);
  js << "{\n";
  js << "  \"artifact\": \"PPG12 Fig.81 data reco vertex period-only "
        "reproduction\",\n";
  js << "  \"scope\": \"data h_D only; iterative MC and SIM vertices omitted by "
        "request\",\n";
  js << "  \"source_data_glob_sdcc\": \"" << kSourceGlob << "\",\n";
  js << "  \"source_config\": "
        "\"ppg12codeGit/efficiencytool/truth_vertex_reweight/config.yaml\",\n";
  js << "  \"source_plot_code\": "
        "\"ppg12codeGit/efficiencytool/truth_vertex_reweight/"
        "plot_comparisons.py\",\n";
  js << "  \"trigger_bit\": 30,\n";
  js << "  \"outputs\": {\n";
  js << "    \"png\": \"" << kOutDir
     << "/ppg12_fig81_data_reco_vertex_periods_only.png\",\n";
  js << "    \"csv\": \"" << kOutDir
     << "/ppg12_fig81_data_reco_vertex_periods_only.csv\"\n";
  js << "  },\n";
  js << "  \"periods\": [\n";
  for (size_t i = 0; i < hists.size(); ++i) {
    const auto &p = inputs[i];
    const auto *h = hists[i].get();
    js << "    {\n";
    js << "      \"period\": \"" << p.period << "\",\n";
    js << "      \"label\": \"" << p.label << "\",\n";
    js << "      \"root_path\": \"" << p.root_path << "\",\n";
    js << "      \"histogram\": \"h_D\",\n";
    js << "      \"run_min\": " << p.run_min << ",\n";
    js << "      \"run_max\": " << p.run_max << ",\n";
    js << "      \"z_cut_cm\": " << p.z_cut << ",\n";
    js << "      \"integral_all_bins\": "
       << h->Integral(0, h->GetNbinsX() + 1) << ",\n";
    js << "      \"mean_cm\": " << h->GetMean() << ",\n";
    js << "      \"rms_cm\": " << h->GetRMS() << ",\n";
    js << "      \"nbins\": " << h->GetNbinsX() << ",\n";
    js << "      \"xmin_cm\": " << h->GetXaxis()->GetXmin() << ",\n";
    js << "      \"xmax_cm\": " << h->GetXaxis()->GetXmax() << "\n";
    js << "    }" << (i + 1 == hists.size() ? "\n" : ",\n");
  }
  js << "  ]\n";
  js << "}\n";
}

}  // namespace

void MakePPG12Fig81DataVertexPeriodsOnly() {
  SetPlainStyle();
  gSystem->mkdir(kOutDir, true);

  std::vector<PeriodInput> inputs = {
      {"1p5mrad", "1.5 mrad data",
       "ppg12codeGit/efficiencytool/truth_vertex_reweight/output/1p5mrad/"
       "reweight.root.photon10_backup",
       kBlue + 1, 83.0, 51274, 54000},
      {"0mrad", "0 mrad data",
       "ppg12codeGit/efficiencytool/truth_vertex_reweight/output/0mrad/"
       "reweight.root.photon10_backup",
       kRed + 1, 144.0, 47289, 51274},
  };

  std::vector<std::unique_ptr<TH1>> hists;
  for (const auto &input : inputs) {
    auto h = LoadHist(input);
    if (!h) {
      return;
    }
    hists.push_back(std::move(h));
  }

  auto c = std::make_unique<TCanvas>("c_fig81_data_vertex_periods",
                                     "PPG12 Fig81 data vertex periods", 1100,
                                     760);
  c->SetLeftMargin(0.13);
  c->SetRightMargin(0.04);
  c->SetTopMargin(0.12);
  c->SetBottomMargin(0.13);
  c->SetLogy();
  c->SetGrid(1, 1);

  auto frame = std::unique_ptr<TH1>(dynamic_cast<TH1 *>(hists[0]->Clone(
      "frame_fig81_data_vertex_periods")));
  frame->Reset("ICESM");
  frame->SetMinimum(2.0e2);
  frame->SetMaximum(5.0e6);
  frame->GetXaxis()->SetRangeUser(-200.0, 200.0);
  frame->GetXaxis()->SetTitle("z_{reco} (cm)");
  frame->GetYaxis()->SetTitle("events");
  frame->GetXaxis()->SetTitleSize(0.055);
  frame->GetYaxis()->SetTitleSize(0.055);
  frame->GetXaxis()->SetLabelSize(0.045);
  frame->GetYaxis()->SetLabelSize(0.045);
  frame->Draw("AXIS");

  for (const auto &h : hists) {
    h->Draw("HIST SAME");
  }

  for (const auto &p : inputs) {
    for (double x : {-p.z_cut, p.z_cut}) {
      auto line = new TLine(x, 2.0e2, x, 5.0e6);
      line->SetLineColor(p.color);
      line->SetLineStyle(3);
      line->SetLineWidth(1);
      line->Draw();
    }
  }

  TLatex tex;
  tex.SetNDC();
  tex.SetTextFont(42);
  tex.SetTextSize(0.033);
  tex.DrawLatex(0.16, 0.955, "Reco vertex: PPG12 data periods");
  tex.SetTextSize(0.030);
  tex.DrawLatex(0.16, 0.855, "#it{sPHENIX} Internal");
  tex.DrawLatex(0.16, 0.800, "p+p #sqrt{s} = 200 GeV");
  tex.DrawLatex(0.16, 0.750, "Photon 4 GeV trigger, data only");

  auto leg = std::make_unique<TLegend>(0.43, 0.13, 0.70, 0.28);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.028);
  for (size_t i = 0; i < inputs.size(); ++i) {
    leg->AddEntry(hists[i].get(), inputs[i].label, "l");
  }
  leg->Draw();

  c->RedrawAxis();
  const std::string png =
      std::string(kOutDir) + "/ppg12_fig81_data_reco_vertex_periods_only.png";
  c->SaveAs(png.c_str());

  WriteCsv(inputs, hists);
  WriteManifest(inputs, hists);

  std::cout << "wrote " << png << std::endl;
}
