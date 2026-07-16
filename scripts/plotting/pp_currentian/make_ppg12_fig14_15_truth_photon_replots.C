#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

namespace
{
constexpr int kBlackLine = 1;
constexpr int kRedLine = 2;
constexpr int kBlueLine = 4;
constexpr int kGreenLine = 8;

void set_style()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTextFont(42);
  gStyle->SetTitleSize(0.055, "XYZ");
  gStyle->SetLabelSize(0.047, "XYZ");
  gStyle->SetTitleOffset(1.15, "X");
  gStyle->SetTitleOffset(1.25, "Y");
  gStyle->SetNdivisions(510, "X");
  gStyle->SetNdivisions(510, "Y");
  TGaxis::SetMaxDigits(3);
}

void draw_sphenix(double x, double y, double size = 0.045)
{
  TLatex t;
  t.SetNDC();
  t.SetTextSize(size);
  t.SetTextFont(42);
  t.DrawLatex(x, y, "#bf{#it{sPHENIX}} Internal");
}

void draw_text(double x, double y, const char *text, double size = 0.04)
{
  TLatex t;
  t.SetNDC();
  t.SetTextSize(size);
  t.SetTextFont(42);
  t.DrawLatex(x, y, text);
}

void style_line(TH1 *h, Color_t color, Style_t marker, int width = 2)
{
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  h->SetMarkerSize(0.85);
  h->SetLineWidth(width);
}

void normalize(TH1 *h)
{
  const double integral = h->Integral();
  if (integral > 0.0)
  {
    h->Scale(1.0 / integral);
  }
}

std::string join_path(const std::string &a, const std::string &b)
{
  return a.empty() || a.back() == '/' ? a + b : a + "/" + b;
}

void ensure_dir(const std::string &path)
{
  gSystem->mkdir(path.c_str(), true);
}

void write_manifest(const std::string &outdir,
                    const std::string &timing_root,
                    const std::string &shape_root)
{
  std::ofstream out(join_path(outdir, "manifest.json"));
  out << "{\n";
  out << "  \"source\": \"PPG12 SDCC ROOT reproduction of IAN Fig.14/Fig.15 truth-photon timing panels\",\n";
  out << "  \"timing_root\": \"" << timing_root << "\",\n";
  out << "  \"shape_root\": \"" << shape_root << "\",\n";
  out << "  \"code_anchors\": [\n";
  out << "    \"ppg12codeGit/efficiencytool/plot_cluster_time.C\",\n";
  out << "    \"ppg12codeGit/plotting/plot_cluster_timing.C\",\n";
  out << "    \"ppg12codeGit/efficiencytool/ShowerShapeCheck.C\",\n";
  out << "    \"ppg12codeGit/plotting/plot_npb_time_bkgsub.C\"\n";
  out << "  ],\n";
  out << "  \"objects\": {\n";
  out << "    \"fig14_left\": [\"h_all_delta_t_mbd_vs_eta\", \"h_npb_delta_t_mbd_vs_eta\"],\n";
  out << "    \"fig14_middle\": \"h_npb_delta_t_mbd_vs_eta eta-slice ProjectionY, Rebin(2), unit-normalized\",\n";
  out << "    \"fig15_right\": \"h_npb_score_vs_time_eta0_pt0, RebinY(4), normalized per NPB-score row\"\n";
  out << "  }\n";
  out << "}\n";
}
} // namespace

void make_ppg12_fig14_15_truth_photon_replots(
    const char *timing_root =
        "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620/reference_roots/fig14_15_truth_photon/cluster_time_analysis_data.root",
    const char *shape_root =
        "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620/reference_roots/fig14_15_truth_photon/data_histoshower_shape_.root",
    const char *outdir =
        "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620/shower_shape_reference_validation/fig14_15_truth_photon")
{
  set_style();
  ensure_dir(outdir);

  TFile f_time(timing_root, "READ");
  TFile f_shape(shape_root, "READ");
  if (f_time.IsZombie() || f_shape.IsZombie())
  {
    Error("make_ppg12_fig14_15_truth_photon_replots", "Could not open input ROOT(s)");
    return;
  }

  auto *h_all = dynamic_cast<TH2D *>(f_time.Get("h_all_delta_t_mbd_vs_eta"));
  auto *h_npb = dynamic_cast<TH2D *>(f_time.Get("h_npb_delta_t_mbd_vs_eta"));
  auto *h_score = dynamic_cast<TH2D *>(f_shape.Get("h_npb_score_vs_time_eta0_pt0"));
  if (!h_all || !h_npb || !h_score)
  {
    Error("make_ppg12_fig14_15_truth_photon_replots", "Missing required histogram object");
    return;
  }

  // Fig. 14 left: cluster-MBD time vs eta with all-cluster and NPB mean profiles.
  {
    TCanvas c("c_fig14_left", "", 820, 620);
    c.SetLeftMargin(0.13);
    c.SetRightMargin(0.15);
    c.SetBottomMargin(0.13);
    c.SetTopMargin(0.04);
    c.SetLogz();

    auto *h_draw = dynamic_cast<TH2D *>(h_all->Clone("h_fig14_left_draw"));
    h_draw->SetDirectory(nullptr);
    h_draw->SetTitle("");
    h_draw->GetXaxis()->SetTitle("Cluster #eta");
    h_draw->GetYaxis()->SetTitle("Cluster - MBD Time [ns]");
    h_draw->GetYaxis()->SetRangeUser(-20.0, 20.0);
    h_draw->GetZaxis()->SetRangeUser(100.0, std::max(1000.0, 1.2 * h_draw->GetMaximum()));
    h_draw->Draw("COLZ");

    auto *prof_all = h_all->ProfileX("prof_fig14_all_mbd", 1, -1, "");
    auto *prof_npb = h_npb->ProfileX("prof_fig14_npb_mbd", 1, -1, "");
    style_line(prof_all, kRedLine, 20, 3);
    style_line(prof_npb, kBlueLine, 21, 3);
    prof_all->Draw("SAME");
    prof_npb->Draw("SAME");

    TLegend leg(0.16, 0.80, 0.54, 0.92);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.038);
    leg.AddEntry(prof_all, "All Clusters (Mean)", "lp");
    leg.AddEntry(prof_npb, "NPB Selection (Mean)", "lp");
    leg.Draw();

    draw_sphenix(0.18, 0.76, 0.04);
    draw_text(0.18, 0.70, "p_{T} > 10 GeV", 0.038);
    c.SaveAs(join_path(outdir, "ppg12_fig14_cluster_mbd_time_eta_reproduction.png").c_str());
    delete h_draw;
    delete prof_all;
    delete prof_npb;
  }

  // Fig. 14 middle: NPB cluster-MBD time projections in eta slices.
  {
    const std::vector<double> eta_edges = {-1.0, -0.35, 0.35, 0.7, 1.0};
    const std::vector<const char *> eta_labels = {
        "-1.0<#eta<-0.35",
        "-0.35<#eta<0.35",
        "0.35<#eta<0.7",
        "0.7<#eta<1.0"};
    const std::vector<Color_t> colors = {kBlackLine, kRed + 1, kBlue - 4, kGreen + 2};

    std::vector<TH1D *> slices;
    double ymax = 0.0;
    for (int islice = 0; islice < 4; ++islice)
    {
      h_npb->GetXaxis()->SetRangeUser(eta_edges[islice], eta_edges[islice + 1]);
      auto *proj = h_npb->ProjectionY(Form("h_fig14_eta_slice_%d", islice));
      proj->SetDirectory(nullptr);
      proj->Rebin(2);
      normalize(proj);
      style_line(proj, colors[islice], 20, 2);
      ymax = std::max(ymax, proj->GetMaximum());
      slices.push_back(proj);
    }
    h_npb->GetXaxis()->SetRange(0, 0);

    TCanvas c("c_fig14_middle", "", 650, 650);
    c.SetLeftMargin(0.15);
    c.SetRightMargin(0.04);
    c.SetBottomMargin(0.13);
    c.SetTopMargin(0.04);

    TH1D frame("h_fig14_middle_frame", "", 80, -20.0, 20.0);
    frame.SetStats(0);
    frame.GetXaxis()->SetTitle("Cluster - MBD Time [ns]");
    frame.GetYaxis()->SetTitle("Normalized Counts");
    frame.GetYaxis()->SetRangeUser(0.0, std::max(0.18, 1.2 * ymax));
    frame.Draw("AXIS");

    for (auto *h : slices)
    {
      h->Draw("HIST SAME");
    }

    draw_sphenix(0.20, 0.92, 0.04);
    draw_text(0.20, 0.84, "NPB selection", 0.04);

    TLegend leg(0.56, 0.70, 0.92, 0.91);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.034);
    for (int islice = 0; islice < 4; ++islice)
    {
      leg.AddEntry(slices[islice], eta_labels[islice], "l");
    }
    leg.Draw();

    c.SaveAs(join_path(outdir, "ppg12_fig14_npb_cluster_mbd_time_eta_slices_reproduction.png").c_str());
    for (auto *h : slices)
    {
      delete h;
    }
  }

  // Fig. 15: NPB score vs cluster-MBD time, normalized independently per NPB bin.
  {
    TCanvas c("c_fig15", "", 700, 610);
    c.SetLeftMargin(0.15);
    c.SetRightMargin(0.16);
    c.SetBottomMargin(0.14);
    c.SetTopMargin(0.04);

    auto *h_norm = dynamic_cast<TH2D *>(h_score->Clone("h_fig15_npb_score_vs_time_norm"));
    h_norm->SetDirectory(nullptr);
    h_norm->RebinY(4);
    for (int iy = 1; iy <= h_norm->GetNbinsY(); ++iy)
    {
      double row_sum = 0.0;
      for (int ix = 1; ix <= h_norm->GetNbinsX(); ++ix)
      {
        row_sum += h_norm->GetBinContent(ix, iy);
      }
      if (row_sum <= 0.0)
      {
        continue;
      }
      for (int ix = 1; ix <= h_norm->GetNbinsX(); ++ix)
      {
        h_norm->SetBinContent(ix, iy, h_norm->GetBinContent(ix, iy) / row_sum);
      }
    }
    h_norm->SetStats(0);
    h_norm->SetTitle("");
    h_norm->GetXaxis()->SetTitle("Cluster-MBD Time [ns]");
    h_norm->GetYaxis()->SetTitle("NPB Score");
    h_norm->GetXaxis()->SetRangeUser(-15.0, 10.0);
    h_norm->GetYaxis()->SetRangeUser(0.0, 1.0);
    h_norm->GetZaxis()->SetRangeUser(0.02, 0.20);
    h_norm->Draw("COLZ");
    draw_sphenix(0.18, 0.88, 0.04);
    draw_text(0.18, 0.81, "p+p #sqrt{s} = 200 GeV", 0.034);
    draw_text(0.18, 0.75, "|#eta^{#gamma}| < 0.7", 0.034);
    draw_text(0.18, 0.69, "10 < E_{T} < 14 GeV", 0.034);
    draw_text(0.18, 0.63, "Normalized per NPB bin", 0.034);
    c.SaveAs(join_path(outdir, "ppg12_fig15_npb_score_vs_cluster_mbd_time_reproduction.png").c_str());
    delete h_norm;
  }

  write_manifest(outdir, timing_root, shape_root);
}
