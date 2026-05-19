#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TStyle.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace
{
struct Point
{
  double x = 0;
  double y = 0;
  double ey = 0;
};

const std::string kRepo = "/Users/patsfan753/Desktop/ThesisAnalysis";
const std::string kOutDir = kRepo + "/dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439/slideReady/basev3e_reco_efficiency_overlay";
const std::string kReferenceCsv = kRepo + "/dataOutput/target80_first_offline/bdt_target80_gated_20260512_001012/analysis_config_etfine_15to35_target80/efficiency_qa/truth_photon_recovery_box_cuts_points.csv";
const std::string kBaseCsv = kOutDir + "/basev3e_target80_reco_efficiency_pt15to35_points.csv";
const std::string kGlobalBdtCsv = kOutDir + "/global32_cent3grid_target80_reco_efficiency_pt15to35_points.csv";
const std::string kReferenceRecoPtCsv = kOutDir + "/reference_boxcuts_reco_candidate_tight_fraction_recoPt15to35_points.csv";
const std::string kBaseRecoPtCsv = kOutDir + "/basev3e_target80_reco_candidate_tight_fraction_recoPt15to35_points.csv";
const std::string kReferencePurityCsv = kOutDir + "/reference_boxcuts_raw_sim_purity_recoPt15to35_points.csv";
const std::string kBasePurityCsv = kOutDir + "/basev3e_target80_raw_sim_purity_recoPt15to35_points.csv";

std::vector<std::string> splitCsv(const std::string &line)
{
  std::vector<std::string> fields;
  std::string field;
  bool quoted = false;
  for (char ch : line)
  {
    if (ch == '"')
    {
      quoted = !quoted;
      continue;
    }
    if (ch == ',' && !quoted)
    {
      fields.push_back(field);
      field.clear();
    }
    else
    {
      field.push_back(ch);
    }
  }
  fields.push_back(field);
  return fields;
}

double toDouble(const std::string &s)
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

std::map<std::string, std::vector<Point>> readReference()
{
  std::ifstream in(kReferenceCsv);
  std::map<std::string, std::vector<Point>> out;
  std::string line;
  std::getline(in, line);
  while (std::getline(in, line))
  {
    auto f = splitCsv(line);
    if (f.size() < 9) continue;
    const std::string &cent = f[3];
    const double x = toDouble(f[4]);
    if (x < 15.0 || x > 35.0) continue;
    out[cent].push_back({x, toDouble(f[7]), toDouble(f[8])});
  }
  return out;
}

std::map<std::string, std::vector<Point>> readBase()
{
  std::ifstream in(kBaseCsv);
  std::map<std::string, std::vector<Point>> out;
  std::string line;
  std::getline(in, line);
  while (std::getline(in, line))
  {
    auto f = splitCsv(line);
    if (f.size() < 7) continue;
    out[f[1]].push_back({toDouble(f[4]), toDouble(f[5]), toDouble(f[6])});
  }
  return out;
}

std::map<std::string, std::vector<Point>> readRecoilCsv(const std::string &path)
{
  std::ifstream in(path);
  std::map<std::string, std::vector<Point>> out;
  std::string line;
  std::getline(in, line);
  while (std::getline(in, line))
  {
    auto f = splitCsv(line);
    if (f.size() < 7) continue;
    out[f[1]].push_back({toDouble(f[4]), toDouble(f[5]), toDouble(f[6])});
  }
  return out;
}

TGraphErrors *makeGraph(const std::vector<Point> &pts, int color, int marker, double shift, double markerSize = 2.05)
{
  auto *g = new TGraphErrors(static_cast<int>(pts.size()));
  for (int i = 0; i < static_cast<int>(pts.size()); ++i)
  {
    g->SetPoint(i, pts[i].x + shift, pts[i].y);
    g->SetPointError(i, 0.0, pts[i].ey);
  }
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(marker);
  g->SetMarkerSize(markerSize);
  g->SetLineWidth(2);
  return g;
}

void drawLabelBlock(const std::string &sampleLine, const std::string &definitionLine)
{
  TLatex lat;
  lat.SetNDC();
  lat.SetTextFont(42);
  lat.SetTextSize(0.044);
  lat.DrawLatex(0.070, 0.965, "#it{#bf{sPHENIX}} Internal");
  lat.SetTextSize(0.029);
  lat.DrawLatex(0.070, 0.918, sampleLine.c_str());
  lat.SetTextSize(0.028);
  lat.DrawLatex(0.070, 0.878, definitionLine.c_str());
}

void writeRows(std::ofstream &csv,
               const std::string &source,
               const std::string &label,
               const std::string &cent,
               const std::vector<Point> &pts)
{
  for (const auto &p : pts)
    csv << source << "," << label << "," << cent << "," << p.x << "," << p.y << "," << p.ey << "\n";
}

void drawComparison(const bool includeGlobalBdt)
{
  gStyle->SetOptStat(0);
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  const auto reference = readReference();
  const auto base = readRecoilCsv(kBaseCsv);
  const auto globalBdt = includeGlobalBdt ? readRecoilCsv(kGlobalBdtCsv) : std::map<std::string, std::vector<Point>>{};
  const std::vector<std::pair<std::string, std::string>> cents = {
      {"0_20", "0-20%"},
      {"20_50", "20-50%"},
      {"50_80", "50-80%"},
  };

  TCanvas c(includeGlobalBdt ? "c_basev3e_reco_eff_with_global" : "c_basev3e_reco_eff",
            includeGlobalBdt ? "c_basev3e_reco_eff_with_global" : "c_basev3e_reco_eff",
            3590, 1220);
  c.Divide(3, 1, 0.006, 0.001);

  const std::string csvName = includeGlobalBdt
                                  ? kOutDir + "/basev3e_target80_vs_reference_vs_global32_reco_efficiency_pt15to35_points.csv"
                                  : kOutDir + "/basev3e_target80_vs_reference_boxcuts_reco_efficiency_pt15to35_points.csv";
  std::ofstream csv(csvName);
  csv << "source,label,centrality,pt_mid,value,error\n";

  for (int ic = 0; ic < 3; ++ic)
  {
    c.cd(ic + 1);
    auto *pad = static_cast<TPad *>(gPad);
    pad->SetTopMargin(0.28);
    pad->SetBottomMargin(0.15);
    pad->SetLeftMargin(ic == 0 ? 0.16 : 0.06);
    pad->SetRightMargin(ic == 2 ? 0.040 : 0.025);
    pad->SetGridx(false);
    pad->SetGridy(false);
    pad->SetTicks(1, 1);

    auto *frame = pad->DrawFrame(14.5, 0.0, 35.5, 0.62);
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle(ic == 1 ? "truth photon p_{T} [GeV]" : "");
    frame->GetYaxis()->SetTitle(ic == 0 ? "Reco efficiency" : "");
    frame->GetXaxis()->SetTitleSize(0.055);
    frame->GetYaxis()->SetTitleSize(0.060);
    frame->GetXaxis()->SetLabelSize(0.046);
    frame->GetYaxis()->SetLabelSize(0.046);
    frame->GetYaxis()->SetTitleOffset(1.18);
    frame->GetXaxis()->SetNdivisions(506);
    frame->GetYaxis()->SetNdivisions(506);

    TLatex title;
    title.SetNDC();
    title.SetTextAlign(22);
    title.SetTextFont(62);
    title.SetTextSize(0.064);
    title.DrawLatex(0.52, 0.76, cents[ic].second.c_str());

    const auto &cent = cents[ic].first;
    auto *gRef = makeGraph(reference.at(cent), kBlack, 20, includeGlobalBdt ? -0.16 : -0.10, 2.15);
    auto *gBase = makeGraph(base.at(cent), kAzure + 2, 20, includeGlobalBdt ? 0.00 : 0.10, 2.20);
    TGraphErrors *gGlobal = nullptr;
    if (includeGlobalBdt) gGlobal = makeGraph(globalBdt.at(cent), kOrange + 7, 25, 0.16, 2.15);

    gRef->Draw("PZ SAME");
    gBase->Draw("PZ SAME");
    if (gGlobal) gGlobal->Draw("PZ SAME");

    writeRows(csv, "reference_box_cuts", "Reference box cuts", cent, reference.at(cent));
    writeRows(csv, "basev3e_centrality_target80", "Base v3E + centrality target-80", cent, base.at(cent));
    if (gGlobal) writeRows(csv, "global32_cent3grid_target80", "32-feature BDT, 8 p_{T} x 3 cent target-80", cent, globalBdt.at(cent));

    if (ic == 2)
    {
      auto *leg = new TLegend(includeGlobalBdt ? 0.13 : 0.24,
                              includeGlobalBdt ? 0.16 : 0.18,
                              0.79,
                              includeGlobalBdt ? 0.38 : 0.32);
      leg->SetBorderSize(0);
      leg->SetFillStyle(1001);
      leg->SetFillColorAlpha(kWhite, 0.86);
      leg->SetTextFont(42);
      leg->SetTextSize(includeGlobalBdt ? 0.030 : 0.034);
      leg->AddEntry(gRef, "Reference box cuts", "pe");
      leg->AddEntry(gBase, "Base v3E + centrality", "pe");
      if (gGlobal) leg->AddEntry(gGlobal, "32-feature BDT, 8 p_{T} #times 3 cent", "pe");
      leg->Draw();
    }
  }

  c.cd();
  drawLabelBlock("Embedded Photon12+20 signal, 15 < truth photon p_{T} < 35 GeV",
                 "Reco efficiency = 1 - truth photons missed after reconstruction and ID selection");
  const std::string out = includeGlobalBdt
                              ? kOutDir + "/basev3e_target80_vs_reference_vs_global32_reco_efficiency_pt15to35_1x3.png"
                              : kOutDir + "/basev3e_target80_vs_reference_boxcuts_reco_efficiency_pt15to35_1x3.png";
  c.SaveAs(out.c_str());
  std::cout << out << "\n";
}

void drawRecoPtTightFraction()
{
  gStyle->SetOptStat(0);
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  const auto reference = readRecoilCsv(kReferenceRecoPtCsv);
  const auto base = readRecoilCsv(kBaseRecoPtCsv);
  const std::vector<std::pair<std::string, std::string>> cents = {
      {"0_20", "0-20%"},
      {"20_50", "20-50%"},
      {"50_80", "50-80%"},
  };

  TCanvas c("c_basev3e_reco_pt_tight_fraction", "c_basev3e_reco_pt_tight_fraction", 3590, 1220);
  c.Divide(3, 1, 0.006, 0.001);

  std::ofstream csv(kOutDir + "/basev3e_target80_vs_reference_boxcuts_reco_candidate_tight_fraction_recoPt15to35_points.csv");
  csv << "source,label,centrality,pt_mid,value,error\n";

  for (int ic = 0; ic < 3; ++ic)
  {
    c.cd(ic + 1);
    auto *pad = static_cast<TPad *>(gPad);
    pad->SetTopMargin(0.28);
    pad->SetBottomMargin(0.15);
    pad->SetLeftMargin(ic == 0 ? 0.16 : 0.06);
    pad->SetRightMargin(ic == 2 ? 0.040 : 0.025);
    pad->SetTicks(1, 1);

    auto *frame = pad->DrawFrame(14.5, 0.0, 35.5, 1.04);
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle(ic == 1 ? "reco photon p_{T} [GeV]" : "");
    frame->GetYaxis()->SetTitle(ic == 0 ? "Tight-ID fraction" : "");
    frame->GetXaxis()->SetTitleSize(0.055);
    frame->GetYaxis()->SetTitleSize(0.060);
    frame->GetXaxis()->SetLabelSize(0.046);
    frame->GetYaxis()->SetLabelSize(0.046);
    frame->GetYaxis()->SetTitleOffset(1.18);
    frame->GetXaxis()->SetNdivisions(506);
    frame->GetYaxis()->SetNdivisions(506);

    TLatex title;
    title.SetNDC();
    title.SetTextAlign(22);
    title.SetTextFont(62);
    title.SetTextSize(0.064);
    title.DrawLatex(0.52, 0.76, cents[ic].second.c_str());

    const auto &cent = cents[ic].first;
    auto *gRef = makeGraph(reference.at(cent), kBlack, 20, -0.10, 2.15);
    auto *gBase = makeGraph(base.at(cent), kAzure + 2, 20, 0.10, 2.20);
    gRef->Draw("PZ SAME");
    gBase->Draw("PZ SAME");

    writeRows(csv, "reference_box_cuts", "Reference box cuts", cent, reference.at(cent));
    writeRows(csv, "basev3e_centrality_target80", "Base v3E + centrality target-80", cent, base.at(cent));

    if (ic == 2)
    {
      auto *leg = new TLegend(0.24, 0.18, 0.79, 0.32);
      leg->SetBorderSize(0);
      leg->SetFillStyle(1001);
      leg->SetFillColorAlpha(kWhite, 0.86);
      leg->SetTextFont(42);
      leg->SetTextSize(0.034);
      leg->AddEntry(gRef, "Reference box cuts", "pe");
      leg->AddEntry(gBase, "Base v3E + centrality", "pe");
      leg->Draw();
    }
  }

  c.cd();
  drawLabelBlock("Embedded Photon12+20 signal, 15 < reco photon p_{T} < 35 GeV",
                 "Tight-ID fraction among truth-matched reconstructed photon candidates");
  const std::string out = kOutDir + "/basev3e_target80_vs_reference_boxcuts_reco_candidate_tight_fraction_recoPt15to35_1x3.png";
  c.SaveAs(out.c_str());
  std::cout << out << "\n";
}

void drawRawSimPurity()
{
  gStyle->SetOptStat(0);
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  const auto reference = readRecoilCsv(kReferencePurityCsv);
  const auto base = readRecoilCsv(kBasePurityCsv);
  const std::vector<std::pair<std::string, std::string>> cents = {
      {"0_20", "0-20%"},
      {"20_50", "20-50%"},
      {"50_80", "50-80%"},
  };

  TCanvas c("c_basev3e_raw_sim_purity", "c_basev3e_raw_sim_purity", 3590, 1220);
  c.Divide(3, 1, 0.006, 0.001);

  std::ofstream csv(kOutDir + "/basev3e_target80_vs_reference_boxcuts_raw_sim_purity_recoPt15to35_points.csv");
  csv << "source,label,centrality,pt_mid,value,error\n";

  for (int ic = 0; ic < 3; ++ic)
  {
    c.cd(ic + 1);
    auto *pad = static_cast<TPad *>(gPad);
    pad->SetTopMargin(0.28);
    pad->SetBottomMargin(0.15);
    pad->SetLeftMargin(ic == 0 ? 0.16 : 0.06);
    pad->SetRightMargin(ic == 2 ? 0.040 : 0.025);
    pad->SetTicks(1, 1);

    auto *frame = pad->DrawFrame(14.5, 0.78, 35.5, 0.96);
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle(ic == 1 ? "reco photon p_{T} [GeV]" : "");
    frame->GetYaxis()->SetTitle(ic == 0 ? "Raw SIM purity" : "");
    frame->GetXaxis()->SetTitleSize(0.055);
    frame->GetYaxis()->SetTitleSize(0.060);
    frame->GetXaxis()->SetLabelSize(0.046);
    frame->GetYaxis()->SetLabelSize(0.046);
    frame->GetYaxis()->SetTitleOffset(1.18);
    frame->GetXaxis()->SetNdivisions(506);
    frame->GetYaxis()->SetNdivisions(506);

    TLatex title;
    title.SetNDC();
    title.SetTextAlign(22);
    title.SetTextFont(62);
    title.SetTextSize(0.064);
    title.DrawLatex(0.52, 0.76, cents[ic].second.c_str());

    const auto &cent = cents[ic].first;
    auto *gRef = makeGraph(reference.at(cent), kBlack, 20, -0.10, 2.15);
    auto *gBase = makeGraph(base.at(cent), kAzure + 2, 20, 0.10, 2.20);
    gRef->Draw("PZ SAME");
    gBase->Draw("PZ SAME");

    writeRows(csv, "reference_box_cuts", "Reference box cuts", cent, reference.at(cent));
    writeRows(csv, "basev3e_centrality_target80", "Base v3E + centrality target-80", cent, base.at(cent));

    if (ic == 2)
    {
      auto *leg = new TLegend(0.21, 0.18, 0.79, 0.32);
      leg->SetBorderSize(0);
      leg->SetFillStyle(1001);
      leg->SetFillColorAlpha(kWhite, 0.86);
      leg->SetTextFont(42);
      leg->SetTextSize(0.034);
      leg->AddEntry(gRef, "Reference box cuts", "pe");
      leg->AddEntry(gBase, "Base v3E + centrality", "pe");
      leg->Draw();
    }
  }

  c.cd();
  drawLabelBlock("Embedded Photon12+20 signal and inclusive-jet background, 15 < reco photon p_{T} < 35 GeV",
                 "Raw SIM purity = S_{tight}^{truth-matched} / (S_{tight}^{truth-matched} + B_{tight}^{inclusive})");
  const std::string out = kOutDir + "/basev3e_target80_vs_reference_boxcuts_raw_sim_purity_recoPt15to35_1x3.png";
  c.SaveAs(out.c_str());
  std::cout << out << "\n";
}
}

void PlotBaseV3ERecoEfficiencyOverlay()
{
  drawComparison(false);
  drawComparison(true);
  drawRecoPtTightFraction();
  drawRawSimPurity();
}
