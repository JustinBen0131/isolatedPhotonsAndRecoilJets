#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"

namespace
{
struct Point
{
  double x = 0.0;
  double ex = 0.0;
  double y = 0.0;
  double ey = 0.0;
};

double Count(TDirectory* dir, const std::string& name)
{
  if (!dir) return 0.0;
  TH1* h = dynamic_cast<TH1*>(dir->Get(name.c_str()));
  return h ? h->GetBinContent(1) : 0.0;
}

double RawPurity(double A, double B, double C, double D)
{
  if (A <= 0.0 || D <= 0.0) return 0.0;
  double Asig = A - B * (C / D);
  if (Asig < 0.0) Asig = 0.0;
  return Asig / A;
}

double RawPurityError(double A, double B, double C, double D)
{
  if (A <= 0.0 || D <= 0.0) return 0.0;

  const double dPdA =  (B * C) / (A * A * D);
  const double dPdB = -(C) / (A * D);
  const double dPdC = -(B) / (A * D);
  const double dPdD =  (B * C) / (A * D * D);

  double var = 0.0;
  if (A > 0.0) var += dPdA * dPdA * A;
  if (B > 0.0) var += dPdB * dPdB * B;
  if (C > 0.0) var += dPdC * dPdC * C;
  if (D > 0.0) var += dPdD * dPdD * D;
  return (var > 0.0) ? std::sqrt(var) : 0.0;
}

std::vector<Point> BuildPurity(TDirectory* dir, const std::string& cent)
{
  const std::vector<int> edges = {10, 12, 14, 16, 18, 20, 22, 24, 26, 35};
  std::vector<Point> pts;
  pts.reserve(edges.size() - 1);

  for (std::size_t i = 0; i + 1 < edges.size(); ++i)
  {
    const int lo = edges[i];
    const int hi = edges[i + 1];
    const std::string suffix = TString::Format("_pT_%d_%d_cent_%s", lo, hi, cent.c_str()).Data();

    const double A = Count(dir, "h_isIsolated_isTight" + suffix);
    const double B = Count(dir, "h_notIsolated_isTight" + suffix);
    const double C = Count(dir, "h_isIsolated_notTight" + suffix);
    const double D = Count(dir, "h_notIsolated_notTight" + suffix);

    if (A <= 0.0 && B <= 0.0 && C <= 0.0 && D <= 0.0) continue;

    Point p;
    p.x = 0.5 * (lo + hi);
    p.ex = 0.5 * (hi - lo);
    p.y = RawPurity(A, B, C, D);
    p.ey = RawPurityError(A, B, C, D);
    pts.push_back(p);

    std::cout << cent << " " << lo << "-" << hi
              << " A=" << A << " B=" << B << " C=" << C << " D=" << D
              << " P=" << p.y << " +/- " << p.ey << "\n";
  }
  return pts;
}

TGraphErrors* MakeGraph(const std::vector<Point>& pts, int color, int marker, const char* name)
{
  auto* g = new TGraphErrors(static_cast<int>(pts.size()));
  g->SetName(name);
  for (int i = 0; i < static_cast<int>(pts.size()); ++i)
    g->SetPoint(i, pts[i].x, pts[i].y), g->SetPointError(i, pts[i].ex, pts[i].ey);
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(marker);
  g->SetMarkerSize(1.25);
  g->SetLineWidth(2);
  return g;
}

void DrawOne(TDirectory* dBox, TDirectory* dNPB, const std::string& cent, const std::string& outDir)
{
  std::cout << "\n[Box-cuts]\n";
  const auto box = BuildPurity(dBox, cent);
  std::cout << "\n[NPB presel]\n";
  const auto npb = BuildPurity(dNPB, cent);

  auto* gBox = MakeGraph(box, kBlack, 20, ("g_box_" + cent).c_str());
  auto* gNPB = MakeGraph(npb, kBlue + 1, 20, ("g_npb_" + cent).c_str());

  TCanvas c(("c_" + cent).c_str(), ("c_" + cent).c_str(), 1100, 820);
  c.SetFillColor(kWhite);
  c.SetLeftMargin(0.13);
  c.SetRightMargin(0.04);
  c.SetBottomMargin(0.13);
  c.SetTopMargin(0.10);
  c.SetTicks(1, 1);

  TH1F frame(("frame_" + cent).c_str(), "", 100, 10.0, 35.0);
  frame.SetDirectory(nullptr);
  frame.SetStats(false);
  frame.SetMinimum(0.0);
  frame.SetMaximum(1.05);
  frame.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
  frame.GetYaxis()->SetTitle("Purity (raw ABCD)");
  frame.GetXaxis()->SetTitleSize(0.050);
  frame.GetYaxis()->SetTitleSize(0.052);
  frame.GetXaxis()->SetLabelSize(0.043);
  frame.GetYaxis()->SetLabelSize(0.043);
  frame.GetXaxis()->SetTitleOffset(1.05);
  frame.GetYaxis()->SetTitleOffset(1.12);
  frame.Draw();

  gBox->Draw("P SAME");
  gNPB->Draw("P SAME");

  TLegend leg(0.18, 0.20, 0.46, 0.35);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.042);
  leg.AddEntry(gBox, "box-cuts", "pe");
  leg.AddEntry(gNPB, "NPB presel", "pe");
  leg.Draw();

  TLatex t;
  t.SetNDC(true);
  t.SetTextFont(42);
  t.SetTextColor(kBlack);
  t.SetTextAlign(13);
  t.SetTextSize(0.044);
  std::string centLabel = cent;
  std::replace(centLabel.begin(), centLabel.end(), '_', '-');
  t.DrawLatex(0.16, 0.885, TString::Format("Au+Au %s%% centrality", centLabel.c_str()).Data());
  t.SetTextSize(0.033);
  t.DrawLatex(0.16, 0.835, "Photon 10 GeV + MBD N&S #geq 2");
  t.DrawLatex(0.16, 0.790, "|v_{z}| < 150 cm, sliding isolation, R_{cone} < 0.4");

  t.SetTextAlign(31);
  t.SetTextSize(0.036);
  t.SetTextFont(62);
  const double sphY = (cent == "50_80") ? 0.350 : 0.725;
  t.DrawLatex(0.93, sphY, "sPHENIX");
  t.SetTextFont(52);
  t.DrawLatex(0.93, sphY - 0.040, "Internal");

  const std::string out = outDir + "/purity_boxCuts_vs_NPBPresel_cent_" + cent + ".png";
  c.SaveAs(out.c_str());

  delete gBox;
  delete gNPB;
}
}

void CompareAuAuBoxToNPB()
{
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(4);
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");

  const std::string boxFile =
      "InputFiles/auau25/RecoilJets_auau_ALL_jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionReference_tightReference_nonTightReference.root";
  const std::string npbFile =
      "InputFiles/auau25/RecoilJets_auau_ALL_jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionVariantA_tightReference_nonTightReference.root";
  const std::string trig = "photon_10_plus_MBD_NS_geq_2_vtx_lt_150";
  const std::string outDir = "dataOutput/auau/compareBoxToNPB";

  gSystem->mkdir(outDir.c_str(), true);

  TFile fBox(boxFile.c_str(), "READ");
  TFile fNPB(npbFile.c_str(), "READ");
  if (fBox.IsZombie() || fNPB.IsZombie())
  {
    std::cerr << "Failed to open input files\n";
    return;
  }

  TDirectory* dBox = fBox.GetDirectory(trig.c_str());
  TDirectory* dNPB = fNPB.GetDirectory(trig.c_str());
  if (!dBox || !dNPB)
  {
    std::cerr << "Missing trigger directory: " << trig << "\n";
    return;
  }

  for (const std::string& cent : {"0_20", "20_50", "50_80"})
    DrawOne(dBox, dNPB, cent, outDir);
}
