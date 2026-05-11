#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TNamed.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <regex>
#include <sstream>
#include <string>

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

bool ExtractJet20RefScalePbPerEvent(const std::string& mergeInfo, double& scale)
{
  scale = 0.0;
  const std::regex jet20Pattern(R"(\[embeddedJet20 Nraw=([0-9.eE+\-]+) sigma_pb=([0-9.eE+\-]+) w=([0-9.eE+\-]+)\])");
  std::smatch match;
  if (!std::regex_search(mergeInfo, match, jet20Pattern) || match.size() < 4) return false;

  const double nRaw = std::stod(match[1].str());
  const double sigmaPb = std::stod(match[2].str());
  if (!(nRaw > 0.0) || !(sigmaPb > 0.0)) return false;

  scale = sigmaPb / nRaw;
  return std::isfinite(scale) && scale > 0.0;
}

std::unique_ptr<TH1> RatioToSmooth(const TH1* h, const TH1* smooth)
{
  if (!h || !smooth) return nullptr;
  std::unique_ptr<TH1> ratio(dynamic_cast<TH1*>(h->Clone("h_embeddedInclusiveJetFinalStitch_ratio")));
  if (!ratio) return nullptr;
  ratio->SetDirectory(nullptr);
  ratio->Reset("ICES");
  if (ratio->GetSumw2N() == 0) ratio->Sumw2();

  for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
  {
    const double den = smooth->GetBinContent(ib);
    const double num = h->GetBinContent(ib);
    if (den <= 0.0 || num <= 0.0) continue;
    ratio->SetBinContent(ib, num / den);
    ratio->SetBinError(ib, h->GetBinError(ib) / den);
  }
  return ratio;
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
  double refScalePbPerEvent = 1.0;
  const bool haveAbsScale = ExtractJet20RefScalePbPerEvent(mergeInfo, refScalePbPerEvent);
  if (!haveAbsScale)
  {
    std::cerr << "[WARN] Could not parse embeddedJet20 scale from MERGE_INFO; plotting relative merged entries." << std::endl;
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

  std::unique_ptr<TH1> smooth(dynamic_cast<TH1*>(h->Clone("h_embeddedInclusiveJetFinalStitch_smooth")));
  smooth->SetDirectory(nullptr);
  smooth->Smooth(2);
  smooth->SetLineColor(kGray + 2);
  smooth->SetLineStyle(2);
  smooth->SetLineWidth(2);
  smooth->SetMarkerSize(0);
  smooth->SetStats(false);

  std::unique_ptr<TH1> ratio = RatioToSmooth(h.get(), smooth.get());
  ratio->SetStats(false);
  ratio->SetLineColor(kBlack);
  ratio->SetMarkerColor(kBlack);
  ratio->SetMarkerStyle(20);
  ratio->SetMarkerSize(0.85);
  ratio->SetLineWidth(2);

  const double xMin = 10.0;
  const double xMax = 45.0;
  h->GetXaxis()->SetRangeUser(xMin, xMax);
  h->GetYaxis()->SetTitle(haveAbsScale ? "#sigma_{eff}/N scaled entries [pb / bin]"
                                       : "weighted merged entries / bin");
  h->GetYaxis()->SetTitleOffset(1.12);
  h->SetMinimum(std::max(1.0e-8, h->GetMaximum() * 1.0e-5));
  h->SetMaximum(std::max(1.0e-6, h->GetMaximum() * 70.0));

  ratio->GetXaxis()->SetRangeUser(xMin, xMax);
  ratio->GetXaxis()->SetTitle("max p_{T}^{jet,truth} [GeV]");
  ratio->GetYaxis()->SetTitle("sum / smooth");
  ratio->GetYaxis()->SetRangeUser(0.85, 1.15);
  ratio->GetYaxis()->SetNdivisions(505);
  ratio->GetYaxis()->SetTitleSize(0.085);
  ratio->GetYaxis()->SetLabelSize(0.075);
  ratio->GetYaxis()->SetTitleOffset(0.50);
  ratio->GetXaxis()->SetTitleSize(0.090);
  ratio->GetXaxis()->SetLabelSize(0.080);

  TCanvas c("c_embeddedInclusiveJet_finalMergedStitchSpectrum",
            "embedded inclusive jet final merged stitch spectrum", 1050, 850);
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
  smooth->Draw("HIST SAME");

  TLegend leg(0.15, 0.17, 0.52, 0.32);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.038);
  leg.AddEntry(h.get(), "Jet12+20 stitched", "lep");
  leg.AddEntry(smooth.get(), "Smoothed reference", "l");
  leg.Draw();

  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(0.040);
  lat.DrawLatex(0.15, 0.86, "Embedded inclusive-jet generator stitching spectrum");
  lat.SetTextSize(0.031);
  lat.DrawLatex(0.15, 0.80, "Final pulled RecoilJets_embeddedJet12plus20_MERGED.root");
  lat.DrawLatex(0.15, 0.75, "Jet12: 12 #leq max p_{T}^{jet,truth} < 20 GeV; Jet20: max p_{T}^{jet,truth} #geq 20 GeV");
  if (haveAbsScale)
  {
    lat.DrawLatex(0.15, 0.70, ("Final stitch scale = " + Sci(refScalePbPerEvent) + " pb per Jet20-weighted event").c_str());
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
  ratio->Draw("E1");
  TLine one(xMin, 1.0, xMax, 1.0);
  one.SetLineColor(kGray + 1);
  one.SetLineStyle(2);
  one.Draw("SAME");

  const std::string outPng = ParentDir(inPath) + "/embeddedInclusiveJet_finalMergedStitchedTruthJetPtSpectrum.png";
  c.SaveAs(outPng.c_str());

  std::cout << "[DONE] Wrote " << outPng << std::endl;
  std::cout << "[INFO] input=" << inPath << std::endl;
  std::cout << "[INFO] mergeInfo=" << mergeInfo << std::endl;
  std::cout << "[INFO] refScalePbPerEvent=" << refScalePbPerEvent << std::endl;
}
