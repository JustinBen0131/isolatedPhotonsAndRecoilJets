#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TH1.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace
{
// IMPORTANT:
//   After RecoilJets_AuAu stitching, PhotonJet12 must use the EXCLUSIVE
//   generator cross section for 12 <= pT_filter^gamma < 20 GeV.
//   Do not use the old inclusive PhotonJet12 value 2719.20189 pb here.
constexpr double kSigmaEmbeddedPhoton12To20_pb = 2661.18552;
constexpr double kSigmaEmbeddedPhoton20Plus_pb = 133.317866;

const std::string kConfigTag =
    "jetMinPt5_7pi_8_vz60_isoR40_fixedIso2GeV_baseVariant_"
    "preselectionReference_tightReference_nonTightReference";

const std::string kFile12 =
    "InputFiles/simEmbedded/RecoilJets_embeddedPhoton12_ALL_" + kConfigTag + ".root";
const std::string kFile20 =
    "InputFiles/simEmbedded/RecoilJets_embeddedPhoton20_ALL_" + kConfigTag + ".root";

const std::string kOutDir =
    "dataOutput/combinedSimOnlyEMBEDDED/" + kConfigTag + "/photonJet12and20merged_SIM";

double ReadEventCount(TFile* f)
{
  if (!f) return 0.0;
  TDirectory* d = f->GetDirectory("SIM");
  if (!d) return 0.0;

  TH1* cnt = dynamic_cast<TH1*>(d->Get("cnt_SIM"));
  if (!cnt) return 0.0;

  const double bin1 = cnt->GetBinContent(1);
  if (bin1 > 0.0) return bin1;

  const double integral = cnt->Integral(0, cnt->GetNbinsX() + 1);
  if (integral > 0.0) return integral;

  return cnt->GetEntries();
}

TH1* CloneFromDir(TDirectory* d, const std::string& name, const std::string& cloneName)
{
  if (!d) return nullptr;
  TH1* h = dynamic_cast<TH1*>(d->Get(name.c_str()));
  if (!h) return nullptr;

  TH1* c = dynamic_cast<TH1*>(h->Clone(cloneName.c_str()));
  if (!c) return nullptr;
  c->SetDirectory(nullptr);
  c->Sumw2();
  return c;
}

std::unique_ptr<TH1> AddABCDSet(TDirectory* d,
                                const std::vector<std::string>& suffixes,
                                const std::string& cloneName)
{
  const std::vector<std::string> regions = {"A", "B", "C", "D"};
  std::unique_ptr<TH1> sum;

  for (const std::string& suffix : suffixes)
  {
    for (const std::string& region : regions)
    {
      const std::string name = "h_pTgamma_ABCD_" + region + suffix;
      std::unique_ptr<TH1> h(CloneFromDir(d, name, cloneName + "_" + region + suffix));
      if (!h) continue;

      if (!sum)
      {
        sum.reset(dynamic_cast<TH1*>(h->Clone(cloneName.c_str())));
        if (!sum) return nullptr;
        sum->SetDirectory(nullptr);
        sum->Reset("ICES");
        sum->Sumw2();
      }
      sum->Add(h.get());
    }
  }

  return sum;
}

std::unique_ptr<TH1> BuildABCDPtSpectrum(TFile* f, const std::string& cloneName)
{
  if (!f) return nullptr;
  TDirectory* d = f->GetDirectory("SIM");
  if (!d) return nullptr;

  // Prefer the inclusive ABCD spectra if available. If only centrality slices
  // exist, fall back to the three standard AuAu centrality groups.
  std::unique_ptr<TH1> inclusive = AddABCDSet(d, {""}, cloneName);
  if (inclusive && inclusive->Integral(0, inclusive->GetNbinsX() + 1) > 0.0)
  {
    return inclusive;
  }

  return AddABCDSet(d, {"_cent_0_20", "_cent_20_50", "_cent_50_80"}, cloneName);
}

void StyleHist(TH1* h, int color, int marker, const char* title)
{
  if (!h) return;
  h->SetTitle(title);
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  h->SetMarkerSize(0.9);
  h->SetLineWidth(2);
}

std::string Sci(double v, int p = 3)
{
  std::ostringstream os;
  os << std::scientific << std::setprecision(p) << v;
  return os.str();
}
}

void MakeEmbeddedPhotonXsecNormQA()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0.5);

  std::unique_ptr<TFile> f12(TFile::Open(kFile12.c_str(), "READ"));
  std::unique_ptr<TFile> f20(TFile::Open(kFile20.c_str(), "READ"));
  if (!f12 || f12->IsZombie())
  {
    std::cerr << "[ERROR] Cannot open " << kFile12 << std::endl;
    return;
  }
  if (!f20 || f20->IsZombie())
  {
    std::cerr << "[ERROR] Cannot open " << kFile20 << std::endl;
    return;
  }

  const double n12 = ReadEventCount(f12.get());
  const double n20 = ReadEventCount(f20.get());
  if (n12 <= 0.0 || n20 <= 0.0)
  {
    std::cerr << "[ERROR] Bad event counts: N12=" << n12 << " N20=" << n20 << std::endl;
    return;
  }

  std::unique_ptr<TH1> h12 = BuildABCDPtSpectrum(f12.get(), "h_embeddedPhoton12_weightedABCDPt");
  std::unique_ptr<TH1> h20 = BuildABCDPtSpectrum(f20.get(), "h_embeddedPhoton20_weightedABCDPt");
  if (!h12 || !h20)
  {
    std::cerr << "[ERROR] Could not build ABCD photon pT spectra from inputs." << std::endl;
    return;
  }

    const double raw12 = h12->Integral(0, h12->GetNbinsX() + 1);
    const double raw20 = h20->Integral(0, h20->GetNbinsX() + 1);

    if (kSigmaEmbeddedPhoton12To20_pb <= 0.0)
    {
      std::cerr << "[ERROR] kSigmaEmbeddedPhoton12To20_pb is not set." << std::endl;
      std::cerr << "        Rerun the Pythia xsec estimate for the exclusive window" << std::endl;
      std::cerr << "        12 <= pT_filter^gamma < 20 GeV and put that value here." << std::endl;
      std::cerr << "        Do NOT use the inclusive PhotonJet12 value 2719.20189 pb." << std::endl;
      return;
    }

    const double w12 = kSigmaEmbeddedPhoton12To20_pb / n12;
    const double w20 = kSigmaEmbeddedPhoton20Plus_pb / n20;

    h12->Scale(w12);
  h20->Scale(w20);

  std::unique_ptr<TH1> hSum(dynamic_cast<TH1*>(h12->Clone("h_embeddedPhoton12plus20_weightedABCDPt")));
  hSum->SetDirectory(nullptr);
  hSum->Add(h20.get());

  std::unique_ptr<TH1> frac12(dynamic_cast<TH1*>(h12->Clone("h_fractionPhoton12")));
  std::unique_ptr<TH1> frac20(dynamic_cast<TH1*>(h20->Clone("h_fractionPhoton20")));
  frac12->SetDirectory(nullptr);
  frac20->SetDirectory(nullptr);
  frac12->Divide(hSum.get());
  frac20->Divide(hSum.get());

  StyleHist(h12.get(), kBlue + 1, 20, "");
  StyleHist(h20.get(), kRed + 1, 21, "");
  StyleHist(hSum.get(), kBlack, 24, "");
  StyleHist(frac12.get(), kBlue + 1, 20, "");
  StyleHist(frac20.get(), kRed + 1, 21, "");

  hSum->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
  hSum->GetYaxis()->SetTitle("#sigma_{eff}/N scaled entries [pb / bin]");
  hSum->GetYaxis()->SetTitleOffset(1.12);

  frac12->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
  frac12->GetYaxis()->SetTitle("sample fraction");
  frac12->GetYaxis()->SetRangeUser(0.0, 1.08);
  frac12->GetYaxis()->SetNdivisions(505);
  frac12->GetYaxis()->SetTitleSize(0.085);
  frac12->GetYaxis()->SetLabelSize(0.075);
  frac12->GetYaxis()->SetTitleOffset(0.50);
  frac12->GetXaxis()->SetTitleSize(0.090);
  frac12->GetXaxis()->SetLabelSize(0.080);

    gSystem->mkdir(kOutDir.c_str(), true);
    const std::string outPng = kOutDir + "/embeddedPhoton12to20_plusPhoton20_xsecNormalizationQA.png";

    TCanvas c("c_embeddedPhotonXsecNormQA", "embedded photon stitched xsec normalization QA", 1050, 850);
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
  const double ymax = std::max({hSum->GetMaximum(), h12->GetMaximum(), h20->GetMaximum()});
  hSum->SetMinimum(std::max(1.0e-8, ymax * 2.0e-5));
  hSum->SetMaximum(std::max(1.0e-6, ymax * 18.0));
  hSum->Draw("E1");
  h12->Draw("E1 SAME");
  h20->Draw("E1 SAME");

    TLegend leg(0.55, 0.63, 0.91, 0.87);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.040);
    leg.AddEntry(hSum.get(), "stitched PhotonJet12 + PhotonJet20", "lep");
    leg.AddEntry(h12.get(), ("PhotonJet12, 12 #leq p_{T,filter}^{#gamma} < 20, #sigma_{eff}=" + Sci(kSigmaEmbeddedPhoton12To20_pb, 3) + " pb").c_str(), "lep");
    leg.AddEntry(h20.get(), ("PhotonJet20, p_{T,filter}^{#gamma} #geq 20, #sigma_{eff}=" + Sci(kSigmaEmbeddedPhoton20Plus_pb, 3) + " pb").c_str(), "lep");
    leg.Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.040);
    lat.DrawLatex(0.15, 0.86, "Embedded photon stitched normalization QA");
    lat.SetTextSize(0.032);
    lat.DrawLatex(0.15, 0.80, "Reference ID; stitched summed ABCD photon spectra");
    lat.DrawLatex(0.15, 0.75, ("w_{12#rightarrow20}=" + Sci(w12) + " pb/event, w_{20+}=" + Sci(w20) + " pb/event").c_str());
    lat.DrawLatex(0.15, 0.70, ("stitched raw ABCD entries: 12#rightarrow20=" + Sci(raw12, 2) + ", 20+=" + Sci(raw20, 2)).c_str());

  bot.cd();
  frac12->Draw("E1");
  frac20->Draw("E1 SAME");
  TLine line(frac12->GetXaxis()->GetXmin(), 0.5, frac12->GetXaxis()->GetXmax(), 0.5);
  line.SetLineColor(kGray + 1);
  line.SetLineStyle(2);
  line.Draw("SAME");

  TLegend leg2(0.58, 0.70, 0.91, 0.93);
  leg2.SetBorderSize(0);
  leg2.SetFillStyle(0);
  leg2.SetTextSize(0.075);
  leg2.AddEntry(frac12.get(), "PhotonJet12 / sum", "lep");
  leg2.AddEntry(frac20.get(), "PhotonJet20 / sum", "lep");
  leg2.Draw();

  c.SaveAs(outPng.c_str());

  std::cout << "[DONE] Wrote " << outPng << std::endl;
  std::cout << "[INFO] N12=" << n12 << " N20=" << n20
            << " w12=" << w12 << " pb/event"
            << " w20=" << w20 << " pb/event" << std::endl;
}
