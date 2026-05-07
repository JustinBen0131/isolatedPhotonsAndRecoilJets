#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TSystem.h>
#include <TStyle.h>

namespace
{
std::string envOrDefault(const char* name, const std::string& fallback)
{
  if (const char* value = std::getenv(name))
  {
    if (*value) return value;
  }
  return fallback;
}

std::string joinPath(const std::string& a, const std::string& b)
{
  if (a.empty()) return b;
  if (a.back() == '/') return a + b;
  return a + "/" + b;
}

TH1* loadHist(TFile& f, const std::string& key)
{
  const std::string name = "h_maxEnergyClus_NewTriggerFilling_perRunCorrected_" + key;
  TDirectory* dir = dynamic_cast<TDirectory*>(f.Get(key.c_str()));
  if (!dir) return nullptr;
  TH1* src = dynamic_cast<TH1*>(dir->Get(name.c_str()));
  if (!src) return nullptr;
  TH1* h = dynamic_cast<TH1*>(src->Clone((name + "_plotClone").c_str()));
  if (h) h->SetDirectory(nullptr);
  return h;
}

TH1* normalizedClone(TH1* src, const char* name)
{
  if (!src) return nullptr;
  TH1* h = dynamic_cast<TH1*>(src->Clone(name));
  if (!h) return nullptr;
  h->SetDirectory(nullptr);
  const double integral = h->Integral(0, h->GetNbinsX() + 1);
  if (integral > 0.0) h->Scale(1.0 / integral);
  return h;
}

void styleHist(TH1* h, int color, int markerStyle)
{
  if (!h) return;
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetLineWidth(3);
  h->SetMarkerStyle(markerStyle);
  h->SetMarkerSize(1.0);
  h->SetStats(false);
}

double smallestPositiveBinContent(TH1* h)
{
  if (!h) return 0.0;
  double best = 0.0;
  for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
  {
    const double y = h->GetBinContent(ib);
    if (!std::isfinite(y) || y <= 0.0) continue;
    if (best <= 0.0 || y < best) best = y;
  }
  return best;
}

double findCrossingX(TH1* h, double targetY, double xLo, double xHi)
{
  if (!h || !std::isfinite(targetY)) return -1.0;

  for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
  {
    const double x = h->GetXaxis()->GetBinCenter(ib);
    if (x < xLo || x > xHi) continue;

    const double y = h->GetBinContent(ib);
    if (!std::isfinite(y)) continue;
    if (y >= targetY) return x;
  }
  return -1.0;
}

double tailMean(TH1* h, double xTailLo)
{
  if (!h) return 0.90;

  const int bLo = h->GetXaxis()->FindBin(xTailLo + 1e-6);
  const int bHi = h->GetXaxis()->FindBin(20.0 - 1e-6);

  double sumW = 0.0;
  double sumWY = 0.0;
  for (int ib = bLo; ib <= bHi; ++ib)
  {
    const double y = h->GetBinContent(ib);
    const double ey = h->GetBinError(ib);
    if (!std::isfinite(y) || y <= 0.0) continue;

    const double w = (std::isfinite(ey) && ey > 0.0) ? 1.0 / (ey * ey) : 1.0;
    sumW += w;
    sumWY += w * y;
  }

  if (sumW <= 0.0) return 0.90;
  return std::max(0.80, std::min(1.0, sumWY / sumW));
}

double fitPlateauConst(TH1* h, double xTailLo)
{
  if (!h) return 0.90;
  xTailLo = std::max(9.0, std::min(19.0, xTailLo));

  const int bLo = h->GetXaxis()->FindBin(xTailLo + 1e-6);
  const int bHi = h->GetXaxis()->FindBin(20.0 - 1e-6);
  if (bHi - bLo + 1 < 3) return tailMean(h, xTailLo);

  static int counter = 0;
  TF1 fPlateau(Form("f_scaledTrig_plateau_tmp_%d", counter++), "pol0", xTailLo, 20.0);
  TFitResultPtr fitRes = h->Fit(&fPlateau, "RQ0NS", "", xTailLo, 20.0);

  double plateau = tailMean(h, xTailLo);
  if (static_cast<int>(fitRes) == 0)
  {
    const double p0 = fPlateau.GetParameter(0);
    if (std::isfinite(p0) && p0 > 0.0) plateau = p0;
  }

  return std::max(0.80, std::min(1.0, plateau));
}

double rangeMean(TH1* h, double xLo, double xHi, double fallback)
{
  if (!h) return fallback;
  const int bLo = h->GetXaxis()->FindBin(xLo + 1e-6);
  const int bHi = h->GetXaxis()->FindBin(xHi - 1e-6);
  double sumWY = 0.0;
  double sumW = 0.0;
  for (int ib = bLo; ib <= bHi; ++ib)
  {
    const double y = h->GetBinContent(ib);
    if (!std::isfinite(y)) continue;
    double err = h->GetBinError(ib);
    if (!std::isfinite(err) || err <= 0.0) err = 0.03;
    const double w = 1.0 / (err * err);
    sumWY += w * y;
    sumW += w;
  }
  if (sumW <= 0.0) return fallback;
  return sumWY / sumW;
}

TF1* buildTurnOnFit(TH1* hRatio, int thresholdGeV, int color,
                    double plateauVal, double floorVal)
{
  const double defaultX0 = (thresholdGeV == 10) ? 5.2 : 6.2;
  const double defaultSlope = (thresholdGeV == 10) ? 1.1 : 1.0;

  const double span = std::max(0.1, plateauVal - floorVal);
  const double x20 = findCrossingX(hRatio, floorVal + 0.20 * span, 1.0, 15.0);
  const double x50 = findCrossingX(hRatio, floorVal + 0.50 * span, 1.0, 15.0);
  const double x80 = findCrossingX(hRatio, floorVal + 0.80 * span, 1.0, 15.0);

  double slopeGuess = defaultSlope;
  if (std::isfinite(x20) && std::isfinite(x80) && x80 > x20)
  {
    slopeGuess = 2.7726 / (x80 - x20);
  }
  slopeGuess = std::max(0.30, std::min(3.0, slopeGuess));

  double x0Guess = defaultX0;
  if (std::isfinite(x50) && x50 > 0.0)
  {
    x0Guess = x50;
  }

  TF1* f = new TF1(Form("f_scaledTrig_turnon_%d", thresholdGeV),
                   "[3]+([0]-[3])/pow(1+exp(-[1]*(x-[2])),[4])",
                   0.0, 20.0);
  if (!f) return nullptr;

  f->SetParNames("Plateau", "Slope", "X0", "Floor", "Shape");
  f->SetParameter(0, plateauVal);
  f->SetParameter(1, slopeGuess);
  f->SetParameter(2, x0Guess);
  f->SetParameter(3, floorVal);
  f->SetParameter(4, 1.0);
  f->FixParameter(0, plateauVal);
  f->SetParLimits(1, 0.20, 3.5);
  f->SetParLimits(2, std::max(1.0, x0Guess - 3.0), std::min(14.0, x0Guess + 3.0));
  f->SetParLimits(3, 0.0, std::min(0.08, floorVal + 0.04));
  f->SetParLimits(4, 0.35, 3.5);
  f->SetNpx(500);
  f->SetLineColor(color);
  f->SetLineStyle(2);
  f->SetLineWidth(3);
  return f;
}

double fitTurnOnAndGetX90(TH1* hRatio, int thresholdGeV, int color,
                          double& plateauVal, std::unique_ptr<TF1>& fitOut)
{
  fitOut.reset();
  if (!hRatio) return -1.0;

  plateauVal = fitPlateauConst(hRatio, 15.0);
  const double floorVal = std::max(0.0, std::min(0.06, rangeMean(hRatio, 1.0, 2.5, 0.01)));
  std::unique_ptr<TF1> f(buildTurnOnFit(hRatio, thresholdGeV, color, plateauVal, floorVal));
  if (!f) return findCrossingX(hRatio, 0.90, 2.0, 20.0);

  const double fitLo = std::max(1.5, findCrossingX(hRatio, floorVal + 0.05 * (plateauVal - floorVal), 1.0, 8.0) - 1.0);
  const double fitHi = std::min(15.0, std::max(10.0, findCrossingX(hRatio, floorVal + 0.96 * (plateauVal - floorVal), 3.0, 18.0) + 3.0));
  TFitResultPtr fitRes = hRatio->Fit(f.get(), "RQSI0", "", fitLo, fitHi);
  int fitStatus = static_cast<int>(fitRes);

  if (fitStatus == 0)
  {
    const double x0 = f->GetParameter(2);
    const double slope = std::fabs(f->GetParameter(1));
    const double refinedTailLo = std::max(12.0, std::min(17.0, x0 + 4.0 / std::max(0.2, slope)));
    plateauVal = fitPlateauConst(hRatio, refinedTailLo);
    f->SetParameter(0, plateauVal);
    f->FixParameter(0, plateauVal);
    fitRes = hRatio->Fit(f.get(), "RQSI0", "", fitLo, fitHi);
    fitStatus = static_cast<int>(fitRes);
  }

  double x90 = -1.0;
  if (fitStatus == 0)
  {
    plateauVal = f->GetParameter(0);
    const double floor = f->GetParameter(3);
    x90 = f->GetX(floor + 0.90 * (plateauVal - floor), 0.0, 20.0);
  }
  if (!std::isfinite(x90) || x90 < 0.0 || x90 > 20.0)
  {
    x90 = findCrossingX(hRatio, 0.90, 2.0, 20.0);
  }

  fitOut = std::move(f);
  return x90;
}

void saveCanvas(TCanvas& c, const std::string& path)
{
  c.SaveAs(path.c_str());
  std::cout << "[WROTE] " << path << std::endl;
}
}

int PlotScaledTriggerStudy(const char* inputPathArg = "", const char* outDirArg = "")
{
  const std::string cfg =
    "jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_"
    "preselectionReference_tightReference_nonTightReference_scaledTriggerStudy";

  const std::string inputPath = inputPathArg && *inputPathArg
    ? inputPathArg
    : envOrDefault("RJ_SCALED_TRIGGER_INPUT",
                   "InputFiles/auau25/RecoilJets_auau_ALL_" + cfg + ".root");

  const std::string outDir = outDirArg && *outDirArg
    ? outDirArg
    : envOrDefault("RJ_SCALED_TRIGGER_PLOT_DIR",
                   "dataOutput/auau/" + cfg + "/scaledTrigQA");

  gSystem->mkdir(outDir.c_str(), true);
  gStyle->SetOptStat(0);

  TFile f(inputPath.c_str(), "READ");
  if (!f.IsOpen() || f.IsZombie())
  {
    std::cerr << "[ERROR] Could not open " << inputPath << std::endl;
    return 1;
  }

  std::unique_ptr<TH1> hBase(loadHist(f, "MBD_NS_geq_2_vtx_lt_150"));
  std::unique_ptr<TH1> hPho10(loadHist(f, "Photon_10"));
  std::unique_ptr<TH1> hPho12(loadHist(f, "Photon_12"));
  if (!hBase || !hPho10 || !hPho12)
  {
    std::cerr << "[ERROR] Missing one or more perRunCorrected scaled-trigger histograms in "
              << inputPath << std::endl;
    return 2;
  }

  constexpr int colorPho10 = kBlue + 1;
  constexpr int colorPho12 = kRed + 1;

  styleHist(hBase.get(), kBlack, 20);
  styleHist(hPho10.get(), colorPho10, 23);
  styleHist(hPho12.get(), colorPho12, 24);

  {
    TCanvas c("c_scaledTrigQA_groupOverlay", "c_scaledTrigQA_groupOverlay", 1200, 800);
    c.cd();
    c.SetLeftMargin(0.14);
    c.SetRightMargin(0.05);
    c.SetBottomMargin(0.14);
    c.SetTopMargin(0.08);
    gPad->SetLogy();
    gPad->SetTicks(1, 1);

    const double maxY = std::max(hBase->GetMaximum(),
                                 std::max(hPho10->GetMaximum(), hPho12->GetMaximum()));
    double minPos = smallestPositiveBinContent(hBase.get());
    const double p10Min = smallestPositiveBinContent(hPho10.get());
    const double p12Min = smallestPositiveBinContent(hPho12.get());
    if (p10Min > 0.0) minPos = (minPos > 0.0) ? std::min(minPos, p10Min) : p10Min;
    if (p12Min > 0.0) minPos = (minPos > 0.0) ? std::min(minPos, p12Min) : p12Min;

    hBase->SetTitle("");
    hBase->GetXaxis()->SetTitle("Max cluster energy [GeV], E_{clus} > 1 GeV");
    hBase->GetYaxis()->SetTitle("Scaled counts");
    hBase->GetXaxis()->SetTitleSize(0.052);
    hBase->GetYaxis()->SetTitleSize(0.052);
    hBase->GetXaxis()->SetTitleOffset(1.08);
    hBase->GetYaxis()->SetTitleOffset(1.10);
    hBase->GetXaxis()->SetLabelSize(0.043);
    hBase->GetYaxis()->SetLabelSize(0.043);
    hBase->GetXaxis()->SetRangeUser(1.0, 20.0);
    hBase->SetMinimum((minPos > 0.0) ? std::max(1.0, 0.5 * minPos) : 1.0);
    hBase->SetMaximum(std::max(10.0, 1.6 * maxY));
    hBase->Draw("hist");
    hPho10->Draw("hist same");
    hPho12->Draw("hist same");

    TLegend leg(0.54, 0.67, 0.90, 0.86);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.034);
    leg.AddEntry(hBase.get(), "MBD N&S #geq 2, |v_{z}| < 150 cm", "l");
    leg.AddEntry(hPho10.get(), "Photon 10, scaled", "l");
    leg.AddEntry(hPho12.get(), "Photon 12, scaled", "l");
    leg.Draw();

    TLatex t;
    t.SetNDC(true);
    t.SetTextFont(42);
    t.SetTextSize(0.036);
    t.DrawLatex(0.57, 0.57, "#it{#bf{sPHENIX}} Internal");
    t.DrawLatex(0.57, 0.52, "Au+Au #sqrt{s_{NN}} = 200 GeV");
    t.DrawLatex(0.57, 0.47, "Scaled-trigger study, 620 runs");

    saveCanvas(c, joinPath(outDir, "hMaxClusterEnergy_groupOverlay.png"));
    saveCanvas(c, joinPath(outDir, "scaledTriggerStudy_maxClusterEnergy_scaledCounts.png"));
  }

  std::unique_ptr<TH1> hPho10Eff(dynamic_cast<TH1*>(hPho10->Clone("hPho10Eff_scaledTrigQA")));
  std::unique_ptr<TH1> hPho12Eff(dynamic_cast<TH1*>(hPho12->Clone("hPho12Eff_scaledTrigQA")));
  if (hPho10Eff) hPho10Eff->SetDirectory(nullptr);
  if (hPho12Eff) hPho12Eff->SetDirectory(nullptr);
  hPho10Eff->Sumw2();
  hPho12Eff->Sumw2();
  hBase->Sumw2();
  hPho10Eff->Divide(hPho10.get(), hBase.get(), 1.0, 1.0, "B");
  hPho12Eff->Divide(hPho12.get(), hBase.get(), 1.0, 1.0, "B");
  styleHist(hPho10Eff.get(), colorPho10, 23);
  styleHist(hPho12Eff.get(), colorPho12, 24);

  double plateau10 = 1.0;
  double plateau12 = 1.0;
  std::unique_ptr<TF1> fit10;
  std::unique_ptr<TF1> fit12;
  const double x90_10 = fitTurnOnAndGetX90(hPho10Eff.get(), 10, colorPho10, plateau10, fit10);
  const double x90_12 = fitTurnOnAndGetX90(hPho12Eff.get(), 12, colorPho12, plateau12, fit12);

  {
    TCanvas c("c_scaledTrigQA_groupTurnOnOverlay", "c_scaledTrigQA_groupTurnOnOverlay", 1200, 800);
    c.cd();
    c.SetLeftMargin(0.14);
    c.SetRightMargin(0.05);
    c.SetBottomMargin(0.14);
    c.SetTopMargin(0.08);
    gPad->SetTicks(1, 1);

    const double maxY = std::max(hPho10Eff->GetMaximum(), hPho12Eff->GetMaximum());
    hPho10Eff->SetTitle("");
    hPho10Eff->GetXaxis()->SetTitle("Max cluster energy [GeV], E_{clus} > 1 GeV");
    hPho10Eff->GetYaxis()->SetTitle("Trigger / MBD efficiency");
    hPho10Eff->GetXaxis()->SetTitleSize(0.052);
    hPho10Eff->GetYaxis()->SetTitleSize(0.052);
    hPho10Eff->GetXaxis()->SetTitleOffset(1.08);
    hPho10Eff->GetYaxis()->SetTitleOffset(1.10);
    hPho10Eff->GetXaxis()->SetLabelSize(0.043);
    hPho10Eff->GetYaxis()->SetLabelSize(0.043);
    hPho10Eff->GetXaxis()->SetRangeUser(1.0, 20.0);
    hPho10Eff->SetMinimum(0.0);
    hPho10Eff->SetMaximum(std::max(1.22, 1.15 * maxY));
    hPho10Eff->Draw("E1");
    hPho12Eff->Draw("E1 same");

    if (fit10) fit10->Draw("same");
    if (fit12) fit12->Draw("same");

    TLine l1(1.0, 1.0, 20.0, 1.0);
    l1.SetLineColor(kGray + 2);
    l1.SetLineStyle(3);
    l1.SetLineWidth(2);
    l1.Draw("same");

    std::unique_ptr<TLine> v10Line;
    std::unique_ptr<TLine> v12Line;
    if (x90_10 > 0.0 && x90_10 <= 20.0)
    {
      const double y90 = fit10 ? fit10->Eval(x90_10) : 0.90;
      v10Line.reset(new TLine(x90_10, 0.0, x90_10, y90));
      v10Line->SetLineColor(colorPho10);
      v10Line->SetLineStyle(2);
      v10Line->SetLineWidth(2);
      v10Line->Draw("same");
    }
    if (x90_12 > 0.0 && x90_12 <= 20.0)
    {
      const double y90 = fit12 ? fit12->Eval(x90_12) : 0.90;
      v12Line.reset(new TLine(x90_12, 0.0, x90_12, y90));
      v12Line->SetLineColor(colorPho12);
      v12Line->SetLineStyle(2);
      v12Line->SetLineWidth(2);
      v12Line->Draw("same");
    }

    const std::string leg10 = (x90_10 > 0.0)
      ? Form("Photon 10 / MBD, 90%%=%.2f", x90_10)
      : std::string("Photon 10 / MBD");
    const std::string leg12 = (x90_12 > 0.0)
      ? Form("Photon 12 / MBD, 90%%=%.2f", x90_12)
      : std::string("Photon 12 / MBD");

    TLegend leg(0.56, 0.19, 0.93, 0.31);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.030);
    leg.AddEntry(hPho10Eff.get(), leg10.c_str(), "pe");
    leg.AddEntry(hPho12Eff.get(), leg12.c_str(), "pe");
    leg.Draw();

    TLatex t;
    t.SetNDC(true);
    t.SetTextFont(42);
    t.SetTextSize(0.030);
    t.SetTextAlign(13);
    t.DrawLatex(0.17, 0.91, "#it{#bf{sPHENIX}} Internal");
    t.DrawLatex(0.17, 0.865, "Au+Au #sqrt{s_{NN}} = 200 GeV");
    t.DrawLatex(0.17, 0.82, "Scaled-trigger study, 620 runs");

    t.SetTextSize(0.030);
    t.SetTextColor(colorPho10);
    t.DrawLatex(0.17, 0.675, Form("Photon 10 plateau = %.3f", plateau10));
    t.SetTextColor(colorPho12);
    t.DrawLatex(0.17, 0.635, Form("Photon 12 plateau = %.3f", plateau12));
    t.SetTextColor(kBlack);

    t.SetTextAlign(33);
    t.SetTextSize(0.030);
    t.DrawLatex(0.90, 0.51, "Fit: generalized sigmoid");
    t.DrawLatex(0.90, 0.465, "floor + plateau turn-on");

    saveCanvas(c, joinPath(outDir, "hMaxClusterEnergy_groupTurnOnOverlay.png"));
    saveCanvas(c, joinPath(outDir, "scaledTriggerStudy_triggerEfficiency.png"));
  }

  return 0;
}
