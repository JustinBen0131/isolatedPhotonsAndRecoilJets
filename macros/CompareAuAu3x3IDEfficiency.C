#include "sPhenixStyle.C"

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace
{
const std::string kOutDir = "dataOutput/auau_bdt3x3_mc_validation/id_efficiency_compare";
const std::string kFullCfg = "preselectionNewPPG12_tightAuAuCentInputBDT_nonTightAuAuBDTComplement_baseVariant";
const std::string k3x3Cfg = "preselectionNewPPG12_tightAuAuCentInput3x3BDT_nonTightAuAuBDTComplement_baseVariant";

struct PtBin
{
  int lo = 0;
  int hi = 0;
};

struct CentBin
{
  std::string suffix;
  std::string label;
};

struct SampleSpec
{
  std::string key;
  std::string label;
  std::string sampleType;
  std::string path;
  int color = kBlack;
  int marker = 20;
  double xShift = 0.0;
};

struct Yield
{
  double tight = 0.0;
  double nonTight = 0.0;
  double tightErr = 0.0;
  double nonTightErr = 0.0;

  double total() const { return tight + nonTight; }
  double eff() const { return total() > 0.0 ? tight / total() : 0.0; }
  double effErr() const
  {
    const double den = total();
    if (den <= 0.0) return 0.0;
    const double eDen = std::hypot(tightErr, nonTightErr);
    const double relT = tight > 0.0 ? tightErr / tight : 0.0;
    const double relD = eDen / den;
    const double err = eff() * std::hypot(relT, relD);
    return std::isfinite(err) ? err : 0.0;
  }
};

const std::vector<PtBin> kPtBins = {
    {5, 8}, {8, 10}, {10, 12}, {12, 14}, {14, 16}, {16, 18},
    {18, 20}, {20, 22}, {22, 24}, {24, 26}, {26, 35}};

const std::vector<CentBin> kCentBins = {
    {"0_20", "0-20% central"},
    {"20_50", "20-50% mid-central"},
    {"50_80", "50-80% peripheral"}};

std::string signalPathFull()
{
  return "dataOutput/combinedSimOnlyEMBEDDED/" + kFullCfg +
         "/photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root";
}

std::string backgroundPathFull()
{
  return "dataOutput/combinedSimOnlyEMBEDDED/" + kFullCfg +
         "/embeddedJet12and20merged_SIM/RecoilJets_embeddedJet12plus20_MERGED.root";
}

std::string signalPath3x3(const std::string& wp)
{
  return "dataOutput/auau_bdt3x3_mc_validation/" + wp +
         "/combinedSimOnlyEMBEDDED/" + k3x3Cfg +
         "/photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root";
}

std::string backgroundPath3x3(const std::string& wp)
{
  return "dataOutput/auau_bdt3x3_mc_validation/" + wp +
         "/combinedSimOnlyEMBEDDED/" + k3x3Cfg +
         "/embeddedJet12and20merged_SIM/RecoilJets_embeddedJet12plus20_MERGED.root";
}

TH1* getHist(TDirectory* d, const std::string& name)
{
  if (!d) return nullptr;
  return dynamic_cast<TH1*>(d->Get(name.c_str()));
}

double integralAndError(TH1* h, double& err)
{
  err = 0.0;
  if (!h) return 0.0;
  return h->IntegralAndError(0, h->GetNbinsX() + 1, err);
}

Yield readYield(TDirectory* d, const PtBin& pt, const CentBin& cent)
{
  const std::string tag =
      "isoR30_pT_" + std::to_string(pt.lo) + "_" + std::to_string(pt.hi) +
      "_cent_" + cent.suffix;
  Yield y;
  y.tight = integralAndError(getHist(d, "h_Eiso_tight_" + tag), y.tightErr);
  y.nonTight = integralAndError(getHist(d, "h_Eiso_nonTight_" + tag), y.nonTightErr);
  return y;
}

Yield addYield(const Yield& a, const Yield& b)
{
  Yield out;
  out.tight = a.tight + b.tight;
  out.nonTight = a.nonTight + b.nonTight;
  out.tightErr = std::hypot(a.tightErr, b.tightErr);
  out.nonTightErr = std::hypot(a.nonTightErr, b.nonTightErr);
  return out;
}

std::unique_ptr<TGraphErrors> buildGraph(TDirectory* d,
                                         const SampleSpec& sample,
                                         const CentBin& cent,
                                         const std::vector<PtBin>& bins,
                                         std::ofstream& csv,
                                         const std::map<std::string, Yield>& baselineByBin)
{
  std::unique_ptr<TGraphErrors> g(new TGraphErrors());
  g->SetName(("g_" + sample.sampleType + "_" + sample.key + "_" + cent.suffix).c_str());
  g->SetMarkerStyle(sample.marker);
  g->SetMarkerSize(1.18);
  g->SetMarkerColor(sample.color);
  g->SetLineColor(sample.color);
  g->SetLineWidth(0);

  for (const auto& pt : bins)
  {
    const Yield y = readYield(d, pt, cent);
    if (y.total() <= 0.0) continue;

    const std::string binKey =
        sample.sampleType + "|" + cent.suffix + "|" +
        std::to_string(pt.lo) + "_" + std::to_string(pt.hi);
    auto it = baselineByBin.find(binKey);
    const double baseEff = (it != baselineByBin.end()) ? it->second.eff() : 0.0;
    const double delta = y.eff() - baseEff;
    const double ratio = baseEff > 0.0 ? y.eff() / baseEff : 0.0;

    const int ip = g->GetN();
    g->SetPoint(ip, 0.5 * (pt.lo + pt.hi) + sample.xShift, y.eff());
    g->SetPointError(ip, 0.0, y.effErr());

    csv << sample.sampleType << "," << sample.key << "," << sample.label << ","
        << cent.suffix << "," << pt.lo << "," << pt.hi << ","
        << std::setprecision(12) << y.tight << "," << y.nonTight << ","
        << y.eff() << "," << y.effErr() << ","
        << baseEff << "," << delta << "," << ratio << "\n";
  }
  return g;
}

void drawHeader(double x, double y)
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(13);
  tx.SetTextSize(0.037);
  tx.DrawLatex(x, y, "#it{#bf{sPHENIX}} Internal");
  tx.SetTextSize(0.029);
  tx.DrawLatex(x, y - 0.070, "Pythia overlay, #sqrt{s_{NN}} = 200 GeV");
}

void drawPanelLabel(const std::string& text)
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(22);
  tx.SetTextSize(0.052);
  tx.DrawLatex(0.50, 0.955, text.c_str());
}

std::string sampleAxisTitle(const std::string& sampleType)
{
  if (sampleType == "signal") return "Signal tight-ID efficiency";
  return "Background tight-ID acceptance";
}

void makePlot(const std::string& sampleType,
              const std::vector<SampleSpec>& samples,
              const std::map<std::string, Yield>& baselineByBin,
              std::ofstream& csv)
{
  TCanvas c(("c_" + sampleType + "_id_eff_3x3").c_str(),
            ("c_" + sampleType + "_id_eff_3x3").c_str(), 2200, 760);
  c.Divide(3, 1, 0.012, 0.0);

  std::vector<std::unique_ptr<TFile>> files;
  std::vector<std::unique_ptr<TGraphErrors>> graphs;
  std::vector<std::unique_ptr<TH1F>> frames;
  std::vector<std::unique_ptr<TLegend>> legends;

  for (size_t ic = 0; ic < kCentBins.size(); ++ic)
  {
    c.cd(ic + 1);
    gPad->SetTicks(1, 1);
    gPad->SetGrid(0, 0);
    gPad->SetLeftMargin(ic == 0 ? 0.16 : 0.07);
    gPad->SetRightMargin(ic == kCentBins.size() - 1 ? 0.035 : 0.015);
    gPad->SetTopMargin(0.32);
    gPad->SetBottomMargin(0.14);

    std::unique_ptr<TH1F> frame(new TH1F(("frame_" + sampleType + "_" + kCentBins[ic].suffix).c_str(), "", 100, 4.6, 35.4));
    frame->SetDirectory(nullptr);
    frame->SetStats(false);
    frame->GetXaxis()->SetTitle("Photon candidate E_{T} [GeV]");
    frame->GetYaxis()->SetTitle(ic == 0 ? sampleAxisTitle(sampleType).c_str() : "");
    frame->GetXaxis()->SetTitleSize(0.044);
    frame->GetYaxis()->SetTitleSize(0.043);
    frame->GetXaxis()->SetLabelSize(0.037);
    frame->GetYaxis()->SetLabelSize(0.037);
    frame->GetYaxis()->SetTitleOffset(ic == 0 ? 1.04 : 1.00);
    if (sampleType == "signal")
    {
      frame->SetMinimum(0.35);
      frame->SetMaximum(1.03);
    }
    else
    {
      frame->SetMinimum(0.0);
      frame->SetMaximum(0.84);
    }
    frame->Draw();
    frames.push_back(std::move(frame));

    drawPanelLabel(kCentBins[ic].label);
    if (ic == 0) drawHeader(0.145, 0.870);

    std::unique_ptr<TLegend> leg;
    if (ic == 0)
      leg.reset(new TLegend(0.50, 0.690, 0.940, 0.850));
    else
      leg.reset(new TLegend(0.18, 0.715, 0.875, 0.885));
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.030);

    for (const auto& sample : samples)
    {
      std::unique_ptr<TFile> f(TFile::Open(sample.path.c_str(), "READ"));
      if (!f || f->IsZombie())
      {
        std::cerr << "[WARN] missing or zombie ROOT file: " << sample.path << "\n";
        continue;
      }
      TDirectory* d = dynamic_cast<TDirectory*>(f->Get("SIM"));
      if (!d)
      {
        std::cerr << "[WARN] missing SIM directory in: " << sample.path << "\n";
        continue;
      }
      auto g = buildGraph(d, sample, kCentBins[ic], kPtBins, csv, baselineByBin);
      if (g && g->GetN() > 0)
      {
        g->Draw("P SAME");
        leg->AddEntry(g.get(), sample.label.c_str(), "p");
        graphs.push_back(std::move(g));
      }
      files.push_back(std::move(f));
    }
    leg->Draw();
    legends.push_back(std::move(leg));
  }

  const std::string png = kOutDir + "/bdt3x3_vs_fullcluster_" + sampleType + "_tight_fraction.png";
  c.SaveAs(png.c_str());
}

void printIntegratedSummary(const std::vector<SampleSpec>& allSamples)
{
  std::cout << "\n================ 3x3 ID EFFICIENCY SUMMARY ================\n";
  std::cout << "Comparison baseline: full-cluster width variables, centrality as input, WP0.50\n";
  std::cout << "3x3 model swaps the #eta/#phi width inputs to the 3x3-window versions.\n\n";

  std::map<std::string, Yield> baselineIntegrated;

  for (const auto& sample : allSamples)
  {
    std::unique_ptr<TFile> f(TFile::Open(sample.path.c_str(), "READ"));
    if (!f || f->IsZombie())
    {
      std::cout << "[WARN] cannot open " << sample.path << "\n";
      continue;
    }
    TDirectory* d = dynamic_cast<TDirectory*>(f->Get("SIM"));
    if (!d)
    {
      std::cout << "[WARN] missing SIM directory in " << sample.path << "\n";
      continue;
    }

    Yield all;
    std::map<std::string, Yield> byCent;
    for (const auto& cent : kCentBins)
    {
      Yield centYield;
      for (const auto& pt : kPtBins)
      {
        centYield = addYield(centYield, readYield(d, pt, cent));
      }
      byCent[cent.suffix] = centYield;
      all = addYield(all, centYield);
    }

    if (sample.key == "full_wp050")
    {
      baselineIntegrated[sample.sampleType + "|all"] = all;
      for (const auto& cent : kCentBins)
      {
        baselineIntegrated[sample.sampleType + "|" + cent.suffix] = byCent[cent.suffix];
      }
    }

    auto printLine = [&](const std::string& centKey, const Yield& y) {
      const std::string baseKey = sample.sampleType + "|" + centKey;
      const double baseEff = baselineIntegrated.count(baseKey) ? baselineIntegrated[baseKey].eff() : y.eff();
      const double delta = y.eff() - baseEff;
      const double rel = baseEff > 0.0 ? 100.0 * delta / baseEff : 0.0;
      std::cout << std::setw(12) << sample.sampleType
                << " | " << std::setw(14) << sample.key
                << " | " << std::setw(6) << centKey
                << " | eff=" << std::fixed << std::setprecision(4) << y.eff()
                << " +/- " << y.effErr()
                << " | delta_vs_full_wp050=" << std::showpos << delta
                << " (" << rel << "%)" << std::noshowpos
                << " | tight=" << std::setprecision(2) << y.tight
                << " nonTight=" << y.nonTight << "\n";
    };

    printLine("all", all);
    for (const auto& cent : kCentBins)
    {
      printLine(cent.suffix, byCent[cent.suffix]);
    }
    std::cout << "\n";
  }
  std::cout << "============================================================\n";
}
}  // namespace

void CompareAuAu3x3IDEfficiency()
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  gSystem->mkdir(kOutDir.c_str(), true);

  const std::vector<SampleSpec> signalSamples = {
      {"full_wp050", "Full-cluster widths, BDT > 0.50", "signal", signalPathFull(), kBlack, 20, -0.16},
      {"3x3_wp050", "3#times3 widths, BDT > 0.50", "signal", signalPath3x3("wp050"), kAzure + 2, 20, 0.00},
      {"3x3_wp080", "3#times3 widths, BDT > 0.80", "signal", signalPath3x3("wp080"), kOrange + 7, 20, 0.16}};

  const std::vector<SampleSpec> backgroundSamples = {
      {"full_wp050", "Full-cluster widths, BDT > 0.50", "background", backgroundPathFull(), kBlack, 20, -0.16},
      {"3x3_wp050", "3#times3 widths, BDT > 0.50", "background", backgroundPath3x3("wp050"), kAzure + 2, 20, 0.00},
      {"3x3_wp080", "3#times3 widths, BDT > 0.80", "background", backgroundPath3x3("wp080"), kOrange + 7, 20, 0.16}};

  std::map<std::string, Yield> baselineByBin;
  auto collectBaseline = [&](const SampleSpec& sample) {
    std::unique_ptr<TFile> f(TFile::Open(sample.path.c_str(), "READ"));
    if (!f || f->IsZombie()) return;
    TDirectory* d = dynamic_cast<TDirectory*>(f->Get("SIM"));
    if (!d) return;
    for (const auto& cent : kCentBins)
    {
      for (const auto& pt : kPtBins)
      {
        const std::string key =
            sample.sampleType + "|" + cent.suffix + "|" +
            std::to_string(pt.lo) + "_" + std::to_string(pt.hi);
        baselineByBin[key] = readYield(d, pt, cent);
      }
    }
  };
  collectBaseline(signalSamples.front());
  collectBaseline(backgroundSamples.front());

  std::ofstream csv(kOutDir + "/bdt3x3_vs_fullcluster_id_efficiency.csv");
  csv << "sample_type,model_key,label,cent,pt_lo,pt_hi,tight_sum,nontight_sum,"
         "tight_fraction,tight_fraction_err,baseline_full_wp050,delta_vs_full_wp050,ratio_vs_full_wp050\n";

  makePlot("signal", signalSamples, baselineByBin, csv);
  makePlot("background", backgroundSamples, baselineByBin, csv);

  std::vector<SampleSpec> allSamples;
  allSamples.insert(allSamples.end(), signalSamples.begin(), signalSamples.end());
  allSamples.insert(allSamples.end(), backgroundSamples.begin(), backgroundSamples.end());
  printIntegratedSummary(allSamples);

  std::cout << "[DONE] CSV: " << kOutDir << "/bdt3x3_vs_fullcluster_id_efficiency.csv\n";
  std::cout << "[DONE] Signal plot: " << kOutDir << "/bdt3x3_vs_fullcluster_signal_tight_fraction.png\n";
  std::cout << "[DONE] Background plot: " << kOutDir << "/bdt3x3_vs_fullcluster_background_tight_fraction.png\n";
}
