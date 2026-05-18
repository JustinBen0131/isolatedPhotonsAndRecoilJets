#include "AnalyzeRecoilJets.h"

#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1F.h"
#include "TKey.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

namespace
{
struct PtIsoHist
{
    std::string name;
    int lo = 0;
    int hi = 0;
    double center = 0.0;
    std::unique_ptr<TH1> hist;
};

struct FitResult
{
    double eff = 0.0;
    double aGeV = 0.0;
    double bPerGeV = 0.0;
    double aErr = 0.0;
    double bErr = 0.0;
    int nPoints = 0;
};

bool FindEfficiencyCut(TH1* hIn, double eff, double& cut, double& cutErr)
{
    cut = 0.0;
    cutErr = 0.0;
    if (!hIn) return false;
    const int nb = hIn->GetNbinsX();
    const double total = hIn->Integral(1, nb);
    if (!(total > 0.0)) return false;

    double running = 0.0;
    int ibCut = nb;
    for (int ib = 1; ib <= nb; ++ib)
    {
        running += hIn->GetBinContent(ib);
        if ((running / total) >= eff)
        {
            ibCut = ib;
            break;
        }
    }

    const double binLo = hIn->GetXaxis()->GetBinLowEdge(ibCut);
    const double binHi = hIn->GetXaxis()->GetBinUpEdge(ibCut);
    const double prev = running - hIn->GetBinContent(ibCut);
    const double binC = hIn->GetBinContent(ibCut);
    const double target = eff * total;
    if (binC > 0.0)
    {
        const double frac = std::min(1.0, std::max(0.0, (target - prev) / binC));
        cut = binLo + frac * (binHi - binLo);
    }
    else
    {
        cut = 0.5 * (binLo + binHi);
    }
    cutErr = 0.5 * (binHi - binLo);
    return std::isfinite(cut);
}

double ReadSIMEventCount(TFile& f)
{
    TDirectory* top = f.GetDirectory(ARJ::kDirSIM.c_str());
    if (!top) top = &f;
    TH1* cnt = dynamic_cast<TH1*>(top->Get(("cnt_" + ARJ::kDirSIM).c_str()));
    return cnt ? cnt->GetBinContent(1) : 0.0;
}

std::map<std::string, PtIsoHist> ReadWeightedIsoHists(const std::vector<std::string>& paths,
                                                       const std::vector<double>& sigmas)
{
    std::map<std::string, PtIsoHist> out;
    const std::regex pat("^h_EisoReco_truthSigMatched_pT_([0-9]+)_([0-9]+)$");

    for (std::size_t i = 0; i < paths.size(); ++i)
    {
        TFile f(paths[i].c_str(), "READ");
        if (f.IsZombie())
        {
            std::cerr << "[ERROR] Cannot open " << paths[i] << "\n";
            continue;
        }
        const double nEvt = ReadSIMEventCount(f);
        if (!(nEvt > 0.0))
        {
            std::cerr << "[ERROR] Bad cnt_SIM in " << paths[i] << "\n";
            continue;
        }
        const double weight = sigmas[i] / nEvt;

        TDirectory* top = f.GetDirectory(ARJ::kDirSIM.c_str());
        if (!top) top = &f;
        TIter next(top->GetListOfKeys());
        while (TKey* key = dynamic_cast<TKey*>(next()))
        {
            const std::string name = key->GetName();
            std::smatch m;
            if (!std::regex_match(name, m, pat)) continue;
            TH1* h = dynamic_cast<TH1*>(top->Get(name.c_str()));
            if (!h || h->GetEntries() <= 0.0) continue;

            const std::string suffix = "_pT_" + m[1].str() + "_" + m[2].str();
            const int lo = std::stoi(m[1].str());
            const int hi = std::stoi(m[2].str());
            std::unique_ptr<TH1> clone(dynamic_cast<TH1*>(h->Clone(("h_merge_" + suffix + "_" + std::to_string(i)).c_str())));
            if (!clone) continue;
            clone->SetDirectory(nullptr);
            if (clone->GetSumw2N() == 0) clone->Sumw2();
            clone->Scale(weight);

            auto& slot = out[suffix];
            if (!slot.hist)
            {
                slot.name = name;
                slot.lo = lo;
                slot.hi = hi;
                slot.center = 0.5 * (lo + hi);
                slot.hist.reset(dynamic_cast<TH1*>(clone->Clone(("h_stitched" + suffix).c_str())));
                slot.hist->SetDirectory(nullptr);
            }
            else
            {
                slot.hist->Add(clone.get());
            }
        }
    }

    return out;
}

void StyleGraph(TGraphErrors& g, int color, int marker)
{
    g.SetLineColor(color);
    g.SetMarkerColor(color);
    g.SetLineWidth(2);
    g.SetMarkerStyle(marker);
    g.SetMarkerSize(1.25);
}

void WriteCoefficients(const std::string& outDir,
                       const std::string& coneLabel,
                       const std::vector<FitResult>& fits,
                       const std::vector<std::string>& inputs,
                       const std::vector<std::string>& histNames)
{
    const std::string csvPath = outDir + "/pp_iso_wp_" + coneLabel + "_coefficients.csv";
    const std::string jsonPath = outDir + "/pp_iso_wp_" + coneLabel + "_coefficients.json";
    {
        std::ofstream csv(csvPath);
        csv << "cone,efficiency,aGeV,bPerGeV,aErr,bErr,nPoints,sideGapGeV,productionDefault\n";
        for (const auto& f : fits)
        {
            csv << coneLabel << "," << f.eff << ","
                << f.aGeV << "," << f.bPerGeV << ","
                << f.aErr << "," << f.bErr << ","
                << f.nPoints << ",1.0," << (std::fabs(f.eff - 0.90) < 1e-9 ? "true" : "false") << "\n";
        }
    }
    {
        std::ofstream js(jsonPath);
        js << "{\n";
        js << "  \"cone\": \"" << coneLabel << "\",\n";
        js << "  \"sideGapGeV\": 1.0,\n";
        js << "  \"productionDefaultEfficiency\": 0.90,\n";
        js << "  \"fitModel\": \"thrReco(pT) = aGeV + bPerGeV * pT\",\n";
        js << "  \"inputs\": [\n";
        for (std::size_t i = 0; i < inputs.size(); ++i)
            js << "    \"" << inputs[i] << "\"" << (i + 1 < inputs.size() ? "," : "") << "\n";
        js << "  ],\n";
        js << "  \"histogramNames\": [\n";
        for (std::size_t i = 0; i < histNames.size(); ++i)
            js << "    \"" << histNames[i] << "\"" << (i + 1 < histNames.size() ? "," : "") << "\n";
        js << "  ],\n";
        js << "  \"fits\": [\n";
        for (std::size_t i = 0; i < fits.size(); ++i)
        {
            const auto& f = fits[i];
            js << "    {\"efficiency\": " << f.eff
               << ", \"aGeV\": " << f.aGeV
               << ", \"bPerGeV\": " << f.bPerGeV
               << ", \"aErr\": " << f.aErr
               << ", \"bErr\": " << f.bErr
               << ", \"nPoints\": " << f.nPoints
               << ", \"productionDefault\": " << (std::fabs(f.eff - 0.90) < 1e-9 ? "true" : "false")
               << "}" << (i + 1 < fits.size() ? "," : "") << "\n";
        }
        js << "  ]\n";
        js << "}\n";
    }
    std::cout << "[WROTE] " << csvPath << "\n";
    std::cout << "[WROTE] " << jsonPath << "\n";
}

std::vector<FitResult> MakeConeFit(const std::string& coneLabel,
                                   const std::vector<std::string>& inputs,
                                   const std::string& outDir)
{
    const std::vector<double> sigmas = {ARJ::kSigmaPhoton5_pb, ARJ::kSigmaPhoton10_pb, ARJ::kSigmaPhoton20_pb};
    auto hists = ReadWeightedIsoHists(inputs, sigmas);
    if (hists.empty())
    {
        std::cerr << "[ERROR] No usable h_EisoReco_truthSigMatched_pT_* histograms for " << coneLabel << "\n";
        return {};
    }

    std::vector<double> x, ex, y70, y80, y90, ey70, ey80, ey90;
    std::vector<std::string> histNames;
    double yMin = std::numeric_limits<double>::max();
    double yMax = -std::numeric_limits<double>::max();

    for (auto& kv : hists)
    {
        auto& slot = kv.second;
        if (!slot.hist || slot.hist->Integral() <= 0.0) continue;
        double c70 = 0.0, c80 = 0.0, c90 = 0.0;
        double e70 = 0.0, e80 = 0.0, e90 = 0.0;
        if (!FindEfficiencyCut(slot.hist.get(), 0.70, c70, e70)) continue;
        if (!FindEfficiencyCut(slot.hist.get(), 0.80, c80, e80)) continue;
        if (!FindEfficiencyCut(slot.hist.get(), 0.90, c90, e90)) continue;
        x.push_back(slot.center);
        ex.push_back(0.0);
        y70.push_back(c70); ey70.push_back(e70);
        y80.push_back(c80); ey80.push_back(e80);
        y90.push_back(c90); ey90.push_back(e90);
        histNames.push_back(slot.name);
        yMin = std::min(yMin, std::min(c70 - e70, std::min(c80 - e80, c90 - e90)));
        yMax = std::max(yMax, std::max(c70 + e70, std::max(c80 + e80, c90 + e90)));
    }

    if (x.size() < 2)
    {
        std::cerr << "[ERROR] Not enough pT bins for " << coneLabel << "\n";
        return {};
    }

    TGraphErrors g90((int)x.size(), x.data(), y90.data(), ex.data(), ey90.data());
    TGraphErrors g80((int)x.size(), x.data(), y80.data(), ex.data(), ey80.data());
    TGraphErrors g70((int)x.size(), x.data(), y70.data(), ex.data(), ey70.data());
    StyleGraph(g90, kMagenta + 1, 20);
    StyleGraph(g80, kGreen + 2, 21);
    StyleGraph(g70, kBlue + 1, 22);

    const double fitXLo = *std::min_element(x.begin(), x.end()) - 0.5;
    const double fitXHi = *std::max_element(x.begin(), x.end()) + 0.5;
    TF1 f90(("f90_" + coneLabel).c_str(), "pol1", fitXLo, fitXHi);
    TF1 f80(("f80_" + coneLabel).c_str(), "pol1", fitXLo, fitXHi);
    TF1 f70(("f70_" + coneLabel).c_str(), "pol1", fitXLo, fitXHi);
    f90.SetLineColor(kMagenta + 1); f90.SetLineWidth(3);
    f80.SetLineColor(kGreen + 2); f80.SetLineWidth(3);
    f70.SetLineColor(kBlue + 1); f70.SetLineWidth(3);
    g90.Fit(&f90, "Q0");
    g80.Fit(&f80, "Q0");
    g70.Fit(&f70, "Q0");

    std::vector<FitResult> fits = {
        {0.70, f70.GetParameter(0), f70.GetParameter(1), f70.GetParError(0), f70.GetParError(1), g70.GetN()},
        {0.80, f80.GetParameter(0), f80.GetParameter(1), f80.GetParError(0), f80.GetParError(1), g80.GetN()},
        {0.90, f90.GetParameter(0), f90.GetParameter(1), f90.GetParError(0), f90.GetParError(1), g90.GetN()}
    };

    TCanvas c(("c_pp_iso_wp_" + coneLabel).c_str(), "", 900, 700);
    c.SetLeftMargin(0.14);
    c.SetRightMargin(0.05);
    c.SetBottomMargin(0.12);
    c.SetTopMargin(0.08);
    c.SetTicks(1, 1);

    const double pad = (std::isfinite(yMin) && std::isfinite(yMax) && yMax > yMin) ? 0.25 * (yMax - yMin) : 0.5;
    TH1F frame(("hFr_" + coneLabel).c_str(), "", 100, fitXLo, fitXHi);
    frame.SetDirectory(nullptr);
    frame.SetStats(false);
    frame.SetMinimum(std::max(0.0, yMin - pad));
    frame.SetMaximum(yMax + pad);
    frame.GetXaxis()->SetTitle("Truth-matched cluster E_{T} [GeV]");
    frame.GetYaxis()->SetTitle("Reco E_{T}^{iso} cutoff [GeV]");
    frame.GetXaxis()->SetTitleSize(0.055);
    frame.GetYaxis()->SetTitleSize(0.055);
    frame.GetXaxis()->SetLabelSize(0.045);
    frame.GetYaxis()->SetLabelSize(0.045);
    frame.GetYaxis()->SetTitleOffset(1.12);
    frame.Draw();

    g90.Draw("PE1 SAME");
    g80.Draw("PE1 SAME");
    g70.Draw("PE1 SAME");
    f90.Draw("SAME");
    f80.Draw("SAME");
    f70.Draw("SAME");

    TLatex info;
    info.SetNDC(true);
    info.SetTextFont(42);
    info.SetTextAlign(13);
    info.SetTextSize(0.030);
    const std::string coneText = (coneLabel == "r30") ? "R = 0.3" : "R = 0.4";
    info.DrawLatex(0.17, 0.88, ("pp PhotonJet5+10+20, " + coneText + " isolation").c_str());
    info.DrawLatex(0.17, 0.83, "Stitched by photon-sample cross section / event count");
    info.DrawLatex(0.17, 0.78, TString::Format("90%%: E_{T}^{iso} = %.3f + %.4f E_{T}", f90.GetParameter(0), f90.GetParameter(1)).Data());
    info.DrawLatex(0.17, 0.73, TString::Format("80%%: E_{T}^{iso} = %.3f + %.4f E_{T}", f80.GetParameter(0), f80.GetParameter(1)).Data());
    info.DrawLatex(0.17, 0.68, TString::Format("70%%: E_{T}^{iso} = %.3f + %.4f E_{T}", f70.GetParameter(0), f70.GetParameter(1)).Data());

    TLatex sph;
    sph.SetNDC(true);
    sph.SetTextFont(42);
    sph.SetTextAlign(33);
    sph.SetTextSize(0.040);
    sph.DrawLatex(0.92, 0.88, "#it{#bf{sPHENIX}} Internal");
    sph.SetTextSize(0.032);
    sph.DrawLatex(0.92, 0.83, "Pythia #sqrt{s}=200 GeV");

    const std::string pngPath = outDir + "/pp_iso_wp_" + coneLabel + "_eff70_80_90.png";
    c.SaveAs(pngPath.c_str());
    std::cout << "[WROTE] " << pngPath << "\n";

    WriteCoefficients(outDir, coneLabel, fits, inputs, histNames);
    return fits;
}
}

void MakePPSlidingIsoRDepFits(const char* inputDir = "dataOutput/ppSlidingIsoFits/input_roots_20260517_rdepiso",
                              const char* outDir = "dataOutput/ppSlidingIsoFits/pp_photonjet5_10_20_rdepiso_20260517")
{
    gSystem->mkdir(outDir, true);
    const std::string dir = inputDir;
    const std::vector<std::string> r30Inputs = {
        dir + "/RecoilJets_photonjet5_ALL_jetMinPt5_7pi_8_vz60_isoR30_fixedIso2GeV.root",
        dir + "/RecoilJets_photonjet10_ALL_jetMinPt5_7pi_8_vz60_isoR30_fixedIso2GeV.root",
        dir + "/RecoilJets_photonjet20_ALL_jetMinPt5_7pi_8_vz60_isoR30_fixedIso2GeV.root"
    };
    const std::vector<std::string> r40Inputs = {
        dir + "/RecoilJets_photonjet5_ALL_jetMinPt5_7pi_8_vz60_isoR40_fixedIso2GeV.root",
        dir + "/RecoilJets_photonjet10_ALL_jetMinPt5_7pi_8_vz60_isoR40_fixedIso2GeV.root",
        dir + "/RecoilJets_photonjet20_ALL_jetMinPt5_7pi_8_vz60_isoR40_fixedIso2GeV.root"
    };

    const auto r30 = MakeConeFit("r30", r30Inputs, outDir);
    const auto r40 = MakeConeFit("r40", r40Inputs, outDir);

    std::ofstream prov(std::string(outDir) + "/provenance.txt");
    prov << "Task: pp R-dependent sliding isolation fits for non-embedded PhotonJet5+10+20 SIM\n";
    prov << "Source: canonical pp SIM ALL files with jetMinPt5_7pi_8_vz60 and fixedIso2GeV, separately for isoR30 and isoR40\n";
    prov << "Histogram discovery regex: ^h_EisoReco_truthSigMatched_pT_([0-9]+)_([0-9]+)$ inside SIM directory\n";
    prov << "Stitch weights: photon sample cross section divided by cnt_SIM event count; constants from AnalyzeRecoilJets.h\n";
    prov << "Fit model: reco isolation threshold = aGeV + bPerGeV * photon ET\n";
    prov << "Production default: 90% signal-retention fit for each cone. sideGapGeV kept at 1.0 unless explicitly changed.\n";
    prov << "\nR30 inputs:\n";
    for (const auto& p : r30Inputs) prov << "  " << p << "\n";
    prov << "\nR40 inputs:\n";
    for (const auto& p : r40Inputs) prov << "  " << p << "\n";
    prov << "\nR30 fit count: " << r30.size() << "\n";
    prov << "R40 fit count: " << r40.size() << "\n";
    prov.close();
    std::cout << "[WROTE] " << std::string(outDir) + "/provenance.txt" << "\n";
}
