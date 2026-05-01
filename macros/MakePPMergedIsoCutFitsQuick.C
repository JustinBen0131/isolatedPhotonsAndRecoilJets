#include "AnalyzeRecoilJets.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace
{
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

void StyleGraph(TGraphErrors& g, int color, int marker, double markerSize)
{
    g.SetLineWidth(2);
    g.SetLineColor(color);
    g.SetMarkerColor(color);
    g.SetMarkerStyle(marker);
    g.SetMarkerSize(markerSize);
}
}

void MakePPMergedIsoCutFitsQuick()
{
    const std::string tag =
        "jetMinPt5_7pi_8_vz60_isoR40_fixedIso2GeV_"
        "preselectionVariantA_tightVariantA_nonTightVariantA";
    const std::string inPath =
        ARJ::kOutputBase + "/combinedSimOnly/" + tag +
        "/photonJet5and10and20merged_SIM/RecoilJets_photonjet5plus10plus20_MERGED.root";
    const std::string outPath =
        ARJ::kOutputBase + "/pp/ppg12Style_isoCutEfficiencyFits_ppMerged_fixedIso2GeV_BDTNPB.png";

    TFile f(inPath.c_str(), "READ");
    if (f.IsZombie())
    {
        std::cerr << "[ERROR] Cannot open merged SIM input: " << inPath << "\n";
        return;
    }

    TDirectory* top = f.GetDirectory(ARJ::kDirSIM.c_str());
    if (!top) top = &f;

    std::vector<double> x, ex;
    std::vector<double> y90, ey90, y80, ey80, y70, ey70;
    double yMin = std::numeric_limits<double>::max();
    double yMax = -std::numeric_limits<double>::max();

    for (int ipt = 1; ipt < ARJ::kNPtBins; ++ipt)
    {
        const ARJ::PtBin& b = ARJ::PtBins()[ipt];
        TH1* hSig = dynamic_cast<TH1*>(top->Get(("h_EisoReco_truthSigMatched" + b.suffix).c_str()));
        if (!hSig || hSig->GetEntries() <= 0.0) continue;

        double c70 = 0.0, e70 = 0.0;
        double c80 = 0.0, e80 = 0.0;
        double c90 = 0.0, e90 = 0.0;
        if (!FindEfficiencyCut(hSig, 0.70, c70, e70)) continue;
        if (!FindEfficiencyCut(hSig, 0.80, c80, e80)) continue;
        if (!FindEfficiencyCut(hSig, 0.90, c90, e90)) continue;

        x.push_back(0.5 * (ARJ::kPtEdges[(std::size_t)ipt] + ARJ::kPtEdges[(std::size_t)ipt + 1]));
        ex.push_back(0.0);
        y90.push_back(c90); ey90.push_back(e90);
        y80.push_back(c80); ey80.push_back(e80);
        y70.push_back(c70); ey70.push_back(e70);

        yMin = std::min(yMin, std::min(c70 - e70, std::min(c80 - e80, c90 - e90)));
        yMax = std::max(yMax, std::max(c70 + e70, std::max(c80 + e80, c90 + e90)));
    }

    if (x.empty())
    {
        std::cerr << "[ERROR] No usable h_EisoReco_truthSigMatched pT bins in " << inPath << "\n";
        return;
    }

    const double pad = (std::isfinite(yMin) && std::isfinite(yMax) && yMax > yMin)
        ? 1.2 * (yMax - yMin)
        : 0.25;

    TGraphErrors g90((int)x.size(), x.data(), y90.data(), ex.data(), ey90.data());
    TGraphErrors g80((int)x.size(), x.data(), y80.data(), ex.data(), ey80.data());
    TGraphErrors g70((int)x.size(), x.data(), y70.data(), ex.data(), ey70.data());
    StyleGraph(g90, kMagenta + 1, 20, 1.4);
    StyleGraph(g80, kGreen + 2, 21, 1.4);
    StyleGraph(g70, kBlue + 1, 22, 1.5);

    const double fitXLo = ARJ::kPtEdges[1];
    const double fitXHi = ARJ::kPtEdges.back();
    TF1 f90("f90_ppMerged_fixedIso2_BDTNPB", "pol1", fitXLo, fitXHi);
    TF1 f80("f80_ppMerged_fixedIso2_BDTNPB", "pol1", fitXLo, fitXHi);
    TF1 f70("f70_ppMerged_fixedIso2_BDTNPB", "pol1", fitXLo, fitXHi);
    f90.SetLineColor(kMagenta + 1); f90.SetLineWidth(3);
    f80.SetLineColor(kGreen + 2);   f80.SetLineWidth(3);
    f70.SetLineColor(kBlue + 1);    f70.SetLineWidth(3);

    TCanvas c("c_ppg12IsoCutFits_ppMerged_fixedIso2_BDTNPB", "", 900, 700);
    c.SetLeftMargin(0.14);
    c.SetRightMargin(0.05);
    c.SetBottomMargin(0.12);
    c.SetTopMargin(0.08);
    c.SetTicks(1, 1);

    TH1F frame("hFr_ppMerged_fixedIso2_BDTNPB", "", 100, ARJ::kPtEdges[1], ARJ::kPtEdges[ARJ::kNPtBins]);
    frame.SetDirectory(nullptr);
    frame.SetStats(0);
    frame.SetMinimum(std::max(0.0, yMin - pad));
    frame.SetMaximum(yMax + pad);
    frame.GetXaxis()->SetTitle("Cluster p_{T} [GeV]");
    frame.GetYaxis()->SetTitle("E_{T}^{iso} Cutoff [GeV]");
    frame.GetXaxis()->SetTitleSize(0.060);
    frame.GetYaxis()->SetTitleSize(0.060);
    frame.GetXaxis()->SetLabelSize(0.050);
    frame.GetYaxis()->SetLabelSize(0.050);
    frame.GetYaxis()->SetTitleOffset(1.05);
    frame.Draw();

    g90.Draw("PE1 SAME");
    g80.Draw("PE1 SAME");
    g70.Draw("PE1 SAME");
    if (g90.GetN() >= 2) { g90.Fit(&f90, "Q0"); f90.Draw("SAME"); }
    if (g80.GetN() >= 2) { g80.Fit(&f80, "Q0"); f80.Draw("SAME"); }
    if (g70.GetN() >= 2) { g70.Fit(&f70, "Q0"); f70.Draw("SAME"); }

    TLegend leg(0.50, 0.62, 0.92, 0.78);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextFont(42);
    leg.SetTextSize(0.030);
    leg.AddEntry(&g90, "90% Efficiency", "ep");
    leg.AddEntry(&g80, "80% Efficiency", "ep");
    leg.AddEntry(&g70, "70% Efficiency", "ep");
    leg.Draw();

    TLatex info;
    info.SetNDC(true);
    info.SetTextFont(42);
    info.SetTextAlign(13);
    info.SetTextSize(0.030);
    info.DrawLatex(0.18, 0.88, "Pythia 5+10+20 merged, #sqrt{s} = 200 GeV");
    if (g90.GetN() >= 2)
        info.DrawLatex(0.18, 0.83,
                       TString::Format("90%%: E_{T}^{iso} = %.3f  + %.3fp_{T}",
                                       f90.GetParameter(0), f90.GetParameter(1)).Data());
    if (g80.GetN() >= 2)
        info.DrawLatex(0.18, 0.78,
                       TString::Format("80%%: E_{T}^{iso} = %.3f  + %.3fp_{T}",
                                       f80.GetParameter(0), f80.GetParameter(1)).Data());
    if (g70.GetN() >= 2)
        info.DrawLatex(0.18, 0.73,
                       TString::Format("70%%: E_{T}^{iso} = %.3f  + %.3fp_{T}",
                                       f70.GetParameter(0), f70.GetParameter(1)).Data());
    info.DrawLatex(0.18, 0.68, "|v_{z}| < 60 cm, |#eta^{#gamma}| < 0.7, #DeltaR_{cone} < 0.4");

    TLatex sph;
    sph.SetNDC(true);
    sph.SetTextFont(42);
    sph.SetTextAlign(33);
    sph.SetTextSize(0.042);
    sph.DrawLatex(0.92, 0.88, "#bf{sPHENIX} #it{Internal}");
    sph.SetTextSize(0.034);
    sph.DrawLatex(0.92, 0.83, "Pythia #sqrt{s} = 200 GeV");

    gSystem->mkdir((ARJ::kOutputBase + "/pp").c_str(), true);
    c.SaveAs(outPath.c_str());
    std::cout << "[WROTE] " << outPath << "\n";
}
