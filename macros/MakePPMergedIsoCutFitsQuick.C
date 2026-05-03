#include "AnalyzeRecoilJets.h"

#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
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

TH1* CloneHistFromTopDir(TFile& f, const std::string& histName, const std::string& cloneName)
{
    TDirectory* top = f.GetDirectory(ARJ::kDirSIM.c_str());
    if (!top) top = &f;
    TH1* h = dynamic_cast<TH1*>(top->Get(histName.c_str()));
    if (!h) return nullptr;

    TH1* c = dynamic_cast<TH1*>(h->Clone(cloneName.c_str()));
    if (!c) return nullptr;
    c->SetDirectory(nullptr);
    if (c->GetSumw2N() == 0) c->Sumw2();
    return c;
}

double ReadSIMEventCount(TFile& f)
{
    TDirectory* top = f.GetDirectory(ARJ::kDirSIM.c_str());
    if (!top) top = &f;
    TH1* cnt = dynamic_cast<TH1*>(top->Get(("cnt_" + ARJ::kDirSIM).c_str()));
    return cnt ? cnt->GetBinContent(1) : 0.0;
}

std::unique_ptr<TH1> MakeSmoothReference(const TH1* h, const std::string& name)
{
    if (!h) return nullptr;
    std::unique_ptr<TH1> ref(dynamic_cast<TH1*>(h->Clone(name.c_str())));
    if (!ref) return nullptr;
    ref->SetDirectory(nullptr);
    ref->Smooth(2);
    for (int ib = 1; ib <= ref->GetNbinsX(); ++ib)
    {
        if (ref->GetBinContent(ib) <= 0.0 && h->GetBinContent(ib) > 0.0)
        {
            ref->SetBinContent(ib, h->GetBinContent(ib));
        }
    }
    return ref;
}

std::unique_ptr<TH1> MakeRatioToSmooth(const TH1* h, const TH1* smooth, const std::string& name)
{
    if (!h || !smooth) return nullptr;
    std::unique_ptr<TH1> ratio(dynamic_cast<TH1*>(h->Clone(name.c_str())));
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

void StyleHist(TH1* h, int color, int marker)
{
    if (!h) return;
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(0.85);
    h->SetLineWidth(2);
}

void DrawPPPhotonStitchTruthSpectrum(const std::string& in5,
                                     const std::string& in10,
                                     const std::string& in20,
                                     const std::string& outPath)
{
    TFile f5(in5.c_str(), "READ");
    TFile f10(in10.c_str(), "READ");
    TFile f20(in20.c_str(), "READ");
    if (f5.IsZombie() || f10.IsZombie() || f20.IsZombie())
    {
        std::cerr << "[ERROR] Cannot open one or more PP photon-jet stitch inputs.\n"
                  << "  " << in5 << "\n"
                  << "  " << in10 << "\n"
                  << "  " << in20 << "\n";
        return;
    }

    std::unique_ptr<TH1> h5(CloneHistFromTopDir(f5, "h_ppPhotonStitch_maxPhotonPt_kept", "h_ppStitch5"));
    std::unique_ptr<TH1> h10(CloneHistFromTopDir(f10, "h_ppPhotonStitch_maxPhotonPt_kept", "h_ppStitch10"));
    std::unique_ptr<TH1> h20(CloneHistFromTopDir(f20, "h_ppPhotonStitch_maxPhotonPt_kept", "h_ppStitch20"));
    if (!h5 || !h10 || !h20)
    {
        std::cerr << "[ERROR] Missing SIM/h_ppPhotonStitch_maxPhotonPt_kept in one or more inputs.\n";
        return;
    }

    h5->Rebin(2);
    h10->Rebin(2);
    h20->Rebin(2);

    const double n5 = ReadSIMEventCount(f5);
    const double n10 = ReadSIMEventCount(f10);
    const double n20 = ReadSIMEventCount(f20);
    if (n5 <= 0.0 || n10 <= 0.0 || n20 <= 0.0)
    {
        std::cerr << "[ERROR] Bad SIM event counts for pp stitching QA: "
                  << "n5=" << n5 << " n10=" << n10 << " n20=" << n20 << "\n";
        return;
    }

    const double ew5 = ARJ::kSigmaPhoton5_pb / n5;
    const double ew10 = ARJ::kSigmaPhoton10_pb / n10;
    const double ew20 = ARJ::kSigmaPhoton20_pb / n20;
    const double w5 = ew5 / ew20;
    const double w10 = ew10 / ew20;
    const double w20 = 1.0;
    h5->Scale(w5);
    h10->Scale(w10);
    h20->Scale(w20);

    StyleHist(h5.get(), kBlue + 1, 20);
    StyleHist(h10.get(), kGreen + 2, 21);
    StyleHist(h20.get(), kRed + 1, 22);

    std::unique_ptr<TH1> hSum(dynamic_cast<TH1*>(h5->Clone("h_ppStitch_sum")));
    hSum->SetDirectory(nullptr);
    hSum->Add(h10.get());
    hSum->Add(h20.get());
    StyleHist(hSum.get(), kBlack, 24);

    std::unique_ptr<TH1> hSmooth = MakeSmoothReference(hSum.get(), "h_ppStitch_smooth");
    std::unique_ptr<TH1> hRatio = MakeRatioToSmooth(hSum.get(), hSmooth.get(), "h_ppStitch_ratio");
    StyleHist(hRatio.get(), kBlack, 20);

    const double xMin = 0.0;
    const double xMax = 60.0;

    TCanvas c("c_ppPhotonStitchTruthSpectrum", "pp photon stitch truth spectrum", 1050, 850);
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
    const double ymax = std::max({h5->GetMaximum(), h10->GetMaximum(), h20->GetMaximum(), hSum->GetMaximum()});
    hSum->SetTitle("");
    hSum->GetXaxis()->SetRangeUser(xMin, xMax);
    hSum->GetXaxis()->SetLabelSize(0.0);
    hSum->GetXaxis()->SetTitleSize(0.0);
    hSum->GetYaxis()->SetTitle("#sigma/#sigma_{20} scaled entries [arb. / bin]");
    hSum->GetYaxis()->SetTitleOffset(1.12);
    hSum->SetMinimum(std::max(1.0e-8, ymax * 2.0e-5));
    hSum->SetMaximum(std::max(1.0e-6, ymax * 85.0));
    hSum->Draw("E1");
    h5->Draw("E1 SAME");
    h10->Draw("E1 SAME");
    h20->Draw("E1 SAME");
    hSmooth->SetLineColor(kGray + 2);
    hSmooth->SetLineStyle(2);
    hSmooth->SetLineWidth(2);
    hSmooth->SetMarkerSize(0);
    hSmooth->Draw("HIST SAME");

    TLegend leg(0.15, 0.15, 0.52, 0.43);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.036);
    leg.AddEntry(h5.get(), "PhotonJet5 stitched", "lep");
    leg.AddEntry(h10.get(), "PhotonJet10 stitched", "lep");
    leg.AddEntry(h20.get(), "PhotonJet20 stitched", "lep");
    leg.AddEntry(hSum.get(), "Combined", "lep");
    leg.AddEntry(hSmooth.get(), "Smoothed reference", "l");
    leg.Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.040);
    lat.DrawLatex(0.15, 0.86, "PP photon generator stitching spectrum");
    lat.SetTextSize(0.031);
    lat.DrawLatex(0.15, 0.80, "Uses SIM/h_ppPhotonStitch_maxPhotonPt_kept");
    lat.DrawLatex(0.15, 0.75, "PhotonJet5: p_{T}^{#gamma}<14; PhotonJet10: 14#leq p_{T}^{#gamma}<22; PhotonJet20: p_{T}^{#gamma}#geq22 GeV");
    lat.DrawLatex(0.15, 0.70, "relative weights: #sigma_{5}/#sigma_{20}, #sigma_{10}/#sigma_{20}, 1");

    TLatex sph;
    sph.SetNDC(true);
    sph.SetTextFont(42);
    sph.SetTextAlign(33);
    sph.SetTextSize(0.042);
    sph.DrawLatex(0.92, 0.58, "#bf{sPHENIX} #it{Internal}");
    sph.SetTextSize(0.034);
    sph.DrawLatex(0.92, 0.53, "Pythia #sqrt{s} = 200 GeV");

    bot.cd();
    hRatio->SetTitle("");
    hRatio->GetXaxis()->SetRangeUser(xMin, xMax);
    hRatio->GetXaxis()->SetTitle("max stored truth p_{T}^{#gamma} [GeV]");
    hRatio->GetYaxis()->SetTitle("sum / smooth");
    hRatio->GetYaxis()->SetRangeUser(0.85, 1.15);
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->GetYaxis()->SetTitleSize(0.085);
    hRatio->GetYaxis()->SetLabelSize(0.075);
    hRatio->GetYaxis()->SetTitleOffset(0.50);
    hRatio->GetXaxis()->SetTitleSize(0.090);
    hRatio->GetXaxis()->SetLabelSize(0.080);
    hRatio->Draw("E1");
    TLine one(xMin, 1.0, xMax, 1.0);
    one.SetLineColor(kGray + 1);
    one.SetLineStyle(2);
    one.Draw("SAME");

    c.SaveAs(outPath.c_str());
    std::cout << "[WROTE] " << outPath << "\n";
}
}

void MakePPMergedIsoCutFitsQuick()
{
    const std::string tag =
        "jetMinPt5_7pi_8_vz60_isoR40_fixedIso2GeV_"
        "preselectionReference_tightReference_nonTightReference";
    const std::string in5 = ARJ::InputSim("photonjet5", tag);
    const std::string in10 = ARJ::InputSim("photonjet10", tag);
    const std::string in20 = ARJ::InputSim("photonjet20", tag);
    const std::string mergedDir = ARJ::OutputCombinedSimOnly(tag, "photonJet5and10and20merged_SIM");
    const std::string inPath = mergedDir + "/RecoilJets_photonjet5plus10plus20_MERGED.root";
    const std::string outPath =
        ARJ::kOutputBase + "/pp/ppg12Style_isoCutEfficiencyFits_ppMerged_fixedIso2GeV_reference.png";
    const std::string stitchOutPath =
        ARJ::kOutputBase + "/pp/ppg12Style_photonJetTruthStitchSpectrum_ppMerged_fixedIso2GeV_reference.png";

    gSystem->mkdir(mergedDir.c_str(), true);
    const bool mergedOk = ARJ::BuildMergedSIMFile_PhotonSlices(
        {in5, in10, in20},
        {ARJ::kSigmaPhoton5_pb, ARJ::kSigmaPhoton10_pb, ARJ::kSigmaPhoton20_pb},
        inPath,
        ARJ::kDirSIM,
        {"photonJet5", "photonJet10", "photonJet20"});
    if (!mergedOk)
    {
        std::cerr << "[ERROR] Failed to build merged photon5+10+20 SIM file.\n";
        return;
    }

    DrawPPPhotonStitchTruthSpectrum(in5, in10, in20, stitchOutPath);

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
