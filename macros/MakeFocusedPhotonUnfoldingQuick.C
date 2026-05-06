#include "AnalyzeRecoilJets.h"

#include "RooUnfold.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

using namespace ARJ;

namespace
{
struct FocusSample
{
    std::string key;
    std::string label;
    std::string triggerLabel;
    std::string cutLabel;
    std::string cutLine1;
    std::string cutLine2;
    std::string simLabel;
    std::string dataPath;
    std::string simPath;
    std::string dataTopDir;
    std::string simTopDir;
    std::string centFolder;
    std::string centSuffix;
    std::string centLabel;
    std::string outDir;
};

struct CorrectionView
{
    std::string folder;
    std::string label;
    std::string recoLegend;
    std::string unfoldedLegend;
    std::string effLegend;
    Color_t recoColor;
    Color_t unfoldColor;
    Style_t recoMarker;
    Style_t unfoldMarker;
};

Dataset OpenDataset(const std::string& label,
                    bool isSim,
                    const std::string& path,
                    const std::string& topDirName,
                    const std::string& centFolder,
                    const std::string& centSuffix,
                    const std::string& centLabel)
{
    Dataset ds;
    ds.label = label;
    ds.isSim = isSim;
    ds.trigger = topDirName;
    ds.topDirName = topDirName;
    ds.inFilePath = path;
    ds.centFolder = centFolder;
    ds.centSuffix = centSuffix;
    ds.centLabel = centLabel;
    ds.file = OpenRequiredFile(path, label);
    ds.topDir = dynamic_cast<TDirectory*>(ds.file->Get(topDirName.c_str()));
    if (!ds.topDir)
    {
        std::cerr << "[FATAL] Missing top-level directory " << topDirName
                  << " in " << path << "\n";
        std::exit(1);
    }
    return ds;
}

TH2D* TransposeTH2Local(const TH2* hIn,
                        const std::string& newName,
                        const std::string& newTitle)
{
    if (!hIn) return nullptr;

    const TAxis* ax = hIn->GetXaxis();
    const TAxis* ay = hIn->GetYaxis();
    const int nx = ax->GetNbins();
    const int ny = ay->GetNbins();

    const bool xVar = (ay->GetXbins() && ay->GetXbins()->GetSize() > 0);
    const bool yVar = (ax->GetXbins() && ax->GetXbins()->GetSize() > 0);

    TH2D* hOut = nullptr;
    if (xVar && yVar)
    {
        hOut = new TH2D(newName.c_str(), newTitle.c_str(),
                        ny, ay->GetXbins()->GetArray(),
                        nx, ax->GetXbins()->GetArray());
    }
    else if (xVar && !yVar)
    {
        hOut = new TH2D(newName.c_str(), newTitle.c_str(),
                        ny, ay->GetXbins()->GetArray(),
                        nx, ax->GetXmin(), ax->GetXmax());
    }
    else if (!xVar && yVar)
    {
        hOut = new TH2D(newName.c_str(), newTitle.c_str(),
                        ny, ay->GetXmin(), ay->GetXmax(),
                        nx, ax->GetXbins()->GetArray());
    }
    else
    {
        hOut = new TH2D(newName.c_str(), newTitle.c_str(),
                        ny, ay->GetXmin(), ay->GetXmax(),
                        nx, ax->GetXmin(), ax->GetXmax());
    }

    hOut->SetDirectory(nullptr);
    hOut->Sumw2();

    for (int ix = 0; ix <= nx + 1; ++ix)
    {
        for (int iy = 0; iy <= ny + 1; ++iy)
        {
            hOut->SetBinContent(iy, ix, hIn->GetBinContent(ix, iy));
            hOut->SetBinError(iy, ix, hIn->GetBinError(ix, iy));
        }
    }

    return hOut;
}

TH1* CloneHist(const TH1* h, const std::string& name)
{
    if (!h) return nullptr;
    TH1* out = dynamic_cast<TH1*>(h->Clone(name.c_str()));
    if (!out) return nullptr;
    out->SetDirectory(nullptr);
    EnsureSumw2(out);
    return out;
}

double HistMaxPositive(const std::vector<TH1*>& hs)
{
    double ymax = 0.0;
    for (const TH1* h : hs)
    {
        if (!h) continue;
        for (int i = 1; i <= h->GetNbinsX(); ++i)
        {
            const double y = h->GetBinContent(i);
            if (std::isfinite(y) && y > ymax) ymax = y;
        }
    }
    return ymax;
}

double HistMaxPositiveInRange(const std::vector<TH1*>& hs, double xmin, double xmax)
{
    double ymax = 0.0;
    for (const TH1* h : hs)
    {
        if (!h) continue;
        for (int i = 1; i <= h->GetNbinsX(); ++i)
        {
            const double lo = h->GetXaxis()->GetBinLowEdge(i);
            const double hi = h->GetXaxis()->GetBinUpEdge(i);
            if (lo < xmin || hi > xmax) continue;
            const double y = h->GetBinContent(i);
            if (std::isfinite(y) && y > ymax) ymax = y;
        }
    }
    return ymax;
}

double HistMinPositive(const std::vector<TH1*>& hs)
{
    double ymin = std::numeric_limits<double>::max();
    for (const TH1* h : hs)
    {
        if (!h) continue;
        for (int i = 1; i <= h->GetNbinsX(); ++i)
        {
            const double y = h->GetBinContent(i);
            if (std::isfinite(y) && y > 0.0 && y < ymin) ymin = y;
        }
    }
    if (ymin == std::numeric_limits<double>::max()) return 0.1;
    return ymin;
}

double HistMinPositiveInRange(const std::vector<TH1*>& hs, double xmin, double xmax)
{
    double ymin = std::numeric_limits<double>::max();
    for (const TH1* h : hs)
    {
        if (!h) continue;
        for (int i = 1; i <= h->GetNbinsX(); ++i)
        {
            const double lo = h->GetXaxis()->GetBinLowEdge(i);
            const double hi = h->GetXaxis()->GetBinUpEdge(i);
            if (lo < xmin || hi > xmax) continue;
            const double y = h->GetBinContent(i);
            if (std::isfinite(y) && y > 0.0 && y < ymin) ymin = y;
        }
    }
    if (ymin == std::numeric_limits<double>::max()) return 0.1;
    return ymin;
}

void StyleHist(TH1* h, Color_t color, Style_t marker, Style_t lineStyle = 1)
{
    if (!h) return;
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(1.15);
    h->SetLineWidth(2);
    h->SetLineStyle(lineStyle);
}

TH1* MakeRatioHist(const TH1* num, const TH1* den, const std::string& name)
{
    TH1* r = CloneHist(num, name);
    if (!r || !den) return r;
    r->Reset("ICES");
    for (int i = 1; i <= r->GetNbinsX(); ++i)
    {
        const double n = num->GetBinContent(i);
        const double d = den->GetBinContent(i);
        const double en = num->GetBinError(i);
        const double ed = den->GetBinError(i);
        if (d <= 0.0 || n <= 0.0) continue;
        const double val = n / d;
        const double rel2 = (en > 0.0 ? (en / n) * (en / n) : 0.0)
                          + (ed > 0.0 ? (ed / d) * (ed / d) : 0.0);
        r->SetBinContent(i, val);
        r->SetBinError(i, val * std::sqrt(rel2));
    }
    return r;
}

TH1* MakePhotonEfficiencyMC(const TH1* truth, const TH1* misses, const std::string& name)
{
    TH1* eff = CloneHist(truth, name);
    if (!eff) return nullptr;
    eff->Reset("ICES");
    for (int i = 1; i <= eff->GetNbinsX(); ++i)
    {
        const double t = truth ? truth->GetBinContent(i) : 0.0;
        const double m = misses ? misses->GetBinContent(i) : 0.0;
        if (t <= 0.0) continue;
        const double e = std::max(0.0, std::min(1.2, 1.0 - m / t));
        eff->SetBinContent(i, e);
        eff->SetBinError(i, std::sqrt(std::max(0.0, e * (1.0 - std::min(e, 1.0)) / t)));
    }
    return eff;
}

void ComputeABCDSignalCounts(Dataset& data,
                             Dataset& sim,
                             int iCanon,
                             bool leakageCorrected,
                             double& A,
                             double& B,
                             double& C,
                             double& D,
                             double& SA,
                             double& eSA)
{
    const PtBin& b = PtBins()[iCanon];
    A = Read1BinCount(data, "h_isIsolated_isTight" + b.suffix);
    B = Read1BinCount(data, "h_notIsolated_isTight" + b.suffix);
    C = Read1BinCount(data, "h_isIsolated_notTight" + b.suffix);
    D = Read1BinCount(data, "h_notIsolated_notTight" + b.suffix);

    SA = 0.0;
    if (A > 0.0 && D > 0.0)
    {
        SA = std::max(0.0, A - B * (C / D));
    }

    if (leakageCorrected)
    {
        double fB = 0.0, fC = 0.0, fD = 0.0;
        TH1* hSig = GetObj<TH1>(sim, "h_sigABCD_MC" + b.suffix, false, true, false);
        const double Asig = hSig ? hSig->GetBinContent(1) : 0.0;
        if (hSig && Asig > 0.0)
        {
            fB = hSig->GetBinContent(2) / Asig;
            fC = hSig->GetBinContent(3) / Asig;
            fD = hSig->GetBinContent(4) / Asig;
        }

        double SAcorr = 0.0;
        if (SolveLeakageCorrectedSA(A, B, C, D, fB, fC, fD, SAcorr))
        {
            SA = std::min(std::max(SAcorr, 0.0), A);
        }
    }

    double varSA = 0.0;
    if (A > 0.0) varSA += A;
    if (D > 0.0)
    {
        const double dSdB = -(C / D);
        const double dSdC = -(B / D);
        const double dSdD = (B * C) / (D * D);
        if (B > 0.0) varSA += dSdB * dSdB * B;
        if (C > 0.0) varSA += dSdC * dSdC * C;
        if (D > 0.0) varSA += dSdD * dSdD * D;
    }
    eSA = (varSA > 0.0) ? std::sqrt(varSA) : 0.0;
}

TH1* MakeLeakageCorrectedReco(Dataset& data,
                              Dataset& sim,
                              const TH1* rawReco,
                              const std::string& name)
{
    TH1* h = CloneHist(rawReco, name);
    if (!h) return nullptr;
    h->Reset("ICES");

    const auto& recoBins = UnfoldRecoPtBins();
    for (const PtBin& ub : recoBins)
    {
        int iCanon = -1;
        if (ub.lo == 8 && ub.hi == 10)
        {
            iCanon = 0;
        }
        else if (ub.lo == 35 && ub.hi == 40)
        {
            iCanon = kNPtBins - 1;
        }
        else
        {
            for (int j = 0; j < kNPtBins; ++j)
            {
                const PtBin& pb = PtBins()[j];
                if (pb.lo == ub.lo && pb.hi == ub.hi)
                {
                    iCanon = j;
                    break;
                }
            }
        }
        if (iCanon < 0) continue;

        double A = 0.0, B = 0.0, C = 0.0, D = 0.0, SA = 0.0, eSA = 0.0;
        ComputeABCDSignalCounts(data, sim, iCanon, true, A, B, C, D, SA, eSA);
        const double cen = 0.5 * (ub.lo + ub.hi);
        const int ib = h->GetXaxis()->FindBin(cen);
        if (ib < 1 || ib > h->GetNbinsX()) continue;
        h->SetBinContent(ib, SA);
        h->SetBinError(ib, eSA);
    }

    return h;
}

TH1* UnfoldPhotonSpectrum(const TH1* reco,
                          TH1* recoSim,
                          TH1* truthSim,
                          TH2* responseTruthXReco,
                          const std::string& name,
                          int nIter = 10)
{
    TH2D* responseRecoXTruth = TransposeTH2Local(responseTruthXReco,
                                                 name + "_response_recoXtruth",
                                                 "");
    RooUnfoldResponse response(recoSim, truthSim, responseRecoXTruth,
                               (name + "_resp").c_str(),
                               (name + "_resp").c_str());
    RooUnfoldBayes unfold(&response, const_cast<TH1*>(reco), nIter);
    TH1* out = dynamic_cast<TH1*>(unfold.Hreco(RooUnfold::kCovariance));
    if (!out)
    {
        delete responseRecoXTruth;
        return nullptr;
    }
    TH1* clone = CloneHist(out, name);
    delete responseRecoXTruth;
    return clone;
}

void DrawYieldOverlay(const FocusSample& s,
                      TH1* rawReco,
                      TH1* rawUnfolded,
                      TH1* corrReco,
                      TH1* corrUnfolded,
                      const std::string& outPath)
{
    StyleHist(rawReco, kRed + 1, 20);
    StyleHist(rawUnfolded, kBlack, 24);
    StyleHist(corrReco, kOrange + 7, 21);
    StyleHist(corrUnfolded, kBlue + 1, 25);

    TH1* rRaw = MakeRatioHist(rawUnfolded, rawReco, s.key + "_ratio_raw");
    TH1* rCorr = MakeRatioHist(corrUnfolded, corrReco, s.key + "_ratio_corr");
    StyleHist(rRaw, kBlack, 24);
    StyleHist(rCorr, kBlue + 1, 25);

    TCanvas c(("c_yield_" + s.key).c_str(), "", 980, 760);
    TPad pTop("pTop", "", 0.0, 0.30, 1.0, 1.0);
    TPad pBot("pBot", "", 0.0, 0.0, 1.0, 0.30);
    pTop.SetBottomMargin(0.025);
    pTop.SetTopMargin(0.28);
    pTop.SetLeftMargin(0.12);
    pTop.SetRightMargin(0.035);
    pBot.SetTopMargin(0.025);
    pBot.SetBottomMargin(0.31);
    pBot.SetLeftMargin(0.12);
    pBot.SetRightMargin(0.035);
    pTop.Draw();
    pBot.Draw();

    pTop.cd();
    pTop.SetLogy();
    TH1* frame = CloneHist(rawReco, s.key + "_yield_frame");
    frame->Reset("ICES");
    const double ymin = std::max(0.05, HistMinPositive({rawReco, rawUnfolded, corrReco, corrUnfolded}) * 0.45);
    const double ymax = std::max(10.0, HistMaxPositive({rawReco, rawUnfolded, corrReco, corrUnfolded}) * 35.0);
    frame->SetTitle("");
    frame->GetXaxis()->SetRangeUser(8.0, 40.0);
    frame->GetYaxis()->SetRangeUser(ymin, ymax);
    frame->GetYaxis()->SetTitle("Counts");
    frame->GetYaxis()->SetTitleSize(0.055);
    frame->GetYaxis()->SetLabelSize(0.045);
    frame->GetXaxis()->SetLabelSize(0.0);
    frame->Draw("axis");

    rawReco->Draw("E1 same");
    rawUnfolded->Draw("E1 same");
    corrReco->Draw("E1 same");
    corrUnfolded->Draw("E1 same");

    TLegend leg(0.48, 0.50, 0.90, 0.76);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.035);
    leg.AddEntry(rawReco, "DATA reco, no purity corr.", "lep");
    leg.AddEntry(rawUnfolded, "Unfolded, no purity corr.", "lep");
    leg.AddEntry(corrReco, "DATA reco, leakage-corrected", "lep");
    leg.AddEntry(corrUnfolded, "Unfolded, leakage-corrected", "lep");
    leg.Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.029);
    lat.DrawLatex(0.145, 0.955, "Photon yield overlay, Bayes it=10");
    lat.DrawLatex(0.145, 0.910, s.triggerLabel.c_str());
    lat.SetTextSize(0.025);
    lat.DrawLatex(0.145, 0.865, s.cutLine1.c_str());
    lat.DrawLatex(0.145, 0.825, s.cutLine2.c_str());
    if (!s.centLabel.empty()) lat.DrawLatex(0.145, 0.785, s.centLabel.c_str());

    pBot.cd();
    TH1* rFrame = CloneHist(rRaw, s.key + "_ratio_frame");
    rFrame->Reset("ICES");
    const double rMax = std::max(1.5, HistMaxPositive({rRaw, rCorr}) * 1.35);
    rFrame->SetTitle("");
    rFrame->GetXaxis()->SetRangeUser(8.0, 40.0);
    rFrame->GetYaxis()->SetRangeUser(0.0, rMax);
    rFrame->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
    rFrame->GetYaxis()->SetTitle("Unfolded / reco");
    rFrame->GetXaxis()->SetTitleSize(0.105);
    rFrame->GetXaxis()->SetLabelSize(0.085);
    rFrame->GetYaxis()->SetTitleSize(0.080);
    rFrame->GetYaxis()->SetLabelSize(0.075);
    rFrame->GetYaxis()->SetTitleOffset(0.55);
    rFrame->GetYaxis()->SetNdivisions(505);
    rFrame->Draw("axis");
    TLine one(8.0, 1.0, 40.0, 1.0);
    one.SetLineStyle(2);
    one.SetLineColor(kGray + 2);
    one.Draw("same");
    rRaw->Draw("E1 same");
    rCorr->Draw("E1 same");

    SaveCanvas(c, outPath);
}

void DrawYieldPair(const FocusSample& s,
                   const CorrectionView& view,
                   TH1* reco,
                   TH1* unfolded,
                   const std::string& outPath)
{
    constexpr double xPlotMin = 10.0;
    constexpr double xPlotMax = 35.0;

    StyleHist(reco, view.recoColor, view.recoMarker);
    StyleHist(unfolded, view.unfoldColor, view.unfoldMarker);

    TH1* ratio = MakeRatioHist(unfolded, reco, s.key + "_" + view.folder + "_ratio");
    StyleHist(ratio, view.unfoldColor, view.unfoldMarker);

    TCanvas c(("c_yield_pair_" + s.key + "_" + view.folder).c_str(), "", 980, 760);
    TPad pTop("pTop", "", 0.0, 0.30, 1.0, 1.0);
    TPad pBot("pBot", "", 0.0, 0.0, 1.0, 0.30);
    pTop.SetBottomMargin(0.025);
    pTop.SetTopMargin(0.28);
    pTop.SetLeftMargin(0.12);
    pTop.SetRightMargin(0.035);
    pBot.SetTopMargin(0.025);
    pBot.SetBottomMargin(0.31);
    pBot.SetLeftMargin(0.12);
    pBot.SetRightMargin(0.035);
    pTop.Draw();
    pBot.Draw();

    pTop.cd();
    pTop.SetLogy();
    TH1* frame = CloneHist(reco, s.key + "_" + view.folder + "_yield_frame");
    frame->Reset("ICES");
    const double ymin = std::max(0.05, HistMinPositiveInRange({reco, unfolded}, xPlotMin, xPlotMax) * 0.45);
    const double ymax = std::max(10.0, HistMaxPositiveInRange({reco, unfolded}, xPlotMin, xPlotMax) * 8.0);
    frame->SetTitle("");
    frame->GetXaxis()->SetRangeUser(xPlotMin, xPlotMax);
    frame->GetYaxis()->SetRangeUser(ymin, ymax);
    frame->GetYaxis()->SetTitle("Counts");
    frame->GetYaxis()->SetTitleSize(0.055);
    frame->GetYaxis()->SetLabelSize(0.045);
    frame->GetXaxis()->SetLabelSize(0.0);
    frame->Draw("axis");

    reco->Draw("E1 same");
    unfolded->Draw("E1 same");

    TLegend leg(0.52, 0.56, 0.90, 0.70);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.037);
    leg.AddEntry(reco, view.recoLegend.c_str(), "lep");
    leg.AddEntry(unfolded, view.unfoldedLegend.c_str(), "lep");
    leg.Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.029);
    lat.DrawLatex(0.145, 0.955, ("Photon yield overlay, " + view.label + ", Bayes it=10").c_str());
    lat.DrawLatex(0.145, 0.910, s.triggerLabel.c_str());
    lat.SetTextSize(0.025);
    lat.DrawLatex(0.145, 0.865, s.cutLine1.c_str());
    lat.DrawLatex(0.145, 0.825, s.cutLine2.c_str());
    if (!s.centLabel.empty()) lat.DrawLatex(0.145, 0.785, s.centLabel.c_str());

    pBot.cd();
    TH1* rFrame = CloneHist(ratio, s.key + "_" + view.folder + "_ratio_frame");
    rFrame->Reset("ICES");
    const double rMax = std::max(1.5, HistMaxPositiveInRange({ratio}, xPlotMin, xPlotMax) * 1.25);
    rFrame->SetTitle("");
    rFrame->GetXaxis()->SetRangeUser(xPlotMin, xPlotMax);
    rFrame->GetYaxis()->SetRangeUser(0.0, rMax);
    rFrame->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
    rFrame->GetYaxis()->SetTitle("Unfolded / reco");
    rFrame->GetXaxis()->SetTitleSize(0.105);
    rFrame->GetXaxis()->SetLabelSize(0.085);
    rFrame->GetYaxis()->SetTitleSize(0.080);
    rFrame->GetYaxis()->SetLabelSize(0.075);
    rFrame->GetYaxis()->SetTitleOffset(0.55);
    rFrame->GetYaxis()->SetNdivisions(505);
    rFrame->Draw("axis");
    TLine one(xPlotMin, 1.0, xPlotMax, 1.0);
    one.SetLineStyle(2);
    one.SetLineColor(kGray + 2);
    one.Draw("same");
    ratio->Draw("E1 same");

    SaveCanvas(c, outPath);
}

void DrawResponseMatrix(const FocusSample& s,
                        TH2* responseTruthXReco,
                        const std::string& outPath)
{
    TH2D* h = TransposeTH2Local(responseTruthXReco,
                                s.key + "_response_reco_vs_truth_draw",
                                "");
    h->SetTitle("");
    h->GetXaxis()->SetTitle("p_{T}^{#gamma, reco} [GeV]");
    h->GetYaxis()->SetTitle("p_{T}^{#gamma, truth} [GeV]");
    h->GetZaxis()->SetTitle("Counts");
    h->GetXaxis()->SetRangeUser(8.0, 40.0);
    h->GetYaxis()->SetRangeUser(5.0, 40.0);
    h->SetMinimum(1.0);

    TCanvas c(("c_response_" + s.key).c_str(), "", 900, 760);
    c.SetLeftMargin(0.13);
    c.SetRightMargin(0.16);
    c.SetTopMargin(0.28);
    c.SetBottomMargin(0.12);
    c.SetLogz();
    h->Draw("COLZ");

    TLine lRecoLo(10.0, 5.0, 10.0, 40.0);
    TLine lRecoHi(35.0, 5.0, 35.0, 40.0);
    TLine lTruthLo(8.0, 10.0, 40.0, 10.0);
    TLine lTruthHi(8.0, 35.0, 40.0, 35.0);
    for (TLine* l : {&lRecoLo, &lRecoHi, &lTruthLo, &lTruthHi})
    {
        l->SetLineColor(kBlack);
        l->SetLineStyle(2);
        l->SetLineWidth(2);
        l->Draw("same");
    }

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.032);
    lat.DrawLatex(0.15, 0.955, "Photon response matrix");
    lat.SetTextSize(0.027);
    lat.DrawLatex(0.15, 0.910, s.simLabel.c_str());
    lat.DrawLatex(0.15, 0.865, s.cutLine1.c_str());
    lat.DrawLatex(0.15, 0.825, s.cutLine2.c_str());
    if (!s.centLabel.empty()) lat.DrawLatex(0.15, 0.785, s.centLabel.c_str());

    SaveCanvas(c, outPath);
    delete h;
}

void DrawEfficiencyDiagnostic(const FocusSample& s,
                              TH1* rawReco,
                              TH1* rawUnfolded,
                              TH1* corrReco,
                              TH1* corrUnfolded,
                              TH1* truthSim,
                              TH1* missesSim,
                              const std::string& outPath)
{
    TH1* effRaw = MakeRatioHist(rawReco, rawUnfolded, s.key + "_eff_raw");
    TH1* effCorr = MakeRatioHist(corrReco, corrUnfolded, s.key + "_eff_corr");
    TH1* effMC = MakePhotonEfficiencyMC(truthSim, missesSim, s.key + "_eff_mc");

    StyleHist(effRaw, kBlack, 20);
    StyleHist(effCorr, kRed + 1, 21);
    StyleHist(effMC, kBlue + 1, 22);

    TCanvas c(("c_eff_" + s.key).c_str(), "", 900, 720);
    c.SetLeftMargin(0.12);
    c.SetRightMargin(0.035);
    c.SetTopMargin(0.23);
    c.SetBottomMargin(0.13);

    TH1* frame = CloneHist(effRaw, s.key + "_eff_frame");
    frame->Reset("ICES");
    frame->SetTitle("");
    frame->GetXaxis()->SetRangeUser(5.0, 40.0);
    frame->GetYaxis()->SetRangeUser(0.0, 1.25);
    frame->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
    frame->GetYaxis()->SetTitle("#epsilon_{#gamma}(p_{T}^{#gamma})");
    frame->GetXaxis()->SetTitleSize(0.052);
    frame->GetYaxis()->SetTitleSize(0.052);
    frame->GetXaxis()->SetLabelSize(0.043);
    frame->GetYaxis()->SetLabelSize(0.043);
    frame->Draw("axis");

    effRaw->Draw("E1 same");
    effCorr->Draw("E1 same");
    effMC->Draw("E1 same");

    TLegend leg(0.40, 0.805, 0.89, 0.985);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.029);
    leg.AddEntry(effRaw, "DATA reco / unfolded, no purity corr.", "lep");
    leg.AddEntry(effCorr, "DATA reco / unfolded, leakage-corrected", "lep");
    leg.AddEntry(effMC, "SIM 1 - misses / truth", "lep");
    leg.Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.030);
    lat.DrawLatex(0.14, 0.975, "Photon efficiency diagnostic");
    lat.SetTextSize(0.025);
    lat.DrawLatex(0.14, 0.930, s.triggerLabel.c_str());
    if (!s.centLabel.empty()) lat.DrawLatex(0.14, 0.890, s.centLabel.c_str());

    SaveCanvas(c, outPath);
}

void DrawEfficiencyDiagnosticPair(const FocusSample& s,
                                  const CorrectionView& view,
                                  TH1* reco,
                                  TH1* unfolded,
                                  TH1* truthSim,
                                  TH1* missesSim,
                                  const std::string& outPath)
{
    TH1* effData = MakeRatioHist(reco, unfolded, s.key + "_" + view.folder + "_eff_data");
    TH1* effMC = MakePhotonEfficiencyMC(truthSim, missesSim, s.key + "_" + view.folder + "_eff_mc");

    StyleHist(effData, view.recoColor, view.recoMarker);
    StyleHist(effMC, kBlue + 1, 22);

    TCanvas c(("c_eff_pair_" + s.key + "_" + view.folder).c_str(), "", 900, 720);
    c.SetLeftMargin(0.12);
    c.SetRightMargin(0.035);
    c.SetTopMargin(0.23);
    c.SetBottomMargin(0.13);

    TH1* frame = CloneHist(effData, s.key + "_" + view.folder + "_eff_frame");
    frame->Reset("ICES");
    frame->SetTitle("");
    frame->GetXaxis()->SetRangeUser(5.0, 40.0);
    frame->GetYaxis()->SetRangeUser(0.0, 1.25);
    frame->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
    frame->GetYaxis()->SetTitle("#epsilon_{#gamma}(p_{T}^{#gamma})");
    frame->GetXaxis()->SetTitleSize(0.052);
    frame->GetYaxis()->SetTitleSize(0.052);
    frame->GetXaxis()->SetLabelSize(0.043);
    frame->GetYaxis()->SetLabelSize(0.043);
    frame->Draw("axis");

    effData->Draw("E1 same");
    effMC->Draw("E1 same");

    TLegend leg(0.40, 0.785, 0.89, 0.910);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.031);
    leg.AddEntry(effData, view.effLegend.c_str(), "lep");
    leg.AddEntry(effMC, "SIM 1 - misses / truth", "lep");
    leg.Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.030);
    lat.DrawLatex(0.14, 0.975, ("Photon efficiency diagnostic, " + view.label).c_str());
    lat.SetTextSize(0.025);
    lat.DrawLatex(0.14, 0.930, s.triggerLabel.c_str());
    if (!s.centLabel.empty()) lat.DrawLatex(0.14, 0.890, s.centLabel.c_str());

    SaveCanvas(c, outPath);
}

void ProcessSample(const FocusSample& s, std::vector<std::string>& outputs)
{
    EnsureDir(s.outDir);

    Dataset data = OpenDataset("DATA", false, s.dataPath, s.dataTopDir,
                               s.centFolder, s.centSuffix, s.centLabel);
    Dataset sim = OpenDataset("SIM", true, s.simPath, s.simTopDir,
                              s.centFolder, s.centSuffix, s.centLabel);

    TH1* rawReco = CloneHist(GetObj<TH1>(data, "h_unfoldRecoPho_pTgamma", true, true, false),
                             s.key + "_raw_reco");
    TH1* recoSim = CloneHist(GetObj<TH1>(sim, "h_unfoldRecoPho_pTgamma", true, true, false),
                             s.key + "_reco_sim");
    TH1* truthSim = CloneHist(GetObj<TH1>(sim, "h_unfoldTruthPho_pTgamma", true, true, false),
                              s.key + "_truth_sim");
    TH1* missesSim = CloneHist(GetObj<TH1>(sim, "h_unfoldTruthPhoMisses_pTgamma", true, true, false),
                               s.key + "_misses_sim");
    TH2* responseTruthXReco = dynamic_cast<TH2*>(
        GetObj<TH2>(sim, "h2_unfoldResponsePho_pTgamma", true, true, false));

    if (!rawReco || !recoSim || !truthSim || !missesSim || !responseTruthXReco)
    {
        std::cerr << "[FATAL] Missing required photon unfolding histogram(s) for "
                  << s.key << "\n";
        std::exit(1);
    }

    TH1* corrReco = MakeLeakageCorrectedReco(data, sim, rawReco, s.key + "_leakage_corr_reco");
    TH1* rawUnfolded = UnfoldPhotonSpectrum(rawReco, recoSim, truthSim, responseTruthXReco,
                                            s.key + "_unfolded_raw", 10);
    TH1* corrUnfolded = UnfoldPhotonSpectrum(corrReco, recoSim, truthSim, responseTruthXReco,
                                             s.key + "_unfolded_leakage_corr", 10);

    if (!corrReco || !rawUnfolded || !corrUnfolded)
    {
        std::cerr << "[FATAL] Failed to produce unfolded spectra for "
                  << s.key << "\n";
        std::exit(1);
    }

    const CorrectionView nonCorrected{
        "nonCorrected",
        "no purity correction",
        "DATA reco, no purity corr.",
        "Unfolded, no purity corr.",
        "DATA reco / unfolded, no purity corr.",
        kRed + 1,
        kBlack,
        20,
        24
    };

    const CorrectionView leakageCorrected{
        "leakageCorrected",
        "leakage-corrected",
        "DATA reco, leakage-corrected",
        "Unfolded, leakage-corrected",
        "DATA reco / unfolded, leakage-corrected",
        kOrange + 7,
        kBlue + 1,
        21,
        25
    };

    for (const CorrectionView& view : {nonCorrected, leakageCorrected})
    {
        const std::string modeDir = JoinPath(s.outDir, view.folder);
        EnsureDir(modeDir);
        const bool isLeakage = (view.folder == "leakageCorrected");
        TH1* reco = isLeakage ? corrReco : rawReco;
        TH1* unfolded = isLeakage ? corrUnfolded : rawUnfolded;

        const std::string yieldPath = JoinPath(modeDir, "01_photon_yield_overlay.png");
        const std::string responsePath = JoinPath(modeDir, "02_photon_response_matrix_reco_vs_truth.png");
        const std::string effPath = JoinPath(modeDir, "03_photon_efficiency_diagnostic.png");

        DrawYieldPair(s, view, reco, unfolded, yieldPath);
        DrawResponseMatrix(s, responseTruthXReco, responsePath);
        DrawEfficiencyDiagnosticPair(s, view, reco, unfolded, truthSim, missesSim, effPath);

        outputs.push_back(yieldPath);
        outputs.push_back(responsePath);
        outputs.push_back(effPath);
    }
}
} // namespace

void MakeFocusedPhotonUnfoldingQuick()
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleFillColor(0);

    const std::string baseOut =
        "dataOutput/focusedPhotonUnfoldingQuick/pp_vs_auau020_referenceCuts";

    const FocusSample pp{
        "pp",
        "pp",
        "pp data: Photon 4 GeV + MBD NS #geq 1",
        "jetMinPt5, #Delta#phi>7#pi/8, |v_{z}|<60 cm, isoR40 fixedIso2GeV, reference ID",
        "jet p_{T}>5 GeV, #Delta#phi>7#pi/8, |v_{z}|<60 cm",
        "isoR40 fixedIso2GeV, reference preselection/tight/nonTight",
        "SIM: photon+jet 5+10+20 merged",
        "InputFiles/pp24/RecoilJets_pp_ALL_jetMinPt5_7pi_8_vz60_isoR40_fixedIso2GeV_preselectionReference_tightReference_nonTightReference.root",
        "dataOutput/combinedSimOnly/jetMinPt5_7pi_8_vz60_isoR40_fixedIso2GeV_preselectionReference_tightReference_nonTightReference/photonJet5and10and20merged_SIM/RecoilJets_photonjet5plus10plus20_MERGED.root",
        "Photon_4_GeV_plus_MBD_NS_geq_1",
        "SIM",
        "",
        "",
        "",
        JoinPath(baseOut, "pp")
    };

    const FocusSample auau020{
        "auau_0_20",
        "AuAu 0-20%",
        "Au+Au data: Photon 10 GeV + MBD NS #geq 2, |v_{z}|<150 cm",
        "jetMinPt5, #Delta#phi>7#pi/8, |v_{z}|<60 cm, isoR40 sliding, baseVariant, reference ID",
        "jet p_{T}>5 GeV, #Delta#phi>7#pi/8, |v_{z}|<60 cm",
        "isoR40 sliding, baseVariant, reference preselection/tight/nonTight",
        "Embedded SIM: photon+jet 12+20 merged",
        "InputFiles/auau25/RecoilJets_auau_ALL_jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionReference_tightReference_nonTightReference.root",
        "InputFiles/simEmbedded/merged/jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionReference_tightReference_nonTightReference/photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root",
        "photon_10_plus_MBD_NS_geq_2_vtx_lt_150",
        "SIM",
        "0_20",
        "_cent_0_20",
        "Centrality: 0-20%",
        JoinPath(baseOut, "auau_0_20")
    };

    std::vector<std::string> outputs;
    ProcessSample(pp, outputs);

    const bool mergeOK = BuildMergedSIMFile_PhotonSlices(
        {
            "InputFiles/simEmbedded/RecoilJets_embeddedPhoton12_ALL_jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionReference_tightReference_nonTightReference.root",
            "InputFiles/simEmbedded/RecoilJets_embeddedPhoton20_ALL_jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionReference_tightReference_nonTightReference.root"
        },
        {kSigmaEmbeddedPhoton12_pb, kSigmaEmbeddedPhoton20_pb},
        auau020.simPath,
        kDirSIM,
        {"embeddedPhoton12", "embeddedPhoton20"});
    if (!mergeOK)
    {
        std::cerr << "[FATAL] Failed to rebuild AuAu embedded photon+jet 12+20 merge.\n";
        std::exit(1);
    }

    ProcessSample(auau020, outputs);

    const std::string indexPath = JoinPath(baseOut, "focused_photon_unfolding_index.txt");
    EnsureParentDirForFile(indexPath);
    std::ofstream index(indexPath);
    index << "Focused photon unfolding quick output\n";
    index << "PP data: " << pp.dataPath << "\n";
    index << "PP sim:  " << pp.simPath << "\n";
    index << "AuAu data: " << auau020.dataPath << "\n";
    index << "AuAu sim:  " << auau020.simPath << "\n\n";
    for (const std::string& out : outputs)
    {
        index << out << "\n";
        std::cout << "[OK] " << out << "\n";
    }
}
