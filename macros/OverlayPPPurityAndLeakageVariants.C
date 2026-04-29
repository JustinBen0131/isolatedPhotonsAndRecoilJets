#include "AnalyzeRecoilJets.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace
{
struct VariantSpec
{
    std::string label;
    std::string tag;
    int color;
    int marker;
};

double ReadBin1(TDirectory* dir, const std::string& name)
{
    if (!dir) return 0.0;
    TH1* h = dynamic_cast<TH1*>(dir->Get(name.c_str()));
    return h ? h->GetBinContent(1) : 0.0;
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

std::string PPPath(const std::string& tag)
{
    return ARJ::kInputBase + "/pp24/RecoilJets_pp_ALL_" + tag + ".root";
}

std::string SimPath(const std::string& sample, const std::string& tag)
{
    return ARJ::kInputBase + "/simPhotonJet/RecoilJets_" + sample + "_ALL_" + tag + ".root";
}

std::string MergedSimPath(const std::string& tag)
{
    return ARJ::kOutputBase + "/combinedSimOnly/" + tag +
           "/photonJet5and10and20merged_SIM/RecoilJets_photonjet5plus10plus20_MERGED.root";
}

bool EnsureMergedSim(const VariantSpec& spec)
{
    const std::string out = MergedSimPath(spec.tag);
    if (!gSystem->AccessPathName(out.c_str()))
    {
        std::cout << "[OK] merged SIM exists for " << spec.label << ": " << out << "\n";
        return true;
    }

    const std::vector<std::string> inputs = {
        SimPath("photonjet5", spec.tag),
        SimPath("photonjet10", spec.tag),
        SimPath("photonjet20", spec.tag)
    };

    for (const auto& path : inputs)
    {
        if (gSystem->AccessPathName(path.c_str()))
        {
            std::cerr << "[ERROR] Missing SIM input for " << spec.label << ": " << path << "\n";
            return false;
        }
    }

    std::cout << "[MERGE] Building 5+10+20 merged SIM for " << spec.label << "\n";
    return ARJ::BuildMergedSIMFile_PhotonSlices(
        inputs,
        {ARJ::kSigmaPhoton5_pb, ARJ::kSigmaPhoton10_pb, ARJ::kSigmaPhoton20_pb},
        out,
        ARJ::kDirSIM,
        {"photonJet5", "photonJet10", "photonJet20"});
}

bool BuildPurityGraph(const VariantSpec& spec,
                      std::vector<double>& x,
                      std::vector<double>& ex,
                      std::vector<double>& y,
                      std::vector<double>& ey)
{
    TFile f(PPPath(spec.tag).c_str(), "READ");
    if (f.IsZombie())
    {
        std::cerr << "[ERROR] Cannot open pp input for " << spec.label << ": " << PPPath(spec.tag) << "\n";
        return false;
    }

    TDirectory* top = f.GetDirectory(ARJ::kTriggerPP.c_str());
    if (!top) top = &f;

    const int n = ARJ::kNPtBins - 1;
    x.assign(n, 0.0);
    ex.assign(n, 0.0);
    y.assign(n, 0.0);
    ey.assign(n, 0.0);

    for (int i = 1; i < ARJ::kNPtBins; ++i)
    {
        const int j = i - 1;
        const auto& b = ARJ::PtBins()[i];
        const std::string& suf = b.suffix;

        const double A = ReadBin1(top, "h_isIsolated_isTight" + suf);
        const double B = ReadBin1(top, "h_notIsolated_isTight" + suf);
        const double C = ReadBin1(top, "h_isIsolated_notTight" + suf);
        const double D = ReadBin1(top, "h_notIsolated_notTight" + suf);

        const double ptLo = ARJ::kPtEdges[(std::size_t)i];
        const double ptHi = ARJ::kPtEdges[(std::size_t)i + 1];
        x[j] = 0.5 * (ptLo + ptHi);
        ex[j] = 0.5 * (ptHi - ptLo);

        if (A > 0.0 && D > 0.0)
        {
            double sig = A - B * (C / D);
            if (sig < 0.0) sig = 0.0;
            y[j] = sig / A;
            ey[j] = RawPurityError(A, B, C, D);
        }
    }

    return true;
}

bool BuildLeakageCGraph(const VariantSpec& spec,
                        std::vector<double>& x,
                        std::vector<double>& ex,
                        std::vector<double>& y,
                        std::vector<double>& ey)
{
    if (!EnsureMergedSim(spec)) return false;

    TFile f(MergedSimPath(spec.tag).c_str(), "READ");
    if (f.IsZombie())
    {
        std::cerr << "[ERROR] Cannot open merged SIM for " << spec.label << ": " << MergedSimPath(spec.tag) << "\n";
        return false;
    }

    TDirectory* top = f.GetDirectory(ARJ::kDirSIM.c_str());
    if (!top) top = &f;

    const int n = ARJ::kNPtBins - 1;
    x.assign(n, 0.0);
    ex.assign(n, 0.0);
    y.assign(n, 0.0);
    ey.assign(n, 0.0);

    for (int i = 1; i < ARJ::kNPtBins; ++i)
    {
        const int j = i - 1;
        const auto& b = ARJ::PtBins()[i];
        const std::string hname = "h_sigABCD_MC" + b.suffix;
        TH1* h = dynamic_cast<TH1*>(top->Get(hname.c_str()));

        const double ptLo = ARJ::kPtEdges[(std::size_t)i];
        const double ptHi = ARJ::kPtEdges[(std::size_t)i + 1];
        x[j] = 0.5 * (ptLo + ptHi);
        ex[j] = 0.5 * (ptHi - ptLo);

        if (!h)
        {
            std::cerr << "[WARN] Missing " << hname << " in merged SIM for " << spec.label << "\n";
            continue;
        }

        const double As = h->GetBinContent(1);
        const double Cs = h->GetBinContent(3);
        y[j] = (As > 0.0) ? (Cs / As) : 0.0;
    }

    return true;
}

void DrawOverlay(const std::vector<VariantSpec>& specs,
                 const std::vector<std::vector<double>>& xs,
                 const std::vector<std::vector<double>>& exs,
                 const std::vector<std::vector<double>>& ys,
                 const std::vector<std::vector<double>>& eys,
                 const std::string& yTitle,
                 const std::string& title,
                 const std::string& outPath,
                 double yMax)
{
    TCanvas c("c_overlay", "c_overlay", 1000, 760);
    c.SetLeftMargin(0.13);
    c.SetRightMargin(0.04);
    c.SetBottomMargin(0.13);
    c.SetTopMargin(0.08);
    c.SetTicks(1, 1);

    TH1F frame("frame", "", 100, ARJ::kPtEdges[1], ARJ::kPtEdges.back());
    frame.SetDirectory(nullptr);
    frame.SetStats(0);
    frame.SetMinimum(0.0);
    frame.SetMaximum(yMax);
    frame.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
    frame.GetYaxis()->SetTitle(yTitle.c_str());
    frame.GetXaxis()->SetTitleSize(0.045);
    frame.GetYaxis()->SetTitleSize(0.045);
    frame.GetXaxis()->SetLabelSize(0.040);
    frame.GetYaxis()->SetLabelSize(0.040);
    frame.Draw();

    std::vector<TGraphErrors*> graphs;
    for (std::size_t i = 0; i < specs.size(); ++i)
    {
        if (xs[i].empty()) continue;
        auto* g = new TGraphErrors((int)xs[i].size(),
                                   const_cast<double*>(xs[i].data()),
                                   const_cast<double*>(ys[i].data()),
                                   const_cast<double*>(exs[i].data()),
                                   const_cast<double*>(eys[i].data()));
        g->SetLineColor(specs[i].color);
        g->SetMarkerColor(specs[i].color);
        g->SetLineWidth(2);
        g->SetMarkerStyle(specs[i].marker);
        g->SetMarkerSize(1.15);
        g->Draw("PE1 SAME");
        graphs.push_back(g);
    }

    TLegend leg(0.58, 0.70, 0.92, 0.89);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextFont(42);
    leg.SetTextSize(0.034);
    for (std::size_t i = 0; i < graphs.size() && i < specs.size(); ++i)
        leg.AddEntry(graphs[i], specs[i].label.c_str(), "pe");
    leg.Draw();

    TLatex text;
    text.SetNDC(true);
    text.SetTextFont(42);
    text.SetTextAlign(13);
    text.SetTextSize(0.038);
    text.DrawLatex(0.16, 0.93, title.c_str());
    text.SetTextSize(0.030);
    text.DrawLatex(0.16, 0.86, "Run24 pp, Photon 4 GeV + MBD NS #geq 1");
    text.DrawLatex(0.16, 0.815, "jetMinPt5, 7#pi/8, |v_{z}| < 60 cm, isoR40, sliding iso");

    c.SaveAs(outPath.c_str());
    std::cout << "[WROTE] " << outPath << "\n";

    for (auto* g : graphs) delete g;
}
}

void OverlayPPPurityAndLeakageVariants()
{
    const std::vector<VariantSpec> specs = {
        {"box-cuts",
         "jetMinPt5_7pi_8_vz60_isoR40_isSliding_preselectionReference_tightReference_nonTightReference",
         kBlack, 20},
        {"BDT with NPB pre",
         "jetMinPt5_7pi_8_vz60_isoR40_isSliding_preselectionVariantA_tightVariantA_nonTightReference",
         kBlue + 1, 21},
        {"BDT only",
         "jetMinPt5_7pi_8_vz60_isoR40_isSliding_preselectionVariantB_tightVariantA_nonTightReference",
         kRed + 1, 22}
    };

    gSystem->mkdir((ARJ::kOutputBase + "/pp").c_str(), true);

    std::vector<std::vector<double>> xs(specs.size()), exs(specs.size());
    std::vector<std::vector<double>> ys(specs.size()), eys(specs.size());
    for (std::size_t i = 0; i < specs.size(); ++i)
    {
        if (!BuildPurityGraph(specs[i], xs[i], exs[i], ys[i], eys[i]))
            std::cerr << "[WARN] Skipping purity graph for " << specs[i].label << "\n";
    }
    DrawOverlay(specs, xs, exs, ys, eys,
                "Purity (raw ABCD)",
                "Raw ABCD purity overlay",
                ARJ::kOutputBase + "/pp/rawPurity_overlay_isSliding_photonIDVariants.png",
                1.05);

    for (auto& v : xs) v.clear();
    for (auto& v : exs) v.clear();
    for (auto& v : ys) v.clear();
    for (auto& v : eys) v.clear();

    for (std::size_t i = 0; i < specs.size(); ++i)
    {
        if (!BuildLeakageCGraph(specs[i], xs[i], exs[i], ys[i], eys[i]))
            std::cerr << "[WARN] Skipping region-C leakage graph for " << specs[i].label << "\n";
    }
    DrawOverlay(specs, xs, exs, ys, eys,
                "Region C signal leakage, f_{C} = C_{sig}/A_{sig}",
                "Truth-signal leakage into ABCD region C",
                ARJ::kOutputBase + "/pp/signalLeakage_regionC_overlay_isSliding_photonIDVariants.png",
                1.05);
}
