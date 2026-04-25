#include "AnalyzeRecoilJets.h"

int GenerateYieldOverlayQuick()
{
    using namespace ARJ;

    auto TriggerLabelFromName = [](const std::string& trigName) -> std::string
    {
        if (trigName == kTriggerPP)
            return "Trigger: Photon 4 GeV + MBD NS #geq 1";
        int photonPt = 0;
        if (std::sscanf(trigName.c_str(), "photon_%d_plus", &photonPt) == 1)
            return TString::Format("Trigger: Photon %d GeV + MBD NS #geq 2, vtx < 150 cm", photonPt).Data();
        if (trigName.find("MBD_NS_geq_2_vtx_lt_150") != std::string::npos)
            return "Trigger: MBD NS #geq 2, vtx < 150 cm";
        return "Trigger: " + trigName;
    };

    const string isoConeLabel = (kAA_IsoConeR == "isoR40")
        ? "#DeltaR_{cone} < 0.4" : "#DeltaR_{cone} < 0.3";
    string isoModeLabel;
    if (kAA_IsoMode == "fixedIso5GeV") isoModeLabel = "E_{T}^{iso} < 5 GeV";
    else                               isoModeLabel = "Sliding iso cut";
    const string vzLabel = TString::Format("|v_{z}| < %d cm", kAA_VzCut).Data();

    TFile* fAA = TFile::Open(InputAuAu().c_str(), "READ");
    if (!fAA || fAA->IsZombie())
    {
        if (fAA) { fAA->Close(); delete fAA; }
        std::cerr << "[FATAL] Cannot open AuAu input: " << InputAuAu() << "\n";
        return 1;
    }

    TFile* fPP = TFile::Open(InputPP(isRun25pp).c_str(), "READ");
    if (!fPP || fPP->IsZombie())
    {
        if (fPP) { fPP->Close(); delete fPP; }
        fAA->Close();
        delete fAA;
        std::cerr << "[FATAL] Cannot open PP input: " << InputPP(isRun25pp) << "\n";
        return 1;
    }

    auto RawSignalYield = [](double A, double B, double C, double D,
                             double& Asig, double& eAsig) -> bool
    {
        Asig = 0.0;
        eAsig = 0.0;
        if (!(A > 0.0) || !(D > 0.0)) return false;

        Asig = A - B * (C / D);
        if (Asig < 0.0) Asig = 0.0;

        const double dSdA = 1.0;
        const double dSdB = -(C / D);
        const double dSdC = -(B / D);
        const double dSdD =  (B * C) / (D * D);

        double var = 0.0;
        if (A > 0.0) var += dSdA * dSdA * A;
        if (B > 0.0) var += dSdB * dSdB * B;
        if (C > 0.0) var += dSdC * dSdC * C;
        if (D > 0.0) var += dSdD * dSdD * D;
        eAsig = (var > 0.0) ? std::sqrt(var) : 0.0;
        return (Asig > 0.0);
    };

    struct SelCent
    {
        int lo;
        int hi;
        int color;
        std::vector<std::string> suffixes;
    };

    std::vector<SelCent> selCents;
    {
        const auto& overlayCentBins = OverlayCentBins();
        const int selCentColors[] = {kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2, kRed + 1};
        for (int i = 0; i < (int)overlayCentBins.size(); ++i)
        {
            SelCent sc;
            sc.lo = overlayCentBins[i].lo;
            sc.hi = overlayCentBins[i].hi;
            sc.color = selCentColors[i % 6];
            sc.suffixes = overlayCentBins[i].suffixes;
            selCents.push_back(sc);
        }
    }

    struct RegionStyle
    {
        string key;
        string label;
        int color;
        int marker;
    };

    const std::vector<RegionStyle> kRegionStyles = {
        {"A", "Region A", kRed + 1, 20},
        {"B", "Region B", kBlue + 1, 21},
        {"C", "Region C", kGreen + 2, 22},
        {"D", "Region D", kMagenta + 1, 33},
    };

    auto RegionYieldFromPP = [&](TDirectory* ppDir, const string& regionKey, int ptIdx) -> double
    {
        if (!ppDir) return 0.0;
        const PtBin& b = PtBins()[ptIdx];
        string histName;
        if      (regionKey == "A") histName = "h_isIsolated_isTight"   + b.suffix;
        else if (regionKey == "B") histName = "h_notIsolated_isTight"  + b.suffix;
        else if (regionKey == "C") histName = "h_isIsolated_notTight"  + b.suffix;
        else if (regionKey == "D") histName = "h_notIsolated_notTight" + b.suffix;
        else return 0.0;
        TH1* h = dynamic_cast<TH1*>(ppDir->Get(histName.c_str()));
        return h ? h->GetBinContent(1) : 0.0;
    };

    auto RegionYieldFromAA = [&](TDirectory* trigDir, const std::vector<std::string>& suffixes,
                                 const string& regionKey, int ptIdx) -> double
    {
        if (!trigDir) return 0.0;
        const PtBin& b = PtBins()[ptIdx];
        double sum = 0.0;
        for (const auto& suf : suffixes)
        {
            string histName;
            if      (regionKey == "A") histName = "h_isIsolated_isTight"   + b.suffix + suf;
            else if (regionKey == "B") histName = "h_notIsolated_isTight"  + b.suffix + suf;
            else if (regionKey == "C") histName = "h_isIsolated_notTight"  + b.suffix + suf;
            else if (regionKey == "D") histName = "h_notIsolated_notTight" + b.suffix + suf;
            else continue;
            TH1* h = dynamic_cast<TH1*>(trigDir->Get(histName.c_str()));
            if (h) sum += h->GetBinContent(1);
        }
        return sum;
    };

    auto DrawRegionOverlay = [&](const string& outPath,
                                 const string& plotTitle,
                                 const string& trigLabel,
                                 const std::vector<string>& regionKeys,
                                 const std::function<double(const string&, int)>& YieldGetter,
                                 bool isPP) -> bool
    {
        static int overlayCounter = 0;
        ++overlayCounter;
        const TString canvasName = TString::Format("c_reg_overlay_quick_%d", overlayCounter);
        const TString frameName = TString::Format("hYieldRegionFrameQuick_%d", overlayCounter);

        TCanvas cReg(canvasName, canvasName, 900,700);
        ApplyCanvasMargins1D(cReg);
        cReg.SetLogy();

        TH1F hFrame(frameName,"",100, kPtEdges.front(), kPtEdges.back());
        hFrame.SetDirectory(nullptr);
        hFrame.SetStats(0);
        hFrame.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
        hFrame.GetYaxis()->SetTitle("Yield");

        TLegend leg(0.18, 0.17, 0.50, 0.31);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.SetTextFont(42);
        leg.SetTextSize(0.033);
        leg.SetNColumns(2);

        std::vector<TGraphErrors*> keep;
        double yMin = std::numeric_limits<double>::max();
        double yMax = 0.0;

        for (const auto& regionKey : regionKeys)
        {
            const auto* style = static_cast<const RegionStyle*>(nullptr);
            for (const auto& rs : kRegionStyles)
            {
                if (rs.key == regionKey)
                {
                    style = &rs;
                    break;
                }
            }
            if (!style) continue;

            std::vector<double> x;
            std::vector<double> ex;
            std::vector<double> y;
            std::vector<double> ey;

            for (int i = 0; i < kNPtBins; ++i)
            {
                const double yield = YieldGetter(regionKey, i);
                if (!(yield > 0.0)) continue;
                const double ptLo = kPtEdges[(std::size_t)i];
                const double ptHi = kPtEdges[(std::size_t)i + 1];
                x.push_back(0.5 * (ptLo + ptHi));
                ex.push_back(0.5 * (ptHi - ptLo));
                y.push_back(yield);
                ey.push_back(std::sqrt(yield));
                yMin = std::min(yMin, std::max(1e-6, yield - std::sqrt(yield)));
                yMax = std::max(yMax, yield + std::sqrt(yield));
            }

            if (x.empty()) continue;

            TGraphErrors* g = new TGraphErrors((int)x.size(), &x[0], &y[0], &ex[0], &ey[0]);
            g->SetLineWidth(2);
            g->SetLineColor(style->color);
            g->SetMarkerStyle(isPP ? style->marker + 4 : style->marker);
            g->SetMarkerSize(1.1);
            g->SetMarkerColor(style->color);
            keep.push_back(g);
            leg.AddEntry(g, style->label.c_str(), "pe");
        }

        if (keep.empty()) return false;

        if (!(yMin < std::numeric_limits<double>::max()) || !(yMin > 0.0)) yMin = 0.5;
        if (!(yMax > 0.0)) yMax = 10.0;

        hFrame.SetMinimum(std::max(0.5, 0.35 * yMin));
        hFrame.SetMaximum(7.0 * yMax);
        hFrame.Draw();

        for (auto* g : keep) g->Draw("P SAME");
        leg.Draw();

        TLatex tTitle;
        tTitle.SetNDC(true);
        tTitle.SetTextFont(42);
        tTitle.SetTextAlign(23);
        tTitle.SetTextSize(0.043);
        tTitle.DrawLatex(0.58, 0.975, plotTitle.c_str());

        TLatex tCuts;
        tCuts.SetNDC(true);
        tCuts.SetTextFont(42);
        tCuts.SetTextAlign(33);
        tCuts.SetTextSize(0.033);
        tCuts.DrawLatex(0.93, 0.90, trigLabel.c_str());
        tCuts.DrawLatex(0.93, 0.85, isoConeLabel.c_str());
        tCuts.DrawLatex(0.93, 0.80, isoModeLabel.c_str());
        tCuts.DrawLatex(0.93, 0.75, vzLabel.c_str());

        SaveCanvas(cReg, outPath);
        std::cout << "[WROTE] " << outPath << "\n";
        for (auto* g : keep) delete g;
        return true;
    };

    for (const auto& trigAA : kTriggersAuAu)
    {
        TDirectory* trigDir = fAA->GetDirectory(trigAA.c_str());
        TDirectory* ppDir = fPP->GetDirectory(kTriggerPP.c_str());
        if (!trigDir) continue;
        if (!ppDir) ppDir = fPP;

        const string outDirYield = JoinPath(JoinPath(OutputAuAu(), trigAA), "yieldOverlays");
        EnsureDir(outDirYield);
        const std::string trigLabelYield = TriggerLabelFromName(trigAA);
        const std::string trigLabelPP = TriggerLabelFromName(kTriggerPP);

        TCanvas cYield("c_yield_raw_centSelect_quick","c_yield_raw_centSelect_quick",900,700);
        ApplyCanvasMargins1D(cYield);
        cYield.SetLogy();

        TH1F hFrameYield("hYieldSelFrameQuick","",100, kPtEdges.front(), kPtEdges.back());
        hFrameYield.SetDirectory(nullptr);
        hFrameYield.SetStats(0);
        hFrameYield.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
        hFrameYield.GetYaxis()->SetTitle("Raw signal yield in region A");

        TLegend legYield(0.19, 0.17, 0.56, 0.31);
        legYield.SetBorderSize(0);
        legYield.SetFillStyle(0);
        legYield.SetTextFont(42);
        legYield.SetTextSize(0.033);
        legYield.SetNColumns(2);

        std::vector<TGraphErrors*> keepYield;
        double yMinYield = std::numeric_limits<double>::max();
        double yMaxYield = 0.0;

        auto UpdateYieldRange = [&](double y, double ey) -> void
        {
            if (!(y > 0.0)) return;
            const double yLo = std::max(1e-6, y - ey);
            const double yHi = y + ey;
            if (yLo > 0.0) yMinYield = std::min(yMinYield, yLo);
            yMaxYield = std::max(yMaxYield, yHi);
        };

        {
            std::vector<double> xPP, exPP, yPP, eyPP;
            for (int i = 0; i < kNPtBins; ++i)
            {
                const PtBin& b = PtBins()[i];
                auto Get1PP = [&](const string& hname)->double {
                    TH1* h = dynamic_cast<TH1*>(ppDir->Get(hname.c_str()));
                    return h ? h->GetBinContent(1) : 0.0;
                };

                const double A = Get1PP("h_isIsolated_isTight"   + b.suffix);
                const double B = Get1PP("h_notIsolated_isTight"  + b.suffix);
                const double C = Get1PP("h_isIsolated_notTight"  + b.suffix);
                const double D = Get1PP("h_notIsolated_notTight" + b.suffix);

                double Asig = 0.0, eAsig = 0.0;
                if (!RawSignalYield(A, B, C, D, Asig, eAsig)) continue;

                const double ptLo = kPtEdges[(std::size_t)i];
                const double ptHi = kPtEdges[(std::size_t)i + 1];
                xPP.push_back(0.5 * (ptLo + ptHi));
                exPP.push_back(0.5 * (ptHi - ptLo));
                yPP.push_back(Asig);
                eyPP.push_back(eAsig);
                UpdateYieldRange(Asig, eAsig);
            }

            if (!xPP.empty())
            {
                TGraphErrors* gPP = new TGraphErrors((int)xPP.size(), &xPP[0], &yPP[0], &exPP[0], &eyPP[0]);
                gPP->SetLineWidth(2);
                gPP->SetLineColor(kRed + 1);
                gPP->SetMarkerStyle(24);
                gPP->SetMarkerSize(1.1);
                gPP->SetMarkerColor(kRed + 1);
                keepYield.push_back(gPP);
                legYield.AddEntry(gPP, "pp", "pe");
            }
        }

        for (const auto& sc : selCents)
        {
            std::vector<double> xAA, exAA, yAA, eyAA;
            for (int i = 0; i < kNPtBins; ++i)
            {
                const PtBin& b = PtBins()[i];
                auto Get1AA = [&](const string& hname)->double {
                    TH1* h = dynamic_cast<TH1*>(trigDir->Get(hname.c_str()));
                    return h ? h->GetBinContent(1) : 0.0;
                };

                double A = 0.0, B = 0.0, C = 0.0, D = 0.0;
                for (const auto& suf : sc.suffixes)
                {
                    A += Get1AA("h_isIsolated_isTight"   + b.suffix + suf);
                    B += Get1AA("h_notIsolated_isTight"  + b.suffix + suf);
                    C += Get1AA("h_isIsolated_notTight"  + b.suffix + suf);
                    D += Get1AA("h_notIsolated_notTight" + b.suffix + suf);
                }

                double Asig = 0.0, eAsig = 0.0;
                if (!RawSignalYield(A, B, C, D, Asig, eAsig)) continue;

                const double ptLo = kPtEdges[(std::size_t)i];
                const double ptHi = kPtEdges[(std::size_t)i + 1];
                xAA.push_back(0.5 * (ptLo + ptHi));
                exAA.push_back(0.5 * (ptHi - ptLo));
                yAA.push_back(Asig);
                eyAA.push_back(eAsig);
                UpdateYieldRange(Asig, eAsig);
            }

            if (xAA.empty()) continue;

            TGraphErrors* gAA = new TGraphErrors((int)xAA.size(), &xAA[0], &yAA[0], &exAA[0], &eyAA[0]);
            gAA->SetLineWidth(2);
            gAA->SetLineColor(sc.color);
            gAA->SetMarkerStyle(20);
            gAA->SetMarkerSize(1.1);
            gAA->SetMarkerColor(sc.color);
            keepYield.push_back(gAA);
            legYield.AddEntry(gAA, TString::Format("AuAu %d-%d%%", sc.lo, sc.hi).Data(), "pe");
        }

        if (!keepYield.empty())
        {
            if (!(yMinYield < std::numeric_limits<double>::max()) || !(yMinYield > 0.0)) yMinYield = 0.5;
            if (!(yMaxYield > 0.0)) yMaxYield = 10.0;

            hFrameYield.SetMinimum(std::max(0.5, 0.40 * yMinYield));
            hFrameYield.SetMaximum(4.5 * yMaxYield);
            hFrameYield.Draw();

            for (auto* g : keepYield) g->Draw("P SAME");
            legYield.Draw();

            TLatex tTitleYield;
            tTitleYield.SetNDC(true);
            tTitleYield.SetTextFont(42);
            tTitleYield.SetTextAlign(23);
            tTitleYield.SetTextSize(0.043);
            tTitleYield.DrawLatex(0.56, 0.975,
                                  "ABCD Raw Signal Yield vs p_{T}^{#gamma} for each centrality, Run3auau");

            TLatex tCutsYield;
            tCutsYield.SetNDC(true);
            tCutsYield.SetTextFont(42);
            tCutsYield.SetTextAlign(23);
            tCutsYield.SetTextSize(0.033);
            tCutsYield.DrawLatex(0.67, 0.90, trigLabelYield.c_str());
            tCutsYield.DrawLatex(0.67, 0.85, isoConeLabel.c_str());
            tCutsYield.DrawLatex(0.67, 0.80, isoModeLabel.c_str());
            tCutsYield.DrawLatex(0.67, 0.75, vzLabel.c_str());

            const string outPath = JoinPath(outDirYield, "yield_rawSignal_centSelect_ppOverlay.png");
            SaveCanvas(cYield, outPath);
            std::cout << "[WROTE] " << outPath << "\n";
        }

        for (auto* g : keepYield) delete g;

        for (std::size_t iCent = 0; iCent < selCents.size(); ++iCent)
        {
            const auto& sc = selCents[iCent];
            const string centDir = JoinPath(outDirYield, TString::Format("cent_%d_%d", sc.lo, sc.hi).Data());
            EnsureDir(centDir);

            DrawRegionOverlay(
                JoinPath(centDir, "yield_regionsABCD_overlay.png"),
                TString::Format("ABCD Region Yields vs p_{T}^{#gamma}, AuAu %d-%d%%", sc.lo, sc.hi).Data(),
                trigLabelYield,
                {"A", "B", "C", "D"},
                [&](const string& regionKey, int ptIdx) { return RegionYieldFromAA(trigDir, sc.suffixes, regionKey, ptIdx); },
                false
            );

            DrawRegionOverlay(
                JoinPath(centDir, "yield_regionsAC_overlay.png"),
                TString::Format("AC Region Yields vs p_{T}^{#gamma}, AuAu %d-%d%%", sc.lo, sc.hi).Data(),
                trigLabelYield,
                {"A", "C"},
                [&](const string& regionKey, int ptIdx) { return RegionYieldFromAA(trigDir, sc.suffixes, regionKey, ptIdx); },
                false
            );

            if (iCent == 0)
            {
                DrawRegionOverlay(
                    JoinPath(centDir, "yield_regionsABCD_overlay_pp.png"),
                    "ABCD Region Yields vs p_{T}^{#gamma}, pp",
                    trigLabelPP,
                    {"A", "B", "C", "D"},
                    [&](const string& regionKey, int ptIdx) { return RegionYieldFromPP(ppDir, regionKey, ptIdx); },
                    true
                );

                DrawRegionOverlay(
                    JoinPath(centDir, "yield_regionsAC_overlay_pp.png"),
                    "AC Region Yields vs p_{T}^{#gamma}, pp",
                    trigLabelPP,
                    {"A", "C"},
                    [&](const string& regionKey, int ptIdx) { return RegionYieldFromPP(ppDir, regionKey, ptIdx); },
                    true
                );
            }
        }
    }

    fPP->Close();
    fAA->Close();
    delete fPP;
    delete fAA;
    return 0;
}
