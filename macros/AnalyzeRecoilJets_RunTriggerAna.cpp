void RunTriggerAna_DoNotScaleMaxClusterEnergy(Dataset& ds)
{
    if (ds.isSim) return;
    if (!ds.file || ds.file->IsZombie()) return;
    
    static std::set<std::string> s_doneOutDirs;
    
    // For AuAu (no centFolder): outBase = OutputAuAu()/<trigger>
    //   -> dirname = OutputAuAu() -> triggerQA lands in auau/triggerQA
    // For AuAu (with centFolder): outBase = OutputAuAu()/<trigger>/<cent>
    //   -> dirname(dirname) = OutputAuAu() -> triggerQA lands in auau/triggerQA
    // For PP: outBase = OutputPP()/<trigger>
    //   -> dirname = OutputPP() -> triggerQA lands in pp/triggerQA
    const std::string baseDir = ds.centFolder.empty()
    ? DirnameFromPath(ds.outBase)
    : DirnameFromPath(DirnameFromPath(ds.outBase));
    const std::string outDir  = JoinPath(baseDir, "triggerQA");
    
    if (s_doneOutDirs.count(outDir)) return;
    s_doneOutDirs.insert(outDir);
    
    EnsureDir(outDir);
    
    const std::string prefix        = "h_maxEnergyClus_NewTriggerFilling_doNotScale_";
    const std::string baselineTag   = "baseline_";
    
    struct RunSpan
    {
        uint64_t lo = 0;
        uint64_t hi = 0;
    };
    
    struct LoadedPair
    {
        std::string basisKey;
        std::string basisLabel;
        std::string groupFolder;
        std::string probeKey;
        std::string probeLabel;
        std::string baselineKey;
        std::string baselineLabel;
        std::string coverageText;
        int         nRunSpans = 0;
        TH1*        hProbe = nullptr;
        TH1*        hBase  = nullptr;
    };
    
    auto datasetBucket = [&]() -> std::string
    {
        const std::string hay = ds.inFilePath + " " + ds.outBase;
        if (hay.find("pp25") != std::string::npos) return "pp25";
        if (hay.find("RecoilJets_pp_ALL.root") != std::string::npos) return "pp24";
        if (hay.find("auau") != std::string::npos) return "auau";
        return "unknown";
    };
    
    auto extractPhotonThresholdGeV = [&](const std::string& key) -> int
    {
        std::smatch m;
        static const std::regex re("Photon_([0-9]+)_GeV");
        if (std::regex_search(key, m, re) && m.size() > 1)
        {
            return std::atoi(m[1].str().c_str());
        }
        return 999;
    };
    
    auto deduceBasisKey = [&](const std::string& probeKey) -> std::string
    {
        if (probeKey.find("_MBD_NS_geq_1_vtx_lt_10") != std::string::npos)   return "MBD_NS_geq_1_vtx_lt_10";
        if (probeKey.find("_MBD_NS_geq_2_vtx_lt_150") != std::string::npos)  return "MBD_NS_geq_2_vtx_lt_150";
        if (probeKey.find("_MBD_NS_geq_2_vtx_lt_10") != std::string::npos)   return "MBD_NS_geq_2_vtx_lt_10";
        if (probeKey.find("_MBD_NS_geq_2") != std::string::npos)             return "MBD_NS_geq_2";
        if (probeKey.find("_MBD_NS_geq_1") != std::string::npos)             return "MBD_NS_geq_1";
        return "";
    };
    
    auto basisLabelFromKey = [&](const std::string& basisKey) -> std::string
    {
        if (basisKey == "MBD_NS_geq_1")            return "MBD N&S >= 1";
        if (basisKey == "MBD_NS_geq_1_vtx_lt_10")  return "MBD N&S >= 1, vtx < 10 cm";
        if (basisKey == "MBD_NS_geq_2")            return "MBD N&S >= 2";
        if (basisKey == "MBD_NS_geq_2_vtx_lt_10")  return "MBD N&S >= 2, vtx < 10 cm";
        if (basisKey == "MBD_NS_geq_2_vtx_lt_150") return "MBD N&S >= 2, vtx < 150 cm";
        return basisKey;
    };
    
    auto groupFolderFromBasis = [&](const std::string& basisKey) -> std::string
    {
        if (basisKey == "MBD_NS_geq_1")            return "commonBasis_MBD_NandS_geq_1";
        if (basisKey == "MBD_NS_geq_1_vtx_lt_10")  return "commonBasis_MBD_NandS_geq_1_vtx_lt_10";
        if (basisKey == "MBD_NS_geq_2")            return "commonBasis_MBD_NS_geq_2";
        if (basisKey == "MBD_NS_geq_2_vtx_lt_10")  return "commonBasis_MBD_NS_geq_2_vtx_lt_10";
        if (basisKey == "MBD_NS_geq_2_vtx_lt_150") return "commonBasis_MBD_NS_geq_2_vtx_lt_150";
        return "commonBasis_UNKNOWN";
    };
    
    auto probeLabelFromKey = [&](const std::string& probeKey) -> std::string
    {
        const int thr = extractPhotonThresholdGeV(probeKey);
        const std::string basisKey = deduceBasisKey(probeKey);
        const std::string basisLabel = basisLabelFromKey(basisKey);
        
        if (thr >= 0 && !basisLabel.empty())
        {
            return TString::Format("Photon %d GeV + %s", thr, basisLabel.c_str()).Data();
        }
        
        return probeKey;
    };
    
    auto colorForProbe = [&](const std::string& probeKey) -> int
    {
        const int thr = extractPhotonThresholdGeV(probeKey);
        if (thr <= 2)  return kOrange+7;
        if (thr == 3)  return kBlue+1;
        if (thr == 4)  return kGreen+2;
        if (thr == 5)  return kMagenta+1;
        if (thr == 6)  return kRed+1;
        if (thr == 8)  return kAzure+2;
        if (thr == 10) return kBlue+1;
        if (thr == 12) return kRed+1;
        if (thr == 14) return kOrange+1;
        if (thr == 18) return kGreen+3;
        if (thr == 20) return kPink+1;
        return kRed+1;
    };
    
    auto markerForProbe = [&](const std::string& probeKey) -> int
    {
        const int thr = extractPhotonThresholdGeV(probeKey);
        if (thr <= 2)  return 20;
        if (thr == 3)  return 21;
        if (thr == 4)  return 22;
        if (thr == 5)  return 33;
        if (thr == 6)  return 29;
        if (thr == 8)  return 34;
        if (thr == 10) return 23;
        if (thr == 12) return 24;
        if (thr == 14) return 25;
        if (thr == 18) return 26;
        if (thr == 20) return 32;
        return 20;
    };
    
    auto getHist = [&](const std::string& dirKey)->TH1*
    {
        if (!ds.file) return nullptr;
        
        TDirectory* dir = ds.file->GetDirectory(dirKey.c_str());
        if (!dir) return nullptr;
        
        const std::string hname = prefix + dirKey;
        return dynamic_cast<TH1*>(dir->Get(hname.c_str()));
    };
    
    auto spansForProbe = [&](const std::string& probeKey)->std::vector<RunSpan>
    {
        std::vector<RunSpan> spans;
        
        const std::string bucket = datasetBucket();
        const int thr = extractPhotonThresholdGeV(probeKey);
        const bool geq1   = (probeKey.find("_MBD_NS_geq_1") != std::string::npos);
        const bool geq2   = (probeKey.find("_MBD_NS_geq_2") != std::string::npos);
        const bool vtx10  = (probeKey.find("_vtx_lt_10")  != std::string::npos);
        const bool vtx150 = (probeKey.find("_vtx_lt_150") != std::string::npos);
        
        auto add = [&](uint64_t lo, uint64_t hi)
        {
            spans.push_back({lo, hi});
        };
        
        if (bucket == "pp24")
        {
            if (geq1 && !vtx10 && !vtx150)
            {
                if (thr >= 2 && thr <= 5)
                {
                    add(47289, 51015);
                    add(51093, 52596);
                    add(52610, 53864);
                }
            }
            else if (geq1 && vtx10)
            {
                if (thr >= 3 && thr <= 5)
                {
                    add(51093, 52596);
                    add(52610, 53864);
                }
            }
        }
        else if (bucket == "pp25")
        {
            if (geq1 && !vtx10 && !vtx150)
            {
                if (thr >= 2 && thr <= 5) add(79269, 81664);
            }
            else if (geq1 && vtx10)
            {
                if (thr == 2) add(79495, 81664);
                if (thr >= 3 && thr <= 5) add(79269, 81664);
            }
        }
        else if (bucket == "auau")
        {
            if (geq2 && !vtx10 && !vtx150)
            {
                if (thr >= 6 && thr <= 12 && (thr % 2 == 0)) add(67599, 68155);
                if (thr >= 2 && thr <= 4) add(78686, 78686);
            }
            else if (geq2 && vtx10)
            {
                if (thr >= 6 && thr <= 12 && (thr % 2 == 0))
                {
                    add(68208, 68220);
                    add(68335, 69616);
                    add(71328, 78572);
                    add(78689, 78954);
                }
                if (thr == 3 || thr == 8 || thr == 10) add(78686, 78686);
            }
            else if (geq2 && vtx150)
            {
                if (thr >= 6 && thr <= 12 && (thr % 2 == 0))
                {
                    add(68208, 68220);
                    add(68335, 69616);
                    add(71328, 78572);
                    add(78689, 78954);
                }
                if (thr == 8 || thr == 10) add(78686, 78686);
            }
        }
        
        return spans;
    };
    
    auto coverageTextForProbe = [&](const std::string& probeKey)->std::string
    {
        const auto spans = spansForProbe(probeKey);
        if (spans.empty())
        {
            return "pair-specific coverage encoded by baseline_<probe>; exact run spans not written to ROOT";
        }
        
        std::ostringstream os;
        for (std::size_t i = 0; i < spans.size(); ++i)
        {
            if (i) os << ", ";
            os << spans[i].lo << "-" << spans[i].hi;
        }
        return os.str();
    };
    
    auto runSpanCountForProbe = [&](const std::string& probeKey)->int
    {
        return static_cast<int>(spansForProbe(probeKey).size());
    };
    
    auto findXAtEff = [&](TH1* h, double target)->double
    {
        if (!h) return -1.0;
        for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
        {
            const double y = h->GetBinContent(ib);
            if (!std::isfinite(y)) continue;
            if (y >= target) return h->GetBinCenter(ib);
        }
        return -1.0;
    };
    
    const std::vector<std::string> groupOrder = {
        "commonBasis_MBD_NandS_geq_1",
        "commonBasis_MBD_NandS_geq_1_vtx_lt_10",
        "commonBasis_MBD_NS_geq_2",
        "commonBasis_MBD_NS_geq_2_vtx_lt_10",
        "commonBasis_MBD_NS_geq_2_vtx_lt_150"
    };
    
    auto groupRank = [&](const std::string& folder) -> int
    {
        for (std::size_t i = 0; i < groupOrder.size(); ++i)
        {
            if (groupOrder[i] == folder) return static_cast<int>(i);
        }
        return 999;
    };
    
    auto collectPairsFromFile = [&](TFile* inFile, const std::string& inFilePath, bool verboseWarnings) -> std::vector<LoadedPair>
    {
        std::vector<LoadedPair> outPairs;
        if (!inFile) return outPairs;
        
        auto getHistFromFile = [&](const std::string& dirKey)->TH1*
        {
            TDirectory* dir = inFile->GetDirectory(dirKey.c_str());
            if (!dir) return nullptr;
            
            const std::string hname = prefix + dirKey;
            return dynamic_cast<TH1*>(dir->Get(hname.c_str()));
        };
        
        TIter nextKey(inFile->GetListOfKeys());
        while (TKey* key = dynamic_cast<TKey*>(nextKey()))
        {
            const std::string dirKey = key->GetName();
            if (dirKey.empty()) continue;
            if (dirKey.find(baselineTag) == 0) continue;
            
            TH1* hProbe = getHistFromFile(dirKey);
            if (!hProbe) continue;
            
            const std::string baselineKey = baselineTag + dirKey;
            TH1* hBase = getHistFromFile(baselineKey);
            if (!hBase)
            {
                if (verboseWarnings)
                {
                    cout << ANSI_BOLD_YEL
                    << "[WARN] Missing doNotScale baseline histogram for probe " << dirKey << "\n"
                    << "       Need: " << baselineKey << "/" << prefix << baselineKey << "\n"
                    << ANSI_RESET << "\n";
                }
                continue;
            }
            
            const std::string basisKey = deduceBasisKey(dirKey);
            if (basisKey.empty())
            {
                if (verboseWarnings)
                {
                    cout << ANSI_BOLD_YEL
                    << "[WARN] Could not deduce common baseline group for probe " << dirKey
                    << ANSI_RESET << "\n";
                }
                continue;
            }
            
            LoadedPair P;
            P.basisKey      = basisKey;
            P.basisLabel    = basisLabelFromKey(basisKey);
            P.groupFolder   = groupFolderFromBasis(basisKey);
            P.probeKey      = dirKey;
            P.probeLabel    = probeLabelFromKey(dirKey);
            P.baselineKey   = baselineKey;
            P.baselineLabel = TString::Format("Baseline raw %s", P.basisLabel.c_str()).Data();
            P.coverageText  = coverageTextForProbe(dirKey);
            P.nRunSpans     = runSpanCountForProbe(dirKey);
            P.hProbe        = hProbe;
            P.hBase         = hBase;
            outPairs.push_back(P);
        }
        
        std::sort(outPairs.begin(), outPairs.end(), [&](const LoadedPair& a, const LoadedPair& b)
                  {
            const int ra = groupRank(a.groupFolder);
            const int rb = groupRank(b.groupFolder);
            if (ra != rb) return ra < rb;
            
            const int ta = extractPhotonThresholdGeV(a.probeKey);
            const int tb = extractPhotonThresholdGeV(b.probeKey);
            if (ta != tb) return ta < tb;
            
            return a.probeKey < b.probeKey;
        });
        
        if (outPairs.empty() && verboseWarnings)
        {
            cout << ANSI_BOLD_YEL
            << "[WARN] No pairwise doNotScale numerator/baseline histogram pairs found in " << inFilePath
            << ANSI_RESET << "\n";
        }
        
        return outPairs;
    };
    
    std::vector<LoadedPair> pairs = collectPairsFromFile(ds.file, ds.inFilePath, true);
    if (pairs.empty())
    {
        return;
    }
    
    std::sort(pairs.begin(), pairs.end(), [&](const LoadedPair& a, const LoadedPair& b)
              {
        const int ra = groupRank(a.groupFolder);
        const int rb = groupRank(b.groupFolder);
        if (ra != rb) return ra < rb;
        
        const int ta = extractPhotonThresholdGeV(a.probeKey);
        const int tb = extractPhotonThresholdGeV(b.probeKey);
        if (ta != tb) return ta < tb;
        
        return a.probeKey < b.probeKey;
    });
    
    std::map<std::string, std::vector<LoadedPair>> groups;
    for (const auto& P : pairs) groups[P.groupFolder].push_back(P);
    
    std::vector<std::string> indexLines;
    indexLines.push_back(string("triggerQA doNotScale summary (") + ds.label + ")");
    indexLines.push_back(string("Input: ") + ds.inFilePath);
    indexLines.push_back("Schema: numerator=<probe>, denominator=baseline_<probe>");
    indexLines.push_back("Efficiency definition: probe live / baseline raw");
    indexLines.push_back("");
    
    auto drawPairOverlay = [&](const LoadedPair& P, const std::string& pairOutDir)
    {
        TH1* hBase = CloneTH1(P.hBase, TString::Format("hBase_pairOverlay_%s", P.probeKey.c_str()).Data());
        TH1* hProbe = CloneTH1(P.hProbe, TString::Format("hProbe_pairOverlay_%s", P.probeKey.c_str()).Data());
        if (!hBase || !hProbe)
        {
            if (hBase) delete hBase;
            if (hProbe) delete hProbe;
            return;
        }
        
        TCanvas c(TString::Format("c_trigAna_pairOverlay_%s", P.probeKey.c_str()).Data(),
                  TString::Format("c_trigAna_pairOverlay_%s", P.probeKey.c_str()).Data(),
                  900, 700);
        c.cd();
        c.SetLeftMargin(0.14);
        c.SetRightMargin(0.05);
        c.SetBottomMargin(0.14);
        c.SetTopMargin(0.08);
        c.SetTicks(1,1);
        c.SetLogy();
        
        const double ymax = std::max(hBase->GetMaximum(), hProbe->GetMaximum());
        double minPos = SmallestPositiveBinContent(hBase);
        const double probeMinPos = SmallestPositiveBinContent(hProbe);
        if (probeMinPos > 0.0)
        {
            if (minPos > 0.0) minPos = std::min(minPos, probeMinPos);
            else              minPos = probeMinPos;
        }
        
        hBase->SetTitle("");
        hBase->GetXaxis()->SetTitle("Maximum Cluster Energy [GeV]");
        hBase->GetYaxis()->SetTitle("Counts");
        hBase->GetXaxis()->SetTitleSize(0.055);
        hBase->GetXaxis()->SetTitleOffset(1.05);
        hBase->GetXaxis()->SetLabelSize(0.045);
        hBase->GetYaxis()->SetTitleSize(0.055);
        hBase->GetYaxis()->SetTitleOffset(1.20);
        hBase->GetYaxis()->SetLabelSize(0.045);
        hBase->GetXaxis()->SetRangeUser(0.0, 20.0);
        hBase->SetMinimum((minPos > 0.0) ? (0.5 * minPos) : 1e-6);
        hBase->SetMaximum((ymax > 0.0) ? (1.25 * ymax) : 1.0);
        hBase->SetLineColor(kBlack);
        hBase->SetLineWidth(4);
        
        hProbe->SetLineColor(colorForProbe(P.probeKey));
        hProbe->SetLineWidth(4);
        
        hBase->Draw("HIST");
        hProbe->Draw("HIST SAME");
        
        TLegend leg(0.42, 0.62, 0.90, 0.88);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.028);
        leg.AddEntry(hBase, P.baselineLabel.c_str(), "l");
        leg.AddEntry(hProbe, TString::Format("Probe live %s", P.probeLabel.c_str()).Data(), "l");
        leg.Draw();
        
        std::vector<std::string> lines = DefaultHeaderLines(ds);
        lines.push_back(P.basisLabel);
        lines.push_back("Pairwise doNotScale overlay");
        lines.push_back(TString::Format("Run spans (%d): %s", P.nRunSpans, P.coverageText.c_str()).Data());
        DrawLatexLines(0.18, 0.30, lines, 0.028, 0.040);
        
        const std::string outPng = JoinPath(pairOutDir, "hMaxClusterEnergy_pairOverlay.png");
        SaveCanvas(c, outPng);
        cout << ANSI_BOLD_GRN << "[WROTE] " << outPng << ANSI_RESET << "\n";
        
        delete hBase;
        delete hProbe;
    };
    
    auto drawPairTurnOn = [&](const LoadedPair& P, const std::string& pairOutDir, double& x95Out)
    {
        x95Out = -1.0;
        
        TH1* hBase = CloneTH1(P.hBase, TString::Format("hBase_pairTurnOn_%s", P.probeKey.c_str()).Data());
        TH1* hProbe = CloneTH1(P.hProbe, TString::Format("hProbe_pairTurnOn_%s", P.probeKey.c_str()).Data());
        if (!hBase || !hProbe)
        {
            if (hBase) delete hBase;
            if (hProbe) delete hProbe;
            return;
        }
        
        TH1* hRatio = CloneTH1(hProbe, TString::Format("hRatio_pairTurnOn_%s", P.probeKey.c_str()).Data());
        if (!hRatio)
        {
            delete hBase;
            delete hProbe;
            return;
        }
        
        EnsureSumw2(hBase);
        EnsureSumw2(hProbe);
        EnsureSumw2(hRatio);
        hRatio->SetDirectory(nullptr);
        hRatio->Divide(hProbe, hBase, 1.0, 1.0, "B");
        x95Out = findXAtEff(hRatio, 0.95);
        
        TCanvas c(TString::Format("c_trigAna_pairTurnOn_%s", P.probeKey.c_str()).Data(),
                  TString::Format("c_trigAna_pairTurnOn_%s", P.probeKey.c_str()).Data(),
                  900, 700);
        c.cd();
        c.SetLeftMargin(0.14);
        c.SetRightMargin(0.05);
        c.SetBottomMargin(0.14);
        c.SetTopMargin(0.08);
        c.SetTicks(1,1);
        
        hRatio->SetTitle("");
        hRatio->GetXaxis()->SetTitle("Maximum Cluster Energy [GeV]");
        hRatio->GetYaxis()->SetTitle("Probe live / baseline raw");
        hRatio->GetXaxis()->SetTitleSize(0.055);
        hRatio->GetXaxis()->SetTitleOffset(1.05);
        hRatio->GetXaxis()->SetLabelSize(0.045);
        hRatio->GetYaxis()->SetTitleSize(0.055);
        hRatio->GetYaxis()->SetTitleOffset(1.20);
        hRatio->GetYaxis()->SetLabelSize(0.045);
        hRatio->GetXaxis()->SetRangeUser(0.0, 20.0);
        hRatio->SetMinimum(0.0);
        hRatio->SetMaximum(1.20);
        hRatio->SetMarkerStyle(markerForProbe(P.probeKey));
        hRatio->SetMarkerSize(1.0);
        hRatio->SetMarkerColor(colorForProbe(P.probeKey));
        hRatio->SetLineColor(colorForProbe(P.probeKey));
        hRatio->SetLineWidth(3);
        hRatio->Draw("E1");
        
        TLine l1(0.0, 1.0, 20.0, 1.0);
        l1.SetLineStyle(2);
        l1.SetLineWidth(2);
        l1.SetLineColor(kBlack);
        l1.Draw("SAME");
        
        TLine l95(0.0, 0.95, 20.0, 0.95);
        l95.SetLineStyle(3);
        l95.SetLineWidth(2);
        l95.SetLineColor(kGray+2);
        l95.Draw("SAME");
        
        TLegend leg(0.18, 0.68, 0.64, 0.88);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.028);
        if (x95Out > 0.0)
        {
            leg.AddEntry(hRatio,
                         TString::Format("%s (95%% = %.2f GeV)", P.probeLabel.c_str(), x95Out).Data(),
                         "pe");
        }
        else
        {
            leg.AddEntry(hRatio, P.probeLabel.c_str(), "pe");
        }
        leg.Draw();
        
        std::vector<std::string> lines = DefaultHeaderLines(ds);
        lines.push_back(P.basisLabel);
        lines.push_back("Pairwise doNotScale turn-on");
        lines.push_back(TString::Format("Run spans (%d): %s", P.nRunSpans, P.coverageText.c_str()).Data());
        DrawLatexLines(0.18, 0.30, lines, 0.028, 0.040);
        
        const std::string outPng = JoinPath(pairOutDir, "hMaxClusterEnergy_pairTurnOn.png");
        SaveCanvas(c, outPng);
        cout << ANSI_BOLD_GRN << "[WROTE] " << outPng << ANSI_RESET << "\n";
        
        delete hBase;
        delete hProbe;
        delete hRatio;
    };
    
    auto drawGroupTurnOnOverlay = [&](const std::vector<LoadedPair>& loaded, const std::string& groupOutDir, const std::string& filenameTag = "", const std::string& perRunLabel = "")
    {
        std::vector<TH1*> ratioHists;
        std::vector<TF1*> fitFuncs;
        std::vector<double> x80s;
        std::vector<double> plateauVals;
        ratioHists.reserve(loaded.size());
        fitFuncs.reserve(loaded.size());
        x80s.reserve(loaded.size());
        plateauVals.reserve(loaded.size());
        
        // ── helper: linear-interpolation crossing finder ──
        auto findCrossingX = [&](TH1* h, double target, double xLo, double xHi) -> double
        {
            if (!h) return -1.0;
            
            const int nb = h->GetNbinsX();
            for (int ib = 1; ib < nb; ++ib)
            {
                const double x1 = h->GetBinCenter(ib);
                const double x2 = h->GetBinCenter(ib + 1);
                if (x2 < xLo || x1 > xHi) continue;
                
                const double y1 = h->GetBinContent(ib);
                const double y2 = h->GetBinContent(ib + 1);
                if (!std::isfinite(y1) || !std::isfinite(y2)) continue;
                
                if (y1 == target) return x1;
                if (y2 == target) return x2;
                
                if ((y1 - target) * (y2 - target) < 0.0 && std::fabs(y2 - y1) > 1e-12)
                {
                    const double frac = (target - y1) / (y2 - y1);
                    if (frac >= 0.0 && frac <= 1.0) return x1 + frac * (x2 - x1);
                }
            }
            
            return -1.0;
        };
        
        // ── helper: weighted mean of tail bins (fallback when pol0 fit has too few bins) ──
        auto fallbackTailMean = [&](TH1* h, double xTailLo) -> double
        {
            if (!h) return 0.94;
            
            const int bLo = h->GetXaxis()->FindBin(xTailLo + 1e-6);
            const int bHi = h->GetXaxis()->FindBin(20.0 - 1e-6);
            
            double sumW = 0.0;
            double sumWY = 0.0;
            
            for (int ib = bLo; ib <= bHi; ++ib)
            {
                const double y  = h->GetBinContent(ib);
                const double ey = h->GetBinError(ib);
                if (!std::isfinite(y) || y <= 0.0) continue;
                
                double w = 1.0;
                if (std::isfinite(ey) && ey > 0.0) w = 1.0 / (ey * ey);
                
                sumW  += w;
                sumWY += w * y;
            }
            
            if (sumW <= 0.0)
            {
                const double ymax = h->GetMaximum();
                const double fallback = (std::isfinite(ymax) && ymax > 0.0) ? ymax : 0.94;
                return std::max(0.80, std::min(1.0, fallback));
            }
            
            return std::max(0.80, std::min(1.0, sumWY / sumW));
        };
        
        // ── helper: fit a pol0 constant to the tail region [xTailLo, 20] ──
        auto fitPlateauConst = [&](TH1* h, double xTailLo) -> double
        {
            if (!h) return 0.94;
            
            xTailLo = std::max(9.0, std::min(19.0, xTailLo));
            
            const int bLo = h->GetXaxis()->FindBin(xTailLo + 1e-6);
            const int bHi = h->GetXaxis()->FindBin(20.0 - 1e-6);
            if (bHi - bLo + 1 < 3) return fallbackTailMean(h, xTailLo);
            
            static int plateauCounter = 0;
            TF1 fPlateau(TString::Format("f_plateau_tmp_%d", plateauCounter++).Data(), "pol0", xTailLo, 20.0);
            fPlateau.SetLineWidth(2);
            
            TFitResultPtr r = h->Fit(&fPlateau, "RQ0NS", "", xTailLo, 20.0);
            const int status = r;
            
            double plateau = fallbackTailMean(h, xTailLo);
            if (status == 0)
            {
                const double p0 = fPlateau.GetParameter(0);
                if (std::isfinite(p0) && p0 > 0.0) plateau = p0;
            }
            
            return std::max(0.80, std::min(1.0, plateau));
        };
        
        // ── helper: build a Gumbel CDF fit function for turn-on curves ──
        //
        // Gumbel CDF:  F(x) = P * exp( -exp( -(x - mu) / beta ) )
        //
        // This is ASYMMETRIC — sharp rise from zero, gradual approach to
        // plateau — which matches trigger turn-on physics far better than
        // the symmetric erf.
        //
        // Parameters:
        //   [0] = P      (plateau, fixed from tail fit)
        //   [1] = mu     (location: F(mu) = P * e^{-1} ≈ 0.368*P)
        //   [2] = beta   (scale: controls steepness)
        //
        // Gumbel CDF math for initial guesses:
        //   At fraction f of plateau: x_f = mu - beta * ln(-ln(f))
        //   50% point: x50 = mu - beta * ln(ln2)    [ln(ln2) ≈ -0.3665]
        //              => mu ≈ x50 + 0.3665*beta
        //   80% point: x80 = mu + 1.500*beta
        //   20% point: x20 = mu - 0.476*beta
        //              => beta ≈ (x80 - x20) / 1.976
        //
        // Drawing range [0, 20]: F(0) ≈ 0 for mu ≈ 5–8, so curve hugs
        // the x-axis at the left edge with no visible "slop."
        auto buildGumbelFit = [&](const LoadedPair& P, TH1* hRatio, double plateauVal) -> TF1*
        {
            const int thr = extractPhotonThresholdGeV(P.probeKey);
            
            // Default guesses if data crossings aren't found
            const double defaultMu =
            (thr == 10) ? 4.8 :
            (thr == 12) ? 6.2 :
            std::max(2.0, 0.55 * ((thr > 0) ? static_cast<double>(thr) : 8.0));
            
            const double defaultBeta =
            (thr == 10) ? 1.1 :
            (thr == 12) ? 1.3 :
            1.2;
            
            // Find data crossing points relative to the plateau
            const double x20rel = findCrossingX(hRatio, 0.20 * plateauVal, 1.0, 15.0);
            const double x50rel = findCrossingX(hRatio, 0.50 * plateauVal, 1.0, 15.0);
            const double x80rel = findCrossingX(hRatio, 0.80 * plateauVal, 1.0, 15.0);
            
            // Estimate beta from 20%–80% rise width
            double betaGuess = defaultBeta;
            if (std::isfinite(x20rel) && std::isfinite(x80rel) && x80rel > x20rel)
            {
                betaGuess = (x80rel - x20rel) / 1.976;
            }
            betaGuess = std::max(0.4, std::min(4.0, betaGuess));
            
            // Estimate mu from 50% crossing: mu = x50 + 0.3665*beta
            double muGuess = defaultMu;
            if (std::isfinite(x50rel) && x50rel > 0.0)
            {
                muGuess = x50rel + 0.3665 * betaGuess;
            }
            
            // Gumbel CDF: P * exp(-exp(-(x-mu)/beta))
            // Drawing range [0, 20]: at x=0 with mu=5, beta=1.1:
            //   F(0) = P * exp(-exp(5/1.1)) = P * exp(-95) ≈ 0
            TF1* fGumbel = new TF1(
                                   TString::Format("f_groupTurnOn_%s", P.probeKey.c_str()).Data(),
                                   "[0]*exp(-exp(-((x-[1])/[2])))",
                                   0.0, 20.0
                                   );
            
            if (!fGumbel) return nullptr;
            
            fGumbel->SetParNames("Plateau", "Mu", "Beta");
            fGumbel->SetParameter(0, plateauVal);
            
            if (generatePerRunTriggerAnaForGoodRunsOnly)
            {
                // Good-runs-only mode: fix plateau at 1.0 for all triggers.
                // The tail estimate is unreliable with few runs, and we want
                // to verify the selected runs produce a true turn-on to unity.
                fGumbel->FixParameter(0, 1.0);
            }
            // Photon 10: fix plateau hard at the tail estimate — the data
            // clearly flattens at ~0.94 and any upward float overshoots.
            // Also restrict beta to be small → sharper knee into the plateau.
            // Photon 12+: allow ±0.02 float as before.
            else if (thr <= 10)
            {
                fGumbel->FixParameter(0, plateauVal);
            }
            else
            {
                // Allow plateau to nudge up slightly (+0.015) but not down
                fGumbel->SetParLimits(0, plateauVal, std::min(1.0, plateauVal + 0.01));
            }
            
            fGumbel->SetParameter(1, muGuess);
            fGumbel->SetParameter(2, betaGuess);
            fGumbel->SetParLimits(1, std::max(1.0, muGuess - 3.0), std::min(18.0, muGuess + 3.0));
            
            if (generatePerRunTriggerAnaForGoodRunsOnly)
            {
                // Good-runs-only: relax beta limits with plateau fixed at 1.0
                if (thr <= 10)
                    fGumbel->SetParLimits(2, std::max(0.3, 0.40 * betaGuess), std::max(3.5, 2.0 * betaGuess));
                else
                    fGumbel->SetParLimits(2, std::max(0.3, 0.40 * betaGuess), std::max(2.0, 1.1 * betaGuess));
            }
            else if (thr <= 10)
            {
                // Tighter beta: force a steeper rise and sharper flattening
                fGumbel->SetParLimits(2, std::max(0.3, 0.50 * betaGuess), std::min(1.7, 1.05 * betaGuess));
            }
            else
            {
                // Tighter beta for Photon 12: sharper knee into plateau
                fGumbel->SetParLimits(2, std::max(0.3, 0.50 * betaGuess), std::min(2.1, 1.15 * betaGuess));
            }
            fGumbel->SetNpx(500);
            fGumbel->SetLineColor(colorForProbe(P.probeKey));
            fGumbel->SetLineStyle(2);
            fGumbel->SetLineWidth(3);
            return fGumbel;
        };
        
        // ── main loop: build ratio histograms + two-pass Gumbel fit per probe ──
        for (const auto& P : loaded)
        {
            TH1* hBase = CloneTH1(P.hBase, TString::Format("hBase_groupTurnOn_%s", P.probeKey.c_str()).Data());
            TH1* hProbe = CloneTH1(P.hProbe, TString::Format("hProbe_groupTurnOn_%s", P.probeKey.c_str()).Data());
            if (!hBase || !hProbe)
            {
                if (hBase) delete hBase;
                if (hProbe) delete hProbe;
                continue;
            }
            
            TH1* hRatio = CloneTH1(hProbe, TString::Format("hRatio_groupTurnOn_%s", P.probeKey.c_str()).Data());
            if (!hRatio)
            {
                delete hBase;
                delete hProbe;
                continue;
            }
            
            EnsureSumw2(hBase);
            EnsureSumw2(hProbe);
            EnsureSumw2(hRatio);
            hRatio->Divide(hProbe, hBase, 1.0, 1.0, "B");
            hRatio->SetDirectory(nullptr);
            hRatio->SetMarkerStyle(20);
            hRatio->SetMarkerSize(1.0);
            hRatio->SetMarkerColor(colorForProbe(P.probeKey));
            hRatio->SetLineColor(colorForProbe(P.probeKey));
            hRatio->SetLineWidth(3);
            
            // ────────────────────────────────────────────────────────────
            // TWO-PASS FIT using Gumbel CDF
            //
            // Pass 1: estimate plateau from [15, 20] GeV.
            //         Fit Gumbel over [2, 19] to get mu, beta.
            //
            // Pass 2: refined plateau from max(15, mu + 4*beta),
            //         re-fix plateau, re-fit Gumbel over [2, 19].
            // ────────────────────────────────────────────────────────────
            double plateauVal = 1.0;
            TF1* fFit = nullptr;
            const double target90 = 0.90;
            double x90 = -1.0;
            if (perRunLabel.empty())
            {
                plateauVal = fitPlateauConst(hRatio, 15.0);
                fFit = buildGumbelFit(P, hRatio, plateauVal);
                
                if (fFit)
                {
                    // Pass 1
                    TFitResultPtr fitRes = hRatio->Fit(fFit, "RQS0", "", 2.0, 19.0);
                    int fitStatus = fitRes;
                    
                    if (fitStatus == 0)
                    {
                        const double mu   = fFit->GetParameter(1);
                        const double beta = std::fabs(fFit->GetParameter(2));
                        
                        // Pass 2: refined plateau from deep into the flat region
                        // Gumbel reaches ~98% of plateau at mu + 4*beta
                        const double refinedTailLo = std::max(15.0, std::min(19.0, mu + 4.0 * beta));
                        if (generatePerRunTriggerAnaForGoodRunsOnly)
                        {
                            plateauVal = 1.0;
                        }
                        else
                        {
                            plateauVal = fitPlateauConst(hRatio, refinedTailLo);
                        }
                        
                        fFit->SetParameter(0, plateauVal);
                        {
                            fFit->FixParameter(0, plateauVal);
                        }
                        fFit->SetParameter(1, mu);
                        
                        fitRes = hRatio->Fit(fFit, "RQS0", "", 2.0, 19.0);
                        fitStatus = fitRes;
                    }
                    
                    // Read back the fitted plateau (may have shifted up from tail estimate)
                    if (fitStatus == 0)
                    {
                        plateauVal = fFit->GetParameter(0);
                    }
                    
                    // Extract 80% crossing from the converged fit
                    if (fitStatus == 0)
                    {
                        // Gumbel CDF analytic: x_f = mu - beta * ln(-ln(f))
                        // For f = 0.80: ln(-ln(0.80)) = ln(0.22314) = -1.4999
                        // So x80 = mu + 1.500 * beta
                        const double mu   = fFit->GetParameter(1);
                        const double beta = std::fabs(fFit->GetParameter(2));
                        // Gumbel CDF analytic: x_f = mu - beta * ln(-ln(f))
                        // For f = 0.90: ln(-ln(0.90)) = ln(0.10536) = -2.2504
                        // So x90 = mu + 2.2504 * beta
                        const double x90Analytic = mu - beta * std::log(-std::log(0.90));
                        if (std::isfinite(x90Analytic) && x90Analytic >= 0.0 && x90Analytic <= 20.0)
                        {
                            x90 = x90Analytic;
                        }
                        else
                        {
                            // Numerical fallback using the fitted plateau
                            const double fitPlat = fFit->GetParameter(0);
                            const double x90Num = fFit->GetX(target90 * fitPlat, 0.0, 20.0);
                            if (std::isfinite(x90Num) && x90Num >= 0.0 && x90Num <= 20.0)
                                x90 = x90Num;
                        }
                    }
                }
                
                if (x90 < 0.0) x90 = findCrossingX(hRatio, target90, 2.0, 20.0);
            } // end if (perRunLabel.empty())
                
            ratioHists.push_back(hRatio);
            fitFuncs.push_back(fFit);
            x80s.push_back(x90);
            plateauVals.push_back(plateauVal);
            
            delete hBase;
            delete hProbe;
        }
        
        if (ratioHists.empty()) return;
        
        // ── draw the overlay canvas ──
        TCanvas c(TString::Format("c_trigAna_groupTurnOn_%s", loaded.front().groupFolder.c_str()).Data(),
                  TString::Format("c_trigAna_groupTurnOn_%s", loaded.front().groupFolder.c_str()).Data(),
                  900, 700);
        c.cd();
        c.SetLeftMargin(0.14);
        c.SetRightMargin(0.05);
        c.SetBottomMargin(0.14);
        c.SetTopMargin(0.08);
        c.SetTicks(1,1);
        
        ratioHists[0]->SetTitle("");
        ratioHists[0]->GetXaxis()->SetTitle("Maximum Cluster Energy [GeV]");
        ratioHists[0]->GetYaxis()->SetTitle("Photon X / MBD NS");
        ratioHists[0]->GetXaxis()->SetTitleSize(0.055);
        ratioHists[0]->GetXaxis()->SetTitleOffset(1.05);
        ratioHists[0]->GetXaxis()->SetLabelSize(0.045);
        ratioHists[0]->GetYaxis()->SetTitleSize(0.055);
        ratioHists[0]->GetYaxis()->SetTitleOffset(1.20);
        ratioHists[0]->GetYaxis()->SetLabelSize(0.045);
        ratioHists[0]->GetXaxis()->SetRangeUser(0.0, 20.0);
        ratioHists[0]->SetMinimum(0.0);
        ratioHists[0]->SetMaximum(1.20);
        ratioHists[0]->Draw("E1");
        
        for (std::size_t i = 1; i < ratioHists.size(); ++i)
        {
            ratioHists[i]->Draw("E1 SAME");
        }
        
        for (std::size_t i = 0; i < fitFuncs.size(); ++i)
        {
            if (fitFuncs[i]) fitFuncs[i]->Draw("SAME");
        }
        
        TLine l1(0.0, 1.0, 20.0, 1.0);
        l1.SetLineStyle(1);
        l1.SetLineWidth(3);
        l1.SetLineColor(kBlack);
        l1.Draw("SAME");
        
        TLine l90(0.0, 0.90, 20.0, 0.90);
        l90.SetLineStyle(3);
        l90.SetLineWidth(2);
        l90.SetLineColor(kGray+2);
        l90.Draw("SAME");

        // Vertical dashed lines at each trigger's 90% efficiency point
        std::vector<TLine*> vertLines;
        for (std::size_t i = 0; i < loaded.size() && i < x80s.size(); ++i)
        {
            if (x80s[i] <= 0.0 || x80s[i] > 20.0) continue;
            TLine* lv = new TLine(x80s[i], 0.0, x80s[i], 0.90);
            lv->SetLineStyle(2);
            lv->SetLineWidth(2);
            lv->SetLineColor(colorForProbe(loaded[i].probeKey));
            lv->Draw("SAME");
            vertLines.push_back(lv);
        }

        TLegend extraLegend(0.55, 0.45, 0.90, 0.55);
        extraLegend.SetBorderSize(0);
        extraLegend.SetFillStyle(0);
        extraLegend.SetTextFont(42);
        extraLegend.SetTextSize(0.040);
        extraLegend.AddEntry((TObject*)nullptr, "#it{#bf{sPHENIX}} Internal", "");
        extraLegend.AddEntry((TObject*)nullptr, "Au+Au #sqrt{s_{NN}} = 200 GeV", "");
        extraLegend.Draw();
        
        TLegend leg(0.18, 0.81, 0.45, 0.89);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.03);
        for (std::size_t i = 0; i < loaded.size() && i < ratioHists.size(); ++i)
        {
            if (x80s[i] > 0.0)
            {
                leg.AddEntry(ratioHists[i],
                             TString::Format("%s (90%%=%.2f GeV)",
                                             loaded[i].probeLabel.c_str(),
                                             x80s[i]).Data(),
                             "pe");
            }
            else
            {
                leg.AddEntry(ratioHists[i],
                             loaded[i].probeLabel.c_str(),
                             "pe");
            }
        }
        leg.Draw();
        
        // ── fit info: bottom-right ──
        if (perRunLabel.empty())
        {
            TLatex tFit;
            tFit.SetNDC(true);
            tFit.SetTextFont(43);
            tFit.SetTextAlign(33);

            double yText = 0.30;

            // Fit model label (larger)
            tFit.SetTextSize(24);
            tFit.SetTextColor(kBlack);
            tFit.DrawLatex(0.9, yText, "Fit: Gumbel CDF");
            yText -= 0.045;

            // Formula (larger)
            tFit.SetTextSize(22);
            tFit.DrawLatex(0.9, yText, "f(x) = P #times exp[#minus exp(#minus(x#minus#mu)/#beta)]");
            yText -= 0.045;

            // 90% efficiency line label (larger)
            tFit.SetTextSize(24);
            tFit.SetTextColor(kGray+2);
            tFit.DrawLatex(0.9, yText, "- - - 90% efficiency");
        }

        // ── plateau annotations: LHS, just under the 0.90 line ──
        if (perRunLabel.empty())
        {
            TLatex tPlat;
            tPlat.SetNDC(true);
            tPlat.SetTextFont(43);
            tPlat.SetTextAlign(13);
            tPlat.SetTextSize(22);

            double yText = 0.72;
            for (std::size_t i = 0; i < loaded.size() && i < plateauVals.size(); ++i)
            {
                const int thr = extractPhotonThresholdGeV(loaded[i].probeKey);
                tPlat.SetTextColor(colorForProbe(loaded[i].probeKey));
                tPlat.DrawLatex(0.16, yText,
                                TString::Format("Photon %d plateau = %.3f", thr, plateauVals[i]).Data());
                yText -= 0.04;
            }
            tPlat.SetTextColor(kBlack);
        }
        
        {
            TLatex tDS;
            tDS.SetNDC(true);
            tDS.SetTextFont(42);
            tDS.SetTextAlign(23);
            tDS.SetTextSize(0.035);
            if (isAuAuOnly || isSimAndDataAUAU)
            {
                const std::string& gf = loaded.front().groupFolder;
                if (gf.find("vtx_lt_10") != std::string::npos)
                    tDS.DrawLatex(0.50, 0.97, "Run3auau Photon 10 and 12 GeV + MBD NS #geq 2, vtx < 10 cm Efficiencies");
                else
                    tDS.DrawLatex(0.50, 0.97, "Run3auau Photon 10 and 12 GeV + MBD NS #geq 2, vtx < 150 cm Efficiencies");
            }
            else
                tDS.DrawLatex(0.50, 0.97, isRun25pp ? "Run25pp" : "Run24pp");
        }
        
        if (!perRunLabel.empty())
        {
            const std::string rlShort =
                (perRunLabel.rfind("000", 0) == 0 && perRunLabel.size() > 3)
                ? perRunLabel.substr(3) : perRunLabel;
            TLatex tRun;
            tRun.SetNDC(true);
            tRun.SetTextFont(42);
            tRun.SetTextAlign(22);
            tRun.SetTextSize(0.10);
            tRun.SetTextColor(kBlack);
            tRun.DrawLatex(0.72, 0.38, rlShort.c_str());
        }
        
        const std::string outPng = JoinPath(groupOutDir, "hMaxClusterEnergy_groupTurnOnOverlay" + filenameTag + ".png");
        SaveCanvas(c, outPng);
        cout << ANSI_BOLD_GRN << "[WROTE] " << outPng << ANSI_RESET << "\n";
        
        for (auto* lv : vertLines) delete lv;
        for (auto* f : fitFuncs) if (f) delete f;
        for (auto* h : ratioHists) delete h;
    };
    auto drawGroupMaxClusterOverlay = [&](const std::vector<LoadedPair>& loaded, const std::string& groupOutDir, const std::string& filenameTag = "", const std::string& perRunLabel = "")
    {
        if (loaded.empty()) return;
        
        std::vector<TH1*> probeHists;
        TH1* hBaseGroup = nullptr;
        
        for (const auto& P : loaded)
        {
            TH1* hB = CloneTH1(P.hBase,  TString::Format("hBase_groupOverlay_%s",  P.probeKey.c_str()).Data());
            TH1* hP = CloneTH1(P.hProbe, TString::Format("hProbe_groupOverlay_%s", P.probeKey.c_str()).Data());
            if (!hB || !hP)
            {
                if (hB) delete hB;
                if (hP) delete hP;
                continue;
            }
            if (!hBaseGroup) { hBaseGroup = hB; }
            else             { delete hB; }
            hP->SetLineColor(colorForProbe(P.probeKey));
            hP->SetLineWidth(4);
            probeHists.push_back(hP);
        }
        
        if (!hBaseGroup || probeHists.empty())
        {
            if (hBaseGroup) delete hBaseGroup;
            for (auto* h : probeHists) delete h;
            return;
        }
        
        double ymax = hBaseGroup->GetMaximum();
        for (auto* h : probeHists) ymax = std::max(ymax, h->GetMaximum());
        
        double minPos = SmallestPositiveBinContent(hBaseGroup);
        for (auto* h : probeHists)
        {
            const double m = SmallestPositiveBinContent(h);
            if (m > 0.0) { if (minPos > 0.0) minPos = std::min(minPos, m); else minPos = m; }
        }
        
        TCanvas c(TString::Format("c_trigAna_groupOverlay_%s", loaded.front().groupFolder.c_str()).Data(),
                  TString::Format("c_trigAna_groupOverlay_%s", loaded.front().groupFolder.c_str()).Data(),
                  900, 700);
        c.cd();
        c.SetLeftMargin(0.14);
        c.SetRightMargin(0.05);
        c.SetBottomMargin(0.14);
        c.SetTopMargin(0.08);
        c.SetTicks(1,1);
        c.SetLogy();
        
        hBaseGroup->SetTitle("");
        hBaseGroup->GetXaxis()->SetTitle("Maximum Cluster Energy [GeV]");
        hBaseGroup->GetYaxis()->SetTitle("Counts");
        hBaseGroup->GetXaxis()->SetTitleSize(0.055);
        hBaseGroup->GetXaxis()->SetTitleOffset(1.05);
        hBaseGroup->GetXaxis()->SetLabelSize(0.045);
        hBaseGroup->GetYaxis()->SetTitleSize(0.055);
        hBaseGroup->GetYaxis()->SetTitleOffset(1.20);
        hBaseGroup->GetYaxis()->SetLabelSize(0.045);
        hBaseGroup->GetXaxis()->SetRangeUser(0.0, 20.0);
        hBaseGroup->SetMinimum((minPos > 0.0) ? (0.5 * minPos) : 1e-6);
        hBaseGroup->SetMaximum((ymax  > 0.0) ? (1.25 * ymax)  : 1.0);
        hBaseGroup->SetLineColor(kBlack);
        hBaseGroup->SetLineWidth(4);
        hBaseGroup->Draw("HIST");
        
        for (auto* h : probeHists) h->Draw("HIST SAME");
        
        TLegend leg(0.4, 0.68, 0.95, 0.86);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.028);
        {
            // Use basisLabel (without "Baseline raw" prefix) for the black curve
            std::string baseLbl = loaded.front().basisLabel;
            auto replGeq = [](std::string s) -> std::string {
                std::size_t p = 0;
                while ((p = s.find(">=", p)) != std::string::npos) { s.replace(p, 2, "#geq"); p += 4; }
                return s;
            };
            leg.AddEntry(hBaseGroup, replGeq(baseLbl).c_str(), "l");
            for (std::size_t i = 0; i < loaded.size() && i < probeHists.size(); ++i)
                leg.AddEntry(probeHists[i], replGeq(loaded[i].probeLabel).c_str(), "l");
        }
        leg.Draw();

        TLegend sphLabel(0.18, 0.18, 0.50, 0.28);
        sphLabel.SetBorderSize(0);
        sphLabel.SetFillStyle(0);
        sphLabel.SetTextFont(42);
        sphLabel.SetTextSize(0.040);
        sphLabel.AddEntry((TObject*)nullptr, "#it{#bf{sPHENIX}} Internal", "");
        sphLabel.AddEntry((TObject*)nullptr, "Au+Au #sqrt{s_{NN}} = 200 GeV", "");
        sphLabel.Draw();

        {
            TLatex tDS;
            tDS.SetNDC(true);
            tDS.SetTextFont(42);
            tDS.SetTextAlign(23);
            tDS.SetTextSize(0.032);
            if (isAuAuOnly || isSimAndDataAUAU)
            {
                const std::string& gf = loaded.front().groupFolder;
                if (gf.find("vtx_lt_10") != std::string::npos)
                    tDS.DrawLatex(0.50, 0.97, "Run3auau MBD NS #geq 2, vtx < 10 cm, and Photon 10/12 Max Cluster Energy Overlay");
                else
                    tDS.DrawLatex(0.50, 0.97, "Run3auau MBD NS #geq 2, vtx < 150 cm, and Photon 10/12 Max Cluster Energy Overlay");
            }
            else
                tDS.DrawLatex(0.50, 0.97, isRun25pp ? "Run25pp" : "Run24pp");
        }
        
        if (!perRunLabel.empty())
        {
            const std::string rlShort =
                (perRunLabel.rfind("000", 0) == 0 && perRunLabel.size() > 3)
                ? perRunLabel.substr(3) : perRunLabel;
            TLatex tRun;
            tRun.SetNDC(true);
            tRun.SetTextFont(42);
            tRun.SetTextAlign(22);
            tRun.SetTextSize(0.10);
            tRun.SetTextColor(kBlack);
            tRun.DrawLatex(0.7, 0.50, rlShort.c_str());
        }
        
        const std::string outPng = JoinPath(groupOutDir, "hMaxClusterEnergy_groupOverlay" + filenameTag + ".png");
        SaveCanvas(c, outPng);
        cout << ANSI_BOLD_GRN << "[WROTE] " << outPng << ANSI_RESET << "\n";
        
        delete hBaseGroup;
        for (auto* h : probeHists) delete h;
    };
    
    struct PerRunTriggerAnaItem
    {
        std::string runLabel;
        std::string filePath;
    };
    
    auto filterPhoton10And12FromGroup = [&](const std::vector<LoadedPair>& input, const std::string& targetGroup) -> std::vector<LoadedPair>
    {
        std::vector<LoadedPair> out;
        for (const auto& P : input)
        {
            if (P.groupFolder != targetGroup) continue;
            
            const int thr = extractPhotonThresholdGeV(P.probeKey);
            if (thr == 10 || thr == 12) out.push_back(P);
        }
        
        std::sort(out.begin(), out.end(), [&](const LoadedPair& a, const LoadedPair& b)
                  {
            const int ta = extractPhotonThresholdGeV(a.probeKey);
            const int tb = extractPhotonThresholdGeV(b.probeKey);
            if (ta != tb) return ta < tb;
            return a.probeKey < b.probeKey;
        });
        
        return out;
    };
    
    auto drawCompactGroupMaxClusterOverlayPad = [&](const std::vector<LoadedPair>& loaded,
                                                    const std::string& runLabel,
                                                    std::vector<TObject*>& keepAlive)
    {
        if (!gPad) return;
        
        const std::string runLabelShort =
            (runLabel.rfind("000", 0) == 0 && runLabel.size() > 3) ? runLabel.substr(3) : runLabel;
        
        gPad->SetLeftMargin(0.14);
        gPad->SetRightMargin(0.03);
        gPad->SetBottomMargin(0.14);
        gPad->SetTopMargin(0.10);
        gPad->SetTicks(1,1);
        gPad->SetLogy();
        
        if (loaded.empty())
        {
            TLatex t;
            t.SetNDC(true);
            t.SetTextFont(42);
            t.SetTextAlign(33);
            t.SetTextSize(0.120);
            t.DrawLatex(0.88, 0.82, runLabelShort.c_str()); 
            t.SetTextAlign(22);
            t.SetTextSize(0.120);
            t.DrawLatex(0.50, 0.50, "MISSING");
            return;
        }
        
        std::vector<TH1*> probeHists;
        TH1* hBaseGroup = nullptr;
        
        for (const auto& P : loaded)
        {
            TH1* hB = CloneTH1(P.hBase,  TString::Format("hBase_compactOverlay_%s_%s",  P.probeKey.c_str(), runLabel.c_str()).Data());
            TH1* hP = CloneTH1(P.hProbe, TString::Format("hProbe_compactOverlay_%s_%s", P.probeKey.c_str(), runLabel.c_str()).Data());
            if (!hB || !hP)
            {
                if (hB) delete hB;
                if (hP) delete hP;
                continue;
            }
            
            if (!hBaseGroup)
            {
                hBaseGroup = hB;
                keepAlive.push_back(hBaseGroup);
            }
            else
            {
                delete hB;
            }
            
            hP->SetLineColor(colorForProbe(P.probeKey));
            hP->SetLineWidth(2);
            probeHists.push_back(hP);
            keepAlive.push_back(hP);
        }
        
        if (!hBaseGroup || probeHists.empty())
        {
            TLatex t;
            t.SetNDC(true);
            t.SetTextFont(42);
            t.SetTextAlign(33);
            t.SetTextSize(0.120);
            t.DrawLatex(0.88, 0.82, runLabelShort.c_str()); 
            t.SetTextAlign(22);
            t.SetTextSize(0.120);
            t.DrawLatex(0.50, 0.50, "MISSING");
            return;
        }
        
        double ymax = hBaseGroup->GetMaximum();
        for (auto* h : probeHists) ymax = std::max(ymax, h->GetMaximum());
        
        double minPos = SmallestPositiveBinContent(hBaseGroup);
        for (auto* h : probeHists)
        {
            const double m = SmallestPositiveBinContent(h);
            if (m > 0.0)
            {
                if (minPos > 0.0) minPos = std::min(minPos, m);
                else              minPos = m;
            }
        }
        
        hBaseGroup->SetTitle("");
        hBaseGroup->GetXaxis()->SetTitle("Max E [GeV]");
        hBaseGroup->GetYaxis()->SetTitle("Counts");
        hBaseGroup->GetXaxis()->SetTitleSize(0.075);
        hBaseGroup->GetYaxis()->SetTitleSize(0.075);
        hBaseGroup->GetXaxis()->SetLabelSize(0.060);
        hBaseGroup->GetYaxis()->SetLabelSize(0.060);
        hBaseGroup->GetYaxis()->SetTitleOffset(0.95);
        hBaseGroup->GetXaxis()->SetRangeUser(0.0, 20.0);
        hBaseGroup->SetMinimum((minPos > 0.0) ? (0.5 * minPos) : 1e-6);
        hBaseGroup->SetMaximum((ymax > 0.0) ? (1.25 * ymax) : 1.0);
        hBaseGroup->SetLineColor(kBlack);
        hBaseGroup->SetLineWidth(2);
        hBaseGroup->Draw("HIST");
        
        for (auto* h : probeHists) h->Draw("HIST SAME");
        
        TLatex t;
        t.SetNDC(true);
        t.SetTextFont(42);
        t.SetTextAlign(33);
        t.SetTextSize(0.120);
        t.SetTextColor(kBlack);
        t.DrawLatex(0.95, 0.88, runLabelShort.c_str());
    };
    
    auto drawCompactGroupTurnOnOverlayPad = [&](const std::vector<LoadedPair>& loaded,
                                                const std::string& runLabel,
                                                std::vector<TObject*>& keepAlive)
    {
        if (!gPad) return;
        
        const std::string runLabelShort =
            (runLabel.rfind("000", 0) == 0 && runLabel.size() > 3) ? runLabel.substr(3) : runLabel;
        
        gPad->SetLeftMargin(0.14);
        gPad->SetRightMargin(0.03);
        gPad->SetBottomMargin(0.14);
        gPad->SetTopMargin(0.10);
        gPad->SetTicks(1,1);
        
        if (loaded.empty())
        {
            TLatex t;
            t.SetNDC(true);
            t.SetTextFont(42);
            t.SetTextAlign(33);
            t.SetTextSize(0.120);
            t.DrawLatex(0.88, 0.82, runLabelShort.c_str()); 
            t.SetTextAlign(22);
            t.SetTextSize(0.120);
            t.DrawLatex(0.50, 0.50, "MISSING");
            return;
        }
        
        std::vector<TH1*> ratioHists;
        for (const auto& P : loaded)
        {
            TH1* hBase = CloneTH1(P.hBase,  TString::Format("hBase_compactTurnOn_%s_%s",  P.probeKey.c_str(), runLabel.c_str()).Data());
            TH1* hProbe = CloneTH1(P.hProbe, TString::Format("hProbe_compactTurnOn_%s_%s", P.probeKey.c_str(), runLabel.c_str()).Data());
            if (!hBase || !hProbe)
            {
                if (hBase) delete hBase;
                if (hProbe) delete hProbe;
                continue;
            }
            
            TH1* hRatio = CloneTH1(hProbe, TString::Format("hRatio_compactTurnOn_%s_%s", P.probeKey.c_str(), runLabel.c_str()).Data());
            if (!hRatio)
            {
                delete hBase;
                delete hProbe;
                continue;
            }
            
            EnsureSumw2(hBase);
            EnsureSumw2(hProbe);
            EnsureSumw2(hRatio);
            hRatio->SetDirectory(nullptr);
            hRatio->Divide(hProbe, hBase, 1.0, 1.0, "B");
            hRatio->SetMarkerStyle(markerForProbe(P.probeKey));
            hRatio->SetMarkerSize(0.70);
            hRatio->SetMarkerColor(colorForProbe(P.probeKey));
            hRatio->SetLineColor(colorForProbe(P.probeKey));
            hRatio->SetLineWidth(2);
            
            keepAlive.push_back(hBase);
            keepAlive.push_back(hProbe);
            keepAlive.push_back(hRatio);
            ratioHists.push_back(hRatio);
        }
        
        if (ratioHists.empty())
        {
            TLatex t;
            t.SetNDC(true);
            t.SetTextFont(42);
            t.SetTextAlign(33);
            t.SetTextSize(0.120);
            t.DrawLatex(0.88, 0.82, runLabelShort.c_str()); 
            t.SetTextAlign(22);
            t.SetTextSize(0.120);
            t.DrawLatex(0.50, 0.50, "MISSING");
            return;
        }
        
        ratioHists[0]->SetTitle("");
        ratioHists[0]->GetXaxis()->SetTitle("Max E [GeV]");
        ratioHists[0]->GetYaxis()->SetTitle("Eff.");
        ratioHists[0]->GetXaxis()->SetTitleSize(0.075);
        ratioHists[0]->GetYaxis()->SetTitleSize(0.075);
        ratioHists[0]->GetXaxis()->SetLabelSize(0.060);
        ratioHists[0]->GetYaxis()->SetLabelSize(0.060);
        ratioHists[0]->GetYaxis()->SetTitleOffset(0.95);
        ratioHists[0]->GetXaxis()->SetRangeUser(0.0, 20.0);
        ratioHists[0]->SetMinimum(0.0);
        ratioHists[0]->SetMaximum(1.15);
        ratioHists[0]->Draw("E1");
        
        for (std::size_t i = 1; i < ratioHists.size(); ++i)
        {
            ratioHists[i]->Draw("E1 SAME");
        }
        
        TLine l1(0.0, 1.0, 20.0, 1.0);
        l1.SetLineStyle(1);
        l1.SetLineWidth(2);
        l1.SetLineColor(kBlack);
        l1.DrawClone("SAME");
        
        TLine l90(0.0, 0.90, 20.0, 0.90);
        l90.SetLineStyle(3);
        l90.SetLineWidth(2);
        l90.SetLineColor(kGray+2);
        l90.DrawClone("SAME");
        
        TLatex t;
        t.SetNDC(true);
        t.SetTextFont(42);
        t.SetTextAlign(33);
        t.SetTextSize(0.120);
        t.SetTextColor(kBlack);
        t.DrawLatex(0.95, 0.88, runLabelShort.c_str());
    };
    
    auto writePerRunTablePages = [&](const std::vector<PerRunTriggerAnaItem>& items,
                                     const std::string& targetGroup,
                                     const std::string& tableOutDir,
                                     bool makeTurnOn,
                                     int nCols, int nRows,
                                     int canvasW, int canvasH)
    {
        if (items.empty()) return;
        
        EnsureDir(tableOutDir);
        const std::size_t nPerPage = static_cast<std::size_t>(nCols * nRows);
        
        for (std::size_t pageStart = 0, pageIdx = 0; pageStart < items.size(); pageStart += nPerPage, ++pageIdx)
        {
            TCanvas c(TString::Format("c_perRunTrigAna_%s_%zu",
                                      makeTurnOn ? "turnOn" : "maxCluster", pageIdx + 1).Data(),
                      "c_perRunTrigAna", canvasW, canvasH);
            c.Divide(nCols, nRows, 0.001, 0.001);
            
            std::vector<TObject*> keepAlive;
            keepAlive.reserve(nPerPage * 8);
            
            for (std::size_t i = 0; i < nPerPage; ++i)
            {
                c.cd(static_cast<int>(i + 1));
                
                const std::size_t idx = pageStart + i;
                if (idx >= items.size())
                {
                    gPad->SetTicks(1,1);
                    TLatex t;
                    t.SetNDC(true);
                    t.SetTextFont(42);
                    t.SetTextAlign(22);
                    t.SetTextSize(0.120);
                    t.DrawLatex(0.50, 0.50, "EMPTY");
                    continue;
                }
                
                const auto& item = items[idx];
                const std::string runLabelShort =
                    (item.runLabel.rfind("000", 0) == 0 && item.runLabel.size() > 3) ? item.runLabel.substr(3) : item.runLabel;
                
                TFile* fRun = TFile::Open(item.filePath.c_str(), "READ");
                if (!fRun || fRun->IsZombie())
                {
                    if (fRun) { fRun->Close(); delete fRun; }
                    gPad->SetTicks(1,1);
                    TLatex t;
                    t.SetNDC(true);
                    t.SetTextFont(42);
                    t.SetTextAlign(33);
                    t.SetTextSize(0.120);
                    t.DrawLatex(0.88, 0.82, runLabelShort.c_str()); 
                    t.SetTextAlign(22);
                    t.SetTextSize(0.110);
                    t.DrawLatex(0.50, 0.50, "BAD FILE");
                    continue;
                }
                
                std::vector<LoadedPair> runPairs = collectPairsFromFile(fRun, item.filePath, false);
                std::vector<LoadedPair> filtered = filterPhoton10And12FromGroup(runPairs, targetGroup);
                
                if (makeTurnOn) drawCompactGroupTurnOnOverlayPad(filtered, item.runLabel, keepAlive);
                else            drawCompactGroupMaxClusterOverlayPad(filtered, item.runLabel, keepAlive);
                
                fRun->Close();
                delete fRun;
            }
            
            const std::string outPng = JoinPath(
                tableOutDir,
                TString::Format("%s_page%03zu.png",
                                makeTurnOn ? "turnOnRunByRunTable" : "maxClusterOverlaysRunByRunTable",
                                pageIdx + 1).Data());
            SaveCanvas(c, outPng);
            cout << ANSI_BOLD_GRN << "[WROTE] " << outPng << ANSI_RESET << "\n";
            
            for (auto* obj : keepAlive) delete obj;
        }
    };
    
    std::set<std::string> writtenGroups;
    for (const auto& folder : groupOrder)
    {
        auto it = groups.find(folder);
        if (it == groups.end()) continue;
        writtenGroups.insert(folder);
        
        const std::string groupOutDir = JoinPath(outDir, folder);
        EnsureDir(groupOutDir);
        
        std::vector<std::string> summary;
        summary.push_back(string("triggerQA doNotScale common-basis summary: ") + folder);
        summary.push_back(string("Input: ") + ds.inFilePath);
        summary.push_back(string("Basis: ") + it->second.front().basisLabel);
        summary.push_back("Schema: numerator=<probe>, denominator=baseline_<probe>");
        summary.push_back("Efficiency = probe live / baseline raw");
        summary.push_back("");
        
        for (const auto& P : it->second)
        {
            const std::string pairOutDir = JoinPath(groupOutDir, P.probeKey);
            EnsureDir(pairOutDir);
            
            double x95 = -1.0;
            drawPairOverlay(P, pairOutDir);
            drawPairTurnOn(P, pairOutDir, x95);
            
            summary.push_back(TString::Format(
                                              "%s | baseline=%s | entries(probe)=%.0f | entries(base)=%.0f | run spans (%d): %s | x95=%s",
                                              P.probeLabel.c_str(),
                                              P.baselineKey.c_str(),
                                              (P.hProbe ? P.hProbe->GetEntries() : 0.0),
                                              (P.hBase  ? P.hBase->GetEntries()  : 0.0),
                                              P.nRunSpans,
                                              P.coverageText.c_str(),
                                              (x95 > 0.0 ? TString::Format("%.4f GeV", x95).Data() : "n/a")
                                              ).Data());
        }
        
        if (isAuAuOnly && (folder == "commonBasis_MBD_NS_geq_2_vtx_lt_150" ||
                           folder == "commonBasis_MBD_NS_geq_2_vtx_lt_10"))
        {
            std::vector<LoadedPair> filteredAuAu;
            for (const auto& P : it->second)
            {
                const int thr = extractPhotonThresholdGeV(P.probeKey);
                if (thr == 10 || thr == 12) filteredAuAu.push_back(P);
            }
            if (!filteredAuAu.empty())
            {
                drawGroupMaxClusterOverlay(filteredAuAu, groupOutDir, "");
                drawGroupTurnOnOverlay(filteredAuAu, groupOutDir, "");
            }
            
            std::vector<LoadedPair> filtered6_8;
            for (const auto& P : it->second)
            {
                const int thr = extractPhotonThresholdGeV(P.probeKey);
                if (thr == 6 || thr == 8) filtered6_8.push_back(P);
            }
            if (!filtered6_8.empty())
            {
                drawGroupMaxClusterOverlay(filtered6_8, groupOutDir, "_photon6_8");
                drawGroupTurnOnOverlay(filtered6_8, groupOutDir, "_photon6_8");
            }
        }
        else
        {
            drawGroupTurnOnOverlay(it->second, groupOutDir, "");
        }

        WriteTextFile(JoinPath(groupOutDir, "summary.txt"), summary);
        indexLines.push_back(folder + " -> " + it->second.front().basisLabel);
        for (const auto& P : it->second)
        {
            indexLines.push_back(TString::Format(
                                                 "  %s | baseline_%s | spans=%d | %s",
                                                 P.probeKey.c_str(),
                                                 P.probeKey.c_str(),
                                                 P.nRunSpans,
                                                 P.coverageText.c_str()).Data());
        }
        indexLines.push_back("");
    }
    
    for (const auto& kv : groups)
    {
        if (writtenGroups.count(kv.first)) continue;
        
        const std::string groupOutDir = JoinPath(outDir, kv.first);
        EnsureDir(groupOutDir);
        
        std::vector<std::string> summary;
        summary.push_back(string("triggerQA doNotScale common-basis summary: ") + kv.first);
        summary.push_back(string("Input: ") + ds.inFilePath);
        summary.push_back(string("Basis: ") + kv.second.front().basisLabel);
        summary.push_back("Schema: numerator=<probe>, denominator=baseline_<probe>");
        summary.push_back("Efficiency = probe live / baseline raw");
        summary.push_back("");
        
        for (const auto& P : kv.second)
        {
            const std::string pairOutDir = JoinPath(groupOutDir, P.probeKey);
            EnsureDir(pairOutDir);
            
            double x95 = -1.0;
            drawPairOverlay(P, pairOutDir);
            drawPairTurnOn(P, pairOutDir, x95);
            
            summary.push_back(TString::Format(
                                              "%s | baseline=%s | entries(probe)=%.0f | entries(base)=%.0f | run spans (%d): %s | x95=%s",
                                              P.probeLabel.c_str(),
                                              P.baselineKey.c_str(),
                                              (P.hProbe ? P.hProbe->GetEntries() : 0.0),
                                              (P.hBase  ? P.hBase->GetEntries()  : 0.0),
                                              P.nRunSpans,
                                              P.coverageText.c_str(),
                                              (x95 > 0.0 ? TString::Format("%.4f GeV", x95).Data() : "n/a")
                                              ).Data());
        }
        
        if (isAuAuOnly && (kv.first == "commonBasis_MBD_NS_geq_2_vtx_lt_150" ||
                           kv.first == "commonBasis_MBD_NS_geq_2_vtx_lt_10"))
        {
            std::vector<LoadedPair> filteredAuAu;
            for (const auto& P : kv.second)
            {
                const int thr = extractPhotonThresholdGeV(P.probeKey);
                if (thr == 10 || thr == 12) filteredAuAu.push_back(P);
            }
            if (!filteredAuAu.empty())
            {
                drawGroupMaxClusterOverlay(filteredAuAu, groupOutDir, "");
                drawGroupTurnOnOverlay(filteredAuAu, groupOutDir, "");
            }
            
            std::vector<LoadedPair> filtered6_8;
            for (const auto& P : kv.second)
            {
                const int thr = extractPhotonThresholdGeV(P.probeKey);
                if (thr == 6 || thr == 8) filtered6_8.push_back(P);
            }
            if (!filtered6_8.empty())
            {
                drawGroupMaxClusterOverlay(filtered6_8, groupOutDir, "_photon6_8");
                drawGroupTurnOnOverlay(filtered6_8, groupOutDir, "_photon6_8");
            }
        }
        else
        {
            drawGroupTurnOnOverlay(kv.second, groupOutDir, "");
        }
        WriteTextFile(JoinPath(groupOutDir, "summary.txt"), summary);
        
        indexLines.push_back(kv.first + " -> " + kv.second.front().basisLabel);
        for (const auto& P : kv.second)
        {
            indexLines.push_back(TString::Format(
                                                 "  %s | baseline_%s | spans=%d | %s",
                                                 P.probeKey.c_str(),
                                                 P.probeKey.c_str(),
                                                 P.nRunSpans,
                                                 P.coverageText.c_str()).Data());
        }
        indexLines.push_back("");
    }
    
    WriteTextFile(JoinPath(outDir, "summary_triggerQA_doNotScale.txt"), indexLines);
    
    if (generatePerRunTriggerAnaForGoodRunsOnly && generatePerRunTriggerAna)
    {
        std::cerr << ANSI_BOLD_RED
                  << "[FATAL] generatePerRunTriggerAnaForGoodRunsOnly and generatePerRunTriggerAna "
                     "cannot both be true. Set generatePerRunTriggerAna = false when using good-runs-only mode."
                  << ANSI_RESET << "\n";
        return;
    }
    
    if (generatePerRunTriggerAnaForGoodRunsOnly && datasetBucket() == "auau")
    {
        if (goodTriggerRuns.empty())
        {
            cout << ANSI_BOLD_YEL
                 << "[WARN] generatePerRunTriggerAnaForGoodRunsOnly is true but goodTriggerRuns is empty. Skipping."
                 << ANSI_RESET << "\n";
        }
        else
        {
            const std::string targetGroup = "commonBasis_MBD_NS_geq_2_vtx_lt_150";
            const std::string perRunInputDir = JoinPath(kInputBase + "/auau25/perRun", CfgTagWithUE());
            const std::string groupOutDir = JoinPath(outDir, targetGroup);
            EnsureDir(groupOutDir);
            
            const std::set<int> goodRunSet(goodTriggerRuns.begin(), goodTriggerRuns.end());
            std::map<std::string, LoadedPair> accPairs;
            int nGoodFound = 0;
            int nGoodMissing = 0;
            std::vector<int> missingRuns;
            
            cout << ANSI_BOLD_GRN
                 << "[GOOD-RUNS-ONLY] Combining " << goodTriggerRuns.size()
                 << " runs for group overlays..." << ANSI_RESET << "\n";
            
            for (int runNum : goodTriggerRuns)
            {
                const std::string runDigits = TString::Format("%08d", runNum).Data();
                const std::string fname = "chunkMerge_run_" + runDigits + ".root";
                const std::string fpath = JoinPath(perRunInputDir, fname);
                
                TFile* fRun = TFile::Open(fpath.c_str(), "READ");
                if (!fRun || fRun->IsZombie())
                {
                    ++nGoodMissing;
                    missingRuns.push_back(runNum);
                    if (fRun) { fRun->Close(); delete fRun; }
                    continue;
                }
                
                std::vector<LoadedPair> runPairs = collectPairsFromFile(fRun, fpath, false);
                std::vector<LoadedPair> filtered = filterPhoton10And12FromGroup(runPairs, targetGroup);
                
                if (filtered.empty())
                {
                    ++nGoodMissing;
                    missingRuns.push_back(runNum);
                    fRun->Close();
                    delete fRun;
                    continue;
                }
                
                ++nGoodFound;
                for (const auto& P : filtered)
                {
                    auto ait = accPairs.find(P.probeKey);
                    if (ait == accPairs.end())
                    {
                        LoadedPair LP;
                        LP.basisKey      = P.basisKey;
                        LP.basisLabel    = P.basisLabel;
                        LP.groupFolder   = P.groupFolder;
                        LP.probeKey      = P.probeKey;
                        LP.probeLabel    = P.probeLabel;
                        LP.baselineKey   = P.baselineKey;
                        LP.baselineLabel = P.baselineLabel;
                        LP.coverageText  = "goodRunsOnly combined";
                        LP.nRunSpans     = 1;
                        LP.hProbe = static_cast<TH1*>(P.hProbe->Clone(
                            TString::Format("goodAcc_probe_%s", P.probeKey.c_str())));
                        LP.hProbe->SetDirectory(nullptr);
                        LP.hBase = static_cast<TH1*>(P.hBase->Clone(
                            TString::Format("goodAcc_base_%s", P.probeKey.c_str())));
                        LP.hBase->SetDirectory(nullptr);
                        accPairs[P.probeKey] = LP;
                    }
                    else
                    {
                        ait->second.hProbe->Add(P.hProbe);
                        ait->second.hBase->Add(P.hBase);
                        ait->second.nRunSpans += 1;
                    }
                }
                
                fRun->Close();
                delete fRun;
            }
            
            if (!accPairs.empty())
            {
                std::vector<LoadedPair> accVec;
                accVec.reserve(accPairs.size());
                for (auto& kv : accPairs) accVec.push_back(kv.second);
                std::sort(accVec.begin(), accVec.end(), [&](const LoadedPair& a, const LoadedPair& b)
                {
                    return extractPhotonThresholdGeV(a.probeKey) < extractPhotonThresholdGeV(b.probeKey);
                });
                
                drawGroupMaxClusterOverlay(accVec, groupOutDir, "");
                drawGroupTurnOnOverlay(accVec, groupOutDir, "");
                
                cout << ANSI_BOLD_GRN
                     << "[GOOD-RUNS-ONLY] Built combined overlays from " << nGoodFound
                     << " good runs (" << nGoodMissing << " missing/empty)"
                     << ANSI_RESET << "\n";
                
                if (!missingRuns.empty())
                {
                    cout << ANSI_BOLD_YEL << "  Missing/empty runs: ";
                    for (std::size_t i = 0; i < missingRuns.size(); ++i)
                    {
                        if (i) cout << ", ";
                        cout << missingRuns[i];
                    }
                    cout << ANSI_RESET << "\n";
                }
                
                for (auto& kv : accPairs)
                {
                    delete kv.second.hProbe;
                    delete kv.second.hBase;
                }
            }
            else
            {
                cout << ANSI_BOLD_RED
                     << "[GOOD-RUNS-ONLY] No usable data found from any of the "
                     << goodTriggerRuns.size() << " listed runs."
                     << ANSI_RESET << "\n";
            }
        }
    }
    
    if (generatePerRunTriggerAnaForGoodRunsOnly && datasetBucket() == "auau")
    {
        if (!goodTriggerRuns_vtx10.empty())
        {
            const std::string targetGroup = "commonBasis_MBD_NS_geq_2_vtx_lt_10";
            const std::string perRunInputDir = JoinPath(kInputBase + "/auau25/perRun", CfgTagWithUE());
            const std::string groupOutDir = JoinPath(outDir, targetGroup);
            EnsureDir(groupOutDir);
            
            std::map<std::string, LoadedPair> accPairs;
            int nGoodFound = 0;
            int nGoodMissing = 0;
            std::vector<int> missingRuns;
            
            cout << ANSI_BOLD_GRN
                 << "[GOOD-RUNS-ONLY vtx10] Combining " << goodTriggerRuns_vtx10.size()
                 << " runs for group overlays..." << ANSI_RESET << "\n";
            
            for (int runNum : goodTriggerRuns_vtx10)
            {
                const std::string runDigits = TString::Format("%08d", runNum).Data();
                const std::string fname = "chunkMerge_run_" + runDigits + ".root";
                const std::string fpath = JoinPath(perRunInputDir, fname);
                
                TFile* fRun = TFile::Open(fpath.c_str(), "READ");
                if (!fRun || fRun->IsZombie())
                {
                    ++nGoodMissing;
                    missingRuns.push_back(runNum);
                    if (fRun) { fRun->Close(); delete fRun; }
                    continue;
                }
                
                std::vector<LoadedPair> runPairs = collectPairsFromFile(fRun, fpath, false);
                std::vector<LoadedPair> filtered = filterPhoton10And12FromGroup(runPairs, targetGroup);
                
                if (filtered.empty())
                {
                    ++nGoodMissing;
                    missingRuns.push_back(runNum);
                    fRun->Close();
                    delete fRun;
                    continue;
                }
                
                ++nGoodFound;
                for (const auto& P : filtered)
                {
                    auto ait = accPairs.find(P.probeKey);
                    if (ait == accPairs.end())
                    {
                        LoadedPair LP;
                        LP.basisKey      = P.basisKey;
                        LP.basisLabel    = P.basisLabel;
                        LP.groupFolder   = P.groupFolder;
                        LP.probeKey      = P.probeKey;
                        LP.probeLabel    = P.probeLabel;
                        LP.baselineKey   = P.baselineKey;
                        LP.baselineLabel = P.baselineLabel;
                        LP.coverageText  = "goodRunsOnly vtx10 combined";
                        LP.nRunSpans     = 1;
                        LP.hProbe = static_cast<TH1*>(P.hProbe->Clone(
                            TString::Format("goodAcc10_probe_%s", P.probeKey.c_str())));
                        LP.hProbe->SetDirectory(nullptr);
                        LP.hBase = static_cast<TH1*>(P.hBase->Clone(
                            TString::Format("goodAcc10_base_%s", P.probeKey.c_str())));
                        LP.hBase->SetDirectory(nullptr);
                        accPairs[P.probeKey] = LP;
                    }
                    else
                    {
                        ait->second.hProbe->Add(P.hProbe);
                        ait->second.hBase->Add(P.hBase);
                        ait->second.nRunSpans += 1;
                    }
                }
                
                fRun->Close();
                delete fRun;
            }
            
            if (!accPairs.empty())
            {
                std::vector<LoadedPair> accVec;
                accVec.reserve(accPairs.size());
                for (auto& kv : accPairs) accVec.push_back(kv.second);
                std::sort(accVec.begin(), accVec.end(), [&](const LoadedPair& a, const LoadedPair& b)
                {
                    return extractPhotonThresholdGeV(a.probeKey) < extractPhotonThresholdGeV(b.probeKey);
                });
                
                drawGroupMaxClusterOverlay(accVec, groupOutDir, "");
                drawGroupTurnOnOverlay(accVec, groupOutDir, "");
                
                cout << ANSI_BOLD_GRN
                     << "[GOOD-RUNS-ONLY vtx10] Built combined overlays from " << nGoodFound
                     << " good runs (" << nGoodMissing << " missing/empty)"
                     << ANSI_RESET << "\n";
                
                for (auto& kv : accPairs)
                {
                    delete kv.second.hProbe;
                    delete kv.second.hBase;
                }
            }
        }
    }
    
    if (generatePerRunTriggerAna && datasetBucket() == "auau")
    {
        const std::string targetGroup = "commonBasis_MBD_NS_geq_2_vtx_lt_150";
        const std::string perRunInputDir = JoinPath(kInputBase + "/auau25/perRun", CfgTagWithUE());
        const std::string perRunGroupOutDir = JoinPath(JoinPath(outDir, targetGroup), "perRun");
        EnsureDir(perRunGroupOutDir);
        
        const bool perRunDirExists = !gSystem->AccessPathName(perRunInputDir.c_str());
        
        std::vector<PerRunTriggerAnaItem> perRunItems;
        std::vector<std::string> perRunSummary;
        perRunSummary.push_back("triggerQA doNotScale per-run summary");
        perRunSummary.push_back(string("Combined input: ") + ds.inFilePath);
        perRunSummary.push_back(string("Per-run input dir: ") + perRunInputDir);
        perRunSummary.push_back(string("Per-run directory exists: ") + (perRunDirExists ? "YES" : "NO"));
        perRunSummary.push_back(string("Target group: ") + targetGroup);
        perRunSummary.push_back("Standalone outputs: hMaxClusterEnergy_groupOverlay.png and hMaxClusterEnergy_groupTurnOnOverlay.png");
        perRunSummary.push_back("");
        
        int nDiscoveredRuns = 0;
        int nOpenableRuns = 0;
        int nBadFiles = 0;
        int nRunsWithAnyDoNotScalePairs = 0;
        int nRunsWithoutAnyDoNotScalePairs = 0;
        int nRunsWithTargetGroup = 0;
        int nRunsWithoutTargetGroup = 0;
        int nRunsWithPhoton10 = 0;
        int nRunsWithPhoton12 = 0;
        int nRunsWithBothPhoton10And12 = 0;
        int nRunsRenderable = 0;
        int nRunsTargetGroupButNoPhoton10Or12 = 0;
        
        const std::set<int> badRunSet(badTriggerRuns.begin(), badTriggerRuns.end());
        std::map<std::string, LoadedPair> accPairs;
        int nBadTriggerRunsSkipped = 0;
        std::map<int, std::vector<std::string>> scaledMissingByPattern;
        
        void* dirp = gSystem->OpenDirectory(perRunInputDir.c_str());
        if (!dirp)
        {
            perRunSummary.push_back("BREAKPOINT: failed before run discovery.");
            perRunSummary.push_back(string("OpenDirectory failed for: ") + perRunInputDir);
            perRunSummary.push_back(string("AccessPathName says exists: ") + (perRunDirExists ? "YES" : "NO"));
            
            cout << ANSI_BOLD_YEL
                 << "[WARN] Could not open per-run triggerQA input directory: "
                 << perRunInputDir
                 << ANSI_RESET << "\n";
        }
        else
        {
            std::vector<PerRunTriggerAnaItem> discovered;
            const char* entry = nullptr;
            while ((entry = gSystem->GetDirEntry(dirp)))
            {
                const std::string fname = entry;
                if (fname == "." || fname == "..") continue;
                if (fname.rfind("chunkMerge_run_", 0) != 0) continue;
                if (fname.size() <= std::string("chunkMerge_run_").size() + std::string(".root").size()) continue;
                if (fname.substr(fname.size() - 5) != ".root") continue;
                
                const std::string runDigits = fname.substr(std::string("chunkMerge_run_").size(),
                                                           fname.size() - std::string("chunkMerge_run_").size() - 5);
                if (runDigits.empty()) continue;
                
                PerRunTriggerAnaItem item;
                item.runLabel = runDigits;
                item.filePath = JoinPath(perRunInputDir, fname);
                discovered.push_back(item);
            }
            gSystem->FreeDirectory(dirp);
            
            std::sort(discovered.begin(), discovered.end(),
                      [&](const PerRunTriggerAnaItem& a, const PerRunTriggerAnaItem& b)
                      {
                          return a.runLabel < b.runLabel;
                      });
            
            nDiscoveredRuns = static_cast<int>(discovered.size());
            perRunSummary.push_back(TString::Format("Discovered chunkMerge_run files: %d", nDiscoveredRuns).Data());
            
            if (discovered.empty())
            {
                perRunSummary.push_back("BREAKPOINT: directory opened, but zero chunkMerge_run_*.root files were found.");
            }
            
            for (const auto& item : discovered)
            {
                TFile* fRun = TFile::Open(item.filePath.c_str(), "READ");
                if (!fRun || fRun->IsZombie())
                {
                    ++nBadFiles;
                    perRunSummary.push_back(string("RUN ") + item.runLabel + " : BAD_FILE | path=" + item.filePath);
                    if (fRun) { fRun->Close(); delete fRun; }
                    continue;
                }
                
                ++nOpenableRuns;
                
                {
                    static const std::string kScaledKeys[3] = {
                        "MBD_NS_geq_2_vtx_lt_150",
                        "photon_10_plus_MBD_NS_geq_2_vtx_lt_150",
                        "photon_12_plus_MBD_NS_geq_2_vtx_lt_150"
                    };
                    int missingFlags = 0;
                    for (int b = 0; b < 3; ++b)
                    {
                        const std::string cntPath = kScaledKeys[b] + "/cnt_" + kScaledKeys[b];
                        TH1* hCnt = dynamic_cast<TH1*>(fRun->Get(cntPath.c_str()));
                        if (!hCnt || hCnt->GetEntries() <= 0)
                            missingFlags |= (1 << b);
                    }
                    if (missingFlags > 0)
                    {
                        const std::string shortLbl =
                            (item.runLabel.rfind("000", 0) == 0 && item.runLabel.size() > 3)
                            ? item.runLabel.substr(3) : item.runLabel;
                        scaledMissingByPattern[missingFlags].push_back(shortLbl);
                    }
                }
                
                std::vector<LoadedPair> runPairs = collectPairsFromFile(fRun, item.filePath, false);
                const int nAllPairs = static_cast<int>(runPairs.size());
                
                if (nAllPairs > 0) ++nRunsWithAnyDoNotScalePairs;
                else               ++nRunsWithoutAnyDoNotScalePairs;
                
                std::vector<LoadedPair> targetPairs;
                for (const auto& P : runPairs)
                {
                    if (P.groupFolder == targetGroup) targetPairs.push_back(P);
                }
                
                const int nTargetPairs = static_cast<int>(targetPairs.size());
                if (nTargetPairs > 0) ++nRunsWithTargetGroup;
                else                  ++nRunsWithoutTargetGroup;
                
                bool hasPhoton10 = false;
                bool hasPhoton12 = false;
                for (const auto& P : targetPairs)
                {
                    const int thr = extractPhotonThresholdGeV(P.probeKey);
                    if (thr == 10) hasPhoton10 = true;
                    if (thr == 12) hasPhoton12 = true;
                }
                
                if (hasPhoton10) ++nRunsWithPhoton10;
                if (hasPhoton12) ++nRunsWithPhoton12;
                if (hasPhoton10 && hasPhoton12) ++nRunsWithBothPhoton10And12;
                
                std::vector<LoadedPair> filtered = filterPhoton10And12FromGroup(runPairs, targetGroup);
                const int nFiltered = static_cast<int>(filtered.size());
                
                if (filtered.empty())
                {
                    if (nAllPairs <= 0)
                    {
                        perRunSummary.push_back(
                            TString::Format("RUN %s : NO_DONOTSCALE_PAIRS | allPairs=%d | targetGroupPairs=%d | photon10=%s | photon12=%s | path=%s",
                                            item.runLabel.c_str(),
                                            nAllPairs,
                                            nTargetPairs,
                                            hasPhoton10 ? "YES" : "NO",
                                            hasPhoton12 ? "YES" : "NO",
                                            item.filePath.c_str()).Data());
                    }
                    else if (nTargetPairs <= 0)
                    {
                        perRunSummary.push_back(
                            TString::Format("RUN %s : TARGET_GROUP_ABSENT | allPairs=%d | targetGroupPairs=%d | photon10=%s | photon12=%s | path=%s",
                                            item.runLabel.c_str(),
                                            nAllPairs,
                                            nTargetPairs,
                                            hasPhoton10 ? "YES" : "NO",
                                            hasPhoton12 ? "YES" : "NO",
                                            item.filePath.c_str()).Data());
                    }
                    else
                    {
                        ++nRunsTargetGroupButNoPhoton10Or12;
                        perRunSummary.push_back(
                            TString::Format("RUN %s : TARGET_GROUP_PRESENT_BUT_NO_PHOTON10_12 | allPairs=%d | targetGroupPairs=%d | photon10=%s | photon12=%s | path=%s",
                                            item.runLabel.c_str(),
                                            nAllPairs,
                                            nTargetPairs,
                                            hasPhoton10 ? "YES" : "NO",
                                            hasPhoton12 ? "YES" : "NO",
                                            item.filePath.c_str()).Data());
                    }
                    
                    fRun->Close();
                    delete fRun;
                    continue;
                }
                
                const std::string runOutDir = JoinPath(perRunGroupOutDir, item.runLabel);
                EnsureDir(runOutDir);
                
                drawGroupMaxClusterOverlay(filtered, runOutDir, "", item.runLabel);
                drawGroupTurnOnOverlay(filtered, runOutDir, "", item.runLabel);
                
                ++nRunsRenderable;
                perRunItems.push_back(item);
                perRunSummary.push_back(
                    TString::Format("RUN %s : RENDERED | allPairs=%d | targetGroupPairs=%d | photon10=%s | photon12=%s | renderPairs=%d | outDir=%s",
                                    item.runLabel.c_str(),
                                    nAllPairs,
                                    nTargetPairs,
                                    hasPhoton10 ? "YES" : "NO",
                                    hasPhoton12 ? "YES" : "NO",
                                    nFiltered,
                                    runOutDir.c_str()).Data());
                
                if (!badRunSet.empty())
                {
                    const int runNumInt = std::stoi(item.runLabel);
                    if (badRunSet.count(runNumInt))
                    {
                        ++nBadTriggerRunsSkipped;
                    }
                    else
                    {
                        for (const auto& P : filtered)
                        {
                            auto ait = accPairs.find(P.probeKey);
                            if (ait == accPairs.end())
                            {
                                LoadedPair LP;
                                LP.basisKey      = P.basisKey;
                                LP.basisLabel    = P.basisLabel;
                                LP.groupFolder   = P.groupFolder;
                                LP.probeKey      = P.probeKey;
                                LP.probeLabel    = P.probeLabel;
                                LP.baselineKey   = P.baselineKey;
                                LP.baselineLabel = P.baselineLabel;
                                LP.coverageText  = "badTrigRun-filtered combined";
                                LP.nRunSpans     = 1;
                                LP.hProbe = static_cast<TH1*>(P.hProbe->Clone(
                                    TString::Format("acc_probe_%s", P.probeKey.c_str())));
                                LP.hProbe->SetDirectory(nullptr);
                                LP.hBase = static_cast<TH1*>(P.hBase->Clone(
                                    TString::Format("acc_base_%s", P.probeKey.c_str())));
                                LP.hBase->SetDirectory(nullptr);
                                accPairs[P.probeKey] = LP;
                            }
                            else
                            {
                                ait->second.hProbe->Add(P.hProbe);
                                ait->second.hBase->Add(P.hBase);
                                ait->second.nRunSpans += 1;
                            }
                        }
                    }
                }
                
                fRun->Close();
                delete fRun;
            }
            
            perRunSummary.push_back("");
            perRunSummary.push_back("[COUNTS]");
            perRunSummary.push_back(TString::Format("Discovered runs: %d", nDiscoveredRuns).Data());
            perRunSummary.push_back(TString::Format("Openable runs: %d", nOpenableRuns).Data());
            perRunSummary.push_back(TString::Format("Bad/unopenable files: %d", nBadFiles).Data());
            perRunSummary.push_back(TString::Format("Runs with any doNotScale pairs: %d", nRunsWithAnyDoNotScalePairs).Data());
            perRunSummary.push_back(TString::Format("Runs with zero doNotScale pairs: %d", nRunsWithoutAnyDoNotScalePairs).Data());
            perRunSummary.push_back(TString::Format("Runs with target group %s: %d", targetGroup.c_str(), nRunsWithTargetGroup).Data());
            perRunSummary.push_back(TString::Format("Runs without target group %s: %d", targetGroup.c_str(), nRunsWithoutTargetGroup).Data());
            perRunSummary.push_back(TString::Format("Runs with Photon 10 in target group: %d", nRunsWithPhoton10).Data());
            perRunSummary.push_back(TString::Format("Runs with Photon 12 in target group: %d", nRunsWithPhoton12).Data());
            perRunSummary.push_back(TString::Format("Runs with both Photon 10 and 12 in target group: %d", nRunsWithBothPhoton10And12).Data());
            perRunSummary.push_back(TString::Format("Runs with target group but no Photon 10/12 pairs: %d", nRunsTargetGroupButNoPhoton10Or12).Data());
            perRunSummary.push_back(TString::Format("Renderable runs (at least one Photon 10/12 pair): %d", nRunsRenderable).Data());
            
            perRunSummary.push_back("");
            if (nDiscoveredRuns <= 0)
                perRunSummary.push_back("BREAKPOINT: directory opened, but no run files were discovered.");
            else if (nOpenableRuns <= 0)
                perRunSummary.push_back("BREAKPOINT: run files were discovered, but none could be opened.");
            else if (nRunsWithAnyDoNotScalePairs <= 0)
                perRunSummary.push_back("BREAKPOINT: run files opened, but none contained any doNotScale numerator/baseline pairs.");
            else if (nRunsWithTargetGroup <= 0)
                perRunSummary.push_back(string("BREAKPOINT: doNotScale pairs exist, but none matched target group ") + targetGroup + ".");
            else if (nRunsRenderable <= 0)
                perRunSummary.push_back("BREAKPOINT: target group exists in some runs, but no Photon 10/12 pairs were present to render.");
            else
                perRunSummary.push_back("BREAKPOINT: none reached; usable per-run outputs were produced.");
            
            if (!badRunSet.empty() && !accPairs.empty())
            {
                std::vector<LoadedPair> accVec;
                accVec.reserve(accPairs.size());
                for (auto& kv : accPairs) accVec.push_back(kv.second);
                std::sort(accVec.begin(), accVec.end(), [&](const LoadedPair& a, const LoadedPair& b)
                {
                    return extractPhotonThresholdGeV(a.probeKey) < extractPhotonThresholdGeV(b.probeKey);
                });
                
                const std::string convergentGroupOutDir = JoinPath(outDir, targetGroup);
                drawGroupMaxClusterOverlay(accVec, convergentGroupOutDir, "");
                drawGroupTurnOnOverlay(accVec, convergentGroupOutDir, "");
                
                const int nGoodRuns = nRunsRenderable - nBadTriggerRunsSkipped;
                cout << "\n" << ANSI_BOLD_GRN
                     << "========== BAD-TRIGGER-RUN FILTER SUMMARY ==========" << ANSI_RESET << "\n"
                     << "  Total renderable runs (with doNotScale data) : " << nRunsRenderable << "\n"
                     << "  Runs in badTriggerRuns vector                : " << static_cast<int>(badTriggerRuns.size()) << "\n"
                     << ANSI_BOLD_RED
                     << "  Runs removed (in vector AND renderable)      : " << nBadTriggerRunsSkipped << ANSI_RESET << "\n"
                     << ANSI_BOLD_GRN
                     << "  Runs used for combined overlays               : " << nGoodRuns << ANSI_RESET << "\n"
                     << ANSI_BOLD_GRN
                     << "====================================================" << ANSI_RESET << "\n\n";
                
                perRunSummary.push_back("");
                perRunSummary.push_back("[BAD-TRIGGER-RUN FILTER SUMMARY]");
                perRunSummary.push_back(TString::Format(
                    "  Total renderable runs: %d", nRunsRenderable).Data());
                perRunSummary.push_back(TString::Format(
                    "  Runs in badTriggerRuns vector: %d", static_cast<int>(badTriggerRuns.size())).Data());
                perRunSummary.push_back(TString::Format(
                    "  Runs removed: %d", nBadTriggerRunsSkipped).Data());
                perRunSummary.push_back(TString::Format(
                    "  Runs used for combined overlays: %d", nGoodRuns).Data());
                
                for (auto& kv : accPairs)
                {
                    delete kv.second.hProbe;
                    delete kv.second.hBase;
                }
                accPairs.clear();
            }
        }
        
        if (!scaledMissingByPattern.empty())
        {
            const std::string T1 = "MBD_NS_geq_2_vtx_lt_150";
            const std::string T2 = "photon_10_plus_MBD_NS_geq_2_vtx_lt_150";
            const std::string T3 = "photon_12_plus_MBD_NS_geq_2_vtx_lt_150";
            
            auto joinVec = [](const std::vector<std::string>& v) -> std::string {
                std::string s;
                for (std::size_t i = 0; i < v.size(); ++i) { if (i) s += ", "; s += v[i]; }
                return s;
            };
            
            auto runsWithBitsMissing = [&](int requiredBits) -> std::vector<std::string> {
                std::vector<std::string> out;
                for (const auto& [pat, runs] : scaledMissingByPattern)
                    if ((pat & requiredBits) == requiredBits)
                        out.insert(out.end(), runs.begin(), runs.end());
                std::sort(out.begin(), out.end());
                return out;
            };
            
            auto runsWithExactPattern = [&](int exact) -> std::vector<std::string> {
                auto it = scaledMissingByPattern.find(exact);
                if (it == scaledMissingByPattern.end()) return {};
                auto sorted = it->second;
                std::sort(sorted.begin(), sorted.end());
                return sorted;
            };
            
            cout << "\n" << ANSI_BOLD_YEL
                 << "========== SCALED-TRIGGER DIAGNOSTIC (cnt_ histogram check) =========="
                 << ANSI_RESET << "\n";
            
            cout << ANSI_BOLD_YEL
                 << "  [EXCLUSIVE PATTERNS] (exactly these missing, others present)"
                 << ANSI_RESET << "\n";
            
            struct PatDesc { int flag; std::string label; };
            const std::vector<PatDesc> excPatterns = {
                {1, "ONLY " + T1},
                {2, "ONLY " + T2},
                {4, "ONLY " + T3},
                {3, T1 + " + " + T2 + " only"},
                {5, T1 + " + " + T3 + " only"},
                {6, T2 + " + " + T3 + " only"},
                {7, "ALL THREE"}
            };
            for (const auto& pd : excPatterns)
            {
                auto runs = runsWithExactPattern(pd.flag);
                if (runs.empty()) continue;
                cout << ANSI_BOLD_RED << "    Missing " << pd.label
                     << " : " << runs.size() << " runs" << ANSI_RESET << "\n"
                     << "      " << joinVec(runs) << "\n";
            }
            
            cout << ANSI_BOLD_YEL
                 << "  [AGGREGATE] (any run missing at least these)"
                 << ANSI_RESET << "\n";
            
            struct AggDesc { int bits; std::string label; };
            const std::vector<AggDesc> aggPatterns = {
                {1, T1},
                {2, T2},
                {4, T3},
                {3, T1 + " AND " + T2},
                {5, T1 + " AND " + T3},
                {6, T2 + " AND " + T3},
                {7, T1 + " AND " + T2 + " AND " + T3}
            };
            for (const auto& ad : aggPatterns)
            {
                auto runs = runsWithBitsMissing(ad.bits);
                if (runs.empty()) continue;
                cout << "    Missing (at least) " << ad.label
                     << " : " << runs.size() << " runs\n"
                     << "      " << joinVec(runs) << "\n";
            }
            
            cout << ANSI_BOLD_YEL
                 << "======================================================================="
                 << ANSI_RESET << "\n\n";
            
            perRunSummary.push_back("");
            perRunSummary.push_back("[SCALED-TRIGGER DIAGNOSTIC]");
            for (const auto& pd : excPatterns)
            {
                auto runs = runsWithExactPattern(pd.flag);
                if (!runs.empty())
                    perRunSummary.push_back("  Missing " + pd.label + " : "
                        + std::to_string(runs.size()) + " runs : " + joinVec(runs));
            }
            for (const auto& ad : aggPatterns)
            {
                auto runs = runsWithBitsMissing(ad.bits);
                if (!runs.empty())
                    perRunSummary.push_back("  Missing (at least) " + ad.label + " : "
                        + std::to_string(runs.size()) + " runs : " + joinVec(runs));
            }
        }
        
        if (!badTriggerRuns.empty())
        {
            std::set<int> renderableRunNums;
            for (const auto& item : perRunItems)
                renderableRunNums.insert(std::stoi(item.runLabel));
            
            std::set<std::string> allMissingShortLabels;
            for (const auto& [pat, runs] : scaledMissingByPattern)
                allMissingShortLabels.insert(runs.begin(), runs.end());
            
            std::vector<int> badFound, badNotFound, badFoundAllTrigsPresent;
            std::map<int, std::vector<int>> badByPattern;
            
            for (int run : badTriggerRuns)
            {
                const std::string shortLbl = std::to_string(run);
                if (!renderableRunNums.count(run))
                {
                    badNotFound.push_back(run);
                    continue;
                }
                badFound.push_back(run);
                
                if (!allMissingShortLabels.count(shortLbl))
                {
                    badFoundAllTrigsPresent.push_back(run);
                    continue;
                }
                
                for (const auto& [pat, runs] : scaledMissingByPattern)
                {
                    for (const auto& r : runs)
                    {
                        if (r == shortLbl)
                        {
                            badByPattern[pat].push_back(run);
                            break;
                        }
                    }
                }
            }
            
            auto joinInts = [](const std::vector<int>& v) -> std::string {
                std::string s;
                for (std::size_t i = 0; i < v.size(); ++i) { if (i) s += ", "; s += std::to_string(v[i]); }
                return s;
            };
            
            const std::string T1 = "MBD_NS_geq_2_vtx_lt_150";
            const std::string T2 = "photon_10_plus_MBD_NS_geq_2_vtx_lt_150";
            const std::string T3 = "photon_12_plus_MBD_NS_geq_2_vtx_lt_150";
            
            struct PatDesc { int flag; std::string label; };
            const std::vector<PatDesc> patDescs = {
                {1, "ONLY " + T1},
                {2, "ONLY " + T2},
                {4, "ONLY " + T3},
                {3, T1 + " + " + T2 + " only"},
                {5, T1 + " + " + T3 + " only"},
                {6, T2 + " + " + T3 + " only"},
                {7, "ALL THREE"}
            };
            
            cout << ANSI_BOLD_YEL
                 << "========== badTriggerRuns VECTOR QA =========="
                 << ANSI_RESET << "\n";
            cout << "  Total runs in badTriggerRuns vector : " << badTriggerRuns.size() << "\n";
            cout << "  Found in renderable per-run files   : " << badFound.size() << "\n";
            cout << "  NOT found (missing or not renderable): " << badNotFound.size() << "\n";
            if (!badNotFound.empty())
                cout << ANSI_BOLD_RED << "    " << joinInts(badNotFound) << ANSI_RESET << "\n";
            
            cout << "  Found but ALL scaled triggers present: " << badFoundAllTrigsPresent.size() << "\n";
            if (!badFoundAllTrigsPresent.empty())
                cout << ANSI_BOLD_YEL << "    (flagged visually, not by cnt_) " << joinInts(badFoundAllTrigsPresent) << ANSI_RESET << "\n";
            
            int nWithMissing = 0;
            for (const auto& [pat, runs] : badByPattern) nWithMissing += static_cast<int>(runs.size());
            cout << "  Found with missing scaled triggers  : " << nWithMissing << "\n";
            
            for (const auto& pd : patDescs)
            {
                auto it = badByPattern.find(pd.flag);
                if (it == badByPattern.end() || it->second.empty()) continue;
                cout << ANSI_BOLD_RED << "    Missing " << pd.label
                     << " : " << it->second.size() << " runs" << ANSI_RESET << "\n"
                     << "      " << joinInts(it->second) << "\n";
            }
            
            cout << ANSI_BOLD_YEL
                 << "==============================================="
                 << ANSI_RESET << "\n\n";
            
            perRunSummary.push_back("");
            perRunSummary.push_back("[badTriggerRuns VECTOR QA]");
            perRunSummary.push_back("  Total in vector: " + std::to_string(badTriggerRuns.size()));
            perRunSummary.push_back("  Found renderable: " + std::to_string(badFound.size()));
            perRunSummary.push_back("  Not found: " + std::to_string(badNotFound.size())
                + (badNotFound.empty() ? "" : " : " + joinInts(badNotFound)));
            perRunSummary.push_back("  All triggers present (visual flag only): " + std::to_string(badFoundAllTrigsPresent.size())
                + (badFoundAllTrigsPresent.empty() ? "" : " : " + joinInts(badFoundAllTrigsPresent)));
            perRunSummary.push_back("  With missing scaled triggers: " + std::to_string(nWithMissing));
            for (const auto& pd : patDescs)
            {
                auto it = badByPattern.find(pd.flag);
                if (it != badByPattern.end() && !it->second.empty())
                    perRunSummary.push_back("    Missing " + pd.label + " : "
                        + std::to_string(it->second.size()) + " : " + joinInts(it->second));
            }
        }
        
        if (!perRunItems.empty())
        {
            const std::string maxTableDir = JoinPath(perRunGroupOutDir, "maxClusterOverlaysRunByRunTables");
            const std::string turnOnTableDir = JoinPath(perRunGroupOutDir, "turnOnRunByRunTables");
            writePerRunTablePages(perRunItems, targetGroup, maxTableDir, false, 12, 8, 5120, 2880);
            writePerRunTablePages(perRunItems, targetGroup, turnOnTableDir, true, 12, 8, 5120, 2880);
            
            if (!badRunSet.empty())
            {
                std::vector<PerRunTriggerAnaItem> badItems;
                for (const auto& item : perRunItems)
                {
                    const int runNum = std::stoi(item.runLabel);
                    if (badRunSet.count(runNum))
                        badItems.push_back(item);
                }
                if (!badItems.empty())
                {
                    const std::string badMaxTableDir = JoinPath(perRunGroupOutDir, "BAD_maxClusterOverlaysRunByRunTables");
                    writePerRunTablePages(badItems, targetGroup, badMaxTableDir, false, 10, 5, 5120, 2880);
                    
                    cout << ANSI_BOLD_YEL
                         << "[BAD-TRIGGER-RUN TABLES] Wrote " << badItems.size()
                         << " bad-run overlays to " << badMaxTableDir
                         << ANSI_RESET << "\n";
                }
            }
        }
        else
        {
            perRunSummary.push_back("No per-run files produced usable photon 10/12 common-basis overlays.");
        }
        
        WriteTextFile(JoinPath(perRunGroupOutDir, "summary.txt"), perRunSummary);
    }
    
    if (generatePerRunTriggerAna && datasetBucket() == "auau")
    {
        const std::string targetGroup = "commonBasis_MBD_NS_geq_2_vtx_lt_10";
        const std::string perRunInputDir = JoinPath(kInputBase + "/auau25/perRun", CfgTagWithUE());
        const std::string perRunGroupOutDir = JoinPath(JoinPath(outDir, targetGroup), "perRun");
        EnsureDir(perRunGroupOutDir);
        
        const std::set<int> badRunSet10(badTriggerRuns_vtx10.begin(), badTriggerRuns_vtx10.end());
        std::map<std::string, LoadedPair> accPairs10;
        int nBadSkipped10 = 0;
        int nRenderable10 = 0;
        std::vector<PerRunTriggerAnaItem> perRunItems10;
        
        void* dirp10 = gSystem->OpenDirectory(perRunInputDir.c_str());
        if (dirp10)
        {
            std::vector<PerRunTriggerAnaItem> discovered;
            const char* entry = nullptr;
            while ((entry = gSystem->GetDirEntry(dirp10)))
            {
                const std::string fname = entry;
                if (fname == "." || fname == "..") continue;
                if (fname.rfind("chunkMerge_run_", 0) != 0) continue;
                if (fname.size() <= std::string("chunkMerge_run_").size() + std::string(".root").size()) continue;
                if (fname.substr(fname.size() - 5) != ".root") continue;
                
                const std::string runDigits = fname.substr(std::string("chunkMerge_run_").size(),
                                                           fname.size() - std::string("chunkMerge_run_").size() - 5);
                if (runDigits.empty()) continue;
                
                PerRunTriggerAnaItem item;
                item.runLabel = runDigits;
                item.filePath = JoinPath(perRunInputDir, fname);
                discovered.push_back(item);
            }
            gSystem->FreeDirectory(dirp10);
            
            std::sort(discovered.begin(), discovered.end(),
                      [&](const PerRunTriggerAnaItem& a, const PerRunTriggerAnaItem& b)
                      { return a.runLabel < b.runLabel; });
            
            for (const auto& item : discovered)
            {
                TFile* fRun = TFile::Open(item.filePath.c_str(), "READ");
                if (!fRun || fRun->IsZombie())
                {
                    if (fRun) { fRun->Close(); delete fRun; }
                    continue;
                }
                
                std::vector<LoadedPair> runPairs = collectPairsFromFile(fRun, item.filePath, false);
                std::vector<LoadedPair> filtered = filterPhoton10And12FromGroup(runPairs, targetGroup);
                
                if (filtered.empty())
                {
                    fRun->Close();
                    delete fRun;
                    continue;
                }
                
                const std::string runOutDir = JoinPath(perRunGroupOutDir, item.runLabel);
                EnsureDir(runOutDir);
                drawGroupMaxClusterOverlay(filtered, runOutDir, "", item.runLabel);
                drawGroupTurnOnOverlay(filtered, runOutDir, "", item.runLabel);
                
                ++nRenderable10;
                perRunItems10.push_back(item);
                
                if (!badRunSet10.empty())
                {
                    const int runNumInt = std::stoi(item.runLabel);
                    if (badRunSet10.count(runNumInt))
                    {
                        ++nBadSkipped10;
                    }
                    else
                    {
                        for (const auto& P : filtered)
                        {
                            auto ait = accPairs10.find(P.probeKey);
                            if (ait == accPairs10.end())
                            {
                                LoadedPair LP;
                                LP.basisKey      = P.basisKey;
                                LP.basisLabel    = P.basisLabel;
                                LP.groupFolder   = P.groupFolder;
                                LP.probeKey      = P.probeKey;
                                LP.probeLabel    = P.probeLabel;
                                LP.baselineKey   = P.baselineKey;
                                LP.baselineLabel = P.baselineLabel;
                                LP.coverageText  = "badTrigRun-filtered vtx10 combined";
                                LP.nRunSpans     = 1;
                                LP.hProbe = static_cast<TH1*>(P.hProbe->Clone(
                                    TString::Format("acc10_probe_%s", P.probeKey.c_str())));
                                LP.hProbe->SetDirectory(nullptr);
                                LP.hBase = static_cast<TH1*>(P.hBase->Clone(
                                    TString::Format("acc10_base_%s", P.probeKey.c_str())));
                                LP.hBase->SetDirectory(nullptr);
                                accPairs10[P.probeKey] = LP;
                            }
                            else
                            {
                                ait->second.hProbe->Add(P.hProbe);
                                ait->second.hBase->Add(P.hBase);
                                ait->second.nRunSpans += 1;
                            }
                        }
                    }
                }
                
                fRun->Close();
                delete fRun;
            }
            
            if (!badRunSet10.empty() && !accPairs10.empty())
            {
                std::vector<LoadedPair> accVec;
                accVec.reserve(accPairs10.size());
                for (auto& kv : accPairs10) accVec.push_back(kv.second);
                std::sort(accVec.begin(), accVec.end(), [&](const LoadedPair& a, const LoadedPair& b)
                {
                    return extractPhotonThresholdGeV(a.probeKey) < extractPhotonThresholdGeV(b.probeKey);
                });
                
                const std::string convergentGroupOutDir = JoinPath(outDir, targetGroup);
                drawGroupMaxClusterOverlay(accVec, convergentGroupOutDir, "");
                drawGroupTurnOnOverlay(accVec, convergentGroupOutDir, "");
                
                const int nGood10 = nRenderable10 - nBadSkipped10;
                cout << "\n" << ANSI_BOLD_GRN
                     << "========== BAD-TRIGGER-RUN FILTER SUMMARY (vtx10) ==========" << ANSI_RESET << "\n"
                     << "  Total renderable runs (with doNotScale data) : " << nRenderable10 << "\n"
                     << "  Runs in badTriggerRuns_vtx10 vector           : " << static_cast<int>(badTriggerRuns_vtx10.size()) << "\n"
                     << ANSI_BOLD_RED
                     << "  Runs removed (in vector AND renderable)       : " << nBadSkipped10 << ANSI_RESET << "\n"
                     << ANSI_BOLD_GRN
                     << "  Runs used for combined overlays                : " << nGood10 << ANSI_RESET << "\n"
                     << ANSI_BOLD_GRN
                     << "=============================================================" << ANSI_RESET << "\n\n";
                
                for (auto& kv : accPairs10)
                {
                    delete kv.second.hProbe;
                    delete kv.second.hBase;
                }
                accPairs10.clear();
            }
        }
        
        if (!perRunItems10.empty())
        {
            const std::string maxTableDir = JoinPath(perRunGroupOutDir, "maxClusterOverlaysRunByRunTables");
            const std::string turnOnTableDir = JoinPath(perRunGroupOutDir, "turnOnRunByRunTables");
            writePerRunTablePages(perRunItems10, targetGroup, maxTableDir, false, 12, 8, 5120, 2880);
            writePerRunTablePages(perRunItems10, targetGroup, turnOnTableDir, true, 12, 8, 5120, 2880);
            
            if (!badRunSet10.empty())
            {
                std::vector<PerRunTriggerAnaItem> badItems;
                for (const auto& item : perRunItems10)
                {
                    const int runNum = std::stoi(item.runLabel);
                    if (badRunSet10.count(runNum))
                        badItems.push_back(item);
                }
                if (!badItems.empty())
                {
                    const std::string badMaxTableDir = JoinPath(perRunGroupOutDir, "BAD_maxClusterOverlaysRunByRunTables");
                    writePerRunTablePages(badItems, targetGroup, badMaxTableDir, false, 10, 5, 5120, 2880);
                }
            }
        }
    }
    
    // -------------------------------------------------------------------
    // NEW: Run24pp vs Run25pp 95% efficiency turn-on point overlay
    //
    // For each of the two common-basis groups:
    //   commonBasis_MBD_NandS_geq_1
    //   commonBasis_MBD_NandS_geq_1_vtx_lt_10
    //
    // Output:
    //   <outDir>/<group>/x95_overlay_run24pp_vs_run25pp.png
    //
    // X-axis: trigger label bins (Photon 2, 3, 4, 5, ...)
    // Y-axis: 95% efficiency turn-on point [GeV]
    // Markers: open circles = Run24pp, closed circles = Run25pp
    // Colors: match colorForProbe per photon threshold
    //
    // FORCEFULLY opens both InputPP() and InputPP25() regardless of toggles.
    // Skipped entirely in AuAu-only mode (no PP files to compare).
    // -------------------------------------------------------------------
    if (!isAuAuOnly)
    {
        const std::vector<std::string> targetGroups = {
            "commonBasis_MBD_NandS_geq_1",
            "commonBasis_MBD_NandS_geq_1_vtx_lt_10"
        };
        
        auto getHistFromFile = [&](TFile* f, const std::string& dirKey) -> TH1*
        {
            if (!f) return nullptr;
            TDirectory* dir = f->GetDirectory(dirKey.c_str());
            if (!dir) return nullptr;
            return dynamic_cast<TH1*>(dir->Get((prefix + dirKey).c_str()));
        };
        
        auto computeX95Map = [&](TFile* f, const std::string& groupTarget) -> std::map<int, double>
        {
            std::map<int, double> result;
            if (!f) return result;
            
            TIter nextKey(f->GetListOfKeys());
            while (TKey* key = dynamic_cast<TKey*>(nextKey()))
            {
                const std::string dirKey = key->GetName();
                if (dirKey.empty() || dirKey.find(baselineTag) == 0) continue;
                
                const std::string bKey = deduceBasisKey(dirKey);
                if (bKey.empty()) continue;
                if (groupFolderFromBasis(bKey) != groupTarget) continue;
                
                TH1* hProbe = getHistFromFile(f, dirKey);
                TH1* hBase  = getHistFromFile(f, baselineTag + dirKey);
                if (!hProbe || !hBase) continue;
                
                TH1* hP = CloneTH1(hProbe, TString::Format("x95ov_p_%s", dirKey.c_str()).Data());
                TH1* hB = CloneTH1(hBase,  TString::Format("x95ov_b_%s", dirKey.c_str()).Data());
                if (!hP || !hB) { if (hP) delete hP; if (hB) delete hB; continue; }
                
                EnsureSumw2(hP);
                EnsureSumw2(hB);
                
                TH1* hR = CloneTH1(hP, TString::Format("x95ov_r_%s", dirKey.c_str()).Data());
                if (!hR) { delete hP; delete hB; continue; }
                EnsureSumw2(hR);
                hR->Divide(hP, hB, 1.0, 1.0, "B");
                
                const int thr = extractPhotonThresholdGeV(dirKey);
                const double x95 = findXAtEff(hR, 0.95);
                if (x95 > 0.0) result[thr] = x95;
                
                delete hP;
                delete hB;
                delete hR;
            }
            return result;
        };
        
        TFile* f24 = TFile::Open(InputPP().c_str(), "READ");
        TFile* f25 = TFile::Open(InputPP25().c_str(), "READ");
        
        const bool ok24 = (f24 && !f24->IsZombie());
        const bool ok25 = (f25 && !f25->IsZombie());
        
        if (ok24 && ok25)
        {
            for (const auto& grp : targetGroups)
            {
                const auto x95_24 = computeX95Map(f24, grp);
                const auto x95_25 = computeX95Map(f25, grp);
                
                std::set<int> allThresholds;
                for (const auto& kv : x95_24) allThresholds.insert(kv.first);
                for (const auto& kv : x95_25) allThresholds.insert(kv.first);
                
                if (allThresholds.empty()) continue;
                
                const int nBins = (int)allThresholds.size();
                std::vector<int> thrs(allThresholds.begin(), allThresholds.end());
                
                TCanvas cOv(TString::Format("c_x95ov_%s", grp.c_str()).Data(),
                            "c_x95ov", 1000, 700);
                cOv.SetLeftMargin(0.14);
                cOv.SetRightMargin(0.05);
                cOv.SetBottomMargin(0.16);
                cOv.SetTopMargin(0.08);
                cOv.SetTicks(1, 1);
                
                TH1F frame(TString::Format("frame_x95ov_%s", grp.c_str()).Data(),
                           "", nBins, 0.5, nBins + 0.5);
                frame.SetTitle("");
                frame.GetXaxis()->SetTitle("Photon Trigger");
                frame.GetYaxis()->SetTitle("95% Efficiency Turn-on [GeV]");
                frame.GetYaxis()->SetTitleOffset(1.3);
                frame.GetXaxis()->SetLabelSize(0.050);
                frame.GetXaxis()->SetTitleSize(0.050);
                frame.GetYaxis()->SetLabelSize(0.045);
                frame.GetYaxis()->SetTitleSize(0.050);
                
                for (int i = 0; i < nBins; ++i)
                    frame.GetXaxis()->SetBinLabel(i + 1, TString::Format("Photon %d", thrs[i]).Data());
                
                double yMax = 0.0;
                for (const auto& kv : x95_24) yMax = std::max(yMax, kv.second);
                for (const auto& kv : x95_25) yMax = std::max(yMax, kv.second);
                
                frame.SetMinimum(0.0);
                frame.SetMaximum(yMax * 1.35);
                frame.Draw("axis");
                
                std::vector<TGraphErrors*> keepG;
                
                TGraphErrors gLeg24(0);
                gLeg24.SetMarkerStyle(24);
                gLeg24.SetMarkerColor(kBlack);
                gLeg24.SetMarkerSize(1.3);
                
                TGraphErrors gLeg25(0);
                gLeg25.SetMarkerStyle(20);
                gLeg25.SetMarkerColor(kBlack);
                gLeg25.SetMarkerSize(1.3);
                
                for (int i = 0; i < nBins; ++i)
                {
                    const int thr = thrs[i];
                    const std::string fakeKey = TString::Format("Photon_%d_GeV_dummy", thr).Data();
                    const int col = colorForProbe(fakeKey);
                    
                    auto it24 = x95_24.find(thr);
                    if (it24 != x95_24.end())
                    {
                        TGraphErrors* g = new TGraphErrors(1);
                        g->SetPoint(0, (double)(i + 1) - 0.12, it24->second);
                        g->SetPointError(0, 0.0, 0.0);
                        g->SetMarkerStyle(24);
                        g->SetMarkerSize(1.5);
                        g->SetMarkerColor(col);
                        g->SetLineColor(col);
                        g->Draw("P same");
                        keepG.push_back(g);
                    }
                    
                    auto it25 = x95_25.find(thr);
                    if (it25 != x95_25.end())
                    {
                        TGraphErrors* g = new TGraphErrors(1);
                        g->SetPoint(0, (double)(i + 1) + 0.12, it25->second);
                        g->SetPointError(0, 0.0, 0.0);
                        g->SetMarkerStyle(20);
                        g->SetMarkerSize(1.5);
                        g->SetMarkerColor(col);
                        g->SetLineColor(col);
                        g->Draw("P same");
                        keepG.push_back(g);
                    }
                }
                
                TLegend leg(0.18, 0.78, 0.50, 0.92);
                leg.SetBorderSize(0);
                leg.SetFillStyle(0);
                leg.SetTextFont(42);
                leg.SetTextSize(0.040);
                leg.AddEntry(&gLeg24, "Run24pp (open)", "p");
                leg.AddEntry(&gLeg25, "Run25pp (closed)", "p");
                leg.Draw();
                
                {
                    TLatex tGrp;
                    tGrp.SetNDC(true);
                    tGrp.SetTextFont(42);
                    tGrp.SetTextAlign(33);
                    tGrp.SetTextSize(0.035);
                    tGrp.DrawLatex(0.93, 0.97, grp.c_str());
                }
                
                const string grpOutDir = JoinPath(outDir, grp);
                EnsureDir(grpOutDir);
                SaveCanvas(cOv, JoinPath(grpOutDir, "x95_overlay_run24pp_vs_run25pp.png"));
                
                cout << ANSI_BOLD_GRN
                << "[WROTE] " << JoinPath(grpOutDir, "x95_overlay_run24pp_vs_run25pp.png")
                << ANSI_RESET << "\n";
                
                for (auto* g : keepG) delete g;
            }
        }
        else
        {
            cout << ANSI_BOLD_YEL
            << "[WARN] x95 overlay: could not open both PP files.\n"
            << "       PP24: " << InputPP() << (ok24 ? " (OK)" : " (MISSING)") << "\n"
            << "       PP25: " << InputPP25() << (ok25 ? " (OK)" : " (MISSING)")
            << ANSI_RESET << "\n";
        }
        
        if (f24) { f24->Close(); delete f24; }
        if (f25) { f25->Close(); delete f25; }
    } // end !isAuAuOnly guard
}
