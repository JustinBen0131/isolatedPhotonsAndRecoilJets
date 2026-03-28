void RunJES3QA(Dataset& ds)
{
    cout << ANSI_BOLD_CYN << "\n==============================\n"
         << "[SECTION 5F] JES3 QA (3D maps) (" << ds.label << ")\n"
         << "==============================" << ANSI_RESET << "\n";

    string outDir = ds.isSim
      ? JoinPath(ds.outBase, "RecoilJetQA/JES3")
      : JoinPath(ds.outBase, "baselineData/RecoilJetQA/JES3");

    EnsureDir(outDir);

    // ---------------------------------------------------------------------------
    // Local helpers (pure refactor: output + filenames + logic unchanged)
    // ---------------------------------------------------------------------------

    auto PrintMeansTableHeader =
      [&](const string& rKey, double R, int wPt, int wN, int wM)
    {
      cout << ANSI_BOLD_CYN
           << "\n[JES3 means] " << ds.label << "  rKey=" << rKey
           << " (R=" << std::fixed << std::setprecision(1) << R << ")\n"
           << ANSI_RESET;

      cout << std::left << std::setw(wPt) << "pTgamma bin"
           << std::right
           << std::setw(wN) << "N(reco)"
           << std::setw(wM) << "<xJ>_re"
           << std::setw(wM) << "<a>_re"
           << std::setw(wN) << "N(tru)"
           << std::setw(wM) << "<xJ>_tr"
           << std::setw(wM) << "<a>_tr"
           << "\n";
      cout << string(wPt + 2*wN + 4*wM, '-') << "\n";
    };

    auto InitSummaryLines =
      [&](const string& rKey, double R)->vector<string>
    {
      vector<string> sumLines;
      sumLines.push_back(string("JES3 summary (") + ds.label + ")");
      sumLines.push_back(string("rKey: ") + rKey + TString::Format("  R=%.1f", R).Data());
      sumLines.push_back("");
      sumLines.push_back("NOTE (requested):");
      sumLines.push_back("  RECO truth-match (MC only): |#eta_{reco}|<0.7, p_{T,reco}>5 GeV, #DeltaR<0.05, best match via CaloRawClusterEval.");
      sumLines.push_back("  TRUTH prompt #gamma (MC): |#eta_{truth}|<0.7, PID=22 final-state, direct/frag via HEPMC history.");
      sumLines.push_back("  TRUTH iso (MC): E_{T}^{iso,truth} = #Sigma E_{T}(#DeltaR<0.3, excl. #gamma).");
      sumLines.push_back("");
      return sumLines;
    };

    auto PickPtAxis =
      [&](TH3* h1, TH3* h2, TH3* h3, TH3* h4,
          TH3* h5, TH3* h6, TH3* h7, TH3* h8)->const TAxis*
    {
      return
        (h1 ? h1->GetXaxis()
        : (h2 ? h2->GetXaxis()
        : (h3 ? h3->GetXaxis()
        : (h4 ? h4->GetXaxis()
        : (h5 ? h5->GetXaxis()
        : (h6 ? h6->GetXaxis()
        : (h7 ? h7->GetXaxis()
        : (h8 ? h8->GetXaxis() : nullptr))))))));
    };

    struct Jes3Dirs
    {
      string rOut;
      string dir2D;
      string dirSumm;

      // xJ(=Y) projections after integrating over alpha(=Z)
      string dirXJProj;

      // Organize JES3 overlays + efficiency dashboards under xJ_fromJES3
      string dirXJProjOverlay;
      string dirXJProjEff;
      string dirXJProjEffLeadMatch;
      string dirXJProjEffLeadJetResponse;
      string dirXJProjEffTagFractions3x3;

      // RECO projections
      string dirXJProjReco;                 // (1) h_JES3_pT_xJ_alpha_<rKey>
      string dirXJProjRecoTruthPhoTagged;   // (7) h_JES3RecoTruthPhoTagged_pT_xJ_alpha_<rKey>
      string dirXJProjRecoTruthTagged;      // (8) h_JES3RecoTruthTagged_pT_xJ_alpha_<rKey>

      // TRUTH projections
      string dirXJProjTruthRecoCond;        // (5) h_JES3Truth_pT_xJ_alpha_<rKey>  (jet-matched truth subset)
      string dirXJProjTruthRecoCondNoJetMatch; // (4) h_JES3TruthRecoCondNoJetMatch_pT_xJ_alpha_<rKey>
      string dirXJProjTruthPure;            // (3) h_JES3TruthPure_pT_xJ_alpha_<rKey>
    };

    auto MakeJes3Dirs =
      [&](const string& rKey, double /*R*/)->Jes3Dirs
    {
      Jes3Dirs D;
      D.rOut     = JoinPath(outDir, rKey);
      D.dir2D    = JoinPath(D.rOut, "2DMaps");
      D.dirSumm  = JoinPath(D.rOut, "summaries");

      D.dirXJProj                       = JoinPath(D.rOut, "xJ_fromJES3");

      // Keep *all* JES3 overlays under xJ_fromJES3/Overlay/
      D.dirXJProjOverlay                = JoinPath(D.dirXJProj, "Overlay");

      // SIM-only efficiency dashboards live under xJ_fromJES3/Efficiency/
      D.dirXJProjEff                    = JoinPath(D.dirXJProj, "Efficiency");
      D.dirXJProjEffLeadMatch           = JoinPath(D.dirXJProjEff, "LeadTruthRecoilMatch");
      D.dirXJProjEffLeadJetResponse     = JoinPath(D.dirXJProjEff, "LeadJetResponse");
      D.dirXJProjEffTagFractions3x3     = JoinPath(D.dirXJProjEff, "JES3_TagFractions_3x3");

      D.dirXJProjReco                   = JoinPath(D.dirXJProj, "RECO");
      D.dirXJProjRecoTruthPhoTagged     = JoinPath(D.dirXJProj, "RECO_truthPhoTagged");
      D.dirXJProjRecoTruthTagged        = JoinPath(D.dirXJProj, "RECO_truthTaggedPhoJet");

      D.dirXJProjTruthRecoCond           = JoinPath(D.dirXJProj, "TRUTH_recoConditioned");
      D.dirXJProjTruthRecoCondNoJetMatch = JoinPath(D.dirXJProj, "TRUTH_recoConditioned_noJetMatch");
      D.dirXJProjTruthPure               = JoinPath(D.dirXJProj, "TRUTH_pure");

      EnsureDir(D.rOut);
      EnsureDir(D.dir2D);
      EnsureDir(D.dirSumm);

        EnsureDir(D.dirXJProj);

      // Avoid cluttering DATA trees with empty Overlay/ and Efficiency/ folders.
      if (ds.isSim)
      {
        EnsureDir(D.dirXJProjOverlay);

        EnsureDir(D.dirXJProjEff);
        EnsureDir(D.dirXJProjEffLeadMatch);
        EnsureDir(D.dirXJProjEffLeadJetResponse);
        EnsureDir(D.dirXJProjEffTagFractions3x3);
      }

      // Always keep RECO outputs for both SIM and DATA
      EnsureDir(D.dirXJProjReco);

      // SIM-only outputs (avoid empty TRUTH/truth-tag folders in DATA)
      if (ds.isSim)
      {
        EnsureDir(D.dirXJProjRecoTruthPhoTagged);
        EnsureDir(D.dirXJProjRecoTruthTagged);

        EnsureDir(D.dirXJProjTruthRecoCond);
        EnsureDir(D.dirXJProjTruthRecoCondNoJetMatch);
        EnsureDir(D.dirXJProjTruthPure);
      }

      return D;
    };

    struct Jes3Hists
    {
      // RECO observables
      TH3* hReco_xJ               = nullptr;  // (1)
      TH3* hReco_j1               = nullptr;  // (2)

      TH3* hRecoTruthPhoTagged_xJ = nullptr;  // (7)
      TH3* hRecoTruthTagged_xJ    = nullptr;  // (8)

      // TRUTH observables
      TH3* hTrut_xJ               = nullptr;  // (5) reco-conditioned truth, jet-matched subset
      TH3* hTrutNoJM_xJ           = nullptr;  // (4) reco-conditioned truth, NO jet match
      TH3* hTrutPure_xJ           = nullptr;  // (3) pure truth
      TH3* hTrut_j1               = nullptr;  // (6)
    };

    auto LoadJes3Hists =
      [&](const string& rKey)->Jes3Hists
    {
      Jes3Hists H;

      // Baseline RECO (exists in SIM and DATA)
      H.hReco_xJ               = GetObj<TH3>(ds, "h_JES3_pT_xJ_alpha_" + rKey, true, true, true);
      H.hReco_j1               = GetObj<TH3>(ds, "h_JES3_pT_jet1Pt_alpha_" + rKey, true, true, true);

      // SIM-only objects: do NOT attempt to load / log-miss these in DATA mode
      if (ds.isSim)
      {
        // TRUTH (pure + reco-conditioned variants)
        H.hTrut_xJ               = GetObj<TH3>(ds, "h_JES3Truth_pT_xJ_alpha_" + rKey, true, true, true);
        H.hTrutNoJM_xJ           = GetObj<TH3>(ds, "h_JES3TruthRecoCondNoJetMatch_pT_xJ_alpha_" + rKey, true, true, true);
        H.hTrutPure_xJ           = GetObj<TH3>(ds, "h_JES3TruthPure_pT_xJ_alpha_" + rKey, true, true, true);
        H.hTrut_j1               = GetObj<TH3>(ds, "h_JES3Truth_pT_jet1Pt_alpha_" + rKey, true, true, true);

        // TRUTH-conditioned RECO (tagged subsets)
        H.hRecoTruthPhoTagged_xJ = GetObj<TH3>(ds, "h_JES3RecoTruthPhoTagged_pT_xJ_alpha_" + rKey, true, true, true);
        H.hRecoTruthTagged_xJ    = GetObj<TH3>(ds, "h_JES3RecoTruthTagged_pT_xJ_alpha_" + rKey, true, true, true);
      }

      return H;
    };


    auto HasAnyJes3 =
      [&](const Jes3Hists& H)->bool
    {
      return (
        H.hReco_xJ || H.hReco_j1 ||
        H.hRecoTruthPhoTagged_xJ || H.hRecoTruthTagged_xJ ||
        H.hTrut_xJ || H.hTrutNoJM_xJ || H.hTrutPure_xJ || H.hTrut_j1
      );
    };


    // ---------------------------------------------------------------------------
    // Per-rKey JES3 (main body)
    // ---------------------------------------------------------------------------

    auto ProcessOneRKey =
      [&](const string& rKey)
    {
      const double R = RFromKey(rKey);

      const Jes3Dirs  D = MakeJes3Dirs(rKey, R);
      const Jes3Hists H = LoadJes3Hists(rKey);

      if (!HasAnyJes3(H))
      {
        cout << ANSI_BOLD_YEL << "[WARN] No JES3 TH3 histograms for " << ds.label << " " << rKey << ANSI_RESET << "\n";
        return;
      }

      // pT axis from available hist
      const TAxis* axPt = PickPtAxis(
          H.hReco_xJ,
          H.hReco_j1,
          H.hRecoTruthPhoTagged_xJ,
          H.hRecoTruthTagged_xJ,
          H.hTrut_xJ,
          H.hTrutNoJM_xJ,
          H.hTrut_j1,
          H.hTrutPure_xJ
      );
      const int nPt = (axPt ? axPt->GetNbins() : 0);

      if (nPt <= 0)
      {
        cout << ANSI_BOLD_YEL << "[WARN] JES3: could not determine pT binning (nPt<=0) for " << ds.label << " " << rKey << ANSI_RESET << "\n";
        return;
      }

      const int wPt = 16;
      const int wN  = 12;
      const int wM  = 12;

      PrintMeansTableHeader(rKey, R, wPt, wN, wM);

      vector<string> sumLines = InitSummaryLines(rKey, R);

      struct BinPack
      {
        string ptLab;

        double nRe = 0, mxRe = 0, maRe = 0;
        double nTr = 0, mxTr = 0, maTr = 0;

        TH1* xJ_re = nullptr;
        TH1* a_re  = nullptr;
        TH1* xJ_tr = nullptr;
        TH1* a_tr  = nullptr;
      };

      auto BuildBinPack =
        [&](int ib)->BinPack
      {
        BinPack P;
        P.ptLab = AxisBinLabel(axPt, ib, "GeV", 0);

        if (H.hReco_xJ)
        {
          P.xJ_re = ProjectY_AtXbin_TH3(H.hReco_xJ, ib, TString::Format("jes3_xJ_re_%s_%d", rKey.c_str(), ib).Data());
          P.a_re  = ProjectZ_AtXbin_TH3(H.hReco_xJ, ib, TString::Format("jes3_a_re_%s_%d", rKey.c_str(), ib).Data());

          if (P.xJ_re)
          {
            P.xJ_re->SetName(TString::Format("h_xJ_re_%s_pTbin%d", rKey.c_str(), ib).Data());
            P.nRe  = P.xJ_re->GetEntries();
            P.mxRe = P.xJ_re->GetMean();
          }
          if (P.a_re)
          {
            P.a_re->SetName(TString::Format("h_a_re_%s_pTbin%d", rKey.c_str(), ib).Data());
            P.maRe = P.a_re->GetMean();
          }
        }

        if (H.hTrut_xJ)
        {
          P.xJ_tr = ProjectY_AtXbin_TH3(H.hTrut_xJ, ib, TString::Format("jes3_xJ_tr_%s_%d", rKey.c_str(), ib).Data());
          P.a_tr  = ProjectZ_AtXbin_TH3(H.hTrut_xJ, ib, TString::Format("jes3_a_tr_%s_%d", rKey.c_str(), ib).Data());

          if (P.xJ_tr)
          {
            P.xJ_tr->SetName(TString::Format("h_xJ_tr_%s_pTbin%d", rKey.c_str(), ib).Data());
            P.nTr  = P.xJ_tr->GetEntries();
            P.mxTr = P.xJ_tr->GetMean();
          }
          if (P.a_tr)
          {
            P.a_tr->SetName(TString::Format("h_a_tr_%s_pTbin%d", rKey.c_str(), ib).Data());
            P.maTr = P.a_tr->GetMean();
          }
        }

        return P;
      };

      auto CleanupBinPack =
        [&](BinPack& P)
      {
        if (P.xJ_re) delete P.xJ_re;
        if (P.a_re)  delete P.a_re;
        if (P.xJ_tr) delete P.xJ_tr;
        if (P.a_tr)  delete P.a_tr;

        P.xJ_re = nullptr;
        P.a_re  = nullptr;
        P.xJ_tr = nullptr;
        P.a_tr  = nullptr;
      };

        auto SaveXJRecoPNGs =
            [&](TH1* xJ_re, int ib, const string& ptLab)
          {
            if (!xJ_re) return;

            // (1) Inclusive in alpha
            const double ptMinGamma = H.hReco_xJ->GetXaxis()->GetBinLowEdge(ib);
            const double ptMaxGamma = H.hReco_xJ->GetXaxis()->GetBinUpEdge(ib);

            DrawAndSave_xJRecoIntegratedAlpha_WithFloors(ds, xJ_re,
                JoinPath(D.dirXJProjReco, TString::Format("xJ_reco_integratedAlpha_pTbin%d.png", ib).Data()),
                ptMinGamma, ptMaxGamma, R, 0.0, false);

            // (1b) NEW: overlayedWithSim (DATA reco vs SIM reco, same rKey/pT bin)
            if (isSimAndDataPP && !ds.isSim)
            {
                const string dirOv = JoinPath(D.dirXJProjReco, "insituCalib");
                EnsureDir(dirOv);

                static std::string s_lastSimPath = "";
                static TFile* s_fSim = nullptr;
                static TDirectory* s_simTopDir = nullptr;

                const std::string simPath = SimInputPathForSample(CurrentSimSample());
                if (!simPath.empty())
                {
                  if (!s_fSim || s_lastSimPath != simPath)
                  {
                    if (s_fSim) { s_fSim->Close(); delete s_fSim; s_fSim = nullptr; }
                    s_simTopDir = nullptr;

                    s_fSim = TFile::Open(simPath.c_str(), "READ");
                      if (s_fSim && !s_fSim->IsZombie())
                      {
                        s_simTopDir = s_fSim->GetDirectory(kDirSIM.c_str());
                        if (!s_simTopDir) s_simTopDir = s_fSim;
                      }
                    s_lastSimPath = simPath;
                  }
                }

                TH3* hSim3 = (s_simTopDir ? dynamic_cast<TH3*>(s_simTopDir->Get(("h_JES3_pT_xJ_alpha_" + rKey).c_str())) : nullptr);
              if (hSim3)
              {
                TH1* xJ_sim_raw = ProjectY_AtXbin_AndAlphaMax_TH3(
                  hSim3, ib, hSim3->GetZaxis()->GetXmax(),
                  TString::Format("jes3_xJ_sim_%s_%d_intAlpha", rKey.c_str(), ib).Data()
                );

                if (xJ_sim_raw)
                {
                  xJ_sim_raw->SetDirectory(nullptr);

                  TH1* xJ_dat = CloneTH1(xJ_re, TString::Format("jes3_xJ_dat_%s_%d_intAlpha", rKey.c_str(), ib).Data());
                  TH1* xJ_sim = CloneTH1(xJ_sim_raw, TString::Format("jes3_xJ_sim_%s_%d_intAlpha_clone", rKey.c_str(), ib).Data());
                  delete xJ_sim_raw;

                  if (xJ_dat && xJ_sim)
                  {
                    EnsureSumw2(xJ_dat);
                    EnsureSumw2(xJ_sim);

                    const double iDat = xJ_dat->Integral(0, xJ_dat->GetNbinsX() + 1);
                    const double iSim = xJ_sim->Integral(0, xJ_sim->GetNbinsX() + 1);
                    if (iDat > 0.0) xJ_dat->Scale(1.0 / iDat);
                    if (iSim > 0.0) xJ_sim->Scale(1.0 / iSim);

                    TCanvas c("c_xJRecoOv","c_xJRecoOv",900,700);
                    ApplyCanvasMargins1D(c);
                    c.SetLogy(false);

                    xJ_dat->SetTitle("");
                    xJ_dat->SetLineWidth(2);
                    xJ_dat->SetLineColor(kGreen + 2);
                    xJ_dat->SetMarkerStyle(20);
                    xJ_dat->SetMarkerSize(1.0);
                    xJ_dat->SetMarkerColor(kGreen + 2);

                    xJ_sim->SetLineWidth(2);
                    xJ_sim->SetLineColor(kOrange + 7);
                    xJ_sim->SetMarkerStyle(20);
                    xJ_sim->SetMarkerSize(1.0);
                    xJ_sim->SetMarkerColor(kOrange + 7);

                    xJ_dat->GetXaxis()->SetTitle("x_{J#gamma}");
                    xJ_dat->GetXaxis()->SetRangeUser(0.0, 2.0);
                    xJ_dat->GetYaxis()->SetTitle("Normalized counts");

                    xJ_dat->Draw("E1");
                    xJ_sim->Draw("E1 same");
                    gPad->Update();

                    const double jetPtMin_GeV = static_cast<double>(kJetPtMin);
                    const string bbLabel = B2BLabel();

                    const double xAbs  = (ptMaxGamma > 0.0) ? (jetPtMin_GeV / ptMaxGamma) : -1.0;
                    const double xFull = (ptMinGamma > 0.0) ? (jetPtMin_GeV / ptMinGamma) : -1.0;

                    const double yMin = gPad->GetUymin();
                    const double yMax = gPad->GetUymax();

                    TLine* lnAbs = new TLine(xAbs,  yMin, xAbs,  yMax);
                    lnAbs->SetLineColor(kBlue + 1);
                    lnAbs->SetLineStyle(2);
                    lnAbs->SetLineWidth(2);

                    TLine* lnFull = new TLine(xFull, yMin, xFull, yMax);
                    lnFull->SetLineColor(kRed + 1);
                    lnFull->SetLineStyle(2);
                    lnFull->SetLineWidth(2);

                    if (xAbs > 0.0)  lnAbs->Draw("same");
                    if (xFull > 0.0) lnFull->Draw("same");

                      TLegend* leg = new TLegend(0.62, 0.70, 0.95, 0.92);
                      leg->SetBorderSize(0);
                      leg->SetFillStyle(0);
                      leg->SetTextFont(42);
                      leg->SetTextSize(0.034);

                      leg->AddEntry(xJ_dat, "DATA (reco)", "ep");
                      leg->AddEntry(xJ_sim, "SIM (reco)",  "ep");
                      if (xAbs > 0.0)
                        leg->AddEntry(lnAbs,
                          TString::Format("x_{J, min}^{abs} = #frac{%.0f}{p_{T, max}^{#gamma}} = %.3f", jetPtMin_GeV, xAbs),
                          "l");
                      if (xFull > 0.0)
                        leg->AddEntry(lnFull,
                          TString::Format("x_{J, min}^{full} = #frac{%.0f}{p_{T, min}^{#gamma}} = %.3f", jetPtMin_GeV, xFull),
                          "l");

                      leg->Draw();

                      {
                        TLatex tCuts;
                        tCuts.SetNDC(true);
                        tCuts.SetTextFont(42);
                        tCuts.SetTextAlign(33);
                        tCuts.SetTextSize(0.038);
                        tCuts.DrawLatex(0.92, 0.62, TString::Format("|#Delta#phi(#gamma,jet)| > %s", bbLabel.c_str()).Data());
                        tCuts.DrawLatex(0.92, 0.54, TString::Format("p_{T}^{jet} > %.0f GeV", jetPtMin_GeV).Data());
                        tCuts.DrawLatex(0.92, 0.46, TString::Format("|v_{z}| < %.0f cm", vzCutCm).Data());
                      }

                    TLatex ttl;
                    ttl.SetNDC(true);
                    ttl.SetTextFont(42);
                    ttl.SetTextSize(0.052);
                    ttl.DrawLatex(0.12, 0.94,
                      TString::Format("RECO x_{J#gamma} (DATA vs SIM), p_{T}^{#gamma} = %.0f - %.0f GeV, R = %.1f",
                        ptMinGamma, ptMaxGamma, R).Data());

                    SaveCanvas(c, JoinPath(dirOv,
                      TString::Format("xJ_reco_integratedAlpha_overlayedWithSim_pTbin%d.png", ib).Data()));

                    delete leg;
                    delete lnFull;
                    delete lnAbs;
                  }

                  if (xJ_dat) delete xJ_dat;
                  if (xJ_sim) delete xJ_sim;
                }
              }
            }

                              
            // (2) NEW: alpha-cut variants (presentation-driven)
            // Adjust cut values freely; these are "reasonable" to see shape evolution.
            const vector<double> alphaMaxCuts = {0.20, 0.30, 0.40, 0.50};

            auto AlphaTag = [&](double aMax)->string
            {
              std::ostringstream s;
              s << std::fixed << std::setprecision(2) << aMax;  // e.g. "0.20"
              string t = s.str();
              std::replace(t.begin(), t.end(), '.', 'p');       // -> "0p20"
              return string("alphaLT") + t;                     // -> "alphaLT0p20"
            };

            const string dirAlphaBase = JoinPath(D.dirXJProjReco, "alphaCuts");
            EnsureDir(dirAlphaBase);

            for (double aMax : alphaMaxCuts)
            {
              const string aTag = AlphaTag(aMax);
              const string aDir = JoinPath(dirAlphaBase, aTag);
              EnsureDir(aDir);

              TH1* xJ_cut = ProjectY_AtXbin_AndAlphaMax_TH3(
                H.hReco_xJ, ib, aMax,
                TString::Format("jes3_xJ_re_%s_%d_%s", rKey.c_str(), ib, aTag.c_str()).Data()
              );

              if (!xJ_cut) continue;

              vector<string> linesCut = {
                "JES3 (RECO): x_{J#gamma}",
                TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data(),
                TString::Format("p_{T}^{#gamma}: %s", ptLab.c_str()).Data(),
                TString::Format("#alpha cut: #alpha < %.2f", aMax).Data(),
                "Projection: integrate only selected #alpha bins"
              };

              DrawAndSaveTH1_Common(ds, xJ_cut,
                JoinPath(aDir, TString::Format("xJ_reco_alphaCut_%s_pTbin%d.png", aTag.c_str(), ib).Data()),
                "x_{J#gamma}", "Counts", linesCut, false, false, 0.0, "E1");

              delete xJ_cut;
            }
        };

        auto SaveXJTruthPNGs =
          [&](TH1* xJ_tr, int ib, const string& ptLab)
        {
          if (!xJ_tr) return;

          vector<string> lines = {
            "JES3 (TRUTH reco-conditioned, jet-matched): x_{J#gamma}^{truth}",
            TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data(),
            TString::Format("p_{T}^{#gamma,truth}: %s", ptLab.c_str()).Data(),
            "Filled with: (tPt, xJt=tj1Pt/tPt, aT=tj2Pt/tPt)",
            "NOTE: this TRUTH TH3 is filled only when reco jet1 matches truth jet1 (ΔR<=0.3)"
          };

          DrawAndSaveTH1_Common(ds, xJ_tr,
            JoinPath(D.dirXJProjTruthRecoCond, TString::Format("xJ_truth_integratedAlpha_recoConditioned_pTbin%d.png", ib).Data()),
            "x_{J#gamma}^{truth}", "Counts", lines, false, false, 0.0, "E1");
        };

        auto SaveXJTruthNoJetMatchPNGs =
          [&](int ib)
        {
          if (!H.hTrutNoJM_xJ) return;

          TH1* xJ_nojm = ProjectY_AtXbin_TH3(
            H.hTrutNoJM_xJ, ib,
            TString::Format("jes3_xJ_trNoJM_%s_%d", rKey.c_str(), ib).Data()
          );

          if (!xJ_nojm) return;

          xJ_nojm->SetName(TString::Format("h_xJ_trNoJM_%s_pTbin%d", rKey.c_str(), ib).Data());

          const string ptLab = AxisBinLabel(H.hTrutNoJM_xJ->GetXaxis(), ib, "GeV", 0);

          vector<string> lines = {
            "JES3 (TRUTH reco-conditioned, NO jet match): x_{J#gamma}^{truth}",
            TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data(),
            TString::Format("p_{T}^{#gamma,truth}: %s", ptLab.c_str()).Data(),
            "Filled with: (tPt, xJt=tj1Pt/tPt, aT=tj2Pt/tPt)",
            "NOTE: reco success + truth-photon matched, but NO reco-jet↔truth-jet matching"
          };

          DrawAndSaveTH1_Common(ds, xJ_nojm,
            JoinPath(D.dirXJProjTruthRecoCondNoJetMatch,
              TString::Format("xJ_truth_integratedAlpha_recoConditioned_noJetMatch_pTbin%d.png", ib).Data()),
            "x_{J#gamma}^{truth}", "Counts", lines, false, false, 0.0, "E1");

          delete xJ_nojm;
        };

        auto SaveXJTruthPurePNGs =
          [&](int ib)
        {
          if (!H.hTrutPure_xJ) return;

          TH1* xJ_pure = ProjectY_AtXbin_TH3(
            H.hTrutPure_xJ, ib,
            TString::Format("jes3_xJ_trPure_%s_%d", rKey.c_str(), ib).Data()
          );

          if (!xJ_pure) return;

          xJ_pure->SetName(TString::Format("h_xJ_trPure_%s_pTbin%d", rKey.c_str(), ib).Data());

          const string ptLab = AxisBinLabel(H.hTrutPure_xJ->GetXaxis(), ib, "GeV", 0);

          vector<string> lines = {
            "JES3 (TRUTH pure): x_{J#gamma}^{truth}",
            TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data(),
            TString::Format("p_{T}^{#gamma,truth}: %s", ptLab.c_str()).Data(),
            "Filled with: (tPt, xJt=tj1Pt/tPt, aT=tj2Pt/tPt)",
            "NOTE: no reco gating / no reco-jet matching"
          };

          DrawAndSaveTH1_Common(ds, xJ_pure,
            JoinPath(D.dirXJProjTruthPure, TString::Format("xJ_truth_integratedAlpha_pureTruth_pTbin%d.png", ib).Data()),
            "x_{J#gamma}^{truth}", "Counts", lines, false, false, 0.0, "E1");

          delete xJ_pure;
        };

        auto SaveXJRecoTruthPhoTaggedPNGs =
          [&](int ib)
        {
          if (!H.hRecoTruthPhoTagged_xJ) return;

          TH1* xJ_tag = ProjectY_AtXbin_TH3(
            H.hRecoTruthPhoTagged_xJ, ib,
            TString::Format("jes3_xJ_rePhoTag_%s_%d", rKey.c_str(), ib).Data()
          );
          if (!xJ_tag) return;

          xJ_tag->SetName(TString::Format("h_xJ_rePhoTag_%s_pTbin%d", rKey.c_str(), ib).Data());

          const string ptLab = AxisBinLabel(H.hRecoTruthPhoTagged_xJ->GetXaxis(), ib, "GeV", 0);

          vector<string> lines = {
            "JES3 (RECO truth-PHOTON-tagged): x_{J#gamma}",
            TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data(),
            TString::Format("p_{T}^{#gamma}: %s", ptLab.c_str()).Data(),
            "Reco fill, but only when reco photon is matched to a truth signal photon",
            "Filled with: (p_{T}^{#gamma}, x_{J#gamma}, #alpha)"
          };

          DrawAndSaveTH1_Common(ds, xJ_tag,
            JoinPath(D.dirXJProjRecoTruthPhoTagged,
              TString::Format("xJ_reco_integratedAlpha_truthPhoTagged_pTbin%d.png", ib).Data()),
            "x_{J#gamma}", "Counts", lines, false, false, 0.0, "E1");

          delete xJ_tag;
        };

        auto SaveXJRecoTruthTaggedPNGs =
          [&](int ib)
        {
          if (!H.hRecoTruthTagged_xJ) return;

          TH1* xJ_tag = ProjectY_AtXbin_TH3(
            H.hRecoTruthTagged_xJ, ib,
            TString::Format("jes3_xJ_reTag_%s_%d", rKey.c_str(), ib).Data()
          );
          if (!xJ_tag) return;

          xJ_tag->SetName(TString::Format("h_xJ_reTag_%s_pTbin%d", rKey.c_str(), ib).Data());

          const string ptLab = AxisBinLabel(H.hRecoTruthTagged_xJ->GetXaxis(), ib, "GeV", 0);

          vector<string> lines = {
            "JES3 (RECO truth-tagged PHOTON+JET): x_{J#gamma}",
            TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data(),
            TString::Format("p_{T}^{#gamma}: %s", ptLab.c_str()).Data(),
            "Reco fill, but only when: truth-signal photon matched AND reco jet1 matches truth jet1 (ΔR<=0.3)",
            "Filled with: (p_{T}^{#gamma}, x_{J#gamma}, #alpha)"
          };

          DrawAndSaveTH1_Common(ds, xJ_tag,
            JoinPath(D.dirXJProjRecoTruthTagged,
              TString::Format("xJ_reco_integratedAlpha_truthTaggedPhoJet_pTbin%d.png", ib).Data()),
            "x_{J#gamma}", "Counts", lines, false, false, 0.0, "E1");

          delete xJ_tag;
        };

      auto SaveJes3Maps2D_ForBin =
        [&](int ib, const string& ptLab)
      {
        // Existing behavior: RECO/TRUTH xJ-alpha maps and jet1Pt-alpha maps

        if (H.hReco_xJ)
        {
          TH2* h2 = ProjectYZ_AtXbin_TH3(H.hReco_xJ, ib, TString::Format("jes3_xJalpha_re_%s_%d", rKey.c_str(), ib).Data());
          if (h2)
          {
            vector<string> lines = {
              "JES3 (RECO): x_{J} vs #alpha",
              TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data(),
              TString::Format("p_{T}^{#gamma}: %s", ptLab.c_str()).Data()
            };
            DrawAndSaveTH2_Common(ds, h2,
              JoinPath(D.dir2D, TString::Format("jes3_reco_xJ_vs_alpha_pTbin%d.png", ib).Data()),
              "x_{J}", "#alpha", "Counts", lines, true);
            delete h2;
          }
        }

        if (H.hTrut_xJ)
        {
          TH2* h2 = ProjectYZ_AtXbin_TH3(H.hTrut_xJ, ib, TString::Format("jes3_xJalpha_tr_%s_%d", rKey.c_str(), ib).Data());
          if (h2)
          {
            vector<string> lines = {
              "JES3 (TRUTH): x_{J} vs #alpha",
              TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data(),
              TString::Format("p_{T}^{#gamma}: %s", ptLab.c_str()).Data()
            };
            DrawAndSaveTH2_Common(ds, h2,
              JoinPath(D.dir2D, TString::Format("jes3_truth_xJ_vs_alpha_pTbin%d.png", ib).Data()),
              "x_{J}", "#alpha", "Counts", lines, true);
            delete h2;
          }
        }

        if (H.hReco_j1)
        {
          TH2* h2 = ProjectYZ_AtXbin_TH3(H.hReco_j1, ib, TString::Format("jes3_j1alpha_re_%s_%d", rKey.c_str(), ib).Data());
          if (h2)
          {
            vector<string> lines = {
              "JES3 (RECO): p_{T}^{jet1} vs #alpha",
              TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data(),
              TString::Format("p_{T}^{#gamma}: %s", ptLab.c_str()).Data()
            };
            DrawAndSaveTH2_Common(ds, h2,
              JoinPath(D.dir2D, TString::Format("jes3_reco_jet1Pt_vs_alpha_pTbin%d.png", ib).Data()),
              "p_{T}^{jet1} [GeV]", "#alpha", "Counts", lines, true);
            delete h2;
          }
        }

        if (H.hTrut_j1)
        {
          TH2* h2 = ProjectYZ_AtXbin_TH3(H.hTrut_j1, ib, TString::Format("jes3_j1alpha_tr_%s_%d", rKey.c_str(), ib).Data());
          if (h2)
          {
            vector<string> lines = {
              "JES3 (TRUTH): p_{T}^{jet1} vs #alpha",
              TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data(),
              TString::Format("p_{T}^{#gamma}: %s", ptLab.c_str()).Data()
            };
            DrawAndSaveTH2_Common(ds, h2,
              JoinPath(D.dir2D, TString::Format("jes3_truth_jet1Pt_vs_alpha_pTbin%d.png", ib).Data()),
              "p_{T}^{jet1} [GeV]", "#alpha", "Counts", lines, true);
            delete h2;
          }
        }
      };

      // For each pT bin: terminal table + summary + plots
      for (int ib = 1; ib <= nPt; ++ib)
      {
        BinPack P = BuildBinPack(ib);

        cout << std::left << std::setw(wPt) << P.ptLab
             << std::right
             << std::setw(wN) << std::fixed << std::setprecision(0) << P.nRe
             << std::setw(wM) << std::fixed << std::setprecision(4) << P.mxRe
             << std::setw(wM) << std::fixed << std::setprecision(4) << P.maRe
             << std::setw(wN) << std::fixed << std::setprecision(0) << P.nTr
             << std::setw(wM) << std::fixed << std::setprecision(4) << P.mxTr
             << std::setw(wM) << std::fixed << std::setprecision(4) << P.maTr
             << "\n";

        sumLines.push_back(TString::Format(
          "pTgamma=%s  Nre=%.0f  <xJ>re=%.6f  <a>re=%.6f  Ntr=%.0f  <xJ>tr=%.6f  <a>tr=%.6f",
          P.ptLab.c_str(), P.nRe, P.mxRe, P.maRe, P.nTr, P.mxTr, P.maTr
        ).Data());

        // Save xJ distributions (integrated over alpha) as individual PNGs
        SaveXJRecoPNGs(P.xJ_re, ib, P.ptLab);

        // (5) TRUTH reco-conditioned (jet-matched): h_JES3Truth_pT_xJ_alpha_<rKey>
        SaveXJTruthPNGs(P.xJ_tr, ib, P.ptLab);

        // (4) TRUTH reco-conditioned (NO jet match): h_JES3TruthRecoCondNoJetMatch_pT_xJ_alpha_<rKey>
        SaveXJTruthNoJetMatchPNGs(ib);

        // (3) TRUTH pure: h_JES3TruthPure_pT_xJ_alpha_<rKey>
        SaveXJTruthPurePNGs(ib);

        // (7) RECO truth-PHOTON-tagged
        SaveXJRecoTruthPhoTaggedPNGs(ib);

        // (8) RECO truth-tagged PHOTON+JET
        SaveXJRecoTruthTaggedPNGs(ib);

        // RECO/TRUTH xJ-alpha maps and jet1Pt-alpha maps
        SaveJes3Maps2D_ForBin(ib, P.ptLab);

        CleanupBinPack(P);
      }

        // -------------------------------------------------------------------------
        // overlayedWithSim: 2x3 table (DATA reco vs SIM reco) for integrated alpha
        // -------------------------------------------------------------------------
        if (isSimAndDataPP && !ds.isSim)
        {
          const string dirOv = JoinPath(D.dirXJProjReco, "insituCalib");
          EnsureDir(dirOv);

          const string dirFits = JoinPath(dirOv, "GaussianFits");
          EnsureDir(dirFits);

          static std::string s_lastSimPath_tbl = "";
          static TFile* s_fSim_tbl = nullptr;
          static TDirectory* s_simTopDir_tbl = nullptr;

          const std::string simPath = SimInputPathForSample(CurrentSimSample());
          if (simPath.empty())
          {
            cout << ANSI_BOLD_RED
                 << "[ERROR] [JES3 overlay table] SimInputPathForSample(CurrentSimSample()) returned EMPTY."
                 << "  rKey=" << rKey
                 << "  hist=h_JES3_pT_xJ_alpha_" << rKey
                 << ANSI_RESET << "\n";
          }
          else
          {
            if (!s_fSim_tbl || s_lastSimPath_tbl != simPath)
            {
              if (s_fSim_tbl) { s_fSim_tbl->Close(); delete s_fSim_tbl; s_fSim_tbl = nullptr; }
              s_simTopDir_tbl = nullptr;

              s_fSim_tbl = TFile::Open(simPath.c_str(), "READ");
              s_lastSimPath_tbl = simPath;

              if (!s_fSim_tbl || s_fSim_tbl->IsZombie())
              {
                cout << ANSI_BOLD_RED
                     << "[ERROR] [JES3 overlay table] Failed to open SIM file:"
                     << "  " << simPath
                     << "  (null or zombie)"
                     << ANSI_RESET << "\n";
              }
              else
              {
                s_simTopDir_tbl = s_fSim_tbl->GetDirectory(kDirSIM.c_str());
                if (!s_simTopDir_tbl)
                {
                  cout << ANSI_BOLD_RED
                       << "[ERROR] [JES3 overlay table] SIM topDir '" << kDirSIM << "' not found in file:"
                       << "  " << simPath
                       << "  -> falling back to file root"
                       << ANSI_RESET << "\n";
                  s_simTopDir_tbl = s_fSim_tbl;
                }
              }
            }
          }

          TH3* hSim3 = nullptr;
          if (s_simTopDir_tbl)
          {
            hSim3 = dynamic_cast<TH3*>(s_simTopDir_tbl->Get(("h_JES3_pT_xJ_alpha_" + rKey).c_str()));
            if (!hSim3)
            {
              cout << ANSI_BOLD_RED
                   << "[ERROR] [JES3 overlay table] Missing SIM TH3 in topDir '" << kDirSIM << "':"
                   << "  h_JES3_pT_xJ_alpha_" << rKey
                   << "  file=" << s_lastSimPath_tbl
                   << ANSI_RESET << "\n";
            }
          }
          else
          {
            cout << ANSI_BOLD_RED
                 << "[ERROR] [JES3 overlay table] s_simTopDir_tbl is NULL."
                 << "  file=" << s_lastSimPath_tbl
                 << ANSI_RESET << "\n";
          }

            if (hSim3 && H.hReco_xJ)
            {
                const int nCols = 3;
                const int nRows = 3;
                const int perPage = nCols * nRows;

                struct InSituPtGroup
                {
                  double ptLo = 0.0;
                  double ptHi = 0.0;
                  string label;
                  string tag;
                };

                const std::vector<InSituPtGroup> inSituPtGroups = {
                  {10.0, 12.0, "10 - 12", "10_12"},
                  {12.0, 14.0, "12 - 14", "12_14"},
                  {14.0, 16.0, "14 - 16", "14_16"},
                  {16.0, 18.0, "16 - 18", "16_18"},
                  {18.0, 20.0, "18 - 20", "18_20"},
                  {20.0, 24.0, "20 - 24", "20_24"},
                  {24.0, 35.0, "24 - 35", "24_35"}
                };

                const int nCalibBins = (int)inSituPtGroups.size();

                auto ProjectGroupedY_IntegratedAlpha_TH3 =
                  [&](const TH3* h3,
                      double ptLo,
                      double ptHi,
                      const std::string& newName)->TH1*
                {
                  if (!h3) return nullptr;

                  int xbinLo = -1;
                  int xbinHi = -1;
                  std::vector<double> w;
                  if (!XaxisOverlapWeights(h3->GetXaxis(), ptLo, ptHi, xbinLo, xbinHi, w)) return nullptr;

                  TH1* sum = nullptr;
                  const double alphaMax = h3->GetZaxis()->GetXmax();

                  for (int xb = xbinLo; xb <= xbinHi; ++xb)
                  {
                    TH1* h1 = ProjectY_AtXbin_AndAlphaMax_TH3(
                      h3, xb, alphaMax,
                      newName + TString::Format("_xb%d", xb).Data()
                    );
                    if (!h1) continue;

                    const int iw = xb - xbinLo;
                    const double ww = (iw >= 0 && iw < (int)w.size()) ? w[iw] : 1.0;

                    if (ww <= 0.0)
                    {
                      delete h1;
                      continue;
                    }

                    h1->Scale(ww);

                    if (!sum)
                    {
                      sum = CloneTH1(h1, newName);
                      if (sum)
                      {
                        sum->Reset("ICES");
                        sum->SetDirectory(nullptr);
                      }
                    }

                    if (sum) sum->Add(h1);
                    delete h1;
                  }

                  return sum;
                };

                auto PaintPadBlack = [&]()
                {
                  if (!gPad) return;
                  gPad->Clear();
                  gPad->SetFillColor(kBlack);
                  gPad->SetFrameFillColor(kBlack);
                  gPad->SetFillStyle(1001);
                  gPad->Range(0.0, 0.0, 1.0, 1.0);
                };

                // ---------------------------------------------------------------------
                // Table 1: overlays only (3x3 layout; last two pads black)
                // ---------------------------------------------------------------------
                TCanvas canTbl(
                  TString::Format("c_tbl_%s_dataVsSim", rKey.c_str()).Data(),
                  "c_tbl_dataVsSim", 1500, 900
                );
                canTbl.Divide(nCols, nRows, 0.001, 0.001);

                std::vector<TH1*> keep;
                keep.reserve(2 * nCalibBins);

                for (int k = 0; k < perPage; ++k)
                {
                  canTbl.cd(k + 1);

                  if (k >= nCalibBins)
                  {
                    PaintPadBlack();
                    continue;
                  }

                  const auto& G = inSituPtGroups[k];

                  gPad->SetLeftMargin(0.14);
                  gPad->SetRightMargin(0.05);
                  gPad->SetTopMargin(0.12);
                  gPad->SetBottomMargin(0.14);

                  const double ptMinGamma = G.ptLo;
                  const double ptMaxGamma = G.ptHi;

                  TH1* hDatRaw = ProjectGroupedY_IntegratedAlpha_TH3(
                    H.hReco_xJ, G.ptLo, G.ptHi,
                    TString::Format("h_tbl_dat_%s_%s", rKey.c_str(), G.tag.c_str()).Data()
                  );
                  TH1* hSimRaw = ProjectGroupedY_IntegratedAlpha_TH3(
                    hSim3, G.ptLo, G.ptHi,
                    TString::Format("h_tbl_sim_%s_%s", rKey.c_str(), G.tag.c_str()).Data()
                  );

                  if (!hDatRaw || !hSimRaw)
                  {
                    if (hDatRaw) delete hDatRaw;
                    if (hSimRaw) delete hSimRaw;
                    continue;
                  }

                  hDatRaw->SetDirectory(nullptr);
                  hSimRaw->SetDirectory(nullptr);

                  EnsureSumw2(hDatRaw);
                  EnsureSumw2(hSimRaw);

                  const double iDat = hDatRaw->Integral(0, hDatRaw->GetNbinsX() + 1);
                  const double iSim = hSimRaw->Integral(0, hSimRaw->GetNbinsX() + 1);
                  if (iDat > 0.0) hDatRaw->Scale(1.0 / iDat);
                  if (iSim > 0.0) hSimRaw->Scale(1.0 / iSim);

                  hDatRaw->SetTitle("");
                  hDatRaw->SetLineWidth(2);
                  hDatRaw->SetLineColor(kGreen + 2);
                  hDatRaw->SetMarkerStyle(20);
                  hDatRaw->SetMarkerSize(1.0);
                  hDatRaw->SetMarkerColor(kGreen + 2);

                  hSimRaw->SetLineWidth(2);
                  hSimRaw->SetLineColor(kOrange + 7);
                  hSimRaw->SetMarkerStyle(20);
                  hSimRaw->SetMarkerSize(1.0);
                  hSimRaw->SetMarkerColor(kOrange + 7);

                  hDatRaw->GetXaxis()->SetTitle("x_{J#gamma}");
                  hDatRaw->GetXaxis()->SetRangeUser(0.0, 2.0);
                  hDatRaw->GetYaxis()->SetTitle("Normalized counts");

                  hDatRaw->Draw("E1");
                  hSimRaw->Draw("E1 same");
                  gPad->Update();

                  const double jetPtMin_GeV = static_cast<double>(kJetPtMin);
                  const string bbLabel = B2BLabel();

                  TLegend* leg = new TLegend(0.70, 0.75, 0.94, 0.88);
                  leg->SetBorderSize(0);
                  leg->SetFillStyle(0);
                  leg->SetTextFont(42);
                  leg->SetTextSize(0.04);
                  leg->AddEntry(hDatRaw, "DATA (reco)", "ep");
                  leg->AddEntry(hSimRaw, "SIM (reco)",  "ep");
                  leg->DrawClone();
                  delete leg;

                  {
                      TLatex tCuts;
                      tCuts.SetNDC(true);
                      tCuts.SetTextFont(42);
                      tCuts.SetTextAlign(33);
                      tCuts.SetTextSize(0.04);
                      tCuts.DrawLatex(0.92, 0.62, TString::Format("|#Delta#phi(#gamma,jet)| > %s", bbLabel.c_str()).Data());
                      tCuts.DrawLatex(0.92, 0.54, TString::Format("p_{T}^{jet} > %.0f GeV", jetPtMin_GeV).Data());
                      tCuts.DrawLatex(0.92, 0.46, TString::Format("|v_{z}| < %.0f cm", vzCutCm).Data());
                  }

                  TLatex ttl;
                  ttl.SetNDC(true);
                  ttl.SetTextFont(42);
                  ttl.SetTextSize(0.043);
                  ttl.DrawLatex(0.12, 0.94,
                      TString::Format("RECO x_{J#gamma} (DATA vs SIM), p_{T}^{#gamma} = %.0f - %.0f GeV, R = %.1f",
                        ptMinGamma, ptMaxGamma, R).Data());

                  keep.push_back(hDatRaw);
                  keep.push_back(hSimRaw);
                }

                SaveCanvas(canTbl, JoinPath(dirOv, "table3x2_overlay_integratedAlpha_overlayedWithSim.png"));

                for (auto* h1 : keep) delete h1;

                // ---------------------------------------------------------------------
                // Table 2: overlays + iterative Gaussian fits + fit text
                // plus mean vs pT plot from Gaussian mean
                // ---------------------------------------------------------------------
                TCanvas canTblFits(
                    TString::Format("c_tbl_%s_dataVsSim_withFits", rKey.c_str()).Data(),
                    "c_tbl_dataVsSim_withFits", 1500, 900
                );
                canTblFits.Divide(nCols, nRows, 0.001, 0.001);

                TCanvas canTblMeans(
                    TString::Format("c_tbl_%s_dataVsSim_withMeans", rKey.c_str()).Data(),
                    "c_tbl_dataVsSim_withMeans", 1500, 900
                );
                canTblMeans.Divide(nCols, nRows, 0.001, 0.001);

                std::vector<TH1*> keepFitsH;
                keepFitsH.reserve(2 * nCalibBins);

                std::vector<TH1*> keepMeansH;
                keepMeansH.reserve(2 * nCalibBins);

                std::vector<TF1*> keepFitFns;
                keepFitFns.reserve(2 * nCalibBins);

                std::vector<double> vPtCtr, vPtErr;
                std::vector<double> vMuDat, vMuDatErr;
                std::vector<double> vMuSim, vMuSimErr;
                std::vector<double> vMeanDat, vMeanDatErr;
                std::vector<double> vMeanSim, vMeanSimErr;
                std::vector<double> vSigDat, vSigDatErr;
                std::vector<double> vSigSim, vSigSimErr;
                std::vector<double> vChi2NdfDat;
                std::vector<double> vChi2NdfSim;

                vPtCtr.reserve(nCalibBins);
                vPtErr.reserve(nCalibBins);
                vMuDat.reserve(nCalibBins);
                vMuDatErr.reserve(nCalibBins);
                vMuSim.reserve(nCalibBins);
                vMuSimErr.reserve(nCalibBins);
                vMeanDat.reserve(nCalibBins);
                vMeanDatErr.reserve(nCalibBins);
                vMeanSim.reserve(nCalibBins);
                vMeanSimErr.reserve(nCalibBins);
                vSigDat.reserve(nCalibBins);
                vSigDatErr.reserve(nCalibBins);
                vSigSim.reserve(nCalibBins);
                vSigSimErr.reserve(nCalibBins);
                vChi2NdfDat.reserve(nCalibBins);
                vChi2NdfSim.reserve(nCalibBins);

                auto FitIterGaus = [&](TH1* h, const std::string& fname, int lcolor) -> TF1*
                {
                      if (!h) return nullptr;
                      if (h->GetNbinsX() <= 0) return nullptr;
                      if (h->Integral(0, h->GetNbinsX() + 1) <= 0.0) return nullptr;

                      const TAxis* ax = h->GetXaxis();
                      if (!ax) return nullptr;

                      const double xLoHard = ax->GetXmin();
                      const double xHiHard = ax->GetXmax();

                      int firstNZ = 1;
                      while (firstNZ <= h->GetNbinsX() && h->GetBinContent(firstNZ) <= 0.0) ++firstNZ;

                      int lastNZ = h->GetNbinsX();
                      while (lastNZ >= 1 && h->GetBinContent(lastNZ) <= 0.0) --lastNZ;

                      if (firstNZ > lastNZ) return nullptr;

                      const double supportLo = h->GetBinLowEdge(firstNZ);
                      const double supportHi = h->GetBinLowEdge(lastNZ) + h->GetBinWidth(lastNZ);

                      int nActiveBins = 0;
                      for (int ib = firstNZ; ib <= lastNZ; ++ib)
                      {
                        if (h->GetBinContent(ib) > 0.0) ++nActiveBins;
                      }

                      const double nEff = h->GetEffectiveEntries();
                      const bool sparseHist     = (nEff < 80.0 || nActiveBins < 16);
                      const bool verySparseHist = (nEff < 35.0 || nActiveBins < 11);

                      auto SmoothedContent = [&](int ib) -> double
                      {
                        if (ib < firstNZ) ib = firstNZ;
                        if (ib > lastNZ)  ib = lastNZ;

                        double sum  = 0.0;
                        double wsum = 0.0;

                        for (int jb = ib - 1; jb <= ib + 1; ++jb)
                        {
                          if (jb < firstNZ || jb > lastNZ) continue;
                          const double w = (jb == ib ? 0.50 : 0.25);
                          sum  += w * std::max(0.0, h->GetBinContent(jb));
                          wsum += w;
                        }

                        return (wsum > 0.0 ? sum / wsum : std::max(0.0, h->GetBinContent(ib)));
                      };

                      int peakBin = firstNZ;
                      double peakY = -1.0;
                      for (int ib = firstNZ; ib <= lastNZ; ++ib)
                      {
                        const double s = SmoothedContent(ib);
                        if (s > peakY)
                        {
                          peakY   = s;
                          peakBin = ib;
                        }
                      }

                      if (peakBin < 1 || peakBin > h->GetNbinsX()) return nullptr;
                      if (!std::isfinite(peakY) || peakY <= 0.0) return nullptr;

                      const double binW = h->GetBinWidth(peakBin);

                      int plateauLo = peakBin;
                      int plateauHi = peakBin;
                      const double plateauFrac = (verySparseHist ? 0.90 : (sparseHist ? 0.92 : 0.94));

                      while (plateauLo > firstNZ && SmoothedContent(plateauLo - 1) >= plateauFrac * peakY) --plateauLo;
                      while (plateauHi < lastNZ  && SmoothedContent(plateauHi + 1) >= plateauFrac * peakY) ++plateauHi;

                      peakBin = plateauLo + (plateauHi - plateauLo) / 2;

                      double peakX = 0.5 * (h->GetBinCenter(plateauLo) + h->GetBinCenter(plateauHi));
                      {
                        double sw  = 0.0;
                        double swx = 0.0;

                        for (int jb = std::max(firstNZ, plateauLo - 1); jb <= std::min(lastNZ, plateauHi + 1); ++jb)
                        {
                          const double y = SmoothedContent(jb);
                          double w = std::max(0.0, y - plateauFrac * peakY);
                          if (w <= 0.0) w = ((plateauLo <= firstNZ + 1) ? 0.06 : 0.10) * y;
                          if (w <= 0.0) continue;

                          sw  += w;
                          swx += w * h->GetBinCenter(jb);
                        }

                        if (sw > 0.0) peakX = swx / sw;
                      }

                      if (peakBin > firstNZ && peakBin < lastNZ)
                      {
                        const double yL = SmoothedContent(peakBin - 1);
                        const double yC = SmoothedContent(peakBin);
                        const double yR = SmoothedContent(peakBin + 1);

                        const double denom = (yL - 2.0 * yC + yR);
                        if (std::isfinite(denom) && std::fabs(denom) > 1.0e-12)
                        {
                          double delta = 0.5 * (yL - yR) / denom;
                          if (delta < -0.40) delta = -0.40;
                          if (delta >  0.40) delta =  0.40;

                          const double parPeakX = h->GetBinCenter(peakBin) + delta * binW;
                          peakX = 0.75 * peakX + 0.25 * parPeakX;
                        }
                      }

                      auto FindCrossX = [&](double frac, bool goLeft, double& xCross) -> bool
                      {
                        const double thr = frac * peakY;

                        if (goLeft)
                        {
                          int ib = peakBin;
                          while (ib > firstNZ && SmoothedContent(ib) > thr) --ib;
                          if (SmoothedContent(ib) > thr) return false;
                          if (ib >= peakBin)
                          {
                            xCross = h->GetBinCenter(peakBin);
                            return true;
                          }

                          const double x1 = h->GetBinCenter(ib);
                          const double x2 = h->GetBinCenter(ib + 1);
                          const double y1 = SmoothedContent(ib);
                          const double y2 = SmoothedContent(ib + 1);

                          if (!std::isfinite(y1) || !std::isfinite(y2) || std::fabs(y2 - y1) < 1.0e-12)
                          {
                            xCross = x2;
                            return true;
                          }

                          const double t = (thr - y1) / (y2 - y1);
                          xCross = x1 + t * (x2 - x1);
                          return std::isfinite(xCross);
                        }

                        int ib = peakBin;
                        while (ib < lastNZ && SmoothedContent(ib) > thr) ++ib;
                        if (SmoothedContent(ib) > thr) return false;
                        if (ib <= peakBin)
                        {
                          xCross = h->GetBinCenter(peakBin);
                          return true;
                        }

                        const double x1 = h->GetBinCenter(ib - 1);
                        const double x2 = h->GetBinCenter(ib);
                        const double y1 = SmoothedContent(ib - 1);
                        const double y2 = SmoothedContent(ib);

                        if (!std::isfinite(y1) || !std::isfinite(y2) || std::fabs(y2 - y1) < 1.0e-12)
                        {
                          xCross = x1;
                          return true;
                        }

                        const double t = (thr - y1) / (y2 - y1);
                        xCross = x1 + t * (x2 - x1);
                        return std::isfinite(xCross);
                      };

                      double xL70 = 0.0, xL55 = 0.0;
                      double xR70 = 0.0, xR55 = 0.0, xR40 = 0.0;

                      const bool hasL70 = FindCrossX(0.70, true,  xL70);
                      const bool hasL55 = FindCrossX(0.55, true,  xL55);
                      const bool hasR70 = FindCrossX(0.70, false, xR70);
                      const bool hasR55 = FindCrossX(0.55, false, xR55);
                      const bool hasR40 = FindCrossX(0.40, false, xR40);

                      auto SigmaFromOneSide = [&](double xCross, double frac) -> double
                      {
                        if (!std::isfinite(xCross) || xCross <= peakX) return 0.0;
                        const double scale = std::sqrt(std::max(1.0e-12, -2.0 * std::log(frac)));
                        if (!std::isfinite(scale) || scale <= 0.0) return 0.0;
                        return (xCross - peakX) / scale;
                      };

                      auto SigmaFromTwoSides = [&](double xLeft, double xRight, double frac) -> double
                      {
                        if (!std::isfinite(xLeft) || !std::isfinite(xRight) || xRight <= xLeft) return 0.0;
                        const double scale = std::sqrt(std::max(1.0e-12, -2.0 * std::log(frac)));
                        if (!std::isfinite(scale) || scale <= 0.0) return 0.0;
                        return 0.5 * (xRight - xLeft) / scale;
                      };

                      double sigNum = 0.0;
                      double sigDen = 0.0;

                      auto AddSigma = [&](double s, double w)
                      {
                        if (!std::isfinite(s) || s <= 0.0 || !std::isfinite(w) || w <= 0.0) return;
                        sigNum += w * s;
                        sigDen += w;
                      };

                      if (hasL70 && hasR70) AddSigma(SigmaFromTwoSides(xL70, xR70, 0.70), 2.0);
                      if (hasL55 && hasR55) AddSigma(SigmaFromTwoSides(xL55, xR55, 0.55), 1.5);

                      if (sigDen <= 0.0)
                      {
                        if (hasR70) AddSigma(SigmaFromOneSide(xR70, 0.70), 2.0);
                        if (hasR55) AddSigma(SigmaFromOneSide(xR55, 0.55), 1.5);
                        if (hasR40) AddSigma(SigmaFromOneSide(xR40, 0.40), 1.0);
                      }

                      double sig = (sigDen > 0.0 ? sigNum / sigDen : 0.0);

                      if (!std::isfinite(sig) || sig <= 0.0)
                      {
                        double sw   = 0.0;
                        double swx  = 0.0;
                        double swx2 = 0.0;

                        const int lo = std::max(firstNZ, peakBin - 1);
                        const int hi = std::min(lastNZ,  peakBin + (sparseHist ? 2 : 1));

                        for (int ib = lo; ib <= hi; ++ib)
                        {
                          const double w = SmoothedContent(ib);
                          if (w <= 0.0) continue;

                          const double x = h->GetBinCenter(ib);
                          sw   += w;
                          swx  += w * x;
                          swx2 += w * x * x;
                        }

                        if (sw > 0.0)
                        {
                          const double muTmp  = swx / sw;
                          const double varTmp = swx2 / sw - muTmp * muTmp;
                          if (varTmp > 0.0) sig = std::sqrt(varTmp);
                        }
                      }

                      if (!std::isfinite(sig) || sig <= 0.0) sig = 2.0 * binW;

                      const double visibleSpan = std::max(3.0 * binW, supportHi - supportLo);
                      const double sigFloor = std::max(1.00 * binW, (verySparseHist ? 0.050 : (sparseHist ? 0.040 : 0.030)) * visibleSpan);
                      if (sig < sigFloor) sig = sigFloor;

                      const double firstNZCenter = h->GetBinCenter(firstNZ);
                      const bool sharpTurnOn = (!(hasL70 && hasL55) || plateauLo <= firstNZ + 1 || ((peakX - firstNZCenter) < 2.5 * binW));

                      double seedLo = sharpTurnOn
                                    ? (hasL70 ? xL70 : peakX - 0.75 * sig)
                                    : (hasL55 ? xL55 : peakX - 1.10 * sig);

                      double seedHi = hasR40
                                    ? xR40
                                    : (hasR55 ? peakX + 1.55 * sig : peakX + 1.80 * sig);

                      if (sparseHist)     seedHi += 0.25 * sig;
                      if (verySparseHist) seedHi += 0.35 * sig;

                      if (!std::isfinite(seedLo)) seedLo = peakX - 0.90 * sig;
                      if (!std::isfinite(seedHi)) seedHi = peakX + 1.80 * sig;

                      if (seedLo < supportLo) seedLo = supportLo;
                      if (seedHi > supportHi) seedHi = supportHi;

                      const double minSeedWidth = (verySparseHist ? 4.8 : (sparseHist ? 4.4 : 4.0)) * binW;
                      if (seedHi - seedLo < minSeedWidth)
                      {
                        const double leftFrac  = sharpTurnOn ? 0.35 : 0.45;
                        const double rightFrac = 1.0 - leftFrac;

                        seedLo = peakX - leftFrac  * minSeedWidth;
                        seedHi = peakX + rightFrac * minSeedWidth;

                        if (seedLo < supportLo)
                        {
                          seedHi += (supportLo - seedLo);
                          seedLo  = supportLo;
                        }
                        if (seedHi > supportHi)
                        {
                          seedLo -= (seedHi - supportHi);
                          seedHi  = supportHi;
                        }
                        if (seedLo < supportLo) seedLo = supportLo;
                      }

                      if (seedHi <= seedLo) return nullptr;

                      TF1* f = new TF1(fname.c_str(), "gaus", seedLo, seedHi);
                      f->SetLineColor(lcolor);
                      f->SetLineWidth(3);
                      f->SetLineStyle(1);
                      f->SetNpx(800);
                      f->SetParameters(peakY, peakX, sig);
                      f->SetParLimits(0, 0.0, 5.0 * std::max(peakY, h->GetMaximum()));

                      double muLoLimit = sharpTurnOn ? std::max(seedLo, peakX - 0.10 * sig)
                                                     : std::max(seedLo, peakX - 0.30 * sig);
                      double muHiLimit = sharpTurnOn ? std::min(seedHi, peakX + 0.30 * sig)
                                                     : std::min(seedHi, peakX + 0.55 * sig);

                      if (muHiLimit <= muLoLimit)
                      {
                        muLoLimit = seedLo;
                        muHiLimit = seedHi;
                      }

                      f->SetParLimits(1, muLoLimit, muHiLimit);
                      f->SetParLimits(2, 0.60 * sigFloor, 0.35 * (xHiHard - xLoHard));

                      int fitStatus = h->Fit(f, "RQ0NWLI");
                      if (fitStatus != 0) fitStatus = h->Fit(f, "RQ0NLI");

                      double muFit  = f->GetParameter(1);
                      double sigFit = std::fabs(f->GetParameter(2));

                      if (!std::isfinite(muFit) || muFit < seedLo || muFit > seedHi) muFit = peakX;
                      if (!std::isfinite(sigFit) || sigFit <= 0.0) sigFit = sig;
                      if (sigFit < sigFloor) sigFit = sigFloor;

                      double finalLo = sharpTurnOn ? std::max(seedLo, muFit - 0.75 * sigFit)
                                                   : std::max(seedLo, muFit - 1.00 * sigFit);
                      double finalHi = std::min(
                        supportHi,
                        std::max(seedHi, muFit + (verySparseHist ? 2.20 * sigFit : (sparseHist ? 1.95 * sigFit : 1.70 * sigFit)))
                      );

                      const double minFinalWidth = (verySparseHist ? 4.4 : (sparseHist ? 4.0 : 3.6)) * binW;
                      if (finalHi - finalLo < minFinalWidth)
                      {
                        const double leftFrac  = sharpTurnOn ? 0.35 : 0.45;
                        const double rightFrac = 1.0 - leftFrac;

                        finalLo = muFit - leftFrac  * minFinalWidth;
                        finalHi = muFit + rightFrac * minFinalWidth;

                        if (finalLo < supportLo)
                        {
                          finalHi += (supportLo - finalLo);
                          finalLo  = supportLo;
                        }
                        if (finalHi > supportHi)
                        {
                          finalLo -= (finalHi - supportHi);
                          finalHi  = supportHi;
                        }
                        if (finalLo < supportLo) finalLo = supportLo;
                      }

                      if (finalLo < seedLo) finalLo = seedLo;
                      if (finalHi > supportHi) finalHi = supportHi;
                      if (finalHi <= finalLo)
                      {
                        finalLo = seedLo;
                        finalHi = seedHi;
                      }

                      f->SetRange(finalLo, finalHi);
                      f->SetParameter(1, muFit);
                      f->SetParameter(2, sigFit);

                      fitStatus = h->Fit(f, "RQ0NWLI");
                      if (fitStatus != 0) fitStatus = h->Fit(f, "RQ0NLI");

                      double fitMu  = f->GetParameter(1);
                      double fitSig = std::fabs(f->GetParameter(2));

                      if (!std::isfinite(fitMu) || !std::isfinite(fitSig) || fitSig <= 0.0)
                      {
                        delete f;
                        return nullptr;
                      }

                      if (fitSig < sigFloor) fitSig = sigFloor;

                      if (sparseHist)
                      {
                        const double refitLo = finalLo;
                        const double refitHi = std::min(
                          supportHi,
                          std::max(finalHi, fitMu + (verySparseHist ? 2.35 * fitSig : 2.05 * fitSig))
                        );

                        if (refitHi > refitLo + 0.5 * binW)
                        {
                          f->SetRange(refitLo, refitHi);
                          f->SetParameter(1, fitMu);
                          f->SetParameter(2, fitSig);

                          fitStatus = h->Fit(f, "RQ0NWLI");
                          if (fitStatus != 0) fitStatus = h->Fit(f, "RQ0NLI");

                          const double refitMu  = f->GetParameter(1);
                          const double refitSig = std::fabs(f->GetParameter(2));

                          if (std::isfinite(refitMu) && std::isfinite(refitSig) && refitSig > 0.0)
                          {
                            fitMu  = refitMu;
                            fitSig = refitSig;
                            finalHi = refitHi;
                          }
                        }
                      }

                      if (!std::isfinite(fitMu) || !std::isfinite(fitSig) || fitSig <= 0.0)
                      {
                        delete f;
                        return nullptr;
                      }

                      f->SetRange(finalLo, finalHi);
                      f->SetLineStyle(1);
                      f->SetLineColor(lcolor);
                      f->SetLineWidth(3);
                      f->SetNpx(800);
                      return f;
                };
                
                auto DrawFitOnly = [&](TF1* f)
                {
                      if (!f) return;

                      f->SetLineStyle(1);
                      f->SetLineWidth(3);
                      f->SetNpx(800);
                      f->Draw("same");
                };
                for (int k = 0; k < perPage; ++k)
                {
                  canTblFits.cd(k + 1);

                  if (k >= nCalibBins)
                  {
                    PaintPadBlack();
                    canTblMeans.cd(k + 1);
                    PaintPadBlack();
                    continue;
                  }

                  const auto& G = inSituPtGroups[k];

                  const double ptMinGamma = G.ptLo;
                  const double ptMaxGamma = G.ptHi;
                  const double ptCtr = 0.5 * (ptMinGamma + ptMaxGamma);
                  const double ptErr = 0.5 * (ptMaxGamma - ptMinGamma);

                  TH1* hDatRaw = ProjectGroupedY_IntegratedAlpha_TH3(
                    H.hReco_xJ, G.ptLo, G.ptHi,
                    TString::Format("h_tbl_fit_dat_%s_%s", rKey.c_str(), G.tag.c_str()).Data()
                  );
                  TH1* hSimRaw = ProjectGroupedY_IntegratedAlpha_TH3(
                    hSim3, G.ptLo, G.ptHi,
                    TString::Format("h_tbl_fit_sim_%s_%s", rKey.c_str(), G.tag.c_str()).Data()
                  );

                  if (!hDatRaw || !hSimRaw)
                  {
                    if (hDatRaw) delete hDatRaw;
                    if (hSimRaw) delete hSimRaw;
                    canTblFits.cd(k + 1);
                    PaintPadBlack();
                    canTblMeans.cd(k + 1);
                    PaintPadBlack();
                    continue;
                  }

                  hDatRaw->SetDirectory(nullptr);
                  hSimRaw->SetDirectory(nullptr);

                  EnsureSumw2(hDatRaw);
                  EnsureSumw2(hSimRaw);

                  const double iDat = hDatRaw->Integral(0, hDatRaw->GetNbinsX() + 1);
                  const double iSim = hSimRaw->Integral(0, hSimRaw->GetNbinsX() + 1);
                  if (iDat > 0.0) hDatRaw->Scale(1.0 / iDat);
                  if (iSim > 0.0) hSimRaw->Scale(1.0 / iSim);

                  TH1* hDatMeanTbl = CloneTH1(hDatRaw, TString::Format("h_tbl_mean_dat_%s_%s", rKey.c_str(), G.tag.c_str()).Data());
                  TH1* hSimMeanTbl = CloneTH1(hSimRaw, TString::Format("h_tbl_mean_sim_%s_%s", rKey.c_str(), G.tag.c_str()).Data());

                  if (hDatMeanTbl) hDatMeanTbl->SetDirectory(nullptr);
                  if (hSimMeanTbl) hSimMeanTbl->SetDirectory(nullptr);

                  const double jetPtMin_GeV = static_cast<double>(kJetPtMin);
                  const string bbLabel = B2BLabel();

                  canTblFits.cd(k + 1);

                  gPad->SetLeftMargin(0.14);
                  gPad->SetRightMargin(0.05);
                  gPad->SetTopMargin(0.12);
                  gPad->SetBottomMargin(0.14);

                  hDatRaw->SetTitle("");
                  hDatRaw->SetLineWidth(2);
                  hDatRaw->SetLineColor(kGreen + 2);
                  hDatRaw->SetMarkerStyle(20);
                  hDatRaw->SetMarkerSize(1.0);
                  hDatRaw->SetMarkerColor(kGreen + 2);

                  hSimRaw->SetLineWidth(2);
                  hSimRaw->SetLineColor(kOrange + 7);
                  hSimRaw->SetMarkerStyle(20);
                  hSimRaw->SetMarkerSize(1.0);
                  hSimRaw->SetMarkerColor(kOrange + 7);

                  hDatRaw->GetXaxis()->SetTitle("x_{J#gamma}");
                  hDatRaw->GetXaxis()->SetRangeUser(0.0, 2.0);
                  hDatRaw->GetYaxis()->SetTitle("Normalized counts");

                  hDatRaw->Draw("E1");
                  hSimRaw->Draw("E1 same");
                  gPad->Update();

                  TF1* fDat = FitIterGaus(hDatRaw, TString::Format("f_tbl_dat_%s_%s", rKey.c_str(), G.tag.c_str()).Data(), kGreen + 2);
                  TF1* fSim = FitIterGaus(hSimRaw, TString::Format("f_tbl_sim_%s_%s", rKey.c_str(), G.tag.c_str()).Data(), kOrange + 7);

                  if (fDat) { DrawFitOnly(fDat); keepFitFns.push_back(fDat); }
                  if (fSim) { DrawFitOnly(fSim); keepFitFns.push_back(fSim); }

                  {
                    TLegend* leg = new TLegend(0.70, 0.75, 0.94, 0.88);
                    leg->SetBorderSize(0);
                    leg->SetFillStyle(0);
                    leg->SetTextFont(42);
                    leg->SetTextSize(0.04);
                    leg->AddEntry(hDatRaw, "DATA (reco)", "ep");
                    leg->AddEntry(hSimRaw, "SIM (reco)",  "ep");
                    leg->DrawClone();
                    delete leg;
                  }

                  {
                    if (fDat && fDat->GetNDF() > 0)
                    {
                      const double mu   = fDat->GetParameter(1);
                      const double sig  = fDat->GetParameter(2);
                      const double chi2 = fDat->GetChisquare();
                      const double ndf  = fDat->GetNDF();

                      vMuDat.push_back(mu);
                      vMuDatErr.push_back(fDat->GetParError(1));

                      vSigDat.push_back(sig);
                      vSigDatErr.push_back(fDat->GetParError(2));

                      vChi2NdfDat.push_back(chi2 / ndf);
                    }
                    else
                    {
                      vMuDat.push_back(-1.0);
                      vMuDatErr.push_back(0.0);

                      vSigDat.push_back(-1.0);
                      vSigDatErr.push_back(0.0);

                      vChi2NdfDat.push_back(-1.0);
                    }

                    if (fSim && fSim->GetNDF() > 0)
                    {
                      const double mu   = fSim->GetParameter(1);
                      const double sig  = fSim->GetParameter(2);
                      const double chi2 = fSim->GetChisquare();
                      const double ndf  = fSim->GetNDF();

                      vMuSim.push_back(mu);
                      vMuSimErr.push_back(fSim->GetParError(1));

                      vSigSim.push_back(sig);
                      vSigSimErr.push_back(fSim->GetParError(2));

                      vChi2NdfSim.push_back(chi2 / ndf);
                    }
                    else
                    {
                      vMuSim.push_back(-1.0);
                      vMuSimErr.push_back(0.0);

                      vSigSim.push_back(-1.0);
                      vSigSimErr.push_back(0.0);

                      vChi2NdfSim.push_back(-1.0);
                    }
                  }

                  const double muDatDraw = (fDat && fDat->GetNDF() > 0) ? fDat->GetParameter(1) : -1.0;
                  const double muSimDraw = (fSim && fSim->GetNDF() > 0) ? fSim->GetParameter(1) : -1.0;

                  gPad->Update();
                  const double yMinFit = gPad->GetUymin();
                  const double yMaxFit = gPad->GetUymax();

                  if (muDatDraw >= 0.0)
                  {
                    TLine lDatFit(muDatDraw, yMinFit, muDatDraw, yMaxFit);
                    lDatFit.SetLineColor(kGreen + 2);
                    lDatFit.SetLineStyle(2);
                    lDatFit.SetLineWidth(2);
                    lDatFit.DrawClone("same");
                  }
                  if (muSimDraw >= 0.0)
                  {
                    TLine lSimFit(muSimDraw, yMinFit, muSimDraw, yMaxFit);
                    lSimFit.SetLineColor(kOrange + 7);
                    lSimFit.SetLineStyle(2);
                    lSimFit.SetLineWidth(2);
                    lSimFit.DrawClone("same");
                  }

                  {
                    TLatex tMean;
                    tMean.SetNDC(true);
                    tMean.SetTextFont(42);
                    tMean.SetTextAlign(33);
                    tMean.SetTextSize(0.040);
                    if (muDatDraw >= 0.0) tMean.DrawLatex(0.92, 0.70, TString::Format("Data <x_{J}> = %.4f", muDatDraw).Data());
                    else                  tMean.DrawLatex(0.92, 0.70, "Data <x_{J}> = N/A");
                    if (muSimDraw >= 0.0) tMean.DrawLatex(0.92, 0.64, TString::Format("Sim <x_{J}> = %.4f", muSimDraw).Data());
                    else                  tMean.DrawLatex(0.92, 0.64, "Sim <x_{J}> = N/A");
                  }

                  {
                    TLatex tCuts;
                    tCuts.SetNDC(true);
                    tCuts.SetTextFont(42);
                    tCuts.SetTextAlign(33);
                    tCuts.SetTextSize(0.04);
                    tCuts.DrawLatex(0.92, 0.56, TString::Format("|#Delta#phi(#gamma,jet)| > %s", bbLabel.c_str()).Data());
                    tCuts.DrawLatex(0.92, 0.50, TString::Format("p_{T}^{jet} > %.0f GeV", jetPtMin_GeV).Data());
                    tCuts.DrawLatex(0.92, 0.44, TString::Format("|v_{z}| < %.0f cm", vzCutCm).Data());
                  }

                  {
                    TLatex ttl;
                    ttl.SetNDC(true);
                    ttl.SetTextFont(42);
                    ttl.SetTextSize(0.043);
                    ttl.DrawLatex(0.12, 0.94,
                      TString::Format("RECO x_{J#gamma} (DATA vs SIM), p_{T}^{#gamma} = %.0f - %.0f GeV, R = %.1f",
                        ptMinGamma, ptMaxGamma, R).Data());
                  }

                  canTblMeans.cd(k + 1);

                  gPad->SetLeftMargin(0.14);
                  gPad->SetRightMargin(0.05);
                  gPad->SetTopMargin(0.12);
                  gPad->SetBottomMargin(0.14);

                  if (hDatMeanTbl && hSimMeanTbl)
                  {
                    hDatMeanTbl->SetTitle("");
                    hDatMeanTbl->SetLineWidth(2);
                    hDatMeanTbl->SetLineColor(kGreen + 2);
                    hDatMeanTbl->SetMarkerStyle(20);
                    hDatMeanTbl->SetMarkerSize(1.0);
                    hDatMeanTbl->SetMarkerColor(kGreen + 2);

                    hSimMeanTbl->SetLineWidth(2);
                    hSimMeanTbl->SetLineColor(kOrange + 7);
                    hSimMeanTbl->SetMarkerStyle(20);
                    hSimMeanTbl->SetMarkerSize(1.0);
                    hSimMeanTbl->SetMarkerColor(kOrange + 7);

                    hDatMeanTbl->GetXaxis()->SetTitle("x_{J#gamma}");
                    hDatMeanTbl->GetXaxis()->SetRangeUser(0.0, 2.0);
                    hDatMeanTbl->GetYaxis()->SetTitle("Normalized counts");

                    hDatMeanTbl->Draw("E1");
                    hSimMeanTbl->Draw("E1 same");
                    gPad->Update();

                    const double meanDatDraw = hDatMeanTbl->GetMean();
                    const double meanSimDraw = hSimMeanTbl->GetMean();

                    const double yMinMean = gPad->GetUymin();
                    const double yMaxMean = gPad->GetUymax();

                    if (std::isfinite(meanDatDraw))
                    {
                      TLine lDatMean(meanDatDraw, yMinMean, meanDatDraw, yMaxMean);
                      lDatMean.SetLineColor(kGreen + 2);
                      lDatMean.SetLineStyle(2);
                      lDatMean.SetLineWidth(2);
                      lDatMean.DrawClone("same");
                    }
                    if (std::isfinite(meanSimDraw))
                    {
                      TLine lSimMean(meanSimDraw, yMinMean, meanSimDraw, yMaxMean);
                      lSimMean.SetLineColor(kOrange + 7);
                      lSimMean.SetLineStyle(2);
                      lSimMean.SetLineWidth(2);
                      lSimMean.DrawClone("same");
                    }

                    {
                      TLegend* leg = new TLegend(0.70, 0.75, 0.94, 0.88);
                      leg->SetBorderSize(0);
                      leg->SetFillStyle(0);
                      leg->SetTextFont(42);
                      leg->SetTextSize(0.04);
                      leg->AddEntry(hDatMeanTbl, "DATA (reco)", "ep");
                      leg->AddEntry(hSimMeanTbl, "SIM (reco)",  "ep");
                      leg->DrawClone();
                      delete leg;
                    }

                    {
                      TLatex tMean;
                      tMean.SetNDC(true);
                      tMean.SetTextFont(42);
                      tMean.SetTextAlign(33);
                      tMean.SetTextSize(0.040);
                      if (std::isfinite(meanDatDraw)) tMean.DrawLatex(0.92, 0.70, TString::Format("Data <x_{J}> = %.4f", meanDatDraw).Data());
                      else                           tMean.DrawLatex(0.92, 0.70, "Data <x_{J}> = N/A");
                      if (std::isfinite(meanSimDraw)) tMean.DrawLatex(0.92, 0.64, TString::Format("Sim <x_{J}> = %.4f", meanSimDraw).Data());
                      else                           tMean.DrawLatex(0.92, 0.64, "Sim <x_{J}> = N/A");
                    }

                    {
                      TLatex tCuts;
                      tCuts.SetNDC(true);
                      tCuts.SetTextFont(42);
                      tCuts.SetTextAlign(33);
                      tCuts.SetTextSize(0.04);
                      tCuts.DrawLatex(0.92, 0.56, TString::Format("|#Delta#phi(#gamma,jet)| > %s", bbLabel.c_str()).Data());
                      tCuts.DrawLatex(0.92, 0.50, TString::Format("p_{T}^{jet} > %.0f GeV", jetPtMin_GeV).Data());
                      tCuts.DrawLatex(0.92, 0.44, TString::Format("|v_{z}| < %.0f cm", vzCutCm).Data());
                    }

                    {
                      TLatex ttl;
                      ttl.SetNDC(true);
                      ttl.SetTextFont(42);
                      ttl.SetTextSize(0.043);
                      ttl.DrawLatex(0.12, 0.94,
                        TString::Format("RECO x_{J#gamma} (DATA vs SIM), p_{T}^{#gamma} = %.0f - %.0f GeV, R = %.1f",
                          ptMinGamma, ptMaxGamma, R).Data());
                    }

                    keepMeansH.push_back(hDatMeanTbl);
                    keepMeansH.push_back(hSimMeanTbl);
                  }
                  else
                  {
                    if (hDatMeanTbl) delete hDatMeanTbl;
                    if (hSimMeanTbl) delete hSimMeanTbl;
                  }

                  vPtCtr.push_back(ptCtr);
                  vPtErr.push_back(ptErr);

                  keepFitsH.push_back(hDatRaw);
                  keepFitsH.push_back(hSimRaw);
                }

                SaveCanvas(canTblFits, JoinPath(dirFits, "table3x3_overlay_integratedAlpha_overlayedWithSim_withFits.png"));
                SaveCanvas(canTblMeans, JoinPath(dirFits, "table3x3_overlay_integratedAlpha_overlayedWithSim_withMeans.png"));

                for (auto* h1 : keepMeansH) delete h1;

              // Cache GetMean values for the already-drawn table bins.
              {
                  const int nTableStored = (int)vPtCtr.size();

                  vMeanDat.clear();
                  vMeanDatErr.clear();
                  vMeanSim.clear();
                  vMeanSimErr.clear();

                  for (int i = 0; i < nTableStored; ++i)
                  {
                    TH1* hDatMean = (2*i     < (int)keepFitsH.size()) ? keepFitsH[2*i]     : nullptr;
                    TH1* hSimMean = (2*i + 1 < (int)keepFitsH.size()) ? keepFitsH[2*i + 1] : nullptr;

                    if (hDatMean)
                    {
                      vMeanDat.push_back(hDatMean->GetMean());
                      vMeanDatErr.push_back(hDatMean->GetMeanError());
                    }
                    else
                    {
                      vMeanDat.push_back(-1.0);
                      vMeanDatErr.push_back(0.0);
                    }

                    if (hSimMean)
                    {
                      vMeanSim.push_back(hSimMean->GetMean());
                      vMeanSimErr.push_back(hSimMean->GetMeanError());
                    }
                    else
                    {
                      vMeanSim.push_back(-1.0);
                      vMeanSimErr.push_back(0.0);
                    }
                  }
                }

                auto SaveGaussianFitsPerPtBin =
                  [&](int ig)
                {
                  if (ig < 0 || ig >= nCalibBins) return;

                  const auto& G = inSituPtGroups[ig];
                  const double ptMinGamma = G.ptLo;
                  const double ptMaxGamma = G.ptHi;

                  TH1* hDatOne = ProjectGroupedY_IntegratedAlpha_TH3(
                    H.hReco_xJ, G.ptLo, G.ptHi,
                    TString::Format("h_gausFit_dat_%s_%s", rKey.c_str(), G.tag.c_str()).Data()
                  );
                  TH1* hSimOne = ProjectGroupedY_IntegratedAlpha_TH3(
                    hSim3, G.ptLo, G.ptHi,
                    TString::Format("h_gausFit_sim_%s_%s", rKey.c_str(), G.tag.c_str()).Data()
                  );

                  if (!hDatOne || !hSimOne)
                  {
                    if (hDatOne) delete hDatOne;
                    if (hSimOne) delete hSimOne;
                    return;
                  }

                  hDatOne->SetDirectory(nullptr);
                  hSimOne->SetDirectory(nullptr);

                  EnsureSumw2(hDatOne);
                  EnsureSumw2(hSimOne);

                  const double iDatOne = hDatOne->Integral(0, hDatOne->GetNbinsX() + 1);
                  const double iSimOne = hSimOne->Integral(0, hSimOne->GetNbinsX() + 1);
                  if (iDatOne > 0.0) hDatOne->Scale(1.0 / iDatOne);
                  if (iSimOne > 0.0) hSimOne->Scale(1.0 / iSimOne);

                  TF1* fDatOne = FitIterGaus(hDatOne, TString::Format("f_gausFit_dat_%s_%s", rKey.c_str(), G.tag.c_str()).Data(), kGreen + 2);
                  TF1* fSimOne = FitIterGaus(hSimOne, TString::Format("f_gausFit_sim_%s_%s", rKey.c_str(), G.tag.c_str()).Data(), kOrange + 7);

                  TCanvas cOne(
                    TString::Format("c_gausFit_%s_%s", rKey.c_str(), G.tag.c_str()).Data(),
                    "c_gausFit", 900, 700
                  );
                  ApplyCanvasMargins1D(cOne);
                  cOne.SetLogy(false);

                  hDatOne->SetTitle("");
                  hDatOne->SetLineWidth(2);
                  hDatOne->SetLineColor(kGreen + 2);
                  hDatOne->SetMarkerStyle(20);
                  hDatOne->SetMarkerSize(1.0);
                  hDatOne->SetMarkerColor(kGreen + 2);

                  hSimOne->SetLineWidth(2);
                  hSimOne->SetLineColor(kOrange + 7);
                  hSimOne->SetMarkerStyle(20);
                  hSimOne->SetMarkerSize(1.0);
                  hSimOne->SetMarkerColor(kOrange + 7);

                  hDatOne->GetXaxis()->SetTitle("x_{J#gamma}");
                  hDatOne->GetXaxis()->SetRangeUser(0.0, 2.0);
                  hDatOne->GetYaxis()->SetTitle("Normalized counts");

                  hDatOne->Draw("E1");
                  hSimOne->Draw("E1 same");
                  if (fDatOne) DrawFitOnly(fDatOne);
                  if (fSimOne) DrawFitOnly(fSimOne);

                  TLegend leg(0.70, 0.75, 0.94, 0.88);
                  leg.SetBorderSize(0);
                  leg.SetFillStyle(0);
                  leg.SetTextFont(42);
                  leg.SetTextSize(0.04);
                  leg.AddEntry(hDatOne, "DATA (reco)", "ep");
                  leg.AddEntry(hSimOne, "SIM (reco)",  "ep");
                  leg.Draw();

                  {
                    const double jetPtMin_GeV = static_cast<double>(kJetPtMin);
                    const string bbLabel = B2BLabel();

                    TLatex tCuts;
                    tCuts.SetNDC(true);
                    tCuts.SetTextFont(42);
                    tCuts.SetTextAlign(33);
                    tCuts.SetTextSize(0.04);
                    tCuts.DrawLatex(0.92, 0.62, TString::Format("|#Delta#phi(#gamma,jet)| > %s", bbLabel.c_str()).Data());
                    tCuts.DrawLatex(0.92, 0.54, TString::Format("p_{T}^{jet} > %.0f GeV", jetPtMin_GeV).Data());
                    tCuts.DrawLatex(0.92, 0.46, TString::Format("|v_{z}| < %.0f cm", vzCutCm).Data());
                  }

                  TLatex ttl;
                  ttl.SetNDC(true);
                  ttl.SetTextFont(42);
                  ttl.SetTextSize(0.043);
                  ttl.DrawLatex(0.12, 0.94,
                    TString::Format("RECO x_{J#gamma} (DATA vs SIM), p_{T}^{#gamma} = %.0f - %.0f GeV, R = %.1f",
                      ptMinGamma, ptMaxGamma, R).Data());

                  const double muDatDraw = (fDatOne && fDatOne->GetNDF() > 0) ? fDatOne->GetParameter(1) : -1.0;
                  const double muSimDraw = (fSimOne && fSimOne->GetNDF() > 0) ? fSimOne->GetParameter(1) : -1.0;

                  TLatex tMean;
                  tMean.SetNDC(true);
                  tMean.SetTextFont(42);
                  tMean.SetTextAlign(13);
                  tMean.SetTextSize(0.040);
                  if (muDatDraw >= 0.0) tMean.DrawLatex(0.16, 0.86, TString::Format("Data <x_{J}> = %.4f", muDatDraw).Data());
                  else                  tMean.DrawLatex(0.16, 0.86, "Data <x_{J}> = N/A");
                  if (muSimDraw >= 0.0) tMean.DrawLatex(0.16, 0.80, TString::Format("Sim <x_{J}> = %.4f", muSimDraw).Data());
                  else                  tMean.DrawLatex(0.16, 0.80, "Sim <x_{J}> = N/A");

                  SaveCanvas(cOne, JoinPath(dirFits,
                    TString::Format("xJ_reco_integratedAlpha_overlayedWithSim_withFits_pTbin%d.png", ig + 1).Data()));

                  if (fDatOne) delete fDatOne;
                  if (fSimOne) delete fSimOne;
                  delete hDatOne;
                  delete hSimOne;
                };

                for (int ig = 0; ig < nCalibBins; ++ig)
                {
                  SaveGaussianFitsPerPtBin(ig);
                }

              if ((int)vPtCtr.size() > 0)
              {
                TCanvas cMean(
                  TString::Format("c_meanVsPt_%s_dataVsSim_withFits", rKey.c_str()).Data(),
                  "c_meanVsPt_dataVsSim_withFits", 900, 700
                );

                TGraphErrors* gDat = new TGraphErrors((int)vPtCtr.size());
                TGraphErrors* gSim = new TGraphErrors((int)vPtCtr.size());

                for (int i = 0; i < (int)vPtCtr.size(); ++i)
                {
                  gDat->SetPoint(i, vPtCtr[i], vMuDat[i]);
                  gDat->SetPointError(i, vPtErr[i], vMuDatErr[i]);

                  gSim->SetPoint(i, vPtCtr[i], vMuSim[i]);
                  gSim->SetPointError(i, vPtErr[i], vMuSimErr[i]);
                }

                gDat->SetTitle("");
                gDat->SetMarkerStyle(20);
                gDat->SetMarkerSize(1.2);
                gDat->SetMarkerColor(kGreen + 2);
                gDat->SetLineColor(kGreen + 2);
                gDat->SetLineWidth(2);

                gSim->SetMarkerStyle(20);
                gSim->SetMarkerSize(1.2);
                gSim->SetMarkerColor(kOrange + 7);
                gSim->SetLineColor(kOrange + 7);
                gSim->SetLineWidth(2);

                gDat->Draw("AP");
                gDat->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                gDat->GetYaxis()->SetTitle("Gaussian mean of x_{J#gamma}");
                gSim->Draw("P same");

                  TLegend* legM = new TLegend(0.62, 0.16, 0.88, 0.28);
                  legM->SetBorderSize(0);
                  legM->SetFillStyle(0);
                  legM->SetTextFont(42);
                  legM->SetTextSize(0.04);
                  legM->AddEntry(gDat, "DATA (reco)", "p");
                  legM->AddEntry(gSim, "SIM (reco)",  "p");
                  legM->Draw();

                  SaveCanvas(cMean, JoinPath(dirOv, "meanVsPt_reco_integratedAlpha_overlayedWithSim_withFits.png"));

                  // ------------------------------------------------------------------
                  // NEW: (DATA/SIM) mean ratio vs pT with proper error propagation
                  // ------------------------------------------------------------------
                  std::vector<double> vPtCtr_R, vPtErr_R, vMuRatio, vMuRatioErr;
                  vPtCtr_R.reserve(vPtCtr.size());
                  vPtErr_R.reserve(vPtErr.size());
                  vMuRatio.reserve(vPtCtr.size());
                  vMuRatioErr.reserve(vPtCtr.size());

                  for (int i = 0; i < (int)vPtCtr.size(); ++i)
                  {
                    if (vMuDat[i] <= 0.0) continue;
                    if (vMuSim[i] <= 0.0) continue;

                    const double r  = vMuDat[i] / vMuSim[i];
                    const double ed = vMuDatErr[i];
                    const double es = vMuSimErr[i];

                    const double relDat = (vMuDat[i] != 0.0) ? (ed / vMuDat[i]) : 0.0;
                    const double relSim = (vMuSim[i] != 0.0) ? (es / vMuSim[i]) : 0.0;

                    const double er = std::fabs(r) * std::sqrt(relDat*relDat + relSim*relSim);

                    vPtCtr_R.push_back(vPtCtr[i]);
                    vPtErr_R.push_back(vPtErr[i]);
                    vMuRatio.push_back(r);
                    vMuRatioErr.push_back(er);
                  }

                  if ((int)vPtCtr_R.size() > 0)
                  {
                    TCanvas cRatio(
                      TString::Format("c_meanRatioVsPt_%s_dataOverSim_withFits", rKey.c_str()).Data(),
                      "c_meanRatioVsPt_dataOverSim_withFits", 900, 700
                    );

                    TGraphErrors* gRatio = new TGraphErrors((int)vPtCtr_R.size());
                    for (int i = 0; i < (int)vPtCtr_R.size(); ++i)
                    {
                      gRatio->SetPoint(i, vPtCtr_R[i], vMuRatio[i]);
                      gRatio->SetPointError(i, vPtErr_R[i], vMuRatioErr[i]);
                    }

                    gRatio->SetTitle("");
                    gRatio->SetMarkerStyle(20);
                    gRatio->SetMarkerSize(1.2);
                    gRatio->SetMarkerColor(kBlack);
                    gRatio->SetLineColor(kBlack);
                    gRatio->SetLineWidth(2);

                    gRatio->Draw("AP");
                    gRatio->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                    gRatio->GetYaxis()->SetTitle("Gaussian mean ratio: DATA / SIM");

                    double xMin = H.hReco_xJ->GetXaxis()->GetBinLowEdge(1);
                    double xMax = H.hReco_xJ->GetXaxis()->GetBinUpEdge(nPt);

                    TLine* lOne = new TLine(xMin, 1.0, xMax, 1.0);
                    lOne->SetLineStyle(2);
                    lOne->SetLineWidth(2);
                    lOne->Draw();

                    SaveCanvas(cRatio, JoinPath(dirOv, "meanRatioVsPt_reco_integratedAlpha_overlayedWithSim_withFits.png"));

                    // ------------------------------------------------------------------
                    // NEW: Figure-8-style panel plot:
                    //   top: mean vs pT (DATA vs SIM)
                    //   bottom: (DATA/SIM) + constant fit (pT > 13 GeV)
                    // ------------------------------------------------------------------
                    const double fitMinUser = 13.0;
                    const double fitMin     = std::max(fitMinUser, xMin);

                    int nFitPts = 0;
                    for (int i = 0; i < (int)vPtCtr_R.size(); ++i)
                    {
                      if (vPtCtr_R[i] >= fitMin) ++nFitPts;
                    }

                    double jesVal = -1.0;
                    double jesErr =  0.0;

                    TF1* fJES = nullptr;
                    if (nFitPts >= 2)
                    {
                      fJES = new TF1(
                        TString::Format("f_meanRatio_pol0_%s", rKey.c_str()).Data(),
                        "pol0", fitMin, xMax
                      );
                      fJES->SetLineColor(kRed + 1);
                      fJES->SetLineWidth(2);
                      fJES->SetLineStyle(2);

                      gRatio->Fit(fJES, "Q0R");

                      jesVal = fJES->GetParameter(0);
                      jesErr = fJES->GetParError(0);
                    }

                    // Determine y-ranges for the panel
                    double yMinTop =  1e9;
                    double yMaxTop = -1e9;
                    for (int i = 0; i < (int)vPtCtr.size(); ++i)
                    {
                      if (vMuDat[i] > 0.0)
                      {
                        yMinTop = std::min(yMinTop, vMuDat[i] - vMuDatErr[i]);
                        yMaxTop = std::max(yMaxTop, vMuDat[i] + vMuDatErr[i]);
                      }
                      if (vMuSim[i] > 0.0)
                      {
                        yMinTop = std::min(yMinTop, vMuSim[i] - vMuSimErr[i]);
                        yMaxTop = std::max(yMaxTop, vMuSim[i] + vMuSimErr[i]);
                      }
                    }
                    if (yMinTop > yMaxTop)
                    {
                      yMinTop = 0.0;
                      yMaxTop = 1.0;
                    }
                    const double padYTop = std::max(0.02, 0.12 * (yMaxTop - yMinTop));
                    yMinTop -= padYTop;
                    yMaxTop += padYTop;

                    double yMinBot =  1e9;
                    double yMaxBot = -1e9;
                    for (int i = 0; i < (int)vMuRatio.size(); ++i)
                    {
                      yMinBot = std::min(yMinBot, vMuRatio[i] - vMuRatioErr[i]);
                      yMaxBot = std::max(yMaxBot, vMuRatio[i] + vMuRatioErr[i]);
                    }
                    yMinBot = std::min(yMinBot, 1.0);
                    yMaxBot = std::max(yMaxBot, 1.0);
                    if (yMinBot > yMaxBot)
                    {
                      yMinBot = 0.8;
                      yMaxBot = 1.2;
                    }
                    const double padYBot = std::max(0.02, 0.20 * (yMaxBot - yMinBot));
                    yMinBot -= padYBot;
                    yMaxBot += padYBot;

                    TCanvas cPanel(
                      TString::Format("c_meanVsPt_withRatioPanel_%s_dataVsSim", rKey.c_str()).Data(),
                      "c_meanVsPt_withRatioPanel_dataVsSim", 900, 900
                    );

                    TPad padTop(
                      TString::Format("padTop_mean_%s", rKey.c_str()).Data(),
                      "padTop", 0.0, 0.30, 1.0, 1.0
                    );
                    TPad padBot(
                      TString::Format("padBot_mean_%s", rKey.c_str()).Data(),
                      "padBot", 0.0, 0.00, 1.0, 0.30
                    );

                    padTop.SetBottomMargin(0.02);
                    padTop.SetLeftMargin(0.14);
                    padTop.SetRightMargin(0.05);
                    padTop.SetTopMargin(0.08);
                    padTop.SetTicks(1,1);

                    padBot.SetTopMargin(0.02);
                    padBot.SetBottomMargin(0.32);
                    padBot.SetLeftMargin(0.14);
                    padBot.SetRightMargin(0.05);
                    padBot.SetTicks(1,1);

                    padTop.Draw();
                    padBot.Draw();

                    padTop.cd();
                    TH1F* hTop = new TH1F(
                      TString::Format("hMeanPanelTop_%s", rKey.c_str()).Data(),
                      "", 100, xMin, xMax
                    );
                    hTop->SetStats(0);
                    hTop->SetMinimum(yMinTop);
                    hTop->SetMaximum(yMaxTop);
                    hTop->SetTitle("");
                    hTop->GetXaxis()->SetLabelSize(0.0);
                    hTop->GetXaxis()->SetTitleSize(0.0);
                    hTop->GetYaxis()->SetTitle("Gaussian mean of x_{J#gamma}");
                    hTop->GetYaxis()->SetTitleSize(0.07);
                    hTop->GetYaxis()->SetLabelSize(0.06);
                    hTop->GetYaxis()->SetTitleOffset(0.90);
                    hTop->Draw();

                    gDat->Draw("P same");
                    gSim->Draw("P same");

                    TLegend* legP = new TLegend(0.62, 0.14, 0.88, 0.30);
                    legP->SetBorderSize(0);
                    legP->SetFillStyle(0);
                    legP->SetTextFont(42);
                    legP->SetTextSize(0.055);
                    legP->AddEntry(gDat, "DATA (reco)", "p");
                    legP->AddEntry(gSim, "SIM (reco)",  "p");
                    legP->Draw();

                    TLatex lat;
                    lat.SetNDC(true);
                    lat.SetTextFont(42);
                    lat.SetTextSize(0.045);
                    lat.DrawLatex(
                        0.16, 0.92,
                        TString::Format("Photon+ Jet 10 + 20 Sim and Run24pp <x_{J#gamma}>, R = %.1f", RFromKey(rKey)).Data()
                    );

                    padBot.cd();
                    TH1F* hBot = new TH1F(
                        TString::Format("hMeanPanelBot_%s", rKey.c_str()).Data(),
                        "", 100, xMin, xMax
                    );
                    hBot->SetStats(0);
                    hBot->SetMinimum(yMinBot);
                    hBot->SetMaximum(yMaxBot);
                    hBot->SetTitle("");
                    hBot->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                    hBot->GetYaxis()->SetTitle("Ratio DATA / SIM");
                    hBot->GetXaxis()->SetTitleSize(0.12);
                    hBot->GetXaxis()->SetLabelSize(0.10);
                    hBot->GetYaxis()->SetTitleSize(0.10);
                    hBot->GetYaxis()->SetLabelSize(0.10);
                    hBot->GetYaxis()->SetTitleOffset(0.55);
                    hBot->Draw();

                    gRatio->Draw("P same");

                    TLine lOneP(xMin, 1.0, xMax, 1.0);
                    lOneP.SetLineStyle(2);
                    lOneP.SetLineWidth(2);
                    lOneP.Draw("same");

                    if (fJES)
                    {
                        fJES->SetLineColor(kRed + 1);
                        fJES->SetLineStyle(2);
                        fJES->SetLineWidth(2);
                        fJES->Draw("same");

                        TLatex lat2;
                        lat2.SetNDC(true);
                        lat2.SetTextFont(42);
                        lat2.SetTextSize(0.10);
                        lat2.SetTextColor(kRed + 1);
                        lat2.SetTextAlign(31);
                        lat2.DrawLatex(
                          0.95, 0.88,
                          TString::Format("in situ JES = %.4f #pm %.4f", jesVal, jesErr).Data()
                        );
                    }

                      SaveCanvas(cPanel, JoinPath(dirFits, "meanVsPt_withRatioPanel_reco_integratedAlpha_overlayedWithSim_withFits.png"));

                      if (fJES)   delete fJES;
                      delete legP;
                      delete hTop;
                      delete hBot;

                      delete lOne;
                      delete gRatio;
                    }

                    // ------------------------------------------------------------------
                    // NEW: sigma vs pT (DATA vs SIM)
                    // ------------------------------------------------------------------
                    TCanvas cSig(
                      TString::Format("c_sigmaVsPt_%s_dataVsSim_withFits", rKey.c_str()).Data(),
                      "c_sigmaVsPt_dataVsSim_withFits", 900, 700
                    );

                    TGraphErrors* gSigDat = new TGraphErrors((int)vPtCtr.size());
                    TGraphErrors* gSigSim = new TGraphErrors((int)vPtCtr.size());

                    for (int i = 0; i < (int)vPtCtr.size(); ++i)
                    {
                      gSigDat->SetPoint(i, vPtCtr[i], vSigDat[i]);
                      gSigDat->SetPointError(i, vPtErr[i], vSigDatErr[i]);

                      gSigSim->SetPoint(i, vPtCtr[i], vSigSim[i]);
                      gSigSim->SetPointError(i, vPtErr[i], vSigSimErr[i]);
                    }

                    gSigDat->SetTitle("");
                    gSigDat->SetMarkerStyle(20);
                    gSigDat->SetMarkerSize(1.2);
                    gSigDat->SetMarkerColor(kGreen + 2);
                    gSigDat->SetLineColor(kGreen + 2);
                    gSigDat->SetLineWidth(2);

                    gSigSim->SetMarkerStyle(20);
                    gSigSim->SetMarkerSize(1.2);
                    gSigSim->SetMarkerColor(kOrange + 7);
                    gSigSim->SetLineColor(kOrange + 7);
                    gSigSim->SetLineWidth(2);

                    gSigDat->Draw("AP");
                    gSigDat->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                    gSigDat->GetYaxis()->SetTitle("Gaussian #sigma of x_{J#gamma}");
                    gSigSim->Draw("P same");

                    TLegend* legS = new TLegend(0.62, 0.16, 0.88, 0.28);
                    legS->SetBorderSize(0);
                    legS->SetFillStyle(0);
                    legS->SetTextFont(42);
                    legS->SetTextSize(0.04);
                    legS->AddEntry(gSigDat, "DATA (reco)", "p");
                    legS->AddEntry(gSigSim, "SIM (reco)",  "p");
                    legS->Draw();

                    SaveCanvas(cSig, JoinPath(dirOv, "sigmaVsPt_reco_integratedAlpha_overlayedWithSim_withFits.png"));

                    delete legS;
                    delete gSigDat;
                    delete gSigSim;

                    // ------------------------------------------------------------------
                    // NEW: chi2/ndf vs pT (DATA vs SIM)
                    // ------------------------------------------------------------------
                    TCanvas cChi(
                      TString::Format("c_chi2NdfVsPt_%s_dataVsSim_withFits", rKey.c_str()).Data(),
                      "c_chi2NdfVsPt_dataVsSim_withFits", 900, 700
                    );

                    TGraphErrors* gChiDat = new TGraphErrors((int)vPtCtr.size());
                    TGraphErrors* gChiSim = new TGraphErrors((int)vPtCtr.size());

                    for (int i = 0; i < (int)vPtCtr.size(); ++i)
                    {
                      gChiDat->SetPoint(i, vPtCtr[i], vChi2NdfDat[i]);
                      gChiDat->SetPointError(i, vPtErr[i], 0.0);

                      gChiSim->SetPoint(i, vPtCtr[i], vChi2NdfSim[i]);
                      gChiSim->SetPointError(i, vPtErr[i], 0.0);
                    }

                    gChiDat->SetTitle("");
                    gChiDat->SetMarkerStyle(20);
                    gChiDat->SetMarkerSize(1.2);
                    gChiDat->SetMarkerColor(kGreen + 2);
                    gChiDat->SetLineColor(kGreen + 2);
                    gChiDat->SetLineWidth(2);

                    gChiSim->SetMarkerStyle(20);
                    gChiSim->SetMarkerSize(1.2);
                    gChiSim->SetMarkerColor(kOrange + 7);
                    gChiSim->SetLineColor(kOrange + 7);
                    gChiSim->SetLineWidth(2);

                    gChiDat->Draw("AP");
                    gChiDat->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                    gChiDat->GetYaxis()->SetTitle("#chi^{2}/ndf of Gaussian fit");
                    gChiSim->Draw("P same");

                    TLegend* legC = new TLegend(0.62, 0.16, 0.88, 0.28);
                    legC->SetBorderSize(0);
                    legC->SetFillStyle(0);
                    legC->SetTextFont(42);
                    legC->SetTextSize(0.04);
                    legC->AddEntry(gChiDat, "DATA (reco)", "p");
                    legC->AddEntry(gChiSim, "SIM (reco)",  "p");
                    legC->Draw();

                    SaveCanvas(cChi, JoinPath(dirOv, "chi2NdfVsPt_reco_integratedAlpha_overlayedWithSim_withFits.png"));

                    delete legC;
                    delete gChiDat;
                    delete gChiSim;

                    if ((int)vMeanDat.size() == (int)vPtCtr.size() && (int)vMeanSim.size() == (int)vPtCtr.size())
                    {
                    // -----------------------------------------------------------------
                    // NEW: Using GetMean two-panel output
                    // -----------------------------------------------------------------
                    TCanvas cMeanGetMeanWithRatio(
                        TString::Format("c_meanVsPt_withRatioPanel_%s_dataVsSim_usingGetMean", rKey.c_str()).Data(),
                        "c_meanVsPt_withRatioPanel_dataVsSim_usingGetMean", 900, 900
                    );

                    TPad pTop(
                        TString::Format("pTop_%s_usingGetMean", rKey.c_str()).Data(),
                        "pTop_usingGetMean", 0.0, 0.30, 1.0, 1.0
                    );
                    TPad pBot(
                        TString::Format("pBot_%s_usingGetMean", rKey.c_str()).Data(),
                        "pBot_usingGetMean", 0.0, 0.00, 1.0, 0.30
                    );

                    pTop.SetLeftMargin(0.14);
                    pTop.SetRightMargin(0.05);
                    pTop.SetTopMargin(0.08);
                    pTop.SetBottomMargin(0.02);

                    pBot.SetLeftMargin(0.14);
                    pBot.SetRightMargin(0.05);
                    pBot.SetTopMargin(0.03);
                    pBot.SetBottomMargin(0.32);

                    pTop.Draw();
                    pBot.Draw();

                    TGraphErrors* gDatMean = new TGraphErrors((int)vPtCtr.size());
                    TGraphErrors* gSimMean = new TGraphErrors((int)vPtCtr.size());

                    for (int i = 0; i < (int)vPtCtr.size(); ++i)
                    {
                        gDatMean->SetPoint(i, vPtCtr[i], vMeanDat[i]);
                        gDatMean->SetPointError(i, vPtErr[i], vMeanDatErr[i]);

                        gSimMean->SetPoint(i, vPtCtr[i], vMeanSim[i]);
                        gSimMean->SetPointError(i, vPtErr[i], vMeanSimErr[i]);
                    }

                    pTop.cd();

                    gDatMean->SetTitle("");
                    gDatMean->SetMarkerStyle(20);
                    gDatMean->SetMarkerSize(1.2);
                    gDatMean->SetMarkerColor(kBlue + 1);
                    gDatMean->SetLineColor(kBlue + 1);
                    gDatMean->SetLineWidth(2);

                    gSimMean->SetMarkerStyle(20);
                    gSimMean->SetMarkerSize(1.2);
                    gSimMean->SetMarkerColor(kRed + 1);
                    gSimMean->SetLineColor(kRed + 1);
                    gSimMean->SetLineWidth(2);

                    const double xMinMeanPlot = H.hReco_xJ->GetXaxis()->GetBinLowEdge(1);
                    const double xMaxMeanPlot = H.hReco_xJ->GetXaxis()->GetBinUpEdge(nPt);

                    gDatMean->Draw("AP");
                    gDatMean->GetXaxis()->SetLimits(xMinMeanPlot, xMaxMeanPlot);
                    gDatMean->GetXaxis()->SetTitle("");
                    gDatMean->GetYaxis()->SetTitle("Mean of x_{J#gamma} (Using GetMean)");
                    gSimMean->Draw("P same");

                    TLegend* legMean = new TLegend(0.25, 0.16, 0.5, 0.28);
                    legMean->SetBorderSize(0);
                    legMean->SetFillStyle(0);
                    legMean->SetTextFont(42);
                    legMean->SetTextSize(0.04);
                    legMean->AddEntry(gDatMean, "DATA (GetMean)", "p");
                    legMean->AddEntry(gSimMean, "SIM (GetMean)",  "p");
                    legMean->Draw();

                    TLatex ttlMean;
                    ttlMean.SetNDC(true);
                    ttlMean.SetTextFont(42);
                    ttlMean.SetTextSize(0.045);
                    ttlMean.DrawLatex(0.16, 0.92, "RECO x_{J#gamma} mean vs p_{T}^{#gamma} Using GetMean");

                    std::vector<double> vRatioMeanPt;
                    std::vector<double> vRatioMeanPtErr;
                    std::vector<double> vRatioMean;
                    std::vector<double> vRatioMeanErr;

                    for (int i = 0; i < (int)vPtCtr.size(); ++i)
                    {
                        if (vMeanSim[i] <= 0.0) continue;

                        const double ratio = vMeanDat[i] / vMeanSim[i];
                        const double fracDat = (vMeanDat[i] != 0.0) ? (vMeanDatErr[i] / vMeanDat[i]) : 0.0;
                        const double fracSim = (vMeanSim[i] != 0.0) ? (vMeanSimErr[i] / vMeanSim[i]) : 0.0;
                        const double ratioErr = ratio * std::sqrt(fracDat * fracDat + fracSim * fracSim);

                        vRatioMeanPt.push_back(vPtCtr[i]);
                        vRatioMeanPtErr.push_back(vPtErr[i]);
                        vRatioMean.push_back(ratio);
                        vRatioMeanErr.push_back(ratioErr);
                    }

                    pBot.cd();

                    TGraphErrors* gRatioMean = new TGraphErrors((int)vRatioMeanPt.size());
                    for (int i = 0; i < (int)vRatioMeanPt.size(); ++i)
                    {
                        gRatioMean->SetPoint(i, vRatioMeanPt[i], vRatioMean[i]);
                        gRatioMean->SetPointError(i, vRatioMeanPtErr[i], vRatioMeanErr[i]);
                    }

                    gRatioMean->SetTitle("");
                    gRatioMean->SetMarkerStyle(20);
                    gRatioMean->SetMarkerSize(1.1);
                    gRatioMean->SetMarkerColor(kBlack);
                    gRatioMean->SetLineColor(kBlack);
                    gRatioMean->SetLineWidth(2);

                    gRatioMean->Draw("AP");
                    gRatioMean->GetXaxis()->SetLimits(xMinMeanPlot, xMaxMeanPlot);
                    gRatioMean->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                    gRatioMean->GetYaxis()->SetTitle("DATA / SIM");
                    gRatioMean->GetXaxis()->SetTitleSize(0.12);
                    gRatioMean->GetXaxis()->SetLabelSize(0.10);
                    gRatioMean->GetYaxis()->SetTitleSize(0.10);
                    gRatioMean->GetYaxis()->SetTitleOffset(0.55);
                    gRatioMean->GetYaxis()->SetLabelSize(0.09);
                    gRatioMean->SetMinimum(0.6);
                    gRatioMean->SetMaximum(1.4);

                    TLine lOne(
                            xMinMeanPlot, 1.0,
                            xMaxMeanPlot, 1.0
                    );
                    lOne.SetLineStyle(2);
                    lOne.SetLineWidth(2);
                    lOne.Draw("same");

                    TF1* fRatioMean = nullptr;
                    double jesMeanVal = -1.0;
                    double jesMeanErr =  0.0;
                    if ((int)vRatioMeanPt.size() > 0)
                    {
                        const double xLoFit = xMinMeanPlot;
                        const double xHiFit = xMaxMeanPlot;

                        fRatioMean = new TF1(
                              TString::Format("fRatioMean_%s", rKey.c_str()).Data(),
                              "pol0", xLoFit, xHiFit
                        );
                        fRatioMean->SetLineColor(kRed + 1);
                        fRatioMean->SetLineStyle(2);
                        fRatioMean->SetLineWidth(2);
                        gRatioMean->Fit(fRatioMean, "Q0R");
                        fRatioMean->Draw("same");

                        jesMeanVal = fRatioMean->GetParameter(0);
                        jesMeanErr = fRatioMean->GetParError(0);

                        TLatex latMean;
                        latMean.SetNDC(true);
                        latMean.SetTextFont(42);
                        latMean.SetTextSize(0.10);
                        latMean.SetTextColor(kRed + 1);
                        latMean.SetTextAlign(31);
                        latMean.DrawLatex(
                              0.95, 0.88,
                              TString::Format("in situ JES = %.4f #pm %.4f", jesMeanVal, jesMeanErr).Data()
                        );
                    }

                    SaveCanvas(cMeanGetMeanWithRatio,
                        JoinPath(dirFits, "meanVsPt_withRatioPanel_reco_integratedAlpha_overlayedWithSim_usingGetMean.png"));

                    if (fRatioMean) delete fRatioMean;
                    delete legMean;
                    delete gDatMean;
                    delete gSimMean;
                    delete gRatioMean;

                    // -----------------------------------------------------------------
                    //  pure mean overlay (DATA Gaussian means + DATA GetMean only)
                    // -----------------------------------------------------------------
                    TCanvas cMeanOverlay(
                            TString::Format("c_meanVsPt_%s_gaussianAndGetMean", rKey.c_str()).Data(),
                            "c_meanVsPt_gaussianAndGetMean", 900, 700
                    );
                    ApplyCanvasMargins1D(cMeanOverlay);

                    TGraphErrors* gDatGausOv = new TGraphErrors((int)vPtCtr.size());
                    TGraphErrors* gDatMeanOv = new TGraphErrors((int)vPtCtr.size());

                    double yMaxOv = 0.0;

                    for (int i = 0; i < (int)vPtCtr.size(); ++i)
                    {
                        gDatGausOv->SetPoint(i, vPtCtr[i], vMuDat[i]);
                        gDatGausOv->SetPointError(i, vPtErr[i], vMuDatErr[i]);

                        gDatMeanOv->SetPoint(i, vPtCtr[i], vMeanDat[i]);
                        gDatMeanOv->SetPointError(i, vPtErr[i], vMeanDatErr[i]);

                        yMaxOv = std::max(yMaxOv, vMuDat[i]   + vMuDatErr[i]);
                        yMaxOv = std::max(yMaxOv, vMeanDat[i] + vMeanDatErr[i]);
                    }

                    gDatGausOv->SetTitle("");
                    gDatGausOv->SetMarkerStyle(20);
                    gDatGausOv->SetMarkerSize(1.2);
                    gDatGausOv->SetMarkerColor(kGreen + 2);
                    gDatGausOv->SetLineColor(kGreen + 2);
                    gDatGausOv->SetLineWidth(2);

                    gDatMeanOv->SetMarkerStyle(20);
                    gDatMeanOv->SetMarkerSize(1.2);
                    gDatMeanOv->SetMarkerColor(kBlue + 1);
                    gDatMeanOv->SetLineColor(kBlue + 1);
                    gDatMeanOv->SetLineWidth(2);

                    const double xMinMeanOverlay = H.hReco_xJ->GetXaxis()->GetBinLowEdge(1);
                    const double xMaxMeanOverlay = H.hReco_xJ->GetXaxis()->GetBinUpEdge(nPt);

                    gDatGausOv->Draw("AP");
                    gDatGausOv->GetXaxis()->SetLimits(xMinMeanOverlay, xMaxMeanOverlay);
                    gDatGausOv->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                    gDatGausOv->GetYaxis()->SetTitle("Mean of x_{J#gamma}");
                    gDatGausOv->SetMinimum(0.0);
                    gDatGausOv->SetMaximum((yMaxOv > 0.0) ? (1.20 * yMaxOv) : 1.0);

                    gDatMeanOv->Draw("P same");

                    TLegend* legOv = new TLegend(0.18, 0.18, 0.46, 0.32);
                    legOv->SetBorderSize(0);
                    legOv->SetFillStyle(0);
                    legOv->SetTextFont(42);
                    legOv->SetTextSize(0.035);
                    legOv->AddEntry(gDatGausOv, "DATA Gaussian mean", "p");
                    legOv->AddEntry(gDatMeanOv, "DATA GetMean",       "p");
                    legOv->Draw();

                    TLatex ttlOv;
                    ttlOv.SetNDC(true);
                    ttlOv.SetTextFont(42);
                    ttlOv.SetTextSize(0.045);
                    ttlOv.DrawLatex(0.16, 0.92, "RECO x_{J#gamma} mean overlay: DATA Gaussian fits and DATA GetMean");

                    SaveCanvas(cMeanOverlay,
                          JoinPath(dirFits, "meanVsPt_reco_integratedAlpha_overlayedWithSim_gaussianAndGetMean.png"));

                    delete legOv;
                    delete gDatGausOv;
                    delete gDatMeanOv;
                  }

                  delete legM;
                  delete gDat;
                  delete gSim;
              }

              for (auto* f : keepFitFns) delete f;
              for (auto* h1 : keepFitsH) delete h1;
          }
      }

      auto Make3x3Table_xJ_FromTH3 =
        [&](const TH3* h3, const string& outBaseDir, const string& tag, bool logy)
      {
          if (!h3) return;

            const int nAll = h3->GetXaxis()->GetNbins();

            const bool wantLast6 = (tag == "RECO");
            const int  n         = wantLast6 ? std::min(6, nAll) : nAll;

            const int perPage = wantLast6 ? 6 : (n <= 6 ? 6 : 9);

            const int nCols = 3;
            const int nRows = (perPage == 6 ? 2 : 3);

            const int firstBin = wantLast6 ? std::max(1, nAll - n + 1) : 1;
            const int lastStartBin = wantLast6 ? firstBin : n;

            int page = 0;
            for (int start = firstBin; start <= lastStartBin; start += perPage)
            {
              ++page;

              TCanvas c(
                TString::Format("c_tbl_xJ_%s_%s_%s_%s_p%d",
                  ds.label.c_str(),
                  rKey.c_str(),
                  tag.c_str(),
                  logy ? "logy" : "lin",
                  page).Data(),
                "c_tbl_xJ", 1
              );

              c.SetCanvasSize(2400, (perPage == 6 ? 1400 : 2100));
              c.Divide(nCols, nRows, 0.002, 0.002);

              vector<TObject*> keep;
              int pad = 0;

              for (int ib = start; ib < start + perPage; ++ib)
              {
                ++pad;
                c.cd(pad);

                gPad->SetLeftMargin(0.13);
                gPad->SetRightMargin(0.04);
                gPad->SetBottomMargin(0.14);
                gPad->SetTopMargin(0.10);
                gPad->SetLogy(logy);

              if (ib > nAll)
              {
                TLatex t;
                t.SetNDC(true);
                t.SetTextFont(42);
                t.SetTextSize(0.06);
                t.DrawLatex(0.20, 0.55, "EMPTY");
                continue;
              }

              TH1* hx = ProjectY_AtXbin_TH3(
                h3, ib,
                TString::Format("jes3_xJ_tbl_%s_%s_%s_%d",
                  rKey.c_str(), tag.c_str(), logy ? "logy" : "lin", ib).Data()
                );

              if (!hx)
              {
                TLatex t;
                t.SetNDC(true);
                t.SetTextFont(42);
                t.SetTextSize(0.06);
                t.DrawLatex(0.15, 0.55, "MISSING");
                continue;
              }

              hx->SetDirectory(nullptr);
              EnsureSumw2(hx);

              hx->SetLineWidth(2);
              hx->SetMarkerStyle(20);
              hx->SetMarkerSize(1.0);

                hx->SetTitle("");
                hx->GetXaxis()->SetTitle((tag == "TRUTH") ? "x_{J#gamma}^{truth}" : "x_{J#gamma}");
                hx->GetXaxis()->SetRangeUser(0.0, 2.0);
                hx->GetYaxis()->SetTitle((ds.isSim && IsWeightedSIMSelected()) ? "Counts / pb^{-1}" : "Counts");

              if (logy)
              {
                const double minPos = SmallestPositiveBinContent(hx);
                hx->SetMinimum((minPos > 0.0) ? (0.5 * minPos) : 1e-6);
              }

                hx->Draw("E1");

                const double jetPtMin_GeV = static_cast<double>(kJetPtMin);

                const double ptMin = h3->GetXaxis()->GetBinLowEdge(ib);
                const double ptMax = h3->GetXaxis()->GetBinUpEdge(ib);

                const double xAbs  = (ptMax > 0.0) ? (jetPtMin_GeV / ptMax) : -1.0;
                const double xFull = (ptMin > 0.0) ? (jetPtMin_GeV / ptMin) : -1.0;

                gPad->Update();
                const double yMin = gPad->GetUymin();
                const double yMax = gPad->GetUymax();

                TLine* lnAbs = new TLine(xAbs,  yMin, xAbs,  yMax);
                lnAbs->SetLineColor(kBlue + 1);
                lnAbs->SetLineStyle(2);
                lnAbs->SetLineWidth(2);

                TLine* lnFull = new TLine(xFull, yMin, xFull, yMax);
                lnFull->SetLineColor(kRed + 1);
                lnFull->SetLineStyle(2);
                lnFull->SetLineWidth(2);

                if (xAbs > 0.0)  lnAbs->Draw("same");
                if (xFull > 0.0) lnFull->Draw("same");

                TLegend* leg = new TLegend(0.44, 0.60, 0.96, 0.90);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextFont(42);
                leg->SetTextSize(0.034);
                leg->SetEntrySeparation(0.22);

                if (xAbs > 0.0)
                  leg->AddEntry(lnAbs,  TString::Format("x_{J, min}^{abs} = #frac{%.0f}{p_{T, max}^{#gamma}} = %.3f", jetPtMin_GeV, xAbs),  "l");
                if (xFull > 0.0)
                  leg->AddEntry(lnFull, TString::Format("x_{J, min}^{full} = #frac{%.0f}{p_{T, min}^{#gamma}} = %.3f", jetPtMin_GeV, xFull), "l");

                leg->Draw();

                  {
                    TLatex tCuts;
                    tCuts.SetNDC(true);
                    tCuts.SetTextFont(42);
                    tCuts.SetTextAlign(33);
                    tCuts.SetTextSize(0.038);
                    tCuts.DrawLatex(0.92, 0.78, "pp trigger = Photon 4 GeV + MBD NS #geq 1");
                    tCuts.DrawLatex(0.92, 0.70, "auau trigger = MBD NS #geq 2 vtx < 150");
                    tCuts.DrawLatex(0.92, 0.62, TString::Format("|#Delta#phi(#gamma,jet)| > %s", B2BLabel().c_str()).Data());
                    tCuts.DrawLatex(0.92, 0.54, TString::Format("p_{T}^{jet} > %.0f GeV", jetPtMin_GeV).Data());
                  }

                {
                  TLatex t;
                  t.SetNDC(true);
                  t.SetTextFont(42);
                  t.SetTextAlign(13);
                  t.SetTextSize(0.052);
                  t.DrawLatex(0.14, 0.98,
                    TString::Format("%s x_{J#gamma}, p_{T}^{#gamma} = %.0f - %.0f GeV, R = %.1f",
                      tag.c_str(),
                      ptMin, ptMax, R).Data());
                }

                keep.push_back(hx);
                keep.push_back(lnAbs);
                keep.push_back(lnFull);
                keep.push_back(leg);
            }

            string outName;
            if (n <= perPage)
            {
              outName = TString::Format("table3x3_xJ_%s_integratedAlpha%s.png",
                tag.c_str(), logy ? "_logy" : "").Data();
            }
            else
            {
              outName = TString::Format("table3x3_xJ_%s_integratedAlpha%s_page%d.png",
                tag.c_str(), logy ? "_logy" : "", page).Data();
            }

            SaveCanvas(c, JoinPath(outBaseDir, outName));

            for (auto* h : keep) delete h;
          }
        };

        auto Make3x3Table_xJ_AuAuPPOverlay_FromTH3 =
          [&](const TH3* h3Au, const string& outBaseDir, const string& tag, bool logy)
          {
            if (!h3Au) return;
            if (!(isAuAuOnly && !ds.isSim && !ds.centSuffix.empty())) return;
            if (tag != "RECO") return;

            // Never produce a log-y overlay table (user requested no *_logy.png)
            if (logy) return;

            static TFile* fPP = nullptr;
            static TDirectory* dirPP = nullptr;

            if (!fPP)
            {
                fPP = TFile::Open(InputPP().c_str(), "READ");
                if (fPP)
                {
                  dirPP = fPP->GetDirectory(kTriggerPP.c_str());
                  if (!dirPP) dirPP = fPP;
                }
            }

            if (!fPP || !dirPP) return;

            TH3* h3PP = dynamic_cast<TH3*>(dirPP->Get(
              TString::Format("h_JES3_pT_xJ_alpha_%s", rKey.c_str()).Data()
            ));
            if (!h3PP) return;

            const auto& recoBinsAll = UnfoldAnalysisRecoPtBins();
            const int n = std::min(9, (int)recoBinsAll.size());
            if (n <= 0) return;

            const int nCols = 3;
            const int nRows = 3;
            const int perPage = nCols * nRows;

            const string overlayDir = JoinPath(outBaseDir, "auau_pp_overlay");
            EnsureDir(overlayDir);

            std::string auauLeg = "AuAu";
            {
                int clo = -1, chi = -1;
                if (std::sscanf(ds.centSuffix.c_str(), "_cent_%d_%d", &clo, &chi) == 2)
                {
                  auauLeg = TString::Format("AuAu (%d-%d%%)", clo, chi).Data();
                }
            }

            auto NormalizeUnitArea = [&](TH1* h)->void
            {
              if (!h) return;
              const double integ = h->Integral(0, h->GetNbinsX() + 1);
              if (integ > 0.0) h->Scale(1.0 / integ);
            };

            auto ProjectGroupedY_IntegratedAlpha_TH3 =
              [&](const TH3* h3,
                  const PtBin& pb,
                  const std::string& newName)->TH1*
            {
              if (!h3) return nullptr;

              int xbinLo = -1;
              int xbinHi = -1;
              std::vector<double> w;
              if (!XaxisOverlapWeights(h3->GetXaxis(), pb.lo, pb.hi, xbinLo, xbinHi, w)) return nullptr;

              TH1* sum = nullptr;
              const double alphaMax = h3->GetZaxis()->GetXmax();

              for (int xb = xbinLo; xb <= xbinHi; ++xb)
              {
                TH1* h1 = ProjectY_AtXbin_AndAlphaMax_TH3(
                  h3, xb, alphaMax,
                  newName + TString::Format("_xb%d", xb).Data()
                );
                if (!h1) continue;

                const int iw = xb - xbinLo;
                const double ww = (iw >= 0 && iw < (int)w.size()) ? w[iw] : 1.0;

                if (ww <= 0.0)
                {
                  delete h1;
                  continue;
                }

                h1->Scale(ww);

                if (!sum)
                {
                  sum = CloneTH1(h1, newName);
                  if (sum)
                  {
                    sum->Reset("ICES");
                    sum->SetDirectory(nullptr);
                  }
                }

                if (sum) sum->Add(h1);
                delete h1;
              }

              return sum;
            };

            auto StyleAuAuPP = [&](TH1* hxAu, TH1* hxPP, bool normalizeShape)->void
            {
              hxAu->SetLineWidth(2);
              hxAu->SetMarkerStyle(20);
              hxAu->SetMarkerSize(1.0);

              hxPP->SetLineWidth(2);
              hxPP->SetLineColor(kRed + 1);
              hxPP->SetMarkerColor(kRed + 1);
              hxPP->SetMarkerStyle(24);
              hxPP->SetMarkerSize(1.0);

              hxAu->SetTitle("");
              hxAu->GetXaxis()->SetTitle("x_{J#gamma}");
              hxAu->GetXaxis()->SetRangeUser(0.0, 2.0);
              hxAu->GetYaxis()->SetTitle(normalizeShape ? "Normalized Counts" : "Counts");

              const double maxAu = hxAu->GetMaximum();
              const double maxPP = hxPP->GetMaximum();
              const double maxY  = std::max(maxAu, maxPP);

              hxAu->SetMinimum(0.0);
              hxAu->SetMaximum((maxY > 0.0) ? (1.25 * maxY) : 1.0);
            };

            auto Draw3x3Table = [&](bool normalizeShape, const string& outName)
            {
              TCanvas c(
                TString::Format("c_tbl_xJ_auauPP_%s_%s_%s",
                  ds.label.c_str(),
                  rKey.c_str(),
                  normalizeShape ? "shape" : "counts").Data(),
                "c_tbl_xJ_auauPP", 1500, 1200
              );

              c.Divide(nCols, nRows, 0.002, 0.002);

              vector<TObject*> keep;
              int pad = 0;

              for (int i = 0; i < perPage; ++i)
              {
                ++pad;
                c.cd(pad);

                gPad->SetLeftMargin(0.13);
                gPad->SetRightMargin(0.04);
                gPad->SetBottomMargin(0.14);
                gPad->SetTopMargin(0.10);
                gPad->SetLogy(false);

                if (i >= n)
                {
                  TLatex t;
                  t.SetNDC(true);
                  t.SetTextFont(42);
                  t.SetTextSize(0.06);
                  t.DrawLatex(0.20, 0.55, "EMPTY");
                  continue;
                }

                const PtBin& pb = recoBinsAll[i];

                TH1* hxAu = ProjectGroupedY_IntegratedAlpha_TH3(
                  h3Au, pb,
                  TString::Format("jes3_xJ_tbl_auau_%s_%s_%s",
                    rKey.c_str(), tag.c_str(), pb.folder.c_str()).Data()
                );

                TH1* hxPP = ProjectGroupedY_IntegratedAlpha_TH3(
                  h3PP, pb,
                  TString::Format("jes3_xJ_tbl_pp_%s_%s_%s",
                    rKey.c_str(), tag.c_str(), pb.folder.c_str()).Data()
                );

                if (!hxAu || !hxPP)
                {
                  TLatex t;
                  t.SetNDC(true);
                  t.SetTextFont(42);
                  t.SetTextSize(0.06);
                  t.DrawLatex(0.15, 0.55, "MISSING");
                  if (hxAu) delete hxAu;
                  if (hxPP) delete hxPP;
                  continue;
                }

                hxAu->SetDirectory(nullptr);
                hxPP->SetDirectory(nullptr);
                EnsureSumw2(hxAu);
                EnsureSumw2(hxPP);

                if (normalizeShape)
                {
                  NormalizeUnitArea(hxAu);
                  NormalizeUnitArea(hxPP);
                }

                StyleAuAuPP(hxAu, hxPP, normalizeShape);

                hxAu->Draw("E1");
                hxPP->Draw("E1 same");

                TLegend* leg = new TLegend(0.18, 0.78, 0.52, 0.90);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextFont(42);
                leg->SetTextSize(0.038);
                leg->SetEntrySeparation(0.22);
                leg->AddEntry(hxPP, "pp", "ep");
                leg->AddEntry(hxAu, auauLeg.c_str(), "ep");
                leg->Draw();

                const double jetPtMin_GeV = static_cast<double>(kJetPtMin);

                {
                  TLatex tCuts;
                  tCuts.SetNDC(true);
                  tCuts.SetTextFont(42);
                  tCuts.SetTextAlign(33);
                  tCuts.SetTextSize(0.038);
                  tCuts.DrawLatex(0.92, 0.62, TString::Format("|#Delta#phi(#gamma,jet)| > %s", B2BLabel().c_str()).Data());
                  tCuts.DrawLatex(0.92, 0.54, TString::Format("p_{T}^{jet} > %.0f GeV", jetPtMin_GeV).Data());
                }

                {
                  TLatex t;
                  t.SetNDC(true);
                  t.SetTextFont(42);
                  t.SetTextAlign(13);
                  t.SetTextSize(0.052);
                  t.DrawLatex(0.14, 0.98,
                    TString::Format("RECO x_{J#gamma}, p_{T}^{#gamma} = %d - %d GeV, R = %.1f",
                      pb.lo, pb.hi, R).Data());
                }

                keep.push_back(hxAu);
                keep.push_back(hxPP);
                keep.push_back(leg);
              }

              SaveCanvas(c, JoinPath(overlayDir, outName));

              for (auto* h : keep) delete h;
            };

            Draw3x3Table(false,
              TString::Format("table2x3_xJ_%s_integratedAlpha_counts.png", tag.c_str()).Data());

            Draw3x3Table(true,
              TString::Format("table2x3_xJ_%s_integratedAlpha.png", tag.c_str()).Data());

            Draw3x3Table(true,
              TString::Format("table3x3_xJ_%s_integratedAlpha_shapeNormalized.png", tag.c_str()).Data());

              // ----------------------------------------
              // (B) INDIVIDUAL PNG for EVERY analysis pT bin
              // ----------------------------------------
            for (int i = 0; i < (int)recoBinsAll.size(); ++i)
            {
                const PtBin& pb = recoBinsAll[i];

                TH1* hxAu = ProjectGroupedY_IntegratedAlpha_TH3(
                  h3Au, pb,
                  TString::Format("jes3_xJ_indiv_auau_%s_%s_%s",
                    rKey.c_str(), tag.c_str(), pb.folder.c_str()).Data()
                );

                TH1* hxPP = ProjectGroupedY_IntegratedAlpha_TH3(
                  h3PP, pb,
                  TString::Format("jes3_xJ_indiv_pp_%s_%s_%s",
                    rKey.c_str(), tag.c_str(), pb.folder.c_str()).Data()
                );

                if (!hxAu || !hxPP)
                {
                  if (hxAu) delete hxAu;
                  if (hxPP) delete hxPP;
                  continue;
                }

                hxAu->SetDirectory(nullptr);
                hxPP->SetDirectory(nullptr);
                EnsureSumw2(hxAu);
                EnsureSumw2(hxPP);

                NormalizeUnitArea(hxAu);
                NormalizeUnitArea(hxPP);

                StyleAuAuPP(hxAu, hxPP, true);

                const double ptMin = pb.lo;
                const double ptMax = pb.hi;

                TCanvas c(
                  TString::Format("c_xJ_auauPP_%s_%s_%d",
                    ds.label.c_str(), rKey.c_str(), i + 1).Data(),
                  "c_xJ_auauPP", 900, 700
                );

                c.SetTopMargin(0.10);
                c.SetBottomMargin(0.14);
                c.SetLeftMargin(0.13);
                c.SetRightMargin(0.05);

                hxAu->Draw("E1");
                hxPP->Draw("E1 same");

                TLegend* leg = new TLegend(0.62, 0.78, 0.92, 0.90);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextFont(42);
                leg->SetTextSize(0.040);
                leg->SetEntrySeparation(0.22);
                leg->AddEntry(hxPP, "pp",   "ep");
                leg->AddEntry(hxAu, auauLeg.c_str(), "ep");
                leg->Draw();

                const double jetPtMin_GeV = static_cast<double>(kJetPtMin);

                {
                  TLatex tCuts;
                  tCuts.SetNDC(true);
                  tCuts.SetTextFont(42);
                  tCuts.SetTextAlign(33);
                  tCuts.SetTextSize(0.038);
                  tCuts.DrawLatex(0.92, 0.62, TString::Format("|#Delta#phi(#gamma,jet)| > %s", B2BLabel().c_str()).Data());
                  tCuts.DrawLatex(0.92, 0.54, TString::Format("p_{T}^{jet} > %.0f GeV", jetPtMin_GeV).Data());
                }

                {
                  TLatex t;
                  t.SetNDC(true);
                  t.SetTextFont(42);
                  t.SetTextAlign(13);
                  t.SetTextSize(0.052);
                  t.DrawLatex(0.14, 0.98,
                    TString::Format("RECO x_{J#gamma}, p_{T}^{#gamma} = %.0f - %.0f GeV, R = %.1f",
                      ptMin, ptMax, R).Data());
                }

                const string outName = TString::Format("xJ_%s_integratedAlpha_auau_pp_pTgamma_%d_%d.png",
                  tag.c_str(), (int)ptMin, (int)ptMax).Data();

                SaveCanvas(c, JoinPath(overlayDir, outName));

                delete leg;
                delete hxAu;
                delete hxPP;
            }
        };

        auto Make3x3Table_xJ_WithWithoutUESub_FromTH3 =
          [&](const TH3* h3NoUE, const string& outBaseDir, const string& tag, bool logy)
        {
          if (!h3NoUE) return;
          if (!(isAuAuOnly && !ds.isSim && !ds.centSuffix.empty())) return;
          if (tag != "RECO") return;

          // Never produce log-y overlays for this comparison
          if (logy) return;

          static TFile* fUE = nullptr;
          static TDirectory* dirUE = nullptr;

          if (!fUE)
          {
              fUE = TFile::Open(InputAuAu().c_str(), "READ");
              if (fUE)
              {
                dirUE = fUE->GetDirectory(ds.trigger.c_str());
                if (!dirUE) dirUE = fUE;
              }
          }

          if (!fUE || !dirUE) return;

          // Fetch the SAME TH3 from the UE-subtracted file, matching centrality
          TH3* h3UE = nullptr;
          {
            const std::string baseName =
              TString::Format("h_JES3_pT_xJ_alpha_%s", rKey.c_str()).Data();

            TObject* obj = dirUE->Get((baseName + ds.centSuffix).c_str());
            if (!obj) obj = dirUE->Get(baseName.c_str());
            h3UE = dynamic_cast<TH3*>(obj);
          }
          if (!h3UE) return;

          const int nAll = h3NoUE->GetXaxis()->GetNbins();
          if (nAll <= 0) return;

          // 2x3 table uses the FIRST 6 pT bins
          const int n = std::min(6, nAll);
          const int perPage = 6;
          const int nCols = 3;
          const int nRows = 2;

          const int firstBin = 1;
          const int lastStartBin = 1;

            const string overlayDir = JoinPath(outBaseDir, "with_withoutUEsub");
            EnsureDir(overlayDir);

            std::string centPar = "";
            {
              int clo = -1, chi = -1;
              if (std::sscanf(ds.centSuffix.c_str(), "_cent_%d_%d", &clo, &chi) == 2)
              {
                centPar = TString::Format(" (%d-%d%%)", clo, chi).Data();
              }
            }

          auto NormalizeUnitArea = [&](TH1* h)->void
          {
            if (!h) return;
            const double integ = h->Integral(0, h->GetNbinsX() + 1);
            if (integ > 0.0) h->Scale(1.0 / integ);
          };

          auto StyleNoUEvsUE = [&](TH1* hNoUE, TH1* hUE)->void
          {
            // without UE sub (baseline file): black
            hNoUE->SetLineWidth(2);
            hNoUE->SetMarkerStyle(20);
            hNoUE->SetMarkerSize(1.0);

            // with UE sub (kInAuAuGoldNew): red
            hUE->SetLineWidth(2);
            hUE->SetLineColor(kRed + 1);
            hUE->SetMarkerColor(kRed + 1);
            hUE->SetMarkerStyle(24);
            hUE->SetMarkerSize(1.0);

            hNoUE->SetTitle("");
            hNoUE->GetXaxis()->SetTitle("x_{J#gamma}");
            hNoUE->GetXaxis()->SetRangeUser(0.0, 2.0);
            hNoUE->GetYaxis()->SetTitle("Normalized Counts");

            const double maxNoUE = hNoUE->GetMaximum();
            const double maxUE   = hUE->GetMaximum();
            const double maxY    = std::max(maxNoUE, maxUE);

            hNoUE->SetMinimum(0.0);
            hNoUE->SetMaximum((maxY > 0.0) ? (1.25 * maxY) : 1.0);
          };

          // -----------------------
          // (A) 2x3 TABLE (first 6)
          // -----------------------
          {
            int page = 0;
            for (int start = firstBin; start <= lastStartBin; start += perPage)
            {
              ++page;

              TCanvas c(
                TString::Format("c_tbl_xJ_withWithoutUE_%s_%s_%s_lin_p%d",
                  ds.label.c_str(),
                  rKey.c_str(),
                  tag.c_str(),
                  page).Data(),
                "c_tbl_xJ_withWithoutUE", 1
              );

              c.SetCanvasSize(2400, 1400);
              c.Divide(nCols, nRows, 0.002, 0.002);

              vector<TObject*> keep;
              int pad = 0;

              for (int k = 0; k < perPage; ++k)
              {
                const int ib = start + k;
                ++pad;
                c.cd(pad);

                gPad->SetLeftMargin(0.13);
                gPad->SetRightMargin(0.04);
                gPad->SetBottomMargin(0.14);
                gPad->SetTopMargin(0.10);
                gPad->SetLogy(false);

                if (ib > n)
                {
                  TLatex t;
                  t.SetNDC(true);
                  t.SetTextFont(42);
                  t.SetTextSize(0.06);
                  t.DrawLatex(0.20, 0.55, "EMPTY");
                  continue;
                }

                TH1* hxNoUE = ProjectY_AtXbin_TH3(
                  h3NoUE, ib,
                  TString::Format("jes3_xJ_tbl_noUE_%s_%s_lin_%d",
                    rKey.c_str(), tag.c_str(), ib).Data()
                );

                TH1* hxUE = ProjectY_AtXbin_TH3(
                  h3UE, ib,
                  TString::Format("jes3_xJ_tbl_withUE_%s_%s_lin_%d",
                    rKey.c_str(), tag.c_str(), ib).Data()
                );

                if (!hxNoUE || !hxUE)
                {
                  TLatex t;
                  t.SetNDC(true);
                  t.SetTextFont(42);
                  t.SetTextSize(0.06);
                  t.DrawLatex(0.15, 0.55, "MISSING");
                  if (hxNoUE) delete hxNoUE;
                  if (hxUE)   delete hxUE;
                  continue;
                }

                hxNoUE->SetDirectory(nullptr);
                hxUE->SetDirectory(nullptr);
                EnsureSumw2(hxNoUE);
                EnsureSumw2(hxUE);

                NormalizeUnitArea(hxNoUE);
                NormalizeUnitArea(hxUE);

                StyleNoUEvsUE(hxNoUE, hxUE);

                hxNoUE->Draw("E1");
                hxUE->Draw("E1 same");

                TLegend* leg = new TLegend(0.18, 0.78, 0.56, 0.90);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextFont(42);
                leg->SetTextSize(0.038);
                leg->SetEntrySeparation(0.22);
                leg->AddEntry(hxUE,   "with UE sub",    "ep");
                leg->AddEntry(hxNoUE, "without UE sub", "ep");
                leg->Draw();

                const double jetPtMin_GeV = static_cast<double>(kJetPtMin);

                const double ptMin = h3NoUE->GetXaxis()->GetBinLowEdge(ib);
                const double ptMax = h3NoUE->GetXaxis()->GetBinUpEdge(ib);

                {
                  TLatex tCuts;
                  tCuts.SetNDC(true);
                  tCuts.SetTextFont(42);
                  tCuts.SetTextAlign(33);
                  tCuts.SetTextSize(0.038);
                  tCuts.DrawLatex(0.92, 0.62, TString::Format("|#Delta#phi(#gamma,jet)| > %s", B2BLabel().c_str()).Data());
                  tCuts.DrawLatex(0.92, 0.54, TString::Format("p_{T}^{jet} > %.0f GeV", jetPtMin_GeV).Data());
                }

                {
                  TLatex t;
                  t.SetNDC(true);
                  t.SetTextFont(42);
                  t.SetTextAlign(13);
                  t.SetTextSize(0.052);
                    t.DrawLatex(0.14, 0.98,
                      TString::Format("RECO x_{J#gamma}, p_{T}^{#gamma} = %.0f - %.0f GeV%s, R = %.1f",
                        ptMin, ptMax, centPar.c_str(), R).Data());
                }

                keep.push_back(hxNoUE);
                keep.push_back(hxUE);
                keep.push_back(leg);
              }

              const string outName = TString::Format("table2x3_xJ_%s_integratedAlpha.png", tag.c_str()).Data();
              SaveCanvas(c, JoinPath(overlayDir, outName));

              for (auto* h : keep) delete h;
            }
          }

          // ----------------------------------------
          // (B) INDIVIDUAL PNG for EVERY pT bin (1..nAll)
          // ----------------------------------------
          for (int ib = 1; ib <= nAll; ++ib)
          {
            TH1* hxNoUE = ProjectY_AtXbin_TH3(
              h3NoUE, ib,
              TString::Format("jes3_xJ_indiv_noUE_%s_%s_lin_%d",
                rKey.c_str(), tag.c_str(), ib).Data()
            );

            TH1* hxUE = ProjectY_AtXbin_TH3(
              h3UE, ib,
              TString::Format("jes3_xJ_indiv_withUE_%s_%s_lin_%d",
                rKey.c_str(), tag.c_str(), ib).Data()
            );

            if (!hxNoUE || !hxUE)
            {
              if (hxNoUE) delete hxNoUE;
              if (hxUE)   delete hxUE;
              continue;
            }

            hxNoUE->SetDirectory(nullptr);
            hxUE->SetDirectory(nullptr);
            EnsureSumw2(hxNoUE);
            EnsureSumw2(hxUE);

            NormalizeUnitArea(hxNoUE);
            NormalizeUnitArea(hxUE);

            StyleNoUEvsUE(hxNoUE, hxUE);

            const double ptMin = h3NoUE->GetXaxis()->GetBinLowEdge(ib);
            const double ptMax = h3NoUE->GetXaxis()->GetBinUpEdge(ib);

            TCanvas c(
              TString::Format("c_xJ_withWithoutUE_%s_%s_%d",
                ds.label.c_str(), rKey.c_str(), ib).Data(),
              "c_xJ_withWithoutUE", 900, 700
            );

            c.SetTopMargin(0.10);
            c.SetBottomMargin(0.14);
            c.SetLeftMargin(0.13);
            c.SetRightMargin(0.05);

            hxNoUE->Draw("E1");
            hxUE->Draw("E1 same");

            TLegend* leg = new TLegend(0.18, 0.78, 0.56, 0.90);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextFont(42);
            leg->SetTextSize(0.040);
            leg->SetEntrySeparation(0.22);
            leg->AddEntry(hxUE,   "with UE sub",    "ep");
            leg->AddEntry(hxNoUE, "without UE sub", "ep");
            leg->Draw();

            const double jetPtMin_GeV = static_cast<double>(kJetPtMin);

            {
              TLatex tCuts;
              tCuts.SetNDC(true);
              tCuts.SetTextFont(42);
              tCuts.SetTextAlign(33);
              tCuts.SetTextSize(0.038);
              tCuts.DrawLatex(0.92, 0.62, TString::Format("|#Delta#phi(#gamma,jet)| > %s", B2BLabel().c_str()).Data());
              tCuts.DrawLatex(0.92, 0.54, TString::Format("p_{T}^{jet} > %.0f GeV", jetPtMin_GeV).Data());
            }

            {
              TLatex t;
              t.SetNDC(true);
              t.SetTextFont(42);
              t.SetTextAlign(13);
              t.SetTextSize(0.052);
              t.DrawLatex(0.14, 0.98,
                TString::Format("RECO x_{J#gamma}, p_{T}^{#gamma} = %.0f - %.0f GeV, R = %.1f",
                  ptMin, ptMax, R).Data());
            }

            const string outName = TString::Format("xJ_%s_integratedAlpha_with_withoutUEsub_pTgamma_%d_%d.png",
              tag.c_str(), (int)ptMin, (int)ptMax).Data();

            SaveCanvas(c, JoinPath(overlayDir, outName));

            delete leg;
            delete hxNoUE;
            delete hxUE;
          }
        };

        // NEW: 3x3 tables for RECO alpha cuts (linear-y only; requested)
        auto Make3x3Table_xJ_FromTH3_AlphaCut =
          [&](const TH3* h3, const string& outBaseDir, const string& tag, double alphaMax)
        {
          if (!h3) return;

        const int n = h3->GetXaxis()->GetNbins();
        const int perPage = 9;

        int page = 0;
        for (int start = 1; start <= n; start += perPage)
        {
          ++page;

          TCanvas c(
            TString::Format("c_tbl_xJ_%s_%s_%s_alphaLT_%.2f_p%d",
              ds.label.c_str(),
              rKey.c_str(),
              tag.c_str(),
              alphaMax,
              page).Data(),
            "c_tbl_xJ_alphaCut", 1500, 1200
          );

          c.Divide(3,3, 0.001, 0.001);

          std::vector<TH1*> keep;
          keep.reserve(perPage);

          for (int k = 0; k < perPage; ++k)
          {
            const int ib = start + k;
            c.cd(k+1);

            gPad->SetLeftMargin(0.14);
            gPad->SetRightMargin(0.05);
            gPad->SetBottomMargin(0.14);
            gPad->SetTopMargin(0.10);
            gPad->SetLogy(false);

            if (ib > n)
            {
              TLatex t;
              t.SetNDC(true);
              t.SetTextFont(42);
              t.SetTextSize(0.06);
              t.DrawLatex(0.20, 0.55, "EMPTY");
              continue;
            }

            TH1* hx = ProjectY_AtXbin_AndAlphaMax_TH3(
              h3, ib, alphaMax,
              TString::Format("jes3_xJ_tbl_%s_%s_alphaLT%.2f_%d",
                rKey.c_str(), tag.c_str(), alphaMax, ib).Data()
            );

            if (!hx)
            {
              TLatex t;
              t.SetNDC(true);
              t.SetTextFont(42);
              t.SetTextSize(0.06);
              t.DrawLatex(0.15, 0.55, "MISSING");
              continue;
            }

            hx->SetDirectory(nullptr);
            EnsureSumw2(hx);

            hx->SetLineWidth(2);
            hx->SetMarkerStyle(20);
            hx->SetMarkerSize(1.0);

              hx->SetTitle("");
              hx->GetXaxis()->SetTitle((tag == "TRUTH") ? "x_{J#gamma}^{truth}" : "x_{J#gamma}");
              hx->GetYaxis()->SetTitle((ds.isSim && IsWeightedSIMSelected()) ? "Counts / pb^{-1}" : "Counts");
              hx->Draw("E1");

            const string ptLab = AxisBinLabel(h3->GetXaxis(), ib, "GeV", 0);

            vector<string> lines;
            lines.push_back(TString::Format("JES3 %s: x_{J#gamma}", tag.c_str()).Data());
            lines.push_back(TString::Format("alpha cut: #alpha < %.2f", alphaMax).Data());
            lines.push_back(TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data());
            lines.push_back(TString::Format("p_{T}^{#gamma}: %s", ptLab.c_str()).Data());
            DrawLatexLines(0.16, 0.90, lines, 0.040, 0.050);

            keep.push_back(hx);
          }

          string outName;
          if (n <= perPage)
          {
            outName = TString::Format("table3x3_xJ_%s_alphaLT%.2f.png", tag.c_str(), alphaMax).Data();
          }
          else
          {
            outName = TString::Format("table3x3_xJ_%s_alphaLT%.2f_page%d.png", tag.c_str(), alphaMax, page).Data();
          }

          SaveCanvas(c, JoinPath(outBaseDir, outName));

          for (auto* h : keep) delete h;
        }
      };

        // integrated-alpha tables (linear + logy)
        Make3x3Table_xJ_FromTH3(H.hReco_xJ,                    D.dirXJProjReco,                   "RECO",  false);
        Make3x3Table_xJ_AuAuPPOverlay_FromTH3(H.hReco_xJ,      D.dirXJProjReco,                   "RECO",  false);

        // =================================================================
        // NEW: xJoverlays — multi-centrality AuAu + PP overlay for first
        //      pT bin (10-12 GeV), output to auau/<trigger>/xJoverlays/
        //      Runs once (static guard) for rKey == "r04" only.
        // =================================================================
        {
          static bool s_xJoverlaysDone = false;
          if (!s_xJoverlaysDone
              && isAuAuOnly && !ds.isSim && !ds.centSuffix.empty()
              && rKey == "r04")
          {
            s_xJoverlaysDone = true;

            // Build trigger-level output: <OutputAuAu()>/<trigger>/xJoverlays/
            const string trigBase = JoinPath(OutputAuAu(), ds.trigger);
            const string xJovDir  = JoinPath(trigBase, "xJoverlays");
            EnsureDir(xJovDir);

            // Open PP file
            TFile* fPPov = TFile::Open(InputPP().c_str(), "READ");
            TDirectory* dirPPov = nullptr;
            if (fPPov)
            {
                dirPPov = fPPov->GetDirectory(kTriggerPP.c_str());
                if (!dirPPov) dirPPov = fPPov;
            }
            // Open AuAu file (baseline, non-UE-subtracted)
            TFile* fAAov = TFile::Open(InputAuAu().c_str(), "READ");
            TDirectory* dirAAov = nullptr;
            if (fAAov)
            {
                dirAAov = fAAov->GetDirectory(ds.trigger.c_str());
                if (!dirAAov) dirAAov = fAAov;
            }

            if (fPPov && dirPPov && fAAov && dirAAov)
            {
              const string h3name = "h_JES3_pT_xJ_alpha_r04";

              // PP TH3
              TH3* h3PP = dynamic_cast<TH3*>(dirPPov->Get(h3name.c_str()));

              // AuAu TH3 per centrality
              struct CentEntry { const char* suffix; const char* label; int color; };
              const CentEntry centEntries[] = {
                { "_cent_0_10",  "AuAu 0-10%",  kBlue + 1  },
                { "_cent_10_20", "AuAu 10-20%", kGreen + 2 },
                { "_cent_20_40", "AuAu 20-40%", kMagenta + 1 },
                { "_cent_40_60", "AuAu 40-60%", kBlack     },
              };
              const int nCent = 4;

              TH3* h3AA[4] = { nullptr, nullptr, nullptr, nullptr };
              for (int ic = 0; ic < nCent; ++ic)
              {
                const string fullName = h3name + centEntries[ic].suffix;
                h3AA[ic] = dynamic_cast<TH3*>(dirAAov->Get(fullName.c_str()));
              }

              // First pT bin: 10-12 GeV
              const auto& recoBinsAll = UnfoldAnalysisRecoPtBins();
              const PtBin* pb10_12 = nullptr;
              for (const auto& pb : recoBinsAll)
              {
                if (pb.lo == 10 && pb.hi == 12) { pb10_12 = &pb; break; }
              }

              if (h3PP && pb10_12)
              {
                // Reuse the same grouped-projection helper already in scope
                auto ProjectXJ = [&](const TH3* h3, const PtBin& pb, const string& name)->TH1*
                {
                  if (!h3) return nullptr;
                  int xbinLo = -1, xbinHi = -1;
                  std::vector<double> w;
                  if (!XaxisOverlapWeights(h3->GetXaxis(), pb.lo, pb.hi, xbinLo, xbinHi, w))
                    return nullptr;

                  TH1* sum = nullptr;
                  const double alphaMax = h3->GetZaxis()->GetXmax();
                  for (int xb = xbinLo; xb <= xbinHi; ++xb)
                  {
                    TH1* h1 = ProjectY_AtXbin_AndAlphaMax_TH3(
                      h3, xb, alphaMax,
                      name + TString::Format("_xb%d", xb).Data()
                    );
                    if (!h1) continue;
                    const int iw = xb - xbinLo;
                    const double ww = (iw >= 0 && iw < (int)w.size()) ? w[iw] : 1.0;
                    if (ww <= 0.0) { delete h1; continue; }
                    h1->Scale(ww);
                    if (!sum)
                    {
                      sum = CloneTH1(h1, name);
                      if (sum) { sum->Reset("ICES"); sum->SetDirectory(nullptr); }
                    }
                    if (sum) sum->Add(h1);
                    delete h1;
                  }
                  return sum;
                };

                auto NormUnit = [](TH1* h) {
                  if (!h) return;
                  const double I = h->Integral(0, h->GetNbinsX() + 1);
                  if (I > 0.0) h->Scale(1.0 / I);
                };

                // Project PP
                TH1* hxPP = ProjectXJ(h3PP, *pb10_12, "xJov_pp_10_12");
                if (hxPP) { hxPP->SetDirectory(nullptr); EnsureSumw2(hxPP); NormUnit(hxPP); }

                // Project AuAu centralities
                TH1* hxAA[4] = { nullptr, nullptr, nullptr, nullptr };
                for (int ic = 0; ic < nCent; ++ic)
                {
                  if (!h3AA[ic]) continue;
                  hxAA[ic] = ProjectXJ(h3AA[ic], *pb10_12,
                    TString::Format("xJov_auau%s_10_12", centEntries[ic].suffix).Data());
                  if (hxAA[ic]) { hxAA[ic]->SetDirectory(nullptr); EnsureSumw2(hxAA[ic]); NormUnit(hxAA[ic]); }
                }

                // Check at least PP + one AuAu exist
                bool anyAA = false;
                for (int ic = 0; ic < nCent; ++ic) if (hxAA[ic]) anyAA = true;

                if (hxPP && anyAA)
                {
                  TCanvas cOv("c_xJov_multiCent", "c_xJov_multiCent", 900, 700);
                  cOv.SetTopMargin(0.08);
                  cOv.SetBottomMargin(0.14);
                  cOv.SetLeftMargin(0.13);
                  cOv.SetRightMargin(0.05);

                  // Style PP: open red circles
                  hxPP->SetLineWidth(2);
                  hxPP->SetLineColor(kRed + 1);
                  hxPP->SetMarkerColor(kRed + 1);
                  hxPP->SetMarkerStyle(24);
                  hxPP->SetMarkerSize(1.1);

                  // Style AuAu: closed circles, color per centrality
                  for (int ic = 0; ic < nCent; ++ic)
                  {
                    if (!hxAA[ic]) continue;
                    hxAA[ic]->SetLineWidth(2);
                    hxAA[ic]->SetLineColor(centEntries[ic].color);
                    hxAA[ic]->SetMarkerColor(centEntries[ic].color);
                    hxAA[ic]->SetMarkerStyle(20);
                    hxAA[ic]->SetMarkerSize(1.0);
                  }

                  // Determine y-max across all curves
                  double ymax = hxPP->GetMaximum();
                  for (int ic = 0; ic < nCent; ++ic)
                    if (hxAA[ic]) ymax = std::max(ymax, hxAA[ic]->GetMaximum());

                  // Use PP as frame
                  hxPP->SetTitle("");
                  hxPP->GetXaxis()->SetTitle("x_{J#gamma}");
                  hxPP->GetXaxis()->SetRangeUser(0.0, 2.0);
                  hxPP->GetYaxis()->SetTitle("Normalized Counts");
                  hxPP->SetMinimum(0.0);
                  hxPP->SetMaximum(ymax * 1.35);

                  hxPP->Draw("E1");
                  for (int ic = 0; ic < nCent; ++ic)
                    if (hxAA[ic]) hxAA[ic]->Draw("E1 same");

                  // Legend in top-right
                  TLegend* legOv = new TLegend(0.60, 0.65, 0.93, 0.92);
                  legOv->SetBorderSize(0);
                  legOv->SetFillStyle(0);
                  legOv->SetTextFont(42);
                  legOv->SetTextSize(0.035);
                  legOv->AddEntry(hxPP, "pp", "ep");
                  for (int ic = 0; ic < nCent; ++ic)
                    if (hxAA[ic]) legOv->AddEntry(hxAA[ic], centEntries[ic].label, "ep");
                  legOv->Draw();

                  // Title annotation
                  {
                    TLatex t;
                    t.SetNDC(true);
                    t.SetTextFont(42);
                    t.SetTextAlign(13);
                    t.SetTextSize(0.042);
                    t.DrawLatex(0.14, 0.98,
                      TString::Format("RECO x_{J#gamma}, p_{T}^{#gamma} = 10-12 GeV, R = 0.4").Data());
                  }

                  SaveCanvas(cOv, JoinPath(xJovDir, "xJ_RECO_multiCent_pp_overlay_pTgamma_10_12.png"));

                  delete legOv;
                }

                // Cleanup
                if (hxPP) delete hxPP;
                for (int ic = 0; ic < nCent; ++ic) if (hxAA[ic]) delete hxAA[ic];
              }
            }

            if (fPPov) { fPPov->Close(); delete fPPov; }
            if (fAAov) { fAAov->Close(); delete fAAov; }
          }
        }

        Make3x3Table_xJ_WithWithoutUESub_FromTH3(H.hReco_xJ,   D.dirXJProjReco,                   "RECO",  false);
        Make3x3Table_xJ_FromTH3(H.hRecoTruthPhoTagged_xJ,      D.dirXJProjRecoTruthPhoTagged,     "RECO",  false);
        Make3x3Table_xJ_FromTH3(H.hRecoTruthTagged_xJ,    D.dirXJProjRecoTruthTagged,        "RECO",  false);

        Make3x3Table_xJ_FromTH3(H.hTrut_xJ,               D.dirXJProjTruthRecoCond,           "TRUTH", false);
        Make3x3Table_xJ_FromTH3(H.hTrutNoJM_xJ,           D.dirXJProjTruthRecoCondNoJetMatch, "TRUTH", false);
        Make3x3Table_xJ_FromTH3(H.hTrutPure_xJ,           D.dirXJProjTruthPure,               "TRUTH", false);

        Make3x3Table_xJ_FromTH3(H.hReco_xJ,               D.dirXJProjReco,                   "RECO",  true);
        Make3x3Table_xJ_FromTH3(H.hRecoTruthPhoTagged_xJ, D.dirXJProjRecoTruthPhoTagged,     "RECO",  true);
        Make3x3Table_xJ_FromTH3(H.hRecoTruthTagged_xJ,    D.dirXJProjRecoTruthTagged,        "RECO",  true);

        Make3x3Table_xJ_FromTH3(H.hTrut_xJ,               D.dirXJProjTruthRecoCond,           "TRUTH", true);
        Make3x3Table_xJ_FromTH3(H.hTrutNoJM_xJ,           D.dirXJProjTruthRecoCondNoJetMatch, "TRUTH", true);
        Make3x3Table_xJ_FromTH3(H.hTrutPure_xJ,           D.dirXJProjTruthPure,               "TRUTH", true);

        // -------------------------------------------------------------------------
        //  Overlays (shape), integrated over alpha
        //
        // Output folders (per rKey) under:
        //   RecoilJetQA/JES3/<rKey>/xJ_fromJES3/
        //     RECO_vs_TRUTHrecoConditioned/                      (legacy: jet-matched truth)
        //     RECO_vs_TRUTHrecoConditioned_noJetMatch/
        //     TRUTHrecoConditioned_noJetMatch_vs_jetMatched/
        //     RECO_vs_RECO_truthPhoTagged/
        //     RECO_truthPhoTagged_vs_truthTaggedPhoJet/
        //
        // Each folder contains:
        //   perPtBin/overlay_pTbinX.png
        //   table3x3_overlay_shape(_pageY).png
        // -------------------------------------------------------------------------

        auto MakeOverlayShape_TH3xJ =
          [&](TH3* hBlack, TH3* hRed,
              const string& ovTag,
              const string& legBlack,
              const string& legRed,
              const vector<string>& headerLines)
        {
          cout << ANSI_BOLD_CYN
               << "\n[JES3 overlay 2-way] ovTag=" << ovTag
               << "  rKey=" << rKey
               << "  ds=" << ds.label
               << ANSI_RESET << "\n";
          cout << "  hBlack=" << (hBlack ? "FOUND" : "MISSING")
               << "  hRed=" << (hRed ? "FOUND" : "MISSING") << "\n";
          if (hBlack)
          {
            cout << "    hBlack name=" << hBlack->GetName()
                 << "  entries=" << hBlack->GetEntries()
                 << "  nPt=" << hBlack->GetXaxis()->GetNbins() << "\n";
          }
          if (hRed)
          {
            cout << "    hRed   name=" << hRed->GetName()
                 << "  entries=" << hRed->GetEntries()
                 << "  nPt=" << hRed->GetXaxis()->GetNbins() << "\n";
          }
          if (!hBlack || !hRed)
          {
            cout << ANSI_BOLD_YEL
                 << "  [SKIP] 2-way overlay missing source TH3"
                 << ANSI_RESET << "\n";
            return;
          }

            const string dirOvBase = JoinPath(D.dirXJProjOverlay, ovTag);
            const string dirOvPer  = JoinPath(dirOvBase, "perPtBin");
            EnsureDir(D.dirXJProjOverlay);
            EnsureDir(dirOvBase);
            EnsureDir(dirOvPer);

          const int nA = hBlack->GetXaxis()->GetNbins();
          const int nB = hRed->GetXaxis()->GetNbins();
          const int nPtOv = std::min(nA, nB);

          // --- Per pT bin overlay PNGs (shape overlay for fair comparison) ---
          for (int ib = 1; ib <= nPtOv; ++ib)
          {
            TH1* hA = ProjectY_AtXbin_TH3(hBlack, ib,
              TString::Format("ov_%s_blk_%s_%d", ovTag.c_str(), rKey.c_str(), ib).Data()
            );
            TH1* hB = ProjectY_AtXbin_TH3(hRed, ib,
              TString::Format("ov_%s_red_%s_%d", ovTag.c_str(), rKey.c_str(), ib).Data()
            );

            if (hA) { hA->SetDirectory(nullptr); EnsureSumw2(hA); }
            if (hB) { hB->SetDirectory(nullptr); EnsureSumw2(hB); }

              if (!hA || !hB || (hA->GetEntries() <= 0.0 && hB->GetEntries() <= 0.0))
              {
                cout << ANSI_BOLD_YEL
                     << "  [SKIP perPtBin] ovTag=" << ovTag
                     << "  rKey=" << rKey
                     << "  ib=" << ib
                     << "  hA=" << (hA ? "FOUND" : "MISSING")
                     << "  hB=" << (hB ? "FOUND" : "MISSING");
                if (hA) cout << "  entriesA=" << hA->GetEntries();
                if (hB) cout << "  entriesB=" << hB->GetEntries();
                cout << ANSI_RESET << "\n";

                if (hA) delete hA;
                if (hB) delete hB;
                continue;
              }

            NormalizeToUnitArea(hA);
            NormalizeToUnitArea(hB);

              // Style: default black vs red, but override for specific legend labels
              int colA = 1;  // default "black"
              int colB = 2;  // default "red"

              // Force special colors by legend label (applies no matter which histogram is A/B)
              if (legBlack == "Reco (#gamma^{truth} + jet^{truth} tagged)") colA = kViolet + 1;  // purple
              if (legBlack == "Truth (#gamma^{reco} + jet^{reco} tagged)")
              {
                // For RECO_vs_TRUTHrecoConditioned we want Truth drawn in BLUE (open circles).
                // Keep PINK for other overlay folders that intentionally use pink.
                colA = (ovTag == "RECO_vs_TRUTHrecoConditioned") ? (kBlue + 1) : (kPink + 7);
              }
              if (legBlack == "Truth (PURE)" || legBlack == "uncond truth") colA = kBlue + 1;  // blue open circles for unconditioned truth

              if (legRed  == "Reco (#gamma^{truth} + jet^{truth} tagged)") colB = kViolet + 1;  // purple
              if (legRed  == "Truth (#gamma^{reco} + jet^{reco} tagged)")
              {
                colB = (ovTag == "RECO_vs_TRUTHrecoConditioned") ? (kBlue + 1) : (kPink + 7);
              }
              if (legRed == "Truth (PURE)" || legRed == "uncond truth") colB = kBlue + 1;    // blue open circles for unconditioned truth

              hA->SetLineWidth(2);
              hA->SetMarkerStyle(24);
              hA->SetMarkerSize(1.00);
              hA->SetLineColor(colA);
              hA->SetMarkerColor(colA);

              hB->SetLineWidth(2);
              hB->SetMarkerStyle(20);
              hB->SetMarkerSize(1.00);
              hB->SetLineColor(colB);
              hB->SetMarkerColor(colB);

            const string ptLab = AxisBinLabel(hBlack->GetXaxis(), ib, "GeV", 0);

            const string outPng = JoinPath(dirOvPer, TString::Format("overlay_pTbin%d.png", ib).Data());

            TCanvas c(
              TString::Format("c_ov_%s_%s_b%d", ds.label.c_str(), rKey.c_str(), ib).Data(),
              "c_ov", 900, 700
            );
            ApplyCanvasMargins1D(c);

            hA->SetTitle("");
            hA->GetXaxis()->SetTitle("x_{J#gamma}");
            hA->GetYaxis()->SetTitle("A.U.");

            const double ymax = std::max(hA->GetMaximum(), hB->GetMaximum());
            hA->SetMaximum(ymax * 1.25);

            hA->Draw("E1");
            hB->Draw("E1 same");

            const std::string bbLabel = B2BLabel();
            const double jetMinPt     = static_cast<double>(kJetPtMin);
            const double vzCut        = vzCutCm;

            const bool isThisSimRecoVsRecoTruthTaggedPhoJet =
                (ds.isSim && (ovTag == "RECO_vs_RECO_truthTaggedPhoJet"));

            if (isThisSimRecoVsRecoTruthTaggedPhoJet)
            {
                // Match the 2x3-table pad look (and clamp x-range)
                hA->GetXaxis()->SetRangeUser(0.0, 2.0);

                // Title + pT label (same placement style as table)
                {
                  TLatex tt;
                  tt.SetNDC(true);
                  tt.SetTextFont(42);

                  // Main title (top-center)
                  tt.SetTextAlign(22);
                  tt.SetTextSize(0.050);
                  tt.DrawLatex(0.52, 0.95,
                    TString::Format("Photon+Jet 10 and 20 Combined Sim (R = %.1f)", R).Data()
                  );

                  // pT label (upper-left)
                  tt.SetTextAlign(13);
                  tt.SetTextSize(0.043);
                  tt.DrawLatex(0.18, 0.89,
                    TString::Format("p_{T}^{#gamma} = %s", ptLab.c_str()).Data()
                  );
                }

                // Legend (top-center-ish, consistent with table readability)
                TLegend leg(0.40, 0.78, 0.72, 0.90);
                leg.SetTextFont(42);
                leg.SetTextSize(0.045);
                leg.SetFillStyle(0);
                leg.SetBorderSize(0);
                leg.AddEntry(hA, legBlack.c_str(), "ep");
                leg.AddEntry(hB, legRed.c_str(),   "ep");
                leg.Draw();

                // Cut text block on RHS middle (like DATA vs SIM overlay example / your table pads)
                {
                  TLatex tCuts;
                  tCuts.SetNDC(true);
                  tCuts.SetTextFont(42);
                  tCuts.SetTextAlign(33);   // right-top anchored
                  tCuts.SetTextSize(0.045);

                  const double tx = 0.92;
                  double ty = 0.63;
                  const double dY = 0.060;

                  tCuts.DrawLatex(tx, ty,
                    TString::Format("|#Delta#phi(#gamma,jet)| > %s", bbLabel.c_str()).Data()
                  );
                  ty -= dY;
                  tCuts.DrawLatex(tx, ty,
                    TString::Format("p_{T}^{jet} > %.0f GeV", jetMinPt).Data()
                  );
                  ty -= dY;
                  tCuts.DrawLatex(tx, ty,
                    TString::Format("|v_{z}| < %.0f cm", std::fabs(vzCut)).Data()
                  );
                }
              }
              else
              {
                // Default: keep existing perPtBin styling for all other overlay folders

                TLegend leg(0.55, 0.75, 0.85, 0.90);
                leg.SetTextFont(42);
                leg.SetTextSize(0.038);
                leg.SetFillStyle(0);
                leg.SetBorderSize(0);
                leg.AddEntry(hA, legBlack.c_str(), "ep");
                leg.AddEntry(hB, legRed.c_str(),   "ep");
                leg.Draw();

                std::string titleLine;
                if (ovTag == "RECO_vs_RECO_truthTaggedPhoJet")
                {
                  titleLine = "RECO vs RECO (#gamma+jet truth tagged)";
                }
                else
                {
                  titleLine = headerLines.empty() ? ovTag : headerLines.front();
                  const std::string pref = "Overlay (shape): ";
                  if (titleLine.rfind(pref, 0) == 0) titleLine = titleLine.substr(pref.size());
                }

                const string ptLabNoUnit = AxisBinLabel(hBlack->GetXaxis(), ib, "", 0);

                TLatex t;
                t.SetNDC(true);
                t.SetTextFont(42);
                t.SetTextAlign(13);
                t.SetTextSize(0.036);
                t.DrawLatex(0.14, 0.9, titleLine.c_str());

                t.SetTextSize(0.032);
                t.DrawLatex(0.14, 0.845,
                  TString::Format("p_{T}^{#gamma}: %s GeV (R = %.1f)", ptLabNoUnit.c_str(), R).Data()
                );
                t.DrawLatex(0.14, 0.795,
                  TString::Format("|#Delta#phi(#gamma,jet)| > %s", bbLabel.c_str()).Data()
                );
                t.DrawLatex(0.14, 0.755,
                  TString::Format("p_{T}^{jet} > %.0f GeV", jetMinPt).Data()
                );
              }

              SaveCanvas(c, outPng);
              cout << ANSI_BOLD_GRN
                   << "  [WROTE perPtBin] " << outPng
                   << ANSI_RESET << "\n";

              delete hA;
              delete hB;
          }

            // --- 3x3 table page of overlays (shape): all available pT bins (up to 9) ---
            const int nCols   = 3;
            const int nRows   = 3;
            const int perPage = nCols * nRows;  // 9

            const int startBin = 1;

            if (startBin <= nPtOv)
            {
              const int page = 1;

              TCanvas c(
                TString::Format("c_tbl_ov_%s_%s_p%d", ovTag.c_str(), rKey.c_str(), page).Data(),
                "c_tbl_ov", 1500, 1200
              );
              c.Divide(nCols, nRows, 0.001, 0.001);

              vector<TH1*> keep;
              keep.reserve(3 * perPage);

              const int nThisPage = std::min(perPage, nPtOv - startBin + 1);
              for (int k = 0; k < nThisPage; ++k)
              {
                const int ib = startBin + k;
                c.cd(k+1);

                gPad->SetLeftMargin(0.14);
                gPad->SetRightMargin(0.05);
                gPad->SetBottomMargin(0.14);
                gPad->SetTopMargin(0.10);

                TH1* hA = ProjectY_AtXbin_TH3(hBlack, ib,
                  TString::Format("tbl_%s_blk_%s_%d_p%d", ovTag.c_str(), rKey.c_str(), ib, page).Data()
                );
                TH1* hB = ProjectY_AtXbin_TH3(hRed, ib,
                  TString::Format("tbl_%s_red_%s_%d_p%d", ovTag.c_str(), rKey.c_str(), ib, page).Data()
                );

                if (hA) { hA->SetDirectory(nullptr); EnsureSumw2(hA); }
                if (hB) { hB->SetDirectory(nullptr); EnsureSumw2(hB); }

                if (!hA || !hB || (hA->GetEntries() <= 0.0 && hB->GetEntries() <= 0.0))
                {
                  if (hA) delete hA;
                  if (hB) delete hB;
                  TLatex t;
                  t.SetNDC(true);
                  t.SetTextFont(42);
                  t.SetTextSize(0.06);
                  t.DrawLatex(0.15, 0.55, "MISSING");
                  continue;
                }

                NormalizeToUnitArea(hA);
                NormalizeToUnitArea(hB);

                  // Style: default black vs red, but override for specific legend labels
                  int colA = 1;  // default "black"
                  int colB = 2;  // default "red"

                  // Force special colors by legend label (applies no matter which histogram is A/B)
                  if (legBlack == "Reco (#gamma^{truth} + jet^{truth} tagged)") colA = kViolet + 1;  // purple
                  if (legBlack == "Truth (#gamma^{reco} + jet^{reco} tagged)")
                  {
                    // For RECO_vs_TRUTHrecoConditioned we want Truth drawn in BLUE (open circles).
                    // Keep PINK for other overlay folders that intentionally use pink.
                    colA = (ovTag == "RECO_vs_TRUTHrecoConditioned") ? (kBlue + 1) : (kPink + 7);
                  }
                  if (legBlack == "Truth (PURE)" || legBlack == "uncond truth") colA = kBlue + 1;  // unconditioned truth in blue

                  if (legRed  == "Reco (#gamma^{truth} + jet^{truth} tagged)") colB = kViolet + 1;  // purple
                  if (legRed  == "Truth (#gamma^{reco} + jet^{reco} tagged)")
                  {
                    colB = (ovTag == "RECO_vs_TRUTHrecoConditioned") ? (kBlue + 1) : (kPink + 7);
                  }
                  if (legRed == "Truth (PURE)" || legRed == "uncond truth") colB = kBlue + 1;       // unconditioned truth in blue

                  hA->SetLineWidth(2);
                  hA->SetMarkerStyle(24);
                  hA->SetMarkerSize(0.95);
                  hA->SetLineColor(colA);
                  hA->SetMarkerColor(colA);

                  hB->SetLineWidth(2);
                  hB->SetMarkerStyle(20);
                  hB->SetMarkerSize(0.95);
                  hB->SetLineColor(colB);
                  hB->SetMarkerColor(colB);

                const double ymax = std::max(hA->GetMaximum(), hB->GetMaximum());
                hA->SetMaximum(ymax * 1.25);

                  hA->SetTitle("");
                  hA->GetXaxis()->SetTitle("x_{J#gamma}");
                  hA->GetYaxis()->SetTitle("A.U.");
                  if (ds.isSim && (ovTag == "RECO_vs_RECO_truthTaggedPhoJet"))
                  {
                    hA->GetXaxis()->SetRangeUser(0.0, 2.0);
                  }
                  hA->Draw("E1");
                  hB->Draw("E1 same");

                const string ptLab = AxisBinLabel(hBlack->GetXaxis(), ib, "GeV", 0);

                  const std::string bbLabel = B2BLabel();
                  const double jetMinPt     = static_cast<double>(kJetPtMin);
                  const double vzCut        = vzCutCm;

                  const bool isThisSimRecoVsRecoTruthTaggedPhoJet =
                    (ds.isSim && (ovTag == "RECO_vs_RECO_truthTaggedPhoJet"));

                  // Title + pT label
                  {
                    TLatex tt;
                    tt.SetNDC(true);
                    tt.SetTextFont(42);

                      if (isThisSimRecoVsRecoTruthTaggedPhoJet)
                      {
                        // Main title (top-center)
                        tt.SetTextAlign(22);
                        tt.SetTextSize(0.050);
                        tt.DrawLatex(0.52, 0.95,
                          TString::Format("Photon+Jet 10 and 20 Combined Sim (R = %.1f)", R).Data()
                        );

                        // pT label (upper-left, smaller & snug in corner)
                        tt.SetTextAlign(13);
                        tt.SetTextSize(0.043);
                        tt.DrawLatex(0.18, 0.89,
                          TString::Format("p_{T}^{#gamma} = %s", ptLab.c_str()).Data()
                        );
                      }
                      else
                      {
                        // Default behavior
                        tt.SetTextAlign(22);
                        tt.SetTextSize(0.060);
                        tt.DrawLatex(0.52, 0.95,
                          TString::Format("p_{T}^{#gamma} = %s  (R=%.1f)", ptLab.c_str(), R).Data()
                        );
                      }
                  }

                  // Legend placement: top-right, but protect long labels in table pads
                  double lx1 = 0.52, ly1 = 0.74, lx2 = 0.92, ly2 = 0.90;
                  double legTextSize = 0.055;
                  double legMargin   = 0.25;   // fraction of box reserved for markers/lines

                  if (isThisSimRecoVsRecoTruthTaggedPhoJet)
                  {
                    // Snug legend in the top-right and tighten vertically
                    lx1 = 0.87; ly1 = 0.78; lx2 = 0.99; ly2 = 0.9;
                    legTextSize = 0.048;
                    legMargin   = 0.20;
                  }

                  // If either legend entry is long, shift the whole legend LEFT and give it more width
                  const size_t maxLegLen = (legBlack.size() > legRed.size()) ? legBlack.size() : legRed.size();

                  if (maxLegLen >= 28)
                  {
                    lx1 = 0.44;  lx2 = 0.92;   // shift left, but not too far
                    legTextSize = 0.045;
                    legMargin   = 0.18;        // more room for text
                  }
                  if (maxLegLen >= 40)
                  {
                    lx1 = 0.40;  lx2 = 0.92;   // extra shift for very long labels
                    legTextSize = 0.040;
                    legMargin   = 0.16;
                  }

                  TLegend leg(lx1, ly1, lx2, ly2);
                  leg.SetTextFont(42);
                  leg.SetTextSize(legTextSize);
                  leg.SetFillStyle(0);
                  leg.SetBorderSize(0);
                  leg.SetMargin(legMargin);
                  leg.SetEntrySeparation(0.08);
                  leg.AddEntry(hA, legBlack.c_str(), "ep");
                  leg.AddEntry(hB, legRed.c_str(),   "ep");
                  leg.DrawClone();

                  // Cut text block: under legend, right-middle (presentation-ready)
                  {
                    TLatex tCuts;
                    tCuts.SetNDC(true);
                    tCuts.SetTextFont(42);
                    tCuts.SetTextAlign(33);   // right-top anchored
                    tCuts.SetTextSize(0.045);

                    const double tx = lx2 - 0.02;  // anchor near legend right edge
                    double ty = ly1 - 0.1;        // start just below legend
                    const double dY = 0.060;

                    tCuts.DrawLatex(tx, ty,
                      TString::Format("|#Delta#phi(#gamma,jet)| > %s", bbLabel.c_str()).Data()
                    );
                    ty -= dY;
                    tCuts.DrawLatex(tx, ty,
                      TString::Format("p_{T}^{jet} > %.0f GeV", jetMinPt).Data()
                    );
                    ty -= dY;
                    tCuts.DrawLatex(tx, ty,
                      TString::Format("|v_{z}| < %.0f cm", std::fabs(vzCut)).Data()
                    );
                }

                keep.push_back(hA);
                keep.push_back(hB);
              }

              const string outName = "table3x3_overlay_shape.png";

                SaveCanvas(c, JoinPath(dirOvBase, outName));
                cout << ANSI_BOLD_GRN
                     << "  [WROTE table] " << JoinPath(dirOvBase, outName)
                     << ANSI_RESET << "\n";

                for (auto* h : keep) delete h;
            }
        };

        // Keep your legacy overlay folder name (this is jet-matched truth in your code):
        MakeOverlayShape_TH3xJ(
          H.hTrut_xJ, H.hReco_xJ,
          "RECO_vs_TRUTHrecoConditioned",
          "Truth (#gamma^{reco} + jet^{reco} tagged)",
          "Reco",
          {"Overlay (shape): RECO vs TRUTH reco-conditioned", "Truth is jet-matched (ΔR<=0.3)"}
        );

        // PURE overlay (no reco conditioning on truth; reco is unconditional JES3 reco)
        MakeOverlayShape_TH3xJ(
          H.hTrutPure_xJ, H.hReco_xJ,
          "RECO_vs_TRUTHpure",
          "uncond truth",
          "uncond reco",
          {"unconditioned truth versus reco"}
        );


        // New overlays covering the missing JES3 variants:
        MakeOverlayShape_TH3xJ(
          H.hTrutNoJM_xJ, H.hReco_xJ,
          "RECO_vs_TRUTHrecoConditioned_noJetMatch",
          "Truth (reco-cond, NO jet match)",
          "Reco",
          {"Overlay (shape): RECO vs TRUTH reco-conditioned", "Truth is NO-jet-match"}
        );

        MakeOverlayShape_TH3xJ(
          H.hTrutNoJM_xJ, H.hTrut_xJ,
          "TRUTHrecoConditioned_noJetMatch_vs_jetMatched",
          "Truth (NO jet match)",
          "Truth (jet-matched)",
          {"Overlay (shape): TRUTH reco-conditioned comparison"}
        );

        MakeOverlayShape_TH3xJ(
          H.hRecoTruthPhoTagged_xJ, H.hReco_xJ,
          "RECO_vs_RECO_truthPhoTagged",
          "Reco (#gamma^{truth} tagged)",
          "Reco (all)",
          {"Overlay (shape): RECO baseline vs truth-#gamma tagged"}
        );

        MakeOverlayShape_TH3xJ(
          H.hRecoTruthPhoTagged_xJ, H.hRecoTruthTagged_xJ,
          "RECO_truthPhoTagged_vs_truthTaggedPhoJet",
          "Reco (#gamma^{truth} tagged)",
          "Reco (#gamma^{truth} + truth jet)",
          {"Overlay (shape): RECO truth-tagged subsets"}
        );
        
        // NEW: Reco (all) vs Reco (truth-#gamma + jet1 matched)
        MakeOverlayShape_TH3xJ(
          H.hRecoTruthTagged_xJ,  H.hReco_xJ,
          "RECO_vs_RECO_truthTaggedPhoJet",
          "Reco (#gamma^{truth} + jet^{truth} tagged)",
          "Reco",
          {"Overlay (shape): RECO vs doubly truth-tagged RECO", "Reco is jet-matched (#DeltaR<=0.3)"}
        );

        // NEW: RECO_vs_RECO_truthTaggedPhoJet + DATA overlay
        //   Three curves per pad: SIM reco (truth-tagged), SIM reco (all), DATA reco
        //   Only active when isSimAndDataPP = true (uses the default PP data file).
        if (isSimAndDataPP && H.hRecoTruthTagged_xJ && H.hReco_xJ)
        {
             static TFile* fPPdataOv = nullptr;
             static TDirectory* dirPPdataOv = nullptr;

             if (!fPPdataOv)
             {
               fPPdataOv = TFile::Open(InputPP().c_str(), "READ");
               if (fPPdataOv)
               {
                 dirPPdataOv = fPPdataOv->GetDirectory(kTriggerPP.c_str());
                 if (!dirPPdataOv) dirPPdataOv = fPPdataOv;
               }
             }

             TH3* hDataReco_xJ = nullptr;
             if (fPPdataOv && dirPPdataOv)
             {
               hDataReco_xJ = dynamic_cast<TH3*>(dirPPdataOv->Get(
                 TString::Format("h_JES3_pT_xJ_alpha_%s", rKey.c_str()).Data()
               ));
             }

             if (hDataReco_xJ)
             {
               const string ovTagData = "RECO_vs_RECO_truthTaggedPhoJet_data";
               const string dirOvBaseData = JoinPath(D.dirXJProjOverlay, ovTagData);
               const string dirOvPerData  = JoinPath(dirOvBaseData, "perPtBin");
               EnsureDir(D.dirXJProjOverlay);
               EnsureDir(dirOvBaseData);
               EnsureDir(dirOvPerData);

               cout << ANSI_BOLD_CYN
                    << "\n[JES3 overlay 3-way+data] ovTag=" << ovTagData
                    << "  rKey=" << rKey
                    << "  ds=" << ds.label
                    << ANSI_RESET << "\n";
               cout << "  hRecoTruthTagged=" << (H.hRecoTruthTagged_xJ ? "FOUND" : "MISSING")
                    << "  hReco=" << (H.hReco_xJ ? "FOUND" : "MISSING")
                    << "  hDataReco=" << (hDataReco_xJ ? "FOUND" : "MISSING") << "\n";

               const int nPtSim  = std::min(H.hRecoTruthTagged_xJ->GetXaxis()->GetNbins(),
                                            H.hReco_xJ->GetXaxis()->GetNbins());
               const int nPtData = hDataReco_xJ->GetXaxis()->GetNbins();
               const int nPtOvD  = std::min(nPtSim, nPtData);

               const std::string legSimTruthTag = "Sim Reco (#gamma^{truth} + jet^{truth})";
               const std::string legSimReco     = "Sim Reco";
               const std::string legData        = "Run24pp";

               // ---- Per pT bin overlay PNGs ----
               for (int ib = 1; ib <= nPtOvD; ++ib)
               {
                 TH1* hA = ProjectY_AtXbin_TH3(H.hRecoTruthTagged_xJ, ib,
                   TString::Format("ov3d_%s_simTT_%s_%d", ovTagData.c_str(), rKey.c_str(), ib).Data()
                 );
                 TH1* hB = ProjectY_AtXbin_TH3(H.hReco_xJ, ib,
                   TString::Format("ov3d_%s_simR_%s_%d", ovTagData.c_str(), rKey.c_str(), ib).Data()
                 );
                 TH1* hD = ProjectY_AtXbin_TH3(hDataReco_xJ, ib,
                   TString::Format("ov3d_%s_data_%s_%d", ovTagData.c_str(), rKey.c_str(), ib).Data()
                 );

                 if (hA) { hA->SetDirectory(nullptr); EnsureSumw2(hA); }
                 if (hB) { hB->SetDirectory(nullptr); EnsureSumw2(hB); }
                 if (hD) { hD->SetDirectory(nullptr); EnsureSumw2(hD); }

                 if (!hA || !hB || !hD ||
                     (hA->GetEntries() <= 0.0 && hB->GetEntries() <= 0.0 && hD->GetEntries() <= 0.0))
                 {
                   if (hA) delete hA;
                   if (hB) delete hB;
                   if (hD) delete hD;
                   continue;
                 }

                 NormalizeToUnitArea(hA);
                 NormalizeToUnitArea(hB);
                 NormalizeToUnitArea(hD);

                 hA->SetLineWidth(2);
                 hA->SetMarkerStyle(24);
                 hA->SetMarkerSize(1.00);
                 hA->SetLineColor(kViolet + 1);
                 hA->SetMarkerColor(kViolet + 1);

                 hB->SetLineWidth(2);
                 hB->SetMarkerStyle(20);
                 hB->SetMarkerSize(1.00);
                 hB->SetLineColor(kRed);
                 hB->SetMarkerColor(kRed);

                 hD->SetLineWidth(2);
                 hD->SetMarkerStyle(21);
                 hD->SetMarkerSize(1.00);
                 hD->SetLineColor(kBlack);
                 hD->SetMarkerColor(kBlack);

                 const string ptLab = AxisBinLabel(H.hReco_xJ->GetXaxis(), ib, "GeV", 0);

                 const string outPng = JoinPath(dirOvPerData, TString::Format("overlay_pTbin%d.png", ib).Data());

                 TCanvas c(
                   TString::Format("c_ov3d_%s_%s_b%d", ds.label.c_str(), rKey.c_str(), ib).Data(),
                   "c_ov3d", 900, 700
                 );
                 ApplyCanvasMargins1D(c);

                 const double ymax = std::max(std::max(hA->GetMaximum(), hB->GetMaximum()), hD->GetMaximum());

                 hA->SetTitle("");
                 hA->GetXaxis()->SetTitle("x_{J#gamma}");
                 hA->GetYaxis()->SetTitle("A.U.");
                 hA->GetXaxis()->SetRangeUser(0.0, 2.0);
                 hA->SetMaximum(ymax * 1.25);

                 hA->Draw("E1");
                 hB->Draw("E1 same");
                 hD->Draw("E1 same");

                 {
                   TLatex tt;
                   tt.SetNDC(true);
                   tt.SetTextFont(42);

                   tt.SetTextAlign(22);
                   tt.SetTextSize(0.050);
                   tt.DrawLatex(0.52, 0.95,
                     TString::Format("Photon+Jet 10 and 20 Combined Sim + Data (R = %.1f)", R).Data()
                   );

                   tt.SetTextAlign(13);
                   tt.SetTextSize(0.043);
                   tt.DrawLatex(0.18, 0.89,
                     TString::Format("p_{T}^{#gamma} = %s", ptLab.c_str()).Data()
                   );
                 }

                 TLegend leg(0.55, 0.73, 0.92, 0.90);
                 leg.SetTextFont(42);
                 leg.SetTextSize(0.030);
                 leg.SetFillStyle(0);
                 leg.SetBorderSize(0);
                 leg.AddEntry(hA, legSimTruthTag.c_str(), "ep");
                 leg.AddEntry(hB, legSimReco.c_str(),     "ep");
                 leg.AddEntry(hD, legData.c_str(),        "ep");
                 leg.Draw();

                 {
                   const std::string bbLabelD = B2BLabel();
                   const double jetMinPtD     = static_cast<double>(kJetPtMin);

                   TLatex tCuts;
                   tCuts.SetNDC(true);
                   tCuts.SetTextFont(42);
                   tCuts.SetTextAlign(33);
                   tCuts.SetTextSize(0.045);

                   const double tx = 0.92;
                   double ty = 0.63;
                   const double dY = 0.060;

                   tCuts.DrawLatex(tx, ty,
                     TString::Format("|#Delta#phi(#gamma,jet)| > %s", bbLabelD.c_str()).Data()
                   );
                   ty -= dY;
                   tCuts.DrawLatex(tx, ty,
                     TString::Format("p_{T}^{jet} > %.0f GeV", jetMinPtD).Data()
                   );
                   ty -= dY;
                   tCuts.DrawLatex(tx, ty,
                     TString::Format("|v_{z}| < %.0f cm", std::fabs(vzCutCm)).Data()
                   );
                 }

                 SaveCanvas(c, outPng);
                 cout << ANSI_BOLD_GRN
                      << "  [WROTE perPtBin+data] " << outPng
                      << ANSI_RESET << "\n";

                 delete hA;
                 delete hB;
                 delete hD;
               }

               // ---- 3x3 table overlay (shape) ----
               {
                 const int nCols   = 3;
                 const int nRows   = 3;
                 const int perPage = nCols * nRows;

                 if (nPtOvD > 0)
                 {
                   TCanvas c(
                     TString::Format("c_tbl_ov3d_%s_%s", ovTagData.c_str(), rKey.c_str()).Data(),
                     "c_tbl_ov3d", 1500, 1200
                   );
                   c.Divide(nCols, nRows, 0.001, 0.001);

                   vector<TH1*> keep;
                   keep.reserve(3 * perPage);

                   const int nThisPage = std::min(perPage, nPtOvD);
                   for (int k = 0; k < nThisPage; ++k)
                   {
                     const int ib = 1 + k;
                     c.cd(k+1);

                     gPad->SetLeftMargin(0.14);
                     gPad->SetRightMargin(0.05);
                     gPad->SetBottomMargin(0.14);
                     gPad->SetTopMargin(0.10);

                     TH1* hA = ProjectY_AtXbin_TH3(H.hRecoTruthTagged_xJ, ib,
                       TString::Format("tbl3d_%s_simTT_%s_%d", ovTagData.c_str(), rKey.c_str(), ib).Data()
                     );
                     TH1* hB = ProjectY_AtXbin_TH3(H.hReco_xJ, ib,
                       TString::Format("tbl3d_%s_simR_%s_%d", ovTagData.c_str(), rKey.c_str(), ib).Data()
                     );
                     TH1* hD = ProjectY_AtXbin_TH3(hDataReco_xJ, ib,
                       TString::Format("tbl3d_%s_data_%s_%d", ovTagData.c_str(), rKey.c_str(), ib).Data()
                     );

                     if (hA) { hA->SetDirectory(nullptr); EnsureSumw2(hA); }
                     if (hB) { hB->SetDirectory(nullptr); EnsureSumw2(hB); }
                     if (hD) { hD->SetDirectory(nullptr); EnsureSumw2(hD); }

                     if (!hA || !hB || !hD ||
                         (hA->GetEntries() <= 0.0 && hB->GetEntries() <= 0.0 && hD->GetEntries() <= 0.0))
                     {
                       if (hA) delete hA;
                       if (hB) delete hB;
                       if (hD) delete hD;
                       TLatex t;
                       t.SetNDC(true);
                       t.SetTextFont(42);
                       t.SetTextSize(0.06);
                       t.DrawLatex(0.15, 0.55, "MISSING");
                       continue;
                     }

                     NormalizeToUnitArea(hA);
                     NormalizeToUnitArea(hB);
                     NormalizeToUnitArea(hD);

                     hA->SetLineWidth(2);
                     hA->SetMarkerStyle(24);
                     hA->SetMarkerSize(0.95);
                     hA->SetLineColor(kViolet + 1);
                     hA->SetMarkerColor(kViolet + 1);

                     hB->SetLineWidth(2);
                     hB->SetMarkerStyle(20);
                     hB->SetMarkerSize(0.95);
                     hB->SetLineColor(kRed);
                     hB->SetMarkerColor(kRed);

                     hD->SetLineWidth(2);
                     hD->SetMarkerStyle(21);
                     hD->SetMarkerSize(0.95);
                     hD->SetLineColor(kBlack);
                     hD->SetMarkerColor(kBlack);

                     const double ymax = std::max(std::max(hA->GetMaximum(), hB->GetMaximum()), hD->GetMaximum());
                     hA->SetMaximum(ymax * 1.25);

                     hA->SetTitle("");
                     hA->GetXaxis()->SetTitle("x_{J#gamma}");
                     hA->GetYaxis()->SetTitle("A.U.");
                     hA->GetXaxis()->SetRangeUser(0.0, 2.0);
                     hA->Draw("E1");
                     hB->Draw("E1 same");
                     hD->Draw("E1 same");

                     const string ptLab = AxisBinLabel(H.hReco_xJ->GetXaxis(), ib, "GeV", 0);

                     {
                       TLatex tt;
                       tt.SetNDC(true);
                       tt.SetTextFont(42);

                       tt.SetTextAlign(22);
                       tt.SetTextSize(0.050);
                       tt.DrawLatex(0.52, 0.95,
                         TString::Format("Photon+Jet 10 and 20 Combined Sim + Data (R = %.1f)", R).Data()
                       );

                       tt.SetTextAlign(13);
                       tt.SetTextSize(0.043);
                       tt.DrawLatex(0.18, 0.89,
                         TString::Format("p_{T}^{#gamma} = %s", ptLab.c_str()).Data()
                       );
                     }

                     TLegend leg(0.87, 0.73, 0.99, 0.9);
                     leg.SetTextFont(42);
                     leg.SetTextSize(0.040);
                     leg.SetFillStyle(0);
                     leg.SetBorderSize(0);
                     leg.SetMargin(0.20);
                     leg.SetEntrySeparation(0.08);
                     leg.AddEntry(hA, legSimTruthTag.c_str(), "ep");
                     leg.AddEntry(hB, legSimReco.c_str(),     "ep");
                     leg.AddEntry(hD, legData.c_str(),        "ep");
                     leg.DrawClone();

                     {
                       const std::string bbLabelD = B2BLabel();
                       const double jetMinPtD     = static_cast<double>(kJetPtMin);

                       TLatex tCuts;
                       tCuts.SetNDC(true);
                       tCuts.SetTextFont(42);
                       tCuts.SetTextAlign(33);
                       tCuts.SetTextSize(0.045);

                       const double tx = 0.85;
                       double ty = 0.63;
                       const double dY = 0.060;

                       tCuts.DrawLatex(tx, ty,
                         TString::Format("|#Delta#phi(#gamma,jet)| > %s", bbLabelD.c_str()).Data()
                       );
                       ty -= dY;
                       tCuts.DrawLatex(tx, ty,
                         TString::Format("p_{T}^{jet} > %.0f GeV", jetMinPtD).Data()
                       );
                       ty -= dY;
                       tCuts.DrawLatex(tx, ty,
                         TString::Format("|v_{z}| < %.0f cm", std::fabs(vzCutCm)).Data()
                       );
                     }

                     keep.push_back(hA);
                     keep.push_back(hB);
                     keep.push_back(hD);
                   }

                   const string outName = "table3x3_overlay_shape.png";
                   SaveCanvas(c, JoinPath(dirOvBaseData, outName));
                   cout << ANSI_BOLD_GRN
                        << "  [WROTE table+data] " << JoinPath(dirOvBaseData, outName)
                        << ANSI_RESET << "\n";

                   for (auto* h : keep) delete h;
                 }
               }
             }
             else
             {
               cout << ANSI_BOLD_YEL
                    << "  [SKIP] RECO_vs_RECO_truthTaggedPhoJet_data: missing DATA TH3 for rKey=" << rKey
                    << ANSI_RESET << "\n";
             }
        }

        // NEW: Reco (truth-#gamma + truth-jet tagged) vs Truth (reco-#gamma + reco-jet tagged)
        //   Purple: Reco (#gamma^{truth} + jet^{truth} tagged)  -> H.hRecoTruthTagged_xJ
        //   Pink  : Truth (#gamma^{reco} + jet^{reco} tagged)  -> H.hTrut_xJ  (your legacy jet-matched truth)
        MakeOverlayShape_TH3xJ(
          H.hRecoTruthTagged_xJ,  H.hTrut_xJ,
          "RECO_truthTaggedPhoJet_vs_TRUTHrecoConditioned",
          "Reco (#gamma^{truth} + jet^{truth} tagged)",
          "Truth (#gamma^{reco} + jet^{reco} tagged)",
          {"Overlay (shape): RECO truth-tagged vs TRUTH reco-conditioned", "Both are jet-matched (#DeltaR<=0.3)"}
        );

        
        auto MakeOverlayShape_TH3xJ_3way =
          [&](TH3* hReco, TH3* hTruth, TH3* hRecoTruth,
              const string& ovTag,
              const string& legReco,
              const string& legTruth,
              const string& legRecoTruth,
              const vector<string>& headerLines)
        {
          cout << ANSI_BOLD_CYN
               << "\n[JES3 overlay 3-way] ovTag=" << ovTag
               << "  rKey=" << rKey
               << "  ds=" << ds.label
               << ANSI_RESET << "\n";
          cout << "  hReco=" << (hReco ? "FOUND" : "MISSING")
               << "  hTruth=" << (hTruth ? "FOUND" : "MISSING")
               << "  hRecoTruth=" << (hRecoTruth ? "FOUND" : "MISSING") << "\n";
          if (hReco)
          {
            cout << "    hReco      name=" << hReco->GetName()
                 << "  entries=" << hReco->GetEntries()
                 << "  nPt=" << hReco->GetXaxis()->GetNbins() << "\n";
          }
          if (hTruth)
          {
            cout << "    hTruth     name=" << hTruth->GetName()
                 << "  entries=" << hTruth->GetEntries()
                 << "  nPt=" << hTruth->GetXaxis()->GetNbins() << "\n";
          }
          if (hRecoTruth)
          {
            cout << "    hRecoTruth name=" << hRecoTruth->GetName()
                 << "  entries=" << hRecoTruth->GetEntries()
                 << "  nPt=" << hRecoTruth->GetXaxis()->GetNbins() << "\n";
          }
          if (!hReco || !hTruth || !hRecoTruth)
          {
            cout << ANSI_BOLD_YEL
                 << "  [SKIP] 3-way overlay missing source TH3"
                 << ANSI_RESET << "\n";
            return;
          }

          // --------------------------------------------------------------------------
          // KINEMATIC FLOOR (requested)
          //   jet pT cut = 5 GeV  =>  x_{Jγ} >= 5 / pTγ(event)
          // For a pTγ BIN, we draw a single representative vertical line at:
          //   x_floor = 5 / <pTγ>   using the bin center.
          // --------------------------------------------------------------------------
          const double jetPtMin_GeV = 5.0;

          const string dirOvRoot = JoinPath(D.dirXJProj, "Overlay");
          const string dirOvBase = JoinPath(dirOvRoot, ovTag);
          const string dirOvPer  = JoinPath(dirOvBase, "perPtBin");
          EnsureDir(dirOvRoot);
          EnsureDir(dirOvBase);
          EnsureDir(dirOvPer);

          const int nA = hReco->GetXaxis()->GetNbins();
          const int nB = hTruth->GetXaxis()->GetNbins();
          const int nC = hRecoTruth->GetXaxis()->GetNbins();
          const int nPtOv = std::min(std::min(nA, nB), nC);

          auto XFloorForPtBin = [&](const TAxis* ax, int ib)->double
          {
            if (!ax) return -1.0;
            const double ptC = ax->GetBinCenter(ib);
            if (!(ptC > 0.0)) return -1.0;
            return jetPtMin_GeV / ptC;
          };

          auto DrawKinematicFloorLine = [&](double xFloor, double yMin, double yMax)
          {
            if (!(xFloor > 0.0) || !(yMax > yMin)) return;
            TLine ln(xFloor, yMin, xFloor, yMax);
            ln.SetLineColor(kGray + 2);
            ln.SetLineStyle(2);
            ln.SetLineWidth(2);
            ln.Draw("same");
          };

          auto MakeSafeRatio = [&](const TH1* hNum, const TH1* hDen, const char* name)->TH1*
          {
            if (!hNum || !hDen) return nullptr;

            TH1* r = (TH1*)hNum->Clone(name);
            if (!r) return nullptr;
            r->SetDirectory(nullptr);
            EnsureSumw2(r);

            r->Divide(hDen); // standard propagation

            // Sanitize bins where denominator is 0 (or non-finite)
            for (int i = 1; i <= r->GetNbinsX(); ++i)
            {
              const double den = hDen->GetBinContent(i);
              if (!(den > 0.0))
              {
                r->SetBinContent(i, 0.0);
                r->SetBinError(i,   0.0);
                continue;
              }

              const double v = r->GetBinContent(i);
              const double e = r->GetBinError(i);
              if (!TMath::Finite(v) || !TMath::Finite(e))
              {
                r->SetBinContent(i, 0.0);
                r->SetBinError(i,   0.0);
              }
            }
            return r;
          };

          // Fixed ratio y-range for stable/consistent visual readout in small pads
          const double ratioYMin = 0.0;
          const double ratioYMax = 2.0;

          // ==========================================================================
          // (A) PER pT BIN: overlay + ratio (same perPtBin/ folder)
          // ==========================================================================
          for (int ib = 1; ib <= nPtOv; ++ib)
          {
            TH1* hR = ProjectY_AtXbin_TH3(hReco, ib,
              TString::Format("ov3_%s_reco_%s_%d", ovTag.c_str(), rKey.c_str(), ib).Data()
            );
            TH1* hT = ProjectY_AtXbin_TH3(hTruth, ib,
              TString::Format("ov3_%s_truth_%s_%d", ovTag.c_str(), rKey.c_str(), ib).Data()
            );
            TH1* hB = ProjectY_AtXbin_TH3(hRecoTruth, ib,
              TString::Format("ov3_%s_recotruth_%s_%d", ovTag.c_str(), rKey.c_str(), ib).Data()
            );

            if (hR) { hR->SetDirectory(nullptr); EnsureSumw2(hR); }
            if (hT) { hT->SetDirectory(nullptr); EnsureSumw2(hT); }
            if (hB) { hB->SetDirectory(nullptr); EnsureSumw2(hB); }

              if (!hR || !hT || !hB || (hR->GetEntries() <= 0.0 && hT->GetEntries() <= 0.0 && hB->GetEntries() <= 0.0))
              {
                cout << ANSI_BOLD_YEL
                     << "  [SKIP perPtBin 3-way] ovTag=" << ovTag
                     << "  rKey=" << rKey
                     << "  ib=" << ib
                     << "  hR=" << (hR ? "FOUND" : "MISSING")
                     << "  hT=" << (hT ? "FOUND" : "MISSING")
                     << "  hB=" << (hB ? "FOUND" : "MISSING");
                if (hR) cout << "  entriesReco=" << hR->GetEntries();
                if (hT) cout << "  entriesTruth=" << hT->GetEntries();
                if (hB) cout << "  entriesRecoTruth=" << hB->GetEntries();
                cout << ANSI_RESET << "\n";

                if (hR) delete hR;
                if (hT) delete hT;
                if (hB) delete hB;
                continue;
              }

            // Normalize to unit area (shape overlays, consistent with your existing output)
            NormalizeToUnitArea(hR);
            NormalizeToUnitArea(hT);
            NormalizeToUnitArea(hB);

            // Style: RECO=black, TRUTH=red, RECO(truth-cond)=blue
            hR->SetLineWidth(2);
            hR->SetMarkerStyle(20);
            hR->SetMarkerSize(1.00);
            hR->SetLineColor(kBlack);
            hR->SetMarkerColor(kBlack);

            hT->SetLineWidth(2);
            hT->SetMarkerStyle(24);
            hT->SetMarkerSize(1.00);
            hT->SetLineColor(kRed);
            hT->SetMarkerColor(kRed);

            hB->SetLineWidth(2);
            hB->SetMarkerStyle(21);
            hB->SetMarkerSize(1.00);
            hB->SetLineColor(kBlue);
            hB->SetMarkerColor(kBlue);

            const double ymax = std::max(std::max(hR->GetMaximum(), hT->GetMaximum()), hB->GetMaximum());
            const double yMaxPlot = ymax * 1.25;
            hR->SetMaximum(yMaxPlot);

            hR->SetTitle("");
            hR->GetXaxis()->SetTitle("x_{J#gamma}");
            hR->GetYaxis()->SetTitle("A.U.");

            const string ptLab  = AxisBinLabel(hReco->GetXaxis(), ib, "GeV", 0);
            const string outPng = JoinPath(dirOvPer, TString::Format("overlay_pTbin%d.png", ib).Data());

            const double xFloor = XFloorForPtBin(hReco->GetXaxis(), ib);

            // ---------------------- Overlay canvas ----------------------
            TCanvas c(
              TString::Format("c_ov3_%s_%s_b%d", ds.label.c_str(), rKey.c_str(), ib).Data(),
              "c_ov3", 900, 700
            );
            ApplyCanvasMargins1D(c);

              hR->Draw("E1");
              hT->Draw("E1 same");
              hB->Draw("E1 same");

              // Kinematic floor line:
              //   KEEP for other overlays, but DO NOT draw for truthTaggedPhoJet overlays.
              if (ovTag.find("truthTaggedPhoJet") == std::string::npos)
              {
                DrawKinematicFloorLine(xFloor, 0.0, yMaxPlot);

                if (xFloor > 0.0)
                {
                  TLatex tFloor;
                  tFloor.SetNDC(true);
                  tFloor.SetTextFont(42);
                  tFloor.SetTextAlign(13);
                  tFloor.SetTextSize(0.034);
                  tFloor.DrawLatex(0.16, 0.74,
                    TString::Format("p_{T}^{jet}>%.0f GeV  #Rightarrow  x_{min}#approx%.3f", jetPtMin_GeV, xFloor).Data()
                  );
                }
              }

              // Legend (top-right)
              double lx1 = 0.50, ly1 = 0.72, lx2 = 0.86, ly2 = 0.90;
              double legTextSize = 0.041;
              if (ovTag.find("truthTaggedPhoJet") != std::string::npos)
              {
                lx1 = 0.55; lx2 = 0.86;
                legTextSize = 0.038;
              }

              TLegend leg(lx1, ly1, lx2, ly2);
              leg.SetTextFont(42);
              leg.SetTextSize(legTextSize);
              leg.SetFillStyle(0);
              leg.SetBorderSize(0);
              leg.AddEntry(hR, legReco.c_str(),      "ep");
              leg.AddEntry(hT, legTruth.c_str(),     "ep");
              leg.AddEntry(hB, legRecoTruth.c_str(), "ep");
              leg.Draw();

              // Plot labeling:
              //   For RECO_vs_RECO_truthTaggedPhoJet overlays:
              //     - Use a single title line (no dataset header block)
              //     - Put cuts on RHS middle (like DATA vs SIM overlay style)
              //     - Do NOT show kinematic-floor line/label
              if (ovTag.find("truthTaggedPhoJet") != std::string::npos)
              {
                TLatex tTitle;
                tTitle.SetNDC(true);
                tTitle.SetTextFont(42);
                tTitle.SetTextAlign(13);
                tTitle.SetTextSize(0.050);
                tTitle.DrawLatex(0.14, 0.92,
                  TString::Format("Reco vs Reco (#gamma+jet truth tagged), p_{T}^{#gamma} = %s, R = %.1f", ptLab.c_str(), R).Data()
                );

                vector<string> cutLines;
                cutLines.push_back(TString::Format("|#Delta#phi(#gamma,jet)| > %s", B2BLabel().c_str()).Data());
                cutLines.push_back(TString::Format("p_{T}^{jet} > %.0f GeV", jetPtMin_GeV).Data());
                cutLines.push_back(TString::Format("|v_{z}| < %.0f cm", vzCutCm).Data());
                DrawLatexLines(0.62, 0.62, cutLines, 0.036, 0.045);
              }
              else
              {
                DrawLatexLines(0.14, 0.92, DefaultHeaderLines(ds), 0.034, 0.045);

                vector<string> lines = headerLines;
                lines.push_back(TString::Format("p_{T}^{#gamma}: %s  (R=%.1f)", ptLab.c_str(), R).Data());
                lines.push_back("Integrated over #alpha");
                DrawLatexLines(0.14, 0.84, lines, 0.030, 0.040);
              }

              SaveCanvas(c, outPng);
              cout << ANSI_BOLD_GRN
                   << "  [WROTE perPtBin 3-way] " << outPng
                   << ANSI_RESET << "\n";

            // ---------------------- Ratio canvas (requested) ----------------------
            // Ratios (shape ratios because inputs are unit-area normalized):
            //   (1) Reco(truth-cond) / Reco(all)
            //   (2) Reco(truth-cond) / Truth(reco-cond)
            TH1* rB_over_R = MakeSafeRatio(hB, hR,
              TString::Format("ratio_%s_BoverR_%s_%d", ovTag.c_str(), rKey.c_str(), ib).Data()
            );
            TH1* rB_over_T = MakeSafeRatio(hB, hT,
              TString::Format("ratio_%s_BoverT_%s_%d", ovTag.c_str(), rKey.c_str(), ib).Data()
            );

            if (rB_over_R && rB_over_T)
            {
              rB_over_R->SetTitle("");
              rB_over_R->GetXaxis()->SetTitle("x_{J#gamma}");
              rB_over_R->GetYaxis()->SetTitle("Ratio (unit-area shapes)");
              rB_over_R->SetMinimum(ratioYMin);
              rB_over_R->SetMaximum(ratioYMax);

              rB_over_R->SetLineWidth(2);
              rB_over_R->SetMarkerStyle(20);
              rB_over_R->SetMarkerSize(1.00);
              rB_over_R->SetLineColor(kBlue);
              rB_over_R->SetMarkerColor(kBlue);

              rB_over_T->SetLineWidth(2);
              rB_over_T->SetMarkerStyle(24);
              rB_over_T->SetMarkerSize(1.00);
              rB_over_T->SetLineColor(kRed);
              rB_over_T->SetMarkerColor(kRed);

              const string outRatioPng = JoinPath(dirOvPer, TString::Format("ratio_pTbin%d.png", ib).Data());

              TCanvas cr(
                TString::Format("c_ratio_%s_%s_b%d", ds.label.c_str(), rKey.c_str(), ib).Data(),
                "c_ratio", 900, 700
              );
              ApplyCanvasMargins1D(cr);

              rB_over_R->Draw("E1");
              rB_over_T->Draw("E1 same");

              // Unity line
              const double xMin = rB_over_R->GetXaxis()->GetXmin();
              const double xMax = rB_over_R->GetXaxis()->GetXmax();
              TLine unity(xMin, 1.0, xMax, 1.0);
              unity.SetLineColor(kGray + 2);
              unity.SetLineStyle(2);
              unity.SetLineWidth(2);
              unity.Draw("same");

              // Kinematic floor line
              DrawKinematicFloorLine(xFloor, ratioYMin, ratioYMax);

              // Legend
              TLegend legr(0.50, 0.72, 0.88, 0.90);
              legr.SetTextFont(42);
              legr.SetTextSize(0.038);
              legr.SetFillStyle(0);
              legr.SetBorderSize(0);
              legr.AddEntry(rB_over_R, "Reco(truth-cond) / Reco(all)", "ep");
              legr.AddEntry(rB_over_T, "Reco(truth-cond) / Truth(reco-cond)", "ep");
              legr.Draw();

              DrawLatexLines(0.14, 0.92, DefaultHeaderLines(ds), 0.034, 0.045);

              vector<string> rlines = headerLines;
              rlines.push_back(TString::Format("p_{T}^{#gamma}: %s  (R=%.1f)", ptLab.c_str(), R).Data());
              rlines.push_back("Ratios of unit-area normalized x_{J#gamma} spectra");
              rlines.push_back("Integrated over #alpha");
              DrawLatexLines(0.14, 0.84, rlines, 0.030, 0.040);

              SaveCanvas(cr, outRatioPng);
            }

            if (rB_over_R) delete rB_over_R;
            if (rB_over_T) delete rB_over_T;

            delete hR;
            delete hT;
            delete hB;
          }

            // ==========================================================================
            // (B) 2x3 TABLE PAGE: overlays (shape) WITH floor line (last 6 pT bins; skip first pT bin)
            // ==========================================================================
            const int nCols   = 3;
            const int nRows   = 2;
            const int perPage = nCols * nRows;  // 6
            int page = 0;

            const int firstBinToTable = 2;  // skip the first pT bin
            const int startBin = std::max(firstBinToTable, nPtOv - perPage + 1);

            if (startBin <= nPtOv)
            {
              ++page;

              TCanvas c(
                TString::Format("c_tbl_ov3_%s_%s_p%d", ovTag.c_str(), rKey.c_str(), page).Data(),
                "c_tbl_ov3", 1500, 900
              );
              c.Divide(nCols, nRows, 0.001, 0.001);

              vector<TH1*> keep;
              keep.reserve(3 * perPage);

              const int nThisPage = std::min(perPage, nPtOv - startBin + 1);
              for (int k = 0; k < nThisPage; ++k)
              {
                const int ib = startBin + k;
                c.cd(k+1);

                gPad->SetLeftMargin(0.14);
                gPad->SetRightMargin(0.05);
                gPad->SetBottomMargin(0.14);
                gPad->SetTopMargin(0.10);

                TH1* hR = ProjectY_AtXbin_TH3(hReco, ib,
                  TString::Format("tbl3_%s_reco_%s_%d_p%d", ovTag.c_str(), rKey.c_str(), ib, page).Data()
                );
                TH1* hT = ProjectY_AtXbin_TH3(hTruth, ib,
                  TString::Format("tbl3_%s_truth_%s_%d_p%d", ovTag.c_str(), rKey.c_str(), ib, page).Data()
                );
                TH1* hB = ProjectY_AtXbin_TH3(hRecoTruth, ib,
                  TString::Format("tbl3_%s_recotruth_%s_%d_p%d", ovTag.c_str(), rKey.c_str(), ib, page).Data()
                );

                if (hR) { hR->SetDirectory(nullptr); EnsureSumw2(hR); }
                if (hT) { hT->SetDirectory(nullptr); EnsureSumw2(hT); }
                if (hB) { hB->SetDirectory(nullptr); EnsureSumw2(hB); }

                if (!hR || !hT || !hB || (hR->GetEntries() <= 0.0 && hT->GetEntries() <= 0.0 && hB->GetEntries() <= 0.0))
                {
                  if (hR) delete hR;
                  if (hT) delete hT;
                  if (hB) delete hB;
                  TLatex t;
                  t.SetNDC(true);
                  t.SetTextFont(42);
                  t.SetTextSize(0.06);
                  t.DrawLatex(0.15, 0.55, "MISSING");
                  continue;
                }

                NormalizeToUnitArea(hR);
                NormalizeToUnitArea(hT);
                NormalizeToUnitArea(hB);

                hR->SetLineWidth(2);
                hR->SetMarkerStyle(20);
                hR->SetMarkerSize(0.92);
                hR->SetLineColor(kBlack);
                hR->SetMarkerColor(kBlack);

                hT->SetLineWidth(2);
                hT->SetMarkerStyle(24);
                hT->SetMarkerSize(0.92);
                hT->SetLineColor(kRed);
                hT->SetMarkerColor(kRed);

                hB->SetLineWidth(2);
                hB->SetMarkerStyle(21);
                hB->SetMarkerSize(0.92);
                hB->SetLineColor(kBlue);
                hB->SetMarkerColor(kBlue);

                const double ymax = std::max(std::max(hR->GetMaximum(), hT->GetMaximum()), hB->GetMaximum());
                const double yMaxPlot = ymax * 1.25;
                hR->SetMaximum(yMaxPlot);

                hR->SetTitle("");
                hR->GetXaxis()->SetTitle("x_{J#gamma}");
                hR->GetYaxis()->SetTitle("A.U.");
                hR->Draw("E1");
                hT->Draw("E1 same");
                hB->Draw("E1 same");

                const double xFloor = XFloorForPtBin(hReco->GetXaxis(), ib);
                DrawKinematicFloorLine(xFloor, 0.0, yMaxPlot);

                const string ptLab = AxisBinLabel(hReco->GetXaxis(), ib, "GeV", 0);

                TLatex tt;
                tt.SetNDC(true);
                tt.SetTextFont(42);
                tt.SetTextAlign(22);
                tt.SetTextSize(0.060);
                tt.DrawLatex(0.52, 0.95,
                  TString::Format("p_{T}^{#gamma} = %s  (R=%.1f)", ptLab.c_str(), R).Data()
                );

                // Legend (top-right)
                double lx1 = 0.50, ly1 = 0.70, lx2 = 0.86, ly2 = 0.90;
                double legTextSize = 0.048;
                if (ovTag.find("truthTaggedPhoJet") != std::string::npos)
                {
                  lx1 = 0.50; lx2 = 0.88;
                  legTextSize = 0.045;
                }

                TLegend leg(lx1, ly1, lx2, ly2);
                leg.SetTextFont(42);
                leg.SetTextSize(legTextSize);
                leg.SetFillStyle(0);
                leg.SetBorderSize(0);
                leg.AddEntry(hR, legReco.c_str(),      "ep");
                leg.AddEntry(hT, legTruth.c_str(),     "ep");
                leg.AddEntry(hB, legRecoTruth.c_str(), "ep");
                leg.DrawClone();

                keep.push_back(hR);
                keep.push_back(hT);
                keep.push_back(hB);
              }

              const string outName = "table2x3_overlay_shape.png";

              SaveCanvas(c, JoinPath(dirOvBase, outName));

              for (auto* h : keep) delete h;
            }

          // ==========================================================================
          // (C) 3x3 TABLE PAGES: ratios (requested) in SAME ovTag folder
          //    Output:
          //      table3x3_ratio.png  (or ..._pageN.png)
          //    Each pad shows BOTH:
          //      Reco(truth-cond)/Reco(all)  and  Reco(truth-cond)/Truth(reco-cond)
          //    plus unity line + kinematic floor line
          // ==========================================================================
          page = 0;

          for (int start = 1; start <= nPtOv; start += perPage)
          {
            ++page;

            TCanvas c(
              TString::Format("c_tbl_ratio_%s_%s_p%d", ovTag.c_str(), rKey.c_str(), page).Data(),
              "c_tbl_ratio", 1500, 1200
            );
            c.Divide(3,3, 0.001, 0.001);

            vector<TH1*> keep;
            keep.reserve(5 * perPage);

            for (int k = 0; k < perPage; ++k)
            {
              const int ib = start + k;
              c.cd(k+1);

              gPad->SetLeftMargin(0.14);
              gPad->SetRightMargin(0.05);
              gPad->SetBottomMargin(0.14);
              gPad->SetTopMargin(0.10);

              if (ib > nPtOv)
              {
                TLatex t;
                t.SetNDC(true);
                t.SetTextFont(42);
                t.SetTextSize(0.06);
                t.DrawLatex(0.20, 0.55, "EMPTY");
                continue;
              }

              TH1* hR = ProjectY_AtXbin_TH3(hReco, ib,
                TString::Format("tblR_%s_reco_%s_%d_p%d", ovTag.c_str(), rKey.c_str(), ib, page).Data()
              );
              TH1* hT = ProjectY_AtXbin_TH3(hTruth, ib,
                TString::Format("tblR_%s_truth_%s_%d_p%d", ovTag.c_str(), rKey.c_str(), ib, page).Data()
              );
              TH1* hB = ProjectY_AtXbin_TH3(hRecoTruth, ib,
                TString::Format("tblR_%s_recotruth_%s_%d_p%d", ovTag.c_str(), rKey.c_str(), ib, page).Data()
              );

              if (hR) { hR->SetDirectory(nullptr); EnsureSumw2(hR); }
              if (hT) { hT->SetDirectory(nullptr); EnsureSumw2(hT); }
              if (hB) { hB->SetDirectory(nullptr); EnsureSumw2(hB); }

              if (!hR || !hT || !hB || (hR->GetEntries() <= 0.0 && hT->GetEntries() <= 0.0 && hB->GetEntries() <= 0.0))
              {
                if (hR) delete hR;
                if (hT) delete hT;
                if (hB) delete hB;
                TLatex t;
                t.SetNDC(true);
                t.SetTextFont(42);
                t.SetTextSize(0.06);
                t.DrawLatex(0.15, 0.55, "MISSING");
                continue;
              }

              NormalizeToUnitArea(hR);
              NormalizeToUnitArea(hT);
              NormalizeToUnitArea(hB);

              TH1* rB_over_R = MakeSafeRatio(hB, hR,
                TString::Format("tbl_ratio_%s_BoverR_%s_%d_p%d", ovTag.c_str(), rKey.c_str(), ib, page).Data()
              );
              TH1* rB_over_T = MakeSafeRatio(hB, hT,
                TString::Format("tbl_ratio_%s_BoverT_%s_%d_p%d", ovTag.c_str(), rKey.c_str(), ib, page).Data()
              );

              if (!rB_over_R || !rB_over_T)
              {
                if (hR) delete hR;
                if (hT) delete hT;
                if (hB) delete hB;
                if (rB_over_R) delete rB_over_R;
                if (rB_over_T) delete rB_over_T;

                TLatex t;
                t.SetNDC(true);
                t.SetTextFont(42);
                t.SetTextSize(0.06);
                t.DrawLatex(0.15, 0.55, "MISSING");
                continue;
              }

              rB_over_R->SetTitle("");
              rB_over_R->GetXaxis()->SetTitle("x_{J#gamma}");
              rB_over_R->GetYaxis()->SetTitle("Ratio");
              rB_over_R->SetMinimum(ratioYMin);
              rB_over_R->SetMaximum(ratioYMax);

              rB_over_R->SetLineWidth(2);
              rB_over_R->SetMarkerStyle(20);
              rB_over_R->SetMarkerSize(0.92);
              rB_over_R->SetLineColor(kBlue);
              rB_over_R->SetMarkerColor(kBlue);

              rB_over_T->SetLineWidth(2);
              rB_over_T->SetMarkerStyle(24);
              rB_over_T->SetMarkerSize(0.92);
              rB_over_T->SetLineColor(kRed);
              rB_over_T->SetMarkerColor(kRed);

              rB_over_R->Draw("E1");
              rB_over_T->Draw("E1 same");

              // Unity line
              const double xMin = rB_over_R->GetXaxis()->GetXmin();
              const double xMax = rB_over_R->GetXaxis()->GetXmax();
              TLine unity(xMin, 1.0, xMax, 1.0);
              unity.SetLineColor(kGray + 2);
              unity.SetLineStyle(2);
              unity.SetLineWidth(2);
              unity.Draw("same");

              // Kinematic floor line
              const double xFloor = XFloorForPtBin(hReco->GetXaxis(), ib);
              DrawKinematicFloorLine(xFloor, ratioYMin, ratioYMax);

              const string ptLab = AxisBinLabel(hReco->GetXaxis(), ib, "GeV", 0);

              TLatex tt;
              tt.SetNDC(true);
              tt.SetTextFont(42);
              tt.SetTextAlign(22);
              tt.SetTextSize(0.060);
              tt.DrawLatex(0.52, 0.95,
                TString::Format("p_{T}^{#gamma} = %s  (R=%.1f)", ptLab.c_str(), R).Data()
              );

              // Legend (bottom-left to keep the title area clear)
              TLegend leg(0.16, 0.18, 0.92, 0.34);
              leg.SetTextFont(42);
              leg.SetTextSize(0.048);
              leg.SetFillStyle(0);
              leg.SetBorderSize(0);
              leg.AddEntry(rB_over_R, "Reco(truth-cond) / Reco(all)", "ep");
              leg.AddEntry(rB_over_T, "Reco(truth-cond) / Truth(reco-cond)", "ep");
              leg.DrawClone();

              // Keep + cleanup
              keep.push_back(hR);
              keep.push_back(hT);
              keep.push_back(hB);
              keep.push_back(rB_over_R);
              keep.push_back(rB_over_T);
            }

            const string outName =
              (nPtOv <= perPage)
                ? "table3x3_ratio.png"
                : TString::Format("table3x3_ratio_page%d.png", page).Data();

              SaveCanvas(c, JoinPath(dirOvBase, outName));
              cout << ANSI_BOLD_GRN
                   << "  [WROTE ratio/table] " << JoinPath(dirOvBase, outName)
                   << ANSI_RESET << "\n";

              for (auto* h : keep) delete h;
          }
        };


        MakeOverlayShape_TH3xJ_3way(
          H.hReco_xJ,
          H.hTrut_xJ,
          H.hRecoTruthTagged_xJ,
          "RECO_vs_TRUTHrecoConditioned_vs_RECO_truthTaggedPhoJet",
          "Reco (all)",
          "Truth (reco #gamma+jet cond)",
          "Reco (truth #gamma+jet cond)",
          {"Overlay (shape): RECO vs TRUTH vs RECO truth-conditioned"}
        );

        // RECO alpha-cut tables into subfolders (requested)
        // + in the *parent* alphaCuts folder, write a 3x3 table overlaying ALL alpha cuts per pT bin.
        if (H.hReco_xJ)
        {
          const vector<double> alphaMaxCuts = {0.20, 0.30, 0.40, 0.50};

          auto AlphaTag = [&](double aMax)->string
          {
            std::ostringstream s;
            s << std::fixed << std::setprecision(2) << aMax;
            string t = s.str();
            std::replace(t.begin(), t.end(), '.', 'p');
            return string("alphaLT") + t;
          };

          const string dirAlphaBase = JoinPath(D.dirXJProjReco, "alphaCuts");
          EnsureDir(dirAlphaBase);

          // (A) Existing behavior: per-alpha-cut folders with their own 3x3 tables
          for (double aMax : alphaMaxCuts)
          {
            const string aDir = JoinPath(dirAlphaBase, AlphaTag(aMax));
            EnsureDir(aDir);
            Make3x3Table_xJ_FromTH3_AlphaCut(H.hReco_xJ, aDir, "RECO", aMax);
          }

          // (B) NEW: One 3x3 table in the PARENT alphaCuts folder overlaying all alpha cuts
          // Marker requirement: circles, two open + two filled; all different colors.
          auto Make3x3Table_xJ_RECO_AlphaCutsOverlay =
            [&](const TH3* h3, const string& outBaseDir)
          {
            if (!h3) return;

            struct CutStyle
            {
              double aMax;
              int    color;
              int    marker;
            };

            // Two open circles (24) + two filled circles (20), all different colors
            const vector<CutStyle> styles =
            {
              {0.20, 2, 24}, // open circle, red
              {0.30, 4, 20}, // filled circle, blue
              {0.40, 6, 24}, // open circle, magenta
              {0.50, 1, 20}  // filled circle, black
            };

            const int n = h3->GetXaxis()->GetNbins();
            const int perPage = 9;
            int page = 0;

            for (int start = 1; start <= n; start += perPage)
            {
              ++page;

              TCanvas c(
                TString::Format("c_tbl_xJ_%s_%s_RECO_alphaCutsOverlay_p%d",
                  ds.label.c_str(), rKey.c_str(), page).Data(),
                "c_tbl_xJ_alphaCutsOverlay", 1500, 1200
              );

              c.Divide(3,3, 0.001, 0.001);

              std::vector<TH1*> keep;
              keep.reserve(perPage * styles.size());

              for (int k = 0; k < perPage; ++k)
              {
                const int ib = start + k;
                c.cd(k+1);

                gPad->SetLeftMargin(0.14);
                gPad->SetRightMargin(0.05);
                gPad->SetBottomMargin(0.14);
                gPad->SetTopMargin(0.10);
                gPad->SetLogy(false);

                if (ib > n)
                {
                  TLatex t;
                  t.SetNDC(true);
                  t.SetTextFont(42);
                  t.SetTextSize(0.06);
                  t.DrawLatex(0.20, 0.55, "EMPTY");
                  continue;
                }

                std::vector<TH1*> hs;
                hs.reserve(styles.size());

                double ymax = 0.0;

                for (const auto& st : styles)
                {
                  const string aTag = AlphaTag(st.aMax);

                  TH1* hx = ProjectY_AtXbin_AndAlphaMax_TH3(
                    h3, ib, st.aMax,
                    TString::Format("jes3_xJ_tbl_%s_RECO_alphaOv_%s_b%d",
                      rKey.c_str(), aTag.c_str(), ib).Data()
                  );

                  if (!hx)
                  {
                    hs.push_back(nullptr);
                    continue;
                  }

                  hx->SetDirectory(nullptr);
                  EnsureSumw2(hx);

                  hx->SetTitle("");
                  hx->SetLineWidth(2);
                  hx->SetLineColor(st.color);
                  hx->SetMarkerStyle(st.marker);
                  hx->SetMarkerSize(0.95);
                  hx->SetMarkerColor(st.color);
                  hx->SetFillStyle(0);

                  hx->GetXaxis()->SetTitle("x_{J#gamma}");
                  hx->GetYaxis()->SetTitle((ds.isSim && IsWeightedSIMSelected()) ? "Counts / pb^{-1}" : "Counts");

                  ymax = std::max(ymax, hx->GetMaximum());

                  hs.push_back(hx);
                  keep.push_back(hx);
                }

                TH1* first = nullptr;
                for (auto* h : hs) { if (h) { first = h; break; } }

                if (!first)
                {
                  TLatex t;
                  t.SetNDC(true);
                  t.SetTextFont(42);
                  t.SetTextSize(0.06);
                  t.DrawLatex(0.15, 0.55, "MISSING");
                  continue;
                }

                first->SetMaximum((ymax > 0.0) ? (ymax * 1.25) : 1.0);

                bool drawn = false;
                for (auto* h : hs)
                {
                  if (!h) continue;
                  if (!drawn) { h->Draw("E1"); drawn = true; }
                  else        { h->Draw("E1 same"); }
                }

                const string ptLab = AxisBinLabel(h3->GetXaxis(), ib, "GeV", 0);

                TLatex tt;
                tt.SetNDC(true);
                tt.SetTextFont(42);
                tt.SetTextAlign(13);
                tt.SetTextSize(0.055);
                tt.DrawLatex(0.16, 0.88,
                  TString::Format("p_{T}^{#gamma}: %s", ptLab.c_str()).Data()
                );

                TLegend leg(0.48, 0.56, 0.93, 0.92);
                leg.SetTextFont(42);
                leg.SetTextSize(0.048);
                leg.SetFillStyle(0);
                leg.SetBorderSize(0);

                for (size_t j = 0; j < styles.size(); ++j)
                {
                  if (!hs[j]) continue;
                  leg.AddEntry(hs[j],
                    TString::Format("#alpha < %.2f", styles[j].aMax).Data(),
                    "ep"
                  );
                }
                leg.DrawClone();
              }

              string outName;
              if (n <= perPage)
              {
                outName = "table3x3_xJ_RECO_alphaCutsOverlay.png";
              }
              else
              {
                outName = TString::Format("table3x3_xJ_RECO_alphaCutsOverlay_page%d.png", page).Data();
              }

              SaveCanvas(c, JoinPath(outBaseDir, outName));

              for (auto* h : keep) delete h;
            }
          };

            Make3x3Table_xJ_RECO_AlphaCutsOverlay(H.hReco_xJ, dirAlphaBase);

            // -----------------------------------------------------------------------------
            // (B2) NEW: same 3x3 alpha overlay, but only THREE curves:
            //      "no alpha cut", "alpha < 0.20", "alpha < 0.50"
            //      Output goes to the SAME directory as the 4-curve overlay:
            //        <...>/xJ_fromJES3/RECO/alphaCuts/
            // -----------------------------------------------------------------------------
            auto Make3x3Table_xJ_RECO_AlphaCutsOverlay_3Curves =
              [&](const TH3* h3, const string& outBaseDir)
            {
              if (!h3) return;

              struct CutStyle
              {
                double aMax;
                int    color;
                int    marker;
                string label;
                string tag;
              };

              // 3 curves: no cut (integrate full alpha axis), alpha<0.20, alpha<0.50
              const vector<CutStyle> styles =
              {
                {999.0, 1, 20, "no #alpha cut", "noCut"},   // full alpha integration
                {0.20,  2, 24, "#alpha < 0.20", "a020"},    // red open
                {0.50,  4, 20, "#alpha < 0.50", "a050"}     // blue filled
              };

              const int n = h3->GetXaxis()->GetNbins();
              const int perPage = 9;
              int page = 0;

              for (int start = 1; start <= n; start += perPage)
              {
                ++page;

                TCanvas c(
                  TString::Format("c_tbl_xJ_%s_%s_RECO_alphaCuts3_p%d",
                    ds.label.c_str(), rKey.c_str(), page).Data(),
                  "c_tbl_xJ_alphaCuts3", 1500, 1200
                );

                c.Divide(3,3, 0.001, 0.001);

                std::vector<TH1*> keep;
                keep.reserve(perPage * styles.size());

                for (int k = 0; k < perPage; ++k)
                {
                  const int ib = start + k;
                  c.cd(k+1);

                  gPad->SetLeftMargin(0.14);
                  gPad->SetRightMargin(0.05);
                  gPad->SetBottomMargin(0.14);
                  gPad->SetTopMargin(0.10);
                  gPad->SetLogy(false);

                  if (ib > n)
                  {
                    TLatex t;
                    t.SetNDC(true);
                    t.SetTextFont(42);
                    t.SetTextSize(0.06);
                    t.DrawLatex(0.20, 0.55, "EMPTY");
                    continue;
                  }

                  std::vector<TH1*> hs;
                  hs.reserve(styles.size());

                  double ymax = 0.0;

                  for (const auto& st : styles)
                  {
                    TH1* hx = ProjectY_AtXbin_AndAlphaMax_TH3(
                      h3, ib, st.aMax,
                      TString::Format("jes3_xJ_tbl_%s_RECO_alphaOv3_%s_b%d",
                        rKey.c_str(), st.tag.c_str(), ib).Data()
                    );

                    if (!hx)
                    {
                      hs.push_back(nullptr);
                      continue;
                    }

                    hx->SetDirectory(nullptr);
                    EnsureSumw2(hx);

                    hx->SetTitle("");
                    hx->SetLineWidth(2);
                    hx->SetLineColor(st.color);
                    hx->SetMarkerStyle(st.marker);
                    hx->SetMarkerSize(0.95);
                    hx->SetMarkerColor(st.color);
                    hx->SetFillStyle(0);

                    hx->GetXaxis()->SetTitle("x_{J#gamma}");
                    hx->GetYaxis()->SetTitle((ds.isSim && IsWeightedSIMSelected()) ? "Counts / pb^{-1}" : "Counts");

                    ymax = std::max(ymax, hx->GetMaximum());

                    hs.push_back(hx);
                    keep.push_back(hx);
                  }

                  TH1* first = nullptr;
                  for (auto* h : hs) { if (h) { first = h; break; } }

                  if (!first)
                  {
                    TLatex t;
                    t.SetNDC(true);
                    t.SetTextFont(42);
                    t.SetTextSize(0.06);
                    t.DrawLatex(0.15, 0.55, "MISSING");
                    continue;
                  }

                  first->SetMaximum((ymax > 0.0) ? (ymax * 1.25) : 1.0);

                  bool drawn = false;
                  for (auto* h : hs)
                  {
                    if (!h) continue;
                    if (!drawn) { h->Draw("E1"); drawn = true; }
                    else        { h->Draw("E1 same"); }
                  }

                  const string ptLab = AxisBinLabel(h3->GetXaxis(), ib, "GeV", 0);

                  TLatex tt;
                  tt.SetNDC(true);
                  tt.SetTextFont(42);
                  tt.SetTextAlign(13);
                  tt.SetTextSize(0.055);
                  tt.DrawLatex(0.16, 0.88,
                    TString::Format("p_{T}^{#gamma}: %s", ptLab.c_str()).Data()
                  );

                  TLegend leg(0.48, 0.62, 0.93, 0.92);
                  leg.SetTextFont(42);
                  leg.SetTextSize(0.048);
                  leg.SetFillStyle(0);
                  leg.SetBorderSize(0);

                  for (size_t j = 0; j < styles.size(); ++j)
                  {
                    if (!hs[j]) continue;
                    leg.AddEntry(hs[j], styles[j].label.c_str(), "ep");
                  }
                  leg.DrawClone();
                }

                string outName;
                if (n <= perPage)
                {
                  outName = "table3x3_xJ_RECO_alphaCutsOverlay_3curves.png";
                }
                else
                {
                  outName = TString::Format("table3x3_xJ_RECO_alphaCutsOverlay_3curves_page%d.png", page).Data();
                }

                SaveCanvas(c, JoinPath(outBaseDir, outName));

                for (auto* h : keep) delete h;
              }
            };

            Make3x3Table_xJ_RECO_AlphaCutsOverlay_3Curves(H.hReco_xJ, dirAlphaBase);

            // -----------------------------------------------------------------------------
            // (C) NEW: 3x3 table per pT bin of "How much does each alpha cut remove?"
            // For each pT^gamma bin, compute:
            //   f_removed(aMax) = N(alpha > aMax) / N(total)
            // using the RECO JES3 TH3 (integrated over xJ).
            // Output goes to the SAME dir as the xJ alpha-cuts overlay:
            //   <...>/xJ_fromJES3/RECO/alphaCuts/
            // -----------------------------------------------------------------------------
            auto Make3x3Table_AlphaCutRemovedFractions =
              [&](const TH3* h3, const string& outBaseDir)
            {
              if (!h3) return;

              const vector<double> aCuts = alphaMaxCuts;

              const int n = h3->GetXaxis()->GetNbins();
              const int perPage = 9;
              int page = 0;

              for (int start = 1; start <= n; start += perPage)
              {
                ++page;

                TCanvas c(
                  TString::Format("c_tbl_alphaRemoved_%s_%s_p%d",
                    ds.label.c_str(), rKey.c_str(), page).Data(),
                  "c_tbl_alphaRemoved", 1500, 1200
                );
                c.Divide(3,2,0.001,0.001);

                std::vector<TObject*> keep;
                keep.reserve(perPage * 3);

                for (int k = 0; k < perPage; ++k)
                {
                  const int ib = start + k;
                  c.cd(k+1);

                  gPad->SetLeftMargin(0.14);
                  gPad->SetRightMargin(0.05);
                  gPad->SetBottomMargin(0.14);
                  gPad->SetTopMargin(0.10);
                  gPad->SetLogy(false);
                  gPad->SetTicks(1,1);

                  if (ib > n)
                  {
                    TLatex t;
                    t.SetNDC(true);
                    t.SetTextFont(42);
                    t.SetTextSize(0.06);
                    t.DrawLatex(0.20, 0.55, "EMPTY");
                    continue;
                  }

                  TH1* hA = ProjectZ_AtXbin_TH3(
                    h3, ib,
                    TString::Format("jes3_alphaProj_%s_%d", rKey.c_str(), ib).Data()
                  );

                  if (!hA)
                  {
                    TLatex t;
                    t.SetNDC(true);
                    t.SetTextFont(42);
                    t.SetTextSize(0.06);
                    t.DrawLatex(0.15, 0.55, "MISSING");
                    continue;
                  }

                  hA->SetDirectory(nullptr);
                  EnsureSumw2(hA);

                  const int nAlphaBins = hA->GetNbinsX();
                  const double nTot = hA->Integral(1, nAlphaBins);

                  std::vector<double> xs;
                  std::vector<double> ys;
                  xs.reserve(aCuts.size());
                  ys.reserve(aCuts.size());

                  for (double aMax : aCuts)
                  {
                    double nAbove = 0.0;
                    for (int iz = 1; iz <= nAlphaBins; ++iz)
                    {
                      const double aC = hA->GetXaxis()->GetBinCenter(iz);
                      if (aC > aMax) nAbove += hA->GetBinContent(iz);
                    }
                    const double fRemoved = (nTot > 0.0) ? (nAbove / nTot) : 0.0;
                    xs.push_back(aMax);
                    ys.push_back(fRemoved);
                  }

                  // Frame for clear axes in small pads
                  TH1F* frame = new TH1F(
                    TString::Format("h_frame_alphaRemoved_%s_%s_%d", ds.label.c_str(), rKey.c_str(), ib).Data(),
                    "",
                    1, 0.15, 0.55
                  );
                  frame->SetDirectory(nullptr);
                  frame->SetStats(0);
                  frame->SetMinimum(0.0);
                  frame->SetMaximum(1.0);
                  frame->GetXaxis()->SetTitle("#alpha_{max}");
                  frame->GetYaxis()->SetTitle("f(#alpha > #alpha_{max})");
                  frame->GetXaxis()->SetTitleSize(0.060);
                  frame->GetYaxis()->SetTitleSize(0.060);
                  frame->GetXaxis()->SetLabelSize(0.050);
                  frame->GetYaxis()->SetLabelSize(0.050);
                  frame->GetYaxis()->SetTitleOffset(1.10);
                  frame->Draw("hist");

                  TGraph* g = new TGraph((int)xs.size());
                  for (int j = 0; j < (int)xs.size(); ++j) g->SetPoint(j, xs[j], ys[j]);
                  g->SetLineWidth(2);
                  g->SetLineColor(kBlack);
                  g->SetMarkerStyle(20);
                  g->SetMarkerSize(1.0);
                  g->SetMarkerColor(kBlack);
                  g->Draw("PL same");

                  // Print values near points
                  TLatex tt;
                  tt.SetTextFont(42);
                  tt.SetTextAlign(12);
                  tt.SetTextSize(0.050);
                  for (int j = 0; j < (int)xs.size(); ++j)
                  {
                    const double xP = xs[j];
                    const double yP = ys[j];
                    tt.DrawLatex(xP + 0.005, std::min(yP + 0.04, 0.95), TString::Format("%.3f", yP).Data());
                  }

                  const string ptLab = AxisBinLabel(h3->GetXaxis(), ib, "GeV", 0);
                  TLatex ptl;
                  ptl.SetNDC(true);
                  ptl.SetTextFont(42);
                  ptl.SetTextAlign(13);
                  ptl.SetTextSize(0.055);
                  ptl.DrawLatex(0.16, 0.88, TString::Format("p_{T}^{#gamma}: %s", ptLab.c_str()).Data());

                  keep.push_back(hA);
                  keep.push_back(frame);
                  keep.push_back(g);
                }

                string outName;
                if (n <= perPage)
                  outName = "table3x3_alphaCutRemovedFractions_RECO.png";
                else
                  outName = TString::Format("table3x3_alphaCutRemovedFractions_RECO_page%d.png", page).Data();

                SaveCanvas(c, JoinPath(outBaseDir, outName));

                for (auto* o : keep) delete o;
              }
            };

            Make3x3Table_AlphaCutRemovedFractions(H.hReco_xJ, dirAlphaBase);

            // -----------------------------------------------------------------------------
            // (D) NEW: 3x3 table of h_jet2Pt_<rKey><slice> with log-y + P(no jet2)
            // Since online fills jet2Pt=0 when jet2 doesn't exist, bin1 (~[0,0.5)) measures no-jet2.
            // Output goes into the SAME alphaCuts folder as the xJ overlay.
            // -----------------------------------------------------------------------------
            auto Make3x3Table_Jet2Pt_LogZeroFrac =
              [&](const string& outBaseDir)
            {
              TCanvas c(
                TString::Format("c_tbl_jet2Pt_zeroFrac_%s_%s", ds.label.c_str(), rKey.c_str()).Data(),
                "c_tbl_jet2Pt_zeroFrac", 1500, 1200
              );
              c.Divide(3,2,0.001,0.001);

              std::vector<TH1*> keep;
              keep.reserve(kNPtBins);

              for (int i = 0; i < kNPtBins; ++i)
              {
                const PtBin& b = PtBins()[i];
                c.cd(i+1);

                gPad->SetLeftMargin(0.14);
                gPad->SetRightMargin(0.05);
                gPad->SetBottomMargin(0.14);
                gPad->SetTopMargin(0.10);
                gPad->SetLogy(true);
                gPad->SetTicks(1,1);

                TH1* h = GetObj<TH1>(ds, string("h_jet2Pt_") + rKey + b.suffix, false, false, false);
                if (!h)
                {
                  TLatex t;
                  t.SetNDC(true);
                  t.SetTextFont(42);
                  t.SetTextSize(0.06);
                  t.DrawLatex(0.15, 0.55, "MISSING h_jet2Pt");
                  continue;
                }

                TH1* hc = CloneTH1(h, TString::Format("h_jet2Pt_clone_%s_%s_%d", ds.label.c_str(), rKey.c_str(), i).Data());
                if (!hc) continue;

                hc->SetDirectory(nullptr);
                EnsureSumw2(hc);

                hc->SetTitle("");
                hc->SetLineWidth(2);
                hc->SetLineColor(kBlack);
                hc->SetMarkerStyle(20);
                hc->SetMarkerSize(0.85);
                hc->SetMarkerColor(kBlack);

                hc->GetXaxis()->SetTitle("p_{T}^{jet2} [GeV]");
                hc->GetYaxis()->SetTitle((ds.isSim && IsWeightedSIMSelected()) ? "Counts / pb^{-1}" : "Counts");

                // log-y safety
                hc->SetMinimum(0.5);

                hc->Draw("E1");

                const double nTot = hc->Integral(1, hc->GetNbinsX());
                const double n0   = hc->GetBinContent(1); // [0,0.5): your "no jet2 -> 0" fills
                const double f0   = (nTot > 0.0) ? (n0 / nTot) : 0.0;

                TLatex t0;
                t0.SetNDC(true);
                t0.SetTextFont(42);
                t0.SetTextAlign(13);
                t0.SetTextSize(0.055);
                t0.DrawLatex(0.16, 0.80, TString::Format("P(no jet2) #approx %.3f", f0).Data());

                TLatex ptl;
                ptl.SetNDC(true);
                ptl.SetTextFont(42);
                ptl.SetTextAlign(13);
                ptl.SetTextSize(0.055);
                ptl.DrawLatex(0.16, 0.88, TString::Format("p_{T}^{#gamma}: %s", b.label.c_str()).Data());

                keep.push_back(hc);
              }

              SaveCanvas(c, JoinPath(outBaseDir, "table3x3_jet2Pt_logy_zeroFrac.png"));

              for (auto* h : keep) delete h;
            };

            Make3x3Table_Jet2Pt_LogZeroFrac(dirAlphaBase);

            // -----------------------------------------------------------------------------
            // (E) NEW: mean <N_recoil jets> vs pTgamma (sanity check: if ~1, alpha often = 0)
            // Output goes into the SAME alphaCuts folder as the xJ overlay.
            // -----------------------------------------------------------------------------
            {
              if (TProfile* pN = GetObj<TProfile>(ds, "p_nRecoilJets_vs_pTgamma_" + rKey, false, false, false))
              {
                TProfile* pc = (TProfile*)pN->Clone(
                  TString::Format("p_nRecoilJets_clone_%s_%s", ds.label.c_str(), rKey.c_str()).Data()
                );
                if (pc)
                {
                  pc->SetDirectory(nullptr);
                  pc->SetTitle("");
                  pc->SetLineWidth(2);
                  pc->SetMarkerStyle(20);
                  pc->SetMarkerSize(0.95);

                  TCanvas cN(
                    TString::Format("c_meanNrecoil_%s_%s", ds.label.c_str(), rKey.c_str()).Data(),
                    "c_meanNrecoil", 900, 700
                  );
                  ApplyCanvasMargins1D(cN);

                  pc->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                  pc->GetYaxis()->SetTitle("<N_{recoil jets}>");
                  pc->SetMinimum(0.0);
                  pc->Draw("E1");

                  DrawLatexLines(0.14,0.92, DefaultHeaderLines(ds), 0.034, 0.045);
                  DrawLatexLines(0.14,0.82,
                    { "Mean number of recoil jets vs p_{T}^{#gamma}",
                      rKey + TString::Format(" (R=%.1f)", R).Data() },
                    0.030, 0.040
                  );

                  SaveCanvas(cN, JoinPath(dirAlphaBase, "profile_meanNrecoilJets_vs_pTgamma.png"));
                  delete pc;
                }
              }
            }

            // -----------------------------------------------------------------------------
            // (F) NEW: 3x3 table of <xJ>(alpha) profiles (xJ-alpha correlation diagnostic)
            // Uses TH3 h_JES3_pT_xJ_alpha_<rKey>:
            //   - Project (xJ,alpha) at fixed pT bin
            //   - ProfileY => mean(xJ) vs alpha
            // Output goes into the SAME alphaCuts folder as the xJ overlay.
            // -----------------------------------------------------------------------------
            auto Make3x3Table_ProfileMeanxJ_vs_Alpha =
              [&](const TH3* h3, const string& outBaseDir)
            {
              if (!h3) return;

              const int n = h3->GetXaxis()->GetNbins();
              const int perPage = 9;
              int page = 0;

              for (int start = 1; start <= n; start += perPage)
              {
                ++page;

                TCanvas c(
                  TString::Format("c_tbl_meanxJ_vs_alpha_%s_%s_p%d",
                    ds.label.c_str(), rKey.c_str(), page).Data(),
                  "c_tbl_meanxJ_vs_alpha", 1500, 1200
                );
                c.Divide(3,2,0.001,0.001);

                std::vector<TObject*> keep;
                keep.reserve(perPage * 2);

                for (int k = 0; k < perPage; ++k)
                {
                  const int ib = start + k;
                  c.cd(k+1);

                  gPad->SetLeftMargin(0.14);
                  gPad->SetRightMargin(0.05);
                  gPad->SetBottomMargin(0.14);
                  gPad->SetTopMargin(0.10);
                  gPad->SetLogy(false);
                  gPad->SetTicks(1,1);

                  if (ib > n)
                  {
                    TLatex t;
                    t.SetNDC(true);
                    t.SetTextFont(42);
                    t.SetTextSize(0.06);
                    t.DrawLatex(0.20, 0.55, "EMPTY");
                    continue;
                  }

                  // Project YZ at this pT bin: Y=xJ, Z=alpha -> TH2 with X=xJ, Y=alpha
                  TH2* h2 = ProjectYZ_AtXbin_TH3(
                    h3, ib,
                    TString::Format("jes3_xJ_vs_alpha_%s_%d", rKey.c_str(), ib).Data()
                  );
                  if (!h2)
                  {
                    TLatex t;
                    t.SetNDC(true);
                    t.SetTextFont(42);
                    t.SetTextSize(0.06);
                    t.DrawLatex(0.15, 0.55, "MISSING");
                    continue;
                  }
                  h2->SetDirectory(nullptr);
                  EnsureSumw2(h2);

                  // ProfileY: mean X (xJ) vs Y (alpha)
                  TProfile* p = h2->ProfileY(
                    TString::Format("p_meanxJ_vs_alpha_%s_%d", rKey.c_str(), ib).Data()
                  );
                  if (!p)
                  {
                    delete h2;
                    continue;
                  }
                  p->SetDirectory(nullptr);

                  p->SetTitle("");
                  p->SetLineWidth(2);
                  p->SetLineColor(kBlack);
                  p->SetMarkerStyle(20);
                  p->SetMarkerSize(0.85);
                  p->SetMarkerColor(kBlack);

                  p->GetXaxis()->SetTitle("#alpha");
                  p->GetYaxis()->SetTitle("<x_{J#gamma}>");
                  p->SetMinimum(0.0);
                  p->SetMaximum(1.2);

                  p->Draw("E1");

                  const string ptLab = AxisBinLabel(h3->GetXaxis(), ib, "GeV", 0);
                  TLatex ptl;
                  ptl.SetNDC(true);
                  ptl.SetTextFont(42);
                  ptl.SetTextAlign(13);
                  ptl.SetTextSize(0.055);
                  ptl.DrawLatex(0.16, 0.88, TString::Format("p_{T}^{#gamma}: %s", ptLab.c_str()).Data());

                  keep.push_back(h2);
                  keep.push_back(p);
                }

                string outName;
                if (n <= perPage)
                  outName = "table3x3_meanxJ_vs_alphaProfile_RECO.png";
                else
                  outName = TString::Format("table3x3_meanxJ_vs_alphaProfile_RECO_page%d.png", page).Data();

                SaveCanvas(c, JoinPath(outBaseDir, outName));

                for (auto* o : keep) delete o;
              }
            };

            Make3x3Table_ProfileMeanxJ_vs_Alpha(H.hReco_xJ, dirAlphaBase);
        }

        // =============================================================================
        // SIM-only: JES3 efficiency dashboard (leading-jet focused)
        //
        // Outputs live in:
        //   <...>/<rKey>/xJ_fromJES3/Efficiency/...
        //
        // This block:
        //   - prints verbose terminal summaries (including explicit output paths)
        //   - writes a summary text file next to the plots
        //   - only runs in SIM (truth-tagged + truth response hists exist)
        // =============================================================================
        if (ds.isSim)
        {
          // Make sure output dirs exist (MakeJes3Dirs also creates these, but keep this robust)
          EnsureDir(D.dirXJProjEff);
          EnsureDir(D.dirXJProjEffLeadMatch);
          EnsureDir(D.dirXJProjEffLeadJetResponse);
          EnsureDir(D.dirXJProjEffTagFractions3x3);

          cout << ANSI_BOLD_CYN
               << "\n[JES3 EfficiencyDashboard] " << ds.label << "  (" << rKey << ", R=" << R << ")"
               << ANSI_RESET << "\n";
          cout << "  Overlay root   : " << D.dirXJProjOverlay << "\n";
          cout << "  Efficiency root: " << D.dirXJProjEff << "\n";
          cout << "    - LeadTruthRecoilMatch : " << D.dirXJProjEffLeadMatch << "\n";
          cout << "    - LeadJetResponse      : " << D.dirXJProjEffLeadJetResponse << "\n";
          cout << "    - TagFractions_3x3     : " << D.dirXJProjEffTagFractions3x3 << "\n";

          vector<string> effSummary;
          effSummary.push_back("JES3 Efficiency Dashboard (SIM only)");
          effSummary.push_back("Dataset: " + ds.label);
          effSummary.push_back("rKey: " + rKey + TString::Format("  R=%.1f", R).Data());
          effSummary.push_back("Overlay root: " + D.dirXJProjOverlay);
          effSummary.push_back("Efficiency root: " + D.dirXJProjEff);
          effSummary.push_back("");

          // ---------------------------------------------------------------------------
          // (1) Lead truth recoil match efficiency vs truth photon pT
          //
          // Inputs (filled online in RecoilJets.cc):
          //   h_leadTruthRecoilMatch_den_pTgammaTruth_<rKey>
          //   h_leadTruthRecoilMatch_num_pTgammaTruth_<rKey>
          //   h_leadTruthRecoilMatch_missA_pTgammaTruth_<rKey>
          //   h_leadTruthRecoilMatch_missB_pTgammaTruth_<rKey>
          // ---------------------------------------------------------------------------
          {
              const string hDenName    = "h_leadTruthRecoilMatch_den_pTgammaTruth_" + rKey;
              const string hNumName    = "h_leadTruthRecoilMatch_num_pTgammaTruth_" + rKey;
              const string hMissAName  = "h_leadTruthRecoilMatch_missA_pTgammaTruth_" + rKey;
              const string hMissA1Name = "h_leadTruthRecoilMatch_missA1_pTgammaTruth_" + rKey;
              const string hMissA2Name = "h_leadTruthRecoilMatch_missA2_pTgammaTruth_" + rKey;
              const string hMissBName  = "h_leadTruthRecoilMatch_missB_pTgammaTruth_" + rKey;

              TH1* hDen    = GetObj<TH1>(ds, hDenName,    false, false, false);
              TH1* hNum    = GetObj<TH1>(ds, hNumName,    false, false, false);
              TH1* hMissA  = GetObj<TH1>(ds, hMissAName,  false, false, false);
              TH1* hMissA1 = GetObj<TH1>(ds, hMissA1Name, false, false, false);
              TH1* hMissA2 = GetObj<TH1>(ds, hMissA2Name, false, false, false);
              TH1* hMissB  = GetObj<TH1>(ds, hMissBName,  false, false, false);


            if (!hDen || !hNum || !hMissA || !hMissB)
            {
              cout << ANSI_BOLD_YEL << "  [LeadTruthRecoilMatch] Missing one or more inputs; skipping.\n" << ANSI_RESET;
              if (!hDen)   cout << "    - MISSING: " << hDenName   << "\n";
              if (!hNum)   cout << "    - MISSING: " << hNumName   << "\n";
              if (!hMissA) cout << "    - MISSING: " << hMissAName << "\n";
              if (!hMissB) cout << "    - MISSING: " << hMissBName << "\n";

              effSummary.push_back("LeadTruthRecoilMatch: MISSING INPUTS (skipped)");
              effSummary.push_back("");
            }
            else
            {
                // Build the single efficiency curve:
                //   eps_lead(pTgammaTruth) = NUM / DEN
                //
                // Definition:
                //   DEN: truth-leading recoil jet exists (fiducial truth definition)
                //   NUM: that truth-leading recoil jet is matched to the chosen reco recoil jet
                //
                // Output (exact):
                //   <...>/<rKey>/xJ_fromJES3/Efficiency/LeadTruthRecoilMatch/
                //     leadRecoilJetMatchEff_vs_pTgammaTruth_<rKey>.png

                TH1* denC = CloneTH1(hDen, "h_leadRecoilJetMatchEff_den_clone_" + rKey);
                TH1* numC = CloneTH1(hNum, "h_leadRecoilJetMatchEff_num_clone_" + rKey);
                EnsureSumw2(denC);
                EnsureSumw2(numC);

                std::unique_ptr<TH1> hEff( CloneTH1(numC, "h_leadRecoilJetMatchEff_" + rKey) );
                if (hEff)
                {
                  EnsureSumw2(hEff.get());
                  hEff->Divide(numC, denC, 1.0, 1.0, "B");
                }

                // Integrated diagnostic (keep this for sanity)
                // JES3-only convention: ignore truth pT < 15 GeV (underflow support for unfolding)
                const double minTruthPtJES3 = 15.0;
                int firstBinJES3 = 1;
                if (denC && denC->GetXaxis())
                {
                  firstBinJES3 = denC->GetXaxis()->FindBin(minTruthPtJES3 + 1e-6); // lands in 15-17 when 15 is an edge
                  if (firstBinJES3 < 1) firstBinJES3 = 1;
                }
                const int lastBinJES3 = (denC ? denC->GetNbinsX() : 0);

                const double denInt   = (denC   ? denC->Integral(firstBinJES3, lastBinJES3) : 0.0);
                const double numInt   = (numC   ? numC->Integral(firstBinJES3, lastBinJES3) : 0.0);
                const double missAInt = (hMissA ? hMissA->Integral(firstBinJES3, lastBinJES3) : 0.0);
                const double missBInt = (hMissB ? hMissB->Integral(firstBinJES3, lastBinJES3) : 0.0);

                const double effInt = SafeDivide(numInt,   denInt, 0.0);
                const double aFrac  = SafeDivide(missAInt, denInt, 0.0);
                const double bFrac  = SafeDivide(missBInt, denInt, 0.0);

                cout << "  [LeadRecoilJetMatchEff] Integrated over pTgammaTruth bins:\n";
                cout << "    DEN=" << denInt << "  NUM=" << numInt
                     << "  MISS_A=" << missAInt << "  MISS_B=" << missBInt << "\n";
                cout << "    eps=NUM/DEN=" << effInt
                     << "  MISS_A/DEN=" << aFrac
                     << "  MISS_B/DEN=" << bFrac
                     << "  (sum=" << (effInt + aFrac + bFrac) << ")\n";

                effSummary.push_back("LeadRecoilJetMatchEff (single-curve):");
                effSummary.push_back("  Output dir: " + D.dirXJProjEffLeadMatch);
                effSummary.push_back(TString::Format("  Integrated: eps=%.4f  missA/DEN=%.4f  missB/DEN=%.4f  sum=%.4f",
                  effInt, aFrac, bFrac, effInt + aFrac + bFrac).Data());

                  if (hEff)
                  {
                  TCanvas c(TString::Format("c_leadRecoilJetMatchEff_%s_%s", ds.label.c_str(), rKey.c_str()).Data(),
                            "c_leadRecoilJetMatchEff", 900, 720);
                  ApplyCanvasMargins1D(c);

                    hEff->SetTitle("");
                    hEff->GetXaxis()->SetTitle("p_{T}^{#gamma,truth} [GeV]");
                    hEff->GetYaxis()->SetTitle("#varepsilon_{lead} = P(reco jet matches truth lead recoil)");

                    // JES3-only convention: do not display truth pT < 15 GeV
                    const double minTruthPtJES3_plot = 15.0;
                    hEff->GetXaxis()->SetRangeUser(minTruthPtJES3_plot, hEff->GetXaxis()->GetXmax());

                    // Auto-scale y-axis to the data (avoid lots of empty space),
                    // but ONLY over the displayed truth pT region (>=15 GeV).
                    double yMin =  1e9;
                    double yMax = -1e9;
                    for (int ib = 1; ib <= hEff->GetNbinsX(); ++ib)
                    {
                      if (hEff->GetXaxis()->GetBinLowEdge(ib) < minTruthPtJES3_plot) continue;

                        const double v = hEff->GetBinContent(ib);
                        const double e = hEff->GetBinError(ib);
                        if (v <= 0.0) continue; // ignore empty/unfilled bins

                        yMin = std::min(yMin, v - e);
                        yMax = std::max(yMax, v + e);
                    }
                      if (yMin > yMax) { yMin = 0.0; yMax = 1.0; } // fallback
                      const double pad = 0.06;
                      const double yLo = std::max(0.0, yMin - pad);
                      const double yHi = std::min(1.05, yMax + pad);
                      hEff->GetYaxis()->SetRangeUser(yLo, yHi);

                    // Draw as a TGraphErrors with x-errors forced to 0:
                    // this guarantees "vertical-only" statistical errors and no histogram bin-width segments.
                    std::unique_ptr<TGraphErrors> gEff(new TGraphErrors());
                    gEff->SetName(TString::Format("g_leadRecoilJetMatchEff_%s", rKey.c_str()).Data());

                    int ip = 0;
                      for (int ib = 1; ib <= hEff->GetNbinsX(); ++ib)
                      {
                          const double y  = hEff->GetBinContent(ib);
                          const double ey = hEff->GetBinError(ib);

                          const double x  = hEff->GetXaxis()->GetBinCenter(ib);
                          gEff->SetPoint(ip, x, y);
                          gEff->SetPointError(ip, 0.0, ey); // xerr=0.0  -> NO horizontal error bars
                          ++ip;
                      }

                  // Error bars are drawn using the graph's LINE attributes.
                  gEff->SetLineWidth(2);
                  gEff->SetLineColor(1);

                  gEff->SetMarkerStyle(20);
                  gEff->SetMarkerSize(1.10);
                  gEff->SetMarkerColor(1);

                  // Use the histogram only as an axis frame
                  hEff->SetLineWidth(0);
                  hEff->SetMarkerSize(0);
                  hEff->Draw("AXIS");

                      // Draw points + vertical errors only:
                      //  P = markers, E = error bars, Z = no end-caps
                      gEff->Draw("PEZ SAME");

                  // Optional reference line at 1 (start at 15 GeV)
                  const double xMin = minTruthPtJES3_plot;
                  const double xMax = hEff->GetXaxis()->GetXmax();
                  TLine one(xMin, 1.0, xMax, 1.0);
                  one.SetLineStyle(2);
                  one.SetLineWidth(2);
                  one.Draw("same");

                  const double Rval =
                        (rKey == "r02") ? 0.2 :
                        (rKey == "r04") ? 0.4 :
                        (rKey == "r06") ? 0.6 : -1.0;

                  const string bbLabel = B2BLabel();
                  const double jetPtMin_GeV = static_cast<double>(kJetPtMin);

                  DrawLatexLines(0.14, 0.92,
                        {TString::Format("Lead recoil-jet match efficiency (R = %.1f)", Rval).Data()},
                        0.030, 0.040
                  );
                  DrawLatexLines(0.14, 0.86, DefaultHeaderLines(ds), 0.030, 0.040);

                  TLatex tCuts;
                  tCuts.SetNDC(true);
                  tCuts.SetTextFont(42);
                  tCuts.SetTextAlign(33);
                  tCuts.SetTextSize(0.028);
                  tCuts.DrawLatex(0.92, 0.78, TString::Format("|#Delta#phi(#gamma,jet)| > %s", bbLabel.c_str()).Data());
                  tCuts.DrawLatex(0.92, 0.70, TString::Format("Reco (p_{T}^{jet} > %.0f GeV)", jetPtMin_GeV).Data());


                  const string outPng = JoinPath(
                      D.dirXJProjEffLeadMatch,
                      TString::Format("leadRecoilJetMatchEff_vs_pTgammaTruth_%s.png", rKey.c_str()).Data()
                  );
                  SaveCanvas(c, outPng);

                  cout << "    Wrote: " << outPng << "\n";
                  effSummary.push_back("  - " + outPng);

                  // ------------------------------------------------------------
                  // NEW: overlay efficiency curves for r02 and r04 on one canvas
                  //   - r02: red filled circles
                  //   - r04: blue filled circles
                  //
                  // Output (exact):
                  //   <...>/xJ_fromJES3/Efficiency/LeadTruthRecoilMatch/
                  //     leadRecoilJetMatchEff_vs_pTgammaTruth_overlay_r02_r04.png
                  // ------------------------------------------------------------
                  if (rKey == "r04")
                  {
                      const std::string rk02 = "r02";
                      const std::string rk04 = "r04";

                      const string hDen02Name = "h_leadTruthRecoilMatch_den_pTgammaTruth_" + rk02;
                      const string hNum02Name = "h_leadTruthRecoilMatch_num_pTgammaTruth_" + rk02;

                      const string hDen04Name = "h_leadTruthRecoilMatch_den_pTgammaTruth_" + rk04;
                      const string hNum04Name = "h_leadTruthRecoilMatch_num_pTgammaTruth_" + rk04;

                      TH1* hDen02 = GetObj<TH1>(ds, hDen02Name, false, false, false);
                      TH1* hNum02 = GetObj<TH1>(ds, hNum02Name, false, false, false);
                      TH1* hDen04 = GetObj<TH1>(ds, hDen04Name, false, false, false);
                      TH1* hNum04 = GetObj<TH1>(ds, hNum04Name, false, false, false);

                      if (!hDen02 || !hNum02 || !hDen04 || !hNum04)
                      {
                        cout << ANSI_BOLD_YEL << "  [LeadRecoilJetMatchEff Overlay] Missing r02/r04 inputs; skipping overlay.\n" << ANSI_RESET;
                        if (!hDen02) cout << "    - MISSING: " << hDen02Name << "\n";
                        if (!hNum02) cout << "    - MISSING: " << hNum02Name << "\n";
                        if (!hDen04) cout << "    - MISSING: " << hDen04Name << "\n";
                        if (!hNum04) cout << "    - MISSING: " << hNum04Name << "\n";
                      }
                      else
                      {
                        // Build binomial efficiencies: NUM/DEN (do NOT modify on-file hists)
                        TH1* den02C = CloneTH1(hDen02, "h_leadRecoilJetMatchEffOverlay_den02_clone");
                        TH1* num02C = CloneTH1(hNum02, "h_leadRecoilJetMatchEffOverlay_num02_clone");
                        TH1* den04C = CloneTH1(hDen04, "h_leadRecoilJetMatchEffOverlay_den04_clone");
                        TH1* num04C = CloneTH1(hNum04, "h_leadRecoilJetMatchEffOverlay_num04_clone");
                        EnsureSumw2(den02C); EnsureSumw2(num02C);
                        EnsureSumw2(den04C); EnsureSumw2(num04C);

                        std::unique_ptr<TH1> hEff02( CloneTH1(num02C, "h_leadRecoilJetMatchEffOverlay_eff02") );
                        std::unique_ptr<TH1> hEff04( CloneTH1(num04C, "h_leadRecoilJetMatchEffOverlay_eff04") );
                        if (hEff02) { EnsureSumw2(hEff02.get()); hEff02->Divide(num02C, den02C, 1.0, 1.0, "B"); }
                        if (hEff04) { EnsureSumw2(hEff04.get()); hEff04->Divide(num04C, den04C, 1.0, 1.0, "B"); }

                          // JES3-only convention: do not display truth pT < 15 GeV
                          const double minTruthPtJES3_plot = 15.0;

                          auto MakeGraphNoXerr = [&](TH1* h, int mStyle, double mSize, int lColor, int mColor) -> std::unique_ptr<TGraphErrors>
                          {
                            std::unique_ptr<TGraphErrors> g(new TGraphErrors());
                            int ip = 0;
                            for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
                            {
                              if (h->GetXaxis()->GetBinLowEdge(ib) < minTruthPtJES3_plot) continue;

                              const double y  = h->GetBinContent(ib);
                              const double ey = h->GetBinError(ib);
                              if (y <= 0.0) continue;

                              const double x  = h->GetXaxis()->GetBinCenter(ib);
                              g->SetPoint(ip, x, y);
                              g->SetPointError(ip, 0.0, ey); // vertical-only errors
                              ++ip;
                            }
                            g->SetLineWidth(2);
                            g->SetLineColor(lColor);
                            g->SetMarkerStyle(mStyle);
                            g->SetMarkerSize(mSize);
                            g->SetMarkerColor(mColor);
                            return g;
                          };

                        auto g02 = (hEff02 ? MakeGraphNoXerr(hEff02.get(), 20, 1.05, 2, 2) : nullptr); // r02 red
                        auto g04 = (hEff04 ? MakeGraphNoXerr(hEff04.get(), 20, 1.05, 4, 4) : nullptr); // r04 blue

                        if (g02 && g04)
                        {
                          // Compute y-range over BOTH curves
                          double yMin =  1e9;
                          double yMax = -1e9;
                          auto AccumRange = [&](TH1* h)
                          {
                            for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
                            {
                                if (h->GetXaxis()->GetBinLowEdge(ib) < minTruthPtJES3_plot) continue;
                                const double v = h->GetBinContent(ib);
                                const double e = h->GetBinError(ib);
                                if (v <= 0.0) continue;
                              yMin = std::min(yMin, v - e);
                              yMax = std::max(yMax, v + e);
                            }
                          };
                          AccumRange(hEff02.get());
                          AccumRange(hEff04.get());
                          if (yMin > yMax) { yMin = 0.0; yMax = 1.0; }
                          const double pad = 0.06;
                          const double yLo = std::max(0.0, yMin - pad);
                          const double yHi = std::min(1.05, yMax + pad);

                          TCanvas cOv(TString::Format("c_leadRecoilJetMatchEff_overlay_%s", ds.label.c_str()).Data(),
                                      "c_leadRecoilJetMatchEff_overlay", 900, 720);
                          ApplyCanvasMargins1D(cOv);

                          // Use r04 hist as an axis frame
                          hEff04->SetTitle("");
                          hEff04->GetXaxis()->SetTitle("p_{T}^{#gamma,truth} [GeV]");
                          hEff04->GetYaxis()->SetTitle("#varepsilon_{lead} = P(reco jet matches truth lead recoil)");
                          hEff04->GetYaxis()->SetRangeUser(yLo, yHi);
                          hEff04->GetXaxis()->SetRangeUser(minTruthPtJES3_plot, hEff04->GetXaxis()->GetXmax());
                          hEff04->SetLineWidth(0);
                          hEff04->SetMarkerSize(0);
                          hEff04->Draw("AXIS");

                          // Draw both curves with stat error bars
                          g02->Draw("PE SAME");
                          g04->Draw("PE SAME");

                          // Reference line at 1
                          const double xMin = minTruthPtJES3_plot;
                          const double xMax = hEff04->GetXaxis()->GetXmax();
                          TLine one(xMin, 1.0, xMax, 1.0);
                          one.SetLineStyle(2);
                          one.SetLineWidth(2);
                          one.Draw("same");

                          // Header + title
                          DrawLatexLines(0.14, 0.94, DefaultHeaderLines(ds), 0.030, 0.040);
                          DrawLatexLines(0.14, 0.86,
                            {"Lead recoil-jet match efficiency overlay (r02 vs r04)"},
                            0.030, 0.040
                          );

                          // Legend (smaller, shifted left)
                          TLegend leg(0.36, 0.76, 0.70, 0.90);
                          leg.SetTextFont(42);
                          leg.SetTextSize(0.028);
                          leg.SetFillStyle(0);
                          leg.SetBorderSize(0);
                          leg.AddEntry(g02.get(), "r02", "ep");
                          leg.AddEntry(g04.get(), "r04", "ep");
                          leg.Draw();

                          const string outOvPng = JoinPath(
                            D.dirXJProjEffLeadMatch,
                            "leadRecoilJetMatchEff_vs_pTgammaTruth_overlay_r02_r04.png"
                          );
                          SaveCanvas(cOv, outOvPng);

                          cout << "    Wrote: " << outOvPng << "\n";
                          effSummary.push_back("  - " + outOvPng);
                        }

                        delete den02C;
                        delete num02C;
                        delete den04C;
                        delete num04C;
                      }
                    }

                    // ------------------------------------------------------------
                    // NEW: missA and missB probabilities vs pTgammaTruth
                    //   fA = MISS_A / DEN  (matched, but not to the chosen/leading reco recoil jet)
                    //   fB = MISS_B / DEN  (no reco match)
                    //
                    // Output (exact):
                    //   <...>/<rKey>/xJ_fromJES3/Efficiency/LeadTruthRecoilMatch/
                    //     leadRecoilJetMatchMisses_vs_pTgammaTruth_<rKey>.png
                    // ------------------------------------------------------------
                    {
                      // Build binomial fraction histograms (do NOT modify on-file histograms)
                      TH1* missA_C = CloneTH1(hMissA, "h_leadTruthRecoilMatch_missA_cloneForFrac_" + rKey);
                      TH1* missB_C = CloneTH1(hMissB, "h_leadTruthRecoilMatch_missB_cloneForFrac_" + rKey);
                      EnsureSumw2(missA_C);
                      EnsureSumw2(missB_C);

                      std::unique_ptr<TH1> hFracA( CloneTH1(missA_C, "h_leadTruthRecoilMatch_fracA_" + rKey) );
                      std::unique_ptr<TH1> hFracB( CloneTH1(missB_C, "h_leadTruthRecoilMatch_fracB_" + rKey) );

                      if (hFracA && hFracB)
                      {
                        EnsureSumw2(hFracA.get());
                        EnsureSumw2(hFracB.get());
                        hFracA->Divide(missA_C, denC, 1.0, 1.0, "B");
                        hFracB->Divide(missB_C, denC, 1.0, 1.0, "B");

                          // JES3-only convention: do not display truth pT < 15 GeV
                          const double minTruthPtJES3_plot = 15.0;

                          // Auto-scale y-axis for misses only, but ONLY over displayed truth pT (>=15 GeV)
                          double yMinM =  1e9;
                          double yMaxM = -1e9;
                          for (int ib = 1; ib <= hFracA->GetNbinsX(); ++ib)
                          {
                            if (hFracA->GetXaxis()->GetBinLowEdge(ib) < minTruthPtJES3_plot) continue;

                            const double a  = hFracA->GetBinContent(ib);
                            const double ea = hFracA->GetBinError(ib);
                            const double b  = hFracB->GetBinContent(ib);
                            const double eb = hFracB->GetBinError(ib);

                            if (a > 0.0) { yMinM = std::min(yMinM, a - ea); yMaxM = std::max(yMaxM, a + ea); }
                            if (b > 0.0) { yMinM = std::min(yMinM, b - eb); yMaxM = std::max(yMaxM, b + eb); }
                          }
                          if (yMinM > yMaxM) { yMinM = 0.0; yMaxM = 0.20; } // fallback
                          const double padM = 0.02;
                          const double yLoM = 0.0;
                          const double yHiM = std::min(1.0, yMaxM + padM);

                          // Build TGraphErrors with xerr=0 (vertical-only statistical errors)
                          auto MakeGraphNoXerr =
                            [&](TH1* h, int mStyle, double mSize, int lColor, int mColor) -> std::unique_ptr<TGraphErrors>
                          {
                            std::unique_ptr<TGraphErrors> g(new TGraphErrors());
                            int ip = 0;
                            for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
                            {
                              if (h->GetXaxis()->GetBinLowEdge(ib) < minTruthPtJES3_plot) continue;

                              const double y  = h->GetBinContent(ib);
                              const double ey = h->GetBinError(ib);
                              if (y <= 0.0) continue;

                              const double x  = h->GetXaxis()->GetBinCenter(ib);

                              g->SetPoint(ip, x, y);
                              g->SetPointError(ip, 0.0, ey); // NO horizontal error bars
                              ++ip;
                            }
                            g->SetLineWidth(2);          // makes vertical errors visible
                            g->SetLineColor(lColor);
                            g->SetMarkerStyle(mStyle);
                            g->SetMarkerSize(mSize);
                            g->SetMarkerColor(mColor);
                            return g;
                          };

                          // Use FILLED CIRCLE markers for BOTH series (different colors distinguish A vs B)
                          auto gA = MakeGraphNoXerr(hFracA.get(), 20, 1.05, 2, 2); // red filled circles
                          auto gB = MakeGraphNoXerr(hFracB.get(), 20, 1.05, 4, 4); // blue filled circles

                          TCanvas cM(TString::Format("c_leadTruthMatchMisses_%s_%s", ds.label.c_str(), rKey.c_str()).Data(),
                                     "c_leadTruthMatchMisses", 900, 720);
                          ApplyCanvasMargins1D(cM);

                          // Use hFracA as the axis frame
                          hFracA->SetTitle("");
                          hFracA->GetXaxis()->SetTitle("p_{T}^{#gamma,truth} [GeV]");
                          hFracA->GetYaxis()->SetTitle("Miss probability (fraction of DEN)");
                          hFracA->GetXaxis()->SetRangeUser(minTruthPtJES3_plot, hFracA->GetXaxis()->GetXmax());
                          hFracA->GetYaxis()->SetRangeUser(yLoM, yHiM);

                          // Draw axis only
                          hFracA->SetLineWidth(0);
                          hFracA->SetMarkerSize(0);
                          hFracA->Draw("AXIS");

                        // Draw graphs WITH statistical error bars (xerr=0 from builder)
                        if (gA) gA->Draw("PE SAME");
                        if (gB) gB->Draw("PE SAME");

                        // 1) Legend: smaller + shifted LEFT so it stays fully on canvas
                        TLegend legM(0.38, 0.78, 0.73, 0.90);
                        legM.SetTextFont(42);
                        legM.SetTextSize(0.028);
                        legM.SetFillStyle(0);
                        legM.SetBorderSize(0);
                        if (gA) legM.AddEntry(gA.get(), "MissA / DEN (matched to non-leading reco recoil jet)", "ep");
                        if (gB) legM.AddEntry(gB.get(), "MissB / DEN (no reco match)", "ep");
                        legM.Draw();

                        // ------------------------------------------------------------
                        // NEW: Closure check block (middle-right)
                        //   DEN = NUM + MissA + MissB  =>  eps + fA + fB = 1
                        // ------------------------------------------------------------
                        {
                            // JES3-only convention: ignore truth pT < 15 GeV for closure readout (match displayed region)
                            const double minTruthPtJES3_plot = 15.0;
                            int firstBinJES3 = 1;
                            if (denC && denC->GetXaxis())
                            {
                              firstBinJES3 = denC->GetXaxis()->FindBin(minTruthPtJES3_plot + 1e-6);
                              if (firstBinJES3 < 1) firstBinJES3 = 1;
                            }
                            const int lastBinJES3 = (denC ? denC->GetNbinsX() : 0);

                            const double denIntC   = (denC    ? denC->Integral(firstBinJES3, lastBinJES3) : 0.0);
                            const double numIntC   = (numC    ? numC->Integral(firstBinJES3, lastBinJES3) : 0.0);
                            const double missAIntC = (missA_C ? missA_C->Integral(firstBinJES3, lastBinJES3) : 0.0);
                            const double missBIntC = (missB_C ? missB_C->Integral(firstBinJES3, lastBinJES3) : 0.0);

                          const double epsIntC = (denIntC > 0.0 ? numIntC   / denIntC : 0.0);
                          const double fAIntC  = (denIntC > 0.0 ? missAIntC / denIntC : 0.0);
                          const double fBIntC  = (denIntC > 0.0 ? missBIntC / denIntC : 0.0);
                          const double sumIntC = epsIntC + fAIntC + fBIntC;

                          double maxAbsDev = 0.0;
                          if (denC && numC && missA_C && missB_C)
                          {
                              const int nb = denC->GetNbinsX();
                              for (int ib = firstBinJES3; ib <= nb; ++ib)
                            {
                              const double den = denC->GetBinContent(ib);
                              if (den <= 0.0) continue;

                              const double sumBin =
                                (numC->GetBinContent(ib) +
                                 missA_C->GetBinContent(ib) +
                                 missB_C->GetBinContent(ib)) / den;

                              maxAbsDev = std::max(maxAbsDev, std::fabs(1.0 - sumBin));
                            }
                          }

                          TLatex tc;
                          tc.SetNDC(true);
                          tc.SetTextFont(42);
                          tc.SetTextAlign(13); // left/top

                          const double x0 = 0.55;
                          const double y0 = 0.62;
                          const double dy = 0.040;

                          tc.SetTextSize(0.030);
                          tc.DrawLatex(x0, y0, "#bf{Closure:}  #varepsilon_{lead} + f_{A} + f_{B} = 1");

                          tc.SetTextSize(0.028);
                          tc.DrawLatex(x0, y0 - 1.0*dy,
                            TString::Format("Integrated: %.3f + %.3f + %.3f = %.3f",
                                            epsIntC, fAIntC, fBIntC, sumIntC).Data());

                          tc.DrawLatex(x0, y0 - 2.0*dy,
                            TString::Format("max |1-(NUM+A+B)/DEN| = %.2e", maxAbsDev).Data());
                        }

                        // 2) TLatex header + subtitle: use DefaultHeaderLines(ds) so SIM sample label matches header toggles
                        {
                          const auto hdr = DefaultHeaderLines(ds);
                          const std::string dsLine = hdr.empty() ? std::string("Dataset: (unknown)") : hdr.front();

                          TLatex t;
                          t.SetNDC(true);
                          t.SetTextFont(42);
                          t.SetTextAlign(31); // right-aligned

                          // dataset header (top-right)
                          t.SetTextSize(0.034);
                          t.DrawLatex(0.89, 0.965, dsLine.c_str());

                          // plot label (just below header, top-right)
                          t.SetTextSize(0.032);
                          t.DrawLatex(0.89, 0.925,
                            TString::Format("Lead truth recoil miss probabilities (%s)", rKey.c_str()).Data()
                          );
                        }

                        const string outMissPng = JoinPath(
                          D.dirXJProjEffLeadMatch,
                          TString::Format("leadRecoilJetMatchMisses_vs_pTgammaTruth_%s.png", rKey.c_str()).Data()
                        );
                        SaveCanvas(cM, outMissPng);

                        cout << "    Wrote: " << outMissPng << "\n";
                        effSummary.push_back("  - " + outMissPng);
                      }

                        delete missA_C;
                        delete missB_C;

                        // ------------------------------------------------------------
                        // NEW: MissA subtype probabilities vs pTgammaTruth
                        //   fA1 = MissA1 / DEN
                        //   fA2 = MissA2 / DEN
                        //
                        // Output:
                        //   leadRecoilJetMatchMissA_subtypes_vs_pTgammaTruth_<rKey>.png
                        // ------------------------------------------------------------
                        if (hMissA1 && hMissA2)
                        {
                          TH1* missA1_C = CloneTH1(hMissA1, "h_leadTruthRecoilMatch_missA1_cloneForFrac_" + rKey);
                          TH1* missA2_C = CloneTH1(hMissA2, "h_leadTruthRecoilMatch_missA2_cloneForFrac_" + rKey);
                          EnsureSumw2(missA1_C);
                          EnsureSumw2(missA2_C);

                          std::unique_ptr<TH1> hFracA1( CloneTH1(missA1_C, "h_leadTruthRecoilMatch_fracA1_" + rKey) );
                          std::unique_ptr<TH1> hFracA2( CloneTH1(missA2_C, "h_leadTruthRecoilMatch_fracA2_" + rKey) );

                          if (hFracA1 && hFracA2)
                          {
                            EnsureSumw2(hFracA1.get());
                            EnsureSumw2(hFracA2.get());
                            hFracA1->Divide(missA1_C, denC, 1.0, 1.0, "B");
                            hFracA2->Divide(missA2_C, denC, 1.0, 1.0, "B");

                            // JES3-only convention: do not display truth pT < 15 GeV
                            const double minTruthPtJES3_plot = 15.0;

                            // Auto-scale y-axis over displayed truth pT (>=15 GeV)
                            double yMinS =  1e9;
                            double yMaxS = -1e9;
                            for (int ib = 1; ib <= hFracA1->GetNbinsX(); ++ib)
                            {
                              if (hFracA1->GetXaxis()->GetBinLowEdge(ib) < minTruthPtJES3_plot) continue;

                              const double a1  = hFracA1->GetBinContent(ib);
                              const double ea1 = hFracA1->GetBinError(ib);
                              const double a2  = hFracA2->GetBinContent(ib);
                              const double ea2 = hFracA2->GetBinError(ib);

                              if (a1 > 0.0) { yMinS = std::min(yMinS, a1 - ea1); yMaxS = std::max(yMaxS, a1 + ea1); }
                              if (a2 > 0.0) { yMinS = std::min(yMinS, a2 - ea2); yMaxS = std::max(yMaxS, a2 + ea2); }
                            }
                            if (yMinS > yMaxS) { yMinS = 0.0; yMaxS = 0.20; } // fallback
                            const double padS = 0.02;
                            const double yLoS = 0.0;
                            const double yHiS = std::min(1.0, yMaxS + padS);

                            auto MakeGraphNoXerr =
                              [&](TH1* h, int mStyle, double mSize, int lColor, int mColor) -> std::unique_ptr<TGraphErrors>
                            {
                              std::unique_ptr<TGraphErrors> g(new TGraphErrors());
                              int ip = 0;
                              for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
                              {
                                if (h->GetXaxis()->GetBinLowEdge(ib) < minTruthPtJES3_plot) continue;

                                const double y  = h->GetBinContent(ib);
                                const double ey = h->GetBinError(ib);
                                if (y <= 0.0) continue;

                                const double x  = h->GetXaxis()->GetBinCenter(ib);

                                g->SetPoint(ip, x, y);
                                g->SetPointError(ip, 0.0, ey); // NO horizontal error bars
                                ++ip;
                              }

                              g->SetMarkerStyle(mStyle);
                              g->SetMarkerSize(mSize);
                              g->SetLineColor(lColor);
                              g->SetMarkerColor(mColor);
                              g->SetLineWidth(2);
                              return g;
                            };

                            auto gA1 = MakeGraphNoXerr(hFracA1.get(), 20, 1.05, 2, 2); // red
                            auto gA2 = MakeGraphNoXerr(hFracA2.get(), 20, 1.05, 4, 4); // blue

                            TCanvas cS(TString::Format("c_leadTruthMatchMissA_Subtypes_%s_%s", ds.label.c_str(), rKey.c_str()).Data(),
                                       "c_leadTruthMatchMissA_Subtypes", 900, 720);
                            ApplyCanvasMargins1D(cS);

                            // Axis frame: use hFracA1
                            hFracA1->SetTitle("");
                            hFracA1->GetXaxis()->SetTitle("p_{T}^{#gamma,truth} [GeV]");
                            hFracA1->GetYaxis()->SetTitle("MissA subtype probability (fraction of DEN)");
                            hFracA1->GetXaxis()->SetRangeUser(minTruthPtJES3_plot, hFracA1->GetXaxis()->GetXmax());
                            hFracA1->GetYaxis()->SetRangeUser(yLoS, yHiS);

                            hFracA1->SetLineWidth(0);
                            hFracA1->SetMarkerSize(0);
                            hFracA1->Draw("AXIS");

                            if (gA1) gA1->Draw("PE SAME");
                            if (gA2) gA2->Draw("PE SAME");

                            TLegend legS(0.38, 0.78, 0.73, 0.90);
                            legS.SetTextFont(42);
                            legS.SetTextSize(0.028);
                            legS.SetFillStyle(0);
                            legS.SetBorderSize(0);
                            if (gA1) legS.AddEntry(gA1.get(), "MissA1 / DEN (truth-matched reco passes recoil cuts)", "ep");
                            if (gA2) legS.AddEntry(gA2.get(), "MissA2 / DEN (truth-matched reco fails recoil cuts)", "ep");
                            legS.Draw();

                            {
                              const auto hdr = DefaultHeaderLines(ds);
                              const std::string dsLine = hdr.empty() ? std::string("Dataset: (unknown)") : hdr.front();

                              TLatex t;
                              t.SetNDC(true);
                              t.SetTextFont(42);
                              t.SetTextAlign(31); // right-aligned

                              t.SetTextSize(0.034);
                              t.DrawLatex(0.89, 0.965, dsLine.c_str());

                              t.SetTextSize(0.032);
                              t.DrawLatex(0.89, 0.925,
                                TString::Format("MissA subtype probabilities (%s)", rKey.c_str()).Data()
                              );
                            }

                            const string outMissASubPng = JoinPath(
                              D.dirXJProjEffLeadMatch,
                              TString::Format("leadRecoilJetMatchMissA_subtypes_vs_pTgammaTruth_%s.png", rKey.c_str()).Data()
                            );
                            SaveCanvas(cS, outMissASubPng);

                            cout << "    Wrote: " << outMissASubPng << "\n";
                            effSummary.push_back("  - " + outMissASubPng);
                          }

                          delete missA1_C;
                          delete missA2_C;
                        }

                        // ------------------------------------------------------------
                        // NEW: MissA composition (within MissA): MissA1/MissA and MissA2/MissA
                        //
                        // Output:
                        //   leadRecoilJetMatchMissA_composition_vs_pTgammaTruth_<rKey>.png
                        // ------------------------------------------------------------
                        if (hMissA && hMissA1 && hMissA2)
                        {
                          TH1* missA_C2  = CloneTH1(hMissA,  "h_leadTruthRecoilMatch_missA_cloneForComp_"  + rKey);
                          TH1* missA1_C2 = CloneTH1(hMissA1, "h_leadTruthRecoilMatch_missA1_cloneForComp_" + rKey);
                          TH1* missA2_C2 = CloneTH1(hMissA2, "h_leadTruthRecoilMatch_missA2_cloneForComp_" + rKey);
                          EnsureSumw2(missA_C2);
                          EnsureSumw2(missA1_C2);
                          EnsureSumw2(missA2_C2);

                          std::unique_ptr<TH1> hCompA1( CloneTH1(missA1_C2, "h_leadTruthRecoilMatch_compA1_" + rKey) );
                          std::unique_ptr<TH1> hCompA2( CloneTH1(missA2_C2, "h_leadTruthRecoilMatch_compA2_" + rKey) );

                          if (hCompA1 && hCompA2)
                          {
                            EnsureSumw2(hCompA1.get());
                            EnsureSumw2(hCompA2.get());
                            hCompA1->Divide(missA1_C2, missA_C2, 1.0, 1.0, "B");
                            hCompA2->Divide(missA2_C2, missA_C2, 1.0, 1.0, "B");

                            const double minTruthPtJES3_plot = 15.0;

                            auto MakeGraphNoXerr =
                              [&](TH1* h, int mStyle, double mSize, int lColor, int mColor) -> std::unique_ptr<TGraphErrors>
                            {
                              std::unique_ptr<TGraphErrors> g(new TGraphErrors());
                              int ip = 0;
                              for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
                              {
                                if (h->GetXaxis()->GetBinLowEdge(ib) < minTruthPtJES3_plot) continue;

                                const double y  = h->GetBinContent(ib);
                                const double ey = h->GetBinError(ib);
                                if (y < 0.0) continue;

                                const double x  = h->GetXaxis()->GetBinCenter(ib);

                                g->SetPoint(ip, x, y);
                                g->SetPointError(ip, 0.0, ey);
                                ++ip;
                              }

                              g->SetMarkerStyle(mStyle);
                              g->SetMarkerSize(mSize);
                              g->SetLineColor(lColor);
                              g->SetMarkerColor(mColor);
                              g->SetLineWidth(2);
                              return g;
                            };

                            auto gC1 = MakeGraphNoXerr(hCompA1.get(), 20, 1.05, 2, 2); // red
                            auto gC2 = MakeGraphNoXerr(hCompA2.get(), 20, 1.05, 4, 4); // blue

                            TCanvas cC(TString::Format("c_leadTruthMatchMissA_Comp_%s_%s", ds.label.c_str(), rKey.c_str()).Data(),
                                       "c_leadTruthMatchMissA_Comp", 900, 720);
                            ApplyCanvasMargins1D(cC);

                            // Axis frame: reuse hCompA1
                            hCompA1->SetTitle("");
                            hCompA1->GetXaxis()->SetTitle("p_{T}^{#gamma,truth} [GeV]");
                            hCompA1->GetYaxis()->SetTitle("Composition within MissA");
                            hCompA1->GetXaxis()->SetRangeUser(minTruthPtJES3_plot, hCompA1->GetXaxis()->GetXmax());
                            hCompA1->GetYaxis()->SetRangeUser(0.0, 1.05);

                            hCompA1->SetLineWidth(0);
                            hCompA1->SetMarkerSize(0);
                            hCompA1->Draw("AXIS");

                            if (gC1) gC1->Draw("PE SAME");
                            if (gC2) gC2->Draw("PE SAME");

                            TLegend legC(0.38, 0.78, 0.73, 0.90);
                            legC.SetTextFont(42);
                            legC.SetTextSize(0.028);
                            legC.SetFillStyle(0);
                            legC.SetBorderSize(0);
                            if (gC1) legC.AddEntry(gC1.get(), "MissA1 / MissA", "ep");
                            if (gC2) legC.AddEntry(gC2.get(), "MissA2 / MissA", "ep");
                            legC.Draw();

                            {
                              const auto hdr = DefaultHeaderLines(ds);
                              const std::string dsLine = hdr.empty() ? std::string("Dataset: (unknown)") : hdr.front();

                              TLatex t;
                              t.SetNDC(true);
                              t.SetTextFont(42);
                              t.SetTextAlign(31);

                              t.SetTextSize(0.034);
                              t.DrawLatex(0.89, 0.965, dsLine.c_str());

                              t.SetTextSize(0.032);
                              t.DrawLatex(0.89, 0.925,
                                TString::Format("MissA composition (%s)", rKey.c_str()).Data()
                              );
                            }

                            const string outMissACompPng = JoinPath(
                              D.dirXJProjEffLeadMatch,
                              TString::Format("leadRecoilJetMatchMissA_composition_vs_pTgammaTruth_%s.png", rKey.c_str()).Data()
                            );
                            SaveCanvas(cC, outMissACompPng);

                            cout << "    Wrote: " << outMissACompPng << "\n";
                            effSummary.push_back("  - " + outMissACompPng);
                          }

                          delete missA_C2;
                          delete missA1_C2;
                          delete missA2_C2;
                        }
                      }

                }
                else
                {
                  cout << ANSI_BOLD_YEL << "  [LeadRecoilJetMatchEff] Failed to build efficiency curve; skipping plot.\n" << ANSI_RESET;
                  effSummary.push_back("  - FAILED to build efficiency curve (skipped)");
                }

                // ---------------------------------------------------------------------------
                // (1b) NEW: LeadTruthRecoilMatch diagnostics (SIM-only)
                //
                // Inputs (filled online in RecoilJets.cc, split by NUM/MissA/MissB):
                //   (A1) h2_leadTruthRecoilMatch_pTrecoJet1_vs_pTtruthLead_{num,missA,missB}_<rKey>
                //   (A2) h2_leadTruthRecoilMatch_pTrecoJet1_vs_pTrecoTruthMatch_{num,missA}_<rKey>
                //   (B3) h2_leadTruthRecoilMatch_dphiRecoJet1_{num,missA,missB}_pTgammaTruth_<rKey>
                //   (B4) h2_leadTruthRecoilMatch_dRRecoJet1_vs_truthLead_{num,missA,missB}_pTgammaTruth_<rKey>
                //   (C5) h2_leadTruthRecoilMatch_xJRecoJet1_vs_dphiRecoJet1_{num,missA,missB}_<rKey>
                //
                // Outputs:
                //   <...>/<rKey>/xJ_fromJES3/Efficiency/LeadTruthRecoilMatch/Diagnostics/
                //     Pt/*.png, Angle/*.png, XJvsDphi/*.png
                // ---------------------------------------------------------------------------
                {
                  const string dirDiag = JoinPath(D.dirXJProjEffLeadMatch, "Diagnostics");
                  const string dirPt   = JoinPath(dirDiag, "Pt");
                  const string dirAng  = JoinPath(dirDiag, "Angle");
                  const string dirXJ2D = JoinPath(dirDiag, "XJvsDphi");
                  const string dirED   = JoinPath(dirDiag, "eventDisplay");
                  EnsureDir(dirDiag);
                  EnsureDir(dirPt);
                  EnsureDir(dirAng);
                  EnsureDir(dirXJ2D);
                  EnsureDir(dirED);

                  effSummary.push_back("LeadTruthRecoilMatch Diagnostics:");
                  effSummary.push_back("  Output dir: " + dirDiag);

                  // ---------------------------------------------------------------------------
                  // (ED) EventDisplay (offline) from EventDisplayTree
                  //
                  // Purpose
                  //   - Read the EventDisplayTree payload written by RecoilJets (diagnostics-only)
                  //   - Render tower E_T maps in a “save3D”-style presentation (TH2 COLZ + LEGO2)
                  //   - Select one pseudo-random entry per category for this rKey:
                  //       NUM  (truth-matched)
                  //       MissA (wrong jet)
                  //       MissB (no reco match)
                  //
                  // Output (per rKey):
                  //   <...>/<rKey>/xJ_fromJES3/Efficiency/LeadTruthRecoilMatch/Diagnostics/eventDisplay/
                  //     NUM/eventDisplay_NUM_<rKey>_runXXXX_evtYYYYYY.png
                  //     MissA/eventDisplay_MissA_<rKey>_runXXXX_evtYYYYYY.png
                  //     MissB/eventDisplay_MissB_<rKey>_runXXXX_evtYYYYYY.png
                  //
                  // NOTE
                  //   - This block is analysis-only; it does not alter selection/physics logic.
                  //   - Added: strong branch/type/bind diagnostics + payload summaries + “why blank?” clues.
                  // ---------------------------------------------------------------------------
                  {
                      const string dirED      = JoinPath(dirDiag, "eventDisplay");
                      const string dirED_NUM  = JoinPath(dirED, "NUM");
                      const string dirED_MA   = JoinPath(dirED, "MissA");
                      const string dirED_MB   = JoinPath(dirED, "MissB");
                      const string dirED_MTJ  = JoinPath(dirED, "MatchedTruthJets");
                      EnsureDir(dirED_NUM);
                      EnsureDir(dirED_MA);
                      EnsureDir(dirED_MB);
                      EnsureDir(dirED_MTJ);

                      auto PrintED =
                        [&](const string& tag, const string& msg, const string& color)->void
                      {
                        cout << color << tag << ANSI_RESET << " " << msg << "\n";
                      };

                      auto PrintEDKV =
                        [&](const string& tag,
                            const string& key,
                            const string& val,
                            const string& color)->void
                      {
                        cout << color << tag << ANSI_RESET << " " << key << "=" << val << "\n";
                      };

                      TTree* tED = (ds.topDir ? dynamic_cast<TTree*>(ds.topDir->Get("EventDisplayTree")) : nullptr);
                      if (!tED && ds.file)
                      {
                        tED = dynamic_cast<TTree*>(ds.file->Get("EventDisplayTree"));
                      }

                      if (!tED)
                      {
                        PrintED("  [EventDisplay][MISSING]",
                                "EventDisplayTree not found (searched ds.topDir then ds.file). topDirName=\"" + ds.topDirName + "\"",
                                ANSI_BOLD_YEL);
                        effSummary.push_back("  EventDisplayTree: MISSING in input ROOT (no eventDisplay PNGs generated)");
                      }
                      else
                      {
                        PrintED("  [EventDisplay][FOUND]",
                                "Tree=\"" + string(tED->GetName()) + "\" entries=" + std::to_string((long long)tED->GetEntries()) +
                                  "  outDir=" + dirED + "  (topDirName=\"" + ds.topDirName + "\")",
                                ANSI_BOLD_CYN);

                        cout << ANSI_BOLD_RED
                             << "  [EventDisplay][TOTAL ENTRIES] "
                             << ANSI_RESET
                             << "inFile=\"" << ds.inFilePath << "\""
                             << "  tree=\"" << tED->GetName() << "\""
                             << "  entries=" << (long long)tED->GetEntries()
                             << "\n";

                        effSummary.push_back("  EventDisplayTree: generating eventDisplay PNGs in " + dirED);

                        // -------------------------------------------------------------------
                        // Branch buffers
                        // -------------------------------------------------------------------
                        int b_run = 0;
                        int b_evt = 0;
                        float b_vz = 0.0f;

                        std::string* b_rKey = nullptr;
                        int b_cat = -1;

                        float b_ptGammaTruth = 0.0f;

                        float b_sel_eta = 0.0f;
                        float b_sel_phi = 0.0f;

                        float b_best_eta = 0.0f;
                        float b_best_phi = 0.0f;

                        float b_truth_eta = 0.0f;
                        float b_truth_phi = 0.0f;

                        std::vector<float>* b_sel_etaTower = nullptr;
                        std::vector<float>* b_sel_phiTower = nullptr;
                        std::vector<float>* b_sel_etTower  = nullptr;

                        std::vector<float>* b_best_etaTower = nullptr;
                        std::vector<float>* b_best_phiTower = nullptr;
                        std::vector<float>* b_best_etTower  = nullptr;

                        // -------------------------------------------------------------------
                        // (ED) Branch existence/type checks + bind return codes
                        // -------------------------------------------------------------------
                        auto DumpBranchInfo =
                          [&](const string& brName)->bool
                        {
                          TBranch* br = tED->GetBranch(brName.c_str());
                          if (!br)
                          {
                            PrintED("    [ED][BRANCH][MISSING]", "\"" + brName + "\"", ANSI_BOLD_YEL);
                            return false;
                          }

                          const string cls   = (br->GetClassName() ? br->GetClassName() : "");
                          const string title = (br->GetTitle() ? br->GetTitle() : "");
                          cout << ANSI_BOLD_GRN << "    [ED][BRANCH][OK]" << ANSI_RESET
                               << " \"" << brName << "\""
                               << "  class=\"" << cls << "\""
                               << "  title=\"" << title << "\"\n";
                          return true;
                        };

                        auto BindBranch =
                          [&](const string& brName, void* addr)->int
                        {
                          const int rc = tED->SetBranchAddress(brName.c_str(), addr);
                          if (rc < 0)
                          {
                            cout << ANSI_BOLD_RED << "    [ED][BIND][FAIL]" << ANSI_RESET
                                 << " rc=" << rc << "  \"" << brName << "\"\n";
                          }
                          else
                          {
                            cout << ANSI_BOLD_GRN << "    [ED][BIND][OK]" << ANSI_RESET
                                 << "   rc=" << rc << "  \"" << brName << "\"\n";
                          }
                          return rc;
                        };

                        // Require the branches we depend on (so the failure mode is explicit)
                        const std::vector<string> reqBranches = {
                          "run","evt","vz","rKey","cat","ptGammaTruth",
                          "sel_eta","sel_phi","best_eta","best_phi","truthLead_eta","truthLead_phi",
                          "sel_etaTower","sel_phiTower","sel_etTower",
                          "best_etaTower","best_phiTower","best_etTower"
                        };

                        bool allPresent = true;
                        for (const auto& b : reqBranches)
                        {
                          allPresent = DumpBranchInfo(b) && allPresent;
                        }

                        // Bind (even if some are missing; rc will tell us exactly what broke)
                        BindBranch("run", &b_run);
                        BindBranch("evt", &b_evt);
                        BindBranch("vz",  &b_vz);

                        BindBranch("rKey", &b_rKey);
                        BindBranch("cat",  &b_cat);

                        BindBranch("ptGammaTruth", &b_ptGammaTruth);

                        BindBranch("sel_eta", &b_sel_eta);
                        BindBranch("sel_phi", &b_sel_phi);

                        BindBranch("best_eta", &b_best_eta);
                        BindBranch("best_phi", &b_best_phi);

                        BindBranch("truthLead_eta", &b_truth_eta);
                        BindBranch("truthLead_phi", &b_truth_phi);

                        BindBranch("sel_etaTower", &b_sel_etaTower);
                        BindBranch("sel_phiTower", &b_sel_phiTower);
                        BindBranch("sel_etTower",  &b_sel_etTower);

                        BindBranch("best_etaTower", &b_best_etaTower);
                        BindBranch("best_phiTower", &b_best_phiTower);
                        BindBranch("best_etTower",  &b_best_etTower);

                        cout << "    [ED][PTR] sel(phi/eta/et)="
                             << (void*)b_sel_phiTower << "/" << (void*)b_sel_etaTower << "/" << (void*)b_sel_etTower
                             << "  best(phi/eta/et)="
                             << (void*)b_best_phiTower << "/" << (void*)b_best_etaTower << "/" << (void*)b_best_etTower
                             << "\n";

                        if (!allPresent)
                        {
                          PrintED("    [ED][WARN]",
                                  "One or more required branches are missing. EventDisplay payloads may be empty or unusable.",
                                  ANSI_BOLD_YEL);
                        }

                        // -------------------------------------------------------------------
                        // (ED) Sanity probe: entry 0, plus first entry matching rKey if available
                        // -------------------------------------------------------------------
                        auto PrintEntrySummary =
                          [&](const string& tag)->void
                        {
                          const string rk = (b_rKey ? *b_rKey : string("<null>"));
                          cout << ANSI_BOLD_CYN << "    [ED][ENTRY]" << ANSI_RESET
                               << " " << tag
                               << "  run=" << b_run
                               << "  evt=" << b_evt
                               << "  vz=" << b_vz
                               << "  rKey=" << rk
                               << "  cat=" << b_cat
                               << "  ptGammaTruth=" << b_ptGammaTruth
                               << "  sel(phi,eta)=(" << b_sel_phi << "," << b_sel_eta << ")"
                               << "  best(phi,eta)=(" << b_best_phi << "," << b_best_eta << ")"
                               << "  truth(phi,eta)=(" << b_truth_phi << "," << b_truth_eta << ")"
                               << "\n";

                          cout << "    [ED][ENTRY] " << tag
                               << "  sel sizes: phi="
                               << (b_sel_phiTower ? b_sel_phiTower->size() : 0)
                               << " eta=" << (b_sel_etaTower ? b_sel_etaTower->size() : 0)
                               << " et="  << (b_sel_etTower  ? b_sel_etTower->size()  : 0)
                               << "  best sizes: phi="
                               << (b_best_phiTower ? b_best_phiTower->size() : 0)
                               << " eta=" << (b_best_etaTower ? b_best_etaTower->size() : 0)
                               << " et="  << (b_best_etTower  ? b_best_etTower->size()  : 0)
                               << "\n";
                        };

                        if (tED->GetEntries() > 0)
                        {
                          const Long64_t nb0 = tED->GetEntry(0);
                          cout << "    [ED][SANITY] GetEntry(0) bytes=" << nb0 << "\n";
                          PrintEntrySummary("entry=0");
                        }

                        // If entry 0 isn’t for this rKey, find the first entry for this rKey and print it
                        bool printedFirstForThisRKey = false;
                        const Long64_t nEntScanMax = std::min<Long64_t>(tED->GetEntries(), 5000); // bounded scan
                        for (Long64_t i = 0; i < nEntScanMax; ++i)
                        {
                          tED->GetEntry(i);
                          if (!b_rKey) continue;
                          if (*b_rKey != rKey) continue;
                          cout << "    [ED][SANITY] First entry for rKey=\"" << rKey << "\" appears at i=" << i << "\n";
                          PrintEntrySummary("firstForRKey");
                          printedFirstForThisRKey = true;
                          break;
                        }
                        if (!printedFirstForThisRKey)
                        {
                          PrintED("    [ED][SANITY]",
                                  "No entries for this rKey found in first " + std::to_string((long long)nEntScanMax) + " entries (will still try full index build).",
                                  ANSI_BOLD_YEL);
                        }

                        // -------------------------------------------------------------------
                        // Utility: wrap phi to [-pi,pi]
                        // -------------------------------------------------------------------
                        auto WrapPhi =
                          [&](float phi)->float
                        {
                          while (phi <= -M_PI) phi += 2.0f*(float)M_PI;
                          while (phi >   M_PI) phi -= 2.0f*(float)M_PI;
                          return phi;
                        };

                        // -------------------------------------------------------------------
                        // (ED) Tower payload dump (what is actually in the vectors)
                        // -------------------------------------------------------------------
                        auto DumpTowerTriplet =
                          [&](const string& tag,
                              const std::vector<float>* vphi,
                              const std::vector<float>* veta,
                              const std::vector<float>* vet,
                              double xMin, double xMax,
                              double yMin, double yMax)->void
                        {
                          cout << ANSI_BOLD_CYN << "    [ED][TOWERS]" << ANSI_RESET << " " << tag << "\n";

                          if (!vphi || !veta || !vet)
                          {
                            PrintED("      [ED][TOWERS][MISSING]", "one or more tower vector pointers are null", ANSI_BOLD_RED);
                            return;
                          }

                          const size_t nPhi = vphi->size();
                          const size_t nEta = veta->size();
                          const size_t nEt  = vet->size();
                          const size_t n    = std::min(nPhi, std::min(nEta, nEt));

                          cout << "      sizes: phi=" << nPhi << " eta=" << nEta << " et=" << nEt
                               << "  using n=" << n << "\n";

                          if (n == 0)
                          {
                            PrintED("      [ED][TOWERS][EMPTY]",
                                    "vectors are empty (no towers to draw) — plots will be blank except markers/axes",
                                    ANSI_BOLD_YEL);
                            return;
                          }

                          float phiRawMin =  1e9f, phiRawMax = -1e9f;
                          float phiWMin   =  1e9f, phiWMax   = -1e9f;
                          float etaMin    =  1e9f, etaMax    = -1e9f;
                          float etMin     =  1e9f, etMax     = -1e9f;

                          double sumEt = 0.0;
                          size_t nEtPos = 0;
                          size_t nEtNeg = 0;
                          size_t nPhiOutRaw = 0;
                          size_t nEtaOut    = 0;

                          for (size_t i = 0; i < n; ++i)
                          {
                            const float phiRaw = (*vphi)[i];
                            const float phiW   = WrapPhi(phiRaw);
                            const float eta    = (*veta)[i];
                            const float et     = (*vet)[i];

                            phiRawMin = std::min(phiRawMin, phiRaw);
                            phiRawMax = std::max(phiRawMax, phiRaw);
                            phiWMin   = std::min(phiWMin,   phiW);
                            phiWMax   = std::max(phiWMax,   phiW);
                            etaMin    = std::min(etaMin, eta);
                            etaMax    = std::max(etaMax, eta);
                            etMin     = std::min(etMin, et);
                            etMax     = std::max(etMax, et);

                            sumEt += et;
                            if (et >= 0.0f) ++nEtPos;
                            else            ++nEtNeg;

                            if (phiRaw < xMin || phiRaw > xMax) ++nPhiOutRaw;
                            if (eta    < yMin || eta    > yMax) ++nEtaOut;
                          }

                          cout << "      et: sum=" << sumEt
                               << "  min=" << etMin
                               << "  max=" << etMax
                               << "  nPos=" << nEtPos
                               << "  nNeg=" << nEtNeg
                               << "\n";

                          cout << "      phi raw[min,max]=[" << phiRawMin << "," << phiRawMax << "]"
                               << "  outOfRange(raw)=" << nPhiOutRaw << "/" << n
                               << "  target=[" << xMin << "," << xMax << "]\n";

                          cout << "      phi wrapped[min,max]=[" << phiWMin << "," << phiWMax << "]"
                               << "  target=[" << xMin << "," << xMax << "]\n";

                          cout << "      eta[min,max]=[" << etaMin << "," << etaMax << "]"
                               << "  outOfRange(eta)=" << nEtaOut << "/" << n
                               << "  target=[" << yMin << "," << yMax << "]\n";

                          const size_t nTop = 8;
                          std::vector<size_t> top;
                          top.reserve(nTop);

                          for (size_t i = 0; i < n; ++i)
                          {
                            const float et = (*vet)[i];

                            size_t pos = 0;
                            while (pos < top.size() && et < (*vet)[top[pos]]) ++pos;

                            if (top.size() < nTop)
                            {
                              top.insert(top.begin() + pos, i);
                            }
                            else if (pos < nTop)
                            {
                              top.insert(top.begin() + pos, i);
                              top.pop_back();
                            }
                          }

                          cout << "      top towers (rank: phiRaw  phiW  eta  et)\n";
                          for (size_t ir = 0; ir < top.size(); ++ir)
                          {
                            const size_t i = top[ir];
                            const float phiRaw = (*vphi)[i];
                            const float phiW   = WrapPhi(phiRaw);
                            const float eta    = (*veta)[i];
                            const float et     = (*vet)[i];

                            cout << "        " << (ir + 1) << ": "
                                 << phiRaw << "  " << phiW << "  " << eta << "  " << et << "\n";
                          }
                        };

                        // Fill helper (phi-wrapped into [-pi,pi])
                        auto FillJetHist =
                          [&](TH2F* h,
                              const std::vector<float>* vphi,
                              const std::vector<float>* veta,
                              const std::vector<float>* vet)->void
                        {
                          if (!h || !vphi || !veta || !vet) return;
                          const size_t n = std::min(vphi->size(), std::min(veta->size(), vet->size()));
                          for (size_t i = 0; i < n; ++i)
                          {
                            const float phiW = WrapPhi((*vphi)[i]);
                            h->Fill(phiW, (*veta)[i], (*vet)[i]);
                          }
                        };

                        auto StyleHistLikeSave3D =
                          [&](TH2F* h)->void
                        {
                          if (!h) return;
                          h->SetStats(0);
                          h->SetContour(99);
                          h->SetTitle("");

                          h->GetXaxis()->SetTitle("#phi");
                          h->GetYaxis()->SetTitle("#eta");
                          h->GetZaxis()->SetTitle("Tower E_{T} [GeV]");

                          h->GetXaxis()->SetTitleOffset(1.6);
                          h->GetYaxis()->SetTitleOffset(2.0);
                          h->GetZaxis()->SetTitleOffset(1.4);
                        };

                        auto StylePadLikeSave3D =
                          [&](TPad* pad)->void
                        {
                          if (!pad) return;
                          pad->SetLeftMargin(0.14);
                          pad->SetBottomMargin(0.14);
                          pad->SetRightMargin(0.32);
                          pad->SetTopMargin(0.06);
                        };

                        auto AdjustPaletteLikeSave3D =
                          [&](TH2F* h)->void
                        {
                          if (!h) return;
                          gPad->Update();

                          // Local include (keeps this change self-contained)
                          #include <TPaletteAxis.h>

                          if (TPaletteAxis* pal = (TPaletteAxis*)h->GetListOfFunctions()->FindObject("palette"))
                          {
                            pal->SetX1NDC(0.72);
                            pal->SetX2NDC(0.86);
                            pal->SetY1NDC(0.15);
                            pal->SetY2NDC(0.90);
                          }
                        };

                        auto DrawMarkersAndHeader =
                          [&](float jetPhi, float jetEta,
                              float truthPhi, float truthEta,
                              const string& headerLine,
                              const string& subLine)->void
                        {
                          const float jetPhiW   = WrapPhi(jetPhi);
                          const float truthPhiW = WrapPhi(truthPhi);

                          // Markers
                          TMarker mReco(jetPhiW, jetEta, 29);
                          mReco.SetMarkerSize(1.8);
                          mReco.Draw();

                          TMarker mTruth(truthPhiW, truthEta, 24);
                          mTruth.SetMarkerSize(1.4);
                          mTruth.Draw();

                          // Header
                          TLatex tex;
                          tex.SetNDC();
                          tex.SetTextSize(0.040);
                          tex.DrawLatex(0.14, 0.96, headerLine.c_str());

                          tex.SetTextSize(0.034);
                          tex.DrawLatex(0.14, 0.91, subLine.c_str());
                        };

                        auto DrawPanelSave3DStyle =
                          [&](TPad* pad,
                              TH2F* h,
                              const string& drawOpt,
                              float jetPhi, float jetEta,
                              float truthPhi, float truthEta,
                              const string& headerLine,
                              const string& subLine)->void
                        {
                          if (!pad || !h) return;

                          pad->cd();
                          StylePadLikeSave3D(pad);
                          StyleHistLikeSave3D(h);

                          h->Draw(drawOpt.c_str());
                          AdjustPaletteLikeSave3D(h);
                          DrawMarkersAndHeader(jetPhi, jetEta, truthPhi, truthEta, headerLine, subLine);
                        };

                        // Build index lists by category for this rKey
                        std::vector<Long64_t> idxByCat[3];
                        const Long64_t nEnt = tED->GetEntries();

                        for (Long64_t ient = 0; ient < nEnt; ++ient)
                        {
                          tED->GetEntry(ient);
                          if (!b_rKey) continue;
                          if (*b_rKey != rKey) continue;
                          if (b_cat < 0 || b_cat > 2) continue;
                          idxByCat[b_cat].push_back(ient);
                        }

                        cout << ANSI_BOLD_CYN << "  [EventDisplay][INDEX]" << ANSI_RESET
                             << " rKey=\"" << rKey << "\""
                             << "  NUM=" << idxByCat[0].size()
                             << "  MissA=" << idxByCat[1].size()
                             << "  MissB=" << idxByCat[2].size()
                             << "  (tree entries=" << nEnt << ")"
                             << "\n";

                        auto PickIndex = [&](int icat)->Long64_t
                        {
                          if (icat < 0 || icat > 2) return -1;
                          if (idxByCat[icat].empty()) return -1;

                          // Pseudo-random but stable per run: seed mixes category size and total entries
                          const unsigned int seed =
                            0xC0FFEEu ^ (unsigned int)(idxByCat[icat].size() * 131u) ^ (unsigned int)(nEnt * 17u);

                          std::mt19937 rng(seed);
                          std::uniform_int_distribution<size_t> uni(0, idxByCat[icat].size() - 1);
                          return idxByCat[icat][uni(rng)];
                        };

                        auto SaveNUMorMissB =
                          [&](const string& catName, int icat, const string& outDir)->void
                        {
                          const Long64_t pick = PickIndex(icat);
                          if (pick < 0)
                          {
                            PrintED("  [EventDisplay][SKIP]",
                                    catName + " has no entries for rKey=\"" + rKey + "\" (icat=" + std::to_string(icat) + ")",
                                    ANSI_BOLD_YEL);
                            return;
                        }

                        tED->GetEntry(pick);

                        const string outPng =
                            JoinPath(outDir,
                                     "eventDisplay_" + catName + "_" + rKey +
                                     "_run" + std::to_string(b_run) +
                                     "_evt" + std::to_string(b_evt) + ".png");

                        cout << ANSI_BOLD_CYN << "  [EventDisplay][DO]" << ANSI_RESET
                               << " " << catName
                               << " pick=" << pick
                               << "  run=" << b_run
                               << "  evt=" << b_evt
                               << "  vz=" << b_vz
                               << "  ptGammaTruth=" << b_ptGammaTruth
                               << "  -> " << outPng
                               << "\n";

                        const float selPhiW     = WrapPhi(b_sel_phi);
                        const float truthPhiW   = WrapPhi(b_truth_phi);
                        const float dPhiSelTruth = WrapPhi(selPhiW - truthPhiW);

                        cout << "    [ED][KIN] sel(phi,eta)=(" << b_sel_phi << "," << b_sel_eta << ")  wrappedPhi=" << selPhiW << "\n";
                        cout << "    [ED][KIN] truth(phi,eta)=(" << b_truth_phi << "," << b_truth_eta << ")  wrappedPhi=" << truthPhiW << "\n";
                        cout << "    [ED][KIN] dphi(sel-truth) wrapped=" << dPhiSelTruth << "  |dphi|=" << std::fabs(dPhiSelTruth) << "\n";

                        TH2F hColz("hED_colz", "", 64, -M_PI, M_PI, 48, -1.1, 1.1);
                        TH2F hLego("hED_lego", "", 64, -M_PI, M_PI, 48, -1.1, 1.1);

                        DumpTowerTriplet(catName + string(": sel"),
                                           b_sel_phiTower, b_sel_etaTower, b_sel_etTower,
                                           hColz.GetXaxis()->GetXmin(), hColz.GetXaxis()->GetXmax(),
                                           hColz.GetYaxis()->GetXmin(), hColz.GetYaxis()->GetXmax());

                        FillJetHist(&hColz, b_sel_phiTower, b_sel_etaTower, b_sel_etTower);
                        FillJetHist(&hLego, b_sel_phiTower, b_sel_etaTower, b_sel_etTower);

                        cout << "    [ED][HIST] hColz: entries=" << hColz.GetEntries()
                               << " sumW=" << hColz.GetSumOfWeights()
                               << " maxBin=" << hColz.GetMaximum()
                               << " integral=" << hColz.Integral()
                               << "\n";

                        TCanvas c("cED", "", 1600, 800);
                        c.Divide(2, 1, 0.0, 0.0);

                        const string header = "EventDisplay " + catName + "  " + rKey;
                        const string sub    = "run " + std::to_string(b_run) +
                                                "  evt " + std::to_string(b_evt) +
                                                "  v_{z}=" + std::to_string((int)std::round(b_vz)) + " cm" +
                                                "  pT_{#gamma}^{truth}=" + std::to_string((int)std::round(b_ptGammaTruth)) + " GeV";

                        DrawPanelSave3DStyle((TPad*)c.cd(1), &hColz, "COLZ",  selPhiW, b_sel_eta, truthPhiW, b_truth_eta, header, sub);
                        DrawPanelSave3DStyle((TPad*)c.cd(2), &hLego, "LEGO2", selPhiW, b_sel_eta, truthPhiW, b_truth_eta, header, sub);

                        c.SaveAs(outPng.c_str());
                        PrintED("  [EventDisplay][WROTE]", outPng, ANSI_BOLD_GRN);

                        // -----------------------------------------------------------------
                        // (ED) Publication-ready 3D-only panel for truth-matched recoil jets (NUM)
                        // -----------------------------------------------------------------
                        if (icat == 0)
                        {
                              const string outPng3D =
                                JoinPath(dirED_MTJ,
                                         "matchedTruthJets3D_" + rKey +
                                         "_run" + std::to_string(b_run) +
                                         "_evt" + std::to_string(b_evt) + ".png");

                              TCanvas c3D("cED_MatchedTruthJets3D", "", 1400, 1000);

                              const string header3D = "Matched Recoil Truth Jets, Photon 10 + 20 GeV Sim";
                              const string sub3D    = rKey + "  run " + std::to_string(b_run) +
                                                      "  evt " + std::to_string(b_evt) +
                                                      "  v_{z}=" + std::to_string((int)std::round(b_vz)) + " cm" +
                                                      "  pT_{#gamma}^{truth}=" + std::to_string((int)std::round(b_ptGammaTruth)) + " GeV";

                              DrawPanelSave3DStyle((TPad*)c3D.cd(), &hLego, "LEGO2",
                                                   selPhiW, b_sel_eta, truthPhiW, b_truth_eta,
                                                   header3D, sub3D);

                              c3D.SaveAs(outPng3D.c_str());
                              PrintED("  [EventDisplay][WROTE]", outPng3D, ANSI_BOLD_GRN);
                            }
                        };

                        auto SaveMissA =
                          [&](const string& outDir)->void
                        {
                          const Long64_t pick = PickIndex(1);
                          if (pick < 0)
                          {
                            PrintED("  [EventDisplay][SKIP]",
                                    "MissA has no entries for rKey=\"" + rKey + "\"",
                                    ANSI_BOLD_YEL);
                            return;
                          }

                          tED->GetEntry(pick);

                          const string outPng =
                            JoinPath(outDir,
                                     "eventDisplay_MissA_" + rKey +
                                     "_run" + std::to_string(b_run) +
                                     "_evt" + std::to_string(b_evt) + ".png");

                          cout << ANSI_BOLD_CYN << "  [EventDisplay][DO]" << ANSI_RESET
                               << " MissA"
                               << " pick=" << pick
                               << "  run=" << b_run
                               << "  evt=" << b_evt
                               << "  vz=" << b_vz
                               << "  ptGammaTruth=" << b_ptGammaTruth
                               << "  -> " << outPng
                               << "\n";

                          const float selPhiW     = WrapPhi(b_sel_phi);
                          const float bestPhiW    = WrapPhi(b_best_phi);
                          const float truthPhiW   = WrapPhi(b_truth_phi);

                          const float dPhiSelTruth  = WrapPhi(selPhiW  - truthPhiW);
                          const float dPhiBestTruth = WrapPhi(bestPhiW - truthPhiW);
                          const float dPhiSelBest   = WrapPhi(selPhiW  - bestPhiW);

                          cout << "    [ED][KIN] sel(phi,eta)=(" << b_sel_phi  << "," << b_sel_eta  << ")  wrappedPhi=" << selPhiW  << "\n";
                          cout << "    [ED][KIN] best(phi,eta)=(" << b_best_phi << "," << b_best_eta << ")  wrappedPhi=" << bestPhiW << "\n";
                          cout << "    [ED][KIN] truth(phi,eta)=(" << b_truth_phi<< "," << b_truth_eta<< ")  wrappedPhi=" << truthPhiW << "\n";
                          cout << "    [ED][KIN] dphi(sel-truth)=" << dPhiSelTruth  << " |dphi|=" << std::fabs(dPhiSelTruth)
                               << "  dphi(best-truth)=" << dPhiBestTruth << " |dphi|=" << std::fabs(dPhiBestTruth)
                               << "  dphi(sel-best)="  << dPhiSelBest   << " |dphi|=" << std::fabs(dPhiSelBest)
                               << "\n";

                          TH2F hSelColz("hED_sel_colz", "", 64, -M_PI, M_PI, 48, -1.1, 1.1);
                          TH2F hSelLego("hED_sel_lego", "", 64, -M_PI, M_PI, 48, -1.1, 1.1);
                          TH2F hBestColz("hED_best_colz", "", 64, -M_PI, M_PI, 48, -1.1, 1.1);
                          TH2F hBestLego("hED_best_lego", "", 64, -M_PI, M_PI, 48, -1.1, 1.1);

                          DumpTowerTriplet("MissA: sel",
                                           b_sel_phiTower, b_sel_etaTower, b_sel_etTower,
                                           hSelColz.GetXaxis()->GetXmin(), hSelColz.GetXaxis()->GetXmax(),
                                           hSelColz.GetYaxis()->GetXmin(), hSelColz.GetYaxis()->GetXmax());

                          DumpTowerTriplet("MissA: best",
                                           b_best_phiTower, b_best_etaTower, b_best_etTower,
                                           hBestColz.GetXaxis()->GetXmin(), hBestColz.GetXaxis()->GetXmax(),
                                           hBestColz.GetYaxis()->GetXmin(), hBestColz.GetYaxis()->GetXmax());

                          FillJetHist(&hSelColz,  b_sel_phiTower,  b_sel_etaTower,  b_sel_etTower);
                          FillJetHist(&hSelLego,  b_sel_phiTower,  b_sel_etaTower,  b_sel_etTower);
                          FillJetHist(&hBestColz, b_best_phiTower, b_best_etaTower, b_best_etTower);
                          FillJetHist(&hBestLego, b_best_phiTower, b_best_etaTower, b_best_etTower);

                          cout << "    [ED][HIST] hSelColz:  entries=" << hSelColz.GetEntries()
                               << " sumW=" << hSelColz.GetSumOfWeights()
                               << " maxBin=" << hSelColz.GetMaximum()
                               << " integral=" << hSelColz.Integral()
                               << "\n";
                          cout << "    [ED][HIST] hBestColz: entries=" << hBestColz.GetEntries()
                               << " sumW=" << hBestColz.GetSumOfWeights()
                               << " maxBin=" << hBestColz.GetMaximum()
                               << " integral=" << hBestColz.Integral()
                               << "\n";

                          TCanvas c("cED_MissA", "", 1600, 1400);
                          c.Divide(2, 2, 0.0, 0.0);

                          const string header = "EventDisplay MissA  " + rKey;
                          const string sub    = "run " + std::to_string(b_run) +
                                                "  evt " + std::to_string(b_evt) +
                                                "  v_{z}=" + std::to_string((int)std::round(b_vz)) + " cm" +
                                                "  pT_{#gamma}^{truth}=" + std::to_string((int)std::round(b_ptGammaTruth)) + " GeV";

                          DrawPanelSave3DStyle((TPad*)c.cd(1), &hSelColz,  "COLZ",  selPhiW,  b_sel_eta,  truthPhiW, b_truth_eta, header, sub + "  (selected jet)");
                          DrawPanelSave3DStyle((TPad*)c.cd(2), &hSelLego,  "LEGO2", selPhiW,  b_sel_eta,  truthPhiW, b_truth_eta, header, sub + "  (selected jet)");
                          DrawPanelSave3DStyle((TPad*)c.cd(3), &hBestColz, "COLZ",  bestPhiW, b_best_eta, truthPhiW, b_truth_eta, header, sub + "  (truth-matched reco jet)");
                          DrawPanelSave3DStyle((TPad*)c.cd(4), &hBestLego, "LEGO2", bestPhiW, b_best_eta, truthPhiW, b_truth_eta, header, sub + "  (truth-matched reco jet)");

                          c.SaveAs(outPng.c_str());
                          PrintED("  [EventDisplay][WROTE]", outPng, ANSI_BOLD_GRN);

                          // -----------------------------------------------------------------
                          // (ED) Publication-ready 3D-only panel for the truth-matched recoil jet in MissA
                          // -----------------------------------------------------------------
                          {
                              const string outPng3D =
                                JoinPath(dirED_MTJ,
                                         "matchedTruthJets3D_MissA_" + rKey +
                                         "_run" + std::to_string(b_run) +
                                         "_evt" + std::to_string(b_evt) + ".png");

                              TCanvas c3D("cED_MatchedTruthJets3D_MissA", "", 1400, 1000);

                              const string header3D = "Matched Recoil Truth Jets, Photon 10 + 20 GeV Sim";
                              const string sub3D    = rKey + "  run " + std::to_string(b_run) +
                                                      "  evt " + std::to_string(b_evt) +
                                                      "  v_{z}=" + std::to_string((int)std::round(b_vz)) + " cm" +
                                                      "  pT_{#gamma}^{truth}=" + std::to_string((int)std::round(b_ptGammaTruth)) + " GeV" +
                                                      "  (truth-matched reco jet)";

                              DrawPanelSave3DStyle((TPad*)c3D.cd(), &hBestLego, "LEGO2",
                                                   bestPhiW, b_best_eta, truthPhiW, b_truth_eta,
                                                   header3D, sub3D);

                              c3D.SaveAs(outPng3D.c_str());
                              PrintED("  [EventDisplay][WROTE]", outPng3D, ANSI_BOLD_GRN);
                            }
                          };

                          SaveNUMorMissB("NUM",   0, dirED_NUM);
                          SaveMissA(dirED_MA);
                          SaveNUMorMissB("MissB", 2, dirED_MB);

                          // -------------------------------------------------------------------
                          //  Side-by-side 3D-only (LEGO2) comparison panel:
                          //   NUM (selected jet)  vs  MissA (selected jet)
                          // Produces one additional PNG under dirED (not inside NUM/MissA/MissB).
                          // -------------------------------------------------------------------
                          {
                            const Long64_t pickNUM = PickIndex(0);

                            Long64_t pickMA = -1;
                            if (!idxByCat[1].empty())
                            {
                              for (const Long64_t ient : idxByCat[1])
                              {
                                tED->GetEntry(ient);
                                if (!b_best_etTower) continue;
                                if (!b_best_etTower->empty())
                                {
                                  pickMA = ient;
                                  break;
                                }
                              }
                              if (pickMA < 0) pickMA = PickIndex(1);
                            }

                            if (pickNUM >= 0 && pickMA >= 0)
                            {
                              // -----------------------------
                              // Load NUM entry (selected jet)
                              // -----------------------------
                              tED->GetEntry(pickNUM);

                              const int   num_run = b_run;
                              const int   num_evt = b_evt;
                              const float num_vz  = b_vz;
                              const float num_pt  = b_ptGammaTruth;

                              const float num_selPhiW   = WrapPhi(b_sel_phi);
                              const float num_selEta    = b_sel_eta;
                              const float num_truthPhiW = WrapPhi(b_truth_phi);
                              const float num_truthEta  = b_truth_eta;

                              TH2F hNumLego("hED_num_lego", "", 64, -M_PI, M_PI, 48, -1.1, 1.1);
                              FillJetHist(&hNumLego, b_sel_phiTower, b_sel_etaTower, b_sel_etTower);

                              // -----------------------------
                              // Load MissA entry (truth-matched reco jet)
                              // -----------------------------
                              tED->GetEntry(pickMA);

                              const int   ma_run = b_run;
                              const int   ma_evt = b_evt;
                              const float ma_vz  = b_vz;
                              const float ma_pt  = b_ptGammaTruth;

                              const float ma_bestPhiW  = WrapPhi(b_best_phi);
                              const float ma_bestEta   = b_best_eta;
                              const float ma_truthPhiW = WrapPhi(b_truth_phi);
                              const float ma_truthEta  = b_truth_eta;

                              TH2F hMALego("hED_missa_lego", "", 64, -M_PI, M_PI, 48, -1.1, 1.1);
                              FillJetHist(&hMALego, b_best_phiTower, b_best_etaTower, b_best_etTower);

                              const string outPng =
                                JoinPath(dirED,
                                         "eventDisplay_NUM_vs_MissA_3D_" + rKey +
                                         "_NUMrun" + std::to_string(num_run) + "_evt" + std::to_string(num_evt) +
                                         "_MArun"  + std::to_string(ma_run)  + "_evt" + std::to_string(ma_evt) + ".png");

                              cout << ANSI_BOLD_CYN << "  [EventDisplay][DO]" << ANSI_RESET
                                   << " NUM vs MissA (3D-only)"
                                   << "  NUM(run,evt)=(" << num_run << "," << num_evt << ")"
                                   << "  MissA(run,evt)=(" << ma_run << "," << ma_evt << ")"
                                   << "  -> " << outPng
                                   << "\n";

                              TCanvas c2("cED_NUM_vs_MissA_3D", "", 1800, 800);
                              c2.Divide(2, 1, 0.0, 0.0);

                              const string header = "EventDisplay 3D  NUM vs MissA  " + rKey;

                              const string subNUM =
                                "NUM: run " + std::to_string(num_run) +
                                "  evt " + std::to_string(num_evt) +
                                "  v_{z}=" + std::to_string((int)std::round(num_vz)) + " cm" +
                                "  pT_{#gamma}^{truth}=" + std::to_string((int)std::round(num_pt)) + " GeV";

                              const string subMA =
                                "MissA: run " + std::to_string(ma_run) +
                                "  evt " + std::to_string(ma_evt) +
                                "  v_{z}=" + std::to_string((int)std::round(ma_vz)) + " cm" +
                                "  pT_{#gamma}^{truth}=" + std::to_string((int)std::round(ma_pt)) + " GeV";

                              DrawPanelSave3DStyle((TPad*)c2.cd(1), &hNumLego, "LEGO2", num_selPhiW, num_selEta, num_truthPhiW, num_truthEta, header, subNUM);
                              DrawPanelSave3DStyle((TPad*)c2.cd(2), &hMALego,  "LEGO2", ma_bestPhiW, ma_bestEta, ma_truthPhiW, ma_truthEta, header, subMA);

                              c2.SaveAs(outPng.c_str());
                              PrintED("  [EventDisplay][WROTE]", outPng, ANSI_BOLD_GRN);
                            }
                            else
                            {
                              PrintED("  [EventDisplay][SKIP]",
                                      "NUM vs MissA (3D-only) skipped: missing NUM and/or MissA entry for rKey=\"" + rKey + "\"",
                                      ANSI_BOLD_YEL);
                            }
                          }
                      }
                  }


                  auto Save2DColz =
                    [&](TH2* hIn,
                        const string& outPng,
                        const string& mainTitle,
                        const string& xTitle,
                        const string& yTitle,
                        bool drawDiag)->bool
                  {
                    if (!hIn) return false;

                    std::unique_ptr<TH2> h( CloneTH2(hIn, TString::Format("h2c_%s_%s", rKey.c_str(), mainTitle.c_str()).Data()) );
                    if (!h) return false;
                    EnsureSumw2(h.get());

                    TCanvas c(TString::Format("c_diag2D_%s_%s_%s", ds.label.c_str(), rKey.c_str(), mainTitle.c_str()).Data(),
                              "c_diag2D", 950, 780);
                    ApplyCanvasMargins2D(c);

                    h->SetTitle("");
                    h->GetXaxis()->SetTitle(xTitle.c_str());
                    h->GetYaxis()->SetTitle(yTitle.c_str());
                    h->Draw("colz");

                    if (drawDiag)
                    {
                      const double xmin = h->GetXaxis()->GetXmin();
                      const double xmax = h->GetXaxis()->GetXmax();
                      const double ymin = h->GetYaxis()->GetXmin();
                      const double ymax = h->GetYaxis()->GetXmax();
                      const double lo = std::max(xmin, ymin);
                      const double hi = std::min(xmax, ymax);
                      TLine diag(lo, lo, hi, hi);
                      diag.SetLineStyle(2);
                      diag.SetLineWidth(2);
                      diag.Draw("same");
                    }

                    DrawLatexLines(0.14, 0.94, DefaultHeaderLines(ds), 0.028, 0.038);
                    DrawLatexLines(0.14, 0.86,
                      {TString::Format("%s (%s)", mainTitle.c_str(), rKey.c_str()).Data()},
                      0.028, 0.038
                    );

                    SaveCanvas(c, outPng);
                    cout << "  [LeadTruthRecoilMatchDiag] Wrote: " << outPng << "\n";
                    effSummary.push_back("  - " + outPng);
                    return true;
                  };

                  // -------------------------
                  // Load new TH2 diagnostics
                  // -------------------------
                  const string hA1_num   = "h2_leadTruthRecoilMatch_pTrecoJet1_vs_pTtruthLead_num_"   + rKey;
                  const string hA1_mA    = "h2_leadTruthRecoilMatch_pTrecoJet1_vs_pTtruthLead_missA_" + rKey;
                  const string hA1_mA1   = "h2_leadTruthRecoilMatch_pTrecoJet1_vs_pTtruthLead_missA1_" + rKey;
                  const string hA1_mA2   = "h2_leadTruthRecoilMatch_pTrecoJet1_vs_pTtruthLead_missA2_" + rKey;
                  const string hA1_mB    = "h2_leadTruthRecoilMatch_pTrecoJet1_vs_pTtruthLead_missB_" + rKey;

                  const string hA2_num   = "h2_leadTruthRecoilMatch_pTrecoJet1_vs_pTrecoTruthMatch_num_"   + rKey;
                  const string hA2_mA    = "h2_leadTruthRecoilMatch_pTrecoJet1_vs_pTrecoTruthMatch_missA_" + rKey;

                  const string hB3_num   = "h2_leadTruthRecoilMatch_dphiRecoJet1_num_pTgammaTruth_"   + rKey;
                  const string hB3_mA    = "h2_leadTruthRecoilMatch_dphiRecoJet1_missA_pTgammaTruth_" + rKey;
                  const string hB3_mB    = "h2_leadTruthRecoilMatch_dphiRecoJet1_missB_pTgammaTruth_" + rKey;

                  const string hB4_num   = "h2_leadTruthRecoilMatch_dRRecoJet1_vs_truthLead_num_pTgammaTruth_"   + rKey;
                  const string hB4_mA    = "h2_leadTruthRecoilMatch_dRRecoJet1_vs_truthLead_missA_pTgammaTruth_" + rKey;
                  const string hB4_mB    = "h2_leadTruthRecoilMatch_dRRecoJet1_vs_truthLead_missB_pTgammaTruth_" + rKey;

                  const string hC5_num   = "h2_leadTruthRecoilMatch_xJRecoJet1_vs_dphiRecoJet1_num_"   + rKey;
                  const string hC5_mA    = "h2_leadTruthRecoilMatch_xJRecoJet1_vs_dphiRecoJet1_missA_" + rKey;
                  const string hC5_mB    = "h2_leadTruthRecoilMatch_xJRecoJet1_vs_dphiRecoJet1_missB_" + rKey;

                  TH2* HA1n  = GetObj<TH2>(ds, hA1_num, false, false, false);
                  TH2* HA1a  = GetObj<TH2>(ds, hA1_mA,  false, false, false);
                  TH2* HA1a1 = GetObj<TH2>(ds, hA1_mA1, false, false, false);
                  TH2* HA1a2 = GetObj<TH2>(ds, hA1_mA2, false, false, false);
                  TH2* HA1b  = GetObj<TH2>(ds, hA1_mB,  false, false, false);

                  TH2* HA2n = GetObj<TH2>(ds, hA2_num, false, false, false);
                  TH2* HA2a = GetObj<TH2>(ds, hA2_mA,  false, false, false);

                  TH2* HB3n = GetObj<TH2>(ds, hB3_num, false, false, false);
                  TH2* HB3a = GetObj<TH2>(ds, hB3_mA,  false, false, false);
                  TH2* HB3b = GetObj<TH2>(ds, hB3_mB,  false, false, false);

                  TH2* HB4n = GetObj<TH2>(ds, hB4_num, false, false, false);
                  TH2* HB4a = GetObj<TH2>(ds, hB4_mA,  false, false, false);
                  TH2* HB4b = GetObj<TH2>(ds, hB4_mB,  false, false, false);

                  TH2* HC5n = GetObj<TH2>(ds, hC5_num, false, false, false);
                  TH2* HC5a = GetObj<TH2>(ds, hC5_mA,  false, false, false);
                  TH2* HC5b = GetObj<TH2>(ds, hC5_mB,  false, false, false);

                    // -------------------------
                    // (A1) pT(truth lead recoil) vs pT(reco recoilJet1) by class
                    // -------------------------
                    if (HA1n) Save2DColz(HA1n,
                      JoinPath(dirPt, TString::Format("leadTruthRecoilMatch_pTrecoJet1_vs_pTtruthLead_NUM_%s.png", rKey.c_str()).Data()),
                      "pT(recoJet1) vs pT(truth lead recoil) [NUM]",
                      "p_{T}^{truth lead recoil} [GeV]",
                      "p_{T}^{recoilJet_{1},reco} [GeV]",
                      true);

                    if (HA1a) Save2DColz(HA1a,
                      JoinPath(dirPt, TString::Format("leadTruthRecoilMatch_pTrecoJet1_vs_pTtruthLead_MissA_%s.png", rKey.c_str()).Data()),
                      "pT(recoJet1) vs pT(truth lead recoil) [MissA]",
                      "p_{T}^{truth lead recoil} [GeV]",
                      "p_{T}^{recoilJet_{1},reco} [GeV]",
                      true);

                    if (HA1b) Save2DColz(HA1b,
                      JoinPath(dirPt, TString::Format("leadTruthRecoilMatch_pTrecoJet1_vs_pTtruthLead_MissB_%s.png", rKey.c_str()).Data()),
                      "pT(recoJet1) vs pT(truth lead recoil) [MissB]",
                      "p_{T}^{truth lead recoil} [GeV]",
                      "p_{T}^{recoilJet_{1},reco} [GeV]",
                      true);

                    // -------------------------
                    // (A1b) NEW: ProfileX overlay (1D): <pT(recoJet1)> vs pT(truth lead recoil)
                    //   Overlay: NUM (black) vs MissA (red) vs MissB (blue) + y=x line
                    // Output:
                    //   .../LeadTruthRecoilMatch/Diagnostics/Pt/
                    //     leadTruthRecoilMatch_profile_pTrecoJet1_vs_pTtruthLead_NUM_MissA_MissB_<rKey>.png
                    // -------------------------
                    if (HA1n && HA1a && HA1b)
                    {
                        TProfile* pN_tmp = HA1n->ProfileX(TString::Format("p_prof_A1_NUM_tmp_%s", rKey.c_str()).Data());
                        TProfile* pA_tmp = HA1a->ProfileX(TString::Format("p_prof_A1_MissA_tmp_%s", rKey.c_str()).Data());
                        TProfile* pB_tmp = HA1b->ProfileX(TString::Format("p_prof_A1_MissB_tmp_%s", rKey.c_str()).Data());

                        std::unique_ptr<TProfile> pN(pN_tmp ? static_cast<TProfile*>(pN_tmp->Clone(TString::Format("p_prof_A1_NUM_%s", rKey.c_str()).Data())) : nullptr);
                        std::unique_ptr<TProfile> pA(pA_tmp ? static_cast<TProfile*>(pA_tmp->Clone(TString::Format("p_prof_A1_MissA_%s", rKey.c_str()).Data())) : nullptr);
                        std::unique_ptr<TProfile> pB(pB_tmp ? static_cast<TProfile*>(pB_tmp->Clone(TString::Format("p_prof_A1_MissB_%s", rKey.c_str()).Data())) : nullptr);

                        if (pN_tmp) { pN_tmp->SetDirectory(nullptr); delete pN_tmp; }
                        if (pA_tmp) { pA_tmp->SetDirectory(nullptr); delete pA_tmp; }
                        if (pB_tmp) { pB_tmp->SetDirectory(nullptr); delete pB_tmp; }

                        if (pN) pN->SetDirectory(nullptr);
                        if (pA) pA->SetDirectory(nullptr);
                        if (pB) pB->SetDirectory(nullptr);



                      if (pN && pA && pB)
                      {
                        // Autoscale y-axis using all three profiles (only populated bins)
                        double yMin =  1e9;
                        double yMax = -1e9;

                        auto AccumRange = [&](TProfile* p)
                        {
                          if (!p) return;
                          for (int ib = 1; ib <= p->GetNbinsX(); ++ib)
                          {
                            const double y  = p->GetBinContent(ib);
                            const double ey = p->GetBinError(ib);
                            if (y <= 0.0) continue;
                            yMin = std::min(yMin, y - ey);
                            yMax = std::max(yMax, y + ey);
                          }
                        };

                        AccumRange(pN.get());
                        AccumRange(pA.get());
                        AccumRange(pB.get());

                        if (yMin > yMax) { yMin = 0.0; yMax = 60.0; } // fallback
                        const double pad = 2.0;
                        const double yLo = std::max(0.0, yMin - pad);
                        const double yHi = std::min(60.0, yMax + pad);

                        TCanvas cP(TString::Format("c_prof_A1_%s_%s", ds.label.c_str(), rKey.c_str()).Data(),
                                   "c_prof_A1", 900, 720);
                        ApplyCanvasMargins1D(cP);

                          // Style
                          pN->SetTitle("");
                          pN->GetXaxis()->SetTitle("p_{T}^{truth lead recoil} [GeV]");
                          pN->GetYaxis()->SetTitle("< p_{T}^{recoilJet_{1},reco} > [GeV]");
                          pN->GetYaxis()->SetRangeUser(0.0, yHi);

                        pN->SetLineWidth(2);
                        pA->SetLineWidth(2);
                        pB->SetLineWidth(2);

                        pN->SetMarkerStyle(20); pN->SetMarkerSize(1.00);
                        pA->SetMarkerStyle(20); pA->SetMarkerSize(1.00);
                        pB->SetMarkerStyle(20); pB->SetMarkerSize(1.00);

                        pN->SetLineColor(kBlack);  pN->SetMarkerColor(kBlack);
                        pA->SetLineColor(kRed+1);  pA->SetMarkerColor(kRed+1);
                        pB->SetLineColor(kBlue+1); pB->SetMarkerColor(kBlue+1);

                        // Draw
                        pN->Draw("E1");
                        pA->Draw("E1 SAME");
                        pB->Draw("E1 SAME");

                          // ------------------------------------------------------------
                          // NEW: MissB horizontal weighted-mean line + N(<5 GeV) counter
                          //   - yConst = weighted mean of MissB profile points
                          //   - nBelow5 = number of populated MissB bins with <pT> < 5 GeV
                          // ------------------------------------------------------------
                          double sumW  = 0.0;
                          double sumWY = 0.0;
                          int    nBelow5 = 0;
                          int    nPop = 0;

                          for (int ib = 1; ib <= pB->GetNbinsX(); ++ib)
                          {
                            const double y  = pB->GetBinContent(ib);
                            const double ey = pB->GetBinError(ib);
                            if (!(y > 0.0) || !(ey > 0.0)) continue;
                            ++nPop;
                            if (y < 5.0) ++nBelow5;
                            const double w = 1.0 / (ey * ey);
                            sumW  += w;
                            sumWY += w * y;
                          }
                          const double yConst = (sumW > 0.0 ? (sumWY / sumW) : 0.0);

                          // Draw horizontal mean line (blue dashed), only if defined
                          TLine* hConstLine = nullptr;
                          if (yConst > 0.0)
                          {
                            const double xmin = pN->GetXaxis()->GetXmin();
                            const double xmax = pN->GetXaxis()->GetXmax();
                            hConstLine = new TLine(xmin, yConst, xmax, yConst);
                            hConstLine->SetLineStyle(2);
                            hConstLine->SetLineWidth(2);
                            hConstLine->SetLineColor(kBlue+1);
                            hConstLine->Draw("same");
                          }

                          // Legend (keep bottom-right, add mean-line entry)
                          TLegend leg(0.58, 0.23, 0.9, 0.44);
                          leg.SetTextFont(42);
                          leg.SetTextSize(0.028);
                          leg.SetFillStyle(0);
                          leg.SetBorderSize(0);
                          leg.AddEntry(pN.get(), "NUM (truth-matched)", "ep");
                          leg.AddEntry(pA.get(), "MissA (wrong jet)", "ep");
                          leg.AddEntry(pB.get(), "MissB (no reco match)", "ep");
                          if (hConstLine)
                          {
                            leg.AddEntry(hConstLine,
                              TString::Format("MissB flat mean").Data(),
                              "l"
                            );
                          }
                          leg.Draw();

                          // Header + title (as before)
                          DrawLatexLines(0.14, 0.94, DefaultHeaderLines(ds), 0.030, 0.040);
                          DrawLatexLines(0.14, 0.86,
                            {TString::Format("<p_{T}^{recoilJet,reco}> vs p_{T}^{truth lead recoil}  (%s, R=%.1f)", rKey.c_str(), R).Data()},
                            0.030, 0.040
                          );

                          // Default SIM cuts — print UNDER the title block (2 lines)
                          DrawLatexLines(0.14, 0.80,
                            {
                              TString::Format("#Delta#phi(#gamma,jet) > %s", B2BLabel().c_str()).Data(),
                              TString::Format("p_{T}^{jet} > %.0f GeV", static_cast<double>(kJetPtMin)).Data()
                            },
                            0.028, 0.036
                          );

                          // MissB mean — place on the RIGHT, below the cut block (no overlap)
                          if (yConst > 0.0)
                          {
                            TLatex t;
                            t.SetNDC(true);
                            t.SetTextFont(42);
                            t.SetTextAlign(31); // right-aligned
                            t.SetTextSize(0.032);
                            t.DrawLatex(0.39, 0.71, TString::Format("MissB mean = %.2f GeV", yConst).Data());
                          }

                          const string outProf = JoinPath(
                            dirPt,
                            TString::Format("leadTruthRecoilMatch_profile_pTrecoJet1_vs_pTtruthLead_NUM_MissA_MissB_%s.png", rKey.c_str()).Data()
                          );
                          SaveCanvas(cP, outProf);

                          cout << "  [LeadTruthRecoilMatchDiag] Wrote: " << outProf << "\n";
                          effSummary.push_back("  - " + outProf);


                          // -------------------------
                          // (A1c) NEW: ProfileX overlay (1D): <pT(recoJet1)> vs pT(truth lead recoil) within MissA subtypes
                          //   Overlay: MissA (black) vs MissA1 (red) vs MissA2 (blue) + y=x line
                          // Output:
                          //   .../LeadTruthRecoilMatch/Diagnostics/Pt/
                          //     leadTruthRecoilMatch_profile_pTrecoJet1_vs_pTtruthLead_MissA_MissA1_MissA2_<rKey>.png
                          // -------------------------
                          if (HA1a && HA1a1 && HA1a2)
                          {
                              TProfile* pA_tmp  = HA1a ->ProfileX(TString::Format("p_prof_A1_MissA_tmp_%s",  rKey.c_str()).Data());
                              TProfile* pA1_tmp = HA1a1->ProfileX(TString::Format("p_prof_A1_MissA1_tmp_%s", rKey.c_str()).Data());
                              TProfile* pA2_tmp = HA1a2->ProfileX(TString::Format("p_prof_A1_MissA2_tmp_%s", rKey.c_str()).Data());

                              std::unique_ptr<TProfile> pA (pA_tmp  ? static_cast<TProfile*>(pA_tmp ->Clone(TString::Format("p_prof_A1_MissA_%s",  rKey.c_str()).Data())) : nullptr);
                              std::unique_ptr<TProfile> pA1(pA1_tmp ? static_cast<TProfile*>(pA1_tmp->Clone(TString::Format("p_prof_A1_MissA1_%s", rKey.c_str()).Data())) : nullptr);
                              std::unique_ptr<TProfile> pA2(pA2_tmp ? static_cast<TProfile*>(pA2_tmp->Clone(TString::Format("p_prof_A1_MissA2_%s", rKey.c_str()).Data())) : nullptr);

                              if (pA_tmp)  { pA_tmp->SetDirectory(nullptr);  delete pA_tmp; }
                              if (pA1_tmp) { pA1_tmp->SetDirectory(nullptr); delete pA1_tmp; }
                              if (pA2_tmp) { pA2_tmp->SetDirectory(nullptr); delete pA2_tmp; }

                              if (pA)  pA->SetDirectory(nullptr);
                              if (pA1) pA1->SetDirectory(nullptr);
                              if (pA2) pA2->SetDirectory(nullptr);



                            if (pA && pA1 && pA2)
                            {
                              // Autoscale y-axis using all three profiles (only populated bins)
                              double yMin =  1e9;
                              double yMax = -1e9;

                              auto AccumRange = [&](TProfile* p)
                              {
                                if (!p) return;
                                for (int ib = 1; ib <= p->GetNbinsX(); ++ib)
                                {
                                  const double y  = p->GetBinContent(ib);
                                  const double ey = p->GetBinError(ib);
                                  if (y <= 0.0) continue;
                                  yMin = std::min(yMin, y - ey);
                                  yMax = std::max(yMax, y + ey);
                                }
                              };

                              AccumRange(pA.get());
                              AccumRange(pA1.get());
                              AccumRange(pA2.get());

                              if (yMin > yMax) { yMin = 0.0; yMax = 60.0; } // fallback
                              const double pad = 2.0;
                              const double yLo = std::max(0.0, yMin - pad);
                              const double yHi = yMax + pad;

                              TCanvas cP2(TString::Format("c_prof_A1_MissA_Subtypes_%s_%s", ds.label.c_str(), rKey.c_str()).Data(),
                                          "c_prof_A1_MissA_Subtypes", 900, 720);
                              ApplyCanvasMargins1D(cP2);

                              // Axis frame: use pA
                              pA->SetTitle("");
                              pA->GetXaxis()->SetTitle("p_{T}^{truth lead recoil} [GeV]");
                              pA->GetYaxis()->SetTitle("<p_{T}^{recoilJet_{1},reco}> [GeV]");
                              pA->GetYaxis()->SetRangeUser(yLo, yHi);

                              // Draw axis only
                              pA->SetLineWidth(0);
                              pA->SetMarkerSize(0);
                              pA->Draw("AXIS");

                              // Style
                              pA->SetLineColor(kBlack);
                              pA->SetMarkerColor(kBlack);
                              pA->SetMarkerStyle(20);
                              pA->SetMarkerSize(1.05);

                              pA1->SetLineColor(kRed+1);
                              pA1->SetMarkerColor(kRed+1);
                              pA1->SetMarkerStyle(20);
                              pA1->SetMarkerSize(1.05);

                              pA2->SetLineColor(kBlue+1);
                              pA2->SetMarkerColor(kBlue+1);
                              pA2->SetMarkerStyle(20);
                              pA2->SetMarkerSize(1.05);

                              pA ->Draw("PE SAME");
                              pA1->Draw("PE SAME");
                              pA2->Draw("PE SAME");

                              // y=x reference
                              const double xMin = pA->GetXaxis()->GetXmin();
                              const double xMax = pA->GetXaxis()->GetXmax();
                              const double lo = std::max(xMin, yLo);
                              const double hi = std::min(xMax, yHi);
                              if (hi > lo)
                              {
                                TLine diag(lo, lo, hi, hi);
                                diag.SetLineStyle(2);
                                diag.SetLineWidth(2);
                                diag.Draw("same");
                              }

                              // Header + title
                              DrawLatexLines(0.14, 0.94, DefaultHeaderLines(ds), 0.030, 0.040);
                              DrawLatexLines(0.14, 0.86,
                                {TString::Format("<pT(recoJet1)> vs pT(truth lead recoil) (MissA subtypes) (%s)", rKey.c_str()).Data()},
                                0.030, 0.040
                              );

                              // Legend
                              TLegend leg2(0.38, 0.76, 0.70, 0.90);
                              leg2.SetTextFont(42);
                              leg2.SetTextSize(0.028);
                              leg2.SetFillStyle(0);
                              leg2.SetBorderSize(0);
                              leg2.AddEntry(pA.get(),  "MissA",  "ep");
                              leg2.AddEntry(pA1.get(), "MissA1", "ep");
                              leg2.AddEntry(pA2.get(), "MissA2", "ep");
                              leg2.Draw();

                              const string outProf2 = JoinPath(
                                dirPt,
                                TString::Format("leadTruthRecoilMatch_profile_pTrecoJet1_vs_pTtruthLead_MissA_MissA1_MissA2_%s.png", rKey.c_str()).Data()
                              );
                              SaveCanvas(cP2, outProf2);

                              cout << "  [LeadTruthRecoilMatchDiag] Wrote: " << outProf2 << "\n";
                              effSummary.push_back("  - " + outProf2);
                            }
                          }

                          // ------------------------------------------------------------------
                          // Terminal diagnostics for MissB profile (blue)
                          //   - Print per truth-pT bin: <pT(reco)> ± err
                          //   - Check whether values stay above 5 GeV
                          //   - Compute weighted mean (constant) and chi2/ndf vs flat hypothesis
                          // ------------------------------------------------------------------
                          {
                            // ANSI helpers (use what your file already uses elsewhere)
                            const std::string B = ANSI_BOLD_CYN;
                            const std::string Y = ANSI_BOLD_YEL;
                            const std::string Rr = ANSI_BOLD_RED;
                            const std::string G = ANSI_BOLD_GRN;
                            const std::string N = ANSI_RESET;

                            auto F2 = [](double x) {
                              std::ostringstream os; os.setf(std::ios::fixed); os << std::setprecision(2) << x; return os.str();
                            };
                            auto F3 = [](double x) {
                              std::ostringstream os; os.setf(std::ios::fixed); os << std::setprecision(3) << x; return os.str();
                            };

                            cout << B << "\n[MissB pT floor diagnostic] " << N
                                 << "(profile: <pT^{reco}> vs pT^{truth lead recoil})  "
                                 << "rKey=" << rKey << "  R=" << std::fixed << std::setprecision(1) << R << "\n";

                            cout << "  Expectation: selected reco jet pT >= 5 GeV (analysis jet cut)\n";
                            cout << "  Values below 5 GeV would indicate a selection/fill mismatch.\n\n";

                            // Table header
                            cout << "  "
                                 << std::left << std::setw(14) << "truth pT bin"
                                 << std::right << std::setw(14) << "<pT> [GeV]"
                                 << std::right << std::setw(14) << "err [GeV]"
                                 << std::right << std::setw(10) << "N"
                                 << "\n";
                            cout << "  " << std::string(52, '-') << "\n";

                            double minMean = 1e9;
                            int    minBin  = -1;

                            // Weighted-mean (constant) computation
                            double sumW  = 0.0;
                            double sumWY = 0.0;
                            double chi2  = 0.0;
                            int    nUsed = 0;

                            // First pass: accumulate weighted mean
                            for (int ib = 1; ib <= pB->GetNbinsX(); ++ib)
                            {
                              const double y  = pB->GetBinContent(ib);
                              const double ey = pB->GetBinError(ib);
                              const double Nbin = pB->GetBinEntries(ib);

                              if (!(y > 0.0) || !(ey > 0.0)) continue; // skip empty bins
                              const double w = 1.0 / (ey * ey);
                              sumW  += w;
                              sumWY += w * y;
                              ++nUsed;
                            }

                            const double yConst = (sumW > 0.0 ? (sumWY / sumW) : 0.0);

                            // Second pass: print rows + chi2
                            for (int ib = 1; ib <= pB->GetNbinsX(); ++ib)
                            {
                              const double y  = pB->GetBinContent(ib);
                              const double ey = pB->GetBinError(ib);
                              const double Nbin = pB->GetBinEntries(ib);

                              // Bin label from x-axis edges
                              const double xlo = pB->GetXaxis()->GetBinLowEdge(ib);
                              const double xhi = pB->GetXaxis()->GetBinUpEdge(ib);

                              if (!(y > 0.0) || !(ey > 0.0))
                              {
                                // still print bin with N=0 for clarity
                                std::ostringstream bl;
                                bl << (int)std::lround(xlo) << "-" << (int)std::lround(xhi);

                                cout << "  "
                                       << std::left  << std::setw(14) << bl.str()
                                       << std::right << std::setw(14) << "--"
                                       << std::right << std::setw(14) << "--"
                                       << std::right << std::setw(10) << "0"
                                       << "\n";

                                continue;
                              }

                              // track min
                              if (y < minMean)
                              {
                                minMean = y;
                                minBin  = ib;
                              }

                              // chi2 contribution
                              const double pull = (y - yConst) / ey;
                              chi2 += pull * pull;
                              ++nUsed;

                              const std::string flag =
                                (y < 5.0 ? (Rr + "  <5!" + N) :
                                 (y < 6.0 ? (Y + "  ~thr" + N) :
                                  (G + "   OK" + N)));

                                std::ostringstream bl;
                                bl << (int)std::lround(xlo) << "-" << (int)std::lround(xhi);

                                cout << "  "
                                     << std::left  << std::setw(14) << bl.str()
                                     << std::right << std::setw(14) << F2(y)
                                     << std::right << std::setw(14) << F2(ey)
                                     << std::right << std::setw(10) << (int)std::lround(Nbin)
                                     << flag
                                     << "\n";

                            }

                            // Summary line
                            const int ndf = std::max(0, nUsed - 1);
                            const double chi2ndf = (ndf > 0 ? chi2 / ndf : 0.0);

                            cout << "  " << std::string(52, '-') << "\n";
                            if (minBin > 0)
                            {
                              const double xlo = pB->GetXaxis()->GetBinLowEdge(minBin);
                              const double xhi = pB->GetXaxis()->GetBinUpEdge(minBin);
                                cout << "  Min <pT> = " << F2(minMean) << " GeV"
                                     << " in truth bin " << (int)std::lround(xlo) << "-" << (int)std::lround(xhi) << " GeV\n";

                            }
                            cout << "  Flat-mean (weighted) <pT>_const = " << F2(yConst) << " GeV"
                                 << "   chi2/ndf = " << F2(chi2ndf) << "\n";

                            if (minMean < 5.0)
                            {
                              cout << Rr << "  WARNING: MissB profile dips below 5 GeV -> investigate selection/fill consistency." << N << "\n";
                            }
                            else
                            {
                              cout << G << "  OK: MissB profile stays >= 5 GeV (consistent with analysis jet pT cut)." << N << "\n";
                            }
                            cout << "\n";
                          }

                      }
                    }

                  // -------------------------
                  // (A2) pT(reco match to truth lead) vs pT(selected recoJet1) (NUM & MissA)
                  // -------------------------
                  if (HA2n) Save2DColz(HA2n,
                    JoinPath(dirPt, TString::Format("leadTruthRecoilMatch_pTrecoJet1_vs_pTrecoTruthMatch_NUM_%s.png", rKey.c_str()).Data()),
                    "pT(recoJet1) vs pT(reco truth-match) [NUM]",
                    "p_{T}^{reco truth-match} [GeV]",
                    "p_{T}^{recoilJet_{1},reco} [GeV]",
                    true);

                  if (HA2a) Save2DColz(HA2a,
                    JoinPath(dirPt, TString::Format("leadTruthRecoilMatch_pTrecoJet1_vs_pTrecoTruthMatch_MissA_%s.png", rKey.c_str()).Data()),
                    "pT(recoJet1) vs pT(reco truth-match) [MissA]",
                    "p_{T}^{reco truth-match} [GeV]",
                    "p_{T}^{recoilJet_{1},reco} [GeV]",
                    true);

                  // -------------------------
                  // (B4) ΔR(recoJet1, truth lead recoil) vs pTγ,truth by class
                  // -------------------------
                  if (HB4n) Save2DColz(HB4n,
                    JoinPath(dirAng, TString::Format("leadTruthRecoilMatch_dRRecoJet1_vs_truthLead_NUM_pTgammaTruth_%s.png", rKey.c_str()).Data()),
                    "#DeltaR(recoJet1, truth lead recoil) vs p_{T}^{#gamma,truth} [NUM]",
                    "p_{T}^{#gamma,truth} [GeV]",
                    "#DeltaR(recoilJet_{1}^{reco}, truth lead recoil)",
                    false);

                  if (HB4a) Save2DColz(HB4a,
                    JoinPath(dirAng, TString::Format("leadTruthRecoilMatch_dRRecoJet1_vs_truthLead_MissA_pTgammaTruth_%s.png", rKey.c_str()).Data()),
                    "#DeltaR(recoJet1, truth lead recoil) vs p_{T}^{#gamma,truth} [MissA]",
                    "p_{T}^{#gamma,truth} [GeV]",
                    "#DeltaR(recoilJet_{1}^{reco}, truth lead recoil)",
                    false);

                  if (HB4b) Save2DColz(HB4b,
                    JoinPath(dirAng, TString::Format("leadTruthRecoilMatch_dRRecoJet1_vs_truthLead_MissB_pTgammaTruth_%s.png", rKey.c_str()).Data()),
                    "#DeltaR(recoJet1, truth lead recoil) vs p_{T}^{#gamma,truth} [MissB]",
                    "p_{T}^{#gamma,truth} [GeV]",
                    "#DeltaR(recoilJet_{1}^{reco}, truth lead recoil)",
                    false);

                    // -------------------------
                    // (B3) NEW: 1D overlay — back-to-back tight fraction vs pT^{gamma,truth}
                    //   f_BB(C)(i) = N(|dphi| > dphi0) / N(all), computed per pT^{gamma,truth} bin
                    //   C ∈ {NUM, MissA, MissB}
                    //   Binomial errors via TH1::Divide(...,"B") (same logic as MissA/DEN).
                    //
                    // Inputs:
                    //   h2_leadTruthRecoilMatch_dphiRecoJet1_{num,missA,missB}_pTgammaTruth_<rKey>
                    //
                    // Output:
                    //   <...>/LeadTruthRecoilMatch/Diagnostics/Angle/
                    //     leadTruthRecoilMatch_fBackToBackGT2p8_vs_pTgammaTruth_overlay_NUM_MissA_MissB_<rKey>.png
                    // -------------------------
                    if (HB3n && HB3a && HB3b)
                    {
                      const double dphi0 = 2.8; // very back-to-back

                      auto BuildBinomialFracHist =
                        [&](TH2* h2, const std::string& baseName) -> std::unique_ptr<TH1>
                      {
                        if (!h2) return nullptr;

                        const TAxis* ax = h2->GetXaxis();
                        const int nbx = ax->GetNbins();

                        // Build TH1 with same x-binning as the TH2 x-axis
                        std::unique_ptr<TH1> hAll;
                        std::unique_ptr<TH1> hBB;

                        if (ax->GetXbins() && ax->GetXbins()->GetSize() > 0)
                        {
                          const double* edges = ax->GetXbins()->GetArray();
                          hAll.reset(new TH1D((baseName + "_all").c_str(), "", nbx, edges));
                          hBB .reset(new TH1D((baseName + "_bb").c_str(),  "", nbx, edges));
                        }
                        else
                        {
                          hAll.reset(new TH1D((baseName + "_all").c_str(), "", nbx, ax->GetXmin(), ax->GetXmax()));
                          hBB .reset(new TH1D((baseName + "_bb").c_str(),  "", nbx, ax->GetXmin(), ax->GetXmax()));
                        }

                        if (!hAll || !hBB) return nullptr;

                        hAll->Sumw2();
                        hBB->Sumw2();

                        const int nby = h2->GetYaxis()->GetNbins();
                        int yCut = h2->GetYaxis()->FindBin(dphi0);
                        if (yCut < 1) yCut = 1;
                        if (yCut > nby) yCut = nby;

                          for (int ix = 1; ix <= nbx; ++ix)
                          {
                            double nAll  = 0.0;
                            double e2All = 0.0;

                            double nBB   = 0.0;
                            double e2BB  = 0.0;

                            // Sum over Y bins manually so we also propagate Sumw2
                            for (int iy = 1; iy <= nby; ++iy)
                            {
                              const double v  = h2->GetBinContent(ix, iy);
                              const double ev = h2->GetBinError(ix, iy);   // sqrt(sumw2) if Sumw2 is on
                              nAll  += v;
                              e2All += ev * ev;

                              if (iy >= yCut)
                              {
                                nBB  += v;
                                e2BB += ev * ev;
                              }
                            }

                            hAll->SetBinContent(ix, nAll);
                            hAll->SetBinError(ix, std::sqrt(e2All));

                            hBB->SetBinContent(ix, nBB);
                            hBB->SetBinError(ix, std::sqrt(e2BB));
                          }

                          std::unique_ptr<TH1> hFrac( (TH1*)hBB->Clone((baseName + "_frac").c_str()) );
                          if (!hFrac) return nullptr;
                          hFrac->Sumw2();

                          // Ratio with proper error propagation (works for weighted and unweighted)
                          // (If you *only* want strict binomial for unweighted integer counts, keep "B" for that case.)
                          hFrac->Divide(hBB.get(), hAll.get(), 1.0, 1.0);
                          hFrac->SetDirectory(nullptr);

                          return hFrac;
                      };

                      // Build fraction histograms (NUM/MissA/MissB)
                      auto hF_N = BuildBinomialFracHist(HB3n, "h_fBB_NUM_"   + rKey);
                      auto hF_A = BuildBinomialFracHist(HB3a, "h_fBB_MissA_" + rKey);
                      auto hF_B = BuildBinomialFracHist(HB3b, "h_fBB_MissB_" + rKey);

                      if (hF_N && hF_A && hF_B)
                      {
                          // Convert to TGraphErrors with xerr=0 (vertical-only stat errors)
                          // JES3-only convention: do not display truth pT < 15 GeV
                          const double minTruthPtJES3_plot = 15.0;

                          auto MakeGraphNoXerr = [&](TH1* h, int color, int mStyle) -> std::unique_ptr<TGraphErrors>
                          {
                            std::unique_ptr<TGraphErrors> g(new TGraphErrors());
                            int ip = 0;
                            for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
                            {
                              if (h->GetXaxis()->GetBinLowEdge(ib) < minTruthPtJES3_plot) continue;

                              const double y  = h->GetBinContent(ib);
                              const double ey = h->GetBinError(ib);
                              const double x  = h->GetXaxis()->GetBinCenter(ib);

                              if (ey <= 0.0 && y <= 0.0) continue;

                              g->SetPoint(ip, x, y);
                              g->SetPointError(ip, 0.0, ey);
                              ++ip;
                            }
                            g->SetLineWidth(2);
                            g->SetLineColor(color);
                            g->SetMarkerStyle(mStyle);
                            g->SetMarkerSize(1.00);
                            g->SetMarkerColor(color);
                            return g;
                          };

                        auto gN = MakeGraphNoXerr(hF_N.get(), kBlack, 20);
                        auto gA = MakeGraphNoXerr(hF_A.get(), kRed+1, 20);
                        auto gB = MakeGraphNoXerr(hF_B.get(), kBlue+1, 20);

                        TCanvas cF(TString::Format("c_fBB_%s_%s", ds.label.c_str(), rKey.c_str()).Data(),
                                   "c_fBB", 900, 720);
                        ApplyCanvasMargins1D(cF);

                        // Use NUM hist as axis frame
                        hF_N->SetTitle("");
                        hF_N->GetXaxis()->SetTitle("p_{T}^{#gamma,truth} [GeV]");
                        hF_N->GetYaxis()->SetTitle(TString::Format("f_{BB}(|#Delta#phi|>%.1f) = N_{BB}/N_{all}", dphi0).Data());
                        hF_N->GetYaxis()->SetRangeUser(0.0, 1.05);
                        hF_N->GetXaxis()->SetRangeUser(minTruthPtJES3_plot, hF_N->GetXaxis()->GetXmax());

                        hF_N->SetLineWidth(0);
                        hF_N->SetMarkerSize(0);
                        hF_N->Draw("AXIS");

                        if (gN) gN->Draw("PE SAME");
                        if (gA) gA->Draw("PE SAME");
                        if (gB) gB->Draw("PE SAME");

                          // Legend (bottom-right)
                          TLegend leg(0.58, 0.18, 0.90, 0.34);
                          leg.SetTextFont(42);
                          leg.SetTextSize(0.030);
                          leg.SetFillStyle(0);
                          leg.SetBorderSize(0);
                          if (gN) leg.AddEntry(gN.get(), "NUM (truth-matched)", "ep");
                          if (gA) leg.AddEntry(gA.get(), "MissA (wrong jet)", "ep");
                          if (gB) leg.AddEntry(gB.get(), "MissB (no reco match)", "ep");
                          leg.Draw();


                        DrawLatexLines(0.14, 0.94, DefaultHeaderLines(ds), 0.030, 0.040);
                        DrawLatexLines(0.14, 0.86,
                          {TString::Format("Back-to-back tight fraction vs p_{T}^{#gamma,truth}  (%s, R=%.1f)", rKey.c_str(), R).Data()},
                          0.030, 0.040
                        );

                        const string outF = JoinPath(
                          dirAng,
                          TString::Format("leadTruthRecoilMatch_fBackToBackGT2p8_vs_pTgammaTruth_overlay_NUM_MissA_MissB_%s.png", rKey.c_str()).Data()
                        );
                        SaveCanvas(cF, outF);

                        cout << "  [LeadTruthRecoilMatchDiag] Wrote: " << outF << "\n";
                        effSummary.push_back("  - " + outF);
                      }

                      // -------------------------
                      // (B3) Existing: 3x3 table overlay of dphi(recoJet1) shapes (NUM vs MissA vs MissB)
                      // -------------------------
                      auto NormalizeVisible = [](TH1* h)
                      {
                        if (!h) return;
                        const int nb = h->GetNbinsX();
                        const double integral = h->Integral(1, nb);
                        if (integral > 0.0) h->Scale(1.0 / integral);
                      };

                      // Prefer the canonical bins (kNPtBins) and skip the 2 truth underflow bins if present.
                      // With your updated truth edges {5,10,15,17,19,21,23,26,35,40}, the canonical 6 bins are:
                      //   15-17, 17-19, 19-21, 21-23, 23-26, 26-35  -> startBin=3.
                      const int nXBins = HB3n->GetXaxis()->GetNbins();
                      int startBin = 1;
                      int nUse = std::min(kNPtBins, nXBins);
                      if (nXBins >= (kNPtBins + 2)) { startBin = 3; nUse = kNPtBins; }

                      TCanvas cTbl(
                          TString::Format("c_tbl_dphiRecoJet1_byClass_%s_%s", ds.label.c_str(), rKey.c_str()).Data(),
                          "c_tbl_dphiRecoJet1_byClass", 1500, 1200
                      );
                      cTbl.Divide(3,2, 0.001, 0.001);

                      std::vector<TObject*> keep;
                      keep.reserve(kNPtBins * 6);

                      for (int i = 0; i < kNPtBins; ++i)
                      {
                        cTbl.cd(i+1);
                        gPad->SetLeftMargin(0.14);
                        gPad->SetRightMargin(0.05);
                        gPad->SetBottomMargin(0.14);
                        gPad->SetTopMargin(0.10);
                        gPad->SetTicks(1,1);

                        if (i >= nUse)
                        {
                          TLatex t;
                          t.SetNDC(true);
                          t.SetTextFont(42);
                          t.SetTextSize(0.06);
                          t.DrawLatex(0.20, 0.55, "EMPTY");
                          continue;
                        }

                        const int xbin = startBin + i;

                        TH1D* pN = HB3n->ProjectionY(
                          TString::Format("p_dphiRecoJet1_NUM_%s_%s_%d", ds.label.c_str(), rKey.c_str(), i).Data(),
                          xbin, xbin
                        );
                        TH1D* pA = HB3a->ProjectionY(
                          TString::Format("p_dphiRecoJet1_MissA_%s_%s_%d", ds.label.c_str(), rKey.c_str(), i).Data(),
                          xbin, xbin
                        );
                        TH1D* pB = HB3b->ProjectionY(
                          TString::Format("p_dphiRecoJet1_MissB_%s_%s_%d", ds.label.c_str(), rKey.c_str(), i).Data(),
                          xbin, xbin
                        );

                        if (!pN || !pA || !pB || pN->GetEntries() <= 0.0 || pA->GetEntries() <= 0.0 || pB->GetEntries() <= 0.0)
                        {
                          if (pN) delete pN;
                          if (pA) delete pA;
                          if (pB) delete pB;
                          TLatex t;
                          t.SetNDC(true);
                          t.SetTextFont(42);
                          t.SetTextSize(0.06);
                          t.DrawLatex(0.15, 0.55, "MISSING");
                          continue;
                        }

                        pN->SetDirectory(nullptr);
                        pA->SetDirectory(nullptr);
                        pB->SetDirectory(nullptr);

                        NormalizeVisible(pN);
                        NormalizeVisible(pA);
                        NormalizeVisible(pB);

                        pN->SetLineWidth(2);
                        pA->SetLineWidth(2);
                        pB->SetLineWidth(2);

                        pN->SetLineColor(kBlack);
                        pA->SetLineColor(kRed+1);
                        pB->SetLineColor(kBlue+1);

                        const double ymax = std::max(pN->GetMaximum(), std::max(pA->GetMaximum(), pB->GetMaximum()));
                        pN->SetMaximum(ymax * 1.28);

                        pN->SetTitle("");
                        pN->GetXaxis()->SetTitle("|#Delta#phi(#gamma^{truth}, recoilJet_{1}^{reco})| [rad]");
                        pN->GetYaxis()->SetTitle("A.U.");

                        pN->Draw("hist");
                        pA->Draw("hist same");
                        pB->Draw("hist same");

                        // Reference line at π/2
                        {
                          TLine* l = new TLine(TMath::Pi()/2.0, 0.0, TMath::Pi()/2.0, pN->GetMaximum());
                          l->SetLineStyle(2);
                          l->SetLineWidth(2);
                          l->SetLineColor(kGray+2);
                          l->Draw("same");
                          keep.push_back(l);
                        }

                        TLegend* leg = new TLegend(0.48, 0.66, 0.93, 0.90);
                        leg->SetTextFont(42);
                        leg->SetTextSize(0.028);
                        leg->SetFillStyle(0);
                        leg->SetBorderSize(0);
                        leg->AddEntry(pN, "NUM (truth-matched)", "l");
                        leg->AddEntry(pA, "MissA (wrong jet)", "l");
                        leg->AddEntry(pB, "MissB (no reco match)", "l");
                        leg->Draw();

                        const string ptLab = AxisBinLabel(HB3n->GetXaxis(), xbin, "GeV", 0);
                        DrawLatexLines(
                          0.16, 0.90,
                          {
                            TString::Format("p_{T}^{#gamma,truth}: %s", ptLab.c_str()).Data(),
                            "Overlay (shape)"
                          },
                          0.040, 0.050
                        );

                        keep.push_back(pN);
                        keep.push_back(pA);
                        keep.push_back(pB);
                        keep.push_back(leg);
                      }

                      const string outTbl = JoinPath(
                        dirAng,
                        TString::Format("table3x3_overlay_dphiRecoJet1_byClass_shape_%s.png", rKey.c_str()).Data()
                      );
                      SaveCanvas(cTbl, outTbl);
                      cout << "  [LeadTruthRecoilMatchDiag] Wrote: " << outTbl << "\n";
                      effSummary.push_back("  - " + outTbl);

                      for (auto* o : keep) delete o;
                    }
                    else
                    {
                      cout << ANSI_BOLD_YEL << "  [LeadTruthRecoilMatchDiag] Missing one or more dphi-by-class TH2; skipping 3x3 overlay.\n" << ANSI_RESET;
                    }

                  // -------------------------
                  // (C5) xJ(recoJet1) vs dphi(recoJet1) by class (2D)
                  // -------------------------
                  if (HC5n) Save2DColz(HC5n,
                    JoinPath(dirXJ2D, TString::Format("leadTruthRecoilMatch_xJRecoJet1_vs_dphiRecoJet1_NUM_%s.png", rKey.c_str()).Data()),
                    "xJ(recoJet1) vs dphi(recoJet1) [NUM]",
                    "|#Delta#phi(#gamma^{truth}, recoilJet_{1}^{reco})| [rad]",
                    "x_{J}^{reco} = p_{T}^{recoilJet_{1}} / p_{T}^{#gamma}",
                    false);

                  if (HC5a) Save2DColz(HC5a,
                    JoinPath(dirXJ2D, TString::Format("leadTruthRecoilMatch_xJRecoJet1_vs_dphiRecoJet1_MissA_%s.png", rKey.c_str()).Data()),
                    "xJ(recoJet1) vs dphi(recoJet1) [MissA]",
                    "|#Delta#phi(#gamma^{truth}, recoilJet_{1}^{reco})| [rad]",
                    "x_{J}^{reco} = p_{T}^{recoilJet_{1}} / p_{T}^{#gamma}",
                    false);

                  if (HC5b) Save2DColz(HC5b,
                    JoinPath(dirXJ2D, TString::Format("leadTruthRecoilMatch_xJRecoJet1_vs_dphiRecoJet1_MissB_%s.png", rKey.c_str()).Data()),
                    "xJ(recoJet1) vs dphi(recoJet1) [MissB]",
                    "|#Delta#phi(#gamma^{truth}, recoilJet_{1}^{reco})| [rad]",
                    "x_{J}^{reco} = p_{T}^{recoilJet_{1}} / p_{T}^{#gamma}",
                    false);

                  effSummary.push_back("");
                }

                effSummary.push_back("");

                delete denC;
                delete numC;

            }
          }

            // ---------------------------------------------------------------------------
            // (2) Lead jet response QA (truth<->reco)
            //
            // Inputs (filled online in RecoilJets.cc):
            //   h_leadRecoilJetMatch_dR_<rKey>
            //   h2_leadRecoilJetPtResponse_pTtruth_ratio_<rKey>
            //   h2_leadRecoilJet_pTtruth_pTreco_<rKey>
            // ---------------------------------------------------------------------------
          {
            const string hDRName      = "h_leadRecoilJetMatch_dR_" + rKey;
            const string hRespName    = "h2_leadRecoilJetPtResponse_pTtruth_ratio_" + rKey;
            const string hTruthRecoNm = "h2_leadRecoilJet_pTtruth_pTreco_" + rKey;

            TH1* hDR   = GetObj<TH1>(ds, hDRName,      false, false, false);
            TH2* hResp = GetObj<TH2>(ds, hRespName,    false, false, false);
            TH2* hTR   = GetObj<TH2>(ds, hTruthRecoNm, false, false, false);


            effSummary.push_back("LeadJetResponse:");
            effSummary.push_back("  Output dir: " + D.dirXJProjEffLeadJetResponse);

            // (2a) #DeltaR match
            if (hDR)
            {
              TCanvas c(TString::Format("c_leadJetDR_%s_%s", ds.label.c_str(), rKey.c_str()).Data(),
                        "c_leadJetDR", 900, 700);
              ApplyCanvasMargins1D(c);

              hDR->SetTitle("");
              hDR->GetXaxis()->SetTitle("#DeltaR(lead recoil jet^{reco}, lead recoil jet^{truth})");
              hDR->GetYaxis()->SetTitle((ds.isSim && IsWeightedSIMSelected()) ? "Counts / pb^{-1}" : "Counts");
              gPad->SetLogy(true);

              hDR->SetLineWidth(2);
              hDR->Draw("hist");

              DrawLatexLines(0.14, 0.94, DefaultHeaderLines(ds), 0.030, 0.040);
              DrawLatexLines(0.14, 0.86, {TString::Format("Lead jet match #DeltaR (%s)", rKey.c_str()).Data()}, 0.030, 0.040);

              const string outPng = JoinPath(
                D.dirXJProjEffLeadJetResponse,
                TString::Format("leadJet_match_dR_%s.png", rKey.c_str()).Data()
              );
              SaveCanvas(c, outPng);

              cout << "  [LeadJetResponse] Wrote: " << outPng << "\n";
              effSummary.push_back("  - " + outPng);
            }
            else
            {
              cout << ANSI_BOLD_YEL << "  [LeadJetResponse] Missing " << hDRName << ANSI_RESET << "\n";
              effSummary.push_back("  - MISSING: " + hDRName);
            }

            // (2b) pT response: pTtruth vs (pTreco/pTtruth)
            if (hResp)
            {
              // 2D view
              {
                TCanvas c(TString::Format("c_leadJetResp2D_%s_%s", ds.label.c_str(), rKey.c_str()).Data(),
                          "c_leadJetResp2D", 900, 760);
                ApplyCanvasMargins2D(c);

                hResp->SetTitle("");
                hResp->GetXaxis()->SetTitle("p_{T}^{jet,truth} [GeV]");
                hResp->GetYaxis()->SetTitle("p_{T}^{jet,reco} / p_{T}^{jet,truth}");
                hResp->Draw("colz");

                DrawLatexLines(0.14, 0.94, DefaultHeaderLines(ds), 0.028, 0.038);
                DrawLatexLines(0.14, 0.86, {TString::Format("Lead jet pT response 2D (%s)", rKey.c_str()).Data()}, 0.028, 0.038);

                const string outPng = JoinPath(
                  D.dirXJProjEffLeadJetResponse,
                  TString::Format("leadJet_pTtruth_vs_ratio_%s.png", rKey.c_str()).Data()
                );
                SaveCanvas(c, outPng);

                cout << "  [LeadJetResponse] Wrote: " << outPng << "\n";
                effSummary.push_back("  - " + outPng);
              }

              // ProfileX view
              {
                  TProfile* p_tmp = hResp->ProfileX(
                    TString::Format("p_leadJetResp_tmp_%s", rKey.c_str()).Data()
                  );
                  std::unique_ptr<TProfile> p(p_tmp ? static_cast<TProfile*>(p_tmp->Clone(
                    TString::Format("p_leadJetResp_%s", rKey.c_str()).Data()
                  )) : nullptr);
                  if (p_tmp) { p_tmp->SetDirectory(nullptr); delete p_tmp; }
                  if (p) p->SetDirectory(nullptr);
                  if (p)
                  {
                    p->SetDirectory(nullptr);
                    if (gDirectory) gDirectory->Remove(p.get());

                    TCanvas c(TString::Format("c_leadJetRespProf_%s_%s", ds.label.c_str(), rKey.c_str()).Data(),
                              "c_leadJetRespProf", 900, 700);
                    ApplyCanvasMargins1D(c);

                  p->SetTitle("");
                  p->GetXaxis()->SetTitle("p_{T}^{jet,truth} [GeV]");
                  p->GetYaxis()->SetTitle("< p_{T}^{jet,reco} / p_{T}^{jet,truth} >");
                  p->GetYaxis()->SetRangeUser(0.0, 2.0);

                  p->SetLineWidth(2);
                  p->SetMarkerStyle(20);
                  p->SetMarkerSize(0.95);
                  p->Draw("E1");

                  DrawLatexLines(0.14, 0.94, DefaultHeaderLines(ds), 0.030, 0.040);
                  DrawLatexLines(0.14, 0.86, {TString::Format("Lead jet response profile (%s)", rKey.c_str()).Data()}, 0.030, 0.040);

                  const string outPng = JoinPath(
                    D.dirXJProjEffLeadJetResponse,
                    TString::Format("leadJet_pTresponse_profile_%s.png", rKey.c_str()).Data()
                  );
                  SaveCanvas(c, outPng);

                  cout << "  [LeadJetResponse] Wrote: " << outPng << "\n";
                  effSummary.push_back("  - " + outPng);
                }
              }
            }
            else
            {
              cout << ANSI_BOLD_YEL << "  [LeadJetResponse] Missing " << hRespName << ANSI_RESET << "\n";
              effSummary.push_back("  - MISSING: " + hRespName);
            }

              // (2c) pTtruth vs pTreco (2D)
              if (hTR)
              {
                TCanvas c(TString::Format("c_leadJetTR2D_%s_%s", ds.label.c_str(), rKey.c_str()).Data(),
                          "c_leadJetTR2D", 900, 760);
                ApplyCanvasMargins2D(c);

                hTR->SetTitle("");
                hTR->GetXaxis()->SetTitle("p_{T}^{jet,truth} [GeV]");
                hTR->GetYaxis()->SetTitle("p_{T}^{jet,reco} [GeV]");
                hTR->Draw("colz");

                // y=x reference
                const double xmin = hTR->GetXaxis()->GetXmin();
                const double xmax = hTR->GetXaxis()->GetXmax();
                TLine diag(xmin, xmin, xmax, xmax);
                diag.SetLineStyle(2);
                diag.SetLineWidth(2);
                diag.Draw("same");

                DrawLatexLines(0.14, 0.94, DefaultHeaderLines(ds), 0.028, 0.038);
                DrawLatexLines(0.14, 0.86, {TString::Format("Lead jet pT^{truth} vs pT^{reco} (%s)", rKey.c_str()).Data()}, 0.028, 0.038);

                const string outPng = JoinPath(
                  D.dirXJProjEffLeadJetResponse,
                  TString::Format("leadJet_pTtruth_vs_pTreco_%s.png", rKey.c_str()).Data()
                );
                SaveCanvas(c, outPng);

                cout << "  [LeadJetResponse] Wrote: " << outPng << "\n";
                effSummary.push_back("  - " + outPng);
              }
              else
              {
                cout << ANSI_BOLD_YEL << "  [LeadJetResponse] Missing " << hTruthRecoNm << ANSI_RESET << "\n";
                effSummary.push_back("  - MISSING: " + hTruthRecoNm);
              }

              // (2d) NEW: 3x3 pT^{#gamma}-binned table of the analysis-selected recoilJet1 Δphi distribution
              //   - recoilJet1: |Δphi(γ, recoilJet1)| from h_match_dphi_vs_pTgamma_<rKey>
              //   - NOTE: this histogram is filled only when recoilJet1 exists (analysis recoil requirement applied)
              // Output:
              //   <...>/<rKey>/xJ_fromJES3/Efficiency/LeadJetResponse/
              //     table3x3_dphiRecoilJet1_shape_<rKey>.png
              {
                const string hDphiName = "h_match_dphi_vs_pTgamma_" + rKey;
                TH2* hDphi = GetObj<TH2>(ds, hDphiName, false, false, false);

                if (!hDphi)
                {
                  cout << ANSI_BOLD_YEL << "  [LeadJetResponse] Missing " << hDphiName << ANSI_RESET << "\n";
                  effSummary.push_back("  - MISSING: " + hDphiName);
                }
                else
                {
                  auto NormalizeVisible = [](TH1* h)
                  {
                    if (!h) return;
                    const int nb = h->GetNbinsX();
                    const double integral = h->Integral(1, nb); // visible bins only
                    if (integral > 0.0) h->Scale(1.0 / integral);
                  };

                  const int nPtAxis = hDphi->GetXaxis()->GetNbins();
                  const int nPtUse  = std::min(kNPtBins, std::min(nPtAxis, (int)PtBins().size()));

                  TCanvas c(
                      TString::Format("c_tbl_dphiRecoilJet1_%s_%s", ds.label.c_str(), rKey.c_str()).Data(),
                      "c_tbl_dphiRecoilJet1", 1500, 1200
                  );
                  c.Divide(3,2, 0.001, 0.001);

                  std::vector<TObject*> keep;
                  keep.reserve(2 * nPtUse);

                  for (int i = 0; i < kNPtBins; ++i)
                  {
                    c.cd(i+1);
                    gPad->SetLeftMargin(0.14);
                    gPad->SetRightMargin(0.05);
                    gPad->SetBottomMargin(0.14);
                    gPad->SetTopMargin(0.10);
                    gPad->SetTicks(1,1);

                    if (i >= nPtUse)
                    {
                      TLatex t;
                      t.SetNDC(true);
                      t.SetTextFont(42);
                      t.SetTextSize(0.06);
                      t.DrawLatex(0.20, 0.55, "EMPTY");
                      continue;
                    }

                    const PtBin& pb = PtBins()[i];
                    const int xbin = i + 1;

                    TH1D* p = hDphi->ProjectionY(
                      TString::Format("p_tbl_dphiRecoilJet1_%s_%s_%d", ds.label.c_str(), rKey.c_str(), i).Data(),
                      xbin, xbin
                    );

                    if (!p || p->GetEntries() <= 0.0)
                    {
                      if (p) delete p;
                      TLatex t;
                      t.SetNDC(true);
                      t.SetTextFont(42);
                      t.SetTextSize(0.06);
                      t.DrawLatex(0.15, 0.55, "MISSING");
                      continue;
                    }

                    p->SetDirectory(nullptr);
                    NormalizeVisible(p);

                    p->SetTitle("");
                    p->GetXaxis()->SetTitle("|#Delta#phi(#gamma, recoilJet_{1})| [rad]");
                    p->GetYaxis()->SetTitle("A.U.");

                    p->SetLineWidth(2);
                    p->SetLineColor(kBlack);

                    const double ymax = p->GetMaximum();
                    p->SetMaximum(ymax * 1.28);

                    p->Draw("hist");

                    // Draw a reference line at π/2 (analysis recoil requirement)
                    {
                      TLine* l = new TLine(TMath::Pi()/2.0, 0.0, TMath::Pi()/2.0, p->GetMaximum());
                      l->SetLineStyle(2);
                      l->SetLineWidth(2);
                      l->SetLineColor(kGray+2);
                      l->Draw("same");
                      keep.push_back(l);
                    }

                    // Small legend (even for 1 curve) + indicates π/2 line meaning
                    TLegend* leg = new TLegend(0.52, 0.72, 0.93, 0.90);
                    leg->SetTextFont(42);
                    leg->SetTextSize(0.028);
                    leg->SetFillStyle(0);
                    leg->SetBorderSize(0);
                    leg->AddEntry(p, "recoilJet_{1} (analysis-selected)", "l");
                    leg->AddEntry((TObject*)0, "|#Delta#phi| > #pi/2 required", "");
                    leg->Draw();

                    DrawLatexLines(
                      0.16, 0.90,
                      {
                        TString::Format("p_{T}^{#gamma}: %d-%d GeV", pb.lo, pb.hi).Data(),
                        "recoilJet_{1} #Delta#phi (shape)"
                      },
                      0.040, 0.050
                    );

                    keep.push_back(p);
                    keep.push_back(leg);
                  }

                  const string outPng = JoinPath(
                    D.dirXJProjEffLeadJetResponse,
                    TString::Format("table3x3_dphiRecoilJet1_shape_%s.png", rKey.c_str()).Data()
                  );
                  SaveCanvas(c, outPng);

                  cout << "  [LeadJetResponse] Wrote: " << outPng << "\n";
                  effSummary.push_back("  - " + outPng);

                  for (auto* o : keep) delete o;
                }
              }

              effSummary.push_back("");

          }

          // ---------------------------------------------------------------------------
          // (3) JES3 tag-fraction 3x3 tables vs xJ (COUNT ratios, not shape-normalized)
          // ---------------------------------------------------------------------------
          {
            auto AlphaTag = [&](double aMax)->string
            {
              const int aInt = (int)std::lround(aMax * 100.0);
              return TString::Format("alphaLT0p%02d", aInt).Data();
            };

            auto Make3x3CountRatioTable_TH3xJ =
              [&](TH3* hNum, TH3* hDen,
                  const string& outPng,
                  const vector<string>& titleLines,
                  double alphaMax,
                  double yMin,
                  double yMax)
            {
              if (!hNum || !hDen) return;

              TCanvas c(
                TString::Format("c_tbl_cntRatio_%s_%s_%s",
                  ds.label.c_str(), rKey.c_str(), AlphaTag(alphaMax).c_str()
                ).Data(),
                "c_tbl_cntRatio", 1500, 1200
              );
              c.Divide(3,2, 0.001, 0.001);

              vector<TH1*> keep;

              for (int ib = 1; ib <= kNPtBins; ++ib)
              {
                c.cd(ib);

                gPad->SetLeftMargin(0.14);
                gPad->SetRightMargin(0.05);
                gPad->SetBottomMargin(0.14);
                gPad->SetTopMargin(0.10);

                TH1* n = ProjectY_AtXbin_AndAlphaMax_TH3(
                  hNum, ib,
                  alphaMax,
                  TString::Format("cntR_n_%s_%s_b%d", rKey.c_str(), AlphaTag(alphaMax).c_str(), ib).Data()
                );
                TH1* d = ProjectY_AtXbin_AndAlphaMax_TH3(
                  hDen, ib,
                  alphaMax,
                  TString::Format("cntR_d_%s_%s_b%d", rKey.c_str(), AlphaTag(alphaMax).c_str(), ib).Data()
                );

                if (n) { n->SetDirectory(nullptr); EnsureSumw2(n); }
                if (d) { d->SetDirectory(nullptr); EnsureSumw2(d); }

                if (!n || !d || d->Integral(1, d->GetNbinsX()) <= 0.0)
                {
                  if (n) delete n;
                  if (d) delete d;

                  TLatex t;
                  t.SetNDC(true);
                  t.SetTextFont(42);
                  t.SetTextSize(0.06);
                  t.DrawLatex(0.15, 0.55, "MISSING");
                  continue;
                }

                TH1* r = CloneTH1(
                  n,
                  TString::Format("cntR_r_%s_%s_b%d", rKey.c_str(), AlphaTag(alphaMax).c_str(), ib).Data()
                );
                if (!r)
                {
                  delete n;
                  delete d;
                  continue;
                }
                EnsureSumw2(r);
                r->Divide(n, d, 1.0, 1.0, "B");

                r->SetTitle("");
                r->GetXaxis()->SetTitle("x_{J#gamma}");
                r->GetYaxis()->SetTitle("Fraction");
                r->GetYaxis()->SetRangeUser(yMin, yMax);

                r->SetLineWidth(2);
                r->SetMarkerStyle(20);
                r->SetMarkerSize(0.95);
                r->Draw("E1");

                const string ptLab = AxisBinLabel(hDen->GetXaxis(), ib, "GeV", 0);

                TLatex tt;
                tt.SetNDC(true);
                tt.SetTextFont(42);
                tt.SetTextAlign(22);
                tt.SetTextSize(0.060);
                tt.DrawLatex(0.52, 0.95,
                  TString::Format("p_{T}^{#gamma} = %s  (R=%.1f)", ptLab.c_str(), R).Data()
                );

                keep.push_back(n);
                keep.push_back(d);
                keep.push_back(r);
              }

              // Put the common title block in the first pad (keeps the other pads uncluttered)
              c.cd(1);
              DrawLatexLines(0.16, 0.86, titleLines, 0.040, 0.050);

              SaveCanvas(c, outPng);

              for (auto* h : keep) delete h;
            };

            const vector<double> alphaMaxCuts = {0.20, 0.30, 0.40, 0.50};

            // (3a) Jet-tag conditional fraction vs xJ: truthTaggedPhoJet / truthPhoTagged
            if (H.hRecoTruthTagged_xJ && H.hRecoTruthPhoTagged_xJ)
            {
              for (double aMax : alphaMaxCuts)
              {
                const string outPng = JoinPath(
                  D.dirXJProjEffTagFractions3x3,
                  TString::Format("table3x3_tagFrac_jetMatchedOverPhoTagged_%s_%s.png",
                    rKey.c_str(), AlphaTag(aMax).c_str()
                  ).Data()
                );

                Make3x3CountRatioTable_TH3xJ(
                  H.hRecoTruthTagged_xJ, H.hRecoTruthPhoTagged_xJ,
                  outPng,
                  {"Tag-fraction vs x_{J#gamma} (COUNT ratio)",
                   "Numerator: RECO truth-#gamma + jet1 matched",
                   "Denominator: RECO truth-#gamma tagged",
                   TString::Format("#alpha < %.2f", aMax).Data()},
                  aMax,
                  0.0, 1.05
                );

                cout << "  [TagFractions] Wrote: " << outPng << "\n";
                effSummary.push_back("TagFractions_3x3:");
                effSummary.push_back("  - " + outPng);
              }
              effSummary.push_back("");
            }
            else
            {
              cout << ANSI_BOLD_YEL << "  [TagFractions] Missing truthTaggedPhoJet and/or truthPhoTagged TH3; skipping jetMatchedOverPhoTagged.\n" << ANSI_RESET;
              effSummary.push_back("TagFractions_3x3: MISSING inputs (jetMatchedOverPhoTagged skipped)");
              effSummary.push_back("");
            }

            // (3b) Photon-tag fraction vs xJ: truthPhoTagged / recoAll
            if (H.hRecoTruthPhoTagged_xJ && H.hReco_xJ)
            {
              for (double aMax : alphaMaxCuts)
              {
                const string outPng = JoinPath(
                  D.dirXJProjEffTagFractions3x3,
                  TString::Format("table3x3_tagFrac_phoTaggedOverRecoAll_%s_%s.png",
                    rKey.c_str(), AlphaTag(aMax).c_str()
                  ).Data()
                );

                Make3x3CountRatioTable_TH3xJ(
                  H.hRecoTruthPhoTagged_xJ, H.hReco_xJ,
                  outPng,
                  {"Tag-fraction vs x_{J#gamma} (COUNT ratio)",
                   "Numerator: RECO truth-#gamma tagged",
                   "Denominator: RECO (all)",
                   TString::Format("#alpha < %.2f", aMax).Data()},
                  aMax,
                  0.0, 1.05
                );

                cout << "  [TagFractions] Wrote: " << outPng << "\n";
                effSummary.push_back("TagFractions_3x3:");
                effSummary.push_back("  - " + outPng);
              }
              effSummary.push_back("");
            }
            else
            {
              cout << ANSI_BOLD_YEL << "  [TagFractions] Missing truthPhoTagged and/or recoAll TH3; skipping phoTaggedOverRecoAll.\n" << ANSI_RESET;
              effSummary.push_back("TagFractions_3x3: MISSING inputs (phoTaggedOverRecoAll skipped)");
              effSummary.push_back("");
            }

            // (3c) Optional: doubly-tagged fraction vs xJ: truthTaggedPhoJet / recoAll
            if (H.hRecoTruthTagged_xJ && H.hReco_xJ)
            {
              for (double aMax : alphaMaxCuts)
              {
                const string outPng = JoinPath(
                  D.dirXJProjEffTagFractions3x3,
                  TString::Format("table3x3_tagFrac_jetMatchedOverRecoAll_%s_%s.png",
                    rKey.c_str(), AlphaTag(aMax).c_str()
                  ).Data()
                );

                Make3x3CountRatioTable_TH3xJ(
                  H.hRecoTruthTagged_xJ, H.hReco_xJ,
                  outPng,
                  {"Tag-fraction vs x_{J#gamma} (COUNT ratio)",
                   "Numerator: RECO truth-#gamma + jet1 matched",
                   "Denominator: RECO (all)",
                   TString::Format("#alpha < %.2f", aMax).Data()},
                  aMax,
                  0.0, 1.05
                );

                cout << "  [TagFractions] Wrote: " << outPng << "\n";
                effSummary.push_back("TagFractions_3x3:");
                effSummary.push_back("  - " + outPng);
              }
              effSummary.push_back("");
            }
          }

          // Write a human-readable summary right next to the plots
          const string sumTxt = JoinPath(D.dirXJProjEff, "summary_efficiency_dashboard.txt");
          WriteTextFile(sumTxt, effSummary);
          cout << "  Wrote summary: " << sumTxt << "\n";
        }

        WriteTextFile(JoinPath(D.dirSumm, "summary_JES3_means.txt"), sumLines);
    };

    // =============================================================================
    // Per-rKey JES3: keep existing outputs + add:
    //  - 3x3 xJ tables (linear) [EXISTING]
    //  - 3x3 xJ tables (logy)   [NEW, additional]
    // =============================================================================
    for (const auto& rKey : kRKeys)
    {
      ProcessOneRKey(rKey);
    }

    // r02/r04 overlays (SIM only) - keep under xJ_fromJES3/Overlay/
    if (ds.isSim)
    {
      const string dirSharedXJProj        = JoinPath(outDir, "xJ_fromJES3");
      const string dirSharedXJProjOverlay = JoinPath(dirSharedXJProj, "Overlay");
      EnsureDir(dirSharedXJProj);
      EnsureDir(dirSharedXJProjOverlay);

      JES3_R02R04Overlays_MaybeRun(ds, dirSharedXJProjOverlay);

        // NEW: RECO-only Δφ-cut overlay table (π/2 vs 7π/8), r04 only
        // Output: <outDir>/r04/xJ_fromJES3/RECO/table3x2_overlay_integratedAlpha_dPhiCuts.png
        JES3_RecoDeltaPhiCutOverlay_MaybeRun(ds, outDir);

        // NEW: RECO-only p_{T}^{jet,min} compare overlay (3 vs 5 vs 10), r02/r04/r06
        // Output (per radius):
        //   <outDir>/<rKey>/xJ_fromJES3/RECO/table3x2_overlay_integratedAlpha_pTminCompare.png
        JES3_RecoBinningOverlay_MaybeRun(ds, outDir);

        // NEW: Sam vs Justin unsmear overlay (hard-coded inputs) for pTminJet3 + 7pi/8, r02/r04/r06
        // Output (per radius):
        //   <outDir>/<rKey>/xJ_fromJES3/RECO/SamVsJustin_pTminJet3_7piOver8/overlay_SamVsJustin_JES3_RECO_pTgamma_13_15_<rKey>.png
        JES3_SamVsJustinUnsmearOverlay_MaybeRun(ds, outDir);
    }

    // =============================================================================
    // In-situ JES3 residual calibration (DATA vs SIM) using reco-only JES3 TH3s.
    // Refactored: implemented as a helper function defined above RunJES3QA.
    // (Functionality and outputs preserved exactly.)
    // =============================================================================
    JES3_InSituResidualCalibration_MaybeRun(ds);
}
